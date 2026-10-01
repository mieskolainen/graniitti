# Create and submit Condor test suites and measurement icepacks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import json
import math
import os
import pathlib
import re
import shutil
import subprocess
import sys

from core.io.files import dated_directory

from tests.condor import job_runner
from tests.condor.job_runner import POOL_COUNT, ROOT, write_json


# Read Condor resource settings
def resource_settings():
    return json.loads((ROOT / "tests/condor/SETTINGS.json").read_text(encoding="ascii"))


# Parse positive integer command line values
def positive_int(value: str) -> int:
    result = int(value)
    if result <= 0:
        raise argparse.ArgumentTypeError("value must be positive")
    return result


# Parse a worker count large enough to serve every array
def worker_count(value: str) -> int:
    result = positive_int(value)
    if result < POOL_COUNT:
        raise argparse.ArgumentTypeError(f"value must be at least {POOL_COUNT}")
    return result


# Parse suite selection and pass explicit trailing options to condor_submit_dag
def parse_args(argv=None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Submit tests or icepacks to Condor")
    parser.add_argument("suite", nargs="?", choices=("all", "cpp", "pytest", "studies", "icepacks"), default="all")
    parser.add_argument("--run-name", help="unique run name")
    parser.add_argument("--output-root", type=pathlib.Path,
                        help="output root (default: figs/icepack for icepacks, runs/tests otherwise)")
    parser.add_argument("--study-events", type=positive_int, default=5000)
    parser.add_argument("--build-jobs", type=int, choices=range(1, 9), default=resource_settings()["jobs"]["build"]["request_cpus"])
    parser.add_argument("--max-jobs", type=worker_count, default=200,
                        help="maximum materialized workers shared across test arrays")
    parser.add_argument("--shard-waves", type=positive_int, default=4)
    parser.add_argument("--timings", type=pathlib.Path, help="timing history, defaulting to the newest previous run")
    parser.add_argument("--dry-run", action="store_true", help="prepare the run without executing or submitting")
    parser.add_argument("--validate", action="store_true", help="generate the DAGMan submission before submitting")
    parser.add_argument("--name", help="folder relative to icepack/, searched recursively including nonmeasurement cards")
    parser.add_argument("--local", action="store_true", help="run icepacks sequentially on this machine")
    parser.add_argument("--report", type=pathlib.Path, help="refresh an existing icepack result table")
    argv = list(sys.argv[1:] if argv is None else argv)
    split = argv.index("--") if "--" in argv else len(argv)
    args = parser.parse_args(argv[:split])
    args.output_root = args.output_root or ROOT / ("figs/icepack" if args.suite == "icepacks" else "runs/tests")
    if args.build_jobs not in range(1, 9):
        parser.error("build request_cpus must be between 1 and 8")
    args.condor_args = argv[split + 1:]
    if args.suite != "icepacks" and (args.name or args.local or args.report):
        parser.error("--name, --local and --report require the icepacks suite")
    if args.local and (args.validate or args.condor_args):
        parser.error("--local cannot submit or validate a Condor DAG")
    return args


# Validate one run name for safe use in paths
def validate_run_name(value: str) -> str:
    if value in (".", "..") or not re.fullmatch(r"[A-Za-z0-9_.-]+", value):
        raise ValueError("run name may contain only letters, digits, underscores, dots, and hyphens")
    return value


# Compute the configured Pytest collection command
def python_collection_command() -> list[str]:
    return [
        sys.executable,
        "-m",
        "pytest",
        "--collect-only",
        "-q",
        "--run-integration",
        "--run-physics",
        "--LOOPSCREEN",
        "1",
        "-p",
        "tests.condor.pytest_plugin",
    ]


# Collect every exact Python test node through Pytest
def collect_python_tests(run_dir: pathlib.Path) -> list[str]:
    collection_path = run_dir / "manifests" / "python_nodes.json"
    environment = os.environ.copy()
    environment["GRANIITTI_PYTEST_COLLECTION"] = str(collection_path)
    environment.pop("GRANIITTI_PYTEST_TIMINGS", None)
    environment["PYTHONDONTWRITEBYTECODE"] = "1"
    command = python_collection_command()
    result = subprocess.run(
        command,
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    collection_log = result.stdout + result.stderr
    (run_dir / "logs" / "python_collection.out").write_text(collection_log, encoding="utf-8")
    if result.returncode != 0:
        raise RuntimeError(
            f"Python test collection failed, see {run_dir / 'logs' / 'python_collection.out'}"
        )
    with collection_path.open(encoding="utf-8") as handle:
        nodes = json.load(handle)
    if not isinstance(nodes, list) or not all(isinstance(node, str) and node for node in nodes):
        raise RuntimeError("Python test collection produced an invalid node list")
    if len(nodes) != len(set(nodes)):
        raise RuntimeError("Python test collection produced duplicate node IDs")
    return sorted(nodes)


# Compute every canonical study launcher relative to the repository
def discover_studies(root: pathlib.Path = ROOT) -> list[str]:
    study_root = root / "tests" / "physics" / "studies"
    return sorted(path.relative_to(root).as_posix() for path in study_root.rglob("run.sh"))


# Load explicit or latest timing history
def load_timings(output_root: pathlib.Path, explicit_path: pathlib.Path | None) -> dict[str, float]:
    candidates: list[pathlib.Path]
    if explicit_path is not None:
        candidates = [explicit_path.expanduser().resolve()]
    else:
        candidates = sorted(
            output_root.glob("*/timings.json"),
            key=lambda path: path.stat().st_mtime,
            reverse=True,
        )
    if not candidates:
        return {}
    with candidates[0].open(encoding="ascii") as handle:
        raw_timings = json.load(handle)
    timings: dict[str, float] = {}
    for key, value in raw_timings.items():
        duration = float(value)
        if duration >= 0.0 and math.isfinite(duration):
            timings[str(key)] = max(duration, 0.01)
    return timings


# Estimate one Python node duration from history and suite type
def estimate_python_duration(node: str, timings: dict[str, float]) -> float:
    history = timings.get(f"python::{node}")
    if history is not None:
        return history
    path = node.split("::", maxsplit=1)[0]
    if path.startswith("tests/physics/validation/"):
        return 1800.0
    if path.startswith(("tests/technical/integration/", "tests/technical/basic/")):
        return 600.0
    if path.startswith(("tests/technical/inference/", "tests/technical/tuning/")):
        return 30.0
    return 3.0


# Estimate one study duration from history
def estimate_study_duration(path: str, timings: dict[str, float]) -> float:
    return timings.get(f"study::{path}", 1800.0)


# Pack weighted Python nodes into balanced longest first shards
def shard_python_tests(
    nodes: list[str],
    timings: dict[str, float],
    shard_count: int,
) -> list[dict]:
    shards = [{"estimated_seconds": 0.0, "nodes": []} for _ in range(shard_count)]
    aggregate = "tests/physics/validation/test_icepacks.py::test_integrated_xs_table"
    integrated = {aggregate}
    if aggregate in nodes:
        from tests.physics.validation.test_icepacks import ICEPACK_ROOT, INTEGRATED_DATASETS

        integrated.update(
            "tests/physics/validation/test_icepacks.py::test_icepack_event_flow"
            f"[{path.parent.relative_to(ICEPACK_ROOT).as_posix()}]"
            for path in INTEGRATED_DATASETS
        )
    # Keep overlapping generation in one process so the table reuses the reports
    groups = {}
    for node in nodes:
        groups.setdefault(aggregate if node in integrated else node, []).append(node)
    weighted_nodes = sorted(
        ((sum(estimate_python_duration(node, timings) for node in group), sorted(group))
         for group in groups.values()),
        reverse=True,
    )
    for duration, group in weighted_nodes:
        shard = min(shards, key=lambda value: (value["estimated_seconds"], len(value["nodes"])))
        shard["nodes"].extend(group)
        shard["estimated_seconds"] += duration
    return [shard for shard in shards if shard["nodes"]]


# Write balanced Python shard task descriptions
def write_python_tasks(run_dir: pathlib.Path, shards: list[dict]) -> None:
    for index, shard in enumerate(shards):
        relative = pathlib.Path("tasks") / f"python.{index:06d}.json"
        payload = {
            "category": "python",
            "estimated_seconds": shard["estimated_seconds"],
            "item": relative.as_posix(),
            "members": shard["nodes"],
        }
        write_json(run_dir / relative, payload)


# Write one isolated task description for every study
def write_study_tasks(
    run_dir: pathlib.Path,
    studies: list[str],
    timings: dict[str, float],
) -> None:
    for index, study in enumerate(studies):
        relative = pathlib.Path("tasks") / f"study.{index:06d}.json"
        payload = {
            "category": "study",
            "estimated_seconds": estimate_study_duration(study, timings),
            "item": study,
            "members": [study],
        }
        write_json(run_dir / relative, payload)


# Reject manifest values that the Condor item parser cannot represent safely
def validate_manifest_items(items: list[str], category: str) -> None:
    if not items:
        raise RuntimeError(f"no {category} items were discovered")
    invalid = [item for item in items if not item or any(character.isspace() for character in item)]
    if invalid:
        raise RuntimeError(f"{category} manifest paths may not contain whitespace: {invalid}")
    if len(items) != len(set(items)):
        raise RuntimeError(f"duplicate {category} manifest items were discovered")


# Escape one value for a DAG VARS declaration
def dag_quote(value: pathlib.Path | str | int) -> str:
    text = str(value)
    return text.replace("\\", "\\\\").replace('"', '\\"')


# Write a submit description through the shared environment launcher
def write_submit(run_dir, name, module, arguments, *, queue="1", **extra):
    path = run_dir / "condor" / f"{name}.sub"
    path.parent.mkdir(parents=True, exist_ok=True)
    settings = resource_settings()
    fields = {
        "universe": "vanilla",
        "executable": ROOT / "tests/condor/run_job.sh",
        "arguments": f'"{ROOT} {module} {arguments}"',
        "initialdir": ROOT,
        "output": run_dir / "logs" / f"{name}.$(ClusterId).$(ProcId).out",
        "error": run_dir / "logs" / f"{name}.$(ClusterId).$(ProcId).err",
        "log": run_dir / "logs" / f"{name}.condor.log",
        "getenv": "True",
        "should_transfer_files": "NO",
        **settings["defaults"],
        **settings["jobs"][name],
        "notification": "Never",
        **extra,
    }
    path.write_text("".join(f"{key} = {value}\n" for key, value in fields.items()) + f"queue {queue}\n",
                    encoding="ascii")
    return path


# Write the selected test arrays with an optional build and an unconditional final verifier
def write_dag(path, run_dir, array_limits, build_jobs=None, build=True):
    work = write_submit(run_dir, "work", "job_runner", f"work $(test_index) $(test_item) {run_dir}",
                        queue=f"test_index,test_item from {run_dir}/manifests/work_$(POOL).tsv",
                        max_materialize="$(MAX_MATERIALIZE)")
    config = run_dir / "condor/dagman.config"
    config.write_text("DAGMAN_PROHIBIT_MULTI_JOBS = False\n", encoding="ascii")
    lines = [f"CONFIG {config}"]
    names = [f"WORK_{index}" for index in range(len(array_limits))]
    if build:
        job = write_submit(run_dir, "build", "job_runner", f"build 000000 build {run_dir}",
                           **({"request_cpus": build_jobs} if build_jobs is not None else {}))
        lines.append(f"JOB BUILD {job}")
    for index, (name, limit) in enumerate(zip(names, array_limits, strict=True)):
        lines.extend([f"JOB {name} {work}", f'VARS {name} POOL="{index}" MAX_MATERIALIZE="{limit}"',
                      f"RETRY {name} 1 UNLESS-EXIT 64"])
    if build:
        lines.extend([f"PARENT BUILD CHILD {' '.join(names)}", "RETRY BUILD 1 UNLESS-EXIT 64"])
    final = write_submit(run_dir, "verifier", "verify", f"{run_dir} $(DAG_STATUS) $(FAILED_COUNT)")
    path.write_text("\n".join([*lines, f"FINAL VERIFIER {final}"]) + "\n", encoding="ascii")


# Write one job per active icepack and always finish with its measurement table
def write_icepack_dag(path, run_dir, entries):
    job = write_submit(run_dir, "icepacks", "icepacks", f"--run {run_dir} $(index)",
                       **{"+JobDescription": '"$(entry)"'})
    final = write_submit(run_dir, "summary", "icepacks", f"--report {run_dir} --final")
    lines = [f'JOB D{i} {job}\nVARS D{i} index="{i}" entry="{dag_quote(item["entry"])}"'
             for i, item in enumerate(entries) if item["active"]]
    path.write_text("\n".join([*lines, f"FINAL SUMMARY {final}"]) + "\n", encoding="ascii")


# Compute the current Git revision when available
def git_revision() -> str | None:
    result = subprocess.run(
        ["git", "rev-parse", "HEAD"],
        cwd=ROOT,
        check=False,
        capture_output=True,
        text=True,
    )
    return result.stdout.strip() if result.returncode == 0 else None


# Prepare only the selected test categories and keep the total worker budget fixed
def prepare_tests(run_dir, args):
    categories = {"all": ["cpp", "python", "study"], "cpp": ["cpp"],
                  "pytest": ["python"], "studies": ["study"]}[args.suite]
    build = "cpp" in categories
    timings = load_timings(run_dir.parent, args.timings)
    nodes = collect_python_tests(run_dir) if "python" in categories else []
    if "python" in categories and not nodes:
        raise RuntimeError("no Python tests were discovered")
    studies = discover_studies() if "study" in categories else []
    if "study" in categories:
        validate_manifest_items(studies, "study")
    shards = shard_python_tests(nodes, timings, min(len(nodes), args.max_jobs * args.shard_waves))
    write_python_tasks(run_dir, shards)
    write_study_tasks(run_dir, studies, timings)
    pools = POOL_COUNT if build else min(POOL_COUNT, len(shards) + len(studies))
    if shards or studies:
        job_runner.write_work_manifests(run_dir, pools)
    base, remainder = divmod(args.max_jobs, pools)
    limits = [base + int(index < remainder) for index in range(pools)]
    write_json(run_dir / "scheduling_timings.json", timings)
    write_dag(run_dir / "condor/run.dag", run_dir, limits, args.build_jobs, build)
    return {
        "categories": categories, "build": build, "array_limits": limits,
        "build_jobs": args.build_jobs, "max_jobs": args.max_jobs, "loopscreen": 1,
        "python_shards": len(shards), "python_tests": len(nodes), "shard_waves": args.shard_waves,
        "studies": len(studies), "study_events": args.study_events, "timing_history_entries": len(timings),
    }


# Create every suite under one output root with separate configuration and writable outputs
def create_run(args: argparse.Namespace) -> pathlib.Path:
    run_name = validate_run_name(args.run_name) if args.run_name else None
    output_root = args.output_root.expanduser().resolve()
    if any(character.isspace() for character in f"{ROOT}{output_root}"):
        raise ValueError("repository and output paths may not contain whitespace")
    config = {"suite": args.suite, "run_name": run_name, "repo_root": str(ROOT), "git_revision": git_revision()}
    if args.suite == "icepacks":
        from tests.condor import icepacks

        config.update(entries=icepacks.datasets(args.name), steering=icepacks.steering())
    if run_name is None:
        run_dir = dated_directory(output_root, prefix=f"{args.suite}.")
    else:
        run_dir = output_root / run_name
        run_dir.mkdir(parents=True, exist_ok=False)
    config["run_name"] = run_dir.name
    for name in ("condor", "junit", "logs", "manifests", "outputs", "results", "tasks", "tmp", "work"):
        (run_dir / name).mkdir()
    if args.suite == "icepacks":
        write_icepack_dag(run_dir / "condor/run.dag", run_dir, config["entries"])
    else:
        config.update(prepare_tests(run_dir, args))
    write_json(run_dir / "config.json", config)
    return run_dir


# Validate or submit one generated DAG with all scheduler outputs inside the run
def submit_dag(run_dir: pathlib.Path, validate: bool, condor_args=()) -> None:
    dag_executable = shutil.which("condor_submit_dag")
    if dag_executable is None:
        raise RuntimeError("condor_submit_dag is not available on PATH")
    dag_path = run_dir / "condor/run.dag"
    command = [dag_executable, *condor_args]
    if validate:
        submit_executable = shutil.which("condor_submit")
        if submit_executable is None:
            raise RuntimeError("condor_submit is not available on PATH")
        subprocess.run([*command, "-no_submit", str(dag_path)], cwd=dag_path.parent, check=True)
        subprocess.run([submit_executable, f"{dag_path}.condor.sub"], cwd=dag_path.parent, check=True)
    else:
        subprocess.run([*command, str(dag_path)], cwd=dag_path.parent, check=True)


# Generate the selected suite, refresh its report or run measurement icepacks locally
def main() -> int:
    args = parse_args()
    if args.suite == "icepacks":
        from tests.condor import icepacks

        if args.report:
            return icepacks.summary(args.report.expanduser().resolve())
    run_dir = create_run(args)
    config = job_runner.load_config(run_dir)
    print(f"Generated test DAG: {run_dir / 'condor/run.dag'}")
    if args.suite == "icepacks":
        icepacks.summary(run_dir)
        print(f"Run report: {run_dir / 'report.md'}")
    else:
        print(f"Run report: {run_dir / 'report.txt'}")
        print(f"Array limits: {config['array_limits']} (maximum {config['max_jobs']})")
        print(f"Collected {config['python_tests']} Python tests into {config['python_shards']} balanced shards")
    if args.dry_run:
        return 0
    if args.local:
        for index, item in enumerate(config["entries"]):
            if item["active"]:
                icepacks.worker(run_dir, index)
        return icepacks.summary(run_dir, final=True)
    submit_dag(run_dir, args.validate, args.condor_args)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
