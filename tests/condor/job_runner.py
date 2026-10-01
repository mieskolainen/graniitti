# Run one build, CTest, Python, or study array item
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import argparse
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import time
import traceback

ROOT = pathlib.Path(__file__).resolve().parents[2]
POOL_COUNT = 3


# Parse one array item invocation
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("category", choices=("build", "work"))
    parser.add_argument("index")
    parser.add_argument("item")
    parser.add_argument("run_dir", type=pathlib.Path)
    return parser.parse_args()


# Load the immutable run configuration
def load_config(run_dir: pathlib.Path) -> dict:
    with (run_dir / "config.json").open(encoding="ascii") as handle:
        return json.load(handle)


# Write JSON through an atomic rename
def write_json(path: pathlib.Path, payload: dict) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    temporary.write_text(json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="ascii")
    temporary.replace(path)


# Write indexed manifest records through an atomic rename
def write_manifest_records(path: pathlib.Path, records: list[tuple[str, str]]) -> None:
    temporary = path.with_name(f"{path.name}.tmp.{os.getpid()}")
    lines = [f"{index} {item}" for index, item in records]
    temporary.write_text("\n".join(lines) + "\n", encoding="ascii")
    temporary.replace(path)


# Run one command and return its exit code
def run_command(
    command: list[str],
    environment: dict[str, str],
    input_text: str | None = None,
    working_directory: pathlib.Path = ROOT,
) -> int:
    print(f"Running: {' '.join(command)}", flush=True)
    result = subprocess.run(
        command,
        cwd=working_directory,
        env=environment,
        input=input_text,
        text=True,
        check=False,
    )
    return result.returncode


# Discover all configured CTest targets after a successful build
def discover_cpp_tests(environment: dict[str, str]) -> list[str]:
    result = subprocess.run(
        ["ctest", "--test-dir", "build", "--show-only=json-v1"],
        cwd=ROOT,
        env=environment,
        check=False,
        capture_output=True,
        text=True,
    )
    if result.returncode != 0:
        print(result.stdout, end="")
        print(result.stderr, end="", file=sys.stderr)
        raise RuntimeError("CTest discovery failed")
    payload = json.loads(result.stdout)
    tests = sorted(test["name"] for test in payload.get("tests", []))
    if not tests:
        raise RuntimeError("CTest discovery returned no tests")
    if len(tests) != len(set(tests)):
        raise RuntimeError("CTest discovery returned duplicate test names")
    return tests


# Estimate one newly discovered CTest target from timing history
def estimate_cpp_duration(name: str, timings: dict[str, float]) -> float:
    history = timings.get(f"cpp::{name}")
    if history is not None:
        return float(history)
    if name == "cpp.test_neurojac":
        return 600.0
    if name in ("cpp.test_models", "cpp.test_hard_pomeron"):
        return 180.0
    if name.startswith("cpp.clang_tidy."):
        return 30.0
    return 20.0


# Add CTest task descriptions after build time discovery
def write_cpp_tasks(run_dir: pathlib.Path, tests: list[str], timings: dict[str, float]) -> None:
    for index, test in enumerate(tests):
        relative = pathlib.Path("tasks") / f"cpp.{index:06d}.json"
        payload = {
            "category": "cpp",
            "estimated_seconds": estimate_cpp_duration(test, timings),
            "item": test,
            "members": [test],
        }
        write_json(run_dir / relative, payload)


# Load and validate one generic work task description
def load_work_task(path: pathlib.Path, run_dir: pathlib.Path) -> dict:
    resolved = path.resolve()
    task_root = (run_dir / "tasks").resolve()
    if not resolved.is_relative_to(task_root) or not resolved.is_file():
        raise RuntimeError(f"invalid work task path: {path}")
    with resolved.open(encoding="ascii") as handle:
        task = json.load(handle)
    if task.get("category") not in ("cpp", "python", "study"):
        raise RuntimeError(f"invalid work task category: {path}")
    if not isinstance(task.get("item"), str) or not task["item"]:
        raise RuntimeError(f"invalid work task item: {path}")
    members = task.get("members")
    if (
        not isinstance(members, list)
        or not members
        or not all(isinstance(member, str) and member for member in members)
    ):
        raise RuntimeError(f"invalid work task members: {path}")
    duration = float(task.get("estimated_seconds", 0.0))
    if duration < 0.0:
        raise RuntimeError(f"invalid work task estimate: {path}")
    task["estimated_seconds"] = duration
    task["path"] = resolved.relative_to(run_dir).as_posix()
    return task


# Mix all task types into balanced work array manifests
def write_work_manifests(run_dir: pathlib.Path, pool_count: int = POOL_COUNT) -> None:
    tasks = [load_work_task(path, run_dir) for path in sorted((run_dir / "tasks").glob("*.json"))]
    if not tasks:
        raise RuntimeError("no work tasks were prepared")
    pools = [{"load": 0.0, "records": []} for _ in range(pool_count)]
    category_records: dict[str, list[tuple[str, str]]] = {
        "cpp": [],
        "python": [],
        "study": [],
    }
    ordered = sorted(
        tasks, key=lambda task: (task["estimated_seconds"], task["path"]), reverse=True
    )
    for number, task in enumerate(ordered):
        index = f"{number:06d}"
        pool = min(pools, key=lambda value: (value["load"], len(value["records"])))
        pool["records"].append((index, task["path"]))
        pool["load"] += task["estimated_seconds"]
        category_records[task["category"]].append((index, task["item"]))
    for pool_index, pool in enumerate(pools):
        write_manifest_records(run_dir / "manifests" / f"work_{pool_index}.tsv", pool["records"])
    for category, records in category_records.items():
        write_manifest_records(run_dir / "manifests" / f"{category}.tsv", records)
    loads = ", ".join(f"{pool['load']:.1f}" for pool in pools)
    print(f"Prepared {len(tasks)} tasks with estimated pool seconds: {loads}", flush=True)


# Configure and build the C++ test targets
def run_build(
    config: dict,
    run_dir: pathlib.Path,
    environment: dict[str, str],
) -> tuple[int, list[list[str]], list[str]]:
    commands = [
        ["cmake", "-S", ".", "-B", "build", "-DWITH_TEST=ON"],
        ["cmake", "--build", "build", f"-j{config['build_jobs']}"],
    ]
    for command in commands:
        exit_code = run_command(command, environment)
        if exit_code != 0:
            return exit_code, commands, ["build"]
    tests = discover_cpp_tests(environment)
    with (run_dir / "scheduling_timings.json").open(encoding="ascii") as handle:
        timings = json.load(handle)
    write_cpp_tasks(run_dir, tests, timings)
    write_work_manifests(run_dir, len(config["array_limits"]))
    print(f"Discovered {len(tests)} CTest targets", flush=True)
    return 0, commands, ["build"]


# Run one exact CTest target
def run_cpp(item: str, environment: dict[str, str]) -> tuple[int, list[list[str]], list[str]]:
    command = [
        "ctest",
        "--test-dir",
        "build",
        "--output-on-failure",
        "--no-tests=error",
        "-R",
        f"^{re.escape(item)}$",
    ]
    return run_command(command, environment), [command], [item]


# Run every test case assigned to one balanced Python shard
def run_python(
    nodes: list[str],
    index: str,
    config: dict,
    run_dir: pathlib.Path,
    environment: dict[str, str],
) -> tuple[int, list[list[str]], list[str]]:
    junit_path = run_dir / "junit" / f"python.{index}.xml"
    timing_path = run_dir / "junit" / f"python.{index}.timings.json"
    environment.pop("GRANIITTI_PYTEST_COLLECTION", None)
    environment["GRANIITTI_PYTEST_TIMINGS"] = str(timing_path)
    environment["GRANIITTI_TEST_OUTPUT_DIR"] = str(run_dir / "outputs" / "python" / index)
    command = [
        sys.executable,
        "-m",
        "pytest",
        "-q",
        "-s",
        *nodes,
        "--run-integration",
        "--run-physics",
        "--LOOPSCREEN",
        str(config["loopscreen"]),
        "-rP",
        "--junitxml",
        str(junit_path),
        "-p",
        "tests.condor.pytest_plugin",
    ]
    return run_command(command, environment), [command], nodes


# Link existing cache inputs while keeping the destination directory writable
def link_cache_inputs(source: pathlib.Path, destination: pathlib.Path) -> None:
    destination.mkdir()
    for entry in source.iterdir():
        (destination / entry.name).symlink_to(entry, target_is_directory=entry.is_dir())


# Link the test tree while copying the writable physics studies
def link_study_tests(source: pathlib.Path, destination: pathlib.Path) -> None:
    destination.mkdir()
    for entry in source.iterdir():
        target = destination / entry.name
        if entry.name != "physics":
            target.symlink_to(entry, target_is_directory=entry.is_dir())
            continue
        target.mkdir()
        for physics_entry in entry.iterdir():
            physics_target = target / physics_entry.name
            if physics_entry.name == "studies":
                shutil.copytree(physics_entry, physics_target)
            else:
                physics_target.symlink_to(
                    physics_entry,
                    target_is_directory=physics_entry.is_dir(),
                )


# Build an isolated filesystem view for one parallel study
def prepare_study_workspace(run_dir: pathlib.Path, index: str) -> pathlib.Path:
    workspace = run_dir / "work" / "study" / index
    if workspace.is_dir():
        return workspace
    temporary = workspace.with_name(f"{workspace.name}.prepare.{os.getpid()}")
    temporary.mkdir(parents=True, exist_ok=False)
    excluded = {
        ".git",
        "build",
        "eikonal",
        "figs",
        "output",
        "runs",
        "sudakov",
        "tests",
        "tmp",
        "vgrid",
    }
    for source in ROOT.iterdir():
        if source.name in excluded:
            continue
        (temporary / source.name).symlink_to(source, target_is_directory=source.is_dir())
    for name in ("figs", "output", "tmp"):
        (temporary / name).mkdir()
    shutil.copytree(ROOT / "vgrid", temporary / "vgrid")
    link_cache_inputs(ROOT / "eikonal", temporary / "eikonal")
    link_cache_inputs(ROOT / "sudakov", temporary / "sudakov")
    link_study_tests(ROOT / "tests", temporary / "tests")
    temporary.replace(workspace)
    return workspace


# Run one study in generation and analysis mode
def run_study(
    item: str,
    index: str,
    config: dict,
    run_dir: pathlib.Path,
    environment: dict[str, str],
) -> tuple[int, list[list[str]], list[str]]:
    study_path = (ROOT / item).resolve()
    study_root = (ROOT / "tests" / "physics" / "studies").resolve()
    if not study_path.is_relative_to(study_root) or not study_path.is_file():
        raise RuntimeError(f"invalid study path: {item}")
    workspace = prepare_study_workspace(run_dir, index)
    output_dir = run_dir / "outputs" / "study" / index
    output_dir.mkdir(parents=True, exist_ok=True)
    matplotlib_dir = run_dir / "tmp" / "matplotlib" / index
    matplotlib_dir.mkdir(parents=True, exist_ok=True)
    environment.pop("EVENTS", None)
    environment["NEVENTS"] = str(config["study_events"])
    environment["STUDY_ACTION"] = "both"
    environment["GRANIITTI_TEST_OUTPUT_DIR"] = str(output_dir)
    environment["MPLCONFIGDIR"] = str(matplotlib_dir)
    if "process_tables" in study_path.parts:
        environment["LOGDIR"] = str(output_dir / "logs")
        environment["TABLE"] = str(output_dir / "cross_sections.txt")
    command = ["bash", item]
    exit_code = run_command(
        command,
        environment,
        working_directory=workspace,
    )
    return exit_code, [command], [item]


# Dispatch one requested array item
def run_item(
    args: argparse.Namespace,
    config: dict,
    environment: dict[str, str],
) -> tuple[int, list[list[str]], list[str], str, str]:
    if args.category == "build":
        exit_code, commands, members = run_build(config, args.run_dir, environment)
        return exit_code, commands, members, "build", "build"
    task = load_work_task(args.run_dir / args.item, args.run_dir)
    category = task["category"]
    item = task["item"]
    if category == "cpp":
        exit_code, commands, members = run_cpp(item, environment)
    elif category == "python":
        exit_code, commands, members = run_python(
            task["members"], args.index, config, args.run_dir, environment
        )
    else:
        exit_code, commands, members = run_study(
            item, args.index, config, args.run_dir, environment
        )
    return exit_code, commands, members, category, item


# Record test results and stop dependent jobs if the build fails
def main() -> int:
    args = parse_args()
    args.run_dir = args.run_dir.resolve()
    if not args.index.isdigit():
        raise SystemExit("array index must be numeric")
    config = load_config(args.run_dir)
    environment = os.environ.copy()
    environment["MPLBACKEND"] = "Agg"
    environment["PYTHONUNBUFFERED"] = "1"
    temporary_dir = args.run_dir / "tmp" / args.category / args.index
    temporary_dir.mkdir(parents=True, exist_ok=True)
    environment["TMPDIR"] = str(temporary_dir)

    started = time.time()
    exit_code = 70
    commands: list[list[str]] = []
    members: list[str] = [args.item]
    result_category = args.category
    result_item = args.item
    error: str | None = None
    try:
        exit_code, commands, members, result_category, result_item = run_item(
            args, config, environment
        )
    except Exception:
        error = traceback.format_exc()
        print(error, file=sys.stderr, flush=True)
    ended = time.time()
    status = {
        "category": result_category,
        "commands": commands,
        "duration_seconds": ended - started,
        "ended_unix": ended,
        "error": error,
        "exit_code": exit_code,
        "index": args.index,
        "item": result_item,
        "members": members,
        "outcome": "PASS" if exit_code == 0 else "FAIL",
        "started_unix": started,
    }
    status_path = args.run_dir / "results" / result_category / f"{args.index}.json"
    try:
        write_json(status_path, status)
    except Exception:
        traceback.print_exc()
        return 64
    print(f"Recorded {status['outcome']}: {status_path}", flush=True)
    return exit_code if args.category == "build" else 0


if __name__ == "__main__":
    raise SystemExit(main())
