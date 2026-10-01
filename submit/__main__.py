# Prepare and submit resolved icetune campaigns
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import json
import os
import pathlib
import shlex
import shutil
import subprocess

import yaml
from core.io.files import ensure_dir

from submit import CAMPAIGN_DIR, runtime, schedulers
from submit import campaign as campaigns


# Parse one Ray campaign submission
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--campaign", required=True, help="campaign name from the YAML catalog")
    parser.add_argument("--scheduler", choices=sorted(schedulers.SUPPORTED_SCHEDULERS), default="local")
    parser.add_argument("--catalog", type=pathlib.Path, default=CAMPAIGN_DIR / "campaigns.yml")
    parser.add_argument("--run-name", default=None, help="override RUN_NAME")
    parser.add_argument(
        "--workers", type=int, default=None, help="number of Ray worker allocations, excluding the dedicated head")
    for role in ("init", "head", "worker"):
        for resource in ("cpu", "gpu"):
            parser.add_argument(f"--{role}-{resource}", type=int, default=None, help=f"override ray.{role}_{resource}")
    parser.add_argument("--max-concurrent", type=int, default=None, help="override MAX_CONCURRENT")
    parser.add_argument("--num-trials", type=int, default=None, help="override total NUM_TRIALS")
    parser.add_argument("--rand-trials", type=int, default=None, help="override RAND_TRIALS warm-up")
    parser.add_argument("--rngseed", type=int, default=None, help="override optimizer.rngseed")
    parser.add_argument("--no-ampfit-reuse", action="store_true", help="rebuild ampfit banks instead of reusing current or previous banks")
    parser.add_argument("--repo-dir", type=pathlib.Path, default=pathlib.Path.cwd(), help="simulator working directory")
    parser.add_argument("--output-dir", type=pathlib.Path, default=None, help="generated scheduler submission directory")
    parser.add_argument("--storage-path", type=pathlib.Path, default=None,
        help="persistent Ray Tune storage which shared mode exposes on every node", )
    parser.add_argument("--upload-runtime", action="store_true",
        help="upload the checkout through Ray instead of requiring a shared path", )
    parser.add_argument(
        "--lxplus-work-dir", type=pathlib.Path, default=None, help="lxplus work directory, defaults to runs/icetune/<run>"
    )
    parser.add_argument("--address", default=None, help="existing Ray head address for --scheduler external")
    parser.add_argument("--partition", default=None, help="Slurm partition")
    parser.add_argument("--queue", default=None, help="PBS queue")
    parser.add_argument("--set", dest="environment_overrides", action="append", default=[],
        metavar="SECTION.FIELD=VALUE", help="override a campaign run setting, for example ray.request_memory=8GB", )
    parser.add_argument("--no-submit", action="store_true", help="generate and validate without running the scheduler")
    return parser.parse_args()


# Compute the CPU capacity available to this local process
def available_cpu_count() -> int:
    try:
        return max(1, len(os.sched_getaffinity(0)))
    except (AttributeError, OSError):
        return max(1, os.cpu_count() or 1)


# Run the driver-independent icetune fit-configuration check locally
def run_icetune_preflight(*, environment: dict[str, str], repo_dir: pathlib.Path, tunesetup_path: pathlib.Path) -> None:
    command = runtime.icetune_command(environment, cdir=repo_dir, phase="preflight", tunesetup_path=tunesetup_path)
    runtime_environment = os.environ.copy()
    runtime_environment.update(environment)
    try:
        completed = subprocess.run(
            command, cwd=repo_dir, env=runtime_environment, check=False, capture_output=True, text=True)
    except (OSError, ValueError, subprocess.SubprocessError) as exc:
        raise RuntimeError(f"icetune preflight rejected the campaign before submission: {exc}") from exc
    if completed.returncode != 0:
        detail = "\n".join(part.strip() for part in (completed.stdout, completed.stderr) if part.strip())
        raise RuntimeError(f"icetune preflight rejected the campaign before submission:\n{detail}")
    if completed.stdout:
        print(completed.stdout, end="" if completed.stdout.endswith("\n") else "\n")
    print(f"ICETUNE preflight: valid initial configuration for {environment['TUNESETUP']}")


# Resolve campaign options without writing files or requesting allocations
def resolve_campaign(args) -> dict:
    if args.partition is not None and args.scheduler != "slurm":
        raise ValueError("--partition is used only with --scheduler slurm")
    if args.queue is not None and args.scheduler != "pbs":
        raise ValueError("--queue is used only with --scheduler pbs")
    catalog_path = args.catalog.resolve()
    catalog = campaigns.load_campaign_catalog(catalog_path)
    overrides = campaigns.parse_overrides(args.environment_overrides)
    if getattr(args, "no_ampfit_reuse", False):
        if "AMPFIT_REUSE" in overrides:
            raise ValueError("Dedicated option and --set both override: AMPFIT_REUSE")
        overrides["AMPFIT_REUSE"] = "0"
    for name in ("init_cpu", "init_gpu", "head_cpu", "head_gpu", "worker_gpu"):
        value = getattr(args, name, None)
        if value is not None:
            key = f"RAY_{name.upper()}"
            if key in overrides:
                raise ValueError("Dedicated option and --set both override: " + key)
            overrides[key] = str(value)
    workers = args.workers
    if args.scheduler in {"local", "condor"} and workers is None:
        workers = 1
    worker_cpu = args.worker_cpu
    if args.scheduler == "local":
        if worker_cpu is None and "RAY_WORKER_CPU" not in overrides:
            worker_cpu = available_cpu_count()
        local_cpus = worker_cpu if worker_cpu is not None else int(overrides["RAY_WORKER_CPU"])
        campaign = catalog["campaigns"].get(args.campaign, {}).get("ray", {})
        default_cpus = catalog["defaults"]["ray"]["cpu_per_trial"]
        requested_cpus = int(campaign.get("cpu_per_trial", default_cpus))
        if local_cpus > 0 and "CPU_PER_TRIAL" not in overrides and requested_cpus > local_cpus:
            overrides["CPU_PER_TRIAL"] = str(local_cpus)
    return campaigns.resolve(catalog, campaign_name=args.campaign, catalog_path=catalog_path, run_name=args.run_name,
        repo_dir=args.repo_dir,
        workers=workers, worker_cpu=worker_cpu, max_concurrent=args.max_concurrent, num_trials=args.num_trials,
        rand_trials=args.rand_trials, rngseed=args.rngseed, environment_overrides=overrides, )


# Resolve Conda while the submitting user's environment directories are available
def lxplus_conda_environment(environment: dict[str, str]) -> dict[str, str]:
    conda_exe = os.environ.get("CONDA_EXE") or shutil.which("conda")
    if not conda_exe or not os.access(conda_exe, os.X_OK):
        raise ValueError("lxplus requires Conda on PATH or CONDA_EXE pointing to its executable")
    requested = environment.get("ICETUNE_CONDA_ENV")
    if not requested:
        raise ValueError("Set runtime.conda_env in the campaign steering")
    result = subprocess.run(
        [conda_exe, "shell.posix+json", "activate", requested], capture_output=True, text=True, check=False)
    if result.returncode != 0:
        raise ValueError(f"Cannot resolve Conda environment {requested!r}:\n{result.stdout}{result.stderr}")
    exports = json.loads(result.stdout)["vars"]["export"]
    prefix = pathlib.Path(exports.get("CONDA_PREFIX", os.environ.get("CONDA_PREFIX", ""))).absolute()
    if not (prefix / "conda-meta").is_dir() or not os.access(prefix / "bin/python", os.X_OK):
        raise ValueError(f"Invalid Conda environment prefix: {prefix}")
    return {"CONDA_EXE": campaigns.environment_value(str(pathlib.Path(conda_exe).resolve())),
            "ICETUNE_CONDA_PREFIX": campaigns.environment_value(str(prefix))}


# Validate inputs and write a complete scheduler submission before executing it
def prepare_submission(environment: dict, args) -> dict:
    catalog_path = args.catalog.resolve()
    if args.scheduler == "lxplus":
        schedulers.validate_lxplus_environment(environment)
    repo_dir = args.repo_dir.resolve()
    if not repo_dir.is_dir():
        raise ValueError(f"Simulator directory does not exist: {repo_dir}")
    if any(character.isspace() for character in str(repo_dir)):
        raise ValueError("Ray repository path cannot contain whitespace")
    if repo_dir.is_relative_to("/eos") and args.scheduler != "lxplus":
        raise ValueError("Ray trial IO requires a non-EOS repository path")
    lxplus_work_dir = None
    if args.scheduler == "lxplus":
        if args.storage_path is not None:
            raise ValueError(f"lxplus Ray output locations follow {campaigns.shared_output_dir_key(environment)}")
        shared_root = campaigns.shared_output_dir(environment)
        lxplus_work_dir = (args.lxplus_work_dir or shared_root / "runs/icetune" / environment["RUN_NAME"]
        ).resolve()
    storage_path = (args.storage_path or (campaigns.shared_output_dir(environment) / "runs" / "icetune"
            if lxplus_work_dir is not None else repo_dir / "runs" / "icetune")).resolve()
    if any(character.isspace() for character in str(storage_path)):
        raise ValueError("Ray storage path cannot contain whitespace")
    if storage_path.is_relative_to("/eos") and args.scheduler != "lxplus":
        raise ValueError("Ray Tune storage requires a non-EOS filesystem")
    environment = dict(environment)
    if args.scheduler == "lxplus":
        environment.update(lxplus_conda_environment(environment))
    environment["RAY_REPO_DIR"] = str(repo_dir)
    environment["RAY_STORAGE_PATH"] = str(storage_path)
    environment["RAY_UPLOAD_RUNTIME"] = str(
        int(bool(getattr(args, "upload_runtime", False) or args.scheduler == "lxplus")))
    if lxplus_work_dir is not None:
        environment["RAY_LXPLUS_WORK_DIR"] = str(lxplus_work_dir)
    if args.address is not None and args.scheduler != "external":
        raise ValueError(f"--address is not used with --scheduler {args.scheduler}")
    if args.scheduler in {"local", "condor"}:
        if int(environment["RAY_WORKERS"]) != 1:
            raise ValueError(f"{args.scheduler} Ray mode requires RAY_WORKERS=1")
        environment["RAY_ADDRESS"] = "local"
    elif args.scheduler == "external":
        address = (args.address or "").strip()
        if not address or address.lower() == "local":
            raise ValueError("--scheduler external requires a non-local --address")
        environment["RAY_ADDRESS"] = address
    elif args.scheduler == "lxplus":
        environment["RAY_ADDRESS"] = "lxplus-bootstrap"

    output_dir = args.output_dir
    if output_dir is None:
        output_dir = (lxplus_work_dir / "submit" if lxplus_work_dir is not None
            else repo_dir / "runs" / "icetune" / environment["RUN_NAME"] / "submit")
    output_dir = output_dir.resolve()
    if any(character.isspace() for character in str(output_dir)):
        raise ValueError("Ray output path cannot contain whitespace")
    ensure_dir(output_dir)
    environment["RAY_SUBMIT_DIR"] = str(output_dir)
    tunesetup_path = output_dir / f"{environment['RUN_NAME']}.tunesetup.json"
    print(f"Random seed:       {environment['RNGSEED']} ({environment['ALGORITHM']})", flush=True)
    if environment["ALGORITHM"].lower() == "ampfit":
        print("ampfit bank reuse: " + ("enabled" if environment.get("AMPFIT_REUSE", "1") == "1" else "disabled"), flush=True)
    run_icetune_preflight(environment=dict(environment), repo_dir=repo_dir, tunesetup_path=tunesetup_path)
    tuning_definition = json.loads(tunesetup_path.read_text(encoding="utf-8"))
    environment["TUNESETUP"] = "tmp/icetune/tunesetup.json" if args.scheduler == "lxplus" else str(tunesetup_path)
    for algorithm in tuning_definition["optimizer"]:
        if tuning_definition["optimizer"][algorithm] is not None:
            environment[algorithm.upper() + "_CONFIG"] = environment["TUNESETUP"] + "#/optimizer/" + algorithm
    job_environment = dict(environment)
    schedd = runtime_archive = runtime_sha256 = campaign_fingerprint = None
    if args.scheduler == "lxplus":
        schedd, schedd_address = schedulers.discover_cern_schedd()
        job_environment["RAY_LXPLUS_SCHEDD"] = schedd
        job_environment["RAY_LXPLUS_SCHEDD_ADDR"] = schedd_address
        runtime_archive, runtime_sha256 = runtime.prepare_lxplus_runtime(
            source=repo_dir, output_dir=output_dir, run_name=environment["RUN_NAME"], tunesetup_path=tunesetup_path)
        campaign_fingerprint = runtime.campaign_state_fingerprint(environment, runtime_sha256)
        job_environment["RAY_CAMPAIGN_FINGERPRINT"] = campaign_fingerprint
        job_environment["RAY_INIT_RUNTIME_ARCHIVE"] = str(runtime_archive)
        job_environment["RAY_INIT_RUNTIME_SHA256"] = runtime_sha256
        if environment.get("PRECOMPUTE", "1") == "1":
            job_environment["RAY_INIT_STATE_RELATIVE"] = (
                pathlib.PurePosixPath("runs") / "icetune" / environment["RUN_NAME"] / "ray_init.json").as_posix()
    payload = {"schema_version": 1, "catalog": str(catalog_path), "campaign": args.campaign,
        "environment": dict(environment), "tuning_definition": tuning_definition, "partition": args.partition,
        "queue": args.queue, "scheduler": args.scheduler, }
    if campaign_fingerprint is not None:
        payload["campaign_fingerprint"] = campaign_fingerprint
        payload["runtime_sha256"] = runtime_sha256
    payload["fingerprint"] = runtime.resolved_fingerprint(payload)
    resolved_revision = f".{campaign_fingerprint[:16]}" if campaign_fingerprint is not None else ""
    resolved_path = output_dir / f"{environment['RUN_NAME']}{resolved_revision}.resolved.yml"
    runtime.validate_resolved_snapshot(resolved_path, payload)
    resolved_path.write_text(yaml.safe_dump(payload, sort_keys=True), encoding="utf-8")
    if args.scheduler == "lxplus":
        if runtime_archive is None or runtime_sha256 is None:
            raise RuntimeError("Ray lxplus runtime identity was not prepared")
        init_path, steer_path, dag_path = schedulers.write_lxplus_dag(environment=job_environment,
            output_dir=output_dir, runtime_archive=runtime_archive, runtime_sha256=runtime_sha256,
            bank_samples=tuning_definition.get("bank_samples"))
        job_path = dag_path
        command = ["condor_submit_dag", "-update_submit", "-Lockfile", str(dag_path.with_suffix(".lock")),
            str(dag_path), ]
    else:
        init_path = None
        steer_path = None
        job_path, command = schedulers.write_scheduler_job(scheduler=args.scheduler, environment=job_environment,
            output_dir=output_dir, partition=args.partition, queue=args.queue, )

    print(f"Resolved campaign: {resolved_path}")
    if init_path is not None and steer_path is not None:
        print(f"INIT submit:      {init_path}")
        print(f"STEER submit:     {steer_path}")
        print(f"DAG:              {job_path}")
    elif job_path is not None:
        print(f"Scheduler job:     {job_path}")
    if job_path is not None:
        print(f"Scheduler logs:    {output_dir / 'logs'}")
    if schedd is not None:
        print(f"CERN schedd:       {schedd} {job_environment['RAY_LXPLUS_SCHEDD_ADDR']}")
    print("Submit command:    " + " ".join(shlex.quote(item) for item in command))
    return {"environment": job_environment, "command": command, "resolved": resolved_path, "job": job_path}


# Execute only a prepared submission with its resolved environment
def submit(submission: dict) -> int:
    environment = {**os.environ, **submission["environment"]}
    return subprocess.run(submission["command"], check=False, env=environment).returncode


# Resolve, prepare and optionally submit one campaign through the public CLI
def main() -> int:
    args = parse_args()
    submission = prepare_submission(resolve_campaign(args), args)
    return 0 if args.no_submit else submit(submission)


if __name__ == "__main__":
    raise SystemExit(main())
