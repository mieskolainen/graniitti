# Stage icetune runtime inputs and verify immutable run settings
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import json
import os
import pathlib
import re
import shutil
import subprocess
import sys
import tarfile
import tempfile
from importlib.metadata import PackageNotFoundError, distribution

import yaml
from core import resource
from core.io.files import ensure_dir
from core.io.serialize import sha256_file

from submit import CAMPAIGN_DIR

COPY_IGNORE = shutil.ignore_patterns("__pycache__", "*.pyc", "*.pyo", "*._old", "condor_logs")


BASE_PATHS = ("bin/gr", "bin/ampfit", "icepack", "gencard", "HEPData", "modeldata", "VERSION.json",
    "install/setconda_lxplus.sh", "install/setenv.sh", "tests/environment.sh", )


# Copy the compact static runtime into lxplus local scratch
def stage_runtime(*, source: pathlib.Path, target: pathlib.Path) -> pathlib.Path:
    source = source.resolve()
    target = target.resolve()
    if source == target or source in target.parents or target in source.parents:
        raise ValueError("Ray lxplus runtime scratch must be separate from the source checkout")
    ensure_dir(target)
    for relative in BASE_PATHS:
        src = source / relative
        if not src.exists():
            raise FileNotFoundError(f"Ray lxplus runtime input is missing: {src}")
        dst = target / relative
        ensure_dir(dst.parent)
        if src.is_dir():
            shutil.copytree(src, dst, dirs_exist_ok=True, ignore=COPY_IGNORE)
        else:
            shutil.copy2(src, dst)
    shutil.copytree(CAMPAIGN_DIR, target / "submit", dirs_exist_ok=True, ignore=COPY_IGNORE)
    shutil.copytree(resource(""), target / "python/core", dirs_exist_ok=True, ignore=COPY_IGNORE)
    try:
        dist = distribution("core")
    except PackageNotFoundError:
        return target
    for item in dist.files or ():
        if item.parts[0].endswith(".dist-info"):
            output = target / "python" / str(item)
            ensure_dir(output.parent)
            shutil.copy2(dist.locate_file(item), output)
    return target


# Extract the exact content-addressed runtime completed by the INIT DAG node
def extract_lxplus_runtime(*, archive: pathlib.Path, target: pathlib.Path, checksum: str) -> pathlib.Path:
    if re.fullmatch(r"[0-9a-f]{64}", checksum) is None:
        raise ValueError("Ray lxplus runtime SHA256 must be lowercase hexadecimal")
    if not archive.is_file():
        raise FileNotFoundError(f"Ray lxplus runtime archive is missing: {archive}")
    if sha256_file(archive) != checksum:
        raise RuntimeError(f"Ray lxplus runtime checksum mismatch: {archive}")
    ensure_dir(target)
    completed = subprocess.run(
        ["tar", "--zstd", "-xf", str(archive), "-C", str(target)], check=False, capture_output=True, text=True)
    if completed.returncode != 0:
        raise RuntimeError(f"Cannot extract Ray lxplus runtime: {(completed.stderr or completed.stdout).strip()}")
    return target


# Package one staged runtime on disk for transfer to the Ray allocations
def pack_runtime(runtime: pathlib.Path) -> pathlib.Path:
    target = runtime.with_name(runtime.name + ".tar.gz")
    with tarfile.open(target, mode="w:gz") as archive:
        archive.add(runtime, arcname="runtime", recursive=True)
    return target


# Publish an immutable archive for Ray nodes to extract without a GCS upload
def publish_ray_runtime(*, archive: pathlib.Path, directory: pathlib.Path) -> str:
    digest = sha256_file(archive)
    ensure_dir(directory)
    target = directory / f"{digest[:16]}.tar.gz"
    if target.exists():
        if sha256_file(target) != digest:
            raise ValueError(f"Ray runtime archive checksum mismatch: {target}")
    else:
        with tempfile.NamedTemporaryFile(dir=directory, prefix=".runtime-", delete=False) as handle:
            with archive.open("rb") as source:
                shutil.copyfileobj(source, handle)
            handle.flush()
            os.fsync(handle.fileno())
        pathlib.Path(handle.name).replace(target)
    return target.resolve().as_uri()


# Stage the verified simulator bootstrap into the runtime sent to the Ray head
def stage_ray_bootstrap(*, runtime: pathlib.Path, coord_dir: pathlib.Path, run_name: str, simdriver: str) -> str:
    init_path = coord_dir / "init.json"
    state = json.loads(init_path.read_text(encoding="utf-8"))
    bootstrap = state.get("bootstrap")
    sys.path.insert(0, str(runtime / "python"))
    try:
        from importlib import import_module

        from core.tune import init as icetune_init

        from submit.campaign import simdriver_name

        driver = simdriver_name({"SIMDRIVER": simdriver}).lower()
        driver_runtime = import_module(f"core.tune.drivers.{driver}.runtime")
    finally:
        sys.path.pop(0)
    if (state.get("backend") != "ray" or state.get("protocol") != icetune_init.PROTOCOL_NAME
        or state.get("protocol_version") != icetune_init.PROTOCOL_VERSION or state.get("status") != "completed"
        or not isinstance(bootstrap, dict)):
        raise RuntimeError(f"Ray INIT did not publish a valid Ray bootstrap: {init_path}")

    if bootstrap.get("kind") != driver:
        raise RuntimeError(f"Ray INIT bootstrap does not match {simdriver}")
    driver_runtime.stage_bootstrap(cdir=runtime, bootstrap=bootstrap)
    state_target = runtime / "runs" / "icetune" / run_name / "ray_init.json"
    ensure_dir(state_target.parent)
    shutil.copy2(init_path, state_target)
    return str(bootstrap["fingerprint"])


# Create one deterministic archive for the lxplus INIT and Ray allocations
def pack_lxplus_runtime(*, runtime: pathlib.Path, output: pathlib.Path) -> str:
    ensure_dir(output.parent)
    process_environment = os.environ.copy()
    process_environment["ZSTD_CLEVEL"] = "1"
    completed = subprocess.run(["tar", "--zstd", "--sort=name", "--mtime=@0", "--owner=0", "--group=0",
            "--numeric-owner", "--mode=u+rwX,go+rX,go-w", "-cf", str(output), "-C", str(runtime), ".", ],
        check=False, capture_output=True, env=process_environment, text=True, )
    if completed.returncode != 0:
        raise RuntimeError(f"Cannot pack Ray lxplus runtime: {(completed.stderr or completed.stdout).strip()}")
    return sha256_file(output)


# Stage one immutable runtime archive before submitting the lxplus DAG
def prepare_lxplus_runtime(*, source: pathlib.Path, output_dir: pathlib.Path, run_name: str, tunesetup_path: pathlib.Path
) -> tuple[pathlib.Path, str]:
    with tempfile.TemporaryDirectory(prefix=f"icetune-ray-{run_name}-") as stage:
        stage_path = pathlib.Path(stage)
        runtime = stage_runtime(source=source, target=stage_path / "runtime")
        definition = runtime / "tmp" / "icetune" / "tunesetup.json"
        ensure_dir(definition.parent)
        shutil.copy2(tunesetup_path, definition)
        payload = json.loads(tunesetup_path.read_text())
        for relative, source_path in payload["aux_param_space"].get("runtime_files", {}).items():
            target = (runtime / relative).resolve()
            if not target.is_relative_to(runtime.resolve()):
                raise ValueError(f"Runtime input is outside the staged repository: {relative}")
            ensure_dir(target.parent)
            shutil.copy2(source_path, target)
        temporary = stage_path / "runtime.tar.zst"
        checksum = pack_lxplus_runtime(runtime=runtime, output=temporary)
        archive = output_dir / f"{run_name}.runtime-{checksum[:16]}.tar.zst"
        if archive.is_file():
            if sha256_file(archive) != checksum:
                raise RuntimeError(f"Conflicting Ray lxplus runtime archive: {archive}")
        else:
            ensure_dir(archive.parent)
            shutil.copy2(temporary, archive)
    return archive, checksum


# Compute the stable identity of one generated Ray submission
def resolved_fingerprint(payload: dict) -> str:
    identity = {key: payload.get(key) for key in ("campaign", "campaign_fingerprint", "environment", "partition",
            "queue", "runtime_sha256", "scheduler", "tuning_definition", )}
    encoded = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


# Compute the immutable state identity before requesting batch allocations
def campaign_state_fingerprint(environment: dict[str, str], runtime_sha256: str) -> str:
    if re.fullmatch(r"[0-9a-f]{64}", runtime_sha256) is None:
        raise ValueError("Ray lxplus runtime SHA256 must be lowercase hexadecimal")
    execution_keys = {"ICETUNE_LCG_VIEW", "CONDA_EXE", "ICETUNE_CONDA_PREFIX"}
    settings = {key: value for key, value in environment.items()
        if not key.startswith("RAY_") and not key.endswith("_SHARED_OUTPUT_DIR") and key not in execution_keys}
    encoded = json.dumps(
        {"runtime_sha256": runtime_sha256, "settings": settings}, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


# Refuse to reuse a Ray output directory with incompatible settings
def validate_resolved_snapshot(path: pathlib.Path, payload: dict) -> None:
    if not path.is_file():
        return
    existing = yaml.safe_load(path.read_text(encoding="utf-8"))
    if not isinstance(existing, dict) or existing.get("fingerprint") != resolved_fingerprint(existing):
        raise ValueError(f"Malformed resolved Ray campaign snapshot: {path}")
    if existing["fingerprint"] != payload["fingerprint"]:
        raise ValueError(f"Ray output directory contains another campaign: {path}; use --run-name")


# Build the same icetune settings for preflight, BANK, INIT and fitting
def icetune_command(environment, *, cdir, phase="fit", tunesetup_path=None, bank_shard=None, bank_index=None, bank_jobs=None):
    env = environment
    command = [sys.executable, "-m", "core.icetune"]
    if phase == "preflight":
        command += ["--preflight", "--save_tunesetup", str(tunesetup_path)]
    elif phase not in {"sample", "bank", "finalize", "init", "fit"}:
        raise ValueError(f"Unknown icetune phase: {phase}")
    command += ["--backend", "ray", "--cdir", str(cdir)]
    for name in ("simdriver", "tunesetup", "tune_default", "run_name", "num_trials", "rand_trials", "rngseed", "algorithm"):
        command += ["--" + name, env[name.upper()]]
    command += ["--precompute", "false" if env.get("PRECOMPUTE", "1") == "0" else "true"]
    if env.get("SAVE_EVENTS", "0") == "1":
        command += ["--save_events"]
        if phase == "fit" and env.get("RAY_CAMPAIGN_FINGERPRINT"):
            from submit.campaign import shared_output_dir
            from submit.lxplus.outputs import campaign_run_path

            command += ["--events_dir", str(campaign_run_path(shared_output_dir(env), env["RUN_NAME"],
                env["RAY_CAMPAIGN_FINGERPRINT"]) / "events")]
    if env["SIMDRIVER"] == "GRANIITTI":
        command += ["--data_covariance_mode", env["DATA_COVARIANCE_MODE"]]
    algorithm = env["ALGORITHM"].lower()
    if algorithm in {"hebo", "icebo", "ampfit", "lbfgs"}:
        command += [f"--{algorithm}_config", env[f"{algorithm.upper()}_CONFIG"]]
    if algorithm == "ampfit":
        command += ["--ampfit_reuse", "true" if env.get("AMPFIT_REUSE", "1") == "1" else "false"]
        if phase in {"sample", "bank", "finalize", "init"} or (phase == "preflight" and env.get("RAY_LXPLUS_WORK_DIR")):
            from submit.campaign import shared_output_dir

            command += ["--bank_shared_dir", str(shared_output_dir(env))]
    if phase == "preflight":
        return command
    command += ["--restore", "--cost_rho", "quadratic", "--max_concurrent_trials", env["MAX_CONCURRENT"]]
    for name in ("proposal_batch_size", "proposal_batch_fraction", "max_t", "cost", "cost_avg"):
        command += ["--" + name, env[name.upper()]]
    command += ["--plot", "true" if env.get("PLOT", "1") == "1" else "false"]
    if env.get("FULL_OUTPUT", "1") == "1":
        command += ["--pickle_dump"]
    if phase in {"sample", "bank", "finalize", "init"}:
        command += ["--phase", "init" if phase == "init" else "bank", "--cpu_per_trial", env["RAY_INIT_CPU"],
                    "--coord_dir", env["RAY_INIT_COORD_DIR"], "--bootstrap_cache_url", env["RAY_INIT_BOOTSTRAP_CACHE"],
                    "--runtime_sha256", env["RAY_INIT_RUNTIME_SHA256"]]
        from submit.campaign import shared_output_dir

        if phase in {"sample", "bank", "finalize"}:
            command += ["--bank_jobs", str(bank_jobs if bank_jobs is not None else env["RAY_INIT_JOBS"])]
            if bank_index is not None:
                command += ["--bank_index", str(bank_index)]
            if phase in {"sample", "finalize"}:
                command += ["--bank_" + phase]
            else:
                if bank_shard is None:
                    raise ValueError("Ray BANK phase requires a shard index")
                command += ["--bank_shard", str(bank_shard)]
        if phase == "init":
            if algorithm == "ampfit" and int(env.get("RAY_INIT_JOBS", "1")) > 1:
                command += ["--bank_jobs", env["RAY_INIT_JOBS"]]
            if env.get("FULL_OUTPUT", "1") == "1":
                command += ["--pickle_dir", str(shared_output_dir(env) / "runs/icetune" / env["RUN_NAME"] / "results")]
        return command
    address = env.get("RAY_ADDRESS", "local")
    command += ["--address", address, "--cpu_per_trial", env["CPU_PER_TRIAL"], "--gpu_per_trial", env["GPU_PER_TRIAL"],
                "--ray_storage_path", env.get("RAY_STORAGE_PATH", str(pathlib.Path(cdir) / "runs/icetune")),
                "--history_interval_s", env.get("RAY_HISTORY_INTERVAL_S", "60")]
    for name in ("worker_admission_timeout_s", "worker_poll_interval_s", "worker_failure_limit", "worker_recovery_timeout_s",
                 "restore_retry_limit", "trial_retry_limit", "restore_retry_delay_s", "resource_status_interval_s", "status_interval_s"):
        command += ["--ray_" + name, env["RAY_" + name.upper()]]
    if env.get("RAY_UPLOAD_RUNTIME", "0") == "1":
        command += ["--ray_upload_runtime"]
        if env.get("RAY_RUNTIME_URI"):
            command += ["--ray_runtime_uri", env["RAY_RUNTIME_URI"]]
    if env.get("RAY_INIT_STATE_RELATIVE"):
        command += ["--ray_init_state", str(pathlib.Path(cdir) / env["RAY_INIT_STATE_RELATIVE"]),
                    "--runtime_sha256", env["RAY_INIT_RUNTIME_SHA256"]]
    if address.lower() == "local":
        command += ["--ray_worker_cpu", env["RAY_WORKER_CPU"], "--ray_worker_gpu", env["RAY_WORKER_GPU"],
                    "--ray_temp_dir", env["RAY_TMP_DIR"]]
    else:
        command += ["--wait_workers", env.get("RAY_WAIT_WORKERS", env["RAY_WORKERS"]),
                    "--wait_workers_timeout_s", env.get("RAY_STARTUP_TIMEOUT_S", "10800")]
    return command


# Replace the launcher with icetune so scheduler signals reach it directly
if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(description="Launch resolved icetune settings")
    parser.add_argument("--phase", choices=("sample", "bank", "finalize", "init", "fit"), default="fit")
    parser.add_argument("--cdir", required=True)
    for name in ("bank_shard", "bank_index", "bank_jobs"):
        parser.add_argument("--" + name, type=int, default=None)
    args = parser.parse_args()
    command = icetune_command(os.environ, cdir=args.cdir, phase=args.phase, bank_shard=args.bank_shard,
                              bank_index=args.bank_index, bank_jobs=args.bank_jobs)
    os.execv(command[0], command)
