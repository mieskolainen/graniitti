# Render scheduler jobs and CERN initialization dependencies
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import pathlib
import re
import shlex
import subprocess

from core.io.files import ensure_dir

from submit import SHELL_DIR, campaign
from submit.lxplus import nodes

SUPPORTED_SCHEDULERS = {"condor", "external", "local", "lxplus", "pbs", "slurm"}
LXPLUS_MAX_RUNTIME_S = 604800
LXPLUS_RETURN_TIME_S = 900


# Validate CERN lxplus resource and shared output settings
def validate_lxplus_environment(environment: dict[str, str]) -> None:
    shared_root = campaign.shared_output_dir(environment)
    if not shared_root.is_absolute():
        raise ValueError(f"lxplus {campaign.shared_output_dir_key(environment)} must be an absolute POSIX path")
    worker_cpu = int(environment["RAY_WORKER_CPU"])
    if nodes.ray_host_slots(worker_cpu) <= 0:
        raise ValueError("Ray lxplus CPU request exceeds the CERN port range")
    if int(environment["RAY_INIT_MAX_RUNTIME_S"]) > LXPLUS_MAX_RUNTIME_S:
        raise ValueError(f"Ray lxplus INIT maximum runtime is {LXPLUS_MAX_RUNTIME_S} seconds")
    if int(environment["RAY_MAX_RUNTIME_S"]) > LXPLUS_MAX_RUNTIME_S:
        raise ValueError(f"Ray lxplus maximum runtime is {LXPLUS_MAX_RUNTIME_S} seconds")
    lxplus_steering_runtime(environment)


# Compute the steering lifetime needed to start slots and return outputs
def lxplus_steering_runtime(environment: dict[str, str]) -> int:
    startup = int(environment.get("RAY_LXPLUS_STARTUP_TIMEOUT_S", environment["RAY_STARTUP_TIMEOUT_S"]))
    runtime = int(environment["RAY_MAX_RUNTIME_S"]) + startup + LXPLUS_RETURN_TIME_S
    if runtime > LXPLUS_MAX_RUNTIME_S:
        raise ValueError(f"Ray lxplus runtime, startup and output return budget exceeds {LXPLUS_MAX_RUNTIME_S} seconds")
    return runtime


# Compute the current CERN schedd name and full HTCondor address
def discover_cern_schedd() -> tuple[str, str]:
    query = subprocess.run(["condor_q", "-totals"], check=False, capture_output=True, text=True)
    output = "\n".join(part for part in (query.stdout, query.stderr) if part)
    if query.returncode != 0:
        raise RuntimeError(f"Cannot query the current CERN schedd: {output.strip()}")
    match = re.search(r"(?m)^--\s+Schedd:\s+([A-Za-z0-9_.-]+)\s+:", output)
    if match is None:
        raise RuntimeError(f"Cannot read the current CERN schedd name from: {output.strip()}")
    schedd = match.group(1)
    status = subprocess.run(
        ["condor_status", "-schedd", "-constraint", f'Name == "{schedd}"', "-af", "MyAddress"], check=False,
        capture_output=True, text=True, )
    detail = "\n".join(part for part in (status.stdout, status.stderr) if part)
    if status.returncode != 0:
        raise RuntimeError(f"Cannot query the CERN schedd address: {detail.strip()}")
    addresses = [line.strip() for line in status.stdout.splitlines() if line.strip()]
    if len(addresses) != 1 or re.fullmatch(r"<[^<>\s]+>", addresses[0]) is None:
        raise RuntimeError(f"Cannot identify one CERN schedd address from: {detail.strip()}")
    return schedd, addresses[0]


# Compute a scheduler-safe memory spelling
def scheduler_memory(value: str, *, scheduler: str) -> str:
    match = re.fullmatch(r"([1-9][0-9]*)(KB|MB|GB|TB)", value.upper())
    if match is None:
        raise ValueError(f"Invalid RAY_REQUEST_MEMORY={value!r}")
    amount, unit = match.groups()
    if scheduler == "slurm":
        return amount + {"KB": "K", "MB": "M", "GB": "G", "TB": "T"}[unit]
    if scheduler == "pbs":
        return amount + unit.lower()
    return amount + unit


# Convert seconds to a scheduler wall-time value
def scheduler_walltime(seconds: int, *, scheduler: str) -> str:
    if seconds <= 0:
        raise ValueError("Ray maximum runtime must be positive")
    days, remainder = divmod(seconds, 86400)
    hours, remainder = divmod(remainder, 3600)
    minutes, secs = divmod(remainder, 60)
    if scheduler == "slurm" and days:
        return f"{days}-{hours:02d}:{minutes:02d}:{secs:02d}"
    return f"{days * 24 + hours:02d}:{minutes:02d}:{secs:02d}"


# Render reproducible shell exports for one resolved job
def shell_exports(environment: dict[str, str]) -> str:
    return "\n".join(f"export {key}={shlex.quote(str(environment[key]))}" for key in sorted(environment))


# Render one HTCondor allocation containing local Ray and the fitting process
def render_condor_job(environment: dict[str, str], *, log_dir: pathlib.Path) -> str:
    if int(environment["RAY_WORKERS"]) != 1:
        raise ValueError(
            "HTCondor Ray mode uses a single allocation. Use lxplus, PBS or Slurm for multiple allocations")
    gpu_line = ""
    if int(environment["RAY_WORKER_GPU"]) > 0:
        gpu_line = f"request_gpus   = {environment['RAY_WORKER_GPU']}\n"
    rendered_environment = " ".join(f"{key}={environment[key]}" for key in sorted(environment))
    return f'''# Generated by python -m submit
universe       = vanilla
initialdir     = {environment["RAY_REPO_DIR"]}
executable     = {SHELL_DIR}/run_ray.sh
log            = {log_dir}/ICETUNE_ray_$(Cluster).log
output         = {log_dir}/ICETUNE_ray_$(Cluster)_$(Process).out
error          = {log_dir}/ICETUNE_ray_$(Cluster)_$(Process).err
request_cpus   = {environment["RAY_WORKER_CPU"]}
{gpu_line}request_memory = {scheduler_memory(environment["RAY_REQUEST_MEMORY"], scheduler="condor")}
+MaxRuntime    = {environment["RAY_MAX_RUNTIME_S"]}
environment    = "{rendered_environment}"
should_transfer_files = NO
queue 1
'''


# Render one submitted CERN steering process which acquires the Ray slots
def render_lxplus_job(environment: dict[str, str], *, log_dir: pathlib.Path) -> str:
    rendered_environment = " ".join(f"{key}={environment[key]}" for key in sorted(environment))
    steering_runtime = lxplus_steering_runtime(environment)
    return f'''# Generated by python -m submit
universe       = vanilla
initialdir     = {environment["RAY_REPO_DIR"]}
executable     = {SHELL_DIR}/steer_ray_lxplus.sh
log            = {log_dir}/ICETUNE_ray_steer_$(Cluster).log
output         = {log_dir}/ICETUNE_ray_steer_$(Cluster)_$(Process).out
error          = {log_dir}/ICETUNE_ray_steer_$(Cluster)_$(Process).err
request_cpus   = 1
request_memory = {scheduler_memory(environment["RAY_STEER_REQUEST_MEMORY"], scheduler="condor")}
request_disk   = {scheduler_memory(environment["RAY_REQUEST_DISK"], scheduler="condor")}
+MaxRuntime    = {steering_runtime}
MY.WantOS      = "el9"
MY.SendCredential = True
MY.IsDaskWorker = True
+JobBatchName  = "icetune-ray-{environment["RUN_NAME"]}-steer"
environment    = "{rendered_environment}"
want_graceful_removal = True
# Signal 2 lets STEER return state before HTCondor removal
kill_sig       = 2
job_max_vacate_time = {environment["RAY_LXPLUS_SHUTDOWN_TIMEOUT_S"]}
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_output_files = ""
queue 1
'''


# Compute the shared locations used by the serialized Ray initialization
def ray_init_locations(*, shared_root: pathlib.Path, run_name: str, campaign_fingerprint: str
) -> tuple[pathlib.Path, pathlib.Path]:
    if re.fullmatch(r"[0-9a-f]{64}", campaign_fingerprint) is None:
        raise ValueError("Ray campaign fingerprint must be lowercase hexadecimal")
    cache = shared_root / "runs" / "icetune" / "bootstrap_cache"
    coord = shared_root / "runs" / "icetune" / run_name / "ray" / "init" / campaign_fingerprint[:16]
    ensure_dir(cache)
    ensure_dir(coord)
    return cache, coord


# Compute the common INIT and BANK job environment used by the Ray lxplus DAG
def _ray_init_environment(*, environment: dict[str, str], runtime_archive: pathlib.Path, runtime_sha256: str,
    bootstrap_cache: pathlib.Path, coord_dir: pathlib.Path) -> str:
    common_keys = ("CONDA_EXE", "ICETUNE_CONDA_PREFIX", "ALGORITHM", "COST", "COST_AVG", "RAY_INIT_CPU", "RAY_INIT_GPU", "RAY_INIT_JOBS", "DATA_COVARIANCE_MODE",
        "FULL_OUTPUT", "SAVE_EVENTS", "ICEBO_CONFIG", "HEBO_CONFIG", "AMPFIT_CONFIG", "LBFGS_CONFIG", "AMPFIT_REUSE", "MAX_CONCURRENT", "MAX_T", "NUM_TRIALS",
        "PLOT", "PROPOSAL_BATCH_FRACTION", "PROPOSAL_BATCH_SIZE", "RAND_TRIALS", "RNGSEED", "RUN_NAME",
        "SIMDRIVER", "TUNESETUP", "TUNE_DEFAULT", )
    driver = campaign.simdriver_name(environment)
    driver_keys = {
        key for key in environment if key == campaign.shared_output_dir_key(environment) or key.startswith(f"{driver}_")}
    init_environment = {
        key: str(environment[key]) for key in (*common_keys, *sorted(driver_keys)) if key in environment}
    init_environment.update({"RAY_INIT_BOOTSTRAP_CACHE": str(bootstrap_cache),
            "RAY_INIT_COORD_DIR": str(coord_dir), "RAY_INIT_RUNTIME_ARCHIVE": runtime_archive.name,
            "RAY_INIT_RUNTIME_SHA256": runtime_sha256, })
    for key, value in init_environment.items():
        if re.search(r'[\s"\\]', value):
            raise ValueError(f"Unsafe Ray INIT environment value {key}={value!r}")
    return " ".join(f"{key}={value}" for key, value in sorted(init_environment.items()))


# Render an initialization or per-sample preparation job
def render_ray_init_job(*, environment: dict[str, str], runtime_archive: pathlib.Path, runtime_sha256: str,
    init_script: pathlib.Path, log_dir: pathlib.Path, bootstrap_cache: pathlib.Path, coord_dir: pathlib.Path,
    phase="init", arguments="") -> str:
    if phase == "bank":
        environment = environment | {"RAY_INIT_CPU": "$(BankCpu)"}
    rendered_environment = _ray_init_environment(environment=environment, runtime_archive=runtime_archive,
        runtime_sha256=runtime_sha256, bootstrap_cache=bootstrap_cache, coord_dir=coord_dir)
    arguments = f"arguments      = {arguments}\n" if arguments else ""
    return f'''# Generated by python -m submit
universe       = vanilla
initialdir     = {runtime_archive.parent}
executable     = {init_script}
{arguments}log            = {log_dir}/ICETUNE_ray_{phase}_$(Cluster).log
output         = {log_dir}/ICETUNE_ray_{phase}_$(Cluster)_$(Process).out
error          = {log_dir}/ICETUNE_ray_{phase}_$(Cluster)_$(Process).err
request_cpus   = {environment["RAY_INIT_CPU"]}
request_gpus   = {environment["RAY_INIT_GPU"]}
request_memory = {environment["RAY_INIT_REQUEST_MEMORY"]}
request_disk   = {environment["RAY_INIT_REQUEST_DISK"]}
+MaxRuntime    = {environment["RAY_INIT_MAX_RUNTIME_S"]}
MY.WantOS      = "el9"
MY.SendCredential = True
+JobBatchName  = "icetune-ray-{environment["RUN_NAME"]}-{phase}"
environment    = "{rendered_environment}"
should_transfer_files = YES
when_to_transfer_output = ON_EXIT
transfer_input_files = {runtime_archive}
transfer_output_files = ""
queue 1
'''


# Divide the preparation allocations according to event and component counts
def bank_allocations(samples, jobs, cpus):
    counts = [1] * len(samples)
    limits = [max(1, (sample["components"] + cpus - 1) // cpus) for sample in samples]
    for _ in range(max(0, jobs - len(samples))):
        available = [index for index, count in enumerate(counts) if count < limits[index]]
        if not available:
            break
        index = max(available, key=lambda index: samples[index]["nevents"] * samples[index]["components"] / counts[index])
        counts[index] += 1
    return counts


# Quote one DAGMan value
def dag_quote(value: object) -> str:
    return str(value).replace("\\", "\\\\").replace('"', '\\"')


# Write the Ray lxplus DAG with optional INIT before STEER
def write_lxplus_dag(
    *, environment: dict[str, str], output_dir: pathlib.Path, runtime_archive: pathlib.Path, runtime_sha256: str,
    bank_samples: list[dict] | None = None
) -> tuple[pathlib.Path | None, pathlib.Path, pathlib.Path]:
    run_name = environment["RUN_NAME"]
    if environment.get("RAY_INIT_RUNTIME_SHA256") != runtime_sha256:
        raise ValueError("Ray lxplus DAG runtime identity is inconsistent")
    log_dir = output_dir / "logs"
    ensure_dir(log_dir)
    shared_root = campaign.shared_output_dir(environment)
    bootstrap_cache, coord_dir = ray_init_locations(
        shared_root=shared_root, run_name=run_name, campaign_fingerprint=environment["RAY_CAMPAIGN_FINGERPRINT"])
    environment["RAY_INIT_COORD_DIR"] = str(coord_dir)
    revision = environment["RAY_CAMPAIGN_FINGERPRINT"][:16]
    init_path = output_dir / f"{run_name}.{revision}.init.condor"
    steer_path = output_dir / f"{run_name}.{revision}.steer.condor"
    dag_path = output_dir / f"{run_name}.{revision}.dag"
    if environment.get("PRECOMPUTE", "1") == "0":
        steer_path.write_text(render_lxplus_job(environment, log_dir=log_dir), encoding="utf-8")
        dag_path.write_text(
            f"JOB STEER {dag_quote(steer_path)}\nRETRY STEER 0\nABORT-DAG-ON STEER 64 RETURN 1\n", encoding="utf-8")
        return None, steer_path, dag_path
    init_path.write_text(render_ray_init_job(environment=environment, runtime_archive=runtime_archive,
            runtime_sha256=runtime_sha256, init_script=SHELL_DIR / "steer_ray_init_lxplus.sh", log_dir=log_dir,
            bootstrap_cache=bootstrap_cache, coord_dir=coord_dir, ), encoding="utf-8", )
    steer_path.write_text(render_lxplus_job(environment, log_dir=log_dir), encoding="utf-8")
    dag_lines = [f"JOB INIT {dag_quote(init_path)}", f"JOB STEER {dag_quote(steer_path)}"]
    jobs, cpus = int(environment["RAY_INIT_JOBS"]), int(environment["RAY_INIT_CPU"])
    if jobs > 1 and environment.get("ALGORITHM", "").lower() == "ampfit":
        if not bank_samples:
            raise ValueError("Distributed ampfit requires the sample plan from icetune preflight")
        common = dict(environment=environment, runtime_archive=runtime_archive, runtime_sha256=runtime_sha256,
                      log_dir=log_dir, bootstrap_cache=bootstrap_cache, coord_dir=coord_dir)
        sample_path = output_dir / f"{run_name}.{revision}.sample.condor"
        bank_path = output_dir / f"{run_name}.{revision}.bank.condor"
        finalize_path = output_dir / f"{run_name}.{revision}.finalize.condor"
        sample_path.write_text(render_ray_init_job(**common, phase="sample",
            init_script=SHELL_DIR / "steer_ray_sample_lxplus.sh", arguments="$(BankIndex) $(Shards)"), encoding="utf-8")
        bank_path.write_text(render_ray_init_job(**common, phase="bank",
            init_script=SHELL_DIR / "steer_ray_bank_lxplus.sh", arguments="$(BankIndex) $(Shard) $(Shards)"), encoding="utf-8")
        finalize_path.write_text(render_ray_init_job(**common, phase="finalize",
            init_script=SHELL_DIR / "steer_ray_bank_lxplus.sh", arguments="$(BankIndex) finalize $(Shards)"), encoding="utf-8")
        from core.io.readers import clean_filename

        for index, shards in enumerate(bank_allocations(bank_samples, jobs, cpus)):
            entry = bank_samples[index]
            directory = (shared_root / "runs/icetune" / run_name / "results/amplitude"
                         / str(entry["dataset"]) / clean_filename(entry["sample"]))
            if environment.get("AMPFIT_REUSE", "1") == "1" and (directory / "finalized.pkl").is_file():
                continue
            sample = f"SAMPLE_{index}"
            dag_lines += [f"JOB {sample} {dag_quote(sample_path)}",
                          f'VARS {sample} BankIndex="{index}" Shards="{shards}"',
                          f"RETRY {sample} 2 UNLESS-EXIT 64", f"ABORT-DAG-ON {sample} 64 RETURN 1",
                          f"CATEGORY {sample} PREP"]
            banks = []
            for shard in range(shards):
                bank = f"BANK_{index}_{shard}"
                bank_cpu = min(cpus, max(1, (bank_samples[index]["components"] + shards - 1 - shard) // shards))
                banks.append(bank)
                dag_lines += [f"JOB {bank} {dag_quote(bank_path)}",
                              f'VARS {bank} BankIndex="{index}" Shard="{shard}" Shards="{shards}" BankCpu="{bank_cpu}"',
                              f"RETRY {bank} 2 UNLESS-EXIT 64", f"ABORT-DAG-ON {bank} 64 RETURN 1",
                              f"CATEGORY {bank} PREP"]
            dag_lines.append(f"PARENT {sample} CHILD {' '.join(banks)}")
            if shards > 1:
                finalize = f"FINALIZE_{index}"
                dag_lines += [f"JOB {finalize} {dag_quote(finalize_path)}",
                              f'VARS {finalize} BankIndex="{index}" Shards="{shards}"',
                              f"RETRY {finalize} 2 UNLESS-EXIT 64", f"ABORT-DAG-ON {finalize} 64 RETURN 1",
                              f"CATEGORY {finalize} PREP", f"PARENT {' '.join(banks)} CHILD {finalize}",
                              f"PARENT {finalize} CHILD INIT"]
            else:
                dag_lines.append(f"PARENT {banks[0]} CHILD INIT")
        dag_lines.append(f"MAXJOBS PREP {jobs}")
    dag_lines += ["PARENT INIT CHILD STEER", "RETRY INIT 2 UNLESS-EXIT 64", "ABORT-DAG-ON INIT 64 RETURN 1",
                # A STEER retry submits a new head and the entire worker array
                "RETRY STEER 0", "ABORT-DAG-ON STEER 64 RETURN 1", ]
    dag_path.write_text("\n".join(dag_lines) + "\n", encoding="utf-8", )
    return init_path, steer_path, dag_path


# Render one PBS allocation that starts a Ray cluster internally
def render_pbs_job(environment: dict[str, str], *, log_dir: pathlib.Path, queue: str | None) -> str:
    if queue is not None and re.fullmatch(r"[A-Za-z0-9_.-]+", queue) is None:
        raise ValueError(f"Unsafe PBS queue {queue!r}")
    queue_line = f"#PBS -q {queue}\n" if queue else ""
    select = (f"1:ncpus={environment['RAY_HEAD_CPU']}:ngpus={environment['RAY_HEAD_GPU']}"
        f":mem={scheduler_memory(environment['RAY_HEAD_REQUEST_MEMORY'], scheduler='pbs')}"
        f"+{environment['RAY_WORKERS']}:ncpus={environment['RAY_WORKER_CPU']}" f":ngpus={environment['RAY_WORKER_GPU']}"
        f":mem={scheduler_memory(environment['RAY_REQUEST_MEMORY'], scheduler='pbs')}")
    return f"""#!/bin/bash
# Generated by python -m submit
#PBS -N ICETUNE_{environment["RUN_NAME"]}
#PBS -l select={select}
#PBS -l place=scatter
#PBS -l walltime={scheduler_walltime(int(environment["RAY_MAX_RUNTIME_S"]), scheduler="pbs")}
#PBS -o {log_dir}/ICETUNE_ray_pbs.log
#PBS -j oe
#PBS -V
#PBS -S /bin/bash
{queue_line}
set -euo pipefail
{shell_exports(environment)}
cd "${{RAY_REPO_DIR}}"
bash {SHELL_DIR}/steer_ray_allocation.sh pbs
"""


# Render one Slurm allocation that starts a Ray cluster internally
def render_slurm_job(environment: dict[str, str], *, log_dir: pathlib.Path, partition: str | None) -> str:
    if partition is not None and re.fullmatch(r"[A-Za-z0-9_.-]+", partition) is None:
        raise ValueError(f"Unsafe Slurm partition {partition!r}")
    partition_line = f"#SBATCH --partition={partition}\n" if partition else ""
    # Reserve uniform physical hosts large enough for either head or worker allocation
    host_cpu = max(int(environment["RAY_HEAD_CPU"]), int(environment["RAY_WORKER_CPU"]))
    host_gpu = max(int(environment["RAY_HEAD_GPU"]), int(environment["RAY_WORKER_GPU"]))
    host_memory = max((environment["RAY_HEAD_REQUEST_MEMORY"], environment["RAY_REQUEST_MEMORY"]),
        key=lambda value: int(value[:-2]) * 1024 ** {"KB": 1, "MB": 2, "GB": 3, "TB": 4}[value[-2:].upper()], )
    gpu_line = ""
    if host_gpu > 0:
        gpu_line = f"#SBATCH --gpus-per-node={host_gpu}\n"
    return f"""#!/bin/bash
# Generated by python -m submit
#SBATCH --job-name=ICETUNE_{environment["RUN_NAME"]}
#SBATCH --chdir={environment["RAY_REPO_DIR"]}
#SBATCH --nodes={int(environment["RAY_WORKERS"]) + 1}
#SBATCH --ntasks-per-node=1
#SBATCH --cpus-per-task={host_cpu}
{gpu_line}#SBATCH --mem={scheduler_memory(host_memory, scheduler="slurm")}
#SBATCH --time={scheduler_walltime(int(environment["RAY_MAX_RUNTIME_S"]), scheduler="slurm")}
#SBATCH --output={log_dir}/ICETUNE_ray_slurm_%j.out
{partition_line}
set -euo pipefail
{shell_exports(environment)}
bash {SHELL_DIR}/steer_ray_allocation.sh slurm
"""


# Write the scheduler description and return its submission command
def write_scheduler_job(
    *, scheduler: str, environment: dict[str, str], output_dir: pathlib.Path, partition: str | None, queue: str | None
) -> tuple[pathlib.Path | None, list[str]]:
    log_dir = output_dir / "logs"
    ensure_dir(log_dir)
    if scheduler in {"external", "local"}:
        return None, ["bash", str(SHELL_DIR / "run_ray.sh")]
    if scheduler == "lxplus":
        raise ValueError("Ray lxplus submission requires write_lxplus_dag")
    suffix = {"condor": "condor", "pbs": "pbs", "slurm": "slurm"}[scheduler]
    job_path = output_dir / f"{environment['RUN_NAME']}.{suffix}"
    if scheduler == "condor":
        text = render_condor_job(environment, log_dir=log_dir)
        command = ["condor_submit", str(job_path)]
    elif scheduler == "pbs":
        text = render_pbs_job(environment, log_dir=log_dir, queue=queue)
        command = ["qsub", str(job_path)]
    else:
        text = render_slurm_job(environment, log_dir=log_dir, partition=partition)
        command = ["sbatch", str(job_path)]
    job_path.write_text(text, encoding="utf-8")
    job_path.chmod(0o755)
    return job_path, command
