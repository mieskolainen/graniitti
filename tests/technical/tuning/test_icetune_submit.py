# Tests for icetune campaign resolution and submission
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
import argparse
import asyncio
import io
import json
import re
import shutil
import socket
import subprocess
import sys
import tarfile
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from threading import Barrier
from types import SimpleNamespace

import pytest
import yaml
from core.io.serialize import sha256_file
from core.tune.core import ICETUNE_HISTORY_FIGURES

from submit import CAMPAIGN_DIR, campaign_source, schedulers
from submit import __main__ as submit
from submit import runtime as run_runtime
from submit.lxplus import cluster as ray_cluster
from submit.lxplus import common as run_common
from submit.lxplus import main as run_main
from submit.lxplus import nodes as ray_nodes
from submit.lxplus import outputs as run_outputs

ROOT = Path(__file__).resolve().parents[3]
ICETUNE_TEST_DIR = CAMPAIGN_DIR

from core.tune.drivers.graniitti.driver import GraniittiDriver  # noqa: E402
from core.tune.drivers.pandora.driver import PandoraDriver  # noqa: E402
from core.tune.tunesetup import load_tunesetup, resolve_tunesetup, save_tunesetup  # noqa: E402

from submit import campaign as campaign_config  # noqa: E402

CAMPAIGN_FP = "c" * 64


# Resolve shared paths at submission for both simulator drivers
@pytest.mark.parametrize("campaign", ["tune-gpom-res-con", "tune-pandora-v0"])
@pytest.mark.parametrize("shared", [None, "outputs", "~/icetune", "/eos/user/example/tuning"])
def test_campaign_shared_path(tmp_path, monkeypatch, campaign, shared):
    monkeypatch.chdir(tmp_path)
    catalog = campaign_config.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    for repo_dir in (None, tmp_path / "simulator"):
        environment = campaign_config.resolve(catalog, campaign_name=campaign, repo_dir=repo_dir,
            environment_overrides={} if shared is None else {"SHARED_OUTPUT_DIR": shared})
        expected = ((repo_dir or tmp_path) / Path(shared or ".").expanduser()).resolve()
        assert campaign_config.shared_output_dir(environment) == expected
    catalog["defaults"]["runtime"]["shared_output_dir"] = "default-output"
    catalog["campaigns"][campaign].setdefault("runtime", {})["shared_output_dir"] = "campaign-output"
    environment = campaign_config.resolve(catalog, campaign_name=campaign,
        environment_overrides={} if shared is None else {"SHARED_OUTPUT_DIR": shared})
    assert campaign_config.shared_output_dir(environment) == (tmp_path / Path(shared or "campaign-output").expanduser()).resolve()


# Preserve campaign storage overrides and resolve relative paths against --repo-dir
def test_submit_shared_path(submit_args):
    submit_args.repo_dir = submit_args.output_dir / "simulator"
    environment = submit.resolve_campaign(submit_args)
    assert campaign_config.shared_output_dir(environment) == submit_args.repo_dir.resolve()
    submit_args.environment_overrides = ["runtime.shared_output_dir=../shared"]
    environment = submit.resolve_campaign(submit_args)
    assert campaign_config.shared_output_dir(environment) == (submit_args.repo_dir / "../shared").resolve()


# Resolve the Pandora installation before jobs change their working directory
@pytest.mark.parametrize("pfa", [None, "../custom-pandora", "~/PandoraPFA", "/opt/PandoraPFA"])
def test_campaign_pfa_path(tmp_path, monkeypatch, pfa):
    monkeypatch.chdir(tmp_path)
    catalog = campaign_config.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    catalog["campaigns"]["tune-pandora-v0"]["runtime"].pop("pfa_dir")
    for repo_dir in (None, tmp_path / "simulator"):
        environment = campaign_config.resolve(catalog, campaign_name="tune-pandora-v0", repo_dir=repo_dir,
            environment_overrides={} if pfa is None else {"PANDORA_PFA_DIR": pfa})
        expected = ((repo_dir or tmp_path) / Path(pfa or "../PandoraPFA").expanduser()).resolve()
        assert Path(environment["PANDORA_PFA_DIR"]) == expected


# Resolve source tunesetup without confusing saved runtime copies with new definitions
def test_named_tunesetup_ignores_saved_runtimes(tmp_path):
    source = tmp_path / "submit/tunecards/graniitti/tune-gpom-xpom-ampfit.json"
    saved = tmp_path / "runs/icetune/example/runtime/submit/tunecards/graniitti/tune-gpom-xpom-ampfit.json"
    for path in (source, saved):
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text("{}", encoding="utf-8")
    assert resolve_tunesetup(cdir=tmp_path, simdriver="GRANIITTI", name=str(source)) == source
    assert resolve_tunesetup(cdir=tmp_path, simdriver="GRANIITTI", name=str(saved)) == saved


# Check Ray activates the resolved prefix or steering environment before GRANIITTI setup
def test_ray_runtime_environment_order():
    runtime = ICETUNE_TEST_DIR / "shell" / "ray_runtime.sh"
    script = r'''source "$1"
activate_conda_environment() { echo "activate=$1"; }
source() { echo "source=$1"; }
python() { return 0; }
run() { (unset ICETUNE_CONDA_PREFIX ICETUNE_CONDA_ENV ICETUNE_ENV_SETUP; eval "$1"; prepare_icetune_ray_environment "$2"; echo "raylet_wait=${RAY_raylet_start_wait_time_s}"); }
run "export ICETUNE_CONDA_PREFIX=/opt/conda/envs/custom ICETUNE_CONDA_ENV=custom" "$2"
run "export ICETUNE_CONDA_ENV=custom" "$2"
run "export ICETUNE_ENV_SETUP=custom-setup" "$2"'''
    completed = subprocess.run(
        ["bash", "-c", script, "bash", str(runtime), str(ROOT)],
        cwd=ROOT,
        check=True,
        capture_output=True,
        text=True,
    )
    assert completed.stdout.splitlines() == [
        "activate=/opt/conda/envs/custom",
        "source=install/setenv.sh",
        "raylet_wait=600",
        "activate=custom",
        "source=install/setenv.sh",
        "raylet_wait=600",
        "source=custom-setup",
        "source=install/setenv.sh",
        "raylet_wait=600",
    ]


# Check Ray resource overrides recompute or explicitly fix concurrency
def test_ray_resource_overrides_control_concurrency(monkeypatch):
    monkeypatch.setattr(sys, "argv", ["submit_ray.py", "--campaign", "tune-gpom-res-con",
                                     "--workers", "3", "--worker-cpu", "32"])
    args = submit.parse_args()
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    derived = campaign_config.resolve(
        catalog,
        campaign_name=args.campaign,
        workers=args.workers,
        worker_cpu=args.worker_cpu,
    )
    assert derived["MAX_CONCURRENT"] == "12"
    schedulers.validate_lxplus_environment(derived)
    explicit = campaign_config.resolve(
        catalog,
        campaign_name="tune-gpom-res-con",
        workers=3,
        worker_cpu=32,
        max_concurrent=7,
    )
    assert explicit["MAX_CONCURRENT"] == "7"
    gpu = campaign_config.resolve(
        catalog,
        campaign_name="tune-gpom-res-con",
        workers=3,
        worker_cpu=32,
        environment_overrides={"GPU_PER_TRIAL": "1", "RAY_WORKER_GPU": "2"},
    )
    assert gpu["MAX_CONCURRENT"] == "6"
    uneven = campaign_config.resolve(
        catalog,
        campaign_name="tune-gpom-res-con",
        workers=2,
        worker_cpu=10,
        environment_overrides={"CPU_PER_TRIAL": "6"},
    )
    assert uneven["MAX_CONCURRENT"] == "2"
    with pytest.raises(ValueError, match="exceeds RAY_WORKER_CPU"):
        campaign_config.resolve(
            catalog,
            campaign_name="tune-gpom-res-con",
            workers=2,
            worker_cpu=8,
            environment_overrides={"CPU_PER_TRIAL": "12"},
        )
    with pytest.raises(ValueError, match="exceeds RAY_WORKER_GPU"):
        campaign_config.resolve(
            catalog,
            campaign_name="tune-gpom-res-con",
            workers=2,
            environment_overrides={"GPU_PER_TRIAL": "3", "RAY_WORKER_GPU": "2"},
        )
    with pytest.raises(ValueError, match="available Ray trial slots"):
        campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con", workers=1, max_concurrent=9)
    with pytest.raises(ValueError, match="Unsupported ray optimizer"):
        campaign_config.resolve(
            catalog, campaign_name="tune-gpom-res-con", environment_overrides={"ALGORITHM": "broken"}
        )
    without_icebo = dict(derived)
    without_icebo.pop("ICEBO_CONFIG")
    without_icebo["ALGORITHM"] = "hebo"
    campaign_config.validate_environment(without_icebo)
    without_icebo["ALGORITHM"] = "icebo"
    with pytest.raises(ValueError, match="ICEBO_CONFIG"):
        campaign_config.validate_environment(without_icebo,)
    without_hebo = dict(derived)
    without_hebo.pop("HEBO_CONFIG")
    without_hebo["ALGORITHM"] = "hebo"
    with pytest.raises(ValueError, match="HEBO_CONFIG"):
        campaign_config.validate_environment(without_hebo,)
    with pytest.raises(ValueError, match="Set both RAY_MIN_WORKER_PORT"):
        campaign_config.resolve(
            catalog,
            campaign_name="tune-gpom-res-con",
            environment_overrides={"RAY_MIN_WORKER_PORT": "12000"},
        )
    allocated = campaign_config.resolve(
        catalog,
        campaign_name="tune-gpom-res-con",
        environment_overrides={"RAY_WORKER_GPU": "1"},
    )
    assert allocated["GPU_PER_TRIAL"] == "0"


# Check dotted overrides select run controls and reject physics switches
def test_campaign_field_overrides():
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    overrides = campaign_config.parse_overrides([
        "optimizer.num_trials=24", "optimizer.rand_trials=6", "runtime.conda_env=custom",
        "ray.worker_admission_timeout_s=900",
        "ray.trial_retry_limit=0",
        "ray.workers=3", "ray.worker_cpu=8", "ray.cpu_per_trial=2",
    ])
    environment = campaign_config.resolve(
        catalog, campaign_name="tune-gpom-res-con", environment_overrides=overrides,
    )
    assert environment["NUM_TRIALS"] == "24"
    assert environment["RAND_TRIALS"] == "6"
    assert environment["ICETUNE_CONDA_ENV"] == "custom"
    assert environment["RAY_WORKER_ADMISSION_TIMEOUT_S"] == "900"
    assert environment["RAY_TRIAL_RETRY_LIMIT"] == "0"
    assert environment["RAY_WORKERS"] == "3"
    assert environment["RAY_WORKER_CPU"] == "8"
    assert environment["CPU_PER_TRIAL"] == "2"
    assert environment["MAX_CONCURRENT"] == "12"
    with pytest.raises(ValueError, match="RAY_WORKER_ADMISSION_TIMEOUT_S"):
        campaign_config.resolve(
            catalog, campaign_name="tune-gpom-res-con",
            environment_overrides={"RAY_WORKER_ADMISSION_TIMEOUT_S": "0"},
        )
    with pytest.raises(ValueError, match="Unknown campaign field"):
        campaign_config.parse_overrides(["core.analysis.FF_prod=1"])
    with pytest.raises(ValueError, match="Duplicate"):
        campaign_config.parse_overrides(["optimizer.num_trials=24", "optimizer.num_trials=25"])
    with pytest.raises(ValueError, match="Dedicated option and --set"):
        campaign_config.resolve(
            catalog, campaign_name="tune-gpom-res-con", num_trials=30, environment_overrides=overrides,
        )


# Check lxplus allows several trials per allocation and rejects invalid initialization and port settings
@pytest.mark.parametrize("campaign", ["tune-gpom-res-con", "tune-pandora-v0"])
def test_ray_lxplus_resource_limits(campaign):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name=campaign, workers=1)
    schedulers.validate_lxplus_environment(environment)
    with pytest.raises(ValueError, match="exceeds the CERN port range"):
        schedulers.validate_lxplus_environment(dict(environment, RAY_WORKER_CPU="84"))
    packed = campaign_config.resolve(catalog, campaign_name=campaign, workers=2, worker_cpu=8,
        environment_overrides={"CPU_PER_TRIAL": "1", "RAY_INIT_CPU": "1"})
    schedulers.validate_lxplus_environment(packed)
    assert packed["MAX_CONCURRENT"] == "16"
    for init_cpu in ("1", "4"):
        independent = campaign_config.resolve(
            catalog, campaign_name=campaign,
            environment_overrides={"CPU_PER_TRIAL": "2", "RAY_INIT_CPU": init_cpu},
        )
        schedulers.validate_lxplus_environment(independent)
        assert independent["CPU_PER_TRIAL"] == "2"
        assert independent["RAY_INIT_CPU"] == init_cpu
    with pytest.raises(ValueError, match="604800 seconds"):
        schedulers.validate_lxplus_environment(dict(environment, RAY_MAX_RUNTIME_S="604801"))
    with pytest.raises(ValueError, match="runtime, startup and output return budget"):
        schedulers.validate_lxplus_environment(dict(environment, RAY_MAX_RUNTIME_S="603601"))
    for memory in ("0GB", "16", "invalid"):
        with pytest.raises(ValueError, match="RAY_STEER_REQUEST_MEMORY"):
            campaign_config.resolve(catalog, campaign_name=campaign,
                environment_overrides=campaign_config.parse_overrides([f"ray.steer_request_memory={memory}"]))


# Check lxplus captures the schedd address used by the submission host
def test_ray_lxplus_discovers_submit_point(monkeypatch):
    address = "<188.185.149.150:9618?addrs=188.185.149.150-9618&alias=bigbird103.cern.ch&noUDP>"
    calls = []

    # Compute the abbreviated queue banner followed by the full collector address
    def fake_run(command, **kwargs):
        calls.append(command)
        if command[0] == "condor_q":
            return SimpleNamespace(
                returncode=0,
                stdout=(
                    "-- Schedd: bigbird103.cern.ch : "
                    "<188.185.149.150:9618?... @ 08/26/26 00:00:00\n"
                ),
                stderr="",
            )
        return SimpleNamespace(returncode=0, stdout=f"{address}\n", stderr="")

    monkeypatch.setattr(submit.subprocess, "run", fake_run)

    assert schedulers.discover_cern_schedd() == ("bigbird103.cern.ch", address)
    assert calls == [
        ["condor_q", "-totals"],
        [
            "condor_status",
            "-schedd",
            "-constraint",
            'Name == "bigbird103.cern.ch"',
            "-af",
            "MyAddress",
        ],
    ]


# Check generated scheduler descriptions carry the complete resolved budget
@pytest.mark.parametrize("steer_memory", ["7GB", "9000mb"])
def test_ray_scheduler_rendering_resolved_campaign(tmp_path, steer_memory):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con", workers=1,
        environment_overrides=campaign_config.parse_overrides([f"ray.steer_request_memory={steer_memory}"]))
    environment["RAY_REPO_DIR"] = str(ROOT)
    environment["RAY_ADDRESS"] = "local"

    condor = schedulers.render_condor_job(environment, log_dir=tmp_path)
    assert "queue 1" in condor
    assert "request_cpus   = 8" in condor
    assert f"NUM_TRIALS={environment['NUM_TRIALS']}" in condor
    assert f"RAND_TRIALS={environment['RAND_TRIALS']}" in condor
    assert "RAY_ADDRESS=local" in condor

    lxplus_job = schedulers.render_lxplus_job(dict(environment, RAY_WORKERS="3"), log_dir=tmp_path)
    assert "steer_ray_lxplus.sh" in lxplus_job
    assert "request_cpus   = 1" in lxplus_job
    assert f"request_memory = {steer_memory.upper()}" in lxplus_job
    assert "MY.SendCredential = True" in lxplus_job
    assert "MY.IsDaskWorker = True" in lxplus_job
    assert "want_graceful_removal = True" in lxplus_job
    assert "kill_sig       = 2" in lxplus_job
    assert "job_max_vacate_time = 120" in lxplus_job
    assert "on_exit_hold" not in lxplus_job
    assert "stream_output" not in lxplus_job
    assert "stream_error" not in lxplus_job
    assert f'JobBatchName  = "icetune-ray-{environment["RUN_NAME"]}-steer"' in lxplus_job

    multi = dict(environment, RAY_WORKERS="2", RAY_REQUEST_MEMORY="4GB")
    with pytest.raises(ValueError, match="single allocation"):
        schedulers.render_condor_job(multi, log_dir=tmp_path)

    multi.pop("RAY_ADDRESS")
    slurm = schedulers.render_slurm_job(multi, log_dir=tmp_path, partition="compute")
    assert "#SBATCH --nodes=3" in slurm
    assert "#SBATCH --partition=compute" in slurm
    assert f"export NUM_TRIALS={environment['NUM_TRIALS']}" in slurm
    pbs = schedulers.render_pbs_job(multi, log_dir=tmp_path, queue="hx")
    assert "#PBS -l select=1:ncpus=8:ngpus=0:mem=8gb+2:ncpus=8:ngpus=0:mem=4gb" in pbs
    assert "#PBS -q hx" in pbs
    with pytest.raises(ValueError, match="Unsafe Slurm partition"):
        schedulers.render_slurm_job(multi, log_dir=tmp_path, partition="compute\n#SBATCH --nodes=9")
    with pytest.raises(ValueError, match="Unsafe PBS queue"):
        schedulers.render_pbs_job(multi, log_dir=tmp_path, queue="hx\n#PBS -V")


# Check a failed icetune preflight exposes its diagnostic without submitting
def test_preflight_reports_invalid_initial_config(monkeypatch, tmp_path):
    calls = []

    def fake_run(command, **kwargs):
        calls.append((command, kwargs))
        return SimpleNamespace(
            returncode=64,
            stdout="",
            stderr=(
                'Initial fit configuration is invalid: Parameter "coupling" received '
                "out-of-bounds value 2.5; expected [0.0, 2.0]\n"
            ),
        )

    monkeypatch.setattr(submit.subprocess, "run", fake_run)
    environment = {
        "RUN_NAME": "invalid",
        "NUM_TRIALS": "24",
        "RAND_TRIALS": "6",
        "RNGSEED": "1234",
        "DATA_COVARIANCE_MODE": "diagonal",
        "SIMDRIVER": "GRANIITTI",
        "TUNESETUP": "tune-gpom-res-con",
        "TUNE_DEFAULT": "TUNE0",
        "ALGORITHM": "hebo",
        "HEBO_CONFIG": "python/src/core/tune/settings/hebo.json",
    }

    with pytest.raises(RuntimeError, match=r'"coupling".*expected \[0\.0, 2\.0\]'):
        submit.run_icetune_preflight(environment=environment, repo_dir=ROOT, tunesetup_path=tmp_path / "tunesetup.json")

    command, kwargs = calls[0]
    assert command[1:4] == ["-m", "core.icetune", "--preflight"]
    assert command[command.index("--hebo_config") + 1] == environment["HEBO_CONFIG"]
    assert kwargs["check"] is False
    assert kwargs["capture_output"] is True


# Supply independent submission arguments for each scheduler test
@pytest.fixture
def submit_args(tmp_path, monkeypatch):
    args = SimpleNamespace(
        campaign="tune-gpom-res-con",
        scheduler="slurm",
        catalog=ICETUNE_TEST_DIR / "campaigns.yml",
        run_name="ray-test",
        workers=2,
        worker_cpu=16,
        max_concurrent=None,
        num_trials=20,
        rand_trials=5,
        rngseed=None,
        repo_dir=ROOT,
        output_dir=tmp_path,
        storage_path=tmp_path / "storage",
        address=None,
        partition="compute",
        queue=None,
        environment_overrides=[],
        no_submit=True,
    )
    run = subprocess.run

    # Execute real preflight and archive commands and reject scheduler submission
    def run_local(command, **kwargs):
        assert command[0] in (sys.executable, "tar", "/bin/bash") or Path(command[0]).name == "conda"
        return run(command, **kwargs)

    monkeypatch.setattr(submit.subprocess, "run", run_local)
    monkeypatch.setattr(submit, "parse_args", lambda: args)
    return args


# Check Ray dry submission records state without invoking the scheduler
@pytest.mark.parametrize("algorithm", ["hebo", "lbfgs"])
def test_ray_no_submit_writes_resolved_state(tmp_path, monkeypatch, submit_args, algorithm):
    submit_args.environment_overrides = [f"optimizer.algorithm={algorithm}"]
    preflights = []
    preflight = submit.run_icetune_preflight

    # Record the real saved tuning definition used by submission
    def record_preflight(**kwargs):
        preflights.append(kwargs)
        return preflight(**kwargs)

    monkeypatch.setattr(submit, "run_icetune_preflight", record_preflight)

    assert submit.main() == 0
    resolved = yaml.safe_load((tmp_path / "ray-test.resolved.yml").read_text(encoding="utf-8"))
    assert resolved["environment"]["NUM_TRIALS"] == "20"
    assert resolved["environment"]["RAND_TRIALS"] == "5"
    assert resolved["environment"]["MAX_CONCURRENT"] == "4"
    assert (tmp_path / "ray-test.slurm").is_file()
    assert preflights[0]["repo_dir"] == ROOT
    assert Path(preflights[0]["environment"]["TUNESETUP"]).stem == "tune-gpom-res-con"

    # Reuse an unchanged run and reject a saved continuum definition that differs from the source
    assert submit.main() == 0
    tunesetup_path = tmp_path / "ray-test.tunesetup.json"
    saved = json.loads(tunesetup_path.read_text(encoding="utf-8"))
    from core.icetune import resolve_optimizer_settings

    algorithm = resolved["environment"]["ALGORITHM"]
    optimizer = SimpleNamespace(cdir=str(ROOT), algorithm=algorithm,
                                **{algorithm + "_config": resolved["environment"][algorithm.upper() + "_CONFIG"]})
    resolve_optimizer_settings(optimizer, argparse.ArgumentParser())
    assert getattr(optimizer, algorithm + "_settings") == saved["optimizer"][algorithm]
    key = next(
        key for key in saved["param_space"]
        if key.startswith("CON_GP|") and "[321,321]/opposite:" in key and key.endswith("@NORM")
    )
    saved["param_space"].pop(key)
    tunesetup_path.write_text(json.dumps(saved), encoding="utf-8")
    with pytest.raises(RuntimeError, match="different tuning definition"):
        submit.main()
    assert json.loads(tunesetup_path.read_text(encoding="utf-8")) == saved


# Check local submission uses the process CPU allocation and fits one trial
def test_ray_local_submission_detects_available_cpus(tmp_path, monkeypatch, submit_args):
    args = submit_args
    args.scheduler = "local"
    args.run_name = "ray-local-test"
    args.workers = None
    args.worker_cpu = None
    args.partition = None
    monkeypatch.setattr(submit, "available_cpu_count", lambda: 4)

    assert submit.main() == 0
    resolved = yaml.safe_load(
        (tmp_path / "ray-local-test.resolved.yml").read_text(encoding="utf-8")
    )
    environment = resolved["environment"]
    assert environment["RAY_ADDRESS"] == "local"
    assert environment["RAY_WORKER_CPU"] == "4"
    assert environment["CPU_PER_TRIAL"] == "4"
    assert environment["MAX_CONCURRENT"] == "1"
    assert environment["RAY_STORAGE_PATH"] == str((tmp_path / "storage").resolve())


# Check an existing portable Ray cluster can run through the same campaign helper
def test_ray_external_cluster_storage(tmp_path, monkeypatch, submit_args):
    args = submit_args
    args.scheduler = "external"
    args.run_name = "ray-external-test"
    args.workers = 3
    args.storage_path = tmp_path / "shared-storage"
    args.address = "ray-head.example:60010"
    args.partition = None

    assert submit.main() == 0
    resolved = yaml.safe_load(
        (tmp_path / "ray-external-test.resolved.yml").read_text(encoding="utf-8")
    )
    environment = resolved["environment"]
    assert environment["RAY_ADDRESS"] == "ray-head.example:60010"
    assert environment["RAY_WORKERS"] == "3"
    assert environment["MAX_CONCURRENT"] == "6"
    assert environment["RAY_STORAGE_PATH"] == str((tmp_path / "shared-storage").resolve())


# Check lxplus submission paths and portable runtime inputs
@pytest.mark.parametrize("directory", ["default", "work", "submit"])
def test_ray_lxplus_records_portable_runtime(tmp_path, monkeypatch, submit_args, directory):
    shared_output_dir = tmp_path / "graniitti"
    args = submit_args
    args.scheduler = "lxplus"
    args.run_name = "ray-lxplus-test"
    args.workers = 3
    args.worker_cpu = 8
    args.worker_gpu = 2
    args.init_cpu = 2
    args.init_gpu = 1
    args.head_cpu = 4
    args.head_gpu = 3
    args.output_dir = None
    args.storage_path = None
    args.upload_runtime = False
    args.lxplus_work_dir = tmp_path / "work" if directory == "work" else None
    args.output_dir = tmp_path / "submission" if directory == "submit" else None
    args.partition = None
    args.environment_overrides = [f"runtime.shared_output_dir={shared_output_dir}"]
    monkeypatch.setattr(
        schedulers,
        "discover_cern_schedd",
        lambda: ("bigbird17.cern.ch", "<188.185.149.150:9618?noUDP>"),
    )

    # Package the actual campaign and tunesetup inputs through the production archive path
    monkeypatch.setattr(run_runtime, "BASE_PATHS", ("submit/campaigns.yml",))

    assert submit.main() == 0
    work_dir = args.lxplus_work_dir or shared_output_dir / "runs/icetune/ray-lxplus-test"
    output_dir = args.output_dir or work_dir / "submit"
    resolved_path = next(output_dir.glob("ray-lxplus-test.*.resolved.yml"))
    resolved = yaml.safe_load(resolved_path.read_text(encoding="utf-8"))
    environment = resolved["environment"]
    assert environment["RAY_ADDRESS"] == "lxplus-bootstrap"
    for role in ("init", "head", "worker"):
        for field in ("cpu", "gpu"):
            name = f"{role}_{field}"
            assert environment[f"RAY_{name.upper()}"] == str(getattr(args, name))
    assert environment["ICETUNE_LCG_VIEW"].startswith("/cvmfs/sft.cern.ch/lcg/views/LCG_110/")
    assert environment["RAY_WORKERS"] == "3"
    assert environment["RAY_UPLOAD_RUNTIME"] == "1"
    assert environment["GRANIITTI_SHARED_OUTPUT_DIR"] == str(shared_output_dir)
    assert environment["RAY_LXPLUS_WORK_DIR"] == str(work_dir.resolve())
    assert environment["RAY_SUBMIT_DIR"] == str(output_dir.resolve())
    assert environment["RAY_STORAGE_PATH"] == str((shared_output_dir / "runs" / "icetune").resolve())
    assert "RAY_LXPLUS_SCHEDD" not in environment
    assert "RAY_LXPLUS_SCHEDD_ADDR" not in environment
    archive = next(output_dir.glob("*.tar.zst"))
    assert resolved["campaign_fingerprint"] == run_runtime.campaign_state_fingerprint(
        environment,
        sha256_file(archive),
    )
    assert resolved_path.name == (
        f"ray-lxplus-test.{resolved['campaign_fingerprint'][:16]}.resolved.yml"
    )
    revision = resolved["campaign_fingerprint"][:16]
    job = output_dir / f"ray-lxplus-test.{revision}.steer.condor"
    assert job.is_file()
    job_text = job.read_text(encoding="utf-8")
    assert "steer_ray_lxplus.sh" in job_text
    assert f"RAY_SUBMIT_DIR={output_dir}" in job_text
    assert f"{output_dir}/logs/ICETUNE_ray_steer_" in job_text
    assert "RAY_LXPLUS_SCHEDD=bigbird17.cern.ch" in job_text
    assert "RAY_LXPLUS_SCHEDD_ADDR=<188.185.149.150:9618?noUDP>" in job_text
    assert "RAY_LXPLUS_STARTUP_TIMEOUT_S" not in job_text
    assert "RAY_CAMPAIGN_FINGERPRINT=" in job_text
    coord = shared_output_dir / "runs/icetune/ray-lxplus-test/ray/init" / revision
    assert f"RAY_INIT_COORD_DIR={coord}" in job_text
    assert coord.is_dir()
    assert "RAY_HEAD_CPU=4" in job_text and "RAY_HEAD_GPU=3" in job_text
    assert "RAY_WORKER_CPU=8" in job_text and "RAY_WORKER_GPU=2" in job_text
    assert "request_disk   = 40GB" in job_text
    assert "on_exit_hold" not in job_text
    init_job = output_dir / f"ray-lxplus-test.{revision}.init.condor"
    dag = output_dir / f"ray-lxplus-test.{revision}.dag"
    assert init_job.is_file()
    init_text = init_job.read_text(encoding="utf-8")
    assert Path(environment["CONDA_EXE"]).is_file()
    assert f"CONDA_EXE={environment['CONDA_EXE']}" in init_text
    assert f"CONDA_EXE={environment['CONDA_EXE']}" in job_text
    assert Path(environment["ICETUNE_CONDA_PREFIX"]).samefile(sys.prefix)
    assert f"ICETUNE_CONDA_PREFIX={environment['ICETUNE_CONDA_PREFIX']}" in init_text
    assert f"ICETUNE_CONDA_PREFIX={environment['ICETUNE_CONDA_PREFIX']}" in job_text
    assert "request_cpus   = 2" in init_text
    assert "request_gpus   = 1" in init_text
    assert dag.is_file()
    assert "PARENT INIT CHILD STEER" in dag.read_text(encoding="utf-8")


# Reject missing Pandora simulation inputs before constructing any batch jobs
def test_lxplus_rejects_missing_pandora_inputs(tmp_path, monkeypatch, submit_args):
    shared_output_dir = tmp_path / "pandora"
    args = submit_args
    args.campaign = "tune-pandora-v0"
    args.scheduler = "lxplus"
    args.run_name = "pandora-lxplus-test"
    args.workers = 3
    args.worker_cpu = None
    args.output_dir = None
    args.storage_path = None
    args.upload_runtime = False
    args.lxplus_work_dir = None
    args.partition = None
    args.environment_overrides = [
            f"runtime.shared_output_dir={shared_output_dir}",
            f"runtime.pfa_dir={ROOT.parent / 'PandoraPFA'}",
        ]
    monkeypatch.setattr(
        schedulers,
        "discover_cern_schedd",
        lambda: ("bigbird17.cern.ch", "<188.185.149.150:9618?noUDP>"),
    )

    monkeypatch.setenv("PANDORA_INPUT_DIR", str(tmp_path / "missing_inputs"))
    with pytest.raises(RuntimeError, match="preflight rejected.*",):
        submit.main()
    assert not list(shared_output_dir.rglob("*.dag"))
    assert not list(shared_output_dir.rglob("*.condor"))


# Check a submitted steering job recovers its original CERN schedd
def test_ray_lxplus_resolves_submit_point(tmp_path, monkeypatch):
    job_ad = tmp_path / "job.ad"
    job_ad.write_text(
        'ClusterId = 123\nGlobalJobId = "bigbird17.cern.ch#123.0#1770000000"\n',
        encoding="utf-8",
    )
    monkeypatch.setenv("_CONDOR_JOB_AD", str(job_ad))
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "scratch"))
    monkeypatch.delenv("RAY_LXPLUS_SCHEDD", raising=False)
    assert ray_nodes.condor_schedd() == "bigbird17.cern.ch"
    assert ray_nodes.condor_cluster_id() == 123
    ray_tmp = ray_nodes.ray_cluster_tmp("bigbird17.cern.ch", 123)
    assert ray_tmp.parent == Path("/tmp")
    assert re.fullmatch(r"irt-[0-9a-f]{12}", ray_tmp.name)
    monkeypatch.setenv("RAY_LXPLUS_SCHEDD", "bigbird18.cern.ch")
    assert ray_nodes.condor_schedd() == "bigbird18.cern.ch"
    address = "<188.185.149.150:9618?noUDP>"
    monkeypatch.setenv("RAY_LXPLUS_SCHEDD_ADDR", address)
    assert ray_nodes.condor_target() == ("bigbird18.cern.ch", ["-addr", address])
    monkeypatch.setenv("RAY_LXPLUS_SCHEDD", "bad name")
    with pytest.raises(ValueError, match="Invalid HTCondor schedd"):
        ray_nodes.condor_schedd()
    monkeypatch.setenv("RAY_LXPLUS_SCHEDD", "bigbird18.cern.ch")
    monkeypatch.setenv("RAY_LXPLUS_SCHEDD_ADDR", "bad address")
    with pytest.raises(ValueError, match="Invalid HTCondor schedd address"):
        ray_nodes.condor_target()


# Check Ray rejects precreated links at its common host-local session path
def test_lxplus_rejects_linked_session_tmp(tmp_path):
    target = tmp_path / "target"
    target.mkdir()
    linked = tmp_path / "irt-0123456789ab"
    linked.symlink_to(target, target_is_directory=True)

    with pytest.raises(RuntimeError, match="private owned directory"):
        ray_nodes.prepare_ray_tmp(linked)


# Check worker placement uses exact validated HTCondor Machine values
def test_lxplus_builds_machine_exclusions(tmp_path, monkeypatch):
    machine_ad = tmp_path / "machine.ad"
    machine_ad.write_text(
        'Arch = "X86_64"\nMachine = "lxbatch0123.cern.ch"\n',
        encoding="utf-8",
    )
    monkeypatch.setenv("_CONDOR_MACHINE_AD", str(machine_ad))

    assert ray_nodes.condor_machine() == "lxbatch0123.cern.ch"
    assert ray_nodes.machine_requirements(
        ["steer.cern.ch", "lxbatch0123.cern.ch", "steer.cern.ch"]
    ) == ('(TARGET.Machine =!= "steer.cern.ch") && (TARGET.Machine =!= "lxbatch0123.cern.ch")')
    assert ray_nodes.machine_requirements(["steer.cern.ch"]) == (
        '(TARGET.Machine =!= "steer.cern.ch")'
    )
    assert [ray_nodes.ray_host_slots(cpus) for cpus in (8, 16, 32, 64, 66)] == [2, 2, 1, 1, 1]
    assert ray_nodes.ray_host_slots(67) == 0

    machine_ad.write_text('Machine = "bad name"\n', encoding="utf-8")
    with pytest.raises(RuntimeError, match="no valid Machine"):
        ray_nodes.condor_machine()


# Check the Ray resource marker uses the exact worker Condor identity
def test_lxplus_builds_worker_resource_marker(tmp_path, monkeypatch):
    job_ad = tmp_path / "job.ad"
    job_ad.write_text("ClusterId = 245570\nProcId = 17\n", encoding="utf-8")
    monkeypatch.setenv("_CONDOR_JOB_AD", str(job_ad))
    identity = ray_nodes.condor_job_identity()
    assert identity == {
        "cluster_id": 245570,
        "job_id": "245570.17",
        "proc_id": 17,
    }
    assert ray_cluster.ray_worker_marker(identity) == "icetune_worker_245570_17"


# Check port exhaustion is reported as a normal local startup failure
def test_lxplus_port_capacity_failure(monkeypatch):
    # Raise the remote port-capacity diagnostic
    def fake_start(**kwargs):
        del kwargs
        raise RuntimeError(
            "Ray needs 23 free reserved ports in CERN range 10000:10100, but only 4 are available"
        )

    monkeypatch.setattr(ray_nodes, "start_ray_node", fake_start)
    result = ray_nodes.start_ray_worker_node(
        cpus=8,
        temp_dir="/tmp/irt-0123456789ab",
        gpus=0,
        head_address="192.0.2.2:10008",
        resource_marker="icetune_worker_245570_17",
    )

    assert result["status"] == "failed"
    assert "only 4 are available" in result["error"]


# Check each Ray worker starts once without an internal retry delay
def test_ray_lxplus_starts_worker_once(monkeypatch):
    calls = []

    # Compute one valid simulated Ray worker node
    def fake_start(**kwargs):
        assert kwargs["role"] == "worker"
        calls.append(kwargs)
        return {"hostname": "worker.cern.ch", "spill_dir": "/tmp/ray-spill"}

    monkeypatch.setattr(ray_nodes, "start_ray_node", fake_start)

    result = ray_nodes.start_ray_worker_node(
        cpus=8,
        temp_dir="/tmp/irt-0123456789ab",
        gpus=0,
        head_address="192.0.2.2:10008",
        resource_marker="icetune_worker_245570_17",
    )

    assert len(calls) == 1
    assert result["status"] == "started"


# Check slow Ray starts remain in flight while completed workers are collected
def test_lxplus_starts_workers_asynchronously():
    calls = []

    class FakeFuture:
        # Store one simulated task result or transport failure
        def __init__(self, *, value=None, error=None, complete=True):
            self.value = value
            self.error = error
            self.complete = complete
            self.result_calls = 0

        # Report whether this simulated startup has completed
        def done(self):
            return self.complete

        # Compute only this allocation's startup outcome
        def result(self):
            self.result_calls += 1
            assert self.complete
            if self.error is not None:
                raise self.error
            return self.value

    class FakeClient:
        # Submit one task pinned to exactly one Dask allocation
        def submit(self, function, **kwargs):
            worker = kwargs["workers"][0]
            calls.append((function, worker, kwargs.copy()))
            if worker == "submit-failure":
                raise RuntimeError("scheduler rejected task")
            if worker == "task-failure":
                return FakeFuture(error=RuntimeError("worker disappeared"))
            if worker == "slow":
                return FakeFuture(complete=False)
            return FakeFuture(
                value={
                    "node": {
                        "hostname": f"{worker}.cern.ch",
                        "spill_dir": f"/scratch/{worker}/ray-spill",
                    },
                    "status": "started",
                }
            )

    workers = ["healthy-a", "submit-failure", "task-failure", "slow", "healthy-b"]
    starting = {}
    immediate = ray_cluster.submit_ray_worker_starts(
        client=FakeClient(),
        workers=workers,
        markers={
            worker: f"icetune_worker_1_{index}" for index, worker in enumerate(workers, start=1)
        },
        starting=starting,
        cpus=8,
        temp_dir="/tmp/irt-0123456789ab",
        gpus=0,
        head_address="192.0.2.2:10008",
    )
    completed = ray_cluster.collect_ray_worker_starts(starting=starting)
    results = {**immediate, **completed}

    assert results["healthy-a"]["status"] == "started"
    assert results["healthy-b"]["status"] == "started"
    assert results["submit-failure"]["status"] == "failed"
    assert "submission failed" in results["submit-failure"]["error"]
    assert results["task-failure"]["status"] == "failed"
    assert "task failed" in results["task-failure"]["error"]
    assert set(starting) == {"slow"}
    assert starting["slow"].result_calls == 0
    assert [worker for _, worker, _ in calls] == workers
    for function, worker, kwargs in calls:
        assert function is ray_nodes.start_ray_worker_node
        assert kwargs["workers"] == [worker]
        assert kwargs["allow_other_workers"] is False
        assert kwargs["pure"] is False
        assert kwargs["resource_marker"].startswith("icetune_worker_1_")


# Check fractional launch deadlines and normal polling when startup is blocked
@pytest.mark.parametrize("now, pending, starting, expected", [
    (0.0, ["worker"], {}, 1.0),
    (1.0, ["worker"], {}, 0.5),
    (2.0, ["worker"], {}, 0.0),
    (2.0, [], {}, 1.0),
    (2.0, ["worker"], {"busy": object()}, 1.0),
])
def test_ray_lxplus_start_delay(monkeypatch, now, pending, starting, expected):
    monkeypatch.setattr(run_main.time, "monotonic", lambda: now)
    assert ray_cluster.ray_start_delay(
        pending=pending, starting=starting, max_starting=1,
        launch_state={"next_at": 1.5},
    ) == pytest.approx(expected)


# Check the central launcher does not exceed its simultaneous startup limit
def test_ray_lxplus_caps_worker_startups():
    class PendingFuture:
        # Keep one simulated Ray startup unresolved
        def done(self):
            return False

    class FakeClient:
        # Compute the fixed set of live Dask allocations
        def scheduler_info(self):
            return {
                "workers": {
                    "head": {},
                    "slow-a": {},
                    "slow-b": {},
                    "waiting": {},
                }
            }

        # Reject any new startup while the configured limit is occupied
        def submit(self, *args, **kwargs):
            pytest.fail("worker startup exceeded the simultaneous limit")

    pending = ["waiting"]
    starting = {"slow-a": PendingFuture(), "slow-b": PendingFuture()}
    ray_cluster.advance_ray_workers(
        client=FakeClient(),
        head_key="head",
        pending=pending,
        seen={"head"},
        joined={"head"},
        failed=set(),
        starting=starting,
        attempts={},
        errors={},
        markers={},
        launch_state={"next_at": 0.0},
        cpus=8,
        gpus=0,
        temp_dir="/tmp/irt-0123456789ab",
        head_address="192.0.2.2:10008",
        steering_machine="steer.cern.ch",
        start_interval_s=10.0,
        max_starting=2,
        retry_limit=3,
    )

    assert pending == ["waiting"]
    assert set(starting) == {"slow-a", "slow-b"}


# Check one dask-lxplus job renders and owns the complete Ray slot array
def test_lxplus_builds_fixed_condor_array():
    class FakeCernJob:
        # Build one simulated dask-lxplus worker description
        def __init__(self, *args, name=None, **kwargs):
            self.name = name
            self.worker_args = kwargs["worker_extra_args"]
            self.job_header_dict = {
                "batch_name": name,
                "JobBatchName": '"wrong-name"',
                "MY.SendCredential": "True",
                "request_memory": "4GB",
                "transfer_output_files": '""',
            }
            if kwargs.get("log_directory"):
                self.job_header_dict.update(
                    {
                        "error": "worker-$(ClusterId).$(ProcId).err",
                        "log": "worker-$(ClusterId).log",
                        "output": "worker-$(ClusterId).$(ProcId).out",
                    }
                )

        # Compute the normal single worker Queue statement
        def job_script(self):
            return getattr(self, "script", None) or (
                "Executable = /bin/sh\nrequest_memory = 4GB\n"
                "request_cpus = 8\n"
                "MY.SendCredential = True\nshould_transfer_files = YES\n"
                'when_to_transfer_output = ON_EXIT\ntransfer_output_files = ""\n\nQueue\n'
            )

        # Compute the normal single worker HTCondor identifier
        def _job_id_from_submit_output(self, output):
            assert output == "4 job(s) submitted to cluster 244912."
            return "244912.0"

    ArrayJob = ray_cluster.cern_array_job(FakeCernJob)
    job = ArrayJob(
        "tcp://scheduler:10000",
        name="CernCluster-0",
        array_size=3,
        array_name="icetune-ray-run-workers",
        head_cpu=4,
        head_gpu=0,
        head_memory="8GB",
        worker_gpu=0,
        worker_memory="4GB",
        log_directory="/eos/user/e/example/grdev/tmp/run/submit/logs",
        worker_extra_args=["--worker-port", "10000:10100"],
    )

    assert job.worker_args[-2:] == ["--worker-port", "10000:10009"]
    assert job.name == "CernCluster-0-$(ProcId)"
    assert "batch_name" not in job.job_header_dict
    assert job.job_header_dict["transfer_output_files"] == '""'
    assert job.job_header_dict["MY.SendCredential"] == "True"
    assert job.job_header_dict["should_transfer_files"] == "YES"
    assert job.job_header_dict["when_to_transfer_output"] == "ON_EXIT"
    assert job.job_header_dict["JobBatchName"] == '"icetune-ray-run-workers"'
    assert job.job_header_dict["error"] == "worker-$(ClusterId).$(ProcId).err"
    assert job.job_header_dict["log"] == "worker-$(ClusterId).log"
    assert job.job_header_dict["output"] == "worker-$(ClusterId).$(ProcId).out"
    assert "/dev/null" not in job.job_header_dict.values()
    assert job.job_script() == (
        "Executable = /bin/sh\nrequest_memory = 8GB\n"
        "request_cpus = 4\n"
        "MY.SendCredential = True\nshould_transfer_files = YES\n"
        'when_to_transfer_output = ON_EXIT\ntransfer_output_files = ""\n\n'
        "request_gpus = 0\nQueue 1\nrequest_cpus = 8\nrequest_memory = 4GB\nrequest_gpus = 0\n"
        "on_exit_remove = False\non_exit_hold = NumJobStarts >= 3\nQueue 3\n"
    )
    assert job._job_id_from_submit_output("4 job(s) submitted to cluster 244912.") == "244912"

    job.script = (
        'MY.SendCredential = True\nshould_transfer_files = NO\ntransfer_output_files = ""\nQueue\n'
    )
    with pytest.raises(RuntimeError, match="must use an isolated transfer sandbox"):
        job.job_script()
    job.script = (
        "MY.SendCredential = True\nshould_transfer_files = YES\n"
        'transfer_output_files = "scratch"\nQueue\n'
    )
    with pytest.raises(RuntimeError, match="must return no scratch files"):
        job.job_script()
    job.script = 'should_transfer_files = YES\ntransfer_output_files = ""\nQueue\n'
    with pytest.raises(RuntimeError, match="must forward the CERN credential"):
        job.job_script()


# Check repeated Dask executions launch only one driver and preserve its completed state
@pytest.mark.parametrize("returncode", [0, 7])
def test_ray_lxplus_runtime_runs_once(tmp_path, monkeypatch, returncode):
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    monkeypatch.setenv("CUDA_VISIBLE_DEVICES", "GPU-head")
    environment = {"RUN_NAME": "run", "ICETUNE_CONDA_PREFIX": sys.prefix,
                   "CUDA_VISIBLE_DEVICES": "GPU-steer"}
    source = tmp_path / "runtime"
    runner = source / "submit/shell" / "run_ray.sh"
    runner.parent.mkdir(parents=True)
    runner.write_text(
        "#!/bin/bash\n"
        'printf "%s" "$CUDA_VISIBLE_DEVICES" > gpu.txt\n'
        "printf 'trial\\n' >> trial-count\n"
        "printf completed > state.txt\n"
        f"exit {returncode}\n",
        encoding="ascii",
    )
    (source / "state.txt").write_text("initial", encoding="ascii")
    environment["RAY_RUNTIME_URI"] = run_runtime.publish_ray_runtime(
        archive=run_runtime.pack_runtime(source), directory=tmp_path / "shared")
    barrier = Barrier(2)

    # Start two concurrent deliveries of the same Dask task
    def execute():
        barrier.wait(timeout=10)
        return run_main.run_icetune_runtime(environment=environment)

    results = []
    failures = []
    with ThreadPoolExecutor(max_workers=2) as pool:
        futures = [pool.submit(execute) for _ in range(2)]
        for future in futures:
            try:
                results.append(future.result(timeout=20))
            except RuntimeError as exc:
                failures.append(str(exc))
    assert len(results) == 1
    assert results[0]["returncode"] == returncode
    assert len(failures) == 1
    assert "Refusing to repeat the campaign" in failures[0]
    assert not run_common.runtime_pid_path(environment).exists()

    # Reject a lost result replay even after the original driver has exited
    with pytest.raises(RuntimeError, match="Refusing to repeat the campaign"):
        run_main.run_icetune_runtime(environment=environment)
    runtime = run_common.runtime_work_root(environment) / "runtime"
    assert (runtime / "trial-count").read_text() == "trial\n"
    assert (runtime / "state.txt").read_text() == "completed"
    assert (runtime / "gpu.txt").read_text() == "GPU-head"


# Transfer a runtime through a real Dask task using only its shared archive location
def test_lxplus_runtime_dask(tmp_path, monkeypatch):
    distributed = pytest.importorskip("distributed")
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    source = tmp_path / "runtime"
    runner = source / "submit/shell/run_ray.sh"
    runner.parent.mkdir(parents=True)
    runner.write_text("#!/bin/bash\ncp input.bin copied.bin\n", encoding="ascii")
    data = bytes(range(256)) * 4096
    (source / "input.bin").write_bytes(data)
    environment = dict(RUN_NAME="dask", ICETUNE_CONDA_PREFIX=sys.prefix,
                       RAY_RUNTIME_URI=run_runtime.publish_ray_runtime(
                           archive=run_runtime.pack_runtime(source), directory=tmp_path / "shared"))
    with (distributed.LocalCluster(n_workers=1, threads_per_worker=1, processes=True, dashboard_address=None,
                                   local_directory=str(tmp_path / "dask"), memory_limit=0) as cluster,
          distributed.Client(cluster) as client):
        result = client.submit(run_main.run_icetune_runtime, environment=environment, pure=False).result(timeout=30)
    assert result["returncode"] == 0
    copied = run_common.runtime_work_root(environment) / "runtime/copied.bin"
    assert copied.read_bytes() == data


# Check a failed final output transfer cannot run the campaign again
def test_lxplus_failed_snapshot_no_retry(tmp_path, monkeypatch):
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    environment = {"RUN_NAME": "run", "ICETUNE_CONDA_PREFIX": sys.prefix}
    source = tmp_path / "runtime"
    runner = source / "submit/shell" / "run_ray.sh"
    runner.parent.mkdir(parents=True)
    runner.write_text("#!/bin/bash\nprintf 'trial\\n' >> trial-count\n", encoding="ascii")
    environment["RAY_RUNTIME_URI"] = run_runtime.publish_ray_runtime(
        archive=run_runtime.pack_runtime(source), directory=tmp_path / "shared")

    # Fail after the simulator work has completed
    def fail_snapshot(**kwargs):
        raise OSError("final snapshot failed")

    monkeypatch.setattr(run_outputs, "snapshot_final_outputs", fail_snapshot)
    with pytest.raises(OSError, match="final snapshot failed"):
        run_main.run_icetune_runtime(environment=environment)
    with pytest.raises(RuntimeError, match="Refusing to repeat the campaign"):
        run_main.run_icetune_runtime(environment=environment)
    runtime = run_common.runtime_work_root(environment) / "runtime"
    assert (runtime / "trial-count").read_text() == "trial\n"


# Check a graceful stop signals only the recorded head local driver process group
def test_lxplus_requests_scoped_runtime_stop(tmp_path, monkeypatch):
    environment = {"RUN_NAME": "run"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path))
    path = run_common.runtime_pid_path(environment)
    path.parent.mkdir(parents=True)
    path.write_text('{"birth": "1234", "pid": 4321}\n', encoding="utf-8")
    signals = []
    monkeypatch.setattr(run_common, "process_birth", lambda pid: "1234")
    monkeypatch.setattr(run_outputs.os, "getpgid", lambda pid: pid)
    monkeypatch.setattr(run_outputs.os, "killpg", lambda pid, value: signals.append((pid, value)))

    assert run_main.request_runtime_stop(environment=environment) == {
        "pid": 4321,
        "requested": True,
    }
    assert signals == [(4321, run_main.signal.SIGINT)]


# Check graceful cleanup waits for the driver result before releasing allocations
def test_lxplus_graceful_shutdown_orders_cleanup(tmp_path, monkeypatch):
    lifecycle = []

    class FakeFuture:
        # Keep the simulated Ray driver active until it receives the stop request
        def done(self):
            return False

        # Compute one final head result within the graceful timeout
        def result(self, *, timeout):
            lifecycle.append(("driver_result", timeout))
            return {"returncode": 130}

    class FakeClient:
        # Send one simulated scoped interrupt on the Ray head allocation
        def run(self, function, *, environment, workers):
            assert function is run_main.request_runtime_stop
            assert environment["RUN_NAME"] == "run"
            assert workers == ["head"]
            lifecycle.append("driver_stop")
            return {"head": {"pid": 4321, "requested": True}}

        # Close the simulated Dask client after Ray exits
        def close(self):
            lifecycle.append("client_close")

    class FakeCluster:
        # Release the simulated HTCondor allocation array last
        def close(self):
            lifecycle.append("cluster_close")

    monkeypatch.setattr(
        run_outputs,
        "publish_runtime_result",
        lambda **kwargs: lifecycle.append("publish_final") or 130,
    )
    run_main.graceful_ray_shutdown(
        client=FakeClient(),
        cluster=FakeCluster(),
        future=FakeFuture(),
        head_key="head",
        environment={"RUN_NAME": "run"},
        shared_root=tmp_path,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
        figure_revision=None,
        live_revision=None,
        timeout_s=120.0,
    )

    assert lifecycle == [
        "driver_stop",
        ("driver_result", 90.0),
        "publish_final",
        "client_close",
        "cluster_close",
    ]


# Check a natural-exit race still publishes the completed final checkpoint
def test_lxplus_shutdown_exited_driver(tmp_path, monkeypatch):
    lifecycle = []

    class FakeFuture:
        # Simulate a stale completion check immediately before driver exit
        def done(self):
            return False

        # Compute the completed final output after the stop request race
        def result(self, *, timeout):
            lifecycle.append(("driver_result", timeout))
            return {"checkpoint": b"checkpoint", "returncode": 0}

    class FakeClient:
        # Report that the driver exited before the interrupt arrived
        def run(self, function, *, environment, workers):
            assert function is run_main.request_runtime_stop
            lifecycle.append("driver_exited")
            return {"head": {"requested": False, "reason": "driver process already exited"}}

        # Close after publishing the completed result
        def close(self):
            lifecycle.append("client_close")

    class FakeCluster:
        # Release the allocation last
        def close(self):
            lifecycle.append("cluster_close")

    monkeypatch.setattr(
        run_outputs,
        "publish_runtime_result",
        lambda **kwargs: lifecycle.append("publish_final") or 0,
    )
    run_main.graceful_ray_shutdown(
        client=FakeClient(),
        cluster=FakeCluster(),
        future=FakeFuture(),
        head_key="head",
        environment={"RUN_NAME": "run"},
        shared_root=tmp_path,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
        figure_revision=None,
        live_revision=None,
        timeout_s=120.0,
    )

    assert lifecycle == [
        "driver_exited",
        ("driver_result", 90.0),
        "publish_final",
        "client_close",
        "cluster_close",
    ]


# Check an unconfirmed shutdown returns live files without creating a checkpoint
def test_lxplus_shutdown_timeout_skips_checkpoint(tmp_path, monkeypatch):
    lifecycle = []

    class FakeFuture:
        # Keep the simulated Ray driver active through the final wait
        def done(self):
            return False

        # Exhaust the simulated graceful driver wait
        def result(self, *, timeout):
            lifecycle.append(("driver_result", timeout))
            raise TimeoutError("driver still stopping")

    class FakeClient:
        # Record the scoped driver interrupt
        def run(self, function, *, environment, workers):
            assert function is run_main.request_runtime_stop
            lifecycle.append("driver_stop")
            return {"head": {"pid": 4321, "requested": True}}

        # Close after the last live output return
        def close(self):
            lifecycle.append("client_close")

    class FakeCluster:
        # Release the allocation after the Dask client
        def close(self):
            lifecycle.append("cluster_close")

    # Record the bounded final live output return
    def fake_live_outputs(**kwargs):
        assert 0.0 < kwargs["timeout_s"] <= 30.0
        lifecycle.append("live_outputs")

    monkeypatch.setattr(run_outputs, "return_live_outputs", fake_live_outputs)
    monkeypatch.setattr(
        run_outputs,
        "publish_runtime_result",
        lambda **kwargs: pytest.fail("an unavailable final result must not be published"),
    )
    run_main.graceful_ray_shutdown(
        client=FakeClient(),
        cluster=FakeCluster(),
        future=FakeFuture(),
        head_key="head",
        environment={"RUN_NAME": "run"},
        shared_root=tmp_path,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
        figure_revision=None,
        live_revision=None,
        timeout_s=120.0,
    )

    assert lifecycle == [
        "driver_stop",
        ("driver_result", 90.0),
        "live_outputs",
        "client_close",
        "cluster_close",
    ]


# Check transient steering failures remain eligible for DAGMan retries
def test_ray_lxplus_steering_exit_codes():
    assert run_common.steer_exit_code(KeyError("RUN_NAME")) == 64
    assert run_common.steer_exit_code(ValueError("invalid resource value")) == 64
    assert run_common.steer_exit_code(OSError("temporary scheduler failure")) == 1
    assert run_common.steer_exit_code(RuntimeError("Ray worker failed to join")) == 1


# Check Ray INIT is one DAG allocation with the configured initialization resources
@pytest.mark.parametrize("init_cpu", [1, 4])
@pytest.mark.parametrize("init_gpu", [0, 2])
def test_ray_lxplus_renders_init_node(tmp_path, init_cpu, init_gpu):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(
        catalog, campaign_name="tune-gpom-res-con",
        environment_overrides={"CPU_PER_TRIAL": "2", "RAY_INIT_CPU": str(init_cpu),
                               "RAY_INIT_GPU": str(init_gpu)},
    )
    archive = tmp_path / "ray-init-runtime.tar.zst"
    init_script = tmp_path / "runtime" / "submit/shell" / "steer_ray_init_lxplus.sh"
    rendered = schedulers.render_ray_init_job(
        environment=environment,
        runtime_archive=archive,
        runtime_sha256="a" * 64,
        init_script=init_script,
        log_dir=tmp_path / "logs",
        bootstrap_cache=tmp_path / "bootstrap_cache",
        coord_dir=tmp_path / "ray" / "init",
    )

    assert f"executable     = {init_script}" in rendered
    assert f"request_cpus   = {init_cpu}" in rendered
    assert f"RAY_INIT_CPU={init_cpu}" in rendered
    assert f"request_gpus   = {init_gpu}" in rendered
    assert "CPU_PER_TRIAL=" not in rendered
    assert "request_memory = 4GB" in rendered
    assert "request_disk   = 20GB" in rendered
    assert "+MaxRuntime    = 14400" in rendered
    assert f"TUNESETUP={environment['TUNESETUP']}" in rendered
    assert "Queue" not in rendered
    assert rendered.endswith("queue 1\n")
    assert "on_exit_hold" not in rendered


# Check DAGMan owns the INIT to STEER dependency and failure transitions
def test_lxplus_dag_orders_init_before_steer(tmp_path):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con")
    environment.update(
        {
            "GRANIITTI_SHARED_OUTPUT_DIR": str(tmp_path / "shared"),
            "RAY_INIT_RUNTIME_ARCHIVE": str(tmp_path / "runtime.tar.zst"),
            "RAY_INIT_RUNTIME_SHA256": "a" * 64,
            "RAY_CAMPAIGN_FINGERPRINT": CAMPAIGN_FP,
            "RAY_INIT_STATE_RELATIVE": "runs/icetune/CENTRAL_GP_RES_CON/ray_init.json",
            "RAY_LXPLUS_SCHEDD": "bigbird17.cern.ch",
            "RAY_LXPLUS_SCHEDD_ADDR": "<188.185.149.150:9618?noUDP>",
            "RAY_REPO_DIR": str(ROOT),
        }
    )
    archive = Path(environment["RAY_INIT_RUNTIME_ARCHIVE"])
    archive.write_bytes(b"runtime")

    init_path, steer_path, dag_path = schedulers.write_lxplus_dag(
        environment=environment,
        output_dir=tmp_path,
        runtime_archive=archive,
        runtime_sha256="a" * 64,
    )

    dag = dag_path.read_text(encoding="utf-8")
    assert f"JOB INIT {init_path}" in dag
    assert f"JOB STEER {steer_path}" in dag
    assert "PARENT INIT CHILD STEER" in dag
    assert "RETRY INIT 2 UNLESS-EXIT 64" in dag
    assert [line for line in dag.splitlines() if line.startswith("RETRY STEER ")] == ["RETRY STEER 0"]
    assert "ABORT-DAG-ON INIT 64 RETURN 1" in dag
    assert "condor_q" not in steer_path.read_text(encoding="utf-8")
    assert init_path.name == f"{environment['RUN_NAME']}.{CAMPAIGN_FP[:16]}.init.condor"
    assert steer_path.name == f"{environment['RUN_NAME']}.{CAMPAIGN_FP[:16]}.steer.condor"

    first = {path: path.read_bytes() for path in (init_path, steer_path, dag_path)}
    environment["RAY_INIT_RUNTIME_SHA256"] = "b" * 64
    environment["RAY_CAMPAIGN_FINGERPRINT"] = "d" * 64
    second = schedulers.write_lxplus_dag(
        environment=environment,
        output_dir=tmp_path,
        runtime_archive=archive,
        runtime_sha256="b" * 64,
    )
    assert not set(first).intersection(second)
    assert all(path.read_bytes() == payload for path, payload in first.items())


# Check independent sample dependencies, component coverage and preparation resource limits
@pytest.mark.parametrize("state", ["empty", "partial", "complete", "rebuild"])
@pytest.mark.parametrize("jobs,cpus,components", [(4, 2, [12, 48]), (4, 8, [1, 4]), (2, 2, [3, 6, 8])])
def test_lxplus_dag_orders_ampfit_bank_array(tmp_path, jobs, cpus, components, state):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name="tune-gpom-star-cms-ampfit")
    environment.update(RAY_INIT_JOBS=str(jobs), RAY_INIT_CPU=str(cpus))
    samples = [dict(dataset=index, sample="nominal", nevents=100 * (index + 1), components=count)
               for index, count in enumerate(components)]
    environment.update(
        {
            "GRANIITTI_SHARED_OUTPUT_DIR": str(tmp_path / "shared"),
            "RAY_INIT_RUNTIME_ARCHIVE": str(tmp_path / "runtime.tar.zst"),
            "RAY_INIT_RUNTIME_SHA256": "a" * 64,
            "RAY_CAMPAIGN_FINGERPRINT": CAMPAIGN_FP,
            "RAY_INIT_STATE_RELATIVE": "runs/icetune/GP_AMPFIT/ray_init.json",
            "RAY_LXPLUS_SCHEDD": "bigbird17.cern.ch",
            "RAY_LXPLUS_SCHEDD_ADDR": "<188.185.149.150:9618?noUDP>",
            "RAY_REPO_DIR": str(ROOT),
        }
    )
    archive = Path(environment["RAY_INIT_RUNTIME_ARCHIVE"])
    archive.write_bytes(b"runtime")

    completed = [] if state == "empty" else samples[:1] if state == "partial" else samples
    environment["AMPFIT_REUSE"] = "0" if state == "rebuild" else "1"
    for sample in completed:
        directory = (tmp_path / "shared/runs/icetune" / environment["RUN_NAME"] / "results/amplitude"
                     / str(sample["dataset"]) / sample["sample"])
        directory.mkdir(parents=True)
        (directory / "finalized.pkl").touch()

    init_path, steer_path, dag_path = schedulers.write_lxplus_dag(
        environment=environment,
        output_dir=tmp_path,
        runtime_archive=archive,
        runtime_sha256="a" * 64,
        bank_samples=samples,
    )

    dag = dag_path.read_text(encoding="utf-8")
    if state == "complete":
        assert [line.split()[1] for line in dag.splitlines() if line.startswith("JOB ")] == ["INIT", "STEER"]
        assert "PARENT INIT CHILD STEER" in dag
        return
    sample_path = tmp_path / f"{environment['RUN_NAME']}.{CAMPAIGN_FP[:16]}.sample.condor"
    bank_path = tmp_path / f"{environment['RUN_NAME']}.{CAMPAIGN_FP[:16]}.bank.condor"
    counts = schedulers.bank_allocations(samples, jobs, cpus)
    assert sum(counts) <= max(jobs, len(samples))
    assert f"MAXJOBS PREP {jobs}" in dag
    for index, shards in enumerate(counts):
        if index < len(completed) and state != "rebuild":
            assert f"JOB SAMPLE_{index} " not in dag
            assert f"JOB BANK_{index}_" not in dag
            assert f"JOB FINALIZE_{index} " not in dag
            continue
        sample = f"SAMPLE_{index}"
        banks = [f"BANK_{index}_{shard}" for shard in range(shards)]
        assert f"JOB {sample} {sample_path}" in dag
        assert f"PARENT {sample} CHILD {' '.join(banks)}" in dag
        if shards > 1:
            finalize = f"FINALIZE_{index}"
            assert f"PARENT {' '.join(banks)} CHILD {finalize}" in dag
            assert f"PARENT {finalize} CHILD INIT" in dag
            assert f"RETRY {finalize} 2 UNLESS-EXIT 64" in dag
            assert f"ABORT-DAG-ON {finalize} 64 RETURN 1" in dag
            assert f"CATEGORY {finalize} PREP" in dag
            assert f"PARENT {' '.join(banks)} CHILD INIT" not in dag
        else:
            assert f"PARENT {banks[0]} CHILD INIT" in dag
            assert f"JOB FINALIZE_{index} " not in dag
        covered = []
        for shard, bank in enumerate(banks):
            columns = list(range(shard, components[index], shards))
            covered.extend(columns)
            assert columns
            assert f"JOB {bank} {bank_path}" in dag
            assert f"RETRY {bank} 2 UNLESS-EXIT 64" in dag
            assert f"CATEGORY {bank} PREP" in dag
            assert f'BankCpu="{min(cpus, len(columns))}"' in next(
                line for line in dag.splitlines() if line.startswith(f"VARS {bank} "))
        assert sorted(covered) == list(range(components[index]))
    assert f"JOB INIT {init_path}" in dag
    assert f"JOB STEER {steer_path}" in dag
    assert "PARENT INIT CHILD STEER" in dag
    bank = bank_path.read_text(encoding="utf-8")
    assert "arguments      = $(BankIndex) $(Shard) $(Shards)" in bank
    assert "request_cpus   = $(BankCpu)" in bank
    assert "queue 1" in bank
    sample = sample_path.read_text(encoding="utf-8")
    assert "arguments      = $(BankIndex) $(Shards)" in sample
    assert f"request_cpus   = {cpus}" in sample
    assert "queue 1" in sample


# Check a single init job keeps the plain INIT to STEER dependency
def test_lxplus_dag_single_init_job_no_bank(tmp_path):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(
        catalog, campaign_name="tune-gpom-star-cms-ampfit",
        environment_overrides={"RAY_INIT_JOBS": "1"},
    )
    environment.update(
        {
            "GRANIITTI_SHARED_OUTPUT_DIR": str(tmp_path / "shared"),
            "RAY_INIT_RUNTIME_ARCHIVE": str(tmp_path / "runtime.tar.zst"),
            "RAY_INIT_RUNTIME_SHA256": "a" * 64,
            "RAY_CAMPAIGN_FINGERPRINT": CAMPAIGN_FP,
            "RAY_LXPLUS_SCHEDD": "bigbird17.cern.ch",
            "RAY_LXPLUS_SCHEDD_ADDR": "<188.185.149.150:9618?noUDP>",
            "RAY_REPO_DIR": str(ROOT),
        }
    )
    archive = Path(environment["RAY_INIT_RUNTIME_ARCHIVE"])
    archive.write_bytes(b"runtime")

    schedulers.write_lxplus_dag(
        environment=environment,
        output_dir=tmp_path,
        runtime_archive=archive,
        runtime_sha256="a" * 64,
    )

    dag = (tmp_path / f"{environment['RUN_NAME']}.{CAMPAIGN_FP[:16]}.dag").read_text(encoding="utf-8")
    assert "JOB BANK" not in dag
    assert "JOB SAMPLE" not in dag
    assert "PARENT INIT CHILD STEER" in dag
    assert not list(tmp_path.glob("*.bank.condor"))
    assert not list(tmp_path.glob("*.sample.condor"))


# Check elastic submission omits INIT while the general switch can restore it
@pytest.mark.parametrize("campaign", ["tune-elastic-single", "tune-elastic-double", "tune-elastic-triple"])
def test_elastic_submission_skips_precomputation(tmp_path, campaign):
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name=campaign)
    assert environment["PRECOMPUTE"] == "0"
    environment.update({
        "GRANIITTI_SHARED_OUTPUT_DIR": str(tmp_path / "shared"),
        "RAY_INIT_RUNTIME_SHA256": "a" * 64,
        "RAY_CAMPAIGN_FINGERPRINT": CAMPAIGN_FP,
        "RAY_REPO_DIR": str(ROOT),
    })
    init_path, steer_path, dag_path = schedulers.write_lxplus_dag(
        environment=environment, output_dir=tmp_path,
        runtime_archive=tmp_path / "runtime.tar.zst", runtime_sha256="a" * 64,
    )
    assert init_path is None
    assert f"JOB STEER {steer_path}" in dag_path.read_text()
    assert "INIT" not in dag_path.read_text()
    assert not list(tmp_path.glob("*.init.condor"))
    assert "PRECOMPUTE=0" in steer_path.read_text()
    restored = campaign_config.resolve(
        catalog, campaign_name=campaign, environment_overrides={"PRECOMPUTE": "1"},
    )
    assert restored["PRECOMPUTE"] == "1"
    with pytest.raises(ValueError, match="requires PRECOMPUTE=1"):
        campaign_config.resolve(
            catalog, campaign_name=campaign, environment_overrides={"DATA_COVARIANCE_MODE": "full"},
        )


# Check INIT and STEER consume one identical content-addressed runtime
def test_lxplus_init_and_steer_share_runtime(tmp_path, monkeypatch):
    source = tmp_path / "source"
    runtime_input = source / "icetune"
    runtime_input.mkdir(parents=True)
    (runtime_input / "input.txt").write_text("same-runtime\n", encoding="utf-8")
    output_dir = tmp_path / "runs" / "nested" / "submit"
    monkeypatch.setattr(run_runtime, "BASE_PATHS", ("icetune",))

    tunesetup_path = tmp_path / "tunesetup.json"
    tunesetup = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source("tune-gpom-res-con"))
    save_tunesetup(tunesetup=tunesetup, path=tunesetup_path, simdriver="GRANIITTI", tune_default="TUNE0")
    archive, checksum = run_runtime.prepare_lxplus_runtime(
        source=source,
        output_dir=output_dir,
        run_name="run",
        tunesetup_path=tunesetup_path,
    )
    extracted = run_runtime.extract_lxplus_runtime(
        archive=archive,
        target=tmp_path / "extracted",
        checksum=checksum,
    )

    assert (extracted / "tmp/icetune/tunesetup.json").read_bytes() == tunesetup_path.read_bytes()
    assert checksum == sha256_file(archive)
    assert (extracted / "icetune" / "input.txt").read_text() == "same-runtime\n"


# Check the Ray steering process refuses another backend's bootstrap metadata
def test_lxplus_rejects_incompatible_init_metadata(tmp_path):
    runtime = tmp_path / "runtime"
    coord = tmp_path / "ray" / "init"
    coord.mkdir(parents=True)
    (coord / "init.json").write_text(
        json.dumps(
            {
                "backend": "other",
                "bootstrap": {"fingerprint": "physics-fingerprint"},
                "protocol": "core.tune.other.bootstrap",
                "protocol_version": 1,
                "status": "completed",
            }
        ),
        encoding="utf-8",
    )

    with pytest.raises(RuntimeError, match="valid Ray bootstrap"):
        run_runtime.stage_ray_bootstrap(
            runtime=runtime,
            coord_dir=coord,
            run_name="run",
            simdriver="GRANIITTI",
        )


# Check the common INIT staging accepts PANDORA's no-output bootstrap
def test_lxplus_stages_pandora_init_metadata(tmp_path):
    runtime = tmp_path / "runtime"
    coord = tmp_path / "ray" / "init"
    coord.mkdir(parents=True)
    (coord / "init.json").write_text(
        json.dumps(
            {
                "backend": "ray",
                "bootstrap": {
                    "archive_sha256": None,
                    "archive_size": 0,
                    "archive_url": None,
                    "files": [],
                    "fingerprint": "physics-fingerprint",
                    "kind": "pandora",
                },
                "protocol": "core.tune.ray.bootstrap",
                "protocol_version": 1,
                "status": "completed",
            }
        ),
        encoding="utf-8",
    )

    fingerprint = run_runtime.stage_ray_bootstrap(
        runtime=runtime,
        coord_dir=coord,
        run_name="run",
        simdriver="PANDORA",
    )

    assert fingerprint == "physics-fingerprint"
    assert (runtime / "runs" / "icetune" / "run" / "ray_init.json").is_file()


# Check the GRANIITTI head reads HEPData and the saved initial point without generating trials
def test_ray_lxplus_head_loads_init_point(tmp_path, data_card):
    shutil.copytree(ROOT / "modeldata/TUNE0", tmp_path / "modeldata/TUNE0", dirs_exist_ok=True)
    steering = {"tune_default": "TUNE0", "data_covariance_mode": "diagonal", "ampfit_reuse_root": "/eos/shared/banks"}
    tunesetup = SimpleNamespace(datacards=[{"datacard": str(data_card)}])
    args = SimpleNamespace(cdir=str(tmp_path), obs_module="default", pickle_dump=False,
                           ray_init_state=str(tmp_path / "ray_init.json"), run_name="head", runtime_sha256="a" * 64)
    prepared = GraniittiDriver()
    prepared.init_data(run_name=args.run_name, datacards=tunesetup.datacards, obs_module=args.obs_module,
                       cdir=args.cdir, pickle_dump=False)
    fingerprint = prepared.bootstrap_fingerprint(cdir=args.cdir, datacards=tunesetup.datacards,
                                                mc_steer=steering, runtime_sha256=args.runtime_sha256, tunesetup_path=None)
    initial = {"parameter": 0.42}
    Path(args.ray_init_state).write_text(json.dumps({
        "backend": "ray", "bootstrap": {"fingerprint": fingerprint}, "initial_points": initial,
        "protocol": "core.tune.ray.bootstrap", "protocol_version": 1, "status": "completed"}))
    driver = GraniittiDriver()
    steering["ampfit_reuse_root"] = None
    result = driver.initialize_ray_bootstrap(args=args, tunesetup=tunesetup, mc_steer=steering)
    assert result == initial and result is not initial
    assert driver.initialized and driver.data
    assert driver.dataset_paths == prepared.dataset_paths
    assert not list(tmp_path.rglob("*.vgrid"))
    assert not list(tmp_path.rglob("*.hepmc3"))


# Check the Pandora head restores its actual initialization record without reconstruction
def test_lxplus_pandora_head_loads_init_point(tmp_path):
    tunesetup = load_tunesetup(cdir=ROOT, simdriver="PANDORA", name=campaign_source("tune-pandora-v0"))
    driver = PandoraDriver()
    args = SimpleNamespace(cdir=str(tmp_path), obs_module="default", pickle_dump=False,
                           ray_init_state=str(tmp_path / "ray_init.json"), run_name="head", runtime_sha256="a" * 64)
    bootstrap = driver.bootstrap_payload(runtime_sha256=args.runtime_sha256, datacards=tunesetup.datacards)
    initial = {"parameter": 0.42}
    Path(args.ray_init_state).write_text(json.dumps({
        "backend": "ray", "bootstrap": bootstrap, "initial_points": initial,
        "protocol": "core.tune.ray.bootstrap", "protocol_version": 1, "status": "completed"}))
    result = driver.initialize_ray_bootstrap(args=args, tunesetup=tunesetup, mc_steer={})
    assert result == initial and result is not initial
    assert driver.initialized and driver.datacards == tunesetup.datacards
    assert not list(tmp_path.rglob("*.root"))


# Check Ray slots may share workers but never use the steering machine
def test_lxplus_accepts_packed_machine_placement():
    class FakeCluster:
        # Build one mutable simulated Dask worker specification
        def __init__(self):
            self.new_spec = {"options": {"job_extra_directives": {}}}
            self.jobs = 0

        # Record one concurrent scale request
        def scale(self, *, jobs):
            self.jobs = jobs

        # Complete one simulated worker submission update
        async def _correct_state(self):
            return None

        # Run one simulated synchronous cluster operation
        def sync(self, function):
            return asyncio.run(function())

    class FakeClient:
        # Store one simulated CERN cluster
        def __init__(self, cluster, machines):
            self.cluster = cluster
            self.machines = machines

        # Check only the first worker request is awaited
        def wait_for_workers(self, workers, *, timeout):
            assert workers == 1
            assert self.cluster.jobs == 1
            assert timeout == 30.0

        # Compute the complete simulated Dask worker set
        def scheduler_info(self):
            return {"workers": {"worker-0": {}, "worker-1": {}}}

        # Compute the selected physical HTCondor machines
        def run(self, function, *, workers):
            if function is ray_nodes.condor_machine:
                return dict(zip(workers, self.machines, strict=True))
            assert function is ray_nodes.condor_job_identity
            return {
                worker: {"cluster_id": 1, "job_id": f"1.{index}", "proc_id": index}
                for index, worker in enumerate(workers)
            }

    cluster = FakeCluster()
    ray_cluster.submit_ray_slots(cluster=cluster, workers=2, cpus=8, steering_machine="steer.cern.ch")
    workers = ray_cluster.wait_ray_slots(client=FakeClient(cluster, ("worker-a.cern.ch", "worker-a.cern.ch")),
                                        timeout=30.0, steering_machine="steer.cern.ch")
    assert workers == ["worker-0", "worker-1"]
    directives = cluster.new_spec["options"]["job_extra_directives"]
    assert directives["requirements"] == '(TARGET.Machine =!= "steer.cern.ch")'
    assert "rank" not in directives
    assert cluster.new_spec["group"] == ["-0", "-1", "-2"]

    rejected = FakeCluster()
    with pytest.raises(RuntimeError, match="steering machine"):
        ray_cluster.submit_ray_slots(cluster=rejected, workers=2, cpus=8, steering_machine="steer.cern.ch")
        ray_cluster.wait_ray_slots(client=FakeClient(rejected, ("steer.cern.ch", "worker-b.cern.ch")),
                                   timeout=30.0, steering_machine="steer.cern.ch")


# Check nested condor_submit failures leave the Dask reconciliation call
def test_lxplus_propagates_worker_submission_failure():
    events = []

    class FakeCluster:
        # Build one mutable simulated Dask worker specification
        def __init__(self):
            self.new_spec = {"options": {"job_extra_directives": {}}}

        # Record that scale runs inside the synchronous Dask loop call
        def scale(self, *, jobs):
            assert jobs == 1
            asyncio.get_running_loop()
            events.append("scale")

        # Fail one simulated nested HTCondor submission
        async def _correct_state(self):
            events.append("submit")
            raise RuntimeError("condor_submit failed")

        # Run one simulated synchronous cluster operation
        def sync(self, function):
            events.append("sync")
            return asyncio.run(function())

    class FakeClient:
        # Reject any wait after a failed worker submission
        def wait_for_workers(self, workers, *, timeout):
            pytest.fail("failed worker submissions must not enter the startup wait")

    cluster = FakeCluster()
    with pytest.raises(RuntimeError, match="condor_submit failed"):
        ray_cluster.submit_ray_slots(cluster=cluster, workers=2, cpus=8, steering_machine="steer.cern.ch")
    assert events == ["sync", "scale", "submit"]
    assert cluster.new_spec["group"] == ["-0", "-1", "-2"]


# Check the steering process tests direct schedd write access before staging
def test_ray_lxplus_checks_submit_access(monkeypatch):
    address = "<188.185.149.150:9618?noUDP>"
    commands = []

    # Record one successful HTCondor authorization check
    def fake_run(command, **kwargs):
        commands.append(command)
        return SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)
    ray_nodes.validate_condor_write(address)

    assert commands == [["condor_ping", "-address", address, "-quiet", "WRITE"]]


# Check lxplus staging copies only the selected runtime paths
def test_ray_lxplus_stages_compact_runtime(tmp_path, monkeypatch):
    source = tmp_path / "source"
    target = tmp_path / "target"
    (source / "bin").mkdir(parents=True)
    (source / "bin" / "gr").write_text("binary", encoding="utf-8")
    (source / "python").mkdir()
    (source / "python" / "module.py").write_text("VALUE = 1\n", encoding="utf-8")
    (source / "output").mkdir()
    (source / "output" / "large.raw").write_text("excluded", encoding="utf-8")
    monkeypatch.setattr(run_runtime, "BASE_PATHS", ("bin/gr", "python"))

    runtime = run_runtime.stage_runtime(source=source, target=target)

    assert (runtime / "bin" / "gr").read_text(encoding="utf-8") == "binary"
    assert (runtime / "python" / "module.py").is_file()
    assert not (runtime / "output").exists()

    shared = tmp_path / "graniitti"
    prior_state = shared / "runs" / "icetune" / "run" / "tuner.pkl"
    prior_state.parent.mkdir(parents=True)
    prior_state.write_text("state", encoding="utf-8")
    checkpoint_payload = run_outputs.pack_output_tree(prior_state.parent, Path("run"))
    run_outputs.publish_checkpoint(
        payload=checkpoint_payload,
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    staged = run_outputs.stage_previous_outputs(
        shared_root=shared,
        runtime=runtime,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    assert staged == [run_outputs.checkpoint_path(shared, "run", CAMPAIGN_FP)]
    assert (runtime / "runs" / "icetune" / "run" / "tuner.pkl").is_file()
    assert not (runtime / "figs").exists()

    packed = run_runtime.pack_runtime(runtime)
    with tarfile.open(packed, mode="r:gz") as archive:
        assert "runtime/bin/gr" in archive.getnames()


# Check an incompatible campaign cannot observe another campaign checkpoint
def test_lxplus_campaign_state_content_addressed(tmp_path):
    shared = tmp_path / "shared"
    source = tmp_path / "source" / "run"
    source.mkdir(parents=True)
    (source / "tuner.pkl").write_text("campaign-a", encoding="utf-8")
    payload = run_outputs.pack_output_tree(source, Path("run"))
    campaign_a = "a" * 64
    campaign_b = "b" * 64
    run_outputs.publish_checkpoint(
        payload=payload,
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=campaign_a,
    )

    other_runtime = tmp_path / "other"
    assert (
        run_outputs.stage_previous_outputs(
            shared_root=shared,
            runtime=other_runtime,
            run_name="run",
            campaign_fingerprint=campaign_b,
        )
        == []
    )
    assert not (other_runtime / "runs" / "icetune" / "run" / "tuner.pkl").exists()

    matching_runtime = tmp_path / "matching"
    assert run_outputs.stage_previous_outputs(
        shared_root=shared,
        runtime=matching_runtime,
        run_name="run",
        campaign_fingerprint=campaign_a,
    ) == [run_outputs.checkpoint_path(shared, "run", campaign_a)]
    assert (matching_runtime / "runs" / "icetune" / "run" / "tuner.pkl").read_text(
        encoding="utf-8"
    ) == "campaign-a"


# Check Ray binds to the address already reachable by the Dask scheduler
def test_ray_lxplus_uses_dask_worker_address(monkeypatch):
    worker = SimpleNamespace(address="tcp://worker-control.cern.ch:10042")
    monkeypatch.setitem(
        sys.modules,
        "distributed",
        SimpleNamespace(get_worker=lambda: worker),
    )
    monkeypatch.setattr(socket, "gethostbyname", lambda host: "188.184.1.2")
    assert ray_nodes.dask_node_ip() == "188.184.1.2"

    monkeypatch.setattr(socket, "gethostbyname", lambda host: "169.254.3.1")
    with pytest.raises(RuntimeError, match="not routable"):
        ray_nodes.dask_node_ip()


# Check Ray head services use distinct explicit ports from the CERN range
@pytest.mark.parametrize("failure", [None, "exit", "timeout"])
def test_ray_lxplus_node_ports_and_gpus(tmp_path, monkeypatch, failure):
    monkeypatch.setenv("RAY_HEAD_REQUEST_MEMORY", "8GB")
    monkeypatch.setenv("RAY_REQUEST_MEMORY", "4GB")
    monkeypatch.setattr(ray_nodes, "ray_memory", lambda requested: {"memory": 1024**3, "object_store_memory": 128 * 1024**2})
    class FakeSocket:
        # Enter one simulated free port probe
        def bind(self, address):
            self.address = address

        # Enter one simulated socket lifetime
        def __enter__(self):
            return self

        # Close one simulated socket lifetime
        def __exit__(self, *args):
            return None

        # Close one simulated free port probe
        def close(self):
            return None

    commands = []
    environments = []

    # Record one simulated Ray startup command
    def fake_run(command, **kwargs):
        commands.append(command)
        environments.append(kwargs["env"])
        if failure:
            root = Path(next(value.split("=", 1)[1] for value in command if value.startswith("--temp-dir=")))
            logs = root / "session_2026_09_08" / "logs"
            logs.mkdir(parents=True)
            (logs / "gcs_server.out").write_text("GCS diagnostic")
            (logs / "raylet.err").write_bytes(b"x" * 40000 + b"Raylet diagnostic\xff")
            if failure == "timeout":
                raise subprocess.TimeoutExpired(command, 1200, output=b"launcher stdout", stderr=b"launcher stderr")
            return SimpleNamespace(returncode=1, stdout="launcher stdout", stderr="launcher stderr")
        return SimpleNamespace(returncode=0, stderr="", stdout="")

    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path))
    monkeypatch.delenv("ICETUNE_CONDA_PREFIX", raising=False)
    monkeypatch.delenv("ICETUNE_RAY_BIN", raising=False)
    monkeypatch.setenv("RAY_LXPLUS_PORT_DIR", str(tmp_path))
    monkeypatch.setenv("PYTHONHOME", "/cvmfs/lcg")
    monkeypatch.setenv("PYTHONPATH", "/cvmfs/lcg/site-packages")
    monkeypatch.setattr(run_outputs.shutil, "which", lambda executable, path=None: executable)
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(ray_nodes, "dask_node_ip", lambda: "192.0.2.1")
    socket_type = socket.socket
    monkeypatch.setattr(socket, "socket", lambda *args: socket_type(*args) if args[0] == socket.AF_UNIX else FakeSocket())
    monkeypatch.setattr(socket, "getfqdn", lambda: "node.cern.ch")
    monkeypatch.setattr(socket, "gethostbyname", lambda host: "192.0.2.1")

    requested_tmp = tmp_path / "irt-0123456789ab"
    ray_tmp = tmp_path / requested_tmp.name
    if failure:
        with pytest.raises(RuntimeError) as caught:
            ray_nodes.start_ray_node(
                role="head", cpus=8, temp_dir=str(ray_tmp), gpus=1,
                conda_prefix="/opt/conda/envs/graniitti", ray_bin="/opt/conda/envs/graniitti/bin/ray",
            )
        detail = str(caught.value)
        assert "GCS diagnostic" in detail
        assert "Raylet diagnostic" in detail
        assert "launcher stdout" in detail and "launcher stderr" in detail
        logs = ray_nodes.ray_log_tails(ray_tmp / "icetune_head/ray/session_2026_09_08/logs")
        assert len(logs["raylet.err"]) == 32768
        assert logs["raylet.err"].endswith("Raylet diagnostic\ufffd")
        assert "timed out" in detail if failure == "timeout" else "startup failed" in detail
        return
    result = ray_nodes.start_ray_node(
        role="head",
        cpus=8,
        temp_dir=str(ray_tmp),
        gpus=1,
        conda_prefix="/opt/conda/envs/graniitti",
        ray_bin="/opt/conda/envs/graniitti/bin/ray",
    )
    command = commands[0]
    service_options = (
        "--node-manager-port=",
        "--object-manager-port=",
        "--runtime-env-agent-port=",
        "--dashboard-agent-listen-port=",
        "--dashboard-agent-grpc-port=",
        "--metrics-export-port=",
        "--port=",
        "--ray-client-server-port=",
        "--dashboard-port=",
    )
    ports = [
        int(argument.split("=", 1)[1])
        for argument in command
        if argument.startswith(service_options)
    ]
    worker_argument = next(
        argument for argument in command if argument.startswith("--worker-port-list=")
    )
    ports.extend(int(port) for port in worker_argument.split("=", 1)[1].split(","))

    assert result["head_address"] == "192.0.2.1:10016"
    node_tmp = ray_tmp / "icetune_head"
    assert result["ray_tmp"] == str(node_tmp / "ray")
    assert "--ray-client-server-port=10017" in command
    assert "--dashboard-port=10018" in command
    assert f"--temp-dir={node_tmp / 'ray'}" in command
    assert "--object-manager-port=10011" in command
    assert "--num-gpus=0" in command
    assert '--resources={"icetune_head":1}' in command
    assert len(ports) == len(set(ports))
    assert all(10010 <= port <= 10100 for port in ports)
    assert environments[0]["PATH"].startswith("/opt/conda/envs/graniitti/bin:")
    assert environments[0]["TMPDIR"] == str(tmp_path)
    assert environments[0]["XDG_CACHE_HOME"] == str(tmp_path / "cache")
    assert environments[0]["RAY_TMPDIR"] == str(node_tmp)
    assert environments[0]["RAY_raylet_start_wait_time_s"] == "600"
    assert "PYTHONHOME" not in environments[0]
    assert "PYTHONPATH" not in environments[0]


# Check packed Ray nodes retain disjoint host ports and release failed leases
def test_ray_lxplus_packed_node_port_leases(tmp_path, monkeypatch):
    monkeypatch.setenv("RAY_HEAD_REQUEST_MEMORY", "8GB")
    monkeypatch.setenv("RAY_REQUEST_MEMORY", "4GB")
    monkeypatch.setattr(ray_nodes, "ray_memory", lambda requested: {"memory": 1024**3, "object_store_memory": 128 * 1024**2})
    class FakeSocket:
        # Accept one simulated free port probe
        def bind(self, address):
            self.address = address

        # Enter one simulated socket lifetime
        def __enter__(self):
            return self

        # Close one simulated socket lifetime
        def __exit__(self, *args):
            return None

        # Close one simulated free port probe
        def close(self):
            return None

    commands = []
    environments = []
    current_pid = [110001]
    fail = [False]

    # Record one simulated Ray startup
    def fake_run(command, **kwargs):
        commands.append(command)
        environments.append(kwargs["env"])
        return SimpleNamespace(
            returncode=int(fail[0]),
            stderr="simulated startup failure" if fail[0] else "",
            stdout="",
        )

    # Compute every explicit port from one Ray command
    def command_ports(command):
        prefixes = (
            "--node-manager-port=",
            "--object-manager-port=",
            "--runtime-env-agent-port=",
            "--dashboard-agent-listen-port=",
            "--dashboard-agent-grpc-port=",
            "--metrics-export-port=",
            "--port=",
            "--ray-client-server-port=",
            "--dashboard-port=",
        )
        ports = {
            int(argument.split("=", 1)[1]) for argument in command if argument.startswith(prefixes)
        }
        worker_argument = next(
            argument for argument in command if argument.startswith("--worker-port-list=")
        )
        ports.update(int(port) for port in worker_argument.split("=", 1)[1].split(","))
        return ports

    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "scratch"))
    monkeypatch.setenv("ICETUNE_CONDA_PREFIX", "/opt/conda/envs/graniitti")
    monkeypatch.setenv("ICETUNE_RAY_BIN", "/opt/conda/envs/graniitti/bin/ray")
    monkeypatch.setenv("RAY_LXPLUS_PORT_DIR", str(tmp_path / "ports"))
    monkeypatch.setattr(run_outputs.os, "getpid", lambda: current_pid[0])
    monkeypatch.setattr(run_common, "process_birth", lambda pid: str(pid * 10))
    monkeypatch.setattr(run_outputs.shutil, "which", lambda executable, path=None: executable)
    monkeypatch.setattr(subprocess, "run", fake_run)
    monkeypatch.setattr(ray_nodes, "dask_node_ip", lambda: "192.0.2.2")
    socket_type = socket.socket
    monkeypatch.setattr(socket, "socket", lambda *args: socket_type(*args) if args[0] == socket.AF_UNIX else FakeSocket())
    monkeypatch.setattr(socket, "getfqdn", lambda: "node.cern.ch")
    monkeypatch.setattr(socket, "gethostbyname", lambda host: "192.0.2.1")

    requested_tmp = tmp_path / "irt-0123456789ab"
    ray_tmp = requested_tmp
    port_root = tmp_path / "ports"
    port_root.mkdir()
    (port_root / "icetune-ray-ports.txt").write_text(
        "110001 999 0123456789ab 10000\n", encoding="utf-8"
    )
    head = ray_nodes.start_ray_node(role="head", cpus=1, temp_dir=str(requested_tmp))
    current_pid[0] = 110002
    ray_nodes.start_ray_node(
        role="worker",
        cpus=1,
        temp_dir=str(requested_tmp),
        head_address=head["head_address"],
        resource_marker="icetune_worker_244912_1",
    )
    assert command_ports(commands[0]).isdisjoint(command_ports(commands[1]))
    assert any(argument.startswith("--ray-client-server-port=") for argument in commands[1])
    assert f"--temp-dir={ray_tmp / 'icetune_head' / 'ray'}" in commands[0]
    assert not any(argument.startswith("--temp-dir=") for argument in commands[1])
    assert environments[0]["RAY_TMPDIR"] == str(ray_tmp / "icetune_head")
    assert environments[1]["RAY_TMPDIR"] == str(ray_tmp / "icetune_worker_244912_1")
    spill_dir = tmp_path / "scratch" / "ray-spill" / "0123456789ab"
    assert f"--object-spilling-directory={spill_dir}" in commands[0]
    assert f"--object-spilling-directory={spill_dir}" in commands[1]
    assert "--include-log-monitor=False" in commands[1]
    assert '--resources={"icetune_worker_244912_1":1}' in commands[1]

    current_pid[0] = 110003
    monkeypatch.setenv("RAY_LXPLUS_PORT_DIR", str(tmp_path / "another-host"))
    monkeypatch.setattr(ray_nodes, "dask_node_ip", lambda: "192.0.2.3")
    fail[0] = True
    with pytest.raises(RuntimeError, match="simulated startup failure"):
        ray_nodes.start_ray_node(
            role="worker",
            cpus=1,
            temp_dir=str(requested_tmp),
            head_address=head["head_address"],
            resource_marker="icetune_worker_244912_2",
        )
    fail[0] = False
    ray_nodes.start_ray_node(role="worker", cpus=1, temp_dir=str(requested_tmp),
                         head_address=head["head_address"], resource_marker="icetune_worker_244912_2")
    assert command_ports(commands[-1]) == command_ports(commands[-2])
    with (
        pytest.raises(RuntimeError, match="but only"),
        ray_nodes.reserve_ray_ports("192.0.2.1", 102, "0123456789ab"),
    ):
        pytest.fail("an impossible CERN port request must not be entered")


# Check host reservations survive separate Condor temporary directories
def test_lxplus_cross_namespace_port_lease(tmp_path, monkeypatch):
    monkeypatch.setenv("RAY_LXPLUS_PORT_DIR", str(tmp_path))
    try:
        with socket.socket() as probe:
            probe.bind(("127.0.0.1", 0))
        with ray_nodes.reserve_ray_ports("127.0.0.1", 2, "0123456789ab") as (first, monitor):
            assert monitor
        monkeypatch.setenv("RAY_LXPLUS_PORT_DIR", str(tmp_path / "other-job"))
        with ray_nodes.reserve_ray_ports("127.0.0.1", 2, "0123456789ab") as (second, monitor):
            assert monitor
            assert set(first).isdisjoint(second)
    except PermissionError:
        pytest.skip("Port allocation requires local socket support")


# Check a figure copy race cannot suppress changed status and history files
def test_lxplus_live_figure_race_json(tmp_path, monkeypatch):
    environment = {"RUN_NAME": "run"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    runtime = run_common.runtime_work_root(environment) / "runtime"
    figures = runtime / "figs" / "icetune" / "run"
    state = runtime / "runs" / "icetune" / "run"
    figures.mkdir(parents=True)
    state.mkdir(parents=True)
    (figures / "summary.json").write_text('{"trial_id":"best"}\n', encoding="utf-8")
    (state / "status.json").write_text('{"completed":7}\n', encoding="utf-8")
    (state / "history.json").write_text('{"history_schema_version":1}\n', encoding="utf-8")
    pack_outputs = run_outputs.pack_runtime_outputs

    # Fail only the independently packaged figure tree
    def fail_figures(**kwargs):
        if kwargs["figure_dir"] is not None:
            raise RuntimeError("figure changed")
        return pack_outputs(**kwargs)

    monkeypatch.setattr(run_outputs, "pack_runtime_outputs", fail_figures)
    snapshot = run_outputs.snapshot_live_outputs(
        environment=environment,
        figure_revision=None,
    )

    assert snapshot["figure_outputs"] is None
    assert "figures" in snapshot["deferred"]
    with tarfile.open(fileobj=io.BytesIO(snapshot["outputs"]), mode="r:gz") as archive:
        names = set(archive.getnames())
    assert "runs/icetune/run/status.json" in names
    assert "runs/icetune/run/history.json" in names


# Check final checkpoint packaging refuses a live icetune driver record
def test_lxplus_checkpoint_stopped_driver(tmp_path, monkeypatch):
    environment = {"RUN_NAME": "run"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    pid_path = run_common.runtime_pid_path(environment)
    pid_path.parent.mkdir(parents=True)
    pid_path.write_text('{"birth":"1234","pid":4321}\n', encoding="utf-8")

    with pytest.raises(RuntimeError, match="Refusing to checkpoint"):
        run_outputs.snapshot_final_outputs(environment=environment)


# Check every history plot transfers independently of the unchanged best trial
@pytest.mark.parametrize("name", ICETUNE_HISTORY_FIGURES)
def test_lxplus_history_transfer(tmp_path, name):
    """Transfer changed history plots while preserving existing best trial files"""
    figures = tmp_path / "figures"
    figures.mkdir()
    (figures / "summary.json").write_text('{"trial_id":"best-1"}')
    (figures / "best.pdf").write_bytes(b"%PDF-best")
    for filename in ICETUNE_HISTORY_FIGURES:
        signature = b"%PDF-" if filename.endswith(".pdf") else b"\x89PNG\r\n\x1a\n"
        (figures / filename).write_bytes(signature + b"history-1")
    initial, _, revision, deferred = run_outputs.snapshot_live_figures(
        figure_dir=figures, run_name="run", figure_revision=None)
    assert not deferred
    publication = dict(shared_root=tmp_path / "shared", run_name="run", campaign_fingerprint=CAMPAIGN_FP)
    run_outputs.restore_outputs(payload=initial, **publication)
    assert not (publication["shared_root"] / "runs/icetune/run/ray/submit").exists()
    target = run_outputs.campaign_figure_path(publication["shared_root"], "run", CAMPAIGN_FP)
    best_stat = (target / "best.pdf").stat()
    updated = (figures / name).read_bytes() + b"updated"
    (figures / name).write_bytes(updated)
    best, history, revision, deferred = run_outputs.snapshot_live_figures(
        figure_dir=figures, run_name="run", figure_revision=revision)
    assert best is None
    assert history == {name: updated}
    assert not deferred
    run_outputs.publish_history_figures(payload=history, **publication)
    assert (target / name).read_bytes() == updated
    assert (target / "best.pdf").stat() == best_stat
    best, history, _, deferred = run_outputs.snapshot_live_figures(
        figure_dir=figures, run_name="run", figure_revision=revision)
    assert best is None and history is None and not deferred
    (figures / "summary.json").write_text('{"trial_id":"best-2"}')
    best, history, _, deferred = run_outputs.snapshot_live_figures(
        figure_dir=figures, run_name="run", figure_revision=revision)
    assert best is not None and history is None and not deferred


# Check live returns stay incremental and final output contains restorable Tune state
def test_lxplus_live_figures_final_checkpoint(tmp_path, monkeypatch):
    environment = {"COST": "loss", "PLOT": "1", "RUN_NAME": "run"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    runtime = run_common.runtime_work_root(environment) / "runtime"
    figures = runtime / "figs" / "icetune" / "run"
    state = runtime / "runs" / "icetune" / "run"
    figures.mkdir(parents=True)
    state.mkdir(parents=True)
    (figures / "summary.json").write_text('{"trial_id":"best-1"}\n', encoding="utf-8")
    (figures / "cost_evolution.png").write_bytes(b"\x89PNG\r\n\x1a\nhistory-1")
    (state / "tuner.pkl").write_bytes(b"tuner-state")
    (state / "experiment_state.json").write_text("{}\n", encoding="utf-8")
    (state / "status.json").write_text(
        json.dumps(
            {
                "phase": "running",
                "ray_resources": {
                    "available": {"CPU": 8.0},
                    "cluster": {
                        "CPU": 8.0,
                        "icetune_head": 1.0,
                        "icetune_worker_252171_1": 1.0,
                    },
                    "nodes": {},
                },
                "ray_resources_observed_at_unix": 122.0,
                "updated_at_unix": 123.0,
            }
        ),
        encoding="utf-8",
    )
    (state / "history.json").write_text('{\n  "history_schema_version": 1\n}\n', encoding="utf-8")
    failures = state / "failures"
    failures.mkdir()
    (failures / "trial-000007.json").write_text(
        '{\n  "error": "complete simulator traceback"\n}\n', encoding="utf-8"
    )
    results = state / "results"
    results.mkdir()
    result_name = "TUNE_icetune_trial-000001.pkl"
    (results / result_name).write_bytes(b"replica-1")
    nested = state / "CFunc_trial" / "results"
    nested.mkdir(parents=True)
    (nested / result_name).write_bytes(b"replica-1")
    nested_figures = state / "CFunc_trial" / "figs" / "icetune" / "run"
    nested_figures.mkdir(parents=True)
    (nested_figures / "summary.json").write_text("{}\n", encoding="utf-8")

    first = run_outputs.snapshot_live_outputs(
        environment=environment,
        figure_revision=None,
    )
    assert first["outputs"] is not None
    assert first["figure_outputs"] is not None
    with tarfile.open(fileobj=io.BytesIO(first["outputs"]), mode="r:gz") as archive:
        names = set(archive.getnames())
    assert "runs/icetune/run/history.json" in names
    assert "runs/icetune/run/status.json" in names
    assert "runs/icetune/run/failures/trial-000007.json" in names
    assert f"runs/icetune/run/results/{result_name}" not in names
    with tarfile.open(fileobj=io.BytesIO(first["figure_outputs"]), mode="r:gz") as archive:
        assert "figs/icetune/run/summary.json" in archive.getnames()
    assert first["ray_resources"]["cluster"]["icetune_head"] == 1.0
    assert first["ray_resources_observed_at_unix"] == 122.0
    shared = tmp_path / "shared"
    run_outputs.restore_outputs(
        payload=first["outputs"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    run_outputs.restore_outputs(
        payload=first["figure_outputs"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    run_output = run_outputs.campaign_run_path(shared, "run", CAMPAIGN_FP)
    figure_output = run_outputs.campaign_figure_path(shared, "run", CAMPAIGN_FP)
    assert not (run_output / "results" / result_name).exists()
    assert json.loads((run_output / "status.json").read_text())["phase"] == "running"
    assert (run_output / "failures" / "trial-000007.json").is_file()
    run_outputs.restore_outputs(
        payload=first["outputs"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )

    unchanged = run_outputs.snapshot_live_outputs(
        environment=environment,
        figure_revision=first["figure_revision"],
        live_revision=first["live_revision"],
    )
    assert unchanged["outputs"] is None
    assert unchanged["figure_outputs"] is None
    assert unchanged["history_figures"] is None

    (figures / "cost_evolution.png").write_bytes(b"\x89PNG\r\n\x1a\nhistory-2")
    history_only = run_outputs.snapshot_live_outputs(
        environment=environment,
        figure_revision=first["figure_revision"],
        live_revision=first["live_revision"],
    )
    assert history_only["outputs"] is None
    assert history_only["figure_outputs"] is None
    assert history_only["history_figures"] == {"cost_evolution.png": b"\x89PNG\r\n\x1a\nhistory-2"}
    run_outputs.publish_history_figures(
        payload=history_only["history_figures"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    assert (figure_output / "summary.json").is_file()
    assert (figure_output / "cost_evolution.png").read_bytes() == (b"\x89PNG\r\n\x1a\nhistory-2")

    (state / "status.json").write_text('{\n  "phase": "done"\n}\n', encoding="utf-8")
    status_only = run_outputs.snapshot_live_outputs(
        environment=environment,
        figure_revision=history_only["figure_revision"],
        live_revision=first["live_revision"],
    )
    with tarfile.open(fileobj=io.BytesIO(status_only["outputs"]), mode="r:gz") as archive:
        names = set(archive.getnames())
    assert "runs/icetune/run/status.json" in names
    assert "runs/icetune/run/history.json" not in names

    (figures / "summary.json").write_text(
        '{"metrics":{"loss":1.0},"trial_id":"best-22"}\n',
        encoding="utf-8",
    )
    result_descriptor = {
        "filename": result_name,
        "sha256": sha256_file(results / result_name),
        "size": (results / result_name).stat().st_size,
    }
    (state / "experiment_state-000.json").write_text(
        json.dumps(
            {
                "trial_data": [
                    [
                        json.dumps({"status": "TERMINATED", "trial_id": "tied-earlier"}),
                        json.dumps(
                            {
                                "last_result": {
                                    "is_initial": False,
                                    "loss": 1.0,
                                }
                            }
                        ),
                    ],
                    [
                        json.dumps({"status": "TERMINATED", "trial_id": "best-22"}),
                        json.dumps(
                            {
                                "last_result": {
                                    "is_initial": False,
                                    "loss": 1.0,
                                    "trial_pickle": result_descriptor,
                                }
                            }
                        ),
                    ],
                ]
            }
        ),
        encoding="utf-8",
    )
    checkpoint = run_outputs.snapshot_final_outputs(environment=environment)
    assert checkpoint["recovery"] is not None
    with tarfile.open(fileobj=io.BytesIO(checkpoint["outputs"]), mode="r:gz") as archive:
        names = set(archive.getnames())
    assert "figs/icetune/run/cost_evolution.png" in names
    assert f"runs/icetune/run/results/{result_name}" in names
    assert checkpoint["result_names"] == [result_name]
    with tarfile.open(fileobj=io.BytesIO(checkpoint["checkpoint"]), mode="r:gz") as archive:
        names = set(archive.getnames())
    assert "run/tuner.pkl" in names
    assert not any("/results/TUNE_icetune_" in name for name in names)
    assert not any("/figs/icetune/" in name for name in names)


# Check trial pickles leave head scratch only with a matching terminated generation
def test_lxplus_defers_uncommitted_trial_pickle(tmp_path, monkeypatch):
    environment = {"COST": "loss", "PLOT": "0", "RUN_NAME": "run"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    runtime = run_common.runtime_work_root(environment) / "runtime"
    state = runtime / "runs" / "icetune" / "run"
    results = state / "results"
    results.mkdir(parents=True)
    (state / "tuner.pkl").write_bytes(b"tuner-state")
    (state / "status.json").write_text('{"phase":"running"}\n', encoding="utf-8")
    result = results / "TUNE_icetune_trial-000001.pkl"
    result.write_bytes(b"uncommitted-at-crash")
    descriptor = {
        "filename": result.name,
        "sha256": sha256_file(result),
        "size": result.stat().st_size,
    }
    experiment = state / "experiment_state-000.json"
    experiment.write_text(
        json.dumps(
            {
                "trial_data": [
                    [
                        json.dumps({"status": "RUNNING", "trial_id": "trial-000001"}),
                        json.dumps(
                            {
                                "last_result": {
                                    "loss": 1.0,
                                    "trial_pickle": descriptor,
                                }
                            }
                        ),
                    ]
                ]
            }
        ),
        encoding="utf-8",
    )

    crashed = run_outputs.snapshot_final_outputs(environment=environment)
    assert crashed["recovery"] is not None
    assert crashed["result_names"] == []
    with tarfile.open(fileobj=io.BytesIO(crashed["outputs"]), mode="r:gz") as archive:
        assert not any("TUNE_icetune_" in name for name in archive.getnames())

    shared = tmp_path / "shared"
    run_outputs.restore_outputs(
        payload=crashed["outputs"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    run_outputs.publish_recovery_generation(
        payload=crashed["recovery"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    restored = tmp_path / "restored"
    run_outputs.stage_previous_outputs(
        shared_root=shared,
        runtime=restored,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    assert not (restored / "runs" / "icetune" / "run" / "results" / result.name).exists()
    assert not (
        run_outputs.campaign_run_path(shared, "run", CAMPAIGN_FP) / "results" / result.name
    ).exists()

    result.write_bytes(b"retried-and-terminated")
    descriptor = {
        "filename": result.name,
        "sha256": sha256_file(result),
        "size": result.stat().st_size,
    }
    experiment.write_text(
        json.dumps(
            {
                "trial_data": [
                    [
                        json.dumps({"status": "TERMINATED", "trial_id": "trial-000001"}),
                        json.dumps(
                            {
                                "last_result": {
                                    "loss": 1.0,
                                    "trial_pickle": descriptor,
                                }
                            }
                        ),
                    ]
                ]
            }
        ),
        encoding="utf-8",
    )
    committed = run_outputs.snapshot_final_outputs(environment=environment)
    assert committed["recovery"] is not None
    assert committed["result_names"] == [result.name]
    with tarfile.open(fileobj=io.BytesIO(committed["outputs"]), mode="r:gz") as archive:
        assert f"runs/icetune/run/results/{result.name}" in archive.getnames()
    run_outputs.restore_outputs(
        payload=committed["outputs"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    run_outputs.publish_recovery_generation(
        payload=committed["recovery"],
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )
    assert (
        run_outputs.campaign_run_path(shared, "run", CAMPAIGN_FP) / "results" / result.name
    ).read_bytes() == b"retried-and-terminated"


# Check restart selects one complete state and figure generation after partial returns
def test_lxplus_recovery_generation_atomic(tmp_path):
    shared = tmp_path / "shared"
    state = tmp_path / "head" / "run"
    state.mkdir(parents=True)
    (state / "tuner.pkl").write_text("matching-state", encoding="utf-8")
    checkpoint = run_outputs.pack_output_tree(state, Path("run"))

    figures = tmp_path / "head" / "figs" / "icetune" / "run"
    figures.mkdir(parents=True)
    (figures / "summary.json").write_text(
        json.dumps({"trial_id": "matching-best"}),
        encoding="utf-8",
    )
    status = state / "status.json"
    status.write_text(json.dumps({"completed": 1}), encoding="utf-8")
    outputs = run_outputs.pack_runtime_outputs(
        figure_dir=figures,
        result_files=[],
        run_name="run",
        state_files=[status],
    )

    run_output = run_outputs.campaign_run_path(shared, "run", CAMPAIGN_FP)
    figure_output = run_outputs.campaign_figure_path(shared, "run", CAMPAIGN_FP)
    result_dir = run_output / "results"
    result_dir.mkdir(parents=True)
    result = result_dir / "TUNE_icetune_matching-best.pkl"
    result.write_bytes(b"matching-result")
    descriptor = {
        "filename": result.name,
        "sha256": sha256_file(result),
        "size": result.stat().st_size,
    }
    recovery = run_outputs.pack_recovery_generation(
        checkpoint=checkpoint,
        outputs=outputs,
        result_descriptors=[descriptor],
        run_name="run",
    )
    run_outputs.publish_recovery_generation(
        payload=recovery,
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )

    canonical = figure_output
    canonical.mkdir(parents=True)
    (canonical / "summary.json").write_text(
        json.dumps({"trial_id": "uncommitted-newer-best"}),
        encoding="utf-8",
    )
    newer = tmp_path / "newer" / "run"
    newer.mkdir(parents=True)
    (newer / "tuner.pkl").write_text("uncommitted-newer-state", encoding="utf-8")
    run_outputs.publish_checkpoint(
        payload=run_outputs.pack_output_tree(newer, Path("run")),
        shared_root=shared,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )

    runtime = tmp_path / "runtime"
    restored = run_outputs.stage_previous_outputs(
        shared_root=shared,
        runtime=runtime,
        run_name="run",
        campaign_fingerprint=CAMPAIGN_FP,
    )

    assert restored == [run_outputs.recovery_path(shared, "run", CAMPAIGN_FP)]
    assert (runtime / "runs" / "icetune" / "run" / "tuner.pkl").read_text(
        encoding="utf-8"
    ) == "matching-state"
    assert (
        json.loads((runtime / "figs" / "icetune" / "run" / "summary.json").read_text())["trial_id"]
        == "matching-best"
    )
    assert (
        runtime / "runs" / "icetune" / "run" / "results" / result.name
    ).read_bytes() == b"matching-result"
    assert json.loads((canonical / "summary.json").read_text())["trial_id"] == "matching-best"


# Check only the failed Condor process is requeued and original requirements survive
def test_ray_worker_replacement_is_scoped(monkeypatch):
    commands = []
    monkeypatch.setenv("RAY_CONDOR_COMMAND_TIMEOUT_S", "47")
    monkeypatch.setattr(ray_nodes, "condor_target", lambda: ("schedd", ["-addr", "<host:9618>"]))

    # Record allocation controls without contacting Condor
    def run(command, **kwargs):
        assert kwargs["timeout"] == 47
        commands.append(command)
        return SimpleNamespace(stdout='[{"Requirements": "TARGET.Memory >= 4096"}]')

    monkeypatch.setattr(subprocess, "run", run)
    ray_cluster.replace_ray_worker(marker="icetune_worker_12_3", machine="bad.cern.ch", exhausted=False)
    assert [command[0] for command in commands] == ["condor_q", "condor_hold", "condor_qedit", "condor_release"]
    assert all("12.3" in command for command in commands)
    assert "TARGET.Memory >= 4096" in commands[2][-1]
    assert 'TARGET.Machine =!= "bad.cern.ch"' in commands[2][-1]
    commands.clear()
    ray_cluster.replace_ray_worker(marker="icetune_worker_12_3", machine="bad.cern.ch", exhausted=True)
    assert [command[0] for command in commands] == ["condor_q", "condor_hold"]
    with pytest.raises(ValueError):
        ray_cluster.replace_ray_worker(marker="icetune_worker_12_0", machine="head.cern.ch", exhausted=False)


# Check interrupted replacement resumes without contacting the already held worker
@pytest.mark.parametrize("limit", [0, 1, 2])
def test_ray_worker_replacement_retry_limit(tmp_path, monkeypatch, limit):
    monkeypatch.setenv("RAY_WORKER_REPLACEMENT_LIMIT", str(limit))
    monkeypatch.setenv("RAY_WORKER_FAILURE_LIMIT", "3")
    marker = "icetune_worker_12_3"
    state, calls = {}, []
    detail = {"hostname": "bad.cern.ch", "job": {"job_id": "12.3"}, "logs": {}}
    client = SimpleNamespace(run=lambda *a, **k: {"worker": detail})
    fail = [True]

    # Simulate a failed qedit after Condor already held the worker
    def replace(**kwargs):
        calls.append(kwargs)
        if fail[0]:
            fail[0] = False
            raise OSError("qedit unavailable")

    monkeypatch.setattr(ray_cluster, "replace_ray_worker", replace)
    resources = {"resources": {"unhealthy_workers": {marker: {"node_id": "node1", "count": 2, "error": "ports"}}}}
    ray_cluster.recover_ray_workers(client=client, resources=resources, markers={"worker": marker}, recovery=state, log_dir=tmp_path)
    assert not calls
    resources["resources"]["unhealthy_workers"][marker]["count"] = 3
    ray_cluster.recover_ray_workers(client=client, resources=resources, markers={"worker": marker}, recovery=state, log_dir=tmp_path)
    assert not state[marker]["done"]
    ray_cluster.recover_ray_workers(client=None, resources=None, markers={}, recovery=state, log_dir=tmp_path)
    assert state[marker]["done"]
    assert state[marker]["attempts"] == 1
    for node in ("node1", "node2", "node3", "node4"):
        resources["resources"]["unhealthy_workers"][marker]["node_id"] = node
        ray_cluster.recover_ray_workers(client=client, resources=resources, markers={"worker": marker}, recovery=state, log_dir=tmp_path)
    assert len(calls) == limit + 2
    assert calls[-1]["exhausted"]
    assert len(list(tmp_path.glob("*.json"))) == limit + 1


# Check the allocation request and detected cgroup limit both constrain Ray memory
@pytest.mark.parametrize("detected,expected", [(2 * 1024**3, 2 * 1024**3), (128 * 1024**3, 4 * 1024**3)])
def test_ray_allocation_memory_budget(monkeypatch, detected, expected):
    monkeypatch.setitem(sys.modules, "distributed.system", SimpleNamespace(MEMORY_LIMIT=detected))
    monkeypatch.setitem(sys.modules, "dask.utils", SimpleNamespace(parse_bytes=lambda value: 4 * 1024**3))
    budget = ray_nodes.ray_memory("4GB")
    assert budget["memory"] == expected * 3 // 4
    assert budget["object_store_memory"] == expected // 8
    assert sum(budget.values()) < expected


# Check each campaign allocation has independent validated CPU and GPU counts
@pytest.mark.parametrize("role", ["init", "head", "worker"])
@pytest.mark.parametrize("resource,value", [("cpu", "3"), ("gpu", "2")])
def test_campaign_allocation_resources(monkeypatch, role, resource, value):
    monkeypatch.setattr(sys, "argv", ["submit_ray.py", "--campaign", "tune-gpom-res-con",
                                     f"--{role}-{resource}", value])
    assert getattr(submit.parse_args(), f"{role}_{resource}") == int(value)
    catalog = campaign_config.load_campaign_catalog(ICETUNE_TEST_DIR / "campaigns.yml")
    baseline = campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con",
                                             environment_overrides={"CPU_PER_TRIAL": "1"})
    key = f"RAY_{role.upper()}_{resource.upper()}"
    overrides = campaign_config.parse_overrides([f"ray.{role}_{resource}={value}"])
    selected = campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con",
                                             environment_overrides={"CPU_PER_TRIAL": "1", **overrides})
    assert selected[key] == value
    if role != "worker":
        assert selected["MAX_CONCURRENT"] == baseline["MAX_CONCURRENT"]
    for other in ("init", "head", "worker"):
        for count in ("cpu", "gpu"):
            other_key = f"RAY_{other.upper()}_{count.upper()}"
            if other_key != key:
                assert selected[other_key] == baseline[other_key]
    invalid = "0" if resource == "cpu" else "-1"
    with pytest.raises(ValueError, match=key):
        campaign_config.resolve(catalog, campaign_name="tune-gpom-res-con",
                                       environment_overrides={key: invalid})


# Check the real Dask Condor classes keep head and worker resource requests separate
@pytest.mark.parametrize("head_gpu,worker_gpu", [(0, 0), (1, 0), (0, 2), (1, 2)])
@pytest.mark.parametrize("backend", ["condor", "lxplus"])
def test_dask_condor_allocation_resources(head_gpu, worker_gpu, backend):
    if backend == "lxplus":
        job_cls = pytest.importorskip("dask_lxplus").CernCluster.job_cls
    else:
        job_cls = pytest.importorskip("dask_jobqueue.htcondor").HTCondorJob
    job = ray_cluster.cern_array_job(job_cls)(
        "tcp://127.0.0.1:8786", name="gpu-test", cores=8, memory="4GB", disk="4GB",
        python=sys.executable, worker_command="distributed.cli.dask_worker",
        processes=1, array_size=3, array_name="icetune-gpu-test",
        head_cpu=3, head_gpu=head_gpu, head_memory="8GB",
        worker_gpu=worker_gpu, worker_memory="4GB",
        job_extra_directives={"MY.SendCredential": "True", "should_transfer_files": "YES",
                              "transfer_output_files": '""'},
    )
    head, workers = job.job_script().split("Queue 1\n")
    assert re.search(r"(?im)^RequestCpus\s*=\s*3$", head)
    assert re.search(r"(?im)^RequestCpus\s*=\s*MY.DaskWorkerCores$", workers)
    assert job.job_header_dict["MY.DaskWorkerCores"] == 8
    assert re.search(rf"(?im)^request_gpus\s*=\s*{head_gpu}$", head)
    assert re.search(rf"(?im)^request_gpus\s*=\s*{worker_gpu}$", workers)
    assert workers.endswith("Queue 3\n")


# Execute serialized steering functions in an interpreter without the checkout on its import path
def test_lxplus_worker_serialize_without_checkout(tmp_path):
    import cloudpickle

    from submit.lxplus.main import register_worker_modules

    modules = [module for name, module in tuple(sys.modules.items())
               if name.startswith("submit") or name == "core.io.files"]
    register_worker_modules()
    try:
        payload = cloudpickle.dumps((ray_nodes.ray_host_slots, run_common.runtime_work_root,
                                     run_outputs.pack_output_tree, run_runtime.pack_runtime, run_main.request_runtime_stop,
                                     run_main.run_icetune_runtime, ray_nodes.check_worker_runtime))
    finally:
        for module in modules:
            cloudpickle.unregister_pickle_by_value(module)
    path = tmp_path / "steering.pkl"
    path.write_bytes(payload)
    source = tmp_path / "sample"
    source.mkdir()
    (source / "value.txt").write_text("worker output")
    script = '''import cloudpickle, io, pathlib, sys, tarfile
slots, scratch, pack, runtime, stop, run, check = cloudpickle.loads(pathlib.Path(sys.argv[1]).read_bytes())
assert slots(8) > 0
assert callable(run)
assert not stop(environment={"RUN_NAME": "isolated"})["requested"]
assert scratch({"RUN_NAME": "isolated"}).name == "icetune-ray-isolated"
source = pathlib.Path(sys.argv[2])
check({}, {"command": [sys.executable, "-c", "import ctypes; ctypes.CDLL(None)"], "env": {}, "timeout": 10})
try:
    check({str(source / "missing"): {}}, None)
except RuntimeError as exc:
    assert "No such file" in str(exc)
else:
    raise AssertionError("Unreadable dependencies were admitted")
with tarfile.open(fileobj=io.BytesIO(pack(source, pathlib.Path("sample"))), mode="r:gz") as archive:
    assert archive.extractfile("sample/value.txt").read() == b"worker output"
with tarfile.open(runtime(source), mode="r:gz") as archive:
    assert archive.extractfile("runtime/value.txt").read() == b"worker output"
'''
    subprocess.run([sys.executable, "-I", "-c", script, str(path), str(source)], cwd=tmp_path, check=True)


# Propagate the submission seed and reuse switch through the real CLI for every optimizer
@pytest.mark.parametrize("algorithm", sorted(campaign_config.SUPPORTED_ALGORITHMS))
def test_submit_optimizer_seed(submit_args, monkeypatch, algorithm):
    from core.icetune import parse_arguments

    submit_args.rngseed = 9187
    submit_args.no_ampfit_reuse = True
    submit_args.environment_overrides = [f"optimizer.algorithm={algorithm}"]
    environment = submit.resolve_campaign(submit_args)
    assert environment["RNGSEED"] == "9187"
    assert environment["AMPFIT_REUSE"] == "0"
    command = run_runtime.icetune_command(environment, cdir=ROOT, phase="preflight",
                                          tunesetup_path=submit_args.output_dir / "setup.json")
    monkeypatch.setattr(sys, "argv", ["icetune", *command[3:]])
    args = parse_arguments()
    assert args.algorithm == algorithm and args.rngseed == 9187


# Parse the shared launcher command with the real icetune CLI for every phase
@pytest.mark.parametrize("phase", ["preflight", "sample", "bank", "finalize", "init", "fit"])
@pytest.mark.parametrize("driver,campaign", [("GRANIITTI", "tune-gpom-res-con"), ("PANDORA", "tune-pandora-v0"),
                                            ("GRANIITTI", "tune-gpom-star-cms-ampfit")])
@pytest.mark.parametrize("address", ["local", "ray-head:6379"])
def test_shared_icetune_command(tmp_path, monkeypatch, phase, driver, campaign, address):
    from core.icetune import parse_arguments

    catalog = campaign_config.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
    environment = campaign_config.resolve(catalog, campaign_name=campaign, workers=1, rngseed=9187,
                                          environment_overrides={"AMPFIT_REUSE": "0"})
    if phase in {"sample", "bank", "finalize"} and environment["ALGORITHM"] != "ampfit":
        pytest.skip("Amplitude bank phases require ampfit")
    environment.update(RAY_ADDRESS=address, RAY_TMP_DIR=str(tmp_path / "ray"),
                       RAY_INIT_COORD_DIR=str(tmp_path / "coord"), RAY_INIT_BOOTSTRAP_CACHE=str(tmp_path / "cache"),
                       RAY_INIT_RUNTIME_SHA256="a" * 64, RAY_INIT_STATE_RELATIVE="runs/init.json", RAY_UPLOAD_RUNTIME="1",
                       RAY_RUNTIME_URI=(tmp_path / "runtime.tar.gz").as_uri())
    command = run_runtime.icetune_command(environment, cdir=ROOT, phase=phase,
                                          tunesetup_path=tmp_path / "setup.json", bank_shard=1,
                                          bank_index=2 if phase in {"sample", "bank", "finalize"} else None, bank_jobs=3)
    monkeypatch.setattr(sys, "argv", ["icetune", *command[3:]])
    args = parse_arguments()
    if phase in {"sample", "bank", "finalize", "init"}:
        assert args.coord_dir == environment["RAY_INIT_COORD_DIR"]
        position = sys.argv.index("--coord_dir")
        monkeypatch.setattr(sys, "argv", sys.argv[:position] + sys.argv[position + 2:])
        assert Path(parse_arguments().coord_dir) == ROOT / "runs/icetune" / args.run_name / "ray/init"
    assert args.simdriver == driver and args.backend == "ray"
    assert args.num_trials == int(environment["NUM_TRIALS"])
    assert args.rand_trials == int(environment["RAND_TRIALS"])
    assert args.rngseed == 9187
    assert not args.ampfit_reuse if environment["ALGORITHM"] == "ampfit" else args.ampfit_reuse
    if phase == "preflight":
        assert args.preflight and args.save_tunesetup == str(tmp_path / "setup.json")
    elif phase in {"sample", "bank", "finalize", "init"}:
        if phase != "init":
            assert args.phase == "bank" and args.bank_sample == (phase == "sample")
            assert args.bank_finalize == (phase == "finalize")
            assert args.bank_index == 2 and args.bank_jobs == 3
            assert args.bank_shard == (1 if phase == "bank" else None)
            assert args.bank_shared_dir == str(campaign_config.shared_output_dir(environment))
            return
        assert args.phase == "init" and args.cpu_per_trial == int(environment["RAY_INIT_CPU"])
        assert args.bootstrap_cache_url == environment["RAY_INIT_BOOTSTRAP_CACHE"]
    else:
        assert args.cpu_per_trial == int(environment["CPU_PER_TRIAL"]) and args.address == address
        assert args.ray_init_state == str(ROOT / "runs/init.json") and args.ray_upload_runtime
        assert args.ray_runtime_uri == environment["RAY_RUNTIME_URI"]
        assert args.runtime_sha256 == environment["RAY_INIT_RUNTIME_SHA256"]


# Verify Ray installs a prepared runtime locally without uploading it to GCS
def test_ray_prepared_runtime_archive(tmp_path):
    from ray._private.runtime_env.packaging import download_and_unpack_package
    from ray._private.runtime_env.working_dir import upload_working_dir_if_needed
    from ray.runtime_env import RuntimeEnv

    runtime = tmp_path / "runtime"
    bank = runtime / "runs/icetune/test/results/amplitude/0/nominal/amplitudes.bin"
    bank.parent.mkdir(parents=True)
    # A sparse bank exceeds the failed local package limit without costly MC generation
    with bank.open("wb") as handle:
        handle.seek(513 * 1024 * 1024)
        handle.write(b"bank")
    executable = runtime / "bin/python"
    executable.parent.mkdir()
    shutil.copy2(sys.executable, executable)
    packed = run_runtime.pack_runtime(runtime)
    directory = tmp_path / "shared"
    uri = run_runtime.publish_ray_runtime(archive=packed, directory=directory)
    assert run_runtime.publish_ray_runtime(archive=packed, directory=directory) == uri
    env = RuntimeEnv(working_dir=uri)
    assert upload_working_dir_if_needed(env, include_gitignore=False,
                                       scratch_dir=str(tmp_path / "upload"))["working_dir"] == uri
    for node in ("head", "worker"):
        cache = tmp_path / node
        cache.mkdir()
        installed = Path(asyncio.run(download_and_unpack_package(uri, str(cache))))
        copied = installed / bank.relative_to(runtime)
        assert copied.stat().st_size == bank.stat().st_size
        with copied.open("rb") as handle:
            handle.seek(-4, 2)
            assert handle.read() == b"bank"
        assert (installed / "bin/python").stat().st_mode & 0o111


# Check EOS result publication needs no hard links and rejects conflicting trials
def test_lxplus_results_without_hard_links(tmp_path, monkeypatch):
    import errno

    # Reproduce the EOS hard link failure while using the real output publisher
    def no_link(*args, **kwargs):
        raise OSError(errno.EXDEV, "Invalid cross-device link")

    monkeypatch.setattr(run_outputs.os, "link", no_link)
    source, target = tmp_path / "source", tmp_path / "target"
    source.mkdir()
    trial = source / "TUNE_icetune_trial.pkl"
    trial.write_bytes(b"completed trial")
    (source / "data_covariance.json").write_text('{"bins": 2}')
    run_outputs.publish_result_files(source=source, target=target)
    run_outputs.publish_result_files(source=source, target=target)
    assert (target / trial.name).read_bytes() == trial.read_bytes()
    assert (target / "data_covariance.json").read_bytes() == (source / "data_covariance.json").read_bytes()
    trial.write_bytes(b"different trial")
    with pytest.raises(RuntimeError, match="Conflicting"):
        run_outputs.publish_result_files(source=source, target=target)
    assert (target / trial.name).read_bytes() == b"completed trial"


# Check failed worker diagnostics also find the session path inherited from the head
def test_ray_worker_inherited_session_logs(tmp_path):
    root = tmp_path / "icetune_head/ray/session_test/logs"
    root.mkdir(parents=True)
    (root / "raylet.err").write_text("Address already in use")
    detail = ray_nodes.ray_startup_detail(["ray", "start"], tmp_path / "icetune_worker_1_1/ray", "", "")
    assert "Address already in use" in detail


# Reuse compact runtime archives only when their complete checksum still matches
def test_compact_runtime_integrity(tmp_path):
    source = tmp_path / "runtime"
    source.mkdir()
    (source / "input.txt").write_text("runtime input\n")
    packed = run_runtime.pack_runtime(source)
    target = tmp_path / "published"
    uri = run_runtime.publish_ray_runtime(archive=packed, directory=target)
    archive, = target.glob("*.tar.gz")
    digest = sha256_file(archive)
    assert archive.name == f"{digest[:16]}.tar.gz"
    assert run_runtime.publish_ray_runtime(archive=packed, directory=target) == uri
    with tarfile.open(archive) as contents:
        assert contents.extractfile("runtime/input.txt").read() == (source / "input.txt").read_bytes()
    archive.write_bytes(b"conflicting runtime")
    with pytest.raises(ValueError, match="checksum mismatch"):
        run_runtime.publish_ray_runtime(archive=packed, directory=target)
    assert archive.read_bytes() == b"conflicting runtime"


# Pack and restore a real runtime below initially missing submission folders
def test_lxplus_archive_missing_parent(tmp_path):
    source = tmp_path / "runtime"
    source.mkdir()
    (source / "input.txt").write_text("runtime input\n")
    archive = tmp_path / "runs" / "submit" / "runtime.tar.zst"
    checksum = run_runtime.pack_lxplus_runtime(runtime=source, output=archive)
    restored = run_runtime.extract_lxplus_runtime(
        archive=archive, target=tmp_path / "worker" / "runtime", checksum=checksum,
    )
    assert (restored / "input.txt").read_bytes() == (source / "input.txt").read_bytes()


# Transfer standalone L-BFGS outputs without a Ray Tune checkpoint
def test_lbfgs_output_transfer(tmp_path, monkeypatch):
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path))
    environment = {"RUN_NAME": "fit", "ALGORITHM": "lbfgs"}
    state = run_common.runtime_work_root(environment) / "runtime/runs/icetune/fit"
    (state / "results").mkdir(parents=True)
    (state / "history.json").write_text('{"trials": []}')
    (state / "summary.json").write_text('{"algorithm": "lbfgs"}')
    result = state / "results/TUNE_icetune_trial.pkl"
    result.write_bytes(b"trial output")
    snapshot = run_outputs.snapshot_final_outputs(environment=environment)
    assert snapshot["checkpoint"] is None
    assert snapshot["result_names"] == [result.name]
    shared = tmp_path / "shared"
    run_outputs.restore_outputs(payload=snapshot["outputs"], shared_root=shared,
                                run_name="fit", campaign_fingerprint=CAMPAIGN_FP)
    target = run_outputs.campaign_run_path(shared, "fit", CAMPAIGN_FP)
    for relative in ("history.json", "summary.json", "results/" + result.name):
        assert (target / relative).read_bytes() == (state / relative).read_bytes()
