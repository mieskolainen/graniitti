# Pandora batch initialization and steering interpreter checks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

from core import resource
from core.tune.drivers.pandora import runtime
from core.tune.drivers.pandora.driver import PandoraDriver
from core.tune.optimizers.hebo.config import load_hebo_config
from core.tune.tunesetup import load_tunesetup, save_tunesetup

from submit import CAMPAIGN_DIR, campaign, campaign_source, schedulers
from submit import runtime as submission

ROOT = Path(__file__).resolve().parents[3]


# Validate INIT metadata without Conda fitting and parameter editing packages
# STEER runs in the separate CERN LCG interpreter
def test_steering_imports():
    command = """
from core.tune.drivers.pandora.runtime import stage_bootstrap
stage_bootstrap(cdir='.', bootstrap=dict(archive_sha256=None, archive_size=0, archive_url=None, files=[]))
try:
    stage_bootstrap(cdir='.', bootstrap=dict(archive_sha256='unexpected', archive_size=0, archive_url=None, files=[]))
except RuntimeError:
    pass
else:
    raise AssertionError('Pandora accepted generator caches')
"""
    environment = dict(os.environ, PYTHONPATH=str(resource("").parent))
    result = subprocess.run([sys.executable, "-S", "-c", command], env=environment, capture_output=True, text=True)
    assert result.returncode == 0, result.stdout + result.stderr


# Run the real transferred INIT shell and load its published point on the fit head
# Reconstruction files are frozen as archive inputs but no reconstruction is requested
def test_batch_init(tmp_path, pandora_inputs):
    environment = campaign.resolve(campaign.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml"),
        campaign_name="tune-pandora-v0", run_name="PANDORA_V0", repo_dir=ROOT)
    tunesetup = load_tunesetup(cdir=ROOT, simdriver="PANDORA", name=campaign_source("tune-pandora-v0"),
                              tune_default=environment["TUNE_DEFAULT"])
    root = Path(tunesetup.datacards[0]["pandora_dir"])
    checkout = tmp_path / "checkout"
    for relative in submission.BASE_PATHS:
        entry = checkout / relative
        entry.parent.mkdir(parents=True, exist_ok=True)
        entry.symlink_to(ROOT / relative)
    source = submission.stage_runtime(source=checkout, target=tmp_path / "runtime")
    relative = "tmp/icetune/pandora"
    manifest, archive = runtime.freeze_runtime(
        root, ["run/PandoraSettings_v6.xml", "run/run_reco_pandora.py"], source / relative,
        {"setup": "", "dependencies": {}},
    )
    manifest.update(archive=archive.relative_to(source).as_posix(), driver_files={
        name: runtime.file_identity(Path(runtime.__file__).parent / name)
        for name in ("driver.py", "runtime.py", "tunesetup/common.py")
    })
    for card in tunesetup.datacards:
        card["runtime"] = manifest
    tunesetup.aux_param_space["runtime_files"] = {manifest["archive"]: str(archive)}
    tunesetup.optimizer = {"hebo": load_hebo_config(resource("tune/settings/hebo.json"))}
    definition = source / "tmp/icetune/tunesetup.json"
    save_tunesetup(tunesetup=tunesetup, path=definition, simdriver="PANDORA", tune_default=tunesetup.tune_default)
    batch_archive = tmp_path / "runtime.tar.zst"
    checksum = submission.pack_lxplus_runtime(runtime=source, output=batch_archive)
    environment.update(CONDA_EXE=os.environ["CONDA_EXE"], ICETUNE_CONDA_PREFIX=sys.prefix,
                       TUNESETUP="tmp/icetune/tunesetup.json", HEBO_CONFIG="tmp/icetune/tunesetup.json#/optimizer/hebo")
    coord = tmp_path / "coord"
    rendered = schedulers._ray_init_environment(
        environment=environment, runtime_archive=batch_archive, runtime_sha256=checksum,
        bootstrap_cache=tmp_path / "cache", coord_dir=coord,
    )
    job_env = {key: os.environ[key] for key in ("HOME", "USER", "LOGNAME", "PATH") if key in os.environ}
    job_env.update(dict(item.split("=", 1) for item in shlex.split(rendered)), _CONDOR_SCRATCH_DIR=str(tmp_path))
    result = subprocess.run(["bash", str(ROOT / "submit/shell/steer_ray_init_lxplus.sh")],
        cwd=tmp_path, env=job_env, capture_output=True, text=True, timeout=120)
    assert result.returncode == 0, result.stdout + result.stderr
    state = json.loads((coord / "init.json").read_text())
    assert state["status"] == "completed" and state["bootstrap"]["kind"] == "pandora"
    driver = PandoraDriver()
    expected = driver.bootstrap_payload(runtime_sha256=checksum, datacards=tunesetup.datacards)
    assert state["bootstrap"] == expected
    assert state["initial_points"] == driver.get_initial_param(
        tunesetup.param_space, tunesetup.aux_param_space, str(source), tunesetup.tune_default)
    submission.stage_ray_bootstrap(runtime=source, coord_dir=coord, run_name="PANDORA_V0", simdriver="PANDORA")
    args = SimpleNamespace(run_name="PANDORA_V0", runtime_sha256=checksum,
                           ray_init_state=source / "runs/icetune/PANDORA_V0/ray_init.json")
    assert driver.initialize_ray_bootstrap(args=args, tunesetup=tunesetup, mc_steer={}) == state["initial_points"]
    assert not list(tmp_path.rglob("reco_*.root"))
    assert not list(tmp_path.rglob("tuner.pkl"))
