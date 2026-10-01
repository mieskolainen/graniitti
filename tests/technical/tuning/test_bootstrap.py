# Tests for portable icetune bootstrap identities and mismatch diagnostics
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import json
import shutil
from pathlib import Path
from types import SimpleNamespace

import pytest
from core import resource
from core.tune import cache, core
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.optimizers.ampfit.config import load_settings

from submit.lxplus import common, outputs

ROOT = Path(__file__).resolve().parents[3]


# Read real HEPData and native amplitude inputs in an isolated checkout
@pytest.fixture
def bootstrap_case(tmp_path, data_card):
    shutil.copytree(ROOT / "modeldata/TUNE0", tmp_path / "modeldata/TUNE0", dirs_exist_ok=True)
    (tmp_path / "bin").mkdir()
    for name in ("gr", "ampfit"):
        shutil.copy2(ROOT / "bin" / name, tmp_path / "bin" / name)
    shutil.copy2(ROOT / "VERSION.json", tmp_path / "VERSION.json")
    driver = GraniittiDriver()
    args = SimpleNamespace(
        tune_default="TUNE0", tunesetup="tunesetup.json", algorithm="ampfit",
        ampfit_settings=load_settings(resource("tune/settings/ampfit.json")), ampfit_reuse=True,
        bank_shared_dir="/eos/shared/banks", data_covariance_mode="diagonal",
        mc_correlation_events=None, mc_correlation_weighting="card",
    )
    options = dict(cdir=str(tmp_path), datacards=[{"datacard": data_card.relative_to(tmp_path).as_posix()}],
                   mc_steer=driver.build_run_steering(args), runtime_sha256="a" * 64, tunesetup_path=None)
    driver.init_data(run_name="init", datacards=options["datacards"], obs_module="default",
                     cdir=options["cdir"], pickle_dump=False)
    return driver, args, options


# Preserve INIT identity when fit storage is absent, changed or omitted from steering
@pytest.mark.parametrize("location", [None, "/other/shared/banks", "omitted"])
def test_bank_storage_identity(bootstrap_case, location):
    driver, args, options = bootstrap_case
    original = copy.deepcopy(options)
    expected = driver.bootstrap_fingerprint(**options)
    args.bank_shared_dir = location
    steering = driver.build_run_steering(args)
    if location == "omitted":
        steering.pop("ampfit_reuse_root")
    assert driver.bootstrap_fingerprint(**(options | {"mc_steer": steering})) == expected
    assert options == original


# Preserve the identity after copying native inputs into a different worker checkout
# Sort package and checkout records by their portable names, regardless of absolute paths
@pytest.mark.parametrize("changed", [False, True])
def test_staged_inputs(bootstrap_case, tmp_path, changed):
    driver, _, options = bootstrap_case
    inputs = driver.bootstrap_inputs(**options)
    staged = tmp_path / "worker"
    for name in ("bin", "modeldata", "icepack"):
        shutil.copytree(tmp_path / name, staged / name)
    shutil.copy2(tmp_path / "VERSION.json", staged / "VERSION.json")
    restored = GraniittiDriver()
    restored.init_data(run_name="fit", datacards=options["datacards"], obs_module="default",
                       cdir=str(staged), pickle_dump=False)
    if changed:
        with (staged / "bin/ampfit").open("ab") as stream:
            stream.write(b"changed native amplitude binary")
    actual = restored.bootstrap_inputs(**(options | {"cdir": str(staged)}))
    assert (cache.json_fingerprint(inputs) == cache.json_fingerprint(actual)) is not changed
    assert actual["amplitude_files"] == sorted(actual["amplitude_files"], key=lambda record: record["path"])


# Retain checks of physics steering, dataset event counts and the submitted runtime
@pytest.mark.parametrize("changed", ["bank", "covariance", "events", "runtime"])
def test_physics_identity(bootstrap_case, changed):
    driver, _, options = bootstrap_case
    expected = driver.bootstrap_fingerprint(**options)
    varied = copy.deepcopy(options)
    if changed == "bank":
        varied["mc_steer"]["ampfit"]["nodes"] += 1
    elif changed == "covariance":
        varied["mc_steer"]["data_covariance_mode"] = "full"
    elif changed == "events":
        varied["datacards"][0]["nevents"] = 16
    else:
        varied["runtime_sha256"] = "b" * 64
    assert driver.bootstrap_fingerprint(**varied) != expected


# Publish the hashed inputs and retain both sides when a staged runtime is incompatible
# Exercise the real archive writer and INIT reader without generating fit trials
@pytest.mark.parametrize("reuse", [True, False])
@pytest.mark.parametrize("writable", [True, False])
def test_bootstrap_diagnostics(bootstrap_case, tmp_path, monkeypatch, writable, reuse):
    driver, _, options = bootstrap_case
    options["mc_steer"]["ampfit_reuse"] = reuse
    output = tmp_path / "vgrid/grid.dat"
    output.parent.mkdir()
    bootstrap = driver.prepare_bootstrap(
        **options, cache_base_url=str(tmp_path / "cache"), init_force=False,
        initialize=lambda: output.write_bytes(b"cache transport test"), reusable_vgrids=lambda: [output],
    )
    previous = bootstrap
    bootstrap = driver.prepare_bootstrap(
        **options, cache_base_url=str(tmp_path / "cache"), init_force=False,
        initialize=lambda: output.write_bytes(b"rebuilt cache transport test"), reusable_vgrids=lambda: [output],
    )
    assert (bootstrap == previous) is reuse
    assert Path(previous["archive_url"]).is_file()
    assert cache.json_fingerprint(bootstrap["fingerprint_inputs"]) == bootstrap["fingerprint"]
    environment = {"RUN_NAME": "fit"}
    monkeypatch.setenv("_CONDOR_SCRATCH_DIR", str(tmp_path / "condor"))
    state_path = common.runtime_work_root(environment) / "runtime/runs/icetune/fit/ray_init.json"
    state_path.parent.mkdir(parents=True)
    state_path.write_text(json.dumps({
        "backend": "ray", "protocol": "core.tune.ray.bootstrap", "protocol_version": 1,
        "status": "completed", "bootstrap": bootstrap, "initial_points": {"parameter": 0.42},
    }))
    loaded, points = core.load_ray_init_state(path=state_path, fingerprint=bootstrap["fingerprint"])
    assert loaded == bootstrap and points == {"parameter": 0.42}
    runtime_inputs = driver.bootstrap_inputs(**(options | {"runtime_sha256": "b" * 64}))
    runtime_fingerprint = cache.json_fingerprint(runtime_inputs)
    failures = state_path.parent / "failures"
    if not writable:
        failures.write_text("not a directory")
    with pytest.raises(cache.PermanentConfigurationError) as caught:
        core.load_ray_init_state(path=state_path, fingerprint=runtime_fingerprint, fingerprint_inputs=runtime_inputs)
    assert bootstrap["fingerprint"] in str(caught.value) and runtime_fingerprint in str(caught.value)
    if writable:
        comparison, = failures.glob("ray_init_*.json")
        report = json.loads(comparison.read_text())
        assert report["init_inputs"] == bootstrap["fingerprint_inputs"]
        assert report["runtime_inputs"] == runtime_inputs
        snapshot = outputs.snapshot_final_outputs(environment=environment)
        assert snapshot["checkpoint"] is None
        shared = tmp_path / "shared"
        outputs.restore_outputs(payload=snapshot["outputs"], shared_root=shared, run_name="fit", campaign_fingerprint="c" * 64)
        published = outputs.campaign_run_path(shared, "fit", "c" * 64) / "failures" / comparison.name
        assert json.loads(published.read_text()) == report
    assert json.loads(state_path.read_text())["bootstrap"] == bootstrap
