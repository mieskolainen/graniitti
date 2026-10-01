# Cached Ray event partitions against native amplitudes and exact derivatives
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import time

import pytest
import ray
import torch
from core import resource
from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank, bank_files, prepare_bank, stage_bank
from core.tune.drivers.graniitti.ampfit.distributed import Partition, remote_amplitude, start_workers
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.drivers.graniitti.tunesetup import card
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.parameters.space import normalize_continuous_param_space

from tests.technical.integration.test_ampfit import ROOT, dataset_card, kaon_events  # noqa: F401


# Prepare physical GP amplitudes with a nonlinear relative phase
@pytest.fixture
def bank(tmp_path, request):
    events = request.getfixturevalue("kaon_events")
    driver = GraniittiDriver()
    dataset = card(dataset_card(tmp_path, "STAR_1792394/KK"), nevents=16, loopscreen=False,
                   xsmode="sample", swap_process="GP[RES+CON]<F>", integrator="VEGAS")
    driver.init_data(str(tmp_path / "run"), [dataset], "default", str(ROOT), pickle_dump=False)
    names = ["RES|f0_980:GP:FF_prod.Lambda2", "RES|f0_980:GP:g_ls(0,0)@NORM", "RES|f0_980:GP:phi@PHASE"]
    initial = driver.get_initial_param(dict.fromkeys(names, {}), {}, str(ROOT))
    space = {name: {"type": "uniform", "lower": value - abs(value) / 2 - 0.1,
                    "upper": value + abs(value) / 2 + 0.1} for name, value in initial.items()}
    run = driver._run_cards(index=0, datacard=dataset, mc_steer={"tune_default": "TUNE0", "tunesetup_name": tmp_path.name},
                           tunename="TUNE0", cdir=str(ROOT))[0]
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"] | {"batch_events": 3}
    directory = tmp_path / "banks/0/nominal"
    prepare_bank(driver=driver, directory=directory, temporary=tmp_path / "basis", run=run,
                 tune=ROOT / "modeldata/TUNE0", initial=initial, bounds=normalize_continuous_param_space(space),
                 pid=driver.pid[0][0], events=events, controls=controls, cdir=str(ROOT),
                 deadline=time.monotonic() + 120)
    native = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                           cuts=driver.cuts[0], controls=controls)
    return native, directory, driver, controls


# Compare complex outputs, gradients and Hessians after removing access to the shared source
def test_cached_partitions(bank):
    native, directory, driver, controls = bank
    ray.init(num_cpus=4, include_dashboard=False, _temp_dir=str(ROOT / "tmp/r"))
    try:
        actors = [ray.remote(Partition).options(num_cpus=1).remote(
            [("0/nominal", str(directory), start, stop)], 1, 3) for start, stop in ((0, 7), (7, 16))]
        ray.get([actor.ready.remote() for actor in actors])
        workers = [(actor, "0/nominal", start, stop) for actor, (start, stop) in zip(actors, ((0, 7), (7, 16)), strict=True)]
        names = sorted(native.initial)
        theta = torch.tensor([native.initial[name] for name in names], dtype=torch.float64, requires_grad=True)
        parameters = dict(zip(names, theta, strict=True))
        expected = native.amplitude(parameters)
        actual = remote_amplitude(parameters, workers)
        torch.testing.assert_close(actual, expected, rtol=1e-12, atol=1e-12)
        weight = torch.arange(expected.numel(), dtype=torch.float64).reshape(expected.shape) + 1
        first, = torch.autograd.grad((expected.abs().square() * weight).sum(), theta, create_graph=True)
        second, = torch.autograd.grad((actual.abs().square() * weight).sum(), theta, create_graph=True)
        torch.testing.assert_close(second, first, rtol=1e-10, atol=1e-10)
        torch.testing.assert_close(torch.autograd.grad(second.sum(), theta)[0],
                                  torch.autograd.grad(first.sum(), theta)[0], rtol=1e-9, atol=1e-9)
        # Histogram likelihood and curvature use the same contracted complex amplitudes
        def objective(values, workers):
            native.workers = workers
            predictions = native.predict(dict(zip(names, values, strict=True)))
            results = dict(mc=[predictions], data=driver.data, obs=driver.obs, datasets=driver.datasets)
            param = dict(data_covariance_mode="diagonal", cost="chi2", cost_avg="sum", cost_rho="quadratic",
                         mc_steer={"ampfit": controls})
            return driver.trial_costs(results, param)["metrics"]["chi2"]

        local_hessian = torch.autograd.functional.hessian(lambda values: objective(values, None), theta)
        remote_hessian = torch.autograd.functional.hessian(lambda values: objective(values, workers), theta)
        torch.testing.assert_close(remote_hessian, local_hessian, rtol=1e-9, atol=1e-9)
        saved = directory.with_name("nominal._old")
        directory.rename(saved)
        try:
            # Persistent partitions continue to evaluate without any EOS reads
            torch.testing.assert_close(remote_amplitude(parameters, workers), expected, rtol=1e-12, atol=1e-12)
        finally:
            saved.rename(directory)
        for actor in actors:
            ray.kill(actor)
        driver.amplitude_path = str(directory.parents[1])
        pool = start_workers(driver, {"cdir": str(ROOT), "run_name": "unused", "max_t": 60,
                                      "mc_steer": {"ampfit": controls}})
        torch.testing.assert_close(remote_amplitude(parameters, pool["0/nominal"]), expected, rtol=1e-12, atol=1e-12)
        actor = pool["0/nominal"][0][0]
        ray.kill(actor, no_restart=False)
        ray.get(actor.ready.remote(), timeout=60)
        torch.testing.assert_close(remote_amplitude(parameters, pool["0/nominal"]), expected, rtol=1e-12, atol=1e-12)
    finally:
        ray.shutdown()



# Round trip a worker bootstrap containing references instead of native bank arrays
def test_bank_bootstrap(bank, tmp_path):
    from core.tune import cache
    from core.tune.drivers.graniitti.runtime import allowed_bootstrap_path, stage_bootstrap

    native, directory, driver, controls = bank
    root = tmp_path / "runtime"
    relative = "runs/icetune/fit/results/amplitude"
    stage_bank(directory, root / relative / "0/nominal", predictions=False)
    files = bank_files(root / relative)
    assert not any(path.suffix == ".bin" or "parts" in path.parts for path in files)
    records = cache.file_records(files, root=root, allowed=allowed_bootstrap_path)
    temporary = tmp_path / "archive"
    temporary.mkdir()
    bootstrap = cache.publish_content_archive(cache_base_url=str(tmp_path / "cache"), kind="graniitti",
        fingerprint=cache.json_fingerprint(records), root=root, records=records, temporary_dir=temporary,
        schema_version=driver.BOOTSTRAP_SCHEMA_VERSION, label="GRANIITTI")
    worker = tmp_path / "worker"
    stage_bootstrap(cdir=worker, bootstrap=bootstrap)
    restored = AmplitudeBank(driver=driver, directory=worker / relative / "0/nominal", obs=driver.obs[0],
                             pid=driver.pid[0], cuts=driver.cuts[0], controls=controls)
    torch.testing.assert_close(restored.amplitude(native.initial), native.amplitude(native.initial), rtol=0, atol=0)
    assert not list(worker.rglob("*.bin"))
