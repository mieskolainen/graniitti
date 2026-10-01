# Native amplitude preparation and coherent reweighting through the real icetune driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import math
import multiprocessing
import subprocess
import time
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyjson5
import pytest
import torch
from core import resource
from core.io import readers
from core.io.serialize import load_json_file, write_json_file
from core.plot import plot
from core.tune.drivers.graniitti.ampfit.amplitude import (
    AmplitudeBank,
    BankPlan,
    finalize_bank,
    input_card,
    prepare_bank,
)
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.drivers.graniitti.tunesetup import card
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.parameters.space import normalize_continuous_param_space

from submit import campaign_source

ROOT = Path(__file__).resolve().parents[3]


# Retain the original HEPData and fiducial cuts while reducing VEGAS sampling for amplitude closure
def dataset_card(directory, dataset_name, resonance="f0_980"):
    source = ROOT / "icepack/SOFTCEP" / dataset_name / "dataset.json"
    dataset = load_json_file(source, loader=pyjson5.load)
    dataset["reader"] = str((source.parent / dataset["reader"]).resolve())
    for entry in dataset["sets"]:
        for field in ("cuts", "obs"):
            entry[field] = str((source.parent / entry[field]).resolve())
    sample = dataset["samples"][0]
    gencard = load_json_file((source.parent / sample["gencard"]).resolve(), loader=pyjson5.load)
    gencard["GENERIC"]["OUTPUT"] = f"{directory.parent.name}_{directory.name}"
    if dataset_name.endswith("/KK"):
        # Isolate the selected pole and sample physical kaon pairs above threshold
        gencard["SCATTERING"]["RES"] = [resonance]
        gencard["GENCUTS"]["<F>"]["M"][0] *= 1.2
    sample["gencard"] = str(directory / "gencard.json")
    write_json_file(sample["gencard"], gencard)
    sample["parameters"].update({
        "GENERIC.WEIGHTED": True, "INTEGRATOR.min_samples": 2048,
        "INTEGRATOR.max_samples": 8192, "INTEGRATOR.precision": 0.9,
        "INTEGRATOR.VEGAS.ncall": 1000, "INTEGRATOR.VEGAS.rounds": 2,
    })
    path = directory / "dataset.json"
    write_json_file(path, dataset)
    return str(path)


# Generate one common kaon phase space sample for fixed event comparisons across models
@pytest.fixture(scope="module")
def kaon_events(tmp_path_factory):
    directory = tmp_path_factory.mktemp("ampfit_kaon_events")
    driver = GraniittiDriver()
    datacards = [card(dataset_card(directory, "STAR_1792394/KK"), nevents=16, loopscreen=False,
                     xsmode="sample", swap_process="MP[RES+CON]<F>", integrator="VEGAS")]
    steering = dict(tune_default="TUNE0", tunesetup_name=directory.name)
    driver.init_data(str(directory / "run"), datacards, "default", str(ROOT), pickle_dump=False)
    run = driver._run_cards(index=0, datacard=datacards[0], mc_steer=steering, tunename="TUNE0", cdir=str(ROOT))[0]
    command = [str(ROOT / "bin/gr"), "-i", run.inputfile, "-o", run.output, "-m", "TUNE0",
               "-n", "-1", "-w", "true", "-l", "false", "-h", "0", "-c", "1",
               "--set", "NUMERICS.json:NUMERICS_VEGAS.automatic_convergence=false"]
    driver._sample_command(list(command), run=run, cdir=str(ROOT), deadline=time.monotonic() + 600,
                           stage="vgrid_initialization", dataset_index=0, output=run.output)
    output = run.output + "_bank"
    command[command.index("-n") + 1] = str(datacards[0]["nevents"])
    command[command.index("-o") + 1] = output
    command.extend(("-d", run.gridfile))
    driver._sample_command(command, run=run, cdir=str(ROOT), deadline=time.monotonic() + 600,
                           stage="event_generation", dataset_index=0, output=output)
    return ROOT / "output" / f"{output}.hepmc3"


# Copy a native bank during submission preflight only when reuse is enabled
def test_ampfit_preflight_reuse(tmp_path, kaon_events):
    driver = GraniittiDriver()
    space = {"RES|f0_980:GP:phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi}}
    tunesetup = SimpleNamespace(param_space=space, aux_param_space={}, datacards=[
        card(dataset_card(tmp_path, "STAR_1792394/KK"), nevents=16, loopscreen=False,
             xsmode="sample", swap_process="GP[RES+CON]<F>", integrator="VEGAS")])
    args = SimpleNamespace(algorithm="ampfit", cost="chi2", cdir=str(ROOT), tune_default="TUNE0", rngseed=73,
                           run_name=f"bank-preflight-{tmp_path.parent.name}", tunesetup="bank_preflight", obs_module="default",
                           data_covariance_mode="diagonal", mc_correlation_events=None, mc_correlation_weighting="sample",
                           ampfit_settings=load_settings(resource("tune/settings/ampfit.json")),
                           preflight=True, ampfit_reuse=False, bank_shared_dir=str(tmp_path))
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    source = ROOT / "runs/icetune" / args.run_name / "init/ampfit_preflight"
    previous = tmp_path / "runs/icetune/previous/results/amplitude/0/nominal"
    run = driver._run_cards(index=0, datacard=tunesetup.datacards[0], mc_steer=driver.build_run_steering(args),
                            tunename="TUNE0", cdir=str(ROOT))[0]
    initial = driver.get_initial_param(space, {}, str(ROOT))
    prepare_bank(driver=driver, directory=previous, temporary=tmp_path / "basis", run=run, tune=source,
                 initial=initial, bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                 events=kaon_events, controls=args.ampfit_settings["bank"], cdir=str(ROOT), deadline=time.monotonic() + 120)
    # Recover finished shards from the failed full-bank merge without changing native bytes
    metadata = previous / "amplitudes.bin.json"
    metadata.rename(metadata.with_suffix(".json._old"))
    partial = previous / "amplitudes.bin.partial"
    partial.write_bytes(b"interrupted merge")
    plan = BankPlan(driver=driver, run=run, tune=source, initial=initial,
                    bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                    controls=args.ampfit_settings["bank"], cdir=str(ROOT), nevents=16)
    plan.recover(previous)
    assert not plan.mismatch(previous)
    assert partial.read_bytes() == b"interrupted merge"
    assert not (previous / "amplitudes.bin").exists()
    finalize_bank(directory=previous, temporary=tmp_path / "finalize", parameters=initial, covariance=False, ready=True,
                  options=dict(driver=driver, obs=driver.obs[0], pid=driver.pid[0], cuts=driver.cuts[0],
                               controls=args.ampfit_settings["bank"]))
    target = driver.amplitude_directory(cdir=str(tmp_path), run_name=args.run_name) / "0/nominal"
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    assert not target.exists()
    args.ampfit_reuse = True
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    banks = [AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                           cuts=driver.cuts[0], controls=args.ampfit_settings["bank"]) for directory in (previous, target)]
    torch.testing.assert_close(banks[0].coherent(banks[0].coefficients(initial)), banks[1].coherent(banks[1].coefficients(initial)))
    published = {path: path.stat().st_mtime_ns for path in target.glob("amplitudes.bin*")}
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    assert published == {path: path.stat().st_mtime_ns for path in published}
    # Report changed physics even when the requested event count also differs
    key = next(iter(space))
    plan = BankPlan(driver=driver, run=run, tune=source, initial={**initial, key: initial[key] + 0.1},
                    bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                    controls=args.ampfit_settings["bank"], cdir=str(ROOT), nevents=32)
    assert set(plan.mismatch(target)) == {"event count 16 != 32", "settings.initial"}
    tunesetup.datacards[0]["nevents"] = 32
    driver.prepare_tunesetup(tunesetup=tunesetup, args=args)
    assert not (target / "finalized.pkl").exists()
    assert list(target.glob("finalized.*.pkl._old"))
    assert (previous / "finalized.pkl").is_file()
    # Disabling reuse also replaces banks already present in the current run
    args.ampfit_reuse = False
    tunesetup.datacards[0]["nevents"] = 16
    driver.initialize(run_name=args.run_name, tunesetup=tunesetup, mc_steer=driver.build_run_steering(args),
                      obs_module="default", cdir=str(ROOT), init_force=False, pickle_dump=False, max_t=600,
                      bank_mode={"phase": "sample", "jobs": 1, "shared_root": str(tmp_path)})
    assert not (target / "amplitudes.bin").exists()
    assert load_json_file(target / "source.bin.json")["shape"][1] == 16



# Compare complete init preparation, native varied amplitudes and the shared icepack histograms
@pytest.mark.parametrize("resonance", ["f0_980", "f2_1950"])
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
@pytest.mark.parametrize("screened", [False, True])
@pytest.mark.parametrize("dataset_name", ["CMS_2752118/pipi_0p7", "STAR_1792394/KK"])
def test_ampfit_init_and_reweighting(tmp_path, request, model, screened, dataset_name, resonance):
    driver = GraniittiDriver()
    key = f"RES|{resonance}:{model}:phi"
    space = {key: {"type": "uniform", "lower": -math.pi, "upper": math.pi}}
    pole = load_json_file(ROOT / f"modeldata/TUNE0/RES/{resonance}.json", loader=pyjson5.load)["PARAM_RES"]["MODELS"][model]
    space.update({f"RES|{resonance}:{model}:{field}": {"type": "uniform", "lower": 0.8 * pole[field], "upper": 1.2 * pole[field]}
                  for field in ("mass", "width")})
    tunesetup = SimpleNamespace(
        param_space=space, aux_param_space={},
        datacards=[card(dataset_card(tmp_path, dataset_name, resonance), nevents=16, loopscreen=screened, xsmode="sample",
                        swap_process=f"{model}[RES+CON]<F>", integrator="VEGAS")],
    )
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"]
    steering = dict(tune_default="TUNE0", tunesetup_name=tmp_path.name, data_covariance_mode="diagonal", ampfit=controls)
    directory = driver.amplitude_directory(cdir=str(ROOT), run_name=str(tmp_path / "run")) / "0/nominal"
    kaons = dataset_name.endswith("/KK")
    if kaons:
        # Compare physical amplitudes at common four-momenta independently of the proposal model
        events = request.getfixturevalue("kaon_events")
        driver.init_data(str(tmp_path / "run"), tunesetup.datacards, "default", str(ROOT), pickle_dump=False)
        run = driver._run_cards(index=0, datacard=tunesetup.datacards[0], mc_steer=steering,
                                tunename="TUNE0", cdir=str(ROOT))[0]
        initial = driver.get_initial_param(tunesetup.param_space, {}, str(ROOT))
        prepare_bank(driver=driver, directory=directory, temporary=tmp_path / "basis", run=run,
                     tune=ROOT / "modeldata/TUNE0", initial=initial,
                     bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                     events=events, controls=controls,
                     cdir=str(ROOT), deadline=time.monotonic() + 120)
    else:
        initial = driver.initialize(run_name=str(tmp_path / "run"), tunesetup=tunesetup, mc_steer=steering,
                                    obs_module="default", cdir=str(ROOT), init_force=True,
                                    max_t=600, pickle_dump=False, processes=1, rngseed=9182)
    bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                         cuts=driver.cuts[0], controls=controls)
    report = bank.validate_source(bank.amplitude(initial), controls["closure_rtol"], directory / "closure.json")
    assert report["weighted_rms_relative_amplitude"] < controls["closure_rtol"]
    with pytest.raises(ValueError, match="RMS source closure"):
        bank.validate_source(bank.amplitude(initial) * 1.1, controls["closure_rtol"], tmp_path / "wrong_closure.json")
    if kaons:
        reuse = dict(driver=driver, run=run, tune=ROOT / "modeldata/TUNE0",
                     initial=initial, bounds=normalize_continuous_param_space(space),
                     pid=driver.pid[0][0], controls=controls, cdir=str(ROOT), nevents=16)
        plan = BankPlan(**reuse)
        assert not plan.mismatch(directory)
        fixed_phase = {name: limits for name, limits in reuse["bounds"].items() if name != key}
        assert not BankPlan(**{**reuse, "bounds": fixed_phase}).mismatch(directory)
        assert BankPlan(**{**reuse, "nevents": 32}).mismatch(directory)
        assert BankPlan(**{**reuse, "initial": {**initial, key: initial[key] + 0.1}}).mismatch(directory)
        if model == "GP" and not screened:
            # Reuse an earlier run through the real bank API and compare complex amplitudes
            root = directory.parents[4]
            copied = root / "reused/results/amplitude/0/nominal"
            assert BankPlan(**{**reuse, "nevents": 32}).reuse(copied, root) is None
            assert plan.reuse(copied, root) == directory
            reused = AmplitudeBank(driver=driver, directory=copied, obs=driver.obs[0], pid=driver.pid[0],
                                   cuts=driver.cuts[0], controls=controls)
            torch.testing.assert_close(reused.coherent(reused.coefficients(initial)), bank.coherent(bank.coefficients(initial)))
            intensity = directory / "source.bin.intensity"
            intensity.rename(intensity.with_suffix(".intensity._old"))
            assert plan.mismatch(copied)
    reconstructed = bank.coherent(bank.coefficients(initial)).detach().numpy()
    np.testing.assert_allclose(reconstructed, bank.amplitudes[0].numpy(), rtol=1e-9, atol=1e-12)
    if not kaons:
        # Only the matched source model defines physical MC histogram weights
        ordinary = readers.read_hepmc3(str(directory / "events.hepmc3"), obs=driver.obs[0], pid=driver.pid[0], cuts=driver.cuts[0])
        for prediction, sample, obs in zip(bank.predict(initial), ordinary, driver.obs[0], strict=True):
            expected = plot.histmc(sample, obs)
            for name in expected:
                np.testing.assert_allclose(prediction[name]["hdata"].counts_scaled, expected[name]["hdata"].counts_scaled, rtol=1e-9, atol=1e-12)
                np.testing.assert_allclose(prediction[name]["hdata"].errs_scaled, expected[name]["hdata"].errs_scaled, rtol=1e-9, atol=1e-12)

    varied = copy.deepcopy(initial)
    varied[key] += 0.07
    varied[f"RES|{resonance}:{model}:mass"] *= 0.998
    varied[f"RES|{resonance}:{model}:width"] *= 1.02
    direct_tune = tmp_path / "direct_tune"
    driver.create_steering_card(param_space=varied, tunename=str(direct_tune), cdir=str(ROOT))
    run = driver._run_cards(index=0, datacard=tunesetup.datacards[0], mc_steer=steering, tunename="TUNE0", cdir=str(ROOT))[0]
    direct_card = tmp_path / "direct.json"
    write_json_file(direct_card, input_card(driver, run, direct_tune))
    output = tmp_path / "nested output" / "direct.bin"
    command = [str(ROOT / "bin/ampfit"), "--input", str(direct_card), "--events", str(directory / "events.hepmc3"),
               "--output", str(output), "--batch-events", "3", "--closure-rtol", str(controls["closure_rtol"])]
    result = subprocess.run(command, cwd=ROOT, capture_output=True, text=True, timeout=120, check=False)
    assert result.returncode == 0, result.stdout + result.stderr
    expected = np.fromfile(output, dtype=np.complex128).reshape(bank.amplitudes[0].shape)
    np.testing.assert_allclose(bank.coherent(bank.coefficients(varied)).detach().numpy(), expected, rtol=1e-9, atol=1e-12)
    np.testing.assert_allclose(bank.amplitude(varied).detach().numpy(), expected, rtol=1e-9, atol=1e-12)
    # Check bounded validation batches against the native complex amplitudes
    bank.batch_events = 3
    with torch.no_grad():
        np.testing.assert_allclose(bank.amplitude(varied).numpy(), expected, rtol=1e-9, atol=1e-12)
    names = sorted(varied)
    value = torch.tensor([varied[name] for name in names], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(
        lambda theta: bank.intensity(bank.coefficients(dict(zip(names, theta, strict=True)))) / bank.reference, (value,))


# Compare all continuum form factors with native complex amplitudes inside and outside screening
@pytest.mark.parametrize("ff_type", ["exp", "power", "orear", "logexp"])
@pytest.mark.parametrize("screened", [False, True])
@pytest.mark.parametrize("wide", [False, True])
def test_continuum_form_factors(tmp_path, kaon_events, ff_type, screened, wide):
    from core.tune.drivers.graniitti.ampfit.ampcheck import validate_points

    driver = GraniittiDriver()
    prefix = "CON_GP|990:[321,321]:"
    form = driver._continuum_ff_template(ff_type)
    keys = [prefix + "FF_offshell." + field for field in form if field not in {"type", "norm"}]
    keys.append(prefix + "FF_transfer.LambdaInv2")
    if wide:
        keys = keys[-1:]
    fixed = {prefix + "FF_offshell.type": ff_type}
    initial = driver.get_initial_param(dict.fromkeys(keys), fixed, str(ROOT)) | fixed
    space = {key: {"type": "uniform", "lower": initial[key] * 0.9, "upper": initial[key] * 1.1} for key in keys}
    # Isolate transfer interpolation over a wide range without changing off shell factors
    if wide:
        space[keys[-1]]["lower"] = initial[keys[-1]] * 0.01
    tune = tmp_path / "source"
    driver.create_steering_card(param_space=initial, tunename=str(tune), cdir=str(ROOT))
    cards = [card(dataset_card(tmp_path, "STAR_1792394/KK"), nevents=16, loopscreen=screened,
                  xsmode="sample", swap_process="GP[RES+CON]<F>", integrator="VEGAS")]
    driver.init_data(str(tmp_path / "run"), cards, "default", str(ROOT), pickle_dump=False)
    run = driver._run_cards(index=0, datacard=cards[0], mc_steer={"tunesetup_name": tmp_path.name, "tune_default": "TUNE0"},
                           tunename="TUNE0", cdir=str(ROOT))[0]
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"] | {"nodes": 3}
    directory = tmp_path / "bank"
    prepare_bank(driver=driver, directory=directory, temporary=tmp_path / "basis", run=run, tune=tune,
                 initial=initial, bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                 events=kaon_events, controls=controls, cdir=str(ROOT), deadline=time.monotonic() + 600)
    bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                         cuts=driver.cuts[0], controls=controls)
    points = [{key: spec["lower"] for key, spec in space.items()}, {key: initial[key] * 1.04 for key in keys}]
    report = validate_points(driver=driver, bank=bank, run=run, points=points, directory=tmp_path / "check",
                             cdir=str(ROOT), controls=controls, events=16, timeout=600)
    assert report["results"][0]["rms_relative_amplitude"] < 1e-9
    assert report["results"][1]["rms_relative_amplitude"] < 0.01
    values = torch.tensor([points[1][key] for key in keys], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(
        lambda theta: torch.view_as_real(bank.coherent(bank.coefficients(initial | dict(zip(keys, theta, strict=True))))),
        (values,), fast_mode=True)


# Compare screened resonance form factors with native complex amplitudes
@pytest.mark.parametrize("dataset_name,vector", [("CMS_2752118/pipi_0p7", "rho_770"),
                                                ("STAR_1792394/KK", "phi_1020")])
def test_gp_resonance_form_factors(tmp_path, dataset_name, vector):
    from core.tune.drivers.graniitti.ampfit.ampcheck import validate_points

    driver = GraniittiDriver()
    dataset = dataset_card(tmp_path, dataset_name)
    gencard = load_json_file(tmp_path / "gencard.json")
    gencard["SCATTERING"]["RES"] = ["f0_980", vector]
    if vector == "phi_1020":
        # Retain the phi pole in the sampled mass range instead of the scalar threshold test range
        original = load_json_file(ROOT / "icepack/SOFTCEP" / dataset_name / "gencard.json", loader=pyjson5.load)
        gencard["GENCUTS"]["<F>"]["M"] = original["GENCUTS"]["<F>"]["M"]
    write_json_file(tmp_path / "gencard.json", gencard)
    keys = ["RES|f0_980:GP:FF_transfer.LambdaInv2", f"RES|{vector}:GP:FF_prod.Lambda2"]
    source = driver.get_initial_param(dict.fromkeys(keys), {}, str(ROOT))
    space = {key: {"type": "uniform", "lower": value * 0.75, "upper": value * 1.25}
             for key, value in source.items()}
    cards = [card(dataset, nevents=16, loopscreen=True, xsmode="sample",
                  swap_process="GP[RES+CON]<F>", integrator="VEGAS")]
    controls = load_settings(campaign_source("tune-gpom-star-cms-ampfit") + "#/optimizer/ampfit")["bank"]
    steering = dict(tune_default="TUNE0", tunesetup_name=tmp_path.name,
                    data_covariance_mode="diagonal", ampfit=controls)
    tunesetup = SimpleNamespace(param_space=space, aux_param_space={}, datacards=cards)
    driver.initialize(run_name=str(tmp_path / "run"), tunesetup=tunesetup, mc_steer=steering,
                      obs_module="default", cdir=str(ROOT), init_force=True,
                      max_t=600, pickle_dump=False, processes=2, rngseed=9182)
    directory = driver.amplitude_directory(cdir=str(ROOT), run_name=str(tmp_path / "run")) / "0/nominal"
    run = driver._run_cards(index=0, datacard=cards[0], mc_steer=steering,
                           tunename="TUNE0", cdir=str(ROOT))[0]
    reuse = dict(driver=driver, run=run, tune=directory / "tune",
                 initial=source, bounds=normalize_continuous_param_space(space), pid=driver.pid[0][0],
                 controls=controls, cdir=str(ROOT), nevents=16)
    assert BankPlan(**reuse).reusable_events(directory)
    assert not BankPlan(**{**reuse, "nevents": 32}).reusable_events(directory)
    bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                         cuts=driver.cuts[0], controls=controls)
    baseline = bank.coherent(bank.coefficients(source))
    changed = bank.coherent(bank.coefficients(source | {keys[1]: source[keys[1]] * 0.75}))
    assert torch.linalg.vector_norm(changed - baseline) > 1e-10 * torch.linalg.vector_norm(baseline)
    points = [{key: spec["lower"] for key, spec in space.items()},
              {key: value * 1.07 for key, value in source.items()},
              {keys[0]: space[keys[0]]["lower"], keys[1]: source[keys[1]] * 0.77},
              {keys[0]: space[keys[0]]["lower"], keys[1]: source[keys[1]] * 1.19}]
    report = validate_points(driver=driver, bank=bank, run=run, points=points,
                             directory=tmp_path / "check", cdir=str(ROOT), controls=controls,
                             events=8, timeout=120)
    assert report["results"][0]["rms_relative_amplitude"] < 1e-9
    assert report["results"][1]["rms_relative_amplitude"] < 0.01
    # Isolate exact production factors from transfer interpolation by using a grid node
    for result in report["results"][2:]:
        assert result["rms_relative_amplitude"] < 1e-9
    value = torch.tensor([points[1][key] for key in keys], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(
        lambda theta: bank.intensity(bank.coefficients(dict(zip(keys, theta, strict=True)))) / bank.reference,
        (value,))
    varied = bank.predict(points[1])
    for prediction, sample, obs in zip(varied, bank.sample, driver.obs[0], strict=True):
        reference = plot.histmc(sample, obs)
        for name, item in prediction.items():
            np.testing.assert_allclose(item["source_statistics"].errs_scaled, reference[name]["hdata"].errs_scaled,
                                       rtol=1e-12, atol=1e-12)


# Prepare one distributed shard in an independent simulator process
def prepare_shard(settings, shard, shared_root, index, shards):
    GraniittiDriver().initialize(**settings, bank_mode={"phase": "shard", "shard": shard, "jobs": shards,
                                                       "index": index, "shared_root": shared_root})


# Prepare screened banks concurrently and reuse their events for data correlations
@pytest.mark.parametrize("distributed", [False, True])
def test_ampfit_parallel_initialization(tmp_path, distributed):
    driver = GraniittiDriver()
    datacards = []
    for index in range(2):
        directory = tmp_path / str(index)
        directory.mkdir()
        datacards.append(card(dataset_card(directory, "CMS_2752118/pipi_0p7"), nevents=16,
                             loopscreen=True, xsmode="sample", swap_process="XP[RES+CON]<F>", integrator="VEGAS"))
    inactive = load_json_file(datacards[-1]["datacard"])
    inactive["active"] = False
    inactive_path = tmp_path / "inactive.json"
    write_json_file(inactive_path, inactive)
    datacards.append({**datacards[-1], "datacard": str(inactive_path)})
    tunesetup = SimpleNamespace(
        datacards=datacards, aux_param_space={},
        param_space={"RES|f0_980:XP:phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi}},
    )
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"]
    steering = dict(tune_default="TUNE0", tunesetup_name=tmp_path.name, data_covariance_mode="full", ampfit=controls)
    run_name = str(tmp_path / "run")
    settings = dict(run_name=run_name, tunesetup=tunesetup, mc_steer=steering,
                    obs_module="default", cdir=str(ROOT), init_force=True,
                    max_t=600, pickle_dump=False, processes=3, rngseed=9182)
    shared_root = str(tmp_path / "shared")
    if distributed:
        root = driver.amplitude_directory(cdir=shared_root, run_name=run_name)
        shared, generated = {}, {}
        with ProcessPoolExecutor(max_workers=2, mp_context=multiprocessing.get_context("spawn")) as pool:
            for index, shards in enumerate((2, 1)):
                driver.initialize(**settings, bank_mode={"phase": "sample", "index": index, "jobs": shards,
                                                        "shared_root": shared_root})
                directory = root / str(index) / "nominal"
                source = load_json_file(directory / "source.bin.json")
                assert source["indices"] == [0] and source["shards"] == shards
                assert not (directory / "amplitudes.bin").exists()
                shared.update({path: (path.stat().st_mtime_ns, path.read_bytes())
                               for path in directory.rglob("*") if path.is_file()})
                generated.update({path: path.stat().st_mtime_ns for folder in ("vgrid", "output")
                                  for path in (ROOT / folder).glob(f"*{tmp_path.name}*")})
                futures = [pool.submit(prepare_shard, settings, shard, shared_root, index, shards)
                           for shard in range(shards)]
                for future in futures:
                    future.result()
                if index == 0:
                    # The first sample completes its bank before the second sample exists
                    assert not (root / "1").exists()
                parts = [directory / "parts" / f"part_{shard}.bin" for shard in range(shards)]
                columns = [column for part in parts for column in load_json_file(str(part) + ".json")["indices"]]
                assert sorted(columns) == list(range(1, source["columns"] + 1))
                if shards > 1:
                    # A serial retry reproduces complex components without reevaluating the source
                    part = parts[1]
                    values = np.fromfile(part, dtype=np.complex128)
                    pool.submit(prepare_shard, {**settings, "processes": 1}, 1, shared_root, index, shards).result()
                    np.testing.assert_allclose(np.fromfile(part, dtype=np.complex128), values, rtol=1e-12, atol=1e-12)
                    driver.initialize(**settings, bank_mode={"phase": "finalize", "jobs": shards,
                                                            "index": index, "shared_root": shared_root})
                completed = directory / "finalized.pkl"
                assert completed.is_file() and (directory / "closure.json").is_file()
                assert not (directory / "amplitudes.bin").exists()
                assert load_json_file(directory / "amplitudes.bin.json")["parts"] == [
                    "source.bin", *(f"parts/part_{shard}.bin" for shard in range(shards))]
                assert not list(directory.glob("amplitudes.bin.*partial"))
                published = {path: path.stat().st_mtime_ns for path in directory.glob("amplitudes.bin*")}
                published[completed] = completed.stat().st_mtime_ns
                driver.initialize(**settings, bank_mode={"phase": "finalize", "jobs": shards,
                                                        "index": index, "shared_root": shared_root})
                assert published == {path: path.stat().st_mtime_ns for path in published}
                # Retry after metadata publication without rewriting shared complex amplitudes
                completed.rename(completed.with_name("finalized.pkl._old"))
                driver.initialize(**settings, bank_mode={"phase": "finalize", "jobs": shards,
                                                        "index": index, "shared_root": shared_root})
                assert completed.is_file()
                assert not (directory / "amplitudes.bin").exists()
        assert shared == {path: (path.stat().st_mtime_ns, path.read_bytes()) for path in shared}
        assert generated == {path: path.stat().st_mtime_ns for path in generated}
    if distributed:
        # INIT must require completed predictions instead of silently repeating preparation
        completed = root / "0/nominal/finalized.pkl"
        saved = completed.with_suffix(".pkl._old")
        completed.rename(saved)
        with pytest.raises(FileNotFoundError, match="finalized.pkl"):
            driver.initialize(**settings, bank_mode={"phase": "combine", "jobs": 2, "shared_root": shared_root})
        saved.rename(completed)
        # INIT loads completed predictions without repeating the bank reuse checks
        for directory in root.glob("*/nominal"):
            request_path = directory / "request.json"
            request_path.rename(request_path.with_suffix(".json._old"))
    initial = driver.initialize(**settings, bank_mode={"phase": "combine", "jobs": 2,
                                                       "shared_root": shared_root} if distributed else None)
    if distributed:
        local = driver.amplitude_directory(cdir=str(ROOT), run_name=run_name)
        assert not list(local.rglob("parts"))
        assert not list(local.rglob("finalized.pkl"))
    assert driver.data[-1] == []
    samples, predictions = [], []
    for index in range(2):
        directory = driver.amplitude_directory(cdir=str(ROOT), run_name=run_name) / str(index) / "nominal"
        bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[index], pid=driver.pid[index],
                             cuts=driver.cuts[index], controls=controls)
        if distributed:
            shared_bank = AmplitudeBank(driver=driver, directory=root / str(index) / "nominal",
                                       obs=driver.obs[index], pid=driver.pid[index], cuts=driver.cuts[index], controls=controls)
            torch.testing.assert_close(shared_bank.amplitude(initial), bank.amplitude(initial), rtol=0, atol=0)
        correlation = driver.data_covariance_payload["mc_correlation_samples"][index]
        weights = bank.sample[0]["sample_weights"]
        assert correlation["event_count_read"] == len(weights)
        assert correlation["weight_sum"] == pytest.approx(np.sum(weights))
        assert correlation["weight_squared_sum"] == pytest.approx(np.sum(np.square(weights)))
        np.testing.assert_allclose(bank.coherent(bank.coefficients(initial)).detach().numpy(),
                                   bank.amplitudes[0].numpy(), rtol=1e-9, atol=1e-12)
        predictions.append(bank.predict(initial))
        for prediction, sample, obs in zip(predictions[-1], bank.sample, driver.obs[index], strict=True):
            expected = plot.histmc(sample, obs)
            for name in expected:
                np.testing.assert_allclose(prediction[name]["hdata"].counts_scaled,
                                           expected[name]["hdata"].counts_scaled, rtol=1e-9, atol=1e-12)
        samples.append(bank.mass2.numpy())
    assert not np.array_equal(*samples)
    costs = driver.trial_costs(dict(mc=[*predictions, []], data=driver.data, obs=driver.obs, datasets=driver.datasets),
        dict(data_covariance_mode="full", mc_steer=steering, cost="chi2", cost_avg="sum", cost_rho="quadratic"))
    assert costs["cost_arr"]["chi2"][-1] == []


# Differentiate the real prepared bank twice with both data covariance and MC error modes
@pytest.mark.parametrize("covariance_mode", ["diagonal", "full"])
@pytest.mark.parametrize("mc_stat", ["source", "reweighted"])
@pytest.mark.parametrize("cost", ["ratio2", "gaussian"])
def test_ampfit_covariance(tmp_path, kaon_events, covariance_mode, mc_stat, cost):
    from core.numerics import lbfgsb
    from core.stats import cov
    from core.tune.optimizers.ampfit.report import amplitude_objective, trial_covariance

    driver = GraniittiDriver()
    driver.amplitude_path = str(tmp_path / "amplitude")
    run_name = f"covariance-{tmp_path.name}"
    cards = [card(dataset_card(tmp_path, "STAR_1792394/KK"), nevents=16, loopscreen=False,
                  xsmode="sample", swap_process="GP[RES+CON]<F>", integrator="VEGAS")]
    controls = {**load_settings(resource("tune/settings/ampfit.json"))["bank"], "mc_stat": mc_stat}
    steering = dict(tune_default="TUNE0", tunesetup_name=tmp_path.name,
                    data_covariance_mode=covariance_mode, ampfit=controls)
    driver.init_data(run_name, cards, "default", str(ROOT), pickle_dump=False)
    run = driver._run_cards(index=0, datacard=cards[0], mc_steer=steering, tunename="TUNE0", cdir=str(ROOT))[0]
    names = ["RES|f0_980:GP:FF_prod.Lambda2", "RES|f0_980:GP:g_ls(0,0)@NORM", "RES|f0_980:GP:phi@PHASE"]
    initial = driver.get_initial_param(dict.fromkeys(names, {}), {}, str(ROOT))
    space = {name: {"type": "float", "lower": value - abs(value) / 2 - 0.1,
                    "upper": value + abs(value) / 2 + 0.1} for name, value in initial.items()}
    prepare_bank(driver=driver, directory=Path(driver.amplitude_path) / "0/nominal", temporary=tmp_path / "basis",
                 run=run, tune=ROOT / "modeldata/TUNE0", initial=initial, bounds=space, pid=driver.pid[0][0],
                 events=kaon_events, controls=controls, cdir=str(ROOT), deadline=time.monotonic() + 120)
    if covariance_mode == "full":
        driver.data_covariance_payload = cov.build_data_covariance(data=driver.data, obs=driver.obs,
                                                                   mcdata_by_dataset=[None])
    param = dict(cdir=str(ROOT), run_name=run_name, datacards=cards, mc_steer=steering,
                 aux_param_space={}, parameter_space=[{"name": name, **space[name]} for name in names],
                 data_covariance_mode=covariance_mode, cost=cost, cost_avg="sum" if cost == "gaussian" else "dataset-mean", cost_rho="quadratic",
                 optimization={"optimizer": "ampfit", "ampfit": load_settings(resource("tune/settings/ampfit.json"))},
                 rngseed=0, max_t=300, chunksize=10000, processes=1)
    driver.prepare_trial_runtime(param)
    objective = amplitude_objective(driver, param, names)
    values = torch.tensor([initial[name] for name in names], dtype=torch.float64)
    value, gradient = lbfgsb.value_gradient(objective, values)
    hessian, _, _ = lbfgsb.hessian_covariance(objective, values)
    for index in range(len(values)):
        step = torch.zeros_like(values)
        step[index] = 1e-4 * values[index].abs().clamp_min(1.0)
        _, plus = lbfgsb.value_gradient(objective, values + step)
        _, minus = lbfgsb.value_gradient(objective, values - step)
        np.testing.assert_allclose(hessian[:, index], (plus - minus) / (2 * step[index]), rtol=1e-5, atol=1e-8)
    outputs = driver.evaluate_trial_outputs(config=initial, param=param, trial_id="covariance",
                                            tunename=str(tmp_path / "trial"))
    assert outputs["error"] is None
    assert float(value) == pytest.approx(outputs["likelihood"]["objective"]["value"])
    result = trial_covariance(driver=driver, record=outputs, param=param, output_root=tmp_path / "report")
    assert result["objective"]["name"] == ("gaussian" if cost == "gaussian" else "chi2")
    np.testing.assert_allclose(result["hessian"], hessian)
    np.testing.assert_allclose(result["gradient"], gradient)

    if covariance_mode == "diagonal" and mc_stat == "source":
        import hashlib
        import json
        import pickle
        from concurrent.futures import ThreadPoolExecutor

        import ray
        from core.tune.backends.ray import TrialOutputCallback, create_global_state
        from ray.cluster_utils import Cluster

        # Publish on the zero CPU head while workers compute the covariance
        cluster = Cluster()
        cluster.add_node(num_cpus=0, resources={"icetune_head": 1}, include_dashboard=False,
                         temp_dir=str(ROOT / "tmp/r"))
        cluster.add_node(num_cpus=2)
        ray.init(address=cluster.address)
        try:
            publication = {**param, "cdir": str(tmp_path)}
            experiment = tmp_path / "runs" / run_name
            actor = create_global_state(experiment_dir=str(experiment), cdir=str(tmp_path), run_name=run_name,
                                        cost=param["cost"], render_figures=True, require_head_resource=True,
                                        keep_pickles=False)
            callback = TrialOutputCallback(experiment_dir=str(experiment), param=publication, global_state=actor,
                                           simdriver=driver, runtime_param=param)
            raw = pickle.dumps(outputs)
            filename = "TUNE_icetune_covariance.pkl"
            descriptor = dict(filename=filename, relative_path=f"results/{filename}", size=len(raw),
                              sha256=hashlib.sha256(raw).hexdigest())
            staged = ray.get(actor.publish_trial_outputs.remote(
                trial_id=outputs["trial_id"], cost=outputs["metrics"][param["cost"]],
                transfer={"pickle": raw, "figures": None}, descriptor=descriptor,
                rendered_initial=False, rendered_best=False, render_pending=True))
            accepted = ray.get(actor.commit_trial_outputs.remote(outputs["trial_id"], staged["output_token"]))
            callback.render_records[accepted["render_token"]] = accepted["render_record"]
            callback.param["optimization"]["ampfit"]["covariance"] = False
            with ThreadPoolExecutor(max_workers=1) as calls:
                calls.submit(callback.flush, wait=True).result(timeout=360)
            root = callback.figure_dir
            summary = json.loads((root / "summary.json").read_text())
            assert summary["trial_id"] == outputs["trial_id"]
            assert not (root / "covariance").exists()
            assert callback.covariance_worker is None
            callback.param["optimization"]["ampfit"]["covariance"] = True
            callback.covariance_next_at = time.time() + 300
            callback.flush()
            assert callback.output_jobs["covariance"]["future"] is None
            # Finalization bypasses the interval when the accepted best has no uncertainty
            with ThreadPoolExecutor(max_workers=1) as calls:
                calls.submit(callback.flush, wait=True).result(timeout=360)
            covariance_root = root / "covariance"
            covariance = json.loads((covariance_root / "covariance.json").read_text())
            assert summary["trial_id"] == covariance["trial_id"] == outputs["trial_id"]
            assert "error" not in covariance
            np.testing.assert_allclose(covariance["hessian"], hessian)
            for basis in ("optimizer", "physical"):
                parameters = json.loads((covariance_root / basis / "parameters.json").read_text())
                assert parameters["trial_id"] == summary["trial_id"]
            assert callback._next_render() is None
            assert callback.covariance_next_at > time.time()
        finally:
            ray.shutdown()
            cluster.shutdown()


# Resolve independent Regge projections and their mixed products at complex amplitude level
@pytest.mark.parametrize("dataset_name,pdg", [("CMS_2752118/pipi_0p7", 211), ("STAR_1792394/KK", 321)])
@pytest.mark.parametrize("screened", [False, True])
@pytest.mark.parametrize("projections", [[0, 1], [0, 1, 2]])
def test_gp_continuum_rows(tmp_path, dataset_name, pdg, screened, projections):
    from core.tune.drivers.graniitti.ampfit.ampcheck import validate_points
    from core.tune.drivers.graniitti.tunesetup.continuum import _continuum_exchange_ids

    driver = GraniittiDriver()
    general = load_json_file(ROOT / "modeldata/TUNE0/GENERAL.json", loader=pyjson5.load)
    exchanges = _continuum_exchange_ids(general["PARAM_REGGE"]["PARAM_CON"]["GP"], (pdg, -pdg), tune_all_channels=True)
    bases = [f"CON_GP|{exchange}:[{pdg},{pdg}]/opposite:helicity(0,0,{-m})"
             for exchange in exchanges for m in projections]
    keys = [base + suffix for base in bases for suffix in ("@RE", "@IM")]
    initial = driver.get_initial_param(dict.fromkeys(keys), {}, str(ROOT))
    scale = max(abs(value) for value in initial.values())
    space = {key: {"type": "uniform", "lower": -2 * scale, "upper": 2 * scale} for key in keys}
    cards = [card(dataset_card(tmp_path, dataset_name), nevents=16, loopscreen=screened, xsmode="auto",
                  swap_process="GP[RES+CON]<F>", integrator="VEGAS")]
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"]
    steering = dict(tune_default="TUNE0", tunesetup_name=tmp_path.name, data_covariance_mode="diagonal", ampfit=controls)
    tunesetup = SimpleNamespace(param_space=space, aux_param_space={}, datacards=cards)
    driver.initialize(run_name=str(tmp_path / "run"), tunesetup=tunesetup, mc_steer=steering,
                      obs_module="default", cdir=str(ROOT), init_force=True,
                      max_t=600, pickle_dump=False, processes=1, rngseed=9182)
    run = driver._run_cards(index=0, datacard=cards[0], mc_steer=steering, tunename="TUNE0", cdir=str(ROOT))[0]
    directory = driver.amplitude_directory(cdir=str(ROOT), run_name=str(tmp_path / "run")) / "0/nominal"
    bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                         cuts=driver.cuts[0], controls=controls)
    points = [initial]
    for only_nonzero_m in (False, True):
        point = {}
        for index, base in enumerate(bases):
            m0 = base.endswith(",0)")
            point[base + "@RE"] = (0.0 if only_nonzero_m and m0 else
                                   scale * (0.6 if m0 else -0.3) * (1 + 0.1 * index))
            point[base + "@IM"] = (0.0 if only_nonzero_m and m0 else
                                   scale * (0.1 if m0 else 0.2) * (1 + 0.05 * index))
        points.append(point)
    report = validate_points(driver=driver, bank=bank, run=run, points=points, directory=tmp_path / "check",
                             cdir=str(ROOT), controls=controls, events=16, timeout=600)
    # Couplings enter polynomially, so fixed form factors permit numerical precision closure
    assert max(row["rms_relative_amplitude"] for row in report["results"]) < 1e-9
    parameters = points[1]
    torch.testing.assert_close(bank.amplitude(parameters), bank.coherent(bank.coefficients(parameters)))
    values = torch.tensor([parameters[key] for key in keys], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(
        lambda theta: torch.view_as_real(bank.amplitude(dict(zip(keys, theta, strict=True)))), (values,), fast_mode=True)
    assert torch.autograd.gradgradcheck(
        lambda theta: torch.view_as_real(bank.amplitude(dict(zip(keys, theta, strict=True)))), (values,), fast_mode=True)
