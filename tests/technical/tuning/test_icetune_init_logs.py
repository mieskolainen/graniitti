# Test GRANIITTI initialization logging and datacard preparation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import shutil
import types
from pathlib import Path

import numpy as np
import pytest
from core.stats.uncertainty import source_covariance
from core.tune.runtime import process as iceruntime

ROOT = Path(__file__).resolve().parents[3]

json5 = pytest.importorskip("pyjson5")
from core.tune.drivers.graniitti import driver as graniitti_driver


# Copy the normal runtime into an isolated generator work directory
def _prepare_cdir(tmp_path):
    cdir = tmp_path / "graniitti"
    tune0 = cdir / "modeldata" / "TUNE0"
    tune0.parent.mkdir(parents=True)
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune0)
    (cdir / "gencard").mkdir(parents=True, exist_ok=True)
    (cdir / "bin").mkdir()
    for path in ("bin/gr", "VERSION.json", "modeldata/mass_width_2026.mcd"):
        shutil.copy2(ROOT / path, cdir / path)
    return cdir


# Write a complete physical pion-pair input with bounded proposal adaptation
def _write_gencard(cdir: Path, *, output: str = "TESTOUT", energy=(6500, 6500)):
    payload = json5.loads((ROOT / "gencard/test.json").read_text())
    payload["SCATTERING"].update(ENERGY=list(energy), PROCESS="GP[CON]<F> -> pi+ pi-", RES=[])
    payload["GENERIC"].update(OUTPUT=output, INTEGRATOR="VEGAS", WEIGHTED=True, CORES=1)
    payload["INTEGRATOR"].update(min_samples=1000, max_samples=2000, precision=1.0)
    payload["INTEGRATOR"]["VEGAS"].update(ncall=1000, rounds=2)
    payload["GENCUTS"]["<F>"]["M"] = [0.5, 2.0]
    with open(cdir / "gencard" / "test.json", "w", encoding="utf-8") as f:
        json.dump(payload, f)


# Configure the production driver for a real inclusive pion mass histogram
def _make_driver(cdir: Path):
    driver = graniitti_driver.GraniittiDriver()
    driver.initialized = True
    driver.datasets = [
        {
            "active": True,
            "fit": {"normalization": "cross_section"},
            "sets": [{"sample": "nominal"}],
            "samples": [
                {
                    "name": "nominal",
                    "label": "GRANIITTI",
                    "gencard": "gencard/test.json",
                    "parameters": {},
                }
            ],
        }
    ]
    driver.dataset_paths = [str(cdir / "icepack" / "test" / "dataset.json")]
    observable = graniitti_driver.readers.get_observables("default")["M"]
    driver.obs = [[{"M": {**observable, "bins": np.linspace(0.5, 2.0, 4)}}]]
    driver.pid = [[[211, -211]]]
    driver.cuts = [["core.analysis.cuts.default"]]
    return driver


# Resolve one test run through the production filename logic
def _resolved_run(driver, cdir: Path, datacard: dict, *, tunename: str = "TUNE0"):
    return driver._run_card(
        index=0,
        datacard=datacard,
        mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
        tunename=tunename,
        cdir=str(cdir),
    )


# Generate a real VEGAS grid through the normal executable
def _write_valid_vgrid(run, *, integrator: str = "VEGAS", proposal_only: bool = False):
    path = Path(run.gridfile)
    cdir = path.parent.parent
    command = [str(cdir / "bin/gr"), "-i", run.inputfile, "-p", run.process,
               "-g", integrator, "-o", run.output, "-n", "-1" if proposal_only else "0",
               "-l", "false", "-c", "1", "-h", "0", "-m", str(cdir / "modeldata/TUNE0")]
    for setting in run.overrides:
        command.extend(("--set", setting))
    assert iceruntime.execute_cmd(cmd=command, cwd=str(cdir), max_t=120)
    assert path.is_file()
    return path


# Build one valid dataset with intentionally independent plot and fit normalization
def _policy_dataset(*, fit_normalization, plot_normalization, mc_scale=2.5):
    plot = {
        "normalization": plot_normalization,
        "ratio_uncertainty": "combined",
        "stack": False,
        "data_style": "hist",
    }
    if plot_normalization == "unit_density":
        plot["density_uncertainty"] = "scaled"
    return {
        "active": True,
        "type": "HEPDATA_SCALAR",
        "reader": "unused_reader",
        "datapath": "unused_data",
        "samples": [
            {
                "name": "nominal",
                "label": "GRANIITTI",
                "gencard": "gencard/test.json",
                "parameters": {},
            }
        ],
        "plot": plot,
        "fit": {"normalization": fit_normalization},
        "sets": [
            {
                "name": fit_normalization,
                "plotname": fit_normalization,
                "mc_scale": mc_scale,
                "pid": [211, -211],
                "cuts": "unused_cuts",
                "obs": "unused_obs",
                "hist": [],
            }
        ],
    }


# Check real HEPData density and covariance follow fit settings and explicit overrides
@pytest.mark.parametrize(('fit', 'override'), [('cross_section', None), ('unit_density', None),
                                             ('cross_section', True), ('unit_density', False)])
def test_init_data_normalization(tmp_path, data_card, fit, override):
    source = json5.loads(data_card.read_text())
    source['fit']['normalization'] = fit
    source['plot']['normalization'] = 'cross_section' if fit == 'unit_density' else 'unit_density'
    if source['plot']['normalization'] == 'unit_density':
        source['plot']['density_uncertainty'] = 'scaled'
    else:
        source['plot'].pop('density_uncertainty', None)
    data_card.write_text(json.dumps(source))
    steering = {'datacard': str(data_card)}
    if override is not None:
        steering['force_density'] = override
    driver = graniitti_driver.GraniittiDriver()
    driver.init_data(run_name='density', datacards=[steering], obs_module='default', cdir=str(tmp_path))
    density = fit == 'unit_density' or override is True
    assert json5.loads(data_card.read_text()) == source
    for subset in driver.data[0]:
        for item in subset.values():
            histogram = item['hdata']
            if density:
                assert histogram.integral() == pytest.approx(1.)
                widths = np.diff(histogram.bins)
                covariance = sum(source_covariance(error) for error in item['uncertainties'])
                assert widths @ covariance @ widths == pytest.approx(0., abs=1e-12)
                np.testing.assert_allclose(np.diag(covariance), histogram.errs_scaled**2)
            else:
                assert histogram.integral() > 0.
                assert not np.isclose(histogram.integral(), 1.)
                assert np.all(np.diag(histogram.covariance_scaled) > 0.)


# Exercise real generation, HepMC reading and serial or parallel MC normalization
@pytest.mark.parametrize('processes', [1, 2])
@pytest.mark.parametrize('density', [False, True])
def test_compute_mc_density(tmp_path, processes, density):
    cdir = _prepare_cdir(tmp_path)
    card = json5.loads((ROOT / 'gencard/test.json').read_text())
    card['SCATTERING']['PROCESS'] = 'GP[CON]<F> -> pi+ pi-'
    card['SCATTERING']['RES'] = []
    card['GENERIC']['INTEGRATOR'] = 'VEGAS'
    card['GENERIC']['WEIGHTED'] = True
    card['INTEGRATOR'].update(min_samples=1000, max_samples=2000, precision=1.)
    card['INTEGRATOR']['VEGAS'].update(ncall=1000, rounds=2)
    card['GENCUTS']['<F>']['M'] = [.5, 2.]
    (cdir / 'gencard/test.json').write_text(json.dumps(card))
    driver = _make_driver(cdir)
    driver.datasets[0] = _policy_dataset(
        fit_normalization='unit_density' if density else 'cross_section',
        plot_normalization='cross_section' if density else 'unit_density')
    observable = graniitti_driver.readers.get_observables('default')['M']
    driver.obs = [[{'M': {**observable, 'bins': np.linspace(.5, 2., 4)}}]]
    driver.pid = [[[211, -211]]]
    driver.cuts = [['core.analysis.cuts.default']]
    result = driver.compute(
        tunename='TUNE0', datacards=[{'nevents': 32, 'weighted': True, 'loopscreen': False, 'xsmode': 'sample'}],
        mc_steer={'tunesetup_name': 'DENSITY', 'tune_default': 'TUNE0'}, cdir=str(cdir),
        chunksize=16, processes=processes, rngseed=31, max_t=120)
    histogram = result[0][0]['M']['hdata']
    assert np.all(np.isfinite(histogram.counts_scaled))
    widths = np.diff(histogram.bins)
    variance = widths @ histogram.covariance_scaled @ widths
    if density:
        assert histogram.integral() == pytest.approx(1.)
        assert variance == pytest.approx(0., abs=1e-12)
    else:
        assert histogram.integral() > 0. and variance > 0.


# Check initialization preserves the requested normalization and screening modes
def test_prepare_init_datacards_xsmode_loopscreen():
    datacards = [
        {"nevents": 50000, "weighted": True, "loopscreen": True, "xsmode": "sample"},
        {"nevents": 1000, "weighted": False, "loopscreen": False, "xsmode": "reset"},
    ]

    init_datacards = graniitti_driver.prepare_init_datacards(datacards)

    assert init_datacards[0]["nevents"] == -1
    assert init_datacards[0]["loopscreen"] is True
    assert init_datacards[0]["xsmode"] == "sample"
    assert init_datacards[1]["nevents"] == -1
    assert init_datacards[1]["loopscreen"] is False
    assert init_datacards[1]["xsmode"] == "reset"


# Check MC correlation event and generator-weight controls are independent
def test_mc_correlated_sampling():
    datacards = [
        {"nevents": 50000, "weighted": True},
        {"nevents": 1000, "weighted": False},
    ]

    card_mode = graniitti_driver.prepare_mc_correlation_datacards(
        datacards,
        event_count=None,
        weighting="card",
    )
    unweighted = graniitti_driver.prepare_mc_correlation_datacards(
        datacards,
        event_count=25000,
        weighting="unweighted",
    )
    weighted = graniitti_driver.prepare_mc_correlation_datacards(
        datacards,
        event_count=30000,
        weighting="weighted",
    )

    assert [card["nevents"] for card in card_mode] == [50000, 1000]
    assert [card["weighted"] for card in card_mode] == [True, False]
    assert [card["nevents"] for card in unweighted] == [25000, 25000]
    assert [card["weighted"] for card in unweighted] == [False, False]
    assert [card["nevents"] for card in weighted] == [30000, 30000]
    assert [card["weighted"] for card in weighted] == [True, True]
    assert datacards[0]["nevents"] == 50000


# Check that only sample-mode initialization grids enter the immutable worker cache
@pytest.mark.parametrize(
    "tunesetup_name", ["TESTSET", "tmp/icetune/TESTSET.json", "/pool/condor/grdev/tmp/icetune/TESTSET.json"]
)
def test_reusable_vgrid_paths_excludes_reset_mode(tmp_path, tunesetup_name):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir, output="SHARED")
    driver = _make_driver(cdir)

    paths = driver.reusable_vgrid_paths(
        datacards=[
            {
                "loopscreen": True,
                "swap_process": "GP[RES+CON]<F>",
                "integrator": "NEUROJAC",
                "xsmode": "sample",
            }
        ],
        mc_steer={"tune_default": "TUNE0", "tunesetup_name": tunesetup_name},
        cdir=str(cdir),
    )
    reset_paths = driver.reusable_vgrid_paths(
        datacards=[{"loopscreen": True, "xsmode": "reset"}],
        mc_steer={"tune_default": "TUNE0", "tunesetup_name": tunesetup_name},
        cdir=str(cdir),
    )

    expected = _resolved_run(
        driver,
        cdir,
        {
            "loopscreen": True,
            "swap_process": "GP[RES+CON]<F>",
            "integrator": "NEUROJAC",
            "xsmode": "sample",
        },
    )
    assert paths == [expected.gridfile]
    assert Path(paths[0]).parent == cdir / "vgrid"
    assert "__ENERGY_6500_6500__CARD_" in paths[0]
    assert reset_paths == []


# Route several generator models to the same measurement set
def test_sample_indices_accept_model_overlays(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    driver = _make_driver(cdir)
    driver.datasets[0]["sets"] = [
        {"samples": ["impulse", "shadowing"]},
        {"sample": "other"},
    ]

    assert driver._sample_set_indices(index=0, sample_name="impulse") == [0]
    assert driver._sample_set_indices(index=0, sample_name="shadowing") == [0]
    assert driver._sample_set_indices(index=0, sample_name="other") == [1]


# Reject ambiguous fitted predictions before any generator command is constructed
def test_fit_rejects_overlapping_model_predictions(tmp_path):
    driver = _make_driver(_prepare_cdir(tmp_path))
    driver.datasets[0]["samples"] = [{"name": "impulse"}, {"name": "shadowing"}]
    driver.datasets[0]["sets"] = [{"samples": ["impulse", "shadowing"]}]
    with pytest.raises(ValueError, match="one MC sample per fitted set"):
        driver._run_cards(index=0, datacard={}, mc_steer={}, tunename="TUNE0", cdir=str(tmp_path))


# Check routed samples produce distinct grids and carry their generator overrides
def test_multi_sample_runs_route_gen_overrides(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir, output="ROUTED")
    driver = _make_driver(cdir)
    driver.datasets[0]["samples"] = [
        {"name": "low", "label": "Low energy", "gencard": "gencard/test.json",
         "parameters": {"SCATTERING.ENERGY": [15.3, 15.3]}},
        {"name": "high", "label": "High energy", "gencard": "gencard/test.json", "scale": 1.5,
         "parameters": {"SCATTERING.BEAM": ["p+", "p-"], "SCATTERING.ENERGY": [31.15, 31.15]}},
    ]
    driver.datasets[0]["sets"] = [{"sample": "high", "mc_scale": 2.0}, {"sample": "low"}]
    driver.obs[0] *= 2
    driver.pid[0] *= 2
    driver.cuts[0] *= 2
    datacard = {"nevents": 32, "weighted": True, "loopscreen": False, "xsmode": "sample"}
    steering = {"tunesetup_name": "TESTSET", "tune_default": "TUNE0"}
    runs = driver._run_cards(index=0, datacard=datacard, mc_steer=steering, tunename="TUNE0", cdir=str(cdir))
    result = driver.compute(tunename="TUNE0", datacards=[datacard], mc_steer=steering,
                            cdir=str(cdir), processes=1, rngseed=31, max_t=120)
    driver.datasets[0]["samples"][1]["scale"] = 1.0
    driver.datasets[0]["sets"][0]["mc_scale"] = 1.0
    unscaled = driver.compute(tunename="TUNE0", datacards=[datacard], mc_steer=steering,
                              cdir=str(cdir), processes=1, rngseed=31, max_t=120)
    for subset, reference, run, scale in zip(result[0], unscaled[0], runs[::-1], (3.0, 1.0), strict=True):
        expected = reference["M"]["hdata"].integral() * scale
        assert expected > 0.0
        assert subset["M"]["hdata"].integral() == pytest.approx(expected)
        grid = json.loads(Path(run.gridfile).read_text())
        assert grid["BEAM_ENERGY"] == pytest.approx(run.beam_energy)


# Check beam-energy changes always select distinct reusable grid files
def test_vgrid_beam_energy_identity(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    driver = _make_driver(cdir)
    datacard = {"loopscreen": True, "xsmode": "sample"}

    _write_gencard(cdir, energy=(6500, 6500))
    cms = _resolved_run(driver, cdir, datacard)
    _write_gencard(cdir, energy=(100, 100))
    star = _resolved_run(driver, cdir, datacard)

    assert cms.gridfile != star.gridfile
    assert "__ENERGY_6500_6500__" in cms.gridfile
    assert "__ENERGY_100_100__" in star.gridfile


# Check only proposal-only grids may cross the LOOPSCREEN compatibility boundary
def test_vgrid_proposal_screening(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    driver = _make_driver(cdir)
    run = _resolved_run(driver, cdir, {"loopscreen": True, "xsmode": "sample"})
    gridfile = _write_valid_vgrid(run, proposal_only=True)
    payload = json.loads(gridfile.read_text(encoding="utf-8"))
    payload["LOOPSCREEN"] = False

    assert graniitti_driver._vgrid_rebuild_reason(payload, run) is None

    payload["PROPOSAL_ONLY"] = False
    assert "LOOPSCREEN" in graniitti_driver._vgrid_rebuild_reason(payload, run)

    payload["PROPOSAL_ONLY"] = "true"
    assert (
        graniitti_driver._vgrid_rebuild_reason(payload, run) == "PROPOSAL_ONLY metadata is invalid"
    )


# Check stale schemas and wrong beam energies trigger grid rebuilding
@pytest.mark.parametrize("failure_mode", ("schema", "energy"))
def test_compute_rebuilds_incompatible_vgrid(tmp_path, failure_mode):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    driver = _make_driver(cdir)
    datacard = {"nevents": 0, "weighted": True, "loopscreen": False, "xsmode": "sample"}
    run = _resolved_run(driver, cdir, datacard)
    gridfile = _write_valid_vgrid(run)
    payload = json.loads(gridfile.read_text())
    if failure_mode == "schema":
        payload.pop("SCHEMA_VERSION")
    else:
        payload["BEAM_ENERGY"] = [100.0, 100.0]
    gridfile.write_text(json.dumps(payload))
    logs = cdir / "runs/init_logs"
    driver.compute(tunename="TUNE0", datacards=[datacard],
                   mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
                   cdir=str(cdir), init_log_dir=str(logs), processes=1, max_t=120)
    rebuilt = json.loads(gridfile.read_text())
    assert graniitti_driver._vgrid_rebuild_reason(rebuilt, run) is None
    assert rebuilt["VEGAS"]["xmat"]
    log = json.loads(next(logs.glob("*.json")).read_text())
    assert log["returncode"] == 0
    assert log["metadata"]["stage"] == "vgrid_reinit_incompatible"
    assert log["metadata"]["vgrid_rebuild_reason"]


# Check icetune swaps the complete process tag and preserves the final state
def test_process_swap_uses_full_pomeron_tag():
    process = "MP[CON]<F> -> pi+ pi- @RES{f0_980:1}"
    assert (
        graniitti_driver._swap_process_tag(
            process,
            "GP[RES+CON]<F>",
        )
        == "GP[RES+CON]<F> -> pi+ pi- @RES{f0_980:1}"
    )
    with pytest.raises(ValueError, match="full Pomeron process tag"):
        graniitti_driver._swap_process_tag(process, "[GP+GP]")


# Accept an existing output directory through the normal filesystem helper
def test_eos_mkdir_eexist(tmp_path):
    target = tmp_path / "results"
    target.mkdir()
    graniitti_driver.ensure_dir(str(target))
    assert target.is_dir()


# Check that directory creation is delegated to the shared EOS retry helper
def test_eos_tolerant_makedirs_creates_nested_dir(tmp_path):
    target = tmp_path / "results/plots"
    graniitti_driver.ensure_dir(str(target))
    assert target.is_dir()


# Initialize a real proposal and record successful completion in the run directory
def test_init_run_log_path(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    driver = _make_driver(cdir)
    datacard = {"nevents": 0, "weighted": True, "loopscreen": False, "xsmode": "sample", "integrator": "VEGAS"}
    run = _resolved_run(driver, cdir, datacard)
    logs = cdir / "runs/init_logs"
    driver.compute(tunename="TUNE0", datacards=[datacard],
                   mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
                   cdir=str(cdir), init_log_dir=str(logs), processes=1, max_t=120)
    grid = json.loads(Path(run.gridfile).read_text())
    assert grid["PROPOSAL_ONLY"]
    assert grid["INTEGRATOR"] == "VEGAS"
    log = json.loads(next(logs.glob("*.json")).read_text())
    assert log["returncode"] == 0
    assert log["metadata"]["stage"] == "vgrid_init"
    assert log["metadata"]["output"] == run.output
    assert not list((cdir / "output").glob("*.hepmc3"))


# Generate events from a saved proposal and preserve its initialization record
def test_compute_event_generation_init_log(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    driver = _make_driver(cdir)
    datacard = {"nevents": 0, "weighted": True, "loopscreen": False, "xsmode": "sample"}
    steering = {"tunesetup_name": "TESTSET", "tune_default": "TUNE0"}
    logs = cdir / "runs/init_logs"
    kwargs = dict(tunename="TUNE0", mc_steer=steering, cdir=str(cdir), init_log_dir=str(logs), processes=1, max_t=120)
    driver.compute(datacards=[datacard], **kwargs)
    before = {p.name: p.read_bytes() for p in logs.glob("*.json")}
    result = driver.compute(datacards=[{**datacard, "nevents": 32}], rngseed=31, **kwargs)
    assert result[0][0]["M"]["hdata"].integral() > 0.0
    assert before == {p.name: p.read_bytes() for p in logs.glob("*.json")}


# Preserve the simulator root cause and bounded output in command failures
def test_simulator_failure_diagnostics(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    card = cdir / "gencard/test.json"
    payload = json.loads(card.read_text())
    payload["SCATTERING"]["PROCESS"] = "INVALID[CON]<F> -> pi+ pi-"
    card.write_text(json.dumps(payload))
    driver = _make_driver(cdir)
    with pytest.raises(graniitti_driver.GraniittiCommandError) as caught:
        driver.compute(tunename="TUNE0",
                       datacards=[{"nevents": 32, "weighted": True, "loopscreen": False, "xsmode": "sample"}],
                       mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
                       cdir=str(cdir), processes=1, max_t=120)
    failure = caught.value.failure
    assert failure["returncode"] != 0
    assert failure["stage"] == "vgrid_initialization"
    assert failure["root_cause"] and "INVALID" in failure["output_tail"]
    assert failure["argv"][0] == str(cdir / "bin/gr")
    assert not list((cdir / "output").glob("*.hepmc3"))


# Check failed trial payloads keep the original GRANIITTI exception context
def test_trial_penalty_original_exception_context(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    card = cdir / "gencard/test.json"
    payload = json.loads(card.read_text())
    payload["SCATTERING"]["PROCESS"] = "INVALID[CON]<F> -> pi+ pi-"
    card.write_text(json.dumps(payload))
    driver = _make_driver(cdir)
    outputs = driver.evaluate_trial_outputs(
        config={}, param={"aux_param_space": {}, "cdir": str(cdir), "chunksize": 10000,
                          "datacards": [{"nevents": 32, "weighted": True, "loopscreen": False, "xsmode": "sample"}],
                          "max_t": 120, "mc_steer": {"tune_default": "TUNE0", "tunesetup_name": "TESTSET"},
                          "processes": 1}, trial_id="trial-failed", tunename="TUNE_FAILED")
    assert outputs["error_type"] == "GraniittiCommandError"
    assert outputs["failure"]["returncode"] != 0
    assert outputs["stage"] == "simulation"
    assert "INVALID" in outputs["traceback"]
    assert not outputs.get("results")


# Initialize the unscreened proposal used by subsequent screened sampling
def test_init_unscreened_proposal(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    source = cdir / "gencard/test.json"
    payload = json.loads(source.read_text())
    payload["GENERIC"]["WEIGHTED"] = False
    source.write_text(json.dumps(payload))
    driver = _make_driver(cdir)
    datacard = {"nevents": 0, "weighted": True, "loopscreen": True, "xsmode": "sample"}
    run = _resolved_run(driver, cdir, datacard)
    driver.compute(tunename="TUNE0", datacards=[datacard],
                   mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
                   cdir=str(cdir), processes=1, max_t=120)
    grid = json.loads(Path(run.gridfile).read_text())
    assert grid["PROPOSAL_ONLY"] and not grid["LOOPSCREEN"]
    assert graniitti_driver._vgrid_rebuild_reason(grid, run) is None


# Run the temporary tune through its explicit model directory and verify cleanup
def test_temporary_tune_model_path(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    _write_gencard(cdir)
    driver = _make_driver(cdir)
    tunename = "TUNE_icetune_node_a_trial_1_pid143"
    driver.create_steering_card(param_space={}, tunename=tunename, cdir=str(cdir), tune_default="TUNE0")
    model = graniitti_driver._resolve_tune_dir(cdir=str(cdir), tunename=tunename)
    assert (model / "GENERAL.json").is_file()
    assert model != cdir / "modeldata/TUNE0"
    datacard = {"nevents": 0, "weighted": True, "loopscreen": False, "xsmode": "sample"}
    run = _resolved_run(driver, cdir, datacard, tunename=tunename)
    driver.compute(tunename=tunename, datacards=[datacard],
                   mc_steer={"tunesetup_name": "TESTSET", "tune_default": "TUNE0"},
                   cdir=str(cdir), processes=1, max_t=120)
    assert Path(run.gridfile).is_file()
    assert not model.exists()


# Initialize a published dimuon measurement with real grids and optional MC correlations
@pytest.mark.parametrize("mode", ["diagonal", "full"])
def test_init_sampling_and_data_cov(tmp_path, data_card, mode):
    cdir = _prepare_cdir(tmp_path)
    dataset = json5.loads(data_card.read_text())
    dataset["samples"] = dataset["samples"][:1]
    card_path = data_card.parent / "gencard.json"
    card = json5.loads(card_path.read_text())
    card["GENERIC"].update(INTEGRATOR="VEGAS", CORES=1, WEIGHTED=True)
    card["INTEGRATOR"].update(min_samples=1000, max_samples=2000, precision=1.0)
    card["INTEGRATOR"]["VEGAS"].update(ncall=1000, rounds=2)
    card_path.write_text(json.dumps(card))
    data_card.write_text(json.dumps(dataset))
    tunesetup = types.SimpleNamespace(
        datacards=[{"datacard": str(data_card), "nevents": 128, "weighted": True,
                    "loopscreen": False, "xsmode": "sample"}],
        param_space={}, aux_param_space={})
    steering = {"tune_default": "TUNE0", "tunesetup_name": "BOOTSTRAP", "data_covariance_mode": mode,
                "mc_correlation_events": 128, "mc_correlation_weighting": "weighted"}
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.initialize(run_name="bootstrap", tunesetup=tunesetup, mc_steer=steering,
                                obs_module="default", cdir=str(cdir), init_force=False, max_t=120)
    assert initial == {}
    grids = driver.reusable_vgrid_paths(datacards=tunesetup.datacards, mc_steer=steering, cdir=str(cdir))
    assert len(grids) == 1
    grid = json.loads(Path(grids[0]).read_text())
    assert grid["PROPOSAL_ONLY"] and not grid["LOOPSCREEN"]
    events = list((cdir / "output").glob("*.hepmc3"))
    if mode == "diagonal":
        assert driver.data_covariance_payload is None
        assert not events
    else:
        samples = driver.data_covariance_payload["mc_correlation_samples"]
        assert len(samples) == 1
        assert samples[0]["event_count_read"] == 128
        assert samples[0]["generator_weighted"]
        covariance = driver.data_covariance_payload["data_total_covariance"]
        assert np.isfinite(covariance).all()
        np.testing.assert_allclose(covariance, covariance.T, atol=1e-14)
        assert np.linalg.eigvalsh(covariance).min() >= -1e-12 * np.max(np.diag(covariance))
