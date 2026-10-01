# Test GRANIITTI tuning-card helper utilities
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pyjson5
import pytest
from core.tune.drivers.graniitti.tunesetup import card, continuum, resonance

pytestmark = pytest.mark.usefixtures("tuning_context")

# Check a continuum-only dataset builds without resonance parameters
def test_ppbar_continuum_only():
    from core.tune.drivers.graniitti.tunesetup import central

    datasets = [card("icepack/SOFTCEP/STAR_1792394/ppbar/dataset.json", swap_process="GP[CON]<F>",
                     nevents=50000, loopscreen=True, xsmode="sample", integrator="VEGAS")]
    assert resonance.discover_resonances_from_datacards(datasets, model="GP") == []
    _, parameters, fixed = central(model="GP", datacards=datasets, continuum_options={
        "FF_offshell": True, "FF_transfer": True, "ff_type": "exp", "tune_couplings": True,
        "production_geometry": "projective", "pveto": True})
    assert parameters
    assert all(key.startswith("CON_GP|") and ":[2212,2212]" in key for key in parameters | fixed)


# Check elastic and central-production steering can coexist
def test_datacard_supports_mixed_processes():
    datacards = [
        card(
            "icepack/ELASTIC/STAR_1791591/dataset.json",
            nevents=200000,
            loopscreen=True,
            xsmode="reset",
            integrator="vegas",
        ),
        card(
            "icepack/SOFTCEP/STAR_1792394/pipi/dataset.json",
            nevents=100000,
            loopscreen=True,
            xsmode="sample",
            swap_process="GP[RES+CON]<F>",
            integrator="neurojac",
        ),
    ]

    assert "swap_process" not in datacards[0]
    assert datacards[0]["xsmode"] == "reset"
    assert datacards[0]["integrator"] == "VEGAS"
    assert datacards[1]["swap_process"] == "GP[RES+CON]<F>"
    assert datacards[1]["integrator"] == "NEUROJAC"
    assert datacards[1]["xsmode"] == "sample"

    final_states = continuum.discover_final_states_from_datacards(datacards)
    resonances = resonance.discover_resonances_from_datacards(
        datacards,
        model="GP",
    )

    assert (2212, 2212) not in final_states
    assert (211, -211) in final_states
    assert "f0_980" in resonances


# Check an explicit tunesetup density switch is retained in the run card
def test_datacard_preserves_force_density():
    datacard = card(
        "icepack/example/dataset.json",
        nevents=1000,
        loopscreen=False,
        xsmode="header",
        force_density=True,
    )

    assert datacard["force_density"] is True


# Check malformed density overrides fail at tunesetup construction
def test_datacard_rejects_non_boolean_force_density():
    with pytest.raises(TypeError, match="must be boolean"):
        card(
            "icepack/example/dataset.json",
            nevents=1000,
            loopscreen=False,
            xsmode="header",
            force_density=1,
        )


# Check process-specific steering is passed through unchanged
def test_datacard_process_specific_steering():
    datacard = card(
        "icepack/example/dataset.json",
        nevents=1000,
        loopscreen=False,
        xsmode="sample",
        custom_process_option={"mode": "example"},
    )

    assert datacard["custom_process_option"] == {"mode": "example"}


# Check icetune rejects integrators without a reusable adaptation grid
def test_datacard_rejects_unknown_integrator():
    with pytest.raises(ValueError, match="VEGAS or NEUROJAC"):
        card(
            "icepack/example/dataset.json",
            nevents=1000,
            loopscreen=False,
            xsmode="sample",
            integrator="FLAT",
        )




# Check a real continuum coupling has a finite phase variation in Cartesian coordinates
@pytest.mark.parametrize("phase", [0.0, 1.5707963267948966, 2.3])
@pytest.mark.parametrize("magnitude", [0.0, 1.7])
def test_con_cartesian_domains_normalize(phase, magnitude):
    from core.tune.parameters.space import normalize_param_space

    parameters = {}
    continuum._add_continuum_production_block_params(
        parameters, continuum_model="MP", exchange="990", pair_key="[211,-211]", sector="direct",
        block={"basis": "crossed_auto_min_L", "g": [magnitude, phase]},
        mode="mag_phase_cartesian", relative_range=1.0, equal_gp_m=False, geometry="direct",
        active_only=True, projections=None,
    )
    normalized = normalize_param_space(parameters)
    assert len(normalized) == (2 if magnitude > 0.0 else 0)
    for domain in normalized.values():
        assert domain["upper"] - domain["lower"] == pytest.approx(2.0 * magnitude)


# Reject invalid coupling intervals before constructing optimizer parameters
@pytest.mark.parametrize("magnitude,cartesian", [
    ((0.0, 2.0), ((float("nan"), float("nan")), (float("nan"), float("nan")))),
    ((0.0, 2.0), ((0.0, 1.0), (0.0, 0.0))),
    ((0.0, float("inf")), None),
    ((2.0, 1.0), None),
    ((-1.0, 1.0), None),
])
def test_coupling_domains_reject_invalid_intervals(magnitude, cartesian):
    from core.tune.drivers.graniitti.tunesetup.domains import add_magnitude_phase_params

    parameters = {}
    with pytest.raises(ValueError):
        add_magnitude_phase_params(
            parameters, base_key="g", magnitude_key="g[0]", phase_key="g[1]",
            mode="mag_phase_cartesian", magnitude_bounds=magnitude, cartesian_bounds=cartesian,
        )
    assert not parameters


# Build simultaneous studies from distinct source tunes without sharing model state
def test_selected_source_tunes_are_isolated(tmp_path):
    import copy
    import json
    from concurrent.futures import ThreadPoolExecutor
    from pathlib import Path

    from core.io.serialize import load_json_file
    from core.tune.drivers.graniitti.tunesetup.domains import load_model_json
    from core.tune.tunesetup import load_tunesetup, save_tunesetup

    from submit import CAMPAIGN_DIR, campaign_source
    settings = load_json_file(CAMPAIGN_DIR / "tunecards/graniitti/_defaults.json")

    root = Path(__file__).resolve().parents[3]
    original = load_model_json(root / "modeldata/TUNE0/GENERAL.json")
    paths, expected = [], []
    for index, scale in enumerate((1.0, 1.25)):
        source = copy.deepcopy(original)
        odderon = source["PARAM_SOFT"]["MODEL"]["single"]["EXCHANGE"]["O"]
        odderon["on"] = True
        odderon["g"][0][0] *= scale
        directory = tmp_path / str(index)
        directory.mkdir()
        (directory / "GENERAL.json").write_text(json.dumps(source))
        paths.append(directory)
        expected.append([x * odderon["g"][0][0] for x in settings["icetune"]["eikonal"]["EXCHANGE"]["g"]["diagonal_relative_bounds"]])

    # Exercise the actual loader and source model discovery in parallel
    def build(path):
        return load_tunesetup(cdir=root, simdriver="GRANIITTI", tune_default=str(path),
                            name=campaign_source("tune-elastic-single"))

    with ThreadPoolExecutor(2) as pool:
        studies = list(pool.map(build, paths))
    for index, study in enumerate(studies):
        domain = study.param_space["SOFT|MODEL.single:EXCHANGE.O.g[0,0]"]
        assert [domain.lower, domain.upper] == pytest.approx(expected[index])
        saved = tmp_path / f"{index}.json"
        save_tunesetup(tunesetup=study, path=saved, simdriver="GRANIITTI", tune_default=study.tune_default)
        assert load_tunesetup(cdir=root, simdriver="GRANIITTI", name=str(saved)).tune_default == str(paths[index])
        with pytest.raises(ValueError, match="different baseline"):
            load_tunesetup(cdir=root, simdriver="GRANIITTI", name=str(saved), tune_default=str(paths[1 - index]))


# Apply explicit JSON ranges and fixed coordinates through the actual study loader
def test_json_parameter_overrides(tmp_path):
    import json
    from pathlib import Path

    from core.io.serialize import load_json_file
    from core.tune.parameters.space import normalize_param_space
    from core.tune.tunesetup import load_tunesetup, save_tunesetup

    from submit import campaign_source

    root = Path(__file__).resolve().parents[3]
    source = campaign_source("tune-gpom-res")
    original = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)
    bounds = normalize_param_space(original.param_space)
    varied, fixed = sorted(bounds)[:2]
    definition = load_json_file(source, loader=pyjson5.load)
    interval = bounds[varied]
    width = interval["upper"] - interval["lower"]
    selected = {**interval, "lower": interval["lower"] + width / 4, "upper": interval["upper"] - width / 4}
    constant = (bounds[fixed]["lower"] + bounds[fixed]["upper"]) / 2
    definition.update(parameters={varied: selected}, fixed={fixed: constant})
    source = tmp_path / "source.json"
    source.write_text(json.dumps(definition))
    study = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)
    assert normalize_param_space(study.param_space)[varied] == selected
    assert fixed not in study.param_space
    assert study.aux_param_space[fixed] == pytest.approx(constant)
    saved = tmp_path / "resolved.json"
    save_tunesetup(tunesetup=study, path=saved, simdriver="GRANIITTI", tune_default=study.tune_default)
    restored = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=saved)
    assert normalize_param_space(restored.param_space) == normalize_param_space(study.param_space)
    assert restored.aux_param_space == study.aux_param_space
    definition["fixed"][varied] = constant
    source.write_text(json.dumps(definition))
    with pytest.raises(ValueError, match="both varied and fixed"):
        load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)


# Close the source JSON, saved setup and physical model card loop with varied and fixed parameters
@pytest.mark.parametrize("bound", ["lower", "upper"])
def test_json_model_cards_roundtrip(tmp_path, bound):
    import json
    from pathlib import Path

    from core.io.serialize import load_json_file
    from core.tune.drivers.graniitti.driver import GraniittiDriver
    from core.tune.tunesetup import load_tunesetup, save_tunesetup

    from submit import campaign_source

    root = Path(__file__).resolve().parents[3]
    targets = {
        "REGGE|s0": ("GENERAL.json", ["PARAM_REGGE", "s0"]),
        "CON_GP|990:[211,211]:FF_offshell.b": ("CON_GP.json", ["990", "[211,211]", "FF_offshell", "b"]),
        "RES|f0_980:GP:mass": ("RES/f0_980.json", ["PARAM_RES", "MODELS", "GP", "mass"]),
        "DECAY|9010221:[321,-321]:zeta.GP@PHASE": ("DECAYS.json", ["9010221", "[321,-321]", "zeta", "GP"]),
    }
    definition = load_json_file(campaign_source("tune-gpom-res"), loader=pyjson5.load)
    definition.update(parameters={}, fixed={})
    original = {}
    for index, (key, (filename, fields)) in enumerate(targets.items()):
        path = root / "modeldata/TUNE0" / filename
        original[path] = path.read_bytes()
        value = pyjson5.loads(original[path].decode())
        for field in fields:
            value = value[field]
        if index < 2:
            definition["parameters"][key] = {"lower": value * 0.9, "upper": value * 1.1, "type": "float"}
        else:
            definition["fixed"][key] = value * 1.05
    source = tmp_path / "source.json"
    source.write_text(json.dumps(definition, allow_nan=False))
    study = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)
    saved = tmp_path / "resolved.json"
    save_tunesetup(tunesetup=study, path=saved, simdriver="GRANIITTI", tune_default=study.tune_default)
    study = load_tunesetup(cdir=root, simdriver="GRANIITTI", name=saved)
    expected = {**definition["fixed"], **{key: spec[bound] for key, spec in definition["parameters"].items()}}
    driver = GraniittiDriver()
    config = driver.get_initial_param(study.param_space, study.aux_param_space, str(root))
    config.update({key: getattr(study.param_space[key], bound) for key in definition["parameters"]})
    output = tmp_path / "model"
    driver.create_steering_card(param_space={**config, **study.aux_param_space}, tunename=str(output), cdir=str(root))
    for key, (filename, fields) in targets.items():
        value = pyjson5.loads((output / filename).read_text())
        for field in fields:
            value = value[field]
        assert value == pytest.approx(expected[key]), key
    assert all(path.read_bytes() == content for path, content in original.items())


# Reject inconsistent shared model parameters before any fit or event generation
def test_json_shared_parameter_conflict(tmp_path):
    import copy
    import json
    from pathlib import Path

    from core.io.serialize import load_json_file
    from core.tune.tunesetup import load_tunesetup

    from submit import campaign_source

    definition = load_json_file(campaign_source("tune-gpom-res"), loader=pyjson5.load)
    duplicate = copy.deepcopy(definition["models"][0])
    duplicate["production"]["relative_bounds"]["default"][1] *= 2
    definition["models"].append(duplicate)
    source = tmp_path / "conflicting.json"
    source.write_text(json.dumps(definition))
    root = Path(__file__).resolve().parents[3]
    with pytest.raises(ValueError, match="Conflicting shared fit parameters"):
        load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)
    definition["models"] = []
    source.write_text(json.dumps(definition))
    with pytest.raises(ValueError, match="non-empty models"):
        load_tunesetup(cdir=root, simdriver="GRANIITTI", name=source)
