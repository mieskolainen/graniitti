# Test icetune production and GRANIITTI tuning helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import re
import shutil
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core.io.serialize import load_json_file
from core.tune import core as icetune
from core.tune.tunesetup import load_tunesetup, save_tunesetup

from submit import CAMPAIGN_DIR, campaign_source

ROOT = Path(__file__).resolve().parents[3]

from core.tune.drivers.graniitti import driver as graniitti_driver
from core.tune.drivers.graniitti.tunesetup import continuum, domains, eikonal, resonance
from core.tune.parameters import space as parameter_space
from core.tune.parameters import tools

XP_P_POMERON_COUPLING = "CON_XP|995:[2212,2212]/opposite:g_ls(1,1)"
XP_PI_POMERON_COUPLING = "CON_XP|995:[211,211]/opposite:g_ls(2,0)"
MP_PI_POMERON_COUPLING = "CON_MP|995:[211,211]/opposite:g"
MP_K_POMERON_COUPLING = "CON_MP|995:[321,321]/opposite:g"
XP_P_PRODUCTION = "CON_XP|995:[2212,2212]/opposite"


pytestmark = pytest.mark.usefixtures("tuning_context")

# Compute the active eta-style inline block from one resonance model
def _active_model_block(param_res, model_name):
    model = param_res["PARAM_RES"]["MODELS"][model_name]
    pairs = [key for key in model if re.fullmatch(r"\[-?\d+,-?\d+\]", key)]
    assert len(pairs) == 1
    return model[pairs[0]]


# Load one checked-in TUNE0 JSON5 card
def _load_tune0_card(relative_path):
    with open(ROOT / "modeldata" / "TUNE0" / relative_path, encoding="utf-8") as handle:
        return load_json_file(handle.name, loader=graniitti_driver.json5.load)


# Compute one TUNE0 scalar resonance channel in Cartesian coordinates
def _tune0_resonance_channel(resonance, model):
    card = _load_tune0_card(f"RES/{resonance}.json")
    block = _active_model_block(card, model)
    row = block["g"]
    return tools.encode_polar(float(row[0]), float(row[1]))


# Compute one TUNE0 continuum coupling row in Cartesian coordinates
def _tune0_continuum_channel(model, row_index, exchange=None):
    exchange = exchange or {"XP": "995", "MP": "995"}.get(model, "990")
    pair_keys = ["[2212,2212]", "[211,211]", "[321,321]"]
    block = _load_tune0_card(f"CON_{model}.json")[exchange][pair_keys[row_index]]["opposite"]
    if "g" in block:
        row = block["g"]
    else:
        field = "g_ls" if block["basis"] == "crossed_ls" else "helicity"
        row = block[field][0][-2:]
    return tools.encode_polar(float(row[0]), float(row[1]))


# Compute the magnitude or phase key for one coupling representation
def _channel_component_key(base_key, column):
    if base_key.endswith(":g"):
        return f"{base_key}[{column}]"
    return tools.magnitude_key(base_key) if column == 2 else tools.raw_phase_key(base_key)


# Check that a coherent resonance channel uses Cartesian complex coordinates
def _assert_cartesian_channel(param_space, base_key):
    assert tools.complex_re_key(base_key) in param_space
    assert tools.complex_im_key(base_key) in param_space
    assert f"{base_key}[0]" not in param_space
    assert f"{base_key}[1]" not in param_space


# Check that a fixed-strength resonance channel exposes one periodic phase
def _assert_raw_phase_channel_only(param_space, base_key):
    magnitude_column, phase_column = (0, 1) if base_key.endswith(":g") else (2, 3)
    mag_key = _channel_component_key(base_key, magnitude_column)
    phase_key = _channel_component_key(base_key, phase_column)

    assert mag_key not in param_space
    assert phase_key not in param_space
    assert tools.raw_phase_key(phase_key) in param_space
    assert _uniform_bounds(param_space[tools.raw_phase_key(phase_key)]) == pytest.approx(
        (-math.pi, math.pi)
    )
    assert tools.complex_re_key(base_key) not in param_space
    assert tools.complex_im_key(base_key) not in param_space


# Check that one resonance model has exactly one periodic resonance phase
def _assert_resonance_phi(param_space, resonance_name, model):
    key = tools.raw_phase_key(f"RES|{resonance_name}:{model}:phi")
    assert key in param_space
    assert _uniform_bounds(param_space[key]) == pytest.approx((-math.pi, math.pi))
    assert tools.phase_u_key(f"RES|{resonance_name}:{model}:phi") not in param_space
    assert tools.phase_v_key(f"RES|{resonance_name}:{model}:phi") not in param_space


# Compute numeric bounds from a Ray Tune uniform parameter
def _uniform_bounds(spec):
    return spec.lower, spec.upper


# Copy the checked-in tuning defaults into an isolated work directory
@pytest.fixture
def cdir(tmp_path):
    cdir = tmp_path / "graniitti"
    tune0 = cdir / "modeldata" / "TUNE0"
    tune0.parent.mkdir(parents=True)
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune0)
    return cdir


# Read modified source cards through the normal JSON5 loader in an isolated tune
@pytest.fixture
def source_model(cdir):
    directory = cdir / "modeldata" / "TUNE0"

    # Write an input card and invalidate its previous decoded contents
    def write(name, payload):
        (directory / name).write_text(json.dumps(payload), encoding="utf-8")
        domains.load_model_json.cache_clear()

    with domains.context(cdir=ROOT, model_path=directory, settings=domains.settings()):
        yield write
    domains.load_model_json.cache_clear()


# Define signed scalar and tensor LS inputs independently of the source tune
@pytest.fixture
def ls_model(source_model):
    for name, model, labels in (("f0_980", "GP", [(0, 0), (2, 2), (4, 4)]),
                                ("f2_1270", "XP", [(0, 2), (2, 0), (2, 2), (2, 4), (4, 2), (4, 4), (6, 4)])):
        card = _load_tune0_card(f"RES/{name}.json")
        block = _active_model_block(card, model)
        block["basis"] = "g_ls"
        block["g_ls"] = [[orbital, spin, 1.0 / (index + 1), -math.pi if index == 2 else 0.0]
                         for index, (orbital, spin) in enumerate(labels)]
        source_model(f"RES/{name}.json", card)


# Create one generated tuning card from the checked-in defaults
def _create_card(cdir, tunename, params, driver=None):
    driver = driver or graniitti_driver.GraniittiDriver()
    return driver.create_steering_card(
        param_space=params,
        tunename=tunename,
        cdir=str(cdir),
        tune_default="TUNE0",
    )


# Load one generated tuning card with its JSON5 comments and shared references
def _load_generated_card(cdir, tunename, relative_path):
    return load_json_file(cdir / "modeldata" / tunename / relative_path, loader=graniitti_driver.json5.load)


# Check resonance FF_prod scales reach their separate model fields
def test_gp_res_model_ff_prod_round_trip(cdir):
    names = ["f0_980", "f2_1270"]
    params = resonance.res_form_factors(fields={"FF_prod"}, model="GP", resonances=names)
    values = {f"RES|{name}:GP:FF_prod.Lambda2": 0.8 + index for index, name in enumerate(names)}
    assert set(params) == set(values)
    _create_card(cdir, "FF_PROD", values)
    general = _load_generated_card(cdir, "FF_PROD", "GENERAL.json")
    assert "PARAM_RES" not in general["PARAM_REGGE"]
    for name in names:
        card = _load_generated_card(cdir, "FF_PROD", f"RES/{name}.json")
        assert card["PARAM_RES"]["MODELS"]["GP"]["FF_prod"]["Lambda2"] == pytest.approx(
            values[f"RES|{name}:GP:FF_prod.Lambda2"]
        )


# Check each joint tune varies RES and CON against STAR and writes real cards
@pytest.mark.parametrize(
    "model,tunesetup,campaign",
    [
        ("GP", "tune-gpom-res-con", "tune-gpom-res-con"),
        ("TP", "tune-tpom-res-con", "tune-tpom-res-con"),
        ("XP", "tune-xpom-res-con", "tune-xpom-res-con"),
        ("MP", "tune-mpom-res-con-coh", "tune-mpom-res-con-coh"),
        ("MP", "tune-mpom-res-con-incoh", "tune-mpom-res-con-incoh"),
    ],
)
@pytest.mark.parametrize("use_campaign", [False, True])
def test_res_con_tunesetup_round_trip(monkeypatch, cdir, model, tunesetup, campaign, use_campaign):

    from submit import campaign as campaign_config

    if use_campaign:
        catalog = campaign_config.load_campaign_catalog(CAMPAIGN_DIR / "campaigns.yml")
        environment = campaign_config.resolve(
            catalog, campaign_name=campaign,
        )
        for key, value in environment.items():
            monkeypatch.setenv(key, value)
    module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(tunesetup))

    assert all(entry["swap_process"] == f"{model}[RES+CON]<F>" for entry in module.datacards)
    assert all(entry["integrator"] == "VEGAS" for entry in module.datacards)
    assert any(key.startswith("RES|") and ":g" in key for key in module.param_space)
    assert any(
        key.startswith(f"CON_{model}|") and (":g" in key or ":helicity(" in key)
        for key in module.param_space
    )
    discovered = resonance.discover_resonances_from_datacards(module.datacards, model=model)
    targets = resonance.default_kk_branching_ratio_targets(discovered)
    assert {key for key in module.param_space if key.endswith(":BR")} == {
        f"DECAY|{pdg}:[{','.join(map(str, pair))}]:BR" for pdg, pair in targets
    }

    saved_path = cdir / "tunesetup.json"
    save_tunesetup(tunesetup=module, path=saved_path, simdriver="GRANIITTI", tune_default="TUNE0")
    saved = load_tunesetup(cdir=cdir, simdriver="GRANIITTI", name="tunesetup.json")
    assert saved.datacards == module.datacards
    assert saved.aux_param_space == module.aux_param_space
    assert parameter_space.normalize_param_space(saved.param_space) == parameter_space.normalize_param_space(module.param_space)
    save_tunesetup(tunesetup=saved, path=saved_path, simdriver="GRANIITTI", tune_default="TUNE0")
    module.datacards[0]["nevents"] += 1
    with pytest.raises(ValueError, match="different tuning definition"):
        save_tunesetup(tunesetup=module, path=saved_path, simdriver="GRANIITTI", tune_default="TUNE0")
    module = saved

    values = {key: (domain.lower + domain.upper) / 2 for key, domain in module.param_space.items()}
    driver = graniitti_driver.GraniittiDriver()
    _create_card(cdir, tunesetup, values | module.aux_param_space, driver)
    restored = driver.get_initial_param(module.param_space, {}, str(cdir), tune_default=tunesetup)
    assert restored == pytest.approx(values)


# Check every supported TP resonance class uses its implemented tensor basis
@pytest.mark.parametrize(
    "resonance_name,channel,labels",
    [
        ("f0_500", "[995,995]", "0,0 2,2"),
        ("eta", "[995,995]", "1,1 3,3"),
        ("rho_770", "[22,995]", "Gamma0 Gamma2"),
        ("f1_1420", "[995,995]", "2,2 4,4"),
        ("f2_1270", "[995,995]", "0,2 (2,0)-(2,2) (2,0)+(2,2) 2,4 4,2 4,4 6,4"),
    ],
)
def test_tp_res_tensor_names_physical(resonance_name, channel, labels):
    assert resonance.tensor_names(resonance_name, channel)[0] == tuple(
        f"g_tensor({label})" for label in labels.split()
    )


# Check scalar, mixed LS and Gamma TP names update exact resonance rows
def test_tp_res_tensor_names_roundtrip(cdir):
    params = {
        "RES|f0_500:TP:[995,995]:g_tensor(0,0)": 0.31,
        "RES|f2_1270:TP:[995,995]:g_tensor((2,0)-(2,2))": -0.42,
        "RES|rho_770:TP:[22,995]:g_tensor(Gamma2)": 5.1,
    }
    defaults = _create_card(cdir, "TUNE_TEST_TP_RESONANCE_NAMES", params)
    path = cdir / "modeldata" / "TUNE_TEST_TP_RESONANCE_NAMES" / "RES"
    assert set(defaults) == set(params)
    for key, value in params.items():
        _, resonance_name, _, channel, name = tools.substring_split(key, ["|", ":"])
        row = resonance.tensor_names(resonance_name, channel)[0].index(name)
        card = json.loads((path / f"{resonance_name}.json").read_text())
        assert card["PARAM_RES"]["MODELS"]["TP"][channel]["g_tensor"][row] == pytest.approx(value)


# Check the shared TP production switch can keep tensor couplings fixed
def test_tp_resonance_setup_fixed_production():
    params, auxiliary = resonance.setup(
        model="TP", resonances=["f0_980"], production={"mode": "none"}
    )
    assert params == auxiliary == {}
    with pytest.raises(ValueError, match="real tensor couplings"):
        resonance.setup(model="TP", resonances=["f0_980"], production={"mode": "mag_phase_cartesian"})


# Check normalized TP tensor coordinates preserve the real coefficient vector
@pytest.mark.parametrize("geometry", ["projective", "spherical"])
def test_tp_res_tensor_direction_roundtrip(geometry, cdir):
    driver = graniitti_driver.GraniittiDriver()
    param_space, _ = resonance.tensor_setup(
        resonances=["f2_1270"],
        geometry=geometry,
    )
    initial = driver.get_initial_param(param_space, {}, str(cdir))
    tunename = f"TUNE_TEST_TP_{geometry.upper()}"
    original = _create_card(cdir, tunename, initial, driver)
    generated = _load_generated_card(cdir, tunename, "RES/f2_1270.json")
    source = _load_tune0_card("RES/f2_1270.json")
    generated_tensor = generated["PARAM_RES"]["MODELS"]["TP"]["[995,995]"]["g_tensor"]
    source_tensor = source["PARAM_RES"]["MODELS"]["TP"]["[995,995]"]["g_tensor"]

    assert original == pytest.approx(initial)
    assert generated_tensor == pytest.approx(source_tensor)
    suffix_check = (
        tools.is_projective_angle_key if geometry == "projective" else tools.is_spherical_angle_key
    )
    assert any(suffix_check(key) for key in param_space)


# Check the shared resonance phase surface has one circular coordinate per model
def test_res_phi_surface_periodic_no_row_phases():
    for model in ("MP", "XP", "GP", "TP"):
        params, aux = resonance.phase_setup(
            resonances=["f0_500"],
            model=model,
            mode="phase_raw",
        )
        assert aux == {}
        assert set(params) == {tools.raw_phase_key(f"RES|f0_500:{model}:phi")}
        assert _uniform_bounds(next(iter(params.values()))) == pytest.approx((-math.pi, math.pi))
        assert not any("g_ls" in key or "helicity" in key for key in params)

    with pytest.raises(ValueError, match="none.*phase_raw"):
        resonance.phase_setup(
            resonances=["f0_500"],
            model="GP",
            mode="phase_cayley",
        )


# Check that eikonal steering follows the source-card matrix Odderon setup
def test_eikonal_matrix_odderon_params(cdir):
    tunename = "TUNE_TEST_ZERO_ODDERON"

    param, aux = eikonal.setup(model="single", pFF="DPOW")
    assert eikonal.ordered_dpow_key("SOFT|MODEL.single:FF.P.param[0,0]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]) in param
    assert (
        eikonal.descending_coupling_key("SOFT|MODEL.single:EXCHANGE.P.g[0,0]") in param
    )
    assert not any(":EXCHANGE.P.g[1," in key for key in param)
    assert "SOFT|MODEL.single:EXCHANGE.O.alpha[0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O.sign" not in param
    assert "SOFT|MODEL.single:EXCHANGE.O.ff" in aux
    assert "SOFT|MODEL.single:EXCHANGE.O3g.g[0,0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O3g.alpha[0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O3g.alpha[1]" in param
    assert "SOFT|MODEL.single:FF.O3g.param[0,0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O3g.ff" in aux

    _create_card(cdir, tunename, aux)
    general = _load_generated_card(cdir, tunename, "GENERAL.json")

    tune0_general = _load_tune0_card("GENERAL.json")
    assert general["PARAM_SOFT"]["active_model"] == "single"
    generated = general["PARAM_SOFT"]["MODEL"]["single"]
    source = tune0_general["PARAM_SOFT"]["MODEL"]["single"]
    assert generated["EXCHANGE"]["P"]["ff"] == "P"
    assert generated["EXCHANGE"]["O"]["g"] == source["EXCHANGE"]["O"]["g"]
    assert generated["EXCHANGE"]["O3g"]["g"] == source["EXCHANGE"]["O3g"]["g"]
    assert generated["EXCHANGE"]["O"]["sign"] == -1


# Check enabled source Odderon exchange exposes Odderon parameters
def test_eikonal_tunes_enabled_odderon_params(source_model):
    general = json.loads(json.dumps(eikonal._load_general_data()))
    general["PARAM_SOFT"]["MODEL"]["single"]["EXCHANGE"]["O"]["g"] = [[0.5]]
    source_model("GENERAL.json", general)

    param, aux = eikonal.setup(model="single", pFF="DPOW")
    assert "SOFT|MODEL.single:EXCHANGE.O.g[0,0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O.alpha[0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O.alpha[1]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O.sign" not in param
    assert eikonal.ordered_dpow_key("SOFT|MODEL.single:FF.O.param[0,0]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]) in param
    assert eikonal.ordered_dpow_key("SOFT|MODEL.single:FF.O.param[0,1]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]) in param
    assert "SOFT|MODEL.single:FF.O.param[0,2]" in param
    assert aux["SOFT|MODEL.single:EXCHANGE.O.ff"] == "O"


# Check enabled hard three-gluon exchange exposes its elastic tune block
def test_eikonal_tunes_enabled_o3g_params(source_model):
    general = json.loads(json.dumps(eikonal._load_general_data()))
    general["PARAM_SOFT"]["MODEL"]["single"]["EXCHANGE"]["O3g"]["g"] = [[0.28]]
    source_model("GENERAL.json", general)

    param, aux = eikonal.setup(model="single", pFF="DPOW")
    assert "SOFT|MODEL.single:EXCHANGE.O3g.g[0,0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O3g.alpha[0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O3g.alpha[1]" in param
    assert "SOFT|MODEL.single:FF.O3g.param[0,0]" in param
    assert aux["SOFT|MODEL.single:EXCHANGE.O3g.ff"] == "O3g"


# Check secondary Reggeon tuning includes both trajectories and their shared form factor
def test_eikonal_tunes_enabled_reggeon_params():
    param, aux = eikonal.setup(model="double", pFF="DPOW", tune_reggeons=True)
    source = eikonal.source_model_block("double")["EXCHANGE"]

    for exchange in ("R_f2", "R_rho"):
        assert _uniform_bounds(
            param[f"SOFT|MODEL.double:EXCHANGE.{exchange}.alpha[0]"]
        ) == pytest.approx((0.40, 0.70))
        assert _uniform_bounds(
            param[f"SOFT|MODEL.double:EXCHANGE.{exchange}.alpha[1]"]
        ) == pytest.approx((0.70, 1.10))
        for row in range(2):
            source_value = source[exchange]["g"][row][row]
            assert _uniform_bounds(
                param[f"SOFT|MODEL.double:EXCHANGE.{exchange}.g[{row},{row}]"]
            ) == pytest.approx((0.5 * source_value, 1.5 * source_value))
        transition = eikonal.symmetric_coupling_key(
            f"SOFT|MODEL.double:EXCHANGE.{exchange}.g[0,1]"
        )
        assert _uniform_bounds(param[transition]) == pytest.approx((-1.0, 1.0))
        assert aux[f"SOFT|MODEL.double:EXCHANGE.{exchange}.ff"] == "R"

    for row in range(2):
        assert _uniform_bounds(param[f"SOFT|MODEL.double:FF.R.param[{row},0]"]) == pytest.approx(
            (0.01, 10.0)
        )
        assert not any(f":FF.R.param[{row},{column}]" in key for column in (1, 2) for key in param)


# Check elastic tuning fixes every Pomeron transition coupling to zero
def test_eikonal_pomeron_coupling_matrix_diagonal(cdir):
    param, aux = eikonal.setup(model="double", pFF="DPOW")
    raw_key = "SOFT|MODEL.double:EXCHANGE.P.g[0,1]"
    fixed_key = eikonal.symmetric_coupling_key(raw_key)

    assert fixed_key not in param
    assert aux[fixed_key] == pytest.approx(0.0)

    source_path = cdir / "modeldata" / "TUNE0" / "GENERAL.json"
    source = graniitti_driver.json5.loads(source_path.read_text(encoding="utf-8"))
    coupling = source["PARAM_SOFT"]["MODEL"]["double"]["EXCHANGE"]["P"]["g"]
    coupling[0][1] = 2.5
    coupling[1][0] = 2.5
    source_path.write_text(json.dumps(source), encoding="utf-8")

    _create_card(cdir, "TUNE_TEST_DIAGONAL_POMERON", aux)
    generated_path = cdir / "modeldata" / "TUNE_TEST_DIAGONAL_POMERON" / "GENERAL.json"
    generated = json.loads(generated_path.read_text(encoding="utf-8"))
    matrix = generated["PARAM_SOFT"]["MODEL"]["double"]["EXCHANGE"]["P"]["g"]
    assert matrix[0][1] == pytest.approx(0.0)
    assert matrix[1][0] == pytest.approx(0.0)


# Check disabled exchanges are absent from the elastic tuning surface
def test_eikonal_excludes_disabled_exchanges(source_model):
    general = json.loads(json.dumps(eikonal._load_general_data()))
    exchange = general["PARAM_SOFT"]["MODEL"]["single"]["EXCHANGE"]
    exchange["O"]["on"] = False
    exchange["O3g"]["on"] = False
    source_model("GENERAL.json", general)

    param, aux = eikonal.setup(model="single", pFF="DPOW")
    assert not any(":EXCHANGE.O." in key for key in param | aux)
    assert not any(":EXCHANGE.O3g." in key for key in param | aux)


# Check eikonal tuning resolves Pomeron and Odderon slots by exchange name
def test_eikonal_reordered_exchanges(source_model):
    general = json.loads(json.dumps(eikonal._load_general_data()))
    soft = general["PARAM_SOFT"]
    model = soft["MODEL"]["single"]
    soft["EXCHANGE_DEF"] = {name: soft["EXCHANGE_DEF"][name] for name in ("R_f2", "P", "O")}
    model["EXCHANGE"] = {name: model["EXCHANGE"][name] for name in ("R_f2", "P", "O")}
    model["EXCHANGE"]["O"]["g"] = [[0.5]]
    source_model("GENERAL.json", general)

    param, aux = eikonal.setup(model="single", pFF="DPOW")
    assert (
        eikonal.descending_coupling_key("SOFT|MODEL.single:EXCHANGE.P.g[0,0]") in param
    )
    assert "SOFT|MODEL.single:EXCHANGE.O.g[0,0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.P.alpha[0]" in param
    assert "SOFT|MODEL.single:EXCHANGE.O.alpha[0]" in param
    assert aux["SOFT|MODEL.single:EXCHANGE.P.ff"] == "P"
    assert aux["SOFT|MODEL.single:EXCHANGE.O.ff"] == "O"


# Check arbitrary channel tuning follows the explicit forward Pomeron selector
def test_eikonal_renamed_pomeron(source_model):
    general = json.loads(json.dumps(eikonal._load_general_data()))
    soft = general["PARAM_SOFT"]
    generated = json.loads(json.dumps(soft["MODEL"]["single"]))
    generated["GW"] = {
        "theta": [0.0 for _ in range(6)],
        "a_c": [1.0, 0.0, 0.0],
    }
    for exchange in generated["EXCHANGE"].values():
        diagonal = float(exchange["g"][0][0])
        exchange["g"] = [
            [diagonal if row == column else 0.0 for column in range(4)] for row in range(4)
        ]
    for bank in generated["FF"].values():
        bank["param"] = [list(bank["param"][0]) for _ in range(4)]
    soft["EXCHANGE_DEF"]["P_aux"] = json.loads(json.dumps(soft["EXCHANGE_DEF"]["P"]))
    generated["EXCHANGE"]["P_aux"] = json.loads(json.dumps(generated["EXCHANGE"]["P"]))
    generated["EXCHANGE"]["P_aux"]["g"] = [
        [8.0 if row == column else 0.0 for column in range(4)] for row in range(4)
    ]
    generated["EIKONAL"]["excitation_exchanges"] = ["P_aux"]
    soft["MODEL"]["generated_four"] = generated
    source_model("GENERAL.json", general)

    param, aux = eikonal.setup(model="generated_four", pFF="DPOW")
    assert eikonal.source_channel_count("generated_four") == 4
    assert eikonal.source_forward_exchange("generated_four") == "P_aux"
    assert {key for key in param if ":GW.theta[" in key} == {
        f"SOFT|MODEL.generated_four:GW.theta[{index}]" for index in range(3)
    }
    for row in range(4):
        selected = eikonal.descending_coupling_key(
            f"SOFT|MODEL.generated_four:EXCHANGE.P_aux.g[{row},{row}]"
        )
        assert selected in param
        assert (
            eikonal.ordered_dpow_key(f"SOFT|MODEL.generated_four:FF.P.param[{row},0]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1])
            in param
        )
    assert not any(":EXCHANGE.P.g[" in key for key in param)
    assert aux["SOFT|MODEL.generated_four:EXCHANGE.P_aux.ff"] == "P"


# Check that unsupported source-card form-factor families are rejected
def test_eikonal_missing_gkernel():
    with pytest.raises(Exception, match="Unknown FF = GKERNEL"):
        eikonal.setup(model="single", pFF="GKERNEL")


# Check that excited-state rotations remain fixed outside the tune parameter space
def test_eikonal_proton_mixing():
    expected = {
        "single": set(),
        "double": {"SOFT|MODEL.double:GW.theta[0]"},
    }

    for model, expected_keys in expected.items():
        param = eikonal.general(model=model, tune_odderon=False)
        theta_keys = {key for key in param if ":GW.theta[" in key}
        assert theta_keys == expected_keys


# Check eikonal domains use canonical couplings, DPOW scales and a free q coordinate
def test_eikonal_canonical_optimizer_coords():
    for model, rows in (("single", 1), ("double", 2)):
        param, _ = eikonal.setup(model=model, pFF="DPOW")
        source = eikonal.source_model_block(model)["EXCHANGE"]
        assert _uniform_bounds(param[f"SOFT|MODEL.{model}:EIKONAL.q"]) == pytest.approx((0.3, 1.0))

        for row in range(rows):
            coupling_key = eikonal.descending_coupling_key(
                f"SOFT|MODEL.{model}:EXCHANGE.P.g[{row},{row}]"
            )
            expected_bounds = (7.0, 10.0) if row == 0 else (0.5, 1.0)
            assert _uniform_bounds(param[coupling_key]) == pytest.approx(expected_bounds)
            odderon_source = source["O"]["g"][row][row]
            assert _uniform_bounds(
                param[f"SOFT|MODEL.{model}:EXCHANGE.O.g[{row},{row}]"]
            ) == pytest.approx((0.5 * odderon_source, 1.5 * odderon_source))
            o3g_source = source["O3g"]["g"][row][row]
            assert _uniform_bounds(
                param[f"SOFT|MODEL.{model}:EXCHANGE.O3g.g[{row},{row}]"]
            ) == pytest.approx((0.5 * o3g_source, 1.5 * o3g_source))

            lower_key = eikonal.ordered_dpow_key(
                f"SOFT|MODEL.{model}:FF.P.param[{row},0]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]
            )
            fraction_key = eikonal.ordered_dpow_key(
                f"SOFT|MODEL.{model}:FF.P.param[{row},1]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]
            )
            assert _uniform_bounds(param[lower_key]) == pytest.approx((0.01, 5.0))
            assert _uniform_bounds(param[fraction_key]) == pytest.approx((0.0, 1.0))

        if rows > 1:
            assert _uniform_bounds(
                param[
                    eikonal.symmetric_coupling_key(
                        f"SOFT|MODEL.{model}:EXCHANGE.O.g[0,1]"
                    )
                ]
            ) == pytest.approx((-1.0, 1.0))
            assert _uniform_bounds(
                param[
                    eikonal.symmetric_coupling_key(
                        f"SOFT|MODEL.{model}:EXCHANGE.O3g.g[0,1]"
                    )
                ]
            ) == pytest.approx((-1.0, 1.0))


# Check canonical optimizer values round-trip through the unchanged flat card schema
def test_eikonal_coords_roundtrip(cdir):
    driver = graniitti_driver.GraniittiDriver()
    param_space, aux_param_space = eikonal.setup(model="double", pFF="DPOW", tune_reggeons=True)
    initial = driver.get_initial_param(
        param_space=param_space,
        aux_param_space=aux_param_space,
        cdir=str(cdir),
    )
    source = _load_tune0_card("GENERAL.json")["PARAM_SOFT"]["MODEL"]["double"]

    coupling_keys = [
        eikonal.descending_coupling_key(f"SOFT|MODEL.double:EXCHANGE.P.g[{row},{row}]")
        for row in range(2)
    ]
    assert eikonal.decode_descending_couplings(
        [initial[key] for key in coupling_keys]
    ) == pytest.approx([source["EXCHANGE"]["P"]["g"][row][row] for row in range(2)])

    for row in range(2):
        scale_keys = [
            eikonal.ordered_dpow_key(f"SOFT|MODEL.double:FF.P.param[{row},{column}]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1])
            for column in range(2)
        ]
        assert eikonal.decode_ordered_dpow_scales(
            *(initial[key] for key in scale_keys), limit=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]
        ) == pytest.approx(source["FF"]["P"]["param"][row][:2])

    config = {
        coupling_keys[0]: 9.0,
        coupling_keys[1]: 0.8,
        eikonal.symmetric_coupling_key("SOFT|MODEL.double:EXCHANGE.O.g[0,1]"): 0.25,
        eikonal.symmetric_coupling_key("SOFT|MODEL.double:EXCHANGE.O3g.g[0,1]"): -0.15,
        eikonal.symmetric_coupling_key("SOFT|MODEL.double:EXCHANGE.R_f2.g[0,1]"): 0.75,
        eikonal.ordered_dpow_key("SOFT|MODEL.double:FF.P.param[0,0]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]): 1.0,
        eikonal.ordered_dpow_key("SOFT|MODEL.double:FF.P.param[0,1]", upper=domains.settings()["bounds"]["eikonal"]["FF"]["DPOW"]["param"][1][1]): 0.5,
        "SOFT|MODEL.double:FF.R.param[0,0]": 4.5,
    }
    tunename = "TUNE_TEST_CANONICAL_EIKONAL"
    _create_card(cdir, tunename, config | aux_param_space, driver)
    generated = _load_generated_card(cdir, tunename, "GENERAL.json")["PARAM_SOFT"]["MODEL"][
        "double"
    ]
    assert [generated["EXCHANGE"]["P"]["g"][row][row] for row in range(2)] == pytest.approx(
        [9.0, 7.2]
    )
    assert generated["EXCHANGE"]["O"]["g"][0][1] == pytest.approx(0.25)
    assert generated["EXCHANGE"]["O"]["g"][1][0] == pytest.approx(0.25)
    assert generated["EXCHANGE"]["O3g"]["g"][0][1] == pytest.approx(-0.15)
    assert generated["EXCHANGE"]["O3g"]["g"][1][0] == pytest.approx(-0.15)
    assert generated["EXCHANGE"]["R_f2"]["g"][0][1] == pytest.approx(0.75)
    assert generated["EXCHANGE"]["R_f2"]["g"][1][0] == pytest.approx(0.75)
    assert generated["FF"]["P"]["param"][0][:2] == pytest.approx([1.0, 3.0])
    assert generated["FF"]["R"]["param"][0] == pytest.approx([4.5])

    decoded = driver._decoded_card_config(
        config=config,
        path=str(cdir / "modeldata" / tunename),
    )["parameters"]
    assert decoded["SOFT|MODEL.double:EXCHANGE.P.g[0,0]"] == pytest.approx(9.0)
    assert decoded["SOFT|MODEL.double:EXCHANGE.P.g[1,1]"] == pytest.approx(7.2)
    assert decoded["SOFT|MODEL.double:EXCHANGE.O.g[0,1]"] == pytest.approx(0.25)
    assert decoded["SOFT|MODEL.double:EXCHANGE.O.g[1,0]"] == pytest.approx(0.25)
    assert decoded["SOFT|MODEL.double:EXCHANGE.O3g.g[0,1]"] == pytest.approx(-0.15)
    assert decoded["SOFT|MODEL.double:EXCHANGE.O3g.g[1,0]"] == pytest.approx(-0.15)
    assert decoded["SOFT|MODEL.double:EXCHANGE.R_f2.g[0,1]"] == pytest.approx(0.75)
    assert decoded["SOFT|MODEL.double:EXCHANGE.R_f2.g[1,0]"] == pytest.approx(0.75)
    assert decoded["SOFT|MODEL.double:FF.P.param[0,0]"] == pytest.approx(1.0)
    assert decoded["SOFT|MODEL.double:FF.P.param[0,1]"] == pytest.approx(3.0)
    assert decoded["SOFT|MODEL.double:FF.R.param[0,0]"] == pytest.approx(4.5)


# Check partial coupling coordinate groups cannot bypass canonical ordering
def test_eikonal_partial_couplings():
    driver = graniitti_driver.GraniittiDriver()
    key = eikonal.descending_coupling_key("SOFT|MODEL.double:EXCHANGE.P.g[0,0]")

    with pytest.raises(Exception, match="Incomplete descending group"):
        driver._prepare_eikonal_canonical_params({key: 9.0})


# Check the elastic tunesetup share datasets and keep topology-specific settings
def test_elastic_tunesetup_shared_surface():
    expected_datasets = [
        "icepack/ELASTIC/ISR_212895/dataset.json",
        "icepack/ELASTIC/ISR_214689/dataset.json",
        "icepack/ELASTIC/STAR_1791591/dataset.json",
        "icepack/ELASTIC/D0_1117021/dataset.json",
        "icepack/ELASTIC/TOTEM_922651/dataset.json",
        "icepack/ELASTIC/TOTEM_1220862/dataset.json",
        "icepack/ELASTIC/TOTEM_1489188/dataset.json",
        "icepack/ELASTIC/TOTEM_1710340/dataset.json",
    ]
    for name, model, nevents in (
        ("tune-elastic-single", "single", 200000),
        ("tune-elastic-double", "double", 200000),
    ):
        module_name = f"tunesetup.graniitti.{name}"
        module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(module_name.rsplit(".", 1)[-1]))

        assert [entry["datacard"] for entry in module.datacards] == expected_datasets
        assert all(entry["nevents"] == nevents for entry in module.datacards)
        assert all(entry["loopscreen"] is True for entry in module.datacards)
        assert all(entry["xsmode"] == "reset" for entry in module.datacards)
        assert all(entry["integrator"] == "VEGAS" for entry in module.datacards)
        assert module.aux_param_space["SOFT|active_model"] == model
        assert module.aux_param_space[f"SOFT|MODEL.{model}:EXCHANGE.P.ff"] == "P"
        assert module.aux_param_space[f"SOFT|MODEL.{model}:EXCHANGE.O3g.ff"] == "O3g"
        assert f"SOFT|MODEL.{model}:EXCHANGE.O3g.alpha[0]" in module.param_space
        assert f"SOFT|MODEL.{model}:EXCHANGE.O3g.alpha[1]" in module.param_space
        assert not any(":EXCHANGE.O.sign" in key for key in module.param_space)
        for exchange in ("R_f2", "R_rho"):
            assert f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.alpha[0]" in module.param_space
            assert f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.alpha[1]" in module.param_space
            assert f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.ff" in module.aux_param_space
            for row in range(eikonal.source_channel_count(model)):
                assert f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.g[{row},{row}]" in module.param_space
                for column in range(row + 1, eikonal.source_channel_count(model)):
                    key = eikonal.symmetric_coupling_key(
                        f"SOFT|MODEL.{model}:EXCHANGE.{exchange}.g[{row},{column}]"
                    )
                    assert key in module.param_space
        assert any(":FF.R.param[" in key for key in module.param_space)
        odderon_transition = eikonal.symmetric_coupling_key(
            f"SOFT|MODEL.{model}:EXCHANGE.O.g[0,1]"
        )
        o3g_transition = eikonal.symmetric_coupling_key(
            f"SOFT|MODEL.{model}:EXCHANGE.O3g.g[0,1]"
        )
        assert (odderon_transition in module.param_space) is (model != "single")
        assert (o3g_transition in module.param_space) is (model != "single")
        for row in range({"single": 1, "double": 2}[model]):
            assert f"SOFT|MODEL.{model}:EXCHANGE.O3g.g[{row},{row}]" in module.param_space
            assert f"SOFT|MODEL.{model}:FF.O3g.param[{row},0]" in module.param_space


# Check real elastic tuning data resolve exactly one fitted model per measurement
@pytest.mark.parametrize("model", ["single", "double", "triple"])
def test_elastic_tunesetup_run_only_fitted_model(model, tmp_path):
    module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(f"tune-elastic-{model}"))
    driver = graniitti_driver.GraniittiDriver()
    driver.init_data(
        run_name=str(tmp_path / "elastic"), datacards=module.datacards,
        obs_module=None, cdir=str(ROOT), pickle_dump=False,
    )
    mc_steer = {"tunesetup_name": f"elastic_{model}", "tune_default": "TUNE0"}
    total = 0
    for index, datacard in enumerate(module.datacards):
        runs = driver._run_cards(
            index=index, datacard=datacard, mc_steer=mc_steer, tunename="TUNE0", cdir=str(ROOT),
        )
        total += len(runs)
        routed = []
        for run in runs:
            assert run.process == "X[EL]<Q>"
            assert run.loopscreen == "true"
            assert f'GENERAL.json:PARAM_SOFT.active_model="{model}"' in run.overrides
            assert sum("PARAM_SOFT.active_model=" in override for override in run.overrides) == 1
            routed.extend(driver._sample_set_indices(index=index, sample_name=run.sample_name))
        assert sorted(routed) == list(range(len(driver.datasets[index]["sets"])))
    assert total == 14


# Check disabled precomputation loads real elastic data and parameters without generation
@pytest.mark.parametrize("model", ["single", "double"])
def test_elastic_head_without_precomputation(model, tmp_path, monkeypatch):
    import runpy
    from types import SimpleNamespace

    from core.tune.runtime import process as iceruntime

    # Fail if head initialization attempts to run the event generator
    def reject_generation(*args, **kwargs):
        pytest.fail("Disabled precomputation launched the generator on the head")

    monkeypatch.setattr(iceruntime, "execute_cmd", reject_generation)
    module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(f"tune-elastic-{model}"))
    driver = graniitti_driver.GraniittiDriver()
    cli = runpy.run_module("core.icetune")
    args = SimpleNamespace(
        precompute=False, phase="fit", backend="ray", ray_init_state=None,
        run_name=str(tmp_path / "elastic"), obs_module=None, cdir=str(ROOT), pickle_dump=False,
    )
    result = cli["initialize_backend_context"](
        args=args, tunesetup=module, simdriver=driver, mc_steer={"tune_default": "TUNE0"},
    )
    assert driver.initialized
    assert result["initial_points"]
    assert result["bootstrapper"] is None
    assert driver.data_covariance_payload is None


# Check that central MP tunes do not refit the elastic SOFT parameters
def test_mp_disabled_elastic_params():
    modules = [
        "tune-mpom-res-con-coh",
        "tune-mpom-res-con-incoh",
    ]

    for name in modules:
        module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(name))
        soft_parameters = {key for key in module.param_space if key.startswith("SOFT|")}
        assert soft_parameters == set()
        assert not any(key.startswith("SOFT|") for key in module.aux_param_space)


# Check the combined tunesetup joins elastic and STAR central GP surfaces
def test_gp_central_elastic_params():
    names = (
        "tunesetup.graniitti.tune-elastic-single",
        "tunesetup.graniitti.tune-gpom-res",
        "tunesetup.graniitti.tune-elastic-gpom",
    )
    elastic, gp, combined = (load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(name.rsplit(".", 1)[-1])) for name in names)

    assert combined.datacards == [*elastic.datacards, *gp.datacards]
    assert set(combined.param_space) == set(elastic.param_space) | set(gp.param_space)
    assert set(combined.aux_param_space) == set(elastic.aux_param_space) | set(gp.aux_param_space)
    assert not set(elastic.param_space) & set(gp.param_space)
    assert not set(elastic.aux_param_space) & set(gp.aux_param_space)
    assert combined.aux_param_space["SOFT|active_model"] == "single"
    assert combined.aux_param_space["REGGE|use_zeta.GP"] is True
    assert not any(key.startswith("CON_GP|") for key in combined.aux_param_space)


def test_prod_decay_params_route_separate_files(cdir):
    tunename = "TUNE_TEST_PRODUCTION"
    xp_p_sub = f"{XP_P_PRODUCTION}:g_ls(3,1)"
    eta_xp = "RES|eta:XP:g_ls(1,1)"
    f2_gp_lead = "RES|f2_1270:GP:g_ls(0,2)"
    f2_gp_sub = "RES|f2_1270:GP:g_ls(2,0)"

    params = {
        "DECAY|225:[321,-321]:BR": 0.11,
        "DECAY|225:[321,-321]:zeta.GP": 1.25,
        tools.magnitude_key(xp_p_sub): 0.35,
        tools.raw_phase_key(xp_p_sub): -0.45,
        tools.complex_re_key("RES|f0_500:MP:g"): 0.91,
        tools.complex_im_key("RES|f0_500:MP:g"): 0.0,
        tools.complex_re_key(eta_xp): 0.7 * math.cos(0.3),
        tools.complex_im_key(eta_xp): 0.7 * math.sin(0.3),
        tools.complex_re_key(f2_gp_lead): 0.2,
        tools.complex_im_key(f2_gp_lead): 0.4,
        tools.magnitude_key(f2_gp_sub): 0.25,
        tools.raw_phase_key(f2_gp_sub): 0.5,
        tools.complex_re_key(XP_PI_POMERON_COUPLING): 0.82,
        tools.complex_im_key(XP_PI_POMERON_COUPLING): -0.31,
        tools.raw_phase_key("RES|f0_500:MP:phi"): 0.37,
        tools.raw_phase_key("RES|eta:XP:phi"): -0.48,
        tools.raw_phase_key("RES|f2_1270:GP:phi"): 0.59,
        tools.raw_phase_key("RES|f0_500:TP:phi"): -0.71,
    }

    orig = _create_card(cdir, tunename, params)
    branching = _load_generated_card(cdir, tunename, "DECAYS.json")
    production = _load_generated_card(cdir, tunename, "CON_XP.json")
    general = _load_generated_card(cdir, tunename, "GENERAL.json")

    tune0_decays = _load_tune0_card("DECAYS.json")
    tune0_production = _load_tune0_card("CON_XP.json")
    f0_res_re, f0_res_im = _tune0_resonance_channel("f0_500", "MP")
    xcon_pi_re, xcon_pi_im = _tune0_continuum_channel("XP", 1)

    assert orig["DECAY|225:[321,-321]:BR"] == pytest.approx(tune0_decays["225"]["[321,-321]"]["BR"])
    assert orig["DECAY|225:[321,-321]:zeta.GP"] == pytest.approx(
        tune0_decays["225"]["[321,-321]"]["zeta"]["GP"]
    )
    tune0_g_ls = tune0_production["995"]["[2212,2212]"]["opposite"]["g_ls"]
    assert orig[tools.magnitude_key(xp_p_sub)] == pytest.approx(tune0_g_ls[1][2])
    assert orig[tools.raw_phase_key(xp_p_sub)] == pytest.approx(tune0_g_ls[1][3])
    assert orig[tools.complex_re_key("RES|f0_500:MP:g")] == pytest.approx(f0_res_re)
    assert orig[tools.complex_im_key("RES|f0_500:MP:g")] == pytest.approx(f0_res_im)
    assert orig[tools.complex_re_key(XP_PI_POMERON_COUPLING)] == pytest.approx(xcon_pi_re)
    assert orig[tools.complex_im_key(XP_PI_POMERON_COUPLING)] == pytest.approx(xcon_pi_im)
    for resonance_name, model in (
        ("f0_500", "MP"),
        ("eta", "XP"),
        ("f2_1270", "GP"),
        ("f0_500", "TP"),
    ):
        key = tools.raw_phase_key(f"RES|{resonance_name}:{model}:phi")
        tune0 = _load_tune0_card(f"RES/{resonance_name}.json")
        assert orig[key] == pytest.approx(tune0["PARAM_RES"]["MODELS"][model]["phi"])

    assert branching["225"]["[321,-321]"]["BR"] == pytest.approx(0.11)
    assert branching["225"]["[321,-321]"]["zeta"]["GP"] == pytest.approx(1.25)
    assert production["995"]["[2212,2212]"]["opposite"]["g_ls"][0][2:] == pytest.approx(
        tune0_g_ls[0][2:]
    )
    assert production["995"]["[2212,2212]"]["opposite"]["g_ls"][1][2:] == pytest.approx(
        [0.35, -0.45]
    )
    assert production["995"]["[211,211]"]["opposite"]["g_ls"][0][2] == pytest.approx(
        math.hypot(0.82, -0.31)
    )
    assert production["995"]["[211,211]"]["opposite"]["g_ls"][0][3] == pytest.approx(
        tools.canonicalize_phase(math.atan2(-0.31, 0.82))
    )
    assert general["PARAM_REGGE"]["PARAM_CON"]["XP"] == _load_tune0_card("GENERAL.json")["PARAM_REGGE"]["PARAM_CON"]["XP"]

    res_f0_500 = _load_generated_card(cdir, tunename, "RES/f0_500.json")
    res_eta = _load_generated_card(cdir, tunename, "RES/eta.json")
    res_f2 = _load_generated_card(cdir, tunename, "RES/f2_1270.json")

    assert _active_model_block(res_f0_500, "MP")["g"] == pytest.approx([0.91, 0.0])
    eta_xres = _active_model_block(res_eta, "XP")
    assert eta_xres["g_ls"][0][2:] == pytest.approx([0.7, 0.3])
    assert _active_model_block(res_f2, "GP")["g_ls"][0][2:] == pytest.approx(
        [math.hypot(0.2, 0.4), tools.canonicalize_phase(math.atan2(0.4, 0.2))]
    )
    assert _active_model_block(res_f2, "GP")["g_ls"][1][2:] == pytest.approx([0.25, 0.5])
    assert res_f0_500["PARAM_RES"]["MODELS"]["MP"]["phi"] == pytest.approx(0.37)
    assert res_eta["PARAM_RES"]["MODELS"]["XP"]["phi"] == pytest.approx(-0.48)
    assert res_f2["PARAM_RES"]["MODELS"]["GP"]["phi"] == pytest.approx(0.59)
    assert res_f0_500["PARAM_RES"]["MODELS"]["TP"]["phi"] == pytest.approx(-0.71)


def test_xp_gp_coords_route_selected_rows(cdir):
    tunename = "TUNE_TEST_RELATIVE_SHAPES"
    xp_row = "RES|f0_500:XP:g_ls(2,2)"
    gp_row = "RES|f2_1270:GP:g_ls(2,0)"
    params = {
        tools.complex_re_key(xp_row): 0.63,
        tools.complex_im_key(xp_row): -0.20,
        tools.complex_re_key(gp_row): 0.40,
        tools.complex_im_key(gp_row): 0.30,
    }

    _create_card(cdir, tunename, params)
    res_f0_500 = _load_generated_card(cdir, tunename, "RES/f0_500.json")
    res_f2 = _load_generated_card(cdir, tunename, "RES/f2_1270.json")

    xp_block = _active_model_block(res_f0_500, "XP")
    gp_block = _active_model_block(res_f2, "GP")
    assert xp_block["g_ls"][1][2:] == pytest.approx(
        [math.hypot(0.63, -0.20), math.atan2(-0.20, 0.63)]
    )
    assert gp_block["g_ls"][1][2:] == pytest.approx([0.5, math.atan2(0.30, 0.40)])


# Check incoherent simplex angles build the symmetric diagonal density matrix
def test_incoherent_rho_populations(cdir):
    base = "RES|f2_1270:MP:polarization.rho_population"
    angles = [0.61, 0.83]
    params = {tools.simplex_theta_key(base, index): value for index, value in enumerate(angles)}
    params["RES|f2_1270:MP:polarization.mode"] = "rho"
    tunename = "TUNE_TEST_RHO_SIMPLEX"
    _create_card(cdir, tunename, params)
    block = _active_model_block(_load_generated_card(cdir, tunename, "RES/f2_1270.json"), "MP")
    polarization = block["polarization"]
    probabilities = tools.simplex_probabilities_from_angles(angles)
    side = 0.5 * np.asarray(probabilities[1:])
    expected = np.concatenate((side[::-1], probabilities[:1], side))
    np.testing.assert_allclose(np.diag(polarization["rho_mag"]), expected)
    np.testing.assert_allclose(
        np.asarray(polarization["rho_mag"]) - np.diag(np.diag(polarization["rho_mag"])),
        0.0,
    )
    np.testing.assert_allclose(polarization["rho_phase"], 0.0)
    assert polarization["mode"] == "rho"
    assert polarization["a_Jz"]


# Check direct GP steering updates one independent photoproduction coupling
def test_direct_gp_helicity_update_photo_tables(cdir):
    tunename = "TUNE_TEST_HELICITY_COMPACT"

    params = {
        tools.complex_re_key("RES|rho_770:GP:helicity(-1,0)"): 0.3,
        tools.complex_im_key("RES|rho_770:GP:helicity(-1,0)"): -0.3,
    }

    _create_card(cdir, tunename, params)
    res_rho = _load_generated_card(cdir, tunename, "RES/rho_770.json")

    x_helicity = _active_model_block(res_rho, "XP")["helicity"]
    g_helicity = _active_model_block(res_rho, "GP")["helicity"]
    tune0_rho = _load_tune0_card("RES/rho_770.json")
    tune0_x_helicity = _active_model_block(tune0_rho, "XP")["helicity"]
    tune0_g_helicity = _active_model_block(tune0_rho, "GP")["helicity"]

    assert x_helicity == tune0_x_helicity
    assert [row[:2] for row in g_helicity] == [row[:2] for row in tune0_g_helicity]
    expected_pair = [math.hypot(0.3, -0.3), math.atan2(-0.3, 0.3)]
    assert g_helicity[0][2:] == pytest.approx(expected_pair)


# Check XP tuning writes independent absolute LS couplings
def test_xp_complex_roundtrip_absolute_ls_couplings(cdir):
    driver = graniitti_driver.GraniittiDriver()
    param_space, _ = resonance.setup(
        model="XP",
        resonances=["chi_c2"],
        production={"mode": "mag_phase_cartesian"},
    )
    params = driver.get_initial_param(param_space, {}, str(cdir))

    first = "RES|chi_c2:XP:g_ls(0,2)"
    params[tools.complex_re_key(first)] = 0.2 * math.cos(-0.31)
    params[tools.complex_im_key(first)] = 0.2 * math.sin(-0.31)
    second = "RES|chi_c2:XP:g_ls(2,0)"
    params[tools.complex_re_key(second)] = 0.2 * math.cos(0.47)
    params[tools.complex_im_key(second)] = 0.2 * math.sin(0.47)
    tunename = "TUNE_TEST_XP_PHASE_GAUGE"
    _create_card(cdir, tunename, params, driver)
    card = _load_generated_card(cdir, tunename, "RES/chi_c2.json")
    block = _active_model_block(card, "XP")
    assert block["g_ls"][0][2:] == pytest.approx([0.2, -0.31])
    assert block["g_ls"][1][2:] == pytest.approx([0.2, 0.47])
    assert "phase" not in block


# Check half-integer GP helicity names round trip through the continuum card
def test_half_integer_gp_helicity_round_trip(cdir):
    source_path = cdir / "modeldata" / "TUNE0" / "CON_GP.json"
    source = graniitti_driver.json5.loads(source_path.read_text(encoding="utf-8"))
    block = source["990"]["[2212,2212]"]["same"]
    block.clear()
    block.update(basis="crossed_helicity", CP=[True, True], helicity=[[-0.5, -0.5, 0, 1.0, 0.0]])
    source_path.write_text(json.dumps(source), encoding="utf-8")
    tunename = "TUNE_TEST_GP_ALPHA_LS"
    row_name = domains.angular_row_name("helicity", [-0.5, -0.5, 0, 1.0, 0.0])
    assert row_name == "helicity(-1/2,-1/2,0)"
    base = f"CON_GP|990:[2212,2212]/same:{row_name}"
    params = {
        tools.magnitude_key(base): 0.52,
        tools.raw_phase_key(base): -0.25,
    }

    _create_card(cdir, tunename, params)
    production = _load_generated_card(cdir, tunename, "CON_GP.json")

    helicity = production["990"]["[2212,2212]"]["same"]["helicity"]
    tune0_helicity = block["helicity"]
    assert [row[:3] for row in helicity] == [row[:3] for row in tune0_helicity]
    selected = next(row for row in helicity if row[:3] == [-0.5, -0.5, 0])
    assert selected[3] == pytest.approx(0.52)
    assert selected[4] == pytest.approx(-0.25)


# Check GP scalar m couplings remain independently tunable
def test_gp_scalar_m_coupling_round_trip(cdir, source_model):
    source = _load_tune0_card("CON_GP.json")
    defaults = [[0, 0, -2, 1.0, 0.0], [0, 0, -1, 2.0, -math.pi], [0, 0, 0, 3.0, 0.0]]
    source["990"]["[321,321]"]["opposite"]["helicity"] = defaults
    source_model("CON_GP.json", source)
    params, _ = continuum.setup(final_state_pdgs=[[211, -211], [321, -321]], 
        REGGE=False,
        continuum_models="GP",
        production_mode="magnitude",
        production_relative_range=0.25,
        equal_gp_m=False,
    )
    coupling_keys = sorted(
        key
        for key in params
        if key.startswith("CON_GP|990:") and ":helicity(0,0," in key and key.endswith("@MAG")
    )
    assert len(coupling_keys) == sum(
        row[-2] > 0.0 for pair in ("[211,211]", "[321,321]")
        for row in source["990"][pair]["opposite"]["helicity"]
    )

    values = graniitti_driver.GraniittiDriver().get_initial_param(params, {}, str(cdir))
    for key in coupling_keys:
        magnitude = values[key]
        assert _uniform_bounds(params[key]) == pytest.approx((0.75 * magnitude, 1.25 * magnitude))

    selected = next(key for key in coupling_keys if ":[321,321]/opposite:helicity(0,0,-1)" in key)
    values[selected] = 3.25
    _create_card(cdir, "TUNE_TEST_GP_EQUAL_M", values)
    exchange = selected.split("|", maxsplit=1)[1].split(":", maxsplit=1)[0]
    pair_key = selected.split(":", maxsplit=2)[1].split("/", maxsplit=1)[0]
    rows = _load_generated_card(cdir, "TUNE_TEST_GP_EQUAL_M", "CON_GP.json")[exchange][
        pair_key
    ]["opposite"]["helicity"]
    assert [row[3] for row in rows] == pytest.approx([defaults[0][3], 3.25, defaults[2][3]])
    assert [row[4] for row in rows] == pytest.approx([0.0, -math.pi, 0.0])


# Check GP continuum LS rows retain their explicit Regge projection
def test_gp_con_ls_projection_roundtrip(cdir):
    source_path = cdir / "modeldata" / "TUNE0" / "CON_GP.json"
    with source_path.open(encoding="utf-8") as stream:
        source = load_json_file(stream.name, loader=graniitti_driver.json5.load)
    block = source["990"]["[211,211]"]["opposite"]
    block["basis"] = "crossed_ls"
    block.pop("helicity")
    block["g_ls"] = [[0, 0, 0, 1.0, 0.1], [2, 0, 2, 0.5, -0.2]]
    source_path.write_text(json.dumps(source), encoding="utf-8")

    first = "CON_GP|990:[211,211]/opposite:g_ls(0,0,0)"
    second = "CON_GP|990:[211,211]/opposite:g_ls(2,0,2)"
    params = {
        tools.magnitude_key(first): 0.7,
        tools.raw_phase_key(first): -0.4,
        tools.complex_re_key(second): 0.3,
        tools.complex_im_key(second): 0.4,
    }
    tunename = "TUNE_TEST_GP_LS_PROJECTION"
    _create_card(cdir, tunename, params)
    rows = _load_generated_card(cdir, tunename, "CON_GP.json")["990"]["[211,211]"][
        "opposite"
    ]["g_ls"]

    assert [row[:3] for row in rows] == [[0, 0, 0], [2, 0, 2]]
    assert rows[0][3:] == pytest.approx([0.7, -0.4])
    assert rows[1][3:] == pytest.approx([0.5, math.atan2(0.4, 0.3)])


# Check photon GP continuum rows keep the fixed-spin four-column schema
def test_gp_photon_con_helicity_roundtrip(cdir):
    source_path = cdir / "modeldata" / "TUNE0" / "CON_GP.json"
    with source_path.open(encoding="utf-8") as stream:
        source = load_json_file(stream.name, loader=graniitti_driver.json5.load)
    source["22"] = {
        "[211,211]": {
            "opposite": {
                "basis": "crossed_helicity",
                "CP": [True, True],
                "helicity": [[0, 0, 1.0, 0.1]],
            }
        }
    }
    source_path.write_text(json.dumps(source), encoding="utf-8")

    base = "CON_GP|22:[211,211]/opposite:helicity(0,0)"
    params = {
        tools.magnitude_key(base): 0.6,
        tools.raw_phase_key(base): -0.3,
    }
    tunename = "TUNE_TEST_GP_PHOTON_HELICITY"
    _create_card(cdir, tunename, params)
    row = _load_generated_card(cdir, tunename, "CON_GP.json")["22"]["[211,211]"][
        "opposite"
    ]["helicity"][0]

    assert row == pytest.approx([0, 0, 0.6, -0.3])


# Check explicit production, transfer and decay form factors round trip independently
def test_explicit_res_form_factor_roundtrip(cdir):
    driver = graniitti_driver.GraniittiDriver()
    params = {
        "RES|f0_500:GP:FF_prod.Lambda2": 3.25,
        "RES|f0_500:TP:FF_transfer.LambdaInv2": 2.5,
        "DECAY|9000221:[211,-211]:FF_decay.TP.Lambda2": 2.75,
    }
    source = _load_tune0_card("RES/f0_500.json")["PARAM_RES"]["MODELS"]
    decays = _load_tune0_card("DECAYS.json")
    initial = driver.get_initial_param(params, {}, str(cdir))
    assert initial == pytest.approx({
        "RES|f0_500:GP:FF_prod.Lambda2": source["GP"]["FF_prod"]["Lambda2"],
        "RES|f0_500:TP:FF_transfer.LambdaInv2": 1.0 / source["TP"]["FF_transfer"]["Lambda2"],
        "DECAY|9000221:[211,-211]:FF_decay.TP.Lambda2": decays["9000221"]["[211,-211]"]["FF_decay"]["TP"]["Lambda2"],
    })
    _create_card(cdir, "EXPLICIT_FF", params, driver)
    tuned = _load_generated_card(cdir, "EXPLICIT_FF", "RES/f0_500.json")["PARAM_RES"]["MODELS"]
    decay = _load_generated_card(cdir, "EXPLICIT_FF", "DECAYS.json")
    assert tuned["GP"]["FF_prod"]["Lambda2"] == pytest.approx(3.25)
    assert tuned["TP"]["FF_transfer"]["Lambda2"] == pytest.approx(0.4)
    assert decay["9000221"]["[211,-211]"]["FF_decay"]["TP"]["Lambda2"] == pytest.approx(2.75)
    assert tuned["XP"] == source["XP"]


# Check exponential resonance transfer steering preserves its declared form
def test_exponential_res_transfer_roundtrip(cdir):
    params = resonance.res_form_factors("TP", resonances=["f1_1420"], fields={"FF_transfer"})
    key = "RES|f1_1420:TP:FF_transfer.b"
    assert set(params) == {key}
    source = _load_tune0_card("RES/f1_1420.json")["PARAM_RES"]["MODELS"]["TP"]["FF_transfer"]
    _create_card(cdir, "EXP_TRANSFER", {key: 1.2 * source["b"]})
    tuned = _load_generated_card(cdir, "EXP_TRANSFER", "RES/f1_1420.json")
    assert tuned["PARAM_RES"]["MODELS"]["TP"]["FF_transfer"] == {**source, "b": pytest.approx(1.2 * source["b"])}


# Check manually supplied inactive and zero transfer parameters are rejected
@pytest.mark.parametrize(
    ("parameter", "value"),
    [
        ("RES|f0_500:GP:FF_transfer.LambdaInv2", 2.0),
        ("RES|f0_500:TP:FF_transfer.LambdaInv2", 0.0),
        ("CON_GP|990:[211,211]:FF_transfer.LambdaInv2", 2.0),
    ],
)
def test_rejects_inactive_or_zero_transfer_params(parameter, value, cdir):
    if parameter.startswith("RES|f0_500:GP:"):
        path = cdir / "modeldata/TUNE0/RES/f0_500.json"
        card = load_json_file(path, loader=graniitti_driver.json5.load)
        card["PARAM_RES"]["MODELS"]["GP"]["FF_transfer"] = {"type": "none"}
        path.write_text(json.dumps(card), encoding="utf-8")
    if parameter.startswith("CON_GP|"):
        path = cdir / "modeldata/TUNE0/CON_GP.json"
        card = graniitti_driver.json5.loads(path.read_text(encoding="utf-8"))
        card["990"]["[211,211]"]["FF_transfer"] = {"type": "none"}
        path.write_text(json.dumps(card), encoding="utf-8")
    with pytest.raises(ValueError):
        _create_card(cdir, "TUNE_TEST_INVALID_TRANSFER", {parameter: value})


# Check active inverse transfer scales write power forms
def test_con_transfer_inverse_scale_roundtrip(cdir):
    driver = graniitti_driver.GraniittiDriver()
    source_path = cdir / "modeldata" / "TUNE0" / "CON_GP.json"
    source = graniitti_driver.json5.loads(source_path.read_text(encoding="utf-8"))
    source["990"]["[211,211]"]["FF_transfer"] = {
        "type": "power", "norm": "zero", "Lambda2": 0.5, "n": 1.0
    }
    source["990"]["[321,321]"]["FF_transfer"] = {
        "type": "power", "norm": "zero", "Lambda2": 2.0, "n": 1.0
    }
    source_path.write_text(json.dumps(source), encoding="utf-8")
    params = {
        "CON_GP|990:[211,211]:FF_transfer.LambdaInv2": 2.5,
        "CON_GP|990:[321,321]:FF_transfer.LambdaInv2": 4.0,
    }
    tunename = "TUNE_TEST_GP_TRANSFER"

    initial = driver.get_initial_param(params, {}, str(cdir))
    expected = {
        "CON_GP|990:[211,211]:FF_transfer.LambdaInv2": 2.0,
        "CON_GP|990:[321,321]:FF_transfer.LambdaInv2": 0.5,
    }
    assert initial == pytest.approx(expected)
    _create_card(cdir, tunename, params, driver)
    production = _load_generated_card(cdir, tunename, "CON_GP.json")
    pion = production["990"]["[211,211]"]["FF_transfer"]
    kaon = production["990"]["[321,321]"]["FF_transfer"]
    assert pion == {
        "type": "power",
        "norm": "zero",
        "Lambda2": pytest.approx(0.4),
        "n": pytest.approx(1.0),
    }
    assert kaon == {
        "type": "power",
        "norm": "zero",
        "Lambda2": pytest.approx(0.25),
        "n": pytest.approx(1.0),
    }


def test_probe_tune_name_unique():
    driver = graniitti_driver.GraniittiDriver()
    tunename_a = driver.trial_tunename(trial_id="probe", pid=143, purpose="probe")
    tunename_b = driver.trial_tunename(trial_id="probe", pid=143, purpose="probe")

    assert tunename_a != tunename_b


def test_default_tmp_token_is_unique():
    token_a = icetune.default_tmp_token()
    token_b = icetune.default_tmp_token()

    assert token_a != token_b


def test_failed_card_cleanup(monkeypatch, cdir):
    driver = graniitti_driver.GraniittiDriver()
    tunename = "TUNE_icetune_stage_fail_pid143"
    original_copyfile = graniitti_driver.shutil.copyfile

    def flaky_copyfile(src, dst, *args, **kwargs):
        if Path(src).name == "GENERAL.json":
            raise OSError("simulated copy failure")
        return original_copyfile(src, dst, *args, **kwargs)

    monkeypatch.setattr(graniitti_driver.shutil, "copyfile", flaky_copyfile)

    with pytest.raises(OSError, match="simulated copy failure"):
        driver.create_steering_card(
            param_space={},
            tunename=tunename,
            cdir=str(cdir),
            tune_default="TUNE0",
        )

    assert not (cdir / "tmp" / tunename).exists()
    leftover = list((cdir / "tmp").glob(f".{tunename}.tmp.*"))
    assert leftover == []


def test_stale_tune_cleanup(cdir):
    tunename = "TUNE_icetune_stale_pid143"
    stale_dir = cdir / "tmp" / tunename
    stale_stage = cdir / "tmp" / f".{tunename}.tmp.stale"
    stale_dir.mkdir(parents=True)
    stale_stage.mkdir(parents=True)
    (stale_dir / "stale.txt").write_text("old", encoding="utf-8")
    (stale_stage / "stale.txt").write_text("old-stage", encoding="utf-8")

    _create_card(cdir, tunename, {"DECAY|225:[321,-321]:zeta.GP": 1.25})

    assert (cdir / "tmp" / tunename / "GENERAL.json").exists()
    assert not (cdir / "tmp" / tunename / "stale.txt").exists()
    assert not stale_stage.exists()


# Check Cayley scalar phases and Cartesian observable couplings route correctly
def test_encoded_obs_phases_couplings_route_cards(cdir):
    tunename = "TUNE_TEST_PHASE_ENCODING"

    params = {
        tools.phase_u_key("DECAY|225:[321,-321]:zeta.GP"): 0.4,
        tools.phase_v_key("DECAY|225:[321,-321]:zeta.GP"): -0.7,
        tools.phase_u_key("RES|eta:XP:g_ls(1,1)"): 0.3,
        tools.phase_v_key("RES|eta:XP:g_ls(1,1)"): 0.9,
        tools.complex_re_key(MP_PI_POMERON_COUPLING): -0.4,
        tools.complex_im_key(MP_PI_POMERON_COUPLING): 0.2,
    }

    orig = _create_card(cdir, tunename, params)
    branching = _load_generated_card(cdir, tunename, "DECAYS.json")
    res_eta = _load_generated_card(cdir, tunename, "RES/eta.json")
    production = _load_generated_card(cdir, tunename, "CON_MP.json")

    assert branching["225"]["[321,-321]"]["zeta"]["GP"] == pytest.approx(
        tools.decode_phase_components(0.4, -0.7)
    )
    phase_eta = tools.decode_phase_components(0.3, 0.9)
    mag_con, phase_con = tools.decode_cartesian(-0.4, 0.2)
    assert _active_model_block(res_eta, "XP")["g_ls"][0][3] == pytest.approx(phase_eta)
    assert production["995"]["[211,211]"]["opposite"]["g"][0] == pytest.approx(mag_con)
    assert production["995"]["[211,211]"]["opposite"]["g"][1] == pytest.approx(phase_con)

    tune0_decays = _load_tune0_card("DECAYS.json")
    default_zeta = tools.encode_phase(tune0_decays["225"]["[321,-321]"]["zeta"]["GP"])
    assert orig[tools.phase_u_key("DECAY|225:[321,-321]:zeta.GP")] == pytest.approx(default_zeta[0])
    assert orig[tools.phase_v_key("DECAY|225:[321,-321]:zeta.GP")] == pytest.approx(default_zeta[1])

    con_pi_re, con_pi_im = _tune0_continuum_channel("MP", 1)
    assert orig[tools.complex_re_key(MP_PI_POMERON_COUPLING)] == pytest.approx(con_pi_re)
    assert orig[tools.complex_im_key(MP_PI_POMERON_COUPLING)] == pytest.approx(con_pi_im)


def test_get_initial_param_encoded_phase_defaults(cdir):
    driver = graniitti_driver.GraniittiDriver()

    param_space = {
        tools.phase_u_key("DECAY|225:[321,-321]:zeta.GP"): 0.0,
        tools.phase_v_key("DECAY|225:[321,-321]:zeta.GP"): 0.0,
        tools.complex_re_key(XP_PI_POMERON_COUPLING): 0.0,
        tools.complex_im_key(XP_PI_POMERON_COUPLING): 0.0,
    }

    orig = driver.get_initial_param(
        param_space=param_space,
        aux_param_space={},
        cdir=str(cdir),
        tune_default="TUNE0",
    )

    tune0_decays = _load_tune0_card("DECAYS.json")
    zeta_default = tools.encode_phase(tune0_decays["225"]["[321,-321]"]["zeta"]["GP"])
    assert orig[tools.phase_u_key("DECAY|225:[321,-321]:zeta.GP")] == pytest.approx(zeta_default[0])
    assert orig[tools.phase_v_key("DECAY|225:[321,-321]:zeta.GP")] == pytest.approx(zeta_default[1])
    xcon_pi_re, xcon_pi_im = _tune0_continuum_channel("XP", 1)
    assert orig[tools.complex_re_key(XP_PI_POMERON_COUPLING)] == pytest.approx(xcon_pi_re)
    assert orig[tools.complex_im_key(XP_PI_POMERON_COUPLING)] == pytest.approx(xcon_pi_im)


# Check publication metadata contains every optimizer value in generated-card coordinates
def test_optimizer_encoding_decode(cdir):
    driver = graniitti_driver.GraniittiDriver()
    tunename = "TUNE_TEST_CARD_CONFIG"
    channel_base = "RES|f0_500:GP:g_ls(0,0)"
    vector_base = "DECAY|9080225:[113,113]:alpha_ls"
    config = {
        tools.raw_phase_key("DECAY|225:[321,-321]:zeta.GP"): 0.25,
        tools.complex_re_key(channel_base): 0.3,
        tools.complex_im_key(channel_base): -0.4,
        "DECAY|9080225:[113,113]:alpha_ls(2,0)@PROJECTIVE": 0.5,
        "DECAY|9080225:[113,113]:alpha_ls(2,2)@PROJECTIVE": -0.25,
        "DECAY|9080225:[113,113]:alpha_ls(4,2)@PROJECTIVE": 0.35,
    }

    _create_card(cdir, tunename, config, driver)
    card_config = driver._decoded_card_config(
        config=config,
        path=str(cdir / "modeldata" / tunename),
    )

    parameters = card_config["parameters"]
    assert parameters["DECAY|225:[321,-321]:zeta.GP"] == pytest.approx(0.25)
    assert parameters["RES|f0_500:GP:g_ls(0,0)@MAG"] == pytest.approx(0.5)
    assert parameters["RES|f0_500:GP:g_ls(0,0)@PHASE"] == pytest.approx(
        math.atan2(-0.4, 0.3)
    )
    expected = tools.projective_vector_from_angles([0.5, -0.25, 0.35])
    table = card_config["tables"][vector_base]
    for row, coefficient in enumerate(expected):
        row_name = domains.angular_row_name("alpha_ls", table["rows"][row])
        physical_base = f"DECAY|9080225:[113,113]:{row_name}"
        assert parameters[tools.magnitude_key(physical_base)] == pytest.approx(abs(coefficient))
        expected_phase = tools.real_phase(coefficient)
        assert parameters[tools.raw_phase_key(physical_base)] == pytest.approx(expected_phase)

    assert table["basis"] == "alpha_ls"
    assert table["representation"] == "signed_real_projective"
    assert table["phase_encoding"] == (
        "0 for nonnegative coefficients and -pi for negative coefficients"
    )
    assert [row[2] for row in table["rows"]] == pytest.approx([abs(value) for value in expected])
    assert sum(row[2] ** 2 for row in table["rows"]) == pytest.approx(1.0)
    assert not any(":g_ls[" in key or ":alpha_ls[" in key for key in parameters)


# Check relative LS norm bounds use the baseline coordinates and preserve directions
@pytest.mark.parametrize("res", ["f0_500", "f0_980", "f0_1500", "f0_1710", "f2_1270", "f2_1525", "f2_1950", "eta"])
@pytest.mark.parametrize("override", [False, True])
@pytest.mark.parametrize("model", ["GP", "XP"])
def test_relative_production_norm_bounds(cdir, res, override, model):
    production = {"mode": "magnitude", "geometry": "projective", "active_only": True}
    absolute, _ = resonance.setup(model=model, resonances=[res], production=production)
    scales = {"default": [0.5, 1.5], res: [0.1, 1.5]} if override else [0.5, 1.5]
    relative, auxiliary = resonance.setup(
        model=model, resonances=[res], production={**production, "relative_bounds": scales})
    initial = graniitti_driver.GraniittiDriver().get_initial_param(relative, auxiliary, str(cdir))
    assert relative.keys() == absolute.keys()
    for key, domain in relative.items():
        if tools.is_projective_angle_key(key):
            assert (domain.lower, domain.upper) == (absolute[key].lower, absolute[key].upper)
        else:
            assert (domain.lower, domain.upper) == pytest.approx(((0.1 if override else 0.5) * initial[key], 1.5 * initial[key]))


# Check MP strengths use relative bounds without changing spin coordinates
@pytest.mark.parametrize("coherence", ["coherent", "incoherent"])
def test_relative_production_mp_bounds(cdir, coherence):
    production = {"mode": "magnitude", "coherence": coherence}
    resonances = ["f2_1270", "f2_1950"]
    absolute, _ = resonance.setup(model="MP", resonances=resonances, production=production)
    relative, auxiliary = resonance.setup(model="MP", resonances=resonances, production={
        **production, "relative_bounds": {"default": [0.5, 1.5], "f2_1950": [0.1, 1.5]}})
    initial = graniitti_driver.GraniittiDriver().get_initial_param(relative, auxiliary, str(cdir))
    assert relative.keys() == absolute.keys()
    for key, domain in relative.items():
        if key.endswith(":g[0]"):
            lower = 0.1 if key.startswith("RES|f2_1950:") else 0.5
            assert (domain.lower, domain.upper) == pytest.approx((lower * initial[key], 1.5 * initial[key]))
        else:
            assert (domain.lower, domain.upper) == (absolute[key].lower, absolute[key].upper)


# Reject invalid relative norm ranges before constructing a fit
@pytest.mark.parametrize("bounds", [[-0.5, 1.5], [1.5, 0.5], [1.0, 1.0], [0.5], [0.5, float("inf")]])
def test_relative_production_norm_invalid(bounds):
    with pytest.raises(ValueError, match="production.relative_bounds"):
        resonance.setup(model="GP", resonances=["f2_1950"], production={
            "mode": "magnitude", "geometry": "projective", "relative_bounds": bounds})


# Reject ambiguous manual and relative production norm bounds
def test_relative_production_norm_conflict():
    with pytest.raises(ValueError, match="not both"):
        resonance.setup(model="GP", resonances=["f2_1950"], production={
            "mode": "magnitude", "geometry": "projective", "relative_bounds": [0.5, 1.5], "max_abs": 6.0})


# Check projective production coordinates preserve the absolute resonance card schema
def test_projective_ls_roundtrip(cdir, ls_model):
    driver = graniitti_driver.GraniittiDriver()
    param_space, _ = resonance.setup(
        model="XP",
        resonances=["f2_1270"],
        production={"mode": "magnitude", "geometry": "projective", "active_only": True},
    )
    param_space = {key: value for key, value in param_space.items() if ":g_ls" in key}
    initial = driver.get_initial_param(param_space, {}, str(cdir))
    norm_key = tools.projective_norm_key("RES|f2_1270:XP:g_ls(0,2)")
    angle_keys = [key for key in param_space if tools.is_projective_angle_key(key)]
    angles = (0.3 * np.sin(np.arange(1, len(angle_keys) + 1))).tolist()
    config = {norm_key: 2.5, **dict(zip(angle_keys, angles, strict=True))}
    source_path = cdir / "modeldata" / "TUNE0" / "RES" / "f2_1270.json"
    source = source_path.read_text(encoding="utf-8")

    tunename = "TUNE_TEST_PROJECTIVE_PRODUCTION"
    original = _create_card(cdir, tunename, config, driver)
    generated = _load_generated_card(cdir, tunename, "RES/f2_1270.json")
    block = generated["PARAM_RES"]["MODELS"]["XP"]["[995,995]"]
    direction = tools.projective_vector_from_angles(angles)
    source_rows = _active_model_block(graniitti_driver.json5.loads(source), "XP")["g_ls"]
    orientation = math.copysign(1.0, math.cos(float(source_rows[0][3])))

    assert source_path.read_text(encoding="utf-8") == source
    assert set(original) == set(initial)
    assert original == pytest.approx(initial)
    assert [row[2] for row in block["g_ls"]] == pytest.approx(
        [2.5 * abs(value) for value in direction]
    )
    assert [row[3] for row in block["g_ls"]] == pytest.approx(
        [tools.real_phase(orientation * value) for value in direction]
    )
    assert all(-math.pi <= row[3] < math.pi for row in block["g_ls"])
    card_config = driver._decoded_card_config(
        config=config,
        path=str(cdir / "modeldata" / tunename),
    )
    vector_base = "RES|f2_1270:XP:g_ls"
    assert vector_base in card_config["tables"]
    assert card_config["tables"][vector_base]["representation"] == "signed_real_projective"
    assert tools.magnitude_key("RES|f2_1270:XP:g_ls(0,2)") in card_config["parameters"]
    assert norm_key not in card_config["parameters"]
    assert not set(angle_keys).intersection(card_config["parameters"])
    restored = driver.get_initial_param(param_space, {}, str(cdir), tune_default=tunename)
    assert restored == pytest.approx(config)


# Check Cartesian resonance couplings are equivalent to norm, phase and projective LS coordinates
def test_gp_cartesian_polar_roundtrip(cdir, ls_model):
    driver = graniitti_driver.GraniittiDriver()
    param_space, _ = resonance.setup(
        model="GP",
        resonances=["f0_980"],
        production={"mode": "mag_phase_cartesian", "geometry": "projective", "active_only": True},
    )
    param_space = {key: value for key, value in param_space.items() if ":g_ls" in key}
    reference = "RES|f0_980:GP:g_ls(0,0)"
    re_key, im_key = tools.projective_coefficient_keys(reference)
    angle_keys = [key for key in param_space if tools.is_projective_angle_key(key)]
    norm = 1.4
    phase = -0.8
    angles = [0.35, -0.25]
    config = {
        re_key: norm * math.cos(phase),
        im_key: norm * math.sin(phase),
        **dict(zip(angle_keys, angles, strict=True)),
    }

    tunename = "TUNE_TEST_CARTESIAN_PROJECTIVE_GP"
    original = _create_card(cdir, tunename, config, driver)
    generated = _load_generated_card(cdir, tunename, "RES/f0_980.json")["PARAM_RES"]["MODELS"]["GP"]
    direction = tools.hemisphere_vector_from_angles(angles)

    assert original == pytest.approx(driver.get_initial_param(param_space, {}, str(cdir)))
    assert generated["phi"] == pytest.approx(phase)
    assert [row[2] for row in generated["[990,990]"]["g_ls"]] == pytest.approx(
        [norm * abs(value) for value in direction]
    )
    assert [row[3] for row in generated["[990,990]"]["g_ls"]] == pytest.approx(
        [tools.real_phase(value) for value in direction]
    )
    restored = driver.get_initial_param(param_space, {}, str(cdir), tune_default=tunename)
    assert restored == pytest.approx(config)


# Check spherical GP coordinates preserve and can reverse the complete LS orientation
def test_gp_spherical_signs_norm(cdir, ls_model):
    driver = graniitti_driver.GraniittiDriver()
    param_space, _ = resonance.setup(
        model="GP",
        resonances=["f0_980"],
        production={"mode": "magnitude", "geometry": "spherical", "active_only": True},
    )
    param_space = {key: value for key, value in param_space.items() if ":g_ls" in key}
    norm_key = tools.projective_norm_key("RES|f0_980:GP:g_ls(0,0)")
    angle_keys = [key for key in param_space if tools.is_spherical_angle_key(key)]
    target = [-0.3657660990230842, 0.13757158886937965, -0.11732778800409273]
    norm = math.sqrt(sum(value * value for value in target))

    for index, coefficients in enumerate((target, [-value for value in target])):
        angles = tools.spherical_angles_from_vector(coefficients)
        config = {norm_key: norm, **dict(zip(angle_keys, angles, strict=True))}
        tunename = f"TUNE_TEST_SPHERICAL_GP_{index}"
        _create_card(cdir, tunename, config, driver)
        rows = _active_model_block(
            _load_generated_card(cdir, tunename, "RES/f0_980.json"), "GP"
        )["g_ls"]
        signed = [
            row[2] if math.cos(tools.canonicalize_phase(row[3])) >= 0.0 else -row[2]
            for row in rows
        ]
        card_config = driver._decoded_card_config(
            config=config,
            path=str(cdir / "modeldata" / tunename),
        )

        assert signed == pytest.approx(coefficients)
        assert card_config["tables"]["RES|f0_980:GP:g_ls"]["representation"] == (
            "signed_real_spherical"
        )

    assert norm_key in param_space
    assert len(angle_keys) == 2
    assert _uniform_bounds(param_space[angle_keys[0]]) == pytest.approx(
        (-tools.PROJECTIVE_ANGLE_LIMIT, tools.PROJECTIVE_ANGLE_LIMIT)
    )
    assert _uniform_bounds(param_space[angle_keys[1]]) == pytest.approx((-math.pi, math.pi))


# Check GP fixes secondary trajectories and exposes physical helicity domains
def test_gp_con_fixed_regge_trajectories():
    params, _ = continuum.setup(final_state_pdgs=[[211, -211], [321, -321]], 
        REGGE=True,
        continuum_models="GP",
        production_mode="magnitude",
    )

    active_model = _load_tune0_card("GENERAL.json")["PARAM_SOFT"]["active_model"]
    prefix = f"SOFT|MODEL.{active_model}"
    expected_trajectory_bounds = {
        f"{prefix}:EXCHANGE.P.alpha[0]": tuple(domains.settings()["bounds"]["eikonal"]["EXCHANGE"]["pomeron"]["alpha"][0]),
        f"{prefix}:EXCHANGE.P.alpha[1]": tuple(domains.settings()["bounds"]["eikonal"]["EXCHANGE"]["pomeron"]["alpha"][1]),
    }
    for key, bounds in expected_trajectory_bounds.items():
        assert _uniform_bounds(params[key]) == pytest.approx(bounds)
    assert not any("EXCHANGE.R_f2.alpha" in key or "EXCHANGE.R_rho.alpha" in key for key in params)
    assert not any(key.startswith("REGGE|") for key in params)

    couplings = {
        key: spec
        for key, spec in params.items()
        if key.startswith("CON_GP|") and key.endswith("@MAG")
    }
    assert couplings
    assert all(":helicity(" in key and ":helicity[" not in key for key in couplings)
    assert all(_uniform_bounds(spec)[0] == 0.0 for spec in couplings.values())
    assert all(_uniform_bounds(spec)[1] > 0.0 for spec in couplings.values())

    rho_params, _ = continuum.setup(final_state_pdgs=[[113, 113]], 
        REGGE=True,
        continuum_models="GP",
        production_mode="none",
    )
    assert f"{prefix}:EXCHANGE.P.alpha[0]" in rho_params
    assert not any("EXCHANGE.R_f2.alpha" in key for key in rho_params)


def test_con_xp_and_gp_tune_direct_prod_rows():
    x_pi_params, _ = continuum.setup(final_state_pdgs=[[211, -211]], 
        REGGE=False,
        continuum_models="XP",
        production_mode="mag_phase_cartesian",
    )
    x_pi_keys = [key for key in x_pi_params if key.startswith("CON_XP|") and ":g_ls(" in key]
    assert x_pi_keys
    assert all(":g_ls(" in key and ":g_ls[" not in key for key in x_pi_keys)
    assert all("[211,211]/opposite:" in key for key in x_pi_keys)
    assert not any(key.startswith("CON_GP|990:") for key in x_pi_params)

    complex_params, _ = continuum.setup(final_state_pdgs=[[211, -211]], 
        REGGE=False,
        continuum_models="XP",
        production_mode="mag_phase_cayley",
    )
    assert any(key.endswith("@MAG") for key in complex_params if key.startswith("CON_XP|"))
    assert any(key.endswith("@PHASE_U") for key in complex_params if key.startswith("CON_XP|"))
    assert any(key.endswith("@PHASE_V") for key in complex_params if key.startswith("CON_XP|"))

    proton_params, _ = continuum.setup(final_state_pdgs=[[2212, -2212]], 
        REGGE=False,
        continuum_models="GP",
        production_mode="mag_phase_cartesian",
    )
    proton_keys = [key for key in proton_params if key.startswith("CON_GP|") and ":helicity(" in key]
    assert proton_keys
    assert all("[2212,2212]/opposite:helicity(" in key for key in proton_keys)
    assert all(":helicity[" not in key for key in proton_keys)
    assert not any(key.startswith("CON_XP|995:") for key in proton_params)
    assert not any(key.startswith("CON_MP|995:") for key in proton_params)

    pi_k_params, _ = continuum.setup(final_state_pdgs=[[211, -211], [321, -321]], 
        REGGE=False,
        continuum_models="GP",
        production_mode="mag_phase_cartesian",
    )
    pi_k_keys = [key for key in pi_k_params if key.startswith("CON_GP|") and ":helicity(" in key]
    assert any("[211,211]/opposite:helicity(" in key for key in pi_k_keys)
    assert any("[321,321]/opposite:helicity(" in key for key in pi_k_keys)


# Configure the real pion card in either supported GP input basis
@pytest.fixture(params=["g_ls", "helicity"])
def gp_pion_basis(cdir, monkeypatch, request):
    model_dir = cdir / "modeldata" / "TUNE0"
    path = model_dir / "CON_GP.json"
    card = graniitti_driver.json5.loads(path.read_text(encoding="utf-8"))
    block = card["990"]["[211,211]"]["opposite"]
    block.clear()
    rows = [[2, 0, -2, 1.7, 0.0], [2, 0, -1, 0.8, -math.pi], [2, 0, 0, 1.2, 0.0]]
    if request.param == "helicity":
        rows = [[0, 0, row[2], row[3] * math.sqrt(2.0 / 3.0), row[4]] for row in rows]
    block.update(basis="crossed_ls" if request.param == "g_ls" else "crossed_helicity", CP=[True, True])
    block[request.param] = rows
    path.write_text(json.dumps(card), encoding="utf-8")
    monkeypatch.setattr(continuum, "model_dir", lambda: model_dir)
    return request.param, block[request.param]


# Check the Pomeron continuum shape reproduces its parity-completed Frobenius normalization
def test_gp_con_projective_norm(cdir, gp_pion_basis):
    field, defaults = gp_pion_basis
    params, _ = continuum.setup(final_state_pdgs=[[211, -211]], 
        REGGE=False,
        FF_offshell=False,
        continuum_models="GP",
        production_mode="magnitude",
        production_geometry="projective",
        production_relative_range=0.5,
    )
    prefix = "CON_GP|990:[211,211]/opposite:"
    norm_key = next(
        key for key in params if key.startswith(prefix) and key.endswith(tools.PROJECTIVE_NORM_SUFFIX)
    )
    angle_keys = [key for key in params if key.startswith(prefix) and tools.is_projective_angle_key(key)]
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(params, {}, str(cdir))

    expected_norm = math.sqrt(sum((1 if int(row[2]) == 0 else 2) * row[3] ** 2 for row in defaults))
    assert initial[norm_key] == pytest.approx(expected_norm)
    assert len(angle_keys) == 2
    trial = {**initial, norm_key: 3.2, angle_keys[0]: 0.31, angle_keys[1]: -0.27}
    tunename = "TUNE_TEST_GP_CON_PROJECTIVE"
    _create_card(cdir, tunename, trial, driver)
    rows = _load_generated_card(cdir, tunename, "CON_GP.json")["990"]["[211,211]"]["opposite"][
        field
    ]
    completed_norm = math.sqrt(2.0 * rows[0][3] ** 2 + 2.0 * rows[1][3] ** 2 + rows[2][3] ** 2)
    restored = driver.get_initial_param(params, {}, str(cdir), tune_default=tunename)

    assert completed_norm == pytest.approx(trial[norm_key])
    assert restored == pytest.approx(trial)


# Check Cartesian GP continuum couplings are equivalent to the complex channel coefficient
def test_gp_con_cartesian_roundtrip(cdir, gp_pion_basis):
    field, defaults = gp_pion_basis
    params, _ = continuum.setup(final_state_pdgs=[[211, -211]], 
        REGGE=False,
        FF_offshell=False,
        continuum_models="GP",
        production_mode="mag_phase_cartesian",
        production_geometry="projective",
        production_relative_range=0.5,
    )
    re_key = next(key for key in params if key.endswith(tools.PROJECTIVE_COEFF_RE_SUFFIX))
    im_key = next(key for key in params if key.endswith(tools.PROJECTIVE_COEFF_IM_SUFFIX))
    angle_keys = [key for key in params if tools.is_projective_angle_key(key)]
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(params, {}, str(cdir))
    norm = 3.0
    phase = 0.7
    trial = {
        **initial,
        re_key: norm * math.cos(phase),
        im_key: norm * math.sin(phase),
        angle_keys[0]: 0.31,
        angle_keys[1]: -0.27,
    }
    tunename = "TUNE_TEST_GP_CON_CARTESIAN"
    _create_card(cdir, tunename, trial, driver)
    rows = _load_generated_card(cdir, tunename, "CON_GP.json")["990"]["[211,211]"]["opposite"][
        field
    ]
    completed_norm = math.sqrt(2.0 * rows[0][3] ** 2 + 2.0 * rows[1][3] ** 2 + rows[2][3] ** 2)
    restored = driver.get_initial_param(params, {}, str(cdir), tune_default=tunename)
    topology = tools.build_parameter_topology(params)

    assert completed_norm == pytest.approx(norm)
    assert tools.canonicalize_phase(rows[0][4]) == pytest.approx(phase)
    assert topology["groups"][0]["kind"] == "cartesian_projective"
    assert restored == pytest.approx(trial)


# Check explicit MP and XP continuum tables preserve complex projective couplings
@pytest.mark.parametrize("model", ["MP", "XP"])
@pytest.mark.parametrize("field,labels", [("g_ls", [[1, 1], [3, 1]]), ("helicity", [[0.5, 0.5], [0.5, -0.5]])])
@pytest.mark.parametrize("mode", ["magnitude", "mag_phase_cartesian"])
def test_continuum_projective_round_trip(cdir, source_model, model, field, labels, mode):
    import torch
    from core.tune.drivers.graniitti.ampfit.coefficients import AmplitudeSteering

    filename = f"CON_{model}.json"
    card = _load_tune0_card(filename)
    block = card["995"]["[2212,2212]"]
    block["opposite"] = {"basis": f"crossed_{'ls' if field == 'g_ls' else field}", "CP": [True, True],
                         field: [[*labels[0], 1.2, 0.0], [*labels[1], 0.8, math.pi]]}
    source_model(filename, card)
    params, _ = continuum.setup(final_state_pdgs=[[2212, -2212]], REGGE=False, FF_offshell=False,
                                continuum_models=model, production_mode=mode, production_geometry="projective")
    base = f"CON_{model}|995:[2212,2212]/opposite:{field}"
    selected = {key: value for key, value in params.items() if key.startswith(base)}
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(selected, {}, str(cdir))
    _create_card(cdir, "PROJECTIVE_INITIAL", initial, driver)
    rows = _load_generated_card(cdir, "PROJECTIVE_INITIAL", filename)["995"]["[2212,2212]"]["opposite"][field]
    np.testing.assert_allclose([row[-2] * np.exp(1j * row[-1]) for row in rows], [1.2, -0.8], atol=1e-12)
    angle = next(key for key in selected if tools.is_projective_angle_key(key))
    trial = {**initial, angle: -0.37}
    norm, phase = 2.1, 0.7 if mode == "mag_phase_cartesian" else 0.0
    if mode == "mag_phase_cartesian":
        re_key = next(key for key in selected if key.endswith(tools.PROJECTIVE_COEFF_RE_SUFFIX))
        im_key = next(key for key in selected if key.endswith(tools.PROJECTIVE_COEFF_IM_SUFFIX))
        trial.update({re_key: norm * math.cos(phase), im_key: norm * math.sin(phase)})
    else:
        trial[next(key for key in selected if tools.is_projective_norm_key(key))] = norm
    _create_card(cdir, "PROJECTIVE_TRIAL", trial, driver)
    rows = _load_generated_card(cdir, "PROJECTIVE_TRIAL", filename)["995"]["[2212,2212]"]["opposite"][field]
    expected = norm * np.asarray(tools.projective_vector_from_angles([trial[angle]])) * np.exp(1j * phase)
    if field == "helicity":
        expected /= math.sqrt(2.0)
    np.testing.assert_allclose([row[-2] * np.exp(1j * row[-1]) for row in rows], expected, atol=1e-12)
    assert driver.get_initial_param(selected, {}, str(cdir), tune_default="PROJECTIVE_TRIAL") == pytest.approx(trial)
    group = tools.build_parameter_topology(selected)["groups"][0]
    assert group["kind"] == ("cartesian_projective" if mode == "mag_phase_cartesian" else "radial_projective")

    # Check the complex couplings and derivatives used by ampfit
    def evaluate(values):
        steering = AmplitudeSteering(driver, cdir / "modeldata/TUNE0", dict(zip(trial, values, strict=True)))
        amplitudes = [steering.coupling(f"{base}[{i},2]", f"{base}[{i},3]", row[-2:]) for i, row in enumerate(rows)]
        return torch.view_as_real(torch.stack(amplitudes))

    values = torch.tensor(list(trial.values()), dtype=torch.float64, requires_grad=True)
    np.testing.assert_allclose(torch.view_as_complex(evaluate(values)).detach().numpy(), expected, atol=1e-12)
    assert torch.autograd.gradcheck(evaluate, (values,), eps=1e-5)


# Check XP roundtrips the active proton production rows through the real driver
def test_xp_con_proton_roundtrip(cdir):
    params, _ = continuum.setup(final_state_pdgs=[[2212, -2212]], 
        REGGE=False,
        continuum_models="XP",
        tune_all_channels=True,
        production_mode="mag_phase_cartesian",
    )
    re_keys = sorted(key for key in params if key.startswith("CON_XP|") and key.endswith("@RE"))
    assert re_keys
    assert all(f"{key[:-3]}@IM" in params for key in re_keys)

    driver = graniitti_driver.GraniittiDriver()
    trial = driver.get_initial_param(params, {}, str(cdir))
    re_key = re_keys[0]
    im_key = f"{re_key[:-3]}@IM"
    trial[re_key] = 0.37
    trial[im_key] = 0.48
    tunename = "TUNE_TEST_XP_PRODUCTION"
    orig = _create_card(cdir, tunename, trial, driver)
    restored = driver.get_initial_param(params, {}, str(cdir), tune_default=tunename)
    assert orig[re_key] != pytest.approx(trial[re_key])
    assert restored[re_key] == pytest.approx(trial[re_key])
    assert restored[im_key] == pytest.approx(trial[im_key])


# Check particle antiparticle pairs route to physical opposite-sector names
@pytest.mark.parametrize(
    ("model", "field"),
    [
        ("MP", "g"),
        ("XP", "g_ls"),
        ("GP", "helicity"),
    ],
)
@pytest.mark.parametrize(
    ("pdg_pair", "pair_key"),
    [
        ((211, -211), "[211,211]"),
        ((-211, 211), "[211,211]"),
        ((321, -321), "[321,321]"),
        ((-321, 321), "[321,321]"),
        ((311, -311), "[311,311]"),
        ((-311, 311), "[311,311]"),
        ((2212, -2212), "[2212,2212]"),
        ((-2212, 2212), "[2212,2212]"),
    ],
)
def test_con_tunes_physical_opposite_sector(model, field, pdg_pair, pair_key):
    params, _ = continuum.setup(
        REGGE=False,
        continuum_models=model,
        production_mode="mag_phase_cartesian",
        final_state_pdgs=[pdg_pair],
    )

    keys = [key for key in params if key.startswith(f"CON_{model}|") and f":{field}" in key]
    assert keys
    assert all(f":{pair_key}/opposite:{field}" in key for key in keys)
    assert not any(f":{pair_key}/same:{field}" in key for key in keys)
    if field in {"g_ls", "helicity"}:
        assert all(f":{field}(" in key and f":{field}[" not in key for key in keys)


# Check Poisson veto tuning is explicit and restricted to supported continuum models
def test_con_poisson_veto_activation():
    with pytest.raises(ValueError, match="not supported by TP"):
        continuum.setup(final_state_pdgs=[[211, -211], [321, -321]], continuum_models="TP", pveto=True)

    params, aux = continuum.setup(final_state_pdgs=[[211, -211], [321, -321], [2212, -2212]], 
        REGGE=False,
        pveto=True,
        continuum_models="XP",
        production_mode="magnitude",
    )

    production = _load_tune0_card("CON_XP.json")
    veto_keys = [key for key in params if key.endswith(":pveto.M0")]
    assert {key.split(":")[1] for key in veto_keys} == {"[211,211]", "[321,321]", "[2212,2212]"}
    assert len(veto_keys) == 3
    for key in veto_keys:
        exchange, pair, _ = key.partition("|")[2].split(":")
        default_m0 = production[exchange][pair]["pveto"]["M0"]
        assert _uniform_bounds(params[key]) == pytest.approx((0.75 * default_m0, 1.50 * default_m0))
        assert _uniform_bounds(params[key.replace(".M0", ".c")]) == pytest.approx((0.0, 1.5))
        assert aux[key.replace(".M0", ".active")] is True


# Include a constant off shell form factor in the allowed exponential slope range
def test_con_exp_form_factor_bounds_include_zero(monkeypatch):
    monkeypatch.setitem(domains.settings()["bounds"]["continuum"]["FF_offshell"]["exp"], "b", [0.0, 1.5])
    params, aux = continuum.form_factor(model="MP", exchange="995", pdg_pair=(211, -211), ff_type="exp")

    assert aux["CON_MP|995:[211,211]:FF_offshell.type"] == "exp"
    assert _uniform_bounds(params["CON_MP|995:[211,211]:FF_offshell.b"]) == pytest.approx((0.0, 1.5))


# Check inactive continuum transfer factors are absent from tuning
def test_con_inactive_transfer_scales(cdir, monkeypatch):
    model_dir = cdir / "modeldata" / "TUNE0"
    for model in ("MP", "XP", "GP"):
        path = model_dir / f"CON_{model}.json"
        card = graniitti_driver.json5.loads(path.read_text(encoding="utf-8"))
        for exchange in card.values():
            for pair in exchange.values():
                pair["FF_transfer"] = {"type": "none"}
        path.write_text(json.dumps(card), encoding="utf-8")
    monkeypatch.setattr(continuum, "model_dir", lambda: model_dir)
    params, _ = continuum.setup(final_state_pdgs=[[211, -211], [321, -321], [2212, -2212]], 
        REGGE=False,
        FF_transfer=True,
        continuum_models=("MP", "XP", "GP"),
        tune_couplings=False,
    )

    keys = [key for key in params if key.endswith(":FF_transfer.LambdaInv2")]
    assert keys == []


# Check that retained continuum bridge factors expose their parameters in natural order
def test_con_bridge_form_factor_params_natural_order():
    power_params, power_aux = continuum.form_factor(
        model="MP", exchange="995", pdg_pair=(211, -211), ff_type="power"
    )
    gkernel_params, gkernel_aux = continuum.form_factor(
        model="MP", exchange="995", pdg_pair=(211, -211), ff_type="gkernel"
    )

    assert power_aux["CON_MP|995:[211,211]:FF_offshell.type"] == "power"
    assert _uniform_bounds(
        power_params["CON_MP|995:[211,211]:FF_offshell.Lambda2"]
    ) == pytest.approx((0.1, 5.0))
    assert _uniform_bounds(power_params["CON_MP|995:[211,211]:FF_offshell.n"]) == pytest.approx(
        (0.1, 3.0)
    )
    assert gkernel_aux["CON_MP|995:[211,211]:FF_offshell.type"] == "gkernel"
    expected_bounds = (
        (0.05, 5.0),
        (0.5, 1.5),
        (0.0, 2.0),
        (0.0, 2.0),
    )
    for field, bounds in zip(("a", "p", "nu", "mu2"), expected_bounds, strict=True):
        key = f"CON_MP|995:[211,211]:FF_offshell.terms[0].{field}"
        assert _uniform_bounds(gkernel_params[key]) == pytest.approx(bounds)


# Check every supported continuum form-factor family remains selectable
@pytest.mark.parametrize(
    ("form_type", "parameter_count"),
    [
        ("exp", 1),
        ("logexp", 2),
        ("power", 2),
        ("orear", 2),
        ("gkernel", 4),
    ],
)
def test_con_form_factor_families_available(
    form_type,
    parameter_count,
):
    params, aux = continuum.form_factor(
        model="MP", exchange="995", pdg_pair=(211, -211), ff_type=form_type
    )

    assert aux["CON_MP|995:[211,211]:FF_offshell.type"] == form_type
    assert len(params) == parameter_count


# Check continuum tuning inherits the active form factor and its initial parameters
@pytest.mark.parametrize("model,exchange", [("MP", "995"), ("XP", "995"), ("GP", "990"), ("GP", "9910")])
def test_con_defaults_to_card_form_factor(cdir, model, exchange):
    production = _load_tune0_card(f"CON_{model}.json")

    param_space = {}
    aux_param_space = {}
    for pdg_pair, pair_key, _entry in (
        ((211, -211), "[211,211]", "[211,-211]"),
        ((321, -321), "[321,321]", "[321,-321]"),
        ((2212, -2212), "[2212,2212]", "[2212,-2212]"),
    ):
        params, aux = continuum.form_factor(model=model, exchange=exchange, pdg_pair=pdg_pair)
        param_space.update(params)
        aux_param_space.update(aux)

        configured_form = production[exchange][pair_key]["FF_offshell"]["type"]
        assert aux[f"CON_{model}|{exchange}:{pair_key}:FF_offshell.type"] == configured_form

    initial = graniitti_driver.GraniittiDriver().get_initial_param(
        param_space=param_space,
        aux_param_space=aux_param_space,
        cdir=str(cdir),
        tune_default="TUNE0",
    )
    bounds = domains.settings()["bounds"]["continuum"]["FF_offshell"]
    for pair_key, _entry in (
        ("[211,211]", "[211,-211]"),
        ("[321,321]", "[321,-321]"),
        ("[2212,2212]", "[2212,-2212]"),
    ):
        form = production[exchange][pair_key]["FF_offshell"]
        key = f"CON_{model}|{exchange}:{pair_key}:FF_offshell.b"
        assert initial[key] == pytest.approx(form["b"])
        assert _uniform_bounds(param_space[key]) == pytest.approx(bounds[form["type"]]["b"])
        if form["type"] == "logexp":
            key = key.removesuffix("b") + "Lambda2"
            assert initial[key] == pytest.approx(form["Lambda2"])
            assert _uniform_bounds(param_space[key]) == pytest.approx(bounds[form["type"]]["Lambda2"])


# Check that removed continuum form-factor aliases fail at tune construction
def test_con_removed_form_factors_rejected():
    for removed in ("exponential", "QEXP", "KAPPA", "STEXP", "SPECMIX"):
        with pytest.raises(Exception, match="Unknown ff_type"):
            continuum.form_factor(model="MP", exchange="995", pdg_pair=(211, -211), ff_type=removed)


def test_con_kaon_proton_update(cdir):
    tunename = "TUNE_TEST_CONTINUUM_KP"

    params = {
        tools.complex_re_key(MP_K_POMERON_COUPLING): 0.33,
        tools.complex_im_key(MP_K_POMERON_COUPLING): -0.44,
        tools.complex_re_key(XP_P_POMERON_COUPLING): -0.55,
        tools.complex_im_key(XP_P_POMERON_COUPLING): 0.66,
    }

    orig = _create_card(cdir, tunename, params)
    con_mp = _load_generated_card(cdir, tunename, "CON_MP.json")
    con_xp = _load_generated_card(cdir, tunename, "CON_XP.json")

    con_k_mag, con_k_phase = tools.decode_cartesian(0.33, -0.44)
    xcon_p_mag, xcon_p_phase = tools.decode_cartesian(-0.55, 0.66)

    assert con_mp["995"]["[321,321]"]["opposite"]["g"] == pytest.approx([con_k_mag, con_k_phase])
    assert con_xp["995"]["[2212,2212]"]["opposite"]["g_ls"][0][2:] == pytest.approx(
        [xcon_p_mag, xcon_p_phase]
    )

    con_k_re, con_k_im = _tune0_continuum_channel("MP", 2)
    xcon_p_re, xcon_p_im = _tune0_continuum_channel("XP", 0)
    assert orig[tools.complex_re_key(MP_K_POMERON_COUPLING)] == pytest.approx(con_k_re)
    assert orig[tools.complex_im_key(MP_K_POMERON_COUPLING)] == pytest.approx(con_k_im)
    assert orig[tools.complex_re_key(XP_P_POMERON_COUPLING)] == pytest.approx(xcon_p_re)
    assert orig[tools.complex_im_key(XP_P_POMERON_COUPLING)] == pytest.approx(xcon_p_im)


def test_con_form_pdg_mapping():
    params, aux = continuum.setup(final_state_pdgs=[[113, 113], [333, 333]], 
        REGGE=False,
        continuum_models="MP",
    )

    all_keys = set(params) | set(aux)
    assert not any("[*]" in key for key in all_keys)
    assert not any(key.startswith("CON_PAIR|") for key in all_keys)



# Check explicit final states exclusively control form-factor selection
@pytest.mark.parametrize(
    ("pdg_pair", "entry"),
    [
        ((211, -211), "[211,-211]"),
        ((321, -321), "[321,-321]"),
        ((2212, -2212), "[2212,-2212]"),
        ((113, 113), "[113,113]"),
        ((333, 333), "[333,333]"),
    ],
)
def test_con_final_state_form_selection(pdg_pair, entry):
    params, aux = continuum.setup(
        REGGE=False,
        continuum_models="MP",
        production_mode="none",
        final_state_pdgs=[pdg_pair],
    )

    form_keys = {key for key in set(params) | set(aux) if key.startswith("CON_MP|") and "FF_offshell" in key}
    assert form_keys
    pair_key = f"[{abs(pdg_pair[0])},{abs(pdg_pair[1])}]"
    assert all(f":{pair_key}:FF_offshell." in key for key in form_keys)


# Check an unlisted continuum pair cannot use a model channel fallback
def test_con_rejects_unlisted_neutron_pair(source_model):
    general = json.loads(json.dumps(continuum._load_general_data()))
    general["PARAM_REGGE"]["PARAM_CON"]["MP"].pop("[2112,-2112]")
    source_model("GENERAL.json", general)

    with pytest.raises(ValueError, match="No explicit PARAM_CON entry"):
        continuum.setup(
            REGGE=False,
            continuum_models="MP",
            production_mode="magnitude",
            final_state_pdgs=[(2112, -2112)],
        )


# Check reversed continuum pair keys cannot introduce order-dependent selection
def test_con_tuner_rejects_duplicate_pair_keys(source_model):
    general = json.loads(json.dumps(continuum._load_general_data()))
    gp = general["PARAM_REGGE"]["PARAM_CON"]["GP"]
    gp["[-211,211]"] = gp["[211,-211]"]
    source_model("GENERAL.json", general)

    with pytest.raises(ValueError, match="duplicates the order-equivalent pair key"):
        continuum.setup(final_state_pdgs=[[211, -211]], 
            REGGE=False,
            continuum_models="GP",
            production_mode="none",
        )


# Check noncanonical direct production pair keys are rejected
def test_con_tuner_rejects_noncanonical_prod_pair(source_model):
    production = json.loads(json.dumps(continuum._load_production_data("GP")))
    production["990"]["[211, 211]"] = production["990"]["[211,211]"]
    source_model("CON_GP.json", production)

    with pytest.raises(ValueError, match="Invalid CON_GP.json pair key"):
        continuum.setup(final_state_pdgs=[[211, -211]], 
            REGGE=False,
            continuum_models="GP",
            production_mode="none",
        )


# Check the tuner uses the same GP analytic CP declaration as runtime
@pytest.mark.parametrize("mode", ["none", "magnitude"])
def test_con_analytic_cp_required(source_model, mode):
    production = json.loads(json.dumps(continuum._load_production_data("GP")))
    production["990"]["[211,211]"]["opposite"]["CP"] = [False, True]
    source_model("CON_GP.json", production)

    with pytest.raises(ValueError, match=r"must declare CP = \[true, true\]"):
        continuum.setup(final_state_pdgs=[[211, -211]], 
            REGGE=False,
            continuum_models="GP",
            production_mode=mode,
        )


# Check the tuner rejects an ambiguous self-conjugate sector declaration
def test_con_tuner_rejects_self_charge_sectors(source_model):
    production = json.loads(json.dumps(continuum._load_production_data("GP")))
    pair_block = production["990"]["[111,111]"]
    pair_block["same"] = json.loads(json.dumps(pair_block["self"]))
    source_model("CON_GP.json", production)

    with pytest.raises(ValueError, match="self sector cannot coexist"):
        continuum.setup(
            REGGE=False,
            continuum_models="GP",
            production_mode="magnitude",
            final_state_pdgs=[(111, 111)],
        )


# Check generated-card updates reject the same ambiguous declaration
def test_con_card_update_rejects_self_charge_sectors():
    driver = graniitti_driver.GraniittiDriver()
    mother = {"[111,111]": {"self": {}, "same": {}}}

    with pytest.raises(Exception, match="mixes self with same or opposite"):
        driver.production_sector_block(mother, "[111,111]/self")


# Check GP continuum tuning requires an explicit selected final-state pair
def test_con_gcon_rejects_missing_explicit_pair():
    with pytest.raises(ValueError, match="No explicit PARAM_CON entry"):
        continuum.setup(
            REGGE=False,
            continuum_models="GP",
            production_mode="magnitude",
            final_state_pdgs=[(3122, -3122)],
        )

    neutral_params, _ = continuum.setup(
        REGGE=False,
        continuum_models="GP",
        production_mode="magnitude",
        final_state_pdgs=[(111, 111)],
    )
    assert "CON_GP|990:[111,111]/self:helicity(0,0,0)@MAG" in neutral_params


# Check explicit final-state steering cannot select an empty surface
def test_con_empty_final_selection():
    with pytest.raises(ValueError, match="must not be empty"):
        continuum.setup(REGGE=False, final_state_pdgs=[])


# Check continuum steering requires two explicit PDGs
@pytest.mark.parametrize("pdg_pair", [(None, 211), (211, None), None, (None, None)])
def test_con_rejects_null_explicit_final_state(pdg_pair):
    with pytest.raises(ValueError, match="explicit PDG pair|two integer PDGs"):
        continuum.setup(REGGE=False, final_state_pdgs=[pdg_pair])


# Check physical tensor continuum coordinates round trip through model cards
def test_tensor_continuum_round_trip(cdir):
    tp, _ = continuum.setup(
        REGGE=False,
        pveto=False,
        ff_type="power",
        continuum_models="TP",
        tune_couplings=True,
        tune_reggeize=False,
        final_state_pdgs=[(211, -211), (321, -321)],
    )
    scalar_keys = [key for key in tp if key.startswith("CON_TP|")]
    assert scalar_keys
    assert any(":[211,211]:g_tensor" in key for key in scalar_keys)
    assert any(":[321,321]:g_tensor" in key for key in scalar_keys)
    vector_tp, _ = continuum.setup(
        REGGE=False,
        continuum_models="TP",
        tune_couplings=True,
        final_state_pdgs=[(113, 113)],
    )
    assert any(key.endswith(":g_tensor(Gamma0)") for key in vector_tp)
    assert any(key.endswith(":g_tensor(Gamma2)") for key in vector_tp)
    assert not any("g_tensor[" in key for key in (*tp, *vector_tp))

    driver = graniitti_driver.GraniittiDriver()
    key = next(key for key in scalar_keys if ":[211,211]:g_tensor" in key)
    vector_key = next(key for key in vector_tp if key.endswith(":g_tensor(Gamma2)"))
    selected = {key: tp[key], vector_key: vector_tp[vector_key]}
    original = driver.get_initial_param(selected, {}, str(cdir))
    updated = 0.75 * _uniform_bounds(selected[key])[1]
    vector_updated = 0.75 * _uniform_bounds(selected[vector_key])[1]
    tunename = "TUNE_TEST_TP_CONTINUUM"
    defaults = _create_card(cdir, tunename, {key: updated, vector_key: vector_updated}, driver)
    restored = driver.get_initial_param(selected, {}, str(cdir), tune_default=tunename)
    assert defaults == pytest.approx(original)
    assert restored[key] == pytest.approx(updated)
    assert restored[vector_key] == pytest.approx(vector_updated)


# Check the shared density switch selects shape fits without central screening
def test_central_density_screening(tmp_path):
    source = campaign_source("tune-gpom-res-con")
    definition = load_json_file(source, loader=graniitti_driver.json5.load)
    for model in definition["models"]:
        for dataset in model["datasets"]:
            dataset.update(force_density=True, loopscreen=False, xsmode="header")
    path = tmp_path / "density.json"
    path.write_text(json.dumps(definition))
    density = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=path)
    standard = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=source)
    assert all(entry["force_density"] and not entry["loopscreen"] and entry["xsmode"] == "header" for entry in density.datacards)
    assert all(entry["loopscreen"] and entry["xsmode"] == "sample" for entry in standard.datacards)
    assert parameter_space.normalize_param_space(density.param_space) == parameter_space.normalize_param_space(standard.param_space)
    assert density.aux_param_space == standard.aux_param_space


# Check the default GP surfaces tune the f0(980) K+K- decay phase
def test_gp_star_tunes_f0_980_kk_zeta_phase(monkeypatch):
    key = tools.raw_phase_key("DECAY|9010221:[321,-321]:zeta.GP")
    name = "tunesetup.graniitti.tune-gpom-res"
    gp = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source(name.rsplit(".", 1)[-1]))
    assert _uniform_bounds(gp.param_space[key]) == pytest.approx((-math.pi, math.pi))


# Check common family switches can replace the GP resonance parameter selection
def test_central_unselected_params_stay_fixed():
    from core.tune.drivers.graniitti import tunesetup as common
    cards = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source("tune-gpom-res")).datacards
    _, parameters, auxiliary = common.central(
        model="GP", datacards=cards, continuum_options={"FF_offshell": True, "ff_type": "power"},
    )
    assert parameters
    assert all(key.startswith("CON_GP|") and ":FF_offshell." in key for key in parameters)
    assert all(key.startswith("CON_GP|") and ":FF_offshell." in key for key in auxiliary)


# Check one GP parameterization switch controls each RES and CON coupling
def test_gp_cartesian_projective_options(ls_model):
    res_mode = con_mode = ("mag_phase_cartesian", "projective")
    from core.tune.drivers.graniitti import tunesetup as common
    res_coordinate, res_geometry = res_mode
    con_coordinate, con_geometry = con_mode
    cards = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source("tune-gpom-res")).datacards
    cards, parameters, auxiliary = common.central(
        model="GP", datacards=cards, exclude=("rho_770", "phi_1020"),
        production={"mode": res_coordinate, "geometry": res_geometry, "active_only": True},
        phi={"mode": "phase_raw"},
        continuum_options={"tune_couplings": True, "production_mode": con_coordinate, "production_geometry": con_geometry},
    )
    gp = SimpleNamespace(datacards=cards, param_space=parameters, aux_param_space=auxiliary)
    reference = "RES|f0_980:GP:g_ls(0,0)"
    re_key, im_key = tools.projective_coefficient_keys(reference)
    topology = tools.build_parameter_topology(gp.param_space)
    group = next(item for item in topology["groups"] if item["base"] == "RES|f0_980:GP:g_ls")

    assert re_key in gp.param_space
    assert im_key in gp.param_space
    assert tools.raw_phase_key("RES|f0_980:GP:phi") not in gp.param_space
    assert tools.raw_phase_key("RES|rho_770:GP:phi") in gp.param_space
    assert tools.raw_phase_key("RES|phi_1020:GP:phi") in gp.param_space
    assert group["kind"] == "cartesian_projective"
    assert group["parameters"][:2] == [re_key, im_key]

    continuum_groups = [
        item for item in topology["groups"] if item["base"].startswith("CON_GP|")
    ]
    assert continuum_groups
    assert all(item["kind"] == "cartesian_projective" for item in continuum_groups)
    assert all(item["parameters"][0].endswith("@COEFF_RE") for item in continuum_groups)
    assert all(item["parameters"][1].endswith("@COEFF_IM") for item in continuum_groups)
    assert all("weights" in item for item in continuum_groups)


# Check every GP coupling coordinate switch emits its physical topology
@pytest.mark.parametrize(
    ("res_mode", "con_mode", "res_kind", "con_kind"),
    [
        (("magnitude", "projective"), ("magnitude", "projective"), "polar_projective", "radial_projective"),
        (("magnitude", "spherical"), ("magnitude", "direct"), "polar_sphere", None),
        (("magnitude", "direct"), ("magnitude", "direct"), "polar_components", None),
        (("mag_phase_cartesian", "direct"), ("mag_phase_cartesian", "direct"), None, None),
    ],
)
# Check optimizer topology with an explicit continuum table containing multiple active rows
def test_gp_parametrization_topologies(gp_pion_basis, ls_model, res_mode, con_mode, res_kind, con_kind):
    from core.tune.drivers.graniitti import tunesetup as common
    res_coordinate, res_geometry = res_mode
    con_coordinate, con_geometry = con_mode
    cards = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source("tune-gpom-res")).datacards
    cards, parameters, auxiliary = common.central(
        model="GP", datacards=cards, exclude=("rho_770", "phi_1020"),
        production={"mode": res_coordinate, "geometry": res_geometry, "active_only": True},
        phi={"mode": "phase_raw"},
        continuum_options={"tune_couplings": True, "production_mode": con_coordinate, "production_geometry": con_geometry},
    )
    gp = SimpleNamespace(datacards=cards, param_space=parameters, aux_param_space=auxiliary)
    topology = tools.build_parameter_topology(gp.param_space)
    res_groups = [group for group in topology["groups"] if group["base"] == "RES|f0_980:GP:g_ls"]
    con_groups = [group for group in topology["groups"] if group["base"].startswith("CON_GP|")]

    assert ([group["kind"] for group in res_groups] if res_groups else None) == (
        [res_kind] if res_kind is not None else None
    )
    assert ([group["kind"] for group in con_groups] if con_groups else None) == (
        [con_kind] * len(con_groups) if con_kind is not None else None
    )


def test_mp_incoherent_res_basis(cdir):
    incoherent = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=campaign_source("tune-mpom-res-con-incoh"))
    tunename = "TUNE_TEST_INCOHERENT_VERTEX_DEFAULTS"

    _create_card(cdir, tunename, incoherent.aux_param_space)
    block = _active_model_block(_load_generated_card(cdir, tunename, "RES/f2_1270.json"), "MP")

    assert block["basis"] == _active_model_block(_load_tune0_card("RES/f2_1270.json"), "MP")["basis"]
    assert block["polarization"]["mode"] == "rho"
    assert block["polarization"]["a_Jz"]
    assert block["polarization"]["rho_mag"]
    assert block["polarization"]["rho_phase"]


# Production dynamics and spin steering are independent in the actual tuning writer
@pytest.mark.parametrize("basis", ["auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity"])
@pytest.mark.parametrize("coherence,mode", [("coherent", "a_Jz"), ("incoherent", "rho"), ("none", "none")])
def test_minpom_setup_independent_spin_steering(cdir, basis, coherence, mode):
    _, auxiliary = resonance.setup(model="MP", resonances=["f2_1270"],
                                   production={"coherence": coherence, "basis": basis})
    tunename = "TUNE_TEST_MP_BASIS"
    _create_card(cdir, tunename, auxiliary)
    block = _active_model_block(_load_generated_card(cdir, tunename, "RES/f2_1270.json"), "MP")
    assert block["basis"] == basis
    assert block["polarization"]["mode"] == mode


@pytest.mark.parametrize(
    ("params", "message"),
    [
        (
            {
                "DECAY|225:[321,-321]:zeta.GP": 0.2,
                tools.phase_u_key("DECAY|225:[321,-321]:zeta.GP"): 0.6,
                tools.phase_v_key("DECAY|225:[321,-321]:zeta.GP"): 0.8,
            },
            "Cannot mix raw and encoded phase inputs",
        ),
        (
            {tools.phase_u_key("DECAY|225:[321,-321]:alpha_ls(2,0)"): 0.6},
            "Unsupported encoded phase target",
        ),
        (
            {tools.phase_u_key("RES|eta:XP:g_ls(1,1)"): 0.6},
            "Incomplete encoded phase pair",
        ),
        (
            {
                tools.raw_phase_key("RES|eta:XP:g_ls(1,1)"): 0.3,
                tools.phase_u_key("RES|eta:XP:g_ls(1,1)"): 0.6,
                tools.phase_v_key("RES|eta:XP:g_ls(1,1)"): 0.8,
            },
            "Cannot mix phase encodings",
        ),
    ],
)
def test_mixed_incomplete_phase_params(params, message, cdir):
    with pytest.raises(Exception, match=message):
        _create_card(cdir, "TUNE_TEST_PHASE_ERROR", params)


# Check model phases and Regge settings change only their selected model
@pytest.mark.parametrize("model", ["MP", "XP", "GP", "TP"])
def test_decay_phase_regge_independence(cdir, model):
    tunename = "MODEL_PARAMS"
    parameters = {tools.raw_phase_key(f"DECAY|225:[321,-321]:zeta.{model}"): 0.37}
    if model != "TP":
        parameters.update({
            f"REGGE|omega.{model}": 0.23,
            f"REGGE|DECAY_BARRIERS.{model}": True,
            f"REGGE|use_zeta.{model}": False,
            f"REGGE|FORWARD_VERTEX.{model}": "unit_residue",
            f"REGGE|PHOTON_VERTEX.{model}": "QED",
            f"REGGE|TU_SIGN.{model}": "positive",
        })
    _create_card(cdir, tunename, parameters)
    tune0 = Path(cdir) / "modeldata" / "TUNE0"
    output = Path(cdir) / "modeldata" / tunename
    before = graniitti_driver.json5.loads((tune0 / "DECAYS.json").read_text())
    after = graniitti_driver.json5.loads((output / "DECAYS.json").read_text())
    before["225"]["[321,-321]"]["zeta"][model] = pytest.approx(0.37)
    assert after == before
    before = graniitti_driver.json5.loads((tune0 / "GENERAL.json").read_text())
    after = graniitti_driver.json5.loads((output / "GENERAL.json").read_text())
    for key, value in parameters.items():
        if key.startswith("REGGE|"):
            field = key.partition("|")[2].partition(".")[0]
            before["PARAM_REGGE"][field][model] = value
    assert after == before


# Fit a shared scalar once and preserve its relative file references in the trial tune
def test_con_fit_follows_arbitrary_shared_scalar(monkeypatch, cdir):
    tune_dir = cdir / "modeldata/TUNE0"
    card_path = tune_dir / "CON_XP.json"
    card = load_json_file(card_path, loader=graniitti_driver.json5.load)
    source = tune_dir / "shared.json"
    source.write_text('{"slope": 0.3}', encoding="utf-8")
    for pairs in card.values():
        if "[211,211]" in pairs:
            pairs["[211,211]"]["FF_offshell"] = {
                "type": "exp", "norm": "pole", "b": {"$ref": "shared.json#/slope"}
            }
    card_path.write_text(json.dumps(card), encoding="utf-8")
    monkeypatch.setattr(continuum, "model_dir", lambda: tune_dir)
    domains.load_model_json.cache_clear()
    params, aux = continuum.setup(final_state_pdgs=[[211, -211]], REGGE=False, tune_couplings=False, continuum_models="XP")
    assert len(params) == 1
    key = next(iter(params))
    assert key.endswith(":FF_offshell.b")
    tunename = "TUNE_TEST_JSON_SHARED_SLOPE"
    _create_card(cdir, tunename, {**aux, key: 1.2}, graniitti_driver.GraniittiDriver())
    trial = cdir / "modeldata" / tunename
    assert load_json_file(source)["slope"] == pytest.approx(0.3)
    assert load_json_file(trial / "shared.json")["slope"] == pytest.approx(1.2)
    resolved = load_json_file(trial / "CON_XP.json")
    raw = json.loads((trial / "CON_XP.json").read_text())
    for exchange, pairs in resolved.items():
        if "[211,211]" in pairs:
            assert pairs["[211,211]"]["FF_offshell"]["b"] == pytest.approx(1.2)
            assert raw[exchange]["[211,211]"]["FF_offshell"]["b"] == {"$ref": "shared.json#/slope"}
    domains.load_model_json.cache_clear()


# Tune nested reggeization controls once per shared meson line
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
def test_continuum_reggeize_controls(cdir, model):
    params, aux = continuum.setup(
        final_state_pdgs=[[211, -211]], continuum_models=model,
        REGGE=False, FF_offshell=False, tune_couplings=False,
        tune_reggeize=True, freeze_scale2_range=(0.4, 2.0),
    )
    assert len(params) == 2
    active = next(key for key in params if key.endswith(":reggeize.active"))
    scale = next(key for key in params if key.endswith(":reggeize.freeze_scale2"))
    assert _uniform_bounds(params[scale]) == pytest.approx((0.4, 2.0))
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(params, aux, str(cdir))
    original = _load_tune0_card(f"CON_{model}.json")
    exchange, pair, _ = scale.split("|")[1].split(":")
    assert initial[scale] == pytest.approx(original[exchange][pair]["reggeize"]["freeze_scale2"])
    name = f"TUNE_TEST_REGGEIZE_{model}"
    _create_card(cdir, name, {**aux, active: False, scale: 0.6}, driver)
    updated = load_json_file(cdir / "modeldata" / name / f"CON_{model}.json")
    for pairs in updated.values():
        if pair in pairs:
            assert pairs[pair]["reggeize"] == {"active": False, "freeze_scale2": 0.6}


# Check campaign width bounds follow the shared physical parameter definition
def test_campaign_width_uses_shared_bounds(tmp_path):
    source = Path(campaign_source("tune-gpom-res"))
    shutil.copytree(source.parent, tmp_path, dirs_exist_ok=True)
    defaults = load_json_file(tmp_path / "_defaults.json")
    defaults["bounds"]["resonance"]["PARAM_RES"]["f2_1950"]["width"] = [0.42, 0.58]
    (tmp_path / "_defaults.json").write_text(json.dumps(defaults))
    card = tmp_path / source.name
    module = load_tunesetup(cdir=ROOT, simdriver="GRANIITTI", name=str(card))
    assert _uniform_bounds(module.param_space["RES|f2_1950:GP:width"]) == pytest.approx((0.42, 0.58))


# Read and write independent pole values through the real tuning driver
@pytest.mark.parametrize("selected", [("GP",), ("MP",), ("XP",), ("TP",), ("GP", "MP", "XP", "TP")])
def test_resonance_poles_are_independent(tmp_path, selected):
    driver = graniitti_driver.GraniittiDriver()
    source = _load_tune0_card("RES/f0_500.json")["PARAM_RES"]["MODELS"]
    parameters = {f"RES|f0_500:{model}:{field}": source[model][field] * (1.1 + 0.1 * index)
                  for index, model in enumerate(selected) for field in ("mass", "width")}
    initial = driver.get_initial_param(parameters, {}, str(ROOT))
    target = tmp_path / "poles"
    original = driver.create_steering_card(param_space=parameters, tunename=str(target), cdir=str(ROOT))
    assert original == pytest.approx(initial)
    generated = load_json_file(target / "RES/f0_500.json")["PARAM_RES"]["MODELS"]
    for model, pole in source.items():
        for field in ("mass", "width"):
            key = f"RES|f0_500:{model}:{field}"
            assert generated[model][field] == pytest.approx(parameters.get(key, pole[field]))
    assert driver.get_initial_param(parameters, {}, str(ROOT), tune_default=str(target)) == pytest.approx(parameters)


# Reject unsupported tensor BW steering before trial generation
def test_tensor_bw_steering_rejected():
    driver = graniitti_driver.GraniittiDriver()
    with pytest.raises(ValueError, match="do not use BW"):
        driver.get_initial_param({"RES|f0_500:TP:BW": 1.0}, {}, str(ROOT))


# Activate a zero GP row while retaining compact parity input and the completed norm
def test_gp_con_zero_row(cdir, gp_pion_basis):
    field, rows = gp_pion_basis
    path = cdir / "modeldata/TUNE0/CON_GP.json"
    data = load_json_file(path, loader=graniitti_driver.json5.load)
    rows[1][-2] = 0.0
    data["990"]["[211,211]"]["opposite"][field] = rows
    path.write_text(json.dumps(data))
    options = dict(final_state_pdgs=[[211, -211]], FF_offshell=False, continuum_models="GP",
                   production_geometry="projective")
    active, _ = continuum.setup(**options)
    params, _ = continuum.setup(**options, production_active_only=False)
    prefix = "CON_GP|990:[211,211]/opposite:"
    angle = prefix + domains.angular_row_name(field, rows[1]) + tools.PROJECTIVE_VECTOR_SUFFIX
    assert angle not in active and angle in params
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(params, {}, str(cdir))
    trial = initial | {angle: 0.3}
    _create_card(cdir, "TUNE_TEST_GP_ZERO", trial, driver)
    generated = _load_generated_card(cdir, "TUNE_TEST_GP_ZERO", "CON_GP.json")["990"]["[211,211]"]["opposite"][field]
    assert [row[:3] for row in generated] == [row[:3] for row in rows]
    assert generated[1][-2] > 0.0
    expected = sum(weight * row[-2]**2 for weight, row in zip((2, 2, 1), rows, strict=True))
    assert sum(weight * row[-2]**2 for weight, row in zip((2, 2, 1), generated, strict=True)) == pytest.approx(expected)
    restored = driver.get_initial_param(params, {}, str(cdir), tune_default="TUNE_TEST_GP_ZERO")
    assert restored == pytest.approx(trial)


# Restrict tuned Regge projections while preserving the other physical rows
@pytest.mark.parametrize("projections", [[0, 1], [0, 1, 2]])
def test_gp_con_m_selection(cdir, gp_pion_basis, projections):
    field, defaults = gp_pion_basis
    params, _ = continuum.setup(final_state_pdgs=[[211, -211]], FF_offshell=False, continuum_models="GP",
                                production_geometry="projective", production_active_only=False, production_m=projections)
    prefix = "CON_GP|990:[211,211]/opposite:"
    selected = {key for key in params if key.startswith(prefix)}
    assert len(selected) == len(projections)
    assert all(any(f",{-m})@" in key for m in projections) for key in selected)
    driver = graniitti_driver.GraniittiDriver()
    initial = driver.get_initial_param(params, {}, str(cdir))
    trial = {key: value * 1.1 if key.endswith("@NORM") else 0.2 for key, value in initial.items()}
    _create_card(cdir, "TUNE_TEST_GP_M", trial, driver)
    rows = _load_generated_card(cdir, "TUNE_TEST_GP_M", "CON_GP.json")["990"]["[211,211]"]["opposite"][field]
    for row, original in zip(rows, defaults, strict=True):
        if abs(row[2]) not in projections:
            assert row == pytest.approx(original)
    restored = driver.get_initial_param(params, {}, str(cdir), tune_default="TUNE_TEST_GP_M")
    assert restored == pytest.approx(trial)
