# Test GRANIITTI resonance tuning helpers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import json
import shutil
from copy import deepcopy

import pytest

pytest.importorskip("ray")
from core.io.serialize import load_json_file
from core.tune.drivers.graniitti.tunesetup import domains, resonance
from core.tune.parameters import tools

from submit import CAMPAIGN_DIR

settings = load_json_file(CAMPAIGN_DIR / "tunecards/graniitti/_defaults.json")

pytestmark = pytest.mark.usefixtures("tuning_context")

# Compute numeric bounds from a Ray Tune uniform parameter
def _uniform_bounds(spec):
    return spec.lower, spec.upper


# Check optional production geometry separates the LS norm and direction
def test_xp_ls_norm_angles():
    p = {}
    rows = [[ell, spin, 1.0, 0.0] for ell, spin in
            ((0, 2), (2, 0), (2, 2), (2, 4), (4, 2), (4, 4), (6, 4))]
    assert resonance._add_production_ls_params(
        p, res="f2_1270", model="XP", rows=rows, max_abs=3.0,
        mode="magnitude", geometry="projective")

    norm_key = tools.projective_norm_key("RES|f2_1270:XP:g_ls(0,2)")
    angle_keys = [key for key in p if tools.is_projective_angle_key(key)]
    assert _uniform_bounds(p[norm_key]) == pytest.approx((0.0, 3.0))
    assert len(angle_keys) == 6
    assert all(
        _uniform_bounds(p[key])
        == pytest.approx((-tools.PROJECTIVE_ANGLE_LIMIT, tools.PROJECTIVE_ANGLE_LIMIT))
        for key in angle_keys
    )
    assert not any(key.endswith("@MAG") and ":XP:g_ls(" in key for key in p)
    topology = tools.build_parameter_topology(p)
    assert topology["groups"][0] == {
        "kind": "projective",
        "base": "RES|f2_1270:XP:g_ls",
        "parameters": angle_keys,
    }


# Check exact-zero LS rows remain outside the normalized production direction
def test_projective_ls_geometry_skips_zero_rows():
    p = {}
    handled = resonance._add_production_ls_params(
        p,
        res="f2_test",
        model="XP",
        rows=[[0, 2, 0.5, 0.0], [2, 0, 0.0, 0.0], [2, 2, 0.2, 0.0]],
        max_abs=3.0,
        mode="magnitude",
        geometry="projective",
    )

    assert handled is True
    assert tools.projective_norm_key("RES|f2_test:XP:g_ls(0,2)") in p
    assert "RES|f2_test:XP:g_ls(2,2)@PROJECTIVE" in p
    assert not any("g_ls(2,0)" in key for key in p)


def test_mp_incoherent_setup_density_matrix_only():
    p, aux = resonance.setup(model="MP", production={"coherence": "incoherent"})

    assert "SPIN|MINPOM_MODE" not in aux
    assert not any(key.startswith("SPIN|MP_") for key in aux)
    assert all(
        value == "rho" for key, value in aux.items() if key.endswith(":MP:polarization.mode")
    )
    assert not any(key.endswith(":MP:basis") for key in aux)
    assert any(tools.is_simplex_theta_key(key) for key in p)
    assert not any("@AJZP" in key for key in p)
    assert not any(key.startswith("CON_") for key in p)


def test_mp_coherent_setup_per_res_basis():
    _, aux = resonance.setup(model="MP", production={"coherence": "coherent"})

    assert "SPIN|MINPOM_MODE" not in aux
    assert not any(key.startswith("SPIN|MP_") for key in aux)
    assert all(
        value == "a_Jz" for key, value in aux.items() if key.endswith(":MP:polarization.mode")
    )
    assert not any(key.endswith(":MP:basis") for key in aux)


def test_frame_steering_exposed_only_for_mp():
    _, mp_aux = resonance.setup(
        model="MP",
        resonances=["f0_980"],
        production={"coherence": "coherent"},
        frame="HX",
    )
    _, xp_aux = resonance.setup(
        model="XP",
        resonances=["f0_980"],
        production={"mode": "none"},
    )
    _, gp_aux = resonance.setup(
        model="GP",
        resonances=["f0_980"],
        production={"mode": "none"},
    )

    assert mp_aux["REGGE|MP_FRAME"] == "HX"
    assert "REGGE|MP_FRAME" not in xp_aux
    assert "REGGE|MP_FRAME" not in gp_aux
    with pytest.raises(ValueError, match="frame applies only"):
        resonance.setup(
            model="XP",
            resonances=["f0_980"],
            production={"mode": "none"},
            frame="HX",
        )


# Check explicit production, decay and transfer tuning selects separate fields
def test_omega_form_factor_tuning():
    p, _ = resonance.setup(
        model="XP",
        resonances=["f0_500", "rho_770"],
        vary=("omega", "FF_prod", "FF_decay"),
    )

    assert "REGGE|omega.XP" in p
    assert _uniform_bounds(p["REGGE|omega.XP"]) == pytest.approx((0.0, 1.0))
    assert _uniform_bounds(p["RES|f0_500:XP:FF_prod.Lambda2"]) == pytest.approx(
        settings["bounds"]["resonance"]["MODELS"]["XP"]["FF_prod"]["Lambda2"]
    )
    assert "REGGE|PARAM_RES:FF_decay.Lambda2" not in p
    assert "REGGE|PARAM_RES:FF_transfer.LambdaInv2" not in p
    assert not any(key.startswith("RES|") and "FF_decay.Lambda2" in key for key in p)

    configured = deepcopy(settings)
    configured["bounds"]["decay"] = {"FF_decay": {"TP": {"Lambda2": [0.6, 3.5]}}}
    with domains.context(cdir=domains.source_root(), model_path=domains.model_dir(), settings=configured):
        tp_mass, _ = resonance.setup(
            model="TP", resonances=["f0_500"], production={"mode": "none"}, vary=("FF_prod", "FF_decay")
        )
    tp_transfer, _ = resonance.setup(
        model="TP", resonances=["f0_500"], production={"mode": "none"}, vary=("FF_transfer",)
    )
    transfer_key = "RES|f0_500:TP:FF_transfer.LambdaInv2"
    assert transfer_key not in tp_mass
    assert set(tp_transfer) == {transfer_key}

    xp_decay, _ = resonance.setup(
        model="XP", resonances=["f0_500"], production={"mode": "none"}, vary=("FF_decay",)
    )
    assert "REGGE|PARAM_RES:FF_prod.Lambda2" not in xp_decay


# Check independent model pole parameters use the configured resonance bounds
def test_f0_mass_width_models():
    limits = settings["bounds"]["resonance"]["PARAM_RES"]["f0_500"]

    setup_kwargs_list = [
        {"model": "XP"},
        {"model": "GP"},
        {"model": "TP"},
        {"model": "MP", "production": {"coherence": "coherent"}},
        {"model": "MP", "production": {"coherence": "incoherent"}},
    ]
    for setup_kwargs in setup_kwargs_list:
        p, _ = resonance.setup(**setup_kwargs, resonances=["f0_500"], vary=("mass_width",))

        prefix = f"RES|f0_500:{setup_kwargs['model']}:"
        assert _uniform_bounds(p[prefix + "mass"]) == pytest.approx(limits["mass"])
        assert _uniform_bounds(p[prefix + "width"]) == pytest.approx(limits["width"])
        assert not any(key.endswith((":mass", ":width")) and not key.startswith(prefix) for key in p)

    p, _ = resonance.setup(model="XP", resonances=["f0_980"])

    assert "RES|f0_500:XP:mass" not in p
    assert "RES|f0_500:XP:width" not in p


def test_gp_cannot_override_the_card_basis():
    with pytest.raises(ValueError, match="applies only to model"):
        resonance.setup(
            model="GP",
            resonances=["f2_1270"],
            production={"basis": "helicity", "mode": "magnitude"},
        )


def test_old_tune_mode_names_are_rejected():
    with pytest.raises(ValueError, match=r"Unknown production\.mode"):
        resonance.setup(
            model="XP",
            resonances=["f0_500"],
            production={"mode": "magnitude_phase"},
        )
    with pytest.raises(ValueError, match=r"Unknown production\.mode"):
        resonance.setup(
            model="GP",
            resonances=["f0_500"],
            production={"mode": "phase"},
        )


# Reject misspelled grouped configuration fields
def test_grouped_setup_unknown_fields():
    with pytest.raises(ValueError, match="Unknown production configuration keys"):
        resonance.setup(
            model="XP",
            resonances=["f0_500"],
            production={"coupling_mode": "magnitude"},
        )


# Reject unknown varied parameter families
def test_grouped_setup_unknown_vary_names():
    with pytest.raises(ValueError, match="Unknown varied resonance parameter families"):
        resonance.setup(
            model="XP",
            resonances=["f0_500"],
            vary=("form_factor",),
        )


# Check physical production and decay bounds reach separate model parameters
@pytest.mark.parametrize("model", ["MP", "XP", "GP", "TP"])
def test_separate_production_decay_bounds(model, tmp_path):
    configured = deepcopy(settings)
    bounds = configured["bounds"]["resonance"]
    bounds["MODELS"][model]["FF_prod"]["Lambda2"] = [1.2, 2.4]
    configured["bounds"]["decay"] = {"FF_decay": {model: {"Lambda2": [2.5, 3.5]}}}
    source = tmp_path / "tune"
    shutil.copytree(domains.model_dir(), source)
    decays = deepcopy(resonance._load_branching_data())
    pdg = str(resonance._load_resonance_param("f0_500")["PDG"])
    for channel, block in decays[pdg].items():
        if resonance._is_decay_channel_key(channel):
            block["FF_decay"][model] = {"type": "gaussian", "norm": "pole", "Lambda2": 1.0}
    (source / "DECAYS.json").write_text(json.dumps(decays))
    with domains.context(cdir=domains.source_root(), model_path=source, settings=configured):
        parameters = resonance.res_form_factors(model, resonances=["f0_500"], fields={"FF_prod", "FF_decay"})
    assert _uniform_bounds(parameters[f"RES|f0_500:{model}:FF_prod.Lambda2"]) == pytest.approx((1.2, 2.4))
    decay = [value for key, value in parameters.items() if key.startswith("DECAY|")]
    assert decay
    assert all(_uniform_bounds(value) == pytest.approx((2.5, 3.5)) for value in decay)
