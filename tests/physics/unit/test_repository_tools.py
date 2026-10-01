# Tests for standalone physics card and coupling tools
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import cmath
import copy
import importlib
import json
import math
import re
import shutil
import subprocess
import sys
from dataclasses import replace
from pathlib import Path

import pytest
import sympy as sp
import sympy.physics.wigner as wigner

from develop.tools.lib import common as tool_common
from develop.tools.lib import pole as pole_math
from develop.tools.lib import soft_exchange
from develop.tools.lib.hera import couplings as hera_couplings
from develop.tools.lib.hera import fit as hera_fit
from develop.tools.lib.hera import inputs as hera_inputs
from develop.tools.lib.hera import model as hera_model
from develop.tools.lib.hera import tune as hera_tune
from icepack._common import hepdata
from icepack.PHOTOPROD._common import hera_reader

ROOT = Path(__file__).resolve().parents[3]
TOOLS = ROOT / "develop" / "tools"


# Require supported measurements for the complete HERA derivation
@pytest.fixture(scope="module")
def hera_channels():
    try:
        return hera_couplings.channels()
    except ValueError as error:
        if str(error) in {
            "Measurement unavailable: a supported HEPData JSON input is required",
            "H1 covariance fit unavailable: a supported HEPData JSON covariance input is required",
        }:
            pytest.skip(str(error))
        raise


# Load the actual tool module with fresh module state for each test
def load_tool(name):
    return importlib.reload(importlib.import_module(f"develop.tools.{name}"))


# Run one repository tool in the same interpreter environment as pytest
def run_tool(name, *arguments, check=True):
    command = [sys.executable, str(TOOLS / f"{name}.py"), *arguments]
    return subprocess.run(command, cwd=ROOT, check=check, capture_output=True, text=True)


# Load the checked-in TUNE0 general card
def load_tune0_general():
    json5 = pytest.importorskip("pyjson5")
    path = ROOT / "modeldata" / "TUNE0" / "GENERAL.json"
    with path.open(encoding="utf-8") as stream:
        return json5.load(stream)


# Check a confirmed push cannot overwrite card edits made during the preview
def test_tool_push_edits_made_during_confirmation(tmp_path):
    push = load_tool("lib.push")
    cards = [tmp_path / "one.json", tmp_path / "two.json"]
    original = '{"g": 1.0, "phase": 0.0}\n'
    edited = '{"g": 1.0, "phase": 0.5}\n'
    for card in cards:
        card.write_text(original, encoding="utf-8")
    updates = {card: [push.ScalarUpdate(("g",), "g", 2.0)] for card in cards}

    # Simulate a card edit while the parameter preview is being reviewed
    def confirm(rows):
        cards[1].write_text(edited, encoding="utf-8")
        return True

    with pytest.raises(ValueError, match="Card changed after the push preview"):
        push.push_json5_updates(updates, confirm=confirm)
    assert cards[0].read_text(encoding="utf-8") == original
    assert cards[1].read_text(encoding="utf-8") == edited


# Apply scalar and array edits to the same card even when its input path is relative
@pytest.mark.parametrize("relative", [False, True])
def test_tool_push_scalar_array(tmp_path, monkeypatch, relative):
    push = load_tool("lib.push")
    monkeypatch.chdir(tmp_path)
    path = tmp_path / "GENERAL.json"
    path.write_text('{"value": 1.0, "rows": [[1,2.0]]}\n')
    card = Path(path.name) if relative else path
    updates = {card: [push.ScalarUpdate(("value",), "value", 3.0)]}
    arrays = {card: {("rows",): [[1,4.0], [2,5.0]]}}
    assert push.push_json5_updates(updates, array_updates=arrays, confirm=lambda rows: True)
    assert json.loads(path.read_text()) == {"value": 3.0, "rows": [[1,4.0], [2,5.0]]}
    assert not push.push_json5_updates(updates, array_updates=arrays, confirm=lambda rows: True)


# Check arbitrary channel resolved Good Walker direction validation
def test_soft_loader_validates_resolved_direction():
    soft = load_tool("lib.soft_exchange")
    general = load_tune0_general()
    general["PARAM_SOFT"]["active_model"] = "double"

    scaled = copy.deepcopy(general)
    scaled["PARAM_SOFT"]["MODEL"]["double"]["GW"]["a_c"] = [3.0]
    configuration = soft.load(scaled)
    assert configuration.channels == 2
    assert soft.direction([3.0, 4.0], 3) == pytest.approx([0.6, 0.8])

    single = copy.deepcopy(general)
    single["PARAM_SOFT"]["active_model"] = "single"
    single["PARAM_SOFT"]["MODEL"]["single"]["GW"]["a_c"] = []
    assert soft.load(single).channels == 1

    invalid_length = copy.deepcopy(general)
    invalid_length["PARAM_SOFT"]["MODEL"]["double"]["GW"]["a_c"] = []
    with pytest.raises(ValueError, match="must contain N-1 coefficients"):
        soft.load(invalid_length)

    nonfinite = copy.deepcopy(general)
    nonfinite["PARAM_SOFT"]["MODEL"]["double"]["GW"]["a_c"] = [math.inf]
    with pytest.raises(ValueError, match=r"GW.a_c\[0\] must be finite"):
        soft.load(nonfinite)

    zero = copy.deepcopy(general)
    zero["PARAM_SOFT"]["MODEL"]["double"]["GW"]["a_c"] = [0.0]
    with pytest.raises(ValueError, match="direction must be nonzero"):
        soft.load(zero)


# A Good Walker direction is unchanged by a finite positive overall scale
@pytest.mark.parametrize("scale", [1e-300, 1.0, 1e300])
def test_gw_direction_scale(scale):
    soft = load_tool("lib.soft_exchange")
    assert soft.direction([3 * scale, -4 * scale], 3) == pytest.approx([0.6, -0.8])


# Check the explicit forward excitation Pomeron remains unambiguous with several screening Pomerons
def test_soft_forward_selection():
    soft = load_tool("lib.soft_exchange")
    general = load_tune0_general()
    active_name = general["PARAM_SOFT"]["active_model"]
    model = general["PARAM_SOFT"]["MODEL"][active_name]

    model["EXCHANGE"]["P_aux"] = copy.deepcopy(model["EXCHANGE"]["P"])
    general["PARAM_SOFT"]["EXCHANGE_DEF"]["P_aux"] = copy.deepcopy(
        general["PARAM_SOFT"]["EXCHANGE_DEF"]["P"]
    )
    model["EIKONAL"]["screening_exchanges"] = ["P", "P_aux"]
    configuration = soft.load(general)
    assert configuration.forward_excitation_exchange == "P"

    selected_auxiliary = copy.deepcopy(general)
    selected_auxiliary["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"]["excitation_exchanges"] = [
        "P_aux"
    ]
    assert (
        soft.load(selected_auxiliary).forward_excitation_exchange
        == "P_aux"
    )

    missing = copy.deepcopy(general)
    del missing["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"]["excitation_exchanges"]
    with pytest.raises(ValueError, match="excitation_exchanges must be a nonempty array"):
        soft.load(missing)

    disabled = copy.deepcopy(general)
    disabled["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P_aux"]["on"] = False
    disabled["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"]["excitation_exchanges"] = ["P_aux"]
    with pytest.raises(ValueError, match="selects disabled exchange P_aux"):
        soft.load(disabled)

    wrong_role = copy.deepcopy(general)
    wrong_role["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"]["excitation_exchanges"] = ["R_f2"]
    with pytest.raises(ValueError, match="must select pomeron exchanges"):
        soft.load(wrong_role)


# Check card-local SOFT and central mapping validation against the C++ schema
def test_soft_exchange_input():
    soft = load_tool("lib.soft_exchange")
    general = load_tune0_general()
    active_name = general["PARAM_SOFT"]["active_model"]

    # Apply one invalid mutation to an otherwise complete checked-in card
    def require_invalid(mutator, exception, message):
        invalid = copy.deepcopy(general)
        mutator(invalid)
        with pytest.raises(exception, match=message):
            soft.load(invalid)

    require_invalid(
        lambda card: card["PARAM_REGGE"]["EXCHANGES"][0].update({"extra": True}),
        ValueError,
        "must contain only role, pdg, pole_spin, and soft_exchange",
    )
    require_invalid(
        lambda card: card["PARAM_REGGE"]["EXCHANGES"][1].update(
            {"soft_exchange": card["PARAM_REGGE"]["EXCHANGES"][0]["soft_exchange"]}
        ),
        ValueError,
        "mapped by multiple rows",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["EXCHANGE_DEF"]["P"].update({"crossing": -1}),
        ValueError,
        "tau and crossing must agree",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["EXCHANGE_DEF"]["P"].update({"tau": -1, "crossing": -1}),
        ValueError,
        "pomeron requires positive signature",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["EXCHANGE_DEF"]["R_f2"].update(
            {"trajectory_mode": "pion_loop"}
        ),
        ValueError,
        "pion_loop for a nonpomeron exchange",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].update({"sign": 0}),
        ValueError,
        r"sign must be -1 or 1",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].update({"on": 1}),
        TypeError,
        "on must be boolean",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].pop("on"),
        TypeError,
        "on must be boolean",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].update(
            {"eta_mode": "unknown"}
        ),
        ValueError,
        "eta_mode has an unknown mode",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].update(
            {"alpha": [2.1, 0.2], "eta_mode": "gamma"}
        ),
        ValueError,
        "eta_mode has an unknown mode",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"].update(
            {"extra": copy.deepcopy(card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"])}
        ),
        ValueError,
        "EXCHANGE has unknown exchange extra",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"].pop("O3g"),
        ValueError,
        "EXCHANGE is missing O3g",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"].update(
            {"excitation_exchanges": ["P", "P"]}
        ),
        ValueError,
        "excitation_exchanges must not contain duplicates",
    )
    require_invalid(
        lambda card: card["PARAM_SOFT"]["MODEL"][active_name]["EXCHANGE"]["P"].update(
            {"ff": "O3g"}
        ),
        ValueError,
        "pomeron cannot use a three gluon form factor",
    )


# Check full eikonal, helicity, and unused form factor bank validation
def test_soft_model_blocks():
    soft = load_tool("lib.soft_exchange")
    general = load_tune0_general()
    general["PARAM_SOFT"]["active_model"] = "double"
    active_name = general["PARAM_SOFT"]["active_model"]

    # Apply one invalid mutation to the active model
    def require_invalid(mutator, exception, message):
        invalid = copy.deepcopy(general)
        mutator(invalid["PARAM_SOFT"]["MODEL"][active_name])
        with pytest.raises(exception, match=message):
            soft.load(invalid)

    for field, value, exception, message in (
        ("unitarization", "power", ValueError, "unitarization must be exp or q_exp"),
        ("q", 0.0, ValueError, r"q must be in \(0,1\]"),
        ("q", True, TypeError, "q must be numeric"),
    ):
        require_invalid(
            lambda model, field=field, value=value: model["EIKONAL"].update({field: value}),
            exception,
            message,
        )

    require_invalid(
        lambda model: model["EIKONAL"].update({"unitarization": "exp", "q": 0.5}),
        ValueError,
        "exp unitarization requires q equal to one",
    )
    require_invalid(
        lambda model: model["EIKONAL"].update({"screening_exchanges": "P"}),
        ValueError,
        "screening_exchanges must be a nonempty array",
    )
    require_invalid(
        lambda model: model["EIKONAL"].update({"screening_exchanges": []}),
        ValueError,
        "screening_exchanges must be a nonempty array",
    )
    require_invalid(
        lambda model: model["EIKONAL"].update({"screening_exchanges": ["P", "P"]}),
        ValueError,
        "must not contain duplicates",
    )
    require_invalid(
        lambda model: model["EIKONAL"].update({"screening_exchanges": ["missing"]}),
        ValueError,
        "selects unknown exchange missing",
    )
    wildcard = copy.deepcopy(general)
    wildcard_model = wildcard["PARAM_SOFT"]["MODEL"][active_name]
    wildcard_model["EIKONAL"]["screening_exchanges"] = ["*"]
    wildcard_model["EIKONAL"]["excitation_exchanges"] = ["*"]
    wildcard_configuration = soft.load(wildcard)
    assert wildcard_configuration.forward_excitation_exchange == "P"

    require_invalid(
        lambda model: model["EIKONAL"].update({"screening_exchanges": ["*", "P"]}),
        ValueError,
        "wildcard must be used alone",
    )
    require_invalid(
        lambda model: (
            model["EXCHANGE"]["P"].update({"on": False}),
            model["EIKONAL"].update({"screening_exchanges": ["P"]}),
        ),
        ValueError,
        "selects disabled exchange P",
    )

    require_invalid(
        lambda model: model["EIKONAL"].update({"helicity": "regge_vertex"}),
        TypeError,
        "EIKONAL.helicity must be a boolean",
    )
    require_invalid(
        lambda model: (
            model["EIKONAL"].update({"helicity": True}),
            model["EXCHANGE"]["P"].pop("helicity"),
        ),
        TypeError,
        "EXCHANGE.P.helicity must be an object",
    )
    require_invalid(
        lambda model: (
            model["EIKONAL"].update({"helicity": True}),
            model["EXCHANGE"]["P"]["helicity"].pop("kappa"),
        ),
        TypeError,
        "helicity.kappa must be numeric",
    )
    require_invalid(
        lambda model: (
            model["EIKONAL"].update({"helicity": True}),
            model["EXCHANGE"]["P"]["helicity"].update({"B_kappa": -0.1}),
        ),
        ValueError,
        "helicity.B_kappa must be nonnegative",
    )

    matrix_helicity = copy.deepcopy(general)
    matrix_model = matrix_helicity["PARAM_SOFT"]["MODEL"][active_name]
    matrix_model["EIKONAL"]["helicity"] = True
    matrix_model["EXCHANGE"]["P"]["helicity"] = {
        "kappa": [[0.1, -0.2], [-0.2, 0.3]],
        "B_kappa": [[0.4, 0.1], [0.1, 0.5]],
    }
    assert soft.load(matrix_helicity).channels == 2
    require_invalid(
        lambda model: (
            model["EIKONAL"].update({"helicity": True}),
            model["EXCHANGE"]["P"]["helicity"].update({"kappa": [[0.1, 0.2], [0.3, 0.4]]}),
        ),
        ValueError,
        "helicity.kappa must be symmetric",
    )
    require_invalid(
        lambda model: (
            model["EIKONAL"].update({"helicity": True}),
            model["EXCHANGE"]["P"]["helicity"].update({"kappa": [[0.1]]}),
        ),
        ValueError,
        "helicity.kappa must use the common N channel count",
    )

    # Add an unused bank so its invalid rows cannot be reached through an exchange
    def add_unused_kernel(model):
        model["FF"]["UNUSED"] = {
            "type": "GKERNEL",
            "param": [[1.0, 1.0, 0.0, 0.0], [1.0, 1.0, 0.0, 0.0]],
        }

    for coordinate, value in enumerate((-1.0, 0.0, -1.0, -1.0)):

        def mutate_kernel(model, coordinate=coordinate, value=value):
            add_unused_kernel(model)
            model["FF"]["UNUSED"]["param"][0][coordinate] = value

        require_invalid(
            mutate_kernel,
            ValueError,
            "GKERNEL parameters are outside their physical ranges",
        )
    require_invalid(
        lambda model: model["FF"].update(
            {
                "UNUSED": {
                    "type": "GKERNEL",
                    "param": [
                        [1.0, 0.5, 0.0, 0.0],
                        [1.0, 1.0, 0.0, 0.0],
                    ],
                }
            }
        ),
        ValueError,
        "GKERNEL has a nonfinite forward slope",
    )
    require_invalid(
        lambda model: model["FF"].update(
            {
                "UNUSED": {
                    "type": "GKERNEL",
                    "param": [
                        [-0.1, 1.0, 1.0, 0.0, 0.0],
                        [0.0, 1.0, 1.0, 0.0, 0.0],
                    ],
                }
            }
        ),
        ValueError,
        "GKERNEL prefactor must be nonnegative",
    )
    require_invalid(
        lambda model: model["FF"].update({"UNUSED": {"type": "UNKNOWN", "param": [[1.0], [1.0]]}}),
        ValueError,
        "unsupported proton form factor",
    )
    require_invalid(
        lambda model: model["FF"].update({"UNUSED": {"type": "3G", "param": [[1.0]]}}),
        ValueError,
        "param must have N rows",
    )


# Select the unique checked-in resonance card matching one HERA production topology
def tune0_resonance_card_for_channel(channel):
    json5 = pytest.importorskip("pyjson5")
    matches = []
    for path in sorted((ROOT / "modeldata" / "TUNE0" / "RES").glob("*.json")):
        with path.open(encoding="utf-8") as stream:
            card = json5.load(stream)
        param_res = card["PARAM_RES"]
        if int(param_res["PDG"]) != channel.pdg:
            continue
        models = param_res["MODELS"]
        if all(
            f"[{first_pdg},{second_pdg}]" in models[model_name]
            for model_name, first_pdg, second_pdg in channel.production_channels
        ):
            matches.append(card)
    assert len(matches) == 1
    return matches[0]


# Check compact scalar-vector serialization
def test_output_roundtrip_compact_arrays():
    common = load_tool("lib.common")
    payload = {"rows": [[113, 40.0, 9.59], [333, 70.0, 7.3]], "enabled": True}
    encoded = common.dumps(payload)

    assert json.loads(encoded) == payload
    assert "[113,40.0,9.59]" in encoded


# Check the raw vector STF norms for all coupled spins
def test_dl_vector_stf_norms():
    # Scalar contraction delta_ij delta_ij = 3, antisymmetric contraction = 2,
    # and the normalized symmetric traceless tensor has unit norm
    assert [pole_math.stf_norm(1, 1, s) for s in range(3)] == pytest.approx(
        [math.sqrt(3), math.sqrt(2), 1.0]
    )


# Check scalar, spinor and vector pole normalization across exchanges and models
@pytest.mark.parametrize("pdg", [211, 2212, 113])
@pytest.mark.parametrize("spin", range(4))
def test_dl_pole_normalization(pdg, spin):
    dl = load_tool("DL_couplings")
    particles = pole_math.particles()
    first, second = particles[pdg], particles[pdg if pdg == 113 else -pdg]
    mother = pole_math.Particle(990, "test", 2 * spin, (-1)**spin, (-1)**spin)
    allowed = sorted(pole_math.ls_rows(mother, first, second, True, True, crossed=True))
    pole = mother, first, second, allowed
    raw = pole_math.ls_norm(mother, first, second, *allowed[0])
    rows = pole_math.helicity_rows(mother, first, second, allowed, allowed[0], p_symmetry=False)
    assert sum(row[2]**2 for row in rows) == pytest.approx(raw**2, rel=1e-11)
    for model in ("MP", "XP", "GP"):
        # A unit imaginary LS seed tests phase preservation as well as magnitude
        block = {"CP": [True, True], "basis": "crossed_ls", "g_ls": [
            [ell, s] + ([0] if model == "GP" else []) + [float(i == 0), math.pi / 2]
            for i, (ell, s) in enumerate(allowed)
        ]}
        for shape in ("preserve", "lowest"):
            out = copy.deepcopy(block)
            for update in dl._coupling_updates((), block, 1.0, pole, model, shape):
                dl._set(out, update.path, update.value)
            actual = sum(pole_math.ls_norm(mother, first, second, *row[:2])**2 * row[-2]**2
                         for row in out["g_ls"])
            # The unit vertex normalization requires sum |H|^2 / (2s+1) = 1
            assert actual == pytest.approx(first.spin2 + 1, rel=2e-8)
            phase = 1j if shape == "preserve" else 1.0
            assert cmath.rect(out["g_ls"][0][-2], out["g_ls"][0][-1]) == pytest.approx(
                phase * math.sqrt(first.spin2 + 1) / raw, rel=2e-8
            )


# Check JSON5 parsing and exact scalar spans used for card updates
def test_dl_card_document():
    dl = load_tool("DL_couplings")
    document = dl.CardDocument("{/*retain*/ value:[1, 2.,], name:'a//b', flag:true,}")
    assert document.data == {"value": [1, 2.0], "name": "a//b", "flag": True}
    assert document.source[slice(*document.spans[("value", 1)])] == "2."


# Check DL coupling factorization and all model specific Regge aliases
def test_dl_exchange_aliases():
    dl = load_tool("DL_couplings")
    coefficients = dl.coefficients()
    residues = dl.reference_couplings(coefficients)
    couplings = dl.derived_couplings(residues)
    neutral_kaon = next(state for state in dl.FINAL_STATES if state.key == "K0")

    assert residues[("pi", "P")] == pytest.approx(4.688992782, rel=1e-9)
    assert neutral_kaon.pdg == (311, -311)
    assert couplings[("n", "P")][0] == pytest.approx(couplings[("p", "P")][0])
    assert couplings[("pi0", "P")][0] == pytest.approx(couplings[("pi", "P")][0])
    assert couplings[("K0", "P")][0] == pytest.approx(couplings[("K", "P")][0])
    for state in ("pi0", "rho", "phi"):
        assert couplings[(state, "C-")][0] == pytest.approx(0.0, abs=0.0)
    assert dl.model_exchange_pdgs("MP", dl.EXCHANGES[0]) == (991, 993, 995)
    assert dl.model_exchange_pdgs("XP", dl.EXCHANGES[0]) == (991, 993, 995)
    assert dl.model_exchange_pdgs("GP", dl.EXCHANGES[0]) == (990,)
    assert set(dl.model_card_block("MP", couplings, 6)["g"]) == {
        "991",
        "993",
        "995",
        "9915",
        "9925",
        "9933",
        "9943",
    }
    assert set(dl.model_card_block("XP", couplings, 6)["g"]) == {
        "991",
        "993",
        "995",
        "9915",
        "9925",
        "9933",
        "9943",
    }
    assert set(dl.model_card_block("GP", couplings, 6)["g"]) == {
        "990",
        "9910",
        "9920",
        "9930",
        "9940",
    }
    for model, exchange in (("MP", "9933"), ("XP", "9933"), ("GP", "9930")):
        rows = dl.model_card_block(model, couplings, 6)["g"][exchange]
        assert any(row[:2] == [311, -311] for row in rows)


# Check every DL exchange uses its mapped bare scale and the DL hadron ratios
def test_dl_couplings_soft_bare_scales_dl_ratios():
    dl = load_tool("DL_couplings")
    coefficients = dl.coefficients()
    beam = {"P": 8.0, "C+": 16.0, "C-": 4.0}
    optical = {"P": 0.98, "C+": 0.75, "C-": 0.64}
    residues = dl.bare_couplings(coefficients, beam)
    born = dl.born_coefficients(residues, optical)

    for exchange in dl.EXCHANGES:
        key = exchange.key
        assert residues[("p", key)] == pytest.approx(beam[key])
        for hadron in coefficients:
            assert residues[(hadron, key)] / residues[("p", key)] == pytest.approx(
                coefficients[hadron][key] / coefficients["p"][key]
            )
            assert born[hadron][key] == pytest.approx(
                optical[key] * residues[(hadron, key)] * residues[("p", key)] / dl.MB_TO_GEV2
            )
    assert coefficients["K"] == pytest.approx({"P": 11.93, "C+": 16.455, "C-": 8.875})

    payload = json.loads(run_tool("DL_couplings", "--format", "json").stdout)
    configured_beam = payload["residue_factorization"]["beam_residues_GeV_minus1"]
    configured_optical = payload["residue_factorization"]["optical_weights"]
    configured_born = payload["residue_factorization"]["bare_Born_coefficients_mb"]
    for row in payload["residue_couplings_GeV_minus1"]:
        key = row["exchange"]
        hadron = row["hadron"]
        assert row["bare_value"] == pytest.approx(
            configured_beam[key] * coefficients[hadron][key] / coefficients["p"][key]
        )
        assert configured_born[hadron][key] == pytest.approx(
            configured_optical[key] * configured_beam[key] * row["bare_value"] / dl.MB_TO_GEV2
        )

    with pytest.raises(ValueError, match="finite and positive"):
        dl.bare_couplings(coefficients, {"P": 0.0, "C+": 1.0, "C-": 1.0})
    with pytest.raises(ValueError, match="optical weight must be finite and positive"):
        dl.born_coefficients(
            residues,
            {"P": 1.0, "C+": 0.0, "C-": 1.0},
        )


# Check neutral kaon amplitudes reproduce the neutron total cross sections by isospin
@pytest.mark.parametrize("sqrt_s", [10.0, 100.0, 1000.0])
def test_dl_neutral_kaon_optical_theorem(sqrt_s):
    dl = load_tool("DL_couplings")
    residues = dl.reference_couplings(dl.coefficients())
    couplings = dl.derived_couplings(residues)
    s = sqrt_s**2
    even = odd = 0j
    for exchange in dl.EXCHANGES:
        key = exchange.key
        alpha = 1 + dl.DL_EPSILON if key == "P" else 1 - dl.DL_ETA
        tau = -1 if key == "C-" else 1
        amplitude = dl.raw_eta_factor(alpha, tau) * s**alpha * residues[("p", key)] * couplings[("K0", key)][0]
        if tau > 0:
            even += amplitude
        else:
            odd += amplitude
    # [REFERENCE: LNS, arXiv:1804.04706, Eq. (3.29)]
    for crossing, neutron_y in ((1, 9.08), (-1, 19.09)):
        sigma = (even + crossing * odd).imag / (s * dl.MB_TO_GEV2)
        expected = dl.DL_MB["K"]["P"] * s**dl.DL_EPSILON + neutron_y * s**(-dl.DL_ETA)
        assert sigma == pytest.approx(expected)
    assert couplings[("K0", "P")][0] == pytest.approx(couplings[("K", "P")][0])


# Check a confirmed DL push changes only matched continuum magnitude tokens
def test_dl_push_updates_direct_con_couplings(tmp_path):
    json5 = pytest.importorskip("pyjson5")
    dl = load_tool("DL_couplings")
    tune_dir = tmp_path / "TUNE0"
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune_dir)
    general_path = tune_dir / "GENERAL.json"
    con_path = tune_dir / "CON_GP.json"
    with con_path.open(encoding="utf-8") as stream:
        source_card = json5.load(stream)
    pair = source_card["990"]["[211,211]"]
    for sector in ("same", "opposite"):
        pair[sector] = {"basis": "crossed_ls", "CP": [True, True],
                        "g_ls": [[2, 0, -2, 1.2, 0.2], [2, 0, -1, 0.8, -0.2], [2, 0, 0, 1.0, 0.0]]}
    pair["FF_transfer"] = {"type": "none"}
    pair["FF_offshell"] = {"type": "power", "norm": "pole", "Lambda2": 1.3, "n": 1.5}
    con_path.write_text(json.dumps(source_card), encoding="utf-8")
    source = con_path.read_text(encoding="utf-8")
    source_ff = source_card["990"]["[211,211]"]["FF_transfer"]
    source_offshell = source_card["990"]["[211,211]"]["FF_offshell"]
    generated_g = {"990": [[211, -211, 5.333, 0.0], [None, None, 4.5, 0.0]]}
    previews = []

    applied = dl.push_cards(
        general_path,
        {"GP": {"g": generated_g}},
        confirm=lambda rows: previews.extend(rows) or True,
        shape="preserve",
    )

    assert applied
    assert con_path.read_text(encoding="utf-8").startswith(source.split("{", 1)[0])
    with con_path.open(encoding="utf-8") as stream:
        card = json5.load(stream)
    for sector in ("same", "opposite"):
        block = card["990"]["[211,211]"][sector]
        rows = block["g_ls"]
        original = source_card["990"]["[211,211]"][sector]["g_ls"]
        anchor = next(row for row in rows if row[2] == 0)
        original_anchor = next(row for row in original if row[2] == 0)
        assert anchor[3] == pytest.approx(5.333 / math.sqrt(2.0 / 3.0))
        scale = anchor[3] / original_anchor[3]
        for before, after in zip(original, rows, strict=True):
            assert after[:3] == before[:3]
            assert float(after[3]) == pytest.approx(scale * float(before[3]))
            assert float(after[4]) == pytest.approx(float(before[4]))
    assert card["990"]["[211,211]"]["FF_transfer"] == source_ff
    assert card["990"]["[211,211]"]["FF_offshell"] == source_offshell
    assert previews


# Check a GP DL push retains the absolute canonical pole coefficient
def test_dl_push_scales_gp_helicity_shape():
    dl = load_tool("DL_couplings")
    card = {
        "990": {
            "[2212,2212]": {
                "opposite": {
                    "basis": "crossed_helicity",
                    "CP": [True, True],
                    "helicity": [
                        [-0.5, -0.5, 0, 3.0, 0.2],
                        [-0.5, 0.5, 1, 4.0, -0.3],
                    ],
                }
            }
        }
    }
    generated = {
        "g": {
            "990": [[2212, -2212, 10.0, 0.0], [None, None, 4.0, 0.0]],
        }
    }

    updates = dl.card_updates(card, generated, "GP", shape="preserve")
    values = {update.path: update.value for update in updates}
    base = ("990", "[2212,2212]", "opposite", "helicity")
    scale = 10.0 * math.sqrt(2.0 / 18.0)
    assert values[(*base, 0, 3)] == pytest.approx(3.0 * scale)
    assert values[(*base, 1, 3)] == pytest.approx(4.0 * scale)


# Check finite m rows scale through the equal spin-one m zero pole anchor
def test_dl_push_scales_equal_spin_gp_shape():
    dl = load_tool("DL_couplings")
    card = {
        "990": {
            "[113,113]": {
                "self": {
                    "basis": "crossed_helicity",
                    "CP": [True, True],
                    "helicity": [[-1, 0, 1, 3.0, 0.2], [-1, 0, 0, 2.0, 0.1]],
                }
            }
        }
    }
    generated = {
        "g": {
            "990": [[113, -113, 10.0, 0.0], [None, None, 4.0, 0.0]],
        }
    }

    updates = dl.card_updates(card, generated, "GP", shape="preserve")
    values = {update.path: update.value for update in updates}
    base = ("990", "[113,113]", "self", "helicity")
    assert values[(*base, 0, 3)] == pytest.approx(7.5 * math.sqrt(3.0))
    assert values[(*base, 1, 3)] == pytest.approx(5.0 * math.sqrt(3.0))


# Check a finite m GP shape without a reduced pole anchor is rejected
def test_dl_gp_missing_anchor():
    dl = load_tool("DL_couplings")
    card = {
        "990": {
            "[113,113]": {
                "self": {
                    "basis": "crossed_helicity",
                    "CP": [True, True],
                    "helicity": [[-1, 0, 1, 3.0, 0.2]],
                }
            }
        }
    }
    generated = {
        "g": {
            "990": [[113, -113, 10.0, 0.0], [None, None, 4.0, 0.0]],
        }
    }

    with pytest.raises(ValueError, match="requires an m=0 pole anchor"):
        dl.card_updates(card, generated, "GP", shape="preserve")


# Check a crossing-related DL coupling seeds each allowed sector at its lowest LS row
def test_dl_push_fills_lowest_ls_row():
    dl = load_tool("DL_couplings")
    card = {"995": {"[2212,2212]": {
        "same": {"basis": "crossed_ls", "CP": [True, True],
                 "g_ls": [[2, 1, 3.0, 0.4], [2, 0, 2.0, 0.3]]},
        "opposite": {"basis": "crossed_ls", "CP": [True, True],
                     "g_ls": [[3, 1, 3.0, 0.4], [1, 1, 2.0, 0.3]]},
    }}}
    generated = {"g": {"995": [[2212, -2212, 8.0, 0.0]]}}
    updates = dl.card_updates(card, generated, "XP", shape="lowest")
    values = {update.path: update.value for update in updates}
    for sector in ("same", "opposite"):
        base = ("995", "[2212,2212]", sector, "g_ls")
        assert values[(*base, 0, 2)] == 0.0
        raw_norm2 = 4.0 / 3.0 if sector == "same" else 10.0 / 3.0
        assert values[(*base, 1, 2)] == pytest.approx(8.0 * math.sqrt(2.0 / raw_norm2))
        assert values[(*base, 0, 3)] == 0.0
        assert values[(*base, 1, 3)] == 0.0
    with pytest.raises(ValueError, match="no supported DL target channels"):
        dl.card_updates({"995": {}}, generated, "XP", shape="lowest")


# Check the canonical spinor vertex in both bases and all three models
@pytest.mark.parametrize("model", ["MP", "XP", "GP"])
@pytest.mark.parametrize("field", ["g_ls", "helicity"])
@pytest.mark.parametrize("shape", ["preserve", "lowest"])
def test_dl_push_spinor_pole_in_both_bases(model, field, shape):
    dl = load_tool("DL_couplings")
    exchange = "990" if model == "GP" else "995"
    # Raw norms are N(1,1)^2 = 10/3 and N(3,1)^2 = 4/5
    rows = [[3, 1, 2.0, 0.2], [1, 1, 3.0, -0.4]] if field == "g_ls" else [
        [-0.5, -0.5, 2.0, 0.2], [-0.5, 0.5, 3.0, -0.4]
    ]
    if model == "GP":
        rows = [[*row[:2], 0, *row[2:]] for row in rows]
        rows.append([*rows[0][:2], -1, 5.0, 0.7])
    block = {"basis": "crossed_ls" if field == "g_ls" else "crossed_helicity", "CP": [True, True], field: rows}
    card = {exchange: {"[2212,2212]": {"opposite": block}}}
    generated = {"g": {exchange: [[2212, -2212, 7.0, 0.0]]}}
    updated = copy.deepcopy(card)
    for update in dl.card_updates(card, generated, model, shape=shape):
        target = updated
        for part in update.path[:-1]:
            target = target[part]
        target[update.path[-1]] = update.value
    actual = updated[exchange]["[2212,2212]"]["opposite"][field]
    if shape == "lowest":
        coefficient = 7.0 * math.sqrt(3.0 / 5.0)
        expected = [0.0, coefficient] if field == "g_ls" else [coefficient * math.sqrt(2.0 / 3.0), coefficient]
        assert [row[-2] for row in actual[:2]] == pytest.approx(expected)
        assert all(abs(row[-1]) < 1e-12 for row in actual)
        if model == "GP":
            assert actual[2][-2] == 0.0
    else:
        norm2 = (4.0 / 5.0) * 2.0**2 + (10.0 / 3.0) * 3.0**2 if field == "g_ls" else 2 * (2.0**2 + 3.0**2)
        scale = 7.0 * math.sqrt(2.0 / norm2)
        assert [row[-2] for row in actual] == pytest.approx([row[-2] * scale for row in rows])
        assert [row[-1] for row in actual] == [row[-1] for row in rows]


# Check DL updates reject an ambiguous self-conjugate family declaration
def test_dl_push_rejects_self_charge_sectors():
    dl = load_tool("DL_couplings")
    sector = {
        "basis": "crossed_helicity",
        "CP": [True, True],
        "helicity": [[0, 0, 0, 1.0, 0.0]],
    }
    card = {
        "990": {
            "[111,111]": {
                "self": copy.deepcopy(sector),
                "same": copy.deepcopy(sector),
            }
        }
    }
    generated = {"g": {"990": [[111, 111, 4.0, 0.0]]}}

    with pytest.raises(ValueError, match="self sector cannot coexist"):
        dl.card_updates(card, generated, "GP", shape="preserve")


# Check declining the DL CLI preview leaves the tune card byte identical
@pytest.mark.parametrize("push_args", [(), ("--push",)])
@pytest.mark.parametrize("reply", ["cancel\n", "", "preserve\nno\n", "lowest\nno\n"])
def test_dl_push_decline_general_card(tmp_path, push_args, reply):
    tune_dir = tmp_path / "TUNE0"
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune_dir)
    general_path = tune_dir / "GENERAL.json"
    before = {path: path.read_bytes() for path in tune_dir.glob("*.json")}
    command = [
        sys.executable,
        str(TOOLS / "DL_couplings.py"),
        *push_args,
        "--tune-general",
        str(general_path),
    ]

    subprocess.run(
        command,
        cwd=ROOT,
        input=reply,
        check=True,
        capture_output=True,
        text=True,
    )

    assert {path: path.read_bytes() for path in tune_dir.glob("*.json")} == before


# Check one CLI confirmation updates real cards with the stated pole strength
@pytest.mark.parametrize(("name", "shape"), [("DL_couplings", "preserve"), ("DL_couplings", "lowest"), ("HERA_couplings", None)])
def test_coupling_push_residues(tmp_path, name, shape, request):
    if name == "HERA_couplings":
        request.getfixturevalue("hera_channels")
    json5 = pytest.importorskip("pyjson5")
    tool = load_tool(name)
    tune_dir = tmp_path / "TUNE0"
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune_dir)
    general_path = tune_dir / "GENERAL.json"
    general = json5.loads(general_path.read_text(encoding="utf-8"))
    soft = general["PARAM_SOFT"]
    matrix = soft["MODEL"][soft["active_model"]]["EXCHANGE"]["P"]["g"]
    for row in matrix:
        row[:] = [1.25 * value for value in row]
    general_path.write_text(json.dumps(general), encoding="utf-8")
    command = [sys.executable, str(TOOLS / f"{name}.py"), "--tune-general", str(general_path)]
    if shape == "lowest":
        command.extend(("--shape", shape))
    reply = "yes\n" if name == "DL_couplings" else "invalid\nyes\n"
    subprocess.run(command, cwd=ROOT, input=reply, check=True, capture_output=True, text=True)

    if name == "DL_couplings":
        beam = tool.beam_couplings(tool.load_config(general_path))
        for model in ("MP", "XP", "GP"):
            card = json5.loads((tune_dir / f"CON_{model}.json").read_text(encoding="utf-8"))
            for exchange in tool.EXCHANGES:
                for pdg in tool.model_exchange_pdgs(model, exchange):
                    for sector in ("same", "opposite"):
                        block = card[str(pdg)]["[2212,2212]"][sector]
                        particles = tool.load_particles(tune_dir=tune_dir)
                        pole = tool._pole_reference(pdg, [2212, 2212], sector, block, particles)
                        physical = tool.pole_tensor_norm(block, pole, model) / math.sqrt(2.0)
                        assert physical == pytest.approx(beam[exchange.key], rel=2e-6)
    else:
        proton = hera_model.load_proton(general_path)
        for channel in hera_couplings.channels(proton):
            _, resonance = hera_tune._resonance_card(tune_dir / "RES", channel)
            for model, first, second in channel.production_channels:
                block = resonance["PARAM_RES"]["MODELS"][model][f"[{first},{second}]"]
                coupling = block["helicity"][0][2]
                forward = coupling**2 * proton.beam_residue_per_gev**2 * hera_model.gev2_to_microbarn() / (16.0 * math.pi)
                assert forward == pytest.approx(
                    channel.forward_dsigma_dt_ub_per_gev2.central * channel.normalization_scale, rel=5e-6
                )

    after = {path: path.read_bytes() for path in tune_dir.rglob("*.json")}
    reply = "yes\n"
    subprocess.run(command, cwd=ROOT, input=reply, check=True, capture_output=True, text=True)
    assert {path: path.read_bytes() for path in tune_dir.rglob("*.json")} == after


# Check DL card generation inherits a changed fitted SOFT beam coupling
def test_dl_push_payload_tracks_soft_beam_residue(tmp_path):
    json5 = pytest.importorskip("pyjson5")
    dl = load_tool("DL_couplings")
    for name in ("first", "second"):
        shutil.copytree(ROOT / "modeldata" / "TUNE0", tmp_path / name)
    first_path = tmp_path / "first/GENERAL.json"
    second_path = tmp_path / "second/GENERAL.json"
    source = (ROOT / "modeldata" / "TUNE0" / "GENERAL.json").read_text(encoding="utf-8")
    first_path.write_text(source, encoding="utf-8")
    mutated = json5.loads(source)
    active = mutated["PARAM_SOFT"]["active_model"]
    coupling = mutated["PARAM_SOFT"]["MODEL"][active]["EXCHANGE"]["P"]["g"]
    coupling[0][0] *= 1.5
    second_path.write_text(json.dumps(mutated), encoding="utf-8")

    first = run_tool(
        "DL_couplings",
        "--format",
        "cards",
        "--model",
        "GP",
        "--tune-general",
        str(first_path),
    ).stdout
    second = run_tool(
        "DL_couplings",
        "--format",
        "cards",
        "--model",
        "GP",
        "--tune-general",
        str(second_path),
    ).stdout

    first_card = json.loads(first)
    second_card = json.loads(second)
    first_pion = first_card["CON"]["GP"]["990"]["[211,211]"]["opposite"]["helicity"][0][-2]
    second_pion = second_card["CON"]["GP"]["990"]["[211,211]"]["opposite"]["helicity"][0][-2]
    assert second_pion > first_pion
    assert second_pion / first_pion == pytest.approx(
        dl.beam_couplings(dl.load_config(second_path))["P"]
        / dl.beam_couplings(dl.load_config(first_path))["P"],
        rel=2.0e-6,
    )


# Check Pomeron and secondary Reggeon phase calculations and payload structure
def test_dl_exchange_phases():
    dl = load_tool("DL_couplings")
    assert dl.raw_eta_factor(1.0, 1) == pytest.approx(1j)
    assert dl.raw_eta_factor(0.0, -1) == pytest.approx(-1j)
    assert dl.raw_eta_factor(0.5, 1) == pytest.approx(complex(-1.0, 1.0))
    assert not math.isfinite(dl.raw_eta_factor(0.0, 1).real)
    assert dl.eta_factor(0.0, 1.08, 1, "raw") == 0j
    with pytest.raises(ValueError, match="unknown mode full_t0"):
        dl.eta_factor(0.7, 1.08, 1, "full_t0")
    assert dl.eta_factor(0.7, 1.08, 1, "rotating_t0") == pytest.approx(
        dl.rotating_eta_factor(1.08, 1)
    )
    with pytest.raises(ValueError, match="unknown mode gamma"):
        dl.eta_factor(0.0, 1.08, 1, "gamma")
    with pytest.raises(ValueError, match="exact integer"):
        soft_exchange.exact_integer(0.5, "PARAM_SOFT.EXCHANGE_DEF.P.tau")
    with pytest.raises(TypeError, match="must be numeric"):
        soft_exchange.finite_number(True, "PARAM_SOFT.EXCHANGE.P.alpha[0]")
    for alpha in (-0.7, 0.3, 1.17, 1.6):
        for tau in (-1, 1):
            expected = -(1 + tau * cmath.exp(-1j * math.pi * alpha)) / math.sin(math.pi * alpha)
            assert dl.raw_eta_factor(alpha, tau) == pytest.approx(expected)
            assert abs(dl.rotating_eta_factor(alpha, tau)) == pytest.approx(1.0)

    configuration = dl.load_config()
    optical = dl.optical_factors(configuration)
    configured = {row["exchange"]: row for row in dl.exchange_output(configuration)}
    assert set(configured) == {"P", "C+", "C-", "O"}
    assert configured["P"]["eta_mode"] == "rotating"
    assert configured["P"]["configured_eta_factor_t0"] == pytest.approx(
        dl.complex_payload(dl.rotating_eta_factor(configured["P"]["a0"], 1))
    )
    assert configured["C+"]["B_GeV_minus2"] > 0.0
    general = load_tune0_general()
    model_name = general["PARAM_SOFT"]["active_model"]
    model = general["PARAM_SOFT"]["MODEL"][model_name]
    for row, exchange in zip(
        general["PARAM_REGGE"]["EXCHANGES"],
        ("P", "C+", "C-", "O"),
        strict=True,
    ):
        residue = soft_exchange.proton_vertex(model, row["soft_exchange"], model_name)
        assert configured[exchange]["beam_residue_t0_GeV_minus1"] == pytest.approx(
            residue.beam_residue_per_gev
        )
        assert configured[exchange]["beam_residue_source"] == residue.source
    assert configured["C+"]["eta_mode"] == "rotating"
    assert configured["C-"]["eta_mode"] == "rotating"
    assert configured["O"]["eta_mode"] == "rotating"
    for exchange in dl.EXCHANGES:
        assert optical[exchange.key] == pytest.approx(
            abs(
                dl.eta_factor(
                    configured[exchange.key]["a0"],
                    configured[exchange.key]["a0"],
                    configured[exchange.key]["tau"],
                    configured[exchange.key]["eta_mode"],
                ).imag
            )
        )

    payload = json.loads(run_tool("DL_couplings", "--format", "json").stdout)
    assert len(payload["configured_exchange_parameters"]) == 4
    assert payload["eta_convention"]["configured_eta_mode"] == [
        row.eta_mode for row in configuration["exchanges"]
    ]
    modes = payload["eta_convention"]["soft_modes"]
    assert modes["model"] == configuration["soft_model"]
    assert modes["eta_mode_P"] == configuration["soft_eta_mode_P"]
    assert modes["eta_mode_O"] == configuration["soft_eta_mode_O"]
    assert modes["3P"]["eta_mode"] == configuration["soft_eta_mode_3P"]
    assert payload["eta_convention"]["photoprod_eta_mode"] == configuration["photoprod_eta_mode"]


# Check the amplitude-level HERA fit point and J/psi steering rows
@pytest.mark.usefixtures("hera_channels")
def test_hera_jpsi_coupling_reconstructs_forward_xs():
    channel = next(row for row in hera_couplings.channels() if row.pdg == 443)
    residue = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    derived = hera_model.derive(channel, residue)
    payload = hera_couplings.channel_output(channel, derived)

    general = load_tune0_general()
    soft = general["PARAM_SOFT"]
    model = soft["MODEL"][soft["active_model"]]
    pomeron = model["EXCHANGE"]["P"]
    expected_residue = soft_exchange.coupling(model, "P")
    bank = pomeron["ff"]
    form_factor = str(model["FF"][bank]["type"])
    expected_parameters = tuple(
        tuple(float(value) for value in row) for row in model["FF"][bank]["param"]
    )
    assert residue.beam_residue_per_gev == pytest.approx(expected_residue)
    assert residue.form_factor == form_factor
    assert residue.parameters == expected_parameters
    proton = soft_exchange.proton_state(model["GW"]["theta"], len(pomeron["g"]))
    slopes = [
        soft_exchange.form_slope(form_factor, tuple(float(value) for value in row))
        for row in model["FF"][bank]["param"]
    ]
    expected_derivative = sum(
        proton[row]
        * float(pomeron["g"][row][column])
        * proton[column]
        * 0.5
        * (slopes[row] + slopes[column])
        for row in range(len(pomeron["g"]))
        for column in range(len(pomeron["g"]))
        if pomeron["transition_ff"] != "diagonal" or row == column
    )
    assert residue.amplitude_slope_per_gev2 == pytest.approx(expected_derivative / expected_residue)
    model["EXCHANGE"]["P"]["g"] = [
        [2.0 * value for value in row] for row in model["EXCHANGE"]["P"]["g"]
    ]
    assert soft_exchange.coupling(model, "P") == pytest.approx(2.0 * expected_residue)

    fitted_forward = channel.forward_dsigma_dt_ub_per_gev2.central * channel.normalization_scale
    assert derived.gamma_pomeron_vector_coupling_per_gev == pytest.approx(
        hera_model.coupling(fitted_forward, residue.beam_residue_per_gev)
    )
    assert derived.reconstructed_forward_dsigma_dt_ub_per_gev2 == pytest.approx(fitted_forward)
    assert derived.exponential_sigma_ub == pytest.approx(
        hera_model.sigma_exp(
            channel.forward_dsigma_dt_ub_per_gev2.central,
            channel.b0_per_gev2.central,
            channel.tmax_gev2,
        )
    )
    assert derived.matched_sigma_ub == pytest.approx(
        fitted_forward * hera_model.profile_moments(residue, channel.vector_cross_section_slope_per_gev2,
                                              channel.tmax_gev2)[0]
    )
    assert derived.vector_cross_section_slope_per_gev2 == pytest.approx(
        channel.vector_cross_section_slope_per_gev2
    )
    assert derived.factorized_cross_slope_per_gev2 == pytest.approx(
        2.0 * residue.amplitude_slope_per_gev2 + channel.vector_cross_section_slope_per_gev2
    )
    assert payload["cards"]["photoprod"] == pytest.approx([
        channel.pdg,
        channel.w0_gev,
        channel.vector_cross_section_slope_per_gev2,
        channel.alpha0.central,
        channel.alpha_prime_per_gev2.central,
    ], abs=5.0e-5)
    for model_name, first_pdg, second_pdg in channel.production_channels:
        model = payload["cards"]["models"][model_name]
        block = model[f"[{first_pdg},{second_pdg}]"]
        assert "coupling" not in block and "reference" not in block
        assert block["basis"] == "helicity"
        coupling = round(derived.gamma_pomeron_vector_coupling_per_gev, 9)
        assert block["helicity"] == [[-1, 0, pytest.approx(coupling), 0.0]]
        if model_name == "MP":
            assert block["polarization"] == {"mode": "none"}
        if model_name in {"MP", "XP"}:
            assert block["Lambda"] == pytest.approx(1.0)
        else:
            assert "Lambda" not in block
        assert set(model) == {f"[{first_pdg},{second_pdg}]"}


# Check a fixed cross section scale acts equally on coupling magnitudes and uncertainties
@pytest.mark.parametrize("scale", [0.5, 1.5])
def test_hera_norm_scale_relative_coupling_errors(scale):
    residue = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    channel = hera_inputs.light()[0]
    original = hera_model.derive(channel, residue)
    changed = hera_model.derive(
        replace(channel, normalization_scale=scale * channel.normalization_scale), residue
    )
    for field in (
        "gamma_pomeron_vector_coupling_per_gev",
        "coupling_stat_up_per_gev", "coupling_stat_down_per_gev",
        "coupling_syst_up_per_gev", "coupling_syst_down_per_gev",
    ):
        assert getattr(changed, field) == pytest.approx(math.sqrt(scale) * getattr(original, field))
    assert changed.reconstructed_forward_dsigma_dt_ub_per_gev2 == pytest.approx(
        scale * original.reconstructed_forward_dsigma_dt_ub_per_gev2
    )


# Check elastic and dissociative HERA cross sections stay fixed when the SOFT proton vertex changes
@pytest.mark.usefixtures("hera_channels")
def test_hera_dissociation_norm():
    residue = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    scaled = replace(residue, beam_residue_per_gev=1.5 * residue.beam_residue_per_gev)
    for channel in hera_couplings.channels():
        original = hera_model.derive(channel, residue)
        changed = hera_model.derive(channel, scaled)
        assert changed.gamma_pomeron_vector_coupling_per_gev == pytest.approx(
            original.gamma_pomeron_vector_coupling_per_gev / 1.5
        )
        diss = hera_couplings.dissociation(channel, changed)
        assert diss["photoprod_diss"] == pytest.approx(
            hera_couplings.dissociation(channel, original)["photoprod_diss"]
        )
        _, W, ratio, _, b, n, _, _ = diss["photoprod_diss"]
        forward = (
            hera_model.forward(
                changed.gamma_pomeron_vector_coupling_per_gev, scaled.beam_residue_per_gev
            )
            * (W / channel.w0_gev) ** changed.forward_w_exponent
            * ratio
        )
        if channel.dlog is not None:
            forward *= hera_model.dlog_factor(W, channel.w0_gev, *channel.dlog)
        if channel.pdg == 333:
            for t, xs, stat, syst_up, syst_down, mass_up, mass_down in diss["data"]:
                error = math.sqrt(stat**2 + max(syst_up, syst_down) ** 2 + max(mass_up, mass_down) ** 2)
                assert hera_fit.dissociation_spectrum(t, forward, b, n) == pytest.approx(xs, abs=error)
        else:
            assert hera_model.sigma_profile(forward, b, n, diss["tmax_GeV2"]) == pytest.approx(
                diss["sigma_ub"], rel=0.002
            )


# Check inherited proton transitions preserve forward ratios and the relative energy dependence
@pytest.mark.parametrize("pdg", [100443, 553, 100553, 200553])
@pytest.mark.parametrize("curve", [None, (11.0, 1.7)])
@pytest.mark.usefixtures("hera_channels")
def test_hera_dissociation_transfer(pdg, curve):
    residue = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    channels = {channel.pdg: channel for channel in hera_couplings.channels()}
    reference = replace(channels[443], dlog=curve)
    parent = hera_couplings.dissociation(reference, hera_model.derive(reference, residue))
    channel = replace(channels[pdg], dlog=curve if pdg == 100443 else None, diss_reference=reference)
    derived = hera_model.derive(channel, residue)
    original = hera_couplings.dissociation(channel, derived)
    varied_channel = replace(channel, normalization_scale=2.0 * channel.normalization_scale,
                             alpha0=replace(channel.alpha0, central=channel.alpha0.central + 0.1))
    varied = hera_couplings.dissociation(varied_channel, hera_model.derive(varied_channel, residue))
    for child in (original, varied):
        assert child["photoprod_diss"][2] == pytest.approx(parent["photoprod_diss"][2])
        assert child["photoprod_diss"][4:] == pytest.approx(parent["photoprod_diss"][4:])
    W0 = original["photoprod_diss"][1]
    assert varied["forward_ub_per_GeV2"] / original["forward_ub_per_GeV2"] == pytest.approx(
        2.0 * (W0 / channel.w0_gev) ** 0.4
    )
    # Reconstruct the pd/el ratio from the two separately normalized energy laws
    def ratio_at_w(state, payload, W):
        _, w0, ratio, delta, *_ = payload["photoprod_diss"]
        elastic = (W / w0)**(4 * (state.alpha0.central - 1))
        if state.dlog is not None:
            elastic *= hera_model.dlog_factor(W, state.w0_gev, *state.dlog) / hera_model.dlog_factor(w0, state.w0_gev, *state.dlog)
        dissociative = ratio * (W / w0)**delta
        for _, scale2, coefficient in payload.get("photoprod_diss_dlog", []):
            dissociative *= hera_model.dlog_factor(W, w0, scale2, coefficient)
        return dissociative / elastic

    for W in (0.7 * W0, 1.8 * W0):
        assert ratio_at_w(channel, original, W) == pytest.approx(ratio_at_w(reference, parent, W), rel=2e-4)
    assert varied["photoprod_diss"][3] - original["photoprod_diss"][3] == pytest.approx(0.4)


# Check every HERA production model stores only the physical coupling magnitude
def test_hera_model_helicity_rows_physical_norm():
    physical = 0.42
    channel = next(row for row in hera_inputs.light() if row.pdg == 113)
    models = hera_couplings.production(channel, physical)

    for model_name, first_pdg, second_pdg in channel.production_channels:
        block = models[model_name][f"[{first_pdg},{second_pdg}]"]
        rows = block["helicity"]
        active = [row[2] for row in rows if row[2] != 0.0]
        assert active == [pytest.approx(physical)]


# Check DL continuum and HERA transitions share one active bare proton scale
@pytest.mark.usefixtures("hera_channels")
def test_dl_hera_share_active_bare_pomeron_residue():
    dl = load_tool("DL_couplings")
    general_path = ROOT / "modeldata" / "TUNE0" / "GENERAL.json"

    configuration = dl.load_config(general_path)
    beam = dl.beam_couplings(configuration)
    bare = dl.bare_couplings(dl.coefficients(), beam)
    proton = hera_model.load_proton(general_path)
    channel = next(row for row in hera_couplings.channels() if row.pdg == 443)
    derived = hera_model.derive(channel, proton)
    payload = hera_couplings.output([(channel, derived)], proton)

    assert bare[("p", "P")] == pytest.approx(proton.beam_residue_per_gev)
    assert payload["amplitude_convention"]["Pomeron_beam_residue_t0_GeV_minus1"] == (
        pytest.approx(bare[("p", "P")])
    )
    fitted_forward = channel.forward_dsigma_dt_ub_per_gev2.central * channel.normalization_scale
    assert hera_model.forward(
        derived.gamma_pomeron_vector_coupling_per_gev,
        bare[("p", "P")],
    ) == pytest.approx(fitted_forward)


# Check SOFT vertex projection implements every transition profile at arbitrary N
def test_soft_residue_transition_form_factor_mode():
    soft = load_tool("lib.soft_exchange")
    model = {
        "GW": {"theta": [math.pi / 4.0]},
        "EXCHANGE": {
            "P": {
                "g": [[2.0, 1.0], [1.0, 4.0]],
                "ff": "P",
                "transition_ff": "geometric",
            }
        },
        "FF": {
            "P": {
                "type": "DPOW",
                "param": [[1.0, 2.0, 1.0], [2.0, 4.0, 1.0]],
            }
        },
    }
    proton = [math.sqrt(0.5), math.sqrt(0.5)]
    slopes = [1.5, 0.75]
    geometric_residue = sum(
        proton[row] * model["EXCHANGE"]["P"]["g"][row][column] * proton[column]
        for row in range(2)
        for column in range(2)
    )
    geometric_derivative = sum(
        proton[row]
        * model["EXCHANGE"]["P"]["g"][row][column]
        * proton[column]
        * (slopes[row] if row == column else 0.5 * (slopes[row] + slopes[column]))
        for row in range(2)
        for column in range(2)
    )

    geometric = soft.proton_vertex(model, "P")
    assert geometric.beam_residue_per_gev == pytest.approx(geometric_residue)
    assert geometric.amplitude_slope_per_gev2 == pytest.approx(
        geometric_derivative / geometric_residue
    )

    model["EXCHANGE"]["P"]["transition_ff"] = "arithmetic"
    arithmetic = soft.proton_vertex(model, "P")
    assert arithmetic == geometric

    model["EXCHANGE"]["P"]["transition_ff"] = "diagonal"
    diagonal = soft.proton_vertex(model, "P")
    diagonal_residue = sum(
        proton[index] ** 2 * model["EXCHANGE"]["P"]["g"][index][index] for index in range(2)
    )
    diagonal_derivative = sum(
        proton[index] ** 2 * model["EXCHANGE"]["P"]["g"][index][index] * slopes[index]
        for index in range(2)
    )
    assert diagonal.beam_residue_per_gev == pytest.approx(diagonal_residue)
    assert diagonal.amplitude_slope_per_gev2 == pytest.approx(
        diagonal_derivative / diagonal_residue
    )


# Check the rho fit against independent integration of the physical differential amplitude
def test_hera_rho_slope():
    import numpy as np
    from scipy.integrate import quad

    residue = hera_model.load_proton(ROOT / "modeldata/TUNE0/GENERAL.json")
    fit = hera_fit.light(residue, 113)
    folder = ROOT / "HEPData/PHOTOPROD/HEPData-ins1798511-v1-json"
    data = hera_reader.read_rho_wt(next(folder.glob("Table12-13,*")),
                                      folder / "Table12-13statisticalcorrelations.json")
    forward, slope, delta, alpha_prime = fit["parameters"]
    w0 = hera_inputs.light()[0].w0_gev
    predicted = [quad(lambda t, w=w: forward * hera_model.proton_profile(t, residue)**2 * np.exp(-slope*t)
                      * (w/w0)**(delta-4*alpha_prime*t), low, high)[0] / (high-low)
                 for (low, high), w in zip(data["binedges"], data["W"], strict=True)]
    np.testing.assert_allclose(fit["prediction"], predicted, rtol=1e-9)
    difference = np.asarray(predicted) - data["y"]
    assert fit["chi2"] == pytest.approx(difference @ np.linalg.solve(data["stat_cov"], difference))
    assert np.linalg.eigvalsh(fit["covariance"]).min() > 0
    assert np.all(np.asarray(fit["syst_up"]) >= 0)


# Check the measured excited-state ratio and integrated bottomonium energy power with actual profiles
@pytest.mark.usefixtures("hera_channels")
def test_hera_excited_states():
    channels = {row.pdg: row for row in hera_couplings.channels()}
    residue = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    measured_parent, measured_child = hera_inputs.charmonium()
    parent, child = channels[443], channels[100443]
    measured_ratio = hera_model.sigma_exp(
        measured_child.forward_dsigma_dt_ub_per_gev2.central, measured_child.b0_per_gev2.central, child.tmax_gev2
    ) / hera_model.sigma_exp(
        measured_parent.forward_dsigma_dt_ub_per_gev2.central, measured_parent.b0_per_gev2.central, child.tmax_gev2
    )
    actual_ratio = child.forward_dsigma_dt_ub_per_gev2.central * hera_model.profile_moments(
        residue, child.vector_cross_section_slope_per_gev2, child.tmax_gev2
    )[0] / (parent.forward_dsigma_dt_ub_per_gev2.central * hera_model.profile_moments(
        residue, parent.vector_cross_section_slope_per_gev2, child.tmax_gev2
    )[0])
    assert actual_ratio == pytest.approx(measured_ratio)
    for measured in hera_inputs.bottomonium():
        channel = channels[measured.pdg]
        derived = hera_model.derive(channel, residue)
        assert derived.matched_sigma_ub == pytest.approx(
            measured.forward_dsigma_dt_ub_per_gev2.central * hera_model.profile_moments(
                residue, measured.vector_cross_section_slope_per_gev2, measured.tmax_gev2)[0]
        )
        _, mean = hera_model.profile_moments(residue, channel.vector_cross_section_slope_per_gev2, channel.tmax_gev2)
        assert derived.forward_w_exponent - derived.shrinkage_per_gev2 * mean == pytest.approx(
            4.0 * (measured.alpha0.central - 1.0) - 4.0 * measured.alpha_prime_per_gev2.central * mean
        )


# Check finite-interval slope matching against an exactly exponential proton form factor
@pytest.mark.parametrize("tmax", [0.7, 2.3, None])
def test_hera_exp_projection(tmax):
    proton = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    proton = replace(proton, form_factor="EXP", parameters=((2.5,),), weights=((1.0,),), transition="diagonal")
    slope, norm = hera_model.match_slope(proton, 6.0, 0.0, tmax)
    assert slope == pytest.approx(6.0 - 2.5)
    assert norm == pytest.approx(1.0)
    with pytest.raises(ValueError, match="broader"):
        hera_model.match_slope(proton, 1.0, 0.0, tmax)


# Check finite-transfer Good Walker projections for every transition prescription
@pytest.mark.parametrize("transition", ["diagonal", "arithmetic", "geometric"])
def test_hera_proton_profile(transition):
    import numpy as np

    proton = hera_model.load_proton(ROOT / "modeldata" / "TUNE0" / "GENERAL.json")
    weights = np.array([[0.3, 0.1], [0.1, 0.5]])
    if transition == "diagonal":
        weights *= np.eye(2)
        weights /= weights.sum()
    proton = replace(proton, form_factor="EXP", parameters=((2.0,), (4.0,)),
                     weights=tuple(tuple(row) for row in weights), transition=transition)
    t = np.array([0.0, 0.1, 0.9, 2.5])
    first, second = np.exp(-t), np.exp(-2.0 * t)
    expected = weights[0, 0] * first + weights[1, 1] * second
    if transition == "arithmetic":
        expected += weights[0, 1] * (first + second)
    elif transition == "geometric":
        expected += 2.0 * weights[0, 1] * np.sqrt(first * second)
    np.testing.assert_allclose(hera_model.proton_profile(t, proton), expected)


# Check HERA push updates trajectory and coupling magnitudes without touching phases
@pytest.mark.usefixtures("hera_channels")
def test_hera_push_tune_card_comments_phases_format(tmp_path):
    json5 = pytest.importorskip("pyjson5")
    tune_dir = tmp_path / "TUNE0"
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune_dir)
    general_path = tune_dir / "GENERAL.json"
    rho_path = tune_dir / "RES" / "rho_770.json"
    general_source = re.sub(r'("photoprod"\s*:\s*\[\s*\[\s*113\s*,\s*)[0-9.]+',
                            r'\g<1>41.0', general_path.read_text(encoding="utf-8"), count=1)
    rho_source = rho_path.read_text(encoding="utf-8")
    card = json5.loads(rho_source)
    models = card["PARAM_RES"]["MODELS"]
    original = {
        "MP": str(models["MP"]["[22,991]"]["helicity"][0][2]),
        "GP": str(models["GP"]["[22,990]"]["helicity"][0][2]),
        "XP": str(models["XP"]["[22,991]"]["helicity"][0][2]),
    }
    mp_old = f"[[-1,0,{original['MP']},0.0]]"
    mp_pos = rho_source.index(mp_old, rho_source.index('"MP":'))
    rho_source = rho_source[:mp_pos] + "[[-1,0,0.5,0.0]]" + rho_source[mp_pos + len(mp_old):]
    gp_start = rho_source.index('"GP": {')
    gp_old = f"[[-1,0,{original['GP']},0.0]]"
    gp_pos = rho_source.index(gp_old, gp_start)
    rho_source = rho_source[:gp_pos] + "[[-1,0,0.6,0.0]]" + rho_source[gp_pos + len(gp_old) :]
    original_phase = "2.47740328385816"
    xp_old = f"[[-1,0,{original['XP']},0.0]]"
    xp_pos = rho_source.index(xp_old, rho_source.index('"XP":'))
    rho_source = rho_source[:xp_pos] + f"[[-1,0,0.7,{original_phase}]]" + rho_source[xp_pos + len(xp_old):]
    assert json5.loads(general_source)["PARAM_REGGE"]["photoprod"][0][1] == pytest.approx(41.0)
    assert '[[-1,0,0.5,0.0]]' in rho_source
    assert "[[-1,0,0.6,0.0]]" in rho_source
    general_path.write_text(general_source, encoding="utf-8")
    rho_path.write_text(rho_source, encoding="utf-8")
    before = {path: path.read_bytes() for path in tune_dir.rglob("*.json")}
    channel = next(row for row in hera_couplings.channels() if row.pdg == 113)
    residue = hera_model.load_proton(general_path)
    selected = [(channel, hera_model.derive(channel, residue))]
    previews = []

    applied = hera_tune.push_cards(
        general_path,
        selected,
        confirm=lambda rows: previews.extend(rows) or True,
    )

    after = {path: path.read_bytes() for path in tune_dir.rglob("*.json")}
    expected = dict(before)
    expected_general = json5.loads(general_source)
    cards = hera_couplings.channel_output(channel, selected[0][1])["cards"]
    for field in ("photoprod", "photoprod_diss"):
        rows = expected_general["PARAM_REGGE"][field]
        index = next(index for index, row in enumerate(rows) if row[0] == channel.pdg)
        rows[index] = cards[field]
    assert json5.loads(after[general_path].decode()) == expected_general
    expected[general_path] = after[general_path]
    physical = selected[0][1].gamma_pomeron_vector_coupling_per_gev
    mp_coupling = round(physical, 9)
    gp_coupling = round(physical, 9)
    xp_coupling = round(physical, 9)
    expected_rho = rho_source.replace("[[-1,0,0.5,0.0]]", f"[[-1,0,{mp_coupling},0.0]]", 1)
    expected_rho = expected_rho.replace("[[-1,0,0.6,0.0]]", f"[[-1,0,{gp_coupling},0.0]]", 1)
    expected_rho = expected_rho.replace(
        f"0.7,{original_phase}",
        f"{xp_coupling},{original_phase}",
        1,
    )
    expected[rho_path] = expected_rho.encode()
    assert applied
    assert after == expected
    assert {row[0] for row in previews} >= {
        "GENERAL.json:PARAM_REGGE.photoprod[113].W0",
        "RES/rho_770.json:PARAM_RES.MODELS.MP.[22,991].helicity[(-1, 0)].magnitude",
        "RES/rho_770.json:PARAM_RES.MODELS.GP.[22,990].helicity[(-1, 0)].magnitude",
        "RES/rho_770.json:PARAM_RES.MODELS.XP.[22,991].helicity[(-1, 0)].magnitude",
    }
    assert original_phase in rho_path.read_text(encoding="utf-8")


# Reject unavailable measurements before the HERA CLI can modify any tune card
@pytest.mark.parametrize("push_args", [(), ("--push",)])
def test_hera_missing_inputs_preserve_tune_cards(tmp_path, push_args):
    tune_dir = tmp_path / "TUNE0"
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune_dir)
    general_path = tune_dir / "GENERAL.json"
    before = {path: path.read_bytes() for path in tune_dir.rglob("*.json")}
    command = [
        sys.executable,
        str(TOOLS / "HERA_couplings.py"),
        *push_args,
        "--channel",
        "rho",
        "--tune-general",
        str(general_path),
    ]

    result = subprocess.run(
        command,
        cwd=ROOT,
        input="yes\n",
        check=False,
        capture_output=True,
        text=True,
    )

    assert result.returncode == 1
    assert "Measurement unavailable: a supported HEPData JSON input is required" in result.stderr
    after = {path: path.read_bytes() for path in tune_dir.rglob("*.json")}
    assert after == before


# Check HERA normalization requires the same unit modulus phase as the C++ reader
@pytest.mark.parametrize("mode", ["raw", "unknown"])
def test_hera_rejects_nonunit_photo_phase(tmp_path, mode):
    general = load_tune0_general()
    general["PARAM_REGGE"]["photoprod_eta_mode"] = mode
    path = tmp_path / "GENERAL.json"
    path.write_text(json.dumps(general), encoding="utf-8")
    with pytest.raises(ValueError, match="photoprod_eta_mode must be rotating_t0 or rotating"):
        hera_model.load_proton(path)


# Check HERA push refuses a local proton normalization override
def test_hera_push_rejects_proton_override(monkeypatch):
    hera = load_tool("HERA_couplings")
    monkeypatch.setattr(sys, "argv", ["HERA_couplings.py", "--push", "--g-p-pp", "8.3"])

    with pytest.raises(ValueError, match="derived from the mapped SOFT exchange"):
        hera.main()


# Check exact Jacob-Wick orthogonality and parity for integer and half integer spins
@pytest.mark.parametrize("spins", [(0, 0, 0), (1, 1, 0), (0, 1, 1), (2, 1, 1),
                                  (1, "1/2", "1/2"), ("1/2", "1/2", 1),
                                  ("3/2", 1, "1/2"), (2, "3/2", "1/2")])
def test_ls_helicity_matrix_is_unitary(spins):
    tool = load_tool("ls_helicity_conversion")
    J, s1, s2 = map(sp.Rational, spins)
    h_rows = tool.helicity_rows(J, s1, s2)
    for parity in (None, -1, 1):
        ls_rows = tool.ls_basis(J, s1, s2, parity, 1, 1)
        if not ls_rows:
            continue
        matrix = tool.transform(J, s1, s2, ls_rows, h_rows)
        n = len(ls_rows)
        assert sp.simplify(matrix.T * matrix) == sp.eye(n)
        if parity is None:
            assert matrix.rows == matrix.cols
            assert sp.simplify(matrix * matrix.T) == sp.eye(n)
        else:
            phase = parity * (-1) ** (s1 + s2 - J)
            for index, row in enumerate(h_rows):
                opposite = h_rows.index(tool.HelicityRow(-row.lambda1, -row.lambda2))
                assert sp.simplify(matrix[opposite, :] - phase * matrix[index, :]) == sp.zeros(1, n)


# Check forbidden angular momentum couplings fail before printing a conversion
@pytest.mark.parametrize("arguments", [
    ("--J", "1", "--s1", "1", "--s2", "0", "--ls", "0:0"),
    ("--J", "1", "--s1", "1", "--s2", "0", "--ls", "4:1"),
    ("--J", "1", "--s1", "1/2", "--s2", "0"),
])
def test_ls_helicity_rejects_forbidden_couplings(arguments):
    result = run_tool("ls_helicity_conversion", *arguments, check=False)
    assert result.returncode > 0, result.stdout + result.stderr
    assert "forbidden LS" in result.stderr or "spins must differ" in result.stderr


# Check an incomplete parity selection cannot silently remove the parity constraint
@pytest.mark.parametrize("flags", [("--p1", "1"), ("--p2", "-1"),
                                  ("--parity", "1"), ("--p1", "1", "--p2", "1")])
def test_ls_helicity_rejects_incomplete_parity(flags):
    result = run_tool("ls_helicity_conversion", "--J", "1", "--s1", "1", "--s2", "0",
                      *flags, check=False)
    assert result.returncode > 0, result.stdout + result.stderr
    assert "--parity requires both --p1 and --p2" in result.stderr


# Check the spin-two singlet uses a single spin representation for every helicity
def test_analytic_jz_scalar_singlet():
    tool = load_tool("analytic_Jz")
    model = tool.Model(0, 1, 2, 1, True, pole_spin=2)
    row = tool.LSRow(0, 0)
    values = [float(tool.ls_component(model, row, (m, m))) for m in range(-2, 3)]
    assert values == pytest.approx([sign / math.sqrt(5) for sign in (1, -1, 1, -1, 1)])


# Check every constructed pure-|Jz| representative remains normalized
def test_analytic_jz_representatives_normalized():
    tool = load_tool("analytic_Jz")
    result = tool.derive(tool.Model(1, 1, 1, 1, True, pole_spin=2))
    found = [sector for sector in result["sectors"] if sector["found"]]

    assert found
    for sector in found:
        assert sum(row[2] ** 2 for row in sector["helicity"]) == pytest.approx(1.0)


# Check pure Jz rows obey GP parity and identical exchange after card expansion
@pytest.mark.parametrize("spin,parity,naturality,identical", [
    (0, 1, 1, True), (1, 1, 1, True), (1, -1, 1, True),
    (2, 1, 1, True), (2, -1, 1, False), (3, 1, 1, True), (1, 1, -1, False),
])
def test_analytic_jz_card_symmetries(spin, parity, naturality, identical):
    tool = load_tool("analytic_Jz")
    result = tool.derive(tool.Model(spin, parity, 1, naturality, identical, pole_spin=2))
    cards = tool.cards([result])[result["J_P"]]
    parity_phase = parity * naturality * (-1) ** spin
    exchange_phase = (-1) ** spin
    found = [sector for sector in result["sectors"] if sector["found"]]
    assert found
    for sector in found:
        assert all(-math.pi <= row[3] < math.pi for row in sector["helicity"])
        dense = {(m1, m2): magnitude * cmath.exp(1j * phase)
                 for m1, m2, magnitude, phase in sector["helicity"]}
        for (m1, m2), value in dense.items():
            assert dense.get((-m1, -m2), 0.0) == pytest.approx(parity_phase * value)
            if identical:
                assert dense.get((m2, m1), 0.0) == pytest.approx(exchange_phase * value)
        expanded = {}
        for m1, m2, magnitude, phase in cards[f"abs_Jz_{sector['abs_Jz']}"]["helicity"]:
            value = magnitude * cmath.exp(1j * phase)
            assert (m1, m2) not in expanded
            expanded[m1, m2] = value
            expanded[-m1, -m2] = parity_phase * value
            if identical:
                expanded[m2, m1] = exchange_phase * value
                expanded[-m2, -m1] = parity_phase * exchange_phase * value
        assert expanded == pytest.approx(dense)


# Check the resonance builder emits one independent GP helicity row
def test_gp_reduced_helicity():
    payload = json.loads(
        run_tool(
            "resonance_card_builder",
            "--model",
            "GP",
            "--fuse",
            "22",
            "990",
            "--spinX2",
            "2",
            "--P",
            "-1",
            "--C",
            "-1",
            "--basis",
            "helicity",
            "--mmax",
            "3",
            "--helicity-row",
            "-1",
            "0",
            "0.42",
            "0.3",
            "--format",
            "card",
        ).stdout
    )

    block = payload["GP"]["[22,990]"]
    assert block["basis"] == "helicity"
    assert block["helicity"] == [[-1, 0, 0.42, 0.3]]
    assert "g_ls" not in block


# Check the resonance builder rejects a repeated GP parity orbit
def test_gp_parity_orbit():
    result = run_tool(
        "resonance_card_builder",
        "--model",
        "GP",
        "--fuse",
        "22",
        "990",
        "--spinX2",
        "2",
        "--P",
        "-1",
        "--C",
        "-1",
        "--basis",
        "helicity",
        "--mmax",
        "3",
        "--helicity-row",
        "-1",
        "0",
        "0.42",
        "0.3",
        "--helicity-row",
        "1",
        "0",
        "0.42",
        "0.3",
        check=False,
    )
    assert result.returncode > 0, result.stdout + result.stderr
    assert "repeat one symmetry orbit" in result.stderr


# Check the scalar GP pole contains every even L=S operator with one active row
def test_res_builder_complete_gp_ls():
    builder = load_tool("resonance_card_builder")
    tune = ROOT / "modeldata" / "TUNE0"
    particles = builder.load_particles(tune)
    _, regge = builder.load_regge(tune)
    pole = builder.gp_pole(particles[990], particles, regge)
    payload = json.loads(run_tool(
        "resonance_card_builder", "--model", "GP", "--fuse", "990", "990",
        "--spinX2", "0", "--P", "1", "--C", "1", "--active-ls", "0", "0",
        "--channel-mag", "0.42", "--channel-phase", "0.3", "--format", "card",
    ).stdout)
    rows = payload["GP"]["[990,990]"]["g_ls"]
    assert [tuple(row[:2]) for row in rows] == [(spin, spin) for spin in range(0, pole.spin2 + 1, 2)]
    assert rows[0][2:] == pytest.approx([0.42, 0.3])
    assert all(abs(row[2]) < 1e-15 for row in rows[1:])


# Check diphoton width normalization retains the complete scalar pole LS table
@pytest.mark.parametrize("model", ["XP", "GP"])
def test_res_builder_complete_diphoton_ls(model):
    payload = json.loads(run_tool(
        "resonance_card_builder", "--model", model, "--fuse", "22", "22",
        "--spinX2", "0", "--P", "1", "--C", "1", "--active-ls", "0", "0",
        "--channel-phase", "4", "--format", "card",
    ).stdout)
    rows = payload[model]["[22,22]"]["g_ls"]
    assert [row[:2] for row in rows] == [[0, 0], [2, 2]]
    assert rows[0][2] is None
    assert -math.pi <= rows[0][3] < math.pi
    assert cmath.exp(1j * rows[0][3]) == pytest.approx(cmath.exp(4j))
    assert rows[1][2:] == [0.0, 0.0]


# Check phase wrapping preserves complex couplings in both production bases
@pytest.mark.parametrize("model,basis", [("XP", "g_ls"), ("XP", "helicity"),
                                        ("GP", "g_ls"), ("GP", "helicity")])
def test_res_builder_canonical_phases(model, basis):
    flags = ["--helicity-row", "-1", "0", "0.5", "4"] if (
        model == "GP" and basis == "helicity"
    ) else ["--active-ls", "0", "1", "--channel-phase", "4"]
    arguments = ["--model", model, "--basis", basis,
        "--fuse", "22", "990" if model == "GP" else "995",
        "--spinX2", "2", "--P", "-1", "--C", "-1", "--format", "card"]
    payload = json.loads(run_tool("resonance_card_builder", *arguments, *flags).stdout)
    reference = json.loads(run_tool("resonance_card_builder", *arguments, *flags[:-1], "0").stdout)
    block = next(iter(payload[model].values()))
    field = "g_ls" if basis == "g_ls" else "helicity"
    rows = block[field]
    original = next(iter(reference[model].values()))[field]
    assert all(-math.pi <= row[3] < math.pi for row in rows)
    for row, initial in zip(rows, original, strict=True):
        actual = row[2] * cmath.exp(1j * row[3])
        expected = initial[2] * cmath.exp(1j * initial[3]) * cmath.exp(4j)
        assert actual == pytest.approx(expected)
        if row[2] <= 0.0:
            assert abs(row[3]) < 1e-15


# Check GP input cannot exceed the configured physical pole spin
@pytest.mark.parametrize("basis", ["g_ls", "helicity"])
def test_res_builder_gp_pole_bounds(basis):
    builder = load_tool("resonance_card_builder")
    tune = ROOT / "modeldata" / "TUNE0"
    particles = builder.load_particles(tune)
    _, regge = builder.load_regge(tune)
    pole = builder.gp_pole(particles[990], particles, regge)
    extra = pole.spin2 + 2 if basis == "g_ls" else pole.spin2 // 2 + 1
    flags = (["--active-ls", str(extra), str(extra)] if basis == "g_ls" else
             ["--mmax", str(extra), "--helicity-row", str(extra), str(extra), "1", "0"])
    result = run_tool(
        "resonance_card_builder", "--model", "GP", "--fuse", "990", "990",
        "--spinX2", "0", "--P", "1", "--C", "1", "--basis", basis,
        *flags, check=False,
    )
    assert result.returncode > 0, result.stdout + result.stderr
    assert "not allowed" in result.stderr or "pole spin" in result.stderr


# Check identical GP helicity steering accepts an independent nonzero orbit
def test_res_builder_identical_gp_helicity():
    payload = json.loads(run_tool(
        "resonance_card_builder", "--model", "GP", "--fuse", "990", "990",
        "--spinX2", "2", "--P", "1", "--C", "1", "--basis", "helicity",
        "--helicity-row", "-1", "0", "0.5", "0.25", "--format", "card",
    ).stdout)
    assert payload["GP"]["[990,990]"]["helicity"] == [[-1, 0, 0.5, 0.25]]


# Check identical exchange removes a diagonal helicity for odd mother spin
def test_res_builder_rejects_bose_zero():
    result = run_tool(
        "resonance_card_builder", "--model", "GP", "--fuse", "990", "990",
        "--spinX2", "2", "--P", "-1", "--C", "1", "--basis", "helicity",
        "--helicity-row", "0", "0", "1", "0", check=False,
    )
    assert result.returncode > 0, result.stdout + result.stderr
    assert "forces this helicity to zero" in result.stderr


# Check the resonance builder consumes the complete validated SOFT snapshot
def test_res_builder_invalid_soft_model(tmp_path):
    builder = load_tool("resonance_card_builder")
    general = load_tune0_general()
    active_name = general["PARAM_SOFT"]["active_model"]
    general["PARAM_SOFT"]["MODEL"][active_name]["EIKONAL"]["helicity"] = "scalar"
    (tmp_path / "GENERAL.json").write_text(json.dumps(general), encoding="utf-8")

    with pytest.raises(TypeError, match="EIKONAL.helicity must be a boolean"):
        builder.load_regge(tmp_path)


# Check the resonance builder emits direct XP covariant couplings
def test_res_builder_emits_covariant_xp_schema():
    payload = json.loads(
        run_tool(
            "resonance_card_builder",
            "--model",
            "XP",
            "--fuse",
            "22",
            "991",
            "--spinX2",
            "2",
            "--P",
            "-1",
            "--C",
            "-1",
            "--active-ls",
            "0",
            "1",
            "--channel-mag",
            "0.42",
            "--channel-phase",
            "0.3",
            "--Lambda",
            "1.5",
            "--format",
            "card",
        ).stdout
    )

    model = payload["XP"]
    block = model["[22,991]"]
    assert set(model) == {"[22,991]"}
    assert block["basis"] == "g_ls"
    assert block["Lambda"] == pytest.approx(1.5)
    assert block["g_ls"] == [[0, 1, 0.42, 0.3], [2, 1, 0.0, 0.0]]
    assert "phase" not in block
    assert "alpha_ls" not in block
    assert "helicity" not in block


# Check the resonance builder emits direct XP helicity input
def test_res_builder_emits_xp_helicity_schema():
    payload = json.loads(
        run_tool(
            "resonance_card_builder",
            "--model",
            "XP",
            "--fuse",
            "22",
            "991",
            "--spinX2",
            "2",
            "--P",
            "-1",
            "--C",
            "-1",
            "--active-ls",
            "0",
            "1",
            "--basis",
            "helicity",
            "--format",
            "card",
        ).stdout
    )

    block = payload["XP"]["[22,991]"]
    assert block["basis"] == "helicity"
    assert "reference" not in block
    assert block["Lambda"] == pytest.approx(1.0)
    assert block["helicity"] == [[-1, 0, 1.0, 0.0]]
    assert "g_ls" not in block


# Check the XP card builder rejects forbidden real-photon quantum numbers
@pytest.mark.parametrize(
    ("spin_x2", "cparity", "message"),
    [
        (2, 1, "cannot produce a spin 1 resonance"),
        (0, -1, "requires C = \\+1"),
    ],
)
def test_diphoton_quantum_numbers(
    spin_x2, cparity, message
):
    command = [
        sys.executable,
        str(TOOLS / "resonance_card_builder.py"),
        "--model",
        "XP",
        "--fuse",
        "22",
        "22",
        "--spinX2",
        str(spin_x2),
        "--P",
        "1",
        "--C",
        str(cparity),
        "--format",
        "card",
    ]
    result = subprocess.run(command, cwd=ROOT, check=False, capture_output=True, text=True)

    assert result.returncode == 1
    assert re.search(message, result.stderr)


# Check transverse photon projection removes LS rows without pole support
def test_res_builder_rejects_xp_photon_scalar_j0():
    command = [
        sys.executable,
        str(TOOLS / "resonance_card_builder.py"),
        "--model",
        "XP",
        "--fuse",
        "22",
        "991",
        "--spinX2",
        "0",
        "--P",
        "1",
        "--C",
        "-1",
        "--format",
        "card",
    ]
    result = subprocess.run(command, cwd=ROOT, check=False, capture_output=True, text=True)

    assert result.returncode == 1
    assert "no allowed LS rows" in result.stderr


# Check each card-producing CLI emits parseable compact JSON
@pytest.mark.parametrize(
    ("name", "arguments", "required_key"),
    [
        ("DL_couplings", ("--format", "cards", "--model", "GP"), "CON"),
        ("HERA_couplings", ("--format", "cards"), "PARAM_REGGE"),
        ("analytic_Jz", ("--J", "0", "--mmax", "1", "--pole-spin", "2", "--format", "cards"), "0+"),
        (
            "resonance_card_builder",
            (
                "--model",
                "GP",
                "--fuse",
                "22",
                "990",
                "--spinX2",
                "2",
                "--P",
                "-1",
                "--C",
                "-1",
                "--basis",
                "helicity",
                "--helicity-row",
                "1",
                "0",
                "1.0",
                "0.0",
                "--format",
                "card",
            ),
            "GP",
        ),
    ],
)
def test_card_tool_output_is_valid_json(name, arguments, required_key, request):
    if name == "HERA_couplings":
        request.getfixturevalue("hera_channels")
    payload = json.loads(run_tool(name, *arguments).stdout)

    assert required_key in payload


# Check the finite power profile integral through the logarithmic limit
@pytest.mark.parametrize("power", [0.5, 1.0 - 1e-8, 1.0, 1.0 + 1e-8, 2.0])
def test_hera_power_profile_integral_unit_power(power):
    integrate = pytest.importorskip("scipy.integrate")
    slope, tmax, forward = 3.2, 0.7, 2.5
    expected = forward * integrate.quad(
        lambda t: (1.0 + slope * t / power) ** (-power), 0.0, tmax
    )[0]
    assert hera_model.sigma_profile(forward, slope, power, tmax) == pytest.approx(
        expected, rel=1e-12
    )


# Check finite cross sections as the slope or measured interval tends to zero
def test_hera_profile_integrals_zero_slope_interval():
    integrate = pytest.importorskip("scipy.integrate")
    forward = 2.5
    for slope in (0.0, 1e-18, 1e-8, 3.2):
        for tmax in (0.0, 1e-18, 0.7):
            exponential = forward * integrate.quad(
                lambda t, b: math.exp(-b * t), 0.0, tmax, args=(slope,)
            )[0]
            assert hera_model.sigma_exp(forward, slope, tmax) == pytest.approx(
                exponential, rel=1e-12, abs=0.0
            )
            for power in (0.0, 0.5, 1.0, 2.0, 1e12):
                expected = exponential if power <= 0.0 else forward * integrate.quad(
                    lambda t, b, n: math.exp(-n * math.log1p(b * t / n)),
                    0.0, tmax, args=(slope, power)
                )[0]
                assert hera_model.sigma_profile(forward, slope, power, tmax) == pytest.approx(
                    expected, rel=1e-12, abs=0.0
                )


# Check parity disabled helicity steering retains both vector spin orientations
def test_res_builder_helicities_without_parity():
    payload = json.loads(run_tool(
        "resonance_card_builder", "--model", "XP", "--fuse", "113", "333",
        "--spinX2", "0", "--P", "1", "--C", "1", "--active-ls", "0", "0",
        "--basis", "helicity", "--no-p", "--format", "card",
    ).stdout)
    block = payload["XP"]["[113,333]"]
    assert block["CP"][1] is False
    rows = {tuple(row[:2]): row[2] for row in block["helicity"]}
    assert set(rows) == {(-1, -1), (0, 0), (1, 1)}
    assert rows[(-1, -1)] == pytest.approx(rows[(1, 1)])
    assert sum(value**2 for value in rows.values()) == pytest.approx(3.0)


# Check GP rejects symmetry flags that its C++ reader cannot accept
@pytest.mark.parametrize("flag", ["--no-c", "--no-p"])
def test_res_builder_rejects_gp_symmetry_flags(flag):
    result = run_tool(
        "resonance_card_builder", "--model", "GP", "--fuse", "990", "990",
        "--spinX2", "0", "--P", "1", "--C", "1", flag, check=False,
    )
    assert result.returncode > 0, result.stdout + result.stderr
    assert "GP requires CP = [true,true]" in result.stderr


# Follow references when physics tools update a shared coupling in another card
def test_tool_push_updates_shared_reference(tmp_path):
    push = load_tool("lib.push")
    source = tmp_path / "source.json"
    alias = tmp_path / "alias.json"
    source.write_text('{"rows": [[0,0,1.0,0.0]]}\n', encoding="utf-8")
    raw = '{"vertex": {"$ref": "source.json#/rows"}}\n'
    alias.write_text(raw, encoding="utf-8")
    updates = {alias: [push.ScalarUpdate(("vertex", 0, 2), "coupling", 2.5)]}
    assert push.push_json5_updates(updates, confirm=lambda rows: True)
    assert alias.read_text(encoding="utf-8") == raw
    assert json.loads(source.read_text())["rows"][0][2] == pytest.approx(2.5)


# Propagate the measured ZEUS normalization uncertainties without inventing a shape covariance
def test_hera_zeus_covariance():
    import numpy as np

    spectra = hera_fit.zeus_spectra()
    for _, data, _, _, covariance, common in spectra:
        shifts = data["y"][:, None] * np.column_stack(list(common.values()))
        np.testing.assert_allclose(covariance, data["stat_cov"] + shifts @ shifts.T)
        assert np.all(np.diag(covariance) > data["y_err"]**2)
        assert all(np.all(value > 0) for value in common.values())
        assert data["bins"][0] == pytest.approx(0.0, abs=1e-14)
        original = hepdata.read_table(data["source"])
        np.testing.assert_allclose(data["x"], [hepdata.number(row["x"][0]["value"]) for row in original["values"]])


# Check the joint fit retains the measured ZEUS normalization constraints
@pytest.mark.usefixtures("hera_channels")
def test_hera_zeus_fit():
    fit = hera_fit.jpsi(hera_model.load_proton(ROOT / "modeldata/TUNE0/GENERAL.json"))
    assert {"zeus_muon", "zeus_electron"}.issubset(fit["nuisance_pulls"])
    assert not fit["energy_used_in_fit"]
    assert fit["ndf"] == sum(sum(row["shape_fit_bins"]) for row in fit["spectra"]) - len(fit["parameters"])
    assert fit["chi2"] == pytest.approx(sum(row["fit_chi2"] for row in fit["spectra"]) + fit["nuisance_chi2"])


# Transpose the original two-dimensional H1 measurement without duplicating or changing values
@pytest.mark.usefixtures("hera_channels")
def test_hera_h1_differential_reader():
    import numpy as np

    spectra = hera_fit.h1_wt_spectra()
    original = hepdata.read_table(spectra[0][1]["source"])
    assert len(spectra) == len(original["values"])
    for (_, data, energy, weight, covariance, common), row in zip(spectra, original["values"], strict=True):
        np.testing.assert_allclose(data["y"], [hepdata.number(value["value"]) for value in row["y"]])
        np.testing.assert_allclose(np.diag(covariance), data["y_err_stat"]**2 + data["y_err_syst"]**2)
        assert energy[0] == pytest.approx(hepdata.number(row["x"][0]["value"]))
        assert weight.sum() == pytest.approx(1.0)
        assert "rher1" in common
        np.testing.assert_allclose(data["binedges"][:-1, 1], data["binedges"][1:, 0])
    for (_, data, _, _, covariance, common), (_, diagonal, *_) in zip(hera_fit.h1_wt_spectra(True), spectra, strict=True):
        np.testing.assert_allclose(np.diag(covariance), diagonal["y_err"]**2)
        shift = data["y"][:, None] * np.column_stack(list(common.values()))
        np.testing.assert_allclose(covariance, data["stat_cov"] + shift @ shift.T)


# Reject unsupported covariance inputs instead of substituting independent errors
@pytest.mark.parametrize("reader", [hera_fit.h1_spectra, hera_fit.h1_wt_spectra,
                                    hera_fit.jpsi_energy, hera_fit.jpsi_dissociation])
def test_hera_missing_covariance(reader):
    with pytest.raises(ValueError, match="H1 covariance fit unavailable: a supported HEPData JSON covariance input is required"):
        reader()




# Reconstruct every fitted integrated cross section from the physical amplitude
@pytest.mark.usefixtures("hera_channels")
def test_hera_jpsi_energy_amplitude_normalization():
    import numpy as np
    from scipy.integrate import quad

    proton = hera_model.load_proton(ROOT / "modeldata/TUNE0/GENERAL.json")
    channel = next(row for row in hera_couplings.channels(proton) if row.pdg == 443)
    derived = hera_model.derive(channel, proton)
    fit = channel.fit
    forward = 1000.0 * hera_model.forward(
        derived.gamma_pomeron_vector_coupling_per_gev, proton.beam_residue_per_gev)
    predicted = []
    for _, energies, _, _, _, tmax in hera_fit.jpsi_energy():
        for w in energies:
            power = (w / channel.w0_gev)**derived.forward_w_exponent
            if channel.dlog is not None:
                power *= hera_model.dlog_factor(w, channel.w0_gev, *channel.dlog)
            slope = channel.vector_cross_section_slope_per_gev2 + 4 * channel.alpha_prime_per_gev2.central * math.log(w / channel.w0_gev)
            integral = quad(lambda t, b=slope: math.exp(-b*t) * hera_model.proton_profile(t, proton)**2,
                            0, np.inf if tmax is None else tmax, epsabs=1e-10)[0]
            predicted.append(forward * power * integral)
    np.testing.assert_allclose(predicted, fit["energy_prediction"], rtol=2e-8)
    covariance = np.asarray(fit["covariance"])
    np.testing.assert_allclose(covariance, covariance.T, rtol=1e-10)
    assert np.all(np.linalg.eigvalsh(covariance) > 0)


# Check the fitted histogram prediction against independent integration of the physical vertex
@pytest.mark.usefixtures("hera_channels")
def test_hera_jpsi_bin_integrals():
    import numpy as np
    from scipy.integrate import quad

    proton = hera_model.load_proton(ROOT / "modeldata/TUNE0/GENERAL.json")
    fit = hera_fit.jpsi(proton)
    forward, slope, delta, alpha_prime = fit["parameters"]
    for report, (_, data, energy, flux, _, _) in zip(
            fit["spectra"], hera_fit.jpsi_spectra(), strict=True):
        predicted = hera_fit.jpsi_spectrum(fit["parameters"], proton, data, energy, flux)
        integrals = [quad(lambda t, energy=energy, flux=flux: forward * float(hera_model.proton_profile(t, proton))**2
                         * math.exp(-slope * t)
                         * np.sum((energy / 100.0)**(delta - 4 * alpha_prime * t) * flux),
                         low, high, epsabs=1e-9)[0] for low, high in data["binedges"]]
        np.testing.assert_allclose(predicted * data["binwidth"], integrals, rtol=1e-9)
        np.testing.assert_allclose(report["prediction"], predicted, rtol=1e-12)


# Reject angular momentum labels outside the SU(2) representation before factorial evaluation
def test_cg_quantum_numbers():
    for values in ((0.5, 0.5, 0.5, -0.5, 0.5, 0.0), (1.0, 1.0, 0.5, -0.5, 1.0, 0.0),
                   (0.3, 0.3, 0.3, -0.3, 0.0, 0.0)):
        assert abs(pole_math.clebsch_gordan(*values)) < 1e-14
    for j1, j2, j in ((0.5, 0.5, 1.0), (1.0, 0.5, 1.5), (2.0, 1.0, 2.0)):
        for m1 in (-j1, j1):
            for m2 in (-j2, j2):
                expected = wigner.clebsch_gordan(j1, j2, j, m1, m2, m1 + m2)
                assert pole_math.clebsch_gordan(j1, j2, m1, m2, j, m1 + m2) == pytest.approx(float(expected))


# Ensure a measurement lookup cannot read a saved copy or choose an ambiguous original table
def test_hera_table_selection(tmp_path, monkeypatch):
    from develop.tools.lib.hera import inputs as inputs

    original = inputs.table(451266, "Table9.json")
    folder = tmp_path / original.parent.relative_to(ROOT)
    folder.mkdir(parents=True)
    shutil.copyfile(original, folder / "Table9.json._old")
    monkeypatch.setattr(inputs.model, "REPOSITORY_ROOT", tmp_path)
    with pytest.raises(ValueError, match="found 0"):
        inputs.table(451266, "Table9*")
    shutil.copyfile(original, folder / original.name)
    assert inputs.measurement(451266, "Table9*", 0).central > 0
    shutil.copyfile(original, folder / "Table9-copy.json")
    with pytest.raises(ValueError, match="found 2"):
        inputs.table(451266, "Table9*")


# Ensure invalid fit values cannot be serialized as nonstandard JSON tokens
def test_tool_json_finite():
    with pytest.raises(ValueError):
        tool_common.dumps({"couplings": [1.0, math.nan]})


# Check derived energy curvature is inserted, updated and removed with the corresponding couplings
@pytest.mark.usefixtures("hera_channels")
def test_hera_push_curvature(tmp_path):
    json5 = pytest.importorskip('pyjson5')
    tune = tmp_path / 'tune'
    shutil.copytree(ROOT / 'modeldata/TUNE0', tune)
    general_path = tune / 'GENERAL.json'
    original = json5.loads(general_path.read_text())
    for field in ('photoprod_dlog', 'photoprod_diss_dlog'):
        original['PARAM_REGGE'].pop(field, None)
    general_path.write_text(json.dumps(original, indent=2))
    proton = hera_model.load_proton(general_path)
    channel = replace(hera_inputs.charmonium()[0], dlog=(11.0, 1.7))
    selected = [(channel, hera_model.derive(channel, proton))]
    before = {path: path.read_bytes() for path in tune.rglob('*.json')}
    assert not hera_tune.push_cards(general_path, selected, confirm=lambda rows: False)
    assert before == {path: path.read_bytes() for path in tune.rglob('*.json')}
    assert hera_tune.push_cards(general_path, selected, confirm=lambda rows: True)
    updated = json5.loads(general_path.read_text())
    assert updated['PARAM_REGGE']['photoprod_dlog'] == [[channel.pdg, *channel.dlog]]
    assert updated['PARAM_SOFT'] == original['PARAM_SOFT']
    assert not hera_tune.push_cards(general_path, selected, confirm=lambda rows: True)
    channel = replace(channel, dlog=None)
    selected = [(channel, hera_model.derive(channel, proton))]
    assert hera_tune.push_cards(general_path, selected, confirm=lambda rows: True)
    assert json5.loads(general_path.read_text())['PARAM_REGGE']['photoprod_dlog'] == []


# Reject updates when a derivation input changes during confirmation
def test_dl_read_dependencies(tmp_path):
    dl = load_tool("DL_couplings")
    general, card = tmp_path / "GENERAL.json", tmp_path / "CON_GP.json"
    general.write_text('{"value":1.0}\n')
    card.write_text('{"g":1.0}\n')
    before = card.read_bytes()
    reader = dl.CardReader({general, card})
    dl.read_card(general, reader)
    updates = {card: [dl.ScalarUpdate(("g",), "g", 2.0)]}

    # Change a read dependency after the update preview
    def confirm(rows):
        general.write_text('{"value":2.0}\n')
        return True

    with pytest.raises(ValueError, match="Card changed after preview"):
        dl.push_json5_updates(updates, reader=reader, confirm=confirm)
    assert card.read_bytes() == before
