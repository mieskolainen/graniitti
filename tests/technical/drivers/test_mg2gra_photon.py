# Regression tests for automated MadGraph photon hard processes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import importlib
import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
MODULE_DIR = ROOT / "develop/MG2GRA"
CARDS_ROOT = ROOT / "MG5cards"
PHOTON_PROCESSES = (
    ("yy_ll", [-11, 11], 1, 16, 4, 1, 0, 2),
    ("yy_uubarg", [2, -2, 21], 1, 32, 4, 1, 1, 2),
    ("yy_uubargg", [2, -2, 21, 21], 2, 64, 8, 2, 2, 2),
)


# Compute one projection-owned standalone photon color sidecar
def photon_color_path(output_layout, name):
    return CARDS_ROOT / output_layout.color_path("photon", name)


# Load one converter module without relying on the working directory
def load_module(name):
    qualified = name if name == "regenerate" else f"modules.{name}"
    sys.path.insert(0, str(MODULE_DIR))
    try:
        module = importlib.import_module(qualified)
    finally:
        sys.path.remove(str(MODULE_DIR))
    return module


# Preserve the generated model mass table through direct and nested kinematics
@pytest.mark.parametrize("decay_chain", [False, True])
@pytest.mark.parametrize("photon", [False, True])
def test_generated_kinematics_model_masses(decay_chain, photon):
    converter = load_module("standalone")
    source = converter.make_setup_kinematics("Amplitude", decay_chain, photon)
    assert "OnShellFinal(final_buffer, mME)" in source
    assert "mME.push_back" not in source
    assert "mME.assign" not in source
    assert "ExternalMass" not in source
    if decay_chain:
        assert "StableDecayLeaves" in source


# Verify every photon sidecar has complete stable-state and color process_data
@pytest.mark.parametrize(
    (
        "name",
        "final_pdgs",
        "ncolor",
        "nhelicity",
        "process_denominator",
        "final_symmetry_factor",
        "alpha_s_power",
        "alpha_qed_power",
    ),
    PHOTON_PROCESSES,
)
def test_photon_metadata_is_complete(
    name,
    final_pdgs,
    ncolor,
    nhelicity,
    process_denominator,
    final_symmetry_factor,
    alpha_s_power,
    alpha_qed_power,
):
    output_layout = load_module("output_layout")
    process_data = json.loads(
        photon_color_path(output_layout, name).read_text()
    )

    assert process_data["incoming_pdgs"] == [22, 22]
    assert process_data["final_pdgs"] == final_pdgs
    assert process_data["ncolor"] == ncolor
    assert process_data["rank"] > 0
    assert len(process_data["helicities"]) == nhelicity
    assert process_data["process_denominator"] == process_denominator
    assert process_data["final_symmetry_factor"] == final_symmetry_factor
    assert process_data["initial_state_denominator"] == 4
    assert (
        process_data["initial_state_denominator"] * process_data["final_symmetry_factor"]
        == process_data["process_denominator"]
    )
    assert process_data["default_alpha_s"] > 0.0
    assert process_data["default_alpha_qed"] > 0.0
    assert process_data["alpha_s_power"] == alpha_s_power
    assert process_data["alpha_qed_power"] == alpha_qed_power
    assert process_data["color_metric_residual"] < 1.0e-12
    assert process_data["factorization_residual"] < 1.0e-12


# Verify the colored proof process carries a complete exact and shower color map
def test_photon_su3_color_flows():
    output_layout = load_module("output_layout")
    process_data = json.loads(
        photon_color_path(output_layout, "yy_uubarg").read_text()
    )

    assert process_data["final_color_representations"] == [3, -3, 8]
    assert process_data["ncolor"] == 1
    assert process_data["rank"] == 1
    assert len(process_data["flow_candidates"]) == 1
    assert len(process_data["flow_candidates"][0]) == 3
    assert process_data["flow_weights"] == [[1.0]]
    assert process_data["projectors"] == [[[2.0, 0.0]]]


# Verify a multi-dimensional color basis is fully factorized and partitioned
def test_photon_projection_rank():
    output_layout = load_module("output_layout")
    process_data = json.loads(
        photon_color_path(output_layout, "yy_uubargg").read_text()
    )

    assert process_data["final_color_representations"] == [3, -3, 8, 8]
    assert process_data["ncolor"] == 2
    assert process_data["rank"] == 2
    assert len(process_data["projectors"]) == 2
    assert all(len(row) == 2 for row in process_data["projectors"])
    assert len(process_data["flow_candidates"]) == 2
    assert process_data["flow_weights"] == [[1.0, 0.0], [0.0, 1.0]]


# Reject type, dimension and non-finite mutations in photon sidecars
@pytest.mark.parametrize(
    "mutation",
    (
        "denominator_type",
        "decay_flag_type",
        "alpha_nan",
        "alpha_qed_nan",
        "alpha_power_negative",
        "alpha_power_mismatch",
        "alpha_qed_power_mismatch",
        "helicity_type",
        "projector_inf",
        "weight_width",
        "weight_negative",
        "weight_partition",
    ),
)
def test_photon_registry_invalid_sidecar_values(tmp_path, mutation):
    registry = load_module("photon_registry")
    output_layout = load_module("output_layout")
    process_data = json.loads(
        photon_color_path(output_layout, "yy_uubarg").read_text()
    )
    if mutation == "denominator_type":
        process_data["process_denominator"] = 4.0
    elif mutation == "decay_flag_type":
        process_data["has_decay_chain"] = 0
    elif mutation == "alpha_nan":
        process_data["default_alpha_s"] = math.nan
    elif mutation == "alpha_qed_nan":
        process_data["default_alpha_qed"] = math.nan
    elif mutation == "alpha_power_negative":
        process_data["alpha_s_power"] = -1
    elif mutation == "alpha_power_mismatch":
        process_data["alpha_s_power"] = 2
    elif mutation == "alpha_qed_power_mismatch":
        process_data["alpha_qed_power"] = 3
    elif mutation == "helicity_type":
        process_data["helicities"][0][0] = True
    elif mutation == "projector_inf":
        process_data["projectors"][0][0][0] = math.inf
    elif mutation == "weight_width":
        process_data["flow_weights"][0].clear()
    elif mutation == "weight_negative":
        process_data["flow_weights"][0][0] = -1.0
    else:
        process_data["flow_weights"][0][0] = 0.5

    directory = tmp_path / output_layout.cards_directory("photon", "yy_uubarg")
    directory.mkdir(parents=True)
    (directory / output_layout.COLOR_DATA).write_text(json.dumps(process_data))
    entry = {
        "name": "yy_uubarg",
        "process": process_data["process"],
        "projection": "photon",
    }

    with pytest.raises(RuntimeError, match="Invalid|Non-unit"):
        registry.load_process_data(tmp_path, entry)


# Verify registry generation is deterministic and rejects ambiguous matches
def test_photon_registry_signatures():
    regenerate = load_module("regenerate")
    registry = load_module("photon_registry")
    manifest = regenerate.load_manifest(MODULE_DIR / "processes.json")

    first = registry.generate_registry(CARDS_ROOT, manifest)
    second = registry.generate_registry(CARDS_ROOT, manifest)
    assert first == second

    entries = [
        (
            entry,
            registry.load_process_data(CARDS_ROOT, entry),
        )
        for entry in manifest["processes"]
        if entry["projection"] == "photon"
    ]
    duplicate = (copy.deepcopy(entries[0][0]), copy.deepcopy(entries[0][1]))
    duplicate[0]["name"] = "duplicate"
    with pytest.raises(RuntimeError, match="Duplicate photon"):
        registry.validate_signatures(entries + [duplicate])

    direct = (
        {"name": "direct", "process": "a a > e+ e-"},
        {"final_pdgs": [-11, 11]},
    )
    cascade = (
        {"name": "cascade", "process": "a a > z, z > e+ e-"},
        {"final_pdgs": [-11, 11]},
    )
    registry.validate_signatures([direct, cascade])

# Accept MG5 bounds, exact orders and mixed orders without guessing amplitude powers
def test_photon_coupling_orders():
    regenerate = load_module("regenerate")
    registry = load_module("photon_registry")
    manifest = regenerate.load_manifest(MODULE_DIR / "processes.json")
    family = copy.deepcopy(registry.photon_families(manifest)[0])
    for processes in (["a a > j j"], ["a a > e+ e- QED==2"],
                      ["a a > u u~ QCD=0", "a a > u u~ g QCD=1"]):
        family["processes"] = processes
        assert registry.photon_families({"families": [family]}) == [family]
