# Regression tests for exact finite-Nc Durham color generation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[3]
MG2GRA = ROOT / "develop" / "MG2GRA"
CARDS_ROOT = ROOT / "MG5cards"
sys.path.insert(0, str(MG2GRA))

from modules import durham_registry, mg5_color, output_layout


# Deserialize one generated complex projector matrix
def _projectors(metadata):
    return [[complex(real, imaginary) for real, imaginary in row] for row in metadata["projectors"]]


# Form P dagger P from projector rows
def _gram(projectors):
    width = len(projectors[0])
    return [
        [sum(row[i].conjugate() * row[j] for row in projectors) for j in range(width)]
        for i in range(width)
    ]


# Load one checked-in Durham metadata sidecar
def _metadata(name):
    path = CARDS_ROOT / output_layout.color_path("durham", name)
    return json.loads(path.read_text())


# Check two complex matrices element by element
def _require_matrix_close(left, right, tolerance=1.0e-12):
    assert len(left) == len(right)
    for left_row, right_row in zip(left, right, strict=False):
        assert len(left_row) == len(right_row)
        for left_value, right_value in zip(left_row, right_row, strict=False):
            assert left_value == pytest.approx(right_value, abs=tolerance)


# Verify the deterministic factorization for a nontrivial complex PSD matrix
def test_pivoted_cholesky_reconstructs_complex_gram():
    seed = [
        [2.0 + 0.0j, 1.0 - 2.0j, -0.5 + 0.25j],
        [0.0 + 0.0j, 3.0 + 0.0j, 1.5 - 0.75j],
    ]
    gram = _gram(seed)
    first = mg5_color.pivoted_cholesky(gram)
    second = mg5_color.pivoted_cholesky(gram)

    assert first == second
    assert len(first) == 2
    _require_matrix_close(_gram(first), gram)


# Verify known analytic SU(3) projectors and the corrected qqbar-g basis order
def test_su3_projectors():
    qqbar = _metadata("gg_uubar")
    qqbar_expected = [[complex(math.sqrt(2.0 / 3.0))] * 2]
    _require_matrix_close(_gram(_projectors(qqbar)), _gram(qqbar_expected))

    root_two = math.sqrt(2.0)
    qqbarg_expected = [
        [
            4.0 / (3.0 * root_two),
            -1.0 / (6.0 * root_two),
            4.0 / (3.0 * root_two),
            -1.0 / (6.0 * root_two),
            4.0 / (3.0 * root_two),
            4.0 / (3.0 * root_two),
        ]
    ]
    qqbarg = _metadata("gg_uubarg")
    _require_matrix_close(_gram(_projectors(qqbarg)), _gram(qqbarg_expected))

    inv_sqrt3 = 1.0 / math.sqrt(3.0)
    inv8_sqrt3 = inv_sqrt3 / 8.0
    d_leading = math.sqrt(5.0 / 27.0)
    d_sublead = d_leading / 8.0
    ggg_expected = [
        [
            inv_sqrt3,
            -inv_sqrt3,
            -inv_sqrt3,
            inv_sqrt3,
            inv_sqrt3,
            -inv_sqrt3,
            -inv8_sqrt3,
            inv8_sqrt3,
            -inv8_sqrt3,
            inv_sqrt3,
            inv8_sqrt3,
            -inv_sqrt3,
            inv8_sqrt3,
            -inv8_sqrt3,
            inv8_sqrt3,
            -inv_sqrt3,
            -inv8_sqrt3,
            inv_sqrt3,
            -inv8_sqrt3,
            inv8_sqrt3,
            -inv8_sqrt3,
            inv_sqrt3,
            inv8_sqrt3,
            -inv_sqrt3,
        ],
        [
            d_leading,
            d_leading,
            d_leading,
            d_leading,
            d_leading,
            d_leading,
            -d_sublead,
            -d_sublead,
            -d_sublead,
            d_leading,
            -d_sublead,
            d_leading,
            -d_sublead,
            -d_sublead,
            -d_sublead,
            d_leading,
            -d_sublead,
            d_leading,
            -d_sublead,
            -d_sublead,
            -d_sublead,
            d_leading,
            -d_sublead,
            d_leading,
        ],
    ]
    ggg = _metadata("gg_ggg")
    _require_matrix_close(_gram(_projectors(ggg)), _gram(ggg_expected))


# Verify every Durham sidecar is internally consistent with the manifest
def test_durham_registry_metadata_complete_finite():
    manifest = json.loads((MG2GRA / "processes.json").read_text())
    entries = [
        entry for entry in manifest["processes"] if entry["projection"] == "durham"
    ]
    assert entries

    matches = set()
    for entry in entries:
        metadata = durham_registry.load_process_data(CARDS_ROOT, entry)
        signature = tuple(metadata["final_pdgs"])
        assert signature not in matches
        matches.add(signature)
        assert metadata["rank"] > 0
        assert metadata["ncolor"] >= metadata["rank"]
        assert math.isfinite(metadata["factorization_residual"])
        assert metadata["factorization_residual"] < 1.0e-11
        assert len(metadata["final_color_representations"]) == len(signature)
        assert all(len(candidate) == len(signature) for candidate in metadata["flow_candidates"])
        if metadata["flow_weights"]:
            for column in range(metadata["ncolor"]):
                total = sum(weights[column] for weights in metadata["flow_weights"])
                assert total == pytest.approx(1.0, abs=1.0e-14)
        else:
            assert all(
                abs(representation) == 1
                for representation in metadata["final_color_representations"]
            )


# Reject type, dimension and non-finite mutations in Durham sidecars
@pytest.mark.parametrize(
    "mutation",
    (
        "integer_type",
        "projector_nan",
        "projector_width",
        "weight_inf",
        "weight_negative",
        "weight_partition",
    ),
)
def test_durham_registry_invalid_sidecar_values(tmp_path, mutation):
    metadata = _metadata("gg_gg")
    if mutation == "integer_type":
        metadata["rank"] = True
    elif mutation == "projector_nan":
        metadata["projectors"][0][0][0] = math.nan
    elif mutation == "projector_width":
        metadata["projectors"][0].pop()
    elif mutation == "weight_inf":
        metadata["flow_weights"][0][0] = math.inf
    elif mutation == "weight_negative":
        metadata["flow_weights"][0][0] = -1.0
    else:
        metadata["flow_weights"][0][0] = 0.5

    path = tmp_path / output_layout.color_path("durham", "gg_gg")
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(metadata))
    entry = {
        "name": "gg_gg",
        "process": metadata["process"],
        "projection": "durham",
    }

    with pytest.raises(RuntimeError, match="Invalid|Non-unit"):
        durham_registry.load_process_data(tmp_path, entry)


# Verify equal stable leaves remain valid when their decay topologies differ
def test_durham_decay_signature():
    direct = (
        {"name": "direct", "process": "g g > u u~"},
        {"final_pdgs": [2, -2]},
    )
    cascade = (
        {"name": "cascade", "process": "g g > z, z > u u~"},
        {"final_pdgs": [2, -2]},
    )
    durham_registry.validate_signatures([direct, cascade])


# Verify matrix_elements retain distinct topologies with identical stable leaves
def test_durham_nested_factories(tmp_path):
    entries = [
        {"name": "direct", "process": "g g > u u~", "projection": "durham"},
        {
            "name": "cascade",
            "process": "g g > z, z > u u~",
            "projection": "durham",
        },
    ]
    for entry in entries:
        path = tmp_path / output_layout.color_path("durham", entry["name"])
        path.parent.mkdir(parents=True)
        metadata = {
            "version": 1,
            "process": entry["process"],
            "incoming_pdgs": [21, 21],
            "final_pdgs": [2, -2],
            "final_color_representations": [3, -3],
            "ncolor": 1,
            "rank": 1,
            "basis_sha256": "0" * 64,
            "factorization_residual": 0.0,
            "projectors": [[[1.0, 0.0]]],
            "flow_candidates": [[[501, 0], [0, 501]]],
            "flow_weights": [[1.0]],
        }
        path.write_text(json.dumps(metadata))

    header, source = durham_registry.generate_registry(tmp_path, {"processes": entries})
    assert 'ProcessMatches(process, "direct", {2, -2})' in source
    assert 'ProcessMatches(process, "cascade", {2, -2})' in source
    assert "AssignDurhamColorFlowCandidate" in header
    assert "CollectStableLeaves" in source


# Verify the generated four-gluon process has the expected finite-Nc dimensions
def test_gggg_rank_and_color_flows():
    metadata = _metadata("gg_gggg")

    assert metadata["final_pdgs"] == [21, 21, 21, 21]
    assert metadata["final_color_representations"] == [8, 8, 8, 8]
    assert metadata["ncolor"] == 120
    assert metadata["rank"] == 8
    assert len(metadata["projectors"]) == 8
    assert len(metadata["flow_candidates"]) == 9
    assert metadata["factorization_residual"] < 1.0e-12


# Verify the generated source owns every marked process through one interface
def test_durham_registry_generation_deterministic():
    manifest = json.loads((MG2GRA / "processes.json").read_text())
    first = durham_registry.generate_registry(CARDS_ROOT, manifest)
    second = durham_registry.generate_registry(CARDS_ROOT, manifest)

    assert first == second
    for entry in manifest["processes"]:
        if entry["projection"] == "durham":
            assert f"AMP_MG5_{entry['name']}" in first[1]


