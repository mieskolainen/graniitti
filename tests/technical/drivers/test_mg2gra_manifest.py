# Tests for the projection-owned MadGraph manifest and card layout
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import copy
import sys
from pathlib import Path
from types import SimpleNamespace

import pytest

ROOT = Path(__file__).resolve().parents[3]
MG2GRA = ROOT / "develop/MG2GRA"
sys.path.insert(0, str(MG2GRA))

import regenerate  # noqa: E402
from modules import process_registry  # noqa: E402


# Load the strict checked-in process manifest
def manifest() -> dict:
    return regenerate.load_manifest(MG2GRA / "processes.json")


# Derive incoming PDGs and public wrapper names from physical manifest data
def test_projection_incoming_wrapper():
    data = manifest()
    incoming = {
        entry["projection"]: process_registry.standalone_incoming_pdg(entry)
        for entry in data["processes"]
    }
    assert incoming == {"durham": 21, "photon": 22}
    wrappers = {
        family["name"]: process_registry.family_wrapper_name(family)
        for family in data["families"]
    }
    assert wrappers["MG5_PP_Z"] == "AMP_MG5_pp_z"
    assert wrappers["MG5_YY_WW"] == "AMP_MG5_yy_ww"


# Reject a projection which disagrees with the MG5 incoming state
def test_projection_incoming_mismatch():
    data = copy.deepcopy(manifest())
    data["processes"][0]["projection"] = "photon"
    with pytest.raises(RuntimeError, match="requires incoming a a"):
        regenerate.validate_manifest(data)


# Add a family using only projection-owned public amplitude inputs
def test_cli_family_add_projection_schema():
    data = {"models": {"sm": {"import": "sm"}}, "processes": [], "families": []}
    args = SimpleNamespace(
        add=None,
        add_family="MG5_TEST",
        remove=None,
        remove_family=None,
        family_process=["p p > mu+ mu-"],
        definition=[],
        projection="parton",
        channel="test",
        model="sm",
        model_import=None,
        complex_mass_scheme=False,
        alpha_charge=None,
        alpha_charge_square=None,
        inverse_alpha_zero=137.03599908,
        mass_override=[],
        mg5_process=None,
    )
    candidate = regenerate.add_registry_entry(args, data)
    assert candidate["families"] == [
        {
            "name": "MG5_TEST",
            "model": "sm",
            "definitions": [],
            "processes": ["p p > mu+ mu-"],
            "projection": "parton",
            "channel": "test",
        }
    ]
