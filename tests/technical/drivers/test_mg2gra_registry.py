# Tests for the generated MadGraph amplitude process registry
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import sys
from pathlib import Path
from types import SimpleNamespace
from typing import Any

import pytest

ROOT = Path(__file__).resolve().parents[3]
MODULE_DIR = ROOT / "develop/MG2GRA"
sys.path.insert(0, str(MODULE_DIR))

import regenerate  # noqa: E402
from modules import mg5_family, output_layout, parton_registry, process_registry  # noqa: E402


# Load the checked-in converter manifest and generated process entries
def registry_datas() -> tuple[dict[str, Any], list[dict[str, Any]]]:
    manifest = json.loads((MODULE_DIR / "processes.json").read_text())
    registry = process_registry.build_registry(regenerate.CARDS_ROOT, manifest)
    return manifest, registry["processes"]


# Compute one generated process by its stable process name
def find_data(entries: list[dict[str, Any]], process_name: str) -> dict[str, Any]:
    matches = [data for data in entries if data["process_name"] == process_name]
    assert len(matches) == 1
    return matches[0]


# Build one synthetic process directly from parsed MG5 syntax
def parsed_process(process_syntax: str) -> dict[str, Any]:
    parsed = process_registry.parse_process_syntax(process_syntax)
    topology = process_registry.numeric_topology(
        parsed, process_registry.standard_particle_selectors()
    )
    return process_registry.process_data(
        "test",
        "test_process",
        process_syntax,
        parsed,
        topology,
        allows_isolated_resonance=not parsed.has_decay_chain,
    )


# Write one strict colorless family channel file
def write_test_family_channels(
    base_dir: Path, family: dict[str, Any], model: dict[str, Any]
) -> None:
    family_name = family["name"]
    channel_path = base_dir / output_layout.channels_path(family["projection"], family_name)
    channel_path.parent.mkdir(parents=True, exist_ok=True)
    family_channels = {
        "schema": 1,
        "family": family_name,
        "generation_setup": process_registry.family_generation_setup_data(
            process_registry.family_generation_setup_from_manifest(family, model)
        ),
        "mg5_source": {"mg5_version": "2.9.27", "source_digest": "0" * 64},
        "subprocesses": [
            {
                "generated_name": "P1_Sigma_sm_aa_epem",
                "channels": [
                    {
                        "initial": [22, 22],
                        "final": [-11, 11],
                        "topology": [
                            {"pdg": -11, "daughters": []},
                            {"pdg": 11, "daughters": []},
                        ],
                        "external_color_representations": [1, 1, 1, 1],
                        "external_color_flows": [[0, 0, 0, 0, 0, 0, 0, 0]],
                    }
                ],
            }
        ],
    }
    channel_path.write_text(json.dumps(family_channels))


# Convert one numeric topology to a stable JSON comparison key
def topology_key(topology: list[dict[str, Any]]) -> str:
    return json.dumps(topology, sort_keys=True, separators=(",", ":"))


# Check MG5 sigmaHat multiplicities match the exact generated channel count
def test_channel_multiplicities():
    source = """double Sigma::sigmaHat()
{
  if(id1 == 2 && id2 == -2)
  {
    return matrix_element[0] * 2;
  }
  else if(id1 == 21 && id2 == 21)
  {
    return matrix_element[1];
  }
  return 0.;
}
"""
    channels = (
        mg5_family.Channel((2, -2), (-11, 11), ()),
        mg5_family.Channel((2, -2), (-13, 13), ()),
        mg5_family.Channel((21, 21), (21, 21), ()),
    )

    assert mg5_family.selected_process_multiplicities(source, "Sigma") == {
        (2, -2): 2,
        (21, 21): 1,
    }
    mg5_family.validate_selected_multiplicities(source, "Sigma", channels)
    with pytest.raises(RuntimeError, match="multiplicities and exact channels disagree"):
        mg5_family.validate_selected_multiplicities(source, "Sigma", channels[:-1])


# Check family wrappers are derived uniquely without redundant manifest fields
@pytest.mark.parametrize(
    ("family_name", "wrapper"),
    (
        ("MG5_TEST", "AMP_MG5_test"),
        ("MG5_YY_WW", "AMP_MG5_yy_ww"),
        ("MG5_BSM_X2", "AMP_MG5_bsm_x2"),
    ),
)
def test_family_wrapper_name_is_derived(family_name: str, wrapper: str):
    assert process_registry.family_wrapper_name({"name": family_name}) == wrapper


# Check generated families select exactly one supported physical projection
def test_family_projection_is_strict():
    family = {"name": "MG5_TEST"}
    with pytest.raises(RuntimeError, match="requires projection"):
        process_registry.family_projection(family)

    family["projection"] = "durham"
    with pytest.raises(RuntimeError, match="requires projection"):
        process_registry.family_projection(family)

    family["projection"] = "parton"
    assert process_registry.family_projection(family) == "parton"
    family["projection"] = "photon"
    assert process_registry.family_projection(family) == "photon"


# Check hard parton amplitudes are accepted only in generated form
def test_hard_parton_amp_fully_generated(tmp_path: Path):
    family = {
        "name": "MG5_TEST",
        "model": "sm",
        "definitions": [],
        "processes": ["p p > mu+ mu-"],
        "projection": "parton",
        "channel": "test",
    }
    include_dir = tmp_path / output_layout.CPP_INCLUDE_ROOT / output_layout.PARTON
    source_dir = tmp_path / output_layout.CPP_SOURCE_ROOT / output_layout.PARTON
    include_dir.mkdir(parents=True)
    source_dir.mkdir(parents=True)
    wrapper = process_registry.family_wrapper_name(family)
    header = include_dir / f"{wrapper}.h"
    source = source_dir / f"{wrapper}.cc"
    header.write_text(parton_registry.amplitude_header(family))
    source.write_text(parton_registry.amplitude_source(family))
    regenerate.validate_family_amplitude(family, tmp_path)

    header.write_text(
        parton_registry.amplitude_header(family).replace(
            "public PartonMG5Process", "public amplitude::ProcessRegistry"
        )
    )
    with pytest.raises(RuntimeError, match="header is stale"):
        regenerate.validate_family_amplitude(family, tmp_path)


# Check a CLI family addition feeds the generated hard process factory
def test_parton_registry_accepts_cli_family_addition(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    manifest = {"models": {"sm": {"import": "sm"}}, "processes": [], "families": []}
    manifest_path = tmp_path / "processes.json"
    manifest_path.write_text(regenerate.dump_manifest(manifest))
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "regenerate.py",
            "--manifest",
            str(manifest_path),
            "--add-family",
            "MG5_TEST",
            "--family-process",
            "p p > mu+ mu-",
            "--projection",
            "parton",
            "--channel",
            "test",
            "--mass-override",
            "5=0",
        ],
    )

    candidate = regenerate.add_registry_entry(regenerate.parse_args(), manifest)
    _, source = parton_registry.generate_registry(candidate)
    assert candidate["families"][0]["projection"] == "parton"
    assert process_registry.family_wrapper_name(candidate["families"][0]) == ("AMP_MG5_test")
    assert candidate["families"][0]["mass_overrides"] == {"5": 0.0}
    parsed = [process_registry.parse_process_syntax("p p > mu+ mu-")]
    assert process_registry.family_process_names(candidate["families"][0], parsed)[0].startswith(
        "test_mupmum_topology_"
    )
    assert '#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_test.h"' in source
    assert 'if (process_family == "MG5_TEST")' in source
    assert "return std::make_unique<AMP_MG5_test>();" in source
    assert '"MG5cards/Parton/MG5_TEST/param_card.dat"' in source


# Check family removal leaves independently model-owned support untouched
def test_cli_family_removal_model_definition():
    manifest = {
        "models": {"sm": {"import": "sm"}},
        "processes": [],
        "families": [
            {
                "name": "MG5_FIRST",
                "model": "sm",
                "definitions": [],
                "processes": ["a a > e+ e-"],
                "projection": "photon",
                "channel": "test",
            }
        ],
    }
    candidate = json.loads(json.dumps(manifest))
    regenerate.remove_registry_entry(
        SimpleNamespace(remove=None, remove_family="MG5_FIRST"), candidate
    )

    assert candidate["families"] == []
    assert candidate["models"]["sm"] == {"import": "sm"}


# Check retirement touches family-owned outputs but not shared model support
def test_family_retirement_shared_model_support(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    base_dir = tmp_path / "develop/MG2GRA"
    cards_dir = tmp_path / output_layout.CARDS_ROOT
    include_dir = tmp_path / output_layout.CPP_INCLUDE_ROOT
    source_dir = tmp_path / output_layout.CPP_SOURCE_ROOT
    family_include = include_dir / output_layout.PARTON / "MG5_TEST"
    family_source = source_dir / output_layout.PARTON / "MG5_TEST"
    wrapper_include = include_dir / output_layout.PARTON / "AMP_MG5_test.h"
    wrapper_source = source_dir / output_layout.PARTON / "AMP_MG5_test.cc"
    model_include = include_dir / output_layout.RUNTIME / "Models/sm"
    family_cards = cards_dir / output_layout.PARTON / "MG5_TEST"
    family_include.mkdir(parents=True)
    family_source.mkdir(parents=True)
    wrapper_include.parent.mkdir(parents=True, exist_ok=True)
    wrapper_source.parent.mkdir(parents=True, exist_ok=True)
    model_include.mkdir(parents=True)
    family_cards.mkdir(parents=True)
    base_dir.mkdir(parents=True)
    (family_include / "Processes.h").write_text("family\n")
    (family_source / "Processes.cc").write_text("family\n")
    wrapper_include.write_text("wrapper\n")
    wrapper_source.write_text("wrapper\n")
    (model_include / "Parameters_sm.h").write_text("shared parameters\n")
    (model_include / "HelAmps_sm.h").write_text("shared HELAS\n")
    (family_cards / output_layout.CHANNEL_DATA).write_text("channels\n")
    (family_cards / output_layout.PARAMETER_CARD).write_text("parameters\n")

    monkeypatch.setattr(regenerate, "ROOT", tmp_path)
    monkeypatch.setattr(regenerate, "BASE_DIR", base_dir)
    monkeypatch.setattr(regenerate, "CARDS_ROOT", cards_dir)
    original = {
        "processes": [],
        "families": [
            {
                "name": "MG5_TEST",
                "model": "sm",
                "projection": "parton",
                "channel": "test",
            }
        ],
    }
    current = {"processes": [], "families": []}
    transaction = regenerate.InstallTransaction(tmp_path / "journal")
    regenerate.retire_removed_outputs(original, current, transaction)

    assert not (family_include / "Processes.h").exists()
    assert not (family_source / "Processes.cc").exists()
    assert not wrapper_include.exists()
    assert not wrapper_source.exists()
    assert not (family_cards / output_layout.CHANNEL_DATA).exists()
    assert not (family_cards / output_layout.PARAMETER_CARD).exists()
    assert (model_include / "Parameters_sm.h").read_text() == "shared parameters\n"
    assert (model_include / "HelAmps_sm.h").read_text() == "shared HELAS\n"


# Check compact source patterns own the exact correlated generated channel sets
def test_family_channel_topologies():
    manifest, entries = registry_datas()
    expected_counts = {
        "MG5_PP_Z": (6, [1] * 6),
        "MG5_PP_ZJ": (33, [11, 11, 11]),
        "MG5_PP_JJ": (66, [66]),
        "MG5_PP_W": (6, [1, 1, 1, 1, 1, 1]),
        "MG5_YY_JJ": (5, [5]),
        "MG5_YY_WW": (9, [1, 1, 1, 1, 1, 1, 1, 1, 1]),
        "MG5_YY_ZJJ": (15, [1] * 15),
    }
    for family_name, (channel_count, process_channel_counts) in expected_counts.items():
        family = next(entry for entry in manifest["families"] if entry["name"] == family_name)
        family_datas = [data for data in entries if data["process_family"] == family_name]
        assert [len(data["channel_topologies"]) for data in family_datas] == (
            process_channel_counts
        )
        registered = [topology for data in family_datas for topology in data["channel_topologies"]]
        generated = process_registry.load_family_channel_topologies(
            regenerate.CARDS_ROOT,
            family,
            manifest["models"][family["model"]],
        )
        assert len(registered) == channel_count
        assert {topology_key(topology) for topology in registered} == {
            topology_key(topology) for topology in generated
        }

    dijet = find_data(entries, "pp_jj")
    stable_pairs = {
        tuple(process_registry.topology_stable_pattern(topology)[0])
        for topology in dijet["channel_topologies"]
    }
    assert (21, 2) in stable_pairs
    assert (2, 21) not in stable_pairs


# Check exact standalone processes expose their one generated channel
def test_standalone_topology():
    manifest, entries = registry_datas()
    standalone_names = {entry["name"] for entry in manifest["processes"]}
    standalone = [data for data in entries if data["process_name"] in standalone_names]
    assert standalone
    for data in standalone:
        if data["matrix_element_form"] == "generated" and data["topology_mode"] == "exact":
            assert data["channel_topologies"] == [data["topology"]]


# Check stale generated process data is rejected before registry generation
def test_process_registry_rejects_stale_process_data(tmp_path: Path):
    process_dir = tmp_path / output_layout.PHOTON / "yy_test"
    process_dir.mkdir(parents=True)
    (process_dir / output_layout.COLOR_DATA).write_text(
        json.dumps({"process": "a a > mu+ mu-", "final_pdgs": [-11, 11]})
    )
    manifest = {
        "processes": [
            {
                "name": "yy_test",
                "model": "sm",
                "process": "a a > e+ e- QCD=0",
                "projection": "photon",
            }
        ],
        "families": [],
    }
    with pytest.raises(RuntimeError, match="Stale generated MG5 process data"):
        process_registry.build_registry(tmp_path, manifest)


# Check standalone process data keeps exact PDG and decay-chain JSON types
@pytest.mark.parametrize(
    ("field", "value", "message"),
    (
        ("final_pdgs", [True, 11], "stable final-state PDGs"),
        ("final_pdgs", [-11, 11.0], "stable final-state PDGs"),
        ("final_pdgs", [0, 11], "stable final-state PDGs"),
        ("final_pdgs", [], "stable final-state PDGs"),
        ("has_decay_chain", 0, "decay-chain flag"),
        ("has_decay_chain", "false", "decay-chain flag"),
    ),
)
def test_malformed_standalone_data(
    tmp_path: Path, field: str, value: Any, message: str
):
    process_dir = tmp_path / output_layout.PHOTON / "yy_test"
    process_dir.mkdir(parents=True)
    process_data = {
        "process": "a a > e+ e- QCD=0",
        "final_pdgs": [-11, 11],
        "has_decay_chain": False,
    }
    process_data[field] = value
    (process_dir / output_layout.COLOR_DATA).write_text(json.dumps(process_data))
    entry = {
        "name": "yy_test",
        "model": "sm",
        "process": "a a > e+ e- QCD=0",
        "projection": "photon",
    }

    with pytest.raises(RuntimeError, match=message):
        process_registry.load_standalone_process_data(tmp_path, entry)


# Check every registered WW channel is a complete stable-leaf amplitude
def test_ww_topology_and_decay_structures():
    _, entries = registry_datas()
    ww = [entry for entry in entries if entry["process_family"] == "MG5_YY_WW"]

    assert len(ww) == 9
    assert {tuple(entry["stable_pdgs"]) for entry in ww} == {
        (-plus, plus + 1, minus, -(minus + 1)) for plus in (11, 13, 15) for minus in (11, 13, 15)
    }
    for entry in ww:
        assert entry["topology_signature"].startswith("w+(")
        assert entry["decay_structure"] == {
            "type": "full",
            "allows_isolated_resonance": False,
        }


# Check C++ generation preserves the complete generated decay amplitude
def test_cpp_decay_structure_decay_amps():
    process = parsed_process("a a > e+ e-")
    process_registry.validate_registry([process])
    generated = process_registry.cpp_process(process)
    assert generated.startswith("      Process{")
    assert "DecayType::Full, true" in generated


# Reject proposal declarations that omit or duplicate generated cascade propagators
@pytest.mark.parametrize(
    "syntax",
    (
        "a a > e+ e-",
        "p p > z j, z > mu+ mu-",
        "p p > t t~, (t > w+ b, w+ > e+ ve), (t~ > w- b~, w- > e- ve~)",
        "a a > h, (h > z z, (z > e+ e-), (z > mu+ mu-))",
    ),
)
@pytest.mark.parametrize("decay_type", ("incoherent", "jacob_wick"))
def test_decay_topology_mismatch(syntax: str, decay_type: str):
    process = parsed_process(syntax)
    process_registry.validate_registry([process])
    decay = process["decay_structure"]
    decay["type"] = decay_type
    with pytest.raises(RuntimeError, match="complete decay amplitude"):
        process_registry.validate_registry([process])


# Check pure-gluon display labels remain compact while syntax remains valid
def test_durham_gluon_labels_and_syntax():
    _, entries = registry_datas()
    two_gluons = find_data(entries, "gg_gg")
    three_gluons = find_data(entries, "gg_ggg")

    assert two_gluons["stable_final_state"] == "gg"
    assert two_gluons["final_state_syntax"] == "g g"
    assert two_gluons["topology_signature"] == "g,g"
    assert three_gluons["stable_final_state"] == "ggg"
    assert three_gluons["final_state_syntax"] == "g g g"
    assert three_gluons["topology_signature"] == "g,g,g"


# Check direct Drell-Yan and explicit Z decay structures differ
def test_fixed_family_direct_decay_chain_structures():
    _, entries = registry_datas()
    direct = find_data(entries, "pp_z_mupmum")
    cascade = find_data(entries, "pp_zj_mupmum")

    assert direct["process_family"] == "MG5_PP_Z"
    assert direct["final_state_syntax"] == "mu+ mu-"
    assert direct["decay_structure"]["type"] == "full"
    assert not direct["decay_structure"]["allows_isolated_resonance"]
    assert direct["matrix_element_form"] == "generated"
    assert direct["topology_mode"] == "exact"
    assert direct["stable_pdgs"] == [-13, 13]
    assert cascade["process_family"] == "MG5_PP_ZJ"
    assert cascade["final_state_syntax"] == "Z > {mu+ mu-} j"
    assert cascade["topology_signature"] == "z(mu+,mu-),j"
    assert cascade["decay_structure"]["type"] == "full"
    assert not cascade["decay_structure"]["allows_isolated_resonance"]
    assert cascade["matrix_element_form"] == "generated"
    assert cascade["topology_mode"] == "flavour_set"
    assert cascade["stable_pdgs"] == [-13, 13, 0]
    assert cascade["topology"][1] == {
        "allowed_pdgs": [-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89],
        "daughters": [],
    }


# Check GRANIITTI tau spelling does not alter the original MG5 syntax
def test_tau_uses_graniitti_process_spelling():
    _, entries = registry_datas()
    direct = find_data(entries, "pp_z_taptam")
    cascade = find_data(entries, "pp_zj_taptam")

    assert direct["stable_final_state"] == "tau+ tau-"
    assert direct["final_state_syntax"] == "tau+ tau-"
    assert direct["process_syntax"] == "p p > ta+ ta- / h"
    assert cascade["stable_final_state"] == "tau+ tau- j"
    assert cascade["final_state_syntax"] == "Z > {tau+ tau-} j"
    assert cascade["process_syntax"] == "p p > z j, z > ta+ ta-"


# Check generated families use their projection-owned wrapper and channel data
def test_family_registry_projected_outputs(tmp_path: Path):
    family = {
        "name": "MG5_TEST_FAMILY",
        "model": "sm",
        "definitions": [],
        "processes": ["a a > e+ e-"],
        "process_names": ["test_family"],
        "projection": "photon",
        "channel": "test",
    }
    model = {"import": "sm"}
    manifest = {"models": {"sm": model}, "processes": [], "families": [family]}
    write_test_family_channels(tmp_path, family, model)
    entries = process_registry.build_registry(tmp_path, manifest)["processes"]
    assert [(data["process_family"], data["process_name"]) for data in entries] == [
        ("MG5_TEST_FAMILY", "test_family")
    ]
    assert process_registry.family_wrapper_name(family) == "AMP_MG5_test_family"
    assert output_layout.channels_path(family["projection"], family["name"]) == (
        "Photon/MG5_TEST_FAMILY/channels.json"
    )
    assert (
        output_layout.parameter_card_path(family["projection"], family["name"])
        == "Photon/MG5_TEST_FAMILY/param_card.dat"
    )
    source = process_registry.registry_source({"processes": entries}, manifest)
    assert (
        'if (process_family == "MG5_TEST_FAMILY" &&\n'
        '      process_name == "test_family") {\n'
        '    return "MG5cards/Photon/MG5_TEST_FAMILY/param_card.dat";'
        in source
    )


# Check process edits cannot reuse stale generated family channels
def test_family_stale_generation(tmp_path: Path):
    family = {
        "name": "MG5_TEST_FAMILY",
        "model": "sm",
        "definitions": [],
        "processes": ["a a > e+ e-"],
        "projection": "photon",
        "channel": "test",
    }
    model = {"import": "sm"}
    manifest = {"models": {"sm": model}, "processes": [], "families": [family]}
    write_test_family_channels(tmp_path, family, model)
    family["processes"][0] += " QCD=0"

    with pytest.raises(RuntimeError, match="generation setup"):
        process_registry.build_registry(tmp_path, manifest)


# Check an unsupported family projection fails before registry generation
def test_family_registry_invalid_projection(tmp_path: Path):
    family = {
        "name": "MG5_TEST_FAMILY",
        "model": "sm",
        "definitions": [],
        "processes": ["a a > e+ e-"],
        "projection": "durham",
        "channel": "test",
    }
    manifest = {"models": {"sm": {"import": "sm"}}, "processes": [], "families": [family]}
    with pytest.raises(RuntimeError, match="requires projection"):
        process_registry.build_registry(tmp_path, manifest)


# Check duplicate process family and process name keys are rejected
def test_process_registry_rejects_duplicate_keys():
    _, entries = registry_datas()
    duplicate = dict(entries[0])
    with pytest.raises(RuntimeError, match="Duplicate MG5 process key"):
        process_registry.validate_registry([*entries, duplicate])


# Check distinct process names cannot ambiguously claim one family topology
def test_duplicate_family_topologies():
    _, entries = registry_datas()
    duplicate = dict(entries[0])
    duplicate["process_name"] = "distinct_name_same_topology"
    with pytest.raises(RuntimeError, match="Duplicate MG5 process topology"):
        process_registry.validate_registry([*entries, duplicate])


# Check standard parenthesized top and W decay chains retain occurrence order
def test_parenthesized_nested_top_decay_chain():
    source = "p p > t t~, (t > w+ b, w+ > e+ ve), (t~ > w- b~, w- > e- ve~)"
    parsed = process_registry.parse_process_syntax(source)
    process = parsed_process(source)

    assert process_registry.stable_tokens(parsed) == (
        "e+",
        "ve",
        "b",
        "e-",
        "ve~",
        "b~",
    )
    assert process["stable_final_state"] == "e+ ve b e- ve~ b~"
    assert process["final_state_syntax"] == ("t > {W+ > {e+ ve} b} t~ > {W- > {e- ve~} b~}")
    assert process["topology_signature"] == ("t(w+(e+,ve),b),t~(w-(e-,ve~),b~)")
    assert process["stable_pdgs"] == [-11, 12, 5, 11, -12, -5]
    assert process["matrix_element_form"] == "generated"
    assert process["topology_mode"] == "exact"
    assert process["topology"][0]["allowed_pdgs"] == [6]
    assert process["topology"][0]["daughters"][0]["allowed_pdgs"] == [24]
    assert process["decay_structure"] == {
        "type": "full",
        "allows_isolated_resonance": False,
    }


# Check identical mothers receive separately scoped decays by occurrence
def test_repeated_mother_decay_scope():
    source = "a a > h, (h > z z, (z > e+ e-), (z > mu+ mu-))"
    parsed = process_registry.parse_process_syntax(source)
    process = parsed_process(source)

    assert process_registry.stable_tokens(parsed) == ("e+", "e-", "mu+", "mu-")
    assert process["stable_final_state"] == "e+ e- mu+ mu-"
    assert process["final_state_syntax"] == ("H > {Z > {e+ e-} Z > {mu+ mu-}}")
    assert process["topology_signature"] == "h(z(e+,e-),z(mu+,mu-))"
    assert process["stable_pdgs"] == [-11, 11, -13, 13]
    assert process["topology"][0]["daughters"][1]["allowed_pdgs"] == [23]
    assert process["decay_structure"]["type"] == "full"


# Check later repeated mothers bind before deeper descendants of earlier decays
def test_repeated_mothers_attach_shallowest_depth():
    source = "a a > h h, (h > h h), (h > e+ e-)"
    parsed = process_registry.parse_process_syntax(source)

    assert process_registry.topology_signature(parsed) == "h(h,h),h(e+,e-)"
    assert process_registry.stable_tokens(parsed) == ("h", "h", "e+", "e-")


# Check generated families cannot claim standalone process family names
@pytest.mark.parametrize("name", ("DURHAM", "PHOTON"))
def test_reserved_family_names(name):
    with pytest.raises(RuntimeError, match="reserved"):
        process_registry.family_process_family({"name": name})


# Check generated families remain isolated in the MG5 process namespace
def test_process_family_requires_mg5_namespace():
    with pytest.raises(RuntimeError, match="MG5_ namespace"):
        process_registry.family_process_family({"name": "CUSTOM"})


# Check explicit aliases become local ordered selector nodes with hard-jet support
def test_mixed_parton_family_aliases():
    source = "g g > g j"
    parsed = process_registry.parse_process_syntax(source)
    selectors = process_registry.standard_particle_selectors()
    aliases = process_registry.parse_family_aliases(
        [
            "define q = u d s c b u~ d~ s~ c~ b~",
            "define j = g q",
        ]
    )
    topology = process_registry.numeric_topology(parsed, selectors, aliases)
    process = process_registry.process_data(
        "test",
        "mixed_gj",
        source,
        parsed,
        topology,
        allows_isolated_resonance=False,
    )

    assert topology[0] == {"allowed_pdgs": [21], "daughters": []}
    assert topology[1]["allowed_pdgs"] == [
        -89,
        -5,
        -4,
        -3,
        -2,
        -1,
        1,
        2,
        3,
        4,
        5,
        21,
        89,
    ]
    assert process["stable_pdgs"] == [21, 0]
    assert process["matrix_element_form"] == "generated"
    assert process["topology_mode"] == "flavour_set"
    generated = process_registry.cpp_process(process)
    assert "AmplitudeTopologyNode{std::vector<int>{21}" in generated
    assert "MatrixElementForm::Generated" in generated
    assert "TopologyMode::FlavourSet" in generated


# Check model PDG overrides resolve internal BSM nodes while process data owns leaves
def test_standalone_pdg_override():
    source = "a a > x, x > e+ e-"
    parsed = process_registry.parse_process_syntax(source)
    manifest = {
        "models": {"bsm": {"particle_pdgs": {"x": 9000001}}},
    }
    selectors = process_registry.model_particle_selectors(manifest, "bsm")
    topology = process_registry.numeric_topology(parsed, selectors)
    assert topology[0]["allowed_pdgs"] == [9000001]
    assert topology[0]["daughters"][0]["allowed_pdgs"] == [-11]

    invalid_manifest = {
        "models": {"bsm": {"particle_pdgs": {"x": [9000001, 9000002]}}},
    }
    with pytest.raises(RuntimeError, match="one nonzero integer PDG code"):
        process_registry.model_particle_selectors(invalid_manifest, "bsm")

    with pytest.raises(RuntimeError, match="Unknown MG5 particle token x"):
        process_registry.numeric_topology(parsed, process_registry.standard_particle_selectors())

    unknown_leaf = process_registry.parse_process_syntax("a a > x x~")
    standalone = process_registry.numeric_topology(
        unknown_leaf,
        process_registry.standard_particle_selectors(),
        stable_leaf_pdgs=[9000001, -9000001],
    )
    assert [node["allowed_pdgs"] for node in standalone] == [
        [9000001],
        [-9000001],
    ]


# Check selector intersections allow strict fallbacks and reject ambiguity
def test_selector_overlap_position():
    selectors = process_registry.standard_particle_selectors()
    aliases = process_registry.parse_family_aliases(["define j = g u d s c b u~ d~ s~ c~ b~"])

    def data(source: str, name: str) -> dict[str, Any]:
        parsed = process_registry.parse_process_syntax(source)
        topology = process_registry.numeric_topology(parsed, selectors, aliases)
        return process_registry.process_data(
            "test", name, source, parsed, topology, allows_isolated_resonance=False
        )

    exact = data("g g > g g", "exact")
    generic = data("g g > j j", "generic")
    process_registry.validate_registry([exact, generic])

    first = data("g g > g j", "first")
    second = data("g g > j g", "second")
    with pytest.raises(RuntimeError, match="Overlapping incomparable"):
        process_registry.validate_registry([first, second])


# Check topology-derived family names are independent of sibling membership
def test_family_process_names_topology_fallback():
    family = {"name": "MG5_COLLIDE"}
    rows = [
        process_registry.parse_process_syntax("a a > h, h > e+ e-"),
        process_registry.parse_process_syntax("a a > z, z > e+ e-"),
    ]
    first = process_registry.family_process_names(family, rows)
    second = process_registry.family_process_names(family, [rows[1]])
    assert first[1] == second[0]
    assert len(set(first)) == 2
    assert all(name.startswith("collide_epem_topology_") for name in first)


# Check explicit process identifiers remain stable when family membership changes
def test_family_process_names_prefer_explicit_ids():
    family = {"name": "MG5_TEST", "process_names": ["electron", "muon"]}
    rows = [
        process_registry.parse_process_syntax("a a > e+ e-"),
        process_registry.parse_process_syntax("a a > mu+ mu-"),
    ]
    assert process_registry.family_process_names(family, rows) == ["electron", "muon"]


# Check generated process data rejects topology inconsistencies
def test_registry_rejects_inconsistent_topology_data():
    process = parsed_process("a a > e+ e-")

    wrong_stable = dict(process)
    wrong_stable["stable_pdgs"] = [0, 11]
    with pytest.raises(RuntimeError, match="stable_pdgs disagree"):
        process_registry.validate_registry([wrong_stable])

    wrong_match = dict(process)
    wrong_match["topology_mode"] = "flavour_set"
    with pytest.raises(RuntimeError, match="topology mode disagrees"):
        process_registry.validate_registry([wrong_match])

    wrong_selector = json.loads(json.dumps(process))
    wrong_selector["topology"][0]["allowed_pdgs"] = [-11, -13]
    with pytest.raises(RuntimeError, match="sorted, unique and nonzero"):
        process_registry.validate_registry([wrong_selector])

    unresolved_generated = json.loads(json.dumps(process))
    unresolved_generated["channel_topologies"] = [
        [
            {"allowed_pdgs": [-13, -11], "daughters": []},
            {"allowed_pdgs": [11], "daughters": []},
        ]
    ]
    with pytest.raises(RuntimeError, match="generated topology 0 is not exact"):
        process_registry.validate_registry([unresolved_generated])

    outside_generated = json.loads(json.dumps(process))
    outside_generated["channel_topologies"] = [
        [
            {"allowed_pdgs": [-13], "daughters": []},
            {"allowed_pdgs": [13], "daughters": []},
        ]
    ]
    with pytest.raises(RuntimeError, match="outside its selector pattern"):
        process_registry.validate_registry([outside_generated])

    duplicate_generated = json.loads(json.dumps(process))
    duplicate_generated["channel_topologies"] = [
        duplicate_generated["topology"],
        duplicate_generated["topology"],
    ]
    with pytest.raises(RuntimeError, match="duplicate generated topology"):
        process_registry.validate_registry([duplicate_generated])


# Check malformed and unmatched parenthesized groups fail before registration
@pytest.mark.parametrize(
    "source, message",
    [
        ("p p > t t~, (t > w+ b", "Unmatched '\\('"),
        ("p p > t t~, t > w+ b)", "Unmatched '\\)'"),
        (
            "p p > t t~, (t > w+ b) trailing",
            "Malformed parenthesized MG5 decay group",
        ),
    ],
)
def test_parenthesized_decay_group_errors(source: str, message: str):
    with pytest.raises(RuntimeError, match=message):
        process_registry.parse_process_syntax(source)
