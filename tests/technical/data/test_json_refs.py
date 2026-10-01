# Generic JSON reference loading and shared parameter edits
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path

import pyjson5
import pytest
from core.io import jsonref
from core.io.serialize import load_json_file, write_json_file
from core.tune.io import read_json

FIXTURES = Path(__file__).resolve().parents[2] / "data" / "json_refs"


# Retry partial documents and temporarily missing referenced files using real filesystem reads
@pytest.mark.parametrize("reference", [False, True])
def test_tuning_json_retries_partial_visibility(tmp_path, reference):
    path = tmp_path / "card.json"
    target = tmp_path / "source.json" if reference else path
    path.write_text('{"value": {"$ref": "source.json#/value"}}' if reference else '{"value":')
    with ThreadPoolExecutor(max_workers=1) as pool:
        result = pool.submit(read_json, str(path), retries=8)
        try:
            with pytest.raises(TimeoutError):
                result.result(timeout=0.01)
        finally:
            target.write_text('{"value": 3}')
        assert result.result(timeout=10) == {"value": 3}


# Read unexpanded test descriptions that intentionally contain invalid references
def cases():
    return json.loads((FIXTURES / "invalid.json").read_text())


# Exercise the same reference fixtures as the C++ library
def test_shared_fixtures():
    expected = json.loads((FIXTURES / "expected.json").read_text())
    assert load_json_file(FIXTURES / "valid.json", loader=pyjson5.load) == expected
    with ThreadPoolExecutor(max_workers=4) as pool:
        values = list(pool.map(lambda _: jsonref.load(FIXTURES / "valid.json", loader=pyjson5.load), range(8)))
    assert all(value == expected for value in values)
    values[0]["object"]["rows"][1][0] = 99
    assert values[0]["array"] == [3, 4]
    assert values[1] == expected


# Reject malformed targets and cyclic reference graphs before model initialization
@pytest.mark.parametrize("case", cases(), ids=lambda case: case["name"])
def test_invalid_references(case):
    with pytest.raises((ValueError, OSError)):
        jsonref.resolve(case["input"], FIXTURES / "input.json", loader=pyjson5.load)


# Keep shared parameter sources when editing a resolved card
@pytest.mark.parametrize("external", [False, True])
def test_shared_edits(tmp_path, external):
    source = tmp_path / "source.json"
    card = tmp_path / "card.json"
    initial = {"source": {"x": 1, "array": [2, 3], "unused": 4}}
    source.write_text(json.dumps(initial))
    raw = {"alias": {"$ref": ("source.json" if external else "") + "#/source"}}
    if not external:
        raw.update(initial)
    card.write_text(json.dumps(raw))
    payload = load_json_file(card)
    payload["alias"]["x"] = 7
    payload["alias"]["array"][0] = 8
    payload["alias"]["new"] = True
    del payload["alias"]["unused"]
    write_json_file(card, payload)
    assert json.loads(card.read_text())["alias"] == raw["alias"]
    resolved = load_json_file(card)
    assert resolved["alias"] == {"x": 7, "array": [8, 3], "new": True}
    assert load_json_file(source if external else card)["source"] == resolved["alias"]


# Reject conflicting values assigned to aliases of the same source
def test_conflicting_edits(tmp_path):
    card = tmp_path / "card.json"
    raw = {"source": 1, "alias": {"$ref": "#/source"}}
    card.write_text(json.dumps(raw))
    payload = load_json_file(card)
    payload.update(source=2, alias=3)
    with pytest.raises(ValueError, match="Conflicting"):
        write_json_file(card, payload)
    assert json.loads(card.read_text()) == raw


# Removing an alias must leave its source value intact
def test_remove_alias(tmp_path):
    card = tmp_path / "card.json"
    card.write_text('{"source": [1,2], "alias": {"$ref": "#/source"}}')
    payload = load_json_file(card)
    del payload["alias"]
    write_json_file(card, payload)
    assert load_json_file(card) == {"source": [1, 2]}


# Replacing a referenced array or object preserves its reference
def test_replace_shared_value(tmp_path):
    card = tmp_path / "card.json"
    card.write_text('{"source": [1,2], "alias": {"$ref": "#/source"}}')
    payload = load_json_file(card)
    payload["alias"] = {"new": [True, None]}
    write_json_file(card, payload)
    assert load_json_file(card) == {"source": {"new": [True, None]}, "alias": {"new": [True, None]}}
    assert json.loads(card.read_text())["alias"] == {"$ref": "#/source"}


# Permit normal output writes to an existing empty file
def test_empty_output(tmp_path):
    path = tmp_path / "empty.json"
    path.touch()
    write_json_file(path, {"x": 1})
    assert load_json_file(path) == {"x": 1}


# A copied trial cannot change a shared parameter that still points outside its tune
def test_trial_reference_sources_stay_inside_tune(tmp_path):
    source = tmp_path / "source.json"
    source.write_text('{"x": 1}')
    trial = tmp_path / "trial"
    trial.mkdir()
    card = trial / "card.json"
    raw = '{"x": {"$ref": "../source.json#/x"}}'
    card.write_text(raw)
    payload = load_json_file(card)
    payload["x"] = 2
    with pytest.raises(ValueError, match="inside the tune"):
        write_json_file(card, payload, reference_root=trial)
    assert load_json_file(source) == {"x": 1}
    assert card.read_text() == raw


# Select saved optimizer settings while resolving references relative to their JSON file
def test_json_file_fragment(tmp_path):
    (tmp_path / "ranges.json").write_text(json.dumps({"limits": [0.25, 0.75]}))
    path = tmp_path / "setup.json"
    path.write_text(json.dumps({"optimizer": {"by/channel": {"bounds": {"$ref": "ranges.json#/limits"}}}}))
    assert load_json_file(str(path) + "#/optimizer/by~1channel") == {"bounds": [0.25, 0.75]}
    with pytest.raises(ValueError):
        load_json_file(str(path) + "#/optimizer/missing")


# Quoted process keys must select the same card values as generator overrides
@pytest.mark.parametrize("field,parts", [
    ('PARAM_NSTAR.MODEL.MP["[22,P]"]', ["PARAM_NSTAR", "MODEL", "MP", "[22,P]"]),
    ("PARAM_NSTAR.MODEL.ygg['[22]']", ["PARAM_NSTAR", "MODEL", "ygg", "[22]"]),
    ('rows[1,2].value', ["rows", 1, 2, "value"]),
    ('model["a.b"][0]', ["model", "a.b", 0]),
])
def test_field_parts(field, parts):
    assert jsonref.field_parts(field) == parts


# Reject ambiguous keys and incomplete array selectors
@pytest.mark.parametrize("field", ["", ".model", "model.", "model..x", "rows[-1]", "rows[0]x", 'model["MP[RES]"', "model[MP[RES]]"])
def test_invalid_field_parts(field):
    with pytest.raises(ValueError):
        jsonref.field_parts(field)
