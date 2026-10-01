# Check raw HEPData downloads and preservation of replaced files
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import pytest

from HEPData import download

ROOT = Path(__file__).resolve().parents[3] / "HEPData"


# Use an original repository table without hardcoded measurement values
@pytest.fixture
def table():
    raw = next(ROOT.glob("PHOTOPROD/HEPData-*-json/*.json")).read_bytes()
    return raw, json.loads(raw)


# Preserve complete response bytes and every older file across repeated downloads
def test_save(table, tmp_path):
    raw, _ = table
    path = tmp_path / "tables" / "table.json"
    temporary = tmp_path / "tmp"
    assert download.save(path, raw, temporary) == "downloaded"
    assert path.read_bytes() == raw
    before = path.stat().st_mtime_ns
    assert download.save(path, raw, temporary) == "unchanged"
    assert path.stat().st_mtime_ns == before
    for count in (1, 2):
        previous = path.read_bytes()
        response = raw + b"\n" * count
        assert download.save(path, response, temporary) == "replaced"
        assert path.read_bytes() == response
        assert path.with_name(path.name + "._old" * count).read_bytes() == previous


# Reject wrong table identities, missing data and non-JSON HTTP responses
def test_check_table(table):
    raw, data = table
    download.check_table(raw, data["doi"])
    with pytest.raises(ValueError):
        download.check_table(raw, data["doi"] + "/other")
    metadata = json.dumps({key: value for key, value in data.items() if key != "values"}).encode()
    for invalid in (metadata, b"{}", b"[]", b"<html>Service unavailable</html>"):
        with pytest.raises(ValueError):
            download.check_table(invalid, data["doi"])


# Resolve category selections and match filenames through the published DOI
def test_records_and_names():
    selected = download.records(ROOT, ["PHOTOPROD"])
    assert selected and all(p.parent == ROOT / "PHOTOPROD" for p in selected)
    assert download.records(ROOT, [str(selected[0]), str(selected[0])]) == [selected[0]]
    with pytest.raises(ValueError):
        download.records(ROOT, ["not_a_record"])
    source = selected[0]
    path = next(source.glob("*.json"))
    data = json.loads(path.read_bytes())
    tables = [{"doi": data["doi"], "processed_name": "Different official name"}]
    assert download.table_names(source, tables) == [path.name]
    with pytest.raises(ValueError):
        download.table_names(source, tables * 2)


# Preserve LaTeX backslashes in published filenames when the filesystem permits them
def test_latex_names():
    paths = [path for path in ROOT.glob("*/HEPData-*-json/*.json") if "\\" in path.name]
    assert paths
    for path in paths:
        data = json.loads(path.read_bytes())
        tables = [{"doi": data["doi"], "processed_name": path.stem}]
        assert download.table_names(path.parent, tables) == [path.name]
