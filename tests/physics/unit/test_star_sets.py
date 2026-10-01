# STAR histogram coverage and grouping by published fiducial selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from collections import Counter
from html import unescape
from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering

from icepack._common.hepdata import read_table, table_values
from icepack.SOFTCEP.STAR_3075716._common import cuts
from tests.technical.support.hepmc import make_pair_event

ROOT = Path(__file__).resolve().parents[3]


# Load both STAR measurements through the actual icepack card loader
@pytest.fixture(params=["1792394", "3075716"])
def star(request):
    record = request.param
    paths = sorted((ROOT / f"icepack/SOFTCEP/STAR_{record}").glob("*/dataset.json"))
    datasets = [steering.load_dataset(str(path), cdir=ROOT) for path in paths]
    tables = {path.name: read_table(path) for path in sorted(
        (ROOT / f"HEPData/CEP/HEPData-ins{record}-v1-json").glob("Figure*.json"))}
    assert datasets and tables
    return record, datasets, tables


# Associate figure panels with their selections [REFERENCE: arXiv:2510.27482, Figures 7-13]
def regions(table):
    figure, panel = re.fullmatch(r"Figure (\d+)([a-f])", table["name"]).groups()
    figure, panel = int(figure), ord(panel) - ord("a")
    if figure == 7:
        return None, None
    if figure in (8, 9, 10):
        return None, panel % 2
    if figure == 11:
        return panel, None
    assert figure in (12, 13)
    return panel % 3, panel // 3


# Check mass and azimuth cut modules against the original 200 GeV table descriptions
def check_cut(selection, region):
    text = unescape(region)
    limits = [float(value) for value in re.findall(r"\d+(?:\.\d+)?", text)]
    if "invariant masses" in text:
        for mass in (bound * factor for bound in limits for factor in (0.8, 1.2)):
            expected = limits[0] < mass < limits[1] if len(limits) == 2 else (
                mass < limits[0] if "<" in text else mass > limits[0])
            assert bool(selection.cut_func(make_pair_event(pid=211, mass=mass))) is expected
    elif len(limits) == 1:
        for phi in (limits[0] * 0.5, limits[0] * 1.5):
            expected = phi < limits[0] if "<" in text else phi > limits[0]
            assert bool(selection.cut_func(make_pair_event(pid=211, dphi=phi))) is expected
    elif not region:
        assert selection.cut_func(make_pair_event(pid=211))


# Include every original differential table exactly once
def test_spectra_coverage(star):
    _, datasets, tables = star
    files = [hist["file"] for dataset, _ in datasets for entry in dataset["sets"] for hist in entry["hist"]]
    assert Counter(files) == Counter({name: 1 for name in tables})


# Group all spectra with identical published cuts into one set and use those cuts
def test_selection_groups(star):
    record, datasets, tables = star
    groups = set()
    for dataset, path in datasets:
        for entry in dataset["sets"]:
            keys = set()
            for hist in entry["hist"]:
                table = tables[hist["file"]]
                if record == "1792394":
                    description = table["description"].split("\n")[0]
                    region = description.split(", for ", 1)[1] if ", for " in description else ""
                else:
                    region = regions(table)
                keys.add((tuple(table["keywords"]["reactions"]), region))
            assert len(keys) == 1, entry["name"]
            key = keys.pop()
            assert key not in groups, entry["name"]
            groups.add(key)
            selection = readers.load_cut_module(steering.resolve_python_reference(
                entry["cuts"], package="core.analysis.cuts", dataset_path=path, cdir=ROOT,
            ))
            if record == "3075716":
                mass, phi = key[1]
                expected = cuts.parameters(mass=mass, phi=phi)
                assert selection.cut_param["phi"] == phi
                for field in ("M", "dPhi"):
                    np.testing.assert_allclose(np.asarray(selection.cut_param[field], dtype=float),
                                               np.asarray(expected[field], dtype=float))
            else:
                check_cut(selection, key[1])


# Read every histogram with its published bins, units, and original measured values
def test_spectra_reading(star):
    _, datasets, tables = star
    for dataset, path in datasets:
        for entry in dataset["sets"]:
            observables, _ = steering.load_observables(entry["obs"], dataset_path=path, cdir=ROOT)
            data, histograms = readers.read_hepdata(
                entry, dataset["datapath"], dataset["type"], observables, cdir=ROOT,
                reader=dataset["reader"], dataset_path=path,
            )
            assert len(data) == len(entry["hist"])
            for hist in entry["hist"]:
                table = tables[hist["file"]]
                result = data[hist["obs"]]
                valid = result["valid"]
                np.testing.assert_allclose(result["y"][valid], [float(row["value"]) for row in table_values(table)])
                edges = np.asarray([[row["x"][0][side] for side in ("low", "high")] for row in table["values"]], dtype=float)
                np.testing.assert_allclose(np.column_stack((result["bins"][:-1], result["bins"][1:]))[valid], edges)
                scale = 1e-12 if "[PB" in table["headers"][1]["name"].upper() else 1e-9
                assert result["scale"] == pytest.approx(scale)
                np.testing.assert_array_equal(histograms[hist["obs"]]["valid"], valid)
