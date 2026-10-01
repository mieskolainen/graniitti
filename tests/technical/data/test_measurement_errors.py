# Covariance checks through every published measurement reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.io import hepdata_reader, steering
from core.stats.uncertainty import covariance_square_root, source_covariance

ROOT = Path(__file__).resolve().parents[3]


# Cover every original JSON measurement, including inactive and alternate tuning cards
def measurement_cards():
    for path in sorted((ROOT / "icepack").rglob("dataset*.json")):
        if "._old" in str(path):
            continue
        dataset, _ = steering.load_dataset(str(path), cdir=ROOT)
        histograms = [hist for entry in dataset["sets"] if entry.get("data", True) for hist in entry["hist"]]
        if histograms and dataset["datapath"].startswith("HEPData/") and all(hist["file"].endswith(".json") for hist in histograms):
            yield path


CARDS = list(measurement_cards())


# Exercise each configured table selection using the actual reader and its declared arguments
def tables(path):
    card, resolved = steering.load_dataset(str(path), cdir=ROOT)
    for entry in card["sets"]:
        if not entry.get("data", True):
            continue
        obs, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
        for histogram in entry["hist"]:
            filename = steering.resolve_data_reference(histogram["file"], datapath=card["datapath"],
                                                        dataset_path=resolved, cdir=ROOT)
            data = hepdata_reader.read(card["reader"], dataset_path=resolved, cdir=str(ROOT), filename=filename,
                file_filter=histogram.get("file_filter"), rebin_factor=histogram.get("rebin_factor"),
                hist=histogram, dataset=entry, all_obs=obs, obs=histogram["obs"])
            if isinstance(data["uncertainties"], dict):
                for index in data["uncertainties"]:
                    local = {key: value[index] for key, value in data.items() if isinstance(value, dict) and index in value}
                    yield filename, local
            else:
                yield filename, data


# All measured errors must remain finite and yield a positive semidefinite covariance
@pytest.mark.parametrize("path", CARDS, ids=lambda path: str(path.parent.relative_to(ROOT / "icepack")))
def test_measurement_covariance(path):
    count = 0
    for filename, data in tables(path):
        count += 1
        assert data["uncertainties"], filename
        matrix = sum(source_covariance(source) for source in data["uncertainties"])
        assert np.isfinite(matrix).all(), filename
        np.testing.assert_allclose(matrix, matrix.T, err_msg=filename)
        covariance_square_root(matrix)
        np.testing.assert_allclose(data["y_err"]**2, np.diag(matrix), err_msg=filename)
        for side in ("up", "down"):
            error = data.get(f"y_err_{side}", data["y_err"])
            assert np.isfinite(error).all(), filename
            assert np.all(error >= 0), filename
    assert count


# Published CMS covariance off-diagonals must survive point sorting and row selection
@pytest.mark.parametrize("paper,table", [("2648536-v2", 2), ("2899343-v1", 3)])
def test_cms_covariance_selection(paper, table):
    from icepack.UPC._common import upc_reader as reader
    path = ROOT / f"HEPData/UPC/HEPData-ins{paper}-json/Table{table}.json"
    selection = {"covariance": "Totalcovariancematrix.json"}
    whole = reader.read(str(path), hist=selection)
    selected = reader.read(str(path), hist={**selection, "rows": [4, 1]})
    matrix = source_covariance(whole["uncertainties"][0])
    assert np.any(matrix - np.diag(np.diag(matrix)))
    valid = np.flatnonzero(selected["valid"])
    covariance = source_covariance(selected["uncertainties"][0])
    np.testing.assert_allclose(covariance[np.ix_(valid, valid)], matrix[np.ix_([1, 4], [1, 4])])
    np.testing.assert_array_equal(selected["valid"], [True, False, True])
    np.testing.assert_allclose(covariance[~selected["valid"]], 0.0)
    np.testing.assert_allclose(selected["binedges"][valid], whole["binedges"][[1, 4]])


# Select cards with actual histogram values per bin, excluding derived point measurements
def bin_cards():
    for path in CARDS:
        card, resolved = steering.load_dataset(str(path), cdir=ROOT)
        for entry in card["sets"]:
            selected = [item for item in entry["hist"] if item.get("differential") is False]
            if not entry.get("data", True) or not selected:
                continue
            obs, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
            if any(obs[item["obs"]].get("kind", "histogram") == "histogram" for item in selected):
                yield path
                break


# Check every configured cross section or count per bin through the shared comparison path
@pytest.mark.parametrize("path", list(bin_cards()), ids=lambda path: str(path.parent.relative_to(ROOT / "icepack")))
def test_bin_cross_section_normalization(path):
    from core import iceplot
    from core.io import readers
    from core.plot import plot

    card, resolved = steering.load_dataset(str(path), cdir=ROOT)
    density = card["plot"]["normalization"] == "unit_density"
    for entry in card["sets"]:
        selected = [item for item in entry["hist"] if item.get("differential") is False]
        if not entry.get("data", True) or not selected:
            continue
        definitions, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
        selected = [item for item in selected if definitions[item["obs"]].get("kind", "histogram") == "histogram"]
        if not selected:
            continue
        data, obs = readers.read_hepdata({**entry, "hist": selected}, card["datapath"], card["type"], definitions,
            cdir=str(ROOT), reader=card["reader"], dataset_path=resolved)
        measured = plot.histhepdata(data, obs, density=density)
        for item in selected:
            name = item["obs"]
            values = data[name]
            valid = values.get("valid", np.ones(len(values["y"]), dtype=bool))
            expected = values["y"] * item["scale"] * 1e12
            scale = entry.get("mc_scale", 1.0)
            weights = expected / (scale * np.asarray(obs[name].get("mc_scale", 1.0)))
            events = dict(data={name: values["x"][valid]}, weights=weights[valid], xsection_pb=weights[valid].sum())
            mc = plot.histmc(events, {name: obs[name]}, scale=scale, density=density)[name]["hdata"]
            reference = measured[name]
            np.testing.assert_allclose(mc.counts_scaled[valid], reference["hdata"].counts_scaled[valid])
            integral = 1.0 if density else expected[valid].sum()
            assert iceplot.histogram_integral(reference["hdata"]) == pytest.approx(integral)
            assert iceplot.histogram_integral(mc) == pytest.approx(integral)
            if not density:
                covariance = sum(source_covariance(source) for source in values["uncertainties"])
                error = np.sqrt(valid @ covariance @ valid) * item["scale"] * 1e12
                assert iceplot.histogram_integral_error(reference["hdata"], reference["uncertainties"]) == pytest.approx(error)
