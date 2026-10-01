# Regression checks for measured intervals and weighted histogram bins
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from pathlib import Path

import numpy as np
import pytest
import torch
from core.io.hepdata_reader import load_dataset_reader
from core.stats import hist
from core.stats.hist import doublet2linear
from core.stats.uncertainty import source_covariance

from icepack._common.hepdata import rebin_common, table_bins

ROOT = Path(__file__).resolve().parents[3]


# Keep physical gaps under changes of units and reject overlapping intervals
@pytest.mark.parametrize("scale", [1e-12, 1.0, 1e12])
def test_interval_edges(scale):
    np.testing.assert_allclose(
        doublet2linear(np.array([[0., 1.], [2., 3.]]) * scale) / scale, [0., 1., 2., 3.]
    )
    np.testing.assert_allclose(
        doublet2linear(np.array([[0., 1.], [np.nextafter(1., 2.), 2.]]) * scale) / scale, [0., 1., 2.]
    )
    np.testing.assert_array_equal(
        doublet2linear(np.array([[0., 1e-9], [2e-9, 1e10]]) * scale), np.array([0., 1e-9, 2e-9, 1e10]) * scale
    )
    with pytest.raises(ValueError, match="overlap|ordered"):
        doublet2linear(np.array([[0., 2e-9], [1e-9, 1e10]]) * scale)
    with pytest.raises(ValueError, match="overlap|ordered"):
        doublet2linear(np.array([[0., 2.], [1., 3.]]) * scale)


# Reject malformed intervals before constructing physical histogram bins
@pytest.mark.parametrize("edges", [[], [[1., 0.]], [[0., np.inf]], [[np.nan, 1.]], [[0., 0.]], [[1., 2.], [0., 1.]]])
def test_invalid_interval_edges(edges):
    with pytest.raises(ValueError):
        doublet2linear(edges)


# Preserve measured STAR yields and all source covariances across acceptance gaps
@pytest.mark.parametrize("factor", [1, 2, 3, 4, 99])
def test_star_rebin_gaps(factor):
    reader = load_dataset_reader("icepack/SOFTCEP/STAR_1792394/KK/dataset.json", cdir=str(ROOT))
    filename = ROOT / "HEPData/CEP/HEPData-ins1792394-v1-json/Figure12(middle,right).json"
    original = reader.read(str(filename))
    data = reader.read(str(filename), rebin_factor=factor)
    old_width = original["binwidth"] * original["valid"]
    new_width = data["binwidth"] * data["valid"]
    assert new_width.sum() == pytest.approx(old_width.sum())
    assert new_width @ data["y"] == pytest.approx(old_width @ original["y"])
    np.testing.assert_allclose(data["binedges"], np.column_stack((data["bins"][:-1], data["bins"][1:])))
    for old_source, new_source in zip(original["uncertainties"], data["uncertainties"], strict=True):
        assert new_width @ source_covariance(new_source) @ new_width == pytest.approx(
            old_width @ source_covariance(old_source) @ old_width
        )
    groups = [np.flatnonzero(np.isclose(original["bins"], edge))[0] for edge in data["bins"]]
    for index, (start, stop) in enumerate(zip(groups[:-1], groups[1:], strict=True)):
        assert np.all(original["valid"][start:stop] == data["valid"][index])


# Keep published representative coordinates for bins that are not merged
@pytest.mark.parametrize("factor", [1, 2])
def test_rebin_representative_points(factor):
    data = {"bins": np.array([0., 1., 3., 6.]), "x": np.array([.4, 1.8, 4.2]),
            "binwidth": np.array([1., 2., 3.]), "y": np.array([10., 4., 2.]), "y_err": np.ones(3)}
    rebinned = rebin_common(copy.deepcopy(data), factor)
    if factor == 1:
        np.testing.assert_array_equal(rebinned["x"], data["x"])
    else:
        np.testing.assert_allclose(rebinned["x"], [1.5, 4.2])


# Keep small weighted bins independent of large weights in other bins on both backends
@pytest.mark.parametrize("backend", ["numpy", "torch"])
def test_weighted_bins_dynamic_range(backend):
    weights = np.array([1e20, 1., -2.])
    if backend == "torch":
        weights = torch.tensor(weights, requires_grad=True)
    counts, errors, _, _ = hist.hist(np.array([.5, 1.5, 2.5]), np.arange(4.), weights=weights)
    if backend == "torch":
        counts[1].backward()
        np.testing.assert_array_equal(weights.grad.numpy(), [0., 1., 0.])
        counts, errors = counts.detach().numpy(), errors.detach().numpy()
    np.testing.assert_allclose(counts, [1e20, 1., -2.])
    np.testing.assert_allclose(errors, [1e20, 1., 2.])


# Match NumPy edge semantics including underflow, overflow and the final upper edge
@pytest.mark.parametrize("backend", ["numpy", "torch"])
def test_histogram_boundaries(backend):
    edges = np.array([0., .1, .3, 1.])
    points = np.array([-np.inf, -1., 0., .1, np.nextafter(.3, 0.), .3, 1., np.nextafter(1., 2.), np.inf, np.nan])
    weights = np.arange(1., len(points) + 1)
    expected = np.histogram(points, bins=edges, weights=weights)[0]
    errors = np.sqrt(np.histogram(points, bins=edges, weights=weights**2)[0])
    indices = hist.bin_indices(points, edges)
    np.testing.assert_array_equal(indices, [-1, -1, 0, 1, 1, 2, 2, -1, -1, -1])
    if backend == "torch":
        weights = torch.tensor(weights)
    counts, observed_errors, _, _ = hist.hist(points, edges, weights=weights)
    np.testing.assert_allclose(np.asarray(counts), expected)
    np.testing.assert_allclose(np.asarray(observed_errors), errors)


# Reject shuffled centers that would otherwise produce misleading positive-width bins
@pytest.mark.parametrize("centers", [[1., 3., 2., 4.], [1., 1., 2.], [1., np.inf], [np.nan, 2.], np.array([1, 3, 2, 4], dtype=np.uint64)])
def test_invalid_centers(centers):
    with pytest.raises(ValueError, match="finite and strictly increasing"):
        hist.center2edgebins(np.asarray(centers))


# Reject incomplete, ambiguous and inconsistent published bin coordinates
@pytest.mark.parametrize("coordinates", [
    [], [{"value": "1"}], [{"value": "1", "low": "0"}],
    [{"low": "0", "high": "1", "value": "2"}],
    [{"low": "0", "high": "1"}, {"value": "2"}],
])
def test_invalid_table_bins(coordinates):
    with pytest.raises(ValueError):
        table_bins({"values": [{"x": [coordinate]} for coordinate in coordinates]})


# Preserve original UPC point intervals and covariance ordering under row selections
@pytest.mark.parametrize("filename, options", [
    ("Table1.json", {"file_filter": "0n0n"}),
    ("Table2.json", {"hist": {"covariance": "Totalcovariancematrix.json"}}),
    ("Table3.json", {}),
])
@pytest.mark.parametrize("rows", [[0], [0, 1], [-1], [-1, -2], [-1, 0], [-1, 0, 1]])
def test_upc_point_selection(filename, options, rows):
    reader = load_dataset_reader("icepack/UPC/PHOTOPROD/CMS_2648536/jpsi_coherent/dataset.json", cdir=str(ROOT))
    path = ROOT / "HEPData/UPC/HEPData-ins2648536-v2-json" / filename
    original = reader.read(str(path), **options)
    rows = np.arange(len(original["y"]))[rows].tolist()
    selected = reader.read(str(path), **{**options, "hist": {**options.get("hist", {}), "rows": rows}})
    valid = np.flatnonzero(selected["valid"])
    indices = np.searchsorted(original["x"], selected["x"][valid])
    np.testing.assert_allclose(selected["binedges"][valid], original["binedges"][indices])
    np.testing.assert_allclose(selected["y"][valid], original["y"][indices])
    for key in original:
        if key.startswith("y_err"):
            np.testing.assert_allclose(selected[key][valid], original[key][indices])
    for before, after in zip(original["uncertainties"], selected["uncertainties"], strict=True):
        matrix = source_covariance(after)
        np.testing.assert_allclose(matrix[np.ix_(valid, valid)], source_covariance(before)[np.ix_(indices, indices)])
        np.testing.assert_allclose(matrix[~selected["valid"]], 0.0)


# Reject fractional and Boolean row indices instead of selecting a different measurement
@pytest.mark.parametrize("rows", [[0.5], [True]])
def test_upc_invalid_row_indices(rows):
    reader = load_dataset_reader("icepack/UPC/PHOTOPROD/CMS_2648536/jpsi_coherent/dataset.json", cdir=str(ROOT))
    with pytest.raises(ValueError, match="integer table row"):
        reader.read(str(ROOT / "HEPData/UPC/HEPData-ins2648536-v2-json/Table2.json"), hist={"rows": rows})


# Keep narrow categories distinct in MC, measured histograms and event covariance assignments
@pytest.mark.parametrize("width", [1e-12, 1e-3, 1.0])
def test_category_intervals(width):
    from core.plot import plot
    from core.stats import cov

    bins = {index: {0: np.array([index, index + 1]) * width, 1: np.array([0., 1.]),
                    2: np.array([0., 1.])} for index in range(2)}
    obs = {"mass": {"bins": bins, "xlim": {index: [0., 1.] for index in bins},
                    "units": {"x": "GeV", "y": "pb"}, "symmetrized_fill": False,
                    "category_labels": ["t1", "t2"], "category_units": ["GeV^2", "GeV^2"]}}
    data = {"mass": {"bins": bins, "xlim": obs["mass"]["xlim"], "scale": 1e-12, "fitw": 1.0,
                     "x": {index: {2: np.array([.5])} for index in bins},
                     "binwidth": {index: {2: np.array([1.])} for index in bins},
                     "y": {index: np.array([float(index + 1)]) for index in bins},
                     "y_err": {index: np.array([.1]) for index in bins}}}
    events = {"data": {"mass": np.array([[width / 2., .5, .5], [3. * width / 2., .5, .5]])},
              "weights": np.array([1., 2.]), "xsection_pb": 3., "event_ids": np.arange(2)}
    mc = plot.histmc(events, obs)
    measured = plot.histhepdata(data, obs)
    assignments = cov.subset_assignments(events, obs)
    assert len(mc) == len(measured) == len(assignments) == 2
    assert mc.keys() == measured.keys() == assignments.keys()
    for index, name in enumerate(mc):
        assert mc[name]["hdata"].counts[0] == pytest.approx(index + 1.)
        assert measured[name]["hdata"].counts_scaled[0] == pytest.approx(index + 1.)
        expected = np.full(2, -1)
        expected[index] = 0
        np.testing.assert_array_equal(assignments[name][1], expected)
