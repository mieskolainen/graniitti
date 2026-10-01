# Checks for the experiment-specific HEPData readers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io import steering
from core.io.hepdata_reader import load_dataset_reader
from core.plot import plot
from core.stats.hist import rebin_partition
from core.stats.uncertainty import source_covariance

from tests.technical.support.hepdata import assert_original, original

ROOT = Path(__file__).resolve().parents[3]
HEPDATA_ROOT = ROOT / "HEPData" / "CEP"
ATLAS_DATA = HEPDATA_ROOT / "HEPData-ins1377585-v1-json"
CMS_DATA = HEPDATA_ROOT / "HEPData-ins2752118-v1-json"
LHCB_1277076_DATA = HEPDATA_ROOT / "HEPData-ins1277076-v1-json"
LHCB_1373746_DATA = HEPDATA_ROOT / "HEPData-ins1373746-v1-json"
LHCB_2825384_DATA = HEPDATA_ROOT / "HEPData-ins2825384-v1-json"

ATLAS_READER = load_dataset_reader(
    "icepack/GAMMA/ATLAS_1377585/mumu/dataset.json",
    cdir=str(ROOT),
)
CMS_READER = load_dataset_reader("icepack/SOFTCEP/CMS_2752118/pipi/dataset.json", cdir=str(ROOT))
LHCB_1277076_READER = load_dataset_reader(
    "icepack/PHOTOPROD/LHCb_1277076/jpsi/dataset.json",
    cdir=str(ROOT),
)
LHCB_1373746_READER = load_dataset_reader(
    "icepack/PHOTOPROD/LHCb_1373746/dataset.json",
    cdir=str(ROOT),
)
LHCB_2825384_READER = load_dataset_reader(
    "icepack/PHOTOPROD/LHCb_2825384/jpsi/dataset.json",
    cdir=str(ROOT),
)
CMS_CUTS = steering.load_python_module(
    str(ROOT / "icepack" / "SOFTCEP" / "CMS_2752118" / "_common" / "cuts.py"),
)


# Integrate one uncertainty category with its full bin covariance
def integrated_category_uncertainty(data, category):
    size = len(data["y"])
    covariance = np.zeros((size, size), dtype=float)
    for source in data["uncertainties"]:
        if source["category"] == category:
            covariance += source_covariance(source)
    widths = np.asarray(data["binwidth"], dtype=float)
    return float(np.sqrt(widths @ covariance @ widths))


# Compute one named uncertainty source from an ordinary histogram
def named_uncertainty(data, name):
    return next(source for source in data["uncertainties"] if source["name"] == name)


# Keep the last measured STAR bin and mask the intervening unmeasured interval
@pytest.mark.parametrize("rebin_factor", [None, 2])
def test_star_1792394_disconnected_mass_bin(rebin_factor):
    reader = load_dataset_reader("icepack/SOFTCEP/STAR_1792394/KK/dataset.json", cdir=str(ROOT))
    path = HEPDATA_ROOT / "HEPData-ins1792394-v1-json" / "Figure12(middle,right).json"
    rows, values, errors = original(path)
    last = rows[-1]["x"][0]
    gap_low, gap_high = float(rows[-2]["x"][0]["high"]), float(last["low"])
    data = reader.read(str(path), rebin_factor=rebin_factor)
    assert len(data["valid"]) == len(data["y"])
    gap = (data["bins"][:-1] < gap_high) & (data["bins"][1:] > gap_low)
    assert not np.any(data["valid"][gap])
    if rebin_factor is None:
        assert_original(data, path)
        assert data["valid"][-1] and not data["valid"][-2]
        assert data["x"][-1] == pytest.approx(float(last.get("value", (float(last["low"]) + float(last["high"])) / 2)))
        for side, index in (("up", 1), ("down", 0)):
            np.testing.assert_allclose(data[f"y_err_{side}"][data["valid"]], np.sqrt(sum(error[:, index]**2 for error in errors.values())))


# Check the integrated LHCb J/psi and psi(2S) measurements
def test_lhcb_1277076_integrated_xs():
    path = LHCB_1277076_DATA / "Table1.json"
    results = LHCB_1277076_READER.read_integrated(str(path))
    rows, values, errors = original(path)
    for index, row in enumerate(rows):
        state = next(state for state, name in LHCB_1277076_READER.STATE_COLUMNS.items() if f" {name} " in row["x"][0]["value"])
        assert results[state]["value"] == pytest.approx(values[index])
        assert results[state]["stat"] == pytest.approx(errors["stat"][index].mean())
        assert results[state]["syst"] == pytest.approx(errors["sys"][index].mean())


# Check both states' published rapidity bins and integrated normalizations
@pytest.mark.parametrize("state,group", [("J/psi", 0), ("psi(2S)", 1)])
def test_lhcb_1277076_rapidity_distribution(state, group):
    path = LHCB_1277076_DATA / "Table2.json"
    data = LHCB_1277076_READER.read(str(path), state=state, region="fiducial")
    assert_original(data, path, group)
    assert np.all(data["y_err"] > data["y_err_uncorrelated"])


# Read LHCb correlated totals from HEPData without inferring their decomposition
@pytest.mark.parametrize("state", ["J/psi", "psi(2S)"])
def test_lhcb_1277076_correlated_totals(state):
    filename = str(LHCB_1277076_DATA / "Table2.json")
    table = json.loads(Path(filename).read_bytes())
    group = next(entry["group"] for entry in table["qualifiers"]["RE"]
                 if f" {LHCB_1277076_READER.STATE_COLUMNS[state]} " in entry["value"])
    points = [next(value for value in row["y"] if value["group"] == group) for row in table["values"]]
    values = np.asarray([point["value"] for point in points], dtype=float)
    statistical = np.asarray([point["errors"][0]["symerror"] for point in points], dtype=float)
    data = LHCB_1277076_READER.read(filename, state=state, region="fiducial")
    expected = np.diag(statistical**2)
    for category, label in (("statistical", f"total correlated statistical uncertainty ({LHCB_1277076_READER.STATE_COLUMNS[state]})"),
                            ("systematic", "total correlated systematic uncertainty")):
        fractions = [next(error["symerror"] for error in point["errors"] if error.get("label") == f"sys,{label}")
                     for point in points]
        error = np.asarray([float(fraction.rstrip("%")) for fraction in fractions]) * values / 100.0
        assert named_uncertainty(data, f"correlated_{category}")["up"] == pytest.approx(error)
        expected += np.outer(error, error)
    assert len(data["uncertainties"]) == 3
    np.testing.assert_allclose(sum(source_covariance(source) for source in data["uncertainties"]), expected)


# Check state selection through the same dataset metadata used by iceplot
def test_lhcb_1277076_dataset_state_selection():
    data = LHCB_1277076_READER.read(
        str(LHCB_1277076_DATA / "Table2.json"),
        dataset={"state": "psi(2S)", "region": "fiducial"},
    )

    assert_original(data, LHCB_1277076_DATA / "Table2.json", group=1)


# Reject ambiguous calls instead of silently selecting one charmonium state
def test_lhcb_1277076_requires_a_state():
    with pytest.raises(ValueError, match="dataset state is required"):
        LHCB_1277076_READER.read(
            str(LHCB_1277076_DATA / "Table2.json"),
            region="fiducial",
        )


# Reject ambiguous calls without an explicit measurement region
def test_lhcb_1277076_requires_a_region():
    with pytest.raises(ValueError, match="dataset region is required"):
        LHCB_1277076_READER.read(
            str(LHCB_1277076_DATA / "Table2.json"),
            state="J/psi",
        )


# Reject acceptance corrections that are absent from the original HEPData tables
@pytest.mark.parametrize("state", ["J/psi", "psi(2S)"])
def test_lhcb_1277076_requires_hepdata_region(state):
    with pytest.raises(ValueError, match="unsupported region"):
        LHCB_1277076_READER.read(str(LHCB_1277076_DATA / "Table2.json"), state=state, region="central_system")


# Check the official LHCb fiducial Upsilon rapidity cross sections
def test_lhcb_1373746_fiducial_density():
    path = LHCB_1373746_DATA / "Table1.json"
    data = LHCB_1373746_READER.read(str(path), state="Upsilon(1S)", region="fiducial")
    _, errors = assert_original(data, path, density=True)
    assert data["y_err_stat"] == pytest.approx(errors["stat"].mean(axis=1))
    assert data["y_err_syst"] == pytest.approx(errors["sys"].mean(axis=1))
    assert {source["name"] for source in data["uncertainties"]} == {"statistical", "systematic"}
    covariance = sum(source_covariance(source) for source in data["uncertainties"])
    np.testing.assert_allclose(covariance, np.diag(np.square(data["y_err"])))


# Check the official acceptance-corrected LHCb Upsilon rapidity table
def test_lhcb_1373746_acceptance():
    path = LHCB_1373746_DATA / "Table2.json"
    data = LHCB_1373746_READER.read(str(path), state="Upsilon(1S)", region="central_system")
    _, errors = assert_original(data, path)
    assert data["y_err"] == pytest.approx(errors["error"].mean(axis=1))
    assert data["uncertainties"][0]["name"] == "combined"
    assert data["uncertainties"][0]["correlation"] == "uncorrelated"


# Reject a region paired with the wrong Upsilon HEPData table
def test_lhcb_1373746_rejects_region_table_mismatch():
    with pytest.raises(ValueError, match="requires Table2.json"):
        LHCB_1373746_READER.read(
            str(LHCB_1373746_DATA / "Table1.json"),
            state="Upsilon(1S)",
            region="central_system",
        )


# Validate the raw unfolded ATLAS dimuon acoplanarity counts
def test_atlas_1377585_acoplanarity_shape():
    path = ATLAS_DATA / "Table6.json"
    data = ATLAS_READER.read(str(path))
    _, errors = assert_original(data, path)
    for name, label in (("stat", "stat"), ("syst", "sys")):
        np.testing.assert_allclose(data[f"y_err_{name}"][data["valid"]], errors[label].mean(axis=1))
    assert np.all(np.isfinite(data["y_err"]))
    assert np.any(~data["valid"])
    np.testing.assert_allclose(data["y"][~data["valid"]], 0.0)


# Validate the raw unfolded ATLAS dielectron acoplanarity counts
def test_atlas_1377585_dielectron_acoplanarity_shape():
    path = ATLAS_DATA / "Table5.json"
    data = ATLAS_READER.read(str(path))
    _, errors = assert_original(data, path)
    for name, label in (("stat", "stat"), ("syst", "sys")):
        np.testing.assert_allclose(data[f"y_err_{name}"][data["valid"]], errors[label].mean(axis=1))
    assert np.all(np.isfinite(data["y_err"]))
    assert np.any(~data["valid"])
    np.testing.assert_allclose(data["y"][~data["valid"]], 0.0)


# Check that ATLAS source covariances reproduce every regularized marginal error
@pytest.mark.parametrize("table", ["Table5.json", "Table6.json"])
def test_atlas_1377585_error_sources_close(table):
    data = ATLAS_READER.read(str(ATLAS_DATA / table))
    sources = {source["name"]: source for source in data["uncertainties"]}
    assert np.sqrt(np.diag(source_covariance(sources["statistical"]))) == pytest.approx(
        data["y_err_stat"]
    )
    assert np.sqrt(np.diag(source_covariance(sources["systematic"]))) == pytest.approx(
        data["y_err_syst"]
    )
    covariance = sum(
        (source_covariance(source) for source in data["uncertainties"]),
        start=np.zeros((len(data["y"]), len(data["y"])), dtype=float),
    )
    assert np.sqrt(np.diag(covariance)) == pytest.approx(data["y_err"])


# Check ATLAS up and down columns are symmetrized before covariance construction
def test_atlas_1377585_sym_errors(tmp_path):
    table = tmp_path / "atlas.json"
    content = json.loads((ATLAS_DATA / "Table5.json").read_bytes())
    content["values"] = content["values"][:1]
    content["values"][0]["y"][0] = {"group": 0, "value": "10.0", "errors": [
        {"label": "stat", "asymerror": {"plus": "2.0", "minus": "-4.0"}},
        {"label": "sys", "asymerror": {"plus": "3.0", "minus": "-5.0"}},
    ]}
    table.write_text(json.dumps(content), encoding="utf-8")

    data = ATLAS_READER.read(str(table))

    assert data["y_err_stat"] == pytest.approx([3.0])
    assert data["y_err_syst"] == pytest.approx([4.0])
    assert data["y_err"] == pytest.approx([5.0])


# Check the ATLAS count distribution becomes a unit normalized density
def test_atlas_1377585_count_density_unit_integral():
    data = ATLAS_READER.read(str(ATLAS_DATA / "Table6.json"))
    data["scale"] = 1.0
    data["fitw"] = 1.0
    observable = {
        "acoplanarity": {
            "xlabel": "$A_{phi}$",
            "ylabel": "$N$",
            "units": {"x": "unit", "y": "unit"},
        }
    }

    plotted = plot.histhepdata(
        hepdata={"acoplanarity": data},
        obs=observable,
        density=True,
        MC_XS_SCALE=1.0,
    )["acoplanarity"]

    histogram = plotted["hdata"]
    assert histogram.integral() == pytest.approx(1.0)
    assert plotted["obs"]["units"]["y"] == "1"
    assert histogram.counts_scaled[0] == pytest.approx(
        data["y"][0] / np.sum(data["y"] * data["binwidth"])
    )


# Compute the CMS pT category key for one pair of category centers
def cms_category_key(dataset, p1_center, p2_center):
    for key, value in dataset["x"].items():
        if value[0] == pytest.approx(p1_center) and value[1] == pytest.approx(p2_center):
            return key
    raise AssertionError(f"Missing CMS category ({p1_center}, {p2_center})")


# Compute one named uncertainty source from a CMS category
def cms_uncertainty_source(dataset, key, name):
    return next(source for source in dataset["uncertainties"][key] if source["name"] == name)


# Compute mirrored invalid dPhi bin centers from the HEPData validity mask
def cms_mirrored_invalid_dphi_bins(dataset):
    invalid = set()
    for key, value in dataset["x"].items():
        p1_center = round(float(value[0]), 3)
        p2_center = round(float(value[1]), 3)
        for center, valid in zip(value[2], dataset["valid"][key], strict=True):
            if valid:
                continue
            phi_center = round(float(center), 6)
            invalid.add((p1_center, p2_center, phi_center))
            if p1_center != p2_center:
                invalid.add((p2_center, p1_center, phi_center))
    return invalid


# Check that dPhi categories keep the published common angular grid
def test_cms_2752118_dphi_padding():
    path = CMS_DATA / "Figures9-10-1115-16-17.json"
    dataset = CMS_READER.read(str(path))
    rows, _, _ = original(path)
    count = len({row["x"][2]["value"] for row in rows})
    expected_edges = np.linspace(0.0, np.pi, count + 1)
    expected_centers = (expected_edges[:-1] + expected_edges[1:]) / 2
    missing = len(dataset["bins"]) * count - len(rows)
    for key in dataset["bins"]:
        np.testing.assert_array_equal(dataset["bins"][key][2], expected_edges)
        np.testing.assert_array_equal(dataset["x"][key][2], expected_centers)
    assert sum(np.count_nonzero(~mask) for mask in dataset["valid"].values()) == missing
    assert sum(np.count_nonzero((dataset["y"][key] == 0) & (dataset["y_err"][key] == 0)) for key in dataset["y"]) == missing


# Check that the event-level dPhi cut rejects exactly the missing HEPData bins
def test_cms_2752118_dphi_acceptance():
    dataset = CMS_READER.read(str(CMS_DATA / "Figures9-10-1115-16-17.json"))

    expected_invalid = cms_mirrored_invalid_dphi_bins(dataset)
    phi_centers = sorted(
        {round(float(center), 6) for value in dataset["x"].values() for center in value[2]}
    )
    pt_pairs = {
        (round(float(value[0]), 3), round(float(value[1]), 3)) for value in dataset["x"].values()
    }
    pt_pairs |= {(p2, p1) for p1, p2 in pt_pairs}
    actual_invalid = {
        (p1, p2, phi)
        for p1, p2 in pt_pairs
        for phi in phi_centers
        if not bool(CMS_CUTS.dPhi_cuts(p1, p2, phi))
    }

    assert actual_invalid == expected_invalid


# Check that mass categories keep the locally published category span
@pytest.mark.parametrize("table", ["Figures19-20-21.json", "Figures22-23.json"])
def test_cms_2752118_mass_categories(table):
    path = CMS_DATA / table
    dataset = CMS_READER.read(str(path))
    rows, _, _ = original(path)
    coordinates = np.asarray([[axis["value"] for axis in row["x"]] for row in rows], dtype=float)
    width = np.median(np.diff(np.unique(coordinates[:, 2])))
    for key, centers in dataset["x"].items():
        selected = np.all(np.isclose(coordinates[:, :2], [centers[0], centers[1]]), axis=1)
        measured = np.sort(coordinates[selected, 2])
        edges = dataset["bins"][key][2]
        np.testing.assert_allclose(centers[2][dataset["valid"][key]], measured, atol=1e-12)
        np.testing.assert_allclose(edges[[0, -1]], measured[[0, -1]] + np.array([-1, 1]) * width / 2, atol=1e-12)
        np.testing.assert_allclose(np.diff(edges), width, atol=1e-12)


# Preserve the original mass intervals and densities after a range selection
def test_cms_2752118_mass_padding():
    filename = str(CMS_DATA / "Figures19-20-21.json")
    source = CMS_READER.read(filename)
    rows, _, _ = original(filename)
    centers = np.unique([float(row["x"][2]["value"]) for row in rows])
    selected = centers[len(centers) // 2:3 * len(centers) // 4]
    dataset = CMS_READER.read(filename, file_filter=f"`m [GeV]` >= {selected[0]} and `m [GeV]` <= {selected[-1]}")
    for key, values in dataset["x"].items():
        old = source["x"][key][2]
        mask = (old >= selected[0] - 1e-12) & (old <= selected[-1] + 1e-12)
        np.testing.assert_allclose(values[2], old[mask])
        np.testing.assert_allclose(dataset["y"][key], source["y"][key][mask])
        np.testing.assert_allclose(dataset["binwidth"][key][2], source["binwidth"][key][2][mask])


# Check that selecting one proton category retains all three published bin widths
@pytest.mark.parametrize(
    "table", ["Figures9-10-1115-16-17.json", "Figures19-20-21.json", "Figures22-23.json"]
)
def test_cms_2752118_category_filter(table):
    filename = str(CMS_DATA / table)
    source = CMS_READER.read(filename)
    rows, _, _ = original(filename)
    p1, p2 = min(tuple(float(axis["value"]) for axis in row["x"][:2]) for row in rows)
    selected = CMS_READER.read(
        filename, file_filter=f"`p_{{1,T}} [GeV]` == {p1} and `p_{{2,T}} [GeV]` == {p2}"
    )
    key = cms_category_key(source, p1, p2)
    assert len(selected["y"]) == 1
    for axis in range(3):
        np.testing.assert_allclose(selected["bins"][0][axis], source["bins"][key][axis])
    for field in ("y", "y_err", "valid"):
        np.testing.assert_array_equal(selected[field][0], source[field][key])


# Check that a single selected center retains its original edges and normalization
@pytest.mark.parametrize("table", ["Figures19-20-21.json", "Figures22-23.json"])
def test_cms_2752118_single_bin_filter(table):
    filename = str(CMS_DATA / table)
    raw = json.loads(Path(filename).read_bytes())
    centers = np.unique([float(row["x"][2]["value"]) for row in raw["values"]])
    center = centers[0]
    width = np.median(np.diff(centers))
    dataset = CMS_READER.read(filename, file_filter=f"`{raw['headers'][2]['name']}` == {center}")
    expected_categories = {tuple(axis["value"] for axis in row["x"][:2]) for row in raw["values"] if np.isclose(float(row["x"][2]["value"]), center)}
    assert len(dataset["bins"]) == len(expected_categories)
    for key in dataset["bins"]:
        np.testing.assert_allclose(dataset["bins"][key][2], center + np.array([-1, 1]) * width / 2, atol=1e-12)
        np.testing.assert_allclose(dataset["binwidth"][key][2], [width])
        assert dataset["valid"][key].tolist() == [True]


# Check rebinning preserves published intervals without crossing acceptance gaps
@pytest.mark.parametrize("factor", [2, 4])
@pytest.mark.parametrize(
    "table",
    ["Figures9-10-1115-16-17.json", "Figures19-20-21.json", "Figures22-23.json"],
)
def test_cms_2752118_rebin_published_intervals(table, factor):
    filename = str(CMS_DATA / table)
    source = CMS_READER.read(filename)
    rebinned = CMS_READER.read(filename, rebin_factor=factor)

    for key in source["valid"]:
        old_valid = source["valid"][key]
        old_edges = source["bins"][key][2]
        new_edges = rebinned["bins"][key][2]
        groups = np.asarray(
            [np.flatnonzero(np.isclose(old_edges, edge))[0] for edge in new_edges]
        )
        new_valid = rebinned["valid"][key]

        for index, (start, stop) in enumerate(zip(groups[:-1], groups[1:], strict=True)):
            group = old_valid[start:stop]
            assert 0 < len(group) <= factor
            assert np.all(group == group[0])
            assert new_valid[index] == group[0]

        old_width = np.sum(np.diff(old_edges)[old_valid])
        new_width = np.sum(np.diff(new_edges)[new_valid])
        assert new_width == pytest.approx(old_width)
        old_integral = np.sum(source["y"][key][old_valid] * np.diff(old_edges)[old_valid])
        new_integral = np.sum(
            rebinned["y"][key][new_valid] * np.diff(new_edges)[new_valid]
        )
        assert new_integral == pytest.approx(old_integral)
        for old_source, new_source in zip(
            source["uncertainties"][key], rebinned["uncertainties"][key], strict=True
        ):
            old_variance = np.diff(old_edges) @ source_covariance(old_source) @ np.diff(old_edges)
            new_variance = np.diff(new_edges) @ source_covariance(new_source) @ np.diff(new_edges)
            assert new_variance == pytest.approx(old_variance)

    rows, _, _ = original(filename)
    p1, p2 = max(tuple(float(axis["value"]) for axis in row["x"][:2]) for row in rows)
    key = cms_category_key(rebinned, p1, p2)
    assert np.any(rebinned["valid"][key])


# Check every short acceptance pattern produces homogeneous complete groups
def test_cms_2752118_rebin_partition_exhaustive():
    for size in range(1, 9):
        for bits in range(1 << size):
            valid = np.asarray([(bits >> index) & 1 for index in range(size)], dtype=bool)
            for factor in range(1, size + 3):
                groups = rebin_partition(valid, factor)
                assert groups[0] == 0
                assert groups[-1] == size
                assert np.all(np.diff(groups) > 0)
                for start, stop in zip(groups[:-1], groups[1:], strict=True):
                    group = valid[start:stop]
                    assert 0 < len(group) <= factor
                    assert np.all(group == group[0])


# Check signed global CMS normalization and efficiency shifts before rebinning
def test_cms_2752118_signed_systematics():
    path = CMS_DATA / "Figures19-20-21.json"
    dataset = CMS_READER.read(str(path))
    rows, quoted, errors = original(path)
    coordinates = np.asarray([[axis["value"] for axis in row["x"]] for row in rows], dtype=float)
    negative_bins = 0
    for key, values in dataset["y"].items():
        selected = np.all(np.isclose(coordinates[:, :2], [dataset["x"][key][axis] for axis in (0, 1)]), axis=1)
        indices = np.flatnonzero(selected)[np.argsort(coordinates[selected, 2])]
        valid = dataset["valid"][key]
        np.testing.assert_allclose(values[valid], quoted[indices])
        negative_bins += np.count_nonzero(values < 0)
        for name, label in (("normalization", "syst (norm.)"), ("efficiency", "syst (effic.)")):
            source = cms_uncertainty_source(dataset, key, name)
            assert source["correlation"] == "collective"
            assert source["effect"] == "multiplicative"
            assert source["scope"] == f"CMS_2752118:{name}"
            expected = np.sign(quoted[indices]) * errors[label][indices].mean(axis=1)
            np.testing.assert_allclose(source["shift"][valid], expected)
            for side in ("up", "down"):
                np.testing.assert_allclose(source[side], np.abs(source["shift"]))
            assert np.all(source["shift"] * values >= 0)
    assert negative_bins == np.count_nonzero(quoted < 0)


# Preserve the published collective efficiency shift and its covariance after rebinning
def test_cms_2752118_efficiency_rebin():
    path = CMS_DATA / "Figures19-20-21.json"
    raw = CMS_READER.read(str(path))
    dataset = CMS_READER.read(str(path), rebin_factor=4)
    for key, values in dataset["y"].items():
        before = cms_uncertainty_source(raw, key, "efficiency")
        after = cms_uncertainty_source(dataset, key, "efficiency")
        original_width = raw["binwidth"][key][2]
        new_width = dataset["binwidth"][key][2]
        assert new_width @ after["shift"] == pytest.approx(original_width @ before["shift"])
        assert new_width @ source_covariance(after) @ new_width == pytest.approx(
            original_width @ source_covariance(before) @ original_width)
        assert dataset["y_err_sexp"][key].shape == values.shape
        np.testing.assert_allclose(dataset["y_err_sexp"][key], np.abs(after["shift"]))


# Check that max(that,uhat) categories keep the locally published span
def test_cms_2752118_tu_categories():
    filename = str(CMS_DATA / "Figures22-23.json")
    dataset = CMS_READER.read(filename)
    rows, _, _ = original(filename)
    counts = {}
    for row in rows:
        category = tuple(float(axis["value"]) for axis in row["x"][:2])
        counts[category] = counts.get(category, 0) + 1
    for key, centers in dataset["x"].items():
        assert np.count_nonzero(dataset["valid"][key]) == counts[(centers[0], centers[1])]


# Check that missing internal max(that,uhat) rows are masked rather than merged
def test_cms_2752118_tu_sparse_bins():
    path = CMS_DATA / "Figures22-23.json"
    dataset = CMS_READER.read(str(path))
    rows, _, _ = original(path)
    coordinates = np.asarray([[axis["value"] for axis in row["x"]] for row in rows], dtype=float)
    missing = 0
    for key, centers in dataset["x"].items():
        selected = np.all(np.isclose(coordinates[:, :2], [centers[0], centers[1]]), axis=1)
        measured = coordinates[selected, 2]
        expected = np.any(np.isclose(centers[2][:, None], measured), axis=1)
        np.testing.assert_array_equal(dataset["valid"][key], expected)
        missing += np.count_nonzero(~expected)
    assert missing > 0


# Check the ten acceptance-corrected J/psi rapidity bins
def test_lhcb_2825384_jpsi_rapidity_distribution():
    path = LHCB_2825384_DATA / "Table4.json"
    data = LHCB_2825384_READER.read(str(path), state="J/psi", region="central_system")
    assert_original(data, path)
    integrated = LHCB_2825384_READER.integrate_rapidity(data)
    assert integrated["value"] == pytest.approx(np.dot(data["y"], data["binwidth"]))
    assert integrated["unit"] == "nb"
    assert integrated["total"] > integrated["stat"]


# Check the ten acceptance-corrected psi(2S) rapidity bins
def test_lhcb_2825384_psi2s_rapidity_distribution():
    path = LHCB_2825384_DATA / "Table5.json"
    data = LHCB_2825384_READER.read(str(path), state="psi(2S)", region="central_system")
    assert_original(data, path)
    integrated = LHCB_2825384_READER.integrate_rapidity(data)
    assert integrated["value"] == pytest.approx(np.dot(data["y"], data["binwidth"]))
    assert integrated["unit"] == "nb"
    assert integrated["total"] > integrated["stat"]


# Check state selection through the same dataset metadata used by iceplot
def test_lhcb_2825384_dataset_state_selection():
    data = LHCB_2825384_READER.read(
        str(LHCB_2825384_DATA / "Table5.json"),
        dataset={"state": "psi(2S)", "region": "central_system"},
    )

    assert data["state"] == "psi(2S)"


# Reject state-table mismatches instead of silently drawing the wrong state
def test_lhcb_2825384_rejects_state_table_mismatch():
    with pytest.raises(ValueError, match="requires Table5.json"):
        LHCB_2825384_READER.read(
            str(LHCB_2825384_DATA / "Table4.json"),
            state="psi(2S)",
            region="central_system",
        )


# Reject fiducial scalar transcriptions outside the original HEPData record
@pytest.mark.parametrize("state", ["J/psi", "psi(2S)"])
def test_lhcb_2825384_requires_hepdata_region(state):
    with pytest.raises(ValueError, match="unsupported region"):
        LHCB_2825384_READER.read(str(LHCB_2825384_DATA / "Table4.json"), state=state, region="fiducial")


# Preserve luminosity shared across the two states directly from the HEPData columns
def test_lhcb_2825384_luminosity_scope():
    sources = [named_uncertainty(LHCB_2825384_READER.read(str(LHCB_2825384_DATA / table),
               state=state, region="central_system"), "luminosity")
               for state, table in (("J/psi", "Table4.json"), ("psi(2S)", "Table5.json"))]
    assert sources[0]["scope"] == sources[1]["scope"]


# Reject ambiguous 13 TeV reader calls without a region
def test_lhcb_2825384_requires_a_region():
    with pytest.raises(ValueError, match="dataset region is required"):
        LHCB_2825384_READER.read(
            str(LHCB_2825384_DATA / "Table4.json"),
            state="J/psi",
        )
