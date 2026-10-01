# Tests for elastic data uncertainty correlations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io.hepdata_reader import load_dataset_reader
from core.stats.uncertainty import source_covariance

from tests.technical.support.hepdata import original

ROOT = Path(__file__).resolve().parents[3]
EL_DATA = ROOT / "HEPData" / "EL"
READERS = {
    name: load_dataset_reader(f"icepack/ELASTIC/{name}/dataset.json", cdir=str(ROOT))
    for name in (
        "ISR_212895",
        "ISR_214689",
        "STAR_1791591",
        "TOTEM_922651",
        "TOTEM_1220862",
        "TOTEM_1489188",
        "TOTEM_1710340",
    )
}
FILES = {
    "ISR_212895": EL_DATA / "HEPData-ins212895-v1-json" / "Table1.json",
    "ISR_214689": EL_DATA / "HEPData-ins214689-v1-json" / "Table2.json",
    "STAR_1791591": EL_DATA / "HEPData-ins1791591-v1-json" / "Figure5.json",
    "TOTEM_922651": EL_DATA / "HEPData-ins922651-v1-json" / "Table1.json",
    "TOTEM_1220862": EL_DATA / "HEPData-ins1220862-v1-json" / "Table1.json",
    "TOTEM_1489188": EL_DATA / "HEPData-ins1489188-v1-json" / "Table1.json",
    "TOTEM_1710340": EL_DATA / "HEPData-ins1710340-v1-json" / "Table1.json",
}


# Read original measurement values and quoted errors independently of the production parser
def published_values(path):
    return json.loads(path.read_bytes())["values"]


# Compute the directly quoted absolute error on one side of a published source
def published_error(rows, index, side):
    errors = [row["y"][0]["errors"][index] for row in rows]
    direction = "plus" if side == "up" else "minus"
    return np.abs(np.asarray([error["symerror"] if "symerror" in error else error["asymerror"][direction]
                              for error in errors], dtype=float))


# Check ISR dip data retain measured bin edges, point errors and absolute normalization errors
@pytest.mark.parametrize("table", ["Table1.json", "Table2.json"])
def test_isr_212895_published_bins_errors(table):
    path = FILES["ISR_212895"].with_name(table)
    rows = published_values(path)
    dataset = READERS["ISR_212895"].read(str(path))
    y = np.asarray([row["y"][0]["value"] for row in rows], dtype=float)
    point_error = published_error(rows, 0, "up")
    _, _, quoted = original(path)
    normalization = quoted["sys" if table == "Table1.json" else "sys_2"][:, 1] / y

    assert dataset["x"] == pytest.approx(np.asarray([point["value"] if "value" in point else 0.5 * (float(point["low"]) + float(point["high"]))
                                                 for point in (row["x"][0] for row in rows)], dtype=float))
    assert dataset["bins"][:-1] == pytest.approx([float(row["x"][0]["low"]) for row in rows])
    assert dataset["bins"][1:] == pytest.approx([float(row["x"][0]["high"]) for row in rows])
    assert dataset["y"] == pytest.approx(y)
    assert uncertainty_source(dataset, "published")["up"] == pytest.approx(point_error)
    assert uncertainty_source(dataset, "published")["category"] == "combined"
    assert systematic_covariance(dataset) == pytest.approx(np.outer(normalization * y, normalization * y))
    assert dataset["y_err"] == pytest.approx(np.hypot(point_error, normalization * y))

    rebinned = READERS["ISR_212895"].read(str(path), rebin_factor=3)
    assert rebinned["binwidth"] @ rebinned["y"] == pytest.approx(dataset["binwidth"] @ y)
    assert systematic_covariance(rebinned) == pytest.approx(np.outer(normalization[0] * rebinned["y"], normalization[0] * rebinned["y"]))


# Check the shared normalization reproduces the published ppbar/pp relative uncertainty
def test_isr_212895_relative_normalization():
    pp, ppbar = (
        READERS["ISR_212895"].read(str(FILES["ISR_212895"].with_name(table)))
        for table in ("Table1.json", "Table2.json")
    )
    pp_scale = uncertainty_source(pp, "normalization")
    ppbar_scale = uncertainty_source(ppbar, "normalization")
    assert pp_scale["scope"] == ppbar_scale["scope"]
    cross_covariance = np.outer(pp_scale["shift"], ppbar_scale["shift"])
    fractional_covariance = cross_covariance / np.outer(pp["y"], ppbar["y"])
    _, pp_y, pp_errors = original(FILES["ISR_212895"])
    _, ppbar_y, ppbar_errors = original(FILES["ISR_212895"].with_name("Table2.json"))
    pp_error = pp_errors["sys"][:, 1] / pp_y
    ppbar_error = ppbar_errors["sys_2"][:, 1] / ppbar_y
    relative = ppbar_errors["sys_1"][:, 1] / ppbar_y
    np.testing.assert_allclose(pp_error[:, None]**2 + ppbar_error[None, :]**2 - 2 * fractional_covariance,
                               np.broadcast_to(relative**2, fractional_covariance.shape))


# Check ISR axis conversion and errors against the original HEPData table
@pytest.mark.parametrize("table", [f"Table{index}.json" for index in range(2, 8)])
def test_isr_214689_published_points_errors(table):
    path = FILES["ISR_214689"].with_name(table)
    rows = published_values(path)
    dataset = READERS["ISR_214689"].read(str(path))
    assert dataset["x"] == pytest.approx(1.0e-3 * np.asarray([row["x"][0]["value"] for row in rows], dtype=float))
    assert dataset["y"] == pytest.approx(np.asarray([row["y"][0]["value"] for row in rows], dtype=float))
    assert dataset["y_err_up"] == pytest.approx(np.abs(published_error(rows, 0, "up")))
    assert dataset["y_err_down"] == pytest.approx(np.abs(published_error(rows, 0, "down")))
    assert len(dataset["uncertainties"]) == 1
    assert dataset["uncertainties"][0]["category"] == "combined"
    covariance = sum(source_covariance(source) for source in dataset["uncertainties"])
    np.testing.assert_allclose(covariance, np.diag(dataset["y_err"]**2))


# Compute one named uncertainty source from a reader result
def uncertainty_source(dataset: dict, name: str) -> dict:
    """Return one uniquely named source"""
    return next(source for source in dataset["uncertainties"] if source["name"] == name)


# Read the symmetric systematic marginal directly from an original JSON table
def published_systematic(name: str) -> np.ndarray:
    """Return the mean absolute up and down systematic error"""
    rows = published_values(FILES[name])
    return 0.5 * (published_error(rows, 1, "up") + published_error(rows, 1, "down"))


# Sum the covariance of all systematic sources in one reader result
def systematic_covariance(dataset: dict) -> np.ndarray:
    """Return the full systematic covariance matrix"""
    covariance = np.zeros((len(dataset["y"]), len(dataset["y"])), dtype=float)
    for source in dataset["uncertainties"]:
        if source["category"] == "systematic":
            covariance += source_covariance(source)
    return covariance


# Preserve original statistical and systematic marginals without external components
@pytest.mark.parametrize("name", ["TOTEM_922651", "TOTEM_1489188", "TOTEM_1710340"])
def test_elastic_hepdata_errors(name):
    dataset = READERS[name].read(str(FILES[name]))
    rows = published_values(FILES[name])
    assert dataset["y"] == pytest.approx([float(row["y"][0]["value"]) for row in rows])
    assert dataset["y_err_syst"] == pytest.approx(published_systematic(name))
    for side in ("up", "down"):
        stat, syst = published_error(rows, 0, side), published_error(rows, 1, side)
        assert dataset[f"y_err_stat_{side}"] == pytest.approx(stat)
        assert uncertainty_source(dataset, "systematic")[side] == pytest.approx(syst)
        assert dataset[f"y_err_{side}"] == pytest.approx(np.hypot(stat, syst))
    assert {source["name"] for source in dataset["uncertainties"]} == {"statistical", "systematic"}
    covariance = systematic_covariance(dataset)
    np.testing.assert_allclose(covariance, np.diag(published_systematic(name)**2))


# Separate common luminosity without counting it twice in published elastic total errors
@pytest.mark.parametrize("name,companion,label", [
    ("STAR_1791591", "Cross-Sections.json", "SysUnc luminosity"),
    ("TOTEM_1220862", "Table4.json", "sys_3"),
])
def test_elastic_luminosity_covariance(name, companion, label):
    path = FILES[name]
    rows = published_values(path)
    dataset = READERS[name].read(str(path))
    normalization = published_values(path.with_name(companion))[0]["y"][0]
    error = next(item for item in normalization["errors"] if item["label"] == label)
    common, residual = [], []
    for side, direction in (("up", "plus"), ("down", "minus")):
        magnitude = error["symerror"] if "symerror" in error else error["asymerror"][direction]
        lumi = dataset["y"] * abs(float(magnitude)) / float(normalization["value"])
        stat, syst = published_error(rows, 0, side), published_error(rows, 1, side)
        common.append(lumi)
        residual.append(np.sqrt(syst**2 - lumi**2))
        np.testing.assert_allclose(uncertainty_source(dataset, "luminosity")[side], lumi)
        np.testing.assert_allclose(uncertainty_source(dataset, "systematic")[side], residual[-1])
        np.testing.assert_allclose(dataset[f"y_err_{side}"], np.hypot(stat, syst))
    shift = np.mean(common, axis=0)
    covariance = np.diag(np.mean(residual, axis=0)**2) + np.outer(shift, shift)
    np.testing.assert_allclose(systematic_covariance(dataset), covariance)
    widths = dataset["binwidth"]
    np.testing.assert_allclose(widths @ covariance @ widths,
                               np.sum((widths * np.mean(residual, axis=0))**2) + (widths @ shift)**2)
    assert uncertainty_source(dataset, "luminosity")["effect"] == "multiplicative"


# Propagate the original covariance through elastic density rebinning
@pytest.mark.parametrize("name", ["STAR_1791591", "TOTEM_1220862", "TOTEM_922651", "TOTEM_1710340"])
def test_elastic_rebin_covariance(name):
    original = READERS[name].read(str(FILES[name]))
    rebinned = READERS[name].read(str(FILES[name]), rebin_factor=3)
    assert rebinned["binwidth"] @ rebinned["y"] == pytest.approx(original["binwidth"] @ original["y"])
    variance = original["binwidth"] @ systematic_covariance(original) @ original["binwidth"]
    rebinned_variance = rebinned["binwidth"] @ systematic_covariance(rebinned) @ rebinned["binwidth"]
    assert rebinned_variance == pytest.approx(variance)


# Check the physical slope on linear and logarithmic grids including the inserted forward point
@pytest.mark.parametrize("q2", [np.linspace(1e-15, 1.0, 100), np.geomspace(1e-15, 1.0, 100)])
def test_slope_momentum_interval(q2):
    from tests.physics.studies.elastic.analysis.analyze_momentum import local_slope

    q2 = np.concatenate(([0.0], q2))
    grid, slope = local_slope(q2[:, None], np.exp(-20.0 * q2))
    # The first logarithmic intervals have unavoidable floating point cancellation
    np.testing.assert_allclose(slope[grid > 1e-8], 20.0, rtol=1e-8)
