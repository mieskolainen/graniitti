# Integrated and differential literature reader physics and input checks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io.hepdata_reader import load_dataset_reader
from core.stats.uncertainty import source_covariance

from icepack._common import reference_reader
from tests.technical.support.hepdata import original

REFERENCES = Path(__file__).resolve().parents[2] / "references"


# Preserve the separately quoted lower and upper STAR systematic errors
def test_star_asymmetric_errors():
    root = REFERENCES.parents[1]
    path = root / "HEPData/CEP/HEPData-ins1792394-v1-json/Table2.json"
    reader = load_dataset_reader("icepack/SOFTCEP/integrated/STAR_2020/dataset.json", cdir=str(root))
    data = reader.read(str(path), file_filter="P P --> P PI+ PI- P", hist={"rows": [1]})
    _, values, quoted = original(path, rows=[1])
    np.testing.assert_allclose(data["y"], values)
    errors = {source["name"]: source for source in data["uncertainties"]}
    for index, side in enumerate(("down", "up")):
        for name in ("stat", "syst"):
            np.testing.assert_allclose(errors[name][side], quoted[f"{name}."][:, index])
        np.testing.assert_allclose(data[f"y_err_{side}"], np.sqrt(sum(pair[:, index]**2 for pair in quoted.values())))
    variance = sum(np.mean(pair, axis=1)**2 for pair in quoted.values())
    np.testing.assert_allclose(sum(source_covariance(source) for source in errors.values()), np.diag(variance))


# Retain the luminosity uncertainty independently of CMS systematic errors
def test_cms_luminosity_error():
    root = REFERENCES.parents[1]
    path = root / "HEPData/CEP/HEPData-ins954992-v1-json/Table1.json"
    reader = load_dataset_reader("icepack/GAMMA/integrated/CMS_2012/dataset.json", cdir=str(root))
    data = reader.read(str(path))
    _, values, errors = original(path)
    label = next(label for label in errors if "luminosity" in label)
    luminosity = next(source for source in data["uncertainties"] if source["name"] == "lumi")
    np.testing.assert_allclose(data["y"], values)
    for index, side in enumerate(("down", "up")):
        np.testing.assert_allclose(luminosity[side], errors[label][:, index])
    assert luminosity["effect"] == "multiplicative"
    np.testing.assert_allclose(data["y_err"], np.sqrt(sum(pair.mean(axis=1)**2 for pair in errors.values())))


# Reconstruct the native SuperChic mass density from the preserved EPS coordinates
@pytest.mark.parametrize("channel", ["gg", "qqbar_one_massless_flavour", "bbbar", "ggg", "g_qqbar_one_massless_flavour"])
def test_superchic_mass_density(channel):
    path = REFERENCES / "superchic2_1508_02718.json"
    key = f"figure3_{channel}"
    row = json.loads(path.read_text())[key]
    native = row["digitization"]
    mapping = native["coordinate_mapping"]
    expected = 10.0 ** np.interp(native["y_plot"], mapping["y_plot"], mapping["log10_dsigma_pb_per_gev"])
    data = reference_reader.read(str(path), key)
    assert data["y"] == pytest.approx(expected)
    assert data["xlim"] == pytest.approx(native["native_mass_range_gev"])
    integral = np.dot(expected, np.diff(row["bins"]))
    assert np.dot(data["y"], data["binwidth"]) == pytest.approx(integral)
    assert not data["uncertainties"]


# Keep theory predictions distinct from measurements with zero experimental uncertainty
def test_integrated_theory():
    path, key = REFERENCES / "superchic2_1508_02718.json", "table1_gg_m75"
    data = reference_reader.read(str(path), key)
    assert np.dot(data["y"], data["binwidth"]) == pytest.approx(json.loads(path.read_text())[key]["measurement"]["value"])
    assert not data["uncertainties"]
    assert data["reference"]


# Reject upper limits and confidence intervals as ordinary Gaussian measurements
@pytest.mark.parametrize("key", ["cms_totem_exclusive_ww_13tev", "oxford_chic1_90cl", "cms_quasiexclusive_ww_8tev"])
def test_non_gaussian_reference(key):
    with pytest.raises(ValueError, match="not a measured value or prediction"):
        reference_reader.read(str(REFERENCES / "exclusive_measurements.json"), key)


# Reject missing selectors and unsupported rebinning through the real reader API
@pytest.mark.parametrize("options,match", [({}, "file_filter is required"),
    ({"file_filter": "absent"}, "unknown reference ID"),
    ({"file_filter": "figure3_gg", "rebin_factor": 2}, "rebinning is not supported")])
def test_reference_selection(options, match):
    with pytest.raises(ValueError, match=match):
        reference_reader.read(str(REFERENCES / "superchic2_1508_02718.json"), **options)


# Reject malformed numerical input before constructing a comparison histogram
@pytest.mark.parametrize("field,value,match", [
    ("type", "other", "unknown type"),
    ("unit", 5, "invalid cross section unit"),
    ("bins", [0.0, 0.0], "invalid bin edges"),
    ("value", [1.0], "invalid measurement values"),
    ("stat", [-0.1, 0.1], "invalid stat uncertainties"),
    ("syst", [0.1], "invalid syst uncertainties"),
])
def test_invalid_reference(tmp_path, field, value, match):
    row = json.loads((REFERENCES / "superchic2_1508_02718.json").read_text())["figure3_gg"]
    target = row["measurement"] if field in {"value", "stat", "syst"} else row
    target[field] = value
    path = tmp_path / "reference.json"
    path.write_text(json.dumps({"test": row}))
    with pytest.raises(ValueError, match=match):
        reference_reader.read(str(path), "test")


# Propagate bin errors with separate statistical and correlated systematic components
def test_differential_covariance(tmp_path):
    row = json.loads((REFERENCES / "superchic2_1508_02718.json").read_text())["figure3_gg"]
    row["bins"] = row["bins"][:3]
    row["measurement"]["value"] = row["measurement"]["value"][:2]
    row["measurement"]["stat"] = [[0.1, 0.2], [0.3, 0.4]]
    row["measurement"]["syst"] = [[0.5, 0.6], [0.7, 0.8]]
    path = tmp_path / "reference.json"
    path.write_text(json.dumps({"test": row}))
    data = reference_reader.read(str(path), "test")
    covariance = sum(source_covariance(source) for source in data["uncertainties"])
    np.testing.assert_allclose(covariance, [[0.40, 0.42], [0.42, 0.58]])
    np.testing.assert_allclose(data["y_err_down"]**2, [0.26, 0.40])
    np.testing.assert_allclose(data["y_err_up"]**2, [0.58, 0.80])


# Reject copied measured values through the real theory reference reader
@pytest.mark.parametrize("filename", ["lhc_tevatron_rhic.json", "exclusive_measurements.json"])
def test_measured_reference_rejected(filename):
    path = REFERENCES / filename
    for key, row in json.loads(path.read_text()).items():
        if row.get("result") == "measurement":
            with pytest.raises(ValueError, match="original experimental data reader"):
                reference_reader.read(str(path), key)
