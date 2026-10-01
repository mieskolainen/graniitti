# Checks for the HERA exclusive J/psi photoproduction data bundles
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.io.hepdata_reader import load_dataset_reader
from core.stats.uncertainty import source_covariance

from icepack.PHOTOPROD._common import hera_reader
from tests.technical.support.hepdata import assert_original

ROOT = Path(__file__).resolve().parents[3]
H1_HIGH_CARD = "icepack/PHOTOPROD/H1_1228913/high_energy/dataset.json"
H1_LOW_CARD = "icepack/PHOTOPROD/H1_1228913/low_energy/dataset.json"
ZEUS_CARD = "icepack/PHOTOPROD/ZEUS_582237/dataset.json"
H1_DATA = ROOT / "HEPData" / "PHOTOPROD" / "HEPData-ins1228913-v1-json"
ZEUS_DATA = ROOT / "HEPData" / "PHOTOPROD" / "HEPData-ins582237-v1-json"

# Check both H1 elastic momentum-transfer measurements
@pytest.mark.parametrize("card,table", [
    (H1_HIGH_CARD, "Table5.json"), (H1_LOW_CARD, "Table7.json"),
    ("icepack/PHOTOPROD/H1_1228913/dissociative_high_energy/dataset.json", "Table6.json"),
    ("icepack/PHOTOPROD/H1_1228913/dissociative_low_energy/dataset.json", "Table8.json"),
])
def test_h1_1228913_elastic_t_spectra(card, table):
    path = H1_DATA / table
    data = load_dataset_reader(card, cdir=str(ROOT)).read(str(path))
    _, errors = assert_original(data, path)
    for side, index in (("down", 0), ("up", 1)):
        np.testing.assert_allclose(data[f"y_err_{side}"], errors["error"][:, index])
    assert data["uncertainties"][0]["category"] == "combined"


# Check the selected ZEUS W bin and its published t intervals
def test_zeus_582237_elastic_t_spectrum_100_gev():
    path = ZEUS_DATA / "Table9.json"
    data = load_dataset_reader(ZEUS_CARD, cdir=str(ROOT)).read(str(path), file_filter="90 TO 110 GEV")
    assert data["group"] == 3
    _, errors = assert_original(data, path, group=3)
    for side, index in (("down", 0), ("up", 1)):
        np.testing.assert_allclose(data[f"y_err_{side}"], errors["error"][:, index])

# Preserve HEPData marginal uncertainties without supplementary text inputs
@pytest.mark.parametrize("table", ["Table5.json", "Table7.json"])
def test_h1_covariance(table):
    data = hera_reader.read(H1_DATA / table)
    _, errors = assert_original(data, H1_DATA / table)
    expected = np.diag(errors["error"].mean(axis=1)**2)
    np.testing.assert_allclose(sum(source_covariance(item) for item in data["uncertainties"]), expected)
