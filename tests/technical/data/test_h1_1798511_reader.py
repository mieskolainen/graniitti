# Checks for the H1 elastic pion pair photoproduction reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
from core.io import readers, steering

ROOT = Path(__file__).resolve().parents[3]
CARD = "icepack/PHOTOPROD/H1_1798511/dataset.json"
DATA = ROOT / "HEPData" / "PHOTOPROD" / "HEPData-ins1798511-v1-json"
SPECTRUM = DATA / "Table17.json"
STAT_CORR = DATA / "Table17_statCorr.json"

# Check the elastic spectrum and complete published covariance structure
def test_h1_1798511_elastic_mass_covariances():
    dataset, resolved = steering.load_dataset(CARD, cdir=ROOT)
    entry = dataset["sets"][0]
    observables, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
    spectra, observables = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], observables,
        cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved)
    data = spectra["M"]
    rows = [row for row in json.loads(SPECTRUM.read_text())["values"] if "value" in row["x"][1]]
    count = len(rows)
    edges = [float(row["x"][0]["low"]) for row in rows] + [float(rows[-1]["x"][0]["high"])]
    assert len(data["y"]) == count
    np.testing.assert_array_equal(observables["M"]["bins"], edges)
    np.testing.assert_array_equal(data["bins"], edges)
    assert len(data["uncertainties"]) == len(rows[0]["y"][0]["errors"])

    for source in data["uncertainties"]:
        covariance = np.asarray(source["covariance"], dtype=float)
        assert covariance.shape == (count, count)
        assert np.all(np.isfinite(covariance))
        assert np.allclose(covariance, covariance.T)

    statistical = data["uncertainties"][0]
    covariance = np.asarray(statistical["covariance"], dtype=float)
    assert statistical["category"] == "statistical"
    assert np.any(np.abs(covariance - np.diag(np.diag(covariance))) > 0.0)
