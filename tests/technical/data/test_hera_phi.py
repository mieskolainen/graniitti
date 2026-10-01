# Validate ZEUS phi spectra, normalization and proton dissociation cuts
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering
from core.plot import plot
from core.stats.hist import center2edgebins
from core.stats.uncertainty import source_covariance

from icepack._common.hepdata import read_table
from icepack.PHOTOPROD._common import dissociative
from tests.technical.data.test_hera_dissociation import photo_event
from tests.technical.support.hepdata import assert_original

ROOT = Path(__file__).resolve().parents[3]


# Preserve original phi values and errors through the full reader and histogram normalization
@pytest.mark.parametrize("entry", ["ZEUS_415642", "ZEUS_508770", "ZEUS_508770/dissociative"])
def test_phi_spectra(entry):
    card, resolved = steering.load_dataset(f"icepack/PHOTOPROD/{entry}/dataset.json", cdir=ROOT)
    for dataset in card["sets"]:
        definitions, _ = steering.load_observables(dataset["obs"], dataset_path=resolved, cdir=ROOT)
        data, observables = readers.read_hepdata(dataset, card["datapath"], card["type"], definitions,
            cdir=str(ROOT), reader=card["reader"], dataset_path=resolved)
        for hist in dataset["hist"]:
            name = hist["obs"]
            result = data[name]
            path = ROOT / card["datapath"] / hist["file"]
            _, errors = assert_original(result, path)
            for side, index in (("down", 0), ("up", 1)):
                np.testing.assert_allclose(result[f"y_err_{side}"],
                                           np.sqrt(sum(error[:, index]**2 for error in errors.values())))
            table = read_table(path)
            if all("low" not in row["x"][0] for row in table["values"]):
                centers = np.asarray([float(row["x"][0]["value"]) for row in table["values"]])
                np.testing.assert_allclose(result["bins"], center2edgebins(centers))
                expected = sum(np.outer(error.mean(axis=1), error.mean(axis=1))
                               if label == "sys_2" or "normalization" in label
                               else np.diag(error.mean(axis=1)**2) for label, error in errors.items())
                np.testing.assert_allclose(sum(source_covariance(source) for source in result["uncertainties"]), expected)
            expected_pb = result["y"] * hist["scale"] * 1e12
            width = result["binwidth"] if hist.get("differential", True) else np.ones_like(expected_pb)
            weights = expected_pb * width / dataset["mc_scale"]
            mc = plot.histmc({"data": {name: result["x"]}, "weights": weights, "xsection_pb": weights.sum()},
                             {name: observables[name]}, scale=dataset["mc_scale"])[name]["hdata"]
            np.testing.assert_allclose(mc.counts_scaled, expected_pb)
            assert mc.integral() == pytest.approx(expected_pb @ width)


# Apply a W dependent dissociation mass bound to exact two body HepMC events
@pytest.mark.parametrize("direction", [-1, 1])
@pytest.mark.parametrize("fragment", [False, True])
def test_dissociation_mass_fraction(direction, fragment):
    event, _ = photo_event(direction, fragment, 12.0)
    values = dissociative.kinematics(event)
    ratio = values["my2"] / values["w"]**2
    event.cut_param.pop("MY_MAX")
    event.cut_param["MY2_W2_MAX"] = 1.1 * ratio
    assert dissociative.accepted(event)
    event.cut_param["MY2_W2_MAX"] = 0.9 * ratio
    assert not dissociative.accepted(event)
