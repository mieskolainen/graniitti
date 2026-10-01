# Scalar MC conversions through the HEPData plotting and tuning consumers
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering
from core.plot import plot
from core.tune.drivers.graniitti.driver import GraniittiDriver

ROOT = Path(__file__).resolve().parents[3]
CARD = "icepack/UPC/PHOTOPROD/CMS_2908607/phi_coherent/dataset.json"


# Convert decay cross sections to unchanged HEPData production densities in both consumers
@pytest.mark.parametrize("card", [
    pytest.param(CARD, id="cms_phi"),
    pytest.param("icepack/PHOTOPROD/ZEUS_582237/dataset.json", id="zeus_jpsi"),
    pytest.param("icepack/UPC/PHOTOPROD/ALICE_2658375/jpsi_incoherent/dataset.json", id="alice_jpsi"),
])
def test_production_conversion(tmp_path, card):
    dataset, path = steering.load_dataset(card, cdir=ROOT)
    entry = dataset["sets"][0]
    definitions, _ = steering.load_observables(entry["obs"], dataset_path=path, cdir=ROOT)
    spectra, observables = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], definitions,
        reader=dataset["reader"], dataset_path=path, cdir=str(ROOT))
    hist = entry["hist"][0]
    name = hist["obs"]
    spectrum = spectra[name]
    expected = spectrum["y"] * hist["scale"] * 1e12
    weights = expected * spectrum["binwidth"] / entry["mc_scale"]
    mcdata = {"data": {name: spectrum["x"]}, "weights": weights, "xsection_pb": weights.sum()}
    driver = GraniittiDriver()
    driver.init_data(run_name=str(tmp_path / "conversion"), datacards=[{"datacard": card}],
                     obs_module="default", cdir=str(ROOT), pickle_dump=False)
    for obs in (observables, driver.obs[0][0]):
        result = plot.histmc(mcdata, {name: obs[name]}, scale=entry["mc_scale"])[name]["hdata"]
        np.testing.assert_allclose(result.counts_scaled, expected)
    np.testing.assert_allclose(driver.data[0][0][name]["hdata"].counts_scaled, expected)
