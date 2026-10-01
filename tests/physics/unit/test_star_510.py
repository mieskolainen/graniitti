# STAR 510 GeV HEPData spectra and physical event selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering
from core.io.cache import freeze

from icepack._common.hepdata import error_pairs, table_values
from icepack.SOFTCEP.STAR_3075716._common import cuts, obs, reader
from tests.technical.support.hepmc import make_pair_event

ROOT = Path(__file__).resolve().parents[3]
DATA = ROOT / "HEPData/CEP/HEPData-ins3075716-v1-json"


# Compare every measured value and asymmetric error with its original HEPData entry
@pytest.mark.parametrize("path", sorted(DATA.glob("*.json")), ids=lambda path: path.stem)
def test_original_spectra(path):
    table = json.loads(path.read_bytes())
    values = table_values(table)
    data = reader.read(str(path))
    valid = data["valid"]
    np.testing.assert_allclose(data["y"][valid], [float(row["value"]) for row in values])
    sources = {source["name"]: source for source in data["uncertainties"]}
    for name, label in (("statistical", "Stat. error"), ("experimental", "Syst. error")):
        expected = np.abs(error_pairs(values, label))
        np.testing.assert_allclose(sources[name]["up"][valid], expected[:, 0])
        np.testing.assert_allclose(sources[name]["down"][valid], expected[:, 1])
    assert set(sources) == {"statistical", "experimental"}
    np.testing.assert_allclose(data["y"][~valid], 0.0)
    assert data["bins"][-1] == pytest.approx(float(table["values"][-1]["x"][0]["high"]))


# Keep the unmeasured proton azimuth interval out of comparisons before and after rebinning
@pytest.mark.parametrize("factor", [None, 2])
def test_azimuth_gap(factor):
    path = DATA / "Figure11a.json"
    rows = json.loads(path.read_bytes())["values"]
    data = reader.read(str(path), rebin_factor=factor)
    gaps = [(float(left["x"][0]["high"]), float(right["x"][0]["low"]))
            for left, right in zip(rows[:-1], rows[1:], strict=True)
            if float(left["x"][0]["high"]) < float(right["x"][0]["low"])]
    assert gaps
    for low, high in gaps:
        overlap = (data["bins"][:-1] < high) & (data["bins"][1:] > low)
        assert not np.any(data["valid"][overlap])
    assert len(data["valid"]) == len(data["y"])


# Read each spectrum through the same card and bin mapping used by iceplot
@pytest.mark.parametrize("channel", ["pipi", "KK", "ppbar"])
def test_dataset_reading(channel):
    dataset, path = steering.load_dataset(f"icepack/SOFTCEP/STAR_3075716/{channel}/dataset.json", cdir=ROOT)
    for entry in dataset["sets"]:
        selection = readers.load_cut_module(steering.resolve_python_reference(
            entry["cuts"], package="core.analysis.cuts", dataset_path=path, cdir=ROOT,
        ))
        freeze(selection.cut_param, strict=True)
        observables, _ = steering.load_observables(entry["obs"], dataset_path=path, cdir=ROOT)
        data, histograms = readers.read_hepdata(
            entry, dataset["datapath"], dataset["type"], observables, cdir=ROOT,
            reader=dataset["reader"], dataset_path=path,
        )
        for spectrum in entry["hist"]:
            result = data[spectrum["obs"]]
            assert result["scale"] == pytest.approx(1e-12 if result["unit"] == "pb" else 1e-9)
            np.testing.assert_array_equal(histograms[spectrum["obs"]]["valid"], result["valid"])


# Select each published mass region using an on-shell pion pair and the real HepMC3 API
@pytest.mark.parametrize("mass,region", [(0.9, 0), (1.0, 1), (1.4, 1), (1.5, 2), (2.0, 2)])
def test_mass_regions(mass, region):
    for index in range(3):
        event = make_pair_event(mass=mass, pid=211, beam_energy=255.0)
        event.cut_param = cuts.parameters(mass=index)
        assert cuts.cut_func(event) == (index == region)


# Combine mass and azimuth selections without counting the 90 degree boundary twice
@pytest.mark.parametrize("phi,region", [(0.0, 0), (89.999999, 0), (90.0, None), (90.000001, 1), (180.0, 1)])
def test_azimuth_regions(phi, region):
    for index in range(2):
        event = make_pair_event(mass=1.2, pid=211, dphi=phi, beam_energy=255.0)
        event.cut_param = cuts.parameters(mass=1, phi=index)
        assert cuts.cut_func(event) == (index == region)


# Exclude the tagged beam protons when projecting a central proton-antiproton pair
def test_central_protons():
    event = make_pair_event(mass=2.4, beam_energy=255.0, forward_pt=0.6)
    assert obs.obs_M["func"](event) == pytest.approx(2.4)
    assert obs.obs_Rap["func"](event) == pytest.approx(0.0, abs=1e-12)
    assert obs.obs_dPhi_pp["func"](event) == pytest.approx(180.0)
