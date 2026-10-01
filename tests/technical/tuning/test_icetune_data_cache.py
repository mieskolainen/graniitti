# Test icetune HEPData cache input fingerprints
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io import steering
from core.tune.drivers.graniitti import driver as graniitti_driver

ROOT = Path(__file__).resolve().parents[3]


# Fingerprint dataset-relative readers, cuts, observables and original tables
def test_dataset_input_records():
    driver = graniitti_driver.GraniittiDriver()
    records = driver._data_input_file_records(
        cdir=ROOT,
        datacards=[{"datacard": "icepack/PHOTOPROD/LHCb_2825384/jpsi/dataset.json"}],
    )
    paths = {record["path"] for record in records}

    assert "icepack/PHOTOPROD/LHCb_2825384/jpsi/dataset.json" in paths
    assert "icepack/PHOTOPROD/LHCb_2825384/reader.py" in paths
    assert "icepack/PHOTOPROD/LHCb_2825384/cuts.py" in paths
    assert "icepack/PHOTOPROD/LHCb_2825384/obs.py" in paths
    assert "HEPData/CEP/HEPData-ins2825384-v1-json/Table4.json" in paths
    assert "icepack/_common/hepdata.py" in paths


# Invalidate serialized measurement histograms after changes to cuts, readers or MC scales
@pytest.mark.parametrize("changed", ["cuts", "reader", "mc_scale"])
def test_cache_input_fingerprint(tmp_path, data_card, changed):
    settings = dict(run_name='fingerprint', datacards=[{'datacard': str(data_card)}], obs_module='default', cdir=str(tmp_path))
    writer = graniitti_driver.GraniittiDriver()
    writer.init_data(**settings)
    matching = graniitti_driver.GraniittiDriver()
    assert matching.load_data_cache(**settings) and matching.initialized
    assert matching.pid == writer.pid
    for source, restored in zip(writer.data[0], matching.data[0], strict=True):
        for name in source:
            np.testing.assert_array_equal(source[name]['hdata'].counts_scaled, restored[name]['hdata'].counts_scaled)
            np.testing.assert_array_equal(source[name]['hdata'].covariance_scaled, restored[name]['hdata'].covariance_scaled)
    if changed == "mc_scale":
        dataset = json.loads(data_card.read_text())
        dataset["sets"][0]["mc_scale"] *= 2.0
        data_card.write_text(json.dumps(dataset))
    else:
        source = data_card.parent.parent / 'cuts.py' if changed == "cuts" else tmp_path / 'icepack/_common/hepdata.py'
        source.write_text(source.read_text() + '\n# Changed measurement processing\n')
    assert not graniitti_driver.GraniittiDriver().load_data_cache(**settings)


# Include explicitly referenced covariance and photon flux tables in the cache fingerprint
@pytest.mark.parametrize("field,card", [
    ("covariance", "icepack/PHOTOPROD/H1_1798511/dataset.json"),
    ("flux", "icepack/PHOTOPROD/H1_1228913/high_energy/dataset.json"),
])
def test_table_input_fingerprint(tmp_path, field, card):
    driver = graniitti_driver.GraniittiDriver()
    cards = [{"datacard": card}]
    dataset, _ = steering.load_dataset(card, cdir=ROOT)
    histogram = next(hist for entry in dataset["sets"] for hist in entry["hist"] if hist.get(field))
    source = ROOT / dataset["datapath"] / histogram[field]
    records = driver._data_input_file_records(cdir=ROOT, datacards=cards)
    assert source.relative_to(ROOT).as_posix() in {record["path"] for record in records}
    table = tmp_path / source.name
    table.write_bytes(source.read_bytes())
    histogram[field] = str(table)
    settings = dict(cdir=ROOT, datacards=cards, datasets=[dataset])
    before = driver._data_input_file_records(**settings)
    assert table in {(ROOT / record["path"]).resolve() for record in before}
    table.write_bytes(table.read_bytes() + b"\n")
    assert before != driver._data_input_file_records(**settings)


# Changes to companion tables must invalidate saved luminosity covariances
@pytest.mark.parametrize("card,filename", [
    ("icepack/UPC/PHOTOPROD/CMS_2648536/jpsi_coherent/dataset.json", "Totalcovariancematrix.json"),
    ("icepack/ELASTIC/STAR_1791591/dataset.json", "Cross-Sections.json"),
    ("icepack/ELASTIC/TOTEM_1220862/dataset.json", "Table4.json"),
    ("icepack/UPC/GAMMA/ATLAS_1811464/gammagamma/dataset.json", "Table11.json"),
])
def test_companion_table_fingerprint(tmp_path, card, filename):
    driver = graniitti_driver.GraniittiDriver()
    dataset, _ = steering.load_dataset(card, cdir=ROOT)
    folder = ROOT / dataset["datapath"]
    for source in folder.glob("*.json"):
        (tmp_path / source.name).write_bytes(source.read_bytes())
    dataset["datapath"] = str(tmp_path)
    settings = dict(cdir=ROOT, datacards=[{"datacard": card}], datasets=[dataset])
    before = driver._data_input_file_records(**settings)
    source = tmp_path / filename
    assert source in {(ROOT / record["path"]).resolve() for record in before}
    source.write_bytes(source.read_bytes() + b"\n")
    assert before != driver._data_input_file_records(**settings)
