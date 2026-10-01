# Checks for the UPC HEPData readers and immutable source tables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pyjson5 as json5
import pytest
from core.io.hepdata_reader import load_dataset_reader
from core.stats.uncertainty import source_covariance

from tests.technical.support.hepdata import assert_original, original

ROOT = Path(__file__).resolve().parents[3]
UPC_DATA = ROOT / "HEPData" / "UPC"

ATLAS_CARD = "icepack/UPC/GAMMA/ATLAS_1832628/mumu/dataset.json"
ATLAS_LBYL_CARD = "icepack/UPC/GAMMA/ATLAS_1811464/gammagamma/dataset.json"
ALICE_COHERENT_CARD = "icepack/UPC/PHOTOPROD/ALICE_1840600/jpsi_coherent/dataset.json"
ALICE_INCOHERENT_CARD = "icepack/UPC/PHOTOPROD/ALICE_2658375/jpsi_incoherent/dataset.json"
CMS_CARD = "icepack/UPC/PHOTOPROD/CMS_2648536/jpsi_coherent/dataset.json"
CMS_BW_CARD = "icepack/UPC/GAMMA/CMS_2861858/ee/dataset.json"
CMS_LBYL_CARD = "icepack/UPC/GAMMA/CMS_2861858/gammagamma/dataset.json"
ALICE_PPB_DIMUON_CARD = "icepack/UPC/GAMMA/ALICE_2654315/mumu/dataset.json"
ALICE_PPB_JPSI_CARD = "icepack/UPC/PHOTOPROD/ALICE_2654315/jpsi/dataset.json"
ALICE_EMD_CARD = "icepack/UPC/EMD/ALICE_2149540/dataset.json"
ATLAS_JPSI_CARD = "icepack/UPC/PHOTOPROD/ATLAS_2966819/jpsi_coherent/dataset.json"
CMS_INCOHERENT_CARD = "icepack/UPC/PHOTOPROD/CMS_2899343/jpsi_incoherent/dataset.json"
CMS_PHI_CARD = "icepack/UPC/PHOTOPROD/CMS_2908607/phi_coherent/dataset.json"
CMS_UPSILON_CARD = "icepack/UPC/PHOTOPROD/CMS_3140573/upsilon_coherent/dataset.json"
ALICE_RHO_CARD = "icepack/UPC/PHOTOPROD/ALICE_1782227/rho_coherent/dataset.json"
ALICE_ENERGY_CARD = "icepack/UPC/PHOTOPROD/ALICE_2666011/jpsi_coherent/dataset.json"
ALICE_MID_JPSI_CARD = "icepack/UPC/PHOTOPROD/ALICE_1840601/jpsi_coherent/dataset.json"
ALICE_MID_PSI2S_CARD = "icepack/UPC/PHOTOPROD/ALICE_1840601/psi2s_coherent/dataset.json"

ATLAS_DATA = UPC_DATA / "HEPData-ins1832628-v1-json"
ATLAS_LBYL_DATA = UPC_DATA / "HEPData-ins1811464-v1-json"
ALICE_COHERENT_DATA = UPC_DATA / "HEPData-ins1840600-v1-json"
ALICE_INCOHERENT_DATA = UPC_DATA / "HEPData-ins2658375-v1-json"
CMS_DATA = UPC_DATA / "HEPData-ins2648536-v2-json"
CMS_2861858_DATA = UPC_DATA / "HEPData-ins2861858-v2-json"
ALICE_PPB_DATA = UPC_DATA / "HEPData-ins2654315-v1-json"
ALICE_EMD_DATA = UPC_DATA / "HEPData-ins2149540-v1-json"

# Compute one named uncertainty source
def named_uncertainty(data, name):
    return next(source for source in data["uncertainties"] if source["name"] == name)


# Compute the single histogram steering entry from one UPC set
def histogram(card, index=0):
    dataset = json5.loads((ROOT / card).read_text(encoding="utf-8"))
    return dataset["sets"][index]["hist"][0]


# Check the ATLAS dimuon mass and rapidity tables
@pytest.mark.parametrize("table", ["Table7.json", "Table1.json"])
def test_atlas_1832628_dimuon_tables(table):
    path = ATLAS_DATA / table
    data = load_dataset_reader(ATLAS_CARD, cdir=str(ROOT)).read(str(path))
    assert_original(data, path)


# Check all ALICE pPb dimuon groups and the nonoverlapping J/psi rapidity rows
def test_alice_2654315_ppb_tables():
    dimuon_reader = load_dataset_reader(ALICE_PPB_DIMUON_CARD, cdir=str(ROOT))
    jpsi_reader = load_dataset_reader(ALICE_PPB_JPSI_CARD, cdir=str(ROOT))
    dimuon = dimuon_reader.read(str(ALICE_PPB_DATA / "Table1.json"), file_filter="2.5-4.0")
    hist = histogram(ALICE_PPB_JPSI_CARD)
    exclusive = jpsi_reader.read(str(ALICE_PPB_DATA / "Table2.json"), hist=hist)
    assert_original(dimuon, ALICE_PPB_DATA / "Table1.json", group=0)
    assert_original(exclusive, ALICE_PPB_DATA / "Table2.json", rows=hist["rows"])


# Check the ALICE single EMD neutron multiplicity cross sections
def test_alice_2149540_emd_table():
    path = ALICE_EMD_DATA / "Table2.json"
    neutron = load_dataset_reader(ALICE_EMD_CARD, cdir=str(ROOT)).read(str(path))
    assert_original(neutron, path)
    assert named_uncertainty(neutron, "stat")["category"] == "statistical"
    assert named_uncertainty(neutron, "sys")["category"] == "systematic"


# Check the ATLAS light-by-light differential and integrated measurements
@pytest.mark.parametrize("table", ["Table3.json", "Table5.json", "Table7.json", "Table11.json"])
def test_atlas_1811464_light_by_light_tables(table):
    path = ATLAS_LBYL_DATA / table
    data = load_dataset_reader(ATLAS_LBYL_CARD, cdir=str(ROOT)).read(str(path), file_filter="Measured")
    _, errors = assert_original(data, path)
    expected = errors["total"] if "total" in errors else np.sqrt(sum(error**2 for error in errors.values()))
    for side, index in (("down", 0), ("up", 1)):
        np.testing.assert_allclose(data[f"y_err_{side}"], expected[:, index])
    assert named_uncertainty(data, "stat")["category"] == "statistical"
    if "total" in errors:
        assert named_uncertainty(data, "systematic_from_total")["category"] == "systematic"


# Check all CMS dielectron Figure 4 and light by light Figure 7 tables
@pytest.mark.parametrize("card,index", [(CMS_BW_CARD, index) for index in range(4)] + [(CMS_LBYL_CARD, index) for index in range(2)])
def test_cms_2861858_particle_level_tables(card, index):
    hist = histogram(card, index)
    path = CMS_2861858_DATA / hist["file"]
    data = load_dataset_reader(card, cdir=str(ROOT)).read(str(path), file_filter=hist["file_filter"])
    assert_original(data, path, group=0)


# Check the extracted CMS photonuclear cross section and gluon suppression points
@pytest.mark.parametrize("table", ["Table2.json", "Table3.json"])
def test_cms_2648536_extracted_observables(table):
    path = CMS_DATA / table
    data = load_dataset_reader(CMS_CARD, cdir=str(ROOT)).read(str(path))
    assert_original(data, path)


# Check both ALICE coherent momentum distributions and correlations
def test_alice_1840600_coherent_tables():
    reader = load_dataset_reader(ALICE_COHERENT_CARD, cdir=str(ROOT))
    abs_t = reader.read(str(ALICE_COHERENT_DATA / "Table1.json"))
    pt2 = reader.read(str(ALICE_COHERENT_DATA / "Table2.json"), hist=histogram(ALICE_COHERENT_CARD))
    assert_original(abs_t, ALICE_COHERENT_DATA / "Table1.json")
    assert_original(pt2, ALICE_COHERENT_DATA / "Table2.json")
    np.testing.assert_array_equal(abs_t["bins"], pt2["bins"])
    assert named_uncertainty(abs_t, "experimental_correlated_syst")["correlation"] == "collective"
    assert named_uncertainty(pt2, "stat")["correlation"] == "uncorrelated"
    assert named_uncertainty(pt2, "uncorrelated_syst")["correlation"] == "uncorrelated"
    assert named_uncertainty(pt2, "correlated_syst")["correlation"] == "collective"


# Keep the coherent UPC density in the published production units
def test_alice_1840600_data_normalization():
    reader = load_dataset_reader(ALICE_COHERENT_CARD, cdir=str(ROOT))
    hist = histogram(ALICE_COHERENT_CARD)
    path = ALICE_COHERENT_DATA / hist["file"]
    data = reader.read(str(path), hist=hist)
    _, errors = assert_original(data, path)
    np.testing.assert_allclose(named_uncertainty(data, "stat")["up"], errors["stat."][:, 1])


# Check the ALICE incoherent photonuclear momentum transfer measurement
def test_alice_2658375_incoherent_table():
    path = ALICE_INCOHERENT_DATA / "Table1.json"
    data = load_dataset_reader(ALICE_INCOHERENT_CARD, cdir=str(ROOT)).read(str(path))
    _, errors = assert_original(data, path)
    correlated = named_uncertainty(data, "experimental_correlated_syst")
    assert correlated["correlation"] == "collective"
    np.testing.assert_allclose(correlated["up"], errors["experimental correlated syst."][:, 1])


# Preserve the published ALICE photonuclear density and uncertainties
def test_alice_2658375_data_normalization():
    path = ALICE_INCOHERENT_DATA / "Table1.json"
    data = load_dataset_reader(ALICE_INCOHERENT_CARD, cdir=str(ROOT)).read(str(path), hist=histogram(ALICE_INCOHERENT_CARD))
    _, errors = assert_original(data, path)
    np.testing.assert_allclose(named_uncertainty(data, "experimental_stat")["up"], errors["experimental stat."][:, 1])


# Convert generated kaon pairs to the published phi rapidity density
def test_cms_phi_decay_and_rapidity_norm():
    reader = load_dataset_reader(CMS_PHI_CARD, cdir=str(ROOT))
    hist = histogram(CMS_PHI_CARD)
    path = ROOT / "HEPData/UPC/HEPData-ins2908607-v1-json/Table1.json"
    raw = reader.read(str(path), file_filter=hist["file_filter"])
    data = reader.read(str(path), file_filter=hist["file_filter"], hist=hist)

    assert data["y"] == pytest.approx(raw["y"])
    assert data["y_err"] == pytest.approx(raw["y_err"])


# Check all CMS coherent neutron multiplicity dependent variables
@pytest.mark.parametrize("selection,group", [("AnAn", 0), ("0n0n", 1), ("0nXn", 2), ("XnXn", 3)])
def test_cms_2648536_coherent_neutron_classes(selection, group):
    path = CMS_DATA / "Table1.json"
    data = load_dataset_reader(CMS_CARD, cdir=str(ROOT)).read(str(path), file_filter=selection)
    assert data["group"] == group
    assert_original(data, path, group)


# Convert folded MC dimuon rates to the published J/psi density without rescaling data
def test_cms_2648536_prod_norm():
    dataset = json5.loads((ROOT / CMS_CARD).read_text(encoding="utf-8"))
    reader = load_dataset_reader(CMS_CARD, cdir=str(ROOT))
    event_sets = [entry for entry in dataset["sets"] if entry.get("mc", True)]

    for index in range(len(event_sets)):
        hist = histogram(CMS_CARD, index)
        data = reader.read(
            str(CMS_DATA / "Table1.json"),
            file_filter=hist["file_filter"],
            hist=hist,
        )
        raw = reader.read(str(CMS_DATA / "Table1.json"), file_filter=hist["file_filter"])
        assert data["y"] == pytest.approx(raw["y"])
        assert data["y_err"] == pytest.approx(raw["y_err"])


# Reject ambiguous access to a table with several dependent variables
def test_cms_2648536_requires_neutron_class():
    reader = load_dataset_reader(CMS_CARD, cdir=str(ROOT))
    with pytest.raises(ValueError, match="file_filter is required"):
        reader.read(str(CMS_DATA / "Table1.json"))




# Fold original ALICE uncertainties while leaving the cached native measurement unchanged
def test_alice_rho_hepdata_folding():
    from icepack.UPC._common import upc_reader as shared
    reader = load_dataset_reader(ALICE_RHO_CARD, cdir=str(ROOT))
    filename = str(UPC_DATA / "HEPData-ins1782227-v1-json/Table1.json")
    raw = shared.read(filename, file_filter="0n0n", hist={"rows": [2, 3, 4]})
    native_values, native_bins = raw["y"].copy(), raw["bins"].copy()
    for _ in range(2):
        folded = reader.read(filename, file_filter="0n0n")
        np.testing.assert_allclose(folded["y"], 2.0 * native_values)
        for side in ("up", "down"):
            np.testing.assert_allclose(folded[f"y_err_{side}"], 2.0 * raw[f"y_err_{side}"])
        np.testing.assert_array_equal(raw["y"], native_values)
        np.testing.assert_array_equal(raw["bins"], native_bins)
        assert len(folded["uncertainties"]) == len(raw["uncertainties"])


# Preserve the signed ALICE migration covariance across classes in matching source bins
@pytest.mark.parametrize("partner", ["Table2.json", "Table3.json", "Table4.json"])
def test_alice_neutron_migration_covariance(partner):
    import json

    from core.stats.cov import assemble_data_covariances

    from icepack.UPC._common import upc_reader

    reader = load_dataset_reader(ALICE_ENERGY_CARD, cdir=str(ROOT))
    spectra, layout, coordinates, shifts = [], [], [], []
    offset = 0
    for subset, filename in enumerate(("Table1.json", partner)):
        path = UPC_DATA / "HEPData-ins2666011-v1-json" / filename
        data = reader.read(str(path))
        native = upc_reader.read(str(path))
        for name in ("y", "y_err", "y_err_up", "y_err_down"):
            np.testing.assert_allclose(data[name], native[name])
        valid = np.flatnonzero(data.get("valid", np.ones(len(data["y"]), dtype=bool)))
        spectra.append({"rap": data})
        layout.append(dict(dataset=0, subset=subset, observable="rap", bins=valid.tolist(),
                           start=offset, stop=offset + len(valid)))
        offset += len(valid)
        rows = json.loads(path.read_text())["values"]
        rows.sort(key=lambda row: float(row["x"][0]["low"]))
        coordinates.append([tuple(float(row["x"][0][key]) for key in ("low", "high")) for row in rows])
        shifts.append(np.array([
            (float(error["asymerror"]["plus"]) - float(error["asymerror"]["minus"])) / 2
            for row in rows for error in row["y"][0]["errors"] if error["label"] == "migr"
        ]))
    _, systematic, _, _ = assemble_data_covariances([spectra], layout)
    same_bin = np.array([[left == right for right in coordinates[1]] for left in coordinates[0]])
    expected = np.outer(*shifts) * same_bin
    n = len(shifts[0])
    np.testing.assert_allclose(systematic[:n, n:], expected)
    assert np.all(expected[same_bin] < 0)
    assert np.all(np.linalg.eigvalsh(systematic) >= -np.finfo(float).eps)


# Retain the original total errors while correlating luminosity across ATLAS spectra
@pytest.mark.parametrize("table", ["Table3.json", "Table5.json", "Table7.json", "Table11.json"])
def test_atlas_luminosity_covariance(table):
    reader = load_dataset_reader(ATLAS_LBYL_CARD, cdir=str(ROOT))
    data = reader.read(str(ATLAS_LBYL_DATA / table), file_filter="Measured")
    integrated = reader.read(str(ATLAS_LBYL_DATA / "Table11.json"), file_filter="Measured")
    _, value, errors = original(ATLAS_LBYL_DATA / "Table11.json")
    expected = data["y"] * errors["lumi"][0, 1] / value[0]
    source = named_uncertainty(data, "luminosity")
    np.testing.assert_allclose(source["up"], expected)
    np.testing.assert_allclose(source_covariance(source), np.outer(expected, expected))
    assert source["scope"] == named_uncertainty(integrated, "luminosity")["scope"]
    assert source["effect"] == "multiplicative"
    _, _, quoted = original(ATLAS_LBYL_DATA / table)
    total = quoted["total"] if "total" in quoted else np.sqrt(sum(error**2 for error in quoted.values()))
    np.testing.assert_allclose(data["y_err_up"], total[:, 1])
    np.testing.assert_allclose(data["y_err_down"], total[:, 0])
