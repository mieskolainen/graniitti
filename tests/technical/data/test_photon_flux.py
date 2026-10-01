# Original H1 photon flux normalization through iceplot and icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering
from core.plot import plot
from core.tune.drivers.graniitti.driver import GraniittiDriver

from icepack.PHOTOPROD._common.hera_reader import photon_flux

ROOT = Path(__file__).resolve().parents[3]
# Integrated ep measurements require no photon flux conversion
CARDS = [path for path in sorted((ROOT / "icepack/PHOTOPROD").glob("H1_*/**/dataset.json"))
         if steering.load_dataset(str(path), cdir=ROOT)[0]["type"] != "INTEGRATED_HEPDATA"]


# Extract fluxes independently from the original JSON field labels
def original_flux(table, rows):
    columns = [i for i, header in enumerate(table["headers"][:table["x_count"]]) if "phi" in header["name"].lower()]
    if columns:
        return np.asarray([float(row["x"][columns[0]]["value"]) for row in rows])
    return next((float(entries[0]["value"]) for name, entries in table["qualifiers"].items() if "phi" in name.lower()), None)


# Reconstruct every published spectrum from ep weights with one photon flux conversion
@pytest.mark.parametrize("card", CARDS, ids=lambda path: str(path.parent.relative_to(ROOT / "icepack")))
def test_photon_flux_normalization(card, tmp_path):
    dataset, resolved = steering.load_dataset(str(card), cdir=ROOT)
    driver = GraniittiDriver()
    driver.init_data(run_name=str(tmp_path / "flux"), datacards=[{"datacard": str(card)}],
                     obs_module="default", cdir=str(ROOT), pickle_dump=False)
    for index, entry in enumerate(dataset["sets"]):
        definitions, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
        spectra, observables = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], definitions,
            cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved)
        for hist in entry["hist"]:
            name = hist["obs"]
            path = ROOT / dataset["datapath"] / hist["file"]
            table = json.loads(path.read_text())
            rows = [table["values"][i] for i in hist.get("rows", range(len(table["values"])))]
            headers = [column["name"] for column in table["headers"][:table["x_count"]]]
            if hist.get("file_filter") in {"elastic", "pd_t", "pd_W"}:
                field = "value" if hist["file_filter"] == "elastic" else "low"
                target = next((i for i, name in enumerate(headers) if name.startswith("$m_Y$")), 1)
                rows = [row for row in rows if field in row["x"][target]]
            if "flux" in hist:
                flux_table = json.loads((path.parent / hist["flux"]).read_text())
                flux = np.sum(original_flux(flux_table, flux_table["values"]))
            else:
                flux = original_flux(table, rows)
            if flux is None:
                # Calculated fluxes are checked independently against the analytic Q2 integral
                flux = 1.0 / np.asarray(observables[name]["mc_scale"])
            for column, header in enumerate(headers):
                if header.startswith("$t$") and any(axis.startswith(r"$m_{\pi\pi}$") for axis in headers):
                    flux = flux * np.asarray([float(row["x"][column]["high"]) - float(row["x"][column]["low"])
                                              for row in rows])
            native = np.asarray([float(row["y"][0]["value"]) for row in rows])
            expected_pb = native * hist["scale"] * 1e12
            np.testing.assert_allclose(spectra[name]["y"], native)
            width = spectra[name]["binwidth"] if hist.get("differential", True) else np.ones_like(native)
            scale = entry.get("mc_scale", 1.0)
            weights = expected_pb * width * flux / scale
            mcdata = {"data": {name: spectra[name]["x"]}, "weights": weights, "xsection_pb": weights.sum()}
            for obs in (observables, driver.obs[0][index]):
                np.testing.assert_allclose(obs[name]["mc_scale"], 1.0 / flux)
                result = plot.histmc(mcdata, {name: obs[name]}, scale=scale)[name]["hdata"]
                np.testing.assert_allclose(result.counts_scaled, expected_pb)
                np.testing.assert_allclose(result.errs_scaled, np.abs(expected_pb))
            np.testing.assert_allclose(driver.data[0][index][name]["hdata"].counts_scaled, expected_pb)


# Reject an extra card factor when the reader already supplies the photon flux conversion
def test_duplicate_photon_flux_scale():
    card = "icepack/PHOTOPROD/H1_1798511/dataset.json"
    dataset, resolved = steering.load_dataset(card, cdir=ROOT)
    entry = dataset["sets"][0]
    entry["hist"][0]["mc_scale"] = 2.0
    obs, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
    with pytest.raises(ValueError, match="both reader and card"):
        readers.read_hepdata(entry, dataset["datapath"], dataset["type"], obs,
            cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved)


# Reject invalid original flux fields during input parsing
@pytest.mark.parametrize("value", ["0", "-1", "nan", "inf"])
@pytest.mark.parametrize("source", ["column", "qualifier"])
def test_invalid_photon_flux(value, source):
    folder = ROOT / "HEPData/PHOTOPROD"
    path = folder / ("HEPData-ins1228913-v1-json/Table1.json" if source == "column"
                     else "HEPData-ins1798511-v1-json/Table17.json")
    table = json.loads(path.read_text())
    point = table["values"][0]["x"][1] if source == "column" else table["qualifiers"][r"$\Phi_{\gamma/e}$"][0]
    point["value"] = value
    with pytest.raises(ValueError):
        photon_flux(table)
