# Check immutable literature inputs used by the icepacks
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from core.io import readers, steering
from core.plot import plot
from pyHepMC3 import HepMC3 as h

from icepack.SOFTCEP.integrated.CMS_2017 import cuts as cms_cuts
from tests.physics.validation.test_icepacks import INTEGRATED_DATASETS

ROOT = Path(__file__).resolve().parents[3]
INTEGRATED_ROOTS = (
    ROOT / "icepack" / "SOFTCEP" / "integrated",
    ROOT / "icepack" / "GAMMA" / "integrated",
)


# Check the full reader and plotting normalization against the published cross section in pb
@pytest.mark.parametrize("path", sorted(path for root in INTEGRATED_ROOTS for path in root.glob("*/dataset.json")))
def test_integrated_publication_xs_plot_units(path):
    dataset, resolved = steering.load_dataset(str(path), cdir=ROOT)
    if not dataset.get("active"):
        pytest.skip("No original HEPData measurement configured")
    for entry in dataset["sets"]:
        if not entry.get("data", True):
            continue
        observables, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
        data, observables = readers.read_hepdata(
            entry, dataset["datapath"], dataset["type"], observables,
            cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved,
        )
        histograms = plot.histhepdata(data, observables)
        for histogram in entry["hist"]:
            measured = data[histogram["obs"]]
            units_pb = {"mb": 1e9, "ub": 1e6, "nb": 1e3, "pb": 1.0, "fb": 1e-3}
            expected = float(measured["y"] @ measured["binwidth"]) * units_pb[measured["unit"]]
            assert histograms[histogram["obs"]]["hdata"].integral() == pytest.approx(expected)


# Count extra stable pions from forward fragmentation only when they enter the CMS fiducial region
@pytest.mark.parametrize("extra_pt,accepted", [(None, True), (0.05, True), (0.3, False)])
def test_cms_dipion_stable_particle_selection(extra_pt, accepted):
    event = h.GenEvent()
    vertex = h.GenVertex()
    central_energy = 0.0
    for pt in ([0.3] if extra_pt is None else [0.3, extra_pt]):
        for sign in (1, -1):
            energy = math.hypot(pt, 0.1396)
            central_energy += energy
            vertex.add_particle_out(h.GenParticle(h.FourVector(sign * pt, 0, 0, energy), sign * 211, 1))
    for sign in (1, -1):
        for energy, incoming in ((3500.0, True), (3500.0 - central_energy / 2, False)):
            momentum = h.FourVector(0, 0, sign * math.sqrt(energy**2 - 0.9383**2), energy)
            particle = h.GenParticle(momentum, 2212, 4 if incoming else 1)
            (vertex.add_particle_in if incoming else vertex.add_particle_out)(particle)
    event.add_vertex(vertex)
    assert cms_cuts.cut_func(SimpleNamespace(evt=event, cut_param=cms_cuts.cut_param)) is accepted


# Check the integrated report preserves each configured reader value and its physical units
@pytest.mark.parametrize("path", INTEGRATED_DATASETS)
def test_integrated_table_selected_reference(path):
    from tests.physics.validation.test_icepacks import integrated_table_row

    dataset, resolved = steering.load_dataset(str(path), cdir=ROOT)
    for index, entry in enumerate(dataset["sets"]):
        histogram = next(item for item in entry["hist"] if item["obs"] == "cross_section")
        obs, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=ROOT)
        data, _ = readers.read_hepdata(entry, dataset["datapath"], dataset["type"], obs,
                  cdir=str(ROOT), reader=dataset["reader"], dataset_path=resolved)
        measured = data["cross_section"]
        value_pb = float(measured["y"][0]) * float(histogram["scale"]) * 1e12
        sample = dict(data_integral=value_pb, integral=value_pb, integral_error=0.01 * value_pb)
        report = {"sets": [dict(observables=[dict(observable="cross_section", samples=[sample])])] * len(dataset["sets"])}
        row = integrated_table_row(dataset, report, index, resolved)
        scale = {"mb": 1e9, "ub": 1e6, "nb": 1e3, "pb": 1.0, "fb": 1e-3}[row["unit"]]
        assert row["cross_section"] * scale == pytest.approx(value_pb)
        assert row["ratio_to_data"] == pytest.approx(1.0)
        for side in ("up", "down"):
            error = math.sqrt(sum(row[field][side]**2 for field in
                    ("statistical_error", "systematic_error", "luminosity_error") if row[field] is not None))
            assert error * scale == pytest.approx(measured[f"y_err_{side}"][0] * histogram["scale"] * 1e12)
