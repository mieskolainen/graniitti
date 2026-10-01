# Validate H1 proton dissociation data, photon flux conversion and event reconstruction
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core.io import readers, steering
from core.kinematics.vec4 import vec4
from core.plot import plot
from core.stats.uncertainty import source_covariance
from pyHepMC3 import HepMC3 as hepmc3

from icepack.PHOTOPROD._common import diss_reader, dissociative
from tests.technical.support.hepdata import assert_original

DATA = Path(__file__).resolve().parents[3] / "HEPData" / "PHOTOPROD"
JPSI = DATA / "HEPData-ins1228913-v1-json"
RHO = DATA / "HEPData-ins1798511-v1-json"
# Check published gamma p values and the inverse photon flux conversion through the real histogram API
@pytest.mark.parametrize(
    ("directory", "table", "card"),
    [(JPSI, 2, "H1_1228913/dissociative_high_energy"),
     (JPSI, 4, "H1_1228913/dissociative_low_energy"),
     (RHO, 7, "H1_1798511/dissociative")],
)
def test_published_w_flux_integral(directory, table, card):
    path = directory / f"Table{table}.json"
    raw = json.loads(path.read_text())
    rows = diss_reader.rho_rows(raw) if directory == RHO else raw["values"]
    flux_column = 3 if directory == RHO else 1
    flux = np.asarray([float(row["x"][flux_column]["value"]) for row in rows])
    sigma = np.asarray([float(row["y"][0]["value"]) for row in rows])
    root = DATA.parents[1]
    dataset, resolved = steering.load_dataset(str(root / "icepack/PHOTOPROD" / card / "dataset.json"), cdir=root)
    entry = next(entry for entry in dataset["sets"] if entry["hist"][0]["obs"] == "w")
    observables, _ = steering.load_observables(entry["obs"], dataset_path=resolved, cdir=root)
    data, observables = readers.read_hepdata(
        dataset=entry, datapath=dataset["datapath"], datatype=dataset["type"], all_obs=observables,
        cdir=str(root), reader=dataset["reader"], dataset_path=resolved,
    )
    result = data["w"]
    assert result["y"] == pytest.approx(sigma)
    expected_pb = sigma * result["scale"] * 1e12
    ep_weights = flux * expected_pb
    mc = plot.histmc(
        mcdata={"data": {"w": result["x"]}, "weights": ep_weights, "xsection_pb": ep_weights.sum()},
        obs=observables,
    )["w"]["hdata"]
    assert mc.counts_scaled == pytest.approx(expected_pb)
    assert np.dot(mc.counts_scaled, flux) == pytest.approx(ep_weights.sum())
    if directory == JPSI:
        error = np.asarray([float(row["y"][0]["errors"][0]["symerror"]) for row in rows])
        assert result["y_err_up"] == pytest.approx(error)
        covariance = sum(source_covariance(source) for source in result["uncertainties"])
        assert result["y_err"] == pytest.approx(np.sqrt(np.diag(covariance)))
        np.testing.assert_allclose(covariance, np.diag(error**2))


# Check the differential t tails and retain all rho correlated uncertainty sources
def test_dissociative_t_spectra_and_rho_cov():
    for name in ("Table6.json", "Table8.json"):
        path = JPSI / name
        assert_original(diss_reader.read(str(path)), path)
    path = RHO / "Table10.json"
    rows = json.loads(path.read_text())["values"]
    selection = [index for index, row in enumerate(rows) if "low" in row["x"][1]]
    rho = diss_reader.read(str(path), file_filter="pd_t", hist={"covariance": "Table10_statCorr.json"})
    assert_original(rho, path, rows=selection)
    assert len(rho["uncertainties"]) == len(rows[selection[0]]["y"][0]["errors"])
    stat = source_covariance(rho["uncertainties"][0])
    indices = {int(rows[index]["x"][2]["value"]): i for i, index in enumerate(selection)}
    for row in json.loads((RHO / "Table10_statCorr.json").read_text())["values"]:
        first, second = (int(axis["value"]) for axis in row["x"])
        if first in indices and second in indices:
            i, j = indices[first], indices[second]
            assert stat[i, j] / np.sqrt(stat[i, i] * stat[j, j]) == pytest.approx(float(row["y"][0]["value"]))
    covariance = sum(source_covariance(source) for source in rho["uncertainties"])
    assert np.linalg.eigvalsh(covariance).min() > 0.0


# Construct an exact gamma p two-body configuration in either beam direction
def photo_event(direction, fragment, target_mass, *, mass=0.77, decay_pdg=211, decay_mass=0.13957):
    event = hepmc3.GenEvent(hepmc3.Units.GEV, hepmc3.Units.MM)
    w, mp = 50.0, 0.938272
    photon_p = (w * w - mp * mp) / (2.0 * w)
    energy = (w * w + mass * mass - target_mass * target_mass) / (2.0 * w)
    momentum = np.sqrt(energy * energy - mass * mass)
    central = np.asarray([momentum * np.sin(0.02), 0.0, direction * momentum * np.cos(0.02), energy])
    target = np.asarray([0.0, 0.0, 0.0, w]) - central
    photon = np.asarray([0.0, 0.0, direction * photon_p, photon_p])

    # Build HepMC particles with their physical production and fragmentation vertices
    def particle(vector, pid, status):
        return hepmc3.GenParticle(hepmc3.FourVector(*vector), pid, status)

    # Construct on-shell two-body daughters and boost them with the physical parent
    def daughters(parent, masses, pids):
        total = vec4(*parent)
        energy = (total.m2 + masses[0] ** 2 - masses[1] ** 2) / (2.0 * total.m)
        momentum = np.sqrt(energy**2 - masses[0] ** 2)
        first = vec4(0.0, momentum, 0.0, energy)
        second = vec4(0.0, -momentum, 0.0, total.m - energy)
        first.boost(b=total, sign=1)
        second.boost(b=total, sign=1)
        return [particle([p.x, p.y, p.z, p.t], pid, 1) for p, pid in zip([first, second], pids, strict=True)]

    lepton = np.asarray([0.0, 0.0, direction * 10.0, 10.0])
    production = hepmc3.GenVertex()
    production.add_particle_in(particle(lepton + photon, -11, 4))
    production.add_particle_in(particle([0.0, 0.0, -direction * photon_p, w - photon_p], 2212, 4))
    production.add_particle_out(particle(lepton, -11, 1))
    resonance = particle(central, 90, 81)
    excitation = particle(target, 90210, 81 if fragment else 1)
    production.add_particle_out(resonance)
    production.add_particle_out(excitation)
    event.add_vertex(production)
    decay = hepmc3.GenVertex()
    decay.add_particle_in(resonance)
    for daughter in daughters(central, [decay_mass, decay_mass], [decay_pdg, -decay_pdg]):
        decay.add_particle_out(daughter)
    event.add_vertex(decay)
    if fragment:
        forward = hepmc3.GenVertex()
        forward.add_particle_in(excitation)
        for daughter in daughters(target, [0.939565, 0.13957], [2112, 211]):
            forward.add_particle_out(daughter)
        event.add_vertex(forward)
    cuts = {"W": [20.0, 80.0], "Q2_MAX": 2.5, "ABS_T_MAX": 1.5, "MY_MAX": 10.0, "M": [0.2792, 1.53]}
    transfer = photon - central
    expected_t = abs(transfer[3] ** 2 - np.dot(transfer[:3], transfer[:3]))
    return SimpleNamespace(evt=event, pid=[decay_pdg, -decay_pdg], cut_param=cuts), expected_t


# Reconstruct the whole proton system independently of fragmentation and beam ordering
@pytest.mark.parametrize("direction", [-1, 1])
@pytest.mark.parametrize("fragment", [False, True])
@pytest.mark.parametrize("target_mass", [2.0, 12.0])
def test_dissociative_reconstruction(direction, fragment, target_mass):
    event, transfer = photo_event(direction, fragment, target_mass)
    values = dissociative.kinematics(event)
    assert values["w"] == pytest.approx(50.0)
    assert values["abs_t"] == pytest.approx(transfer)
    assert values["mass"] == pytest.approx(0.77)
    assert values["my2"] == pytest.approx(target_mass**2)
    assert dissociative.accepted(event) == (target_mass < 10.0)
