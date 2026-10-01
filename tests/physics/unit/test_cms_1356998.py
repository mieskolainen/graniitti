# Test the CMS 7 TeV diffraction generation and Rivet data routing
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path
from types import SimpleNamespace

import pytest
from core.io.hepdata_reader import load_dataset_reader
from pyHepMC3 import HepMC3 as h

from icepack.MINBIAS.CMS_1356998 import obs

ROOT = Path(__file__).resolve().parents[3]
DATASET_PATH = ROOT / "icepack" / "MINBIAS" / "CMS_1356998" / "dataset.json"
REFERENCE_PATH = ROOT / "HEPData" / "Rivet_1356998" / "CMS_2015_I1356998.yoda"


# Construct conserving pp -> X p or X Y events with physical X,Y -> p gamma decays
def diffractive_event(final, status=2):
    energy = 7000.0
    masses = [mass for _, mass in final]
    momentum = math.sqrt((energy**2 - sum(masses)**2) * (energy**2 - (masses[0] - masses[1])**2)) / (2 * energy)
    vertex = h.GenVertex()
    decays = []
    for sign, (pid, mass) in zip((1, -1), final, strict=True):
        vertex.add_particle_in(h.GenParticle(
            h.FourVector(0.0, 0.0, sign * math.sqrt((energy / 2)**2 - 0.9382720813**2), energy / 2), 2212, 4
        ))
        energy_out = math.hypot(momentum, mass)
        particle = h.GenParticle(h.FourVector(0.0, 0.0, sign * momentum, energy_out), pid, status if pid == 90210 else 1)
        vertex.add_particle_out(particle)
        if pid == 90210:
            decay = h.GenVertex()
            decay.add_particle_in(particle)
            q = (mass**2 - 0.9382720813**2) / (2 * mass)
            for px, rest_energy, daughter_pid in ((q, mass - q, 2212), (-q, q, 22)):
                decay.add_particle_out(h.GenParticle(
                    h.FourVector(px, 0.0, sign * momentum * rest_energy / mass, energy_out * rest_energy / mass),
                    daughter_pid, 1,
                ))
            decays.append(decay)
    event = h.GenEvent()
    event.add_vertex(vertex)
    for decay in decays:
        event.add_vertex(decay)
    return SimpleNamespace(evt=event, pid=[90210])


# Check the fixed CMS beam orientation through real SD and DD event projections
@pytest.mark.parametrize('project,final', [
    (obs.proj_log10_xi_sd, [(90210, 100.0), (2212, 0.9382720813)]),
    (obs.proj_log10_xi_sd, [(90210, 100.0), (90210, 2.0)]),
    (obs.proj_log10_xi_dd_forward, [(90210, 100.0), (90210, 5.0)]),
])
def test_forward_gap_fixed_beam_orientation(project, final):
    assert project(diffractive_event(final)) == pytest.approx(math.log10(100.0**2 / 7000.0**2))
    assert project(diffractive_event(final[::-1])) == obs.UNDERFLOW


# Require masses reconstructed from stable particles to be independent of generator intermediate status
def test_gap_reco_stable_four_vectors():
    final = [(90210, 50.0), (90210, 20.0)]
    expected_gap = -math.log(50.0**2 * 20.0**2 / (7000.0**2 * obs.PROTON_MASS_GEV**2))
    for status in (2, 81):
        event = diffractive_event(final, status=status)
        systems = obs.diffractive_systems(event)
        assert [system.m for system in systems] == pytest.approx([50.0, 20.0])
        assert obs.proj_delta_eta_dd_central(event) == pytest.approx(expected_gap)
        assert (systems[0] + systems[1]).m2 == pytest.approx(obs.collision_s(event))


# Reject YODA measurements through the production reader API
@pytest.mark.parametrize("region", ["sd_forward", "dd_forward", "dd_central"])
def test_measurement_unavailable(region):
    reader = load_dataset_reader(str(DATASET_PATH), cdir=str(ROOT))
    with pytest.raises(ValueError, match="HEPData JSON"):
        reader.read(str(REFERENCE_PATH), region=region)
