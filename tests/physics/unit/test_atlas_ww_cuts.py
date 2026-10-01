# Check the ATLAS WW prompt lepton selection with conserving HepMC3 events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from types import SimpleNamespace

import pytest
from core.kinematics.vec4 import hepmc2vec4, vec4
from pyHepMC3 import HepMC3 as h

from icepack.GAMMA.integrated.ATLAS_2021 import cuts


# Construct pp -> pp W+W- with on shell W decays and an optional leptonic tau decay
def ww_event(tau=False):
    event = h.GenEvent()
    vertex = h.GenVertex()
    w_mass = 80.4
    for sign in (1, -1):
        for energy, incoming in ((6500.0, True), (6500.0 - w_mass, False)):
            particle = h.GenParticle(h.FourVector(0, 0, sign * math.sqrt(energy**2 - 0.9383**2), energy), 2212, 4 if incoming else 1)
            (vertex.add_particle_in if incoming else vertex.add_particle_out)(particle)
    event.add_vertex(vertex)
    for pid, mass, axis in ((-15 if tau else -11, 1.7769 if tau else 0.000511, 0), (13, 0.1057, 1)):
        boson = h.GenParticle(h.FourVector(0, 0, 0, w_mass), 24 if pid < 0 else -24, 2)
        vertex.add_particle_out(boson)
        decay = h.GenVertex()
        decay.add_particle_in(boson)
        q = (w_mass**2 - mass**2) / (2 * w_mass)
        momentum = vec4(q if axis == 0 else 0, q if axis == 1 else 0, 0, w_mass - q)
        lepton = h.GenParticle(h.FourVector(*momentum), pid, 2 if abs(pid) == 15 else 1)
        decay.add_particle_out(lepton)
        recoil = vec4(0, 0, 0, w_mass) - momentum
        decay.add_particle_out(h.GenParticle(h.FourVector(*recoil), 16 if pid == -15 else (12 if pid == -11 else -14), 1))
        event.add_vertex(decay)
        if abs(pid) == 15:
            tau_decay = h.GenVertex()
            tau_decay.add_particle_in(lepton)
            q = (mass**2 - 0.000511**2) / (2 * mass)
            for daughter_pid, rest in ((-11, vec4(q, 0, 0, mass - q)), (12, vec4(-q / 2, 0, 0, q / 2)), (-16, vec4(-q / 2, 0, 0, q / 2))):
                rest.boost(momentum, sign=1)
                tau_decay.add_particle_out(h.GenParticle(h.FourVector(*rest), daughter_pid, 1))
            event.add_vertex(tau_decay)
    return SimpleNamespace(evt=event, cut_param=cuts.cut_param)


# Accept direct W leptons despite the incoming proton ancestry and exclude W -> tau feed down
@pytest.mark.parametrize("tau", (False, True))
def test_prompt_ww_selection_proton_beams(tau):
    event = ww_event(tau)
    for vertex in event.evt.vertices():
        incoming = sum((hepmc2vec4(particle.momentum()) for particle in vertex.particles_in()), vec4())
        outgoing = sum((hepmc2vec4(particle.momentum()) for particle in vertex.particles_out()), vec4())
        assert tuple(incoming) == pytest.approx(tuple(outgoing), abs=1e-8)
    assert cuts.cut_func(event) is not tau
