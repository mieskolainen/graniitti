# Test ATLAS 7 TeV lepton cuts using real HepMC3 momenta
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.kinematics.vec4 import vec4

from icepack.GAMMA.ATLAS_1377585.ee import cuts as ee
from icepack.GAMMA.ATLAS_1377585.mumu import cuts as mumu
from tests.technical.support.hepmc import decay_event


# Test both physical lepton masses through the same published selection cases
@pytest.mark.parametrize('cuts,pid,mass,pt_min', [(ee, 11, 0.000511, 12.0), (mumu, 13, 0.105658, 10.0)])
@pytest.mark.parametrize('system_pt,eta,first_pt,accepted', [
    (1.499, 0.0, 20.0, True), (1.5, 0.0, 20.0, False),
    (0.5, 2.4 - 1e-10, 20.0, True), (0.5, 2.4 + 1e-10, 20.0, False),
    (0.0, 0.0, None, False),
    (0.0, 0.0, 37.5, False), (0.0, 0.0, 53.5, True),
])
def test_lepton_fiducial_selection(cuts, pid, mass, pt_min, system_pt, eta, first_pt, accepted):
    first_pt = pt_min if first_pt is None else first_pt
    particles = []
    for pdg, px in ((pid, first_pt), (-pid, system_pt - first_pt)):
        pz = px * math.sinh(eta)
        particles.append((pdg, vec4(px, 0.0, pz, math.sqrt(px**2 + pz**2 + mass**2)), 1))
    assert bool(cuts.cut_func(decay_event(particles, cut_param=cuts.cut_param))) is accepted


# Require opposite charges independently of the lepton record order
@pytest.mark.parametrize("cuts,pid,mass", [(ee, 11, 0.000511), (mumu, 13, 0.105658)])
@pytest.mark.parametrize("charges", [(1, -1), (-1, 1), (1, 1), (-1, -1)])
def test_lepton_pair_charge(cuts, pid, mass, charges):
    particles = [(sign * pid, vec4(px, 0.0, 0.0, math.hypot(px, mass)), 1)
                 for sign, px in zip(charges, (20.0, -20.0), strict=True)]
    if charges[0] == charges[1]:
        particles.extend((charges[0] * 211, vec4(0.0, 0.0, 0.0, 0.139570), 1) for _ in range(2))
    event = decay_event(particles, cut_param=cuts.cut_param)
    assert bool(cuts.cut_func(event)) is (charges[0] != charges[1])


