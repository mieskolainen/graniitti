# Test LHCb charmonium and bottomonium cuts with real HepMC3 decay events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.kinematics.vec4 import vec4

from icepack.PHOTOPROD.LHCb_1277076 import cuts as tev7
from icepack.PHOTOPROD.LHCb_1373746 import cuts as upsilon
from icepack.PHOTOPROD.LHCb_2825384 import cuts as tev13
from tests.technical.support.hepmc import decay_event


# Apply the same published acceptance checks to all three measurements
@pytest.fixture(params=[tev7, upsilon, tev13], ids=['7TeV', '7and8TeV', '13TeV'])
def cuts(request):
    return request.param


# Construct massive particles with consistent pseudorapidity and four-momentum
def particle(pid, eta, pt=1.0, status=1):
    mass = 0.105658 if abs(pid) == 13 else 0.13957039
    pz = pt * math.sinh(eta)
    return pid, vec4(pt, 0.0, pz, math.sqrt(pt**2 + pz**2 + mass**2)), status


# Accept final opposite-sign muons and ignore nonfinal muons or unrelated pions
@pytest.mark.parametrize('extra', [[], [particle(13, 3.0, status=2), particle(211, 3.0)]])
def test_accepts_final_dimuon_pair(cuts, extra):
    particles = [particle(13, 2.1), particle(-13, 4.4)] + extra
    assert cuts.fiducial_cut(decay_event(particles, pid=[13, -13], cut_param=cuts.cut_param))


# Reject missing, additional and same-sign final muons
@pytest.mark.parametrize('pids', [[], [13], [13, -13, 13], [13, 13], [-13, -13]])
def test_rejects_wrong_muon_multiplicity_or_charge(cuts, pids):
    event = decay_event([particle(pid, 3.0) for pid in pids], pid=[13, -13], cut_param=cuts.cut_param)
    assert not cuts.fiducial_cut(event)


# Resolve both sides of each eta boundary without relying on inverse-function roundoff
@pytest.mark.parametrize('eta,accepted', [(2.0 - 1e-10, False), (2.0 + 1e-10, True),
                                        (4.5 - 1e-10, True), (4.5 + 1e-10, False)])
@pytest.mark.parametrize('pid', [13, -13])
def test_muon_eta_boundaries(cuts, eta, accepted, pid):
    event = decay_event([particle(pid, eta), particle(-pid, 3.0)], cut_param=cuts.cut_param)
    assert bool(cuts.fiducial_cut(event)) is accepted


# Test central rapidity independently of muon eta using a boosted massive pair
@pytest.mark.parametrize('rapidity,accepted', [(2.0 - 1e-10, False), (2.0 + 1e-10, True),
                                             (4.5 - 1e-10, True), (4.5 + 1e-10, False)])
def test_central_rapidity_boundaries(cuts, rapidity, accepted):
    particles = []
    for pid, px in ((13, 1.0), (-13, -1.0)):
        momentum = vec4(px, 0.0, 0.0, math.hypot(px, 0.105658))
        momentum.boost(vec4(0.0, 0.0, math.sinh(rapidity), math.cosh(rapidity)), sign=1)
        particles.append((pid, momentum, 1))
    event = decay_event(particles, cut_param=cuts.cut_param)
    assert bool(cuts.central_system_cut(event)) is accepted
    if not accepted:
        assert not cuts.fiducial_cut(event)


# Corrected measurements have no trigger pT threshold and separate central from muon acceptance
@pytest.mark.parametrize('etas,pt,fiducial', [((2.1, 4.4), 0.1, True),
                                            ((3.0, 3.5), 0.1, True), ((1.0, 5.0), 1.0, False)])
def test_corrected_acceptance(cuts, etas, pt, fiducial):
    event = decay_event([particle(13, etas[0], pt), particle(-13, etas[1], pt)], cut_param=cuts.cut_param)
    assert cuts.central_system_cut(event)
    assert bool(cuts.fiducial_cut(event)) is fiducial
