# Test ATLAS UPC dimuon cuts and observables with physical HepMC3 momenta
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.kinematics.vec4 import vec4

from icepack.UPC.GAMMA.ATLAS_1832628.mumu import cuts, cuts_mass, cuts_rapidity, obs
from tests.technical.support.hepmc import decay_event


# Construct a massive dimuon pair with optional resolved photon recoil
def pair_event(selection=cuts, *, mass=None, rapidity=0.0, pt=0.0, muon_pt=None, photon=False, phi=0.0):
    mass = mass or 1.5 * cuts.cut_param["mass_min"]
    energy = mass / 2.0
    momentum = math.sqrt(energy**2 - 0.105658**2)
    transverse = momentum if muon_pt is None else muon_pt
    longitudinal = math.sqrt(max(0.0, momentum**2 - transverse**2))
    mt = math.hypot(mass, pt)
    system = vec4(pt, 0.0, mt * math.sinh(rapidity), mt * math.cosh(rapidity))
    particles = []
    for sign in (-1, 1):
        p4 = vec4(0.0, sign * transverse, sign * longitudinal, energy)
        p4.boost(system, sign=1)
        p4.rotateZ(phi)
        particles.append((sign * 13, p4, 1))
    if photon:
        recoil = vec4(-pt, 0.0, 0.0, abs(pt))
        recoil.rotateZ(phi)
        particles.append((22, recoil, 1))
    return decay_event(particles, pid=[13, -13], cut_param=selection.cut_param)


# Check every pair cut on both sides of its fiducial boundary
@pytest.mark.parametrize("selection,variable,bound,direction", [
    (cuts, "mass", "mass_min", 1),
    (cuts, "mass", "mass_max", -1),
    (cuts, "pt", "system_pt_max", -1),
    (cuts, "muon_pt", "muon_pt_min", 1),
    (cuts_mass, "rapidity", "abs_rap_max", -1),
    (cuts_rapidity, "mass", "mass_max", -1),
])
@pytest.mark.parametrize("side", [-1, 1])
def test_fiducial_boundaries(selection, variable, bound, direction, side):
    event = pair_event(selection, **{variable: selection.cut_param[bound] + side * 1.0e-5})
    assert bool(selection.cut_func(event)) is (side == direction)


# Check both muon pseudorapidity limits independently of the pair rapidity
@pytest.mark.parametrize("sign", [-1, 1])
@pytest.mark.parametrize("side", [-1, 1])
def test_muon_eta(sign, side):
    eta = sign * (cuts.cut_param["muon_abs_eta_max"] + side * 1.0e-5)
    pt = cuts.cut_param["mass_min"]
    particles = []
    for charge in (-1, 1):
        p4 = vec4()
        p4.setPtEtaPhiM(pt, eta, 0.0 if charge > 0 else math.pi, 0.105658)
        particles.append((charge * 13, p4, 1))
    event = decay_event(particles, cut_param=cuts.cut_param)
    assert bool(cuts.cut_func(event)) is (side < 0)


# A photon can balance the central system while the final muon pair fails its cut
@pytest.mark.parametrize("phi", [0.0, 0.73, -2.4])
@pytest.mark.parametrize("rapidity", [-0.4, 0.4])
def test_radiated_pair_cuts_and_observables(phi, rapidity):
    pt = 1.5 * cuts.cut_param["system_pt_max"]
    event = pair_event(pt=pt, photon=True, phi=phi, rapidity=rapidity)
    assert cuts.dimuon(event).pt == pytest.approx(pt)
    assert obs.proj_abs_rap(event) == pytest.approx(abs(rapidity))
    assert obs.proj_mass(event) == pytest.approx(1.5 * cuts.cut_param["mass_min"])
    mass = obs.proj_mass(event)
    acoplanarity = 2.0 * math.atan(pt / math.sqrt(mass**2 - 4.0 * 0.105658**2)) / math.pi
    assert obs.proj_acoplanarity(event) == pytest.approx(acoplanarity)
    assert not cuts.cut_func(event)
