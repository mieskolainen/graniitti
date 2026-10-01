# Test dilepton cuts on radiated final states with real HepMC3 events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.analysis import obs
from core.kinematics.vec4 import vec4
from pyHepMC3 import HepMC3 as h

from icepack.GAMMA.ATLAS_1377585.ee import cuts as pp_ee
from icepack.GAMMA.ATLAS_1377585.mumu import cuts as pp_mumu
from icepack.GAMMA.mumu_with_pythia import cuts as shower_mumu
from icepack.HARDPOM.z_with_pythia import cuts as shower_z
from icepack.UPC.GAMMA.ALICE_2654315.mumu import cuts as alice
from icepack.UPC.GAMMA.CMS_2861858.ee import cuts as cms
from tests.technical.support.hepmc import decay_event


# Construct an opposite-sign lepton pair with a photon balancing its transverse recoil
def radiated_event(selection, pid, first_pt, second_pt, *, rapidity=0.0, phi=0.0, photon=0.0):
    mass = 0.000511 if pid == 11 else 0.105658
    particles = []
    boost = vec4(0.0, 0.0, math.sinh(rapidity), math.cosh(rapidity))
    for pdg, px in ((pid, first_pt), (-pid, -second_pt), (22, second_pt - first_pt)):
        if pdg == 22 and math.isclose(px, 0.0, abs_tol=1e-12):
            continue
        p4 = vec4(px, 0.0, 0.0, abs(px) if pdg == 22 else math.hypot(px, mass))
        p4.rotateZ(phi)
        p4.boost(boost, sign=1)
        particles.append((pdg, p4, 1))
    if photon > 0.0:
        p4 = vec4(0.0, 0.0, photon, photon)
        p4.boost(boost, sign=1)
        particles.append((22, p4, 1))
    event = decay_event(particles, pid=[pid, -pid], cut_param=selection.cut_param)
    parent = next(p for p in event.evt.particles() if p.pid() == 990)
    system = sum((p4 for _, p4, _ in particles), vec4())
    production = h.GenVertex()
    production.add_particle_out(parent)
    for sign in (-1, 1):
        p4 = vec4(0.0, 0.0, sign * system.m / 2.0, system.m / 2.0)
        p4.boost(system, sign=1)
        production.add_particle_in(h.GenParticle(h.FourVector(*p4), 22, 4))
    event.evt.add_vertex(production)
    return event


# Reject a radiated pair even when the pair plus photon has zero transverse momentum
@pytest.mark.parametrize('selection,pid,pt_min,pt_max', [
    (pp_ee, 11, 'ELECTRON_PT_MIN', 'SYSTEM_PT_MAX'),
    (pp_mumu, 13, 'MUON_PT_MIN', 'SYSTEM_PT_MAX'),
    (shower_mumu, 13, 'MUON_PT_MIN', 'PT_MUMU_MAX'),
    (cms, 11, 'electron_et_min', 'system_pt_max'),
])
@pytest.mark.parametrize('phi', [0.0, 0.73, -2.4])
@pytest.mark.parametrize('side', [-1, 1])
def test_photon_recoil(selection, pid, pt_min, pt_max, phi, side):
    pt = 2.0 * selection.cut_param[pt_min]
    recoil = selection.cut_param[pt_max] * (1.0 + side * 0.01)
    event = radiated_event(selection, pid, pt + recoil, pt, phi=phi)
    assert obs.proj_1D_Pt(event) == pytest.approx(recoil)
    assert bool(selection.cut_func(event)) is (side < 0)


# ALICE measures the final muon-pair mass even when photons raise the full central mass
@pytest.mark.parametrize('bound,direction', [(0, 1), (1, -1)])
@pytest.mark.parametrize('side', [-1, 1])
def test_alice_mass_migration(bound, direction, side):
    mass = alice.cut_param['mass'][bound] * (1.0 + side * 0.01)
    pt = math.sqrt((mass / 2.0)**2 - 0.105658**2)
    event = radiated_event(alice, 13, pt, pt, rapidity=sum(alice.cut_param['rap']) / 2.0, photon=mass)
    assert obs.proj_1D_M(event) == pytest.approx(mass)
    assert bool(alice.cut_func(event)) is (side == direction)


# Showered Z selection accepts migration from a parent mass above the fiducial window
def test_z_mass_migration():
    mass = shower_z.cut_param['Z_MASS']
    pt = math.sqrt((mass / 2.0)**2 - 0.105658**2)
    energy = 2.0 * (shower_z.cut_param['M_MUMU_MAX'] - mass)
    event = radiated_event(shower_z, 13, pt, pt, photon=energy)
    parent = next(p for p in event.evt.particles() if p.pid() == 990)
    assert parent.momentum().m() > shower_z.cut_param['M_MUMU_MAX']
    assert obs.proj_1D_M(event) == pytest.approx(mass)
    assert shower_z.cut_func(event)
