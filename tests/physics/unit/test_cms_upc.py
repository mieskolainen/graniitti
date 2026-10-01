# Test CMS UPC exclusivity with conserving HepMC3 final states
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.kinematics.vec4 import vec4

from icepack.UPC.GAMMA.CMS_2861858.ee import cuts as ee
from icepack.UPC.GAMMA.CMS_2861858.gammagamma import cuts as yy
from tests.technical.support.hepmc import decay_event


# Construct an on-shell final particle at a specified transverse momentum and direction
def particle(pid, pt, eta=0.0, phi=0.0, mass=0.0):
    p4 = vec4()
    p4.setPtEtaPhiM(pt, eta, phi, mass)
    return pid, p4, 1


# Construct a selected pair and extra particles with an exactly conserving parent vertex
def event(selection, extra=(), eta=0.0, phi=0.0):
    pids = [11, -11] if selection is ee else [22, 22]
    pt = 2 * selection.cut_param['electron_et_min' if selection is ee else 'photon_pt_min']
    mass = 0.000511 if selection is ee else 0.0
    pair = [particle(pid, pt, eta, phi + index * math.pi, mass) for index, pid in enumerate(pids)]
    return decay_event(pair + list(extra), pid=pids, cut_param=selection.cut_param)


# Exercise both measured final states through their actual cuts.py modules
@pytest.fixture(params=[ee, yy], ids=['ee', 'gammagamma'])
def selection(request):
    return request.param


# Reject additional tracks by charge, including antiparticles and nuclear fragments
@pytest.mark.parametrize('pid', [13, -13, 211, -321, 2212, -3222, 1000020040])
@pytest.mark.parametrize('above', [False, True])
def test_charged_veto(selection, pid, above):
    pt = selection.cut_param['charged_pt_min'] * (1.01 if above else 0.99)
    sample = event(selection, [particle(pid, pt, phi=math.pi / 2, mass=0.14)])
    assert bool(selection.cut_func(sample)) is not above


# Apply central track and neutral calorimeter pseudorapidity acceptance
@pytest.mark.parametrize('pid,threshold,eta_limit', [
    (211, 'charged_pt_min', 'charged_abs_eta_max'), (22, 'neutral_et_min', 'neutral_abs_eta_max'),
])
@pytest.mark.parametrize('inside', [False, True])
def test_veto_acceptance(selection, pid, threshold, eta_limit, inside):
    param = selection.cut_param
    extra = particle(pid, param[threshold] * 1.1, param[eta_limit] * (0.99 if inside else 1.01), math.pi / 2)
    assert bool(selection.cut_func(event(selection, [extra]))) is not inside


# Reject isolated photons above the neutral transverse-energy threshold
@pytest.mark.parametrize('above', [False, True])
def test_neutral_threshold(selection, above):
    pt = selection.cut_param['neutral_et_min'] * (1.01 if above else 0.99)
    assert bool(selection.cut_func(event(selection, [particle(22, pt, phi=math.pi / 2)]))) is not above


# Respect candidate ECAL windows in the barrel and endcap, including wrapped azimuth
@pytest.mark.parametrize('endcap', [False, True])
@pytest.mark.parametrize('axis', ['eta', 'phi'])
@pytest.mark.parametrize('inside', [False, True])
def test_ecal_window(selection, endcap, axis, inside):
    param = selection.cut_param
    eta = param['ecal_barrel_abs_eta_max'] * (1.2 if endcap else 0.5)
    phi = math.pi - 0.02
    width = param['electron_dphi_max'][endcap] if selection is ee else param['photon_dphi_max']
    shift = (0.99 if inside else 1.01) * (param['neutral_deta_max'] if axis == 'eta' else width)
    extra = particle(22, 1.5 * param['neutral_et_min'], eta + (shift if axis == 'eta' else 0),
                     phi + (shift if axis == 'phi' else 0))
    assert bool(selection.cut_func(event(selection, [extra], eta=eta, phi=phi))) is inside


# Neutral hadrons remain vetoed inside the ECAL windows while neutrinos are invisible
@pytest.mark.parametrize('pid,visible', [(130, True), (2112, True), (-2112, True), (12, False), (-14, False), (16, False)])
def test_neutral_hadrons(selection, pid, visible):
    pt = 1.5 * selection.cut_param['neutral_et_min']
    extra = particle(pid, pt, phi=0.01, mass=0.94 if visible else 0)
    assert bool(selection.cut_func(event(selection, [extra]))) is not visible


# Apply the published ZDC OR, retaining one-sided neutron activity
@pytest.mark.parametrize('counts', [(0, 0), (3, 0), (0, 3), (2, 3), (3, 2), (3, 3)])
def test_zdc_neutrons(selection, counts):
    param = selection.cut_param
    energy, eta, mass = param['zdc_energy_max'] / 2.5, param['zdc_abs_eta_min'] * 1.1, 0.939565
    pt = math.sqrt(energy**2 - mass**2) / math.cosh(eta)
    extra = [particle(2112, pt, side * eta, mass=mass)
             for side, count in zip((-1, 1), counts, strict=True) for _ in range(count)]
    assert bool(selection.cut_func(event(selection, extra))) is (min(counts) < 3)


# Resolve the ZDC energy threshold independently of neutron multiplicity and eta coverage
@pytest.mark.parametrize('inside', [False, True])
@pytest.mark.parametrize('above', [False, True])
def test_zdc_energy_acceptance(selection, inside, above):
    param = selection.cut_param
    energy = param['zdc_energy_max'] * (1.01 if above else 0.99)
    eta, mass = param['zdc_abs_eta_min'] * (1.01 if inside else 0.99), 0.939565
    pt = math.sqrt(energy**2 - mass**2) / math.cosh(eta)
    extra = [particle(2112, pt, side * eta, mass=mass) for side in (-1, 1)]
    assert bool(selection.cut_func(event(selection, extra))) is not (inside and above)
