# Tests for exclusive dimuon rapidity gap and void observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pytest
from core.io import steering
from pyHepMC3 import HepMC3 as h
from scipy.integrate import quad

from tests.technical.support.hepmc import make_particle

ROOT = Path(__file__).resolve().parents[3]

from core.analysis import obs
from core.kinematics.vec4 import vec4

HARDPOM_Z = ROOT / "icepack" / "HARDPOM" / "z_with_pythia"
diffz_cuts = steering.load_python_module(
    str(HARDPOM_Z / "cuts.py"),
)
diffz_obs = steering.load_python_module(
    str(HARDPOM_Z / "obs.py"),
)


# Compare the boundary integral with direct quadrature, including nearly uniform weights
@pytest.mark.parametrize("beta", [1.0e-16, 0.1, 2.0])
@pytest.mark.parametrize("pt", [[0.7, 2.0, 0.3], [1.0e20, 2.0, 0.3]])
def test_gap_flow_matches_boundary_integral(beta, pt):
    eta = np.array([-1.7, 0.2, 1.4])
    pt = np.array(pt)
    low, high, q0 = -2.5, 2.5, 0.8
    # Integrate the defining Laplace weight without subtracting accumulated penalties
    numerator = quad(lambda y: np.exp(beta * (y - high) - pt[eta > y].sum() / q0),
                     low, high, points=eta, epsabs=1.0e-12)[0]
    denominator = quad(lambda y: np.exp(beta * (y - high)), low, high)[0]
    actual = diffz_obs.energy_gap_flow_score(eta, pt, 1, low, high, 0.0, beta, q0)
    assert actual == pytest.approx(numerator / denominator, rel=1.0e-12, abs=1.0e-14)


# Compute a massive four-vector from transverse coordinates
def make_p4(pt, eta, phi, mass=0.105658):
    px = pt * math.cos(phi)
    py = pt * math.sin(phi)
    pz = pt * math.sinh(eta)
    energy = math.sqrt(px * px + py * py + pz * pz + mass * mass)
    return vec4(px, py, pz, energy)


# Compute one particle record compatible with the diffZ projector helpers
def make_record(pid, pt, eta, phi=0.0, is_final=True, status=None):
    if status is None:
        status = 1 if is_final else 2
    return {
        "pid": pid,
        "status": status,
        "p4": make_p4(pt=pt, eta=eta, phi=phi, mass={13: 0.105658, 211: 0.139570, 321: 0.493677, 2212: 0.93827, 22: 0.0}[abs(pid)]),
    }


# Compute one incoming proton beam record on the requested side
def make_beam_record(side):
    energy = 6500.0
    proton_mass = 0.93827
    pz = side * math.sqrt(energy * energy - proton_mass * proton_mass)
    return {
        "pid": 2212,
        "status": 4,
        "p4": vec4(0.0, 0.0, pz, energy),
    }


# Compute one stable forward proton with explicit longitudinal momentum
def make_forward_proton_record(side, pz_abs, pt=0.2):
    proton_mass = 0.93827
    energy = math.sqrt(pz_abs * pz_abs + pt * pt + proton_mass * proton_mass)
    pz = side * pz_abs
    return {
        "pid": 2212,
        "status": 1,
        "p4": vec4(pt, 0.0, pz, energy),
    }


# Build a conserving HepMC event with explicit beams, final states and recoil
def make_event(records):
    event = h.GenEvent()
    vertex = h.GenVertex()
    recoil = vec4()
    for record in records:
        momentum = record['p4']
        particle = make_particle(*momentum, record['pid'], record['status'])
        if record['status'] == 4:
            vertex.add_particle_in(particle)
            recoil += momentum
        else:
            vertex.add_particle_out(particle)
            recoil -= momentum
    assert recoil.e > 0 and recoil.m2 > 0
    vertex.add_particle_out(make_particle(*recoil, 990, 2))
    event.add_vertex(vertex)
    return SimpleNamespace(evt=event, cut_param=dict(diffz_cuts.cut_param))


# Compute a selected Z pair plus controlled charged and ignored particles
def make_selected_z_event(extra=()):
    selected_minus = make_record(13, 46.0, 0.2, 0.0)
    selected_plus = make_record(-13, 45.0, -0.2, math.pi - 0.02)
    extra_hard_plus = make_record(-13, 70.0, 2.3, 1.0)
    extra_soft_minus = make_record(13, 5.0, -2.0, -0.5)
    right_pion = make_record(211, 2.0, 2.0, 0.4)
    neutral = make_record(22, 100.0, 2.2, 0.7)
    nonfinal = make_record(211, 100.0, -2.2, -0.7, is_final=False)
    zero_pt = make_record(321, 0.0, 1.5, 0.0)
    outside = make_record(2212, 1.0, 5.0, 0.0)
    return make_event(
        [
            selected_minus,
            selected_plus,
            extra_hard_plus,
            extra_soft_minus,
            right_pion,
            neutral,
            nonfinal,
            zero_pt,
            outside,
            make_beam_record(-1),
            make_beam_record(1),
            *extra,
        ]
    )


# Compute a symmetry-transformed copy of one HepMC event
def transform_event(event, phi_rotation=0.0, reflect_eta=False, exchange_muon_charge=False):
    cosine = math.cos(phi_rotation)
    sine = math.sin(phi_rotation)
    transformed = []
    for record in obs.proj_event_particles(event):
        if record['pid'] == 990:
            continue
        clone = dict(record)
        momentum = record["p4"]
        px = cosine * momentum.px - sine * momentum.py
        py = sine * momentum.px + cosine * momentum.py
        pz = -momentum.pz if reflect_eta else momentum.pz
        clone["p4"] = vec4(px, py, pz, momentum.e)
        if exchange_muon_charge and abs(record["pid"]) == 13:
            clone["pid"] = -record["pid"]
        transformed.append(clone)
    output = make_event(transformed)
    output.cut_param = dict(event.cut_param)
    return output


# Verify an empty accepted side has unit energy gap flow score
def test_empty_side_has_unit_score():
    score = diffz_obs.energy_gap_flow_score(
        eta=np.asarray([-2.0]),
        pt=np.asarray([1.0]),
        side=1,
        eta_min=0.0,
        eta_max=4.9,
        pt_min=0.0,
        beta=1.0,
        q0=1.0,
    )
    assert score == pytest.approx(1.0)


# Verify the sorted implementation against the one-particle Eq. 35 result
def test_one_particle_score_matches_eq35():
    eta = 2.0
    pt = 1.5
    eta_max = 4.9
    edge_weight = math.exp(eta - eta_max)
    low_weight = math.exp(-eta_max)
    normalization = 1.0 - low_weight
    expected = (math.exp(-pt) * (edge_weight - low_weight) + 1.0 - edge_weight) / normalization

    score = diffz_obs.energy_gap_flow_score(
        eta=np.asarray([eta]),
        pt=np.asarray([pt]),
        side=1,
        eta_min=0.0,
        eta_max=eta_max,
        pt_min=0.0,
        beta=1.0,
        q0=1.0,
    )
    assert score == pytest.approx(expected, abs=1.0e-14)


# Verify the multiplicity score against the one-particle Laplace-void result
def test_multiplicity_matches_eq24():
    eta = 2.0
    pt = 1.5
    p0 = 1.0
    eta_max = 4.9
    edge_weight = math.exp(eta - eta_max)
    low_weight = math.exp(-eta_max)
    normalization = 1.0 - low_weight
    survival = p0 / (p0 + pt)
    expected = (survival * (edge_weight - low_weight) + 1.0 - edge_weight) / normalization

    score = diffz_obs.multiplicity_gap_flow_score(
        eta=np.asarray([eta]),
        pt=np.asarray([pt]),
        side=1,
        eta_min=0.0,
        eta_max=eta_max,
        pt_min=0.0,
        beta=1.0,
        p0=p0,
    )
    assert score == pytest.approx(expected, abs=1.0e-14)


# Verify the multiplicity score against a two-particle boundary integral
def test_two_multiplicity_matches_sorted_integral():
    eta = np.asarray([0.7, 2.1])
    pt = np.asarray([0.5, 3.0])
    eta_min = 0.0
    eta_max = 3.0
    beta = 0.4
    p0 = 1.0
    boundaries = np.asarray([eta_min, eta[0], eta[1], eta_max])
    z = np.exp(beta * (boundaries - eta_max))
    survival = np.asarray([(p0 / (p0 + pt[0])) * (p0 / (p0 + pt[1])), p0 / (p0 + pt[1]), 1.0])
    normalization = -math.expm1(-beta * (eta_max - eta_min))
    expected = np.dot(survival, np.diff(z)) / normalization

    observed = diffz_obs.multiplicity_gap_flow_score(eta, pt, 1, eta_min, eta_max, 0.0, beta, p0)
    assert observed == pytest.approx(expected, abs=1.0e-14)


# Verify collinear fragmentation increases the multiplicity penalty
def test_multiplicity_score_responds_particle_count():
    single = diffz_obs.multiplicity_gap_flow_score(
        np.asarray([2.0]), np.asarray([2.0]), 1, 0.0, 4.9, 0.0, 1.0, 1.0
    )
    split = diffz_obs.multiplicity_gap_flow_score(
        np.asarray([2.0, 2.0]), np.asarray([1.0, 1.0]), 1, 0.0, 4.9, 0.0, 1.0, 1.0
    )
    assert split < single


# Verify a particle suppresses only the side containing its pseudorapidity
def test_particle_affects_only_rapidity_side():
    eta = np.asarray([2.0])
    pt = np.asarray([3.0])
    left = diffz_obs.energy_gap_flow_score(eta, pt, -1, 0.0, 4.9, 0.0, 1.0, 1.0)
    right = diffz_obs.energy_gap_flow_score(eta, pt, 1, 0.0, 4.9, 0.0, 1.0, 1.0)
    assert left == pytest.approx(1.0)
    assert right < 1.0


# Verify eta reflection interchanges the left and right scores
def test_eta_reflection_swaps_side_scores():
    eta = np.asarray([-3.0, -0.7, 1.2, 2.4])
    pt = np.asarray([0.4, 2.0, 1.3, 3.1])
    left = diffz_obs.energy_gap_flow_score(eta, pt, -1, 0.0, 4.9, 0.0, 1.0, 1.0)
    right = diffz_obs.energy_gap_flow_score(eta, pt, 1, 0.0, 4.9, 0.0, 1.0, 1.0)
    reflected_left = diffz_obs.energy_gap_flow_score(-eta, pt, -1, 0.0, 4.9, 0.0, 1.0, 1.0)
    reflected_right = diffz_obs.energy_gap_flow_score(-eta, pt, 1, 0.0, 4.9, 0.0, 1.0, 1.0)
    assert reflected_left == pytest.approx(right)
    assert reflected_right == pytest.approx(left)


# Verify increasing side transverse momentum suppresses the energy gap flow score
def test_score_decreases_with_particle_pt():
    eta = np.asarray([3.0])
    low_pt = diffz_obs.energy_gap_flow_score(eta, np.asarray([0.2]), 1, 0.0, 4.9, 0.0, 1.0, 1.0)
    high_pt = diffz_obs.energy_gap_flow_score(eta, np.asarray([5.0]), 1, 0.0, 4.9, 0.0, 1.0, 1.0)
    assert high_pt < low_pt


# Verify the kernel applies the strict charged-particle transverse-momentum cut
def test_score_applies_pt_min():
    eta = np.asarray([2.0])
    pt = np.asarray([0.1])
    rejected = diffz_obs.energy_gap_flow_score(eta, pt, 1, 0.0, 4.9, 0.1, 1.0, 1.0)
    accepted = diffz_obs.energy_gap_flow_score(eta, pt, 1, 0.0, 4.9, 0.09, 1.0, 1.0)
    assert rejected == pytest.approx(1.0)
    assert accepted < 1.0


# Verify the full oriented scan resolves a four-unit forward gap
def test_full_scan_resolves_four_unit_gap():
    eta = np.asarray([-1.5])
    pt = np.asarray([100.0])
    score = diffz_obs.energy_gap_flow_score(
        eta=eta,
        pt=pt,
        side=1,
        eta_min=-2.5,
        eta_max=2.5,
        pt_min=0.1,
        beta=1.0,
        q0=1.0,
    )
    effective_gap = diffz_obs.energy_gap_flow_effective_gap(score, -2.5, 2.5, 1.0)
    assert effective_gap == pytest.approx(4.0, abs=1.0e-12)


# Verify the effective gap transform maps score boundaries to the rapidity span
def test_effective_gap_transform_boundaries():
    assert diffz_obs.energy_gap_flow_effective_gap(0.0, -2.5, 2.5, 0.1) == pytest.approx(0.0)
    assert diffz_obs.energy_gap_flow_effective_gap(1.0, -2.5, 2.5, 0.1) == pytest.approx(5.0)


# Verify the effective-gap transform is monotonic in the raw score
def test_effective_gap_transform_monotonic():
    scores = np.linspace(0.0, 1.0, 11)
    gaps = np.asarray(
        [diffz_obs.energy_gap_flow_effective_gap(score, -2.5, 2.5, 0.1) for score in scores]
    )
    assert np.all(np.diff(gaps) > 0.0)


# Verify each full-range side receives particles from both hemispheres
def test_full_scan_both_hemispheres_per_side():
    eta = np.asarray([2.0])
    pt = np.asarray([3.0])
    left = diffz_obs.energy_gap_flow_score(eta, pt, -1, -2.5, 2.5, 0.1, 1.0, 1.0)
    right = diffz_obs.energy_gap_flow_score(eta, pt, 1, -2.5, 2.5, 0.1, 1.0, 1.0)
    assert right < left < 1.0


# Verify each forward gap uses the nearest accepted particle across the event
def test_forward_gap_acceptance():
    eta = np.asarray([-1.0, 0.4, 1.7, 2.2])
    pt = np.full(eta.size, 1.0)
    left = diffz_obs.forward_gap_size(eta, pt, -1, 0.0, 2.5, 0.1)
    right = diffz_obs.forward_gap_size(eta, pt, 1, 0.0, 2.5, 0.1)
    cross_central = diffz_obs.forward_gap_size(eta, pt, -1, 1.2, 2.5, 0.1)

    assert left == pytest.approx(1.5)
    assert right == pytest.approx(0.3)
    assert cross_central == pytest.approx(1.3)


# Verify an event without accepted charged particles spans the full detector
def test_empty_forward_gap_spans_full_acceptance():
    eta = np.asarray([], dtype=np.float64)
    pt = np.asarray([], dtype=np.float64)
    assert diffz_obs.forward_gap_size(eta, pt, -1, -2.5, 2.5, 0.1) == pytest.approx(5.0)
    assert diffz_obs.forward_gap_size(eta, pt, 1, -2.5, 2.5, 0.1) == pytest.approx(5.0)


# Verify the conventional gap applies its charged-particle fiducial cuts
def test_forward_gap_applies_pt_eta_thresholds():
    eta = np.asarray([2.8, 2.4, 1.5])
    pt = np.asarray([100.0, 0.1, 0.6])

    assert diffz_obs.forward_gap_size(eta, pt, 1, -2.5, 2.5, 0.5) == pytest.approx(1.0)
    assert diffz_obs.forward_gap_size(eta, pt, 1, -2.5, 2.5, 0.05) == pytest.approx(0.1)


# Verify the largest gap includes internal intervals and both fiducial edges
def test_largest_eta_gap_includes_all_intervals():
    eta = np.asarray([-1.0, -0.2, 2.0])
    pt = np.ones(eta.size)
    observed = diffz_obs.largest_pseudorapidity_gap(eta, pt, -2.5, 2.5, 0.5)
    assert observed == pytest.approx(2.2)


# Verify an empty accepted event has one gap across the full fiducial span
def test_empty_largest_eta_gap_spans_acceptance():
    eta = np.asarray([-3.0, 2.4])
    pt = np.asarray([10.0, 0.5])
    observed = diffz_obs.largest_pseudorapidity_gap(eta, pt, -2.5, 2.5, 0.5)
    assert observed == pytest.approx(5.0)


# Verify eta reflection leaves the largest symmetric-fiducial gap unchanged
def test_largest_eta_gap_reflection_sym():
    eta = np.asarray([-2.1, -0.4, 1.2, 2.3])
    pt = np.asarray([1.0, 2.0, 3.0, 4.0])
    original = diffz_obs.largest_pseudorapidity_gap(eta, pt, -2.5, 2.5, 0.5)
    reflected = diffz_obs.largest_pseudorapidity_gap(-eta, pt, -2.5, 2.5, 0.5)
    assert reflected == pytest.approx(original)


# Verify the two-boundary Laplace score against its one-particle closed form
def test_particle_energy_gap():
    eta_min = -2.5
    eta_max = 2.5
    eta = 0.7
    pt = 1.8
    q0 = 0.9
    d0 = eta - eta_min
    d1 = eta_max - eta
    length = eta_max - eta_min
    expected = (d0**2 + d1**2 + 2.0 * d0 * d1 * math.exp(-pt / q0)) / length**2

    observed = diffz_obs.two_boundary_energy_gap_score(
        np.asarray([eta]), np.asarray([pt]), eta_min, eta_max, 0.0, q0
    )
    assert observed == pytest.approx(expected, abs=1.0e-14)


# Verify empty accepted activity gives unit two-boundary Laplace score
def test_empty_energy_gap_unity():
    observed = diffz_obs.two_boundary_energy_gap_score(
        np.asarray([], dtype=np.float64), np.asarray([], dtype=np.float64), -2.5, 2.5, 0.3, 2.5
    )
    assert observed == pytest.approx(1.0)


# Verify reflection leaves the symmetric two-boundary score unchanged
def test_energy_gap_reflection_sym():
    eta = np.asarray([-2.1, -0.4, 1.2, 2.3])
    pt = np.asarray([1.0, 2.0, 3.0, 4.0])
    original = diffz_obs.two_boundary_energy_gap_score(eta, pt, -2.5, 2.5, 0.5, 2.0)
    reflected = diffz_obs.two_boundary_energy_gap_score(-eta, pt, -2.5, 2.5, 0.5, 2.0)
    assert reflected == pytest.approx(original)


# Verify increasing accepted transverse momentum suppresses the Laplace score
def test_energy_gap_decreases_with_pt():
    eta = np.asarray([0.4])
    low_pt = diffz_obs.two_boundary_energy_gap_score(eta, np.asarray([0.5]), -2.5, 2.5, 0.1, 1.0)
    high_pt = diffz_obs.two_boundary_energy_gap_score(eta, np.asarray([5.0]), -2.5, 2.5, 0.1, 1.0)
    assert high_pt < low_pt


# Verify additive transverse energy makes same-eta collinear splitting invariant
def test_energy_gap_collinear_splitting_invariant():
    single = diffz_obs.two_boundary_energy_gap_score(
        np.asarray([0.4]), np.asarray([2.0]), -2.5, 2.5, 0.1, 1.0
    )
    split = diffz_obs.two_boundary_energy_gap_score(
        np.asarray([0.4, 0.4]), np.asarray([0.75, 1.25]), -2.5, 2.5, 0.1, 1.0
    )
    assert split == pytest.approx(single, abs=1.0e-14)


# Verify hard and transparent Laplace limits recover gap concentration and unity
def test_energy_gap_limits():
    eta = np.asarray([-1.0, -0.2, 2.0])
    pt = np.ones(eta.size)
    hard_expected = (1.5**2 + 0.8**2 + 2.2**2 + 0.5**2) / 5.0**2
    hard = diffz_obs.two_boundary_energy_gap_score(eta, pt, -2.5, 2.5, 0.5, 1.0e-9)
    transparent = diffz_obs.two_boundary_energy_gap_score(eta, pt, -2.5, 2.5, 0.5, 1.0e12)
    assert hard == pytest.approx(hard_expected, abs=1.0e-14)
    assert transparent == pytest.approx(1.0, abs=1.0e-11)


# Verify the two-boundary kernel applies its strict transverse-momentum cut
def test_energy_gap_applies_pt_min():
    eta = np.asarray([0.4])
    pt = np.asarray([0.5])
    rejected = diffz_obs.two_boundary_energy_gap_score(eta, pt, -2.5, 2.5, 0.5, 1.0)
    accepted = diffz_obs.two_boundary_energy_gap_score(eta, pt, -2.5, 2.5, 0.49, 1.0)
    assert rejected == pytest.approx(1.0)
    assert accepted < 1.0


# Verify a nonpositive Laplace transverse-energy scale is rejected
@pytest.mark.parametrize("q0", [0.0, -1.0])
def test_energy_gap_invalid_q0(q0):
    with pytest.raises(ValueError, match="q0"):
        diffz_obs.two_boundary_energy_gap_score(
            np.asarray([0.4]), np.asarray([1.0]), -2.5, 2.5, 0.1, q0
        )


# Verify the exact light-cone loss and selection of the most beam-like proton
@pytest.mark.parametrize('momenta', [(5000.,), (5000., 6200.)])
def test_proton_loss_matches_definition(momenta):
    beam = make_beam_record(1)['p4']
    protons = [make_forward_proton_record(1, pz_abs=pz, pt=.2)['p4'] for pz in momenta]
    expected = 1. - max(p.e + p.pz for p in protons) / (beam.e + beam.pz)
    observed = diffz_obs.leading_proton_lightcone_loss(
        eta=np.array([p.eta for p in protons]), pt=np.array([p.pt for p in protons]),
        pz=np.array([p.pz for p in protons]), energy=np.array([p.e for p in protons]),
        side=1, beam_lightcone=beam.e + beam.pz, eta_min=6., pt_max=2.)
    assert observed == pytest.approx(expected, abs=1.e-14)


# Verify a side without an accepted forward proton has unit momentum loss
def test_proton_loss_missing_candidate_unity():
    observed = diffz_obs.leading_proton_lightcone_loss(
        eta=np.asarray([5.9, 7.0]),
        pt=np.asarray([0.2, 2.0]),
        pz=np.asarray([100.0, 100.0]),
        energy=np.asarray([101.0, 101.0]),
        side=1,
        beam_lightcone=13000.0,
        eta_min=6.0,
        pt_max=2.0,
    )
    assert observed == pytest.approx(1.0)


# Verify side reflection preserves the forward-proton loss
def test_proton_loss_reflection_sym():
    eta = np.asarray([8.0])
    pt = np.asarray([0.3])
    pz = np.asarray([6200.0])
    energy = np.sqrt(pz * pz + pt * pt + 0.93827**2)
    right = diffz_obs.leading_proton_lightcone_loss(eta, pt, pz, energy, 1, 13000.0, 6.0, 2.0)
    left = diffz_obs.leading_proton_lightcone_loss(-eta, pt, -pz, energy, -1, 13000.0, 6.0, 2.0)
    assert left == pytest.approx(right)


# Verify the event projector selects the strongest accepted beam-side proton
def test_forward_proton_projector_minimum_side_loss():
    proton = make_forward_proton_record(1, pz_abs=5000.0, pt=0.2)
    event = make_selected_z_event([proton])
    beam = make_beam_record(1)["p4"]
    expected = 1.0 - (proton["p4"].e + proton["p4"].pz) / (beam.e + beam.pz)

    left, right = diffz_obs.forward_proton_lightcone_losses(event)
    assert left == pytest.approx(1.0)
    assert right == pytest.approx(expected)
    assert diffz_obs.proj_xi_forward_proton_min(event) == pytest.approx(expected)
    assert diffz_obs.proj_d_forward_proton(event) == pytest.approx(-math.log10(expected))


# Verify exact Z records are removed while additional muons remain
def test_selected_z_records_removed_extra_muons():
    event = make_selected_z_event()
    pair = diffz_cuts.select_mumu_pair(event)
    assert pair is not None
    assert [record['pid'] for record in pair['records']] == [13, -13]
    assert [record['pt'] for record in pair['records']] == pytest.approx([46.0, 45.0])

    expected_eta = np.asarray([2.3, -2.0, 2.0])
    expected_pt = np.asarray([70.0, 5.0, 2.0])
    expected_left = diffz_obs.energy_gap_flow_score(
        expected_eta,
        expected_pt,
        -1,
        event.cut_param["EGF_CH_ETA_MIN"],
        event.cut_param["EGF_CH_ETA_MAX"],
        event.cut_param["EGF_CH_PT_MIN"],
        event.cut_param["EGF_BETA"],
        event.cut_param["EGF_Q0"],
    )
    expected_right = diffz_obs.energy_gap_flow_score(
        expected_eta,
        expected_pt,
        1,
        event.cut_param["EGF_CH_ETA_MIN"],
        event.cut_param["EGF_CH_ETA_MAX"],
        event.cut_param["EGF_CH_PT_MIN"],
        event.cut_param["EGF_BETA"],
        event.cut_param["EGF_Q0"],
    )
    expected_multiplicity_left = diffz_obs.multiplicity_gap_flow_score(
        expected_eta,
        expected_pt,
        -1,
        event.cut_param["EGF_CH_ETA_MIN"],
        event.cut_param["EGF_CH_ETA_MAX"],
        event.cut_param["EGF_CH_PT_MIN"],
        event.cut_param["EGF_BETA"],
        event.cut_param["MGF_P0"],
    )
    expected_multiplicity_right = diffz_obs.multiplicity_gap_flow_score(
        expected_eta,
        expected_pt,
        1,
        event.cut_param["EGF_CH_ETA_MIN"],
        event.cut_param["EGF_CH_ETA_MAX"],
        event.cut_param["EGF_CH_PT_MIN"],
        event.cut_param["EGF_BETA"],
        event.cut_param["MGF_P0"],
    )

    left, right = diffz_obs.charged_energy_gap_flow_scores(event)
    assert left == pytest.approx(expected_left)
    assert right == pytest.approx(expected_right)
    assert diffz_obs.proj_n_charged_central(event) == 3
    assert diffz_obs.proj_g_egf_max(event) == pytest.approx(max(left, right))

    multiplicity_left, multiplicity_right = diffz_obs.charged_multiplicity_gap_flow_scores(event)
    assert multiplicity_left == pytest.approx(expected_multiplicity_left)
    assert multiplicity_right == pytest.approx(expected_multiplicity_right)
    assert diffz_obs.proj_g_mgf_max(event) == pytest.approx(
        max(multiplicity_left, multiplicity_right)
    )

    effective_left, effective_right = diffz_obs.charged_energy_gap_flow_effective_gaps(event)
    assert effective_left == pytest.approx(diffz_obs.proj_deta_egf_left(event))
    assert effective_right == pytest.approx(diffz_obs.proj_deta_egf_right(event))
    assert diffz_obs.proj_deta_egf_max(event) == pytest.approx(max(effective_left, effective_right))

    gap_left, gap_right = diffz_obs.charged_forward_gap_sizes(event)
    assert gap_left == pytest.approx(0.5)
    assert gap_right == pytest.approx(0.2)
    assert diffz_obs.proj_deta_f_max(event) == pytest.approx(0.5)
    assert diffz_obs.proj_deta_largest(event) == pytest.approx(4.0)
    expected_G2 = diffz_obs.two_boundary_energy_gap_score(
        expected_eta,
        expected_pt,
        event.cut_param["EGF_CH_ETA_MIN"],
        event.cut_param["EGF_CH_ETA_MAX"],
        event.cut_param["EGF_CH_PT_MIN"],
        event.cut_param["G2_EF_Q0"],
    )
    assert diffz_obs.proj_G2_ef(event) == pytest.approx(expected_G2)
    assert diffz_obs.proj_d_forward_proton(event) == pytest.approx(0.0)

    assert diffz_obs.proj_gap_roc_scores(event) == pytest.approx(
        (
            max(left, right),
            max(multiplicity_left, multiplicity_right),
            0.5,
            4.0,
            expected_G2,
            3,
            0.0,
        )
    )


# Verify forward and reversed ROC steering uses opposite class and cut directions
def test_gap_roc_class_direction_steering():
    forward = diffz_obs.obs_gap_roc["roc"]
    reversed_ = diffz_obs.obs_gap_roc_reversed["roc"]

    assert forward["directions"] == [
        "higher",
        "higher",
        "higher",
        "higher",
        "higher",
        "lower",
        "higher",
    ]
    assert reversed_["directions"] == [
        "lower",
        "lower",
        "lower",
        "lower",
        "lower",
        "higher",
        "lower",
    ]
    assert forward["comparisons"][0]["signal"] == 4
    assert forward["comparisons"][0]["background"] == 2
    assert reversed_["comparisons"][0]["signal"] == 2
    assert reversed_["comparisons"][0]["background"] == 4
    assert len(diffz_obs.proj_gap_roc_scores(make_selected_z_event())) == len(
        forward["score_labels"]
    )
    assert forward["score_colors"] == reversed_["score_colors"]
    assert forward["comparisons"][0]["linestyle"] == "-"
    assert forward["comparisons"][1]["linestyle"] == "--"
    assert reversed_["comparisons"][0]["linestyle"] == "-"
    assert reversed_["comparisons"][1]["linestyle"] == "--"
    assert forward["label_template"] == reversed_["label_template"]


# Check gap reconstruction under rotations, beam reflection and muon charge exchange
@pytest.mark.parametrize(('rotation', 'reflection', 'charge'), [(.73, False, False), (0., True, False), (0., False, True)])
def test_gap_event_symmetries(rotation, reflection, charge):
    event = make_selected_z_event()
    transformed = transform_event(event, rotation, reflection, charge)
    for project in (diffz_obs.charged_energy_gap_flow_scores, diffz_obs.charged_multiplicity_gap_flow_scores,
                    diffz_obs.charged_forward_gap_sizes):
        expected = project(event)
        assert project(transformed) == pytest.approx(expected[::-1] if reflection else expected)
    for project in (diffz_obs.proj_G2_ef, diffz_obs.proj_deta_largest, diffz_obs.proj_n_charged_central):
        assert project(transformed) == pytest.approx(project(event))
    assert diffz_cuts.cut_func(transformed) == diffz_cuts.cut_func(event)


# Verify the event multiplicity cut counts all stable charged final states
def test_event_selection_total_charged_multiplicity():
    event = make_selected_z_event()
    event.cut_param["N_CH_MAX"] = 7
    assert diffz_cuts.cut_func(event)

    event.cut_param["N_CH_MAX"] = 6
    assert not diffz_cuts.cut_func(event)


# Verify every gap-flow constituent and kernel parameter comes from cut_param
def test_gap_flow_obs_follow_event_params():
    event = make_selected_z_event()
    event.cut_param["EGF_CH_ETA_MIN"] = 2.1
    event.cut_param["EGF_CH_ETA_MAX"] = 2.4
    event.cut_param["EGF_CH_PT_MIN"] = 3.0
    event.cut_param["EGF_BETA"] = 0.7
    event.cut_param["EGF_Q0"] = 2.0
    event.cut_param["G2_EF_Q0"] = 0.4
    event.cut_param["MGF_P0"] = 1.3

    expected = diffz_obs.energy_gap_flow_score(
        eta=np.asarray([2.3]),
        pt=np.asarray([70.0]),
        side=1,
        eta_min=2.1,
        eta_max=2.4,
        pt_min=3.0,
        beta=0.7,
        q0=2.0,
    )
    left, right = diffz_obs.charged_energy_gap_flow_scores(event)
    assert left == pytest.approx(1.0)
    assert right == pytest.approx(expected)
    assert diffz_obs.proj_n_charged_central(event) == 1

    expected_multiplicity = diffz_obs.multiplicity_gap_flow_score(
        eta=np.asarray([2.3]),
        pt=np.asarray([70.0]),
        side=1,
        eta_min=2.1,
        eta_max=2.4,
        pt_min=3.0,
        beta=0.7,
        p0=1.3,
    )
    multiplicity_left, multiplicity_right = diffz_obs.charged_multiplicity_gap_flow_scores(event)
    assert multiplicity_left == pytest.approx(1.0)
    assert multiplicity_right == pytest.approx(expected_multiplicity)

    expected_effective = diffz_obs.energy_gap_flow_effective_gap(expected, 2.1, 2.4, 0.7)
    assert diffz_obs.proj_deta_egf_left(event) == pytest.approx(0.3)
    assert diffz_obs.proj_deta_egf_right(event) == pytest.approx(expected_effective)

    gap_left, gap_right = diffz_obs.charged_forward_gap_sizes(event)
    assert gap_left == pytest.approx(0.3)
    assert gap_right == pytest.approx(0.1)

    expected_G2 = diffz_obs.two_boundary_energy_gap_score(
        eta=np.asarray([2.3]), pt=np.asarray([70.0]), eta_min=2.1, eta_max=2.4, pt_min=3.0, q0=0.4
    )
    assert diffz_obs.proj_G2_ef(event) == pytest.approx(expected_G2)
