# Tests for cascade decay angular observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import pytest
from core.analysis import obs
from core.kinematics.vec4 import vec4

from tests.technical.support.hepmc import decay_event


# Verify angle conversions against radians and degrees
def test_angle_conversions():
    assert obs.rad2deg(math.pi) == pytest.approx(180.0)
    assert obs.deg2rad(180.0) == pytest.approx(math.pi)


# Verify the startup preflight can compile and execute a representative kernel
def test_numba_runtime_preflight():
    assert obs.preflight_numba_runtime() is None


# Construct on-shell daughters from the two-body phase-space momentum
def two_body_rest(mother_mass, m1, m2, theta, phi):
    x = mother_mass**2
    y = m1**2
    z = m2**2
    pnorm = math.sqrt(max(0.0, x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z))) / (
        2.0 * mother_mass
    )
    px = pnorm * math.sin(theta) * math.cos(phi)
    py = pnorm * math.sin(theta) * math.sin(phi)
    pz = pnorm * math.cos(theta)
    return (
        vec4(px, py, pz, math.sqrt(m1 * m1 + pnorm * pnorm)),
        vec4(-px, -py, -pz, math.sqrt(m2 * m2 + pnorm * pnorm)),
    )


# Boost a copy of a rest-frame daughter to the supplied frame
def boost_from_rest(p, system):
    out = p.copy()
    out.boost(b=system, sign=1)
    return out


# Construct a four-pion cascade with unequal parent masses
def make_four_pion_event():
    mX = 2.4
    mA = 0.77
    mB = 0.82
    mpi = 0.13957
    X_lab = vec4(0.42, -0.28, 0.64, math.sqrt(mX * mX + 0.42**2 + (-0.28) ** 2 + 0.64**2))

    A_in_X, B_in_X = two_body_rest(mX, mA, mB, 0.82, -0.47)
    A_lab = boost_from_rest(A_in_X, X_lab)
    B_lab = boost_from_rest(B_in_X, X_lab)

    Ap_in_A, Am_in_A = two_body_rest(mA, mpi, mpi, 1.13, 0.71)
    Bp_in_B, Bm_in_B = two_body_rest(mB, mpi, mpi, 0.96, -1.21)

    Ap_lab = boost_from_rest(Ap_in_A, A_lab)
    Am_lab = boost_from_rest(Am_in_A, A_lab)
    Bp_lab = boost_from_rest(Bp_in_B, B_lab)
    Bm_lab = boost_from_rest(Bm_in_B, B_lab)

    return decay_event(
        [
            (211, Ap_lab, 1),
            (211, Bp_lab, 1),
            (-211, Am_lab, 1),
            (-211, Bm_lab, 1),
        ]
    )


# Derive helicity cosines from Lorentz scalars, without invoking a frame transform
@pytest.mark.parametrize("transformed", [False, True])
def test_4body_helicity_matches_invariant_projection(transformed):
    event = make_four_pion_event()
    if transformed:
        records = obs.proj_central_particle_records(event)
        for record in records:
            record["p4"].rotateY(0.63)
            record["p4"].boost(vec4(0.3, -0.2, 0.4, 1.0), sign=1)
        event = decay_event([(r["pid"], r["p4"], 1) for r in records])
    pions_in_hx, systems, particles = obs.proj_1D_4body_helicity(event)
    X = sum(sum(pair) for pair in particles)
    assert [system.m for system in systems] == pytest.approx([0.77, 0.82])
    for index, (A, pair) in enumerate(zip(systems, particles, strict=True)):
        a = pair[0]
        # In the A rest frame the helicity axis is opposite the grandmother momentum
        numerator = a * X - (a * A) * (X * A) / A.m2
        denominator = math.sqrt(((a * A)**2 / A.m2 - a.m2) * ((X * A)**2 / A.m2 - X.m2))
        expected = numerator / denominator
        assert pions_in_hx[index][0].costheta == pytest.approx(expected, abs=1e-12)
        assert pions_in_hx[index][1].costheta == pytest.approx(-expected, abs=1e-12)
        project = obs.proj_1D_4body_cos1 if index == 0 else obs.proj_1D_4body_cos2
        assert project(event) == pytest.approx(expected, abs=1e-12)
        if index == 0:
            assert obs.proj_1D_4body_cos1_dotprod(event) == pytest.approx(expected, abs=1e-12)


# Propagate a mass hypothesis through every paired observable without cache contamination
def test_4body_pairing_mass_argument():
    from core.io.steering import select_histogram_observables

    event = make_four_pion_event()
    functions = [obs.proj_1D_4body_M_A, obs.proj_1D_4body_M_B, obs.proj_1D_4body_DeltaM_AB,
                 obs.proj_1D_4body_cos1, obs.proj_1D_4body_cos2, obs.proj_1D_4body_Deltacos_12,
                 obs.proj_1D_4body_cos1_dotprod, obs.proj_1D_4body_phi12]
    results = []
    for target in [0.8, 3.0, 0.8]:
        helicities, systems, _ = obs.proj_1D_4body_helicity(event, target_m=target)
        a, b = [pair[0] for pair in helicities]
        phi = math.degrees(math.atan2(math.sin(a.phi + b.phi), math.cos(a.phi + b.phi)))
        expected = [systems[0].m, systems[1].m, systems[0].m - systems[1].m,
                    a.costheta, b.costheta, a.costheta - b.costheta, a.costheta, phi]
        definitions = {str(i): {'func': function} for i, function in enumerate(functions)}
        selected = select_histogram_observables({'hist': [
            {'obs': name, 'args': {'target_m': target}} for name in definitions]}, definitions)
        values = [item['func'](event) for item in selected.values()]
        assert values == pytest.approx(expected, abs=1e-12)
        results.append(values)
    assert results[0] == pytest.approx(results[2])
    assert results[0][:2] != pytest.approx(results[1][:2])
