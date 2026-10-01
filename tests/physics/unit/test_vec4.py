# Tests for the standalone Lorentz-vector implementation
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from decimal import Decimal, localcontext

import numpy as np
import pytest
from core.kinematics.vec4 import vec4


# Check rapidity is invariant under unit changes and odd under longitudinal reflection
@pytest.mark.parametrize("scale", [1.0e-200, 1.0, 8.0e307])
@pytest.mark.parametrize("rapidity", [-1.0, 0.0, 1.0])
def test_rapidity_units_longitudinal_reflection(scale, rapidity):
    vector = vec4(scale, 0.0, scale * np.sinh(rapidity), scale * np.cosh(rapidity))
    assert vector.rapidity == pytest.approx(rapidity, abs=5.0e-14)
    vector.setZ(-vector.pz)
    assert vector.rapidity == pytest.approx(-rapidity, abs=5.0e-14)


# Check exact forward massless limits and the undefined zero four-vector
def test_rapidity_lightlike_limits():
    assert np.isposinf(vec4(0.0, 0.0, 4.0, 4.0).rapidity)
    assert np.isneginf(vec4(0.0, 0.0, -4.0, 4.0).rapidity)
    assert np.isnan(vec4().rapidity)
    assert np.isnan(vec4(0.0, 0.0, 2.0, 1.0).rapidity)


# Check the PDG longitudinal boost law for a massive momentum
def test_rapidity_longitudinal_boost_law():
    vector = vec4()
    vector.setPt2RapPhiM2(2.0, 1.2, 0.4, 0.7)
    boost = vec4(0.0, 0.0, np.sinh(0.6), np.cosh(0.6))
    vector.boost(boost)
    assert vector.rapidity == pytest.approx(0.6)


# Validate construction and arithmetic component by component
@pytest.mark.parametrize('actual,expected', [
    (vec4(1, 2, 3, 4), (1, 2, 3, 4)), (vec4(), (0, 0, 0, 0)),
    (vec4(1, 2, 3, 4) + vec4(4, 3, 2, 1), (5, 5, 5, 5)),
    (vec4(1, 2, 3, 4) - vec4(4, 3, 2, 1), (-3, -1, 1, 3)),
    (vec4(1, 2, 3, 4) * 2, (2, 4, 6, 8)),
])
def test_vec4_arithmetic(actual, expected):
    assert (actual.x, actual.y, actual.z, actual.t) == expected


# Validate in-place scaling updates the private vector storage
def test_vec4_in_place_scale():
    v = vec4(1, 2, 3, 4)
    v.scale(2)
    assert tuple(v) == (2, 4, 6, 8)


# Validate the Minkowski dot product
def test_vec4_dot_product():
    v1 = vec4(1, 0, 0, 1)
    v2 = vec4(0, 1, 0, 1)
    dot_product = v1 * v2
    assert dot_product == 1


# Validate the invariant mass squared
def test_vec4_mass_squared():
    v = vec4(1, 2, 3, 4)
    assert v.m2 == 16 - (1**2 + 2**2 + 3**2)


# Validate the invariant mass
def test_vec4_mass():
    v = vec4(1, 2, 3, 4)
    assert v.m == np.sqrt(v.m2)


# Validate transverse momentum
def test_vec4_transverse_momentum():
    v = vec4(3, 4, 0, 5)
    assert v.pt == 5


# Validate the spatial opening angle
def test_vec4_angle():
    v1 = vec4(1, 0, 0, 1)
    v2 = vec4(0, 1, 0, 1)
    angle = v1.angle(v2)
    assert np.isclose(angle, np.pi / 2)


# Validate pseudorapidity
def test_vec4_eta():
    v = vec4(3, 4, 5, 7)

    pT = np.sqrt(v.x**2 + v.y**2)
    theta = np.arctan2(pT, v.z)
    expected_eta = -np.log(np.tan(theta / 2))

    assert np.isclose(v.eta, expected_eta, atol=1e-6), f"Expected: {expected_eta}, Got: {v.eta}"


# Validate polar angles including axial limits
@pytest.mark.parametrize('components,expected', [
    ((3, 4, 5, 7), np.pi / 4), ((3, 4, 0, 5), np.pi / 2), ((0, 0, 5, 5), 0.0),
])
def test_vec4_theta(components, expected):
    assert vec4(*components).theta == pytest.approx(expected, abs=1e-6)


# Validate boosts in both directions
def test_vec4_boost():
    for sign in [-1, 1]:
        v = vec4(0, 0, 1, 2)
        boost_vector = vec4(0, 0, 0.1, np.sqrt(1.01))

        beta = boost_vector.pz / boost_vector.e
        gamma = 1.0 / np.sqrt(1 - beta**2)

        expected_pz = gamma * (v.pz - sign * beta * v.t)
        expected_t = gamma * (v.t - sign * beta * v.pz)

        v.boost(boost_vector, -sign)

        assert np.isclose(v.z, expected_pz, atol=1e-6), f"Expected: {expected_pz}, Got: {v.z}"
        assert np.isclose(v.t, expected_t, atol=1e-6), f"Expected: {expected_t}, Got: {v.t}"


# Validate each rotation axis including the unchanged components and energy
@pytest.mark.parametrize('method,components,expected', [
    ('rotateX', (0, 1, 0, 1), (0, 0, 1, 1)),
    ('rotateY', (1, 0, 0, 1), (0, 0, -1, 1)),
    ('rotateZ', (1, 0, 0, 1), (0, 1, 0, 1)),
])
def test_vec4_rotation(method, components, expected):
    vector = vec4(*components)
    getattr(vector, method)(np.pi / 2)
    assert tuple(vector) == pytest.approx(expected, abs=1e-6)


# Preserve finite forward pseudorapidity and its invariance under momentum rescaling
@pytest.mark.parametrize("scale", [1.0e-100, 1.0, 1.0e100])
@pytest.mark.parametrize("direction", [-1, 1])
def test_eta_forward_scale_invariance(scale, direction):
    p = vec4(scale, 0.0, direction * 1.0e9 * scale, 1.0e9 * scale)
    assert p.eta == pytest.approx(direction * np.arcsinh(1.0e9))
    assert vec4(0.0, 0.0, direction, 1.0).eta == direction * np.inf


# Check boost inverses and invariant mass independently of the boost momentum scale
@pytest.mark.parametrize("scale", [1.0e-100, 1.0, 1.0e100])
def test_boost_scale_invariance(scale):
    original = vec4(0.2, -0.3, 0.4, 1.5)
    p = original.copy()
    b = vec4(0.0, 0.0, 3.0 * scale, 5.0 * scale)
    p.boost(b)
    np.testing.assert_allclose(tuple(p), [0.2, -0.3, 1.25 * (0.4 - 0.6 * 1.5), 1.25 * (1.5 - 0.6 * 0.4)])
    assert p.m2 == pytest.approx(original.m2)
    p.boost(b, sign=1)
    np.testing.assert_allclose(tuple(p), tuple(original), atol=1.0e-14)


# Reject boost momenta with no physical massive rest frame
@pytest.mark.parametrize("b", [vec4(), vec4(0., 0., 1., 1.), vec4(0., 0., 2., 1.)])
def test_boost_requires_timelike_momentum(b):
    with pytest.raises(ValueError):
        vec4(0.0, 0.0, 0.0, 1.0).boost(b)


# Compare an LHC boost with the exact mass of the supplied binary floating-point momentum
@pytest.mark.parametrize("scale", [1.0e-100, 1.0, 1.0e100])
@pytest.mark.parametrize("sign", [-1, 1])
def test_boost_lhc_rest_particle(scale, sign):
    b = vec4(1700.0 * scale, -900.0 * scale,
             np.sqrt(6500.0**2 - 1700.0**2 - 900.0**2) * scale,
             np.sqrt(6500.0**2 + 0.9382720813**2) * scale)
    with localcontext() as context:
        context.prec = 80
        x, y, z, energy = (Decimal.from_float(float(component)) for component in b)
        mass = (energy * energy - x * x - y * y - z * z).sqrt()
        expected = [float(sign * component / mass) for component in (x, y, z)]
        expected.append(float(energy / mass))
    particle = vec4(0.0, 0.0, 0.0, 1.0)
    particle.boost(b, sign)
    np.testing.assert_allclose(tuple(particle), expected, rtol=3.0e-10, atol=0.0)


# Leave the original momentum untouched when a boost cannot produce finite components
@pytest.mark.parametrize("energy", [np.inf, np.nan, np.finfo(float).max])
def test_boost_nonfinite_result_rejected_atomic(energy):
    particle = vec4(0.0, 0.0, 0.0, energy)
    with pytest.raises(ValueError):
        particle.boost(vec4(0.0, 0.0, 0.9, 1.0), sign=1)
    np.testing.assert_equal(tuple(particle), (0.0, 0.0, 0.0, energy))


# Preserve shared vectors while isolating boosts and mutable momentum components
def test_deepcopy_momenta():
    original = vec4(0.3, 0.4, 0.5, 2.0)
    repeated = copy.deepcopy([original, original])
    assert repeated[0] is repeated[1] and repeated[0] is not original
    repeated[0].boost(b=vec4(0.2, -0.1, 0.3, 1.5), sign=1)
    np.testing.assert_allclose(repeated[0].m2, original.m2)
    np.testing.assert_allclose(tuple(original), [0.3, 0.4, 0.5, 2.0])
    momenta = vec4(np.array([0.3, 0.4]), 0.0, 0.0, np.array([1.0, 2.0]))
    independent = copy.deepcopy(momenta)
    independent.x[0] += 1.0
    np.testing.assert_allclose(momenta.x, [0.3, 0.4])
