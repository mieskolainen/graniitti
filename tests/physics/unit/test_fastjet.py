# Unit tests for Numba anti-kT jet clustering
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

import numpy as np
import pytest
from core.kinematics import fastjet


# Compute one massless Cartesian four-vector
def massless_four_vector(pt, eta, phi):
    return pt * math.cos(phi), pt * math.sin(phi), pt * math.sinh(eta), pt * math.cosh(eta)


# Convert particle tuples into component arrays
def component_arrays(particles):
    values = np.asarray(particles, dtype=np.float64)
    return values[:, 0], values[:, 1], values[:, 2], values[:, 3]


# Check that nearby particles are merged with E-scheme recombination
def test_antikt_merges_nearby_particles():
    particles = [
        massless_four_vector(30.0, 0.0, 0.0),
        massless_four_vector(10.0, 0.1, 0.1),
        massless_four_vector(25.0, 0.0, math.pi),
    ]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.6, 5.0, 5.0)

    assert jets.shape == (2, 4)
    assert labels[0] == labels[1]
    assert labels[2] != labels[0]
    assert np.isclose(jets[0, 0], particles[0][0] + particles[1][0])
    assert np.isclose(jets[0, 1], particles[0][1] + particles[1][1])


# Check anti-kT assigns a soft particle by the hard jet distance rather than angle alone
# [REFERENCE: https://arxiv.org/abs/0802.1189]
def test_antikt_soft_particle_prefers_harder_jet():
    particles = np.asarray([massless_four_vector(100.0, 0.0, 0.0),
                            massless_four_vector(10.0, 0.0, 1.0),
                            massless_four_vector(1.0, 0.0, 0.7)])
    # The soft particle is nearer to the 10 GeV jet, but its anti-kT distance favors the 100 GeV jet
    assert 0.7**2 / 100.0**2 < 0.3**2 / 10.0**2
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.8)
    np.testing.assert_array_equal(labels, [0, 1, 0])
    np.testing.assert_allclose(jets, [particles[0] + particles[2], particles[1]], rtol=1e-12, atol=1e-12)


# Check that azimuthal wrapping merges particles across the pi boundary
def test_antikt_wraps_azimuth():
    particles = [massless_four_vector(20.0, 0.0, math.pi - 0.05), massless_four_vector(10.0, 0.0, -math.pi + 0.05)]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.4, 0.0, 5.0)

    assert jets.shape == (1, 4)
    assert labels.tolist() == [0, 0]


# Check nearest-neighbour lookup across a rapidity tile boundary
def test_antikt_merges_across_rapidity_tile_boundary():
    particles = [massless_four_vector(20.0, 0.09, 0.0), massless_four_vector(10.0, 0.11, 0.0)]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.05, 0.0, 5.0)

    assert jets.shape == (1, 4)
    assert labels.tolist() == [0, 0]


# Check transverse momentum and pseudorapidity jet selection
def test_antikt_applies_fiducial_jet_selection():
    particles = [
        massless_four_vector(30.0, 0.0, 0.0),
        massless_four_vector(15.0, 0.0, 2.0),
        massless_four_vector(40.0, 3.0, -2.0),
    ]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.4, 20.0, 2.5)

    assert jets.shape == (1, 4)
    assert labels.tolist() == [0, -1, -1]


# Check invalid clustering parameters return an empty result
@pytest.mark.parametrize(
    "radius,pt_min,abs_eta_max",
    [
        (0.0, 20.0, 2.5),
        (-0.4, 20.0, 2.5),
        (math.nan, 20.0, 2.5),
        (math.inf, 20.0, 2.5),
        (1.0e-200, 20.0, 2.5),
        (1.0e200, 20.0, 2.5),
        (0.4, math.nan, 2.5),
        (0.4, -1.0, 2.5),
        (0.4, 20.0, math.nan),
        (0.4, 20.0, 0.0),
    ],
)
def test_antikt_invalid_params_empty_result(radius, pt_min, abs_eta_max):
    particles = [massless_four_vector(30.0, 0.0, 0.0)]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), radius, pt_min, abs_eta_max)

    assert jets.shape == (0, 4)
    assert labels.size == 0


# Check constituent conservation and covariance under rotations, boosts and particle order
@pytest.mark.parametrize("massive", [False, True])
@pytest.mark.parametrize("radius", [0.05, 0.4, 0.6, 1.2, 2.2, 4.0])
def test_antikt_conserves_momentum_symmetries(massive, radius):
    for seed in range(8):
        random = np.random.default_rng(seed)
        particles = np.asarray([massless_four_vector(
            random.uniform(0.5, 50.0), random.uniform(-3.0, 3.0), random.uniform(-math.pi, math.pi))
            for _ in range(24)])
        if massive:
            particles[:, 3] = np.hypot(particles[:, 3], random.uniform(0.0, 20.0, len(particles)))
        jets, labels = fastjet.cluster_antikt(*component_arrays(particles), radius, 0.0, 10.0)
        assert np.all(labels >= 0)
        assert len(set(labels)) == len(jets)
        for index, jet in enumerate(jets):
            np.testing.assert_allclose(particles[labels == index].sum(axis=0), jet, rtol=1e-12, atol=1e-12)
        angle, rapidity = 0.73, 0.8
        rotation = np.array([[math.cos(angle), -math.sin(angle)], [math.sin(angle), math.cos(angle)]])
        boost = np.array([[math.cosh(rapidity), math.sinh(rapidity)], [math.sinh(rapidity), math.cosh(rapidity)]])
        transformed = particles.copy()
        expected = jets.copy()
        for momenta in (transformed, expected):
            momenta[:, :2] = momenta[:, :2] @ rotation.T
            momenta[:, 2:] = momenta[:, 2:] @ boost.T
        permutation = random.permutation(len(particles))
        result, assigned = fastjet.cluster_antikt(*component_arrays(transformed[permutation]), radius, 0.0, 10.0)
        np.testing.assert_allclose(result, expected, rtol=1e-10, atol=1e-10)
        np.testing.assert_array_equal(assigned[np.argsort(permutation)], labels)


# Check massive clustering uses rapidity and E-scheme four-momentum recombination
# [REFERENCE: https://arxiv.org/abs/0802.1189]
def test_antikt_massive_rapidity_distance():
    particles = np.asarray([(10.0, 0.0, 30.0, math.sqrt(10.0**2 + 30.0**2 + 40.0**2)),
                            massless_four_vector(5.0, 0.75, 0.0)])
    assert abs(math.asinh(3.0) - 0.75) > 0.4
    assert abs(math.atanh(particles[0, 2] / particles[0, 3]) - 0.75) < 0.4
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.4)
    np.testing.assert_array_equal(labels, [0, 0])
    np.testing.assert_allclose(jets, [particles.sum(axis=0)], rtol=1e-12)


# Compare massive E-scheme momenta and constituents with the independent FastJet library
@pytest.mark.parametrize("radius", [0.05, 0.4, 1.2, 4.0])
def test_antikt_matches_fastjet(radius):
    fj = pytest.importorskip("fastjet")
    random = np.random.default_rng(614)
    for _ in range(12):
        particles = np.asarray([massless_four_vector(random.uniform(0.1, 100.0),
                               random.uniform(-4.0, 4.0), random.uniform(-math.pi, math.pi))
                               for _ in range(40)])
        particles[:, 3] = np.hypot(particles[:, 3], random.uniform(0.0, 20.0, len(particles)))
        inputs = [fj.PseudoJet(*p) for p in particles]
        for index, particle in enumerate(inputs):
            particle.set_user_index(index)
        sequence = fj.ClusterSequence(inputs, fj.JetDefinition(fj.antikt_algorithm, radius))
        expected = fj.sorted_by_pt(sequence.inclusive_jets())
        jets, labels = fastjet.cluster_antikt(*component_arrays(particles), radius)
        np.testing.assert_allclose(jets, [[p.px(), p.py(), p.pz(), p.e()] for p in expected], rtol=1e-11, atol=1e-11)
        for index, jet in enumerate(expected):
            assert set(np.flatnonzero(labels == index)) == {p.user_index() for p in jet.constituents()}


# Check a collinear split and a soft emission preserve the resolved hard jets
@pytest.mark.parametrize("split", [0.1, 0.5, 0.9])
def test_antikt_soft_collinear_safety(split):
    particles = np.asarray([massless_four_vector(30.0, 0.2, 0.1), massless_four_vector(20.0, -0.4, 2.7)])
    jets, _ = fastjet.cluster_antikt(*component_arrays(particles), 0.4, 1.0)
    collinear = np.vstack((split * particles[0], (1.0 - split) * particles[0], particles[1]))
    soft = np.asarray(massless_four_vector(1e-10, 0.25, 0.15))
    result, labels = fastjet.cluster_antikt(*component_arrays(np.vstack((collinear, soft))), 0.4, 1.0)
    np.testing.assert_allclose(result, jets, rtol=1e-10, atol=1e-10)
    assert labels[0] == labels[1] == labels[3]
    assert labels[2] != labels[0]


# Check component arrays cannot silently discard particles or mix dimensions
@pytest.mark.parametrize("component", range(4))
@pytest.mark.parametrize("shape", [(1,), (2, 1)])
def test_antikt_invalid_component_arrays(component, shape):
    arrays = [np.ones(2) for _ in range(4)]
    arrays[component] = np.ones(shape)
    with pytest.raises(ValueError, match="component"):
        fastjet.cluster_antikt(*arrays, 0.4)


# Check soft constituents and rapidities are invariant under a common momentum scale
@pytest.mark.parametrize("scale", [1.0e-14, 1.0e-100])
def test_antikt_preserves_soft_constituents(scale):
    particles = np.asarray(
        [
            massless_four_vector(30.0, 1.0, 0.0),
            massless_four_vector(10.0, 1.1, 0.1),
            massless_four_vector(25.0, -1.0, 2.0),
        ]
    )
    particles[:, 3] = np.hypot(particles[:, 3], 5.0)
    expected, expected_labels = fastjet.cluster_antikt(*component_arrays(particles), 0.4)
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles * scale), 0.4)

    np.testing.assert_array_equal(labels, expected_labels)
    np.testing.assert_allclose(jets / scale, expected, rtol=1.0e-12, atol=1.0e-12)
    for particle in particles:
        expected_rapidity = math.atanh(particle[2] / particle[3])
        assert fastjet.jet_rapidity(*(particle * scale)) == pytest.approx(expected_rapidity)


# Check E-scheme recombination retains the energy of a massive particle at rest
def test_antikt_merges_massive_particle_rest():
    particles = [(30.0, 0.0, 0.0, 30.0), (0.0, 0.0, 0.0, 2.0)]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 0.4)

    np.testing.assert_array_equal(labels, [0, 0])
    np.testing.assert_allclose(jets, [[30.0, 0.0, 0.0, 32.0]])


# Check transverse cancellation cannot turn a forward jet into a central jet
def test_antikt_rejects_jet_along_beam():
    particles = [(10.0, 0.0, 10.0, math.sqrt(200.0)), (-10.0, 0.0, 10.0, math.sqrt(200.0))]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 4.0, 0.0, 2.5)

    assert jets.shape == (0, 4)
    np.testing.assert_array_equal(labels, [-1, -1])
    assert math.isinf(fastjet.jet_eta(0.0, 0.0, 20.0))
    assert fastjet.jet_eta(0.0, 0.0, -20.0) < 0.0


# Check large angular scales do not overflow distances for a vanishing jet transverse momentum
def test_antikt_massive_rest_jet_large_radius():
    particles = [(0.0, 0.0, 0.0, 2.0)]
    jets, labels = fastjet.cluster_antikt(*component_arrays(particles), 1.0e6)

    np.testing.assert_allclose(jets, particles)
    np.testing.assert_array_equal(labels, [0])
