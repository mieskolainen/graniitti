// Covariant XP pole propagator numerator tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

// C++
#include <array>
#include <cmath>
#include <complex>
#include <limits>

// Own
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace {

// Compute one compact particle with the requested pole quantum numbers
gra::MParticle Particle(int pdg, int spin_x2, double mass) {
  gra::MParticle particle;
  particle.pdg = pdg;
  particle.spinX2 = spin_x2;
  particle.mass = mass;
  return particle;
}

// Require one square complex matrix to equal the identity
void RequireIdentity(const gra::MMatrix<std::complex<double>> &matrix,
                     double tolerance = 1.0e-10) {
  REQUIRE(matrix.size_row() == matrix.size_col());
  for (std::size_t row = 0; row < matrix.size_row(); ++row) {
    for (std::size_t column = 0; column < matrix.size_col(); ++column) {
      const std::complex<double> expected = row == column ? 1.0 : 0.0;
      REQUIRE(std::abs(matrix[row][column] - expected) < tolerance);
    }
  }
}

} // namespace

TEST_CASE("XP scalar pole numerator has unit residue",
          "[gra::xpom][propagator]") {
  const auto scalar = Particle(211, 0, 0.13957);
  const auto numerator =
      gra::xpom::Numerator(scalar, gra::M4Vec(0.13, -0.08, 0.21, 0.31));
  CHECK(numerator.type == gra::xpom::PoleType::Scalar);
  CHECK(numerator.helicities == std::vector<double>{0.0});
  RequireIdentity(numerator.matrix);
}

TEST_CASE("XP Dirac and Proca pole numerators have canonical residues",
          "[gra::xpom][propagator]") {
  SECTION("Dirac") {
    const auto proton = Particle(2212, 1, 0.9382720813);
    const gra::M4Vec on_shell(0.17, -0.11, 0.39,
                              std::sqrt(0.17 * 0.17 + 0.11 * 0.11 +
                                        0.39 * 0.39 +
                                        proton.mass * proton.mass));
    const auto pole = gra::xpom::Numerator(proton, on_shell);
    CHECK(pole.type == gra::xpom::PoleType::Dirac);
    CHECK(pole.helicities == std::vector<double>{-0.5, 0.5});
    RequireIdentity(pole.matrix);

    auto off_shell = on_shell;
    off_shell.SetE(on_shell.E() + 0.23);
    const auto continued = gra::xpom::Numerator(proton, off_shell);
    CHECK(continued.matrix.IsFinite());
    CHECK(continued.matrix.FrobNorm2() > 0.0);

    const auto reduced = gra::xpom::ReducedNumerator(proton, off_shell);
    CHECK(reduced.type == gra::xpom::PoleType::Dirac);
    CHECK(reduced.helicities == continued.helicities);
    RequireIdentity(reduced.matrix);
    CHECK((continued.matrix - reduced.matrix).FrobNorm2() > 1.0e-8);
  }

  SECTION("Proca") {
    const auto vector = Particle(113, 2, 0.77526);
    const gra::M4Vec on_shell(-0.14, 0.09, 0.31,
                              std::sqrt(0.14 * 0.14 + 0.09 * 0.09 +
                                        0.31 * 0.31 +
                                        vector.mass * vector.mass));
    const auto pole = gra::xpom::Numerator(vector, on_shell);
    CHECK(pole.type == gra::xpom::PoleType::Proca);
    CHECK(pole.helicities == std::vector<double>{-1.0, 0.0, 1.0});
    RequireIdentity(pole.matrix);
  }
}

TEST_CASE("XP photon pole current obeys the Ward identity",
          "[gra::xpom][propagator]") {
  const auto photon = Particle(gra::PDG::PDG_gamma, 2, 0.0);
  const gra::M4Vec momentum(0.21, -0.17, 0.43,
                            std::sqrt(0.21 * 0.21 + 0.17 * 0.17 + 0.43 * 0.43));
  const auto numerator = gra::xpom::Numerator(photon, momentum);
  CHECK(numerator.type == gra::xpom::PoleType::Photon);
  CHECK(numerator.helicities == std::vector<double>{-1.0, 1.0});
  RequireIdentity(numerator.matrix);

  constexpr std::array<double, 4> metric = {1.0, -1.0, -1.0, -1.0};
  for (const int helicity : {-1, 1}) {
    const auto field = gra::xpom::FieldStrength(momentum, helicity);
    for (std::size_t nu = 0; nu < 4; ++nu) {
      std::complex<double> ward = 0.0;
      for (std::size_t mu = 0; mu < 4; ++mu) {
        ward += metric[mu] * momentum[mu] * field[mu][nu];
      }
      REQUIRE(std::abs(ward) < 1.0e-10);
    }
  }
}

// Distinguish the projected Proca numerator from the unit reduced pole coefficient
TEST_CASE("XP Proca projection is distinct from reduced pole sewing",
          "[gra::xpom][propagator][normalization][covariance]") {
  const auto vector = Particle(113, 2, 1.0);
  const std::array<std::array<double, 3>, 3> directions = {
      std::array<double, 3>{0.0, 0.0, 1.0}, {1.0, 0.0, 0.0}, {0.36, -0.48, 0.8}};
  for (const auto &direction : directions) {
    for (const double p : {0.0, 0.4, 1.0, 2.0}) {
      const double pole_energy = std::sqrt(p * p + vector.mass * vector.mass);
      for (const double energy : {0.0, -0.3, 0.6, pole_energy}) {
        CAPTURE(direction, p, energy);
        const gra::M4Vec momentum(p * direction[0], p * direction[1], p * direction[2], energy);
        const auto numerator = gra::xpom::Numerator(vector, momentum);
        // epsilon_0.k = p (E - E_pole) / m in the physical pole basis
        const double overlap = p * (energy - pole_energy) / (vector.mass * vector.mass);
        gra::MMatrix<std::complex<double>> expected(3, 3, "eye");
        expected[1][1] += overlap * overlap;
        CHECK((numerator.matrix - expected).FrobNorm2() < 1.0e-20);
        const auto reduced = gra::xpom::ReducedNumerator(vector, momentum);
        RequireIdentity(reduced.matrix);
      }
    }
  }
}

// Count non-finite pole projections as sampling failures even for finite input components
TEST_CASE("XP pole numerator numerical failures use amplitude bookkeeping",
          "[gra::xpom][propagator][failure]") {
  const auto vector = Particle(113, 2, 1.0);
  const auto proton = Particle(2212, 1, 1.0);
  const double infinity = std::numeric_limits<double>::infinity();
  const double nan = std::numeric_limits<double>::quiet_NaN();
  for (const auto &particle : {vector, proton}) {
    for (const double energy : {infinity, nan}) {
      REQUIRE_THROWS_AS(gra::xpom::Numerator(particle, gra::M4Vec(0.0, 0.0, 1.0, energy)),
                        gra::AmplitudeFailure);
      REQUIRE_THROWS_AS(gra::xpom::ReducedNumerator(particle, gra::M4Vec(0.0, 0.0, 1.0, energy)),
                        gra::AmplitudeFailure);
    }
  }
  // The Proca k^mu k^nu term overflows while the physical pole basis remains finite
  REQUIRE_FALSE(gra::xpom::Numerator(vector, gra::M4Vec(0.0, 0.0, 1.0, 1.0e200)).matrix.IsFinite());
  RequireIdentity(gra::xpom::ReducedNumerator(vector, gra::M4Vec(0.0, 0.0, 1.0, 1.0e200)).matrix);
}

// Compare reduced expressions to covariant spinor and polarization contractions
TEST_CASE("XP pole projections retain the absolute covariant residues",
          "[gra::xpom][propagator][normalization][physics]") {
  const gra::MDirac dirac("DIRAC");
  constexpr std::array<double, 4> metric = {1.0, -1.0, -1.0, -1.0};
  for (const double mass : {0.5, 1.0, 2.0}) {
    const gra::M4Vec pole(0.3, -0.4, 0.7, std::sqrt(0.74 + mass * mass));
    for (const double energy : {-0.3, 0.0, 0.6, pole.E()}) {
      auto momentum = pole;
      momentum.SetE(energy);
      const auto proton = Particle(2212, 1, mass);
      const auto vector = Particle(113, 2, mass);
      const auto fermion = gra::xpom::Numerator(proton, momentum).matrix;
      const auto proca = gra::xpom::Numerator(vector, momentum).matrix;
      const auto slash = dirac.FSlash(momentum) + gra::MMatrix<std::complex<double>>::IdentityMatrix(4) * mass;
      const std::array<int, 2> spin = {-1, 1};
      for (const auto i : indices(spin)) {
        const auto left = dirac.Bar(dirac.uHelDirac(pole, spin[i]));
        for (const auto j : indices(spin)) {
          const auto right = dirac.uHelDirac(pole, spin[j]);
          const auto expected = slash.BilinearForm(left, right) / (4.0 * mass * mass);
          CHECK(std::abs(fermion[i][j] - expected) < 1.0e-11);
        }
      }
      gra::MMatrix<std::complex<double>> numerator(4, 4, 0.0);
      for (const auto mu : indices(metric)) {
        for (const auto nu : indices(metric)) {
          numerator[mu][nu] = (mu == nu ? -metric[mu] : 0.0) + momentum[mu] * momentum[nu] / (mass * mass);
        }
      }
      const std::array<int, 3> helicities = {-1, 0, 1};
      for (const auto i : indices(helicities)) {
        const auto first = dirac.EpsMassiveSpin1(pole, helicities[i]);
        for (const auto j : indices(helicities)) {
          const auto second = dirac.EpsMassiveSpin1(pole, helicities[j]);
          std::array<std::complex<double>, 4> left{};
          std::array<std::complex<double>, 4> right{};
          for (const auto mu : indices(metric)) {
            left[mu] = metric[mu] * std::conj(first(mu));
            right[mu] = metric[mu] * second(mu);
          }
          CHECK(std::abs(proca[i][j] - numerator.BilinearForm(left, right)) < 1.0e-11);
        }
      }
    }
  }
}

// Keep unit pole coefficients when explicit boosted projections lose the mass scale
TEST_CASE("XP massive pole residues remain stable at large boosts",
          "[gra::xpom][propagator][normalization][covariance]") {
  for (const double mass : {0.5, 1.0, 2.0}) {
    for (const double boost : {1.0e3, 1.0e4, 1.0e6, 1.0e9}) {
      for (const gra::M4Vec direction : {gra::M4Vec(0.0, 0.0, 1.0, 0.0), gra::M4Vec(0.36, -0.48, 0.8, 0.0)}) {
        auto momentum = direction * (mass * boost);
        const double p = std::hypot(momentum.Px(), momentum.Py(), momentum.Pz());
        momentum.SetE(std::hypot(p, mass));
        CAPTURE(mass, boost);
        RequireIdentity(gra::xpom::Numerator(Particle(2212, 1, mass), momentum).matrix);
        RequireIdentity(gra::xpom::Numerator(Particle(113, 2, mass), momentum).matrix);
      }
    }
  }
}

// Reject invalid massive pole input consistently in full and reduced projections
TEST_CASE("XP massive pole projections validate the nominal mass",
          "[gra::xpom][propagator][params]") {
  for (const double mass : {0.0, -1.0, std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    for (const int spin : {1, 2}) {
      const auto particle = Particle(spin == 1 ? 2212 : 113, spin, mass);
      const gra::M4Vec momentum(0.0, 0.0, 1.0, 2.0);
      REQUIRE_THROWS_AS(gra::xpom::Numerator(particle, momentum), std::invalid_argument);
      REQUIRE_THROWS_AS(gra::xpom::ReducedNumerator(particle, momentum), std::invalid_argument);
    }
  }
}

// Reject non-finite photon tensors before they enter coherent amplitude sums
TEST_CASE("XP photon field strength rejects invalid generated momenta",
          "[gra::xpom][propagator][photon][failure]") {
  for (const double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    REQUIRE_THROWS_AS(gra::xpom::FieldStrength(gra::M4Vec(0.0, 0.0, 1.0, invalid), 1), gra::AmplitudeFailure);
    REQUIRE_THROWS_AS(gra::xpom::FieldStrength(gra::M4Vec(invalid, 0.0, 1.0, 2.0), -1), gra::AmplitudeFailure);
  }
}

// Sew every physical spin once for either orientation of the exchanged momentum
TEST_CASE("XP reduced spin sums preserve crossing and absolute normalization",
          "[gra::xpom][propagator][normalization][crossing]") {
  for (const auto &particle : {Particle(211, 0, 0.5), Particle(2212, 1, 1.0), Particle(-2212, 1, 1.0),
                               Particle(113, 2, 0.5), Particle(333, 2, 2.0), Particle(22, 2, 0.0)}) {
    const std::size_t states = particle.pdg == 22 ? 2 : particle.spinX2 + 1;
    for (const double energy : {0.0, 0.3, std::hypot(1.0, particle.mass)}) {
      const gra::M4Vec momentum(0.36, -0.48, 0.8, energy);
      for (const double orientation : {-1.0, 1.0}) {
        CAPTURE(particle.pdg, energy, orientation);
        const auto numerator = gra::xpom::ReducedNumerator(particle, momentum * orientation);
        const auto metric = gra::xpom::PoleMetric(particle, momentum * orientation);
        REQUIRE(numerator.helicities.size() == states);
        REQUIRE(metric.projections == numerator.helicities);
        RequireIdentity(numerator.matrix);
        RequireIdentity(metric.matrix);
        // A unit pole coefficient sums spin states without averaging by 2s+1 or inserting 2m
        CHECK(numerator.matrix.FrobNorm2() == Approx(static_cast<double>(states)).epsilon(1.0e-12));
      }
    }
  }
}
