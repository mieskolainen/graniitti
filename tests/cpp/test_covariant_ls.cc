// Covariant XP central LS tensor tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

// C++
#include <cmath>
#include <complex>
#include <limits>

// Own
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace {

// Check one complex value component by component
void RequireComplexNear(const std::complex<double>& actual, const std::complex<double>& expected,
                        const double tolerance) {
  CHECK(std::real(actual) == Approx(std::real(expected)).margin(tolerance));
  CHECK(std::imag(actual) == Approx(std::imag(expected)).margin(tolerance));
}

// Compute a compact particle carrying the requested integer spin and parity
gra::MParticle Particle(int pdg, int spin, int parity = 1, int c_parity = 1) {
  gra::MParticle particle;
  particle.pdg    = pdg;
  particle.spinX2 = 2 * spin;
  particle.P      = parity;
  particle.C      = c_parity;
  return particle;
}

}  // namespace

// Check every tensor-tensor pole operator against the Jacob-Wick LS equation
TEST_CASE("Fixed-spin resonance poles obey Jacob-Wick algebra for every parity",
          "[gra::spin][covariant-ls][resonance][parity][regression]") {
  const auto       tensor        = Particle(995, 2, 1, 1);
  constexpr double momentum      = 0.43;
  constexpr double exchange_mass = 0.91;
  const double     energy        = std::sqrt(exchange_mass * exchange_mass + momentum * momentum);
  gra::LORENTZSCALAR lts;
  auto& q1 = lts.q1_in_X;
  q1 = gra::M4Vec(0.0, 0.0, momentum, energy);
  std::size_t      tested = 0;

  for (int spin = 0; spin <= 6; ++spin) {
    for (const int parity : {-1, 1}) {
      const auto mother = Particle(910000 + 10 * spin + parity, spin, parity, 1);
      const auto operators =
          gra::spin::CanonicalPoleOperators(mother, tensor, tensor, true, true, true, gra::spin::VertexContext::Auto);
      REQUIRE_FALSE(operators.empty());

      for (const auto& op : operators) {
        CAPTURE(spin, parity, op.coupling.l, op.coupling.two_s);
        const std::complex<double> coupling =
            std::polar(0.67 + 0.01 * static_cast<double>(tested), -0.31 + 0.02 * static_cast<double>(tested));
        const auto vertex =
            gra::spin::PreparePoleLS(mother, tensor, tensor, {{op.coupling.l, op.coupling.two_s, coupling}}, 1.0, true,
                                     true, true, gra::spin::VertexContext::Auto, 0.0, false);
        const auto   actual     = gra::xpom::Fusion(lts, vertex);
        const double total_spin = 0.5 * static_cast<double>(op.coupling.two_s);
        const double jw =
            std::sqrt((2.0 * static_cast<double>(op.coupling.l) + 1.0) / (2.0 * static_cast<double>(spin) + 1.0));
        const double raw =
            gra::spin::RawPoleLSNormalization(4, 4, 2 * spin, op.coupling.l, static_cast<int>(op.coupling.two_s));
        CHECK(raw == Approx(op.raw_normalization).epsilon(2.0e-13));

        for (std::size_t row = 0; row < vertex.helicity.lambda_values.size_row(); ++row) {
          const double               m1      = vertex.helicity.lambda_values[row][0];
          const double               m2      = vertex.helicity.lambda_values[row][1];
          const double               lambda  = m1 - m2;
          const std::complex<double> reduced = coupling * raw * jw *
                                               gra::wigner::CG(static_cast<double>(op.coupling.l), total_spin, 0.0,
                                                               lambda, static_cast<double>(spin), lambda) *
                                               gra::wigner::CG(2.0, 2.0, m1, -m2, total_spin, lambda);

          for (const auto& column : indices(vertex.helicity.Jz_values)) {
            const double               projection = vertex.helicity.Jz_values[column];
            const std::complex<double> expected =
                std::abs(lambda + projection) < 1.0e-12
                    ? gra::spin::JacobWickSecondLegReversalPhase(static_cast<double>(spin), projection) * reduced
                    : 0.0;
            RequireComplexNear(actual[row][column], expected, 4.0e-12);
          }
        }
        ++tested;
      }
    }
  }
  CHECK(tested > 50);
}

TEST_CASE("Covariant LS event failures reject invalid generated momenta", "[gra::spin][covariant-ls][failure]") {
  const auto   spin_two = Particle(900010, 2);
  const auto   first    = Particle(800006, 0);
  const auto   second   = Particle(800007, 0);
  const auto   vertex   = gra::spin::PreparePoleLS(spin_two, first, second, {{2, 0, 1.0}}, 1.0);
  const double nan      = std::numeric_limits<double>::quiet_NaN();
  const double maximum  = std::numeric_limits<double>::max();

  REQUIRE_THROWS_AS(gra::spin::PoleLSReduced(vertex, nan), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(gra::spin::PoleLSReduced(vertex, -0.1), gra::AmplitudeFailure);
  REQUIRE_FALSE(gra::spin::PoleLSReduced(vertex, maximum).IsFinite());

  gra::LORENTZSCALAR lts;
  lts.q1_in_X = gra::M4Vec(nan, 0.0, 0.4, 0.8);
  REQUIRE_THROWS_AS(gra::xpom::Fusion(lts, vertex), gra::AmplitudeFailure);
  lts.q1_in_X = gra::M4Vec(maximum, maximum, maximum, 0.8);
  REQUIRE_THROWS_AS(gra::xpom::Fusion(lts, vertex), gra::AmplitudeFailure);
}

TEST_CASE("Photon STF vertices use the transverse pole projection", "[gra::spin][covariant-ls]") {
  const auto mother = Particle(900006, 1, -1, -1);
  const auto photon = Particle(22, 1, -1, -1);
  const auto scalar = Particle(991, 0, 1, 1);
  const auto vertex = gra::spin::PreparePoleLS(mother, photon, scalar, {{0, 2, 1.0}}, 1.0);
  REQUIRE(vertex.leg1_transverse_pole);
  REQUIRE_FALSE(vertex.leg2_transverse_pole);

  const auto reduced = gra::spin::PoleLSReduced(vertex, 0.4);
  REQUIRE(std::abs(reduced[0][0]) == Approx(1.0));
  REQUIRE(std::abs(reduced[1][0]) == Approx(0.0).margin(1e-13));
  REQUIRE(std::abs(reduced[2][0]) == Approx(1.0));

  gra::LORENTZSCALAR lts;
  auto& q1 = lts.q1_in_X;
  q1 = gra::M4Vec(0.0, 0.0, 0.4, 0.7);
  const auto       central = gra::xpom::Fusion(lts, vertex);
  REQUIRE(central.size_row() == 2);
  REQUIRE(central.size_col() == 3);
  REQUIRE(central.FrobNorm2() == Approx(2.0));
}

TEST_CASE("Covariant LS central rows follow the common azimuth section", "[gra::spin][covariant-ls]") {
  const auto   mother = Particle(900007, 1);
  const auto   vector = Particle(993, 1);
  const auto   vertex = gra::spin::PreparePoleLS(mother, vector, vector, {{2, 4, {0.7, -0.2}}}, 1.0, true, false, true);
  gra::LORENTZSCALAR lts;
  auto& q1 = lts.q1_in_X;
  q1 = gra::M4Vec(0.31, -0.17, 0.42, 0.83);
  const auto   reference = gra::xpom::Fusion(lts, vertex);
  const double angle     = 0.63;
  q1.RotateZ(angle);
  const auto rotated = gra::xpom::Fusion(lts, vertex);

  for (std::size_t row = 0; row < rotated.size_row(); ++row) {
    const double lambda = vertex.helicity.lambda_values[row][0] - vertex.helicity.lambda_values[row][1];
    for (std::size_t col = 0; col < rotated.size_col(); ++col) {
      const auto phase    = std::exp(std::complex<double>(0.0, (-vertex.helicity.Jz_values[col] - lambda) * angle));
      const auto expected = reference[row][col] * phase;
      REQUIRE(rotated[row][col].real() == Approx(expected.real()).margin(1e-12));
      REQUIRE(rotated[row][col].imag() == Approx(expected.imag()).margin(1e-12));
    }
  }
}

TEST_CASE("Covariant LS production uses the Cartesian STF spin metric", "[gra::spin][covariant-ls]") {
  const auto       mother = Particle(900008, 1, -1);
  const auto       first  = Particle(800003, 0);
  const auto       second = Particle(800004, 0);
  const auto       vertex = gra::spin::PreparePoleLS(mother, first, second, {{1, 0, 1.0}}, 1.0, true, false, true);
  gra::LORENTZSCALAR lts;
  auto& q1 = lts.q1_in_X;
  q1 = gra::M4Vec(0.0, 0.0, 0.4, 0.7);
  const auto       reduced = gra::spin::PoleLSReduced(vertex, 0.4);
  const auto       central = gra::xpom::Fusion(lts, vertex);
  REQUIRE(central.size_row() == 1);
  REQUIRE(central.size_col() == 3);
  REQUIRE(std::abs(central[0][0]) == Approx(0.0).margin(1e-13));
  REQUIRE(central[0][1].real() == Approx(-reduced[0][0].real()));
  REQUIRE(central[0][1].imag() == Approx(-reduced[0][0].imag()));
  REQUIRE(std::abs(central[0][2]) == Approx(0.0).margin(1e-13));
}

TEST_CASE("Covariant LS central matrix follows the analytic STF rotation", "[gra::spin][covariant-ls]") {
  const auto mother = Particle(900009, 2, 1);
  const auto vector = Particle(800005, 1, -1);
  const auto vertex = gra::spin::PreparePoleLS(mother, vector, vector, {{2, 4, {0.6, -0.25}}}, 0.9, true, false, true);
  gra::LORENTZSCALAR lts;
  auto& q1 = lts.q1_in_X;
  q1 = gra::M4Vec(0.29, -0.21, 0.37, 0.83);
  const gra::M4Vec q2(-0.29, 0.21, -0.37, 0.91);
  const auto       relative = (q1 - q2) * 0.5;
  const auto       reduced  = gra::spin::PoleLSReduced(vertex, relative.P3mod());
  const auto       central  = gra::xpom::Fusion(lts, vertex);

  for (std::size_t row = 0; row < central.size_row(); ++row) {
    const double      lambda = vertex.helicity.lambda_values[row][0] - vertex.helicity.lambda_values[row][1];
    const std::size_t i1     = vertex.helicity.lambda_idx[row][0];
    const std::size_t i2     = vertex.helicity.lambda_idx[row][1];
    for (std::size_t col = 0; col < central.size_col(); ++col) {
      const double M        = vertex.helicity.Jz_values[col];
      const auto   expected = gra::spin::JacobWickSecondLegReversalPhase(vertex.helicity.J, M) *
                              gra::wigner::d(relative.Theta(), lambda, -M, vertex.helicity.J) *
                              std::exp(std::complex<double>(0.0, -(M + lambda) * relative.Phi())) * reduced[i1][i2];
      REQUIRE(central[row][col].real() == Approx(expected.real()).margin(1e-12));
      REQUIRE(central[row][col].imag() == Approx(expected.imag()).margin(1e-12));
    }
  }
}
