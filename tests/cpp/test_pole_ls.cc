// Shared physical-pole LS algebra and continuum vertex tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <future>
#include <limits>
#include <type_traits>
#include <vector>

// Own
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Regge/MReggeMulti.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MPoleLS.h"
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

// Check two complex matrices component by component
void RequireMatrixNear(const gra::MMatrix<std::complex<double>>& actual,
                       const gra::MMatrix<std::complex<double>>& expected, const double tolerance) {
  REQUIRE(actual.size_row() == expected.size_row());
  REQUIRE(actual.size_col() == expected.size_col());
  for (std::size_t row = 0; row < actual.size_row(); ++row) {
    for (std::size_t column = 0; column < actual.size_col(); ++column) {
      RequireComplexNear(actual[row][column], expected[row][column], tolerance);
    }
  }
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

// Compute a compact particle carrying one doubled spin label
gra::MParticle ParticleX2(int pdg, int spin_x2, int parity = 1, int c_parity = 1) {
  gra::MParticle particle;
  particle.pdg    = pdg;
  particle.spinX2 = spin_x2;
  particle.P      = parity;
  particle.C      = c_parity;
  return particle;
}

}  // namespace

// Check immutable continuum vertices against fresh LS evaluation and independent worker reads
TEST_CASE("Prepared continuum residues own their couplings and share only const data",
          "[gra::spin][pole-ls][cache][threading]") {
  const auto exchange = Particle(800400, 2);
  const auto scalar   = Particle(800401, 0);
  auto       vertex = gra::spin::PreparePoleLS(exchange, scalar, scalar, {{2, 0, {0.7, -0.4}}}, 1.0, true, false, true,
                                               gra::spin::VertexContext::SubTUChannelExchange);
  const gra::spin::PoleResidue prepared(vertex);
  const auto                   copy = prepared;
  static_assert(std::is_same_v<decltype(prepared.Pole()), const gra::spin::PoleLS&>);
  static_assert(std::is_same_v<decltype(prepared.Reduced()), const gra::spin::EvaluatedPoleSubvertex&>);
  REQUIRE(prepared.Ready());
  REQUIRE(&prepared.Pole() == &copy.Pole());
  const auto   reference = gra::spin::EvaluateCrossedPole(vertex);
  const double density   = gra::spin::LeadingPoleDensity(vertex, vertex.Lambda);
  RequireMatrixNear(prepared.Reduced().frame, reference.frame, 1.0e-14);
  REQUIRE(prepared.Density() == Approx(density).epsilon(1.0e-13));

  const std::complex<double> scale{0.3, 0.8};
  vertex.terms[0].coefficient *= scale;
  const gra::spin::PoleResidue changed(vertex);
  RequireMatrixNear(changed.Reduced().frame, reference.frame * scale, 1.0e-13);
  REQUIRE(changed.Density() == Approx(density * std::norm(scale)).epsilon(1.0e-13));
  RequireMatrixNear(prepared.Reduced().frame, reference.frame, 1.0e-14);

  std::vector<std::future<gra::MMatrix<std::complex<double>>>> workers;
  for (int i = 0; i < 4; ++i) {
    workers.push_back(std::async(std::launch::async, [copy] { return copy.Reduced().frame * copy.Density(); }));
  }
  for (auto& worker : workers) { RequireMatrixNear(worker.get(), reference.frame * density, 1.0e-14); }
  vertex.derivative_factor = true;
  REQUIRE_THROWS_AS(gra::spin::PoleResidue(vertex), std::invalid_argument);
}

// Compare projection before sewing with the full tensor for boson and fermion spin bases
TEST_CASE("Projected subvertices preserve the full complex helicity contraction", "[gra::spin][pole-ls][contraction]") {
  const auto exchange = Particle(800410, 2);
  for (const int spin_x2 : {0, 1, 2}) {
    const auto first     = ParticleX2(800411, spin_x2);
    const auto second    = ParticleX2(-800411, spin_x2);
    const auto operators = gra::spin::CanonicalPoleOperators(exchange, first, second, true, false, false);
    std::vector<gra::spin::LSTerm> terms;
    for (const auto& i : indices(operators)) {
      terms.push_back(
          {operators[i].coupling.l, operators[i].coupling.two_s, std::polar(0.3, 0.17 * static_cast<double>(i))});
    }
    const auto vertex = gra::spin::PreparePoleLS(exchange, first, second, terms, 1.0, true, false, false);
    const auto up     = gra::spin::EvaluatePoleSubvertex(vertex, gra::M4Vec(0.2, -0.3, 0.8, 1.1),
                                                         gra::M4Vec(0.1, 0.4, 1.2, 0.3), false);
    const auto dn     = gra::spin::EvaluatePoleSubvertex(vertex, gra::M4Vec(-0.2, 0.3, -0.8, 1.1),
                                                         gra::M4Vec(-0.1, -0.4, 1.2, 0.3), true);
    for (const std::size_t rows : {2U, 4U}) {
      gra::MMatrix<std::complex<double>> upper(rows, up.frame.size_col(), 0.0);
      gra::MMatrix<std::complex<double>> lower(rows, dn.frame.size_col(), 0.0);
      for (std::size_t i = 0; i < rows; ++i) {
        for (std::size_t j = 0; j < upper.size_col(); ++j) {
          if (j % 2 == 0) { upper[i][j] = std::polar(0.4 + 0.1 * i, 0.37 * (j + i)); }
          lower[i][j] = std::polar(0.7 + 0.2 * i, -0.29 * (j + i));
        }
      }
      for (const bool swap : {false, true}) {
        CAPTURE(spin_x2, rows, swap);
        const auto central   = gra::spin::Subchannel(up, dn, first, second, swap);
        const auto full      = gra::spin::Contract(upper, lower, central);
        const auto projected = gra::spin::ProjectedSubchannel(up, dn, upper, lower, first, second, swap);
        RequireMatrixNear(projected, full, 3.0e-11);
      }
    }
  }
}

TEST_CASE("Raw STF scalar contractions keep every spin component", "[gra::spin][pole-ls]") {
  const auto scalar = Particle(900000, 0);
  for (std::size_t j = 0; j <= 4; ++j) {
    CAPTURE(j);
    const auto              exchange = Particle(900100 + static_cast<int>(j), static_cast<int>(j));
    const gra::spin::LSTerm term{0, 0, 1.0};
    const auto              vertex = gra::spin::PreparePoleLS(scalar, exchange, exchange, {term}, 1.0);
    REQUIRE(vertex.raw_normalization[0] == Approx(std::sqrt(2.0 * static_cast<double>(j) + 1.0)));

    const auto reduced = gra::spin::PoleLSReduced(vertex, 0.37);
    REQUIRE(reduced.size_row() == 2 * j + 1);
    REQUIRE(reduced.size_col() == 2 * j + 1);
    for (std::size_t m1 = 0; m1 < reduced.size_row(); ++m1) {
      for (std::size_t m2 = 0; m2 < reduced.size_col(); ++m2) {
        if (m1 == m2) {
          const int    m        = static_cast<int>(m1) - static_cast<int>(j);
          const double expected = ((static_cast<int>(j) - m) % 2 == 0) ? 1.0 : -1.0;
          REQUIRE(std::real(reduced[m1][m2]) == Approx(expected));
          REQUIRE(std::imag(reduced[m1][m2]) == Approx(0.0).margin(1e-13));
        } else {
          REQUIRE(std::abs(reduced[m1][m2]) == Approx(0.0).margin(1e-13));
        }
      }
    }
  }
}

TEST_CASE("Raw STF normalization follows fixed Cartesian contractions", "[gra::spin][pole-ls]") {
  REQUIRE(gra::spin::RawSTFCouplingNormalization(1, 1, 1) == Approx(std::sqrt(2.0)));
  REQUIRE(gra::spin::RawSTFCouplingNormalization(1, 2, 2) == Approx(std::sqrt(3.0 / 2.0)));
  REQUIRE(gra::spin::RawSTFCouplingNormalization(2, 2, 2) == Approx(std::sqrt(7.0 / 12.0)));
  REQUIRE(gra::spin::RawLSOperatorNormalization(0, 0, 2, 2, 0) == Approx(std::sqrt(2.0 / 3.0)));
  REQUIRE(gra::spin::RawLSOperatorNormalization(1, 2, 2, 1, 2) == Approx(std::sqrt(15.0 / 4.0)));
  REQUIRE_THROWS_AS(gra::spin::RawSTFCouplingNormalization(85, 85, 0), std::invalid_argument);
}

// Keep physical two-body C selection in decays and fusion, with separate crossed vertices
TEST_CASE("Pole LS C parity distinguishes physical pairs from crossed residues",
          "[gra::spin][pole-ls][crossing][parity][regression]") {
  const auto                           mother     = Particle(20223, 1, 1, 1);
  const auto                           proton     = ParticleX2(2212, 1, 1, 0);
  const auto                           antiproton = ParticleX2(-2212, 1, -1, 0);
  const std::vector<gra::spin::LSTerm> terms{{1, 0, 1.0}};
  // The singlet P wave has C = -1, whereas the physical mother has C = +1
  for (const auto context : {gra::spin::VertexContext::Auto, gra::spin::VertexContext::CrossedBeamLeg,
                             gra::spin::VertexContext::SubTUChannelExchange}) {
    CHECK_THROWS_AS(gra::spin::PreparePoleLS(mother, proton, antiproton, terms, 1.0, false, true, true, context),
                    std::invalid_argument);
    if (context == gra::spin::VertexContext::Auto) {
      CHECK_THROWS_AS(gra::spin::PreparePoleLS(mother, proton, antiproton, terms, 1.0, true, true, true, context),
                      std::invalid_argument);
    } else {
      const auto vertex = gra::spin::PreparePoleLS(mother, proton, antiproton, terms, 1.0, true, true, true, context);
      const auto reference =
          gra::spin::PreparePoleLS(mother, proton, antiproton, terms, 1.0, true, false, true, context);
      const auto reduced = gra::spin::PoleLSReduced(vertex, 1.0);
      CHECK(reduced.FrobNorm2() > 0.0);
      RequireMatrixNear(reduced, gra::spin::PoleLSReduced(reference, 1.0), 2.0e-12);
    }
  }
}

TEST_CASE("A nonphoton crossed pole is one angle-free Regge residue", "[gra::spin][pole-ls][continuum][regression]") {
  const auto first_scalar  = Particle(211, 0, -1);
  const auto second_scalar = Particle(-211, 0, -1);

  for (std::size_t spin = 0; spin <= 6; ++spin) {
    CAPTURE(spin);
    const int               parity   = spin % 2 == 0 ? 1 : -1;
    const auto              exchange = Particle(900100 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const gra::spin::LSTerm term{
        spin, 0, std::polar(0.72 + 0.03 * static_cast<double>(spin), -0.18 + 0.05 * static_cast<double>(spin))};
    const auto vertex  = gra::spin::PreparePoleLS(exchange, first_scalar, second_scalar, {term}, 1.0, true, false, true,
                                                  gra::spin::VertexContext::SubTUChannelExchange);
    const auto reduced = gra::spin::PoleLSHelicity(vertex, vertex.Lambda);
    const auto crossed = gra::spin::EvaluateCrossedPole(vertex);
    REQUIRE(crossed.helicity.Jz_values == std::vector<double>{0.0});
    REQUIRE(crossed.frame.size_row() == 1);
    REQUIRE(crossed.frame.size_col() == 1);
    RequireComplexNear(crossed.frame[0][0], reduced.T[0][0], 2.0e-13);
  }
}

// Keep the timelike fixed-pole Wigner theorem separate from crossed Regge use
TEST_CASE("A physical spin-two pole retains its Jacob-Wick angular frame",
          "[gra::spin][pole-ls][physical-pole][wigner][regression]") {
  const auto tensor   = Particle(910002, 2, 1);
  const auto first    = Particle(211, 0, -1);
  const auto second   = Particle(-211, 0, -1);
  const auto vertex   = gra::spin::PreparePoleLS(tensor, first, second, {{2, 0, {0.73, -0.19}}}, 1.0, true, false, true,
                                                 gra::spin::VertexContext::Auto, 0.0, false);
  const auto helicity = gra::spin::PoleLSHelicity(vertex, 1.0);
  const std::size_t zero = gra::spin::SpinProjectionIndex(0.0, helicity.J, "physical spin-two zero projection");
  const gra::M4Vec  axis(0.0, 0.0, 1.0, 1.2);
  const auto        frame = [&](const double theta) {
    const gra::M4Vec final(std::sin(theta), 0.0, std::cos(theta), 1.4);
    return gra::spin::VirtualSubchannelFrame(helicity, final, axis, false);
  };
  const auto forward    = frame(0.0);
  const auto transverse = frame(0.5 * gra::math::PI);
  const auto node       = frame(std::acos(1.0 / std::sqrt(3.0)));
  REQUIRE(std::abs(forward[0][zero]) > 0.0);
  RequireComplexNear(transverse[0][zero], -0.5 * forward[0][zero], 2.0e-12);
  RequireComplexNear(node[0][zero], 0.0, 2.0e-12);
}

// Check the crossed second daughter row and phase for non-scalar helicities
TEST_CASE("Reduced crossed residues retain the Jacob-Wick second-leg map",
          "[gra::spin][pole-ls][continuum][crossing][regression]") {
  const auto exchange  = Particle(910102, 2, 1);
  const auto first     = ParticleX2(910011, 2, 1);
  const auto second    = ParticleX2(910012, 2, 1);
  const auto operators = gra::spin::CanonicalPoleOperators(exchange, first, second, true, false, false,
                                                           gra::spin::VertexContext::SubTUChannelExchange);
  REQUIRE_FALSE(operators.empty());
  std::vector<gra::spin::LSTerm> terms;
  terms.reserve(operators.size());
  for (const auto& i : indices(operators)) {
    terms.push_back({operators[i].coupling.l, operators[i].coupling.two_s,
                     std::polar(0.41 + 0.03 * static_cast<double>(i), -0.37 + 0.11 * static_cast<double>(i))});
  }
  const auto vertex    = gra::spin::PreparePoleLS(exchange, first, second, terms, 1.0, true, false, false,
                                                  gra::spin::VertexContext::SubTUChannelExchange, 0.0, false);
  const auto auxiliary = gra::spin::PoleLSHelicity(vertex, 1.0);
  const auto physical  = gra::spin::EvaluateCrossedPole(vertex);
  REQUIRE(physical.frame.size_col() == 1);
  REQUIRE(physical.frame.size_row() == auxiliary.lambda_values.size_row());
  for (std::size_t row = 0; row < auxiliary.lambda_values.size_row(); ++row) {
    const double lambda1 = auxiliary.lambda_values[row][0];
    const double lambda2 = auxiliary.lambda_values[row][1];
    std::size_t  crossed = auxiliary.lambda_values.size_row();
    for (std::size_t candidate = 0; candidate < auxiliary.lambda_values.size_row(); ++candidate) {
      if (std::abs(auxiliary.lambda_values[candidate][0] - lambda1) < 1.0e-12 &&
          std::abs(auxiliary.lambda_values[candidate][1] + lambda2) < 1.0e-12) {
        crossed = candidate;
        break;
      }
    }
    REQUIRE(crossed < auxiliary.lambda_values.size_row());
    const std::size_t i1    = auxiliary.lambda_idx[crossed][0];
    const std::size_t i2    = auxiliary.lambda_idx[crossed][1];
    const double      phase = gra::spin::JacobWickSecondLegReversalPhase(auxiliary.s2, lambda2);
    RequireComplexNear(physical.frame[row][0], phase * auxiliary.T[i1][i2], 3.0e-12);
  }
}

TEST_CASE("MP and XP share every reduced scalar-pair pole kernel", "[gra::spin][pole-ls][continuum]") {
  gra::MDecayBranch first;
  gra::MDecayBranch second;
  first.p   = Particle(211, 0, -1);
  second.p  = Particle(-211, 0, -1);
  first.p4  = gra::M4Vec(0.31, -0.18, 0.42, 0.62);
  second.p4 = gra::M4Vec(-0.31, 0.18, -0.42, 0.62);
  const gra::M4Vec upper(0.13, -0.07, 0.21, 0.18);
  const gra::M4Vec lower(-0.09, 0.11, -0.17, 0.16);

  for (std::size_t spin = 0; spin <= 4; ++spin) {
    CAPTURE(spin);
    const int               parity   = spin % 2 == 0 ? 1 : -1;
    const auto              exchange = Particle(800100 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const gra::spin::LSTerm term{spin, 0, 1.0};
    gra::ReggeContinuumPole cache;
    cache.pole_operator[0] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, first.p, second.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    cache.pole_operator[1] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, second.p, first.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));

    const auto mp = gra::rspin::PairKernel(cache, first, second, upper, lower);
    REQUIRE(mp.size_row() == 1);
    REQUIRE(mp.size_col() == 1);
    REQUIRE(mp.IsFinite());
    const auto upper_pole = gra::spin::EvaluateCrossedPole(cache.pole_operator[0].Pole());
    const auto lower_pole = gra::spin::EvaluateCrossedPole(cache.pole_operator[1].Pole());
    RequireComplexNear(mp[0][0], upper_pole.frame[0][0] * lower_pole.frame[0][0], 2.0e-13);

    auto       rotated_first  = first;
    auto       rotated_second = second;
    gra::M4Vec rotated_upper  = upper;
    gra::M4Vec rotated_lower  = lower;
    rotated_first.p4.RotateZ(0.73);
    rotated_second.p4.RotateZ(0.73);
    rotated_upper.RotateZ(0.73);
    rotated_lower.RotateZ(0.73);
    const auto rotated = gra::rspin::PairKernel(cache, rotated_first, rotated_second, rotated_upper, rotated_lower);
    RequireMatrixNear(rotated, mp, 2.0e-13);
  }
}

// Check the full physical-pole pair kernel independently of Regge reduction
TEST_CASE("MP and XP share every physical scalar-pair pole kernel", "[gra::spin][pole-ls][physical-pole][continuum]") {
  gra::MDecayBranch first;
  gra::MDecayBranch second;
  first.p   = Particle(211, 0, -1);
  second.p  = Particle(-211, 0, -1);
  first.p4  = gra::M4Vec(0.31, -0.18, 0.42, 0.62);
  second.p4 = gra::M4Vec(-0.31, 0.18, -0.42, 0.62);
  const gra::M4Vec upper(0.13, -0.07, 0.21, 0.18);
  const gra::M4Vec lower(-0.09, 0.11, -0.17, 0.16);

  for (std::size_t spin = 0; spin <= 4; ++spin) {
    CAPTURE(spin);
    const int               parity   = spin % 2 == 0 ? 1 : -1;
    const auto              exchange = Particle(800200 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const gra::spin::LSTerm term{spin, 0, 1.0};
    gra::ReggeContinuumPole pole;
    pole.pole_operator[0] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, first.p, second.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    pole.pole_operator[1] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, second.p, first.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    for (auto& residue : pole.pole_operator) {
      auto vertex                    = residue.Pole();
      vertex.helicity.exchange_basis = gra::ExchangeBasisType::HelicityTransport;
      residue                        = gra::spin::PoleResidue(vertex);
    }

    const auto mp = gra::rspin::PairKernel(pole, first, second, upper, lower);
    REQUIRE(mp.size_row() == 2 * spin + 1);
    REQUIRE(mp.size_col() == 2 * spin + 1);
    REQUIRE(mp.IsFinite());
  }
}

// Check the Feynman-bilinear lower vertex for every representative pole spin
TEST_CASE("Crossed Regge pole vertices use one bilinear lower residue",
          "[gra::spin][pole-ls][continuum][phase][regression]") {
  const auto                 first        = Particle(211, 0, -1);
  const auto                 second       = Particle(-211, 0, -1);
  constexpr double           phase        = 0.43;
  const std::complex<double> phase_factor = std::polar(1.0, phase);

  for (std::size_t spin = 0; spin <= 6; ++spin) {
    CAPTURE(spin);
    const int                  parity   = spin % 2 == 0 ? 1 : -1;
    const auto                 exchange = Particle(810100 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const std::complex<double> coupling =
        std::polar(0.73 + 0.02 * static_cast<double>(spin), 0.19 + 0.07 * static_cast<double>(spin));
    const gra::spin::LSTerm term{spin, 0, coupling};
    const auto upper_vertex = gra::spin::PreparePoleLS(exchange, first, second, {term}, 1.0, true, false, true,
                                                       gra::spin::VertexContext::SubTUChannelExchange);
    const auto lower_vertex = gra::spin::PreparePoleLS(exchange, second, first, {term}, 1.0, true, false, true,
                                                       gra::spin::VertexContext::SubTUChannelExchange);
    const auto upper        = gra::spin::EvaluateCrossedPole(upper_vertex);
    const auto lower        = gra::spin::EvaluateCrossedPole(lower_vertex);
    REQUIRE(upper.helicity.Jz_values == std::vector<double>{0.0});
    REQUIRE(lower.helicity.Jz_values == std::vector<double>{0.0});

    const auto kernel =
        gra::spin::Subchannel(upper, lower, first, second, false, true);
    REQUIRE(kernel.size_row() == 1);
    REQUIRE(kernel.size_col() == 1);
    const std::complex<double> bilinear = upper.frame[0][0] * lower.frame[0][0];
    RequireComplexNear(kernel[0][0], bilinear, 4.0e-13);
    REQUIRE(std::abs(kernel[0][0] - upper.frame[0][0] * std::conj(lower.frame[0][0])) > 1.0e-6);

    auto phased_vertex = lower_vertex;
    phased_vertex.terms[0].coefficient *= phase_factor;
    const auto phased_lower  = gra::spin::EvaluateCrossedPole(phased_vertex);
    const auto phased_kernel = gra::spin::Subchannel(upper, phased_lower, first, second, false, true);
    RequireComplexNear(phased_kernel[0][0], phase_factor * kernel[0][0], 6.0e-13);

    const std::vector<std::complex<double>> upper_source = {std::polar(0.81, 0.13)};
    const std::vector<std::complex<double>> lower_source = {std::polar(0.76, -0.17)};
    const auto projected = gra::spin::ProjectedSubchannel(upper, lower, upper_source, lower_source, first, second, false);
    REQUIRE(projected.size() == 1);
    const std::complex<double> upper_projection = (upper.frame * upper_source)[0];
    const std::complex<double> lower_projection = (lower.frame * lower_source)[0];
    const std::complex<double> expected         = upper_projection * lower_projection;
    RequireComplexNear(projected[0], expected, 6.0e-13);

    {
      const auto metric = gra::spin::ExchangeMetric(gra::ExchangeBasisType::ReducedRegge, {0.0}, static_cast<double>(spin));
      REQUIRE(metric.size_row() == 1);
      REQUIRE(metric.size_col() == 1);
      RequireComplexNear(metric[0][0], 1.0, 2.0e-13);
    }
    const auto gp_metric = gra::spin::ExchangeMetric(gra::ExchangeBasisType::ReggeHelicity, {0.0}, static_cast<double>(spin));
    REQUIRE(gp_metric.size_row() == 1);
    REQUIRE(gp_metric.size_col() == 1);
    RequireComplexNear(gp_metric[0][0], 1.0, 2.0e-13);
  }
}

// Check one full spherical metric and the bilinear lower pole sewing
TEST_CASE("Physical pole lower vertices use one spherical metric and bilinear sewing",
          "[gra::spin][pole-ls][physical-pole][parity][phase][regression]") {
  const auto       first  = Particle(211, 0, -1);
  const auto       second = Particle(-211, 0, -1);
  const gra::M4Vec upper_final(0.29, 0.17, -0.22, 0.72);
  const gra::M4Vec lower_final(-0.16, -0.25, 0.19, 0.70);
  const gra::M4Vec upper_parent(0.18, -0.12, 0.31, 0.65);
  const gra::M4Vec lower_parent(-0.21, 0.14, 0.27, 0.61);

  for (std::size_t spin = 0; spin <= 6; ++spin) {
    CAPTURE(spin);
    const int                  parity   = spin % 2 == 0 ? 1 : -1;
    const auto                 exchange = Particle(810200 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const std::complex<double> coupling =
        std::polar(0.73 + 0.02 * static_cast<double>(spin), 0.19 + 0.07 * static_cast<double>(spin));
    const gra::spin::LSTerm term{spin, 0, coupling};
    auto upper_vertex = gra::spin::PreparePoleLS(exchange, first, second, {term}, 1.0, true, false, true,
                                                 gra::spin::VertexContext::SubTUChannelExchange);
    auto lower_vertex = gra::spin::PreparePoleLS(exchange, second, first, {term}, 1.0, true, false, true,
                                                 gra::spin::VertexContext::SubTUChannelExchange);
    upper_vertex.helicity.exchange_basis = gra::ExchangeBasisType::HelicityTransport;
    lower_vertex.helicity.exchange_basis = gra::ExchangeBasisType::HelicityTransport;
    const auto  upper       = gra::spin::EvaluatePoleSubvertex(upper_vertex, upper_final, upper_parent, false);
    const auto  lower       = gra::spin::EvaluatePoleSubvertex(lower_vertex, lower_final, lower_parent, true);
    const auto& projections = lower.helicity.Jz_values;
    REQUIRE(projections.size() == 2 * spin + 1);

    const auto mp_metric = gra::spin::ExchangeMetric(gra::ExchangeBasisType::HelicityTransport, projections, static_cast<double>(spin));
    const auto xp_metric = gra::spin::ExchangeMetric(gra::ExchangeBasisType::HelicityTransport, projections, static_cast<double>(spin));
    for (const auto& row : indices(projections)) {
      for (const auto& column : indices(projections)) {
        const bool                 crossed = std::abs(projections[row] + projections[column]) < 1.0e-13;
        const std::complex<double> expected =
            crossed ? gra::spin::JacobWickSecondLegReversalPhase(static_cast<double>(spin), projections[row]) : 0.0;
        RequireComplexNear(mp_metric[row][column], expected, 2.0e-13);
        RequireComplexNear(xp_metric[row][column], expected, 2.0e-13);
      }
    }

    const gra::M4Vec lower_relative = lower_final + lower_parent * 0.5;
    const auto       lower_raw = gra::spin::VirtualSubchannelFrame(lower.helicity, lower_relative, lower_parent, true);
    const auto       lower_expected = lower_raw * mp_metric;
    REQUIRE(lower.frame.size_row() == lower_expected.size_row());
    REQUIRE(lower.frame.size_col() == lower_expected.size_col());
    REQUIRE(lower.frame.size_row() == 1);
    bool has_asymmetric_columns = spin == 0;
    for (const auto& column : indices(projections)) {
      RequireComplexNear(lower.frame[0][column], lower_expected[0][column], 3.0e-13);
      const auto opposite = std::find_if(
          projections.cbegin(), projections.cend(),
          [projection = projections[column]](const double value) { return std::abs(value + projection) < 1.0e-13; });
      REQUIRE(opposite != projections.cend());
      const std::size_t          opposite_column = static_cast<std::size_t>(opposite - projections.cbegin());
      const std::complex<double> component_expected =
          gra::spin::JacobWickSecondLegReversalPhase(static_cast<double>(spin), projections[column]) *
          lower_raw[0][opposite_column];
      RequireComplexNear(lower.frame[0][column], component_expected, 3.0e-13);
      has_asymmetric_columns =
          has_asymmetric_columns || std::abs(lower_raw[0][column] - lower_raw[0][opposite_column]) > 1.0e-6;
    }
    REQUIRE(has_asymmetric_columns);

    const auto kernel =
        gra::spin::Subchannel(upper, lower, first, second, false, true);
    REQUIRE(kernel.size_row() == projections.size() * projections.size());
    REQUIRE(kernel.size_col() == 1);
    bool distinguishes_hermitian = false;
    for (const auto& upper_spin : indices(projections)) {
      for (const auto& lower_spin : indices(projections)) {
        const std::size_t          row       = upper_spin * projections.size() + lower_spin;
        const std::complex<double> expected  = upper.frame[0][upper_spin] * lower.frame[0][lower_spin];
        const std::complex<double> hermitian = upper.frame[0][upper_spin] * std::conj(lower.frame[0][lower_spin]);
        RequireComplexNear(kernel[row][0], expected, 4.0e-13);
        distinguishes_hermitian = distinguishes_hermitian || std::abs(kernel[row][0] - hermitian) > 1.0e-6;
      }
    }
    REQUIRE(distinguishes_hermitian);

    std::vector<std::complex<double>> upper_source(projections.size());
    std::vector<std::complex<double>> lower_source(projections.size());
    for (const auto& projection : indices(projections)) {
      upper_source[projection] =
          std::polar(0.81 + 0.03 * static_cast<double>(projection), 0.13 * static_cast<double>(projection + 1));
      lower_source[projection] =
          std::polar(0.76 + 0.04 * static_cast<double>(projection), -0.17 * static_cast<double>(projection + 1));
    }
    const auto projected = gra::spin::ProjectedSubchannel(upper, lower, upper_source, lower_source, first, second, false);
    REQUIRE(projected.size() == 1);
    const std::complex<double> upper_projection = (upper.frame * upper_source)[0];
    const std::complex<double> lower_projection = (lower.frame * lower_source)[0];
    const std::complex<double> expected         = upper_projection * lower_projection;
    RequireComplexNear(projected[0], expected, 6.0e-13);
    REQUIRE(std::abs(projected[0] - upper_projection * std::conj(lower_projection)) > 1.0e-6);
  }
}

// Check physical photon transport and its single lower spherical metric
TEST_CASE("Photon crossing supports a restricted transverse parent basis",
          "[gra::spin][pole-ls][continuum][photon][regression]") {
  const auto photon = Particle(22, 1, -1, -1);
  const auto first  = Particle(211, 0, -1);
  const auto second = Particle(-211, 0, -1);
  const auto vertex = gra::spin::PreparePoleLS(photon, second, first, {{1, 0, {0.71, 0.23}}}, 1.0, true, false, true,
                                               gra::spin::VertexContext::SubTUChannelExchange);
  REQUIRE(vertex.helicity.exchange_basis == gra::ExchangeBasisType::HelicityTransport);
  REQUIRE(vertex.helicity.Jz_values.size() == 2);
  CHECK(vertex.helicity.Jz_values[0] == Approx(-1.0));
  CHECK(vertex.helicity.Jz_values[1] == Approx(1.0));

  const gra::M4Vec final(-0.16, -0.25, 0.19, 0.70);
  const gra::M4Vec parent(-0.21, 0.14, 0.27, 0.61);
  const auto       lower = gra::spin::EvaluatePoleSubvertex(vertex, final, parent, true);
  REQUIRE(lower.frame.size_col() == 2);
  const auto       metric = gra::spin::ExchangeMetric(gra::ExchangeBasisType::HelicityTransport, lower.helicity.Jz_values, 1.0);
  const gra::M4Vec relative = final + parent * 0.5;
  const auto       raw      = gra::spin::VirtualSubchannelFrame(lower.helicity, relative, parent, true);
  const auto       expected = raw * metric;
  REQUIRE(lower.frame.size_row() == 1);
  for (const auto& column : indices(lower.helicity.Jz_values)) {
    RequireComplexNear(lower.frame[0][column], expected[0][column], 3.0e-13);
  }

  constexpr double rotation       = 0.57;
  gra::M4Vec       rotated_final  = final;
  gra::M4Vec       rotated_parent = parent;
  rotated_final.RotateZ(rotation);
  rotated_parent.RotateZ(rotation);
  const auto rotated_lower = gra::spin::EvaluatePoleSubvertex(vertex, rotated_final, rotated_parent, true);
  const auto upper         = gra::spin::EvaluatePoleSubvertex(vertex, final, parent, false);
  const auto rotated_upper = gra::spin::EvaluatePoleSubvertex(vertex, rotated_final, rotated_parent, false);
  for (const auto& column : indices(lower.helicity.Jz_values)) {
    const double m = lower.helicity.Jz_values[column];
    RequireComplexNear(rotated_lower.frame[0][column], std::exp(gra::math::zi * m * rotation) * lower.frame[0][column],
                       4.0e-13);
    RequireComplexNear(rotated_upper.frame[0][column], std::exp(-gra::math::zi * m * rotation) * upper.frame[0][column],
                       4.0e-13);
  }
}

// Check that every nonphoton fixed pole uses the reduced ladder contraction
TEST_CASE("MP and XP ladder boundaries close reduced pair kernels",
          "[gra::spin][pole-ls][continuum][rotation][regression]") {
  gra::MDecayBranch first;
  gra::MDecayBranch second;
  first.p   = Particle(211, 0, -1);
  second.p  = Particle(-211, 0, -1);
  first.p4  = gra::M4Vec(0.36, -0.24, 0.31, 0.64);
  second.p4 = gra::M4Vec(-0.36, 0.24, -0.31, 0.64);
  const gra::M4Vec             upper_transfer(0.53, -0.41, 0.21, 0.18);
  const gra::M4Vec             lower_transfer(-0.47, 0.38, -0.17, 0.16);
  const gra::spin::ForwardSpec forward  = {gra::ForwardVertexMode::HelicityResidue, 1.0,
                                           gra::ExchangeBasisType::ReducedRegge};
  constexpr double             rotation = 0.63;

  for (std::size_t spin = 0; spin <= 6; ++spin) {
    CAPTURE(spin);
    const int               parity   = spin % 2 == 0 ? 1 : -1;
    const auto              exchange = Particle(820100 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const gra::spin::LSTerm term{
        spin, 0, std::polar(0.79 + 0.02 * static_cast<double>(spin), -0.16 + 0.05 * static_cast<double>(spin))};
    gra::ReggeContinuumPole pole;
    pole.pole_operator[0] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, first.p, second.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    pole.pole_operator[1] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, second.p, first.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));

    for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
      CAPTURE(model);
      const auto upper_basis = gra::rspin::LocalBasis(pole, 0);
      const auto lower_basis = gra::rspin::LocalBasis(pole, 1);
      REQUIRE(upper_basis == std::vector<double>{0.0});
      REQUIRE(lower_basis == std::vector<double>{0.0});

      // Evaluate one rotated public boundary and pair-kernel contraction
      const auto evaluate = [&](const double angle) {
        gra::MDecayBranch rotated_first  = first;
        gra::MDecayBranch rotated_second = second;
        gra::M4Vec        rotated_upper  = upper_transfer;
        gra::M4Vec        rotated_lower  = lower_transfer;
        rotated_first.p4.RotateZ(angle);
        rotated_second.p4.RotateZ(angle);
        rotated_upper.RotateZ(angle);
        rotated_lower.RotateZ(angle);
        const auto kernel = gra::rspin::PairKernel(pole, rotated_first, rotated_second, rotated_upper, -rotated_lower);
        const auto upper  = gra::spin::ExchangeBoundary(upper_basis, rotated_upper, false, forward);
        const auto lower  = gra::spin::ExchangeBoundary(lower_basis, rotated_lower, true, forward);
        REQUIRE(upper == std::vector<std::complex<double>>{1.0});
        REQUIRE(lower == std::vector<std::complex<double>>{1.0});
        const std::vector<gra::MMatrix<std::complex<double>>> kernels = {kernel};
        const std::vector<gra::MMatrix<std::complex<double>>> metrics;
        return gra::ladder::Contract(upper, kernels, metrics, lower);
      };

      const std::complex<double> reference = evaluate(0.0);
      const std::complex<double> rotated   = evaluate(rotation);
      REQUIRE(std::isfinite(std::real(reference)));
      REQUIRE(std::isfinite(std::imag(reference)));
      REQUIRE(std::abs(reference) > 0.0);
      const double tolerance = 2.0e-11 * std::max({1.0, std::abs(reference), std::abs(rotated)});
      RequireComplexNear(rotated, reference, tolerance);
    }
  }

  const auto analytic_lower = gra::spin::ExchangeBoundary({0.0}, lower_transfer, true, forward);
  REQUIRE(analytic_lower.size() == 1);
  RequireComplexNear(analytic_lower[0], 1.0, 2.0e-13);
}

// Check full physical-pole boundary transport in a non-axial ladder
TEST_CASE("MP and XP physical-pole ladder boundaries close non-axial kernels",
          "[gra::spin][pole-ls][physical-pole][rotation][regression]") {
  gra::MDecayBranch first;
  gra::MDecayBranch second;
  first.p   = Particle(211, 0, -1);
  second.p  = Particle(-211, 0, -1);
  first.p4  = gra::M4Vec(0.36, -0.24, 0.31, 0.64);
  second.p4 = gra::M4Vec(-0.36, 0.24, -0.31, 0.64);
  const gra::M4Vec             upper_transfer(0.53, -0.41, 0.21, 0.18);
  const gra::M4Vec             lower_transfer(-0.47, 0.38, -0.17, 0.16);
  const gra::spin::ForwardSpec forward  = {gra::ForwardVertexMode::HelicityResidue, 1.0,
                                           gra::ExchangeBasisType::HelicityTransport};
  constexpr double             rotation = 0.63;

  for (std::size_t spin = 0; spin <= 6; ++spin) {
    CAPTURE(spin);
    const int               parity   = spin % 2 == 0 ? 1 : -1;
    const auto              exchange = Particle(820200 + static_cast<int>(spin), static_cast<int>(spin), parity);
    const gra::spin::LSTerm term{
        spin, 0, std::polar(0.79 + 0.02 * static_cast<double>(spin), -0.16 + 0.05 * static_cast<double>(spin))};
    gra::ReggeContinuumPole pole;
    pole.pole_operator[0] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, first.p, second.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    pole.pole_operator[1] = gra::spin::PoleResidue(gra::spin::PreparePoleLS(
        exchange, second.p, first.p, {term}, 1.0, true, false, true, gra::spin::VertexContext::SubTUChannelExchange));
    for (auto& residue : pole.pole_operator) {
      auto vertex                    = residue.Pole();
      vertex.helicity.exchange_basis = gra::ExchangeBasisType::HelicityTransport;
      residue                        = gra::spin::PoleResidue(vertex);
    }

    for (const auto model : {gra::ReggeProductionModel::MP, gra::ReggeProductionModel::XP}) {
      CAPTURE(model);
      const auto upper_basis = gra::rspin::LocalBasis(pole, 0);
      const auto lower_basis = gra::rspin::LocalBasis(pole, 1);
      REQUIRE(upper_basis.size() == 2 * spin + 1);
      REQUIRE(lower_basis.size() == 2 * spin + 1);

      // Evaluate one rotated public boundary and pair-kernel contraction
      const auto evaluate = [&](const double angle) {
        gra::MDecayBranch rotated_first  = first;
        gra::MDecayBranch rotated_second = second;
        gra::M4Vec        rotated_upper  = upper_transfer;
        gra::M4Vec        rotated_lower  = lower_transfer;
        rotated_first.p4.RotateZ(angle);
        rotated_second.p4.RotateZ(angle);
        rotated_upper.RotateZ(angle);
        rotated_lower.RotateZ(angle);
        const auto kernel = gra::rspin::PairKernel(pole, rotated_first, rotated_second, rotated_upper, -rotated_lower);
        const auto upper  = gra::spin::ExchangeBoundary(upper_basis, rotated_upper, false, forward);
        const auto lower  = gra::spin::ExchangeBoundary(lower_basis, rotated_lower, true, forward);
        const auto upper_expected = gra::spin::ExchangeHelicityResidue(upper_basis, rotated_upper, false, forward);
        const auto lower_expected = gra::spin::ExchangeHelicityResidue(lower_basis, rotated_lower, true, forward);
        REQUIRE(upper.size() == upper_expected.size());
        REQUIRE(lower.size() == lower_expected.size());
        for (const auto& projection : indices(upper_basis)) {
          RequireComplexNear(upper[projection], upper_expected[projection], 2.0e-13);
          RequireComplexNear(lower[projection], lower_expected[projection], 2.0e-13);
        }
        const std::vector<gra::MMatrix<std::complex<double>>> kernels = {kernel};
        const std::vector<gra::MMatrix<std::complex<double>>> metrics;
        return gra::ladder::Contract(upper, kernels, metrics, lower);
      };

      const std::complex<double> reference = evaluate(0.0);
      const std::complex<double> rotated   = evaluate(rotation);
      REQUIRE(std::isfinite(std::real(reference)));
      REQUIRE(std::isfinite(std::imag(reference)));
      REQUIRE(std::abs(reference) > 0.0);
      const double tolerance = 2.0e-11 * std::max({1.0, std::abs(reference), std::abs(rotated)});
      RequireComplexNear(rotated, reference, tolerance);
    }
  }
}

// Keep the analytically continued GP vertex as a bilinear contraction
TEST_CASE("Analytic crossed lower residues remain bilinear", "[gra::spin][pole-ls][GP][continuum][regression]") {
  gra::HELMatrix upper_input;
  gra::HELMatrix lower_input;
  gra::gpom::InitCrossed(upper_input, 0.0, 0.0, 0, "analytic upper dual test");
  gra::gpom::InitCrossed(lower_input, 0.0, 0.0, 0, "analytic lower dual test");
  upper_input.coupling_basis = gra::CouplingBasis::Helicity;
  lower_input.coupling_basis = gra::CouplingBasis::Helicity;
  upper_input.T[0][0]        = std::polar(0.83, 0.31);
  lower_input.T[0][0]        = std::polar(1.17, -0.46);
  upper_input.T_set[0][0]    = true;
  lower_input.T_set[0][0]    = true;
  const auto upper           = gra::gpom::Crossed(upper_input, 1.23, gra::M4Vec{}, false);
  const auto lower           = gra::gpom::Crossed(lower_input, 1.67, gra::M4Vec{}, false);
  REQUIRE(upper.helicity.UsesReggeDomain());
  REQUIRE(lower.helicity.UsesReggeDomain());

  const auto first  = Particle(211, 0, -1);
  const auto second = Particle(-211, 0, -1);
  const auto kernel =
      gra::spin::Subchannel(upper, lower, first, second, false, true);
  REQUIRE(kernel.size_row() == 1);
  REQUIRE(kernel.size_col() == 1);
  const std::complex<double> bilinear  = upper.frame[0][0] * lower.frame[0][0];
  const std::complex<double> hermitian = upper.frame[0][0] * std::conj(lower.frame[0][0]);
  RequireComplexNear(kernel[0][0], bilinear, 2.0e-13);
  REQUIRE(std::abs(kernel[0][0] - hermitian) > 1.0e-3);

  const std::vector<std::complex<double>> upper_source = {std::polar(0.91, 0.27)};
  const std::vector<std::complex<double>> lower_source = {std::polar(1.08, -0.39)};
  const auto projected = gra::spin::ProjectedSubchannel(upper, lower, upper_source, lower_source, first, second, false);
  REQUIRE(projected.size() == 1);
  const std::complex<double> projected_bilinear = (upper.frame * upper_source)[0] * (lower.frame * lower_source)[0];
  const std::complex<double> projected_hermitian =
      (upper.frame * upper_source)[0] * std::conj((lower.frame * lower_source)[0]);
  RequireComplexNear(projected[0], projected_bilinear, 2.0e-13);
  REQUIRE(std::abs(projected[0] - projected_hermitian) > 1.0e-3);
}

TEST_CASE("Crossed LS evaluates continuous J in the m zero sector", "[gra::spin][pole-ls][GP][continuum]") {
  gra::HELMatrix analytic;
  gra::gpom::InitCrossed(analytic, 0.0, 0.0, 0, "crossed LS GP test");
  analytic.coupling_basis = gra::CouplingBasis::LS;
  analytic.m_ls.resize(1);
  analytic.m_ls[0].Set(2, 0, 1.0);
  gra::gpom::InitCrossedLS(analytic, 2);

  const auto pole = gra::gpom::Crossed(analytic, 2.0, gra::M4Vec{}, false);
  REQUIRE(pole.frame.size_row() == 1);
  REQUIRE(pole.frame.size_col() == 1);
  RequireComplexNear(pole.helicity.T[0][0], std::sqrt(2.0 / 3.0), 2.0e-13);
  RequireComplexNear(pole.frame[0][0], pole.helicity.T[0][0], 2.0e-13);

  const auto continued = gra::gpom::Crossed(analytic, 1.37, gra::M4Vec{}, false);
  REQUIRE(continued.helicity.T.IsFinite());
  RequireComplexNear(continued.frame[0][0], continued.helicity.T[0][0], 2.0e-13);
}

TEST_CASE("Crossed helicity keeps the direct m zero residue", "[gra::spin][pole-ls][GP][continuum]") {
  gra::HELMatrix analytic;
  gra::gpom::InitCrossed(analytic, 0.0, 0.0, 0, "crossed helicity GP test");
  analytic.coupling_basis            = gra::CouplingBasis::Helicity;
  const std::complex<double> residue = std::polar(1.7, 0.23);
  analytic.T[0][0]                   = residue;
  analytic.T_set[0][0]               = true;

  const auto value = gra::gpom::Crossed(analytic, 2.31, gra::M4Vec{}, false);
  REQUIRE(value.frame.size_row() == 1);
  REQUIRE(value.frame.size_col() == 1);
  RequireComplexNear(value.helicity.T[0][0], residue, 2.0e-13);
  RequireComplexNear(value.frame[0][0], residue, 2.0e-13);
}

TEST_CASE("GP crossed residues reject the resonance fusion m basis", "[gra::spin][pole-ls][GP][continuum]") {
  gra::HELMatrix  fusion;
  const gra::json channel  = {{"CP", {true, true}}, {"basis", "helicity"}, {"helicity", {{0, 0, 1.0, 0.0}}}};
  const auto      exchange = Particle(990, 2, 1);
  const auto      scalar   = Particle(900001, 0, 1);
  REQUIRE_THROWS_AS(gra::ParseTwoBodyCouplings(
                        fusion, channel, "GP crossed boundary test", "GP", exchange, {scalar, scalar}, {}, true,
                        gra::spin::VertexContext::SubTUChannelExchange, gra::regge::Signature::Positive, true, 2),
                    std::invalid_argument);
}

TEST_CASE("Raw STF LS operators enforce quantum numbers at setup", "[gra::spin][pole-ls]") {
  const auto              vector       = Particle(993, 1, -1);
  const auto              pseudoscalar = Particle(900001, 0, -1);
  const gra::spin::LSTerm allowed{1, 2, 1.0};
  REQUIRE_NOTHROW(gra::spin::PreparePoleLS(pseudoscalar, vector, vector, {allowed}, 1.0));

  const auto scalar = Particle(900002, 0, 1);
  REQUIRE_THROWS_AS(gra::spin::PreparePoleLS(scalar, vector, vector, {allowed}, 1.0), std::invalid_argument);

  const auto              spin_one = Particle(900003, 1, 1);
  const gra::spin::LSTerm bose_forbidden{1, 0, 1.0};
  REQUIRE_THROWS_AS(gra::spin::PreparePoleLS(spin_one, vector, vector, {bose_forbidden}, 1.0, true, false, false),
                    std::invalid_argument);

  const auto spin_two = Particle(900004, 2, 1);
  const auto tensor   = Particle(995, 2, 1);
  const auto mixed = gra::spin::PreparePoleLS(spin_two, vector, tensor, {{1, 4, {0.7, -0.2}}}, 0.8, true, false, true);
  REQUIRE(mixed.raw_normalization[0] == Approx(std::sqrt(15.0 / 4.0)));
}

TEST_CASE("Pole LS derivative factors preserve absolute coefficients", "[gra::spin][pole-ls]") {
  const auto spin_two = Particle(900005, 2);
  const auto first    = Particle(800001, 0);
  const auto second   = Particle(800002, 0);
  const auto unit     = gra::spin::PreparePoleLS(spin_two, first, second, {{2, 0, 1.0}}, 0.5);
  const auto triple   = gra::spin::PreparePoleLS(spin_two, first, second, {{2, 0, 3.0}}, 0.5);
  const auto disabled = gra::spin::PreparePoleLS(spin_two, first, second, {{2, 0, 1.0}}, 0.5, true, true, true,
                                                 gra::spin::VertexContext::Auto, 0.0, false);

  const auto at_scale          = gra::spin::PoleLSReduced(unit, 0.5);
  const auto at_double         = gra::spin::PoleLSReduced(unit, 1.0);
  const auto at_zero           = gra::spin::PoleLSReduced(unit, 0.0);
  const auto coefficient_three = gra::spin::PoleLSReduced(triple, 0.5);
  const auto disabled_low      = gra::spin::PoleLSReduced(disabled, 0.2);
  const auto disabled_high     = gra::spin::PoleLSReduced(disabled, 4.0);
  REQUIRE(std::abs(at_double[0][0] / at_scale[0][0]) == Approx(4.0));
  REQUIRE(std::abs(coefficient_three[0][0] / at_scale[0][0]) == Approx(3.0));
  REQUIRE(std::abs(at_zero[0][0]) == Approx(0.0).margin(1e-13));
  REQUIRE_FALSE(disabled.derivative_factor);
  REQUIRE(std::abs(disabled_high[0][0] - disabled_low[0][0]) == Approx(0.0).margin(1.0e-13));
}

TEST_CASE("Crossed continuum pole residues do not acquire decay barriers", "[gra::spin][pole-ls][continuum]") {
  const auto spin_two = Particle(995, 2);
  const auto first    = Particle(211, 0, -1);
  const auto second   = Particle(-211, 0, -1);
  const auto vertex   = gra::spin::PreparePoleLS(spin_two, first, second, {{2, 0, 1.0}}, 1.0, true, true, true,
                                                 gra::spin::VertexContext::SubTUChannelExchange);

  const auto low  = gra::spin::PoleLSReduced(vertex, 0.2);
  const auto high = gra::spin::PoleLSReduced(vertex, 4.0);
  REQUIRE(low.size_row() == 1);
  REQUIRE(low.size_col() == 1);
  REQUIRE(std::abs(low[0][0]) > 0.0);
  REQUIRE(std::abs(high[0][0] - low[0][0]) == Approx(0.0).margin(1.0e-13));
}

// Check coherent LS cancellations without rescaling the remaining helicities
TEST_CASE("Canonical pole cancellations are never normalized away", "[gra::spin][pole-ls][normalization]") {
  const auto scalar   = Particle(900015, 0, 1);
  const auto tensor   = Particle(995, 2, 1);
  const auto s_wave   = gra::spin::PreparePoleLS(scalar, tensor, tensor, {{0, 0, 1.0}}, 1.0);
  const auto d_wave   = gra::spin::PreparePoleLS(scalar, tensor, tensor, {{2, 4, 1.0}}, 1.0);
  const auto s_matrix = gra::spin::PoleLSReduced(s_wave, 1.0);
  const auto d_matrix = gra::spin::PoleLSReduced(d_wave, 1.0);
  REQUIRE(std::abs(s_matrix[2][2]) > 1.0e-12);
  REQUIRE(std::abs(d_matrix[2][2]) > 1.0e-12);

  const std::complex<double> cancellation = -s_matrix[2][2] / d_matrix[2][2];
  const auto coherent = gra::spin::PreparePoleLS(scalar, tensor, tensor, {{0, 0, 1.0}, {2, 4, cancellation}}, 1.0);
  const auto combined = gra::spin::PoleLSReduced(coherent, 1.0);
  REQUIRE(std::abs(combined[2][2]) < 1.0e-12);
  REQUIRE(coherent.terms[0].coefficient == std::complex<double>(1.0, 0.0));
  for (std::size_t i = 0; i < combined.size_row(); ++i) {
    for (std::size_t j = 0; j < combined.size_col(); ++j) {
      const auto expected = s_matrix[i][j] + cancellation * d_matrix[i][j];
      REQUIRE(std::real(combined[i][j]) == Approx(std::real(expected)).margin(1.0e-13));
      REQUIRE(std::imag(combined[i][j]) == Approx(std::imag(expected)).margin(1.0e-13));
    }
  }
}

TEST_CASE("Pole LS pole operators support spin-half pairs", "[gra::spin][pole-ls]") {
  const auto tensor = Particle(995, 2);
  const auto proton = ParticleX2(2212, 1);
  const auto vertex = gra::spin::PreparePoleLS(tensor, proton, proton, {{2, 0, 1.0}}, 1.0, true, false, true);
  REQUIRE(vertex.raw_normalization[0] == Approx(2.0 / std::sqrt(3.0)));

  const auto reduced = gra::spin::PoleLSReduced(vertex, 1.0);
  REQUIRE(reduced.size_row() == 2);
  REQUIRE(reduced.size_col() == 2);
  REQUIRE(reduced.IsFinite());
  REQUIRE(reduced.FrobNorm2() == Approx(4.0 / 3.0));
}

// Check zero physical densities at threshold and at a coherent photon LS node
TEST_CASE("Physical pole density preserves threshold zeros and photon LS nodes",
          "[gra::spin][pole-ls][normalization][regression]") {
  const auto tensor       = Particle(995, 2);
  const auto pseudoscalar = Particle(221, 0, -1);
  const auto odd          = gra::spin::PreparePoleLS(pseudoscalar, tensor, tensor, {{1, 2, 1.0}}, 1.0);
  CHECK(gra::spin::LeadingPoleDensity(odd, 0.0) == Approx(0.0).margin(1.0e-24));
  CHECK(gra::spin::LeadingPoleDensity(odd, 0.73) > 0.0);

  const auto       rho         = Particle(113, 1, -1, -1);
  const auto       photon      = Particle(22, 1, -1, -1);
  const auto       scalar      = Particle(991, 0);
  const auto       s_wave      = gra::spin::PreparePoleLS(rho, photon, scalar, {{0, 2, 1.0}}, 1.0);
  const auto       d_wave      = gra::spin::PreparePoleLS(rho, photon, scalar, {{2, 2, 1.0}}, 1.0);
  constexpr double momentum    = 0.73;
  const auto       s_matrix    = gra::spin::PoleLSReduced(s_wave, momentum);
  const auto       d_matrix    = gra::spin::PoleLSReduced(d_wave, momentum);
  const auto       coefficient = -s_matrix[0][0] / d_matrix[0][0];
  const auto       node        = gra::spin::PreparePoleLS(rho, photon, scalar, {{0, 2, 1.0}, {2, 2, coefficient}}, 1.0);
  CHECK(gra::spin::PoleLSReduced(node, momentum).FrobNorm2() < 1.0e-24);
  CHECK(gra::spin::LeadingPoleDensity(node, momentum) == Approx(0.0).margin(1.0e-24));
  CHECK(gra::spin::LeadingPoleDensity(node, momentum + 0.1) > 0.0);
}
