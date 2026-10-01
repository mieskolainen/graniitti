// Check soft CEP fusion algebra with explicit couplings and no model cards
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <complex>
#include <vector>

#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "catch.hpp"

using gra::aux::indices;

namespace {

// Define a test particle by its spin, parity and charge conjugation
gra::MParticle Particle(int pdg, int spin) {
  gra::MParticle particle;
  particle.pdg    = pdg;
  particle.spinX2 = 2 * spin;
  particle.P = particle.C = 1;
  return particle;
}

// Compare complete complex amplitudes without discarding their phases
void CheckAmplitude(const gra::HelAmp& value, const gra::HelAmp& expected) {
  REQUIRE(value.size_row() == expected.size_row());
  REQUIRE(value.size_col() == expected.size_col());
  REQUIRE(value.IsFinite());
  CHECK((value - expected).FrobNorm2() <= 1.0e-22 * expected.FrobNorm2() + 1.0e-26);
}

// Check coupling linearity and the norm identity for independently evaluated amplitudes
template <typename Evaluate>
void CheckCoherence(const std::vector<gra::spin::LSTerm>& terms, Evaluate evaluate) {
  auto second = terms;
  for (const auto i : indices(second)) { second[i].coefficient = std::polar(0.7, -0.31 * (i + 1)); }
  const auto first_amp  = evaluate(terms);
  const auto second_amp = evaluate(second);
  REQUIRE(first_amp.FrobNorm2() > 0.0);
  REQUIRE(second_amp.FrobNorm2() > 0.0);
  std::array<gra::HelAmp, 2> amplitudes;
  for (const auto i : indices(amplitudes)) {
    const std::complex<double> phase(0.0, i == 0 ? 1.0 : -1.0);
    auto                       combined = terms;
    for (const auto j : indices(combined)) { combined[j].coefficient += phase * second[j].coefficient; }
    amplitudes[i] = evaluate(combined);
    CheckAmplitude(amplitudes[i], first_amp + second_amp * phase);
  }
  CHECK(amplitudes[0].FrobNorm2() + amplitudes[1].FrobNorm2() ==
        Approx(2.0 * (first_amp.FrobNorm2() + second_amp.FrobNorm2())).epsilon(1.0e-11));
}

}  // namespace

// Exercise every allowed LS operator without assuming the contents of a fitted tune
TEST_CASE("Soft CEP fusion preserves complex couplings", "[softcep][normalization][phase][rotation]") {
  const auto exchange = Particle(800000, 2);
  for (const int spin : {0, 2}) {
    CAPTURE(spin);
    const auto resonance = Particle(800001, spin);
    const auto operators = gra::spin::CanonicalPoleOperators(resonance, exchange, exchange);
    REQUIRE_FALSE(operators.empty());
    std::vector<gra::spin::LSTerm> terms;
    for (const auto i : indices(operators)) {
      terms.push_back({operators[i].coupling.l, operators[i].coupling.two_s, std::polar(0.4, 0.23 * (i + 1))});
    }
    gra::LORENTZSCALAR lts;
    lts.process.MP_FRAME = "CM";
    lts.q1_in_X          = gra::M4Vec(0.2, -0.3, 0.7, 0.0);
    for (const auto fusion : {gra::mpom::Fusion, gra::xpom::Fusion}) {
      CAPTURE(fusion == gra::mpom::Fusion ? "MP" : "XP");
      CheckCoherence(terms, [&](const auto& coefficients) {
        const auto pole      = gra::spin::PreparePoleLS(resonance, exchange, exchange, coefficients, 1.0);
        const auto amplitude = fusion(lts, pole);
        auto       rotated   = lts;
        rotated.q1_in_X.RotateZ(0.61);
        auto expected = amplitude;
        for (std::size_t row = 0; row < expected.size_row(); ++row) {
          const double lambda = pole.helicity.lambda_values[row][0] - pole.helicity.lambda_values[row][1];
          for (std::size_t col = 0; col < expected.size_col(); ++col) {
            expected[row][col] *= std::polar(1.0, -(pole.helicity.Jz_values[col] + lambda) * 0.61);
          }
        }
        CheckAmplitude(fusion(rotated, pole), expected);
        return amplitude;
      });
    }
    const auto pole = gra::spin::PreparePoleLS(resonance, exchange, exchange, terms, 1.0);
    CheckAmplitude(gra::mpom::Fusion(lts, pole), gra::xpom::Fusion(lts, pole));
    for (const auto alpha : {std::array{2.0, 2.0}, std::array{1.1, 0.9}}) {
      CAPTURE(alpha);
      CheckCoherence(terms, [&](const auto& coefficients) {
        gra::HELMatrix hel;
        hel.coupling_basis  = gra::CouplingBasis::LS;
        hel.analytic_Lambda = 1.0;
        hel.analytic_MMAX   = exchange.spinX2 / 2;
        for (const auto& term : coefficients) { hel.alpha_ls.Set(term.l, term.two_s, term.coefficient); }
        gra::gpom::InitResonanceLS(hel, spin, exchange.spinX2, exchange.spinX2);
        auto source         = gra::gpom::detail::PoleSource(hel, false, false);
        source.basis.alpha1 = alpha[0];
        source.basis.alpha2 = alpha[1];
        return gra::gpom::Fusion(hel, source, lts.q1_in_X.P3mod(), true);
      });
    }
  }
}
