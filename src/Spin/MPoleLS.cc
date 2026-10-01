// Shared physical pole LS algebra and immutable continuum vertices
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityNorm.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::spin {

namespace {

constexpr double kTolerance = 1e-12;

// Largest irreducible rank supported by the finite Wigner factorial backend
constexpr std::size_t kMaximumSTFRank = 84;

// Compute one finite integer spin encoded by spinX2
std::size_t IntegerSpin(const MParticle& particle, const std::string& context) {
  if (particle.spinX2 < 0 || particle.spinX2 % 2 != 0) {
    throw std::invalid_argument(context + ": the pole mother requires integer spin");
  }
  const std::size_t rank = static_cast<std::size_t>(particle.spinX2 / 2);
  if (rank > kMaximumSTFRank) {
    throw std::invalid_argument(context + ": STF rank exceeds the supported Wigner range");
  }
  return rank;
}

// Compute one supported non-negative twice-spin label
int ParticleSpinX2(const MParticle& particle, const std::string& context) {
  if (particle.spinX2 < 0 || particle.spinX2 > 2 * static_cast<int>(kMaximumSTFRank)) {
    throw std::invalid_argument(context + ": unsupported spin label");
  }
  return particle.spinX2;
}

// Convert one supported rank to the signed representation used in setup
int SignedRank(std::size_t rank, const std::string& context) {
  if (rank > kMaximumSTFRank) {
    throw std::invalid_argument(context + ": STF rank exceeds the supported Wigner range");
  }
  return static_cast<int>(rank);
}

// Compute whether three integer ranks satisfy the SU(2) triangle rule
bool Triangle(std::size_t a, std::size_t b, std::size_t c) { return c <= a + b && c >= ((a > b) ? a - b : b - a); }

// Compute one normalized spin-j polarization in the symmetric occupation basis
std::map<std::array<int, 3>, double> PolarizationOccupations(int j, int m) {
  std::map<std::array<int, 3>, double> state;
  state[{j, 0, 0}] = 1.0;
  for (int current_m = j; current_m > m; --current_m) {
    const double denominator = std::sqrt(static_cast<double>((j + current_m) * (j - current_m + 1)));
    std::map<std::array<int, 3>, double> lowered;
    for (const auto& [occupation, coefficient] : state) {
      const int plus  = occupation[0];
      const int zero  = occupation[1];
      const int minus = occupation[2];
      if (plus > 0) {
        lowered[{plus - 1, zero + 1, minus}] +=
            coefficient * std::sqrt(2.0 * static_cast<double>(plus * (zero + 1))) / denominator;
      }
      if (zero > 0) {
        lowered[{plus, zero - 1, minus + 1}] +=
            coefficient * std::sqrt(2.0 * static_cast<double>(zero * (minus + 1))) / denominator;
      }
    }
    state = std::move(lowered);
  }
  return state;
}

// Compute one ordered Cartesian component of a normalized symmetric tensor
double OrderedPolarizationComponent(int rank, int projection, const std::array<int, 3>& occupation) {
  const auto state = PolarizationOccupations(rank, projection);
  const auto found = state.find(occupation);
  if (found == state.end()) { return 0.0; }
  // Use reentrant Gamma evaluation because tensor normalizations enter concurrent amplitude calculations
  const double log_ordered_weight =
      0.5 *
      (math::LogGamma(static_cast<double>(occupation[0] + 1)) + math::LogGamma(static_cast<double>(occupation[1] + 1)) +
       math::LogGamma(static_cast<double>(occupation[2] + 1)) - math::LogGamma(static_cast<double>(rank + 1)));
  return found->second * std::exp(log_ordered_weight);
}

// Compute the overlap of z_<i1...iL> with the normalized m=0 STF tensor
double OrbitalSTFProjection(std::size_t l) {
  const double ld = static_cast<double>(l);
  return std::exp(0.5 * (ld * std::log(2.0) + 2.0 * math::LogGamma(ld + 1.0) - math::LogGamma(2.0 * ld + 1.0)));
}

// Compute physical pole helicities for photons and all projections otherwise
std::vector<double> VertexHelicities(const MParticle& particle) {
  if (particle.pdg == PDG::PDG_gamma) {
    if (particle.spinX2 != 2) { throw std::invalid_argument("PreparePoleLS: a photon must have spin one"); }
    return {-1.0, 1.0};
  }
  return SpinProjections(0.5 * static_cast<double>(particle.spinX2));
}

// Compute one valid unnormalized Jacob-Wick LS basis matrix
HelAmp ReducedLSBasis(int j1x2, int j2x2, int Jx2, std::size_t l, int two_s, bool photon1, bool photon2) {
  const SpinRep     rep1(j1x2, "ReducedLSBasis leg 1");
  const SpinRep     rep2(j2x2, "ReducedLSBasis leg 2");
  const SpinBasis   basis1 = photon1 ? SpinBasis(rep1, {-2, 2}, "ReducedLSBasis photon leg 1") : SpinBasis(rep1);
  const SpinBasis   basis2 = photon2 ? SpinBasis(rep2, {-2, 2}, "ReducedLSBasis photon leg 2") : SpinBasis(rep2);
  const TwoBodySpin spin(SpinRep(Jx2, "ReducedLSBasis mother"), basis1, basis2);
  return spin.JacobWick(l, two_s, "ReducedLSBasis");
}

// Validate angular momentum, discrete symmetries and identical-boson exchange
void ValidateTerm(const MParticle& mother, const MParticle& leg1, const MParticle& leg2, const LSTerm& term, int j1x2,
                  int j2x2, int Jx2, bool production_mode, bool C_symmetry, bool P_symmetry, VertexContext context) {
  const TwoBodySpin spin(SpinRep(Jx2, "ValidateTerm mother"), SpinBasis(SpinRep(j1x2, "ValidateTerm leg 1")),
                         SpinBasis(SpinRep(j2x2, "ValidateTerm leg 2")));
  if (!spin.AllowsLS(term.l, static_cast<int>(term.two_s))) {
    throw std::invalid_argument("PreparePoleLS: LS row violates angular momentum coupling");
  }
  const int orbital_parity = (term.l % 2 == 0) ? 1 : -1;
  if (P_symmetry && mother.P != leg1.P * leg2.P * orbital_parity) {
    throw std::invalid_argument("PreparePoleLS: LS row violates parity conservation");
  }
  const double s = 0.5 * static_cast<double>(term.two_s);
  // Crossed vertices do not describe a physical two-body C eigenstate
  const bool crossed =
      production_mode && (context == VertexContext::CrossedBeamLeg || context == VertexContext::SubTUChannelExchange);
  if (!TwoBodyCParityAllowed(mother, leg1, leg2, term.l, s, C_symmetry && !crossed)) {
    throw std::invalid_argument("PreparePoleLS: LS row violates C parity conservation");
  }
  if (PhysicalIdenticalPair(leg1, leg2, production_mode, context) &&
      ((leg1.spinX2 % 2 == 0 && !BoseSymmetry(static_cast<int>(term.l), static_cast<int>(s))) ||
       (leg1.spinX2 % 2 != 0 && !FermiSymmetry(static_cast<int>(term.l), static_cast<int>(s))))) {
    throw std::invalid_argument(
        "PreparePoleLS: LS row violates "
        "identical-particle symmetry");
  }
}

}  // namespace

// Compute the raw Cartesian STF coupling relative to normalized SU(2) CG
// coefficients
//
// The canonical raw map contracts fixed index pairs with delta_ij, uses
// unit-weight symmetrization followed by STF projection, and does not sum over
// equivalent index choices. When a+b-c is odd, the remaining pair is coupled
// with epsilon_ijk. Each irreducible map is multiplied by one constant phase
// so its highest nonzero spherical component has the real Condon-Shortley CG
// sign. Thus an epsilon map includes the phase that removes the Cartesian i
// from e_+ cross e_0 = i e_+. Card phases refer to this rephased operator
double RawSTFCouplingNormalization(std::size_t rank1, std::size_t rank2, std::size_t output_rank) {
  if (!Triangle(rank1, rank2, output_rank)) {
    throw std::invalid_argument("RawSTFCouplingNormalization: ranks violate the triangle rule");
  }
  const std::size_t difference   = rank1 + rank2 - output_rank;
  const std::size_t contractions = difference / 2;
  const bool        epsilon      = difference % 2 != 0;
  if (rank2 < contractions + static_cast<std::size_t>(epsilon)) {
    throw std::logic_error("RawSTFCouplingNormalization: invalid contraction count");
  }

  const int                rank1_i       = SignedRank(rank1, "RawSTFCouplingNormalization");
  const int                rank2_i       = SignedRank(rank2, "RawSTFCouplingNormalization");
  const int                output_rank_i = SignedRank(output_rank, "RawSTFCouplingNormalization");
  const int                projection    = output_rank_i - rank1_i;
  const std::array<int, 3> occupation    = {static_cast<int>(rank2 - contractions - static_cast<std::size_t>(epsilon)),
                                            epsilon ? 1 : 0, static_cast<int>(contractions)};
  const double             cartesian     = OrderedPolarizationComponent(rank2_i, projection, occupation);
  const double             cg =
      wigner::CG(static_cast<double>(rank1), static_cast<double>(rank2), static_cast<double>(rank1),
                 static_cast<double>(projection), static_cast<double>(output_rank), static_cast<double>(output_rank));
  if (std::abs(cartesian) <= kTolerance || std::abs(cg) <= kTolerance) {
    throw std::logic_error("RawSTFCouplingNormalization: zero highest-weight projection");
  }
  const double normalization = std::abs(cartesian / cg);
  if (!std::isfinite(normalization) || normalization <= 0.0) {
    throw std::overflow_error("RawSTFCouplingNormalization: non-finite normalization");
  }
  return normalization;
}

// Compute the complete raw STF normalization relative to the JW LS matrix
double RawLSOperatorNormalization(std::size_t j1, std::size_t j2, std::size_t J, std::size_t l, std::size_t s) {
  if (!Triangle(j1, j2, s) || !Triangle(l, s, J)) {
    throw std::invalid_argument("RawLSOperatorNormalization: ranks violate an LS triangle rule");
  }
  const double jacob_wick_conversion =
      std::sqrt((2.0 * static_cast<double>(J) + 1.0) / (2.0 * static_cast<double>(l) + 1.0));
  const double normalization = RawSTFCouplingNormalization(j1, j2, s) * RawSTFCouplingNormalization(l, s, J) *
                               OrbitalSTFProjection(l) * jacob_wick_conversion;
  if (!std::isfinite(normalization) || normalization <= 0.0) {
    throw std::overflow_error("RawLSOperatorNormalization: non-finite normalization");
  }
  return normalization;
}

// Compute the raw pole-operator normalization for boson or spinor leg pairs
double RawPoleLSNormalization(int j1x2, int j2x2, int Jx2, std::size_t l, int two_s) {
  if (Jx2 < 0 || Jx2 % 2 != 0 || two_s < 0 || two_s % 2 != 0) {
    throw std::invalid_argument("RawPoleLSNormalization: integer mother and coupled spin required");
  }
  const std::size_t J                 = static_cast<std::size_t>(Jx2 / 2);
  const std::size_t s                 = static_cast<std::size_t>(two_s / 2);
  double            leg_normalization = 0.0;
  if (j1x2 % 2 == 0 && j2x2 % 2 == 0) {
    leg_normalization =
        RawSTFCouplingNormalization(static_cast<std::size_t>(j1x2 / 2), static_cast<std::size_t>(j2x2 / 2), s);
  } else if (j1x2 == 1 && j2x2 == 1) {
    // Canonical epsilon and Pauli spinor bilinears are sqrt(2) times CG rows
    leg_normalization = std::sqrt(2.0);
  } else {
    throw std::invalid_argument("RawPoleLSNormalization: unsupported mixed spinor-tensor legs");
  }
  const double jacob_wick_conversion =
      std::sqrt((2.0 * static_cast<double>(J) + 1.0) / (2.0 * static_cast<double>(l) + 1.0));
  const double normalization =
      leg_normalization * RawSTFCouplingNormalization(l, s, J) * OrbitalSTFProjection(l) * jacob_wick_conversion;
  if (!std::isfinite(normalization) || normalization <= 0.0) {
    throw std::overflow_error("RawPoleLSNormalization: non-finite normalization");
  }
  return normalization;
}

// Enumerate every independent canonical physical-pole operator
std::vector<CanonicalPoleOperator> CanonicalPoleOperators(const MParticle& mother, const MParticle& leg1,
                                                          const MParticle& leg2, bool production_mode, bool C_symmetry,
                                                          bool P_symmetry, VertexContext context) {
  HELMatrix metadata;
  metadata.BR         = 1.0;
  metadata.C_symmetry = C_symmetry;
  metadata.P_symmetry = P_symmetry;
  metadata.alpha_ls.Clear();
  const auto allowed =
      AllowedLSCouplings(metadata, mother, leg1, leg2, production_mode, "CanonicalPoleOperators", context);
  std::vector<CanonicalPoleOperator> output;
  output.reserve(allowed.size());
  for (const auto& row : allowed) {
    const auto pole = ReducedLSBasis(leg1.spinX2, leg2.spinX2, mother.spinX2, row.l, static_cast<int>(row.two_s),
                                     leg1.pdg == PDG::PDG_gamma, leg2.pdg == PDG::PDG_gamma);
    if (pole.FrobNorm2() <= kTolerance * kTolerance) { continue; }
    const double exponent =
        static_cast<double>(row.l) + 0.5 * leg1.spinX2 + 0.5 * leg2.spinX2 - 0.5 * static_cast<double>(row.two_s);
    if (std::abs(exponent - std::round(exponent)) > kTolerance) {
      throw std::invalid_argument("CanonicalPoleOperators: non-integral leg-exchange phase");
    }
    const double phase = std::llround(exponent) % 2 == 0 ? 1.0 : -1.0;
    output.push_back(
        {row, RawPoleLSNormalization(leg1.spinX2, leg2.spinX2, mother.spinX2, row.l, static_cast<int>(row.two_s)),
         phase});
  }
  if (output.empty()) { throw std::invalid_argument("CanonicalPoleOperators: no physical-pole operator is allowed"); }
  return output;
}

// Validate a complete canonical pole coefficient table
void ValidateCanonicalPoleTerms(const MParticle& mother, const MParticle& leg1, const MParticle& leg2,
                                const std::vector<LSTerm>& terms, bool production_mode, bool C_symmetry,
                                bool P_symmetry, VertexContext context) {
  const auto allowed = CanonicalPoleOperators(mother, leg1, leg2, production_mode, C_symmetry, P_symmetry, context);
  std::map<std::pair<std::size_t, std::size_t>, std::complex<double>> supplied;
  for (const auto& term : terms) {
    if (!supplied.emplace(std::make_pair(term.l, term.two_s), term.coefficient).second) {
      throw std::invalid_argument("ValidateCanonicalPoleTerms: duplicate LS operator");
    }
  }
  if (supplied.size() != allowed.size()) {
    throw std::invalid_argument("ValidateCanonicalPoleTerms: explicit table is not complete");
  }
  for (const auto& entry : allowed) {
    if (!supplied.contains({entry.coupling.l, entry.coupling.two_s})) {
      throw std::invalid_argument(
          "ValidateCanonicalPoleTerms: explicit table omits an allowed "
          "operator");
    }
  }
}

// Apply the canonical raw pole normalization to LS coefficients
void ApplyCanonicalPoleLSCoefficients(HELMatrix& helicity, const MParticle& mother, const MParticle& leg1,
                                      const MParticle& leg2) {
  if (helicity.alpha_ls.Empty()) {
    throw std::invalid_argument("ApplyCanonicalPoleLSCoefficients: invalid LS metadata");
  }
  for (LSTerm& term : helicity.alpha_ls) {
    term.coefficient *=
        RawPoleLSNormalization(leg1.spinX2, leg2.spinX2, mother.spinX2, term.l, static_cast<int>(term.two_s));
  }
}

// Prepare immutable pole-normalized STF LS matrices and helicity metadata
PoleLS PreparePoleLS(const MParticle& mother, const MParticle& leg1, const MParticle& leg2,
                     const std::vector<LSTerm>& terms, double Lambda, bool production_mode, bool C_symmetry,
                     bool P_symmetry, VertexContext context, double coupling_min, bool derivative_factor) {
  if (terms.empty()) { throw std::invalid_argument("PreparePoleLS: at least one LS term is required"); }
  if (!std::isfinite(Lambda) || Lambda <= 0.0) {
    throw std::invalid_argument("PreparePoleLS: Lambda must be finite and positive");
  }
  if (!std::isfinite(coupling_min) || coupling_min < 0.0) {
    throw std::invalid_argument("PreparePoleLS: coupling threshold must be nonnegative");
  }
  const std::size_t J    = IntegerSpin(mother, "PreparePoleLS mother");
  const int         Jx2  = static_cast<int>(2 * J);
  const int         j1x2 = ParticleSpinX2(leg1, "PreparePoleLS leg 1");
  const int         j2x2 = ParticleSpinX2(leg2, "PreparePoleLS leg 2");

  PoleLS vertex;
  vertex.rows                 = static_cast<std::size_t>(j1x2 + 1);
  vertex.cols                 = static_cast<std::size_t>(j2x2 + 1);
  vertex.Lambda               = Lambda;
  vertex.context              = context;
  vertex.derivative_factor    = derivative_factor && context != VertexContext::SubTUChannelExchange;
  vertex.leg1_transverse_pole = leg1.pdg == PDG::PDG_gamma;
  vertex.leg2_transverse_pole = leg2.pdg == PDG::PDG_gamma;
  vertex.helicity.C_symmetry  = C_symmetry;
  vertex.helicity.P_symmetry  = P_symmetry;
  InitTwoBodyBasis(vertex.helicity, 0.5 * static_cast<double>(Jx2), 0.5 * static_cast<double>(j1x2),
                   0.5 * static_cast<double>(j2x2), VertexHelicities(mother), VertexHelicities(leg1),
                   VertexHelicities(leg2), "PreparePoleLS");
  vertex.helicity.exchange_basis = context == VertexContext::SubTUChannelExchange && mother.pdg != PDG::PDG_gamma
                                       ? ExchangeBasisType::ReducedRegge
                                       : ExchangeBasisType::HelicityTransport;
  vertex.dual_spin.reserve(vertex.helicity.Jz_values.size());
  vertex.spin_metric.reserve(vertex.helicity.Jz_values.size());
  for (const double Jz : vertex.helicity.Jz_values) {
    if (!std::isfinite(Jz)) { throw std::logic_error("PreparePoleLS: non-finite resonance spin"); }
    vertex.dual_spin.push_back(-Jz);
    vertex.spin_metric.push_back(JacobWickSecondLegReversalPhase(vertex.helicity.J, Jz));
  }

  vertex.rotation = std::make_shared<const wigner::Rotation>(vertex.helicity.jw_rotation.difference,
                                                             vertex.dual_spin, vertex.helicity.J);

  std::set<std::pair<std::size_t, std::size_t>> seen;
  bool                                          nonzero = false;
  for (const auto& term : terms) {
    if (!std::isfinite(term.coefficient.real()) || !std::isfinite(term.coefficient.imag()) ||
        term.l > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2) ||
        term.two_s > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("PreparePoleLS: invalid LS term");
    }
    if (!seen.emplace(term.l, term.two_s).second) {
      throw std::invalid_argument("PreparePoleLS: duplicate LS operator");
    }
    ValidateTerm(mother, leg1, leg2, term, j1x2, j2x2, Jx2, production_mode, C_symmetry, P_symmetry, context);
    if (production_mode && std::abs(term.coefficient) <= coupling_min) { continue; }
    const double raw_normalization = RawPoleLSNormalization(j1x2, j2x2, Jx2, term.l, static_cast<int>(term.two_s));
    auto basis = ReducedLSBasis(j1x2, j2x2, Jx2, term.l, static_cast<int>(term.two_s), vertex.leg1_transverse_pole,
                                vertex.leg2_transverse_pole);
    if (!basis.IsFinite() || basis.FrobNorm2() <= kTolerance * kTolerance) {
      throw std::invalid_argument(
          "PreparePoleLS: LS operator has no physical pole "
          "helicity support");
    }
    if (vertex.reduced_basis.empty()) {
      vertex.rows = basis.size_row();
      vertex.cols = basis.size_col();
    } else if (basis.size_row() != vertex.rows || basis.size_col() != vertex.cols) {
      throw std::logic_error("PreparePoleLS: inconsistent LS operator dimensions");
    }
    basis *= raw_normalization;
    if (!basis.IsFinite()) { throw std::overflow_error("PreparePoleLS: non-finite normalized LS operator"); }
    vertex.terms.push_back(term);
    vertex.raw_normalization.push_back(raw_normalization);
    vertex.reduced_basis.push_back(std::move(basis));
    nonzero = nonzero || std::norm(term.coefficient) > 0.0;
  }
  if (!nonzero && !production_mode) { throw std::invalid_argument("PreparePoleLS: no active LS coefficient"); }
  vertex.ready = true;
  return vertex;
}

// Average the squared coherent pole amplitudes over the lowest active helicity sector
double LeadingPoleDensity(const PoleLS& vertex, double momentum) {
  const auto        reduced = PoleLSReduced(vertex, momentum);
  const std::size_t rows    = vertex.helicity.lambda_values.size_row();
  HelAmp            reference(rows, 1, 0.0);
  for (std::size_t row = 0; row < rows; ++row) {
    const std::size_t i1 = vertex.helicity.lambda_idx[row][0];
    const std::size_t i2 = vertex.helicity.lambda_idx[row][1];
    if (i1 >= reduced.size_row() || i2 >= reduced.size_col()) {
      throw std::invalid_argument(
          "LeadingPoleDensity: helicity index outside the reduced "
          "matrix");
    }
    reference[row][0] = reduced[i1][i2];
  }
  return ForwardHelicityDensity(reference, vertex.helicity.lambda_values,
                                vertex.leg1_transverse_pole ? ForwardLegType::RealPhoton : ForwardLegType::Hadron,
                                vertex.leg2_transverse_pole ? ForwardLegType::RealPhoton : ForwardLegType::Hadron);
}

// Sum canonical pole operators with the momentum rule fixed by vertex role
HelAmp PoleLSReduced(const PoleLS& vertex, double momentum) {
  if (!std::isfinite(momentum) || momentum < 0.0) {
    throw AmplitudeFailure("PoleLSReduced: invalid generated relative momentum");
  }
  const double scaled_momentum = vertex.derivative_factor ? momentum / vertex.Lambda : 1.0;
  if (!std::isfinite(scaled_momentum)) {
    throw AmplitudeFailure("PoleLSReduced: generated derivative factor is non-finite");
  }
  HelAmp reduced(vertex.rows, vertex.cols, 0.0);
  for (std::size_t i = 0; i < vertex.terms.size(); ++i) {
    const auto& term  = vertex.terms[i];
    const auto& basis = vertex.reduced_basis[i];
    // A timelike fusion operator may contain its derivative momentum power
    // A crossed continuum vertex has no two-body decay threshold, so its DL
    // reduced vertex amplitude does not acquire this frame-dependent momentum power
    double derivative = 1.0;
    if (vertex.derivative_factor && term.l > 0) {
      derivative = math::IntegerPower(scaled_momentum, static_cast<unsigned int>(term.l));
    }
    const std::complex<double> scale = term.coefficient * derivative;
    reduced.AddScaled(basis, scale);
  }
  return reduced;
}

// Evaluate the raw LS operators into one reduced-helicity matrix cache
HELMatrix PoleLSHelicity(const PoleLS& vertex, double momentum) {
  HELMatrix helicity      = vertex.helicity;
  helicity.T              = PoleLSReduced(vertex, momentum);
  helicity.T_set          = helicity.T.Transform([](const auto& value) { return std::abs(value) > kTolerance; });
  helicity.coupling_basis = CouplingBasis::Helicity;
  return helicity;
}

// Evaluate one fixed pole as a dimensionless crossed Regge vertex
EvaluatedPoleSubvertex EvaluateCrossedPole(const PoleLS& vertex) {
  auto helicity = PoleLSHelicity(vertex, vertex.Lambda);

  HelAmp        residue(helicity.lambda_values.size_row(), 1, 0.0);
  MMatrix<bool> active(helicity.lambda_values.size_row(), 1, false);
  for (std::size_t row = 0; row < helicity.lambda_values.size_row(); ++row) {
    const std::size_t i1 = helicity.lambda_idx[row][0];
    const std::size_t i2 = helicity.lambda_idx[row][1];

    residue[row][0] = helicity.T[i1][i2];
    active[row][0]  = helicity.T_set[i1][i2];
  }
  helicity.T     = residue;
  helicity.T_set = std::move(active);
  helicity.T_active.clear();
  helicity.Jz_values   = {0.0};
  helicity.jw_rotation = {};
  const auto frame     = ReducedCrossedFrame(helicity);
  return {std::move(helicity), frame};
}

// Evaluate one canonical fixed-spin pole subvertex in its local frame
EvaluatedPoleSubvertex EvaluatePoleSubvertex(const PoleLS& vertex, const M4Vec& final_in_X,
                                             const M4Vec& parent_dir_in_X, bool second_exchange_daughter) {
  const double parent_sign = second_exchange_daughter ? 0.5 : -0.5;
  const M4Vec  relative    = final_in_X + parent_dir_in_X * parent_sign;
  auto         helicity    = PoleLSHelicity(vertex, relative.P3mod());
  return EvaluatePoleSubvertex(std::move(helicity), final_in_X, parent_dir_in_X, second_exchange_daughter);
}

// Prepare only quantities independent of every screening momentum shift
PoleResidue::PoleResidue(const PoleLS& vertex) {
  if (!vertex.ready || vertex.context != VertexContext::SubTUChannelExchange || vertex.derivative_factor) {
    throw std::invalid_argument("PoleResidue requires a prepared continuum vertex without derivative momentum powers");
  }
  const double density = LeadingPoleDensity(vertex, vertex.Lambda);
  if (!std::isfinite(density) || density < 0.0) {
    throw std::invalid_argument("PoleResidue has a non-finite or negative leading pole density");
  }
  // Photon angles are evaluated at the current momentum, so only its helicity tensor is fixed
  auto reduced = vertex.helicity.exchange_basis == ExchangeBasisType::ReducedRegge
                     ? EvaluateCrossedPole(vertex)
                     : EvaluatedPoleSubvertex{PoleLSHelicity(vertex, vertex.Lambda), {}};
  // Copy the input so retained references to its coefficients cannot modify the prepared vertex
  data = std::make_shared<const Data>(Data{vertex, std::move(reduced), density});
}

}  // namespace gra::spin
