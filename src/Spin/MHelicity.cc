// Fundamental spin bases, couplings and density matrix algebra
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// In general, functions here use helicities for fermions which are
// not normalized by the spin vector norm
// That is, functions use [-1/2, 1/2] not for example [-1, 1]
//
// Notation:
// - lambda = lambda1 - lambda2
// - The transition matrix is implemented as
//     f_{lambda1,lambda2;Jz}(theta,phi) = D^{J*}_{Jz,lambda}(phi,theta,0)
//     T_{lambda1,lambda2}
//   with wigner::D(theta,phi,lambda,Jz,J) = d^J_{Jz,lambda}(theta) exp(+i Jz
//   phi)
// - This differs from references that write D^J_{lambda,Jz}; the convention is
// internally consistent

// C++ standard
#include <algorithm>
#include <cmath>
#include <complex>
#include <iostream>
#include <sstream>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "rang.hpp"

using gra::aux::indices;

using gra::math::msqrt;
using gra::math::zi;

namespace gra {
namespace spin {

// Used for checking eigenvalues
const double epsilon = 1e-6;

constexpr double kSpinLabelTolerance = 1e-9;

// Infer a daughter spin from explicit helicity metadata when the cache lacks s1/s2
// The maximum listed |lambda| fixes the helicity basis spin
//
double SpinFromHelicityColumn(const gra::HELMatrix &hel, std::size_t column) {
  if (column >= hel.lambda_values.size_col()) {
    throw std::invalid_argument("MHelicity: Requested helicity column outside metadata bounds");
  }

  double spin = 0.0;
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) {
    spin = std::max(spin, std::abs(hel.lambda_values[i][column]));
  }
  return SpinLabelX2(spin, "MHelicity::SpinFromHelicityColumn") / 2.0;
}

// Find a row with the requested two-helicity assignment
// Crossed-row remapping needs the exact lambda1,lambda2 partner
//
std::size_t FindHelicityPairRow(const gra::HELMatrix &hel, double lambda1, double lambda2) {
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) {
    if (std::abs(hel.lambda_values[i][0] - lambda1) < kSpinLabelTolerance &&
        std::abs(hel.lambda_values[i][1] - lambda2) < kSpinLabelTolerance) {
      return i;
    }
  }

  throw std::invalid_argument("MHelicity: Could not find crossed helicity partner row");
}

// Convert auxiliary crossed-antiparticle rows into physical second-leg helicity rows
// Production cards use crossed auxiliary rows but amplitudes use physical rows
//
MMatrix<std::complex<double>> CrossedSecondRowsToPhysical(const gra::HELMatrix                &hel,
                                                          const MMatrix<std::complex<double>> &aux) {
  const double             s2 = (hel.s2 >= 0.0) ? hel.s2 : SpinFromHelicityColumn(hel, 1);
  std::vector<std::size_t> permutation(aux.size_row());
  std::vector<double>      phases(aux.size_row());
  for (std::size_t row = 0; row < hel.lambda_values.size_row(); ++row) {
    const double lambda1          = hel.lambda_values[row][0];
    const double physical_lambda2 = hel.lambda_values[row][1];
    permutation[row]              = FindHelicityPairRow(hel, lambda1, -physical_lambda2);
    phases[row]                   = JacobWickSecondLegReversalPhase(s2, physical_lambda2);
  }
  return aux.SelectRows(permutation).LeftDiagonalProduct(phases);
}

// Rotate helicity columns between azimuthal spin sections
// A spin projection m picks up exp(i c m phi) under the rotation
//
void ApplyColumnHelicityPhase(MMatrix<std::complex<double>> &f, const gra::HELMatrix &hel, double phi,
                              double coefficient) {
  for (std::size_t col = 0; col < f.size_col(); ++col) {
    f.ScaleColumn(col, std::exp(gra::math::zi * coefficient * hel.Jz_values[col] * phi));
  }
}

// Build one reduced crossed Regge frame without an angular rotation
MMatrix<std::complex<double>> ReducedCrossedFrame(const gra::HELMatrix &hel) {
  return CrossedSecondRowsToPhysical(hel, gra::spin::DirectFrame(hel, "MHelicity::ReducedCrossedFrame"));
}

// Build one angle-free crossed Regge-helicity frame with its spin section
MMatrix<std::complex<double>> ReggeCrossedFrame(const gra::HELMatrix &hel, const gra::M4Vec &axis,
                                                const bool second_exchange_daughter) {
  MMatrix<std::complex<double>> out =
      CrossedSecondRowsToPhysical(hel, gra::spin::DirectFrame(hel, "MHelicity::ReggeCrossedFrame"));
  ApplyColumnHelicityPhase(out, hel, axis.Phi(), second_exchange_daughter ? 1.0 : -1.0);
  // Convert the Fourier vertex to the dual spherical basis of the forward source
  // The lower axis is reversed relative to its collider transfer momentum
  for (const auto &column : indices(hel.Jz_values)) {
    const int m = static_cast<int>(std::llround(hel.Jz_values[column]));
    out.ScaleColumn(column, SphericalSection(second_exchange_daughter ? -m : m));
  }
  return out;
}

// Build one virtual subchannel vertex in the timelike central-system frame
// An off-shell t/u pair has no real Jacob-Wick rest frame
MMatrix<std::complex<double>> VirtualSubchannelFrame(const gra::HELMatrix &hel, gra::M4Vec final_in_X,
                                                     const gra::M4Vec &parent_dir_in_X, bool second_exchange_daughter) {
  if (hel.exchange_basis == ExchangeBasisType::ReducedRegge) { return ReducedCrossedFrame(hel); }
  if (hel.UsesReggeDomain()) {
    throw std::invalid_argument(
        "MHelicity::VirtualSubchannelFrame requires a reduced Regge "
        "residue or a physical pole tensor");
  }

  gra::kinematics::RotateHelicityAxes(final_in_X, parent_dir_in_X);
  MMatrix<std::complex<double>> out =
      CrossedSecondRowsToPhysical(hel, gra::spin::fDecayMatrix(hel, final_in_X.Theta(), final_in_X.Phi()));
  // Physical pole columns use one ket section before the lower metric is applied
  ApplyColumnHelicityPhase(out, hel, parent_dir_in_X.Phi(), -1.0);
  return out;
}

// Evaluate one prepared physical-pole helicity tensor in its local frame
EvaluatedPoleSubvertex EvaluatePoleSubvertex(HELMatrix helicity, const gra::M4Vec &final_in_X,
                                             const gra::M4Vec &parent_dir_in_X, bool second_exchange_daughter) {
  const double     parent_sign = second_exchange_daughter ? 0.5 : -0.5;
  const gra::M4Vec relative    = final_in_X + parent_dir_in_X * parent_sign;
  auto             frame       = VirtualSubchannelFrame(helicity, relative, parent_dir_in_X, second_exchange_daughter);
  if (second_exchange_daughter) {
    std::vector<std::size_t> crossed;
    std::vector<double>      metric;
    crossed.reserve(helicity.Jz_values.size());
    metric.reserve(helicity.Jz_values.size());
    for (const double projection : helicity.Jz_values) {
      const auto found =
          std::find_if(helicity.Jz_values.cbegin(), helicity.Jz_values.cend(),
                       [projection](const double value) { return std::abs(value + projection) < kSpinLabelTolerance; });
      if (found == helicity.Jz_values.cend()) {
        throw std::invalid_argument(
            "MHelicity::EvaluatePoleSubvertex requires a complete "
            "physical pole basis");
      }
      crossed.push_back(static_cast<std::size_t>(found - helicity.Jz_values.cbegin()));
      metric.push_back(JacobWickSecondLegReversalPhase(helicity.J, projection));
    }
    frame = frame.SelectColumns(crossed).RightDiagonalProduct(metric);
  }
  return {std::move(helicity), std::move(frame)};
}

// Compute the physical spin-half helicity transitions
std::vector<std::pair<double, double>> SpinHalfTransitions(bool drop_flip) {
  std::vector<std::pair<double, double>> transitions;
  const auto                             helicities = BinaryHelicityLabelsX2();
  for (const int lambda_in_x2 : helicities) {
    for (const int lambda_out_x2 : helicities) {
      const double lambda_in  = lambda_in_x2 / 2.0;
      const double lambda_out = lambda_out_x2 / 2.0;
      if (drop_flip && std::abs(lambda_in - lambda_out) > kSpinLabelTolerance) { continue; }
      transitions.push_back({lambda_in, lambda_out});
    }
  }
  return transitions;
}

// Compute one spherical forward helicity barrier and azimuthal factor
std::complex<double> HelicityFactor(int m, double delta_lambda, double qt, double phi, double exchange_mass,
                                    double flip_mass, bool reverse_phase, bool use_exchange_barrier) {
  if (!(exchange_mass > 0.0) || !(flip_mass > 0.0)) {
    throw std::invalid_argument("gra::spin::HelicityFactor: mass scales must be positive");
  }

  const double exchange_barrier =
      use_exchange_barrier ? math::IntegerPower(qt / exchange_mass, static_cast<unsigned int>(std::abs(m))) : 1.0;
  const double flip_barrier = std::pow(qt / flip_mass, std::abs(delta_lambda));
  const double phase_sign   = reverse_phase ? -1.0 : 1.0;
  const double angle        = phase_sign * static_cast<double>(m) * phi;
  // Complete the Condon-Shortley section of a real spherical source
  // The real spherical source obeys R_{-m} = (-1)^m R_m^* on either ordered exchange leg
  return SphericalSection(m) * exchange_barrier * flip_barrier * std::complex<double>(std::cos(angle), std::sin(angle));
}

// Compute the external helicity phase of one ordered collider beam source
std::complex<double> ForwardHelicitySectionPhase(const double lambda_in, const double lambda_out,
                                                 const double collider_transfer_azimuth,
                                                 const bool   second_exchange_daughter) {
  const double delta_lambda  = lambda_in - lambda_out;
  const double phase_sign = second_exchange_daughter ? -1.0 : 1.0;
  return std::exp(zi * phase_sign * delta_lambda * collider_transfer_azimuth);
}

// Compute the complete external section factor of a factorized spin-half source
std::complex<double> SpinHalfForwardHelicitySectionFactor(const double lambda_in, const double lambda_out,
                                                          const double collider_transfer_azimuth,
                                                          const bool   second_exchange_daughter) {
  double reference_sign = 1.0;
  if (std::abs(lambda_in - lambda_out) > kSpinLabelTolerance) {
    reference_sign = (second_exchange_daughter ? -2.0 : 2.0) * lambda_in;
  }
  return reference_sign *
         ForwardHelicitySectionPhase(lambda_in, lambda_out, collider_transfer_azimuth, second_exchange_daughter);
}

// Compute the Jacob-Wick phase from reversing the ordered second leg
double JacobWickSecondLegReversalPhase(const double spin, const double physical_helicity) {
  static_cast<void>(SpinLabelX2(spin, "gra::spin::JacobWickSecondLegReversalPhase spin"));
  if (!std::isfinite(physical_helicity) || std::abs(physical_helicity) > spin + kSpinLabelTolerance) {
    throw std::invalid_argument("gra::spin::JacobWickSecondLegReversalPhase: invalid helicity");
  }
  const double exponent = spin - physical_helicity;
  const double rounded  = std::round(exponent);
  if (std::abs(exponent - rounded) > kSpinLabelTolerance) {
    throw std::invalid_argument("gra::spin::JacobWickSecondLegReversalPhase: non-integral exponent");
  }
  return (std::llround(rounded) % 2 == 0) ? 1.0 : -1.0;
}

// Compute a validated direct helicity basis matrix
MMatrix<std::complex<double>> DirectFrame(const gra::HELMatrix &hel, const std::string &context) {
  if (hel.UsesReggeDomain()) {
    if (hel.analytic_MMAX < 0 || hel.Jz_values.size() != static_cast<std::size_t>(2 * hel.analytic_MMAX + 1)) {
      throw std::invalid_argument(context + ": direct trajectory basis has invalid dimensions");
    }
  }
  if (hel.T.size_row() != hel.lambda_values.size_row() || hel.T.size_col() != hel.Jz_values.size()) {
    throw std::invalid_argument(context + ": direct helicity matrix has invalid dimensions");
  }
  return hel.T;
}

// Compute true when a final particle has only the two transverse vector helicities
// Real photons and gluons have no physical lambda=0 polarization
bool IsTransverseMasslessVector(const gra::MParticle &particle) {
  const bool massless_vector = particle.pdg == PDG::PDG_gamma || particle.pdg == PDG::PDG_gluon;
  if (massless_vector && particle.spinX2 != 2) {
    throw std::invalid_argument("gra::spin::IsTransverseMasslessVector: photon/gluon spinX2 must be 2");
  }
  return massless_vector;
}

// Compute physical final-state helicities for one external particle
// Massive particles use SU(2), while real gamma/g use lambda=+-1
std::vector<double> FinalStateHelicities(const gra::MParticle &particle, const std::string &context) {
  if (IsTransverseMasslessVector(particle)) {
    const auto labels = BinaryHelicityLabelsX2();
    return {static_cast<double>(labels[0]), static_cast<double>(labels[1])};
  }
  return SpinRep(particle.spinX2, context).Projections();
}

// Compute the physical final-state helicity dimension of one particle
std::size_t FinalStateHelicityCount(const gra::MParticle &particle, const std::string &context) {
  return IsTransverseMasslessVector(particle) ? 2 : SpinRep(particle.spinX2, context).Dim();
}

// Compute the compact final-state index of one physical helicity
std::size_t FinalStateHelicityIndex(const gra::MParticle &particle, double helicity, const std::string &context) {
  const auto helicities = FinalStateHelicities(particle, context);
  for (const auto &i : indices(helicities)) {
    if (std::abs(helicities[i] - helicity) < kSpinLabelTolerance) { return i; }
  }
  throw std::invalid_argument(context + ": helicity is not in the physical final-state basis");
}

// Fill pair labels and full SU(2) matrix indices from explicit helicity lists
// Compact physical rows may point to non-contiguous entries of the
// spin matrix
void FillHelicityPairProjections(MMatrix<double> &lambda_values, MMatrix<std::size_t> &lambda_idx,
                                 const std::vector<double> &helicities1, const std::vector<double> &helicities2,
                                 double s1, double s2, const std::string &context) {
  const std::size_t count = helicities1.size() * helicities2.size();
  lambda_values           = MMatrix<double>(count, 2, 0.0);
  lambda_idx              = MMatrix<std::size_t>(count, 2, 0);

  std::size_t row = 0;
  for (const double lambda1 : helicities1) {
    for (const double lambda2 : helicities2) {
      lambda_values[row][0] = lambda1;
      lambda_values[row][1] = lambda2;
      lambda_idx[row][0]    = SpinProjectionIndex(lambda1, s1, context + " leg 1");
      lambda_idx[row][1]    = SpinProjectionIndex(lambda2, s2, context + " leg 2");
      ++row;
    }
  }
}

// Initialize one two-body helicity basis and its Jacob-Wick row cache
void InitTwoBodyBasis(gra::HELMatrix &hel, double J, double s1, double s2, const std::vector<double> &parent_helicities,
                      const std::vector<double> &helicities1, const std::vector<double> &helicities2,
                      const std::string &context) {
  hel.J         = J;
  hel.s1        = s1;
  hel.s2        = s2;
  hel.Jz_values = parent_helicities;
  FillHelicityPairProjections(hel.lambda_values, hel.lambda_idx, helicities1, helicities2, s1, s2, context);
  InitJWRotation(hel);
}

// Compute whether one daughter is an external final particle at this JW vertex
// Decay daughters and the first sub-t/u daughter are physical outgoing
// states
bool IsFinalVertexLeg(bool production_mode, VertexContext context, std::size_t leg) {
  if (!production_mode) { return true; }
  return context == VertexContext::SubTUChannelExchange && leg == 0;
}

// Compute whether one vertex leg is a physical external massless-vector state
// Central production photons are transverse source states
bool IsPhysicalMasslessVectorLeg(const gra::MParticle &particle, bool production_mode, VertexContext context,
                                 std::size_t leg) {
  if (IsFinalVertexLeg(production_mode, context, leg)) { return true; }
  return production_mode && context == VertexContext::Auto && particle.pdg == PDG::PDG_gamma;
}

// Compute the helicity list used for one leg of a JW vertex
// Physical massless vectors use only transverse helicities
std::vector<double> VertexLegHelicities(const gra::MParticle &particle, bool production_mode, VertexContext context,
                                        std::size_t leg, const std::string &label) {
  if (IsPhysicalMasslessVectorLeg(particle, production_mode, context, leg) ||
      IsFinalVertexLeg(production_mode, context, leg)) {
    return FinalStateHelicities(particle, label);
  }
  return SpinProjections(particle.spinX2 / 2.0);
}

// Project an LS matrix onto the physical external massless-vector subspace
// Only decay tensors are normalized after the physical projection
// [REFERENCE: Landau, Dokl. Akad. Nauk SSSR 60 (1948) 207]
// [REFERENCE: Yang, Phys. Rev. 77 (1950) 242]
void ProjectFinalMasslessVectorHelicities(MMatrix<std::complex<double>> &T, const gra::MParticle &p1,
                                          const gra::MParticle &p2, bool production_mode, VertexContext context,
                                          const std::string &label) {
  const bool restrict1 = IsPhysicalMasslessVectorLeg(p1, production_mode, context, 0) && IsTransverseMasslessVector(p1);
  const bool restrict2 = IsPhysicalMasslessVectorLeg(p2, production_mode, context, 1) && IsTransverseMasslessVector(p2);
  if (!restrict1 && !restrict2) { return; }

  for (std::size_t i = 0; i < T.size_row(); ++i) {
    const double lambda1 = -0.5 * p1.spinX2 + static_cast<double>(i);
    for (std::size_t j = 0; j < T.size_col(); ++j) {
      const double lambda2 = -0.5 * p2.spinX2 + static_cast<double>(j);
      if ((restrict1 && std::abs(lambda1) < kSpinLabelTolerance) ||
          (restrict2 && std::abs(lambda2) < kSpinLabelTolerance)) {
        T[i][j] = 0.0;
      }
    }
  }

  const double norm2 = T.FrobNorm2();
  if (norm2 < 1e-20) {
    throw std::invalid_argument(label + ": no physical transverse massless-vector helicity amplitude");
  }
  if (!production_mode) { T *= 1.0 / std::sqrt(norm2); }
}

// Compute a compact integer or half-integer spin label for diagnostics
// Labels are printed as helicity quantum numbers
//
std::string FormatSpinLabel(double value) {
  const double rounded = std::round(value);
  if (std::abs(value - rounded) < 1e-12) { return std::to_string(static_cast<int>(std::llround(rounded))); }
  std::ostringstream out;
  out << value;
  return out.str();
}

namespace {

// Identify an LS row rejected by the physical selection rules
//
class LSCouplingRejected final : public std::invalid_argument {
 public:
  using std::invalid_argument::invalid_argument;
};

// Compute integer total spin needed by identical-particle statistics tests
// Bose and Fermi exchange rules use integer LS exponents
//
int IntegerSpin(double spin, const std::string &context) {
  const double rounded = std::round(spin);
  if (std::abs(spin - rounded) > kSpinLabelTolerance) {
    throw std::invalid_argument(context + ": total spin is not integer");
  }
  return static_cast<int>(std::llround(rounded));
}

// Build one unnormalized Jacob-Wick LS basis matrix
// All selection rules are applied once per LS row before recoupling
MMatrix<std::complex<double>> LSHelicityBasis(const gra::HELMatrix &prototype, const gra::MParticle &p,
                                              const gra::MParticle &p1, const gra::MParticle &p2, std::size_t l,
                                              std::size_t two_s, bool production_mode, VertexContext context) {
  const double                  s = 0.5 * static_cast<double>(two_s);
  const TwoBodySpin             spin(SpinRep(p.spinX2, "gra::spin::LSHelicityBasis mother"),
                                     SpinBasis(SpinRep(p1.spinX2, "gra::spin::LSHelicityBasis daughter 1")),
                                     SpinBasis(SpinRep(p2.spinX2, "gra::spin::LSHelicityBasis daughter 2")));
  MMatrix<std::complex<double>> basis(static_cast<std::size_t>(p1.spinX2 + 1), static_cast<std::size_t>(p2.spinX2 + 1),
                                      0.0);

  if (!spin.AllowsLS(l, static_cast<int>(two_s))) { return basis; }
  const bool crossed =
      production_mode && (context == VertexContext::CrossedBeamLeg || context == VertexContext::SubTUChannelExchange);
  if (prototype.P_symmetry && p1.P * p2.P * ((l % 2 == 0) ? 1 : -1) != p.P) { return basis; }
  if (!TwoBodyCParityAllowed(p, p1, p2, l, s, prototype.C_symmetry && !crossed)) { return basis; }

  const bool identical    = PhysicalIdenticalPair(p1, p2, production_mode, context);
  const bool boson_pair   = p1.spinX2 % 2 == 0 && p2.spinX2 % 2 == 0;
  const bool fermion_pair = p1.spinX2 % 2 != 0 && p2.spinX2 % 2 != 0;
  if (identical && boson_pair &&
      !BoseSymmetry(static_cast<int>(l), IntegerSpin(s, "gra::spin::LSHelicityBasis Bose symmetry"))) {
    return basis;
  }
  if (identical && fermion_pair &&
      !FermiSymmetry(static_cast<int>(l), IntegerSpin(s, "gra::spin::LSHelicityBasis Fermi symmetry"))) {
    return basis;
  }

  return spin.JacobWick(l, static_cast<int>(two_s), "gra::spin::LSHelicityBasis");
}

// Build one physical LS basis row after external-helicity projection
MMatrix<std::complex<double>> PhysicalLSHelicityBasis(const gra::HELMatrix &prototype, const gra::MParticle &p,
                                                      const gra::MParticle &p1, const gra::MParticle &p2, std::size_t l,
                                                      std::size_t two_s, bool production_mode, VertexContext context) {
  auto basis = LSHelicityBasis(prototype, p, p1, p2, l, two_s, production_mode, context);
  if (basis.FrobNorm2() <= 1e-24) { return {}; }
  try {
    ProjectFinalMasslessVectorHelicities(basis, p1, p2, production_mode, context, "gra::spin::PhysicalLSHelicityBasis");
  } catch (const std::invalid_argument &) { return {}; }
  return basis;
}

// Visit every allowed physical LS matrix in orbital order
template <typename Visitor>
void ForEachLS(const HELMatrix &hc, const MParticle &p, const MParticle &p1, const MParticle &p2, bool production_mode,
               VertexContext context, Visitor &&visit) {
  const SpinRep     parent(p.spinX2, "LS mother"), first(p1.spinX2, "LS first leg"), second(p2.spinX2, "LS second leg");
  const std::size_t max_two_s = static_cast<std::size_t>(first.X2()) + second.X2();
  const std::size_t max_l     = (parent.X2() + max_two_s) / 2;
  for (std::size_t l = 0; l <= max_l; ++l) {
    for (std::size_t two_s = static_cast<std::size_t>(std::abs(first.X2() - second.X2())); two_s <= max_two_s;
         two_s += 2) {
      auto matrix = PhysicalLSHelicityBasis(hc, p, p1, p2, l, two_s, production_mode, context);
      if (!matrix.isEmpty()) { visit(LSCoupling{l, two_s}, std::move(matrix)); }
    }
  }
}

// Store the allowed direct-helicity subspace generated by all valid LS rows
// Direct H(lambda1,lambda2) must lie in the JW LS span
//
struct DirectHelicitySubspace {
  std::vector<std::vector<std::complex<double>>> basis;
  std::vector<std::vector<std::complex<double>>> raw_basis;
  std::vector<MMatrix<std::complex<double>>>     raw_matrices;
  std::vector<LSCoupling>                        raw_rows;
  MMatrix<bool>                                  support;
};

}  // namespace

namespace {

// Compute a compact LS-row label for diagnostics
//
std::string FormatLSCouplingLabel(const LSCoupling &row) {
  return "(" + std::to_string(row.l) + "," + gra::aux::Spin2XtoString(static_cast<int>(row.two_s)) + ")";
}

// Compute the equivalent LS coefficients of one direct-helicity vector
//
std::vector<std::complex<double>> EquivalentLSCoefficients(const DirectHelicitySubspace            &subspace,
                                                           const std::vector<std::complex<double>> &helicity_vector) {
  if (subspace.raw_basis.size() != subspace.raw_rows.size()) {
    throw std::invalid_argument("gra::spin::EquivalentLSCoefficients: LS row bookkeeping mismatch");
  }
  if (subspace.raw_basis.empty()) { return {}; }

  const auto basis = gra::MMatrix<std::complex<double>>::FromColumns(subspace.raw_basis);
  if (basis.size_row() != helicity_vector.size()) {
    throw std::invalid_argument(
        "gra::spin::EquivalentLSCoefficients: basis "
        "vector dimension mismatch");
  }
  // Resolve only the independent physical helicities and choose the minimum-norm LS coefficients
  const auto projection = gra::MMatrix<std::complex<double>>::FromColumns(subspace.basis).Dagger();
  return (projection * basis).PseudoInverse(0.0) * (projection * helicity_vector);
}

// Print the direct-helicity coefficients as a compact table
//
void PrintDirectHelicityCoefficientTable(const gra::HELMatrix &hc, double s1, double s2, const std::string &title) {
  std::vector<std::vector<std::string>> rows;
  for (std::size_t i = 0; i < hc.T.size_row(); ++i) {
    for (std::size_t j = 0; j < hc.T.size_col(); ++j) {
      if (std::abs(hc.T[i][j]) <= 1e-12) { continue; }
      const double lambda1 = -s1 + static_cast<double>(i);
      const double lambda2 = -s2 + static_cast<double>(j);
      rows.push_back({FormatSpinLabel(lambda1), FormatSpinLabel(lambda2), gra::aux::ToString(std::abs(hc.T[i][j]), 3),
                      gra::aux::ToString(std::arg(hc.T[i][j]), 3)});
    }
  }
  if (rows.empty()) { return; }
  std::cout << title << std::endl;
  gra::aux::PrintTable({"lambda1", "lambda2", "|H|", "arg(H)"}, rows);
  std::cout << std::endl;
}

// Print the LS coefficients equivalent to one direct-helicity vector
//
void PrintEquivalentLSCoefficientTable(const DirectHelicitySubspace            &subspace,
                                       const std::vector<std::complex<double>> &helicity_vector,
                                       const std::string                       &title) {
  const auto                            coefficients = EquivalentLSCoefficients(subspace, helicity_vector);
  std::vector<std::vector<std::string>> rows;
  rows.reserve(coefficients.size());
  for (const auto &i : indices(coefficients)) {
    rows.push_back({FormatLSCouplingLabel(subspace.raw_rows[i]), gra::aux::ToString(std::abs(coefficients[i]), 3),
                    gra::aux::ToString(std::arg(coefficients[i]), 3),
                    std::abs(coefficients[i]) > 1e-12 ? "yes" : "no"});
  }
  if (rows.empty()) { return; }
  std::cout << title << std::endl;
  gra::aux::PrintTable({"(l,s)", "|alpha_eq|", "arg(alpha_eq)", "active"}, rows);
  std::cout << std::endl;
}

// Compute allowed helicity support and basis vectors from the JW LS construction
// Each accepted (l,S) row is recoupled with CG coefficients
//
DirectHelicitySubspace BuildDirectHelicitySubspace(const gra::HELMatrix &hc, const gra::MParticle &p,
                                                   const gra::MParticle &p1, const gra::MParticle &p2,
                                                   bool production_mode, VertexContext context) {
  DirectHelicitySubspace out;
  out.support = MMatrix<bool>(SpinStateCount(p1.spinX2 / 2.0, "LS first leg"),
                              SpinStateCount(p2.spinX2 / 2.0, "LS second leg"), false);
  ForEachLS(hc, p, p1, p2, production_mode, context, [&](LSCoupling row, auto matrix) {
    const auto raw = matrix.Flatten();
    for (std::size_t i = 0; i < matrix.size_row(); ++i) {
      for (std::size_t j = 0; j < matrix.size_col(); ++j) {
        if (std::abs(matrix[i][j]) > 1e-10) { out.support[i][j] = true; }
      }
    }
    out.raw_basis.push_back(raw);
    out.raw_matrices.push_back(std::move(matrix));
    out.raw_rows.push_back(row);
    gra::AddOrthonormal(out.basis, raw, 1e-12);
  });

  return out;
}

// Check whether the two daughters form an ordered particle-antiparticle pair
// Charged two-body C parity is defined for particle-antiparticle pairs
//
bool IsParticleAntiparticlePair(const gra::MParticle &p1, const gra::MParticle &p2) { return p1.pdg == -p2.pdg; }

// Evaluate the C-parity of a two-body LS state when it is well-defined
// Particle-antiparticle states use C = (-1)^(l+S)
//
bool TwoBodyCParity(const gra::MParticle &p1, const gra::MParticle &p2, std::size_t l, double s, int &Ctot) {
  if (IsParticleAntiparticlePair(p1, p2)) {
    const double exponent = static_cast<double>(l) + s;
    if (std::abs(exponent - std::round(exponent)) > 1e-9) { return false; }
    Ctot = (static_cast<long long>(std::llround(exponent)) % 2 == 0) ? 1 : -1;
    return true;
  }

  if (IsValidCParity(p1.C) && IsValidCParity(p2.C)) {
    Ctot = p1.C * p2.C;
    return true;
  }

  return false;
}

}  // namespace

// Check if an integer stores a defined C-parity eigenvalue
// Only neutral C eigenstates carry C = +-1
//
bool IsValidCParity(int C) { return C == -1 || C == 1; }

// Compute whether the two physical legs are identical in the current vertex context
// Crossed production vertices use C-eigen metadata for neutral legs
//
bool PhysicalIdenticalPair(const gra::MParticle &p1, const gra::MParticle &p2, bool production_mode,
                           VertexContext context) {
  const bool crossed =
      production_mode && (context == VertexContext::CrossedBeamLeg || context == VertexContext::SubTUChannelExchange);
  return crossed ? (p1.pdg == p2.pdg && IsValidCParity(p1.C) && IsValidCParity(p2.C)) : (p1.pdg == p2.pdg);
}

// Compute true when one two-body LS row passes the requested C-parity constraint
// C filtering is applied only when the two-body C value is defined
//
bool TwoBodyCParityAllowed(const gra::MParticle &mother, const gra::MParticle &p1, const gra::MParticle &p2,
                           std::size_t l, double s, bool c_conservation) {
  if (!c_conservation || !IsValidCParity(mother.C)) { return true; }

  int Ctot = 0;
  if (!TwoBodyCParity(p1, p2, l, s, Ctot)) { return true; }
  return Ctot == mother.C;
}

// Compute the flattened row-major direct-helicity coordinate index
// This labels one H(lambda1,lambda2) coordinate
//
std::size_t DirectHelicityCoordinateIndex(double lambda1, double lambda2, double s1, double s2,
                                          const std::string &context) {
  const std::size_t n2 = SpinStateCount(s2, context + " leg 2");
  const std::size_t i1 = SpinProjectionIndex(lambda1, s1, context + " leg 1");
  const std::size_t i2 = SpinProjectionIndex(lambda2, s2, context + " leg 2");
  return i1 * n2 + i2;
}

// Compute the direct-helicity leg-exchange phase for swapped table lookup
// Two-body helicity exchange gives (-1)^(s1+s2-J)
//
double DirectHelicityLegExchangePhase(const gra::MParticle &mother, const gra::MParticle &leg1,
                                      const gra::MParticle &leg2, const std::string &context) {
  const double J        = mother.spinX2 / 2.0;
  const double s1       = leg1.spinX2 / 2.0;
  const double s2       = leg2.spinX2 / 2.0;
  const double exponent = s1 + s2 - J;
  const double rounded  = std::round(exponent);
  if (std::abs(exponent - rounded) > kSpinLabelTolerance) {
    throw std::invalid_argument(context + ": non-integral direct-helicity leg-exchange phase");
  }
  return (std::llround(rounded) % 2 == 0) ? 1.0 : -1.0;
}

}  // namespace spin

namespace spin {

// Compute the active LS rows in minimal orbital-angular-momentum order
// Allowed rows satisfy angular momentum, P, C and statistics rules
//
std::vector<LSCoupling> AllowedLSCouplings(const gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1,
                                           const gra::MParticle &p2, bool production_mode,
                                           const std::string &source_name, VertexContext context) {
  (void)source_name;
  std::vector<LSCoupling> out;
  ForEachLS(hc, p, p1, p2, production_mode, context, [&](LSCoupling row, const auto &) { out.push_back(row); });
  return out;
}

// Normalize one non-zero direct helicity vector in place
// Reduced helicity vectors are unit-normalized before storage
//
void NormalizeDirectHelicityVector(std::vector<std::complex<double>> &v) {
  const double norm2 = gra::SquaredNorm(v);
  if (norm2 <= 1e-24) { throw std::invalid_argument("gra::spin::NormalizeDirectHelicityVector: zero vector"); }
  v = gra::NormalizedL2(v);
}

// Fix an arbitrary direct-vector phase by making the first non-zero element real positive
// The overall phase of a single vertex basis vector is unobservable
//
void CanonicalizeDirectHelicityVectorPhase(std::vector<std::complex<double>> &v) {
  for (const auto &x : v) {
    if (std::abs(x) <= 1e-12) { continue; }
    const std::complex<double> phase = std::conj(x) / std::abs(x);
    gra::Scale(v, phase);
    return;
  }
}

// Add one vector to an orthonormal direct-helicity basis
// This accumulates independent JW tensor structures
//
void AddDirectHelicityBasisVector(std::vector<std::vector<std::complex<double>>> &basis,
                                  const std::vector<std::complex<double>> &candidate, double tol) {
  gra::AddOrthonormal(basis, candidate, tol);
}

// Compute the orthonormal JW direct-helicity basis generated by all allowed LS rows
// Direct cards are checked against the same LS tensor algebra
//
DirectHelicityBasis BuildJWDirectHelicityBasis(const gra::HELMatrix &hc, const gra::MParticle &p,
                                               const gra::MParticle &p1, const gra::MParticle &p2, bool production_mode,
                                               const std::string &source_name, const std::string &verbose_label,
                                               VertexContext context) {
  (void)source_name;
  (void)verbose_label;
  const auto local = BuildDirectHelicitySubspace(hc, p, p1, p2, production_mode, context);
  return {local.basis, local.support};
}

// Compute the raw LS-generated direct-helicity vectors before orthonormalization
// Raw vectors preserve coordinate phases for parity-pair inference
//
std::vector<std::vector<std::complex<double>>> BuildJWDirectHelicityRawBasis(
    const gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1, const gra::MParticle &p2,
    bool production_mode, const std::string &source_name, const std::string &verbose_label, VertexContext context) {
  (void)source_name;
  (void)verbose_label;
  const auto local = BuildDirectHelicitySubspace(hc, p, p1, p2, production_mode, context);
  return local.raw_basis;
}

// Convert one direct-helicity matrix into its complete JW LS coefficients
DirectHelicityLSExpansion DirectHelicityToLSCoefficients(const gra::HELMatrix &helicity, const gra::MParticle &mother,
                                                         const gra::MParticle &leg1, const gra::MParticle &leg2,
                                                         bool production_mode, VertexContext context) {
  if (helicity.T.isEmpty() || !helicity.T.IsFinite()) {
    throw std::invalid_argument("gra::spin::DirectHelicityToLSCoefficients: invalid helicity matrix");
  }
  const auto subspace = BuildDirectHelicitySubspace(helicity, mother, leg1, leg2, production_mode, context);
  const auto target   = helicity.T.Flatten();
  DirectHelicityLSExpansion out;
  out.rows         = subspace.raw_rows;
  out.coefficients = EquivalentLSCoefficients(subspace, target);
  if (out.rows.size() != out.coefficients.size()) {
    throw std::invalid_argument("gra::spin::DirectHelicityToLSCoefficients: LS dimensions differ");
  }
  std::vector<std::complex<double>> residual = target;
  for (const auto &i : indices(out.coefficients)) {
    gra::AddScaled(residual, subspace.raw_basis[i], -out.coefficients[i]);
  }
  out.residual_norm2 = gra::SquaredNorm(residual);
  if (!std::isfinite(out.residual_norm2)) {
    throw std::invalid_argument("gra::spin::DirectHelicityToLSCoefficients: non-finite residual");
  }
  return out;
}

// Compute the fixed phase relating two direct-helicity coordinates in a basis
// All nonzero LS tensors must give the same coordinate phase
//
DirectHelicityPhase DirectHelicityCoordinatePhase(const std::vector<std::vector<std::complex<double>>> &basis_vectors,
                                                  std::size_t from_index, std::size_t to_index, double tolerance) {
  DirectHelicityPhase out;
  const double        active_tol2 = tolerance * tolerance;

  for (const auto &v : basis_vectors) {
    if (from_index >= v.size() || to_index >= v.size()) {
      throw std::invalid_argument(
          "gra::spin::DirectHelicityCoordinatePhase: "
          "coordinate outside basis vector");
    }
    const std::complex<double> a        = v[from_index];
    const std::complex<double> b        = v[to_index];
    const bool                 a_active = gra::math::abs2(a) > active_tol2;
    const bool                 b_active = gra::math::abs2(b) > active_tol2;
    if (!a_active && !b_active) { continue; }
    if (a_active != b_active) { return {}; }

    std::complex<double> ratio     = b / a;
    const double         ratio_abs = std::abs(ratio);
    if (ratio_abs < tolerance) { return {}; }
    ratio /= ratio_abs;

    if (!out.found) {
      out.phase = ratio;
      out.found = true;
    } else if (std::abs(out.phase - ratio) > tolerance) {
      return {};
    }
  }

  return out;
}

}  // namespace spin

namespace spin {

// Project a direct helicity vector into one orthonormal basis
// This removes components outside the allowed JW tensor space
//
std::vector<std::complex<double>> ProjectDirectHelicityVector(
    const std::vector<std::vector<std::complex<double>>> &basis, const std::vector<std::complex<double>> &seed) {
  return gra::ProjectOrthonormal(basis, seed);
}

// Build a flat direct-helicity target projected into an allowed JW basis
// Flat-helicity mode chooses equal phase support after LS projection
//
std::vector<std::complex<double>> BuildFlatDirectHelicityVector(const DirectHelicityBasis &basis, double J, double s1,
                                                                double s2, const std::string &context) {
  const std::size_t                 n1 = SpinStateCount(s1, context + " leg 1");
  const std::size_t                 n2 = SpinStateCount(s2, context + " leg 2");
  std::vector<std::complex<double>> target(n1 * n2, 0.0);
  std::vector<bool>                 assigned(n1 * n2, false);

  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      const double lambda1 = -s1 + static_cast<double>(i);
      const double lambda2 = -s2 + static_cast<double>(j);
      if (std::abs(lambda1 - lambda2) > J + kSpinLabelTolerance) { continue; }
      const std::size_t seed_index = i * n2 + j;
      if (assigned[seed_index]) { continue; }

      std::vector<std::complex<double>> seed(n1 * n2, 0.0);
      seed[seed_index] = 1.0;
      auto projected   = ProjectDirectHelicityVector(basis.basis, seed);
      if (gra::SquaredNorm(projected) <= 1e-20) { continue; }
      NormalizeDirectHelicityVector(projected);
      CanonicalizeDirectHelicityVectorPhase(projected);

      for (const auto &k : indices(projected)) {
        if (std::abs(projected[k]) <= 1e-12 || assigned[k]) { continue; }
        target[k]   = projected[k] / std::abs(projected[k]);
        assigned[k] = true;
      }
    }
  }

  bool has_target = false;
  for (const auto &x : target) {
    if (std::abs(x) > 1e-12) {
      has_target = true;
      break;
    }
  }
  if (!has_target) { throw std::invalid_argument(context + ": could not build a flat direct-helicity target"); }

  auto selected = ProjectDirectHelicityVector(basis.basis, target);
  if (gra::SquaredNorm(selected) <= 1e-20) {
    throw std::invalid_argument(context +
                                ": flat direct-helicity target is orthogonal "
                                "to the allowed LS subspace");
  }
  NormalizeDirectHelicityVector(selected);
  CanonicalizeDirectHelicityVectorPhase(selected);
  return selected;
}

namespace {

// Initialize direct helicity dimensions and projection labels in a HELMatrix
// The reduced matrix stores H(lambda1,lambda2)
//
void PrepareDirectHelicityMatrix(gra::HELMatrix &hc, double J, const gra::MParticle &p1, const gra::MParticle &p2,
                                 bool production_mode, VertexContext context) {
  const double      s1 = p1.spinX2 / 2.0;
  const double      s2 = p2.spinX2 / 2.0;
  const std::size_t n1 = SpinStateCount(s1, "gra::spin::PrepareDirectHelicityMatrix leg 1");
  const std::size_t n2 = SpinStateCount(s2, "gra::spin::PrepareDirectHelicityMatrix leg 2");
  hc.coupling_basis    = gra::CouplingBasis::Helicity;
  InitTwoBodyBasis(hc, J, s1, s2, SpinProjections(J),
                   VertexLegHelicities(p1, production_mode, context, 0, "gra::spin::ValidateDirectTMatrix daughter 1"),
                   VertexLegHelicities(p2, production_mode, context, 1, "gra::spin::ValidateDirectTMatrix daughter 2"),
                   "gra::spin::ValidateDirectTMatrix");
  hc.T     = MMatrix<std::complex<double>>(n1, n2, 0.0);
  hc.T_set = MMatrix<bool>(n1, n2, false);
}

}  // namespace

// Store a flattened direct-helicity vector into a prepared HELMatrix
// Flattened amplitude vectors are reshaped into helicity rows
//
void StoreDirectHelicityVector(gra::HELMatrix &hc, const std::vector<std::complex<double>> &vector) {
  const std::size_t n1 = hc.T.size_row();
  const std::size_t n2 = hc.T.size_col();
  if (vector.size() != n1 * n2) {
    throw std::invalid_argument("gra::spin::StoreDirectHelicityVector: vector size mismatch");
  }
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      hc.T[i][j]     = vector[i * n2 + j];
      hc.T_set[i][j] = std::abs(hc.T[i][j]) > 1e-12;
    }
  }
}

// Build one automatic two-body central coupling in LS or direct-helicity basis
// Automatic modes choose a relative shape in the allowed subspace
//
gra::HELMatrix BuildAutomaticCentralCoupling(const gra::MParticle              &mother,
                                             const std::vector<gra::MParticle> &vertex_legs,
                                             AutoCentralCouplingMode mode, bool production_mode,
                                             const std::string &source_name, const std::string &verbose_label,
                                             bool verbose_output, bool p_symmetry, bool c_symmetry,
                                             VertexContext context) {
  if (vertex_legs.size() != 2) {
    throw std::invalid_argument("gra::spin::BuildAutomaticCentralCoupling: expected two vertex legs");
  }

  gra::HELMatrix hc;
  hc.BR         = 1.0;
  hc.P_symmetry = p_symmetry;
  hc.C_symmetry = c_symmetry;
  hc.alpha_ls.Clear();

  const auto allowed_ls =
      AllowedLSCouplings(hc, mother, vertex_legs[0], vertex_legs[1], production_mode, source_name, context);
  if (allowed_ls.empty()) {
    throw std::invalid_argument(
        "gra::spin::BuildAutomaticCentralCoupling: no "
        "allowed LS coupling for PDG = " +
        std::to_string(mother.pdg) + " from " + source_name);
  }

  if (mode == AutoCentralCouplingMode::MinL || mode == AutoCentralCouplingMode::MinS ||
      mode == AutoCentralCouplingMode::FlatLS) {
    if (mode == AutoCentralCouplingMode::MinL || mode == AutoCentralCouplingMode::MinS) {
      // Select the requested primary LS quantum number and break ties with the
      // other
      const auto row =
          *std::min_element(allowed_ls.cbegin(), allowed_ls.cend(), [mode](const auto &left, const auto &right) {
            if (mode == AutoCentralCouplingMode::MinL) {
              return std::tie(left.l, left.two_s) < std::tie(right.l, right.two_s);
            }
            return std::tie(left.two_s, left.l) < std::tie(right.two_s, right.l);
          });
      hc.alpha_ls.Set(row.l, row.two_s, 1.0);
    } else {
      for (const auto &row : allowed_ls) { hc.alpha_ls.Set(row.l, row.two_s, 1.0); }
    }
    InitTMatrix(hc, mother, vertex_legs[0], vertex_legs[1], production_mode, source_name, false, verbose_output,
                verbose_label, context);
    return hc;
  }

  const double J     = mother.spinX2 / 2.0;
  const double s1    = vertex_legs[0].spinX2 / 2.0;
  const double s2    = vertex_legs[1].spinX2 / 2.0;
  const auto   basis = BuildJWDirectHelicityBasis(hc, mother, vertex_legs[0], vertex_legs[1], production_mode,
                                                  source_name, verbose_label, context);
  if (basis.basis.empty()) {
    throw std::invalid_argument(
        "gra::spin::BuildAutomaticCentralCoupling: no "
        "allowed direct-helicity basis for PDG = " +
        std::to_string(mother.pdg) + " from " + source_name);
  }

  std::vector<std::complex<double>> selected;
  if (mode != AutoCentralCouplingMode::FlatHelicity) {
    throw std::invalid_argument(
        "gra::spin::BuildAutomaticCentralCoupling: "
        "unsupported automatic coupling mode");
  }
  selected = BuildFlatDirectHelicityVector(basis, J, s1, s2, "gra::spin::BuildAutomaticCentralCoupling flat helicity");

  PrepareDirectHelicityMatrix(hc, J, vertex_legs[0], vertex_legs[1], production_mode, context);
  StoreDirectHelicityVector(hc, selected);
  ValidateDirectTMatrix(hc, mother, vertex_legs[0], vertex_legs[1], production_mode, source_name, false, verbose_output,
                        verbose_label, context);
  return hc;
}

namespace {

// Store the common physical description of one two-body vertex
struct VertexDescription {
  bool        crossed    = false;
  bool        subchannel = false;
  bool        physical_c = false;
  std::string process;
  std::string role0;
  std::string role1;
  std::string role2;
  std::string arrow;
  std::string coupling;
  std::string steering;
  std::string descriptor;
};

// Build the common physical description of one two-body vertex
VertexDescription DescribeVertex(const std::string &caller, const gra::HELMatrix &hel, const gra::MParticle &mother,
                                 const gra::MParticle &leg1, const gra::MParticle &leg2, int Jx2, int j1x2, int j2x2,
                                 bool production_mode, const std::string &source_name, VertexContext context,
                                 bool include_identical) {
  VertexDescription out;
  out.crossed    = production_mode && context == VertexContext::CrossedBeamLeg;
  out.subchannel = production_mode && context == VertexContext::SubTUChannelExchange;
  out.physical_c = hel.C_symmetry && IsValidCParity(mother.C) && !(out.crossed || out.subchannel);
  out.process    = out.crossed       ? "crossed beam-leg t-channel"
                   : out.subchannel  ? "sub t/u-channel exchange"
                   : production_mode ? "production"
                                     : "decay";
  out.role0      = out.crossed       ? "exchange state"
                   : out.subchannel  ? "external exchange"
                   : production_mode ? "produced state"
                                     : "mother state";
  out.role1      = out.crossed       ? "beam particle"
                   : out.subchannel  ? "outgoing particle"
                   : production_mode ? "leg 1"
                                     : "daughter 1";
  out.role2      = out.crossed       ? "scattered beam"
                   : out.subchannel  ? "internal off-shell"
                   : production_mode ? "leg 2"
                                     : "daughter 2";
  out.arrow      = production_mode ? "[state 0 <- state 1 + state 2]" : "[state 0 -> state 1 + state 2]";
  out.coupling   = production_mode ? "production" : "decay";
  out.steering   = source_name.empty() ? (production_mode ? "central continuum card" : "DECAYS.json") : source_name;
  const std::string compact_arrow = production_mode ? " <- " : " -> ";
  out.descriptor = caller + ":: [" + gra::aux::Spin2XtoString(Jx2) + "^" + gra::aux::ParityToString(mother.P) +
                   compact_arrow + gra::aux::Spin2XtoString(j1x2) + "^" + gra::aux::ParityToString(leg1.P) + " " +
                   gra::aux::Spin2XtoString(j2x2) + "^" + gra::aux::ParityToString(leg2.P) + "]";
  if (include_identical) {
    const bool identical = PhysicalIdenticalPair(leg1, leg2, production_mode, context);
    out.descriptor += " (identical = " + std::string(identical ? "true" : "false") + ")";
  }
  out.descriptor += " [PDG: " + std::to_string(mother.pdg) + compact_arrow + std::to_string(leg1.pdg) + " " +
                    std::to_string(leg2.pdg) + "]";
  return out;
}

// Print the common particle content and symmetry constraints of one vertex
void PrintVertexDescription(const std::string &caller, const VertexDescription &info, const gra::HELMatrix &hel,
                            const gra::MParticle &mother, const gra::MParticle &leg1, const gra::MParticle &leg2,
                            const std::string &verbose_label, bool direct) {
  const auto row = [](const std::string &state, const std::string &role, bool crossed, const gra::MParticle &particle) {
    return std::vector<std::string>{state,
                                    role,
                                    crossed ? "yes" : "no",
                                    crossed ? "antiparticle data" : "physical data",
                                    crossed ? "lambda_c=-lambda_phys" : "identity",
                                    particle.spinX2 % 2 == 0 ? "boson" : "fermion",
                                    gra::aux::NullableSpin2XtoString(particle.spinX2),
                                    gra::aux::ParityToString(particle.P),
                                    gra::aux::ParityToString(particle.C),
                                    std::to_string(particle.pdg),
                                    particle.name};
  };
  std::cout << std::endl << caller << ": " << info.process << " helicity decomposition" << std::endl << std::endl;
  if (!verbose_label.empty()) { std::cout << "Vertex: " << verbose_label << std::endl << std::endl; }
  std::cout << info.arrow << std::endl;
  auto rows =
      std::vector<std::vector<std::string>>{row("0", info.role0, false, mother), row("1", info.role1, false, leg1),
                                            row("2", info.role2, info.crossed || info.subchannel, leg2)};
  if (info.subchannel && !rows[2][10].empty()) { rows[2][10] += "*"; }
  gra::aux::PrintTable({"ID", "Role", "Crossed", "PDG basis", "Helicity map", "Type", "J", "P", "C", "PDG", "Name"},
                       rows);
  std::cout << std::endl;
  if (info.crossed || info.subchannel) {
    std::cout << "Crossed leg 2:" << std::endl
              << "  particle data: antiparticle" << std::endl
              << "  rows: lambda_c=-lambda_phys, phase "
                 "(-1)^(s2-lambda_phys)"
              << std::endl
              << std::endl;
  }
  if (hel.P_symmetry) {
    std::cout << (direct ? "Parity constraint: direct H must lie in the "
                           "allowed JW parity subspace"
                         : "Parity constraint: P_0 = P_1 x P_2 x (-1)^l")
              << std::endl;
  }
  if (info.physical_c) {
    std::cout << (direct ? "C-parity constraint: direct H must lie in the "
                           "allowed JW C subspace"
                         : "C-parity constraint: C_0 = C(two-body state)")
              << std::endl;
  }
  std::cout << std::endl;
}

// Store validated spin labels and dimensions for one two-body vertex
struct VertexSpin {
  double      J;
  double      s1;
  double      s2;
  int         Jx2;
  int         s1x2;
  int         s2x2;
  std::size_t n1;
  std::size_t n2;
};

// Initialize the three irreducible representations of one two-body vertex
VertexSpin InitVertexSpin(const gra::MParticle &mother, const gra::MParticle &leg1, const gra::MParticle &leg2,
                          const std::string &context) {
  const SpinRep parent(mother.spinX2, context + " mother");
  const SpinRep first(leg1.spinX2, context + " daughter 1");
  const SpinRep second(leg2.spinX2, context + " daughter 2");
  return {parent.Spin(), first.Spin(), second.Spin(), parent.X2(), first.X2(), second.X2(), first.Dim(), second.Dim()};
}

}  // namespace

// Initialize helicity decay amplitude matrix T with (2s1 + 1) x (2s2 + 1) dimensions
// This is the Jacob-Wick LS -> helicity recoupling
//
// where s1, s2 are daughter spins (0, 1/2, 1, ...)
//
// T_{\lambda_1,\lambda_2}
//   = \sum_{ls} \alpha_{ls} x
//          < J\lambda|ls0\lambda> <s\lambda|s1s2\lambda1,-\lambda2>
//
void InitTMatrix(gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1, const gra::MParticle &p2,
                 bool production_mode, const std::string &source_name, bool require_complete_l2s, bool verbose,
                 const std::string &verbose_label, VertexContext context) {
  hc.gp_orbital                                  = {};
  hc.jw_rotation                                 = {};
  const auto [J, s1, s2, J2, s1x2, s2x2, n1, n2] = InitVertexSpin(p, p1, p2, "gra::spin::InitTMatrix");

  const auto         info = DescribeVertex("gra::spin::InitTMatrix", hc, p, p1, p2, J2, s1x2, s2x2, production_mode,
                                           source_name, context, true);
  const std::string &coupling_label     = info.coupling;
  const std::string &steering_label     = info.steering;
  const std::string &process_descriptor = info.descriptor;

  // -------------------------------------------------------------------
  // Helicity decay amplitude matrix

  hc.T = MMatrix<std::complex<double>>(static_cast<unsigned int>(n1), static_cast<unsigned int>(n2), 0.0);
  std::vector<LSCoupling> missing_ls;
  unsigned int            nonzero = 0;

  // Construct T-matrix (2s1 + 1) x (2s2 + 1) elements
  if (verbose) {
    PrintVertexDescription("gra::spin::InitTMatrix", info, hc, p, p1, p2, verbose_label, false);
    std::cout << "gra::spin::InitTMatrix: SU(2) decomposition"
              << " [lambda = lambda1 - lambda2]" << std::endl;
  }

  std::vector<std::vector<std::string>> su2_rows;
  struct ActiveLSBasis {
    std::size_t                   l     = 0;
    std::size_t                   two_s = 0;
    MMatrix<std::complex<double>> matrix;
  };
  std::vector<ActiveLSBasis> bases;
  const std::size_t          max_two_s = static_cast<std::size_t>(s1x2 + s2x2);
  for (std::size_t two_s = 0; two_s <= max_two_s; ++two_s) {
    const double s = 0.5 * static_cast<double>(two_s);
    if (!wigner::IsInt(s1 + s2 + s)) { continue; }
    const std::size_t max_l = static_cast<std::size_t>(std::floor(J + s + kSpinLabelTolerance));
    for (std::size_t l = 0; l <= max_l; ++l) {
      MMatrix<std::complex<double>> basis = LSHelicityBasis(hc, p, p1, p2, l, two_s, production_mode, context);
      if (basis.FrobNorm2() <= 1e-24) { continue; }
      ++nonzero;
      if (require_complete_l2s && !hc.alpha_ls.Contains(l, two_s)) { missing_ls.push_back({l, two_s}); }
      bases.push_back({l, two_s, std::move(basis)});
    }
  }

  for (const LSTerm &term : hc.alpha_ls) {
    const bool active = std::any_of(bases.cbegin(), bases.cend(), [&term](const ActiveLSBasis &basis) {
      return basis.l == term.l && basis.two_s == term.two_s;
    });
    if (!active) {
      const std::string str = ": Forbidden helicity " + coupling_label + " alpha_ls-coupling from " + steering_label +
                              ": (l,s) = (" + std::to_string(term.l) + "," +
                              gra::aux::Spin2XtoString(static_cast<int>(term.two_s)) + ")";
      throw LSCouplingRejected(process_descriptor + str);
    }
  }

  const double sum_alphasq = hc.alpha_ls.Norm2();
  if (sum_alphasq < 1e-12) {
    throw LSCouplingRejected(process_descriptor + ": All active helicity " + coupling_label +
                             " alpha_ls-couplings are zero "
                             "(sum_{ls} |alpha_ls|^2 = 0)");
  }
  if (!production_mode && std::abs(sum_alphasq - 1.0) > 1e-6) {
    if (verbose) {
      std::cout << "Renormalizing input: sum_{ls} |alpha_{ls}|^2 = " << sum_alphasq << " => 1" << std::endl
                << std::endl;
    }
    hc.alpha_ls.Scale(1.0 / msqrt(sum_alphasq));
  }

  hc.ls_components.clear();
  for (const auto &basis : bases) {
    const LSTerm *term = hc.alpha_ls.Find(basis.l, basis.two_s);
    if (term == nullptr || std::fpclassify(math::abs2(term->coefficient)) == FP_ZERO) { continue; }
    const std::complex<double>    alpha    = term->coefficient;
    MMatrix<std::complex<double>> physical = basis.matrix;
    try {
      ProjectFinalMasslessVectorHelicities(physical, p1, p2, production_mode, context, process_descriptor);
    } catch (const std::invalid_argument &) { continue; }
    hc.T += physical * alpha;
    hc.ls_components.push_back({basis.l, basis.two_s, alpha, std::move(physical)});
    if (!verbose) { continue; }
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n2; ++j) {
        if (std::abs(basis.matrix[i][j]) < 1e-12) { continue; }
        su2_rows.push_back(
            {"(" + std::to_string(basis.l) + "," + gra::aux::Spin2XtoString(static_cast<int>(basis.two_s)) + ")",
             FormatSpinLabel(-s1 + static_cast<double>(i)), FormatSpinLabel(-s2 + static_cast<double>(j)),
             gra::aux::ToString(std::real(basis.matrix[i][j]), 3), gra::aux::ToString(std::abs(alpha), 3),
             gra::aux::ToString(std::arg(alpha), 3)});
      }
    }
  }

  // Real external photons and gluons carry only lambda=+-1
  // The zero-norm check enforces Landau-Yang and related helicity selection
  // rules
  const double projected_norm2 = hc.T.FrobNorm2();
  ProjectFinalMasslessVectorHelicities(hc.T, p1, p2, production_mode, context, process_descriptor);
  if (!production_mode && projected_norm2 > 0.0) {
    const double scale = std::sqrt(hc.T.FrobNorm2() / projected_norm2);
    for (auto &component : hc.ls_components) { component.alpha *= scale; }
  }

  const std::string str0 = process_descriptor;

  if (nonzero == 0) {
    const std::string str = ": Process is impossible (spin-parity-C-statistics) [.P_symmetry = " +
                            (hc.P_symmetry ? std::string("true") : std::string("false")) +
                            ", .C_symmetry = " + (hc.C_symmetry ? std::string("true]") : std::string("false]"));
    throw LSCouplingRejected(str0 + str);
  }

  if (verbose) {
    gra::aux::PrintTable({"(l,s)", "lambda1", "lambda2", "JW", "|a_ls|", "arg(a_ls)"}, su2_rows);
    std::cout << std::endl;
  }

  if (verbose) {
    std::cout << std::endl;
    std::cout << "T matrix (2s1 + 1) x (2s2 + 1):" << std::endl;
    hc.T.PrintSeparate();
    PrintDirectHelicityCoefficientTable(hc, s1, s2, "Equivalent direct-helicity coefficients from LS input:");
  }

  hc.alpha_ls.RemoveBelow(0.0);

  if (verbose) {
    std::cout << (production_mode ? "Production coupling verification:" : "Normalization verification:") << std::endl;
    if (production_mode) {
      std::cout << "sum |g_{LS}|^2 = " << hc.alpha_ls.Norm2() << ", sum |H_{lambda1,lambda2}|^2 = " << hc.T.FrobNorm2();
    } else {
      std::cout << "sum |alpha_{ls}|^2 = " << hc.alpha_ls.Norm2()
                << "  <=> sum |T_{lambda1,lambda2}|^2 = " << hc.T.FrobNorm2();
    }
    std::cout << std::endl << std::endl;
  }

  // Print missing but allowed alpha_ls rows without failing initialization
  if (!missing_ls.empty()) {
    std::string str = "";
    for (const LSCoupling &row : missing_ls) {
      str += "(" + std::to_string(row.l) + "," + gra::aux::Spin2XtoString(static_cast<int>(row.two_s)) + ") ";
    }

    const std::string middle = ", WARNING: Following helicity " + coupling_label +
                               " alpha_ls-couplings "
                               "missing from " +
                               steering_label +
                               ": (l,s) "
                               "= ";
    std::cout << rang::fg::red << str0 + middle + str << rang::fg::reset << std::endl;
  }

  // ------------------------------------------------------------
  // Init indexing structures

  InitTwoBodyBasis(hc, J, s1, s2, SpinProjections(J),
                   VertexLegHelicities(p1, production_mode, context, 0, "gra::spin::InitTMatrix daughter 1"),
                   VertexLegHelicities(p2, production_mode, context, 1, "gra::spin::InitTMatrix daughter 2"),
                   "gra::spin::InitTMatrix");
  // ------------------------------------------------------------
}

// Validate direct reduced helicity amplitudes against the allowed Jacob-Wick subspace
// Explicit H(lambda1,lambda2) rows must be expressible by legal LS tensors
//
void ValidateDirectTMatrix(gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1,
                           const gra::MParticle &p2, bool production_mode, const std::string &source_name,
                           bool require_complete_helicity, bool verbose, const std::string &verbose_label,
                           VertexContext context) {
  hc.gp_orbital                                  = {};
  hc.jw_rotation                                 = {};
  const auto [J, s1, s2, J2, s1x2, s2x2, n1, n2] = InitVertexSpin(p, p1, p2, "gra::spin::ValidateDirectTMatrix");

  const auto info = DescribeVertex("gra::spin::ValidateDirectTMatrix", hc, p, p1, p2, J2, s1x2, s2x2, production_mode,
                                   source_name, context, false);
  const std::string &coupling_label     = info.coupling;
  const std::string &steering_label     = info.steering;
  const std::string &process_descriptor = info.descriptor;

  if (verbose) { PrintVertexDescription("gra::spin::ValidateDirectTMatrix", info, hc, p, p1, p2, verbose_label, true); }

  if (hc.T.size_row() != n1 || hc.T.size_col() != n2) {
    throw std::invalid_argument(process_descriptor + ": direct helicity T matrix has invalid dimensions");
  }
  if (hc.T_set.isEmpty()) {
    hc.T_set = MMatrix<bool>(n1, n2, false);
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n2; ++j) {
        if (std::abs(hc.T[i][j]) > 1e-12) { hc.T_set[i][j] = true; }
      }
    }
  }
  if (hc.T_set.size_row() != n1 || hc.T_set.size_col() != n2) {
    throw std::invalid_argument(process_descriptor + ": direct helicity T_set matrix has invalid dimensions");
  }

  InitTwoBodyBasis(hc, J, s1, s2, SpinProjections(J),
                   VertexLegHelicities(p1, production_mode, context, 0, "gra::spin::ValidateDirectTMatrix daughter 1"),
                   VertexLegHelicities(p2, production_mode, context, 1, "gra::spin::ValidateDirectTMatrix daughter 2"),
                   "gra::spin::ValidateDirectTMatrix");

  const auto subspace = BuildDirectHelicitySubspace(hc, p, p1, p2, production_mode, context);
  if (subspace.basis.empty()) {
    throw std::invalid_argument(process_descriptor + ": no allowed Jacob-Wick LS subspace exists");
  }

  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      if (hc.T_set[i][j] && !subspace.support[i][j]) {
        const double lambda1 = -s1 + static_cast<double>(i);
        const double lambda2 = -s2 + static_cast<double>(j);
        throw std::invalid_argument(process_descriptor + ": Forbidden helicity " + coupling_label + " amplitude from " +
                                    steering_label + ": (lambda1,lambda2) = (" + FormatSpinLabel(lambda1) + "," +
                                    FormatSpinLabel(lambda2) + ")");
      }
    }
  }

  if (require_complete_helicity) {
    std::vector<std::pair<double, double>> missing;
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n2; ++j) {
        if (subspace.support[i][j] && !hc.T_set[i][j]) {
          missing.push_back({-s1 + static_cast<double>(i), -s2 + static_cast<double>(j)});
        }
      }
    }
    if (!missing.empty()) {
      std::ostringstream message;
      message << process_descriptor << ": helicity " << coupling_label << " table from " << steering_label
              << " omits independent allowed coordinates";
      for (const auto &row : missing) {
        message << " (" << FormatSpinLabel(row.first) << "," << FormatSpinLabel(row.second) << ")";
      }
      throw std::invalid_argument(message.str());
    }
  }

  const double norm2_before = hc.T.FrobNorm2();
  if (!production_mode && norm2_before < 1e-12) {
    throw std::invalid_argument(process_descriptor + ": All active helicity " + coupling_label +
                                " amplitudes are zero (sum |H|^2 = 0)");
  }

  std::vector<std::complex<double>> residual = hc.T.Flatten();
  for (const auto &basis_vector : subspace.basis) {
    const std::complex<double> c = gra::InnerProduct(basis_vector, residual);
    gra::AddScaled(residual, basis_vector, -c);
  }
  const double residual2 = gra::SquaredNorm(residual);
  if (residual2 > 1e-10 * std::max(1.0, norm2_before)) {
    std::ostringstream ss;
    ss << process_descriptor << ": direct helicity " << coupling_label << " amplitudes from " << steering_label
       << " are not compatible with the allowed Jacob-Wick LS subspace"
       << " (residual norm^2 = " << residual2 << ")";
    throw std::invalid_argument(ss.str());
  }

  if (!production_mode && std::abs(norm2_before - 1.0) > 1e-6) {
    const double inv_norm = 1.0 / std::sqrt(norm2_before);
    hc.T *= inv_norm;
    if (verbose) {
      std::cout << "Renormalizing input: sum |H_{lambda1,lambda2}|^2 = " << norm2_before << " => 1" << std::endl
                << std::endl;
    }
  }

  const auto coefficients = EquivalentLSCoefficients(subspace, hc.T.Flatten());
  if (subspace.raw_matrices.size() != coefficients.size()) {
    throw std::logic_error("gra::spin::ValidateDirectTMatrix: LS basis cache size mismatch");
  }
  hc.ls_components.clear();
  for (const auto &i : indices(coefficients)) {
    if (std::abs(coefficients[i]) <= 1e-12) { continue; }
    hc.ls_components.push_back(
        {subspace.raw_rows[i].l, subspace.raw_rows[i].two_s, coefficients[i], subspace.raw_matrices[i]});
  }

  if (verbose) {
    PrintEquivalentLSCoefficientTable(subspace, hc.T.Flatten(),
                                      "Equivalent LS coefficients from direct-helicity input:");

    std::vector<std::vector<std::string>> rows;
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n2; ++j) {
        if (!hc.T_set[i][j] && std::abs(hc.T[i][j]) <= 1e-12) { continue; }
        const double lambda1 = -s1 + static_cast<double>(i);
        const double lambda2 = -s2 + static_cast<double>(j);
        rows.push_back({FormatSpinLabel(lambda1), FormatSpinLabel(lambda2), gra::aux::ToString(std::abs(hc.T[i][j]), 3),
                        gra::aux::ToString(std::arg(hc.T[i][j]), 3), subspace.support[i][j] ? "yes" : "no",
                        hc.T_set[i][j] ? "yes" : "no"});
      }
    }
    std::cout << std::endl;
    std::cout << "gra::spin::ValidateDirectTMatrix: direct helicity basis" << std::endl;
    if (!verbose_label.empty()) { std::cout << "Vertex: " << verbose_label << std::endl; }
    gra::aux::PrintTable({"lambda1", "lambda2", "|H|", "arg(H)", "allowed", "set"}, rows);
    std::cout << std::endl;
    std::cout << "T matrix (direct reduced helicity basis, rows lambda1 x lambda2):" << std::endl;
    hc.T.PrintSeparate();
    std::cout << (production_mode ? "Production coupling verification:" : "Normalization verification:") << std::endl;
    if (production_mode) {
      std::cout << "sum |H_{lambda1,lambda2}|^2 = " << hc.T.FrobNorm2()
                << ", allowed JW subspace dimension = " << subspace.basis.size();
    } else {
      std::cout << "sum |H_{lambda1,lambda2}|^2 = " << hc.T.FrobNorm2()
                << ", allowed JW subspace dimension = " << subspace.basis.size();
    }
    std::cout << std::endl << std::endl;
  }
}

// Compute the identical-boson exchange selection rule
// Bose symmetry requires l-S to be even
bool BoseSymmetry(int l, int s) {
  // l - s must be even for symmetric wavefunction
  return (l - s) % 2 == 0;
}

// Compute the identical-fermion exchange selection rule
// Fermi symmetry requires l+S to be even
bool FermiSymmetry(int l, int s) {
  // l + s must be even for antisymmetric wavefunction
  return (l + s) % 2 == 0;
}

// -------------------------------------------------------

// Compute the pole-normalized helicity matrix with LS threshold factors
MMatrix<std::complex<double>> DecayLSHelicityMatrix(const gra::HELMatrix &hel, double q_ratio, bool barrier) {
  if (!std::isfinite(q_ratio) || q_ratio < 0.0) {
    throw AmplitudeFailure("DecayLSHelicityMatrix: q/q0 must be finite and non-negative");
  }
  if (!barrier || std::abs(q_ratio - 1.0) <= 1e-14) { return hel.T; }

  MMatrix<std::complex<double>> out(hel.T.size_row(), hel.T.size_col(), 0.0);
  for (const auto &component : hel.ls_components) {

    out += component.matrix * (component.alpha * std::pow(q_ratio, component.l));
  }
  return out;
}

// Compute the angle-integrated LS threshold intensity
double DecayLSIntensity(const gra::HELMatrix &hel, double q_ratio, bool barrier) {

  if (!barrier) { return 1.0; }

  const double pole_norm = hel.T.FrobNorm2();

  // The physical helicity trace retains overlaps of transverse LS operators
  return DecayLSHelicityMatrix(hel, q_ratio, true).FrobNorm2() / pole_norm;
}

// Precompute the finite Jacob-Wick helicity row map
void InitJWRotation(gra::HELMatrix &hel) {
  const std::size_t rows = hel.lambda_values.size_row();
  const std::size_t cols = hel.Jz_values.size();
  if (cols == 0 || hel.lambda_values.size_col() != 2 || hel.lambda_idx.size_col() != 2 ||
      hel.lambda_idx.size_row() != rows) {
    throw std::invalid_argument("InitJWRotation: invalid helicity projection dimensions");
  }
  const double J = hel.J >= 0.0 ? hel.J : static_cast<double>(cols - 1) / 2.0;
  if (J < 0.0) { throw std::invalid_argument("InitJWRotation: negative mother spin"); }

  std::vector<std::size_t> mother_indices;
  mother_indices.reserve(cols);
  for (const double Jz : hel.Jz_values) {
    const std::size_t index = SpinProjectionIndex(Jz, J, "InitJWRotation mother projection");
    if (std::find(mother_indices.cbegin(), mother_indices.cend(), index) != mother_indices.cend()) {
      throw std::invalid_argument("InitJWRotation: duplicate mother spin projection");
    }
    mother_indices.push_back(index);
  }

  JWRotationBasis basis;
  basis.row.resize(rows);
  basis.difference.reserve(rows);
  for (const auto &row : indices(basis.row)) {
    const double lambda = hel.lambda_values[row][0] - hel.lambda_values[row][1];
    if (!std::isfinite(lambda)) { throw std::invalid_argument("InitJWRotation: non-finite helicity difference"); }
    const auto found = std::find(basis.difference.cbegin(), basis.difference.cend(), lambda);
    if (found == basis.difference.cend()) {
      basis.row[row] = basis.difference.size();
      basis.difference.push_back(lambda);
    } else {
      basis.row[row] = static_cast<std::size_t>(found - basis.difference.cbegin());
    }
  }
  basis.rotation  = std::make_shared<const wigner::Rotation>(basis.difference, hel.Jz_values, J);
  basis.ready     = true;
  hel.jw_rotation = std::move(basis);
}

// Construct the event-dependent Jacob-Wick angular decay matrix
// The decay matrix is f = D* T in the mother spin basis
// [REFERENCE: Jacob, Wick, On the general theory of particles with spin, 1959]
// [REFERENCE: Amsler, Bizot, Simulation of angular distributions, 1983]
MMatrix<std::complex<double>> fDecayMatrix(const gra::HELMatrix &hel, double theta, double phi) {
  const std::size_t rows = hel.lambda_values.size_row();
  const std::size_t cols = hel.Jz_values.size();

  // Construct transition amplitude matrix
  MMatrix<std::complex<double>> f(rows, cols);

  const double J = (hel.J >= 0.0) ? hel.J : static_cast<double>(cols - 1) / 2.0;

  if (!std::isfinite(theta) || !std::isfinite(phi)) {
    throw AmplitudeFailure("fDecayMatrix: non-finite generated decay angle");
  }

  // The scalar representation is the identity and has no angular dependence
  if (std::fpclassify(J) == FP_ZERO) {
    for (std::size_t i = 0; i < rows; ++i) {
      const std::size_t i1 = hel.lambda_idx[i][0];
      const std::size_t i2 = hel.lambda_idx[i][1];
      const double lambda = hel.lambda_values[i][0] - hel.lambda_values[i][1];
      if (std::abs(lambda) <= kSpinLabelTolerance) { f[i][0] = hel.T[i1][i2]; }
    }
    return f;
  }

  std::vector<std::complex<double>> phase(cols);
  for (std::size_t j = 0; j < cols; ++j) { phase[j] = std::exp(std::complex<double>(0.0, hel.Jz_values[j] * phi)); }

  const MMatrix<double> rotation = hel.jw_rotation.rotation->Evaluate(theta);

  // Rows = final state spin projections
  for (std::size_t i = 0; i < rows; ++i) {
    const std::size_t i1 = hel.lambda_idx[i][0];
    const std::size_t i2 = hel.lambda_idx[i][1];
    const std::complex<double> T = hel.T[i1][i2];

    // Columns = initial state polarizations
    for (std::size_t j = 0; j < cols; ++j) {
      // Wigner D * Helicity amplitude
      f[i][j] = rotation(hel.jw_rotation.row[i], j) * phase[j] * T;
    }
  }

  return f;
}

// Compute a source density averaged over the incoming spin states
// The row trace is divided by the physical beam spin multiplicity
double SourceSpinAveragedDensity(const MMatrix<std::complex<double>> &matrix, std::size_t incoming_spin_states,
                                 const std::string &context) {
  if (incoming_spin_states == 0) {
    throw std::invalid_argument(context + ": incoming spin-state count should be positive");
  }
  const double rho = matrix.FrobNorm2() / static_cast<double>(incoming_spin_states);
  if (!std::isfinite(rho)) { throw AmplitudeFailure(context + ": source has non-finite spin density"); }
  return rho;
}

// Compute the unpolarized density carried by one exchange helicity column
// The density is averaged over the physical incoming spin states
double ForwardSourceColumnSpinAveragedDensity(const MMatrix<std::complex<double>> &matrix, std::size_t column,
                                              std::size_t incoming_spin_states, const std::string &context) {
  if (column >= matrix.size_col()) {
    throw std::invalid_argument(context + ": reference helicity column is outside source matrix");
  }
  if (incoming_spin_states == 0) {
    throw std::invalid_argument(context + ": incoming spin-state count should be positive");
  }
  const double density = matrix.ColumnNorm2(column) / static_cast<double>(incoming_spin_states);
  if (!std::isfinite(density)) {
    throw AmplitudeFailure(context + ": reference helicity column has non-finite spin density");
  }
  return density;
}

// Compute the von Neumann entropy measuring spin-state mixedness
// S = -Tr(rho log rho) = -sum_i lambda_i log lambda_i
double VonNeumannEntropy(const MMatrix<std::complex<double>> &rho) {
  try {
    return rho.SelfAdjointEntropy(epsilon);
  } catch (const std::exception &) { return -1.0; }
}

// Density matrix properties:
//
//  1. Tr[rho^2] <= 1     (equality only for a pure-state projector)
//  2. rho^\dagger = rho  (hermiticity)
//  3. Tr[rho] = 1        (normalization)
//  4. rho >= 0           (positivity, eigenvalues greater or equal to zero)
//
// Test positivity of an integer- or half-integer-spin density matrix,
// - Also checks the trivial fact that rho
//   has right dimensions given J, that is (2J+1)
// - Also checks that the normalization is ok
//
//
// The density matrix rho must be Hermitian, normalized and positive semidefinite
//
bool Positivity(const MMatrix<std::complex<double>> &rho, double J) {
  const std::size_t n        = rho.size_row();
  const std::size_t expected = SpinStateCount(J, "gra::spin::Positivity spin");

  // Dimensions
  if (n != expected || rho.size_col() != expected) {
    throw std::invalid_argument(
        "gra::spin::Positivity: Density matrix must be "
        "square with dimension 2J+1");
  }

  // Normalization
  const std::complex<double> tracerho  = rho.Trace();
  constexpr double           TRACE_TOL = 1e-6;
  if (std::abs(std::imag(tracerho)) > TRACE_TOL || std::abs(std::real(tracerho) - 1.0) > TRACE_TOL) {
    std::string str =
        "gra::spin::Positivity: Density matrix is not properly "
        "normalized (Tr[rho] = " +
        std::to_string(std::real(tracerho)) + " + i " + std::to_string(std::imag(tracerho)) + ") (should be 1) !";
    rho.PrintSeparate();
    throw std::invalid_argument(str);
  }

  constexpr double HERM_TOL = 1e-10;
  if (!rho.IsHermitian(HERM_TOL)) {
    rho.PrintSeparate();
    throw std::invalid_argument("gra::spin::Positivity: Density matrix is not Hermitian");
  }

  try {
    for (const double eigenvalue : rho.SelfAdjointEigenvalues(HERM_TOL)) {
      if (eigenvalue < -1e-10) {
        rho.PrintSeparate();
        throw std::invalid_argument(
            "gra::spin::Positivity: Eigenvalues of "
            "the density matrix are not >= 0 !");
      }
    }
  } catch (const std::runtime_error &) {
    rho.PrintSeparate();
    throw std::invalid_argument("gra::spin::Positivity: Eigenvalue decomposition failed");
  }
  return true;
}

// Decompose a positive density matrix into weighted pure-state spin vectors
// The decomposition is rho = sum |psi_k><psi_k| with sqrt eigenvalue in each psi_k
//
std::vector<std::vector<std::complex<double>>> SpectralPureStates(const MMatrix<std::complex<double>> &rho) {
  constexpr double tol = 1e-10;
  try {
    return rho.PositiveSpectralVectors(tol);
  } catch (const std::domain_error &) {
    throw std::invalid_argument(
        "gra::spin::SpectralPureStates: Spin density "
        "matrix has negative eigenvalues");
  } catch (const std::exception &) {
    throw std::invalid_argument(
        "gra::spin::SpectralPureStates: Failed to diagonalize spin density "
        "matrix");
  }
}

// Generate a random density matrix of spin weights
// Reflection averaging commutes with the parity of the production tensor
MMatrix<std::complex<double>> RandomRho(int spinX2, bool parity, MRandom &rng) {
  if (spinX2 < 0) { throw std::invalid_argument("gra::spin::RandomRho: negative spinX2"); }
  const std::size_t             n = static_cast<std::size_t>(spinX2 + 1);
  MMatrix<std::complex<double>> random(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) { random[i][j] = std::complex<double>(rng.G(0.0, 1.0), rng.G(0.0, 1.0)); }
  }

  MMatrix<std::complex<double>> rho = random * random.Dagger();
  if (parity) {
    MMatrix<std::complex<double>> parity_image(n, n, 0.0);
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) { parity_image[i][j] = rho[n - 1 - i][n - 1 - j]; }
    }
    rho = (rho + parity_image) * 0.5;
  }

  const double trace = std::real(rho.Trace());
  if (!(trace > 0.0) || !std::isfinite(trace)) {
    throw std::runtime_error("gra::spin::RandomRho: invalid random trace");
  }
  rho *= 1.0 / trace;
  (void)Positivity(rho, spinX2 / 2.0);
  return rho;
}

}  // namespace spin
}  // namespace gra
