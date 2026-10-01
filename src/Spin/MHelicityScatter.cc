// Production (scattering) side of helicity amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Spin/MHelicityScatter.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {

namespace spin {

namespace {

constexpr double kIntegerTolerance = 1e-9;

// Compute true if a floating spin projection is numerically an integer label
// Regge exchange helicities m are integer labels
//
bool IsIntegerProjection(double value) { return std::abs(value - std::round(value)) < kIntegerTolerance; }

// Compute one integer spin projection with validation
// Analytic and Regge vertices use integer m powers and phases
//
int IntegerProjection(double value, const std::string& context) {
  if (!IsIntegerProjection(value)) { throw std::invalid_argument(context + ": non-integer exchange helicity label"); }
  return static_cast<int>(std::llround(value));
}

// Compute selected physical helicity transitions from a HELMatrix row list
// Selected rows define incoming and outgoing beam helicities
//
std::vector<std::pair<double, double>> PhotonTransitionsFromRows(const HELMatrix&                hel,
                                                                 const std::vector<std::size_t>& rows) {
  std::vector<std::pair<double, double>> out;
  out.reserve(rows.size());
  for (const std::size_t row : rows) { out.push_back({hel.lambda_values[row][0], hel.lambda_values[row][1]}); }
  return out;
}

// Build separable forward factors in the selected exchanged-state basis
// The reduced Regge basis has one column while physical poles retain helicity
//
ForwardSourceFactors ForwardExchangeResidueFactors(const HELMatrix& hel, const M4Vec& incoming_beam,
                                                   const M4Vec& outgoing_system, bool second_exchange_daughter,
                                                   const std::vector<std::size_t>& rows, ForwardSpec forward) {
  const M4Vec  q = incoming_beam - outgoing_system;
  const double collider_transfer_phi =
      (second_exchange_daughter ? incoming_beam - outgoing_system : outgoing_system - incoming_beam).Phi();
  // Matrix form R = f r^T with f_(lambda_i,lambda_f) = F_(lambda_i-lambda_f)
  // Component form R_(lambda_i,lambda_f,m) = F_(lambda_i-lambda_f) R_m
  // lambda_i,f are beam helicities, m is exchange helicity and F contains only
  // the external helicity section R_m contains the barrier and azimuthal phase
  auto column_factor = forward.basis == ExchangeBasisType::ReducedRegge
                           ? HelVec{1.0}
                           : ExchangeHelicityResidue(hel.Jz_values, q, second_exchange_daughter, forward);

  // The fitted SOFT coupling supplies the physical proton flip magnitude
  HelVec row_factor(rows.size());
  for (const auto& i : indices(rows)) {
    const std::size_t row        = rows[i];
    const double      lambda_in  = hel.lambda_values[row][0];
    const double      lambda_out = hel.lambda_values[row][1];
    row_factor[i] = spin::SpinHalfForwardHelicitySectionFactor(lambda_in, lambda_out, collider_transfer_phi,
                                                               second_exchange_daughter);
  }

  return {std::move(row_factor), std::move(column_factor)};
}

}  // namespace

// Transform a pair into the requested production spin frame
void ProductionFrame(std::vector<M4Vec>& p, const LORENTZSCALAR& lts, const std::string& frame,
                     const M4Vec& helicity_dir) {
  if (frame == "HX") {
    kinematics::HXframe(p, helicity_dir);
    return;
  }

  if (frame == "CS") {
    kinematics::CSframe(p, lts.pfinal[0], lts.pbeam1, lts.pbeam2);
    kinematics::BoostToRestFrame(p, "MHelicityScatter::ProductionFrame CS");
    return;
  }

  if (frame == "CM") {
    kinematics::BoostToRestFrame(p, "MHelicityScatter::ProductionFrame CM");
    return;
  }

  throw std::invalid_argument("MHelicityScatter: Unsupported production spin frame '" + frame +
                              "' for JW production amplitudes; "
                              "allowed values are CS, HX, CM");
}

// Transport CM unit axes with the same rotations used for the root decay
HelAmp ProductionRotation(const LORENTZSCALAR& lts, const std::string& frame, double J) {
  std::vector<M4Vec> axes = {M4Vec(1.0, 0.0, 0.0, 0.0), M4Vec(0.0, 0.0, 1.0, 0.0)};
  if (frame == "HX") {
    for (auto& axis : axes) { kinematics::RotateHelicityAxes(axis, lts.pfinal[0]); }
  } else if (frame == "CS") {
    const auto upper = kinematics::BoostToRestFrame(lts.pbeam1, lts.pfinal[0], "ProductionRotation upper beam");
    const auto lower = kinematics::BoostToRestFrame(lts.pbeam2, lts.pfinal[0], "ProductionRotation lower beam");
    kinematics::LorentzFrame(axes, upper, lower, axes, frame, -1);
  } else if (frame != "CM") {
    throw std::invalid_argument("ProductionRotation: unsupported production spin frame '" + frame + "'");
  }
  return SpinFrame(axes[0], axes[1], J);
}

// Build the exchange-helicity part of one fixed-spin forward vertex
HelVec ExchangeHelicityResidue(const std::vector<double>& projections, const M4Vec& transfer,
                               bool second_exchange_daughter, ForwardSpec forward) {
  const double exchange_scale = math::msqrt(forward.s0);
  HelVec       residue(projections.size(), 0.0);
  for (const auto& column : indices(projections)) {
    const int m = IntegerProjection(projections[column], "MHelicityScatter::ExchangeHelicityResidue projection");
    residue[column] =
        spin::HelicityFactor(m, 0.0, transfer.Pt(), transfer.Phi(), exchange_scale, exchange_scale,
                             second_exchange_daughter, forward.mode == ForwardVertexMode::HelicityResidue);
  }
  return residue;
}

// Convert full or compact leg rows with negative helicity first to canonical
std::size_t CanonicalProtonPairSpinLayout::FromNegativeFirstLegRows(std::size_t upper_row, std::size_t lower_row,
                                                                    std::size_t upper_rows, std::size_t lower_rows) {
  if (upper_row >= upper_rows || lower_row >= lower_rows) {
    throw std::invalid_argument("CanonicalProtonPairSpinLayout: leg row is outside the basis");
  }
  if (upper_rows == 4 && lower_rows == 4) {
    return HardRow(upper_row / 2, lower_row / 2, upper_row % 2, lower_row % 2);
  }
  if (upper_rows == 2 && lower_rows == 2) { return CompactNoFlipHardRow(upper_row, lower_row); }
  throw std::invalid_argument(
      "CanonicalProtonPairSpinLayout: expected two compact or four full "
      "rows per proton leg");
}

// Build the producer row permutation for a full or compact proton pair
CanonicalProtonPairSpinLayout::KroneckerRowMap CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
    std::size_t upper_rows, std::size_t lower_rows) {
  if (upper_rows != lower_rows || (upper_rows != 2 && upper_rows != 4)) {
    throw std::invalid_argument(
        "CanonicalProtonPairSpinLayout: expected two compact or four full "
        "rows per proton leg");
  }
  KroneckerRowMap rows;
  rows.count = upper_rows * lower_rows;
  for (std::size_t upper_row = 0; upper_row < upper_rows; ++upper_row) {
    for (std::size_t lower_row = 0; lower_row < lower_rows; ++lower_row) {
      const std::size_t source_row = upper_row * lower_rows + lower_row;
      rows.destination[source_row] = FromNegativeFirstLegRows(upper_row, lower_row, upper_rows, lower_rows);
    }
  }
  return rows;
}

namespace {

// Compute true if this exchange branch is the EPA photon source, not a Regge vertex
// Photons use QED/EPA source matrices instead of Regge vertices
//
bool IsPhotonExchangeBranch(const MDecayBranch& branch) { return branch.p.pdg == PDG::PDG_gamma; }

// Select the physical forward no-flip beam rows after crossed-row remapping
// No-flip means incoming and outgoing physical helicities are equal
//
std::vector<std::size_t> HelicityConservingRows(const HELMatrix& hel) {
  if (hel.lambda_values.size_col() != 2 || hel.lambda_idx.size_col() != 2 ||
      hel.lambda_idx.size_row() != hel.lambda_values.size_row()) {
    throw std::invalid_argument(
        "MHelicityScatter: Beam-leg helicity metadata "
        "has inconsistent dimensions");
  }

  const double             s1      = hel.s1 >= 0.0 ? hel.s1 : SpinFromHelicityColumn(hel, 0);
  const auto               count   = SpinStateCount(s1, "HelicityConservingRows incoming spin");
  const auto               missing = hel.lambda_values.size_row();
  std::vector<std::size_t> rows(count, missing);
  for (std::size_t i = 0; i < missing; ++i) {
    if (std::abs(hel.lambda_values[i][0] - hel.lambda_values[i][1]) > kIntegerTolerance) { continue; }
    const auto incoming = SpinProjectionIndex(hel.lambda_values[i][0], s1, "HelicityConservingRows");
    if (rows[incoming] != missing) { throw std::invalid_argument("Duplicate forward no-flip helicity row"); }
    rows[incoming] = i;
  }
  if (std::find(rows.begin(), rows.end(), missing) != rows.end()) {
    throw std::invalid_argument("Missing forward no-flip helicity row");
  }
  return rows;
}

// Select every row in the explicit beam-leg helicity basis
// Full forward bases retain both no-flip and spin-flip transitions
//
std::vector<std::size_t> AllHelicityRows(const HELMatrix& hel) {
  std::vector<std::size_t> rows;
  rows.reserve(hel.lambda_values.size_row());
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) { rows.push_back(i); }
  return rows;
}

}  // namespace

// Contract a pair of continuum subchannel vertices into one final-state helicity matrix
// Internal exchange helicities are summed with lambda_ex_up + lambda_ex_dn = 0
//
// Compute one projection index or the basis dimension when it is absent
std::size_t ProjectionIndex(const std::vector<double>& basis, double projection) {
  const auto found = std::find_if(basis.cbegin(), basis.cend(), [projection](double value) {
    return std::abs(value - projection) < kIntegerTolerance;
  });
  return static_cast<std::size_t>(found - basis.cbegin());
}

// Compute the internal helicity contraction weight
std::complex<double> InternalWeight(const InternalHelicityMetric* metric, double upper, double lower) {
  if (metric == nullptr) { return std::abs(upper + lower) < kIntegerTolerance ? 1.0 : 0.0; }
  const std::size_t i = ProjectionIndex(metric->projections, upper);
  const std::size_t j = ProjectionIndex(metric->projections, -lower);
  if (i == metric->projections.size() || j == metric->projections.size()) { return 0.0; }
  return metric->matrix[i][j];
}

// Traverse the physical row contractions shared by full and projected operators
template <typename Action>
void ForEachSubchannelTerm(const HELMatrix& upper, const HELMatrix& lower, const InternalHelicityMetric* metric,
                           const MParticle& left_particle, const MParticle& right_particle, bool swap_final_order,
                           bool skip_unphysical, Action&& action) {
  if (metric != nullptr && (metric->matrix.size_row() != metric->projections.size() ||
                            metric->matrix.size_col() != metric->projections.size())) {
    throw std::invalid_argument("MHelicityScatter::Subchannel: invalid internal helicity metric");
  }
  const auto left  = FinalStateHelicities(left_particle, "MHelicityScatter::Subchannel left helicities");
  const auto right = FinalStateHelicities(right_particle, "MHelicityScatter::Subchannel right helicities");
  for (std::size_t i = 0; i < upper.lambda_values.size_row(); ++i) {
    for (std::size_t j = 0; j < lower.lambda_values.size_row(); ++j) {
      const auto weight = InternalWeight(metric, upper.lambda_values[i][1], lower.lambda_values[j][1]);
      if (std::fpclassify(std::abs(weight)) == FP_ZERO) { continue; }
      const double      lambda_left  = swap_final_order ? lower.lambda_values[j][0] : upper.lambda_values[i][0];
      const double      lambda_right = swap_final_order ? upper.lambda_values[i][0] : lower.lambda_values[j][0];
      const std::size_t left_index   = ProjectionIndex(left, lambda_left);
      const std::size_t right_index  = ProjectionIndex(right, lambda_right);
      if (left_index == left.size() || right_index == right.size()) {
        if (skip_unphysical) { continue; }
        throw std::invalid_argument("MHelicityScatter::Subchannel: final helicity is not physical");
      }
      action(i, j, left_index * right.size() + right_index, weight);
    }
  }
}

// Sew helicity rows of either full or projected pole frames
HelAmp SewSubchannel(const HELMatrix& upper, const HELMatrix& lower, const HelAmp& upper_frame,
                     const HelAmp& lower_frame, const MParticle& left, const MParticle& right, bool swap,
                     bool skip_unphysical, const InternalHelicityMetric* metric) {
  if (upper_frame.size_row() != upper.lambda_values.size_row() ||
      lower_frame.size_row() != lower.lambda_values.size_row()) {
    throw std::invalid_argument("Subchannel: helicity metadata dimensions disagree");
  }
  const auto n_final =
      FinalStateHelicityCount(left, "Subchannel left") * FinalStateHelicityCount(right, "Subchannel right");
  HelAmp out(upper_frame.size_col() * lower_frame.size_col(), n_final, 0.0);
  ForEachSubchannelTerm(upper, lower, metric, left, right, swap, skip_unphysical,
                        [&](std::size_t i, std::size_t j, std::size_t final, std::complex<double> weight) {
                          out.AddColumnKroneckerProduct(final, upper_frame.Row(i), lower_frame.Row(j), weight);
                        });
  return out;
}

// Sew two evaluated vertices with an optional internal pole numerator
HelAmp Subchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower, const MParticle& left,
                  const MParticle& right, bool swap, bool skip_unphysical, const InternalHelicityMetric* metric) {
  return SewSubchannel(upper.helicity, lower.helicity, upper.frame, lower.frame, left, right, swap, skip_unphysical,
                       metric);
}

// Project exchange helicities before sewing the internal pole
HelVec ProjectedSubchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower,
                           const HelVec& upper_source, const HelVec& lower_source, const MParticle& left,
                           const MParticle& right, bool swap, bool skip_unphysical,
                           const InternalHelicityMetric* metric) {
  const auto up = upper.frame * upper_source;
  const auto dn = lower.frame * lower_source;
  if (up.size() != upper.helicity.lambda_values.size_row() || dn.size() != lower.helicity.lambda_values.size_row()) {
    throw std::invalid_argument("ProjectedSubchannel: helicity metadata dimensions disagree");
  }
  HelVec out(FinalStateHelicityCount(left, "Subchannel left") * FinalStateHelicityCount(right, "Subchannel right"),
             0.0);
  ForEachSubchannelTerm(upper.helicity, lower.helicity, metric, left, right, swap, skip_unphysical,
                        [&](std::size_t i, std::size_t j, std::size_t final, std::complex<double> weight) {
                          out[final] += weight * up[i] * dn[j];
                        });
  return out;
}

// Project forward sources before sewing into canonical proton rows
HelAmp ProjectedSubchannel(const EvaluatedPoleSubvertex& upper, const EvaluatedPoleSubvertex& lower,
                           const HelAmp& upper_source, const HelAmp& lower_source, const MParticle& left,
                           const MParticle& right, bool swap, bool skip_unphysical,
                           const InternalHelicityMetric* metric) {
  const auto destination =
      CanonicalProtonPairSpinLayout::KroneckerDestinationRows(upper_source.size_row(), lower_source.size_row());
  return SewSubchannel(upper.helicity, lower.helicity, upper.frame * upper_source.Transpose(),
                       lower.frame * lower_source.Transpose(), left, right, swap, skip_unphysical, metric)
      .PermuteRows(destination);
}

// Build an incoherent spin-independent production operator
HelAmp Blind(std::size_t initial_states, std::size_t spin_states, std::complex<double> scale) {
  if (initial_states == 0 || spin_states == 0) { return HelAmp(); }

  HelAmp                     out(initial_states * spin_states, spin_states, 0.0);
  const std::complex<double> weight = scale / math::msqrt(static_cast<double>(spin_states));
  for (std::size_t i = 0; i < initial_states; ++i) {
    for (std::size_t h = 0; h < spin_states; ++h) { out[i * spin_states + h][h] = weight; }
  }
  return out;
}

// Build a spin-independent operator with exact forward transition sections
HelAmp Blind(const HelVec& upper_transition, const HelVec& lower_transition, std::size_t spin_states,
             std::complex<double> scale) {
  const auto destination_rows =
      CanonicalProtonPairSpinLayout::KroneckerDestinationRows(upper_transition.size(), lower_transition.size());
  const HelVec pair_transition = MappedKroneckerProduct(upper_transition, lower_transition, destination_rows);
  HelAmp       out             = Blind(pair_transition.size(), spin_states, scale);
  for (const auto& row : indices(pair_transition)) {
    for (std::size_t h = 0; h < spin_states; ++h) { out.ScaleRow(row * spin_states + h, pair_transition[row]); }
  }
  return out;
}

// Select physical rows from one forward helicity source
std::vector<std::size_t> Rows(const MDecayBranch& branch, bool drop_flip) {
  return drop_flip ? HelicityConservingRows(branch.hel) : AllHelicityRows(branch.hel);
}

// Build selected rows of one forward photon or exchanged-state source
HelAmp Forward(const LORENTZSCALAR& lts, const MDecayBranch& branch, const M4Vec& incoming, const M4Vec& outgoing,
               bool second, const std::vector<std::size_t>& rows, const std::string& photon_vertex,
               ForwardSpec forward) {
  if (IsPhotonExchangeBranch(branch)) {
    return qed::PhotonSourceMatrixTransitions(lts, second ? 2 : 1, PhotonTransitionsFromRows(branch.hel, rows), {-1, 1},
                                              photon_vertex, second);
  }
  const auto factors = ForwardExchangeResidueFactors(branch.hel, incoming, outgoing, second, rows, forward);
  return OuterProduct(factors.beam_transition, factors.exchange_helicity);
}

// Compute exact fixed-spin factors or no value for photon sources
std::optional<ForwardSourceFactors> ForwardFactors(const MDecayBranch& branch, const M4Vec& incoming,
                                                   const M4Vec& outgoing, bool second,
                                                   const std::vector<std::size_t>& rows, ForwardSpec forward) {
  if (IsPhotonExchangeBranch(branch)) { return std::nullopt; }
  return ForwardExchangeResidueFactors(branch.hel, incoming, outgoing, second, rows, forward);
}

// Multiply the physical production tensor by coherent central-spin weights
// The setup-time vertex reference and density matrix fix the scale
//
HelAmp Steer(const HelAmp& central, const HelVec& spin_steering) {
  if (central.size_col() != spin_steering.size()) {
    throw std::invalid_argument(
        "MHelicityScatter::Steer: central matrix and spin steering dimensions "
        "differ");
  }

  const double steering_norm2 = SquaredNorm(spin_steering);
  if (!(steering_norm2 > 0.0) || !std::isfinite(steering_norm2)) {
    throw std::invalid_argument("MHelicityScatter::Steer: spin steering has zero or non-finite norm");
  }

  return central.RightDiagonalProduct(spin_steering);
}

// Matrix form C = c kron(U,D) X with canonical proton-pair hard rows
// Component form C_(i1,i2,f1,f2,x) = c sum_(mn)
// U_(i1,f1,m) D_(i2,f2,n) X_(m,n,x)
// Raw U and D rows and persistent pair bits use the canonical negative-first
// order
//
HelAmp Contract(const HelAmp& upper, const HelAmp& lower, const HelAmp& central, std::complex<double> scale) {
  const std::size_t central_rows = upper.size_col() * lower.size_col();
  if (central.size_row() != central_rows) {
    throw std::invalid_argument(
        "MHelicityScatter::Contract: central "
        "basis dimension mismatch");
  }

  const auto destination_rows =
      CanonicalProtonPairSpinLayout::KroneckerDestinationRows(upper.size_row(), lower.size_row());
  HelAmp out = upper.KroneckerMultiply(lower, central, destination_rows);
  out *= scale;
  return out;
}

// Contract C = c (f_u kron f_d) [(r_u kron r_d)^T X] exactly once
// Fixed-spin beam transitions factor from exchange helicity vertices
HelAmp Contract(const ForwardSourceFactors& upper, const ForwardSourceFactors& lower, const HelAmp& central,
                std::complex<double> scale) {
  const std::size_t central_rows = upper.exchange_helicity.size() * lower.exchange_helicity.size();
  if (central.size_row() != central_rows) {
    throw std::invalid_argument(
        "MHelicityScatter::Contract: central "
        "basis dimension mismatch");
  }

  const auto destination_rows = CanonicalProtonPairSpinLayout::KroneckerDestinationRows(upper.beam_transition.size(),
                                                                                        lower.beam_transition.size());
  HelAmp     out              = central.FactorizedKroneckerMultiply(upper.beam_transition, upper.exchange_helicity,
                                                                    lower.beam_transition, lower.exchange_helicity, destination_rows);
  out *= scale;
  return out;
}

// Expand projected subchannels into canonical proton-pair hard rows
HelPair Contract(const ForwardSourceFactors& upper, const ForwardSourceFactors& lower, const HelVec& sub_t,
                 const HelVec& sub_u) {
  if (sub_t.size() != sub_u.size()) {
    throw std::invalid_argument(
        "MHelicityScatter::Contract: projected continuum "
        "dimensions disagree");
  }
  const auto   destination_rows = CanonicalProtonPairSpinLayout::KroneckerDestinationRows(upper.beam_transition.size(),
                                                                                          lower.beam_transition.size());
  const HelVec beam_transition = MappedKroneckerProduct(upper.beam_transition, lower.beam_transition, destination_rows);
  return {OuterProduct(beam_transition, sub_t), OuterProduct(beam_transition, sub_u)};
}

// Build the spherical propagator metric in one retained helicity basis
HelAmp ExchangeMetric(ExchangeBasisType basis, const std::vector<double>& projections, double pole_spin) {
  if (projections.empty()) { throw std::invalid_argument("ExchangeMetric: empty helicity basis"); }
  if (!std::isfinite(pole_spin)) { throw AmplitudeFailure("ExchangeMetric: non-finite generated spin"); }
  if (basis != ExchangeBasisType::HelicityTransport) {
    if (projections.size() != 1 || std::abs(projections[0]) > 1.0e-9) {
      throw std::invalid_argument("ExchangeMetric: single-residue basis must contain only zero");
    }
    return HelAmp(1, 1, 1.0);
  }
  HelAmp metric(projections.size(), projections.size(), 0.0);
  for (const auto& row : indices(projections)) {
    const int  m        = IntegerProjection(projections[row], "ExchangeMetric");
    const auto opposite = std::find_if(projections.cbegin(), projections.cend(),
                                       [m](double value) { return std::abs(value + static_cast<double>(m)) < 1.0e-9; });
    if (opposite == projections.cend()) { throw std::invalid_argument("ExchangeMetric: incomplete helicity basis"); }
    const std::size_t column = static_cast<std::size_t>(opposite - projections.cbegin());
    metric(row, column)      = spin::JacobWickSecondLegReversalPhase(pole_spin, m);
  }
  return metric;
}

// Build one boundary exchange-helicity vertex without proton spin factors
HelVec ExchangeBoundary(const std::vector<double>& projections, const M4Vec& transfer, bool second_exchange_daughter,
                        ForwardSpec forward) {
  if (forward.basis != ExchangeBasisType::HelicityTransport) {
    if (projections.size() != 1 || std::abs(projections[0]) > 1.0e-9) {
      throw std::invalid_argument("ExchangeBoundary: single-residue basis must contain only zero");
    }
    return {1.0};
  }
  return spin::ExchangeHelicityResidue(projections, transfer, second_exchange_daughter, forward);
}

}  // namespace spin
}  // namespace gra
