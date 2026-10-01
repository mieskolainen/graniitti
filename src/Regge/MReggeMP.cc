// Minimal Pomeron central production amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeMP.h"

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::mpom {

// Evaluate the complete MP fusion tensor in the selected spin frame
HelAmp Fusion(const LORENTZSCALAR& lts, const spin::PoleLS& vertex) {
  const auto   helicity    = spin::PoleLSHelicity(vertex, lts.q1_in_X.P3mod());
  const double central_phi = lts.q1_in_X.Phi();
  const HelAmp central     = spin::fDecayMatrix(helicity, lts.q1_in_X.Theta(), central_phi);

  // Complete the CM Jacob-Wick section carried by the collider exchange vertices
  HelVec phase(central.size_row());
  for (const auto row : indices(phase)) {
    const double exchange_helicity = helicity.lambda_values[row][0] - helicity.lambda_values[row][1];
    phase[row]                     = std::exp(-math::zi * exchange_helicity * central_phi);
  }
  if (vertex.dual_spin.size() != central.size_col() || vertex.spin_metric.size() != central.size_col()) {
    throw std::invalid_argument("MReggeMP::Fusion: invalid produced-spin metric");
  }
  std::vector<std::size_t> dual;
  dual.reserve(vertex.dual_spin.size());
  for (const double projection : vertex.dual_spin) {
    dual.push_back(spin::SpinProjectionIndex(projection, helicity.J, "MReggeMP::Fusion produced spin"));
  }
  // C_CM = E F G has exchange-helicity rows and resonance-spin columns, with exchange phase E
  // Decay-like F carries m, paired with M by the spherical metric G_{m M} = (-1)^(J-M) delta_{m,-M}
  // These dual projections m = -M couple fusion and decay spin-J indices to total spin zero
  // Wigner D = D(R) rotates magnetic projections, leaving reduced LS couplings invariant
  // In a common external helicity basis, C_R = C_CM D^T and B_R = D^* B_CM
  // T denotes transpose, * conjugation. Unitarity gives (D^dagger D)^* = D^T D^* = I
  // Hence the fusion amplitude C_R B_R = C_CM B_CM is independent of the spin axes
  const auto exchange = HelAmp::DiagonalMatrix(phase);
  const auto metric =
      HelAmp::IdentityMatrix(central.size_col()).SelectColumns(dual).RightDiagonalProduct(vertex.spin_metric);
  const auto resonance = metric * spin::ProductionRotation(lts, lts.process.MP_FRAME, helicity.J).Transpose();

  // C_R = E F G D(R)^T, with exchange = E, central = F, resonance = G D(R)^T
  return exchange * central * resonance;
}

// Build the MP production matrices with C_filtered = C S^T
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, const double s0,
                              const std::array<ForwardLegState, 2>* forward_state) {
  const auto* filter = res.UsesUnrestrictedSpinBasis() ? nullptr : &res.MP.filter;
  return rspin::Resonance(lts, res, Fusion, s0, filter, forward_state);
}

// Build the MP crossed continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, const double s0) { return rspin::Continuum(lts, s0); }

}  // namespace gra::mpom
