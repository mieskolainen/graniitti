// Generated event invariants and forward-leg kinematics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Kinematics/MEventKinematics.h"

// C++
#include <array>
#include <cmath>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MSetup.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::msqrt;

namespace gra {

namespace {

// Compute a local roundoff bound for one reconstructed lightlike transfer
double LightlikeTransferRoundoffTolerance(const gra::M4Vec &transfer, const gra::M4Vec &beam,
                                          const gra::M4Vec &forward) {
  const double invariant_scale   = transfer.E() * transfer.E() + transfer.P3mod2();
  const double subtraction_scale = std::abs(transfer.E()) * (std::abs(beam.E()) + std::abs(forward.E())) +
                                   std::abs(transfer.Px()) * (std::abs(beam.Px()) + std::abs(forward.Px())) +
                                   std::abs(transfer.Py()) * (std::abs(beam.Py()) + std::abs(forward.Py())) +
                                   std::abs(transfer.Pz()) * (std::abs(beam.Pz()) + std::abs(forward.Pz()));
  return 128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, invariant_scale + subtraction_scale);
}

// Canonicalize one transfer known to be exactly lightlike
bool CanonicalizeLightlikeTransfer(gra::M4Vec &transfer, double &virtuality, const gra::M4Vec &beam,
                                   const gra::M4Vec &forward) {
  const double tolerance = LightlikeTransferRoundoffTolerance(transfer, beam, forward);
  if (!(transfer.E() > 0.0) || !(transfer.P3mod2() > 0.0) || std::abs(virtuality) > tolerance) { return false; }
  transfer.SetE(transfer.P3mod());
  virtuality = 0.0;
  return true;
}

}  // namespace

// Resolve and validate one generated forward leg from event kinematics
ForwardLegState ResolveForwardLegState(const LORENTZSCALAR &lts, const ForwardBeamLeg leg) {
  if (leg != ForwardBeamLeg::Upper && leg != ForwardBeamLeg::Lower) {
    throw std::invalid_argument("ResolveForwardLegState: invalid generated beam leg");
  }
  const bool        upper       = leg == ForwardBeamLeg::Upper;
  const std::size_t final_index = upper ? 1 : 2;
  if (lts.pfinal.size() <= final_index) {
    throw AmplitudeFailure("ResolveForwardLegState: generated forward system is absent");
  }

  ForwardLegState state;
  state.leg = leg;
  state.final_state =
      (upper ? lts.excite1 : lts.excite2) ? ForwardFinalState::InclusiveExcitation : ForwardFinalState::Elastic;
  state.has_xi  = upper ? lts.has_xi1 : lts.has_xi2;
  state.xi      = upper ? lts.xi1 : lts.xi2;
  state.t       = upper ? lts.t1 : lts.t2;
  state.qt      = upper ? lts.qt1 : lts.qt2;
  state.mass2   = lts.pfinal[final_index].M2();
  state.emitter = upper ? lts.beam1 : lts.beam2;
  if (lts.model_cache != nullptr) { state.structure = lts.model_cache->Tune().Structure(); }
  state.is_nuclear = nuclear::IsNuclearPDG(state.emitter.pdg);
  const auto upc   = lts.upc_event != nullptr ? lts.upc_event : lts.upc_model;
  if (upc != nullptr) {
    const std::size_t index = upper ? 0 : 1;
    state.emission          = upc->Param().emission[index];
    state.target            = upc->Param().target[index];
    state.neutron           = upc->Param().neutron[index];
  }
  state.upc_model = upc;
  state.incoming  = upper ? lts.pbeam1 : lts.pbeam2;
  state.outgoing  = lts.pfinal[final_index];
  state.transfer  = upper ? lts.q1 : lts.q2;

  // Remove only independent EMD excitation, preserving hard knockout invariant masses
  const double emd = lts.forward_emd[upper ? 0 : 1];
  if (state.is_nuclear && emd > 0.0) {
    const auto& opposite = upper ? lts.pbeam2 : lts.pbeam1;
    const double mass2 = math::pow2(state.emitter.mass), product = state.incoming.DotM(opposite);
    const double other_mass2 = math::pow2(upper ? lts.beam2.mass : lts.beam1.mass);
    const auto lightlike = opposite - state.incoming * (other_mass2 / (product + std::sqrt(product * product - mass2 * other_mass2)));
    const double hard_mass2 = math::pow2(std::sqrt(state.mass2) - emd);
    state.transfer += lightlike * ((state.mass2 - hard_mass2) / (2.0 * state.outgoing.DotM(lightlike)));
    state.outgoing = state.incoming - state.transfer;
    state.mass2 = hard_mass2;
    state.t = state.transfer.M2();
  }

  if (!std::isfinite(state.t) || state.t > 1.0e-12 || !std::isfinite(state.qt) || state.qt < 0.0 ||
      !std::isfinite(state.mass2) || state.mass2 <= 0.0 || !std::isfinite(state.transfer.M2())) {
    throw AmplitudeFailure("ResolveForwardLegState: invalid generated forward kinematics");
  }
  if (state.has_xi && !std::isfinite(state.xi)) {
    throw AmplitudeFailure("ResolveForwardLegState: invalid forward momentum loss");
  }
  return state;
}

// Build and check Lorentz scalars
// Input (Nf) is the number of final states
//
bool kinematics::SetLorentzScalars(MProcessState &state, unsigned int Nf, bool loop_active,
                                   TransferVirtualityPolicy transfer_policy) {
  // Example Nf = 4 gives 8 Lorentz scalars:
  //
  //        t1
  // ---------------->                 state.lts.pfinal[1]
  //         |   s13
  //         |------------->           state.lts.pfinal[3]
  //         |
  // s   that,uhat       shat (s34)    state.lts.pfinal[0]
  //         |
  //         |------------->           state.lts.pfinal[4]
  //         |   s24
  // ---------------->                 state.lts.pfinal[2]
  //        t2

  // ------------------------------------------------------------------
  // s-type Lorentz scalars -->

  const int         offset          = 3;  // central system indexing starts from offset
  const std::size_t scalar_capacity = sizeof(state.lts.ss) / sizeof(state.lts.ss[0]);
  if (Nf >= scalar_capacity || Nf >= state.lts.pfinal.size()) {
    throw std::invalid_argument(
        "SetLorentzScalars: final-state "
        "multiplicity exceeds scalar storage");
  }

  // Reject non-finite generated momenta before constructing invariants and rest frames
  if (!gra::AllFinite(state.lts.pbeam1.Contravariant()) || !gra::AllFinite(state.lts.pbeam2.Contravariant())) {
    return false;
  }
  for (std::size_t i = 0; i <= Nf; ++i) {
    if (!gra::AllFinite(state.lts.pfinal[i].Contravariant())) { return false; }
  }

  // The accepted Born point defines the masses preserved by screening shifts
  if (!loop_active) { state.lts.forward_mass2 = {state.lts.pfinal[1].M2(), state.lts.pfinal[2].M2()}; }

  // Diagonal (Note! do not remove, this is used by some processes)
  //
  for (std::size_t i = 0; i <= Nf; ++i) {
    state.lts.ss[i][i] = state.lts.pfinal[i].M2();
    if (!std::isfinite(state.lts.ss[i][i])) { return false; }
  }

  // Compute upper triangle values
  // [ -, o, o, o ]
  // [  , -, o, o ]
  // [  ,  , -, o ]
  // [  ,  ,  , - ]
  //
  for (std::size_t i = 0; i <= Nf; ++i) {
    for (std::size_t j = i + 1; j <= Nf; ++j) {
      state.lts.ss[i][j] = (state.lts.pfinal[i] + state.lts.pfinal[j]).M2();
      if (!std::isfinite(state.lts.ss[i][j]) || state.lts.ss[i][j] < 0) { return false; }
    }
  }
  // Copy the i<->j permutated to the lower triangle
  // (faster than calculating twice)
  for (std::size_t i = 1; i <= Nf; ++i) {
    for (std::size_t j = 0; j < i; ++j) { state.lts.ss[i][j] = state.lts.ss[j][i]; }
  }

  // DEBUG PRINT
  /*
  for (std::size_t i = 0; i <= Nf; ++i) {  // start at 0
    for (std::size_t j = 0; j <= Nf; ++j) {
      printf("%0.3f (%d,%d) ", state.lts.ss[i][j], i, j);
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
  */

  // Sub invariants w.r.t central system
  state.lts.s1 = state.lts.ss[1][0];
  state.lts.s2 = state.lts.ss[2][0];
  if (state.lts.s1 < 0 || state.lts.s2 < 0) { return false; }

  // ------------------------------------------------------------------
  // Propagator vectors

  std::array<M4Vec, 2> transfer = {state.lts.pbeam1 - state.lts.pfinal[1], state.lts.pbeam2 - state.lts.pfinal[2]};
  
  if (transfer_policy == TransferVirtualityPolicy::CollinearLightlike) {

    // Preserve soft collinear photons without subtracting beam and remnant momenta
    const double e1 = 0.5 * state.lts.pfinal[0].LightconePos();
    const double e2 = 0.5 * state.lts.pfinal[0].LightconeNeg();
    transfer = {M4Vec(0.0, 0.0, e1, e1), M4Vec(0.0, 0.0, -e2, e2)};
  }
  if (state.lts.upc_model != nullptr &&
      !nuclear::CloseTransfers({state.lts.beam1.pdg, state.lts.beam2.pdg}, {state.lts.pbeam1, state.lts.pbeam2},
                               {state.lts.pfinal[1], state.lts.pfinal[2]}, state.lts.pfinal[0], transfer)) {
    return false;
  }
  M4Vec q1 = transfer[0];
  M4Vec q2 = transfer[1];

  // ------------------------------------------------------------------
  // t-type Lorentz scalars -->

  double t1 = q1.M2();
  double t2 = q2.M2();

  if (!std::isfinite(t1) || !std::isfinite(t2)) { return false; }

  // Apply the configured virtuality condition independently to both legs
  const bool lightlike1 = transfer_policy == TransferVirtualityPolicy::CollinearLightlike ||
                          (transfer_policy == TransferVirtualityPolicy::HardDiffractive && !state.lts.hard_diff1);
  const bool lightlike2 = transfer_policy == TransferVirtualityPolicy::CollinearLightlike ||
                          (transfer_policy == TransferVirtualityPolicy::HardDiffractive && !state.lts.hard_diff2);
  const bool valid1 =
      lightlike1 ? (CanonicalizeLightlikeTransfer(q1, t1, state.lts.pbeam1, state.lts.pfinal[1]) ||
                    (transfer_policy == TransferVirtualityPolicy::HardDiffractive && t1 < 0.0)) : t1 <= 0.0;
  const bool valid2 =
      lightlike2 ? (CanonicalizeLightlikeTransfer(q2, t2, state.lts.pbeam2, state.lts.pfinal[2]) ||
                    (transfer_policy == TransferVirtualityPolicy::HardDiffractive && t2 < 0.0)) : t2 <= 0.0;
  if (!valid1 || !valid2) { return false; }
  state.lts.q1 = q1;
  state.lts.q2 = q2;
  state.lts.t1 = t1;
  state.lts.t2 = t2;

  // Boost propagators to the system X-rest frame after virtuality handling
  M4Vec        q1_in_X = state.lts.q1;
  M4Vec        q2_in_X = state.lts.q2;
  const double M0      = state.lts.pfinal[0].M();

  gra::kinematics::LorentzBoost(state.lts.pfinal[0], M0, q1_in_X, -1);
  gra::kinematics::LorentzBoost(state.lts.pfinal[0], M0, q2_in_X, -1);
  if (!kinematics::FiniteFourVector(q1_in_X) || !kinematics::FiniteFourVector(q2_in_X)) { return false; }
  state.lts.q1_in_X = q1_in_X;
  state.lts.q2_in_X = q2_in_X;

  // For 2-body central processes
  if (state.lts.decaytree.size() == 2) {
    state.lts.t_hat = (state.lts.q1 - state.lts.decaytree[0].p4).M2();  // note q1 on both
    state.lts.u_hat = (state.lts.q1 - state.lts.decaytree[1].p4).M2();  // in t and u!
  }

  // Compute (only) for higher multiplicity processes
  if (Nf > 4) {
    for (std::size_t i = offset; i <= Nf; ++i) { state.lts.tt_1[i] = (state.lts.q1 - state.lts.pfinal[i]).M2(); }
    for (std::size_t i = offset; i <= Nf; ++i) {
      for (std::size_t j = offset; j <= Nf; ++j) {
        state.lts.tt_xy[i][j] = (state.lts.q1 - state.lts.pfinal[i] - state.lts.pfinal[j]).M2();
      }
    }
    for (std::size_t i = offset; i <= Nf; ++i) { state.lts.tt_2[i] = (state.lts.q2 - state.lts.pfinal[i]).M2(); }
  }

  // Incoming exchange fractions follow from the complete forward-system recoil
  state.lts.x1 = kinematics::LongitudinalMomentumLoss(state.lts.pbeam1, state.lts.pfinal[1], true);
  state.lts.x2 = kinematics::LongitudinalMomentumLoss(state.lts.pbeam2, state.lts.pfinal[2], false);

  // Generic noncollinear processes expose the same losses as physical forward
  // xi
  state.lts.has_xi1 = state.phase_space_class != "P";
  state.lts.has_xi2 = state.phase_space_class != "P";
  state.lts.xi1     = state.lts.has_xi1 ? state.lts.x1 : 0.0;
  state.lts.xi2     = state.lts.has_xi2 ? state.lts.x2 : 0.0;

  // Bjorken-x [0,1] (this Lorentz invariant expression
  // gives 1 always for elastic central production forward leg)
  state.lts.xbj1 = state.lts.t1 / (2 * (state.lts.pbeam1 * state.lts.q1));
  state.lts.xbj2 = state.lts.t2 / (2 * (state.lts.pbeam2 * state.lts.q2));

  // Propagator pt
  state.lts.qt1 = state.lts.q1.Pt();
  state.lts.qt2 = state.lts.q2.Pt();

  // No need to update these when loop integrating
  // (central system momentum not impacted)
  if (loop_active == false) {
    // System variables
    state.lts.m2    = state.lts.pfinal[0].M2();
    state.lts.s_hat = state.lts.m2;
    state.lts.Y     = state.lts.pfinal[0].Rap();
    state.lts.Pt    = state.lts.pfinal[0].Pt();

    // For 2-body central processes
    if (state.lts.decaytree.size() == 2) {
      // Boost daughters to the system X-rest frame
      std::vector<M4Vec> p  = {state.lts.decaytree[0].p4, state.lts.decaytree[1].p4};
      const double       M0 = state.lts.pfinal[0].M();

      for (const auto &i : gra::aux::indices(p)) {
        gra::kinematics::LorentzBoost(state.lts.pfinal[0], M0, p[i],
                                      -1);  // Note the minus sign
      }
      state.lts.d0_in_X = p[0];
      state.lts.d1_in_X = p[1];
    }
  }

  return true;
}

// Reconstruct exact on-shell forward momenta at fixed central four-momentum
// p1T -> p1T-kT, p2T -> p2T+kT and p1+p2+pX = P1+P2
bool kinematics::RebuildScreeningKinematics(LORENTZSCALAR &lts, const std::array<double, 2> &upper_transverse,
                                            const std::array<double, 2> &lower_transverse,
                                            const bool                   central_system_active) {
  if (lts.pfinal.size() < 3 || lts.pfinal_orig.size() < 3) { return false; }
  const M4Vec  incoming       = lts.pbeam1 + lts.pbeam2;
  const double collision_mass = incoming.M();
  if (!std::isfinite(collision_mass) || !(collision_mass > 0.0)) { return false; }

  M4Vec central_cm = central_system_active ? lts.pfinal[0] : M4Vec();
  if (central_system_active) { LorentzBoost(incoming, collision_mass, central_cm, -1); }
  M4Vec upper_cm = lts.pfinal_orig[1];
  M4Vec lower_cm = lts.pfinal_orig[2];
  upper_cm.SetPxPy(upper_transverse[0], upper_transverse[1]);
  lower_cm.SetPxPy(lower_transverse[0], lower_transverse[1]);

  // Keep the sampled invariant masses of the two outgoing beam systems fixed at every kT node
  const double upper_mass2 = lts.forward_mass2[0];
  const double lower_mass2 = lts.forward_mass2[1];
  if (!std::isfinite(upper_mass2) || !std::isfinite(lower_mass2) || upper_mass2 < 0.0 || lower_mass2 < 0.0) {
    return false;
  }
  const double upper_pz = SolvePz(std::sqrt(upper_mass2), std::sqrt(lower_mass2), upper_cm.Pt(), lower_cm.Pt(),
                                  central_cm.Pz(), central_cm.E(), incoming.M2());
  const double lower_pz = -(central_cm.Pz() + upper_pz);
  if (!std::isfinite(upper_pz) || upper_pz < 0.0 || lower_pz > 0.0) { return false; }
  upper_cm.SetPzE(upper_pz, msqrt(upper_mass2 + math::pow2(upper_cm.Pt()) + math::pow2(upper_pz)));
  lower_cm.SetPzE(lower_pz, msqrt(lower_mass2 + math::pow2(lower_cm.Pt()) + math::pow2(lower_pz)));

  LorentzBoost(incoming, collision_mass, upper_cm, 1);
  LorentzBoost(incoming, collision_mass, lower_cm, 1);
  const M4Vec final_sum = upper_cm + lower_cm + (central_system_active ? lts.pfinal[0] : M4Vec());
  if (!math::CheckEMC(incoming - final_sum)) { return false; }
  lts.pfinal[1] = upper_cm;
  lts.pfinal[2] = lower_cm;
  return true;
}


}  // namespace gra
