// Generated MG5 amplitude for gamma-gamma WW production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/Processes.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {

// Construct all generated subprocesses
AMP_MG5_yy_ww::AMP_MG5_yy_ww()
    : PhotonMG5Process(amplitude::Processes("MG5_YY_WW"), 0.0) {
  std::vector<mg5::Subprocess<MG5_YY_WW::ProcessBase>> subprocesses;
  MG5_YY_WW::BuildSubprocesses(
      subprocesses,
      aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_WW/param_card.dat"));
  mg5::ValidateSubprocessInitialStates(
      subprocesses, mg5::PhotonInitialStates(), "AMP_MG5_yy_ww");
  subprocess_sum.SetSubprocesses(std::move(subprocesses));
}

// Compute the Born matrix-element power of alpha_s
int AMP_MG5_yy_ww::AlphaSPower() const noexcept { return subprocess_sum.AlphaSPower(); }

// Compute the Born matrix-element power of alpha_QED
int AMP_MG5_yy_ww::AlphaQEDPower() const noexcept { return subprocess_sum.AlphaQEDPower(); }

// Evaluate the summed matrix element squared
mg5helas::MatrixElementEvaluation AMP_MG5_yy_ww::Evaluate(
    LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) {
  has_final_state_color_ = false;
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }
  const bool epa_hard = coherent_epa &&
                        mg5helas::HasEPAHardBeamState(lts);
  const auto status = subprocess_sum.Prepare(lts, alpha_s, epa_hard);
  if (!mg5helas::EvaluationSucceeded(status)) {
    return {status, 0.0};
  }
  const auto *channel = subprocess_sum.MatchingChannel(
      lts, {PDG::PDG_gamma, PDG::PDG_gamma});
  if (channel == nullptr) {
    return {mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};
  }
  has_final_state_color_ = mg5::ChannelHasFinalStateColor(*channel);
  double amp2 = 0.0;
  const auto amplitude_status =
      subprocess_sum.CalcPreparedPhotonHelicityAmp2(lts, coherent_epa, amp2);
  if (!mg5helas::EvaluationSucceeded(amplitude_status)) {
    has_final_state_color_ = false;
  }
  return {amplitude_status, amp2};
}

// Sample one generated final-state color flow
bool AMP_MG5_yy_ww::SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) {
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {
    return false;
  }
  if (!has_final_state_color_) {
    return ClearHardColorFlow(lts);
  }
  const auto *channel = subprocess_sum.MatchingChannel(
      lts, {PDG::PDG_gamma, PDG::PDG_gamma});
  if (channel == nullptr) {
    return false;
  }
  return SampleHardColorFlow(
      lts, random, lts.hard_color_flows,
      channel->external_color_flows.size());
}

// Compute whether the selected event channel has colored stable final states
bool AMP_MG5_yy_ww::HasFinalStateColor() const {
  return has_final_state_color_;
}

// Compute the number of generated subprocesses
std::size_t AMP_MG5_yy_ww::SubprocessCount() const {
  return subprocess_sum.SubprocessCount();
}

}  // namespace gra
