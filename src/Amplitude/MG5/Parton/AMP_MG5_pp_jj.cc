// Generated MG5 amplitude for parton jj production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_jj.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_JJ/Processes.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {

// Construct all generated subprocesses
AMP_MG5_pp_jj::AMP_MG5_pp_jj()
    : PartonMG5Process(amplitude::Processes("MG5_PP_JJ")) {
  std::vector<mg5::Subprocess<MG5_PP_JJ::ProcessBase>> subprocesses;
  MG5_PP_JJ::BuildSubprocesses(
      subprocesses,
      aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_JJ/param_card.dat"));
  mg5::ValidateSubprocessInitialStates(
      subprocesses, mg5::MasslessQCDInitialStates(), "AMP_MG5_pp_jj");
  subprocess_sum.SetSubprocesses(std::move(subprocesses));
}

// Prepare event-local momenta and final-state channels
mg5helas::EvaluationStatus AMP_MG5_pp_jj::Prepare(LORENTZSCALAR &lts,
                                                double alpha_s) {
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {
    return mg5helas::EvaluationStatus::AmplitudeFailure;
  }
  return subprocess_sum.Prepare(lts, alpha_s);
}

// Evaluate one prepared event and return all event-local results
PartonMG5Evaluation AMP_MG5_pp_jj::EvaluatePrepared(LORENTZSCALAR &lts,
                                                double alpha_s) {
  PartonMG5Evaluation result;
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {
    result.status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return result;
  }
  if (!subprocess_sum.PreparedStateMatches(lts, alpha_s)) {
    result.status = Prepare(lts, alpha_s);
  }
  if (!result.Valid()) {
    return result;
  }
  result.status = subprocess_sum.CalcPreparedHelicityAmp2(lts, result.amp2);
  result.contributing_subprocesses =
      subprocess_sum.ContributingSubprocesses();
  result.color_flows = subprocess_sum.LastHardColorFlows();
  return result;
}

// Compute the number of generated subprocesses
std::size_t AMP_MG5_pp_jj::SubprocessCount() const {
  return subprocess_sum.SubprocessCount();
}

}  // namespace gra
