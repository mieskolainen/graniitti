#ifndef GRANIITTI_AMPLITUDE_MG5_YY_JJ_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_YY_JJ_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_bbx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_ddx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_uux.h"

namespace MG5_YY_JJ {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_bbx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_ddx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_JJ_P1_Sigma_sm_lepton_masses_aa_uux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
}


// Compute generated subprocess channels in BuildProcesses order
inline std::vector<std::vector<gra::mg5::Channel>> SubprocessChannels() {
  using gra::AmplitudeTopologyNode;
  using gra::mg5helas::ColorFlowLeg;
  using gra::mg5::Channel;
  return {
      {Channel{{22, 22}, {5, -5}, {AmplitudeTopologyNode{{5}, {}}, AmplitudeTopologyNode{{-5}, {}}}, {1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {1, -1}, {AmplitudeTopologyNode{{1}, {}}, AmplitudeTopologyNode{{-1}, {}}}, {1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}, Channel{{22, 22}, {3, -3}, {AmplitudeTopologyNode{{3}, {}}, AmplitudeTopologyNode{{-3}, {}}}, {1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {2, -2}, {AmplitudeTopologyNode{{2}, {}}, AmplitudeTopologyNode{{-2}, {}}}, {1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}, Channel{{22, 22}, {4, -4}, {AmplitudeTopologyNode{{4}, {}}, AmplitudeTopologyNode{{-4}, {}}}, {1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}}
  };
}

// Build generated subprocesses with their exact channels
inline void BuildSubprocesses(std::vector<gra::mg5::Subprocess<ProcessBase>> &out,
                              const std::string &param_card) {
  std::vector<std::unique_ptr<ProcessBase>> matrix_elements;
  BuildProcesses(matrix_elements, param_card);
  auto channels = SubprocessChannels();
  if (matrix_elements.size() != channels.size()) {
    throw std::invalid_argument(
        "MG5_YY_JJ::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_YY_JJ

#endif
