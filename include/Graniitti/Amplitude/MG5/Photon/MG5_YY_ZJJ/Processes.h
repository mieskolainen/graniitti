#ifndef GRANIITTI_AMPLITUDE_MG5_YY_ZJJ_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_YY_ZJJ_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P10_Sigma_sm_lepton_masses_aa_epemccx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P11_Sigma_sm_lepton_masses_aa_mupmumccx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P12_Sigma_sm_lepton_masses_aa_taptamccx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P13_Sigma_sm_lepton_masses_aa_epembbx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P14_Sigma_sm_lepton_masses_aa_mupmumbbx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P15_Sigma_sm_lepton_masses_aa_taptambbx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P1_Sigma_sm_lepton_masses_aa_epemuux.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P2_Sigma_sm_lepton_masses_aa_mupmumuux.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P3_Sigma_sm_lepton_masses_aa_taptamuux.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P4_Sigma_sm_lepton_masses_aa_epemddx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P5_Sigma_sm_lepton_masses_aa_mupmumddx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P6_Sigma_sm_lepton_masses_aa_taptamddx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P7_Sigma_sm_lepton_masses_aa_epemssx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P8_Sigma_sm_lepton_masses_aa_mupmumssx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/MG5_YY_ZJJ_P9_Sigma_sm_lepton_masses_aa_taptamssx.h"

namespace MG5_YY_ZJJ {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P10_Sigma_sm_lepton_masses_aa_epemccx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P11_Sigma_sm_lepton_masses_aa_mupmumccx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P12_Sigma_sm_lepton_masses_aa_taptamccx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P13_Sigma_sm_lepton_masses_aa_epembbx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P14_Sigma_sm_lepton_masses_aa_mupmumbbx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P15_Sigma_sm_lepton_masses_aa_taptambbx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P1_Sigma_sm_lepton_masses_aa_epemuux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P2_Sigma_sm_lepton_masses_aa_mupmumuux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P3_Sigma_sm_lepton_masses_aa_taptamuux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P4_Sigma_sm_lepton_masses_aa_epemddx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P5_Sigma_sm_lepton_masses_aa_mupmumddx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P6_Sigma_sm_lepton_masses_aa_taptamddx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P7_Sigma_sm_lepton_masses_aa_epemssx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P8_Sigma_sm_lepton_masses_aa_mupmumssx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_ZJJ_P9_Sigma_sm_lepton_masses_aa_taptamssx>();
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
      {Channel{{22, 22}, {-11, 11, 4, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{4}, {}}, AmplitudeTopologyNode{{-4}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-13, 13, 4, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{4}, {}}, AmplitudeTopologyNode{{-4}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-15, 15, 4, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{4}, {}}, AmplitudeTopologyNode{{-4}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-11, 11, 5, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{5}, {}}, AmplitudeTopologyNode{{-5}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-13, 13, 5, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{5}, {}}, AmplitudeTopologyNode{{-5}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-15, 15, 5, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{5}, {}}, AmplitudeTopologyNode{{-5}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-11, 11, 2, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{2}, {}}, AmplitudeTopologyNode{{-2}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-13, 13, 2, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{2}, {}}, AmplitudeTopologyNode{{-2}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-15, 15, 2, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{2}, {}}, AmplitudeTopologyNode{{-2}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-11, 11, 1, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{1}, {}}, AmplitudeTopologyNode{{-1}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-13, 13, 1, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{1}, {}}, AmplitudeTopologyNode{{-1}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-15, 15, 1, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{1}, {}}, AmplitudeTopologyNode{{-1}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-11, 11, 3, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{3}, {}}, AmplitudeTopologyNode{{-3}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-13, 13, 3, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{3}, {}}, AmplitudeTopologyNode{{-3}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}},
      {Channel{{22, 22}, {-15, 15, 3, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{3}, {}}, AmplitudeTopologyNode{{-3}, {}}}, {1, 1, 1, 1, 3, -3}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}}}}}
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
        "MG5_YY_ZJJ::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_YY_ZJJ

#endif
