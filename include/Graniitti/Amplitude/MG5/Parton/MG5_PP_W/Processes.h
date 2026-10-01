#ifndef GRANIITTI_AMPLITUDE_MG5_PP_W_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_PP_W_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P1_Sigma_sm_lepton_masses_udx_epve.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P2_Sigma_sm_lepton_masses_dux_emvex.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P3_Sigma_sm_lepton_masses_udx_mupvm.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P4_Sigma_sm_lepton_masses_dux_mumvmx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P5_Sigma_sm_lepton_masses_udx_tapvt.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/MG5_PP_W_P6_Sigma_sm_lepton_masses_dux_tamvtx.h"

namespace MG5_PP_W {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_PP_W_P1_Sigma_sm_lepton_masses_udx_epve>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_W_P2_Sigma_sm_lepton_masses_dux_emvex>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_W_P3_Sigma_sm_lepton_masses_udx_mupvm>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_W_P4_Sigma_sm_lepton_masses_dux_mumvmx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_W_P5_Sigma_sm_lepton_masses_udx_tapvt>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_W_P6_Sigma_sm_lepton_masses_dux_tamvtx>();
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
      {Channel{{2, -1}, {-11, 12}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 2}, {-11, 12}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -3}, {-11, 12}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 4}, {-11, 12}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -2}, {11, -12}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 1}, {11, -12}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -4}, {11, -12}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 3}, {11, -12}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -1}, {-13, 14}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 2}, {-13, 14}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -3}, {-13, 14}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 4}, {-13, 14}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -2}, {13, -14}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 1}, {13, -14}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -4}, {13, -14}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 3}, {13, -14}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -1}, {-15, 16}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 2}, {-15, 16}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -3}, {-15, 16}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 4}, {-15, 16}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -2}, {15, -16}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 1}, {15, -16}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -4}, {15, -16}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 3}, {15, -16}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}}
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
        "MG5_PP_W::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_PP_W

#endif
