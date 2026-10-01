#ifndef GRANIITTI_AMPLITUDE_MG5_YY_WW_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_YY_WW_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P1_Sigma_sm_lepton_masses_aa_epveemvex.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P2_Sigma_sm_lepton_masses_aa_epvemumvmx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P3_Sigma_sm_lepton_masses_aa_epvetamvtx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P4_Sigma_sm_lepton_masses_aa_mupvmemvex.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P5_Sigma_sm_lepton_masses_aa_mupvmmumvmx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P6_Sigma_sm_lepton_masses_aa_mupvmtamvtx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P7_Sigma_sm_lepton_masses_aa_tapvtemvex.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P8_Sigma_sm_lepton_masses_aa_tapvtmumvmx.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/MG5_YY_WW_P9_Sigma_sm_lepton_masses_aa_tapvttamvtx.h"

namespace MG5_YY_WW {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_YY_WW_P1_Sigma_sm_lepton_masses_aa_epveemvex>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P2_Sigma_sm_lepton_masses_aa_epvemumvmx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P3_Sigma_sm_lepton_masses_aa_epvetamvtx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P4_Sigma_sm_lepton_masses_aa_mupvmemvex>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P5_Sigma_sm_lepton_masses_aa_mupvmmumvmx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P6_Sigma_sm_lepton_masses_aa_mupvmtamvtx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P7_Sigma_sm_lepton_masses_aa_tapvtemvex>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P8_Sigma_sm_lepton_masses_aa_tapvtmumvmx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_YY_WW_P9_Sigma_sm_lepton_masses_aa_tapvttamvtx>();
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
      {Channel{{22, 22}, {-11, 12, 11, -12}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-11, 12, 13, -14}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-11, 12, 15, -16}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{12}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-13, 14, 11, -12}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-13, 14, 13, -14}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-13, 14, 15, -16}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{14}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-15, 16, 11, -12}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{11}, {}}, AmplitudeTopologyNode{{-12}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-15, 16, 13, -14}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{13}, {}}, AmplitudeTopologyNode{{-14}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{22, 22}, {-15, 16, 15, -16}, {AmplitudeTopologyNode{{24}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{16}, {}}}}, AmplitudeTopologyNode{{-24}, {AmplitudeTopologyNode{{15}, {}}, AmplitudeTopologyNode{{-16}, {}}}}}, {1, 1, 1, 1, 1, 1}, {{ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}}
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
        "MG5_YY_WW::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_YY_WW

#endif
