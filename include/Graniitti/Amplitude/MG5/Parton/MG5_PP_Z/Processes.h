#ifndef GRANIITTI_AMPLITUDE_MG5_PP_Z_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_PP_Z_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P1_Sigma_sm_lepton_masses_bbx_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P1_Sigma_sm_lepton_masses_ddx_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P1_Sigma_sm_lepton_masses_uux_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P2_Sigma_sm_lepton_masses_bbx_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P2_Sigma_sm_lepton_masses_ddx_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P2_Sigma_sm_lepton_masses_uux_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P3_Sigma_sm_lepton_masses_bbx_taptam.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P3_Sigma_sm_lepton_masses_ddx_taptam.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P3_Sigma_sm_lepton_masses_uux_taptam.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P4_Sigma_sm_lepton_masses_bbx_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P4_Sigma_sm_lepton_masses_ddx_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P4_Sigma_sm_lepton_masses_uux_epem.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P5_Sigma_sm_lepton_masses_bbx_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P5_Sigma_sm_lepton_masses_ddx_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P5_Sigma_sm_lepton_masses_uux_mupmum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P6_Sigma_sm_lepton_masses_bbx_taptam.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P6_Sigma_sm_lepton_masses_ddx_taptam.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/MG5_PP_Z_P6_Sigma_sm_lepton_masses_uux_taptam.h"

namespace MG5_PP_Z {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_PP_Z_P1_Sigma_sm_lepton_masses_bbx_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P1_Sigma_sm_lepton_masses_ddx_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P1_Sigma_sm_lepton_masses_uux_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P2_Sigma_sm_lepton_masses_bbx_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P2_Sigma_sm_lepton_masses_ddx_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P2_Sigma_sm_lepton_masses_uux_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P3_Sigma_sm_lepton_masses_bbx_taptam>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P3_Sigma_sm_lepton_masses_ddx_taptam>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P3_Sigma_sm_lepton_masses_uux_taptam>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P4_Sigma_sm_lepton_masses_bbx_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P4_Sigma_sm_lepton_masses_ddx_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P4_Sigma_sm_lepton_masses_uux_epem>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P5_Sigma_sm_lepton_masses_bbx_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P5_Sigma_sm_lepton_masses_ddx_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P5_Sigma_sm_lepton_masses_uux_mupmum>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P6_Sigma_sm_lepton_masses_bbx_taptam>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P6_Sigma_sm_lepton_masses_ddx_taptam>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_Z_P6_Sigma_sm_lepton_masses_uux_taptam>();
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
      {Channel{{5, -5}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-11, 11}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{5, -5}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-13, 13}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{5, -5}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-15, 15}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{5, -5}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-11, 11}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{5, -5}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-13, 13}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{5, -5}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-5, 5}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{1, -1}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-1, 1}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{3, -3}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-3, 3}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}},
      {Channel{{2, -2}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-2, 2}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{4, -4}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {3, -3, 1, 1}, {{ColorFlowLeg{1, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}, Channel{{-4, 4}, {-15, 15}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}}, {-3, 3, 1, 1}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}}}}}
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
        "MG5_PP_Z::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_PP_Z

#endif
