#ifndef GRANIITTI_AMPLITUDE_MG5_PP_ZJ_PROCESSES_H
#define GRANIITTI_AMPLITUDE_MG5_PP_ZJ_PROCESSES_H

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_bbx_epemg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_ddx_epemg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gb_epemb.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gbx_epembx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gd_epemd.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gdx_epemdx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gu_epemu.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gux_epemux.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_uux_epemg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_bbx_mupmumg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_ddx_mupmumg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gb_mupmumb.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gbx_mupmumbx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gd_mupmumd.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gu_mupmumu.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gux_mupmumux.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_uux_mupmumg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_bbx_taptamg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_ddx_taptamg.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gb_taptamb.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gbx_taptambx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gd_taptamd.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gdx_taptamdx.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gu_taptamu.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gux_taptamux.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg.h"

namespace MG5_PP_ZJ {

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_bbx_epemg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_ddx_epemg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gb_epemb>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gbx_epembx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gd_epemd>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gdx_epemdx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gu_epemu>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_gux_epemux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P1_Sigma_sm_lepton_masses_uux_epemg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_bbx_mupmumg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_ddx_mupmumg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gb_mupmumb>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gbx_mupmumbx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gd_mupmumd>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gdx_mupmumdx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gu_mupmumu>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_gux_mupmumux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P2_Sigma_sm_lepton_masses_uux_mupmumg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_bbx_taptamg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_ddx_taptamg>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gb_taptamb>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gbx_taptambx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gd_taptamd>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gdx_taptamdx>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gu_taptamu>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_gux_taptamux>();
    proc->initProc(param_card);
    out.push_back(std::move(proc));
  }
  {
    auto proc = std::make_unique<MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg>();
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
      {Channel{{5, -5}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-5, 5}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{1, -1}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-1, 1}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{3, -3}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-3, 3}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{21, 5}, {-11, 11, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{5, 21}, {-11, 11, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -5}, {-11, 11, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-5, 21}, {-11, 11, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 1}, {-11, 11, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{1, 21}, {-11, 11, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 3}, {-11, 11, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{3, 21}, {-11, 11, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -1}, {-11, 11, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-1, 21}, {-11, 11, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -3}, {-11, 11, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-3, 21}, {-11, 11, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 2}, {-11, 11, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{2, 21}, {-11, 11, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 4}, {-11, 11, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{4, 21}, {-11, 11, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -2}, {-11, 11, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-2, 21}, {-11, 11, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -4}, {-11, 11, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-4, 21}, {-11, 11, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{2, -2}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-2, 2}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{4, -4}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-4, 4}, {-11, 11, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-11}, {}}, AmplitudeTopologyNode{{11}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{5, -5}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-5, 5}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{1, -1}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-1, 1}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{3, -3}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-3, 3}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{21, 5}, {-13, 13, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{5, 21}, {-13, 13, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -5}, {-13, 13, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-5, 21}, {-13, 13, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 1}, {-13, 13, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{1, 21}, {-13, 13, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 3}, {-13, 13, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{3, 21}, {-13, 13, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -1}, {-13, 13, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-1, 21}, {-13, 13, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -3}, {-13, 13, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-3, 21}, {-13, 13, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 2}, {-13, 13, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{2, 21}, {-13, 13, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 4}, {-13, 13, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{4, 21}, {-13, 13, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -2}, {-13, 13, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-2, 21}, {-13, 13, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -4}, {-13, 13, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-4, 21}, {-13, 13, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{2, -2}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-2, 2}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{4, -4}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-4, 4}, {-13, 13, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-13}, {}}, AmplitudeTopologyNode{{13}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{5, -5}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-5, 5}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{1, -1}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-1, 1}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{3, -3}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-3, 3}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}},
      {Channel{{21, 5}, {-15, 15, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{5, 21}, {-15, 15, 5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{5}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -5}, {-15, 15, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-5, 21}, {-15, 15, -5}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-5}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 1}, {-15, 15, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{1, 21}, {-15, 15, 1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{1}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 3}, {-15, 15, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{3, 21}, {-15, 15, 3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{3}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -1}, {-15, 15, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-1, 21}, {-15, 15, -1}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-1}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -3}, {-15, 15, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-3, 21}, {-15, 15, -3}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-3}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{21, 2}, {-15, 15, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{2, 21}, {-15, 15, 2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{2}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{21, 4}, {-15, 15, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {8, 3, 1, 1, 3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}, Channel{{4, 21}, {-15, 15, 4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{4}, {}}}, {3, 8, 1, 1, 3}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{1, 0}}}}},
      {Channel{{21, -2}, {-15, 15, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-2, 21}, {-15, 15, -2}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-2}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{21, -4}, {-15, 15, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {8, -3, 1, 1, -3}, {{ColorFlowLeg{1, 2}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}, Channel{{-4, 21}, {-15, 15, -4}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{-4}, {}}}, {-3, 8, 1, 1, -3}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{1, 2}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 2}}}}},
      {Channel{{2, -2}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-2, 2}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{4, -4}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {3, -3, 1, 1, 8}, {{ColorFlowLeg{2, 0}, ColorFlowLeg{0, 1}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}, Channel{{-4, 4}, {-15, 15, 21}, {AmplitudeTopologyNode{{23}, {AmplitudeTopologyNode{{-15}, {}}, AmplitudeTopologyNode{{15}, {}}}}, AmplitudeTopologyNode{{21}, {}}}, {-3, 3, 1, 1, 8}, {{ColorFlowLeg{0, 1}, ColorFlowLeg{2, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{0, 0}, ColorFlowLeg{2, 1}}}}}
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
        "MG5_PP_ZJ::BuildSubprocesses: subprocess and channel counts disagree");
  }
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {
    out.push_back({std::move(matrix_elements[i]), std::move(channels[i])});
  }
}

}  // namespace MG5_PP_ZJ

#endif
