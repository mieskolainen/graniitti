// Regge Amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEPARAM_H
#define MREGGEPARAM_H

// C++
#include <array>
#include <cstddef>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <string>
#include <tuple>
#include <vector>

// Own
#include "Graniitti/MModelTune.h"
#include "Graniitti/Photon/MPhotoDiss.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MFormFactor.h"
#include "Graniitti/Regge/MReggeProductionModel.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "json.hpp"

namespace gra {

// Regge amplitude parameters
namespace regge {

struct VertexParam {
  int                    first  = 0;
  int                    second = 0;
  std::array<FFParam, 2> transfer;
  std::array<FFParam, 2> offshell;
};

// Store one continuum vertex form-factor pair
struct VertexForm {
  FFParam transfer;
  FFParam offshell;
};

// Store one continuum zero-secondary-production veto
struct VetoParam {
  bool active = false;
  double M0 = 0.0;
  double c = 0.0;
};

// Store meson reggeization and the freezing virtuality in GeV^2
struct ReggeizeParam {
  bool active = false;
  double freeze_scale2 = 0.0;
};

struct PairParam {
  std::vector<int>         pdg;
  ReggeizeParam             reggeize;
  VetoParam                veto;
  std::map<int, VertexForm> forms;
  std::vector<VertexParam> channels;
};

enum class PermType { Auto, Charged, All };

// Store the result of one secondary-Reggeon vertex quantum-number check
struct VertexCheck {
  bool        applies            = false;
  bool        allowed            = true;
  int         representative_pdg = 0;
  std::string reason;
};

// Double-logarithmic correction anchored in value and energy derivative at W0
struct PhotoDLog {
  double scale2 = 0.0;
  double c = 0.0;
  double root0 = 0.0;
};

// Store one vector-meson photoproduction trajectory and meson-side cross-section slope
struct PhotoParam {
  double                        W0        = 0.0;
  double                        B_gammaPV = 0.0;  // d ln(|F_gammaPV|^2) / dt at t = 0 [GeV^-2]
  double                        a0        = 0.0;
  double                        ap        = 0.0;
  std::optional<flux::PhotoDissParam> diss;
  std::optional<PhotoDLog> dlog;
  std::vector<PhotoDLog> diss_dlog;
};

// Store one exchanged-meson Regge trajectory anchored at its physical pole
struct MesonTraj {
  double spin = 0.0;
  double ap   = 0.0;
};

// Store one typed central trajectory alias row and its SOFT exchange mapping
struct ExchangeParam {
  SoftExchangeRole role = SoftExchangeRole::Reggeon;
  std::vector<int> pdg;
  int              pole_spinX2 = -1;
  std::string      soft_exchange_name;
  SoftExchangeId   soft_exchange{0};
};

using Topology = std::vector<int>;

// Continuum specific Regge parameters
struct MReggeConParam {
  bool                                 multiregge_transfer_ext        = true;
  bool                                 multiregge_transfer_int        = false;
  bool                                 multiregge_secondary_exchanges = true;
  std::map<int, std::vector<Topology>> multiregge_topologies;
  PermType                             permutations = PermType::Auto;
  std::vector<PairParam>               MP;
  std::vector<PairParam>               XP;
  std::vector<PairParam>               GP;
  std::vector<int>                     form_pdgs;
};

struct Param {
  SoftModelPtr               soft_model;
  std::vector<ExchangeParam> exchanges;
  std::map<int, std::size_t> PDG_TO_INDEX;
  std::map<int, int>         PDG_TO_C;
  std::size_t                pomeron_trajectory = std::numeric_limits<std::size_t>::max();

  double  s0                 = 0.0;
  double  gp_alpha_min       = std::numeric_limits<double>::quiet_NaN();
  EtaMode photoprod_eta_mode = EtaMode::Rotating;
  std::map<ReggeProductionModel, bool> use_zeta;


  std::map<int, PhotoParam> photo_channels;
  std::map<int, MesonTraj>  meson_trajectories;

  std::map<ReggeProductionModel, double> omega;
  MReggeConParam con;
};

using ParamPtr = std::shared_ptr<const Param>;

std::string FormatPDGVector(const std::vector<int> &pdg);

// Validate continuum propagation and zero-secondary-production settings
void ValidateContinuumControls(const nlohmann::json &block, const std::string &label);

// Collect the current decay-tree final-state PDG codes
std::vector<int> FinalPDGs(const LORENTZSCALAR &lts);

// Convert one central amplitude index to a decay-tree slot
std::size_t DecayIndex(int amplitude_index, const LORENTZSCALAR &lts, const std::string &context);

// Compute the C-induced sign of one parsed continuum vertex
double VertexSign(const VertexParam &row, const MDecayBranch &left, const MDecayBranch &right, const Param &param);

std::size_t TrajectoryIndex(const Param &param, int pdg);

// Compute the mapped Regge trajectory at one generated transfer
double      Alpha(const Param &param, int pdg, double transfer);
int         CParity(const Param &param, int pdg);
int         PairCParity(const Param &param, int first_pdg, int second_pdg);

// Compute the canonical fixed-spin pole matching one trajectory alias
const MParticle &PoleRepresentative(const Param &param, const MPDG &pdg_table, int pdg);

// Compute whether one mapped trajectory contains a secondary Reggeon pole
bool IsSecondaryReggeonTrajectory(const Param &param, int pdg);

// Check charge and JPC compatibility of an f2, a2, rho or omega Reggeon vertex
VertexCheck CheckVertex(const Param &param, const gra::MPDG &pdg_table, int exchange_pdg, const std::vector<int> &pair_pdgs);

// Compute the channel-specific photoproduction parameters for one central state
const PhotoParam &Photo(const Param &param, int pdg);

// Compute the pole-anchored Regge trajectory for one exchanged meson
const MesonTraj &MesonTrajectory(const Param &param, int pdg);

// Resolve the configured continuum amplitude permutation construction
int PermCount(const std::vector<gra::MDecayBranch> &decaytree, std::size_t n_central, PermType mode);

// Compute the exchange C-parity sign induced by an antiparticle vertex
double AntiparticleSign(const Param &param, int exchange_pdg, int particle_pdg);

// Compute only an explicit continuum entry for a final-state pair
const PairParam *FindPair(const Param &param, const std::vector<int> &pdg, ReggeProductionModel model);
const PairParam &Pair(const Param &param, const std::vector<int> &pdg, ReggeProductionModel model);

// Compute whether pair entries contain one connected serial block
bool ContinuumLadderEntriesConnect(const Param &param, const std::vector<const PairParam *> &entries, ReggeProductionModel model, bool pomeron_only = false);

// Compute whether one continuum vertex is enabled in direct multi-Regge ladders
bool VertexAllowed(const Param &param, const VertexParam &vertex);

// Resolve direct multi-body continuum ladder permutations supported by a model
std::vector<std::vector<int>> LadderPermutations(const gra::LORENTZSCALAR &lts, const Param &param, std::size_t n_central, ReggeProductionModel model);

// Validate cached direct multi-body ladder permutations
void CheckLadderPermutations(const std::vector<std::vector<int>> &permutations, std::size_t n_central, const std::string &context);

// Validate cached serial and parallel multi-Regge topology partitions
void CheckMultiReggeTopologies(const std::vector<Topology> &topologies, std::size_t n_central, const std::string &context);

Param ReadParam(const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table, const MModelTune &tune);

// Construct one shared immutable Regge parameter block
ParamPtr ReadParamPtr(const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table, const MModelTune &tune);

// Compute one run owned immutable Regge parameter block
ParamPtr GetParam(MModelCache &cache, const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table);

}  // namespace regge
}  // namespace gra

#endif
