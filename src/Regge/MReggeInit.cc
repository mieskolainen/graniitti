// Regge process initialization and prepared amplitude data
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeMPInit.h"
#include "Graniitti/Regge/MReggeXPInit.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityNorm.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tensor/MTensorInit.h"
#include "Graniitti/Tensor/MTensorPomeron.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

// Compute (-1)^n for an integer-valued exponent
int IntegerPhaseSign(double exponent, const std::string &context) {
  constexpr double kComparisonTolerance = 1e-9;
  const double rounded = std::round(exponent);
  if (std::abs(exponent - rounded) > kComparisonTolerance) { throw std::invalid_argument(context + " has non-integral naturality phase"); }
  return (std::llround(rounded) % 2 == 0) ? 1 : -1;
}

// Fix the physical transverse norm to the measured diphoton partial width
// An isolated decay also removes the resonance line normalization
double GammaGammaScale(const PARAM_RES &res, double norm2, bool isolated_decay) {
  if (!(norm2 > 0.0) || !std::isfinite(norm2)) { throw std::invalid_argument("Amplitude initialization: invalid diphoton pole norm"); }
  double coupling = resonance::GammaGammaResonanceCoupling(res.p, res.modelparam.empty() ? MODELPARAM : res.modelparam);
  if (isolated_decay) { coupling *= math::msqrt(2.0 * res.p.mass * res.p.width); }
  return coupling / math::msqrt(norm2);
}

// Dispatch one ordered continuum pair to the selected production model
ReggeContinuumPole BuildContinuumPair(MProcessSetup &setup, const regge::Param &param, const ReggeProductionModel model, const MParticle &upper_exchange, const MParticle &lower_exchange, const MParticle &first, const MParticle &second,
                                      const std::string &context) {
  switch (model) {
    case ReggeProductionModel::MP:
      return mpom::PreparePair(setup, upper_exchange, lower_exchange, first, second);
    case ReggeProductionModel::XP:
      return xpom::PreparePair(setup, upper_exchange, lower_exchange, first, second);
    case ReggeProductionModel::GP:
      return gpom::PreparePair(setup, param, upper_exchange, lower_exchange, first, second, context);
    default:
      throw std::invalid_argument("Amplitude initialization: continuum pair requires MP, XP or GP");
  }
}

// Prepare every local exchange-helicity kernel used by a direct ladder
void SetupContinuumLadderPole(MProcessSetup &setup, const regge::Param &param, const ReggeProductionModel model) {
  if (model != ReggeProductionModel::MP && model != ReggeProductionModel::XP && model != ReggeProductionModel::GP) { throw std::invalid_argument("Amplitude initialization: unsupported continuum ladder spin model"); }
  const std::string model_name = ReggeProductionModelName(model);
  const std::string context    = "Amplitude initialization: " + model_name + " continuum ladder";
  auto             &pole       = setup.lts.process.CONT_LADDER_POLE;
  for (const auto &permutation : setup.lts.process.CONT_LADDER_PERMUTATIONS) {
    for (std::size_t pair = 0; pair < permutation.size(); pair += 2) {
      const std::size_t first_index  = regge::DecayIndex(permutation.at(pair), setup.lts, context);
      const std::size_t second_index = regge::DecayIndex(permutation.at(pair + 1), setup.lts, context);
      const auto       &first        = setup.lts.decaytree[first_index].p;
      const auto       &second       = setup.lts.decaytree[second_index].p;
      const auto       *entry        = regge::FindPair(param, {first.pdg, second.pdg}, model);
      if (entry == nullptr) { throw std::invalid_argument("Amplitude initialization: missing continuum ladder entry"); }
      for (const auto &row : entry->channels) {
        if (!regge::VertexAllowed(param, row)) { continue; }
        const ReggeContinuumPoleKey key = {first.pdg, second.pdg, row.first, row.second};
        if (pole.contains(key)) { continue; }
        const auto &upper_exchange = setup.lts.PDG.FindByPDG(row.first);
        const auto &lower_exchange = setup.lts.PDG.FindByPDG(row.second);
        if (upper_exchange.pdg == PDG::PDG_gamma || lower_exchange.pdg == PDG::PDG_gamma) {
          throw std::invalid_argument(
              "Amplitude initialization: continuum ladders do not "
              "support photon exchange");
        }
        auto prepared = BuildContinuumPair(setup, param, model, upper_exchange, lower_exchange, first, second, model_name + " continuum ladder");
        if (model == ReggeProductionModel::GP) {
          for (auto &vertex : prepared.gp_vertex) {
            gpom::PrepareLadder(vertex, model_name + " continuum ladder m=0");
          }
        }
        pole.emplace(key, std::move(prepared));
      }
    }
  }
  if (pole.empty()) { throw std::invalid_argument("Amplitude initialization: continuum ladder pole data is empty"); }
}


// Parse the forward beam helicity source mode during process initialization
ForwardVertexMode ParseForwardVertexMode(const std::string &mode) {
  if (mode == "helicity_residue") { return ForwardVertexMode::HelicityResidue; }
  if (mode == "unit_residue") { return ForwardVertexMode::UnitResidue; }
  throw std::invalid_argument(
      "Amplitude initialization: PARAM_REGGE.FORWARD_VERTEX should "
      "be 'helicity_residue' or 'unit_residue'");
}

// Compute the injected immutable SOFT model after validating process setup
const SoftModelPtr &RequireSetupSoftModel(const MProcessSetup &setup) {
  if (!setup.soft_model) { throw std::invalid_argument("Amplitude initialization requires a SOFT model"); }
  return setup.soft_model;
}

// Compute the injected immutable tune after validating process setup
const MModelTune &RequireSetupModelTune(const MProcessSetup &setup) {
  if (!setup.model_tune) { throw std::invalid_argument("Amplitude initialization requires a model tune"); }
  return *setup.model_tune;
}

// Compute the run owned model cache required by process initialization
MModelCache &RequireSetupModelCache(const MProcessSetup &setup) {
  if (setup.lts.model_cache == nullptr) { throw std::invalid_argument("Process initialization requires run owned model stores"); }
  return *setup.lts.model_cache;
}

// Compute the support restriction for one active photoproduction channel
std::string PhotoProductionError(ReggeProductionModel model, nuclear::CollisionType collision,
                               const MParticle &resonance, const std::array<int, 2> &exchange,
                               const MPDG &pdg, const nlohmann::json &general) {
  if ((exchange[0] == PDG::PDG_gamma) == (exchange[1] == PDG::PDG_gamma)) {
    return "photoproduction requires exactly one photon exchange";
  }
  const int target = exchange[exchange[0] == PDG::PDG_gamma ? 1 : 0];
  const bool ion = collision == nuclear::CollisionType::PA || collision == nuclear::CollisionType::EA ||
                   collision == nuclear::CollisionType::AA;
  const bool spin3 = model == ReggeProductionModel::TP && !ion && resonance.spinX2 == 6;
  if ((model == ReggeProductionModel::TP || ion) &&
      ((!spin3 && resonance.spinX2 != 2) || resonance.P != -1 || resonance.C != -1 || resonance.chargeX3 != 0)) {
    return "Tensor photoproduction requires neutral 1-- or 3-- states, nuclear production requires 1--";
  }
  if (model == ReggeProductionModel::TP) {
    const int isospin = pdg.FindByPDG(target).isospinX2;
    if (ion && isospin != 0 && isospin != 2) {
      return "Tensor nuclear photoproduction requires isoscalar or isovector target exchanges";
    }
  } else {
    const auto &regge = general.at("PARAM_REGGE");
    const bool pomeron = std::any_of(regge.at("EXCHANGES").begin(), regge.at("EXCHANGES").end(), [&](const auto &row) {
      const auto &aliases = row.at("pdg");
      return row.at("role") == "pomeron" && std::find(aliases.begin(), aliases.end(), target) != aliases.end();
    });
    if (!pomeron) { return "Regge photoproduction requires a Pomeron target exchange"; }
    const auto &photo = regge.at("photoprod");
    if (std::none_of(photo.begin(), photo.end(), [&](const auto &row) { return row.at(0) == resonance.pdg; })) {
      return "Regge photoproduction has no photoprod parameters for this resonance";
    }
  }
  return {};
}

// Compute true when one configured production channel contains a photon on a beam leg
// QED current validation is needed only for photon-emitting beam sides
bool UsesPhotonProductionLeg(const LORENTZSCALAR &lts, std::size_t leg) {
  const auto tree_has_photon = [leg](const std::vector<MDecayBranch> &tree) { return leg < tree.size() && tree[leg].p.pdg == PDG::PDG_gamma; };
  for (const auto &[name, resonance] : lts.process.RESONANCES) {
    (void)name;
    for (const auto &production : resonance.production) {
      if (tree_has_photon(production.tree)) { return true; }
    }
  }
  if (lts.process.ROOT_RES_ACTIVE) {
    for (const auto &production : lts.process.ROOT_RES.production) {
      if (tree_has_photon(production.tree)) { return true; }
    }
  }
  for (const auto &tree : lts.process.CONT_PRODUCTIONTREE) {
    if (tree_has_photon(tree)) { return true; }
  }
  return false;
}

// Validate every physical beam emitter of a direct non-collinear photon mode
void ValidateDirectPhotonEmitters(const MProcessSetup &setup) {
  if (setup.istate != "yy" && setup.istate != "ygg") { return; }
  for (std::size_t leg = 0; leg < 2; ++leg) {
    if (setup.istate == "ygg" && !gra::flux::SupportsPhotoDirection(setup.lts, static_cast<int>(leg + 1))) { continue; }
    gra::flux::ValidatePhotonEmitter(setup.lts, static_cast<int>(leg + 1), "Amplitude initialization " + setup.istate + "[" + setup.channel + "]");
  }
}

// Classify the active subprocess channel into branching setup modes
BranchingProcessFlags ClassifyBranchingProcess(const MProcessSetup &setup) {
  BranchingProcessFlags flags;
  flags.model              = setup.info.model;
  flags.continuum          = setup.info.continuum;
  flags.resonance          = setup.info.resonance == ResonanceType::Required ||
                             (setup.info.resonance == ResonanceType::Optional && !setup.lts.process.RESONANCES.empty());
  flags.jw_helicity_algebra = setup.info.jw_helicity_algebra;
  return flags;
}

// Compute true when a process-root resonance should set up its decay amplitude
bool RootResonanceDecayRequested(const LORENTZSCALAR &lts) { return lts.process.root_decay_mode != RootDecayMode::None && lts.process.root_resonance_pdg != 0; }

// Compute true when a process-root resonance decay is only an isolated proposal
bool RootResonanceDecayIsolated(const LORENTZSCALAR &lts) { return lts.process.root_decay_mode == RootDecayMode::Isolated; }




// Remove decay-amplitude normalization while retaining isolated branching ratio
// data
void ApplyIsolatedDecayNormalization(HELMatrix &hel) { hel.g_decay = std::complex<double>(1.0, 0.0); }

// Compute true when every tensor-Pomeron resonance uses the generic axial path
bool IsTensorAxialResonanceSet(const MProcessSetup &setup, const std::map<std::string, PARAM_RES> &resonances) {
  return setup.istate == "TP" && setup.channel == "RES" && !resonances.empty() &&
         std::all_of(resonances.begin(), resonances.end(), [](const auto &entry) { return entry.second.p.spinX2 == 2 && entry.second.p.P == 1 && entry.second.p.C == 1; });
}

// Validate flat continuum final-state shape before pole construction
void ValidateContinuumFinalState(const BranchingProcessFlags &flags, const std::vector<MDecayBranch> &tree) {
  if (!flags.continuum || decay::HasCascadedDecay(tree)) { return; }

  const std::size_t n_central = tree.size();
  if (n_central != 2 && n_central != 4 && n_central != 6) {
    throw std::invalid_argument(
        "Amplitude initialization: direct continuum final "
        "states require exactly 2, 4 or 6 "
        "central particles; got " +
        std::to_string(n_central));
  }
}

// Reject cascades when a Regge continuum selects a direct 4/6-body amplitude
void ValidateReggeContinuumTopology(const MProcessSetup &setup, const std::vector<MDecayBranch> &tree) {
  if ((setup.istate != "MP" && setup.istate != "XP" && setup.istate != "GP") || setup.channel != "CON") { return; }
  if (tree.size() == 2) { return; }
  if ((tree.size() == 4 || tree.size() == 6) && !decay::HasCascadedDecay(tree)) { return; }
  throw std::invalid_argument("Amplitude initialization: " + setup.istate +
                              "[CON] uses direct matrix elements for the 4/6-body "
                              "ladder and does not support cascades in those topologies");
}

// Validate every decay vertex whose branching ratio enters a physical JW
// amplitude
void ValidateTwoBodyDecayBranch(const MDecayBranch &branch) {
  if (branch.legs.empty()) { return; }
  if (branch.legs.size() != 2) {
    throw std::invalid_argument(
        "Amplitude initialization: physical helicity "
        "decay vertices must be two-body; PDG " +
        std::to_string(branch.p.pdg) + " has " + std::to_string(branch.legs.size()) + " daughters");
  }
  for (const auto &leg : branch.legs) { ValidateTwoBodyDecayBranch(leg); }
}

// Validate binary helicity cascades while retaining the TP axial multi-body root fallback
void ValidatePhysicalDecayArities(const MProcessSetup &setup, const BranchingProcessFlags &flags) {
  const auto &lts = setup.lts;
  if (RootResonanceDecayIsolated(lts)) { return; }
  const bool axial_binary = lts.decaytree.size() == 2 && IsTensorAxialResonanceSet(setup, lts.process.RESONANCES);
  if (!flags.jw_helicity_algebra && !axial_binary) { return; }
  if ((flags.resonance || RootResonanceDecayRequested(lts)) && lts.decaytree.size() != 2) {
    throw std::invalid_argument(
        "Amplitude initialization: physical resonance amplitudes "
        "require a "
        "two-body top-level "
        "decay; use an explicit cascade or &> isolated phase-space sampling");
  }
  for (const auto &branch : lts.decaytree) { ValidateTwoBodyDecayBranch(branch); }
}

// Compute true for the tensor-Pomeron vector-pair cascade topology
bool IsTensorVectorPairCascade(const std::vector<MDecayBranch> &tree) { return tree.size() == 2 && !tree[0].legs.empty() && !tree[1].legs.empty() && tree[0].legs.size() == 2 && tree[1].legs.size() == 2; }

// Validate one tensor-Pomeron vector branch as a direct V to two stable
// pseudoscalars decay
void ValidateTensorVectorCascadeBranch(const MDecayBranch &branch, const std::string &context) {
  if (branch.p.spinX2 != 2 || branch.legs.size() != 2) { throw std::invalid_argument(context + ": cascade branches must be spin-1 vectors with two daughters"); }
  for (const auto &daughter : branch.legs) {
    if (daughter.p.spinX2 != 0 || !daughter.legs.empty()) { throw std::invalid_argument(context + ": vector cascades require two stable pseudoscalar daughters"); }
  }
}

// Validate tensor-Pomeron resonance and continuum topology before coupling
// construction
void ValidateTensorTopology(const MProcessSetup &setup, const std::vector<MDecayBranch> &tree, const std::map<std::string, PARAM_RES> &resonances) {
  const bool tensor_channel = setup.istate == "TP";
  if (!tensor_channel) { return; }
  if (IsTensorAxialResonanceSet(setup, resonances)) { return; }
  if (tree.size() != 2) { throw std::invalid_argument("Amplitude initialization: " + setup.channel + " requires exactly two top-level central branches"); }

  const bool left_cascade  = !tree[0].legs.empty();
  const bool right_cascade = !tree[1].legs.empty();
  if (setup.channel == "RES+CON" && (left_cascade || right_cascade)) {
    throw std::invalid_argument(
        "Amplitude initialization: TP[RES+CON] supports "
        "only undecayed direct two-body states");
  }
  if (left_cascade != right_cascade) {
    throw std::invalid_argument(
        "Amplitude initialization: tensor vector-pair "
        "cascades require both branches to decay");
  }
  if (left_cascade) {
    ValidateTensorVectorCascadeBranch(tree[0], "Amplitude initialization");
    ValidateTensorVectorCascadeBranch(tree[1], "Amplitude initialization");
  }
}

// Reject isolated modes whose matrix element contains daughter-level dynamics
void ValidateIsolatedTopology(const MProcessSetup &setup, const BranchingProcessFlags &flags, const LORENTZSCALAR &lts, const std::map<std::string, PARAM_RES> &resonances) {
  if (!RootResonanceDecayIsolated(lts)) { return; }
  if (flags.continuum) {
    throw std::invalid_argument(
        "Amplitude initialization: &> requires a production-only "
        "resonance "
        "amplitude; " +
        setup.channel + " contains a daughter-dependent continuum matrix element");
  }
  if (flags.model == ReggeProductionModel::TP && (flags.resonance || flags.continuum) && !IsTensorAxialResonanceSet(setup, resonances)) {
    throw std::invalid_argument("Amplitude initialization: &> is not supported for " + setup.channel + " because its tensor production and decay vertices are intertwined");
  }
}

// Print one named setup section banner
void PrintBranchingSection(const std::string &title) {
  std::cout << std::endl;
  aux::PrintBar(".");
  std::cout << rang::fg::magenta << title << rang::fg::reset << std::endl;
  aux::PrintBar(".");
  std::cout << std::endl;
}

class CParityLookup {
 public:
  explicit CParityLookup(const MPDG &pdg_table) : pdg_table_(pdg_table) {}

  // Compute the C-parity for one exchange PDG from the PDG table
  int ParticleC(int pdg) {
    const auto cached = cache_.find(pdg);
    if (cached != cache_.end()) { return cached->second; }

    int c = 0;
    try {
      c = pdg_table_.FindByPDG(pdg).C;
    } catch (const std::exception &e) { throw std::invalid_argument("Amplitude initialization: exchange PDG " + std::to_string(pdg) + " is not in the PDG table: " + e.what()); }
    if (c != -1 && c != 1) { throw std::invalid_argument("Amplitude initialization: exchange PDG " + std::to_string(pdg) + " has undefined C-parity in the PDG table"); }
    cache_[pdg] = c;
    return c;
  }

  // Compute the product C-parity for one two-exchange channel
  int PairC(const std::vector<int> &channel) {
    if (channel.size() != 2) {
      throw std::invalid_argument(
          "Amplitude initialization: production channel "
          "must have two exchange PDGs");
    }
    return ParticleC(channel[0]) * ParticleC(channel[1]);
  }

  // Compute the product C parity for one typed two exchange channel
  int PairC(const std::array<int, 2> &channel) { return ParticleC(channel[0]) * ParticleC(channel[1]); }

 private:
  const MPDG        &pdg_table_;
  std::map<int, int> cache_;
};

// Compute whether one process uses MP, XP or GP Regge spin steering
bool UsesReggeSpinDefaults(const ReggeProductionModel model) { return model == ReggeProductionModel::MP || model == ReggeProductionModel::XP || model == ReggeProductionModel::GP; }

// Read complete model specific Regge spin steering tables at initialization
void ValidateSpinTables(const nlohmann::json &regge) {
  for (const std::string key : {"DERIVATIVE_FACTOR", "DECAY_BARRIERS", "FORWARD_VERTEX", "PHOTON_VERTEX", "TU_SIGN"}) {
    const auto &table = regge.at(key);
    if (!table.is_object() || table.size() != 3) { throw std::invalid_argument("Amplitude initialization: PARAM_REGGE." + key + " must contain MP, XP and GP"); }
    for (const std::string model : {"MP", "XP", "GP"}) {
      if (!table.contains(model)) { throw std::invalid_argument("Amplitude initialization: PARAM_REGGE." + key + " is missing " + model); }
      const auto &value = table.at(model);
      if (key == "DERIVATIVE_FACTOR" || key == "DECAY_BARRIERS") {
        if (!value.is_boolean()) { throw std::invalid_argument("Amplitude initialization: PARAM_REGGE." + key + "." + model + " must be boolean"); }
      } else if (key == "FORWARD_VERTEX") {
        (void)ParseForwardVertexMode(value.get<std::string>());
      } else if (key == "PHOTON_VERTEX") {
        qed::ValidatePhotonMode(value.get<std::string>(), "Amplitude initialization PARAM_REGGE.PHOTON_VERTEX." + model);
      } else {
        const auto sign = value.get<std::string>();
        if (sign != "auto" && sign != "positive" && sign != "negative") { throw std::invalid_argument("Amplitude initialization: PARAM_REGGE.TU_SIGN." + model + " must be auto, positive or negative"); }
      }
    }
  }
}

// Apply shared and model-specific spin steering defaults from GENERAL.json
void ApplySpinDefaults(MProcessSetup &setup, const ReggeProductionModel model) {
  bool default_spingen        = true;
  bool default_spindec        = true;
  bool default_forward_noflip = false;

  const MModelTune &tune = RequireSetupModelTune(setup);
  const auto       &spin = tune.General("PARAM_SPIN");
  default_spingen        = spin.at("SPINGEN");
  default_spindec        = spin.at("SPINDEC");
  default_forward_noflip = spin.at("FORWARD_NOFLIP");
  const bool decay_sym   = spin.at("DECAY_SYM");

  if (!setup.spingen_user) { setup.lts.process.SPINGEN = default_spingen; }
  if (!setup.spindec_user) { setup.lts.process.SPINDEC = default_spindec; }
  setup.lts.process.FORWARD_NOFLIP = default_forward_noflip;
  setup.lts.amplitude.DECAY_SYM    = decay_sym;
  if (setup.lts.process.root_decay_mode == RootDecayMode::Isolated) { setup.lts.amplitude.DECAY_SYM = false; }
  if (model != ReggeProductionModel::MP && setup.lts.process.MP_FRAME != "null") { throw std::invalid_argument("Amplitude initialization: @MP_FRAME applies only to MP processes"); }

  if (!UsesReggeSpinDefaults(model)) { return; }
  const auto       &regge          = tune.General("PARAM_REGGE");
  ValidateSpinTables(regge);
  const std::string name           = ReggeProductionModelName(model);
  const std::string forward_vertex = regge.at("FORWARD_VERTEX").at(name);
  const std::string photon_vertex  = regge.at("PHOTON_VERTEX").at(name);
  const auto       &mmax           = tune.Numerics("NUMERICS_REGGE").at("MMAX");
  if (!mmax.is_number_integer() || mmax.get<long double>() < 0.0L || mmax.get<long double>() > (std::numeric_limits<int>::max() - 1) / 2) {
    throw std::invalid_argument("Amplitude initialization: NUMERICS_REGGE.MMAX must be a supported non-negative integer");
  }
  const int         MMAX           = mmax.get<int>();
  const std::string tu_sign        = regge.at("TU_SIGN").at(name);
  if (model == ReggeProductionModel::MP) {
    const std::string default_mp_frame = regge.at("MP_FRAME");
    if (default_mp_frame != "CS" && default_mp_frame != "HX" && default_mp_frame != "CM") {
      throw std::invalid_argument(
          "Amplitude initialization: PARAM_REGGE.MP_FRAME should be 'CS', "
          "'HX' or 'CM'");
    }
    if (setup.lts.process.MP_FRAME == "null") { setup.lts.process.MP_FRAME = default_mp_frame; }
    if (setup.lts.process.MP_FRAME != "CS" && setup.lts.process.MP_FRAME != "HX" && setup.lts.process.MP_FRAME != "CM") { throw std::invalid_argument("Amplitude initialization: @MP_FRAME should be 'CS', 'HX' or 'CM'"); }
  }
  setup.lts.process.DECAY_BARRIER     = regge.at("DECAY_BARRIERS").at(name);
  setup.lts.process.DERIVATIVE_FACTOR = regge.at("DERIVATIVE_FACTOR").at(name).get<bool>();
  setup.lts.process.FORWARD_VERTEX    = ParseForwardVertexMode(forward_vertex);
  setup.lts.process.PHOTON_VERTEX     = photon_vertex;
  if (!setup.mmax_user) { setup.lts.process.MMAX = MMAX; }
  (void)gpom::AnalyticMIndex(0, setup.lts.process.MMAX, "Amplitude initialization: effective NUMERICS_REGGE.MMAX");
  setup.lts.process.TU_SIGN = tu_sign;
}

// Warn when an MRegge running width has no pole-normalized two-body profile
void WarnRunningWidthFallback(const MProcessSetup &setup, const PARAM_RES &res, const std::vector<MParticle> &daughter, const BranchingProcessFlags &flags) {
  if (res.BW != BreitWigner::RunningWidth || !(flags.ResonanceModel(ReggeProductionModel::MP) || flags.ResonanceModel(ReggeProductionModel::XP) || flags.ResonanceModel(ReggeProductionModel::GP))) { return; }

  const bool two_body  = daughter.size() == 2;
  const bool pole_open = two_body && gra::kinematics::DecayMomentum(res.p.mass, daughter[0].mass, daughter[1].mass) > 0.0;
  const bool ls_ready  = !setup.lts.process.DECAY_BARRIER || !res.hel_decay.ls_components.empty();
  if (two_body && pole_open && ls_ready) { return; }

  gra::aux::PrintWarning();
  std::cout << rang::fg::red << "Amplitude initialization: resonance " << res.p.name
            << " has no pole-normalized two-body LS running-width profile, "
               "using fixed-width behavior"
            << rang::fg::reset << std::endl;
}

// Build one resonance decay helicity structure with isolated-decay fallback
void PrepareResonanceDecay(MProcessSetup &setup, PARAM_RES &res, const std::vector<MParticle> &daughter, const BranchingProcessFlags &flags, const bool isolated) {
  try {
    res.hel_decay = setup.ProcessHelicityStructure(res.p, daughter, false, !isolated);
  } catch (const MissingHelicityData &error) {
    if (!isolated) { throw; }
    gra::aux::PrintWarning();
    std::cout << rang::fg::red
              << "Amplitude initialization: &> isolated resonance decay data "
                 "not found; manual BR factor is not inferred automatically: "
              << error.what() << rang::fg::reset << std::endl;
    res.hel_decay = HELMatrix();
    res.hel_decay.alpha_ls.Clear();
  }
  if (isolated) { ApplyIsolatedDecayNormalization(res.hel_decay); }
  WarnRunningWidthFallback(setup, res, daughter, flags);
  res.production_model = flags.model;
}

// Setup resonance decay helicity structures and production channel selections
void SetupResonanceDecayStructures(MProcessSetup &setup, const BranchingProcessFlags &flags) {
  const bool has_cascaded_branch = std::any_of(setup.lts.decaytree.cbegin(), setup.lts.decaytree.cend(), [](const auto &branch) { return !branch.legs.empty(); });
  if (!flags.resonance && !has_cascaded_branch) { return; }

  PrintBranchingSection("Amplitude initialization: Decay spin structure");

  if (flags.resonance) {
    std::vector<MParticle> daughter;
    daughter.reserve(setup.lts.decaytree.size());
    for (const auto &branch : setup.lts.decaytree) { daughter.push_back(branch.p); }
    const bool isolated = RootResonanceDecayIsolated(setup.lts);
    for (auto &[name, res] : setup.lts.process.RESONANCES) {
      (void)name;
      PrepareResonanceDecay(setup, res, daughter, flags, isolated);
    }
  }

  for (const auto &i : indices(setup.lts.decaytree)) {
    if (!setup.lts.decaytree[i].legs.empty()) { setup.ProcessHelicityTree(setup.lts.decaytree[i], false, !setup.isolate); }
  }
}

// Store one ordered beam side production path
struct ProductionPath {
  std::array<int, 2> exchange;
  std::size_t        source = 0;
};

// Compute the resonance production channels of one active model
const std::vector<RES_PRODUCTION_CHANNEL> &ProductionChannels(const PARAM_RES &RES, ReggeProductionModel model) {
  if (model == ReggeProductionModel::MP) { return RES.MP.channels; }
  if (model == ReggeProductionModel::XP) { return RES.XP.channels; }
  if (model == ReggeProductionModel::GP) { return RES.GP.channels; }
  throw std::invalid_argument("Amplitude initialization: resonance production requires MP, XP or GP");
}

// Select and order resonance channels compatible with the central C parity
std::vector<ProductionPath> SelectProductionChannels(const PARAM_RES &RES, ReggeProductionModel model,
                                                    CParityLookup &exchange_c) {
  if (model != ReggeProductionModel::MP && model != ReggeProductionModel::XP && model != ReggeProductionModel::GP) { return {}; }
  const auto                                &channels = ProductionChannels(RES, model);
  std::vector<ProductionPath>                output;
  std::map<std::pair<int, int>, std::size_t> seen;
  for (const auto &i : indices(channels)) {
    const auto &exchange = channels[i].exchange;
    std::cout << rang::fg::yellow << "Production input[" << i << "]: [" << exchange[0] << ", " << exchange[1] << "]" << rang::fg::reset << std::endl;
    const auto key                  = std::minmax(exchange[0], exchange[1]);
    const auto [position, inserted] = seen.emplace(key, i);
    if (!inserted) {
      throw std::invalid_argument("Amplitude initialization: duplicate unordered production channel [" + std::to_string(exchange[0]) + ", " + std::to_string(exchange[1]) + "] duplicates production[" + std::to_string(position->second) +
                                  "]");
    }
    if (RES.p.C != 0) {
      const int pair_c = exchange_c.PairC(exchange);
      if (pair_c != RES.p.C) {
        throw std::invalid_argument("Amplitude initialization: production channel violates resonance C parity");
      }
    }
    output.push_back({exchange, i});
    if (exchange[0] != exchange[1]) { output.push_back({{exchange[1], exchange[0]}, i}); }
  }
  return output;
}

// Build one exchange to beam scattering branch
MDecayBranch ProductionBranch(const MPDG &pdg, int exchange, int beam) {
  MDecayBranch branch;
  branch.p     = pdg.FindByPDG(exchange);
  branch.depth = 0;
  branch.name  = std::to_string(exchange) + "#0";
  branch.legs.resize(2);
  for (auto &leg : branch.legs) {
    leg.p     = pdg.FindByPDG(beam);
    leg.depth = 1;
    leg.name  = std::to_string(beam) + "#1";
  }
  return branch;
}

// Build one complete runtime resonance production channel
RES_PRODUCTION BuildProduction(MProcessSetup &setup, PARAM_RES &RES, const regge::ParamPtr &param, const ProductionPath &path, ReggeProductionModel model) {
  RES_PRODUCTION output;
  output.tree = {ProductionBranch(setup.lts.PDG, path.exchange[0], setup.lts.beam1.pdg), ProductionBranch(setup.lts.PDG, path.exchange[1], setup.lts.beam2.pdg)};
  std::vector<MParticle> legs;
  for (auto &branch : output.tree) {
    legs.push_back(branch.p);
    if (model != ReggeProductionModel::GP) { setup.ProcessHelicityTree(branch, true, true); }
  }
  const auto &channel = ProductionChannels(RES, model).at(path.source);
  if (model == ReggeProductionModel::GP) {
    if (!param) { throw std::logic_error("Amplitude initialization: missing GP Regge model"); }
    output.hel = gpom::PrepareResonance(RES.p, legs, channel, *param, setup.lts.PDG, setup.lts.process.MMAX, setup.model_tune->Global().coupling_min);
    if (channel.width_derived) {
      gpom::NormalizeGammaGamma(output.hel, RES, legs, setup.lts.process.DERIVATIVE_FACTOR,
                                  RootResonanceDecayIsolated(setup.lts));
    }
    return output;
  }
  output.pole = (model == ReggeProductionModel::MP ? mpom::PrepareResonance : xpom::PrepareResonance)(setup, RES, channel, legs);
  if (channel.width_derived) { rspin::NormalizeGammaGamma(*output.pole, RES, RootResonanceDecayIsolated(setup.lts)); }
  output.hel = output.pole->helicity;
  if (model == ReggeProductionModel::XP) { xpom::PrintPole("XP resonance pole operator", RES.p, legs, *output.pole); }
  return output;
}

// Build the immutable resonance production plan
void BuildResonancePlan(MProcessSetup &setup, const BranchingProcessFlags &flags, CParityLookup &exchange_c) {
  if (!flags.jw_helicity_algebra || !flags.resonance) { return; }

  regge::ParamPtr regge_param;
  if (flags.ResonanceModel(ReggeProductionModel::GP)) {
    std::vector<int> final_pdgs;
    final_pdgs.reserve(setup.lts.decaytree.size());
    for (const auto &branch : setup.lts.decaytree) { final_pdgs.push_back(branch.p.pdg); }
    regge_param = regge::GetParam(RequireSetupModelCache(setup), final_pdgs, setup.lts.PDG);
  }

  PrintBranchingSection("Amplitude initialization: Production spin structure");
  for (auto &[name, RES] : setup.lts.process.RESONANCES) {
    if (flags.model == ReggeProductionModel::MP && !RES.UsesMPPolarizationMode()) { throw std::invalid_argument("Amplitude initialization: unsupported MP spin basis " + RES.spin_basis); }
    if (!RES.C_from_card && RES.p.C == 0) {
      try {
        RES.p.C = setup.lts.PDG.FindByPDG(RES.p.pdg).C;
      } catch (...) {
        // Charged or user-defined resonances may not have a C eigenvalue in the
        // particle table
      }
    }
    if (flags.model == ReggeProductionModel::MP) { mpom::PreparePolarization(RES, setup.lts); }
    RES.production.clear();
    const auto selected = SelectProductionChannels(RES, flags.model, exchange_c);
    const auto collision = nuclear::ClassifyCollision(setup.lts.beam1.pdg, setup.lts.beam2.pdg);
    for (const auto &k : indices(selected)) {
      const bool active = ProductionChannels(RES, flags.model)[selected[k].source].Active(setup.model_tune->Global().coupling_min);
      if (active && collision != nuclear::CollisionType::PP && collision != nuclear::CollisionType::Invalid) {
        const auto error = PhotoProductionError(flags.model, collision, RES.p, selected[k].exchange,
                                                setup.lts.PDG, setup.model_tune->General());
        if (!error.empty()) {
          throw std::invalid_argument("Amplitude initialization: " + setup.istate + "[" + setup.channel + "] " +
                                       nuclear::CollisionName(collision) + " " + name + ": " + error);
        }
      }
      std::cout << rang::fg::yellow << "Production channel[" << k << "]: [" << selected[k].exchange[0] << ", " << selected[k].exchange[1] << "]" << rang::fg::reset << std::endl;
      auto production = BuildProduction(setup, RES, regge_param, selected[k], flags.model);
      if (active) { RES.production.push_back(std::move(production)); }
    }
    if (flags.model == ReggeProductionModel::MP && !selected.empty()) { mpom::PrintPolarization(RES, setup.lts.process.QMETRICS); }

    (void)name;
    std::cout << std::endl;
  }
  // Omit zero resonances only after validating every supplied production channel
  std::erase_if(setup.lts.process.RESONANCES, [](const auto &entry) { return entry.second.production.empty(); });
}

// Build one immutable direct continuum ladder plan
void BuildContinuumLadderPlan(MProcessSetup &setup, const ReggeProductionModel model, const std::size_t n_central) {
  setup.lts.process.CONT_LADDER_POLE.clear();
  if (setup.lts.decaytree.size() != n_central) { throw std::invalid_argument("Amplitude initialization: continuum plan multiplicity mismatch"); }

  const SoftModelPtr &soft_model = RequireSetupSoftModel(setup);
  std::vector<int>    final_pdgs;
  final_pdgs.reserve(n_central);
  for (const auto &branch : setup.lts.decaytree) {
    if (branch.p.spinX2 != 0) {
      throw std::invalid_argument("Amplitude initialization: continuum ladders require spin-zero final states");
    }
    final_pdgs.push_back(branch.p.pdg);
  }
  const auto param         = regge::GetParam(RequireSetupModelCache(setup), final_pdgs, setup.lts.PDG);
  const auto topology_bank = param->con.multiregge_topologies.find(static_cast<int>(n_central));
  if (topology_bank == param->con.multiregge_topologies.end() || topology_bank->second.empty()) {
    throw std::invalid_argument("Amplitude initialization: PARAM_CON.MULTI has no partitions for " + std::to_string(n_central) + " central particles");
  }
  const bool parallel                     = std::any_of(topology_bank->second.cbegin(), topology_bank->second.cend(), [](const auto &topology) { return topology.size() > 1; });
  setup.lts.process.MULTIREGGE_TOPOLOGIES = topology_bank->second;
  // The parallel subladder bank stores diagonal production vertices
  const SoftExchangeId pomeron = param->exchanges.at(param->pomeron_trajectory).soft_exchange;
  if (parallel && !soft_model->Exchange(pomeron).CouplingMatrix().IsDiagonal()) {
    throw std::invalid_argument(
        "Amplitude initialization: parallel multi-Regge topologies require "
        "diagonal production Pomeron residues");
  }
  if (parallel && setup.excitation != 0) {
    throw std::invalid_argument(
        "Amplitude initialization: parallel multi-Regge topologies require "
        "elastic forward hadrons");
  }
  setup.lts.process.CONT_LADDER_PERMUTATIONS = regge::LadderPermutations(setup.lts, *param, n_central, model);
  if (setup.lts.process.CONT_LADDER_PERMUTATIONS.empty()) {
    const std::string model_name = ReggeProductionModelName(model);
    throw std::invalid_argument("Amplitude initialization: no model-supported " + model_name +
                                " direct continuum ladder "
                                "permutation");
  }
  SetupContinuumLadderPole(setup, *param, model);
}

using ReggeBeamExchangePair = std::pair<int, int>;

// Compute whether one beam exchange can produce the selected inelastic leg
bool SupportsReggeForwardExcitation(const regge::Param &param, const int exchange_pdg) {
  if (exchange_pdg == PDG::PDG_gamma) { return true; }
  const std::size_t    trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  const SoftExchangeId exchange   = param.exchanges.at(trajectory).soft_exchange;
  return param.soft_model->Exchange(exchange).Role() == SoftExchangeRole::Pomeron;
}

// Resolve the explicit pair entries of one direct continuum ladder permutation
std::vector<const regge::PairParam *> ReggeLadderEntries(const MProcessSetup &setup, const regge::Param &param, const std::vector<int> &permutation, const ReggeProductionModel model) {
  const std::string                     model_name     = ReggeProductionModelName(model);
  constexpr int                         central_offset = 3;
  std::vector<const regge::PairParam *> entries;
  entries.reserve(permutation.size() / 2);
  for (std::size_t pair = 0; pair < permutation.size(); pair += 2) {
    const int first  = permutation.at(pair) - central_offset;
    const int second = permutation.at(pair + 1) - central_offset;
    if (first < 0 || second < 0 || static_cast<std::size_t>(first) >= setup.lts.decaytree.size() || static_cast<std::size_t>(second) >= setup.lts.decaytree.size()) {
      throw std::invalid_argument("Amplitude initialization: invalid " + model_name + " ladder index");
    }
    const auto *entry = regge::FindPair(param, {setup.lts.decaytree[first].p.pdg, setup.lts.decaytree[second].p.pdg}, model);
    if (entry == nullptr) { throw std::invalid_argument("Amplitude initialization: " + model_name + " ladder entry disappeared"); }
    entries.push_back(entry);
  }
  return entries;
}

// Extend one connected continuum exchange chain to its lower beam exchange
void AppendConnectedReggeLadderPairs(const regge::Param &param, const std::vector<const regge::PairParam *> &entries, const std::size_t entry_index, const int upper_exchange, const int previous_exchange,
                                     const std::size_t previous_trajectory, const ReggeProductionModel model, std::vector<ReggeBeamExchangePair> &pairs) {
  for (const auto &row : entries.at(entry_index)->channels) {
    if (!regge::VertexAllowed(param, row)) { continue; }
    const bool connected = model == ReggeProductionModel::GP ? regge::TrajectoryIndex(param, row.first) == previous_trajectory : row.first == previous_exchange;
    if (!connected) { continue; }
    if (entry_index + 1 == entries.size()) {
      pairs.emplace_back(upper_exchange, row.second);
      continue;
    }
    AppendConnectedReggeLadderPairs(param, entries, entry_index + 1, upper_exchange, row.second, regge::TrajectoryIndex(param, row.second), model, pairs);
  }
}

// Collect the outer beam exchange pairs of every connected continuum ladder
void AppendReggeLadderBeamPairs(const MProcessSetup &setup, const regge::Param &param, const ReggeProductionModel model, std::vector<ReggeBeamExchangePair> &pairs) {
  const std::string model_name = ReggeProductionModelName(model);
  for (const auto &permutation : setup.lts.process.CONT_LADDER_PERMUTATIONS) {
    const auto entries = ReggeLadderEntries(setup, param, permutation, model);
    if (entries.size() < 2) { throw std::invalid_argument("Amplitude initialization: incomplete " + model_name + " ladder"); }
    for (const auto &row : entries.front()->channels) {
      if (!regge::VertexAllowed(param, row)) { continue; }
      AppendConnectedReggeLadderPairs(param, entries, 1, row.first, row.second, regge::TrajectoryIndex(param, row.second), model, pairs);
    }
  }
}

// Collect all initialized MRegge beam exchange orderings
std::vector<ReggeBeamExchangePair> ReggeProductionOrderings(const MProcessSetup &setup, const BranchingProcessFlags &flags, const regge::Param &param) {
  std::vector<ReggeBeamExchangePair> pairs;
  if (flags.resonance) {
    for (const auto &[name, resonance] : setup.lts.process.RESONANCES) {
      (void)name;
      for (const auto &production : resonance.production) {
        const auto &tree = production.tree;
        if (tree.size() != 2) {
          throw std::invalid_argument(
              "Amplitude initialization: MRegge resonance "
              "production ordering must have two exchanges");
        }
        pairs.emplace_back(tree[0].p.pdg, tree[1].p.pdg);
      }
    }
  }
  if (flags.continuum && setup.lts.decaytree.size() == 2) {
    for (const auto &channel : setup.lts.process.CONT_PRODUCTION) {
      if (channel.size() != 2) {
        throw std::invalid_argument(
            "Amplitude initialization: MRegge continuum "
            "production ordering must have two exchanges");
      }
      pairs.emplace_back(channel[0], channel[1]);
    }
  }
  if (flags.continuum && setup.lts.decaytree.size() > 2) { AppendReggeLadderBeamPairs(setup, param, flags.model, pairs); }
  return pairs;
}

// Require one supported production ordering for every sampled excitation side
void ValidateReggeExcitationOrderings(const MProcessSetup &setup, const BranchingProcessFlags &flags) {
  if (setup.excitation == 0) { return; }
  std::vector<int> final_pdgs;
  final_pdgs.reserve(setup.lts.decaytree.size());
  for (const auto &branch : setup.lts.decaytree) { final_pdgs.push_back(branch.p.pdg); }
  const auto param    = regge::GetParam(RequireSetupModelCache(setup), final_pdgs, setup.lts.PDG);
  const auto pairs    = ReggeProductionOrderings(setup, flags, *param);
  if (flags.continuum && setup.lts.process.DISSOCIATION == DissociationType::Hera) {
    throw std::invalid_argument("PARAM_NSTAR.MODEL: hera has no hadronic continuum profile");
  }
  if (flags.continuum && setup.lts.process.PHOTO_DISSOCIATION == DissociationType::Hera) {
    for (const auto &channel : setup.lts.process.CONT_PRODUCTION) {
      if (std::count(channel.begin(), channel.end(), PDG::PDG_gamma) == 1) {
        throw std::invalid_argument("PARAM_NSTAR.MODEL: hera requires a vector resonance, not a photoproduction continuum");
      }
    }
  }
  const auto supports = [&](const bool upper, const bool lower) {
    return std::any_of(pairs.begin(), pairs.end(), [&](const ReggeBeamExchangePair &pair) { return (!upper || SupportsReggeForwardExcitation(*param, pair.first)) && (!lower || SupportsReggeForwardExcitation(*param, pair.second)); });
  };
  if (setup.excitation == 1 && nuclear::IsProton(setup.lts.beam1.pdg) && !supports(true, false)) {
    throw std::invalid_argument(
        "Amplitude initialization: NSTARS=1 has no MRegge "
        "production ordering for upper forward excitation");
  }
  if (setup.excitation == 1 && nuclear::IsProton(setup.lts.beam2.pdg) && !supports(false, true)) {
    throw std::invalid_argument(
        "Amplitude initialization: NSTARS=1 has no MRegge "
        "production ordering for lower forward excitation");
  }
  if (setup.excitation == 2 && !supports(true, true)) {
    throw std::invalid_argument(
        "Amplitude initialization: NSTARS=2 has no MRegge "
        "production ordering for double forward excitation");
  }
  for (const auto &[name, resonance] : setup.lts.process.RESONANCES) {
    for (const auto &production : resonance.production) {
      const auto &tree = production.tree;
      const bool photo = std::any_of(tree.begin(), tree.end(), [](const auto &leg) { return leg.p.pdg == PDG::PDG_gamma; });
      const auto type = photo ? setup.lts.process.PHOTO_DISSOCIATION : setup.lts.process.DISSOCIATION;
      if (type != DissociationType::Hera) { continue; }
      if (resonance.p.spinX2 != 2) { throw std::invalid_argument("PARAM_NSTAR.MODEL: hera requires vector meson photoproduction"); }
      if (tree.size() != 2) { throw std::invalid_argument("PARAM_NSTAR.MODEL: hera requires gamma-Pomeron production"); }
      const bool upper_gamma = tree[0].p.pdg == PDG::PDG_gamma;
      const bool lower_gamma = tree[1].p.pdg == PDG::PDG_gamma;
      const int exchange = upper_gamma ? tree[1].p.pdg : tree[0].p.pdg;
      if (upper_gamma == lower_gamma || regge::TrajectoryIndex(*param, exchange) != param->pomeron_trajectory) {
        throw std::invalid_argument("PARAM_NSTAR.MODEL: hera requires gamma-Pomeron production");
      }
      const int target_pdg = upper_gamma ? setup.lts.beam2.pdg : setup.lts.beam1.pdg;
      if (nuclear::IsProton(target_pdg) && !regge::Photo(*param, resonance.p.pdg).diss) {
        throw std::invalid_argument("Amplitude initialization: photoprod_diss has no proton-dissociative row for " + name);
      }
    }
  }
}

// Setup process-root resonance decay helicity structures for physical or
// isolated arrows
void SetupRootResonanceDecayStructure(MProcessSetup &setup, const BranchingProcessFlags &flags) {
  setup.lts.process.ROOT_RES_ACTIVE = false;
  setup.lts.process.ROOT_RES        = gra::PARAM_RES();

  if (!RootResonanceDecayRequested(setup.lts)) { return; }

  const int      pdg = setup.lts.process.root_resonance_pdg;
  gra::PARAM_RES res;
  res.p = setup.lts.PDG.FindByPDG(pdg);

  std::vector<MParticle> daughter;
  for (const auto &branch : setup.lts.decaytree) { daughter.push_back(branch.p); }

  PrintBranchingSection("Amplitude initialization: Root decay spin structure");
  PrepareResonanceDecay(setup, res, daughter, flags, RootResonanceDecayIsolated(setup.lts));

  setup.lts.process.ROOT_RES        = res;
  setup.lts.process.ROOT_RES_ACTIVE = true;
}

// Initialize shared physical or isolated decay state before sampling
void InitializeDecayState(MProcessSetup &setup, const BranchingProcessFlags &flags) {
  ValidatePhysicalDecayArities(setup, flags);
  if (setup.isolate && setup.flat_mass2) {
    throw std::invalid_argument(
        "Process initialization: &> isolated decay sampling requires the "
        "default Breit-Wigner mass proposal; @FLATMASS2:true is not supported");
  }
  setup.symmetry_factor = decay::FinalStateSymmetryFactor(setup.lts.decaytree, setup.isolate);
  if (flags.resonance && setup.lts.process.RESONANCES.empty()) {
    throw std::invalid_argument("Amplitude initialization: RES requires at least one active resonance");
  }
  if (flags.jw_helicity_algebra || flags.model != ReggeProductionModel::None) {
    SetupResonanceDecayStructures(setup, flags);
  }
  SetupRootResonanceDecayStructure(setup, flags);
}

// Initialize physical Regge or tensor branching state
void InitializePhysicalProcessState(MProcessSetup &setup, const BranchingProcessFlags &flags) {
  setup.lts.process.REGGE_MODEL = flags.model;
  if (setup.lts.decaytree.empty()) {
    std::cout << rang::fg::red << "Physical process initialization: decaytree.size() == 0 !" << rang::fg::reset << std::endl;
    return;
  }

  SetupTensorProcessModel(setup);
  ValidateContinuumFinalState(flags, setup.lts.decaytree);
  ValidateReggeContinuumTopology(setup, setup.lts.decaytree);
  ValidateTensorTopology(setup, setup.lts.decaytree, setup.lts.process.RESONANCES);
  ValidateIsolatedTopology(setup, flags, setup.lts, setup.lts.process.RESONANCES);
  if (setup.istate == "TP" && setup.channel == "CON" && !decay::HasCascadedDecay(setup.lts.decaytree) && setup.lts.decaytree.size() != 2) {
    throw std::invalid_argument(
        "Amplitude initialization: TP[CON] does not "
        "support direct multi-body final states; "
        "use vector cascades");
  }
  CParityLookup exchange_c(setup.lts.PDG);

  ApplySpinDefaults(setup, flags.model);
  InitializeDecayState(setup, flags);
  BuildResonancePlan(setup, flags, exchange_c);
  setup.lts.process.CONT_LADDER_PERMUTATIONS.clear();
  setup.lts.process.CONT_LADDER_POLE.clear();
  setup.lts.process.MULTIREGGE_TOPOLOGIES.clear();
}

// Initialize process-independent root decay and phase-space state
void InitializeGenericProcessState(MProcessSetup &setup) {
  if (setup.lts.decaytree.empty()) { return; }
  const auto flags = ClassifyBranchingProcess(setup);
  ApplySpinDefaults(setup, ReggeProductionModel::None);
  ValidateDirectPhotonEmitters(setup);
  InitializeDecayState(setup, flags);
}

// Prepare the common MP, XP and GP density spaces before any event is evaluated
void ConfigureReggeQMetrics(MProcessSetup &setup, const ReggeProductionModel model) {
  auto &lts = setup.lts;
  if (!lts.process.QMETRICS) { return; }
  const auto &numerics = setup.model_tune->Numerics("NUMERICS_SPIN");
  std::vector<spin::SpinSpec> specs;
  std::string note;
  const bool cascade = decay::HasCascadedDecay(lts.decaytree);
  if (lts.decaytree.size() != 2 || lts.upc_model != nullptr) {
    note = "Integrated spin metrics require a two-body central vertex without nuclear projection";
  } else if (!lts.process.SPINGEN) {
    note = "Production spin correlations are disabled";
  } else {
    if (lts.process.SPINDEC && lts.process.root_decay_mode != RootDecayMode::Isolated) {
      spin::SpinSpec pair;
      pair.name = cascade ? "intermediate_pair" : "final_pair";
      pair.type = cascade ? spin::SpinType::IntermediatePair : spin::SpinType::FinalPair;
      pair.frame = "Pair CM helicities";
      pair.measure = cascade
          ? "central pair spin state before decays, with final-state cuts applied"
          : "coherent full process, final phase space and fiducial cuts";
      for (std::size_t leg = 0; leg < 2; ++leg) {
        pair.pdg[leg] = lts.decaytree[leg].p.pdg;
        for (const double h : spin::FinalStateHelicities(lts.decaytree[leg].p, "spin metrics")) {
          pair.helicity[leg].push_back(static_cast<int>(std::llround(2.0 * h)));
        }
      }
      specs.push_back(std::move(pair));
    } else { note = "Final-pair metrics are unavailable with disabled or isolated spin decay"; }
    for (auto &[name, res] : lts.process.RESONANCES) {
      spin::SpinSpec parent;
      parent.name = name;
      parent.frame = model == ReggeProductionModel::MP ? lts.process.MP_FRAME : "CM";
      parent.measure = "production reference, decay operator omitted, selected final phase-space measure";
      parent.pdg[0] = res.p.pdg;
      for (int m = -res.p.spinX2; m <= res.p.spinX2; m += 2) { parent.helicity[0].push_back(m); }
      parent.helicity[1] = {0};
      res.spin_index = specs.size();
      specs.push_back(std::move(parent));
    }
  }
  const bool protons = setup.excitation == 0 && nuclear::IsProton(lts.beam1.pdg) && nuclear::IsProton(lts.beam2.pdg);
  if (!protons && !specs.empty()) { note += " Forward proton metrics require two intact protons"; }
  lts.qmetrics.Configure(std::move(specs), numerics.at("metric_groups"), numerics.at("metric_max_dimension"),
                             std::move(note), protons ? std::array<int, 2>{lts.beam1.pdg, lts.beam2.pdg}
                                                       : std::array<int, 2>{});
}

// Initialize Regge resonance and continuum branching structures
void MRegge::InitializeBranching(MProcessSetup &setup, ReggeProductionModel model, MReggeMode mode) {
  const BranchingProcessFlags flags = ClassifyBranchingProcess(setup);
  if (flags.model != model || model == ReggeProductionModel::TP || model == ReggeProductionModel::None) { throw std::invalid_argument("MRegge::InitializeBranching: production model mismatch"); }
  InitializePhysicalProcessState(setup, flags);
  if (flags.continuum) {
    if (setup.lts.decaytree.size() == 2) { BuildContinuum2Plan(setup, model); }
    if (setup.lts.decaytree.size() == 4) { BuildContinuum4Plan(setup, model); }
    if (setup.lts.decaytree.size() == 6) { BuildContinuum6Plan(setup, model); }
  }
  for (std::size_t leg = 0; leg < 2; ++leg) {
    if (!UsesPhotonProductionLeg(setup.lts, leg) || !gra::flux::SupportsPhotoDirection(setup.lts, static_cast<int>(leg + 1))) { continue; }
    gra::flux::ValidatePhotonEmitter(setup.lts, static_cast<int>(leg + 1), "Amplitude initialization PARAM_REGGE.PHOTON_VERTEX=" + setup.lts.process.PHOTON_VERTEX);
  }
  ValidateReggeExcitationOrderings(setup, flags);
  ValidateInitializedProcess(mode, setup.lts, setup.istate + "[" + setup.channel + "]");
  InitializeContinuumInterference(setup.lts);
  ConfigureReggeQMetrics(setup, model);
}

// Convert exchange PDG pairs into beam-side topology strings
// The topology strings are consumed by MPDG::TokenizeProcess
std::vector<std::string> MakeContinuumTopologyStrings(const std::vector<std::vector<int>> &production, int beam1_pdg, int beam2_pdg) {
  std::vector<std::string> channels;
  channels.reserve(production.size());
  const std::string beam1 = std::to_string(beam1_pdg);
  const std::string beam2 = std::to_string(beam2_pdg);

  for (const auto &entry : production) {
    if (entry.size() != 2) {
      throw std::invalid_argument(
          "Amplitude initialization: Continuum exchange "
          "channel must have two PDGs");
    }
    channels.push_back(std::to_string(entry[0]) + " > {" + beam1 + " " + beam1 + "} " + std::to_string(entry[1]) + " > {" + beam2 + " " + beam2 + "}");
  }
  return channels;
}

// Build the immutable two-body continuum production plan
void BuildContinuum2Plan(MProcessSetup &setup, const ReggeProductionModel model) {
  if (model != ReggeProductionModel::MP && model != ReggeProductionModel::XP && model != ReggeProductionModel::GP) { throw std::invalid_argument("Amplitude initialization: unsupported continuum production model"); }
  std::vector<int> refdecay;
  refdecay.reserve(setup.lts.decaytree.size());
  for (const auto &branch : setup.lts.decaytree) { refdecay.push_back(branch.p.pdg); }
  const auto              regge_param = regge::GetParam(RequireSetupModelCache(setup), refdecay, setup.lts.PDG);
  const regge::PairParam &pair        = regge::Pair(*regge_param, refdecay, model);

  setup.lts.process.CONT_PRODUCTION.clear();
  setup.lts.process.CONT_PRODUCTION.reserve(pair.channels.size());
  for (const auto &channel : pair.channels) { setup.lts.process.CONT_PRODUCTION.push_back({channel.first, channel.second}); }
  const std::vector<std::string> channels = MakeContinuumTopologyStrings(setup.lts.process.CONT_PRODUCTION, setup.lts.beam1.pdg, setup.lts.beam2.pdg);

  if (setup.lts.process.CONT_PRODUCTION.size() != channels.size()) {
    throw std::invalid_argument(
        "Amplitude initialization: Continuum "
        "exchange channel expansion mismatch");
  }

  setup.lts.process.CONT_PRODUCTIONTREE.clear();
  setup.lts.process.CONTINUUM_POLE.clear();
  setup.lts.process.CONTINUUM_GP.clear();

  std::string central_state = FormatPDGList(refdecay);
  central_state.front()     = '{';
  central_state.back()      = '}';

  for (const auto &k : indices(channels)) {
    std::cout << rang::fg::yellow << "Continuum production channel[" << k << "]: topology = " << channels[k] << ", central state = " << central_state << rang::fg::reset << std::endl;
    std::vector<MDecayBranch> tree;
    setup.lts.PDG.TokenizeProcess(channels[k], 0, tree);
    // GP uses analytic beam-side proton vertices instead of fixed-spin branch
    // helicity tables
    if (model != ReggeProductionModel::GP) {
      for (const auto &i : indices(tree)) {
        if (tree[i].legs.size() != 0) { setup.ProcessHelicityTree(tree[i], true, true); }
      }
    }
    setup.lts.process.CONT_PRODUCTIONTREE.push_back(tree);
    if (model == ReggeProductionModel::GP) {
      setup.lts.process.CONTINUUM_GP.push_back(gpom::PrepareContinuum(setup, setup.lts.process.CONT_PRODUCTIONTREE.back(), *regge_param));
    } else {
      auto pole = (model == ReggeProductionModel::MP ? mpom::PrepareContinuum : xpom::PrepareContinuum)(setup, tree);
      if (model == ReggeProductionModel::XP) {
        const std::vector<MParticle> upper_legs = {setup.lts.decaytree[0].p, setup.lts.decaytree[1].p};
        const std::vector<MParticle> lower_legs = {setup.lts.decaytree[1].p, setup.lts.decaytree[0].p};
        xpom::PrintPole("XP continuum upper pole operator", setup.lts.process.CONT_PRODUCTIONTREE.back()[0].p,
                            upper_legs, pole[0].Pole());
        xpom::PrintPole("XP continuum lower pole operator", setup.lts.process.CONT_PRODUCTIONTREE.back()[1].p,
                            lower_legs, pole[1].Pole());
      }
      setup.lts.process.CONTINUUM_POLE.push_back(std::move(pole));
    }
    if (model == ReggeProductionModel::GP) { gpom::PrintContinuum(setup.lts, setup.lts.process.CONT_PRODUCTION[k], setup.lts.process.CONT_PRODUCTIONTREE[k]); }
  }
}

// Build the immutable four-body continuum production plan
void BuildContinuum4Plan(MProcessSetup &setup, const ReggeProductionModel model) { BuildContinuumLadderPlan(setup, model, 4); }

// Build the immutable six-body continuum production plan
void BuildContinuum6Plan(MProcessSetup &setup, const ReggeProductionModel model) { BuildContinuumLadderPlan(setup, model, 6); }

}  // namespace gra
