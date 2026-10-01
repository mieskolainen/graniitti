// Regge process interfaces and common exchange kernels
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

using gra::math::msqrt;

namespace gra {

// Compute the beam particle attached to one forward leg
//
const gra::MParticle &ReggeBeamParticleForLeg(const gra::LORENTZSCALAR &lts, int leg) {
  if (leg == 1) { return lts.beam1; }
  if (leg == 2) { return lts.beam2; }
  throw std::invalid_argument("MRegge: beam leg should be 1 or 2");
}

namespace {

// Validate one generated Regge subenergy and momentum transfer
void ValidateReggeEventKinematics(const double s, const double t, const std::string &context) {
  if (!std::isfinite(s) || !(s > 0.0) || !std::isfinite(t)) { throw AmplitudeFailure(context + ": invalid generated s or t"); }
}

// Compute true when every top-level branch is a stable direct final state
bool ReggeDirectFinalState(const std::vector<gra::MDecayBranch> &tree) {
  return std::all_of(tree.begin(), tree.end(), [](const auto &branch) { return branch.legs.empty(); });
}

// Match the decay-tree multiplicities accepted by one Regge amplitude mode
bool ReggeProcessAccepts(const std::vector<gra::MDecayBranch> &tree, gra::MReggeMode mode) {
  switch (mode) {
    case gra::MReggeMode::Generic:
    case gra::MReggeMode::Soft:
      return true;
    case gra::MReggeMode::Resonance:
      return tree.size() >= 2;
    case gra::MReggeMode::ContinuumTwoBody:
    case gra::MReggeMode::ResonanceContinuumTwoBody:
      return tree.size() == 2;
    case gra::MReggeMode::ContinuumTwoFourSixBody:
      if (tree.size() == 2) { return true; }
      return (tree.size() == 4 || tree.size() == 6) && ReggeDirectFinalState(tree);
  }
  return false;
}

// Compute a compact final-state pattern for one Regge amplitude mode
std::string ReggeProcessPattern(gra::MReggeMode mode) {
  switch (mode) {
    case gra::MReggeMode::Generic:
    case gra::MReggeMode::Soft:
      return "event Regge final state";
    case gra::MReggeMode::Resonance:
      return "resonance decays";
    case gra::MReggeMode::ContinuumTwoBody:
    case gra::MReggeMode::ResonanceContinuumTwoBody:
      return "2 hadrons";
    case gra::MReggeMode::ContinuumTwoFourSixBody:
      return "2/4/6 hadrons and pair decays";
  }
  throw std::invalid_argument("MRegge::ProcessDefinitionFor: unknown process mode");
}

// Require complete two-body continuum caches built from the selected model
void ValidateReggeContinuumCache(const gra::LORENTZSCALAR &lts, const std::string &process_name) {
  if (lts.process.CONT_PRODUCTION.empty()) { throw std::invalid_argument("MRegge process " + process_name + " has no model-supported continuum exchange channel"); }
  for (const auto &channel : lts.process.CONT_PRODUCTION) {
    if (channel.size() != 2) { throw std::invalid_argument("MRegge process " + process_name + " has a malformed continuum exchange channel"); }
  }
  const bool        gp              = !lts.process.CONTINUUM_GP.empty();
  const std::size_t vertex_channels = gp ? lts.process.CONTINUUM_GP.size() : lts.process.CONTINUUM_POLE.size();
  if (lts.process.CONT_PRODUCTIONTREE.size() != lts.process.CONT_PRODUCTION.size() || vertex_channels != lts.process.CONT_PRODUCTION.size()) {
    throw std::invalid_argument("MRegge process " + process_name + " has inconsistent continuum channel caches");
  }
  for (const auto &i : indices(lts.process.CONT_PRODUCTION)) {
    const std::size_t vertices = gp ? lts.process.CONTINUUM_GP[i].size() : lts.process.CONTINUUM_POLE[i].size();
    if (lts.process.CONT_PRODUCTIONTREE[i].size() != 2 || vertices != 4) { throw std::invalid_argument("MRegge process " + process_name + " has malformed continuum production or helicity data"); }
    if (!gp && std::any_of(lts.process.CONTINUUM_POLE[i].cbegin(), lts.process.CONTINUUM_POLE[i].cend(),
                           [](const spin::PoleResidue &residue) { return !residue.Ready(); })) {
      throw std::invalid_argument("MRegge process " + process_name + " has an unprepared continuum pole residue");
    }
  }
}

}  // namespace

// Validate all Regge caches and steering fixed during initialization
void MRegge::ValidateInitializedProcess(MReggeMode mode, const LORENTZSCALAR &lts, const std::string &process_name) {
  if (lts.decaytree.empty()) { return; }

  if (mode == gra::MReggeMode::Resonance) {
    if (lts.process.RESONANCES.empty()) { throw std::invalid_argument("MRegge process " + process_name + " requires an active resonance"); }
    return;
  }

  const bool two_body_continuum = mode == gra::MReggeMode::ContinuumTwoBody || mode == gra::MReggeMode::ResonanceContinuumTwoBody || (mode == gra::MReggeMode::ContinuumTwoFourSixBody && lts.decaytree.size() == 2);
  if (two_body_continuum) {
    ValidateReggeContinuumCache(lts, process_name);
    if (mode == gra::MReggeMode::ResonanceContinuumTwoBody && lts.process.RESONANCES.empty()) { throw std::invalid_argument("MRegge process " + process_name + " requires an active resonance"); }
    return;
  }

  if (mode == gra::MReggeMode::ContinuumTwoFourSixBody && (lts.decaytree.size() == 4 || lts.decaytree.size() == 6)) {
    if (!lts.process.CONT_LADDER_PERMUTATIONS.empty()) { regge::CheckLadderPermutations(lts.process.CONT_LADDER_PERMUTATIONS, lts.decaytree.size(), "MRegge process " + process_name); }
    if (!lts.process.MULTIREGGE_TOPOLOGIES.empty()) { regge::CheckMultiReggeTopologies(lts.process.MULTIREGGE_TOPOLOGIES, lts.decaytree.size(), "MRegge process " + process_name); }
  }
}

// Build one immutable Regge process definition for a selected process mode
std::shared_ptr<const amplitude::ProcessDefinition> MRegge::ProcessDefinitionFor(MReggeMode mode, const std::string &process_name) {
  return std::make_shared<amplitude::AnalyticProcess>(
      "MREGGE", process_name, ReggeProcessPattern(mode), DecayStructure{}, [mode](const std::vector<MDecayBranch> &tree) { return ReggeProcessAccepts(tree, mode); },
      [mode](const LORENTZSCALAR &lts) {
        if (mode == MReggeMode::Generic || mode == MReggeMode::Soft) { return DecayStructure{}; }
        if (mode == MReggeMode::ContinuumTwoFourSixBody && lts.decaytree.size() > 2) {
          return DecayStructure{DecayType::Full};
        }
        return decay::JacobWickStructure(lts);
      });
}

// Initialize immutable Regge parameters before worker process copies
void MRegge::InitializeParameters(const MProcessSetup &setup) {
  if (!setup.model_tune) { throw std::invalid_argument("MRegge initialization requires a model tune"); }
  MModelCache &cache = RequireModelCache(setup.lts.model_cache, setup.model_tune, "MRegge initialization");
  (void)regge::GetParam(cache, regge::FinalPDGs(setup.lts), setup.lts.PDG);
  (void)GetReggeNumerics(cache);
}

namespace {

// Build the reusable polar transform from the Regge production quadrature
std::unique_ptr<math::MPolarFourier> BuildParallelFourierTransform(const MReggeNumerics &numerics) {
  const auto &rule = numerics.ParallelRule();
  if (rule.kt.empty() || rule.radial_weight.size() != rule.kt.size() || rule.measure_weight.size_row() != rule.kt.size() || rule.measure_weight.size_col() == 0 || rule.kt_x.size_row() != rule.kt.size() ||
      rule.kt_x.size_col() != rule.measure_weight.size_col() || rule.kt_y.size_row() != rule.kt.size() || rule.kt_y.size_col() != rule.measure_weight.size_col()) {
    throw std::invalid_argument("MRegge parallel Fourier transform has an invalid production rule");
  }

  math::PolarMeasureRule momentum;
  momentum.node           = rule.kt;
  momentum.measure_weight = rule.radial_weight;

  const std::size_t azimuth_count = rule.measure_weight.size_col();
  if (azimuth_count > std::numeric_limits<unsigned int>::max()) { throw std::overflow_error("MRegge parallel Fourier transform rule is too large"); }
  const std::size_t max_harmonic = (azimuth_count - 1) / 2;
  return std::make_unique<math::MPolarFourier>(math::BuildPolarFourierTransform(std::move(momentum), static_cast<unsigned int>(azimuth_count), numerics.ParallelImpactCount(), numerics.ParallelImpactMax(), max_harmonic));
}

// Compute whether the selected topology bank contains simultaneous ladders
bool RequiresParallelProduction(const LORENTZSCALAR &lts, const regge::Param &param) {
  const std::size_t central_count = lts.decaytree.size();
  if (central_count != 4 && central_count != 6) { return false; }
  const auto configured = param.con.multiregge_topologies.find(static_cast<int>(central_count));
  if (configured == param.con.multiregge_topologies.end()) { return false; }
  const auto &topologies = lts.process.MULTIREGGE_TOPOLOGIES.empty() ? configured->second : lts.process.MULTIREGGE_TOPOLOGIES;
  return std::any_of(topologies.cbegin(), topologies.cend(), [](const auto &topology) { return topology.size() > 1; });
}

// Compute whether the selected six-body topology bank needs a Fourier product
bool RequiresParallelFourier(const LORENTZSCALAR &lts, const regge::Param &param) {
  if (lts.decaytree.size() != 6) { return false; }
  const auto configured = param.con.multiregge_topologies.find(6);
  if (configured == param.con.multiregge_topologies.end()) { return false; }
  const auto &topologies = lts.process.MULTIREGGE_TOPOLOGIES.empty() ? configured->second : lts.process.MULTIREGGE_TOPOLOGIES;
  return std::any_of(topologies.cbegin(), topologies.cend(), [](const auto &topology) { return topology.size() == 3; });
}

}  // namespace

// Class constructor
MRegge::MRegge(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune_handle, std::shared_ptr<const amplitude::ProcessDefinition> definition)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(RequireModelCache(lts.model_cache, model_tune_handle, "MRegge").TunePtr()),
      soft_model(model_tune->Soft()),
      param_handle(regge::GetParam(*lts.model_cache, regge::FinalPDGs(lts), lts.PDG)),
      param(*param_handle),
      continuum_plan(lts.process.CONT_LADDER_POLE.empty() || lts.process.REGGE_MODEL == ReggeProductionModel::None
                         ? regge::ContinuumPlan{}
                         : lts.decaytree.size() == 4 ? regge::BuildContinuum4Plan(lts, param, lts.process.REGGE_MODEL)
                                                     : lts.decaytree.size() == 6 ? regge::BuildContinuum6Plan(lts, param, lts.process.REGGE_MODEL) : regge::ContinuumPlan{}),
      parallel_required(RequiresParallelProduction(lts, param)),
      regge_numerics(GetReggeNumerics(*lts.model_cache)),
      parallel_fourier_required(RequiresParallelFourier(lts, param)),
      triple_pomeron_root(soft_model->TriplePomeronCouplingRoot(param.exchanges.at(param.pomeron_trajectory).soft_exchange)) {
  if (lts.process.CONT_TU_SIGN.empty()) { InitializeContinuumInterference(lts); }
  if (!parallel_required) { return; }
  if (parallel_fourier_required) { parallel_fourier = BuildParallelFourierTransform(*regge_numerics); }
}

// Compute one mapped SOFT propagator without beam vertices
std::complex<double> MRegge::PropagatorForExchange(const double s, const double t, const SoftExchangeId exchange) const {
  ValidateReggeEventKinematics(s, t, "MRegge::PropagatorForExchange");
  return PreparedPropagator(s, t, exchange);
}

// Compute one prepared mapped SOFT propagator without repeated lookup
std::complex<double> MRegge::PreparedPropagator(const double s, const double t, const SoftExchangeId exchange) const {
  const auto                &soft  = soft_model->Exchange(exchange);
  const double               alpha = soft_model->Alpha(exchange, t);
  const std::complex<double> pole =
      soft.Eta() == EtaMode::Raw ? regge::Pole(s, param.s0, {alpha, 0.0}, soft.Signature(), regge::Rim::Lower) : regge::EtaFactor(alpha, soft.Alpha0(), soft.Signature(), soft.Eta()) * std::pow(s / param.s0, alpha);

  return static_cast<double>(soft.ResidueSign()) * pole;
}

// Compute the propagator of the explicitly mapped Pomeron trajectory
std::complex<double> MRegge::PomeronKernel(const double s, const double t) const { return PropagatorForExchange(s, t, param.exchanges.at(param.pomeron_trajectory).soft_exchange); }

// Compute the propagator of one central trajectory PDG alias
std::complex<double> MRegge::ExchangeKernel(const double s, const double t, const int exchange_pdg) const {
  const std::size_t trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  return PropagatorForExchange(s, t, param.exchanges.at(trajectory).soft_exchange);
}

}  // namespace gra
