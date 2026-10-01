// Two-body Regge resonance and continuum amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeCEP.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {
namespace {

// Count unstable lines connecting reduced vertices in one branch
std::size_t CascadeLineCount(const MDecayBranch &branch) {
  if (branch.legs.empty()) { return 0; }
  std::size_t lines = 1;
  for (const auto &leg : branch.legs) { lines += CascadeLineCount(leg); }
  return lines;
}

// Count unstable lines connecting reduced vertices in one decay tree
std::size_t CascadeLineCount(const std::vector<MDecayBranch> &tree) {
  std::size_t lines = 0;
  for (const auto &branch : tree) { lines += CascadeLineCount(branch); }
  return lines;
}

// Store event-local factors shared by ordered resonance production rows
struct ResonanceProductionContext {
  int                  central_pdg           = 0;
  std::complex<double> subenergy_suppression = 0.0;
  std::complex<double> line_shape            = 0.0;
  double               production_ff         = 1.0;
  regge::FFParam       ff_transfer;
};

// Build the pole-normalized two-body LS decay profile
gra::MReggeDecayState BuildResonanceDecayKinematics(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &resonance) {
  gra::MReggeDecayState out;
  if (lts.decaytree.size() != 2 || !(lts.m2 > 0.0) || !(resonance.p.mass > 0.0)) { return out; }

  const double mass_a = lts.decaytree[0].p4.M();
  const double mass_b = lts.decaytree[1].p4.M();
  const double q0     = gra::kinematics::DecayMomentum(resonance.p.mass, mass_a, mass_b);
  if (!(q0 > 0.0)) { return out; }

  const double mass = msqrt(lts.m2);
  const double q    = gra::kinematics::DecayMomentum(mass, mass_a, mass_b);
  out.q_ratio       = q > 0.0 ? q / q0 : 0.0;
  if (lts.process.DECAY_BARRIER && resonance.hel_decay.ls_components.empty()) {
    out.q_ratio = 1.0;
    return out;
  }

  const double intensity = gra::spin::DecayLSIntensity(resonance.hel_decay, out.q_ratio, lts.process.DECAY_BARRIER);

  // The running profile is W(s) = M0/sqrt(s) times q/q0 times the orthogonal LS intensity
  out.running_profile = resonance.p.mass / mass * out.q_ratio * intensity;
  out.valid           = true;
  return out;
}

// Build the MP, XP or GP decay matrix with the selected LS barrier policy
MMatrix<std::complex<double>> BuildReggeResonanceDecayMatrix(gra::LORENTZSCALAR &lts, const gra::PARAM_RES &resonance, const gra::MReggeDecayState &decay, const std::string &frame, const bool cascade = true) {
  if (!lts.process.DECAY_BARRIER || resonance.hel_decay.ls_components.empty() || lts.process.root_decay_mode == gra::RootDecayMode::Isolated || lts.decaytree.size() != 2) { return gra::spin::ResonanceDecayMatrix(lts, resonance, resonance.hel_decay, frame, cascade); }

  if (!lts.process.SPINDEC) { return gra::spin::ResonanceDecayMatrix(lts, resonance, frame) * std::sqrt(gra::spin::DecayLSIntensity(resonance.hel_decay, decay.q_ratio, true)); }

  gra::HELMatrix scaled = resonance.hel_decay;
  scaled.T              = gra::spin::DecayLSHelicityMatrix(resonance.hel_decay, decay.q_ratio, true);
  return gra::spin::ResonanceDecayMatrix(lts, resonance, scaled, frame, cascade);
}

// Count the strong exchange legs in one ordered two-leg production row
std::size_t HadronicLegCount(const int upper_pdg, const int lower_pdg) { return static_cast<std::size_t>(upper_pdg != PDG::PDG_gamma) + static_cast<std::size_t>(lower_pdg != PDG::PDG_gamma); }

// Compute the selected MP, XP or GP resonance model settings
const RES_PRODUCTION_MODEL &ResonanceModelForm(const PARAM_RES &resonance, const ReggeProductionModel model) {
  if (model == ReggeProductionModel::MP) { return resonance.MP; }
  if (model == ReggeProductionModel::XP) { return resonance.XP; }
  if (model == ReggeProductionModel::GP) { return resonance.GP; }
  throw std::invalid_argument("MRegge resonance form factors require MP, XP or GP");
}

// Collect the event-local context shared by all ordered resonance rows
ResonanceProductionContext BuildResonanceProductionContext(gra::LORENTZSCALAR &lts, const gra::PARAM_RES &resonance, const regge::Param &param, const ReggeProductionModel model, const gra::MReggeDecayState &decay) {
  regge::CheckMass2(lts.m2, "MRegge resonance production");
  ResonanceProductionContext out;
  const RES_MODEL_FORM &form  = ResonanceModelForm(resonance, model);
  const auto           &prod  = form.ff_prod;
  const auto           &decay_ff = resonance.hel_decay.ff_decay;
  out.ff_transfer             = form.ff_transfer;
  out.central_pdg             = resonance.p.pdg;
  out.subenergy_suppression   = std::pow(param.s0 / lts.m2, param.omega.at(model));
  out.line_shape              = gra::resonance::LineShape(lts.m2, resonance, decay.valid ? decay.running_profile : 1.0) * regge::MassFF(lts.m2, pow2(resonance.p.mass), decay_ff);
  out.production_ff           = regge::MassFF(lts.m2, pow2(resonance.p.mass), prod);
  return out;
}

// One forward Regge source before the two-leg Cartesian product
std::vector<PairSource> ProductionSources(const gra::LORENTZSCALAR &lts, const gra::MRegge &regge, const regge::Param &param, const std::vector<gra::MDecayBranch> &tree, const ResonanceProductionContext &production,
                                          ReggeSourceCache &source_cache) {
  const int                  up            = tree[0].p.pdg;
  const int                  dn            = tree[1].p.pdg;
  const bool                 upper_gamma   = up == PDG::PDG_gamma;
  const bool                 lower_gamma   = dn == PDG::PDG_gamma;
  const bool                 one_gamma     = upper_gamma != lower_gamma;
  const bool                 photo_pomeron = one_gamma && regge::TrajectoryIndex(param, upper_gamma ? dn : up) == param.pomeron_trajectory;
  const double               upper_s       = source_cache.Subenergy(ForwardBeamLeg::Upper);
  const double               lower_s       = source_cache.Subenergy(ForwardBeamLeg::Lower);
  const auto                &upper         = source_cache.Sources(ForwardBeamLeg::Upper, up, upper_s, photo_pomeron && !upper_gamma, production.central_pdg);
  const auto                &lower         = source_cache.Sources(ForwardBeamLeg::Lower, dn, lower_s, photo_pomeron && !lower_gamma, production.central_pdg);
  const std::complex<double> scale         = !upper_gamma && !lower_gamma ? production.subenergy_suppression : std::complex<double>(1.0, 0.0);
  return PairSources(upper, lower, regge.SoftModelHandle()->GoodWalker(), ProtonSpinLayout(lts), scale);
}

// Accumulate ordered resonance production channels coherently by sector
void SumRes(gra::LORENTZSCALAR &lts, const gra::MRegge &regge, const regge::Param &param, const std::vector<MMatrix<std::complex<double>>> &channel_prod, const std::vector<gra::RES_PRODUCTION> &production_channels,
            const MReggeDecayState &decay, const ResonanceProductionContext &production, std::complex<double> common_pref, const std::size_t coherence_group, ReggeSourceCache &source_cache,
            const std::optional<std::size_t> channel_filter, const std::optional<std::size_t> spin_index) {

  const auto pair_index = lts.qmetrics.Active() ? lts.qmetrics.Find("intermediate_pair") : std::nullopt;
  for (const auto &i : indices(channel_prod)) {
    if (channel_filter.has_value() && i != *channel_filter) { continue; }
    const auto &tree = production_channels[i].tree;
    const int         upper_pdg     = tree[0].p.pdg;
    const int         lower_pdg     = tree[1].p.pdg;
    const std::size_t hadronic_legs = HadronicLegCount(upper_pdg, lower_pdg);
    double            form_factor   = production.production_ff;
    if (hadronic_legs == 2) { form_factor *= regge::TransferFF(lts.t1, production.ff_transfer) * regge::TransferFF(lts.t2, production.ff_transfer); }
    const auto central = channel_prod[i].MultiplyScaled(decay.matrix, common_pref * form_factor);
    const auto pairs   = ProductionSources(lts, regge, param, tree, production, source_cache);
    AddGoodWalker(lts, pairs, central, coherence_group, regge.SoftModelHandle());
    if (pair_index) {
      AddGoodWalker(lts, pairs, channel_prod[i].MultiplyScaled(decay.pair, common_pref * form_factor),
                    coherence_group, regge.SoftModelHandle(), *pair_index);
    }
    if (lts.qmetrics.Active() && spin_index) {
      AddGoodWalker(lts, pairs, channel_prod[i] * (common_pref * form_factor), coherence_group,
                    regge.SoftModelHandle(), *spin_index);
    }
  }
}

// Identify ordered particle-antiparticle two-body final states
bool ParticleAntiparticlePair(const gra::MDecayBranch &a, const gra::MDecayBranch &b) { return a.p.pdg == -b.p.pdg; }

// Compute the crossed Bose/Fermi sign of a particle-antiparticle pair
// [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:1801.03902, Eqs. (2.10)--(2.13)]
int ParticleAntiparticleTUSign(const gra::LORENTZSCALAR &lts) {
  if (lts.decaytree.size() != 2 || !ParticleAntiparticlePair(lts.decaytree[0], lts.decaytree[1])) {
    throw std::invalid_argument(
        "MRegge::Con: particle-antiparticle t/u sign requested for another "
        "final state");
  }
  if (lts.decaytree[0].p.spinX2 != lts.decaytree[1].p.spinX2) { throw std::invalid_argument("MRegge::Con: particle-antiparticle pair has unequal spins"); }
  return (lts.decaytree[0].p.spinX2 % 2 == 0) ? 1 : -1;
}

// Resolve C for a two-body final state with individually defined C eigenvalues
bool FixedTwoBodyCParity(const gra::LORENTZSCALAR &lts, int &C_final) {
  if (lts.decaytree.size() != 2) { return false; }
  const auto &a = lts.decaytree[0];
  const auto &b = lts.decaytree[1];
  if (ParticleAntiparticlePair(a, b)) { return false; }
  if (!spin::IsValidCParity(a.p.C) || !spin::IsValidCParity(b.p.C)) { return false; }
  C_final = a.p.C * b.p.C;
  return true;
}

// Compute true for identical top-level continuum particles, before later cascade
// decays
bool IdenticalTopLevelPair(const gra::LORENTZSCALAR &lts) { return lts.decaytree.size() == 2 && lts.decaytree[0].p.pdg == lts.decaytree[1].p.pdg; }

// Compute the Bose/Fermi sign required by top-level identical continuum
// particles
int TopLevelStatisticsSign(const gra::LORENTZSCALAR &lts) {
  if (!IdenticalTopLevelPair(lts)) {
    throw std::invalid_argument(
        "MRegge::Con: top-level statistics sign "
        "requested for non-identical states");
  }
  return (lts.decaytree[0].p.spinX2 % 2 == 0) ? 1 : -1;
}

// Compute whether an exact top-level Bose/Fermi exchange sign exists
bool IdenticalTopLevelStatisticsProjector(const gra::LORENTZSCALAR &lts, int &sign) {
  if (!IdenticalTopLevelPair(lts)) { return false; }
  sign = TopLevelStatisticsSign(lts);
  return true;
}

// Format one top-level continuum final-state pair for diagnostics
std::string ContinuumFinalStateLabel(const gra::LORENTZSCALAR &lts) {
  if (lts.decaytree.size() != 2) { return "<non two-body>"; }
  return std::to_string(lts.decaytree[0].p.pdg) + " " + std::to_string(lts.decaytree[1].p.pdg);
}

// Resolve the continuum sign once, with explicit steering overriding the statistics choice
// Fixed-C channel selection is independent of the relative t/u sign
double ResolveContinuumTUSign(const gra::LORENTZSCALAR &lts, const regge::Param &param, const std::vector<int> &exchange_pair) {
  const auto &mode = lts.process.TU_SIGN;
  if (mode != "auto" && mode != "positive" && mode != "negative") {
    throw std::invalid_argument("MRegge::Con: Unknown PARAM_REGGE.TU_SIGN = " + mode);
  }
  if (exchange_pair.size() == 2 && lts.decaytree.size() == 2) {
    const int C_exchange = regge::PairCParity(param, exchange_pair[0], exchange_pair[1]);
    int C_final = 0;
    if (FixedTwoBodyCParity(lts, C_final) && C_exchange != C_final) {
      throw std::invalid_argument("MRegge::Con: PARAM_REGGE.TU_SIGN=" + mode + " cannot make C=" + std::to_string(C_exchange) + " exchange produce fixed-C final state [" + ContinuumFinalStateLabel(lts) + "] with C=" + std::to_string(C_final));
    }
  }
  if (mode == "positive") { return 1.0; }
  if (mode == "negative") { return -1.0; }

  int statistics_sign = 1;
  if (IdenticalTopLevelStatisticsProjector(lts, statistics_sign)) { return statistics_sign; }
  if (lts.decaytree.size() == 2 && ParticleAntiparticlePair(lts.decaytree[0], lts.decaytree[1])) { return ParticleAntiparticleTUSign(lts); }
  return 1.0;
}

// Common pair profile and scalar factors for the continuum t and u channels
struct ContinuumTUFactors {
  std::vector<PairSource> pair;
  std::complex<double>    t = 0.0;
  std::complex<double>    u = 0.0;
};

// Select the pre-resolved continuum vertex for one production-plan channel
const regge::VertexParam &ContinuumVertex(const regge::PairParam &entry, const std::size_t channel,
                                         const int first, const int second, regge::VertexParam &photon) {
  if (channel < entry.channels.size() && entry.channels[channel].first == first &&
      entry.channels[channel].second == second) {
    return entry.channels[channel];
  }
  const auto match = std::find_if(entry.channels.cbegin(), entry.channels.cend(), [&](const auto &vertex) {
    return vertex.first == first && vertex.second == second;
  });
  if (match != entry.channels.cend()) { return *match; }

  photon.first  = first;
  photon.second = second;
  for (const std::size_t leg : {std::size_t{0}, std::size_t{1}}) {
    const int exchange = leg == 0 ? first : second;
    const auto &form = entry.forms.at(exchange);
    photon.transfer[leg] = exchange == PDG::PDG_gamma ? regge::FFParam{} : form.transfer;
    photon.offshell[leg] = form.offshell;
  }
  return photon;
}

// Compute the pair-space continuum Regge factors
ContinuumTUFactors ME4ContinuumTU(const gra::MRegge &regge, gra::LORENTZSCALAR &lts, const regge::Param &param, const std::vector<int> &exchange_channel, const std::size_t channel,
                                 const ReggeProductionModel model, ReggeSourceCache &source_cache) {

  const double               M2_t      = pow2(lts.decaytree[0].p.mass);
  const double               M2_u      = pow2(lts.decaytree[1].p.mass);
  const int                  A         = exchange_channel[0];
  const int                  B         = exchange_channel[1];
  const std::vector<int>     pair_pdgs = {lts.decaytree[0].p.pdg, lts.decaytree[1].p.pdg};
  const auto                &entry     = regge::Pair(param, pair_pdgs, model);
  regge::VertexParam         photon_vertex;
  const auto                &vertex    = ContinuumVertex(entry, channel, A, B, photon_vertex);
  const std::complex<double> t_vertex = regge::PairExchange(param, entry, vertex, lts.decaytree[0], lts.decaytree[1], lts.t_hat, M2_t, lts.t1, lts.t2);
  const std::complex<double> u_vertex = regge::PairExchange(param, entry, vertex, lts.decaytree[1], lts.decaytree[0], lts.u_hat, M2_u, lts.t1, lts.t2);

  const auto        &upper       = source_cache.Profile(ForwardBeamLeg::Upper, A);
  const auto        &lower       = source_cache.Profile(ForwardBeamLeg::Lower, B);
  const auto        &good_walker = regge.SoftModelHandle()->GoodWalker();
  ContinuumTUFactors out;
  out.pair = PairSources(upper, lower, good_walker, ProtonSpinLayout(lts), 1.0);
  // Each two-body t or u graph joins two reduced vertices through one line
  // The resulting SewingSign(1) = -1 is common to both graph orderings
  const double sewing = SewingSign(1);
  out.t               = sewing * source_cache.Kernel(ForwardBeamLeg::Upper, A, lts.ss[1][3], false, 0) * source_cache.Kernel(ForwardBeamLeg::Lower, B, lts.ss[2][4], false, 0) * t_vertex;
  out.u               = sewing * source_cache.Kernel(ForwardBeamLeg::Upper, A, lts.ss[1][4], false, 0) * source_cache.Kernel(ForwardBeamLeg::Lower, B, lts.ss[2][3], false, 0) * u_vertex;
  return out;
}

// Sum continuum channels with common propagators, statistics, decay and Good-Walker factors
template <typename Contract>
void SumCon(const MRegge &regge, LORENTZSCALAR &lts, const regge::Param &param,
            ReggeProductionModel model, ReggeSourceCache &source_cache, Contract contract) {
  const auto strong_veto = regge::Veto(regge::Pair(param, {lts.decaytree[0].p.pdg, lts.decaytree[1].p.pdg}, model).veto, msqrt(lts.s_hat));
  const double sewing = SewingSign(CascadeLineCount(lts.decaytree));
  for (const auto i : indices(lts.process.CONT_PRODUCTION)) {
    const auto &exchange = lts.process.CONT_PRODUCTION[i];
    const auto factors = ME4ContinuumTU(regge, lts, param, exchange, i, model, source_cache);
    const double veto = HadronicLegCount(exchange[0], exchange[1]) == 2 ? strong_veto : 1.0;
    AddGoodWalker(lts, factors.pair, contract(i, factors, sewing * veto), 0, regge.SoftModelHandle());
  }
}

}  // namespace

// Cache the continuum t/u interference signs fixed by model steering
void MRegge::InitializeContinuumInterference(LORENTZSCALAR &lts) {
  lts.process.CONT_TU_SIGN.clear();
  if (lts.decaytree.size() != 2 || lts.process.CONT_PRODUCTION.empty()) { return; }
  if (lts.model_cache == nullptr) { throw std::invalid_argument("MRegge::InitializeContinuumInterference: missing model cache"); }
  const auto param_handle = regge::GetParam(*lts.model_cache, regge::FinalPDGs(lts), lts.PDG);
  lts.process.CONT_TU_SIGN.reserve(lts.process.CONT_PRODUCTION.size());
  for (const auto &exchange_pair : lts.process.CONT_PRODUCTION) { lts.process.CONT_TU_SIGN.push_back(ResolveContinuumTUSign(lts, *param_handle, exchange_pair)); }
}

namespace reggeamp {

// Compute the root decay frame sewn to one Regge production basis
std::string RootFrame(const gra::LORENTZSCALAR &lts, const ReggeProductionModel model) {
  if (model == ReggeProductionModel::MP) { return lts.process.MP_FRAME; }
  if (model == ReggeProductionModel::XP || model == ReggeProductionModel::GP) { return "CM"; }
  throw std::invalid_argument("MRegge::Amp2: unsupported root decay frame");
}

// Resolve one transfer-independent decay matrix for the active event
const MReggeDecayState *Decay(LORENTZSCALAR &lts, const std::string &name, const PARAM_RES &resonance, const std::string &root_frame, MReggeDecayState &local) {
  if (lts.screening.active && lts.amplitude.central.Active()) { return &lts.amplitude.regge_decay.Get(lts.amplitude.central, name, "MRegge::EvalRes"); }
  local        = BuildResonanceDecayKinematics(lts, resonance);
  // Keep external helicities in the same CM basis as the continuum
  local.matrix = BuildReggeResonanceDecayMatrix(lts, resonance, local, "CM");
  if (lts.qmetrics.Active() && lts.qmetrics.Find("intermediate_pair")) {
    local.pair = BuildReggeResonanceDecayMatrix(lts, resonance, local, "CM", false);
  }
  if (root_frame != "CM" && lts.process.SPINDEC && lts.process.root_decay_mode != RootDecayMode::Isolated && lts.decaytree.size() == 2) {
    const auto rotation = spin::ProductionRotation(lts, root_frame, 0.5 * resonance.p.spinX2).Conj();
    local.matrix = rotation * local.matrix;
    if (!local.pair.isEmpty()) { local.pair = rotation * local.pair; }
  }
  if (lts.amplitude.central.Active()) { return &lts.amplitude.regge_decay.Store(lts.amplitude.central, name, std::move(local), "MRegge::EvalRes"); }
  return &local;
}

// Evaluate one selected resonance implementation with shared GP work
void EvalResonance(const gra::MRegge &regge, gra::LORENTZSCALAR &lts, const gra::regge::Param &param, const std::string &name, gra::PARAM_RES &resonance, gra::ReggeProductionModel model, gra::gpom::AmpCache &gp_amp,
                   ReggeSourceCache &source_cache, const std::size_t coherence_group, const std::optional<std::size_t> channel_filter, const MReggeDecayState *prepared_decay, const std::array<ForwardLegState, 2> *forward_state) {

  const std::string                root_frame = RootFrame(lts, model);
  gra::MReggeDecayState            local_decay;
  const gra::MReggeDecayState     *decay      = prepared_decay != nullptr ? prepared_decay : Decay(lts, name, resonance, root_frame, local_decay);
  const ResonanceProductionContext production = BuildResonanceProductionContext(lts, resonance, param, model, *decay);
  // The direct production-decay resonance line gives SewingSign(1) = -1
  // Every additional unstable cascade line supplies one further minus sign
  const double sewing = SewingSign(1 + CascadeLineCount(lts.decaytree));
  const std::complex<double> decay_coupling =
      param.use_zeta.at(model) ? resonance.hel_decay.g_decay
                     : std::complex<double>(std::abs(resonance.hel_decay.g_decay), 0.0);
  const std::complex<double> common =
      sewing * decay_coupling * production.line_shape * std::polar(1.0, ResonanceModelForm(resonance, model).phi);

  const auto prod = model == ReggeProductionModel::MP
                        ? mpom::Resonance(lts, resonance, param.s0, forward_state)
                    : model == ReggeProductionModel::XP ? xpom::Resonance(lts, resonance, param.s0, forward_state)
                                                        : gpom::Resonance(lts, param, resonance, &gp_amp);
  SumRes(lts, regge, param, prod, resonance.production, *decay, production, common, coherence_group,
         source_cache, channel_filter, resonance.spin_index);
}

}  // namespace reggeamp

namespace {

// Evaluate one selected continuum implementation
void EvalContinuum2(const gra::MRegge &regge, gra::LORENTZSCALAR &lts, const gra::regge::Param &param, gra::ReggeProductionModel model, ReggeSourceCache &source_cache, gra::gpom::AmpCache &gp_amp) {
  const bool analytic       = model == ReggeProductionModel::GP;
  const auto decay          = gra::spin::ContinuumDecayMatrix(lts, "CM");
  const bool compact_scalar = !lts.process.SPINGEN && !analytic && lts.hamp.metadata.spin_basis == ScreeningSpinBasis::ProtonIdentity && lts.hamp.metadata.spin_rows == 1;
  if (compact_scalar) {
    SumCon(regge, lts, param, model, source_cache, [&](std::size_t i, const ContinuumTUFactors &f, double scale) {
      const auto pole = rspin::ContinuumPoleScales(lts.process.CONTINUUM_POLE[i]);
      return decay * ((pole.t * f.t + lts.process.CONT_TU_SIGN[i] * pole.u * f.u) * scale);
    });
    return;
  }
  const auto pair_index = lts.qmetrics.Active() ? lts.qmetrics.Find("intermediate_pair") : std::nullopt;
  const auto channels = model == ReggeProductionModel::GP ? gpom::Continuum(lts, param, &gp_amp) : model == ReggeProductionModel::MP ? mpom::Continuum(lts, param.s0) : xpom::Continuum(lts, param.s0);
  SumCon(regge, lts, param, model, source_cache, [&](std::size_t i, const ContinuumTUFactors &f, double scale) {
    auto production = channels[i].first * f.t;
    production.AddScaled(channels[i].second, lts.process.CONT_TU_SIGN[i] * f.u);
    if (pair_index) { AddGoodWalker(lts, f.pair, production * scale, 0, regge.SoftModelHandle(), *pair_index); }
    return production.MultiplyScaled(decay, scale);
  });
}

}  // namespace

// Evaluate one complete MP, XP or GP central process
double MRegge::Amp2(gra::LORENTZSCALAR &lts, ReggeProductionModel model, MReggeMode mode) {
  lts.proton_good_walker.reset();
  const bool add_res = mode == MReggeMode::Resonance || mode == MReggeMode::ResonanceContinuumTwoBody;
  const bool add_con = mode == MReggeMode::ContinuumTwoBody || mode == MReggeMode::ContinuumTwoFourSixBody || mode == MReggeMode::ResonanceContinuumTwoBody;

  const std::size_t n = lts.decaytree.size();
  if (add_con && n != 2) {
    const std::complex<double> amplitude = n == 4 ? EvalContinuum4(lts, model) : EvalContinuum6(lts, model);
    return std::norm(amplitude);
  }

  const ForwardLegState upper_state = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState lower_state = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  gpom::AmpCache        gp_amp;
  ReggeSourceCache      source_cache(lts, *this, param, upper_state, lower_state);
  if (add_con) { EvalContinuum2(*this, lts, param, model, source_cache, gp_amp); }
  if (add_res) {
    std::size_t incoherent_group = 1;
    std::map<int, std::size_t> spin_groups;
    const bool blind_spin = !lts.process.SPINGEN || !lts.process.SPINDEC || lts.process.root_decay_mode == RootDecayMode::Isolated || n != 2;
    for (auto &entry : lts.process.RESONANCES) {
      const bool density = model == ReggeProductionModel::MP && entry.second.UsesDensitySpinBasis();
      std::size_t group = 0;
      if (density) {
        group = incoherent_group;
        ++incoherent_group;
      } else if (blind_spin) {
        // Angular averaging separates J sectors and removes continuum interference
        const auto [position, inserted] = spin_groups.try_emplace(entry.second.p.spinX2, incoherent_group);
        group = position->second;
        if (inserted) { ++incoherent_group; }
      }
      reggeamp::EvalResonance(*this, lts, param, entry.first, entry.second, model, gp_amp, source_cache, group);
    }
  }
  if (!lts.proton_good_walker) { throw AmplitudeFailure("MRegge::Amp2: no applicable Good Walker source"); }
  if (GoodWalkerOnly(lts, "MRegge::Amp2")) { return 0.0; }
  ProjectGoodWalkerBorn(lts);
  return BornNorm(lts);
}

}  // namespace gra
