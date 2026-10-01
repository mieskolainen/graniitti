// Regge photoproduction beam contraction for proton and nuclear collisions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggePhoto.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <map>
#include <optional>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Regge/MReggeCEP.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {
namespace {

// Build elementary gamma-proton subenergies for both directions
std::array<double, 2> PhotoSubenergies(const LORENTZSCALAR &lts, const std::array<ForwardLegState, 2> &elementary, const std::array<bool, 2> &direction) {
  std::array<double, 2> subenergy = {1.0, 1.0};
  if (direction[1]) { subenergy[0] = (lts.q2 + elementary[0].incoming).M2(); }
  if (direction[0]) { subenergy[1] = (lts.q1 + elementary[1].incoming).M2(); }
  for (const auto &i : indices(subenergy)) {
    if (direction[1 - i] && (!std::isfinite(subenergy[i]) || !(subenergy[i] > 0.0))) {
      throw AmplitudeFailure(
          "EvalReggePhoto: invalid gamma-proton "
          "subenergy");
    }
  }
  return subenergy;
}

// Derive the nuclear shadowing profile from one elementary Regge channel
flux::PhotoTargetProfile TargetProfile(const ReggePhotoAmplitude &term, const ReggePhotoWidths &width_ee) {
  const auto found = width_ee.find(term.central_pdg);
  if (!(term.w2 > 0.0)) { throw AmplitudeFailure("EvalReggePhoto: invalid generated gamma-proton subenergy"); }
  const std::complex<double> reduced   = term.target_forward / term.w2;
  const double               fv_over_e = std::sqrt(qed::alpha_QED() * term.central_mass / (3.0 * found->second));
  flux::PhotoTargetProfile   profile;
  profile.sigma_eff = fv_over_e * std::abs(reduced.imag()) * PDG::GeV2mb;
  profile.slope     = term.target_slope;
  profile.eta       = term.target_eta;
  profile.x         = term.central_mass * term.central_mass / term.w2;
  profile.scale2    = term.central_mass * term.central_mass / 4.0;
  if (!std::isfinite(profile.sigma_eff) || profile.sigma_eff < 0.0) { throw AmplitudeFailure("EvalReggePhoto: invalid vector-nucleon profile"); }
  return profile;
}

// Expand one elementary direction over explicit source and target sectors
//
nuclear::PhotoTerm ExpandDirection(const LORENTZSCALAR &lts, const ReggePhotoAmplitude &term, const ForwardLegState &upper, const ForwardLegState &lower, const ReggePhotoWidths &width_ee) {
  nuclear::PhotoTerm     output;
  const ForwardLegState &photon = term.photon_leg == ForwardBeamLeg::Upper ? upper : lower;
  const ForwardLegState &target = term.photon_leg == ForwardBeamLeg::Upper ? lower : upper;
  if (target.IsNuclear() && !(std::abs(term.target_elastic) > 0.0)) { return output; }
  const auto source = flux::PhotoSourceScalarAmplitudes(lts, photon);
  if (source.empty()) { return output; }

  std::vector<std::complex<double>> target_ratio = {1.0};
  if (target.IsNuclear()) {
    const flux::PhotoTargetProfile profile                                = TargetProfile(term, width_ee);
    auto                           factors                                = flux::ResolvePhotoTarget(target, profile, term.target_elastic, term.target_elastic);
    target_ratio                                                          = factors.factor;
    output.photo_current.at(static_cast<std::size_t>(target.Index() - 1)) = std::move(factors.current);
    for (auto &factor : target_ratio) { factor /= term.target_elastic; }
  }

  output.amplitude.reserve(term.amplitude.size() * source.size() * target_ratio.size());
  for (const auto &hard : term.amplitude) {
    for (const auto &emitter : source) {
      for (const auto &transition : target_ratio) { output.amplitude.push_back(hard * emitter * transition); }
    }
  }
  return output;
}

// Compute the elementary helicity-conserving vertex in the same kinematics as production
std::complex<double> PhotoVertex(const LORENTZSCALAR &lts, const regge::Param &param, const RES_PRODUCTION &production,
                                 ReggeProductionModel model, double t) {
  if (model == ReggeProductionModel::MP || model == ReggeProductionModel::XP) {
    return (model == ReggeProductionModel::MP ? mpom::PhotoCoupling : xpom::PhotoCoupling)(production);
  }
  return gpom::PhotoCoupling(lts, param, production, t);
}

// Elementary gamma-proton profile derived from one production model
struct ReggePhotoProfile {
  std::complex<double> forward = 0.0;
  double               slope   = 0.0;
  double               eta     = 0.0;
};

// Derive the vector-nucleon profile from a model-specific gamma-proton amplitude
ReggePhotoProfile ElementaryProfile(const MRegge &regge, const regge::Param &param, const LORENTZSCALAR &lts,
                                    const PARAM_RES &resonance, const RES_PRODUCTION &production, ReggeProductionModel model,
                                    double w2, int exchange_pdg, double dt) {
  const std::size_t trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  if (trajectory != param.pomeron_trajectory) { throw AmplitudeFailure("MRegge::PhotoTargetProfile: photo channel requires a Pomeron target"); }
  const SoftExchangeId exchange = param.exchanges[trajectory].soft_exchange;
  const auto kernel = [&](double t) { return regge.SoftModelHandle()->PhysicalResidue(exchange, t) * regge.PhotoKernel(w2, t, resonance.p.pdg); };
  const auto forward = PhotoVertex(lts, param, production, model, 0.0) * kernel(0.0);
  const auto nearby = PhotoVertex(lts, param, production, model, -dt) * kernel(-dt);
  ReggePhotoProfile profile;
  profile.forward = forward;
  // A vanishing forward amplitude gives zero VMD absorption and retains the finite-transfer channel
  const double norm0 = std::norm(forward);
  profile.slope = norm0 > 0.0 ? std::log(norm0 / std::norm(nearby)) / dt
                              : std::log(std::norm(kernel(0.0)) / std::norm(kernel(-dt))) / dt;
  const auto reduced = forward / w2;
  const double phase_scale = std::max(std::abs(reduced.real()), std::abs(reduced.imag()));
  if (std::abs(reduced.imag()) > std::numeric_limits<double>::epsilon() * phase_scale) { profile.eta = reduced.real() / reduced.imag(); }
  if (!std::isfinite(profile.slope) || !(profile.slope > 0.0)) { throw AmplitudeFailure("MRegge::PhotoTargetProfile: invalid vector-nucleon profile"); }
  return profile;
}

}  // namespace

// Pomeron kernel without vertex couplings for vector meson photoproduction
//
// K(V) = F_gammaPV(t) (W_yp^2 / W0^2)^(alpha(t) - 1), where
// F_gammaPV(t) = exp(B_gammaPV t / 2) and the Good Walker proton vertex is
// supplied separately by ReggeSourceCache
//
// [REFERENCE: H1 Collaboration, arxiv.org/abs/2005.14471]
// [REFERENCE: ZEUS Collaboration, arxiv.org/abs/hep-ex/0201043]
// [REFERENCE: ZEUS Collaboration, arxiv.org/abs/hep-ex/9601009]
// [REFERENCE: Ewerz, Maniatis, Nachtmann, arxiv.org/abs/1309.3478]

// Compute the anchored log amplitude and local trajectory correction for one double-log term
// [REFERENCE: arXiv:1307.7099, Eqs. (5), (17), fixed-scale evolution and local derivative dispersion approximation]
std::array<double, 2> PhotoEvolution(const regge::PhotoDLog &term, double s, double s0) {
  if (!(s > term.scale2)) { throw AmplitudeFailure("PhotoEvolution: subenergy is below the double-log scale"); }
  const double root = std::sqrt(std::log(s / term.scale2));
  const double delta = std::log(s / s0) / (root + term.root0);
  // The anchored log amplitude is quadratic in root-root0, avoiding cancellation near W0
  const double factor = -0.25 * term.c * delta / term.root0;
  return {factor * delta, factor / root};
}

// Compute the channel-specific photoproduction kernel without vertex couplings
std::complex<double> MRegge::PhotoKernel(const double s, const double t, const int central_pdg) const {
  if (!std::isfinite(s) || !(s > 0.0) || !std::isfinite(t)) { throw AmplitudeFailure("MRegge::PhotoKernel: invalid generated s or t"); }
  const auto &channel = regge::Photo(param, central_pdg);
  const double log_s = std::log(s / math::pow2(channel.W0));
  double log_amplitude = (channel.a0 - 1.0 + channel.ap * t) * log_s;
  double alpha0 = channel.a0;
  if (channel.dlog) {
    const auto evolution = PhotoEvolution(*channel.dlog, s, math::pow2(channel.W0));
    log_amplitude += evolution[0];
    alpha0 += evolution[1];
  }
  const auto &pomeron = soft_model->Exchange(param.exchanges.at(param.pomeron_trajectory).soft_exchange);
  const auto eta = regge::EtaFactor(alpha0 + channel.ap * t, alpha0, pomeron.Signature(), param.photoprod_eta_mode);
  return static_cast<double>(pomeron.ResidueSign()) * eta * s * std::exp(log_amplitude) *
         form::ExpSlopeAmplitude(channel.B_gammaPV, t);
}

// Compute the selected target transition from the elastic forward amplitude
// HERA uses its fitted W and t dependence with a density in dM_Y^2
// The smooth mass continuum extends to threshold without resolving individual N* resonances
// [REFERENCE: H1 Collaboration, arXiv:2005.14471, Eqs. (25), (41) and Tables 8, 11]
// [REFERENCE: H1 Collaboration, arXiv:1304.5162, Sections 2.2, 3.1, 3.2 and Tables 2, 3]
std::complex<double> MRegge::PhotoDissKernel(const double s, const double t, const double mass2, const int central_pdg, DissociationType type) const {
  if (type == DissociationType::Soft) {
    const auto exchange = param.exchanges[param.pomeron_trajectory].soft_exchange;
    return PhotoKernel(s, 0.0, central_pdg) * soft_model->ForwardExcitationFactor(exchange, t, mass2);
  }
  const auto& channel = regge::Photo(param, central_pdg);
  const auto& diss = *channel.diss;
  const double factor = flux::PhotoDissFactor(diss, s, t, mass2);
  const auto& pomeron = soft_model->Exchange(param.exchanges.at(param.pomeron_trajectory).soft_exchange);
  double alpha = 1.0 + diss.delta / 4.0;
  double log_amplitude = (channel.a0 - 1.0) * std::log(math::pow2(diss.W0 / channel.W0));
  if (channel.dlog) { log_amplitude += PhotoEvolution(*channel.dlog, math::pow2(diss.W0), math::pow2(channel.W0))[0]; }
  for (const auto &term : channel.diss_dlog) {
    const auto evolution = PhotoEvolution(term, s, math::pow2(diss.W0));
    log_amplitude += evolution[0];
    alpha += evolution[1];
  }
  const auto eta = regge::EtaFactor(alpha, alpha, pomeron.Signature(), param.photoprod_eta_mode);
  return static_cast<double>(pomeron.ResidueSign()) * eta * s * std::exp(log_amplitude) * factor;
}

// Build the immutable VMD data required by selected nuclear photo channels
ReggePhotoWidths BuildPhotoPlan(const LORENTZSCALAR &lts) {
  ReggePhotoWidths width_ee;
  if (lts.upc_model == nullptr) { return width_ee; }
  for (const auto &[name, resonance] : lts.process.RESONANCES) {
    static_cast<void>(name);
    const bool photo = std::any_of(resonance.production.cbegin(), resonance.production.cend(), [](const auto &production) {
      const auto &tree = production.tree;
      return tree.size() == 2 && ((tree[0].p.pdg == PDG::PDG_gamma) != (tree[1].p.pdg == PDG::PDG_gamma));
    });
    if (photo) { width_ee.emplace(resonance.p.pdg, resonance::ElectronicPartialWidth(resonance.p, resonance.modelparam.empty() ? gra::MODELPARAM : resonance.modelparam)); }
  }
  return width_ee;
}

// Resolve selected model-specific photo channels with caller supplied beam states
std::vector<ReggePhotoAmplitude> MRegge::PhotoAmplitudes(gra::LORENTZSCALAR &lts, const ReggeProductionModel model, const std::array<ForwardLegState, 2> &state, const std::array<double, 2> &subenergy,
                                                         const std::array<bool, 2> &direction) const {
  if (model != ReggeProductionModel::MP && model != ReggeProductionModel::XP && model != ReggeProductionModel::GP) { throw std::invalid_argument("MRegge::PhotoAmplitudes requires MP, XP or GP"); }
  lts.proton_good_walker.reset();
  gpom::AmpCache                   gp_amp;
  ReggeSourceCache                 source_cache(lts, *this, param, state[0], state[1], subenergy, {1.0, 1.0});
  std::vector<ReggePhotoAmplitude> output;
  for (auto &[name, resonance] : lts.process.RESONANCES) {
    MReggeDecayState        local_decay;
    const MReggeDecayState *decay = reggeamp::Decay(lts, name, resonance, reggeamp::RootFrame(lts, model), local_decay);
    for (const auto &channel : indices(resonance.production)) {
      const auto &production  = resonance.production[channel];
      const auto &tree        = production.tree;
      const bool  upper_gamma = tree.size() == 2 && tree[0].p.pdg == PDG::PDG_gamma;
      const bool  lower_gamma = tree.size() == 2 && tree[1].p.pdg == PDG::PDG_gamma;
      if (upper_gamma == lower_gamma || !direction[upper_gamma ? 0 : 1]) { continue; }

      ReggePhotoAmplitude term;
      term.photon_leg                     = upper_gamma ? ForwardBeamLeg::Upper : ForwardBeamLeg::Lower;
      const ForwardBeamLeg   target_leg   = upper_gamma ? ForwardBeamLeg::Lower : ForwardBeamLeg::Upper;
      const ForwardLegState &target       = state[target_leg == ForwardBeamLeg::Upper ? 0 : 1];
      const int              exchange_pdg = upper_gamma ? tree[1].p.pdg : tree[0].p.pdg;
      const std::size_t      trajectory   = regge::TrajectoryIndex(param, exchange_pdg);
      const SoftExchangeId   exchange     = param.exchanges[trajectory].soft_exchange;
      term.central_pdg                    = resonance.p.pdg;
      term.central_mass                   = resonance.p.mass;
      term.w2                             = source_cache.Subenergy(target_leg);
      const ReggePhotoProfile    profile  = ElementaryProfile(*this, param, lts, resonance, production, model, term.w2, exchange_pdg, regge_numerics->PhotoStep());
      const std::complex<double> shifted  = target.IsExcited()
          ? SoftModelHandle()->PhysicalResidue(exchange, 0.0) * PhotoDissKernel(term.w2, target.t, target.mass2, term.central_pdg, lts.process.PHOTO_DISSOCIATION)
          : SoftModelHandle()->PhysicalResidue(exchange, target.t) * PhotoKernel(term.w2, target.t, term.central_pdg);
      term.target_forward = profile.forward;
      term.target_elastic = shifted / (SoftModelHandle()->PhysicalResidue(exchange, 0.0) * PhotoKernel(term.w2, 0.0, term.central_pdg));
      term.target_slope   = profile.slope;
      term.target_eta     = profile.eta;

      lts.proton_good_walker.reset();
      reggeamp::EvalResonance(*this, lts, param, name, resonance, model, gp_amp, source_cache, 0, channel, decay, &state);
      if (!lts.proton_good_walker) { continue; }
      term.amplitude = ProjectGoodWalker(*lts.proton_good_walker);
      output.push_back(std::move(term));
    }
  }
  lts.proton_good_walker.reset();
  return output;
}

// Contract model-specific elementary directions with proton or nuclear beams
double EvalReggePhoto(MRegge &regge, LORENTZSCALAR &lts, ReggeProductionModel model, const ReggePhotoWidths &width_ee) {
  const auto upper = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const auto lower = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const std::array<bool, 2> active = {flux::SupportsPhotoTarget(lower), flux::SupportsPhotoTarget(upper)};
  const std::array<ForwardLegState, 2> elementary = {
      flux::ElementaryPhotoTarget(lts, upper), flux::ElementaryPhotoTarget(lts, lower)};
  const auto subenergy = PhotoSubenergies(lts, elementary, active);
  std::vector<nuclear::PhotoTerm> terms;
  for (const auto &term : regge.PhotoAmplitudes(lts, model, elementary, subenergy, active)) {
    auto expanded = ExpandDirection(lts, term, upper, lower, width_ee);
    if (expanded.amplitude.empty()) { continue; }
    expanded.direction = term.photon_leg == ForwardBeamLeg::Upper ? nuclear::PhotoDirection::Upper : nuclear::PhotoDirection::Lower;
    terms.push_back(std::move(expanded));
  }
  flux::SumPhotoTerms(lts, std::move(terms));
  return BornNorm(lts);
}

}  // namespace gra
