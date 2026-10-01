// Exclusive Z photoproduction amplitude using the GRANIITTI Sudakov-Shuvaev UGD
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MPhotoZ.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MException.h"

using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

namespace {

struct PhotoZFinalState {
  int  fermion_index     = 0;
  int  antifermion_index = 1;
  int  fermion_pdg       = 0;
  bool is_quark          = false;
  int  color_count       = 1;
};

struct PhotoZEWParameters {
  double q2         = 0.0;
  double alpha_qed  = 0.0;
  double e_qed      = 0.0;
  double sin2thetaW = 0.0;
  double z_mass     = 0.0;
  double z_width    = 0.0;
};

struct PhotoZProductionAmplitudes {
  std::complex<double> gamma = 0.0;
  std::complex<double> z     = 0.0;
};

struct PhotoZTargetChannels {
  std::array<PhotoZProductionAmplitudes, 2> amplitude{};
  std::optional<nuclear::PhotoCurrent>      current;
  std::size_t                               count = 1;
};

struct PhotoZProductionChannels {
  std::array<PhotoZProductionAmplitudes, 4> amplitude{};
  std::size_t                               count = 1;
};

struct PhotoZDecayCurrents {
  MDirac::Current gamma{};
  MDirac::Current z{};
};

// Compute true for charged leptons supported by the direct PhotoZ syntax
bool IsPhotoZLepton(int apdg) { return apdg == 11 || apdg == 13 || apdg == 15; }

// Compute true for light quarks supported by the direct PhotoZ syntax
bool IsPhotoZQuark(int apdg) { return apdg >= 1 && apdg <= 5; }

// Compute true for one direct same-flavour PhotoZ fermion pair
bool PhotoZProcessAccepts(const std::vector<MDecayBranch> &tree) {
  if (tree.size() != 2 || !tree[0].legs.empty() || !tree[1].legs.empty()) { return false; }
  const int first_pdg = tree[0].p.pdg;
  const int apdg      = std::abs(first_pdg);
  return first_pdg == -tree[1].p.pdg && (IsPhotoZLepton(apdg) || IsPhotoZQuark(apdg));
}

// Compute alpha_QED according to the existing GRANIITTI steering convention
double PhotoZAlphaQED(const gra::LORENTZSCALAR &lts, double q2) {
  constexpr double alpha_mg = 1.0 / 132.507;
  return gra::qed::alpha_QED(q2, lts.model_cache->Tune().Structure().QED_alpha, alpha_mg);
}

// Compute electroweak parameters from the particle table and model tune
PhotoZEWParameters BuildPhotoZEWParameters(const gra::LORENTZSCALAR &lts, double q2) {
  PhotoZEWParameters ew;
  ew.q2        = q2;
  ew.alpha_qed = PhotoZAlphaQED(lts, q2);
  ew.e_qed     = msqrt(4.0 * gra::math::PI * ew.alpha_qed);

  const gra::MParticle &z = lts.PDG.FindByPDG(23);
  const gra::MParticle &w = lts.PDG.FindByPDG(24);
  ew.z_mass               = z.mass;
  ew.z_width              = z.width;
  ew.sin2thetaW           = 1.0 - pow2(w.mass / z.mass);
  return ew;
}

// Compute the compact colour-dipole nucleon cross section in millibarns
// sigma_dip = pi^2 alpha_s xg/(3 Q2), converted from GeV^-2 to mb
// [REFERENCE: Coelho and Goncalves, Nucl. Phys. B956 (2020) 115013]
double PhotoZDipoleCrossSection(const gra::LORENTZSCALAR &lts, const double x, const double scale2) {
  try {
    const double alpha_s = lts.GlobalSudakovPtr->AlphaS_Q2(scale2);
    const double xg      = lts.GlobalSudakovPtr->xg_xQ2(x, scale2);
    return math::PIPI / (3.0 * scale2) * alpha_s * xg * PDG::GeV2mb;
  } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("MPhotoZ::PhotoZDipoleCrossSection: ") + error.what());
  }
}

// Compute the fermion ordering of the initialized neutral-current final state
PhotoZFinalState ResolvePhotoZFinalState(const gra::LORENTZSCALAR &lts) {
  const int pdg0 = lts.decaytree[0].p.pdg;
  const int apdg = std::abs(pdg0);
  const bool is_quark = IsPhotoZQuark(apdg);
  PhotoZFinalState state;
  state.fermion_index     = (pdg0 > 0) ? 0 : 1;
  state.antifermion_index = 1 - state.fermion_index;
  state.fermion_pdg       = std::abs(pdg0);
  state.is_quark          = is_quark;
  state.color_count       = is_quark ? 3 : 1;
  return state;
}

// Compute the timelike propagator for a virtual photon or Z boson
// D_gamma = 1/q2, D_Z = 1/(q2-M_Z^2+i M_Z Gamma_Z)
std::complex<double> NeutralBosonPropagator(bool z_boson, const PhotoZEWParameters &ew) {
  if (!z_boson) { return 1.0 / ew.q2; }
  return 1.0 / std::complex<double>(ew.q2 - pow2(ew.z_mass), ew.z_mass * ew.z_width);
}

// Compute the gamma* and Z* impact amplitudes with the CSS colour normalization
// [REFERENCE: Cisek, Schafer, Szczurek, arXiv:0906.1739, Eqs. (2.1)-(2.3)]
PhotoZProductionAmplitudes PhotoZCSSImpactAmplitudes(const gra::LORENTZSCALAR &lts, const gra::MPhotoQCDParam &param,
                                                     const gra::MPhotoQCDNumerics &num, const PhotoZEWParameters &ew,
                                                     double w2) {
  if (!(w2 > ew.q2)) { return {}; }

  const double xg = gra::MPhotoQCD::XGluon(param, ew.q2, w2);
  if (!(xg > 0.0 && xg < 1.0)) { return {}; }

  const double sw           = std::sqrt(ew.sin2thetaW);
  const double cw           = std::sqrt(1.0 - ew.sin2thetaW);
  const double mu           = gra::MPhotoQCD::HardScale(param, ew.q2);
  const double prefactor    = 2.0 * ew.alpha_qed / gra::math::PI;

  PhotoZProductionAmplitudes out;
  for (const int qpdg : {1, 2, 3, 4, 5}) {
    const gra::qed::NeutralCurrentCouplings c = gra::qed::FermionNeutralCurrentCouplings(qpdg, ew.sin2thetaW);
    const gra::MParticle                   &q = lts.PDG.FindByPDG(qpdg);
    const std::complex<double> flavour = gra::MPhotoQCD::CSSFlavorIntegral(lts, param, num, xg, q.mass, ew.q2, mu);
    const double               gamma_coupling = pow2(c.chargeX1);
    const double               z_coupling     = c.chargeX1 * c.gV / (sw * cw);
    out.gamma += prefactor * gamma_coupling * flavour;
    out.z += prefactor * z_coupling * flavour;
  }
  return out;
}

// Compute one photon-Pomeron production direction before spin phase
PhotoZTargetChannels DirectionTargetAmplitudes(const gra::LORENTZSCALAR &lts, const gra::MPhotoQCDParam &param,
                                               const gra::MPhotoQCDNumerics &num, const PhotoZEWParameters &ew,
                                               bool photon_from_upper, const SoftModel &soft_model) {
  const gra::M4Vec     &q_photon = photon_from_upper ? lts.q1 : lts.q2;
  const ForwardLegState target_state =
      ResolveForwardLegState(lts, photon_from_upper ? ForwardBeamLeg::Lower : ForwardBeamLeg::Upper);
  PhotoZTargetChannels out;

  if (!gra::flux::SupportsPhotoTarget(target_state)) { return out; }

  const gra::M4Vec p_target     = gra::flux::PhotoTargetMomentum(target_state);
  const double     w2           = (q_photon + p_target).M2();
  const double     target_mass2 = target_state.IsExcited() ? target_state.mass2 : p_target.M2();
  if (!(target_mass2 > 0.0) || !(w2 > pow2(msqrt(ew.q2) + msqrt(target_mass2)))) {
    auto target = gra::flux::ZeroPhotoTarget(target_state);
    out.count   = target.factor.size();
    out.current = std::move(target.current);
    return out;
  }

  const PhotoZProductionAmplitudes impact  = PhotoZCSSImpactAmplitudes(lts, param, num, ew, w2);
  const double                     elastic = gra::MPhotoQCD::TSlopeFactor(param, w2, target_state.t);
  const double                     hadron = gra::MPhotoQCD::TargetTransitionFactor(target_state, param, w2, soft_model);
  gra::flux::PhotoTargetProfile    profile;
  profile.slope                     = gra::MPhotoQCD::TSlope(param, w2);
  profile.eta                       = gra::MPhotoQCD::RealPartFactor(param).real();
  profile.x                         = param.xg_coefficient * ew.q2 / w2;
  profile.scale2                    = std::max(param.mu_MIN * param.mu_MIN, param.mu_over_m * param.mu_over_m * ew.q2);
  if (target_state.IsNuclear()) { profile.sigma_eff = PhotoZDipoleCrossSection(lts, profile.x, profile.scale2); }
  auto target                       = gra::flux::ResolvePhotoTarget(target_state, profile, elastic, hadron);
  out.current                       = std::move(target.current);
  const std::complex<double> common = gra::MPhotoQCD::RealPartFactor(param) * w2;
  out.count                         = std::min(out.amplitude.size(), target.factor.size());
  for (std::size_t i = 0; i < out.count; ++i) {
    out.amplitude[i].gamma = common * impact.gamma * target.factor[i];
    out.amplitude[i].z     = common * impact.z * target.factor[i];
  }
  return out;
}

// Couple explicit emitter sectors to all target-sector amplitudes
PhotoZProductionChannels DirectionProductionAmplitudes(const gra::LORENTZSCALAR   &lts,
                                                       const PhotoZTargetChannels &target, const int photon_leg,
                                                       const int lambda) {
  const ForwardLegState state =
      ResolveForwardLegState(lts, photon_leg == 1 ? ForwardBeamLeg::Upper : ForwardBeamLeg::Lower);
  const auto               source = gra::flux::PhotoSourceAmplitudes(lts, state, lambda);
  PhotoZProductionChannels out;
  out.count           = source.size() * target.count;
  std::size_t channel = 0;
  for (const std::complex<double> photon : source) {
    for (std::size_t transition = 0; transition < target.count; ++transition) {
      out.amplitude[channel++] = PhotoZProductionAmplitudes{photon * target.amplitude[transition].gamma,
                                                            photon * target.amplitude[transition].z};
    }
  }
  return out;
}

// Compute decay currents that are independent of the parent boson helicity
std::vector<PhotoZDecayCurrents> BuildPhotoZDecayCurrentCache(const gra::MDecayBranch  &fermion,
                                                              const gra::MDecayBranch  &anti,
                                                              const PhotoZFinalState   &state,
                                                              const PhotoZEWParameters &ew, const gra::MDirac &dirac) {
  std::vector<PhotoZDecayCurrents> currents;
  currents.reserve(4);
  constexpr auto helicities = spin::BinaryHelicityLabelsX2();
  for (const int hf : helicities) {
    for (const int ha : helicities) {
      PhotoZDecayCurrents cache;
      cache.gamma = gra::qed::PhotonFermionCurrent(dirac, fermion.p4, anti.p4, state.fermion_pdg, hf, ha, ew.e_qed);
      cache.z =
          gra::qed::ZFermionCurrent(dirac, fermion.p4, anti.p4, state.fermion_pdg, hf, ha, ew.e_qed, ew.sin2thetaW);
      currents.push_back(cache);
    }
  }
  return currents;
}

// Compute all PhotoZ helicity amplitudes for one validated final state
std::vector<std::complex<double>> BuildPhotoZHelicityAmplitudes(gra::LORENTZSCALAR &lts, const PhotoZFinalState &state,
                                                                const gra::MPhotoQCDParam    &param,
                                                                const gra::MPhotoQCDNumerics &num,
                                                                const PhotoZEWParameters &ew, const gra::MDirac &dirac,
                                                                const SoftModel &soft_model) {
  const gra::MDecayBranch   &fermion             = lts.decaytree[state.fermion_index];
  const gra::MDecayBranch   &anti                = lts.decaytree[state.antifermion_index];
  const gra::M4Vec           boson               = fermion.p4 + anti.p4;
  const ForwardLegState      upper_state         = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState      lower_state         = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const bool                 separate_directions = gra::flux::SplitPhotoDirections(lts, upper_state, lower_state);
  const std::complex<double> gamma_prop          = NeutralBosonPropagator(false, ew);
  const std::complex<double> z_prop              = NeutralBosonPropagator(true, ew);
  const std::vector<PhotoZDecayCurrents> decay_currents = BuildPhotoZDecayCurrentCache(fermion, anti, state, ew, dirac);
  PhotoZTargetChannels upper_target           = DirectionTargetAmplitudes(lts, param, num, ew, true, soft_model);
  PhotoZTargetChannels lower_target           = DirectionTargetAmplitudes(lts, param, num, ew, false, soft_model);
  lts.screening.photo.current[1]              = std::move(upper_target.current);
  lts.screening.photo.current[0]              = std::move(lower_target.current);
  const std::size_t upper_count     = gra::flux::PhotoSourceCount(upper_state) * upper_target.count;
  const std::size_t lower_count     = gra::flux::PhotoSourceCount(lower_state) * lower_target.count;
  const std::size_t channel_count   = separate_directions ? upper_count + lower_count : 1;

  std::vector<std::complex<double>> amplitudes;
  amplitudes.assign(static_cast<std::size_t>(4 * state.color_count) * channel_count, 0.0);
  const auto polarization = MPhotoQCD::TransverseStates(lts, boson);

  // The CSS high-energy impact factor is used in its transverse SCHC
  // approximation
  constexpr auto helicities = spin::BinaryHelicityLabelsX2();
  for (const int lambda : helicities) {
    const FTensor::Tensor1<std::complex<double>, 4> eps   = polarization[spin::BinaryHelicityIndexX2(lambda)];
    const PhotoZProductionChannels                  upper = DirectionProductionAmplitudes(lts, upper_target, 1, lambda);
    const PhotoZProductionChannels                  lower = DirectionProductionAmplitudes(lts, lower_target, 2, lambda);

    std::size_t row = 0;
    for (const PhotoZDecayCurrents &current : decay_currents) {
      const std::complex<double> decay_gamma = gra::MPhotoQCD::ContractCurrentPolarization(current.gamma, eps);
      const std::complex<double> decay_z     = gra::MPhotoQCD::ContractCurrentPolarization(current.z, eps);
      for (int color = 0; color < state.color_count; ++color) {
        if (!separate_directions) {
          amplitudes[row++] += (upper.amplitude[0].gamma + lower.amplitude[0].gamma) * gamma_prop * decay_gamma +
                               (upper.amplitude[0].z + lower.amplitude[0].z) * z_prop * decay_z;
        } else {
          // Store emitter sectors outside target sectors in each direction
          for (std::size_t i = 0; i < upper.count; ++i) {
            amplitudes[row++] += upper.amplitude[i].gamma * gamma_prop * decay_gamma +
                                 upper.amplitude[i].z * z_prop * decay_z;
          }
          for (std::size_t i = 0; i < lower.count; ++i) {
            amplitudes[row++] += lower.amplitude[i].gamma * gamma_prop * decay_gamma +
                                 lower.amplitude[i].z * z_prop * decay_z;
          }
        }
      }
    }
  }
  return gra::flux::CompletePhotoInitialSpinStates(lts, amplitudes);
}

}  // namespace

// Build the immutable direct neutral-current process definition
std::shared_ptr<const amplitude::ProcessDefinition> MPhotoZ::ProcessDefinitionFor() {
  return std::make_shared<amplitude::AnalyticProcess>("MPHOTOZ", "ygg_Z", "charged lepton or quark pair (u/d/s/c/b)",
                                                      DirectDecayStructure(), PhotoZProcessAccepts,
                                                      [](const LORENTZSCALAR &) { return DirectDecayStructure(); });
}

// Initialize photoproduction and Sudakov parameters before worker copies
void MPhotoZ::InitializeParameters(MProcessSetup &setup) {
  if (!PhotoZProcessAccepts(setup.lts.decaytree)) {
    throw std::invalid_argument("MPhotoZ::InitializeParameters: unsupported direct fermion pair");
  }
  if (!setup.model_tune) { throw std::invalid_argument("MPhotoZ::InitializeParameters: missing model tune"); }
  MModelCache &cache = RequireModelCache(setup.lts.model_cache, setup.model_tune, "MPhotoZ initialization");
  (void)GetPhotoQCDParam(cache, "PARAM_PHOTOZ");
  (void)GetPhotoQCDNumerics(cache, "NUMERICS_PHOTOZ");
  MPhotoQCD::EnsureSudakov(setup.lts, setup.soft_model, "MPhotoZ::InitializeParameters");
}

// Construct the photoproduction amplitude with one immutable process definition
MPhotoZ::MPhotoZ(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune_snapshot,
                 std::shared_ptr<const amplitude::ProcessDefinition> definition)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(RequireModelCache(lts.model_cache, model_tune_snapshot, "MPhotoZ").TunePtr()),
      soft_model(model_tune->Soft()),
      dirac("DIRAC") {
  if (!PhotoZProcessAccepts(lts.decaytree)) {
    throw std::invalid_argument("MPhotoZ: unsupported direct fermion pair");
  }
  MModelCache &cache = *lts.model_cache;
  param              = GetPhotoQCDParam(cache, "PARAM_PHOTOZ");
  numerics           = GetPhotoQCDNumerics(cache, "NUMERICS_PHOTOZ");
  const auto ew      = BuildPhotoZEWParameters(lts, pow2(lts.PDG.FindByPDG(23).mass));
  if (!std::isfinite(ew.sin2thetaW) || ew.sin2thetaW <= 0.0 || ew.sin2thetaW >= 1.0) {
    throw std::invalid_argument("MPhotoZ::BuildPhotoZEWParameters: invalid weak mixing angle");
  }
  EnsureSudakov(lts);
}

// Acquire shared read-only Sudakov and UGD tables
void MPhotoZ::EnsureSudakov(gra::LORENTZSCALAR &lts) const {
  gra::MPhotoQCD::EnsureSudakov(lts, soft_model, "MPhotoZ::EnsureSudakov");
}

// Evaluate gamma p to neutral-current f fbar kinematics
double MPhotoZ::Amp2(gra::LORENTZSCALAR &lts) const {
  // Preserve both exact photon-emitter and target forward systems through LHE
  // conversion
  lts.exact_forward_photon_kinematics = true;
  const PhotoZFinalState state        = ResolvePhotoZFinalState(lts);
  const double           q2           = gra::MPhotoQCD::CentralMass2(lts);
  if (!(q2 > 0.0) || !std::isfinite(q2)) {
    throw AmplitudeFailure("MPhotoZ::Amp2: central invariant mass must be positive");
  }

  const PhotoZEWParameters ew = BuildPhotoZEWParameters(lts, q2);
  lts.hard_color_flows.clear();
  lts.hamp = BuildPhotoZHelicityAmplitudes(lts, state, *param, *numerics, ew, dirac, *soft_model);
  gra::flux::ConfigurePhotoLayout(lts);
  return gra::flux::PhotoInitialSpinAverage(lts) * gra::SquaredNorm(lts.hamp);
}

// Assign shower-compatible color flow for quark-pair final states
void MPhotoZ::SampleColorFlow(gra::LORENTZSCALAR &lts) const {
  const PhotoZFinalState state = ResolvePhotoZFinalState(lts);
  for (auto &branch : lts.decaytree) { branch.p.color_flow.clear(); }
  if (!state.is_quark) { return; }

  lts.decaytree[state.fermion_index].p.color_flow.flow1     = 501;
  lts.decaytree[state.antifermion_index].p.color_flow.flow2 = 501;
}

// Compute the on-shell photon-proton to Z-proton total cross section in nb
// sigma = |A(t=0)|^2/[16 pi W^4 B(W)]
double MPhotoZ::GammaPTotalCrossSectionNb(gra::LORENTZSCALAR &lts, double W) const {
  EnsureSudakov(lts);
  if (!(W > 0.0) || !std::isfinite(W)) {
    throw std::invalid_argument("MPhotoZ::GammaPTotalCrossSectionNb: W must be positive");
  }

  const double w2 = pow2(W);
  const double q2 = pow2(lts.PDG.FindByPDG(23).mass);
  if (!(w2 > pow2(msqrt(q2) + gra::PDG::mp))) { return 0.0; }

  const PhotoZEWParameters         ew        = BuildPhotoZEWParameters(lts, q2);
  const PhotoZProductionAmplitudes impact    = PhotoZCSSImpactAmplitudes(lts, *param, *numerics, ew, w2);
  const std::complex<double>       amplitude = MPhotoQCD::RealPartFactor(*param) * w2 * impact.z;
  const double                     slope     = MPhotoQCD::TSlope(*param, w2);
  if (!(slope > 0.0)) { return 0.0; }

  const double dsigma_dt0 = std::norm(amplitude) / (16.0 * math::PI * pow2(w2));
  return dsigma_dt0 / slope * PDG::GeV2barn * 1.0e9;
}

}  // namespace gra
