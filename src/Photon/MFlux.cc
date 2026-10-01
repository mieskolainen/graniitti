// Photon and other fluxes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <stdexcept>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MPhoton.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MPhotoQCD.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace flux {

using gra::aux::indices;
using math::abs2;
using math::msqrt;

namespace {

// Check whether a flux momentum fraction is physical
bool ValidMomentumFraction(double x) { return x > 0.0 && x < 1.0; }

// Check whether a cross-section density is finite and positive
bool PositiveFinite(double value) { return std::isfinite(value) && value > 0.0; }

// Compute the explicit photon-emission sectors of one forward leg
std::vector<nuclear::CoherenceType> EmissionSectors(const ForwardLegState &state) {
  const auto sector = state.upc_model != nullptr ? state.upc_model->SourceSector(state.Index(), state.emission)
                                                : state.emission;
  return state.IsNuclear() ? nuclear::CoherenceSectors(sector)
                           : std::vector<nuclear::CoherenceType>{state.emission};
}

// Factor one positive diagonal photon density into real eigen-sources
// rho_T = diag(a_parallel^2,a_perpendicular^2)
TransversePhotonSource DensityEigenSource(const TransversePhotonFlux &density) {
  return {std::sqrt(density.parallel), std::sqrt(density.perpendicular)};
}

// Compute one explicit transverse photon source without inferring its phase
TransversePhotonSource PhotonSectorSource(const ForwardLegState &state, const nuclear::CoherenceType emission,
                                          const TransversePhotonFlux &density) {
  if (!state.IsNuclear() || emission != nuclear::CoherenceType::Coherent) { return DensityEigenSource(density); }
  if (state.upc_model == nullptr) { return {}; }
  const nuclear::MPhoton *photon = state.upc_model->Photon(state.Index());
  if (photon == nullptr) { return {}; }
  return {photon->CoherentAmp(state.xi, state.t, state.qt), 0.0};
}

// Compute the scalar source used by a helicity-averaged hard amplitude
double ScalarSectorSource(const ForwardLegState &state, const nuclear::CoherenceType emission,
                          const PhotonFluxSector &sector) {
  if (state.IsNuclear() && emission == nuclear::CoherenceType::Coherent) { return sector.source.parallel; }
  return std::hypot(sector.source.parallel, sector.source.perpendicular);
}

// Compute the explicit photonuclear target sectors of one forward leg
std::vector<nuclear::CoherenceType> TargetSectors(const ForwardLegState &state) {
  const auto sector = state.upc_model != nullptr ? state.upc_model->SourceSector(state.Index(), state.target)
                                                : state.target;
  return state.IsNuclear() ? nuclear::CoherenceSectors(sector)
                           : std::vector<nuclear::CoherenceType>{state.target};
}

// Evaluate the Drees-Zeppenfeld form-factor integral without cancellation
// I = ln(A)-11/6+3/A-3/(2A^2)+1/(3A^3), A = 1+0.71/Qmin^2
double DZIntegralBracket(double Q2min) {
  if (!(Q2min > 0.0) || !std::isfinite(Q2min)) { return 0.0; }

  constexpr double dipole_scale2 = 0.71;
  const double     w             = dipole_scale2 / (Q2min + dipole_scale2);
  if (w < 1.0e-3) {
    double sum   = 0.0;
    double power = math::pow4(w);
    for (int order = 4; order <= 12; ++order) {
      sum += power / static_cast<double>(order);
      power *= w;
    }
    return sum;
  }
  return std::log1p(dipole_scale2 / Q2min) - w - math::pow2(w) / 2.0 - math::pow3(w) / 3.0;
}

// Set all stored helicity amplitudes to zero
void ZeroHelicityAmplitudes(gra::LORENTZSCALAR &lts) {
  for (const auto &i : indices(lts.hamp)) { lts.hamp[i] = 0.0; }
}

// Store the ordered nuclear emission sectors for screening
void StoreEmissionSectorMetadata(gra::LORENTZSCALAR                                       &lts,
                                 const std::array<std::vector<nuclear::CoherenceType>, 2> &sector,
                                 const std::size_t                                         rows) {
  lts.hamp.layout.nuclear_type        = nuclear::ScreenType::Fusion;
  lts.hamp.layout.epa_sector_resolved = true;
  lts.hamp.layout.epa_rows_per_sector = rows;
  for (const auto &leg : indices(sector)) {
    lts.hamp.layout.epa_sector_count[leg] = static_cast<std::uint8_t>(sector[leg].size());
    lts.hamp.layout.epa_sector_type[leg]  = {255, 255};
    for (const auto &i : indices(sector[leg])) {
      lts.hamp.layout.epa_sector_type[leg][i] = static_cast<std::uint8_t>(sector[leg][i]);
    }
  }
}

// Compute amplitude weights for the resolved sectors of an inclusive nuclear EPA state
// w_ij = a_i a_j / sqrt[(sum_k n_k)(sum_l n_l)]
EPASectorWeights ResolveInclusiveEmissionSectors(gra::LORENTZSCALAR &lts, const ForwardLegState &upper,
                                                 const ForwardLegState &lower) {
  EPASectorWeights weights;
  if (lts.hamp.layout.epa_sector_resolved) { return weights; }
  const auto upper_sector = EmissionSectors(upper);
  const auto lower_sector = EmissionSectors(lower);
  StoreEmissionSectorMetadata(lts, {upper_sector, lower_sector}, 1);

  const auto upper_flux  = ForwardPhotonFluxSectors(upper);
  const auto lower_flux  = ForwardPhotonFluxSectors(lower);
  double     upper_total = 0.0;
  double     lower_total = 0.0;
  for (const auto &sector : upper_flux) { upper_total += sector.density.Trace(); }
  for (const auto &sector : lower_flux) { lower_total += sector.density.Trace(); }
  if (!(upper_total > 0.0) || !(lower_total > 0.0)) {
    weights.amplitude_fraction = {0.0};
    return weights;
  }

  const double norm = std::sqrt(upper_total * lower_total);
  weights.amplitude_fraction.clear();
  weights.amplitude_fraction.reserve(upper_sector.size() * lower_sector.size());
  for (const auto &i : indices(upper_sector)) {
    const double upper_source = ScalarSectorSource(upper, upper_sector[i], upper_flux[i]);
    for (const auto &j : indices(lower_sector)) {
      const double lower_source = ScalarSectorSource(lower, lower_sector[j], lower_flux[j]);
      weights.amplitude_fraction.push_back(upper_source * lower_source / norm);
    }
  }
  return weights;
}

}  // namespace

// Compute ordered forward photon densities without exposing their source model
std::vector<PhotonFluxSector> ForwardPhotonFluxSectors(const ForwardLegState &state) {
  const auto                    sector = EmissionSectors(state);
  std::vector<PhotonFluxSector> output(sector.size());
  for (const auto &i : indices(sector)) {
    output[i].density =
        state.IsNuclear() ? NuclearPhotonFluxTransverse(state, sector[i]) : ForwardPhotonFluxTransverse(state);
    output[i].source = PhotonSectorSource(state, sector[i], output[i].density);
    output[i].code   = static_cast<std::uint8_t>(sector[i]);
  }
  return output;
}

// Validate one configured photon emitter through the common source layer
void ValidatePhotonEmitter(const LORENTZSCALAR &lts, const int leg, const std::string &context) {
  if (leg != 1 && leg != 2) { throw std::invalid_argument(context + ": photon leg should be 1 or 2"); }
  const MParticle &beam = leg == 1 ? lts.beam1 : lts.beam2;
  if (nuclear::IsNuclearPDG(beam.pdg)) {
    if (lts.upc_model == nullptr || lts.upc_model->Photon(leg) == nullptr) {
      throw std::invalid_argument(context + ": missing nuclear photon source");
    }
    return;
  }
  qed::ValidateEmitter(lts, leg, context);
}

// Compute whether one forward state supports the elementary photon target
bool SupportsPhotoTarget(const ForwardLegState &state) {
  if (state.IsNuclear()) {
    return state.upc_model != nullptr && state.upc_model->Nucleus(state.Index()) != nullptr &&
           state.upc_model->Photo(state.Index()) != nullptr;
  }
  return std::abs(state.emitter.pdg) == PDG::PDG_p;
}

// Compute whether one ordered photon beam direction has a physical target
bool SupportsPhotoDirection(const LORENTZSCALAR &lts, const int photon_leg) {
  if (photon_leg != 1 && photon_leg != 2) {
    throw std::invalid_argument("SupportsPhotoDirection: photon leg should be one or two");
  }
  const int        target_leg = photon_leg == 1 ? 2 : 1;
  const MParticle &target     = target_leg == 1 ? lts.beam1 : lts.beam2;
  if (nuclear::IsNuclearPDG(target.pdg)) {
    return lts.upc_model != nullptr && lts.upc_model->Photo(target_leg) != nullptr;
  }
  return nuclear::IsProton(target.pdg);
}

// Compute the hadron or constituent-nucleon momentum seen by the hard amplitude
M4Vec PhotoTargetMomentum(const ForwardLegState &state) {
  if (!state.IsNuclear()) { return state.incoming; }
  if (state.upc_model == nullptr) { return {}; }
  const nuclear::MNucleus *nucleus = state.upc_model->Nucleus(state.Index());
  if (nucleus == nullptr || nucleus->A() == 0) { return {}; }
  return state.incoming / static_cast<double>(nucleus->A());
}

// Construct zero target amplitudes with the configured sectors and current bank
PhotoTargetFactors ZeroPhotoTarget(const ForwardLegState &state) {
  PhotoTargetFactors output;
  output.factor.assign(TargetSectors(state).size(), 0.0);
  if (state.IsNuclear() && state.upc_model != nullptr) {
    if (const auto *bank = state.upc_model->Bank(state.Index()); bank != nullptr) {
      output.current.emplace();
      output.current->sample.assign(bank->Size(), 0.0);
    }
  }
  return output;
}

// Resolve hadron and nuclear target factors outside the hard amplitude module
PhotoTargetFactors ResolvePhotoTarget(const ForwardLegState &state, const PhotoTargetProfile &profile,
                                      const std::complex<double> elastic_factor,
                                      const std::complex<double> hadron_factor) {
  PhotoTargetFactors output;
  if (!state.IsNuclear()) {
    output.factor = {hadron_factor};
    return output;
  }

  const auto sector = TargetSectors(state);
  output.factor.assign(sector.size(), 0.0);
  if (state.upc_model == nullptr) {
    return output;
  }
  const nuclear::MPhoto   *photo   = state.upc_model->Photo(state.Index());
  const nuclear::MNucleus *nucleus = state.upc_model->Nucleus(state.Index());
  if (photo == nullptr || nucleus == nullptr || nucleus->A() == 0) { return output; }
  // Convert the reduced nucleon current to the external nuclear invariant
  // amplitude: M_A = [2 q.P_A / (2 q.p_N)] J_A M_N = A J_A M_N
  const double          target_norm = static_cast<double>(nucleus->A());
  nuclear::PhotoProfile nuclear_profile;
  nuclear_profile.sigma_eff = profile.sigma_eff;
  nuclear_profile.slope     = profile.slope;
  nuclear_profile.eta       = profile.eta;
  nuclear_profile.x         = profile.x;
  nuclear_profile.scale2    = profile.scale2;
  nuclear_profile.isospin   = profile.isospin;
  nuclear::PhotoTransition transition;
  try {
    const M3Vec                    transfer = nuclear::RestTransfer(state.incoming, state.transfer, state.emitter.mass);
    const nuclear::PhotonDirection direction = nuclear::TargetPhotonDirection(state.Index());
    transition = photo->Factors(nuclear_profile, transfer[0], transfer[1], transfer[2], direction,
                                state.upc_model->Bank(state.Index()), state.upc_model->HotSpot(state.Index()),
                                state.upc_model->HasSamples(), state.target);
  } catch (const AmplitudeFailure &) { throw; } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("ResolvePhotoTarget: ") + error.what());
  }
  output.current = std::move(transition.current);
  for (const auto &i : indices(sector)) {
    if (sector[i] == nuclear::CoherenceType::Coherent) {
      output.factor[i] = target_norm * elastic_factor * transition.coherent;
    } else if (sector[i] == nuclear::CoherenceType::Incoherent) {
      output.factor[i] = target_norm * elastic_factor * transition.incoherent;
    }
  }
  return output;
}

// Complete the elementary photoproduction amplitudes over initial spin states
std::vector<std::complex<double>> CompletePhotoInitialSpinStates(const LORENTZSCALAR                     &lts,
                                                                 const std::vector<std::complex<double>> &amplitudes) {
  return lts.upc_model != nullptr ? amplitudes : MPhotoQCD::InitialProtonSpinCopies(amplitudes);
}

// Compute the initial-state spin average for one photoproduction process
double PhotoInitialSpinAverage(const LORENTZSCALAR &lts) { return lts.upc_model == nullptr ? 0.25 : 1.0; }

// Compute the number of resolved photon-source sectors for one forward leg
std::size_t PhotoSourceCount(const ForwardLegState &state) { return EmissionSectors(state).size(); }

// Compute whether two photon directions need distinct amplitude channels
bool SplitPhotoDirections(const LORENTZSCALAR &lts, const ForwardLegState &upper, const ForwardLegState &lower) {
  const bool sampled = (upper.upc_model != nullptr && upper.upc_model->HasSamples()) ||
                       (lower.upc_model != nullptr && lower.upc_model->HasSamples());
  const bool config_survival = lts.upc_model != nullptr && lts.upc_model->HadronicConvolution() &&
                               lts.upc_model->Param().survival == nuclear::SurvivalType::MCGGCF;
  const bool emd = lts.upc_model && lts.upc_model->Param().reaction && lts.upc_model->Param().additional_emd;
  return sampled || emd || config_survival || upper.IsExcited() || lower.IsExcited() ||
         (upper.IsNuclear() &&
          (upper.emission != nuclear::CoherenceType::Coherent || upper.target != nuclear::CoherenceType::Coherent)) ||
         (lower.IsNuclear() &&
          (lower.emission != nuclear::CoherenceType::Coherent || lower.target != nuclear::CoherenceType::Coherent));
}

// Store the canonical direction and nuclear-sector order for screening
void ConfigurePhotoLayout(LORENTZSCALAR &lts) {
  auto &layout = lts.hamp.layout;
  layout.nuclear_type          = lts.upc_model == nullptr ? nuclear::ScreenType::Scalar : nuclear::ScreenType::Photo;
  layout.photo_sector_resolved = false;
  layout.photo_channel_count   = 0;
  layout.photo_channel         = {};
  if (lts.upc_model == nullptr) { return; }
  const ForwardLegState upper = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState lower = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  if (SplitPhotoDirections(lts, upper, lower)) { layout.ConfigurePhotoChannels(*upper.upc_model); }
}

// Compute one explicit EPA source-sector amplitude with its helicity phase
std::complex<double> PhotoSourceAmplitude(const LORENTZSCALAR &lts, const ForwardLegState &state,
                                          const nuclear::CoherenceType emission, const int helicity) {
  const double recoil = std::sqrt(std::max(0.0, 1.0 - state.xi));
  if (!state.IsNuclear()) { return recoil * qed::TransversePhotonSourceAmplitude(lts, state, helicity); }
  const TransversePhotonFlux density = NuclearPhotonFluxTransverse(state, emission);
  PhotonFluxSector           sector;
  sector.density      = density;
  sector.source       = PhotonSectorSource(state, emission, density);
  const double source = ScalarSectorSource(state, emission, sector);
  return recoil * nuclear::PhotonAmp(source, state.xi, qed::Charge(lts, state.Index()), state.transfer.Phi(), helicity);
}

// Compute all resolved photon-source amplitudes for one helicity
std::vector<std::complex<double>> PhotoSourceAmplitudes(const LORENTZSCALAR &lts, const ForwardLegState &state,
                                                        const int helicity) {
  const auto                        sector = EmissionSectors(state);
  std::vector<std::complex<double>> amplitudes;
  amplitudes.reserve(sector.size());
  for (const auto emission : sector) { amplitudes.push_back(PhotoSourceAmplitude(lts, state, emission, helicity)); }
  return amplitudes;
}

// Compute resolved scalar photon-source factors without helicity phases
std::vector<std::complex<double>> PhotoSourceScalarAmplitudes(const LORENTZSCALAR &lts, const ForwardLegState &state) {
  // d qT^2 = (1-xi) d Q^2 converts the EPA density to the exact one-photon recoil phase space
  // [REFERENCE: Frixione et al., arXiv:hep-ph/9310350, Eqs. (17)-(20)]
  const double recoil = std::sqrt(std::max(0.0, 1.0 - state.xi));
  if (!state.IsNuclear()) { return {recoil * qed::PhotonSourceAmplitude(lts, state, lts.process.PHOTON_VERTEX)}; }
  const auto                        sector = EmissionSectors(state);
  std::vector<std::complex<double>> amplitudes;
  amplitudes.reserve(sector.size());
  const double charge = qed::Charge(lts, state.Index());
  if (!std::isfinite(charge) || !(std::abs(charge) > 0.0) || !(state.xi > 0.0)) {
    return std::vector<std::complex<double>>(sector.size(), 0.0);
  }
  const double charge_phase = charge / std::abs(charge);
  for (const auto emission : sector) {
    PhotonFluxSector resolved;
    resolved.density    = NuclearPhotonFluxTransverse(state, emission);
    resolved.source     = PhotonSectorSource(state, emission, resolved.density);
    const double source = ScalarSectorSource(state, emission, resolved);
    amplitudes.push_back(recoil * charge_phase * source / std::sqrt(state.xi));
  }
  return amplitudes;
}

// Apply the common EPA normalization and ordered sector weights to amplitudes
bool ApplyEPAAmplitudeWeights(std::vector<std::complex<double>> &amplitudes, const double common,
                              const EPASectorWeights &weights) {
  if (common < 0.0 || weights.amplitude_fraction.empty()) { return false; }

  std::vector<std::complex<double>> weighted;
  weighted.reserve(amplitudes.size() * weights.amplitude_fraction.size());
  for (const auto &amplitude : amplitudes) {
    for (const double factor : weights.amplitude_fraction) {
      const std::complex<double> value = common * factor * amplitude;
      weighted.push_back(value);
    }
  }
  amplitudes.swap(weighted);
  return true;
}

// Compute the transverse EPA density eigenvalues of a charged spin-half emitter
// eigenvalues = (n_parallel,n_perpendicular)
TransversePhotonFlux ElasticSpinHalfFluxTransverse(double xi, double t, double pt, double mass,
                                                   double electric_form_factor, double magnetic_form_factor,
                                                   double charge) {
  if (!std::isfinite(xi) || !std::isfinite(t) || t > 0.0 || !std::isfinite(pt)) { return {}; }
  const double pt2 = math::pow2(pt);
  if (!ValidMomentumFraction(xi)) { return {}; }

  const double xi2   = math::pow2(xi);
  const double mass2 = math::pow2(mass);
  const double Q2    = -t;
  const double denom = pt2 + xi2 * mass2;
  if (!(denom > 0.0)) { return {}; }
  const double electric_part =
      (4.0 * mass2 * math::pow2(electric_form_factor) + Q2 * math::pow2(magnetic_form_factor)) / (4.0 * mass2 + Q2);
  const double magnetic_part      = math::pow2(magnetic_form_factor);
  // Cancel pt2 before evaluating the finite magnetic density at pt = 0
  const double delta              = pt2 / denom;
  const double common             = math::pow2(charge) * 16.0 * math::PIPI * qed::alpha_QED() / (math::PI * xi * denom);
  const double electric           = common * (1.0 - xi) * delta * electric_part;
  const double magnetic_per_state = common * (xi2 / 4.0) * magnetic_part;
  return {electric + magnetic_per_state, magnetic_per_state};
}

// Compute the elastic transverse EPA density of a proton
TransversePhotonFlux CohFluxTransverse(double xi, double t, double pt, const form::ParamStore &structure) {
  if (!std::isfinite(t) || t > 0.0) { return {}; }
  const double Q2 = -t;
  return ElasticSpinHalfFluxTransverse(xi, t, pt, PDG::mp, form::G_E(Q2, structure), form::G_M(Q2, structure), 1.0);
}

// Compute the scalar coherent EPA flux
// n_gamma = n_parallel+n_perpendicular
double CohFlux(double xi, double t, double pt, const form::ParamStore &structure) {
  return CohFluxTransverse(xi, t, pt, structure).Trace();
}

// Compute the inelastic transverse EPA density of a dissociated proton
TransversePhotonFlux IncohFluxTransverse(double xi, double t, double pt, double M2, const form::ParamStore &structure) {
  constexpr double proton_mass2 = math::pow2(PDG::mp);
  if (!std::isfinite(xi) || !std::isfinite(t) || t > 0.0 || !std::isfinite(pt) || !std::isfinite(M2)) { return {}; }
  const double pt2 = math::pow2(pt);
  if (!ValidMomentumFraction(xi)) { return {}; }

  const double xi2   = math::pow2(xi);
  const double Q2    = -t;
  const double denom = Q2 + M2 - proton_mass2;
  if (!(denom > 0.0)) { return {}; }
  const double xbj = Q2 / denom;
  if (!(xbj > 0.0) || xbj > 1.0) { return {}; }
  const double delta_denom = pt2 + xi * (M2 - proton_mass2) + xi2 * proton_mass2;
  if (!(delta_denom > 0.0)) { return {}; }

  const double delta    = pt2 / delta_denom;
  const double common   = 16.0 * math::PIPI * qed::alpha_QED() / (math::PI * xi * delta_denom);
  const double electric = common * (1.0 - xi) * delta * form::F2xQ2(xbj, Q2, structure) / denom;
  const double magnetic_per_state =
      common * (xi2 / (4.0 * math::pow2(xbj))) * (2.0 * xbj * form::F1xQ2(xbj, Q2, structure)) / denom;
  return {electric + magnetic_per_state, magnetic_per_state};
}

// Compute the scalar inelastic EPA flux
// n_gamma^inel = n_parallel+n_perpendicular
double IncohFlux(double xi, double t, double pt, double M2, const form::ParamStore &structure) {
  return IncohFluxTransverse(xi, t, pt, M2, structure).Trace();
}

// Compute the momentum-transfer integrated Drees-Zeppenfeld photon flux
// n_DZ(x) = alpha[1+(1-x)^2] I(Qmin^2)/(2 pi x)
double DZFlux(double x) {
  if (!ValidMomentumFraction(x)) { return 0.0; }
  const double Q2min = math::pow2(PDG::mp * x) / (1.0 - x);
  const double flux  = qed::alpha_QED() / (2.0 * math::PI * x) * (1.0 + math::pow2(1.0 - x)) * DZIntegralBracket(Q2min);
  return std::isfinite(flux) && flux > 0.0 ? flux : 0.0;
}

// Compute the coherent or inclusive photon flux for one generated forward leg
double ForwardPhotonFlux(const ForwardLegState &state) { return ForwardPhotonFluxTransverse(state).Trace(); }

// Compute one selected nuclear photon-emission density of a nuclear leg
TransversePhotonFlux NuclearPhotonFluxTransverse(const ForwardLegState &state, const nuclear::CoherenceType emission) {
  if (!state.IsNuclear() || state.upc_model == nullptr || !state.has_xi || !ValidMomentumFraction(state.xi)) {
    return {};
  }
  const nuclear::MPhoton* photon = state.upc_model->Photon(state.Index());
  if (photon == nullptr) { return {}; }
  try {
    nuclear::PhotonDensity density;
    if (const auto* bank = state.upc_model->Bank(state.Index()); bank != nullptr) {
      const M3Vec q = nuclear::RestTransfer(state.incoming, state.transfer, state.emitter.mass);
      density       = photon->Density(emission, state.xi, state.t, state.qt, q, bank);
    } else {
      density = photon->Density(emission, state.xi, state.t, state.qt);
    }
    return {density.parallel, density.perpendicular};
  } catch (const AmplitudeFailure&) { throw; } catch (const std::exception& error) {
    throw AmplitudeFailure(std::string("NuclearPhotonFluxTransverse: ") + error.what());
  }
}

// Compute the transverse density used to normalize the active photon sources
TransversePhotonFlux ForwardPhotonFluxTransverse(const ForwardLegState& state) {
  if (!state.has_xi || !ValidMomentumFraction(state.xi)) { return {}; }
  if (state.IsNuclear()) {
    const auto sector = state.upc_model != nullptr ? state.upc_model->SourceSector(state.Index(), state.emission)
                                                  : state.emission;
    return NuclearPhotonFluxTransverse(state, sector);
  }
  if (state.IsExcited()) { return IncohFluxTransverse(state.xi, state.t, state.qt, state.mass2, state.structure); }
  if (std::abs(state.emitter.pdg) == PDG::PDG_p) {
    return CohFluxTransverse(state.xi, state.t, state.qt, state.structure);
  }
  const int pdg = std::abs(state.emitter.pdg);
  if ((pdg == 11 || pdg == 13 || pdg == 15) && state.emitter.spinX2 == 1) {
    return ElasticSpinHalfFluxTransverse(state.xi, state.t, state.qt, state.emitter.mass, 1.0, 1.0,
                                         state.emitter.chargeX3 / 3.0);
  }
  return {};
}

// Compute the exact EPA conversion from the pp phase space to dx1 dx2 and hard
// flux
//
// The exact light-cone Jacobian replaces dx1 dx2 = dM2 dY / s, while
// (q1 + q2)^2 replaces x1 x2 s in the on-shell gamma-gamma Moller flux
// [REFERENCE: Luszczak, Schafer and Szczurek, JHEP 05 (2018) 064, Eq. (2.9)]
double ExactktEPAPhaseSpaceFactor(const gra::LORENTZSCALAR &lts) {
  if (!lts.has_xi1 || !lts.has_xi2 || !ValidMomentumFraction(lts.xi1) || !ValidMomentumFraction(lts.xi2) ||
      !(lts.s > 0.0) || lts.pfinal.size() < 3) {
    return 0.0;
  }

  const M4Vec &forward1 = lts.pfinal[1];
  const M4Vec &forward2 = lts.pfinal[2];
  if (!(forward1.E() > 0.0) || !(forward2.E() > 0.0)) { return 0.0; }

  const double beam_plus    = lts.pbeam1.LightconePos();
  const double beam_minus   = lts.pbeam2.LightconeNeg();
  const double retained1    = 1.0 - lts.xi1;
  const double retained2    = 1.0 - lts.xi2;
  const double beam_product = beam_plus * beam_minus;
  if (!(beam_product > 0.0) || !(retained1 > 0.0) || !(retained2 > 0.0)) { return 0.0; }

  const double beam_mass1      = lts.beam1.mass;
  const double beam_mass2      = lts.beam2.mass;
  const double forward_mass1sq = lts.forward_mass2[0];
  const double forward_mass2sq = lts.forward_mass2[1];
  if (!std::isfinite(beam_mass1) || !std::isfinite(beam_mass2) || beam_mass1 < 0.0 || beam_mass2 < 0.0 ||
      !std::isfinite(forward_mass1sq) || !std::isfinite(forward_mass2sq) || forward_mass1sq < 0.0 ||
      forward_mass2sq < 0.0) {
    return 0.0;
  }
  const double mt1sq = forward_mass1sq + forward1.Pt2();
  const double mt2sq = forward_mass2sq + forward2.Pt2();
  const double lightcone_determinant =
      beam_product - mt1sq * mt2sq / (beam_product * math::pow2(retained1) * math::pow2(retained2));
  if (!std::isfinite(lightcone_determinant) || std::abs(lightcone_determinant) <= 1.0e-24) { return 0.0; }
  const double dx_jacobian = 1.0 / std::abs(lightcone_determinant);

  const double velocity_difference = std::abs(forward1.Pz() / forward1.E() - forward2.Pz() / forward2.E());
  const double beta                = kinematics::beta12(lts.s, beam_mass1, beam_mass2);
  const double hard_mass2          = (lts.q1 + lts.q2).M2();
  if (!(velocity_difference > 0.0) || !(beta > 0.0) || !(hard_mass2 > 0.0)) { return 0.0; }

  const double factor =
      2.0 * lts.s * beta * forward1.E() * forward2.E() * velocity_difference * dx_jacobian / hard_mass2;
  return std::isfinite(factor) && factor > 0.0 ? factor : 0.0;
}

// Convert the incoming proton flux to the massless collinear photon flux
// F_pp/F_gg = s beta_pp/shat
// [REFERENCE: Particle Data Group, Review of Particle Physics, Kinematics]
double CollinearPhotonPhaseSpaceFactor(const gra::LORENTZSCALAR &lts) {
  if (!ValidMomentumFraction(lts.x1) || !ValidMomentumFraction(lts.x2) || !(lts.s > 0.0) || !(lts.s_hat > 0.0)) {
    return 0.0;
  }

  const double beam_mass1 = lts.beam1.mass;
  const double beam_mass2 = lts.beam2.mass;
  if (!std::isfinite(beam_mass1) || !std::isfinite(beam_mass2) || beam_mass1 < 0.0 || beam_mass2 < 0.0) { return 0.0; }
  const double beta   = kinematics::beta12(lts.s, beam_mass1, beam_mass2);
  const double factor = lts.s * beta / lts.s_hat;
  return std::isfinite(factor) && factor > 0.0 ? factor : 0.0;
}

// Set the photon hard scales and alpha_s from the configured PDF
bool SetPhotonAlphaS(gra::LORENTZSCALAR &lts) {
  const double Q2 = lts.s_hat / 4.0;
  lts.muF         = Q2 > 0.0 ? math::msqrt(Q2) : 0.0;
  lts.muR         = lts.muF;
  lts.scalup      = lts.muF;
  lts.alphaQCD    = 0.0;
  if (!(Q2 > 0.0) || !std::isfinite(Q2)) { return false; }

  try {
    if (!lts.GlobalPdfPtr->hasAlphaS() || !lts.GlobalPdfPtr->inRangeQ2(Q2)) { return false; }
    const double alpha_s = lts.GlobalPdfPtr->alphasQ2(Q2);
    if (!(alpha_s > 0.0) || !std::isfinite(alpha_s)) { return false; }
    lts.alphaQCD = alpha_s;
  } catch (...) {
    lts.alphaQCD = 0.0;
    return false;
  }
  return true;
}

// Apply transverse EPA currents with one supplied outer phase-space conversion
EPAWeight ApplyktEPAcurrents(double amp2, gra::LORENTZSCALAR &lts, const double phase_space) {
  // kT-EPA carries physical non-collinear forward particles into the LHE record
  lts.exact_forward_photon_kinematics = true;
  if (!lts.has_xi1 || !lts.has_xi2 || !ValidMomentumFraction(lts.xi1) ||
      !ValidMomentumFraction(lts.xi2)) {
    ZeroHelicityAmplitudes(lts);
    return {};
  }

  const ForwardLegState upper      = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState lower      = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const double          gammaflux1 = ForwardPhotonFlux(upper);
  const double          gammaflux2 = ForwardPhotonFlux(lower);

  if (!PositiveFinite(gammaflux1) || !PositiveFinite(gammaflux2)) {
    ZeroHelicityAmplitudes(lts);
    return {};
  }

  if (!PositiveFinite(phase_space)) {
    ZeroHelicityAmplitudes(lts);
    return {};
  }

  // Combine all
  const double fluxes = gammaflux1 * gammaflux2 * phase_space;
  if (fluxes <= 0.0) {
    ZeroHelicityAmplitudes(lts);
    return {};
  }
  const double tot = fluxes * amp2;

  // --------------------------------------------------------------------
  // Apply fluxes to helicity amplitudes
  EPAWeight result;
  result.amp2            = tot;
  result.amplitude_scale = msqrt(fluxes);
  result.sector          = ResolveInclusiveEmissionSectors(lts, upper, lower);
  if (!ApplyEPAAmplitudeWeights(lts.hamp, result.amplitude_scale, result.sector)) {
    ZeroHelicityAmplitudes(lts);
    return {};
  }
  // --------------------------------------------------------------------

  // ** Save for HepMC output **
  lts.id1     = PDG::PDG_gamma;
  lts.id2     = PDG::PDG_gamma;
  lts.pdf_xf1 = lts.xi1 * gammaflux1;
  lts.pdf_xf2 = lts.xi2 * gammaflux2;
  lts.muF     = math::msqrt(std::max(lts.s_hat, 0.0)) / 2.0;
  lts.muR     = lts.muF;
  lts.scalup  = lts.muF;

  return result;
}

// Apply non-collinear EPA fluxes at cross section level
// Use with full 2 -> N kinematics
EPAWeight ApplyktEPAfluxes(double amp2, gra::LORENTZSCALAR &lts) {
  return ApplyktEPAcurrents(amp2, lts, ExactktEPAPhaseSpaceFactor(lts));
}

// Apply Gamma-Gamma collinear Drees-Zeppenfeld (coherent flux) at cross section
// level Use with collinear kinematics
double ApplyDZfluxes(double amp2, gra::LORENTZSCALAR &lts) {
  // Collinear DZ photons use the ordinary partonic LHE view
  lts.exact_forward_photon_kinematics = false;
  if (!ValidMomentumFraction(lts.x1) || !ValidMomentumFraction(lts.x2)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }

  // Evaluate gamma pdfs
  const double f1 = DZFlux(lts.x1);
  const double f2 = DZFlux(lts.x2);
  if (!PositiveFinite(f1) || !PositiveFinite(f2)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }

  const double phasespace = CollinearPhotonPhaseSpaceFactor(lts);
  if (!PositiveFinite(phasespace)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  const double fluxes = f1 * f2 * phasespace;
  if (fluxes <= 0.0) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  const double tot = fluxes * amp2;

  // --------------------------------------------------------------------
  // Apply fluxes to helicity amplitudes
  if (!ApplyEPAAmplitudeWeights(lts.hamp, msqrt(fluxes), EPASectorWeights{})) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  // --------------------------------------------------------------------

  // ** Save for HepMC output **
  lts.id1     = PDG::PDG_gamma;
  lts.id2     = PDG::PDG_gamma;
  lts.pdf_xf1 = lts.x1 * f1;
  lts.pdf_xf2 = lts.x2 * f2;
  lts.muF     = math::msqrt(std::max(lts.s_hat, 0.0)) / 2.0;
  lts.muR     = lts.muF;
  lts.scalup  = lts.muF;

  return tot;
}

// Apply Gamma-Gamma LUX-pdf (use at \mu > 10 GeV) at cross section level
// Use with collinear kinematics
// [REFERENCE: Manohar, Nason, Salam and Zanderighi, arXiv:1708.01256]
double ApplyLUXfluxes(double amp2, gra::LORENTZSCALAR &lts) {
  // Collinear LUX photons use the ordinary partonic LHE view
  lts.exact_forward_photon_kinematics = false;
  if (!ValidMomentumFraction(lts.x1) || !ValidMomentumFraction(lts.x2)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }

  // pdf factorization scale
  const double Q2 = lts.s_hat / 4.0;
  if (!(Q2 > 0.0)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  if (!lts.GlobalPdfPtr->inRangeXQ2(lts.x1, Q2) || !lts.GlobalPdfPtr->inRangeXQ2(lts.x2, Q2)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }

  // Evaluate gamma pdfs
  double f1 = 0.0;
  double f2 = 0.0;
  try {
    // Divide x out
    f1 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x1, Q2) / lts.x1;
    f2 = lts.GlobalPdfPtr->xfxQ2(PDG::PDG_gamma, lts.x2, Q2) / lts.x2;
  } catch (...) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  if (!PositiveFinite(f1) || !PositiveFinite(f2)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }

  const double phasespace = CollinearPhotonPhaseSpaceFactor(lts);
  if (!PositiveFinite(phasespace)) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  const double fluxes = f1 * f2 * phasespace;
  if (fluxes <= 0.0) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  const double tot = fluxes * amp2;

  // --------------------------------------------------------------------
  // Apply fluxes to helicity amplitudes
  if (!ApplyEPAAmplitudeWeights(lts.hamp, msqrt(fluxes), EPASectorWeights{})) {
    ZeroHelicityAmplitudes(lts);
    return 0.0;
  }
  // --------------------------------------------------------------------

  // ** Save for HepMC output **
  lts.id1     = PDG::PDG_gamma;
  lts.id2     = PDG::PDG_gamma;
  lts.pdf_xf1 = lts.x1 * f1;
  lts.pdf_xf2 = lts.x2 * f2;
  lts.muF     = math::msqrt(Q2);
  lts.muR     = lts.muF;
  lts.scalup  = lts.muF;

  return tot;
}

// Compute the elementary elastic nucleon state in the nuclear transfer Breit frame
ForwardLegState ElementaryPhotoTarget(const LORENTZSCALAR &lts, const ForwardLegState &state) {
  ForwardLegState elementary = state;
  if (!state.IsNuclear()) { return elementary; }
  const int proton_pdg   = state.emitter.pdg < 0 ? -PDG::PDG_p : PDG::PDG_p;
  elementary.emitter     = lts.PDG.FindByPDG(proton_pdg);
  // Equal-mass nucleon momenta in the transfer Breit frame preserve the Ward identities
  const M4Vec average = state.incoming - state.transfer * ((state.incoming * state.transfer) / state.t);
  const double mass2 = elementary.emitter.mass * elementary.emitter.mass;
  const double norm2 = average.M2();
  if (!(state.t < 0.0) || !(norm2 > 0.0)) { throw AmplitudeFailure("ElementaryPhotoTarget: invalid nuclear transfer"); }
  const M4Vec center = average * std::sqrt((mass2 - state.t / 4.0) / norm2);
  elementary.incoming    = center + state.transfer / 2.0;
  elementary.outgoing    = center - state.transfer / 2.0;
  elementary.mass2       = elementary.emitter.mass * elementary.emitter.mass;
  elementary.final_state = ForwardFinalState::Elastic;
  elementary.is_nuclear  = false;
  elementary.upc_model.reset();
  return elementary;
}

// Sum physical photon directions and preserve each nuclear target current for screening
void SumPhotoTerms(LORENTZSCALAR &lts, std::vector<nuclear::PhotoTerm> terms) {
  lts.screening.photo.term.clear();
  for (auto &current : lts.screening.photo.current) { current.reset(); }
  const auto upper = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const auto lower = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const std::array<bool, 2> active = {SupportsPhotoTarget(lower), SupportsPhotoTarget(upper)};
  std::array<std::size_t, 2> channel_count{};
  if (lts.upc_model == nullptr) {
    for (const auto &direction : indices(active)) { channel_count[direction] = active[direction] ? 1 : 0; }
  } else {
    for (const auto &channel : nuclear::PhotoChannels(*upper.upc_model)) {
      ++channel_count[channel.direction == nuclear::PhotoDirection::Upper ? 0 : 1];
    }
  }
  const bool split = SplitPhotoDirections(lts, upper, lower);
  std::array<std::vector<std::complex<double>>, 2> direction_sum;
  for (const auto &term : terms) {
    const std::size_t direction = term.direction == nuclear::PhotoDirection::Upper ? 0 : 1;
    if (term.amplitude.empty() || channel_count[direction] == 0) { continue; }
    if (term.amplitude.size() % channel_count[direction] != 0) {
      throw AmplitudeFailure("SumPhotoTerms: inconsistent nuclear sector layout");
    }
    auto &sum = direction_sum[direction];
    if (sum.empty()) {
      sum = term.amplitude;
    } else {
      if (sum.size() != term.amplitude.size()) { throw AmplitudeFailure("SumPhotoTerms: incompatible helicity rows"); }
      gra::AddScaled(sum, term.amplitude, 1.0);
    }
  }
  std::size_t hard_count = 0;
  for (const auto &direction : indices(direction_sum)) {
    if (!direction_sum[direction].empty()) {
      hard_count = direction_sum[direction].size() / channel_count[direction];
      break;
    }
  }
  if (hard_count == 0) { throw AmplitudeFailure("SumPhotoTerms: no physical photon direction"); }

  lts.hamp.clear();
  if (!split) {
    for (const auto &direction : indices(channel_count)) {
      if (active[direction] && channel_count[direction] != 1) { throw AmplitudeFailure("SumPhotoTerms: unresolved coherent direction layout"); }
    }
    lts.hamp.assign(hard_count, 0.0);
    for (const auto &h : indices(lts.hamp)) {
      for (const auto &direction : indices(direction_sum)) {
        if (!direction_sum[direction].empty()) { lts.hamp[h] += direction_sum[direction][h]; }
      }
    }
    flux::ConfigurePhotoLayout(lts);
  } else {
    const std::size_t sector_count = channel_count[0] + channel_count[1];
    lts.hamp.reserve(hard_count * sector_count);
    for (std::size_t h = 0; h < hard_count; ++h) {
      for (const auto &direction : indices(direction_sum)) {
        for (std::size_t sector = 0; sector < channel_count[direction]; ++sector) {
          const auto &amplitude = direction_sum[direction];
          lts.hamp.push_back(amplitude.empty() ? std::complex<double>(0.0, 0.0) : amplitude[h * channel_count[direction] + sector]);
        }
      }
    }
    flux::ConfigurePhotoLayout(lts);

    lts.screening.photo.term.reserve(terms.size());
    for (auto &term : terms) {
      const std::size_t direction          = term.direction == nuclear::PhotoDirection::Upper ? 0U : 1U;
      const std::size_t direction_channels = channel_count[direction];
      if (direction_channels == 0 || term.amplitude.size() != hard_count * direction_channels) { throw AmplitudeFailure("SumPhotoTerms: photonuclear term layout changed"); }
      nuclear::PhotoTerm resolved;
      resolved.direction = term.direction;
      resolved.amplitude.assign(hard_count * sector_count, 0.0);
      resolved.photo_current   = std::move(term.photo_current);
      const std::size_t offset = direction == 0 ? 0 : channel_count[0];
      for (std::size_t h = 0; h < hard_count; ++h) {
        for (std::size_t channel = 0; channel < direction_channels; ++channel) { resolved.amplitude[h * sector_count + offset + channel] = term.amplitude[h * direction_channels + channel]; }
      }
      lts.screening.photo.term.push_back(std::move(resolved));
    }
  }
}

}  // namespace flux
}  // namespace gra
