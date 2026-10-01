// Shared pp and UPC screening of process amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <span>
#include <vector>

#include "Graniitti/Eikonal/MProtonScreen.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {

namespace {

// Limit central event caches to one complete Born and screening evaluation
class CentralCacheGuard {
 public:
  // Start a fresh cache generation for this amplitude evaluation
  explicit CentralCacheGuard(MAmplitudeEventState &state) : state(state) { state.BeginCentral(); }

  // Invalidate all cached central data when the evaluation leaves scope
  ~CentralCacheGuard() { state.AbortCentral(); }

  // Prevent a second guard from owning the same cache lifetime
  CentralCacheGuard(const CentralCacheGuard &)            = delete;
  CentralCacheGuard &operator=(const CentralCacheGuard &) = delete;

 private:
  MAmplitudeEventState &state;
};

// Compute whether one beam removes the need for a hadronic survival loop
bool HasChargedLeptonBeam(const LORENTZSCALAR &lts) {
  const auto charged_lepton = [](const MParticle &beam) {
    const int pdg = std::abs(beam.pdg);
    return pdg == 11 || pdg == 13 || pdg == 15;
  };
  return charged_lepton(lts.beam1) || charged_lepton(lts.beam2);
}

// Compute the checked photonuclear channel layout stored by one amplitude
std::vector<nuclear::PhotoChannel> UPCPhotoChannels(const ScreeningLayout &layout) {
  if (!layout.photo_sector_resolved) { return {}; }
  if (layout.photo_channel_count == 0 || layout.photo_channel_count > layout.photo_channel.size()) {
    throw AmplitudeFailure("UPCPhotoChannels: invalid channel count");
  }
  return {layout.photo_channel.begin(), layout.photo_channel.begin() + layout.photo_channel_count};
}

// Decode one explicit EPA source sector stored in screening metadata
nuclear::CoherenceType UPCEmissionSector(const std::uint8_t type) {
  if (type == static_cast<std::uint8_t>(nuclear::CoherenceType::Coherent)) { return nuclear::CoherenceType::Coherent; }
  if (type == static_cast<std::uint8_t>(nuclear::CoherenceType::Incoherent)) {
    return nuclear::CoherenceType::Incoherent;
  }
  throw AmplitudeFailure("UPCEmissionSector: unresolved EPA sector");
}

// Map process screening states to the nuclear convolution layout
nuclear::ScreenLayout MakeUPCScreenLayout(const ScreeningLayout &event_layout, const nuclear::MUPC &upc) {
  nuclear::ScreenLayout layout;
  layout.photo = UPCPhotoChannels(event_layout);
  if (event_layout.nuclear_type == nuclear::ScreenType::Fusion) {
    layout.fusion.rows = event_layout.epa_rows_per_sector;
    if (!event_layout.epa_sector_resolved) {
      throw AmplitudeFailure("UPCScreenLayout: unresolved photon-fusion sectors");
    }
    for (std::size_t leg = 0; leg < 2; ++leg) {
      const std::size_t count = event_layout.epa_sector_count[leg];
      if (count == 0 || count > event_layout.epa_sector_type[leg].size()) {
        throw AmplitudeFailure("UPCScreenLayout: invalid photon-fusion sector count");
      }
      for (std::size_t sector = 0; sector < count; ++sector) {
        layout.fusion.sector[leg].push_back(UPCEmissionSector(event_layout.epa_sector_type[leg][sector]));
      }
    }
    if (upc.HasSamples()) { layout.type = nuclear::ScreenType::Fusion; }
  } else if (event_layout.nuclear_type == nuclear::ScreenType::Photo) {
    if (layout.photo.empty()) {
      const auto &param = upc.Param();
      for (std::size_t leg = 0; leg < 2; ++leg) {
        if (upc.Nucleus(static_cast<int>(leg + 1)) != nullptr &&
            (param.emission[leg] != nuclear::CoherenceType::Coherent ||
             param.target[leg] != nuclear::CoherenceType::Coherent)) {
          throw AmplitudeFailure("UPCScreenLayout: unresolved photonuclear final sector");
        }
      }
      layout.fixed.valid = true;
    } else if (upc.HasSamples()) {
      layout.type = nuclear::ScreenType::Photo;
    }
  }
  return layout;
}

// Extract the amplitude rows while leaving hard color tags outside screening
std::vector<std::vector<std::complex<double>>> HardColorAmplitudes(const std::vector<mg5helas::HardColorFlow> &flows) {
  std::vector<std::vector<std::complex<double>>> amplitude;
  amplitude.reserve(flows.size());
  for (const auto &flow : flows) { amplitude.push_back(flow.amplitudes); }
  return amplitude;
}

// Compare the ordered external color assignments and amplitude row counts
bool SameHardColorLayout(const std::vector<mg5helas::HardColorFlow> &born,
                         const std::vector<mg5helas::HardColorFlow> &shifted) {
  if (born.size() != shifted.size()) { return false; }
  for (const auto &flow : indices(born)) {
    if (born[flow].channel != shifted[flow].channel ||
        born[flow].amplitudes.size() != shifted[flow].amplitudes.size() ||
        born[flow].external.size() != shifted[flow].external.size()) {
      return false;
    }
    for (const auto &leg : indices(born[flow].external)) {
      if (born[flow].external[leg].color != shifted[flow].external[leg].color ||
          born[flow].external[leg].anticolor != shifted[flow].external[leg].anticolor) {
        return false;
      }
    }
  }
  return true;
}

// Move the amplitude data produced at the current kinematic point into the UPC convolution
nuclear::ScreenPoint TakeUPCScreenPoint(LORENTZSCALAR &lts, std::vector<std::vector<std::complex<double>>> flow) {
  nuclear::ScreenPoint point;
  point.amplitude.assign(lts.hamp.begin(), lts.hamp.end());
  point.flow     = std::move(flow);
  const bool emd = lts.upc_model->Param().reaction && lts.upc_model->Param().additional_emd;
  for (const auto& leg : indices(point.transfer)) {
    const bool upper = leg == 0;
    const auto transfer = emd ? ResolveForwardLegState(lts, upper ? ForwardBeamLeg::Upper : ForwardBeamLeg::Lower).transfer
                              : (upper ? lts.q1 : lts.q2);
    point.transfer[leg] = nuclear::RestTransfer(upper ? lts.pbeam1 : lts.pbeam2, transfer,
                                               upper ? lts.beam1.mass : lts.beam2.mass);
  }
  point.photo_current.swap(lts.screening.photo.current);
  point.photo_terms.swap(lts.screening.photo.term);
  return point;
}

// Require agreement between the selected representation and its event payload
void RequireAmplitudeRepresentation(const ScreeningMetadata                        &metadata,
                                    const std::optional<ProtonGoodWalkerAmplitude> &payload) {
  const bool expects_pair = metadata.amplitude_type == ScreeningAmplitudeType::GoodWalker;
  if (expects_pair != payload.has_value()) {
    throw AmplitudeFailure("RequireAmplitudeRepresentation: Good Walker amplitude is inconsistent");
  }
}

// Compute the event Good Walker class selected by the forward systems
GWFinalClass EventGoodWalkerClass(const LORENTZSCALAR &lts) {
  if (lts.excite1 && lts.excite2) { return GWFinalClass::DoubleDissociation; }
  if (lts.excite1) { return GWFinalClass::SingleDissociation1; }
  if (lts.excite2) { return GWFinalClass::SingleDissociation2; }
  return GWFinalClass::Elastic;
}


}  // namespace

// Reconstruct the complete event and evaluate its selected hard amplitude
//
bool MProcess::ComputeScreeningNode(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) {
  evaluation_status = mg5helas::EvaluationStatus::Success;
  if (!LoopKinematics(p1p, p2p)) {
    return false;
  }
  state.lts.proton_good_walker.reset();
  state.lts.screening.BeginNode();
  state.lts.hamp.layout.Clear();
  const double amp2 = EvaluateBareAmplitude();
  if (evaluation_status == mg5helas::EvaluationStatus::KinematicsFailure) {
    // Consume the rejected node without changing the status of the completed convolution
    evaluation_status = mg5helas::EvaluationStatus::Success;
    return false;
  }
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    throw AmplitudeFailure("MProcess::EvaluateScreeningNode: amplitude evaluation failed");
  }
  ValidateAmplitude(amp2);
  if (state.lts.hamp.metadata.amplitude_type != ScreeningAmplitudeType::GoodWalker && state.lts.hamp.empty()) {
    throw AmplitudeFailure("MProcess::EvaluateScreeningNode: empty amplitude");
  }
  return true;
}

// Screen one Good Walker pair-space amplitude with the pp eikonal
//
double MProcess::ScreenPPGoodWalker(ProtonGoodWalkerAmplitude born, const ScreeningMetadata &metadata,
                                    const std::array<double, 2> &upper_transverse,
                                    const std::array<double, 2> &lower_transverse,
                                    const MEikonal::LoopConst   &loop_const) {
  eikonal::MProtonGoodWalkerScreen screen(born, metadata, loop_const, eikonal.SoftModelHandle(),
                                          eikonal.GetChannelCount());
  ScreeningLoopGuard               loop_guard(*this);
  const std::size_t                n_phi = loop_const.node_weight.size_col();

  for (const auto &i : indices(loop_const.kt2)) {
    for (std::size_t j = 0; j < n_phi; ++j) {
      const double                kt_x          = loop_const.kt_x[i][j];
      const double                kt_y          = loop_const.kt_y[i][j];
      const std::array<double, 2> upper_shifted = {upper_transverse[0] - kt_x, upper_transverse[1] - kt_y};
      const std::array<double, 2> lower_shifted = {lower_transverse[0] + kt_x, lower_transverse[1] + kt_y};
      if (!EvaluateScreeningNode(upper_shifted, lower_shifted)) { continue; }
      if (!state.lts.proton_good_walker) {
        throw AmplitudeFailure(
            "MProcess::ScreenPPGoodWalker: shifted Good Walker source is "
            "absent");
      }
      screen.Add(i, j, *state.lts.proton_good_walker);
    }
  }

  if (!RestoreScreeningBornState()) {
    throw AmplitudeFailure(
        "MProcess::ScreenPPGoodWalker: failed to restore Born "
        "kinematics");
  }
  auto *metrics = state.lts.qmetrics.Active() ? &state.lts.qmetrics : nullptr;
  if (metrics != nullptr) { metrics->ClearEvent(); }
  state.lts.hamp                = screen.Result(metrics);
  const AmplitudeWeights weight = ProcPtr.NormalizeScreenedAmplitude(state.lts.hamp, true);

  PostScreeningAmplitude(weight.helicity_weight);

  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    throw AmplitudeFailure("MProcess::ScreenPPGoodWalker: amplitude finalization failed");
  }
  state.lts.proton_good_walker.emplace(std::move(born));
  loop_guard.Accept();
  return weight.total_weight;
}

// Compute |A_Born + integral d2kT S(kT) A(p1T-kT,p2T+kT)|^2 with process normalization
double MProcess::ScreenedAmplitudeSquared(bool include_screening) {
  if (state.lts.screening.active) {
    // This indicates an internal event-flow violation and should never occur during sampling
    throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: nested screening loop state");
  }
  CentralCacheGuard       cache_guard(state.lts.amplitude);
  const ScreeningMetadata metadata = state.lts.hamp.metadata;
  if (state.lts.proton_good_walker) {
    throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: stale Good Walker amplitude");
  }

  // Elastic scattering always uses the physical Coulomb and nuclear matrix
  if (metadata.amplitude_type == ScreeningAmplitudeType::ElasticCNI) {
    if (!std::isfinite(metadata.amplitude_normalization) || metadata.amplitude_normalization <= 0.0) {
      throw AmplitudeFailure(
          "MProcess::ScreenedAmplitudeSquared: invalid elastic initial-state "
          "average");
    }
    if (!eikonal.HasElasticCNI()) {
      throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: elastic CNI is not initialized");
    }
    if (state.lts.pfinal.size() < 3) {
      throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: elastic final state is missing");
    }

    ProtonHelicityMatrix matrix;
    try {
      matrix = eikonal.PhysicalElasticHelicityMatrix(state.lts.pbeam1, state.lts.pbeam2, state.lts.pfinal[1],
                                                     state.lts.pfinal[2]);
    } catch (const std::exception &error) {
      throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: elastic CNI evaluation failed: " +
                             std::string(error.what()));
    } catch (...) { throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: elastic CNI evaluation failed"); }

    state.lts.hamp.assign(matrix.begin(), matrix.end());
    const AmplitudeWeights weight = ProcPtr.NormalizeAmplitude(state.lts.hamp);
    PostScreeningAmplitude(weight.helicity_weight);
    if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
      throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: amplitude finalization failed");
    }
    return weight.total_weight;
  }

  const bool excitation = state.lts.upc_excitation != nullptr;
  const bool convolution = state.lts.upc_model != nullptr ? (excitation || state.lts.upc_model->HadronicConvolution())
                                                          : state.screening && !HasChargedLeptonBeam(state.lts);

  // Compute the bare amplitude when convolution is disabled for this evaluation
  if ((!include_screening && !excitation) || !convolution) {
    EvaluateAmplitude();
    RequireAmplitudeRepresentation(metadata, state.lts.proton_good_walker);
    const ScreeningMetadata evaluated_metadata = state.lts.hamp.metadata;
    const ScreeningLayout   evaluated_layout   = state.lts.hamp.layout;
    if (state.lts.upc_event != nullptr &&
        (evaluated_layout.photo_sector_resolved ||
         (evaluated_layout.nuclear_type == nuclear::ScreenType::Fusion && state.lts.upc_event->HasSamples()))) {
      try {
        const auto&               upc    = *state.lts.upc_event;
        const auto                layout = MakeUPCScreenLayout(evaluated_layout, upc);
        const nuclear::MUPCScreen screen(
            upc, layout, TakeUPCScreenPoint(state.lts, HardColorAmplitudes(state.lts.hard_color_flows)));
        return FinalizeUPCAmplitude(screen.Result(), layout, evaluated_metadata);
      } catch (const AmplitudeFailure&) { throw; } catch (const std::exception& error) {
        throw AmplitudeFailure(std::string("MProcess::ScreenedAmplitudeSquared: ") + error.what());
      }
    }
    if (evaluated_layout.photo_sector_resolved) {
      try {
        state.lts.hamp = nuclear::CombinePhotoChannels(state.lts.hamp, UPCPhotoChannels(evaluated_layout));
      } catch (const std::exception &error) {
        throw AmplitudeFailure(std::string("MProcess::ScreenedAmplitudeSquared: ") + error.what());
      }
      state.lts.hamp.Configure(evaluated_metadata);
    }
    const AmplitudeWeights weight = ProcPtr.NormalizeAmplitude(state.lts.hamp);
    if (state.lts.upc_event != nullptr) {
      const nuclear::ScreenLayout layout = MakeUPCScreenLayout(evaluated_layout, *state.lts.upc_event);
      SetUPCFinalWeights(nuclear::FinalSectorWeights(layout, weight.helicity_weight));
    }
    PostScreeningAmplitude(weight.helicity_weight);
    if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
      throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: amplitude finalization failed");
    }
    return weight.total_weight;
  }

  // Evaluate the complete Born amplitude banks once before the convolution
  const double born_amp2 = EvaluateAmplitude();
  RequireAmplitudeRepresentation(metadata, state.lts.proton_good_walker);

  if (state.lts.hamp.size() == 0) {
    if (std::fpclassify(born_amp2) == FP_ZERO) { return 0.0; }
    throw AmplitudeFailure(
        "MProcess::ScreenedAmplitudeSquared: positive amplitude has "
        "no screening components");
  }

  // Save the accepted Born momenta for exact restoration after shifted evaluations
  state.lts.pfinal_orig = state.lts.pfinal;

  // Keep the Born transverse momenta fixed as the integration origin
  const std::array<double, 2> p1T = {state.lts.pfinal[1].Px(), state.lts.pfinal[1].Py()};
  const std::array<double, 2> p2T = {state.lts.pfinal[2].Px(), state.lts.pfinal[2].Py()};

  // Nuclear UPC screening uses the optical or sampled Glauber kernel
  const auto &screening_upc = state.lts.upc_event != nullptr ? state.lts.upc_event : state.lts.upc_model;
  if (screening_upc != nullptr) {
    if (screening_upc->HadronicConvolution() && screening_upc->Param().survival == nuclear::SurvivalType::MCGGCF &&
        !screening_upc->HasSamples()) {
      throw AmplitudeFailure(
          "MProcess::ScreenedAmplitudeSquared: nuclear fluctuation was not "
          "sampled for the event");
    }
    if (metadata.amplitude_type != ScreeningAmplitudeType::Physical) {
      throw AmplitudeFailure(
          "MProcess::ScreenedAmplitudeSquared: nuclear UPC screening requires a "
          "physical helicity amplitude");
    }
    return ScreenUPCAmplitude(p1T, p2T, *screening_upc);
  }

  // Proton-proton screening uses the coupled-channel eikonal kernel
  const auto &loop_const = eikonal.GetLoopConst(state.lts.s);

  if (metadata.amplitude_type == ScreeningAmplitudeType::GoodWalker) {
    ProtonGoodWalkerAmplitude born = std::move(*state.lts.proton_good_walker);
    state.lts.proton_good_walker.reset();
    return ScreenPPGoodWalker(std::move(born), metadata, p1T, p2T, loop_const);
  }

  return ScreenPPAmplitude(p1T, p2T, loop_const);
}


// Screen physical helicity amplitudes with the pp eikonal
//
double MProcess::ScreenPPAmplitude(const std::array<double, 2> &p1T, const std::array<double, 2> &p2T,
                                   const MEikonal::LoopConst &loop_const) {
  const std::vector<std::complex<double>> hamp_0             = state.lts.hamp;
  const ScreeningMetadata                 screening_metadata = state.lts.hamp.metadata;
  const ScreeningLayout                   screening_layout   = state.lts.hamp.layout;
  const auto                              hard_color_flows_0 = state.lts.hard_color_flows;

  std::vector<std::vector<std::complex<double>>> bank_source = {hamp_0};
  for (const auto &flow : hard_color_flows_0) { bank_source.push_back(flow.amplitudes); }
  eikonal::MProtonScreen screen(bank_source, screening_metadata, loop_const, EventGoodWalkerClass(state.lts));

  const int    born_id1       = state.lts.id1;
  const int    born_id2       = state.lts.id2;
  const double born_pdf_xf1   = state.lts.pdf_xf1;
  const double born_pdf_xf2   = state.lts.pdf_xf2;
  const double born_muF       = state.lts.muF;
  const double born_muR       = state.lts.muR;
  const double born_scalup    = state.lts.scalup;
  const bool   born_exact_epa = state.lts.exact_forward_photon_kinematics;

  // Recontract the prepared EPA tensor with node-local sources, or rebuild the full amplitude
  const auto evaluate_node = [&](const double kx, const double ky) {
    const std::array<double, 2> p1T_new = {p1T[0] - kx, p1T[1] - ky};
    const std::array<double, 2> p2T_new = {p2T[0] + kx, p2T[1] + ky};
    if (!EvaluateScreeningNode(p1T_new, p2T_new)) { return false; }
    if (state.lts.hamp.size() != hamp_0.size() || !(state.lts.hamp.layout == screening_layout) ||
        !SameHardColorLayout(hard_color_flows_0, state.lts.hard_color_flows)) {
      throw AmplitudeFailure("MProcess::ScreenPPAmplitude: shifted amplitude layout changed");
    }
    return true;
  };

  // View node-local physical and projected rows without allocating in the loop
  std::vector<std::span<const std::complex<double>>> node_source(bank_source.size());
  const auto                                         view_node = [&]() {
    node_source.front() = state.lts.hamp;
    for (const auto &flow : indices(state.lts.hard_color_flows)) {
      node_source[flow + 1] = state.lts.hard_color_flows[flow].amplitudes;
    }
  };

  ScreeningLoopGuard loop_guard(*this);
  const std::size_t  n_phi = loop_const.physical_screening_weight.size_col();
  for (const auto &i : indices(loop_const.kt2)) {
    for (std::size_t j = 0; j < n_phi; ++j) {
      if (!evaluate_node(loop_const.kt_x[i][j], loop_const.kt_y[i][j])) { continue; }
      view_node();
      screen.Add(i, j, node_source);
    }
  }

  // Restore the accepted Born point without solving its longitudinal momenta
  // again
  if (!RestoreScreeningBornState()) {
    throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: failed to restore Born kinematics");
  }

  // Keep event-record flux information at the accepted Born point
  state.lts.id1                             = born_id1;
  state.lts.id2                             = born_id2;
  state.lts.pdf_xf1                         = born_pdf_xf1;
  state.lts.pdf_xf2                         = born_pdf_xf2;
  state.lts.muF                             = born_muF;
  state.lts.muR                             = born_muR;
  state.lts.scalup                          = born_scalup;
  state.lts.exact_forward_photon_kinematics = born_exact_epa;

  // ------------------------------------------------------------
  // Final amplitude (squared)

  // Publish the physical bank and keep color metadata outside the convolution
  auto screened  = screen.Result();
  state.lts.hamp = std::move(screened.front());

  state.lts.hard_color_flows = hard_color_flows_0;
  for (auto &flow : state.lts.hard_color_flows) { flow.screened_weight.reset(); }
  for (const auto &c : indices(hard_color_flows_0)) {
    state.lts.hard_color_flows[c].amplitudes = std::move(screened[c + 1]);
  }

  const bool             dense_proton_spin = state.lts.hamp.size() != hamp_0.size();
  const AmplitudeWeights weight = ProcPtr.NormalizeScreenedAmplitude(state.lts.hamp, dense_proton_spin);
  PostScreeningAmplitude(weight.helicity_weight);
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    throw AmplitudeFailure("MProcess::ScreenedAmplitudeSquared: amplitude finalization failed");
  }
  loop_guard.Accept();
  return weight.total_weight;
}


// Sample nuclear fluctuations once for the active UPC event evaluation
//
void MProcess::PrepareUPCEvent(const bool include_screening) {
  if (state.lts.upc_model == nullptr || state.lts.upc_event != nullptr) { return; }
  const auto &param          = state.lts.upc_model->Param();
  const bool  convolve       = include_screening && state.lts.upc_model->HadronicConvolution();
  bool        current_sample = false;
  for (const auto &i : indices(param.emission)) {
    const int leg = static_cast<int>(i + 1);
    current_sample |=
        (state.lts.upc_model->Photon(leg) != nullptr && param.emission[i] != nuclear::CoherenceType::Coherent) ||
        (state.lts.upc_model->Photo(leg) != nullptr && param.target[i] != nuclear::CoherenceType::Coherent);
  }
  const bool config_survival =
      convolve && state.lts.upc_model->Glauber() != nullptr && param.survival == nuclear::SurvivalType::MCGGCF;
  if (param.config.count == 0 || (!config_survival && !current_sample)) {
    state.lts.upc_event = state.lts.upc_excitation
                              ? state.lts.upc_model->WithExcitation(state.lts.upc_excitation, include_screening)
                              : state.lts.upc_model;
    return;
  }
  try {
    // Isolate variable-cost nuclear sampling from the phase-space RNG stream
    MRandom nuclear_random;
    nuclear_random.SetSeed(static_cast<std::uint32_t>(state.random.rng()));
    const std::size_t count = config_survival ? 0 : param.current_count;
    state.lts.upc_event     = state.lts.upc_model->Sample(nuclear_random, convolve, count);
    if (state.lts.upc_excitation) {
      state.lts.upc_event = state.lts.upc_event->WithExcitation(state.lts.upc_excitation, include_screening);
    }
  } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("MProcess::PrepareUPCEvent: ") + error.what());
  }
}

// Normalize a resolved nuclear projection and finalize event weights
double MProcess::FinalizeUPCAmplitude(const nuclear::ScreenResult& result, const nuclear::ScreenLayout& layout,
                                      const ScreeningMetadata& metadata) {
  if (state.lts.hard_color_flows.size() != result.color_norm.size()) {
    throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: screened hard color flow count changed");
  }
  for (const auto& i : indices(result.color_norm)) {
    if (!std::isfinite(result.color_norm[i]) || result.color_norm[i] < 0.0) {
      throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: invalid hard color flow probability");
    }
    state.lts.hard_color_flows[i].screened_weight = result.color_norm[i];
  }
  state.lts.hamp.assign(result.amplitude.begin(), result.amplitude.end());
  state.lts.hamp.Configure(metadata);
  if (!gra::AllFinite(state.lts.hamp) || !gra::AllFinite(result.photo)) {
    throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: non-finite amplitude");
  }
  const AmplitudeWeights weight = ProcPtr.NormalizeAmplitude(result.helicity_norm, metadata);
  SetUPCFinalWeights(nuclear::FinalSectorWeights(layout, weight.helicity_weight));
  state.upc_photo = result.photo;
  PostScreeningAmplitude(weight.helicity_weight);
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: amplitude finalization failed");
  }
  return weight.total_weight;
}

// Screen one heavy-ion UPC amplitude through the nuclear screening
//
double MProcess::ScreenUPCAmplitude(const std::array<double, 2> &p1T, const std::array<double, 2> &p2T,
                                    const nuclear::MUPC &upc) {
  try {
    const ScreeningMetadata     metadata        = state.lts.hamp.metadata;
    const ScreeningLayout       event_layout    = state.lts.hamp.layout;
    const nuclear::ScreenLayout layout          = MakeUPCScreenLayout(event_layout, upc);
    const auto                  born_hard_flows = state.lts.hard_color_flows;

    const nuclear::ScreenPoint born = TakeUPCScreenPoint(state.lts, HardColorAmplitudes(born_hard_flows));
    nuclear::MUPCScreen    screen(upc, layout, born);
    const int                  born_id1       = state.lts.id1;
    const int                  born_id2       = state.lts.id2;
    const double               born_pdf_xf1   = state.lts.pdf_xf1;
    const double               born_pdf_xf2   = state.lts.pdf_xf2;
    const double               born_muF       = state.lts.muF;
    const double               born_muR       = state.lts.muR;
    const double               born_scalup    = state.lts.scalup;
    const bool                 born_exact_epa = state.lts.exact_forward_photon_kinematics;

    ScreeningLoopGuard loop_guard(*this);
    for (const auto& node : upc.Nodes(born.transfer)) {
      const std::array<double, 2> p1T_new = {p1T[0] - node.kx, p1T[1] - node.ky};
      const std::array<double, 2> p2T_new = {p2T[0] + node.kx, p2T[1] + node.ky};
      if (!EvaluateScreeningNode(p1T_new, p2T_new)) { continue; }
      if (state.lts.hamp.empty() || !(state.lts.hamp.layout == event_layout)) {
        throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: shifted amplitude layout changed");
      }
      if (state.lts.hard_color_flows.size() != born_hard_flows.size()) {
        throw AmplitudeFailure("MProcess::ScreenUPCAmplitude: shifted hard color flow layout changed");
      }
      screen.Add(node, TakeUPCScreenPoint(state.lts, HardColorAmplitudes(state.lts.hard_color_flows)));
    }

    if (!RestoreScreeningBornState()) {
      throw AmplitudeFailure(
          "MProcess::ScreenUPCAmplitude: failed to restore Born "
          "kinematics");
    }
    state.lts.hard_color_flows = born_hard_flows;

    // Keep event-record flux information at the accepted Born point
    state.lts.id1                             = born_id1;
    state.lts.id2                             = born_id2;
    state.lts.pdf_xf1                         = born_pdf_xf1;
    state.lts.pdf_xf2                         = born_pdf_xf2;
    state.lts.muF                             = born_muF;
    state.lts.muR                             = born_muR;
    state.lts.scalup                          = born_scalup;
    state.lts.exact_forward_photon_kinematics = born_exact_epa;

    const double total = FinalizeUPCAmplitude(screen.Result(), layout, metadata);
    loop_guard.Accept();
    return total;
  } catch (const AmplitudeFailure &) { throw; } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("MProcess::ScreenUPCAmplitude: ") + error.what());
  }
}

}  // namespace gra
