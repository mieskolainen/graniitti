// Hard diffraction with factorized or adaptive central phase space
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <functional>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MHardDiffraction.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/Attribute.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

namespace gra {

using gra::math::CheckEMC;
using gra::math::msqrt;
using gra::math::pow2;
using gra::PDG::GeV2barn;

// Construct the process list
MHardDiffraction::MHardDiffraction() { Initialize(); }

// Construct one selected hard-diffraction process
MHardDiffraction::MHardDiffraction(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  Initialize(process.find("<C>") != std::string::npos ? "C" : "F");
  InitHistograms();
  SetProcess(process, syntax, std::move(tune));

  M4Vec zerovec(0, 0, 0, 0);
  state.lts.pfinal.assign(11, zerovec);

  std::cout << "MHardDiffraction:: [Constructor done]" << std::endl;
}

// Destroy one selected hard-diffraction process
MHardDiffraction::~MHardDiffraction() {}

// Initialize supported hard-diffraction process names
void MHardDiffraction::Initialize(const std::string &mode) {
  const std::vector<std::string> supported = {"IPp", "IPIP"};
  state.phase_space_class                  = mode;
  ProcPtr                                  = MSubProc(supported, state.phase_space_class);
  const MSubProc other(supported, mode == "F" ? "C" : "F");
  ProcPtr.ProcessRegistry.insert(other.ProcessRegistry.begin(), other.ProcessRegistry.end());
  state.lts.central_phase_space_mode       = CentralPhaseSpaceMode::HardDiffraction;
}

// Initialize process-specific phase-space dimensionality
void MHardDiffraction::FinalizeProcessConfiguration() {
  if (state.lts.decaytree.size() > 8) {
    throw std::invalid_argument(
        "MHardDiffraction::FinalizeProcessConfiguration: direct "
        "central multiplicity cannot exceed 8");
  }

  // Hard diffraction currently has no consistent excited-forward event model
  if (state.excitation != 0) {
    throw std::invalid_argument(
        "MHardDiffraction::FinalizeProcessConfiguration: NSTARS is not "
        "supported for IPp "
        "or IPIP hard diffraction");
  }

  if (GetSoftModel() == nullptr) {
    throw std::invalid_argument(
        "MHardDiffraction::FinalizeProcessConfiguration: missing SOFT model "
        "snapshot");
  }
  if (state.lts.model_cache == nullptr) {
    throw std::invalid_argument(
        "MHardDiffraction::FinalizeProcessConfiguration: missing run owned "
        "model caches");
  }
  hard_pomeron_pdf = state.lts.model_cache->hard_pomeron.GetHardPomeronPDF(GetSoftModel());

  if (ProcPtr.ISTATE == "IPp") {
    if (!(ProtonRemnantMass() > PDG::mp)) {
      throw std::invalid_argument("MHardDiffraction::FinalizeProcessConfiguration: remnant_mass must exceed the proton mass");
    }
    for (const int flavour : HardPartonFlavours()) {
      const auto ids = ProtonRemnantIDs(flavour, PDG::PDG_p);
      for (const int id : ids) {
        if (!IsDiquark(id) && state.lts.PDG.PDG_table.count(id) == 0) {
          throw std::invalid_argument("MHardDiffraction::FinalizeProcessConfiguration: missing remnant PDG data");
        }
      }
      if (!(ProtonRemnantMass() > HardParticle(ids[0]).mass + HardParticle(ids[1]).mass)) {
        throw std::invalid_argument("MHardDiffraction::FinalizeProcessConfiguration: remnant_mass is below a constituent threshold");
      }
    }
  }

  SetTechnicalBoundaries(state.gcuts, state.excitation);
  if (ProcPtr.ISTATE == "IPp" && state.lts.GlobalPdfPtr != nullptr) {
    hard_pomeron_pdf->ValidateAlphaS(*state.lts.GlobalPdfPtr,
                                    hard_pomeron_pdf->FactorizationQ2(pow2(state.gcuts.M_min)),
                                    hard_pomeron_pdf->FactorizationQ2(pow2(state.gcuts.M_max)));
  }


  const bool double_diff = (ProcPtr.ISTATE == "IPIP");
  ProcPtr.LIPSDIM        = double_diff ? 8 : 5;
  if (state.phase_space_class == "C" && state.lts.decaytree.size() > 1) {
    if (!state.lts.PS_active) {
      throw std::invalid_argument("Hard diffraction <C> requires the central phase-space weight, use <F> for normalized decays");
    }
    central_coordinates.resize(3 * state.lts.decaytree.size() - 4);
    ProcPtr.LIPSDIM += central_coordinates.size();
  }
}

// Update loop kinematics for the eikonal screening integration
bool MHardDiffraction::LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) {
  if (!kinematics::RebuildScreeningKinematics(state.lts, p1p, p2p, true)) { return false; }
  if (!kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true,
                                     kinematics::TransferVirtualityPolicy::HardDiffractive)) {
    return false;
  }
  SetTaggedForwardXi();

  return BuildHardRemnantBranches(false);
}

// Recompute hard diffraction kinematics from the restored Born four-momenta
bool MHardDiffraction::RefreshBornKinematics() {
  if (!kinematics::SetLorentzScalars(state, state.lts.decaytree.size() + 2, true,
                                     kinematics::TransferVirtualityPolicy::HardDiffractive)) {
    return false;
  }
  SetTaggedForwardXi();
  return BuildHardRemnantBranches(false);
}

// Compute one event weight
double MHardDiffraction::ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) {
  hard_parton_channels.clear();
  hard_channel_component_counts.clear();
  hard_component_channel_indices.clear();
  hard_color_flows.clear();
  hard_selected_channel_index = 0;
  double W                    = 0.0;

  if (parton::PrepareFinalStateProposal(state.lts, state.random)) { CalculateSymmetryFactor(); }

  PreparePhaseSpacePoint(HardRandomKin(randvec), aux);

  if (aux.Valid()) {
    const double MatESQ = GetAmp2(aux.include_screening, aux);

    double exact = 0.0;
    DecayWidthPS(exact);
    if (!aux.adaptation_mode) {
      state.lts.DW_sum_exact.AddLogWeight(gra::kinematics::MCW(exact), aux.log_inverse_density);
      state.lts.DW_sum.AddLogWeight(state.lts.DW, aux.log_inverse_density);
    }

    double C_space = 1.0;
    if (state.lts.decaytree.size() != 0 && state.lts.PS_active) {
      C_space = state.lts.DW.Integral();
      C_space *= CascadePS();
    }

    W = DissociationCrossSectionFactor() * DecaySymmetryCompensationFactor() * C_space * (1.0 / state.symmetry_factor) *
        MatESQ * HardIntegralVolume() * HardPhaseSpaceWeight() * parton::ProposalWeight(state.lts) * GeV2barn /
        kinematics::MollerFlux(state.lts.pbeam1, state.lts.pbeam2);
  }

  return W;
}

// Apply common fiducial cuts
bool MHardDiffraction::FiducialCuts() const { return CommonCuts(); }

// Apply hard-diffraction generation cuts to the sampled hard system
bool MHardDiffraction::HardGenerationCuts() const {
  const double M = state.lts.pfinal[0].M();
  const double Y = state.lts.pfinal[0].Rap();
  if (!std::isfinite(M) || !std::isfinite(Y)) { return false; }
  if (M < state.gcuts.M_min || M > state.gcuts.M_max) { return false; }
  if (Y < state.gcuts.Y_min || Y > state.gcuts.Y_max) { return false; }

  if (state.lts.hard_diff1 && !(state.lts.diff_xi1 >= state.gcuts.XI_min && state.lts.diff_xi1 <= state.gcuts.XI_max)) {
    return false;
  }
  if (state.lts.hard_diff2 && !(state.lts.diff_xi2 >= state.gcuts.XI_min && state.lts.diff_xi2 <= state.gcuts.XI_max)) {
    return false;
  }

  return true;
}

// Save one event to HepMC
bool MHardDiffraction::BuildEventRecord(HepMC3::GenEvent &evt) {
  if ((!state.lts.hard_diff1 && (state.lts.id1 == 0 || ProtonRemnantIDs(state.lts.id1, state.lts.beam1.pdg)[0] == 0)) ||
      (!state.lts.hard_diff2 && (state.lts.id2 == 0 || ProtonRemnantIDs(state.lts.id2, state.lts.beam2.pdg)[0] == 0))) {
    return false;
  }
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, leaves); }
  for (const auto *leaf : leaves) {
    if (std::abs(leaf->p.pdg) == PDG::PDG_hard_jet) { return false; }
  }
  if (!MProcess::BuildEventRecord(evt)) { return false; }
  const int sides = (state.lts.hard_diff1 ? 1 : 0) | (state.lts.hard_diff2 ? 2 : 0);
  evt.add_attribute("graniitti_hard_diffraction", std::make_shared<HepMC3::IntAttribute>(sides));
  // Cross the net remnant color into each incoming hard parton for the optional ISR view
  const std::array<MDecayBranch *, 2> forward = {&state.lts.decayforward1, &state.lts.decayforward2};
  for (const auto &side : indices(forward)) {
    leaves.clear();
    CollectStableLeaves(*forward[side], leaves);
    std::map<int, int> balance;
    for (const auto *leaf : leaves) {
      if (leaf->p.color_flow.flow1 != 0) { ++balance[leaf->p.color_flow.flow1]; }
      if (leaf->p.color_flow.flow2 != 0) { --balance[leaf->p.color_flow.flow2]; }
    }
    for (const auto &[tag, count] : balance) {
      if (count == 0) { continue; }
      if (std::abs(count) != 1) { return false; }
      const std::string name = "graniitti_hard" + std::to_string(side + 1) + (count < 0 ? "_flow1" : "_flow2");
      if (evt.attribute<HepMC3::IntAttribute>(name)) { return false; }
      evt.add_attribute(name, std::make_shared<HepMC3::IntAttribute>(tag));
    }
  }
  return true;
}

// Build the incoming hard-parton basis for the current phase-space point
std::vector<MHardDiffraction::HardPartonChannel> MHardDiffraction::BuildHardPartonChannels(double Q2) {
  const auto flavours = HardPartonFlavours();

  std::vector<std::pair<int, double>> density1;
  std::vector<std::pair<int, double>> density2;
  density1.reserve(flavours.size());
  density2.reserve(flavours.size());

  for (const auto id1 : flavours) {
    const double f1 = HardPartonDensityForSide(id1, true, Q2);
    if (std::isfinite(f1) && f1 > 0.0) { density1.emplace_back(id1, f1); }
  }

  for (const auto id2 : flavours) {
    const double f2 = HardPartonDensityForSide(id2, false, Q2);
    if (std::isfinite(f2) && f2 > 0.0) { density2.emplace_back(id2, f2); }
  }

  std::vector<HardPartonChannel> channels;
  channels.reserve(density1.size() * density2.size());
  for (const auto &leg1 : density1) {
    for (const auto &leg2 : density2) {
      channels.push_back({leg1.first, leg2.first, leg1.second, leg2.second, state.lts.diff_xhard1 * leg1.second,
                          state.lts.diff_xhard2 * leg2.second});
    }
  }
  return channels;
}

// Recalculate channel densities without changing the cached channel ordering
std::vector<MHardDiffraction::HardPartonChannel> MHardDiffraction::RefreshHardPartonChannelDensities(
    const std::vector<HardPartonChannel> &channels, double Q2) {
  std::vector<HardPartonChannel> refreshed;
  refreshed.reserve(channels.size());

  for (const auto &channel : channels) {
    const double f1      = HardPartonDensityForSide(channel.id1, true, Q2);
    const double f2      = HardPartonDensityForSide(channel.id2, false, Q2);
    const double safe_f1 = (std::isfinite(f1) && f1 > 0.0) ? f1 : 0.0;
    const double safe_f2 = (std::isfinite(f2) && f2 > 0.0) ? f2 : 0.0;
    refreshed.push_back(
        {channel.id1, channel.id2, safe_f1, safe_f2, state.lts.diff_xhard1 * safe_f1, state.lts.diff_xhard2 * safe_f2});
  }

  return refreshed;
}

// Compute the current diffractive proton momentum transfer for one beam side
double MHardDiffraction::CurrentDiffractiveT(bool first_side) const {
  const bool   hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  const double sampled_t = first_side ? state.lts.diff_t1 : state.lts.diff_t2;
  if (!hard_diff || !state.lts.screening.active) { return sampled_t; }

  const MDecayBranch &branch = first_side ? state.lts.decayforward1 : state.lts.decayforward2;
  if (branch.legs.empty()) { return sampled_t; }

  const M4Vec  beam     = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const M4Vec  transfer = beam - branch.legs.front().p4;
  const double t        = transfer.M2();
  return std::isfinite(t) ? t : sampled_t;
}

// Compute the PDF or DPDF density for one incoming beam side
double MHardDiffraction::HardPartonDensityForSide(int pid, bool first_side, double Q2) {
  const bool hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  if (!hard_diff) {
    const double xhard = first_side ? state.lts.diff_xhard1 : state.lts.diff_xhard2;
    const int beam_pdg = first_side ? state.lts.beam1.pdg : state.lts.beam2.pdg;
    return ProtonPartonDensity(beam_pdg < 0 && pid != PDG::PDG_gluon ? -pid : pid, xhard, Q2);
  }

  const double xi   = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
  const double beta = first_side ? state.lts.diff_beta1 : state.lts.diff_beta2;
  const double t    = CurrentDiffractiveT(first_side);
  try {
    return hard_pomeron_pdf->DiffractiveDensity(pid, xi, beta, t, Q2);
  } catch (...) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return 0.0;
  }
}

// Evaluate the PDF-weighted hard subprocess amplitude
double MHardDiffraction::EvaluateBareAmplitude() {
  const double Q2         = hard_pomeron_pdf->FactorizationQ2(state.lts.pfinal[0].M2());
  const double hard_mass2 = state.lts.s * state.lts.diff_xhard1 * state.lts.diff_xhard2;
  // Convert the outer proton flux to the sampled massless parton flux
  // F_beam/F_parton = F_beam/(2 shat)
  const double flux_ratio = std::isfinite(hard_mass2) && hard_mass2 > 0.0
                                ? kinematics::MollerFlux(state.lts.pbeam1, state.lts.pbeam2) / (2.0 * hard_mass2)
                                : 0.0;
  if (!(std::isfinite(Q2) && Q2 > 0.0) || !(std::isfinite(flux_ratio) && flux_ratio > 0.0)) {
    evaluation_status = mg5helas::EvaluationStatus::KinematicsFailure;
    state.lts.hamp.clear();
    return 0.0;
  }

  state.lts.alphaQCD = AlphaQCD(Q2);
  if (!mg5helas::EvaluationSucceeded(evaluation_status) || !std::isfinite(state.lts.alphaQCD) ||
      state.lts.alphaQCD <= 0.0) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    state.lts.hamp.clear();
    return 0.0;
  }
  state.lts.muF    = msqrt(Q2);
  state.lts.muR    = state.lts.muF;
  state.lts.scalup = state.lts.muF;
  ProcPtr.PrepareBareAmplitude(state.lts);
  if (!mg5helas::EvaluationSucceeded(ProcPtr.EvaluationStatus())) {
    evaluation_status = ProcPtr.EvaluationStatus();
    state.lts.hamp.clear();
    return 0.0;
  }

  std::vector<HardPartonChannel> channels;
  if (state.lts.screening.active && !hard_parton_channels.empty()) {
    channels = RefreshHardPartonChannelDensities(hard_parton_channels, Q2);
  } else {
    channels             = BuildHardPartonChannels(Q2);
    hard_parton_channels = channels;
  }
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    state.lts.hamp.clear();
    return 0.0;
  }

  std::vector<std::complex<double>> channel_hamp;
  std::vector<std::size_t>          channel_component_counts(channels.size(), 0);
  std::vector<std::size_t>          component_channel_indices;
  // Keep raw MG5 flow amplitudes aligned with their incoming parton channels
  std::vector<mg5helas::HardColorFlow> color_flows;

  double weighted_amp2 = 0.0;
  for (const auto &channel_index : gra::aux::indices(channels)) {
    const auto &channel = channels[channel_index];
    state.lts.id1       = channel.id1;
    state.lts.id2       = channel.id2;
    const double amp2   = ProcPtr.GetPreparedBareAmplitude2(state.lts);
    if (!mg5helas::EvaluationSucceeded(ProcPtr.EvaluationStatus())) {
      evaluation_status = ProcPtr.EvaluationStatus();
      state.lts.hamp.clear();
      state.lts.hard_color_flows.clear();
      return 0.0;
    }
    if (state.lts.screening.active && state.lts.hamp.empty() &&
        channel_index < hard_channel_component_counts.size()) {
      state.lts.hamp.assign(hard_channel_component_counts[channel_index], 0.0);
    }
    channel_component_counts[channel_index] = state.lts.hamp.size();
    const double density2                   = channel.f1 * channel.f2 * flux_ratio;
    const double density_scale              = (std::isfinite(density2) && density2 > 0.0) ? msqrt(density2) : 0.0;

    if (amp2 > 0.0 && state.lts.hamp.empty()) {
      evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
      state.lts.hard_color_flows.clear();
      return 0.0;
    }
    for (const auto amplitude : state.lts.hamp) {
      const std::complex<double> value      = density_scale * amplitude;
      const double               component2 = std::norm(value);
      if (std::isfinite(component2)) {
        channel_hamp.push_back(value);
        weighted_amp2 += component2;
      } else {
        evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
        state.lts.hamp.clear();
        state.lts.hard_color_flows.clear();
        return 0.0;
      }
      component_channel_indices.push_back(channel_index);
    }

    // Weight each exact generated flow by the same PDF density as its helicity
    // amplitude
    const auto &generated_flows = state.lts.hard_color_flows;
    if (!generated_flows.empty()) {
      if (!BuildHardRemnantBranches(false)) {
        evaluation_status = mg5helas::EvaluationStatus::KinematicsFailure;
        state.lts.hamp.clear();
        state.lts.hard_color_flows.clear();
        return 0.0;
      }
      for (const auto &flow : generated_flows) {
        std::vector<MColorFlow> candidate;
        if (!BuildGeneratedHardColorFlowCandidate(flow.external, candidate)) {
          evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
          state.lts.hamp.clear();
          state.lts.hard_color_flows.clear();
          return 0.0;
        }
        std::vector<std::complex<double>> weighted_flow;
        weighted_flow.reserve(flow.amplitudes.size());
        for (const auto amplitude : flow.amplitudes) {
          const std::complex<double> weighted = density_scale * amplitude;
          if (!std::isfinite(weighted.real()) || !std::isfinite(weighted.imag())) {
            evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
            state.lts.hamp.clear();
            state.lts.hard_color_flows.clear();
            return 0.0;
          }
          weighted_flow.push_back(weighted);
        }
        mg5helas::ExternalColorFlow assignment;
        assignment.reserve(candidate.size());
        for (const auto &leg : candidate) { assignment.push_back({leg.flow1, leg.flow2}); }
        color_flows.push_back({std::move(weighted_flow), std::move(assignment), std::nullopt, channel_index});
      }
    }
  }
  state.lts.hamp = std::move(channel_hamp);

  if (state.lts.screening.active) {
    if (channel_component_counts != hard_channel_component_counts ||
        component_channel_indices != hard_component_channel_indices ||
        !HardColorFlowsMatch(color_flows, hard_color_flows)) {
      evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
      state.lts.hamp.clear();
      state.lts.hard_color_flows.clear();
      return 0.0;
    }
  } else {
    hard_channel_component_counts  = std::move(channel_component_counts);
    hard_component_channel_indices = std::move(component_channel_indices);
    hard_color_flows               = color_flows;
  }

  state.lts.hard_color_flows = std::move(color_flows);

  if (channels.empty()) {
    state.lts.id1     = 0;
    state.lts.id2     = 0;
    state.lts.pdf_xf1 = 0.0;
    state.lts.pdf_xf2 = 0.0;
  }

  return weighted_amp2;
}

// Sample one incoming hard-parton channel from final component weights
void MHardDiffraction::PostScreeningAmplitude(const std::vector<double> &component_amp_squared) {
  for (const auto value : component_amp_squared) {
    if (!std::isfinite(value) || value < 0.0) {
      evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
      state.lts.hamp.clear();
      return;
    }
  }
  if (!SelectHardPartonChannel(component_amp_squared)) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
  }
}

// Finalize the selected hard channel and both remnant branches atomically
bool MHardDiffraction::FinalizeAmplitudeEventState() { return BuildHardRemnantBranches(true); }

// Publish tagged leading proton losses without confusing them with beta times
// xi
void MHardDiffraction::SetTaggedForwardXi() {
  state.lts.has_xi1 = state.lts.hard_diff1;
  state.lts.has_xi2 = state.lts.hard_diff2;
  state.lts.xi1     = state.lts.hard_diff1 ? state.lts.diff_xi1 : 0.0;
  state.lts.xi2     = state.lts.hard_diff2 ? state.lts.diff_xi2 : 0.0;
}

// Map one unit coordinate logarithmically onto a positive interval
// x = xmin exp[u log(xmax/xmin)], dx/du = x log(xmax/xmin)
bool MHardDiffraction::MapLog(double unit, double min, double max, double &value, double &jacobian) {
  value    = 0.0;
  jacobian = 0.0;
  if (!std::isfinite(unit) || unit < 0.0 || unit > 1.0 || !std::isfinite(min) || !std::isfinite(max) || !(min > 0.0) ||
      !(min < max)) {
    return false;
  }

  const double log_ratio = std::log(max / min);
  value                  = min * std::exp(unit * log_ratio);
  jacobian               = value * log_ratio;
  return std::isfinite(value) && std::isfinite(jacobian) && jacobian > 0.0;
}

// Map one unit coordinate with density proportional to exp(Bx)
// p(x) = B exp(Bx)/[exp(B xmax)-exp(B xmin)], dx/du = 1/p(x)
bool MHardDiffraction::MapExp(double unit, double min, double max, double slope, double &value, double &jacobian) {
  value    = 0.0;
  jacobian = 0.0;
  if (!std::isfinite(unit) || unit < 0.0 || unit > 1.0 || !std::isfinite(min) || !std::isfinite(max) ||
      !std::isfinite(slope) || !(min < max)) {
    return false;
  }

  const double span = max - min;
  if (std::abs(slope * span) < 1.0e-8) {
    value    = min + unit * span;
    jacobian = span;
    return true;
  }

  if (slope > 0.0) {
    const double ratio = std::exp(-slope * span);
    const double mix   = ratio + unit * (1.0 - ratio);
    value              = max + std::log(mix) / slope;
    jacobian           = (1.0 - ratio) / (slope * mix);
  } else {
    const double ratio = std::exp(slope * span);
    const double mix   = 1.0 + unit * (ratio - 1.0);
    value              = min + std::log(mix) / slope;
    jacobian           = (ratio - 1.0) / (slope * mix);
  }
  return std::isfinite(value) && std::isfinite(jacobian) && jacobian > 0.0;
}

// Map collinear hard fractions through mass squared and rapidity
// x1 = sqrt(M2/s)e^Y, x2 = sqrt(M2/s)e^-Y, |J| = 1/s
bool MHardDiffraction::MapHardLong(double mass_unit, double rapidity_unit, double s, double mass2_min, double mass2_max,
                                   const std::array<double, 2> &xhard1_range, const std::array<double, 2> &xhard2_range,
                                   double xi_product, HardLongPoint &point) {
  point = {};
  if (!std::isfinite(s) || !(s > 0.0) || !std::isfinite(xi_product) || !(xi_product > 0.0) ||
      !std::isfinite(rapidity_unit) || rapidity_unit < 0.0 || rapidity_unit > 1.0 || !std::isfinite(xhard1_range[0]) ||
      !std::isfinite(xhard1_range[1]) || !std::isfinite(xhard2_range[0]) || !std::isfinite(xhard2_range[1]) ||
      xhard1_range[0] < 0.0 || xhard2_range[0] < 0.0 || !(xhard1_range[0] < xhard1_range[1]) ||
      !(xhard2_range[0] < xhard2_range[1]) || xhard1_range[1] > 1.0 || xhard2_range[1] > 1.0) {
    return false;
  }

  const double reachable_min = s * xhard1_range[0] * xhard2_range[0];
  const double reachable_max = s * xhard1_range[1] * xhard2_range[1];
  const double lower         = std::max(mass2_min, reachable_min);
  const double upper         = std::min(mass2_max, reachable_max);
  if (!MapLog(mass_unit, lower, upper, point.mass2, point.mass_jacobian)) { return false; }

  const double ratio        = std::sqrt(point.mass2 / s);
  double       rapidity_min = std::log(ratio / xhard2_range[1]);
  double       rapidity_max = std::log(xhard1_range[1] / ratio);
  if (xhard1_range[0] > 0.0) { rapidity_min = std::max(rapidity_min, std::log(xhard1_range[0] / ratio)); }
  if (xhard2_range[0] > 0.0) { rapidity_max = std::min(rapidity_max, std::log(ratio / xhard2_range[0])); }
  if (!std::isfinite(rapidity_min) || !std::isfinite(rapidity_max) || !(rapidity_min < rapidity_max)) { return false; }

  const double rapidity_jacobian = rapidity_max - rapidity_min;
  point.rapidity                 = rapidity_min + rapidity_unit * rapidity_jacobian;
  point.xhard1                   = ratio * std::exp(point.rapidity);
  point.xhard2                   = ratio * std::exp(-point.rapidity);
  point.jacobian                 = point.mass_jacobian * rapidity_jacobian / (s * xi_product);
  return std::isfinite(point.xhard1) && std::isfinite(point.xhard2) && std::isfinite(point.jacobian) &&
         point.xhard1 > 0.0 && point.xhard2 > 0.0 && point.jacobian > 0.0;
}

// Sample one diffractive external leg and return its xi-t Jacobian
bool MHardDiffraction::SampleDiffractiveSide(bool first_side, double xi_unit, double t_unit, double phi_unit,
                                             double &jacobian) {
  jacobian = 0.0;
  if (!std::isfinite(phi_unit) || phi_unit < 0.0 || phi_unit > 1.0) { return false; }

  const auto  &xi_range    = hard_pomeron_pdf->XiRange();
  const double xi_min      = std::max(xi_range[0], state.gcuts.XI_min);
  const double xi_max      = std::min(xi_range[1], state.gcuts.XI_max);
  double       xi          = 0.0;
  double       xi_jacobian = 0.0;
  if (!MapLog(xi_unit, xi_min, xi_max, xi, xi_jacobian)) { return false; }

  const M4Vec &beam           = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const auto  &t_range        = hard_pomeron_pdf->TRange();
  const double physical_t_max = -xi * xi * std::max(beam.M2(), 0.0) / (1.0 - xi);
  const double t_max          = std::min(t_range[1], physical_t_max);
  double       t              = 0.0;
  double       t_jacobian     = 0.0;
  if (!MapExp(t_unit, t_range[0], t_max, hard_pomeron_pdf->FluxTSlope(xi), t, t_jacobian)) { return false; }

  const double phi = 2.0 * math::PI * phi_unit;
  jacobian         = xi_jacobian * t_jacobian;
  if (!std::isfinite(jacobian) || !(jacobian > 0.0)) { return false; }

  if (first_side) {
    state.lts.hard_diff1  = true;
    state.lts.diff_xi1    = xi;
    state.lts.diff_beta1  = 0.0;
    state.lts.diff_t1     = t;
    state.lts.diff_phi1   = phi;
    state.lts.diff_xhard1 = 0.0;
  } else {
    state.lts.hard_diff2  = true;
    state.lts.diff_xi2    = xi;
    state.lts.diff_beta2  = 0.0;
    state.lts.diff_t2     = t;
    state.lts.diff_phi2   = phi;
    state.lts.diff_xhard2 = 0.0;
  }
  return true;
}

// Compute the hard-fraction support for one current beam side
std::array<double, 2> MHardDiffraction::HardXRange(bool first_side) const {
  const bool hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  if (!hard_diff) {
    double     x_max      = 1.0 - 1.0e-12;
    const bool other_diff = first_side ? state.lts.hard_diff2 : state.lts.hard_diff1;
    // Keep the ordinary remnant causal after the exact mass recoil rescaling
    if (other_diff) {
      const double xi        = first_side ? state.lts.diff_xi2 : state.lts.diff_xi1;
      const double t         = first_side ? state.lts.diff_t2 : state.lts.diff_t1;
      const double lightcone = state.lts.pbeam1.LightconePos() * state.lts.pbeam2.LightconeNeg();
      if (!(xi > 0.0) || !(state.lts.s > 0.0) || !(lightcone > 0.0)) { return {1.0, 0.0}; }
      x_max = std::min(x_max, (lightcone + t / xi) / state.lts.s);
    }
    return {1.0e-12, x_max};
  }

  const double xi         = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
  const auto  &beta_range = hard_pomeron_pdf->BetaRange();
  return {xi * beta_range[0], xi * beta_range[1]};
}

// Store one sampled hard fraction as beta or ordinary proton x
bool MHardDiffraction::SetHardX(bool first_side, double xhard) {
  if (!std::isfinite(xhard) || !(xhard > 0.0) || !(xhard < 1.0)) { return false; }

  const bool hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  double     beta      = 0.0;
  if (hard_diff) {
    const double xi         = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
    const auto  &beta_range = hard_pomeron_pdf->BetaRange();
    beta                    = xhard / xi;
    const double tolerance  = 32.0 * std::numeric_limits<double>::epsilon() * beta_range[1];
    if (beta < beta_range[0] - tolerance || beta > beta_range[1] + tolerance) { return false; }
    beta  = std::clamp(beta, beta_range[0], beta_range[1]);
    xhard = xi * beta;
  }

  if (first_side) {
    state.lts.diff_xhard1 = xhard;
    state.lts.diff_beta1  = beta;
  } else {
    state.lts.diff_xhard2 = xhard;
    state.lts.diff_beta2  = beta;
  }
  return true;
}

// Sample both hard fractions through mass squared and rapidity
bool MHardDiffraction::SampleHardLong(double mass_unit, double rapidity_unit, double mass_sum, double mass_max,
                                      double &jacobian) {
  jacobian                  = 0.0;
  const double physical_min = std::max(state.gcuts.M_min, mass_sum);
  const double physical_max = std::min(state.gcuts.M_max, mass_max);
  if (!std::isfinite(physical_min) || !std::isfinite(physical_max) || !(physical_min < physical_max)) { return false; }

  const auto xhard1_range = HardXRange(true);
  const auto xhard2_range = HardXRange(false);
  double     xi_product   = 1.0;
  if (state.lts.hard_diff1) { xi_product *= state.lts.diff_xi1; }
  if (state.lts.hard_diff2) { xi_product *= state.lts.diff_xi2; }

  HardLongPoint point;
  if (!MapHardLong(mass_unit, rapidity_unit, state.lts.s, pow2(physical_min), pow2(physical_max), xhard1_range,
                   xhard2_range, xi_product, point)) {
    return false;
  }
  if (!SetHardX(true, point.xhard1) || !SetHardX(false, point.xhard2)) { return false; }
  // The rapidity range depends only on the total mass and cancels between decay histories
  state.lts.central_phase_space_mass_cut_min =
      std::max(state.gcuts.M_min, msqrt(state.lts.s * xhard1_range[0] * xhard2_range[0]));
  state.lts.central_phase_space_mass_max =
      std::min(physical_max, msqrt(state.lts.s * xhard1_range[1] * xhard2_range[1]));
  state.lts.central_phase_space_mass_margin = 0.0;
  state.lts.central_phase_space_generated_jacobian = point.mass_jacobian;
  jacobian = point.jacobian;
  return true;
}

// Build one random hard-diffraction phase-space point
bool MHardDiffraction::HardRandomKin(const std::vector<double> &randvec) {
  hard_integral_volume = 0.0;
  state.lts.id1 = state.lts.id2 = 0;
  remnant_angles = {state.random.U(0.0, 1.0), state.random.U(0.0, 1.0)};

  const bool double_ip = (ProcPtr.ISTATE == "IPIP");

  if ((!double_ip && randvec.size() < 5) || (double_ip && randvec.size() < 8)) { return false; }
  const std::size_t offset = double_ip ? 8 : 5;
  if (randvec.size() < offset + central_coordinates.size()) { return false; }
  std::copy_n(randvec.begin() + offset, central_coordinates.size(), central_coordinates.begin());

  double external_jacobian = 1.0;
  if (double_ip) {
    double first_jacobian  = 0.0;
    double second_jacobian = 0.0;
    if (!SampleDiffractiveSide(true, randvec[0], randvec[2], randvec[3], first_jacobian) ||
        !SampleDiffractiveSide(false, randvec[4], randvec[6], randvec[7], second_jacobian)) {
      return false;
    }
    external_jacobian = first_jacobian * second_jacobian;
  } else {
    const bool   first_side = randvec[0] < 0.5;
    const double xi_unit    = first_side ? 2.0 * randvec[0] : 2.0 * randvec[0] - 1.0;
    if (!SampleDiffractiveSide(first_side, xi_unit, randvec[2], randvec[3], external_jacobian)) { return false; }
  }

  double M_sum = 0.0;
  double M_max = state.gcuts.M_max;
  PrepareDecaySymmetryProposal();
  if (!PrepareCentralBranchMasses(M_sum, &M_max)) { return false; }

  double       longitudinal_jacobian = 0.0;
  const double rapidity_unit         = double_ip ? randvec[5] : randvec[4];
  if (!SampleHardLong(randvec[1], rapidity_unit, M_sum, M_max, longitudinal_jacobian)) { return false; }
  hard_integral_volume = external_jacobian * longitudinal_jacobian;

  return HardBuildKin(state.lts.diff_xhard1, state.lts.diff_xhard2);
}

// Build t-dependent hard-scattering kinematics
bool MHardDiffraction::HardBuildKin(double xhard1, double xhard2) {
  if (!(xhard1 > 0.0) || !(xhard2 > 0.0)) { return false; }

  const M4Vec beamsum = state.lts.pbeam1 + state.lts.pbeam2;

  M4Vec q1;
  M4Vec q2;
  if (!BuildHardPair(xhard1, xhard2, q1, q2)) { return false; }

  state.lts.pfinal[1] = state.lts.pbeam1 - q1;
  state.lts.pfinal[2] = state.lts.pbeam2 - q2;
  state.lts.pfinal[0] = q1 + q2;
  if (!HardGenerationCuts()) { return false; }

  if (!CheckEMC(beamsum - (state.lts.pfinal[1] + state.lts.pfinal[2] + state.lts.pfinal[0]))) { return false; }

  std::vector<double> masses;
  double              mass_sum = 0.0;
  for (const auto &i : gra::aux::indices(state.lts.decaytree)) {
    masses.push_back(state.lts.decaytree[i].m_offshell);
    mass_sum += state.lts.decaytree[i].m_offshell;
  }
  if (state.lts.pfinal[0].M2() <= pow2(mass_sum)) { return false; }

  std::vector<M4Vec>   products;
  const bool           UNWEIGHT = !state.lts.PS_active;
  gra::kinematics::MCW w;
  if (state.lts.decaytree.size() == 1) {
    if (!SetSingleCentralRootKinematics()) { return false; }
  } else if (!central_coordinates.empty()) {
    w = gra::kinematics::NBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products,
                                         central_coordinates);
  } else if (state.lts.decaytree.size() == 2) {
    w = gra::kinematics::TwoBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products,
                                           state.random);
  } else if (state.lts.decaytree.size() == 3) {
    w = gra::kinematics::ThreeBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products, UNWEIGHT,
                                             state.random);
  } else if (state.lts.decaytree.size() > 3) {
    w = gra::kinematics::NBodyPhaseSpace(state.lts.pfinal[0], state.lts.pfinal[0].M(), masses, products, UNWEIGHT,
                                         state.random);
  }

  if (state.lts.decaytree.size() != 1) {
    if (w.GetW() < 0) { return false; }
    state.lts.DW = w;

    const unsigned int offset = 3;
    for (const auto &i : gra::aux::indices(state.lts.decaytree)) {
      state.lts.decaytree[i].p4    = products[i];
      state.lts.pfinal[i + offset] = products[i];
    }
  }

  for (const auto &i : gra::aux::indices(state.lts.decaytree)) {
    if (!ConstructDecayKinematics(state.lts.decaytree[i])) { return false; }
  }
  if (!ApplyDecaySymmetryProposal()) { return false; }

  const unsigned int Nf = state.lts.decaytree.size() + 2;
  if (!kinematics::SetLorentzScalars(state, Nf, false, kinematics::TransferVirtualityPolicy::HardDiffractive)) {
    return false;
  }
  SetTaggedForwardXi();

  return BuildHardRemnantBranches(false);
}

// Build one hard parton pair with the exact sampled DPDF mass
bool MHardDiffraction::BuildHardPair(double xhard1, double xhard2, M4Vec &q1, M4Vec &q2) const {
  M4Vec hard1;
  M4Vec hard2;
  if (!(xhard1 > 0.0) || !(xhard2 > 0.0) || xhard1 >= 1.0 || xhard2 >= 1.0) { return false; }

  const double mass2 = state.lts.s * xhard1 * xhard2;
  if (!std::isfinite(mass2) || !(mass2 > 0.0)) { return false; }

  if (state.lts.hard_diff1 && state.lts.hard_diff2) {
    if (!BuildHardPartonForSide(true, hard1) || !BuildHardPartonForSide(false, hard2)) { return false; }
    if (!ScaleDoubleDiffractivePair(mass2, state.lts.diff_beta1, state.lts.diff_beta2, hard1, hard2)) { return false; }
  } else if (state.lts.hard_diff1 != state.lts.hard_diff2) {
    M4Vec &qdiff = state.lts.hard_diff1 ? hard1 : hard2;
    M4Vec &qcol  = state.lts.hard_diff1 ? hard2 : hard1;
    if (!BuildHardPartonForSide(state.lts.hard_diff1, qdiff) ||
        !BuildProtonHardParton(!state.lts.hard_diff1, qdiff, mass2, qcol)) { return false; }
  } else { return false; }

  const double closure = std::abs((hard1 + hard2).M2() - mass2);
  if (!AcceptHardParton(true, hard1) || !AcceptHardParton(false, hard2) || closure >= 1.0e-9 * std::max(1.0, mass2)) {
    return false;
  }
  q1 = hard1;
  q2 = hard2;
  return true;
}

// Scale two diffractive hard legs symmetrically to the exact DPDF mass
bool MHardDiffraction::ScaleDoubleDiffractivePair(double mass2, double beta1, double beta2, M4Vec &q1, M4Vec &q2) {
  if (!std::isfinite(mass2) || !(mass2 > 0.0) || !std::isfinite(beta1) || !std::isfinite(beta2) || !(beta1 > 0.0) ||
      !(beta2 > 0.0) || beta1 > 1.0 || beta2 > 1.0) {
    return false;
  }

  const double base_mass2 = (q1 + q2).M2();
  if (!std::isfinite(base_mass2) || !(base_mass2 > 0.0)) { return false; }
  double       scale     = std::sqrt(mass2 / base_mass2);
  const double tolerance = 128.0 * std::numeric_limits<double>::epsilon();
  if (!std::isfinite(scale) || scale < 1.0 - tolerance) { return false; }
  scale                      = std::max(scale, 1.0);
  const auto remnant_support = [&](double beta) {
    return scale * beta < 1.0 || (std::abs(scale - 1.0) <= tolerance && std::abs(beta - 1.0) <= tolerance);
  };
  if (!remnant_support(beta1) || !remnant_support(beta2)) { return false; }

  /*
   * beta and xhard remain the collinear DPDF integration variables. The common
   * Lorentz-scalar factor only assigns the finite-t recoil of the generated
   * event and leaves the hard-system rapidity and orientation unchanged.
   */
  M4Vec        scaled1 = q1 * scale;
  M4Vec        scaled2 = q2 * scale;
  const double closure = std::abs((scaled1 + scaled2).M2() - mass2);
  if (!(scaled1.E() > 0.0) || !(scaled2.E() > 0.0) || scaled1.M2() > 1.0e-10 || scaled2.M2() > 1.0e-10 ||
      closure >= 1.0e-9 * std::max(1.0, mass2)) {
    return false;
  }
  q1 = scaled1;
  q2 = scaled2;
  return true;
}

// Build one incoming hard parton for a selected beam side
bool MHardDiffraction::BuildHardPartonForSide(bool first_side, M4Vec &parton) const {
  const bool  hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  const M4Vec beam      = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const bool  plus_side = first_side;

  if (!hard_diff) { return false; }

  const double xi   = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
  const double beta = first_side ? state.lts.diff_beta1 : state.lts.diff_beta2;
  const double t    = first_side ? state.lts.diff_t1 : state.lts.diff_t2;
  const double phi  = first_side ? state.lts.diff_phi1 : state.lts.diff_phi2;

  M4Vec leading;
  if (!kinematics::BuildForwardParticleXiT(beam, xi, t, phi, plus_side, leading)) { return false; }
  return BuildDiffractiveHardParton(beam, beam - leading, xi, beta, plus_side, parton);
}

// Build the tagged particle and physical exchange remnant on one side
bool MHardDiffraction::BuildTaggedRemnant(bool first_side, const M4Vec &parton, M4Vec &leading, M4Vec &remnant) const {
  const M4Vec  beam = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const double xi   = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
  const double t    = first_side ? state.lts.diff_t1 : state.lts.diff_t2;
  const double phi  = first_side ? state.lts.diff_phi1 : state.lts.diff_phi2;
  if (!kinematics::BuildForwardParticleXiT(beam, xi, t, phi, first_side, leading) || !AcceptRemnantMomentum(leading)) {
    return false;
  }
  remnant = beam - leading - parton;
  if (remnant.E() > 1.0e-12) { return AcceptRemnantMomentum(remnant); }
  return std::abs(remnant.E()) < 1.0e-12 && std::abs(remnant.Px()) < 1.0e-12 && std::abs(remnant.Py()) < 1.0e-12 &&
         std::abs(remnant.Pz()) < 1.0e-12;
}

// Check one hard parton against its physical beam-remnant support
bool MHardDiffraction::AcceptHardParton(bool first_side, const M4Vec &parton) const {
  if (!(parton.E() > 0.0)) { return false; }
  const bool hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  if (!hard_diff) {
    const M4Vec &beam = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
    return AcceptRemnantMomentum(beam - parton);
  }
  M4Vec leading;
  M4Vec remnant;
  return BuildTaggedRemnant(first_side, parton, leading, remnant);
}

// Build one four-vector from light-cone components and transverse momentum
M4Vec MHardDiffraction::LightConeVector(double plus, double minus, double px, double py) const {
  return M4Vec(px, py, 0.5 * (plus - minus), 0.5 * (plus + minus));
}

// Build one spacelike hard parton and a physical lightlike Pomeron remnant
bool MHardDiffraction::BuildDiffractiveHardParton(const M4Vec &beam, const M4Vec &exchange, double xi, double beta,
                                                  bool plus_side, M4Vec &parton) const {
  if (!(xi > 0.0) || !(beta > 0.0) || beta > 1.0) { return false; }

  /*
   * A spacelike Pomeron cannot split into two future-directed on-shell objects
   * while conserving four-momentum. Build the DPDF light-cone split with
   * q^2 = beta t and a lightlike remnant. Double diffraction subsequently
   * rescales both hard legs by one Lorentz scalar so the hard mass is exact.
   */
  if (std::is_eq(beta <=> 1.0)) {
    parton = exchange;
  } else {
    const double remnant_px  = (1.0 - beta) * exchange.Px();
    const double remnant_py  = (1.0 - beta) * exchange.Py();
    const double remnant_pt2 = remnant_px * remnant_px + remnant_py * remnant_py;
    M4Vec        remnant;

    if (plus_side) {
      const double plus = (1.0 - beta) * exchange.LightconePos();
      if (!(plus > 0.0)) { return false; }
      remnant = LightConeVector(plus, remnant_pt2 / plus, remnant_px, remnant_py);
    } else {
      const double minus = (1.0 - beta) * exchange.LightconeNeg();
      if (!(minus > 0.0)) { return false; }
      remnant = LightConeVector(remnant_pt2 / minus, minus, remnant_px, remnant_py);
    }
    parton = exchange - remnant;
  }

  const double x =
      plus_side ? parton.LightconePos() / beam.LightconePos() : parton.LightconeNeg() / beam.LightconeNeg();
  const double virtuality = beta * exchange.M2();
  return std::isfinite(parton.E()) && parton.E() > 0.0 && std::abs(x - beta * xi) < 1.0e-8 &&
         std::abs(parton.M2() - virtuality) < 1.0e-8 * std::max(1.0, std::abs(virtuality));
}

// Resolve the ordinary proton remnant mass for this process
double MHardDiffraction::ProtonRemnantMass() const {
  return hard_pomeron_pdf->RemnantMass();
}

// Solve (beam + qdiff - remnant)^2 = mass2 with a forward massive remnant
bool MHardDiffraction::BuildProtonHardParton(bool first_side, const M4Vec &qdiff, double mass2, M4Vec &parton) const {
  const M4Vec &beam = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const M4Vec total = beam + qdiff;
  const long double large = first_side ? total.LightconePos() : total.LightconeNeg();
  const long double small = first_side ? total.LightconeNeg() : total.LightconePos();
  const long double mr2 = pow2(ProtonRemnantMass());
  const long double mt2 = static_cast<long double>(mass2) + total.Pt2();
  const long double b = large * small + mr2 - mt2;
  const long double disc = b * b - 4 * large * small * mr2;
  if (!(large > 0.0L) || !(small > 0.0L) || !(b > 0.0L) || !(disc > 0.0L)) { return false; }

  // The larger root keeps the remnant in the ordinary beam hemisphere
  const long double recoil = (b + std::sqrt(disc)) / (2 * small);
  const long double beam_large = first_side ? beam.LightconePos() : beam.LightconeNeg();
  const long double beam_small = first_side ? beam.LightconeNeg() : beam.LightconePos();
  const double qlarge = static_cast<double>(beam_large - recoil);
  const double qsmall = static_cast<double>(beam_small - mr2 / recoil);
  if (!(qlarge > 0.0) || !(qlarge < beam_large) || !(qsmall < 0.0)) { return false; }
  parton = first_side ? LightConeVector(qlarge, qsmall, 0.0, 0.0)
                      : LightConeVector(qsmall, qlarge, 0.0, 0.0);
  return std::isfinite(parton.E()) && parton.E() > 0.0 && AcceptRemnantMomentum(beam - parton);
}

// Split one ordinary remnant with a normalized isotropic two-body distribution
bool MHardDiffraction::BuildProtonRemnant(bool first_side, MDecayBranch &branch) const {
  const int beam_pdg = first_side ? state.lts.beam1.pdg : state.lts.beam2.pdg;
  const int parton_id = first_side ? state.lts.id1 : state.lts.id2;
  const auto ids = ProtonRemnantIDs(parton_id, beam_pdg);
  if (ids[0] == 0 || !AcceptRemnantMomentum(branch.p4)) { return false; }
  if (parton_id == 0) {
    branch.legs.push_back(MakeRemnantBranch(ids[0], branch.p4));
    return true;
  }

  // Transport the same Born split through every screening momentum
  const M4Vec &parent = state.lts.screening.active
                           ? state.lts.pfinal_orig[first_side ? 1 : 2] : branch.p4;
  const double mass = parent.M();
  const double m1 = HardParticle(ids[0]).mass;
  const double m2 = HardParticle(ids[1]).mass;
  if (!(mass > m1 + m2)) { return false; }
  M4Vec axis = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  kinematics::LorentzBoost(parent, mass, axis, -1);
  if (!(axis.E() > 0.0) || !(axis.P3mod() > 0.0)) { return false; }
  axis.SetE(0.0);
  axis /= axis.P3mod();
  const auto &forward = state.lts.screening.active ? state.lts.pfinal_orig : state.lts.pfinal;
  M4Vec plane = (first_side ? state.lts.pbeam2 : state.lts.pbeam1) - forward[first_side ? 2 : 1];
  kinematics::LorentzBoost(parent, mass, plane, -1);
  plane.SetE(0.0);
  plane -= axis * plane.Dot3(axis);
  if (!(plane.P3mod() > 1.0e-12)) {
    plane = std::abs(axis.Px()) < 0.9 ? M4Vec(1, 0, 0, 0) : M4Vec(0, 1, 0, 0);
    plane -= axis * plane.Dot3(axis);
  }
  plane /= plane.P3mod();
  const double cosine = 2.0 * remnant_angles[0] - 1.0;
  const double sine = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
  const double phi = 2.0 * math::PI * remnant_angles[1];
  const double p = kinematics::DecayMomentum(mass, m1, m2);
  M4Vec first = (axis * cosine + (plane * std::cos(phi) + axis.Cross3(plane) * std::sin(phi)) * sine) * p;
  first.SetE(std::hypot(p, m1));
  M4Vec second = first * -1.0;
  second.SetE(std::hypot(p, m2));
  kinematics::LorentzBoost(branch.p4, mass, first, 1);
  kinematics::LorentzBoost(branch.p4, mass, second, 1);
  if (!AcceptRemnantMomentum(first) || !AcceptRemnantMomentum(second) ||
      !CheckEMC(branch.p4 - first - second)) { return false; }
  branch.legs.push_back(MakeRemnantBranch(ids[0], first));
  branch.legs.push_back(MakeRemnantBranch(ids[1], second));
  return true;
}

// Rebuild deterministic forward remnant decay branches
bool MHardDiffraction::BuildHardRemnantBranches(bool assign_colors) {
  if (state.lts.pfinal.size() < 3 || (state.lts.screening.active && state.lts.pfinal_orig.size() < 3)) { return false; }
  std::array<MDecayBranch, 2> branches;
  branches[0].p4 = state.lts.pfinal[1];
  branches[1].p4 = state.lts.pfinal[2];
  branches[0].p  = HardParticle(PDG::PDG_NSTAR);
  branches[1].p  = HardParticle(PDG::PDG_NSTAR);

  const bool ok =
      state.lts.screening.active ? BuildLoopHardRemnantBranches(branches) : BuildNominalHardRemnantBranches(branches);
  if (!ok) { return false; }

  state.lts.decayforward1 = std::move(branches[0]);
  state.lts.decayforward2 = std::move(branches[1]);
  return !assign_colors || AssignHardColorFlow();
}

// Rebuild nominal t-dependent forward remnant branches
bool MHardDiffraction::BuildNominalHardRemnantBranches(std::array<MDecayBranch, 2> &branches) {
  return BuildNominalHardRemnantSide(true, branches[0]) && BuildNominalHardRemnantSide(false, branches[1]);
}

// Rebuild one nominal forward remnant side
bool MHardDiffraction::BuildNominalHardRemnantSide(bool first_side, MDecayBranch &branch) {
  const M4Vec beam      = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const M4Vec parent    = first_side ? state.lts.pfinal[1] : state.lts.pfinal[2];
  const bool  hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  const int   beam_pdg  = first_side ? state.lts.beam1.pdg : state.lts.beam2.pdg;
  const int   parton_id = first_side ? state.lts.id1 : state.lts.id2;

  if (hard_diff) {
    M4Vec leading;
    M4Vec remnant;
    if (!BuildTaggedRemnant(first_side, beam - parent, leading, remnant)) { return false; }
    branch.legs.push_back(MakeRemnantBranch(beam_pdg, leading));
    if (remnant.E() > 1.0e-12) { branch.legs.push_back(MakeRemnantBranch(RemnantPartnerPDG(parton_id), remnant)); }
  } else {
    if (!BuildProtonRemnant(first_side, branch)) { return false; }
  }

  M4Vec forward_sum(0, 0, 0, 0);
  for (const auto &leg : branch.legs) { forward_sum += leg.p4; }
  return CheckEMC(parent - forward_sum);
}

// Rebuild loop-deformed forward remnant branches by momentum fractions
bool MHardDiffraction::BuildLoopHardRemnantBranches(std::array<MDecayBranch, 2> &branches) {
  return BuildLoopHardRemnantSide(true, branches[0]) && BuildLoopHardRemnantSide(false, branches[1]);
}

// Rebuild one loop-deformed forward remnant side
bool MHardDiffraction::BuildLoopHardRemnantSide(bool first_side, MDecayBranch &branch) {
  const M4Vec beam      = first_side ? state.lts.pbeam1 : state.lts.pbeam2;
  const M4Vec parent    = first_side ? state.lts.pfinal[1] : state.lts.pfinal[2];
  const bool  hard_diff = first_side ? state.lts.hard_diff1 : state.lts.hard_diff2;
  const int   beam_pdg  = first_side ? state.lts.beam1.pdg : state.lts.beam2.pdg;
  const int   parton_id = first_side ? state.lts.id1 : state.lts.id2;

  if (!hard_diff) {
    return BuildProtonRemnant(first_side, branch);
  }

  const double xi          = first_side ? state.lts.diff_xi1 : state.lts.diff_xi2;
  const double t           = first_side ? state.lts.diff_t1 : state.lts.diff_t2;
  const double phi         = first_side ? state.lts.diff_phi1 : state.lts.diff_phi2;
  const M4Vec  parent_orig = first_side ? state.lts.pfinal_orig[1] : state.lts.pfinal_orig[2];

  M4Vec leading_orig;
  if (!kinematics::BuildForwardParticleXiT(beam, xi, t, phi, first_side, leading_orig)) { return false; }

  M4Vec leading = leading_orig;
  M4Vec remnant = parent_orig - leading_orig;
  if (!AcceptRemnantMomentum(leading)) { return false; }

  // At beta = 1 the exact leading baryon is the complete forward system
  if (remnant.E() <= 1.0e-12) {
    leading = parent;
    remnant = M4Vec(0.0, 0.0, 0.0, 0.0);
    if (!AcceptRemnantMomentum(leading)) { return false; }
    branch.legs.push_back(MakeRemnantBranch(beam_pdg, leading));
    return CheckEMC(parent - leading);
  }
  if (!AcceptRemnantMomentum(remnant)) { return false; }

  // Transport the physical Born two-body forward system to the loop-deformed
  // parent
  kinematics::LorentzBoost(parent_orig, parent_orig.M(), leading, -1);
  kinematics::LorentzBoost(parent_orig, parent_orig.M(), remnant, -1);
  kinematics::LorentzBoost(parent, parent.M(), leading, 1);
  kinematics::LorentzBoost(parent, parent.M(), remnant, 1);
  if (!AcceptRemnantMomentum(leading) || !AcceptRemnantMomentum(remnant)) { return false; }
  if (!CheckEMC(parent - (leading + remnant))) { return false; }

  branch.legs.push_back(MakeRemnantBranch(beam_pdg, leading));
  if (remnant.E() > 1.0e-12) { branch.legs.push_back(MakeRemnantBranch(RemnantPartnerPDG(parton_id), remnant)); }
  return true;
}

// Assign shower-compatible color tags to hard partons and remnants
bool MHardDiffraction::AssignHardColorFlow() {
  if (!hard_color_flows.empty()) {
    if (!HardColorFlowsMatch(state.lts.hard_color_flows, hard_color_flows)) { return false; }
    const auto selected =
        mg5helas::SelectHardColorFlow(state.lts.hard_color_flows, state.random, hard_selected_channel_index);
    if (selected.has_value()) {
      std::vector<MColorFlow> assignment;
      assignment.reserve(state.lts.hard_color_flows[*selected].external.size());
      for (const auto &leg : state.lts.hard_color_flows[*selected].external) {
        assignment.push_back({leg.color, leg.anticolor});
      }
      return ApplyHardColorFlowCandidate(assignment);
    }
    return false;
  }

  // Construct the analytic color assignment for non-MG5 hard subprocesses
  struct ColorSlot {
    int                *tag   = nullptr;
    const MDecayBranch *owner = nullptr;
  };

  // Close the hard parton colors before resolving either composite remnant
  std::array<MDecayBranch, 2> crossed;
  const std::array<int, 2> incoming = {state.lts.id1, state.lts.id2};
  for (const auto &side : indices(crossed)) {
    crossed[side].p.pdg = incoming[side] == PDG::PDG_gluon ? incoming[side] : -incoming[side];
  }
  std::vector<MDecayBranch *> leaves = {&crossed[0], &crossed[1]};
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, leaves); }

  std::vector<ColorSlot> color_slots;
  std::vector<ColorSlot> anticolor_slots;
  for (auto *leaf : leaves) {
    leaf->p.color_flow.clear();
    const int pdg = leaf->p.pdg;
    if (!IsColoredParton(pdg)) { continue; }

    if (pdg == PDG::PDG_gluon) {
      color_slots.push_back({&leaf->p.color_flow.flow1, leaf});
      anticolor_slots.push_back({&leaf->p.color_flow.flow2, leaf});
    } else if (IsDiquark(pdg)) {
      if (pdg > 0) {
        anticolor_slots.push_back({&leaf->p.color_flow.flow2, leaf});
      } else {
        color_slots.push_back({&leaf->p.color_flow.flow1, leaf});
      }
    } else if (pdg > 0) {
      color_slots.push_back({&leaf->p.color_flow.flow1, leaf});
    } else {
      anticolor_slots.push_back({&leaf->p.color_flow.flow2, leaf});
    }
  }

  if (color_slots.size() != anticolor_slots.size()) { return false; }

  std::vector<int>                                      anticolor_match(anticolor_slots.size(), -1);
  std::function<bool(std::size_t, std::vector<bool> &)> match_color = [&](std::size_t        color_index,
                                                                          std::vector<bool> &seen) {
    for (const auto &anti_index : indices(anticolor_slots)) {
      if (seen[anti_index] || color_slots[color_index].owner == anticolor_slots[anti_index].owner) { continue; }
      seen[anti_index] = true;
      if (anticolor_match[anti_index] < 0 || match_color(static_cast<std::size_t>(anticolor_match[anti_index]), seen)) {
        anticolor_match[anti_index] = static_cast<int>(color_index);
        return true;
      }
    }
    return false;
  };

  for (const auto &color_index : indices(color_slots)) {
    std::vector<bool> seen(anticolor_slots.size(), false);
    if (!match_color(color_index, seen)) { return false; }
  }

  constexpr int first_tag = 501;
  for (const auto &anti_index : indices(anticolor_slots)) {
    const std::size_t color_index    = static_cast<std::size_t>(anticolor_match[anti_index]);
    const int         tag            = first_tag + static_cast<int>(color_index);
    *color_slots[color_index].tag    = tag;
    *anticolor_slots[anti_index].tag = tag;
  }
  mg5helas::ExternalColorFlow external;
  external.reserve(leaves.size());
  for (const auto &i : indices(leaves)) {
    const auto &flow = leaves[i]->p.color_flow;
    external.push_back(i < 2 ? mg5helas::ColorFlowLeg{flow.flow2, flow.flow1}
                             : mg5helas::ColorFlowLeg{flow.flow1, flow.flow2});
  }
  std::vector<MColorFlow> candidate;
  return BuildGeneratedHardColorFlowCandidate(external, candidate) && ApplyHardColorFlowCandidate(candidate);
}

// Enumerate leading-color assignments in stable-leaf order
std::vector<std::vector<MColorFlow>> MHardDiffraction::BuildHardColorFlowCandidates() {
  struct Slot {
    std::size_t owner = 0;
  };

  std::vector<MDecayBranch *> leaves;
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, leaves); }
  CollectStableLeaves(state.lts.decayforward1, leaves);
  CollectStableLeaves(state.lts.decayforward2, leaves);

  std::vector<Slot> colors;
  std::vector<Slot> anticolors;
  for (const auto &i : indices(leaves)) {
    const int pdg = leaves[i]->p.pdg;
    if (!IsColoredParton(pdg)) { continue; }
    if (pdg == PDG::PDG_gluon) {
      colors.push_back({i});
      anticolors.push_back({i});
    } else if (IsDiquark(pdg)) {
      if (pdg > 0) {
        anticolors.push_back({i});
      } else {
        colors.push_back({i});
      }
    } else if (pdg > 0) {
      colors.push_back({i});
    } else {
      anticolors.push_back({i});
    }
  }
  if (colors.size() != anticolors.size()) { return {}; }

  std::vector<std::vector<MColorFlow>> candidates;
  std::vector<int>                     match(colors.size(), -1);
  std::vector<bool>                    used(anticolors.size(), false);
  std::function<void(std::size_t)>     enumerate = [&](std::size_t color_index) {
    if (color_index == colors.size()) {
      std::vector<MColorFlow>               candidate(leaves.size());
      std::vector<std::vector<std::size_t>> graph(leaves.size());
      for (const auto &i : indices(colors)) {
        const std::size_t anti_index  = static_cast<std::size_t>(match[i]);
        const std::size_t color_owner = colors[i].owner;
        const std::size_t anti_owner  = anticolors[anti_index].owner;
        const int         tag         = 501 + static_cast<int>(i);
        candidate[color_owner].flow1  = tag;
        candidate[anti_owner].flow2   = tag;
        graph[color_owner].push_back(anti_owner);
        graph[anti_owner].push_back(color_owner);
      }

      bool has_quark_end = false;
      for (const auto &i : indices(leaves)) {
        const int pdg = leaves[i]->p.pdg;
        has_quark_end = has_quark_end || (IsColoredParton(pdg) && pdg != PDG::PDG_gluon);
      }
      if (has_quark_end) {
        std::vector<bool> seen(leaves.size(), false);
        for (const auto &root : indices(leaves)) {
          if (seen[root] || !IsColoredParton(leaves[root]->p.pdg)) { continue; }
          bool                     component_gluon = false;
          bool                     component_quark = false;
          std::vector<std::size_t> stack           = {root};
          seen[root]                               = true;
          while (!stack.empty()) {
            const std::size_t node = stack.back();
            stack.pop_back();
            component_gluon = component_gluon || leaves[node]->p.pdg == PDG::PDG_gluon;
            component_quark = component_quark || leaves[node]->p.pdg != PDG::PDG_gluon;
            for (const auto next : graph[node]) {
              if (!seen[next]) {
                seen[next] = true;
                stack.push_back(next);
              }
            }
          }
          // Leading-color gluons lie on an open quark line when quark line ends
          // exist
          if (component_gluon && !component_quark) { return; }
        }
      }
      candidates.push_back(std::move(candidate));
      return;
    }

    for (const auto &anti : indices(anticolors)) {
      if (used[anti] || colors[color_index].owner == anticolors[anti].owner) { continue; }
      used[anti]         = true;
      match[color_index] = static_cast<int>(anti);
      enumerate(color_index + 1);
      used[anti] = false;
    }
  };
  enumerate(0);
  return candidates;
}

// Convert one generated symbolic color pair to event color tags
MColorFlow MHardDiffraction::ConvertGeneratedColorFlowLeg(const mg5helas::ColorFlowLeg &leg,
                                                          std::map<int, int>           &tag_map) const {
  const auto convert_tag = [&tag_map](int tag) {
    if (tag == 0) { return 0; }
    const auto [entry, inserted] = tag_map.try_emplace(std::abs(tag), 501 + static_cast<int>(tag_map.size()));
    return (tag > 0 ? 1 : -1) * entry->second;
  };
  return {convert_tag(leg.color), convert_tag(leg.anticolor)};
}

// Compute true when one event color pair matches the particle representation
bool MHardDiffraction::HardColorFlowMatchesParticle(const MColorFlow &flow, const MParticle &particle) const {
  const int pdg = particle.pdg;
  if (particle.color == 6) { return flow.flow1 >= 501 && flow.flow2 <= -501; }
  if (particle.color == -6) { return flow.flow1 <= -501 && flow.flow2 >= 501; }
  if (particle.color == 8 || pdg == PDG::PDG_gluon) { return flow.flow1 >= 501 && flow.flow2 >= 501 && flow.flow1 != flow.flow2; }
  if (IsDiquark(pdg)) { return pdg > 0 ? flow.flow1 == 0 && flow.flow2 >= 501 : flow.flow1 >= 501 && flow.flow2 == 0; }
  if (std::abs(particle.color) == 3 || (std::abs(pdg) >= 1 && std::abs(pdg) <= 6)) {
    return (particle.color != 0 ? particle.color > 0 : pdg > 0) ? flow.flow1 >= 501 && flow.flow2 == 0 : flow.flow1 == 0 && flow.flow2 >= 501;
  }
  return flow.empty();
}

// Distribute the crossed color across a singlet baryon, diquark or resolved gluon remnant
bool MHardDiffraction::MapRemnantColor(const std::vector<MDecayBranch *> &leaves, const MColorFlow &crossed,
                                      std::size_t offset, int &next_tag, std::vector<MColorFlow> &candidate) const {
  std::vector<std::size_t> colored;
  for (const auto &i : indices(leaves)) {
    if (IsColoredParton(leaves[i]->p.pdg)) { colored.push_back(i); }
  }
  if (crossed.empty()) { return colored.empty(); }
  if (colored.size() == 1) {
    const auto i = colored.front();
    if (!HardColorFlowMatchesParticle(crossed, leaves[i]->p)) { return false; }
    candidate[offset + i] = crossed;
    return true;
  }
  if (colored.size() != 2) { return false; }

  const auto i = colored[0];
  const auto j = colored[1];
  if (crossed.flow1 != 0 && crossed.flow2 != 0) {
    // A gluon removed from a baryon leaves a triplet and an antitriplet
    const MColorFlow color = {crossed.flow1, 0};
    const MColorFlow anti = {0, crossed.flow2};
    const bool first_color = HardColorFlowMatchesParticle(color, leaves[i]->p);
    candidate[offset + i] = first_color ? color : anti;
    candidate[offset + j] = first_color ? anti : color;
  } else {
    // Insert the remnant gluon on the open diquark color line
    if (!IsDiquark(leaves[i]->p.pdg) || leaves[j]->p.pdg != PDG::PDG_gluon) { return false; }
    const int tag = next_tag++;
    candidate[offset + i] = crossed.flow1 != 0 ? MColorFlow{tag, 0} : MColorFlow{0, tag};
    candidate[offset + j] = crossed.flow1 != 0 ? MColorFlow{crossed.flow1, tag} : MColorFlow{tag, crossed.flow2};
  }
  return HardColorFlowMatchesParticle(candidate[offset + i], leaves[i]->p) &&
         HardColorFlowMatchesParticle(candidate[offset + j], leaves[j]->p);
}

// Validate one complete candidate without changing the event color state
bool MHardDiffraction::ValidateHardColorFlowCandidate(const std::vector<MColorFlow> &candidate) {
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, leaves); }
  CollectStableLeaves(state.lts.decayforward1, leaves);
  CollectStableLeaves(state.lts.decayforward2, leaves);
  if (candidate.size() != leaves.size()) { return false; }

  std::map<int, std::array<std::size_t, 2>> tag_counts;
  for (const auto &i : indices(leaves)) {
    if (!HardColorFlowMatchesParticle(candidate[i], leaves[i]->p)) { return false; }
    if (candidate[i].flow1 != 0) { ++tag_counts[std::abs(candidate[i].flow1)][candidate[i].flow1 < 0 ? 1 : 0]; }
    if (candidate[i].flow2 != 0) { ++tag_counts[std::abs(candidate[i].flow2)][candidate[i].flow2 > 0 ? 1 : 0]; }
  }
  for (const auto &entry : tag_counts) {
    if (entry.first < 501 || entry.second[0] != 1 || entry.second[1] != 1) { return false; }
  }
  return true;
}

// Map one exact generated external flow onto central and remnant leaves
bool MHardDiffraction::BuildGeneratedHardColorFlowCandidate(const mg5helas::ExternalColorFlow &external_flow,
                                                            std::vector<MColorFlow>           &candidate) {
  candidate.clear();
  if (external_flow.empty()) {
    const auto candidates = BuildHardColorFlowCandidates();
    if (candidates.size() != 1 || !ValidateHardColorFlowCandidate(candidates.front())) { return false; }
    candidate = candidates.front();
    return true;
  }

  std::vector<MDecayBranch *> central;
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, central); }
  if (external_flow.size() != central.size() + 2) { return false; }

  std::vector<MDecayBranch *> forward1;
  std::vector<MDecayBranch *> forward2;
  CollectStableLeaves(state.lts.decayforward1, forward1);
  CollectStableLeaves(state.lts.decayforward2, forward2);

  std::map<int, int>      tag_map;
  std::vector<MColorFlow> converted;
  converted.reserve(external_flow.size());
  for (const auto &leg : external_flow) {
    converted.push_back(ConvertGeneratedColorFlowLeg(leg, tag_map));
  }

  candidate.assign(central.size() + forward1.size() + forward2.size(), {});
  for (const auto &i : indices(central)) { candidate[i] = converted[i + 2]; }

  const std::array<std::vector<MDecayBranch *> *, 2> forward = {&forward1, &forward2};
  std::size_t offset = central.size();
  int next_tag = 501 + static_cast<int>(tag_map.size());
  for (const auto &side : indices(forward)) {
    const MColorFlow crossed = {converted[side].flow2, converted[side].flow1};
    if (!MapRemnantColor(*forward[side], crossed, offset, next_tag, candidate)) {
      candidate.clear();
      return false;
    }
    offset += forward[side]->size();
  }

  if (!ValidateHardColorFlowCandidate(candidate)) {
    candidate.clear();
    return false;
  }
  return true;
}

// Compute true when two candidate tables contain identical ordered tags
bool MHardDiffraction::HardColorFlowsMatch(const std::vector<mg5helas::HardColorFlow> &first,
                                           const std::vector<mg5helas::HardColorFlow> &second) const {
  if (first.size() != second.size()) { return false; }
  for (const auto &row : indices(first)) {
    if (first[row].channel != second[row].channel || first[row].external.size() != second[row].external.size()) {
      return false;
    }
    for (const auto &leg : indices(first[row].external)) {
      if (first[row].external[leg].color != second[row].external[leg].color ||
          first[row].external[leg].anticolor != second[row].external[leg].anticolor) {
        return false;
      }
    }
  }
  return true;
}

// Apply one generated color-flow candidate to all stable leaves
bool MHardDiffraction::ApplyHardColorFlowCandidate(const std::vector<MColorFlow> &candidate) {
  if (DecayTreeHasColoredIntermediate(state.lts.decaytree) ||
      DecayBranchHasColoredIntermediate(state.lts.decayforward1) ||
      DecayBranchHasColoredIntermediate(state.lts.decayforward2)) {
    return false;
  }
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : state.lts.decaytree) { CollectStableLeaves(branch, leaves); }
  CollectStableLeaves(state.lts.decayforward1, leaves);
  CollectStableLeaves(state.lts.decayforward2, leaves);
  if (candidate.size() != leaves.size() || !ValidateHardColorFlowCandidate(candidate)) { return false; }
  for (auto &branch : state.lts.decaytree) { ClearDecayBranchColorFlow(branch); }
  ClearDecayBranchColorFlow(state.lts.decayforward1);
  ClearDecayBranchColorFlow(state.lts.decayforward2);
  for (const auto &i : indices(leaves)) { leaves[i]->p.color_flow = candidate[i]; }
  return true;
}

// Collect stable leaves from one mutable decay branch
void MHardDiffraction::CollectStableLeaves(MDecayBranch &branch, std::vector<MDecayBranch *> &leaves) const {
  if (branch.legs.empty()) {
    leaves.push_back(&branch);
    return;
  }
  for (auto &leg : branch.legs) { CollectStableLeaves(leg, leaves); }
}

// Compute true for partons carrying QCD color
bool MHardDiffraction::IsColoredParton(int pdg) const {
  return pdg == PDG::PDG_gluon || (std::abs(pdg) >= 1 && std::abs(pdg) <= 5) || IsDiquark(pdg);
}

// Compute true for colored diquark remnant ids
bool MHardDiffraction::IsDiquark(int pdg) const {
  const int apdg = std::abs(pdg);
  const int spin = apdg % 10;
  return apdg >= 1000 && apdg < 6000 && (spin == 1 || spin == 3);
}

// Compute a standard colored remnant partner for an extracted parton
int MHardDiffraction::RemnantPartnerPDG(int pdg) const {
  if (pdg == PDG::PDG_gluon) { return PDG::PDG_gluon; }
  if (std::abs(pdg) >= 1 && std::abs(pdg) <= 5) { return -pdg; }
  return PDG::PDG_fragment;
}

// Compute a minimal flavour-conserving remnant with effective light-quark recombination
// [REFERENCE: T. Sjostrand et al., hep-ph/0603175, section 11, https://arxiv.org/abs/hep-ph/0603175]
std::array<int, 2> MHardDiffraction::ProtonRemnantIDs(int pdg, int beam_pdg) const {
  // Before amplitude evaluation the incoming flavour is not yet selected
  if (pdg == 0) { return {PDG::PDG_fragment, 0}; }
  if (std::abs(beam_pdg) != PDG::PDG_p) { return {}; }

  const int sign = (beam_pdg > 0) ? 1 : -1;
  const int id   = sign * pdg;

  if (pdg == PDG::PDG_gluon) { return {sign * 2, sign * 2101}; }
  if (id == 2) { return {sign * 2101, PDG::PDG_gluon}; }
  if (id == 1) { return {sign * 2203, PDG::PDG_gluon}; }
  // Recombine the sea companion with a valence constituent of the color octet core
  if (id == 3) { return {sign * 321, sign * 2101}; }
  if (id == 4) { return {-sign * 421, sign * 2101}; }
  if (id == 5) { return {sign * 521, sign * 2101}; }
  if (id == -1) { return {sign * 2112, sign * 2}; }
  if (id == -2) { return {sign * 2212, sign * 2}; }
  if (id == -3) { return {sign * 3122, sign * 2}; }
  if (id == -4) { return {sign * 4122, sign * 2}; }
  if (id == -5) { return {sign * 5122, sign * 2}; }
  return {};
}

// Build one stable remnant branch
MDecayBranch MHardDiffraction::MakeRemnantBranch(int pdg, const M4Vec &p4) const {
  MDecayBranch branch;
  branch.p  = HardParticle(pdg);
  branch.p4 = p4;
  branch.m_offshell = p4.M();
  return branch;
}

// Compute true when one remnant four-vector is usable in an event record
bool MHardDiffraction::AcceptRemnantMomentum(const M4Vec &p4) const {
  // Permit numerical lightlike remnants but reject genuinely spacelike
  // event-record particles
  return std::isfinite(p4.E()) && std::isfinite(p4.Px()) && std::isfinite(p4.Py()) && std::isfinite(p4.Pz()) &&
         p4.E() > 0.0 && p4.M2() > -1.0e-8;
}

// Resolve a physical particle definition or an internal remnant container
MParticle MHardDiffraction::HardParticle(int pdg) const {
  if (pdg == 0 || std::abs(pdg) == PDG::PDG_NSTAR || pdg == PDG::PDG_fragment || IsDiquark(pdg)) {
    MParticle particle;
    particle.pdg  = pdg;
    particle.name = "hard-remnant";
    if (IsDiquark(pdg)) {
      const int sign = pdg > 0 ? 1 : -1;
      particle.mass = 2.0 * PDG::mp / 3.0;
      particle.chargeX3 = sign * (std::abs(pdg) == 2203 ? 4 : 1);
      particle.color = -3 * sign;
      particle.spinX2 = std::abs(pdg) % 10 - 1;
    }
    return particle;
  }

  try {
    return state.lts.PDG.FindByPDG(pdg);
  } catch (const std::exception &) {
    throw PhaseSpaceFailure("MHardDiffraction::HardParticle: missing PDG data for " + std::to_string(pdg));
  }
}

// Compute the sampled hard-diffraction integral volume
double MHardDiffraction::HardIntegralVolume() const { return hard_integral_volume; }

// Compute the hard-scattering phase-space Jacobian
double MHardDiffraction::HardPhaseSpaceWeight() const { return 1.0; }

// Calculate pure phase-space decay width
void MHardDiffraction::DecayWidthPS(double &exact) const { exact = CentralDecayWidthPS(); }

// Compute the supported incoming hard-parton flavours
std::vector<int> MHardDiffraction::HardPartonFlavours() const {
  // Expand configured quark species to both signs while retaining each species
  // only once
  std::vector<int> flavours;
  for (const int pid : hard_pomeron_pdf->PartonFlavours()) {
    if (pid == PDG::PDG_gluon) {
      if (std::find(flavours.begin(), flavours.end(), pid) == flavours.end()) { flavours.push_back(pid); }
    } else {
      const int q = std::abs(pid);
      if (std::find(flavours.begin(), flavours.end(), q) == flavours.end()) {
        flavours.push_back(q);
        flavours.push_back(-q);
      }
    }
  }
  return flavours;
}

// Select and publish one hard-parton channel to state.lts metadata
bool MHardDiffraction::SelectHardPartonChannel(const std::vector<double> &component_amp2) {
  const std::size_t born_component_count = hard_component_channel_indices.size();
  if ((born_component_count == 0 && !component_amp2.empty()) ||
      (born_component_count != 0 && component_amp2.size() % born_component_count != 0)) {
    return false;
  }

  std::vector<double> channel_amp2(hard_parton_channels.size(), 0.0);
  for (const auto &i : gra::aux::indices(component_amp2)) {
    const double weight = component_amp2[i];
    if (std::isfinite(weight) && weight > 0.0) {
      // Dense proton spin and Good Walker blocks keep Born components fastest
      const std::size_t channel = hard_component_channel_indices[i % born_component_count];
      if (channel >= channel_amp2.size()) { return false; }
      channel_amp2[channel] += weight;
    }
  }

  double total = 0.0;
  for (const auto weight : channel_amp2) {
    if (std::isfinite(weight) && weight > 0.0) { total += weight; }
  }

  if (!(total > 0.0) || hard_parton_channels.empty()) {
    state.lts.id1               = 0;
    state.lts.id2               = 0;
    state.lts.pdf_xf1           = 0.0;
    state.lts.pdf_xf2           = 0.0;
    hard_selected_channel_index = 0;
    return true;
  }

  const double target   = state.random.U(0.0, total);
  double       sum      = 0.0;
  std::size_t  selected = hard_parton_channels.size() - 1;
  for (const auto &i : gra::aux::indices(channel_amp2)) {
    const double weight = channel_amp2[i];
    if (!std::isfinite(weight) || weight <= 0.0) { continue; }
    sum += weight;
    if (target <= sum) {
      selected = i;
      break;
    }
  }

  const auto &channel         = hard_parton_channels[selected];
  hard_selected_channel_index = selected;
  state.lts.id1               = channel.id1;
  state.lts.id2               = channel.id2;
  state.lts.pdf_xf1           = channel.pdf_xf1;
  state.lts.pdf_xf2           = channel.pdf_xf2;
  return true;
}

// Compute a proton parton density f_i/p(x,Q2)
double MHardDiffraction::ProtonPartonDensity(int pid, double x, double Q2) {
  if (!(x > 0.0) || x > 1.0) { return 0.0; }
  try {
    const double xfx = state.lts.GlobalPdfPtr->xfxQ2(pid, x, Q2);
    const double f   = xfx / x;
    return std::isfinite(f) ? f : 0.0;
  } catch (...) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return 0.0;
  }
}

// Compute alpha_s for the current hard scale
double MHardDiffraction::AlphaQCD(double Q2) {
  try {
    return hard_pomeron_pdf->AlphaS(Q2);
  } catch (...) {
    evaluation_status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return 0.0;
  }
}

// Print the hard-diffraction setup
void MHardDiffraction::PrintInit(bool silent) const {
  if (!silent) {
    PrintSetup();

    const std::string proton1 = (ProcPtr.ISTATE == "IPIP") ? "-----------DPDF-------->" : "--------DPDF/pdf------>";
    const std::string proton2 = (ProcPtr.ISTATE == "IPIP") ? "-----------DPDF-------->" : "--------pdf/DPDF======>";
    const std::vector<std::string> feynmangraph = {"||          ", "||          ", "xx--------->", "||          ",
                                                   "||          "};

    std::cout << proton1 << std::endl;
    for (const auto &row : feynmangraph) {
      if (state.screening) {
        std::cout << rang::fg::red << "     **    " << rang::style::reset;
      } else {
        std::cout << rang::fg::red << "           " << rang::style::reset;
      }
      std::cout << row << std::endl;
    }
    std::cout << proton2 << std::endl;
    std::cout << std::endl;
    std::cout << rang::style::bold << "Hard Pomeron PDF:" << rang::style::reset << std::endl;
    if (hard_pomeron_pdf != nullptr) { hard_pomeron_pdf->PrintSummary(std::cout); }
    std::cout << std::endl;
    PrintFiducialCuts();
  }
}

}  // namespace gra
