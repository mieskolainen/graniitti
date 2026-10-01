// Abstract process class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <mutex>
#include <thread>
#include <stdexcept>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Process/MEventRecord.h"
#include "Graniitti/Process/MForwardExcitation.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Tech/MException.h"

// HepMC3
#include "HepMC3/GenEvent.h"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

// Serialize each diagnostic across workers and keep logging failures nonfatal
void MProcess::PrintFailure(const char *type, const char *message) const noexcept {
  if (!debug) { return; }
  try {
    std::lock_guard<std::mutex> lock(gra::g_mutex);
    std::cerr << "[debug][" << std::this_thread::get_id() << "][" << type << "] " << message << std::endl;
  } catch (...) {}
}

// Separate invalid numerical weights from zero physical amplitude support
void MProcess::BookkeepAmplitudeWeight(double &W, MEventWeightState &aux) {
  if (aux.technical_failure) {
    W = 0.0;
    return;
  }
  if (!CheckInfNan(W)) {
    PrintFailure("TechnicalFailure", "MProcess::BookkeepAmplitudeWeight: non-finite event weight");
    aux.technical_failure = true;
    return;
  }
  if (W < 0.0) {
    PrintFailure("TechnicalFailure", "MProcess::BookkeepAmplitudeWeight: negative event weight");
    W                     = 0.0;
    aux.technical_failure = true;
    return;
  }
  if (math::IsZero(W)) { aux.amplitude_ok = false; }
}

// Report rejected screening nodes at the common process boundary
bool MProcess::EvaluateScreeningNode(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p) {
  const bool success = ComputeScreeningNode(p1p, p2p);
  if (!success) {
    PrintFailure("KinematicsFailure", "MProcess::EvaluateScreeningNode: loop or amplitude kinematics failed");
  }
  return success;
}

// Evaluate one process weight and apply mode-independent failure bookkeeping
double MProcess::EventWeight(const std::vector<double> &randvec, MEventWeightState &aux) {
  BeginEventEvaluation(aux);

  const unsigned int dimension = GetdLIPSDim();
  if (randvec.size() < dimension || !std::all_of(randvec.cbegin(), randvec.cbegin() + dimension, [](const double unit) {
        return std::isfinite(unit) && unit >= 0.0 && unit <= 1.0;
      })) {
    aux.kinematics_ok     = false;
    aux.technical_failure = true;
    PrintFailure("KinematicsFailure", "MProcess::EventWeight: invalid phase-space sampling coordinates");
    if (aux.qmetrics) { state.lts.qmetrics.Observe(0.0, aux.log_inverse_density, aux.sample_index); }
    ResetRejectedEventState();
    return 0.0;
  }

  double weight = 0.0;
  bool kinematics_reported = false;
  try {
    const double radiative_weight =
        state.radiative.GenerateISR(randvec, ProcPtr.LIPSDIM, state.lts.beam1, state.lts.beam2, state.lts.pbeam1,
                                    state.lts.pbeam2, state.lts.s, state.lts.sqrt_s);
    if (!std::isfinite(radiative_weight) || !(radiative_weight > 0.0)) {
      PrintFailure("KinematicsFailure", "MProcess::EventWeight: ISR kinematics or weight failed");
      kinematics_reported = true;
      aux.kinematics_ok = false;
      weight            = 0.0;
    } else {
      weight = radiative_weight * ComputeEventWeight(randvec, aux);
      if (aux.qmetrics) { aux.phase_weight *= radiative_weight; }
      BookkeepAmplitudeWeight(weight, aux);
      if (aux.Valid()) {
        // Finalize the accepted amplitude exactly once before FSR and cuts
        if (!SampleUPCFinalState()) { throw AmplitudeFailure("MProcess::EventWeight: UPC final-state sampling failed"); }
        if (!FinalizeAmplitudeEventState()) {
          throw AmplitudeFailure("MProcess::EventWeight: process final-state sampling failed");
        }
        amplitude_event_state_finalized = true;
        ClassifyFinalState(aux);
        if (aux.Valid() && state.nuclear_final.has_value()) {
          state.nuclear_event.emplace(HepMC3::Units::GEV, HepMC3::Units::MM);
          if (!BuildEventRecord(*state.nuclear_event)) {
            throw PhaseSpaceFailure("MProcess::EventWeight: nuclear hard event construction failed");
          }
          weight *= state.nuclear_final->Complete(*state.nuclear_event, state.random);
          BookkeepAmplitudeWeight(weight, aux);
        }
        if (!aux.Valid()) { weight = 0.0; }
      }
    }
  } catch (const AmplitudeFailure &error) {
    PrintFailure("AmplitudeFailure", error.what());
    aux.RecordAmplitudeFailure(error.what());
    weight = 0.0;
  } catch (const PhaseSpaceFailure &error) {
    PrintFailure("PhaseSpaceFailure", error.what());
    kinematics_reported = true;
    aux.kinematics_ok     = false;
    aux.technical_failure = true;
    weight                = 0.0;
  } catch (...) {
    ResetRejectedEventState();
    throw;
  }

  // Fill persistent histograms only after post-radiation classification
  if (!aux.kinematics_ok && !kinematics_reported) {
    PrintFailure("KinematicsFailure", "MProcess::EventWeight: phase-space construction returned false");
  }
  if (!aux.adaptation_mode) {
    const double totalweight = gra::statistics::ImportanceWeight(weight, aux.log_inverse_density);
    FillHistograms(totalweight, state.lts);
  }

  if (aux.qmetrics) {
    state.lts.qmetrics.Observe(aux.Valid() ? aux.phase_weight : 0.0, aux.log_inverse_density, aux.sample_index);
  }
  if (!aux.Valid()) { ResetRejectedEventState(); }
  return weight;
}

// Construct one event record and clear all event-local state on every exit
bool MProcess::EventRecord(HepMC3::GenEvent &evt) {
  try {
    if (!amplitude_event_state_finalized) {
      throw AmplitudeFailure("MProcess::EventRecord: event amplitude was not finalized by EventWeight");
    }
    bool success = true;
    if (state.nuclear_event.has_value()) {
      evt = *state.nuclear_event;
    } else {
      success = BuildEventRecord(evt);
    }
    state.nuclear_event.reset();
    if (!success) { PrintFailure("EventRecordFailure", "MProcess::EventRecord: event construction returned false"); }
    forward::ClearBranches(state.lts);
    amplitude_event_state_finalized = false;
    if (success) {
      ResetTransientAmplitudeState();
    } else {
      ResetFailedAmplitudeState();
    }
    return success;
  } catch (const AmplitudeFailure &error) {
    PrintFailure("AmplitudeFailure", error.what());
    ResetRejectedEventState();
    return false;
  } catch (const PhaseSpaceFailure &error) {
    PrintFailure("PhaseSpaceFailure", error.what());
    ResetRejectedEventState();
    return false;
  } catch (...) {
    ResetRejectedEventState();
    throw;
  }
}

// Add configured QED FSR photons to one accepted event record
bool MProcess::ApplyRadiation(HepMC3::GenEvent &evt) noexcept {
  if (state.radiative.GetConfig().fsr == radiative::Mode::YFS && !state.radiative.GetFSR().prepared) {
    PrintFailure("RadiativeFailure", "MProcess::ApplyRadiation: FSR state was not prepared");
    return false;
  }
  std::exception_ptr failure;
  const bool success = state.radiative.Apply(evt, state.random, debug ? &failure : nullptr);
  if (!success) {
    if (failure) {
      try {
        std::rethrow_exception(failure);
      } catch (const std::exception &error) {
        PrintFailure("RadiativeFailure", error.what());
      } catch (...) {
        PrintFailure("RadiativeFailure", "MProcess::ApplyRadiation: unknown radiation failure");
      }
    } else {
      PrintFailure("RadiativeFailure", "MProcess::ApplyRadiation: ISR or FSR event construction returned false");
    }
  }
  return success;
}

// Construct the common central-production event record
bool MProcess::BuildEventRecord(HepMC3::GenEvent &evt) { return record::WriteCentral(state, evt); }

// Amplitude squared with optional eikonal screening
//
double MProcess::GetAmp2(bool include_screening, MEventWeightState &aux) {
  const bool decay_sym_input = state.lts.amplitude.DECAY_SYM;

  double amp2 = 0.0;
  try {
    ProcPtr.DecayStructureFor(state.lts);
    state.lts.amplitude.DECAY_SYM = state.lts.process.root_decay_mode != gra::RootDecayMode::Isolated &&
                                   state.lts.decay_structure.Coherent();
    amp2 = EvaluateAmplitudeBoundary(include_screening, aux);
  } catch (...) {
    state.lts.amplitude.DECAY_SYM = decay_sym_input;
    throw;
  }
  state.lts.amplitude.DECAY_SYM = decay_sym_input;

  return amp2;
}

// Apply one generic exception and status boundary around an amplitude
// evaluation
double MProcess::EvaluateAmplitudeBoundary(bool include_screening, MEventWeightState &aux) {
  evaluation_status     = mg5helas::EvaluationStatus::Success;
  aux.amplitude_failure = false;

  double amp2 = 0.0;
  try {
    PrepareUPCEvent(include_screening);
    amp2 = (state.flat_amplitude == 0) ? ScreenedAmplitudeSquared(include_screening) : GetFlatAmp2(state.lts);
  } catch (const AmplitudeFailure &error) {
    PrintFailure("AmplitudeFailure", error.what());
    amp2                  = 0.0;
    evaluation_status     = mg5helas::EvaluationStatus::AmplitudeFailure;
    aux.RecordAmplitudeFailure(error.what());
    ResetFailedAmplitudeState();
  } catch (const PhaseSpaceFailure &) { throw; }

  // Generated kinematic and amplitude failures are technical, not physics cuts
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) { aux.technical_failure = true; }
  return amp2;
}

// Evaluate an analytic or generated amplitude with shared failure bookkeeping
double MProcess::EvaluateAmplitude() {
  state.lts.proton_good_walker.reset();
  state.lts.screening.BeginNode();
  state.lts.hamp.layout.Clear();
  const double value = EvaluateBareAmplitude();
  if (!mg5helas::EvaluationSucceeded(evaluation_status)) {
    if (evaluation_status == mg5helas::EvaluationStatus::KinematicsFailure) {
      throw PhaseSpaceFailure("MProcess::EvaluateAmplitude: amplitude kinematics failed");
    }
    throw AmplitudeFailure("MProcess::EvaluateAmplitude: amplitude evaluation failed");
  }

  ValidateAmplitude(value);
  return value;
}

// Validate generated amplitudes before normalization or screening accumulation
void MProcess::ValidateAmplitude(double value) const {
  if (!std::isfinite(value) || value < 0.0) {
    throw AmplitudeFailure("MProcess::EvaluateAmplitude: invalid amplitude squared");
  }
  if (!gra::AllFinite(state.lts.hamp)) {
    throw AmplitudeFailure("MProcess::EvaluateAmplitude: non-finite helicity amplitude");
  }
  const bool good_walker = state.lts.hamp.metadata.amplitude_type == ScreeningAmplitudeType::GoodWalker;
  if (good_walker != state.lts.proton_good_walker.has_value()) {
    throw AmplitudeFailure("MProcess::ValidateAmplitude: inconsistent Good Walker amplitude");
  }
  if (good_walker) {
    if (state.lts.proton_good_walker->components.empty()) {
      throw AmplitudeFailure("MProcess::ValidateAmplitude: empty Good Walker amplitude");
    }
    for (const auto &component : state.lts.proton_good_walker->components) {
      if (component.source.size_row() == 0 || component.source.size_col() == 0 || !component.source.IsFinite()) {
        if (component.spin_index) { continue; }  // Diagnostic failures are recorded by MQMetrics
        throw AmplitudeFailure("MProcess::ValidateAmplitude: invalid Good Walker source");
      }
    }
  }
  for (const auto &flow : state.lts.hard_color_flows) {
    if (!gra::AllFinite(flow.amplitudes) ||
        (flow.screened_weight.has_value() && (!std::isfinite(*flow.screened_weight) || *flow.screened_weight < 0.0))) {
      throw AmplitudeFailure("MProcess::EvaluateAmplitude: invalid color-flow amplitude or weight");
    }
  }
}

// Compute "flat/special" matrix element squared |A|^2
// for evaluating the phase space and other special purposes
double MProcess::GetFlatAmp2(const gra::LORENTZSCALAR &input) const {
  double W = 1.0;
  if (state.flat_amplitude >= 1 && state.flat_amplitude <= 3 && !state.model_tune) {
    throw AmplitudeFailure("MProcess::GetFlatAmp2: missing model tune");
  }
  const double B = state.model_tune ? state.model_tune->Flat().B : 0.0;

  // Ansatz: |A|^2 ~ exp(Bt1) exp(Bt2)
  if (state.flat_amplitude == 1) {
    W = pow2(input.s) * form::ExpSlopeWeight(B, input.t1) * form::ExpSlopeWeight(B, input.t2);
  }
  // Ansatz: |A|^2 ~ exp(Bt1) exp(Bt2) / sqrt(shat)
  else if (state.flat_amplitude == 2) {
    W = pow2(input.s) * form::ExpSlopeWeight(B, input.t1) * form::ExpSlopeWeight(B, input.t2) /
        msqrt(input.s_hat);
  }
  // Ansatz: |A|^2 ~ exp(Bt1) exp(Bt2) / shat
  else if (state.flat_amplitude == 3) {
    W = pow2(input.s) * form::ExpSlopeWeight(B, input.t1) * form::ExpSlopeWeight(B, input.t2) /
        input.s_hat;
  }
  // Constant
  else if (state.flat_amplitude == 4) {
    W = 1.0;
  } else {
    // Throw an error, unknown mode
    std::string str =
        "MProcess::GetFlatAmp2: Unknown state.flat_amplitude (|A|^2 is 0 = "
        "off, 1 = "
        "exp{b(t1+t2)}, 2 = exp{b(t1+t2)}/shat^{1/2}, 3 = exp{b(t1+t2)}/shat, "
        "4 = 1.0) "
        ": input was " +
        std::to_string(state.flat_amplitude);
    throw std::invalid_argument(str);
  }
  return W;
}

// Clear event-local color and parton metadata before another sampled point
void MProcess::ResetEventRecordState() noexcept {
  state.nuclear_event.reset();
  state.lts.upc_excitation.reset();
  if (state.nuclear_final.has_value()) { state.nuclear_final->Reset(); }
  for (auto &branch : state.lts.decaytree) { ClearDecayBranchColorFlow(branch); }
  ClearDecayBranchColorFlow(state.lts.decayforward1);
  ClearDecayBranchColorFlow(state.lts.decayforward2);
  state.lts.id1                             = 0;
  state.lts.id2                             = 0;
  state.lts.pdf_xf1                         = 0.0;
  state.lts.pdf_xf2                         = 0.0;
  state.lts.muF                             = 0.0;
  state.lts.muR                             = 0.0;
  state.lts.scalup                          = 0.0;
  state.lts.alphaQCD                        = 0.0;
  state.lts.parton_proposal_weight          = 1.0;
  state.lts.exact_forward_photon_kinematics = false;
  amplitude_event_state_finalized           = false;
  state.multipomeron_chain.clear();
  state.multipomeron_impact_parameter = 0.0;
}

// Initialize all generic event-local state before one phase-space evaluation
void MProcess::BeginEventEvaluation(MEventWeightState &aux) noexcept {
  aux.ResetStatus();
  if (aux.qmetrics) { state.lts.qmetrics.Begin(true); }

  state.radiative.ResetISR(state.lts.pbeam1, state.lts.pbeam2, state.lts.s, state.lts.sqrt_s);
  state.radiative.ResetFSR();

  evaluation_status = mg5helas::EvaluationStatus::Success;
  state.lts.upc_event.reset();
  state.lts.amplitude.AbortCentral();
  state.lts.screening.Clear();
  forward::Reset(state.lts);
  ResetEventRecordState();
  ResetUPCFinalState();
  ResetTransientAmplitudeState();
}

// Clear generic state belonging to a rejected or interrupted event
void MProcess::ResetRejectedEventState() noexcept {
  state.radiative.ResetFSR();
  ResetFailedAmplitudeState();
  forward::Reset(state.lts);
  ResetEventRecordState();
  ResetUPCFinalState();
  state.lts.upc_event.reset();
}

// Prepare common state after one physical phase-space construction
void MProcess::PreparePhaseSpacePoint(bool kinematics_ok, MEventWeightState &aux) {
  aux.kinematics_ok = kinematics_ok;
  if (aux.kinematics_ok && !CEPForwardFragment()) {
    throw PhaseSpaceFailure("MProcess::PreparePhaseSpacePoint: forward excitation failed");
  }
}

// Apply FSR and the final fiducial and veto classification
void MProcess::ClassifyFinalState(MEventWeightState &aux) {
  state.radiative.PrepareFSR(state.lts.decaytree, state.random);
  const auto &central = state.radiative.FiducialTree(state.lts.decaytree);
  aux.fidcuts_ok      = FiducialCuts();
  aux.vetocuts_ok     = aux.fidcuts_ok && state.vetocuts.Pass(state.lts, central);
}

// Apply common fiducial cuts for central and factorized processes
bool MProcess::CommonCuts() const {
  if (!state.fcuts.active) { return true; }
  const auto &central = state.radiative.FiducialTree(state.lts.decaytree);
  if (!UserCut(state.usercuts, state.lts, central)) { return false; }
  if (state.phase_space_class != "P" && !state.fcuts.PassForward(state.lts)) { return false; }
  return state.fcuts.PassCentralSystem(state.lts) && state.fcuts.PassCentralParticles(central) &&
         state.fcuts.PassSelectedParticles(central);
}

// Apply the common central-production fiducial cuts
bool MProcess::FiducialCuts() const { return CommonCuts(); }

// Restore the accepted Born four-momenta and process kinematics
bool MProcess::RestoreBornKinematics() {
  if (state.lts.pfinal_orig.size() < 3 || state.lts.pfinal_orig.size() != state.lts.pfinal.size()) { return false; }
  state.lts.pfinal = state.lts.pfinal_orig;
  return RefreshBornKinematics();
}

// Restore the Born point without allowing cleanup to hide an amplitude failure
bool MProcess::RestoreScreeningBornState() noexcept {
  state.lts.screening.Clear();
  bool restored = false;
  try {
    restored = RestoreBornKinematics();
  } catch (const std::exception &error) {
    PrintFailure("KinematicsFailure", error.what());
    restored = false;
  } catch (...) {
    PrintFailure("KinematicsFailure", "MProcess::RestoreScreeningBornState: unknown restoration failure");
    restored = false;
  }
  if (!restored && state.lts.pfinal_orig.size() == state.lts.pfinal.size()) {
    try {
      state.lts.pfinal = state.lts.pfinal_orig;
    } catch (...) { restored = false; }
  }
  state.lts.screening.Clear();
  return restored;
}

// Roll back loop kinematics and clear every incomplete amplitude output
void MProcess::AbortScreeningLoopState() noexcept {
  RestoreScreeningBornState();
  ResetTransientAmplitudeState();
}

// Clear amplitude outputs after their event-local consumers have finished
void MProcess::ResetTransientAmplitudeState() noexcept {
  state.lts.hamp.clear();
  state.lts.proton_good_walker.reset();
  state.lts.hard_color_flows.clear();
  if (!amplitude_event_state_finalized) { state.lts.upc_event.reset(); }
}

// Clear incomplete amplitude outputs without restoring an unrelated Born point
void MProcess::ResetFailedAmplitudeState() noexcept {
  state.lts.amplitude.AbortCentral();
  if (state.lts.screening.active) {
    AbortScreeningLoopState();
    return;
  }
  state.lts.screening.Clear();
  ResetTransientAmplitudeState();
}

// Validate and normalize orthogonal UPC final-sector weights
void MProcess::SetUPCFinalWeights(const nuclear::FinalWeights &weight) {
  double total = 0.0;
  for (const double value : weight) {
    if (!std::isfinite(value) || value < 0.0) {
      throw AmplitudeFailure("MProcess::SetUPCFinalWeights: invalid final-sector weight");
    }
    total += value;
  }
  ResetUPCFinalState();
  if (!(total > 0.0) || !std::isfinite(total)) { return; }
  for (const auto &i : indices(weight)) { state.upc_final_weight[i] = weight[i] / total; }
}

// Sample the event-local UPC final sector at most once
bool MProcess::SampleUPCFinalState() {
  if (state.lts.upc_model == nullptr || state.upc_final.valid) { return true; }
  double total = 0.0;
  for (const double value : state.upc_final_weight) { total += value; }
  if (!(total > 0.0)) { return false; }
  state.upc_final = nuclear::SampleFinalState(state.upc_final_weight, state.random.U(0.0, 1.0));
  return state.upc_final.valid;
}

// Finalize UPC, color-flow and process-specific event state atomically
bool MProcess::FinalizeEventState() { return SampleUPCFinalState() && FinalizeAmplitudeEventState(); }

// Clear event-local UPC final-sector probabilities and selection
void MProcess::ResetUPCFinalState() noexcept {
  state.upc_final_weight.fill(0.0);
  state.upc_photo = {};
  state.upc_final = {};
}

// Compute that the base process has no additional Born quantities to refresh
bool MProcess::RefreshBornKinematics() { return true; }

}  // namespace gra
