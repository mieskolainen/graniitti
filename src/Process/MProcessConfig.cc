// Process configuration and initialization operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <cstdlib>
#include <iostream>
#include <map>
#include <numeric>
#include <stdexcept>
#include <string>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MIncoherent.h"
#include "Graniitti/Nuclear/MSetup.h"
#include "Graniitti/Photon/MPhotoVM.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Process/MForwardExcitation.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

namespace gra {

using math::msqrt;
using math::pow2;

namespace {

// Require forward particle selections to use the completed HepMC3 record
void ValidateNuclearVeto(const VETOCUT &cuts) {
  if (!cuts.active) { return; }
  for (const auto &cut : cuts.cuts) {
    if (cut.source_forward) {
      throw std::invalid_argument("Nuclear final states require forward particle cuts in icepacks");
    }
  }
}

// Validate every beam combination against the selected process
void ValidateInitialState(const MProcessState &state, const MSubProc &process) {
  const std::array<MParticle, 2> beam      = {state.lts.beam1, state.lts.beam2};
  const bool                     proton1   = nuclear::IsProton(beam[0].pdg);
  const bool                     proton2   = nuclear::IsProton(beam[1].pdg);
  const bool                     ion1      = nuclear::IsNuclearPDG(beam[0].pdg);
  const bool                     ion2      = nuclear::IsNuclearPDG(beam[1].pdg);
  const auto                     collision = nuclear::ClassifyCollision(beam[0].pdg, beam[1].pdg);

  if (!process.SelectedProcess()->Info().beams.Supports(beam[0].pdg, beam[1].pdg)) {
    throw std::invalid_argument("MProcess::PrepareRun: " + process.ISTATE + "[" + process.CHANNEL +
                                "] does not support beam PDGs " + std::to_string(beam[0].pdg) + ", " +
                                std::to_string(beam[1].pdg));
  }
  nuclear::ValidateBeam(beam[0].pdg, beam[0].chargeX3);
  nuclear::ValidateBeam(beam[1].pdg, beam[1].chargeX3);

  const bool nonproton = !proton1 || !proton2;
  if (!nonproton) { return; }

  if (state.excitation != 0 && !((collision == nuclear::CollisionType::EP || collision == nuclear::CollisionType::PA) && state.excitation == 1)) {
    throw std::invalid_argument("MProcess::PrepareRun: NSTARS requires pp beams or NSTARS=1 with ep or pA beams");
  }
  if ((ion1 || ion2) && state.lts.upc_model == nullptr) {
    throw std::invalid_argument("MProcess::PrepareRun: nuclear beams require NUCLEAR steering");
  }
}

// Compute whether either beam has no hadronic overlap profile
bool HasLeptonBeam(const MProcessState &state) {
  return nuclear::IsChargedLepton(state.lts.beam1.pdg) || nuclear::IsChargedLepton(state.lts.beam2.pdg);
}

// Reject forward decay weights that cannot enter the generated cross section
void ValidateBeamFragmentation(const MProcessState &state) {
  if (state.beamfrag == BeamFragType::FewBody && state.nstar_param->fewbody.body_br[1] > 0.0 &&
      !state.nstar_param->fewbody.unweight) {
    throw std::invalid_argument("PARAM_NSTAR.FEWBODY.unweight must be true for three body decays because forward decay weights are not propagated");
  }
}

// Validate screening options which cannot modify the selected amplitude
void ValidateScreeningConfiguration(const MProcessState &state) {
  if (!state.screening) { return; }
  if (HasLeptonBeam(state)) {
    throw std::invalid_argument("MProcess::PrepareRun: LOOPSCREEN requires two hadronic beams");
  }
  if (state.flat_amplitude != 0) {
    throw std::invalid_argument("MProcess::PrepareRun: LOOPSCREEN cannot be combined with @FLATAMP:N for N > 0");
  }
}

// Compute true when a decay branch contains a generic final-state parton alias
bool HasGenericPartonAlias(const MDecayBranch &branch) {
  if (std::abs(branch.p.pdg) == PDG::PDG_hard_jet) { return true; }
  return std::any_of(branch.legs.begin(), branch.legs.end(),
                     [](const MDecayBranch &leg) { return HasGenericPartonAlias(leg); });
}

// Compute true when a decay tree contains a generic final-state parton alias
bool HasGenericPartonAlias(const std::vector<MDecayBranch> &tree) {
  return std::any_of(tree.begin(), tree.end(),
                     [](const MDecayBranch &branch) { return HasGenericPartonAlias(branch); });
}

// Validate charge conservation at every explicit physical decay vertex
void ValidateDecayCharge(const std::vector<MDecayBranch> &tree) {
  for (const auto &branch : tree) {
    if (branch.legs.empty()) { continue; }
    ValidateDecayCharge(branch.legs);
    // Generic parton charges are resolved by the generated process selection
    const bool fixed = std::abs(branch.p.pdg) != PDG::PDG_hard_jet &&
                       std::all_of(branch.legs.begin(), branch.legs.end(), [](const MDecayBranch &leg) {
                         return std::abs(leg.p.pdg) != PDG::PDG_hard_jet;
                       });
    const int charge = std::accumulate(branch.legs.begin(), branch.legs.end(), 0,
                                        [](int sum, const MDecayBranch &leg) { return sum + leg.p.chargeX3; });
    if (fixed && branch.p.chargeX3 != charge) {
      throw std::invalid_argument("MProcess::SetDecayMode: decay of " + branch.p.name +
                                  " violates charge conservation");
    }
  }
}

// Apply process-local external masses to stable outgoing quarks
void ApplyFinalStatePartonMasses(MDecayBranch &branch, const std::map<int, double> &masses) {
  if (branch.legs.empty()) {
    const auto mass = masses.find(std::abs(branch.p.pdg));
    if (mass != masses.end()) { branch.p.mass = mass->second; }
    return;
  }
  for (auto &leg : branch.legs) { ApplyFinalStatePartonMasses(leg, masses); }
}

// Apply process-local external masses to the central decay tree
void ApplyFinalStatePartonMasses(std::vector<MDecayBranch> &tree, const std::map<int, double> &masses) {
  for (auto &branch : tree) { ApplyFinalStatePartonMasses(branch, masses); }
}

// Reject Tensor excitation when no channel-resolved hard source is available
void ValidateTensorExcitationScreening(const MProcessState &state, const MSubProc &process,
                                       const std::size_t channels) {
  if (state.screening && state.excitation > 0 && process.ISTATE == "TP" && channels > 1) {
    throw std::invalid_argument(
        "MProcess::PrepareRun: multichannel Tensor forward excitation needs "
        "a channel-resolved hard source");
  }
}

}  // namespace

// Calculate the process-independent QFT symmetry factor
void MProcess::CalculateSymmetryFactor() {
  state.symmetry_factor = decay::FinalStateSymmetryFactor(state.lts.decaytree, GetISOLATE());
}

// Set decaymode
void MProcess::SetDecayMode(std::string str) {
  // Clear decaytree
  state.lts.decaytree.clear();

  // ------------------------------------------------------------------
  // 0. Check if contains only spaces (or is empty)
  std::string           check   = str;
  std::string::iterator end_pos = std::remove(check.begin(), check.end(), ' ');
  check.erase(end_pos, check.end());

  // Empty final states require quasielastic or factorized luminosity kinematics
  const bool factorized = state.lts.central_phase_space_mode == CentralPhaseSpaceMode::Factorized;
  if (check.empty() && state.process.find("<Q>") == std::string::npos && !factorized) {
    throw std::invalid_argument(
        "MProcess::SetDecayMode: DECAY string input is "
        "empty (can be only for <Q> or "
        "factorized phase space)!");
  }

  // ------------------------------------------------------------------
  // 1. Syntax check, we must have same amount of left { and right } brackets
  unsigned int arrows = 0;
  unsigned int L      = 0;
  unsigned int R      = 0;
  for (const auto &i : indices(str)) {
    if (str[i] == '{') ++L;
    if (str[i] == '}') ++R;
    if (str[i] == '>') ++arrows;
  }
  if (L != R) {
    std::string strerr = "MProces::SetDecayMode: ERROR: Decay tree has " + std::to_string(L) + " left brackets { and " +
                         std::to_string(R) + " right brackets } !";
    throw std::invalid_argument(strerr);
  }
  if (arrows != L) {
    std::string strerr =
        "MProces::SetDecayMode: ERROR: Decay tree syntax with "
        "arrows > not equal to { "
        "} brackets!";
    throw std::invalid_argument(strerr);
  }

  // ===================================================================
  // ** Read decay by recursion **
  state.lts.PDG.TokenizeProcess(str, 0, state.lts.decaytree);
  ValidateDecayCharge(state.lts.decaytree);
  // ===================================================================

  const auto processes = ProcPtr.Processes();
  amplitude::OrderMG5ProcessSyntax(processes, state.lts.decaytree, str);
  const bool has_mg5 = std::any_of(processes.begin(), processes.end(), [](const auto &process) {
    return amplitude::IsGenerated(process.matrix_element_form);
  });
  if (has_mg5 && !HasGenericPartonAlias(state.lts.decaytree) &&
      !ProcPtr.MatchProcess(state.lts.decaytree).has_value()) {
    throw std::invalid_argument("Unsupported final state '" + str +
                                "': no analytic or one-to-one MG5 PDG ordering matches");
  }

  parton::ConfigureFinalStateMasses(state.lts, ProcPtr);
  ApplyFinalStatePartonMasses(state.lts.decaytree, state.lts.final_state_parton_masses);

  if (state.lts.final_state_partons.IsConfigured() && !HasGenericPartonAlias(state.lts.decaytree)) {
    throw std::invalid_argument(
        "MProcess::SetDecayMode: @j requires at least "
        "one final-state j or ~j alias");
  }
  parton::ConfigureFinalStateProposal(state.lts, ProcPtr);

  // Save decaymode string
  state.decay_mode = str;

  // Calculate symmetry factor
  CalculateSymmetryFactor();

  // Preserve generic aliases so each event can propose an on-shell physical
  // mode
  if (HasGenericPartonAlias(state.lts.decaytree)) {
    state.lts.generic_parton_decaytree = state.lts.decaytree;
  } else {
    state.lts.generic_parton_decaytree.clear();
  }
  state.lts.parton_proposal_weight = 1.0;

  // Remind the user only for analytic factorized decay amplitudes
  const bool factorized_decay_warning = state.lts.decaytree.size() > 2 &&
                                        state.lts.process.root_decay_mode != RootDecayMode::Isolated &&
                                        ProcPtr.UsesJWHelicityAlgebra();
  if (factorized_decay_warning) {
    gra::aux::PrintWarning();
    std::cout << "Reminder: Resonance decay |matrix element 1->K|^2 is "
                 "non-factorizable from the phase space for K = " +
                     std::to_string(state.lts.decaytree.size()) + " > 2!"
              << std::endl;
    std::cout << "Use &> arrow for fully separating 2->3 [+] 1->K" << std::endl;
    std::cout << std::endl;
  }
}

// Load process-local immutable helicity steering cards
void MProcess::SetHelicityConfig(MModelTunePtr tune) { helicity_config.Load(tune); }

// Construct the process-independent amplitude initialization context
MProcessSetup MProcess::CreateProcessSetup() {
  if (state.soft_model == nullptr) { throw std::logic_error("MProcess::CreateProcessSetup: SOFT model is not set"); }
  if (!helicity_config.IsLoaded()) {
    throw std::logic_error("MProcess::CreateProcessSetup: helicity cards are not loaded");
  }

  const auto structure = [this](const MParticle &particle, const std::vector<MParticle> &legs, bool production_mode,
                                bool strict_mode, const std::string &verbose_label, bool verbose_output,
                                gra::spin::VertexContext context) {
    return ProcessHelicityStructure(particle, legs, production_mode, strict_mode, verbose_label, verbose_output,
                                    context);
  };
  const auto tree = [this](MDecayBranch &branch, bool production_mode, bool strict_mode) {
    ProcessHelicityTree(branch, production_mode, strict_mode);
  };
  const auto pole_operator = [this](ReggeProductionModel model, const MParticle &particle,
                                    const std::vector<MParticle> &legs, gra::spin::VertexContext context) {
    return ProcessPoleOperatorStructure(model, particle, legs, context);
  };

  return MProcessSetup(state.lts, state.model_tune, ProcPtr.ISTATE, ProcPtr.CHANNEL, ProcPtr.SelectedProcess()->Info(), state.flat_mass2, GetISOLATE(),
                       state.excitation, state.spingen_user, state.spindec_user, state.mmax_user, state.symmetry_factor,
                       structure, tree, pole_operator);
}

// Initialize the selected physical amplitude with worker-local services
void MProcess::InitializeProcessAmplitude() {
  ProcPtr.BindRandom(state.random);
  if (!state.qmetrics_user) {
    state.lts.process.QMETRICS = state.model_tune->General("PARAM_SPIN").at("QMETRICS").get<bool>();
  }
  MProcessSetup setup = CreateProcessSetup();
  ProcPtr.InitializeAmplitude(setup);
  if (state.lts.process.QMETRICS &&
      (!state.lts.qmetrics.Configured() || state.flat_amplitude != 0 ||
       (state.phase_space_class != "F" && state.phase_space_class != "C"))) {
    const auto &numerics = state.model_tune->Numerics("NUMERICS_SPIN");
    state.lts.qmetrics.Configure({}, numerics.at("metric_groups"), numerics.at("metric_max_dimension"),
                                     "Spin metrics are unavailable for this process");
  }

  const ScreeningMetadata &screening = state.lts.hamp.metadata;
  if (screening.spin_basis == ScreeningSpinBasis::Unset || !std::isfinite(screening.amplitude_normalization) ||
      screening.amplitude_normalization <= 0.0 || screening.spin_rows == 0) {
    throw std::invalid_argument("MProcess::InitializeProcessAmplitude: invalid screening metadata");
  }
}

// No additional eikonal momentum-table range is required by default
double MProcess::EikonalMaxKT2() const { return 0.0; }

// Initialize the complete process runtime before worker copies are made
void MProcess::PrepareRun() {
  ValidateInitialState(state, ProcPtr);
  ValidateBeamFragmentation(state);
  ValidateScreeningConfiguration(state);
  InitializeProcessAmplitude();
  FinalizeProcessConfiguration();
  ValidateTensorExcitationScreening(state, ProcPtr, state.soft_model->GoodWalker().ChannelCount());

  const bool elastic_cni         = ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "EL";
  const bool soft_nondiffractive = ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "ND";
  const bool nuclear_upc         = state.lts.upc_model != nullptr;
  const bool lepton_beam         = HasLeptonBeam(state);
  if (!nuclear_upc && !lepton_beam && (state.screening || elastic_cni || soft_nondiffractive) &&
      !eikonal.IsInitialized()) {
    eikonal.S3Constructor(state.lts.s, GetInitialState(), false, 0, 0, EikonalMaxKT2());
  }
  if (state.screening && !nuclear_upc && !lepton_beam) { (void)eikonal.GetLoopConst(state.lts.s); }
  if (elastic_cni) { eikonal.InitializeElasticCNI(state.gcuts.q_t_abs_min, state.gcuts.q_t_abs_max); }
}

// Initialize the complete process runtime with an external eikonal
void MProcess::PrepareRun(const MEikonal &eikonal) {
  ValidateInitialState(state, ProcPtr);
  ValidateBeamFragmentation(state);
  ValidateScreeningConfiguration(state);
  if (state.lts.upc_model != nullptr || HasLeptonBeam(state)) {
    throw std::invalid_argument(
        "MProcess::PrepareRun: external pp eikonal is invalid for pA, AA, "
        "ep, eA or ee beams");
  }
  InitializeProcessAmplitude();
  FinalizeProcessConfiguration();
  SetEikonal(eikonal);
  ValidateTensorExcitationScreening(state, ProcPtr, this->eikonal.GetChannelCount());

  if (state.screening) { (void)this->eikonal.GetLoopConst(state.lts.s); }
  if (ProcPtr.ISTATE == "X" && ProcPtr.CHANNEL == "EL") {
    this->eikonal.InitializeElasticCNI(state.gcuts.q_t_abs_min, state.gcuts.q_t_abs_max);
  }
}

// Construct one helicity vertex from the process-local steering cards
HELMatrix MProcess::ProcessHelicityStructure(const MParticle &particle, const std::vector<MParticle> &legs,
                                             bool production_mode, bool strict_mode, const std::string &verbose_label,
                                             bool verbose_output, gra::spin::VertexContext context) const {
  return helicity_config.ProcessHelicityStructure(particle, legs, state.lts, ProcPtr, state.model_tune,
                                                  state.decay_mode, production_mode, strict_mode, verbose_label,
                                                  verbose_output, context);
}

// Construct helicity vertices recursively through one decay branch
void MProcess::ProcessHelicityTree(MDecayBranch &branch, bool production_mode, bool strict_mode) {
  helicity_config.ProcessHelicityTree(branch, state.lts, ProcPtr, state.model_tune, state.decay_mode, production_mode,
                                      strict_mode);
}

// Construct one canonical MP or XP pole operator from its continuum card
spin::PoleLS MProcess::ProcessPoleOperatorStructure(ReggeProductionModel model, const MParticle &particle,
                                                               const std::vector<MParticle> &legs,
                                                               gra::spin::VertexContext      context) const {
  return helicity_config.ProcessPoleOperatorStructure(model, particle, legs, state.lts, ProcPtr, state.model_tune,
                                                      context);
}

// Set isolated phase-space sampling
void MProcess::SetISOLATE(bool value) { state.lts.PS_active = !value; }

// Compute whether isolated phase-space sampling is active
bool MProcess::GetISOLATE() const { return !state.lts.PS_active; }

// Set root decay semantics selected by the process arrow
void MProcess::SetRootDecayMode(RootDecayMode mode) { state.lts.process.root_decay_mode = mode; }

// Set flat mass-squared sampling in decay cascades
void MProcess::SetFLATMASS2(bool value) {
  aux::PrintNotice();
  std::cout << rang::fg::red
            << "MProcess::SetFLATMASS2: Set flat in mass^2 sampling in decay trees: " << rang::fg::reset << std::endl;
  state.flat_mass2      = value;
  state.flat_mass2_user = true;
}

// Check whether flat mass-squared sampling is active
bool MProcess::GetFLATMASS2() const { return state.flat_mass2; }

// Set the offshell sampling window in decay cascades
void MProcess::SetOFFSHELL(double value) {
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument("MProcess::SetOFFSHELL: input must be finite and nonnegative");
  }
  aux::PrintNotice();
  std::cout << rang::fg::red << "MProcess::SetOFFSHELL: Set number of decay widths in decay trees: " << value
            << rang::fg::reset << std::endl;
  state.offshell_widths      = value;
  state.offshell_widths_user = true;
}

// Compute the offshell sampling window
double MProcess::GetOFFSHELL() const { return state.offshell_widths; }

// Set the default offshell sampling window
void MProcess::SetDefaultOFFSHELL(double value) {
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument("MProcess::SetDefaultOFFSHELL: input must be finite and nonnegative");
  }
  state.offshell_widths = value;
}

// Set the minimum width for numerical offshell sampling
void MProcess::SetWIDTHMIN(double value) {
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument("MProcess::SetWIDTHMIN: input must be finite and nonnegative");
  }
  state.width_min = value;
}

// Compute the minimum width for numerical offshell sampling
double MProcess::GetWIDTHMIN() const { return state.width_min; }

// Set generation spin correlations
void MProcess::SetSPINGEN(bool value) {
  std::cout << rang::fg::red
            << "MProcess::SetSPINGEN: Set generation 2->1 spin correlations: " << (value ? "true" : "false")
            << rang::fg::reset << std::endl;
  state.lts.process.SPINGEN = value;
  state.spingen_user        = true;
}

// Set decay spin correlations
void MProcess::SetSPINDEC(bool value) {
  std::cout << rang::fg::red << "MProcess::SetSPINDEC: Set decay 1->2 spin correlations: " << (value ? "true" : "false")
            << rang::fg::reset << std::endl;
  state.lts.process.SPINDEC = value;
  state.spindec_user        = true;
}

// Enable initialization references and integrated spin metrics
void MProcess::SetQMetrics(const bool value) {
  state.lts.process.QMETRICS = value;
  state.qmetrics_user = true;
}

// Set the MP central spin basis frame
void MProcess::SetMPFrame(const std::string &frame) {
  std::cout << rang::fg::red << "MProcess::SetMPFrame: Set MP central spin basis frame: " << frame << rang::fg::reset
            << std::endl;
  state.lts.process.MP_FRAME = frame;
}

// Set the maximum retained analytic Regge helicity
void MProcess::SetMMAX(int value) {
  std::cout << rang::fg::red << "MProcess::SetMMAX: Set maximum analytic Regge helicity: " << value << rang::fg::reset
            << std::endl;
  state.lts.process.MMAX = value;
  state.mmax_user        = true;
}

// Compute the configured initial state
std::vector<MParticle> MProcess::GetInitialState() const { return {state.lts.beam1, state.lts.beam2}; }

// Compute the phase-space dimension
unsigned int MProcess::GetdLIPSDim() const {
  return ProcPtr.LIPSDIM + state.radiative.ISRDim(state.lts.beam1, state.lts.beam2);
}

// Set eikonal screening
void MProcess::SetScreening(bool value) { state.screening = value; }

// Check whether eikonal screening is active
bool MProcess::GetScreening() const { return state.screening; }

// Set a validated external eikonal
void MProcess::SetEikonal(const MEikonal &eikonal) {
  if (state.soft_model == nullptr || eikonal.SoftModelHandle() != state.soft_model) {
    throw std::invalid_argument("MProcess::SetEikonal: SOFT model mismatch");
  }
  if (!eikonal.IsInitialized()) {
    throw std::invalid_argument("MProcess::SetEikonal: external eikonal is not initialized");
  }
  if (!std::isfinite(state.lts.s) || state.lts.s <= 0.0 || state.lts.beam1.pdg == 0 || state.lts.beam2.pdg == 0) {
    throw std::invalid_argument("MProcess::SetEikonal: process initial state is not initialized");
  }
  if (!std::is_eq(eikonal.InitializedMandelstamS() <=> state.lts.s)) {
    throw std::invalid_argument("MProcess::SetEikonal: exact Mandelstam s mismatch");
  }

  const std::array<MParticle, 2> beams         = {state.lts.beam1, state.lts.beam2};
  const std::vector<MParticle>  &eikonal_beams = eikonal.InitialState();
  if (eikonal_beams.size() != beams.size()) {
    throw std::invalid_argument("MProcess::SetEikonal: external eikonal needs exactly two beams");
  }
  for (const auto &beam : indices(beams)) {
    if (eikonal_beams[beam].pdg != beams[beam].pdg || !std::is_eq(eikonal_beams[beam].mass <=> beams[beam].mass) ||
        eikonal_beams[beam].chargeX3 != beams[beam].chargeX3) {
      throw std::invalid_argument("MProcess::SetEikonal: ordered beam " + std::to_string(beam + 1) + " mismatch");
    }
  }
  this->eikonal = eikonal;
}

// Compute the coupled eikonal
const MEikonal &MProcess::GetEikonal() const { return eikonal; }

// Set the immutable SOFT model shared by amplitudes and screening
void MProcess::SetSoftModel(SoftModelPtr model) {
  if (model == nullptr) { throw std::invalid_argument("MProcess::SetSoftModel: null SOFT model"); }
  const int extra_loop_kt  = eikonal.Numerics.extra_NumberLoopKT;
  const int extra_loop_phi = eikonal.Numerics.extra_NumberLoopPHI;
  state.soft_model         = std::move(model);
  eikonal.Numerics.SetLoopDiscretization(extra_loop_kt, extra_loop_phi);
}

// Access the immutable SOFT model
const SoftModelPtr &MProcess::GetSoftModel() const noexcept { return state.soft_model; }

// Set the immutable complete model tune shared by amplitudes and screening
void MProcess::SetModelTune(MModelTunePtr tune) {
  if (tune == nullptr) { throw std::invalid_argument("MProcess::SetModelTune: null model tune"); }
  state.model_tune      = std::move(tune);
  state.lts.model_cache = std::make_shared<MModelCache>(state.model_tune);
  SetSoftModel(state.model_tune->Soft());
  const int extra_loop_kt  = eikonal.Numerics.extra_NumberLoopKT;
  const int extra_loop_phi = eikonal.Numerics.extra_NumberLoopPHI;
  eikonal                  = MEikonal(state.model_tune);
  eikonal.Numerics.SetLoopDiscretization(extra_loop_kt, extra_loop_phi);
  state.nstar_param = GetNstarParam(*state.lts.model_cache);
}

// Access the immutable complete model tune
const MModelTunePtr &MProcess::GetModelTune() const noexcept { return state.model_tune; }

// Set the process PDF set
void MProcess::SetLHAPDF(const std::string &name) {
  std::cout << "MProcess::SetLHAPDF: " << name << std::endl;
  state.lts.LHAPDFSET    = name;
  state.lts.GlobalPdfPtr = nullptr;
}

// Set generation cuts
void MProcess::SetGenCuts(const GENCUT &cuts) { state.gcuts = cuts; }

// Set fiducial cuts
void MProcess::SetFidCuts(const FIDCUT &cuts) { state.fcuts = cuts; }

// Set the custom fiducial-cut selector
void MProcess::SetUserCuts(std::int64_t id) { state.usercuts = id; }

// Set veto cuts
void MProcess::SetVetoCuts(const VETOCUT &cuts) {
  if (state.nuclear_final.has_value()) { ValidateNuclearVeto(cuts); }
  state.vetocuts = cuts;
}

// Select and validate beam remnant fragmentation before event generation
void MProcess::SetBeamFrag(const std::string &mode) {
  const std::map<std::string, BeamFragType> modes = {
      {"none", BeamFragType::None}, {"fewbody", BeamFragType::FewBody},
      {"cylinder", BeamFragType::Cylinder}, {"diquark", BeamFragType::Diquark}};
  const auto selected = modes.find(mode);
  if (selected == modes.end()) {
    throw std::invalid_argument("SCATTERING.BEAMFRAG must be none, fewbody, cylinder or diquark, got " + mode);
  }
  state.beamfrag = selected->second;
}

// Set forward excitation
void MProcess::SetExcitation(int value) {
  if (value < 0 || value > 2) { throw std::invalid_argument("MProcess::SetExcitation: input must be 0, 1 or 2"); }
  state.excitation = value;
  if (value > 0) {
    aux::PrintWarning();
    std::cout << rang::fg::red << "MProcess::SetExcitation: proton excitation is under construction" << rang::fg::reset
              << std::endl;
  }
}

// Set the flat matrix-element debug mode
void MProcess::SetFLATAMP(int value) {
  if (value < 0 || value > 4) { throw std::invalid_argument("MProcess::SetFLATAMP: input must be between 0 and 4"); }
  state.flat_amplitude = value;
}

// Set input resonances
void MProcess::SetResonances(const std::map<std::string, PARAM_RES> &resonances) {
  state.lts.process.RESONANCES = resonances;
}

// Compute input resonances
const std::map<std::string, PARAM_RES> &MProcess::GetResonances() const { return state.lts.process.RESONANCES; }

// Print the shared process and active subprocess setup
void MProcess::PrintSetup() const {
  std::cout << std::endl;
  std::cout << rang::style::bold << "Process setup:" << rang::style::reset << std::endl << std::endl;
  std::cout << "- Random seed:      " << state.random.GetSeed() << std::endl;
  std::cout << "- Initial state:    [" << state.lts.beam1.name << " " << state.lts.beam2.name << "]" << std::endl;
  const double a1 = static_cast<double>(nuclear::BeamMassNumber(state.lts.beam1.pdg));
  const double a2 = static_cast<double>(nuclear::BeamMassNumber(state.lts.beam2.pdg));
  printf("- Energies:         [%0.1f %0.1f] GeV\n", state.lts.pbeam1.E() / a1, state.lts.pbeam2.E() / a2);
  printf("- Constituent CMS:   %0.1f GeV\n", (state.lts.pbeam1 / a1 + state.lts.pbeam2 / a2).M());
  std::cout << "- Process:          " << state.process << rang::fg::green << " [";

  for (const auto &str : ProcPtr.GetProcessDescriptor(state.process)) { std::cout << str << " | "; }
  std::cout << "]" << rang::fg::reset << std::endl << std::endl;

  // Subprocess setup
  std::cout << rang::style::bold << "Subprocess setup:" << rang::style::reset << std::endl << std::endl;
  const bool proton_pair = nuclear::IsProton(state.lts.beam1.pdg) && nuclear::IsProton(state.lts.beam2.pdg);
  if (proton_pair) { std::cout << "- Pomeron loop screening:  " << std::boolalpha << state.screening << std::endl; }
  if (state.lts.upc_model != nullptr) { nuclear::PrintUPCSetup(*state.lts.upc_model, std::cout); }

  // All other than inclusive processes
  if (ProcPtr.ISTATE != "X") {
    std::cout << "- Final state:             " << state.decay_mode << std::endl;
    if (proton_pair) {
      std::cout << "- Proton N* excitation:    " << std::boolalpha << state.excitation;

      if (state.excitation == 0) { std::cout << rang::fg::green << "  <elastic>" << rang::fg::reset << std::endl; }
      if (state.excitation == 1) { std::cout << rang::fg::green << "  <single>" << rang::fg::reset << std::endl; }
      if (state.excitation == 2) { std::cout << rang::fg::green << "  <double>" << rang::fg::reset << std::endl; }
    }
  }
  if (state.flat_amplitude != 0) { std::cout << "- Flat amplitude mode:     " << state.flat_amplitude << std::endl; }

  std::cout << std::endl << std::endl;
}


// Set beam particle 4-vectors with nuclear energies steered per nucleon
void MProcess::SetInitialState(const std::vector<std::string> &beam, const std::vector<double> &energy) {
  if (beam.size() != 2) { throw std::invalid_argument("MProcess::SetInitialState: Input BEAM vector not dim 2!"); }
  if (energy.size() != 2) { throw std::invalid_argument("MProcess::SetInitialState: Input ENERGY vector not dim 2!"); }
  for (const double value : energy) {
    if (!std::isfinite(value) || value < 0.0) {
      throw std::invalid_argument("MProcess::SetInitialState: beam energies must be finite and nonnegative");
    }
  }

  state.lts.beam1 = state.lts.PDG.FindByPDGName(beam[0]);
  state.lts.beam2 = state.lts.PDG.FindByPDGName(beam[1]);

  const double input1 = nuclear::FullBeamEnergy(state.lts.beam1.pdg, energy[0]);
  const double input2 = nuclear::FullBeamEnergy(state.lts.beam2.pdg, energy[1]);

  // Preserve the requested on-shell energies, including a stationary target
  SetBeamEnergies(input1, input2);
}

// Set the heavy-ion UPC model after both beam particles are known
void MProcess::SetNuclear(const nuclear::UPCParam &param, std::unique_ptr<nuclear::MReaction> model) {
  const std::array<MParticle, 2> beam       = {state.lts.beam1, state.lts.beam2};
  const std::array<M4Vec, 2>     momentum   = {state.lts.pbeam1, state.lts.pbeam2};
  nuclear::UPCParam              configured = param;
  configured.glauber.profile = {};
  configured.survival_eikonal.clear();
  configured.s_nn = 0.0;
  if (model) { configured.reaction = nuclear::ReactionType::External; }
  if (state.model_tune == nullptr) { throw std::logic_error("MProcess::SetNuclear: model tune is not set"); }
  configured.structure_param = state.model_tune->Structure();
  const auto is_hadron       = [](const MParticle &particle) {
    return nuclear::IsNuclearPDG(particle.pdg) || std::abs(particle.pdg) == PDG::PDG_p;
  };
  const bool hadronic_pair = is_hadron(beam[0]) && is_hadron(beam[1]);
  if (hadronic_pair && state.screening) {
    const std::array<double, 2> mass_number = {static_cast<double>(nuclear::BeamMassNumber(beam[0].pdg)),
                                               static_cast<double>(nuclear::BeamMassNumber(beam[1].pdg))};
    // A proton has A = 1 while each ion momentum is reduced per nucleon
    const M4Vec  constituent1 = momentum[0] / mass_number[0];
    const M4Vec  constituent2 = momentum[1] / mass_number[1];
    const double s_nn         = (constituent1 + constituent2).M2();
    if (!std::isfinite(s_nn) || !(s_nn > 4.0 * PDG::mp * PDG::mp)) {
      throw std::invalid_argument("MProcess::SetNuclear: invalid nucleon-nucleon energy");
    }
    const MParticle proton = state.lts.PDG.FindByPDG(PDG::PDG_p);
    const auto tune = configured.survival == nuclear::SurvivalType::Optical
                          ? state.model_tune : nuclear::GGCFTune(state.model_tune, configured.ggcf);
    configured.survival_eikonal = tune->Soft()->ActiveModel();
    configured.s_nn = s_nn;
    MEikonal        nucleon_eikonal(tune);
    nucleon_eikonal.S3Constructor(s_nn, {proton, proton}, true, static_cast<int>(configured.glauber.profile_nodes), 2);
    configured.glauber.profile = nuclear::BuildNNProfile(nucleon_eikonal, configured.sigma_nn);
  } else if (configured.sigma_nn.has_value()) {
    throw std::invalid_argument("MProcess::SetNuclear: NUCLEAR::sigma_NN requires active hadronic survival");
  }
  nuclear::UPCMode mode;
  mode.screening = state.screening;
  mode.emission = {false, false};
  mode.target   = {false, false};
  mode.current  = {false, false};
  if (ProcPtr.ISTATE == "yy") {
    for (const auto &leg : indices(beam)) {
      mode.emission[leg] = nuclear::IsNuclearPDG(beam[leg].pdg);
      mode.current[leg]  = mode.emission[leg] && configured.emission[leg] != nuclear::CoherenceType::Coherent;
    }
  } else if (ProcPtr.ISTATE == "ygg" || ProcPtr.ISTATE == "GP" || ProcPtr.ISTATE == "TP") {
    for (const auto &leg : indices(configured.emission)) {
      if (!nuclear::IsNuclearPDG(beam[leg].pdg)) { continue; }
      const std::size_t other = 1 - leg;
      mode.emission[leg]      = is_hadron(beam[other]);
      mode.target[leg]        = nuclear::IsEPAEmitter(beam[other].pdg);
      mode.current[leg]       = (mode.emission[leg] && configured.emission[leg] != nuclear::CoherenceType::Coherent) ||
                          (mode.target[leg] && configured.target[leg] != nuclear::CoherenceType::Coherent);
    }
  }
  if (configured.photo_model == nuclear::PhotoModel::LTA && (mode.target[0] || mode.target[1])) {
    // The input-scale gluon response supports charmonium, without nuclear DGLAP evolution
    const bool jmrt  = ProcPtr.ISTATE == "ygg" && (ProcPtr.CHANNEL == "jpsi" || ProcPtr.CHANNEL == "psi(2S)");
    bool       charm = jmrt || (ProcPtr.ISTATE == "GP" && !state.lts.process.RESONANCES.empty());
    for (const auto &[name, res] : state.lts.process.RESONANCES) { charm &= res.p.pdg == 443 || res.p.pdg == 100443; }
    if (!charm || state.phase_space_class != "F") {
      throw std::invalid_argument("MProcess::SetNuclear: input-scale LTA requires charmonium in <F> phase space");
    }
    const double scale =
        jmrt ? pow2(ReadPhotoVMNumerics(*state.model_tune)->HeavyQuarkMass(4)) : pow2(state.gcuts.M_max / 2.0);
    const double mass2 = jmrt ? 4.0 * scale : pow2(state.gcuts.M_min);
    auto         cuts  = state.gcuts;
    SetTechnicalBoundaries(cuts, state.excitation);
    const double mt = std::hypot(cuts.M_max, 2.0 * cuts.forward_pt_max);
    for (const auto &leg : indices(mode.target)) {
      if (!mode.target[leg]) { continue; }
      const auto  &shadow   = configured.photo[leg].shadow;
      const auto  &target   = momentum[leg];
      const auto  &source   = momentum[1 - leg];
      const double sign     = leg == 0 ? -1.0 : 1.0;
      const double rapidity = leg == 0 ? -state.gcuts.Y_min : state.gcuts.Y_max;
      const double plus = target.E() + sign * target.Pz(), minus = target.E() - sign * target.Pz();
      const double w2 = (target.M2() / nuclear::BeamMassNumber(beam[leg].pdg) +
                         (mt * std::exp(rapidity) + plus) * minus + (source.E() - sign * source.Pz()) * plus) /
                        nuclear::BeamMassNumber(beam[leg].pdg);
      if (!(mass2 > 0.0) || scale > shadow.q2_max || mass2 / w2 < shadow.x_min) {
        throw std::invalid_argument(
            "MProcess::SetNuclear: mass, rapidity or transverse cuts exceed the LTA input domain");
      }
    }
  }
  auto upc             = nuclear::BuildUPC(beam, momentum, state.excitation, configured, mode);
  state.upc_param      = upc->Param();
  state.upc_configured = true;
  state.lts.upc_model  = std::move(upc);
  state.nuclear_final.reset();
  if (configured.reaction.has_value()) {
    ValidateNuclearVeto(state.vetocuts);
    if (state.phase_space_class != "F" && state.phase_space_class != "C") {
      throw std::invalid_argument("MProcess::SetNuclear: nuclear final states require <F> or <C> phase space");
    }
    if (model) {
      state.nuclear_final.emplace(std::move(model), configured.neutron);
    } else if (configured.reaction == nuclear::ReactionType::Internal) {
      state.nuclear_final.emplace(std::make_unique<nuclear::MIncoherent>(state.lts.upc_model), configured.neutron);
    } else {
      state.nuclear_final.emplace(configured.final_library, nlohmann::json::parse(configured.final_config), state.upc_param);
    }
  }
}

// Set the same nuclear reaction interface for an embedded external backend
void MProcess::SetNuclearFinalState(std::unique_ptr<nuclear::MReaction> model) {
  if (!HasNuclear() || (state.phase_space_class != "F" && state.phase_space_class != "C")) {
    throw std::invalid_argument("MProcess::SetNuclearFinalState: requires nuclear <F> or <C> production");
  }
  if (!model) { throw std::invalid_argument("MProcess::SetNuclearFinalState: null model"); }
  SetNuclear(state.upc_param, std::move(model));
}

// Compute whether a heavy-ion UPC model is active
bool MProcess::HasNuclear() const noexcept { return state.lts.upc_model != nullptr; }

// Set finite on-shell beam momenta after the beam particles are specified
void MProcess::SetBeamEnergies(double E1, double E2) {
  if (state.lts.beam1.pdg == 0 || state.lts.beam2.pdg == 0) {
    throw std::invalid_argument("MProcess::SetBeamEnergies: Beam PDG particles not set yet!");
  }

  const double m1 = state.lts.beam1.mass;
  const double m2 = state.lts.beam2.mass;
  if (!std::isfinite(m1) || !std::isfinite(m2) || m1 < 0.0 || m2 < 0.0 ||
      !std::isfinite(E1) || !std::isfinite(E2) || !(E1 > 0.0) || !(E2 > 0.0) || E1 < m1 || E2 < m2) {
    throw std::invalid_argument("MProcess::SetBeamEnergies: finite beam energies must be at least their masses");
  }

  // Form the new state before replacing the nominal ISR beams
  const M4Vec beam1(0, 0, std::sqrt((E1 - m1) * (E1 + m1)), E1);
  const M4Vec beam2(0, 0, -std::sqrt((E2 - m2) * (E2 + m2)), E2);
  const double s = (beam1 + beam2).M2();
  if (!std::isfinite(s) || !(s > 0.0) || std::sqrt(s) < m1 + m2) {
    throw std::invalid_argument("MProcess::SetBeamEnergies: invalid invariant collision energy");
  }
  state.lts.pbeam1 = beam1;
  state.lts.pbeam2 = beam2;
  state.lts.s = s;
  state.lts.sqrt_s = std::sqrt(s);
  state.radiative.SetBeams(beam1, beam2);

  printf(
      "MProcess::SetBeamEnergies: beam: [%s, %s], steering energy = [%0.1f, "
      "%0.1f] GeV\n",
      state.lts.beam1.name.c_str(), state.lts.beam2.name.c_str(), nuclear::SteeringBeamEnergy(state.lts.beam1.pdg, E1),
      nuclear::SteeringBeamEnergy(state.lts.beam2.pdg, E2));
}

// Set the QED initial and final state radiation configuration
void MProcess::SetRadiative(const radiative::Config &config) {
  state.radiative.Configure(config);
  state.radiative.ValidateApplicability(state.lts.beam1, state.lts.beam2, state.lts.decaytree);
  state.radiative.SetBeams(state.lts.pbeam1, state.lts.pbeam2);
}


// Validate fiducial cuts against the configured phase-space class
void MProcess::ValidateFiducialCuts(const gra::FIDCUT &cuts, std::int64_t usercuts) const {
  if (!cuts.active && usercuts != 0) {
    throw std::invalid_argument("MProcess::ValidateFiducialCuts: USERCUTS requires active FIDCUTS");
  }
  if (!cuts.active) { return; }
  if (!IsKnownUserCut(usercuts)) {
    throw std::invalid_argument("MProcess::ValidateFiducialCuts: unknown USERCUTS ID " + std::to_string(usercuts));
  }

  for (const auto &cut : cuts.pdg_cuts) {
    if (cut.pdg.size() != cut.pdg_abs.size()) {
      throw std::invalid_argument("MProcess::ValidateFiducialCuts: malformed PDG selector");
    }
  }

  const bool central_cuts = cuts.HasCentralParticleCuts() || cuts.HasCentralSystemCuts() || !cuts.pdg_cuts.empty();
  if (state.phase_space_class == "Q" && central_cuts) {
    throw std::invalid_argument(
        "MProcess::ValidateFiducialCuts: <Q> does not "
        "support CENTRAL fiducial cuts");
  }
  if (state.phase_space_class == "P" && cuts.HasForwardCuts()) {
    throw std::invalid_argument(
        "MProcess::ValidateFiducialCuts: <P> does not "
        "support FORWARD fiducial cuts");
  }
  if (!cuts.forward_M_active) { return; }

  if (state.phase_space_class == "Q" && ProcPtr.CHANNEL != "SD" && ProcPtr.CHANNEL != "DD") {
    throw std::invalid_argument(
        "MProcess::ValidateFiducialCuts: FORWARD M "
        "requires <Q> SD or DD excitation");
  }
  if (state.phase_space_class != "Q" && state.excitation == 0) {
    throw std::invalid_argument("MProcess::ValidateFiducialCuts: FORWARD M requires NSTARS excitation");
  }
}


// This is called last by the initialization routines
// as the last step before event generation
void MProcess::SetTechnicalBoundaries(gra::GENCUT &cuts, unsigned int excitation) {
  if (excitation != 0) { (void)forward::MassBounds(state, cuts); }
  if (cuts.forward_pt_min < 0.0) {  // Not set yet by the USER
    cuts.forward_pt_min = 0.0;
  }

  if (cuts.forward_pt_max < 0.0) {  // Not set yet by the USER

    if (excitation == 0) {  // Elastic forward protons, default values
      cuts.forward_pt_max = 2.5;
    }

    else if (excitation == 1) {  // Single excitation
      cuts.forward_pt_max = 50.0;
    } else if (excitation == 2) {  // Double excitation
      cuts.forward_pt_max = 100.0;
    }
  }
}


// Set process
void MProcess::SetProcess(std::string &str, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune) {
  // SET IT HERE!
  state.process = str;

  // Call this always first
  if (tune != nullptr) {
    state.lts.PDG.ReadParticleData(MPDG::DataFile(), tune->Directory());
    SetModelTune(std::move(tune));
  } else {
    state.lts.PDG.ReadParticleData();
  }

  // Read the process-local outgoing parton species without changing PDF
  // flavours
  state.lts.final_state_partons = MPartonFlavour();
  std::size_t selector_count    = 0;
  for (const auto &command : syntax) {
    if (command.id != "j") { continue; }
    ++selector_count;
    if (selector_count > 1) { throw std::invalid_argument("@Syntax error: duplicate @j final-state selector"); }
    state.lts.final_state_partons.Configure(command.values);
  }

  state.lts.particle_mass_overrides.clear();
  state.lts.particle_width_overrides.clear();

  // @SYNTAX Read and set new PDG input
  for (const auto &i : indices(syntax)) {
    if (syntax[i].id == "PDG") {
      // Take target string
      std::string pdg_;
      if (syntax[i].target.size() == 1) {
        pdg_ = syntax[i].target[0];
      } else if (syntax[i].target.size() == 0) {
        throw std::invalid_argument("@Syntax error: invalid PDG[] without any target []");
      } else if (syntax[i].target.size() > 1) {
        throw std::invalid_argument("@Syntax error: invalid PDG[] with multiple targets inside []");
      }
      // Conversion to integer
      int pdg = 0;

      try {
        pdg = aux::ParseInt(pdg_, "@PDG target");
      } catch (const std::exception &error) {
        throw std::invalid_argument("@Syntax error: invalid PDG[] target number '" + pdg_ +
                                    "', string to int conversion fails: " + error.what());
      }

      // Try to find the particle from PDG table, will throw exception if fails
      MParticle p = state.lts.PDG.FindByPDG(pdg);

      // Check if it has an anti-particle
      MParticle p_anti;
      bool      found_anti = false;
      try {
        p_anti     = state.lts.PDG.FindByPDG(-pdg);
        found_anti = true;
      } catch (...) {}  // do nothing,

      // Set new properties
      for (const auto &x : syntax[i].arg) {
        if (x.first == "M") {
          p.mass = aux::ParseDouble(x.second, "@PDG mass");
          if (!std::isfinite(p.mass) || p.mass < 0.0) {
            throw std::invalid_argument("@Syntax error: @PDG mass should be finite and nonnegative");
          }
          state.lts.particle_mass_overrides[std::abs(pdg)] = p.mass;
          if (found_anti) { p_anti.mass = p.mass; }
        } else if (x.first == "W") {
          p.width = aux::ParseDouble(x.second, "@PDG width");
          if (!std::isfinite(p.width) || p.width < 0.0) {
            throw std::invalid_argument("@Syntax error: @PDG width should be finite and nonnegative");
          }
          state.lts.particle_width_overrides[std::abs(pdg)] = p.width;
          p.tau = p.width > 0.0 ? PDG::hbar / p.width : 0.0;
          if (found_anti) {
            p_anti.width = p.width;
            p_anti.tau   = p.tau;
          }
        } else {
          throw std::invalid_argument("@Syntax error: unknown @PDG property '" + x.first + "', expected M or W");
        }
      }
      // Set new modified to the PDG table
      state.lts.PDG.PDG_table[pdg] = p;
      if (found_anti) { state.lts.PDG.PDG_table[-pdg] = p_anti; }

      std::cout << rang::fg::red
                << "MProcess::SetProcess: New particle properties set "
                   "with @PDG[number]{key:val} syntax:"
                << rang::fg::reset << std::endl;
      p.print();
    }
  }

  // Remove whitespace
  std::string::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
  str.erase(end_pos, str.end());

  // Parse commandline string
  std::string istate  = "";
  std::string channel = "";
  std::string mc      = "";
  ParseCMD(str, istate, channel, mc);

  // Setup subprocess
  ProcPtr.Initialize(istate, channel);
  state.lts.process.root_resonance_pdg = ProcPtr.RootResonancePDG(state.lts);

  // Check do we find the process
  if (ProcPtr.ProcessRegistry.count(str)) {
    // fine
  } else {
    throw std::invalid_argument("MProcess::SetProcess: Unknown state.process: " + str);
  }
}

// "MP[RES]<C>" is an example of valid string to be parsed
//
void MProcess::ParseCMD(const std::string &str, std::string &first, std::string &second, std::string &third) const {
  // First and Second
  bool found = false;
  ;
  for (const auto &i : indices(str)) {
    if (str[i] != '[' && !found) {
      first += str[i];
      continue;
    } else if (str[i] == '[' && !found) {
      found = true;
      continue;
    }
    if (found && str[i] != ']') { second += str[i]; }
    if (found && str[i] == ']') { break; }
  }

  // Third
  int mark1 = 0;
  int mark2 = 0;
  for (const auto &i : indices(str)) {
    if (str[i] == '<') {
      mark1 = i;
      continue;
    }
    if (str[i] == '>') {
      mark2 = i;
      break;  // First >, then break
    }
  }
  third = str.substr(mark1 + 1, mark2 - mark1 - 1);
}


namespace {

// Print one active PDG fiducial range
void PrintPDGFiducialRange(const char *name, const char *unit, const FIDCUTRANGE &range) {
  if (!range.active) { return; }
  printf("  - %-5s [min, max] = [%0.2f, %0.2f]%s \n", name, range.min, range.max, unit);
}

// Print one configured scalar fiducial range
void PrintConfiguredFiducialRange(const char *name, const char *unit, bool active, double min, double max) {
  if (!active) { return; }
  printf("- %-5s [min, max] = [%0.2f, %0.2f]%s \n", name, min, max, unit);
}

}  // namespace

// Print fiducial cuts set by user
void MProcess::PrintFiducialCuts() const {
  gra::aux::PrintBar("-");
  std::cout << std::endl;
  std::cout << rang::style::bold;
  std::cout << "Fiducial cuts:" << std::endl << std::endl;
  std::cout << rang::style::reset;

  if (state.fcuts.active == true) {
    std::cout << "Central final states" << std::endl;
    if (state.fcuts.HasCentralParticleCuts()) {
      PrintConfiguredFiducialRange("Eta", "", state.fcuts.particle_eta_active, state.fcuts.eta_min,
                                   state.fcuts.eta_max);
      PrintConfiguredFiducialRange("Pt", " GeV", state.fcuts.particle_pt_active, state.fcuts.pt_min,
                                   state.fcuts.pt_max);
      PrintConfiguredFiducialRange("Et", " GeV", state.fcuts.particle_Et_active, state.fcuts.Et_min,
                                   state.fcuts.Et_max);
      PrintConfiguredFiducialRange("Rap", "", state.fcuts.particle_rap_active, state.fcuts.rap_min,
                                   state.fcuts.rap_max);
    } else {
      std::cout << "- no cuts" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "Central system" << std::endl;
    if (state.fcuts.HasCentralSystemCuts()) {
      PrintConfiguredFiducialRange("M", " GeV", state.fcuts.system_M_active, state.fcuts.M_min, state.fcuts.M_max);
      PrintConfiguredFiducialRange("Pt", " GeV", state.fcuts.system_Pt_active, state.fcuts.Pt_min, state.fcuts.Pt_max);
      PrintConfiguredFiducialRange("Rap", "", state.fcuts.system_Rap_active, state.fcuts.Y_min, state.fcuts.Y_max);
    } else {
      std::cout << "- no cuts" << std::endl;
    }
    std::cout << std::endl;

    std::cout << "PDG based" << std::endl;
    if (state.fcuts.pdg_cuts.empty()) {
      std::cout << "- no cuts" << std::endl;
    } else {
      for (const auto &cut : state.fcuts.pdg_cuts) {
        const char *matching = cut.pdg.size() == 1 ? "all matches" : "exact multiplicity";
        std::cout << "- " << cut.SelectorString() << " (" << matching << ")" << std::endl;
        PrintPDGFiducialRange("M", " GeV", cut.M);
        PrintPDGFiducialRange("Rap", "", cut.Rap);
        PrintPDGFiducialRange("Eta", "", cut.Eta);
        PrintPDGFiducialRange("Pt", " GeV", cut.Pt);
        PrintPDGFiducialRange("Et", " GeV", cut.Et);
      }
    }
    std::cout << std::endl;

    std::cout << "Forward kinematics" << std::endl;
    if (state.fcuts.HasForwardCuts()) {
      PrintConfiguredFiducialRange("M", " GeV", state.fcuts.forward_M_active, state.fcuts.forward_M_min,
                                   state.fcuts.forward_M_max);
      PrintConfiguredFiducialRange("|t_i|", " GeV^2", state.fcuts.forward_t_active, state.fcuts.forward_t_min,
                                   state.fcuts.forward_t_max);
      PrintConfiguredFiducialRange("|t1|", " GeV^2", state.fcuts.forward_t1.active,
                                   state.fcuts.forward_t1.min, state.fcuts.forward_t1.max);
      PrintConfiguredFiducialRange("|t2|", " GeV^2", state.fcuts.forward_t2.active,
                                   state.fcuts.forward_t2.min, state.fcuts.forward_t2.max);
      PrintConfiguredFiducialRange("Xi_i", "", state.fcuts.forward_xi_active, state.fcuts.forward_xi_min,
                                   state.fcuts.forward_xi_max);
      PrintConfiguredFiducialRange("dPhi", " deg", state.fcuts.forward_dPhi_active, state.fcuts.forward_dPhi_min,
                                   state.fcuts.forward_dPhi_max);
    } else {
      std::cout << "- no cuts" << std::endl;
    }
  } else {
    std::cout << "- Not active" << std::endl;
    std::cout << std::endl;
    std::cout << "PDG based" << std::endl;
    std::cout << "- no cuts" << std::endl;
  }
  gra::aux::PrintBar("-");
  std::cout << std::endl;
  std::cout << rang::style::bold;
  std::cout << "Extra custom fiducial cuts:" << std::endl << std::endl;
  std::cout << rang::style::reset;
  std::cout << rang::fg::red;
  if (state.fcuts.active && state.usercuts != 0) {
    std::cout << "- Active ID = " << state.usercuts << " (see MUserCuts.cc)" << std::endl;
  } else {
    std::cout << "- Not active" << std::endl;
  }
  std::cout << rang::fg::reset;

  gra::aux::PrintBar("-");
  std::cout << std::endl;
  std::cout << rang::style::bold;
  std::cout << "VETO cuts:" << std::endl << std::endl;
  std::cout << rang::style::reset;

  if (state.vetocuts.active == true) {
    for (const auto &i : indices(state.vetocuts.cuts)) {
      const auto &domain = state.vetocuts.cuts[i];
      std::string sources;
      if (domain.source_forward) { sources = "forward"; }
      if (domain.source_central) {
        if (!sources.empty()) { sources += ", "; }
        sources += "central";
      }
      const char *charge = "any";
      if (domain.charge == gra::VetoCharge::Charged) {
        charge = "charged";
      } else if (domain.charge == gra::VetoCharge::Neutral) {
        charge = "neutral";
      }

      std::cout << "DOMAIN: " << i << std::endl;
      printf("- Eta   [min, max] = [%0.2f, %0.2f] \n", domain.eta_min, domain.eta_max);
      printf("- Pt    [min, max] = [%0.2f, %0.2f] GeV \n", domain.pt_min, domain.pt_max);
      std::cout << "- Sources             = [" << sources << "]" << std::endl;
      std::cout << "- Particle charge     = " << charge << std::endl;
      std::cout << std::endl;
    }
  } else {
    std::cout << "- Not active" << std::endl;
  }

  gra::aux::PrintBar("-");
  std::cout << std::endl;
  std::cout << rang::style::bold;
  std::cout << "Central system final state tree:" << std::endl << std::endl;
  std::cout << rang::style::reset;
  if (state.lts.process.ROOT_RES_ACTIVE && state.lts.process.root_decay_mode != gra::RootDecayMode::None) {
    if (state.lts.process.ROOT_RES.hel_decay.BR_set) {
      printf("Root [%s] PDG=%d, BR=%0.3E \n", state.lts.process.ROOT_RES.p.name.c_str(),
             state.lts.process.ROOT_RES.p.pdg, state.lts.process.ROOT_RES.hel_decay.BR);
    } else {
      printf("Root [%s] PDG=%d, BR=not set [&> manual] \n", state.lts.process.ROOT_RES.p.name.c_str(),
             state.lts.process.ROOT_RES.p.pdg);
    }
  }
  if (state.lts.process.root_decay_mode == gra::RootDecayMode::Isolated) {
    for (const auto &resonance : state.lts.process.RESONANCES) {
      if (resonance.second.hel_decay.BR_set) {
        printf("Central resonance [%s] PDG=%d, BR=%0.3E \n", resonance.second.p.name.c_str(), resonance.second.p.pdg,
               resonance.second.hel_decay.BR);
      } else {
        printf("Central resonance [%s] PDG=%d, BR=not set [&> manual] \n", resonance.second.p.name.c_str(),
               resonance.second.p.pdg);
      }
    }
  }
  // Print out decaytree recursively
  for (const auto &i : indices(state.lts.decaytree)) { decay::PrintTree(state.lts.decaytree[i]); }

  std::cout << std::endl;
  if (GetISOLATE()) {
    std::cout << "Isolated &> decay treatment:" << std::endl
              << "- Applied decay coupling:       g_decay = 1" << std::endl
              << "- Final-state symmetry factor:  1/S = " << 1.0 / state.symmetry_factor << std::endl
              << "- Decay phase-space volume:     not included in "
                 "cross-section weight"
              << std::endl
              << "- Daughter kinematics:          according to phase-space" << std::endl
              << "- Branching fractions:          not applied to event weights" << std::endl;
  } else {
    std::cout << "Final state symmetry factor (1/S = 1/" << state.symmetry_factor << ") applied at cross section level"
              << std::endl;
  }

  gra::aux::PrintBar("-");
  std::cout << std::endl;
}


}  // namespace gra
