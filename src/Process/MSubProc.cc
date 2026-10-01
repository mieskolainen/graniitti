// (Sub)-Processes and Amplitude containers
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <complex>
#include <random>
#include <string>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Tensor/MTensorInit.h"

// Libraries
#include "rang.hpp"

namespace gra {

// Compute beam support for the configured photoproduction continuum
bool PROC_403_TP_PHOTO::SupportsCollision(nuclear::CollisionType collision, const nlohmann::json &general, const MPDG &pdg) const {
  return Info().beams.Supports(collision) && TensorPhotoBeamError(collision, general, pdg).empty();
}

// Compute the excitation prescriptions implemented by this process
std::vector<DissociationType> MProc::DissociationModels() const {
  if (ISTATE == "yy") { return {DissociationType::Structure}; }
  if (ISTATE == "X" && (CHANNEL == "SD" || CHANNEL == "DD")) { return {DissociationType::TripleRegge}; }
  if ((ISTATE == "MP" || ISTATE == "XP" || ISTATE == "GP" || ISTATE == "TP") &&
      (CHANNEL == "RES" || CHANNEL == "RES+CON")) {
    return {DissociationType::Soft, DissociationType::Hera};
  }
  if (IncludesParameterSet(parameter_sets, MParameterSet::PhotoVM)) {
    return {DissociationType::Soft, DissociationType::Hera};
  }
  if (ISTATE == "MP" || ISTATE == "XP" || ISTATE == "GP" || ISTATE == "TP" || ISTATE == "gg" || ISTATE == "ygg") {
    return {DissociationType::Soft};
  }
  return {};
}

// Validate every configured assignment and select the active excitation prescription
void InitializeDissociation(MProcessSetup &setup, const MProc &active,
                            const std::vector<std::shared_ptr<MProc>> &processes) {
  auto &cache = RequireModelCache(setup.lts.model_cache, setup.model_tune, "Forward excitation");
  const auto param = GetNstarParam(cache);
  for (const auto &[key, rule] : param->model) {
    for (const auto type : {rule.hadron, rule.photo}) {
      if (key == "ygg" && type == DissociationType::None) { continue; }
      const bool supported = std::any_of(processes.begin(), processes.end(), [&](const auto &process) {
        if (key != process->ISTATE) { return false; }
        const auto models = process->DissociationModels();
        return std::find(models.begin(), models.end(), type) != models.end();
      });
      if (!supported) { throw std::invalid_argument("PARAM_NSTAR.MODEL: unknown class or unsupported prescription for " + key); }
    }
  }
  const auto models = active.DissociationModels();
  const auto rule = param->Model(active.ISTATE);
  const auto type = active.ISTATE == "ygg" ? rule.photo : rule.hadron;
  if (models.empty()) {
    if (setup.excitation != 0) { throw std::invalid_argument(active.Label() + " requires NSTARS=0"); }
  } else if (type == DissociationType::None ||
             (setup.excitation != 0 && std::find(models.begin(), models.end(), type) == models.end())) {
    throw std::invalid_argument("PARAM_NSTAR.MODEL: missing or unsupported prescription for " + active.Label());
  }
  setup.lts.process.DISSOCIATION = rule.hadron;
  setup.lts.process.PHOTO_DISSOCIATION = rule.photo;
}

// Add all processes
std::vector<std::shared_ptr<MProc>> MSubProc::CreateAllProcesses() const {
  std::vector<std::shared_ptr<MProc>> v;

  v.push_back(std::shared_ptr<MProc>{new PROC_001_QED_YY_HIGGS()});
  v.push_back(std::shared_ptr<MProc>{new PROC_002_QED_YY_MONOPOLIUM0()});
  v.push_back(std::shared_ptr<MProc>{new PROC_003_QED_YY_EPA()});
  v.push_back(std::shared_ptr<MProc>{new PROC_004_QED_YY_QED()});
  v.push_back(std::make_shared<PROC_005_QED_YY_YY>());
  v.push_back(std::shared_ptr<MProc>{new PROC_020_QED_YY_FLUX()});
  v.push_back(std::shared_ptr<MProc>{new PROC_030_QED_YY_DZ_EPA()});
  v.push_back(std::shared_ptr<MProc>{new PROC_031_QED_YY_DZ_FLUX()});
  v.push_back(std::shared_ptr<MProc>{new PROC_040_QED_YY_LUX_EPA()});
  v.push_back(std::shared_ptr<MProc>{new PROC_050_PHOTO_YGG_Z()});
  v.push_back(std::shared_ptr<MProc>{new PROC_051_PHOTO_YGG_VM("jpsi", 51)});
  v.push_back(std::shared_ptr<MProc>{new PROC_051_PHOTO_YGG_VM("psi(2S)", 52)});
  v.push_back(std::shared_ptr<MProc>{new PROC_051_PHOTO_YGG_VM("Upsilon(1S)", 53)});
  v.push_back(std::shared_ptr<MProc>{new PROC_051_PHOTO_YGG_VM("Upsilon(2S)", 54)});
  v.push_back(std::shared_ptr<MProc>{new PROC_051_PHOTO_YGG_VM("Upsilon(3S)", 55)});

  // Expose every generated photon family through all supported EPA fluxes
  int photon_index = 0;
  for (const auto &info : PhotonMG5ProcessInfos()) {
    v.push_back(
        std::make_shared<MGeneratedPhotonProc>("yy", info.channel,
                                               ProcessDescriptor{"Generated gamma-gamma " + info.channel, "kT-EPA",
                                                                 "yy", "MG5 amplitudes", 6 + photon_index},
                                               info.process_family));
    v.push_back(std::make_shared<MGeneratedPhotonProc>(
        "yy_DZ", info.channel,
        ProcessDescriptor{"Generated gamma-gamma " + info.channel, "Collinear Drees-Zeppenfeld EPA", "yy",
                          "MG5 amplitudes", 32 + photon_index},
        info.process_family));
    v.push_back(std::make_shared<MGeneratedPhotonProc>(
        "yy_LUX", info.channel,
        ProcessDescriptor{"Generated gamma-gamma " + info.channel, "Collinear LUX-PDF", "yy", "MG5 amplitudes",
                          41 + photon_index},
        info.process_family));
    ++photon_index;
  }

  // Add one compact MP, XP or GP process registration
  const auto add = [&v](const std::string &istate, const std::string &channel, const std::string &title,
                        const std::string &model, int id, ReggeProductionModel production_model, MReggeMode mode) {
    v.push_back(std::make_shared<MReggeCentralProc>(
        istate, channel, ProcessDescriptor{title, model, "card", "", id},
        production_model, mode));
  };
  add("MP", "RES", "Parametric resonance", "M-Pomeron", 100, ReggeProductionModel::MP, MReggeMode::Resonance);
  add("MP", "CON", "Hadron continuum 2/4/6-body", "M-Pomeron", 101, ReggeProductionModel::MP,
      MReggeMode::ContinuumTwoFourSixBody);
  add("MP", "RES+CON", "Hadron resonances + continuum 2-body", "M-Pomeron", 102, ReggeProductionModel::MP,
      MReggeMode::ResonanceContinuumTwoBody);
  add("XP", "RES", "Parametric resonance", "X-Pomeron", 200, ReggeProductionModel::XP, MReggeMode::Resonance);
  add("XP", "CON", "Hadron continuum 2/4/6-body", "X-Pomeron", 201, ReggeProductionModel::XP,
      MReggeMode::ContinuumTwoFourSixBody);
  add("XP", "RES+CON", "Hadron resonances + continuum 2-body", "X-Pomeron", 202, ReggeProductionModel::XP,
      MReggeMode::ResonanceContinuumTwoBody);
  add("GP", "RES", "Analytic Regge helicity amplitudes", "G-Pomeron", 300, ReggeProductionModel::GP, MReggeMode::Resonance);
  add("GP", "CON", "Hadron continuum 2/4/6-body", "G-Pomeron", 301, ReggeProductionModel::GP,
      MReggeMode::ContinuumTwoFourSixBody);
  add("GP", "RES+CON", "Hadron resonances + continuum 2-body", "G-Pomeron", 302, ReggeProductionModel::GP,
      MReggeMode::ResonanceContinuumTwoBody);
  add("TP", "RES", "Parametric resonance", "T-Pomeron", 400, ReggeProductionModel::TP, MReggeMode::Resonance);
  add("TP", "CON", "Hadron continuum 2-body / vector cascades", "T-Pomeron", 401, ReggeProductionModel::TP,
      MReggeMode::ContinuumTwoBody);
  add("TP", "RES+CON", "Hadron resonances + continuum 2-body", "T-Pomeron", 402, ReggeProductionModel::TP,
      MReggeMode::ResonanceContinuumTwoBody);
  v.push_back(std::make_shared<PROC_403_TP_PHOTO>());
  v.push_back(std::shared_ptr<MProc>{new PROC_600_SOFT_EL()});
  v.push_back(std::shared_ptr<MProc>{new PROC_601_SOFT_SD()});
  v.push_back(std::shared_ptr<MProc>{new PROC_602_SOFT_DD()});
  v.push_back(std::shared_ptr<MProc>{new PROC_603_SOFT_ND()});
  v.push_back(std::make_shared<MDurhamResonanceProc>("chic(0)", "gg_chic0", 10441, 700));
  v.push_back(std::make_shared<MDurhamResonanceProc>("chic(1)", "gg_chic1", 20443, 701));
  v.push_back(std::make_shared<MDurhamResonanceProc>("chic(2)", "gg_chic2", 445, 702));
  v.push_back(std::make_shared<MDurhamContinuumProc>("QCD", MDurhamMode::Parton, 703));
  v.push_back(std::make_shared<MDurhamContinuumProc>("MM", MDurhamMode::MesonPair, 704));
  v.push_back(std::make_shared<MDurhamContinuumProc>("yy", MDurhamMode::PhotonPair, 705));
  v.push_back(std::shared_ptr<MProc>{new PROC_704_DURHAM_FLUX()});
  // Expose every generated hard family through both SD and DD diffraction
  int hard_id = 798;
  for (const std::string initial : {"IPp", "IPIP"}) {
    const bool single = initial == "IPp";
    for (const auto &info : PartonMG5ProcessInfos()) {
      v.push_back(std::make_shared<MGeneratedPartonProc>(
          initial, info.channel,
          ProcessDescriptor{std::string(single ? "Single" : "Double") + " hard diffractive " + info.channel,
                            single ? "Pomeron PDF x proton PDF" : "Two Pomeron PDFs", initial,
                            "MG5 amplitudes", hard_id++},
          info.process_family));
    }
  }

  return v;
}

// Build the process registry used for the first initialization
MSubProc::MSubProc(const std::vector<std::string> &istate, const std::string &mc) {
  for (const auto &i : aux::indices(istate)) {
    ConstructDescriptions(InitialStateSelection{istate[i]}, PhaseSpaceSelection{mc});
  }
}

// Set spesific initial state and channel
void MSubProc::Initialize(const std::string &istate, const std::string &channel) {
  ISTATE           = istate;
  CHANNEL          = channel;
  process_prepared = false;
  pr.reset();
}

// Construct textual descriptions for one initial-state and phase-space
// selection
void MSubProc::ConstructDescriptions(InitialStateSelection istate, PhaseSpaceSelection mc) {
  const std::vector<std::shared_ptr<MProc>> p = CreateAllProcesses();

  for (const auto &i : aux::indices(p)) {
    if (p[i]->ISTATE == istate.value) {
      ProcessRegistry.insert(
          std::make_pair(p[i]->ISTATE + "[" + p[i]->CHANNEL + "]<" + mc.value + ">", p[i]->DESCRIPTION));
    }
  }
}

// Compute the selected physical process for initialization and queries
std::shared_ptr<MProc> MSubProc::SelectedProcess() const {
  if (pr != nullptr && pr->ISTATE == ISTATE && pr->CHANNEL == CHANNEL) { return pr; }
  for (const auto &process : CreateAllProcesses()) {
    if (process->ISTATE == ISTATE && process->CHANNEL == CHANNEL) { return process; }
  }
  throw std::invalid_argument("MSubProc::SelectedProcess: Unknown ISTATE = " + ISTATE + " or CHANNEL = " + CHANNEL);
}

// Activate the selected process
void MSubProc::ActivateProcess() {
  pr = SelectedProcess();
  if (model_tune != nullptr) { pr->BindModelTune(model_tune); }
  if (random != nullptr) { pr->BindRandom(*random); }
}

// Compute concrete native and generated final states joined to valid commands
std::vector<std::vector<std::string>> MSubProc::SupportedFinalStateRows() const {
  std::vector<std::vector<std::string>>     rows;
  const std::vector<std::shared_ptr<MProc>> processes = CreateAllProcesses();

  for (const auto &process : processes) {
    const std::vector<amplitude::Process> supported = process->Processes();
    if (supported.empty()) { continue; }

    const std::string command_prefix = process->ISTATE + "[" + process->CHANNEL + "]<";
    for (const auto &entry : ProcessRegistry) {
      if (entry.first.compare(0, command_prefix.size(), command_prefix) != 0) { continue; }
      for (const auto &final_state : supported) {
        if (final_state.topology.empty()) { continue; }
        const bool stable = std::all_of(final_state.topology.begin(), final_state.topology.end(),
                                        [](const auto &node) { return node.daughters.empty(); });
        rows.push_back({entry.first + " -> " + final_state.final_state_syntax,
                        process->Info().beams.Names(),
                        stable ? final_state.final_state_syntax : final_state.stable_final_state,
                        final_state.process_syntax});
      }
    }
  }

  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

// Check if process exists
bool MSubProc::ProcessExist(const std::string &str) const {
  if (ProcessRegistry.find(str) != ProcessRegistry.end()) { return true; }
  return false;
}

// Compute process description string
std::vector<std::string> MSubProc::GetProcessDescriptor(const std::string &str) const {
  if (!ProcessExist(str)) {
    throw std::invalid_argument("MSubProc::GetProcessDescriptor: Process by name " + str + " does not exist");
  }
  return ProcessRegistry.find(str)->second.AsVector();
}

// Bind the process-local random number generator to the active subprocess
void MSubProc::BindRandom(gra::MRandom &rng) {
  random = &rng;
  if (pr != nullptr) { pr->BindRandom(rng); }
}

// Initialize the selected physical amplitude for this process copy
void MSubProc::InitializeAmplitude(gra::MProcessSetup &setup) {
  if (setup.model_tune == nullptr || setup.soft_model == nullptr) {
    throw std::invalid_argument("MSubProc::InitializeAmplitude: null model tune");
  }
  if (model_tune != nullptr && model_tune != setup.model_tune) {
    throw std::invalid_argument("MSubProc::InitializeAmplitude: model tune identity changed");
  }
  if (soft_model != nullptr && soft_model != setup.soft_model) {
    throw std::invalid_argument("MSubProc::InitializeAmplitude: SOFT model identity changed");
  }
  model_tune = setup.model_tune;
  soft_model = setup.soft_model;
  if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) {
    pr.reset();
    ActivateProcess();
  }
  if (random != nullptr) { pr->BindRandom(*random); }
  pr->BindModelTune(model_tune);
  InitializeDissociation(setup, *pr, CreateAllProcesses());
  pr->InitializeAmplitude(setup);
  // Validate resolved decay amplitudes before sampling, after branching data exist
  if (pr->MatchProcess(setup.lts.decaytree).has_value()) {
    setup.lts.decay_structure = pr->DecayStructureFor(setup.lts);
  }
  process_prepared = true;
}

// Sample a process-specific event color flow
bool MSubProc::SampleColorFlow(gra::LORENTZSCALAR &lts) { return pr == nullptr || pr->SampleColorFlow(lts); }

// Compute the root resonance PDG code advertised by the selected process
int MSubProc::RootResonancePDG(const gra::LORENTZSCALAR &lts) const {
  return SelectedProcess()->RootResonancePDG(lts);
}

// Prepare reusable amplitude state for one phase-space point
void MSubProc::PrepareBareAmplitude(gra::LORENTZSCALAR &lts) {
  RequirePreparedProcess();
  DecayStructureFor(lts);
  pr->InitializeWorkerAmplitude(lts);
  pr->PrepareAmp2(lts);
}

// Resolve and store the active amplitude decay structure
DecayStructure MSubProc::DecayStructureFor(gra::LORENTZSCALAR &lts) {
  if (!process_prepared) { lts.decay_structure = {}; }
  if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) {
    pr.reset();
    ActivateProcess();
  }
  if (!process_prepared) { lts.decay_structure = pr->DecayStructureFor(lts); }
  if (lts.decay_structure.type == DecayType::JacobWickCoherent ||
      lts.decay_structure.type == DecayType::JacobWickIncoherent) {
    lts.decay_structure.type = decay::JacobWickStructure(lts).type;
  }
  return lts.decay_structure;
}

// Match the active amplitude process against one decay tree
std::optional<amplitude::Process> MSubProc::MatchProcess(const std::vector<MDecayBranch> &decaytree) {
  if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) {
    pr.reset();
    ActivateProcess();
  }
  return pr->MatchProcess(decaytree);
}

// Compute an optional amplitude constraint for one generic-parton mode
std::optional<bool> MSubProc::AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) {
  if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) {
    pr.reset();
    ActivateProcess();
  }
  return pr->AcceptsPartonMode(decaytree);
}

// Compute processes for the selected process and activate it when needed
std::vector<amplitude::Process> MSubProc::Processes() {
  if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) {
    pr.reset();
    ActivateProcess();
  }
  return pr->Processes();
}

// Amplitude function
double MSubProc::GetBareAmplitude2(gra::LORENTZSCALAR &lts) {
  RequirePreparedProcess();
  DecayStructureFor(lts);
  pr->InitializeWorkerAmplitude(lts);
  pr->ResetEvaluationStatus();
  lts.muF    = 0.0;
  lts.muR    = 0.0;
  lts.scalup = 0.0;
  return pr->Amp2(lts);
}

// Amplitude function using previously prepared amplitude state
double MSubProc::GetPreparedBareAmplitude2(gra::LORENTZSCALAR &lts) {
  RequirePreparedProcess();
  DecayStructureFor(lts);
  pr->InitializeWorkerAmplitude(lts);
  pr->ResetEvaluationStatus();
  return pr->PreparedAmp2(lts);
}

// Compute the status of the most recent event-local matrix-element evaluation
mg5helas::EvaluationStatus MSubProc::EvaluationStatus() const {
  return pr == nullptr ? mg5helas::EvaluationStatus::AmplitudeFailure : pr->EvaluationStatus();
}

}  // namespace gra
