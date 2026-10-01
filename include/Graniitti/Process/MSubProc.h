// (Sub)-Processes and Amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSUBPROC_H
#define MSUBPROC_H

// C++
#include <algorithm>
#include <complex>
#include <cstdint>
#include <functional>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Amplitude/Photon/AMP_yy_yy.h"
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MGamma.h"
#include "Graniitti/Photon/MPhotoVM.h"
#include "Graniitti/Photon/MPhotoZ.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/QCD/MDurham.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggePhoto.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tensor/MTensorPomeron.h"

using gra::aux::indices;

namespace gra {

using math::abs2;
using math::msqrt;

// Collect stable final-state PDGs in decay-tree momentum order
inline std::vector<int> PhotonStablePDGs(const std::vector<MDecayBranch> &tree) {
  const auto       leaves = mg5::StableDecayLeaves(tree);
  std::vector<int> pdgs;
  pdgs.reserve(leaves.size());
  for (const auto *leaf : leaves) { pdgs.push_back(leaf->p.pdg); }
  return pdgs;
}

// Structured process descriptor used by help tables and setup summaries
struct ProcessDescriptor {
  std::string process;
  std::string model;
  std::string channels;
  std::string comments;
  int         display_order = 9999;

  // Compute the vector form used in setup summaries
  std::vector<std::string> AsVector() const {
    std::vector<std::string> out = {process, model, channels};
    if (!comments.empty()) { out.push_back(comments); }
    return out;
  }
};

// Immutable parameter families prepared by a process before worker copies
enum class MParameterSet : std::uint32_t {
  None         = 0,
  Regge        = 1U << 0,
  Tensor       = 1U << 1,
  MonopolePair = 1U << 2,
  Monopolium   = 1U << 3,
  Durham       = 1U << 4,
  PhotoZ       = 1U << 5,
  PhotoVM      = 1U << 6,
  PhotonPDF    = 1U << 7,
  ProtonPDF    = 1U << 8
};

// Combine immutable parameter families into one declarative process mask
constexpr MParameterSet operator|(const MParameterSet lhs, const MParameterSet rhs) noexcept {
  return static_cast<MParameterSet>(static_cast<std::uint32_t>(lhs) | static_cast<std::uint32_t>(rhs));
}

// Compute whether one immutable parameter family is present in a process mask
constexpr bool IncludesParameterSet(const MParameterSet mask, const MParameterSet value) noexcept {
  return (static_cast<std::uint32_t>(mask) & static_cast<std::uint32_t>(value)) != 0U;
}

// Helper for copy constructors of unique_ptr
template <class T>
std::unique_ptr<T> copy_unique(const std::unique_ptr<T> &source) {
  return source ? std::make_unique<T>(*source) : nullptr;
}

// Build generated photon processes and the generic physical pair selector
inline std::shared_ptr<const amplitude::ProcessDefinition> PhotonContinuumProcesses() {
  auto                                                             generated = amplitude::Processes("PHOTON");
  std::vector<std::shared_ptr<const amplitude::ProcessDefinition>> definitions;
  definitions.push_back(std::make_shared<amplitude::ProcessRegistry>(std::move(generated)));
  definitions.push_back(MGamma::ProcessDefinitionFor(MGammaMode::FermionPair));
  return std::make_shared<amplitude::ProcessAlternatives>(std::move(definitions));
}

// Build an explicit process for an inline unit amplitude
inline std::shared_ptr<const amplitude::ProcessDefinition> UnitAmplitudeProcessDefinitionFor(
    const std::string &process_name) {
  return std::make_shared<amplitude::AnalyticProcess>(
      "UNIT", process_name, "unit-amplitude final state", DecayStructure{},
      [](const std::vector<MDecayBranch> &) { return true; },
      [](const LORENTZSCALAR &) {
        return DecayStructure{};
      });
}

// Abstract process base class
class MProc : public amplitude::ProcessFamily {
 public:
  // Build process commands used by event channel selection and help tables
  MProc(const std::string &i, const std::string &j, const ProcessDescriptor &k,
        std::shared_ptr<const amplitude::ProcessDefinition> definition, ScreeningMetadata screening,
        MParameterSet parameter_sets_in = MParameterSet::None, ProcessInfo info_in = {})
      : amplitude::ProcessFamily(std::move(definition)),
        ISTATE(i),
        CHANNEL(j),
        DESCRIPTION(k),
        screening_definition(screening),
        parameter_sets(parameter_sets_in), info(std::move(info_in)) {}
  std::string       ISTATE;
  std::string       CHANNEL;
  ProcessDescriptor DESCRIPTION;

  // Process class containers
  std::unique_ptr<MGamma>           Gamma               = nullptr;
  std::unique_ptr<MDurham>          Durham              = nullptr;
  std::unique_ptr<MRegge>           Regge               = nullptr;
  std::unique_ptr<MPhotoZ>          PhotoZ              = nullptr;
  std::unique_ptr<MPhotoVM>         PhotoVM             = nullptr;
  std::unique_ptr<MTensorPomeron>   Tensor              = nullptr;
  std::unique_ptr<PhotonMG5Process> PhotonMatrixElement = nullptr;
  std::optional<amplitude::Process> PhotonProcess;

  // Bind the immutable SOFT model shared by every process worker
  void BindSoftModel(const SoftModelPtr &model) {
    if (model == nullptr) { throw std::invalid_argument("MProc::BindSoftModel: null SOFT model"); }
    if (soft_model != nullptr && soft_model != model) {
      throw std::invalid_argument("MProc::BindSoftModel: SOFT model identity changed");
    }
    soft_model = model;
  }

  // Bind the immutable complete model tune
  void BindModelTune(const MModelTunePtr &tune) {
    if (tune == nullptr) { throw std::invalid_argument("MProc::BindModelTune: null model tune"); }
    if (model_tune != nullptr && model_tune != tune) {
      throw std::invalid_argument("MProc::BindModelTune: model tune identity changed");
    }
    model_tune = tune;
    BindSoftModel(model_tune->Soft());
  }

  // Select and initialize the gamma-gamma amplitude before sampling
  void InitGamma(gra::LORENTZSCALAR &lts) {
    if (CHANNEL != "EPA" ||
        (AssertN(2, lts.decaytree.size()) && AssertLeptonQuarkMonopolePair(PDGlist(lts)))) {
      if (Gamma == nullptr) { Gamma = make_unique<MGamma>(lts, GetModelTune(), Family()); }
      return;
    }
    auto &cache = RequireModelCache(lts.model_cache, GetModelTune(), Label());
    PhotonProcess = GeneratedPhotonProcess(lts);
    if (!PhotonProcess.has_value()) { throw std::invalid_argument(Label() + ": no generated photon process for the decay tree"); }
    PhotonMatrixElement = CreatePhotonMG5Process(*PhotonProcess);
    if (PhotonMatrixElement == nullptr) { throw std::invalid_argument(Label() + ": generated photon matrix element is unavailable"); }
    if (PhotonMatrixElement->AlphaSPower() > 0) {
      if (lts.GlobalPdfPtr == nullptr) {
        lts.GlobalPdfPtr = cache.pdf.GetPDF(lts.LHAPDFSET, 0);
      }
      if (!lts.GlobalPdfPtr->hasAlphaS()) { throw std::invalid_argument(Label() + ": PDF has no alpha_s evolution"); }
    }
  }

  // Initialize the Durham amplitude with the bound random generator
  void InitDurham(gra::LORENTZSCALAR &lts) {
    if (random == nullptr) { throw std::invalid_argument("MProc::InitDurham: Random number generator is not bound"); }
    if (Durham == nullptr) { Durham = make_unique<MDurham>(lts, GetModelTune(), *random, Family()); }
  }

  // Initialize one Tensor runtime
  void InitTensor(gra::LORENTZSCALAR &lts, MTensorPomeronMode mode) {
    if (Tensor == nullptr) { Tensor = make_unique<MTensorPomeron>(lts, GetModelTune(), Family(), mode); }
  }
  void InitRegge(gra::LORENTZSCALAR &lts) {
    if (Regge == nullptr) { Regge = make_unique<MRegge>(lts, GetModelTune(), Family()); }
  }

  // Lazily initialize the exclusive Z photoproduction matrix_element
  void InitPhotoZ(gra::LORENTZSCALAR &lts) {
    if (PhotoZ == nullptr) { PhotoZ = make_unique<MPhotoZ>(lts, GetModelTune(), Family()); }
  }

  // Lazily initialize the exclusive vector-meson photoproduction matrix_element
  void InitPhotoVM(gra::LORENTZSCALAR &lts) {
    if (PhotoVM == nullptr) { PhotoVM = make_unique<MPhotoVM>(lts, GetModelTune(), CHANNEL, Family()); }
  }

  // ---------------------------------------------------------------------

  virtual ~MProc() {}  // Needs to be virtual

  // Initialize the selected physical amplitude before worker copies are made
  void InitializeAmplitude(gra::MProcessSetup &setup) {
    if (process_initialized) { return; }
    BindModelTune(setup.model_tune);
    InitializeParameterSets(setup);
    ValidateRuntime();
    InitializeBranching(setup);
    setup.lts.hamp.Configure(Screening(setup));
    process_initialized = true;
  }

  // Construct only the mutable amplitude runtime for one worker copy
  void InitializeWorkerAmplitude(gra::LORENTZSCALAR &lts) {
    if (runtime_initialized) { return; }
    InitializeProcess(lts);
    runtime_initialized = true;
  }

  // Check whether fixed physical process state has been initialized
  bool IsProcessInitialized() const { return process_initialized; }

  // Check whether the mutable amplitude runtime has been initialized
  bool IsRuntimeInitialized() const { return runtime_initialized; }

  virtual double Amp2(gra::LORENTZSCALAR &lts) = 0;

  // Convert external-helicity amplitudes to process-normalized physical weights
  virtual AmplitudeWeights NormalizeAmplitude(const MHelicityAmplitudes &amplitude) const {
    std::vector<double> helicity_norm;
    helicity_norm.reserve(amplitude.size());
    for (const auto &value : amplitude) { helicity_norm.push_back(std::norm(value)); }
    return NormalizeProjectedAmplitude(helicity_norm, amplitude.metadata);
  }

  // Normalize external-helicity norms supplied by a resolved screening projection
  virtual AmplitudeWeights NormalizeProjectedAmplitude(const std::vector<double> &helicity_norm,
                                                       const ScreeningMetadata   &metadata) const {
    const double normalization = metadata.amplitude_normalization;
    AmplitudeWeights result;
    result.helicity_weight.reserve(helicity_norm.size());
    for (const double norm : helicity_norm) {
      if (norm < 0.0) {
        throw AmplitudeFailure("MProc::NormalizeAmplitude: invalid helicity norm");
      }
      const double weight = normalization * norm;
      result.helicity_weight.push_back(weight);
      result.total_weight += weight;
    }
    if (!std::isfinite(result.total_weight)) {
      throw AmplitudeFailure("MProc::NormalizeAmplitude: non-finite physical weight");
    }
    return result;
  }

  // Normalize a proton-screened amplitude in its resulting spin representation
  virtual AmplitudeWeights NormalizeScreenedAmplitude(const MHelicityAmplitudes &amplitude,
                                                      const bool                 dense_proton_spin) const {
    (void)dense_proton_spin;
    return NormalizeAmplitude(amplitude);
  }

  // Reset event-local matrix-element evaluation status
  void ResetEvaluationStatus() { evaluation_status = mg5helas::EvaluationStatus::Success; }

  // Compute event-local matrix-element evaluation status
  mg5helas::EvaluationStatus EvaluationStatus() const { return evaluation_status; }

  // Prepare reusable amplitude state for one phase-space point
  virtual void PrepareAmp2(gra::LORENTZSCALAR &lts) { (void)lts; }
  // Evaluate with previously prepared reusable amplitude state
  virtual double PreparedAmp2(gra::LORENTZSCALAR &lts) { return Amp2(lts); }
  // Compute the capabilities declared by this process implementation
  const ProcessInfo &Info() const { return info; }

  // Compute the excitation prescriptions implemented by this process
  std::vector<DissociationType> DissociationModels() const;

  // Compute collision support with the model cards used for process help
  virtual bool SupportsCollision(nuclear::CollisionType collision, const nlohmann::json &, const MPDG &) const {
    return info.beams.Supports(collision);
  }

  // Compute the root resonance PDG code for analytic resonance processes
  virtual int RootResonancePDG(const gra::LORENTZSCALAR &lts) const {
    (void)lts;
    return 0;
  }

  // Bind the process-local random number generator before lazy amplitude
  // construction
  void BindRandom(gra::MRandom &rng) {
    if (random != &rng) {
      random = &rng;
      Durham.reset();
      runtime_initialized = false;
    }
  }

  // Sample a process-specific event color flow after the final amplitude is
  // known
  bool SampleColorFlow(gra::LORENTZSCALAR &lts) {
    if (Durham != nullptr) {
      if (!Durham->SampleColorFlow(lts)) {
        SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
        return false;
      }
      return true;
    }
    if (PhotoZ != nullptr) {
      PhotoZ->SampleColorFlow(lts);
      return true;
    }
    if (PhotoVM != nullptr) {
      PhotoVM->SampleColorFlow(lts);
      return true;
    }
    if (Gamma != nullptr) {
      if (!Gamma->SampleColorFlow(lts)) {
        SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
        return false;
      }
      return true;
    }
    const bool generated = PhotonMatrixElement != nullptr && PhotonProcess.has_value();
    const bool color_ready = !generated || !PhotonMatrixElement->HasFinalStateColor()
                                 ? ClearHardColorFlow(lts)
                                 : random != nullptr && PhotonMatrixElement->SampleColorFlow(lts, *random);
    if (!color_ready) {
      SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
      return false;
    }
    return true;
  }

  // Evaluate the initialized generated photon matrix element with exact finite-Nc color
  double GeneratedGammaGammaCON(gra::LORENTZSCALAR &lts, bool coherent_epa) {
    PhotonMG5Process &matrix_element = *PhotonMatrixElement;
    if (matrix_element.AlphaSPower() > 0) {
      if (!std::isfinite(lts.alphaQCD) || lts.alphaQCD < 0.0 ||
          (!(lts.alphaQCD > 0.0) && !flux::SetPhotonAlphaS(lts))) {
        SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
        lts.hamp.clear();
        return 0.0;
      }
    }
    const auto result = matrix_element.Evaluate(lts, lts.alphaQCD, coherent_epa);
    SetEvaluationStatus(result.status);
    if (!result.Valid()) { return 0.0; }
    return result.amp2;
  }

  // Evaluate one Durham hard kernel and propagate its event-local status
  double EvaluateDurhamQCD(gra::LORENTZSCALAR &lts, const std::string &process) {
    const double amp2 = Durham->DurhamQCD(lts, process);
    SetEvaluationStatus(Durham->EvaluationStatus());
    return amp2;
  }

  // Compute the exact generated photon process for one decay tree
  std::optional<amplitude::Process> GeneratedPhotonProcess(const gra::LORENTZSCALAR &lts) const {
    const auto process = MatchProcess(lts.decaytree);
    if (!process.has_value() || !HasPhotonMG5Process(*process)) { return std::nullopt; }
    return process;
  }

  // Apply EPA sector weights to every event-local generated color row
  bool WeightHardColorFlows(LORENTZSCALAR &lts, double common, const flux::EPASectorWeights &weights) {
    auto flows = lts.hard_color_flows;
    for (auto &flow : flows) {
      flow.screened_weight.reset();
      if (!flux::ApplyEPAAmplitudeWeights(flow.amplitudes, common, weights)) { return false; }
    }
    lts.hard_color_flows = std::move(flows);
    return true;
  }

  // Evaluate shifted photon sources with the node-local EPA phase space
  double EPAHardAmp2(gra::LORENTZSCALAR &lts) {
    if (!lts.epa_hard.Ready()) {
      SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
      return 0.0;
    }
    const auto result = mg5helas::ContractEPAHard(lts, lts.epa_hard);
    SetEvaluationStatus(result.status);
    if (!result.Valid()) { return 0.0; }
    const double phase_space = flux::ExactktEPAPhaseSpaceFactor(lts);
    if (!(phase_space > 0.0)) {
      SetEvaluationStatus(mg5helas::EvaluationStatus::KinematicsFailure);
      return 0.0;
    }
    const flux::EPAWeight weight = flux::ApplyktEPAcurrents(result.amp2, lts, phase_space);
    if (!lts.hard_color_flows.empty()) {
      if (!WeightHardColorFlows(lts, weight.amplitude_scale, weight.sector)) {
        SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
        return 0.0;
      }
    }
    return weight.amp2;
  }

  // Apply kT-EPA normalization to the hard and shower-flow amplitudes
  double ApplyktEPA(const double amp2, gra::LORENTZSCALAR &lts) {
    const flux::EPAWeight weight = flux::ApplyktEPAfluxes(amp2, lts);
    if (!lts.hard_color_flows.empty()) {
      if (!WeightHardColorFlows(lts, weight.amplitude_scale, weight.sector)) {
        SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
        return 0.0;
      }
    }
    return weight.amp2;
  }

  // Evaluate the initialized continuum amplitude with the selected photon flux
  double GammaGammaCON(gra::LORENTZSCALAR &lts, bool coherent_epa) {
    if (PhotonMatrixElement != nullptr) { return GeneratedGammaGammaCON(lts, coherent_epa); }
    const auto result = Gamma->yyffbar(lts, coherent_epa);
    SetEvaluationStatus(result.status);
    return result.Valid() ? result.amp2 : 0.0;
  }

  // Access the immutable SOFT model bound to this process worker
  const SoftModelPtr &GetSoftModel() const noexcept { return soft_model; }

  // Access the immutable complete model tune bound to this process worker
  const MModelTunePtr &GetModelTune() const noexcept { return model_tune; }

  // Compute subprocess label used in diagnostics
  std::string Label() const { return ISTATE + "[" + CHANNEL + "]"; }

  // Throw for process entries whose external matrix-element matrix_element is
  // absent
  double ThrowMissingMatrixElement(const std::string &matrix_element) const {
    throw std::invalid_argument(Label() + " requires " + matrix_element +
                                " matrix elements, but the matrix_element is "
                                "not linked in this checkout");
  }

  // Require that a resonance call produced a non-empty helicity vector
  void RequireHelicityAmplitudes(const std::vector<std::complex<double>> &hamp, const std::string &context) const {
    if (hamp.empty()) { throw AmplitudeFailure(context + " produced no resonance helicity amplitudes"); }
  }

  // Assert the final state list
  bool AssertN(const std::vector<int> &reference, const std::vector<int> &input) const {
    if (reference.size() != input.size()) { return false; }
    return std::is_permutation(reference.begin(), reference.end(), input.begin());
  }

  // Assert the number of final states
  bool AssertN(int reference, int input) const { return (reference == input) ? true : false; }

  // Assert a five-flavour Durham q q~ state with an optional final gluon
  bool AssertDurhamQuarkState(const std::vector<int> &input, bool with_gluon) const {
    for (const int flavour : {1, 2, 3, 4, 5}) {
      const std::vector<int> reference =
          with_gluon ? std::vector<int>{flavour, -flavour, PDG::PDG_gluon} : std::vector<int>{flavour, -flavour};
      if (AssertN(reference, input)) { return true; }
    }
    return false;
  }

  // Assert lepton or quark pair
  bool AssertLeptonQuarkMonopolePair(const std::vector<int> &input) const {
    // PDG identifiers
    static const MMatrix<int> X = {{11, -11}, {13, -13}, {15, -15}, {1, -1}, {2, -2},
                                   {3, -3},   {4, -4},   {5, -5},   {6, -6}, {PDG::PDG_monopole, -PDG::PDG_monopole}};

    for (std::size_t i = 0; i < X.size_row(); ++i) {
      if (AssertN({X[i][0], X[i][1]}, input)) { return true; }
    }
    return false;
  }

  // Get PDG list
  std::vector<int> PDGlist(const LORENTZSCALAR &lts) const {
    std::vector<int> L(lts.decaytree.size(), 0);
    for (const auto &i : aux::indices(lts.decaytree)) { L[i] = lts.decaytree[i].p.pdg; }
    return L;
  }

  void ThrowUnknownFinalState() const {
    throw std::invalid_argument(ISTATE + "[" + CHANNEL + "]" + " with unsupported final state");
  }

 protected:
  // Prepare every immutable parameter family declared by this process route
  void InitializeParameterSets(gra::MProcessSetup &setup) const {
    if (IncludesParameterSet(parameter_sets, MParameterSet::Regge)) { MRegge::InitializeParameters(setup); }
    if (IncludesParameterSet(parameter_sets, MParameterSet::Tensor)) { MTensorPomeron::InitializeParameters(setup); }
    if (IncludesParameterSet(parameter_sets, MParameterSet::MonopolePair)) {
      MGamma::InitializeParameters(setup, MGammaMode::FermionPair);
    }
    if (IncludesParameterSet(parameter_sets, MParameterSet::Monopolium)) {
      MGamma::InitializeParameters(setup, MGammaMode::Monopolium);
    }
    if (IncludesParameterSet(parameter_sets, MParameterSet::Durham)) { MDurham::InitializeParameters(setup); }
    if (IncludesParameterSet(parameter_sets, MParameterSet::PhotoZ)) { MPhotoZ::InitializeParameters(setup); }
    if (IncludesParameterSet(parameter_sets, MParameterSet::PhotoVM)) { MPhotoVM::InitializeParameters(setup); }

    const bool needs_photon_pdf = IncludesParameterSet(parameter_sets, MParameterSet::PhotonPDF);
    const bool needs_proton_pdf = IncludesParameterSet(parameter_sets, MParameterSet::ProtonPDF);
    if (needs_photon_pdf || needs_proton_pdf) {
      if (setup.lts.LHAPDFSET.empty() || setup.lts.LHAPDFSET == "null") {
        throw std::invalid_argument(Label() + " requires a configured LHAPDF set");
      }
      if (setup.lts.model_cache == nullptr) {
        throw std::invalid_argument("Amplitude initialization requires run owned model caches");
      }
      setup.lts.GlobalPdfPtr = setup.lts.model_cache->pdf.GetPDF(setup.lts.LHAPDFSET, 0);
    }
    if (needs_photon_pdf && !setup.lts.GlobalPdfPtr->hasFlavor(PDG::PDG_gamma)) {
      throw std::invalid_argument(Label() + ": LHAPDF set '" + setup.lts.LHAPDFSET + "' has no photon distribution");
    }
  }

  // Initialize process-specific parameters and fixed amplitude structures
  virtual void InitializeProcess(gra::LORENTZSCALAR &lts) = 0;

  // Validate one optional external amplitude runtime before worker copies
  virtual void ValidateRuntime() const {}

  // Initialize process-independent root decay and phase-space state
  virtual void InitializeBranching(gra::MProcessSetup &setup) {
    InitializeGenericProcessState(setup);
  }

  // Resolve the complete immutable screening layout once after branching setup
  virtual ScreeningMetadata Screening(const gra::MProcessSetup &setup) const {
    ScreeningMetadata layout = screening_definition;
    if (ISTATE == "ygg" && setup.lts.upc_model != nullptr) {
      layout                = {ScreeningSpinBasis::Scalar, ProtonScreeningMode::None, 1.0, layout.amplitude_type};
      layout.spin_rows      = 1;
      layout.forward_noflip = true;
      return layout;
    }
    if (layout.spin_basis == ScreeningSpinBasis::ProtonHelicity) {
      if (layout.amplitude_type == ScreeningAmplitudeType::ElasticCNI) {
        layout.forward_noflip = false;
      } else if (IncludesParameterSet(parameter_sets, MParameterSet::Tensor) && ISTATE == "TP") {
        layout.forward_noflip = GetTensorParam(*setup.lts.model_cache, setup.lts.PDG, setup.lts.process.RESONANCES)->FORWARD_NOFLIP;
      } else {
        layout.forward_noflip = setup.lts.process.FORWARD_NOFLIP;
      }
      layout.spin_rows = layout.forward_noflip ? 4 : 16;
    }
    return layout;
  }

  // Store one event-local matrix-element failure without throwing
  void SetEvaluationStatus(mg5helas::EvaluationStatus status) { evaluation_status = status; }

  gra::MRandom              *random              = nullptr;  // Non-owning process-local random number generator
  mg5helas::EvaluationStatus evaluation_status   = mg5helas::EvaluationStatus::Success;
  bool                       process_initialized = false;
  bool                       runtime_initialized = false;
  const ScreeningMetadata    screening_definition;
  MModelTunePtr              model_tune;
  SoftModelPtr               soft_model;
  MParameterSet              parameter_sets = MParameterSet::None;
  const ProcessInfo          info;
};

// Route one generated parton amplitude family through its stable matrix element
// ID
class MGeneratedPartonProc : public MProc {
 public:
  // Construct one generated parton process descriptor and process definition
  MGeneratedPartonProc(const std::string &istate, const std::string &channel, const ProcessDescriptor &description,
                       const std::string &process_family)
      : MProc(istate, channel, description,
              std::make_shared<amplitude::ProcessRegistry>(amplitude::Processes(process_family)),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::TripleRegge, 1.0},
              istate == "IPp" ? MParameterSet::ProtonPDF : MParameterSet::None, {.jw_helicity_algebra = false}),
        parton_process_family(process_family) {}

  // Destroy the event-local generated parton matrix element
  ~MGeneratedPartonProc() override = default;

  // Construct the selected generated parton matrix element
  void InitializeProcess(gra::LORENTZSCALAR &lts) override {
    (void)lts;
    parton_process = CreatePartonMG5Process(parton_process_family);
    if (parton_process == nullptr) {
      throw std::invalid_argument(Label() + " has no generated parton matrix element for " + parton_process_family);
    }
  }

  // Prepare one generated parton matrix element for the current event
  void PrepareAmp2(gra::LORENTZSCALAR &lts) override {
    parton_evaluation = {};
    lts.hard_color_flows.clear();
    parton_evaluation.status = parton_process->Prepare(lts, lts.alphaQCD);
    SetEvaluationStatus(parton_evaluation.status);
  }

  // Evaluate one generated parton matrix element for the current event
  double Amp2(gra::LORENTZSCALAR &lts) override {
    PrepareAmp2(lts);
    if (!parton_evaluation.Valid()) { return 0.0; }
    return PreparedAmp2(lts);
  }

  // Evaluate one prepared generated parton matrix element for the current event
  double PreparedAmp2(gra::LORENTZSCALAR &lts) override {
    parton_evaluation = parton_process->EvaluatePrepared(lts, lts.alphaQCD);
    SetEvaluationStatus(parton_evaluation.status);
    lts.hard_color_flows = parton_evaluation.color_flows;
    return parton_evaluation.Valid() ? parton_evaluation.amp2 : 0.0;
  }

 private:
  const std::string                 parton_process_family;
  std::unique_ptr<PartonMG5Process> parton_process;
  PartonMG5Evaluation               parton_evaluation;
};

// Route one generated gamma-gamma family through a selected photon flux
class MGeneratedPhotonProc : public MProc {
 public:
  // Construct one generated photon process descriptor and process definition
  MGeneratedPhotonProc(const std::string &istate, const std::string &channel, const ProcessDescriptor &description,
                       const std::string &process_family)
      : MProc(istate, channel, description,
              std::make_shared<amplitude::ProcessRegistry>(amplitude::Processes(process_family)),
              {ScreeningSpinBasis::ProtonIdentity,
               istate == "yy" ? ProtonScreeningMode::ForwardExcitation : ProtonScreeningMode::None, 0.25},
              istate == "yy_LUX" ? MParameterSet::PhotonPDF : MParameterSet::None,
              {.jw_helicity_algebra = false, .beams = istate == "yy" ? nuclear::BeamSupport::EPA() : nuclear::BeamSupport{}}),
        photon_process_family(process_family) {}

  // Destroy the event-local generated photon matrix element
  ~MGeneratedPhotonProc() override = default;

  // Construct the selected generated photon matrix element
  void InitializeProcess(LORENTZSCALAR &lts) override {
    PhotonMatrixElement = CreatePhotonMG5Process(photon_process_family);
    if (PhotonMatrixElement == nullptr) {
      throw std::invalid_argument(Label() + " has no generated photon matrix_element for " + photon_process_family);
    }
    PhotonProcess = PhotonMatrixElement->MatchProcess(lts.decaytree);
    if (!PhotonProcess.has_value()) { throw std::invalid_argument(Label() + ": incompatible generated photon decay tree"); }
  }

  // Evaluate one generated photon family with its selected flux
  double Amp2(LORENTZSCALAR &lts) override {
    const bool coherent_epa = ISTATE == "yy";
    const auto result       = PhotonMatrixElement->Evaluate(lts, lts.alphaQCD, coherent_epa);
    SetEvaluationStatus(result.status);
    if (!result.Valid()) { return 0.0; }
    if (ISTATE == "yy") { return ApplyktEPA(result.amp2, lts); }
    if (ISTATE == "yy_DZ") { return flux::ApplyDZfluxes(result.amp2, lts); }
    if (ISTATE == "yy_LUX") { return flux::ApplyLUXfluxes(result.amp2, lts); }
    SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
    return 0.0;
  }

 private:
  const std::string photon_process_family;
};

// -----------------------------------------------------------------------
// Gamma-Gamma processes
// -----------------------------------------------------------------------

class PROC_001_QED_YY_HIGGS : public MProc {
 public:
  PROC_001_QED_YY_HIGGS()
      : MProc("yy", "Higgs", {"SM Higgs", "kT-EPA", "yy", "", 1}, MGamma::ProcessDefinitionFor(MGammaMode::Higgs),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25}, MParameterSet::None,
              {.beams = nuclear::BeamSupport::EPA(),
               .example_decay = RootDecayMode::Isolated}) {}
  // Compute the SM Higgs root resonance PDG code
  virtual int RootResonancePDG(const gra::LORENTZSCALAR &lts) const {
    (void)lts;
    return 25;
  }
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    const double amp2 = Gamma->yyHiggs(lts);
    if (!std::isfinite(amp2) || amp2 < 0.0 || lts.hamp.empty()) {
      SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
      return 0.0;
    }
    return ApplyktEPA(amp2, lts);
  }
};
class PROC_002_QED_YY_MONOPOLIUM0 : public MProc {
 public:
  PROC_002_QED_YY_MONOPOLIUM0()
      : MProc("yy", "monopolium(0)", {"Monopolium (J=0)", "kT-EPA", "yy", "", 2},
              MGamma::ProcessDefinitionFor(MGammaMode::Monopolium),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25},
              MParameterSet::Monopolium,
              {.beams = nuclear::BeamSupport::EPA(),
               .example_decay = RootDecayMode::Isolated}) {}
  // Compute the default monopolium root resonance PDG code
  virtual int RootResonancePDG(const gra::LORENTZSCALAR &lts) const {
    (void)lts;
    return PDG::PDG_monopolium;
  }
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    const double amp2 = Gamma->yyMP(lts);
    if (!std::isfinite(amp2) || amp2 < 0.0 || lts.hamp.empty()) {
      SetEvaluationStatus(mg5helas::EvaluationStatus::AmplitudeFailure);
      return 0.0;
    }
    return ApplyktEPA(amp2, lts);
  }
};
class PROC_003_QED_YY_EPA : public MProc {
 public:
  PROC_003_QED_YY_EPA()
      : MProc("yy", "EPA", {"Continuum l+l-, qqbar, monopolepair", "kT-EPA", "yy", "MG5 amplitudes", 3},
              PhotonContinuumProcesses(), {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::None, 0.25},
              MParameterSet::MonopolePair | MParameterSet::ProtonPDF,
              {.jw_helicity_algebra = false, .beams = nuclear::BeamSupport::EPA()}) {}
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    double amp2 = GammaGammaCON(lts, true);
    return ApplyktEPA(amp2, lts);
  }
};
class PROC_004_QED_YY_QED : public MProc {
 public:
  // Construct the proton QED process from its fermion pair amplitude definition
  PROC_004_QED_YY_QED()
      : MProc("yy", "QED", {"Continuum l+l-, qqbar", "Full QED", "yy", "", 4},
              MTensorPomeron::ProcessDefinitionFor(MTensorPomeronMode::QED),
              {ScreeningSpinBasis::ProtonHelicity, ProtonScreeningMode::ForwardExcitation, 0.25},
              MParameterSet::Tensor, {.beams = {{nuclear::CollisionType::PP}, false}}) {}

  // Initialize the proton QED amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitTensor(lts, MTensorPomeronMode::QED); }

  // Evaluate the full proton QED fermion pair amplitude
  double Amp2(gra::LORENTZSCALAR &lts) override { return Tensor->ME4(lts, TensorContinuumMode::QED); }
};

// Route the analytic charged Standard Model light-by-light amplitude through
// kT-EPA
class PROC_005_QED_YY_YY : public MProc {
 public:
  // Construct the exact analytic gamma gamma to gamma gamma process
  PROC_005_QED_YY_YY()
      : MProc("yy", "yy", {"SM light-by-light scattering", "kT-EPA", "yy", "Analytic charged SM loops", 5},
              std::make_shared<amplitude::ProcessRegistry>(LightByLightProcesses()),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25}, MParameterSet::None,
              {.jw_helicity_algebra = false, .beams = nuclear::BeamSupport::EPA()}) {}

  // Construct one worker-local analytic amplitude
  void InitializeProcess(LORENTZSCALAR &lts) override {
    PhotonMatrixElement = std::make_unique<AMP_yy_yy>(lts.model_cache->Tune().SM());
  }

  // Evaluate the charged-SM loop and apply both incoming photon sources
  double Amp2(LORENTZSCALAR &lts) override {
    const auto result = PhotonMatrixElement->Evaluate(lts, lts.alphaQCD, true);
    SetEvaluationStatus(result.status);
    if (!result.Valid()) { return 0.0; }
    return ApplyktEPA(result.amp2, lts);
  }
};

// -----------------------------------------------------------------------
// Inclusive processes
// -----------------------------------------------------------------------

class PROC_600_SOFT_EL : public MProc {
 public:
  PROC_600_SOFT_EL()
      : MProc("X", "EL", {"Elastic", "Eikonal CNI", "y, card", "Use with screening loop on", 600},
              MRegge::ProcessDefinitionFor(MReggeMode::Soft, "soft_elastic"),
              {ScreeningSpinBasis::ProtonHelicity, ProtonScreeningMode::Elastic, 0.25,
               ScreeningAmplitudeType::ElasticCNI},
              MParameterSet::Regge, {.example_decay = RootDecayMode::None}) {}
  // Initialize the selected Regge amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitRegge(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) { return abs2(Regge->ME2(lts, MReggeInclusive::EL)); }
};
class PROC_601_SOFT_SD : public MProc {
 public:
  PROC_601_SOFT_SD()
      : MProc("X", "SD", {"Single Diffractive", "Triple Pomeron", "card", "Use with screening loop on", 601},
              MRegge::ProcessDefinitionFor(MReggeMode::Soft, "soft_single_diffraction"),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::TripleRegge, 1.0,
               ScreeningAmplitudeType::GoodWalker},
              MParameterSet::Regge, {.example_decay = RootDecayMode::None}) {}
  // Initialize the selected Regge amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitRegge(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) { return abs2(Regge->ME2(lts, MReggeInclusive::SD)); }
};
class PROC_602_SOFT_DD : public MProc {
 public:
  PROC_602_SOFT_DD()
      : MProc("X", "DD", {"Double Diffractive", "Triple Pomeron", "card", "Use with screening loop on", 602},
              MRegge::ProcessDefinitionFor(MReggeMode::Soft, "soft_double_diffraction"),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::TripleRegge, 1.0,
               ScreeningAmplitudeType::GoodWalker},
              MParameterSet::Regge, {.example_decay = RootDecayMode::None}) {}
  // Initialize the selected Regge amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitRegge(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) { return abs2(Regge->ME2(lts, MReggeInclusive::DD)); }
};
class PROC_603_SOFT_ND : public MProc {
 public:
  PROC_603_SOFT_ND()
      : MProc("X", "ND",
              {"Non-Diffractive", "N-cut soft Pomerons", "card", "Unit normalized, see 'minbias' under bin/", 603},
              UnitAmplitudeProcessDefinitionFor("soft_nondiffractive"),
              {ScreeningSpinBasis::Scalar, ProtonScreeningMode::None, 1.0}, MParameterSet::None,
              {.example_decay = RootDecayMode::None}) {}
  // Initialize the unit soft non-diffractive amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { (void)lts; }

  virtual double Amp2(gra::LORENTZSCALAR &lts) { return 1.0; }
};

// -----------------------------------------------------------------------
// Pomeron-Pomeron processes
// -----------------------------------------------------------------------

// Map one common central-process mode to the independent TP implementation
inline MTensorPomeronMode TensorPomeronModeFor(MReggeMode mode) {
  switch (mode) {
    case MReggeMode::Resonance:
      return MTensorPomeronMode::Resonance;
    case MReggeMode::ContinuumTwoBody:
    case MReggeMode::ContinuumTwoFourSixBody:
      return MTensorPomeronMode::Continuum;
    case MReggeMode::ResonanceContinuumTwoBody:
      return MTensorPomeronMode::ResonanceContinuum;
    case MReggeMode::Generic:
    case MReggeMode::Soft:
      break;
  }
  throw std::invalid_argument("TensorPomeronModeFor: unsupported process mode");
}

// Select the process definition owned by one production model
inline std::shared_ptr<const amplitude::ProcessDefinition> CentralProcessDefinitionFor(
    ReggeProductionModel model, MReggeMode mode, const std::string &process_name) {
  if (model == ReggeProductionModel::TP) { return MTensorPomeron::ProcessDefinitionFor(TensorPomeronModeFor(mode)); }
  if (model == ReggeProductionModel::MP || model == ReggeProductionModel::XP || model == ReggeProductionModel::GP) {
    return MRegge::ProcessDefinitionFor(mode, process_name);
  }
  throw std::invalid_argument("CentralProcessDefinitionFor: production model is not selected");
}

// Select the immutable parameter family owned by one production model
inline MParameterSet CentralParameterSetFor(ReggeProductionModel model) {
  if (model == ReggeProductionModel::TP) { return MParameterSet::Tensor; }
  if (model == ReggeProductionModel::MP || model == ReggeProductionModel::XP || model == ReggeProductionModel::GP) {
    return MParameterSet::Regge;
  }
  throw std::invalid_argument("CentralParameterSetFor: production model is not selected");
}

// Compute capabilities from the registered production model and amplitude mode
inline ProcessInfo ReggeProcessInfo(ReggeProductionModel production_model, MReggeMode mode) {
  using enum nuclear::CollisionType;
  ProcessInfo info;
  info.model = production_model;
  info.resonance = mode == MReggeMode::Resonance || mode == MReggeMode::ResonanceContinuumTwoBody
                       ? ResonanceType::Required : ResonanceType::None;
  info.continuum = mode == MReggeMode::ContinuumTwoBody || mode == MReggeMode::ContinuumTwoFourSixBody ||
                   mode == MReggeMode::ResonanceContinuumTwoBody;
  info.jw_helicity_algebra = production_model != ReggeProductionModel::TP;
  info.beams.antiprotons = production_model != ReggeProductionModel::TP;
  if (mode == MReggeMode::Resonance) {
    info.beams.collisions = production_model == ReggeProductionModel::MP || production_model == ReggeProductionModel::XP
                                ? std::vector{PP, EP} : std::vector{PP, EP, PA, EA, AA};
  }
  return info;
}

class MReggeCentralProc : public MProc {
 public:
  // Construct one typed MP, XP, GP or TP central process
  MReggeCentralProc(const std::string &istate, const std::string &channel, const ProcessDescriptor &description,
                    ReggeProductionModel production_model, MReggeMode mode)
      : MProc(istate, channel, description,
              CentralProcessDefinitionFor(production_model, mode,
                                          istate + "_" + (channel == "RES+CON" ? "RES_CON" : channel)),
              {ScreeningSpinBasis::ProtonHelicity, ProtonScreeningMode::ForwardExcitation, 0.25,
               ScreeningAmplitudeType::GoodWalker},
              CentralParameterSetFor(production_model), ReggeProcessInfo(production_model, mode)),
        production_model(production_model),
        mode(mode) {}

  // Initialize only the runtime owned by the selected production model
  void InitializeProcess(gra::LORENTZSCALAR &lts) override {
    if (production_model == ReggeProductionModel::TP) {
      InitTensor(lts, TensorPomeronModeFor(mode));
      return;
    }
    regge_photo = UsesExternalPhotoBeams(lts);
    InitRegge(lts);
  }

  // Initialize the selected central Regge branching structures
  void InitializeBranching(gra::MProcessSetup &setup) override {
    if (production_model == ReggeProductionModel::TP) {
      MTensorPomeron::InitializeBranching(setup, TensorPomeronModeFor(mode));
      return;
    }
    MRegge::InitializeBranching(setup, production_model, mode);
    InitPhotoPlan(setup.lts);
  }

  // Resolve the complete Regge screening layout from the selected physics mode
  ScreeningMetadata Screening(const gra::MProcessSetup &setup) const override {
    ScreeningMetadata layout = MProc::Screening(setup);
    if (production_model == ReggeProductionModel::TP && setup.excitation != 0) {
      // Forward dissociation retains its explicit one-channel source model
      layout.amplitude_type = ScreeningAmplitudeType::Physical;
    }
    const bool asymmetric_photo = UsesExternalPhotoBeams(setup.lts);
    if (asymmetric_photo) {
      if (mode != MReggeMode::Resonance) {
        throw std::invalid_argument("Regge asymmetric photoproduction requires the RES channel");
      }
      layout.amplitude_type = ScreeningAmplitudeType::Physical;
    }
    const std::size_t central_particles = setup.lts.decaytree.size();
    const bool        scalar_pair       = !setup.lts.process.SPINGEN && central_particles == 2 &&
                             setup.lts.decaytree[0].p.spinX2 == 0 && setup.lts.decaytree[1].p.spinX2 == 0;
    const bool blind_hadronic =
        scalar_pair && production_model != ReggeProductionModel::GP && setup.lts.process.FORWARD_NOFLIP &&
        std::all_of(setup.lts.process.CONT_PRODUCTION.cbegin(), setup.lts.process.CONT_PRODUCTION.cend(),
                    [](const auto &channel) {
                      return channel.size() == 2 && channel[0] != PDG::PDG_gamma && channel[1] != PDG::PDG_gamma;
                    });
    if (mode == MReggeMode::ContinuumTwoFourSixBody && blind_hadronic) {
      // Spin-disabled scalar sources share one proton identity row
      layout.spin_basis              = ScreeningSpinBasis::ProtonIdentity;
      layout.amplitude_normalization = 1.0;
      layout.spin_rows               = 1;
      layout.forward_noflip          = true;
    }
    if (mode == MReggeMode::ContinuumTwoFourSixBody && (central_particles == 4 || central_particles == 6)) {
      layout.spin_basis              = ScreeningSpinBasis::ProtonIdentity;
      layout.amplitude_normalization = 1.0;
      layout.spin_rows               = 1;
      layout.forward_noflip          = true;
    }
    return layout;
  }

  // Evaluate one typed Regge central process
  double Amp2(gra::LORENTZSCALAR &lts) override {
    if (production_model == ReggeProductionModel::TP) { return TensorAmp2(lts); }
    if (regge_photo) { return EvalReggePhoto(*Regge, lts, production_model, lts.process.PHOTO_WIDTH_EE); }
    return Regge->Amp2(lts, production_model, mode);
  }

 private:
  // Initialize and validate the immutable VMD inputs before event sampling
  void InitPhotoPlan(LORENTZSCALAR &lts) {
    regge_photo = UsesExternalPhotoBeams(lts);
    if (!regge_photo) {
      lts.process.PHOTO_WIDTH_EE.clear();
      return;
    }
    if (mode != MReggeMode::Resonance) {
      throw std::invalid_argument("Regge asymmetric photoproduction requires the RES channel");
    }
    lts.process.PHOTO_WIDTH_EE = BuildPhotoPlan(lts);
  }

  // Compute whether MP, XP or GP needs external photo beam contraction
  bool UsesExternalPhotoBeams(const LORENTZSCALAR &lts) const {
    const bool regge_model = production_model == ReggeProductionModel::MP ||
                             production_model == ReggeProductionModel::XP ||
                             production_model == ReggeProductionModel::GP || production_model == ReggeProductionModel::TP;
    return regge_model && (nuclear::IsChargedLepton(lts.beam1.pdg) || nuclear::IsChargedLepton(lts.beam2.pdg) ||
                           nuclear::IsNuclearPDG(lts.beam1.pdg) || nuclear::IsNuclearPDG(lts.beam2.pdg));
  }

  // Add two coherent Tensor Good Walker amplitudes after exact layout validation
  static void MergeGoodWalker(std::optional<ProtonGoodWalkerAmplitude>       &target,
                              const std::optional<ProtonGoodWalkerAmplitude> &addition, const ScreeningMetadata &layout,
                              const std::string &context) {
    const bool expects_pair = layout.amplitude_type == ScreeningAmplitudeType::GoodWalker;
    if (!expects_pair) {
      if (target || addition) { throw AmplitudeFailure(context + ": unexpected Tensor Good Walker payload"); }
      return;
    }
    if (!target || !addition) { throw AmplitudeFailure(context + ": incomplete Tensor Good Walker sum"); }

    auto       &sum  = *target;
    const auto &term = *addition;
    if (sum.model == nullptr || term.model == nullptr || sum.model != term.model || sum.channel_count == 0 ||
        sum.channel_count != term.channel_count || sum.model->GoodWalker().ChannelCount() != sum.channel_count ||
        sum.components.empty() || sum.components.size() != term.components.size()) {
      throw AmplitudeFailure(context + ": incompatible Tensor Good Walker layouts");
    }

    const std::size_t pair_dimension = sum.model->GoodWalker().PairDimension();
    std::vector<bool> matched(sum.components.size(), false);
    for (const auto &added : term.components) {
      const auto found = std::find_if(sum.components.begin(), sum.components.end(), [&](const auto &component) {
        return component.coherence_group == added.coherence_group &&
               component.upper_sector == added.upper_sector && component.lower_sector == added.lower_sector;
      });
      if (found == sum.components.end()) {
        throw AmplitudeFailure(context + ": Tensor Good Walker component is missing");
      }
      const std::size_t index     = static_cast<std::size_t>(std::distance(sum.components.begin(), found));
      auto             &component = *found;
      if (matched[index] || layout.spin_rows == 0 || component.source.size_row() == 0 ||
          component.source.size_row() % layout.spin_rows != 0 || component.source.size_col() != pair_dimension ||
          component.source.size_row() != added.source.size_row() ||
          component.source.size_col() != added.source.size_col()) {
        throw AmplitudeFailure(context + ": incompatible Tensor Good Walker sources");
      }
      matched[index] = true;
      component.source += added.source;
    }
    if (std::find(matched.begin(), matched.end(), false) != matched.end()) {
      throw AmplitudeFailure(context + ": Tensor Good Walker component is unmatched");
    }
  }

  // Evaluate the independent TP resonance implementation
  double TensorResonanceAmp2(gra::LORENTZSCALAR &lts) {
    const double amp2 = UsesExternalPhotoBeams(lts) ? Tensor->PhotoProduction(lts, false) : Tensor->ME3(lts);
    RequireHelicityAmplitudes(lts.hamp, Label());
    return amp2;
  }

  // Evaluate the independent TP continuum implementation
  double TensorContinuumAmp2(gra::LORENTZSCALAR &lts) {
    return lts.decaytree[0].legs.empty() ? Tensor->ME4(lts, TensorContinuumMode::TensorPomeron) : Tensor->ME6(lts);
  }
  // Evaluate the coherent independent TP continuum plus resonance amplitude
  double TensorCoherentAmp2(gra::LORENTZSCALAR &lts) {

    // 1. Evaluate continuum matrix element -> helicity amplitudes to lts.hamp
    const double continuum_amp2 = Tensor->ME4(lts, TensorContinuumMode::TensorPomeron);
    if (lts.process.RESONANCES.empty()) { return continuum_amp2; }
    std::vector<std::complex<double>>              continuum_hamp = std::move(lts.hamp);
    const std::optional<ProtonGoodWalkerAmplitude> continuum_pair = std::move(lts.proton_good_walker);

    // 2. Evaluate resonance matrix elements -> helicity amplitudes to lts.hamp
    // We loop over resonances inside ME3
    Tensor->ME3(lts);
    if (lts.hamp.size() != continuum_hamp.size()) {
      throw AmplitudeFailure(
          "TP[RES+CON]: Continuum and resonance "
          "helicity bases are incompatible");
    }
    std::transform(lts.hamp.begin(), lts.hamp.end(), continuum_hamp.begin(), lts.hamp.begin(), std::plus<>{});
    MergeGoodWalker(lts.proton_good_walker, continuum_pair, lts.hamp.metadata, "TP[RES+CON]");

    return gra::SquaredNorm(lts.hamp) / 4.0;
  }

  // Dispatch one independent TP amplitude without exposing another process type
  double TensorAmp2(gra::LORENTZSCALAR &lts) {
    switch (mode) {
      case MReggeMode::Resonance:
        return TensorResonanceAmp2(lts);
      case MReggeMode::ContinuumTwoBody:
      case MReggeMode::ContinuumTwoFourSixBody:
        return TensorContinuumAmp2(lts);
      case MReggeMode::ResonanceContinuumTwoBody:
        return TensorCoherentAmp2(lts);
      case MReggeMode::Generic:
      case MReggeMode::Soft:
        break;
    }
    throw std::invalid_argument("TP process has an unsupported amplitude mode");
  }

  ReggeProductionModel production_model;
  MReggeMode           mode;
  bool                 regge_photo = false;
};

// Gauge-restored charged-meson photo and electroproduction process
class PROC_403_TP_PHOTO : public MProc {
 public:
  // Compute beam support for the configured photoproduction continuum
  bool SupportsCollision(nuclear::CollisionType collision, const nlohmann::json &general, const MPDG &pdg) const override;

  // Construct the explicit Tensor Pomeron photoproduction channel
  PROC_403_TP_PHOTO()
      : MProc("TP", "PHOTO",
              {"Charged-meson photo and electroproduction", "T-Pomeron", "card", "arXiv:2508.06334", 403},
              MTensorPomeron::ProcessDefinitionFor(MTensorPomeronMode::Photo),
              {ScreeningSpinBasis::ProtonHelicity, ProtonScreeningMode::ForwardExcitation, 0.25,
               ScreeningAmplitudeType::Physical},
              MParameterSet::Tensor,
              {.model = ReggeProductionModel::TP, .resonance = ResonanceType::Optional, .jw_helicity_algebra = false,
               .beams = nuclear::BeamSupport::Photo(false)}) {}

  // Initialize the shared Tensor runtime for photoproduction
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitTensor(lts, MTensorPomeronMode::Photo); }

  // Initialize meson and optional vector-resonance branching data
  void InitializeBranching(gra::MProcessSetup &setup) override {
    MTensorPomeron::InitializeBranching(setup, MTensorPomeronMode::Photo);
  }

  // Evaluate the coherent Drell-Soding and configured resonance amplitude
  double Amp2(gra::LORENTZSCALAR &lts) override {
    return Tensor->MEPhoto(lts);
  }
};

// -----------------------------------------------------------------------
// Durham QCD processes
// -----------------------------------------------------------------------

class MDurhamResonanceProc : public MProc {
 public:
  // Construct one Durham charmonium resonance channel
  MDurhamResonanceProc(const std::string &channel, const std::string &process, int pdg, int order)
      : MProc("gg", channel, {"QCD resonance " + channel, "Durham QCD", "gg", "PDF sensitive process", order},
              MDurham::ProcessDefinitionFor(MDurhamMode::Resonance, process),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 1.0}, MParameterSet::Durham),
        root_pdg(pdg) {}
  // Compute the selected charmonium PDG code
  int RootResonancePDG(const gra::LORENTZSCALAR &lts) const override {
    (void)lts;
    return root_pdg;
  }
  // Initialize the selected Durham resonance amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitDurham(lts); }
  // Evaluate the selected Durham resonance amplitude
  double Amp2(gra::LORENTZSCALAR &lts) override { return EvaluateDurhamQCD(lts, CHANNEL); }

 private:
  int root_pdg;
};

class MDurhamContinuumProc : public MProc {
 public:
  // Construct one Durham continuum family with its own final-state definition
  MDurhamContinuumProc(const std::string &channel, MDurhamMode mode, int order)
      : MProc("gg", channel, {"Durham continuum", "Durham QCD", "gg", "", order},
              MDurham::ProcessDefinitionFor(mode, "gg_" + channel),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 1.0},
              MParameterSet::Durham),
        mode(mode) {}
  // Initialize the selected Durham amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitDurham(lts); }

  // Evaluate the hard family selected and validated during initialization
  double Amp2(gra::LORENTZSCALAR &lts) override {
    return EvaluateDurhamQCD(lts, mode == MDurhamMode::MesonPair ? "MMbar" : "MG5");
  }

 private:
  MDurhamMode mode;
};

class PROC_704_DURHAM_FLUX : public MProc {
 public:
  PROC_704_DURHAM_FLUX()
      : MProc("gg", "FLUX", {"Durham flux with |A|^2 = 1", "Durham QCD", "gg", "Test process", 706},
              MDurham::ProcessDefinitionFor(MDurhamMode::Flux, "gg_flux"),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 1.0},
              MParameterSet::Durham, {.example_decay = RootDecayMode::None}) {}
  // Initialize the selected Durham amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitDurham(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) { return EvaluateDurhamQCD(lts, "FLUX"); }
};

// -----------------------------------------------------------------------
// Gamma-Gamma processes
// -----------------------------------------------------------------------

class PROC_040_QED_YY_LUX_EPA : public MProc {
 public:
  PROC_040_QED_YY_LUX_EPA()
      : MProc("yy_LUX", "EPA", {"Continuum l+l-, qqbar, monopolepair", "Collinear LUX-PDF", "yy", "MG5 amplitudes", 40},
              PhotonContinuumProcesses(),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25},
              MParameterSet::MonopolePair | MParameterSet::PhotonPDF, {.jw_helicity_algebra = false}) {}
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    double amp2 = GammaGammaCON(lts, false);
    return flux::ApplyLUXfluxes(amp2, lts);
  }
};

class PROC_030_QED_YY_DZ_EPA : public MProc {
 public:
  PROC_030_QED_YY_DZ_EPA()
      : MProc("yy_DZ", "EPA",
              {"Continuum l+l-, qqbar, monopolepair", "Collinear Drees-Zeppenfeld EPA", "yy", "MG5 amplitudes", 30},
              PhotonContinuumProcesses(), {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::None, 0.25},
              MParameterSet::MonopolePair, {.jw_helicity_algebra = false}) {}
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    double amp2 = GammaGammaCON(lts, false);
    return flux::ApplyDZfluxes(amp2, lts);
  }
};

class PROC_020_QED_YY_FLUX : public MProc {
 public:
  PROC_020_QED_YY_FLUX()
      : MProc("yy", "FLUX", {"kt-EPA flux with |A|^2 = 1", "kT-EPA", "yy", "Test process", 20},
              MGamma::ProcessDefinitionFor(MGammaMode::Flux),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 1.0}, MParameterSet::None,
              {.beams = nuclear::BeamSupport::EPA(), .example_decay = RootDecayMode::None}) {}
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    constexpr double amp2 = 1.0;
    lts.hamp.assign(1, std::complex<double>(1.0, 0.0));
    return ApplyktEPA(amp2, lts);
  }
};

class PROC_031_QED_YY_DZ_FLUX : public MProc {
 public:
  PROC_031_QED_YY_DZ_FLUX()
      : MProc("yy_DZ", "FLUX", {"DZ flux with |A|^2 = 1", "Collinear Drees-Zeppenfeld EPA", "yy", "Test process", 31},
              MGamma::ProcessDefinitionFor(MGammaMode::Flux),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::None, 1.0}, MParameterSet::None,
              {.example_decay = RootDecayMode::None}) {}
  // Initialize the selected gamma-gamma amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitGamma(lts); }

  // Apply collinear photon fluxes to a unit helicity amplitude
  virtual double Amp2(gra::LORENTZSCALAR &lts) {
    constexpr double amp2 = 1.0;
    lts.hamp.assign(1, std::complex<double>(1.0, 0.0));
    return flux::ApplyDZfluxes(amp2, lts);
  }
};

class PROC_050_PHOTO_YGG_Z : public MProc {
 public:
  // Construct the ygg to Z process descriptor
  PROC_050_PHOTO_YGG_Z()
      : MProc("ygg", "Z", {"Photoproduced Z", "Dipole/kT-factorized", "ygg", "", 50}, MPhotoZ::ProcessDefinitionFor(),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25,
               ScreeningAmplitudeType::Physical, 4},
              MParameterSet::PhotoZ,
              {.jw_helicity_algebra = false, .beams = nuclear::BeamSupport::Photo()}) {}
  // Initialize the exclusive Z photoproduction amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitPhotoZ(lts); }

  // Evaluate the photoproduced Z amplitude
  virtual double Amp2(gra::LORENTZSCALAR &lts) { return PhotoZ->Amp2(lts); }
};

class PROC_051_PHOTO_YGG_VM : public MProc {
 public:
  // Construct one ygg heavy-vector process descriptor
  PROC_051_PHOTO_YGG_VM(const std::string &channel, int display_order)
      : MProc("ygg", channel, {"Photoproduced vector meson", "JMRT", "ygg", "", display_order},
              MPhotoVM::ProcessDefinitionFor(channel),
              {ScreeningSpinBasis::ProtonIdentity, ProtonScreeningMode::ForwardExcitation, 0.25,
               ScreeningAmplitudeType::Physical, 4},
              MParameterSet::PhotoVM,
              {.jw_helicity_algebra = false, .beams = nuclear::BeamSupport::Photo()}) {}
  // Initialize the selected vector-meson photoproduction amplitude
  void InitializeProcess(gra::LORENTZSCALAR &lts) override { InitPhotoVM(lts); }

  // Evaluate the photoproduced heavy-vector amplitude
  virtual double Amp2(gra::LORENTZSCALAR &lts) { return PhotoVM->Amp2(lts); }
};

// Umbrella class
class MSubProc {
 public:
  MSubProc(const std::vector<std::string> &istate, const std::string &mc);
  void Initialize(const std::string &istate, const std::string &channel);

  // ---------------------------------------------------------------------
  MSubProc() {}
  ~MSubProc() {}

  // Processes are copied into per-thread process objects. Do not copy
  // the activated MProc instance: generated amplitudes own mutable buffers
  // and should be initialized independently in each thread. The copied
  // subprocess keeps only the selected process description and lazily rebuilds
  // pr on first GetBareAmplitude2 call
  MSubProc(MSubProc const &other)
      : ISTATE(other.ISTATE),
        CHANNEL(other.CHANNEL),
        LIPSDIM(other.LIPSDIM),
        ProcessRegistry(other.ProcessRegistry),
        random(nullptr),
        pr(nullptr),
        process_prepared(other.process_prepared),
        model_tune(other.model_tune),
        soft_model(other.soft_model) {}
  MSubProc(MSubProc &&) = default;

  MSubProc &operator=(MSubProc const &other) {
    if (this == &other) { return *this; }
    ISTATE           = other.ISTATE;
    CHANNEL          = other.CHANNEL;
    LIPSDIM          = other.LIPSDIM;
    ProcessRegistry  = other.ProcessRegistry;
    process_prepared = other.process_prepared;
    model_tune       = other.model_tune;
    soft_model       = other.soft_model;
    random           = nullptr;
    pr.reset();
    return *this;
  }
  MSubProc &operator=(MSubProc &&) = default;
  // ---------------------------------------------------------------------

  std::string                              ISTATE;           // "MP","XP","GP","TP","yy","gg" etc
  std::string                              CHANNEL;          // "RES","CON","RES+CON" etc
  unsigned int                             LIPSDIM = 0;      // Lorentz Invariant Phase Space Dimension
  std::map<std::string, ProcessDescriptor> ProcessRegistry;  // Process descriptions

  void BindRandom(gra::MRandom &rng);

  // Access the immutable SOFT model bound to this subprocess
  const SoftModelPtr &GetSoftModel() const noexcept { return soft_model; }

  // Access the immutable complete model tune bound to this subprocess
  const MModelTunePtr &GetModelTune() const noexcept { return model_tune; }

  // Initialize the selected physical amplitude for this process copy
  void InitializeAmplitude(gra::MProcessSetup &setup);

  bool SampleColorFlow(gra::LORENTZSCALAR &lts);
  // Evaluate the active process at one screening-loop kinematic point
  double ScreeningAmp2(gra::LORENTZSCALAR &lts) {
    return lts.epa_hard.Ready() ? EPAHardAmp2(lts) : GetBareAmplitude2(lts);
  }
  // Normalize external-helicity amplitudes through their process screening definition
  AmplitudeWeights NormalizeAmplitude(const MHelicityAmplitudes &amplitude) {
    return ActiveProcess().NormalizeAmplitude(amplitude);
  }
  // Normalize projected external-helicity norms through their process screening definition
  AmplitudeWeights NormalizeAmplitude(const std::vector<double> &helicity_norm, const ScreeningMetadata &metadata) {
    return ActiveProcess().NormalizeProjectedAmplitude(helicity_norm, metadata);
  }
  // Normalize a proton-screened amplitude through its owning process
  AmplitudeWeights NormalizeScreenedAmplitude(const MHelicityAmplitudes &amplitude, const bool dense_proton_spin) {
    return ActiveProcess().NormalizeScreenedAmplitude(amplitude, dense_proton_spin);
  }
  // Evaluate shifted photon sources with the node-local EPA phase space
  double EPAHardAmp2(gra::LORENTZSCALAR &lts) { return pr == nullptr ? 0.0 : pr->EPAHardAmp2(lts); }
  // Prepare reusable amplitude state for one phase-space point
  void PrepareBareAmplitude(gra::LORENTZSCALAR &lts);
  // Resolve and store the active amplitude decay structure
  DecayStructure DecayStructureFor(gra::LORENTZSCALAR &lts);
  // Match the active amplitude process against one decay tree
  std::optional<amplitude::Process> MatchProcess(const std::vector<MDecayBranch> &decaytree);
  // Compute an optional amplitude constraint for one generic-parton mode
  std::optional<bool> AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree);
  // Compute processes for the selected process and activate it when needed
  std::vector<amplitude::Process> Processes();
  double                          GetBareAmplitude2(gra::LORENTZSCALAR &lts);
  // Evaluate with previously prepared reusable amplitude state
  double GetPreparedBareAmplitude2(gra::LORENTZSCALAR &lts);
  // Compute the status of the most recent event-local matrix-element evaluation
  mg5helas::EvaluationStatus EvaluationStatus() const;
  int                        RootResonancePDG(const gra::LORENTZSCALAR &lts) const;
  bool                       ProcessExist(const std::string &process) const;
  // Compute concrete native and generated final states for this phase space
  std::vector<std::vector<std::string>> SupportedFinalStateRows() const;
  std::vector<std::shared_ptr<MProc>>   CreateAllProcesses() const;
  // Compute the selected physical process for initialization and queries
  std::shared_ptr<MProc> SelectedProcess() const;
  std::vector<std::string>              GetProcessDescriptor(const std::string &process) const;
  // Compute whether the selected process uses Jacob-Wick decay amplitudes
  bool UsesJWHelicityAlgebra() const { return SelectedProcess()->Info().jw_helicity_algebra; }

 private:
  struct InitialStateSelection {
    const std::string &value;
  };
  struct PhaseSpaceSelection {
    const std::string &value;
  };

  void ConstructDescriptions(InitialStateSelection istate, PhaseSpaceSelection mc);
  void ActivateProcess();

  // Resolve the selected process before applying its physical normalization
  MProc &ActiveProcess() {
    if (pr == nullptr || pr->ISTATE != ISTATE || pr->CHANNEL != CHANNEL) { ActivateProcess(); }
    return *pr;
  }

  // Require master process initialization before a worker amplitude is used
  void RequirePreparedProcess() const {
    if (!process_prepared) { throw std::logic_error("MSubProc: physical process state is not initialized"); }
  }

  gra::MRandom *random = nullptr;  // Non-owning process-local random number generator

  // Active mutable amplitude runtime owned by one process copy
  std::shared_ptr<MProc> pr               = nullptr;
  bool                   process_prepared = false;
  MModelTunePtr          model_tune;
  SoftModelPtr           soft_model;
};

}  // namespace gra

#endif
