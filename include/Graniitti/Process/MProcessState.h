// Process configuration and worker-local mutable state
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPROCESSSTATE_H
#define MPROCESSSTATE_H

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Tensor algebra
#include "FTensor.hpp"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MFinal.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/PDF/MSudakov.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Photon/MRadiative.h"
#include "Graniitti/QCD/MPartonFlavour.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Regge/MReggeProductionModel.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Sampling/MCW.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MQMetrics.h"
#include "Graniitti/Spin/MHELMatrix.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {

namespace mg5helas {

// Store one fixed-basis complex hard-amplitude component
struct HelicityComponent {
  std::array<int, 2>                incoming = {0, 0};
  std::vector<int>                  outgoing;
  std::size_t                       color = 0;
  std::complex<double>              value = 0.0;
  std::vector<std::complex<double>> flow_values;
};

// Store one external-leg color and anticolor pair
struct ColorFlowLeg {
  int color     = 0;
  int anticolor = 0;
};

using ExternalColorFlow = std::vector<ColorFlowLeg>;

// Pair one color amplitude row with its external color assignment
struct HardColorFlow {
  std::vector<std::complex<double>> amplitudes;
  ExternalColorFlow                 external;
  std::optional<double>             screened_weight;
  std::optional<std::size_t>        channel;

  // Compute the explicit screened weight or the coherent row norm
  double Weight() const { return screened_weight.has_value() ? *screened_weight : gra::SquaredNorm(amplitudes); }
};

// Sample one finite hard color-flow weight, optionally within one subprocess channel
std::optional<std::size_t> SelectHardColorFlow(const std::vector<HardColorFlow> &flow, MRandom &random,
                                               std::optional<std::size_t> channel = std::nullopt);

// Store the transfer-independent hard rest-frame transformation
struct EPAHardFrame {
  M4Vec                hard;
  double               mass = 0.0;
  std::array<M4Vec, 2> beam;
  std::array<M4Vec, 2> incoming;
  MMatrix<double>      rotation;
  bool                 valid = false;
};

// Store one color projection of the same hard photon tensor
struct EPAHardColor {
  std::vector<HelicityComponent> components;
  ExternalColorFlow              external;
};

// Own the immutable hard photon tensor prepared by the current Born event
struct EPAHardTensor {
  EPAHardFrame                   frame;
  std::vector<HelicityComponent> amplitude;
  std::vector<EPAHardColor>      color;
  double                         normalization = 1.0;

  // Clear the tensor before preparing a new Born event
  void Clear() {
    frame = {};
    amplitude.clear();
    color.clear();
    normalization = 1.0;
  }

  // Compute whether the Born event prepared a contractible tensor
  bool Ready() const { return frame.valid && !amplitude.empty() && normalization > 0.0; }
};

}  // namespace mg5helas

// Root decay operator semantics selected by the process arrow
enum class RootDecayMode { None, Physical, Isolated };

// Central phase-space map used by proposal-density transformations
enum class CentralPhaseSpaceMode { Unknown, Factorized, Central, Collinear, HardDiffraction };

// Forward proton spin content stored in the amplitude array
enum class ScreeningSpinBasis { Unset, Scalar, ProtonIdentity, ProtonHelicity };

// Good Walker state content of the forward process
enum class ProtonScreeningMode { Elastic, ForwardExcitation, TripleRegge, None };

// Select the amplitude representation consumed by screening
enum class ScreeningAmplitudeType { Physical, ElasticCNI, GoodWalker };

// One prepared hard proton transition with pair order (--,-+,+-,++)
// Negative proton helicity is bit 0 and positive proton helicity is bit 1
struct ScreeningSpinTransition {
  std::uint8_t initial      = 0;
  std::uint8_t intermediate = 0;
  std::uint8_t source_row   = 0;
};

// Fixed screening information selected by one physical process class
struct ScreeningMetadata {
  // Construct one complete process-level screening definition
  constexpr ScreeningMetadata(ScreeningSpinBasis     spin_basis_in              = ScreeningSpinBasis::Unset,
                              ProtonScreeningMode    proton_mode_in             = ProtonScreeningMode::None,
                              double                 amplitude_normalization_in = 1.0,
                              ScreeningAmplitudeType amplitude_type_in          = ScreeningAmplitudeType::Physical,
                              std::size_t spin_rows_in = 1, bool forward_noflip_in = false)
      : spin_basis(spin_basis_in),
        proton_mode(proton_mode_in),
        amplitude_normalization(amplitude_normalization_in),
        spin_rows(spin_rows_in),
        forward_noflip(forward_noflip_in),
        amplitude_type(amplitude_type_in) {}

  ScreeningSpinBasis                      spin_basis              = ScreeningSpinBasis::Unset;
  ProtonScreeningMode                     proton_mode             = ProtonScreeningMode::None;
  double                                  amplitude_normalization = 1.0;
  std::size_t                             spin_rows               = 1;
  bool                                    forward_noflip          = false;
  std::array<ScreeningSpinTransition, 16> spin_transition;
  std::size_t                             spin_transition_count = 0;
  ScreeningAmplitudeType                  amplitude_type;

  // Prepare the compact hard spin-row map once during process initialization
  void PrepareSpinTransitions() {
    spin_transition_count = 0;
    const auto add = [this](const std::size_t initial, const std::size_t intermediate, const std::size_t source_row) {
      spin_transition[spin_transition_count++] = {static_cast<std::uint8_t>(initial),
                                                  static_cast<std::uint8_t>(intermediate),
                                                  static_cast<std::uint8_t>(source_row)};
    };
    if (spin_basis == ScreeningSpinBasis::ProtonIdentity) {
      if (spin_rows != 1 && spin_rows != 4) {
        throw std::invalid_argument("ScreeningMetadata: proton identity basis needs 1 or 4 rows");
      }
      for (std::size_t initial = 0; initial < 4; ++initial) { add(initial, initial, spin_rows == 1 ? 0 : initial); }
    } else if (spin_basis == ScreeningSpinBasis::ProtonHelicity) {
      if (spin_rows != (forward_noflip ? 4 : 16)) {
        throw std::invalid_argument("ScreeningMetadata: proton helicity row count is inconsistent");
      }
      for (std::size_t initial = 0; initial < 4; ++initial) {
        const std::size_t count = forward_noflip ? 1 : 4;
        for (std::size_t offset = 0; offset < count; ++offset) {
          const std::size_t intermediate = forward_noflip ? initial : offset;
          const std::size_t source_row =
              forward_noflip ? initial : spin::PairHelicityTransitionIndex(initial, intermediate);
          add(initial, intermediate, source_row);
        }
      }
    }
  }
};

// Event-local resolved EPA and photonuclear amplitude rows
struct ScreeningLayout {
  nuclear::ScreenType                        nuclear_type          = nuclear::ScreenType::Scalar;
  bool                                       epa_sector_resolved   = false;
  std::size_t                                epa_rows_per_sector   = 0;
  std::array<std::uint8_t, 2>                epa_sector_count      = {0, 0};
  std::array<std::array<std::uint8_t, 2>, 2> epa_sector_type       = {std::array<std::uint8_t, 2>{255, 255},
                                                                      std::array<std::uint8_t, 2>{255, 255}};
  bool                                       photo_sector_resolved = false;
  std::size_t                                photo_channel_count   = 0;
  std::array<nuclear::PhotoChannel, 8>       photo_channel{};

  // Compare the complete discrete amplitude layout used by screening
  constexpr bool operator==(const ScreeningLayout &other) const noexcept = default;

  // Clear all resolved source-sector rows before amplitude evaluation
  void Clear() noexcept { *this = {}; }

  // Prepare the physical upper-then-lower photoproduction channel order
  void ConfigurePhotoChannels(const nuclear::MUPC &upc) {
    nuclear_type          = nuclear::ScreenType::Photo;
    photo_sector_resolved = false;
    photo_channel_count   = 0;
    photo_channel         = {};
    const auto channel    = nuclear::PhotoChannels(upc);
    if (channel.size() > photo_channel.size()) {
      throw std::invalid_argument("ScreeningLayout::ConfigurePhotoChannels: too many channels");
    }
    std::copy(channel.begin(), channel.end(), photo_channel.begin());
    photo_channel_count   = channel.size();
    photo_sector_resolved = photo_channel_count > 0;
  }
};

// Contiguous helicity amplitudes with process-defined screening information
class MHelicityAmplitudes : public std::vector<std::complex<double>> {
 public:
  using std::vector<std::complex<double>>::vector;
  using std::vector<std::complex<double>>::operator=;

  // Store immutable process screening information without changing amplitudes
  void Configure(const ScreeningMetadata &input) {
    metadata = input;
    metadata.PrepareSpinTransitions();
  }

  ScreeningMetadata metadata;
  ScreeningLayout   layout;
};

// Store normalized external-helicity weights and their physical sum
struct AmplitudeWeights {
  std::vector<double> helicity_weight;
  double              total_weight = 0.0;
};

// Own photonuclear amplitude data produced at one screening node
struct PhotoScreenState {
  std::array<std::optional<nuclear::PhotoCurrent>, 2> current;
  std::vector<nuclear::PhotoTerm>                     term;

  // Clear node-local photonuclear currents and direction terms
  void Clear() noexcept {
    for (auto &value : current) { value.reset(); }
    term.clear();
  }
};

// Own Durham quantities reused between projections at one screening node
struct DurhamScreenState {
  bool  jet_cuts_cached = false;
  bool  jet_cuts_pass   = false;
  bool  hard_cached     = false;
  M4Vec hard_k1;
  M4Vec hard_k2;

  // Start one shifted hard evaluation while retaining central-system jet cuts
  void BeginNode() noexcept {
    hard_cached = false;
    hard_k1     = {};
    hard_k2     = {};
  }

  // Clear every node-local Durham result
  void Clear() noexcept { *this = {}; }
};

// Own the complete process-local state of one screening convolution
struct ScreeningEventState {
  bool              active = false;
  PhotoScreenState  photo;
  DurhamScreenState durham;

  // Start one Born or shifted amplitude evaluation with empty node-local caches
  void BeginNode() noexcept {
    photo.Clear();
    durham.BeginNode();
  }

  // Clear every process cache when a screening loop ends
  void Clear() noexcept {
    active = false;
    photo.Clear();
    durham.Clear();
  }
};

// Orthogonal forward sector carried by one Good Walker pair component
enum class ProtonGoodWalkerSector { Elastic, TripleResolved, TripleInclusive, InelasticEPA, PhotoDiss };

// Compute the physical Good Walker basis carried by one forward sector
inline GoodWalkerFinalBasis SectorFinalBasis(const ProtonGoodWalkerSector sector) {
  switch (sector) {
    case ProtonGoodWalkerSector::Elastic:
      return GoodWalkerFinalBasis::Proton;
    case ProtonGoodWalkerSector::TripleResolved:
      return GoodWalkerFinalBasis::Excited;
    case ProtonGoodWalkerSector::TripleInclusive:
    case ProtonGoodWalkerSector::InelasticEPA:
    case ProtonGoodWalkerSector::PhotoDiss:
      return GoodWalkerFinalBasis::Complete;
  }
  throw std::invalid_argument("SectorFinalBasis: unknown sector");
}

// One coherent rectangular source in the Good Walker pair space
struct ProtonGoodWalkerComponent {
  std::size_t                   coherence_group         = 0;
  ProtonGoodWalkerSector        upper_sector            = ProtonGoodWalkerSector::Elastic;
  ProtonGoodWalkerSector        lower_sector            = ProtonGoodWalkerSector::Elastic;
  MMatrix<std::complex<double>> source;
  std::optional<std::size_t> spin_index = std::nullopt;
};

// Event-local arbitrary-channel Good Walker amplitude
struct ProtonGoodWalkerAmplitude {
  SoftModelPtr                           model;
  std::size_t                            channel_count = 0;
  std::vector<ProtonGoodWalkerComponent> components;
};

// Select the physical decay construction independently of the sampling proposal
// None omits physical decay amplitudes, Full owns the complete stable-particle amplitude
// JacobWick types distinguish coherent history sums from unsymmetrized or spin-blind cascades
enum class DecayType { None, JacobWickIncoherent, JacobWickCoherent, Full };

// Declare amplitude ownership and support for isolated production
struct DecayStructure {
  DecayType type = DecayType::None;
  bool allows_isolated_resonance = true;

  // Compute whether the physical amplitude includes all indistinguishable histories
  constexpr bool Coherent() const { return type == DecayType::JacobWickCoherent || type == DecayType::Full; }

  // Compare physical decay declarations
  constexpr bool operator==(const DecayStructure &) const = default;
};

// Select the forward beam helicity source convention
enum class ForwardVertexMode { HelicityResidue, UnitResidue };

// Key one local continuum pole pair by ordered final and exchange PDGs
using ReggeContinuumPoleKey = std::array<int, 4>;

// Store the two prepared pole operators of one local continuum pair
struct ReggeContinuumPole {
  std::array<spin::PoleResidue, 2>       pole_operator;
  std::array<HELMatrix, 2>               gp_vertex;
};

// Store prepared process data copied from the initialized master process
struct MProcessModelState {
  // Resonance and root-decay amplitude data
  std::map<std::string, PARAM_RES> RESONANCES;
  PARAM_RES                        ROOT_RES;
  bool                             ROOT_RES_ACTIVE = false;

  // Resolved subprocess root identity and decay semantics
  int           root_resonance_pdg = 0;
  RootDecayMode root_decay_mode    = RootDecayMode::None;

  // Forward excitation prescription resolved during process initialization
  DissociationType DISSOCIATION = DissociationType::Soft;
  DissociationType PHOTO_DISSOCIATION = DissociationType::Soft;

  // Spin and production-vertex steering
  std::string       MP_FRAME                = "null";
  bool              SPINGEN                 = true;
  bool              SPINDEC                 = true;
  bool              QMETRICS             = false;
  bool              DECAY_BARRIER           = true;
  bool              DERIVATIVE_FACTOR       = false;
  bool              FORWARD_NOFLIP          = true;
  ForwardVertexMode FORWARD_VERTEX          = ForwardVertexMode::HelicityResidue;
  std::string       PHOTON_VERTEX           = "EPA";
  int               MMAX                    = 0;
  std::string       TU_SIGN                 = "auto";
  ReggeProductionModel REGGE_MODEL          = ReggeProductionModel::None;

  // Electronic widths prepared for nuclear vector meson photoproduction
  std::map<int, double> PHOTO_WIDTH_EE;

  // Prepared Regge continuum production channels
  std::vector<std::vector<int>>          CONT_PRODUCTION;
  std::vector<std::vector<MDecayBranch>> CONT_PRODUCTIONTREE;
  // Canonical physical pole operators used by MP and XP continuum
  std::vector<std::vector<spin::PoleResidue>> CONTINUUM_POLE;
  // Analytic GP continuum subvertices
  std::vector<std::vector<HELMatrix>>                 CONTINUUM_GP;
  std::vector<double>                                 CONT_TU_SIGN;
  std::vector<std::vector<int>>                       CONT_LADDER_PERMUTATIONS;
  std::map<ReggeContinuumPoleKey, ReggeContinuumPole> CONT_LADDER_POLE;
  std::vector<std::vector<int>>                       MULTIREGGE_TOPOLOGIES;

  // Prepared Tensor Pomeron particle content
  bool               TENSOR_MODEL_READY = false;
  std::vector<int>   TENSOR_PSEUDOSCALAR_PDGS;
  std::vector<int>   TENSOR_BARYON_PDGS;
  std::map<int, int> TENSOR_VECTOR_DECAY_PDGS;
};

// Identify one worker-local event whose invariant data may be cached
class MEventCacheScope {
 public:
  MEventCacheScope() = default;

  // Copy no active event identity into another process state
  MEventCacheScope(const MEventCacheScope &) noexcept {}

  // Move no active event identity into another process state
  MEventCacheScope(MEventCacheScope &&) noexcept {}

  // Reset the destination instead of copying an event identity
  MEventCacheScope &operator=(const MEventCacheScope &other) noexcept {
    if (this != &other) {
      active = false;
      id     = 0;
    }
    return *this;
  }

  // Reset the destination instead of moving an event identity
  MEventCacheScope &operator=(MEventCacheScope &&other) noexcept {
    if (this != &other) {
      active = false;
      id     = 0;
    }
    return *this;
  }

  // Start one new event generation and invalidate every previous binding
  void Begin() noexcept {
    active = true;
    id     = id == std::numeric_limits<std::uint64_t>::max() ? 1 : id + 1;
  }

  // Invalidate the current event after a failed or rejected calculation
  void Abort() noexcept { active = false; }

  // Check whether one event generation is currently active
  bool Active() const noexcept { return active; }

  // Compute the worker-local identity of the current event generation
  std::uint64_t Id() const noexcept { return id; }

  // Reject cache access outside an established event generation
  void Require(const std::string &context) const {
    if (!active) { throw std::logic_error(context + ": no active event cache scope"); }
  }

 private:
  bool          active = false;
  std::uint64_t id     = 0;
};

// Store typed values which are valid only within one explicit event scope
template <typename Key, typename Value>
class MEventCache {
 public:
  MEventCache() = default;

  // Copy no event binding because a scope has object-local identity
  MEventCache(const MEventCache &) noexcept {}

  // Move no event binding because a scope has object-local identity
  MEventCache(MEventCache &&) noexcept {}

  // Clear event data when assigning across scope-owning objects
  MEventCache &operator=(const MEventCache &other) noexcept {
    if (this != &other) { Clear(); }
    return *this;
  }

  // Clear event data when moving across scope-owning objects
  MEventCache &operator=(MEventCache &&other) noexcept {
    if (this != &other) { Clear(); }
    return *this;
  }

  // Insert one value exactly once under the active event generation
  const Value &Store(const MEventCacheScope &scope, const Key &key, Value value, const std::string &context) {
    scope.Require(context);
    Bind(scope);
    const auto [entry, inserted] = values.emplace(key, std::move(value));
    if (!inserted) { throw std::logic_error(context + ": duplicate event cache key"); }
    return entry->second;
  }

  // Compute one value only when the cache belongs to the active event generation
  const Value &Get(const MEventCacheScope &scope, const Key &key, const std::string &context) const {
    scope.Require(context);
    if (owner != &scope || id != scope.Id()) { throw std::logic_error(context + ": event cache generation mismatch"); }
    const auto entry = values.find(key);
    if (entry == values.end()) { throw std::logic_error(context + ": event cache key is absent"); }
    return entry->second;
  }

  // Clear all values and remove their event-generation binding
  void Clear() noexcept {
    owner = nullptr;
    id    = 0;
    values.clear();
  }

 private:
  // Bind an empty cache to one generation or replace an older generation
  void Bind(const MEventCacheScope &scope) {
    if (owner == &scope && id == scope.Id()) { return; }
    values.clear();
    owner = &scope;
    id    = scope.Id();
  }

  const MEventCacheScope *owner = nullptr;
  std::uint64_t           id    = 0;
  std::map<Key, Value>    values;
};

// Store one central resonance decay result for a complete Born event
struct MReggeDecayState {
  double                        q_ratio         = 1.0;
  double                        running_profile = 1.0;
  bool                          valid           = false;
  MMatrix<std::complex<double>> matrix;
  MMatrix<std::complex<double>> pair;
};

// Classify the covariant Tensor resonance amplitude
enum class TensorResonanceType { Scalar, Pseudoscalar, Vector, AxialVector, Tensor, Spin3 };

// Store one central Tensor resonance decay block for a complete Born event
struct TensorDecayState {
  TensorResonanceType                                       type       = TensorResonanceType::Scalar;
  std::complex<double>                                      propagator = 0.0;
  std::vector<std::complex<double>>                         scalar;
  FTensor::Tensor1<std::complex<double>, 4>                 vector;
  FTensor::Tensor3<double, 4, 4, 4>                         spin3;
  MMatrix<std::complex<double>>                             axial;
  std::vector<FTensor::Tensor2<std::complex<double>, 4, 4>> tensor;
};

// Store mutable amplitude data for one sampled phase-space point
struct MAmplitudeEventState {
  // Start one central event and invalidate all previous typed cache entries
  void BeginCentral() noexcept {
    central.Begin();
    continuum_decay.Clear();
    regge_decay.Clear();
    tensor_decay.Clear();
    durham_decay.Clear();
  }

  // Invalidate every central cache after an event-level failure
  void AbortCentral() noexcept {
    central.Abort();
    continuum_decay.Clear();
    regge_decay.Clear();
    tensor_decay.Clear();
    durham_decay.Clear();
  }

  bool                                                    DECAY_SYM = false;
  MEventCacheScope                                        central;
  MEventCache<std::string, MMatrix<std::complex<double>>> continuum_decay;
  MEventCache<std::string, MReggeDecayState>              regge_decay;
  MEventCache<std::string, TensorDecayState>              tensor_decay;
  MEventCache<std::string, MMatrix<std::complex<double>>> durham_decay;
};

// Lorentz scalars and other common kinematic variables
class LORENTZSCALAR {
 public:
  LORENTZSCALAR() {}
  ~LORENTZSCALAR() {}

  // Event-local screening state owned by the active subprocess
  ScreeningEventState screening;

  // --------------------------------------------------------------------
  // Particle Database
  MPDG PDG;

  // Process-local physical species selected for generic final-state partons
  MPartonFlavour final_state_partons;

  // Process-local external quark masses required by the selected amplitude
  std::map<int, double> final_state_parton_masses;

  // Preserve explicit particle requests for generated model initialization
  std::map<int, double> particle_mass_overrides;
  std::map<int, double> particle_width_overrides;

  // Parsed alias tree restored before each physical parton proposal
  std::vector<MDecayBranch> generic_parton_decaytree;

  // Physical stable-PDG modes accepted by the active amplitude definition
  std::vector<std::vector<int>> final_state_parton_modes;

  // Inverse discrete probability 1/q(mode) of the physical parton proposal
  double parton_proposal_weight = 1.0;

  // --------------------------------------------------------------------
  // PDF access
  //
  // These caches are shared_ptr because LORENTZSCALAR is copied into
  // temporary work objects and thread-local process copies. LHAPDF members and
  // Sudakov/Shuvaev tables are acquired from process wide stores and treated as
  // read only
  //
  // Initialization is serialized at the call sites. After initialization the
  // cached objects must be treated as read only. shared_ptr makes the lifetime
  // safe across copies. These objects must not carry event-local mutable state
  std::shared_ptr<const MSudakov> GlobalSudakovPtr = nullptr;

  // Expensive immutable tables shared only by workers of this run
  MModelCachePtr model_cache;

  // Immutable heavy-ion UPC physics shared by worker-local process copies
  std::shared_ptr<const nuclear::MUPC> upc_model;

  // Nuclear fluctuations sampled once for the active event evaluation
  std::shared_ptr<const nuclear::MUPC> upc_event;

  // Orthogonal EMD channel sampled before the hard kinematics and fixed throughout screening
  std::shared_ptr<const nuclear::ExcitationChannel> upc_excitation;

  std::string                        LHAPDFSET    = "null";
  std::shared_ptr<const LHAPDF::PDF> GlobalPdfPtr = nullptr;

  // Explicit hard-diffraction variables for f_i/p^D(xi,beta,t,Q2) bookkeeping
  // Side 1 / 2
  bool   hard_diff1  = false;  // Side 1 uses a diffractive parton density
  bool   hard_diff2  = false;  //
  double diff_xi1    = 0.0;    // Leading hadron loss xi
  double diff_xi2    = 0.0;    //
  double diff_beta1  = 0.0;    // Parton fraction beta in the side 1 exchange
  double diff_beta2  = 0.0;    //
  double diff_t1     = 0.0;    // Leading hadron transfer t in GeV^2 on side 1
  double diff_t2     = 0.0;    //
  double diff_phi1   = 0.0;    // Leading hadron azimuth in radians on side 1
  double diff_phi2   = 0.0;    //
  double diff_xhard1 = 0.0;    // Hard fraction xi1 beta1 on side 1
  double diff_xhard2 = 0.0;    //
  // --------------------------------------------------------------------

  // Prepared amplitude state initialized before worker processes are copied
  MProcessModelState process;

  // Mutable amplitude state rebuilt for each sampled phase-space point
  MAmplitudeEventState amplitude;

  // Bose-Einstein decay symmetrizations
  std::string                           decay_symmetry_topology_key;
  std::vector<std::vector<std::size_t>> decay_symmetry_assignments;
  std::vector<double>                   decay_symmetry_statistics_signs;
  std::size_t                           decay_symmetry_proposal_index       = 0;
  bool                                  decay_symmetry_proposal_active      = false;
  double                                decay_symmetry_proposal_phase_space = 0.0;

  // Central phase space support used by crossed proposal channels
  CentralPhaseSpaceMode central_phase_space_mode                    = CentralPhaseSpaceMode::Unknown;
  double                central_phase_space_rap_min_cm              = 0.0;
  double                central_phase_space_rap_max_cm              = 0.0;
  double                central_phase_space_mass_cut_min            = 0.0;
  double                central_phase_space_mass_max                = 0.0;
  double                central_phase_space_mass_margin             = 0.0;
  double                central_phase_space_generated_jacobian = 0.0;
  double                central_phase_space_transverse_radius2 = 0.0;

  // Helicity amplitudes returned by amplitude functions,
  // used in screening loop etc
  MHelicityAmplitudes hamp;

  // Worker-local spin integration sums and event scratch
  spin::MQMetrics qmetrics;

  // Transfer-independent photon tensor and hard-process color-flow output
  mg5helas::EPAHardTensor              epa_hard;
  std::vector<mg5helas::HardColorFlow> hard_color_flows;

  // Pair-space hard sources used by arbitrary-channel screening
  std::optional<ProtonGoodWalkerAmplitude> proton_good_walker;

  // Central system decaytree
  std::vector<MDecayBranch> decaytree;

  // Forward system decaytree
  MDecayBranch decayforward1;
  MDecayBranch decayforward2;

  // Integral weight container (central system phase space)
  gra::kinematics::MCW DW;

  // Active in the factorized phase space product (by default, true)
  // This is controlled by the spesific amplitudes
  bool PS_active = true;

  // State which decay factors are included in the matrix element
  DecayStructure decay_structure;

  // Sum containers
  gra::kinematics::MCWSUM DW_sum;
  gra::kinematics::MCWSUM DW_sum_exact;

  // Initial states
  gra::MParticle beam1;
  gra::MParticle beam2;

  // Four momenta of initial and final states
  M4Vec              pbeam1;
  M4Vec              pbeam2;
  std::vector<M4Vec> pfinal = std::vector<M4Vec>(11);
  std::vector<M4Vec> pfinal_orig;

  // Basic Lorentz scalars
  double s      = 0.0;
  double t      = 0.0;
  double u      = 0.0;
  double sqrt_s = 0.0;

  // Sub-energies^2
  double s1 = 0.0;
  double s2 = 0.0;

  // 4-momentum transfer squared
  double t1 = 0.0;
  double t2 = 0.0;

  // Sub-Mandelstam variables
  double s_hat = 0.0;
  double t_hat = 0.0;
  double u_hat = 0.0;

  // Central system
  double m2 = 0.0;
  double Y  = 0.0;
  double Pt = 0.0;

  M4Vec d0_in_X;  // 2-body decay (direct) products in (central system) X-frame
  M4Vec d1_in_X;

  // Slot 0 stores the central sum in addition to
  // 2 forward and 8 direct central particles
  double ss[11][11]    = {{0.0}};
  double tt_1[11]      = {0.0};
  double tt_2[11]      = {0.0};
  double tt_xy[11][11] = {{0.0}};

  // Longitudinal momentum fractions
  double x1      = 0.0;    // Side 1 exchange fraction used in the hard amplitude or density
  double x2      = 0.0;    //
  double xi1     = 0.0;    // Tagged forward momentum loss on side 1
  double xi2     = 0.0;    //
  bool   has_xi1 = false;  // True when xi1 describes a physical forward state
  bool   has_xi2 = false;  //
  double xbj1    = 0.0;    // DIS Bjorken variable -Q1^2/(2 p1 dot q1) on side 1
  double xbj2    = 0.0;    //

  // Forward proton excitation tag
  bool excite1 = false;
  bool excite2 = false;
  // Forward external-system invariant mass squares before decay
  std::array<double, 2> forward_mass2 = {-1.0, -1.0};
  std::array<double, 2> forward_emd{};  // Independent EMD excitation in GeV

  // Propagators from proton1 (up) and proton2 (down)
  M4Vec q1;
  M4Vec q2;
  M4Vec q1_in_X;  // In (central system) X-frame
  M4Vec q2_in_X;

  // Propagator pt
  double qt1 = 0.0;
  double qt2 = 0.0;

  // ** HepMC output **
  int    id1     = 0;    // Incoming parton or photon PDG id on beam side 1
  int    id2     = 0;    // Incoming parton or photon PDG id on beam side 2
  double pdf_xf1 = 0.0;  // HepMC PDF value x1 f1(x1,Q2), not Feynman xF
  double pdf_xf2 = 0.0;  // HepMC PDF value x2 f2(x2,Q2), not Feynman xF

  // Exact non-collinear photon exchange with physical forward-state kinematics
  bool exact_forward_photon_kinematics = false;

  // Hard process scales in GeV
  double muF    = 0.0;  // PDF factorization scale
  double muR    = 0.0;  // QCD renormalization scale
  double scalup = 0.0;  // LHE shower starting or veto scale

  // Running couplings
  double alphaQCD = 0.0;
};

// Per-sample auxiliary data owned by one worker thread
struct MEventWeightState {
  // Logarithm of the inverse phase-space sampling density
  double log_inverse_density = 0.0;

  // Integration-only density collection and amplitude-free event measure
  bool qmetrics = false;
  std::uint64_t sample_index = 0;
  double phase_weight = 0.0;

  // Proposal adaptation controls
  bool adaptation_mode   = false;
  bool include_screening = true;

  // Configure one proposal adaptation evaluation
  void ConfigureAdaptation(bool fast_adaptation) {
    adaptation_mode   = true;
    include_screening = !fast_adaptation;
  }

  // Sequential physics acceptance flags
  bool amplitude_ok  = true;
  bool kinematics_ok = true;
  bool fidcuts_ok    = true;
  bool vetocuts_ok   = true;

  // Technical failure independent of the sequential physics acceptance
  bool technical_failure = false;

  // Event-local amplitude evaluation failed after reaching the amplitude stage
  bool amplitude_failure = false;

  // First diagnostic raised by the event-local amplitude boundary
  std::string amplitude_failure_message;

  // Forced acceptance of the event
  bool forced_accept = false;

  // Reset event-local acceptance and failure bookkeeping
  void ResetStatus() noexcept {
    phase_weight      = 0.0;
    amplitude_ok      = true;
    kinematics_ok     = true;
    fidcuts_ok        = true;
    vetocuts_ok       = true;
    technical_failure = false;
    amplitude_failure = false;
    forced_accept     = false;
    amplitude_failure_message.clear();
  }

  // Record one amplitude failure while preserving its first diagnostic
  void RecordAmplitudeFailure(const std::string &message) {
    amplitude_ok      = false;
    technical_failure = true;
    amplitude_failure = true;
    if (amplitude_failure_message.empty()) { amplitude_failure_message = message; }
  }

  // Compute combined event validity flag
  bool Valid() const {
    return amplitude_ok && kinematics_ok && fidcuts_ok && vetocuts_ok && !amplitude_failure && !technical_failure;
  }
};

// Multipomeron kinematics for one sampled screening chain
struct MultipomeronKinematics {
  M4Vec p1i;
  M4Vec p2i;

  M4Vec q1;
  M4Vec q2;

  M4Vec k;

  M4Vec p1f;
  M4Vec p2f;

  M4Vec p3;
  M4Vec p4;
};

// Generator cut (default parameters set here)
struct GENCUT {
  // Direct daughter laboratory rapidities, class <C>
  double rap_min = -9.0;
  double rap_max = 9.0;

  // Central system laboratory rapidity for factorized and partonic production
  double Y_min = -9.0;
  double Y_max = 9.0;

  // Central system mass, classes <C> and <F> and similar
  double M_min = 0.0;
  double M_max = 0.0;

  // Both <C> and <F> class forward legs
  double forward_pt_min = -1.0;  // Keep at -1 for user setup trigger
  double forward_pt_max = -1.0;  // Keep at -1 for user setup trigger

  // ---------------------------------------

  // Quasi-Elastic phase space <Q> or forward excitation
  double XI_min = 0.0;
  double XI_max = 1.0;

  // Quasielastic minimum sampled |t| in GeV^2, negative uses zero
  double q_t_abs_min = -1.0;

  // Quasielastic maximum sampled |t| in GeV^2, negative uses loop bound
  double q_t_abs_max = -1.0;
};

// Fiducial observable range
struct FIDCUTRANGE {
  bool   active = false;
  double min    = 0.0;
  double max    = 0.0;

  // Compute true when one value is inside this active inclusive range
  bool Contains(double value) const;
};

// PDG-selected fiducial cuts for one particle or composite system
struct FIDPDGCUT {
  std::vector<int> pdg;

  // True entries match by absolute PDG id instead of signed PDG id
  std::vector<bool> pdg_abs;

  FIDCUTRANGE M;
  FIDCUTRANGE Rap;
  FIDCUTRANGE Eta;
  FIDCUTRANGE Pt;
  FIDCUTRANGE Et;

  // Compute true when at least one observable range is configured
  bool HasActiveRange() const { return M.active || Rap.active || Eta.active || Pt.active || Et.active; }

  // Compute this PDG selector in steering-card syntax
  std::string SelectorString() const;
};

// Fiducial cuts (default parameters set here)
struct FIDCUT {
  bool active = false;

  // "Central particles"
  bool   particle_eta_active = false;
  double eta_min             = -30.0;
  double eta_max             = 30.0;

  bool   particle_rap_active = false;
  double rap_min             = -30.0;
  double rap_max             = 30.0;

  bool   particle_pt_active = false;
  double pt_min             = 0.0;
  double pt_max             = 1000000.0;

  bool   particle_Et_active = false;
  double Et_min             = 0.0;
  double Et_max             = 1000000.0;

  // "Central system"
  bool   system_M_active = false;
  double M_min           = 0.0;
  double M_max           = 1000000.0;

  bool   system_Rap_active = false;
  double Y_min             = -30.0;
  double Y_max             = 30.0;

  bool   system_Pt_active = false;
  double Pt_min           = 0.0;
  double Pt_max           = 1000000.0;

  // "Forward system"
  bool   forward_t_active = false;
  double forward_t_min    = 0.0;
  double forward_t_max    = 1000000.0;

  // Independent beam transfers allow a photon Q2 cut without cutting the hadron transfer
  FIDCUTRANGE forward_t1;
  FIDCUTRANGE forward_t2;

  bool   forward_M_active = false;
  double forward_M_min    = 0.0;
  double forward_M_max    = 1000000.0;

  // Fractional longitudinal momentum loss of each forward final leg
  bool   forward_xi_active = false;
  double forward_xi_min    = 0.0;
  double forward_xi_max    = 1.0;

  // Absolute azimuthal separation of the two forward final legs in degrees
  bool   forward_dPhi_active = false;
  double forward_dPhi_min    = 0.0;
  double forward_dPhi_max    = 180.0;

  // Single-PDG cuts apply to every match, while multi-PDG cuts require one
  // exact selector multiset and reject missing or excess matching leaves
  std::vector<FIDPDGCUT> pdg_cuts;

  // Compute true when any generic central particle observable is configured
  bool HasCentralParticleCuts() const {
    return particle_eta_active || particle_rap_active || particle_pt_active || particle_Et_active;
  }

  // Compute true when any central-system observable is configured
  bool HasCentralSystemCuts() const { return system_M_active || system_Rap_active || system_Pt_active; }

  // Compute true when any forward-system observable is configured
  bool HasForwardCuts() const {
    return forward_t_active || forward_t1.active || forward_t2.active || forward_M_active || forward_xi_active || forward_dPhi_active;
  }

  // Compute true when the forward event kinematics pass all configured cuts
  bool PassForward(const LORENTZSCALAR &lts) const;

  // Compute true when the two forward legs pass the azimuthal-separation cut
  bool PassForwardDeltaPhi(const LORENTZSCALAR &lts) const;

  // Compute true when the available forward legs pass the momentum-loss cut
  bool PassForwardXi(const LORENTZSCALAR &lts) const;

  // Compute true when the central system passes all configured cuts
  bool PassCentralSystem(const LORENTZSCALAR &lts) const;

  // Compute true when every stable central particle passes the common cuts
  bool PassCentralParticles(const std::vector<MDecayBranch> &tree) const;

  // Compute true when all PDG-selected particle and system cuts pass
  bool PassSelectedParticles(const std::vector<MDecayBranch> &tree) const;
};

// Particle charge category selected by one veto domain
enum class VetoCharge { Any, Charged, Neutral };

// Kinematic and particle-source selection for one veto domain
struct VETODOMAIN {
  double eta_min = -30.0;
  double eta_max = 30.0;

  double pt_min = 0.0;
  double pt_max = 1000000.0;

  bool       source_forward = true;
  bool       source_central = true;
  VetoCharge charge         = VetoCharge::Any;
};

// Complete event-veto configuration
struct VETOCUT {
  bool                    active = false;
  std::vector<VETODOMAIN> cuts;

  // Compute true when no stable forward or central particle enters a veto domain
  bool Pass(const LORENTZSCALAR &lts) const;

  // Test veto domains with an explicit central momentum tree
  bool Pass(const LORENTZSCALAR &lts, const std::vector<MDecayBranch> &central) const;
};

// Complete configuration and mutable data for one process worker
struct MProcessState {
  // Kinematic and amplitude data
  LORENTZSCALAR lts;

  // Event-local multipomeron screening chain
  std::vector<MultipomeronKinematics> multipomeron_chain;
  double                              multipomeron_impact_parameter = 0.0;

  // Generation, fiducial and veto cuts
  GENCUT  gcuts;
  FIDCUT  fcuts;
  VETOCUT vetocuts;

  // Process-local random numbers
  MRandom random;

  // Process-local QED ISR and FSR service
  radiative::MYFS radiative;

  // Immutable SOFT exchange and Good Walker snapshot
  SoftModelPtr soft_model;

  // Immutable complete model tune shared by all workers
  MModelTunePtr model_tune;

  // Validated heavy-ion UPC controls
  nuclear::UPCParam upc_param;
  bool              upc_configured = false;

  // Event-local orthogonal nuclear final-sector probabilities and selection
  nuclear::FinalWeights upc_final_weight{};
  std::array<double, 3> upc_photo{};
  nuclear::FinalState   upc_final;

  // Own the nuclear backend per worker and the completed HepMC3 event per point
  std::optional<nuclear::MFinal> nuclear_final;
  std::optional<HepMC3::GenEvent> nuclear_event;

  // Immutable forward excitation and fragmentation parameters
  MNstarParamPtr nstar_param;

  // Cross-section statistical factor for identical final states
  double symmetry_factor = 0.0;

  // Steering parameters
  std::string  process;
  std::string  phase_space_class;
  std::string  decay_mode;
  bool         screening      = false;
  int          excitation     = 0;
  BeamFragType beamfrag       = BeamFragType::None;
  std::int64_t usercuts       = 0;
  int          flat_amplitude = 0;

  // Phase-space control
  bool   flat_mass2           = false;
  bool   flat_mass2_user      = false;
  double offshell_widths      = 5.0;
  bool   offshell_widths_user = false;
  double width_min            = 1e-12;
  bool   spingen_user         = false;
  bool   spindec_user         = false;
  bool   qmetrics_user     = false;
  bool   mmax_user            = false;
};

}  // namespace gra

#endif
