// Nuclear UPC steering and scalar screening convolution
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARUPC_H
#define MNUCLEARUPC_H

#include <array>
#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Math/MPolarQuadrature.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MBreakup.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MGlauber.h"
#include "Graniitti/Nuclear/MGGCF.h"
#include "Graniitti/Nuclear/MMass.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MPhoto.h"
#include "Graniitti/Nuclear/MPhoton.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra::nuclear {

// Select a lepton, proton or nuclear UPC beam
enum class BeamType { Lepton, Proton, Nucleus };

// Select an observed neutron class
struct NeutronSelection {
  enum class Type { Any, Equal, NotEqual, Greater, GreaterEqual, Less, LessEqual };
  Type type = Type::Any;
  std::size_t count = 0;

  // Check the emitted neutron multiplicity against the parsed selection
  bool Accept(std::size_t n) const;
};

// Select the spatial resolution of nuclear currents
enum class StructureType { Smooth, Nucleon, Hotspot };

// Select the hadronic survival calculation
enum class SurvivalType { Optical, OpticalGGCF, MCGGCF };

// Select an optional exclusive nuclear reaction calculation
enum class ReactionType { Internal, External };

// Select the UPC services required by one process runtime
struct UPCMode {
  std::array<bool, 2> emission  = {true, true};
  std::array<bool, 2> target    = {true, true};
  std::array<bool, 2> current   = {true, true};
  bool                screening = true;
};

// Store the UPC impact profile convolution controls
struct ConvolutionParam {
  double       b_max           = 0.0;  // Impact parameter transform limit [fm]
  unsigned int smooth_b_nodes  = 0;    // Smooth impact parameter nodes
  unsigned int sample_b_nodes  = 0;    // Sampled impact parameter nodes
  unsigned int b_phi_nodes     = 0;    // Sampled impact azimuthal nodes
  unsigned int smooth_kt_nodes = 0;    // Smooth radial lookup nodes
  unsigned int sample_kt_nodes = 0;    // Sampled-profile radial transform nodes
};

// Parse one lowercase beam name
BeamType ParseBeamType(const std::string &name);

// Parse one lowercase nuclear-sector name
CoherenceType ParseCoherenceType(const std::string &name);

// Parse one lowercase neutron-class name
NeutronSelection ParseNeutronSelection(const std::string &name);

// Parse one lowercase nuclear-structure name
StructureType ParseStructureType(const std::string &name);

// Parse one lowercase survival-model name
SurvivalType ParseSurvivalType(const std::string &name);

// Parse one lowercase photonuclear model name
PhotoModel ParsePhotoModel(const std::string &name);

// Compute one lowercase beam name
std::string BeamName(BeamType type);

// Compute one lowercase nuclear-sector name
std::string CoherenceName(CoherenceType type);

// Compute one lowercase neutron-class name
std::string NeutronName(NeutronSelection type);

// Compute one lowercase nuclear-structure name
std::string StructureName(StructureType type);

// Compute one lowercase survival-model name
std::string SurvivalName(SurvivalType type);

// Compute one lowercase photonuclear model name
std::string PhotoModelName(PhotoModel model);

// Store complete two-leg UPC physics and numerical controls
struct UPCParam {
  std::optional<ReactionType>  reaction;                // Empty keeps inclusive forward systems
  bool                         additional_emd = false;  // Additional photon exchange independent of hard breakup
  std::string                  final_library;           // Optional HepMC3 nuclear reaction library
  std::string                  final_config = "{}";     // External nuclear model configuration
  MassParam                    mass;                    // Evaluated masses and explicit theoretical completion
  unsigned int                 decay_nodes   = 0;       // Continuum decay integration intervals
  unsigned int                 decay_steps   = 0;       // Maximum sequential emissions
  unsigned int                 impulse_nodes = 0;       // Fermi sphere integration nodes
  double                       emd_focus     = 0.0;     // Survival-focused absorption-history proposal fraction
  bool                         emd_condition = false;   // Condition coherent Xn legs on photon absorption
  std::array<CoherenceType, 2> emission      = {CoherenceType::Coherent, CoherenceType::Coherent};  // EPA sectors
  std::array<CoherenceType, 2> target        = {CoherenceType::Coherent, CoherenceType::Coherent};  // Target sectors
  PhotoModel                   photo_model   = PhotoModel::Glauber;                   // Photonuclear target model
  std::array<NeutronSelection, 2>   neutron       = {};  // Neutron classes
  StructureType                structure     = StructureType::Smooth;                 // Nuclear current structure
  SurvivalType                 survival      = SurvivalType::OpticalGGCF;             // Hadronic survival model
  std::optional<double>        sigma_nn;           // Inelastic NN cross section override [mb]
  GGCFParam                    ggcf;               // Named eikonal for GGCF survival
  std::string                  survival_eikonal;   // Resolved NN profile model
  double                       s_nn = 0.0;         // Nucleon-nucleon energy squared [GeV^2]
  GlauberParam                 glauber;            // Glauber-Gribov controls
  GeometryParam                geometry;           // Nuclear density controls
  std::array<BreakupParam, 2>  breakup;            // Derived leg breakup models
  EMDParam                     emd;                // Responses selected by target isotope
  std::array<PhotoParam, 2>    photo;              // Photonuclear target controls
  ConfigParam                  config;             // Nuclear configuration controls
  std::size_t                  current_count = 0;  // Current-only configurations per nuclear leg
  math::PolarParam             loop;               // Nuclear screening loop quadrature
  ConvolutionParam             convolution;        // Impact profile convolution controls
  form::ParamStore             structure_param;    // Proton structure parameters
  std::string                  table_fingerprint;  // Steering and external-table identity
};

// Store one precomputed two-dimensional screening-loop node
struct LoopNode {
  double                            kt     = 0.0;
  double                            phi    = 0.0;
  double                            kx     = 0.0;
  double                            ky     = 0.0;
  std::complex<double>              weight = {0.0, 0.0};
  std::vector<std::complex<double>> sample_weight;
};

// Store deterministic quadrature and transform tables for config screening
struct ConfigGrid {
  std::vector<double>               b_node;
  std::vector<double>               b_weight;
  std::vector<std::complex<double>> b_unit;
  std::vector<std::complex<double>> b_phase;
  std::vector<double>               k_node;
  std::vector<double>               bessel;
};

// Store one importance-normalized orthogonal excitation amplitude
struct ExcitationChannel {
  double born = 1.0;  // Large-b amplitude, zero for every non-vacuum absorption history
  double radius = 0.0;  // Largest physical excitation support in fm
  std::vector<double> momentum, bare, screened;
  std::vector<double> impact, amplitude;  // Unscreened sqrt(P(h|b)/q(h)) for this orthogonal history
};

// Own one immutable scalar UPC survival convolution
class MUPC {
 public:
  using SectorRatios = std::array<std::vector<std::complex<double>>, 2>;

  // Construct and precompute one heterogeneous UPC convolution
  MUPC(std::array<BeamType, 2> type, std::array<std::shared_ptr<const MNucleus>, 2> nucleus, UPCParam param,
       std::array<std::shared_ptr<const MConfigBank>, 2> bank = {}, UPCMode mode = {});

  // Construct and precompute one nuclear screening convolution
  MUPC(const MNucleus &beam1, const MNucleus &beam2, UPCParam param, std::shared_ptr<const MConfigBank> bank1 = nullptr,
       std::shared_ptr<const MConfigBank> bank2 = nullptr, UPCMode mode = {});

  // Compute the validated UPC controls
  const UPCParam &Param() const { return param_; }

  // Compute one checked beam profile type
  BeamType Type(int leg) const;

  // Compute one checked beam nucleus
  const MNucleus *Nucleus(int leg) const;

  // Compute one optional immutable configuration bank
  const MConfigBank *Bank(int leg) const;

  // Compute one optional event-local hotspot bank
  const MHotSpot *HotSpot(int leg) const;

  // Compute the finite-range Glauber model
  const MGlauber *Glauber() const { return glauber_.get(); }

  // Compute whether both beams are hadrons
  bool HasHadronicPair() const;

  // Compute whether nuclear hadronic screening is enabled
  bool Screening() const { return screening_; }

  // Compute whether hadronic survival requires an impact-space convolution
  bool HadronicConvolution() const;

  // Compute whether survival or a resolved excitation requires convolution
  bool Convolution() const { return HadronicConvolution() || channel_ != nullptr; }

  // Compute the asymptotic excitation amplitude multiplying the Born term
  double BornWeight() const { return channel_ ? channel_->born : 1.0; }

  // Attach a worker-local excitation channel without changing shared nuclear tables
  std::shared_ptr<const MUPC> WithExcitation(std::shared_ptr<const ExcitationChannel> channel, bool survival) const;

  // Compute one optional preconstructed nuclear photon source
  const MPhoton *Photon(int leg) const;

  // Compute one optional preconstructed photonuclear target model
  const MPhoto *Photo(int leg) const;

  // Compute one electromagnetic breakup model
  const MBreakup *Breakup(int leg) const;

  // Compute event-local screening nodes clustered at shifted EPA currents
  std::vector<LoopNode> Nodes(const std::array<M3Vec, 2> &transfer) const;

  // Compute whether a configuration-current ensemble is available
  bool HasSamples() const;

  // Include mean and fluctuation sources before projecting a sampled nuclear final state
  CoherenceType SourceSector(int leg, CoherenceType selected) const;

  // Compute the complete Cartesian configuration sample count
  std::size_t SampleCount() const { return sample_count_; }

  // Compute the independent upper and lower bank sizes
  const std::array<std::size_t, 2> &SampleShape() const { return sample_shape_; }

  // Compute one leg's configuration index in the paired ensemble
  std::size_t SampleIndex(std::size_t sample, int leg) const;

  // Compute normalized charge currents for one explicit emission sector
  std::vector<std::complex<double>> EmissionRatios(int leg, CoherenceType emission, const M3Vec &q) const;

  // Compute coherent and incoherent charge-current ratios in that order
  SectorRatios EmissionComponents(int leg, const M3Vec &q) const;

  // Compute normalized matter currents for one explicit target sector
  std::vector<std::complex<double>> TargetRatios(int leg, CoherenceType target, const M3Vec &q) const;

  // Compute target ratios from one event-local configuration-current bank
  std::vector<std::complex<double>> TargetCurrentRatios(int leg, CoherenceType target,
                                                        const PhotoCurrent &current) const;

  // Compute the momentum-space hadronic profile in GeV^-2
  double ProfileTransform(double kt) const;

  // Compute the selected no-additional-interaction amplitude
  double SurvivalAmp(double b) const;

  // Sample one complete event-local nuclear fluctuation ensemble
  std::shared_ptr<const MUPC> Sample(MRandom &random, bool prepare_convolution = true, std::size_t count = 0) const;

 private:
  std::shared_ptr<const ExcitationChannel>             channel_;
  UPCParam                                             param_;
  std::array<BeamType, 2>                              type_;
  std::array<std::shared_ptr<const MNucleus>, 2>       nucleus_;
  std::shared_ptr<const MGlauber>                      glauber_;
  std::array<std::shared_ptr<const MBreakup>, 2>       breakup_;
  std::array<std::shared_ptr<const MConfigBank>, 2>    bank_;
  std::array<std::shared_ptr<const MConfigSampler>, 2> sampler_;
  std::array<std::optional<MHotSpot>, 2>               hotspot_;
  bool                                                 screening_ = true;
  std::shared_ptr<const ConfigGrid>                    config_grid_;
  std::array<std::shared_ptr<const MPhoton>, 2>        photon_;
  std::array<std::shared_ptr<const MPhoto>, 2>         photo_;
  std::vector<double>                                  b_node_;
  std::vector<double>                                  b_weight_;
  std::vector<double>                                  profile_;
  std::vector<std::vector<std::complex<double>>>       sample_profile_;
  std::array<std::size_t, 2>                           sample_shape_ = {1, 1};
  std::size_t                                          sample_count_ = 0;
  std::vector<double>                                  loop_k_;
  std::vector<double>                                  loop_profile_;
  std::vector<std::complex<double>>                    loop_harmonic_;

  // Construct one event-local runtime sharing all deterministic model tables
  MUPC(const MUPC &model, std::array<std::shared_ptr<const MConfigBank>, 2> bank,
       std::array<std::optional<MHotSpot>, 2> hotspot, bool prepare_convolution);

  // Validate all UPC controls and optional configuration banks
  void Validate() const;

  // Validate only the new event configuration banks against the immutable beam model
  void ValidateBanks() const;

  // Prepare the impact-parameter and momentum-space convolution
  void Prepare();

  // Construct the full physics and numerical disk-cache key
  std::string CacheKey() const;

  // Load and validate one compressed JSON screening table
  bool ReadCache(const std::string &filename, const std::string &key);

  // Save one complete compressed JSON screening table
  void WriteCache(const std::string &filename, const std::string &key) const;

  // Load one deterministic config-screening grid
  bool ReadConfigCache(const std::string &filename, const std::string &key);

  // Save one deterministic config-screening grid
  void WriteConfigCache(const std::string &filename, const std::string &key) const;

  // Prepare deterministic config-screening quadrature and transform tables
  void PrepareConfigGrid();

  // Prepare the dense radial transform used by smooth focused quadrature
  void PrepareScalarTransform();

  // Interpolate the smooth radial profile transform
  double ScalarTransform(double kt) const;

  // Prepare the complete Cartesian configuration ensemble
  void PrepareSamples();

  // Expand one leg-local current bank over the paired ensemble
  std::vector<std::complex<double>> ExpandRatios(const std::vector<std::complex<double>> &ratio, int leg) const;

  // Compute one charge or matter current bank
  std::vector<std::complex<double>> Currents(int leg, const M3Vec &q, bool charge) const;

  // Normalize and expand one leg-local current bank
  std::vector<std::complex<double>> Ratios(int leg, CoherenceType type,
                                           const std::vector<std::complex<double>> &current) const;
};

}  // namespace gra::nuclear

#endif
