// Nuclear UPC amplitude screening and Good-Walker projection
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MUPCSCREEN_H
#define MUPCSCREEN_H

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Nuclear/MUPC.h"

namespace gra::nuclear {

// Select one photon direction in a photonuclear amplitude
enum class PhotoDirection : std::uint8_t { Upper, Lower };

// Store one sampled orthogonal two-leg nuclear final state
struct FinalState {
  bool                         valid = false;
  std::array<CoherenceType, 2> leg   = {CoherenceType::Coherent, CoherenceType::Coherent};
};

using FinalWeights = std::array<double, 4>;

// Compute the compact upper-lower final-sector index
std::size_t FinalPairIndex(CoherenceType upper, CoherenceType lower);

// Store one explicit photon-emitter and photonuclear-target channel
struct PhotoChannel {
  PhotoDirection direction = PhotoDirection::Upper;
  CoherenceType  emission  = CoherenceType::Coherent;
  CoherenceType  target    = CoherenceType::Coherent;

  // Compute the upper-lower final pair from the photon direction
  std::size_t Pair() const {
    return direction == PhotoDirection::Upper ? FinalPairIndex(emission, target) : FinalPairIndex(target, emission);
  }

  // Compare two photonuclear channels field by field
  constexpr bool operator==(const PhotoChannel& other) const noexcept = default;
};

// Compute the explicit coherent and incoherent sectors for one choice
std::vector<CoherenceType> CoherenceSectors(CoherenceType sector);

// Compute the canonical upper-then-lower photonuclear channel order
std::vector<PhotoChannel> PhotoChannels(const UPCParam& param);

// Compute only physical photonuclear channels of one active UPC runtime
std::vector<PhotoChannel> PhotoChannels(const MUPC& upc);

// Coherently combine photon directions ending in the same final sector
std::vector<std::complex<double>> CombinePhotoChannels(const std::vector<std::complex<double>>& amplitude,
                                                       const std::vector<PhotoChannel>&         channel);

using NuclearGoodWalkerWeights = std::array<std::array<double, 2>, 2>;

// Project one Cartesian two-leg eigen-amplitude into all resolved final sectors
NuclearGoodWalkerWeights NuclearGoodWalkerProject(const std::vector<std::complex<double>>& amplitude,
                                                  const std::array<std::size_t, 2>&        shape);

// Select the nuclear amplitude convolution used for one hard process
enum class ScreenType { Scalar, Fusion, Photo };

// Store the resolved two-photon source layout
struct FusionLayout {
  std::size_t                               rows = 0;
  std::array<std::vector<CoherenceType>, 2> sector;
};

// Store the complete process-independent nuclear amplitude layout
struct ScreenLayout {
  ScreenType                type = ScreenType::Scalar;
  FusionLayout              fusion;
  std::vector<PhotoChannel> photo;
  FinalState                fixed;
};

// Aggregate hard-channel weights into the four orthogonal final sectors
FinalWeights FinalSectorWeights(const ScreenLayout& layout, const std::vector<double>& helicity_norm);

// Select one orthogonal final sector from finite nonnegative weights
FinalState SampleFinalState(const FinalWeights& weight, double unit);

// Store one photoproduction direction term with its target currents
struct PhotoTerm {
  PhotoDirection                             direction = PhotoDirection::Upper;
  std::vector<std::complex<double>>          amplitude;
  std::array<std::optional<PhotoCurrent>, 2> photo_current;
};

// Store one hard amplitude and its two nuclear-rest-frame momentum transfers
struct ScreenPoint {
  std::vector<std::complex<double>>              amplitude;
  std::vector<std::vector<std::complex<double>>> flow;
  std::array<M3Vec, 2>                           transfer{};
  std::array<std::optional<PhotoCurrent>, 2>     photo_current;
  std::vector<PhotoTerm>                         photo_terms;
};

// Store screened amplitudes and their raw external-helicity and color-projection norms
struct ScreenResult {
  std::array<double, 3>             photo{};  // |A_upper|^2, |A_lower|^2, 2 Re(A_upper* A_lower) in the coherent sector
  std::vector<std::complex<double>> amplitude;
  std::vector<double>               helicity_norm;
  std::vector<double>               color_norm;
  std::vector<FinalWeights>         color_sector;
};

// Accumulate one event-local nuclear UPC screening convolution
class MUPCScreen {
 public:
  // Initialize one convolution or a bare projection with unit Born weight
  MUPCScreen(const MUPC& upc, ScreenLayout layout, ScreenPoint born);

  // Add one shifted hard amplitude at a precomputed Glauber node
  void Add(const LoopNode& node, const ScreenPoint& point);

  // Project the accumulated ensemble into resolved physical sectors
  ScreenResult Result() const;

 private:
  const MUPC&                                                 upc_;
  ScreenLayout                                                layout_;
  std::size_t                                                 flow_count_    = 0;
  std::size_t                                                 born_size_     = 0;
  std::size_t                                                 sample_count_  = 0;
  std::size_t                                                 channel_count_ = 0;
  std::size_t                                                 hard_count_    = 0;
  bool                                                        photo_terms_   = false;
  std::array<bool, 2>                                         photo_current_{};
  std::array<std::array<bool, 4>, 2>                         selected_{};
  std::vector<std::vector<std::complex<double>>>              born_bank_;
  std::complex<double>                                        scalar_weight_ = 0.0;
  std::vector<std::vector<std::complex<double>>>              scalar_bank_;
  std::vector<std::vector<std::vector<std::complex<double>>>> sample_bank_;
  std::array<std::vector<std::complex<double>>, 2>            photo_bank_;

  // Validate a complete point before changing the accumulated amplitude
  void ValidatePoint(const ScreenPoint& point) const;

  // Initialize the scalar screening convolution
  void PrepareScalar(const ScreenPoint& born);

  // Initialize the resolved two-photon configuration convolution
  void PrepareFusion(const ScreenPoint& born);

  // Initialize the resolved photonuclear configuration convolution
  void PreparePhoto(const ScreenPoint& born);

  // Add one scalar Glauber contribution
  void AddScalar(const LoopNode& node, const ScreenPoint& point);

  // Accumulate Born and shifted two-photon currents with scalar or sample weights
  void AccumulateFusionPoint(const ScreenPoint& point, std::complex<double> common,
                             const std::vector<std::complex<double>>& weight = {});

  // Compute all sample ratios for one resolved photonuclear current
  std::vector<std::vector<std::complex<double>>> PhotoRatios(const std::array<M3Vec, 2>&                       transfer,
                                                             const std::array<std::optional<PhotoCurrent>, 2>& current,
                                                             std::optional<PhotoDirection> direction) const;

  // Accumulate one photonuclear point with scalar or sample-dependent weights
  void AccumulatePhotoPoint(const ScreenPoint& point, std::complex<double> common,
                            const std::vector<std::complex<double>>& weight = {});

  // Accumulate one photonuclear amplitude into final-sector sample blocks
  void AccumulatePhoto(const std::vector<std::complex<double>>&              amplitude,
                       const std::vector<std::vector<std::complex<double>>>& ratio, std::complex<double> common,
                       const std::vector<std::complex<double>>& weight, std::size_t bank);

  // Compute the selected scalar or resolved Good-Walker result
  ScreenResult ScalarResult() const;
  ScreenResult EnsembleResult() const;
};

}  // namespace gra::nuclear

#endif
