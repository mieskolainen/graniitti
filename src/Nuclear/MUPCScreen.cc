// Nuclear UPC amplitude screening and Good-Walker projection
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MUPCScreen.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <numeric>
#include <string>
#include <utility>

#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::nuclear {

namespace {

// Compute the source-sector pair of one fusion amplitude row
std::array<std::size_t, 2> FusionPair(const FusionLayout& layout, const std::size_t component) {
  const std::size_t lower_rows = layout.rows * layout.sector[1].size();
  const std::size_t block      = layout.rows * layout.sector[0].size() * lower_rows;
  const std::size_t slot       = component % block;
  return {(slot / lower_rows) / layout.rows, (slot % lower_rows) / layout.rows};
}

// Compute the hard and photon-polarization row shared by all nuclear source sectors
std::size_t FusionRow(const FusionLayout& layout, const std::size_t component) {
  const std::size_t lower_rows = layout.rows * layout.sector[1].size();
  const std::size_t block      = layout.rows * layout.sector[0].size() * lower_rows;
  return (component / block) * layout.rows * layout.rows +
         ((component % block) / lower_rows) % layout.rows * layout.rows + component % layout.rows;
}

// Compute only the explicit photon-emission ratios used by one leg
std::array<std::vector<std::complex<double>>, 2> EmissionRatios(const MUPC& upc, const int leg,
                                                                const std::vector<CoherenceType>& sector,
                                                                const M3Vec&                      q) {
  if (sector.empty() || sector.size() > 2) { throw AmplitudeFailure("EmissionRatios: invalid sector count"); }
  std::array<std::vector<std::complex<double>>, 2> ratio;
  if (sector.size() == 1) {
    ratio[0] = upc.EmissionRatios(leg, sector[0], q);
    return ratio;
  }

  const auto component = upc.EmissionComponents(leg, q);
  for (const auto& i : indices(sector)) {
    if (sector[i] == CoherenceType::Coherent) {
      ratio[i] = component[0];
    } else if (sector[i] == CoherenceType::Incoherent) {
      ratio[i] = component[1];
    } else {
      throw AmplitudeFailure("EmissionRatios: unresolved sector");
    }
  }
  return ratio;
}

// Compute one representative channel for each physical final-sector pair
std::array<std::size_t, 4> PhotoRepresentatives(const std::vector<PhotoChannel>& channel) {
  std::array<std::size_t, 4> representative;
  representative.fill(std::numeric_limits<std::size_t>::max());
  for (const auto& i : indices(channel)) {
    const std::size_t pair = channel[i].Pair();
    if (representative[pair] == std::numeric_limits<std::size_t>::max()) { representative[pair] = i; }
  }
  return representative;
}

// Enumerate photon directions before their mean and fluctuation sources
std::vector<PhotoChannel> PhotoSources(const std::array<CoherenceType, 2>& emission,
                                       const std::array<CoherenceType, 2>& target) {
  std::vector<PhotoChannel> channel;
  channel.reserve(8);
  for (const auto direction : {PhotoDirection::Upper, PhotoDirection::Lower}) {
    const std::size_t emitter = direction == PhotoDirection::Upper ? 0U : 1U;
    for (const auto source : CoherenceSectors(emission[emitter])) {
      for (const auto transition : CoherenceSectors(target[1U - emitter])) {
        channel.push_back({direction, source, transition});
      }
    }
  }
  return channel;
}

// Compute the helicity-summed coherent photon direction density
std::array<double, 3> PhotoInterference(const std::array<std::vector<std::complex<double>>, 2>& amplitude) {
  return {gra::SquaredNorm(amplitude[0]), gra::SquaredNorm(amplitude[1]),
          2.0 * std::real(gra::InnerProduct(amplitude[0], amplitude[1]))};
}

// Validate all supplied target currents against their nuclear configuration banks
void ValidateCurrents(const std::array<std::optional<PhotoCurrent>, 2>& current, const MUPC& upc) {
  for (const auto& leg : indices(current)) {
    if (!current[leg].has_value()) { continue; }
    const auto* bank  = upc.Bank(static_cast<int>(leg + 1));
    const auto& value = *current[leg];
    if (bank == nullptr || value.sample.size() != bank->Size() || !gra::AllFinite(value.sample) ||
        !std::isfinite(value.stat.mean.real()) || !std::isfinite(value.stat.mean.imag()) ||
        !std::isfinite(value.stat.second) || !std::isfinite(value.stat.variance) || value.stat.second < 0.0 ||
        value.stat.variance < 0.0) {
      throw AmplitudeFailure("MUPCScreen: invalid nuclear target current");
    }
  }
}

// Validate that separately shadowed terms reproduce one coherent amplitude
void ValidatePhotoTerms(const ScreenPoint& point, const std::size_t size, const std::vector<PhotoChannel>& channel,
                        const bool required, const MUPC& upc) {
  if (point.photo_terms.empty()) {
    if (required) { throw AmplitudeFailure("MUPCScreen: photonuclear term bank is absent"); }
    return;
  }
  if (!required) { throw AmplitudeFailure("MUPCScreen: photonuclear term layout changed"); }
  if (channel.empty() || size % channel.size() != 0) {
    throw AmplitudeFailure("MUPCScreen: photonuclear term channel layout changed");
  }

  std::vector<std::complex<double>> sum(size, 0.0);
  std::vector<double>               scale(size, 0.0);
  for (const auto& term : point.photo_terms) {
    if (term.amplitude.size() != size || !gra::AllFinite(term.amplitude)) {
      throw AmplitudeFailure("MUPCScreen: photonuclear term amplitude layout changed");
    }
    if (term.direction != PhotoDirection::Upper && term.direction != PhotoDirection::Lower) {
      throw AmplitudeFailure("MUPCScreen: invalid photonuclear direction");
    }
    ValidateCurrents(term.photo_current, upc);
    for (const auto& h : indices(term.amplitude)) {
      if (channel[h % channel.size()].direction != term.direction && std::abs(term.amplitude[h]) > 0.0) {
        throw AmplitudeFailure("MUPCScreen: photonuclear term occupies the wrong direction");
      }
      sum[h] += term.amplitude[h];
      scale[h] += std::abs(term.amplitude[h]);
    }
    const std::size_t emitter = term.direction == PhotoDirection::Upper ? 0U : 1U;
    const std::size_t target  = 1U - emitter;
    if (term.photo_current[emitter].has_value()) {
      throw AmplitudeFailure("MUPCScreen: photonuclear term carries an emitter-side target current");
    }
    if (upc.Bank(static_cast<int>(target + 1)) != nullptr && !term.photo_current[target].has_value()) {
      throw AmplitudeFailure("MUPCScreen: photonuclear term is missing its nuclear target current");
    }
  }
  if (!gra::AllFinite(sum) || !gra::AllFinite(scale)) {
    throw AmplitudeFailure("MUPCScreen: non-finite photonuclear term sum");
  }
  for (const auto& h : indices(sum)) {
    const double tolerance =
        128.0 * std::numeric_limits<double>::epsilon() * std::max(scale[h], std::abs(point.amplitude[h]));
    if (std::abs(sum[h] - point.amplitude[h]) > tolerance) {
      throw AmplitudeFailure("MUPCScreen: photonuclear terms do not sum to the amplitude");
    }
  }
}

}  // namespace

// Compute the compact upper-lower final-sector index
std::size_t FinalPairIndex(const CoherenceType upper, const CoherenceType lower) {
  if ((upper != CoherenceType::Coherent && upper != CoherenceType::Incoherent) ||
      (lower != CoherenceType::Coherent && lower != CoherenceType::Incoherent)) {
    throw AmplitudeFailure("FinalPairIndex: final sector must be resolved");
  }
  return 2 * static_cast<std::size_t>(upper) + static_cast<std::size_t>(lower);
}

// Project one Cartesian two-leg eigen-amplitude into all resolved final sectors
// W_CC = |<A_ij>|^2, W_IC = <|<A_ij>_j - <A_ij>|^2>_i
// W_CI = <|<A_ij>_i - <A_ij>|^2>_j
// W_II = <|A_ij - <A_ij>_j - <A_ij>_i + <A_ij>|^2>_ij
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
NuclearGoodWalkerWeights NuclearGoodWalkerProject(const std::vector<std::complex<double>>& amplitude,
                                                  const std::array<std::size_t, 2>&        shape) {
  if (shape[0] == 0 || shape[1] == 0 || shape[0] > std::numeric_limits<std::size_t>::max() / shape[1] ||
      amplitude.size() != shape[0] * shape[1] || !gra::AllFinite(amplitude)) {
    throw AmplitudeFailure("NuclearGoodWalkerProject: invalid Cartesian ensemble");
  }
  std::vector<std::complex<double>> upper(shape[0], 0.0);
  std::vector<std::complex<double>> lower(shape[1], 0.0);
  std::complex<double>              mean = 0.0;
  for (std::size_t i = 0; i < shape[0]; ++i) {
    for (std::size_t j = 0; j < shape[1]; ++j) {
      const std::complex<double> value = amplitude[i * shape[1] + j];
      upper[i] += value;
      lower[j] += value;
      mean += value;
    }
  }
  const double inverse_upper = 1.0 / static_cast<double>(shape[0]);
  const double inverse_lower = 1.0 / static_cast<double>(shape[1]);
  gra::Scale(upper, inverse_lower);
  gra::Scale(lower, inverse_upper);
  mean *= inverse_upper * inverse_lower;

  NuclearGoodWalkerWeights result{};
  result[0][0] = std::norm(mean);
  for (const auto& value : upper) { result[1][0] += std::norm(value - mean) * inverse_upper; }
  for (const auto& value : lower) { result[0][1] += std::norm(value - mean) * inverse_lower; }
  for (std::size_t i = 0; i < shape[0]; ++i) {
    for (std::size_t j = 0; j < shape[1]; ++j) {
      const std::complex<double> residual = amplitude[i * shape[1] + j] - upper[i] - lower[j] + mean;
      result[1][1] += std::norm(residual) * inverse_upper * inverse_lower;
    }
  }
  if (!gra::AllFinite(result[0]) || !gra::AllFinite(result[1])) {
    throw AmplitudeFailure("NuclearGoodWalkerProject: non-finite projection");
  }
  return result;
}

// Compute the explicit coherent and incoherent sectors for one choice
std::vector<CoherenceType> CoherenceSectors(const CoherenceType sector) {
  if (sector == CoherenceType::Inclusive) { return {CoherenceType::Coherent, CoherenceType::Incoherent}; }
  return {sector};
}

// Compute the canonical upper-then-lower photonuclear channel order
std::vector<PhotoChannel> PhotoChannels(const UPCParam& param) { return PhotoSources(param.emission, param.target); }

// Compute only physical photonuclear channels of one active UPC runtime
std::vector<PhotoChannel> PhotoChannels(const MUPC& upc) {
  const auto can_emit = [&upc](const int leg) {
    return upc.Type(leg) != BeamType::Nucleus || upc.Photon(leg) != nullptr;
  };
  const auto can_target = [&upc](const int leg) {
    return upc.Type(leg) == BeamType::Proton || (upc.Type(leg) == BeamType::Nucleus && upc.Photo(leg) != nullptr);
  };
  const bool upper    = can_emit(1) && can_target(2);
  const bool lower    = can_emit(2) && can_target(1);
  auto       emission = upc.Param().emission;
  auto       target   = upc.Param().target;
  for (const auto& leg : indices(emission)) {
    emission[leg] = upc.SourceSector(static_cast<int>(leg + 1), emission[leg]);
    target[leg]   = upc.SourceSector(static_cast<int>(leg + 1), target[leg]);
  }
  auto channel = PhotoSources(emission, target);
  channel.erase(std::remove_if(channel.begin(), channel.end(),
                               [&](const auto& item) {
                                 return (item.direction == PhotoDirection::Upper && !upper) ||
                                        (item.direction == PhotoDirection::Lower && !lower);
                               }),
                channel.end());
  return channel;
}

// Coherently combine photon directions ending in the same final sector
// A_f = sum_{channels c -> f} A_c
std::vector<std::complex<double>> CombinePhotoChannels(const std::vector<std::complex<double>>& amplitude,
                                                       const std::vector<PhotoChannel>&         channel) {
  if (amplitude.empty() || channel.empty() || amplitude.size() % channel.size() != 0 || !gra::AllFinite(amplitude)) {
    throw AmplitudeFailure("CombinePhotoChannels: invalid amplitude layout");
  }
  const auto                        representative = PhotoRepresentatives(channel);
  std::vector<std::complex<double>> combined(amplitude.size(), 0.0);
  std::vector<bool>                 occupied(amplitude.size(), false);
  for (const auto& h : indices(amplitude)) {
    const std::size_t hard   = h / channel.size();
    const std::size_t source = h % channel.size();
    const std::size_t pair   = channel[source].Pair();
    const std::size_t output = hard * channel.size() + representative[pair];
    if (pair != 0 && std::abs(amplitude[h]) > 0.0) {
      if (occupied[output]) {
        throw AmplitudeFailure("CombinePhotoChannels: incoherent interference requires configuration currents");
      }
      occupied[output] = true;
    }
    combined[output] += amplitude[h];
  }
  if (!gra::AllFinite(combined)) { throw AmplitudeFailure("CombinePhotoChannels: non-finite amplitude sum"); }
  return combined;
}

// Aggregate amplitude components into the four orthogonal final sectors
// W_f = sum_{components h in f} W_h
FinalWeights FinalSectorWeights(const ScreenLayout& layout, const std::vector<double>& helicity_norm) {
  FinalWeights weight{};
  double       total_weight = 0.0;
  for (const double value : helicity_norm) {
    if (!std::isfinite(value) || value < 0.0) {
      throw AmplitudeFailure("FinalSectorWeights: channel weights must be finite and nonnegative");
    }
    total_weight += value;
  }

  if (layout.fixed.valid) {
    weight[FinalPairIndex(layout.fixed.leg[0], layout.fixed.leg[1])] = total_weight;
    return weight;
  }

  if (!layout.photo.empty()) {
    if (helicity_norm.size() % layout.photo.size() != 0) {
      throw AmplitudeFailure("FinalSectorWeights: invalid photonuclear component stride");
    }
    for (const auto& h : indices(helicity_norm)) {
      const auto& channel = layout.photo[h % layout.photo.size()];
      weight[channel.Pair()] += helicity_norm[h];
    }
    return weight;
  }

  if (layout.fusion.rows == 0) { return weight; }
  const std::size_t upper_count = layout.fusion.sector[0].size();
  const std::size_t lower_count = layout.fusion.sector[1].size();
  const std::size_t lower_rows  = layout.fusion.rows * lower_count;
  const std::size_t upper_rows  = layout.fusion.rows * upper_count;
  const std::size_t block       = upper_rows * lower_rows;
  if (upper_count == 0 || lower_count == 0 || block == 0 || helicity_norm.size() % block != 0) {
    throw AmplitudeFailure("FinalSectorWeights: invalid photon-fusion component stride");
  }
  for (const auto& h : indices(helicity_norm)) {
    const std::size_t slot  = h % block;
    const std::size_t upper = (slot / lower_rows) / layout.fusion.rows;
    const std::size_t lower = (slot % lower_rows) / layout.fusion.rows;
    weight[FinalPairIndex(layout.fusion.sector[0][upper], layout.fusion.sector[1][lower])] += helicity_norm[h];
  }
  return weight;
}

// Select one orthogonal final sector from finite nonnegative weights
// P(f) = W_f / sum_g W_g
FinalState SampleFinalState(const FinalWeights& weight, const double unit) {
  if (!std::isfinite(unit) || unit < 0.0 || unit > 1.0) {
    throw AmplitudeFailure("SampleFinalState: invalid random coordinate");
  }
  double      total    = 0.0;
  std::size_t last     = 0;
  bool        positive = false;
  for (const auto& i : indices(weight)) {
    if (!std::isfinite(weight[i]) || weight[i] < 0.0) {
      throw AmplitudeFailure("SampleFinalState: weights must be finite and nonnegative");
    }
    total += weight[i];
    if (weight[i] > 0.0) {
      last     = i;
      positive = true;
    }
  }
  if (!positive || !(total > 0.0) || !std::isfinite(total)) { return {}; }

  const double target     = unit * total;
  double       cumulative = 0.0;
  std::size_t  selected   = last;
  for (const auto& i : indices(weight)) {
    cumulative += weight[i];
    if (target < cumulative) {
      selected = i;
      break;
    }
  }
  FinalState state;
  state.valid  = true;
  state.leg[0] = selected / 2 == 0 ? CoherenceType::Coherent : CoherenceType::Incoherent;
  state.leg[1] = selected % 2 == 0 ? CoherenceType::Coherent : CoherenceType::Incoherent;
  return state;
}

// Initialize one convolution from its immutable model and Born point
MUPCScreen::MUPCScreen(const MUPC& upc, ScreenLayout layout, ScreenPoint born)
    : upc_(upc),
      layout_(std::move(layout)),
      flow_count_(born.flow.size()),
      born_size_(born.amplitude.size()),
      sample_count_(upc.SampleCount()) {
  if (born_size_ == 0) { throw AmplitudeFailure("MUPCScreen: empty Born amplitude"); }
  photo_terms_ = !born.photo_terms.empty();
  for (const auto& leg : indices(photo_current_)) { photo_current_[leg] = born.photo_current[leg].has_value(); }
  ValidatePoint(born);

  switch (layout_.type) {
    case ScreenType::Scalar:
      PrepareScalar(born);
      return;
    case ScreenType::Fusion:
      PrepareFusion(born);
      return;
    case ScreenType::Photo:
      PreparePhoto(born);
      return;
  }
  throw AmplitudeFailure("MUPCScreen: invalid screening type");
}

// Validate a complete point before changing the accumulated amplitude
void MUPCScreen::ValidatePoint(const ScreenPoint& point) const {
  if (point.amplitude.size() != born_size_ || !gra::AllFinite(point.amplitude) || point.flow.size() != flow_count_) {
    throw AmplitudeFailure("MUPCScreen: amplitude layout changed");
  }
  for (const auto& flow : point.flow) {
    if (flow.size() != born_size_ || !gra::AllFinite(flow)) {
      throw AmplitudeFailure("MUPCScreen: invalid color amplitude");
    }
  }
  for (const auto& q : point.transfer) {
    if (!std::isfinite(q[0]) || !std::isfinite(q[1]) || !std::isfinite(q[2])) {
      throw AmplitudeFailure("MUPCScreen: non-finite momentum transfer");
    }
  }
  for (const auto& leg : indices(photo_current_)) {
    if (point.photo_current[leg].has_value() != photo_current_[leg]) {
      throw AmplitudeFailure("MUPCScreen: nuclear target current layout changed");
    }
  }
  ValidateCurrents(point.photo_current, upc_);
  if (layout_.type == ScreenType::Photo || (layout_.type == ScreenType::Scalar && !layout_.photo.empty())) {
    ValidatePhotoTerms(point, born_size_, layout_.photo, photo_terms_, upc_);
  } else if (!point.photo_terms.empty()) {
    throw AmplitudeFailure("MUPCScreen: unexpected photonuclear terms");
  }
}

// Add one shifted hard amplitude at a precomputed Glauber node
void MUPCScreen::Add(const LoopNode& node, const ScreenPoint& point) {
  ValidatePoint(point);
  if (!std::isfinite(node.kt) || !std::isfinite(node.phi) || !std::isfinite(node.kx) || !std::isfinite(node.ky) ||
      !std::isfinite(node.weight.real()) || !std::isfinite(node.weight.imag()) ||
      (!node.sample_weight.empty() && node.sample_weight.size() != sample_count_) ||
      !gra::AllFinite(node.sample_weight)) {
    throw AmplitudeFailure("MUPCScreen::Add: invalid convolution weight");
  }
  switch (layout_.type) {
    case ScreenType::Scalar:
      AddScalar(node, point);
      return;
    case ScreenType::Fusion:
      AccumulateFusionPoint(point, node.weight, node.sample_weight);
      return;
    case ScreenType::Photo:
      AccumulatePhotoPoint(point, node.weight, node.sample_weight);
      return;
  }
  throw AmplitudeFailure("MUPCScreen::Add: invalid screening type");
}

// Project the accumulated ensemble into resolved physical sectors
ScreenResult MUPCScreen::Result() const {
  switch (layout_.type) {
    case ScreenType::Scalar:
      return ScalarResult();
    case ScreenType::Fusion:
      return EnsembleResult();
    case ScreenType::Photo:
      return EnsembleResult();
  }
  throw AmplitudeFailure("MUPCScreen::Result: invalid screening type");
}

// Initialize the scalar screening convolution
// A_scr^(0) = W_Born A(0)
void MUPCScreen::PrepareScalar(const ScreenPoint& born) {
  born_bank_ = {born.amplitude};
  born_bank_.insert(born_bank_.end(), born.flow.begin(), born.flow.end());
  scalar_bank_ = born_bank_;
  for (auto& bank : scalar_bank_) { gra::Scale(bank, upc_.BornWeight()); }
  if (!layout_.photo.empty() && layout_.photo != PhotoChannels(upc_)) {
    throw AmplitudeFailure("MUPCScreen: noncanonical photonuclear layout");
  }
}

// Initialize the resolved two-photon configuration convolution
// A_c,h(0) = sum_ab A_ab,h(0) R_a,c(0) R_b,c(0)
void MUPCScreen::PrepareFusion(const ScreenPoint& born) {
  if (sample_count_ == 0 || !upc_.HasSamples() || layout_.fusion.rows == 0) {
    throw AmplitudeFailure("MUPCScreen: missing fusion ensemble");
  }
  for (std::size_t leg = 0; leg < 2; ++leg) {
    if (layout_.fusion.sector[leg] !=
        CoherenceSectors(upc_.SourceSector(static_cast<int>(leg + 1), upc_.Param().emission[leg]))) {
      throw AmplitudeFailure("MUPCScreen: noncanonical fusion sectors");
    }
  }
  const std::size_t upper_rows = layout_.fusion.rows * layout_.fusion.sector[0].size();
  const std::size_t lower_rows = layout_.fusion.rows * layout_.fusion.sector[1].size();
  const std::size_t block      = upper_rows * lower_rows;
  if (block == 0 || born_size_ % block != 0) { throw AmplitudeFailure("MUPCScreen: invalid fusion amplitude stride"); }

  hard_count_ = born_size_ / (layout_.fusion.sector[0].size() * layout_.fusion.sector[1].size());
  for (const auto upper : CoherenceSectors(upc_.Param().emission[0])) {
    for (const auto lower : CoherenceSectors(upc_.Param().emission[1])) {
      selected_[0][FinalPairIndex(upper, lower)] = true;
    }
  }
  sample_bank_.assign(1 + born.flow.size(), std::vector<std::vector<std::complex<double>>>(
                                                sample_count_, std::vector<std::complex<double>>(hard_count_, 0.0)));
  AccumulateFusionPoint(born, upc_.BornWeight());
}

// Initialize the resolved photonuclear configuration convolution
// A_c,h(0) = sum_ab A_ab,h(0) R_emit,a,c(0) R_target,b,c(0)
void MUPCScreen::PreparePhoto(const ScreenPoint& born) {
  if (sample_count_ == 0 || !upc_.HasSamples() || layout_.photo.empty() || layout_.photo != PhotoChannels(upc_)) {
    throw AmplitudeFailure("MUPCScreen: invalid photonuclear ensemble");
  }
  channel_count_ = layout_.photo.size();
  if (born_size_ % channel_count_ != 0) { throw AmplitudeFailure("MUPCScreen: invalid photonuclear amplitude stride"); }
  // Direction-dependent target shadowing requires the same decomposition for every color projection
  if (photo_terms_ && !born.flow.empty()) {
    throw AmplitudeFailure("MUPCScreen: direction-resolved photonuclear color amplitudes are absent");
  }
  hard_count_ = born_size_ / channel_count_;
  for (const auto& channel : PhotoChannels(upc_.Param())) {
    selected_[channel.direction == PhotoDirection::Upper ? 0U : 1U][channel.Pair()] = true;
  }
  if (upc_.Param().reaction && upc_.Param().additional_emd) {
    for (auto& amplitude : photo_bank_) { amplitude.resize(hard_count_, 0.0); }
  }
  const std::size_t pair_count = 4;  // CC, CI, IC and II sectors
  sample_bank_.assign(1 + born.flow.size(),
                      std::vector<std::vector<std::complex<double>>>(
                          sample_count_, std::vector<std::complex<double>>(hard_count_ * pair_count, 0.0)));
  AccumulatePhotoPoint(born, upc_.BornWeight());
}

// Add one scalar Glauber contribution
// A_scr += w(k) [A(k) - A(0)]
void MUPCScreen::AddScalar(const LoopNode& node, const ScreenPoint& point) {
  scalar_weight_ += node.weight;
  for (const auto& h : indices(point.amplitude)) {
    scalar_bank_[0][h] += node.weight * (point.amplitude[h] - born_bank_[0][h]);
  }
  for (const auto& flow : indices(point.flow)) {
    for (const auto& h : indices(point.flow[flow])) {
      scalar_bank_[flow + 1][h] += node.weight * (point.flow[flow][h] - born_bank_[flow + 1][h]);
    }
  }
}

// Reconstruct each eigen-amplitude before applying the configuration-dependent kernel
// A_c,h += w_c(k) sum_ab A_ab,h(k) R_a,c(k) R_b,c(k), R_C = 1, R_I = delta J / sqrt(Var J)
void MUPCScreen::AccumulateFusionPoint(const ScreenPoint& point, const std::complex<double> common,
                                       const std::vector<std::complex<double>>& weight) {
  const auto upper = EmissionRatios(upc_, 1, layout_.fusion.sector[0], point.transfer[0]);
  const auto lower = EmissionRatios(upc_, 2, layout_.fusion.sector[1], point.transfer[1]);
  for (std::size_t sample = 0; sample < sample_count_; ++sample) {
    const std::complex<double> kernel = weight.empty() ? common : weight[sample];
    for (const auto& h : indices(point.amplitude)) {
      const auto [i, j]                 = FusionPair(layout_.fusion, h);
      const std::complex<double> factor = upper[i][sample] * lower[j][sample] * kernel;
      const std::size_t          row    = FusionRow(layout_.fusion, h);
      sample_bank_[0][sample][row] += factor * point.amplitude[h];
      for (const auto& flow : indices(point.flow)) {
        sample_bank_[flow + 1][sample][row] += factor * point.flow[flow][h];
      }
    }
  }
}

// Compute all sample ratios for the resolved photonuclear channels
// R_c = R_emit,c R_target,c
std::vector<std::vector<std::complex<double>>> MUPCScreen::PhotoRatios(
    const std::array<M3Vec, 2>& transfer, const std::array<std::optional<PhotoCurrent>, 2>& current,
    const std::optional<PhotoDirection> direction) const {
  using CurrentBank  = std::array<std::array<std::vector<std::complex<double>>, 2>, 2>;
  using CurrentReady = std::array<std::array<bool, 2>, 2>;
  CurrentBank                                    emission_bank;
  CurrentBank                                    target_bank;
  CurrentReady                                   emission_ready{};
  CurrentReady                                   target_ready{};
  std::vector<std::vector<std::complex<double>>> ratio(channel_count_,
                                                       std::vector<std::complex<double>>(sample_count_, 0.0));

  // Compute the bank index of one resolved sector
  const auto sector_index = [](const CoherenceType sector) { return sector == CoherenceType::Coherent ? 0U : 1U; };

  for (const auto& channel : indices(layout_.photo)) {
    const auto& item = layout_.photo[channel];
    if (direction.has_value() && item.direction != *direction) { continue; }
    const std::size_t emitter_leg = item.direction == PhotoDirection::Upper ? 0U : 1U;
    const std::size_t target_leg  = 1U - emitter_leg;
    const std::size_t source      = sector_index(item.emission);
    const std::size_t transition  = sector_index(item.target);
    if (!emission_ready[emitter_leg][source]) {
      emission_bank[emitter_leg][source] =
          upc_.EmissionRatios(static_cast<int>(emitter_leg + 1), item.emission, transfer[emitter_leg]);
      emission_ready[emitter_leg][source] = true;
    }
    if (!target_ready[target_leg][transition]) {
      const auto& target_current = current[target_leg];
      if (target_current.has_value()) {
        target_bank[target_leg][transition] =
            upc_.TargetCurrentRatios(static_cast<int>(target_leg + 1), item.target, *target_current);
      } else {
        target_bank[target_leg][transition] =
            upc_.TargetRatios(static_cast<int>(target_leg + 1), item.target, transfer[target_leg]);
      }
      target_ready[target_leg][transition] = true;
    }
    const auto& emitter = emission_bank[emitter_leg][source];
    const auto& target  = target_bank[target_leg][transition];
    for (std::size_t sample = 0; sample < sample_count_; ++sample) {
      ratio[channel][sample] = emitter[sample] * target[sample];
    }
  }
  return ratio;
}

// Accumulate separately shadowed terms before the Good-Walker projection
void MUPCScreen::AccumulatePhotoPoint(const ScreenPoint& point, const std::complex<double> common,
                                      const std::vector<std::complex<double>>& weight) {
  if (!photo_terms_) {
    AccumulatePhoto(point.amplitude, PhotoRatios(point.transfer, point.photo_current, std::nullopt), common, weight, 0);
  } else {
    for (const auto& term : point.photo_terms) {
      AccumulatePhoto(term.amplitude, PhotoRatios(point.transfer, term.photo_current, term.direction), common, weight,
                      0);
    }
  }
  if (!point.flow.empty()) {
    const auto ratio = PhotoRatios(point.transfer, point.photo_current, std::nullopt);
    for (const auto& flow : indices(point.flow)) { AccumulatePhoto(point.flow[flow], ratio, common, weight, flow + 1); }
  }
}

// Reconstruct each photon direction before the final nuclear projections
// A_c,f += w_c sum_ab A_ab R_a,c R_b,c for every selected final sector f
void MUPCScreen::AccumulatePhoto(const std::vector<std::complex<double>>&              amplitude,
                                 const std::vector<std::vector<std::complex<double>>>& ratio,
                                 const std::complex<double> common, const std::vector<std::complex<double>>& weight,
                                 const std::size_t bank) {
  if (bank >= sample_bank_.size()) { throw AmplitudeFailure("MUPCScreen::AccumulatePhoto: invalid bank"); }
  const std::size_t pair_count = 4;  // CC, CI, IC and II sectors
  for (const auto& h : indices(amplitude)) {
    const std::size_t channel   = h % channel_count_;
    const std::size_t hard      = h / channel_count_;
    const auto        direction = layout_.photo[channel].direction == PhotoDirection::Upper ? 0U : 1U;
    for (std::size_t sample = 0; sample < sample_count_; ++sample) {
      const auto kernel = weight.empty() ? common : weight[sample];
      const auto value  = kernel * ratio[channel][sample] * amplitude[h];
      for (const auto& pair : indices(selected_[direction])) {
        if (selected_[direction][pair]) { sample_bank_[bank][sample][hard * pair_count + pair] += value; }
      }
      if (bank == 0 && selected_[direction][0] && !photo_bank_[0].empty()) {
        photo_bank_[direction][hard] += value / static_cast<double>(sample_count_);
      }
    }
  }
}

// Compute the scalar screening result
// A_scr = A(0) + sum_k w(k) A(k)
// W_h = N_avg |A_scr,h|^2
ScreenResult MUPCScreen::ScalarResult() const {
  ScreenResult result;
  auto         bank = scalar_bank_;
  for (const auto& i : indices(bank)) {
    for (const auto& h : indices(bank[i])) { bank[i][h] += scalar_weight_ * born_bank_[i][h]; }
  }
  if (!layout_.photo.empty() && upc_.Param().reaction && upc_.Param().additional_emd) {
    std::array<std::vector<std::complex<double>>, 2> direction;
    for (auto& amplitude : direction) { amplitude.resize(born_size_ / layout_.photo.size(), 0.0); }
    for (const auto& h : indices(bank[0])) {
      const auto& channel = layout_.photo[h % layout_.photo.size()];
      if (channel.Pair() != 0) { continue; }
      direction[channel.direction == PhotoDirection::Upper ? 0U : 1U][h / layout_.photo.size()] += bank[0][h];
    }
    result.photo = PhotoInterference(direction);
  }
  if (!layout_.photo.empty()) {
    for (auto& amplitude : bank) { amplitude = CombinePhotoChannels(amplitude, layout_.photo); }
  }
  result.amplitude = std::move(bank.front());
  result.helicity_norm.resize(result.amplitude.size(), 0.0);
  for (const auto& h : indices(result.amplitude)) { result.helicity_norm[h] = std::norm(result.amplitude[h]); }
  result.color_norm.resize(bank.size() - 1, 0.0);
  result.color_sector.resize(result.color_norm.size());
  for (const auto& i : indices(result.color_norm)) {
    result.color_norm[i] = gra::SquaredNorm(bank[i + 1]);
    std::vector<double> weight(bank[i + 1].size(), 0.0);
    for (const auto& h : indices(weight)) { weight[h] = std::norm(bank[i + 1][h]); }
    result.color_sector[i] = FinalSectorWeights(layout_, weight);
  }
  return result;
}

// Project either fusion or photonuclear ensembles through the same Good-Walker decomposition
// Only the coherent sector has one complex amplitude, excited sectors retain their norms
ScreenResult MUPCScreen::EnsembleResult() const {
  ScreenResult result;
  result.amplitude.resize(born_size_, 0.0);
  result.helicity_norm.resize(born_size_, 0.0);
  result.color_norm.resize(sample_bank_.size() - 1, 0.0);
  result.color_sector.resize(result.color_norm.size());
  const bool photo          = layout_.type == ScreenType::Photo;
  const auto representative = PhotoRepresentatives(layout_.photo);
  if (photo) { result.photo = PhotoInterference(photo_bank_); }
  std::vector<std::complex<double>> ensemble(sample_count_, 0.0);
  for (std::size_t h = 0; h < born_size_; ++h) {
    std::size_t pair = 0, component = h;
    if (photo) {
      const auto channel = h % channel_count_;
      pair               = layout_.photo[channel].Pair();
      if (representative[pair] != channel) { continue; }
      component = (h / channel_count_) * 4 + pair;
    } else {
      const auto [upper, lower] = FusionPair(layout_.fusion, h);
      pair                      = FinalPairIndex(layout_.fusion.sector[0][upper], layout_.fusion.sector[1][lower]);
      if (!selected_[0][pair]) { continue; }
      component = FusionRow(layout_.fusion, h);
    }
    for (const auto& bank : indices(sample_bank_)) {
      for (const auto& sample : indices(ensemble)) { ensemble[sample] = sample_bank_[bank][sample][component]; }
      const auto value = NuclearGoodWalkerProject(ensemble, upc_.SampleShape())[pair / 2][pair % 2];
      if (bank == 0) {
        result.helicity_norm[h] = value;
        if (pair == 0) {
          result.amplitude[h] = std::accumulate(ensemble.begin(), ensemble.end(), std::complex<double>{}) /
                                static_cast<double>(sample_count_);
        }
      } else {
        result.color_norm[bank - 1] += value;
        result.color_sector[bank - 1][pair] += value;
      }
    }
  }
  return result;
}

}  // namespace gra::nuclear
