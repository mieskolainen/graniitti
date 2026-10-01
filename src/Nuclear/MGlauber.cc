// Finite-range optical, GGCF and configuration Glauber survival
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MGlauber.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {
namespace {

// Compute the exact interpolation support of one finite NN profile
double ProfileSupport(const NNProfile& profile) {
  std::size_t last = 0;
  for (const auto& i : indices(profile.inelastic)) {
    if (profile.inelastic[i] > 0.0) { last = i; }
  }
  return profile.b_node[std::min(last + 1, profile.b_node.size() - 1)];
}

// Visit nucleons in one transverse interval, treating a proton as a center at zero
// Stop once every fluctuation state has zero survival probability
template <class F>
bool VisitNucleons(const MConfig* config, double xmin, double xmax, F&& visit) {
  if (config == nullptr) { return xmin > 0.0 || xmax < 0.0 || visit(std::array<double, 3>{}); }
  const auto& nucleons = config->Nucleons();
  const auto& order    = config->TransverseOrder();
  const auto  begin =
      std::lower_bound(order.begin(), order.end(), xmin, [&](std::size_t i, double x) { return nucleons[i].x[0] < x; });
  const auto end =
      std::upper_bound(begin, order.end(), xmax, [&](double x, std::size_t i) { return x < nucleons[i].x[0]; });
  for (auto it = begin; it != end; ++it) {
    if (!visit(nucleons[*it].x)) { return false; }
  }
  return true;
}

}  // namespace

// Construct and tabulate one heterogeneous hadronic overlap
MGlauber::MGlauber(std::array<HadronType, 2> type, std::array<std::shared_ptr<const MNucleus>, 2> nucleus,
                   GlauberParam param)
    : type_(type), nucleus_(std::move(nucleus)), param_(param) {
  param_.fluctuation.omega = param_.profile.omega;
  Validate();
  PrepareProfile();
  PrepareSigma();
  Prepare();
}

// Construct and tabulate one finite-range nuclear overlap
MGlauber::MGlauber(const MNucleus& beam1, const MNucleus& beam2, GlauberParam param)
    : MGlauber({HadronType::Nucleus, HadronType::Nucleus},
               {std::make_shared<const MNucleus>(beam1), std::make_shared<const MNucleus>(beam2)}, param) {}

// Validate all physical and numerical controls
void MGlauber::Validate() const {
  for (const auto& i : indices(type_)) {
    if (type_[i] == HadronType::Nucleus && nucleus_[i] == nullptr) {
      throw std::invalid_argument("MGlauber: a nuclear beam requires a nuclear model");
    }
    if (type_[i] == HadronType::Proton && nucleus_[i] != nullptr) {
      throw std::invalid_argument("MGlauber: a proton beam must not carry a nuclear model");
    }
  }
  if (!std::isfinite(param_.profile.sigma) || !std::isfinite(param_.profile.omega) || !std::isfinite(param_.b_max) ||
      !std::isfinite(param_.q_max) || !(param_.profile.sigma > 0.0) || param_.profile.omega < 0.0 ||
      !(param_.b_max > 0.0) || !(param_.q_max > 0.0) || param_.profile.fingerprint.empty() ||
      param_.profile.b_node.size() < 2 || param_.profile.b_node.size() != param_.profile.inelastic.size()) {
    throw std::invalid_argument("MGlauber: invalid physical controls");
  }
  if (param_.b_nodes < 32 || param_.b_nodes > 4096 || param_.q_nodes < 32 || param_.q_nodes > 4096 ||
      param_.profile_nodes < 32 || param_.profile_nodes > 16384 || param_.profile_nodes % 2 != 0 ||
      param_.fluctuation.nodes == 0 || param_.fluctuation.nodes > 256) {
    throw std::invalid_argument("MGlauber: invalid quadrature node count");
  }
  for (const auto& i : indices(param_.profile.b_node)) {
    if (!std::isfinite(param_.profile.b_node[i]) || !std::isfinite(param_.profile.inelastic[i]) ||
        param_.profile.b_node[i] < 0.0 || param_.profile.inelastic[i] < 0.0 || param_.profile.inelastic[i] > 1.0 ||
        (i > 0 && !(param_.profile.b_node[i] > param_.profile.b_node[i - 1]))) {
      throw std::invalid_argument("MGlauber: invalid pp eikonal profile");
    }
  }
  const double sigma =
      10.0 * 2.0 * math::PI * math::LinearRadialIntegral(param_.profile.b_node, param_.profile.inelastic);
  if (std::abs(sigma - param_.profile.sigma) > 1.0e-8 * param_.profile.sigma) {
    throw std::invalid_argument("MGlauber: pp profile normalization disagrees with sigma_NN");
  }
}

// Detect an equally spaced pp profile for direct interpolation lookup
void MGlauber::PrepareProfile() {
  const auto& node = param_.profile.b_node;
  profile_step_    = (node.back() - node.front()) / static_cast<double>(node.size() - 1);
  profile_uniform_ = profile_step_ > 0.0;
  for (const auto& i : indices(node)) {
    const double expected = node.front() + static_cast<double>(i) * profile_step_;
    // The epsilon multiplier only recognizes an exactly uniform input grid
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(expected));
    if (std::abs(node[i] - expected) > tolerance) {
      profile_uniform_ = false;
      break;
    }
  }
  if (profile_uniform_) { profile_inverse_step_ = 1.0 / profile_step_; }
}

// Prepare the moment-calibrated Gamma fluctuation and pair rules
// s_ij = f_i f_j with <f> = 1 and Var[f] = omega_N from one-sided diffraction
// [REFERENCE: Alvioli and Strikman, Phys. Lett. B722 (2013) 347]
void MGlauber::PrepareSigma() {
  sigma_scale_     = CrossSectionRule(param_.fluctuation);
  sigma_scale_max_ = *std::max_element(sigma_scale_.begin(), sigma_scale_.end());
  sigma_scale_max_ *= sigma_scale_max_;
  for (const auto& upper : indices(sigma_scale_)) {
    for (std::size_t lower = upper; lower < sigma_scale_.size(); ++lower) {
      pair_scale_.push_back(sigma_scale_[upper] * sigma_scale_[lower]);
      pair_radius_scale_.push_back(std::sqrt(pair_scale_.back()));
      pair_weight_.push_back(upper == lower ? 1.0 : 2.0);
    }
  }
}

// Compute the normalized inelastic nucleon form factor
// F_NN(q) = 2 pi / sigma_NN int b db P_inel(b) J_0(qb / hbarc)
double MGlauber::NNForm(const double q) const {
  constexpr double hbarc         = PDG::GeV2fm;
  constexpr double inverse_sqrt3 = 0.57735026918962576451;
  double           integral      = 0.0;
  for (std::size_t i = 1; i < param_.profile.b_node.size(); ++i) {
    const double lower      = param_.profile.b_node[i - 1];
    const double upper      = param_.profile.b_node[i];
    const double center     = 0.5 * (lower + upper);
    const double half_width = 0.5 * (upper - lower);
    for (const double sign : {-1.0, 1.0}) {
      const double b        = center + sign * half_width * inverse_sqrt3;
      const double fraction = (b - lower) / (upper - lower);
      const double probability =
          param_.profile.inelastic[i - 1] * (1.0 - fraction) + param_.profile.inelastic[i] * fraction;
      integral += half_width * b * probability * std::cyl_bessel_j(0.0, q * b / hbarc);
    }
  }
  return 2.0 * math::PI * integral / (0.1 * param_.profile.sigma);
}

// Prepare the finite-range overlap interpolation table
// T_12(b) = A_1 A_2 / (2 pi hbarc^2) int q dq F_1 F_2 F_NN J_0(qb / hbarc)
// [REFERENCE: Glauber, Lectures in Theoretical Physics 1 (1959) 315]
void MGlauber::Prepare() {
  const auto          q_rule = math::GaussLegendreRule(param_.q_nodes, 0.0, param_.q_max);
  std::vector<double> kernel(q_rule.first.size(), 0.0);
  ParallelFor(kernel.size(), [&](const std::size_t i) {
    const double q = q_rule.first[i];
    kernel[i]      = q_rule.second[i] * q * MatterForm(1, q) * MatterForm(2, q) * NNForm(q);
  });

  b_node_.resize(static_cast<std::size_t>(param_.b_nodes) + 1, 0.0);
  overlap_.resize(b_node_.size(), 0.0);
  constexpr double hbarc = PDG::GeV2fm;
  const double     scale = static_cast<double>(MassNumber(1)) * MassNumber(2) / (2.0 * math::PI * hbarc * hbarc);
  ParallelFor(b_node_.size(), [&](const std::size_t ib) {
    const double b  = param_.b_max * static_cast<double>(ib) / static_cast<double>(param_.b_nodes);
    b_node_[ib]     = b;
    double integral = 0.0;
    for (const auto& iq : indices(kernel)) {
      integral += kernel[iq] * std::cyl_bessel_j(0.0, q_rule.first[iq] * b / hbarc);
    }
    overlap_[ib] = std::max(0.0, scale * integral);
  });

  for (std::size_t i = 1; i < overlap_.size(); ++i) { overlap_[i] = std::min(overlap_[i], overlap_[i - 1]); }
}

// Compute one checked beam nucleus
const MNucleus* MGlauber::Nucleus(const int leg) const {
  if (leg == 1) { return nucleus_[0].get(); }
  if (leg == 2) { return nucleus_[1].get(); }
  throw std::out_of_range("MGlauber::Nucleus: leg must be one or two");
}

// Compute whether one checked beam uses the proton profile
bool MGlauber::IsProton(const int leg) const {
  if (leg < 1 || leg > 2) { throw std::out_of_range("MGlauber::IsProton: leg must be one or two"); }
  return type_[static_cast<std::size_t>(leg - 1)] == HadronType::Proton;
}

// Compute one beam mass number for overlap normalization
unsigned int MGlauber::MassNumber(const int leg) const {
  if (IsProton(leg)) { return 1; }
  return Nucleus(leg)->A();
}

// Compute one normalized beam matter form factor
// F_m(q) = int rho_m(r) exp(i q.r / hbarc) d^3r
double MGlauber::MatterForm(const int leg, const double q) const {
  if (IsProton(leg)) { return 1.0; }
  return Nucleus(leg)->MatterDensity().Form(q);
}

// Compute the finite-range nuclear overlap in fm^-2
// T_12(b) = int d^2s T_1(s) T_2(b - s) folded with P_NN
double MGlauber::Overlap(const double b) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MGlauber::Overlap: impact parameter must be finite and nonnegative");
  }
  if (b >= param_.b_max) { return 0.0; }
  return math::LinearInterpolateValidatedGrid(b_node_, overlap_, b);
}

// Compute the optical no-additional-interaction amplitude
// S_opt(b) = exp[-sigma_NN T_12(b) / 2]
double MGlauber::OpticalAmp(const double b) const {
  const double opacity = 0.05 * param_.profile.sigma * Overlap(b);
  return std::exp(-opacity);
}

// Compute the optical no-additional-interaction probability
// P_opt(b) = |S_opt(b)|^2
double MGlauber::OpticalProb(const double b) const {
  const double amplitude = OpticalAmp(b);
  return amplitude * amplitude;
}

// Compute the GGCF no-additional-interaction amplitude
// S_GGCF(b) = <exp[-s_ij sigma_NN T_12(b) / 2]>_ij
// [REFERENCE: Alvioli and Strikman, Phys. Lett. B722 (2013) 347]
double MGlauber::GGCFAmp(const double b) const {
  const double opacity   = 0.05 * param_.profile.sigma * Overlap(b);
  double       amplitude = 0.0;
  for (const auto& state : indices(pair_scale_)) {
    amplitude += pair_weight_[state] * std::exp(-opacity * pair_scale_[state]);
  }
  const double states = static_cast<double>(sigma_scale_.size());
  return amplitude / (states * states);
}

// Compute the GGCF no-additional-interaction probability
// P_GGCF(b) = |S_GGCF(b)|^2
double MGlauber::GGCFProb(const double b) const {
  const double amplitude = GGCFAmp(b);
  return amplitude * amplitude;
}

// Compute the scaled inelastic nucleon probability at one separation
// P_NN(b,s) = P_NN(b / sqrt(s),1)
double MGlauber::NNProb(const double b, const double sigma_scale) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MGlauber::NNProb: separation must be finite and nonnegative");
  }
  if (!std::isfinite(sigma_scale) || !(sigma_scale > 0.0)) {
    throw std::invalid_argument("MGlauber::NNProb: cross-section scale must be positive");
  }
  const double scaled_b = b / std::sqrt(sigma_scale);
  if (scaled_b >= param_.profile.b_node.back()) { return 0.0; }
  if (scaled_b <= param_.profile.b_node.front()) { return param_.profile.inelastic.front(); }
  return std::clamp(math::LinearInterpolateValidatedGrid(param_.profile.b_node, param_.profile.inelastic, scaled_b),
                    0.0, 1.0);
}

// Compute one prevalidated fluctuation-pair interaction probability
// P_ij(b) = P_NN(b / sqrt(s_ij),1)
double MGlauber::PairProb(const double b, const std::size_t state) const {
  const double scaled_b = b / pair_radius_scale_[state];
  const auto&  node     = param_.profile.b_node;
  const auto&  value    = param_.profile.inelastic;
  if (scaled_b >= node.back()) { return 0.0; }
  if (scaled_b <= node.front()) { return value.front(); }
  if (!profile_uniform_) { return std::clamp(math::LinearInterpolateValidatedGrid(node, value, scaled_b), 0.0, 1.0); }

  const double coordinate = (scaled_b - node.front()) * profile_inverse_step_;
  std::size_t  lower      = std::min(static_cast<std::size_t>(coordinate), node.size() - 2);
  while (lower > 0 && scaled_b < node[lower]) { --lower; }
  while (lower + 1 < node.size() - 1 && scaled_b >= node[lower + 1]) { ++lower; }
  const std::size_t upper    = lower + 1;
  const double      fraction = (scaled_b - node[lower]) / (node[upper] - node[lower]);
  return std::clamp(value[lower] * (1.0 - fraction) + value[upper] * fraction, 0.0, 1.0);
}

// Compute the finite-range nucleon no-interaction amplitude at one sigma
// S_NN(b,sigma) = sqrt[1 - P_NN(b,sigma / sigma_0)]
// This defines the positive scalar no-inelastic channel of the absorptive
// Glauber model and is not a reconstruction of the complex elastic pp S matrix
double MGlauber::NNAmp(const double b, const double sigma) const {
  if (!std::isfinite(b) || b < 0.0 || !std::isfinite(sigma) || !(sigma > 0.0)) {
    throw std::invalid_argument("MGlauber::NNAmp: impact parameter and sigma must be physical");
  }
  const double probability = NNProb(b, sigma / param_.profile.sigma);
  return std::sqrt(1.0 - probability);
}

// Accumulate one NN pair into the surviving fluctuation probabilities
void MGlauber::Attenuate(double radius, std::size_t state, std::span<double> probability, std::size_t& active) const {
  for (const auto& i : indices(probability)) {
    if (!(probability[i] > 0.0)) { continue; }
    probability[i] *= 1.0 - PairProb(radius, state < pair_scale_.size() ? state : i);
    if (probability[i] < std::numeric_limits<double>::min()) {
      probability[i] = 0.0;
      --active;
    }
  }
}

// Average no-interaction amplitudes over the requested fluctuation states
// S = sum_s w_s sqrt(P_s) / sum_s w_s
double MGlauber::Average(std::span<const double> probability, std::size_t state) const {
  if (state < pair_scale_.size()) { return std::sqrt(probability.front()); }
  double amplitude = 0.0;
  for (const auto& i : indices(probability)) { amplitude += pair_weight_[i] * std::sqrt(probability[i]); }
  return amplitude / math::pow2(static_cast<double>(sigma_scale_.size()));
}

// Compute one fixed-geometry amplitude averaged over fluctuation pairs
// S_c(b) = <prod_ij sqrt[1 - P_ij(|b + r_j - r_i|)]>_states
// [REFERENCE: Glauber, Lectures in Theoretical Physics 1 (1959) 315]
double MGlauber::PairAmp(const double bx, const double by, const MConfig* config1, const MConfig* config2,
                         const std::size_t state) const {
  if (!std::isfinite(bx) || !std::isfinite(by)) {
    throw std::invalid_argument("MGlauber::PairAmp: impact vector must be finite");
  }
  const std::size_t count1 = config1 == nullptr ? 1 : config1->Nucleons().size();
  const std::size_t count2 = config2 == nullptr ? 1 : config2->Nucleons().size();
  if (count1 > std::numeric_limits<std::size_t>::max() / count2) {
    throw std::invalid_argument("MGlauber::PairAmp: geometry is too large");
  }
  const bool                  fixed         = state < pair_scale_.size();
  const std::size_t           state_count   = fixed ? 1 : pair_scale_.size();
  const double                radius_scale  = fixed ? pair_radius_scale_[state] : std::sqrt(sigma_scale_max_);
  const double                radius_max    = ProfileSupport(param_.profile) * radius_scale;
  const double                radius2_max   = radius_max * radius_max;
  const std::array<double, 4> proton_bounds = {0.0, 0.0, 0.0, 0.0};
  const auto&                 bounds1       = config1 == nullptr ? proton_bounds : config1->TransverseBounds();
  const auto&                 bounds2       = config2 == nullptr ? proton_bounds : config2->TransverseBounds();
  const double                xmin          = bounds1[0] - bounds2[1];
  const double                xmax          = bounds1[1] - bounds2[0];
  const double                ymin          = bounds1[2] - bounds2[3];
  const double                ymax          = bounds1[3] - bounds2[2];
  const double                dx_box        = bx < xmin ? xmin - bx : (bx > xmax ? bx - xmax : 0.0);
  const double                dy_box        = by < ymin ? ymin - by : (by > ymax ? by - ymax : 0.0);
  if (dx_box * dx_box + dy_box * dy_box > radius2_max) { return 1.0; }
  constexpr std::size_t           local_count = 64;  // Allocation only, not physics
  std::array<double, local_count> local{};
  std::vector<double>             heap(state_count > local_count ? state_count : 0, 1.0);
  std::span<double> probability = heap.empty() ? std::span<double>(local.data(), state_count) : std::span<double>(heap);
  std::fill(probability.begin(), probability.end(), 1.0);
  std::size_t active_count = state_count;

  VisitNucleons(config1, bounds2[0] + bx - radius_max, bounds2[1] + bx + radius_max, [&](const auto& x1) {
    const double center = x1[0] - bx;
    return VisitNucleons(config2, center - radius_max, center + radius_max, [&](const auto& x2) {
      const double dx        = x1[0] - x2[0] - bx;
      const double dy        = x1[1] - x2[1] - by;
      const double distance2 = dx * dx + dy * dy;
      if (distance2 <= radius2_max) { Attenuate(std::sqrt(distance2), state, probability, active_count); }
      return active_count > 0;
    });
  });
  return Average(probability, state);
}

// Validate one pair of optional configuration banks
void MGlauber::ValidateBanks(const MConfigBank* bank1, const MConfigBank* bank2) const {
  const std::array<const MConfigBank*, 2> bank = {bank1, bank2};
  for (const auto& i : indices(type_)) {
    if (type_[i] == HadronType::Nucleus &&
        (bank[i] == nullptr || bank[i]->Nucleus().ID().pdg != nucleus_[i]->ID().pdg)) {
      throw std::invalid_argument("MGlauber: nuclear configuration bank does not match");
    }
    if (type_[i] == HadronType::Proton && bank[i] != nullptr) {
      throw std::invalid_argument("MGlauber: proton leg must not carry a nuclear bank");
    }
  }
}

// Select a validated nuclear configuration or a point proton on each leg
std::array<const MConfig*, 2> MGlauber::Configs(const MConfigBank* bank1, const MConfigBank* bank2, std::size_t config1,
                                                std::size_t config2) const {
  ValidateBanks(bank1, bank2);
  if (config1 >= (bank1 ? bank1->Size() : 1) || config2 >= (bank2 ? bank2->Size() : 1)) {
    throw std::out_of_range("MGlauber::ConfigPairAmp: configuration is out of range");
  }
  return {bank1 ? &bank1->At(config1) : nullptr, bank2 ? &bank2->At(config2) : nullptr};
}

// Compute the number of Cartesian configuration-pair samples
std::size_t MGlauber::SampleCount(const MConfigBank* bank1, const MConfigBank* bank2) const {
  ValidateBanks(bank1, bank2);
  const std::size_t count1 = bank1 == nullptr ? 1 : bank1->Size();
  const std::size_t count2 = bank2 == nullptr ? 1 : bank2->Size();
  if (count1 > std::numeric_limits<std::size_t>::max() / count2) {
    throw std::invalid_argument("MGlauber::SampleCount: Cartesian ensemble is too large");
  }
  return count1 * count2;
}

// Compute one leg cross section averaged over the other leg [mb]
// <sigma_ij>_j = sigma_NN f_i
double MGlauber::Sigma(const std::size_t node) const { return param_.profile.sigma * sigma_scale_.at(node); }

// Compute the physical pair cross section at two independent fluctuation nodes [mb]
double MGlauber::Sigma(const std::size_t upper, const std::size_t lower) const {
  return param_.profile.sigma * sigma_scale_.at(upper) * sigma_scale_.at(lower);
}

// Compute one paired configuration survival amplitude
// S_c(b) = <prod_ij sqrt[1 - P_ij]>_states
double MGlauber::SampleAmp(const double b, const MConfigBank* bank1, const MConfigBank* bank2,
                           const std::size_t sample) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument(
        "MGlauber::SampleAmp: impact parameter must be finite and "
        "nonnegative");
  }
  return SampleAmp(b, 0.0, bank1, bank2, sample);
}

// Compute one paired configuration amplitude at a transverse vector
// S_c(b_x,b_y) = <prod_ij sqrt[1 - P_ij]>_states
double MGlauber::SampleAmp(const double bx, const double by, const MConfigBank* bank1, const MConfigBank* bank2,
                           const std::size_t sample) const {
  if (!std::isfinite(bx) || !std::isfinite(by)) {
    throw std::invalid_argument("MGlauber::SampleAmp: impact vector must be finite");
  }
  const std::size_t count = SampleCount(bank1, bank2);
  if (sample >= count) { throw std::out_of_range("MGlauber::SampleAmp: sample is out of range"); }
  const std::size_t count2  = bank2 == nullptr ? 1 : bank2->Size();
  const std::size_t config1 = sample / count2;
  const std::size_t config2 = sample % count2;
  return ConfigPairAmp(bx, by, bank1, bank2, config1, config2);
}

// Compute one configuration-pair amplitude averaged over GGCF fluctuations
// S_c(b) = <S_c(b,s_ij)>_states
double MGlauber::ConfigPairAmp(const double bx, const double by, const MConfigBank* bank1, const MConfigBank* bank2,
                               const std::size_t config1, const std::size_t config2) const {
  if (!std::isfinite(bx) || !std::isfinite(by)) {
    throw std::invalid_argument("MGlauber::ConfigPairAmp: impact vector must be finite");
  }
  const auto [first, second] = Configs(bank1, bank2, config1, config2);
  return PairAmp(bx, by, first, second, pair_scale_.size());
}

// Compute the packed symmetric GGCF pair-state index
std::size_t MGlauber::PairState(const std::size_t sigma1, const std::size_t sigma2) const {
  const std::size_t count = sigma_scale_.size();
  if (sigma1 >= count || sigma2 >= count) {
    throw std::out_of_range("MGlauber::PairState: eigenstate is out of range");
  }
  const std::size_t upper = std::min(sigma1, sigma2);
  const std::size_t lower = std::max(sigma1, sigma2);
  return upper * (2 * count - upper + 1) / 2 + lower - upper;
}

// Compute one configuration-pair amplitude at fixed GGCF eigenstates
// S_c(b,s_1,s_2) = prod_ij sqrt[1 - P_ij(b,s_1 s_2)]
double MGlauber::ConfigPairAmp(const double bx, const double by, const MConfigBank* bank1, const MConfigBank* bank2,
                               const std::size_t config1, const std::size_t config2, const std::size_t sigma1,
                               const std::size_t sigma2) const {
  if (!std::isfinite(bx) || !std::isfinite(by)) {
    throw std::invalid_argument("MGlauber::ConfigPairAmp: impact vector must be finite");
  }
  const auto [first, second] = Configs(bank1, bank2, config1, config2);
  return PairAmp(bx, by, first, second, PairState(sigma1, sigma2));
}

// Compute one GGCF-averaged configuration amplitude at many impact vectors
// S_c(b_k) = <S_c(b_k,s_ij)>_states
std::vector<double> MGlauber::ConfigPairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                                            const MConfigBank* bank2, const std::size_t config1,
                                            const std::size_t config2) const {
  return PairAmp(impact, bank1, bank2, config1, config2, pair_scale_.size());
}

// Compute one fixed-state configuration amplitude at many impact vectors
// S_c(b_k,s_1,s_2) = prod_ij sqrt[1 - P_ij(b_k,s_1 s_2)]
std::vector<double> MGlauber::ConfigPairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                                            const MConfigBank* bank2, const std::size_t config1,
                                            const std::size_t config2, const std::size_t sigma1,
                                            const std::size_t sigma2) const {
  return PairAmp(impact, bank1, bank2, config1, config2, PairState(sigma1, sigma2));
}

// Compute fixed or GGCF-averaged amplitudes using one spatial pair grid
// S_c(b_k) = <prod_ij sqrt[1 - P_ij(|b_k + r_j - r_i|)]>_states
std::vector<double> MGlauber::PairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                                      const MConfigBank* bank2, const std::size_t config1, const std::size_t config2,
                                      const std::size_t state) const {
  const auto [first, second] = Configs(bank1, bank2, config1, config2);
  for (const auto& b : impact) {
    if (!gra::AllFinite(b)) { throw std::invalid_argument("MGlauber::ConfigPairAmp: impact vector must be finite"); }
  }

  const std::array<double, 3> proton        = {0.0, 0.0, 0.0};
  const std::array<double, 4> proton_bounds = {0.0, 0.0, 0.0, 0.0};
  const auto&                 bounds1       = first == nullptr ? proton_bounds : first->TransverseBounds();
  const auto&                 bounds2       = second == nullptr ? proton_bounds : second->TransverseBounds();
  const std::array<double, 4> bounds = {bounds1[0] - bounds2[1], bounds1[1] - bounds2[0], bounds1[2] - bounds2[3],
                                        bounds1[3] - bounds2[2]};

  const bool        fixed        = state < pair_scale_.size();
  const std::size_t state_count  = fixed ? 1 : pair_scale_.size();
  const double      radius_scale = fixed ? pair_radius_scale_[state] : std::sqrt(sigma_scale_max_);
  const double      radius       = ProfileSupport(param_.profile) * radius_scale;
  const double      radius2      = radius * radius;
  constexpr long    cell_span    = 4;  // Exact lookup acceleration only
  const double      cell_width   = radius / static_cast<double>(cell_span);
  const std::size_t nx =
      std::max<std::size_t>(1, static_cast<std::size_t>(std::floor((bounds[1] - bounds[0]) / cell_width)) + 1);
  const std::size_t ny =
      std::max<std::size_t>(1, static_cast<std::size_t>(std::floor((bounds[3] - bounds[2]) / cell_width)) + 1);
  if (nx > std::numeric_limits<std::size_t>::max() / ny) {
    throw std::invalid_argument("MGlauber::ConfigPairAmp: grid is too large");
  }
  std::vector<std::vector<std::array<double, 2>>> cell(nx * ny);
  const auto insert = [&](const std::array<double, 3>& x1, const std::array<double, 3>& x2) {
    const std::array<double, 2> point = {x1[0] - x2[0], x1[1] - x2[1]};
    const std::size_t           ix    = std::min(nx - 1, static_cast<std::size_t>((point[0] - bounds[0]) / cell_width));
    const std::size_t           iy    = std::min(ny - 1, static_cast<std::size_t>((point[1] - bounds[2]) / cell_width));
    cell[iy * nx + ix].push_back(point);
  };
  if (first == nullptr && second == nullptr) {
    insert(proton, proton);
  } else if (first == nullptr) {
    for (const auto& item : second->Nucleons()) { insert(proton, item.x); }
  } else if (second == nullptr) {
    for (const auto& item : first->Nucleons()) { insert(item.x, proton); }
  } else {
    for (const auto& upper : first->Nucleons()) {
      for (const auto& lower : second->Nucleons()) { insert(upper.x, lower.x); }
    }
  }

  std::vector<double> amplitude(impact.size(), 1.0);
  std::vector<double> pair_amplitude(state_count, 1.0);
  for (const auto& i : indices(impact)) {
    const auto&  b      = impact[i];
    const double dx_box = b[0] < bounds[0] ? bounds[0] - b[0] : (b[0] > bounds[1] ? b[0] - bounds[1] : 0.0);
    const double dy_box = b[1] < bounds[2] ? bounds[2] - b[1] : (b[1] > bounds[3] ? b[1] - bounds[3] : 0.0);
    if (dx_box * dx_box + dy_box * dy_box > radius2) { continue; }
    const long ix = static_cast<long>(std::floor((b[0] - bounds[0]) / cell_width));
    const long iy = static_cast<long>(std::floor((b[1] - bounds[2]) / cell_width));
    std::fill(pair_amplitude.begin(), pair_amplitude.end(), 1.0);
    std::size_t active_count = state_count;
    for (long ring = 0; ring <= cell_span && active_count > 0; ++ring) {
      for (long dy = -ring; dy <= ring && active_count > 0; ++dy) {
        for (long dx = -ring; dx <= ring && active_count > 0; ++dx) {
          if (std::max(std::abs(dx), std::abs(dy)) != ring) { continue; }
          const long x = ix + dx;
          const long y = iy + dy;
          if (x < 0 || x >= static_cast<long>(nx) || y < 0 || y >= static_cast<long>(ny)) { continue; }
          for (const auto& point : cell[static_cast<std::size_t>(y) * nx + static_cast<std::size_t>(x)]) {
            const double dx        = point[0] - b[0];
            const double dy        = point[1] - b[1];
            const double distance2 = dx * dx + dy * dy;
            if (distance2 > radius2) { continue; }
            const double distance = std::sqrt(distance2);
            Attenuate(distance, state, pair_amplitude, active_count);
          }
        }
      }
    }
    amplitude[i] = Average(pair_amplitude, state);
  }
  return amplitude;
}

// Compute the configuration Glauber survival amplitude
// S_MC(b) = <S_c(b)>_configs
// [REFERENCE: Glauber, Lectures in Theoretical Physics 1 (1959) 315]
double MGlauber::ConfigAmp(const double b, const MConfigBank* bank1, const MConfigBank* bank2) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument(
        "MGlauber::ConfigAmp: impact parameter must be finite and "
        "nonnegative");
  }
  ValidateBanks(bank1, bank2);
  const std::size_t count1    = bank1 == nullptr ? 1 : bank1->Size();
  const std::size_t count2    = bank2 == nullptr ? 1 : bank2->Size();
  double            amplitude = 0.0;
  for (std::size_t i = 0; i < count1; ++i) {
    for (std::size_t j = 0; j < count2; ++j) { amplitude += ConfigPairAmp(b, 0.0, bank1, bank2, i, j); }
  }
  return amplitude / static_cast<double>(count1 * count2);
}

// Compute the configuration Glauber survival probability
// P_MC(b) = |S_MC(b)|^2
double MGlauber::ConfigProb(const double b, const MConfigBank* bank1, const MConfigBank* bank2) const {
  const double amplitude = ConfigAmp(b, bank1, bank2);
  return amplitude * amplitude;
}

}  // namespace gra::nuclear
