// Correlated nuclear proton and neutron configurations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MConfig.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <random>
#include <stdexcept>
#include <unordered_map>
#include <utility>
#include <vector>

#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra::nuclear {
namespace {

// Draw one isotropic center at a sampled radius
// x = r (sqrt(1 - u^2) cos(phi), sqrt(1 - u^2) sin(phi), u)
std::array<double, 3> DrawCenter(const double r, MRandom &random) {
  const double cosine = random.U(-1.0, 1.0);
  const double sine   = std::sqrt(std::max(0.0, 1.0 - cosine * cosine));
  const double phi    = random.U(0.0, 2.0 * math::PI);
  return {r * sine * std::cos(phi), r * sine * std::sin(phi), r * cosine};
}

// Preserve proton and neutron mean-square radii after independent-particle recentering
// R_s^2 = (1 - 2/A) U_s + sum_i U_i/A^2, U_s = A [R_s^2 - R_m^2/(A - 1)]/(A - 2)
std::array<double, 2> RadiusScales(const MNucleus &nucleus) {
  const double a       = static_cast<double>(nucleus.A());
  const double z       = static_cast<double>(nucleus.Z());
  const double n       = static_cast<double>(nucleus.N());
  const double charge2 = math::pow2(nucleus.ChargeDensity().Rms());
  const double matter2 = math::pow2(nucleus.MatterDensity().Rms());
  if (nucleus.A() < 2) { throw std::invalid_argument("MConfigSampler: recentering requires at least two nucleons"); }
  if (nucleus.A() == 2) {
    if (nucleus.Z() > 0 && std::abs(charge2 - matter2) > 1.0e-10 * matter2) {
      throw std::invalid_argument("MConfigSampler: two recentered nucleons require equal proton and matter radii");
    }
    return {std::sqrt(2.0), std::sqrt(2.0)};
  }
  const std::array<double, 2> radius2 = {charge2, n > 0.0 ? (a * matter2 - z * charge2) / n : charge2};
  std::array<double, 2>       scale   = {1.0, 1.0};
  for (const auto &i : indices(scale)) {
    if ((i == 0 && nucleus.Z() == 0) || (i == 1 && nucleus.N() == 0)) { continue; }
    const double value = a * (radius2[i] - matter2 / (a - 1.0)) / ((a - 2.0) * radius2[i]);
    if (!std::isfinite(value) || !(value > 0.0)) {
      throw std::invalid_argument("MConfigSampler: nuclear radii cannot be preserved by recentering");
    }
    scale[i] = std::sqrt(value);
  }
  return scale;
}

using SeparationCell = std::array<long long, 3>;

// Hash one three-dimensional hard-core cell
struct SeparationCellHash {
  std::size_t operator()(const SeparationCell &cell) const noexcept {
    std::size_t value = 0;
    for (const long long coordinate : cell) {
      const std::size_t part = std::hash<long long>{}(coordinate);
      value ^= part + 0x9e3779b9U + (value << 6U) + (value >> 2U);
    }
    return value;
  }
};

// Index exact neighboring searches for the hard-core nucleon constraint
class SeparationGrid {
 public:
  // Construct one cell list with the hard-core distance as its cell width
  explicit SeparationGrid(const double d_min) : step_(d_min), distance2_(d_min * d_min) {}

  // Compute whether one proposal respects every hard-core pair constraint
  // |r_i - r_j| >= d_min for every accepted pair
  bool Accept(const std::array<double, 3> &center, const std::vector<Nucleon> &nucleon,
              const std::size_t moved = std::numeric_limits<std::size_t>::max()) const {
    if (std::fpclassify(step_) == FP_ZERO) { return true; }
    const SeparationCell cell = Cell(center);
    for (long long dx = -1; dx <= 1; ++dx) {
      for (long long dy = -1; dy <= 1; ++dy) {
        for (long long dz = -1; dz <= 1; ++dz) {
          const SeparationCell neighbor = {cell[0] + dx, cell[1] + dy, cell[2] + dz};
          const auto           found    = bin_.find(neighbor);
          if (found == bin_.cend()) { continue; }
          for (const std::size_t index : found->second) {
            if (index != moved && gra::SquaredNorm(gra::Subtract(center, nucleon[index].x)) < distance2_) {
              return false;
            }
          }
        }
      }
    }
    return true;
  }

  // Insert one accepted nucleon into its exact hard-core cell
  void Add(const std::size_t index, const std::array<double, 3> &center) {
    if (std::fpclassify(step_) != FP_ZERO) { bin_[Cell(center)].push_back(index); }
  }

  // Move one accepted nucleon between exact hard-core cells
  void Move(const std::size_t index, const std::array<double, 3> &previous, const std::array<double, 3> &center) {
    if (std::fpclassify(step_) == FP_ZERO) { return; }
    const SeparationCell old_cell = Cell(previous);
    const SeparationCell new_cell = Cell(center);
    if (old_cell == new_cell) { return; }
    auto found = bin_.find(old_cell);
    if (found == bin_.end()) { throw std::logic_error("SeparationGrid: previous cell is missing"); }
    auto item = std::find(found->second.begin(), found->second.end(), index);
    if (item == found->second.end()) { throw std::logic_error("SeparationGrid: nucleon index is missing"); }
    found->second.erase(item);
    if (found->second.empty()) { bin_.erase(found); }
    bin_[new_cell].push_back(index);
  }

 private:
  double                                                                           step_      = 0.0;
  double                                                                           distance2_ = 0.0;
  std::unordered_map<SeparationCell, std::vector<std::size_t>, SeparationCellHash> bin_;

  // Compute the signed cell containing one nucleon center
  SeparationCell Cell(const std::array<double, 3> &center) const {
    return {static_cast<long long>(std::floor(center[0] / step_)),
            static_cast<long long>(std::floor(center[1] / step_)),
            static_cast<long long>(std::floor(center[2] / step_))};
  }
};

// Recenter one configuration at its nucleon center of mass
// r_i -> r_i - sum_j r_j / A
void Recenter(std::vector<Nucleon> &nucleon) {
  std::array<double, 3> center = {0.0, 0.0, 0.0};
  for (const auto &item : nucleon) { gra::AddScaled(center, item.x, 1.0); }
  gra::Scale(center, 1.0 / static_cast<double>(nucleon.size()));
  for (auto &item : nucleon) { item.x = gra::Subtract(item.x, center); }
}

// Draw one center from the declared density for a selected nucleon type
// r = c_type F_type^-1(u) with separate proton and neutron radius corrections
std::array<double, 3> DrawNucleon(const MNucleus &nucleus, const MConfigSampler &sampler, const NucleonType type,
                                  const std::array<double, 2> &scale, MRandom &random) {
  const double radius = type == NucleonType::Proton ? nucleus.ChargeDensity().Radius(random.U(0.0, 1.0))
                                                    : sampler.NeutronRadius(random.U(0.0, 1.0));
  return DrawCenter(scale[type == NucleonType::Proton ? 0 : 1] * radius, random);
}

// Equilibrate the symmetric hard-core distribution by independence updates
// P({r_i}) proportional to prod_i rho_i(r_i) prod_{i<j} Theta(r_ij - d_min)
void Equilibrate(const MNucleus &nucleus, const ConfigParam &param, const MConfigSampler &sampler,
                 const std::array<double, 2> &scale, std::vector<Nucleon> &nucleon, MRandom &random) {
  if (std::fpclassify(param.d_min) == FP_ZERO || param.sweeps == 0) { return; }
  std::uniform_int_distribution<std::size_t> select(0, nucleon.size() - 1);
  SeparationGrid                             grid(param.d_min);
  for (const auto &i : indices(nucleon)) { grid.Add(i, nucleon[i].x); }
  for (std::size_t sweep = 0; sweep < param.sweeps; ++sweep) {
    for (std::size_t proposal = 0; proposal < nucleon.size(); ++proposal) {
      const std::size_t moved  = select(random.rng);
      const auto        center = DrawNucleon(nucleus, sampler, nucleon[moved].type, scale, random);
      if (grid.Accept(center, nucleon, moved)) {
        grid.Move(moved, nucleon[moved].x, center);
        nucleon[moved].x = center;
      }
    }
  }
}

// Sample one complete proton and neutron configuration
// P({r_i}) proportional to prod_i rho_i(r_i) prod_{i<j} Theta(r_ij - d_min)
std::vector<Nucleon> DrawConfig(const MNucleus &nucleus, const ConfigParam &param, const MConfigSampler &sampler,
                                MRandom &random) {
  const auto               scale = RadiusScales(nucleus);
  std::vector<NucleonType> type(nucleus.A(), NucleonType::Neutron);
  std::fill_n(type.begin(), nucleus.Z(), NucleonType::Proton);
  std::shuffle(type.begin(), type.end(), random.rng);

  std::vector<Nucleon> nucleon;
  nucleon.reserve(nucleus.A());
  SeparationGrid grid(param.d_min);
  for (const NucleonType current_type : type) {
    bool accepted = false;
    for (std::size_t trial = 0; trial < param.max_trials; ++trial) {
      const auto center = DrawNucleon(nucleus, sampler, current_type, scale, random);
      if (grid.Accept(center, nucleon)) {
        nucleon.push_back({center, current_type});
        grid.Add(nucleon.size() - 1, center);
        accepted = true;
        break;
      }
    }
    if (!accepted) { throw std::runtime_error("MConfigBank: failed to sample the requested nucleon correlation"); }
  }
  Equilibrate(nucleus, param, sampler, scale, nucleon, random);
  Recenter(nucleon);
  return nucleon;
}

}  // namespace

// Compute the first two moments of one finite configuration-current bank
// <J> = sum_c J_c / N_c, Var_u[J] = N_c/(N_c - 1) [<|J|^2> - |<J>|^2]
CurrentStat CurrentMoments(const std::vector<std::complex<double>> &current) {
  const auto  moment = statistics::ComplexMoments(current);
  CurrentStat stat{moment.mean, moment.second,
                   current.size() > 1
                       ? moment.variance * static_cast<double>(current.size()) / static_cast<double>(current.size() - 1)
                       : 0.0};
  stat.second = std::norm(stat.mean) + stat.variance;
  if (!std::isfinite(stat.second)) { throw std::runtime_error("CurrentMoments: non-finite current moment"); }
  return stat;
}

// Construct one validated immutable nucleon configuration
MConfig::MConfig(const unsigned int a, const unsigned int z, std::vector<Nucleon> nucleon)
    : a_(a), z_(z), nucleon_(std::move(nucleon)) {
  if (a_ == 0 || z_ > a_ || nucleon_.size() != a_) {
    throw std::invalid_argument("MConfig: inconsistent A, Z or nucleon count");
  }
  proton_.reserve(z_);
  for (const auto &i : indices(nucleon_)) {
    if (nucleon_[i].type == NucleonType::Proton) { proton_.push_back(i); }
  }
  if (proton_.size() != z_) { throw std::invalid_argument("MConfig: inconsistent proton count"); }
  for (const auto &item : nucleon_) {
    if (!gra::AllFinite(item.x)) { throw std::invalid_argument("MConfig: non-finite nucleon center"); }
  }
  transverse_order_.resize(nucleon_.size());
  for (const auto &i : indices(transverse_order_)) { transverse_order_[i] = i; }
  std::sort(transverse_order_.begin(), transverse_order_.end(), [&](const std::size_t first, const std::size_t second) {
    if (nucleon_[first].x[0] < nucleon_[second].x[0]) { return true; }
    if (nucleon_[second].x[0] < nucleon_[first].x[0]) { return false; }
    return first < second;
  });
  transverse_bounds_ = {nucleon_[transverse_order_.front()].x[0], nucleon_[transverse_order_.back()].x[0],
                        nucleon_.front().x[1], nucleon_.front().x[1]};
  for (const auto &item : nucleon_) {
    transverse_bounds_[2] = std::min(transverse_bounds_[2], item.x[1]);
    transverse_bounds_[3] = std::max(transverse_bounds_[3], item.x[1]);
  }
  PreparePhotoGeometry();
}

// Tabulate exact downstream pair geometry for both photon directions
// b_ij^2 = (x_j - x_i)^2 + (y_j - y_i)^2 with sign(z_j - z_i) fixed
void MConfig::PreparePhotoGeometry() {
  for (auto &offset : downstream_offset_) { offset.assign(nucleon_.size() + 1, 0); }
  for (auto &radius2 : downstream_r2_) {
    radius2.clear();
    radius2.reserve(nucleon_.size() * (nucleon_.size() - 1) / 2);
  }
  for (const auto &produced : indices(nucleon_)) {
    const auto &source = nucleon_[produced].x;
    for (const auto &other : indices(nucleon_)) {
      if (other == produced) { continue; }
      const auto  &center = nucleon_[other].x;
      const double dz     = center[2] - source[2];
      if (std::fpclassify(dz) == FP_ZERO) { continue; }
      const double      dx        = center[0] - source[0];
      const double      dy        = center[1] - source[1];
      const std::size_t direction = dz > 0.0 ? 0 : 1;
      downstream_r2_[direction].push_back(dx * dx + dy * dy);
    }
    for (const auto &direction : indices(downstream_offset_)) {
      downstream_offset_[direction][produced + 1] = downstream_r2_[direction].size();
    }
  }
}

// Compute exact downstream transverse separations for one produced nucleon
std::span<const double> MConfig::DownstreamR2(const std::size_t produced, const bool positive_z) const {
  if (produced >= nucleon_.size()) { throw std::out_of_range("MConfig::DownstreamR2: nucleon is out of range"); }
  const std::size_t direction = positive_z ? 0 : 1;
  const std::size_t first     = downstream_offset_[direction][produced];
  const std::size_t count     = downstream_offset_[direction][produced + 1] - first;
  return std::span<const double>(downstream_r2_[direction]).subspan(first, count);
}

// Compute one selected configuration phase sum
// J(q) = sum_i exp(i q.r_i / hbarc)
std::complex<double> MConfig::Current(const double qx, const double qy, const double qz, const bool charge_only) const {
  const std::array<double, 3> q = {qx, qy, qz};
  if (!gra::AllFinite(q)) { throw std::invalid_argument("MConfig::Current: momentum must be finite"); }
  constexpr double     hbarc   = PDG::GeV2fm;
  std::complex<double> current = {0.0, 0.0};
  if (charge_only) {
    for (const std::size_t i : proton_) {
      const double phase = gra::BilinearProduct(q, nucleon_[i].x) / hbarc;
      current += std::polar(1.0, phase);
    }
    return current;
  }
  for (const auto &item : nucleon_) {
    const double phase = gra::BilinearProduct(q, item.x) / hbarc;
    current += std::polar(1.0, phase);
  }
  return current;
}

// Compute the coherent phase sum over proton centers for q in GeV
// J_ch(q) = sum_{i in protons} exp(i q.r_i / hbarc)
std::complex<double> MConfig::ChargeCurrent(const double qx, const double qy, const double qz) const {
  return Current(qx, qy, qz, true);
}

// Compute the coherent phase sum over all nucleon centers for q in GeV
// J_m(q) = sum_{i in nucleons} exp(i q.r_i / hbarc)
std::complex<double> MConfig::MatterCurrent(const double qx, const double qy, const double qz) const {
  return Current(qx, qy, qz, false);
}

// Construct all deterministic radial sampling tables
MConfigSampler::MConfigSampler(const MNucleus &nucleus, ConfigParam param) : nucleus_(nucleus), param_(param) {
  Prepare();
}

// Validate configuration controls and prepare the neutron radial sampler
// rho_n(r) = [A rho_m(r) - Z rho_ch(r)] / N
void MConfigSampler::Prepare() {
  if (nucleus_.ID().lambda != 0) {
    throw std::invalid_argument("MConfigSampler: hypernuclear configurations are not supported");
  }
  if (param_.count == 0 || param_.max_trials == 0 || param_.sweeps > 10000 || !std::isfinite(param_.d_min) ||
      param_.d_min < 0.0 || param_.neutron_nodes < 32 || !std::isfinite(param_.density_rel_tol) ||
      !std::isfinite(param_.negative_norm_tol) || !std::isfinite(param_.norm_tol) || !(param_.density_rel_tol > 0.0) ||
      !(param_.negative_norm_tol > 0.0) || !(param_.norm_tol > 0.0)) {
    throw std::invalid_argument("MConfigSampler: invalid configuration controls");
  }
  static_cast<void>(RadiusScales(nucleus_));

  const unsigned int n = nucleus_.N();
  if (n == 0) { return; }
  const MDensity   &charge    = nucleus_.ChargeDensity();
  const MDensity   &matter    = nucleus_.MatterDensity();
  const double      r_max     = std::max(charge.Param().r_max, matter.Param().r_max);
  const std::size_t intervals = param_.neutron_nodes;
  const double      step      = r_max / static_cast<double>(intervals);
  const double      a         = static_cast<double>(nucleus_.A());
  const double      z         = static_cast<double>(nucleus_.Z());
  const double      inverse_n = 1.0 / static_cast<double>(n);

  neutron_r_.resize(intervals + 1, 0.0);
  std::vector<double> density(intervals + 1, 0.0);
  double              positive_peak = 0.0;
  double              negative_peak = 0.0;
  for (std::size_t i = 0; i <= intervals; ++i) {
    const double r     = static_cast<double>(i) * step;
    const double value = (a * matter.Rho(r) - z * charge.Rho(r)) * inverse_n;
    if (!std::isfinite(value)) { throw std::invalid_argument("MConfigSampler: non-finite derived neutron density"); }
    neutron_r_[i] = r;
    density[i]    = value;
    positive_peak = std::max(positive_peak, value);
    negative_peak = std::max(negative_peak, -value);
  }
  if (!(positive_peak > 0.0)) {
    throw std::invalid_argument("MConfigSampler: derived neutron density is not positive");
  }

  double signed_norm   = 0.0;
  double negative_norm = 0.0;
  neutron_cdf_.resize(intervals + 1, 0.0);
  for (std::size_t i = 1; i <= intervals; ++i) {
    const double previous = 4.0 * math::PI * neutron_r_[i - 1] * neutron_r_[i - 1] * density[i - 1];
    const double current  = 4.0 * math::PI * neutron_r_[i] * neutron_r_[i] * density[i];
    signed_norm += 0.5 * step * (previous + current);
    negative_norm += 0.5 * step * (std::max(0.0, -previous) + std::max(0.0, -current));
    neutron_cdf_[i] = neutron_cdf_[i - 1] + 0.5 * step * (std::max(0.0, previous) + std::max(0.0, current));
  }
  if (negative_peak > param_.density_rel_tol * positive_peak || negative_norm > param_.negative_norm_tol) {
    throw std::invalid_argument(
        "MConfigSampler: charge and matter densities imply a negative neutron "
        "density");
  }
  if (!std::isfinite(signed_norm) || std::abs(signed_norm - 1.0) > param_.norm_tol) {
    throw std::invalid_argument("MConfigSampler: derived neutron density is not normalized");
  }

  const double cdf_norm = neutron_cdf_.back();
  if (!std::isfinite(cdf_norm) || !(cdf_norm > 0.0) || std::abs(cdf_norm - 1.0) > param_.norm_tol) {
    throw std::invalid_argument("MConfigSampler: derived neutron radial sampler is not normalized");
  }
  for (double &value : neutron_cdf_) { value /= cdf_norm; }
  neutron_cdf_.back() = 1.0;
}

// Map one uniform variate to the normalized neutron radial density
// u = 4 pi int_0^r r'^2 rho_n(r') dr'
double MConfigSampler::NeutronRadius(const double u) const {
  if (u < 0.0 || u > 1.0) { throw std::invalid_argument("MConfigSampler::NeutronRadius: variate must be in [0,1]"); }
  if (neutron_cdf_.empty()) { throw std::logic_error("MConfigSampler::NeutronRadius: nucleus has no neutrons"); }
  return math::LinearInterpolateValidatedGrid(neutron_cdf_, neutron_r_, u);
}

// Draw one correlated proton and neutron configuration
// P({r_i}) proportional to prod_i rho_i(r_i) prod_{i<j} Theta(r_ij - d_min)
MConfig MConfigSampler::Draw(MRandom &random) const {
  return MConfig(nucleus_.A(), nucleus_.Z(), DrawConfig(nucleus_, param_, *this, random));
}

// Sample a correlated configuration bank from one prepared density model
MConfigBank::MConfigBank(const MConfigSampler &sampler, MRandom &random)
    : MConfigBank(sampler, random, sampler.Param().count) {}

// Sample a bounded event-local subset from one prepared density model
MConfigBank::MConfigBank(const MConfigSampler &sampler, MRandom &random, const std::size_t count)
    : nucleus_(sampler.Nucleus()), param_(sampler.Param()) {
  if (count == 0 || count > param_.count) { throw std::invalid_argument("MConfigBank: invalid event sample count"); }
  param_.count = count;
  config_.reserve(count);
  for (std::size_t i = 0; i < count; ++i) { config_.push_back(sampler.Draw(random)); }
}

// Compute one fixed configuration with checked indexing
const MConfig &MConfigBank::At(const std::size_t index) const { return config_.at(index); }

// Compute current moments from all fixed configurations
// <J> = sum_c J_c / N_c, Var_u[J] = N_c/(N_c - 1) [<|J|^2> - |<J>|^2]
CurrentStat MConfigBank::Stat(const double qx, const double qy, const double qz, const bool charge_only) const {
  std::vector<std::complex<double>> current(config_.size());
  for (const auto &i : indices(current)) {
    current[i] = charge_only ? config_[i].ChargeCurrent(qx, qy, qz) : config_[i].MatterCurrent(qx, qy, qz);
  }
  return CurrentMoments(current);
}

// Compute charge-current ensemble moments for q in GeV
// <J_ch> = sum_c J_ch,c / N_c, Var_u[J_ch] = N_c/(N_c - 1) [<|J_ch|^2> - |<J_ch>|^2]
CurrentStat MConfigBank::ChargeStat(const double qx, const double qy, const double qz) const {
  return Stat(qx, qy, qz, true);
}

// Compute matter-current ensemble moments for q in GeV
// <J_m> = sum_c J_m,c / N_c, Var_u[J_m] = N_c/(N_c - 1) [<|J_m|^2> - |<J_m>|^2]
CurrentStat MConfigBank::MatterStat(const double qx, const double qy, const double qz) const {
  return Stat(qx, qy, qz, false);
}

}  // namespace gra::nuclear
