// Coherent and incoherent photonuclear target transitions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MPhoto.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Nuclear/MCache.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tech/MJsonZip.h"
#include "json.hpp"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {
namespace {

// Compute the current persistent photonuclear geometry representation version
constexpr std::size_t PhotoCacheVersion() { return 2; }

// Resolve density gradients and the largest Fourier phase on a Gauss rule
unsigned int PhotoNodes(unsigned int minimum, double momentum, double length, double scale) {
  const double count = std::ceil(0.5 * scale * std::abs(momentum) * length / PDG::GeV2fm);
  if (!std::isfinite(count) || count > 4096.0) {
    throw AmplitudeFailure("MPhoto: momentum exceeds the resolved quadrature");
  }
  return std::max(minimum, static_cast<unsigned int>(count));
}

// Compute the checked array index of one photon propagation direction
std::size_t DirectionIndex(const PhotonDirection direction) {
  switch (direction) {
    case PhotonDirection::PositiveZ:
      return 0;
    case PhotonDirection::NegativeZ:
      return 1;
  }
  throw std::invalid_argument("MPhoto: invalid photon direction");
}

// Validate one channel and energy dependent elementary profile
// Gamma_N(0) = sigma_eff / (4 pi B hbarc^2), Gamma_N(0) <= 2 / (1 + eta^2)
void ValidateProfile(const PhotoProfile& profile) {
  if (!std::isfinite(profile.sigma_eff) || !std::isfinite(profile.slope) || !std::isfinite(profile.eta) ||
      !std::isfinite(profile.omega) || profile.sigma_eff < 0.0 || !(profile.slope > 0.0) || profile.omega < 0.0) {
    throw AmplitudeFailure("MPhoto: invalid elementary profile");
  }
  constexpr double hbarc  = PDG::GeV2fm;
  const double     gamma0 = 0.1 * profile.sigma_eff / (4.0 * math::PI * profile.slope * hbarc * hbarc);
  if (gamma0 > 2.0 / (1.0 + profile.eta * profile.eta)) {
    throw AmplitudeFailure("MPhoto: elementary profile violates impact unitarity");
  }
}

}  // namespace

// Compute the incident photon direction for one target beam leg
PhotonDirection TargetPhotonDirection(const int target_leg) {
  if (target_leg == 1) { return PhotonDirection::NegativeZ; }
  if (target_leg == 2) { return PhotonDirection::PositiveZ; }
  throw std::out_of_range("TargetPhotonDirection: target leg must be one or two");
}

// Construct one immutable photonuclear target model
MPhoto::MPhoto(const MNucleus& target, PhotoParam param, const PhotoModel model)
    : MPhoto(std::make_shared<const MNucleus>(target), param, model) {}

// Construct one target transition sharing an immutable nuclear model
MPhoto::MPhoto(std::shared_ptr<const MNucleus> target, PhotoParam param, const PhotoModel model)
    : target_(std::move(target)), param_(param), model_(model) {
  if (target_ == nullptr) { throw std::invalid_argument("MPhoto: null nuclear model"); }
  Prepare();
}

// Validate controls and construct the radial integration rule
void MPhoto::Prepare() {
  sigma_scale_ = CrossSectionRule(param_.fluctuation);
  if (param_.b_nodes < 32 || param_.b_nodes > 4096) {
    throw std::invalid_argument("MPhoto: b_nodes must be in [32,4096]");
  }
  if (param_.z_nodes < 32 || param_.z_nodes > 4096) {
    throw std::invalid_argument("MPhoto: z_nodes must be in [32,4096]");
  }
  if (!std::isfinite(param_.phase_scale) || param_.phase_scale < 1.0) {
    throw std::invalid_argument("MPhoto: phase_scale must be finite and at least one");
  }
  if (model_ == PhotoModel::Impulse) { return; }
  {
    if (!std::isfinite(param_.table.qt_max) || !(param_.table.qt_max > 0.0) || !std::isfinite(param_.table.qz_max) ||
        !(param_.table.qz_max > 0.0) || param_.table.qt_nodes < 4 || param_.table.qt_nodes > 4097 ||
        param_.table.qz_nodes < 4 || param_.table.qz_nodes > 4097 || param_.table.series_terms < 4 ||
        param_.table.series_terms > 256 || !std::isfinite(param_.table.series_abs_tol) ||
        !(param_.table.series_abs_tol > 0.0) || !(param_.table.series_abs_tol < 1.0)) {
      throw std::invalid_argument("MPhoto: invalid optical response controls");
    }
    constexpr std::size_t response_limit = 50000000;
    const std::size_t     plane          = static_cast<std::size_t>(param_.table.qt_nodes) * param_.table.qz_nodes;
    if (plane > response_limit / param_.table.series_terms) {
      throw std::invalid_argument("MPhoto: optical response table is too large");
    }
  }
  const auto& density = target_->MatterDensity().Param();
  const auto& charge = target_->ChargeDensity().Param();
  r_max_ = std::max(density.r_max, charge.r_max);
  const auto nodes = std::max(density.nodes, charge.nodes);
  try {
    param_.b_nodes = PhotoNodes(std::max(param_.b_nodes, nodes), param_.table.qt_max, r_max_, param_.phase_scale);
    param_.z_nodes = PhotoNodes(std::max(param_.z_nodes, nodes), param_.table.qz_max, 2.0 * r_max_, param_.phase_scale);
  } catch (const AmplitudeFailure& error) { throw std::invalid_argument(error.what()); }
  auto rule   = math::GaussLegendreRule(param_.b_nodes, 0.0, r_max_);
  b_node_     = std::move(rule.first);
  b_weight_   = std::move(rule.second);
  auto z_rule = math::GaussLegendreRule(param_.z_nodes, -r_max_, r_max_);
  z_node_     = std::move(z_rule.first);
  z_weight_   = std::move(z_rule.second);

  if (model_ == PhotoModel::LTA) { shadow_model_ = std::make_unique<MShadow>(param_.shadow); }

  const std::string key      = CacheKey();
  const std::string filename = NuclearCacheFilename("PHOTOPROD_" + std::to_string(target_->ID().pdg), key);
  MCacheLock        lock(filename + ".lock");
  const bool        loaded = ReadCache(filename, key);
  if (loaded) {
    std::cout << "Loaded photonuclear geometry cache: " << filename << std::endl;
  } else {
    PrepareOpticalKernel();
    PrepareResponse();
    if (lock.Acquired()) {
      try {
        WriteCache(filename, key);
        std::cout << "Saved photonuclear geometry cache: " << filename << std::endl;
      } catch (const std::exception& error) {
        std::cerr << "WARNING: MPhoto cache not saved: " << error.what() << std::endl;
      }
    }
  }
  depth_max_ = *std::max_element(optical_depth_.cbegin(), optical_depth_.cend());
  thickness_.resize(b_node_.size());
  for (const auto& ib : indices(b_node_)) { thickness_[ib] = target_->MatterDensity().Thick(b_node_[ib]); }
  if (model_ == PhotoModel::Glauber) {
    longitudinal_.assign(b_node_.size() * param_.table.qz_nodes, 0.0);
    charge_longitudinal_.assign(longitudinal_.size(), 0.0);
    ParallelFor(param_.table.qz_nodes, [&](const std::size_t iq) {
      const double qz = iq * param_.table.qz_max / (param_.table.qz_nodes - 1);
      for (const auto& iz : indices(z_node_)) {
        const double phase = std::cos(qz * z_node_[iz] / PDG::GeV2fm);
        for (const auto& ib : indices(b_node_)) {
          longitudinal_[ib * param_.table.qz_nodes + iq] += density_weight_[ib * z_node_.size() + iz] * phase;
          charge_longitudinal_[ib * param_.table.qz_nodes + iq] += charge_weight_[ib * z_node_.size() + iz] * phase;
        }
      }
    });
  }
}

// Tabulate the common radial and matter-density integration weights
void MPhoto::PrepareDensityKernel() {
  const std::size_t b_size = b_node_.size();
  const std::size_t z_size = z_node_.size();
  const std::size_t size   = b_size * z_size;
  radial_weight_.resize(b_size, 0.0);
  density_weight_.resize(size, 0.0);
  charge_weight_.resize(size, 0.0);
  ParallelFor(
      b_size,
      [&](const std::size_t ib) {
        radial_weight_[ib]      = 2.0 * math::PI * b_weight_[ib] * b_node_[ib];
        const std::size_t first = ib * z_size;
        for (const auto& iz : indices(z_node_)) {
          const std::size_t index = first + iz;
          const double radius = std::hypot(b_node_[ib], z_node_[iz]);
          density_weight_[index] = z_weight_[iz] * target_->MatterDensity().Rho(radius);
          charge_weight_[index] = z_weight_[iz] * target_->ChargeDensity().Rho(radius);
        }
      },
      "Photonuclear density kernel");
}

// Integrate outgoing thickness directly on the optical quadrature
// D(b,z) = int_z^infinity rho_m(sqrt(b^2+s^2)) ds
void MPhoto::PrepareOpticalKernel() {
  PrepareDensityKernel();
  optical_depth_.assign(density_weight_.size(), 0.0);
  const auto  rule    = math::GaussLegendreRule(param_.z_nodes, 0.0, 1.0);
  const auto& density = target_->MatterDensity();
  ParallelFor(
      b_node_.size(),
      [&](std::size_t ib) {
        const double b = b_node_[ib];
        if (b >= r_max_) { return; }
        const double limit = std::sqrt(r_max_ * r_max_ - b * b);
        if (model_ == PhotoModel::LTA) {
          const double thickness = density.Thick(b);
          for (const auto& iz : indices(z_node_)) { optical_depth_[ib * z_node_.size() + iz] = thickness; }
          return;
        }
        for (const auto& iz : indices(z_node_)) {
          const double lower  = std::clamp(z_node_[iz], -limit, limit);
          const double length = limit - lower;
          double       depth  = 0.0;
          for (const auto& k : indices(rule.first)) {
            depth += length * rule.second[k] * density.Rho(std::hypot(b, lower + length * rule.first[k]));
          }
          optical_depth_[ib * z_node_.size() + iz] = depth;
        }
      },
      "Photonuclear optical kernel");
}

// Tabulate momentum transforms of the optical-depth moment basis
// R_n(q) = int d^2b dz rho_{m,p}(r) D(b,z)^n J_0(q_T b/hbarc) exp(i q_z z/hbarc)
void MPhoto::PrepareResponse() {
  const std::size_t qt_size = param_.table.qt_nodes, qz_size = param_.table.qz_nodes;
  const std::size_t terms = param_.table.series_terms;
  MMatrix<double>   bessel(b_node_.size(), qt_size);
  ParallelFor(qt_size, [&](const std::size_t iq) {
    const double qt = iq * param_.table.qt_max / static_cast<double>(qt_size - 1);
    for (const auto& ib : indices(b_node_)) {
      bessel(ib, iq) = radial_weight_[ib] * std::cyl_bessel_j(0.0, qt * b_node_[ib] / PDG::GeV2fm);
    }
  });

  // TensorProductProjection conjugates the basis, giving exp(+i qz z) in the response
  MMatrix<std::complex<double>> phase(z_node_.size(), qz_size);
  ParallelFor(qz_size, [&](const std::size_t iq) {
    const double qz = iq * param_.table.qz_max / static_cast<double>(qz_size - 1);
    for (const auto& iz : indices(z_node_)) { phase(iz, iq) = std::polar(1.0, -qz * z_node_[iz] / PDG::GeV2fm); }
  });

  const auto project = [&](const std::vector<double>& weight, std::vector<std::complex<double>>& table) {
    table.resize(terms * qt_size * qz_size);
    ParallelFor(terms, [&](const std::size_t order) {
      std::vector<double> density(weight);
      for (const auto& index : indices(density)) { density[index] *= std::pow(optical_depth_[index], order); }
      const auto response = TensorProductProjection(density, bessel, phase);
      std::copy(response.begin(), response.end(), table.begin() + order * qt_size * qz_size);
    }, "Photonuclear optical response table");
  };
  project(density_weight_, response_);
  project(charge_weight_, charge_response_);
}

// Construct the exact nuclear-geometry cache key
std::string MPhoto::CacheKey() const {
  const auto&        density = target_->MatterDensity().Param();
  std::ostringstream stream;
  stream << PhotoCacheVersion() << ';' << static_cast<int>(model_) << ';' << target_->ID().pdg << ';' << std::hexfloat
         << r_max_ << ';' << param_.phase_scale << ';' << param_.b_nodes << ';' << param_.z_nodes << ';'
         << param_.table.qt_max << ';' << param_.table.qz_max << ';' << param_.table.qt_nodes << ';'
         << param_.table.qz_nodes << ';' << param_.table.series_terms << ';' << param_.table.series_abs_tol << ';'
         << density.radius << ';' << density.skin << ';' << density.r_max << ';' << density.nodes << ';'
         << density.cdf_nodes << ';' << density.form_q_max << ';' << density.form_nodes << ';' << density.form_abs_tol
         << ';';
  const auto& charge = target_->ChargeDensity().Param();
  stream << charge.radius << ';' << charge.skin << ';' << charge.r_max << ';' << charge.nodes << ';'
         << charge.cdf_nodes << ';' << charge.form_q_max << ';' << charge.form_nodes << ';' << charge.form_abs_tol;
  return stream.str();
}

// Load one exact kinematics-independent optical kernel
bool MPhoto::ReadCache(const std::string& filename, const std::string& key) {
  try {
    nlohmann::json cache;
    if (!ReadNuclearCache(filename, PhotoCacheVersion(), "photoprod_geometry", key, cache)) { return false; }
    const std::size_t optical_size = static_cast<std::size_t>(param_.b_nodes) * param_.z_nodes;
    const std::size_t response_size =
        static_cast<std::size_t>(param_.table.series_terms) * param_.table.qt_nodes * param_.table.qz_nodes;
    return MJsonZip::DecompressVector(cache.at("radial_weight"), radial_weight_, param_.b_nodes) &&
           MJsonZip::DecompressVector(cache.at("density_weight"), density_weight_, optical_size) &&
           MJsonZip::DecompressVector(cache.at("charge_weight"), charge_weight_, optical_size) &&
           MJsonZip::DecompressComplexVector(cache.at("charge_response"), charge_response_, response_size) &&
           MJsonZip::DecompressVector(cache.at("optical_depth"), optical_depth_, optical_size) &&
           MJsonZip::DecompressComplexVector(cache.at("response"), response_, response_size);
  } catch (const std::exception&) { return false; }
}

// Save one exact kinematics-independent optical kernel
void MPhoto::WriteCache(const std::string& filename, const std::string& key) const {
  auto cache = NuclearCache(PhotoCacheVersion(), "photoprod_geometry", key);
  cache.update({{"radial_weight", MJsonZip::CompressVector(radial_weight_)},
                {"density_weight", MJsonZip::CompressVector(density_weight_)},
                {"charge_weight", MJsonZip::CompressVector(charge_weight_)},
                {"charge_response", MJsonZip::CompressComplexVector(charge_response_)},
                {"optical_depth", MJsonZip::CompressVector(optical_depth_)},
                {"response", MJsonZip::CompressComplexVector(response_)}});
  PublishNuclearCache(filename, cache);
}

// Compute configured or event-specific cross-section eigenvalues
// <s> = 1 and Var[s] = omega
// The Gamma ansatz fixes two moments, with an assumed small-sigma tail and higher moments
const std::vector<double>& MPhoto::CrossSectionScales(const PhotoProfile& profile, std::vector<double>& scratch) const {
  const double scale = std::max({1.0, std::abs(profile.omega), std::abs(param_.fluctuation.omega)});
  if (std::abs(profile.omega - param_.fluctuation.omega) <= 16.0 * std::numeric_limits<double>::epsilon() * scale) {
    return sigma_scale_;
  }
  FluctuationParam fluctuation = param_.fluctuation;
  fluctuation.omega            = profile.omega;
  try {
    scratch = CrossSectionRule(fluctuation);
  } catch (const std::exception& error) { throw AmplitudeFailure(error.what()); }
  return scratch;
}

// Compute the coherent optical mean from cached outgoing-depth Fourier moments
// For a spherical density, J_-(qt,qz) = J_+(qt,-qz)
CurrentStat MPhoto::OpticalStat(const PhotoProfile& profile, double qt, double qz, PhotonDirection direction) const {
  qz = DirectionIndex(direction) == 0 ? qz : -qz;
  if (model_ != PhotoModel::Impulse) { ValidateProfile(profile); }
  if (!std::isfinite(qt) || !std::isfinite(qz) || qt < 0.0) { throw AmplitudeFailure("MPhoto: invalid momentum"); }
  if (model_ != PhotoModel::Glauber || !(profile.sigma_eff > 0.0) || qt > param_.table.qt_max ||
      std::abs(qz) > param_.table.qz_max) {
    return OpticalStatDirect(profile, qt, qz);
  }
  CurrentStat                       stat;
  const double                      a     = static_cast<double>(target_->A());
  const double                      sigma = 0.1 * profile.sigma_eff;
  std::vector<double>               sigma_scratch;
  const std::vector<double>&        sigma_scale = CrossSectionScales(profile, sigma_scratch);
  std::vector<std::complex<double>> opacity(sigma_scale.size());
  double                            max_opacity = 0.0;
  for (const auto& i : indices(opacity)) {
    opacity[i]  = 0.5 * sigma * sigma_scale[i] * a * std::complex<double>(1.0, -profile.eta);
    max_opacity = std::max(max_opacity, std::abs(opacity[i]));
  }
  std::size_t terms = 4;
  for (; terms <= param_.table.series_terms; ++terms) {
    const double mean_tail = a * math::ExponentialTailBound(max_opacity * depth_max_, terms);
    if (mean_tail <= param_.table.series_abs_tol) { break; }
  }
  if (terms > param_.table.series_terms) { return OpticalStatDirect(profile, qt, qz); }

  const auto moments = Moments(qt, qz, terms, profile.isospin);

  std::vector<std::complex<double>> mean_coeff(sigma_scale.begin(), sigma_scale.end());
  const double                      inverse = 1.0 / static_cast<double>(opacity.size());
  for (std::size_t order = 0; order < terms; ++order) {
    std::complex<double> mean_weight = 0.0;
    for (const auto& coefficient : mean_coeff) { mean_weight += coefficient; }
    mean_weight *= inverse;
    stat.mean += a * mean_weight * moments[order];

    const double divisor = static_cast<double>(order + 1);
    for (const auto& i : indices(mean_coeff)) { mean_coeff[i] *= -opacity[i] / divisor; }
  }
  stat.second = std::norm(stat.mean) + stat.variance;
  return stat;
}

// Compute the source density for neutral isoscalar or isovector production
// J_0 = J_p + J_n, J_1 = J_p - J_n with A rho_m = Z rho_p + N rho_n
double MPhoto::Source(double matter, double charge, PhotoIsospin isospin) const {
  return isospin == PhotoIsospin::Isoscalar ? matter : 2.0 * target_->Z() / target_->A() * charge - matter;
}

// Interpolate the production density moments on the common optical-depth basis
std::vector<std::complex<double>> MPhoto::Moments(double qt, double qz, std::size_t terms, PhotoIsospin isospin) const {
  const auto t = math::UniformCubicWeights(param_.table.qt_nodes, 0.0, param_.table.qt_max, qt);
  const auto z = math::UniformCubicWeights(param_.table.qz_nodes, 0.0, param_.table.qz_max, std::abs(qz));
  std::vector<std::complex<double>> moments(terms, 0.0);
  for (const auto& order : indices(moments)) {
    for (const auto& it : indices(t.index)) {
      for (const auto& iz : indices(z.index)) {
        const auto index = (order * param_.table.qt_nodes + t.index[it]) * param_.table.qz_nodes + z.index[iz];
        const auto response = isospin == PhotoIsospin::Isoscalar ? response_[index]
            : 2.0 * target_->Z() / target_->A() * charge_response_[index] - response_[index];
        moments[order] += t.weight[it] * z.weight[iz] * response;
      }
    }
    if (qz < 0.0) { moments[order] = std::conj(moments[order]); }
  }
  // The zero-opacity term uses the independently validated spherical transform
  const double q = std::hypot(qt, qz);
  moments[0] = Source(target_->MatterDensity().Form(q), target_->ChargeDensity().Form(q), isospin);
  return moments;
}

// Compute direct optical current moments with a normalized zero-opacity term
// Omega_s = A sigma_eff s (1 - i eta) / 2
// <J_I> = A int rho_I exp(i q.r/hbarc) <s exp[-Omega_s D]>_s d^3r
// rho_0 = rho_m, rho_1 = 2 Z rho_p / A - rho_m
// [REFERENCE: Frankfurt et al., Phys. Lett. B752 (2016) 51, Eq. (10)]
CurrentStat MPhoto::OpticalStatDirect(const PhotoProfile& profile, const double qt, const double qz) const {
  const double a = static_cast<double>(target_->A());
  CurrentStat  stat;
  if (model_ == PhotoModel::LTA) { return LeadingStat(profile, qt, qz); }
  if (model_ == PhotoModel::Impulse || std::fpclassify(profile.sigma_eff) == FP_ZERO) {
    const double q = std::hypot(qt, qz), z = target_->Z(), n = target_->N();
    const double matter = target_->MatterDensity().Form(q), proton = target_->ChargeDensity().Form(q);
    const double neutron = n > 0.0 ? (a * matter - z * proton) / n : 0.0;
    stat.mean = a * Source(matter, proton, profile.isospin);
    stat.variance = std::max(0.0, a - z * proton * proton - n * neutron * neutron);
    stat.second       = std::norm(stat.mean) + stat.variance;
    return stat;
  }

  constexpr double                  hbarc = PDG::GeV2fm;
  const double                      sigma = 0.1 * profile.sigma_eff;
  std::vector<double>               sigma_scratch;
  const std::vector<double>&        sigma_scale = CrossSectionScales(profile, sigma_scratch);
  std::vector<std::complex<double>> opacity(sigma_scale.size(), 0.0);
  for (const auto& i : indices(sigma_scale)) {
    opacity[i] = 0.5 * sigma * sigma_scale[i] * a * std::complex<double>(1.0, -profile.eta);
  }
  // Refine only transfers whose phases exceed the prepared geometry resolution
  const auto nb     = PhotoNodes(param_.b_nodes, qt, r_max_, param_.phase_scale);
  const auto nz     = PhotoNodes(param_.z_nodes, qz, 2.0 * r_max_, param_.phase_scale);
  const bool refine = nb > param_.b_nodes || nz > param_.z_nodes;
  const auto b_rule =
      refine ? math::GaussLegendreRule(nb, 0.0, r_max_) : std::pair<std::vector<double>, std::vector<double>>{};
  const auto z_rule =
      refine ? math::GaussLegendreRule(nz, -r_max_, r_max_) : std::pair<std::vector<double>, std::vector<double>>{};
  const auto                        depth_rule = refine ? math::GaussLegendreRule(param_.z_nodes, 0.0, 1.0)
                                                        : std::pair<std::vector<double>, std::vector<double>>{};
  const auto&                       b_node     = refine ? b_rule.first : b_node_;
  const auto&                       z_node     = refine ? z_rule.first : z_node_;
  const auto&                       density    = target_->MatterDensity();
  std::vector<std::complex<double>> z_phase(z_node.size(), 0.0);
  for (const auto& iz : indices(z_node)) { z_phase[iz] = std::polar(1.0, qz * z_node[iz] / hbarc); }
  stat.mean = a * Source(density.Form(std::hypot(qt, qz)), target_->ChargeDensity().Form(std::hypot(qt, qz)),
                         profile.isospin);
  for (const auto& ib : indices(b_node)) {
    const double         b            = b_node[ib];
    const double         bessel       = std::cyl_bessel_j(0.0, qt * b / hbarc);
    std::complex<double> longitudinal = 0.0;
    const std::size_t    first        = ib * z_node.size();
    const double         area         = refine ? 2.0 * math::PI * b * b_rule.second[ib] : radial_weight_[ib];
    const double         limit        = std::sqrt(std::max(0.0, r_max_ * r_max_ - b * b));
    for (const auto& iz : indices(z_node)) {
      const std::size_t index = first + iz;
      double            depth = 0.0;
      if (refine) {
        const double lower = std::clamp(z_node[iz], -limit, limit), length = limit - lower;
        for (const auto& k : indices(depth_rule.first)) {
          depth += length * depth_rule.second[k] * density.Rho(std::hypot(b, lower + length * depth_rule.first[k]));
        }
      } else {
        depth = optical_depth_[index];
      }
      const double radius = std::hypot(b, z_node[iz]);
      const double matter = refine ? z_rule.second[iz] * density.Rho(radius) : density_weight_[index];
      const double charge = refine ? z_rule.second[iz] * target_->ChargeDensity().Rho(radius) : charge_weight_[index];
      const double weight = Source(matter, charge, profile.isospin);
      std::complex<double> attenuation = 0.0;
      for (const auto& eigenstate : indices(opacity)) {
        attenuation += sigma_scale[eigenstate] * std::exp(-opacity[eigenstate] * depth);
      }
      attenuation /= static_cast<double>(sigma_scale.size());
      longitudinal += weight * z_phase[iz] * (attenuation - 1.0);
    }
    stat.mean += a * area * bessel * longitudinal;
  }
  stat.second = std::norm(stat.mean) + stat.variance;
  if (!std::isfinite(stat.mean.real()) || !std::isfinite(stat.mean.imag()) || !std::isfinite(stat.second) ||
      !std::isfinite(stat.variance)) {
    throw AmplitudeFailure("MPhoto::OpticalStat: non-finite current moment");
  }
  return stat;
}

// Compute long-coherence optical fluctuations with elastic rescattering and finite-A subtraction
// The longitudinal form factor does not supply finite-coherence rescattering correlations
// [REFERENCE: Guzey, Strikman and Zhalov, Eur. Phys. J. C 74 (2014) 2942, Eqs. (13)-(14)]
double MPhoto::GlauberVariance(const PhotoProfile& profile, double qt, double qz) const {
  std::vector<double> scratch;
  const auto&         scale = CrossSectionScales(profile, scratch);
  const double        a = target_->A(), sigma = 0.1 * profile.sigma_eff;
  const auto          c = 0.5 * a * sigma * std::complex<double>(1.0, -profile.eta);
  // The eigenstate profile has transverse slope B_s = s B
  const double elastic   = a * sigma * sigma * (1.0 + profile.eta * profile.eta) /
                           (8.0 * math::PI * profile.slope * PDG::GeV2fm * PDG::GeV2fm);
  const bool   tabulated = qt <= param_.table.qt_max && std::abs(qz) <= param_.table.qz_max;
  const auto b = tabulated
                     ? std::pair<std::vector<double>, std::vector<double>>{}
                     : math::GaussLegendreRule(PhotoNodes(param_.b_nodes, qt, r_max_, param_.phase_scale), 0.0, r_max_);
  const auto z = tabulated ? std::pair<std::vector<double>, std::vector<double>>{}
                           : math::GaussLegendreRule(PhotoNodes(param_.z_nodes, qz, 2.0 * r_max_, param_.phase_scale),
                                                     -r_max_, r_max_);
  const auto stencil          = math::UniformCubicWeights(param_.table.qz_nodes, 0.0, param_.table.qz_max,
                                                          std::min(std::abs(qz), param_.table.qz_max));
  const auto&          radius = tabulated ? b_node_ : b.first;
  double               local  = 1.0;
  std::complex<double> mean = target_->MatterDensity().Form(std::hypot(qt, qz));
  std::complex<double> proton = target_->ChargeDensity().Form(std::hypot(qt, qz));
  for (const auto& ib : indices(radius)) {
    const double t     = tabulated ? thickness_[ib] : target_->MatterDensity().Thick(radius[ib]);
    const double area  = tabulated ? radial_weight_[ib] : 2.0 * math::PI * radius[ib] * b.second[ib];
    double       phase = 0.0, charge_phase = 0.0, second = 0.0;
    if (tabulated) {
      for (const auto& i : indices(stencil.index)) {
        phase += stencil.weight[i] * longitudinal_[ib * param_.table.qz_nodes + stencil.index[i]];
        charge_phase += stencil.weight[i] * charge_longitudinal_[ib * param_.table.qz_nodes + stencil.index[i]];
      }
    } else {
      for (const auto& iz : indices(z.first)) {
        const double r = std::hypot(radius[ib], z.first[iz]);
        const double weight = z.second[iz] * std::cos(qz * z.first[iz] / PDG::GeV2fm);
        phase += weight * target_->MatterDensity().Rho(r);
        charge_phase += weight * target_->ChargeDensity().Rho(r);
      }
    }
    std::complex<double> first = 0.0;
    for (const double s : scale) {
      first += s * std::exp(-c * s * t);
      for (const double u : scale) {
        second += s * u * std::exp((-c * s - std::conj(c) * u + elastic * s * u / (s + u)) * t).real();
      }
    }
    first /= scale.size();
    second /= scale.size() * scale.size();
    const auto transform = area * math::BesselJ0(qt * radius[ib] / PDG::GeV2fm) * (first - 1.0);
    mean += transform * phase;
    proton += transform * charge_phase;
    local += area * t * (second - 1.0);
  }
  const double z_count = target_->Z(), n_count = target_->N();
  const auto neutron = n_count > 0.0 ? (a * mean - z_count * proton) / n_count : std::complex<double>{};
  return std::max(0.0, a * local - z_count * std::norm(proton) - n_count * std::norm(neutron));
}

// Compute the LTA density response with the finite-A coherent subtraction
// g(T) = (1-r) T + 2r [1-exp(-sigma_3 A T/2)]/(sigma_3 A)
// w(T) = dg/dT = 1-r+r exp(-sigma_3 A T/2)
// Var[J] = A int rho_m w^2 - Z |int rho_p w exp(i q.r)|^2 - N |int rho_n w exp(i q.r)|^2
// This is the linear density response, with sigma_3^in approximated by sigma_3
// [REFERENCE: Guzey, Strikman and Zhalov, Eur. Phys. J. C 74 (2014) 2942, Eqs. (9), (15)-(16)]
CurrentStat MPhoto::LeadingStat(const PhotoProfile& profile, const double qt, const double qz) const {
  const auto           xs = shadow_model_->CrossSections(profile.x, profile.scale2);
  const double         a = target_->A(), opacity = 0.05 * xs.sigma3 * a, r = xs.ratio;
  CurrentStat          stat;
  std::complex<double> response = 0.0, proton = 0.0;
  double               local    = 0.0;
  std::size_t          terms    = 4;
  for (; terms <= param_.table.series_terms; ++terms) {
    if (a * math::ExponentialTailBound(opacity * depth_max_, terms) <= param_.table.series_abs_tol) { break; }
  }
  if (qt <= param_.table.qt_max && std::abs(qz) <= param_.table.qz_max && terms <= param_.table.series_terms) {
    const auto moments = Moments(qt, qz, terms);
    const auto vector = Moments(qt, qz, terms, PhotoIsospin::Isovector);
    double     coefficient = 1.0;
    for (const auto& n : indices(moments)) {
      const double direct = n == 0 ? 1.0 - r : 0.0;
      stat.mean += a * (direct + r * coefficient / (n + 1.0)) * moments[n];
      response += (direct + r * coefficient) * moments[n];
      if (target_->Z() > 0) { proton += (direct + r * coefficient) * (moments[n] + vector[n]) * a / (2.0 * target_->Z()); }
      coefficient *= -opacity / (n + 1.0);
    }
    // The second moment has no momentum phase, so a single radial sum avoids a doubled-opacity series
    local = 1.0;
    for (const auto& ib : indices(thickness_)) {
      const double delta = r * std::expm1(-opacity * thickness_[ib]);
      local += radial_weight_[ib] * thickness_[ib] * delta * (2.0 + delta);
    }
  } else {
    // Direct evaluation is reserved for transfers outside the interpolation table
    const auto& density = target_->MatterDensity();
    const auto  b = math::GaussLegendreRule(PhotoNodes(param_.b_nodes, qt, r_max_, param_.phase_scale), 0.0, r_max_);
    const auto  z =
        math::GaussLegendreRule(PhotoNodes(param_.z_nodes, qz, 2.0 * r_max_, param_.phase_scale), -r_max_, r_max_);
    const double form = density.Form(std::hypot(qt, qz));
    stat.mean         = a * form;
    response          = form;
    proton            = target_->ChargeDensity().Form(std::hypot(qt, qz));
    local             = 1.0;
    std::vector<std::complex<double>> phase(z.first.size());
    for (const auto& iz : indices(phase)) { phase[iz] = std::polar(1.0, qz * z.first[iz] / PDG::GeV2fm); }
    for (const auto& ib : indices(b.first)) {
      const double         thickness    = density.Thick(b.first[ib]);
      const double         x            = opacity * thickness;
      const double         coherent     = 1.0 - r + r * (x > 1.0e-12 ? -std::expm1(-x) / x : 1.0);
      const double         weight       = 1.0 - r + r * std::exp(-x);
      std::complex<double> longitudinal = 0.0, charge = 0.0;
      for (const auto& iz : indices(z.first)) {
        const double radius = std::hypot(b.first[ib], z.first[iz]);
        longitudinal += z.second[iz] * density.Rho(radius) * phase[iz];
        charge += z.second[iz] * target_->ChargeDensity().Rho(radius) * phase[iz];
      }
      const double area      = 2.0 * math::PI * b.first[ib] * b.second[ib];
      const auto   transform = area * std::cyl_bessel_j(0.0, qt * b.first[ib] / PDG::GeV2fm) * longitudinal;
      stat.mean += a * (coherent - 1.0) * transform;
      response += (weight - 1.0) * transform;
      proton += (weight - 1.0) * area * std::cyl_bessel_j(0.0, qt * b.first[ib] / PDG::GeV2fm) * charge;
      local += area * thickness * (weight - 1.0) * (weight + 1.0);
    }
  }
  const double z_count = target_->Z(), n_count = target_->N();
  const auto neutron = n_count > 0.0 ? (a * response - z_count * proton) / n_count : std::complex<double>{};
  stat.variance = std::max(0.0, a * local - z_count * std::norm(proton) - n_count * std::norm(neutron));
  stat.second   = std::norm(stat.mean) + stat.variance;
  return stat;
}

// Interpolate the full nuclear thickness for the LTA density response
// Both propagation directions cross the same spherical target thickness
double MPhoto::Thickness(double b) const {
  if (b >= r_max_) { return 0.0; }
  if (b > b_node_.back()) { return target_->MatterDensity().Thick(b); }
  const auto            upper = std::upper_bound(b_node_.cbegin(), b_node_.cend(), b);
  const auto            count = static_cast<std::size_t>(std::distance(b_node_.cbegin(), upper));
  const auto            first = std::min(count > 1 ? count - 2 : 0, b_node_.size() - 4);
  std::array<double, 4> node, value;
  for (const auto& i : indices(node)) {
    node[i]  = b_node_[first + i];
    value[i] = thickness_[first + i];
  }
  return std::max(0.0, math::CubicLagrangeInterpolate(node, value, b));
}

// Factor the elementary profile into exact attenuation coefficients
// Gamma_s(b) = sigma_eff (1 - i eta) exp[-b^2/(2 B s hbarc^2)] / (4 pi B hbarc^2)
// w_LTA = 1 - r + r exp[-sigma_3 A T(b)/2]
// [REFERENCE: Guzey, Strikman and Zhalov, Eur. Phys. J. C 74 (2014) 2942, Eq. (16)]
MPhoto::AttenuationKernel MPhoto::PrepareAttenuation(const PhotoProfile& profile) const {
  if (model_ != PhotoModel::Impulse) { ValidateProfile(profile); }
  AttenuationKernel kernel;
  kernel.isospin = profile.isospin;
  if (model_ == PhotoModel::Impulse) { return kernel; }
  constexpr double hbarc = PDG::GeV2fm;
  if (model_ == PhotoModel::LTA) {
    if (shadow_model_ == nullptr) { throw std::logic_error("MPhoto: leading twist model is unavailable"); }
    const ShadowXS xs = shadow_model_->CrossSections(profile.x, profile.scale2);
    kernel.strength   = 0.05 * xs.sigma3 * target_->A();
    kernel.direct     = 1.0 - xs.ratio;
    kernel.rescatter  = xs.ratio;
    return kernel;
  }
  if (std::fpclassify(profile.sigma_eff) == FP_ZERO) { return kernel; }
  const double slope = profile.slope * hbarc * hbarc;
  kernel.strength    = 0.1 * profile.sigma_eff * std::complex<double>(1.0, -profile.eta) / (4.0 * math::PI * slope);
  std::vector<double>        scale_scratch;
  const std::vector<double>& scale = CrossSectionScales(profile, scale_scratch);
  kernel.production                = scale;
  kernel.inverse_width.resize(scale.size(), 0.0);
  for (const auto& i : indices(scale)) { kernel.inverse_width[i] = 1.0 / (2.0 * slope * scale[i]); }
  return kernel;
}

// Factor J_c(q) = sum_i exp(i q.r_i) w_ci with q-independent geometry in w_ci
// w_ci = prod_{j downstream i} [1 - Gamma_s(b_ij)]
// [REFERENCE: Frankfurt et al., Phys. Rev. C93 (2016) 055202]
std::complex<double> MPhoto::ConfigAttenuation(const MConfig& config, const std::size_t produced,
                                               const double inverse_width, const PhotonDirection direction,
                                               const AttenuationKernel& kernel) const {
  if (kernel.inverse_width.empty()) { return {1.0, 0.0}; }
  std::complex<double> attenuation = {1.0, 0.0};
  const bool           positive_z  = DirectionIndex(direction) == 0;
  for (const double radius2 : config.DownstreamR2(produced, positive_z)) {
    const double profile = std::exp(-radius2 * inverse_width);
    attenuation *= std::complex<double>(1.0, 0.0) - kernel.strength * profile;
    if (std::norm(attenuation) < std::numeric_limits<double>::min()) { return {0.0, 0.0}; }
  }
  return attenuation;
}

// Compute coherent and incoherent transitions from one full momentum transfer
// C(q) = (<J(q)>_optical, sqrt[Var[J(q)]_config])
PhotoTransition MPhoto::Factors(const PhotoProfile& profile, const double qx, const double qy, const double qz,
                                const PhotonDirection direction, const MConfigBank* bank, const MHotSpot* hotspot,
                                const bool with_current, const CoherenceType sector) const {
  CurrentStat optical = OpticalStat(profile, std::hypot(qx, qy), qz, direction);
  if (bank == nullptr && sector != CoherenceType::Coherent &&
      model_ == PhotoModel::Glauber && profile.sigma_eff > 0.0) {
    optical.variance = GlauberVariance(profile, std::hypot(qx, qy), qz);
  }
  if (bank == nullptr) { return {optical.mean, std::sqrt(std::max(0.0, optical.variance))}; }
  auto        sample = ShadowCurrents(profile, *bank, qx, qy, qz, direction, hotspot);
  CurrentStat stat   = CurrentMoments(sample);
  if (!with_current) { return {optical.mean, std::sqrt(std::max(0.0, stat.variance))}; }
  // Fix the coherent mean without changing centered configuration fluctuations
  for (auto& current : sample) { current += optical.mean - stat.mean; }
  stat.mean   = optical.mean;
  stat.second = std::norm(stat.mean) + stat.variance;
  return {optical.mean, std::sqrt(std::max(0.0, stat.variance)), PhotoCurrent{stat, std::move(sample)}};
}

// Compute nucleon-center and hotspot shadowed currents in that order
// J_c(q) = sum_i tau_i exp(i q.r_i/hbarc) <s prod_{j downstream i} [1 - Gamma_s(b_ij)]>_s
// tau_p = 1, tau_n = +1 for isoscalar and -1 for isovector production
// Attenuation keeps the produced vector state fixed, with no interchannel rescattering
// [REFERENCE: Frankfurt, Strikman and Zhalov, Phys. Lett. B540 (2002) 220]
std::array<std::complex<double>, 2> MPhoto::ShadowPair(const MConfig& config, const std::size_t sample, const double qx,
                                                       const double qy, const double qz,
                                                       const PhotonDirection direction, const MHotSpot* hotspot,
                                                       const AttenuationKernel& kernel) const {
  const std::array<double, 3> q = {qx, qy, qz};
  static_cast<void>(DirectionIndex(direction));
  constexpr double                    hbarc   = PDG::GeV2fm;
  std::array<std::complex<double>, 2> current = {0.0, 0.0};
  for (const auto& nucleon : indices(config.Nucleons())) {
    const auto&          center      = config.Nucleons()[nucleon];
    const double         phase       = gra::BilinearProduct(q, center.x) / hbarc;
    std::complex<double> attenuation = 0.0;
    if (model_ == PhotoModel::LTA) {
      attenuation = std::exp(-kernel.strength.real() * Thickness(std::hypot(center.x[0], center.x[1])));
    } else if (kernel.inverse_width.empty()) {
      attenuation = 1.0;
    } else {
      for (const auto& i : indices(kernel.inverse_width)) {
        attenuation +=
            kernel.production[i] * ConfigAttenuation(config, nucleon, kernel.inverse_width[i], direction, kernel);
      }
      attenuation /= static_cast<double>(kernel.inverse_width.size());
    }
    const double sign = kernel.isospin == PhotoIsospin::Isoscalar || center.type == NucleonType::Proton ? 1.0 : -1.0;
    const std::complex<double> contribution = sign * std::polar(1.0, phase) * (kernel.direct + kernel.rescatter * attenuation);
    current[0] += contribution;
    current[1] += hotspot != nullptr ? contribution * hotspot->Factor(sample, nucleon, qx, qy) : contribution;
  }
  return current;
}

// Compute one nucleon-center shadowed configuration current
// J_c(q) = sum_i exp(i q.r_i/hbarc) w_ci
std::complex<double> MPhoto::ShadowCurrent(const PhotoProfile& profile, const MConfig& config, const double qx,
                                           const double qy, const double qz, const PhotonDirection direction) const {
  if (config.Nucleons().size() != target_->A() || !gra::AllFinite(std::array<double, 3>{qx, qy, qz})) {
    throw std::invalid_argument("MPhoto::ShadowCurrent: configuration or momentum does not match");
  }
  const AttenuationKernel kernel = PrepareAttenuation(profile);
  return ShadowPair(config, 0, qx, qy, qz, direction, nullptr, kernel)[0];
}

// Compute the physical hotspot currents of one fixed bank
// J_c^phys = J_c^N + J_c^H - <H> J_c^N
// [REFERENCE: Good and Walker, Phys. Rev. 120 (1960) 1857]
std::vector<std::complex<double>> MPhoto::ShadowCurrents(const PhotoProfile& profile, const MConfigBank& bank,
                                                         const double qx, const double qy, const double qz,
                                                         const PhotonDirection direction,
                                                         const MHotSpot*       hotspot) const {
  if (bank.Nucleus().ID().pdg != target_->ID().pdg || bank.Size() == 0 ||
      !gra::AllFinite(std::array<double, 3>{qx, qy, qz})) {
    throw std::invalid_argument("MPhoto::ShadowCurrents: configuration bank does not match");
  }
  std::vector<std::complex<double>> current(bank.Size(), 0.0);
  if (hotspot != nullptr && (hotspot->ConfigCount() != bank.Size() || hotspot->NucleonCount() != target_->A())) {
    throw std::invalid_argument("MPhoto::ShadowCurrents: hotspot bank does not match");
  }
  const double            hotspot_mean  = hotspot != nullptr ? hotspot->Mean(qx, qy) : 1.0;
  const AttenuationKernel kernel        = PrepareAttenuation(profile);
  std::complex<double>    residual_mean = 0.0;
  for (const auto& sample : indices(current)) {
    const auto pair = ShadowPair(bank.At(sample), sample, qx, qy, qz, direction, hotspot, kernel);
    current[sample] = pair[0];
    if (hotspot != nullptr) {
      // Add H - <H> so the ensemble coherent current stays unchanged
      const std::complex<double> residual = pair[1] - hotspot_mean * pair[0];
      current[sample] += residual;
      residual_mean += residual;
    }
  }
  if (hotspot != nullptr) {
    residual_mean /= static_cast<double>(current.size());
    for (auto& value : current) { value -= residual_mean; }
  }
  if (!gra::AllFinite(current)) { throw AmplitudeFailure("MPhoto::ShadowCurrents: non-finite current"); }
  return current;
}

// Compute the shadowed current moments of one fixed configuration bank
// <J> = <J_c^N>, Var_u[J] = N_c/(N_c - 1) [<|J_c^phys|^2> - |<J_c^phys>|^2]
CurrentStat MPhoto::ShadowStat(const PhotoProfile& profile, const MConfigBank& bank, const double qx, const double qy,
                               const double qz, const PhotonDirection direction, const MHotSpot* hotspot) const {
  return CurrentMoments(ShadowCurrents(profile, bank, qx, qy, qz, direction, hotspot));
}

}  // namespace gra::nuclear
