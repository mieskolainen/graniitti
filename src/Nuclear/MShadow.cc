// Leading twist nuclear gluon shadowing
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MShadow.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <sstream>
#include <stdexcept>
#include <utility>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MCache.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonZip.h"
#include "json.hpp"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {
namespace {

// Compute the persistent leading twist table version
constexpr std::size_t ShadowCacheVersion() { return 2; }

// Validate one perturbatively matched inclusive and diffractive PDF pair
void ValidatePDFPair(const ShadowParam &param, const LHAPDF::PDF &pdf, const LHAPDF::PDF &dpdf) {
  if (!pdf.hasFlavor(PDG::PDG_gluon) || !dpdf.hasFlavor(PDG::PDG_gluon) || !pdf.hasAlphaS() || !dpdf.hasAlphaS()) {
    throw std::invalid_argument("MShadow: PDFs require gluons and alpha_s");
  }

  const auto &pdf_info  = pdf.info();
  const auto &dpdf_info = dpdf.info();
  // LTA consistency requires matching PDF and alpha_s evolution orders
  // Small differences between fitted alpha_s(MZ) values remain a PDF input uncertainty
  const int pdf_order        = pdf_info.get_entry_as<int>("OrderQCD", -1);
  const int dpdf_order       = dpdf_info.get_entry_as<int>("OrderQCD", -1);
  const int pdf_alpha_order  = pdf_info.get_entry_as<int>("AlphaS_OrderQCD", -1);
  const int dpdf_alpha_order = dpdf_info.get_entry_as<int>("AlphaS_OrderQCD", -1);
  if (pdf_order < 0 || pdf_order != dpdf_order || pdf_alpha_order < 0 || pdf_alpha_order != dpdf_alpha_order) {
    throw std::invalid_argument("MShadow: inclusive and diffractive PDFs use incompatible QCD evolution");
  }

  // Require the validated (x, Q2) domain to lie inside the common LHAPDF grid
  const double beta_min   = param.x_min / param.x_max;
  const double q2_min     = std::max(math::pow2(pdf_info.get_entry_as<double>("QMin", 0.0)),
                                     math::pow2(dpdf_info.get_entry_as<double>("QMin", 0.0)));
  const double q2_max     = std::min(math::pow2(pdf_info.get_entry_as<double>("QMax", 0.0)),
                                     math::pow2(dpdf_info.get_entry_as<double>("QMax", 0.0)));
  const double pdf_x_min  = pdf_info.get_entry_as<double>("XMin", 1.0);
  const double pdf_x_max  = pdf_info.get_entry_as<double>("XMax", 0.0);
  const double dpdf_x_min = dpdf_info.get_entry_as<double>("XMin", 1.0);
  const double dpdf_x_max = dpdf_info.get_entry_as<double>("XMax", 0.0);
  const double tolerance  = 1.0e-12;
  if (param.x_min < pdf_x_min || param.x_max > pdf_x_max || beta_min < dpdf_x_min || dpdf_x_max < 1.0 - tolerance ||
      param.q2_min < q2_min * (1.0 - tolerance) || param.q2_max > q2_max * (1.0 + tolerance)) {
    throw std::invalid_argument("MShadow: model domain exceeds the common PDF grid");
  }
}

}  // namespace

// Construct immutable PDF and interpolation tables
MShadow::MShadow(ShadowParam param) : param_(std::move(param)) { Prepare(); }

// Validate controls and prepare the persistent interpolation table
void MShadow::Prepare() {
  if (param_.pdf_set.empty() || param_.dpdf_set.empty() || param_.pdf_member < 0 || param_.dpdf_member < 0 ||
      !std::isfinite(param_.alpha0) || !(param_.alpha0 > 1.0) || !std::isfinite(param_.dpdf_diss) ||
      param_.dpdf_diss < 1.0 || !std::isfinite(param_.alpha_prime) || param_.alpha_prime < 0.0 ||
      !std::isfinite(param_.flux_b) || !(param_.flux_b > 0.0) || !std::isfinite(param_.b_diff) ||
      !(param_.b_diff > 0.0) || !std::isfinite(param_.x_min) || !std::isfinite(param_.x_max) ||
      !std::isfinite(param_.q2_min) || !std::isfinite(param_.q2_max) || !std::isfinite(param_.sigma3_anchor_x) ||
      !std::isfinite(param_.sigma3_fade_x) || !std::isfinite(param_.sigma3_power) || !(param_.x_min > 0.0) ||
      !(param_.x_max > param_.x_min) || !(param_.q2_min > 0.0) || !(param_.q2_max > param_.q2_min) ||
      param_.sigma3_anchor_x < param_.x_min || !(param_.sigma3_fade_x > param_.sigma3_anchor_x) ||
      !(param_.sigma3_fade_x < param_.x_max) || !(param_.sigma3_power > 0.0) || param_.x_nodes < 17 ||
      param_.scale_nodes < 9 || param_.integral_nodes < 16) {
    throw std::invalid_argument("MShadow: invalid leading twist controls");
  }
  pdf_  = pdf_store_.GetPDF(param_.pdf_set, param_.pdf_member);
  dpdf_ = pdf_store_.GetPDF(param_.dpdf_set, param_.dpdf_member);
  ValidatePDFPair(param_, *pdf_, *dpdf_);

  log_x_      = math::linspace(std::log(param_.x_min), std::log(param_.x_max), param_.x_nodes);
  log_scale2_ = math::linspace(std::log(param_.q2_min), std::log(param_.q2_max), param_.scale_nodes);

  const std::string key      = CacheKey();
  const std::string filename = NuclearCacheFilename("LTA", key);
  MCacheLock        lock(filename + ".lock");
  const bool        loaded = ReadCache(filename, key);
  if (loaded) {
    std::cout << "Loaded leading twist shadowing cache: " << filename << std::endl;
    return;
  }

  sigma2_.resize(log_x_.size() * log_scale2_.size(), 0.0);
  ParallelFor(
      sigma2_.size(),
      [&](const std::size_t index) {
        const std::size_t ix = index / log_scale2_.size();
        const std::size_t iq = index % log_scale2_.size();
        sigma2_[index]       = Sigma2(std::exp(log_x_[ix]), std::exp(log_scale2_[iq]));
      },
      "Leading twist shadowing table");
  if (!gra::AllFinite(sigma2_)) { throw std::runtime_error("MShadow: non-finite shadowing table"); }
  if (lock.Acquired()) {
    WriteCache(filename, key);
    std::cout << "Saved leading twist shadowing cache: " << filename << std::endl;
  }
}

// Evaluate the two nucleon effective cross section directly
// sigma_2 = 16 pi / [(1 + eta^2) xg(x)] int dx_P beta g_D(4)(beta,x_P,t_min)
// GKG18 absorbs the flux normalization into g_P, dpdf_diss removes proton dissociation
// [REFERENCE: arXiv:1106.2091, Eq. (52), Sec. 5.1.1; arXiv:1802.01363, Eq. (8), Sec. III B]
double MShadow::Sigma2(const double x, const double scale2) const {
  if (x >= param_.x_max) { return 0.0; }
  const double proton = pdf_->xfxQ2(21, x, scale2);
  if (!(proton > 0.0) || !std::isfinite(proton)) {
    throw std::runtime_error("MShadow: invalid inclusive gluon density");
  }
  const auto rule     = math::GaussLegendreRule(param_.integral_nodes, std::log(x), std::log(param_.x_max));
  double     integral = 0.0;
  for (const auto &i : indices(rule.first)) {
    const double xp         = std::exp(rule.first[i]);
    const double beta       = x / xp;
    const double t_min      = -math::pow2(PDG::mp * xp) / (1.0 - xp);
    const double slope      = param_.flux_b - 2.0 * param_.alpha_prime * std::log(xp);
    const double flux       = std::pow(xp, 1.0 - 2.0 * param_.alpha0) * std::exp(slope * t_min) / param_.dpdf_diss;
    const double beta_gluon = dpdf_->xfxQ2(21, beta, scale2);
    integral += rule.second[i] * xp * flux * beta_gluon;
  }
  const double eta        = std::tan(0.5 * math::PI * (param_.alpha0 - 1.0));
  const double sigma_gev2 = 16.0 * math::PI * integral / ((1.0 + eta * eta) * proton);
  const double sigma_mb   = sigma_gev2 * PDG::GeV2mb;
  if (!std::isfinite(sigma_mb) || sigma_mb < 0.0) {
    throw std::runtime_error("MShadow: invalid two nucleon cross section");
  }
  return sigma_mb;
}

// Interpolate the two nucleon effective cross section
double MShadow::Interpolate(const double x, const double scale2) const {
  const auto        x_position = math::LocateCell(log_x_, std::log(x));
  const auto        q_position = math::LocateCell(log_scale2_, std::log(scale2));
  const std::size_t stride     = log_scale2_.size();
  const auto        value      = [&](const std::size_t ix, const std::size_t iq) { return sigma2_[ix * stride + iq]; };
  const double      lower =
      value(x_position.index, q_position.index) +
      q_position.fraction * (value(x_position.index, q_position.index + 1) - value(x_position.index, q_position.index));
  const double upper = value(x_position.index + 1, q_position.index) +
                       q_position.fraction * (value(x_position.index + 1, q_position.index + 1) -
                                              value(x_position.index + 1, q_position.index));
  return lower + x_position.fraction * (upper - lower);
}

// Compute effective cross sections at one gluon x and scale
// sigma_3 = sigma_2 at small x, a power continuation at intermediate x, then a linear fade
// r = sigma_2 / sigma_3, sigma_3^in = sigma_3 - sigma_3^2 / (16 pi B_diff)
// [REFERENCE: Frankfurt, Guzey and Strikman, Phys. Rept. 512 (2012) 255]
ShadowXS MShadow::CrossSections(const double x, const double scale2) const {
  if (!std::isfinite(x) || !std::isfinite(scale2) || x < param_.x_min || !(x < 1.0) || !(scale2 > 0.0) ||
      scale2 > param_.q2_max) {
    throw std::invalid_argument("MShadow: x or scale is outside the model domain");
  }
  // Freeze low scale LTA moments at the configured input scale
  const double q2 = std::max(scale2, param_.q2_min);
  ShadowXS     out;
  // The diffractive integral vanishes above x_max, leaving the impulse response
  if (x >= param_.x_max) { return out; }
  out.sigma2 = Interpolate(x, q2);
  // Continue the FGS10_H gluon moment from its black-disk lower bound
  const double anchor = Interpolate(param_.sigma3_anchor_x, q2);
  if (x <= param_.sigma3_anchor_x) {
    out.sigma3 = out.sigma2;
  } else if (x <= param_.sigma3_fade_x) {
    out.sigma3 = anchor * std::pow(param_.sigma3_anchor_x / x, param_.sigma3_power);
  } else {
    out.sigma3 = anchor * std::pow(param_.sigma3_anchor_x / param_.sigma3_fade_x, param_.sigma3_power) *
                 (param_.x_max - x) / (param_.x_max - param_.sigma3_fade_x);
  }
  if (!std::isfinite(out.sigma3) || !(out.sigma3 > 0.0) || out.sigma2 > out.sigma3) {
    std::ostringstream message;
    message << "MShadow: unphysical shadowing cross sections at x=" << x << ", scale2=" << scale2
            << " GeV2, sigma2=" << out.sigma2 << " mb, sigma3=" << out.sigma3 << " mb";
    throw std::runtime_error(message.str());
  }
  out.ratio     = out.sigma2 / out.sigma3;
  out.sigma3_in = out.sigma3 * (1.0 - out.sigma3 / (16.0 * math::PI * param_.b_diff * PDG::GeV2mb));
  // Elastic subtraction can round to zero in the continuous impulse limit
  if (!std::isfinite(out.sigma3_in) || !(out.sigma3_in > 0.0) || out.sigma3_in > out.sigma3) {
    throw std::runtime_error("MShadow: unphysical inelastic cross section");
  }
  return out;
}

// Construct the complete table identity
std::string MShadow::CacheKey() const {
  std::ostringstream stream;
  stream << ShadowCacheVersion() << ';' << param_.pdf_set << ';' << param_.pdf_member << ';' << param_.dpdf_set << ';'
         << param_.dpdf_member << ';' << LHAPDF::version() << ';' << pdf_->dataversion() << ';' << dpdf_->dataversion()
         << ';' << std::hexfloat << param_.alpha0 << ';' << param_.alpha_prime << ';' << param_.flux_b << ';'
         << param_.dpdf_diss << ';' << param_.b_diff << ';' << param_.x_min << ';' << param_.x_max << ';'
         << param_.q2_min << ';' << param_.q2_max << ';' << param_.sigma3_anchor_x << ';' << param_.sigma3_fade_x << ';'
         << param_.sigma3_power << ';' << param_.x_nodes << ';' << param_.scale_nodes << ';' << param_.integral_nodes
         << ';';
  return stream.str();
}

// Load one compressed leading twist table
bool MShadow::ReadCache(const std::string &filename, const std::string &key) {
  try {
    nlohmann::json cache;
    return ReadNuclearCache(filename, ShadowCacheVersion(), "leading_twist_shadowing", key, cache) &&
           MJsonZip::DecompressVector(cache.at("sigma2"), sigma2_, log_x_.size() * log_scale2_.size());
  } catch (const std::exception &) { return false; }
}

// Save one compressed leading twist table
void MShadow::WriteCache(const std::string &filename, const std::string &key) const {
  auto cache      = NuclearCache(ShadowCacheVersion(), "leading_twist_shadowing", key);
  cache["sigma2"] = MJsonZip::CompressVector(sigma2_);
  PublishNuclearCache(filename, cache);
}

}  // namespace gra::nuclear
