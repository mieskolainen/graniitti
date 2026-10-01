// Spherical nuclear charge and matter densities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MDensity.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <stdexcept>
#include <utility>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;
using gra::math::ParallelFor;

namespace gra::nuclear {

// Construct and normalize one density on a fixed radial quadrature
MDensity::MDensity(DensityParam param) : param_(param) {
  Validate();
  Prepare();
}

// Compute the unnormalized two-parameter Fermi profile
// f(r) = 1 / [1 + exp((r - R) / a)]
double MDensity::Profile(const double r) const {
  const double x = (r - param_.radius) / param_.skin;
  if (x >= 0.0) {
    const double inverse = std::exp(-x);
    return inverse / (1.0 + inverse);
  }
  return 1.0 / (1.0 + std::exp(x));
}

// Validate all density parameters and the finite tail limit
void MDensity::Validate() {
  if (!std::isfinite(param_.radius) || !std::isfinite(param_.skin) || !std::isfinite(param_.r_max) ||
      !(param_.radius > 0.0) || !(param_.skin > 0.0)) {
    throw std::invalid_argument("MDensity: radius and skin must be finite and positive");
  }
  if (param_.nodes < 32 || param_.nodes > 4096) { throw std::invalid_argument("MDensity: nodes must be in [32,4096]"); }
  if (param_.cdf_nodes < 32) { throw std::invalid_argument("MDensity: invalid radial sampling controls"); }
  if (!std::isfinite(param_.form_q_max) || !(param_.form_q_max > 0.0) || param_.form_nodes < 4 ||
      param_.form_nodes > 262145 || !std::isfinite(param_.form_abs_tol) || !(param_.form_abs_tol > 0.0) ||
      param_.form_abs_tol > 1.0e-3) {
    throw std::invalid_argument("MDensity: invalid form-factor controls");
  }
  if (!(param_.r_max > param_.radius)) { throw std::invalid_argument("MDensity: r_max must exceed the radius"); }
}

// Prepare normalization, moments and inverse radial sampling data
// N = 4 pi int r^2 f(r) dr, <r^2> = int r^4 f(r) dr / int r^2 f(r) dr
void MDensity::Prepare() {
  auto rule = math::GaussLegendreRule(param_.nodes, 0.0, param_.r_max);
  r_node_   = std::move(rule.first);
  r_weight_ = std::move(rule.second);

  double radial_norm = 0.0;
  double radial_r2   = 0.0;
  for (const auto &i : indices(r_node_)) {
    const double r2    = r_node_[i] * r_node_[i];
    const double value = r_weight_[i] * r2 * Profile(r_node_[i]);
    radial_norm += value;
    radial_r2 += value * r2;
  }
  norm_ = 4.0 * math::PI * radial_norm;
  if (!std::isfinite(norm_) || !(norm_ > 0.0) || !std::isfinite(radial_r2) || !(radial_r2 > 0.0)) {
    throw std::invalid_argument("MDensity: density normalization failed");
  }
  rms_ = std::sqrt(radial_r2 / radial_norm);
  PrepareForm();

  const std::size_t intervals = param_.cdf_nodes;
  cdf_node_.resize(intervals + 1, 0.0);
  cdf_.resize(intervals + 1, 0.0);
  const double step     = param_.r_max / static_cast<double>(intervals);
  double       previous = 0.0;
  for (const auto &i : indices(cdf_)) {
    if (i == 0) { continue; }
    const double r     = static_cast<double>(i) * step;
    const double value = r * r * Profile(r);
    cdf_node_[i]       = r;
    cdf_[i]            = cdf_[i - 1] + 0.5 * step * (previous + value);
    previous           = value;
  }
  const double cdf_norm = cdf_.back();
  if (!std::isfinite(cdf_norm) || !(cdf_norm > 0.0)) {
    throw std::invalid_argument("MDensity: radial sampling table failed");
  }
  for (double &value : cdf_) { value /= cdf_norm; }
  cdf_.back() = 1.0;
}

// Compute the direct radial-quadrature form factor
// F(q) = 4 pi / N int r^2 f(r) j_0(qr / hbarc) dr
double MDensity::FormDirect(const double q) const {
  if (q <= param_.form_q_max) { return FormRule(q, form_r_node_, form_r_weight_); }
  const double inverse = PDG::GeV2fm / q, factor = 4.0 * math::PI / norm_;
  // Twice integrate r f(r) by parts, with int |(r f)''| dr <= 2 + R/(2 a)
  const double bound = factor * std::pow(inverse, 3) * (2.0 + param_.radius / (2.0 * param_.skin));
  if (bound <= param_.form_abs_tol) {
    const double tail = Profile(param_.r_max), phase = q * param_.r_max / PDG::GeV2fm;
    return factor * inverse * inverse * tail *
           (-param_.r_max * std::cos(phase) +
            inverse * (1.0 - param_.r_max * (1.0 - tail) / param_.skin) * std::sin(phase));
  }
  // Reuse the prepared rule on shorter intervals to resolve every Fourier phase
  const std::size_t   intervals = static_cast<std::size_t>(std::ceil(1.5 * q / param_.form_q_max));
  const double        scale     = 1.0 / static_cast<double>(intervals);
  std::vector<double> node(form_r_node_.size()), weight(form_r_weight_);
  for (double &value : weight) { value *= scale; }
  double integral = 0.0;
  for (std::size_t j = 0; j < intervals; ++j) {
    for (const auto &i : indices(node)) { node[i] = (form_r_node_[i] + j * param_.r_max) * scale; }
    integral += FormRule(q, node, weight);
  }
  return integral;
}

// Evaluate the form factor with one explicit radial rule
// F(q) = 4 pi / N int r^2 f(r) sin(qr / hbarc) / (qr / hbarc) dr
double MDensity::FormRule(const double q, const std::vector<double> &node, const std::vector<double> &weight) const {
  if (node.empty() || node.size() != weight.size()) {
    throw std::logic_error("MDensity::FormRule: radial rule is invalid");
  }
  double integral = 0.0;
  for (const auto &i : indices(node)) {
    const double r      = node[i];
    const double phase  = q * r / PDG::GeV2fm;
    const double bessel = std::fpclassify(phase) == FP_ZERO ? 1.0 : std::sin(phase) / phase;
    integral += weight[i] * r * r * Profile(r) * bessel;
  }
  const double value = 4.0 * math::PI * integral / norm_;
  if (!std::isfinite(value)) { throw std::runtime_error("MDensity::FormRule: non-finite form factor"); }
  return value;
}

// Compute a radial node count resolving one oscillatory form factor
std::size_t MDensity::FormNodeCount(const double q, const double scale) const {
  const double phase    = std::abs(q) * param_.r_max / PDG::GeV2fm;
  const double required = std::ceil(scale * phase) + 32.0;
  if (!std::isfinite(required) || required > 4096.0) {
    throw std::invalid_argument("MDensity: form-factor momentum exceeds the resolved range");
  }
  return std::max<std::size_t>(param_.nodes, static_cast<std::size_t>(required));
}

// Interpolate the cached form factor on its uniform momentum grid
double MDensity::FormCached(const double q) const {
  // Interpolate F - 1 in q^2 near the origin to preserve the even charge-radius expansion
  if (q < form_step_) {
    const std::array<double, 4> node{0.0, 1.0, 4.0, 9.0};
    const std::array<double, 4> value{0.0, form_[1] - 1.0, form_[2] - 1.0, form_[3] - 1.0};
    return 1.0 + math::CubicLagrangeInterpolate(node, value, math::pow2(q / form_step_));
  }
  const auto stencil = math::UniformCubicWeights(form_.size(), 0.0, param_.form_q_max, q);
  double     value   = 0.0;
  for (const auto &i : indices(stencil.index)) { value += stencil.weight[i] * form_[stencil.index[i]]; }
  return value;
}

// Prepare and validate the immutable form-factor table
void MDensity::PrepareForm() {
  const std::size_t form_count = FormNodeCount(param_.form_q_max, 0.5);
  auto              form_rule  = math::GaussLegendreRule(static_cast<unsigned int>(form_count), 0.0, param_.r_max);
  form_r_node_                 = std::move(form_rule.first);
  form_r_weight_               = std::move(form_rule.second);

  form_step_ = param_.form_q_max / static_cast<double>(param_.form_nodes - 1);
  form_.resize(param_.form_nodes, 0.0);
  ParallelFor(form_.size(), [&](const std::size_t i) { form_[i] = FormDirect(static_cast<double>(i) * form_step_); });
  form_.front() = 1.0;

  const std::size_t   check_count = FormNodeCount(param_.form_q_max, 0.75);
  const auto          check_rule  = math::GaussLegendreRule(static_cast<unsigned int>(check_count), 0.0, param_.r_max);
  std::vector<double> error(form_.size() - 1, 0.0);
  ParallelFor(error.size(), [&](const std::size_t i) {
    const double q = (static_cast<double>(i) + 0.5) * form_step_;
    error[i]       = std::abs(FormCached(q) - FormRule(q, check_rule.first, check_rule.second));
  });
  const double max_error = *std::max_element(error.cbegin(), error.cend());
  if (!std::isfinite(max_error) || max_error > param_.form_abs_tol) {
    throw std::invalid_argument("MDensity: form-factor interpolation tolerance not reached");
  }
}

// Compute the normalized three-dimensional density in fm^-3
// rho(r) = f(r) / N with int rho(r) d^3r = 1
double MDensity::Rho(const double r) const {
  if (!std::isfinite(r) || r < 0.0) {
    throw std::invalid_argument("MDensity::Rho: radius must be finite and nonnegative");
  }
  if (r >= param_.r_max) { return 0.0; }
  return Profile(r) / norm_;
}

// Compute the normalized spherical form factor for momentum in GeV
// F(q) = int rho(r) exp(i q.r / hbarc) d^3r
double MDensity::Form(const double q) const {
  if (!std::isfinite(q)) { throw std::invalid_argument("MDensity::Form: momentum must be finite"); }
  const double q_abs = std::abs(q);
  if (std::fpclassify(q_abs) == FP_ZERO) { return 1.0; }
  if (q_abs <= param_.form_q_max) { return FormCached(q_abs); }
  return FormDirect(q_abs);
}

// Compute the normalized transverse thickness in fm^-2
// T(b) = int rho(sqrt(b^2 + z^2)) dz
double MDensity::Thick(const double b) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MDensity::Thick: impact parameter must be finite and nonnegative");
  }
  if (b >= param_.r_max) { return 0.0; }
  const double z_max    = std::sqrt(param_.r_max * param_.r_max - b * b);
  const double scale    = z_max / param_.r_max;
  double       integral = 0.0;
  for (const auto &i : indices(r_node_)) {
    const double z = scale * r_node_[i];
    const double r = std::sqrt(b * b + z * z);
    integral += scale * r_weight_[i] * Profile(r);
  }
  const double value = 2.0 * integral / norm_;
  if (!std::isfinite(value) || value < 0.0) { throw std::runtime_error("MDensity::Thick: non-finite thickness"); }
  return value;
}

// Compute the charge fraction inside one transverse cylinder
// C(b) = 4 pi [int_0^b r^2 rho(r) dr + b^2 int_0^zmax z rho(hypot(b,z))/(hypot(b,z)+z) dz]
double MDensity::Cylinder(const double b) const {
  if (!std::isfinite(b) || b < 0.0) {
    throw std::invalid_argument("MDensity::Cylinder: radius must be finite and nonnegative");
  }
  if (b >= param_.r_max) { return 1.0; }
  const double radial_scale = b / param_.r_max;
  const double z_max        = std::sqrt((param_.r_max - b) * (param_.r_max + b));
  const double z_scale      = z_max / param_.r_max;
  const double profile      = Profile(b);
  // Integrate the constant cap density analytically to resolve arbitrarily small cylinders
  double integral = profile * b * b * (param_.r_max + z_max * z_max / (param_.r_max + z_max) - b) / 3.0;
  for (const auto& i : indices(r_node_)) {
    const double r = radial_scale * r_node_[i], z = z_scale * r_node_[i], radius = std::hypot(b, z);
    // Resolve the spherical cap with z = sqrt(r^2 - b^2), retaining C(b) proportional to b^2
    integral += r_weight_[i] *
                (radial_scale * r * r * Profile(r) + z_scale * b * b * z / (radius + z) * (Profile(radius) - profile));
  }
  const double fraction = 4.0 * math::PI * integral / norm_;
  if (!std::isfinite(fraction)) { throw std::runtime_error("MDensity::Cylinder: non-finite charge fraction"); }
  return std::clamp(fraction, 0.0, 1.0);
}

// Map one uniform variate to the normalized radial density
// u = 4 pi int_0^r r'^2 rho(r') dr'
double MDensity::Radius(const double u) const {
  if (u < 0.0 || u > 1.0) { throw std::invalid_argument("MDensity::Radius: variate must be in [0,1]"); }
  return math::LinearInterpolateValidatedGrid(cdf_, cdf_node_, u);
}

}  // namespace gra::nuclear
