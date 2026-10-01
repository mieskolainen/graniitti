// Generic probability distribution functions and moment matching
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MDistribution.h"
#include "Graniitti/Math/MSpecialFunctions.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>

namespace gra::math {

// Compute the standard normal quantile using a rational approximation
// z = Phi^{-1}(p)
double NormalQuantile(const double probability) {
  constexpr double low = 0.02425;
  constexpr std::array<double, 6> a = {
      -3.969683028665376e+01, 2.209460984245205e+02,  -2.759285104469687e+02,
      1.383577518672690e+02,  -3.066479806614716e+01, 2.506628277459239e+00};
  constexpr std::array<double, 5> b = {
      -5.447609879822406e+01, 1.615858368580409e+02, -1.556989798598866e+02,
      6.680131188771972e+01, -1.328068155288572e+01};
  constexpr std::array<double, 6> c = {
      -7.784894002430293e-03, -3.223964580411365e-01, -2.400758277161838e+00,
      -2.549732539343734e+00, 4.374664141464968e+00,  2.938163982698783e+00};
  constexpr std::array<double, 4> d = {
      7.784695709041462e-03, 3.224671290700398e-01, 2.445134137142996e+00,
      3.754408661907416e+00};
  if (!(probability > 0.0) || !(probability < 1.0)) {
    throw std::invalid_argument("NormalQuantile: probability is invalid");
  }
  if (probability < low) {
    const double q = std::sqrt(-2.0 * std::log(probability));
    return (((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
            c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  if (probability > 1.0 - low) {
    const double q = std::sqrt(-2.0 * std::log1p(-probability));
    return -(((((c[0] * q + c[1]) * q + c[2]) * q + c[3]) * q + c[4]) * q +
             c[5]) /
           ((((d[0] * q + d[1]) * q + d[2]) * q + d[3]) * q + 1.0);
  }
  const double q = probability - 0.5;
  const double r = q * q;
  return (((((a[0] * r + a[1]) * r + a[2]) * r + a[3]) * r + a[4]) * r + a[5]) *
         q /
         (((((b[0] * r + b[1]) * r + b[2]) * r + b[3]) * r + b[4]) * r + 1.0);
}

// Validate Gamma distribution iteration controls once at construction
GammaNumerics::GammaNumerics(double tolerance, unsigned int max_iter, double normal_min)
    : cdf_tol(tolerance), cdf_max_iter(max_iter), normal_shape_min(normal_min) {
  if (!std::isfinite(cdf_tol) || !(cdf_tol > 0.0) || !(cdf_tol < 1.0) ||
      cdf_max_iter == 0 || !std::isfinite(normal_shape_min) || !(normal_shape_min > 0.0)) {
    throw std::invalid_argument("Gamma distribution controls are invalid");
  }
}

// Compute the regularized lower incomplete Gamma function
// P(a,x) = gamma(a,x)/Gamma(a)
double GammaCDF(const double shape, const double value,
                const GammaNumerics &numerics) {
  if (!(shape > 0.0) || value < 0.0 || !std::isfinite(shape) ||
      !std::isfinite(value)) {
    throw std::invalid_argument("GammaCDF: invalid input");
  }
  if (std::fpclassify(value) == FP_ZERO) {
    return 0.0;
  }
  if (value < shape + 1.0) {
    // Use reentrant Gamma evaluation when probability tables are prepared concurrently
    const double exponent = -value + shape * std::log(value) - gra::math::LogGamma(shape + 1.0);
    double term = 1.0;
    double       sum           = term;
    double       shifted_shape = shape;
    for (unsigned int iteration = 0; iteration < numerics.cdf_max_iter; ++iteration) {
      shifted_shape += 1.0;
      term *= value / shifted_shape;
      sum += term;
      if (std::abs(term) <= numerics.cdf_tol * std::abs(sum)) {
        const double cdf = std::exp(std::log(sum) + exponent);
        if (!std::isfinite(cdf)) { throw std::runtime_error("GammaCDF: non-finite series"); }
        return std::clamp(cdf, 0.0, 1.0);
      }
    }
  } else {
    const double exponent = -value + shape * std::log(value) - gra::math::LogGamma(shape);
    const double floor    = std::numeric_limits<double>::min() / numerics.cdf_tol;
    double       offset   = value + 1.0 - shape;
    double       c        = 1.0 / floor;
    double       d        = 1.0 / std::max(std::abs(offset), floor);
    double       fraction = d;
    for (unsigned int iteration = 0; iteration < numerics.cdf_max_iter; ++iteration) {
      const double index       = static_cast<double>(iteration) + 1.0;
      const double coefficient = -index * (index - shape);
      offset += 2.0;
      d = coefficient * d + offset;
      if (std::abs(d) < floor) { d = std::copysign(floor, d); }
      c = offset + coefficient / c;
      if (std::abs(c) < floor) { c = std::copysign(floor, c); }
      d                   = 1.0 / d;
      const double change = d * c;
      fraction *= change;
      if (std::abs(change - 1.0) <= numerics.cdf_tol) {
        const double cdf = 1.0 - std::exp(exponent) * fraction;
        if (!std::isfinite(cdf)) { throw std::runtime_error("GammaCDF: non-finite continued fraction"); }
        return std::clamp(cdf, 0.0, 1.0);
      }
    }
  }
  throw std::runtime_error("GammaCDF: iteration did not converge");
}

// Compute one Gamma quantile with unit scale
// x = P^{-1}(a,p), with x ~= a[1 - 1/(9a) + Phi^{-1}(p)/(3sqrt(a))]^3 at large a
double GammaQuantile(const double probability, const double shape, const GammaNumerics &numerics) {
  if (!(probability > 0.0) || !(probability < 1.0) || !(shape > 0.0) ||
      !std::isfinite(shape)) {
    throw std::invalid_argument("GammaQuantile: invalid input");
  }
  if (shape >= numerics.normal_shape_min) {
    const double z = NormalQuantile(probability);
    const double base =
        1.0 - 1.0 / (9.0 * shape) + z / (3.0 * std::sqrt(shape));
    const double approximate = shape * base * base * base;
    if (approximate > 0.0 && std::isfinite(approximate)) { return approximate; }
  }
  double lower = 0.0;
  double upper = std::max(1.0, shape + 12.0 * std::sqrt(shape) + 12.0);
  while (GammaCDF(shape, upper, numerics) < probability) {
    upper *= 2.0;
  }
  for (unsigned int iteration = 0; iteration < numerics.cdf_max_iter;
       ++iteration) {
    const double middle = lower + 0.5 * (upper - lower);
    if (GammaCDF(shape, middle, numerics) < probability) {
      lower = middle;
    } else {
      upper = middle;
    }
    if (upper - lower <= numerics.cdf_tol * std::max(std::abs(upper), std::numeric_limits<double>::denorm_min())) {
      return lower + 0.5 * (upper - lower);
    }
    if (std::nextafter(lower, upper) >= upper) { return lower + 0.5 * (upper - lower); }
  }
  throw std::runtime_error("GammaQuantile: iteration did not converge");
}

}  // namespace gra::math
