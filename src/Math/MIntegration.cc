// Numerical integration rules and sampled integrals
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MIntegration.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"

namespace gra::math {

// Compute the composite Simpson or trapezoidal integration weight
// w_i = (1,4,2,...,4,1)/3 for even N and (1/2,1,...,1,1/2) for odd N
double CompositeWeight(const std::size_t index, const std::size_t intervals) {
  if (intervals == 0 || index > intervals) { throw std::invalid_argument("math::CompositeWeight: index outside grid"); }
  if (intervals % 2 == 0) {
    if (index == 0 || index == intervals) { return 1.0 / 3.0; }
    return ((index % 2) == 0 ? 2.0 : 4.0) / 3.0;
  }
  return (index == 0 || index == intervals) ? 0.5 : 1.0;
}

// Integrate x f(x) exactly when f is represented by linear grid interpolation
// I_i = Delta x_i [f_i(2x_i + x_{i+1}) + f_{i+1}(x_i + 2x_{i+1})]/6
double LinearRadialIntegral(const std::vector<double> &node, const std::vector<double> &value) {
  ValidateInterpolationGrid(node);
  if (value.size() != node.size() ||
      !std::all_of(value.cbegin(), value.cend(), [](const double entry) { return std::isfinite(entry); })) {
    throw std::invalid_argument("math::LinearRadialIntegral: non-finite or inconsistent grid");
  }
  double integral = 0.0;
  for (std::size_t i = 1; i < node.size(); ++i) {
    const double lower = node[i - 1];
    const double upper = node[i];
    integral += (upper - lower) * (value[i - 1] * (2.0 * lower + upper) + value[i] * (lower + 2.0 * upper)) / 6.0;
  }
  return integral;
}

// Compute midpoint periodic trapezoid nodes and scaled weights
// x_j = a + (j + 1/2)(b - a)/N, w_j = (b - a)/N
std::pair<std::vector<double>, std::vector<double>> PeriodicTrapzRule(const unsigned int n, const double a,
                                                                      const double b) {
  if (n == 0 || !std::isfinite(a) || !std::isfinite(b) || !(a < b)) {
    throw std::invalid_argument("math::PeriodicTrapzRule: invalid rule");
  }
  const double        step = (b - a) / static_cast<double>(n);
  std::vector<double> node(n, 0.0);
  std::vector<double> weight(n, step);
  for (std::size_t j = 0; j < n; ++j) { node[j] = a + (static_cast<double>(j) + 0.5) * step; }
  return {node, weight};
}

// Compute Gauss Legendre nodes and scaled weights on a finite interval
// x_i = (a+b)/2 + (b-a)z_i/2, w_i = (b-a)/[(1-z_i^2)P'_N(z_i)^2]
std::pair<std::vector<double>, std::vector<double>> GaussLegendreRule(const unsigned int n, const double a,
                                                                      const double b) {
  if (n == 0 || !std::isfinite(a) || !std::isfinite(b) || !(a < b)) {
    throw std::invalid_argument("math::GaussLegendreRule: invalid rule");
  }
  std::vector<double> node(n, 0.0);
  std::vector<double> weight(n, 0.0);
  const double        middle    = 0.5 * (b + a);
  const double        half      = 0.5 * (b - a);
  const unsigned int  roots     = (n + 1) / 2;
  constexpr double    tolerance = 1.0e-15;
  for (unsigned int i = 0; i < roots; ++i) {
    double z          = std::cos(PI * (static_cast<double>(i) + 0.75) / (static_cast<double>(n) + 0.5));
    double previous   = 0.0;
    double derivative = 0.0;
    do {
      double current = 1.0;
      double lower   = 0.0;
      for (unsigned int j = 1; j <= n; ++j) {
        const double lower2 = lower;
        lower               = current;
        current             = ((2.0 * j - 1.0) * z * lower - (j - 1.0) * lower2) / j;
      }
      derivative = n * (z * current - lower) / (z * z - 1.0);
      previous   = z;
      z          = previous - current / derivative;
    } while (std::abs(z - previous) > tolerance);
    node[i]           = middle - half * z;
    node[n - 1 - i]   = middle + half * z;
    weight[i]         = 2.0 * half / ((1.0 - z * z) * derivative * derivative);
    weight[n - 1 - i] = weight[i];
  }
  return {node, weight};
}

// Compute the closed grid Simpson 1/3 weight vector
// w = (1,4,2,4,...,2,4,1)
std::vector<double> Simpson13Weight(const unsigned int intervals) {
  if (intervals < 2 || intervals % 2 != 0) {
    throw std::invalid_argument("math::Simpson13Weight: interval count must be positive and even");
  }
  std::vector<double> weight(intervals + 1, 2.0);
  for (std::size_t i = 1; i < intervals; i += 2) { weight[i] = 4.0; }
  weight.front() = 1.0;
  weight.back()  = 1.0;
  return weight;
}

// Compute the closed grid Simpson 3/8 weight vector
// w_i = 1 at the boundaries, 2 for interior i mod 3 = 0, and 3 otherwise
std::vector<double> Simpson38Weight(const unsigned int intervals) {
  if (intervals < 3 || intervals % 3 != 0) {
    throw std::invalid_argument("math::Simpson38Weight: interval count must be a positive multiple of 3");
  }
  std::vector<double> weight(intervals + 1, 3.0);
  for (std::size_t i = 3; i < intervals; i += 3) { weight[i] = 2.0; }
  weight.front() = 1.0;
  weight.back()  = 1.0;
  return weight;
}

// Compute the closed grid Boole weight vector
// w = (7,32,12,32,14,...,32,12,32,7)
std::vector<double> BooleWeight(const unsigned int intervals) {
  if (intervals < 4 || intervals % 4 != 0) {
    throw std::invalid_argument("math::BooleWeight: interval count must be a positive multiple of 4");
  }
  std::vector<double> weight(intervals + 1, 14.0);
  for (std::size_t j = 0; j < intervals; j += 4) {
    weight[j + 1] = 32.0;
    weight[j + 2] = 12.0;
    weight[j + 3] = 32.0;
  }
  weight.front() = 7.0;
  weight.back()  = 7.0;
  return weight;
}

// Compute the tensor Simpson 1/3 weight matrix
// W_{ij} = w_i w_j
MMatrix<double> Simpson13Weight2D(const unsigned int rows, const unsigned int cols) {
  return gra::OuterProduct(Simpson13Weight(rows), Simpson13Weight(cols));
}

// Compute the tensor Simpson 3/8 weight matrix
// W_{ij} = w_i w_j
MMatrix<double> Simpson38Weight2D(const unsigned int rows, const unsigned int cols) {
  return gra::OuterProduct(Simpson38Weight(rows), Simpson38Weight(cols));
}

// Compute the tensor Boole weight matrix
// W_{ij} = w_i w_j
MMatrix<double> BooleWeight2D(const unsigned int rows, const unsigned int cols) {
  return gra::OuterProduct(BooleWeight(rows), BooleWeight(cols));
}

}  // namespace gra::math
