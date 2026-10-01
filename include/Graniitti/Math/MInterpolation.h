// Interpolation grid utilities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MINTERPOLATION_H
#define MINTERPOLATION_H

#include <algorithm>
#include <array>
#include <cmath>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MFloat.h"

namespace gra::math {

// Store one lower interpolation index and local coordinate
struct InterpolationCell {
  std::size_t index    = 0;
  double      fraction = 0.0;
};

// Validate one finite strictly increasing interpolation grid
void ValidateInterpolationGrid(const std::vector<double> &node);

// Locate one value on a finite strictly increasing interpolation grid
InterpolationCell LocateCell(const std::vector<double> &node, double value);

// Linearly interpolate on a finite nondecreasing grid whose ordering was already validated
// Repeated interior knots select their last value, allowing inversion of CDF plateaus
// f(x) = (1-t)f_i + t f_{i+1}, t = (x-x_i)/(x_{i+1}-x_i)
template <typename T>
inline T LinearInterpolateValidatedGrid(const std::vector<double> &node, const std::vector<T> &value, const double x) {
  if (!std::isfinite(x) || node.empty() || value.size() != node.size()) {
    throw std::invalid_argument("math::LinearInterpolate: empty or inconsistent grid");
  }
  if (x <= node.front()) { return value.front(); }
  if (x >= node.back()) { return value.back(); }
  const auto        upper    = std::upper_bound(node.cbegin(), node.cend(), x);
  const std::size_t hi       = std::distance(node.cbegin(), upper);
  const std::size_t lo       = hi - 1;
  const double      fraction = (x - node[lo]) / (node[hi] - node[lo]);
  T                 output   = value[lo] * (1.0 - fraction);
  if constexpr (requires { output.AddScaled(value[hi], fraction); }) {
    output.AddScaled(value[hi], fraction);
  } else {
    output += value[hi] * fraction;
  }
  return output;
}

// Linearly interpolate a value table on one validated ordered real grid
template <typename T>
inline T LinearInterpolate(const std::vector<double> &node, const std::vector<T> &value, const double x) {
  ValidateInterpolationGrid(node);
  return LinearInterpolateValidatedGrid(node, value, x);
}

// Store indices and weights for cubic interpolation on a uniform grid
struct UniformCubicStencil {
  std::array<std::size_t, 4> index{};
  std::array<double, 4>      weight{};
};

// Construct the weights of a four node cubic Lagrange stencil
std::array<double, 4> CubicLagrangeWeights(const std::array<double, 4> &node, double x);

// Construct one cubic interpolation stencil on a uniform finite grid
UniformCubicStencil UniformCubicWeights(std::size_t count, double minimum, double maximum, double x);

// Interpolate one value with precomputed cubic Lagrange weights
// f(x) = sum_{i=0}^3 w_i(x) f_i
template <typename T>
inline T CubicLagrangeWeightedSum(const std::array<T, 4> &value, const std::array<double, 4> &weight) {
  T output = value[0] * weight[0];
  for (std::size_t i = 1; i < value.size(); ++i) {
    if constexpr (requires { output.AddScaled(value[i], weight[i]); }) {
      output.AddScaled(value[i], weight[i]);
    } else {
      output += value[i] * weight[i];
    }
  }
  return output;
}

// Interpolate one value with a four node cubic Lagrange stencil
// f(x) = sum_i f_i product_{j != i}(x-x_j)/(x_i-x_j)
template <typename T>
inline T CubicLagrangeInterpolate(const std::array<double, 4> &node, const std::array<T, 4> &value, const double x) {
  return CubicLagrangeWeightedSum(value, CubicLagrangeWeights(node, x));
}

}  // namespace gra::math

#endif
