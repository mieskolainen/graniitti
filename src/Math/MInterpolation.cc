// Interpolation grid utilities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MInterpolation.h"

#include <algorithm>
#include <cmath>
#include <iterator>
#include <stdexcept>

namespace gra::math {

// Validate one finite strictly increasing interpolation grid
void ValidateInterpolationGrid(const std::vector<double> &node) {
  if (node.empty() ||
      !std::all_of(node.cbegin(), node.cend(), [](const double value) { return std::isfinite(value); })) {
    throw std::invalid_argument("math::LinearInterpolate: grid nodes must be finite and nonempty");
  }
  for (std::size_t i = 1; i < node.size(); ++i) {
    if (!(node[i] > node[i - 1])) {
      throw std::invalid_argument("math::LinearInterpolate: grid nodes must be strictly increasing");
    }
  }
}

// Locate one value on a finite strictly increasing interpolation grid
// t = (x-x_i)/(x_{i+1}-x_i), clipped to the first or last cell
InterpolationCell LocateCell(const std::vector<double> &node, const double value) {
  if (!std::isfinite(value)) { throw std::invalid_argument("math::LocateCell: invalid grid or value"); }
  ValidateInterpolationGrid(node);
  if (node.size() < 2) { throw std::invalid_argument("math::LocateCell: grid has one node"); }
  if (value <= node.front()) { return {0, 0.0}; }
  if (value >= node.back()) { return {node.size() - 2, 1.0}; }
  const auto        upper = std::upper_bound(node.cbegin(), node.cend(), value);
  const std::size_t index = static_cast<std::size_t>(std::distance(node.cbegin(), upper) - 1);
  return {index, (value - node[index]) / (node[index + 1] - node[index])};
}

// Construct the weights of a four node cubic Lagrange stencil
// w_i(x) = product_{j != i}(x-x_j)/(x_i-x_j)
std::array<double, 4> CubicLagrangeWeights(const std::array<double, 4> &node, const double x) {
  if (!std::isfinite(x) ||
      !std::all_of(node.cbegin(), node.cend(), [](const double coordinate) { return std::isfinite(coordinate); })) {
    throw std::invalid_argument("math::CubicLagrangeWeights: non-finite coordinate");
  }
  std::array<double, 4> weight{};
  for (std::size_t i = 0; i < node.size(); ++i) {
    weight[i] = 1.0;
    for (std::size_t j = 0; j < node.size(); ++j) {
      if (i == j) { continue; }
      const double separation = node[i] - node[j];
      if (math::IsZero(separation)) { throw std::invalid_argument("math::CubicLagrangeWeights: repeated grid node"); }
      weight[i] *= (x - node[j]) / separation;
    }
  }
  return weight;
}

// Construct one cubic interpolation stencil on a uniform finite grid
// x_i = x_min + i(x_max-x_min)/(N-1) with four adjacent Lagrange nodes
UniformCubicStencil UniformCubicWeights(const std::size_t count, const double minimum, const double maximum,
                                        const double x) {
  if (count < 4 || !std::isfinite(minimum) || !std::isfinite(maximum) || !(maximum > minimum) || !std::isfinite(x) ||
      x < minimum || x > maximum) {
    throw std::invalid_argument("math::UniformCubicWeights: invalid grid");
  }
  const double          step     = (maximum - minimum) / static_cast<double>(count - 1);
  const double          position = (x - minimum) / step;
  const std::size_t     cell     = std::min(static_cast<std::size_t>(position), count - 2);
  const std::size_t     first    = std::min(cell > 0 ? cell - 1 : 0, count - 4);
  UniformCubicStencil   stencil;
  std::array<double, 4> node{};
  for (std::size_t i = 0; i < node.size(); ++i) {
    stencil.index[i] = first + i;
    node[i]          = minimum + static_cast<double>(first + i) * step;
  }
  stencil.weight = CubicLagrangeWeights(node, x);
  return stencil;
}

}  // namespace gra::math
