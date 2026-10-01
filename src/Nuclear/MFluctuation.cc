// Positive cross-section fluctuation quadrature
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MFluctuation.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MDistribution.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::nuclear {

using gra::aux::indices;

// Construct positive equal-weight eigenvalues with exact first two moments
// <s> = 1 and <s^2> - <s>^2 = omega
// [REFERENCE: Alvioli and Strikman, Phys. Lett. B722 (2013) 347]
std::vector<double> CrossSectionRule(const FluctuationParam &param) {
  const math::GammaNumerics numerics{param.cdf_tol, param.cdf_max_iter, param.normal_shape_min};
  if (!std::isfinite(param.omega) || param.omega < 0.0 || param.nodes == 0) {
    throw std::invalid_argument("CrossSectionRule: invalid controls");
  }
  if (std::fpclassify(param.omega) == FP_ZERO) { return {1.0}; }
  if (param.nodes < 2 || !(param.omega < static_cast<double>(param.nodes - 1))) {
    throw std::invalid_argument("CrossSectionRule: nodes cannot resolve the requested omega");
  }

  std::vector<double> log_node(param.nodes, 0.0);
  const double        shape = 1.0 / param.omega;
  for (const auto &i : indices(log_node)) {
    const double probability = (static_cast<double>(i) + 0.5) / static_cast<double>(param.nodes);
    const double value       = param.omega * math::GammaQuantile(probability, shape, numerics);
    if (!std::isfinite(value) || !(value > 0.0)) { throw std::runtime_error("CrossSectionRule: invalid eigenvalue"); }
    log_node[i] = std::log(value);
  }

  double lower = 0.0;
  double upper = 1.0;
  while (statistics::PoweredVariance(log_node, upper) < param.omega) {
    upper *= 2.0;
    if (!std::isfinite(upper)) { throw std::runtime_error("CrossSectionRule: moment match did not converge"); }
  }
  for (unsigned int iteration = 0; iteration < param.cdf_max_iter; ++iteration) {
    const double middle = 0.5 * (lower + upper);
    if (statistics::PoweredVariance(log_node, middle) < param.omega) {
      lower = middle;
    } else {
      upper = middle;
    }
    if (upper - lower <= param.cdf_tol * std::max(1.0, std::abs(upper))) { break; }
  }
  if (upper - lower > param.cdf_tol * std::max(1.0, std::abs(upper))) {
    throw std::runtime_error("CrossSectionRule: moment match did not converge");
  }

  const double        power   = 0.5 * (lower + upper);
  const double        maximum = power * *std::max_element(log_node.begin(), log_node.end());
  std::vector<double> scale(param.nodes, 0.0);
  for (const auto &i : indices(scale)) { scale[i] = std::exp(power * log_node[i] - maximum); }
  const double mean = gra::Sum(scale) / static_cast<double>(scale.size());
  for (double &value : scale) { value /= mean; }
  return scale;
}

}  // namespace gra::nuclear
