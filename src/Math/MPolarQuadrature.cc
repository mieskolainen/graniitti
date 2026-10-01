// Polar tensor-product quadrature
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MPolarQuadrature.h"

#include <cmath>
#include <stdexcept>
#include <utility>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra::math {

// Compute the interval multiple required by one radial rule
unsigned int PolarIntervalMultiple(const PolarParam &param) {
  if (param.radial_integrator == "GL") {
    return 1;
  }
  if (param.radial_integrator == "1/3") {
    return 2;
  }
  if (param.radial_integrator == "3/8") {
    return 3;
  }
  if (param.radial_integrator == "Boole") {
    return 4;
  }
  throw std::invalid_argument(
      "PolarIntervalMultiple: unknown radial integrator " +
      param.radial_integrator);
}

// Compute the number of radial nodes required by one polar rule
unsigned int PolarNodeCount(const PolarParam &param) {
  if (param.radial_integrator == "GL") {
    return param.radial_intervals;
  }
  PolarIntervalMultiple(param);
  return param.radial_intervals + 1;
}

// Compute the closed Newton-Cotes coefficient of one radial rule
// c = 1/3, 3/8, or 2/45 for Simpson 1/3, Simpson 3/8, or Boole
double PolarCoefficient(const PolarParam &param) {
  if (param.radial_integrator == "1/3") {
    return 1.0 / 3.0;
  }
  if (param.radial_integrator == "3/8") {
    return 3.0 / 8.0;
  }
  if (param.radial_integrator == "Boole") {
    return 2.0 / 45.0;
  }
  if (param.radial_integrator == "GL") {
    return 1.0;
  }
  throw std::invalid_argument("PolarCoefficient: unknown radial integrator " +
                              param.radial_integrator);
}

namespace {

// Fill the complete transformed radial measure
// measure_i = weight_i jac_i
void FillMeasure(PolarRule1D &rule) {
  rule.measure.resize(rule.node.size(), 0.0);
  for (const auto &i : indices(rule.node)) {
    rule.measure[i] = rule.weight[i] * rule.jac[i];
  }
}

} // namespace

// Construct transformed radial nodes, weights and measure Jacobians
// int r dr f(r) = int du J(u) f[r(u)] with J = r, r^2, or 2 Delta r u r
PolarRule1D PolarRadialRule(const PolarParam &param) {
  const unsigned int multiple = PolarIntervalMultiple(param);
  if (param.radial_intervals == 0 || !std::isfinite(param.r_min) ||
      !std::isfinite(param.r_max) || param.r_min < 0.0 ||
      !(param.r_max > param.r_min) ||
      param.radial_intervals % multiple != 0 ||
      (param.radial_map == RadialMap::Log && !(param.r_min > 0.0)) ||
      (param.radial_map == RadialMap::Square &&
       param.radial_integrator != "GL")) {
    throw std::invalid_argument("PolarRadialRule: invalid radial controls");
  }

  PolarRule1D rule;
  rule.node.assign(PolarNodeCount(param), 0.0);
  rule.weight.assign(rule.node.size(), 0.0);
  rule.jac.assign(rule.node.size(), 0.0);
  if (param.radial_map == RadialMap::Log) {
    rule.map_min = std::log(param.r_min);
    rule.step = (std::log(param.r_max) - rule.map_min) /
                param.radial_intervals;
  } else if (param.radial_map == RadialMap::Square) {
    rule.step = 1.0 / param.radial_intervals;
  } else {
    rule.step = (param.r_max - param.r_min) / param.radial_intervals;
  }

  if (param.radial_integrator == "GL") {
    std::pair<std::vector<double>, std::vector<double>> raw;
    if (param.radial_map == RadialMap::Log) {
      raw = GaussLegendreRule(param.radial_intervals, rule.map_min,
                              std::log(param.r_max));
    } else if (param.radial_map == RadialMap::Square) {
      raw = GaussLegendreRule(param.radial_intervals, 0.0, 1.0);
    } else {
      raw = GaussLegendreRule(param.radial_intervals, param.r_min,
                              param.r_max);
    }
    for (const auto &i : indices(rule.node)) {
      if (param.radial_map == RadialMap::Log) {
        rule.node[i] = std::exp(raw.first[i]);
        rule.weight[i] = raw.second[i];
        rule.jac[i] = pow2(rule.node[i]);
      } else if (param.radial_map == RadialMap::Square) {
        const double unit = raw.first[i];
        const double range = param.r_max - param.r_min;
        rule.node[i] = param.r_min + range * pow2(unit);
        rule.weight[i] = raw.second[i] * 2.0 * range * unit;
        rule.jac[i] = rule.node[i];
      } else {
        rule.node[i] = raw.first[i];
        rule.weight[i] = raw.second[i];
        rule.jac[i] = rule.node[i];
      }
    }
    FillMeasure(rule);
    return rule;
  }

  std::vector<double> raw_weight;
  if (param.radial_integrator == "1/3") {
    raw_weight = Simpson13Weight(param.radial_intervals);
  } else if (param.radial_integrator == "3/8") {
    raw_weight = Simpson38Weight(param.radial_intervals);
  } else if (param.radial_integrator == "Boole") {
    raw_weight = BooleWeight(param.radial_intervals);
  } else {
    throw std::invalid_argument("PolarRadialRule: unknown radial integrator " +
                                param.radial_integrator);
  }
  const double scale = rule.step * PolarCoefficient(param);
  for (const auto &i : indices(rule.node)) {
    const double raw = param.radial_map == RadialMap::Log
                           ? rule.map_min + i * rule.step
                           : param.r_min + i * rule.step;
    const double radius =
        param.radial_map == RadialMap::Log ? std::exp(raw) : raw;
    rule.node[i] = radius;
    rule.weight[i] = scale * raw_weight[i];
    rule.jac[i] =
        param.radial_map == RadialMap::Log ? pow2(radius) : radius;
  }
  FillMeasure(rule);
  return rule;
}

// Construct periodic azimuthal nodes and weights
// phi_j = 2pi(j+1/2)/N, w_j = 2pi/N
PolarRule1D PolarAzimuthRule(const PolarParam &param) {
  if (param.azimuth_nodes == 0 || param.azimuth_integrator != "Trap") {
    throw std::invalid_argument(
        "PolarAzimuthRule: invalid azimuthal controls");
  }
  auto raw = PeriodicTrapzRule(param.azimuth_nodes, 0.0, 2.0 * PI);
  PolarRule1D rule;
  rule.step = raw.second.front();
  rule.node = std::move(raw.first);
  rule.weight = std::move(raw.second);
  rule.jac.assign(rule.node.size(), 1.0);
  rule.measure = rule.weight;
  return rule;
}

} // namespace gra::math
