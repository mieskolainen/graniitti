// Polar quadrature with circular cuts and matching scales
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MPolarCut.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra::math {

// Validate controls and precompute unit interval rules for each smooth segment
PolarCutRule::PolarCutRule(const PolarParam& input) : param(input) {
  (void)PolarRadialRule(param);
  if (param.azimuth_integrator != "GL" && param.azimuth_integrator != "Trap") {
    throw std::invalid_argument("PolarCutRule: unknown azimuth integrator");
  }
  auto unit       = param;
  unit.r_min      = 0.0;
  unit.r_max      = 1.0;
  unit.radial_map = RadialMap::Linear;
  radial          = PolarRadialRule(unit);
  const auto rule = param.azimuth_integrator == "GL" ? GaussLegendreRule(param.azimuth_nodes, 0.0, 1.0)
                                                     : PeriodicTrapzRule(param.azimuth_nodes, 0.0, 1.0);
  azimuth.node    = rule.first;
  azimuth.weight  = rule.second;
}

namespace {

// Split radius where the allowed azimuth changes at circular tangencies and intersections
std::vector<double> PolarRadii(const PolarParam& param, const std::vector<std::array<double, 2>>& centres,
                               double radius, const std::vector<std::array<double, 3>>& splits) {
  std::vector<double> bounds = {param.r_min, param.r_max};
  const auto          add    = [&](double r) {
    if (r > param.r_min && r < param.r_max) { bounds.push_back(r); }
  };
  for (const auto& i : indices(centres)) {
    const auto&  a = centres[i];
    const double d = std::hypot(a[0], a[1]);
    add(std::abs(d - radius));
    add(d + radius);
    add(d / 2.0);
    for (std::size_t j = 0; j < i; ++j) {
      const auto&  b        = centres[j];
      const double dx       = b[0] - a[0];
      const double dy       = b[1] - a[1];
      const double distance = std::hypot(dx, dy);
      if (!(distance > 0.0) || !(distance < 2.0 * radius)) { continue; }
      const double span = std::sqrt(radius * radius - distance * distance / 4.0) / distance;
      for (const double sign : {-1.0, 1.0}) {
        add(std::hypot((a[0] + b[0]) / 2.0 + sign * span * dy, (a[1] + b[1]) / 2.0 - sign * span * dx));
      }
    }
  }
  for (const auto& circle : splits) {
    const double d = std::hypot(circle[0], circle[1]);
    add(std::abs(d - circle[2]));
    add(d + circle[2]);
  }
  std::sort(bounds.begin(), bounds.end());
  return bounds;
}

// Split azimuth at each virtuality cutoff and change of the nearest gluon scale
std::vector<double> PolarAngles(const std::vector<std::array<double, 2>>& centres, double radius, double r,
                                double phi0, const std::vector<std::array<double, 3>>& splits) {
  std::vector<double> angles = {0.0, PI / 2.0, PI, 3.0 * PI / 2.0, 2.0 * PI};
  const auto          add    = [&](double angle) { angles.push_back(std::fmod(angle - phi0 + 4.0 * PI, 2.0 * PI)); };
  for (const auto& centre : centres) {
    const double d = std::hypot(centre[0], centre[1]);
    if (!(d > 0.0)) { continue; }
    const double axis = std::atan2(centre[1], centre[0]);
    for (const double cosine : {(r * r + d * d - radius * radius) / (2.0 * r * d), d / (2.0 * r)}) {
      if (!(std::abs(cosine) < 1.0)) { continue; }
      const double angle = std::acos(cosine);
      add(axis - angle);
      add(axis + angle);
    }
  }
  for (const auto& circle : splits) {
    const double d = std::hypot(circle[0], circle[1]);
    if (!(d > 0.0)) { continue; }
    const double cosine = (r * r + d * d - circle[2] * circle[2]) / (2.0 * r * d);
    if (!(std::abs(cosine) < 1.0)) { continue; }
    const double axis  = std::atan2(circle[1], circle[0]);
    const double angle = std::acos(cosine);
    add(axis - angle);
    add(axis + angle);
  }
  std::sort(angles.begin(), angles.end());
  return angles;
}

}  // namespace

// Integrate common radial rings and smooth angular intervals outside the cut disks
std::vector<PolarPoint> PolarCutRule::Nodes(const std::vector<std::array<double, 2>>& centres, double radius,
                                            double phi0, const std::vector<std::array<double, 3>>& splits) const {
  phi0 = WrapAngle(phi0);
  const auto              bounds      = PolarRadii(param, centres, radius, splits);
  const bool              logarithmic = param.radial_map == RadialMap::Log;
  std::vector<PolarPoint> nodes;
  nodes.reserve(bounds.size() * radial.node.size() * azimuth.node.size() * 4);
  unsigned int ring = 0;
  for (std::size_t b = 1; b < bounds.size(); ++b) {
    const double lower = bounds[b - 1];
    const double upper = bounds[b];
    if (!(upper > lower)) { continue; }
    const double width = logarithmic ? std::log(upper / lower) : upper - lower;
    const auto&  rule  = radial;
    for (const auto& i : indices(rule.node)) {
      double u   = rule.node[i];
      double jac = 1.0;
      // Remove square-root tangencies at interior radial boundaries
      if (b > 1 && b + 1 < bounds.size()) {
        jac = PI / 2.0 * std::sin(PI * u);
        u   = (1.0 - std::cos(PI * u)) / 2.0;
      } else if (b > 1) {
        jac = 2.0 * u;
        u *= u;
      } else if (b + 1 < bounds.size()) {
        jac = 2.0 * (1.0 - u);
        u   = 1.0 - pow2(1.0 - u);
      }
      double r = lower + width * u;
      if (logarithmic) {
        r = lower * std::exp(width * u);
        jac *= r;
      } else if (param.radial_map == RadialMap::Square) {
        r = lower + width * u * u;
        jac *= 2.0 * u;
      }
      jac *= r * width * rule.weight[i];
      ++ring;
      if (!(r > 0.0) || !(jac > 0.0)) { continue; }
      const auto angles = PolarAngles(centres, radius, r, phi0, splits);
      for (std::size_t a = 1; a < angles.size(); ++a) {
        const double span   = angles[a] - angles[a - 1];
        const double middle = phi0 + (angles[a] + angles[a - 1]) / 2.0;
        if (!(span > 0.0) || std::any_of(centres.begin(), centres.end(), [&](const auto& centre) {
              return pow2(r * std::cos(middle) - centre[0]) + pow2(r * std::sin(middle) - centre[1]) < pow2(radius);
            })) {
          continue;
        }
        for (const auto& j : indices(azimuth.node)) {
          const double phi = phi0 + angles[a - 1] + span * azimuth.node[j];
          nodes.push_back({r * std::cos(phi), r * std::sin(phi), jac * span * azimuth.weight[j], ring});
        }
      }
    }
  }
  return nodes;
}

}  // namespace gra::math
