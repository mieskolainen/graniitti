// Polar quadrature with circular cuts and matching scales
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPOLARCUT_H
#define MPOLARCUT_H

#include <array>
#include "Graniitti/Math/MPolarQuadrature.h"

namespace gra::math {

// Store one transverse node with its complete d^2q weight
struct PolarPoint {
  double       x      = 0.0;
  double       y      = 0.0;
  double       weight = 0.0;
  unsigned int radial = 0;
};

// Store unit interval rules for integration between circular boundaries
struct PolarCutRule {
  PolarParam  param;
  PolarRule1D radial;
  PolarRule1D azimuth;

  // Validate controls and precompute the unit interval rules
  explicit PolarCutRule(const PolarParam& input);

  // Construct nodes outside equal cut disks, with optional circle splits {x, y, r}
  std::vector<PolarPoint> Nodes(const std::vector<std::array<double, 2>>& centres, double radius,
                                double phi0 = 0.0, const std::vector<std::array<double, 3>>& splits = {}) const;
};

}  // namespace gra::math

#endif
