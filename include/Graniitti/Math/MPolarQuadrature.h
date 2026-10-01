// Polar tensor-product quadrature
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPOLARQUADRATURE_H
#define MPOLARQUADRATURE_H

#include <string>
#include <vector>

namespace gra::math {

// Select the radial coordinate transformation
enum class RadialMap { Linear, Log, Square };

// Store one polar quadrature definition
struct PolarParam {
  std::string radial_integrator;
  std::string azimuth_integrator;
  RadialMap radial_map = RadialMap::Linear;
  double r_min = 0.0;
  double r_max = 0.0;
  unsigned int radial_intervals = 0;
  unsigned int azimuth_nodes = 0;
};

// Store one transformed one-dimensional quadrature rule
struct PolarRule1D {
  std::vector<double> node;
  std::vector<double> weight;
  std::vector<double> jac;
  std::vector<double> measure;
  double step = 0.0;
  double map_min = 0.0;
};

// Compute the interval multiple required by one radial rule
unsigned int PolarIntervalMultiple(const PolarParam &param);

// Compute the number of radial nodes required by one polar rule
unsigned int PolarNodeCount(const PolarParam &param);

// Compute the closed Newton-Cotes coefficient of one radial rule
double PolarCoefficient(const PolarParam &param);

// Construct transformed radial nodes, weights and measure Jacobians
PolarRule1D PolarRadialRule(const PolarParam &param);

// Construct periodic azimuthal nodes and weights
PolarRule1D PolarAzimuthRule(const PolarParam &param);

} // namespace gra::math

#endif
