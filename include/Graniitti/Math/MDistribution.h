// Generic probability distribution functions and moment matching
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATHDISTRIBUTION_H
#define MMATHDISTRIBUTION_H

namespace gra::math {

// Store iterative controls for Gamma distribution functions
struct GammaNumerics {
  // Validate immutable iteration controls before evaluating a distribution
  GammaNumerics(double tolerance, unsigned int max_iter, double normal_min);
  const double cdf_tol;
  const unsigned int cdf_max_iter;
  const double normal_shape_min;
};

// Compute the standard normal quantile
double NormalQuantile(double probability);

// Compute the regularized lower incomplete Gamma function
double GammaCDF(double shape, double value, const GammaNumerics &numerics);

// Compute one Gamma quantile with unit scale
double GammaQuantile(double probability, double shape,
                     const GammaNumerics &numerics);

} // namespace gra::math

#endif
