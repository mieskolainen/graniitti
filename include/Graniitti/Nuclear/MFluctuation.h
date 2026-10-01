// Positive cross-section fluctuation quadrature
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARFLUCTUATION_H
#define MNUCLEARFLUCTUATION_H

#include <vector>

namespace gra::nuclear {

// Store cross-section fluctuation physics and numerical controls
struct FluctuationParam {
  double       omega            = 0.0;  // Relative cross-section variance
  unsigned int nodes            = 0;    // Equal-weight eigenstate nodes
  double       cdf_tol          = 0.0;  // Gamma CDF relative tolerance
  unsigned int cdf_max_iter     = 0;    // Gamma CDF iteration limit
  double       normal_shape_min = 0.0;  // Normal-quantile Gamma shape boundary
};

// Construct positive equal-weight eigenvalues with exact first two moments
std::vector<double> CrossSectionRule(const FluctuationParam &param);

}  // namespace gra::nuclear

#endif
