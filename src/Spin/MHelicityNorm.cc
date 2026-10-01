// Forward helicity normalization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Spin/MHelicityNorm.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra::spin {

using gra::aux::indices;

namespace {

// Compute one integer twice-helicity grade
std::size_t HelicityGrade(double helicity, ForwardLegType leg, double tolerance) {
  if (!std::isfinite(helicity)) { throw std::invalid_argument("ForwardHelicityDensity: non-finite helicity"); }
  if (leg == ForwardLegType::RealPhoton) {
    if (std::abs(std::abs(helicity) - 1.0) > tolerance) {
      throw std::invalid_argument(
          "ForwardHelicityDensity: photon helicity is not "
          "transverse");
    }
    return 0;
  }
  const double absolute       = std::abs(helicity);
  const double twice_helicity = 2.0 * absolute;
  const double rounded        = std::round(twice_helicity);
  if (std::abs(twice_helicity - rounded) > tolerance || rounded < 0.0) {
    throw std::invalid_argument(
        "ForwardHelicityDensity: hadronic helicity grade is "
        "not integer or half-integer");
  }
  return static_cast<std::size_t>(rounded);
}

// Compute the physical incoming-helicity averaging weight
double IncomingWeight(ForwardLegType upper, ForwardLegType lower) {
  return (upper == ForwardLegType::RealPhoton ? 0.5 : 1.0) * (lower == ForwardLegType::RealPhoton ? 0.5 : 1.0);
}

}  // namespace

// Compute the physical reduced-helicity norm without basis padding
// ||H||^2 = sum_(physical lambda1,lambda2) |H_lambda1,lambda2|^2
double ReducedHelicityNorm2(const HELMatrix &helicity) {
  double norm2 = 0.0;
  for (std::size_t row = 0; row < helicity.lambda_values.size_row(); ++row) {
    if (helicity.lambda_idx.size_col() != 2 || row >= helicity.lambda_idx.size_row()) {
      throw std::invalid_argument("ReducedHelicityNorm2: invalid helicity metadata");
    }
    const std::size_t i1 = helicity.lambda_idx[row][0];
    const std::size_t i2 = helicity.lambda_idx[row][1];
    if (i1 >= helicity.T.size_row() || i2 >= helicity.T.size_col()) {
      throw std::invalid_argument("ReducedHelicityNorm2: helicity index outside T");
    }
    norm2 += gra::math::abs2(helicity.T[i1][i2]);
  }
  if (!(norm2 > 0.0) || !std::isfinite(norm2)) {
    throw std::invalid_argument("ReducedHelicityNorm2: physical tensor norm is invalid");
  }
  return norm2;
}

// Compute the leading forward density of one coherent central-helicity tensor
// rho = sum_leading |M|^2 averaged over physical incoming helicities
double ForwardHelicityDensity(const MMatrix<std::complex<double>> &tensor, const MMatrix<double> &incoming_helicities,
                              ForwardLegType upper, ForwardLegType lower, double tolerance) {
  if (!(tolerance > 0.0) || !std::isfinite(tolerance)) {
    throw std::invalid_argument("ForwardHelicityDensity: invalid tolerance");
  }
  if (tensor.isEmpty() || tensor.size_row() != incoming_helicities.size_row() || incoming_helicities.size_col() != 2) {
    throw std::invalid_argument("ForwardHelicityDensity: invalid tensor metadata");
  }
  if (!tensor.IsFinite()) { throw AmplitudeFailure("ForwardHelicityDensity: non-finite amplitude"); }

  std::size_t  leading_grade    = std::numeric_limits<std::size_t>::max();
  double       weighted_density = 0.0;
  double       total_weight     = 0.0;
  const double incoming_weight  = IncomingWeight(upper, lower);
  for (std::size_t row = 0; row < tensor.size_row(); ++row) {
    const double row_density = gra::SquaredNorm(tensor.Row(row));
    if (row_density <= tolerance * tolerance) { continue; }
    const std::size_t grade = HelicityGrade(incoming_helicities[row][0], upper, tolerance) +
                              HelicityGrade(incoming_helicities[row][1], lower, tolerance);
    if (grade > leading_grade) { continue; }
    if (grade < leading_grade) {
      leading_grade    = grade;
      weighted_density = total_weight = 0.0;
    }
    weighted_density += incoming_weight * row_density;
    total_weight += incoming_weight;
  }
  // A threshold zero or coherent LS cancellation has zero physical density
  if (leading_grade == std::numeric_limits<std::size_t>::max()) { return 0.0; }
  if (!(weighted_density > 0.0) || !(total_weight > 0.0) || !std::isfinite(weighted_density) ||
      !std::isfinite(total_weight)) {
    throw AmplitudeFailure("ForwardHelicityDensity: leading density is invalid");
  }

  return weighted_density / total_weight;
}

}  // namespace gra::spin
