// Shared form factor parameter types
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFORMFACTOR_H
#define MFORMFACTOR_H

#include <cmath>
#include <span>
#include <string>
#include <vector>

#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Tech/MException.h"

namespace gra::regge {

enum class FFType {
  None,
  Exponential,
  Power,
  Orear,
  Gaussian,
  Vector,
  GKernel,
  Dirac,
  LogExp
};

enum class FFNorm { Zero, Pole };

// Store one normalized form factor and its family parameters
struct FFParam {
  FFType type = FFType::None;
  FFNorm norm = FFNorm::Zero;
  std::vector<double> param;
  form::ParamStore structure;
};

// Evaluate the curvature matched log squared profile with validated b and Lambda2
// F = exp[-b*Lambda2*(y + y^2/2)], y = log(1+x/Lambda2), intended for x >= 0
// x and Lambda2 are in GeV^2, b is in GeV^-2
inline double LogExpFF(double x, double b, double Lambda2) {
  if (x <= -Lambda2) {
    throw AmplitudeFailure("LogExpFF: real domain requires x > -Lambda2");
  }
  if (std::fpclassify(b) == FP_ZERO || std::fpclassify(x) == FP_ZERO) { return 1.0; }
  const double z = x / Lambda2;
  // Avoid overflow in x/Lambda2 while retaining the initial derivative near zero
  const double y = std::isfinite(z) ? std::log1p(z)
                                   : std::log(x) - std::log(Lambda2);
  if (std::fpclassify(y) == FP_ZERO) { return std::exp(-b * x); }
  const double phi = b * (Lambda2 * (y * (1.0 + 0.5 * y)));
  return std::exp(-phi);
}

// Evaluate generalized attenuation factors from validated [a,p,nu,mu2] rows
inline double GKernel(double x, std::span<const double> params, const std::string& context, double product = 1.0) {
  for (std::size_t i = 0; i < params.size(); i += 4) {
    const double a = params[i], p = params[i + 1], nu = params[i + 2], mu2 = params[i + 3];
    if (x + mu2 < 0.0) { throw AmplitudeFailure(context + ": negative kernel domain"); }
    const double y = a * (std::pow(x + mu2, p) - std::pow(mu2, p));
    if (!std::isfinite(y)) { throw AmplitudeFailure(context + ": non-finite attenuation"); }
    if (nu <= 1.0e-12) {
      product *= std::exp(-y);
    } else {
      const double base = 1.0 + nu * y;
      if (!std::isfinite(base) || base <= 0.0) { throw AmplitudeFailure(context + ": nonpositive kernel base"); }
      product *= std::exp(-std::log1p(nu * y) / nu);
    }
  }
  if (!std::isfinite(product)) { throw AmplitudeFailure(context + ": non-finite kernel product"); }
  return product;
}

} // namespace gra::regge

#endif
