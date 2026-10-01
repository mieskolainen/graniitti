// Spherical waves and the complete Rayleigh Sommerfeld propagator
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_SOMMERFELD_H
#define PROGRAM_SOMMERFELD_H

#include <complex>
#include <stdexcept>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMath.h"

namespace gra::program {

// Compute the spherical source field with the exp(+i omega t) convention
inline std::complex<double> u_point(const M4Vec& x, const M4Vec& source, double k) {
  const double r = (x - source).P3mod();
  if (!(r > 0.0) || !std::isfinite(r) || !std::isfinite(k) || k <= 0.0) {
    throw std::invalid_argument("sommerfeld: invalid source distance or wavenumber");
  }
  return std::exp(math::zi * k * (x.T() - r)) / r;
}

// Compute the full normal derivative of exp(-ikR)/R on the aperture
// [REFERENCE: https://doi.org/10.1364/JOSAA.21.000510]
inline std::complex<double> RS_integrand(const M4Vec& x, const M4Vec& detector, const M4Vec& source,
                                         const M4Vec& normal, double k) {
  const M4Vec  delta = x - detector;
  const double r     = delta.P3mod();
  if (!(r > 0.0) || !std::isfinite(r)) {
    throw std::invalid_argument("sommerfeld: detector must lie off the aperture");
  }
  return (-math::zi * k - 1.0 / r) / (2.0 * math::PI) * u_point(x, source, k) * std::exp(-math::zi * k * r) / r *
         normal.Dot3(delta) / r;
}

}  // namespace gra::program
#endif
