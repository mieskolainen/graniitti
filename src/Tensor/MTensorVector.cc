// Dispersive rho and omega propagation into charged pions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorVector.h"

#include <cmath>
#include <limits>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Tech/MException.h"

namespace gra::tensor {

// Compute the subtracted P-wave dispersion integral and its absorptive part
// [REFERENCE: Bolz et al., arXiv:1409.8483, Eqs. (B.10)-(B.14)]
std::complex<double> VectorLoop(const double s, const double mass) {
  const double x = s / (4.0 * math::pow2(mass));
  double real = 0.0;
  double imag = 0.0;
  if (std::abs(x) < 0.25) {
    // Expand x * 2F1(1,1;7/2;x)/(480*pi^2) without cancelling inverse powers of x
    double term = 1.0;
    double sum = term;
    std::size_t n = 1;
    do {
      term *= x * n / (n + 2.5);
      sum += term;
      ++n;
    } while (std::abs(term) > std::numeric_limits<double>::epsilon() * std::abs(sum));
    real = x * sum / (480.0 * math::pow2(math::PI));
  } else if (x > 1.0 || x < 0.0) {
    const double xi = std::sqrt(1.0 - 1.0 / x);
    real = (1.0 / 3.0 + math::pow2(xi) - 0.5 * math::pow3(xi) *
            (std::log(std::abs(x)) + 2.0 * std::log1p(xi))) / (96.0 * math::pow2(math::PI));
    if (x > 1.0) { imag = s * math::pow3(xi) / (192.0 * math::PI); }
  } else {
    const double xi = std::sqrt(1.0 / x - 1.0);
    real = (1.0 / 3.0 - math::pow2(xi) + math::pow3(xi) * std::atan2(1.0, xi)) /
           (96.0 * math::pow2(math::PI));
  }
  return {s * real, imag};
}

// Compute the coupled inverse with the pion and kaon cuts and the omega total width
// [REFERENCE: Bolz et al., arXiv:1409.8483, Eqs. (B.6)-(B.17)]
MMatrix<std::complex<double>> RhoOmega::Inverse(const double s) const {
  const double rho2 = math::pow2(mass[0]);
  const double omega2 = math::pow2(mass[1]);
  const auto pi = VectorLoop(s, pion);
  const auto kk = VectorLoop(s, kaon);
  const auto pi_rho = pi - s / rho2 * VectorLoop(rho2, pion).real();
  const auto kk_rho = kk - s / rho2 * VectorLoop(rho2, kaon).real();
  const double kk_omega = kk.real() - s / omega2 * VectorLoop(omega2, kaon).real();
  MMatrix<std::complex<double>> out(2, 2, 0.0);
  out[0][0] = s - rho2 + math::pow2(g[0]) * (pi_rho + 0.5 * kk_rho);
  out[1][1] = s - omega2 + math::pow2(g[0] / 2.0) * kk_omega + math::zi * s / mass[1] * width;
  out[0][1] = out[1][0] = s * b + g[0] * g[1] * pi_rho;
  return out;
}

// Translate a failed spectral solve into the event amplitude bookkeeping
MMatrix<std::complex<double>> RhoOmega::Propagator(const double s) const {
  try {
    return Inverse(s).Solve(MMatrix<std::complex<double>>::IdentityMatrix(2));
  } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("RhoOmega::Propagator: ") + error.what());
  }
}

}  // namespace gra::tensor
