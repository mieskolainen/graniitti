// Fourier evolution of the periodic Kuramoto Sivashinsky equation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_PDEBENCH_H
#define PROGRAM_PDEBENCH_H

#include <Eigen/unsupported/Eigen/FFT>
#include <complex>
#include <limits>
#include <stdexcept>
#include <valarray>
#include <vector>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MFFT.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::program {

// Validate the Fourier lattice and time evolution before allocating arrays
inline void ValidateKS(int length, int points, double dt, int steps, int nplot, bool eigen) {
  if (length <= 0 || points < 4 || points % 2 != 0 || (!eigen && (points & (points - 1)) != 0) || !std::isfinite(dt) ||
      dt <= 0.0 || steps < 0 || nplot <= 0) {
    throw std::invalid_argument("pdebench: invalid length, Fourier grid, time step or plotting interval");
  }
}

class KS {
 public:
  // Initialize the spectral field and nonlinear history for CNAB2
  // [REFERENCE: https://online.kitp.ucsb.edu/online/transturb17/gibson/html/5-kuramoto-sivashinksy.html]
  KS(const std::valarray<std::complex<double>>& field, int length, double dt, bool eigen)
      : eigen_(eigen),
        dt_(dt),
        u_(field),
        A_(field.size()),
        B_(field.size()),
        G_(field.size()),
        previous_(field.size()) {
    if (field.size() > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("pdebench: Fourier grid exceeds integer range");
    }
    const int n = static_cast<int>(field.size());
    ValidateKS(length, n, dt, 0, 1, eigen);
    for (const auto& i : gra::aux::indices(field)) {
      if (!std::isfinite(field[i].real()) || !std::isfinite(field[i].imag()) || std::abs(field[i].imag()) > 1e-12) {
        throw std::invalid_argument("pdebench: initial field must be finite and real");
      }
      const int    j           = static_cast<int>(i);
      const int    mode        = j < n / 2 ? j : j - n;
      const double k           = 2.0 * math::PI * mode / length;
      const double linear      = k * k - math::pow4(k);
      const double denominator = 1.0 - 0.5 * dt * linear;
      if (std::abs(denominator) < 1e-12) { throw std::invalid_argument("pdebench: singular Crank Nicolson step"); }
      A_[i] = 1.0 + 0.5 * dt * linear;
      B_[i] = 1.0 / denominator;
      G_[i] = j == n / 2 ? 0.0 : std::complex<double>(0.0, -0.5 * k);
    }
    previous_ = field * field;
    Transform(previous_, false);
    previous_ *= G_;
    Transform(u_, false);
  }

  // Advance the field using the current nonlinear term and its previous value
  void Step() {
    auto nonlinear = Field();
    for (const auto& i : gra::aux::indices(nonlinear)) { nonlinear[i] = math::pow2(nonlinear[i].real()); }
    Transform(nonlinear, false);
    nonlinear *= G_;
    u_        = B_ * (A_ * u_ + 1.5 * dt_ * nonlinear - 0.5 * dt_ * previous_);
    previous_ = nonlinear;
    for (const auto& i : gra::aux::indices(u_)) {
      if (!std::isfinite(u_[i].real()) || !std::isfinite(u_[i].imag())) {
        throw std::runtime_error("pdebench: non-finite field, reduce the time step");
      }
    }
  }

  // Compute the field in coordinate space
  std::valarray<std::complex<double>> Field() const {
    auto field = u_;
    Transform(field, true);
    return field;
  }

 private:
  // Apply the selected FFT with identical forward and inverse normalization
  void Transform(std::valarray<std::complex<double>>& field, bool inverse) const {
    if (eigen_) {
      Eigen::FFT<double>                fft;
      const auto                        input = gra::valarray2vector(field);
      std::vector<std::complex<double>> output;
      if (inverse) {
        fft.inv(output, input);
      } else {
        fft.fwd(output, input);
      }
      field = gra::vector2valarray(output);
    } else if (inverse) {
      MFFT::ifft(field);
    } else {
      MFFT::fft(field);
    }
  }

  bool                                eigen_;
  double                              dt_;
  std::valarray<std::complex<double>> u_, A_, B_, G_, previous_;
};

}  // namespace gra::program
#endif
