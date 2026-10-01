// Simple in-place Fast Fourier Transform with radix-2 using std::valarray
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Normalization is 1 to 1, that is: ifft( fft(x) ) = x

#ifndef MFFT_H
#define MFFT_H

// C++
#include <complex>
#include <stdexcept>
#include <string>
#include <valarray>

#include "Graniitti/Math/MMath.h"

namespace gra {

namespace MFFT {
// In-place Cooley-Tukey FFT
// X_k = sum_{n=0}^{N-1} x_n exp(-2 pi i kn/N)
template <typename T>
void fft(std::valarray<std::complex<T>> &x) {
  const std::size_t N = x.size();
  if (N <= 1) { return; }  // Trivial case/recursion ends, X = x
  if ((N & (N - 1)) != 0) {
    throw std::invalid_argument("ERROR: MFFT::fft: Input x.size() = " + std::to_string(N) +
                                " not a power of 2!");
  }

  // Radix-2 step
  std::valarray<std::complex<T>> E = x[std::slice(0, N / 2, 2)];
  std::valarray<std::complex<T>> O = x[std::slice(1, N / 2, 2)];

  // Even and odd part via recursion
  MFFT::fft(E);
  MFFT::fft(O);

  for (std::size_t k = 0; k < N / 2; ++k) {
    const std::complex<T> t = std::exp(std::complex<T>(0, -2.0 * math::PI * k / N));
    x[k]                    = E[k] + O[k] * t;
    x[k + N / 2]            = E[k] - O[k] * t;
  }
}

// In-place Cooley-Tukey IFFT
// x_n = N^{-1} sum_{k=0}^{N-1} X_k exp(+2 pi i kn/N)
template <typename T>
void ifft(std::valarray<std::complex<T>> &x) {
  const std::size_t N = x.size();
  if (N <= 1) { return; }
  if ((N & (N - 1)) != 0) {
    throw std::invalid_argument("ERROR: MFFT::ifft: Input x.size() = " + std::to_string(N) + " not a power of 2!");
  }
  x = x.apply(std::conj);  // Conjugate
  MFFT::fft(x);            // FFT
  x = x.apply(std::conj);  // Conjugate
  x /= N;                  // Normalize
}
}  // namespace MFFT

}  // namespace gra

#endif
