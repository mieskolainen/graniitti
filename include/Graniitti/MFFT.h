// Simple in-place Fast Fourier Transform with radix-2 using std::valarray
//
// Normalization is 1 to 1, that is: ifft( fft(x) ) = x
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFFT_H
#define MFFT_H

// C++
#include <complex>
#include <valarray>

namespace gra {
#ifndef M__PI
#define M__PI 3.141592653589793238462643383279502884197169399375105820974944L
#endif

namespace MFFT {
// In-place Cooley–Tukey FFT
template <typename T>
void fft(std::valarray<std::complex<T>> &x) {
  const std::size_t N = x.size();
  if (N <= 1) {
    return;
  }  // Trivial case/recursion ends, X = x
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
    const std::complex<T> t = std::exp(std::complex<T>(0, -2.0 * M__PI * k / N));
    x[k] = E[k] + O[k] * t;
    x[k + N / 2] = E[k] - O[k] * t;
  }
}

// In-place Cooley-Tukey IFFT
template <typename T>
void ifft(std::valarray<std::complex<T>> &x) {
  x = x.apply(std::conj);  // Conjugate
  MFFT::fft(x);            // FFT
  x = x.apply(std::conj);  // Conjugate
  x /= x.size();           // Normalize
}
}

}  // gra namespace ends

#endif