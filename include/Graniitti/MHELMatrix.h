// Helicity coupling structures [HEADER ONLY file]
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELMATRIX_H
#define MHELMATRIX_H

// C++
#include <complex>
#include <random>
#include <valarray>
#include <vector>

#include "Graniitti/MMatrix.h"

namespace gra {

// Decay helicity amplitude information
struct HELMatrix {
  void InitAlphaToZero() {
    const unsigned int N = 20;
    alpha                = MMatrix<std::complex<double>>(N, N, 0.0);
    alpha_set            = MMatrix<bool>(N, N, true);
  }

  // Decay amplitude matrix (calculated by SU(2) decomposition routines)
  gra::MMatrix<std::complex<double>> T;

  // Helicity amplitude decay ls-couplings (constant)
  gra::MMatrix<std::complex<double>> alpha;
  gra::MMatrix<bool>                 alpha_set;

  // Parity conservation
  bool P_conservation = true;

  // --------------------------------------------------------------------
  // Tensor Pomeron couplings
  std::vector<double> g_decay_tensor;  // Decay couplings
  // --------------------------------------------------------------------

  // Branching Ratio to a particular decay
  double BR      = 1.0;
  double g_decay = 1.0;  // Equivalent decay coupling
};

}  // namespace gra

#endif