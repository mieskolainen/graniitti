// Dispersive rho and omega propagation into charged pions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORVECTOR_H
#define MTENSORVECTOR_H

#include <array>
#include <complex>

#include "Graniitti/Math/MMatrix.h"

namespace gra::tensor {

// Compute s R(s,m^2) + i I(s,m^2) for a two-pseudoscalar P-wave cut
std::complex<double> VectorLoop(double s, double mass);

// Store the common pion vertices and pole inputs of the coupled vector propagator
struct RhoOmega {
  std::array<double, 2> mass{};
  std::array<double, 2> g{};
  double width = 0.0;
  double pion = 0.0;
  double kaon = 0.0;
  double b = 0.0;

  // Compute the rho or omega row after species validation at initialization
  std::size_t Index(int pdg) const { return pdg == 113 ? 0 : 1; }

  // Compute the inverse transverse propagator with on-shell real subtractions
  MMatrix<std::complex<double>> Inverse(double s) const;

  // Compute the coupled transverse propagator for timelike vector decay
  MMatrix<std::complex<double>> Propagator(double s) const;
};

}  // namespace gra::tensor

#endif
