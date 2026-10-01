// Spin density and bipartite entanglement diagnostics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSPINDENSITY_H
#define MSPINDENSITY_H

#include "Graniitti/Spin/MHELMatrix.h"

namespace gra::spin {

// Store purity and von Neumann entropy in bits for one normalized spin density
struct SpinMetrics {
  double purity  = 0.0;
  double entropy = 0.0;
};

// Store daughter entropies and negativity for an ordered pair of spin spaces
struct PairMetrics {
  SpinMetrics pair;
  double      entropy1       = 0.0;
  double      entropy2       = 0.0;
  double      negativity     = 0.0;
  double      log_negativity = 0.0;
};

// Compute purity and entropy after normalizing a positive spin density to unit trace
SpinMetrics DensityMetrics(const MMatrix<std::complex<double>>& rho);

// Compute normalized bipartite metrics for rho_(i a),(j b) in the n1 x n2 product basis
PairMetrics EntanglementMetrics(const MMatrix<std::complex<double>>& rho, std::size_t n1, std::size_t n2);

// Compute the conditional daughter density F rho_parent F^dagger / Tr at fixed decay angles
MMatrix<std::complex<double>> DecayDensity(const HELMatrix& hel, const MMatrix<std::complex<double>>& parent,
                                           double theta, double phi);

// Print reduced-vertex pair metrics for an unpolarized mother at fixed daughter direction
void PrintDecayEntanglement(const HELMatrix& hel);

}  // namespace gra::spin

#endif
