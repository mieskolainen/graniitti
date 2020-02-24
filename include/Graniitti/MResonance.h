// Resonance container class
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MRESONANCE_H
#define MRESONANCE_H

// C++
#include <complex>
#include <map>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MHELMatrix.h"
#include "Graniitti/MParticle.h"

namespace gra {

// Resonance parameters
class PARAM_RES {
 public:
  PARAM_RES() {}
  void PrintParam(double sqrts) const {
    std::cout << rang::style::bold << "Custom resonance parameters:" << rang::style::reset
              << std::endl
              << std::endl;

    printf("- PDG ID:      %d \n", p.pdg);
    printf("- Mass M0:     %0.5f GeV \n", p.mass);
    printf("- Width W:     %0.5f GeV \n", p.width);
    printf("- Spin Jx2:    %d \n", p.spinX2);
    printf("- Parity P:    %d \n", p.P);
    std::cout << std::endl;

    printf("<Production> \n");
    printf("- Effective vertex constant g_i:  %0.1E x exp(i x %0.1f) \n", std::abs(g), std::arg(g));
    printf("- Form factor parameter: %0.2f \n", g_FF);

    std::cout << std::endl;
    printf("<Decay> \n");
    printf("- BR:          %0.3E\n", hel.BR);
    printf("- Effective vertex constant g_f:  %0.3E", hel.g_decay);
    std::cout << std::endl << std::endl;

    if (p.spinX2 != 0) {
      std::cout << "- Polarization Lorentz frame: " << FRAME << std::endl;
      std::cout << std::endl;
      std::cout << "- Spin polarization density matrix [rho]:" << std::endl << std::endl;

      // Print elements
      for (std::size_t i = 0; i < rho.size_row(); ++i) {
        for (std::size_t j = 0; j < rho.size_col(); ++j) {
          std::string delim = (j < rho.size_col() - 1) ? ", " : "";
          printf("%6.3f+i%6.3f%s ", std::real(rho[i][j]), std::imag(rho[i][j]), delim.c_str());
        }
        std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
      // printf("Density matrix von Neumann entropy S = %0.3f
      // \n\n", MSpin::VonNeumannEntropy(rho));
    }
  }

  // Particle class
  MParticle p;

  // -------------------------------------------------------------------
  // Tensor Pomeron couplings
  std::vector<double> g_Tensor;  // Production couplings
  // -------------------------------------------------------------------

  // (Complex) production coupling constant
  std::complex<double> g = 0.0;

  // Form factor parameter
  double g_FF = 0.0;

  // Breit-Wigner type
  int BW = 0;

  // 2->1 (generation spin correlation), 1->2 (decay spin correlations)
  bool SPINGEN = true;
  bool SPINDEC = true;

  // Lorentz frame of the spin distribution
  std::string FRAME = "null";

  // Maximum sliding pomeron helicity
  int JMAX = 0;

  // Spin-density matrix (constant)
  MMatrix<std::complex<double>> rho;

  // --------------------------------------------------------------------
  // Helicity and decay amplitude information
  HELMatrix hel;
  // --------------------------------------------------------------------
};

}  // namespace gra

#endif
