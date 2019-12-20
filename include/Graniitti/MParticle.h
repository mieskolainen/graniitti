// Particle and branch objects [HEADER ONLY file]
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPARTICLE_H
#define MPARTICLE_H

// C++
#include <complex>
#include <random>
#include <valarray>
#include <vector>

#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MHELMatrix.h"
#include "Graniitti/MCW.h"

namespace gra {


// Particle class
struct MParticle {
  // These values are read from the PDG-file
  std::string name;
  int         pdg      = 0;  // PDG code
  int         chargeX3 = 0;  // Q x 3
  int         spinX2   = 0;  // J x 2
  int         color    = 0;  // QCD color code

  double mass  = 0.0;
  double width = 0.0;
  double tau   = 0.0;  // hbar / width

  // J^PC
  int          P    = 1;      // Default even P-parity
  int          C    = 0;      // Default zero C-parity (e.g. fermions do not have)
  unsigned int L    = 0;      // For Mesons/Baryons
  bool         glue = false;  // Glueball state

  void setPCL(int _P, int _C, unsigned int _L) {
    P = _P;
    C = _C;
    L = _L;
  }

  // Width cut (in Breit-Wigner sampling)
  double wcut = 0.0;

  void print() const {
    std::cout << " NAME:   " << name << std::endl;
    std::cout << " ID:     " << pdg << std::endl;
    std::cout << " M:      " << mass << std::endl;
    std::cout << " W:      " << width << std::endl;
    std::cout << " J^PC:   " << aux::Spin2XtoString(spinX2) << "^" << aux::ParityToString(P)
              << aux::ParityToString(C) << std::endl
              << std::endl;
  }
};


// Recursive decay tree branch
struct MDecayBranch {
  MDecayBranch() {
    f = MMatrix<std::complex<double>>(1, 1, 1.0);  // Unit element
  }

  // Offshell mass picked event by event
  double m_offshell = 0.0;

  MParticle                 p;               // PDG particle
  M4Vec                     p4;              // 4-momentum
  std::vector<MDecayBranch> legs;            // Daughters
  M4Vec                     decay_position;  // Decay 4-position
  gra::HELMatrix            hel;             // Decay helicity information

  // Used with helicity amplitudes
  MMatrix<std::complex<double>> f;

  // MC weight container
  gra::kinematics::MCW W;

  // Active in the factorized phase space product (by default, no)
  // This is controlled by the spesific amplitudes
  bool PS_active = false;

  // Decay tree current level
  int depth = 0;
};

}  // namespace gra

#endif