// Gamma amplitudes
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGAMMA_H
#define MGAMMA_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAmplitudes.h"
#include "Graniitti/MKinematics.h"

namespace gra {

// Monopole wave functions and parameters
namespace PARAM_MONOPOLE {

extern bool printed;      // Printing called
extern bool initialized;  // For lazy initialization

extern int         En;        // Bound state energy level
extern double      M0;        // Monopole mass
extern double      Gamma0;    // Monopolium width
extern std::string coupling;  // Coupling scenarios
extern int         gn;        // Monopole charge n = 1,2,3,...

double EnergyMP(double n);
double GammaMP(double n, double alpha_g);
double PsiMP(double n);

void PrintParameters(double sqrts);
void ReadParameters(const std::string &modelfile);

}  // namespace PARAM_MONOPOLE


class MGamma : public MAmplitudes {
 public:
  MGamma(gra::LORENTZSCALAR &lts, const std::string &modelfile);
  ~MGamma() {}

  // yy->resonance X
  double yyX(gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const;

  // yy->lepton pair, or monopole antimonopole amplitude
  double yyffbar(gra::LORENTZSCALAR &lts);

  // yy->SM Higgs
  double yyHiggs(gra::LORENTZSCALAR &lts) const;

  // yy->monopolium
  double yyMP(gra::LORENTZSCALAR &lts) const;

 protected:
};

}  // namespace gra

#endif
