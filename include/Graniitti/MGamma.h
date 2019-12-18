// Gamma amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGAMMA_H
#define MGAMMA_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MAmplitudes.h"

namespace gra {

class MGamma : public MAmplitudes {
 public:
  MGamma();
  ~MGamma();

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
