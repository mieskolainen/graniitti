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

// MadGraph
#include "Graniitti/Amplitude/MAmpMG5_yy_ll.h"
#include "Graniitti/Amplitude/MAmpMG5_yy_ww.h"

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MKinematics.h"

namespace gra {
// "Functionoid class"
// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MGamma {
 public:
  MGamma() {}
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
  // MADGRAPH amplitudes added here
  MAmpMG5_yy_ll AmpMG5_yy_ll;
  MAmpMG5_yy_ww AmpMG5_yy_ww;
};

}  // namespace gra

#endif
