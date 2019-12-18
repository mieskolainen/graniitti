// All external amplitudes (matrix elements) from MadGraph etc. collected here
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MAMPLITUDES_H
#define MAMPLITUDES_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MKinematics.h"

// MadGraph
#include "Graniitti/Amplitude/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/AMP_MG5_yy_ll.h"
#include "Graniitti/Amplitude/AMP_MG5_yy_ww_evev.h"
#include "Graniitti/Amplitude/AMP_MG5_gg_gg.h"
#include "Graniitti/Amplitude/AMP_MG5_gg_ggg.h"


namespace gra {

// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MAmplitudes {

public:
  
  MAmplitudes() {}
  ~MAmplitudes() {}

  AMP_MG5_yy_ww AmpMG5_yy_ww;
  AMP_MG5_yy_ll AmpMG5_yy_ll;
  AMP_MG5_yy_ww_evev AmpMG5_yy_ww_evev;
  AMP_MG5_gg_gg  AmpMG5_gg_gg;
  AMP_MG5_gg_ggg AmpMG5_gg_ggg;

protected:

};

}  // namespace gra

#endif
