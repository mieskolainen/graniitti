// Photon and other fluxes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFLUX_H
#define MFLUX_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MKinematics.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MGlobals.h"

namespace gra {
namespace flux {


  double ApplyktEPAfluxes(double amp2, gra::LORENTZSCALAR& lts);
  double ApplyDZfluxes(double amp2, gra::LORENTZSCALAR& lts);
  double ApplyLUXfluxes(double amp2, gra::LORENTZSCALAR& lts);


} // namespace flux
} // namespace gra


#endif