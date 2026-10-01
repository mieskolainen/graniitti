// MP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEMPINIT_H
#define MREGGEMPINIT_H

#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Spin/MPoleLS.h"

namespace gra::mpom {

// Prepare one resonance fusion vertex
spin::PoleLS PrepareResonance(MProcessSetup& setup, const PARAM_RES& res, const RES_PRODUCTION_CHANNEL& channel, const std::vector<MParticle>& legs);
// Prepare the spin filter in the selected production frame
void PreparePolarization(PARAM_RES& res, const LORENTZSCALAR& lts);
// Prepare four ordered continuum vertices
std::vector<spin::PoleResidue> PrepareContinuum(MProcessSetup& setup, const std::vector<MDecayBranch>& tree);
// Prepare a local continuum pair
ReggeContinuumPole PreparePair(MProcessSetup& setup, const MParticle& upper, const MParticle& lower, const MParticle& first, const MParticle& second);
// Print the selected MP polarization
void PrintPolarization(const PARAM_RES& res, bool metrics_enabled);

}  // namespace gra::mpom

#endif
