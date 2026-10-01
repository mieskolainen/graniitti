// XP production vertex preparation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEXPINIT_H
#define MREGGEXPINIT_H

#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Spin/MPoleLS.h"

namespace gra::xpom {

// Prepare one resonance fusion vertex
spin::PoleLS PrepareResonance(MProcessSetup &setup, const PARAM_RES &res, const RES_PRODUCTION_CHANNEL &channel, const std::vector<MParticle> &legs);
// Prepare four ordered continuum vertices
std::vector<spin::PoleResidue> PrepareContinuum(MProcessSetup &setup, const std::vector<MDecayBranch> &tree);
// Prepare a local continuum pair
ReggeContinuumPole PreparePair(MProcessSetup &setup, const MParticle &upper, const MParticle &lower, const MParticle &first, const MParticle &second);
// Print the canonical XP pole operator
void PrintPole(const std::string &label, const MParticle &mother, const std::vector<MParticle> &legs, const spin::PoleLS &vertex);

}  // namespace gra::xpom

#endif
