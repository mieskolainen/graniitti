// HepMC event-record construction shared by physical processes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEVENTRECORD_H
#define MEVENTRECORD_H

#include "Graniitti/Process/MProcessState.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

namespace gra::record {

// Write the common central-production event topology
bool WriteCentral(MProcessState &state, HepMC3::GenEvent &evt);

// Attach standard HepMC color-flow attributes to a generated particle
void AttachColorFlow(const MParticle &particle, const HepMC3::GenParticlePtr &generated);

// Sample lifetimes and write a decay branch with an optional prompt root
void WriteBranch(MDecayBranch &branch, const HepMC3::GenParticlePtr &mother, HepMC3::GenEvent &evt,
                 MRandom &random, bool displaced = true);

}  // namespace gra::record

#endif
