// Cascade phase-space proposal densities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDECAYSAMPLING_H
#define MDECAYSAMPLING_H

#include "Graniitti/Process/MProcessState.h"

namespace gra::decay {

// Compute the invariant-mass proposal density for one internal branch
double MassDensity(const MDecayBranch &branch);

// Compute the mixture density relative to the generated central and decay measure
double MixtureDensity(LORENTZSCALAR &lts, const std::vector<MDecayBranch> &tree);

}  // namespace gra::decay

#endif
