// Custom user defined cuts
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MUSERCUTS_H
#define MUSERCUTS_H

// C++
#include <cstdint>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {
// Compute whether an integer selects an implemented user cut
bool IsKnownUserCut(std::int64_t id) noexcept;

// User cuts (return false for events not passing the cuts)
bool UserCut(std::int64_t id, const gra::LORENTZSCALAR &lts) noexcept;

// User cuts with an explicit post-radiation central momentum tree
bool UserCut(std::int64_t id, const gra::LORENTZSCALAR &lts,
             const std::vector<MDecayBranch> &central) noexcept;

}  // namespace gra

#endif
