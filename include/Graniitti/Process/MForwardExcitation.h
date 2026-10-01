// Forward proton excitation and fragmentation state
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFORWARDEXCITATION_H
#define MFORWARDEXCITATION_H

#include <array>

#include "Graniitti/Process/MProcessState.h"

namespace gra::forward {

// Compute and validate the invariant mass interval of an excited proton
std::array<double, 2> MassBounds(const MProcessState &state, const GENCUT &cuts);

// Clear generated forward branches after event serialization
void ClearBranches(LORENTZSCALAR &lts) noexcept;

// Clear all event-local forward excitation and diffractive variables
void Reset(LORENTZSCALAR &lts) noexcept;

}  // namespace gra::forward

#endif
