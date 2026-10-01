// Nuclear UPC model construction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARSETUP_H
#define MNUCLEARSETUP_H

// C++
#include <array>
#include <memory>
#include <optional>
#include <ostream>

// Own
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Particle/MParticle.h"

namespace gra::nuclear {

// Build one validated UPC runtime from physical beam states
std::shared_ptr<const MUPC> BuildUPC(const std::array<MParticle, 2> &beam, const std::array<M4Vec, 2> &momentum,
                                    int excitation, const UPCParam &param, UPCMode mode = {});

// Print the active nuclear UPC setup fields
void PrintUPCSetup(const MUPC &upc, std::ostream &stream);

// Derive the full elementary inelastic profile from one pp eikonal table
NNProfile BuildNNProfile(const MEikonal &eikonal, std::optional<double> sigma_nn);

}  // namespace gra::nuclear

#endif
