// Tensor Pomeron process initialization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORINIT_H
#define MTENSORINIT_H

#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Process/MProcessSetup.h"

namespace gra {

// Compute the beam restriction for the configured Tensor photo continuum
std::string TensorPhotoBeamError(nuclear::CollisionType collision, const nlohmann::json &general, const MPDG &pdg);

// Resolve the tune-specific Tensor Pomeron process model
void SetupTensorProcessModel(MProcessSetup &setup);

} // namespace gra

#endif
