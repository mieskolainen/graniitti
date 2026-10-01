// Regge photoproduction beam contraction for proton and nuclear collisions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEPHOTO_H
#define MREGGEPHOTO_H

#include <map>

#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MRegge.h"

namespace gra {

// Store selected VMD widths owned by the Regge photo adapter
using ReggePhotoWidths = std::map<int, double>;

// Build the immutable VMD data required by selected nuclear photo channels
ReggePhotoWidths BuildPhotoPlan(const LORENTZSCALAR &lts);

// Contract model-specific elementary directions with proton or nuclear beams
double EvalReggePhoto(MRegge &regge, LORENTZSCALAR &lts, ReggeProductionModel model, const ReggePhotoWidths &width_ee);

}  // namespace gra

#endif
