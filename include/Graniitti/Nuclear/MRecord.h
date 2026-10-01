// Nuclear UPC event record metadata
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARRECORD_H
#define MNUCLEARRECORD_H

#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "HepMC3/GenEvent.h"

namespace gra::nuclear {

// Attach complete UPC model metadata to one HepMC event
void AttachRecord(const MUPC &upc, const FinalState &final, HepMC3::GenEvent &event);

}  // namespace gra::nuclear

#endif
