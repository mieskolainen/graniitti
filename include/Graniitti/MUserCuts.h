// Custom user defined cuts
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MUSERCUTS_H
#define MUSERCUTS_H

// C++
#include <vector>

// Own
#include "Graniitti/MKinematics.h"

namespace gra {
// User cuts (return false for events not passing the cuts)
bool UserCut(int id, const gra::LORENTZSCALAR &lts);

}  // gra namespace ends

#endif