// Nuclear UPC steering parser
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARSTEERING_H
#define MNUCLEARSTEERING_H

#include <string>

// Own
#include "Graniitti/Nuclear/MUPC.h"

// Libraries
#include "json.hpp"

namespace gra::nuclear {

// Parse selectors and complete nuclear model-card blocks
UPCParam ReadUPCSteering(const nlohmann::json &block, const nlohmann::json &physics, const nlohmann::json &numerics,
                         const std::string &tune);

}  // namespace gra::nuclear

#endif
