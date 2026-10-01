// Particle parameters evaluated by generated MG5 models
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_MODEL_H
#define AMP_MG5_MODEL_H

#include <cmath>
#include <limits>
#include <map>
#include <stdexcept>
#include <string>
#include <utility>

#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"

namespace gra::mg5 {

// Store the real pole parameters used by the generated model
struct ParticleParam {
  double mass = 0.0;
  double width = 0.0;
  bool signed_mass = false;
};

using ParticleMap = std::map<int, ParticleParam>;

// Reject card entries incompatible with the evaluated UFO model
inline void ValidateModel(const SLHAReader &card, const ParticleMap &particles) {
  for (const auto &[pdg, particle] : particles) {
    for (const auto &[block, value] : {std::pair{"mass", particle.mass},
                                      std::pair{"decay", particle.width}}) {
      if (!std::isfinite(value) || value < 0.0) {
        throw std::invalid_argument("MG5: invalid model " + std::string(block) + " for PDG " + std::to_string(pdg));
      }
      const auto entry = card.find_block_entry(block, pdg);
      // Generated cards print dependent parameters with seven significant digits
      const double tolerance = 5.1e-7 * std::abs(value) + 64.0 * std::numeric_limits<double>::epsilon();
      // [REFERENCE: hep-ph/0311123, section 2.4.1, signed Majorana masses]
      const bool signed_value = particle.signed_mass &&
          (std::string(block) == "mass" || card.find_block_entry("mass", pdg).value_or(0.0) < 0.0);
      const double card_value = entry.has_value() ? (signed_value ? std::abs(*entry) : *entry) : value;
      if (entry.has_value() && (!std::isfinite(card_value) || std::abs(card_value - value) > tolerance)) {
        throw std::invalid_argument("MG5: card " + std::string(block) + " for PDG " + std::to_string(pdg) +
                                    " conflicts with the generated model. Change the independent model parameters or regenerate");
      }
    }
  }
}

}  // namespace gra::mg5

#endif
