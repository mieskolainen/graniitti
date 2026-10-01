// Tensor Pomeron process classification and decay validation
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORPROCESS_H
#define MTENSORPROCESS_H

#include <string>
#include <vector>

#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {

// Select the Tensor Pomeron amplitude family exposed through the process
enum class MTensorPomeronMode {
  Generic,
  Resonance,
  Continuum,
  ResonanceContinuum,
  Photo,
  QED
};

namespace tensor {

// Classify resonance quantum numbers supported by Tensor amplitudes
TensorResonanceType ClassifyTensorResonance(const MParticle &particle);

// Match the decay-tree structure accepted by one Tensor Pomeron mode
bool ProcessAccepts(const std::vector<MDecayBranch> &tree,
                    MTensorPomeronMode mode);

// Compute the stable process identity of one Tensor Pomeron process mode
std::string ProcessName(MTensorPomeronMode mode);

// Compute the final-state pattern exposed by one Tensor Pomeron mode
std::string ProcessPattern(MTensorPomeronMode mode);

// Resolve the decay structure for one accepted Tensor Pomeron topology
DecayStructure ProcessDecayStructure(MTensorPomeronMode mode,
                                     const LORENTZSCALAR &lts);

} // namespace tensor
} // namespace gra

#endif
