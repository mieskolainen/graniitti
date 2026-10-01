// Good-Walker fluctuation input for nuclear survival
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARGGCF_H
#define MNUCLEARGGCF_H

#include <string>

#include "Graniitti/Eikonal/MEikonal.h"

namespace gra::nuclear {

// Store the GGCF survival eikonal selection
struct GGCFParam {
  std::string eikonal;
};

// Store the forward total cross section and one-projectile diffractive variance
struct GWMoments {
  double sigma = 0.0;  // Total NN cross section [mb]
  double omega = 0.0;  // One-projectile forward dissociation / elastic
};

// Integrate complex forward amplitudes with the other proton kept elastic
GWMoments ForwardGW(const MEikonalMatrix &runtime, unsigned int leg = 0);

// Select and validate the configured multichannel eikonal snapshot
MModelTunePtr GGCFTune(const MModelTunePtr &tune, const GGCFParam &param);

}  // namespace gra::nuclear
#endif
