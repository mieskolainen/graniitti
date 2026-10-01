// HERA proton dissociation profiles for vector meson photoproduction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPHOTODISS_H
#define MPHOTODISS_H

#include <map>

#include "json.hpp"

namespace gra::flux {

// Store the dissociative density relative to the elastic forward normalization at W0
struct PhotoDissParam {
  double W0        = 0.0;
  double ratio     = 0.0;
  double delta     = 0.0;
  double b         = 0.0;
  double n         = 0.0;
  double epsilon   = 0.0;
  double mass_max  = 0.0;
  double mass_norm = 0.0;
};

// Read and normalize the measured proton dissociation profiles
std::map<int, PhotoDissParam> ReadPhotoDiss(const nlohmann::json& rows);

// Compute the relative amplitude density in dM_Y^2 with the fitted W and t dependence
double PhotoDissFactor(const PhotoDissParam& param, double w2, double t, double mass2);

}  // namespace gra::flux

#endif
