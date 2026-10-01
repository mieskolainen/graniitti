// Free particle lattice paths and slit boundary conditions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_PATHMARK_H
#define PROGRAM_PATHMARK_H

#include <cmath>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MMath.h"

namespace gra::program {

struct PathParam {
  int           N       = 4;
  double        dt      = 1.0;
  double        m       = 1.5;
  int           k_slit  = 2;
  double        slit    = 0.02;
  unsigned long samples = 10000000;

  // Validate the source, two slit slices and detector before sampling
  void Validate(double width) const {
    if (N < 4 || k_slit < 2 || k_slit >= N - 1 || samples == 0 || !std::isfinite(dt) || dt <= 0.0 ||
        !std::isfinite(m) || m <= 0.0 || !std::isfinite(slit) || slit <= 0.0 || !std::isfinite(width) ||
        slit >= width) {
      throw std::invalid_argument("pathmark: invalid lattice, slit, mass, time step or sample count");
    }
  }
};

// Compute the free particle action with fixed time spacing
inline double S_lat(const std::vector<double>& x, const PathParam& param) {
  double action = 0.0;
  for (std::size_t j = 0; j + 1 < x.size(); ++j) { action += param.m / (2.0 * param.dt) * math::pow2(x[j + 1] - x[j]); }
  return action;
}

}  // namespace gra::program
#endif
