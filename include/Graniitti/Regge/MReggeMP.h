// Minimal Pomeron central production amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEMP_H
#define MREGGEMP_H

#include <complex>
#include <utility>
#include <vector>

#include "Graniitti/Regge/MReggeMPXP.h"

namespace gra::mpom {

// Share the common finite-spin pole operations
using rspin::LocalBasis;
using rspin::PairKernel;
using rspin::PhotoCoupling;

// Evaluate one MP resonance pole tensor
HelAmp Fusion(const LORENTZSCALAR& lts, const spin::PoleLS& vertex);

// Build the MP resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, double s0, const std::array<ForwardLegState, 2>* forward_state = nullptr);

// Build the MP crossed continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, double s0);

}  // namespace gra::mpom

#endif
