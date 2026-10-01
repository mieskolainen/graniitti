// XP pole amplitudes and covariant numerator projections
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEXP_H
#define MREGGEXP_H

#include <complex>
#include <utility>
#include <vector>

#include "Graniitti/Regge/MReggeMPXP.h"

namespace gra::xpom {

// Share the common finite-spin pole operations
using rspin::LocalBasis;
using rspin::PairKernel;
using rspin::PhotoCoupling;

// Evaluate one XP resonance pole tensor
HelAmp Fusion(const LORENTZSCALAR& lts, const spin::PoleLS& vertex);

// Build the XP resonance production matrices
std::vector<HelAmp> Resonance(const LORENTZSCALAR& lts, const PARAM_RES& res, double s0, const std::array<ForwardLegState, 2>* forward_state = nullptr);

// Build the XP crossed continuum production matrices
std::vector<HelPair> Continuum(const LORENTZSCALAR& lts, double s0);

// Identify the covariant numerator carried by one internal XP field
enum class PoleType { Scalar, Dirac, Proca, Photon };

// Store one pole numerator in the physical helicity basis
struct PoleNumerator {
  PoleType            type = PoleType::Scalar;
  std::vector<double> helicities;
  HelAmp              matrix;
};

// Classify one supported scalar, fermion, vector or photon pole field
PoleType ClassifyPole(const MParticle& particle);

// Project a covariant numerator into the positive-energy pole basis for diagnostics
PoleNumerator Numerator(const MParticle& particle, const M4Vec& momentum);

// Compute the unit pole coefficient supported by reduced continuum vertices
PoleNumerator ReducedNumerator(const MParticle& particle, const M4Vec& momentum);

// Compute the antisymmetric photon field-strength tensor k^mu eps^nu-k^nu eps^mu
HelAmp FieldStrength(const M4Vec& momentum, int helicity);

// Express the reduced XP pole coefficient in the shared continuum helicity metric
spin::InternalHelicityMetric PoleMetric(const MParticle& particle, const M4Vec& momentum);


}  // namespace gra::xpom

#endif
