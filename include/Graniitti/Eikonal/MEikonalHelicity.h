// Coupled-channel elastic proton helicity amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEIKONAL_HELICITY_H
#define MEIKONAL_HELICITY_H

#include <array>
#include <complex>
#include <cstddef>
#include <stdexcept>

// Own
#include "Graniitti/Spin/MHelicity.h"

namespace gra {

class MEikonal;

// Proton-proton amplitudes in the negative-first collider helicity convention
//   phi1 = <++|T|++>, phi2 = <++|T|-->, phi3 = <+-|T|+->,
//   phi4 = <+-|T|-+>, phi5 = <++|T|+->
// <++|T|-+> = -phi5_first in the same external-state phase section
struct ElasticHelicityAmplitudes {
  std::complex<double> phi1 = 0.0;
  std::complex<double> phi2 = 0.0;
  std::complex<double> phi3 = 0.0;
  std::complex<double> phi4 = 0.0;
  std::complex<double> phi5 = 0.0;
  std::complex<double> phi5_first = 0.0;
};

// Named independent amplitudes of the elastic proton helicity operator
enum class ElasticHelicityComponent { Phi1, Phi2, Phi3, Phi4, Phi5, Phi5First };

// One row-major proton pair helicity matrix transition in (--,-+,+-,++) order
struct ProtonHelicityTransition {
  ElasticHelicityComponent component = ElasticHelicityComponent::Phi1;
  double sign = 1.0;
  int azimuth_harmonic = 0;
};

// Compute the collider harmonic for one row-major proton pair transition
constexpr int ProtonPairHelicityHarmonic(const std::size_t row,
                                         const std::size_t col) {
  const auto initial = spin::BinaryPairHelicityLabelsX2(col);
  const auto final = spin::BinaryPairHelicityLabelsX2(row);
  return spin::ColliderSpinHalfHelicityHarmonic(initial[0], initial[1],
                                                final[0], final[1]);
}

// Compute the row-major table with signed collider helicity reciprocity
constexpr std::array<ProtonHelicityTransition, 16>
CanonicalProtonHelicityTransitions() noexcept {
  using Component = ElasticHelicityComponent;
  return {{{Component::Phi1, 1.0, ProtonPairHelicityHarmonic(0, 0)},
           {Component::Phi5, -1.0, ProtonPairHelicityHarmonic(0, 1)},
           {Component::Phi5First, 1.0, ProtonPairHelicityHarmonic(0, 2)},
           {Component::Phi2, 1.0, ProtonPairHelicityHarmonic(0, 3)},
           {Component::Phi5, 1.0, ProtonPairHelicityHarmonic(1, 0)},
           {Component::Phi3, 1.0, ProtonPairHelicityHarmonic(1, 1)},
           {Component::Phi4, 1.0, ProtonPairHelicityHarmonic(1, 2)},
           {Component::Phi5First, 1.0, ProtonPairHelicityHarmonic(1, 3)},
           {Component::Phi5First, -1.0, ProtonPairHelicityHarmonic(2, 0)},
           {Component::Phi4, 1.0, ProtonPairHelicityHarmonic(2, 1)},
           {Component::Phi3, 1.0, ProtonPairHelicityHarmonic(2, 2)},
           {Component::Phi5, -1.0, ProtonPairHelicityHarmonic(2, 3)},
           {Component::Phi2, 1.0, ProtonPairHelicityHarmonic(3, 0)},
           {Component::Phi5First, -1.0, ProtonPairHelicityHarmonic(3, 1)},
           {Component::Phi5, 1.0, ProtonPairHelicityHarmonic(3, 2)},
           {Component::Phi1, 1.0, ProtonPairHelicityHarmonic(3, 3)}}};
}

// Compute the defining reference-plane matrix index for one named amplitude
constexpr std::size_t
ElasticHelicityReferenceIndex(const ElasticHelicityComponent component) {
  switch (component) {
  case ElasticHelicityComponent::Phi1:
    return spin::PairHelicityMatrixIndex(3, 3);
  case ElasticHelicityComponent::Phi2:
    return spin::PairHelicityMatrixIndex(3, 0);
  case ElasticHelicityComponent::Phi3:
    return spin::PairHelicityMatrixIndex(2, 2);
  case ElasticHelicityComponent::Phi4:
    return spin::PairHelicityMatrixIndex(2, 1);
  case ElasticHelicityComponent::Phi5:
    return spin::PairHelicityMatrixIndex(3, 2);
  case ElasticHelicityComponent::Phi5First:
    return spin::PairHelicityMatrixIndex(3, 1);
  }
  throw std::invalid_argument(
      "ElasticHelicityReferenceIndex: unknown component");
}

// Compute the defining reference-plane transition for one named amplitude
constexpr ProtonHelicityTransition
ElasticHelicityReferenceTransition(const ElasticHelicityComponent component) {
  return CanonicalProtonHelicityTransitions()[ElasticHelicityReferenceIndex(
      component)];
}

using ProtonHelicityMatrix = std::array<std::complex<double>, 16>;

// Construct the complete pair-helicity matrix at one transfer azimuth
ProtonHelicityMatrix
BuildProtonHelicityMatrix(const ElasticHelicityAmplitudes &amplitude,
                          double azimuth = 0.0);

// Rotate one complete reference-plane helicity matrix to a transfer azimuth
ProtonHelicityMatrix
RotateProtonHelicityMatrix(const ProtonHelicityMatrix &reference,
                           double azimuth);

// Screened crossing decomposition.  The amplitude for the constructed beam
// state is A_even + sigma_beam A_odd, with sigma_beam=+1 for pp and -1 for
// p pbar. This decomposition is performed after matrix exponentiation
struct ElasticCrossingHelicityAmplitudes {
  ElasticHelicityAmplitudes even;
  ElasticHelicityAmplitudes odd;
};

// Project the coupled Good-Walker matrix onto one exclusive physical final
// state.  f1=f2=0 is elastic pp -> pp (or ppbar -> ppbar)
ElasticHelicityAmplitudes GetElasticHelicityAmplitudes(const MEikonal &eikonal,
                                                       double kt2,
                                                       std::size_t f1 = 0,
                                                       std::size_t f2 = 0);

// Compute the separately screened C-even and C-odd helicity amplitudes
ElasticCrossingHelicityAmplitudes
GetElasticCrossingHelicityAmplitudes(const MEikonal &eikonal, double kt2,
                                     std::size_t f1 = 0, std::size_t f2 = 0);

// Spin-averaged amplitude squared with the scalar normalization preserved
double
UnpolarizedHelicityAmpSquared(const ElasticHelicityAmplitudes &amplitude);

} // namespace gra

#endif
