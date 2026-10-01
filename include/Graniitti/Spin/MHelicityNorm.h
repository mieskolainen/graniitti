// Forward helicity normalization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELICITYNORM_H
#define MHELICITYNORM_H

// C++
#include <complex>

// Own
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Spin/MHELMatrix.h"

namespace gra::spin {

// Identify the forward source measure of one exchanged leg
enum class ForwardLegType { Hadron, RealPhoton };

// Compute the physical reduced-helicity norm without basis padding
double ReducedHelicityNorm2(const HELMatrix &helicity);

// Compute the leading forward density of one coherent central helicity tensor
double ForwardHelicityDensity(const MMatrix<std::complex<double>> &tensor,
                              const MMatrix<double> &incoming_helicities,
                              ForwardLegType upper, ForwardLegType lower,
                              double tolerance = 1e-12);

} // namespace gra::spin

#endif
