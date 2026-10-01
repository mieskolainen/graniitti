// Covariant spin-three decay and effective photoproduction current
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef GRANIITTI_MTENSORSPIN3_H
#define GRANIITTI_MTENSORSPIN3_H

#include "Graniitti/Kinematics/M4Vec.h"
#include "FTensor.hpp"

namespace gra::tensor {

// Compute the symmetric transverse traceless rank-three scalar-pair decay tensor
FTensor::Tensor3<double, 4, 4, 4> Spin3Decay(const M4Vec &positive, const M4Vec &negative);

// Compute the conserved photon current of (d^nu R^{abc}) F_{a nu} P_{bc}
FTensor::Tensor3<double, 4, 4, 4> Spin3Current(
    const M4Vec &q, const M4Vec &p, const FTensor::Tensor3<double, 4, 4, 4> &decay);

}  // namespace gra::tensor
#endif
