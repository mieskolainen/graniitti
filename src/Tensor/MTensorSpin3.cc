// Covariant spin-three decay and effective photoproduction current
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorSpin3.h"

namespace gra::tensor {

// Project the decay momentum onto the massive spin-three representation, all indices down
FTensor::Tensor3<double, 4, 4, 4> Spin3Decay(const M4Vec &positive, const M4Vec &negative) {
  const auto p = positive + negative;
  const auto difference = positive - negative;
  const auto r = difference - p * ((difference * p) / p.M2());
  const auto metric = [&p](int a, int b) {
    return (a == b ? (a == 0 ? 1.0 : -1.0) : 0.0) - (p % a) * (p % b) / p.M2();
  };
  FTensor::Tensor3<double, 4, 4, 4> out;
  for (int a = 0; a < 4; ++a) {
    for (int b = 0; b < 4; ++b) {
      for (int c = 0; c < 4; ++c) {
        out(a, b, c) = (r % a) * (r % b) * (r % c) - r.M2() / 5.0 *
            (metric(a, b) * (r % c) + metric(a, c) * (r % b) + metric(b, c) * (r % a));
      }
    }
  }
  return out;
}

// Contract the field strength with the transverse spin-three decay tensor
FTensor::Tensor3<double, 4, 4, 4> Spin3Current(
    const M4Vec &q, const M4Vec &p, const FTensor::Tensor3<double, 4, 4, 4> &decay) {
  FTensor::Tensor3<double, 4, 4, 4> out;
  for (int b = 0; b < 4; ++b) {
    for (int c = 0; c < 4; ++c) {
      double contraction = 0.0;
      for (int a = 0; a < 4; ++a) { contraction += q[a] * decay(a, b, c); }
      for (int mu = 0; mu < 4; ++mu) { out(mu, b, c) = (q * p) * decay(mu, b, c) - (p % mu) * contraction; }
    }
  }
  return out;
}

}  // namespace gra::tensor
