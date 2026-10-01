// Gauge-restored charged-meson photoproduction in the Tensor Pomeron model
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORPHOTO_H
#define MTENSORPHOTO_H

#include <array>
#include <complex>

#include "Graniitti/Tensor/MTensorPomeron.h"

namespace gra {

using TensorPhotoCurrent = std::array<FTensor::Tensor1<std::complex<double>, 4>, 4>;

// Contract a conserved target current with physical photon sources and nuclear transitions
nuclear::PhotoTerm TensorPhotoTerm(const MTensorPomeron& tensor, const MTensorPomeronParam& param,
                                   const LORENTZSCALAR& lts, const ForwardLegState& target,
                                   const TensorPhotoCurrent& current, const flux::PhotoTargetProfile& profile,
                                   bool continuum);

// Evaluate the unequal-subenergy Drell-Soding amplitude
class MTensorPhoto {
 public:
  using Current = std::array<std::complex<double>, 4>;

  // Bind shared Tensor vertices and immutable photoproduction steering
  MTensorPhoto(const MTensorPomeron& tensor, MTensorPomeronParamPtr parameters);

  // Compute one gamma-star proton Drell-Soding current for fixed target spins
  Current DSCurrent(const LORENTZSCALAR& lts, bool photon_from_upper, std::size_t initial_helicity,
                    std::size_t final_helicity) const;

  // Build conserved meson currents in both directions with the common beam contraction
  std::vector<nuclear::PhotoTerm> Terms(const LORENTZSCALAR& lts) const;

  // Evaluate both physical photon directions in the collider process
  double Amp2(LORENTZSCALAR& lts) const;

 private:
  const MTensorPomeron&      tensor;
  MTensorPomeronParamPtr     parameter_handle;
  const MTensorPomeronParam& parameter;
};

}  // namespace gra

#endif
