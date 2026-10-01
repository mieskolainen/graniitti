// Tensor Pomeron forward source decomposition
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSORFORWARD_H
#define MTENSORFORWARD_H

#include <array>
#include <cstddef>
#include <vector>

#include "Graniitti/Kinematics/MEventKinematics.h"

namespace gra {

// Select one ordered Tensor forward interaction
enum class TensorForwardMechanism {
  PomeronPomeron,
  GammaPomeron,
  PomeronGamma,
  GammaGamma
};

// Select the source carried by one generated forward leg
enum class TensorForwardSource {
  Elastic,
  InclusivePomeron,
  PhotonParallel,
  PhotonPerpendicular
};

// Store orthogonal forward-source components of one Tensor amplitude
class TensorForwardSourceBank {
public:
  // Build the common PP, gamma-P and P-gamma source basis
  static TensorForwardSourceBank Hadronic(const ForwardLegState &upper,
                                          const ForwardLegState &lower);

  // Build the Cartesian gamma-gamma photon-density source basis
  static TensorForwardSourceBank PhotonFusion(const ForwardLegState &upper,
                                              const ForwardLegState &lower);

  // Compute the number of incoherent source components
  std::size_t Size() const noexcept { return labels.size(); }

  // Compute component indices for one ordered interaction
  const std::vector<std::size_t> &
  Components(TensorForwardMechanism mechanism) const;

  // Compute the upper and lower source labels of one component
  const std::array<TensorForwardSource, 2> &Label(std::size_t component) const;

private:
  // Compute the checked storage index for one forward mechanism
  static std::size_t MechanismIndex(TensorForwardMechanism mechanism);

  // Insert one source label and return its stable component index
  std::size_t Add(const std::array<TensorForwardSource, 2> &label);

  std::vector<std::array<TensorForwardSource, 2>> labels;
  std::array<std::vector<std::size_t>, 4> mechanism_components;
};

} // namespace gra

#endif
