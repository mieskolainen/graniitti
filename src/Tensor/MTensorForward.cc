// Tensor Pomeron forward source decomposition
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorForward.h"

#include <algorithm>
#include <stdexcept>

namespace gra {

// Compute the checked storage index for one forward mechanism
std::size_t TensorForwardSourceBank::MechanismIndex(
    const TensorForwardMechanism mechanism) {
  switch (mechanism) {
  case TensorForwardMechanism::PomeronPomeron:
    return 0;
  case TensorForwardMechanism::GammaPomeron:
    return 1;
  case TensorForwardMechanism::PomeronGamma:
    return 2;
  case TensorForwardMechanism::GammaGamma:
    return 3;
  }
  throw std::invalid_argument(
      "TensorForwardSourceBank::MechanismIndex: unknown mechanism");
}

// Insert one Tensor forward label while preserving its first component index
std::size_t
TensorForwardSourceBank::Add(const std::array<TensorForwardSource, 2> &label) {
  const auto found = std::find(labels.begin(), labels.end(), label);
  if (found != labels.end()) {
    return static_cast<std::size_t>(std::distance(labels.begin(), found));
  }
  labels.push_back(label);
  return labels.size() - 1;
}

// Build the common PP, gamma-P and P-gamma source basis
TensorForwardSourceBank
TensorForwardSourceBank::Hadronic(const ForwardLegState &upper,
                                  const ForwardLegState &lower) {
  TensorForwardSourceBank bank;
  const auto pomeron_source = [](const ForwardLegState &state) {
    return state.IsExcited() ? TensorForwardSource::InclusivePomeron
                             : TensorForwardSource::Elastic;
  };
  const auto photon_source = [](const ForwardLegState &state,
                                const std::size_t eigen_source) {
    if (!state.IsExcited()) {
      return TensorForwardSource::Elastic;
    }
    return eigen_source == 0 ? TensorForwardSource::PhotonParallel
                             : TensorForwardSource::PhotonPerpendicular;
  };

  bank.mechanism_components[MechanismIndex(
      TensorForwardMechanism::PomeronPomeron)] = {
      bank.Add({pomeron_source(upper), pomeron_source(lower)})};

  auto &gamma_pomeron = bank.mechanism_components[MechanismIndex(
      TensorForwardMechanism::GammaPomeron)];
  const std::size_t upper_sources = upper.IsExcited() ? 2 : 1;
  for (std::size_t source = 0; source < upper_sources; ++source) {
    gamma_pomeron.push_back(
        bank.Add({photon_source(upper, source), pomeron_source(lower)}));
  }

  auto &pomeron_gamma = bank.mechanism_components[MechanismIndex(
      TensorForwardMechanism::PomeronGamma)];
  const std::size_t lower_sources = lower.IsExcited() ? 2 : 1;
  for (std::size_t source = 0; source < lower_sources; ++source) {
    pomeron_gamma.push_back(
        bank.Add({pomeron_source(upper), photon_source(lower, source)}));
  }
  return bank;
}

// Build the Cartesian gamma-gamma photon-density source basis
TensorForwardSourceBank
TensorForwardSourceBank::PhotonFusion(const ForwardLegState &upper,
                                      const ForwardLegState &lower) {
  TensorForwardSourceBank bank;
  const auto photon_source = [](const ForwardLegState &state,
                                const std::size_t eigen_source) {
    if (!state.IsExcited()) {
      return TensorForwardSource::Elastic;
    }
    return eigen_source == 0 ? TensorForwardSource::PhotonParallel
                             : TensorForwardSource::PhotonPerpendicular;
  };
  auto &gamma_gamma = bank.mechanism_components[MechanismIndex(
      TensorForwardMechanism::GammaGamma)];
  const std::size_t upper_sources = upper.IsExcited() ? 2 : 1;
  const std::size_t lower_sources = lower.IsExcited() ? 2 : 1;
  for (std::size_t upper_source = 0; upper_source < upper_sources;
       ++upper_source) {
    for (std::size_t lower_source = 0; lower_source < lower_sources;
         ++lower_source) {
      gamma_gamma.push_back(bank.Add({photon_source(upper, upper_source),
                                      photon_source(lower, lower_source)}));
    }
  }
  return bank;
}

// Compute component indices for one ordered interaction
const std::vector<std::size_t> &TensorForwardSourceBank::Components(
    const TensorForwardMechanism mechanism) const {
  return mechanism_components[MechanismIndex(mechanism)];
}

// Compute the upper and lower source labels of one component
const std::array<TensorForwardSource, 2> &
TensorForwardSourceBank::Label(const std::size_t component) const {
  if (component >= labels.size()) {
    throw std::out_of_range(
        "TensorForwardSourceBank::Label: component outside source bank");
  }
  return labels[component];
}

} // namespace gra
