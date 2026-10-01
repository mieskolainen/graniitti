// Good-Walker fluctuation input for nuclear survival
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MGGCF.h"

#include <cmath>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::nuclear {
using gra::aux::indices;

// Integrate Gamma_f0(b) before the incoherent sum over orthogonal final states
// omega_N = [d sigma(pp -> Xp)/dt] / [d sigma_el/dt] at t = 0
// This covers the configured low-mass GW states and assumes predominantly imaginary amplitudes
GWMoments ForwardGW(const MEikonalMatrix &runtime, const unsigned int leg) {
  if (leg > 1) { throw std::invalid_argument("ForwardGW: invalid projectile leg"); }
  const auto &b = runtime.ImpactParameterNodes();
  std::vector<ProtonHelicityMatrix> forward(runtime.ChannelCount());
  for (const auto &channel : indices(forward)) {
    const auto f1 = leg == 0 ? channel : 0, f2 = leg == 1 ? channel : 0;
    auto lower = runtime.ImpactHelicityMatrix(b.front(), f1, f2);
    for (std::size_t i = 1; i < b.size(); ++i) {
      const auto upper = runtime.ImpactHelicityMatrix(b[i], f1, f2);
      for (const auto &spin : indices(lower)) {
        if (CanonicalProtonHelicityTransitions()[spin].azimuth_harmonic != 0) { continue; }
        forward[channel][spin] += math::PI * (b[i] - b[i - 1]) * (b[i - 1] * lower[spin] + b[i] * upper[spin]);
      }
      lower = upper;
    }
  }
  const auto &elastic = forward.front();
  const double norm = gra::SquaredNorm(elastic);
  double diss = 0.0;
  for (std::size_t channel = 1; channel < forward.size(); ++channel) { diss += gra::SquaredNorm(forward[channel]); }
  const double sigma = 0.5 * PDG::GeV2mb * std::imag(elastic[0] + elastic[5] + elastic[10] + elastic[15]);
  if (!(norm > 0.0) || !(sigma > 0.0) || !std::isfinite(diss) || !std::isfinite(norm) || !std::isfinite(sigma)) {
    throw std::runtime_error("ForwardGW: invalid forward diffraction amplitudes");
  }
  return {sigma, diss / norm};
}

// Resolve the named GGCF model independently of the active production model
MModelTunePtr GGCFTune(const MModelTunePtr &tune, const GGCFParam &param) {
  const auto selected = tune->WithSoft(param.eikonal);
  if (selected->Soft()->GoodWalker().ChannelCount() < 2) {
    throw std::invalid_argument("PARAM_NUCLEAR.GGCF.eikonal requires a multichannel model");
  }
  return selected;
}

}  // namespace gra::nuclear
