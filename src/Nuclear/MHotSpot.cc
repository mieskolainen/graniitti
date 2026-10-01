// Subnucleon gluon-density fluctuations for incoherent photoproduction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MHotSpot.h"

#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>

#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra::nuclear {

// Validate the shared hotspot geometry and gluon-strength inputs
void HotSpotParam::Validate() const {
  if (count == 0 || !std::isfinite(b_center) || !std::isfinite(b_profile) || !std::isfinite(strength_sigma) ||
      b_center < 0.0 || b_profile < 0.0 || strength_sigma < 0.0) {
    throw std::invalid_argument("MHotSpot: invalid hotspot controls");
  }
}

// Sample one hotspot bank aligned with a fixed nucleon configuration bank
MHotSpot::MHotSpot(const MConfigBank &bank, HotSpotParam param, MRandom &random)
    : param_(param), config_count_(bank.Size()), nucleon_count_(bank.Nucleus().A()) {
  Prepare(random);
}

// Validate controls and sample all immutable hotspot configurations
// b_i ~ N(0, hbarc^2 B_center), w_i = exp(-sigma_w^2 / 2 + sigma_w z_i)
// [REFERENCE: Mantysaari and Schenke, Phys. Rev. Lett. 117 (2016) 052301]
void MHotSpot::Prepare(MRandom &random) {
  param_.Validate();
  const std::size_t count = static_cast<std::size_t>(param_.count);
  if (config_count_ > std::numeric_limits<std::size_t>::max() / nucleon_count_ ||
      config_count_ * nucleon_count_ > std::numeric_limits<std::size_t>::max() / count) {
    throw std::overflow_error("MHotSpot: configuration size overflow");
  }

  constexpr double hbarc          = PDG::GeV2fm;
  const double     center_sigma   = hbarc * std::sqrt(param_.b_center);
  const double     strength_shift = -0.5 * param_.strength_sigma * param_.strength_sigma;
  spot_.reserve(config_count_ * nucleon_count_ * count);
  for (std::size_t sample = 0; sample < config_count_; ++sample) {
    for (std::size_t nucleon = 0; nucleon < nucleon_count_; ++nucleon) {
      for (std::size_t hotspot = 0; hotspot < count; ++hotspot) {
        Spot spot;
        spot.x = random.G(0.0, center_sigma);
        spot.y = random.G(0.0, center_sigma);
        const double strength_fluctuation =
            std::fpclassify(param_.strength_sigma) == FP_ZERO ? 0.0 : param_.strength_sigma * random.G(0.0, 1.0);
        spot.strength = std::exp(strength_shift + strength_fluctuation);
        if (!std::isfinite(spot.strength)) { throw std::runtime_error("MHotSpot: non-finite gluon strength"); }
        spot_.push_back(spot);
      }
    }
  }
}

// Compute the analytic ensemble-mean transverse hotspot form factor
// <H(q_T)> = exp[-(B_center + B_profile) q_T^2 / 2]
double MHotSpot::Mean(const double qx, const double qy) const {
  if (!std::isfinite(qx) || !std::isfinite(qy)) {
    throw std::invalid_argument("MHotSpot::Mean: momentum must be finite");
  }
  const double qt2 = qx * qx + qy * qy;
  return std::exp(-0.5 * (param_.b_center + param_.b_profile) * qt2);
}

// Compute one sampled nucleon gluon-current factor
// H(q_T) = exp(-B_profile q_T^2 / 2) sum_i w_i exp(i q_T.b_i / hbarc) / N_h
// [REFERENCE: Mantysaari and Schenke, Phys. Rev. Lett. 117 (2016) 052301]
std::complex<double> MHotSpot::Factor(const std::size_t sample, const std::size_t nucleon, const double qx,
                                      const double qy) const {
  if (!std::isfinite(qx) || !std::isfinite(qy)) {
    throw std::invalid_argument("MHotSpot::Factor: momentum must be finite");
  }
  if (sample >= config_count_ || nucleon >= nucleon_count_) {
    throw std::out_of_range("MHotSpot::Factor: sample is out of range");
  }
  constexpr double     hbarc   = PDG::GeV2fm;
  const double         qt2     = qx * qx + qy * qy;
  std::complex<double> current = 0.0;
  const std::size_t    first   = (sample * nucleon_count_ + nucleon) * param_.count;
  for (std::size_t hotspot = 0; hotspot < param_.count; ++hotspot) {
    const Spot  &spot  = spot_[first + hotspot];
    const double phase = (qx * spot.x + qy * spot.y) / hbarc;
    current += std::polar(spot.strength, phase);
  }
  const double profile = std::exp(-0.5 * param_.b_profile * qt2);
  return profile * current / static_cast<double>(param_.count);
}

}  // namespace gra::nuclear
