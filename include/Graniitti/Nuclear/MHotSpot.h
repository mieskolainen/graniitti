// Subnucleon gluon-density fluctuations for incoherent photoproduction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARHOTSPOT_H
#define MNUCLEARHOTSPOT_H

#include <complex>
#include <cstddef>
#include <vector>

#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra::nuclear {

// Store transverse hotspot geometry and gluon-strength controls
struct HotSpotParam {
  unsigned int count          = 0;    // Hotspots per nucleon
  double       b_center       = 0.0;  // Hotspot-center width [GeV^-2]
  double       b_profile      = 0.0;  // Single-hotspot width [GeV^-2]
  double       strength_sigma = 0.0;  // Lognormal width

  // Validate hotspot inputs before constructing event configurations
  void Validate() const;
};

// Own one immutable subnucleon Good-Walker configuration bank
class MHotSpot {
 public:
  // Sample one hotspot bank aligned with a fixed nucleon configuration bank
  MHotSpot(const MConfigBank &bank, HotSpotParam param, MRandom &random);

  // Compute the validated hotspot controls
  const HotSpotParam &Param() const { return param_; }

  // Compute the number of nuclear configurations
  std::size_t ConfigCount() const { return config_count_; }

  // Compute the number of nucleons per configuration
  std::size_t NucleonCount() const { return nucleon_count_; }

  // Compute the analytic ensemble-mean transverse hotspot form factor
  double Mean(double qx, double qy) const;

  // Compute one sampled nucleon gluon-current factor
  std::complex<double> Factor(std::size_t sample, std::size_t nucleon, double qx, double qy) const;

 private:
  // Store one transverse hotspot center and normalized gluon strength
  struct Spot {
    double x        = 0.0;
    double y        = 0.0;
    double strength = 1.0;
  };

  HotSpotParam      param_;
  std::size_t       config_count_  = 0;
  std::size_t       nucleon_count_ = 0;
  std::vector<Spot> spot_;

  // Validate controls and sample all immutable hotspot configurations
  void Prepare(MRandom &random);
};

}  // namespace gra::nuclear

#endif
