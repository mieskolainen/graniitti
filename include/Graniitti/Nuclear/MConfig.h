// Correlated nuclear proton and neutron configurations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARCONFIG_H
#define MNUCLEARCONFIG_H

#include <array>
#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MTypes.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra::nuclear {

// Compute current moments with an unbiased, centered configuration variance
CurrentStat CurrentMoments(const std::vector<std::complex<double>>& current);

// Own one fixed proton and neutron configuration
class MConfig {
 public:
  // Construct one validated immutable nucleon configuration
  MConfig(unsigned int a, unsigned int z, std::vector<Nucleon> nucleon);

  // Compute all nucleon centers
  const std::vector<Nucleon> &Nucleons() const { return nucleon_; }

  // Compute nucleon indices sorted by transverse x coordinate
  const std::vector<std::size_t> &TransverseOrder() const { return transverse_order_; }

  // Compute transverse coordinate bounds as xmin, xmax, ymin and ymax
  const std::array<double, 4> &TransverseBounds() const { return transverse_bounds_; }

  // Compute exact downstream transverse separations for one produced nucleon
  std::span<const double> DownstreamR2(std::size_t produced, bool positive_z) const;

  // Compute the coherent phase sum over proton centers for q in GeV
  std::complex<double> ChargeCurrent(double qx, double qy, double qz = 0.0) const;

  // Compute the coherent phase sum over all nucleon centers for q in GeV
  std::complex<double> MatterCurrent(double qx, double qy, double qz = 0.0) const;

 private:
  unsigned int                            a_ = 0;
  unsigned int                            z_ = 0;
  std::vector<Nucleon>                    nucleon_;
  std::vector<std::size_t>                proton_;
  std::vector<std::size_t>                transverse_order_;
  std::array<double, 4>                   transverse_bounds_{};
  std::array<std::vector<std::size_t>, 2> downstream_offset_;
  std::array<std::vector<double>, 2>      downstream_r2_;

  // Compute one selected configuration phase sum
  std::complex<double> Current(double qx, double qy, double qz, bool charge_only) const;

  // Tabulate exact downstream pair geometry for both photon directions
  void PreparePhotoGeometry();
};

// Own immutable density samplers shared by event-local nuclear configurations
class MConfigSampler {
 public:
  // Construct all deterministic radial sampling tables
  MConfigSampler(const MNucleus &nucleus, ConfigParam param);

  // Compute the nuclear density model used for sampling
  const MNucleus &Nucleus() const { return nucleus_; }

  // Compute the validated configuration controls
  const ConfigParam &Param() const { return param_; }

  // Draw one correlated proton and neutron configuration
  MConfig Draw(MRandom &random) const;

  // Map one uniform variate to the normalized neutron radial density
  double NeutronRadius(double u) const;

 private:
  MNucleus            nucleus_;
  ConfigParam         param_;
  std::vector<double> neutron_r_;
  std::vector<double> neutron_cdf_;

  // Validate configuration controls and prepare the neutron radial sampler
  void Prepare();
};

// Own one event-local configuration bank for Good-Walker averages
class MConfigBank {
 public:
  // Sample a correlated configuration bank from one prepared density model
  MConfigBank(const MConfigSampler &sampler, MRandom &random);

  // Sample a bounded event-local subset from one prepared density model
  MConfigBank(const MConfigSampler &sampler, MRandom &random, std::size_t count);

  // Compute the validated sampling controls
  const ConfigParam &Param() const { return param_; }

  // Compute the nuclear density model used for sampling
  const MNucleus &Nucleus() const { return nucleus_; }

  // Compute the number of fixed configurations
  std::size_t Size() const { return config_.size(); }

  // Compute one fixed configuration with checked indexing
  const MConfig &At(std::size_t index) const;

  // Compute charge-current ensemble moments for q in GeV
  CurrentStat ChargeStat(double qx, double qy, double qz = 0.0) const;

  // Compute matter-current ensemble moments for q in GeV
  CurrentStat MatterStat(double qx, double qy, double qz = 0.0) const;

 private:
  MNucleus             nucleus_;
  ConfigParam          param_;
  std::vector<MConfig> config_;

  // Compute current moments from all fixed configurations
  CurrentStat Stat(double qx, double qy, double qz, bool charge_only) const;
};

}  // namespace gra::nuclear

#endif
