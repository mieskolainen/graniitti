// Finite-range optical, GGCF and configuration Glauber survival
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARGLAUBER_H
#define MNUCLEARGLAUBER_H

#include <array>
#include <cstddef>
#include <memory>
#include <span>
#include <string>
#include <vector>

#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MFluctuation.h"
#include "Graniitti/Nuclear/MNucleus.h"

namespace gra::nuclear {

// Select a hadronic Glauber beam profile
enum class HadronType { Proton, Nucleus };

// Store the energy-dependent inelastic nucleon profile from the pp eikonal
struct NNProfile {
  double              sigma = 0.0;  // Integrated inelastic cross section [mb]
  double              omega = 0.0;  // One-projectile Good-Walker relative variance
  std::vector<double> b_node;       // Impact-parameter nodes [fm]
  std::vector<double> inelastic;    // Inelastic interaction probability
  std::string         fingerprint;  // Complete pp eikonal table identity
};

// Store finite-range Glauber inputs and numerical controls
struct GlauberParam {
  NNProfile        profile;              // Computed elementary pp inelastic profile
  double           b_max         = 0.0;  // Impact-parameter table limit [fm]
  double           q_max         = 0.0;  // Overlap transform limit [GeV]
  unsigned int     b_nodes       = 0;    // Impact-parameter table nodes
  unsigned int     q_nodes       = 0;    // Fourier-Bessel momentum nodes
  unsigned int     profile_nodes = 0;    // pp eikonal impact-profile intervals
  FluctuationParam fluctuation;          // Good-Walker cross-section eigenstates
};

// Own one immutable finite-range nucleus-nucleus Glauber model
class MGlauber {
 public:
  // Construct and tabulate one heterogeneous hadronic overlap
  MGlauber(std::array<HadronType, 2> type, std::array<std::shared_ptr<const MNucleus>, 2> nucleus, GlauberParam param);

  // Construct and tabulate one finite-range nuclear overlap
  MGlauber(const MNucleus& beam1, const MNucleus& beam2, GlauberParam param);

  // Compute the validated Glauber controls
  const GlauberParam& Param() const { return param_; }

  // Compute one checked beam nucleus
  const MNucleus* Nucleus(int leg) const;

  // Compute whether one checked beam uses the proton profile
  bool IsProton(int leg) const;

  // Compute the finite-range nuclear overlap in fm^-2
  double Overlap(double b) const;

  // Compute the optical no-additional-interaction amplitude
  double OpticalAmp(double b) const;

  // Compute the optical no-additional-interaction probability
  double OpticalProb(double b) const;

  // Compute the GGCF no-additional-interaction amplitude
  double GGCFAmp(double b) const;

  // Compute the GGCF no-additional-interaction probability
  double GGCFProb(double b) const;

  // Compute the finite-range nucleon no-interaction amplitude at one sigma
  double NNAmp(double b, double sigma) const;

  // Compute the configuration Glauber survival amplitude
  double ConfigAmp(double b, const MConfigBank* bank1, const MConfigBank* bank2) const;

  // Compute the configuration Glauber survival probability
  double ConfigProb(double b, const MConfigBank* bank1, const MConfigBank* bank2) const;

  // Compute the number of Cartesian configuration-pair samples
  std::size_t SampleCount(const MConfigBank* bank1, const MConfigBank* bank2) const;

  // Compute one paired configuration survival amplitude
  double SampleAmp(double b, const MConfigBank* bank1, const MConfigBank* bank2, std::size_t sample) const;

  // Compute one paired configuration amplitude at a transverse vector
  double SampleAmp(double bx, double by, const MConfigBank* bank1, const MConfigBank* bank2, std::size_t sample) const;

  // Compute one configuration-pair amplitude averaged over GGCF fluctuations
  double ConfigPairAmp(double bx, double by, const MConfigBank* bank1, const MConfigBank* bank2, std::size_t config1,
                       std::size_t config2) const;

  // Compute one GGCF-averaged configuration amplitude at many impact vectors
  std::vector<double> ConfigPairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                                    const MConfigBank* bank2, std::size_t config1, std::size_t config2) const;

  // Compute one configuration-pair amplitude at fixed GGCF eigenstates
  double ConfigPairAmp(double bx, double by, const MConfigBank* bank1, const MConfigBank* bank2, std::size_t config1,
                       std::size_t config2, std::size_t sigma1, std::size_t sigma2) const;

  // Compute one fixed-state configuration amplitude at many impact vectors
  std::vector<double> ConfigPairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                                    const MConfigBank* bank2, std::size_t config1, std::size_t config2,
                                    std::size_t sigma1, std::size_t sigma2) const;

  // Compute the number of cross-section fluctuation nodes
  std::size_t SigmaCount() const { return sigma_scale_.size(); }

  // Compute one leg cross section averaged over the other leg [mb]
  double Sigma(std::size_t node) const;

  // Compute the physical pair cross section at two independent fluctuation nodes [mb]
  double Sigma(std::size_t upper, std::size_t lower) const;

 private:
  std::array<HadronType, 2>                      type_;
  std::array<std::shared_ptr<const MNucleus>, 2> nucleus_;
  GlauberParam                                   param_;
  std::vector<double>                            b_node_;
  std::vector<double>                            overlap_;
  std::vector<double>                            sigma_scale_;
  std::vector<double>                            pair_scale_;
  std::vector<double>                            pair_radius_scale_;
  std::vector<double>                            pair_weight_;
  double                                         sigma_scale_max_      = 1.0;
  double                                         profile_step_         = 0.0;
  double                                         profile_inverse_step_ = 0.0;
  bool                                           profile_uniform_      = false;

  // Validate all physical and numerical controls
  void Validate() const;

  // Prepare the finite-range overlap interpolation table
  void Prepare();

  // Detect an equally spaced pp profile for direct interpolation lookup
  void PrepareProfile();

  // Prepare the moment-calibrated Gamma fluctuation and pair rules
  void PrepareSigma();

  // Compute one fixed-geometry amplitude averaged over fluctuation pairs
  double PairAmp(double bx, double by, const MConfig* config1, const MConfig* config2, std::size_t state) const;

  // Compute fixed or GGCF-averaged amplitudes using one spatial pair grid
  std::vector<double> PairAmp(const std::vector<std::array<double, 2>>& impact, const MConfigBank* bank1,
                              const MConfigBank* bank2, std::size_t config1, std::size_t config2,
                              std::size_t state) const;

  // Compute the packed symmetric GGCF pair-state index
  std::size_t PairState(std::size_t sigma1, std::size_t sigma2) const;

  // Compute the scaled inelastic nucleon probability at one separation
  double NNProb(double b, double sigma_scale) const;

  // Compute one prevalidated fluctuation-pair interaction probability
  double PairProb(double b, std::size_t state) const;

  // Accumulate one NN pair into the surviving fluctuation probabilities
  void Attenuate(double radius, std::size_t state, std::span<double> probability, std::size_t& active) const;

  // Average no-interaction amplitudes over the requested fluctuation states
  double Average(std::span<const double> probability, std::size_t state) const;

  // Select a validated nuclear configuration or a point proton on each leg
  std::array<const MConfig*, 2> Configs(const MConfigBank* bank1, const MConfigBank* bank2, std::size_t config1,
                                        std::size_t config2) const;

  // Compute the normalized inelastic nucleon form factor
  double NNForm(double q) const;

  // Validate one pair of optional configuration banks
  void ValidateBanks(const MConfigBank* bank1, const MConfigBank* bank2) const;

  // Compute one beam mass number for overlap normalization
  unsigned int MassNumber(int leg) const;

  // Compute one normalized beam matter form factor
  double MatterForm(int leg, double q) const;
};

}  // namespace gra::nuclear

#endif
