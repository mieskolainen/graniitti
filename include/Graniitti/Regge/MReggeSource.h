// Regge forward sources and Good Walker pair amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGESOURCE_H
#define MREGGESOURCE_H

#include <array>
#include <complex>
#include <cstdint>
#include <deque>
#include <optional>
#include <utility>
#include <vector>

#include "Graniitti/Regge/MRegge.h"

namespace gra {

// Store one forward source before the two-leg Cartesian product
struct LegSource {
  ProtonGoodWalkerSector            sector = ProtonGoodWalkerSector::Elastic;
  std::vector<std::complex<double>> nonflip;
  std::vector<std::complex<double>> flip;
};

// Store one pair-space source with its orthogonal forward sectors
struct PairSource {
  ProtonGoodWalkerSector        upper_sector = ProtonGoodWalkerSector::Elastic;
  ProtonGoodWalkerSector        lower_sector = ProtonGoodWalkerSector::Elastic;
  MMatrix<std::complex<double>> amplitude;
};

// Reuse forward Regge sources only within one fixed amplitude call
class ReggeSourceCache {
 public:
  // Bind the cache to the complete forward state of one amplitude call
  ReggeSourceCache(const LORENTZSCALAR &lts, const MRegge &regge, const regge::Param &param, const ForwardLegState &upper, const ForwardLegState &lower);

  // Bind explicit subenergies and optional photon-source amplitudes
  ReggeSourceCache(const LORENTZSCALAR &lts, const MRegge &regge, const regge::Param &param, const ForwardLegState &upper, const ForwardLegState &lower, std::array<double, 2> subenergy,
                   std::array<std::optional<std::complex<double>>, 2> photon_source);

  // Compute one cached Good Walker profile
  const std::vector<LegSource> &Profile(ForwardBeamLeg leg, int exchange_pdg);

  // Compute one cached scalar Regge kernel
  std::complex<double> Kernel(ForwardBeamLeg leg, int exchange_pdg, double s_forward, bool photoproduction, int central_pdg);

  // Compute one cached profile with its scalar kernel applied
  const std::vector<LegSource> &Sources(ForwardBeamLeg leg, int exchange_pdg, double s_forward, bool photoproduction, int central_pdg);

  // Compute one caller supplied forward subenergy
  double Subenergy(ForwardBeamLeg leg) const;

 private:
  using ProfileKey = std::pair<ForwardBeamLeg, int>;

  // Store a discrete cache identity without floating point comparisons
  struct KernelKey {
    ForwardBeamLeg leg             = ForwardBeamLeg::Upper;
    int            exchange_pdg    = 0;
    std::uint64_t  subenergy       = 0;
    bool           photoproduction = false;
    int            central_pdg     = 0;

    bool operator==(const KernelKey &) const = default;
  };

  // Compute the immutable forward state selected at construction
  const ForwardLegState &State(ForwardBeamLeg leg) const;

  // Canonicalize only the inputs used by one scalar Regge kernel
  KernelKey Key(ForwardBeamLeg leg, int exchange_pdg, double s_forward, bool photoproduction, int central_pdg) const;

  const LORENTZSCALAR                                       &lts;
  const MRegge                                              &regge;
  const regge::Param                                        &param;
  ForwardLegState                                            upper;
  ForwardLegState                                            lower;
  std::array<double, 2>                                      subenergy;
  std::array<std::optional<std::complex<double>>, 2>         photon_source;
  std::size_t                                                max_profiles = 0;
  std::vector<std::pair<ProfileKey, std::vector<LegSource>>> profiles;
  std::vector<std::pair<KernelKey, std::complex<double>>>    kernels;
  std::deque<std::pair<KernelKey, std::vector<LegSource>>>   sources;
};

// Compute (-1)^N when N internal lines connect reduced M subamplitudes
// This follows after removing the overall i from the full iM graph
double SewingSign(std::size_t internal_lines);

// Apply one common factor to every proton Good Walker component
void ScaleGoodWalker(LORENTZSCALAR &lts, double scale, const std::string &context);

// Compute the Born norm with the configured initial spin average
double BornNorm(const LORENTZSCALAR &lts);

// Retain only the proton Good Walker source during shifted screening
bool GoodWalkerOnly(LORENTZSCALAR &lts, const std::string &context);

// Initialize or validate the event-local proton Good Walker amplitude
void PrepareGoodWalker(LORENTZSCALAR &lts, const SoftModelPtr &model);

// Compute the hard proton spin layout used by Good Walker sources
ScreeningMetadata ProtonSpinLayout(const LORENTZSCALAR &lts);

// Build the Good Walker profile carried by one forward leg
std::vector<LegSource> LegProfile(const LORENTZSCALAR &lts, const ForwardLegState &state, const MRegge &regge, const regge::Param &param, int exchange_pdg, std::optional<std::complex<double>> photon_source = std::nullopt);

// Compute the kernel and beam-sign factor carried by one forward leg
std::complex<double> LegKernel(const LORENTZSCALAR &lts, const ForwardLegState &state, const MRegge &regge, const regge::Param &param, int exchange_pdg, double s_forward, bool photoproduction, int central_pdg);

// Build all orthogonal sources carried by one forward leg
std::vector<LegSource> LegSources(const LORENTZSCALAR &lts, const ForwardLegState &state, const MRegge &regge, const regge::Param &param, int exchange_pdg, double s_forward, bool photoproduction, int central_pdg);

// Form scalar-spin upper and lower source combinations
std::vector<PairSource> PairSources(const std::vector<LegSource> &upper, const std::vector<LegSource> &lower, const GoodWalkerSpace &good_walker, std::complex<double> scale);

// Form upper and lower source combinations in the proton spin basis
std::vector<PairSource> PairSources(const std::vector<LegSource> &upper, const std::vector<LegSource> &lower, const GoodWalkerSpace &good_walker, const ScreeningMetadata &spin_layout, std::complex<double> scale);

// Add one central helicity matrix to compatible pair sources
void AddGoodWalker(LORENTZSCALAR &lts, const std::vector<PairSource> &pairs, const MMatrix<std::complex<double>> &central, std::size_t group, const SoftModelPtr &model);

// Add one production spin source only when density collection is active
void AddGoodWalker(LORENTZSCALAR &lts, const std::vector<PairSource> &pairs, const MMatrix<std::complex<double>> &central, std::size_t group, const SoftModelPtr &model, std::size_t spin_index);

// Project every proton Good Walker source onto its physical final-state sector
std::vector<std::complex<double>> ProjectGoodWalker(const ProtonGoodWalkerAmplitude &state, spin::MQMetrics *metrics = nullptr, double normalization = 1.0, std::size_t proton_rows = 0);

}  // namespace gra

#endif
