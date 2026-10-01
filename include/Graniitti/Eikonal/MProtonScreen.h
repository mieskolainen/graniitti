// Universal pp screening of helicity and Good Walker amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPROTONSCREEN_H
#define MPROTONSCREEN_H

#include <complex>
#include <cstddef>
#include <span>
#include <vector>

#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra::eikonal {

// Apply the pp eikonal to one physical amplitude and any color projections
class MProtonScreen {
 public:
  // Expand Born amplitudes into the physical proton and Good Walker basis
  MProtonScreen(const std::vector<std::vector<std::complex<double>>> &born, const ScreeningMetadata &metadata,
                const MEikonal::LoopConst &loop, GWFinalClass final_class);

  // Accumulate S_pp(kT) A_hard(kT) d2kT for every amplitude bank
  void Add(std::size_t radial, std::size_t azimuth,
           const std::vector<std::span<const std::complex<double>>> &amplitude);

  // Compute A_screened = A_Born + A_loop for every amplitude bank
  std::vector<std::vector<std::complex<double>>> Result() const;

 private:
  struct Bank {
    std::vector<std::complex<double>> born;
    std::vector<std::complex<double>> loop;
    std::size_t                       source_size = 0;
    std::size_t                       spectators  = 0;
  };

  ScreeningMetadata                                   metadata_;
  const MEikonal::LoopConst                          &loop_;
  const std::vector<MEikonal::LoopGoodWalkerChannel> *channel_ = nullptr;
  std::vector<Bank>                                   bank_;
  bool                                                good_walker_        = false;
  bool                                                scalar_good_walker_ = false;
  bool                                                helicity_           = false;
};

// Apply the pp eikonal directly in the Good Walker pair space
class MProtonGoodWalkerScreen {
 public:
  // Initialize dense pair-space blocks from the Born amplitude
  MProtonGoodWalkerScreen(const ProtonGoodWalkerAmplitude &born, const ScreeningMetadata &metadata,
                          const MEikonal::LoopConst &loop, const SoftModelPtr &soft_model,
                          std::size_t eikonal_channels);

  // Accumulate S_pp(kT) A_hard(kT) d2kT in proton spin and Good Walker pair space
  void Add(std::size_t radial, std::size_t azimuth, const ProtonGoodWalkerAmplitude &amplitude);

  // Project the screened pair-space amplitude onto physical Good Walker final states
  std::vector<std::complex<double>> Result(spin::MQMetrics *metrics = nullptr) const;

 private:
  const ProtonGoodWalkerAmplitude               &born_;
  ScreeningMetadata                              metadata_;
  const MEikonal::LoopConst                     &loop_;
  std::vector<std::vector<std::complex<double>>> block_;
  std::vector<std::size_t>                       spectators_;
};

}  // namespace gra::eikonal

#endif
