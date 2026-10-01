// Resolved electromagnetic excitation amplitudes and normalized history sampling
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEAREMD_H
#define MNUCLEAREMD_H

#include "Graniitti/Nuclear/MUPC.h"

namespace gra::nuclear {

// Use real scalar EMD amplitudes and trace orthogonal absorption and decay histories
// A_h = Fourier[S(b) sqrt(P(h|b)/q(h)) A(b)], sigma_class = sum_{h in class} q(h) |A_h|^2
class MEMD {
 public:
  struct State {
    std::array<double, 2> excitation{};
    std::shared_ptr<const ExcitationChannel> channel;
  };

  // Prepare the compound Poisson response and long-range Hankel quadrature
  explicit MEMD(std::shared_ptr<const MUPC> upc);

  // Sample h with probability q(h), supplying sqrt(P(h|b)/q(h)) before Fourier transformation
  State Sample(MRandom& random) const;

 private:
  std::shared_ptr<const MUPC> upc_;
  std::array<bool, 2> absorbed_{};
  std::vector<double> radius_, momentum_, survival_, proposal_, log_condition_;
  std::array<std::vector<double>, 2> mean_;
  std::array<MMatrix<double>, 2> log_spectrum_;
  MMatrix<double> kernel_;
};

}  // namespace gra::nuclear

#endif
