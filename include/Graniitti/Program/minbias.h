// Minimum bias component event counts
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef PROGRAM_MINBIAS_H
#define PROGRAM_MINBIAS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Graniitti/Tech/MAux.h"

namespace gra::program {

// Allocate SD, DD and ND events by the largest fractional remainders
inline std::array<int, 3> MinbiasCounts(int events, const std::array<double, 3>& xs) {
  using gra::aux::indices;
  if (events < 0) { throw std::invalid_argument("minbias: negative event count"); }
  for (const double value : xs) {
    if (!std::isfinite(value) || value < 0.0) {
      throw std::invalid_argument("minbias: component cross sections must be finite and nonnegative");
    }
  }
  const long double total = std::accumulate(xs.begin(), xs.end(), 0.0L);
  if (!(total > 0.0L)) { throw std::invalid_argument("minbias: inelastic cross section must be positive"); }
  std::array<int, 3>         counts{};
  std::array<long double, 3> remainder{};
  int                        left = events;
  for (const auto& i : indices(xs)) {
    const long double expected = static_cast<long double>(events) * xs[i] / total;
    counts[i]                  = static_cast<int>(std::floor(expected));
    remainder[i]               = expected - counts[i];
    left -= counts[i];
  }
  if (left < 0 || left > 3) { throw std::runtime_error("minbias: invalid event allocation"); }
  std::array<std::size_t, 3> order{0, 1, 2};
  std::stable_sort(order.begin(), order.end(),
                   [&](std::size_t a, std::size_t b) { return remainder[a] > remainder[b]; });
  for (int i = 0; i < left; ++i) { ++counts[order[i]]; }
  return counts;
}

}  // namespace gra::program
#endif
