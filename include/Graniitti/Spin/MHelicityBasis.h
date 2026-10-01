// Common binary helicity basis order and pair row mappings
// 
// Canonical convention used in GRANIITTI:
// 
// Single helicity: (-,+)
// Pair helicity: (--,-+,+-,++), second label varies fastest
// Physical matrices: rows are final states, columns initial states
// Hard amplitudes: initial pair major, final pair minor
// Higher-spin states: -J,...,+J
// Azimuthal phase:  e^{i(\lambda_{1i} - \lambda_{2i}
//                      - \lambda_{1f} + \lambda_{2f})\phi}
// 
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MHELICITY_BASIS_H
#define MHELICITY_BASIS_H

// C++
#include <array>
#include <cstddef>

namespace gra {
namespace spin {

// Compute the Condon-Shortley section of an integer spherical component
constexpr double SphericalSection(const int m) {
  return m > 0 && m % 2 != 0 ? -1.0 : 1.0;
}

// Compute both binary helicity indices in storage order
constexpr std::array<std::size_t, 2> BinaryHelicityIndices() noexcept {
  return {0, 1};
}

// Compute one negative-first doubled spin-half or transverse helicity label
constexpr int BinaryHelicityLabelX2(const std::size_t index) {
  return index == 0 ? -1 : 1;
}

// Compute the negative-first index of one doubled binary helicity label
constexpr std::size_t BinaryHelicityIndexX2(const int helicity_x2) {
  return helicity_x2 > 0 ? 1 : 0;
}

// Compute both binary helicity labels in negative-first storage order
constexpr std::array<int, 2> BinaryHelicityLabelsX2() {
  return {BinaryHelicityLabelX2(0), BinaryHelicityLabelX2(1)};
}

// Compute one pair index with the second helicity varying fastest
constexpr std::size_t BinaryPairHelicityIndex(const std::size_t first,
                                              const std::size_t second) {
  return 2 * first + second;
}

// Compute one pair index from physical helicity labels
constexpr std::size_t BinaryPairHelicityIndexX2(const int first_x2,
                                                const int second_x2) {
  return BinaryPairHelicityIndex(BinaryHelicityIndexX2(first_x2),
                                 BinaryHelicityIndexX2(second_x2));
}

// Compute negative-first physical helicity labels for one pair index
constexpr std::array<int, 2>
BinaryPairHelicityLabelsX2(const std::size_t pair) {
  return {BinaryHelicityLabelX2(pair / 2), BinaryHelicityLabelX2(pair % 2)};
}

// Compute a row-major operator index with final pair rows and initial columns
constexpr std::size_t PairHelicityMatrixIndex(const std::size_t final_pair,
                                              const std::size_t initial_pair) {
  return 4 * final_pair + initial_pair;
}

// Compute a hard amplitude row with initial pair first and final pair second
constexpr std::size_t
PairHelicityTransitionIndex(const std::size_t initial_pair,
                            const std::size_t final_pair) {
  return 4 * initial_pair + final_pair;
}

} // namespace spin
} // namespace gra

#endif // MHELICITY_BASIS_H
