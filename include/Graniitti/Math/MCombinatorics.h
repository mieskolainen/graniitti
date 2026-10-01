// Combinatorial and discrete mathematical algorithms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MCOMBINATORICS_H
#define MCOMBINATORICS_H

// C++
#include <algorithm>
#include <cstddef>
#include <iterator>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MTensor.h"

namespace gra {
namespace math {

// Index all combinations of numbers in vectors [recursive function]
// Output = S_0 x S_1 x ... x S_{N-1}
//
// -----------------------------------------------------------------------
// Test code:
//
// const std::vector<std::vector<size_t>> allind = {{0,1,2}, {0,1}, {0,1,2,3}};
//
// std::vector<size_t> current;
// std::vector<std::vector<std::size_t>> output;
// math::IndexComb(allind, 0, current, output);
//
// for (std::size_t i = 0; i < output.size(); ++i) {
//   for (std::size_t j = 0; j < output[i].size(); ++j) {
//     std::cout << output[i][j] << " ";
//   }
//   std::cout << std::endl;
// }
// -----------------------------------------------------------------------
template <typename T>
void IndexComb(const std::vector<std::vector<T>> &allind, std::size_t index, std::vector<T> current,
               std::vector<std::vector<T>> &output) {
  if (index >= allind.size()) {
    output.push_back(current);
    return;
  }
  for (std::size_t i = 0; i < allind[index].size(); ++i) {
    std::vector<T> b = current;
    b.push_back(allind[index][i]);
    IndexComb(allind, index + 1, b, output);
  }
}

// Compute the sign of one permutation relative to a reference ordering
// sgn(p) = (-1)^(number of inversions in p)
template <typename T>
int PermutationSign(const std::vector<T> &reference, const std::vector<T> &permuted) {
  if (reference.size() != permuted.size()) { throw std::invalid_argument("math::PermutationSign: size mismatch"); }

  std::vector<std::size_t> rank;
  rank.reserve(permuted.size());
  std::vector<bool> used(reference.size(), false);
  for (const auto &value : permuted) {
    const auto found = std::find(reference.cbegin(), reference.cend(), value);
    if (found == reference.cend()) { throw std::invalid_argument("math::PermutationSign: value mismatch"); }
    const std::size_t index = static_cast<std::size_t>(std::distance(reference.cbegin(), found));
    if (used[index]) { throw std::invalid_argument("math::PermutationSign: repeated value"); }
    used[index] = true;
    rank.push_back(index);
  }

  int sign = 1;
  for (std::size_t i = 0; i < rank.size(); ++i) {
    for (std::size_t j = i + 1; j < rank.size(); ++j) {
      if (rank[i] > rank[j]) { sign = -sign; }
    }
  }
  return sign;
}

// Binomial Coefficient C(n, k)
// C(n,k) = n!/[k!(n-k)!]
constexpr int Cbinom(int n, int k) {
  if (n < 0 || k < 0 || k > n) { return 0; }
  k                = std::min(k, n - k);
  long long result = 1;
  for (int i = 1; i <= k; ++i) {
    result = result * (n - k + i) / i;
    if (result > std::numeric_limits<int>::max()) {
      throw std::overflow_error("math::Cbinom: result does not fit in int");
    }
  }
  return static_cast<int>(result);
}

// N-dim epsilon tensor e_{\mu_1,\mu_2,\mu_3,...,\mu_N}:
// epsilon_{i_1...i_N} = sgn(i_1,...,i_N) for distinct indices and zero otherwise
//  + 1 if even permutation of arguments
//  - 1 if odd permutation of arguments
//    0 otherwise
inline MTensor<int> EpsTensor(std::size_t N) {
  if (N == 0) { throw std::invalid_argument("gra::math::EpsTensor: dimension must be positive"); }

  // Maximum range value for each for-loop
  const std::size_t MAX = N;

  // These hold for-loop index for each nested for loop
  std::vector<std::size_t> ind(N, 0);

  // ------------------------------------------------------------------
  // Permutation tensor definition
  auto permutation = [&]() {
    int value = 1;
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = i + 1; j < N; ++j) {
        // Even permutation +1, Odd permutation -1, Otherwise 0
        if (ind[i] > ind[j]) {
          value = -value;
        } else if (ind[i] == ind[j]) {
          return 0;
        }
      }
    }
    return value;
  };
  // ------------------------------------------------------------------

  const std::vector<std::size_t> dimensions(N, MAX);
  MTensor<int>                   T(dimensions, int(0));

  // Nested for-loop
  std::size_t index = 0;
  while (true) {
    // --------------------------------------------------------------
    // Evaluate the function value
    T(ind) = permutation();
    // --------------------------------------------------------------

    ind[0]++;

    // Carry
    while (ind[index] == MAX) {
      if (index == N - 1) { return T; }

      ind[index++] = 0;
      ind[index]++;
    }
    index = 0;
  }
}

// Construct binary reflect Gray code (BRGC)
// g(n) = n XOR (n >> 1)
// which is the Hamiltonian Path on a N-dim unit hypercube (2^N Boolean vector
// space)
// The operator ^ is Exclusive OR (XOR) and the operator >> is bit shift right
constexpr unsigned int Binary2Gray(unsigned int number) { return number ^ (number >> 1); }

// Convert BRGC (binary reflect Gray code) to a binary number
// n = XOR_{k >= 0} (g >> k)
// Gray code is obtained as XOR (^) with all more significant bits
constexpr unsigned int Gray2Binary(unsigned int number) {
  for (unsigned int mask = number >> 1; mask != 0; mask = mask >> 1) { number = number ^ mask; }
  return number;
}

// Boolean vector to index representation
// n = sum_i b_i 2^(d - 1 - i)
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
inline int Vec2Ind(const std::vector<bool> &vec) {
  if (vec.size() > static_cast<std::size_t>(std::numeric_limits<int>::digits)) {
    throw std::overflow_error("math::Vec2Ind: bit vector does not fit in int");
  }
  unsigned int result = 0;
  for (const bool bit : vec) { result = (result << 1U) | static_cast<unsigned int>(bit); }
  return static_cast<int>(result);
}

// Index to Boolean vector representation
// b_i = [n >> (d - 1 - i)] AND 1
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
// Input: x = index representation
//        d = Boolean vector space dimension
inline std::vector<bool> Ind2Vec(unsigned int ind, unsigned int d) {
  constexpr unsigned int DIGITS = std::numeric_limits<unsigned int>::digits;
  if (d > DIGITS || (d < DIGITS && ind >= (1U << d))) {
    throw std::overflow_error("math::Ind2Vec: index does not fit in requested dimension");
  }
  std::vector<bool> final(d, 0);
  for (unsigned int i = 0; i < d; ++i) { final[d - 1U - i] = ((ind >> i) & 1U) != 0U; }
  return final;
}

// OEIS.org A030109 sequence of left-right bit reversed sequence
// seq(i) is the d-bit reversal of i
// Input: dim = Boolean vector space dimension
inline std::vector<unsigned int> LRsequence(unsigned int dim) {
  if (dim > static_cast<unsigned int>(std::numeric_limits<int>::digits)) {
    throw std::overflow_error("math::LRsequence: dimension does not fit index representation");
  }
  // The sequence
  std::vector<unsigned int> seq(std::size_t{1} << dim, 0);

  for (std::size_t i = 0; i < seq.size(); ++i) {
    // Write i in binary with dimension d
    std::vector<bool> vec = Ind2Vec(i, dim);

    // Reverse the bits
    reverse(vec.begin(), vec.end());

    // Turn back to index representation
    seq[i] = Vec2Ind(vec);
  }

  return seq;
}

// Construct binary matrix in normal binary order
// B_{ij} = [i >> (d - 1 - j)] AND 1
inline std::vector<std::vector<int>> BinaryMatrix(unsigned int d) {
  if (d >= std::numeric_limits<unsigned int>::digits) {
    throw std::overflow_error("math::BinaryMatrix: dimension does not fit row count");
  }
  std::vector<std::vector<int>> B;
  const unsigned int            N = 1U << d;

  // First initialize
  std::vector<int> rvec(d, 0);
  for (std::size_t i = 0; i < N; ++i) { B.push_back(rvec); }

  // Now fill
  for (std::size_t i = 0; i < N; ++i) {
    // Get binary expansion (vector) for this i = 0...2^d-1
    std::vector<bool> binvec = Ind2Vec(i, d);

    for (unsigned int j = 0; j < d; ++j) { B[i][j] = binvec[j]; }
  }
  return B;
}

// Go through all the permutations, N! number of them
// Algorithm from the Knuth's "Art of Computer Programming"
//
// Input example:
//
// std::vector<int> set = {1,2,3};
inline std::vector<std::vector<int>> Permutations(const std::vector<int> &input) {
  std::vector<int> set = input;
  std::sort(set.begin(), set.end());
  std::vector<std::vector<int>> outset;
  do { outset.push_back(set); } while (std::next_permutation(set.begin(), set.end()));
  return outset;
}

// type == 0 : Returns charge conserving amplitude permutations, such as
// pi+pi-pi+pi- final state
// For N = 2,4,6,8 this gives number of 2,16,288,9216 permutations of final
// state amplitudes (legs)
// which is OEIS sequence
// A055546 (absolute coefficients of Cayley-Menger determinant of order N)
//

// type == 1 : Returns fully symmetrized amplitude, such as the case of yyyy (4
// gammas) final state
// This gives 2! 4! 6! (4,24,720) ... (pure factorial) number of permutations
//
inline std::vector<std::vector<int>> GetAmpPerm(int N, int type) {
  if (N < 0 || (type != 0 && type != 1)) {
    throw std::invalid_argument("math::GetAmpPerm: invalid multiplicity or permutation type");
  }
  if (type == 0 && (N == 0 || (N % 2) != 0)) {
    throw std::invalid_argument("math::GetAmpPerm: charged multiplicity must be positive and even");
  }
  // We start indexing at 3 for central system particles relevant for the
  // amplitude
  // indexing: [0 not in use, 1,2 are final state protons]
  int offset = 3;

  std::vector<int> set;
  int              k = offset;
  for (int i = 0; i < N; ++i) {
    set.push_back(k);
    ++k;
  }
  // Get complete permutations
  std::vector<std::vector<int>> outset = Permutations(set);

  if (type == 1) {
    return outset;         // Compute complete permutations
  } else if (type == 0) {  // Compute local charge conserving

    // Create alternating charge series
    // 1,-1,1,-1,...
    std::vector<int> charge;
    for (int i = 0; i < N; ++i) {
      charge.push_back(i % 2 == 0 ? 1 : -1);
    }

    std::vector<std::vector<int>> valid_amplitudes;

    for (std::size_t i = 0; i < outset.size(); ++i) {
      bool ok = true;
      for (std::size_t j = 0; j < outset[0].size() - 1; j += 2) {
        int a = outset[i][j];
        int b = outset[i][j + 1];

        if ((charge.at(a - offset) + charge.at(b - offset)) != 0) ok = false;
      }
      if (ok) { valid_amplitudes.push_back(outset[i]); }
    }
    return valid_amplitudes;
  }
  throw std::logic_error("math::GetAmpPerm: unreachable permutation type");
}


}  // namespace math
}  // namespace gra

#endif
