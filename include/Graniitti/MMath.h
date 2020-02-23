// Mathematical functions [HEADER ONLY file]
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATH_H
#define MMATH_H

// C++
#include <algorithm>
#include <complex>
#include <random>
#include <stdexcept>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MTensor.h"

namespace gra {
namespace math {
// Imaginary unit
constexpr const std::complex<double> zi(0.0, 1.0);

// Geometry
constexpr const double PI   = 3.141592653589793238462643383279502884197169399375105820974944L;
constexpr const double PIPI = PI * PI;

// Index all combinations of numbers in vectors [recursive function]
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

// Moore-Penrose Pseudoinverse for Eigen
template <typename T>
T PseudoInverse(const T &a, double epsilon = std::numeric_limits<double>::epsilon()) {
  Eigen::JacobiSVD<T> svd(a, Eigen::ComputeThinU | Eigen::ComputeThinV);

  const double tolerance =
      epsilon * std::max(a.cols(), a.rows()) * svd.singularValues().array().abs()(0);

  return svd.matrixV() *
         (svd.singularValues().array().abs() > tolerance)
             .select(svd.singularValues().array().inverse(), 0)
             .matrix()
             .asDiagonal() *
         svd.matrixU().adjoint();
}

// Matrix-Vector Multiplication C = AB,
// where A is (n x m)
//       B is (m x 1),
//       which gives C (n x 1)
template <typename T>
inline std::vector<T> vv_vMultiply(const std::vector<std::vector<T>> &A, const std::vector<T> &B) {
  const std::size_t n = A.size();
  const std::size_t m = A[0].size();
  std::vector<T>    C(n, 0.0);  // Init with zero!

  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t k = 0; k < m; ++k) {
      C[i] += A[i][k] * B[k];  // notice plus
    }
  }
  return C;
}

template <typename T>
inline std::vector<T> valarray2vector(const std::valarray<T> &x) {
  std::vector<T> y(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) { y[k] = x[k]; }
  return y;
}
template <typename T>
inline std::valarray<T> vector2valarray(const std::vector<T> &x) {
  std::valarray<T> y(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) { y[k] = x[k]; }
  return y;
}

// Absolute values per vector element
template <typename T>
inline std::vector<T> vabs(const std::vector<T> &x) {
  std::vector<T> y(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) { y[k] = std::abs(x[k]); }
  return y;
}

// MATLAB style, H is std::vector, std::valarray or similar
// use like:
// std::valarray<double> a = linspace<std::valarray>(0.0, 10.0, 16);
//
template <typename T>
inline std::vector<T> linspace(T start, T stop, std::size_t size) {
  if (size == 0) { throw std::invalid_argument("linspace: argument size == 0"); }
  if (size == 1) { return {stop}; }

  std::vector<T> v(size);
  for (std::size_t i = 0; i < size; ++i) { v[i] = start + i * (stop - start) / (size - 1); }
  return v;
}

template <template <typename T> class container_type, class value_type>
inline container_type<value_type> linspace(value_type start, value_type stop, std::size_t size) {
  if (size == 0) { throw std::invalid_argument("linspace: argument size == 0"); }
  if (size == 1) { return {stop}; }

  container_type<value_type> v(size);
  for (std::size_t i = 0; i < size; ++i) { v[i] = start + i * (stop - start) / (size - 1); }
  return v;
}
template <template <typename T> class container_type, class value_type>
inline container_type<value_type> arange(value_type start, value_type step, value_type stop) {
  std::size_t                size = (stop - start) / step;
  container_type<value_type> v(size);
  for (std::size_t i = 0; i < size; ++i) { v[i] = start + step * i; }
  return v;
}

// Cumulative sum
template <typename T>
inline void CumSum(const std::vector<T> &x, std::vector<T> &sumvec) {
  const std::size_t N = x.size();
  sumvec.resize(N, 0.0);
  sumvec[0] = x[0];  // First
  for (std::size_t i = 1; i < N; ++i) { sumvec[i] = sumvec[i - 1] + x[i]; }
}

// Explicit powers for faster evaluation (std::pow is slow)
template <typename T>
constexpr T pow2(T x) {
  return x * x;
}
template <typename T>
constexpr T pow3(T x) {
  return x * pow2(x);
}
template <typename T>
constexpr T pow4(T x) {
  return pow2(pow2(x));
}
template <typename T>
constexpr T pow5(T x) {
  return x * pow4(x);
}
template <typename T>
constexpr T pow6(T x) {
  return pow2(pow3(x));
}
template <typename T>
constexpr T pow7(T x) {
  return x * pow6(x);
}
template <typename T>
constexpr T pow8(T x) {
  return pow2(pow4(x));
}

// Energy-Momentum conservation
inline bool CheckEMC(const M4Vec &diff, double epsilon = 1e-6) {
  if (std::abs(diff.Px()) > epsilon || std::abs(diff.Py()) > epsilon ||
      std::abs(diff.Pz()) > epsilon || diff.E() > epsilon) {
    return false;
  } else {
    return true;
  }
}

// Factorial n! = n*(n-1)*(n-2)*...*1
constexpr double factorial(int n) {
  double x = 1;
  for (int i = 1; i <= n; ++i) { x *= i; }
  return x;
}

// Complex sign function
constexpr double csgn(std::complex<double> z) {
  const double Re = std::real(z);
  const double Im = std::imag(z);

  if (Re > 0.0 || (std::abs(Re) < 1e-12 && Im > 0.0)) {
    return 1.0;
  } else {
    return -1.0;
  }
}

// Sign function
//
// Sign for x=0 is 0
//
constexpr double sign(double x) {
  const double EPS = 1e-9;
  if (std::abs(x) < EPS) {
    return 0.0;
  } else if (x > 0) {
    return 1.0;
  } else {
    return -1.0;
  }
}

// Sgn function
constexpr double sign(double X, double Y) {
  if (Y < 0.0) {
    return (-std::abs(X));
  } else {
    return (std::abs(X));
  }
}

// Sign function
// template <typename T>
// constexpr int sgn(T val) {
//  return (T(0) < val) - (val < T(0));
//}

// Floating point rounding safe square root
inline double msqrt(double x) { return std::sqrt(std::max(0.0, x)); }

// Complex amplitude squared
inline double abs2(const std::complex<double> &M) { return pow2(std::abs(M)); }

// ||v||^2
template <typename T>
inline double vpow2(const std::vector<T> &v) {
  double sum = 0.0;
  for (std::size_t i = 0; i < v.size(); ++i) { sum += abs2(v[i]); }
  return sum;
}

// ||v||
template <typename T>
inline double l2norm(const std::vector<T> &v) {
  double sum = vpow2(v);
  return msqrt(sum);
}

// Template print
template <template <typename T> class container_type, class value_type>
inline void PrintArray(container_type<value_type> x, std::string name) {
  std::cout << "PrintArray: " << name << std::endl;
  for (unsigned int i = 0; i < x.size(); ++i) { std::cout << x[i]; }
  std::cout << std::endl << std::endl;
}

// Degrees to radians
constexpr double Deg2Rad(double deg) { return (deg / 180.0) * gra::math::PI; }

// Radians to degrees
constexpr double Rad2Deg(double rad) { return (rad * 180.0) / gra::math::PI; }

// Binomial Coefficient C(n, k)
constexpr int Cbinom(int n, int k) {
  if (k == 0 || k == n) { return 1; }

  // Recursive function
  return Cbinom(n - 1, k - 1) + Cbinom(n - 1, k);
}

// N-dim epsilon tensor e_{\mu_1,\mu_2,\mu_3,...,\mu_N}:
//  + 1 if even permutation of arguments
//  - 1 if odd permutation of arguments
//    0 otherwise
inline MTensor<int> EpsTensor(std::size_t N) {
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

// Composite trapezoidal integral:
// Has geometric (=fast) convergence for periodic functions:
// - a circle in complex plane
// - interval on real line
// - Hankel contour around (-inf,0]
//
// hstep = (b-a)/N is the discretization size
//
// Indexing: j = 0,1,...,N-1,N,  => length of f = N+1, N = even
//
template <typename T>
inline T CTrapIntegral(const std::vector<T> &f, double hstep) {
  if ((f.size() - 1) % 2 != 0) {  // Must be even
    std::string str = "FATAL ERROR: gra::math::CTrapIntegral N = " + std::to_string(f.size() - 1) +
                      " is not even!";
    throw std::invalid_argument(str);
  }

  const std::size_t N     = f.size() - 1;
  T                 I_sum = 0.0;

  for (std::size_t j = 1; j <= N - 1; ++j) { I_sum += f[j]; }
  T I = hstep * (f[0] + 2 * I_sum + f[N]) / 2;

  return I;
}

// Composite Simpson's rule integral (quadratic interpolation),
// - hstep=(b-a)/N the discretization size
//
// Indexing: j = 0,1,...,N-1,N,  => length of f = N+1, N = even
//
template <typename T>
inline T CSIntegral(const std::vector<T> &f, double hstep) {
  if ((f.size() - 1) % 2 != 0) {  // Must be even
    std::string str =
        "FATAL ERROR: gra::math::CSIntegral N = " + std::to_string(f.size() - 1) + " is not even!";
    throw std::invalid_argument(str);
  }

  const std::size_t N      = f.size() - 1;
  T                 I_even = 0.0;
  T                 I_odd  = 0.0;

  for (std::size_t j = 1; j <= N / 2 - 1; ++j) { I_even += f[2 * j]; }
  for (std::size_t j = 1; j <= N / 2; ++j) { I_odd += f[2 * j - 1]; }
  T I = hstep * (f[0] + 2.0 * I_even + 4.0 * I_odd + f[N]) / 3.0;

  return I;
}

// Composite Simpson's 3/8 rule integral (qubic interpolation),
// - hstep=(b-a)/3N the discretization size
//
// Indexing: j = 0,1,...,3N, => length of f = 3N+1
//
template <typename T>
inline T CS38Integral(const std::vector<T> &f, double hstep) {
  if ((f.size() - 1) % 3 != 0) {  // Must be multiple of 3
    std::string str = "FATAL ERROR: gra::math::CS38Integral N = " +
                      std::to_string((f.size() - 1) / 3) + " is not multiple of 3!";
    throw std::invalid_argument(str);
  }

  const std::size_t N = (f.size() - 1) / 3;

  T I_first  = 0.0;
  T I_second = 0.0;

  for (std::size_t j = 1; j <= N; ++j) { I_first += (f[3 * j - 2] + f[3 * j - 1]); }
  for (std::size_t j = 1; j <= N - 1; ++j) { I_second += f[3 * j]; }
  T I = 3 * hstep * (f[0] + 3 * I_first + 2 * I_second + f[3 * N]) / 8;

  return I;
}

// ---------------------------------------------------------------------
// 2D-Simpson's 1/3 rule, see functions below
//
// f, W with dim[M+1,N+1] where M,N are even
// hstepM, hstepN are discretization step sizes (b-a)/M, (d-c)/N
//
// Based on Fubini's theorem <https://en.wikipedia.org/wiki/Fubini%27s_theorem>

template <typename T>
inline T Simpson13Integral2D(const MMatrix<T> &f, const MMatrix<double> &W, double hstepM,
                             double hstepN) {
  if (((f.size_row() - 1) % 2 != 0) || ((f.size_col() - 1) % 2 != 0)) {  // Must be even
    std::string str =
        "FATAL ERROR: gra::math::Simpson13Integral2D M,N must be "
        "even!";
    throw std::invalid_argument(str);
  }

  T I = 0.0;
  for (std::size_t i = 0; i < f.size_row(); ++i) {
    for (std::size_t j = 0; j < f.size_col(); ++j) { I += W[i][j] * f[i][j]; }
  }
  I *= (hstepM * hstepN) / 9;

  return I;
}

// ---------------------------------------------------------------------
// 2D-Simpson's 3/8 rule, see functions below
//
// f, W with dim[M+1,N+1] where M,N are multiple of 3
// hstepM, hstepN are discretization step sizes (b-a)/M, (d-c)/N
//
// Based on Fubini's theorem <https://en.wikipedia.org/wiki/Fubini%27s_theorem>

template <typename T>
inline T Simpson38Integral2D(const MMatrix<T> &f, const MMatrix<double> &W, double hstepM,
                             double hstepN) {
  if (((f.size_row() - 1) % 3 != 0) || ((f.size_col() - 1) % 3 != 0)) {
    std::string str =
        "FATAL ERROR: gra::math::Simpson38Integral2D M,N must be "
        "multiple of 3!";
    throw std::invalid_argument(str);
  }

  T I = 0.0;
  for (std::size_t i = 0; i < f.size_row(); ++i) {
    for (std::size_t j = 0; j < f.size_col(); ++j) { I += W[i][j] * f[i][j]; }
  }
  I *= 9 * (hstepM * hstepN) / 64;  // 3^2, 9^2
  return I;
}

// Example with M = 6, N = 8:

// 1    4    2    4    2    4    2    4    1
// 4   16    8   16    8   16    8   16    4
// 2    8    4    8    4    8    4    8    2
// 4   16    8   16    8   16    8   16    4
// 2    8    4    8    4    8    4    8    2
// 4   16    8   16    8   16    8   16    4
// 1    4    2    4    2    4    2    4    1

// Simpson's 1/3 rule weight vector dim[N+1]
inline std::vector<double> Simpson13Weight(unsigned int N) {
  if (N % 2 != 0) {  // Must be even
    std::string str =
        "FATAL ERROR: gra::math::SimpsonWeight N = " + std::to_string(N) + " is not even!";
    throw std::invalid_argument(str);
  }
  std::vector<double> W(N + 1, 2.0);  // Init with 2

  // Set 1 as the first and the last
  W[0] = 1.0;
  W[N] = 1.0;

  for (std::size_t i = 1; i < N; i += 2) { W[i] = 4.0; }
  return W;
}

// Simpson's 3/8 rule weight vector dim[N+1]
inline std::vector<double> Simpson38Weight(unsigned int N) {
  if (N % 3 != 0) {  // Must be multiple of 3
    std::string str =
        "FATAL ERROR: gra::math::SimpsonWeight N = " + std::to_string(N) + " is not multiple of 3!";
    throw std::invalid_argument(str);
  }
  std::vector<double> W(N + 1, 3.0);  // Init with 3

  // Set 1 as the first and the last
  W[0] = 1.0;
  W[N] = 1.0;

  for (std::size_t i = 3; i < N; i += 3) { W[i] = 2.0; }
  return W;
}

// Create a weight matrix dim[M+1, N+1] for 2D-Simpson's 1/3 rule
inline MMatrix<double> Simpson13Weight2D(unsigned int M, unsigned int N) {
  MMatrix<double>           W(M + 1, N + 1, 0.0);
  const std::vector<double> U = Simpson13Weight(M);
  const std::vector<double> V = Simpson13Weight(N);

  // Outer (tensor) product: W = U*V^T
  // printf("\n");
  for (std::size_t i = 0; i <= M; ++i) {
    for (std::size_t j = 0; j <= N; ++j) { W[i][j] = U[i] * V[j]; }
  }
  return W;
}

// Create a weight matrix dim[M+1, N+1] for 2D-Simpson's 3/8 rule
inline MMatrix<double> Simpson38Weight2D(unsigned int M, unsigned int N) {
  MMatrix<double>           W(M + 1, N + 1, 0.0);
  const std::vector<double> U = Simpson38Weight(M);
  const std::vector<double> V = Simpson38Weight(N);

  // Outer (tensor) product: W = U*V^T
  // printf("\n");
  for (std::size_t i = 0; i <= M; ++i) {
    for (std::size_t j = 0; j <= N; ++j) { W[i][j] = U[i] * V[j]; }
  }
  return W;
}

// Construct binary reflect Gray code (BRGC)
// which is the Hamiltonian Path on a N-dim unit hypercube (2^N Boolean vector
// space)
// The operator ^ is Exclusive OR (XOR) and the operator >> is bit shift right
constexpr unsigned int Binary2Gray(unsigned int number) { return number ^ (number >> 1); }

// Convert BRGC (binary reflect Gray code) to a binary number
// Gray code is obtained as XOR (^) with all more significant bits
constexpr unsigned int Gray2Binary(unsigned int number) {
  for (unsigned int mask = number >> 1; mask != 0; mask = mask >> 1) { number = number ^ mask; }
  return number;
}

// Boolean vector to index representation
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
inline int Vec2Ind(std::vector<bool> vec) {
  // Swap bit direction (comment this for the reverse ordering)
  std::reverse(vec.begin(), vec.end());

  int retval = 0;
  int i      = 0;

  for (std::vector<bool>::iterator it = vec.begin(); it != vec.end(); ++it, ++i) {
    if (*it) { retval |= 1 << i; }
  }
  return retval;
}

// Index to Boolean vector representation
// [the normal order: 0 ~ 00, 1 ~ 01, 2 ~ 10, 3 ~ 11 ...]
// Input: x = index representation
//        d = Boolean vector space dimension
inline std::vector<bool> Ind2Vec(unsigned int ind, unsigned int d) {
  std::vector<bool> ret;

  while (ind) {
    if (ind & 1)
      ret.push_back(1);
    else
      ret.push_back(0);

    ind >>= 1;
  }
  reverse(ret.begin(), ret.end());

  // Now will the d-dimensional bit vector starting from rightmost bit
  std::vector<bool> final(d, 0);
  for (std::size_t i = 0; i < ret.size(); ++i) { final[d - 1 - i] = ret[ret.size() - 1 - i]; }

  return final;
}

// OEIS.org A030109 sequence of left-right bit reversed sequence
// Input: dim = Boolean vector space dimension
inline std::vector<unsigned int> LRsequence(unsigned int dim) {
  // The sequence
  std::vector<unsigned int> seq(std::pow(2, dim), 0);

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
inline std::vector<std::vector<int>> BinaryMatrix(unsigned int d) {
  std::vector<std::vector<int>> B;
  unsigned int                  N = std::pow(2, d);

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

// Assumes big endian
// Example of how to use: int i = 6; PrintBits(sizeof(i), &i);
// >> 00000000000000000000000000000110
inline void PrintBits(size_t const size, void const *const ptr) {
  unsigned char *b = (unsigned char *)ptr;
  unsigned char  byte;

  // Loop over bytes, e.g. for int32 the size is 4
  for (int i = size - 1; i >= 0; i--) {
    // Loop over bits of a byte
    for (int j = 7; j >= 0; j--) {
      byte = (b[i] >> j) & 1;
      printf("%u", byte);
    }
  }
  puts("");
}


// Go through all the permutations, N! number of them
// Algorithm from the Knuth's "Art of Computer Programming"
//
// Input example:
//
// std::vector<int> set = {1,2,3};
inline std::vector<std::vector<int>> Permutations(std::vector<int> &input) {
  std::size_t N            = input.size();
  int         permutations = (int)factorial(N);

  // Here we save all permutations
  std::vector<std::vector<int>> outset(permutations, std::vector<int>(N));
  // This we modify
  std::vector<int> set = input;

  int temp = 0;
  int j    = 0;
  int l    = 0;

  // Save first
  outset[0] = set;
  int k     = 1;

  while (true) {
    // Find largest j for which set[j] < set[j+1]
    j = -1;
    for (std::size_t i = 0; i < N - 1; ++i) {
      if (set[i] < set[i + 1]) j = i;
    }
    if (j == -1) {
      break;  // No such j, we are done
    }
    for (std::size_t i = j + 1; i < N; ++i) {
      if (set[i] > set[j]) { l = i; }
    }
    temp   = set[j];
    set[j] = set[l];
    set[l] = temp;

    // Reverse from j+1 to the end
    unsigned int swaps = (N - j - 1) / 2;
    for (unsigned int i = 0; i < swaps; ++i) {
      if ((j + 1 + i) < N) {
        temp             = set[(j + 1) + i];
        set[(j + 1) + i] = set[N - 1 - i];
        set[N - 1 - i]   = temp;
      }
    }
    // Save permutation set
    outset[k] = set;
    ++k;
  }
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
    return outset;  // Return complete permutations
  } else {          // Return local charge conserving

    // Create alternating charge series
    // 1,-1,1,-1,...
    std::vector<int> charge;
    for (int i = 0; i < N; ++i) {
      charge.push_back(pow(-1, i));
      //  printf("%d", charge[i]);
    }
    // printf("\n");

    std::vector<std::vector<int>> valid_amplitudes;

    for (std::size_t i = 0; i < outset.size(); ++i) {
      bool ok = true;
      for (std::size_t j = 0; j < outset[0].size() - 1; j += 2) {
        int a = outset[i][j];
        int b = outset[i][j + 1];

        if ((charge.at(a - offset) + charge.at(b - offset)) != 0) ok = false;
        // printf("%d%d", a, b);
      }
      if (ok) { valid_amplitudes.push_back(outset[i]); }
      //  printf(" ok = %d", ok);
      //  printf("\n");
    }
    // printf("Total = %d, Charge conserving = %d \n",
    // outset.size(), valid_amplitudes.size());
    return valid_amplitudes;
  }
}

// Ordinary Legendre polynomials P_l(x) to cross check algorithmic
// implementation
// x usually cos(theta)
inline double LegendrePl(unsigned int l, double x) {
  if (l == 0) {
    return 1;
  } else if (l == 1) {
    return x;
  } else if (l == 2) {
    return 0.5 * (3 * x * x - 1);
  } else if (l == 3) {
    return 0.5 * (5 * x * x * x - 3 * x);
  } else if (l == 4) {
    return (35 * x * x * x * x - 30 * x * x + 3) / 8;
  } else if (l == 5) {
    return (63 * x * x * x * x * x - 70 * x * x * x + 15 * x) / 8;
  } else if (l == 6) {
    return (231 * x * x * x * x * x * x - 315 * x * x * x * x + 105 * x * x - 5) / 16;
  } else if (l == 7) {
    return (429 * x * x * x * x * x * x * x - 693 * x * x * x * x * x + 315 * x * x * x - 35 * x) /
           16;
  } else if (l == 8) {
    return (6435 * x * x * x * x * x * x * x * x - 12012 * x * x * x * x * x * x +
            6930 * x * x * x * x - 1260 * x * x + 35) /
           128;
  } else {
    throw std::invalid_argument("LegendrePl: Not valid l = " + std::to_string(l));
  }
}

// Associated Legendre Polynomial P_l^m(x) solved numerically
//
// Algorithm from:
//
// [REFERENCE: S. Jemma, The Nitty Gritty Details of Spherical Harmonics, SIGRAPH2003]
constexpr double sf_legendre(int l, int m, double x) {
  // Associated Legendre Polynomial
  double pmm = 1.0;
  if (m > 0) {
    const double sx2  = std::sqrt((1.0 - x) * (1.0 + x));
    double       fact = 1.0;
    for (int i = 1; i <= m; ++i) {
      pmm *= (-fact) * sx2;
      fact += 2.0;
    }
  }
  if (l == m) { return pmm; }

  double pmmp1 = x * (2.0 * m + 1.0) * pmm;
  if (l == m + 1) { return pmmp1; }

  double pll = 0.0;
  for (int ll = m + 2; ll <= l; ++ll) {
    pll   = ((2.0 * ll - 1.0) * x * pmmp1 - (ll + m - 1.0) * pmm) / (ll - m);
    pmm   = pmmp1;
    pmmp1 = pll;
  }
  return pll;
}

// -----------------------------------------------------------------------
// Associated Legendre Polynomials with Spherical Harmonics
//
//
// Quantum Mechanics convention
//  Y_l^m(theta,phi)
//  = (-1)^m \sqrt{\frac{2l+1}{4\pi} \frac{(l-m)!}{(l+m)!} P_{lm}(cos(theta))
//  e^{im\phi}
//
// Normalization:
//  d\Omega = sin(theta) d\phi d\theta
//  \int |Y_l^m|^2 d\Omega = 1
//
// Orthogonality:
//  \int_{\theta = 0}^\pi \int_{\phi=0}^{2\pi} Y_l^m Y_l^{m'}^* d\Omega =
//  \delta_{ll'} \delta_{mm'}
//
//
// Relation:
//  Y_l^{m*}(theta,phi) = (-1)^m Y_l^{-m} (theta,phi)
//
// When m = 0, we get ordinary Legendre polynomials
//
// Parity of spherical harmonics:
//       {r, \theta, \phi} -> {r, \pi - \theta, \pi + \phi}
//
//
// Y_l^m(-\vec{r}) = (-1)^l Y_l^m(\vec{r})
// Y_l^m(theta,phi) -> Y_l^m (\pi - \theta, \pi + \phi) = (-1)^l Y_l^m(\theta,
// \phi)
//
//
// Parity = (-1)^l, so parity is even or odd, depending on the quantum number l

// Real valued representation of complex spherical harmonics
//
// Normalization such that \int dOmega |Y_l^m(Omega)|^2 = 1
inline double NReY(const std::complex<double> &Y, int l, int m) {
  if (m < 0) {
    return sqrt(2.0) * std::pow(-1, m) * std::imag(Y);
  } else if (m == 0) {
    return std::real(Y);
  } else {
    return sqrt(2.0) * std::pow(-1, m) * std::real(Y);
  }
}

// Real Basis Laplace Spherical Harmonics
//
// [REFERENCE: https://en.wikipedia.org/wiki/Spherical_harmonics]
inline double Y_real_basis(double costheta, double phi, int l, int m) {
  // Normalization function
  auto K = [](int l, int m) -> double {

    if (m == 0) { return msqrt((2.0 * l + 1.0) / (4.0 * PI)); }  // Speed it up
    double temp =
        ((2.0 * l + 1.0) * factorial(l - m)) / static_cast<double>(4.0 * PI * factorial(l + m));
    return msqrt(temp);
  };
  const double sqrt2 = msqrt(2.0);

  if (m < 0) {
    return std::pow(-1, m) * sqrt2 * K(l, -m) * std::sin(-m * phi) * sf_legendre(l, -m, costheta);
  } else if (m == 0) {
    return K(l, 0) * sf_legendre(l, m, costheta);
  } else {  // m > 0
    return std::pow(-1, m) * sqrt2 * K(l, m) * std::cos(m * phi) * sf_legendre(l, m, costheta);
  }
}

// Complex Basis Laplace Spherical Harmonics
inline std::complex<double> Y_complex_basis(double costheta, double phi, int l, int m) {
  if (m < 0) {
    return 1.0 / msqrt(2.0) *
           (Y_real_basis(costheta, phi, l, std::abs(m)) -
            zi * Y_real_basis(costheta, phi, l, -std::abs(m)));
  } else if (m == 0) {
    return Y_real_basis(costheta, phi, l, 0);
  } else {  // m > 0
    return std::pow(-1, m) / msqrt(2.0) *
           (Y_real_basis(costheta, phi, l, std::abs(m)) +
            zi * Y_real_basis(costheta, phi, l, -std::abs(m)));
  }
}

// Manually written Spherical Harmonics to cross-check algorithmic
// implementations
inline std::complex<double> Y_complex_basis_ref(double costheta, double phi, int l, int m) {
  const double theta = std::acos(costheta);

  if (l == 0 && m == 0) return 0.5 * std::sqrt(1.0 / PI);

  // ------------------------------------------------------------------
  if (l == 1 && m == -1)
    return 0.5 * std::sqrt(3.0 / (2.0 * PI)) * std::exp(-gra::math::zi * phi) * std::sin(theta);

  if (l == 1 && m == 0) return 0.5 * std::sqrt(3.0 / PI) * costheta;

  if (l == 1 && m == 1)
    return -0.5 * std::sqrt(3.0 / (2.0 * PI)) * std::exp(gra::math::zi * phi) * std::sin(theta);

  // ------------------------------------------------------------------
  if (l == 2 && m == -2)
    return 0.25 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2);

  if (l == 2 && m == -1)
    return 0.5 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(-gra::math::zi * phi) * std::sin(theta) *
           costheta;

  if (l == 2 && m == 0) return 0.25 * std::sqrt(5.0 / PI) * (3.0 * std::pow(costheta, 2) - 1.0);

  if (l == 2 && m == 1)
    return -0.5 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(gra::math::zi * phi) * std::sin(theta) *
           costheta;

  if (l == 2 && m == 2)
    return 0.25 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2);

  // ------------------------------------------------------------------

  if (l == 3 && m == -3)
    return 1.0 / 8.0 * std::sqrt(35.0 / PI) * std::exp(-3.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 3);

  if (l == 3 && m == -2)
    return 0.25 * std::sqrt(105.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * costheta;

  if (l == 3 && m == -1)
    return 1.0 / 8.0 * std::sqrt(21.0 / PI) * std::exp(-gra::math::zi * phi) * std::sin(theta) *
           (5.0 * std::pow(costheta, 2) - 1.0);

  if (l == 3 && m == 0)
    return 0.25 * std::sqrt(7.0 / PI) * (5.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 3 && m == 1)
    return -1.0 / 8.0 * std::sqrt(21.0 / PI) * std::exp(gra::math::zi * phi) * std::sin(theta) *
           (5.0 * std::pow(costheta, 2) - 1.0);

  if (l == 3 && m == 2)
    return 0.25 * std::sqrt(105.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * costheta;

  if (l == 3 && m == 3)
    return -1.0 / 8.0 * std::sqrt(35.0 / PI) * std::exp(3.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 3);

  // ------------------------------------------------------------------

  if (l == 4 && m == -4)
    return 3.0 / 16.0 * std::sqrt(35.0 / (2.0 * PI)) * std::exp(-4.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 4);

  if (l == 4 && m == -3)
    return 3.0 / 8.0 * std::sqrt(35.0 / (PI)) * std::exp(-3.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 3) * costheta;

  if (l == 4 && m == -2)
    return 3.0 / 8.0 * std::sqrt(5.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * (7.0 * std::pow(costheta, 2) - 1.0);

  if (l == 4 && m == -1)
    return 3.0 / 8.0 * std::sqrt(5.0 / PI) * std::exp(-gra::math::zi * phi) * std::sin(theta) *
           (7.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 4 && m == 0)
    return 3.0 / 16.0 * std::sqrt(1.0 / PI) *
           (35.0 * std::pow(costheta, 4) - 30.0 * std::pow(costheta, 2) + 3.0);

  if (l == 4 && m == 1)
    return -3.0 / 8.0 * std::sqrt(5.0 / PI) * std::exp(gra::math::zi * phi) * std::sin(theta) *
           (7.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 4 && m == 2)
    return 3.0 / 8.0 * std::sqrt(5.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * (7.0 * std::pow(costheta, 2) - 1.0);

  if (l == 4 && m == 3)
    return -3.0 / 8.0 * std::sqrt(35.0 / (PI)) * std::exp(3.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 3) * costheta;

  if (l == 4 && m == 4)
    return 3.0 / 16.0 * std::sqrt(35.0 / (2.0 * PI)) * std::exp(4.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 4);

  throw std::invalid_argument("Y_complex_basis_ref:: Not supported l = " + std::to_string(l) +
                              ", m = " + std::to_string(m));
}

}  // namespace math
}  // namespace gra

#endif