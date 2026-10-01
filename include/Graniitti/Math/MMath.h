// General mathematical functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATH_H
#define MMATH_H

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <type_traits>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MFloat.h"

namespace gra {
namespace math {

// Imaginary unit
constexpr const std::complex<double> zi(0.0, 1.0);

// Geometry
constexpr const double PI   = 3.141592653589793238462643383279502884197169399375105820974944L;
constexpr const double PIPI = PI * PI;
constexpr const double ZETA2 = PIPI / 6.0;

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

// Compute one nonnegative integer power by repeated multiplication
// value = base^exponent by binary exponentiation
template <typename T>
constexpr T IntegerPower(T base, unsigned int exponent) {
  T value = T{1};
  while (exponent > 0) {
    if ((exponent & 1U) != 0U) { value *= base; }
    exponent >>= 1U;
    if (exponent > 0) { base *= base; }
  }
  return value;
}

// Compute all nonnegative integer powers through one maximum exponent
// powers_n = base^n, n = 0,...,maximum
template <typename T>
std::vector<T> IntegerPowers(const T base, const unsigned int maximum) {
  std::vector<T> powers(static_cast<std::size_t>(maximum) + 1, T{1});
  for (unsigned int exponent = 1; exponent <= maximum; ++exponent) { powers[exponent] = powers[exponent - 1] * base; }
  return powers;
}

// Check energy momentum conservation within one absolute tolerance
// max_mu |Delta p^mu| <= epsilon
inline bool CheckEMC(const M4Vec &diff, double epsilon = 1e-6) {
  if (!std::isfinite(epsilon) || epsilon < 0.0 || !std::isfinite(diff.Px()) || !std::isfinite(diff.Py()) ||
      !std::isfinite(diff.Pz()) || !std::isfinite(diff.E())) {
    return false;
  }
  if (std::abs(diff.Px()) > epsilon || std::abs(diff.Py()) > epsilon || std::abs(diff.Pz()) > epsilon ||
      std::abs(diff.E()) > epsilon) {
    return false;
  } else {
    return true;
  }
}

// Compute the exact arithmetic sign with zero mapped to zero
// sign(x) = 1 for x > 0, 0 for x = 0, and -1 for x < 0
template <typename T>
requires std::is_arithmetic_v<T>
constexpr int sign(T value) { return (T{0} < value) - (value < T{0}); }

// Floating point rounding safe square root
// msqrt(x) = sqrt(max(0,x))
inline double msqrt(double x) { return std::sqrt(std::max(x, 0.0)); }

// Compute the squared magnitude of one real or complex scalar
// abs2(x) = x^2 for real x and x x* for complex x
template <typename T>
inline auto abs2(const T &value) {
  if constexpr (requires {
                  value.real();
                  value.imag();
                }) {
    return std::norm(value);
  } else {
    return value * value;
  }
}

// Degrees to radians
// rad = pi deg/180
constexpr double Deg2Rad(double deg) { return (deg / 180.0) * gra::math::PI; }

// Radians to degrees
// deg = 180 rad/pi
constexpr double Rad2Deg(double rad) { return (rad * 180.0) / gra::math::PI; }

// Wrap one angle to the interval (-pi, pi]
// output = angle mod 2pi in (-pi,pi]
inline double WrapAngle(double angle) {
  constexpr double period = 2.0 * PI;
  double           output = std::fmod(angle + PI, period);
  if (output < 0.0) { output += period; }
  output -= PI;
  return output <= -PI ? PI : output;
}

}  // namespace math
}  // namespace gra

#endif
