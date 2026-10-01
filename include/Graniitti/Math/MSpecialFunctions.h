// Special functions and controlled series approximations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSPECIALFUNCTIONS_H
#define MSPECIALFUNCTIONS_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace math {

// Compute log I_1(x) without exponential overflow for positive x
double LogBesselI1(double x);

// Bound the absolute exponential Taylor tail after a fixed term count
double ExponentialTailBound(double magnitude, std::size_t terms);

// Compute the first-kind Bessel function of integer order zero
double BesselJ0(double x);

// Compute the first-kind Bessel function of integer order one
double BesselJ1(double x);

// Compute the first three first-kind integer-order Bessel functions
std::array<double, 3> BesselJ012(double x);

// Compute the first-kind Bessel function of one nonnegative integer order
double BesselJ(int n, double x);

// Factorial n! = n*(n-1)*(n-2)*...*1
// n! = Gamma(n+1)
inline double factorial(int n) {
  constexpr double f[] = {1.,
                          1.,
                          2.,
                          6.,
                          24.,
                          120.,
                          720.,
                          5040.,
                          40320.,
                          362880.,
                          3628800.,
                          39916800.,
                          479001600.,
                          6227020800.,
                          87178291200.,
                          1307674368000.,
                          20922789888000.,
                          355687428096000.,
                          6402373705728000.,
                          121645100408832000.,
                          2432902008176640000.,
                          51090942171709440000.,
                          1124000727777607680000.,
                          25852016738884976640000.,
                          620448401733239439360000.,
                          15511210043330985984000000.,
                          403291461126605635584000000.,
                          10888869450418352160768000000.,
                          304888344611713860501504000000.,
                          8841761993739701954543616000000.,
                          265252859812191058636308480000000.,
                          8222838654177922817725562880000000.,
                          263130836933693530167218012160000000.,
                          8683317618811886495518194401280000000.,
                          295232799039604140847618609643520000000.,
                          10333147966386144929666651337523200000000.,
                          371993326789901217467999448150835200000000.,
                          13763753091226345046315979581580902400000000.,
                          523022617466601111760007224100074291200000000.,
                          20397882081197443358640281739902897356800000000.,
                          815915283247897734345611269596115894272000000000.};
  if (n < 0) { return 0.0; }
  if (n <= 40) { return f[n]; }
  if (n > 170) { return std::numeric_limits<double>::infinity(); }

  double x = 1;
  for (int i = 1; i <= n; ++i) { x *= i; }
  return x;
}

// Narrow a special function value without hiding non-finite results
inline double FiniteSpecialValue(const long double value) {
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) { throw gra::AmplitudeFailure("Special function: non-finite value"); }
  return result;
}

// Compute log|Gamma(x)| and optionally its sign using storage local to this call
// lgamma writes global signgam state, so concurrent sampling requires the reentrant variants
template <typename T>
inline T LogGamma(const T x, int *sign = nullptr) {
  static_assert(std::is_floating_point_v<T>);
  int local_sign = 1;
  T   value;
  if constexpr (std::is_same_v<T, long double>) {
    value = ::lgammal_r(x, &local_sign);
  } else if constexpr (std::is_same_v<T, float>) {
    value = ::lgammaf_r(x, &local_sign);
  } else {
    value = ::lgamma_r(x, &local_sign);
  }
  if (sign != nullptr) { *sign = local_sign; }
  return value;
}

// Compute Gamma products in logarithms, with exact zeros from denominator poles
// R = product_i Gamma(a_i) / product_j Gamma(b_j)
template <std::size_t N, std::size_t M>
inline long double GammaRatio(const std::array<double, N> &numerator, const std::array<double, M> &denominator) {
  long double logarithm = 0.0L;
  int         sign      = 1;
  // Accumulate a finite Gamma factor using its real sign and logarithmic magnitude
  const auto accumulate = [&](const double x, const int power) {
    if (!std::isfinite(x)) { throw gra::AmplitudeFailure("GammaRatio: non-finite argument"); }
    if (x <= 0.0 && math::IsExactEqual(x, std::round(x))) { return false; }
    // Keep the Gamma sign local to each call because lgamma writes shared signgam state during sampling
    int gamma_sign = 1;
    logarithm += power * LogGamma(static_cast<long double>(x), &gamma_sign);
    sign *= gamma_sign;
    return true;
  };
  for (const double x : denominator) {
    if (!accumulate(x, -1)) { return 0.0; }
  }
  for (const double x : numerator) {
    if (!accumulate(x, 1)) { throw gra::AmplitudeFailure("GammaRatio: numerator pole"); }
  }
  const long double result = sign * std::exp(logarithm);
  if (!std::isfinite(result)) { throw gra::AmplitudeFailure("GammaRatio: non-finite value"); }
  return result;
}

// Compute reciprocal Gamma with exact zeros at nonpositive integer poles
// result = 1/Gamma(x)
inline double ReciprocalGamma(const double x) {
  return FiniteSpecialValue(GammaRatio(std::array<double, 0>{}, std::array{x}));
}

// Bound the remaining Gauss series once all shifted parameters are positive
// |R_n| <= |t_n| r/(1-r), with r bounding every subsequent term ratio
inline long double Hyper2F1Tail(const double a, const double b, const double c, const double z,
                                const int n, const long double term) {
  const long double an = static_cast<long double>(a) + n;
  const long double bn = static_cast<long double>(b) + n;
  const long double cn = static_cast<long double>(c) + n;
  if (!(an > 0.0L && bn > 0.0L && cn > 0.0L)) { return std::numeric_limits<long double>::infinity(); }
  const long double ratio = z * std::max(1.0L, an / cn) * std::max(1.0L, bn / (n + 1.0L));
  return ratio < 1.0L ? std::abs(term) * ratio / (1.0L - ratio) : std::numeric_limits<long double>::infinity();
}

// Continue a Gauss solution from z=1/2 using convergent local Taylor series
// z(1-z)F'' + [c-(a+b+1)z]F' - abF = 0
// [REFERENCE: NIST DLMF 15.10.1, https://dlmf.nist.gov/15.10.E1]
inline double Hyper2F1Continue(const long double a, const long double b, const long double c, const double z,
                               long double value, long double derivative) {
  if (!(z > 0.5 && z < 1.0) || !std::isfinite(value) || !std::isfinite(derivative)) {
    throw gra::AmplitudeFailure("Hyper2F1: invalid continuation state");
  }
  constexpr long double tolerance = 8.0L * std::numeric_limits<long double>::epsilon();
  long double           x         = 0.5L;
  while (x < z) {
    const long double step      = std::min(static_cast<long double>(z) - x, (1.0L - x) / 2.0L);
    const long double quadratic = x * (1.0L - x);
    const long double linear    = c - (a + b + 1.0L) * x;
    long double       previous  = value;
    long double       term      = step * derivative;
    long double       sum       = previous + term;
    long double       slope     = term;
    bool              converged = false;
    for (int n = 0; n < 512; ++n) {
      const long double next =
          ((n + a) * (n + b) * step * step * previous - (n + 1.0L) * (n * (1.0L - 2.0L * x) + linear) * step * term) /
          (quadratic * (n + 2.0L) * (n + 1.0L));
      sum += next;
      slope += (n + 2.0L) * next;
      if (n >= 16 && std::abs(term) + std::abs(next) <= tolerance * std::abs(sum) &&
          (n + 2.0L) * (std::abs(term) + std::abs(next)) <= tolerance * std::max(std::abs(sum), std::abs(slope))) {
        converged = true;
        break;
      }
      previous = term;
      term     = next;
    }
    if (!converged || !std::isfinite(sum) || !std::isfinite(slope)) {
      throw gra::AmplitudeFailure("Hyper2F1: continuation did not converge");
    }
    value      = sum;
    derivative = slope / step;
    x += step;
  }
  const double result = static_cast<double>(value);
  if (!std::isfinite(result)) { throw gra::AmplitudeFailure("Hyper2F1: continuation overflow"); }
  return result;
}

// Evaluate the Gauss hypergeometric function on the physical unit interval
// 2F1(a,b;c;z) = sum_k (a)_k (b)_k z^k / [(c)_k k!]
inline double Hyper2F1(const double a, const double b, const double c, const double z) {
  constexpr int         max_terms = 20000;
  constexpr long double tolerance = 1.0e-15L;
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(z) || z < 0.0 || z > 1.0 ||
      (c <= 0.0 && math::IsExactEqual(c, std::round(c)))) {
    throw gra::AmplitudeFailure("Hyper2F1: invalid argument or denominator pole");
  }
  // Evaluate the exact power law without cancellation in a decaying solution
  if (z < 1.0 && (math::IsExactEqual(a, c) || math::IsExactEqual(b, c))) {
    const long double power = math::IsExactEqual(a, c) ? b : a;
    return FiniteSpecialValue(std::exp(-power * std::log1p(-static_cast<long double>(z))));
  }
  // Identify terminating polynomials before applying the unity convergence condition
  const auto terminal_degree = [](const double value) {
    if (value > 0.0 || value <= -max_terms || !math::IsExactEqual(value, std::round(value))) { return max_terms; }
    return static_cast<int>(-std::llround(value));
  };
  const int terminal = std::min(terminal_degree(a), terminal_degree(b));
  if (terminal == max_terms && math::IsExactEqual(z, 1.0)) {
    const double gap = c - a - b;
    if (!(gap > 0.0)) { throw gra::AmplitudeFailure("Hyper2F1: divergent series at unity"); }
    return FiniteSpecialValue(GammaRatio(std::array{c, gap}, std::array{c - a, c - b}));
  }

  // Factor the decaying solution when Euler's transformation gives a polynomial
  // [REFERENCE: NIST DLMF 15.8.1, https://dlmf.nist.gov/15.8.E1]
  if (terminal == max_terms && z > 0.5 && z < 1.0 &&
      std::min(terminal_degree(c - a), terminal_degree(c - b)) < max_terms) {
    const long double factor =
        std::exp((static_cast<long double>(c) - a - b) * std::log1p(-static_cast<long double>(z)));
    return FiniteSpecialValue(factor * Hyper2F1(c - a, c - b, c, z));
  }

  if (terminal == max_terms && z > 0.5) {
    return Hyper2F1Continue(a, b, c, z, Hyper2F1(a, b, c, 0.5),
                            static_cast<long double>(a) * b / c * Hyper2F1(a + 1.0, b + 1.0, c + 1.0, 0.5));
  }
  const int   last = terminal < max_terms ? terminal : max_terms;
  long double term = 1.0L;
  long double sum  = 1.0L;
  for (int k = 0; k < last; ++k) {
    const long double denominator = (static_cast<long double>(c) + k) * (k + 1.0L);
    if (math::IsZero(denominator)) { throw gra::AmplitudeFailure("Hyper2F1: denominator pole"); }
    term *= (static_cast<long double>(a) + k) * (static_cast<long double>(b) + k) * static_cast<long double>(z) /
            denominator;
    sum += term;
    if (terminal == max_terms && Hyper2F1Tail(a, b, c, z, k + 1, term) <= tolerance * std::abs(sum)) {
      return FiniteSpecialValue(sum);
    }
  }
  if (terminal < max_terms) { return FiniteSpecialValue(sum); }
  throw gra::AmplitudeFailure("Hyper2F1: series did not converge");
}

// Evaluate the regularized Gauss hypergeometric function
// 2F1_reg(a,b;c;z) = sum_k (a)_k (b)_k z^k / [Gamma(c+k) k!]
inline double RegularizedHyper2F1(const double a, const double b, const double c, const double z) {
  constexpr int         max_terms = 20000;
  constexpr long double tolerance = 1.0e-15L;
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(c) || !std::isfinite(z) || z < 0.0 || z > 1.0) {
    throw gra::AmplitudeFailure("RegularizedHyper2F1: invalid argument");
  }
  // Evaluate the exact power law with the regularizing Gamma factor
  if (z < 1.0 && (math::IsExactEqual(a, c) || math::IsExactEqual(b, c))) {
    const long double power  = math::IsExactEqual(a, c) ? b : a;
    const long double factor = GammaRatio(std::array<double, 0>{}, std::array{c});
    if (math::IsZero(factor)) { return 0.0; }
    return FiniteSpecialValue(factor * std::exp(-power * std::log1p(-static_cast<long double>(z))));
  }
  // Identify terminating polynomials before applying the unity convergence condition
  const auto terminal_degree = [](const double value) {
    if (value > 0.0 || value <= -max_terms || !math::IsExactEqual(value, std::round(value))) { return max_terms; }
    return static_cast<int>(-std::llround(value));
  };
  const int terminal = std::min(terminal_degree(a), terminal_degree(b));
  if (terminal == max_terms && math::IsExactEqual(z, 1.0)) {
    const double gap = c - a - b;
    if (!(gap > 0.0)) { throw gra::AmplitudeFailure("RegularizedHyper2F1: divergent series at unity"); }
    return FiniteSpecialValue(GammaRatio(std::array{gap}, std::array{c - a, c - b}));
  }

  // Factor the decaying solution when Euler's transformation gives a polynomial
  // [REFERENCE: NIST DLMF 15.8.1, https://dlmf.nist.gov/15.8.E1]
  if (terminal == max_terms && z > 0.5 && z < 1.0 &&
      std::min(terminal_degree(c - a), terminal_degree(c - b)) < max_terms) {
    const long double factor =
        std::exp((static_cast<long double>(c) - a - b) * std::log1p(-static_cast<long double>(z)));
    return FiniteSpecialValue(factor * RegularizedHyper2F1(c - a, c - b, c, z));
  }

  if (terminal == max_terms && z > 0.5) {
    return Hyper2F1Continue(a, b, c, z, RegularizedHyper2F1(a, b, c, 0.5),
                            static_cast<long double>(a) * b * RegularizedHyper2F1(a + 1.0, b + 1.0, c + 1.0, 0.5));
  }
  int first = 0;
  if (c <= 0.0 && math::IsExactEqual(c, std::round(c))) {
    first = c <= -max_terms ? max_terms : 1 - static_cast<int>(std::llround(c));
  }
  if (terminal < first) { return 0.0; }
  if (first == max_terms) { throw gra::AmplitudeFailure("RegularizedHyper2F1: series starts beyond term limit"); }

  long double numerator       = 1.0L;
  long double factorial_value = 1.0L;
  long double z_power         = 1.0L;
  for (int k = 0; k < first; ++k) {
    numerator *= (static_cast<long double>(a) + k) * (static_cast<long double>(b) + k);
    factorial_value *= static_cast<long double>(k + 1);
    z_power *= static_cast<long double>(z);
  }
  long double term = numerator * z_power * GammaRatio(std::array<double, 0>{}, std::array{c + first}) / factorial_value;
  long double sum  = term;
  const int   last = terminal < max_terms ? terminal : max_terms;
  for (int k = first; k < last; ++k) {
    const long double denominator = (static_cast<long double>(c) + k) * (k + 1.0L);
    if (math::IsZero(denominator)) { throw gra::AmplitudeFailure("RegularizedHyper2F1: denominator pole"); }
    term *= (static_cast<long double>(a) + k) * (static_cast<long double>(b) + k) * static_cast<long double>(z) /
            denominator;
    sum += term;
    if (terminal == max_terms && Hyper2F1Tail(a, b, c, z, k + 1, term) <= tolerance * std::abs(sum)) {
      return FiniteSpecialValue(sum);
    }
  }
  if (terminal < max_terms) { return FiniteSpecialValue(sum); }
  throw gra::AmplitudeFailure("RegularizedHyper2F1: series did not converge");
}

// Bound the remaining 3F2 series by ratios r_k <= k/(k+p), with p = 1 + gap/2
// Nonnegative coefficients after k = n + u prove the bound for every u >= 0
inline long double Hyper3F2Tail(const std::array<double, 3> &a, const std::array<double, 2> &b,
                                const int n, const long double term) {
  const long double gap = static_cast<long double>(b[0]) + b[1] - a[0] - a[1] - a[2];
  if (!(gap > 0.0L) || n <= std::max({-a[0], -a[1], -a[2], -b[0], -b[1]})) {
    return std::numeric_limits<long double>::infinity();
  }
  const long double p  = 1.0L + gap / 2.0L;
  const long double a2 = static_cast<long double>(a[0]) + a[1] + a[2];
  const long double a1 = static_cast<long double>(a[0]) * a[1] + static_cast<long double>(a[0]) * a[2] +
                         static_cast<long double>(a[1]) * a[2];
  const long double a0 = static_cast<long double>(a[0]) * a[1] * a[2];
  const long double c3 = gap / 2.0L;
  const long double c2 = static_cast<long double>(b[0]) * b[1] + b[0] + b[1] - a1 - p * a2;
  const long double c1 = static_cast<long double>(b[0]) * b[1] - a0 - p * a1;
  const long double c0 = -p * a0;
  const long double x  = n;
  if (c2 + 3.0L * x * c3 < 0.0L || c1 + x * (2.0L * c2 + 3.0L * x * c3) < 0.0L ||
      c0 + x * (c1 + x * (c2 + x * c3)) < 0.0L) {
    return std::numeric_limits<long double>::infinity();
  }
  return std::abs(term) * (2.0L * x / gap);
}

// Estimate the 3F2 tail by t_n (n+c)/gap and bound its recurrence residual
// Summing t_k [r_k (1+R(k+1))-R(k)] bounds the error of R(n) t_n
inline std::pair<long double, long double> Hyper3F2Remainder(const std::array<double, 3> &a,
                                                           const std::array<double, 2> &b,
                                                           const int n, const long double term) {
  const long double tail = Hyper3F2Tail(a, b, n, term);
  if (!std::isfinite(tail)) { return {0.0L, tail}; }
  const long double gap = static_cast<long double>(b[0]) + b[1] - a[0] - a[1] - a[2];
  const long double a2 = static_cast<long double>(a[0]) + a[1] + a[2];
  const long double a1 = static_cast<long double>(a[0]) * a[1] + static_cast<long double>(a[0]) * a[2] +
                         static_cast<long double>(a[1]) * a[2];
  const long double a0 = static_cast<long double>(a[0]) * a[1] * a[2];
  const long double b1 = static_cast<long double>(b[0]) * b[1] + b[0] + b[1];
  const long double b0 = static_cast<long double>(b[0]) * b[1];
  const long double b2 = static_cast<long double>(b[0]) + b[1] + 1.0L;
  const long double c = a2 + (a1 - b1) / (gap + 1.0L);
  const long double c3 = a2 + gap + 1.0L - b2;
  const long double c2 = a1 + (c + gap + 1.0L) * a2 - b1 - c * b2;
  const long double c1 = a0 + (c + gap + 1.0L) * a1 - b0 - c * b1;
  const long double c0 = (c + gap + 1.0L) * a0 - c * b0;
  const long double x = n;
  const long double denominator = gap * std::min(1.0L, 1.0L + b[0] / x) *
                                   std::min(1.0L, 1.0L + b[1] / x);
  const long double residual =
      (std::abs(c3) + (std::abs(c2) + (std::abs(c1) + std::abs(c0) / x) / x) / x) / denominator;
  return {term * (x + c) / gap, residual * (std::abs(term) + tail)};
}

// Evaluate regularized 3F2 at unity through denominator Gamma poles
// 3F2_reg(1) = sum_k (a1)_k(a2)_k(a3)_k / [Gamma(b1+k)Gamma(b2+k)k!]
// [REFERENCE: NIST DLMF 16.2.5, https://dlmf.nist.gov/16.2.E5]
double RegularizedHyper3F2Unit(double a1, double a2, double a3, double b1, double b2);

// Ordinary Legendre polynomials P_l(x) to cross check algorithmic
// (l+1)P_{l+1}(x) = (2l+1)xP_l(x) - lP_{l-1}(x)
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
    return (429 * x * x * x * x * x * x * x - 693 * x * x * x * x * x + 315 * x * x * x - 35 * x) / 16;
  } else if (l == 8) {
    return (6435 * x * x * x * x * x * x * x * x - 12012 * x * x * x * x * x * x + 6930 * x * x * x * x - 1260 * x * x +
            35) /
           128;
  } else {
    throw std::invalid_argument("LegendrePl: Not valid l = " + std::to_string(l));
  }
}

// Associated Legendre Polynomial P_l^m(x) solved numerically
// (l-m)P_l^m = (2l-1)xP_{l-1}^m - (l+m-1)P_{l-2}^m
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
//  = sqrt[(2l+1)/(4pi) (l-m)!/(l+m)!] P_l^m(cos(theta)) exp(im phi)
// where P_l^m includes the Condon-Shortley phase (-1)^m
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
// Y_lm^R = sqrt(2)(-1)^m Re Y_l^m for m>0 and -sqrt(2) Im Y_l^m for m<0
//
// Normalization such that \int dOmega |Y_l^m(Omega)|^2 = 1
inline double NReY(const std::complex<double> &Y, int m) {
  const double phase = std::abs(m) % 2 == 0 ? 1.0 : -1.0;
  if (m < 0) {
    return -sqrt(2.0) * std::imag(Y);
  } else if (m == 0) {
    return std::real(Y);
  } else {
    return sqrt(2.0) * phase * std::real(Y);
  }
}

// Compute the normalized associated Legendre function including the Condon-Shortley phase
// S_lm = sqrt[(2l+1)(l-m)!/(4pi(l+m)!)] P_l^m, with 0 <= m <= l
// [REFERENCE: NIST DLMF 14.10.3 and 14.30.1, https://dlmf.nist.gov/14.10.E3 and https://dlmf.nist.gov/14.30.E1]
inline double SphericalLegendre(const int l, const int m, const double x) {
  double current = 1.0 / std::sqrt(4.0 * PI);
  if (m > 0) {
    const double sine = std::sqrt((1.0 - x) * (1.0 + x));
    for (int order = 1; order <= m; ++order) {
      current *= -sine * std::sqrt((2.0 * order + 1.0) / (2.0 * order));
    }
  }
  if (l == m) { return current; }
  double lower = current;
  current *= x * std::sqrt(2.0 * m + 3.0);
  for (int degree = m + 2; degree <= l; ++degree) {
    const double j     = degree;
    const double denom = (j - m) * (j + m);
    const double a     = std::sqrt((4.0 * j * j - 1.0) / denom);
    const double b     = std::sqrt((2.0 * j + 1.0) * (j - 1.0 - m) * (j - 1.0 + m) / ((2.0 * j - 3.0) * denom));
    const double next  = a * x * current - b * lower;
    lower              = current;
    current            = next;
  }
  return current;
}

// Real Basis Laplace Spherical Harmonics
// Y_l0^R = S_l0 and Y_lm^R = sqrt(2)(-1)^|m| S_l|m| times cos or sin
inline double Y_real_basis(double costheta, double phi, int l, int m) {
  const double radial = SphericalLegendre(l, std::abs(m), costheta);
  const double phase  = std::abs(m) % 2 == 0 ? 1.0 : -1.0;
  if (m < 0) {
    return phase * std::sqrt(2.0) * radial * std::sin(-m * phi);
  } else if (m == 0) {
    return radial;
  } else {  // m > 0
    return phase * std::sqrt(2.0) * radial * std::cos(m * phi);
  }
}

// Complex Basis Laplace Spherical Harmonics
// Y_l^m = N_lm P_l^m(cos theta) exp(i m phi) for m >= 0, with negative m from conjugation
inline std::complex<double> Y_complex_basis(double costheta, double phi, int l, int m) {
  const double radial = SphericalLegendre(l, std::abs(m), costheta);
  const double phase  = m < 0 && std::abs(m) % 2 != 0 ? -1.0 : 1.0;
  return phase * radial * std::exp(zi * (m * phi));
}

// Manually written Spherical Harmonics to cross-check algorithmic
// implementations
inline std::complex<double> Y_complex_basis_ref(double costheta, double phi, int l, int m) {
  const double theta = std::acos(costheta);

  if (l == 0 && m == 0) return 0.5 * std::sqrt(1.0 / PI);

  // ------------------------------------------------------------------
  if (l == 1 && m == -1) return 0.5 * std::sqrt(3.0 / (2.0 * PI)) * std::exp(-gra::math::zi * phi) * std::sin(theta);

  if (l == 1 && m == 0) return 0.5 * std::sqrt(3.0 / PI) * costheta;

  if (l == 1 && m == 1) return -0.5 * std::sqrt(3.0 / (2.0 * PI)) * std::exp(gra::math::zi * phi) * std::sin(theta);

  // ------------------------------------------------------------------
  if (l == 2 && m == -2)
    return 0.25 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 2);

  if (l == 2 && m == -1)
    return 0.5 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(-gra::math::zi * phi) * std::sin(theta) * costheta;

  if (l == 2 && m == 0) return 0.25 * std::sqrt(5.0 / PI) * (3.0 * std::pow(costheta, 2) - 1.0);

  if (l == 2 && m == 1)
    return -0.5 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(gra::math::zi * phi) * std::sin(theta) * costheta;

  if (l == 2 && m == 2)
    return 0.25 * std::sqrt(15.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 2);

  // ------------------------------------------------------------------

  if (l == 3 && m == -3)
    return 1.0 / 8.0 * std::sqrt(35.0 / PI) * std::exp(-3.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 3);

  if (l == 3 && m == -2)
    return 0.25 * std::sqrt(105.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 2) *
           costheta;

  if (l == 3 && m == -1)
    return 1.0 / 8.0 * std::sqrt(21.0 / PI) * std::exp(-gra::math::zi * phi) * std::sin(theta) *
           (5.0 * std::pow(costheta, 2) - 1.0);

  if (l == 3 && m == 0) return 0.25 * std::sqrt(7.0 / PI) * (5.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 3 && m == 1)
    return -1.0 / 8.0 * std::sqrt(21.0 / PI) * std::exp(gra::math::zi * phi) * std::sin(theta) *
           (5.0 * std::pow(costheta, 2) - 1.0);

  if (l == 3 && m == 2)
    return 0.25 * std::sqrt(105.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 2) *
           costheta;

  if (l == 3 && m == 3)
    return -1.0 / 8.0 * std::sqrt(35.0 / PI) * std::exp(3.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 3);

  // ------------------------------------------------------------------

  if (l == 4 && m == -4)
    return 3.0 / 16.0 * std::sqrt(35.0 / (2.0 * PI)) * std::exp(-4.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 4);

  if (l == 4 && m == -3)
    return 3.0 / 8.0 * std::sqrt(35.0 / (PI)) * std::exp(-3.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 3) *
           costheta;

  if (l == 4 && m == -2)
    return 3.0 / 8.0 * std::sqrt(5.0 / (2.0 * PI)) * std::exp(-2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * (7.0 * std::pow(costheta, 2) - 1.0);

  if (l == 4 && m == -1)
    return 3.0 / 8.0 * std::sqrt(5.0 / PI) * std::exp(-gra::math::zi * phi) * std::sin(theta) *
           (7.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 4 && m == 0)
    return 3.0 / 16.0 * std::sqrt(1.0 / PI) * (35.0 * std::pow(costheta, 4) - 30.0 * std::pow(costheta, 2) + 3.0);

  if (l == 4 && m == 1)
    return -3.0 / 8.0 * std::sqrt(5.0 / PI) * std::exp(gra::math::zi * phi) * std::sin(theta) *
           (7.0 * std::pow(costheta, 3) - 3.0 * costheta);

  if (l == 4 && m == 2)
    return 3.0 / 8.0 * std::sqrt(5.0 / (2.0 * PI)) * std::exp(2.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 2) * (7.0 * std::pow(costheta, 2) - 1.0);

  if (l == 4 && m == 3)
    return -3.0 / 8.0 * std::sqrt(35.0 / (PI)) * std::exp(3.0 * gra::math::zi * phi) * std::pow(std::sin(theta), 3) *
           costheta;

  if (l == 4 && m == 4)
    return 3.0 / 16.0 * std::sqrt(35.0 / (2.0 * PI)) * std::exp(4.0 * gra::math::zi * phi) *
           std::pow(std::sin(theta), 4);

  throw std::invalid_argument("Y_complex_basis_ref:: Not supported l = " + std::to_string(l) +
                              ", m = " + std::to_string(m));
}

}  // namespace math
}  // namespace gra

#endif
