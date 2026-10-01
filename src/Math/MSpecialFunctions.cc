// Special function implementations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MSpecialFunctions.h"

#include <cmath>
#include <iomanip>
#include <limits>
#include <numeric>
#include <optional>
#include <sstream>
#include <stdexcept>

#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::math {

using gra::aux::indices;

namespace {

// Prepare R(n) = n sum_j c_j/n^j and the polynomial of r_n[1+R(n+1)]-R(n)
// The truncation error is bounded by the recurrence residual times the absolute series tail
// [REFERENCE: J L Willis, Numer Algorithms 59 (2012) 447, https://arxiv.org/abs/1102.3003]
class Hyper3F2Expansion {
 public:
  // Use a quartic inverse power expansion to keep preparation cheaper than direct summation
  enum : std::size_t { order = 4 };

  // Cancel successive inverse powers of the exact rational term recurrence
  Hyper3F2Expansion(const std::array<double, 3>& a, const std::array<double, 2>& b) {
    // Expand the three linear factors in the hypergeometric term ratio
    const auto cubic = [](const long double x, const long double y, const long double z) {
      return std::array<long double, 4>{1.0L, x + y + z, x * y + x * z + y * z, x * y * z};
    };
    const auto        A     = cubic(a[0], a[1], a[2]);
    const auto        B     = cubic(b[0], b[1], 1.0L);
    const auto        abs_A = cubic(std::abs(a[0]), std::abs(a[1]), std::abs(a[2]));
    const auto        abs_B = cubic(std::abs(b[0]), std::abs(b[1]), 1.0L);
    const long double gap   = static_cast<long double>(b[0]) + b[1] - a[0] - a[1] - a[2];

    // Store exact binomial coefficients for the shifted inverse powers
    constexpr auto binomial = [] {
      std::array<std::array<long double, order + 4>, order + 1> out{};
      for (std::size_t n = 0; n <= order; ++n) {
        for (std::size_t k = 0; k <= n; ++k) { out[n][k] = Cbinom(n, k); }
      }
      return out;
    }();
    std::array<long double, 2 * order + 3> rounding{};
    for (std::size_t i = 0; i <= order + 2; ++i) {
      for (std::size_t r = 0; r <= std::min(i, std::size_t{3}); ++r) {
        residual[i + 1] += A[r] * binomial[order - 1][i - r];
        rounding[i + 1] += abs_A[r] * binomial[order - 1][i - r];
      }
    }
    for (std::size_t j = 0; j <= order; ++j) {
      coefficient[j] = residual[j + 1] / (gap + j);
      for (std::size_t i = 1; i <= order + 2 + (j == 0); ++i) {
        long double basis     = 0.0L;
        long double magnitude = 0.0L;
        for (std::size_t r = 0; r <= std::min(i, std::size_t{3}); ++r) {
          basis += A[r] * binomial[order - j][i - r] - B[r] * binomial[order - 1][i - r];
          magnitude += abs_A[r] * binomial[order - j][i - r] + abs_B[r] * binomial[order - 1][i - r];
        }
        // Use the exact linear coefficient without subtracting the expansion orders
        if (i == 1) {
          basis     = -gap - j;
          magnitude = abs_A[1] + abs_B[1] + 1.0L + j;
        }
        residual[j + i] += coefficient[j] * basis;
        rounding[j + i] += std::abs(coefficient[j]) * magnitude;
      }
    }
    for (const auto i : indices(residual)) {
      residual[i] = std::abs(residual[i]) +
                    static_cast<long double>(order + 8) * std::numeric_limits<long double>::epsilon() * rounding[i];
    }
  }

  // Bound the residual for every k >= n after all shifted parameters become positive
  std::pair<long double, long double> Evaluate(const std::array<double, 3>& a, const std::array<double, 2>& b,
                                               const int n, const long double term) const {
    const long double tail_bound = Hyper3F2Tail(a, b, n, term);
    if (!std::isfinite(tail_bound)) { return {0.0L, tail_bound}; }
    const long double x = 1.0L / n;
    // Evaluate the inverse power polynomials by Horner accumulation
    const auto        horner = [x](const long double sum, const long double value) { return sum * x + value; };
    const long double tail   = term * n * std::accumulate(coefficient.rbegin(), coefficient.rend(), 0.0L, horner);
    const long double bound  = std::accumulate(residual.rbegin(), residual.rend() - 1, 0.0L, horner) /
                               (std::min(1.0L, 1.0L + b[0] * x) * std::min(1.0L, 1.0L + b[1] * x));
    return {tail, bound * (std::abs(term) + tail_bound)};
  }

 private:
  std::array<long double, order + 1>     coefficient{};
  std::array<long double, 2 * order + 3> residual{};
};

// Recompute a cancellation dominated partial sum with the widest compiler arithmetic
// Keep the common Gamma factor outside the sum so its error is relative to the final value
// [REFERENCE: GCC additional floating types, https://gcc.gnu.org/onlinedocs/gcc/Floating-Types.html]
std::pair<long double, long double> Hyper3F2Refine(const std::array<double, 3>& a, const std::array<double, 2>& b,
                                                   const int first, const int n, const long double remainder) {
#if defined(__SIZEOF_FLOAT128__)
  using Real           = __float128;
  constexpr int digits = 113;  // IEEE binary128 significand
#else
  using Real           = long double;
  constexpr int digits = std::numeric_limits<Real>::digits;
#endif
  Real term = 1;
  for (int k = 0; k < first; ++k) { term *= (Real(a[0]) + k) * (Real(a[1]) + k) * (Real(a[2]) + k) / (k + 1); }
  Real        sum      = term;
  long double rounding = (8.0L * first + 16.0L) * std::abs(static_cast<long double>(term));
  for (int k = first; k < n; ++k) {
    term *= (Real(a[0]) + k) * (Real(a[1]) + k) * (Real(a[2]) + k) / ((Real(b[0]) + k) * (Real(b[1]) + k) * (k + 1));
    sum += term;
    rounding += 8.0L * (k + 1) * std::abs(static_cast<long double>(term)) + std::abs(static_cast<long double>(sum));
  }
  const Real        tail   = term * Real(remainder);
  const long double factor = GammaRatio(std::array<double, 0>{}, std::array{b[0] + first, b[1] + first});
  const long double value  = static_cast<long double>(sum + tail) * factor;
  rounding += (8.0L * n + 16.0L) * std::abs(static_cast<long double>(tail));
  return {value, std::ldexp(1.0L, 1 - digits) * rounding * std::abs(factor) +
                     16.0L * std::numeric_limits<long double>::epsilon() * std::abs(value)};
}

}  // namespace

// Evaluate regularized 3F2 at unity with a bounded inverse power tail
// [REFERENCE: J L Willis, Numer Algorithms 59 (2012) 447, https://arxiv.org/abs/1102.3003]
double RegularizedHyper3F2Unit(const double a1, const double a2, const double a3, const double b1, const double b2) {
  constexpr int         max_terms = 10000;
  constexpr long double tolerance = 1.0e-14L;
  if (!std::isfinite(a1) || !std::isfinite(a2) || !std::isfinite(a3) || !std::isfinite(b1) || !std::isfinite(b2)) {
    throw gra::AmplitudeFailure("RegularizedHyper3F2Unit: invalid argument");
  }
  // Identify exact terminating upper parameters before applying convergence conditions
  const auto terminal_degree = [](const double value) {
    if (value > 0.0 || value <= -max_terms || !math::IsExactEqual(value, std::round(value))) { return max_terms; }
    return static_cast<int>(-std::llround(value));
  };
  // Skip exact reciprocal Gamma zeros in the regularized series
  const auto pole_start = [](const double value) {
    if (value > 0.0 || !math::IsExactEqual(value, std::round(value))) { return 0; }
    return value <= -max_terms ? max_terms : 1 - static_cast<int>(std::llround(value));
  };
  const int terminal = std::min({terminal_degree(a1), terminal_degree(a2), terminal_degree(a3)});
  const int first    = std::max(pole_start(b1), pole_start(b2));
  if (terminal < first) { return 0.0; }
  if (first == max_terms) { throw gra::AmplitudeFailure("RegularizedHyper3F2Unit: series starts beyond term limit"); }
  if (terminal == max_terms && !(static_cast<long double>(b1) + b2 - a1 - a2 - a3 > 0.0L)) {
    throw gra::AmplitudeFailure("RegularizedHyper3F2Unit: divergent series at unity");
  }
  long double numerator       = 1.0L;
  long double factorial_value = 1.0L;
  for (int k = 0; k < first; ++k) {
    numerator *=
        (static_cast<long double>(a1) + k) * (static_cast<long double>(a2) + k) * (static_cast<long double>(a3) + k);
    factorial_value *= static_cast<long double>(k + 1);
  }
  long double term =
      numerator * GammaRatio(std::array<double, 0>{}, std::array{b1 + first, b2 + first}) / factorial_value;
  long double                      sum       = term;
  long double                      rounding  = (8.0L * first + 16.0L) * std::abs(term);
  const int                        last      = terminal < max_terms ? terminal : max_terms;
  int                              next_tail = first + 1;
  std::optional<Hyper3F2Expansion> expansion;
  for (int k = first; k < last; ++k) {
    const long double denominator =
        (static_cast<long double>(b1) + k) * (static_cast<long double>(b2) + k) * (k + 1.0L);
    term *= (static_cast<long double>(a1) + k) * (static_cast<long double>(a2) + k) *
            (static_cast<long double>(a3) + k) / denominator;
    sum += term;
    // Weight recurrence roundoff by each term, since the large initial terms may cancel
    rounding += 8.0L * (k + 1) * std::abs(term) + std::abs(sum);
    if (terminal < max_terms || k + 1 != next_tail) { continue; }
    const int n = k + 1;
    // Refine only after the truncation bound is small enough and cancellation limits precision
    const auto converge = [&](const long double tail, const long double error) -> std::optional<double> {
      const long double target = tolerance * std::abs(sum + tail);
      if (!(error <= target)) { return std::nullopt; }
      const long double roundoff =
          std::numeric_limits<long double>::epsilon() * (rounding + (8.0L * n + 16.0L) * std::abs(tail));
      if (error + roundoff <= target) { return FiniteSpecialValue(sum + tail); }
      const auto [refined, refined_error] = Hyper3F2Refine({a1, a2, a3}, {b1, b2}, first, n, tail / term);
      if (error + refined_error <= tolerance * std::abs(refined)) { return FiniteSpecialValue(refined); }
      return std::nullopt;
    };
    const auto [tail, error] = expansion ? expansion->Evaluate({a1, a2, a3}, {b1, b2}, n, term)
                                         : Hyper3F2Remainder({a1, a2, a3}, {b1, b2}, n, term);
    if (const auto value = converge(tail, error)) { return *value; }
    if (!expansion && n >= static_cast<int>(Hyper3F2Expansion::order + 1)) {
      expansion.emplace(std::array{a1, a2, a3}, std::array{b1, b2});
      const auto [accelerated, bound] = expansion->Evaluate({a1, a2, a3}, {b1, b2}, n, term);
      if (const auto value = converge(accelerated, bound)) { return *value; }
    }
    next_tail = std::min(2 * next_tail, last);
  }
  if (terminal < max_terms) { return FiniteSpecialValue(sum); }
  std::ostringstream error;
  error << std::setprecision(std::numeric_limits<double>::max_digits10)
        << "RegularizedHyper3F2Unit: series did not converge for (" << a1 << ", " << a2 << ", " << a3 << "; " << b1
        << ", " << b2 << ")";
  throw gra::AmplitudeFailure(error.str());
}

// Match the ordinary evaluation to the large argument expansion at machine precision
// [REFERENCE: NIST DLMF 10.40.1, large argument expansion of the modified Bessel function]
double LogBesselI1(const double x) {
  // Preserve log(x/2) when I_1(x) itself would underflow
  if (x >= 0.0 && x <= std::sqrt(std::numeric_limits<double>::epsilon())) { return std::log(x) - std::log(2.0); }
  if (x < std::log(std::numeric_limits<double>::max()) / 2.0) { return std::log(std::cyl_bessel_i(1.0, x)); }
  double sum = 1.0, term = 1.0;
  for (unsigned int k = 1; std::abs(term) > std::numeric_limits<double>::epsilon() * std::abs(sum); ++k) {
    term *= (pow2(2.0 * k - 1.0) - 4.0) / (8.0 * k * x);
    sum += term;
  }
  return x - 0.5 * (std::log(2.0 * PI) + std::log(x)) + std::log(sum);
}

// Compute the first-kind Bessel function of integer order zero
// J_0(x) = sum_k (-1)^k (x/2)^(2k)/(k!)^2
// [REFERENCE: Abramowitz and Stegun, Handbook of Mathematical Functions, 1965]
double BesselJ0(const double x) {
  if (!std::isfinite(x)) { return std::numeric_limits<double>::quiet_NaN(); }
  if (math::IsZero(x)) { return 1.0; }

  const double ax = std::abs(x);
  if (ax < 8.0) {
    const double y = x * x;
    const double numerator =
        57568490574.0 +
        y * (-13362590354.0 + y * (651619640.7 + y * (-11214424.18 + y * (77392.33017 - y * 184.9052456))));
    const double denominator =
        57568490411.0 + y * (1029532985.0 + y * (9494680.718 + y * (59272.64853 + y * (267.8532712 + y))));
    return numerator / denominator;
  }

  const double z     = 8.0 / ax;
  const double y     = z * z;
  const double phase = ax - 0.785398164;
  const double p = 1.0 + y * (-0.1098628627e-2 + y * (0.2734510407e-4 + y * (-0.2073370639e-5 + y * 0.2093887211e-6)));
  const double q =
      -0.1562499995e-1 + y * (0.1430488765e-3 + y * (-0.6911147651e-5 + y * (0.7621095161e-6 - y * 0.9349451520e-7)));
  return std::sqrt(0.636619772 / ax) * (p * std::cos(phase) - z * q * std::sin(phase));
}

// Compute the first-kind Bessel function of integer order one
// J_1(x) = sum_k (-1)^k (x/2)^(2k+1)/[k!(k+1)!]
// [REFERENCE: Abramowitz and Stegun, Handbook of Mathematical Functions, 1965]
double BesselJ1(const double x) {
  if (!std::isfinite(x)) { return std::numeric_limits<double>::quiet_NaN(); }
  const double ax = std::abs(x);
  if (ax < 8.0) {
    const double y = x * x;
    const double numerator =
        72362614232.0 +
        y * (-7895059235.0 + y * (242396853.1 + y * (-2972611.439 + y * (15704.48260 - y * 30.16036606))));
    const double denominator =
        144725228442.0 + y * (2300535178.0 + y * (18583304.74 + y * (99447.43394 + y * (376.9991397 + y))));
    return x * numerator / denominator;
  }

  const double z     = 8.0 / ax;
  const double y     = z * z;
  const double phase = ax - 3.0 * PI / 4.0;
  const double p     = 1.0 + y * (0.183105e-2 + y * (-0.3516396496e-4 + y * (0.2457520174e-5 - y * 0.240337019e-6)));
  const double q =
      0.04687499995 + y * (-0.2002690873e-3 + y * (0.8449199096e-5 + y * (-0.88228987e-6 + y * 0.105787412e-6)));
  const double value = std::sqrt(0.636619772 / ax) * (std::cos(phase) * p - z * std::sin(phase) * q);
  return x < 0.0 ? -value : value;
}

// Compute the first three first-kind integer-order Bessel functions
// J_2(x) = 2J_1(x)/x - J_0(x)
std::array<double, 3> BesselJ012(const double x) {
  const double j0 = BesselJ0(x);
  const double j1 = BesselJ1(x);
  double       j2 = 0.0;
  if (std::abs(x) <= 1.0) {
    const double x2 = x * x;
    j2              = x2 * (1.0 / 8.0 +
               x2 * (-1.0 / 96.0 +
                     x2 * (1.0 / 3072.0 + x2 * (-1.0 / 184320.0 + x2 * (1.0 / 17694720.0 - x2 / 2477260800.0)))));
  } else {
    j2 = 2.0 * j1 / x - j0;
  }
  return {j0, j1, j2};
}

// Compute the first-kind Bessel function of one nonnegative integer order
// J_{n+1}(x) = 2n J_n(x)/x - J_{n-1}(x)
// [REFERENCE: C. W. Clenshaw, Mathematical Tables, Vol. 5, 1962]
double BesselJ(const int n, const double x) {
  if (n < 0) { throw std::invalid_argument("math::BesselJ: order must be nonnegative"); }
  if (!std::isfinite(x)) { return std::numeric_limits<double>::quiet_NaN(); }
  if (n == 0) { return BesselJ0(x); }
  if (n == 1) { return BesselJ1(x); }
  if (math::IsZero(x)) { return 0.0; }

  // Use the convergent small-argument series without inverse powers of x
  // [REFERENCE: NIST DLMF 10.2.2, https://dlmf.nist.gov/10.2.E2]
  if (std::abs(x) <= 1.0) {
    long double term = 1.0L;
    for (int order = 1; order <= n; ++order) {
      term *= static_cast<long double>(x) / (2.0L * order);
      if (math::IsZero(term)) { return static_cast<double>(term); }
    }
    const long double ratio = -static_cast<long double>(x) * x / 4.0L;
    long double       sum   = term;
    for (unsigned int k = 1; std::abs(term) > std::numeric_limits<double>::epsilon() * std::abs(sum); ++k) {
      term *= ratio / (static_cast<long double>(k) * (static_cast<long double>(n) + k));
      sum += term;
    }
    return static_cast<double>(sum);
  }

  const double inverse_x = 2.0 / x;
  if (std::abs(x) > static_cast<double>(n)) {
    double lower = BesselJ0(x);
    double value = BesselJ1(x);
    for (int order = 1; order < n; ++order) {
      const double upper = order * inverse_x * value - lower;
      lower              = value;
      value              = upper;
    }
    return value;
  }

  constexpr double accuracy = 40.0;
  constexpr double large    = 1.0e10;
  constexpr double small    = 1.0e-10;
  const int        start    = 2 * static_cast<int>((n + std::floor(std::sqrt(accuracy * n))) / 2.0);
  double           upper    = 0.0;
  double           value    = 1.0;
  double           target   = 0.0;
  double           sum      = 0.0;
  bool             add      = false;
  for (int order = start; order > 0; --order) {
    double lower = order * inverse_x * value - upper;
    upper        = value;
    value        = lower;
    if (std::abs(value) > large) {
      value *= small;
      upper *= small;
      target *= small;
      sum *= small;
    }
    if (add) { sum += value; }
    add = !add;
    if (order == n) { target = upper; }
  }
  sum = 2.0 * sum - value;
  return target / sum;
}

// Bound the absolute exponential Taylor tail after a fixed term count
// R_N <= (x^N/N!)/[1 - x/(N+1)] for x/(N+1) < 1
double ExponentialTailBound(const double magnitude, const std::size_t terms) {
  if (!std::isfinite(magnitude) || magnitude < 0.0 || terms == 0) {
    throw std::invalid_argument("ExponentialTailBound: invalid input");
  }
  if (std::fpclassify(magnitude) == FP_ZERO) { return 0.0; }
  double next = 1.0;
  for (std::size_t n = 1; n <= terms; ++n) {
    next *= magnitude / static_cast<double>(n);
    if (!std::isfinite(next)) { return std::numeric_limits<double>::infinity(); }
  }
  const double ratio = magnitude / static_cast<double>(terms + 1);
  if (!(ratio < 1.0)) { return std::numeric_limits<double>::infinity(); }
  return next / (1.0 - ratio);
}

}  // namespace gra::math
