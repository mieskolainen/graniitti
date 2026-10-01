// Random number class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra {
namespace {

// Compute a truncated Breit-Wigner quantile without subtracting angles near pi/2
long double BWQuantile(long double lower, long double upper, long double pole,
                       long double scale, double unit) {
  if (!(unit > 0.0)) { return lower; }
  if (!(unit < 1.0)) { return upper; }
  const long double left = lower - pole;
  const long double right = upper - pole;
  long double offset = 0.0L;
  if (left > 0.0L) {
    const long double angle = (1.0L - unit) * std::atan2(scale, left) +
                              unit * std::atan2(scale, right);
    offset = scale / std::tan(angle);
  } else if (right < 0.0L) {
    const long double angle = (1.0L - unit) * std::atan2(scale, -left) +
                              unit * std::atan2(scale, -right);
    offset = -scale / std::tan(angle);
  } else {
    const long double angle = (1.0L - unit) * std::atan2(left, scale) +
                              unit * std::atan2(right, scale);
    offset = scale * std::tan(angle);
  }
  return std::clamp(pole + offset, lower, upper);
}

// Compute log[P(n)/P(0)] without cancellation in the NBD Poisson limit
long double NBDLogWeight(int n, long double mean, long double shape) {
  long double weight = 0.0L;
  if (n <= 32) {
    for (int j = 0; j < n; ++j) {
      weight += std::log((mean / (j + 1.0L)) * ((shape + j) / (shape + mean)));
    }
    return weight;
  }
  if (shape < 16.0L) {
    // Use reentrant Gamma evaluation because each sampling thread computes multiplicity probabilities
    return gra::math::LogGamma(n + shape) - gra::math::LogGamma(n + 1.0L) - gra::math::LogGamma(shape) -
           n * std::log1p(shape / mean);
  }
  // Evaluate the Stirling remainder in the NBD gamma ratio through x^-9
  // [REFERENCE: NIST DLMF 5.11.1, https://dlmf.nist.gov/5.11.E1]
  const auto correction = [](long double x) {
    const long double inverse = 1.0L / x;
    const long double square = inverse * inverse;
    return inverse * (1.0L / 12.0L + square * (-1.0L / 360.0L + square *
           (1.0L / 1260.0L + square * (-1.0L / 1680.0L + square / 1188.0L))));
  };
  // Remove n log(k) analytically before taking the large shape limit
  const long double remainder =
      (shape + n - 0.5L) * std::log1p(n / shape) - n + correction(shape + n) - correction(shape);
  return n * (std::log(mean) - std::log1p(mean / shape)) - gra::math::LogGamma(n + 1.0L) + remainder;
}

}  // namespace

// Sample Poisson counts with intensity lambda, optionally conditioned on n >= 1
// P(n|lambda) = exp(-lambda) lambda^n/n!, normalized by 1-exp(-lambda) when conditioned
int MRandom::PoissonRandom(double lambda, bool positive) {
  if (!std::isfinite(lambda) || lambda < 0.0 ||
      lambda > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("MRandom::PoissonRandom: mean must be finite, non-negative and representable");
  }
  if (std::fpclassify(lambda) == FP_ZERO) {
    if (positive) { throw std::invalid_argument("MRandom::PoissonRandom: positive count requires positive mean"); }
    return 0;
  }
  // Draw the first arrival conditional on T < 1, then the remaining Poisson process
  if (positive) { lambda = std::max(0.0, lambda + std::log1p(U(0.0, 1.0) * std::expm1(-lambda))); }
  if (!(lambda > 0.0)) { return static_cast<int>(positive); }
  std::poisson_distribution<int> distribution(lambda);
  return static_cast<int>(positive) + distribution(rng);
}

// Exponential random numbers with rate lambda (mean is given by 1/lambda)
// x = -log(1-u)/lambda
double MRandom::ExpRandom(double lambda) {
  if (!std::isfinite(lambda) || lambda <= 0.0) {
    throw std::invalid_argument("MRandom::ExpRandom: rate must be finite and positive");
  }
  const double u = U(0, 1);
  return -std::log1p(-u) / lambda;
}

// Exponential pdf
// f(x) = lambda exp(-lambda x), x >= 0
double MRandom::ExpPdf(double x, double lambda) {
    if (!std::isfinite(lambda) || lambda <= 0.0) {
      throw std::invalid_argument("MRandom::ExpPdf: rate must be finite and positive");
    }
    if (!std::isfinite(x)) {
      throw std::invalid_argument("MRandom::ExpPdf: value must be finite");
    }
    if (x < 0) return 0;
    return lambda * std::exp(-lambda * x);
}

// Bounded [a,b] exponential random numbers
// x = a - log[1-u(1-exp(-lambda(b-a)))]/lambda
double MRandom::ExpBoundedRandom(double a, double b, double lambda) {
    if (!std::isfinite(a) || !std::isfinite(b) || !(a < b) ||
        !std::isfinite(lambda) || lambda <= 0.0) {
      throw std::invalid_argument("MRandom::ExpBoundedRandom: require finite a < b and positive rate");
    }
    const long double u = U(0, 1);
    const long double span = static_cast<long double>(b) - a;
    const long double rate = lambda;
    const long double norm = -std::expm1(-rate * span);
    const long double value = static_cast<long double>(a) - std::log1p(-u * norm) / rate;
    return std::clamp(static_cast<double>(value), a, b);
}

// Bounded [a,b] exponential pdf
// f(x) = lambda exp[-lambda(x-a)]/[1-exp(-lambda(b-a))]
double MRandom::ExpBoundedPdf(double x, double a, double b, double lambda) {
    if (!std::isfinite(a) || !std::isfinite(b) || !(a < b) ||
        !std::isfinite(lambda) || lambda <= 0.0) {
      throw std::invalid_argument("MRandom::ExpBoundedPdf: require finite a < b and positive rate");
    }
    if (!std::isfinite(x)) {
      throw std::invalid_argument("MRandom::ExpBoundedPdf: value must be finite");
    }
    if (x < a || x > b) return 0;
    const long double rate = lambda;
    const long double span = static_cast<long double>(b) - a;
    const long double offset = static_cast<long double>(x) - a;
    return static_cast<double>(rate * std::exp(-rate * offset) / -std::expm1(-rate * span));
}

// Powerlaw random numbers from [a,b] with exponent alpha (e.g. 2 for ~ 1/x^2)
// x = [a^(1-alpha)+u(b^(1-alpha)-a^(1-alpha))]^(1/(1-alpha)), with the log limit at alpha=1
double MRandom::PowerBoundedRandom(double a, double b, double alpha) {
  if (!std::isfinite(a) || !std::isfinite(b) || !std::isfinite(alpha) ||
      a <= 0.0 || b <= a) {
    throw std::invalid_argument(
        "MRandom::PowerBoundedRandom: require finite 0 < a < b and exponent");
  }
  const long double u = U(0, 1);
  const long double log_a = std::log(static_cast<long double>(a));
  const long double log_b = std::log(static_cast<long double>(b));
  const long double span = std::log1p((static_cast<long double>(b) - a) / a);
  const long double power = 1.0L - static_cast<long double>(alpha);
  long double log_x = log_a + u * span;
  if (power < 0.0L) {
    log_x = log_a + std::log1p(u * std::expm1(power * span)) / power;
  } else if (power > 0.0L) {
    // Reflect the uniform draw so the exponential argument remains nonpositive
    log_x = log_b + std::log1p(u * std::expm1(-power * span)) / power;
  }
  return std::clamp(static_cast<double>(std::exp(log_x)), a, b);
}

// Bounded [a,b] powerlaw pdf
// f(x) = (1-alpha)x^(-alpha)/[b^(1-alpha)-a^(1-alpha)], with f=1/[x log(b/a)] at alpha=1
double MRandom::PowerBoundedPdf(double x, double a, double b, double alpha) {
    if (!std::isfinite(x) || !std::isfinite(a) || !std::isfinite(b) ||
        !std::isfinite(alpha) || a <= 0.0 || b <= a) {
      throw std::invalid_argument(
          "MRandom::PowerBoundedPdf: require finite value, 0 < a < b and exponent");
    }
    if (x < a || x > b) return 0;
    const long double log_x = std::log(static_cast<long double>(x));
    const long double span = std::log1p((static_cast<long double>(b) - a) / a);
    const long double power = 1.0L - static_cast<long double>(alpha);
    long double log_pdf = -log_x - std::log(span);
    if (power < 0.0L) {
      log_pdf = std::log(-power) - std::log(-std::expm1(power * span)) +
                power * std::log1p((static_cast<long double>(x) - a) / a) - log_x;
    } else if (power > 0.0L) {
      log_pdf = std::log(power) - std::log(-std::expm1(-power * span)) -
                power * std::log1p((static_cast<long double>(b) - x) / x) - log_x;
    }
    return static_cast<double>(std::exp(log_pdf));
}

// Uniform random numbers from [a,b)
// x = a + (b-a)u
double MRandom::U(double a, double b) {
  if (!std::isfinite(a) || !std::isfinite(b) || a > b) {
    throw std::invalid_argument(
        "MRandom::U: bounds must be finite and ordered");
  }
  if (std::is_eq(a <=> b)) {
    return a;
  }
  const double value = flat(rng);
  return std::min(std::lerp(a, b, value), std::nextafter(b, a));
}

// Gaussian random numbers from normal distribution with mean mu and std sigma
// x = mu + sigma z, z ~ N(0,1)
double MRandom::G(double mu, double sigma) {
  if (!std::isfinite(mu) || !std::isfinite(sigma) || sigma < 0.0) {
    throw std::invalid_argument(
        "MRandom::G: mean and width must be finite with nonnegative width");
  }
  if (std::fpclassify(sigma) == FP_ZERO) {
    return mu;
  }
  const double value = gaussian(rng);
  return mu + sigma * value;
}

// Relativistic Breit-Wigner sampling f(s) ~ 1/[(s-m^2)^2 + m^2Gamma^2]
//
// Input: m0     = Pole mass (GeV)
//        Gamma  = Full width (GeV)
//        LIMIT  = gives m0 +- LIMIT * GAMMA,
//        M_MIN  = optional minimum bound
//
// Compute mass (GeV)
//
double MRandom::RelativisticBWRandom(double m0, double Gamma, double LIMIT, double M_MIN) {
  if (!std::isfinite(m0) || !std::isfinite(Gamma) ||
      !std::isfinite(LIMIT) || !std::isfinite(M_MIN) || m0 < 0.0 ||
      Gamma < 0.0 || LIMIT <= 0.0 || M_MIN < 0.0) {
    throw std::invalid_argument(
        "MRandom::RelativisticBWRandom: invalid mass, width or bounds");
  }
  // No width case
  if (Gamma < 1e-40) {
    if (m0 < M_MIN) {
      throw std::invalid_argument(
          "MRandom::RelativisticBWRandom: pole mass is below minimum bound");
    }
    return m0;
  }

  const double m2    = math::pow2(m0);
  const double mmax   = m0 + LIMIT * Gamma;
  const double mmin   = std::max(M_MIN, std::max(0.0, m0 - LIMIT * Gamma));
  const double m2max  = gra::math::pow2(mmax);
  const double m2min  = gra::math::pow2(mmin);
  if (!(mmax > mmin) || std::fpclassify(m0) == FP_ZERO) {
    throw std::invalid_argument(
        "MRandom::RelativisticBWRandom: empty mass range or zero pole mass");
  }

  const double scale = m0 * Gamma;
  if (!std::isfinite(m2max) || !std::isfinite(m2min) || !std::isfinite(scale) || !(scale > 0.0)) {
    throw std::invalid_argument("MRandom::RelativisticBWRandom: invalid mass interval or width scale");
  }
  const long double mass2 = BWQuantile(m2min, m2max, m2, scale, U(0, 1));
  return std::clamp(static_cast<double>(std::sqrt(mass2)), mmin, mmax);
}

// Relativistic Breit-Wigner mass-squared integral over a finite range
// I = [atan((s-m0^2)/(m0 Gamma))]_{s_min}^{s_max}/(m0 Gamma)
double MRandom::RelativisticBWMass2Integral(double m0, double Gamma,
                                            double m2min, double m2max) {
  if (!std::isfinite(m0) || !std::isfinite(Gamma) ||
      !std::isfinite(m2min) || !std::isfinite(m2max) || m0 <= 0.0 ||
      Gamma <= 0.0 || m2min < 0.0 || m2max < 0.0) {
    throw std::invalid_argument(
        "MRandom::RelativisticBWMass2Integral: mass and width must be positive "
        "and bounds must be finite and nonnegative");
  }
  if (!(m2max > m2min)) { return 0.0; }

  const double scale = m0 * Gamma;
  if (!std::isfinite(scale) || scale <= 0.0) {
    throw std::overflow_error(
        "MRandom::RelativisticBWMass2Integral: mass-width scale overflow");
  }

  const double m02 = math::pow2(m0);
  if (!std::isfinite(m02)) {
    throw std::overflow_error(
        "MRandom::RelativisticBWMass2Integral: pole mass squared overflow");
  }
  // atan(b/a) - atan(c/a) = atan2(a(b-c), a^2+bc) for b > c and a > 0
  const long double left = static_cast<long double>(m2min) - m02;
  const long double right = static_cast<long double>(m2max) - m02;
  const long double span = static_cast<long double>(m2max) - m2min;
  const long double a = scale;
  const double integral = static_cast<double>(std::atan2(a * span, a * a + left * right) / a);
  if (!std::isfinite(integral) || integral < 0.0) {
    throw std::overflow_error(
        "MRandom::RelativisticBWMass2Integral: integral overflow");
  }
  return integral;
}

// Cauchy (non-relativistic Breit-Wigner) sampling
// f(m) proportional to 1/[(m-m0)^2+(Gamma/2)^2]
//
// Input as with RelativisticBWRandom
//
double MRandom::CauchyRandom(double m0, double Gamma, double LIMIT, double M_MIN) {
  if (!std::isfinite(m0) || !std::isfinite(Gamma) ||
      !std::isfinite(LIMIT) || !std::isfinite(M_MIN) || m0 < 0.0 ||
      Gamma < 0.0 || LIMIT <= 0.0 || M_MIN < 0.0) {
    throw std::invalid_argument(
        "MRandom::CauchyRandom: invalid mass, width or bounds");
  }
  if (Gamma < 1e-40) {
    if (m0 < M_MIN) {
      throw std::invalid_argument(
          "MRandom::CauchyRandom: pole mass is below minimum bound");
    }
    return m0;
  }

  const double mmax = m0 + LIMIT * Gamma;
  const double mmin = std::max(std::max(0.0, M_MIN), m0 - LIMIT * Gamma);

  if (!(mmax > mmin)) {
    throw std::invalid_argument("MRandom::CauchyRandom: empty mass range");
  }
  if (!std::isfinite(mmax)) {
    throw std::invalid_argument("MRandom::CauchyRandom: non-finite mass interval");
  }
  return static_cast<double>(BWQuantile(mmin, mmax, m0, 0.5L * Gamma, U(0, 1)));
}

// K-dimensional Dirichlet distribution with parameter vector alpha of length K
// y_i = g_i/sum_j g_j with g_i ~ Gamma(alpha_i,1)
void MRandom::DirRandom(const std::vector<double> &alpha, std::vector<double> &y) {
  if (alpha.empty() ||
      !std::all_of(alpha.begin(), alpha.end(), [](const double value) {
        return std::isfinite(value) && value > 0.0;
      })) {
    throw std::invalid_argument(
        "MRandom::DirRandom: concentration values must be finite and positive");
  }
  const std::size_t K = alpha.size();

  // Gamma(a) = Gamma(a+1) U^(1/a) for a < 1, evaluated in log space
  std::vector<long double> log_y(K);
  for (const auto &i : indices(alpha)) {
    const long double shape = alpha[i];
    std::gamma_distribution<long double> gamma(shape < 1.0L ? shape + 1.0L : shape, 1.0L);
    long double value = 0.0L;
    do { value = gamma(rng); } while (!(value > 0.0L));
    log_y[i] = std::log(value);
    if (shape < 1.0L) {
      log_y[i] += std::log1p(-static_cast<long double>(U(0, 1))) / shape;
    }
  }
  const long double maximum = *std::max_element(log_y.begin(), log_y.end());
  for (long double &value : log_y) { value = std::exp(value - maximum); }
  const long double sum = gra::Sum(log_y);
  y.resize(K);
  for (const auto &i : indices(y)) { y[i] = static_cast<double>(log_y[i] / sum); }
}

// Random sample from NBD distribution with parameters avgN and k
// P(n) = Gamma(n+k)/[Gamma(k)n!] (avgN/(k+avgN))^n (k/(k+avgN))^k
int MRandom::NBDRandom(double avgN, double k, int maxvalue) {
  if (!std::isfinite(avgN) || !std::isfinite(k) || avgN < 0.0 || k <= 0.0 ||
      maxvalue < 0) {
    throw std::invalid_argument(
        "MRandom::NBDRandom: invalid mean, shape or maximum");
  }
  if (std::fpclassify(avgN) == FP_ZERO || maxvalue == 0) {
    return 0;
  }
  const unsigned int MAXTRIAL = 1e7;
  const long double shape = k;
  const long double mean = avgN;
  const int mode = static_cast<int>(std::min(
      static_cast<long double>(maxvalue),
      shape > 1.0L ? std::floor(mean * (1.0L - 1.0L / shape)) : 0.0L));
  const long double log_p = -std::log1p(shape / mean);
  const long double mode_weight = NBDLogWeight(mode, mean, shape);

  // Random integer from [0,NBins-1]
  std::uniform_int_distribution<int> RANDI(0, maxvalue);

  // Acceptance-Rejection
  unsigned int trials = 0;
  while (true) {
    const int    n   = RANDI(rng);
    // Compare to the mode without including the possibly underflowing P(0)
    const long double log_u = std::log1p(-static_cast<long double>(U(0, 1)));
    long double log_ratio = 0.0L;
    if (std::fpclassify(shape - 1.0L) == FP_ZERO) {
      log_ratio = n * log_p;
    } else if (std::abs(n - mode) > 32) {
      log_ratio = NBDLogWeight(n, mean, shape) - mode_weight;
    } else {
      // Resolve nearby probabilities directly with at most 32 recurrence steps
      int j = mode;
      while (j != n && log_ratio > log_u) {
        // Combine the factors before taking logarithms in the Poisson limit
        if (j < n) {
          log_ratio += std::log((mean / (j + 1.0L)) * ((shape + j) / (shape + mean)));
          ++j;
        } else {
          log_ratio += std::log((j / mean) * ((shape + mean) / (shape + (j - 1.0L))));
          --j;
        }
      }
    }
    if (log_u < std::min(log_ratio, 0.0L)) { return n; }
    ++trials;
    if (trials > MAXTRIAL) {
      throw PhaseSpaceFailure("MRandom::NBDRandom: truncated multiplicity sampling exhausted its trials");
    }
  }
}

// Negative binomial distribution with parameters avgN and k
// P(n) = Gamma(n+k)/[Gamma(k)n!] (avgN/(k+avgN))^n (k/(k+avgN))^k
double MRandom::NBDPdf(int n, double avgN, double k) {
  if (!std::isfinite(avgN) || !std::isfinite(k) || avgN < 0.0 || k <= 0.0) {
    throw std::invalid_argument("MRandom::NBDPdf: invalid mean or shape");
  }
  if (n < 0) {
    return 0.0;
  }
  if (std::fpclassify(avgN) == FP_ZERO) {
    return n == 0 ? 1.0 : 0.0;
  }
  const long double shape = k;
  const long double mean = avgN;
  const long double log_zero = -shape * std::log1p(mean / shape);
  return static_cast<double>(std::exp(log_zero + NBDLogWeight(n, mean, shape)));
}

// Random sample from log-distribution with parameter p
int MRandom::LogRandom(double p, int maxvalue) {
  if (!std::isfinite(p) || p <= 0.0 || p >= 1.0 || maxvalue < 1) {
    throw std::invalid_argument("MRandom::LogRandom: require 0 < p < 1 and maxvalue >= 1");
  }

  // Log-distribution support starts at k = 1
  std::uniform_int_distribution<int> RANDI(1, maxvalue);

  // Acceptance-Rejection
  while (true) {
    const int    n   = RANDI(rng);
    const double val = LogPdf(n, p);
    if (U(0, 1) < val) { return n; }
  }
}

// Log-distribution with with parameter p
// P(k) = -p^k/[k log(1-p)], k >= 1
double MRandom::LogPdf(int k, double p) {
  if (!std::isfinite(p) || p <= 0.0 || p >= 1.0) {
    throw std::invalid_argument("MRandom::LogPdf: require 0 < p < 1");
  }
  if (k < 1) { return 0.0; }
  return std::pow(p, k - 1) * (p / -std::log1p(-p)) / static_cast<double>(k);
}

}  // namespace gra
