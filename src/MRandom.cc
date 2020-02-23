// Random number class
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MMath.h"
#include "Graniitti/MRandom.h"

namespace gra {
// Poisson distributed random numbers with mean lambda
int MRandom::PoissonRandom(double lambda) {
  const unsigned int ALGO = 2;

  // Knuth algorithm
  if (ALGO == 1) {
    const double L = std::exp(-lambda);
    int          k = 0;
    double       p = 1;

    do {
      ++k;
      p *= U(0, 1);
    } while (p > L);
    return k - 1;
  }
  // Inversion search algorithm
  if (ALGO == 2) {
    int    k = 0;
    double p = std::exp(-lambda);
    double s = p;
    double u = U(0, 1);

    while (u > s) {
      ++k;
      p = p * lambda / k;
      s += p;
    }
    return k;
  }
}

// Exponential random numbers with rate lambda (mean is given by 1/lambda)
double MRandom::ExpRandom(double lambda) {
  const double R = U(0, 1);
  return -1.0 / lambda * std::log(1.0 - R);
}

// Powerlaw random numbers from [a,b] with exponent alpha (e.g. -2 for ~ 1/x^2)
double MRandom::PowerRandom(double a, double b, double alpha) {
  const double min_pow = std::pow(a, alpha + 1.0);
  const double max_pow = std::pow(b, alpha + 1.0);
  return std::pow((max_pow - min_pow) * U(0, 1) + min_pow, 1.0 / (alpha + 1));
}

// Uniform random numbers from [a,b)
double MRandom::U(double a, double b) {
  const double value = flat(rng);
  return a + (b - a) * value;
}

// Gaussian random numbers from normal distribution with mean mu and std sigma
double MRandom::G(double mu, double sigma) {
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
// Return mass (GeV)
//
double MRandom::RelativisticBWRandom(double m0, double Gamma, double LIMIT, double M_MIN) {
  // No width case
  if (Gamma < 1e-40) { return m0; }

  const double m2    = math::pow2(m0);
  const double m2max = gra::math::pow2(m0 + LIMIT * Gamma);
  const double m2min =
      std::max(std::max(0.0, math::pow2(M_MIN)), gra::math::pow2(m0 - LIMIT * Gamma));

  double m2val = 0.0;
  while (true) {
    const double R = U(0, 1);
    m2val =
        m2 +
        m0 * Gamma *
            std::tan(std::atan2(m2min - m2, m0 * Gamma) +
                     R * (std::atan2(m2max - m2, m0 * Gamma) - std::atan2(m2min - m2, m0 * Gamma)));

    // Note >= handles massless case, otherwise stuck with m = 0
    if (m2val >= 0) { break; }
  }
  return gra::math::msqrt(m2val);
}

// Cauchy (non-relativistic Breit-Wigner) sampling
//
// Input as with RelativisticBWRandom
//
double MRandom::CauchyRandom(double m0, double Gamma, double LIMIT, double M_MIN) {
  if (Gamma < 1e-40) { return m0; }

  const double mmax = m0 + LIMIT * Gamma;
  const double mmin = std::max(std::max(0.0, M_MIN), m0 - LIMIT * Gamma);

  double mval = 0.0;
  while (true) {
    const double R     = U(0, 1);
    const double value = (Gamma * 0.5) * std::tan(gra::math::PI * (R - 0.5)) + m0;

    if (value < mmax && value > mmin) {
      mval = value;
      break;
    }
  }
  return mval;
}

// K-dimensional Dirichlet distribution with parameter vector alpha of length K
void MRandom::DirRandom(const std::vector<double> &alpha, std::vector<double> &y) {
  unsigned int K = alpha.size();

  // Draw K independent samples from Gamma(shape = alpha_i, scale = 1)
  y.resize(K, 0.0);
  for (std::size_t i = 0; i < K; ++i) {
    std::gamma_distribution<double> gammarnd(alpha[i], 1.0);
    y[i] = gammarnd(rng);  // Draw sample
  }
  // Take the sum
  double y_sum = 0.0;
  for (std::size_t i = 0; i < K; ++i) { y_sum += y[i]; }
  // Normalize
  for (std::size_t i = 0; i < K; ++i) { y[i] /= y_sum; }
}

// Negative binomial distribution with parameters avgN and k
double MRandom::NBDpdf(int n, double avgN, double k) {
  // return Cbinom(n+k-1,n)
  return tgamma(n + k) / (tgamma(n + 1) * tgamma(k)) * std::pow(avgN / (k + avgN), n) *
         std::pow(k / (k + avgN), k);
}

// Random sample from NBD distribution with parameters avgN and k
int MRandom::NBDRandom(double avgN, double k, int maxvalue) {
  const unsigned int MAXTRIAL = 1e7;

  // Random integer from [0,NBins-1]
  std::uniform_int_distribution<int> RANDI(1, maxvalue);

  // Acceptance-Rejection
  unsigned int trials = 0;
  while (true) {
    const int    n   = RANDI(rng);
    const double val = NBDpdf(n, avgN, k);
    if (U(0, 1) < val) { return n; }
    ++trials;
    if (trials > MAXTRIAL) {
      return avgN;  // Return the mean value
    }
  }
}

// Log-distribution with with parameter p
double MRandom::Logpdf(int k, double p) {
  return -1.0 / std::log(1 - p) * std::pow(p, k) / static_cast<double>(k);
}

// Random sample from log-distribution with parameter p
int MRandom::LogRandom(double p, int maxvalue) {
  // Random integer from [0,NBins-1]
  std::uniform_int_distribution<int> RANDI(0, maxvalue - 1);

  // Acceptance-Rejection
  while (true) {
    const int    n   = RANDI(rng);
    const double val = Logpdf(n, p);
    if (U(0, 1) < val) { return n; }
  }
}

}  // namespace gra
