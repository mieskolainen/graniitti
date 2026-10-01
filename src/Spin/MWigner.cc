// Wigner rotations and SU(2) angular momentum algebra
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Spin/MSpin.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
namespace wigner {

using gra::aux::indices;

namespace {

constexpr double kTolerance          = 1e-9;
constexpr int    kMaxFiniteFactorial = 170;

// Compute whether two spin labels agree within the numerical tolerance
bool Equal(double x, double y) { return std::abs(x - y) < kTolerance; }

// Compute whether x is an integer or half-integer spin label
bool IsHalfInt(double x) { return IsInt(2.0 * x); }

// Compute doubled spin after validating one SU(2) representation label
int SpinX2(double J, const std::string &context) {
  if (!std::isfinite(J) || J < -kTolerance || !IsHalfInt(J) || J > 0.5 * static_cast<double>(kMaxFiniteFactorial)) {
    throw std::invalid_argument(context + ": J is not a nonnegative integer or half-integer");
  }
  return static_cast<int>(std::llround(2.0 * J));
}

// Compute the sign from one integer phase exponent
// phase = (-1)^n
double Phase(double exponent, const std::string &context) {
  if (!IsInt(exponent) || std::abs(exponent) >= static_cast<double>(std::numeric_limits<long long>::max())) {
    throw std::invalid_argument(context + ": phase exponent is not integer");
  }
  return (std::llround(exponent) % 2 == 0) ? 1.0 : -1.0;
}

// Compute one checked nonnegative integer factorial
double Factorial(double value, const std::string &context) {
  if (value < -kTolerance || !IsInt(value) || value > static_cast<double>(kMaxFiniteFactorial)) {
    throw std::invalid_argument(context + ": invalid factorial argument");
  }
  return gra::math::factorial(static_cast<int>(std::llround(value)));
}

// Prepare one small-d polynomial in integer half-angle powers
// d^J_(mp,m)(theta) = F sum_n c_n cos(theta/2)^(2J-n) sin(theta/2)^n
Polynomial PreparePolynomial(double m, double mp, double J, int J2) {
  Polynomial coefficients;
  coefficients.coefficient.assign(J2 + 1, 0.0);
  if (std::abs(m) > J + kTolerance || std::abs(mp) > J + kTolerance) { return coefficients; }
  coefficients.factor = gra::math::msqrt(Factorial(J + m, "gra::wigner::d") * Factorial(J - m, "gra::wigner::d") *
                                         Factorial(J + mp, "gra::wigner::d") * Factorial(J - mp, "gra::wigner::d"));
  for (int k = 0; k <= J2; ++k) {
    if (J - mp - k < -kTolerance || J + m - k < -kTolerance || k + mp - m < -kTolerance) { continue; }
    const double coefficient = Phase(mp - m + k, "gra::wigner::d") /
                               (Factorial(J - mp - k, "gra::wigner::d") * Factorial(J + m - k, "gra::wigner::d") *
                                Factorial(k + mp - m, "gra::wigner::d") * gra::math::factorial(k));
    const int n = static_cast<int>(std::llround(mp - m + 2.0 * k));
    if (n < 0 || n > J2 || !std::isfinite(coefficient)) {
      throw std::invalid_argument("gra::wigner::d: invalid half-angle coefficient");
    }
    coefficients.coefficient[n] = coefficient;
  }
  return coefficients;
}

// Evaluate one prepared polynomial from shared half-angle powers
double SmallDFromPowers(const Polynomial &polynomial, const std::vector<double> &cosine_power,
                        const std::vector<double> &sine_power) {
  const auto &coefficients = polynomial.coefficient;
  double sum = 0.0;
  for (const auto &n : indices(coefficients)) {
    sum += coefficients[n] * cosine_power[coefficients.size() - 1 - n] * sine_power[n];
  }
  return polynomial.factor * sum;
}

// Compute whether a Wigner 3j tuple satisfies the SU(2) selection rules
bool W3jAllowed(double j1, double j2, double j3, double m1, double m2, double m3) {
  if (j1 < 0.0 || j2 < 0.0 || j3 < 0.0 || !IsHalfInt(j1) || !IsHalfInt(j2) || !IsHalfInt(j3) || !IsInt(j1 + j2 + j3) ||
      !Equal(m1 + m2 + m3, 0.0) || !IsInt(j1 - m1) || !IsInt(j2 - m2) || !IsInt(j3 - m3)) {
    return false;
  }
  return j3 <= j1 + j2 + kTolerance && j3 + kTolerance >= std::abs(j1 - j2) && std::abs(m1) <= j1 + kTolerance &&
         std::abs(m2) <= j2 + kTolerance && std::abs(m3) <= j3 + kTolerance;
}

// Compute Clebsch-Gordan selection rules through the common Wigner 3j tuple
bool CGAllowed(double j1, double j2, double m1, double m2, double j, double m) {
  return IsHalfInt(m1) && IsHalfInt(m2) && IsHalfInt(m) && W3jAllowed(j1, j2, j, m1, m2, -m);
}

// Compute whether one Wigner d row belongs to a finite SU(2) representation
bool DAllowed(double J, double m, double mp) {
  return J >= 0.0 && IsHalfInt(J) && IsHalfInt(m) && IsHalfInt(mp) && IsInt(J - m) && IsInt(J - mp) &&
         std::abs(m) <= J + kTolerance && std::abs(mp) <= J + kTolerance;
}

// Compute log(abs(Gamma(x))) and its real sign away from Gamma poles
std::pair<double, int> SignedLogGamma(double x) {
  if (!std::isfinite(x) || (x <= 0.0 && std::abs(x - std::round(x)) < kTolerance)) {
    return {std::numeric_limits<double>::infinity(), 0};
  }
  // Keep the Gamma sign local because amplitudes are evaluated by concurrent sampling threads
  int          sign  = 1;
  const double value = gra::math::LogGamma(x, &sign);
  if (!std::isfinite(value)) { return {std::numeric_limits<double>::infinity(), 0}; }
  return {value, sign};
}

// Compute the Racah triangle coefficient for one Wigner 3j symbol
// Delta = (j1+j2-j3)!(j1-j2+j3)!(-j1+j2+j3)!/(j1+j2+j3+1)!
double Triangle(double j1, double j2, double j3) {
  return Factorial(j1 + j2 - j3, "gra::wigner::Triangle") * Factorial(j1 - j2 + j3, "gra::wigner::Triangle") *
         Factorial(-j1 + j2 + j3, "gra::wigner::Triangle") / Factorial(j1 + j2 + j3 + 1.0, "gra::wigner::Triangle");
}

// Diagonalize J^2 in the fixed-M product basis to avoid Racah factorial cancellation
double CoupledCG(double j1, double j2, double m1, double J, double M) {
  const double first = std::max(-j1, M - j2);
  const double last = std::min(j1, M + j2);
  const std::size_t n = static_cast<std::size_t>(std::llround(last - first)) + 1;
  MMatrix<double> casimir(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    const double a = first + static_cast<double>(i);
    const double b = M - a;
    casimir[i][i] = j1 * (j1 + 1.0) + j2 * (j2 + 1.0) + 2.0 * a * b;
    if (i + 1 < n) {
      const double ladder = std::sqrt((j1 - a) * (j1 + a + 1.0) * (j2 + b) * (j2 - b + 1.0));
      casimir[i][i + 1] = ladder;
      casimir[i + 1][i] = ladder;
    }
  }
  MMatrix<double> states;
  const auto eigenvalues = casimir.SelfAdjointEigenvalues(1e-13, &states);
  const std::size_t col = static_cast<std::size_t>(std::llround(J - std::max(std::abs(M), std::abs(j1 - j2))));
  const std::size_t row = static_cast<std::size_t>(std::llround(m1 - first));
  if (col >= n || row >= n || std::abs(eigenvalues[col] - J * (J + 1.0)) > 1e-9 * (1.0 + J * (J + 1.0))) {
    throw std::runtime_error("gra::wigner::CoupledCG: invalid angular momentum eigensystem");
  }
  // Propagate the positive largest-m1 phase to a well-resolved coefficient
  // Boundary eigenvector entries can be smaller than eigensolver roundoff
  std::size_t peak = 0;
  for (std::size_t i = 1; i < n; ++i) {
    if (std::abs(states[i][col]) > std::abs(states[peak][col])) { peak = i; }
  }
  double current = 1.0;
  double next = 0.0;
  for (std::size_t i = n - 1; i > peak; --i) {
    const double previous = ((J * (J + 1.0) - casimir[i][i]) * current -
        (i + 1 < n ? casimir[i][i + 1] * next : 0.0)) / casimir[i][i - 1];
    next = current;
    current = previous;
    const double scale = std::max(std::abs(current), std::abs(next));
    if (scale > 1e100) { current /= scale; next /= scale; }
  }
  const double phase = std::copysign(1.0, current) * std::copysign(1.0, states[peak][col]);
  return phase * states[row][col];
}

// Compute selected Wigner rows from a normalized symmetric spinor representation
MMatrix<double> SpinorDRows(double theta, const std::vector<double> &m_values,
                            const std::vector<double> &mp_values, double J) {
  const auto rotation = spin::SpinRotation(spin::SpinHalfRotation(theta, 0.0), J);
  MMatrix<double> output(m_values.size(), mp_values.size(), 0.0);
  const spin::SpinRep rep = spin::SpinRep::FromSpin(J, "gra::wigner::SpinorDRows");
  for (std::size_t i = 0; i < m_values.size(); ++i) {
    for (std::size_t j = 0; j < mp_values.size(); ++j) {
      if (std::abs(m_values[i]) > J + kTolerance || std::abs(mp_values[j]) > J + kTolerance) { continue; }
      output[i][j] = rotation[rep.Index(mp_values[j], "gra::wigner::d")][rep.Index(m_values[i], "gra::wigner::d")].real();
    }
  }
  return output;
}

}  // namespace

// Compute whether x is an integer within the spin-label tolerance
bool IsInt(double x) { return std::isfinite(x) && std::abs(x - std::round(x)) < kTolerance; }

// Compute one SU(2) Clebsch-Gordan coefficient
// <j1 m1,j2 m2|j m> = (-1)^(m+j1-j2) sqrt(2j+1) W3j(j1,j2,j,m1,m2,-m)
double CG(double j1, double j2, double m1, double m2, double j, double m) {
  if (!CGAllowed(j1, j2, m1, m2, j, m)) { return 0.0; }
  return Phase(m + j1 - j2, "gra::wigner::CG") * gra::math::msqrt(2.0 * j + 1.0) * W3j(j1, j2, j, m1, m2, -m);
}

// Compute one Wigner 3j symbol using the Racah formula
double W3j(double j1, double j2, double j3, double m1, double m2, double m3) {
  if (!W3jAllowed(j1, j2, j3, m1, m2, m3)) { return 0.0; }
  if (j1 + j2 + j3 + 1.0 > kMaxFiniteFactorial) {
    throw std::invalid_argument("gra::wigner::W3j: spin sum exceeds the finite factorial range");
  }
  if (j1 + j2 + j3 > 32.0) {
    return Phase(j1 - j2 - m3, "gra::wigner::W3j") * CoupledCG(j1, j2, m1, j3, -m3) / std::sqrt(2.0 * j3 + 1.0);
  }

  const int c1    = static_cast<int>(std::llround(-j2 + m1 + j3));
  const int c2    = static_cast<int>(std::llround(-j1 - m2 + j3));
  const int c3    = static_cast<int>(std::llround(j1 + j2 - j3));
  const int c4    = static_cast<int>(std::llround(j1 - m1));
  const int c5    = static_cast<int>(std::llround(j2 + m2));
  const int lower = std::max({0, -c1, -c2});
  const int upper = std::min({c3, c4, c5});

  double sum = 0.0;
  for (int k = lower; k <= upper; ++k) {
    sum += Phase(static_cast<double>(k), "gra::wigner::W3j") /
           (gra::math::factorial(k) * gra::math::factorial(c1 + k) * gra::math::factorial(c2 + k) *
            gra::math::factorial(c3 - k) * gra::math::factorial(c4 - k) * gra::math::factorial(c5 - k));
  }

  const double result =
      Phase(j1 - j2 - m3, "gra::wigner::W3j") *
      gra::math::msqrt(Triangle(j1, j2, j3) * Factorial(j1 + m1, "gra::wigner::W3j") *
                       Factorial(j1 - m1, "gra::wigner::W3j") * Factorial(j2 + m2, "gra::wigner::W3j") *
                       Factorial(j2 - m2, "gra::wigner::W3j") * Factorial(j3 + m3, "gra::wigner::W3j") *
                       Factorial(j3 - m3, "gra::wigner::W3j")) *
      sum;
  if (!std::isfinite(result)) { throw std::overflow_error("gra::wigner::W3j: coefficient exceeds numerical range"); }
  return result;
}

// Compute Gamma continued Wigner 3j function
// [REFERENCE: A R White, Nucl Phys B67 (1973) 189]
std::complex<double> W3jRegge(const double j1, const double j2, const double j3, const double m1, const double m2,
                              const double m3) {
  if (!gra::AllFinite(std::array{j1, j2, j3, m1, m2, m3})) { throw AmplitudeFailure("wigner::W3jRegge: non-finite Gamma arguments"); }
  if (std::abs(m1 + m2 + m3) > kTolerance) { return 0.0; }
  if (W3jAllowed(j1, j2, j3, m1, m2, m3)) { return W3j(j1, j2, j3, m1, m2, m3); }

  const double                 c1             = -j2 + m1 + j3;
  const double                 c2             = -j1 - m2 + j3;
  const double                 c3             = j1 + j2 - j3;
  const double                 c4             = j1 - m1;
  const double                 c5             = j2 + m2;
  const std::array<double, 10> gamma_argument = {
      j1 + j2 - j3 + 1.0, j1 - j2 + j3 + 1.0, -j1 + j2 + j3 + 1.0, j1 + m1 + 1.0, j1 - m1 + 1.0,
      j2 + m2 + 1.0,      j2 - m2 + 1.0,      j3 + m3 + 1.0,       j3 - m3 + 1.0, j1 + j2 + j3 + 2.0};
  double log_abs       = 0.0;
  int    radicand_sign = 1;
  for (std::size_t index = 0; index + 1 < gamma_argument.size(); ++index) {
    const auto [value, gamma_sign] = SignedLogGamma(gamma_argument[index]);
    if (gamma_sign == 0) { return 0.0; }
    log_abs += 0.5 * value;
    if (gamma_sign < 0) { radicand_sign *= -1; }
  }
  const auto [triangle_den, triangle_sign] = SignedLogGamma(gamma_argument.back());
  if (triangle_sign == 0) { return 0.0; }
  log_abs -= 0.5 * triangle_den;
  if (triangle_sign < 0) { radicand_sign *= -1; }
  // At zero projections Dixon's sum replaces the infinite Racah series
  // Gamma duplication removes its coincident numerator and denominator poles
  // [REFERENCE: NIST DLMF 16.4.4 and 5.5.5, https://dlmf.nist.gov/16.4.E4, https://dlmf.nist.gov/5.5.E5]
  const double hyper =
      math::IsZero(m1) && math::IsZero(m2) && math::IsZero(m3) && j1 + j2 + j3 > -2.0
          ? math::FiniteSpecialValue(std::sqrt(static_cast<long double>(math::PI)) *
                                     std::exp2(static_cast<long double>(c3)) *
                                     math::GammaRatio(std::array{1.0 + 0.5 * (j1 + j2 + j3)},
                                                      std::array{0.5 * (1.0 - c3), 1.0 + 0.5 * (j1 - j2 + j3),
                                                                 1.0 + 0.5 * (-j1 + j2 + j3), 1.0 + j3}))
          : math::RegularizedHyper3F2Unit(-c3, -c4, -c5, c1 + 1.0, c2 + 1.0);
  if (std::fpclassify(hyper) == FP_ZERO) { return 0.0; }
  int outside_sign = hyper < 0.0 ? -1 : 1;
  for (const double value : {c3 + 1.0, c4 + 1.0, c5 + 1.0}) {
    const auto [denominator, gamma_sign] = SignedLogGamma(value);
    if (gamma_sign == 0) { return 0.0; }
    log_abs -= denominator;
    if (gamma_sign < 0) { outside_sign *= -1; }
  }
  log_abs += std::log(std::abs(hyper));
  const std::complex<double> phase      = std::exp(gra::math::zi * gra::math::PI * (j1 - j2 - m3));
  const std::complex<double> root_phase = radicand_sign < 0 ? gra::math::zi : std::complex<double>(1.0, 0.0);
  const std::complex<double> value      = phase * root_phase * static_cast<double>(outside_sign) * std::exp(log_abs);
  return value;
}

// Compute Gamma continued Clebsch-Gordan function
// [REFERENCE: A R White, Nucl Phys B67 (1973) 189]
std::complex<double> CGRegge(const double j1, const double j2, const double m1, const double m2, const double j,
                             const double m) {
  if (std::abs(m1 + m2 - m) > kTolerance) { return 0.0; }
  if (CGAllowed(j1, j2, m1, m2, j, m)) { return CG(j1, j2, m1, m2, j, m); }
  const std::complex<double> phase = std::exp(gra::math::zi * gra::math::PI * (m + j1 - j2));
  return phase * std::sqrt(std::complex<double>(2.0 * j + 1.0, 0.0)) * W3jRegge(j1, j2, j, m1, m2, -m);
}

// Compute one real Wigner element with arguments ordered as m,mp
// d^J_(mp,m)(theta) = <J,mp|exp(-i theta Jy)|J,m>
double d(double theta, double m, double mp, double J) {
  if (!std::isfinite(theta) || !std::isfinite(m) || !std::isfinite(mp)) {
    throw std::invalid_argument("gra::wigner::d: angles and spin projections must be finite");
  }
  const int J2 = SpinX2(J, "gra::wigner::d");
  if (std::abs(m) > J + kTolerance || std::abs(mp) > J + kTolerance) { return 0.0; }
  if (J2 == 0) { return 1.0; }
  if (J2 > 32) { return SpinorDRows(theta, {m}, {mp}, J)[0][0]; }
  const std::vector<double> cosine_power = math::IntegerPowers(std::cos(theta / 2.0), static_cast<unsigned int>(J2));
  const std::vector<double> sine_power   = math::IntegerPowers(std::sin(theta / 2.0), static_cast<unsigned int>(J2));
  return SmallDFromPowers(PreparePolynomial(m, mp, J, J2), cosine_power, sine_power);
}

// Compute Gamma continued Wigner d function
// [REFERENCE: NIST DLMF 15.2, https://dlmf.nist.gov/15.2]
std::complex<double> dRegge(const double theta, double m, double mp, const double J) {
  if (!gra::AllFinite(std::array{theta, m, mp, J})) { throw AmplitudeFailure("wigner::dRegge: non-finite Gamma arguments"); }
  if (!IsInt(m) || !IsInt(mp)) { return 0.0; }
  if (DAllowed(J, m, mp)) { return d(theta, m, mp, J); }
  
  int phase = 1;
  if (mp < m) {
    phase = (std::llround(m - mp) % 2 == 0) ? phase : -phase;
    std::swap(m, mp);
  }
  if (mp + m < 0.0) {
    const double old_m = m;
    m                  = -mp;
    mp                 = -old_m;
  }
  const int delta              = static_cast<int>(std::llround(mp - m));
  const int sum                = static_cast<int>(std::llround(mp + m));
  const auto [log_ap, sign_ap] = SignedLogGamma(J + mp + 1.0);
  const auto [log_am, sign_am] = SignedLogGamma(J - mp + 1.0);
  const auto [log_bp, sign_bp] = SignedLogGamma(J + m + 1.0);
  const auto [log_bm, sign_bm] = SignedLogGamma(J - m + 1.0);
  if (sign_ap == 0 || sign_am == 0 || sign_bp == 0 || sign_bm == 0) { return 0.0; }
  const int    ratio_sign    = sign_ap * sign_am * sign_bp * sign_bm;
  const int    jacobi_sign   = sign_bm * sign_am;
  const double log_prefactor = 0.5 * (log_ap + log_am - log_bp - log_bm) + log_bm - log_am;
  const double sine          = std::sin(theta / 2.0);
  const double cosine        = std::cos(theta / 2.0);
  if (sine < 0.0 || cosine < 0.0) { return 0.0; }
  const double z     = sine * sine;
  const double hyper = gra::math::RegularizedHyper2F1(-J + mp, J + mp + 1.0, static_cast<double>(delta + 1), z);
  const int                  base_phase = delta % 2 == 0 ? 1 : -1;
  const std::complex<double> root_phase = ratio_sign < 0 ? gra::math::zi : std::complex<double>(1.0, 0.0);
  const std::complex<double> value      = static_cast<double>(phase * base_phase * jacobi_sign) * root_phase *
                                     std::exp(log_prefactor) * std::pow(sine, delta) * std::pow(cosine, sum) * hyper;
  return value;
}

// Validate one fixed spin basis and prepare its angular coefficients
Rotation::Rotation(const std::vector<double> &m, const std::vector<double> &mp, double J)
    : spin_x2(SpinX2(J, "gra::wigner::Rotation")), representation(spin_x2, "gra::wigner::Rotation"), rows(m), columns(mp) {
  for (const auto *labels : {&rows, &columns}) {
    for (const double projection : *labels) {
      if (!std::isfinite(projection) || !IsInt(J - projection)) {
        throw std::invalid_argument("gra::wigner::Rotation: invalid spin projection");
      }
    }
  }
  if (spin_x2 > 32) {
    for (const double row : rows) {
      row_indices.push_back(std::abs(row) > J + kTolerance ? representation.Dim() : representation.Index(row, "wigner::Rotation"));
    }
    for (const double column : columns) {
      column_indices.push_back(std::abs(column) > J + kTolerance ? representation.Dim() : representation.Index(column, "wigner::Rotation"));
    }
    return;
  }
  coefficients.reserve(rows.size() * columns.size());
  for (const double row : rows) {
    for (const double column : columns) { coefficients.push_back(PreparePolynomial(row, column, J, spin_x2)); }
  }
}

// Evaluate cached low-spin polynomials or the stable large-spin recursion
MMatrix<double> Rotation::Evaluate(double theta) const {
  if (spin_x2 > 32) {
    const auto rotation = representation.Rotation(spin::SpinHalfRotation(theta, 0.0));
    MMatrix<double> output(rows.size(), columns.size(), 0.0);
    for (const auto &row : indices(rows)) {
      for (const auto &col : indices(columns)) {
        if (row_indices[row] < representation.Dim() && column_indices[col] < representation.Dim()) {
          output[row][col] = rotation[column_indices[col]][row_indices[row]].real();
        }
      }
    }
    return output;
  }
  const auto cosine_power = math::IntegerPowers(std::cos(theta / 2.0), static_cast<unsigned int>(spin_x2));
  const auto sine_power = math::IntegerPowers(std::sin(theta / 2.0), static_cast<unsigned int>(spin_x2));
  MMatrix<double> output(rows.size(), columns.size(), 0.0);
  for (const auto &row : indices(rows)) {
    for (const auto &col : indices(columns)) {
      output[row][col] = SmallDFromPowers(coefficients[row * columns.size() + col], cosine_power, sine_power);
    }
  }
  return output;
}

// Compute real small-d rows for a standalone spin basis
MMatrix<double> dRows(double theta, const std::vector<double> &m_values, const std::vector<double> &mp_values,
                      double J) {
  return Rotation(m_values, mp_values, J).Evaluate(theta);
}

// Compute the conjugate Wigner D element with arguments ordered as m,mp
// D^(J*)_(mp,m)(phi,theta,0) = d^J_(mp,m)(theta) exp(i mp phi)
std::complex<double> D(double theta, double phi, double m, double mp, double J) {
  if (!std::isfinite(phi)) { throw std::invalid_argument("gra::wigner::D: phi must be finite"); }
  return d(theta, m, mp, J) * std::exp(std::complex<double>(0.0, mp * phi));
}

// Compute one Gamma-continued conjugate Wigner D function
// D^(J*)_(mp,m) = d^J_(mp,m)(theta) exp(i mp phi)
std::complex<double> DRegge(const double theta, const double phi, const double m, const double mp, const double J) {
  return dRegge(theta, m, mp, J) * std::exp(std::complex<double>(0.0, mp * phi));
}

// Compute the stored transpose matrix with rows m and columns mp
// matrix_(m,mp) = D^(J*)_(mp,m)
MMatrix<std::complex<double>> DMatrix(double J, double theta, double phi) {
  if (!std::isfinite(phi)) { throw std::invalid_argument("gra::wigner::DMatrix: phi must be finite"); }
  const int           J2 = SpinX2(J, "gra::wigner::DMatrix");
  const std::size_t   n  = static_cast<std::size_t>(J2 + 1);
  std::vector<double> projections(n);
  for (std::size_t index = 0; index < n; ++index) { projections[index] = -J + static_cast<double>(index); }
  const MMatrix<double>         rotation = dRows(theta, projections, projections, J);
  MMatrix<std::complex<double>> matrix(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const double mp = -J + static_cast<double>(j);
      matrix[i][j]    = rotation[i][j] * std::exp(std::complex<double>(0.0, mp * phi));
    }
  }
  return matrix;
}

}  // namespace wigner
}  // namespace gra
