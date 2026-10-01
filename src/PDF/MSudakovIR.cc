// Infrared continuation of the skewed Durham gluon flux
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <stdexcept>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/PDF/MSudakov.h"

using gra::aux::indices;

namespace gra {

using math::pow2;

namespace {

// Compute the weighted squared-exponential kernel before boundary conditioning
// K(s,a) = exp[-w(s+a) -(s-a)^2/(2 ell^2)]
long double WeightedGaussianKernel(long double s, long double anchor, long double length,
                                   long double exponent) {
  const long double difference = (s - anchor) / length;
  return std::exp(-exponent * (s + anchor) - 0.5L * difference * difference);
}

// Compute the weighted squared-exponential kernel conditioned on value and slope at zero
long double ConditionedGaussianKernel(long double s, long double anchor, long double length,
                                      long double exponent) {
  const long double base = WeightedGaussianKernel(s, anchor, length, exponent);
  const long double c0_s = WeightedGaussianKernel(s, 0.0L, length, exponent);
  const long double c0_a = WeightedGaussianKernel(anchor, 0.0L, length, exponent);
  const long double inverse_length2 = 1.0L / (length * length);
  const long double c1_s = c0_s * (-exponent + s * inverse_length2);
  const long double c1_a = c0_a * (-exponent + anchor * inverse_length2);
  const long double inverse00 = 1.0L + exponent * exponent * length * length;
  const long double inverse01 = exponent * length * length;
  const long double inverse11 = length * length;
  const long double projection =
      c0_s * (inverse00 * c0_a + inverse01 * c1_a) +
      c1_s * (inverse01 * c0_a + inverse11 * c1_a);
  return base - projection;
}

// Compute one finite exponential moment
// integral_0^u s^n exp(-r s) ds = gamma(n+1,ru)/r^(n+1)
long double ExponentialMoment(unsigned int power, long double rate,
                              long double upper) {
  if (upper <= 0.0L) { return 0.0L; }
  long double factorial = 1.0L;
  for (unsigned int i = 2; i <= power; ++i) {
    factorial *= static_cast<long double>(i);
  }
  const long double normalization =
      factorial / std::pow(rate, static_cast<int>(power + 1));
  if (std::isinf(upper)) { return normalization; }

  const long double z = rate * upper;
  if (z < 0.05L) {
    long double term =
        std::pow(upper, static_cast<int>(power + 1)) /
        static_cast<long double>(power + 1);
    long double sum = term;
    for (unsigned int order = 0; order < 32; ++order) {
      term *= -z / static_cast<long double>(order + 1) *
              static_cast<long double>(power + order + 1) /
              static_cast<long double>(power + order + 2);
      sum += term;
      if (std::abs(term) <=
          std::numeric_limits<long double>::epsilon() * std::abs(sum)) {
        break;
      }
    }
    return sum;
  }

  long double series = 1.0L;
  long double term = 1.0L;
  for (unsigned int order = 1; order <= power; ++order) {
    term *= z / static_cast<long double>(order);
    series += term;
  }
  return normalization * (1.0L - std::exp(-z) * series);
}

// Compute the positive colour-neutral baseline logarithmic slope
// h0(s) = [1+(c+d s) exp(-lambda s)]^2
long double InfraredBaselineSlope(long double s, long double c, long double d,
                                  long double lambda) {
  const long double profile =
      1.0L + (c + d * s) * std::exp(-lambda * s);
  return profile * profile;
}

// Compute the baseline integral in excess of the colour-neutral unit slope
long double InfraredBaselineExcessIntegral(long double s, long double c,
                                           long double d,
                                           long double lambda) {
  const long double i0 = ExponentialMoment(0, lambda, s);
  const long double i1 = ExponentialMoment(1, lambda, s);
  const long double i0_double = ExponentialMoment(0, 2.0L * lambda, s);
  const long double i1_double = ExponentialMoment(1, 2.0L * lambda, s);
  const long double i2_double = ExponentialMoment(2, 2.0L * lambda, s);
  return 2.0L * c * i0 + 2.0L * d * i1 +
         c * c * i0_double + 2.0L * c * d * i1_double +
         d * d * i2_double;
}

// Compute the conditioned RKHS point-evaluation norm at one infrared distance
long double ConditionedEvaluationNorm(long double s, long double length,
                                      long double exponent) {
  return std::sqrt(std::max(
      0.0L, ConditionedGaussianKernel(s, s, length, exponent)));
}

// Compute the positivity margin for the complete RKHS ball at one s
// m(s) = h0(s) - R sqrt(K_c(s,s))
long double RKHSBallPositivityMargin(long double s, long double c,
                                     long double d, long double lambda,
                                     long double length,
                                     long double exponent,
                                     long double radius) {
  return InfraredBaselineSlope(s, c, d, lambda) -
         radius * ConditionedEvaluationNorm(s, length, exponent);
}

// Refine one bracketed RKHS positivity-margin minimum
long double RefineRKHSBallMinimum(long double left, long double right,
                                  long double c, long double d,
                                  long double lambda, long double length,
                                  long double exponent, long double radius) {
  const long double ratio =
      (std::sqrt(5.0L) - 1.0L) / 2.0L;
  long double x1 = right - ratio * (right - left);
  long double x2 = left + ratio * (right - left);
  long double f1 = RKHSBallPositivityMargin(
      x1, c, d, lambda, length, exponent, radius);
  long double f2 = RKHSBallPositivityMargin(
      x2, c, d, lambda, length, exponent, radius);
  for (unsigned int iteration = 0; iteration < 80; ++iteration) {
    if (f1 <= f2) {
      right = x2;
      x2 = x1;
      f2 = f1;
      x1 = right - ratio * (right - left);
      f1 = RKHSBallPositivityMargin(
          x1, c, d, lambda, length, exponent, radius);
    } else {
      left = x1;
      x1 = x2;
      f1 = f2;
      x2 = left + ratio * (right - left);
      f2 = RKHSBallPositivityMargin(
          x2, c, d, lambda, length, exponent, radius);
    }
  }
  return std::min(f1, f2);
}

// Compute a finite limit beyond which a positive analytic tail bound holds
long double RKHSValidationUpper(long double c, long double d,
                                long double lambda, long double length,
                                long double exponent, long double radius) {
  long double upper =
      std::max({1.0L, 1.0L / lambda, length, 1.0L / exponent});
  for (unsigned int iteration = 0; iteration < 128; ++iteration) {
    const long double profile_bound =
        (std::abs(c) + std::abs(d) * upper) *
        std::exp(-lambda * upper);
    const long double baseline_bound =
        2.0L * profile_bound + profile_bound * profile_bound;
    const long double ball_bound = radius * std::exp(-exponent * upper);
    if (baseline_bound + ball_bound <= 0.25L) { return upper; }
    upper *= 2.0L;
    if (!std::isfinite(upper)) { break; }
  }
  throw std::domain_error(
      "MSudakov: unable to establish the infrared RKHS positivity tail");
}

// Compute the numerically refined full-ball positivity margin over s >= 0
long double RKHSBallMinimum(long double c, long double d,
                            long double lambda, long double length,
                            long double exponent, long double radius) {
  const long double upper =
      RKHSValidationUpper(c, d, lambda, length, exponent, radius);
  const long double shortest_scale =
      std::min({1.0L, 1.0L / lambda, length, 1.0L / exponent});
  const long double target_step = shortest_scale / 64.0L;
  const long double intervals_real = std::ceil(upper / target_step);
  if (!(intervals_real >= 2.0L) || intervals_real > 200000.0L) {
    throw std::domain_error(
        "MSudakov: infrared RKHS positivity validation is under-resolved");
  }
  const std::size_t intervals = static_cast<std::size_t>(intervals_real);
  const long double step = upper / static_cast<long double>(intervals);

  long double previous = RKHSBallPositivityMargin(
      0.0L, c, d, lambda, length, exponent, radius);
  long double current = RKHSBallPositivityMargin(
      step, c, d, lambda, length, exponent, radius);
  long double minimum = std::min(previous, current);
  for (std::size_t i = 1; i < intervals; ++i) {
    const long double next_s =
        static_cast<long double>(i + 1) * step;
    const long double next = RKHSBallPositivityMargin(
        next_s, c, d, lambda, length, exponent, radius);
    minimum = std::min(minimum, next);
    if (current <= previous && current <= next) {
      const long double left =
          static_cast<long double>(i - 1) * step;
      minimum = std::min(
          minimum, RefineRKHSBallMinimum(
                       left, next_s, c, d, lambda, length, exponent, radius));
    }
    previous = current;
    current = next;
  }
  return minimum;
}

// Validate positivity of the baseline and every member in the configured ball
void ValidateInfraredPositivity(const MSudakovModel &model, double c, double d,
                                double lambda) {
  const long double c_ld = c;
  const long double d_ld = d;
  const long double lambda_ld = lambda;
  const long double scale =
      std::max({1.0L, std::abs(c_ld), std::abs(d_ld) / lambda_ld});
  const long double tolerance =
      4096.0L * std::numeric_limits<long double>::epsilon() * scale;
  if (model.Mode() != MSudakovIRMode::IR_RKHS ||
      std::fpclassify(model.rkhs_radius) == FP_ZERO) {
    return;
  }

  const long double margin = RKHSBallMinimum(
      c_ld, d_ld, lambda_ld, model.rkhs_logq2_length,
      model.rkhs_weight_exponent, model.rkhs_radius);
  if (margin < -tolerance) {
    throw std::domain_error(
        "MSudakov: configured RKHS ball does not preserve positive flux");
  }
}

// Compute the weighted squared-exponential kernel integral from zero to upper
long double WeightedGaussianKernelIntegral(long double upper, long double anchor,
                                           long double length, long double exponent) {
  if (upper <= 0.0L) { return 0.0L; }
  const long double root_two = std::sqrt(2.0L);
  const long double shifted = anchor - exponent * length * length;
  const long double lower_argument = -shifted / (root_two * length);
  const long double upper_argument = (upper - shifted) / (root_two * length);
  // Avoid subtracting two erf values rounded to the same Gaussian tail limit
  const long double difference =
      lower_argument >= 0.0L
          ? std::erfc(lower_argument) - std::erfc(upper_argument)
          : (upper_argument <= 0.0L
                 ? std::erfc(-upper_argument) - std::erfc(-lower_argument)
                 : std::erf(upper_argument) - std::erf(lower_argument));
  const long double normalization =
      std::exp(-2.0L * exponent * anchor +
               0.5L * exponent * exponent * length * length) *
      length * std::sqrt(std::acos(-1.0L) / 2.0L);
  return normalization * difference;
}

// Compute the conditioned weighted kernel integral from zero to upper
long double ConditionedGaussianKernelIntegral(long double upper, long double anchor,
                                              long double length, long double exponent) {
  if (upper <= 0.0L) { return 0.0L; }
  const long double base =
      WeightedGaussianKernelIntegral(upper, anchor, length, exponent);
  const long double integral_c0 =
      WeightedGaussianKernelIntegral(upper, 0.0L, length, exponent);
  const long double boundary_c0 =
      std::isinf(upper)
          ? 0.0L
          : WeightedGaussianKernel(upper, 0.0L, length, exponent);
  const long double integral_c1 =
      1.0L - boundary_c0 - 2.0L * exponent * integral_c0;
  const long double c0_a =
      WeightedGaussianKernel(anchor, 0.0L, length, exponent);
  const long double c1_a =
      c0_a * (-exponent + anchor / (length * length));
  const long double inverse00 = 1.0L + exponent * exponent * length * length;
  const long double inverse01 = exponent * length * length;
  const long double inverse11 = length * length;
  const long double projection =
      integral_c0 * (inverse00 * c0_a + inverse01 * c1_a) +
      integral_c1 * (inverse01 * c0_a + inverse11 * c1_a);
  return base - projection;
}

// Compute the squared-exponential correlation kernel in logarithmic x
// K_x(x,a) = exp[-ln^2(x/a)/(2 ell_x^2)]
long double GaussianLogXKernel(long double x, long double anchor, long double length) {
  const long double difference = std::log(x / anchor) / length;
  return std::exp(-0.5L * difference * difference);
}

// Compute a stable factorized RKHS norm for one finite representer member
double WeightedMemberNorm2(const std::vector<double> &anchors,
                           const std::vector<double> &x_anchors,
                           const std::vector<double> &coefficients,
                           double length, double exponent,
                           double logx_length) {
  const std::size_t n = coefficients.size();
  if (anchors.size() != n || x_anchors.size() != n) {
    throw std::invalid_argument("MSudakovModel::MemberNorm2: anchors and coefficients must have equal sizes");
  }
  if (n == 0) { return 0.0; }

  MMatrix<long double> kernel(n, n, 0.0L);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j <= i; ++j) {
      const long double value =
          ConditionedGaussianKernel(anchors[i], anchors[j], length, exponent) *
          GaussianLogXKernel(x_anchors[i], x_anchors[j], logx_length);
      kernel[i][j] = value;
      kernel[j][i] = value;
    }
  }

  MMatrix<long double> lower;
  try {
    constexpr long double tolerance =
        4096.0L * std::numeric_limits<long double>::epsilon();
    lower =
        kernel.CholeskyLower(tolerance, CholeskyMode::StrictPositiveDefinite);
  } catch (const std::invalid_argument &) {
    throw std::invalid_argument(
        "MSudakovModel::MemberNorm2: ill-conditioned member anchors");
  } catch (const std::domain_error &) {
    throw std::invalid_argument(
        "MSudakovModel::MemberNorm2: ill-conditioned member anchors");
  } catch (const std::runtime_error &) {
    throw std::invalid_argument(
        "MSudakovModel::MemberNorm2: ill-conditioned member anchors");
  }

  std::vector<long double> coefficients_ld(coefficients.begin(),
                                           coefficients.end());
  const std::vector<long double> projected =
      lower.LeftMultiply(coefficients_ld);
  return static_cast<double>(gra::SquaredNorm(projected));
}

}  // namespace

// Compute the parsed infrared mode
MSudakovIRMode MSudakovModel::Mode() const {
  if (mode_name == "PERTURBATIVE_ONLY") { return MSudakovIRMode::PERTURBATIVE_ONLY; }
  if (mode_name == "IR_BASELINE") { return MSudakovIRMode::IR_BASELINE; }
  if (mode_name == "IR_RKHS") { return MSudakovIRMode::IR_RKHS; }
  throw std::invalid_argument("MSudakovModel::Mode: unknown mode " + mode_name);
}

// Compute the squared weighted-RKHS norm of the selected functional member
double MSudakovModel::MemberNorm2() const {
  return WeightedMemberNorm2(
      member_s_anchors, member_x_anchors, member_coefficients,
      rkhs_logq2_length, rkhs_weight_exponent, rkhs_logx_length);
}

// Validate the matched-boundary physics parameters
void MSudakovModel::Validate() const {
  (void)Mode();
  if (!std::isfinite(Q0) || Q0 <= 0.0) {
    throw std::invalid_argument("MSudakovModel::Validate: Q0 must be positive");
  }
  if (!std::isfinite(ir_logq2_length) || ir_logq2_length <= 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: ir_logq2_length must be positive");
  }
  if (!std::isfinite(rkhs_logq2_length) || rkhs_logq2_length <= 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: rkhs_logq2_length must be positive");
  }
  if (!std::isfinite(rkhs_weight_exponent) || rkhs_weight_exponent <= 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: rkhs_weight_exponent must be positive");
  }
  if (!std::isfinite(rkhs_logx_length) || rkhs_logx_length <= 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: rkhs_logx_length must be positive");
  }
  if (!std::isfinite(rkhs_radius) || rkhs_radius < 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: rkhs_radius must be non-negative");
  }
  if (member_s_anchors.size() != member_coefficients.size()) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: member s anchors and coefficients must have equal sizes");
  }
  if (member_x_anchors.size() != member_coefficients.size()) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: member x anchors and coefficients must have equal sizes");
  }
  for (const auto &i : indices(member_s_anchors)) {
    const double anchor = member_s_anchors[i];
    if (!std::isfinite(anchor) || anchor <= 0.0) {
      throw std::invalid_argument(
          "MSudakovModel::Validate: functional anchors must be positive");
    }
    const double x_anchor = member_x_anchors[i];
    if (!std::isfinite(x_anchor) || !(x_anchor > 0.0 && x_anchor < 1.0)) {
      throw std::invalid_argument(
          "MSudakovModel::Validate: functional x anchors must lie in (0,1)");
    }
    for (std::size_t j = 0; j < i; ++j) {
      const double s_separation = std::abs(anchor - member_s_anchors[j]);
      const double x_separation =
          std::abs(std::log(x_anchor / member_x_anchors[j]));
      if (s_separation <=
              1e-8 * std::max({1.0, std::abs(anchor),
                               std::abs(member_s_anchors[j])}) &&
          x_separation <= 1e-8) {
        throw std::invalid_argument(
            "MSudakovModel::Validate: functional anchor pairs must be separated");
      }
    }
  }
  for (const double coefficient : member_coefficients) {
    if (!std::isfinite(coefficient)) {
      throw std::invalid_argument(
          "MSudakovModel::Validate: functional coefficients must be finite");
    }
  }
  const double tolerance = 256.0 * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, pow2(rkhs_radius));
  if (MemberNorm2() > pow2(rkhs_radius) + tolerance) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: selected functional member lies outside the weighted-RKHS ball");
  }
  if (Mode() != MSudakovIRMode::IR_RKHS && MemberNorm2() > tolerance) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: non-zero member requires IR_RKHS");
  }
  if (Mode() != MSudakovIRMode::IR_RKHS && rkhs_radius > 0.0) {
    throw std::invalid_argument(
        "MSudakovModel::Validate: non-zero rkhs_radius requires IR_RKHS");
  }
}

// Read matched-boundary physics from one immutable GENERAL JSON text
void MSudakovModel::ConfigureFromJson(
    const std::string &source_file, const std::string &json_text) {
  using json = nlohmann::json;
  try {
    const json j = json::parse(json_text);
    const auto &block = j.at("PARAM_SKEWED_UGD");
    mode_name         = block.at("mode");
    Q0                   = block.at("Q0");
    ir_logq2_length      = block.at("ir_logq2_length");
    rkhs_logq2_length   = block.at("rkhs_logq2_length");
    rkhs_logx_length    = block.at("rkhs_logx_length");
    rkhs_weight_exponent = block.at("rkhs_weight_exponent");
    rkhs_radius = block.at("rkhs_radius");
    member_s_anchors =
        block.at("member_s_anchors").get<std::vector<double>>();
    member_x_anchors =
        block.at("member_x_anchors").get<std::vector<double>>();
    member_coefficients =
        block.at("member_coefficients").get<std::vector<double>>();
    Validate();

    std::cout << "MSudakovModel::ReadParameters: [PARAM_SKEWED_UGD]" << std::endl;
    std::cout << block << std::endl << std::endl;
  } catch (const nlohmann::json::exception &e) {
    throw std::invalid_argument(
        "MSudakovModel::ConfigureFromJson: Error parsing " +
        source_file + " (Check PARAM_SKEWED_UGD): " + e.what());
  }
}

// Compute the tensor-product Green kernel used for one functional member
double MSudakov::FunctionalKernel(double s, double anchor, double x,
                                  double x_anchor) const {
  return static_cast<double>(ConditionedGaussianKernel(
      s, anchor, Model.rkhs_logq2_length, Model.rkhs_weight_exponent) *
      GaussianLogXKernel(x, x_anchor, Model.rkhs_logx_length));
}

// Compute the tensor-product Green kernel integral from the matching point to s
double MSudakov::FunctionalKernelIntegral(double s, double anchor, double x,
                                          double x_anchor) const {
  return static_cast<double>(ConditionedGaussianKernelIntegral(
      s, anchor, Model.rkhs_logq2_length, Model.rkhs_weight_exponent) *
      GaussianLogXKernel(x, x_anchor, Model.rkhs_logx_length));
}

// Compute one selected admissible infrared continuation
//
// For s = log(Q0^2/Q^2), h = f/G and lambda = 1/ell:
//
// h0(s) = [1 + (c + d s) exp(-lambda s)]^2
// G(s)  = G0 exp[-int_0^s (h0 + v) du]
// f(s)  = G(s) [h0(s)+v(s)]
//
// The conditioned RKHS field obeys v(0) = v'(0) = 0 and decays at large s
// [REFERENCE: Martin and Ryskin, Phys. Rev. D64 (2001) 094017, hep-ph/0107149]
MSudakov::MatchedFlux MSudakov::InfraredFlux(double x, double q2, double mu) const {
  if (Model.Mode() == MSudakovIRMode::PERTURBATIVE_ONLY || !(q2 > 0.0)) { return {}; }

  const MatchedFlux boundary = PerturbativeFlux(x, Numerics.q2_MIN, mu);
  if (!(boundary.integrated > 0.0) || !std::isfinite(boundary.integrated)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: positive perturbative boundary value is required");
  }
  const double gamma = boundary.flux / boundary.integrated;
  if (!(gamma > 0.0) || !std::isfinite(gamma)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: positive perturbative boundary slope is required");
  }

  // Continue the infrared profile from the radiating side of mu = Q0
  const double match_mu = std::max(mu, std::nextafter(Numerics.mu_MIN, Numerics.mu_MAX));
  const double boundary_derivative = PerturbativeFluxLogDerivative(x, Numerics.q2_MIN, match_mu);
  const double kappa =
      gamma * gamma - boundary_derivative / boundary.integrated;
  const double s = std::log(Numerics.q2_MIN) - std::log(q2);
  const double lambda = 1.0 / Model.ir_logq2_length;
  const double root_gamma = std::sqrt(gamma);
  const double c = root_gamma - 1.0;
  const double d = kappa / (2.0 * root_gamma) + lambda * c;
  ValidateInfraredPositivity(Model, c, d, lambda);
  const double decay = std::exp(-lambda * s);
  const double profile = 1.0 + (c + d * s) * decay;
  const double h_base = profile * profile;
  double log_shape =
      -s - static_cast<double>(
               InfraredBaselineExcessIntegral(s, c, d, lambda));
  double h = h_base;

  if (Model.Mode() == MSudakovIRMode::IR_RKHS) {
    const double evaluation_norm = std::sqrt(std::max(
        0.0, static_cast<double>(ConditionedGaussianKernel(
                 s, s, Model.rkhs_logq2_length,
                 Model.rkhs_weight_exponent))));
    if (h_base + 1.0e-12 < Model.rkhs_radius * evaluation_norm) {
      throw std::domain_error(
          "MSudakov::InfraredFlux: RKHS ball violates positive infrared flux");
    }
    for (const auto &i : indices(Model.member_coefficients)) {
      h += Model.member_coefficients[i] *
           FunctionalKernel(s, Model.member_s_anchors[i], x,
                            Model.member_x_anchors[i]);
      log_shape -= Model.member_coefficients[i] *
                   FunctionalKernelIntegral(
                       s, Model.member_s_anchors[i], x,
                       Model.member_x_anchors[i]);
    }
  }
  if (!(h >= 0.0) || !std::isfinite(h)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: non-positive functional logarithmic slope");
  }
  if (log_shape < std::log(std::numeric_limits<double>::min())) { return {}; }

  const double integrated = boundary.integrated * std::exp(log_shape);
  return {integrated, integrated * h};
}

// Prepare the infrared continuation parameters and analytic zero limit
void MSudakov::PreparedFlux::PrepareInfrared() {
  infrared_mode = owner->Model.Mode();
  infrared_lambda = 1.0 / owner->Model.ir_logq2_length;
  if (infrared_mode == MSudakovIRMode::PERTURBATIVE_ONLY) {
    infrared_ready = true;
    return;
  }
  const auto boundary = owner->PerturbativeFlux(x, q2_min, mu);
  boundary_integrated = boundary.integrated;
  boundary_flux       = boundary.flux;
  if (!(boundary_integrated > 0.0) || !std::isfinite(boundary_integrated)) {
    throw std::domain_error("MSudakov::PreparedFlux: positive boundary gluon is required");
  }
  const double gamma = boundary_flux / boundary_integrated;
  if (!(gamma > 0.0) || !std::isfinite(gamma)) {
    throw std::domain_error(
        "MSudakov::PreparedFlux::PrepareInfrared: positive boundary slope is required");
  }
  // Continue the infrared profile from the radiating side of mu = Q0
  const double match_mu = std::max(mu, std::nextafter(owner->Numerics.mu_MIN, owner->Numerics.mu_MAX));
  const double boundary_derivative = owner->PerturbativeFluxLogDerivative(x, q2_min, match_mu);
  const double kappa =
      gamma * gamma - boundary_derivative / boundary_integrated;
  const double root_gamma = std::sqrt(gamma);
  infrared_c = root_gamma - 1.0;
  infrared_d =
      kappa / (2.0 * root_gamma) + infrared_lambda * infrared_c;
  boundary_over_q2_min = boundary_integrated / q2_min;
  ValidateInfraredPositivity(
      owner->Model, infrared_c, infrared_d, infrared_lambda);

  double exponent = -static_cast<double>(InfraredBaselineExcessIntegral(
      std::numeric_limits<long double>::infinity(), infrared_c, infrared_d,
      infrared_lambda));
  if (infrared_mode == MSudakovIRMode::IR_RKHS) {
    for (const auto &i : indices(owner->Model.member_coefficients)) {
      const double integral_inf = owner->FunctionalKernelIntegral(
          std::numeric_limits<double>::infinity(),
          owner->Model.member_s_anchors[i], x,
          owner->Model.member_x_anchors[i]);
      exponent -= owner->Model.member_coefficients[i] * integral_inf;
    }
  }
  zero_limit = boundary_over_q2_min * std::exp(exponent);
  infrared_ready = true;
}

// Compute the infrared continuation contribution to f_g/Q2
double MSudakov::PreparedFlux::InfraredOverQ2(double q2) const {
  if (infrared_mode == MSudakovIRMode::PERTURBATIVE_ONLY) {
    return 0.0;
  }
  if (std::fpclassify(q2) == FP_ZERO) { return zero_limit; }
  if (!(boundary_integrated > 0.0) || !std::isfinite(boundary_integrated)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: positive perturbative boundary value is required");
  }
  if (!std::isfinite(infrared_c) || !std::isfinite(infrared_d)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: finite perturbative boundary slope is required");
  }

  const double s = std::log(q2_min) - std::log(q2);
  const double decay = std::exp(-infrared_lambda * s);
  const double profile =
      1.0 + (infrared_c + infrared_d * s) * decay;
  const double h_base = profile * profile;
  double h = h_base;
  // Cancel exp(-s)/Q2 = 1/Q0^2 before exponentiation to preserve the finite limit
  double log_shape =
      -static_cast<double>(InfraredBaselineExcessIntegral(
               s, infrared_c, infrared_d, infrared_lambda));
  if (infrared_mode == MSudakovIRMode::IR_RKHS) {
    const double evaluation_norm = std::sqrt(std::max(
        0.0, static_cast<double>(ConditionedGaussianKernel(
                 s, s, owner->Model.rkhs_logq2_length,
                 owner->Model.rkhs_weight_exponent))));
    if (h_base + 1.0e-12 < owner->Model.rkhs_radius * evaluation_norm) {
      throw std::domain_error(
          "MSudakov::PreparedFlux: RKHS ball violates positive infrared flux");
    }
    for (const auto &i : indices(owner->Model.member_coefficients)) {
      h += owner->Model.member_coefficients[i] *
           owner->FunctionalKernel(
               s, owner->Model.member_s_anchors[i], x,
               owner->Model.member_x_anchors[i]);
      log_shape -= owner->Model.member_coefficients[i] *
                   owner->FunctionalKernelIntegral(
                       s, owner->Model.member_s_anchors[i], x,
                       owner->Model.member_x_anchors[i]);
    }
  }
  if (!(h >= 0.0) || !std::isfinite(h)) {
    throw std::domain_error(
        "MSudakov::InfraredFlux: non-positive functional logarithmic slope");
  }
  return boundary_over_q2_min * std::exp(log_shape) * h;
}

// Compute f_g / Q2 with the analytic color-neutral limit at Q2 = 0
double MSudakov::fg_xQ2MuOverQ2(double x, double q2, double mu) const {
  if (!std::isfinite(q2) || q2 < 0.0) { return 0.0; }
  if (q2 >= Numerics.q2_MIN && q2 > 0.0) {
    return fg_xQ2Mu(x, q2, mu) / q2;
  }
  return PrepareFlux(x, mu).OverQ2(q2);
}

// Compute a read-only clone for coherent functional propagation through observables
std::shared_ptr<const MSudakov> MSudakov::WithMemberCoefficients(
    const std::vector<double> &coefficients) const {
  if (Model.Mode() != MSudakovIRMode::IR_RKHS) {
    throw std::domain_error(
        "MSudakov::WithMemberCoefficients: IR_RKHS mode is not active");
  }
  std::shared_ptr<MSudakov> clone = std::make_shared<MSudakov>(*this);
  clone->Model.member_coefficients = coefficients;
  clone->Model.Validate();
  return clone;
}

}  // namespace gra
