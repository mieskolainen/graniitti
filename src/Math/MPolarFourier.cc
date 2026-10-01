// Non-radial polar Fourier-Bessel transforms and convolutions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

namespace gra::math {
namespace {

// Compute true when one complex value has finite real and imaginary parts
bool FiniteComplex(const std::complex<double> value) {
  return std::isfinite(value.real()) && std::isfinite(value.imag());
}

// Validate one radial quadrature carrying the measure r dr
void ValidatePolarRule(const PolarMeasureRule &rule,
                       const std::string &context) {
  if (rule.node.empty() || rule.node.size() != rule.measure_weight.size()) {
    throw std::invalid_argument(context + ": invalid radial rule dimensions");
  }
  for (const auto &i : indices(rule.node)) {
    if (!std::isfinite(rule.node[i]) || rule.node[i] < 0.0 || !std::isfinite(rule.measure_weight[i]) ||
        rule.measure_weight[i] < 0.0) {
      throw std::invalid_argument(context + ": invalid radial rule value");
    }
    if (i > 0 && !(rule.node[i - 1] < rule.node[i])) {
      throw std::invalid_argument(context + ": radial nodes are not ordered");
    }
  }
  if (std::none_of(rule.measure_weight.begin(), rule.measure_weight.end(),
                   [](const double weight) { return weight > 0.0; })) {
    throw std::invalid_argument(context + ": radial rule has zero measure");
  }
}

// Validate uniform angular nodes spanning one complete period
void ValidateAzimuthNodes(const std::vector<double> &node) {
  if (node.empty()) {
    throw std::invalid_argument(
        "MPolarFourier: angular quadrature must not be empty");
  }
  for (const double phi : node) {
    if (!std::isfinite(phi)) {
      throw std::invalid_argument(
          "MPolarFourier: non-finite angular quadrature node");
    }
  }
  if (node.size() == 1) {
    return;
  }
  const double step = 2.0 * PI / static_cast<double>(node.size());
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() *
                           std::max(1.0, std::abs(node.back()));
  for (std::size_t i = 1; i < node.size(); ++i) {
    if (std::abs((node[i] - node[i - 1]) - step) > tolerance) {
      throw std::invalid_argument(
          "MPolarFourier: angular quadrature is not uniform");
    }
  }
}

// Compute one non-negative integer-order Bessel function of the first kind
// output = J_order(argument)
double IntegerBessel(const std::size_t order, const double argument) {
  if (std::fpclassify(argument) == FP_ZERO) {
    return order == 0 ? 1.0 : 0.0;
  }
  return std::cyl_bessel_j(static_cast<double>(order), argument);
}

// Compute i raised to one non-negative integer power
// output = i^power
std::complex<double> PositiveImaginaryPhase(const std::size_t power) {
  switch (power % 4) {
  case 0:
    return {1.0, 0.0};
  case 1:
    return {0.0, 1.0};
  case 2:
    return {-1.0, 0.0};
  default:
    return {0.0, -1.0};
  }
}

// Compute the matrix column representing one signed harmonic
std::size_t HarmonicColumn(const int harmonic, const std::size_t max_harmonic) {
  if (harmonic < -static_cast<int>(max_harmonic) ||
      harmonic > static_cast<int>(max_harmonic)) {
    throw std::out_of_range("MPolarFourier: harmonic is outside the field");
  }
  return static_cast<std::size_t>(harmonic + static_cast<int>(max_harmonic));
}

// Validate one impact parameter harmonic field against the transform grid
void ValidateHarmonicField(const PolarHarmonicField &field,
                           const std::size_t radial_size) {
  if (field.max_harmonic >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(
        "MPolarFourier: harmonic order exceeds integer range");
  }
  if (field.coefficient.size_row() != radial_size ||
      field.coefficient.size_col() != 2 * field.max_harmonic + 1) {
    throw std::invalid_argument(
        "MPolarFourier: invalid harmonic field dimensions");
  }
}

// Validate one prepared inverse kernel against the transform grid
void ValidateInverseKernel(const PolarInverseKernel &kernel,
                           const std::size_t radial_size) {
  if (kernel.max_harmonic >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      kernel.radial_weight.size_row() != kernel.max_harmonic + 1 ||
      kernel.radial_weight.size_col() != radial_size ||
      kernel.angular_phase.size() != 2 * kernel.max_harmonic + 1) {
    throw std::invalid_argument(
        "MPolarFourier: invalid prepared inverse kernel dimensions");
  }
}

// Apply one validated inverse kernel to one validated harmonic field
// f(q,phi) = sum_m (-i)^|m| exp(i m phi) sum_b w_b J_|m|(qb) F_m(b)
std::complex<double> ApplyInverseKernel(const PolarHarmonicField &field,
                                        const PolarInverseKernel &kernel) {
  if (field.max_harmonic > kernel.max_harmonic) {
    throw std::invalid_argument(
        "MPolarFourier: inverse kernel has insufficient harmonics");
  }

  std::complex<double> out = 0.0;
  for (int harmonic = -static_cast<int>(field.max_harmonic);
       harmonic <= static_cast<int>(field.max_harmonic); ++harmonic) {
    const std::size_t order = static_cast<std::size_t>(std::abs(harmonic));
    const std::size_t field_column =
        HarmonicColumn(harmonic, field.max_harmonic);
    const std::size_t kernel_column =
        HarmonicColumn(harmonic, kernel.max_harmonic);
    std::complex<double> radial = 0.0;
    for (std::size_t ib = 0; ib < kernel.radial_weight.size_col(); ++ib) {
      radial +=
          kernel.radial_weight(order, ib) * field.coefficient(ib, field_column);
    }
    out += kernel.angular_phase[kernel_column] * radial;
  }
  if (!FiniteComplex(out)) {
    throw gra::AmplitudeFailure(
        "MPolarFourier::InverseEvaluate: non-finite transformed field");
  }
  return out;
}

// Compute a zero-initialized harmonic field on one radial grid
PolarHarmonicField ZeroHarmonicField(const std::size_t radial_size,
                                     const std::size_t max_harmonic) {
  PolarHarmonicField out;
  out.max_harmonic = max_harmonic;
  out.coefficient =
      MMatrix<std::complex<double>>(radial_size, 2 * max_harmonic + 1, 0.0);
  return out;
}

} // namespace

// Compute one signed harmonic from a periodic discrete Fourier grid
// m(k) = k for k <= N/2 and k-N otherwise
int SignedHarmonic(const unsigned int mode, const unsigned int count) {
  if (count == 0 || mode >= count ||
      count > static_cast<unsigned int>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument("SignedHarmonic: invalid periodic grid index");
  }
  return mode <= count / 2 ? static_cast<int>(mode)
                           : static_cast<int>(mode) - static_cast<int>(count);
}

// Prepare the Hankel transform by integrating each logarithmic slope analytically
// Integration by parts gives 2 pi [J0(q r0/s)-J0(q r1/s)] / (q^2 log(r1/r0))
MMatrix<double> LogBesselKernel(const std::vector<double>& radius, const std::vector<double>& momentum,
                               const double coordinate_scale) {
  if (radius.size() < 2 || !(coordinate_scale > 0.0) || !std::isfinite(coordinate_scale)) {
    throw std::invalid_argument("LogBesselKernel: invalid radial grid or scale");
  }
  for (const auto& i : indices(radius)) {
    if (!std::isfinite(radius[i]) || !(radius[i] > 0.0) || (i > 0 && !(radius[i] > radius[i - 1]))) {
      throw std::invalid_argument("LogBesselKernel: radius must increase strictly");
    }
  }
  MMatrix<double> kernel(momentum.size(), radius.size(), 0.0);
  for (const auto& k : indices(momentum)) {
    const double q = momentum[k];
    if (!std::isfinite(q) || q < 0.0) { throw std::invalid_argument("LogBesselKernel: invalid momentum"); }
    for (std::size_t i = 1; i < radius.size(); ++i) {
      const double lo = radius[i - 1] / coordinate_scale, hi = radius[i] / coordinate_scale;
      const double log_ratio = std::log(lo / hi);
      double difference = 0.0;
      if (q * hi < 1.0) {
        // Series after dividing by q^2 remains regular at zero momentum
        double term = hi * hi / 4.0;
        for (unsigned int m = 1; ; ++m) {
          const double value = -term * std::expm1(2.0 * m * log_ratio);
          difference += value;
          if (std::abs(value) <= std::numeric_limits<double>::epsilon() * std::abs(difference)) { break; }
          term *= -math::pow2(q * hi / (2.0 * (m + 1)));
        }
      } else {
        difference = (std::cyl_bessel_j(0.0, q * lo) - std::cyl_bessel_j(0.0, q * hi)) / (q * q);
      }
      const double weight = -2.0 * math::PI * difference / log_ratio;
      kernel[k][i - 1] += weight;
      kernel[k][i] -= weight;
    }
  }
  return kernel;
}

// Compute a finite radial order-zero Fourier-Bessel transform
// F(q) = 2pi/s^2 int r dr J_0(qr/s) f(r)
double RadialBesselTransform0(const std::vector<double> &node,
                              const std::vector<double> &weight,
                              const std::vector<double> &field,
                              const double momentum,
                              const double coordinate_scale) {
  if (node.empty() || node.size() != weight.size() ||
      node.size() != field.size() || !std::isfinite(momentum) ||
      !std::isfinite(coordinate_scale) || !(coordinate_scale > 0.0)) {
    throw std::invalid_argument("RadialBesselTransform0: invalid input");
  }
  double integral = 0.0;
  double correction = 0.0;
  for (const auto &i : indices(node)) {
    if (!std::isfinite(node[i]) || node[i] < 0.0 || !std::isfinite(weight[i]) ||
        !std::isfinite(field[i])) {
      throw std::invalid_argument(
          "RadialBesselTransform0: non-finite radial field");
    }
    const double term =
        weight[i] * node[i] *
        std::cyl_bessel_j(0.0, momentum * node[i] / coordinate_scale) *
        field[i];
    const double corrected = term - correction;
    const double total = integral + corrected;
    correction = (total - integral) - corrected;
    integral = total;
  }
  return 2.0 * PI * integral / (coordinate_scale * coordinate_scale);
}

// Construct a linear Gauss-Legendre rule including the radial measure r dr
// measure_weight_i = w_i r_i
PolarMeasureRule GaussLegendrePolarMeasure(const unsigned int node_count,
                                           const double minimum,
                                           const double maximum) {
  if (minimum < 0.0) {
    throw std::invalid_argument(
        "GaussLegendrePolarMeasure: minimum radius must be non-negative");
  }
  const auto rule = GaussLegendreRule(node_count, minimum, maximum);
  PolarMeasureRule out{rule.first, rule.second};
  for (const auto &i : indices(out.measure_weight)) {
    out.measure_weight[i] *= out.node[i];
  }
  ValidatePolarRule(out, "GaussLegendrePolarMeasure");
  return out;
}

// Initialize fixed momentum, impact parameter, and angular quadrature grids
// f_m(k) = N_phi^-1 sum_j f(k,phi_j) exp(-i m phi_j)
MPolarFourier::MPolarFourier(PolarMeasureRule momentum_rule,
                             PolarMeasureRule impact_rule,
                             std::vector<double> azimuth_node,
                             const std::size_t max_harmonic)
    : momentum_rule_(std::move(momentum_rule)),
      impact_rule_(std::move(impact_rule)),
      azimuth_node_(std::move(azimuth_node)), max_harmonic_(max_harmonic) {
  ValidatePolarRule(momentum_rule_, "MPolarFourier momentum rule");
  ValidatePolarRule(impact_rule_, "MPolarFourier impact rule");
  ValidateAzimuthNodes(azimuth_node_);
  if (max_harmonic_ >
          static_cast<std::size_t>(std::numeric_limits<int>::max()) ||
      2 * max_harmonic_ + 1 > azimuth_node_.size()) {
    throw std::invalid_argument(
        "MPolarFourier: too many harmonics for the angular quadrature");
  }

  angular_kernel_ = MMatrix<std::complex<double>>(2 * max_harmonic_ + 1,
                                                  azimuth_node_.size(), 0.0);
  const double angular_weight = 1.0 / azimuth_node_.size();
  for (int harmonic = -static_cast<int>(max_harmonic_);
       harmonic <= static_cast<int>(max_harmonic_); ++harmonic) {
    const std::size_t column = HarmonicColumn(harmonic, max_harmonic_);
    for (const auto &j : indices(azimuth_node_)) {
      angular_kernel_(column, j) =
          angular_weight *
          std::exp(std::complex<double>(0.0, -harmonic * azimuth_node_[j]));
    }
  }

  forward_kernel_.reserve(max_harmonic_ + 1);
  for (std::size_t order = 0; order <= max_harmonic_; ++order) {
    MMatrix<double> kernel(impact_rule_.node.size(), momentum_rule_.node.size(),
                           0.0);
    for (const auto &ib : indices(impact_rule_.node)) {
      for (const auto &ik : indices(momentum_rule_.node)) {
        kernel(ib, ik) = momentum_rule_.measure_weight[ik] *
                         IntegerBessel(order, momentum_rule_.node[ik] *
                                                  impact_rule_.node[ib]);
      }
    }
    forward_kernel_.push_back(std::move(kernel));
  }
}

// Compute the momentum radial quadrature rule
const PolarMeasureRule &MPolarFourier::MomentumRule() const {
  return momentum_rule_;
}

// Compute the impact parameter radial quadrature rule
const PolarMeasureRule &MPolarFourier::ImpactRule() const {
  return impact_rule_;
}

// Compute the uniform angular quadrature nodes
const std::vector<double> &MPolarFourier::AzimuthNodes() const {
  return azimuth_node_;
}

// Compute the largest retained harmonic of each input field
std::size_t MPolarFourier::MaxInputHarmonic() const { return max_harmonic_; }

// Apply F(b)=1/(2 pi) int d2k exp(+i k.b) f(k)
PolarHarmonicField MPolarFourier::Forward(
    const MMatrix<std::complex<double>> &momentum_field) const {
  if (momentum_field.size_row() != momentum_rule_.node.size() ||
      momentum_field.size_col() != azimuth_node_.size()) {
    throw std::invalid_argument(
        "MPolarFourier::Forward: momentum field dimensions disagree");
  }
  if (!gra::AllFinite(momentum_field.Elements())) {
    throw gra::AmplitudeFailure(
        "MPolarFourier::Forward: non-finite momentum field");
  }

  MMatrix<std::complex<double>> momentum_harmonic(momentum_rule_.node.size(),
                                                  2 * max_harmonic_ + 1, 0.0);
  for (const auto &ik : indices(momentum_rule_.node)) {
    for (std::size_t column = 0; column < angular_kernel_.size_row();
         ++column) {
      std::complex<double> value = 0.0;
      for (const auto &j : indices(azimuth_node_)) {
        const auto sample = momentum_field(ik, j);
        value += angular_kernel_(column, j) * sample;
      }
      momentum_harmonic(ik, column) = value;
    }
  }

  auto out = ZeroHarmonicField(impact_rule_.node.size(), max_harmonic_);
  for (int harmonic = -static_cast<int>(max_harmonic_);
       harmonic <= static_cast<int>(max_harmonic_); ++harmonic) {
    const std::size_t order = static_cast<std::size_t>(std::abs(harmonic));
    const std::size_t column = HarmonicColumn(harmonic, max_harmonic_);
    const std::complex<double> phase = PositiveImaginaryPhase(order);
    for (const auto &ib : indices(impact_rule_.node)) {
      std::complex<double> value = 0.0;
      for (const auto &ik : indices(momentum_rule_.node)) {
        value += forward_kernel_[order](ib, ik) * momentum_harmonic(ik, column);
      }
      out.coefficient(ib, column) = phase * value;
    }
  }
  return out;
}

// Multiply two impact parameter fields by discrete angular convolution
// (FG)_m(b) = sum_n F_n(b) G_{m-n}(b)
PolarHarmonicField
MPolarFourier::Multiply(const PolarHarmonicField &first,
                        const PolarHarmonicField &second) const {
  ValidateHarmonicField(first, impact_rule_.node.size());
  ValidateHarmonicField(second, impact_rule_.node.size());
  if (first.max_harmonic >
      std::numeric_limits<std::size_t>::max() - second.max_harmonic) {
    throw std::overflow_error("MPolarFourier::Multiply: harmonic overflow");
  }
  const std::size_t output_max = first.max_harmonic + second.max_harmonic;
  if (output_max > static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::overflow_error(
        "MPolarFourier::Multiply: harmonic order exceeds integer range");
  }
  auto out = ZeroHarmonicField(impact_rule_.node.size(), output_max);
  for (const auto &ib : indices(impact_rule_.node)) {
    for (int first_harmonic = -static_cast<int>(first.max_harmonic);
         first_harmonic <= static_cast<int>(first.max_harmonic);
         ++first_harmonic) {
      const std::size_t first_column = static_cast<std::size_t>(
          first_harmonic + static_cast<int>(first.max_harmonic));
      for (int second_harmonic = -static_cast<int>(second.max_harmonic);
           second_harmonic <= static_cast<int>(second.max_harmonic);
           ++second_harmonic) {
        const std::size_t second_column = static_cast<std::size_t>(
            second_harmonic + static_cast<int>(second.max_harmonic));
        const std::size_t output_column = static_cast<std::size_t>(
            first_harmonic + second_harmonic + static_cast<int>(output_max));
        out.coefficient(ib, output_column) +=
            first.coefficient(ib, first_column) *
            second.coefficient(ib, second_column);
      }
    }
  }
  return out;
}

// Multiply three impact parameter fields by discrete angular convolution
// (FGH)_m = sum_{n,l} F_n G_l H_{m-n-l}
PolarHarmonicField
MPolarFourier::Multiply(const PolarHarmonicField &first,
                        const PolarHarmonicField &second,
                        const PolarHarmonicField &third) const {
  return Multiply(Multiply(first, second), third);
}

// Prepare the inverse kernel once for a fixed output momentum
// K_{m,b} = w_b J_|m|(qb) (-i)^|m| exp(i m phi)
PolarInverseKernel
MPolarFourier::PrepareInverseKernel(const double q, const double azimuth,
                                    const std::size_t max_harmonic) const {
  if (!std::isfinite(q) || q < 0.0 || !std::isfinite(azimuth)) {
    throw std::invalid_argument(
        "MPolarFourier::PrepareInverseKernel: invalid output momentum");
  }
  if (max_harmonic >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(
        "MPolarFourier::PrepareInverseKernel: harmonic order exceeds integer "
        "range");
  }

  PolarInverseKernel kernel;
  kernel.max_harmonic = max_harmonic;
  kernel.radial_weight =
      MMatrix<double>(max_harmonic + 1, impact_rule_.node.size(), 0.0);
  kernel.angular_phase.resize(2 * max_harmonic + 1);
  for (std::size_t order = 0; order <= max_harmonic; ++order) {
    for (const auto &ib : indices(impact_rule_.node)) {
      kernel.radial_weight(order, ib) =
          impact_rule_.measure_weight[ib] *
          IntegerBessel(order, q * impact_rule_.node[ib]);
    }
  }
  for (int harmonic = -static_cast<int>(max_harmonic);
       harmonic <= static_cast<int>(max_harmonic); ++harmonic) {
    const std::size_t order = static_cast<std::size_t>(std::abs(harmonic));
    const std::size_t column = HarmonicColumn(harmonic, max_harmonic);
    kernel.angular_phase[column] =
        std::conj(PositiveImaginaryPhase(order)) *
        std::exp(std::complex<double>(0.0, harmonic * azimuth));
  }
  return kernel;
}

// Evaluate one impact parameter field with a prepared inverse kernel
std::complex<double>
MPolarFourier::InverseEvaluate(const PolarHarmonicField &impact_field,
                               const PolarInverseKernel &kernel) const {
  ValidateHarmonicField(impact_field, impact_rule_.node.size());
  ValidateInverseKernel(kernel, impact_rule_.node.size());
  return ApplyInverseKernel(impact_field, kernel);
}

// Evaluate several impact parameter fields with one prepared inverse kernel
std::vector<std::complex<double>> MPolarFourier::InverseEvaluateBatch(
    const std::vector<PolarHarmonicField> &impact_field,
    const PolarInverseKernel &kernel) const {
  ValidateInverseKernel(kernel, impact_rule_.node.size());
  std::vector<std::complex<double>> output(impact_field.size());
  for (const auto &i : indices(impact_field)) {
    ValidateHarmonicField(impact_field[i], impact_rule_.node.size());
    output[i] = ApplyInverseKernel(impact_field[i], kernel);
  }
  return output;
}

// Evaluate f(q)=1/(2 pi) int d2b exp(-i q.b) F(b) at arbitrary q
std::complex<double>
MPolarFourier::InverseEvaluate(const PolarHarmonicField &impact_field,
                               const double q, const double azimuth) const {
  ValidateHarmonicField(impact_field, impact_rule_.node.size());
  const auto kernel =
      PrepareInverseKernel(q, azimuth, impact_field.max_harmonic);
  return ApplyInverseKernel(impact_field, kernel);
}

// Evaluate the direct two-field transverse convolution at arbitrary q
// int d2k f(k)g(q-k) = 2pi F^{-1}[F(b)G(b)]
std::complex<double>
MPolarFourier::ConvolutionEvaluate(const PolarHarmonicField &first,
                                   const PolarHarmonicField &second,
                                   const double q, const double azimuth) const {
  return 2.0 * PI * InverseEvaluate(Multiply(first, second), q, azimuth);
}

// Evaluate the direct three-field transverse convolution at arbitrary q
// int d2k d2l f(k)g(l)h(q-k-l) = (2pi)^2 F^{-1}[FGH]
std::complex<double>
MPolarFourier::ConvolutionEvaluate(const PolarHarmonicField &first,
                                   const PolarHarmonicField &second,
                                   const PolarHarmonicField &third,
                                   const double q, const double azimuth) const {
  return 4.0 * PIPI *
         InverseEvaluate(Multiply(first, second, third), q, azimuth);
}

// Build one transform from an existing polar momentum rule
MPolarFourier BuildPolarFourierTransform(PolarMeasureRule momentum_rule,
                                         const unsigned int azimuth_count,
                                         const unsigned int impact_count,
                                         const double impact_max,
                                         const std::size_t max_harmonic) {
  if (azimuth_count == 0 || impact_count == 0 || !std::isfinite(impact_max) ||
      impact_max <= 0.0) {
    throw std::invalid_argument(
        "BuildPolarFourierTransform: invalid transform dimensions");
  }
  const auto impact = GaussLegendrePolarMeasure(impact_count, 0.0, impact_max);
  const auto azimuth = PeriodicTrapzRule(azimuth_count, 0.0, 2.0 * PI);
  return {std::move(momentum_rule), impact, azimuth.first, max_harmonic};
}

} // namespace gra::math
