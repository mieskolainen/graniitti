// Unit tests for non-radial polar Fourier-Bessel convolutions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Math/MPolarQuadrature.h"
#include "Graniitti/Tech/MAux.h"
#include "catch.hpp"

using gra::aux::indices;

namespace {

// Store one shifted two-dimensional Gaussian field
struct ShiftedGaussian {
  double width = 1.0;
  double x = 0.0;
  double y = 0.0;
  std::complex<double> coefficient = 1.0;
};

// Sample one shifted Gaussian on the transform momentum grid
gra::MMatrix<std::complex<double>>
SampleGaussian(const gra::math::MPolarFourier &transform,
               const ShiftedGaussian &gaussian) {
  const auto &radial = transform.MomentumRule().node;
  const auto &azimuth = transform.AzimuthNodes();
  gra::MMatrix<std::complex<double>> out(radial.size(), azimuth.size(), 0.0);
  for (const auto &i : indices(radial)) {
    for (const auto &j : indices(azimuth)) {
      const double x = radial[i] * std::cos(azimuth[j]);
      const double y = radial[i] * std::sin(azimuth[j]);
      const double dx = x - gaussian.x;
      const double dy = y - gaussian.y;
      out(i, j) = gaussian.coefficient *
                  std::exp(-gaussian.width * (dx * dx + dy * dy));
    }
  }
  return out;
}

// Evaluate the analytic convolution of shifted two-dimensional Gaussians
template <std::size_t N>
std::complex<double>
GaussianConvolution(const std::array<ShiftedGaussian, N> &gaussian,
                    const double qx, const double qy) {
  double inverse_width_sum = 0.0;
  double center_x = 0.0;
  double center_y = 0.0;
  double width_product = 1.0;
  std::complex<double> coefficient = 1.0;
  for (const auto &field : gaussian) {
    inverse_width_sum += 1.0 / field.width;
    center_x += field.x;
    center_y += field.y;
    width_product *= field.width;
    coefficient *= field.coefficient;
  }
  const double effective_width = 1.0 / inverse_width_sum;
  const double dx = qx - center_x;
  const double dy = qy - center_y;
  const double normalization =
      std::pow(gra::math::PI, static_cast<double>(N - 1)) * effective_width /
      width_product;
  return coefficient * normalization *
         std::exp(-effective_width * (dx * dx + dy * dy));
}

// Construct one transform resolving non-collinear shifted Gaussian fields
gra::math::MPolarFourier GaussianTransform() {
  const auto momentum = gra::math::GaussLegendrePolarMeasure(56, 0.0, 6.0);
  return gra::math::BuildPolarFourierTransform(momentum, 64, 72, 11.0, 20);
}

} // namespace

TEST_CASE("Polar Fourier-Bessel transform reproduces shifted Gaussian "
          "convolutions",
          "[MPolarFourier][convolution]") {
  const auto transform = GaussianTransform();
  const std::array<ShiftedGaussian, 3> gaussian = {
      ShiftedGaussian{1.15, 0.42, -0.18, {1.10, 0.25}},
      ShiftedGaussian{0.83, -0.27, 0.36, {0.72, -0.31}},
      ShiftedGaussian{1.37, 0.16, 0.29, {0.91, 0.18}}};
  std::array<gra::math::PolarHarmonicField, 3> field;
  for (const auto &i : indices(gaussian)) {
    field[i] = transform.Forward(SampleGaussian(transform, gaussian[i]));
  }

  SECTION("Two-field normalization includes one direct convolution measure") {
    const double qx = 0.31;
    const double qy = -0.22;
    const double q = std::hypot(qx, qy);
    const double phi = std::atan2(qy, qx);
    const std::array<ShiftedGaussian, 2> pair = {gaussian[0], gaussian[1]};
    const auto expected = GaussianConvolution(pair, qx, qy);
    const auto value =
        transform.ConvolutionEvaluate(field[0], field[1], q, phi);
    REQUIRE(value.real() == Approx(expected.real()).epsilon(2.0e-7));
    REQUIRE(value.imag() == Approx(expected.imag()).epsilon(2.0e-7));
  }

  SECTION(
      "Three-field normalization includes two direct convolution measures") {
    const double qx = -0.24;
    const double qy = 0.39;
    const double q = std::hypot(qx, qy);
    const double phi = std::atan2(qy, qx);
    const auto expected = GaussianConvolution(gaussian, qx, qy);
    const auto value =
        transform.ConvolutionEvaluate(field[0], field[1], field[2], q, phi);
    REQUIRE(value.real() == Approx(expected.real()).epsilon(3.0e-7));
    REQUIRE(value.imag() == Approx(expected.imag()).epsilon(3.0e-7));
  }

  SECTION("One prepared inverse kernel evaluates a field batch") {
    const double qx = 0.28;
    const double qy = -0.34;
    const double q = std::hypot(qx, qy);
    const double phi = std::atan2(qy, qx);
    const std::array<ShiftedGaussian, 2> first_pair = {gaussian[0],
                                                       gaussian[1]};
    const std::array<ShiftedGaussian, 2> second_pair = {gaussian[1],
                                                        gaussian[2]};
    const std::vector<gra::math::PolarHarmonicField> product = {
        transform.Multiply(field[0], field[1]),
        transform.Multiply(field[1], field[2])};
    const auto kernel =
        transform.PrepareInverseKernel(q, phi, product.front().max_harmonic);
    const auto value = transform.InverseEvaluateBatch(product, kernel);
    const std::array<std::complex<double>, 2> expected = {
        GaussianConvolution(first_pair, qx, qy) / (2.0 * gra::math::PI),
        GaussianConvolution(second_pair, qx, qy) / (2.0 * gra::math::PI)};

    REQUIRE(value.size() == expected.size());
    for (const auto &i : indices(expected)) {
      REQUIRE(value[i].real() == Approx(expected[i].real()).epsilon(2.0e-7));
      REQUIRE(value[i].imag() == Approx(expected[i].imag()).epsilon(2.0e-7));
      const auto scalar = transform.InverseEvaluate(product[i], q, phi);
      REQUIRE(value[i].real() == Approx(scalar.real()).epsilon(1.0e-13));
      REQUIRE(value[i].imag() == Approx(scalar.imag()).epsilon(1.0e-13));
    }
  }
}

TEST_CASE("Polar Fourier-Bessel quadrature validates angular resolution",
          "[MPolarFourier][validation]") {
  const auto momentum = gra::math::GaussLegendrePolarMeasure(4, 0.0, 2.0);
  const auto impact = gra::math::GaussLegendrePolarMeasure(5, 0.0, 3.0);
  const auto azimuth =
      gra::math::PeriodicTrapzRule(8, 0.0, 2.0 * gra::math::PI);
  REQUIRE_THROWS_AS(
      gra::math::MPolarFourier(momentum, impact, azimuth.first, 4),
      std::invalid_argument);
}

// Check the radial measure at the origin for closed Newton-Cotes rules
TEST_CASE("Polar Fourier accepts zero measure at the radial origin", "[MPolarFourier][quadrature]") {
  for (const std::string integrator : {"1/3", "3/8", "Boole"}) {
    const gra::math::PolarParam       param{integrator, "Trap", gra::math::RadialMap::Linear, 0.0, 1.0, 12, 8};
    const auto                        radial = gra::math::PolarRadialRule(param);
    const gra::math::PolarMeasureRule momentum{radial.node, radial.measure};
    REQUIRE_NOTHROW(gra::math::BuildPolarFourierTransform(momentum, 8, 8, 5.0, 2));
    const gra::math::PolarMeasureRule        impact{{0.0, 1.0}, {0.0, 0.5}};
    const auto                               azimuth = gra::math::PeriodicTrapzRule(8, 0.0, 2.0 * gra::math::PI);
    const gra::math::MPolarFourier           transform(momentum, impact, azimuth.first, 2);
    const gra::MMatrix<std::complex<double>> constant(radial.node.size(), 8, 1.0);
    const auto                               field = transform.Forward(constant);
    REQUIRE(field.coefficient(0, 2).real() == Approx(0.5).epsilon(0.0).margin(1.0e-14));
    REQUIRE(field.coefficient(0, 2).imag() == Approx(0.0).margin(1.0e-14));
    auto invalid                   = momentum;
    invalid.measure_weight.front() = -1.0;
    REQUIRE_THROWS_AS(gra::math::BuildPolarFourierTransform(invalid, 8, 8, 5.0, 2), std::invalid_argument);
    std::fill(invalid.measure_weight.begin(), invalid.measure_weight.end(), 0.0);
    REQUIRE_THROWS_AS(gra::math::BuildPolarFourierTransform(invalid, 8, 8, 5.0, 2), std::invalid_argument);
  }
}

// Resolve a smooth profile across many decades without truncating its long tail
TEST_CASE("Logarithmic Hankel kernels reproduce Gaussian amplitudes", "[fourier][EMD]") {
  const std::vector<double> momentum = {0.0, 1.0e-12, 0.03, 0.2, 0.7, 2.0};
  const double scale = 1.7;
  for (const double width : {0.01, 1.0, 100.0}) {
    std::vector<double> radius(4097), field(radius.size());
    const double minimum = 1.0e-6 / std::sqrt(width), maximum = 20.0 / std::sqrt(width);
    for (const auto& i : indices(radius)) {
      radius[i] = minimum * std::exp(std::log(maximum / minimum) * i / (radius.size() - 1));
      field[i] = std::exp(-width * radius[i] * radius[i]);
    }
    const auto result = gra::math::LogBesselKernel(radius, momentum, scale) * field;
    for (const auto& i : indices(momentum)) {
      const double norm = gra::math::PI / (width * scale * scale);
      const double expected = norm * std::exp(-momentum[i] * momentum[i] / (4.0 * width * scale * scale));
      CHECK(result[i] == Approx(expected).margin(2.0e-5 * norm));
    }
  }
}
