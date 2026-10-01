// Unit tests for MMatrix, MTensor, M4Vec class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

//#define CATCH_CONFIG_MAIN

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <span>
#include <sstream>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MDistribution.h"
#include "Graniitti/Math/MFFT.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Math/MTensor.h"
#include "Graniitti/Math/MTransport.h"
#include "Graniitti/Tech/MAux.h"
#include "catch.hpp"

using namespace gra;
using gra::aux::indices;

constexpr double EPS = 1e-10;

// Compare a single tensor-axis contraction to the full Kronecker operator
TEST_CASE("MMatrix transforms one column tensor factor", "[MMatrix][LinearAlgebra]") {
  const MMatrix<std::complex<double>> op = {{{0.3, 0.2}, {-0.4, 0.7}},
                                           {{0.1, -0.5}, {0.8, 0.3}},
                                           {{-0.2, 0.6}, {0.4, -0.1}}};
  MMatrix<std::complex<double>> amplitude(3, 12, 0.0);
  for (std::size_t i = 0; i < amplitude.size_row(); ++i) {
    for (std::size_t j = 0; j < amplitude.size_col(); ++j) {
      amplitude[i][j] = {0.1 * (i + j), 0.2 * static_cast<double>(j) - 0.3 * i};
    }
  }
  const auto identity2 = MMatrix<std::complex<double>>::IdentityMatrix(2);
  const auto identity3 = MMatrix<std::complex<double>>::IdentityMatrix(3);
  const auto explicit_op = identity2.Kronecker(op).Kronecker(identity3);
  const auto actual = amplitude.TransformColumnAxis(op, 2, 3);
  CHECK((actual - amplitude * explicit_op.Transpose()).FrobNorm() < 1e-12);
  REQUIRE_THROWS_AS(amplitude.TransformColumnAxis(op, 3, 3), std::invalid_argument);
  REQUIRE_THROWS_AS(amplitude.TransformColumnAxis(op, 0, 3), std::invalid_argument);
}

// Check ordered tensor traces and transposes on complex nonsymmetric factors
TEST_CASE("MMatrix partial tensor operations preserve factor order", "[MMatrix][LinearAlgebra]") {
  using Complex = std::complex<double>;
  using gra::TensorFactor;
  const MMatrix<Complex> a = {{{0.2, 0.3}, {-0.4, 0.7}}, {{0.8, -0.1}, {0.5, -0.2}}};
  const MMatrix<Complex> b = {{{0.1, 0.2}, {0.3, -0.4}, {-0.6, 0.5}},
                             {{0.7, -0.2}, {0.4, 0.6}, {0.8, -0.3}},
                             {{-0.2, -0.7}, {0.9, 0.1}, {0.5, -0.4}}};
  const auto product = a.Kronecker(b);
  CHECK((product.PartialTrace(2, 3, TensorFactor::Second) - a * b.Trace()).FrobNorm() < 1e-14);
  CHECK((product.PartialTrace(2, 3, TensorFactor::First) - b * a.Trace()).FrobNorm() < 1e-14);
  CHECK((product.PartialTranspose(2, 3, TensorFactor::First) - a.Transpose().Kronecker(b)).FrobNorm() < 1e-14);
  CHECK((product.PartialTranspose(2, 3, TensorFactor::Second) - a.Kronecker(b.Transpose())).FrobNorm() < 1e-14);
  for (const auto factor : {TensorFactor::First, TensorFactor::Second}) {
    CHECK((product.PartialTranspose(2, 3, factor).PartialTranspose(2, 3, factor) - product).FrobNorm() < 1e-14);
    REQUIRE_THROWS_AS(product.PartialTrace(0, 3, factor), std::invalid_argument);
    REQUIRE_THROWS_AS(product.PartialTranspose(2, 2, factor), std::invalid_argument);
    REQUIRE_THROWS_AS(product.PartialTranspose(std::numeric_limits<std::size_t>::max(), 2, factor), std::invalid_argument);
    REQUIRE_THROWS_AS(MMatrix<Complex>(2, 3).PartialTrace(1, 2, factor), std::invalid_argument);
  }
  CHECK((product.PartialTranspose(2, 3, TensorFactor::First).PartialTranspose(2, 3, TensorFactor::Second) -
         product.Transpose()).FrobNorm() < 1e-14);
}

TEST_CASE("MFloat exact predicates preserve IEEE states", "[MFloat]") {
  const double positive_zero = 0.0;
  const double negative_zero = -0.0;
  const double subnormal     = std::numeric_limits<double>::denorm_min();
  const double nan           = std::numeric_limits<double>::quiet_NaN();
  const double infinity      = std::numeric_limits<double>::infinity();

  REQUIRE(gra::math::IsZero(positive_zero));
  REQUIRE(gra::math::IsZero(negative_zero));
  REQUIRE(gra::math::IsExactEqual(positive_zero, negative_zero));
  REQUIRE_FALSE(gra::math::IsZero(subnormal));
  REQUIRE_FALSE(gra::math::IsExactEqual(nan, nan));
  REQUIRE(gra::math::IsExactEqual(infinity, infinity));
}

TEST_CASE("WrapAngle preserves the signed angular convention", "[MMath]") {
  REQUIRE(gra::math::WrapAngle(0.0) == Approx(0.0).margin(1.0e-15));
  REQUIRE(gra::math::WrapAngle(gra::math::PI) == Approx(gra::math::PI).margin(1.0e-15));
  REQUIRE(gra::math::WrapAngle(-gra::math::PI) == Approx(gra::math::PI).margin(1.0e-15));
  REQUIRE(gra::math::WrapAngle(3.0 * gra::math::PI) == Approx(gra::math::PI).margin(1.0e-15));
  REQUIRE(gra::math::WrapAngle(-1.5 * gra::math::PI) == Approx(0.5 * gra::math::PI).margin(1.0e-15));
  REQUIRE(gra::math::sign(-7) == -1);
  REQUIRE(gra::math::sign(0) == 0);
  REQUIRE(gra::math::sign(4) == 1);
}

TEST_CASE("MFFT obeys the discrete Fourier transform and inverse", "[MFFT]") {
  using Complex                  = std::complex<double>;
  std::valarray<Complex> impulse = {Complex(1.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0), Complex(0.0, 0.0)};
  gra::MFFT::fft(impulse);
  for (const Complex value : impulse) {
    REQUIRE(value.real() == Approx(1.0).margin(1.0e-15));
    REQUIRE(value.imag() == Approx(0.0).margin(1.0e-15));
  }

  std::valarray<Complex> signal = {Complex(0.3, -0.2), Complex(-1.1, 0.7), Complex(2.0, 0.4),   Complex(0.5, -0.9),
                                   Complex(-0.6, 1.3), Complex(0.8, 0.2),  Complex(-0.4, -0.7), Complex(1.2, 0.1)};
  const std::valarray<Complex> reference = signal;
  gra::MFFT::fft(signal);
  gra::MFFT::ifft(signal);
  for (std::size_t i = 0; i < signal.size(); ++i) {
    REQUIRE(signal[i].real() == Approx(reference[i].real()).margin(2.0e-14));
    REQUIRE(signal[i].imag() == Approx(reference[i].imag()).margin(2.0e-14));
  }

  std::valarray<Complex>       invalid(Complex(1.0, 2.0), 3);
  const std::valarray<Complex> unchanged = invalid;
  REQUIRE_THROWS_AS(gra::MFFT::ifft(invalid), std::invalid_argument);
  for (std::size_t i = 0; i < invalid.size(); ++i) {
    REQUIRE(gra::math::IsExactEqual(invalid[i].real(), unchanged[i].real()));
    REQUIRE(gra::math::IsExactEqual(invalid[i].imag(), unchanged[i].imag()));
  }
}

TEST_CASE("MDistribution quantiles invert normalized probability laws", "[MDistribution]") {
  const gra::math::GammaNumerics numerics{1.0e-13, 10000, 1.0e12};

  REQUIRE(gra::math::NormalQuantile(0.5) == Approx(0.0).margin(1.0e-15));
  REQUIRE(gra::math::NormalQuantile(0.975) == Approx(1.95996398454005).margin(5.0e-9));
  for (const double probability : {1.0e-5, 0.1, 0.9, 1.0 - 1.0e-5}) {
    REQUIRE(gra::math::NormalQuantile(1.0 - probability) ==
            Approx(-gra::math::NormalQuantile(probability)).margin(2.0e-10));
  }

  for (const double shape : {0.2, 1.0, 3.7, 25.0}) {
    for (const double probability : {1.0e-5, 0.1, 0.5, 0.9, 1.0 - 1.0e-5}) {
      const double quantile = gra::math::GammaQuantile(probability, shape, numerics);
      REQUIRE(gra::math::GammaCDF(shape, quantile, numerics) == Approx(probability).margin(2.0e-12));
    }
  }

  const double value = 1.7;
  REQUIRE(gra::math::GammaCDF(1.0, value, numerics) == Approx(1.0 - std::exp(-value)).margin(2.0e-14));
  REQUIRE(gra::math::GammaCDF(2.0, value, numerics) == Approx(1.0 - std::exp(-value) * (1.0 + value)).margin(2.0e-14));
  REQUIRE(gra::math::GammaQuantile(0.73, 1.0, numerics) == Approx(-std::log(1.0 - 0.73)).margin(2.0e-12));

  const gra::math::GammaNumerics approximate{1.0e-13, 10000, 10.0};
  const double                   approximate_quantile = gra::math::GammaQuantile(0.9, 25.0, approximate);
  REQUIRE(gra::math::GammaCDF(25.0, approximate_quantile, numerics) == Approx(0.9).margin(2.0e-4));

  REQUIRE_THROWS_AS(gra::math::NormalQuantile(0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::GammaNumerics(0.0, 10000, 1.0e12), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::GammaNumerics(1.0e-13, 0, 1.0e12), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::GammaNumerics(1.0e-13, 10000, 0.0), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::GammaCDF(-1.0, 1.0, numerics), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::GammaQuantile(1.0, 1.0, numerics), std::invalid_argument);
}

// ---------------------------------------------------------
// MMatrix

TEST_CASE("MMatrix dagger obeys involution and the Frobenius trace identity", "[MMatrix]") {
  using Complex = std::complex<double>;
  const MMatrix<Complex> matrix{{Complex(0.4, -0.3), Complex(1.2, 0.7)}, {Complex(-0.8, 0.2), Complex(0.1, -1.1)}};

  const MMatrix<Complex> dagger        = matrix.Dagger();
  const MMatrix<Complex> restored      = dagger.Dagger();
  double                 element_norm2 = 0.0;
  for (std::size_t i = 0; i < matrix.size_row(); ++i) {
    for (std::size_t j = 0; j < matrix.size_col(); ++j) {
      element_norm2 += std::norm(matrix(i, j));
      REQUIRE(restored(i, j).real() == Approx(matrix(i, j).real()).margin(1e-12));
      REQUIRE(restored(i, j).imag() == Approx(matrix(i, j).imag()).margin(1e-12));
    }
  }

  const Complex trace = (matrix * dagger).Trace();
  REQUIRE(trace.real() == Approx(element_norm2).margin(1e-12));
  REQUIRE(trace.imag() == Approx(0.0).margin(1e-12));
}

TEST_CASE("MMatrix Basic Initialization", "[MMatrix]") {
  gra::MMatrix<double> mat(3, 3);

  REQUIRE(mat.size_row() == 3);
  REQUIRE(mat.size_col() == 3);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) { REQUIRE(gra::math::IsZero(mat(i, j))); }
  }
}

TEST_CASE("MMatrix Initialization with value", "[MMatrix]") {
  gra::MMatrix<double> mat(2, 3, 5.0);

  REQUIRE(mat.size_row() == 2);
  REQUIRE(mat.size_col() == 3);

  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 3; ++j) { REQUIRE(gra::math::IsExactEqual(mat(i, j), 5.0)); }
  }
}

TEST_CASE("MMatrix Identity matrix", "[MMatrix]") {
  gra::MMatrix<double> mat(3, 3, "eye");

  REQUIRE(mat.size_row() == 3);
  REQUIRE(mat.size_col() == 3);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      if (i == j) {
        REQUIRE(gra::math::IsExactEqual(mat(i, j), 1.0));
      } else {
        REQUIRE(gra::math::IsZero(mat(i, j)));
      }
    }
  }
}

TEST_CASE("MMatrix Minkowski matrix", "[MMatrix]") {
  gra::MMatrix<double> mat(4, 4, "minkowski");

  REQUIRE(mat.size_row() == 4);
  REQUIRE(mat.size_col() == 4);

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      if (i == j) {
        REQUIRE(gra::math::IsExactEqual(mat(i, j), i > 0 ? -1.0 : 1.0));
      } else {
        REQUIRE(gra::math::IsZero(mat(i, j)));
      }
    }
  }
}

TEST_CASE("MMatrix Addition", "[MMatrix]") {
  gra::MMatrix<double> mat1(2, 2, 3.0);
  gra::MMatrix<double> mat2(2, 2, 2.0);

  gra::MMatrix<double> result = mat1 + mat2;

  REQUIRE(result.size_row() == 2);
  REQUIRE(result.size_col() == 2);

  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) { REQUIRE(gra::math::IsExactEqual(result(i, j), 5.0)); }
  }
}

TEST_CASE("MMatrix Subtraction", "[MMatrix]") {
  gra::MMatrix<double> mat1(2, 2, 5.0);
  gra::MMatrix<double> mat2(2, 2, 3.0);

  gra::MMatrix<double> result = mat1 - mat2;

  REQUIRE(result.size_row() == 2);
  REQUIRE(result.size_col() == 2);

  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) { REQUIRE(gra::math::IsExactEqual(result(i, j), 2.0)); }
  }
}

TEST_CASE("MMatrix Scalar Multiplication", "[MMatrix]") {
  gra::MMatrix<double> mat(2, 2, 3.0);
  gra::MMatrix<double> result = mat * 2.0;

  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) { REQUIRE(gra::math::IsExactEqual(result(i, j), 6.0)); }
  }
}

TEST_CASE("MMatrix Scalar Division", "[MMatrix]") {
  gra::MMatrix<double> mat(2, 2, 6.0);
  gra::MMatrix<double> result = mat / 2.0;

  for (std::size_t i = 0; i < 2; ++i) {
    for (std::size_t j = 0; j < 2; ++j) { REQUIRE(gra::math::IsExactEqual(result(i, j), 3.0)); }
  }
}

TEST_CASE("MMatrix Transpose", "[MMatrix]") {
  gra::MMatrix mat        = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  gra::MMatrix transposed = mat.Transpose();
  REQUIRE(transposed.size_row() == 3);
  REQUIRE(transposed.size_col() == 2);
  REQUIRE(gra::math::IsExactEqual(transposed(0, 0), 1.0));
  REQUIRE(gra::math::IsExactEqual(transposed(1, 0), 2.0));
  REQUIRE(gra::math::IsExactEqual(transposed(2, 0), 3.0));
  REQUIRE(gra::math::IsExactEqual(transposed(0, 1), 4.0));
  REQUIRE(gra::math::IsExactEqual(transposed(1, 1), 5.0));
  REQUIRE(gra::math::IsExactEqual(transposed(2, 1), 6.0));
}

TEST_CASE("MMatrix Multiplication", "[MMatrix]") {
  gra::MMatrix mat1   = {{1, 2, 3}, {4, 5, 6}};
  gra::MMatrix mat2   = {{7, 8}, {9, 10}, {11, 12}};
  gra::MMatrix result = mat1 * mat2;
  REQUIRE(result.size_row() == 2);
  REQUIRE(result.size_col() == 2);
  REQUIRE(result(0, 0) == 58);
  REQUIRE(result(0, 1) == 64);
  REQUIRE(result(1, 0) == 139);
  REQUIRE(result(1, 1) == 154);
}

TEST_CASE("MMatrix Vector Multiplication", "[MMatrix]") {
  gra::MMatrix        mat    = {{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  std::vector<double> vec    = {1.0, 2.0, 3.0};
  std::vector<double> result = mat * vec;
  REQUIRE(result.size() == 2);
  REQUIRE(gra::math::IsExactEqual(result[0], 14.0));
  REQUIRE(gra::math::IsExactEqual(result[1], 32.0));
}

// Verify mixed scalar matrix algebra and reusable container operations
TEST_CASE("MMatrix generic linear algebra operations", "[MMatrix]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<double> matrix{{1.0, 2.0}, {3.0, 4.0}};

  const std::vector<Complex> vector{{0.0, 1.0}, {1.0, 0.0}};
  const auto                 product = matrix * vector;
  REQUIRE(product[0] == Complex(2.0, 1.0));
  REQUIRE(product[1] == Complex(4.0, 3.0));

  const auto scaled = matrix * Complex(0.0, 2.0);
  REQUIRE(scaled(0, 1) == Complex(0.0, 4.0));
  REQUIRE(scaled(1, 1) == Complex(0.0, 8.0));

  const auto diagonal = matrix.RightDiagonalProduct(std::vector<Complex>{{1.0, 0.0}, {0.0, 1.0}});
  REQUIRE(diagonal(1, 0) == Complex(3.0, 0.0));
  REQUIRE(diagonal(1, 1) == Complex(0.0, 4.0));
  REQUIRE(matrix.RightDiagonalFrobNorm2(std::vector<Complex>{{1.0, 0.0}, {0.0, 1.0}}) ==
          Approx(diagonal.FrobNorm2()).margin(1e-14));
  const gra::MMatrix<double> scaled_norm_matrix{{1.0e200}};
  REQUIRE(scaled_norm_matrix.RightDiagonalFrobNorm2(std::vector<double>{1.0e-200}) == Approx(1.0));
  auto                      scaled_columns = matrix;
  const std::vector<double> column_scale{2.0, 3.0};
  scaled_columns.ScaleColumns(column_scale);
  REQUIRE(scaled_columns.Flatten() == matrix.RightDiagonalProduct(column_scale).Flatten());
  const auto left_diagonal = matrix.LeftDiagonalProduct(std::vector<double>{2.0, 3.0});
  REQUIRE(left_diagonal(0, 1) == Approx(4.0));
  REQUIRE(left_diagonal(1, 0) == Approx(9.0));
  REQUIRE(matrix.ColumnNorm2(1) == Approx(20.0));

  const auto                block = matrix.Submatrix(0, 1, 2, 1);
  const std::vector<double> expected_block{2.0, 4.0};
  REQUIRE(block.Flatten() == expected_block);
  const auto support = matrix.Transform([](double value) { return value > 2.0; });
  REQUIRE_FALSE(support(0, 1));
  REQUIRE(support(1, 0));

  const gra::MMatrix<Complex> hermitian{{Complex(1.0, 0.0), Complex(0.0, 1.0)},
                                        {Complex(0.0, -1.0), Complex(2.0, 0.0)}};
  REQUIRE(hermitian.IsHermitian(1e-12));

  std::vector<Complex> accumulated(2, 0.0);
  gra::AddScaled(accumulated, vector, Complex(0.0, 1.0));
  REQUIRE(accumulated[0] == Complex(-1.0, 0.0));
  REQUIRE(accumulated[1] == Complex(0.0, 1.0));

  std::vector<Complex> matrix_accumulated{{1.0, 0.0}, {0.0, -1.0}};
  matrix.MultiplyAdd(std::span<const Complex>(vector), std::span<Complex>(matrix_accumulated), 2.0);
  REQUIRE(matrix_accumulated[0] == Complex(5.0, 2.0));
  REQUIRE(matrix_accumulated[1] == Complex(8.0, 5.0));

  const auto columns = gra::MMatrix<double>::FromColumns(std::vector<std::vector<double>>{{1.0, 2.0}, {3.0, 4.0}});
  REQUIRE(columns(0, 1) == Approx(3.0));
  REQUIRE(columns(1, 0) == Approx(2.0));

  const gra::MMatrix<bool> mask{{true, false}, {false, true}};
  REQUIRE(matrix.MaskedSquaredNorm(mask) == Approx(17.0));

  const auto outer = gra::OuterProduct(std::vector<double>{1.0, 2.0}, std::vector<Complex>{{0.0, 1.0}});
  REQUIRE(outer(1, 0) == Complex(0.0, 2.0));

  gra::MMatrix<Complex> accumulated_outer(2, 1, Complex(1.0, 0.0));
  accumulated_outer.AddOuterProduct(std::vector<double>{1.0, 2.0}, std::vector<Complex>{{0.0, 1.0}}, Complex(2.0, 0.0));
  REQUIRE(accumulated_outer(0, 0) == Complex(1.0, 2.0));
  REQUIRE(accumulated_outer(1, 0) == Complex(1.0, 4.0));

  REQUIRE(hermitian.IsFinite());
  REQUIRE(gra::AllFinite(vector));
  auto nonfinite  = hermitian;
  nonfinite(0, 0) = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
  REQUIRE_FALSE(nonfinite.IsFinite());
  REQUIRE_FALSE(nonfinite.IsApprox(nonfinite));
  REQUIRE_FALSE(gra::AllFinite(nonfinite.Flatten()));

  const std::array<double, 3> first_vector  = {1.0, 2.0, 3.0};
  const std::array<double, 3> second_vector = {4.0, 5.0, 6.0};
  REQUIRE((gra::Add(first_vector, second_vector) == std::array<double, 3>{5.0, 7.0, 9.0}));
  REQUIRE((gra::Subtract(second_vector, first_vector) == std::array<double, 3>{3.0, 3.0, 3.0}));
  REQUIRE((gra::Negated(first_vector) == std::array<double, 3>{-1.0, -2.0, -3.0}));
  REQUIRE((gra::CrossProduct(first_vector, second_vector) == std::array<double, 3>{-3.0, 6.0, -3.0}));

  const auto interpolated = gra::math::LinearInterpolate(std::vector<double>{0.0, 2.0},
                                                         std::vector<gra::MMatrix<double>>{matrix, matrix * 3.0}, 0.5);
  REQUIRE(interpolated(1, 1) == Approx(6.0));
}

// Compare factorized and explicitly materialized Kronecker contractions
TEST_CASE("MMatrix factorized Kronecker contraction preserves row mappings", "[MMatrix][LinearAlgebra]") {
  using Complex = std::complex<double>;
  const std::vector<Complex>       upper_rows{{0.7, 0.2}, {-0.3, 0.5}};
  const std::vector<Complex>       upper_columns{{1.1, -0.4}, {0.2, 0.8}};
  const std::vector<Complex>       lower_rows{{-0.6, 0.1}, {0.9, -0.7}};
  const std::vector<Complex>       lower_columns{{0.4, 0.3}, {-1.2, 0.2}};
  const gra::MMatrix<Complex>      central{{Complex(0.0, 0.0), Complex(0.3, -0.2), Complex(0.0, 0.0)},
                                      {Complex(-0.7, 0.1), Complex(0.0, 0.0), Complex(0.5, 0.4)},
                                      {Complex(0.2, -0.9), Complex(0.0, 0.0), Complex(0.0, 0.0)},
                                      {Complex(0.0, 0.0), Complex(1.3, 0.6), Complex(-0.1, 0.8)}};
  const std::array<std::size_t, 4> destination{3, 1, 2, 0};

  const auto upper            = gra::OuterProduct(upper_rows, upper_columns);
  const auto lower            = gra::OuterProduct(lower_rows, lower_columns);
  const auto explicit_product = upper.KroneckerMultiply(lower, central, destination);
  const auto factorized_product =
      central.FactorizedKroneckerMultiply(upper_rows, upper_columns, lower_rows, lower_columns, destination);
  REQUIRE(factorized_product.IsApprox(explicit_product, 1e-13));

  const gra::MMatrix<Complex> empty_central(4, 0);
  const auto                  empty_product =
      empty_central.FactorizedKroneckerMultiply(upper_rows, upper_columns, lower_rows, lower_columns, destination);
  REQUIRE(empty_product.size_row() == 4);
  REQUIRE(empty_product.size_col() == 0);

  const auto mapped_rows = gra::MappedKroneckerProduct(upper_rows, lower_rows, destination);
  for (std::size_t upper_row = 0; upper_row < upper_rows.size(); ++upper_row) {
    for (std::size_t lower_row = 0; lower_row < lower_rows.size(); ++lower_row) {
      const std::size_t source = upper_row * lower_rows.size() + lower_row;
      REQUIRE(std::abs(mapped_rows[destination[source]] - upper_rows[upper_row] * lower_rows[lower_row]) <= 1.0e-14);
    }
  }
  REQUIRE_THROWS_AS(gra::MappedKroneckerProduct(upper_rows, lower_rows, std::array<std::size_t, 4>{0, 0, 2, 3}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::MappedKroneckerProduct(upper_rows, lower_rows, std::array<std::size_t, 3>{0, 1, 2}),
                    std::invalid_argument);
}

// Check exact diagonal structure and allocation free cached contraction
TEST_CASE("MMatrix exact diagonal operations", "[MMatrix][LinearAlgebra]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> matrix =
      gra::MMatrix<Complex>::DiagonalMatrix({Complex(2.0, 1.0), Complex(-1.0, 0.5), Complex(0.0, 0.0)});

  REQUIRE(matrix.IsDiagonal());
  const std::vector<Complex> diagonal = matrix.GetDiag();
  const std::vector<Complex> expected_diagonal{Complex(2.0, 1.0), Complex(-1.0, 0.5), Complex(0.0, 0.0)};
  REQUIRE(diagonal == expected_diagonal);

  auto perturbed  = matrix;
  perturbed(0, 1) = Complex(0.0, std::numeric_limits<double>::epsilon());
  REQUIRE_FALSE(perturbed.IsDiagonal());
  REQUIRE_FALSE(gra::MMatrix<double>(2, 3, 0.0).IsDiagonal());

  const std::vector<double>  finite_source{3.0, 4.0, 5.0};
  const std::vector<Complex> initial_target{Complex(1.0, 0.0), Complex(-2.0, 1.0), Complex(7.0, -1.0)};
  const Complex              scale(0.5, -1.0);
  auto                       dense_target = initial_target;
  matrix.MultiplyAdd(std::span<const double>(finite_source), std::span<Complex>(dense_target), scale);
  auto diagonal_target = initial_target;
  gra::DiagonalMultiplyAdd(std::span<const Complex>(diagonal), std::span<const double>(finite_source),
                           std::span<Complex>(diagonal_target), scale);
  REQUIRE(diagonal_target == dense_target);
  REQUIRE(diagonal_target[0] == Complex(7.0, -4.5));
  REQUIRE(diagonal_target[1] == Complex(-2.0, 6.0));
  REQUIRE(diagonal_target[2] == initial_target[2]);

  const std::vector<double> nonfinite_source{3.0, 4.0, std::numeric_limits<double>::quiet_NaN()};
  auto                      sparse_target = initial_target;
  gra::DiagonalMultiplyAdd(std::span<const Complex>(diagonal), std::span<const double>(nonfinite_source),
                           std::span<Complex>(sparse_target), scale);
  REQUIRE(sparse_target[2] == initial_target[2]);

  REQUIRE_THROWS_AS(
      gra::DiagonalMultiplyAdd(std::span<const Complex>(diagonal), std::span<const double>(finite_source).first(2),
                               std::span<Complex>(diagonal_target), scale),
      std::invalid_argument);
}

TEST_CASE("MMatrix Move Semantics", "[MMatrix]") {
  gra::MMatrix mat1 = {{3.0, 3.0}, {3.0, 3.0}};
  gra::MMatrix mat2 = std::move(mat1);
  REQUIRE(mat2.size_row() == 2);
  REQUIRE(mat2.size_col() == 2);
  REQUIRE(gra::math::IsExactEqual(mat2(0, 0), 3.0));

  mat1 = gra::MMatrix<double>{{1.0}};
  REQUIRE(mat1(0, 0) == Approx(1.0));
}

TEST_CASE("MMatrix Copy Semantics", "[MMatrix]") {
  gra::MMatrix<double> original{{1.0, 2.0}, {3.0, 4.0}};

  SECTION("Copy constructor deep-copies data") {
    gra::MMatrix<double> copy(original);
    original(0, 0) = 9.0;

    REQUIRE(copy.size_row() == 2);
    REQUIRE(copy.size_col() == 2);
    REQUIRE(copy(0, 0) == Approx(1.0));
    REQUIRE(copy(1, 1) == Approx(4.0));
  }

  SECTION("Copy assignment resizes and deep-copies data") {
    gra::MMatrix<double> assigned(1, 3, -1.0);
    assigned       = original;
    original(0, 1) = 8.0;

    REQUIRE(assigned.size_row() == 2);
    REQUIRE(assigned.size_col() == 2);
    REQUIRE(assigned(0, 1) == Approx(2.0));
    REQUIRE(assigned(1, 0) == Approx(3.0));
  }

  SECTION("Self copy-assignment leaves data unchanged") {
    original = original;

    REQUIRE(original.size_row() == 2);
    REQUIRE(original.size_col() == 2);
    REQUIRE(original(0, 0) == Approx(1.0));
    REQUIRE(original(1, 1) == Approx(4.0));
  }
}

TEST_CASE("MMatrix Trace", "[MMatrix]") {
  gra::MMatrix mat = {{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}};
  REQUIRE(gra::math::IsExactEqual(mat.Trace(), 6.0));
}

TEST_CASE("MMatrix flat copies and views retain row-major order", "[MMatrix]") {
  gra::MMatrix              mat  = {{1.0, 2.0}, {3.0, 4.0}};
  const std::vector<double> flat = mat.Flatten();
  REQUIRE(flat.size() == 4);
  REQUIRE(gra::math::IsExactEqual(flat[0], 1.0));
  REQUIRE(gra::math::IsExactEqual(flat[1], 2.0));
  REQUIRE(gra::math::IsExactEqual(flat[2], 3.0));
  REQUIRE(gra::math::IsExactEqual(flat[3], 4.0));

  auto elements = mat.Elements();
  REQUIRE(elements.size() == 4);
  elements[2]          = 8.0;
  const auto immutable = std::as_const(mat).Elements();
  REQUIRE(gra::math::IsExactEqual(immutable[2], 8.0));
}

TEST_CASE("MMatrix Error handling - out of bounds access", "[MMatrix]") {
  gra::MMatrix<double> mat(3, 3);
  REQUIRE_THROWS_AS(mat(3, 3), std::out_of_range);
}

TEST_CASE("MMatrix Error handling - addition with size mismatch", "[MMatrix]") {
  gra::MMatrix<double> mat1(3, 3, 1.0);
  gra::MMatrix<double> mat2(2, 3, 2.0);

  REQUIRE_THROWS_AS(mat1 + mat2, std::invalid_argument);
  REQUIRE_THROWS_AS(mat1 - mat2, std::invalid_argument);
}

TEST_CASE("MMatrix initializer-list validation", "[MMatrix]") {
  REQUIRE_THROWS_AS((gra::MMatrix<int>{{1, 2}, {3}}), std::invalid_argument);
  REQUIRE_THROWS_AS((gra::MMatrix<int>{std::vector<int>{1, 2}, std::vector<int>{3}}), std::invalid_argument);
}

TEST_CASE("MMatrix Complex Initialization and Conjugate Transpose", "[MMatrix]") {
  using Complex = std::complex<double>;
  gra::MMatrix<Complex> mat(2, 2);

  mat(0, 0) = Complex(1.0, 1.0);
  mat(0, 1) = Complex(2.0, -1.0);
  mat(1, 0) = Complex(3.0, 4.0);
  mat(1, 1) = Complex(4.0, -2.0);

  gra::MMatrix<Complex> conj_transposed = mat.ConjTranspose();

  REQUIRE(conj_transposed.size_row() == 2);
  REQUIRE(conj_transposed.size_col() == 2);

  // Check conjugate transpose
  REQUIRE(conj_transposed(0, 0) == std::conj(mat(0, 0)));
  REQUIRE(conj_transposed(0, 1) == std::conj(mat(1, 0)));
  REQUIRE(conj_transposed(1, 0) == std::conj(mat(0, 1)));
  REQUIRE(conj_transposed(1, 1) == std::conj(mat(1, 1)));
}

TEST_CASE("MMatrix Matrix-Vector Multiplication with Complex", "[MMatrix]") {
  using Complex = std::complex<double>;
  gra::MMatrix<Complex> mat(2, 2);
  mat(0, 0) = Complex(1.0, 1.0);
  mat(0, 1) = Complex(2.0, -1.0);
  mat(1, 0) = Complex(3.0, 4.0);
  mat(1, 1) = Complex(4.0, -2.0);

  std::vector<Complex> vec    = {Complex(1.0, 0.0), Complex(0.0, 1.0)};
  std::vector<Complex> result = mat * vec;

  REQUIRE(result.size() == 2);
  REQUIRE(result[0] == Complex(2.0, 3.0));
  REQUIRE(result[1] == Complex(5.0, 8.0));
}

TEST_CASE("MMatrix Move Assignment", "[MMatrix]") {
  gra::MMatrix<double> mat1(3, 3, 1.0);
  gra::MMatrix<double> mat2(2, 2, 2.0);

  mat2 = std::move(mat1);

  REQUIRE(mat2.size_row() == 3);
  REQUIRE(mat2.size_col() == 3);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) { REQUIRE(gra::math::IsExactEqual(mat2(i, j), 1.0)); }
  }

  mat1 = gra::MMatrix<double>{{4.0, 5.0}};
  REQUIRE(mat1(0, 1) == Approx(5.0));
}

TEST_CASE("MMatrix Identity Matrix Addition and Subtraction", "[MMatrix]") {
  gra::MMatrix<double> identity(3, 3, "eye");
  gra::MMatrix<double> mat(3, 3, 2.0);

  // Addition
  gra::MMatrix<double> result_add = mat + identity;

  REQUIRE(result_add.size_row() == 3);
  REQUIRE(result_add.size_col() == 3);

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      if (i == j) {
        REQUIRE(gra::math::IsExactEqual(result_add(i, j),
                                        3.0));  // Diagonal elements increase
      } else {
        REQUIRE(gra::math::IsExactEqual(result_add(i, j), 2.0));  // Off diagonal elements remain unchanged
      }
    }
  }

  // Subtraction
  gra::MMatrix<double> result_sub = mat - identity;

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      if (i == j) {
        REQUIRE(gra::math::IsExactEqual(result_sub(i, j),
                                        1.0));  // Diagonal elements decrease
      } else {
        REQUIRE(gra::math::IsExactEqual(result_sub(i, j), 2.0));  // Off diagonal elements remain unchanged
      }
    }
  }
}

TEST_CASE("MMatrix Scalar Operations with Complex Numbers", "[MMatrix]") {
  using Complex = std::complex<double>;
  gra::MMatrix<Complex> mat(2, 2);
  mat(0, 0) = Complex(1.0, 1.0);
  mat(0, 1) = Complex(2.0, -1.0);
  mat(1, 0) = Complex(3.0, 4.0);
  mat(1, 1) = Complex(4.0, -2.0);

  // Scalar multiplication
  gra::MMatrix<Complex> scaled = mat * Complex(2.0, 0.0);

  REQUIRE(scaled.size_row() == 2);
  REQUIRE(scaled.size_col() == 2);

  REQUIRE(scaled(0, 0) == Complex(2.0, 2.0));
  REQUIRE(scaled(0, 1) == Complex(4.0, -2.0));
  REQUIRE(scaled(1, 0) == Complex(6.0, 8.0));
  REQUIRE(scaled(1, 1) == Complex(8.0, -4.0));

  // Scalar division
  gra::MMatrix<Complex> divided = mat / Complex(2.0, 0.0);

  REQUIRE(divided.size_row() == 2);
  REQUIRE(divided.size_col() == 2);

  REQUIRE(divided(0, 0) == Complex(0.5, 0.5));
  REQUIRE(divided(0, 1) == Complex(1.0, -0.5));
  REQUIRE(divided(1, 0) == Complex(1.5, 2.0));
  REQUIRE(divided(1, 1) == Complex(2.0, -1.0));
}

TEST_CASE("MMatrix compound assignment and unary minus", "[MMatrix]") {
  gra::MMatrix<double> mat1{{1.0, 2.0}, {3.0, 4.0}};
  gra::MMatrix<double> mat2{{0.5, 1.5}, {2.5, 3.5}};

  mat1 += mat2;
  REQUIRE(mat1(0, 0) == Approx(1.5));
  REQUIRE(mat1(0, 1) == Approx(3.5));
  REQUIRE(mat1(1, 0) == Approx(5.5));
  REQUIRE(mat1(1, 1) == Approx(7.5));

  mat1 -= mat2;
  REQUIRE(mat1(0, 0) == Approx(1.0));
  REQUIRE(mat1(0, 1) == Approx(2.0));
  REQUIRE(mat1(1, 0) == Approx(3.0));
  REQUIRE(mat1(1, 1) == Approx(4.0));

  const gra::MMatrix<double> neg = -mat1;
  REQUIRE(neg(0, 0) == Approx(-1.0));
  REQUIRE(neg(0, 1) == Approx(-2.0));
  REQUIRE(neg(1, 0) == Approx(-3.0));
  REQUIRE(neg(1, 1) == Approx(-4.0));
}

TEST_CASE("MMatrix compound assignment dimension mismatch errors", "[MMatrix]") {
  gra::MMatrix<double> lhs(2, 2, 1.0);
  gra::MMatrix<double> rhs(2, 3, 1.0);

  REQUIRE_THROWS_AS(lhs += rhs, std::invalid_argument);
  REQUIRE_THROWS_AS(lhs -= rhs, std::invalid_argument);
}

TEST_CASE("MMatrix reductions and norms", "[MMatrix]") {
  gra::MMatrix<double>       mat{{1.0, 2.0}, {3.0, 4.0}};
  const gra::MMatrix<double> weights{{0.5, 1.0}, {1.5, 2.0}};

  REQUIRE(mat.Sum() == Approx(10.0));
  REQUIRE(mat.ElementwiseProductSum(weights) == Approx(15.0));
  REQUIRE(mat.FrobNorm2() == Approx(30.0));
  REQUIRE(mat.FrobNorm() == Approx(std::sqrt(30.0)));
  REQUIRE(mat.Tr() == Approx(5.0));
  REQUIRE_THROWS_AS(mat.ElementwiseProductSum(gra::MMatrix<double>(1, 2, 1.0)), std::invalid_argument);
}

// Check stable Frobenius norms over the full floating-point range
TEST_CASE("MMatrix Frobenius norm has stable scaling", "[MMatrix][LinearAlgebra]") {
  const gra::MMatrix<double> large{{3.0e200, 4.0e200}};
  REQUIRE(large.FrobNorm() == Approx(5.0e200).epsilon(1.0e-14));

  const gra::MMatrix<double> tiny{{3.0e-200, 4.0e-200}};
  REQUIRE(tiny.FrobNorm() == Approx(5.0e-200).margin(5.0e-212));
  REQUIRE(tiny.FrobNorm() > 0.0);

  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> complex{{Complex(3.0, 4.0), Complex(0.0, 12.0)}};
  REQUIRE(complex.FrobNorm() == Approx(13.0).margin(1.0e-12));

  gra::MMatrix<double> infinite{{1.0, 0.0}};
  infinite(0, 1) = std::numeric_limits<double>::infinity();
  REQUIRE(std::isinf(infinite.FrobNorm()));

  gra::MMatrix<double> not_a_number{{1.0, 0.0}};
  not_a_number(0, 1) = std::numeric_limits<double>::quiet_NaN();
  REQUIRE(std::isnan(not_a_number.FrobNorm()));
}

TEST_CASE("MMatrix calculus operations", "[MMatrix]") {
  const gra::MMatrix<double> coefficient{{2.0, 1.0}, {1.0, 3.0}};
  const gra::MMatrix<double> source{{1.0}, {2.0}};
  const gra::MMatrix<double> solution = coefficient.Solve(source);
  REQUIRE(solution(0, 0) == Approx(0.2));
  REQUIRE(solution(1, 0) == Approx(0.6));
  const gra::MMatrix<double> residual = coefficient * solution - source;
  REQUIRE(residual.FrobNorm() == Approx(0.0).margin(1e-12));

  const gra::MMatrix<double> generator{{std::log(2.0), 0.0}, {0.0, std::log(3.0)}};
  const gra::MMatrix<double> exponential = generator.Exp();
  REQUIRE(exponential.IsDiagonal());
  REQUIRE(exponential(0, 0) == Approx(2.0).epsilon(EPS));
  REQUIRE(exponential(1, 1) == Approx(3.0).epsilon(EPS));
  REQUIRE(exponential(0, 1) == Approx(0.0).margin(EPS));
  const gra::MMatrix<double> exp_inverse  = (-generator).Exp();
  const gra::MMatrix<double> exp_identity = exponential * exp_inverse;
  REQUIRE(exp_identity(0, 0) == Approx(1.0).margin(1e-10));
  REQUIRE(exp_identity(1, 1) == Approx(1.0).margin(1e-10));
  REQUIRE(exp_identity(0, 1) == Approx(0.0).margin(1e-10));
  REQUIRE(exp_identity(1, 0) == Approx(0.0).margin(1e-10));

  const gra::MMatrix<double> nilpotent{{0.0, 1.0}, {0.0, 0.0}};
  const gra::MMatrix<double> nilpotent_exponential = nilpotent.Exp();
  REQUIRE(nilpotent_exponential.IsApprox(gra::MMatrix<double>{{1.0, 1.0}, {0.0, 1.0}}, 1.0e-12));

  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> complex_generator =
      gra::MMatrix<Complex>::DiagonalMatrix({Complex(std::log(2.0), M_PI / 2.0), Complex(std::log(3.0), -M_PI)});
  const gra::MMatrix<Complex> complex_exponential = complex_generator.Exp();
  REQUIRE(complex_exponential.IsDiagonal());
  REQUIRE(complex_exponential(0, 0).real() == Approx(0.0).margin(1.0e-12));
  REQUIRE(complex_exponential(0, 0).imag() == Approx(2.0).margin(1.0e-12));
  REQUIRE(complex_exponential(1, 1).real() == Approx(-3.0).margin(1.0e-12));
  REQUIRE(complex_exponential(1, 1).imag() == Approx(0.0).margin(1.0e-12));

  const gra::MMatrix<double> left{{1.0, 2.0}, {3.0, 4.0}};
  const gra::MMatrix<double> right{{0.0, 5.0}, {6.0, 7.0}};
  const gra::MMatrix<double> product = left.Kronecker(right);
  REQUIRE(product.size_row() == 4);
  REQUIRE(product.size_col() == 4);
  REQUIRE(product(3, 2) == Approx(24.0));
  REQUIRE(product(3, 3) == Approx(28.0));

  gra::MMatrix<double> block_matrix(4, 4, 0.0);
  block_matrix.SetBlock(1, 0, left);
  const gra::MMatrix<double> block = block_matrix.Submatrix(2, 0, 2, 2);
  REQUIRE(block(0, 0) == Approx(1.0));
  REQUIRE(block(1, 1) == Approx(4.0));

  const gra::MMatrix<std::complex<double>> A{{2.0, 0.0}, {0.0, 1.0}};
  REQUIRE(A.MaxSingularValue() == Approx(2.0).epsilon(1e-8));
  const auto eigen_A     = A.ToEigen();
  const auto converted_A = gra::MMatrix<std::complex<double>>::FromEigen(eigen_A);
  REQUIRE(converted_A.IsApprox(A, 1.0e-14));
  const auto square_root = A.PrincipalPower(0.5);
  REQUIRE(square_root.IsDiagonal());
  REQUIRE((square_root * square_root).IsApprox(A, 1.0e-12));

  const gra::MMatrix<Complex> upper_triangular{{Complex(4.0, 0.0), Complex(1.0, 0.0)},
                                               {Complex(0.0, 0.0), Complex(9.0, 0.0)}};
  const auto                  upper_square_root = upper_triangular.PrincipalPower(0.5);
  REQUIRE_FALSE(upper_square_root.IsDiagonal());
  REQUIRE((upper_square_root * upper_square_root).IsApprox(upper_triangular, 1.0e-11));

  const auto zero_diagonal = gra::MMatrix<Complex>::DiagonalMatrix({Complex(1.0, 0.0), Complex(0.0, 0.0)});
  REQUIRE_THROWS(zero_diagonal.PrincipalPower(0.5));
  const auto cut_diagonal = gra::MMatrix<Complex>::DiagonalMatrix({Complex(1.0, 0.0), Complex(-1.0, 0.0)});
  REQUIRE_THROWS(cut_diagonal.PrincipalPower(0.5));

  const gra::MMatrix<double> symmetric{{2.0, 1.0}, {1.0, 2.0}};
  gra::MMatrix<double>       eigenvectors;
  const auto                 eigenvalues = symmetric.SelfAdjointEigenvalues(1.0e-12, &eigenvectors);
  REQUIRE(eigenvalues.size() == 2);
  REQUIRE(eigenvalues[0] + eigenvalues[1] == Approx(4.0).margin(1.0e-12));
  REQUIRE(eigenvalues[0] * eigenvalues[1] == Approx(3.0).margin(1.0e-12));
  REQUIRE((eigenvectors.Transpose() * eigenvectors).IsApprox(gra::MMatrix<double>::IdentityMatrix(2), 1.0e-12));
}

// Check singular values and rectangular least-squares solves
TEST_CASE("MMatrix SVD and least-squares operations", "[MMatrix][LinearAlgebra]") {
  const gra::MMatrix<double> real_design{{1.0, 0.0}, {0.0, 1.0}, {1.0, 1.0}};
  const std::vector<double>  real_expected = {2.0, -1.0};
  const std::vector<double>  real_rhs      = real_design * real_expected;
  const std::vector<double>  real_solution = real_design.SolveLeastSquares(real_rhs);

  REQUIRE(real_solution.size() == real_expected.size());
  REQUIRE(real_solution[0] == Approx(real_expected[0]).margin(1.0e-12));
  REQUIRE(real_solution[1] == Approx(real_expected[1]).margin(1.0e-12));

  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> complex_design{{Complex(1.0, 1.0), Complex(0.0, 0.0)},
                                             {Complex(0.0, 0.0), Complex(2.0, -1.0)},
                                             {Complex(1.0, 0.0), Complex(1.0, 0.0)}};
  const std::vector<Complex>  complex_expected = {Complex(0.5, -0.25), Complex(-1.0, 0.75)};
  const std::vector<Complex>  complex_rhs      = complex_design * complex_expected;
  const std::vector<Complex>  complex_solution = complex_design.SolveLeastSquares(complex_rhs);

  REQUIRE(complex_solution.size() == complex_expected.size());
  for (std::size_t index = 0; index < complex_solution.size(); ++index) {
    REQUIRE(complex_solution[index].real() == Approx(complex_expected[index].real()).margin(1.0e-12));
    REQUIRE(complex_solution[index].imag() == Approx(complex_expected[index].imag()).margin(1.0e-12));
  }

  const gra::MMatrix<double> rank_deficient{{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}};
  const std::vector<double>  rank_rhs      = {3.0, 6.0, 9.0};
  const std::vector<double>  rank_solution = rank_deficient.SolveLeastSquares(rank_rhs);
  const std::vector<double>  rank_residual = gra::Subtract(rank_deficient * rank_solution, rank_rhs);
  REQUIRE(gra::SquaredNorm(rank_residual) == Approx(0.0).margin(1.0e-20));

  const gra::MMatrix<double> singular_matrix{{0.0, 4.0, 0.0}, {0.0, 0.0, 2.0}};
  const std::vector<double>  singular_values = singular_matrix.SingularValues();
  REQUIRE(singular_values.size() == 2);
  REQUIRE(singular_values[0] == Approx(4.0).margin(1.0e-12));
  REQUIRE(singular_values[1] == Approx(2.0).margin(1.0e-12));
  REQUIRE(singular_values[0] >= singular_values[1]);

  gra::MMatrix<double> right_vectors;
  const auto rectangular_values = singular_matrix.SingularValues(&right_vectors);
  REQUIRE(right_vectors.size_row() == 3);
  REQUIRE(right_vectors.size_col() == 2);
  REQUIRE((right_vectors.Dagger() * right_vectors).IsApprox(gra::MMatrix<double>::IdentityMatrix(2), 1.0e-12));
  const auto columns = singular_matrix * right_vectors;
  for (std::size_t i = 0; i < rectangular_values.size(); ++i) {
    REQUIRE(gra::SquaredNorm(columns.Column(i)) == Approx(rectangular_values[i] * rectangular_values[i]).margin(1.0e-12));
  }

  // Resolve survival amplitudes whose squared singular values underflow
  const gra::MMatrix<Complex> absorbing{{Complex(0.0, 0.8), Complex(0.0)},
                                         {Complex(0.0), Complex(0.0, 1.0e-200)}};
  gra::MMatrix<Complex> absorbing_vectors;
  const auto absorbing_values = absorbing.SingularValues(&absorbing_vectors);
  REQUIRE(absorbing_values[0] == Approx(0.8).margin(1.0e-12));
  REQUIRE(absorbing_values[1] > 0.0);
  REQUIRE(absorbing_values[1] / 1.0e-200 == Approx(1.0).epsilon(1.0e-12));
  REQUIRE((absorbing_vectors.Dagger() * absorbing_vectors).IsApprox(gra::MMatrix<Complex>::IdentityMatrix(2), 1.0e-12));
  const auto absorbed = absorbing * absorbing_vectors.Column(1);
  REQUIRE(std::abs(absorbed[0]) <= 1.0e-210);
  REQUIRE(std::abs(absorbed[1]) / 1.0e-200 == Approx(1.0).epsilon(1.0e-12));

  REQUIRE_THROWS_AS(real_design.SolveLeastSquares({1.0, 2.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(real_design.SolveLeastSquares({1.0, std::numeric_limits<double>::quiet_NaN(), 2.0}),
                    std::invalid_argument);
}

// Check real and complex self-adjoint diagonalization
TEST_CASE("MMatrix self-adjoint eigensystems", "[MMatrix][LinearAlgebra]") {
  const gra::MMatrix<double> real_matrix{{2.0, 1.0}, {1.0, 2.0}};
  gra::MMatrix<double>       real_vectors;
  const std::vector<double>  real_values = real_matrix.SelfAdjointEigenvalues(1.0e-12, &real_vectors);
  gra::MMatrix<double>       real_diagonal(2, 2, 0.0);
  real_diagonal(0, 0) = real_values[0];
  real_diagonal(1, 1) = real_values[1];

  REQUIRE(real_values[0] == Approx(1.0).margin(1.0e-12));
  REQUIRE(real_values[1] == Approx(3.0).margin(1.0e-12));
  REQUIRE(real_values[0] <= real_values[1]);
  REQUIRE((real_vectors * real_diagonal * real_vectors.Dagger()).IsApprox(real_matrix, 1.0e-12));

  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> complex_matrix{{Complex(2.0, 0.0), Complex(1.0, 1.0)},
                                             {Complex(1.0, -1.0), Complex(3.0, 0.0)}};
  gra::MMatrix<Complex>       complex_vectors;
  const std::vector<double>   complex_values = complex_matrix.SelfAdjointEigenvalues(1.0e-12, &complex_vectors);
  gra::MMatrix<Complex>       complex_diagonal(2, 2, Complex{});
  complex_diagonal(0, 0) = Complex(complex_values[0], 0.0);
  complex_diagonal(1, 1) = Complex(complex_values[1], 0.0);

  REQUIRE(complex_values[0] == Approx(1.0).margin(1.0e-12));
  REQUIRE(complex_values[1] == Approx(4.0).margin(1.0e-12));
  REQUIRE(complex_values[0] <= complex_values[1]);
  REQUIRE((complex_vectors * complex_diagonal * complex_vectors.Dagger()).IsApprox(complex_matrix, 1.0e-12));
  REQUIRE((complex_vectors.Dagger() * complex_vectors).IsApprox(gra::MMatrix<Complex>::IdentityMatrix(2), 1.0e-12));

  const gra::MMatrix<Complex> density{{Complex(0.25, 0.0), Complex{}}, {Complex{}, Complex(0.75, 0.0)}};
  const auto                  pure_states = density.PositiveSpectralVectors(1.0e-12);
  gra::MMatrix<Complex>       reconstructed(2, 2, Complex{});
  for (const auto &state : pure_states) {
    reconstructed.AddOuterProduct(state, gra::Conjugated(state), Complex(1.0, 0.0));
  }
  REQUIRE(reconstructed.IsApprox(density, 1.0e-12));
  REQUIRE(density.SelfAdjointEntropy(1.0e-12) ==
          Approx(-0.25 * std::log(0.25) - 0.75 * std::log(0.75)).margin(1.0e-12));

  const gra::MMatrix<double> tiny_matrix{{1.0e-200, 0.0}, {0.0, 4.0e-200}};
  const std::vector<double>  tiny_values = tiny_matrix.SelfAdjointEigenvalues(1.0e-12);
  REQUIRE(tiny_values[0] == Approx(1.0e-200).margin(1.0e-212));
  REQUIRE(tiny_values[1] == Approx(4.0e-200).margin(1.0e-212));

  const gra::MMatrix<Complex> nonhermitian{{Complex(1.0, 0.0), Complex(0.0, 1.0)},
                                           {Complex(0.0, 1.0), Complex(2.0, 0.0)}};
  REQUIRE_THROWS_AS(nonhermitian.SelfAdjointEigenvalues(1.0e-12), std::invalid_argument);
  const gra::MMatrix<double> tiny_nonsymmetric{{1.0e-200, 1.0e-200}, {0.0, 1.0e-200}};
  REQUIRE_THROWS_AS(tiny_nonsymmetric.SelfAdjointEigenvalues(1.0e-12), std::invalid_argument);
  REQUIRE_THROWS_AS(real_matrix.SelfAdjointEigenvalues(-1.0), std::invalid_argument);
  const gra::MMatrix<double> indefinite{{-0.1, 0.0}, {0.0, 1.1}};
  REQUIRE_THROWS_AS(indefinite.SelfAdjointEntropy(1.0e-12), std::domain_error);
  REQUIRE_THROWS_AS(indefinite.PositiveSpectralVectors(1.0e-12), std::domain_error);
}

// Check Moore-Penrose identities and singular-value truncation
TEST_CASE("MMatrix pseudoinverse diagnostics and identities", "[MMatrix][LinearAlgebra]") {
  const gra::MMatrix<double>    matrix{{1.0, 2.0}, {0.0, 1.0}, {1.0, 0.0}};
  gra::PseudoInverseDiagnostics diagnostics;
  const gra::MMatrix<double>    inverse = matrix.PseudoInverse(0.0, &diagnostics);

  REQUIRE(diagnostics.numerical_rank == 2);
  REQUIRE(diagnostics.retained_rank == 2);
  REQUIRE(std::isfinite(diagnostics.condition_number));
  REQUIRE((matrix * inverse * matrix).IsApprox(matrix, 1.0e-11));
  REQUIRE((inverse * matrix * inverse).IsApprox(inverse, 1.0e-11));
  REQUIRE((matrix * inverse).IsHermitian(1.0e-11));
  REQUIRE((inverse * matrix).IsHermitian(1.0e-11));

  const gra::MMatrix<double>    cutoff_matrix{{1.0, 0.0}, {0.0, 1.0e-4}};
  gra::PseudoInverseDiagnostics cutoff_diagnostics;
  const gra::MMatrix<double>    cutoff_inverse = cutoff_matrix.PseudoInverse(1.0e-3, &cutoff_diagnostics);
  REQUIRE(cutoff_diagnostics.numerical_rank == 2);
  REQUIRE(cutoff_diagnostics.retained_rank == 1);
  REQUIRE(cutoff_diagnostics.condition_number == Approx(1.0e4));
  REQUIRE(cutoff_inverse(0, 0) == Approx(1.0).margin(1.0e-12));
  REQUIRE(cutoff_inverse(1, 1) == Approx(0.0).margin(0.0));

  const gra::MMatrix<double>    zero_matrix(2, 3, 0.0);
  gra::PseudoInverseDiagnostics zero_diagnostics;
  const gra::MMatrix<double>    zero_inverse = zero_matrix.PseudoInverse(0.0, &zero_diagnostics);
  REQUIRE(zero_inverse.size_row() == 3);
  REQUIRE(zero_inverse.size_col() == 2);
  REQUIRE(zero_inverse.FrobNorm() == Approx(0.0).margin(0.0));
  REQUIRE(zero_diagnostics.numerical_rank == 0);
  REQUIRE(zero_diagnostics.retained_rank == 0);
  REQUIRE(std::isinf(zero_diagnostics.condition_number));

  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> complex_matrix{{Complex(1.0, 1.0), Complex(0.0, 0.0)},
                                             {Complex(0.0, 0.0), Complex(2.0, -1.0)},
                                             {Complex(1.0, 0.0), Complex(1.0, 0.0)}};
  const gra::MMatrix<Complex> complex_inverse = complex_matrix.PseudoInverse(0.0);
  REQUIRE((complex_matrix * complex_inverse * complex_matrix).IsApprox(complex_matrix, 1.0e-11));
  REQUIRE((complex_inverse * complex_matrix * complex_inverse).IsApprox(complex_inverse, 1.0e-11));

  REQUIRE_THROWS_AS(matrix.PseudoInverse(-1.0), std::invalid_argument);
}

// Check stable positive-semidefinite square roots across matrix scales
TEST_CASE("MMatrix principal positive-semidefinite square root", "[MMatrix][LinearAlgebra]") {
  const gra::MMatrix<double> zero_matrix(3, 3, 0.0);
  const gra::MMatrix<double> zero_root = zero_matrix.PrincipalPositiveSemidefiniteSquareRoot();
  REQUIRE(zero_root.size_row() == 3);
  REQUIRE(zero_root.size_col() == 3);
  REQUIRE(zero_root.FrobNorm() == Approx(0.0).margin(0.0));

  const gra::MMatrix<double> tiny_matrix{{1.0e-200, 0.0}, {0.0, 4.0e-200}};
  const gra::MMatrix<double> tiny_root = tiny_matrix.PrincipalPositiveSemidefiniteSquareRoot();
  REQUIRE(tiny_root(0, 0) == Approx(1.0e-100).margin(1.0e-112));
  REQUIRE(tiny_root(1, 1) == Approx(2.0e-100).margin(2.0e-112));
  REQUIRE(tiny_root(0, 1) == Approx(0.0).margin(1.0e-212));
  const gra::MMatrix<double> tiny_reconstructed = tiny_root.Dagger() * tiny_root;
  REQUIRE(tiny_reconstructed(0, 0) == Approx(1.0e-200).margin(1.0e-212));
  REQUIRE(tiny_reconstructed(1, 1) == Approx(4.0e-200).margin(4.0e-212));

  const gra::MMatrix<double> numerical_negative{{-1.0e-15, 0.0}, {0.0, 1.0}};
  const gra::MMatrix<double> clamped_root = numerical_negative.PrincipalPositiveSemidefiniteSquareRoot();
  REQUIRE(clamped_root(0, 0) == Approx(0.0).margin(1.0e-15));
  REQUIRE(clamped_root(1, 1) == Approx(1.0).margin(1.0e-12));

  // A complex projector is its own square root, including its numerical null eigenspace
  for (const std::size_t dimension : {5, 9, 13}) {
    std::vector<std::complex<double>> state(dimension);
    for (const auto i : indices(state)) { state[i] = std::polar(1.0, 0.31 * i * i); }
    gra::Scale(state, 1.0 / std::sqrt(gra::SquaredNorm(state)));
    const auto projector = gra::RankOneProjector(state);
    const auto root = projector.PrincipalPositiveSemidefiniteSquareRoot();
    REQUIRE((root - projector).FrobNorm() < 1e-12);
  }

  const gra::MMatrix<double> indefinite{{-1.0e-4, 0.0}, {0.0, 1.0}};
  REQUIRE_THROWS_AS(indefinite.PrincipalPositiveSemidefiniteSquareRoot(), std::domain_error);

  const gra::MMatrix<double> nonsymmetric{{1.0, 0.1}, {0.0, 1.0}};
  REQUIRE_THROWS_AS(nonsymmetric.PrincipalPositiveSemidefiniteSquareRoot(), std::invalid_argument);
}

// Check strict and semidefinite natural-order Cholesky factors
TEST_CASE("MMatrix long-double Cholesky modes", "[MMatrix][LinearAlgebra]") {
  const long double               pivot_tolerance = 4096.0L * std::numeric_limits<long double>::epsilon();
  const gra::MMatrix<long double> strict_matrix{{1.0L, 0.0L}, {0.0L, 2.0L * pivot_tolerance}};
  const gra::MMatrix<long double> strict_lower =
      strict_matrix.CholeskyLower(pivot_tolerance, gra::CholeskyMode::StrictPositiveDefinite);
  REQUIRE(static_cast<double>(strict_lower(1, 1)) ==
          Approx(static_cast<double>(std::sqrt(2.0L * pivot_tolerance))).epsilon(1.0e-12));
  REQUIRE((strict_lower * strict_lower.Dagger()).IsApprox(strict_matrix, static_cast<double>(8.0L * pivot_tolerance)));

  const gra::MMatrix<long double> boundary_matrix{{1.0L, 0.0L}, {0.0L, pivot_tolerance}};
  REQUIRE_THROWS_AS(boundary_matrix.CholeskyLower(pivot_tolerance, gra::CholeskyMode::StrictPositiveDefinite),
                    std::domain_error);

  const gra::MMatrix<long double> natural_matrix{{4.0L, 2.0L}, {2.0L, 5.0L}};
  const gra::MMatrix<long double> natural_lower =
      natural_matrix.CholeskyLower(pivot_tolerance, gra::CholeskyMode::AllowSemidefinite);
  REQUIRE(static_cast<double>(natural_lower(0, 0)) == Approx(2.0).margin(1.0e-15));
  REQUIRE(static_cast<double>(natural_lower(1, 0)) == Approx(1.0).margin(1.0e-15));
  REQUIRE(static_cast<double>(natural_lower(1, 1)) == Approx(2.0).margin(1.0e-15));
  REQUIRE(
      (natural_lower * natural_lower.Dagger()).IsApprox(natural_matrix, static_cast<double>(8.0L * pivot_tolerance)));

  const gra::MMatrix<long double> rank_one{{1.0L, 1.0L}, {1.0L, 1.0L}};
  const gra::MMatrix<long double> semidefinite_lower =
      rank_one.CholeskyLower(pivot_tolerance, gra::CholeskyMode::AllowSemidefinite);
  REQUIRE(static_cast<double>(semidefinite_lower(0, 0)) == Approx(1.0).margin(1.0e-15));
  REQUIRE(static_cast<double>(semidefinite_lower(1, 0)) == Approx(1.0).margin(1.0e-15));
  REQUIRE(static_cast<double>(semidefinite_lower(1, 1)) == Approx(0.0).margin(0.0));
  REQUIRE((semidefinite_lower * semidefinite_lower.Dagger())
              .IsApprox(rank_one, static_cast<double>(8.0L * pivot_tolerance)));
  REQUIRE_THROWS_AS(rank_one.CholeskyLower(pivot_tolerance, gra::CholeskyMode::StrictPositiveDefinite),
                    std::domain_error);

  const gra::MMatrix<long double> leading_zero{{0.0L, 0.0L}, {0.0L, 1.0L}};
  const gra::MMatrix<long double> leading_zero_lower =
      leading_zero.CholeskyLower(pivot_tolerance, gra::CholeskyMode::AllowSemidefinite);
  REQUIRE(static_cast<double>(leading_zero_lower(0, 0)) == Approx(0.0).margin(0.0));
  REQUIRE(static_cast<double>(leading_zero_lower(1, 1)) == Approx(1.0).margin(1.0e-15));
  REQUIRE((leading_zero_lower * leading_zero_lower.Dagger())
              .IsApprox(leading_zero, static_cast<double>(8.0L * pivot_tolerance)));

  const gra::MMatrix<long double> indefinite{{1.0L, 2.0L}, {2.0L, 1.0L}};
  REQUIRE_THROWS_AS(indefinite.CholeskyLower(pivot_tolerance, gra::CholeskyMode::AllowSemidefinite), std::domain_error);

  const gra::MMatrix<long double> nonsymmetric{{1.0L, 0.1L}, {0.0L, 1.0L}};
  REQUIRE_THROWS_AS(nonsymmetric.CholeskyLower(pivot_tolerance, gra::CholeskyMode::AllowSemidefinite),
                    std::invalid_argument);
}

// Check compact row and column selection without implicit reordering
TEST_CASE("MMatrix row and column selection and row replacement", "[MMatrix][Indexing]") {
  const gra::MMatrix<double> source{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}, {7.0, 8.0, 9.0}};
  const std::vector<double>  column = source.Column(1);
  REQUIRE((column == std::vector<double>{2.0, 5.0, 8.0}));

  const gra::MMatrix<double> selected = source.SelectColumns({2, 0});
  const gra::MMatrix<double> expected_selected{{3.0, 1.0}, {6.0, 4.0}, {9.0, 7.0}};
  REQUIRE(selected.IsApprox(expected_selected));

  const gra::MMatrix<double> selected_rows = source.SelectRows({2, 0});
  const gra::MMatrix<double> expected_rows{{7.0, 8.0, 9.0}, {1.0, 2.0, 3.0}};
  REQUIRE(selected_rows.IsApprox(expected_rows));

  gra::MMatrix<double>       target(4, 3, 0.0);
  const gra::MMatrix<double> replacement{{10.0, 11.0, 12.0}, {20.0, 21.0, 22.0}};
  target.SetRows({3, 1}, replacement);
  REQUIRE(target(3, 0) == Approx(10.0));
  REQUIRE(target(3, 2) == Approx(12.0));
  REQUIRE(target(1, 0) == Approx(20.0));
  REQUIRE(target(1, 2) == Approx(22.0));
  REQUIRE(target(0, 0) == Approx(0.0));

  gra::MMatrix<double> self_alias{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  self_alias.SetRows({2, 0, 1}, self_alias);
  const gra::MMatrix<double> expected_permutation{{3.0, 4.0}, {5.0, 6.0}, {1.0, 2.0}};
  REQUIRE(self_alias.IsApprox(expected_permutation));

  REQUIRE_THROWS_AS(source.Column(3), std::out_of_range);
  REQUIRE_THROWS_AS(source.SelectColumns({0, 3}), std::out_of_range);
  REQUIRE_THROWS_AS(source.SelectRows({0, 3}), std::out_of_range);
  REQUIRE_THROWS_AS(target.SetRows({1, 1}, replacement), std::invalid_argument);
  REQUIRE_THROWS_AS(target.SetRows({0, 4}, replacement), std::out_of_range);
  REQUIRE_THROWS_AS(target.SetRows({0}, gra::MMatrix<double>(2, 3, 0.0)), std::invalid_argument);
}

// Check C++20 span support for vector Hadamard products
TEST_CASE("HadamardProduct accepts spans", "[MMatrix][Containers]") {
  const std::array<double, 4>   first  = {1.0, -2.0, 3.0, -4.0};
  const std::array<double, 4>   second = {0.5, 2.0, -1.0, -0.25};
  const std::span<const double> first_span(first);
  const std::span<const double> second_span(second);
  const std::vector<double>     product = gra::HadamardProduct(first_span, second_span);
  REQUIRE((product == std::vector<double>{0.5, -4.0, -3.0, 1.0}));

  const std::span<const double> short_span(second.data(), 3);
  REQUIRE_THROWS_AS(gra::HadamardProduct(first_span, short_span), std::invalid_argument);
}

TEST_CASE("MMatrix conjugation helpers", "[MMatrix]") {
  using Complex = std::complex<double>;
  gra::MMatrix<Complex> mat{{Complex(1.0, 2.0), Complex(-3.0, 4.0)}, {Complex(5.0, -6.0), Complex(-7.0, -8.0)}};

  const gra::MMatrix<Complex> conj = mat.Conj();
  REQUIRE(conj(0, 0) == std::conj(mat(0, 0)));
  REQUIRE(conj(0, 1) == std::conj(mat(0, 1)));
  REQUIRE(conj(1, 0) == std::conj(mat(1, 0)));
  REQUIRE(conj(1, 1) == std::conj(mat(1, 1)));
}

TEST_CASE("MMatrix Dagger is conjugate transpose alias", "[MMatrix]") {
  using Complex = std::complex<double>;
  gra::MMatrix<Complex> mat{{Complex(1.0, 2.0), Complex(3.0, -4.0), Complex(-5.0, 6.0)},
                            {Complex(7.0, 0.5), Complex(-8.0, 1.0), Complex(9.0, -2.0)}};

  const auto dagger         = mat.Dagger();
  const auto conj_transpose = mat.ConjTranspose();

  REQUIRE(dagger.size_row() == 3);
  REQUIRE(dagger.size_col() == 2);
  for (std::size_t i = 0; i < dagger.size_row(); ++i) {
    for (std::size_t j = 0; j < dagger.size_col(); ++j) { REQUIRE(dagger(i, j) == conj_transpose(i, j)); }
  }
}

TEST_CASE("MMatrix row indexing agrees with element access", "[MMatrix]") {
  gra::MMatrix<int> mat{{1, 2}, {3, 4}};
  mat[1][0] = 9;
  REQUIRE(mat(1, 0) == 9);

  const gra::MMatrix<int> &const_mat = mat;
  REQUIRE(const_mat[0][1] == 2);
  REQUIRE(const_mat(1, 0) == 9);
}

TEST_CASE("MMatrix constructor and square-matrix error paths", "[MMatrix]") {
  REQUIRE_THROWS_AS(gra::MMatrix<double>(2, 2, "bad-init"), std::invalid_argument);

  const std::size_t maximum = std::numeric_limits<std::size_t>::max();
  REQUIRE_THROWS_AS(gra::MMatrix<double>(maximum, 2), std::length_error);

  gra::MMatrix<double> resized(1, 1, 3.0);
  REQUIRE_THROWS_AS(resized.Resize(maximum, 2), std::length_error);
  REQUIRE(resized.size_row() == 1);
  REQUIRE(resized.size_col() == 1);
  REQUIRE(resized(0, 0) == Approx(3.0));

  gra::MMatrix<double> nonsquare(2, 3, 1.0);
  REQUIRE_THROWS_AS(nonsquare.Trace(), std::invalid_argument);
  REQUIRE_THROWS_AS(nonsquare.GetDiag(), std::invalid_argument);
}

TEST_CASE("MMatrix multiplication dimension mismatch errors", "[MMatrix]") {
  gra::MMatrix<double> mat{{1.0, 2.0}, {3.0, 4.0}};
  std::vector<double>  bad_vec{1.0, 2.0, 3.0};
  REQUIRE_THROWS_AS(mat * bad_vec, std::invalid_argument);

  gra::MMatrix<double> lhs(2, 3, 1.0);
  gra::MMatrix<double> rhs(2, 2, 1.0);
  REQUIRE_THROWS_AS(lhs * rhs, std::invalid_argument);
}

TEST_CASE("MMatrix Zero Size Initialization", "[MMatrix]") {
  gra::MMatrix<double> mat;

  REQUIRE(mat.size_row() == 0);
  REQUIRE(mat.size_col() == 0);
  REQUIRE(mat.isEmpty() == true);
}

TEST_CASE("MMatrix non-empty state reports false", "[MMatrix]") {
  gra::MMatrix<double> mat(1, 2, 0.0);

  REQUIRE(mat.size_row() == 1);
  REQUIRE(mat.size_col() == 2);
  REQUIRE(mat.isEmpty() == false);
}

TEST_CASE("MMatrix Diagonal Extraction", "[MMatrix]") {
  gra::MMatrix<double> mat(3, 3);
  mat(0, 0) = 1.0;
  mat(1, 1) = 2.0;
  mat(2, 2) = 3.0;

  std::vector<double> diag = mat.GetDiag();

  REQUIRE(diag.size() == 3);
  REQUIRE(gra::math::IsExactEqual(diag[0], 1.0));
  REQUIRE(gra::math::IsExactEqual(diag[1], 2.0));
  REQUIRE(gra::math::IsExactEqual(diag[2], 3.0));
}

TEST_CASE("MMatrix Element Access Out of Range", "[MMatrix]") {
  gra::MMatrix<double>        mat(3, 3);
  const gra::MMatrix<double> &const_mat = mat;

  // Check for out-of-bounds exception
  REQUIRE_THROWS_AS(mat(3, 0), std::out_of_range);
  REQUIRE_THROWS_AS(mat(0, 3), std::out_of_range);
  REQUIRE_THROWS_AS(const_mat(3, 0), std::out_of_range);
  REQUIRE_THROWS_AS(const_mat(0, 3), std::out_of_range);
}

TEST_CASE("MMatrix Initialization with Nested Initializer Lists", "[MMatrix]") {
  gra::MMatrix<int> mat = {{1, 2}, {3, 4}};

  REQUIRE(mat.size_row() == 2);
  REQUIRE(mat.size_col() == 2);

  REQUIRE(mat(0, 0) == 1);
  REQUIRE(mat(0, 1) == 2);
  REQUIRE(mat(1, 0) == 3);
  REQUIRE(mat(1, 1) == 4);
}

TEST_CASE("MMatrix sparse multiplication skips exact zeros without changing result", "[MMatrix]") {
  gra::MMatrix<double>       sparse{{0.0, 2.0, 0.0}, {3.0, 0.0, 4.0}};
  gra::MMatrix<double>       dense{{5.0, 6.0}, {7.0, 8.0}, {9.0, 10.0}};
  const gra::MMatrix<double> product = sparse * dense;

  REQUIRE(product.size_row() == 2);
  REQUIRE(product.size_col() == 2);
  REQUIRE(product(0, 0) == Approx(14.0));
  REQUIRE(product(0, 1) == Approx(16.0));
  REQUIRE(product(1, 0) == Approx(51.0));
  REQUIRE(product(1, 1) == Approx(58.0));

  const std::vector<double> vector{11.0, 12.0, 13.0};
  const std::vector<double> matvec = sparse * vector;
  REQUIRE(matvec.size() == 2);
  REQUIRE(matvec[0] == Approx(24.0));
  REQUIRE(matvec[1] == Approx(85.0));
}

TEMPLATE_TEST_CASE("MMatrix:: Initialization with initialization list", "[MMatrix][template]", int) {
  const gra::MMatrix<TestType> A = {{1, 2, 3}, {4, 5, 6}, {7, 8, 9}};

  SECTION("Two-fold initialization list") {
    REQUIRE(A[0][1] == 2);
    REQUIRE(A[2][0] == 7);
    REQUIRE(A[1][2] == 6);
  }

  const std::vector<TestType>  row0 = {1, 2, 3};
  const std::vector<TestType>  row1 = {4, 5, 6};
  const std::vector<TestType>  row2 = {7, 8, 9};
  const gra::MMatrix<TestType> B    = {row0, row1, row2};

  SECTION("One-fold initialization list with vector") {
    REQUIRE(B[0][1] == 2);
    REQUIRE(B[2][0] == 7);
    REQUIRE(B[1][2] == 6);
  }

  const gra::MMatrix<TestType> C = A * B;

  SECTION("Matrix x Matrix multiplication") {
    REQUIRE(C[0][0] == 30);
    REQUIRE(C[0][1] == 36);
    REQUIRE(C[2][1] == 126);
  }

  gra::MMatrix<TestType> D = {{1, 2, 3, 4}, {5, 6, 7, 8}};
  D                        = D.Transpose();

  SECTION("Matrix Transpose") {
    REQUIRE(D[0][1] == 5);
    REQUIRE(D[1][1] == 6);
    REQUIRE(D[2][1] == 7);
    REQUIRE(D[3][1] == 8);
  }
}

// ---------------------------------------------------------
// MTensor

TEST_CASE("MTensor initialization and basic operations", "[MTensor]") {
  using namespace gra;

  SECTION("Default constructor") {
    MTensor<double> tensor;
    REQUIRE(tensor.empty());
    REQUIRE(tensor.rank() == 0);
    REQUIRE(tensor.elements() == 0);
    REQUIRE(tensor.raw_data() == nullptr);
  }

  SECTION("Initialization with dimensions") {
    std::vector<std::size_t> dims = {4, 4, 4, 4, 4, 4};
    MTensor<double>          tensor(dims);
    REQUIRE(tensor.size(0) == 4);
    REQUIRE(tensor.size(1) == 4);
  }

  SECTION("Initialization with dimensions and default value") {
    std::vector<std::size_t> dims = {4, 4};
    MTensor<double>          tensor(dims, 2.5);
    REQUIRE(gra::math::IsExactEqual(tensor({0, 0}), 2.5));
    REQUIRE(gra::math::IsExactEqual(tensor({3, 3}), 2.5));
  }

  SECTION("Element assignment and access") {
    std::vector<std::size_t> dims = {4, 4, 4, 4, 4, 4};
    MTensor<double>          tensor(dims, 0.0);
    tensor({0, 3, 2, 0, 1, 2}) = 1.0;
    REQUIRE(gra::math::IsExactEqual(tensor({0, 3, 2, 0, 1, 2}), 1.0));
  }
}

TEST_CASE("MTensor copy constructor and assignment operator", "[MTensor]") {
  using namespace gra;

  SECTION("Copy constructor") {
    std::vector<std::size_t> dims = {3, 3, 3};
    MTensor<double>          tensor1(dims, 5.0);
    MTensor<double>          tensor2 = tensor1;

    tensor1({0, 0, 0}) = -2.0;

    REQUIRE(gra::math::IsExactEqual(tensor2({0, 0, 0}), 5.0));
    REQUIRE(gra::math::IsExactEqual(tensor2({2, 2, 2}), 5.0));
  }

  SECTION("Assignment operator") {
    std::vector<std::size_t> dims1 = {3, 3, 3};
    MTensor<double>          tensor1(dims1, 7.0);
    std::vector<std::size_t> dims2 = {4, 4, 4};
    MTensor<double>          tensor2(dims2, 1.0);

    tensor2 = tensor1;

    tensor1({2, 1, 0}) = -3.0;

    REQUIRE(gra::math::IsExactEqual(tensor2({0, 0, 0}), 7.0));
    REQUIRE(gra::math::IsExactEqual(tensor2({2, 1, 0}), 7.0));
    REQUIRE(tensor2.size(0) == 3);
    REQUIRE(tensor2.size(1) == 3);
  }
}

TEST_CASE("MTensor move, const access and empty-state handling", "[MTensor]") {
  using namespace gra;

  SECTION("Const access returns const reference") {
    using ConstRef =
        decltype(std::declval<const MTensor<double> &>()(std::declval<const std::vector<std::size_t> &>()));
    REQUIRE((std::is_same_v<ConstRef, const double &>));

    MTensor<double> tensor({2, 2}, 0.0);
    tensor({1, 1})                      = 4.5;
    const MTensor<double> &const_tensor = tensor;
    REQUIRE(gra::math::IsExactEqual(const_tensor({1, 1}), 4.5));
  }

  SECTION("Move constructor transfers ownership") {
    MTensor<double> tensor({2, 3}, 1.5);
    tensor({1, 2}) = 9.0;

    MTensor<double> moved = std::move(tensor);
    // MTensor explicitly resets the source volume during a move
    REQUIRE(tensor.empty());  // NOLINT(bugprone-use-after-move)
    REQUIRE(moved.size(0) == 2);
    REQUIRE(moved.size(1) == 3);
    REQUIRE(gra::math::IsExactEqual(moved({1, 2}), 9.0));
  }

  SECTION("Move assignment transfers ownership") {
    MTensor<double> tensor({2, 2, 2}, 0.0);
    tensor({1, 0, 1}) = 7.0;

    MTensor<double> moved({1}, 3.0);
    moved = std::move(tensor);
    // MTensor explicitly resets the source volume during a move
    REQUIRE(tensor.empty());  // NOLINT(bugprone-use-after-move)
    REQUIRE(moved.rank() == 3);
    REQUIRE(gra::math::IsExactEqual(moved({1, 0, 1}), 7.0));
  }

  SECTION("Assignment from empty tensor resets state") {
    MTensor<double> tensor({2, 2}, 5.0);
    tensor = MTensor<double>();
    REQUIRE(tensor.empty());
    REQUIRE(tensor.rank() == 0);
    REQUIRE_THROWS_AS(tensor({}), std::invalid_argument);
  }

  SECTION("Empty tensor access throws") {
    MTensor<double> tensor;
    REQUIRE(tensor.empty());
    REQUIRE_THROWS_AS(tensor({}), std::invalid_argument);
  }
}

TEST_CASE("MTensor raw storage follows its public row-major indexing contract", "[MTensor]") {
  using namespace gra;

  MTensor<double> tensor({2, 3, 4}, 0.0);

  double *base = tensor.raw_data();
  REQUIRE(base != nullptr);

  REQUIRE(&tensor({0, 0, 0}) == base);
  REQUIRE(&tensor({0, 0, 1}) == base + 1);
  REQUIRE(&tensor({0, 0, 3}) == base + 3);

  REQUIRE(&tensor({0, 1, 0}) == base + 4);
  REQUIRE(&tensor({0, 2, 0}) == base + 8);

  REQUIRE(&tensor({1, 0, 0}) == base + 12);
  REQUIRE(&tensor({1, 2, 3}) == base + 23);
}

TEST_CASE("MTensor index validation", "[MTensor]") {
  using namespace gra;

  SECTION("Out-of-bounds indexing throws exception") {
    std::vector<std::size_t> dims = {3, 3, 3};
    MTensor<double>          tensor(dims, 0.0);

    REQUIRE_THROWS_AS(tensor({3, 0, 0}), std::invalid_argument);
    REQUIRE_THROWS_AS(tensor({0, 3, 0}), std::invalid_argument);
    REQUIRE_THROWS_AS(tensor({0, 0, 3}), std::invalid_argument);
  }

  SECTION("Incorrect rank throws exception") {
    std::vector<std::size_t> dims = {3, 3, 3};
    MTensor<double>          tensor(dims, 0.0);

    REQUIRE_THROWS_AS(tensor({1, 2}),
                      std::invalid_argument);  // Fewer dimensions
    REQUIRE_THROWS_AS(tensor({1, 2, 3, 4}),
                      std::invalid_argument);  // More dimensions
    REQUIRE_THROWS_AS(tensor.size(3), std::out_of_range);
  }
}

// ---------------------------------------------------------
// M4Vec

TEST_CASE("Constructor and Initialization", "[M4Vec]") {
  // Test default constructor
  M4Vec vec1;
  REQUIRE(vec1.Px() == Approx(0.0).epsilon(EPS));
  REQUIRE(vec1.Py() == Approx(0.0).epsilon(EPS));
  REQUIRE(vec1.Pz() == Approx(0.0).epsilon(EPS));
  REQUIRE(vec1.E() == Approx(0.0).epsilon(EPS));

  // Test constructor with values
  M4Vec vec2(1.0, 2.0, 3.0, 4.0);
  REQUIRE(vec2.Px() == Approx(1.0).epsilon(EPS));
  REQUIRE(vec2.Py() == Approx(2.0).epsilon(EPS));
  REQUIRE(vec2.Pz() == Approx(3.0).epsilon(EPS));
  REQUIRE(vec2.E() == Approx(4.0).epsilon(EPS));
}

TEST_CASE("Setters and Getters", "[M4Vec]") {
  M4Vec vec;

  // Test individual setters and getters
  vec.SetPx(5.0);
  vec.SetPy(6.0);
  vec.SetPz(7.0);
  vec.SetE(8.0);

  REQUIRE(vec.Px() == Approx(5.0).epsilon(EPS));
  REQUIRE(vec.Py() == Approx(6.0).epsilon(EPS));
  REQUIRE(vec.Pz() == Approx(7.0).epsilon(EPS));
  REQUIRE(vec.E() == Approx(8.0).epsilon(EPS));

  // Test Set method
  vec.Set(9.0, 10.0, 11.0, 12.0);
  REQUIRE(vec.Px() == Approx(9.0).epsilon(EPS));
  REQUIRE(vec.Py() == Approx(10.0).epsilon(EPS));
  REQUIRE(vec.Pz() == Approx(11.0).epsilon(EPS));
  REQUIRE(vec.E() == Approx(12.0).epsilon(EPS));
}

TEST_CASE("Invariants and Physics Quantities", "[M4Vec]") {
  // Create a 4-vector with components (px, py, pz, E)
  M4Vec vec(1.0, 2.0, 3.0, 10.0);

  // Test Invariant mass (E^2 - |p|^2)
  // |p|^2 = px^2 + py^2 + pz^2 = 1^2 + 2^2 + 3^2 = 14
  // Invariant mass M^2 = E^2 - |p|^2 = 10^2 - 14 = 100 - 14 = 86
  REQUIRE(vec.Invariant() == Approx(86.0).epsilon(EPS));
  REQUIRE(vec.M2() == Approx(86.0).epsilon(EPS));
  REQUIRE(vec.M() == Approx(sqrt(86.0)).epsilon(EPS));

  // Test transverse momentum
  // Pt^2 = px^2 + py^2 = 1^2 + 2^2 = 5
  REQUIRE(vec.Pt2() == Approx(5.0).epsilon(EPS));
  REQUIRE(vec.Pt() == Approx(sqrt(5.0)).epsilon(EPS));

  // Test pseudorapidity and rapidity from their logarithmic definitions
  const double momentum          = sqrt(14.0);
  const double expected_eta      = 0.5 * log((momentum + 3.0) / (momentum - 3.0));
  const double expected_rapidity = 0.5 * log((10.0 + 3.0) / (10.0 - 3.0));
  REQUIRE(vec.Eta() == Approx(expected_eta).epsilon(EPS));
  REQUIRE(vec.Rap() == Approx(expected_rapidity).epsilon(EPS));

  // Test gamma (E/m) and beta (|p|/E)
  REQUIRE(vec.Gamma() == Approx(10.0 / sqrt(86.0)).epsilon(EPS));
  REQUIRE(vec.Beta() == Approx(momentum / 10.0).epsilon(EPS));
}

TEST_CASE("Operators and Algebra", "[M4Vec]") {
  M4Vec vec1(1.0, 2.0, 3.0, 4.0);
  M4Vec vec2(2.0, 1.0, 0.5, 3.0);

  // Test addition
  M4Vec sum = vec1 + vec2;
  REQUIRE(sum.Px() == Approx(3.0).epsilon(EPS));
  REQUIRE(sum.Py() == Approx(3.0).epsilon(EPS));
  REQUIRE(sum.Pz() == Approx(3.5).epsilon(EPS));
  REQUIRE(sum.E() == Approx(7.0).epsilon(EPS));

  // Test subtraction
  M4Vec diff = vec1 - vec2;
  REQUIRE(diff.Px() == Approx(-1.0).epsilon(EPS));
  REQUIRE(diff.Py() == Approx(1.0).epsilon(EPS));
  REQUIRE(diff.Pz() == Approx(2.5).epsilon(EPS));
  REQUIRE(diff.E() == Approx(1.0).epsilon(EPS));

  // Test scalar multiplication
  M4Vec scaled = vec1 * 2.0;
  REQUIRE(scaled.Px() == Approx(2.0).epsilon(EPS));
  REQUIRE(scaled.Py() == Approx(4.0).epsilon(EPS));
  REQUIRE(scaled.Pz() == Approx(6.0).epsilon(EPS));
  REQUIRE(scaled.E() == Approx(8.0).epsilon(EPS));

  // Test scalar division
  M4Vec divided = vec1 / 2.0;
  REQUIRE(divided.Px() == Approx(0.5).epsilon(EPS));
  REQUIRE(divided.Py() == Approx(1.0).epsilon(EPS));
  REQUIRE(divided.Pz() == Approx(1.5).epsilon(EPS));
  REQUIRE(divided.E() == Approx(2.0).epsilon(EPS));
}

TEST_CASE("Dot and Cross Products", "[M4Vec]") {
  M4Vec vec1(1.0, 0.0, 0.0, 5.0);
  M4Vec vec2(0.0, 1.0, 0.0, 4.0);

  // Test Minkowski dot product
  REQUIRE(vec1.DotM(vec2) == Approx(20.0).epsilon(EPS));

  // Test 3-vector dot product
  REQUIRE(vec1.Dot3(vec2) == Approx(0.0).epsilon(EPS));

  // Test cross product
  M4Vec cross = vec1.Cross3(vec2);
  REQUIRE(cross.Px() == Approx(0.0).epsilon(EPS));
  REQUIRE(cross.Py() == Approx(0.0).epsilon(EPS));
  REQUIRE(cross.Pz() == Approx(1.0).epsilon(EPS));
  REQUIRE(cross.E() == Approx(0.0).epsilon(EPS));
  REQUIRE(cross.Dot3(vec1) == Approx(0.0).margin(EPS));
  REQUIRE(cross.Dot3(vec2) == Approx(0.0).margin(EPS));

  const auto extended = vec1.Contravariant<long double>();
  STATIC_REQUIRE(std::is_same_v<typename decltype(extended)::value_type, long double>);
  REQUIRE(extended[0] == 5.0L);
  REQUIRE(extended[1] == 1.0L);
  REQUIRE(extended[2] == 0.0L);
  REQUIRE(extended[3] == 0.0L);

  const std::array<long double, 4> extended_vec2 = {4.0L, 0.0L, 1.0L, 0.0L};
  REQUIRE(gra::MinkowskiProduct(extended, extended_vec2) == 20.0L);
  const std::array<double, 4> ordinary_vec2 = {4.0, 0.0, 1.0, 0.0};
  auto                        mixed_product = gra::MinkowskiProduct(ordinary_vec2, extended);
  STATIC_REQUIRE(std::is_same_v<decltype(mixed_product), long double>);
  REQUIRE(mixed_product == 20.0L);
  REQUIRE_THROWS_AS(gra::MinkowskiProduct(std::array<double, 3>{}, extended_vec2), std::invalid_argument);
}

TEST_CASE("Indexing Operators", "[M4Vec]") {
  M4Vec vec(1.0, 2.0, 3.0, 4.0);

  // Test contravariant indexing
  REQUIRE(vec ^ 0 == Approx(4.0).epsilon(EPS));  // E
  REQUIRE(vec ^ 1 == Approx(1.0).epsilon(EPS));  // px
  REQUIRE(vec ^ 2 == Approx(2.0).epsilon(EPS));  // py
  REQUIRE(vec ^ 3 == Approx(3.0).epsilon(EPS));  // pz

  // Test covariant indexing (with metric)
  REQUIRE(vec % 0 == Approx(4.0).epsilon(EPS));   // E
  REQUIRE(vec % 1 == Approx(-1.0).epsilon(EPS));  // -px
  REQUIRE(vec % 2 == Approx(-2.0).epsilon(EPS));  // -py
  REQUIRE(vec % 3 == Approx(-3.0).epsilon(EPS));  // -pz

  double contraction = 0.0;
  for (std::size_t mu = 0; mu < 4; ++mu) { contraction += (vec ^ mu) * (vec % mu); }
  REQUIRE(contraction == Approx(vec.M2()).margin(EPS));
}

TEST_CASE("Comparison Operators", "[M4Vec]") {
  M4Vec vec1(1.0, 2.0, 3.0, 4.0);
  M4Vec vec2(1.0, 2.0, 3.0, 4.0);
  M4Vec vec3(1.0, 2.0, 3.0, 5.0);

  // Test equality
  REQUIRE(vec1 == vec2);
  REQUIRE(vec1 != vec3);
}

TEST_CASE("Physics Quantities", "[M4Vec]") {
  M4Vec vec(3.0, 4.0, 0.0, 5.0);

  // Test transverse mass
  REQUIRE(vec.Mt2() == Approx(25.0).epsilon(EPS));
  REQUIRE(vec.Mt() == Approx(5.0).epsilon(EPS));

  // Test azimuthal angle
  REQUIRE(vec.Phi() == Approx(atan2(4.0, 3.0)).epsilon(EPS));
}

// Basic operations of 4-vectors
TEST_CASE("M4Vec: Basic kinematic operations", "[M4Vec]") {
  const double EPS = 1e-5;
  const double mpi = 0.139;  // GeV

  M4Vec p(0, 0, 0, mpi);
  SECTION("Invariant operators") {
    REQUIRE(p.M() == Approx(mpi).epsilon(EPS));
    REQUIRE(p.M2() == Approx(mpi * mpi).epsilon(EPS));
  }

  M4Vec k(1.0, 2.0, 3.0, mpi);
  SECTION("3-Momentum components") {
    REQUIRE(k.Px() == Approx(1.0).epsilon(EPS));
    REQUIRE(k.Py() == Approx(2.0).epsilon(EPS));
    REQUIRE(k.Pz() == Approx(3.0).epsilon(EPS));
  }

  SECTION("Bracket indexing: k[1] == k.Px()") {
    REQUIRE(k[1] == Approx(k.Px()).epsilon(EPS));
    REQUIRE(k[2] == Approx(k.Py()).epsilon(EPS));
    REQUIRE(k[3] == Approx(k.Pz()).epsilon(EPS));
    REQUIRE(k[0] == Approx(k.E()).epsilon(EPS));
  }

  SECTION("Covariant indexing with %% operator: ") {
    for (std::size_t i = 1; i <= 3; ++i) { REQUIRE((k % i) == Approx(-k[i]).epsilon(EPS)); }
    REQUIRE((k % 0) == Approx(k[0]).epsilon(EPS));
  }
}

TEST_CASE("M4Vec: Setter aliases and helper accessors", "[M4Vec]") {
  M4Vec vec;

  vec.SetX(1.0);
  vec.SetY(2.0);
  vec.SetZ(3.0);
  vec.SetT(4.0);
  REQUIRE(vec.X() == Approx(1.0).epsilon(EPS));
  REQUIRE(vec.Y() == Approx(2.0).epsilon(EPS));
  REQUIRE(vec.Z() == Approx(3.0).epsilon(EPS));
  REQUIRE(vec.T() == Approx(4.0).epsilon(EPS));

  vec.SetPxPy(5.0, 6.0);
  REQUIRE(vec.Px() == Approx(5.0).epsilon(EPS));
  REQUIRE(vec.Py() == Approx(6.0).epsilon(EPS));

  vec.SetPxPyPz(7.0, 8.0, 9.0);
  REQUIRE(vec.Px() == Approx(7.0).epsilon(EPS));
  REQUIRE(vec.Py() == Approx(8.0).epsilon(EPS));
  REQUIRE(vec.Pz() == Approx(9.0).epsilon(EPS));

  vec.SetPzE(10.0, 11.0);
  REQUIRE(vec.Pz() == Approx(10.0).epsilon(EPS));
  REQUIRE(vec.E() == Approx(11.0).epsilon(EPS));

  vec.SetP3({1.0, -2.0, 3.5});
  const auto p3 = vec.P3();
  REQUIRE(p3.size() == 3);
  REQUIRE(p3[0] == Approx(1.0).epsilon(EPS));
  REQUIRE(p3[1] == Approx(-2.0).epsilon(EPS));
  REQUIRE(p3[2] == Approx(3.5).epsilon(EPS));

  vec.SetPxPyPzM(3.0, 4.0, 12.0, 5.0);
  REQUIRE(vec.E() == Approx(13.9283882772).epsilon(1e-10));

  const auto beta = vec.BetaVector();
  REQUIRE(beta.size() == 3);
  REQUIRE(beta[0] == Approx(vec.Px() / vec.E()).epsilon(EPS));
  REQUIRE(beta[1] == Approx(vec.Py() / vec.E()).epsilon(EPS));
  REQUIRE(beta[2] == Approx(vec.Pz() / vec.E()).epsilon(EPS));
}

TEST_CASE("M4Vec: Derived quantities and phase helpers", "[M4Vec]") {
  M4Vec vec(3.0, 4.0, 0.0, 5.0);

  REQUIRE(vec.Perp() == Approx(vec.Pt()).epsilon(EPS));
  REQUIRE(vec.Perp2() == Approx(vec.Pt2()).epsilon(EPS));
  REQUIRE(vec.Et() == Approx(vec.Mt()).epsilon(EPS));
  REQUIRE(vec.Et2() == Approx(vec.Mt2()).epsilon(EPS));
  REQUIRE(vec.Theta() == Approx(M_PI / 2.0).epsilon(1e-12));
  REQUIRE(vec.CosTheta() == Approx(0.0).margin(1e-12));
  REQUIRE(vec.LightconePos() == Approx(5.0).epsilon(EPS));
  REQUIRE(vec.LightconeNeg() == Approx(5.0).epsilon(EPS));
  REQUIRE(vec.ComplexPt().real() == Approx(3.0).epsilon(EPS));
  REQUIRE(vec.ComplexPt().imag() == Approx(4.0).epsilon(EPS));
  REQUIRE(vec.ExpCPhi().real() == Approx(0.6).epsilon(EPS));
  REQUIRE(vec.ExpCPhi().imag() == Approx(0.8).epsilon(EPS));

  M4Vec ref(4.0, -3.0, 0.0, 6.0);
  REQUIRE(vec.DotPt(ref) == Approx(0.0).epsilon(EPS));
  REQUIRE(vec.DeltaPhi(ref) == Approx(M_PI / 2.0).margin(1e-12));
  REQUIRE(vec.DeltaPhiAbs(ref) == Approx(M_PI / 2.0).margin(1e-12));
}

TEST_CASE("M4Vec: Mutating operators and safety checks", "[M4Vec]") {
  M4Vec vec(1.0, -2.0, 3.0, 4.0);
  M4Vec rhs(0.5, 1.5, -1.0, 2.0);

  vec += rhs;
  REQUIRE(vec == M4Vec(1.5, -0.5, 2.0, 6.0));

  vec -= rhs;
  REQUIRE(vec == M4Vec(1.0, -2.0, 3.0, 4.0));

  vec *= 2.0;
  REQUIRE(vec == M4Vec(2.0, -4.0, 6.0, 8.0));

  vec /= 4.0;
  REQUIRE(vec == M4Vec(0.5, -1.0, 1.5, 2.0));

  const auto neg = -vec;
  REQUIRE(neg == M4Vec(-0.5, 1.0, -1.5, -2.0));

  vec.Flip3();
  REQUIRE(vec == M4Vec(-0.5, 1.0, -1.5, 2.0));

  vec[0] = 9.0;
  vec[1] = 8.0;
  vec[2] = 7.0;
  vec[3] = 6.0;
  REQUIRE(vec.E() == Approx(9.0).epsilon(EPS));
  REQUIRE(vec.Px() == Approx(8.0).epsilon(EPS));
  REQUIRE(vec.Py() == Approx(7.0).epsilon(EPS));
  REQUIRE(vec.Pz() == Approx(6.0).epsilon(EPS));
  REQUIRE_THROWS_AS(vec[4], std::out_of_range);
  REQUIRE_THROWS_AS((vec ^ 4), std::out_of_range);
  REQUIRE_THROWS_AS((vec % 4), std::out_of_range);
}

TEST_CASE("M4Vec: Lorentz boost and propagation edge cases", "[M4Vec]") {
  const M4Vec vec(1.0, 2.0, 3.0, 10.0);

  const gra::M3Vec beta{0.21, -0.13, 0.08};
  const M4Vec      boosted  = vec.LorentzBoost(beta);
  const M4Vec      restored = boosted.LorentzBoost(beta, -1);
  REQUIRE(boosted.M2() == Approx(vec.M2()).margin(1e-11));
  REQUIRE(restored.Px() == Approx(vec.Px()).margin(1e-11));
  REQUIRE(restored.Py() == Approx(vec.Py()).margin(1e-11));
  REQUIRE(restored.Pz() == Approx(vec.Pz()).margin(1e-11));
  REQUIRE(restored.E() == Approx(vec.E()).margin(1e-11));

  const M4Vec invalid = vec.LorentzBoost({1.1, 0.0, 0.0});
  REQUIRE(invalid == M4Vec(0.0, 0.0, 0.0, -1.0));

  const M4Vec rest(0.0, 0.0, 0.0, 2.0);
  const M4Vec pos = rest.PropagatePosition(3.0, 2.0);
  REQUIRE(pos.Px() == Approx(0.0).epsilon(EPS));
  REQUIRE(pos.Py() == Approx(0.0).epsilon(EPS));
  REQUIRE(pos.Pz() == Approx(0.0).epsilon(EPS));
  REQUIRE(pos.E() == Approx(3.0 * 2.99792458E8 * 2.0).epsilon(1e-12));
}

// ---------------------------------------------------------
// Canonical tensor and vector algebra

TEST_CASE("FlatTensorTranspose exchanges flattened tensor indices", "[MMatrix]") {
  std::size_t A = 2, B = 3, C = 5;

  // Create a (A x B) x C matrix representing a 3D tensor
  MMatrix<int> input({{1, 2, 3, 4, 5},
                      {6, 7, 8, 9, 10},
                      {11, 12, 13, 14, 15},
                      {16, 17, 18, 19, 20},
                      {21, 22, 23, 24, 25},
                      {26, 27, 28, 29, 30}});

  // Expected output of size (A x C) x B
  MMatrix<int> expected({{1, 6, 11},
                         {2, 7, 12},
                         {3, 8, 13},
                         {4, 9, 14},
                         {5, 10, 15},
                         {16, 21, 26},
                         {17, 22, 27},
                         {18, 23, 28},
                         {19, 24, 29},
                         {20, 25, 30}});

  // Perform transpose
  MMatrix<int> result = input.FlatTensorTranspose(A, B, C);

  // Verify that the result matches the expected output
  REQUIRE(result.IsApprox(expected, 0.0));
}

TEST_CASE("Kronecker maps matrix basis indices to the tensor product", "[MMatrix]") {
  gra::MMatrix<double> A{{1, 2}, {3, 4}};
  gra::MMatrix<double> B{{0, 5}, {6, 7}};
  gra::MMatrix<double> C = A.Kronecker(B);

  gra::MMatrix<double> expectedC{{0, 5, 0, 10}, {6, 7, 12, 14}, {0, 15, 0, 20}, {18, 21, 24, 28}};

  REQUIRE(C.IsApprox(expectedC));

  gra::MMatrix<double> sparseA{{0, 2}, {3, 0}};
  gra::MMatrix<double> sparseC = sparseA.Kronecker(B);
  gra::MMatrix<double> expectedSparseC{{0, 0, 0, 10}, {0, 0, 12, 14}, {0, 15, 0, 0}, {18, 21, 0, 0}};
  REQUIRE(sparseC.IsApprox(expectedSparseC));
}

TEST_CASE("KroneckerProduct preserves vector basis ordering", "[MMatrix]") {
  const std::vector<double> first    = {1.0, 2.0};
  const std::vector<double> second   = {3.0, 4.0, 5.0};
  const std::vector<double> expected = {3.0, 4.0, 5.0, 6.0, 8.0, 10.0};

  REQUIRE(gra::KroneckerProduct(first, second) == expected);
}

// Check an allocation-free Kronecker-vector contraction and its basis order
TEST_CASE("KroneckerBilinearProduct contracts the lower index fastest", "[MMatrix][Kronecker]") {
  using Complex                       = std::complex<double>;
  const std::array<double, 2>  upper  = {2.0, -1.0};
  const std::array<double, 3>  lower  = {0.5, 3.0, -2.0};
  const std::array<Complex, 6> tensor = {Complex(1.0, 0.2), Complex(-0.4, 0.3), Complex(0.7, -0.1),
                                         Complex(0.2, 0.5), Complex(-1.0, 0.4), Complex(0.6, 0.8)};
  const Complex expected = upper[0] * (lower[0] * tensor[0] + lower[1] * tensor[1] + lower[2] * tensor[2]) +
                           upper[1] * (lower[0] * tensor[3] + lower[1] * tensor[4] + lower[2] * tensor[5]);

  const std::span<const double> upper_view(upper);
  const Complex                 actual = gra::KroneckerBilinearProduct(upper_view, lower, tensor);
  REQUIRE(std::abs(actual - expected) == Approx(0.0).margin(1e-14));
  REQUIRE_THROWS_AS(gra::KroneckerBilinearProduct(upper, lower, std::array<Complex, 5>{}), std::invalid_argument);
}

// Check mixed scalar products, exact row maps and fused matrix multiplication
TEST_CASE("Kronecker supports mixed scalars and mapped fused multiplication", "[MMatrix][Kronecker]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<double>  upper{{1.0, 2.0}, {-0.5, 3.0}};
  const gra::MMatrix<Complex> lower{{Complex(0.4, 0.2), Complex(-1.0, 0.3)}, {Complex(2.0, -0.1), Complex(0.5, 0.7)}};
  const auto                  tensor = upper.Kronecker(lower);
  static_assert(std::is_same_v<std::remove_cvref_t<decltype(tensor(0, 0))>, Complex>);

  const std::array<std::size_t, 4> destination_rows = {3, 1, 2, 0};
  const auto                       mapped           = upper.Kronecker(lower, destination_rows);
  for (std::size_t source = 0; source < destination_rows.size(); ++source) {
    for (std::size_t col = 0; col < tensor.size_col(); ++col) {
      REQUIRE(std::abs(mapped(destination_rows[source], col) - tensor(source, col)) == Approx(0.0).margin(1e-14));
    }
  }

  const gra::MMatrix<double> right{{1.0, -0.2}, {0.3, 2.0}, {-1.0, 0.5}, {0.7, -0.4}};
  const auto                 fused = upper.KroneckerMultiply(lower, right, destination_rows);
  const auto expected              = mapped * gra::MMatrix<Complex>{{1.0, -0.2}, {0.3, 2.0}, {-1.0, 0.5}, {0.7, -0.4}};
  REQUIRE(fused.IsApprox(expected, 1e-14));
  const auto empty = upper.KroneckerMultiply(lower, gra::MMatrix<double>(tensor.size_col(), 0), destination_rows);
  REQUIRE(empty.size_row() == destination_rows.size());
  REQUIRE(empty.size_col() == 0);

  REQUIRE_THROWS_AS(upper.Kronecker(lower, std::array<std::size_t, 4>{0, 1, 1, 3}), std::invalid_argument);
  REQUIRE_THROWS_AS(upper.Kronecker(lower, std::array<std::size_t, 4>{0, 1, 2, 4}), std::invalid_argument);
  REQUIRE_THROWS_AS(upper.Kronecker(lower, std::array<int, 4>{0, 1, 2, -1}), std::invalid_argument);
  REQUIRE_THROWS_AS(upper.Kronecker(lower, std::array<std::size_t, 3>{0, 1, 2}), std::invalid_argument);
  REQUIRE_THROWS_AS(upper.KroneckerMultiply(lower, gra::MMatrix<double>(3, 1), destination_rows),
                    std::invalid_argument);
}

// Check fused tensor contraction entrywise with complex nonsymmetric matrices
TEST_CASE("KroneckerMultiply preserves the exact tensor contraction order", "[MMatrix][Kronecker]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> left{{{1.0, 0.3}, {-0.2, 0.7}, {0.5, -0.4}}, {{-0.8, 0.1}, {0.6, -0.9}, {0.2, 0.5}}};
  const gra::MMatrix<Complex> right_factor{
      {{0.4, -0.6}, {1.2, 0.1}}, {{-0.3, 0.8}, {0.7, -0.2}}, {{0.9, 0.5}, {-0.4, -0.7}}, {{0.1, -0.9}, {0.3, 0.6}}};
  gra::MMatrix<Complex> vertex(6, 5, Complex{});
  for (std::size_t row = 0; row < vertex.size_row(); ++row) {
    for (std::size_t col = 0; col < vertex.size_col(); ++col) {
      vertex[row][col] =
          Complex(0.11 * static_cast<double>(1 + row) - 0.03 * col, -0.07 * row + 0.05 * static_cast<double>(1 + col));
    }
  }

  const auto actual = left.KroneckerMultiply(right_factor, vertex);
  REQUIRE(actual.size_row() == 8);
  REQUIRE(actual.size_col() == 5);
  for (std::size_t left_row = 0; left_row < left.size_row(); ++left_row) {
    for (std::size_t right_row = 0; right_row < right_factor.size_row(); ++right_row) {
      for (std::size_t output_col = 0; output_col < vertex.size_col(); ++output_col) {
        Complex expected = 0.0;
        for (std::size_t left_col = 0; left_col < left.size_col(); ++left_col) {
          for (std::size_t right_col = 0; right_col < right_factor.size_col(); ++right_col) {
            expected += left[left_row][left_col] * right_factor[right_row][right_col] *
                        vertex[left_col * right_factor.size_col() + right_col][output_col];
          }
        }
        REQUIRE(std::abs(actual[left_row * right_factor.size_row() + right_row][output_col] - expected) < 1e-13);
      }
    }
  }
  REQUIRE(actual.IsApprox(left.Kronecker(right_factor) * vertex, 1e-13));
  REQUIRE_THROWS_AS(left.KroneckerMultiply(right_factor, gra::MMatrix<Complex>(5, 1, 0.0)), std::invalid_argument);
}

// Check in-place accumulation of one tensor-product vector into a column
TEST_CASE("AddColumnKroneckerProduct uses the lower index fastest", "[MMatrix][Kronecker]") {
  gra::MMatrix<double>        matrix(6, 3, 0.0);
  const std::array<double, 2> first{2.0, -3.0};
  const std::array<double, 3> second{5.0, 7.0, -11.0};
  matrix.AddColumnKroneckerProduct(1, first, second);
  matrix.AddColumnKroneckerProduct(1, first, second);
  REQUIRE(matrix.Column(1) == std::vector<double>{20.0, 28.0, -44.0, -30.0, -42.0, 66.0});
  REQUIRE_THROWS_AS(matrix.AddColumnKroneckerProduct(3, first, second), std::out_of_range);
  REQUIRE_THROWS_AS(matrix.AddColumnKroneckerProduct(0, first, std::array<double, 2>{1.0, 2.0}), std::invalid_argument);
}

// Compare complex row permutations with an explicit tensor product
TEST_CASE("Row permutations preserve complex tensor columns", "[MMatrix][Kronecker]") {
  using Complex                                = std::complex<double>;
  const std::array<Complex, 2>     first       = {Complex{0.3, -0.7}, Complex{-0.2, 0.4}};
  const std::array<double, 3>      second      = {0.0, 1.7, -0.9};
  const std::array<std::size_t, 6> destination = {3, 0, 5, 2, 1, 4};
  const Complex                    scale{0.6, -0.5};
  gra::MMatrix<Complex>            matrix(6, 2, 0.0);
  matrix.AddColumnKroneckerProduct(1, first, second, scale);
  matrix.AddColumnKroneckerProduct(1, first, second, 2.0 * scale);
  const auto mapped = matrix.PermuteRows(destination);
  const auto tensor = gra::KroneckerProduct(first, second);
  for (std::size_t i = 0; i < tensor.size(); ++i) {
    CHECK(std::abs(mapped[destination[i]][1] - 3.0 * scale * tensor[i]) < 1.0e-14);
    CHECK(std::abs(mapped[i][0]) < 1.0e-14);
  }
  REQUIRE(gra::MMatrix<Complex>(6, 0).PermuteRows(destination).size_col() == 0);
  REQUIRE_THROWS_AS(matrix.PermuteRows(std::array<std::size_t, 6>{}), std::invalid_argument);
  REQUIRE_THROWS_AS(matrix.PermuteRows(std::array<std::size_t, 2>{}), std::invalid_argument);
}

// Check inactive sparse tensor rows are never read by the fused product
TEST_CASE("KroneckerMultiply skips exact-zero tensor coefficients", "[MMatrix][Kronecker]") {
  const double                     nan              = std::numeric_limits<double>::quiet_NaN();
  const std::array<std::size_t, 1> destination_rows = {0};

  const gra::MMatrix<double> zero_upper{{2.0, 0.0}};
  const gra::MMatrix<double> dense_lower{{3.0, 4.0}};
  const gra::MMatrix<double> upper_inactive_nan{{5.0}, {6.0}, {nan}, {nan}};
  const auto upper_result = zero_upper.KroneckerMultiply(dense_lower, upper_inactive_nan, destination_rows);
  REQUIRE(upper_result(0, 0) == Approx(78.0));
  REQUIRE(std::isfinite(upper_result(0, 0)));

  const gra::MMatrix<double> dense_upper{{2.0, 4.0}};
  const gra::MMatrix<double> zero_lower{{3.0, 0.0}};
  const gra::MMatrix<double> lower_inactive_nan{{5.0}, {nan}, {6.0}, {nan}};
  const auto lower_result = dense_upper.KroneckerMultiply(zero_lower, lower_inactive_nan, destination_rows);
  REQUIRE(lower_result(0, 0) == Approx(102.0));
  REQUIRE(std::isfinite(lower_result(0, 0)));

  const gra::MMatrix<double> subnormal_upper{{std::numeric_limits<double>::denorm_min()}};
  const gra::MMatrix<double> finite_lower{{0.25}};
  const gra::MMatrix<double> underflow_inactive_nan{{nan}};
  const auto                 underflow_result =
      subnormal_upper.KroneckerMultiply(finite_lower, underflow_inactive_nan, destination_rows);
  REQUIRE(gra::math::IsZero(underflow_result(0, 0)));
  REQUIRE(std::isfinite(underflow_result(0, 0)));
}

TEST_CASE("Kronecker rejects overflowing matrix dimensions", "[MMatrix][validation]") {
  const gra::MMatrix<double> maximum_rows(std::numeric_limits<std::size_t>::max(), 0);
  const gra::MMatrix<double> two_rows(2, 0);

  REQUIRE_THROWS_AS(maximum_rows.Kronecker(two_rows), std::length_error);
}

TEST_CASE("Kronecker preserves trace and dagger product identities", "[MMatrix]") {
  using cdouble = std::complex<double>;

  const gra::MMatrix<cdouble> A{{cdouble(1.0, 0.2), cdouble(-0.4, 0.7)}, {cdouble(0.3, -0.5), cdouble(2.0, -0.1)}};
  const gra::MMatrix<cdouble> B{{cdouble(-0.2, 0.6), cdouble(1.2, -0.3)}, {cdouble(0.8, 0.4), cdouble(0.5, 0.9)}};

  const auto AB = A.Kronecker(B);
  REQUIRE(AB.Trace().real() == Approx((A.Trace() * B.Trace()).real()).margin(1e-12));
  REQUIRE(AB.Trace().imag() == Approx((A.Trace() * B.Trace()).imag()).margin(1e-12));

  const auto lhs = AB.Dagger();
  const auto rhs = A.Dagger().Kronecker(B.Dagger());
  REQUIRE(lhs.size_row() == rhs.size_row());
  REQUIRE(lhs.size_col() == rhs.size_col());

  for (std::size_t i = 0; i < lhs.size_row(); ++i) {
    for (std::size_t j = 0; j < lhs.size_col(); ++j) {
      REQUIRE(lhs[i][j].real() == Approx(rhs[i][j].real()).margin(1e-12));
      REQUIRE(lhs[i][j].imag() == Approx(rhs[i][j].imag()).margin(1e-12));
    }
  }
}

TEST_CASE("MeshGrid broadcasts each coordinate over the Cartesian product", "[MMatrix]") {
  std::vector<double> x = {1.0, 2.0, 3.0};
  std::vector<double> y = {4.0, 5.0};
  const auto [X, Y]     = gra::MeshGrid(x, y);

  gra::MMatrix<double> expectedX{{1.0, 2.0, 3.0}, {1.0, 2.0, 3.0}};
  gra::MMatrix<double> expectedY{{4.0, 4.0, 4.0}, {5.0, 5.0, 5.0}};

  REQUIRE(X.IsApprox(expectedX));
  REQUIRE(Y.IsApprox(expectedY));
}

TEST_CASE("Diagonal products apply independent row and column weights", "[MMatrix]") {
  std::vector<double>  x = {1.0, 2.0};
  gra::MMatrix<double> A{{3.0, 4.0}, {5.0, 6.0}};
  std::vector<double>  y = {7.0, 8.0};

  const gra::MMatrix<double> C = A.LeftDiagonalProduct(x).RightDiagonalProduct(y);

  gra::MMatrix<double> expectedC{{21.0, 32.0}, {70.0, 96.0}};

  REQUIRE(C.IsApprox(expectedC));
}

TEST_CASE("DiagonalMatrix embeds a vector on the matrix diagonal", "[MMatrix]") {
  std::vector<double>  x = {1.0, 2.0, 3.0};
  gra::MMatrix<double> A = gra::MMatrix<double>::DiagonalMatrix(x);

  gra::MMatrix<double> expectedA{{1.0, 0.0, 0.0}, {0.0, 2.0, 0.0}, {0.0, 0.0, 3.0}};

  REQUIRE(A.IsApprox(expectedA));
}

TEST_CASE("NormalizedL2 preserves vector direction", "[MMatrix]") {
  std::vector<double> x      = {3.0, 4.0};
  std::vector<double> unit_x = gra::NormalizedL2(x);

  std::vector<double> expected_x = {0.6, 0.8};

  REQUIRE(unit_x.size() == expected_x.size());
  REQUIRE(unit_x[0] == Approx(expected_x[0]).margin(1e-12));
  REQUIRE(unit_x[1] == Approx(expected_x[1]).margin(1e-12));
  REQUIRE(gra::MMatrix<double>::IdentityMatrix(2).BilinearForm(unit_x, unit_x) == Approx(1.0).margin(1e-12));
}

TEST_CASE("NormalizedSum preserves positive vector ratios", "[MMatrix]") {
  const std::vector<double> input      = {2.0, 3.0, 5.0};
  const std::vector<double> normalized = gra::NormalizedSum(input);

  REQUIRE(gra::Sum(normalized) == Approx(1.0).margin(1e-12));
  REQUIRE(normalized[0] == Approx(0.2).margin(1e-12));
  REQUIRE(normalized[1] == Approx(0.3).margin(1e-12));
  REQUIRE(normalized[2] == Approx(0.5).margin(1e-12));
}

TEST_CASE("Orthonormal vector operators build and apply one basis", "[MMatrix]") {
  using Complex = std::complex<double>;
  std::vector<std::vector<Complex>> basis;
  REQUIRE(gra::AddOrthonormal(basis, std::vector<Complex>{1.0, 1.0}, 1e-12));
  REQUIRE(gra::AddOrthonormal(basis, std::vector<Complex>{1.0, -1.0}, 1e-12));
  REQUIRE_FALSE(gra::AddOrthonormal(basis, std::vector<Complex>{2.0, 0.0}, 1e-12));
  const std::vector<Complex> source    = {Complex(0.2, 0.4), Complex(-0.3, 0.1)};
  const auto                 projected = gra::ProjectOrthonormal(basis, source);
  REQUIRE(projected.size() == source.size());
  for (const auto &i : gra::aux::indices(source)) { REQUIRE(std::abs(projected[i] - source[i]) < 1e-12); }
}

TEST_CASE("BilinearForm agrees with explicit contraction", "[MMatrix]") {
  std::vector<double>  a = {1.0, 2.0};
  gra::MMatrix<double> B{{3.0, 4.0}, {5.0, 6.0}};
  std::vector<double>  c = {7.0, 8.0};

  double result          = B.BilinearForm(a, c);
  double expected_result = 219.0;

  REQUIRE(result == Approx(expected_result).epsilon(1e-9));
}

// Check mixed real-matrix complex-vector contraction without conjugation
TEST_CASE("BilinearForm supports mixed nonconjugating Lorentz contractions", "[MMatrix][complex]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<double>   metric{{1.0, 0.0}, {0.0, -1.0}};
  const std::array<Complex, 2> left      = {Complex(1.0, 1.0), Complex(2.0, -1.0)};
  const std::array<Complex, 2> right     = {Complex(3.0, 2.0), Complex(-1.0, 4.0)};
  const Complex                expected  = left[0] * right[0] - left[1] * right[1];
  const Complex                hermitian = std::conj(left[0]) * right[0] - std::conj(left[1]) * right[1];

  const Complex actual = metric.BilinearForm(left, right);
  REQUIRE(std::abs(actual - expected) == Approx(0.0).margin(1e-15));
  REQUIRE(std::abs(actual - hermitian) > 1.0);
}

// Check vector and matrix bilinear forms with independent coefficient masks
TEST_CASE("Masked bilinear products preserve nonconjugating semantics", "[MMatrix][complex]") {
  using Complex                                = std::complex<double>;
  const std::array<Complex, 3> first           = {Complex(1.0, 1.0), Complex(2.0, -1.0), Complex(-1.0, 2.0)};
  const std::array<double, 3>  second          = {2.0, 3.0, 4.0};
  const std::array<bool, 3>    common_mask     = {true, false, true};
  const Complex                vector_expected = first[0] * second[0] + first[2] * second[2];
  REQUIRE(std::abs(gra::MaskedBilinearProduct(first, second, common_mask) - vector_expected) ==
          Approx(0.0).margin(1e-14));
  const std::array<double, 3> second_with_inactive_nan = {2.0, std::numeric_limits<double>::quiet_NaN(), 4.0};
  REQUIRE(std::abs(gra::MaskedBilinearProduct(first, second_with_inactive_nan, common_mask) - vector_expected) ==
          Approx(0.0).margin(1e-14));

  const gra::MMatrix<double>   matrix{{1.0, 2.0}, {3.0, 4.0}, {5.0, 6.0}};
  const std::array<double, 3>  left            = {0.5, -1.0, 2.0};
  const std::array<Complex, 2> right           = {Complex(0.2, 0.3), Complex(-0.4, 0.7)};
  const std::array<bool, 3>    row_mask        = {true, false, true};
  const std::array<bool, 2>    column_mask     = {false, true};
  const Complex                matrix_expected = left[0] * matrix(0, 1) * right[1] + left[2] * matrix(2, 1) * right[1];
  REQUIRE(std::abs(matrix.MaskedBilinearForm(left, right, row_mask, column_mask) - matrix_expected) ==
          Approx(0.0).margin(1e-14));
  auto matrix_with_inactive_nan  = matrix;
  matrix_with_inactive_nan(1, 1) = std::numeric_limits<double>::quiet_NaN();
  matrix_with_inactive_nan(0, 0) = std::numeric_limits<double>::quiet_NaN();
  REQUIRE(std::abs(matrix_with_inactive_nan.MaskedBilinearForm(left, right, row_mask, column_mask) - matrix_expected) ==
          Approx(0.0).margin(1e-14));

  REQUIRE_THROWS_AS(gra::MaskedBilinearProduct(first, second, std::array<bool, 2>{true, false}), std::invalid_argument);
  REQUIRE_THROWS_AS(matrix.MaskedBilinearForm(left, right, std::array<bool, 2>{true, false}, column_mask),
                    std::invalid_argument);
}

TEST_CASE("Normalization helpers reject zero and nonfinite input", "[MMatrix][validation]") {
  REQUIRE_THROWS_AS(gra::NormalizedL2(std::vector<double>{0.0, 0.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::NormalizedL2(std::vector<double>{1.0, std::numeric_limits<double>::quiet_NaN()}),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::NormalizedSum(std::vector<double>{-1.0, 1.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::NormalizedSum(std::vector<double>{1.0, std::numeric_limits<double>::infinity()}),
                    std::invalid_argument);
}

TEST_CASE("Bilinear and distance helpers reject dimension mismatch", "[MMatrix][validation]") {
  const std::vector<double>  short_vector = {1.0};
  const std::vector<double>  long_vector  = {2.0, 3.0};
  const gra::MMatrix<double> matrix{{1.0, 0.0}, {0.0, 1.0}};

  REQUIRE_THROWS_AS(gra::BilinearProduct(short_vector, long_vector), std::invalid_argument);
  REQUIRE_THROWS_AS(matrix.LeftMultiply(short_vector), std::invalid_argument);
  REQUIRE_THROWS_AS(matrix.BilinearForm(short_vector, long_vector), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::L1Distance(short_vector, long_vector), std::invalid_argument);
}

TEST_CASE("LinearInterpolate validates ordered finite grids", "[MMath][validation]") {
  const std::vector<double> nodes  = {0.0, 1.0, 2.0};
  const std::vector<double> values = {0.0, 1.0, 2.0};
  REQUIRE_NOTHROW(gra::math::ValidateInterpolationGrid(nodes));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(nodes, values, 1.5) == Approx(1.5).margin(1e-15));
  REQUIRE_THROWS_AS(gra::math::LinearInterpolate(std::vector<double>{0.0, 2.0, 1.0}, values, 0.5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::LinearInterpolate(
                        std::vector<double>{0.0, std::numeric_limits<double>::quiet_NaN(), 2.0}, values, 0.5),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::LinearInterpolate(std::vector<double>{0.0, 1.0}, std::vector<double>{0.0, 1.0},
                                                 std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
}

TEST_CASE("Validated interpolation inverts CDF plateaus and small probabilities", "[MMath][interpolation]") {
  const std::vector<double> cdf    = {0.0, 1.0e-20, 0.5, 0.5, 1.0};
  const std::vector<double> radius = {0.0, 1.0, 2.0, 3.0, 4.0};
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 0.0) == Approx(0.0));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 5.0e-21) == Approx(0.5));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 0.25) == Approx(1.5));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 0.5) == Approx(3.0));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 0.75) == Approx(3.5));
  REQUIRE(gra::math::LinearInterpolateValidatedGrid(cdf, radius, 1.0) == Approx(4.0));
}

TEST_CASE("LinearRadialIntegral exactly integrates linear radial fields", "[MMath][integration]") {
  const std::vector<double> nodes = {0.0, 0.3, 1.1, 2.0};
  std::vector<double>       values(nodes.size(), 0.0);
  for (const auto &i : indices(nodes)) { values[i] = 2.0 + 3.0 * nodes[i]; }
  REQUIRE(gra::math::LinearRadialIntegral(nodes, values) == Approx(12.0).margin(2.0e-14));
  REQUIRE_THROWS_AS(gra::math::LinearRadialIntegral(nodes, std::vector<double>{1.0}), std::invalid_argument);
}

TEST_CASE("CubicLagrangeInterpolate reproduces complex cubic fields", "[MMath][interpolation]") {
  const std::array<double, 4> nodes = {-1.3, -0.2, 0.7, 2.1};
  const auto                  field = [](const double x) {
    return std::complex<double>(1.2 - 0.4 * x + 0.7 * x * x - 0.3 * x * x * x,
                                -0.8 + 0.6 * x - 0.2 * x * x + 0.5 * x * x * x);
  };
  const std::array<std::complex<double>, 4> values  = {field(nodes[0]), field(nodes[1]), field(nodes[2]),
                                                       field(nodes[3])};
  const double                              x       = 0.31;
  const auto                                weights = gra::math::CubicLagrangeWeights(nodes, x);
  REQUIRE(std::accumulate(weights.cbegin(), weights.cend(), 0.0) == Approx(1.0).margin(2.0e-15));
  REQUIRE(std::abs(gra::math::CubicLagrangeWeightedSum(values, weights) - field(x)) < 2.0e-14);
  REQUIRE(std::abs(gra::math::CubicLagrangeInterpolate(nodes, values, x) - field(x)) < 2.0e-14);

  auto repeated = nodes;
  repeated[2]   = repeated[1];
  REQUIRE_THROWS_AS(gra::math::CubicLagrangeWeights(repeated, x), std::invalid_argument);
}

TEST_CASE("Interpolation cells and uniform cubic stencils preserve grids", "[MMath][interpolation]") {
  const std::vector<double> nodes = {0.0, 1.0, 3.0};
  const auto                below = gra::math::LocateCell(nodes, -2.0);
  REQUIRE(below.index == 0);
  REQUIRE(below.fraction == Approx(0.0));
  const auto inside = gra::math::LocateCell(nodes, 2.0);
  REQUIRE(inside.index == 1);
  REQUIRE(inside.fraction == Approx(0.5));
  const auto above = gra::math::LocateCell(nodes, 4.0);
  REQUIRE(above.index == 1);
  REQUIRE(above.fraction == Approx(1.0));
  REQUIRE_THROWS_AS(gra::math::LocateCell(std::vector<double>{1.0}, 1.0), std::invalid_argument);

  constexpr std::size_t count   = 7;
  constexpr double      minimum = -1.0;
  constexpr double      maximum = 2.0;
  constexpr double      x       = 0.37;
  const auto            stencil = gra::math::UniformCubicWeights(count, minimum, maximum, x);
  const double          step    = (maximum - minimum) / (count - 1);
  const auto            field   = [](const double coordinate) {
    return 0.4 - 0.8 * coordinate + 0.2 * coordinate * coordinate + 0.7 * coordinate * coordinate * coordinate;
  };
  std::array<double, 4> value{};
  for (std::size_t i = 0; i < value.size(); ++i) {
    REQUIRE(stencil.index[i] < count);
    value[i] = field(minimum + stencil.index[i] * step);
  }
  REQUIRE(gra::math::CubicLagrangeWeightedSum(value, stencil.weight) == Approx(field(x)).margin(2.0e-14));
  REQUIRE_THROWS_AS(gra::math::UniformCubicWeights(3, minimum, maximum, x), std::invalid_argument);
}

// Compare logarithmic evaluations across the overflow boundary to extended precision
TEST_CASE("Log Bessel I1 preserves finite high excitation densities", "[MMath][Bessel]") {
  for (const double x : {1.0e-10, 0.1, 10.0, 354.0, 355.0, 700.0, 1000.0, 4000.0}) {
    const double reference = std::log(std::cyl_bessel_i(1.0L, static_cast<long double>(x)));
    CHECK(gra::math::LogBesselI1(x) == Approx(reference).margin(2.0e-12));
  }
}

// Check independent Bessel values with explicit absolute tolerances
TEST_CASE("BesselJ012 shares stable integer-order evaluations", "[MMath][Bessel]") {
  for (const double x : {0.0, 1.0e-7, -1.0e-4, 0.02, 0.7, 1.0, 2.0, 10.0}) {
    const auto value = gra::math::BesselJ012(x);
    REQUIRE(gra::math::IsExactEqual(value[0], gra::math::BesselJ0(x)));
    REQUIRE(gra::math::IsExactEqual(value[1], gra::math::BesselJ1(x)));
    REQUIRE(value[2] == Approx(gra::math::BesselJ(2, x)).margin(2.0e-10));
    if (x >= 0.0) {
      REQUIRE(value[0] == Approx(std::cyl_bessel_j(0.0, x)).epsilon(0.0).margin(8.0e-9));
      REQUIRE(value[1] == Approx(std::cyl_bessel_j(1.0, x)).epsilon(0.0).margin(8.0e-9));
      REQUIRE(value[2] == Approx(std::cyl_bessel_j(2.0, x)).epsilon(0.0).margin(8.0e-9));
    }
  }
  for (int order = 0; order <= 8; ++order) {
    const double x        = 3.2;
    const double positive = gra::math::BesselJ(order, x);
    const double negative = gra::math::BesselJ(order, -x);
    const double parity   = order % 2 == 0 ? 1.0 : -1.0;
    REQUIRE(positive == Approx(std::cyl_bessel_j(static_cast<double>(order), x)).margin(8.0e-9));
    REQUIRE(negative == Approx(parity * positive).margin(2.0e-14));
  }
  REQUIRE_THROWS_AS(gra::math::BesselJ(-1, 0.4), std::invalid_argument);
  REQUIRE(std::isnan(gra::math::BesselJ(2, std::numeric_limits<double>::infinity())));
}

TEST_CASE("MMatrix IsApprox rejects nonfinite tolerance", "[MMatrix][validation]") {
  const gra::MMatrix<double> matrix{{1.0}};
  REQUIRE_FALSE(matrix.IsApprox(matrix, std::numeric_limits<double>::quiet_NaN()));
}

TEST_CASE("MMatrix IsHermitian rejects nonfinite inputs", "[MMatrix][validation]") {
  using Complex = std::complex<double>;
  const gra::MMatrix<Complex> hermitian{{Complex(1.0, 0.0), Complex(0.0, 1.0)},
                                        {Complex(0.0, -1.0), Complex(2.0, 0.0)}};
  REQUIRE_FALSE(hermitian.IsHermitian(std::numeric_limits<double>::quiet_NaN()));

  auto nonfinite  = hermitian;
  nonfinite(0, 0) = Complex(std::numeric_limits<double>::quiet_NaN(), 0.0);
  REQUIRE_FALSE(nonfinite.IsHermitian(1.0e-12));
}

// Check that small complex fluctuations survive a large coherent mean
TEST_CASE("Complex moments preserve centered fluctuations and phases", "[statistics]") {
  const std::vector<std::complex<double>> sample{{1.0e8 - 1.0, -2.0e8 - 2.0}, {1.0e8 + 1.0, -2.0e8 + 2.0}};
  const auto                              moments = statistics::ComplexMoments(sample);
  REQUIRE(moments.variance == Approx(5.0).epsilon(0.0).margin(1.0e-14));
  auto rotated = sample;
  gra::Scale(rotated, std::complex<double>(0.0, 1.0));
  REQUIRE(statistics::ComplexMoments(rotated).variance == Approx(moments.variance).epsilon(0.0).margin(1.0e-14));
  REQUIRE(statistics::ComplexMoments({{1.0e8, -2.0e8}, {1.0e8, -2.0e8}}).variance == Approx(0.0).margin(1.0e-14));
}

// Check inverse powers and logarithmic ranges without materializing large values
TEST_CASE("Powered variance preserves inverse powers and logarithmic shifts", "[statistics]") {
  const std::vector<double> logs{0.0, std::log(2.0), std::log(4.0)};
  REQUIRE(statistics::PoweredVariance(logs, -1.0) == Approx(2.0 / 7.0).epsilon(1.0e-13));
  REQUIRE(statistics::PoweredVariance({-400.0, 400.0}, -1.0) == Approx(1.0).epsilon(1.0e-14));
  REQUIRE(statistics::PoweredVariance({-400.0, 400.0}, 1.0) == Approx(1.0).epsilon(1.0e-14));
  REQUIRE(statistics::PoweredVariance({-400.0, 400.0}, 0.0) == Approx(0.0).margin(1.0e-14));
  REQUIRE(statistics::PoweredVariance({600.0, 1400.0}, -1.0) == Approx(1.0).epsilon(1.0e-14));
  REQUIRE_THROWS_AS(statistics::PoweredVariance({std::numeric_limits<double>::quiet_NaN()}, 0.0),
                    std::invalid_argument);
}

// Check finite Gamma tails and explicit failure when the root budget is insufficient
TEST_CASE("Gamma quantiles enforce root convergence and positive support", "[MDistribution]") {
  const math::GammaNumerics limited{1.0e-13, 32, 1.0e12};
  REQUIRE(math::GammaCDF(1.0, std::log(2.0), limited) == Approx(0.5).epsilon(1.0e-12));
  REQUIRE_THROWS_AS(math::GammaQuantile(0.5, 1.0, limited), std::runtime_error);
  const math::GammaNumerics exact{1.0e-13, 10000, 1.0e12};
  REQUIRE(math::GammaQuantile(0.5, 1.0, exact) == Approx(std::log(2.0)).epsilon(0.0).margin(1.0e-13));
  for (const auto shape : {10.0, 100.0}) {
    const math::GammaNumerics approximate{1.0e-13, 10000, shape};
    const double              probability = shape < 50.0 ? 1.0e-30 : 1.0e-300;
    const double              quantile    = math::GammaQuantile(probability, shape, approximate);
    REQUIRE(std::isfinite(quantile));
    REQUIRE(quantile > 0.0);
    REQUIRE(math::GammaCDF(shape, quantile, exact) / probability == Approx(1.0).epsilon(2.0e-10));
  }
}

// Check terminating polynomials and exact nonpolynomial identities near unity
TEST_CASE("Gauss functions preserve polynomials and continuation", "[MMath][hypergeometric]") {
  for (const double z : {0.0, 0.5, 0.9999, 1.0}) {
    REQUIRE(math::Hyper2F1(-1.0, 3.0, 1.0, z) == Approx(1.0 - 3.0 * z).epsilon(0.0).margin(1.0e-14));
    REQUIRE(math::RegularizedHyper2F1(-1.0, 3.0, 1.0, z) == Approx(1.0 - 3.0 * z).epsilon(0.0).margin(1.0e-14));
    REQUIRE(math::RegularizedHyper2F1(-1.0, 2.0, 0.0, z) == Approx(-2.0 * z).epsilon(0.0).margin(1.0e-14));
  }
  for (const double z : {0.5, 0.75, 0.99, 0.9999, 1.0 - 1.0e-12}) {
    const double logarithm = -std::log1p(-z) / z;
    REQUIRE(math::Hyper2F1(1.0, 1.0, 2.0, z) == Approx(logarithm).epsilon(2.0e-12));
    REQUIRE(math::RegularizedHyper2F1(1.0, 1.0, 2.0, z) == Approx(logarithm).epsilon(2.0e-12));
    const double power = std::pow(1.0 - z, -0.25);
    REQUIRE(math::Hyper2F1(0.25, 3.0, 3.0, z) == Approx(power).epsilon(2.0e-12));
    REQUIRE(math::RegularizedHyper2F1(0.25, 3.0, 3.0, z) == Approx(power / 2.0).epsilon(2.0e-12));
  }
  REQUIRE_THROWS_AS(math::Hyper2F1(1000.0, 1000.0, 1.0, 0.75), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(math::RegularizedHyper2F1(1000.0, 1000.0, 1.0, 0.75), gra::AmplitudeFailure);
}

// Check scalar and empty tensor states through the copy assignment API
TEST_CASE("Tensor copy assignment distinguishes a scalar from empty storage", "[MTensor]") {
  const std::vector<std::size_t> scalar_index;
  const MTensor<double>          scalar(scalar_index, 7.0);
  const MTensor<double>          empty;
  MTensor<double>                target;
  target = scalar;
  REQUIRE(target.elements() == 1);
  REQUIRE(target(scalar_index) == Approx(7.0));
  target = empty;
  REQUIRE(target.empty());
  REQUIRE(target.raw_data() == nullptr);
  REQUIRE_THROWS_AS(target(scalar_index), std::invalid_argument);
  target = scalar;
  REQUIRE(target.elements() == 1);
  REQUIRE(target(scalar_index) == Approx(7.0));
  MTensor<double> shaped({2, 3}, 1.0);
  shaped = empty;
  REQUIRE(shaped.empty());
  REQUIRE(shaped.rank() == 0);
  shaped = scalar;
  REQUIRE(shaped(scalar_index) == Approx(7.0));
}

// Check signed differences and both container overloads of a closed sequence
TEST_CASE("Linspace preserves descending integer grids", "[MMath][MAlgorithms]") {
  const std::vector<int> expected{3, 2, 1, 0};
  REQUIRE(math::linspace(3, 0, 4) == expected);
  REQUIRE(math::linspace<std::vector>(3, 0, 4) == expected);
  REQUIRE(math::linspace(3U, 0U, 4) == std::vector<unsigned int>{3, 2, 1, 0});
  REQUIRE(math::linspace(2, -3, 4) == std::vector<int>{2, 1, -1, -3});
  REQUIRE(math::linspace(-3, 2, 4) == std::vector<int>{-3, -2, 0, 2});
  const int minimum = std::numeric_limits<int>::min();
  const int maximum = std::numeric_limits<int>::max();
  REQUIRE(math::linspace(minimum, maximum, 3) == std::vector<int>{minimum, -1, maximum});
  REQUIRE(math::linspace(maximum, minimum, 3) == std::vector<int>{maximum, 0, minimum});
  const auto array = math::linspace<std::valarray>(3, 0, 4);
  for (const auto &i : indices(expected)) { REQUIRE(array[i] == expected[i]); }
  const auto extended = math::linspace(0.0L, 1.0L, 4);
  REQUIRE(std::abs(extended[1] - 1.0L / 3.0L) <= 4.0L * std::numeric_limits<long double>::epsilon());
}

// Check that changing the common units does not change a linear solution
TEST_CASE("Matrix solve preserves common scaling", "[MMatrix]") {
  const MMatrix<double> matrix{{2.0, 1.0}, {1.0, 3.0}};
  const MMatrix<double> expected{{1.0, -2.0}, {3.0, 4.0}};
  const auto            source = matrix * expected;
  for (const double scale : {1.0e-200, 1.0e-20, 1.0, 1.0e200}) {
    REQUIRE((matrix * scale).Solve(source * scale).IsApprox(expected, 1.0e-12));
    const MMatrix<double> singular{{scale, scale}, {scale, scale}};
    REQUIRE_THROWS_AS(singular.Solve(source * scale), std::runtime_error);
  }
}

// Check kernel scale invariance and exactly zero marginal support
TEST_CASE("Sinkhorn preserves common kernel scaling", "[MMath][transport]") {
  const std::vector<double> p{0.5, 0.5};
  const std::vector<double> q{0.2, 0.3, 0.5};
  MMatrix<double>           coupling;
  for (const double scale : {1.0e-200, 1.0, 1.0e200}) {
    const MMatrix<double> kernel(2, 3, scale);
    REQUIRE(opt::SinkHorn(coupling, kernel, p, q, 5) == Approx(0.0).margin(1.0e-14));
    REQUIRE(coupling.IsApprox(gra::OuterProduct(p, q), 1.0e-14));
  }
  const MMatrix<double> sparse{{1.0, 2.0}, {0.0, 0.0}};
  REQUIRE(opt::SinkHorn(coupling, sparse, {1.0, 0.0}, {0.2, 0.8}, 5) == Approx(0.0).margin(1.0e-14));
  REQUIRE(coupling.IsApprox(MMatrix<double>{{0.2, 0.8}, {0.0, 0.0}}, 1.0e-14));
  REQUIRE_THROWS_AS(opt::SinkHorn(coupling, sparse, p, {0.2, 0.8}, 5), std::runtime_error);
}

// Check fluctuations below the precision of uncentered exponential moments
TEST_CASE("Powered variance resolves very small powers", "[statistics]") {
  for (const double power : {1.0e-20, -1.0e-20, 1.0e-8, -1.0e-8}) {
    const double expected = std::pow(std::tanh(power / 2.0), 2);
    REQUIRE(statistics::PoweredVariance({0.0, 1.0}, power) / expected == Approx(1.0).epsilon(1.0e-13));
    REQUIRE(statistics::PoweredVariance({100.0, 101.0}, power) / expected == Approx(1.0).epsilon(1.0e-13));
  }
}

// Check logarithmic reductions independently of transport balancing
TEST_CASE("Log product sums preserve finite and zero support", "[MMatrix]") {
  const double              infinity = std::numeric_limits<double>::infinity();
  const std::vector<double> a{1000.0, -1000.0};
  const std::vector<double> b{-1000.0, 1000.0};
  REQUIRE(static_cast<double>(gra::LogProductSum(a, b)) == Approx(std::log(2.0)).epsilon(1.0e-14));
  const std::vector<double> zero{-infinity, -infinity};
  const auto                value = gra::LogProductSum(a, zero);
  REQUIRE(std::isinf(value));
  REQUIRE(std::signbit(value));
  REQUIRE_THROWS_AS(gra::LogProductSum(a, std::vector<double>{0.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::LogProductSum(a, std::vector<double>{infinity, 0.0}), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::LogProductSum(a, std::vector<double>{std::numeric_limits<double>::quiet_NaN(), 0.0}),
                    std::invalid_argument);
}

// Check that positive kernel support survives an extreme range of matrix entries
TEST_CASE("Sinkhorn retains subnormal kernel support", "[MMath][transport]") {
  MMatrix<double> coupling;
  for (const double scale : {1.0e-200, 1.0, 1.0e200}) {
    const MMatrix<double> kernel{{scale, 4.0 * scale}, {9.0 * scale, 16.0 * scale}};
    REQUIRE(opt::SinkHorn(coupling, kernel, {0.5, 0.5}, {0.5, 0.5}, 100) == Approx(0.0).margin(1.0e-13));
    REQUIRE(coupling.IsApprox(MMatrix<double>{{0.2, 0.3}, {0.3, 0.2}}, 1.0e-13));
  }
  for (const double small : {1.0e-200, std::numeric_limits<double>::denorm_min()}) {
    const MMatrix<double> kernel{{small, 0.0}, {0.0, 1.0e200}};
    REQUIRE(opt::SinkHorn(coupling, kernel, {0.5, 0.5}, {0.5, 0.5}, 5) == Approx(0.0).margin(1.0e-13));
    REQUIRE(coupling.IsApprox(MMatrix<double>{{0.5, 0.0}, {0.0, 0.5}}, 1.0e-13));
  }
}

// Check elimination without overflowing intermediate differences
TEST_CASE("Matrix solve normalizes large finite coefficients", "[MMatrix]") {
  const MMatrix<double> matrix{{1.0e308, 1.0e308}, {-1.0e308, 1.0e308}};
  REQUIRE(matrix.Solve(matrix).IsApprox(MMatrix<double>{{1.0, 0.0}, {0.0, 1.0}}, 1.0e-14));
  const MMatrix<double> expected{{1.0, -2.0}, {3.0, 4.0}};
  const MMatrix<double> mixed{{2.0e200, 1.0e200}, {1.0e-100, 3.0e-100}};
  REQUIRE(mixed.Solve(mixed * expected).IsApprox(expected, 1.0e-13));
  const MMatrix<double> diagonal{{1.0e200, 0.0}, {0.0, 1.0e-10}};
  REQUIRE(diagonal.Solve(diagonal * expected).IsApprox(expected, 1.0e-13));
}

// Check Gamma support below the smallest normal floating point value
TEST_CASE("Gamma functions resolve subnormal tails", "[MDistribution]") {
  const math::GammaNumerics exact{1.0e-13, 10000, 1.0e12};
  REQUIRE(math::GammaCDF(1.0e-310, 0.5, exact) == Approx(1.0).epsilon(1.0e-13));
  const double probability = 1.0e-310;
  REQUIRE(math::GammaQuantile(probability, 1.0, exact) / probability == Approx(1.0).epsilon(2.0e-12));
  REQUIRE(math::IsZero(math::GammaQuantile(1.0e-10, 1.0e-3, exact)));
}

// Check Gamma ratios, small regularized values and decaying Gauss solutions
TEST_CASE("Gauss functions retain small values and signal overflow", "[MMath][hypergeometric]") {
  for (const double x : {-2.5, -1.5, -0.5, 0.5, 5.0}) {
    REQUIRE(math::ReciprocalGamma(x) * std::tgamma(x) == Approx(1.0).epsilon(1.0e-14));
  }
  REQUIRE(math::Hyper2F1(1.0, 1.0, 200.0, 1.0) == Approx(199.0 / 198.0).epsilon(1.0e-13));
  const double reciprocal = math::ReciprocalGamma(172.0);
  REQUIRE(reciprocal > 0.0);
  REQUIRE(std::log(reciprocal) + gra::math::LogGamma(172.0) == Approx(0.0).margin(3.0e-13));
  for (const double z : {0.25, 0.5, 0.99}) {
    const double expected = math::Hyper2F1(1.0, 1.0, 20.0, z) / std::tgamma(20.0);
    REQUIRE(math::RegularizedHyper2F1(1.0, 1.0, 20.0, z) / expected == Approx(1.0).epsilon(2.0e-12));
  }
  for (const double z : {0.99, 1.0 - 1.0e-12}) {
    const double expected = std::pow(1.0 - z, 2.5);
    REQUIRE(math::Hyper2F1(-2.5, 3.0, 3.0, z) / expected == Approx(1.0).epsilon(1.0e-13));
    REQUIRE(math::RegularizedHyper2F1(-2.5, 3.0, 3.0, z) / (expected / 2.0) == Approx(1.0).epsilon(1.0e-13));
    const double polynomial = std::pow(1.0 - z, 1.5) * (1.0 - 5.5 * z / 3.0);
    REQUIRE(math::Hyper2F1(-2.5, 4.0, 3.0, z) / polynomial == Approx(1.0).epsilon(1.0e-13));
    REQUIRE(math::RegularizedHyper2F1(-2.5, 4.0, 3.0, z) / (polynomial / 2.0) == Approx(1.0).epsilon(1.0e-13));
  }
  REQUIRE_THROWS_AS(math::Hyper2F1(1000.0, 1000.0, 1.0, 0.5), gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(math::RegularizedHyper2F1(1000.0, 1000.0, 1.0, 0.5), gra::AmplitudeFailure);
}

// Check all Gamma precisions and independent signs during concurrent evaluation
TEST_CASE("Log Gamma keeps signs local to concurrent calls", "[MMath][threading]") {
  const std::array<double, 8> x{-2.5, -1.5, -0.5, 0.5, 1.5, 5.0, 10.0, 20.0};
  std::array<double, 8> error{};
  // Compare logarithms and signs with the independent thread-safe tgamma function
  const auto difference = []<typename T>(const T value) {
    const T gamma = std::tgamma(value);
    int sign = 0;
    const T logarithm = math::LogGamma(value, &sign);
    if (sign != (gamma < T(0) ? -1 : 1)) { return std::numeric_limits<double>::infinity(); }
    return static_cast<double>(std::abs(logarithm - std::log(std::abs(gamma))) / std::numeric_limits<T>::epsilon());
  };
  math::ParallelFor(error.size(), [&](const std::size_t worker) {
    for (std::size_t iteration = 0; iteration < 1000; ++iteration) {
      const double value = x[(worker + iteration) % x.size()];
      error[worker] = std::max({error[worker], difference(static_cast<float>(value)), difference(value),
                                difference(static_cast<long double>(value))});
    }
  });
  for (const double value : error) { REQUIRE(value < 128.0); }
}

// Check matrix reshaping preserves storage, ordering and valid dimensions
TEST_CASE("Matrix reshape preserves row-major entries", "[MMatrix][reshape]") {
  gra::MMatrix<double> matrix{{1.0, 2.0, 3.0}, {4.0, 5.0, 6.0}};
  const auto          *storage = matrix[0];
  const auto           entries = matrix.Flatten();
  matrix.Reshape(3, 2);
  REQUIRE(matrix.size_row() == 3);
  REQUIRE(matrix.size_col() == 2);
  REQUIRE(matrix[0] == storage);
  REQUIRE(matrix.Flatten() == entries);
  REQUIRE(matrix[2][0] == Approx(5.0));
  REQUIRE_THROWS_AS(matrix.Reshape(2, 2), std::invalid_argument);
  REQUIRE_THROWS_AS(matrix.Reshape(std::numeric_limits<std::size_t>::max(), 2), std::invalid_argument);
  REQUIRE(matrix.size_row() == 3);
  REQUIRE(matrix.Flatten() == entries);
  gra::MMatrix<double> empty;
  empty.Reshape(0, 4);
  REQUIRE(empty.size_col() == 4);
}
