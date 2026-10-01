// Test mathematical limits and angular momentum identities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MPolarFourier.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Tech/MAux.h"
#include "catch.hpp"

using gra::aux::indices;

// Check spherical harmonic normalization, parity and azimuthal rotation at large angular momentum
// [REFERENCE: NIST DLMF 14.30.6, 14.30.7 and 14.30.9, https://dlmf.nist.gov/14.30]
TEST_CASE("Spherical harmonics retain normalization beyond the factorial range", "[MMath][spherical][physics]") {
  for (const int l : {0, 1, 4, 16, 85, 86, 100, 180, 256}) {
    for (const double x : {-0.83, 0.0, 0.42, 1.0}) {
      CAPTURE(l, x);
      const double phi   = 0.37;
      const double alpha = 0.43;
      double       norm  = 0.0;
      for (int m = -l; m <= l; ++m) {
        CAPTURE(m);
        const auto y         = gra::math::Y_complex_basis(x, phi, l, m);
        const auto conjugate = gra::math::Y_complex_basis(x, phi, l, -m);
        const auto parity    = gra::math::Y_complex_basis(-x, phi + gra::math::PI, l, m);
        const auto rotated   = gra::math::Y_complex_basis(x, phi + alpha, l, m);
        REQUIRE(std::isfinite(y.real()));
        REQUIRE(std::isfinite(y.imag()));
        CHECK(std::abs(conjugate - (std::abs(m) % 2 == 0 ? 1.0 : -1.0) * std::conj(y)) < 2.0e-12);
        CHECK(std::abs(parity - (l % 2 == 0 ? 1.0 : -1.0) * y) < 2.0e-12);
        CHECK(std::abs(rotated - std::exp(gra::math::zi * (m * alpha)) * y) < 2.0e-12);
        CHECK(gra::math::NReY(y, m) == Approx(gra::math::Y_real_basis(x, phi, l, m)).margin(2.0e-13));
        norm += std::norm(y);
      }
      CHECK(norm == Approx((2.0 * l + 1.0) / (4.0 * gra::math::PI)).epsilon(2.0e-12));
    }
  }
}

// Check rotational invariance of the addition theorem and integrated partial wave norms
// [REFERENCE: NIST DLMF 14.30.8 and 14.30.9, https://dlmf.nist.gov/14.30]
TEST_CASE("Spherical harmonics obey the addition theorem and angular integral", "[MMath][spherical][physics]") {
  const double x      = 0.37;
  const double y      = -0.42;
  const double phi    = 0.29;
  const double psi    = 1.14;
  const double cosine = x * y + std::sqrt((1.0 - x * x) * (1.0 - y * y)) * std::cos(phi - psi);
  for (const int l : {4, 32, 86, 100, 180}) {
    std::complex<double> sum{};
    for (int m = -l; m <= l; ++m) {
      sum += std::conj(gra::math::Y_complex_basis(x, phi, l, m)) * gra::math::Y_complex_basis(y, psi, l, m);
    }
    const double expected = (2.0 * l + 1.0) / (4.0 * gra::math::PI) * gra::math::sf_legendre(l, 0, cosine);
    CHECK(sum.real() == Approx(expected).margin(2.0e-12));
    CHECK(std::abs(sum.imag()) < 2.0e-12);
  }
  const auto [node, weight] = gra::math::GaussLegendreRule(256, -1.0, 1.0);
  for (const std::array<int, 2> lm : {std::array{86, 86}, {100, 70}, {180, 0}, {180, 179}}) {
    CAPTURE(lm[0], lm[1]);
    double integral = 0.0;
    for (const auto &i : indices(node)) {
      integral += 2.0 * gra::math::PI * weight[i] * std::norm(gra::math::Y_complex_basis(node[i], phi, lm[0], lm[1]));
    }
    CHECK(integral == Approx(1.0).epsilon(2.0e-11));
  }
}

// Check the small-argument limit and parity without forming inverse powers
// [REFERENCE: NIST DLMF 10.2.2, https://dlmf.nist.gov/10.2.E2]
TEST_CASE("Bessel functions preserve representable small-argument limits", "[MMath][Bessel]") {
  for (const double x : {std::numeric_limits<double>::denorm_min(), 1.0e-310, 1.0e-100, 1.0e-10}) {
    CHECK(gra::math::LogBesselI1(x) == Approx(std::log(x) - std::log(2.0)).margin(2.0e-13));
  }
  for (const double x : {1.0e-10, 1.0e-50, 1.0e-100, 1.0e-150, 1.0e-310}) {
    for (const int n : {2, 3, 5, 10}) {
      CAPTURE(x, n);
      const long double leading  = std::pow(static_cast<long double>(x) / 2.0L, n) / std::tgamma(n + 1.0L);
      const double      expected = static_cast<double>(leading);
      const double      value    = gra::math::BesselJ(n, x);
      REQUIRE(std::isfinite(value));
      if (expected > 0.0) {
        CHECK(value / expected == Approx(1.0).epsilon(3.0e-14));
      } else {
        CHECK(gra::math::IsZero(value));
      }
      CHECK(gra::math::BesselJ(n, -x) == Approx((n % 2 == 0 ? 1.0 : -1.0) * value).epsilon(3.0e-14));
    }
  }
  for (const int n : {2, 3, 8, 20}) {
    for (const double x : {0.1, 0.7, 1.0}) {
      CHECK(gra::math::BesselJ(n, x) == Approx(std::cyl_bessel_j(static_cast<double>(n), x)).epsilon(3.0e-13));
    }
  }
}

// Check scale independent normalization of real and complex state vectors
TEST_CASE("Vector normalization preserves finite directions at numerical limits", "[MMatrix][normalization]") {
  for (const double scale : {1.0e-310, 1.0e-300, 1.0, 1.0e300, 1.0e307}) {
    const auto real = gra::NormalizedL2(std::array{3.0 * scale, 4.0 * scale});
    CHECK(real[0] == Approx(0.6).epsilon(3.0e-13));
    CHECK(real[1] == Approx(0.8).epsilon(3.0e-13));
    const auto probability = gra::NormalizedSum(std::array{3.0 * scale, 4.0 * scale});
    CHECK(probability[0] == Approx(3.0 / 7.0).epsilon(3.0e-13));
    CHECK(probability[1] == Approx(4.0 / 7.0).epsilon(3.0e-13));
    const auto complex = gra::NormalizedL2(std::vector<std::complex<double>>{{3.0 * scale, 4.0 * scale}, {0.0, 0.0}});
    CHECK(complex[0].real() == Approx(0.6).epsilon(3.0e-13));
    CHECK(complex[0].imag() == Approx(0.8).epsilon(3.0e-13));
  }
  const double maximum     = std::numeric_limits<double>::max();
  const auto   probability = gra::NormalizedSum(std::array{maximum, maximum});
  CHECK(probability[0] == Approx(0.5).epsilon(1.0e-14));
  const auto complex = gra::NormalizedL2(std::array<std::complex<double>, 1>{{{maximum, maximum}}});
  CHECK(complex[0].real() == Approx(1.0 / std::sqrt(2.0)).epsilon(1.0e-14));
  const auto single = gra::NormalizedL2(std::array<std::complex<float>, 1>{{{3.0F, 4.0F}}});
  CHECK(single[0].real() == Approx(0.6).epsilon(1.0e-6));
  CHECK(single[0].imag() == Approx(0.8).epsilon(1.0e-6));
}

// Check analytic Gauss identities near poles and near terminating parameters
TEST_CASE("Gauss functions preserve continuous parameters near integers", "[MMath][hypergeometric]") {
  for (const double c : {-1.0e-13, 1.0e-13}) {
    CHECK(gra::math::Hyper2F1(1.0, c + 1.0, c, 0.5) == Approx(2.0 + 2.0 / c).epsilon(2.0e-13));
  }
  const double a        = -1.0e-13;
  const double z        = 0.75;
  const double expected = std::pow(1.0 - z, -a - 1.0) * (1.0 + (a - 1.0) * z);
  CHECK(gra::math::Hyper2F1(a, 2.0, 1.0, z) == Approx(expected).epsilon(1.0e-14));
  CHECK(gra::math::RegularizedHyper2F1(a, 2.0, 1.0, z) == Approx(expected).epsilon(1.0e-14));
}

// Keep true polynomial zeros distinct from divergent or unresolved special functions
TEST_CASE("Special function failures cannot become physical zeros", "[MMath][hypergeometric]") {
  const double nan = std::numeric_limits<double>::quiet_NaN();
  CHECK(std::isnan(gra::math::msqrt(nan)));
  CHECK(gra::math::msqrt(-1.0e-15) == Approx(0.0));
  CHECK_THROWS_AS(gra::math::ReciprocalGamma(nan), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::math::Hyper2F1(1.0, 1.0, 2.0, 1.0), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::math::RegularizedHyper2F1(1.0, 1.0, 2.0, 1.0), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::math::RegularizedHyper3F2Unit(nan, 1.0, 1.0, 2.0, 2.0), gra::AmplitudeFailure);
  CHECK_THROWS_AS(gra::math::RegularizedHyper3F2Unit(1.0, 1.0, 1.0, 1.0, 1.0), gra::AmplitudeFailure);
  CHECK(gra::math::IsZero(gra::math::RegularizedHyper3F2Unit(-1.0, 1.0, 1.0, -2.0, 1.0)));
  CHECK(gra::math::RegularizedHyper3F2Unit(-1.0, 1.0, 1.0, 0.0, 1.0) == Approx(-1.0));
}

// Check slow convergent sums with exact zeta and Gauss values
// [REFERENCE: NIST DLMF 15.4.20, https://dlmf.nist.gov/15.4.E20]
TEST_CASE("Regularized 3F2 accelerates small convergence gaps", "[MMath][hypergeometric]") {
  CHECK(gra::math::RegularizedHyper3F2Unit(1.0, 1.0, 1.0, 2.0, 2.0) == Approx(gra::math::PIPI / 6.0).epsilon(2.0e-14));
  for (const double gap : {0.015625, 0.125, 0.5, 1.5}) {
    for (const double a : {-1.375, 0.375, 1.125}) {
      const double      b = 0.75;
      const double      c = 0.6875;
      const double      d = a + c + gap;
      const long double expected =
          std::tgamma(static_cast<long double>(gap)) /
          (std::tgamma(static_cast<long double>(b)) * std::tgamma(static_cast<long double>(d - a)) *
           std::tgamma(static_cast<long double>(d - c)));
      CAPTURE(a, b, c, d, gap);
      CHECK(gra::math::RegularizedHyper3F2Unit(a, b, c, b, d) ==
            Approx(static_cast<double>(expected)).epsilon(2.0e-13));
    }
  }
}

// Check generic noninteger parameter differences against Dixon's exact sum
// [REFERENCE: NIST DLMF 16.4.4, https://dlmf.nist.gov/16.4.E4]
TEST_CASE("Regularized 3F2 reproduces Dixon sums near the convergence boundary", "[MMath][hypergeometric]") {
  for (const long double a : {-0.375L, 0.375L, 1.125L, 2.25L}) {
    for (const long double b : {0.3125L, 0.6875L}) {
      for (const long double gap : {0.015625L, 0.125L, 0.5L, 1.5L, 4.0L}) {
        const long double c = 1.0L + a / 2.0L - b - gap / 2.0L;
        const long double expected =
            std::tgamma(1.0L + a / 2.0L) * std::tgamma(gap / 2.0L) /
            (std::tgamma(1.0L + a) * std::tgamma(1.0L + a / 2.0L - b) *
             std::tgamma(1.0L + a / 2.0L - c) * std::tgamma(1.0L + a - b - c));
        CAPTURE(a, b, c, gap);
        CHECK(gra::math::RegularizedHyper3F2Unit(a, b, c, 1.0L + a - b, 1.0L + a - c) ==
              Approx(static_cast<double>(expected)).epsilon(2.0e-13));
      }
    }
  }
}

// Compare with independent finite beta sums evaluated at 100 decimal digits
// Coincident poles use beta derivatives, and dyadic inputs preserve exact integer differences
// [REFERENCE: M A Shpot and H M Srivastava, Appl Math Comput 259 (2015) 819, https://arxiv.org/abs/1411.2455]
TEST_CASE("Regularized 3F2 matches continued GP coefficients at coincident trajectories",
          "[MMath][hypergeometric][physics]") {
  const std::array<std::array<double, 6>, 23> points{{
      {-1.96875, -3.0625, -2.90625, -1.90625, -2.0625, 0.2169816848755289},
      {-1.96875, -1.0625, -0.90625, 0.09375, -0.0625, -1.9217541902798971},
      {0.03125, 0.9375, 0.09375, 4.09375, 2.9375, 0.0783346932786244},
      {2.03125, -3.0625, -0.90625, 2.09375, 3.9375, 0.28866722824721014},
      {-2.125, -3.0625, -3.0625, -2.0625, -2.0625, -1.96722761576142},
      {-2.125, -1.0625, -1.0625, -0.0625, -0.0625, -2.208094512887298},
      {-0.125, 0.9375, -0.0625, 3.9375, 2.9375, 0.09542662416013538},
      {1.875, -3.0625, -1.0625, 1.9375, 3.9375, 0.3343145538955419},
      {-2.1250000009313226, -3.0625, -3.0625000009313226, -2.0625000009313226, -2.0625, -1.9672276341655741},
      {-2.1250000009313226, -1.0625, -1.0625000009313226, -0.06250000093132257, -0.0625, -2.2080945142415382},
      {-0.12500000093132257, 0.9375, -0.06250000093132257, 3.9374999990686774, 2.9375, 0.09542662427160252},
      {1.8749999990686774, -3.0625, -1.0625000009313226, 1.9374999990686774, 3.9375, 0.3343145541697444},
      {-0.40625, -1.3125, -3.09375, -2.09375, -0.3125, 0.03215003516271602},
      {-0.40625, 0.6875, -1.09375, -0.09375, 1.6875, 0.08577599832963172},
      {1.59375, 2.6875, -0.09375, 3.90625, 4.6875, 0.012022571034615665},
      {3.59375, -1.3125, -1.09375, 1.90625, 5.6875, 0.02157510450390786},
      {1.625, -0.6875, -1.6875, -0.6875, 0.3125, -0.1769077313746816},
      {1.625, 1.3125, 0.3125, 1.3125, 2.3125, 2.0075887590253045},
      {3.625, 3.3125, 1.3125, 5.3125, 5.3125, 0.0017281672094344304},
      {5.625, -0.6875, 0.3125, 3.3125, 6.3125, 0.001674451716941692},
      // This sampled GP point cancels initial terms more than 1000 times larger than the sum
      {-0.16972876305578044, -2.0877717302113643, -2.0819570328444161,
       0.91804296715558387, 0.91222826978863569, 0.00147938284995920096},
      // This sampled GP point requires more precision under the roundoff bound
      {-0.16992487261318701, -2.0816025307909674, -2.0883223418222197,
       0.91167765817778035, 0.91839746920903265, 0.000162978227198925524},
      // Resolve a true fusion node with cancellation of twelve decimal digits
      {-0.16995052861057047, -2.0849609375, -2.0849895911105705,
       0.9150104088894295, 0.9150390625, 1.2553947016262305e-12},
  }};
  for (const auto& point : points) {
    std::array<double, 3> upper{point[0], point[1], point[2]};
    std::sort(upper.begin(), upper.end());
    do {
      CAPTURE(upper, point[3], point[4]);
      CHECK(gra::math::RegularizedHyper3F2Unit(upper[0], upper[1], upper[2], point[3], point[4]) ==
            Approx(point[5]).epsilon(2.0e-14));
      CHECK(gra::math::RegularizedHyper3F2Unit(upper[0], upper[1], upper[2], point[4], point[3]) ==
            Approx(point[5]).epsilon(2.0e-14));
    } while (std::next_permutation(upper.begin(), upper.end()));
  }
}

// Check a small regularized 3F2 through an independently summed Gauss identity
TEST_CASE("Regularized 3F2 bounds its full tail relative to the result", "[MMath][hypergeometric]") {
  for (const double b : {2.25, 3.5, 10.0, 20.0, 100.0, 172.0}) {
    CAPTURE(b);
    const long double expected = (b - 1.0L) / (b - 2.0L) / std::tgamma(static_cast<long double>(b));
    const double      value    = gra::math::RegularizedHyper3F2Unit(1.0, 1.0, 1.0, b, 1.0);
    REQUIRE(value > 0.0);
    CHECK(static_cast<double>(static_cast<long double>(value) / expected) == Approx(1.0).epsilon(3.0e-13));
  }
}

// Check the GP scalar fusion series against partial fractions and the Euler beta integral
TEST_CASE("Regularized 3F2 resolves scalar fusion across trajectory arguments", "[MMath][hypergeometric][physics]") {
  for (const long double a : {-0.1L, -1.4L, -2.17L, -3.05L}) {
    for (const long double b : {-3.1L, -1.0L - 1.0e-9L, -1.0L + 1.0e-9L, 0.1L, 1.3L}) {
      for (const long double delta : {0.0001L, 0.17L, 1.4L}) {
        const long double c = b + delta;
        CAPTURE(a, b, c);
        const long double expected = std::tgamma(1.0L - a) / (c - b) *
            (1.0L / (std::tgamma(c) * std::tgamma(b + 1.0L - a)) -
             1.0L / (std::tgamma(b) * std::tgamma(c + 1.0L - a)));
        const double value = gra::math::RegularizedHyper3F2Unit(a, b, c, b + 1.0L, c + 1.0L);
        CHECK(value == Approx(static_cast<double>(expected)).epsilon(3.0e-12));
      }
    }
  }
}

// Check the scalar fusion tail bound against its independent beta integral
TEST_CASE("Regularized 3F2 resolves the continued scalar fusion tail", "[MMath][hypergeometric][physics]") {
  const long double a = -1.4L;
  const long double b = -2.3L;
  const long double c = -3.1L;
  const long double expected = std::tgamma(1.0L - a) / (c - b) *
      (1.0L / (std::tgamma(c) * std::tgamma(b + 1.0L - a)) -
       1.0L / (std::tgamma(b) * std::tgamma(c + 1.0L - a)));
  const double value = gra::math::RegularizedHyper3F2Unit(a, b, c, b + 1.0L, c + 1.0L);
  CHECK(value == Approx(static_cast<double>(expected)).epsilon(3.0e-13));
  long double term = 1.0L / (std::tgamma(b + 1.0L) * std::tgamma(c + 1.0L));
  long double sum = term;
  for (int k = 0; k < 64; ++k) {
    term *= (a + k) * (b + k) * (c + k) / ((b + 1.0L + k) * (c + 1.0L + k) * (k + 1.0L));
    sum += term;
    if (k + 1 == 16 || k + 1 == 64) {
      const auto [tail, error] = gra::math::Hyper3F2Remainder(
          {static_cast<double>(a), static_cast<double>(b), static_cast<double>(c)},
          {static_cast<double>(b + 1.0L), static_cast<double>(c + 1.0L)}, k + 1, term);
      CHECK(std::abs(expected - sum - tail) <= error + 1.0e-15L * std::abs(expected));
      CHECK(error < std::abs(expected - sum));
    }
  }
}

// Keep nonfinite screening fields inside the recoverable amplitude failure path
TEST_CASE("Polar Fourier numerical failures are recoverable during sampling", "[MPolarFourier][physics]") {
  const gra::math::MPolarFourier transform({{1.0}, {1.0}}, {{1.0}, {1.0}}, {0.0}, 0);
  for (const double value : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    const gra::MMatrix<std::complex<double>> matrix(1, 1, std::complex<double>{value, 0.0});
    CHECK_THROWS_AS(transform.Forward(matrix), gra::AmplitudeFailure);
    const gra::math::PolarHarmonicField field{0, matrix};
    CHECK_THROWS_AS(transform.InverseEvaluate(field, 0.0, 0.0), gra::AmplitudeFailure);
  }
}
