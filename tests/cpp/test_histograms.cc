// Unit tests for one-dimensional and two-dimensional histograms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "catch.hpp"

#include <complex>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Tech/MH1.h"
#include "Graniitti/Tech/MH2.h"

// Check physical variable-bin coordinates and safe public indexing
TEST_CASE("MH1 variable bins use their physical centers", "[MH1]") {
  gra::MH1<double> histogram({{10.0, 12.0}, {12.0, 16.0}});
  histogram.Fill(11.0, 1.0);
  histogram.Fill(14.0, 3.0);

  CHECK(histogram.GetBinXVal(0) == Approx(11.0));
  CHECK(histogram.GetBinXVal(1) == Approx(14.0));
  CHECK(histogram.GetMean(1) == Approx(13.25));
  CHECK_THROWS_AS(histogram.GetBinXVal(-1), std::out_of_range);
  CHECK_THROWS_AS(histogram.GetBinXVal(2), std::out_of_range);
}

// Check the documented closed upper boundary and complex-weight uncertainty
TEST_CASE("MH1 includes its upper boundary and computes norm errors", "[MH1]") {
  gra::MH1<double> closed(2, 0.0, 2.0);
  closed.Fill(2.0);
  CHECK(closed.GetBinCount(1) == 1);
  CHECK(closed.SumBinCounts() == 1);

  gra::MH1<std::complex<double>> complex_histogram(1, 0.0, 1.0);
  complex_histogram.Fill(0.5, {3.0, 4.0});
  CHECK(complex_histogram.GetBinError(0) == Approx(5.0));

  gra::MH1<double> signed_histogram(1, 0.0, 1.0);
  signed_histogram.Fill(0.5, -2.0);
  CHECK(signed_histogram.GetBinWeight(0) == Approx(-2.0));
  CHECK(signed_histogram.GetPositiveDefinite(0) == Approx(2.0));
}

// Check autorange for signed coordinates and degenerate input samples
TEST_CASE("MH1 autorange retains signed and repeated samples", "[MH1]") {
  SECTION("signed coordinates") {
    gra::MH1<double> histogram(4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(-4.0, 1.0);
    histogram.Fill(-2.0, 1.0);
    CHECK(histogram.FillCount() == 2);
    CHECK(histogram.SumBinCounts() == 2);
  }

  SECTION("repeated coordinates") {
    gra::MH1<double> histogram(4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(3.0);
    histogram.Fill(3.0);
    int bins = 0;
    double minimum = 0.0;
    double maximum = 0.0;
    histogram.GetBounds(bins, minimum, maximum);
    CHECK(maximum > minimum);
    CHECK(histogram.SumBinCounts() == 2);
  }

  SECTION("clear discards buffered samples") {
    gra::MH1<double> histogram(4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(-10.0);
    histogram.Clear();
    histogram.Fill(2.0);
    histogram.Fill(4.0);
    CHECK(histogram.FillCount() == 2);
    CHECK(histogram.SumBinCounts() == 2);
  }

  SECTION("empty buffers receive finite bounds") {
    gra::MH1<double> histogram(4);
    histogram.FlushBuffer();
    int bins = 0;
    double minimum = 0.0;
    double maximum = 0.0;
    histogram.GetBounds(bins, minimum, maximum);
    CHECK(bins == 4);
    CHECK(minimum == Approx(-0.5));
    CHECK(maximum == Approx(0.5));
  }
}

// Check that arithmetic cannot combine physically different binning
TEST_CASE("Histogram arithmetic requires identical binning", "[MH1][MH2]") {
  const gra::MH1<double> first(2, 0.0, 1.0);
  const gra::MH1<double> second(2, 1.0, 2.0);
  CHECK_THROWS_AS(first + second, std::domain_error);

  const gra::MH2 first2(2, 0.0, 1.0, 2, 0.0, 1.0);
  const gra::MH2 second2(2, 0.0, 2.0, 2, 0.0, 1.0);
  CHECK_THROWS_AS(first2 + second2, std::domain_error);
}

// Check uncertainty propagation for independent histogram subtraction
TEST_CASE("MH1 subtraction adds independent statistical variances", "[MH1]") {
  gra::MH1<double> first(1, 0.0, 1.0);
  gra::MH1<double> second(1, 0.0, 1.0);
  first.Fill(0.5, 3.0);
  second.Fill(0.5, 2.0);
  const auto difference = first - second;
  CHECK(difference.GetBinWeight(0) == Approx(1.0));
  CHECK(difference.GetBinError(0) == Approx(std::sqrt(13.0)));
  CHECK(difference.GetBinCount(0) == 2);
}

// Check independent product and ratio variances with unequal sample counts
TEST_CASE("MH1 products and ratios propagate statistical variances", "[MH1]") {
  gra::MH1<double> first(1, 0.0, 1.0);
  gra::MH1<double> second(1, 0.0, 1.0);
  first.Fill(0.5);
  second.Fill(0.5);
  second.Fill(0.5);
  const auto ratio = first / second;
  CHECK(ratio.GetBinWeight(0) == Approx(0.5));
  CHECK(ratio.GetBinError(0) == Approx(std::sqrt(0.375)));
  CHECK(ratio.GetBinCount(0) == first.GetBinCount(0));
  const auto product = first * second;
  CHECK(product.GetBinWeight(0) == Approx(2.0));
  CHECK(product.GetBinError(0) == Approx(std::sqrt(6.0)));
  CHECK(product.GetBinCount(0) == first.GetBinCount(0));

  gra::MH1<double> empty(1, 0.0, 1.0);
  CHECK((first / empty).GetBinWeight(0) == Approx(0.0));
  CHECK((first / empty).GetBinError(0) == Approx(0.0));

  gra::MH1<std::complex<double>> complex_first(1, 0.0, 1.0);
  gra::MH1<std::complex<double>> complex_second(1, 0.0, 1.0);
  complex_first.Fill(0.5, {0.0, 1.0});
  complex_second.Fill(0.5, {0.0, 1.0});
  complex_second.Fill(0.5, {0.0, 1.0});
  CHECK((complex_first / complex_second).GetBinError(0) == Approx(std::sqrt(0.375)));
  CHECK((complex_first * complex_second).GetBinError(0) == Approx(std::sqrt(6.0)));
  CHECK((complex_first / complex_second).GetdOdXError(0) == Approx(std::sqrt(0.375)));
  CHECK((complex_first * complex_second).GetdOdXError(0) == Approx(4.0 * std::sqrt(6.0)));
}

// Check ratio uncertainties when the squared bin sum exceeds double range
TEST_CASE("MH1 ratio variances retain extreme weight scales", "[MH1]") {
  for (const double scale : {7e153, 1e-150}) {
    gra::MH1<double> first(1, 0.0, 1.0);
    gra::MH1<double> second(1, 0.0, 1.0);
    first.Fill(0.5, scale);
    second.Fill(0.5, scale);
    second.Fill(0.5, scale);
    const auto ratio = first / second;
    CHECK(ratio.GetBinWeight(0) == Approx(0.5));
    CHECK(ratio.GetBinError(0) == Approx(std::sqrt(0.375)));
  }

  gra::MH1<double> large(1, 0.0, 1.0);
  gra::MH1<double> small(1, 0.0, 1.0);
  large.Fill(0.5, 7e153);
  large.Fill(0.5, 7e153);
  small.Fill(0.5, 1e-154);
  const auto product = large * small;
  CHECK(product.GetBinWeight(0) == Approx(1.4));
  CHECK(product.GetBinError(0) == Approx(std::sqrt(2.94)));
}

// Check two-dimensional upper boundaries and empty statistics
TEST_CASE("MH2 handles boundaries and empty statistics", "[MH2]") {
  gra::MH2 histogram(2, 0.0, 2.0, 2, -1.0, 1.0);
  const auto [mean, error] = histogram.WeightMeanAndError();
  CHECK(mean == Approx(0.0));
  CHECK(error == Approx(0.0));
  histogram.Fill(2.0, 1.0, 2.0);
  CHECK(histogram.GetBinCount(1, 1) == 1);
  CHECK(histogram.SumBinCounts() == 1);

  gra::MH2 automatic(2, 2);
  CHECK_NOTHROW(automatic.Clear());
  CHECK_THROWS_AS(automatic.SetLogX(), std::invalid_argument);

  automatic.FlushBuffer();
  int xbins = 0;
  int ybins = 0;
  double xmin = 0.0;
  double xmax = 0.0;
  double ymin = 0.0;
  double ymax = 0.0;
  automatic.GetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
  CHECK(xbins == 2);
  CHECK(ybins == 2);
  CHECK(xmin == Approx(-0.5));
  CHECK(xmax == Approx(0.5));
  CHECK(ymin == Approx(-0.5));
  CHECK(ymax == Approx(0.5));
}

// Check two-dimensional autorange with signed weights and repeated samples
TEST_CASE("MH2 autorange retains signed and repeated samples", "[MH2]") {
  SECTION("signed coordinates and weights") {
    gra::MH2 histogram(4, 4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(-4.0, -8.0, -1.0);
    histogram.Fill(-2.0, -4.0, -1.0);
    CHECK(histogram.FillCount() == 2);
    CHECK(histogram.SumBinCounts() == 2);
  }

  SECTION("repeated coordinates") {
    gra::MH2 histogram(4, 4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(3.0, 5.0);
    histogram.Fill(3.0, 5.0);
    int xbins = 0;
    int ybins = 0;
    double xmin = 0.0;
    double xmax = 0.0;
    double ymin = 0.0;
    double ymax = 0.0;
    histogram.GetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
    CHECK(xmax > xmin);
    CHECK(ymax > ymin);
    CHECK(histogram.SumBinCounts() == 2);
  }

  SECTION("clear discards buffered samples") {
    gra::MH2 histogram(4, 4);
    histogram.SetAutoBuffSize(2);
    histogram.Fill(-10.0, -20.0);
    histogram.Clear();
    histogram.Fill(2.0, 4.0);
    histogram.Fill(4.0, 8.0);
    CHECK(histogram.FillCount() == 2);
    CHECK(histogram.SumBinCounts() == 2);
  }
}

// Check log bins against the same physical variable bins and integrated density
TEST_CASE("Logarithmic histograms use physical widths and moments", "[MH1][MH2][normalization]") {
  gra::MH1<double> logarithmic(2, 1.0, 100.0);
  logarithmic.SetLogX();
  gra::MH1<double> variable({{1.0, 10.0}, {10.0, 100.0}});
  for (auto *histogram : {&logarithmic, &variable}) {
    histogram->Fill(5.0, 1.0);
    histogram->Fill(50.0, 3.0);
  }
  CHECK(logarithmic.GetMean(1) == Approx(variable.GetMean(1)));
  CHECK(logarithmic.GetMean(2) == Approx(variable.GetMean(2)));
  for (int bin = 0; bin < 2; ++bin) {
    CHECK(logarithmic.GetdOdX(bin) == Approx(variable.GetdOdX(bin)));
    CHECK(logarithmic.GetdOdXError(bin) == Approx(variable.GetdOdXError(bin)));
  }
  CHECK(9.0 * logarithmic.GetdOdX(0) + 90.0 * logarithmic.GetdOdX(1) == Approx(2.0));

  gra::MH2 both(2, 1.0, 100.0, 2, 1.0, 100.0);
  both.SetLogXY();
  both.Fill(5.0, 50.0, 1.0);
  both.Fill(50.0, 5.0, 3.0);
  CHECK(both.GetMeanX(1) == Approx(variable.GetMean(1)));
  CHECK(both.GetMeanX(2) == Approx(variable.GetMean(2)));
  CHECK(both.GetMeanY(1) == Approx((55.0 + 3.0 * 5.5) / 4.0));
  CHECK(both.GetMeanY(2) == Approx((55.0 * 55.0 + 3.0 * 5.5 * 5.5) / 4.0));
}

// Check logarithmic bounds cannot be reset into an undefined real domain
TEST_CASE("Histogram resets retain valid logarithmic axes", "[MH1][MH2][regression]") {
  gra::MH1<double> first(2, 1.0, 100.0);
  first.SetLogX();
  first.Fill(5.0);
  CHECK_THROWS_AS(first.ResetBounds(2, 0.0, 100.0), std::invalid_argument);
  CHECK(first.GetMean(1) == Approx(5.5));
  CHECK(first.FillCount() == 1);
  first.ResetBounds(2, 10.0, 1000.0);
  first.Fill(50.0);
  CHECK(first.GetMean(1) == Approx(55.0));

  SECTION("autorange resets clear logarithmic binning") {
    first.ResetBounds(2);
    first.SetAutoBuffSize(2);
    first.Fill(-4.0);
    first.Fill(-2.0);
    CHECK(first.SumBinCounts() == 2);
    CHECK(std::isfinite(first.GetMean(1)));
    CHECK(first.GetMean(1) < 0.0);
  }
  SECTION("explicit variable bins replace logarithmic binning") {
    first.ResetBounds({{-4.0, -2.0}, {-2.0, 0.0}});
    first.Fill(-3.0);
    CHECK(first.GetMean(1) == Approx(-3.0));
  }

  gra::MH2 second(2, 1.0, 100.0, 2, 1.0, 100.0);
  second.SetLogXY();
  second.Fill(5.0, 50.0);
  CHECK_THROWS_AS(second.ResetBounds(2, -1.0, 100.0, 2, 1.0, 100.0), std::invalid_argument);
  CHECK_THROWS_AS(second.ResetBounds(2, 1.0, 100.0, 2, 0.0, 100.0), std::invalid_argument);
  CHECK(second.FillCount() == 1);
  CHECK(second.GetMeanX(1) == Approx(5.5));
  CHECK(second.GetMeanY(1) == Approx(55.0));
}

// Check a failed binning change preserves the stored physical coordinates
TEST_CASE("Logarithmic binning changes are atomic and precede filling", "[MH1][MH2][regression]") {
  gra::MH1<double> first(2, 1.0, 100.0);
  first.Fill(5.0);
  const double mean = first.GetMean(1);
  CHECK_THROWS_AS(first.SetLogX(), std::logic_error);
  CHECK(first.GetMean(1) == Approx(mean));

  gra::MH2 mixed(2, 1.0, 100.0, 2, -1.0, 1.0);
  CHECK_THROWS_AS(mixed.SetLogXY(), std::invalid_argument);
  mixed.Fill(5.0, 0.5);
  CHECK(mixed.GetMeanX(1) == Approx(mean));
  CHECK_THROWS_AS(mixed.SetLogX(), std::logic_error);

  gra::MH2 positive(2, 1.0, 100.0, 2, 1.0, 100.0);
  positive.Fill(5.0, 5.0);
  CHECK_THROWS_AS(positive.SetLogXY(), std::logic_error);
  CHECK_THROWS_AS(positive.SetLogY(), std::logic_error);
  CHECK(positive.GetMeanX(1) == Approx(mean));
  CHECK(positive.GetMeanY(1) == Approx(mean));
  positive.Clear();
  CHECK_NOTHROW(positive.SetLogXY());
  positive.Fill(5.0, 5.0);
  CHECK_NOTHROW(positive.SetLogXY());
  CHECK(positive.GetMeanX(1) == Approx(5.5));
}

// Check roundoff near the upper bound cannot discard an in-range event
TEST_CASE("Histogram indices retain points adjacent to the upper bound", "[MH1][MH2][normalization]") {
  for (const bool logarithmic : {false, true}) {
    const double low = logarithmic ? 1.0 : -1.0;
    const double high = logarithmic ? 100.0 : 1.0;
    gra::MH1<double> first(50, low, high);
    gra::MH2 second(50, low, high, 50, low, high);
    if (logarithmic) { first.SetLogX(); second.SetLogXY(); }
    const double x = std::nextafter(high, low);
    first.Fill(x);
    second.Fill(x, x);
    CHECK(first.GetBinCount(49) == 1);
    CHECK(second.GetBinCount(49, 49) == 1);
    CHECK(first.SumBinCounts() == first.FillCount());
    CHECK(second.SumBinCounts() == second.FillCount());
  }
}

// Check signed-bin statistics and entropy use consistent absolute probabilities
TEST_CASE("Histogram extrema and entropy handle signed and large weights", "[MH1][MH2][statistics]") {
  gra::MH1<double> first(2, 0.0, 2.0);
  gra::MH2 second(2, 0.0, 2.0, 1, 0.0, 1.0);
  first.Fill(0.5, -2.0);
  first.Fill(1.5, 1.0);
  second.Fill(0.5, 0.5, -2.0);
  second.Fill(1.5, 0.5, 1.0);
  const auto probability = first.GetProbDensity();
  const double entropy = -probability[0] * std::log2(probability[0]) - probability[1] * std::log2(probability[1]);
  CHECK(second.ShannonEntropy() == Approx(entropy));
  second.Clear();
  second.Fill(0.5, 0.5, -2.0);
  second.Fill(1.5, 0.5, -1.0);
  CHECK(second.GetMaxWeight() == Approx(-1.0));
  CHECK(second.GetMinWeight() == Approx(-2.0));
  first.Clear();
  second.Clear();
  for (const double x : {0.5, 1.5}) {
    first.Fill(x, 2.0e130);
    second.Fill(x, 0.5, 2.0e130);
  }
  CHECK(first.GetMinWeight() == Approx(2.0e130));
  CHECK(second.GetMinWeight() == Approx(2.0e130));
}

// Check signed differential weights integrate to the mean Monte Carlo weight
TEST_CASE("Signed MH1 distributions retain cancellations and moments", "[MH1][normalization]") {
  gra::MH1<double> first({{0.0, 1.0}, {1.0, 3.0}});
  gra::MH2 second(2, 0.0, 2.0, 1, 0.0, 1.0);
  first.Fill(0.5, -2.0);
  first.Fill(2.0, 3.0);
  first.Fill(4.0, 5.0);
  CHECK(first.GetdOdX(0) == Approx(-2.0 / 3.0));
  CHECK(first.GetdOdX(0) + 2.0 * first.GetdOdX(1) == Approx(first.WeightMeanAndError().first));
  CHECK(first.GetMean(1) == Approx(5.0));
  CHECK(first.GetdOdXError(0) == Approx(2.0 / 3.0));
  std::vector<std::vector<double>> bins, values;
  first.VectorOutput(bins, values);
  CHECK(bins == std::vector<std::vector<double>>{{0.0, 1.0}, {1.0, 3.0}});
  CHECK(values[0][0] == Approx(-2.0));
  CHECK(values[0][1] == Approx(2.0));

  first.ResetBounds(2, 0.0, 2.0);
  for (const auto &[x, weight] : std::vector<std::pair<double, double>>{{0.5, -2.0}, {1.5, 3.0}}) {
    first.Fill(x, weight);
    second.Fill(x, 0.5, weight);
  }
  CHECK(first.GetMean(1) == Approx(second.GetMeanX(1)));
  CHECK(first.GetMean(2) == Approx(second.GetMeanX(2)));
}

// Check shared bin boundaries retain the same events after explicit-edge conversion
TEST_CASE("MH1 bin representations agree on shared boundaries", "[MH1][normalization]") {
  for (const bool logarithmic : {false, true}) {
    gra::MH1<double> first(2, 1.0, logarithmic ? 100.0 : 3.0);
    if (logarithmic) { first.SetLogX(); }
    gra::MH1<double> second({{first.GetBinXVal(0, -1), first.GetBinXVal(0, 1)},
                             {first.GetBinXVal(1, -1), first.GetBinXVal(1, 1)}});
    const double boundary = first.GetBinXVal(0, 1);
    for (const double x : {1.0, std::nextafter(boundary, 1.0), boundary,
                           std::nextafter(boundary, 100.0), first.GetBinXVal(1, 1)}) {
      first.Fill(x);
      second.Fill(x);
    }
    CHECK(first.GetCounts() == second.GetCounts());
    CHECK(first.SumBinCounts() == 5);
  }
}

// Check symmetric autorange preserves both beam orientations without losing support
TEST_CASE("Histogram symmetric autorange is invariant under reflection", "[MH1][MH2][symmetry]") {
  for (const double sign : {-1.0, 1.0}) {
    gra::MH1<double> first(8);
    gra::MH2 second(8, 8);
    first.SetAutoSymmetry(true);
    second.SetAutoSymmetry({true, true});
    for (const double x : {2.0, 4.0}) {
      first.Fill(sign * x);
      second.Fill(sign * x, -sign * x);
    }
    first.FlushBuffer();
    second.FlushBuffer();
    int nx = 0, ny = 0;
    double xmin = 0.0, xmax = 0.0, ymin = 0.0, ymax = 0.0;
    first.GetBounds(nx, xmin, xmax);
    CHECK(xmin == Approx(-xmax));
    CHECK(first.SumBinCounts() == 2);
    second.GetBounds(nx, xmin, xmax, ny, ymin, ymax);
    CHECK(xmin == Approx(-xmax));
    CHECK(ymin == Approx(-ymax));
    CHECK(second.SumBinCounts() == 2);
  }
}

// Check coherent densities retain their normalization under bin splitting and repeated sampling
TEST_CASE("Complex MH1 densities normalize amplitudes before squaring", "[MH1][normalization][phase]") {
  for (const auto phase : {std::complex<double>(1.0, 0.0), std::polar(1.0, 0.73)}) {
    gra::MH1<std::complex<double>> coarse(1, 0.0, 2.0);
    gra::MH1<std::complex<double>> fine(2, 0.0, 2.0);
    for (int repetition = 0; repetition < 2; ++repetition) {
      for (const double x : {0.25, 0.75, 1.25, 1.75}) {
        const auto weight = phase * std::complex<double>(2.0, 2.0);
        coarse.Fill(x, weight);
        fine.Fill(x, weight);
      }
      CHECK(coarse.GetdOdX(0) == Approx(2.0));
      CHECK(fine.GetdOdX(0) == Approx(2.0));
      CHECK(fine.GetdOdX(1) == Approx(2.0));
    }
  }
  gra::MH1<std::complex<double>> variable({{0.0, 1.0}, {1.0, 3.0}});
  for (const double x : {0.5, 1.5, 2.5}) { variable.Fill(x, {3.0, 0.0}); }
  const auto probability = variable.GetProbDensity();
  CHECK(probability[0] == Approx(1.0 / 3.0));
  CHECK(probability[1] == Approx(2.0 / 3.0));
  CHECK(variable.GetMean(1) == Approx(1.5));
  CHECK(variable.GetdOdX(0) == Approx(1.0));
  CHECK(variable.GetdOdX(1) == Approx(1.0));
}

// Check coherent-norm errors against event-by-event real projections and global phase rotation
TEST_CASE("Complex MH1 errors include both amplitude quadratures", "[MH1][statistics][phase]") {
  for (const auto phase : {std::complex<double>(1.0, 0.0), std::polar(1.0, 0.73)}) {
    gra::MH1<std::complex<double>> histogram(1, 0.0, 2.0);
    const std::vector<std::complex<double>> weights{{3.0, 4.0}, {1.0, -2.0}};
    const auto amplitude = phase * (weights[0] + weights[1]);
    double variance = 0.0;
    for (const auto weight : weights) {
      histogram.Fill(0.5, phase * weight);
      variance += std::pow(2.0 * std::real(std::conj(amplitude) * phase * weight), 2);
    }
    std::vector<std::vector<double>> bins, values;
    histogram.VectorOutput(bins, values);
    CHECK(values[0][0] == Approx(20.0));
    CHECK(values[0][1] == Approx(std::sqrt(variance)));
    CHECK(histogram.GetBinError(0) == Approx(std::sqrt(30.0)));
    CHECK(histogram.GetdOdXError(0) == Approx(std::sqrt(variance) / 16.0));
    const auto sum = histogram + histogram;
    CHECK(sum.GetdOdX(0) == Approx(histogram.GetdOdX(0)));
    CHECK(sum.GetdOdXError(0) == Approx(histogram.GetdOdXError(0) / std::sqrt(2.0)));
  }
}
