// 2D histogram class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <iostream>
#include <limits>
#include <numeric>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH2.h"

// Libraries
#include "rang.hpp"

namespace gra {
// Constructor
MH2::MH2(int xbins, double xmin, double xmax, int ybins, double ymin,
         double ymax, std::string namestr) {
  name = namestr;
  ResetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
  FILLBUFF = false;
}

// Constructor with only number of bins
MH2::MH2(int xbins, int ybins, std::string namestr) {
  if (xbins <= 0 || ybins <= 0) {
    throw std::invalid_argument("MH2: autorange bin counts must be positive");
  }
  name = namestr;
  XBINS = xbins;
  YBINS = ybins;
  weights = MMatrix<double>(XBINS, YBINS, 0.0);
  weights2 = weights;
  counts = MMatrix<long long int>(XBINS, YBINS, 0);
  FILLBUFF = true;
}

// Empty constructor
MH2::MH2() {
  XBINS = 50; // Default
  YBINS = 50;
  weights = MMatrix<double>(XBINS, YBINS, 0.0);
  weights2 = weights;
  counts = MMatrix<long long int>(XBINS, YBINS, 0);
  FILLBUFF = true;
}

// Destructor
MH2::~MH2() = default;

// Reset histogram storage to explicit finite bounds
void MH2::ResetBounds(int xbins, double xmin, double xmax, int ybins,
                      double ymin, double ymax) {
  if (xbins <= 0 || ybins <= 0 || !std::isfinite(xmin) ||
      !std::isfinite(xmax) || !std::isfinite(ymin) || !std::isfinite(ymax) ||
      !(xmax > xmin) || !(ymax > ymin) ||
      (LOGX && !(xmin > 0.0)) || (LOGY && !(ymin > 0.0))) {
    throw std::invalid_argument(
        "MH2::ResetBounds: require positive bins, finite ordered bounds and positive logarithmic bounds");
  }
  XMIN = xmin;
  XMAX = xmax;
  XBINS = xbins;

  YMIN = ymin;
  YMAX = ymax;
  YBINS = ybins;

  // Init
  weights = MMatrix<double>(XBINS, YBINS, 0.0);
  weights2 = weights;
  counts = MMatrix<long long int>(XBINS, YBINS, 0);

  Clear();          // Call also this!
  FILLBUFF = false; // No autorange, explicit bounds provided
}

// Print the histogram and its binned statistics
void MH2::Print() const {
  if (!(fills > 0)) { // No fills
    std::cout << "MH2::Print: <" << name << "> Fills = " << fills << std::endl;
    return;
  }

  // Histogram name
  std::cout << "MH2::Print: <" << name << ">" << std::endl;

  std::vector<std::string> ascziart = {" ", ".", ":", "-", "=",
                                       "+", "*", "#", "%", "@"};
  const double maxw = GetMaxWeight();

  std::cout << "          |"; // Empty top left corner
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    std::cout << "=";
  }
  std::cout << "| " << std::endl;

  for (int j = YBINS - 1; j > -1;
       --j) { // Must be int, we index down to negative
    const double binwidth = (YMAX - YMIN) / YBINS;
    printf("%9.2E |", binwidth * (j + 1) - binwidth / 2.0 + YMIN);
    for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
      const double ratio = (maxw > 0) ? weights[i][j] / maxw : 0.0;
      const double w = std::isfinite(ratio) ? std::clamp(ratio, 0.0, 1.0) : 0.0;
      const int ind = std::clamp(static_cast<int>(std::round(w * 9)), 0, 9);
      std::cout << rang::fg::yellow << ascziart[ind] << rang::fg::reset;
    }
    std::cout << "|";
    std::cout << std::endl;
  }
  std::cout << "          |"; // Empty bottom left corner
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    std::cout << "=";
  }
  std::cout << "|" << std::endl;

  const double binwidth = (XMAX - XMIN) / XBINS;
  for (int k = -1; k < 50; ++k) {
    std::cout << "           ";
    int empty = 0;
    for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
      const double binvalue = binwidth * (i + 1) - binwidth / 2.0 + XMIN;
      const std::string let = gra::aux::ToString(std::abs(binvalue), 2);
      const std::string sgn = (binvalue >= 0.0) ? "+" : "-";
      if (k == -1) { // print value sign (-+)
        std::cout << sgn;
        continue;
      }
      if (k < static_cast<int>(let.length())) {
        std::cout << let[k];
      } else {
        std::cout << " ";
        empty++;
      }
    }
    std::cout << std::endl;
    if (empty == XBINS) {
      break;
    } // the end
  }

  // Print statistics
  std::cout << "<binned statistics> " << std::endl;
  std::pair<double, double> valerr = WeightMeanAndError();
  printf(" <W> = %0.3E +- %0.3E  [F = %lld | X: U/O = %lld/%lld, Y: U/O = "
         "%lld/%lld ]\n",
         valerr.first, valerr.second, fills, underflow[0], overflow[0],
         underflow[1], overflow[1]);

  const double X_mean = GetMeanX(1);
  const double X_sqmean = GetMeanX(2);

  const double Y_mean = GetMeanY(1);
  const double Y_sqmean = GetMeanY(2);

  printf(" <X> = %0.3f, <X^2> = %0.3f, <X^2> - <X>^2 = %0.3f \n", X_mean,
         X_sqmean, X_sqmean - std::pow(X_mean, 2));
  printf(" <Y> = %0.3f, <Y^2> = %0.3f, <Y^2> - <Y>^2 = %0.3f \n", Y_mean,
         Y_sqmean, Y_sqmean - std::pow(Y_mean, 2));

  std::cout << std::endl;
}

// Compute the mean event weight and its standard error
// <w> = sum_i w_i/N, delta<w> = sqrt[(<w^2>-<w>^2)/N]
std::pair<double, double> MH2::WeightMeanAndError() const {
  if (fills <= 0) {
    return {0.0, 0.0};
  }
  const double N =
      fills; // Need to use number of total fills here, not counts in bins
  const double val = SumWeights() / N;
  const double err2 = SumWeights2() / N - gra::math::pow2(val);
  const double err = gra::math::msqrt(err2 / N);

  return {val, err};
}

// Compute whether both indices address a physical bin
bool MH2::ValidBin(int xbin, int ybin) const {
  if (xbin >= 0 && xbin < XBINS && ybin >= 0 && ybin < YBINS) {
    return true;
  }
  return false;
}

// Compute whether two histograms use identical physical bin boundaries
bool MH2::CompatibleBinning(const MH2 &rhs) const {
  return XBINS == rhs.XBINS && YBINS == rhs.YBINS &&
         std::is_eq(XMIN <=> rhs.XMIN) && std::is_eq(XMAX <=> rhs.XMAX) &&
         std::is_eq(YMIN <=> rhs.YMIN) && std::is_eq(YMAX <=> rhs.YMAX) &&
         LOGX == rhs.LOGX && LOGY == rhs.LOGY;
}

// Fill one unweighted histogram entry
void MH2::Fill(double xvalue, double yvalue) {
  // Call weighted fill with weight 1.0
  Fill(xvalue, yvalue, 1.0);
}

// Fill one weighted histogram entry
void MH2::Fill(double xvalue, double yvalue, double weight) {
  if (std::isnan(weight)) {
    std::cout << "MH2:Fill: Warning: NaN weight" << std::endl;
    return;
  }
  if (std::isinf(weight)) {
    std::cout << "MH2:Fill: Warning: Inf weight" << std::endl;
    return;
  }

  if (!FILLBUFF) { // Normal filling

    ++fills;

    // Find out bins
    const int xbin = GetIdx(xvalue, XMIN, XMAX, XBINS, LOGX);
    const int ybin = GetIdx(yvalue, YMIN, YMAX, YBINS, LOGY);

    if (xbin == -3 || ybin == -3) {
      nanflow += 1;
    }

    if (xbin == -1) {
      underflow[0] += 1;
    }
    if (ybin == -1) {
      underflow[1] += 1;
    }

    if (xbin == -2) {
      overflow[0] += 1;
    }
    if (ybin == -2) {
      overflow[1] += 1;
    }

    if (ValidBin(xbin, ybin)) {
      weights[xbin][ybin] += weight;
      weights2[xbin][ybin] += weight * weight;
      counts[xbin][ybin] += 1;
    }
  } else { // Autorange initialization

    buff_values.push_back({xvalue, yvalue});
    buff_weights.push_back(weight);

    if (buff_values.size() >= static_cast<std::size_t>(AUTOBUFFSIZE)) {
      FlushBuffer();
    }
  }
}

// Clear the histogram data while preserving bin boundaries
void MH2::Clear() {
  if (FILLBUFF) {
    buff_values.clear();
    buff_weights.clear();
  }
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      weights[i][j] = 0;
      weights2[i][j] = 0;
      counts[i][j] = 0;
    }
  }
  fills = 0;
  underflow = {0, 0};
  overflow = {0, 0};
  nanflow = 0;
}

// Automatic histogram range algorithm
void MH2::FlushBuffer() {
  if (!FILLBUFF) {
    return;
  }
  if (buff_values.empty()) {
    ResetBounds(XBINS, -0.5, 0.5, YBINS, -0.5, 0.5);
    return;
  }
  FILLBUFF = false; // no more filling buffer

  std::vector<double> min = {0.0, 0.0};
  std::vector<double> max = {0.0, 0.0};

  // Loop over dimensions
  for (std::size_t dim = 0; dim < 2; ++dim) {
    std::vector<double> finite_values;
    finite_values.reserve(buff_values.size());
    double mu = 0.0;
    double sumW = 0.0;
    for (std::size_t i = 0; i < buff_values.size(); ++i) {
      if (!std::isfinite(buff_values[i][dim])) {
        continue;
      }
      finite_values.push_back(buff_values[i][dim]);
      const double magnitude = std::abs(buff_weights[i]);
      if (std::isfinite(magnitude)) {
        mu += buff_values[i][dim] * magnitude;
        sumW += magnitude;
      }
    }
    if (finite_values.empty()) {
      min[dim] = -0.5;
      max[dim] = 0.5;
      continue;
    }
    if (sumW > 0.0) {
      mu /= sumW;
    } else {
      mu = std::accumulate(finite_values.begin(), finite_values.end(), 0.0) /
           finite_values.size();
    }

    // Variance
    double var = 0.0;
    for (std::size_t i = 0; i < buff_values.size(); ++i) {
      const double magnitude = std::abs(buff_weights[i]);
      if (std::isfinite(buff_values[i][dim]) && std::isfinite(magnitude)) {
        var += magnitude * std::pow(buff_values[i][dim] - mu, 2);
      }
    }
    if (sumW > 0.0) {
      var /= sumW;
    } else {
      var = 0.0;
      for (const double value : finite_values) {
        var += std::pow(value - mu, 2);
      }
      var /= finite_values.size();
    }

    // Find minimum and maximum
    const auto [minimum, maximum] =
        std::minmax_element(finite_values.begin(), finite_values.end());
    const double minval = *minimum;
    const double maxval = *maximum;

    // Set new histogram bounds
    const double sigma = std::sqrt(std::abs(var));
    double xmin = mu - 2.5 * sigma;
    double xmax = mu + 2.5 * sigma;

    // A numerical failure may happen with variance calculation, then use this
    if (!std::isfinite(xmin) || !std::isfinite(xmax)) {
      xmin = minval;
      xmax = maxval;
    }

    // If symmetric setup set by user
    if (AUTOSYMMETRY[dim]) {
      const double val = std::max(std::abs(xmin), std::abs(xmax));
      xmin = -val;
      xmax = val;
    }

    // We have only positive values, such as invariant mass
    if (!AUTOSYMMETRY[dim] && minval > 0.0) {
      xmin = std::max(0.0, xmin);
    }

    if (!(xmax > xmin)) {
      const double half_width =
          std::max(1e-9, std::max(1.0, std::abs(mu)) * 1e-6);
      xmin = mu - half_width;
      xmax = mu + half_width;
    }

    min[dim] = xmin;
    max[dim] = xmax;
  }

  // New histogram bounds
  ResetBounds(XBINS, min[0], max[0], YBINS, min[1], max[1]);

  // Fill buffered events
  for (std::size_t i = 0; i < buff_values.size(); ++i) {
    Fill(buff_values[i][0], buff_values[i][1], buff_weights[i]);
  }

  // Clear buffers
  buff_values.clear();
  buff_weights.clear();
}

// Compute the binned moment of the x coordinate
// <x^p> = sum_ij w_ij x_i^p/sum_ij w_ij
double MH2::GetMeanX(int power) const {
  double sum = 0.0;
  double norm = 0.0; // Normalization
  const double minimum = LOGX ? std::log(XMIN) : XMIN;
  const double step = ((LOGX ? std::log(XMAX) : XMAX) - minimum) / XBINS;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    const double lower = minimum + i * step;
    const double upper = minimum + (i + 1) * step;
    const double center = LOGX ? 0.5 * (std::exp(lower) + std::exp(upper)) : 0.5 * (lower + upper);
    const double value = std::pow(center, power);
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      const double weight = GetBinWeight(i, j);
      sum += weight * value;
      norm += weight;
    }
  }
  if (!(std::fpclassify(norm) == FP_ZERO) && std::isfinite(norm)) {
    return sum / norm;
  } else {
    return 0.0;
  }
}

// Compute the binned moment of the y coordinate
// <y^p> = sum_ij w_ij y_j^p/sum_ij w_ij
double MH2::GetMeanY(int power) const {
  double sum = 0.0;
  double norm = 0.0; // Normalization
  const double minimum = LOGY ? std::log(YMIN) : YMIN;
  const double step = ((LOGY ? std::log(YMAX) : YMAX) - minimum) / YBINS;
  for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
    const double lower = minimum + j * step;
    const double upper = minimum + (j + 1) * step;
    const double center = LOGY ? 0.5 * (std::exp(lower) + std::exp(upper)) : 0.5 * (lower + upper);
    const double value = std::pow(center, power);
    for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
      const double weight = GetBinWeight(i, j);
      sum += weight * value;
      norm += weight;
    }
  }
  if (!(std::fpclassify(norm) == FP_ZERO) && std::isfinite(norm)) {
    return sum / norm;
  } else {
    return 0.0;
  }
}

// Compute the sum of all in-range bin weights
double MH2::SumWeights() const { return weights.Sum(); }

// Compute the sum of squared in-range event weights
double MH2::SumWeights2() const { return weights2.Sum(); }

// Compute the number of in-range histogram fills
long long int MH2::SumBinCounts() const { return counts.Sum(); }

// Compute the largest bin weight
double MH2::GetMaxWeight() const {
  double maxval = -std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      if (weights[i][j] > maxval) {
        maxval = weights[i][j];
      }
    }
  }
  return maxval;
}

// Compute the smallest bin weight
double MH2::GetMinWeight() const {
  double minval = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      if (weights[i][j] < minval) {
        minval = weights[i][j];
      }
    }
  }
  return minval;
}

// Compute the event count in one bin
long long int MH2::GetBinCount(int xbin, int ybin) const {
  if (ValidBin(xbin, ybin)) {
    return counts[xbin][ybin];
  } else {
    return 0;
  }
}

// Get weight of the bin
double MH2::GetBinWeight(int xbin, int ybin) const {
  if (ValidBin(xbin, ybin)) {
    return weights[xbin][ybin];
  } else {
    return 0;
  }
}

// Get bin indices (i,j) corresponding to value (xvalue,yvalue)
void MH2::GetBinIdx(double xvalue, double yvalue, int &xbin, int &ybin) const {
  // Find out bins
  xbin = GetIdx(xvalue, XMIN, XMAX, XBINS, LOGX);
  ybin = GetIdx(yvalue, YMIN, YMAX, YBINS, LOGY);
}

// Get table/histogram index for linearly or base-10 logarithmically spaced
// bins
// Gives exact uniform filling within bin boundaries
//
// In the logarithmic case, MINVAL and MAXVAL > 0, naturally
//
//
// Underflow returns -1
// Overflow  returns -2
// nan/inf   returns -3
//
int MH2::GetIdx(double value, double minval, double maxval, int nbins,
                bool logbins) const {
  if (std::isnan(value) || std::isinf(value)) {
    return -3;
  }
  if (value < minval) {
    return -1;
  } // Underflow
  if (value > maxval) {
    return -2;
  } // Overflow
  if (std::is_eq(value <=> maxval)) {
    return nbins - 1;
  }
  // Widen the coordinate differences and logarithms before dividing
  const long double lower = minval;
  const long double upper = maxval;
  const long double coordinate = value;
  const long double span = logbins ? std::log10(upper) - std::log10(lower) : upper - lower;
  const long double offset = logbins ? std::log10(coordinate) - std::log10(lower) : coordinate - lower;
  if (!(span > 0.0L) || !std::isfinite(span) || !std::isfinite(offset)) { return -3; }
  const long double bin = std::floor(nbins * (offset / span));
  return static_cast<int>(std::clamp(bin, 0.0L, static_cast<long double>(nbins - 1)));
}

// Compute the Shannon entropy of the histogram in bits
// H = -sum_ij p_ij log2(p_ij)
// (useful for validating statistical properties, for example)
double MH2::ShannonEntropy() const {
  // Match the absolute-weight probability convention used by MH1
  const auto values = weights.Elements();
  const long double sum = std::accumulate(values.begin(), values.end(), 0.0L,
                                          [](long double total, double value) { return total + std::abs(value); });
  if (!(sum > 0.0L)) { return 0.0; }
  long double entropy = 0.0L;
  for (const auto value : values) {
    const long double p = std::abs(value) / sum;
    if (p > 0.0L) { entropy -= p * std::log2(p); }
  }
  return static_cast<double>(entropy);
}

} // namespace gra
