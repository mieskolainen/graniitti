// Templated 1D-histogram class with real or complex weights
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <iostream>
#include <iterator>
#include <limits>
#include <numeric>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MH1.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;

namespace gra {
// Constructor
template <class T>
MH1<T>::MH1(int xbins, double xmin, double xmax, std::string namestr) {
  name = namestr;
  // Init
  ResetBounds(xbins, xmin, xmax);
  FILLBUFF = false;
}

// Constructor with only number of bins
template <class T> MH1<T>::MH1(int xbins, std::string namestr) {
  name = namestr;
  ResetBounds(xbins);
}

// Constructor with varying bin sizes [N x 2]
template <class T>
MH1<T>::MH1(std::vector<std::vector<double>> edges, std::string namestr) {
  name = namestr;

  // Init
  ResetBounds(edges);
  FILLBUFF = false;
}

// Empty constructor
template <class T> MH1<T>::MH1() { ResetBounds(50); }

// Destructor
template <class T> MH1<T>::~MH1() = default;

// Extract physical bin bounds and bin values with statistical errors
template <class T>
void MH1<T>::VectorOutput(std::vector<std::vector<double>> &bins,
                          std::vector<std::vector<double>> &vals) const {
  bins.resize(XBINS);
  vals.resize(XBINS);
  for (const auto i : indices(bins)) {
    bins[i] = {GetBinXVal(i, -1), GetBinXVal(i, 1)};
    vals[i] = {BinValue(i), BinError(i)};
  }
}

// Raw output to a file (or to the screen if file = stdout)
template <class T> void MH1<T>::RawOutput(FILE *file) const {
  if (file == nullptr) {
    throw std::invalid_argument("MH1::RawOutput: file pointer is null");
  }
  fprintf(file, "#binlow,binhigh,value,error,value/binwidth,error/binwidth \n");

  for (std::size_t idx = 0; idx < static_cast<unsigned int>(XBINS); ++idx) {
    // This will take care of amplitude versus amplitude squared
    const double value = BinValue(idx);
    const double value_err = BinError(idx);

    const double lower = GetBinXVal(idx, -1);
    const double upper = GetBinXVal(idx, 1);
    const double width = upper - lower;

    fprintf(file, "%0.6E,%0.6E,%0.6E,%0.6E,%0.6E,%0.6E\n", lower, upper, value,
            value_err, value / width, value_err / width);
  }
}

// Print screen
template <class T> void MH1<T>::Print(double width) const {
  if (!std::isfinite(width) || width <= 0.0) {
    throw std::invalid_argument(
        "MH1::Print: width must be finite and positive");
  }
  if (!(fills > 0)) { // No fills
    std::cout << "MH1::Print: <" << name << "> Fills = " << fills << std::endl;
    return;
  }

  // Histogram name
  std::cout << "MH1::Print: <" << name << ">" << std::endl;

  const double columns = XBINS * width;
  if (!std::isfinite(columns) ||
      columns > static_cast<double>(std::numeric_limits<long long>::max())) {
    throw std::length_error("MH1::Print: requested width is too large");
  }
  const std::size_t N =
      std::max<std::size_t>(2, static_cast<std::size_t>(std::llround(columns)));
  const double maximum = GetMaxWeight();
  const double maxvisual =
      std::isfinite(maximum) && maximum > 0.0 ? 1.1 * maximum : 1.0;

  // Top left corner
  std::cout << "                       |";
  for (std::size_t i = 0; i < N; ++i) {
    std::cout << "=";
  }
  std::cout << "| " << std::endl;

  for (std::size_t idx = 0; idx < static_cast<unsigned int>(XBINS); ++idx) {
    printf("(%9.2E, %9.2E] |", GetBinXVal(idx, -1), GetBinXVal(idx, 1));

    // This will take care of amplitude versus amplitude squared
    const double value = BinValue(idx);
    const double value_err = BinError(idx);

    // Visualization values scaled between [0,maxw]
    const double w =
        std::isfinite(value) ? std::clamp(value / maxvisual, 0.0, 1.0) : 0.0;

    // Y-AXIS index
    const std::size_t ind = std::round(w * (N - 1));

    const int Nc = GetBinCount(idx);
    if (Nc > 0) { // non-zero bin

      // Error on value (only on double histograms)
      const double relative_err =
          value > 0.0 && std::isfinite(value_err) ? value_err / value : 0.0;
      const std::size_t U = static_cast<std::size_t>(
          std::clamp(std::ceil(w * (N - 1) * (1.0 + relative_err)), 0.0,
                     static_cast<double>(N - 1)));
      const std::size_t D = static_cast<std::size_t>(
          std::clamp(std::floor(w * (N - 1) * (1.0 - relative_err)), 0.0,
                     static_cast<double>(N - 1)));
      // Print now
      for (std::size_t j = 0; j < N; ++j) {
        if (j == ind) {
          std::cout << rang::bg::magenta << "+" << rang::bg::reset;
        } else if (j >= D && j != ind && j <= U) {
          std::cout << rang::bg::magenta << " " << rang::bg::reset;
        } else {
          std::cout << " ";
        }
      }
      // No events in the bin
    } else {
      for (std::size_t j = 0; j < N; ++j) {
        std::cout << " ";
      }
    }

    printf("| %0.1E +- %0.1E [dO/dx] = %0.1E \n", value, value_err,
           GetdOdX(idx));
  }

  // Empty bottom left corner
  std::cout << "                       |";

  for (std::size_t i = 0; i < N; ++i) {
    std::cout << "=";
  }
  std::cout << "| " << std::endl;

  // Print statistics
  std::cout << "<binned statistics>" << std::endl;
  std::pair<double, double> valerr = WeightMeanAndError();
  printf(" <W> = %0.3E +- %0.3E  [F = %lld | U/O = %lld/%lld]\n", valerr.first,
         valerr.second, fills, underflow, overflow);

  const double mean = GetMean(1);
  const double sqmean = GetMean(2);
  printf(" <X> = %0.3f, <X^2> = %0.3f, <X^2> - <X>^2 = %0.3f\n", mean, sqmean,
         sqmean - std::pow(mean, 2));

  std::cout << std::endl;
}

// Compute <w>=sum_i w_i/N and delta<w>=sqrt[(<|w|^2>-|<w>|^2)/N]
template <class T>
std::pair<double, double> MH1<T>::WeightMeanAndError() const {
  if (fills <= 0) {
    return {0.0, 0.0};
  }
  const double N =
      fills; // Need to use number of total fills here, not counts in bins
  double val = 0.0;
  double second_moment = 0.0;
  if constexpr (std::is_same_v<T, double>) {
    val = SumWeights() / N;
    second_moment = SumWeights2() / N;
  } else {
    val = std::abs(SumWeights()) / N;
    second_moment = std::abs(SumWeights2()) / N;
  }
  const double err2 = second_moment - gra::math::pow2(val);
  const double err = gra::math::msqrt(err2 / N);

  return {val, err};
}

// Get valarrays containing X values and corresponding density values
template <class T>
void MH1<T>::GetXPositiveDefinite(std::valarray<double> &x,
                                  std::valarray<double> &y) const {
  x.resize(XBINS);
  y.resize(XBINS);
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    x[i] = GetBinXVal(i, 0); // Get center value
    y[i] = GetPositiveDefinite(i);
  }
}

// Compute a histogram moment from the binned values
// <x^p> = sum_i w_i x_i^p/sum_i w_i
template <class T> double MH1<T>::GetMean(int power) const {
  double sum = 0.0;
  double norm = 0.0; // Normalization

  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    const double X_value = GetBinXVal(i, 0);

    const double weight = BinIntegral(i);
    sum += weight * std::pow(X_value, power);
    norm += weight;
  }
  if (!(std::fpclassify(norm) == FP_ZERO) && std::isfinite(norm)) {
    return sum / norm;
  } else {
    return 0.0;
  }
}

// Compute whether the index addresses a physical bin
template <class T> bool MH1<T>::ValidBin(int idx) const {
  if (idx >= 0 && idx < XBINS) {
    return true;
  }
  return false;
}

// Compute whether two histograms use identical physical bin boundaries
template <class T> bool MH1<T>::CompatibleBinning(const MH1<T> &rhs) const {
  const bool same_edges =
      binedges.size() == rhs.binedges.size() &&
      std::equal(binedges.cbegin(), binedges.cend(), rhs.binedges.cbegin(),
                 [](const auto &left, const auto &right) {
                   return left.size() == right.size() &&
                          std::equal(left.cbegin(), left.cend(), right.cbegin(),
                                     [](const double first,
                                        const double second) {
                                       return std::is_eq(first <=> second);
                                     });
                 });
  return XBINS == rhs.XBINS && std::is_eq(XMIN <=> rhs.XMIN) &&
         std::is_eq(XMAX <=> rhs.XMAX) && LOGX == rhs.LOGX && same_edges;
}

// Unweighted fill
template <class T> void MH1<T>::Fill(double xvalue) {
  // Call the weighted with weight 1.0
  Fill(xvalue, 1.0);
}

// Weighted fill
template <class T> void MH1<T>::Fill(double xvalue, T weight) {
  const double magnitude = std::abs(weight);
  if (!std::isfinite(magnitude)) {
    nanflow += 1;
    return;
  }
  if (!FILLBUFF) { // Normal filling

    fills += 1;

    // Find out bin
    int idx = 0;
    GetBinIdx(xvalue, idx);

    if (idx == -3) {
      nanflow += 1;
    }
    if (idx == -1) {
      underflow += 1;
    }
    if (idx == -2) {
      overflow += 1;
    }

    if (ValidBin(idx)) {
      weights[idx] += weight;
      weights2[idx] += magnitude * magnitude;
      if constexpr (std::is_same_v<T, std::complex<double>>) { weights_sq[idx] += weight * weight; }
      counts[idx] += 1;
    }
  } else { // Autorange initialization
    buff_values.push_back(xvalue);
    buff_weights.push_back(weight);

    if (buff_values.size() >= static_cast<std::size_t>(AUTOBUFFSIZE)) {
      FlushBuffer();
    }
  }
}

// Automatic histogram range algorithm
template <class T> void MH1<T>::FlushBuffer() {
  if (!FILLBUFF) {
    return;
  }
  if (buff_values.empty()) {
    ResetBounds(XBINS, -0.5, 0.5);
    return;
  }
  FILLBUFF = false; // no more filling buffer

  std::vector<double> finite_values;
  finite_values.reserve(buff_values.size());
  double mu = 0.0;
  double sumW = 0.0;
  for (std::size_t i = 0; i < buff_values.size(); ++i) {
    if (!std::isfinite(buff_values[i])) {
      continue;
    }
    finite_values.push_back(buff_values[i]);
    const double magnitude = std::abs(buff_weights[i]);
    if (std::isfinite(magnitude)) {
      mu += buff_values[i] * magnitude;
      sumW += magnitude;
    }
  }
  if (finite_values.empty()) {
    ResetBounds(XBINS, -0.5, 0.5);
    for (std::size_t i = 0; i < buff_values.size(); ++i) {
      Fill(buff_values[i], buff_weights[i]);
    }
    buff_values.clear();
    buff_weights.clear();
    return;
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
    if (std::isfinite(buff_values[i]) && std::isfinite(magnitude)) {
      var += magnitude * std::pow(buff_values[i] - mu, 2);
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

  // Minimum and maximum
  auto it1 = std::min_element(finite_values.begin(), finite_values.end());
  auto it2 = std::max_element(finite_values.begin(), finite_values.end());
  const double minval = *it1;
  const double maxval = *it2;

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
  if (AUTOSYMMETRY) {
    const double val = std::max(std::abs(xmin), std::abs(xmax));
    xmin = -val;
    xmax = val;
  }

  // We have only positive values, such as invariant mass
  if (!AUTOSYMMETRY && minval > 0.0) {
    xmin = std::max(0.0, xmin);
  }

  if (!(xmax > xmin)) {
    const double half_width =
        std::max(1e-9, std::max(1.0, std::abs(mu)) * 1e-6);
    xmin = mu - half_width;
    xmax = mu + half_width;
  }

  ResetBounds(XBINS, xmin, xmax);

  // Fill buffered events
  for (std::size_t i = 0; i < buff_values.size(); ++i) {
    Fill(buff_values[i], buff_weights[i]);
  }

  // Clear buffers
  buff_values.clear();
  buff_weights.clear();
}

// Reset histogram completely
template <class T> void MH1<T>::ResetBounds(int xbins) {
  if (xbins <= 0) {
    throw std::invalid_argument(
        "MH1<T>::ResetBounds: bin count must be positive");
  }
  LOGX = false;
  XMIN = 0.0;
  XMAX = 0.0;
  XBINS = xbins;
  binedges.clear();
  weights.assign(static_cast<std::size_t>(XBINS), T{});
  weights2.assign(static_cast<std::size_t>(XBINS), T{});
  counts.assign(static_cast<std::size_t>(XBINS), 0);
  buff_values.clear();
  buff_weights.clear();
  Clear();
  FILLBUFF = true;
}

// Reset histogram completely
template <class T>
void MH1<T>::ResetBounds(std::vector<std::vector<double>> edges) {
  if (edges.empty()) {
    throw std::invalid_argument(
        "MH1<T>::ResetBounds: input edge array is empty " +
        std::to_string(edges.size()));
  }
  if (edges.size() >
      static_cast<std::size_t>(std::numeric_limits<int>::max())) {
    throw std::length_error(
        "MH1<T>::ResetBounds: bin count is not representable");
  }
  for (std::size_t i = 0; i < edges.size(); ++i) {
    if (edges[i].size() != 2) {
      throw std::invalid_argument(
          "MH1<T>::ResetBounds: every element of the edge array should size 2 "
          "(minimum, maximum)");
    }
    if (!std::isfinite(edges[i][0]) || !std::isfinite(edges[i][1]) ||
        edges[i][1] <= edges[i][0] ||
        (i > 0 && edges[i][0] < edges[i - 1][1])) {
      throw std::invalid_argument(
          "MH1<T>::ResetBounds: non-monotonic binedges at index = " +
          std::to_string(i));
    }
  }

  LOGX = false;
  ResetBounds(static_cast<int>(edges.size()), edges.front()[0],
              edges.back()[1]);
  binedges = std::move(edges);
}

// Reset histogram completely
template <class T>
void MH1<T>::ResetBounds(int xbins, double xmin, double xmax) {
  if (xbins <= 0 || !std::isfinite(xmin) || !std::isfinite(xmax) ||
      !(xmax > xmin) || (LOGX && !(xmin > 0.0))) {
    throw std::invalid_argument(
        "MH1<T>::ResetBounds: require positive bins, finite ordered bounds and positive logarithmic bounds");
  }
  XMIN = xmin;
  XMAX = xmax;
  XBINS = xbins;
  binedges.clear();

  // Init
  std::vector<T> null(static_cast<std::size_t>(XBINS), T{});
  weights = null;
  weights2 = null;
  counts = std::vector<long long int>(static_cast<std::size_t>(XBINS), 0);

  Clear();          // Call also this!
  FILLBUFF = false; // No autorange, explicit bounds provided
}

// Clear the histogram data but keep the bounds
template <class T> void MH1<T>::Clear() {
  if (FILLBUFF) {
    buff_values.clear();
    buff_weights.clear();
  }
  if constexpr (std::is_same_v<T, std::complex<double>>) { weights_sq.assign(XBINS, T{}); }
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    weights[i] = 0.0;
    weights2[i] = 0.0;
    counts[i] = 0;
  }
  fills = 0;
  underflow = 0;
  overflow = 0;
  nanflow = 0;
}

// Sum over all histogram bin weights
template <class T> T MH1<T>::SumWeights() const { return gra::Sum(weights); }

// Sum over all histogram bin weights squared
template <class T> T MH1<T>::SumWeights2() const { return gra::Sum(weights2); }

// Sum the number of bin counts (not the same as fills)
template <class T> long long int MH1<T>::SumBinCounts() const {
  return gra::Sum(counts);
}

// Get maximum histogram bin weight, for complex return |w|^2
template <class T> double MH1<T>::GetMaxWeight() const {
  double maxval = 0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    if (GetPositiveDefinite(i) > maxval) {
      maxval = GetPositiveDefinite(i);
    }
  }
  return maxval;
}

// Get minimum histogram bin weight, for complex return |w|^2
template <class T> double MH1<T>::GetMinWeight() const {
  double minval = std::numeric_limits<double>::infinity();
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    if (GetPositiveDefinite(i) < minval) {
      minval = GetPositiveDefinite(i);
    }
  }
  return minval;
}

// Get number of event fills in the bin
template <class T> long long int MH1<T>::GetBinCount(int idx) const {
  if (ValidBin(idx)) {
    return counts[idx];
  } else {
    return 0;
  }
}

// Get weight of the bin, for complex return complex number
template <class T> T MH1<T>::GetBinWeight(int idx) const {
  if (ValidBin(idx)) {
    return weights[idx];
  } else {
    return 0.0;
  }
}

// Get weight of the bin, for complex return complex number
template <class T> T MH1<T>::GetBinWeight2(int idx) const {
  if (ValidBin(idx)) {
    return weights2[idx];
  } else {
    return 0.0;
  }
}

// Compute the statistical error estimate of one bin
// delta w_bin = sqrt(sum_i |w_i|^2)
template <class T> double MH1<T>::GetBinError(int idx) const {
  if (!ValidBin(idx)) {
    return 0.0;
  }
  double err = 0;
  if (GetBinCount(idx) > 0) {
    err = std::sqrt(std::abs(GetBinWeight2(idx))); // [\sum_i |w_i|^2]^{1/2}
  }
  return err;
}

// Compute the signed density or the squared coherent amplitude density
// dO/dx = w_bin/(N_fill Delta x), or |w_bin/(N_fill Delta x)|^2 for amplitudes
template <class T> double MH1<T>::GetdOdX(int idx) const {
  if (!ValidBin(idx)) {
    return 0.0;
  }
  const double norm = fills * GetBinWidth(idx);
  if (norm > 0.0) {
    if constexpr (std::is_same_v<T, double>) {
      return GetBinWeight(idx) / norm;
    } else {
      return std::norm(GetBinWeight(idx) / norm);
    }
  } else {
    return 0.0;
  }
}

// Compute the uncertainty as a differential distribution
// Apply the same linear or quadratic normalization as the bin observable
template <class T> double MH1<T>::GetdOdXError(int idx) const {
  if (!ValidBin(idx)) {
    return 0.0;
  }
  const double norm = fills * GetBinWidth(idx);
  if (norm > 0.0) {
    if constexpr (std::is_same_v<T, double>) {
      return BinError(idx) / norm;
    } else {
      return BinError(idx) / norm / norm;
    }
  } else {
    return 0.0;
  }
}

// Compute bin normalization
template <class T> double MH1<T>::GetBinWidth(int idx) const {
  if (!ValidBin(idx)) {
    throw std::out_of_range("MH1::GetBinWidth: bin index outside histogram");
  }
  if (binedges.empty()) {
    return GetBinXVal(idx, 1) - GetBinXVal(idx, -1);
  } else {
    return binedges[idx][1] - binedges[idx][0];
  }
}

// Compute signed real bin content or the coherent complex-weight norm
template <class T> double MH1<T>::BinValue(int idx) const {
  if constexpr (std::is_same_v<T, double>) {
    return GetBinWeight(idx);
  } else {
    return std::norm(GetBinWeight(idx));
  }
}

// Integrate the bin observable up to its common sample normalization
// A coherent amplitude integral contributes |w_bin|^2 / Delta x to the density integral
template <class T> double MH1<T>::BinIntegral(int idx) const {
  if constexpr (std::is_same_v<T, double>) {
    return GetBinWeight(idx);
  } else {
    if (FILLBUFF || !ValidBin(idx)) { return 0.0; }
    return BinValue(idx) / GetBinWidth(idx);
  }
}

// Propagate independent event-weight fluctuations to the coherent norm at first order
// Var(|A|^2) = 2 (|A|^2 sum_i |w_i|^2 + Re[(A*)^2 sum_i w_i^2])
template <class T> double MH1<T>::BinError(int idx) const {
  if constexpr (std::is_same_v<T, double>) {
    return GetBinError(idx);
  } else {
    if (!ValidBin(idx)) { return 0.0; }
    using Wide = std::complex<long double>;
    const Wide amplitude = weights[idx];
    const long double variance = 2.0L * (std::norm(amplitude) * std::real(weights2[idx]) +
        std::real(math::pow2(std::conj(amplitude)) * Wide(weights_sq[idx])));
    return std::sqrt(std::max(0.0L, variance));
  }
}

// Compute |w| (double) or |w|^2 (complex case)
template <class T> double MH1<T>::GetPositiveDefinite(int idx) const {
  if (ValidBin(idx)) {
    if constexpr (std::is_same_v<T, double>) {
      return std::abs(weights[idx]);
    } else {
      return std::norm(weights[idx]);
    }
  } else {
    return 0.0;
  }
}

// Get bin index (idx) corresponding to value (xvalue)
template <class T> void MH1<T>::GetBinIdx(double xvalue, int &idx) {
  // Equal width binning
  if (binedges.empty()) {
    idx = ComputeIdx(xvalue, XMIN, XMAX, XBINS, LOGX);

    // Variable width binning
  } else {
    idx = ComputeIdx(xvalue, binedges);
  }
}

// Get bin value in units of X for a given bin index
// boundary = -1,0,1 (lower, center, upper)
template <class T> double MH1<T>::GetBinXVal(int idx, int boundary) const {
  if (!ValidBin(idx)) {
    throw std::out_of_range("MH1::GetBinXVal: bin index outside histogram");
  }
  if (boundary < -1 || boundary > 1) {
    throw std::invalid_argument(
        "MH1::GetBinXVal: bin boundary must be -1, 0 or 1");
  }
  // Non-equal width binning
  if (!binedges.empty()) {
    if (boundary == -1)
      return binedges[idx][0];
    if (boundary == 0)
      return 0.5 * (binedges[idx][0] + binedges[idx][1]);
    if (boundary == 1)
      return binedges[idx][1];
  }

  // Equal binning
  if (LOGX) {
    const double log10step = (std::log10(XMAX) - std::log10(XMIN)) / XBINS;
    if (boundary == -1) {
      return std::pow(10, std::log10(XMIN) + idx * log10step);
    } else if (boundary == 0) {
      return (std::pow(10, std::log10(XMIN) + idx * log10step) +
              std::pow(10, std::log10(XMIN) + (idx + 1) * log10step)) /
             2;
    } else if (boundary == 1) {
      return std::pow(10, std::log10(XMIN) + (idx + 1) * log10step);
    }
  } else {
    const double binwidth = (XMAX - XMIN) / XBINS;
    const double value = XMIN + (idx + 1) * binwidth;
    if (boundary == -1) {
      return value - binwidth;
    } else if (boundary == 0) {
      return value - binwidth / 2.0;
    } else if (boundary == 1) {
      return value;
    }
  }
  throw std::logic_error("MH1::GetBinXVal: unreachable bin boundary");
}

// Compute bin index for variable width binning
//
template <class T>
int MH1<T>::ComputeIdx(double value,
                       const std::vector<std::vector<double>> &edges) const {
  if (std::isnan(value) || std::isinf(value)) {
    return -3;
  }

  // Underflow
  if (value < edges[0][0]) {
    return -1;
  }
  // Overflow
  else if (value > edges[edges.size() - 1][1]) {
    return -2;
  }
  // Binary search
  else {
    /*
    // Linear search
    for (std::size_t i = 0; i < binedges.size(); ++i) {
      if (binedges[i][0] < xvalue && xvalue <= binedges[i][1]) {
        idx = i;
      }
    }
    */
    const auto upper = std::upper_bound(
        edges.begin(), edges.end(), value,
        [](const double coordinate, const std::vector<double> &edge) {
          return coordinate < edge[0];
        });
    const auto it = std::prev(upper);
    if (value > (*it)[1]) {
      return -3;
    }
    return static_cast<int>(std::distance(edges.begin(), it));
  }
}

// Compute table/histogram index for linearly or base-10 logarithmically spaced
// bins. Computes exact uniform filling within bin boundaries
//
// In the logarithmic case, MINVAL and MAXVAL > 0, naturally
//
//
// Underflow     returns -1
// Overflow      returns -2
// Invalid input returns -3
//
template <class T>
int MH1<T>::ComputeIdx(double value, double minval, double maxval, int nbins,
                       bool logbins) const {
  if (std::isnan(value) || std::isinf(value)) {
    return -3;
  }
  if (value < minval) {
    return -1;
  } // underflow
  if (value > maxval) {
    return -2;
  } // overflow
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

// Instantiate (necessary for compilation)
template class MH1<double>;
template class MH1<std::complex<double>>;

} // namespace gra
