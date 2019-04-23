// 2D histogram class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MH2.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"

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
  name = namestr;
  XBINS = xbins;
  YBINS = ybins;
  FILLBUFF = true;
}

// Empty constructor
MH2::MH2() {
  XBINS = 50; // Default
  YBINS = 50;
  FILLBUFF = true;
}

// Destructor
MH2::~MH2() {}

void MH2::ResetBounds(int xbins, double xmin, double xmax, int ybins,
                      double ymin, double ymax) {
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
}

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
      const double w = (maxw > 0) ? weights[i][j] / maxw : 0;
      const int ind = std::round(w * 9);
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
  printf(" <W> = %0.3E +- %0.3E \n", valerr.first, valerr.second);
  printf(" [F = %lld | X: U/O = %lld/%lld, Y: U/O = %lld/%lld] \n", fills,
         underflow[0], overflow[0], underflow[1], overflow[1]);

  std::cout << std::endl;
}

std::pair<double, double> MH2::WeightMeanAndError() const {
  const double N =
      fills; // Need to use number of total fills here, not counts in bins
  const double val = SumWeights() / N;
  const double err2 = SumWeights2() / N - gra::math::pow2(val);
  const double err = gra::math::msqrt(err2 / N);

  return {val, err};
}

bool MH2::ValidBin(int xbin, int ybin) const {
  if (xbin >= 0 && xbin < XBINS && ybin >= 0 && ybin < YBINS) {
    return true;
  }
  return false;
}

void MH2::Fill(double xvalue, double yvalue) {
  // Call weighted fill with weight 1.0
  Fill(xvalue, yvalue, 1.0);
}

void MH2::Fill(double xvalue, double yvalue, double weight) {
  if (!FILLBUFF) { // Normal filling

    ++fills;

    // Find out bins
    const int xbin = GetIdx(xvalue, XMIN, XMAX, XBINS, LOGX);
    const int ybin = GetIdx(yvalue, YMIN, YMAX, YBINS, LOGY);

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

    if (buff_values.size() > static_cast<unsigned int>(AUTOBUFFSIZE)) {
      FlushBuffer();
    }
  }
}

void MH2::Clear() {
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      weights[i][j] = 0;
      weights2[i][j] = 0;
      counts[i][j] = 0;
      fills = 0;
      underflow = {0, 0};
      overflow = {0, 0};
    }
  }
}

// Automatic histogram range algorithm
void MH2::FlushBuffer() {
  if (FILLBUFF && buff_values.size() > 0) {
    FILLBUFF = false; // no more filling buffer

    std::vector<double> min = {0.0, 0.0};
    std::vector<double> max = {0.0, 0.0};

    // Loop over dimensions
    for (std::size_t dim = 0; dim < 2; ++dim) {
      // Find out mean
      double mu = 0;
      double sumW = 0;
      for (std::size_t i = 0; i < buff_values.size(); ++i) {
        mu += std::abs(buff_values[i][dim] * buff_weights[i]);
        sumW += std::abs(buff_weights[i]);
      }
      if (sumW > 0) {
        mu /= sumW;
      }

      // Variance
      double var = 0;
      for (std::size_t i = 0; i < buff_values.size(); ++i) {
        var +=
            std::abs(buff_weights[i]) * std::pow(buff_values[i][dim] - mu, 2);
      }
      if (sumW > 0) {
        var /= sumW;
      }

      // Find minimum
      double minval = 1e64;
      for (std::size_t i = 0; i < buff_values.size(); ++i) {
        if (buff_values[i][dim] < minval) {
          minval = buff_values[i][dim];
        }
      }

      // Set new histogram bounds
      const double std = std::sqrt(std::abs(var));
      double xmin = mu - 2.5 * std;
      double xmax = mu + 2.5 * std;

      // If symmetric setup set by user
      if (AUTOSYMMETRY[dim]) {
        double val = (std::abs(xmin) + std::abs(xmax)) / 2.0;
        xmin = -val;
        xmax = val;
      }

      // We have only positive values, such as invariant mass
      if (minval > 0.0) {
        xmin = std::max(0.0, xmin);
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
}

double MH2::SumWeights() const {
  double sum = 0.0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      sum += weights[i][j];
    }
  }
  return sum;
}

double MH2::SumWeights2() const {
  double sum = 0.0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      sum += weights2[i][j];
    }
  }
  return sum;
}

long long int MH2::SumBinCounts() const {
  long long int sum = 0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      sum += counts[i][j];
    }
  }
  return sum;
}

double MH2::GetMaxWeight() const {
  double maxval = 0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      if (weights[i][j] > maxval) {
        maxval = weights[i][j];
      }
    }
  }
  return maxval;
}

double MH2::GetMinWeight() const {
  double minval = 1e128;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      if (weights[i][j] < minval) {
        minval = weights[i][j];
      }
    }
  }
  return minval;
}

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
// Gives exact uniform filling within bin boundaries.
//
// In the logarithmic case, MINVAL and MAXVAL > 0, naturally.
//
//
// Underflow returns -1
// Overflow  returns -2
int MH2::GetIdx(double value, double minval, double maxval, int nbins,
                bool logbins) const {
  if (value < minval) {
    return -1;
  } // Underflow
  if (value > maxval) {
    return -2;
  } // Overflow
  int idx = 0;

  // Logarithmic binning
  if (logbins) {
    // Check do we have a non-negative input
    idx = value > 0
              ? std::floor(nbins * (std::log10(value) - std::log10(minval)) /
                           (std::log10(maxval) - std::log10(minval)))
              : -1;
    // Linear binning
  } else {
    const double BINWIDTH = (maxval - minval) / nbins;
    idx = std::floor((value - minval) / BINWIDTH);
  }
  return idx;
}

// Return Shannon entropy of the histogram in bits (base-2 log)
// (useful for validating statistical properties, for example)
double MH2::ShannonEntropy() const {
  double sum = 0.0;
  for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
    for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
      sum += weights[i][j];
    }
  }
  if (sum > 0) {
    double S = 0.0;
    for (std::size_t i = 0; i < static_cast<unsigned int>(XBINS); ++i) {
      for (std::size_t j = 0; j < static_cast<unsigned int>(YBINS); ++j) {
        if (weights[i][j] > 0) {
          S += weights[i][j] / sum * std::log2(weights[i][j] / sum);
        }
      }
    }
    return -S;
  } else {
    return 0.0;
  }
}

} // gra namespace
