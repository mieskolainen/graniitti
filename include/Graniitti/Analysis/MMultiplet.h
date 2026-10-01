// Multiplet of ROOT histograms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MULTIPLETCLASS_H
#define MULTIPLETCLASS_H

// C++
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

// ROOT
#include "TCanvas.h"
#include "TH1.h"
#include "TH2.h"
#include "TProfile.h"

// Own
#include "Graniitti/Tech/MAux.h"

namespace gra {
// Histogram boundaries
class h1Bound {
public:
  h1Bound(unsigned int _N, double _min, double _max) {
    N = _N;
    min = _min;
    max = _max;
  }
  ~h1Bound() = default;

  unsigned int N = 0;
  double min = 0.0;
  double max = 0.0;
};

class h1Multiplet {
public:
  h1Multiplet(const std::string &name, const std::string &labeltext, int N,
              double minval, double maxval,
              const std::vector<std::string> &legendtext);
  h1Multiplet(const h1Multiplet &) = delete;
  h1Multiplet &operator=(const h1Multiplet &) = delete;
  ~h1Multiplet() {
    for (const auto &i : gra::aux::indices(h)) {
      delete h[i];
    }
  }
  void MultiFill(const std::vector<double> &x) {
    if (x.size() != h.size()) {
      throw std::invalid_argument(
          "h1Multiplet::MultiFill: input size mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i]);
    }
  }
  void MultiFill(const std::vector<double> &x,
                 const std::vector<double> &weight) {
    if (x.size() != h.size() || weight.size() != h.size()) {
      throw std::invalid_argument(
          "h1Multiplet::MultiFill: input size mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i], weight[i]);
    }
  }
  // Normalize all histograms with one event-weight factor per input sample
  void NormalizeAll(const std::vector<double> &normalization,
                    const std::vector<double> &multiplier);

  // Plot and save 1D-histogram Multiplet
  std::vector<double> SaveFig(const std::string &fullpath,
                              bool RATIOPLOT = true) const;

  int N_;
  double minval_;
  double maxval_;
  std::string legendposition_;
  std::string name_;

  // Histogram pointers here
  std::vector<TH1D *> h;
  std::vector<std::string> legendtext_;
};

class h2Multiplet {
public:
  h2Multiplet(const std::string &name, const std::string &labeltext, int N1,
              double minval1, double maxval1, int N2, double minval2,
              double maxval2, const std::vector<std::string> &legendtext);
  h2Multiplet(const h2Multiplet &) = delete;
  h2Multiplet &operator=(const h2Multiplet &) = delete;

  ~h2Multiplet() {
    for (const auto &i : gra::aux::indices(h)) {
      delete h[i];
    }
  }
  void MultiFill(const std::vector<std::vector<double>> &x) {
    if (x.size() != h.size() ||
        !std::all_of(x.begin(), x.end(),
                     [](const auto &value) { return value.size() == 2; })) {
      throw std::invalid_argument(
          "h2Multiplet::MultiFill: input shape mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i][0], x[i][1]);
    }
  }
  void MultiFill(const std::vector<std::vector<double>> &x,
                 const std::vector<double> &weight) {
    if (x.size() != h.size() || weight.size() != h.size() ||
        !std::all_of(x.begin(), x.end(),
                     [](const auto &value) { return value.size() == 2; })) {
      throw std::invalid_argument(
          "h2Multiplet::MultiFill: input shape mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i][0], x[i][1], weight[i]);
    }
  }
  // Normalize all histograms with one event-weight factor per input sample
  void NormalizeAll(const std::vector<double> &normalization,
                    const std::vector<double> &multiplier);

  // Plot and save histogram Multiplet
  double SaveFig(const std::string &fullpath, bool RATIOPLOT = true) const;

  int N1_;
  double minval1_;
  double maxval1_;

  int N2_;
  double minval2_;
  double maxval2_;

  std::string name_;

  // Histogram pointers here
  std::vector<TH2D *> h;
  std::vector<std::string> legendtext_;
};

class hProfMultiplet {
public:
  hProfMultiplet(const std::string &name, const std::string &labeltext, int N,
                 double minval1, double maxval1, double minval2, double maxval2,
                 const std::vector<std::string> &legendtext);
  hProfMultiplet(const hProfMultiplet &) = delete;
  hProfMultiplet &operator=(const hProfMultiplet &) = delete;

  ~hProfMultiplet() {
    for (const auto &i : gra::aux::indices(h)) {
      delete h[i];
    }
  }
  void MultiFill(const std::vector<std::vector<double>> &x) {
    if (x.size() != h.size() ||
        !std::all_of(x.begin(), x.end(),
                     [](const auto &value) { return value.size() == 2; })) {
      throw std::invalid_argument(
          "hProfMultiplet::MultiFill: input shape mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i][0], x[i][1]);
    }
  }
  void MultiFill(const std::vector<std::vector<double>> &x,
                 const std::vector<double> &weight) {
    if (x.size() != h.size() || weight.size() != h.size() ||
        !std::all_of(x.begin(), x.end(),
                     [](const auto &value) { return value.size() == 2; })) {
      throw std::invalid_argument(
          "hProfMultiplet::MultiFill: input shape mismatch");
    }
    for (const auto &i : gra::aux::indices(h)) {
      h[i]->Fill(x[i][0], x[i][1], weight[i]);
    }
  }

  // Plot and save histogram Multiplet
  double SaveFig(const std::string &fullpath, bool RATIOPLOT = true) const;

  int N_;
  double minval1_;
  double maxval1_;

  double minval2_;
  double maxval2_;

  std::string name_;
  std::string legendposition_;

  // Histogram pointers here
  std::vector<TProfile *> h;
  std::vector<std::string> legendtext_;
};

} // namespace gra

#endif
