// 2D histogram class
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MH2_H
#define MH2_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/MMatrix.h"

namespace gra {
class MH2 {
 public:
  MH2(int xbins, double xmin, double xmax, int ybins, double ymin, double ymax,
      std::string namestr = "noname");

  MH2(int xbins, int ybins, std::string namestr = "noname");
  MH2();
  ~MH2();

  void ResetBounds(int xbins, double xmin, double xmax, int ybins, double ymin, double ymax);

  // Get full histogram data
  MMatrix<double>        GetWeights() const { return weights; }
  MMatrix<double>        GetWeights2() const { return weights2; }
  MMatrix<long long int> GetCounts() const { return counts; }

  void                      Fill(double xvalue, double yvalue);
  void                      Fill(double xvalue, double yvalue, double weight);
  void                      Clear();
  std::pair<double, double> WeightMeanAndError() const;

  double GetMeanX(int power) const;
  double GetMeanY(int power) const;

  double        SumWeights() const;
  double        SumWeights2() const;
  long long int SumBinCounts() const;
  long long int FillCount() const { return fills; }

  double        GetBinWeight(int xbin, int ybin) const;
  long long int GetBinCount(int xbin, int ybin) const;
  void          GetBinIdx(double xvalue, double yvalue, int &xbin, int &ybin) const;
  double        GetMaxWeight() const;
  double        GetMinWeight() const;
  void          Print() const;
  double        ShannonEntropy() const;

  // Set logarithmic binning
  // (user needs to take care that XMIN and XMAX > 0)
  void SetLogX() { LOGX = true; }
  void SetLogY() { LOGX = true; }
  void SetLogXY() {
    LOGX = true;
    LOGY = true;
  }

  // Overload + operator to add two histograms
  MH2 operator+(const MH2 &rhs) {
    if ((this->XBINS != rhs.XBINS) || (this->YBINS != rhs.YBINS)) {
      throw std::domain_error("MH2 + operator: Histograms with different number of bins");
    }
    MH2 h(this->XBINS, this->XMIN, this->XMAX, this->YBINS, this->YMIN, this->YMAX, this->name);

    h.fills     = this->fills + rhs.fills;
    h.underflow = {this->underflow[0] + rhs.underflow[0], this->underflow[1] + rhs.underflow[1]};
    h.overflow  = {this->overflow[0] + rhs.overflow[0], this->overflow[1] + rhs.overflow[1]};
    h.nanflow   = this->nanflow + rhs.nanflow;

    // DATA
    h.weights  = this->weights + rhs.weights;
    h.weights2 = this->weights2 + rhs.weights2;
    h.counts   = this->counts + rhs.counts;

    return h;
  }

  // AUTOBUFFSIZE for autorange buffer size
  void SetAutoBuffSize(int n) { AUTOBUFFSIZE = n; }
  void FlushBuffer();

  // Symmetric bounds
  void SetAutoSymmetry(const std::vector<bool> &in) {
    if (in.size() != 2) {
      throw std::invalid_argument("MH2::SetAutoSymmetry: Input should be size 2 boolean vector");
    }
    AUTOSYMMETRY = in;
  }

  void GetBounds(int &xbins, double &xmin, double &xmax, int &ybins, double &ymin,
                 double &ymax) const {
    xbins = XBINS;
    xmin  = XMIN;
    xmax  = XMAX;
    ybins = YBINS;
    ymin  = YMIN;
    ymax  = YMAX;
  }

  void FuseBuffer(const MH2 &rhs) {
    buff_values.insert(buff_values.end(), rhs.buff_values.begin(), rhs.buff_values.end());
    buff_weights.insert(buff_weights.end(), rhs.buff_weights.begin(), rhs.buff_weights.end());
  }

  // Keep it public for buffer fusion
  std::vector<std::vector<double>> buff_values;
  std::vector<double>              buff_weights;

 private:
  std::string name;  // Histogram name

  // -----------------------------------------------------------
  // For autorange

  bool              FILLBUFF     = false;
  int               AUTOBUFFSIZE = 100000;  // Default AUTOBUFFSIZE
  std::vector<bool> AUTOSYMMETRY = {false, false};
  // -----------------------------------------------------------

  // Boundary conditions
  double XMIN  = 0.0;
  double XMAX  = 0.0;
  int    XBINS = 0;

  double YMIN  = 0.0;
  double YMAX  = 0.0;
  int    YBINS = 0;

  // Number of underflow and overflow counts
  long long int              fills     = 0;
  std::vector<long long int> overflow  = {0, 0};
  std::vector<long long int> underflow = {0, 0};
  long long int              nanflow   = 0;

  // Logarithmic binning
  bool LOGX = false;
  bool LOGY = false;
  bool ValidBin(int xbin, int ybin) const;
  int  GetIdx(double value, double minval, double maxval, int nbins, bool logbins) const;

  // Weights (in unweighted case weights = counts)
  MMatrix<double> weights;
  MMatrix<double> weights2;

  // Counts
  MMatrix<long long int> counts;
};

}  // namespace gra

#endif