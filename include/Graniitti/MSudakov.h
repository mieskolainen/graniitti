// Shuvaev PDF and Sudakov suppression class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSUDAKOV_H
#define MSUDAKOV_H

// C++
#include <complex>
#include <memory>
#include <vector>

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MTimer.h"

// Libraries
#include "json.hpp"

namespace gra {

struct MSudakovNumerics {
  double       q2_MIN = 0.0;
  double       q2_MAX = 0.0;
  double       M_MIN  = 0.0;
  double       M_MAX  = 0.0;
  double       x_MIN  = 0.0;
  const double x_MAX  = 1.0 - 1E-9;

  // Numerical integral discretization [SET HERE]
  const unsigned int SudakovIntegralN = 1000;
  const unsigned int ShuvaevIntegralN = 400;

  // Logarithmic stepping true/false
  std::vector<bool> SUDA_log_ON;
  std::vector<bool> SHUV_log_ON;

  // Number of discrete intervals
  std::vector<unsigned int> SUDA_N;
  std::vector<unsigned int> SHUV_N;

  // CMS energy
  double sqrts = 0.0;

  bool DEBUG = false;

  void ReadParameters() {
    // Read and parse
    using json                  = nlohmann::json;
    const std::string inputfile = gra::aux::GetBasePath(2) + "/modeldata/" + "NUMERICS.json";
    const std::string data      = gra::aux::GetInputData(inputfile);
    json              j;

    try {
      j = json::parse(data);

      // JSON block identifier
      const std::string XID = "NUMERICS_SUDAKOV";

      q2_MIN = j.at(XID).at("q2_MIN");
      q2_MAX = j.at(XID).at("q2_MAX");
      M_MIN  = j.at(XID).at("M_MIN");
      // M_MAX = j.at(XID).at("M_MAX"); // Post-Setup
      // x_MIN = j.at(XID).at("x_MIN"); // Post-Setup
      // x_MAX = j.at(XID).at("x_MAX"); // Set above

      // Logarithmic stepping true/false [assign needed for <cast>]
      std::vector<bool> xx = j.at(XID).at("SUDA").at("log_ON");
      SUDA_log_ON          = xx;
      std::vector<bool> yy = j.at(XID).at("SHUV").at("log_ON");
      SHUV_log_ON          = yy;

      // Number of node points [assign needed for <cast>]
      std::vector<unsigned int> nn = j.at(XID).at("SUDA").at("N");
      SUDA_N                       = nn;
      std::vector<unsigned int> mm = j.at(XID).at("SHUV").at("N");
      SHUV_N                       = mm;

      // Discretization
      // ShuvaevIntegralN = j.at(XID).at("ShuvaevIntegralN");
      // SudakovIntegralN = j.at(XID).at("SudakovIntegralN");

      DEBUG = j.at(XID).at("DEBUG");
    } catch (...) {
      std::string str = "MSudakovNumerics::ReadParameters: Error parsing " + inputfile +
                        " (Check for extra/missing commas)";
      throw std::invalid_argument(str);
    }
  }
};


// Interpolation container
class IArray2D {
 public:
  IArray2D(){};
  ~IArray2D(){};

  std::string  name[2];
  double       MIN[2]   = {0};
  double       MAX[2]   = {0};
  unsigned int N[2]     = {0};
  double       STEP[2]  = {0};
  bool         islog[2] = {false};

  // Setup discretization
  void Set(unsigned int VAR, std::string _name, double _min, double _max, double _N,
           double _logarithmic) {
    // Out of index
    if (VAR > 1) {
      throw std::invalid_argument("IArray2D::Set: Error: VAR = " + std::to_string(VAR) + " > 1");
    }

    // AT least two intervals
    if (_N < 2) {
      throw std::invalid_argument("IArray2D::Set: Error: N = " + std::to_string(_N) + " < 2");
    }

    // Sanity
    if (_min >= _max) {
      throw std::invalid_argument("IArray2D::Set: Error: Variable " + _name + " MIN = " +
                                  std::to_string(_min) + " >= MAX = " + std::to_string(_max));
    }

    // Sanity
    if (_logarithmic && _min < 1e-9) {
      throw std::invalid_argument(
          "IArray2D::Set: Error: Variable " + _name +
          " is using logarithmic stepping with boundary MIN = " + std::to_string(_min));
    }
    islog[VAR] = _logarithmic;

    // Logarithm taken here!
    MIN[VAR] = islog[VAR] ? std::log(_min) : _min;
    MAX[VAR] = islog[VAR] ? std::log(_max) : _max;
    N[VAR]   = _N;

    STEP[VAR] = (MAX[VAR] - MIN[VAR]) / N[VAR];
    name[VAR] = _name;
  }

  // Call this last
  void InitArray() {
    // Note N+1 !
    F = std::vector<std::vector<std::vector<double>>>(
        N[0] + 1, std::vector<std::vector<double>>(N[1] + 1, std::vector<double>(4, 0.0)));
  }

  // N+1 per dimension!
  std::vector<std::vector<std::vector<double>>> F;

  // CMS energy (needed for boundary conditions)
  double sqrts = 0.0;

  std::string GetHashString() const {
    std::string str = std::to_string(islog[0]) + std::to_string(islog[1]) + std::to_string(MIN[0]) +
                      std::to_string(MIN[1]) + std::to_string(MAX[0]) + std::to_string(MAX[1]) +
                      std::to_string(N[0]) + std::to_string(N[1]) + std::to_string(sqrts);
    return str;
  }

  bool                      WriteArray(const std::string &filename, bool overwrite) const;
  bool                      ReadArray(const std::string &filename);
  std::pair<double, double> Interpolate2D(double A, double B) const;
};

// Sudakov suppression and skewed pdf
//
// Take care when copying this class - there is a pointer to LHAPDF
class MSudakov {
 public:
  MSudakov();
  ~MSudakov();

  void Init(double sqrts_in, const std::string &PDFSET, bool init_arrays = true);
  void InitLHAPDF(const std::string &PDFSET);
  std::pair<double, double> Shuvaev_H(double q2, double x);
  std::pair<double, double> Sudakov_T(double qt2, double M);

  double fg_xQ2M(double x, double q2, double M) const;
  double AlphaS_Q2(double q2) const;
  double NumFlavor(double q2) const;
  double xg_xQ2(double x, double Q2) const;
  void   TestPDF() const;

  bool initialized = false;

 private:
  int init_trials = 0;

  void   InitArrays();
  double diff_xg_xQ2_wrt_Q2(double x, double q2) const;
  double AP_gg(double delta) const;
  double AP_qg(double delta, double qt2) const;
  void   CalculateArray(IArray2D &arr, std::pair<double, double> (MSudakov::*f)(double, double));

  std::string  PDFSETNAME;
  LHAPDF::PDF *PdfPtr = nullptr;

  IArray2D veto;  // Sudakov veto
  IArray2D spdf;  // Shuvaev pdf

  MSudakovNumerics Numerics;  // Numerics
};

}  // namespace gra

#endif