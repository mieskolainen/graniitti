// Eikonal proton density and Screening amplitude class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEIKONAL_H
#define MEIKONAL_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MPDG.h"

// External
#include "rang.hpp"

namespace gra {

// Numerical integration parameters
struct MEikonalNumerics {

  static constexpr double    MinKT2 = 1E-6;
  static constexpr double    MaxKT2 = 25.0;
  unsigned int NumberKT2 = 0;
  bool         logKT2    = false;

  static constexpr double    MinBT  = 1E-6;
  static constexpr double    MaxBT  = 10.0 / PDG::GeV2fm;
  unsigned int NumberBT  = 0;
  bool         logBT     = false;

  static constexpr double       FBIntegralMinKT = 1E-9;
  static constexpr double       FBIntegralMaxKT = 30.0;
  static constexpr unsigned int FBIntegralN     = 10000;
  static constexpr double       MinLoopKT       = 1E-4;

  double          MaxLoopKT  = 1.75;

  unsigned int NumberLoopKT  = 15;  // Number of kt steps  (default minimum)
  unsigned int NumberLoopPHI = 12;  // Number of phi steps (default minimum)

  // User setup (ND can be negative, to get below the default)
  void SetLoopDiscretization(int ND) {
    NumberLoopKT  = std::max(3, 3 * ND + (int)NumberLoopKT);
    NumberLoopPHI = std::max(3, 3 * ND + (int)NumberLoopPHI);
  
    if (ND < 0) { std::cout << rang::fg::red; }
    std::cout << "MEikonalNumerics::SetLoopDiscretization: ND = " << ND << std::endl;
    std::cout << "-> NumberLoopKT  = " << NumberLoopKT  << std::endl;
    std::cout << "-> NumberLoopPHI = " << NumberLoopPHI << std::endl;
    std::cout << rang::fg::reset;
  }

  // Unique hash
  std::string GetHashString() {
    std::string str =
        std::to_string(MinKT2) + std::to_string(MaxKT2) +
        std::to_string(NumberKT2) + std::to_string(logKT2) +
        std::to_string(MinBT) + std::to_string(MaxBT) +
        std::to_string(NumberBT) + std::to_string(logBT) +
        std::to_string(FBIntegralMinKT) +
        std::to_string(FBIntegralMaxKT) +
        std::to_string(FBIntegralN);
    return str;
  }

  // Read parameters from file
  void ReadParameters() {

    // Read and parse
    using json = nlohmann::json;

    const std::string inputfile = gra::aux::GetBasePath(2) + "/modeldata/" + "NUMERICS.json";
    const std::string data      = gra::aux::GetInputData(inputfile);
    json              j;

    try {
      j = json::parse(data);

      // JSON block identifier
      const std::string XID = "NUMERICS_EIKONAL";

      NumberKT2 = j.at(XID).at("NumberKT2");
      logKT2    = j.at(XID).at("logKT2");
      NumberBT  = j.at(XID).at("NumberBT");
      logBT     = j.at(XID).at("logBT");

    } catch (...) {
      std::string str = "ReadParameters: Error parsing " + inputfile +
                        " (Check for extra/missing commas)";
      throw std::invalid_argument(str);
    }
  }

};


// Interpolation container
class IArray1D {
 public:
  IArray1D(){};
  ~IArray1D(){};

  std::string  name;
  double       MIN   = 0;
  double       MAX   = 0;
  unsigned int N     = 0;
  double       STEP  = 0;
  bool         islog = false;

  // Setup discretization
  void Set(std::string _name, double _min, double _max, double _N, double _logarithmic) {
    // AT least two intervals
    if (_N < 2) {
      throw std::invalid_argument("IArray1D::Set: Error: N = " + std::to_string(_N) + " < 2");
    }

    // Check sanity
    if (_min >= _max) {
      throw std::invalid_argument("IArray1D::Set: Error: Variable " + _name + " MIN = " +
                                  std::to_string(_min) + " >= MAX = " + std::to_string(_max));
    }

    // Check sanity
    if (_logarithmic && _min < 1e-9) {
      throw std::invalid_argument("IArray1D::Set: Error: Variable " + _name +
                                  " is using logarithmic stepping with boundary MIN = " +
                                  std::to_string(_min));
    }

    islog = _logarithmic;

    // Logarithm taken here!
    MIN = islog ? std::log(_min) : _min;
    MAX = islog ? std::log(_max) : _max;
    N   = _N;

    STEP = (MAX - MIN) / N;
    name = _name;
  }

  // Call this last
  void InitArray() {
    // Note N+1 !
    F = MMatrix<std::complex<double>>(N + 1, 2, 0.0);
  }

  // N+1 size!
  MMatrix<std::complex<double>> F;

  // CMS energy (needed for boundary conditions)
  double sqrts = 0.0;

  std::string GetHashString() const {
    std::string str = std::to_string(islog) + std::to_string(MIN) + std::to_string(MAX) +
                      std::to_string(N) + std::to_string(sqrts);
    return str;
  }

  bool WriteArray(const std::string &filename, bool overwrite) const;
  bool ReadArray(const std::string &filename);
  std::complex<double> Interpolate1D(double a) const;

  static const unsigned int X = 0;
  static const unsigned int Y = 1;
};

class MEikonal {
 public:
  MEikonal();
  ~MEikonal();

  // Construct outside
  void S3Constructor(double s_in, const std::vector<gra::MParticle> &initialstate_in,
                     bool onlydensity = false, int NumberBT = 0, int NumberKT2 = 0);
  
  // Get eikonal and amplitude values
  // std::complex<double> S3DensityInterpolator(double bt);
  // std::complex<double> S3ScreeningInterpolator(double kt2);

  // Initialization already done
  bool IsInitialized() const {
    if (S3INIT == true) { return true; }
    return false;
  }

  // Get total cross sections
  void GetTotXS(double &tot, double &el, double &in) const;

  static std::complex<double> SingleAmpElastic(double s, double t, int type);
  std::complex<double> S3Density(double bt) const;
  std::complex<double> S3Screening(double kt2) const;

  // Get random number of cut Pomerons
  template <typename T>
  void S3GetRandomCutsBt(unsigned int &m, double &bt, T &rng) {
    const double STEP = (Numerics.MaxBT - Numerics.MinBT) / Numerics.NumberBT;

    // Numerical integral loop over impact parameter (b_t) space
    // C++11, thread_local is also static
    thread_local std::uniform_real_distribution<double>      flat(0, 1);
    thread_local std::uniform_int_distribution<unsigned int> randbt(0, Numerics.NumberBT);
    thread_local std::uniform_int_distribution<unsigned int> randm(1, MCUT - 1);  // 1,2,...

    // Acceptance-Rejection
    unsigned int n = 0;
    while (true) {
      // Draw random impact parameter flat
      n = randbt(rng);

      const int value = randm(rng);
      if (flat(rng) < P_array[value][n]) {
        m = value;  // m cut Pomerons
        break;
      }
    }
    bt = Numerics.MinBT + n * STEP;
  }

  // Get random number of cut Pomerons
  template <typename T>
  unsigned int S3GetRandomCuts(T &rng) {
    // Random integer from [1,NBins-1]
    // C++11, thread_local is also static
    thread_local std::uniform_real_distribution<double>      flat(0, 1);
    thread_local std::uniform_int_distribution<unsigned int> RANDI(1, P_cut.size() - 1);
    
    // Acceptance-Rejection
    while (true) {
      const unsigned int m = RANDI(rng);
      if (flat(rng) < P_cut[m]) {
        return m;  // m cut Pomerons
      }
    }
  }

  // Arrays
  IArray1D MBT;
  IArray1D MSA;

  // Parameters
  MEikonalNumerics Numerics;
  
 private:
  static const unsigned int MCUT = 25;  // Maximum number of cut Pomerons

  // S3 soft survival initialized
  bool S3INIT = false;

  void S3CalculateArray(IArray1D &add, std::complex<double> (MEikonal::*f)(double) const);
  void S3CalcXS();
  void S3InitCutPomerons();

  // Mandelstam s
  double s = 0.0;

  // Initial state
  std::vector<gra::MParticle> INITIALSTATE;

  // Integrated eikonal based total cross sections
  double sigma_tot  = 0.0;
  double sigma_el   = 0.0;
  double sigma_inel = 0.0;

  // Cut Pomeron probabilities
  std::vector<double>              P_cut;
  std::vector<std::vector<double>> P_array;

  // Two-Channel eikonal based cross sections
  double sigma_diff[4] = {0.0};

  static const unsigned int X = 0;
  static const unsigned int Y = 1;
};

}  // gra namespace ends

#endif