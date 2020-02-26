// Shuvaev PDF and Sudakov suppression class
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <future>
#include <iostream>
#include <memory>
#include <vector>
#include <fstream>

// C file processing
#include <sys/stat.h>
#include <sys/types.h>
#include <unistd.h>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MSudakov.h"
#include "Graniitti/MTimer.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

// Libraries
#include "rang.hpp"

namespace gra {

using aux::indices;
using math::msqrt;
using math::pow2;
using math::zi;

// Constructor
MSudakov::MSudakov() {}

// Destructor
MSudakov::~MSudakov() { delete PdfPtr; }

// Init
void MSudakov::Init(double _sqrts, const std::string &PDFSET, bool init_arrays) {
  // First this
  Numerics.ReadParameters();

  // Setup energy dependent boundaries
  Numerics.sqrts = _sqrts;
  Numerics.M_MAX = _sqrts;

  // Collinear case: M^2 = x_1x_2s, let x_1 = 1.0 gives minimum x_2 = M^2/s
  Numerics.x_MIN = math::pow2(Numerics.M_MIN / Numerics.sqrts);

  // InitLHAPDF("CT10nlo");
  InitLHAPDF(PDFSET);

  if (init_arrays) { InitArrays(); }
}

void MSudakov::InitArrays() {
  // Init Sudakov variables
  {
    std::cout << "Initializing <Sudakov> array:" << std::endl;

    // q2,M
    veto.sqrts = Numerics.sqrts;  // FIRST THIS
    veto.Set(0, "q2", Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SUDA_N[0],
             Numerics.SUDA_log_ON[0]);
    veto.Set(1, "M", Numerics.M_MIN, Numerics.M_MAX, Numerics.SUDA_N[1], Numerics.SUDA_log_ON[1]);
    veto.InitArray();  // Initialize (call last!)

    const unsigned long hash     = gra::aux::djb2hash(veto.GetHashString());
    const std::string   filename = gra::aux::GetBasePath(2) + "/eikonal/SUDA_" + gra::aux::ToString(veto.sqrts, 0) + "_" +
                                 PDFSETNAME + "_" + std::to_string(hash);

    // Try to read pre-calculated
    bool ok = veto.ReadArray(filename);

    while (!ok) {  // Problem, re-calculate
      // Pointer to member function: ReturnType
      // (ClassType::*)(ParameterTypes...)
      std::pair<double, double> (MSudakov::*f)(double, double) = &MSudakov::Sudakov_T;
      CalculateArray(veto, f);
      veto.WriteArray(filename, true);
      ok = veto.ReadArray(filename);
    }
  }

  // Init Shuvaev PDF variables
  {
    std::cout << "Initializing <Shuvaev> array:" << std::endl;

    // q2,x
    spdf.sqrts = Numerics.sqrts;  // FIRST THIS
    spdf.Set(0, "q2", Numerics.q2_MIN, Numerics.q2_MAX, Numerics.SHUV_N[0],
             Numerics.SHUV_log_ON[0]);
    spdf.Set(1, "x", Numerics.x_MIN, Numerics.x_MAX, Numerics.SHUV_N[1], Numerics.SHUV_log_ON[1]);
    spdf.InitArray();  // Initialize (call last!)

    const unsigned long hash     = gra::aux::djb2hash(spdf.GetHashString());
    const std::string   filename = gra::aux::GetBasePath(2) + "/eikonal/SHUV_" + gra::aux::ToString(spdf.sqrts, 0) + "_" +
                                 PDFSETNAME + "_" + std::to_string(hash);

    // Try to read pre-calculated
    bool ok = spdf.ReadArray(filename);

    while (!ok) {  // Problem, re-calculate
      // Pointer to member function: ReturnType
      // (ClassType::*)(ParameterTypes...)
      std::pair<double, double> (MSudakov::*f)(double, double) = &MSudakov::Shuvaev_H;
      CalculateArray(spdf, f);
      spdf.WriteArray(filename, true);
      ok = spdf.ReadArray(filename);
    }
  }
  // --------------------------------------------------
  // Finally
  initialized = true;
  std::cout << std::endl;

  if (Numerics.DEBUG) { TestPDF(); }
}

// Return gluon xg(x,Q2) from LHAPDF
double MSudakov::xg_xQ2(double x, double q2) const {
  const int pid = 21;  // gluon
  try {
    return PdfPtr->xfxQ2(pid, x, q2);
  } catch (...) {
    std::string str =
        "MSudakov::xg_xQ2: Problem with x = " + std::to_string(x) + ", q2 = " + std::to_string(q2);
    throw std::invalid_argument(str);
  }
}

// PDF access
//
// [REFERENCE: LHAPDF6, https://arxiv.org/abs/1412.7420]
void MSudakov::InitLHAPDF(const std::string &pdfname) {
  PDFSETNAME = pdfname;

  // LHAPDF init
  try {
    PdfPtr = LHAPDF::mkPDF(pdfname, 0);
  } catch (...) {
    ++init_trials;

    std::string str = "MSudakov::InitLHAPDF: Trials = " + std::to_string(init_trials) +
                      " :: Problem with reading a pdfset '" + pdfname + "'";

    if (init_trials >= 2) {  // Too many failures
      throw std::invalid_argument(str);
    }

    std::cout << str << std::endl;

    aux::AutoDownloadLHAPDF(pdfname);
    InitLHAPDF(pdfname);  // try again
  }
}

// PDF test routine
void MSudakov::TestPDF() const {
  const double MINLOGX = std::log10(Numerics.x_MIN);
  const double MAXLOGX = std::log10(Numerics.x_MAX);
  const int    NX      = 5;  // Number of points - 1
  const double stepX   = (MAXLOGX - MINLOGX) / NX;

  const double MINLOGQ2 = std::log10(Numerics.q2_MIN);
  const double MAXLOGQ2 = std::log10(Numerics.q2_MAX);
  const int    NQ2      = 5;  // Number of points - 1
  const double stepQ2   = (MAXLOGQ2 - MINLOGQ2) / NQ2;

  const double MINLOGM = std::log10(Numerics.M_MIN);
  const double MAXLOGM = std::log10(Numerics.M_MAX);
  const int    NM      = 5;  // Number of points - 1
  const double stepM   = (MAXLOGM - MINLOGM) / NM;

  // Test loop
  for (std::size_t i = 0; i < NM + 1; ++i) {
    const double log10M = MINLOGM + i * stepM;
    const double M      = std::pow(10, log10M);

    printf("[M = %0.1f GeV] : alpha_s(Q = M GeV) = %0.3f \n\n", M, PdfPtr->alphasQ2(M * M));

    for (std::size_t j = 0; j < NX + 1; ++j) {
      const double log10x = MINLOGX + j * stepX;
      const double x      = std::pow(10, log10x);

      printf("x = %0.5E \n", x);
      for (std::size_t k = 0; k < NQ2 + 1; ++k) {
        const double log10q2 = MINLOGQ2 + k * stepQ2;
        const double q2      = std::pow(10, log10q2);

        // Normal gluon pdf
        const double xf = xg_xQ2(x, q2);

        // Durham flux
        const double hxf = fg_xQ2M(x, q2, M);

        printf(
            "(x = %0.3E, q2 = %0.2f, M = %0.1f) : [gluon pdf: xg(x,q2), "
            "Durham flux: fg(x,q2,M)] = (%0.2f,%0.2f) \n",
            x, q2, M, xf, hxf);
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }
}

// Access QCD coupling alpha_s(Q^2) from LHAPDF
double MSudakov::AlphaS_Q2(double q2) const {
  try {
    return PdfPtr->alphasQ2(q2);
  } catch (...) {
    std::string str = "MSudakov::AlphaS_Q2: Problem with q2 = " + std::to_string(q2);
    throw std::invalid_argument(str);
  }
}

// Calculate differentiation dxg/dQ2 via "Richardson's extrapolation"
// f'(x) = [4*D_0(h) + D_0(2h)] / 3 + O(h^4)
//
// Requires 4 evaluations of densities
//
// [REFERENCE: https://en.wikipedia.org/wiki/Richardson_extrapolation]
double MSudakov::diff_xg_xQ2_wrt_Q2(double x, double q2) const {
  // Do not take h less than 1E-4 with 64-bit double
  const double h   = 1E-4;
  const double hX2 = 2 * h;

  // Calculate central differences
  const double D0_A = (xg_xQ2(x, q2 + h) - xg_xQ2(x, q2 - h)) / (2 * h);
  const double D0_B = (xg_xQ2(x, q2 + hX2) - xg_xQ2(x, q2 - hX2)) / (2 * hX2);

  return (4.0 * D0_A - D0_B) / 3.0;
}

// Durham flux (skewed gluon pdf)
// ~ Shuvaev transformed gluon pdf x Sudakov suppression
//
// f_g(x,x',qt^2,\mu)
// = \frac{\partial}{\partial \ln Q_t^2} [H_g(x/2,x/2,Q_t^2)\sqrt{T(Q_t^2,\mu)}]
//
double MSudakov::fg_xQ2M(double x, double q2, double M) const {
  // Calculate Shuvaev transformation
  // std::pair<double,double> out1 = Shuvaev_H(q2, x);
  std::pair<double, double> out1 = spdf.Interpolate2D(q2, x);
  const double              Hg   = out1.first;
  const double              dHg  = out1.second;

  // Calculate Sudakov veto
  // std::pair<double,double> out2 = Sudakov_T(q2, M);
  std::pair<double, double> out2 = veto.Interpolate2D(q2, M);
  const double              Tg   = out2.first;
  const double              dTg  = out2.second;

  // Chain rule's: d/dln(q^2) [ ... ]
  double total = 0.0;

  if (Tg > 1e-15) {
    total = dHg * math::msqrt(Tg) + Hg * dTg / (2.0 * math::msqrt(Tg));
  } else {
    // printf("MSudakov::fg_xQ2M: Sudakov factor Tg = %0.2E [x = %0.1E, q2 =
    // %0.3f
    // GeV^2, M = %0.3f GeV] \n", Tg, x, q2, M);
  }

  // And because our functions return partial derivatives w.r.t. q^2 (not ln q^2)
  // we transform:
  // dln(q^2)/dq^2 = 1/q^2 <=> dln(q^2) = dq^2/q^2 =>
  // d[...]/dln(q^2) = d[...]/dq^2 * q^2, applied below:
  //
  total *= q2;

  // Truncate negative values to give proper zero
  total = (total < 0.0) ? 0.0 : total;

  return total;
}

// Calculate Shuvaev integral transform
// ----------------------------------------------------------------------
// Identity:
//
// H_g(x, \xi -> 0) = xg(x)
//
// ----------------------------------------------------------------------
//
// Hg  = numerical integral transformed from standard pdf
// dHg = dHg/dq^2 (differentiated numerically)
//
// [REFERENCE: Harland-Lang, https://arxiv.org/abs/1306.6661]
//
std::pair<double, double> MSudakov::Shuvaev_H(double q2, double x) {
  const double y_MIN  = x / 4.0;
  const double y_MAX  = 1.0;
  const double y_STEP = (y_MAX - y_MIN) / Numerics.ShuvaevIntegralN;

  double Hg  = 0.0;
  double dHg = 0.0;

  // Check that we are within valid domain (take into account floating points)
  const double EPS = 1e-5;
  if (x >= Numerics.x_MIN * (1 - EPS) && x <= Numerics.x_MAX * (1 + EPS) &&
      q2 >= Numerics.q2_MIN * (1 - EPS) && q2 <= Numerics.q2_MAX * (1 + EPS)) {
    // N+1!
    std::vector<double> fA(Numerics.ShuvaevIntegralN + 1, 0.0);
    std::vector<double> fB(Numerics.ShuvaevIntegralN + 1, 0.0);

    for (const auto &i : indices(fA)) {
      const double y        = y_MIN + i * y_STEP;
      const double argument = x / (4.0 * y);

      // H_g(x/2,x/2,Q^2) = 4x/\pi \int_{x/4}^1 dy
      // y^{1/2}(1-y)^{1/2}g(x/4y),q^2)
      //
      // => take into account that LHAPDF provides xg(), not g(), gives:
      const double factor = math::msqrt(math::pow3(y) * (1 - y));
      fA[i]               = factor * xg_xQ2(argument, q2);
      fB[i]               = factor * diff_xg_xQ2_wrt_Q2(argument, q2);
    }
    const double norm = 16.0 / math::PI;
    Hg                = norm * math::CSIntegral(fA, y_STEP);
    dHg               = norm * math::CSIntegral(fB, y_STEP);
  } else {
    // Fatal error
    throw std::invalid_argument("MSudakov::Shuvaev_H(q2,x): Input arguments out of domain: q2 = " +
                                std::to_string(q2) + ", x = " + std::to_string(x));
  }
  return {Hg, dHg};
}

// Calculate Sudakov real radiation suppression (veto) factor via direct
// numerical integral.
//
// T(Q_t^2, \mu) = exp( -\int_{Q_t^2}^{\mu^2} dk_t^2 \frac{1}{k_t^2}
// \frac{\alpha_s(k_t^2)}{2\pi}
//                       \int_0^{1-DELTA} dz[zP_gg(z) + \sum_q P_qg(z)] ),
//
// where DELTA = k_t/mu
//
// Tg  = numerical integral value [range [0,1)]
// dTg = dTg/dq^2 (analytic derivative)
//
// [REFERENCE: Coughlin, Forshaw, https://arxiv.org/abs/0912.3280v2]
//
std::pair<double, double> MSudakov::Sudakov_T(double qt2, double mu) {
  // Scale functor
  auto DeltaScale = [](double kt2, double M) -> double { return math::msqrt(kt2) / M; };

  // [OUTER INTEGRAL]: k_t^2 integral upper and lower bounds:
  const double MAX = std::log(math::pow2(mu));
  const double MIN = std::log(qt2);

  // Discretization
  const double STEP = (MAX - MIN) / Numerics.SudakovIntegralN;

  // [OUTER kt^2-integral]:
  // Change of a variable: u := ln(kt^2) [du = (du/dkt2) dkt2]
  // this cancels 1/kt^2 from the integrand and numerical convergence is
  // much improved.
  //
  // [INNER z-integral is analytic]:
  std::vector<double> f(Numerics.SudakovIntegralN + 1, 0.0);
  for (const auto &i : indices(f)) {
    const double u     = MIN + i * STEP;
    const double kt2   = std::exp(u);
    const double delta = DeltaScale(kt2, mu);
    f[i]               = AlphaS_Q2(kt2) * (AP_gg(delta) + AP_qg(delta, kt2));
  }
  const double integral = math::CSIntegral(f, STEP) / (2.0 * math::PI);

  // Final expression and its derivative
  const double Tg    = std::exp(-integral);
  const double delta = DeltaScale(qt2, mu);
  const double dTg =
      Tg * AlphaS_Q2(qt2) / (2.0 * math::PI * qt2) * (AP_gg(delta) + AP_qg(delta, qt2));
  return {Tg, dTg};
}

// Altarelli-Parisi splitting function definite integral over z
//
// see standard literature, e.g.
// [REFERENCE: B.R. Webber, CERN lectures 08,
// www.hep.phy.cam.ac.uk/theory/webber/QCDlect3.pdf]
// [REFERENCE: R.K. Ellis, W.J. Stirling, B.R. Webber, QCD and Collider Physics]
// [REFERENCE: http://pdg.lbl.gov/2018/download/db2018.pdf]
//
// \int z P_gg(z) dz
// \int z *{ 2 * CA * int [z/(1-z)_{+} + (1-z)/z + z*(1-z)] +
//           1/6*(11*CA - 4*Nf*TR)*deltafunc(1-z) } dz, z = 0 ... (1-delta) =
//           (below)
//
// First part: ...
// Second part: \int 1/6*(11*C - 4*N*T)*DiracDelta[1-x] dx, x = 0 ... 1 - delta,
// this vanishes with positive delta
//
// where + is the "plus-description"
//
double MSudakov::AP_gg(double delta) const {
  const double CA = 3.0;  // Structure constant
  return 2.0 * CA *
         (std::log(1.0 / delta) -
          math::pow2(1.0 - delta) * (3.0 * math::pow2(delta) - 2.0 * delta + 11.0) / 12.0);
}

// Altarelli-Parisi splitting function definite integral over z,
//
// \int \sum_q P_qz(z) dz
//
// with sum over active quark flavors with mass threshold qt2
// \sum_q TR * \int (z^2 + (1-z)^2) dz, z = 0 ... (1-delta)
//
double MSudakov::AP_qg(double delta, double qt2) const {
  const double TR = 0.5;  // Structure constant
  return TR * (-2.0 * math::pow3(delta) / 3.0 + math::pow2(delta) - delta + 2.0 / 3.0) *
         NumFlavor(qt2);
}

// Return the number of quark flavors at scale q^2
//
double MSudakov::NumFlavor(double q2) const {
  const double m_charm  = 1.275;  // Gev, PDG-2018 (default definition)
  const double m_bottom = 4.18;   // Gev

  if (q2 < math::pow2(m_charm)) {
    return 3.0;
  } else if (q2 < math::pow2(m_bottom)) {
    return 4.0;
  } else {
    return 5.0;
  }
}

// Constructs interpolation array values
void MSudakov::CalculateArray(IArray2D &arr,
                              std::pair<double, double> (MSudakov::*f)(double, double)) {
  MTimer timer;
  for (const auto &i : indices(arr.F)) {
    const double a = arr.MIN[0] + i * arr.STEP[0];

    // Transform input to linear if log stepping, for the function
    const double var1 = (arr.islog[0]) ? std::exp(a) : a;
    gra::aux::PrintProgress(i / static_cast<double>(arr.N[0] + 1));

    for (const auto &j : indices(arr.F[i])) {
      const double b = arr.MIN[1] + j * arr.STEP[1];

      // Transform input to linear if log stepping, for the function
      const double var2 = (arr.islog[1]) ? std::exp(b) : b;

      // Call function being pointed to
      const std::pair<double, double> output = (this->*f)(var1, var2);

      arr.F[i][j][0] = a;
      arr.F[i][j][1] = b;
      arr.F[i][j][2] = std::abs(output.first) < 1e-64 ? 0 : output.first;    // Underflow protection
      arr.F[i][j][3] = std::abs(output.second) < 1e-64 ? 0 : output.second;  //
    }
  }
  // Progressbar clearing
  gra::aux::ClearProgress();
  printf("MSudakov::CalculateArray: Time elapsed %0.1f sec \n", timer.ElapsedSec());
}

// Write the array to a file
bool IArray2D::WriteArray(const std::string &filename, bool overwrite) const {
  // Do not write if file exists already
  if (gra::aux::FileExist(filename) && !overwrite) {
    // std::cout << "- Found pre-calculated" << std::endl;
    return true;
  }
  // -------------------------------------------------------------

  std::ofstream file;
  file.open(filename);
  if (!file.is_open()) {
    std::string str = "IArray2D::WriteArray: Fatal IO-error with: " + filename;
    throw std::invalid_argument(str);
  }

  std::cout << "IArray2D::WriteArray: ";
  unsigned int line_number = 0;

  try {
    for (const auto &i : indices(F)) {
      for (const auto &j : indices(F[i])) {
        // Write to file
        file << std::setprecision(15) << F[i][j][0] << "," << F[i][j][1] << "," << F[i][j][2] << ","
             << F[i][j][3] << std::endl;

        ++line_number;
      }
    }
  } catch (...) {
    throw std::invalid_argument("IArray2D:WriteArray: Error in file " + filename + " at line " +
                                std::to_string(line_number));
  }

  file.close();
  return true;
}

// Read the array from a file
bool IArray2D::ReadArray(const std::string &filename) {
  std::ifstream file;
  file.open(filename);
  if (!file.is_open()) {
    std::string str = "IArray2D::ReadArray: Fatal IO-error with: " + filename;
    return false;
  }

  std::string  line;
  unsigned int fills       = 0;
  unsigned int line_number = 0;
  std::cout << "IArray2D::ReadArray: ";

  try {
    for (const auto &i : indices(F)) {
      for (const auto &j : indices(F[i])) {
        // Read every line from the stream
        getline(file, line);

        std::istringstream       stream(line);
        std::vector<std::string> columns;
        std::string              element;

        // Get every line element (4 of them) separated by separator
        int k = 0;
        while (getline(stream, element, ',')) {
          F[i][j][k] = std::stod(element);  // string to double
          ++k;
          ++fills;
        }
        ++line_number;
      }
    }

  } catch (...) {
    throw std::invalid_argument("IArray2D:ReadArray: Error in file " + filename + " at line " +
                                std::to_string(line_number));
  }

  file.close();

  if (fills != 4 * (N[0] + 1) * (N[1] + 1)) {
    std::string str = "Corrupted file: " + filename;
    std::cout << str << std::endl;
    return false;
  }
  std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;
  return true;
}

// Standard 2D-bilinear interpolation f(a,b) = Z, and derivative dZ
//
//
// [REFERENCE: https://en.wikipedia.org/wiki/Bilinear_interpolation]
std::pair<double, double> IArray2D::Interpolate2D(double a, double b) const {
  const double EPS = 1e-5;

  if (a < MIN[0]) { a = MIN[0]; }  // Truncate before (possible) logarithm
  if (b < MIN[1]) { b = MIN[1]; }  // Truncate before (possible) logarithm

  // Logarithmic stepping or not
  if (islog[0]) { a = std::log(a); }
  if (islog[1]) { b = std::log(b); }

  if (a > MAX[0] * (1 + EPS) || b > MAX[1] * (1 + EPS)) {
    printf(
        "Interpolate2D(%s,%s) Input out of grid domain: "
        "%s = %0.3f [%0.3f, %0.3f], %s = %0.3f [%0.3f, %0.3f] \n",
        name[0].c_str(), name[1].c_str(), name[0].c_str(), a, MIN[0], MAX[0], name[1].c_str(), b,
        MIN[1], MAX[1]);
  }

  // Get indices
  int i = std::floor((a - MIN[0]) / STEP[0]);
  int j = std::floor((b - MIN[1]) / STEP[1]);

  // Lower boundary protection
  if (i < 0) { i = 0; }  // Int needed for this (not unsigned int)
  if (j < 0) { j = 0; }

  // Upper boundary protection
  if (i >= (int)N[0]) { i = N[0] - 1; }  // We got N+1 elements in F
  if (j >= (int)N[1]) { j = N[1] - 1; }

  // Aux variables for readability
  const unsigned int X = 0;
  const unsigned int Y = 1;

  const double x1 = F[i][j][X];
  const double x2 = F[i + 1][j][X];
  const double y1 = F[i][j][Y];
  const double y2 = F[i][j + 1][Y];

  const double xstep = STEP[0];
  const double ystep = STEP[1];
  const double xval  = a;
  const double yval  = b;

  double values[2] = {0.0};

  // 2 == Z, 3 == dZ
  for (std::size_t C = 2; C <= 3; ++C) {
    const double Q11 = F[i][j][C];
    const double Q12 = F[i][j + 1][C];
    const double Q21 = F[i + 1][j][C];
    const double Q22 = F[i + 1][j + 1][C];

    // Interpolated valued
    values[C - 2] =
        ((y2 - yval) / ystep) * ((x2 - xval) / xstep * Q11 + (xval - x1) / xstep * Q21) +
        ((yval - y1) / ystep) * ((x2 - xval) / xstep * Q12 + (xval - x1) / xstep * Q22);
  }

  return {values[0], values[1]};
}

}  // namespace gra
