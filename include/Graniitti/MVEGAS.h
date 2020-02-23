// VEGAS integrator class grid routines
//
// Future step:
// Factorize all VEGAS functions out from MGraniitti to here (use function pointers etc.)
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MVEGAS_H
#define MVEGAS_H

#include <vector>


namespace gra {

// Vegas MC default parameters
struct VEGASPARAM {
  unsigned int BINS   = 128;  // Maximum number of bins per dimension (EVEN NUMBER!)
  double       LAMBDA = 1.5;  // Regularization parameter

  // Initialization
  unsigned int NCALL     = 20000;  // Number of calls per iteration
  unsigned int ITER      = 15;     // Number of iterations
  double       CHI2MAX   = 10.0;   // Maximum chi2 in initialization
  double       PRECISION = 0.01;   // Maximum relative error of cross section integral
  int          DEBUG     = -1;     // Debug mode

  // User cannot set these
  unsigned int MAXFDIM = 100;      // Maximum integral dimension
  double       EPS     = 1.0e-30;  // Epsilon parameter
};


// Vegas MC adaptation data
struct VEGASData {
  // VEGAS initialization function
  void Init(unsigned int init, const VEGASPARAM &param) {
    // First initialization: Create the grid and initial data
    if (init == 0) {
      ClearAll(param);
      InitGridDependent(param);
    }
    // Use the previous grid but NOT integral data
    else if (init == 1) {
      ClearIntegral();
      InitGridDependent(param);
    }
    // Use the previous grid and its integral data
    else if (init == 2) {
      InitGridDependent(param);
    } else {
      throw std::invalid_argument("VEGASData::Init: Unknown init parameter = " +
                                  std::to_string(init));
    }
  }

  // Create number of calls per thread, they need to sum to calls
  std::vector<unsigned int> GetLocalCalls(int calls, int N_threads) {
    std::vector<unsigned int> localcalls(N_threads, 0.0);
    int                       sum = 0;
    for (int k = 0; k < N_threads; ++k) {
      localcalls[k] = std::floor(calls / N_threads);
      sum += localcalls[k];
    }
    // Add remainder to the thread number[0]
    localcalls[0] += calls - sum;
    return localcalls;
  }

  // Initialize
  void InitGridDependent(const VEGASPARAM &param) {
    // Create grid spacing
    for (std::size_t j = 0; j < FDIM; ++j) { dxvec[j] = region[j + FDIM] - region[j]; }

    // If binning parameter changed from previous call
    if (param.BINS != BINS_prev) {
      for (std::size_t i = 0; i < std::max(param.BINS, BINS_prev); ++i) { rvec[i] = 1.0; }
      for (std::size_t j = 0; j < FDIM; ++j) {
        Rebin(BINS_prev / static_cast<double>(param.BINS), j, param);
      }
      BINS_prev = param.BINS;
    }
  }

  // VEGAS rebinning function (algorithm adapted from Numerical Recipes)
  void Rebin(double ac, unsigned int j, const VEGASPARAM &param) {
    unsigned int k  = 0;
    double       dr = 0.0;
    double       zn = 0.0;
    double       zo = 0.0;

    for (std::size_t i = 0; i < param.BINS - 1; ++i) {
      while (ac > dr) { dr += rvec[(++k) - 1]; }
      if (k > 1) { zo = xmat[k - 2][j]; }
      zn = xmat[k - 1][j];
      dr -= ac;
      xcache[i] = zn - (zn - zo) * dr / rvec[k - 1];
    }

    for (std::size_t i = 0; i < param.BINS - 1; ++i) { xmat[i][j] = xcache[i]; }

    xmat[param.BINS - 1][j] = 1.0;
  }

  // VEGAS grid optimizing function (algorithm adapted from Numerical Recipes)
  void OptimizeGrid(const VEGASPARAM &param) {
    double zo = 0;
    double zn = 0;
    double ac = 0;

    for (std::size_t j = 0; j < FDIM; ++j) {
      zo          = f2mat[0][j];
      zn          = f2mat[1][j];
      f2mat[0][j] = (zo + zn) / 2.0;
      dcache[j]   = f2mat[0][j];

      for (std::size_t i = 2; i < param.BINS; ++i) {
        ac              = zo + zn;
        zo              = zn;
        zn              = f2mat[i][j];
        f2mat[i - 1][j] = (ac + zn) / 3.0;
        dcache[j] += f2mat[i - 1][j];
      }
      f2mat[param.BINS - 1][j] = (zo + zn) / 2.0;
      dcache[j] += f2mat[param.BINS - 1][j];
    }

    for (std::size_t j = 0; j < FDIM; ++j) {
      ac = 0.0;
      for (std::size_t i = 0; i < param.BINS; ++i) {
        f2mat[i][j] = (f2mat[i][j] < param.EPS) ? param.EPS : f2mat[i][j];
        rvec[i]     = std::pow((1.0 - f2mat[i][j] / dcache[j]) /
                               (std::log(dcache[j]) - std::log(f2mat[i][j]) + param.EPS),
                           param.LAMBDA);
        // Floating point protection (integrand close to 0)
        if (std::isnan(rvec[i]) || std::isinf(rvec[i])) { rvec[i] = param.EPS; }
        ac += rvec[i];
      }
      Rebin(ac / param.BINS, j, param);
    }
  }

  // Initialize sampling region [0,1] x [0,1] x ... x [0,1]
  void InitRegion(unsigned int fdim) {
    FDIM = fdim;
    region.resize(2 * FDIM, 0.0);

    // Lower bound [zero]
    for (std::size_t i = 0; i < FDIM; ++i) { region[i] = 0.0; }
    // Upper bound [one]
    for (std::size_t i = FDIM; i < 2 * FDIM; ++i) { region[i] = 1.0; }
  }

  // Clear integrated data
  void ClearIntegral() {
    sumdata = 0.0;
    sumchi2 = 0.0;
    sweight = 0.0;
  }

  // Full init
  void ClearAll(const VEGASPARAM &param) {
    // Vectors [MAXFDIM]
    dcache = std::vector<double>(param.MAXFDIM, 0.0);
    dxvec  = std::vector<double>(param.MAXFDIM, 0.0);

    // Vectors [BINS]
    rvec   = std::vector<double>(param.BINS, 0.0);
    xcache = std::vector<double>(param.BINS, 0.0);

    // Matrices [BINS x MAXFDIM]
    fmat  = std::vector<std::vector<double>>(param.BINS, std::vector<double>(param.MAXFDIM, 0.0));
    f2mat = std::vector<std::vector<double>>(param.BINS, std::vector<double>(param.MAXFDIM, 0.0));
    xmat  = std::vector<std::vector<double>>(param.BINS, std::vector<double>(param.MAXFDIM, 0.0));

    // Init with 1!
    for (std::size_t j = 0; j < FDIM; ++j) { xmat[0][j] = 1.0; }

    // VEGAS integraldata
    sumdata = 0.0;
    sumchi2 = 0.0;
    sweight = 0.0;

    // Previous binning
    BINS_prev = 1;
  }

  // VEGAS scalars
  unsigned int BINS_prev = 0;
  unsigned int FDIM      = 0;

  double fsum  = 0.0;
  double f2sum = 0.0;

  double sumdata = 0.0;
  double sumchi2 = 0.0;
  double sweight = 0.0;

  // Vectors
  std::vector<double> region;

  // Vectors
  std::vector<double> dcache;
  std::vector<double> dxvec;

  // Vectors
  std::vector<double> rvec;
  std::vector<double> xcache;

  // Matrices
  std::vector<std::vector<double>> fmat;
  std::vector<std::vector<double>> f2mat;
  std::vector<std::vector<double>> xmat;

};  // struct VEGASData

}  // namespace gra

#endif