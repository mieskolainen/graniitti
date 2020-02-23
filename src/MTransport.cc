// Optimal Transport Algorithms
//
//
// (c) 2017-2020 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++ standard
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// Own
#include "Graniitti/MMatOper.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MTransport.h"

// Eigen
#include <Eigen/Dense>

namespace gra {
namespace opt {


// Input:  n, m,    Kernel dimensions
//         lambda,  Entropic regularization
//
// Output: K (n x m) Convolution kernel matrix
//
void ConvKernel(std::size_t n, std::size_t m, double lambda, MMatrix<double>& K) {
  std::vector<double> x = math::linspace(0.0, n - 1.0, n);
  matoper::ScaleVector(x, 1.0 / (n - 1.0));

  std::vector<double> y = math::linspace(0.0, m - 1.0, m);
  matoper::ScaleVector(y, 1.0 / (m - 1.0));

  MMatrix<double> X;
  MMatrix<double> Y;
  matoper::MeshGrid(x, y, X, Y);

  // Convolution matrix
  K = MMatrix<double>(X.size_row(), X.size_col());
  for (std::size_t i = 0; i < X.size_row(); ++i) {
    for (std::size_t j = 0; j < X.size_col(); ++j) {
      K[i][j] = std::exp(-math::pow2(X[i][j] - Y[i][j]) / lambda);
    }
  }
}

// Input:  lambda,    Entropic regularization
//         C (n x m), Cost matrix
//
// Output: K (n x m) Gibbs (Gaussian convolution) kernel matrix
//
void GibbsKernel(double lambda, const MMatrix<double>& C, MMatrix<double>& K) {
  // Gibbs Kernel: Entropy regularized distance matrix
  K = MMatrix<double>(C.size_row(), C.size_col());
  for (std::size_t i = 0; i < C.size_row(); ++i) {
    for (std::size_t j = 0; j < C.size_col(); ++j) { K[i][j] = std::exp(-C[i][j] / lambda); }
  }
}

// Sinkhorn-Knopp non-linear but convex optimization algorithm
//
// Input:
//
// Kernel matrix between elements of p and q              (n x m)
// Probability density p (histogram / dirac point mass)   (n x 1)   with \sum = 1
// Probability density q (histogram / dirac point mass)   (m x 1)   with \sum = 1
// Iterations                                             (scalar)
//
// Output:
//
// Optimal Transport Matrix Pi                            (n x m)
// Distance measure                                       (scalar)
//
double SinkHorn(MMatrix<double>& Pi, const MMatrix<double>& K, std::vector<double>& p,
                std::vector<double>& q, std::size_t iter) {
  const double SAFE_EPS = 1e-64;

  std::cout << "opt::SinkHorn optimization:" << std::endl;

  // ============================================================
  // Test the normalization
  auto fsum = [](const std::vector<double>& x) {
    double value = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) { value += x[i]; }
    return value;
  };

  const double EPS = 1e-5;
  if (std::abs(fsum(p) - 1) > EPS) {
    const std::string str = "SinkHorn:: Input p elements sum != 1";
    throw std::invalid_argument(str);
  }
  if (std::abs(fsum(q) - 1) > EPS) {
    const std::string str = "SinkHorn:: Input q elements sum != 1";
    throw std::invalid_argument(str);
  }
  // ============================================================

  const std::size_t n = K.size_row();
  const std::size_t m = K.size_col();

  // Transposed version
  MMatrix<double> KT = K.Transpose();

  // Initialization
  std::vector<double> u(n, 0.0);  // no need for initialization
  std::vector<double> v(m, 1.0);  // I

  // || Pi * 1 - p || or || Pi^T 1 - q ||
  std::vector<double> ones_p(m, 1.0);
  std::vector<double> ones_q(n, 1.0);

  auto resnorm_p = [&](const MMatrix<double>& Pi) {
    double                    value = 0.0;
    const std::vector<double> lhs   = Pi * ones_p;
    for (std::size_t i = 0; i < lhs.size(); ++i) { value += std::abs(lhs[i] - p[i]); }
    return value;
  };

  auto resnorm_q = [&](const MMatrix<double>& PiT) {
    double                    value = 0.0;
    const std::vector<double> lhs   = PiT * ones_q;
    for (std::size_t i = 0; i < lhs.size(); ++i) { value += std::abs(lhs[i] - q[i]); }
    return value;
  };

  // Sinkhorn iterations
  double p_cost = 0.0;
  double q_cost = 0.0;


  MMatrix<double> PiT;

  for (std::size_t k = 0; k < iter; ++k) {
    const std::vector<double> Kv = K * v;  // >= 0
    for (std::size_t i = 0; i < n; ++i) { u[i] = p[i] / std::max(SAFE_EPS, Kv[i]); }

    const std::vector<double> KTu = KT * u;  // >= 0
    for (std::size_t i = 0; i < m; ++i) { v[i] = q[i] / std::max(SAFE_EPS, KTu[i]); }

    // Calculate metrics
    if ((k + 1) % (iter / 10) == 0 || k == 0 || k == iter - 1) {
      // Transport coupling matrix
      Pi  = matoper::diagAdiag(u, K, v);
      PiT = Pi.Transpose();

      // Transport cost
      p_cost = resnorm_p(Pi);
      q_cost = resnorm_q(PiT);

      printf(
          "iter = %4lu / %4lu : p_cost = %0.5E, log(p_cost) = %4.1f, q_cost = %0.5E, log(q_cost) = "
          "%4.1f \n",
          k + 1, iter, p_cost, std::log(p_cost), q_cost, std::log(q_cost));
    }
  }

  return p_cost;
}

}  // namespace opt
}  // namespace gra
