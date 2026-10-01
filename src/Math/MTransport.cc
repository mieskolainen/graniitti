// Entropy regularized optimal transport algorithms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MTransport.h"

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra {
namespace opt {

// Input:  n, m,    Kernel dimensions
//         lambda,  Entropic regularization
//
// Output: K (n x m) Convolution kernel matrix
// K_ij = exp[-(x_i-y_j)^2/lambda]
//
void ConvKernel(std::size_t n, std::size_t m, double lambda, MMatrix<double> &K) {
  if (n < 2 || m < 2 || !std::isfinite(lambda) || lambda <= 0.0) {
    throw std::invalid_argument("opt::ConvKernel: dimensions must be at least two and lambda positive");
  }
  std::vector<double> x = math::linspace(0.0, n - 1.0, n);
  gra::Scale(x, 1.0 / (n - 1.0));

  std::vector<double> y = math::linspace(0.0, m - 1.0, m);
  gra::Scale(y, 1.0 / (m - 1.0));

  auto [Y, X] = gra::MeshGrid(y, x);

  // Convolution matrix
  K = MMatrix<double>(X.size_row(), X.size_col());
  for (std::size_t i = 0; i < X.size_row(); ++i) {
    for (std::size_t j = 0; j < X.size_col(); ++j) { K[i][j] = std::exp(-math::pow2(X[i][j] - Y[i][j]) / lambda); }
  }
}

// Input:  lambda,    Entropic regularization
//         C (n x m), Cost matrix
//
// Output: K (n x m) Gibbs (Gaussian convolution) kernel matrix
// K_ij = exp(-C_ij/lambda)
//
void GibbsKernel(double lambda, const MMatrix<double> &C, MMatrix<double> &K) {
  if (!std::isfinite(lambda) || lambda <= 0.0 || C.size_row() == 0 || C.size_col() == 0 || !C.IsFinite()) {
    throw std::invalid_argument("opt::GibbsKernel: require a finite cost matrix and positive lambda");
  }
  for (std::size_t row = 0; row < C.size_row(); ++row) {
    for (std::size_t col = 0; col < C.size_col(); ++col) {
      if (C(row, col) < 0.0) { throw std::invalid_argument("opt::GibbsKernel: costs must be nonnegative"); }
    }
  }
  // Gibbs Kernel: Entropy regularized distance matrix
  K = C.Transform([lambda](const double value) { return std::exp(-value / lambda); });
}

// Sinkhorn-Knopp non-linear but convex optimization algorithm
// u <- p/(Kv), v <- q/(K^T u), Pi = diag(u) K diag(v)
//
// Input:
//
// Kernel matrix between elements of p and q              (n x m)
// Probability density p (histogram / dirac point mass)   (n x 1)   with \sum =
// 1 Probability density q (histogram / dirac point mass)   (m x 1)   with \sum
// = 1 Iterations                                             (scalar)
//
// Output:
//
// Optimal Transport Matrix Pi                            (n x m)
// Source-marginal L1 convergence residual                (scalar)
//
double SinkHorn(MMatrix<double> &Pi, const MMatrix<double> &K, const std::vector<double> &p,
                const std::vector<double> &q, std::size_t iter) {
  const std::size_t n        = K.size_row();
  const std::size_t m        = K.size_col();
  if (iter == 0 || n == 0 || m == 0 || p.size() != n || q.size() != m || !K.IsFinite()) {
    throw std::invalid_argument("opt::SinkHorn: invalid dimensions or input");
  }
  for (std::size_t row = 0; row < n; ++row) {
    for (std::size_t col = 0; col < m; ++col) {
      if (K(row, col) < 0.0) { throw std::invalid_argument("opt::SinkHorn: kernel elements must be nonnegative"); }
    }
  }
  const auto valid_probability = [](const double value) { return std::isfinite(value) && value >= 0.0; };
  if (!std::all_of(p.begin(), p.end(), valid_probability) || !std::all_of(q.begin(), q.end(), valid_probability)) {
    throw std::invalid_argument("opt::SinkHorn: probability elements must be finite and nonnegative");
  }

  // Test the normalization
  const double EPS = 1e-5;
  if (!std::isfinite(gra::Sum(p)) || std::abs(gra::Sum(p) - 1.0) > EPS) {
    const std::string str = "SinkHorn:: Input p elements sum != 1";
    throw std::invalid_argument(str);
  }
  if (!std::isfinite(gra::Sum(q)) || std::abs(gra::Sum(q) - 1.0) > EPS) {
    const std::string str = "SinkHorn:: Input q elements sum != 1";
    throw std::invalid_argument(str);
  }
  // ============================================================

  std::cout << "opt::SinkHorn optimization:" << std::endl;

  // Balance in logarithms so that every representable positive entry retains support
  const auto kernel = K.Transform([](const double value) { return std::log(static_cast<long double>(value)); });
  const auto KT     = kernel.Transpose();

  // Compute one marginal scaling without assigning artificial kernel support
  const auto ratio = [](const double probability, const long double denominator) {
    if (math::IsZero(probability)) { return -std::numeric_limits<long double>::infinity(); }
    if (!std::isfinite(denominator)) {
      throw std::runtime_error("opt::SinkHorn: marginal has no finite kernel support");
    }
    const long double value = std::log(static_cast<long double>(probability)) - denominator;
    if (!std::isfinite(value)) { throw std::runtime_error("opt::SinkHorn: marginal scaling overflow"); }
    return value;
  };

  // Initialization
  std::vector<long double> u(n, 0.0L);
  std::vector<long double> v(m, 0.0L);

  // || Pi * 1 - p || or || Pi^T 1 - q ||
  std::vector<double> ones_p(m, 1.0);
  std::vector<double> ones_q(n, 1.0);

  auto resnorm_p = [&](const MMatrix<double> &Pi) { return gra::L1Distance(Pi * ones_p, p); };

  auto resnorm_q = [&](const MMatrix<double> &PiT) { return gra::L1Distance(PiT * ones_q, q); };

  // Sinkhorn iterations
  double p_residual = 0.0;
  double q_residual = 0.0;

  MMatrix<double> PiT;

  const std::size_t report_stride = std::max<std::size_t>(1, iter / 10);
  for (std::size_t k = 0; k < iter; ++k) {
    for (const auto &i : indices(u)) { u[i] = ratio(p[i], gra::LogProductSum(kernel.Row(i), v)); }
    for (const auto &i : indices(v)) { v[i] = ratio(q[i], gra::LogProductSum(KT.Row(i), u)); }

    // Calculate metrics
    if ((k + 1) % report_stride == 0 || k == 0 || k == iter - 1) {
      // Transport coupling matrix
      Pi.Resize(n, m);
      for (std::size_t i = 0; i < n; ++i) {
        for (std::size_t j = 0; j < m; ++j) { Pi(i, j) = static_cast<double>(std::exp(u[i] + kernel(i, j) + v[j])); }
      }
      if (!Pi.IsFinite()) { throw std::runtime_error("opt::SinkHorn: non-finite coupling"); }
      PiT = Pi.Transpose();

      // Marginal convergence residuals
      p_residual = resnorm_p(Pi);
      q_residual = resnorm_q(PiT);

      printf(
          "iter = %4zu / %4zu : p_residual = %0.5E, log(p_residual) = "
          "%4.1f, q_residual = %0.5E, log(q_residual) = "
          "%4.1f \n",
          k + 1, iter, p_residual, std::log(p_residual), q_residual, std::log(q_residual));
    }
  }

  return p_residual;
}

}  // namespace opt
}  // namespace gra
