// Numerical integration rules and sampled integrals
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MINTEGRATION_H
#define MINTEGRATION_H

#include <cmath>
#include <cstddef>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Math/MMatrix.h"

namespace gra::math {

// Compute the composite Simpson or trapezoidal integration weight
double CompositeWeight(std::size_t index, std::size_t intervals);

// Integrate x f(x) exactly when f is represented by linear grid interpolation
double LinearRadialIntegral(const std::vector<double> &node, const std::vector<double> &value);

// Integrate closed uniform samples with the composite trapezoidal rule
// I = h[f_0/2 + sum_{i=1}^{N-1} f_i + f_N/2]
template <typename T>
inline T CSTrapzIntegral(const std::vector<T> &f, const double step) {
  if (f.size() < 2) { throw std::invalid_argument("gra::math::CSTrapzIntegral: at least two samples required"); }
  if (!std::isfinite(step)) { throw std::invalid_argument("gra::math::CSTrapzIntegral: step must be finite"); }
  const std::size_t intervals = f.size() - 1;
  T                 sum       = 0.0;
  for (std::size_t j = 1; j < intervals; ++j) { sum += f[j]; }
  return step * (f.front() + 2 * sum + f.back()) / 2.0;
}

// Integrate one uniform period without a duplicated boundary sample
// I = h sum_{i=0}^{N-1} f_i
template <typename T>
inline T PeriodicTrapzIntegral(const std::vector<T> &f, const double step) {
  if (f.empty()) { throw std::invalid_argument("gra::math::PeriodicTrapzIntegral: no samples"); }
  if (!std::isfinite(step)) { throw std::invalid_argument("gra::math::PeriodicTrapzIntegral: step must be finite"); }
  T integral = 0.0;
  for (const auto &value : f) { integral += value; }
  return integral * step;
}

// Compute midpoint periodic trapezoid nodes and scaled weights
std::pair<std::vector<double>, std::vector<double>> PeriodicTrapzRule(unsigned int n, double a, double b);

// Integrate closed uniform samples with the composite Simpson 1/3 rule
// I = h/3 [f_0 + f_N + 4 sum_odd f_i + 2 sum_even f_i]
template <typename T>
inline T CS13Integral(const std::vector<T> &f, const double step) {
  if (f.size() < 3) { throw std::invalid_argument("gra::math::CS13Integral: at least three samples required"); }
  if ((f.size() - 1) % 2 != 0) { throw std::invalid_argument("gra::math::CS13Integral: interval count must be even"); }
  if (!std::isfinite(step)) { throw std::invalid_argument("gra::math::CS13Integral: step must be finite"); }
  const std::size_t panels = (f.size() - 1) / 2;
  T                 even   = 0.0;
  T                 odd    = 0.0;
  for (std::size_t j = 1; j < panels; ++j) { even += f[2 * j]; }
  for (std::size_t j = 1; j <= panels; ++j) { odd += f[2 * j - 1]; }
  return step * (f.front() + 2.0 * even + 4.0 * odd + f.back()) / 3.0;
}

// Integrate closed uniform samples with the composite Simpson 3/8 rule
// I = 3h/8 [f_0 + f_N + 3 sum_{i mod 3 != 0} f_i + 2 sum_{i mod 3 = 0} f_i]
template <typename T>
inline T CS38Integral(const std::vector<T> &f, const double step) {
  if (f.size() < 4) { throw std::invalid_argument("gra::math::CS38Integral: at least four samples required"); }
  if ((f.size() - 1) % 3 != 0) {
    throw std::invalid_argument("gra::math::CS38Integral: interval count must be a multiple of 3");
  }
  if (!std::isfinite(step)) { throw std::invalid_argument("gra::math::CS38Integral: step must be finite"); }
  const std::size_t panels      = (f.size() - 1) / 3;
  T                 nonmultiple = 0.0;
  T                 multiple    = 0.0;
  for (std::size_t j = 1; j <= panels; ++j) { nonmultiple += f[3 * j - 2] + f[3 * j - 1]; }
  for (std::size_t j = 1; j < panels; ++j) { multiple += f[3 * j]; }
  return 3.0 * step * (f.front() + 3.0 * nonmultiple + 2.0 * multiple + f.back()) / 8.0;
}

// Integrate closed uniform samples with the composite Boole rule
// I_panel = 2h/45 (7f_0 + 32f_1 + 12f_2 + 32f_3 + 7f_4)
template <typename T>
inline T CSBooleIntegral(const std::vector<T> &f, const double step) {
  if (f.size() < 5) { throw std::invalid_argument("gra::math::CSBooleIntegral: at least five samples required"); }
  if ((f.size() - 1) % 4 != 0) {
    throw std::invalid_argument("gra::math::CSBooleIntegral: interval count must be a multiple of 4");
  }
  if (!std::isfinite(step)) { throw std::invalid_argument("gra::math::CSBooleIntegral: step must be finite"); }
  T integral = 0.0;
  for (std::size_t j = 0; j + 4 < f.size(); j += 4) {
    integral += 7.0 * f[j] + 32.0 * f[j + 1] + 12.0 * f[j + 2] + 32.0 * f[j + 3] + 7.0 * f[j + 4];
  }
  return integral * 2.0 * step / 45.0;
}

// Integrate a closed rectangular grid with tensor Simpson 1/3 weights
// I = h_x h_y/9 sum_{ij} w_i w_j f_{ij}, w = (1,4,2,...,4,1)
template <typename T>
inline T Simpson13Integral2D(const MMatrix<T> &f, const MMatrix<double> &weight, const double row_step,
                             const double col_step) {
  if (f.size_row() < 3 || f.size_col() < 3) {
    throw std::invalid_argument("gra::math::Simpson13Integral2D: each dimension needs three samples");
  }
  if (f.size_row() != weight.size_row() || f.size_col() != weight.size_col()) {
    throw std::invalid_argument("gra::math::Simpson13Integral2D: data and weight dimensions differ");
  }
  if ((f.size_row() - 1) % 2 != 0 || (f.size_col() - 1) % 2 != 0) {
    throw std::invalid_argument("gra::math::Simpson13Integral2D: interval counts must be even");
  }
  if (!std::isfinite(row_step) || !std::isfinite(col_step)) {
    throw std::invalid_argument("gra::math::Simpson13Integral2D: steps must be finite");
  }
  T integral = f.ElementwiseProductSum(weight);
  integral *= row_step * col_step / 9.0;
  return integral;
}

// Integrate a closed rectangular grid with tensor Simpson 3/8 weights
// I = 9h_x h_y/64 sum_{ij} w_i w_j f_{ij}, w = (1,3,3,2,...,3,3,1)
template <typename T>
inline T Simpson38Integral2D(const MMatrix<T> &f, const MMatrix<double> &weight, const double row_step,
                             const double col_step) {
  if (f.size_row() < 4 || f.size_col() < 4) {
    throw std::invalid_argument("gra::math::Simpson38Integral2D: each dimension needs four samples");
  }
  if (f.size_row() != weight.size_row() || f.size_col() != weight.size_col()) {
    throw std::invalid_argument("gra::math::Simpson38Integral2D: data and weight dimensions differ");
  }
  if ((f.size_row() - 1) % 3 != 0 || (f.size_col() - 1) % 3 != 0) {
    throw std::invalid_argument("gra::math::Simpson38Integral2D: interval counts must be multiples of 3");
  }
  if (!std::isfinite(row_step) || !std::isfinite(col_step)) {
    throw std::invalid_argument("gra::math::Simpson38Integral2D: steps must be finite");
  }
  T integral = f.ElementwiseProductSum(weight);
  integral *= 9.0 * row_step * col_step / 64.0;
  return integral;
}

// Integrate a closed rectangular grid with tensor Boole weights
// I = 4h_x h_y/2025 sum_{ij} w_i w_j f_{ij}, w = (7,32,12,32,14,...,7)
template <typename T>
inline T BooleIntegral2D(const MMatrix<T> &f, const MMatrix<double> &weight, const double row_step,
                         const double col_step) {
  if (f.size_row() < 5 || f.size_col() < 5) {
    throw std::invalid_argument("gra::math::BooleIntegral2D: each dimension needs five samples");
  }
  if (f.size_row() != weight.size_row() || f.size_col() != weight.size_col()) {
    throw std::invalid_argument("gra::math::BooleIntegral2D: data and weight dimensions differ");
  }
  if ((f.size_row() - 1) % 4 != 0 || (f.size_col() - 1) % 4 != 0) {
    throw std::invalid_argument("gra::math::BooleIntegral2D: interval counts must be multiples of 4");
  }
  if (!std::isfinite(row_step) || !std::isfinite(col_step)) {
    throw std::invalid_argument("gra::math::BooleIntegral2D: steps must be finite");
  }
  T integral = f.ElementwiseProductSum(weight);
  integral *= 4.0 * row_step * col_step / 2025.0;
  return integral;
}

// Compute Gauss Legendre nodes and scaled weights on a finite interval
std::pair<std::vector<double>, std::vector<double>> GaussLegendreRule(unsigned int n, double a, double b);

// Integrate one positive variable with logarithmic Gauss Legendre nodes
// int_a^b f(x) dx = int_log(a)^log(b) exp(u) f(exp(u)) du
template <typename Evaluator>
inline double LogGaussIntegral(const unsigned int n, const double minimum, const double maximum, Evaluator evaluator) {
  if (n == 0 || !std::isfinite(minimum) || !std::isfinite(maximum) || !(minimum > 0.0) || !(maximum > minimum)) {
    throw std::invalid_argument("gra::math::LogGaussIntegral: invalid rule");
  }
  const auto [node, weight] = GaussLegendreRule(n, std::log(minimum), std::log(maximum));
  double integral           = 0.0;
  for (std::size_t i = 0; i < node.size(); ++i) {
    const double value = std::exp(node[i]);
    integral += weight[i] * value * evaluator(value);
  }
  return integral;
}

// Integrate one positive variable using the logarithmic measure dlog(x)
// int_a^b f(x) dlog(x) = int_log(a)^log(b) f(exp(u)) du
template <typename Evaluator>
inline double LogMeasureGaussIntegral(const unsigned int n, const double minimum, const double maximum,
                                      Evaluator evaluator) {
  if (n == 0 || !std::isfinite(minimum) || !std::isfinite(maximum) || !(minimum > 0.0) || !(maximum > minimum)) {
    throw std::invalid_argument("gra::math::LogMeasureGaussIntegral: invalid rule");
  }
  const auto [node, weight] = GaussLegendreRule(n, std::log(minimum), std::log(maximum));
  double integral           = 0.0;
  for (std::size_t i = 0; i < node.size(); ++i) { integral += weight[i] * evaluator(std::exp(node[i])); }
  return integral;
}

// Compute the closed grid Simpson 1/3 weight vector
std::vector<double> Simpson13Weight(unsigned int intervals);

// Compute the closed grid Simpson 3/8 weight vector
std::vector<double> Simpson38Weight(unsigned int intervals);

// Compute the closed grid Boole weight vector
std::vector<double> BooleWeight(unsigned int intervals);

// Compute the tensor Simpson 1/3 weight matrix
MMatrix<double> Simpson13Weight2D(unsigned int rows, unsigned int cols);

// Compute the tensor Simpson 3/8 weight matrix
MMatrix<double> Simpson38Weight2D(unsigned int rows, unsigned int cols);

// Compute the tensor Boole weight matrix
MMatrix<double> BooleWeight2D(unsigned int rows, unsigned int cols);

}  // namespace gra::math

#endif
