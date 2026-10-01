// Finite sample statistics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MStatistics.h"

#include <algorithm>
#include <cmath>
#include <complex>
#include <stdexcept>

namespace gra::statistics {

// Compute the population moments of one finite complex sample
// mean = N^{-1} sum_i z_i, variance = Var(Re z) + Var(Im z)
ComplexStat ComplexMoments(const std::vector<std::complex<double>> &sample) {
  if (sample.empty()) { throw std::invalid_argument("statistics::ComplexMoments: empty sample"); }
  ComplexStat stat;
  RunningMoments real;
  RunningMoments imaginary;
  CompensatedSum second;
  for (const auto &value : sample) {
    if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
      throw std::invalid_argument("statistics::ComplexMoments: non-finite sample");
    }
    real.Add(value.real());
    imaginary.Add(value.imag());
    second.Add(std::norm(std::complex<long double>(value.real(), value.imag())));
  }
  const long double count = static_cast<long double>(sample.size());
  stat.mean               = {static_cast<double>(real.Mean()), static_cast<double>(imaginary.Mean())};
  stat.second             = static_cast<double>(second.Value() / count);
  stat.variance           = static_cast<double>(std::max(0.0L, (real.M2() + imaginary.M2()) / count));
  return stat;
}

// Compute the unbiased variance from finite complex sample moments
// s^2 = N/(N-1) [E(|z|^2) - |E(z)|^2]
double UnbiasedComplexVariance(const double second, const std::complex<double> mean, const std::size_t count) {
  if (!std::isfinite(second) || !std::isfinite(mean.real()) || !std::isfinite(mean.imag()) || second < 0.0 ||
      count == 0) {
    throw std::invalid_argument("statistics::UnbiasedComplexVariance: invalid sampled moment");
  }
  if (count == 1) { return 0.0; }
  const double population = std::max(0.0, second - std::norm(mean));
  return population * static_cast<double>(count) / static_cast<double>(count - 1);
}

// Compute the relative variance of positive log values raised to one power
// relative variance = E[exp(2p log x)]/E[exp(p log x)]^2 - 1
double PoweredVariance(const std::vector<double> &log_node, const double power) {
  if (log_node.empty() || !std::isfinite(power)) {
    throw std::invalid_argument("statistics::PoweredVariance: invalid input");
  }
  if (!std::all_of(log_node.begin(), log_node.end(), [](const double value) { return std::isfinite(value); })) {
    throw std::invalid_argument("statistics::PoweredVariance: non-finite input");
  }
  if (math::IsZero(power)) { return 0.0; }
  const auto [minimum, maximum] = std::minmax_element(log_node.begin(), log_node.end());
  const long double shift       = power < 0.0 ? *minimum : *maximum;
  CompensatedSum    first;
  RunningMoments    centered;
  for (const double value : log_node) {
    const long double exponent = static_cast<long double>(power) * (value - shift);
    first.Add(std::exp(exponent));
    centered.Add(std::expm1(exponent));
  }
  const long double count = static_cast<long double>(log_node.size());
  const long double mean  = first.Value() / count;
  return static_cast<double>(std::max(0.0L, centered.M2() / count) / (mean * mean));
}

// Compute self-normalized weighted-mean variance from scaled weight sums
// Var(mean) = [sum_i w_i^2(y_i-mean)^2/(sum_i w_i)^2] N_eff/(N_eff-1)
double WeightedMeanVariance(const long double scaled_weight_sum, const long double scaled_weight2_sum,
                            const long double scaled_centered_sum) {
  if (std::fpclassify(scaled_weight_sum) == FP_ZERO || !(scaled_weight2_sum > 0.0L) || !(scaled_centered_sum >= 0.0L)) {
    return 0.0;
  }
  const long double effective_entries = scaled_weight_sum * scaled_weight_sum / scaled_weight2_sum;
  const long double correction = effective_entries > 1.0L ? effective_entries / (effective_entries - 1.0L) : 1.0L;
  const long double variance   = scaled_centered_sum / (scaled_weight_sum * scaled_weight_sum) * correction;
  return static_cast<double>(std::max(variance, 0.0L));
}

// Compute delete-group jackknife covariance around the replica mean
// Cov = (R-1)/R sum_r (theta_r - mean)(theta_r - mean)^T
MMatrix<double> JackknifeCovariance(const std::vector<std::vector<double>> &replicas) {
  if (replicas.empty()) { return {}; }
  const std::size_t dimension = replicas.front().size();
  for (const auto &replica : replicas) {
    if (replica.size() != dimension ||
        !std::all_of(replica.cbegin(), replica.cend(), [](const double value) { return std::isfinite(value); })) {
      throw std::invalid_argument("statistics::JackknifeCovariance: invalid replica");
    }
  }

  std::vector<double> mean(dimension, 0.0);
  for (const auto &replica : replicas) { gra::AddScaled(mean, replica, 1.0); }
  gra::Scale(mean, 1.0 / static_cast<double>(replicas.size()));

  MMatrix<double> covariance(dimension, dimension, 0.0);
  if (replicas.size() == 1) { return covariance; }
  const double factor = static_cast<double>(replicas.size() - 1) / static_cast<double>(replicas.size());
  for (const auto &replica : replicas) {
    const std::vector<double> residual = gra::Subtract(replica, mean);
    covariance.AddOuterProduct(residual, residual, factor);
  }
  return covariance;
}

}  // namespace gra::statistics
