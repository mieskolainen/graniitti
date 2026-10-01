// Stable statistical accumulation utilities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSTATISTICS_H
#define MSTATISTICS_H

#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <limits>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMatrix.h"

namespace gra {
namespace statistics {

// Store the first two moments and population variance of a complex sample
struct ComplexStat {
  std::complex<double> mean = 0.0;
  double second = 0.0;
  double variance = 0.0;
};

// Compute the population moments of one finite complex sample
ComplexStat ComplexMoments(const std::vector<std::complex<double>> &sample);

// Compute the unbiased variance from finite complex sample moments
double UnbiasedComplexVariance(double second, std::complex<double> mean,
                               std::size_t count);

// Compute the relative variance of positive log values raised to one power
double PoweredVariance(const std::vector<double> &log_node, double power);

// Compute self-normalized weighted-mean variance from scaled weight sums
double WeightedMeanVariance(long double scaled_weight_sum, long double scaled_weight2_sum,
                            long double scaled_centered_sum);

// Compute delete-group jackknife covariance around the replica mean
MMatrix<double> JackknifeCovariance(const std::vector<std::vector<double>> &replicas);

// Combine one finite integrand and inverse density at floating-point limits
// w = f exp(log(1/q))
inline double LogRangeImportanceWeight(double integrand,
                                       double log_inverse_density) {
  const long double log_magnitude =
      std::log(std::abs(static_cast<long double>(integrand))) +
      static_cast<long double>(log_inverse_density);
  const long double maximum = std::log(
      static_cast<long double>(std::numeric_limits<double>::max()));
  const long double minimum = std::log(
      static_cast<long double>(std::numeric_limits<double>::denorm_min()));
  if (log_magnitude > maximum) {
    return std::copysign(std::numeric_limits<double>::infinity(), integrand);
  }
  if (log_magnitude < minimum) {
    return std::copysign(0.0, integrand);
  }
  return std::copysign(
      static_cast<double>(std::exp(log_magnitude)), integrand);
}

// Combine one integrand with a logarithmic inverse proposal density
// w = f/q = f exp(log(1/q))
inline double ImportanceWeight(double integrand,
                               double log_inverse_density) {
  if (math::IsZero(integrand)) {
    return 0.0;
  }
  if (!std::isfinite(integrand) || std::isnan(log_inverse_density)) {
    return integrand * std::numeric_limits<double>::quiet_NaN();
  }
  if (std::isinf(log_inverse_density) &&
      std::signbit(log_inverse_density)) {
    return std::copysign(0.0, integrand);
  }
  if (std::isinf(log_inverse_density) &&
      !std::signbit(log_inverse_density)) {
    return std::copysign(std::numeric_limits<double>::infinity(), integrand);
  }
  if (math::IsZero(log_inverse_density)) {
    return integrand;
  }

  const double inverse_density = std::exp(log_inverse_density);
  if (std::isfinite(inverse_density) &&
      inverse_density >= std::numeric_limits<double>::min()) {
    return integrand * inverse_density;
  }
  return LogRangeImportanceWeight(integrand, log_inverse_density);
}

// Importance weight together with a directly representable proposal density
struct ImportanceSample {
  double weight = 0.0;
  double proposal_density = 0.0;

  // Compute whether the proposal density is safe for direct accumulation
  bool HasMaterializedProposalDensity() const {
    return std::isfinite(proposal_density) &&
           proposal_density >= std::numeric_limits<double>::min();
  }
};

// Evaluate one integration sample while materializing the proposal only once
// weight = integrand/proposal_density
inline ImportanceSample EvaluateImportanceSample(
    double integrand, double log_inverse_density) {
  if (math::IsZero(log_inverse_density) && std::isfinite(integrand)) {
    return {integrand, 1.0};
  }

  const double proposal_density = std::exp(-log_inverse_density);
  if (std::isfinite(integrand) && !std::isnan(log_inverse_density) &&
      std::isfinite(proposal_density) &&
      proposal_density >= std::numeric_limits<double>::min()) {
    return {integrand / proposal_density, proposal_density};
  }
  return {ImportanceWeight(integrand, log_inverse_density),
          proposal_density};
}

// Portable mantissa and exponent representation of one long double
struct ExtendedFloatParts {
  double mantissa = 0.0;
  int    exponent = 0;
};

// Split one finite long double without narrowing its exponent range
// value = mantissa 2^exponent with |mantissa| in [1/2,1)
inline ExtendedFloatParts SplitExtendedFloat(long double value) {
  if (!std::isfinite(value)) {
    throw std::invalid_argument(
        "SplitExtendedFloat: non-finite input");
  }
  int exponent = 0;
  const long double mantissa = std::frexp(value, &exponent);
  double encoded_mantissa = static_cast<double>(mantissa);
  if (math::IsExactEqual(std::abs(encoded_mantissa), 1.0)) {
    encoded_mantissa = std::copysign(0.5, encoded_mantissa);
    ++exponent;
  }
  return {encoded_mantissa, exponent};
}

// Restore one finite long double from its mantissa and exponent
// value = mantissa 2^exponent
inline long double JoinExtendedFloat(const ExtendedFloatParts &parts) {
  if (!std::isfinite(parts.mantissa) ||
      std::abs(parts.mantissa) >= 1.0 ||
      (!math::IsZero(parts.mantissa) &&
       std::abs(parts.mantissa) < 0.5)) {
    throw std::invalid_argument(
        "JoinExtendedFloat: invalid mantissa");
  }
  const long double value =
      std::ldexp(static_cast<long double>(parts.mantissa), parts.exponent);
  if (!std::isfinite(value)) {
    throw std::invalid_argument(
        "JoinExtendedFloat: exponent outside long-double range");
  }
  return value;
}

// Compensated long-double sum with support for exact common rescaling
class CompensatedSum {
 public:
  // Reset the accumulated sum
  void Reset() {
    sum_        = 0.0L;
    correction_ = 0.0L;
  }

  // Add one finite value using the Neumaier recurrence
  // sum' = sum + value with the lost low part accumulated in correction
  void Add(long double value) {
    const long double updated = sum_ + value;
    if (std::abs(sum_) >= std::abs(value)) {
      correction_ += (sum_ - updated) + value;
    } else {
      correction_ += (value - updated) + sum_;
    }
    sum_ = updated;
  }

  // Multiply the complete compensated state by one finite factor
  void Scale(long double factor) {
    sum_ *= factor;
    correction_ *= factor;
  }

  // Compute the compensated sum
  // value = sum + correction
  long double Value() const { return sum_ + correction_; }

 private:
  long double sum_        = 0.0L;
  long double correction_ = 0.0L;
};

// Welford moments for finite scalar observations
class RunningMoments {
 public:
  // Reset all accumulated moments
  void Reset() {
    mean_   = 0.0L;
    m2_     = 0.0L;
    count_  = 0;
    finite_ = true;
  }

  // Add one observation using the Welford recurrence
  // mean_n = mean_{n-1} + delta/n, M2_n = M2_{n-1} + delta(x_n-mean_n)
  void Add(long double value) {
    if (count_ == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("RunningMoments::Add: count overflow");
    }
    ++count_;
    if (!std::isfinite(value)) {
      finite_ = false;
      return;
    }
    const long double delta = value - mean_;
    mean_ += delta / static_cast<long double>(count_);
    m2_ += delta * (value - mean_);
  }

  // Merge an independent moment state using the Chan recurrence
  // M2 = M2_a + M2_b + (mean_b-mean_a)^2 n_a n_b/(n_a+n_b)
  void Merge(const RunningMoments &other) {
    if (other.count_ == 0) {
      return;
    }
    if (count_ == 0) {
      *this = other;
      return;
    }
    if (other.count_ >
        std::numeric_limits<std::size_t>::max() - count_) {
      throw std::overflow_error("RunningMoments::Merge: count overflow");
    }

    const std::size_t merged_count = count_ + other.count_;
    const long double left_count =
        static_cast<long double>(count_);
    const long double right_count =
        static_cast<long double>(other.count_);
    const long double total_count =
        static_cast<long double>(merged_count);
    const long double delta = other.mean_ - mean_;
    mean_ += delta * right_count / total_count;
    m2_ += other.m2_ +
           delta * delta * left_count * right_count / total_count;
    count_ = merged_count;
    finite_ = finite_ && other.finite_;
  }

  // Scale every accumulated observation by one finite factor
  // mean' = c mean, M2' = c^2 M2
  void Scale(long double factor) {
    if (!std::isfinite(factor)) {
      throw std::invalid_argument("RunningMoments::Scale: invalid factor");
    }
    mean_ *= factor;
    m2_ *= factor * factor;
    finite_ = finite_ && std::isfinite(mean_) && std::isfinite(m2_);
  }

  // Restore one serialized moment state after validating its invariants
  void Restore(std::size_t count, long double mean, long double m2,
               bool finite) {
    if (!std::isfinite(mean) || !std::isfinite(m2) || m2 < 0.0L ||
        (count == 0 && (!math::IsZero(mean) || !math::IsZero(m2)))) {
      throw std::invalid_argument(
          "RunningMoments::Restore: invalid moment state");
    }
    count_  = count;
    mean_   = mean;
    m2_     = m2;
    finite_ = finite;
  }

  // Compute the number of observations
  std::size_t Count() const { return count_; }

  // Compute the running mean
  long double Mean() const { return mean_; }

  // Compute the centered sum of squares
  long double M2() const { return m2_; }

  // Check whether every observation was finite
  bool IsFinite() const { return finite_; }

 private:
  long double mean_   = 0.0L;
  long double m2_     = 0.0L;
  std::size_t count_  = 0;
  bool        finite_ = true;
};

// Scale-normalized moments for nonnegative samples supplied directly or in
// log space
class ScaledPositiveMoments {
 public:
  // Reset all accumulated moments and positive extrema
  void Reset() {
    scaled_moments_.Reset();
    positive_count_ = 0;
    log_scale_ = 0.0L;
    min_log_value_ = 0.0L;
    materialized_scale_ = 0.0L;
    materialized_minimum_ = 0.0L;
  }

  // Add one nonnegative observation without throwing on invalid sample values
  void Add(double value) {
    if (!(value >= 0.0) || !std::isfinite(value)) {
      AddInvalid();
      return;
    }
    if (math::IsZero(value)) {
      scaled_moments_.Add(0.0L);
      return;
    }
    AddFiniteValue(static_cast<long double>(value));
  }

  // Add one positive observation represented by its natural logarithm
  void AddLogPositive(long double log_value) {
    if (!std::isfinite(log_value)) {
      AddInvalid();
      return;
    }
    AddFiniteLog(log_value);
  }

  // Restore one serialized scaled-moment state after validating its invariants
  void Restore(std::size_t count, std::size_t positive_count,
               long double scaled_mean, long double scaled_m2,
               long double log_scale, long double min_log_value,
               bool finite) {
    const long double count_value = static_cast<long double>(count);
    const long double positive_value =
        static_cast<long double>(positive_count);
    const long double tolerance =
        256.0L * std::numeric_limits<double>::epsilon();
    const bool invalid_normalized_mean =
        finite && positive_count > 0 &&
        (scaled_mean + tolerance < 1.0L / count_value ||
         scaled_mean > positive_value / count_value + tolerance);
    const long double maximum_m2 =
        finite && positive_count > 0
            ? count_value * scaled_mean *
                  std::max(1.0L - scaled_mean, 0.0L)
            : 0.0L;
    const bool invalid_normalized_m2 =
        finite && scaled_m2 > maximum_m2 + tolerance * count_value;
    if (positive_count > count ||
        (positive_count == 0 &&
         (!math::IsZero(scaled_mean) || !math::IsZero(scaled_m2) ||
          !math::IsZero(log_scale) || !math::IsZero(min_log_value))) ||
        (positive_count > 0 &&
         (!(scaled_mean > 0.0L) || !std::isfinite(log_scale) ||
          !std::isfinite(min_log_value) || min_log_value > log_scale)) ||
        invalid_normalized_mean || invalid_normalized_m2) {
      throw std::invalid_argument(
          "ScaledPositiveMoments::Restore: invalid moment state");
    }
    scaled_moments_.Restore(count, scaled_mean, scaled_m2, finite);
    positive_count_ = positive_count;
    log_scale_ = log_scale;
    min_log_value_ = min_log_value;
    RestoreMaterializedExtrema();
  }

  // Compute the total number of observations including zeros
  std::size_t Count() const { return scaled_moments_.Count(); }

  // Compute the number of strictly positive observations
  std::size_t PositiveCount() const { return positive_count_; }

  // Compute the scale-normalized running mean
  long double ScaledMean() const { return scaled_moments_.Mean(); }

  // Compute the scale-normalized centered sum of squares
  long double ScaledM2() const { return scaled_moments_.M2(); }

  // Compute the logarithm of the largest positive observation
  long double LogScale() const { return log_scale_; }

  // Compute the logarithm of the smallest positive observation
  long double MinLogValue() const { return min_log_value_; }

  // Compute the physical sample mean
  // mean = exp(log_scale) scaled_mean
  long double Mean() const {
    if (!IsFinite()) {
      return std::numeric_limits<long double>::quiet_NaN();
    }
    if (positive_count_ == 0) {
      return 0.0L;
    }
    if (materialized_scale_ > 0.0L) {
      return materialized_scale_ * ScaledMean();
    }
    return std::exp(log_scale_ + std::log(ScaledMean()));
  }

  // Compute the sample standard deviation divided by the absolute mean
  // relative standard deviation = sqrt[M2/(N-1)]/|mean|
  long double RelativeStandardDeviation() const {
    if (!IsFinite()) {
      return std::numeric_limits<long double>::quiet_NaN();
    }
    if (Count() < 2 || positive_count_ == 0 || math::IsZero(ScaledM2())) {
      return 0.0L;
    }
    const long double variance =
        std::max(ScaledM2(), 0.0L) /
        static_cast<long double>(Count() - 1);
    return std::sqrt(variance) / std::abs(ScaledMean());
  }

  // Compute the effective sample size divided by the observation count
  // N_eff/N = N mean^2/[M2 + N mean^2]
  double EffectiveSampleFraction() const {
    if (!IsFinite() || Count() == 0 || positive_count_ == 0) {
      return 0.0;
    }
    const long double count = static_cast<long double>(Count());
    const long double squared_mean = ScaledMean() * ScaledMean();
    const long double denominator = ScaledM2() + count * squared_mean;
    if (!(denominator > 0.0L)) {
      return 0.0;
    }
    const long double fraction = count * squared_mean / denominator;
    return static_cast<double>(std::clamp(fraction, 0.0L, 1.0L));
  }

  // Compute the smallest strictly positive observation
  long double MinimumPositive() const {
    if (!IsFinite()) {
      return std::numeric_limits<long double>::quiet_NaN();
    }
    if (positive_count_ == 0) {
      return 0.0L;
    }
    return materialized_minimum_ > 0.0L ? materialized_minimum_
                                        : std::exp(min_log_value_);
  }

  // Compute the largest strictly positive observation
  long double Maximum() const {
    if (!IsFinite()) {
      return std::numeric_limits<long double>::quiet_NaN();
    }
    if (positive_count_ == 0) {
      return 0.0L;
    }
    return materialized_scale_ > 0.0L ? materialized_scale_
                                      : std::exp(log_scale_);
  }

  // Check whether every accumulated observation was finite and nonnegative
  bool IsFinite() const { return scaled_moments_.IsFinite(); }

 private:
  // Mark one invalid observation while preserving nonthrowing sampling behavior
  void AddInvalid() {
    scaled_moments_.Add(std::numeric_limits<long double>::quiet_NaN());
  }

  // Add one finite value without transcendental normalization
  void AddFiniteValue(long double value) {
    if (positive_count_ == 0) {
      materialized_scale_ = value;
      materialized_minimum_ = value;
      log_scale_ = std::log(value);
      min_log_value_ = log_scale_;
      AddNormalizedPositive(1.0L);
      return;
    }
    if (!(materialized_scale_ > 0.0L)) {
      AddFiniteLog(std::log(value));
      return;
    }
    if (value > materialized_scale_) {
      scaled_moments_.Scale(materialized_scale_ / value);
      materialized_scale_ = value;
      log_scale_ = std::log(value);
    }
    if (value < materialized_minimum_) {
      materialized_minimum_ = value;
      min_log_value_ = std::log(value);
    }
    AddNormalizedPositive(value / materialized_scale_);
  }

  // Add one finite positive logarithmic observation using a common scale
  void AddFiniteLog(long double log_value) {
    materialized_scale_ = 0.0L;
    materialized_minimum_ = 0.0L;
    long double normalized = 1.0L;
    if (positive_count_ == 0) {
      log_scale_ = log_value;
      min_log_value_ = log_value;
    } else {
      if (log_value > log_scale_) {
        const long double ratio = std::exp(log_scale_ - log_value);
        scaled_moments_.Scale(ratio);
        log_scale_ = log_value;
      } else {
        normalized = std::exp(log_value - log_scale_);
      }
      min_log_value_ = std::min(min_log_value_, log_value);
    }
    AddNormalizedPositive(normalized);
  }

  // Add one normalized positive observation and update its count
  void AddNormalizedPositive(long double normalized) {
    scaled_moments_.Add(normalized);
    ++positive_count_;
  }

  // Restore extrema caches when both values fit in long double
  void RestoreMaterializedExtrema() {
    materialized_scale_ = 0.0L;
    materialized_minimum_ = 0.0L;
    if (positive_count_ == 0) {
      return;
    }
    const long double scale = std::exp(log_scale_);
    const long double minimum = std::exp(min_log_value_);
    if (scale > 0.0L && minimum > 0.0L && std::isfinite(scale) &&
        std::isfinite(minimum)) {
      materialized_scale_ = scale;
      materialized_minimum_ = minimum;
    }
  }

  RunningMoments scaled_moments_;
  std::size_t    positive_count_ = 0;
  long double   log_scale_ = 0.0L;
  long double   min_log_value_ = 0.0L;
  long double   materialized_scale_ = 0.0L;
  long double   materialized_minimum_ = 0.0L;
};

// Scale-normalized signed-weight sums and squared-weight sums
class ScaledWeightSums {
 public:
  // Reset all accumulated weight sums
  void Reset() {
    scaled_sum_.Reset();
    scaled_square_sum_.Reset();
    scaled_absolute_sum_.Reset();
    scale_  = 0.0L;
    count_  = 0;
    finite_ = true;
  }

  // Add one signed weight without squaring its physical scale
  // sum w_i = scale sum(w_i/scale), sum w_i^2 = scale^2 sum(w_i/scale)^2
  void Add(double weight) {
    if (count_ == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error("ScaledWeightSums::Add: count overflow");
    }
    ++count_;
    if (!std::isfinite(weight)) {
      finite_ = false;
      return;
    }

    const long double value     = static_cast<long double>(weight);
    const long double magnitude = std::abs(value);
    if (math::IsZero(magnitude)) {
      return;
    }
    if (magnitude > scale_) {
      if (scale_ > 0.0L) {
        const long double ratio = scale_ / magnitude;
        scaled_sum_.Scale(ratio);
        scaled_absolute_sum_.Scale(ratio);
        scaled_square_sum_.Scale(ratio * ratio);
      }
      scale_ = magnitude;
    }

    const long double normalized = value / scale_;
    scaled_sum_.Add(normalized);
    scaled_absolute_sum_.Add(std::abs(normalized));
    scaled_square_sum_.Add(normalized * normalized);
  }

  // Restore one serialized scale-normalized state
  void Restore(std::size_t count, long double scale, long double scaled_sum,
               long double scaled_square_sum,
               long double scaled_absolute_sum, bool finite) {
    if (!std::isfinite(scale) || !std::isfinite(scaled_sum) ||
        !std::isfinite(scaled_square_sum) ||
        !std::isfinite(scaled_absolute_sum) || scale < 0.0L ||
        scaled_square_sum < 0.0L || scaled_absolute_sum < 0.0L ||
        (count == 0 &&
         (!math::IsZero(scale) || !math::IsZero(scaled_sum) ||
          !math::IsZero(scaled_square_sum) ||
          !math::IsZero(scaled_absolute_sum))) ||
        (math::IsZero(scale) &&
         (!math::IsZero(scaled_sum) || !math::IsZero(scaled_square_sum) ||
          !math::IsZero(scaled_absolute_sum)))) {
      throw std::invalid_argument(
          "ScaledWeightSums::Restore: invalid weight state");
    }
    Reset();
    count_  = count;
    scale_  = scale;
    finite_ = finite;
    scaled_sum_.Add(scaled_sum);
    scaled_square_sum_.Add(scaled_square_sum);
    scaled_absolute_sum_.Add(scaled_absolute_sum);
  }

  // Compute the number of accumulated weights
  std::size_t Count() const { return count_; }

  // Compute the common absolute weight scale
  long double Scale() const { return scale_; }

  // Compute the signed sum divided by the common scale
  long double ScaledSum() const { return scaled_sum_.Value(); }

  // Compute the squared-weight sum divided by the squared common scale
  long double ScaledSquareSum() const {
    return scaled_square_sum_.Value();
  }

  // Compute the absolute-weight sum divided by the common scale
  long double ScaledAbsoluteSum() const {
    return scaled_absolute_sum_.Value();
  }

  // Compute the physical signed-weight sum in long-double range
  long double Sum() const { return scale_ * ScaledSum(); }

  // Compute the physical arithmetic mean
  long double Mean() const {
    if (!finite_) {
      return std::numeric_limits<long double>::quiet_NaN();
    }
    return count_ > 0
               ? scale_ * (ScaledSum() / static_cast<long double>(count_))
               : 0.0L;
  }

  // Compute the physical squared-weight sum in long-double range
  long double SquareSum() const {
    return scale_ * scale_ * ScaledSquareSum();
  }

  // Compute the scale-invariant effective sample size
  // N_eff = (sum_i w_i)^2/sum_i w_i^2
  long double EffectiveSampleSize() const {
    const long double square_sum = ScaledSquareSum();
    if (!(square_sum > 0.0L) || !finite_) {
      return 0.0L;
    }
    const long double sum = ScaledSum();
    return std::clamp(sum * sum / square_sum, 0.0L,
                      static_cast<long double>(count_));
  }

  // Compute the effective sample size divided by the observation count
  double EffectiveSampleFraction() const {
    if (count_ == 0) {
      return 0.0;
    }
    const long double fraction =
        EffectiveSampleSize() / static_cast<long double>(count_);
    return static_cast<double>(std::clamp(fraction, 0.0L, 1.0L));
  }

  // Compute whether the signed sum exceeds a relative cancellation threshold
  bool HasSignificantSignedSum(long double relative_tolerance) const {
    if (!(relative_tolerance >= 0.0L) ||
        !std::isfinite(relative_tolerance)) {
      throw std::invalid_argument(
          "ScaledWeightSums::HasSignificantSignedSum: invalid tolerance");
    }
    const long double absolute_sum = ScaledAbsoluteSum();
    return finite_ && absolute_sum > 0.0L &&
           std::abs(ScaledSum()) > relative_tolerance * absolute_sum;
  }

  // Check whether every accumulated weight was finite
  bool IsFinite() const { return finite_; }

 private:
  CompensatedSum scaled_sum_;
  CompensatedSum scaled_square_sum_;
  CompensatedSum scaled_absolute_sum_;
  long double    scale_  = 0.0L;
  std::size_t    count_  = 0;
  bool           finite_ = true;
};

// Accumulate weighted vector sums and their Poisson covariance with compensation
class WeightedVectorSums {
 public:
  // Allocate one fixed-dimensional moment accumulator
  explicit WeightedVectorSums(std::size_t dimension)
      : sum_(dimension), covariance_(dimension, dimension) {}

  // Add w b and w^2 b b^T without narrowing intermediate products
  template <typename Container>
  void Add(const Container &basis, double weight) {
    if (static_cast<std::size_t>(basis.size()) != sum_.size() ||
        !std::isfinite(weight) || !gra::AllFinite(basis)) {
      throw std::invalid_argument("WeightedVectorSums::Add: invalid dimensions or values");
    }
    weights_.Add(weight);
    for (std::size_t i = 0; i < sum_.size(); ++i) {
      const long double value = static_cast<long double>(weight) * basis[i];
      sum_[i].Add(value);
      for (std::size_t j = 0; j <= i; ++j) {
        covariance_[i][j].Add(value * (static_cast<long double>(weight) * basis[j]));
      }
    }
  }

  // Compute the weighted vector sum in double precision
  std::vector<double> Sum() const {
    std::vector<double> result(sum_.size());
    for (std::size_t i = 0; i < result.size(); ++i) { result[i] = static_cast<double>(sum_[i].Value()); }
    if (!gra::AllFinite(result)) { throw std::overflow_error("WeightedVectorSums::Sum: overflow"); }
    return result;
  }

  // Compute the symmetric Poisson covariance in double precision
  MMatrix<double> Covariance() const {
    MMatrix<double> result(sum_.size(), sum_.size());
    for (std::size_t i = 0; i < sum_.size(); ++i) {
      for (std::size_t j = 0; j <= i; ++j) {
        result[i][j] = result[j][i] = static_cast<double>(covariance_[i][j].Value());
      }
    }
    if (!result.IsFinite()) { throw std::overflow_error("WeightedVectorSums::Covariance: overflow"); }
    return result;
  }

  // Compute the event weight summary
  const ScaledWeightSums &Weights() const { return weights_; }

 private:
  std::vector<CompensatedSum> sum_;
  MMatrix<CompensatedSum> covariance_;
  ScaledWeightSums weights_;
};

// Stable inverse-variance moments for independent scalar estimates
class WeightedEstimateMoments {
 public:
  // Reset all accumulated estimates
  void Reset() {
    mean_           = 0.0L;
    scaled_m2_      = 0.0L;
    scaled_weight_  = 0.0L;
    max_log_weight_ = 0.0L;
    samples_        = 0;
  }

  // Add one finite estimate and its positive variance
  // mean = sum_i x_i/sigma_i^2 divided by sum_i 1/sigma_i^2
  void Add(long double value, long double variance) {
    if (!std::isfinite(value) || !std::isfinite(variance) ||
        !(variance > 0.0L)) {
      throw std::invalid_argument(
          "WeightedEstimateMoments::Add: invalid estimate");
    }
    if (samples_ == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error(
          "WeightedEstimateMoments::Add: count overflow");
    }

    const long double log_weight = -std::log(variance);
    if (samples_ == 0) {
      mean_           = value;
      scaled_weight_  = 1.0L;
      scaled_m2_      = 0.0L;
      max_log_weight_ = log_weight;
      samples_        = 1;
      return;
    }

    long double sample_weight = 1.0L;
    if (log_weight > max_log_weight_) {
      const long double factor = std::exp(max_log_weight_ - log_weight);
      scaled_weight_ *= factor;
      scaled_m2_ *= factor;
      max_log_weight_ = log_weight;
    } else {
      sample_weight = std::exp(log_weight - max_log_weight_);
    }

    const long double merged_weight = scaled_weight_ + sample_weight;
    const long double delta         = value - mean_;
    mean_ += delta * sample_weight / merged_weight;
    scaled_m2_ +=
        delta * delta * scaled_weight_ * sample_weight / merged_weight;
    scaled_weight_ = merged_weight;
    ++samples_;
  }

  // Restore one serialized inverse-variance moment state
  void Restore(std::size_t samples, long double mean, long double scaled_m2,
               long double scaled_weight, long double max_log_weight) {
    if (!std::isfinite(mean) || !std::isfinite(scaled_m2) ||
        !std::isfinite(scaled_weight) || !std::isfinite(max_log_weight) ||
        scaled_m2 < 0.0L || scaled_weight < 0.0L ||
        (samples == 0 &&
         (!math::IsZero(mean) || !math::IsZero(scaled_m2) ||
          !math::IsZero(scaled_weight) || !math::IsZero(max_log_weight))) ||
        (samples > 0 && !(scaled_weight > 0.0L))) {
      throw std::invalid_argument(
          "WeightedEstimateMoments::Restore: invalid estimate state");
    }
    samples_        = samples;
    mean_           = mean;
    scaled_m2_      = scaled_m2;
    scaled_weight_  = scaled_weight;
    max_log_weight_ = max_log_weight;
  }

  // Compute the inverse-variance weighted mean
  long double Mean() const { return mean_; }

  // Compute the number of independent estimates
  std::size_t Samples() const { return samples_; }

  // Compute the scaled centered sum used for serialization
  long double ScaledM2() const { return scaled_m2_; }

  // Compute the scaled inverse-variance sum used for serialization
  long double ScaledWeight() const { return scaled_weight_; }

  // Compute the logarithm of the common inverse-variance scale
  long double MaxLogWeight() const { return max_log_weight_; }

  // Compute the variance of the inverse-variance weighted mean
  // Var(mean) = 1/sum_i (1/sigma_i^2)
  long double VarianceOfMean() const {
    if (samples_ == 0 || !(scaled_weight_ > 0.0L)) {
      return std::numeric_limits<long double>::infinity();
    }
    return std::exp(-max_log_weight_) / scaled_weight_;
  }

  // Compute reduced chi2 without materializing the inverse-variance scale
  // chi2/ndf = sum_i (x_i-mean)^2/sigma_i^2 /(N-1)
  double ReducedChi2() const {
    if (samples_ < 2 || !(scaled_m2_ > 0.0L)) {
      return 0.0;
    }
    const long double log_chi2 =
        std::log(scaled_m2_) + max_log_weight_ -
        std::log(static_cast<long double>(samples_ - 1));
    const long double maximum =
        std::log(static_cast<long double>(
            std::numeric_limits<double>::max()));
    const long double minimum =
        std::log(static_cast<long double>(
            std::numeric_limits<double>::denorm_min()));
    if (log_chi2 > maximum) {
      return std::numeric_limits<double>::infinity();
    }
    if (log_chi2 < minimum) {
      return 0.0;
    }
    return static_cast<double>(std::exp(log_chi2));
  }

 private:
  long double mean_           = 0.0L;
  long double scaled_m2_      = 0.0L;
  long double scaled_weight_  = 0.0L;
  long double max_log_weight_ = 0.0L;
  std::size_t samples_        = 0;
};

// Clip one vector by its Euclidean norm without materializing that norm
// x' = x min[1, L/sqrt(sum_i x_i^2)]
inline void ClipEuclideanNorm(std::vector<double> &values,
                              double maximum_norm) {
  if (!(maximum_norm > 0.0) || !std::isfinite(maximum_norm)) {
    throw std::invalid_argument(
        "ClipEuclideanNorm: invalid maximum norm");
  }

  double scale       = 0.0;
  double square_sum  = 1.0;
  bool   has_nonzero = false;
  for (double value : values) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument(
          "ClipEuclideanNorm: non-finite vector component");
    }
    const double magnitude = std::abs(value);
    if (math::IsZero(magnitude)) {
      continue;
    }
    has_nonzero = true;
    if (magnitude > scale) {
      const double ratio = scale / magnitude;
      square_sum = 1.0 + square_sum * ratio * ratio;
      scale = magnitude;
    } else {
      const double ratio = magnitude / scale;
      square_sum += ratio * ratio;
    }
  }
  if (!has_nonzero) {
    return;
  }

  const double normalized_limit = maximum_norm / std::sqrt(square_sum);
  if (scale <= normalized_limit) {
    return;
  }
  for (double &value : values) {
    value = (value / scale) * normalized_limit;
  }
}

}  // namespace statistics
}  // namespace gra

#endif
