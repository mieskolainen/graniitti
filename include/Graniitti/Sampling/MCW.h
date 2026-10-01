// Monte Carlo Weight Objects [HEADER ONLY file]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MCW_H
#define MCW_H

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <complex>
#include <limits>
#include <random>
#include <stdexcept>
#include <valarray>
#include <vector>

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MStatistics.h"


namespace gra {
namespace kinematics {

// A simple Monte Carlo weight (container)
class MCW {
 public:
  MCW() = default;

  // Construct one single-observation Monte Carlo weight
  explicit MCW(double observation) { moments.Add(observation); }

  // Construct one Monte Carlo weight from raw moments
  MCW(double w, double w2, double n) {
    RestoreRawMoments(w, w2, n);
  }

  // Add two independent Monte Carlo moment states
  MCW operator+(const MCW &obj) const {
    MCW res = *this;
    res += obj;
    return res;
  }

  // Merge one independent Monte Carlo moment state
  MCW &operator+=(const MCW &obj) {
    if ((moments.Count() > 0 && !math::IsZero(obj.empty_sum)) ||
        (obj.moments.Count() > 0 && !math::IsZero(empty_sum))) {
      throw std::logic_error("MCW::operator+=: invalid empty weight state");
    }
    moments.Merge(obj.moments);
    if (moments.Count() == 0) {
      empty_sum += obj.empty_sum;
      empty_square_sum += obj.empty_square_sum;
    } else {
      empty_sum = 0.0;
      empty_square_sum = 0.0;
    }
    return *this;
  }

  // Scale all weights without changing their observation count
  // w_i <- scale w_i
  MCW operator*(double scale) const {
    if (!std::isfinite(scale)) {
      throw std::invalid_argument("MCW::operator*: invalid scale");
    }
    MCW result = *this;
    result.moments.Scale(static_cast<long double>(scale));
    result.empty_sum *= scale;
    result.empty_square_sum *= scale * scale;
    return result;
  }

  // Push one new observation into the stable moments
  void Push(double w) { moments.Add(w); }

  // MC estimate of the integral
  // I_hat = N^-1 sum_i w_i
  double Integral() const {
    if (!moments.IsFinite()) { return std::numeric_limits<double>::quiet_NaN(); }
    return moments.Count() > 0 ? static_cast<double>(moments.Mean()) : 0.0;
  }

  // Compute the unbiased variance of the Monte Carlo mean
  // Var(I_hat) = sum_i (w_i-I_hat)^2/[N(N-1)]
  double IntegralError2() const {
    if (!moments.IsFinite()) { return std::numeric_limits<double>::quiet_NaN(); }
    if (moments.Count() < 2) {
      return 0.0;
    }
    const long double count =
        static_cast<long double>(moments.Count());
    return static_cast<double>(moments.M2() / (count * (count - 1.0L)));
  }

  // Compute the standard error of the Monte Carlo mean
  double IntegralError() const { return std::sqrt(IntegralError2()); }

  // Compute the raw sum of observations
  // W = sum_i w_i
  double GetW() const {
    if (!moments.IsFinite()) { return std::numeric_limits<double>::quiet_NaN(); }
    if (moments.Count() == 0) {
      return empty_sum;
    }
    return static_cast<double>(moments.Mean() *
                               static_cast<long double>(moments.Count()));
  }

  // Compute the raw squared-observation sum
  // W2 = sum_i w_i^2
  double GetW2() const {
    if (!moments.IsFinite()) { return std::numeric_limits<double>::quiet_NaN(); }
    if (moments.Count() == 0) {
      return empty_square_sum;
    }
    const long double count =
        static_cast<long double>(moments.Count());
    return static_cast<double>(moments.M2() +
                               count * moments.Mean() * moments.Mean());
  }

  // Compute the observation count
  double GetN() const { return static_cast<double>(moments.Count()); }

 private:
  // Restore stable moments from an integer count and raw sums
  void RestoreRawMoments(double sum, double square_sum, double count_value) {
    if (!std::isfinite(sum) || !std::isfinite(square_sum) ||
        !std::isfinite(count_value) || count_value < 0.0 ||
        !math::IsExactEqual(std::floor(count_value), count_value) ||
        static_cast<long double>(count_value) >
            static_cast<long double>(
                std::numeric_limits<std::size_t>::max())) {
      throw std::invalid_argument("MCW: invalid raw moments");
    }
    const std::size_t count = static_cast<std::size_t>(count_value);
    if (count == 0) {
      empty_sum = sum;
      empty_square_sum = square_sum;
      return;
    }
    const long double mean =
        static_cast<long double>(sum) / static_cast<long double>(count);
    long double m2 = 0.0L;
    if (count > 1) {
      m2 = static_cast<long double>(square_sum) -
           static_cast<long double>(count) * mean * mean;
      const long double tolerance =
          16.0L * std::numeric_limits<double>::epsilon() *
          std::max(std::abs(static_cast<long double>(square_sum)),
                   static_cast<long double>(count) * mean * mean);
      if (m2 < -tolerance) {
        throw std::invalid_argument("MCW: inconsistent raw moments");
      }
      m2 = std::max(m2, 0.0L);
    }
    moments.Restore(count, mean, m2, true);
  }

  statistics::RunningMoments moments;
  double empty_sum = 0.0;
  double empty_square_sum = 0.0;
};


// A weighted sum of MCW container integral values (use with VEGAS, for example)
class MCWSUM {
 public:
  MCWSUM() = default;

  // Add one estimate using a finite positive linear weight
  void Add(const MCW &x, double weight) {
    if (math::IsZero(weight)) {
      return;
    }
    if (!(weight > 0.0) || !std::isfinite(weight)) {
      throw std::invalid_argument("MCWSUM::Add: invalid weight");
    }
    AddLogWeight(x, std::log(weight));
  }

  // Add one estimate using a logarithmic positive weight
  void AddLogWeight(const MCW &x, double log_weight) {
    if (std::isnan(log_weight) ||
        (std::isinf(log_weight) && !std::signbit(log_weight))) {
      throw std::invalid_argument("MCWSUM::AddLogWeight: invalid weight");
    }
    if (std::isinf(log_weight) && std::signbit(log_weight)) {
      return;
    }
    const double integral = x.Integral();
    const double error2 = x.IntegralError2();
    if (!std::isfinite(integral) || !std::isfinite(error2) || error2 < 0.0) {
      throw std::invalid_argument("MCWSUM::AddLogWeight: invalid estimate");
    }

    if (samples == 0) {
      max_log_weight = static_cast<long double>(log_weight);
    } else if (static_cast<long double>(log_weight) > max_log_weight) {
      const long double ratio =
          std::exp(max_log_weight - static_cast<long double>(log_weight));
      scaled_weight_sum.Scale(ratio);
      scaled_integral_sum.Scale(ratio);
      scaled_error2_sum.Scale(ratio * ratio);
      max_log_weight = static_cast<long double>(log_weight);
    }
    const long double normalized =
        std::exp(static_cast<long double>(log_weight) - max_log_weight);
    scaled_weight_sum.Add(normalized);
    scaled_integral_sum.Add(normalized * integral);
    scaled_error2_sum.Add(normalized * normalized * error2);
    ++samples;
  }

  // Compute the stable weighted integral mean
  // I = sum_i a_i I_i / sum_i a_i
  double Integral() const {
    const long double weight_sum = scaled_weight_sum.Value();
    if (!(weight_sum > 0.0L)) {
      return 0.0;
    }
    return static_cast<double>(scaled_integral_sum.Value() / weight_sum);
  }

  // Compute the stable variance of the weighted integral mean
  // Var(I) = sum_i a_i^2 Var(I_i)/(sum_i a_i)^2
  double IntegralError2() const {
    const long double weight_sum = scaled_weight_sum.Value();
    if (!(weight_sum > 0.0L)) {
      return 0.0;
    }
    return static_cast<double>(scaled_error2_sum.Value() /
                               (weight_sum * weight_sum));
  }

  // Compute the standard error of the weighted integral mean
  double IntegralError() const { return gra::math::msqrt(IntegralError2()); }

 private:
  statistics::CompensatedSum scaled_weight_sum;
  statistics::CompensatedSum scaled_integral_sum;
  statistics::CompensatedSum scaled_error2_sum;
  long double max_log_weight = 0.0L;
  std::size_t samples = 0;
};

}  // namespace kinematics
}  // namespace gra

#endif
