// Generic numerical construction algorithms
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MALGORITHMS_H
#define MALGORITHMS_H

// C++
#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstddef>
#include <exception>
#include <iostream>
#include <limits>
#include <mutex>
#include <stdexcept>
#include <stop_token>
#include <string_view>
#include <thread>
#include <type_traits>
#include <utility>
#include <vector>

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
namespace math {

// Fill a closed arithmetic sequence without signed or unsigned overflow
// x_i = start + trunc[i(stop-start)/(N-1)] for integral values
template <typename Container, typename T>
inline void FillLinspace(Container &value, const T start, const T stop) {
  if (value.size() == 0) { throw std::invalid_argument("linspace: argument size == 0"); }
  if (value.size() == 1) {
    value[0] = start;
    return;
  }
  const std::size_t intervals = value.size() - 1;
  value[0]                    = start;
  value[intervals]            = stop;
  if constexpr (std::is_integral_v<T> && !std::is_same_v<T, bool>) {
    using Unsigned            = std::make_unsigned_t<T>;
    using Wide                = std::common_type_t<Unsigned, std::size_t>;
    const bool     descending = stop < start;
    const Unsigned distance   = descending ? Unsigned(start) - Unsigned(stop) : Unsigned(stop) - Unsigned(start);
    const Wide     step       = static_cast<Wide>(distance) / intervals;
    const Wide     remainder  = static_cast<Wide>(distance) % intervals;
    Wide           carry      = 0;
    Unsigned       current    = static_cast<Unsigned>(start);
    for (std::size_t i = 1; i < intervals; ++i) {
      const bool extra         = carry >= intervals - remainder;
      carry                    = extra ? carry - (intervals - remainder) : carry + remainder;
      const Unsigned increment = static_cast<Unsigned>(step + static_cast<Wide>(extra));
      current                  = descending ? current - increment : current + increment;
      value[i]                 = static_cast<T>(current);
    }
  } else {
    for (std::size_t i = 1; i < intervals; ++i) {
      if constexpr (std::is_floating_point_v<T>) {
        const T fraction = static_cast<T>(i) / static_cast<T>(intervals);
        value[i] = std::lerp(start, stop, fraction);
      } else {
        const double fraction = static_cast<double>(i) / static_cast<double>(intervals);
        value[i] = start + (stop - start) * fraction;
      }
    }
  }
}

// MATLAB style, H is std::vector, std::valarray or similar
// x_i = start + i(stop - start)/(N - 1), i = 0,...,N - 1
// use like:
// std::valarray<double> a = linspace<std::valarray>(0.0, 10.0, 16);
//
template <typename T>
inline std::vector<T> linspace(T start, T stop, std::size_t size) {
  if (size == 0) { throw std::invalid_argument("linspace: argument size == 0"); }
  if (size == 1) { return {start}; }

  // Resize separately to avoid GCC 11's spurious free-nonheap-object warning
  std::vector<T> v;
  v.resize(size);
  FillLinspace(v, start, stop);
  return v;
}

// Construct a closed arithmetic sequence in one selected container
// x_i = start + i(stop - start)/(N - 1), i = 0,...,N - 1
template <template <typename, typename...> class container_type, class value_type>
inline container_type<value_type> linspace(value_type start, value_type stop, std::size_t size) {
  if (size == 0) { throw std::invalid_argument("linspace: argument size == 0"); }
  if (size == 1) { return {start}; }

  container_type<value_type> v(size);
  FillLinspace(v, start, stop);
  return v;
}

// Construct a half-open arithmetic progression without integral promotion
// x_i = start + i step for all i with x_i before stop in the step direction
template <template <typename T> class container_type, class value_type>
inline container_type<value_type> arange(value_type start, value_type step, value_type stop) {
  static_assert(std::is_arithmetic_v<value_type>, "arange requires an arithmetic value type");
  if constexpr (std::is_floating_point_v<value_type>) {
    if (!std::isfinite(start) || !std::isfinite(step) || !std::isfinite(stop)) {
      throw std::invalid_argument("arange: arguments must be finite");
    }
  }
  if (math::IsZero(step)) { throw std::invalid_argument("arange: step must be nonzero"); }

  if (step > 0 ? !(start < stop) : !(start > stop)) { return container_type<value_type>(0); }
  std::size_t size = 0;
  if constexpr (std::is_integral_v<value_type>) {
    // Count integral samples exactly even across the signed range
    using Unsigned = std::make_unsigned_t<decltype(+start)>;
    const bool descending = step < 0;
    const Unsigned distance = descending ? Unsigned(start) - Unsigned(stop) : Unsigned(stop) - Unsigned(start);
    const Unsigned stride = descending ? Unsigned(0) - Unsigned(step) : Unsigned(step);
    const Unsigned count = distance / stride + Unsigned(distance % stride != 0);
    if (std::cmp_greater(count, std::numeric_limits<std::size_t>::max())) {
      throw std::length_error("arange: output size is not representable");
    }
    size = static_cast<std::size_t>(count);
  } else {
    const long double span =
        (static_cast<long double>(stop) - static_cast<long double>(start)) / static_cast<long double>(step);
    const long double count = std::max(1.0L, std::ceil(span));
    if (!std::isfinite(count) || count >= std::ldexp(1.0L, std::numeric_limits<std::size_t>::digits)) {
      throw std::length_error("arange: output size is not representable");
    }
    size = static_cast<std::size_t>(count);
    // Exclude samples rounded onto or beyond the open boundary
    while (size > 0) {
      const value_type last = std::fma(step, static_cast<value_type>(size - 1), start);
      if (step > 0 ? last < stop : last > stop) { break; }
      --size;
    }
  }
  container_type<value_type> v(size);
  if constexpr (std::is_integral_v<value_type>) {
    value_type current = start;
    for (std::size_t i = 0; i < size; ++i) {
      v[i] = current;
      if (i + 1 == size) { continue; }
      bool overflow = step > 0 && current > std::numeric_limits<value_type>::max() - step;
      if constexpr (std::is_signed_v<value_type>) {
        overflow = overflow || (step < 0 && current < std::numeric_limits<value_type>::min() - step);
      }
      if (overflow) { throw std::overflow_error("arange: element is not representable"); }
      current += step;
    }
  } else {
    for (std::size_t i = 0; i < size; ++i) { v[i] = std::fma(step, static_cast<value_type>(i), start); }
  }
  return v;
}

// Execute independent nodes with deterministic indexed output
template <typename Function>
void ParallelFor(const std::size_t count, const Function &function, const std::string_view progress_label = {}) {
  if (count == 0) { return; }
  const bool show_progress = !progress_label.empty();
  const bool terminal      = aux::IsTerminal();
  if (show_progress) {
    std::cout << "  " << progress_label << std::endl;
    if (terminal) {
      aux::PrintProgress(0.0);
    } else {
      std::cout << "  " << progress_label << " [0%]" << std::endl;
    }
  }
  const std::size_t thread_count =
      std::min<std::size_t>(count, std::max<std::size_t>(1, std::thread::hardware_concurrency()));
  std::atomic<std::size_t>  next{0};
  std::atomic<std::size_t>  completed{0};
  std::atomic<unsigned int> displayed_percent{0};
  std::atomic<bool>         failed{false};
  std::mutex                error_mutex;
  std::mutex                progress_mutex;
  std::exception_ptr        error;
  // Stop and join constructed workers if a later thread cannot be created
  std::vector<std::jthread> workers;
  workers.reserve(thread_count);
  for (std::size_t thread = 0; thread < thread_count; ++thread) {
    workers.emplace_back([&](const std::stop_token stop) {
      while (!stop.stop_requested() && !failed.load(std::memory_order_acquire)) {
        const std::size_t index = next.fetch_add(1, std::memory_order_relaxed);
        if (index >= count) { return; }
        try {
          function(index);
          if (show_progress) {
            const std::size_t  done    = completed.fetch_add(1, std::memory_order_relaxed) + 1;
            const double       ratio   = done / static_cast<double>(count);
            const auto         percent = static_cast<unsigned int>(100.0 * ratio);
            const unsigned int step    = terminal ? 1U : 10U;
            if (percent >= displayed_percent.load(std::memory_order_relaxed) + step || percent == 100U) {
              std::lock_guard<std::mutex> lock(progress_mutex);
              if (percent >= displayed_percent.load(std::memory_order_relaxed) + step || percent == 100U) {
                displayed_percent.store(percent, std::memory_order_relaxed);
                if (terminal) {
                  aux::PrintProgress(ratio);
                } else {
                  std::cout << "  " << progress_label << " [" << percent << "%]" << std::endl;
                }
              }
            }
          }
        } catch (...) {
          std::lock_guard<std::mutex> lock(error_mutex);
          if (!error) { error = std::current_exception(); }
          failed.store(true, std::memory_order_release);
        }
      }
    });
  }
  for (auto &worker : workers) { worker.join(); }
  if (show_progress) {
    if (!error && terminal) { aux::PrintProgress(1.0); }
    if (terminal) { aux::ClearProgress(); }
  }
  if (error) { std::rethrow_exception(error); }
  if (show_progress) { std::cout << "  " << progress_label << " [DONE]" << std::endl; }
}


}  // namespace math
}  // namespace gra

#endif
