// Floating point classification helpers
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFLOAT_H
#define MFLOAT_H

// C++
#include <cmath>
#include <compare>
#include <complex>
#include <concepts>
#include <type_traits>

namespace gra {
namespace math {

namespace detail {

// Identify one standard complex scalar type
template <typename T> struct IsComplex : std::false_type {};

template <typename T>
struct IsComplex<std::complex<T>> : std::true_type {};

} // namespace detail

// Compute whether one floating point value is exactly zero
// result = [fpclassify(value) = FP_ZERO]
template <typename T>
requires std::is_floating_point_v<T>
inline bool IsZero(const T value) {
  return std::fpclassify(value) == FP_ZERO;
}

// Compute whether one integral value is exactly zero
template <typename T>
requires std::is_integral_v<T>
constexpr bool IsZero(const T value) {
  return value == T{0};
}

// Compute whether both components of one complex value are exactly zero
// result = [Re(value) = 0 and Im(value) = 0]
template <typename T>
inline bool IsZero(const std::complex<T> &value) {
  return IsZero(value.real()) && IsZero(value.imag());
}

// Preserve exact zero support for custom equality comparable scalar types
template <typename T>
requires(!std::is_arithmetic_v<std::remove_cv_t<T>> &&
         !detail::IsComplex<std::remove_cv_t<T>>::value &&
         requires(const T &value) {
           { value == T{} } -> std::convertible_to<bool>;
         })
inline bool IsZero(const T &value) {
  return value == T{};
}

// Compute whether two floating point values are exactly numerically equal
// result = [(first <=> second) is equal]
template <typename T>
requires std::is_floating_point_v<T>
inline bool IsExactEqual(const T first, const T second) {
  return std::is_eq(first <=> second);
}

} // namespace math
} // namespace gra

#endif // MFLOAT_H
