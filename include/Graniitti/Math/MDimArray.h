// Multidimensional fixed size arrays via template alias technique
//
//
// Usage:
// MultiDimArray<int, 3, 2> arr {1, 2, 3, 4, 5, 6};
// assert(arr[1][1] == 4);
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDIMARRAY_H
#define MDIMARRAY_H

#include <array>
#include <cstddef>
#include <complex>

namespace gra {
// Nest the rightmost index innermost for row-major storage
//
template <typename T, std::size_t D1, std::size_t D2, std::size_t... DN>
struct GetArray {
  using type = std::array<typename GetArray<T, D2, DN...>::type, D1>;
};

template <typename T, std::size_t D1, std::size_t D2>
struct GetArray<T, D1, D2> {
  using type = std::array<std::array<T, D2>, D1>;
};

template <typename T, std::size_t D1, std::size_t D2, std::size_t... DN>
using MDimArray = typename GetArray<T, D1, D2, DN...>::type;

}  // namespace gra

#endif
