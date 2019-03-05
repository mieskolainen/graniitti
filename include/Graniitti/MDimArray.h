// Multidimensional fixed size arrays via template alias technique
//
//
// Usage:
// MultiDimArray<int, 3, 2> arr {1, 2, 3, 4, 5, 6};
// assert(arr[1][1] == 4);
// 
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MDIMARRAY_H
#define MDIMARRAY_H

#include <complex>
#include <array>


namespace gra {

// Note the index order below, in order to get the row major format
// 
template <typename T, size_t D1, size_t D2, size_t... DN>
struct GetArray {
    using type = std::array<typename GetArray<T, D2, DN...>::type, D1>;
};

template <typename T, size_t D1, size_t D2>
struct GetArray<T, D1, D2> {
    using type = std::array<std::array<T, D2>, D1>;
};

template <typename T, size_t D1, size_t D2, size_t... DN>
using MDimArray = typename GetArray<T, D1, D2, DN...>::type;


} // gra namespace ends

#endif