// Minimal tensor class
//
// Example: (rank-6 tensor with dim-4 per dimension)
//
// MTensor<double> tensor({4,4,4,4,4,4}, 0.0);
// tensor({0,3,2,0,1,2}) = 1.0;  // Brace indices
// tensor(0,3,2,0,1,2)   = 1.0;  // Variadic indices
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MTENSOR_H
#define MTENSOR_H

#include <algorithm>
#include <array>
#include <cstddef>
#include <initializer_list>
#include <limits>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

namespace gra {

namespace detail {

// C++14-compatible conjunction over integral index types
template <typename... Ts> struct all_integral;

template <> struct all_integral<> : std::true_type {};

template <typename T, typename... Ts>
struct all_integral<T, Ts...>
    : std::integral_constant<
          bool,
          std::is_integral<typename std::decay<T>::type>::value &&
              all_integral<Ts...>::value> {};

} // namespace detail

template <typename T> class MTensor {
public:
  MTensor() noexcept : volume(0), capacity(0), data(nullptr) {}

  // Construct dimensions from a vector
  MTensor(const std::vector<std::size_t> &newdim)
      : volume(0), capacity(0), data(nullptr) {
    Allocate(newdim);

    if (data != nullptr) {
      std::fill_n(data, volume, T());
    }
  }

  // Construct dimensions from a brace list without a temporary vector
  MTensor(std::initializer_list<std::size_t> newdim)
      : volume(0), capacity(0), data(nullptr) {
    Allocate(newdim);

    if (data != nullptr) {
      std::fill_n(data, volume, T());
    }
  }

  // Construct vector dimensions and fill with a value
  MTensor(const std::vector<std::size_t> &newdim, T value)
      : volume(0), capacity(0), data(nullptr) {
    Allocate(newdim);

    if (data != nullptr) {
      std::fill_n(data, volume, value);
    }
  }

  // Construct brace-list dimensions and fill with a value
  MTensor(std::initializer_list<std::size_t> newdim, T value)
      : volume(0), capacity(0), data(nullptr) {
    Allocate(newdim);

    if (data != nullptr) {
      std::fill_n(data, volume, value);
    }
  }

  ~MTensor() { delete[] data; }

  // Copy constructor
  MTensor(const MTensor &a)
      : dim(a.dim),
        strides(a.strides),
        volume(a.volume),
        capacity(a.volume),
        data(nullptr) {

    if (volume != 0) {
      data = new T[volume];
      Copy(a);
    }
  }

  // Copy values while reusing allocated storage when possible
  MTensor &operator=(const MTensor &rhs) {
    if (this != &rhs) {
      if (rhs.volume == 0) {
        auto newdim     = rhs.dim;
        auto newstrides = rhs.strides;
        dim.swap(newdim);
        strides.swap(newstrides);
        volume = 0;
        return *this;
      }
      ReSize(rhs.dim);
      Copy(rhs);
    }

    return *this;
  }

  // Move constructor
  MTensor(MTensor &&other) noexcept
      : dim(std::move(other.dim)),
        strides(std::move(other.strides)),
        volume(other.volume),
        capacity(other.capacity),
        data(other.data) {

    other.volume   = 0;
    other.capacity = 0;
    other.data     = nullptr;
  }

  // Move assignment
  MTensor &operator=(MTensor &&other) noexcept {
    if (this != &other) {
      delete[] data;

      dim      = std::move(other.dim);
      strides  = std::move(other.strides);
      volume   = other.volume;
      capacity = other.capacity;
      data     = other.data;

      other.volume   = 0;
      other.capacity = 0;
      other.data     = nullptr;
    }

    return *this;
  }

  // Index tensor elements with a vector of coordinates
  T &operator()(const std::vector<std::size_t> &ind) {
    return data[Index(ind)];
  }

  const T &operator()(const std::vector<std::size_t> &ind) const {
    return data[Index(ind)];
  }
  
  // Index tensor({i,j,...}) without allocating a temporary vector
  T &operator()(std::initializer_list<std::size_t> ind) {
    return data[Index(ind)];
  }

  const T &operator()(std::initializer_list<std::size_t> ind) const {
    return data[Index(ind)];
  }

  // Index tensor(i,j,...) with integral coordinates stored on the stack
  template <
      typename... Indices,
      typename std::enable_if<
          (sizeof...(Indices) > 0) &&
              detail::all_integral<Indices...>::value,
          int>::type = 0>
  T &operator()(Indices... indices) {

    const std::array<std::size_t, sizeof...(Indices)> ind{{
        static_cast<std::size_t>(indices)...}};

    return data[Index(ind)];
  }

  template <
      typename... Indices,
      typename std::enable_if<
          (sizeof...(Indices) > 0) &&
              detail::all_integral<Indices...>::value,
          int>::type = 0>
  const T &operator()(Indices... indices) const {

    const std::array<std::size_t, sizeof...(Indices)> ind{{
        static_cast<std::size_t>(indices)...}};

    return data[Index(ind)];
  }

  // SIZE INFORMATION

  std::size_t size(std::size_t ind) const {
    return dim.at(ind);
  }

  std::size_t rank() const noexcept {
    return dim.size();
  }

  bool empty() const noexcept {
    return volume == 0;
  }

  // Number of logical tensor elements
  std::size_t elements() const noexcept {
    return volume;
  }

  // Access contiguous row-major storage, or nullptr when logically empty
  T *raw_data() noexcept {
    return (volume == 0) ? nullptr : data;
  }

  const T *raw_data() const noexcept {
    return (volume == 0) ? nullptr : data;
  }

private:

  void ValidateIndexRank(std::size_t index_rank) const {

    if (volume == 0) {
      throw std::invalid_argument(
          "MTensor:: Error: Cannot index an empty tensor");
    }

    if (index_rank != dim.size()) {
      throw std::invalid_argument(
          "MTensor:: Error: Index vector with rank = " +
          std::to_string(index_rank) +
          " c.f. Tensor has rank " +
          std::to_string(dim.size()));
    }
  }

  // Dynamic std::vector indexing
  std::size_t Index(
      const std::vector<std::size_t> &ind) const {

    ValidateIndexRank(ind.size());

    std::size_t offset = 0;

    for (std::size_t i = 0; i < ind.size(); ++i) {

      if (ind[i] >= dim[i]) {
        throw std::invalid_argument(
            "MTensor:: Error: Input index " +
            std::to_string(i) +
            " over bounds: " +
            std::to_string(ind[i]) +
            " >= " +
            std::to_string(dim[i]));
      }

      offset += strides[i] * ind[i];
    }

    return offset;
  }

  // initializer_list indexing
  std::size_t Index(
      std::initializer_list<std::size_t> ind) const {

    ValidateIndexRank(ind.size());

    std::size_t offset = 0;
    std::size_t axis   = 0;

    for (const std::size_t index : ind) {

      if (index >= dim[axis]) {
        throw std::invalid_argument(
            "MTensor:: Error: Input index " +
            std::to_string(axis) +
            " over bounds: " +
            std::to_string(index) +
            " >= " +
            std::to_string(dim[axis]));
      }

      offset += strides[axis] * index;
      ++axis;
    }

    return offset;
  }

  // Compile-time-sized stack index
  template <std::size_t N>
  std::size_t Index(
      const std::array<std::size_t, N> &ind) const {

    ValidateIndexRank(N);

    std::size_t offset = 0;

    for (std::size_t i = 0; i < N; ++i) {

      if (ind[i] >= dim[i]) {
        throw std::invalid_argument(
            "MTensor:: Error: Input index " +
            std::to_string(i) +
            " over bounds: " +
            std::to_string(ind[i]) +
            " >= " +
            std::to_string(dim[i]));
      }

      offset += strides[i] * ind[i];
    }

    return offset;
  }

  void Copy(const MTensor &a) {

    if (volume == 0) {
      return;
    }

    std::copy_n(a.data, volume, data);
  }

  // Reshape with storage reuse, preparing allocation and metadata before modifying the tensor
  void ReSize(
      const std::vector<std::size_t> &newdim) {
    if (dim == newdim && volume != 0) { return; }

    // Prepare first
    std::vector<std::size_t> newdims = newdim;

    const std::size_t newvolume =
        Prod(newdims);

    std::vector<std::size_t> newstrides =
        ComputeStrides(newdims, newvolume);

    T *newdata = nullptr;

    // Allocate only if existing capacity is insufficient
    if (newvolume > capacity) {
      newdata = new T[newvolume];
    }

    // Commit after successful preparation/allocation
    if (newdata != nullptr) {
      delete[] data;

      data     = newdata;
      capacity = newvolume;
    }

    dim.swap(newdims);
    strides.swap(newstrides);
    volume = newvolume;
  }

  // Compute volume with overflow checks, treating rank zero as one scalar
  static std::size_t Prod(
      const std::vector<std::size_t> &x) {

    std::size_t product = 1;

    for (const std::size_t value : x) {

      if (value == 0) {
        return 0;
      }

      if (product >
          std::numeric_limits<std::size_t>::max() / value) {

        throw std::overflow_error(
            "MTensor:: Error: Tensor volume overflow");
      }

      product *= value;
    }

    return product;
  }

  static std::size_t Prod(
      std::initializer_list<std::size_t> x) {

    std::size_t product = 1;

    for (const std::size_t value : x) {

      if (value == 0) {
        return 0;
      }

      if (product >
          std::numeric_limits<std::size_t>::max() / value) {

        throw std::overflow_error(
            "MTensor:: Error: Tensor volume overflow");
      }

      product *= value;
    }

    return product;
  }

  // Allocation
  void Allocate(
      const std::vector<std::size_t> &newdim) {

    std::vector<std::size_t> newdims =
        newdim;

    const std::size_t newvolume =
        Prod(newdims);

    std::vector<std::size_t> newstrides =
        ComputeStrides(newdims, newvolume);

    T *newdata =
        (newvolume != 0)
            ? new T[newvolume]
            : nullptr;

    dim.swap(newdims);
    strides.swap(newstrides);

    volume   = newvolume;
    capacity = newvolume;
    data     = newdata;
  }

  void Allocate(
      std::initializer_list<std::size_t> newdim) {

    const std::size_t newvolume =
        Prod(newdim);

    std::vector<std::size_t> newdims(
        newdim.begin(),
        newdim.end());

    std::vector<std::size_t> newstrides =
        ComputeStrides(newdims, newvolume);

    T *newdata =
        (newvolume != 0)
            ? new T[newvolume]
            : nullptr;

    dim      = std::move(newdims);
    strides  = std::move(newstrides);
    volume   = newvolume;
    capacity = newvolume;
    data     = newdata;
  }

  // Strides

  // Row-major:
  //
  // offset =
  //
  //   i0 * (d1*d2*...*dN)
  // + i1 * (d2*d3*...*dN)
  // + ...
  // + iN
  //
  static std::vector<std::size_t>
  ComputeStrides(
      const std::vector<std::size_t> &newdim,
      const std::size_t newvolume) {

    std::vector<std::size_t> out(
        newdim.size(), 1);

    if (newdim.empty()) {
      return out;
    }

    // No valid element can ever be indexed if one dimension is zero
    // Avoid computing meaningless suffix products in that case
    if (newvolume == 0) {
      std::fill(out.begin(), out.end(), 0);
      return out;
    }

    for (std::size_t i = newdim.size();
         i-- > 1;) {

      // Because the complete non-zero volume has already been checked for
      // overflow, every suffix product is finite as well
      out[i - 1] =
          out[i] * newdim[i];
    }

    return out;
  }

  std::vector<std::size_t> dim;
  std::vector<std::size_t> strides;

  // Logical number of elements
  std::size_t volume;

  // Allocated number of elements. May exceed volume after a reshape to a
  // smaller tensor, allowing subsequent assignments to reuse storage
  std::size_t capacity;

  T *data;
};

} // namespace gra

#endif
