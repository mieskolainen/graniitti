// Reusable templated matrix and container linear algebra
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATRIX_H
#define MMATRIX_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdio>
#include <initializer_list>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <ranges>
#include <span>
#include <stdexcept>
#include <string>
#include <string_view>
#include <type_traits>
#include <utility>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Math/MFloat.h"

// Eigen
#include <Eigen/Dense>

namespace gra {

// Select whether a Cholesky factor requires strict or semidefinite pivots
enum class CholeskyMode { StrictPositiveDefinite, AllowSemidefinite };

// Select one factor of an ordered bipartite tensor product
enum class TensorFactor { First, Second };

// Store the numerical and retained ranks of one Moore-Penrose inverse
struct PseudoInverseDiagnostics {
  std::size_t numerical_rank   = 0;
  std::size_t retained_rank    = 0;
  double      condition_number = 0.0;
};

// Dynamically sized matrix
template <typename T>
class MMatrix {
 public:
  MMatrix() : rows(0), cols(0) { data = nullptr; }

  MMatrix(std::size_t r, std::size_t c) {
    Allocate(r, c);
    if (data != nullptr) {
      std::fill(data, data + rows * cols, T());  // Zero/default initialization
    }
  }

  // Construct a square matrix with default initialized entries
  explicit MMatrix(std::size_t n) : MMatrix(n, n) {}
  MMatrix(std::size_t r, std::size_t c, T value) {
    Allocate(r, c);
    if (data != nullptr) { std::fill(data, data + rows * cols, T(value)); }
  }

  // Special init
  MMatrix(std::size_t r, std::size_t c, const std::string &special) {
    if (special != "eye" && special != "minkowski") {
      throw std::invalid_argument("MMatrix: Unknown initialization string:" + special);
    }

    Allocate(r, c);

    if (special == "eye") {
      Identity();
    } else {
      Minkowski();
    }
  }

  // Set matrix to identity [diagonal = 1, otherwise 0]
  void Identity() {
    const T zero = T(0);
    const T one  = T(1);
    for (std::size_t i = 0; i < rows; ++i) {
      T *row_data = data + cols * i;
      for (std::size_t j = 0; j < cols; ++j) { row_data[j] = (i == j) ? one : zero; }
    }
  }
  // Set matrix to Minkowski metric (+,-,-,-)
  void Minkowski() {
    const T zero  = T(0);
    const T plus  = T(1);
    const T minus = T(-1);
    for (std::size_t i = 0; i < rows; ++i) {
      T *row_data = data + cols * i;
      for (std::size_t j = 0; j < cols; ++j) {
        const T diag = (i > 0 && j > 0) ? minus : plus;
        row_data[j]  = (i == j) ? diag : zero;
      }
    }
  }

  // For initializing with a = {row0-vector, row1-vector, ...}
  // where each row is std::vector<T>.
  //
  // This constructor is need such that plain nested brace
  // initialization prefers the non-template initializer_list overload below
  //
  // Row-major order!
  template <typename U = T, typename std::enable_if<std::is_same<U, T>::value, int>::type = 0>
  MMatrix(std::initializer_list<std::vector<U>> list) {
    rows = list.size();
    cols = GetRectangularCols(list);
    Allocate(rows, cols);

    std::size_t i = 0;
    for (const auto &v : list) {
      for (std::size_t j = 0; j < cols; ++j) { data[cols * i + j] = v[j]; }
      ++i;
    }
  }

  // Initialize one matrix from fixed size row arrays
  template <std::size_t N>
  MMatrix(std::initializer_list<std::array<T, N>> list) {
    rows = list.size();
    cols = N;
    Allocate(rows, cols);

    std::size_t i = 0;
    for (const auto &v : list) {
      for (std::size_t j = 0; j < N; ++j) { data[cols * i + j] = v[j]; }
      ++i;
    }
  }

  // For initializing directly with nested brace lists:
  //
  //   MMatrix<double> A{{1.0, 2.0}, {3.0, 4.0}};
  //
  // Row-major order!
  MMatrix(std::initializer_list<std::initializer_list<T>> list) {
    rows = list.size();
    cols = GetRectangularCols(list);
    Allocate(rows, cols);

    std::size_t i = 0;
    for (const auto &v : list) {
      std::size_t j = 0;
      for (const auto &w : v) {
        data[cols * i + j] = w;
        ++j;
      }
      ++i;
    }
  }

  ~MMatrix() { delete[] data; }

  // Copy constructor
  MMatrix(const MMatrix &a) {
    Allocate(a.rows, a.cols);
    Copy(a);
  }

  // Assignment operator
  MMatrix &operator=(const MMatrix &rhs) {
    if (this != &rhs) {
      if (rows != rhs.rows || cols != rhs.cols) { Resize(rhs.rows, rhs.cols); }
      Copy(rhs);
    }
    return *this;
  }

  // Resize this matrix without initializing the new entries
  void Resize(const std::size_t r, const std::size_t c) {
    if (rows == r && cols == c) { return; }
    T                *new_data = nullptr;
    const std::size_t count    = CheckedElementCount(r, c);
    if (count != 0) { new_data = new T[count]; }
    delete[] data;
    rows = r;
    cols = c;
    data = new_data;
  }

  // Move constructor
  MMatrix(MMatrix &&other) noexcept : rows(other.rows), cols(other.cols), data(other.data) {
    other.data = nullptr;  // Leave the other object in a valid state
    other.rows = 0;
    other.cols = 0;
  }

  // Move assignment operator
  MMatrix &operator=(MMatrix &&other) noexcept {
    if (this != &other) {
      delete[] data;
      rows = other.rows;
      cols = other.cols;
      data = other.data;

      // Leave 'other' in a valid empty state
      other.data = nullptr;
      other.rows = 0;
      other.cols = 0;
    }
    return *this;
  }

  // For indexing with [i][j]
  // Row-major order!
  T *operator[](const std::size_t &row) {
    if (row >= rows) { throw std::out_of_range("MMatrix:: row index over matrix dimensions"); }
    return data + cols * row;
  }
  const T *operator[](const std::size_t &row) const {
    if (row >= rows) { throw std::out_of_range("MMatrix:: row index over matrix dimensions"); }
    return data + cols * row;
  }

  // Row-major order!
  T &operator()(std::size_t i, std::size_t j) {
    if (i >= rows || j >= cols) { throw std::out_of_range("MMatrix:: Index over matrix dimensions!"); }
    return data[cols * i + j];
  }
  const T &operator()(std::size_t i, std::size_t j) const {
    if (i >= rows || j >= cols) { throw std::out_of_range("MMatrix:: Index over matrix dimensions!"); }
    return data[cols * i + j];
  }

  // ------------------------------------------------------------------
  // Add to the left
  // A_{ij} <- A_{ij} + B_{ij}
  MMatrix &operator+=(const MMatrix &rhs) {
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      throw std::invalid_argument("MMatrix:: operator+=: Dimension mismatch");
    }
    for (std::size_t i = 0; i < rows * cols; ++i) { data[i] += rhs.data[i]; }
    return *this;
  }
  // Subtract to the left
  // A_{ij} <- A_{ij} - B_{ij}
  MMatrix &operator-=(const MMatrix &rhs) {
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      throw std::invalid_argument("MMatrix:: operator-=: Dimension mismatch");
    }
    for (std::size_t i = 0; i < rows * cols; ++i) { data[i] -= rhs.data[i]; }
    return *this;
  }

  // ------------------------------------------------------------------

  // Compute negated matrix
  // output_{ij} = -A_{ij}
  MMatrix operator-() const {
    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = -data[i]; }
    return out;
  }

  // Add two matrices
  // output_{ij} = A_{ij} + B_{ij}
  MMatrix operator+(const MMatrix &rhs) const {
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      const std::string str =
          "MMatrix:: Matrix + Matrix with invalid dimensions: "
          "Left matrix (" +
          std::to_string(rows) + "x" + std::to_string(cols) +
          "), "
          "Right matrix (" +
          std::to_string(rhs.size_row()) + "x" + std::to_string(rhs.size_col()) + ")";
      throw std::invalid_argument(str);
    }

    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = data[i] + rhs.data[i]; }

    return out;
  }

  // Subtract two matrices
  // output_{ij} = A_{ij} - B_{ij}
  MMatrix operator-(const MMatrix &rhs) const {
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      const std::string str =
          "MMatrix:: Matrix - Matrix with invalid dimensions: "
          "Left matrix (" +
          std::to_string(rows) + "x" + std::to_string(cols) +
          "), "
          "Right matrix (" +
          std::to_string(rhs.size_row()) + "x" + std::to_string(rhs.size_col()) + ")";
      throw std::invalid_argument(str);
    }

    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = data[i] - rhs.data[i]; }

    return out;
  }
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // We do not allow addition or substraction by scalar, for dimensional
  // "safety" reasons (can lead to accidental results), thus only multiplicative
  // operations

  // Multiply by a scalar
  // output_{ij} = A_{ij} c
  MMatrix operator*(const T &rhs) const {
    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = data[i] * rhs; }

    return out;
  }

  // Multiply this matrix by one compatible scalar type
  // output_{ij} = A_{ij} c
  template <typename Scalar>
  auto operator*(const Scalar &rhs) const -> MMatrix<decltype(std::declval<T>() * std::declval<Scalar>())> {
    using R        = decltype(std::declval<T>() * std::declval<Scalar>());
    MMatrix<R> out = MMatrix<R>::Uninitialized(rows, cols);
    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = data[i] * rhs; }
    return out;
  }

  // Divide by a scalar
  // output_{ij} = A_{ij}/c
  MMatrix operator/(const T &rhs) const {
    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = data[i] / rhs; }

    return out;
  }

  // Multiply this matrix by a scalar
  // A_{ij} <- A_{ij} c
  MMatrix &operator*=(const T &rhs) {
    for (std::size_t i = 0; i < rows * cols; ++i) { data[i] *= rhs; }
    return *this;
  }

  // Accumulate A <- A + c B, where A is this matrix, B is rhs and c is scale
  template <typename Scalar>
  void AddScaled(const MMatrix &rhs, const Scalar &scale) {
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      throw std::invalid_argument("MMatrix::AddScaled: Dimension mismatch");
    }
    for (std::size_t i = 0; i < rows * cols; ++i) { data[i] += rhs.data[i] * scale; }
  }
  // ------------------------------------------------------------------

  // Multiply this matrix by a vector with a compatible scalar type
  // y_i = sum_j A_{ij} x_j
  template <typename U>
  auto operator*(const std::vector<U> &rhs) const -> std::vector<decltype(std::declval<T>() * std::declval<U>())> {
    using R = decltype(std::declval<T>() * std::declval<U>());
    if (rhs.size() != cols) {
      const std::string err = "MMatrix:: mat * vec product with invalid dimensions: " + std::to_string(cols) + " x " +
                              std::to_string(rhs.size());
      throw std::invalid_argument(err);
    }

    std::vector<R> out(rows, R(0));
    for (std::size_t i = 0; i < rows; ++i) {
      const T *lhs_row = data + cols * i;
      for (std::size_t j = 0; j < cols; ++j) {
        const T lhs = lhs_row[j];
        if (IsZero(lhs)) { continue; }
        out[i] += lhs * rhs[j];
      }
    }
    return out;
  }

  // Accumulate target <- target + scale * this * source without allocation
  // y_i <- y_i + c sum_j A_{ij} x_j
  template <typename U, typename V, typename Scalar>
  void MultiplyAdd(const std::span<const U> source, const std::span<V> target, const Scalar &scale) const {
    if (source.size() != cols || target.size() != rows) {
      throw std::invalid_argument("MMatrix::MultiplyAdd: vector dimensions disagree");
    }
    for (std::size_t row = 0; row < rows; ++row) {
      const T *matrix_row  = data + cols * row;
      V        accumulator = target[row];
      bool     has_nonzero = false;
      for (std::size_t col = 0; col < cols; ++col) {
        const T value = matrix_row[col];
        if (IsZero(value)) { continue; }
        accumulator += scale * value * source[col];
        has_nonzero = true;
      }
      if (has_nonzero) { target[row] = accumulator; }
    }
  }

  // Multiply one square matrix by one fixed size vector
  // y_i = sum_j A_{ij} x_j
  template <std::size_t N>
  std::array<T, N> operator*(const std::array<T, N> &rhs) const {
    if (rows != N || cols != N) {
      const std::string err = "MMatrix:: mat * array product with invalid dimensions: " + std::to_string(rows) + " x " +
                              std::to_string(cols) + " versus " + std::to_string(N);
      throw std::invalid_argument(err);
    }

    std::array<T, N> out{};
    for (std::size_t i = 0; i < N; ++i) {
      const T *lhs_row = data + cols * i;
      for (std::size_t j = 0; j < N; ++j) {
        const T lhs = lhs_row[j];
        if (IsZero(lhs)) { continue; }
        out[i] += lhs * rhs[j];
      }
    }
    return out;
  }

  // Multiply one dynamic row vector by this matrix without conjugation
  // y_j = sum_i x_i A_{ij}
  template <typename Container>
  auto LeftMultiply(const Container &row) const {
    if (row.size() != rows) { throw std::invalid_argument("MMatrix::LeftMultiply: vector dimensions disagree"); }
    using U = typename Container::value_type;
    using R = decltype(std::declval<U>() * std::declval<T>());
    std::vector<R> output(cols, R{});
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) { output[j] += row[i] * data[i * cols + j]; }
    }
    return output;
  }

  // Multiply one fixed row vector by this square matrix without conjugation
  // y_j = sum_i x_i A_{ij}
  template <typename U, std::size_t N>
  auto LeftMultiply(const std::array<U, N> &row) const {
    if (rows != N || cols != N) { throw std::invalid_argument("MMatrix::LeftMultiply: array dimensions disagree"); }
    using R = decltype(std::declval<U>() * std::declval<T>());
    std::array<R, N> output{};
    for (std::size_t i = 0; i < N; ++i) {
      for (std::size_t j = 0; j < N; ++j) { output[j] += row[i] * data[i * cols + j]; }
    }
    return output;
  }

  // Compute one nonconjugating row-matrix-column bilinear form
  // output = sum_{ij} x_i A_{ij} y_j
  template <typename FirstContainer, typename SecondContainer>
  auto BilinearForm(const FirstContainer &left, const SecondContainer &right) const {
    if (left.size() != rows || right.size() != cols) {
      throw std::invalid_argument("MMatrix::BilinearForm: vector dimensions disagree");
    }
    using First  = typename FirstContainer::value_type;
    using Second = typename SecondContainer::value_type;
    using R      = decltype(std::declval<First>() * std::declval<T>() * std::declval<Second>());
    R result{};
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) { result += left[i] * data[i * cols + j] * right[j]; }
    }
    return result;
  }

  // Compute one nonconjugating bilinear form restricted by independent masks
  // output = sum_{i in rows,j in cols} x_i A_{ij} y_j
  template <typename FirstContainer, typename SecondContainer, typename RowMask, typename ColumnMask>
  auto MaskedBilinearForm(const FirstContainer &left, const SecondContainer &right, const RowMask &row_mask,
                          const ColumnMask &column_mask) const {
    if (left.size() != rows || right.size() != cols || row_mask.size() != rows || column_mask.size() != cols) {
      throw std::invalid_argument("MMatrix::MaskedBilinearForm: vector dimensions disagree");
    }
    using First  = typename FirstContainer::value_type;
    using Second = typename SecondContainer::value_type;
    using R      = decltype(std::declval<First>() * std::declval<T>() * std::declval<Second>());
    R result{};
    for (std::size_t row = 0; row < rows; ++row) {
      if (!row_mask[row]) { continue; }
      for (std::size_t col = 0; col < cols; ++col) {
        if (column_mask[col]) { result += left[row] * data[row * cols + col] * right[col]; }
      }
    }
    return result;
  }

  // Compute one Hermitian matrix element left dagger times this times right
  // output = sum_{ij} x_i^* A_{ij} y_j
  template <typename LeftContainer, typename RightContainer>
  auto MatrixElement(const LeftContainer &left, const RightContainer &right) const {
    if (left.size() != rows || right.size() != cols) {
      throw std::invalid_argument("MMatrix::MatrixElement: vector dimensions disagree");
    }
    using Left  = typename LeftContainer::value_type;
    using Right = typename RightContainer::value_type;
    using R     = decltype(Conjugate(std::declval<Left>()) * std::declval<T>() * std::declval<Right>());
    R result{};
    for (std::size_t i = 0; i < rows; ++i) {
      const auto left_value = Conjugate(left[i]);
      for (std::size_t j = 0; j < cols; ++j) { result += left_value * data[i * cols + j] * right[j]; }
    }
    return result;
  }

  // Accumulate one scaled outer product into this matrix
  // A_{ij} <- A_{ij} + c x_i y_j
  template <typename FirstContainer, typename SecondContainer, typename Scalar>
  void AddOuterProduct(const FirstContainer &first, const SecondContainer &second, const Scalar &scale) {
    if (!std::cmp_equal(first.size(), rows) || !std::cmp_equal(second.size(), cols)) {
      throw std::invalid_argument("MMatrix::AddOuterProduct: vector dimensions disagree");
    }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) { data[row * cols + col] += scale * first[row] * second[col]; }
    }
  }

  // Accumulate a Kronecker vector into one matrix column
  // A_{iN+j,k} <- A_{iN+j,k} + x_i y_j
  template <typename FirstContainer, typename SecondContainer>
  void AddColumnKroneckerProduct(const std::size_t column, const FirstContainer &first, const SecondContainer &second) {
    AddColumnKroneckerProduct(column, first, second, T(1));
  }

  // Accumulate one scaled Kronecker vector into one matrix column
  // A_{iN+j,k} <- A_{iN+j,k} + c x_i y_j
  template <typename FirstContainer, typename SecondContainer, typename Scalar>
  void AddColumnKroneckerProduct(const std::size_t column, const FirstContainer &first, const SecondContainer &second,
                                 const Scalar &scale) {
    if (column >= cols) { throw std::out_of_range("MMatrix::AddColumnKroneckerProduct: column outside matrix"); }
    if (!second.empty() && first.size() > std::numeric_limits<std::size_t>::max() / second.size()) {
      throw std::invalid_argument("MMatrix::AddColumnKroneckerProduct: dimension overflow");
    }
    if (!std::cmp_equal(first.size() * second.size(), rows)) {
      throw std::invalid_argument("MMatrix::AddColumnKroneckerProduct: dimensions disagree");
    }
    for (std::size_t first_index = 0; first_index < first.size(); ++first_index) {
      for (std::size_t second_index = 0; second_index < second.size(); ++second_index) {
        data[(first_index * second.size() + second_index) * cols + column] +=
            scale * first[first_index] * second[second_index];
      }
    }
  }

  // Permute rows by their destination indices without changing column order
  // output_{p(i),j} = A_{i,j}
  template <typename IndexContainer>
  MMatrix PermuteRows(const IndexContainer &destination_rows) const {
    ValidateDestinationRows(destination_rows, rows, "MMatrix::PermuteRows");
    MMatrix out = MMatrix::Uninitialized(rows, cols);
    if (cols == 0) { return out; }
    for (std::size_t row = 0; row < rows; ++row) {
      std::copy_n(data + row * cols, cols, out.data + static_cast<std::size_t>(destination_rows[row]) * cols);
    }
    return out;
  }

  // Matrix [this] * Matrix [rhs] multiplication
  // C_{ij} = sum_k A_{ik} B_{kj}
  // Optimized matrix multiplication with exact-zero skips for sparse helicity
  // matrices
  MMatrix operator*(const MMatrix &rhs) const {
    const std::size_t n = rows;
    const std::size_t m = cols;
    const std::size_t p = rhs.size_col();

    if (cols != rhs.size_row()) {
      const std::string str =
          "MMatrix:: Matrix * Matrix with invalid dimensions: "
          "Left matrix (" +
          std::to_string(rows) + "x" + std::to_string(cols) +
          "), "
          "Right matrix (" +
          std::to_string(rhs.size_row()) + "x" + std::to_string(rhs.size_col()) + ")";
      throw std::invalid_argument(str);
    }
    MMatrix<T> C(n, p, T(0));

    // Loop reordering to improve cache locality
    for (std::size_t i = 0; i < n; ++i) {
      T       *c_row   = C[i];
      const T *lhs_row = data + cols * i;
      for (std::size_t k = 0; k < m; ++k) {
        const T tmp = lhs_row[k];
        if (IsZero(tmp)) { continue; }
        const T *rhs_row = rhs.data + rhs.cols * k;
        for (std::size_t j = 0; j < p; ++j) { c_row[j] += tmp * rhs_row[j]; }
      }
    }

    return C;
  }

  // Evaluate C = c A B, where A is this matrix, B is rhs and c is scale
  // C_{ij} = c sum_k A_{ik} B_{kj}
  MMatrix MultiplyScaled(const MMatrix &rhs, const T &scale) const {
    MMatrix<T> C = (*this) * rhs;
    C *= scale;
    return C;
  }

  // Solve this matrix times X equals rhs by pivoted Gaussian elimination
  // AX = B
  MMatrix Solve(MMatrix rhs) const {
    if (rows != cols || rows != rhs.rows) { throw std::invalid_argument("MMatrix::Solve: incompatible dimensions"); }
    if (!IsFinite() || !rhs.IsFinite()) { throw std::invalid_argument("MMatrix::Solve: non-finite matrix"); }
    MMatrix           coefficients(*this);
    const std::size_t n             = rows;
    const std::size_t right_columns = rhs.cols;
    const double      pivot_floor   = 64.0 * std::numeric_limits<double>::epsilon();
    // Normalize each equation so that its units do not determine the pivot threshold
    for (std::size_t row = 0; row < n; ++row) {
      double scale = 0.0;
      for (std::size_t col = 0; col < n; ++col) {
        scale = std::max(scale, static_cast<double>(std::abs(coefficients(row, col))));
      }
      if (!(scale > 0.0) || !std::isfinite(scale)) {
        throw std::runtime_error("MMatrix::Solve: invalid matrix scale");
      }
      for (std::size_t col = 0; col < n; ++col) { coefficients(row, col) /= T(scale); }
      for (std::size_t col = 0; col < right_columns; ++col) { rhs(row, col) /= T(scale); }
    }
    if (!rhs.IsFinite()) { throw std::runtime_error("MMatrix::Solve: source scaling overflow"); }

    for (std::size_t pivot_column = 0; pivot_column < n; ++pivot_column) {
      std::size_t pivot_row = pivot_column;
      double      pivot_abs = std::abs(coefficients(pivot_column, pivot_column));
      for (std::size_t row = pivot_column + 1; row < n; ++row) {
        const double candidate = std::abs(coefficients(row, pivot_column));
        if (candidate > pivot_abs) {
          pivot_row = row;
          pivot_abs = candidate;
        }
      }
      if (!std::isfinite(pivot_abs) || !(pivot_abs > pivot_floor)) {
        throw std::runtime_error("MMatrix::Solve: singular matrix");
      }
      if (pivot_row != pivot_column) {
        for (std::size_t col = pivot_column; col < n; ++col) {
          std::swap(coefficients(pivot_column, col), coefficients(pivot_row, col));
        }
        for (std::size_t col = 0; col < right_columns; ++col) {
          std::swap(rhs(pivot_column, col), rhs(pivot_row, col));
        }
      }

      for (std::size_t row = pivot_column + 1; row < n; ++row) {
        const T factor                  = coefficients(row, pivot_column) / coefficients(pivot_column, pivot_column);
        coefficients(row, pivot_column) = T(0);
        for (std::size_t col = pivot_column + 1; col < n; ++col) {
          coefficients(row, col) -= factor * coefficients(pivot_column, col);
        }
        for (std::size_t col = 0; col < right_columns; ++col) { rhs(row, col) -= factor * rhs(pivot_column, col); }
      }
    }

    MMatrix solution = Uninitialized(n, right_columns);
    for (std::size_t reverse_row = n; reverse_row-- > 0;) {
      for (std::size_t col = 0; col < right_columns; ++col) {
        T value = rhs(reverse_row, col);
        for (std::size_t upper_col = reverse_row + 1; upper_col < n; ++upper_col) {
          value -= coefficients(reverse_row, upper_col) * solution(upper_col, col);
        }
        solution(reverse_row, col) = value / coefficients(reverse_row, reverse_row);
      }
    }
    if (!solution.IsFinite()) { throw std::runtime_error("MMatrix::Solve: non-finite solution"); }
    return solution;
  }

  // Compute the scaling and squaring Pade matrix exponential
  // exp(A) = sum_{k=0}^infinity A^k/k!
  MMatrix Exp() const {
    if (rows != cols) { throw std::invalid_argument("MMatrix::Exp requires a square matrix"); }
    if (rows == 0) { return MMatrix(); }
    if (IsDiagonal() && IsFinite()) {
      MMatrix result(rows, cols, T{});
      for (std::size_t index = 0; index < rows; ++index) {
        result.data[index * cols + index] = std::exp(data[index * cols + index]);
      }
      return result;
    }

    constexpr double theta13 = 5.371920351148152;
    const double     norm    = NormOne();
    const int        scaling =
        (norm <= theta13 || math::IsZero(norm)) ? 0 : static_cast<int>(std::ceil(std::log2(norm / theta13)));
    const double  scale    = std::ldexp(1.0, scaling);
    const MMatrix scaled   = (*this) / T(scale);
    const MMatrix squared  = scaled * scaled;
    const MMatrix fourth   = squared * squared;
    const MMatrix sixth    = fourth * squared;
    const MMatrix identity = IdentityMatrix(rows);

    constexpr std::array<double, 14> coefficient = {64764752532480000.0,
                                                    32382376266240000.0,
                                                    7771770303897600.0,
                                                    1187353796428800.0,
                                                    129060195264000.0,
                                                    10559470521600.0,
                                                    670442572800.0,
                                                    33522128640.0,
                                                    1323241920.0,
                                                    40840800.0,
                                                    960960.0,
                                                    16380.0,
                                                    182.0,
                                                    1.0};

    const MMatrix numerator_inner =
        sixth * (sixth * T(coefficient[13]) + fourth * T(coefficient[11]) + squared * T(coefficient[9])) +
        sixth * T(coefficient[7]) + fourth * T(coefficient[5]) + squared * T(coefficient[3]) +
        identity * T(coefficient[1]);
    const MMatrix numerator = scaled * numerator_inner;
    const MMatrix denominator =
        sixth * (sixth * T(coefficient[12]) + fourth * T(coefficient[10]) + squared * T(coefficient[8])) +
        sixth * T(coefficient[6]) + fourth * T(coefficient[4]) + squared * T(coefficient[2]) +
        identity * T(coefficient[0]);

    MMatrix result = (denominator - numerator).Solve(denominator + numerator);
    for (int step = 0; step < scaling; ++step) { result = result * result; }
    return result;
  }

  // Compute the mixed scalar Kronecker product with rhs
  // (A tensor B)_{(i,k),(j,l)} = A_{ij} B_{kl}
  template <typename U>
  auto Kronecker(const MMatrix<U> &rhs) const {
    using R                       = std::decay_t<decltype(std::declval<T>() * std::declval<U>())>;
    const std::size_t output_rows = CheckedElementCount(rows, rhs.rows);
    const std::size_t output_cols = CheckedElementCount(cols, rhs.cols);
    MMatrix<R>        output      = MMatrix<R>::Uninitialized(output_rows, output_cols);
    for (std::size_t first_row = 0; first_row < rows; ++first_row) {
      for (std::size_t first_col = 0; first_col < cols; ++first_col) {
        for (std::size_t second_row = 0; second_row < rhs.rows; ++second_row) {
          for (std::size_t second_col = 0; second_col < rhs.cols; ++second_col) {
            output(first_row * rhs.rows + second_row, first_col * rhs.cols + second_col) =
                data[first_row * cols + first_col] * rhs.data[second_row * rhs.cols + second_col];
          }
        }
      }
    }
    return output;
  }

  // Compute a Kronecker product with conceptual rows mapped to destinations
  // C_{p(i,k),(j,l)} = A_{ij} B_{kl}
  template <typename U, typename IndexContainer>
  auto Kronecker(const MMatrix<U> &rhs, const IndexContainer &destination_rows) const {
    using R                       = std::decay_t<decltype(std::declval<T>() * std::declval<U>())>;
    const std::size_t output_rows = CheckedElementCount(rows, rhs.rows);
    const std::size_t output_cols = CheckedElementCount(cols, rhs.cols);
    ValidateDestinationRows(destination_rows, output_rows, "MMatrix::Kronecker");
    MMatrix<R> output = MMatrix<R>::Uninitialized(output_rows, output_cols);
    for (std::size_t first_row = 0; first_row < rows; ++first_row) {
      for (std::size_t second_row = 0; second_row < rhs.rows; ++second_row) {
        const std::size_t source_row = first_row * rhs.rows + second_row;
        const std::size_t output_row = static_cast<std::size_t>(destination_rows[source_row]);
        for (std::size_t first_col = 0; first_col < cols; ++first_col) {
          for (std::size_t second_col = 0; second_col < rhs.cols; ++second_col) {
            output(output_row, first_col * rhs.cols + second_col) =
                data[first_row * cols + first_col] * rhs.data[second_row * rhs.cols + second_col];
          }
        }
      }
    }
    return output;
  }

  // Multiply a row-mapped Kronecker product by one matrix without forming it
  // C_{p(i,k),m} = sum_{jl} A_{ij} B_{kl} R_{(j,l),m}
  template <typename U, typename V, typename IndexContainer>
  auto KroneckerMultiply(const MMatrix<U> &rhs, const MMatrix<V> &right, const IndexContainer &destination_rows) const {
    using R                       = std::decay_t<decltype(std::declval<T>() * std::declval<U>() * std::declval<V>())>;
    const std::size_t output_rows = CheckedElementCount(rows, rhs.rows);
    const std::size_t tensor_cols = CheckedElementCount(cols, rhs.cols);
    if (right.rows != tensor_cols) {
      throw std::invalid_argument("MMatrix::KroneckerMultiply: matrix dimensions disagree");
    }
    ValidateDestinationRows(destination_rows, output_rows, "MMatrix::KroneckerMultiply");
    MMatrix<R> output(output_rows, right.cols, R{});
    if (right.cols == 0) { return output; }
    for (std::size_t first_row = 0; first_row < rows; ++first_row) {
      for (std::size_t second_row = 0; second_row < rhs.rows; ++second_row) {
        const std::size_t source_row = first_row * rhs.rows + second_row;
        R *output_data = output.data + static_cast<std::size_t>(destination_rows[source_row]) * output.cols;
        for (std::size_t first_col = 0; first_col < cols; ++first_col) {
          const T &first_value = data[first_row * cols + first_col];
          if (math::IsZero(first_value)) { continue; }
          for (std::size_t second_col = 0; second_col < rhs.cols; ++second_col) {
            const U &second_value = rhs.data[second_row * rhs.cols + second_col];
            if (math::IsZero(second_value)) { continue; }
            const auto coefficient = first_value * second_value;
            if (math::IsZero(coefficient)) { continue; }
            const V *right_data = right.data + (first_col * rhs.cols + second_col) * right.cols;
            for (std::size_t output_col = 0; output_col < right.cols; ++output_col) {
              const V &right_value = right_data[output_col];
              if (math::IsZero(right_value)) { continue; }
              output_data[output_col] += coefficient * right_value;
            }
          }
        }
      }
    }
    return output;
  }

  // Contract two factorized matrices with this matrix without forming them
  // C_{p(i,k),m} = a_i b_k sum_{jl} c_j d_l A_{(j,l),m}
  template <typename FirstRows, typename FirstColumns, typename SecondRows, typename SecondColumns,
            typename IndexContainer>
  auto FactorizedKroneckerMultiply(const FirstRows &first_rows, const FirstColumns &first_columns,
                                   const SecondRows &second_rows, const SecondColumns &second_columns,
                                   const IndexContainer &destination_rows) const {
    using FirstRow          = typename FirstRows::value_type;
    using FirstColumn       = typename FirstColumns::value_type;
    using SecondRow         = typename SecondRows::value_type;
    using SecondColumn      = typename SecondColumns::value_type;
    using ColumnCoefficient = std::decay_t<decltype(std::declval<FirstColumn>() * std::declval<SecondColumn>())>;
    using Contracted        = std::decay_t<decltype(std::declval<ColumnCoefficient>() * std::declval<T>())>;
    using RowCoefficient    = std::decay_t<decltype(std::declval<FirstRow>() * std::declval<SecondRow>())>;
    using R                 = std::decay_t<decltype(std::declval<RowCoefficient>() * std::declval<Contracted>())>;

    const std::size_t output_rows = CheckedElementCount(first_rows.size(), second_rows.size());
    const std::size_t tensor_rows = CheckedElementCount(first_columns.size(), second_columns.size());
    if (rows != tensor_rows) {
      throw std::invalid_argument("MMatrix::FactorizedKroneckerMultiply: matrix dimensions disagree");
    }
    ValidateDestinationRows(destination_rows, output_rows, "MMatrix::FactorizedKroneckerMultiply");
    if (cols == 0) { return MMatrix<R>(output_rows, 0); }

    std::vector<Contracted> contracted(cols, Contracted{});
    for (std::size_t first_col = 0; first_col < first_columns.size(); ++first_col) {
      for (std::size_t second_col = 0; second_col < second_columns.size(); ++second_col) {
        const ColumnCoefficient coefficient = first_columns[first_col] * second_columns[second_col];
        if (math::IsZero(coefficient)) { continue; }
        const T *matrix_row = data + (first_col * second_columns.size() + second_col) * cols;
        for (std::size_t col = 0; col < cols; ++col) {
          if (math::IsZero(matrix_row[col])) { continue; }
          contracted[col] += coefficient * matrix_row[col];
        }
      }
    }

    MMatrix<R> output = MMatrix<R>::Uninitialized(output_rows, cols);
    for (std::size_t first_row = 0; first_row < first_rows.size(); ++first_row) {
      for (std::size_t second_row = 0; second_row < second_rows.size(); ++second_row) {
        const std::size_t    source_row  = first_row * second_rows.size() + second_row;
        R                   *output_row  = output.data + static_cast<std::size_t>(destination_rows[source_row]) * cols;
        const RowCoefficient coefficient = first_rows[first_row] * second_rows[second_row];
        for (std::size_t col = 0; col < cols; ++col) { output_row[col] = coefficient * contracted[col]; }
      }
    }
    return output;
  }

  // Multiply a Kronecker product by one matrix without forming it
  // C_{(i,k),m} = sum_{jl} A_{ij} B_{kl} R_{(j,l),m}
  template <typename U, typename V>
  auto KroneckerMultiply(const MMatrix<U> &rhs, const MMatrix<V> &right) const {
    using R                       = std::decay_t<decltype(std::declval<T>() * std::declval<U>() * std::declval<V>())>;
    const std::size_t output_rows = CheckedElementCount(rows, rhs.rows);
    const std::size_t tensor_cols = CheckedElementCount(cols, rhs.cols);
    if (right.rows != tensor_cols) {
      throw std::invalid_argument("MMatrix::KroneckerMultiply: matrix dimensions disagree");
    }
    MMatrix<R> output(output_rows, right.cols, R{});
    if (right.cols == 0) { return output; }
    for (std::size_t first_row = 0; first_row < rows; ++first_row) {
      for (std::size_t second_row = 0; second_row < rhs.rows; ++second_row) {
        R *output_data = output.data + (first_row * rhs.rows + second_row) * output.cols;
        for (std::size_t first_col = 0; first_col < cols; ++first_col) {
          const T &first_value = data[first_row * cols + first_col];
          if (math::IsZero(first_value)) { continue; }
          for (std::size_t second_col = 0; second_col < rhs.cols; ++second_col) {
            const U &second_value = rhs.data[second_row * rhs.cols + second_col];
            if (math::IsZero(second_value)) { continue; }
            const auto coefficient = first_value * second_value;
            if (math::IsZero(coefficient)) { continue; }
            const V *right_data = right.data + (first_col * rhs.cols + second_col) * right.cols;
            for (std::size_t output_col = 0; output_col < right.cols; ++output_col) {
              const V &right_value = right_data[output_col];
              if (math::IsZero(right_value)) { continue; }
              output_data[output_col] += coefficient * right_value;
            }
          }
        }
      }
    }
    return output;
  }

  // ------------------------------------------------------------------

  // Sum all elements
  // output = sum_{ij} A_{ij}
  T Sum() const {
    T sum(T(0));
    for (std::size_t i = 0; i < rows * cols; ++i) { sum += data[i]; }
    return sum;
  }

  // Compute the nonconjugating sum of elementwise products with rhs
  // output = sum_{ij} A_{ij} B_{ij}
  template <typename U>
  auto ElementwiseProductSum(const MMatrix<U> &rhs) const -> decltype(std::declval<T>() * std::declval<U>()) {
    using Result = decltype(std::declval<T>() * std::declval<U>());
    if (rows != rhs.rows || cols != rhs.cols) {
      throw std::invalid_argument("MMatrix::ElementwiseProductSum: dimensions disagree");
    }
    Result result{};
    for (std::size_t index = 0; index < rows * cols; ++index) { result += data[index] * rhs.data[index]; }
    return result;
  }

  bool isEmpty() const {
    if ((size_row() == 0) || (size_col() == 0)) {
      return true;
    } else {
      return false;
    }
  }

  // Compute the Frobenius norm with scaled accumulation
  // ||A||_F = sqrt(sum_{ij} |A_{ij}|^2)
  double FrobNorm() const {
    long double scale = 0.0L;
    long double sum   = 1.0L;
    for (std::size_t index = 0; index < rows * cols; ++index) {
      const long double magnitude = [&]() {
        if constexpr (std::is_arithmetic_v<T>) {
          return std::abs(static_cast<long double>(data[index]));
        } else {
          return static_cast<long double>(std::abs(data[index]));
        }
      }();
      if (!std::isfinite(magnitude)) { return static_cast<double>(magnitude); }
      if (math::IsZero(magnitude)) { continue; }
      if (scale < magnitude) {
        const long double ratio = scale / magnitude;
        sum                     = 1.0L + sum * ratio * ratio;
        scale                   = magnitude;
      } else {
        const long double ratio = magnitude / scale;
        sum += ratio * ratio;
      }
    }
    return math::IsZero(scale) ? 0.0 : static_cast<double>(scale * std::sqrt(sum));
  }

  // Frobenius norm squared
  // ||A||_F^2 = sum_{ij} |A_{ij}|^2
  double FrobNorm2() const {
    double sum = 0.0;
    for (std::size_t i = 0; i < rows * cols; ++i) { sum += Abs2(data[i]); }
    return sum;
  }

  // Compute the squared Frobenius norm after right diagonal multiplication
  // ||A D||_F^2 = sum_{ij} |A_{ij} d_j|^2
  template <typename Diagonal>
  double RightDiagonalFrobNorm2(const Diagonal &diagonal) const {
    if (!std::cmp_equal(diagonal.size(), cols)) {
      throw std::invalid_argument("MMatrix::RightDiagonalFrobNorm2: dimensions disagree");
    }
    double sum = 0.0;
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) { sum += Abs2(data[row * cols + col] * diagonal[col]); }
    }
    return sum;
  }

  // Compute the squared norm of entries selected by a boolean mask
  // output = sum_{ij: mask_{ij}} |A_{ij}|^2
  double MaskedSquaredNorm(const MMatrix<bool> &mask) const {
    if (rows != mask.size_row() || cols != mask.size_col()) {
      throw std::invalid_argument("MMatrix::MaskedSquaredNorm: mask dimensions disagree");
    }
    double result = 0.0;
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        if (mask(row, col)) { result += Abs2(data[row * cols + col]); }
      }
    }
    return result;
  }

  // Scale entries selected by a boolean mask
  template <typename Scalar>
  void ScaleMasked(const MMatrix<bool> &mask, const Scalar &scale) {
    if (rows != mask.size_row() || cols != mask.size_col()) {
      throw std::invalid_argument("MMatrix::ScaleMasked: mask dimensions disagree");
    }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        if (mask(row, col)) { data[row * cols + col] *= scale; }
      }
    }
  }

  // Compute true when every real or complex matrix entry is finite
  bool IsFinite() const {
    for (std::size_t i = 0; i < rows * cols; ++i) {
      if constexpr (std::is_arithmetic_v<T>) {
        if (!std::isfinite(data[i])) { return false; }
      } else if constexpr (requires(const T &value) {
                             value.real();
                             value.imag();
                           }) {
        if (!std::isfinite(data[i].real()) || !std::isfinite(data[i].imag())) { return false; }
      } else {
        static_assert(std::is_arithmetic_v<T>, "MMatrix::IsFinite requires real or complex entries");
      }
    }
    return true;
  }

  // Compute whether this is an exactly diagonal square matrix
  bool IsDiagonal() const {
    if (rows != cols) { return false; }
    for (std::size_t row = 0; row < rows; ++row) {
      const T *matrix_row = data + cols * row;
      for (std::size_t col = 0; col < row; ++col) {
        if (!IsZero(matrix_row[col])) { return false; }
      }
      for (std::size_t col = row + 1; col < cols; ++col) {
        if (!IsZero(matrix_row[col])) { return false; }
      }
    }
    return true;
  }

  // Compute true when another matrix agrees entrywise within a tolerance
  // max_{ij} |A_{ij} - B_{ij}| <= tolerance
  bool IsApprox(const MMatrix &other, const double tolerance = 1.0e-9) const {
    if (rows != other.rows || cols != other.cols || !std::isfinite(tolerance) || tolerance < 0.0) { return false; }
    for (std::size_t index = 0; index < rows * cols; ++index) {
      const double difference = [&]() {
        if constexpr (std::is_arithmetic_v<T>) {
          return static_cast<double>(
              std::abs(static_cast<long double>(data[index]) - static_cast<long double>(other.data[index])));
        } else {
          return static_cast<double>(std::abs(data[index] - other.data[index]));
        }
      }();
      if (!std::isfinite(difference) || difference > tolerance) { return false; }
    }
    return true;
  }

  // Compute the squared Euclidean norm of one matrix column
  // output = sum_i |A_{ij}|^2
  double ColumnNorm2(const std::size_t column) const {
    if (column >= cols) { throw std::out_of_range("MMatrix::ColumnNorm2: column outside matrix"); }
    double sum = 0.0;
    for (std::size_t row = 0; row < rows; ++row) { sum += Abs2(data[row * cols + column]); }
    return sum;
  }

  // Compute one matrix column as a contiguous vector
  std::vector<T> Column(const std::size_t column) const {
    if (column >= cols) { throw std::out_of_range("MMatrix::Column: column outside matrix"); }
    std::vector<T> output(rows);
    for (std::size_t row = 0; row < rows; ++row) { output[row] = data[row * cols + column]; }
    return output;
  }

  // Compute selected matrix columns in the requested order
  MMatrix SelectColumns(const std::vector<std::size_t> &indices) const {
    MMatrix output = Uninitialized(rows, indices.size());
    for (std::size_t output_col = 0; output_col < indices.size(); ++output_col) {
      const std::size_t source_col = indices[output_col];
      if (source_col >= cols) { throw std::out_of_range("MMatrix::SelectColumns: column outside matrix"); }
      for (std::size_t row = 0; row < rows; ++row) {
        output.data[row * output.cols + output_col] = data[row * cols + source_col];
      }
    }
    return output;
  }

  // Transform one tensor factor of the column space without forming a Kronecker matrix
  MMatrix TransformColumnAxis(const MMatrix &op, std::size_t before, std::size_t after) const {
    if (before == 0 || after == 0 || op.cols == 0 || cols % after != 0 ||
        cols / after % op.cols != 0 || cols / after / op.cols != before) {
      throw std::invalid_argument("MMatrix::TransformColumnAxis: dimensions disagree");
    }
    if (op.rows > std::numeric_limits<std::size_t>::max() / before / after) {
      throw std::length_error("MMatrix::TransformColumnAxis: output dimensions overflow");
    }
    MMatrix output(rows, before * op.rows * after, T(0));
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t prefix = 0; prefix < before; ++prefix) {
        for (std::size_t out = 0; out < op.rows; ++out) {
          for (std::size_t in = 0; in < op.cols; ++in) {
            for (std::size_t suffix = 0; suffix < after; ++suffix) {
              output[row][(prefix * op.rows + out) * after + suffix] +=
                  op[out][in] * data[row * cols + (prefix * op.cols + in) * after + suffix];
            }
          }
        }
      }
    }
    return output;
  }

  // Compute selected matrix rows in the requested order
  MMatrix SelectRows(const std::vector<std::size_t> &indices) const {
    MMatrix output = Uninitialized(indices.size(), cols);
    for (std::size_t output_row = 0; output_row < indices.size(); ++output_row) {
      const std::size_t source_row = indices[output_row];
      if (source_row >= rows) { throw std::out_of_range("MMatrix::SelectRows: row outside matrix"); }
      std::copy_n(data + source_row * cols, cols, output.data + output_row * cols);
    }
    return output;
  }

  // Replace selected matrix rows from a compact source matrix
  void SetRows(const std::vector<std::size_t> &indices, const MMatrix &source) {
    if (source.rows != indices.size() || source.cols != cols) {
      throw std::invalid_argument("MMatrix::SetRows: dimensions disagree");
    }
    std::vector<bool> assigned(rows, false);
    for (const std::size_t target_row : indices) {
      if (target_row >= rows) { throw std::out_of_range("MMatrix::SetRows: row outside matrix"); }
      if (assigned[target_row]) { throw std::invalid_argument("MMatrix::SetRows: duplicate row"); }
      assigned[target_row] = true;
    }
    if (cols == 0) { return; }
    if (this == &source) {
      const MMatrix source_copy(source);
      SetRows(indices, source_copy);
      return;
    }
    for (std::size_t source_row = 0; source_row < indices.size(); ++source_row) {
      const std::size_t target_row = indices[source_row];
      std::copy_n(source.data + source_row * source.cols, cols, data + target_row * cols);
    }
  }

  // Scale one matrix column in place
  // A_{ij} <- c A_{ij} at fixed j
  template <typename Scalar>
  void ScaleColumn(const std::size_t column, const Scalar &scale) {
    if (column >= cols) { throw std::out_of_range("MMatrix::ScaleColumn: column outside matrix"); }
    for (std::size_t row = 0; row < rows; ++row) { data[row * cols + column] *= scale; }
  }

  // Scale all matrix columns by one diagonal vector in place
  // A_{ij} <- A_{ij} d_j
  template <typename Diagonal>
  void ScaleColumns(const Diagonal &diagonal) {
    if (!std::cmp_equal(diagonal.size(), cols)) {
      throw std::invalid_argument("MMatrix::ScaleColumns: dimensions disagree");
    }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) { data[row * cols + col] *= diagonal[col]; }
    }
  }

  // Scale one matrix row in place
  // A_{ij} <- c A_{ij} at fixed i
  template <typename Scalar>
  void ScaleRow(const std::size_t row, const Scalar &scale) {
    if (row >= rows) { throw std::out_of_range("MMatrix::ScaleRow: row outside matrix"); }
    T *row_data = data + row * cols;
    for (std::size_t col = 0; col < cols; ++col) { row_data[col] *= scale; }
  }

  // Compute this matrix multiplied on the left by one diagonal vector
  // output_{ij} = d_i A_{ij}
  template <typename U>
  auto LeftDiagonalProduct(const std::vector<U> &diagonal) const
      -> MMatrix<decltype(std::declval<U>() * std::declval<T>())> {
    using R = decltype(std::declval<U>() * std::declval<T>());
    if (diagonal.size() != rows) { throw std::invalid_argument("MMatrix::LeftDiagonalProduct: dimensions disagree"); }
    MMatrix<R> out = MMatrix<R>::Uninitialized(rows, cols);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        out.data[row * cols + col] = diagonal[row] * data[row * cols + col];
      }
    }
    return out;
  }

  // Compute this matrix multiplied on the right by one diagonal vector
  // output_{ij} = A_{ij} d_j
  template <typename U>
  auto RightDiagonalProduct(const std::vector<U> &diagonal) const
      -> MMatrix<decltype(std::declval<T>() * std::declval<U>())> {
    using R = decltype(std::declval<T>() * std::declval<U>());
    if (diagonal.size() != cols) { throw std::invalid_argument("MMatrix::RightDiagonalProduct: dimensions disagree"); }
    MMatrix<R> out = MMatrix<R>::Uninitialized(rows, cols);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        out.data[row * cols + col] = data[row * cols + col] * diagonal[col];
      }
    }
    return out;
  }

  // Compute true when this square matrix is Hermitian within a tolerance
  // max_{ij}|A_{ij} - A_{ji}^*| <= tolerance
  bool IsHermitian(const double tolerance) const {
    if (rows != cols || !std::isfinite(tolerance) || tolerance < 0.0) { return false; }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = row; col < cols; ++col) {
        const auto difference = std::abs(data[row * cols + col] - Conjugate(data[col * cols + row]));
        if (!std::isfinite(difference) || difference > tolerance) { return false; }
      }
    }
    return true;
  }

  // Apply one scalar function independently to every matrix entry
  template <typename Function>
  auto Transform(Function function) const -> MMatrix<std::decay_t<decltype(function(std::declval<const T &>()))>> {
    using R        = std::decay_t<decltype(function(std::declval<const T &>()))>;
    MMatrix<R> out = MMatrix<R>::Uninitialized(rows, cols);
    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = function(data[i]); }
    return out;
  }

  // Compute the induced matrix one norm
  // ||A||_1 = max_j sum_i |A_{ij}|
  double NormOne() const {
    double norm = 0.0;
    for (std::size_t j = 0; j < cols; ++j) {
      double column = 0.0;
      for (std::size_t i = 0; i < rows; ++i) { column += std::abs((*this)(i, j)); }
      norm = std::max(norm, column);
    }
    return norm;
  }

  // Convert this matrix to dynamically sized Eigen storage
  auto ToEigen() const {
    Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic> output(static_cast<Eigen::Index>(rows),
                                                            static_cast<Eigen::Index>(cols));
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        output(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = data[row * cols + col];
      }
    }
    return output;
  }

  // Compute the singular values in descending order and optional right singular vectors
  // A = U diag(sigma_i) V^dagger
  std::vector<double> SingularValues(MMatrix *right_vectors = nullptr) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    if (rows == 0 || cols == 0) {
      if (right_vectors != nullptr) { *right_vectors = MMatrix(); }
      return {};
    }
    if (!IsFinite()) { throw std::invalid_argument("MMatrix::SingularValues: matrix is not finite"); }
    using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const unsigned int options = right_vectors != nullptr ? static_cast<unsigned int>(Eigen::ComputeThinV) : 0U;
    const Eigen::JacobiSVD<EigenMatrix> decomposition(ToEigen(), options);
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::SingularValues: decomposition failed");
    }
    const auto          singular = decomposition.singularValues();
    std::vector<double> output(static_cast<std::size_t>(singular.size()));
    for (Eigen::Index index = 0; index < singular.size(); ++index) {
      output[static_cast<std::size_t>(index)] = static_cast<double>(singular(index));
    }
    if (right_vectors != nullptr) { *right_vectors = FromEigen(decomposition.matrixV()); }
    return output;
  }

  // Compute the largest singular value of this matrix
  // ||A||_2 = max_i sigma_i
  double MaxSingularValue() const requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    const std::vector<double> singular = SingularValues();
    return singular.empty() ? 0.0 : singular.front();
  }

  // Solve one rectangular least-squares system with pivoted QR
  // X = argmin_X ||AX-B||_F
  std::vector<T> SolveLeastSquares(const std::vector<T> &rhs) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    if (rows == 0 || cols == 0 || rhs.size() != rows || !IsFinite()) {
      throw std::invalid_argument("MMatrix::SolveLeastSquares: invalid matrix or source dimensions");
    }
    for (const T &value : rhs) {
      if constexpr (std::is_same_v<T, double>) {
        if (!std::isfinite(value)) { throw std::invalid_argument("MMatrix::SolveLeastSquares: source is not finite"); }
      } else if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
        throw std::invalid_argument("MMatrix::SolveLeastSquares: source is not finite");
      }
    }

    using EigenMatrix                                    = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    using EigenVector                                    = Eigen::Matrix<T, Eigen::Dynamic, 1>;
    const EigenMatrix                             matrix = ToEigen();
    const Eigen::Map<const EigenVector>           source(rhs.data(), static_cast<Eigen::Index>(rhs.size()));
    const Eigen::ColPivHouseholderQR<EigenMatrix> decomposition(matrix);
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::SolveLeastSquares: decomposition failed");
    }
    const EigenVector solution = decomposition.solve(source);
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::SolveLeastSquares: solve failed");
    }
    std::vector<T> output(static_cast<std::size_t>(solution.size()));
    for (Eigen::Index index = 0; index < solution.size(); ++index) {
      const T value = solution(index);
      if constexpr (std::is_same_v<T, double>) {
        if (!std::isfinite(value)) { throw std::runtime_error("MMatrix::SolveLeastSquares: solution is not finite"); }
      } else if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
        throw std::runtime_error("MMatrix::SolveLeastSquares: solution is not finite");
      }
      output[static_cast<std::size_t>(index)] = value;
    }
    return output;
  }

  // Compute the Moore-Penrose inverse with an explicit relative singular cut
  // A^+ = V diag(1/sigma_i) U^dagger for retained sigma_i
  MMatrix PseudoInverse(const double relative_singular_cut, PseudoInverseDiagnostics *diagnostics = nullptr) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    if (rows == 0 || cols == 0 || !IsFinite() || !std::isfinite(relative_singular_cut) || relative_singular_cut < 0.0) {
      throw std::invalid_argument("MMatrix::PseudoInverse: invalid matrix or singular cut");
    }
    if (diagnostics != nullptr) {
      *diagnostics = {};
    }

    using EigenMatrix                          = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const EigenMatrix                   matrix = ToEigen();
    const Eigen::JacobiSVD<EigenMatrix> decomposition(matrix, Eigen::ComputeThinU | Eigen::ComputeThinV);
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::PseudoInverse: decomposition failed");
    }
    const auto singular = decomposition.singularValues();
    if (singular.size() == 0) { throw std::runtime_error("MMatrix::PseudoInverse: decomposition is empty"); }
    if (!(singular(0) > 0.0)) {
      if (diagnostics != nullptr) { diagnostics->condition_number = std::numeric_limits<double>::infinity(); }
      return MMatrix(cols, rows, T{});
    }

    const double numerical_cut =
        std::numeric_limits<double>::epsilon() * static_cast<double>(std::max(rows, cols)) * singular(0);
    const double applied_cut      = std::max(numerical_cut, relative_singular_cut * singular(0));
    auto         inverse_singular = singular;
    inverse_singular.setZero();
    std::size_t numerical_rank = 0;
    std::size_t retained_rank  = 0;
    for (Eigen::Index index = 0; index < singular.size(); ++index) {
      if (singular(index) > numerical_cut) { ++numerical_rank; }
      if (singular(index) > applied_cut) {
        inverse_singular(index) = 1.0 / singular(index);
        ++retained_rank;
      }
    }
    if (diagnostics != nullptr) {
      diagnostics->numerical_rank   = numerical_rank;
      diagnostics->retained_rank    = retained_rank;
      diagnostics->condition_number = singular(0) / singular(singular.size() - 1);
    }
    return FromEigen(decomposition.matrixV() * inverse_singular.asDiagonal() * decomposition.matrixU().adjoint());
  }

  // Diagonalize one finite self-adjoint matrix in ascending eigenvalue order
  // A V = V diag(lambda_i)
  std::vector<double> SelfAdjointEigenvalues(const double relative_tolerance, MMatrix *eigenvectors = nullptr) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    if (rows == 0 || rows != cols || !IsFinite() || !std::isfinite(relative_tolerance) || relative_tolerance < 0.0) {
      throw std::invalid_argument("MMatrix::SelfAdjointEigenvalues: invalid matrix or tolerance");
    }
    const MMatrix adjoint         = Dagger();
    const double  norm            = FrobNorm();
    const double  difference_norm = ((*this) - adjoint).FrobNorm();
    const double  residual        = math::IsZero(norm) ? difference_norm : difference_norm / norm;
    if (!std::isfinite(residual) || residual > relative_tolerance) {
      throw std::invalid_argument("MMatrix::SelfAdjointEigenvalues: matrix is not self-adjoint");
    }
    MMatrix self_adjoint = (*this) * T(0.5);
    self_adjoint.AddScaled(adjoint, T(0.5));
    using EigenMatrix = Eigen::Matrix<T, Eigen::Dynamic, Eigen::Dynamic>;
    const Eigen::SelfAdjointEigenSolver<EigenMatrix> decomposition(self_adjoint.ToEigen());
    if (decomposition.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::SelfAdjointEigenvalues: decomposition failed");
    }
    const auto          values = decomposition.eigenvalues();
    std::vector<double> output(static_cast<std::size_t>(values.size()));
    for (Eigen::Index index = 0; index < values.size(); ++index) {
      output[static_cast<std::size_t>(index)] = static_cast<double>(values(index));
    }
    if (eigenvectors != nullptr) { *eigenvectors = FromEigen(decomposition.eigenvectors()); }
    return output;
  }

  // Compute the entropy minus Tr(A log A) of a positive self-adjoint matrix
  double SelfAdjointEntropy(const double relative_tolerance) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    double entropy = 0.0;
    for (const double value : SelfAdjointEigenvalues(relative_tolerance)) {
      if (value < -relative_tolerance) {
        throw std::domain_error("MMatrix::SelfAdjointEntropy: matrix is not positive semidefinite");
      }
      if (value > relative_tolerance) { entropy -= value * std::log(value); }
    }
    return entropy;
  }

  // Decompose a positive self-adjoint matrix into weighted pure vectors
  std::vector<std::vector<T>> PositiveSpectralVectors(const double relative_tolerance) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    MMatrix                     eigenvectors;
    const auto                  eigenvalues = SelfAdjointEigenvalues(relative_tolerance, &eigenvectors);
    std::vector<std::vector<T>> states;
    for (std::size_t index = 0; index < eigenvalues.size(); ++index) {
      if (eigenvalues[index] < -relative_tolerance) {
        throw std::domain_error(
            "MMatrix::PositiveSpectralVectors: matrix is "
            "not positive semidefinite");
      }
      if (eigenvalues[index] <= relative_tolerance) { continue; }
      auto         state  = eigenvectors.Column(index);
      const double weight = std::sqrt(eigenvalues[index]);
      for (auto &value : state) { value *= weight; }
      states.push_back(std::move(state));
    }
    return states;
  }

  // Compute the principal square root of one positive semidefinite matrix
  // A^(1/2) = V diag(sqrt(lambda_i)) V^dagger
  MMatrix PrincipalPositiveSemidefiniteSquareRoot(
      const double relative_tolerance = 64.0 * std::numeric_limits<double>::epsilon()) const
      requires(std::is_same_v<T, double> || std::is_same_v<T, std::complex<double>>) {
    if (rows == 0 || rows != cols || !IsFinite() || !std::isfinite(relative_tolerance) || relative_tolerance < 0.0) {
      throw std::invalid_argument("MMatrix::PrincipalPositiveSemidefiniteSquareRoot: invalid input");
    }
    const double norm = FrobNorm();
    if (math::IsZero(norm)) { return MMatrix(rows, cols, T{}); }
    const double  dimension            = static_cast<double>(std::max<std::size_t>(1, rows));
    const double  eigenvalue_tolerance = relative_tolerance * dimension * norm;
    const MMatrix adjoint              = Dagger();
    if (((*this) - adjoint).FrobNorm() > eigenvalue_tolerance) {
      throw std::invalid_argument(
          "MMatrix::PrincipalPositiveSemidefiniteSquareRoot: matrix is not "
          "self-adjoint");
    }
    MMatrix target = (*this) * T(0.5);
    target.AddScaled(adjoint, T(0.5));
    MMatrix                   eigenvectors;
    const std::vector<double> eigenvalues =
        target.SelfAdjointEigenvalues(relative_tolerance * dimension, &eigenvectors);
    MMatrix diagonal(rows, cols, T{});
    for (std::size_t index = 0; index < rows; ++index) {
      if (eigenvalues[index] < -eigenvalue_tolerance) {
        throw std::domain_error(
            "MMatrix::PrincipalPositiveSemidefiniteSquareRoot: matrix is not "
            "positive semidefinite");
      }
      // Resolve numerical null eigenvalues symmetrically before taking their square roots
      diagonal(index, index) = eigenvalues[index] > eigenvalue_tolerance ? T(std::sqrt(eigenvalues[index])) : T{};
    }
    const MMatrix root                     = eigenvectors * diagonal * eigenvectors.Dagger();
    const double  reconstruction_tolerance = 4.0 * relative_tolerance * dimension * norm;
    const double  reconstruction_error     = (root.Dagger() * root - target).FrobNorm();
    if (!root.IsFinite() || !std::isfinite(reconstruction_error) || reconstruction_error > reconstruction_tolerance) {
      throw std::runtime_error(
          "MMatrix::PrincipalPositiveSemidefiniteSquareRoot: "
          "reconstruction failed");
    }
    return root;
  }

  // Compute a natural-order lower Cholesky factor for one real matrix
  // A = L L^T
  // The pivot tolerance is relative above unit scale and absolute below it
  MMatrix CholeskyLower(const T pivot_tolerance, const CholeskyMode mode) const requires(std::is_floating_point_v<T>) {
    if (rows == 0 || rows != cols || !IsFinite() || !std::isfinite(pivot_tolerance) || pivot_tolerance < T(0) ||
        (mode != CholeskyMode::StrictPositiveDefinite && mode != CholeskyMode::AllowSemidefinite)) {
      throw std::invalid_argument("MMatrix::CholeskyLower: invalid matrix or tolerance");
    }
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = row + 1; col < cols; ++col) {
        const T tolerance =
            pivot_tolerance * std::max({T(1), std::abs((*this)(row, col)), std::abs((*this)(col, row))});
        if (std::abs((*this)(row, col) - (*this)(col, row)) > tolerance) {
          throw std::invalid_argument("MMatrix::CholeskyLower: matrix is not symmetric");
        }
      }
    }

    MMatrix lower(rows, cols, T(0));
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col <= row; ++col) {
        T value = (*this)(row, col) / T(2) + (*this)(col, row) / T(2);
        for (std::size_t inner = 0; inner < col; ++inner) { value -= lower(row, inner) * lower(col, inner); }
        if (!std::isfinite(value)) { throw std::domain_error("MMatrix::CholeskyLower: pivot is not finite"); }
        if (row == col) {
          const T diagonal_tolerance = pivot_tolerance * std::max(T(1), std::abs((*this)(row, row)));
          if (mode == CholeskyMode::StrictPositiveDefinite) {
            if (!(value > diagonal_tolerance)) {
              throw std::domain_error("MMatrix::CholeskyLower: matrix is not positive definite");
            }
            lower(row, col) = std::sqrt(value);
          } else {
            if (value < -diagonal_tolerance) {
              throw std::domain_error(
                  "MMatrix::CholeskyLower: matrix is not positive "
                  "semidefinite");
            }
            lower(row, col) = std::sqrt(std::max(T(0), value));
          }
        } else if (mode == CholeskyMode::StrictPositiveDefinite || lower(col, col) > pivot_tolerance) {
          lower(row, col) = value / lower(col, col);
        } else {
          const T residual_tolerance =
              pivot_tolerance * std::max({T(1), std::abs((*this)(row, col)), std::abs((*this)(col, row))});
          if (std::abs(value) > residual_tolerance) {
            throw std::domain_error(
                "MMatrix::CholeskyLower: matrix is not positive "
                "semidefinite");
          }
        }
      }
    }

    T maximum_entry = T(0);
    T maximum_error = T(0);
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        const T target        = (*this)(row, col) / T(2) + (*this)(col, row) / T(2);
        T       reconstructed = T(0);
        for (std::size_t inner = 0; inner <= std::min(row, col); ++inner) {
          reconstructed += lower(row, inner) * lower(col, inner);
        }
        maximum_entry = std::max(maximum_entry, std::abs(target));
        maximum_error = std::max(maximum_error, std::abs(reconstructed - target));
      }
    }
    const T verification_relative = std::max(pivot_tolerance, T(8) * std::numeric_limits<T>::epsilon());
    const T verification_tolerance =
        T(8) * static_cast<T>(std::max<std::size_t>(1, rows)) * verification_relative * std::max(T(1), maximum_entry);
    if (!std::isfinite(maximum_error) || maximum_error > verification_tolerance) {
      throw std::runtime_error("MMatrix::CholeskyLower: reconstruction failed");
    }
    return lower;
  }

  // Raise this complex matrix to a real power on the principal branch
  // A^p = V diag(exp[p Log(lambda_i)]) V^-1
  template <typename U = T>
  requires(!std::is_arithmetic_v<U> && std::is_same_v<U, T>) MMatrix
      PrincipalPower(const typename U::value_type exponent)
  const {
    using Real        = typename U::value_type;
    using EigenMatrix = Eigen::Matrix<U, Eigen::Dynamic, Eigen::Dynamic>;
    if (rows == 0 || rows != cols || !std::isfinite(exponent)) {
      throw std::invalid_argument("MMatrix::PrincipalPower: invalid matrix power input");
    }

    constexpr double branch_tolerance = 1.0e-12;
    const auto       principal_power  = [exponent](const U value) {
      const Real scale = std::max(Real(1), std::abs(value));
      if (std::abs(value) <= Real(branch_tolerance) * scale) {
        throw std::runtime_error("MMatrix::PrincipalPower: matrix has an eigenvalue at zero");
      }
      if (value.real() <= Real(0) && std::abs(value.imag()) <= Real(branch_tolerance) * scale) {
        throw std::runtime_error("MMatrix::PrincipalPower: principal branch cut crossed");
      }
      return std::exp(exponent * std::log(value));
    };
    if (IsDiagonal() && IsFinite()) {
      MMatrix output(rows, cols, U{});
      for (std::size_t index = 0; index < rows; ++index) {
        output.data[index * cols + index] = principal_power(data[index * cols + index]);
      }
      return output;
    }

    const EigenMatrix                            matrix = ToEigen();
    const Eigen::ComplexEigenSolver<EigenMatrix> solver(matrix);
    if (solver.info() != Eigen::Success) {
      throw std::runtime_error("MMatrix::PrincipalPower: eigensystem construction failed");
    }

    const EigenMatrix                   eigenvectors = solver.eigenvectors();
    const Eigen::JacobiSVD<EigenMatrix> svd(eigenvectors);
    const auto                          singular = svd.singularValues();
    if (singular.size() == 0 || !(singular(0) > Real(0)) || !(singular(singular.size() - 1) > Real(0))) {
      throw std::runtime_error("MMatrix::PrincipalPower: singular eigenvector matrix");
    }

    constexpr double condition_max = 1.0e12;
    const double     condition     = static_cast<double>(singular(0) / singular(singular.size() - 1));
    if (!std::isfinite(condition) || condition > condition_max) {
      throw std::runtime_error("MMatrix::PrincipalPower: ill-conditioned eigenvector matrix");
    }

    const Eigen::FullPivLU<EigenMatrix> decomposition(eigenvectors);
    if (!decomposition.isInvertible()) {
      throw std::runtime_error("MMatrix::PrincipalPower: defective matrix power input");
    }

    EigenMatrix value_diag = EigenMatrix::Zero(matrix.rows(), matrix.cols());
    EigenMatrix power_diag = EigenMatrix::Zero(matrix.rows(), matrix.cols());
    for (Eigen::Index index = 0; index < matrix.rows(); ++index) {
      const U value            = solver.eigenvalues()(index);
      value_diag(index, index) = value;
      power_diag(index, index) = principal_power(value);
    }

    const EigenMatrix inverse       = decomposition.inverse();
    const EigenMatrix reconstructed = eigenvectors * value_diag * inverse;
    const Real        residual      = (reconstructed - matrix).norm() / std::max(Real(1), matrix.norm());
    if (!std::isfinite(residual) || residual > Real(1.0e-9)) {
      throw std::runtime_error("MMatrix::PrincipalPower: eigensystem reconstruction failed");
    }
    return FromEigen(eigenvectors * power_diag * inverse);
  }

  // Compute a square identity matrix
  static MMatrix IdentityMatrix(const std::size_t n) { return MMatrix(n, n, "eye"); }

  // Compute a square matrix with entries from one diagonal vector
  static MMatrix DiagonalMatrix(const std::vector<T> &diagonal) {
    MMatrix output(diagonal.size(), diagonal.size(), T{});
    for (std::size_t index = 0; index < diagonal.size(); ++index) {
      output.data[index * output.cols + index] = diagonal[index];
    }
    return output;
  }

  // Build a row major matrix from equally sized column containers
  template <typename ColumnContainer>
  static MMatrix FromColumns(const ColumnContainer &columns) {
    if (columns.empty()) { return MMatrix(); }
    const std::size_t column_rows = columns.front().size();
    for (const auto &column : columns) {
      if (column.size() != column_rows) {
        throw std::invalid_argument("MMatrix::FromColumns: column dimensions disagree");
      }
    }
    MMatrix output = Uninitialized(column_rows, columns.size());
    for (std::size_t col = 0; col < columns.size(); ++col) {
      for (std::size_t row = 0; row < column_rows; ++row) { output.data[row * output.cols + col] = columns[col][row]; }
    }
    return output;
  }

  // Convert one finite Eigen matrix expression to this matrix type
  template <typename Derived>
  static MMatrix FromEigen(const Eigen::MatrixBase<Derived> &input) {
    static_assert(std::is_same_v<typename Derived::Scalar, T>, "MMatrix::FromEigen requires identical scalar types");
    MMatrix output = Uninitialized(static_cast<std::size_t>(input.rows()), static_cast<std::size_t>(input.cols()));
    for (Eigen::Index row = 0; row < input.rows(); ++row) {
      for (Eigen::Index col = 0; col < input.cols(); ++col) {
        output.data[static_cast<std::size_t>(row) * output.cols + static_cast<std::size_t>(col)] = input(row, col);
      }
    }
    if (!output.IsFinite()) { throw std::runtime_error("MMatrix::FromEigen: non-finite matrix element"); }
    return output;
  }

  // Construct an arbitrary real orthogonal mixing matrix
  // O = product_{i<j} R_{ij}(theta_{ij}), O O^T = I
  static MMatrix MixingReal(const std::vector<double> &theta,
                            const std::size_t          dimension) requires(std::is_same_v<T, double>) {
    if (dimension == 0) { throw std::invalid_argument("MMatrix::MixingReal: matrix dimension must be positive"); }
    if (dimension > 1 && dimension > std::numeric_limits<std::size_t>::max() / (dimension - 1)) {
      throw std::length_error("MMatrix::MixingReal: dimension overflow");
    }
    const std::size_t expected = dimension * (dimension - 1) / 2;
    if (theta.size() != expected) {
      throw std::invalid_argument("MMatrix::MixingReal: expected " + std::to_string(expected) + " angles for N = " +
                                  std::to_string(dimension) + ", received " + std::to_string(theta.size()));
    }

    MMatrix     mixing      = IdentityMatrix(dimension);
    std::size_t angle_index = 0;
    for (std::size_t first = 0; first + 1 < dimension; ++first) {
      for (std::size_t second = first + 1; second < dimension; ++second) {
        const double angle = theta[angle_index++];
        if (!std::isfinite(angle)) { throw std::invalid_argument("MMatrix::MixingReal: mixing angles must be finite"); }
        MMatrix      rotation    = IdentityMatrix(dimension);
        const double cosine      = std::cos(angle);
        const double sine        = std::sin(angle);
        rotation[first][first]   = cosine;
        rotation[first][second]  = sine;
        rotation[second][first]  = -sine;
        rotation[second][second] = cosine;
        mixing                   = rotation * mixing;
      }
    }
    return mixing;
  }

  // Trace over one factor of A_(i a),(j b) in the ordered n1 x n2 product basis
  MMatrix PartialTrace(std::size_t n1, std::size_t n2, TensorFactor traced) const {
    ValidateTensorShape(n1, n2);
    const bool first = traced == TensorFactor::First;
    const std::size_t kept = first ? n2 : n1, summed = first ? n1 : n2;
    MMatrix out(kept, kept, T{});
    for (std::size_t i = 0; i < kept; ++i) {
      for (std::size_t j = 0; j < kept; ++j) {
        for (std::size_t k = 0; k < summed; ++k) {
          out[i][j] += first ? (*this)[k * n2 + i][k * n2 + j] : (*this)[i * n2 + k][j * n2 + k];
        }
      }
    }
    return out;
  }

  // Transpose one factor of A_(i a),(j b) without complex conjugation
  MMatrix PartialTranspose(std::size_t n1, std::size_t n2, TensorFactor transposed) const {
    ValidateTensorShape(n1, n2);
    MMatrix out(rows, cols, UninitializedTag{});
    for (std::size_t i = 0; i < n1; ++i) {
      for (std::size_t j = 0; j < n1; ++j) {
        for (std::size_t a = 0; a < n2; ++a) {
          for (std::size_t b = 0; b < n2; ++b) {
            out[i * n2 + a][j * n2 + b] = transposed == TensorFactor::First
                ? (*this)[j * n2 + a][i * n2 + b] : (*this)[i * n2 + b][j * n2 + a];
          }
        }
      }
    }
    return out;
  }

  // Get transposed matrix
  // output_{ij} = A_{ji}
  MMatrix Transpose() const {
    MMatrix out(cols, rows, UninitializedTag{});

    for (std::size_t i = 0; i < rows; ++i) {
      const T *src_row = data + cols * i;

      for (std::size_t j = 0; j < cols; ++j) { out.data[rows * j + i] = src_row[j]; }
    }

    return out;
  }

  // Exchange the second and third indices of one flattened rank-three tensor
  // output_{(i,k),j} = input_{(i,j),k}
  MMatrix FlatTensorTranspose(const std::size_t first, const std::size_t second, const std::size_t third) const {
    if ((second != 0 && first > std::numeric_limits<std::size_t>::max() / second) ||
        (third != 0 && first > std::numeric_limits<std::size_t>::max() / third) || rows != first * second ||
        cols != third) {
      throw std::invalid_argument("MMatrix::FlatTensorTranspose: dimensions disagree");
    }
    MMatrix output(first * third, second, UninitializedTag{});
    for (std::size_t i = 0; i < first; ++i) {
      for (std::size_t j = 0; j < second; ++j) {
        for (std::size_t k = 0; k < third; ++k) {
          output.data[(i * third + k) * second + j] = data[(i * second + j) * third + k];
        }
      }
    }
    return output;
  }

  // Get conjugate transposed matrix (dagger)
  // A^dagger_{ij} = A_{ji}^*
  MMatrix Dagger() const { return ConjTranspose(); }

  // Get conjugate transposed matrix (dagger)
  // A^dagger_{ij} = A_{ji}^*
  MMatrix ConjTranspose() const {
    MMatrix out(cols, rows, UninitializedTag{});

    for (std::size_t i = 0; i < rows; ++i) {
      const T *src_row = data + cols * i;

      for (std::size_t j = 0; j < cols; ++j) { out.data[rows * j + i] = Conjugate(src_row[j]); }
    }

    return out;
  }

  // Get matrix with conjugate elements (no transpose)
  // output_{ij} = A_{ij}^*
  MMatrix Conj() const {
    MMatrix out(rows, cols, UninitializedTag{});

    for (std::size_t i = 0; i < rows * cols; ++i) { out.data[i] = Conjugate(data[i]); }

    return out;
  }

  // Get diagonal vector
  std::vector<T> GetDiag() const {
    if (rows != cols) {
      throw std::invalid_argument("MMatrix::GetDiag: Only defined for square matrices, rows = " + std::to_string(rows) +
                                  " , cols = " + std::to_string(cols));
    }
    std::vector<T> diagonal(rows);
    for (std::size_t i = 0; i < rows; ++i) { diagonal[i] = data[i * cols + i]; }
    return diagonal;
  }

  // Compute one rectangular matrix block
  MMatrix Submatrix(const std::size_t first_row, const std::size_t first_col, const std::size_t row_count,
                    const std::size_t col_count) const {
    if (first_row > rows || first_col > cols || row_count > rows - first_row || col_count > cols - first_col) {
      throw std::out_of_range("MMatrix::Submatrix: block outside matrix");
    }
    MMatrix out = Uninitialized(row_count, col_count);
    if (col_count == 0) { return out; }
    for (std::size_t row = 0; row < row_count; ++row) {
      std::copy_n(data + (first_row + row) * cols + first_col, col_count, out.data + row * col_count);
    }
    return out;
  }

  // Set one square block with a scalar coefficient
  template <typename Scalar>
  void SetBlock(const std::size_t block_row, const std::size_t block_col, const MMatrix &block,
                const Scalar &coefficient) {
    const std::size_t dimension = block.rows;
    if (block.rows != block.cols) { throw std::invalid_argument("MMatrix::SetBlock: block must be square"); }
    const std::size_t first_row = CheckedElementCount(block_row, dimension);
    const std::size_t first_col = CheckedElementCount(block_col, dimension);
    if (first_row > rows || dimension > rows - first_row || first_col > cols || dimension > cols - first_col) {
      throw std::invalid_argument("MMatrix::SetBlock: incompatible block dimensions");
    }
    for (std::size_t row = 0; row < dimension; ++row) {
      for (std::size_t col = 0; col < dimension; ++col) {
        data[(first_row + row) * cols + first_col + col] = T(coefficient) * block.data[row * dimension + col];
      }
    }
  }

  // Set one square block with unit coefficient
  void SetBlock(const std::size_t block_row, const std::size_t block_col, const MMatrix &block) {
    SetBlock(block_row, block_col, block, T(1));
  }

  // Reshape the existing row-major entries without copying or changing their order
  void Reshape(std::size_t r, std::size_t c) {
    if ((c != 0 && r > std::numeric_limits<std::size_t>::max() / c) || r * c != rows * cols) {
      throw std::invalid_argument("MMatrix::Reshape: element counts disagree");
    }
    rows = r;
    cols = c;
  }

  // Flatten the matrix into a single vector
  std::vector<T> Flatten() const {
    std::vector<T> flat(rows * cols);
    std::copy_n(data, rows * cols, flat.begin());
    return flat;
  }

  // Get trace
  // Tr(A) = sum_i A_{ii}
  T Trace() const {
    if (rows != cols) {
      throw std::invalid_argument("MMatrix::Trace: Only defined for square matrices, rows = " + std::to_string(rows) +
                                  " , cols = " + std::to_string(cols));
    }
    T sum = T(0);
    for (std::size_t i = 0; i < rows; ++i) { sum += data[cols * i + i]; }
    return sum;
  }
  T Tr() const { return Trace(); }

  void Print(const std::string &name = "") const {
    std::cout << "MMatrix::Print: " << name << " [" << rows << " x " << cols << "]" << std::endl;
    std::cout << std::setprecision(4);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) { std::cout << (*this)[i][j] << "\t"; }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  // Print real and imaginary matrix parts separately
  void PrintSeparate() const {
    std::cout << "Re:" << std::endl;
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        const std::string delimiter = col + 1 < cols ? ", " : "";
        std::printf("%6.3f%s ", std::real(data[row * cols + col]), delimiter.c_str());
      }
      std::cout << std::endl;
    }
    std::cout << std::endl << "Im:" << std::endl;
    for (std::size_t row = 0; row < rows; ++row) {
      for (std::size_t col = 0; col < cols; ++col) {
        const std::string delimiter = col + 1 < cols ? ", " : "";
        std::printf("%6.3f%s ", std::imag(data[row * cols + col]), delimiter.c_str());
      }
      std::cout << std::endl;
    }
    std::cout << std::endl;
  }

  // Overload the << operator for printing the matrix
  friend std::ostream &operator<<(std::ostream &os, const MMatrix<T> &matrix) {
    os << "MMatrix << cout: "
       << " [" << matrix.rows << " x " << matrix.cols << "]" << std::endl;
    os << std::setprecision(4);
    for (std::size_t i = 0; i < matrix.rows; ++i) {
      for (std::size_t j = 0; j < matrix.cols; ++j) {
        os << matrix[i][j] << "\t";  // Print each element with tab separation
      }
      os << std::endl;  // Newline after each row
    }
    os << std::endl;
    return os;  // Support chained stream output (e.g., cout << matrix1 <<
                // matrix2)
  }

  // Size operators
  std::size_t size_row() const { return rows; }
  std::size_t size_col() const { return cols; }

  // Compute one mutable matrix row without copying its elements
  std::span<T> Row(const std::size_t row) { return std::span<T>((*this)[row], cols); }

  // Compute one immutable matrix row without copying its elements
  std::span<const T> Row(const std::size_t row) const { return std::span<const T>((*this)[row], cols); }

  // Compute every mutable matrix element in row-major order without copying
  std::span<T> Elements() { return std::span<T>(data, rows * cols); }

  // Compute every immutable matrix element in row-major order without copying
  std::span<const T> Elements() const { return std::span<const T>(data, rows * cols); }

 private:
  template <typename U>
  friend class MMatrix;

  // Validate an ordered bipartite shape without overflowing the dimension product
  void ValidateTensorShape(std::size_t n1, std::size_t n2) const {
    if (rows == 0 || rows != cols || n1 == 0 || n2 == 0 || rows % n1 != 0 || rows / n1 != n2) {
      throw std::invalid_argument("MMatrix: bipartite dimensions disagree");
    }
  }

  static MMatrix Uninitialized(std::size_t r, std::size_t c) { return MMatrix(r, c, UninitializedTag{}); }

  struct UninitializedTag {};
  MMatrix(std::size_t r, std::size_t c, UninitializedTag) { Allocate(r, c); }

  // Copy data from a to this matrix after Resize
  void Copy(const MMatrix &a) {
    if (rows == 0 || cols == 0) { return; }
    std::copy_n(a.data, rows * cols, data);
  }

  void Allocate(std::size_t r, std::size_t c) {
    const std::size_t count = CheckedElementCount(r, c);
    rows                    = r;
    cols                    = c;
    if (count == 0) {
      data = nullptr;
      return;
    }
    data = new T[count];
  }

  // Compute a checked matrix element count for two dimensions
  static std::size_t CheckedElementCount(const std::size_t r, const std::size_t c) {
    if (r != 0 && c > std::numeric_limits<std::size_t>::max() / r) {
      throw std::length_error("MMatrix: matrix dimensions overflow");
    }
    return r * c;
  }

  // Validate a source-to-destination row map as one exact permutation
  template <typename IndexContainer>
  static void ValidateDestinationRows(const IndexContainer &destination_rows, const std::size_t row_count,
                                      const std::string_view context) {
    if (!std::cmp_equal(destination_rows.size(), row_count)) {
      throw std::invalid_argument(std::string(context) + ": destination row dimensions disagree");
    }
    for (std::size_t source_row = 0; source_row < row_count; ++source_row) {
      const auto destination = destination_rows[source_row];
      if (std::cmp_less(destination, 0) || !std::cmp_less(destination, row_count)) {
        throw std::invalid_argument(std::string(context) + ": destination row outside permutation");
      }
      for (std::size_t previous = 0; previous < source_row; ++previous) {
        if (destination == destination_rows[previous]) {
          throw std::invalid_argument(std::string(context) + ": destination rows are not a permutation");
        }
      }
    }
  }

  template <typename Container>
  std::size_t GetRectangularCols(const Container &list) {
    if (list.size() == 0) { return 0; }

    const std::size_t expected_cols = list.begin()->size();
    for (const auto &row : list) {
      if (row.size() != expected_cols) {
        throw std::invalid_argument("MMatrix: initializer-list rows must all have the same length");
      }
    }
    return expected_cols;
  }

  // Helper function to handle conjugate for different types
  template <typename U = T>
  typename std::enable_if<!std::is_arithmetic<U>::value, U>::type Conjugate(const U &value) const {
    return std::conj(value);  // For complex numbers
  }

  // Preserve real values under complex conjugation
  template <typename U = T>
  typename std::enable_if<std::is_arithmetic<U>::value, U>::type Conjugate(const U &value) const {
    return value;  // For real numbers (identity function)
  }

  // Compute true for exact arithmetic zero values used by sparse matrix products
  static bool IsZero(const T &value) { return math::IsZero(value); }

  // Compute squared absolute value without a sqrt round trip
  template <typename U = T>
  static typename std::enable_if<std::is_arithmetic<U>::value, double>::type Abs2(const U &value) {
    const double x = static_cast<double>(value);
    return x * x;
  }

  // Compute squared complex absolute value without a sqrt round trip
  template <typename U = T>
  static typename std::enable_if<!std::is_arithmetic<U>::value, double>::type Abs2(const U &value) {
    return std::norm(value);
  }

  std::size_t rows;
  std::size_t cols;
  T          *data;
};

// Generic vector and matrix algebra

namespace detail {

template <typename T>
struct IsComplex : std::false_type {};
template <typename T>
struct IsComplex<std::complex<T>> : std::true_type {};

template <typename T>
struct IsFloatingScalar : std::is_floating_point<std::decay_t<T>> {};
template <typename T>
struct IsFloatingScalar<std::complex<T>> : std::is_floating_point<T> {};

// Compute the scalar conjugate used by Hermitian vector operations
// Conjugate(z) = z^*
template <typename T>
inline auto Conjugate(const T &value) {
  if constexpr (IsComplex<std::decay_t<T>>::value) {
    return std::conj(value);
  } else {
    return value;
  }
}

// Compute whether one real or complex scalar is finite
template <typename T>
inline bool IsFinite(const T &value) {
  if constexpr (IsComplex<std::decay_t<T>>::value) {
    return std::isfinite(value.real()) && std::isfinite(value.imag());
  } else {
    return std::isfinite(value);
  }
}

// Compute the squared magnitude of one real or complex scalar
// SquaredMagnitude(z) = z z^*
template <typename T>
inline auto SquaredMagnitude(const T &value) {
  if constexpr (IsComplex<std::decay_t<T>>::value) {
    return std::norm(value);
  } else {
    return value * value;
  }
}

}  // namespace detail

// Compute the Euclidean or Hermitian inner product of two iterator ranges
// <x,y> = sum_i x_i^* y_i
template <typename FirstIterator, typename SecondIterator>
inline auto InnerProduct(FirstIterator first, const FirstIterator last, SecondIterator second) {
  using First  = typename std::iterator_traits<FirstIterator>::value_type;
  using Second = typename std::iterator_traits<SecondIterator>::value_type;
  using Result = decltype(detail::Conjugate(std::declval<First>()) * std::declval<Second>());
  Result result{};
  for (; first != last; ++first, ++second) { result += detail::Conjugate(*first) * *second; }
  return result;
}

// Compute the Euclidean or Hermitian inner product of two containers
// <x,y> = sum_i x_i^* y_i
template <typename FirstContainer, typename SecondContainer>
inline auto InnerProduct(const FirstContainer &first, const SecondContainer &second) {
  if (first.size() != second.size()) { throw std::invalid_argument("gra::InnerProduct: vector dimensions disagree"); }
  return InnerProduct(std::begin(first), std::end(first), std::begin(second));
}

// Compute log(sum_i exp(a_i + b_i)) while preserving zero support as -infinity
template <typename FirstContainer, typename SecondContainer>
inline long double LogProductSum(const FirstContainer &first, const SecondContainer &second) {
  if (first.size() != second.size()) { throw std::invalid_argument("gra::LogProductSum: dimensions disagree"); }
  long double maximum = -std::numeric_limits<long double>::infinity();
  for (std::size_t i = 0; i < first.size(); ++i) {
    const long double a = first[i];
    const long double b = second[i];
    if (std::isnan(a) || std::isnan(b) || (std::isinf(a) && !std::signbit(a)) || (std::isinf(b) && !std::signbit(b))) {
      throw std::invalid_argument("gra::LogProductSum: invalid logarithm");
    }
    maximum = std::max(maximum, a + b);
  }
  if (!std::isfinite(maximum)) { return maximum; }
  long double sum = 0.0L;
  for (std::size_t i = 0; i < first.size(); ++i) {
    sum += std::exp((static_cast<long double>(first[i]) + second[i]) - maximum);
  }
  return maximum + std::log(sum);
}

// Compute the squared Euclidean norm of one vector-like container
// ||x||_2^2 = sum_i |x_i|^2
template <typename Container>
inline auto SquaredNorm(const Container &x) {
  using T      = typename Container::value_type;
  using Result = decltype(detail::SquaredMagnitude(std::declval<T>()));
  Result result{};
  for (const auto &value : x) { result += detail::SquaredMagnitude(value); }
  return result;
}

// Compute the sum of one vector-like container
// output = sum_i x_i
template <typename Container>
inline auto Sum(const Container &x) {
  using T = typename Container::value_type;
  T result{};
  for (const auto &value : x) { result += value; }
  return result;
}

// Compute true when every real or complex container entry is finite
template <typename Container>
inline bool AllFinite(const Container &x) {
  return std::all_of(std::begin(x), std::end(x), [](const auto &value) { return detail::IsFinite(value); });
}

// Accumulate a cached diagonal matrix product without dense matrix traversal
// y_i <- y_i + c d_i x_i
template <typename Diagonal, typename U, typename V, typename Scalar>
inline void DiagonalMultiplyAdd(const std::span<const Diagonal> diagonal, const std::span<const U> source,
                                const std::span<V> target, const Scalar &scale) {
  if (diagonal.size() != source.size() || diagonal.size() != target.size()) {
    throw std::invalid_argument("gra::DiagonalMultiplyAdd: vector dimensions disagree");
  }
  for (std::size_t index = 0; index < diagonal.size(); ++index) {
    const Diagonal value = diagonal[index];
    if (math::IsZero(value)) { continue; }
    target[index] += scale * value * source[index];
  }
}

// Scale one vector-like container in place
// x_i <- c x_i
template <typename Container, typename Scalar>
inline void Scale(Container &x, const Scalar &scale) {
  std::for_each(std::begin(x), std::end(x), [&scale](auto &value) { value *= scale; });
}

// Compute the additive inverse of one vector-like container
// output_i = -x_i
template <typename Container>
requires(!std::ranges::view<Container>) inline Container Negated(const Container &input) {
  Container output = input;
  Scale(output, -1.0);
  return output;
}

// Compute the cross product of two fixed three-vectors
// (a cross b)_i = epsilon_{ijk} a_j b_k
template <typename T>
inline std::array<T, 3> CrossProduct(const std::array<T, 3> &left, const std::array<T, 3> &right) {
  return {left[1] * right[2] - left[2] * right[1], left[2] * right[0] - left[0] * right[2],
          left[0] * right[1] - left[1] * right[0]};
}

// Compute an L2 normalized copy of one vector-like container
// output = x/sqrt(sum_i |x_i|^2)
template <typename Container>
requires(!std::ranges::view<Container>) inline Container NormalizedL2(const Container &x) {
  using Value = typename Container::value_type;
  static_assert(detail::IsFloatingScalar<Value>::value, "gra::NormalizedL2 requires floating or complex values");
  using Real = decltype(std::real(Value{}));
  Real scale = Real(0);
  for (const auto &value : x) {
    if (!detail::IsFinite(value)) { throw std::invalid_argument("gra::NormalizedL2 requires finite values"); }
    scale = std::max({scale, std::abs(std::real(value)), std::abs(std::imag(value))});
  }
  if (!(scale > Real(0))) {
    throw std::invalid_argument("gra::NormalizedL2 requires a finite nonzero norm");
  }
  Container output = x;
  for (auto &value : output) { value /= scale; }
  Scale(output, Real(1) / std::sqrt(SquaredNorm(output)));
  return output;
}

// Compute a positive-sum normalized copy of one vector-like container
// output_i = x_i/sum_j x_j
template <typename Container>
requires(!std::ranges::view<Container>) inline Container NormalizedSum(const Container &x) {
  using Value = typename Container::value_type;
  static_assert(std::is_floating_point_v<Value>, "gra::NormalizedSum requires real floating values");
  Value scale = Value(0);
  for (const auto value : x) {
    if (!std::isfinite(value)) { throw std::invalid_argument("gra::NormalizedSum requires finite values"); }
    scale = std::max(scale, std::abs(value));
  }
  if (!(scale > Value(0))) { throw std::invalid_argument("gra::NormalizedSum requires a positive sum"); }
  Container  output = x;
  for (auto &value : output) { value /= scale; }
  const auto total = Sum(output);
  if (!std::isfinite(total) || total <= 0.0) {
    throw std::invalid_argument("gra::NormalizedSum requires a finite positive sum");
  }
  for (auto &value : output) { value /= total; }
  return output;
}

// Accumulate one scaled vector-like container into another
// target_i <- target_i + c source_i
template <typename TargetContainer, typename SourceContainer, typename Scalar>
inline void AddScaled(TargetContainer &target, const SourceContainer &source, const Scalar &scale) {
  if (target.size() != source.size()) { throw std::invalid_argument("gra::AddScaled: vector dimensions disagree"); }
  std::transform(std::begin(target), std::end(target), std::begin(source), std::begin(target),
                 [&scale](const auto &left, const auto &right) { return left + scale * right; });
}

// Project one vector onto an orthonormal vector basis
template <typename Basis, typename Container>
inline Container ProjectOrthonormal(const Basis &basis, const Container &source) {
  Container output(source.size(), typename Container::value_type{});
  for (const auto &vector : basis) { AddScaled(output, vector, InnerProduct(vector, source)); }
  return output;
}

// Add one linearly independent vector to an orthonormal vector basis
template <typename Basis, typename Container>
inline bool AddOrthonormal(Basis &basis, const Container &candidate, double tolerance) {
  if (!std::isfinite(tolerance) || tolerance < 0.0) {
    throw std::invalid_argument("gra::AddOrthonormal requires a finite nonnegative tolerance");
  }
  Container vector = candidate;
  for (const auto &unit : basis) { AddScaled(vector, unit, -InnerProduct(unit, vector)); }
  const double norm2 = SquaredNorm(vector);
  if (!std::isfinite(norm2)) { throw std::invalid_argument("gra::AddOrthonormal requires finite values"); }
  if (norm2 <= tolerance * tolerance) { return false; }
  Scale(vector, 1.0 / std::sqrt(norm2));
  basis.push_back(std::move(vector));
  return true;
}

// Compute the elementwise sum of two equally sized containers
// output_i = x_i + y_i
template <typename Container>
requires(!std::ranges::view<Container>) inline Container Add(const Container &left, const Container &right) {
  Container output = left;
  AddScaled(output, right, 1.0);
  return output;
}

// Compute the elementwise difference of two equally sized containers
// output_i = x_i - y_i
template <typename Container>
requires(!std::ranges::view<Container>) inline Container Subtract(const Container &left, const Container &right) {
  Container output = left;
  AddScaled(output, right, -1.0);
  return output;
}

// Compute the nonconjugating bilinear product of two equal-size containers
// output = sum_i x_i y_i
template <typename FirstContainer, typename SecondContainer>
inline auto BilinearProduct(const FirstContainer &first, const SecondContainer &second) {
  if (first.size() != second.size()) {
    throw std::invalid_argument("gra::BilinearProduct: vector dimensions disagree");
  }
  using First  = typename FirstContainer::value_type;
  using Second = typename SecondContainer::value_type;
  using Result = decltype(std::declval<First>() * std::declval<Second>());
  Result result{};
  for (std::size_t index = 0; index < first.size(); ++index) { result += first[index] * second[index]; }
  return result;
}

// Compute the Minkowski bilinear product with metric (+,-,-,-)
// x.y = x^0 y^0 - x^1 y^1 - x^2 y^2 - x^3 y^3
template <typename FirstContainer, typename SecondContainer>
inline auto MinkowskiProduct(const FirstContainer &first, const SecondContainer &second) {
  if (first.size() != 4 || second.size() != 4) {
    throw std::invalid_argument("gra::MinkowskiProduct: Lorentz vectors must have four components");
  }
  using First  = typename FirstContainer::value_type;
  using Second = typename SecondContainer::value_type;
  using Result = decltype(std::declval<First>() * std::declval<Second>());
  if constexpr (std::is_floating_point_v<Result>) {
    Result result = std::fma(static_cast<Result>(first[0]), static_cast<Result>(second[0]),
                             -static_cast<Result>(first[1]) * static_cast<Result>(second[1]));
    result        = std::fma(-static_cast<Result>(first[2]), static_cast<Result>(second[2]), result);
    return std::fma(-static_cast<Result>(first[3]), static_cast<Result>(second[3]), result);
  } else {
    Result result = first[0] * second[0];
    for (std::size_t index = 1; index < 4; ++index) { result -= first[index] * second[index]; }
    return result;
  }
}

// Compute a nonconjugating bilinear product restricted by one common mask
// output = sum_{i: mask_i} x_i y_i
template <typename FirstContainer, typename SecondContainer, typename Mask>
inline auto MaskedBilinearProduct(const FirstContainer &first, const SecondContainer &second, const Mask &mask) {
  if (first.size() != second.size() || first.size() != mask.size()) {
    throw std::invalid_argument("gra::MaskedBilinearProduct: vector dimensions disagree");
  }
  using First  = typename FirstContainer::value_type;
  using Second = typename SecondContainer::value_type;
  using Result = decltype(std::declval<First>() * std::declval<Second>());
  Result result{};
  for (std::size_t index = 0; index < first.size(); ++index) {
    if (mask[index]) { result += first[index] * second[index]; }
  }
  return result;
}

// Compute the entrywise complex conjugate of one vector-like container
// output_i = x_i^*
template <typename Container>
requires(!std::ranges::view<Container>) inline Container Conjugated(const Container &input) {
  Container output = input;
  for (auto &value : output) { value = detail::Conjugate(value); }
  return output;
}

// Compute the L1 distance of two equal-size real or complex containers
// d_1(x,y) = sum_i |x_i-y_i|
template <typename FirstContainer, typename SecondContainer>
inline double L1Distance(const FirstContainer &first, const SecondContainer &second) {
  if (first.size() != second.size()) { throw std::invalid_argument("gra::L1Distance: vector dimensions disagree"); }
  double result = 0.0;
  for (std::size_t index = 0; index < first.size(); ++index) { result += std::abs(first[index] - second[index]); }
  return result;
}

// Compute the Kronecker product of two vector-like containers
// output_{iN+j} = x_i y_j
template <typename FirstContainer, typename SecondContainer>
inline auto KroneckerProduct(const FirstContainer &first, const SecondContainer &second) {
  using Result = decltype(first[0] * second[0]);
  if (!second.empty() && first.size() > std::numeric_limits<std::size_t>::max() / second.size()) {
    throw std::invalid_argument("gra::KroneckerProduct: dimension overflow");
  }
  std::vector<Result> output(first.size() * second.size());
  for (std::size_t i = 0; i < first.size(); ++i) {
    for (std::size_t j = 0; j < second.size(); ++j) { output[i * second.size() + j] = first[i] * second[j]; }
  }
  return output;
}

// Contract a Kronecker vector with one flat tensor without allocating it
// output = sum_{ij} x_i y_j T_{iN+j}
template <typename FirstContainer, typename SecondContainer, typename TensorContainer>
inline auto KroneckerBilinearProduct(const FirstContainer &first, const SecondContainer &second,
                                     const TensorContainer &tensor) {
  if (!second.empty() && first.size() > std::numeric_limits<std::size_t>::max() / second.size()) {
    throw std::invalid_argument("gra::KroneckerBilinearProduct: dimension overflow");
  }
  if (first.size() * second.size() != tensor.size()) {
    throw std::invalid_argument("gra::KroneckerBilinearProduct: tensor dimensions disagree");
  }
  using First  = typename FirstContainer::value_type;
  using Second = typename SecondContainer::value_type;
  using Tensor = typename TensorContainer::value_type;
  using Result = decltype(std::declval<First>() * std::declval<Second>() * std::declval<Tensor>());
  Result result{};
  for (std::size_t first_index = 0; first_index < first.size(); ++first_index) {
    for (std::size_t second_index = 0; second_index < second.size(); ++second_index) {
      result += first[first_index] * second[second_index] * tensor[first_index * second.size() + second_index];
    }
  }
  return result;
}

// Compute the elementwise product of two vector-like containers
// output_i = x_i y_i
template <typename FirstContainer, typename SecondContainer>
inline auto HadamardProduct(const FirstContainer &first, const SecondContainer &second) {
  using Result = decltype(first[0] * second[0]);
  if (first.size() != second.size()) {
    throw std::invalid_argument("gra::HadamardProduct: vector dimensions disagree");
  }
  std::vector<Result> output(first.size());
  std::transform(std::begin(first), std::end(first), std::begin(second), output.begin(),
                 [](const auto &left, const auto &right) { return left * right; });
  return output;
}

// Compute the outer product of two vector-like containers
// output_{ij} = x_i y_j
template <typename FirstContainer, typename SecondContainer>
inline auto OuterProduct(const FirstContainer &first, const SecondContainer &second) {
  using Result = decltype(first[0] * second[0]);
  MMatrix<Result> output(first.size(), second.size(), Result{});
  for (std::size_t i = 0; i < first.size(); ++i) {
    for (std::size_t j = 0; j < second.size(); ++j) { output(i, j) = first[i] * second[j]; }
  }
  return output;
}

// Compute a row-mapped Kronecker product of two vector-like containers
// output_{p(iN+j)} = x_i y_j
template <typename FirstContainer, typename SecondContainer, typename IndexContainer>
inline auto MappedKroneckerProduct(const FirstContainer &first, const SecondContainer &second,
                                   const IndexContainer &destination_rows) {
  using First  = typename FirstContainer::value_type;
  using Second = typename SecondContainer::value_type;
  using Result = std::decay_t<decltype(std::declval<First>() * std::declval<Second>())>;

  if (!second.empty() && first.size() > std::numeric_limits<std::size_t>::max() / second.size()) {
    throw std::invalid_argument("gra::MappedKroneckerProduct: dimension overflow");
  }
  const std::size_t output_size = first.size() * second.size();
  if (!std::cmp_equal(destination_rows.size(), output_size)) {
    throw std::invalid_argument("gra::MappedKroneckerProduct: destination dimensions disagree");
  }

  std::vector<Result> output(output_size);
  std::vector<bool>   assigned(output_size, false);
  for (std::size_t first_index = 0; first_index < first.size(); ++first_index) {
    for (std::size_t second_index = 0; second_index < second.size(); ++second_index) {
      const std::size_t source      = first_index * second.size() + second_index;
      const std::size_t destination = static_cast<std::size_t>(destination_rows[source]);
      if (destination >= output_size || assigned[destination]) {
        throw std::invalid_argument(
            "gra::MappedKroneckerProduct: destination rows are not a "
            "permutation");
      }
      assigned[destination] = true;
      output[destination]   = first[first_index] * second[second_index];
    }
  }
  return output;
}

// Compute Cartesian coordinate matrices for two one-dimensional grids
// X_{ij} = x_j, Y_{ij} = y_i
template <typename T>
inline auto MeshGrid(const std::vector<T> &x, const std::vector<T> &y) {
  MMatrix<T> X(y.size(), x.size());
  MMatrix<T> Y(y.size(), x.size());
  for (std::size_t row = 0; row < y.size(); ++row) {
    for (std::size_t col = 0; col < x.size(); ++col) {
      X(row, col) = x[col];
      Y(row, col) = y[row];
    }
  }
  return std::make_pair(std::move(X), std::move(Y));
}

// Compute upper_basis^dagger source conjugate(lower_basis) in row-major order
// output_{ab} = sum_{ik} U_{ia}^* source_{ik} L_{kb}^*
template <typename Source, typename U, typename V>
inline auto TensorProductProjection(const Source &source, const MMatrix<U> &upper_basis,
                                    const MMatrix<V> &lower_basis) {
  using T                      = typename Source::value_type;
  using Intermediate           = decltype(std::declval<T>() * detail::Conjugate(std::declval<V>()));
  using Result                 = decltype(detail::Conjugate(std::declval<U>()) * std::declval<Intermediate>());
  const std::size_t upper_rows = upper_basis.size_row();
  const std::size_t lower_rows = lower_basis.size_row();
  if ((!lower_rows && upper_rows) ||
      (lower_rows && upper_rows > std::numeric_limits<std::size_t>::max() / lower_rows) ||
      source.size() != upper_rows * lower_rows) {
    throw std::invalid_argument("gra::TensorProductProjection: dimensions disagree");
  }

  const std::size_t         lower_cols = lower_basis.size_col();
  std::vector<Intermediate> lower_projection(upper_rows * lower_cols, Intermediate{});
  for (std::size_t i = 0; i < upper_rows; ++i) {
    Intermediate *projected_row = lower_projection.data() + i * lower_cols;
    const T      *source_row    = source.data() + i * lower_rows;
    for (std::size_t k = 0; k < lower_rows; ++k) {
      const V *basis_row = lower_basis[k];
      for (std::size_t b = 0; b < lower_cols; ++b) {
        projected_row[b] += source_row[k] * detail::Conjugate(basis_row[b]);
      }
    }
  }

  std::vector<Result> output(upper_basis.size_col() * lower_cols, Result{});
  for (std::size_t i = 0; i < upper_rows; ++i) {
    const U            *basis_row     = upper_basis[i];
    const Intermediate *projected_row = lower_projection.data() + i * lower_cols;
    for (std::size_t a = 0; a < upper_basis.size_col(); ++a) {
      Result    *output_row  = output.data() + a * lower_cols;
      const auto coefficient = detail::Conjugate(basis_row[a]);
      for (std::size_t b = 0; b < lower_cols; ++b) { output_row[b] += coefficient * projected_row[b]; }
    }
  }
  return output;
}

// Compute the rank one projector formed from one normalized vector
// P_{ij} = v_i v_j^*
template <typename Container>
inline auto RankOneProjector(const Container &vector) {
  using T      = typename Container::value_type;
  using Result = decltype(std::declval<T>() * detail::Conjugate(std::declval<T>()));
  MMatrix<Result> projector(vector.size(), vector.size(), Result{});
  for (std::size_t i = 0; i < vector.size(); ++i) {
    for (std::size_t j = 0; j < vector.size(); ++j) { projector(i, j) = vector[i] * detail::Conjugate(vector[j]); }
  }
  return projector;
}

// Convert between vector and valarray storage
template <typename T>
inline std::vector<T> valarray2vector(const std::valarray<T> &x) {
  std::vector<T> y(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) { y[k] = x[k]; }
  return y;
}
template <typename T>
inline std::valarray<T> vector2valarray(const std::vector<T> &x) {
  std::valarray<T> y(x.size());
  for (std::size_t k = 0; k < x.size(); ++k) { y[k] = x[k]; }
  return y;
}

// Cumulative sum
// sumvec_i = sum_{j=0}^i x_j
template <typename T>
inline void CumSum(const std::vector<T> &x, std::vector<T> &sumvec) {
  const std::size_t N = x.size();
  sumvec.resize(N, 0.0);
  if (N == 0) { return; }
  sumvec[0] = x[0];  // First
  for (std::size_t i = 1; i < N; ++i) { sumvec[i] = sumvec[i - 1] + x[i]; }
}

// Template print
template <template <typename T> class container_type, class value_type>
inline void PrintArray(container_type<value_type> x, std::string name) {
  std::cout << "PrintArray: " << name << std::endl;
  for (unsigned int i = 0; i < x.size(); ++i) { std::cout << x[i]; }
  std::cout << std::endl << std::endl;
}

// Multiply a matrix by a scalar from the left
// output_{ij} = c A_{ij}
template <typename T>
MMatrix<T> operator*(const T &lhs, MMatrix<T> rhs) {
  rhs *= lhs;
  return rhs;
}

}  // namespace gra

#endif
