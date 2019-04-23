// Minimal templated matrix class [HEADER ONLY class]
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATRIX_H
#define MMATRIX_H

#include <initializer_list>
#include <iomanip>
#include <iostream>

namespace gra {

// Constant size matrix
template <typename T> class MMatrix {
public:
  MMatrix() : rows(0), cols(0) { data = nullptr; }

  MMatrix(std::size_t r, std::size_t c) {
    rows = r;
    cols = c;
    data = new T[rows * cols];
    std::fill(data, data + rows * cols, T()); // No initialization
  }
  MMatrix(std::size_t r, std::size_t c, T value) {
    rows = r;
    cols = c;
    data = new T[rows * cols];
    std::fill(data, data + rows * cols, T(value)); // Initialization
  }
  MMatrix(std::size_t r, std::size_t c, const std::string &special) {
    rows = r;
    cols = c;
    data = new T[rows * cols];

    if (special == "eye") {
      std::fill(data, data + rows * cols, T(0.0));
      Identity();
    } else if (special == "minkowski") {
      std::fill(data, data + rows * cols, T(0.0));
      Minkowski();
    } else {
      throw std::invalid_argument("MMatrix: Unknown initialization string:" +
                                  special);
    }
  }
  // Set matrix to identity
  void Identity() {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        this->operator()(i, j) = (i == j) ? 1.0 : 0.0;
      }
    }
  }
  // Set matrix to Minkowski metric
  void Minkowski() {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        double sign = (i > 0 && j > 0) ? -1.0 : 1.0;
        this->operator()(i, j) = (i == j) ? sign * 1.0 : 0.0;
      }
    }
  }

  // For initializing with a = {row0-vector, row1-vector, ...}
  // where each row is std::vector<T>
  //
  // Row-major order!
  MMatrix(std::initializer_list<std::vector<T>> list) {
    rows = list.size();
    cols = list.begin()->size();
    data = new T[rows * cols];

    std::size_t i = 0;
    for (const auto &v : list) {
      for (std::size_t j = 0; j < cols; ++j) {
        data[cols * i + j] = v[j];
      }
      ++i;
    }
  }

  // For initializing with a = { {}, {}, {} }
  //
  // Row-major order!
  MMatrix(std::initializer_list<std::initializer_list<T>> list) {
    rows = list.size();
    cols = list.begin()->size();
    data = new T[rows * cols];

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

  // Copy constructor
  MMatrix(const MMatrix &a) {
    rows = a.rows;
    cols = a.cols;
    data = new T[rows * cols];

    // Copy all elements
    Copy(a);
  }

  ~MMatrix() {
    // delete dynamically allocated memory
    delete[] data;
  }

  // Assignment operator
  MMatrix &operator=(const MMatrix &rhs) {
    if (data != rhs.data && rhs.data != nullptr) {
      ReSize(rhs.rows, rhs.cols);
      Copy(rhs);
    }
    return *this;
  }
  // For indexing with [i][j]
  // Row-major order!
  T *operator[](const std::size_t &row) { return data + cols * row; }
  const T *operator[](const std::size_t &row) const {
    return data + cols * row;
  }

  // Row-major order!
  T &operator()(std::size_t i, std::size_t j) {
    if (i >= rows || j >= cols) {
      throw std::out_of_range("MMatrix:: Index over matrix dimensions!");
    }
    return data[cols * i + j];
  }
  const T &operator()(std::size_t i, std::size_t j) const {
    if (i >= rows || j >= cols) {
      throw std::out_of_range("MMatrix:: Index over matrix dimensions!");
    }
    return data[cols * i + j];
  }

  // ------------------------------------------------------------------
  // Add to the left
  MMatrix &operator+=(const MMatrix &rhs) {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        this->operator()(i, j) += rhs[i][j];
      }
    }
    return *this;
  }
  // Subtract to the left
  MMatrix &operator-=(const MMatrix &rhs) {
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        this->operator()(i, j) -= rhs[i][j];
      }
    }
    return *this;
  }

  // ------------------------------------------------------------------

  // Return negated matrix
  MMatrix operator-() const {
    MMatrix out(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i][j] = -this->operator()(i, j);
      }
    }
    return out;
  }

  // Add two matrices
  MMatrix operator+(const MMatrix &rhs) const {
    MMatrix out(rows, cols);
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      throw std::invalid_argument(
          "MMatrix:: Matrix + Matrix with invalid dimensions");
    }

    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i][j] = this->operator()(i, j) + rhs[i][j];
      }
    }
    return out;
  }
  // Subtract two matrices
  MMatrix operator-(const MMatrix &rhs) const {
    MMatrix out(rows, cols);
    if (rows != rhs.size_row() || cols != rhs.size_col()) {
      throw std::invalid_argument(
          "MMatrix:: Matrix - Matrix with invalid dimensions");
    }

    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i][j] = this->operator()(i, j) - rhs[i][j];
      }
    }
    return out;
  }
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // We do not allow addition or substraction by scalar, for dimensional
  // "safety" reasons (can lead to accidental results), thus only multiplicative
  // operations.

  // Multiply by a scalar
  MMatrix operator*(const T &rhs) const {
    MMatrix out(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i][j] = this->operator()(i, j) * rhs;
      }
    }
    return out;
  }

  // Divide by a scalar
  MMatrix operator/(const T &rhs) const {
    MMatrix out(rows, cols);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i][j] = this->operator()(i, j) / rhs;
      }
    }
    return out;
  }
  // ------------------------------------------------------------------

  // Matrix [this] * Vector [rhs] multiplication
  std::vector<T> operator*(const std::vector<T> &rhs) const {
    if (rhs.size() != cols) {
      throw std::invalid_argument(
          "MMatrix:: matrix * vector product with invalid dimension on rhs");
    }

    std::vector<T> out(rows, 0.0); // Init!
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[i] += this->operator()(i, j) * rhs[j];
      }
    }
    return out;
  }

  // Matrix [this] * Matrix [rhs] multiplication
  MMatrix operator*(const MMatrix &rhs) const {
    const std::size_t n = rows;
    const std::size_t m = cols;
    const std::size_t p = rhs.size_col();

    MMatrix<T> C(n, p, 0.0); // Note initialization to 0.0!

    if (cols != rhs.size_row()) {
      throw std::invalid_argument(
          "MMatrix:: matrix * matrix product with invalid dimensions");
    }

    // Over [i,j] of C
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < p; ++j) {
        for (std::size_t k = 0; k < m; ++k) {
          C[i][j] += this->operator()(i, k) * rhs[k][j]; // notice plus
        }
      }
    }
    return C;
  }

  // ------------------------------------------------------------------

  // Frobenius norm
  double FrobNorm() const {
    double sum = 0.0;
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        const double value = std::abs(this->operator()(i, j));
        sum += value * value; // |a_ij|^2
      }
    }
    return sum;
  }

  // Get transposed matrix
  MMatrix Transpose() const {
    MMatrix out(cols, rows);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[j][i] = this->operator()(i, j);
      }
    }
    return out;
  }
  // Get conjugate transposed matrix (dagger)
  MMatrix ConjTranspose() const {
    MMatrix out(cols, rows);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        out[j][i] = std::conj(this->operator()(i, j));
      }
    }
    return out;
  }
  MMatrix Dagger() const { return ConjTranspose(); }

  // Get diagonal vector
  std::vector<T> GetDiag() const {
    std::vector<T> d(rows);
    for (std::size_t i = 0; i < rows; ++i) {
      d[i] = this->operator()(i, i);
    }
    return d;
  }

  // Get trace
  T Trace() const {
    T sum = 0.0;
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        sum += this->operator()(i, j);
      }
    }
    return sum;
  }
  T Tr() const { return Trace(); }

  void Print(const std::string &name = "") const {
    std::cout << "MMatrix::Print: " << name << std::endl;
    std::cout << std::setprecision(4);
    for (std::size_t i = 0; i < rows; ++i) {
      for (std::size_t j = 0; j < cols; ++j) {
        std::cout << this->operator()(i, j) << "\t";
      }
      std::cout << std::endl;
    }
  }

  // Size operators
  std::size_t size_row() const { return rows; }
  std::size_t size_col() const { return cols; }

private:
  // Copy data from a to *this (after ReSize)
  void Copy(const MMatrix &a) {
    T *p = data + rows * cols;
    T *q = a.data + rows * cols;
    while (p > data) {
      *--p = *--q;
    }
  }

  // Re-Allocate
  void ReSize(std::size_t r, std::size_t c) {
    if (data != nullptr) {
      delete[] data;
    }
    rows = r;
    cols = c;
    data = new T[rows * cols];
  }

  std::size_t rows;
  std::size_t cols;
  T *data;
};

} // gra namespace ends

#endif