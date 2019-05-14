// Simple matrix/vector operations not included in MMatrix. [HEADER ONLY file]
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MMATOPER_H
#define MMATOPER_H

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"

namespace gra {
namespace matoper {

// MATLAB style meshgrid
//
// X is a matrix with each row    == x
// Y is a matrix with each column == y
template <typename T>
inline void MeshGrid(const std::vector<T>& x, const std::vector<T>& y, MMatrix<T>& X, MMatrix<T>& Y) {
  
  X = MMatrix<T>(y.size(), x.size());
  Y = MMatrix<T>(y.size(), x.size());

  for (std::size_t i = 0; i < X.size_row(); ++i) {
    for (std::size_t j = 0; j < X.size_col(); ++j) {
      X[i][j] = x[j];
    }
  }
  for (std::size_t i = 0; i < Y.size_row(); ++i) {
    for (std::size_t j = 0; j < Y.size_col(); ++j) {
      Y[i][j] = y[i];
    }
  }
}

// diag(x) * A * diag(y) product
template <typename T>
inline MMatrix<T> diagAdiag(const std::vector<T>& x, const MMatrix<T>& A, const std::vector<T>& y) {

  // Over [i,j] of A
  const std::size_t n = A.size_row();
  const std::size_t m = A.size_col();
  
  MMatrix<T> C(n, m);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < m; ++j) {
        C[i][j] = x[i] * A[i][j] * y[j];
    }
  }
  return C;
}

// Return diagonal matrix constructed from a vector
template <typename T>
inline MMatrix<T> Diag(const std::vector<T>& x) {

  MMatrix<T> A(x.size(), x.size(), 0.0);
  for (std::size_t i = 0; i < x.size(); ++i) {
    A[i][i] = x[i];
  }
  return A;
}

// Return L2-unit vector (|x| = 1)
template <typename T>
inline std::vector<T> Unit(const std::vector<T> &x) {
  double norm = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) { norm += gra::math::abs2(x[k]); }
  norm               = gra::math::msqrt(norm);  // l2-norm

  std::vector<T> y = x;
  for (std::size_t k = 0; k < y.size(); ++k) {
    y[k] /= norm;  // normalize
  }
  return y;
}

// Return unit sum normalized (\sum_i x_i = 1)
template <typename T>
inline std::vector<T> Normalized(const std::vector<T> &x) {
  double sum = 0.0;
  for (std::size_t k = 0; k < x.size(); ++k) { sum += x[k]; }
  std::vector<T> y = x;
  if (sum > 0) {
  for (std::size_t k = 0; k < y.size(); ++k) {
    y[k] /= sum;  // normalize
  }
  }
  return y;
}

// Multiply vector by scalar
template <typename T, typename T2>
inline void ScaleVector(std::vector<T> &A, T2 scale) {
  for (std::size_t i = 0; i < A.size(); ++i) { A[i] *= scale; }
}

// Add a constant matrix: A + (full ones) * c
template <typename T, typename T2>
inline void AddConstMat(MMatrix<T> &A, T2 c) {
  for (std::size_t i = 0; i < A.size_row(); ++i) {
    for (std::size_t j = 0; j < A.size_col(); ++j) { A[i][j] += c; }
  }
}

// Vector-Matrix Multiplication c = aB,
// where a is (1 x m)
//       B is (m x p),
//       which gives c (1 x p)
template <typename T>
inline std::vector<T> VecMatMultiply(const std::vector<T> &A, const MMatrix<T> &B) {
  const std::size_t m = A.size();
  const std::size_t p = B.size_col();

  std::vector<T> C(p, 0.0);  // Init with zero!

  // Over [j] of C
  for (std::size_t j = 0; j < p; ++j) {
    for (std::size_t k = 0; k < m; ++k) {
      C[j] += A[k] * B[k][j];  // notice plus
    }
  }

  return C;
}

// Vector-Vector multiplication (strictly not a dot product, which would have
// complex conjugation for the other)
template <typename T>
inline T VecVecMultiply(const std::vector<T> &A, const std::vector<T> &B) {
  T out = 0.0;
  for (std::size_t i = 0; i < A.size(); ++i) {
    out += A[i] * B[i];  // notice plus
  }
  return out;
}

// Vector Conjugate Transpose (\dagger)
// X is (n x 1)
// Y is (1 x n)
template <typename T>
inline std::vector<T> VecDagger(const std::vector<T> &X) {
  const unsigned int n = X.size();
  std::vector<T>     Y(n, 0.0);

  // Conjugate element by element
  for (std::size_t i = 0; i < n; ++i) { Y[i] = std::conj(X[i]); }
  return Y;
}

// Cross product for two 3-vectors
template <typename T>
inline std::vector<T> Cross(const std::vector<T> &a, const std::vector<T> &b) {
  std::vector<T> z = {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2],
                      a[0] * b[1] - a[1] * b[0]};
  return z;
}

// Negate vector
template <typename T>
inline std::vector<T> Negat(const std::vector<T> &a) {
  std::vector<T> z(a.size(), 0.0);
  for (std::size_t k = 0; k < a.size(); ++k) { z[k] = -a[k]; }
  return z;
}

// Subtract two vectors
template <typename T>
inline std::vector<T> Minus(const std::vector<T> &a, const std::vector<T> &b) {
  std::vector<T> z(a.size(), 0.0);
  for (std::size_t k = 0; k < a.size(); ++k) { z[k] = a[k] - b[k]; }
  return z;
}

// Add two vectors
template <typename T>
inline std::vector<T> Plus(const std::vector<T> &a, const std::vector<T> &b) {
  std::vector<T> z(a.size(), 0.0);
  for (std::size_t k = 0; k < a.size(); ++k) { z[k] = a[k] + b[k]; }
  return z;
}

// Outer (tensor) product of two vectors: u (x) v = u v^T
template <typename T>
inline MMatrix<T> OuterProd(const std::vector<T> &u, const std::vector<T> &v) {
  const unsigned int m = u.size();
  const unsigned int n = v.size();
  MMatrix<T>         out(m, n);

  for (std::size_t i = 0; i < m; ++i) {
    for (std::size_t j = 0; j < n; ++j) { out[i][j] = u[i] * v[j]; }
  }
  return out;
}

// Kronecker tensor product of two matrices: A(x)B = C
//
// Dyadic product, e.g. 2 x R^3 vectors span 6-D space
// cf. full R^3x3 matrix has 9 components (like product state vs entanglement)
//
template <typename T>
inline MMatrix<T> TensorProd(const MMatrix<T> &A, const MMatrix<T> &B) {
  const unsigned int rowA = A.size_row();
  const unsigned int rowB = B.size_row();
  const unsigned int colA = A.size_col();
  const unsigned int colB = B.size_col();

  MMatrix<T> C(rowA * rowB, colA * colB);

  for (std::size_t i = 0; i < rowA; ++i) {
    for (std::size_t j = 0; j < colA; ++j) {
      for (std::size_t k = 0; k < rowB; ++k) {
        for (std::size_t l = 0; l < colB; ++l) {
          C[i * rowB + k][j * colB + l] = A[i][j] * B[k][l];
        }
      }
    }
  }
  return C;
}

// Print out vector elements
template <typename T>
inline void PrintVector(const std::vector<T> &A, const std::string name = "") {
  // Print elements
  std::cout << name << " : ";
  for (std::size_t i = 0; i < A.size(); ++i) {
    printf("%6.3f+i%6.3f | ", std::real(A[i]), std::imag(A[i]));
  }
  std::cout << std::endl;
}

// Print out matrix elements
template <typename T>
inline void PrintMatrix(const MMatrix<T> &A, const std::string name = "") {
  // Print elements
  std::cout << name << " : ";
  for (std::size_t i = 0; i < A.size_row(); ++i) {
    for (std::size_t j = 0; j < A.size_col(); ++j) {
      const std::string delim = (j < A.size_col() - 1) ? ", " : "";
      printf("%6.3f+i%6.3f%s ", std::real(A[i][j]), std::imag(A[i][j]), delim.c_str());
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

// Print out matrix elements as two separate matrix
template <typename T>
inline void PrintMatrixSeparate(const MMatrix<T> &A) {
  // Print real part
  std::cout << "Re:" << std::endl;
  for (std::size_t i = 0; i < A.size_row(); ++i) {
    for (std::size_t j = 0; j < A.size_col(); ++j) {
      const std::string delim = (j < A.size_col() - 1) ? ", " : "";
      printf("%6.3f%s ", std::real(A[i][j]), delim.c_str());
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;

  // Print imaginary part
  std::cout << "Im:" << std::endl;
  for (std::size_t i = 0; i < A.size_row(); ++i) {
    for (std::size_t j = 0; j < A.size_col(); ++j) {
      const std::string delim = (j < A.size_col() - 1) ? ", " : "";
      printf("%6.3f%s ", std::imag(A[i][j]), delim.c_str());
    }
    std::cout << std::endl;
  }
  std::cout << std::endl;
}

}  // matoper namespace ends
}  // gra namespace ends

#endif