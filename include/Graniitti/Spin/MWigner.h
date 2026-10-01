// Wigner rotations and SU(2) angular momentum algebra
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MWIGNER_H
#define MWIGNER_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Spin/MSpin.h"

namespace gra {
namespace wigner {

// Store one Wigner half-angle polynomial and its factorial normalization
struct Polynomial {
  double factor = 1.0;
  std::vector<double> coefficient;
};

// Store fixed Wigner coefficients for repeated angular evaluations
class Rotation {
 public:
  // Validate spin labels and prepare the half-angle polynomial coefficients
  Rotation(const std::vector<double> &m, const std::vector<double> &mp, double J);
  // Evaluate the prepared rows at one generated polar angle
  MMatrix<double> Evaluate(double theta) const;

 private:
  int spin_x2 = 0;
  spin::SpinRep representation;
  std::vector<double> rows, columns;
  std::vector<std::size_t> row_indices, column_indices;
  std::vector<Polynomial> coefficients;
};

// Compute whether x is an integer within the spin-label tolerance
bool IsInt(double x);

// Compute one SU(2) Clebsch-Gordan coefficient
double CG(double j1, double j2, double m1, double m2, double j, double m);

// Compute one Wigner 3j symbol
double W3j(double j1, double j2, double j3, double m1, double m2, double m3);

// Compute one Gamma continued Wigner 3j function
std::complex<double> W3jRegge(double j1, double j2, double j3, double m1,
                              double m2, double m3);

// Compute one Gamma continued Clebsch-Gordan function
std::complex<double> CGRegge(double j1, double j2, double m1, double m2,
                             double j, double m);

// Compute one real Wigner element d^J_(mp,m) with arguments ordered as m,mp
double d(double theta, double m, double mp, double J);

// Compute one Gamma continued Wigner d function
std::complex<double> dRegge(double theta, double m, double mp, double J);

// Compute real Wigner rows for multiple m and mp labels at one polar angle
MMatrix<double> dRows(double theta, const std::vector<double> &m_values,
                      const std::vector<double> &mp_values, double J);

// Compute D^(J*)_(mp,m)(phi,theta,0) with arguments ordered as m,mp
std::complex<double> D(double theta, double phi, double m, double mp, double J);

// Compute one Gamma continued conjugate Wigner D function
std::complex<double> DRegge(double theta, double phi, double m, double mp,
                            double J);

// Compute the stored transpose matrix with rows m and columns mp
MMatrix<std::complex<double>> DMatrix(double J, double theta, double phi);

} // namespace wigner
} // namespace gra

#endif
