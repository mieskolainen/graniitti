// Finite spin representations and projection bases
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MSPIN_H
#define MSPIN_H

#include <cstddef>
#include <complex>
#include <string>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMatrix.h"

namespace gra::spin {

// Compute R_z(phi) R_y(theta) in the negative-first spin-half basis
MMatrix<std::complex<double>> SpinHalfRotation(double theta, double phi);

// Compute the positive SL(2,C) boost from rest to a timelike momentum
MMatrix<std::complex<double>> SpinHalfBoost(const M4Vec &p);

// Compute the canonical Wigner rotation from an SL(2,C) frame map and a local momentum
MMatrix<std::complex<double>> SpinHalfWigner(const MMatrix<std::complex<double>> &frame, const M4Vec &p);

// Compute the symmetric spin-J representation of a unitary spin-half map
MMatrix<std::complex<double>> SpinRotation(const MMatrix<std::complex<double>> &u, double J);

// Compute the spin rotation from the rotated Cartesian x and z unit axes
// The Euler angles choose the SU(2) lift of the SO(3) frame
MMatrix<std::complex<double>> SpinFrame(const M4Vec &x, const M4Vec &z, double J);

// Store one finite irreducible SU(2) representation using doubled spin labels
class SpinRep {
 public:
  // Construct one representation from a nonnegative doubled spin
  explicit SpinRep(int spin_x2, const std::string &context);

  // Construct one representation from an integer or half-integer spin
  static SpinRep FromSpin(double spin, const std::string &context);

  // Compute the doubled spin label
  int X2() const noexcept;

  // Compute the physical spin
  double Spin() const noexcept;

  // Compute the representation dimension 2J+1
  std::size_t Dim() const noexcept;

  // Evaluate the spin representation of one generated unitary spin-half map
  MMatrix<std::complex<double>> Rotation(const MMatrix<std::complex<double>> &u) const;

  // Compute the ordered doubled projections -2J,-2J+2,...,+2J
  std::vector<int> ProjectionsX2() const;

  // Compute the ordered physical projections -J,-J+1,...,+J
  std::vector<double> Projections() const;

  // Compute the basis index of one doubled projection
  std::size_t IndexX2(int projection_x2, const std::string &context) const;

  // Compute the basis index of one integer or half-integer projection
  std::size_t Index(double projection, const std::string &context) const;

 private:
  int spin_x2 = 0;
};

// Select the physical projections of one finite spin representation
class SpinBasis {
 public:
  // Construct the complete projection basis of one representation
  explicit SpinBasis(SpinRep rep);

  // Construct a physical subspace from doubled projection labels
  SpinBasis(SpinRep rep, std::vector<int> projections_x2,
            const std::string &context);

  // Compute the underlying finite spin representation
  const SpinRep &Rep() const noexcept;

  // Compute the active doubled projection labels
  const std::vector<int> &ProjectionsX2() const noexcept;

  // Compute whether one doubled projection is active
  bool ContainsX2(int projection_x2) const noexcept;

 private:
  SpinRep rep;
  std::vector<int> projections_x2;
};

// Store the spin spaces of one parent and two ordered legs
class TwoBodySpin {
 public:
  // Construct one ordered two-body spin space
  TwoBodySpin(SpinRep parent, SpinBasis leg1, SpinBasis leg2);

  // Compute whether L and S can couple the three representations
  bool AllowsLS(std::size_t l, int two_s) const noexcept;

  // Compute one Jacob-Wick LS operator in the full leg basis
  MMatrix<std::complex<double>> JacobWick(std::size_t l, int two_s,
                                          const std::string &context) const;

 private:
  SpinRep parent;
  SpinBasis leg1;
  SpinBasis leg2;
};

// Compute doubled spin after validating an integer or half-integer label
int SpinLabelX2(double spin, const std::string &context);

// Compute the dimension of one finite spin representation
std::size_t SpinStateCount(double spin, const std::string &context);

// Compute one spin-projection basis index
std::size_t SpinProjectionIndex(double projection, double spin,
                                const std::string &context);

// Compute the ordered physical projections of one spin representation
std::vector<double> SpinProjections(double spin);

}  // namespace gra::spin

#endif
