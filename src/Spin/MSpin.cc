// Finite spin representations and projection bases
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Spin/MSpin.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MException.h"

namespace gra::spin {
namespace {

constexpr double kSpinTolerance = 1e-9;

// Compute whether three doubled spins satisfy the SU(2) triangle rules
bool TriangleX2(int a, int b, int c) noexcept {
  return a >= 0 && b >= 0 && c >= 0 && c <= a + b && c >= std::abs(a - b) && (a + b + c) % 2 == 0;
}

// Compute one finite doubled quantum number from a floating label
int DoubledLabel(double value, const std::string &context) {
  if (!std::isfinite(value)) { throw std::invalid_argument(context + ": quantum number is not finite"); }
  const double doubled = 2.0 * value;
  const double rounded = std::round(doubled);
  if (std::abs(doubled - rounded) > kSpinTolerance ||
      std::abs(rounded) > static_cast<double>(std::numeric_limits<int>::max())) {
    throw std::invalid_argument(context + ": quantum number is not integer or half-integer");
  }
  return static_cast<int>(std::llround(rounded));
}

}  // namespace

// Compute R_z(phi) R_y(theta) without losing the sign of a 2pi rotation
MMatrix<std::complex<double>> SpinHalfRotation(double theta, double phi) {
  if (!std::isfinite(theta) || !std::isfinite(phi)) {
    throw AmplitudeFailure("SpinHalfRotation: non-finite generated angle");
  }
  const double c = std::cos(theta / 2.0);
  const double s = std::sin(theta / 2.0);
  const auto phase = std::exp(std::complex<double>(0.0, phi / 2.0));
  return {{phase * c, phase * s}, {-std::conj(phase) * s, std::conj(phase) * c}};
}

// Compute B(p) = ((E+m)I + sigma.p)/sqrt(2m(E+m))
MMatrix<std::complex<double>> SpinHalfBoost(const M4Vec &p) {
  if (!(p.E() > 0.0) || !std::isfinite(p.E())) {
    throw AmplitudeFailure("SpinHalfBoost: momentum is not future timelike");
  }
  // Scale out E before forming the invariant or the boost normalization
  const M4Vec q(p.Px() / p.E(), p.Py() / p.E(), p.Pz() / p.E(), 1.0);
  const double mass2 = q.M2();
  if (!(mass2 > 0.0) || !std::isfinite(mass2)) {
    throw AmplitudeFailure("SpinHalfBoost: momentum is not future timelike");
  }
  const double mass = std::sqrt(mass2);
  const double e = 1.0 + mass;
  const std::complex<double> transverse(q.Px(), q.Py());
  MMatrix<std::complex<double>> boost = {{e - q.Pz(), transverse}, {std::conj(transverse), e + q.Pz()}};
  boost *= 1.0 / std::sqrt(2.0 * mass * e);
  return boost;
}

// Extract the SU(2) factor of frame B(p) without subtracting a large inverse boost
// For L in SL(2,C), the unitary polar factor is proportional to L + adj(L)^dagger
MMatrix<std::complex<double>> SpinHalfWigner(const MMatrix<std::complex<double>> &frame, const M4Vec &p) {
  if (frame.size_row() != 2 || frame.size_col() != 2) {
    throw std::invalid_argument("SpinHalfWigner: expected an SL(2,C) frame");
  }
  const auto lorentz = frame * SpinHalfBoost(p);
  const auto a = lorentz[0][0] + std::conj(lorentz[1][1]);
  const auto b = lorentz[0][1] - std::conj(lorentz[1][0]);
  const double norm = std::hypot(std::abs(a), std::abs(b));
  if (!(norm > 0.0) || !std::isfinite(norm)) {
    throw AmplitudeFailure("SpinHalfWigner: undefined generated rotation");
  }
  return {{a / norm, b / norm}, {-std::conj(b) / norm, std::conj(a) / norm}};
}

// Couple symmetric spinors recursively with the normalized j x 1/2 CG coefficients
MMatrix<std::complex<double>> SpinRep::Rotation(const MMatrix<std::complex<double>> &u) const {
  const int two_j = spin_x2;
  if (u.size_row() != 2 || u.size_col() != 2) {
    throw std::invalid_argument("SpinRotation: expected a spin-half matrix");
  }
  if (!u.IsFinite()) { throw AmplitudeFailure("SpinRotation: non-finite generated spin map"); }
  const MMatrix<std::complex<double>> identity(2, 2, "eye");
  if ((u.Dagger() * u - identity).FrobNorm() > 1e-8) {
    throw AmplitudeFailure("SpinRotation: generated spin map is not unitary");
  }
  MMatrix<std::complex<double>> rotation(1, 1, 1.0);
  for (int n = 1; n <= two_j; ++n) {
    MMatrix<std::complex<double>> next(n + 1, n + 1, 0.0);
    for (int i = 0; i <= n; ++i) {
      for (int j = 0; j <= n; ++j) {
        if (i < n && j < n) { next[i][j] += std::sqrt(double(n - i) * (n - j)) / n * rotation[i][j] * u[0][0]; }
        if (i < n && j > 0) { next[i][j] += std::sqrt(double(n - i) * j) / n * rotation[i][j - 1] * u[0][1]; }
        if (i > 0 && j < n) { next[i][j] += std::sqrt(double(i) * (n - j)) / n * rotation[i - 1][j] * u[1][0]; }
        if (i > 0 && j > 0) { next[i][j] += std::sqrt(double(i) * j) / n * rotation[i - 1][j - 1] * u[1][1]; }
      }
    }
    rotation = std::move(next);
  }
  return rotation;
}

// Evaluate a standalone spin rotation after validating its representation
MMatrix<std::complex<double>> SpinRotation(const MMatrix<std::complex<double>> &u, double J) {
  return SpinRep::FromSpin(J, "SpinRotation").Rotation(u);
}

// Resolve the final Euler rotation about z from the transported x axis
MMatrix<std::complex<double>> SpinFrame(const M4Vec &x, const M4Vec &z, double J) {
  M4Vec transverse = x;
  transverse.RotateZ(-z.Phi());
  transverse.RotateY(-z.Theta());
  return SpinRotation(SpinHalfRotation(z.Theta(), z.Phi()) * SpinHalfRotation(0.0, transverse.Phi()), J);
}

// Construct one representation from a nonnegative doubled spin
SpinRep::SpinRep(int spin_x2, const std::string &context) : spin_x2(spin_x2) {
  if (spin_x2 < 0) { throw std::invalid_argument(context + ": negative doubled spin"); }
}

// Construct one representation from an integer or half-integer spin
SpinRep SpinRep::FromSpin(double spin, const std::string &context) {
  const int spin_x2 = DoubledLabel(spin, context);
  if (spin_x2 < 0) { throw std::invalid_argument(context + ": negative spin"); }
  return SpinRep(spin_x2, context);
}

// Compute the doubled spin label
int SpinRep::X2() const noexcept { return spin_x2; }

// Compute the physical spin
double SpinRep::Spin() const noexcept { return 0.5 * static_cast<double>(spin_x2); }

// Compute the representation dimension 2J+1
std::size_t SpinRep::Dim() const noexcept { return static_cast<std::size_t>(spin_x2) + 1; }

// Compute the ordered doubled projections -2J,-2J+2,...,+2J
std::vector<int> SpinRep::ProjectionsX2() const {
  std::vector<int> projections;
  projections.reserve(Dim());
  for (std::size_t i = 0; i < Dim(); ++i) {
    const long long projection_x2 = -static_cast<long long>(spin_x2) + 2 * static_cast<long long>(i);
    projections.push_back(static_cast<int>(projection_x2));
  }
  return projections;
}

// Compute the ordered physical projections -J,-J+1,...,+J
std::vector<double> SpinRep::Projections() const {
  std::vector<double> projections;
  projections.reserve(Dim());
  for (const int projection_x2 : ProjectionsX2()) { projections.push_back(0.5 * static_cast<double>(projection_x2)); }
  return projections;
}

// Compute the basis index of one doubled projection
std::size_t SpinRep::IndexX2(int projection_x2, const std::string &context) const {
  if (projection_x2 < -spin_x2 || projection_x2 > spin_x2 || (projection_x2 + spin_x2) % 2 != 0) {
    throw std::invalid_argument(context + ": spin projection outside spin basis");
  }
  return static_cast<std::size_t>((projection_x2 + spin_x2) / 2);
}

// Compute the basis index of one integer or half-integer projection
std::size_t SpinRep::Index(double projection, const std::string &context) const {
  return IndexX2(DoubledLabel(projection, context), context);
}

// Construct the complete projection basis of one representation
SpinBasis::SpinBasis(SpinRep rep) : rep(rep), projections_x2(this->rep.ProjectionsX2()) {}

// Construct a physical subspace from doubled projection labels
SpinBasis::SpinBasis(SpinRep rep, std::vector<int> projections_x2, const std::string &context)
    : rep(rep), projections_x2(std::move(projections_x2)) {
  if (this->projections_x2.empty()) { throw std::invalid_argument(context + ": empty spin basis"); }
  for (const int projection_x2 : this->projections_x2) { static_cast<void>(this->rep.IndexX2(projection_x2, context)); }
  std::sort(this->projections_x2.begin(), this->projections_x2.end());
  if (std::adjacent_find(this->projections_x2.begin(), this->projections_x2.end()) != this->projections_x2.end()) {
    throw std::invalid_argument(context + ": duplicate spin projection");
  }
}

// Compute the underlying finite spin representation
const SpinRep &SpinBasis::Rep() const noexcept { return rep; }

// Compute the active doubled projection labels
const std::vector<int> &SpinBasis::ProjectionsX2() const noexcept { return projections_x2; }

// Compute whether one doubled projection is active
bool SpinBasis::ContainsX2(int projection_x2) const noexcept {
  return std::binary_search(projections_x2.begin(), projections_x2.end(), projection_x2);
}

// Construct one ordered two-body spin space
TwoBodySpin::TwoBodySpin(SpinRep parent, SpinBasis leg1, SpinBasis leg2)
    : parent(parent), leg1(std::move(leg1)), leg2(std::move(leg2)) {}

// Compute whether L and S can couple the three representations
bool TwoBodySpin::AllowsLS(std::size_t l, int two_s) const noexcept {
  if (l > static_cast<std::size_t>(std::numeric_limits<int>::max() / 2)) { return false; }
  return TriangleX2(leg1.Rep().X2(), leg2.Rep().X2(), two_s) && TriangleX2(2 * static_cast<int>(l), two_s, parent.X2());
}

// Compute one Jacob-Wick LS operator in the full leg basis
MMatrix<std::complex<double>> TwoBodySpin::JacobWick(std::size_t l, int two_s, const std::string &context) const {
  if (!AllowsLS(l, two_s)) { throw std::invalid_argument(context + ": forbidden LS coupling"); }
  MMatrix<std::complex<double>> op(leg1.Rep().Dim(), leg2.Rep().Dim(), 0.0);
  const double                  J    = parent.Spin();
  const double                  s1   = leg1.Rep().Spin();
  const double                  s2   = leg2.Rep().Spin();
  const double                  s    = 0.5 * static_cast<double>(two_s);
  const double                  norm = std::sqrt((2.0 * static_cast<double>(l) + 1.0) / (2.0 * J + 1.0));

  for (const int m1x2 : leg1.Rep().ProjectionsX2()) {
    if (!leg1.ContainsX2(m1x2)) { continue; }
    const std::size_t i  = leg1.Rep().IndexX2(m1x2, context + " leg 1");
    const double      m1 = 0.5 * static_cast<double>(m1x2);
    for (const int m2x2 : leg2.Rep().ProjectionsX2()) {
      if (!leg2.ContainsX2(m2x2)) { continue; }
      const std::size_t j      = leg2.Rep().IndexX2(m2x2, context + " leg 2");
      const double      m2     = 0.5 * static_cast<double>(m2x2);
      const double      lambda = m1 - m2;
      if (std::abs(lambda) > J + kSpinTolerance) { continue; }
      op[i][j] =
          norm * wigner::CG(static_cast<double>(l), s, 0.0, lambda, J, lambda) * wigner::CG(s1, s2, m1, -m2, s, lambda);
    }
  }
  return op;
}

// Compute doubled spin after validating an integer or half-integer label
int SpinLabelX2(double spin, const std::string &context) { return SpinRep::FromSpin(spin, context).X2(); }

// Compute the dimension of one finite spin representation
std::size_t SpinStateCount(double spin, const std::string &context) { return SpinRep::FromSpin(spin, context).Dim(); }

// Compute one spin-projection basis index
std::size_t SpinProjectionIndex(double projection, double spin, const std::string &context) {
  return SpinRep::FromSpin(spin, context + " spin").Index(projection, context);
}

// Compute the ordered physical projections of one spin representation
std::vector<double> SpinProjections(double spin) {
  return SpinRep::FromSpin(spin, "gra::spin::SpinProjections").Projections();
}

}  // namespace gra::spin
