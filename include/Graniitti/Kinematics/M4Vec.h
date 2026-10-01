// 4-vectors [HEADER ONLY class]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Standard metric (+,-,-,-) and MC initialization convention (px,py,pz,e)
//
// Operators ^ and % use Lorentz ordering (E,px,py,pz)
// Operator ^ gives contravariant components and % gives covariant components
//
//  p^0 =  E
//  p^1 =  px
//  p^2 =  py
//  p^3 =  pz
//
//  p%0 =  E
//  p%1 = -px
//  p%2 = -py
//  p%3 = -pz

#ifndef M4VEC_H
#define M4VEC_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <numbers>
#include <stdexcept>

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMatrix.h"

namespace gra {

using M3Vec = std::array<double, 3>;

// ----------------------------------------------------------------------
// C++17 allows class template argument default deduction without brackets
//
// M4Vec a;
// M4Vec<> b;
// M4Vec<double> c;
//
// a,b,c are all valid
//
// template <typename T = double>
//
// Complex 4-vectors could be implemented this way, TBD
// ----------------------------------------------------------------------

class M4Vec {
public:
  // All default to zero
  M4Vec() { k = {0.0, 0.0, 0.0, 0.0}; }

  // Initialize (note the internal order)
  M4Vec(double x, double y, double z, double t) { k = {t, x, y, z}; }

  // These work fine with std::array
  M4Vec(const M4Vec &) = default;
  M4Vec &operator=(const M4Vec &) = default;

  // SET methods
  void SetPx(double v) { k[X_] = v; }
  void SetPy(double v) { k[Y_] = v; }
  void SetPz(double v) { k[Z_] = v; }
  void SetE(double v) { k[E_] = v; }

  void SetX(double v) { k[X_] = v; }
  void SetY(double v) { k[Y_] = v; }
  void SetZ(double v) { k[Z_] = v; }
  void SetT(double v) { k[E_] = v; }

  void Set(double x, double y, double z, double t) {
    k[X_] = x;
    k[Y_] = y;
    k[Z_] = z;
    k[E_] = t;
  }

  // Set an on-shell four-momentum with E = sqrt(px^2+py^2+pz^2+m^2)
  void SetPxPyPzM(double x, double y, double z, double m) {
    k[X_] = x;
    k[Y_] = y;
    k[Z_] = z;
    k[E_] = std::hypot(P3mod(), m);
  }

  void SetPxPyPz(double x, double y, double z) {
    k[X_] = x;
    k[Y_] = y;
    k[Z_] = z;
  }

  void SetPxPy(double x, double y) {
    k[X_] = x;
    k[Y_] = y;
  }

  void SetPzE(double z, double e) {
    k[Z_] = z;
    k[E_] = e;
  }

  void SetP3(const M3Vec &vec) {
    k[X_] = vec[0];
    k[Y_] = vec[1];
    k[Z_] = vec[2];
  }

  // Apply metric tensor: eta_{\mu\nu} k^\nu = k_\mu
  void Flip3() {
    k[X_] = -k[X_];
    k[Y_] = -k[Y_];
    k[Z_] = -k[Z_];
  }

  // GET methods
  double Px() const { return k[X_]; }
  double Py() const { return k[Y_]; }
  double Pz() const { return k[Z_]; }
  double E() const { return k[E_]; }

  double X() const { return k[X_]; }
  double Y() const { return k[Y_]; }
  double Z() const { return k[Z_]; }
  double T() const { return k[E_]; }

  // Compute contravariant components in (E, px, py, pz) order
  template <typename T = double> std::array<T, 4> Contravariant() const {
    return {static_cast<T>(E()), static_cast<T>(Px()), static_cast<T>(Py()),
            static_cast<T>(Pz())};
  }

  // Compute a boosted four-vector for a finite velocity with beta^2 < 1
  // E' = gamma(E + beta.p), p' = p + [(gamma-1)(beta.p)/beta^2 + gamma E] beta
  M4Vec LorentzBoost(const M3Vec &betavec, int sign = 1) const {
    // Apply Lorentz boost by beta vector (boost direction 3-momentum divided by
    // energy) sign (+1,-1) allows to flip the boost direction
    //
    if (sign != -1 && sign != 1) { return M4Vec(0.0, 0.0, 0.0, -1.0); }
    const long double bx = sign * static_cast<long double>(betavec[0]);
    const long double by = sign * static_cast<long double>(betavec[1]);
    const long double bz = sign * static_cast<long double>(betavec[2]);
    const long double b2 = bx * bx + by * by + bz * bz;
    if (!std::isfinite(b2) || b2 >= 1.0L || !gra::AllFinite(k)) {
      return M4Vec(0.0, 0.0, 0.0, -1.0);
    }

    const long double dot = bx * k[X_] + by * k[Y_] + bz * k[Z_];
    const long double gamma = 1.0L / std::sqrt(1.0L - b2);
    // (gamma-1)/beta^2 = gamma^2/(gamma+1), including the limit at rest
    const long double gamma2 = gamma * gamma / (gamma + 1.0L);
    const long double scale = gamma2 * dot + gamma * k[E_];
    const M4Vec result(static_cast<double>(k[X_] + scale * bx),
                       static_cast<double>(k[Y_] + scale * by),
                       static_cast<double>(k[Z_] + scale * bz),
                       static_cast<double>(gamma * (k[E_] + dot)));
    return gra::AllFinite(result.k) ? result : M4Vec(0.0, 0.0, 0.0, -1.0);
  }

  // Rotate the spatial components with the stated theta and phi convention
  // R(theta,phi) = R_z(phi) R_y(theta)
  void Rotate(double theta, double phi) {
    const double c1 = std::cos(theta);
    const double s1 = std::sin(theta);
    const double c2 = std::cos(phi);
    const double s2 = std::sin(phi);

    const MMatrix<double> rotation = {
        {c1 * c2, -s2, s1 * c2}, {c1 * s2, c2, s1 * s2}, {-s1, 0.0, c1}};
    Rotate(rotation);
  }

  // Rotate the spatial components with one 3 by 3 rotation matrix
  void Rotate(const MMatrix<double> &rotation) { SetP3(rotation * P3()); }

  // Apply an active counterclockwise rotation around the x axis
  void RotateX(double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double py = cosine * Py() - sine * Pz();
    const double pz = sine * Py() + cosine * Pz();
    SetPy(py);
    SetPz(pz);
  }

  // Apply an active counterclockwise rotation around the y axis
  void RotateY(double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double px = cosine * Px() + sine * Pz();
    const double pz = -sine * Px() + cosine * Pz();
    SetPx(px);
    SetPz(pz);
  }

  // Apply an active counterclockwise rotation around the z axis
  void RotateZ(double angle) {
    const double cosine = std::cos(angle);
    const double sine = std::sin(angle);
    const double px = cosine * Px() - sine * Py();
    const double py = sine * Px() + cosine * Py();
    SetPx(px);
    SetPy(py);
  }

  // Construct the spatial rotation from this direction to the target direction
  // Rodrigues rotation with a perpendicular axis for exactly opposite directions
  MMatrix<double> RotationTo(const M3Vec &target) const {
    const MMatrix<double> identity(3, 3, "eye");
    const M3Vec source = P3();
    const double source_norm = std::hypot(source[0], source[1], source[2]);
    const double target_norm = std::hypot(target[0], target[1], target[2]);
    if (!(source_norm > 0.0) || !(target_norm > 0.0)) {
      return identity;
    }

    const M3Vec source_unit = {source[0] / source_norm, source[1] / source_norm,
                               source[2] / source_norm};
    const M3Vec target_unit = {target[0] / target_norm, target[1] / target_norm,
                               target[2] / target_norm};
    const double cosine =
        std::clamp(gra::InnerProduct(source_unit, target_unit), -1.0, 1.0);

    M3Vec axis = gra::CrossProduct(source_unit, target_unit);
    const double sine = std::hypot(axis[0], axis[1], axis[2]);
    if (!(sine > 0.0) || (cosine < 0.0 && sine < 8.0 * std::numeric_limits<double>::epsilon())) {
      if (cosine >= 0.0) { return identity; }
      M3Vec reference = {1.0, 0.0, 0.0};
      if (std::abs(source_unit[1]) <= std::abs(source_unit[0]) &&
          std::abs(source_unit[1]) <= std::abs(source_unit[2])) {
        reference = {0.0, 1.0, 0.0};
      } else if (std::abs(source_unit[2]) <= std::abs(source_unit[0])) {
        reference = {0.0, 0.0, 1.0};
      }

      const M3Vec axis =
          gra::NormalizedL2(gra::CrossProduct(source_unit, reference));
      const MMatrix<double> projector = {
          {axis[0] * axis[0], axis[0] * axis[1], axis[0] * axis[2]},
          {axis[1] * axis[0], axis[1] * axis[1], axis[1] * axis[2]},
          {axis[2] * axis[0], axis[2] * axis[1], axis[2] * axis[2]}};
      return projector * 2.0 - identity;
    }

    for (auto &component : axis) { component /= sine; }
    // Remove rounding along the source direction, which matters near pi
    gra::AddScaled(axis, source_unit, -gra::InnerProduct(axis, source_unit));
    axis = gra::NormalizedL2(axis);
    const MMatrix<double> cross = {{0.0, -axis[2], axis[1]},
                                   {axis[2], 0.0, -axis[0]},
                                   {-axis[1], axis[0], 0.0}};
    return identity + cross * sine + (cross * cross) * (1.0 - cosine);
  }

  // Rotate this spatial direction onto the target while preserving its norm
  void RotateTo(const M3Vec &target) { Rotate(RotationTo(target)); }

  // Compute 3-vector
  M3Vec P3() const { return {k[X_], k[Y_], k[Z_]}; }

  // Compute Lorentz boost 3-vector, i.e. 3-momentum divided by energy
  // beta = p/E
  M3Vec BetaVector() const {
    if (math::IsZero(k[E_])) {
      return {0.0, 0.0, 0.0};
    }
    return {k[X_] / k[E_], k[Y_] / k[E_], k[Z_] / k[E_]};
  }

  // ALGEBRA methods

  // Space-time invariants
  // M^2 = E^2 - |p|^2
  template <typename T = double>
  T Invariant() const {
    const auto p = Contravariant<T>();
    return MinkowskiProduct(p, p);
  }
  double M2() const { return Invariant(); }
  double M() const { return (M2() > 0.0) ? msqrt(M2()) : -msqrt(-M2()); }

  // Compute a physical invariant mass, tolerating null roundoff and rejecting spacelike momenta
  bool PhysicalMass(double &mass) const {
    const auto p = Contravariant<long double>();
    const long double spatial2 = p[1] * p[1] + p[2] * p[2] + p[3] * p[3];
    long double mass2 = MinkowskiProduct(p, p);
    const long double scale = std::max(p[0] * p[0], spatial2);
    const long double tolerance = 64.0L * std::numeric_limits<double>::epsilon() * scale;
    if (!std::isfinite(mass2) || mass2 < -tolerance) {
      mass = 0.0;
      return false;
    }
    if (std::abs(mass2) <= tolerance) { mass2 = 0.0; }
    mass = static_cast<double>(std::sqrt(mass2));
    return std::isfinite(mass);
  }

  // gamma = E/m = 1/sqrt(1-v^2/c^2) = 1/sqrt(1-beta^2)
  double Gamma() const { return (M() > 0.0) ? E() / M() : -1.0; }
  double Beta() const {
    return !math::IsZero(E()) ? P3mod() / E() : -1.0;
  }

  // Transverse 2-vector norm and norm^2
  // pT^2 = px^2 + py^2
  double Perp() const { return Pt(); }
  double Perp2() const { return Pt2(); }
  double Pt() const { return std::hypot(Px(), Py()); }
  double Pt2() const { return k[X_] * k[X_] + k[Y_] * k[Y_]; }

  // Total 3-vector norm and norm^2
  // |p|^2 = px^2 + py^2 + pz^2
  double P3mod() const { return std::hypot(Px(), Py(), Pz()); }
  double P3mod2() const {
    return k[X_] * k[X_] + k[Y_] * k[Y_] + k[Z_] * k[Z_];
  }

  // Transverse mass (invariant under boost in z-direction)
  // mT^2 = M^2 + pT^2 = E^2 - pz^2
  double Mt() const { return msqrt(Mt2()); }
  double Mt2() const {
    const long double energy = E(), pz = Pz();
    return static_cast<double>((energy - pz) * (energy + pz));
  }

  // Coincides with transverse mass for single particle
  double Et() const { return Mt(); }
  double Et2() const { return Mt2(); }

  // Angles
  // phi = atan2(py,px), theta = atan2(pT,pz)
  double Phi() const {
    return (math::IsZero(Px()) && math::IsZero(Py()))
               ? 0.0
               : std::atan2(Py(), Px());
  } // y / x, range [-PI,PI]
  double Theta() const {
    return (math::IsZero(Px()) && math::IsZero(Py()) &&
            math::IsZero(Pz()))
               ? 0.0
               : std::atan2(Pt(), Pz());
  } // |Pt| / z
  double CosTheta() const { return std::cos(Theta()); }

  // Pseudorapidity and rapidity (boost) in z-direction
  // eta = (1/2)log[(|p|+pz)/(|p|-pz)], y = (1/2)log[(E+pz)/(E-pz)]
  double Eta() const {
    return static_cast<double>(std::asinh(static_cast<long double>(Pz()) / Pt()));
  }
  double Rap() const {
    if (std::abs(Pz()) < 0.5 * std::abs(E())) { return std::atanh(Pz() / E()); }
    const long double energy = E(), pz = Pz();
    return static_cast<double>(0.5L * std::log((energy + pz) / (energy - pz)));
  }

  // Lightcone variable: k_+ = E + p_z
  double LightconePos() const { return E() + Pz(); }

  // Lightcone variable: k_- = E - p_z
  double LightconeNeg() const { return E() - Pz(); }

  // SPINOR-HELICITY CO-VARIABLES

  // pT^C = px + i py

  std::complex<double> ComplexPt() const {
    return Px() + std::complex<double>(0, 1) * Py();
  }
  // exp(i phi) pT/mT = (px+i py)/sqrt[(E+pz)(E-pz)]
  std::complex<double> ExpCPhi() const {
    return ComplexPt() / msqrt(LightconePos() * LightconeNeg());
  }

  // 2-BODY ALGEBRA

  // Minkowski 4-product
  // p.q = p^0 q^0 - p_vec.q_vec
  double DotM(const M4Vec &rhs) const {
    return E() * rhs.E() -
           (Px() * rhs.Px() + Py() * rhs.Py() + Pz() * rhs.Pz());
  }

  // 3-vector dot product
  // p_vec.q_vec = px qx + py qy + pz qz
  double Dot3(const M4Vec &rhs) const {
    return Px() * rhs.Px() + Py() * rhs.Py() + Pz() * rhs.Pz();
  }

  // Transverse 2-vector dot product
  // pT.qT = px qx + py qy
  double DotPt(const M4Vec &rhs) const {
    return Px() * rhs.Px() + Py() * rhs.Py();
  }

  // 3-vector cross product (return vector with 0 energy/time)
  // (p cross q)_i = epsilon_{ijk} p_j q_k
  M4Vec Cross3(const M4Vec &rhs) const {
    M4Vec a(Py() * rhs.Pz() - Pz() * rhs.Py(),
            Pz() * rhs.Px() - Px() * rhs.Pz(),
            Px() * rhs.Py() - Py() * rhs.Px(), 0.0);
    return a;
  }

  // Azimuth angle difference between [-PI,PI]
  double DeltaPhi(const M4Vec &v) const {
    double D = Phi() - v.Phi();
    while (D >= std::numbers::pi) {
      D -= 2.0 * std::numbers::pi;
    }
    while (D < -std::numbers::pi) {
      D += 2.0 * std::numbers::pi;
    }
    return D;
  }

  // Azimuth angle between [0,PI]
  double DeltaPhiAbs(const M4Vec &v) const { return std::abs(DeltaPhi(v)); }

  // OPERATORS
  double operator[](size_t mu) const {
    ValidateIndex(mu);
    return k[mu];
  }
  double &operator[](size_t mu) {
    ValidateIndex(mu);
    return k[mu];
  }

  // Access operator in normal contravariant (upper index) indexing
  double operator^(size_t mu) const {
    ValidateIndex(mu);
    return k[mu];
  }

  // Access operator with simultaneous lowering with metric (covariant index)
  double operator%(size_t mu) const {
    ValidateIndex(mu);
    if (mu == E_) {
      return k[E_];
    } else {
      return -k[mu];
    }
  }

  // Minkowski scalar product
  double operator*(const M4Vec &rhs) const { return DotM(rhs); }

  // 4-vector + 4-vector
  M4Vec operator+(const M4Vec &rhs) const {
    return M4Vec(k[X_] + rhs.k[X_], k[Y_] + rhs.k[Y_], k[Z_] + rhs.k[Z_],
                 k[E_] + rhs.k[E_]);
  }
  M4Vec operator-(const M4Vec &rhs) const {
    return M4Vec(k[X_] - rhs.k[X_], k[Y_] - rhs.k[Y_], k[Z_] - rhs.k[Z_],
                 k[E_] - rhs.k[E_]);
  }

  // Flip sign of all components
  M4Vec operator-() const { return M4Vec(-k[X_], -k[Y_], -k[Z_], -k[E_]); }

  // 4-vector */ scalar
  M4Vec operator*(const double rhs) const {
    return M4Vec(k[X_] * rhs, k[Y_] * rhs, k[Z_] * rhs, k[E_] * rhs);
  }
  M4Vec operator/(const double rhs) const {
    return M4Vec(k[X_] / rhs, k[Y_] / rhs, k[Z_] / rhs, k[E_] / rhs);
  }

  // Comparison
  bool operator==(const M4Vec &rhs) const {
    const double EPS = 1e-10;
    return std::abs(k[E_] - rhs.k[E_]) < EPS &&
           std::abs(k[X_] - rhs.k[X_]) < EPS &&
           std::abs(k[Y_] - rhs.k[Y_]) < EPS &&
           std::abs(k[Z_] - rhs.k[Z_]) < EPS;
  }
  bool operator!=(const M4Vec &rhs) const { return !(*this == rhs); }

  // 4-vector +-= 4-vector
  void operator+=(const M4Vec &rhs) {
    k[E_] += rhs.k[E_];
    k[X_] += rhs.k[X_];
    k[Y_] += rhs.k[Y_];
    k[Z_] += rhs.k[Z_];
  }
  void operator-=(const M4Vec &rhs) {
    k[E_] -= rhs.k[E_];
    k[X_] -= rhs.k[X_];
    k[Y_] -= rhs.k[Y_];
    k[Z_] -= rhs.k[Z_];
  }

  // 4-vector */= scalar
  void operator*=(const double rhs) {
    k[E_] *= rhs;
    k[X_] *= rhs;
    k[Y_] *= rhs;
    k[Z_] *= rhs;
  }
  void operator/=(const double rhs) {
    k[E_] /= rhs;
    k[X_] /= rhs;
    k[Y_] /= rhs;
    k[Z_] /= rhs;
  }

  void Print(const std::string name = "") const {
    std::cout << "M4Vec::" << name << " Px (X): " << Px()
              << ", Py (Y): " << Py() << ", Pz (Z): " << Pz()
              << ", E (T): " << E() << ", M (S): " << M()
              << ", theta: " << Theta() << ", phi: " << Phi() << std::endl;
  }

  // Overload the << operator for M4Vec
  friend std::ostream &operator<<(std::ostream &os, const M4Vec &vec) {
    os << "M4Vec: Px (X): " << vec.Px() << ", Py (Y): " << vec.Py()
       << ", Pz (Z): " << vec.Pz() << ", E (T): " << vec.E()
       << ", M (S): " << vec.M() << ", theta: " << vec.Theta()
       << ", phi: " << vec.Phi();
    return os;
  }

  // Particle 4-position starting starting propagation from (0,0,0,0)
  // x^mu = (gamma tau_0 c, gamma tau_0 beta_vec c) times scale
  // - p is the particle 4-momentum in the lab frame
  // - tau0 is the particle flight time in its rest frame
  M4Vec PropagatePosition(double tau0, double scale) const {
    // Flight time in the lab frame
    const double gamma = Gamma();
    const double tau = gamma * tau0;

    // Velocity
    const double beta = Beta();

    // End point 4-position
    if (math::IsZero(P3mod2())) {
      return M4Vec(0.0, 0.0, 0.0, tau * c * scale);
    }

    return M4Vec(tau * c * beta * Px() / P3mod() * scale,
                 tau * c * beta * Py() / P3mod() * scale,
                 tau * c * beta * Pz() / P3mod() * scale, tau * c * scale);
  }

private:
  // speed of light, c = [m/s] (EXACT/DEFINITION)
  static constexpr const double c = 2.99792458E8;

  // Indices
  static const int E_ = 0;
  static const int X_ = 1;
  static const int Y_ = 2;
  static const int Z_ = 3;

  // Safe sqrt
  double msqrt(double x) const { return std::sqrt(std::max(x, 0.0)); }
  void ValidateIndex(size_t mu) const {
    if (mu > 3) {
      throw std::out_of_range("M4Vec: index out of bounds");
    }
  }

  // 4-vector (faster than std::vector)
  std::array<double, 4> k{};
};

} // namespace gra

#endif
