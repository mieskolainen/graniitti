// Dirac spinors, gamma algebra and polarization states
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

// Tensor algebra
#include "FTensor.hpp"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

using FTensor::Tensor1;
using FTensor::Tensor2;

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;

namespace gra {

namespace {

// Compute a nonnegative fermion mass while tolerating mass-shell roundoff
// m = sqrt(max(0,E^2-|p|^2))
double FermionMass(const M4Vec &p, const std::string &context) {
  const double E  = p.E();
  const double p3 = p.P3mod();
  if (!std::isfinite(E) || !std::isfinite(p3)) { throw AmplitudeFailure(context + ": non-finite momentum"); }

  const double mass2     = std::fma(E, E, -p3 * p3);
  const double scale     = E * E + p3 * p3 + 1.0;
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  if (mass2 < -tolerance || !(E > 0.0)) {
    throw AmplitudeFailure(context + ": invalid positive-energy timelike/null momentum");
  }
  return msqrt(std::max(0.0, mass2));
}

// Build a Dirac-basis particle helicity spinor with a known on-shell mass
MDirac::Spinor DiracParticleHelicitySpinor(const M4Vec &p, int helicity, double mass, const std::string &context) {
  const double denom = p.E() + mass;
  if (!(denom > 0.0) || !std::isfinite(denom)) {
    throw AmplitudeFailure(context + ": invalid positive-energy momentum");
  }
  const double               n      = msqrt(denom);  // u^dagger u = 2E
  const double               ratio  = p.P3mod() / denom;
  const double               theta2 = p.Theta() / 2.0;
  const double               phi    = p.Phi();
  const double               c      = std::cos(theta2);
  const double               s      = std::sin(theta2);
  const std::complex<double> phase  = std::exp(zi * phi);

  if (helicity == 1) { return {n * c, n * phase * s, n * ratio * c, n * ratio * phase * s}; }
  return {-n * s, n * phase * c, n * ratio * s, -n * ratio * phase * c};
}

// Build a Dirac-basis antiparticle helicity spinor with a known on-shell mass
MDirac::Spinor DiracAntiparticleHelicitySpinor(const M4Vec &p, int helicity, double mass, const std::string &context) {
  const MDirac::Spinor particle = DiracParticleHelicitySpinor(p, -helicity, mass, context);
  return {particle[2], particle[3], particle[0], particle[1]};
}

// Transverse spin-1 helicity basis ordered as {-1,+1}.  Construct both states
// together so theta/phi and their trigonometric functions are evaluated once
//
std::array<Tensor1<std::complex<double>, 4>, 2> TransverseSpin1Basis(const M4Vec &k) {
  const double     theta    = k.Theta();
  const double     phi      = k.Phi();
  const double     ct       = std::cos(theta);
  const double     st       = std::sin(theta);
  const double     cp       = std::cos(phi);
  const double     sp       = std::sin(phi);
  constexpr double invsqrt2 = 0.707106781186547524400844362104849039;

  std::array<Tensor1<std::complex<double>, 4>, 2> eps;
  for (const auto &i : indices(eps)) {
    const int    m   = (i == 0) ? -1 : 1;
    const double neg = -static_cast<double>(m);
    const double pos = static_cast<double>(m);

    eps[i](0) = 0.0;
    eps[i](1) = (neg * ct * cp + zi * sp) * invsqrt2;
    eps[i](2) = (neg * ct * sp - zi * cp) * invsqrt2;
    eps[i](3) = pos * st * invsqrt2;
  }
  return eps;
}

// Massive spin-1 basis ordered as {-1,0,+1}.
std::array<Tensor1<std::complex<double>, 4>, 3> MassiveSpin1Basis(const M4Vec &k, const std::string &context) {
  const double M = k.M();
  const double E = k.E();
  if (!(M > 0.0) || !std::isfinite(M) || !std::isfinite(E)) {
    throw AmplitudeFailure(context + ": requires finite timelike momentum with M > 0");
  }

  const double     theta    = k.Theta();
  const double     phi      = k.Phi();
  const double     ct       = std::cos(theta);
  const double     st       = std::sin(theta);
  const double     cp       = std::cos(phi);
  const double     sp       = std::sin(phi);
  constexpr double invsqrt2 = 0.707106781186547524400844362104849039;

  std::array<Tensor1<std::complex<double>, 4>, 3> eps;

  // lambda = -1
  eps[0](0) = 0.0;
  eps[0](1) = (ct * cp + zi * sp) * invsqrt2;
  eps[0](2) = (ct * sp - zi * cp) * invsqrt2;
  eps[0](3) = -st * invsqrt2;

  // lambda = 0
  eps[1](0) = k.P3mod() / M;
  eps[1](1) = E * st * cp / M;
  eps[1](2) = E * st * sp / M;
  eps[1](3) = E * ct / M;

  // lambda = +1
  eps[2](0) = 0.0;
  eps[2](1) = (-ct * cp + zi * sp) * invsqrt2;
  eps[2](2) = (-ct * sp - zi * cp) * invsqrt2;
  eps[2](3) = st * invsqrt2;

  return eps;
}

// Lightlike reference vector opposite to the spatial momentum.  With this
// choice the projected reference-spinor construction is a true helicity state
// (up to an irrelevant global phase), and p.l = E + |p| is never collinearly
// singular for a positive-energy moving particle
M4Vec OppositeLightlikeReference(const M4Vec &p, const std::string &context) {
  const double p3 = p.P3mod();
  if (!(p3 > 0.0) || !std::isfinite(p3)) {
    throw AmplitudeFailure(context + ": helicity reference is undefined for |p| = 0");
  }

  const double invp = 1.0 / p3;
  const double lx   = -p[1] * invp;
  const double ly   = -p[2] * invp;
  const double lz   = -p[3] * invp;
  const double lE   = std::sqrt(lx * lx + ly * ly + lz * lz);

  // (px,py,pz,E)
  return M4Vec(lx, ly, lz, lE);
}

// Apply conjugation and Lorentz-index placement to a polarization basis
template <std::size_t N>
void ApplyPolarizationConvention(std::array<Tensor1<std::complex<double>, 4>, N> &states, const std::string &type,
                                 bool index_up, const std::string &context) {
  const bool conjugate = type == "conj";
  if (!conjugate && type != "none") { throw std::invalid_argument(context + ": type must be none or conj"); }
  for (auto &state : states) {
    if (conjugate) {
      for (std::size_t mu = 0; mu < 4; ++mu) { state(mu) = std::conj(state(mu)); }
    }
    if (!index_up) {
      for (std::size_t mu = 1; mu < 4; ++mu) { state(mu) = -state(mu); }
    }
  }
}

// Transform the real and imaginary parts of one complex Lorentz current
template <typename Transform>
MDirac::Current TransformCurrent(const MDirac::Current &current, Transform transform) {
  M4Vec real(-current[1].real(), -current[2].real(), -current[3].real(), current[0].real());
  M4Vec imag(-current[1].imag(), -current[2].imag(), -current[3].imag(), current[0].imag());
  transform(real);
  transform(imag);
  return {{{real.E(), imag.E()}, {-real.Px(), -imag.Px()}, {-real.Py(), -imag.Py()}, {-real.Pz(), -imag.Pz()}}};
}

// Assemble one four-component operator from two-component spin blocks
MMatrix<std::complex<double>> SpinorBlocks(const MMatrix<std::complex<double>> &a,
                                           const MMatrix<std::complex<double>> &b,
                                           const MMatrix<std::complex<double>> &c,
                                           const MMatrix<std::complex<double>> &d) {
  MMatrix<std::complex<double>> out(4, 4, 0.0);
  out.SetBlock(0, 0, a);
  out.SetBlock(0, 1, b);
  out.SetBlock(1, 0, c);
  out.SetBlock(1, 1, d);
  return out;
}

}  // namespace

// Store basis-independent immutable Dirac algebra
struct MDirac::CommonAlgebra {
  MMatrix<Complex> S_basis = {std::array<Complex, 4>{1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0), 0.0},
                              std::array<Complex, 4>{0.0, 1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)},
                              std::array<Complex, 4>{-1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0), 0.0},
                              std::array<Complex, 4>{0.0, -1.0 / std::sqrt(2.0), 0.0, 1.0 / std::sqrt(2.0)}};
  MMatrix<Complex> sigma_x = {std::array<Complex, 2>{0.0, 1.0}, std::array<Complex, 2>{1.0, 0.0}};
  MMatrix<Complex> sigma_y = {std::array<Complex, 2>{0.0, -zi}, std::array<Complex, 2>{zi, 0.0}};
  MMatrix<Complex> sigma_z = {std::array<Complex, 2>{1.0, 0.0}, std::array<Complex, 2>{0.0, -1.0}};
  std::array<MMatrix<Complex>, 3> angular_momentum;
  std::array<std::size_t, 4>      LI          = {0, 1, 2, 3};
  std::array<int, 2>              SPINORSTATE = spin::BinaryHelicityLabelsX2();
  MMatrix<double>                 g           = MMatrix<double>(4, 4, "minkowski");
  MMatrix<Complex>                I4          = MMatrix<Complex>(4, 4, "eye");

  // Construct the immutable basis-independent matrices once
  CommonAlgebra() : angular_momentum{sigma_x * 0.5, sigma_y * 0.5, sigma_z * 0.5} {}
};

// Store immutable gamma matrices for one representation
struct MDirac::BasisAlgebra {
  std::string                                    BASIS;
  std::array<MMatrix<Complex>, 5>                gamma_up;
  std::array<MMatrix<Complex>, 5>                gamma_lo;
  std::array<std::array<MMatrix<Complex>, 4>, 4> sigma_up;
  std::array<std::array<MMatrix<Complex>, 4>, 4> sigma_lo;
  MMatrix<Complex>                               charge_conjugation;
  MMatrix<Complex>                               right_projector;
  MMatrix<Complex>                               left_projector;

  // Construct all constant matrices for one gamma representation
  explicit BasisAlgebra(const std::string &basis, const CommonAlgebra &common) : BASIS(basis) {
    const MMatrix<Complex> zero(2, 2, 0.0);
    const MMatrix<Complex> identity2(2, 2, "eye");
    const auto             y0_chiral = SpinorBlocks(zero, identity2, identity2, zero);
    const auto             y0_dirac  = SpinorBlocks(identity2, zero, zero, -identity2);
    const auto             y5_chiral = SpinorBlocks(-identity2, zero, zero, identity2);
    const auto            &y5_dirac  = y0_chiral;

    const auto y1_up = SpinorBlocks(zero, common.sigma_x, -common.sigma_x, zero);
    const auto y2_up = SpinorBlocks(zero, common.sigma_y, -common.sigma_y, zero);
    const auto y3_up = SpinorBlocks(zero, common.sigma_z, -common.sigma_z, zero);

    const MMatrix<Complex> &y0 = BASIS == "D" ? y0_dirac : y0_chiral;
    const MMatrix<Complex> &y5 = BASIS == "D" ? y5_dirac : y5_chiral;

    gamma_up = {y0, y1_up, y2_up, y3_up, y5};
    gamma_lo = {y0, -y1_up, -y2_up, -y3_up, y5};

    for (std::size_t mu = 0; mu < 4; ++mu) {
      for (std::size_t nu = 0; nu < 4; ++nu) {
        sigma_up[mu][nu] = (gamma_up[mu] * gamma_up[nu] - gamma_up[nu] * gamma_up[mu]) * (zi / 2.0);
        sigma_lo[mu][nu] = (gamma_lo[mu] * gamma_lo[nu] - gamma_lo[nu] * gamma_lo[mu]) * (zi / 2.0);
      }
    }
    const MMatrix<Complex> identity4(4, 4, "eye");
    charge_conjugation = -gamma_up[2] * gamma_up[0] * zi;
    right_projector    = (identity4 + gamma_up[4]) * 0.5;
    left_projector     = (identity4 - gamma_up[4]) * 0.5;
  }
};

// Compute the immutable basis-independent algebra
const MDirac::CommonAlgebra &MDirac::CommonConstants() {
  static const CommonAlgebra constants;
  return constants;
}

// Compute the immutable algebra for the selected gamma representation
const MDirac::BasisAlgebra &MDirac::BasisConstants(const std::string &basis) {
  if (basis == "DIRAC") {
    static const BasisAlgebra constants("D", CommonConstants());
    return constants;
  }
  if (basis == "CHIRAL") {
    static const BasisAlgebra constants("C", CommonConstants());
    return constants;
  }
  throw std::invalid_argument("MDirac: unknown gamma basis: " + basis);
}

MDirac::MDirac() : MDirac("DIRAC") {}

MDirac::MDirac(const std::string &basis)
    : common_algebra(CommonConstants()),
      basis_algebra(BasisConstants(basis)),
      S_basis(common_algebra.S_basis),
      sigma_x(common_algebra.sigma_x),
      sigma_y(common_algebra.sigma_y),
      sigma_z(common_algebra.sigma_z),
      gamma_up(basis_algebra.gamma_up),
      gamma_lo(basis_algebra.gamma_lo),
      sigma_up(basis_algebra.sigma_up),
      sigma_lo(basis_algebra.sigma_lo),
      LI(common_algebra.LI),
      SPINORSTATE(common_algebra.SPINORSTATE),
      g(common_algebra.g),
      I4(common_algebra.I4),
      BASIS(basis_algebra.BASIS) {}

// Compute chirality projectors valid in every gamma-matrix representation
// P_R,L = (1 +/- gamma5)/2
// Right handed
const MMatrix<std::complex<double>> &MDirac::PR() const { return basis_algebra.right_projector; }
// Left handed
const MMatrix<std::complex<double>> &MDirac::PL() const { return basis_algebra.left_projector; }

// Compute a two-component Weyl spinor in the Jacob-Wick azimuth chart
// sigma.p chi_h = h|p| chi_h
//
MDirac::WeylSpinor MDirac::XiSpinor(const M4Vec &p, int helicity) const {
  const double theta2 = p.Theta() / 2.0;
  const double phi    = p.Phi();
  const double c      = std::cos(theta2);
  const double s      = std::sin(theta2);

  if (helicity == 1) { return {c, s * std::exp(zi * phi)}; }
  return {-s * std::exp(-zi * phi), c};
}

// ----------------------------------------------------------------------
// Helicity eigenstate spinors for particles in the covariant azimuth chart
//
// <@@ DEFINED IN CHIRAL GAMMA MATRIX REPRESENTATION @@>
//
MDirac::Spinor MDirac::uHelChiral(const M4Vec &p, int helicity) const {

  const double E  = p.E();
  const double p3 = p.P3mod();
  const double m  = FermionMass(p, "MDirac::uHelChiral");

  if (!(E + m > 0.0)) { throw AmplitudeFailure("MDirac::uHelChiral: invalid positive-energy timelike/null momentum"); }

  const double               theta2 = p.Theta() / 2.0;
  const double               phi    = p.Phi();
  const double               c      = std::cos(theta2);
  const double               ss     = std::sin(theta2);
  const std::complex<double> phase  = std::exp(zi * phi);

  // E + m - |p| suffers a severe cancellation for E >> m.  Since p.M() is
  // obtained from the same four-vector, use E-|p| = m^2/(E+|p|)
  const double Ep = E + p3;
  if (!(Ep > 0.0) || !std::isfinite(Ep)) { throw AmplitudeFailure("MDirac::uHelChiral: invalid E + |p|"); }
  const double neg = std::fpclassify(m) == FP_ZERO ? 0.0 : m + (m * m) / Ep;
  const double pos = E + m + p3;
  const double N   = 1.0 / (std::sqrt(2.0) * msqrt(E + m));  // u^dagger u = 2E

  if (helicity == 1) { return {c * neg * N, ss * phase * neg * N, c * pos * N, ss * phase * pos * N}; }
  return {-ss * pos * N, c * phase * pos * N, -ss * neg * N, c * phase * neg * N};
}

// Helicity eigenstate spinors for Anti-Particles
//
// <@@ DEFINED IN CHIRAL GAMMA MATRIX REPRESENTATION @@>
//
MDirac::Spinor MDirac::vHelChiral(const M4Vec &p, int helicity) const {

  // Flip the helicity, then use the particle solution with the standard
  // charge-conjugate permutation for this representation/convention
  const Spinor v = uHelChiral(p, -helicity);
  return {-v[0], -v[1], v[2], v[3]};
}
// -----------------------------------------------------------------------

// -----------------------------------------------------------------------
// Helicity eigenstate spinors
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
// [REFERENCE: Thomson, Modern Particle Physics, Cambridge University Press]
// [https://www.hep.phy.cam.ac.uk/~thomson/partIIIparticles/handouts/Handout_2_2011.pdf]
//
// In high energy limit E >> m, helicity eigenstates are eigenstates of gamma^5
// gamma^5 u_+ = +u_+
// gamma^5 u_- = -u_-
// gamma^5 v_+ = -v_+
// gamma^5 v_- = +v_-
//
// Remember: In the limit E >> m (only then)
// -> left and right handed chiral states == helicity states
//
MDirac::Spinor MDirac::uHelDirac(const M4Vec &p, int helicity) const {
  return DiracParticleHelicitySpinor(p, helicity, FermionMass(p, "MDirac::uHelDirac"), "MDirac::uHelDirac");
}

// Helicity eigenstate spinors for Anti-Fermions
//
// <@@ DEFINED IN DIRAC GAMMA-MATRIX REPRESENTATION @@>
//
MDirac::Spinor MDirac::vHelDirac(const M4Vec &p, int helicity) const {
  return DiracAntiparticleHelicitySpinor(p, helicity, FermionMass(p, "MDirac::vHelDirac"), "MDirac::vHelDirac");
}

// Compute the external state phase converting an elastic current to the collider
// section
MDirac::Complex MDirac::ElasticSpinHalfColliderCurrentPhase(const FermionKind fermion_kind, const int beam_leg,
                                                            const double lambda_in, const double lambda_out,
                                                            const double transfer_azimuth) {
  const int sign = fermion_kind == FermionKind::Particle ? 1 : -1;

  const double helicity_in  = 2.0 * lambda_in;
  const double helicity_out = 2.0 * lambda_out;
  if (beam_leg == 1) { return std::exp(zi * 0.5 * (sign - helicity_out) * transfer_azimuth); }
  return -sign * helicity_in * std::exp(zi * 0.5 * (sign + helicity_out) * transfer_azimuth);
}

// -----------------------------------------------------------------------
// Dirac-Spinor "spin eigenstate" fermions: (\gamma^\mu p_\mu - m)u = 0
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
MDirac::Spinor MDirac::uDirac(const M4Vec &p, int spin) const {
  const double E     = p.E();
  const double m     = FermionMass(p, "MDirac::uDirac");
  const double denom = E + m;
  if (!(denom > 0.0)) { throw AmplitudeFailure("MDirac::uDirac: invalid positive-energy timelike/null momentum"); }
  const double N = msqrt(denom);  // u^dagger u = 2E

  if (spin == 1) { return {N, 0.0, N * p[3] / denom, N * (p[1] + zi * p[2]) / denom}; }
  return {0.0, N, N * (p[1] - zi * p[2]) / denom, -N * p[3] / denom};
}

// Dirac-Spinor "spin eigenstate" for anti-fermions: (\gamma^\mu p_\mu + m)v = 0
//
// <@@ DEFINED IN DIRAC GAMMA MATRIX REPRESENTATION @@>
//
MDirac::Spinor MDirac::vDirac(const M4Vec &p, int spin) const {

  // Flip the spin, then use the u-particle solution permutated
  spin           = -spin;
  const Spinor v = uDirac(p, spin);
  return {v[2], v[3], v[0], v[1]};
}
// -----------------------------------------------------------------------

// Compute Hermitian angular momentum operators J_i = 1/2 sigma_i, i = 1,2,3
const MMatrix<std::complex<double>> &MDirac::J_operator(unsigned int i) const {
  if (i >= 1 && i <= 3) { return common_algebra.angular_momentum[i - 1]; }

  throw std::invalid_argument("MDirac::J_operator: i is invalid (not 1,2,3)");
}

// Transformation bilinears:
//
// Vector current: \bar{u}_1 \gamma^\mu u_2
// Axial current:  \bar{u}_1 \gamma^\mu \gamma_5 \u_2
// Tensor current: \bar{u}_1 \sigma^{\mu\nu} \gamma_5 \u_2
// Pseudoscalar current:  \bar{u}_1 \gamma_5 \u_2
//
//
// Adjoing gamma matrix: gamma^\mu\dagger = gamma^0 gamma^\mu gamma^0
//
//
// Charge conjugation operator matrix is basis dependent
// u = C \bar{v}^T <=> \bar{v}^T = C^{-1} u
//
//
// Basis independent definition:
// C^{-1} \gamma_\mu C = -(\gamma_\mu)^T
// equivalently C (\gamma_\mu)^T C^{-1} = -\gamma_\mu
// C^\dagger = C^{-1}
// C^T = -C
//
const MMatrix<std::complex<double>> &MDirac::C_up() const {
  if (BASIS == "D" || BASIS == "C") { return basis_algebra.charge_conjugation; }
  throw std::invalid_argument("MDirac::C_up: unknown gamma basis");
}

// --------------------------------------------------------------------
// Basic notes:
//
// 1. Fermion helicity is conserved only in the massless/chirality-conserving
//    limit. Massive fermions and helicity-flip interactions need both states
// 2. Gauge invariance permits convenient polarization representatives
// Individual
//    gauge-dependent pieces need not be invariant; a physical gauge-invariant
//    amplitude must be
// 3. Sum incoherently only over mutually orthogonal, unobserved external states
//    (and an incoherent initial density matrix). Intermediate helicity
//    amplitudes generally interfere and must be summed at amplitude level
//
// --------------------------------------------------------------------
// Conventions:
//
// Spin-1/2:
// incoming particle:      u(p)
// outgoing particle:      \bar{u}(p)
// incoming anti-particle: \bar{v}(p)
// outgoing anti-particle: v(p)
//
// Propagators:
// spin 1:                 -ig_\mu\nu / q^2
// spin 1/2:                i(\gamma^\mu q_\mu + m)/(q^2 - m^2)
//
// Vertex factor:
//                         ie\gamma^\mu
//
// Matrix element: -iM     product of rules
// --------------------------------------------------------------------

// Photon propagator: -ig_{\mu\nu} / q^2
//
// Remember:
// \sum_{4 virtual polarizations} \eps_\mu^\lambda (\eps_{\nu}^\lambda)^* =
// -g_{\mu\nu}
//
Tensor2<std::complex<double>, 4, 4> MDirac::iD_y(const double q2) const {
  Tensor2<std::complex<double>, 4, 4> T;

  for (const auto &u : LI) {
    for (const auto &v : LI) { T(u, v) = -zi * g[u][v] / q2; }
  }
  return T;
}

// Internal Fermion propagator (matrix with spinor indices, no Lorentz indices)
//
// Input as contravariant (upper) index 4-vector
//
MMatrix<std::complex<double>> MDirac::iD_F(const M4Vec &q, double m) const {
  return (FSlash(q) + I4 * m) * (zi / (q.M2() - pow2(m)));
}

// Polarization vector eps^{(m)\mu}(k) for massless spin-1 with helicity m =
// -1,1
//
// Spatial dependence only on the direction
//
// [REFERENCE: http://scipp.ucsc.edu/~haber/ph218/polsum.pdf]
// [t,x,y,z] order convention!
//
Tensor1<std::complex<double>, 4> MDirac::EpsSpin1(const M4Vec &k, int m) const {
  const auto eps = TransverseSpin1Basis(k);
  return eps[(m + 1) / 2];
}

// Polarization vector eps^{(m)\mu}(k) for massive spin-1 with helicity m =
// -1,0,1
//
// m = +-1 state depends only on the direction of the momentum
// m = 0   state depends also on magnitude
//
// [t,x,y,z] order convention!
//
// Should obey sum:
// \sum_{\lambda = -1,0,1} \eps_\mu(k, \lambda) \eps_\nu(k, \lambda)*
// = -g_{\mu\nu} + k_\mu k_\nu / m^2
//
Tensor1<std::complex<double>, 4> MDirac::EpsMassiveSpin1(const M4Vec &k, int m) const {
  const auto eps = MassiveSpin1Basis(k, "MDirac::EpsMassiveSpin1");
  return eps[m + 1];
}

// Polarization tensor eps^{(m)\mu\nu}(k) for massive spin-2 with helicity m =
// -2,-1,0,1,2
//
// [t,x,y,z] convention!
//
// Should obey sum: (eps_{\mu\nu}^{(m)}(k))^* eps^{(n)\mu\nu}(k) ) = \delta_{mn}
//
Tensor2<std::complex<double>, 4, 4> MDirac::EpsMassiveSpin2(const M4Vec &k, int m) const {

  // Exact 1 x 1 -> 2 Clebsch-Gordan decomposition.  These coefficients are
  // constants, so do not call a generic CG implementation 144 times/state
  const auto  eps1 = MassiveSpin1Basis(k, "MDirac::EpsMassiveSpin2");
  const auto &em   = eps1[0];
  const auto &e0   = eps1[1];
  const auto &ep   = eps1[2];

  constexpr double invsqrt2 = 0.707106781186547524400844362104849039;
  constexpr double invsqrt6 = 0.408248290463863016366214012450981899;

  Tensor2<std::complex<double>, 4, 4> epsmat;
  for (const auto &mu : LI) {
    for (const auto &nu : LI) {
      if (m == 2) {
        epsmat(mu, nu) = ep(mu) * ep(nu);
      } else if (m == 1) {
        epsmat(mu, nu) = (ep(mu) * e0(nu) + e0(mu) * ep(nu)) * invsqrt2;
      } else if (m == 0) {
        epsmat(mu, nu) = (ep(mu) * em(nu) + em(mu) * ep(nu) + 2.0 * e0(mu) * e0(nu)) * invsqrt6;
      } else if (m == -1) {
        epsmat(mu, nu) = (em(mu) * e0(nu) + e0(mu) * em(nu)) * invsqrt2;
      } else {  // m == -2
        epsmat(mu, nu) = em(mu) * em(nu);
      }
    }
  }
  return epsmat;
}

// Adjoint Dirac spinor: \bar{u} = u^dagger * gamma^0
MDirac::Spinor MDirac::Bar(const Spinor &spinor) const {
  // First conjugate elements, then a matrix product with gamma^0 matrix
  return gamma_up[0].LeftMultiply(gra::Conjugated(spinor));
}

// Feynman slash matrix operator: \slash{a} = \gamma_\mu a^\mu = \gamma^\mu
// a_\mu
//
// Input assumed contravariant (upper) index 4-vector
MMatrix<std::complex<double>> MDirac::FSlash(const M4Vec &a) const {
  MMatrix<std::complex<double>> aslash(4, 4, 0.0);  // Init with zero!
  for (const auto &mu : LI) { aslash += gamma_up[mu] * (a % mu); }
  return aslash;
}

// ----------------------------------------------------------------------
// [REFERENCE: https://arxiv.org/pdf/hep-ph/0110108.pdf]

// Spinor product: s_\lambda(p1,p2)
std::complex<double> MDirac::sProd(const M4Vec &p1, const M4Vec &p2, int helicity) const {
  return gra::BilinearProduct(Bar(uGauge(p1, helicity)), uGauge(p2, -helicity));
}

namespace {

// Construct a helicity spinor from a lightlike reference opposite to its momentum
MDirac::Spinor GaugeSpinor(const MDirac &dirac, const M4Vec &p, int helicity, int sign,
                           MDirac::Spinor (MDirac::*spinor)(const M4Vec &, int) const) {
  if (std::fpclassify(p.P3mod()) == FP_ZERO) { return (dirac.*spinor)(p, helicity); }
  const M4Vec  reference = OppositeLightlikeReference(p, "MDirac::GaugeSpinor");
  const double norm      = 2.0 * (p * reference);
  if (!(norm > 0.0) || !std::isfinite(norm)) { throw AmplitudeFailure("MDirac::GaugeSpinor: invalid normalization"); }
  const double mass = FermionMass(p, "MDirac::GaugeSpinor");
  return ((dirac.FSlash(p) * static_cast<double>(sign) + dirac.I4 * mass) / msqrt(norm)) *
         (dirac.*spinor)(reference, -helicity);
}

}  // namespace

// Compute the particle helicity spinor from its massless reference
MDirac::Spinor MDirac::uGauge(const M4Vec &p, int helicity) const {
  return GaugeSpinor(*this, p, helicity, 1, BASIS == "D" ? &MDirac::uHelDirac : &MDirac::uHelChiral);
}

// Compute the antiparticle helicity spinor from its massless reference
MDirac::Spinor MDirac::vGauge(const M4Vec &p, int helicity) const {
  return GaugeSpinor(*this, p, helicity, -1, BASIS == "D" ? &MDirac::vHelDirac : &MDirac::vHelChiral);
}

// Construct spin-1/2 helicity spinors (-1,1) [indexing with 0,1]
std::array<MDirac::Spinor, 2> MDirac::SpinorStates(const M4Vec &p, const std::string &type) const {
  if (type != "u" && type != "ubar" && type != "v" && type != "vbar") {
    throw std::invalid_argument("MDirac::SpinorStates: type must be u, ubar, v or vbar");
  }

  std::array<Spinor, 2> spinor{};
  for (const auto &m : indices(spinor)) {
    const int h = SPINORSTATE[m];
    if (BASIS == "D") {
      if (type == "u") {
        spinor[m] = uHelDirac(p, h);
      } else if (type == "ubar") {
        spinor[m] = Bar(uHelDirac(p, h));
      } else if (type == "v") {
        spinor[m] = vHelDirac(p, h);
      } else {
        spinor[m] = Bar(vHelDirac(p, h));
      }
    } else {
      if (type == "u") {
        spinor[m] = uHelChiral(p, h);
      } else if (type == "ubar") {
        spinor[m] = Bar(uHelChiral(p, h));
      } else if (type == "v") {
        spinor[m] = vHelChiral(p, h);
      } else {
        spinor[m] = Bar(vHelChiral(p, h));
      }
    }
  }
  return spinor;
}

// Construct Massless Spin-1 polarization vectors (m = -1,1) [indexing with 0,1]
//
// Input p as contravariant (upper index)
//
std::array<Tensor1<std::complex<double>, 4>, 2> MDirac::MasslessSpin1States(const M4Vec &p, const std::string &type,
                                                                            bool INDEX_UP) const {
  auto eps = TransverseSpin1Basis(p);
  ApplyPolarizationConvention(eps, type, INDEX_UP, "MDirac::MasslessSpin1States");
  return eps;
}

// Construct Massive Spin-1 polarization vectors (m = -1,0,1) [indexing with
// 0,1,2]
//
// Input p as contravariant (upper index)
//
std::array<Tensor1<std::complex<double>, 4>, 3> MDirac::MassiveSpin1States(const M4Vec &p, const std::string &type,
                                                                           bool INDEX_UP) const {
  auto eps = MassiveSpin1Basis(p, "MDirac::MassiveSpin1States");
  ApplyPolarizationConvention(eps, type, INDEX_UP, "MDirac::MassiveSpin1States");
  return eps;
}

// Compute a lower-index Dirac-Pauli current
// The momentum transfer convention is q = outgoing - incoming and the result has a lower Lorentz index
// Antiparticles follow vbar(in) Gamma(q) v(out) with charge kept external
MDirac::Current MDirac::DiracPauliCurrent(const gra::M4Vec &incoming, const gra::M4Vec &outgoing,
                                          const gra::M4Vec &vertex_q, const FermionKind fermion_kind, double lambda_in,
                                          double lambda_out, double mass, double F1, double F2) const {


  const int                  h_in            = lambda_in > 0.0 ? 1 : -1;
  const int                  h_out           = lambda_out > 0.0 ? 1 : -1;
  const std::complex<double> pauli_prefactor = gra::math::zi * F2 / (2.0 * mass);

  Spinor left;
  Spinor right;
  switch (fermion_kind) {
    case FermionKind::Particle:
      left  = DiracParticleHelicitySpinor(outgoing, h_out, mass, "MDirac::DiracPauliCurrent outgoing");
      right = DiracParticleHelicitySpinor(incoming, h_in, mass, "MDirac::DiracPauliCurrent incoming");
      break;
    case FermionKind::Antiparticle:
      left  = DiracAntiparticleHelicitySpinor(incoming, h_in, mass, "MDirac::DiracPauliCurrent incoming");
      right = DiracAntiparticleHelicitySpinor(outgoing, h_out, mass, "MDirac::DiracPauliCurrent outgoing");
      break;
    default:
      throw std::invalid_argument("MDirac::DiracPauliCurrent: unknown fermion kind");
  }

  // Express both spinors in the selected gamma basis before taking the adjoint
  if (BASIS == "C") {
    const auto to_chiral = S_basis.Dagger();
    left = to_chiral * left;
    right = to_chiral * right;
  }
  left = Bar(left);

  // Evaluate the Dirac and Pauli bilinears without matrix temporaries
  Current current{};
  for (const auto &mu : LI) {
    const std::complex<double> dirac = gamma_lo[mu].BilinearForm(left, right);

    std::complex<double> pauli = 0.0;
    for (const auto &nu : LI) {
      if (nu == mu || std::fpclassify(vertex_q[nu]) == FP_ZERO) { continue; }
      pauli += vertex_q[nu] * sigma_lo[mu][nu].BilinearForm(left, right);
    }
    current[mu] = F1 * dirac + pauli_prefactor * pauli;
  }
  return current;
}

// Compute the complete outgoing-incoming Dirac-Pauli helicity basis
std::array<MDirac::Current, 4> MDirac::DiracPauliCurrentBasis(const M4Vec &incoming, const M4Vec &outgoing,
                                                              const M4Vec &vertex_q, const FermionKind fermion_kind,
                                                              const double mass, const double F1,
                                                              const double F2) const {

  constexpr std::array<int, 2> helicity = {-1, 1};
  std::array<Spinor, 2>        left;
  std::array<Spinor, 2>        right;
  for (std::size_t h = 0; h < 2; ++h) {
    if (fermion_kind == FermionKind::Particle) {
      left[h] =
          DiracParticleHelicitySpinor(outgoing, helicity[h], mass, "MDirac::DiracPauliCurrentBasis outgoing");
      right[h] = DiracParticleHelicitySpinor(incoming, helicity[h], mass, "MDirac::DiracPauliCurrentBasis incoming");
    } else if (fermion_kind == FermionKind::Antiparticle) {
      left[h] =
          DiracAntiparticleHelicitySpinor(incoming, helicity[h], mass, "MDirac::DiracPauliCurrentBasis incoming");
      right[h] =
          DiracAntiparticleHelicitySpinor(outgoing, helicity[h], mass, "MDirac::DiracPauliCurrentBasis outgoing");
    } else {
      throw std::invalid_argument("MDirac::DiracPauliCurrentBasis: unknown fermion type");
    }
    // Use the same gamma representation for every external helicity state
    if (BASIS == "C") {
      const auto to_chiral = S_basis.Dagger();
      left[h] = to_chiral * left[h];
      right[h] = to_chiral * right[h];
    }
    left[h] = Bar(left[h]);
  }

  const std::complex<double> pauli_prefactor = math::zi * F2 / (2.0 * mass);
  std::array<Current, 4>     basis{};
  for (std::size_t out = 0; out < 2; ++out) {
    for (std::size_t in = 0; in < 2; ++in) {
      const Spinor &bra     = fermion_kind == FermionKind::Particle ? left[out] : left[in];
      const Spinor &ket     = fermion_kind == FermionKind::Particle ? right[in] : right[out];
      Current      &current = basis[2 * out + in];
      for (const auto &mu : LI) {
        std::complex<double> pauli = 0.0;
        for (const auto &nu : LI) {
          if (nu != mu) { pauli += vertex_q[nu] * sigma_lo[mu][nu].BilinearForm(bra, ket); }
        }
        current[mu] = F1 * gamma_lo[mu].BilinearForm(bra, ket) + pauli_prefactor * pauli;
      }
    }
  }
  return basis;
}

// Lorentz boost one lower-index complex current as a physical four-vector
MDirac::Current MDirac::BoostCurrent(const Current &current, const M4Vec &boost, const double mass, const int sign) {
  return TransformCurrent(current, [&](M4Vec &part) { kinematics::LorentzBoost(boost, mass, part, sign); });
}

// Rotate one lower-index complex current with a spatial rotation matrix
MDirac::Current MDirac::RotateCurrent(const Current &current, const MMatrix<double> &rotation) {
  if (rotation.size_row() != 3 || rotation.size_col() != 3) {
    throw std::invalid_argument("MDirac::RotateCurrent: rotation matrix must be 3 by 3");
  }
  return TransformCurrent(current, [&](M4Vec &part) { part.Rotate(rotation); });
}

}  // namespace gra
