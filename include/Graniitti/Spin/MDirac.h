// Dirac gamma algebra routines
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDIRAC_H
#define MDIRAC_H

// C++
#include <array>
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"

// FTensor
#include "FTensor.hpp"

namespace gra {
class MDirac {
  struct CommonAlgebra;
  struct BasisAlgebra;

  // Compute immutable algebra shared by every MDirac instance
  static const CommonAlgebra &CommonConstants();

  // Compute immutable gamma algebra for one representation
  static const BasisAlgebra &BasisConstants(const std::string &basis);

  const CommonAlgebra &common_algebra;
  const BasisAlgebra &basis_algebra;

public:
  using Complex = std::complex<double>;
  using WeylSpinor = std::array<Complex, 2>;
  using Spinor = std::array<Complex, 4>;
  using Current = std::array<Complex, 4>;

  // Select particle or antiparticle external spinors for an elastic current
  enum class FermionKind { Particle, Antiparticle };

  MDirac();
  MDirac(const std::string &basis);
  ~MDirac() = default;

  // Covariant polarization representatives without a JW little group phase
  FTensor::Tensor1<std::complex<double>, 4> EpsSpin1(const M4Vec &k,
                                                     int m) const;
  FTensor::Tensor1<std::complex<double>, 4> EpsMassiveSpin1(const M4Vec &k,
                                                            int m) const;
  FTensor::Tensor2<std::complex<double>, 4, 4> EpsMassiveSpin2(const M4Vec &k,
                                                               int m) const;

  // Spin/spinor state collectors
  std::array<Spinor, 2> SpinorStates(const M4Vec &p,
                                     const std::string &type) const;
  std::array<FTensor::Tensor1<std::complex<double>, 4>, 2>
  MasslessSpin1States(const M4Vec &p, const std::string &type,
                      bool INDEX_UP = true) const;
  std::array<FTensor::Tensor1<std::complex<double>, 4>, 3>
  MassiveSpin1States(const M4Vec &p, const std::string &type,
                     bool INDEX_UP = true) const;

  // Compute the two-component Jacob-Wick helicity representative
  WeylSpinor XiSpinor(const M4Vec &p, int helicity) const;

  // Compute four-component helicity spinors in the covariant azimuth chart
  Spinor uHelChiral(const M4Vec &p, int helicity) const;
  Spinor vHelChiral(const M4Vec &p, int helicity) const;

  Spinor uHelDirac(const M4Vec &p, int helicity) const;
  Spinor vHelDirac(const M4Vec &p, int helicity) const;

  // Compute the external state phase converting an elastic current to the
  // collider section
  static Complex ElasticSpinHalfColliderCurrentPhase(FermionKind fermion_kind,
                                                     int beam_leg,
                                                     double lambda_in,
                                                     double lambda_out,
                                                     double transfer_azimuth);

  // Compute a lower-index Dirac-Pauli current with supplied form factors
  Current DiracPauliCurrent(const M4Vec &incoming, const M4Vec &outgoing,
                            const M4Vec &vertex_q, FermionKind fermion_kind,
                            double lambda_in, double lambda_out, double mass,
                            double F1, double F2) const;

  // Compute all outgoing-incoming elastic spin-half transitions
  std::array<Current, 4> DiracPauliCurrentBasis(
      const M4Vec &incoming, const M4Vec &outgoing, const M4Vec &vertex_q,
      FermionKind fermion_kind, double mass, double F1, double F2) const;

  // Lorentz boost one lower-index complex current as a physical four-vector
  static Current BoostCurrent(const Current &current, const M4Vec &boost,
                              double mass, int sign);

  // Rotate one lower-index complex current with a spatial rotation matrix
  static Current RotateCurrent(const Current &current,
                               const MMatrix<double> &rotation);

  // Dirac adjoint
  Spinor Bar(const Spinor &spinor) const;

  // Feynman slash matrix operator
  MMatrix<std::complex<double>> FSlash(const M4Vec &a) const;

  // Propagators
  FTensor::Tensor2<std::complex<double>, 4, 4> iD_y(const double q2) const;
  MMatrix<std::complex<double>> iD_F(const M4Vec &q, double m) const;

  // Dirac spinors
  Spinor uDirac(const M4Vec &p, int spin) const;
  Spinor vDirac(const M4Vec &p, int spin) const;

  // Spinor-Helicity style methods
  std::complex<double> sProd(const M4Vec &p1, const M4Vec &p2,
                             int helicity) const;
  Spinor uGauge(const M4Vec &p, int helicity) const;
  Spinor vGauge(const M4Vec &p, int helicity) const;

  // Charge conjugate operator
  const MMatrix<std::complex<double>> &C_up() const;

  // Right and left chiral projectors
  const MMatrix<std::complex<double>> &PR() const;
  const MMatrix<std::complex<double>> &PL() const;

  // Angular Momentum operators
  const MMatrix<std::complex<double>> &J_operator(unsigned int i) const;

  // Transform matrix from the Weyl to Dirac basis by an operator sandwich:
  //
  // \gamma^\mu_Dirac = S \gamma_{Weyl}^\mu S^\dagger
  //
  // S^{-1} = S^\dagger
  //
  // Spinor Transformation:
  //
  // u_Dirac = S u_Weyl
  // v_Dirac = S v_Weyl
  //
  // S = 1/\sqrt{2}(1 - y^5y^0)
  //
  const MMatrix<std::complex<double>> &S_basis;

  // Pauli matrices in the standard Z-basis
  const MMatrix<std::complex<double>> &sigma_x;
  const MMatrix<std::complex<double>> &sigma_y;
  const MMatrix<std::complex<double>> &sigma_z;

  // Contravariant (up) and covariant (lo) gamma matrix set
  const std::array<MMatrix<std::complex<double>>, 5> &gamma_up;
  const std::array<MMatrix<std::complex<double>>, 5> &gamma_lo;

  // Contravariant (up) and covariant (lo) \sigma_{\mu\nu} matrices
  const std::array<std::array<MMatrix<std::complex<double>>, 4>, 4> &sigma_up;
  const std::array<std::array<MMatrix<std::complex<double>>, 4>, 4> &sigma_lo;

  // Indices etc
  const std::array<std::size_t, 4> &LI;
  const std::array<int, 2> &SPINORSTATE;

  // Minkowski metric tensor (+,-,-,-) in convention (t,px,py,pz)
  const MMatrix<double> &g;

  // Identity matrix
  const MMatrix<std::complex<double>> &I4;

protected:
  const std::string &BASIS; // D for Dirac, C for Chiral
};

} // namespace gra

#endif
