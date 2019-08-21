// Spin polarization functions
//
// (Wigner-D functions, Clebsch-Gordans, Arbitrary Spin-Density matrix)
//
//
// In general, functions here use helicities for fermions which are
// not normalized by the spin vector norm.
// That is, functions use [-1/2, 1/2] not for example [-1, 1].
//
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++ standard
#include <algorithm>
#include <complex>
#include <iostream>
#include <random>
#include <vector>

// Eigen
#include <Eigen/Dense>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MSpin.h"

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;

namespace gra {
namespace spin {
// Used for checking eigenvalues
const double epsilon = 1e-6;

// Initialize helicity decay amplitude matrix T (2s1 + 1) x (2s2 + 1),
//
// where s1, s2 are daughter spins (0, 1/2, 1, ...)
//
// T_{\lambda_1,\lambda_2}
//   = \sum_{ls} \alpha_{ls} x
//          < J\lambda|ls0\lambda> <s\lambda|s1s2\lambda1,-\lambda2>
//
void InitTMatrix(gra::HELMatrix &hc, const gra::MParticle &p, const gra::MParticle &p1,
                 const gra::MParticle &p2) {
  const double J = p.spinX2 / 2.0;
  const int    P = p.P;

  const double s1 = p1.spinX2 / 2.0;
  const int    P1 = p1.P;

  const double s2 = p2.spinX2 / 2.0;
  const int    P2 = p2.P;

  // Boson or Fermion
  const bool is_boson1 = (cint(s1)) ? true : false;
  const bool is_boson2 = (cint(s2)) ? true : false;
  const bool identical = (p1.pdg == p2.pdg) ? true : false;  // e.g. pi0 pi0 vs pi+pi-

  // -------------------------------------------------------------------
  // Helicity decay amplitude matrix

  hc.T = MMatrix<std::complex<double>>(static_cast<unsigned int>(2 * s1 + 1),
                                       static_cast<unsigned int>(2 * s2 + 1), 0.0);

  MMatrix<bool> need_to_set(20, 20, false);
  bool          tag_set = false;
  unsigned int  nonzero = 0;

  // Construct T-matrix (2s1 + 1) x (2s2 + 1) elements
  std::cout << std::endl;
  std::cout << "gra::spin:: "
            << " Particle 1: " << (is_boson1 ? "boson" : "fermion")
            << " , Particle 2: " << (is_boson2 ? "boson" : "fermion") << std::endl;

  std::cout << "gra::spin::InitTMatrix: Calculating SU(2) decomposition [lambda = lambda1 - lambda2]: " << std::endl;
  
  for (int s = 0; s <= static_cast<int>(s1 + s2); ++s) {
    for (int l = 0; l <= static_cast<int>(J + s); ++l) {
      for (int i = 0; i < static_cast<int>(2 * s1 + 1); ++i) {
        for (int j = 0; j < static_cast<int>(2 * s2 + 1); ++j) {
          const double lambda1 = -s1 + static_cast<double>(i);  // From negative to positive
          const double lambda2 = -s2 + static_cast<double>(j);  // From negative to positive

          // Jacob-Wick
          const double lambda = lambda1 - lambda2;

          // Check angular momentum conservation for z-axis {
          if (!(abs(lambda1 - lambda2) <= J)) {
            continue;  // Not conserved
          }
          // <J \lambda | ls0\lambda>
          const double cg1 = gra::spin::ClebschGordan(
              static_cast<double>(l), static_cast<double>(s), 0.0, lambda, J, lambda);

          // <s \lambda | s1s2\lambda1,-\lambda2>
          const double cg2 =
              gra::spin::ClebschGordan(s1, s2, lambda1, -lambda2, static_cast<double>(s), lambda);

          // Angular momentum conserved (should be, by construction
          // here)
          if (!(std::abs(lambda) <= J)) { continue; }  // not conserved

          // If parity conserving decay
          if (hc.P_conservation) {
            // Parity of the final state (hold for both
            // bosons/fermions)
            const int P_tot = P1 * P2 * std::pow(-1, l);
            if (P_tot != P) { continue; }  // not conserved
          }

          // Bose-Symmetry for identical Boson-Pairs (symmetric
          // wavefunction)
          if (is_boson1 && is_boson2 && identical) {
            if (!BoseSymmetry(l, s)) { continue; }  // not right
          }

          // Fermi-Symmetry for identical Fermion-Pairs (antisymmetric
          // wavefunction)
          if (!is_boson1 && !is_boson2 && identical) {
            if (!FermiSymmetry(l, s)) { continue; }  // not right
          }

          // CG-coupling product is non-zero
          if (std::abs(cg1 * cg2) > 1e-15) {
            // Sum
            hc.T[i][j] += hc.alpha[l][s] * cg1 * cg2;

            // Has not been set, throw error next
            if (hc.alpha_set[l][s] == false) {
              need_to_set[l][s] = true;
              tag_set           = true;
            }
            printf(
                "[l=%d,s=%d] :: lambda1 = %4.1f, lambda2 = "
                "%4.1f : cg1*cg2 = %6.3f \n",
                l, s, lambda1, lambda2, cg1 * cg2);
            ++nonzero;
          }
        }
      }
    }
  }

  const std::string str0 = "gra::spin::InitTMatrix:: [" + gra::aux::Spin2XtoString(2 * J) + "^" +
                           gra::aux::ParityToString(P) + " -> " + gra::aux::Spin2XtoString(2 * s1) +
                           "^" + gra::aux::ParityToString(P1) + " " +
                           gra::aux::Spin2XtoString(2 * s2) + "^" + gra::aux::ParityToString(P2) +
                           "]" + " (identical = " + std::string(identical ? "true" : "false") + ")";

  if (nonzero == 0) {
    const std::string str = ": Decay is impossible (spin-parity-statistics) [.P_conservation = " +
                            (hc.P_conservation ? std::string("true]") : std::string("false]"));
    throw std::invalid_argument(str0 + str);
  }
  gra::matoper::PrintMatrixSeparate(hc.T);

  // Do we have missing values?
  if (tag_set) {
    std::string str = "";
    for (std::size_t l = 0; l < need_to_set.size_row(); ++l) {
      for (std::size_t s = 0; s < need_to_set.size_col(); ++s) {
        if (need_to_set[l][s] == true) {
          str += "(" + std::to_string(l) + "," + std::to_string(s) + ") ";
        }
      }
    }

    const std::string middle =
        ", Following helicity decay alpha_ls-couplings "
        "missing from BRANCHING: (l,s) "
        "= ";
    throw std::invalid_argument(str0 + middle + str);
  }
}

// |l-s| <= J <= l + s, parity P = (-1)^l
// l is orbital angular momentum
// s is total spin
//
// Bose-Einstein statistics requires l-s to be even
// for identical Boson pairs.
bool BoseSymmetry(int l, int s) {
  const int number = l - s;
  if (number % 2 == 0) { return true; }
  return false;
}

// Fermi-Dirac statistics requires l+s to be even
// for identical Fermion pairs.
bool FermiSymmetry(int l, int s) {
  const int number = l + s;
  if (number % 2 == 0) { return true; }
  return false;
}

// Recursive decay chain with spin correlations:
// A -> B + C, B -> B1 + B2, C -> C1 + C2
//
// TBD. Write down a fully recursive (infinite length) 2->2 chain version
//
std::complex<double> SpinAmp(const gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) {
  const unsigned int B = 0;
  const unsigned int C = 1;

  MMatrix<std::complex<double>> fA;
  MMatrix<std::complex<double>> fB;
  MMatrix<std::complex<double>> fC;

  // This function handles only 2-body decays
  if (lts.decaytree.size() != 2) { return resonance.g_decay; }

  // No spin, no need
  if (resonance.p.spinX2 == 0) {
    return resonance.hc.T[0][0] * resonance.g_decay;  // (l,s) = (0,0)
  }

  // Transition amplitude matrix
  // A -> B + C (in A rest frame)
  if (lts.decaytree.size() == 2) {
    // Get daughter spins and parities
    const double s1 = lts.decaytree[B].p.spinX2 / 2.0;
    const double s2 = lts.decaytree[C].p.spinX2 / 2.0;

    // Boost daughter to the A-rest frame
    M4Vec boosted_daughter = lts.decaytree[B].p4;
    gra::kinematics::LorentzBoost(lts.pfinal[0], lts.pfinal[0].M(), boosted_daughter, -1);

    fA = fMatrix(resonance.hc.T, resonance.p.spinX2 / 2.0, s1, s2, boosted_daughter.Theta(),
                 boosted_daughter.Phi());
  } else {
    return 1.0;
  }

  // Transition amplitude matrix
  // B -> B1 + B2
  if (lts.decaytree[B].legs.size() == 2) {
    const unsigned int B1 = 0;
    const unsigned int B2 = 1;

    // Get daughter spins and parities
    const double s1 = lts.decaytree[B].legs[B1].p.spinX2 / 2.0;
    const double s2 = lts.decaytree[B].legs[B2].p.spinX2 / 2.0;

    // Constant 1.0 for now (solid for "Magic Decays", i.e., which have T == 1)
    const MMatrix<std::complex<double>> T((unsigned int)(2 * s1 + 1), (unsigned int)(2 * s2 + 1),
                                          1.0);

    // Rotate and boost daughters to the frame where z-spanned by B-flight
    // direction
    // (B-helicity rest frame)
    std::vector<M4Vec> daughter = {lts.decaytree[B].legs[B1].p4, lts.decaytree[B].legs[B2].p4};
    gra::kinematics::HEframe(daughter);

    fB = fMatrix(T, lts.decaytree[B].p.spinX2 / 2.0, s1, s2, daughter[B1].Theta(),
                 daughter[B1].Phi());
  }

  // Transition amplitude matrix
  // C -> C1 + C2
  if (lts.decaytree[C].legs.size() == 2) {
    const unsigned int C1 = 0;
    const unsigned int C2 = 1;

    // Get daughter spins and parities
    const double s1 = lts.decaytree[C].legs[C1].p.spinX2 / 2.0;
    const double s2 = lts.decaytree[C].legs[C2].p.spinX2 / 2.0;

    // Constant 1.0 for now (solid for "Magic Decays", i.e., which have T == 1)
    const MMatrix<std::complex<double>> T((unsigned int)(2 * s1 + 1), (unsigned int)(2 * s2 + 1),
                                          1.0);

    // Rotate and boost daughters to the frame where z-spanned by C-flight
    // direction
    // (C-helicity rest frame)
    std::vector<M4Vec> daughter = {lts.decaytree[C].legs[C1].p4, lts.decaytree[C].legs[C2].p4};
    gra::kinematics::HEframe(daughter);

    fC = fMatrix(T, lts.decaytree[C].p.spinX2 / 2.0, s1, s2, daughter[C1].Theta(),
                 daughter[C1].Phi());
  }

  // Total transition amplitude matrix as a tensor product
  MMatrix<std::complex<double>> f;
  if (fB.size_row() > 0 && fC.size_row() > 0) {
    const MMatrix<std::complex<double>> prod = gra::matoper::TensorProd(fB, fC);
    f                                        = prod * fA;  // Matrix x Matrix product
  } else {
    f = fA;  // No recursion
  }

  // ------------------------------------------------------------------
  // Construct the D-matrix for an event-by-event spin space density operator
  // rotation
  double theta_rotation = 0.0;
  double phi_rotation   = 0.0;

  // Spin polarization density matrix defined:

  // In Helicity frame
  if (resonance.hc.FRAME == "HE") {
    theta_rotation = lts.pfinal[0].Theta();
    phi_rotation   = lts.pfinal[0].Phi();

    // In Non-Rotated frame
  } else if (resonance.hc.FRAME == "SR") {
    theta_rotation = 0;
    phi_rotation   = 0;

    // In Collins-Soper frame
  } else if (resonance.hc.FRAME == "CS") {
    theta_rotation = lts.pfinal[0].Theta();
    phi_rotation   = -lts.pfinal[0].Phi();
  } else {
    // Throw exception
    std::string str =
        "gra::spin::SpinAmp2: Unknown polarization Lorentz rest FRAME chosen: "
        "\"" +
        resonance.hc.FRAME +
        "\" (valid are: \"SR\" (no rotation), \"HE\" (helicity), \"CS\" "
        "(Collins-Soper) )";
    throw std::invalid_argument(str);
  }

  // ------------------------------------------------------------------
  // Rotation does mixing of spin states. N.B. Eigenvalues do not change in
  // rotation.
  const MMatrix<std::complex<double>> D =
      gra::spin::DMatrix(resonance.p.spinX2 / 2.0, theta_rotation, phi_rotation);

  // rho_rot = D*rho*D^dagger
  const MMatrix<std::complex<double>> rho_ROT = D * resonance.hc.rho * D.Dagger();

  // Weight (amplitude squared) of the event by the density matrix formalism:
  // Tr[f*rho*f^dagger]
  const MMatrix<std::complex<double>> fRfd = f * rho_ROT * f.Dagger();

  const double A = msqrt(std::real(fRfd.Trace()));

  // (2J+1) is normalization to match production coupling with spin-0
  // (same cross section obtained as with spin-0 if density matrix == I ==
  // totally
  // unpolarized)
  const double NORM = msqrt(resonance.p.spinX2 + 1.0);

  return NORM * A * resonance.g_decay;
}

// Create spin projections -J <= M <= J (number of 2J+1)
// negative to positive
//
// Spin-1/2:  [-1/2, 1/2]
// Spin-1:    [-1, 0, 1]
// Spin-3/2:  [-3/2, -1/2, 1/2, 3/2]
// Spin-2:    [-2, -1, 0, 1, 2]
//
std::vector<double> SpinProjections(double J) {
  std::vector<double> M(static_cast<int>(2 * J + 1), 0.0);
  double              m = -J;
  for (std::size_t i = 0; i < M.size(); ++i) {
    M[i] += m;
    m += 1.0;
  }
  return M;
}

// Returns transition amplitude matrix of size [(2s1 + 1)x(2s2 + 1)] x (2J + 1)
//
// f_{\lambda_1\lambda_2,M}(\theta,\phi) = D_{\lambda M}^J (\theta,\phi)
// T_{\lambda_1, \lambda_2}
//
// theta,phi defined in the decay rest frame
//
// [REFERENCE: Jacob, Wick, On the general theory of particles with spin, 1959]
// [REFERENCE: Amsler, Bizot, Simulation of angular distributions, 1983]
MMatrix<std::complex<double>> fMatrix(const MMatrix<std::complex<double>> &T, double J, double s1,
                                      double s2, double theta, double phi) {
  if (T.size_row() != static_cast<unsigned int>(2 * s1 + 1) ||
      T.size_col() != static_cast<unsigned int>(2 * s2 + 1)) {
    throw std::invalid_argument(
        "gra::spin::fMatrix: Input T matrix "
        "(2s1+1)x(2s2+1) with wrong dimensions: " +
        std::to_string(T.size_row()) + " x " + std::to_string(T.size_col()));
  }
  // Number of final state configurations
  const unsigned int N =
      static_cast<unsigned int>(2 * s1 + 1) * static_cast<unsigned int>(2 * s2 + 1);

  // Initial state projections
  const std::vector<double> M = SpinProjections(J);

  // -------------------------------------------------------------------
  // Construct final state helicity configurations

  // For indexing T matrix
  MMatrix<double>       lambda_values(N, 2, 0.0);
  MMatrix<unsigned int> lambda_index(N, 2, 0);

  unsigned int i = 0;
  for (int mu = 0; mu < static_cast<int>(2 * s1 + 1); ++mu) {
    for (int nu = 0; nu < static_cast<int>(2 * s2 + 1); ++nu) {
      lambda_values[i][0] = -s1 + static_cast<double>(mu);
      lambda_values[i][1] = -s2 + static_cast<double>(nu);

      lambda_index[i][0] = mu;
      lambda_index[i][1] = nu;
      ++i;
    }
  }

  // -------------------------------------------------------------------

  // Construct transition amplitude matrix
  MMatrix<std::complex<double>> f(N, 2 * J + 1, 0.0);

  // Rows = final state spin projections
  for (std::size_t i = 0; i < N; ++i) {
    // Final state helicities
    const double lambda1 = lambda_values[i][0];
    const double lambda2 = lambda_values[i][1];

    // Total lambda by definition
    const double lambda = lambda1 - lambda2;

    // Columns = initial state polarizations
    for (std::size_t j = 0; j < 2 * J + 1; ++j) {
      // Wigner D * Helicity amplitude
      f[i][j] = WignerD(theta, phi, lambda, M[j], J) * T[lambda_index[i][0]][lambda_index[i][1]];
    }
  }
  // std::cout << "f:::" << std::endl;
  // gra::matoper::PrintMatrixSeparate(f);
  // exit(1);

  return f;
}

// Construct Wigner D-matrix for a spin-space operator rotation
MMatrix<std::complex<double>> DMatrix(double J, double theta_rotation, double phi_rotation) {
  // Create spin projections
  const std::vector<double> M = SpinProjections(J);

  const unsigned int n = static_cast<unsigned int>(2 * J + 1);

  MMatrix<std::complex<double>> D(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      D[i][j] = WignerD(theta_rotation, phi_rotation, M[i], M[j], J);
    }
  }
  return D;
}

// SU(2) Clebsch-Gordan coefficients
double ClebschGordan(double j1, double j2, double m1, double m2, double j, double m) {
  // Selection rules
  if (!CGrules(j1, j2, m1, m2, j, m)) { return 0.0; }

  // Use Wigner-3j symbol to evaluate, note minus on m!
  return std::pow(-1.0, m + j1 - j2) * msqrt(2 * j + 1) * W3j(j1, j2, j, m1, m2, -m);
}

// Check if half-integer
bool chalfint(double x) {
  if (std::abs(2 * x - std::floor(2 * x)) < 1E-9) { return true; }
  return false;
}

// Check if integer
bool cint(double x) {
  if (std::abs(x - std::floor(x)) < 1E-9) { return true; }
  return false;
}

// Floating point equality comparison
bool isequal(double x, double y) {
  if (std::abs(x - y) < 1E-9) { return true; }
  return false;
}

// CG-selection rules
double CGrules(double j1, double j2, double m1, double m2, double j, double m) {
  const bool PRINT_ON = false;

  if (!chalfint(j1) || !chalfint(j2) || !chalfint(j)) {
    if (PRINT_ON) printf("CGrules:: (j1,j2,j) not integer or half-integer\n");

    return false;
  }
  if (!cint(j1 + j2 + j)) {  // "Integer perimeter rule"
    if (PRINT_ON) printf("CGrules:: j1 + j2 + j not integer\n");

    return false;
  }
  if (!chalfint(m1) || !chalfint(m2) || !chalfint(m)) {
    if (PRINT_ON) printf("CGrules:: (m1,m2,m) not integer or half-integer\n");

    return false;
  }
  if (!isequal(m1 + m2, m)) {
    if (PRINT_ON) printf("CGrules:: m1 + m2 != m\n");

    return false;
  }
  if (!isequal(j1 - m1, std::floor(j1 - m1))) {
    if (PRINT_ON) printf("CGrules:: parity of 2j_1 and 2m_1 does not match\n");

    return false;
  }
  if (!isequal(j2 - m2, std::floor(j2 - m2))) {
    if (PRINT_ON) printf("CGrules:: parity of 2j_2 and 2m_2 does not match\n");

    return false;
  }
  if (!isequal(j - m, std::floor(j - m))) {
    if (PRINT_ON) printf("CGrules:: parity of 2j and 2m does not match\n");

    return false;
  }
  if ((j1 + j2) < j || j < std::abs(j1 - j2)) {
    if (PRINT_ON) printf("CGrules:: j over valid domain\n");

    return false;
  }
  if (std::abs(m1) > j1) {
    if (PRINT_ON) printf("CGrules:: |m1| > j1\n");

    return false;
  }
  if (std::abs(m2) > j2) {
    if (PRINT_ON) printf("CGrules:: |m2| > j2\n");

    return false;
  }
  if (std::abs(m) > j) {
    if (PRINT_ON) printf("CGrules:: |m| > j\n");

    return false;
  }
  return true;
}

// Wigner 3j symbol
//
// (j1 j2 j3)
// (m1 m1 m3)
//
// http://mathworld.wolfram.com/Wigner3j-Symbol.html

double W3j(double j1, double j2, double j3, double m1, double m2, double m3) {
  // Check selection rules
  if (!W3jrules(j1, j2, j3, m1, m2, m3)) { return 0.0; }

  // Create coefficient structure
  const int c1 = -j2 + m1 + j3;
  const int c2 = -j1 - m2 + j3;
  const int c3 = j1 + j2 - j3;
  const int c4 = j1 - m1;
  const int c5 = j2 + m2;

  // Boundaries for which factorials are non-negative
  const int lower = std::max(std::max(0, -c1), -c2);
  const int upper = std::min(std::min(c3, c4), c5);
  double    fsum  = 0;

  // Evaluate sum term
  for (int k = lower; k <= upper; ++k) {
    fsum += std::pow(-1.0, k) / (gra::math::factorial(k) * gra::math::factorial(c1 + k) *
                                 gra::math::factorial(c2 + k) * gra::math::factorial(c3 - k) *
                                 gra::math::factorial(c4 - k) * gra::math::factorial(c5 - k));
  }

  // Evaluate by Racah formula
  const double wc = std::pow(-1.0, j1 - j2 - m3) *
                    msqrt(TriangleCoeff(j1, j2, j3) * gra::math::factorial(j1 + m1) *
                          gra::math::factorial(j1 - m1) * gra::math::factorial(j2 + m2) *
                          gra::math::factorial(j2 - m2) * gra::math::factorial(j3 + m3) *
                          gra::math::factorial(j3 - m3)) *
                    fsum;
  return wc;
}

// Triangle coefficient
double TriangleCoeff(double j1, double j2, double j3) {
  return gra::math::factorial(j1 + j2 - j3) * gra::math::factorial(j1 - j2 + j3) *
         gra::math::factorial(-j1 + j2 + j3) / gra::math::factorial(j1 + j2 + j3 + 1);
}

// Selection rules
bool W3jrules(double j1, double j2, double j3, double m1, double m2, double m3) {
  const bool PRINT_ON = false;

  if (!cint(j1 + j2 + j3)) {
    if (PRINT_ON) { std::cout << "W3jSR:: j1 + j2 + j3 not integer\n"; }

    return false;
  }
  if (!isequal(m1 + m2 + m3, 0)) {
    if (PRINT_ON) { std::cout << "W3jSR:: m1 + m2 + m3 != 0\n"; }

    return false;
  }
  if (!isequal(j1 - m1, std::floor(j1 - m1))) {
    if (PRINT_ON) { std::cout << "W3jSR:: parity of 2j_1 and 2m_1 does not match\n"; }

    return false;
  }
  if (!isequal(j2 - m2, std::floor(j2 - m2))) {
    if (PRINT_ON) { std::cout << "W3jSR:: parity of 2j_2 and 2m_2 does not match\n"; }

    return false;
  }
  if (!isequal(j3 - m3, std::floor(j3 - m3))) {
    if (PRINT_ON) { std::cout << "W3jSR:: parity of 2j_3 and 2m_3 does not match\n"; }

    return false;
  }
  if ((j1 + j2) < j3 || j3 < std::abs(j1 - j2)) {
    if (PRINT_ON) { std::cout << "W3jSR:: j over valid domain\n"; }

    return false;
  }
  if (std::abs(m1) > j1) {
    if (PRINT_ON) { std::cout << "W3jSR:: |m1| > j1\n"; }

    return false;
  }
  if (std::abs(m2) > j2) {
    if (PRINT_ON) { std::cout << "W3jSR:: |m2| > j2\n"; }

    return false;
  }
  if (std::abs(m3) > j3) {
    if (PRINT_ON) { std::cout << "W3jSR:: |m3| > j3\n"; }

    return false;
  }
  return true;
}

// Calculate quantum mechanical von Neumann entropy of a density matrix
double VonNeumannEntropy(const MMatrix<std::complex<double>> &rho) {
  const unsigned int N = rho.size_col();
  Eigen::MatrixXcd   A(N, N);

  // Collect matrix elements
  for (std::size_t i = 0; i < N; ++i) {
    for (std::size_t j = 0; j < N; ++j) { A(i, j) = rho[i][j]; }
  }

  // Solve eigenvalues
  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXcd> es(A);
  if (es.info() != Eigen::Success) { return -1.0; }

  // Calculate von Neumann Entropy
  double entropy = 0.0;
  for (std::size_t i = 0; i < N; ++i) {
    std::complex<double> lambda = es.eigenvalues().col(0)[i];
    const double         eigval = std::real(lambda);  // Take real part for numerics
    if (eigval > epsilon) { entropy += eigval * std::log(eigval); }
  }
  entropy = -entropy;
  return entropy;
}

// Wigner small-d (real valued convention, z-y-z convention)
//
// Wigner small d_{mm"}^J(theta) ~ rotation around y_z-axis
//
// [REFERENCE: en.wikipedia.org/wiki/Wigner_D-matrix#Wigner_(small)_d-matrix]
//
// m, mp, J all in integers (not half-integers)
//
double WignerSmalld(double theta, double m, double mp, double J) {
  const double factor =
      gra::math::msqrt(gra::math::factorial(J + m) * gra::math::factorial(J - m) *
                       gra::math::factorial(J + mp) * gra::math::factorial(J - mp));
  double s = 0.0;

  // k = 0,1,...,2J
  for (int k = 0; k <= 2 * J; ++k) {
    double value = 0;
    double val1  = 0.0;
    double val2  = 0.0;

    // Ignore terms with negative factorials
    if (J - mp - k >= 0 && J + m - k >= 0 && k + mp - m >= 0 && k >= 0) {
      val1 = std::pow(-1, mp - m + k) /
             (gra::math::factorial(J - mp - k) * gra::math::factorial(J + m - k) *
              gra::math::factorial(k + mp - m) * gra::math::factorial(k));

      val2 = std::pow(std::cos(theta / 2.0), 2 * J + m - mp - 2 * k) *
             std::pow(std::sin(theta / 2.0), mp - m + 2 * k);

      value = val1 * val2;
    }
    s += value;
  }
  // Fatal error
  if (std::isnan(s)) {
    std::string str = "gra::spin::WignerSmalld: NaN value with (theta,m,mp,J) = " +
                      std::to_string(theta) + " " + std::to_string(m) + " " + std::to_string(mp) +
                      " " + std::to_string(J);
    throw std::invalid_argument(str);
  }
  return factor * s;
}

// Wigner large-D
// Wigner large D(theta,phi) = exp(iJ_y theta) exp(iJ_z phi)
// D = @(theta,phi,m,mp,J) wignersmalld(theta,m,mp,J)*exp(1i*mp*phi);
//
// m, mp, J all in integers (not half-integers)
//
std::complex<double> WignerD(double theta, double phi, double m, double mp, double J) {
  return WignerSmalld(theta, m, mp, J) * std::exp(std::complex<double>(0, 1) * mp * phi);
}

// Density matrix properties:
//
//  1. rho^2 = rho        (projector)
//  2. rho^\dagger = rho  (hermiticity)
//  3. Tr[rho] = 1        (normalization)
//  4. rho >= 0           (positivity, eigenvalues greater or equal to zero)
//
// Test positivity of the density matrix,
// - Also checks the trivial fact that rho
//   has right dimensions given J, that is (2J+1)
// - Also checks that the normalization is ok.
//
//
// TBD, write down generalization for any spin (e.g. done with Eigen library
// eigenvalues)

bool Positivity(const MMatrix<std::complex<double>> &rho, unsigned int J) {
  const unsigned int n = rho.size_row();

  // Dimensions
  if ((n - 1) / 2 != J) {
    std::string str =
        "gra::spin::Positivity: Density matrix does not have a valid "
        "dimension (2J+1) !";
    throw std::invalid_argument(str);
    return false;
  }

  // Normalization, abs taken for comparison
  const double tracerho = std::abs(rho.Trace());
  if (std::abs(tracerho - 1.0) > 1e-2) {
    std::string str =
        "gra::spin::Positivity: Density matrix is not properly "
        "normalized (Tr[rho] = " +
        std::to_string(tracerho) + ") (should be 1) !";
    gra::matoper::PrintMatrixSeparate(rho);
    throw std::invalid_argument(str);
    return false;
  }

  const double NEG_EPS = -1e6;  // Negative for floating point reasons

  // Spin-1
  if (J == 1) {
    // Component as:
    //
    //    rho_i = [0.5*(1-a), b+1i*c, d;
    //             b-1i*c,   a, -b+1i*c;
    //             d, -b-1i*c, 0.5*(1-a)];

    // Test positivity of rho (ANALYTICAL)
    const double a = std::real(rho[1][1]);
    const double b = std::real(rho[0][1]);
    const double c = std::imag(rho[0][1]);
    const double d = std::real(rho[0][2]);

    const double pc1 = -4 * b * b - 4 * c * c - 2 * d * d + a - (3 * a * a) / 2.0 + 1 / 2.0;
    const double pc2 =
        -(3 * (2 * d - a + 1) * (2 * a * d - a + a * a + 4 * b * b + 4 * c * c)) / 2.0;

    if (pc1 >= NEG_EPS && pc2 >= NEG_EPS) {
      return true;
    } else {
      printf("pc1 = %0.5f, pc2 = %0.5f \n", pc1, pc2);
      std::string str =
          "gra::spin::Positivity: Eigenvalues of the density "
          "matrix are not >= 0 !";
      gra::matoper::PrintMatrixSeparate(rho);
      throw std::invalid_argument(str);
      return false;
    }
  }

  // Spin-2
  else if (J == 2) {
    // Components as:
    //
    //    rho_i  = [a,       f+1i*c,   g+1i*d,      h+1i*e,      m;
    //             f-1i*c,   b,        j+1i*k,      l, -h+1i*e;
    //             g-1i*d,   j-1i*k,   1-2*(a+b),  -j+1i*k, g-1i*d;
    //             h-1i*e,   l,       -j-1i*k,      b, -f+1i*c;
    //             m,       -h-1i*e,   g+1i*d,     -f-1i*c, a];

    // Test positivity of rho (ANALYTICAL)
    const double a = std::real(rho[0][0]);
    const double b = std::real(rho[1][1]);
    const double c = std::imag(rho[0][1]);
    const double d = std::imag(rho[0][2]);
    const double e = std::imag(rho[0][3]);
    const double f = std::real(rho[0][1]);
    const double g = std::real(rho[0][2]);
    const double h = std::real(rho[0][3]);
    const double j = std::real(rho[1][2]);
    const double k = std::imag(rho[1][2]);
    const double l = std::real(rho[1][3]);
    const double m = std::real(rho[0][4]);

    const double pc1 = -6 * a * a - 8 * a * b + 4 * a - 6 * b * b + 4 * b - 4 * c * c - 4 * d * d -
                       4 * e * e - 4 * f * f - 4 * g * g - 4 * h * h - 4 * j * j - 4 * k * k -
                       2 * l * l - 2 * m * m;
    const double pc2 = -12 * a * a * a - 48 * a * a * b + 6 * a * a - 48 * a * b * b + 24 * a * b +
                       12 * a * c * c - 12 * a * d * d + 12 * a * e * e + 12 * a * f * f -
                       12 * a * g * g + 12 * a * h * h - 24 * a * j * j - 24 * a * k * k +
                       12 * a * m * m - 12 * b * b * b + 6 * b * b + 12 * b * c * c -
                       24 * b * d * d + 12 * b * e * e + 12 * b * f * f - 24 * b * g * g +
                       12 * b * h * h - 12 * b * j * j - 12 * b * k * k + 12 * b * l * l -
                       12 * c * c + 24 * c * d * j + 24 * c * e * l - 24 * c * e * m -
                       24 * c * g * k + 12 * d * d * m - 24 * d * e * j + 24 * d * f * k -
                       24 * d * h * k - 12 * e * e + 24 * e * g * k - 12 * f * f + 24 * f * g * j +
                       24 * f * h * l - 24 * f * h * m + 12 * g * g * m - 24 * g * h * j -
                       12 * h * h - 12 * j * j * l - 12 * k * k * l - 6 * l * l - 6 * m * m;
    const double pc3 =
        48 * a * a * b - 168 * a * a * b * b - 96 * a * a * a * b + 96 * a * a * c * c +
        96 * a * a * e * e + 96 * a * a * f * f + 96 * a * a * h * h - 48 * a * a * j * j -
        48 * a * a * k * k + 72 * a * a * l * l - 96 * a * b * b * b + 48 * a * b * b +
        144 * a * b * c * c - 96 * a * b * d * d + 144 * a * b * e * e + 144 * a * b * f * f -
        96 * a * b * g * g + 144 * a * b * h * h - 96 * a * b * j * j - 96 * a * b * k * k +
        96 * a * b * l * l + 96 * a * b * m * m - 48 * a * c * c + 96 * a * c * d * j -
        96 * a * c * e * l + 192 * a * c * e * m - 96 * a * c * g * k - 96 * a * d * e * j +
        96 * a * d * f * k - 96 * a * d * h * k - 48 * a * e * e + 96 * a * e * g * k -
        48 * a * f * f + 96 * a * f * g * j - 96 * a * f * h * l + 192 * a * f * h * m -
        96 * a * g * h * j - 48 * a * h * h - 96 * a * j * j * l - 96 * a * k * k * l -
        48 * a * l * l + 96 * b * b * c * c - 48 * b * b * d * d + 96 * b * b * e * e +
        96 * b * b * f * f - 48 * b * b * g * g + 96 * b * b * h * h + 72 * b * b * m * m -
        48 * b * c * c + 96 * b * c * d * j - 192 * b * c * e * l + 96 * b * c * e * m -
        96 * b * c * g * k + 96 * b * d * d * m - 96 * b * d * e * j + 96 * b * d * f * k -
        96 * b * d * h * k - 48 * b * e * e + 96 * b * e * g * k - 48 * b * f * f +
        96 * b * f * g * j - 192 * b * f * h * l + 96 * b * f * h * m + 96 * b * g * g * m -
        96 * b * g * h * j - 48 * b * h * h - 48 * b * m * m + 24 * c * c * c * c +
        48 * c * c * d * d - 48 * c * c * e * e + 48 * c * c * f * f + 48 * c * c * g * g +
        48 * c * c * h * h + 48 * c * c * j * j + 48 * c * c * k * k + 48 * c * c * l * m +
        96 * c * d * d * e + 96 * c * d * j * l - 96 * c * d * j * m - 192 * c * e * f * h +
        96 * c * e * g * g + 96 * c * e * j * j + 96 * c * e * k * k + 96 * c * e * l -
        96 * c * e * m - 96 * c * g * k * l + 96 * c * g * k * m + 48 * d * d * e * e +
        48 * d * d * f * f + 96 * d * d * f * h + 48 * d * d * h * h + 48 * d * d * l * l -
        96 * d * e * j * l + 96 * d * e * j * m + 96 * d * f * k * l - 96 * d * f * k * m -
        96 * d * h * k * l + 96 * d * h * k * m + 24 * e * e * e * e + 48 * e * e * f * f +
        48 * e * e * g * g + 48 * e * e * h * h + 48 * e * e * j * j + 48 * e * e * k * k +
        48 * e * e * l * m + 96 * e * g * k * l - 96 * e * g * k * m + 24 * f * f * f * f +
        48 * f * f * g * g - 48 * f * f * h * h + 48 * f * f * j * j + 48 * f * f * k * k +
        48 * f * f * l * m + 96 * f * g * g * h + 96 * f * g * j * l - 96 * f * g * j * m +
        96 * f * h * j * j + 96 * f * h * k * k + 96 * f * h * l - 96 * f * h * m +
        48 * g * g * h * h + 48 * g * g * l * l - 96 * g * h * j * l + 96 * g * h * j * m +
        24 * h * h * h * h + 48 * h * h * j * j + 48 * h * h * k * k + 48 * h * h * l * m +
        48 * j * j * m * m + 48 * k * k * m * m + 24 * l * l * m * m;
    const double pc4 =
        120 *
        (c * c + 2 * c * e + e * e + f * f + 2 * f * h + h * h - a * b - a * l + b * m + l * m) *
        (a * l - 2 * c * e - a * b - 2 * f * h - b * m + l * m + 2 * a * b * b + 2 * a * a * b -
         2 * a * c * c - 2 * b * c * c - 2 * a * e * e + 2 * b * d * d - 2 * a * f * f -
         2 * b * e * e - 2 * b * f * f - 2 * a * h * h + 2 * b * g * g - 2 * b * h * h +
         2 * a * j * j + 2 * a * k * k - 2 * a * a * l + 2 * b * b * m - 2 * d * d * l -
         2 * g * g * l + 2 * j * j * m + 2 * k * k * m + c * c + e * e + f * f + h * h +
         4 * a * c * e + 4 * b * c * e - 2 * a * b * l + 4 * a * f * h + 2 * a * b * m +
         4 * b * f * h - 4 * c * d * j + 4 * d * e * j + 4 * c * g * k - 4 * d * f * k +
         4 * d * h * k - 4 * e * g * k - 4 * f * g * j + 4 * g * h * j - 2 * a * l * m -
         2 * b * l * m);

    if (pc1 >= NEG_EPS && pc2 >= NEG_EPS && pc3 >= NEG_EPS && pc4 >= NEG_EPS) {
      return true;
    } else {
      printf(
          "pc1 = %0.5f, pc2 = %0.5f, pc3 = %0.5f, pc4 = "
          "%0.5f \n",
          pc1, pc2, pc3, pc4);
      std::string str =
          "gra::spin::Positivity: Eigenvalues of the density "
          "matrix are not >= 0 !";
      gra::matoper::PrintMatrixSeparate(rho);
      throw std::invalid_argument(str);
      return false;
    }
  } else {
    std::string str =
        "gra::spin::Positivity: Not a valid spin (only J = 1 and 2 "
        "accepted) !";
    gra::matoper::PrintMatrixSeparate(rho);
    throw std::invalid_argument(str);
    return false;
  }
}

// Generate random density matrix for the resonance
// Create random (under Haar / Hilbert-Schmidt measure) spin-density matrices
// for any spin,
// and for spin 1 and 2 the possibility to enforce parity conservation
//
//
// TBD, write down the generalization of parity conservation constaint for any
// spin.
MMatrix<std::complex<double>> RandomRho(unsigned int J, bool parity, MRandom &rng) {
  const std::complex<double>    zi(0, 1);  // imag unit
  MMatrix<std::complex<double>> rho;

  // Draw random density matrices until valid found, takes a couple
  // usually
  while (true) {
    // Create random complex matrix (n x n)
    const unsigned int            n = 2 * J + 1;
    MMatrix<std::complex<double>> r(n, n);

    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        r[i][j] = rng.G(0, 1) + zi * rng.G(0, 1);  // real and imag part from
                                                   // normal distribution
      }
    }
    // Take product (outer product per vector): r*r^\dagger
    rho = r * r.Dagger();

    // Normalize the trace to 1 (real just for conversion to double)
    const double scale = 1.0 / std::real(rho.Trace());
    rho                = rho * scale;

    // If parity conservation required, then symmetrize
    if (parity) {
      if (J == 1) {
        // 1.1 Purely real elements by parity and
        // hermiticity
        rho[0][2] = std::real(rho[0][2]);

        // 1.2 Parity constrained
        rho[1][2] = -std::real(rho[0][1]) + zi * std::imag(rho[0][1]);

        // 2. Lower triangle = hermiticity constrained
        for (std::size_t i = 0; i < n; ++i) {
          for (std::size_t j = 0; j <= i; ++j) { rho[i][j] = std::conj(rho[j][i]); }
        }

        // 3. Diagonal parity constrained
        rho[2][2] = std::real(rho[0][0]);
      }
      if (J == 2) {
        // 1.1 Purely real elements by parity and
        // hermiticity
        rho[0][4] = std::real(rho[0][4]);
        rho[1][3] = std::real(rho[1][3]);

        // 1.2 Parity constrained
        rho[1][4] = -std::real(rho[0][3]) + zi * std::imag(rho[0][3]);
        rho[2][3] = -std::real(rho[1][2]) + zi * std::imag(rho[1][2]);
        rho[2][4] = std::conj(rho[0][2]);
        rho[3][4] = -std::real(rho[0][1]) + zi * std::imag(rho[0][1]);

        // 2. Lower triangle = hermiticity constrained
        for (std::size_t i = 0; i < n; ++i) {
          for (std::size_t j = 0; j <= i; ++j) { rho[i][j] = std::conj(rho[j][i]); }
        }

        // 3. Diagonal parity constrained
        rho[3][3] = std::real(rho[1][1]);
        rho[4][4] = std::real(rho[0][0]);
      }

      // Re-Normalize trace to 1 (real just for conversion to double)
      const double newscale = 1.0 / std::real(rho.Trace());
      rho                   = rho * newscale;
      ;
    }

    // printmatrix(rho);
    // Check that the matrix is valid density matrix in terms of
    // positivity requirement
    try {
      if (Positivity(rho, J)) { break; }
    } catch (...) {
      continue;  // not valid density matrix
    }
  }

  gra::matoper::PrintMatrixSeparate(rho);

  return rho;
}

}  // spin namespace
}  // gra namespace
