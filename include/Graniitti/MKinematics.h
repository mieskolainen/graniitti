// Standard relativistic kinematics and related structures [HEADER ONLY file]
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MKINEMATICS_H
#define MKINEMATICS_H

// C++
#include <complex>
#include <random>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MMatOper.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MSudakov.h"

namespace gra {
namespace kinematics {


// Case m3 = m4
inline double SolvePz3_A(double m3, double pt3, double pt4, double pz5, double E5, double s) {
  const double t2  = E5 * E5;
  const double t3  = pt3 * pt3;
  const double t4  = pt4 * pt4;
  const double t5  = pz5 * pz5;
  const double t6  = m3 * m3;
  const double t7  = std::sqrt(s);
  const double t8  = s * t2 * 6.0;
  const double t9  = std::pow(s, 3.0 / 2.0);
  const double t10 = t2 * t2;
  const double t11 = t3 * t3;
  const double t12 = t4 * t4;
  const double t13 = t5 * t5;
  const double t14 = s * s;
  const double t15 = t5 * t6 * 4.0;
  const double t16 = t3 * t5 * 2.0;
  const double t17 = t4 * t5 * 2.0;
  const double t18 = E5 * t6 * t7 * 8.0;
  const double t19 = E5 * t3 * t7 * 4.0;
  const double t20 = E5 * t4 * t7 * 4.0;
  const double t21 = E5 * t5 * t7 * 4.0;
  const double t22 = t8 + t10 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 + t19 + t20 + t21 -
                     E5 * t9 * 4.0 - s * t3 * 2.0 - s * t4 * 2.0 - s * t5 * 2.0 - s * t6 * 4.0 -
                     t2 * t3 * 2.0 - t2 * t4 * 2.0 - t2 * t5 * 2.0 - t3 * t4 * 2.0 - t2 * t6 * 4.0 -
                     E5 * t2 * t7 * 4.0;
  const double t23 = std::sqrt(t22);
  const double t0  = pz5 * (-1.0 / 2.0) -
                    (E5 * t23 * (1.0 / 2.0) - t7 * t23 * (1.0 / 2.0) +
                     pz5 * (t3 * (1.0 / 2.0) - t4 * (1.0 / 2.0))) /
                        (s + t2 - t5 - E5 * t7 * 2.0);

  return t0;
}

/*
// Old version for reference
inline double SolvePz3_B(double m1, double m2, double pt1, double pt2, double
pz, double pE, double
s) {

  using gra::math::msqrt;
  using gra::math::pow2;
  using gra::math::pow3;
  using gra::math::pow4;
  using gra::math::pow5;

  const double sqrt_s = msqrt(s);
  const double p1z =
          (-(pow2(m1) * pow2(pE) * pz) + pow2(m2) * pow2(pE) * pz -
           pow4(pE) * pz - pow2(pE) * pow2(pt1) * pz +
           pow2(pE) * pow2(pt2) * pz + pow2(m1) * pow3(pz) -
           pow2(m2) * pow3(pz) + 2.0 * pow2(pE) * pow3(pz) +
           pow2(pt1) * pow3(pz) - pow2(pt2) * pow3(pz) - pow5(pz) -
           2.0 * pow2(m1) * pE * pz * sqrt_s +
           2.0 * pow2(m2) * pE * pz * sqrt_s -
           2.0 * pE * pow2(pt1) * pz * sqrt_s +
           2.0 * pE * pow2(pt2) * pz * sqrt_s - pow2(m1) * pz * s +
           pow2(m2) * pz * s + 2.0 * pow2(pE) * pz * s - pow2(pt1) * pz * s +
           pow2(pt2) * pz * s + 2.0 * pow3(pz) * s - pz * pow2(s) +
           sqrt(pow2(pE - sqrt_s) *
                        pow2(pow2(pE) - pow2(pz) + 2.0 * pE * sqrt_s + s) *
                        (pow4(m1) + pow4(m2) + pow4(pE) -
                         2.0 * pow2(pE) * pow2(pt1) + pow4(pt1) -
                         2.0 * pow2(pE) * pow2(pt2) - 2.0 * pow2(pt1) * pow2(pt2) +
                         pow4(pt2) - 2.0 * pow2(pE) * pow2(pz) +
                         2.0 * pow2(pt1) * pow2(pz) + 2.0 * pow2(pt2) * pow2(pz) +
                         pow4(pz) - 4.0 * pow3(pE) * sqrt_s +
                         4.0 * pE * pow2(pt1) * sqrt_s +
                         4.0 * pE * pow2(pt2) * sqrt_s +
                         4.0 * pE * pow2(pz) * sqrt_s + 6.0 * pow2(pE) * s -
                         2.0 * pow2(pt1) * s - 2.0 * pow2(pt2) * s -
                         2.0 * pow2(pz) * s - 4.0 * pE * std::pow(s, 1.5) + pow2(s) -
                         2.0 * pow2(m2) * (pow2(pE) + pow2(pt1) - pow2(pt2) -
                                                           pow2(pz) - 2.0 * pE * sqrt_s + s) -
                         2.0 * pow2(m1) *
                                 (pow2(m2) + pow2(pE) - pow2(pt1) + pow2(pt2) -
                                  pow2(pz) - 2.0 * pE * sqrt_s + s)))) /
          (2.0 * (pow4(pE) + pow2(pow2(pz) - s) -
                          2.0 * pow2(pE) * (pow2(pz) + s)));

          return p1z;
}
*/

// Case m3, m4 same or different
inline double SolvePz3_B(double m3, double m4, double pt3, double pt4, double pz5, double E5,
                         double s) {
  const double t2  = E5 * E5;
  const double t3  = m3 * m3;
  const double t4  = m4 * m4;
  const double t5  = pt3 * pt3;
  const double t6  = pt4 * pt4;
  const double t7  = pz5 * pz5;
  const double t8  = std::sqrt(s);
  const double t9  = s * t2 * 6.0;
  const double t10 = std::pow(s, 3.0 / 2.0);
  const double t11 = t2 * t2;
  const double t12 = t3 * t3;
  const double t13 = t4 * t4;
  const double t14 = t5 * t5;
  const double t15 = t6 * t6;
  const double t16 = t7 * t7;
  const double t17 = s * s;
  const double t18 = t3 * t5 * 2.0;
  const double t19 = t4 * t6 * 2.0;
  const double t20 = t3 * t7 * 2.0;
  const double t21 = t4 * t7 * 2.0;
  const double t22 = t5 * t7 * 2.0;
  const double t23 = t6 * t7 * 2.0;
  const double t24 = E5 * t3 * t8 * 4.0;
  const double t25 = E5 * t4 * t8 * 4.0;
  const double t26 = E5 * t5 * t8 * 4.0;
  const double t27 = E5 * t6 * t8 * 4.0;
  const double t28 = E5 * t7 * t8 * 4.0;
  const double t29 = t9 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 + t19 + t20 + t21 + t22 +
                     t23 + t24 + t25 + t26 + t27 + t28 - E5 * t10 * 4.0 - s * t3 * 2.0 -
                     s * t4 * 2.0 - s * t5 * 2.0 - s * t6 * 2.0 - s * t7 * 2.0 - t2 * t3 * 2.0 -
                     t2 * t4 * 2.0 - t2 * t5 * 2.0 - t3 * t4 * 2.0 - t2 * t6 * 2.0 - t2 * t7 * 2.0 -
                     t3 * t6 * 2.0 - t4 * t5 * 2.0 - t5 * t6 * 2.0 - E5 * t2 * t8 * 4.0;
  const double t30 = std::sqrt(t29);
  const double t0 =
      pz5 * (-1.0 / 2.0) -
      (E5 * t30 * (1.0 / 2.0) - t8 * t30 * (1.0 / 2.0) +
       pz5 * (t3 * (1.0 / 2.0) - t4 * (1.0 / 2.0) + t5 * (1.0 / 2.0) - t6 * (1.0 / 2.0))) /
          (s + t2 - t7 - E5 * t8 * 2.0);

  return t0;
}

/*
%% MATLAB symbolic code to generate solutions
%% Solve the non-linear system
syms E3 E4 E5 pz3 pz4 pz5 m3 m4 pt3 pt4 s


assume(E3  > 0);
assume(E4  > 0);
assume(pt3 > 0);
assume(pt4 > 0);
assume(m3  > 0);
assume(m4  > 0);
assume(s   > 0);
assume(pzin,'real');

% Use CMS-frame
eqs = [0       == pz3 + pz4 + pz5
           s^(1/2) == E3 + E4 + E5
           E3^2    == m3^2 + pz3^2 + pt3^2
           E4^2    == m4^2 + pz4^2 + pt4^2].';

S = solve(eqs, [pz3 pz4 E3 E4]);

%% The answer for p3z, rest by substitution
S = simplify(S.pz3, 25);

polybranch = 2; % Polynomial branch
S = S(polybranch);


str = 'const double';

sol_A = simplify(subs(S, [m4], [m3]), 'Steps', 25);
ccode(sol_A, 'File', 'temp.c'); readwritetext('temp.c','solution_A.c', str);

sol_B = simplify(S, 'Steps', 25);
ccode(sol_B, 'File', 'temp.c'); readwritetext('temp.c','solution_B.c', str);
*/

/*
function readwritetext(old_file, new_file, str)
fid_old = fopen(old_file,'r');
fid_new = fopen(new_file,'w');

tline = fgetl(fid_old);
while ischar(tline)
        fprintf(fid_new, '%s%s\n', str, tline);
        tline = fgetl(fid_old);
end
fclose(fid_old);
fclose(fid_new);
end
*/

// This function solves pz3 component (in center-of-momentum frame)
// in p1 + p2 -> p3 + {p5} + p4
//
//
inline double SolvePz(double m3, double m4, double pt3, double pt4, double pz5, double E5,
                      double s) {
  constexpr double EPS = 1e-9;

  if (std::abs(m3 - m4) < EPS) {  // identical mass case (elastic-X-elastic)
    return SolvePz3_A(m3, pt3, pt4, pz5, E5, s);
  } else {
    return SolvePz3_B(m3, m4, pt3, pt4, pz5, E5, s);
  }
}

// Center of mass momentum^2
// for e.g. 1+2 (p_{cm}^2) -> 3+4 (p_{cm}'^2)
//
// dsigma/dOmega (s,theta) = 1/(64*pi^2*s) p_{cm}'/p_{cm} |M|^2,
// where dOmega = dcostheta dphi, if no initial state polarization, no phi dep.
//
// Then integrating phi: dcostheta dphi / (16*pi^2) = dcostheta/(8*pi), gives
// dsigma/dt (s,t) = 1/(16*pi*s^2 * p_{cm}^2) |M(s,t)|^2
//
inline double pcm2(double s, double m1, double m2) {
  return 1 / (4.0 * s) * (s - gra::math::pow2(m1 + m2)) * (s - gra::math::pow2(m1 - m2));
}

// 2->2 scattering invariants (p1 + p2  -> p3 + p4)
//
// -s/2(1 - cos(theta*))
inline double mandelstam_t(const M4Vec &p1, const M4Vec &p3) { return -2.0 * (p1 * p3); }
// -s/2(1 + cos(theta*))
inline double mandelstam_u(const M4Vec &p2, const M4Vec &p4) { return -2.0 * (p2 * p4); }

// 1-body Lorentz invariant phase space volume integral for energy < MAX (GeV)
// \int dPS_1 = \int d^3p / (2E(2pi)^3)
//            = \int_0^{4\pi} d^2\Omega/(2pi)^3 \int_0^{\sqrt{MAX^2-m^2}} p^2
//            dp/(2E)
//            = 4\pi / (2(2\pi)^3) \sqrt{MAX^2 - m^2}^3/(3 MAX)
//
// with E^2 = p^2 + m^2
//
inline double dPhi1(double MAX, double m) {
  using gra::math::pow2;

  return pow2(MAX) / (12 * pow2(gra::math::PI)) * std::pow(1.0 - pow2(m) / pow2(MAX), 3.0 / 2.0);
}

// Two-body Lorentz invariant phase space volume, evaluated in the rest frame
// p_norm is the result given by DecayMomentum()
inline double dPhi2(double E, double p_norm) {
  // Sampling box volume cos(theta) ~ [-1, 1], phi ~ [0,2pi]
  const double V = 2.0 * 2.0 * gra::math::PI;
  return V * p_norm / (16.0 * gra::math::pow2(gra::math::PI) * E);
}

// Kallen triangle function (square root applied here)
// symmetric under interchange of x <-> y <-> z:
// lambda^(1/2) = ( (x - y - z)^2 - 4yz )^(1/2)
//
// 1 / lambda^(1/2)(s,m_a^2,m_b^2) -> 1/(2s) when s >> m_a^2, m_b^2
//
// Note x,y,z are Lorentz scalars (ENERGY squared variables)
inline double SqrtKallenLambda(double x, double y, double z) {
  // Hyperboloid formula
  return gra::math::msqrt(x * x + y * y + z * z - 2.0 * (x * y + x * z + y * z));
  //    return msqrt( pow2( x - y - z) - 4.0*y*z ); // numerically unstable
}

inline double beta12(double s, double m1, double m2) {
  return SqrtKallenLambda(s, gra::math::pow2(m1), gra::math::pow2(m2)) / s;
}

// Kinematic (momentum) bound for phase space (decays)
// This is same as: SqrtKallenLambda(x^2, y^2, z^2) / (2x)
//
// Note x,y,z [mother, daughter1, daughter2] in units of ENERGY (not ENERGY^2)
inline double DecayMomentum(double x, double y, double z) {
  // return SqrtKallenLambda(x*x, y*y, z*z) / (2*x);
  // numerically stable
  return gra::math::msqrt((x - y - z) * (x + y + z) * (x + y - z) * (x - y + z)) / (2.0 * x);
}

// Massive two body phase space integral volume PS^{2}(M^2; m_A^2, m_B^2)
inline double PS2Massive(double M2, double mA2, double mB2) {
  return SqrtKallenLambda(M2, mA2, mB2) / (8.0 * gra::math::PI * M2);
}

// Massless phase space integral volume PS^{(n)} = \int dPS^{(n)}
// where M2 (GeV^2) is the total CMS energy^2 of the n-body system.
// Use this e.g. as a (debug) reference of other routines
// Gives 1/(8*PI) for n = 2 case, regardless of M^2
//
// tgamma is C++/11 math Gamma function
inline double PSnMassless(double M2, int n) {
  // const double denom = tgamma(n) * tgamma(n - 1.0);
  const double denom = gra::math::factorial(n - 1) * gra::math::factorial(n - 2);
  const double E     = gra::math::msqrt(M2);
  const double PS =
      (1.0 / (2.0 * std::pow(4.0 * gra::math::PI, 2 * n - 3))) * (std::pow(E, 2 * n - 4) / denom);

  return PS;
}

// Exact 2-body partial decay width,
// is not dependent on momentum => phase space and matrix element factorize
//
// Input: Mother mass M0^2
//        Decay product 1 mass^2
//        Decay product 2 mass^2
//        Decay matrix element squared |M|^2 value
//        Final state symmetry factor S
//
// [REFERENCE: Alwall et al, https://arxiv.org/abs/1402.1178]
inline double PDW2body(double M0_2, double m1_2, double m2_2, double MatElem2, double S_factor) {
  return SqrtKallenLambda(M0_2, m1_2, m2_2) * MatElem2 /
         (16.0 * gra::math::PI * S_factor * math::pow3(math::msqrt(M0_2)));
}

// Relativistic velocity of a particle with m^2
// in the system with invariant shat
inline double Beta(double m2, double shat) { return gra::math::msqrt(1.0 - 4.0 * m2 / shat); }

// 2->2 Scattering E1 + E2 -> E3 + E4 angle as a function of (s,t,m_i^2)
// [REFERENCE: Kallen, Elementary Particle Physics, Addison-Wesley, 1964]
inline double CosthetaStar(double s, double t, double E1, double E2, double E3, double E4) {
  const double numer = s * s + s * (2.0 * t - E1 - E2 - E3 - E4) + (E1 - E2) * (E3 - E4);
  const double denom = SqrtKallenLambda(s, E1, E2) * SqrtKallenLambda(s, E3, E4);

  return numer / denom;
}

// 2->2 Scattering E1 + E2 -> E3 + E4 Mandelstam t-limits (physical boundaries)
// [REFERENCE: Kallen, Elementary Particle Physics, Addison-Wesley, 1964]
// [REFERENCE: Sjostrand, Mrenna, Skands, https://arxiv.org/abs/hep-ph/0603175]
inline void Two2TwoLimit(double s, double E1, double E2, double E3, double E4, double &tmin,
                         double &tmax) {
  // Triangle function
  const double KL12 = SqrtKallenLambda(s, E1, E2);
  const double KL34 = SqrtKallenLambda(s, E3, E4);

  // Geometry
  const double A = s - (E1 + E2 + E3 + E4) + (E1 - E2) * (E3 - E4) / s;
  const double B = KL12 * KL34 / s;
  const double C = (E3 - E1) * (E4 - E2) + (E1 + E4 - E2 - E3) * (E1 * E4 - E2 * E3) / s;
  tmin           = -0.5 * (A + B);
  tmax           = C / tmin;
}

// Standard conventions:
//
// E^2 = (mc^2)^2 + |pc|^2 = gamma*mc^2, p = gamma*m*v, m = constant!
// gamma = E/m = 1/sqrt(1-v^2/c^2) = 1/sqrt(1-beta^2)
// beta  = v/c = p/E => gamma*beta = p/m

//
// d^3p / (2E(2\pi)^3) = |\partial(p_x,p_y,p_z)/\partial(p_t,\phi,y)| dp_t d\phi
//   dy / (2E(2\pi)^3) = p_t^3 d_pt d\phi dy / (2E^2(2\pi)^3)
//

// A simple Monte Carlo weight (container)
class MCW {
 public:
  MCW() {
    W  = 0.0;
    W2 = 0.0;
    N  = 0.0;
  }
  MCW(double w, double w2, double n) {
    W  = w;
    W2 = w2;
    N  = n;
  }
  // + operator
  MCW operator+(const MCW &obj) {
    MCW res;
    res.W  = W + obj.W;
    res.W2 = W2 + obj.W2;
    res.N  = N + obj.N;
    return res;
  }
  // += operator
  MCW &operator+=(const MCW &obj) {
    this->W += obj.W;
    this->W2 += obj.W2;
    this->N += obj.N;
    return *this;  // return by reference
  }
  // * operator (does scale W and W^2 but not N)
  MCW operator*(double scale) { return MCW(GetW() * scale, GetW2() * scale * scale, GetN()); }

  // Push new weight
  void Push(double w) {
    W += w;
    W2 += w * w;
    N += 1.0;
  }
  // MC estimate of the integral
  double Integral() const { return (N > 0.0) ? W / N : 0.0; }
  // MC estimate of the integral error squared (standard error of the mean)
  double IntegralError2() const {
    if (N > 0.0) {
      return (W2 / N - gra::math::pow2(W / N)) / N;
    } else {
      return 0.0;
    }
  }
  double IntegralError() const { return gra::math::msqrt(IntegralError2()); }
  double GetW() const { return W; }
  double GetW2() const { return W2; }
  double GetN() const { return N; }
  void SetW(double w) { W = w; }
  void SetW2(double w2) { W2 = w2; }
  void SetN(double n) { N = n; }

 private:
  double W;   // sum of weights
  double W2;  // sum of weights^2
  double N;   // trials
};

// A weighted sum of MCW container integral values (use with VEGAS, for example)
class MCWSUM {
 public:
  MCWSUM() {}

  void Add(const MCW &x, double weight) {
    wsum += weight;
    wintsum += weight * x.Integral();
    error2sum += (weight * weight) * x.IntegralError2();
  }
  double Integral() const {
    if (wsum > 0.0) {
      return wintsum / wsum;
    } else {
      return 0.0;
    }
  }
  double IntegralError2() const {
    if (wsum > 0.0) {
      return error2sum / gra::math::pow2(wsum);
    } else {
      return 0.0;
    }
  }
  double IntegralError() const { return gra::math::msqrt(IntegralError2()); }

 private:
  // \sum_i weight_i
  double wsum = 0.0;

  // \sum_i weight_i integral_i
  double wintsum = 0.0;

  // \sum_i weight_i^2 error_i^2
  double error2sum = 0.0;
};

// Generic Lorentz boost into the direction given by 4-momentum of "boost"
//
// A. sign = -1 for boost to the rest frame (out of lab, for example)
// B. sign =  1 out from the rest frame (back to lab, for example)
//
template <typename T>
inline void LorentzBoost(const T &boost, double M0, T &p, int sign) {
  std::valarray<double> p3 = {p.Px(), p.Py(), p.Pz()};

  // Mother beta-vector
  const std::valarray<double> beta = {sign * boost.Px() / boost.E(), sign * boost.Py() / boost.E(),
                                      sign * boost.Pz() / boost.E()};
  // Mother gamma-factor
  const double gamma = boost.E() / M0;

  // Apply transforms
  const double kappa1 = (beta * p3).sum();  // Inner product
  const double kappa2 = gamma * (gamma * kappa1 / (1.0 + gamma) + p.E());

  // Vector sum
  p3 += kappa2 * beta;

  p = T(p3[0], p3[1], p3[2], gamma * (p.E() + kappa1));
}

// Flat variables in spherical coordinates
template <typename T2>
inline void FlatIsotropic(double &costheta, double &sintheta, double &phi, T2 &rng) {
  costheta = rng.U(-1.0, 1.0);  // cos(theta) flat [-1,1]
  sintheta = gra::math::msqrt(1.0 - gra::math::pow2(costheta));
  phi      = rng.U(0.0, 2.0 * gra::math::PI);  // phi flat [0,2pi]
}

// Isotropic decay 0 -> 1 + 2 in spherical coordinates with flat cos(theta),phi
template <typename T1, typename T2>
inline void Isotropic(double pnorm, T1 &p1, T1 &p2, double m1, double m2, T2 &rng) {
  using gra::math::pow2;
  using gra::math::msqrt;

  double costheta, sintheta, phi;
  FlatIsotropic(costheta, sintheta, phi, rng);

  // Jacobian of spherical coordinates
  const std::valarray<double> k = {pnorm * sintheta * std::cos(phi),
                                   pnorm * sintheta * std::sin(phi), pnorm * costheta};

  // Energies by on-shell condition
  const std::valarray<double> e = {msqrt(pow2(m1) + pow2(pnorm)), msqrt(pow2(m2) + pow2(pnorm))};

  // Back-to-back
  p1 = T1(k[0], k[1], k[2], e[0]);
  p2 = T1(-k[0], -k[1], -k[2], e[1]);
}

// 2-Body isotropic phase space decay in the rest frame
//
//           m0
//          /
//         /
// M0 ----*----- m1
//
// For floating point / efficiency reasons, the mother mass needs to be provided
// outside
//
template <typename T1, typename T2>
inline MCW TwoBodyPhaseSpace(const T1 &mother, double M0, const std::vector<double> &m,
                             std::vector<T1> &p, T2 &rng) {
  // Two particles
  p.resize(2);

  // First isotropic decay
  const double pnorm = DecayMomentum(M0, m[0], m[1]);
  Isotropic(pnorm, p[0], p[1], m[0], m[1], rng);

  // Boost daughters to the lab frame
  const int sign = 1;
  for (const auto &i : {0, 1}) { LorentzBoost(mother, M0, p[i], sign); }

  // return phase space weight
  const double weight = dPhi2(M0, pnorm);
  return MCW(weight, gra::math::pow2(weight), 1);
}

// 3-Body isotropic (Dalitz) phase space decay in the rest frame;
// by factorizing the 3-body phase space into two 2-body
//
//          m[0]       m[1]
//          /         /
//         /         /
// M0 ----*---m12---*----- m[2]
//
// dphi_3(pM; p0, p1, p2) ~ dM_12^2 dphi_2 (pM; p0, p12) dphi_2 (p12; p1, p2)
//
// Return: phase space weight for a valid fragmentation and -1.0 for a
// kinematically impossible.
// When unweight == true, This function returns also number of trials for proper
// MC error treatment
// (MCW)
//
template <typename T1, typename T2>
inline MCW ThreeBodyPhaseSpace(const T1 &mother, double M0, const std::vector<double> &m,
                               std::vector<T1> &p, bool unweight, T2 &rng) {
  p.resize(3);  // Three final states
  T1 p12;

  // Phase space boundaries [min,max]
  const std::vector<double> m12bound = {m[1] + m[2], M0 - m[0]};
  double w_max = DecayMomentum(M0, m[0], m12bound[0]) * DecayMomentum(m12bound[1], m[1], m[2]);
  if (unweight == false) { w_max = 0; }

  // Accceptance-Rejection
  const unsigned int  MAXTRIAL = 1e8;
  std::vector<double> pnorm    = {0.0, 0.0};
  double              m12      = 0;
  MCW                 x;
  do {
    m12      = rng.U(m12bound[0], m12bound[1]);  // Flat mass (in GeV, not GeV^2)
    pnorm[0] = DecayMomentum(M0, m[0], m12);
    pnorm[1] = DecayMomentum(m12, m[1], m[2]);

    const double w = (m12 / gra::math::PI) * dPhi2(M0, pnorm[0]) * dPhi2(m12, pnorm[1]);
    x.Push(w);
    if (x.GetN() > MAXTRIAL) { return MCW(-1, 0, 0); }  // Impossible kinematics
  } while ((pnorm[0] * pnorm[1]) < rng.U(0.0, w_max));

  // pM -> p0 + p12 in the pM r.f. and p12 -> p1 + p2 in the p12 r.f.
  Isotropic(pnorm[0], p[0], p12, m[0], m12, rng);
  Isotropic(pnorm[1], p[1], p[2], m[1], m[2], rng);

  // Boost p[1] and p[2] out from the p12 rest frame to the pM rest frame
  const int sign = 1;  // Boost sign
  for (const auto &i : {1, 2}) { LorentzBoost(p12, m12, p[i], sign); }

  // Boost p[0], p[1] and p[2] to the lab frame
  for (const auto &i : {0, 1, 2}) { LorentzBoost(mother, M0, p[i], sign); }

  // Phasespace weight [note m versus m^2 jacobian, taken into account]
  // [m12 volume (GeV)] x [dm_12] x [dphi2] x [dphi2]
  const double volume = (m12bound[1] - m12bound[0]);
  x                   = x * volume;  // operator overloaded

  // Return phase space weight
  return x;
}

// Decay setup for NBodyPhaseSpace
// Uses the "sorting algorithm", see description in F. James, CERN/68
// bool unweighted for (un)weighted operation
//
// Return: W for a valid and -1.0 for a kinematically impossible
//
template <typename T1, typename T2>
inline MCW NBodySetup(const T1 &mother, double M0, const std::vector<double> &m,
                      std::vector<double> &M_eff, std::vector<double> &pnorm, bool unweight,
                      T2 &rng) {
  // Decay multiplicity (1 -> N decay)
  const unsigned int  N = m.size();
  std::vector<double> randvec(N - 2, 0.0);

  // Random variables functor
  auto fillrandom = [&](double &value) -> void { value = rng.U(0, 1); };

  // Effective (intermediate) masses
  M_eff[0]     = m[0];  // First daughter
  M_eff[N - 1] = M0;    // Mother

  // Cumulative sum of daughter masses
  std::vector<double> sumvec;
  gra::math::CumSum(m, sumvec);

  // Sampling mass interval
  const double DELTA = M_eff[N - 1] - sumvec[N - 1];

  // Generate effective masses functor, where sumvec[i] defines the minimum
  auto calcmass = [&]() -> void {
    for (std::size_t i = 1; i < N - 1; ++i) { M_eff[i] = sumvec[i] + randvec[i - 1] * DELTA; }
  };

  // -------------------------------------------------------------------
  // Recursively find maximum weight for Acceptance-Rejection
  double              w_max    = 1.0;
  std::vector<double> m_minmax = {0.0, DELTA + m[N - 1]};
  for (std::size_t i = 1; i < N; ++i) {
    m_minmax[0] += m[i - 1];
    m_minmax[1] += m[i];
    w_max *= DecayMomentum(m_minmax[1], m_minmax[0], m[i]);
  }

  // Functor to calculate decay momentum to pnorm[], and return product
  auto calcmomentum = [&]() -> double {
    double w = 1.0;
    for (std::size_t i = 1; i < N; ++i) {
      pnorm[i] = DecayMomentum(M_eff[i], M_eff[i - 1], m[i]);
      w *= pnorm[i];
    }
    return w;
  };

  // -------------------------------------------------------------------
  // Acceptance-Rejection

  if (unweight == false) { w_max = 0; }

  const unsigned int MAXTRIAL = 1e8;
  MCW                x;
  double             w = 0.0;
  do {
    // 1. Ordered random numbers
    std::for_each(randvec.begin(), randvec.end(), fillrandom);
    std::sort(randvec.begin(), randvec.end());  // Ascending
    // std::sort(randvec.begin(), randvec.end(), [](const double a, const double
    // b) {
    // return a > b; }); // Descending

    // 2. Calculate effective masses
    calcmass();

    // 3. Calculate momentum norm values
    w = calcmomentum();

    x.Push(w);
    if (x.GetN() > MAXTRIAL) { return MCW(-1, 0, 0); }  // Impossible kinematics
  } while (w < rng.U(0.0, w_max));

  // Normalize the phase space integral
  const double volume = 1.0 / (2.0 * std::pow(2.0 * gra::math::PI, 2 * N - 3)) *
                        std::pow(DELTA, N - 2) / gra::math::factorial(N - 2) / M0;
  x = x * volume;  // operator overloaded

  // Return phase space weight
  return x;
}

// Flat N-body phase space by implementing a number of (N-2) recursive 2-body
// decays.
// Input as [mother 4-vector, vector of daughter masses, vector of final
// 4-vectors]
//
// [REFERENCE: F. James, https://cds.cern.ch/record/275743/files/CERN-68-15.pdf]
//
//                          m_n
//                         *             m_{n-1}
//  N-bodies              *             *              m_{n-2}
//       --------------- *             *              *
//       ---------------A------------ *              *
//  M_N  ----------------------------B------------- *
//       ------------------------------------------C
//                      M_{N-1}      M_{N-2}        *
//                                                   *
//                                                    *
//                                                     M_{N-3}
//
// For faster alternatives, investigate:
// M.M. Block, Monte Carlo phase space evaluation, Comp. Phys. Commun.
// 69, 459 (1992)
// S. Platzer, RAMBO on diet, https://arxiv.org/abs/1308.2922
//
// Return: weight for a valid fragmentation and -1.0 for a kinematically
// impossible
//
// When unweight == true, This function returns also number of trials for proper
// MC error treatment
// (MCW obj)
//
template <typename T1, typename T2>
inline MCW NBodyPhaseSpace(const T1 &mother, double M0, const std::vector<double> &m,
                           std::vector<T1> &p, bool unweight, T2 &rng) {
  using gra::math::pow2;
  using gra::math::msqrt;

  const int N = m.size();

  // Generate effective masses and decay momentum
  std::vector<double> M_eff(N, 0.0);  // Effective masses
  std::vector<double> pnorm(N, 0.0);  // Decay momentum weights
  const MCW           x = NBodySetup(mother, M0, m, M_eff, pnorm, unweight, rng);
  if (x.GetW() < 0) { return x; }  // Impossible kinematics

  // Setup decay daughters size
  p.resize(N);

  // Recursively go down from N-body to 1-particle
  T1 k = mother;
  for (std::size_t i = N - 1; i >= 1; --i) {
    double costheta, sintheta, phi;
    FlatIsotropic(costheta, sintheta, phi, rng);
    const double pt = pnorm[i] * sintheta;

    p[i] = T1(pt * std::cos(phi), pt * std::sin(phi), pnorm[i] * costheta,
              msqrt(pow2(m[i]) + pow2(pnorm[i])));

    // Boost to the lab frame
    const int sign = 1;  // positive
    LorentzBoost(k, k.M(), p[i], sign);

    // Subtract the momentum
    k -= p[i];
  }
  // Finally, we are left with only one daughter
  p[0] = k;

  return x;  // Phase-space weight
}

// SO(3) Rotations

// Rotation matrix from 3-vector a to 3-vector b (both need to be unit vectors with ||x|| = 1)
//
inline MMatrix<double> RotOnetoOne(const std::vector<double>& a, const std::vector<double>& b) {

  if (a.size() != 3 || b.size() != 3) {
    throw std::invalid_argument("kinematics::RotOnetoOne: Input should be 3-vectors");
  }
  // Dot product
  const double c = matoper::VecVecMultiply(a,b);

  // Cross product
  const std::vector<double> v = matoper::Cross(a,b);

  // Skew symmetric cross product matrix
  const MMatrix<double> Vx = {{ 0,    -v[2],  v[1]},
                              { v[2],     0, -v[0]},
                              {-v[1],  v[0],    0}};
  // Identity 
  const MMatrix<double> I(3,3,"eye");

  // Rotation matrix, with singularity at c = -1 (a and b back-to-back)
  const MMatrix<double> R = I + Vx + (Vx*Vx) * (1.0/(1.0 + c));

  return R;
}

// Rotate 4-vector spatial part by (theta,phi)
template <typename T>
inline void Rotate(T &p, double theta, double phi) {
  const double c1 = std::cos(theta);
  const double s1 = std::sin(theta);
  const double c2 = std::cos(phi);
  const double s2 = std::sin(phi);

  // 3x3 Rotation matrix
  const MMatrix<double> R = {{c1 * c2, -s2, s1 * c2}, {c1 * s2, c2, s1 * s2}, {-s1, 0.0, c1}};

  // Rotate
  const std::vector<double> p3    = {p.Px(), p.Py(), p.Pz()};
  const std::vector<double> p3new = R * p3;

  p = T(p3new[0], p3new[1], p3new[2], p.E());
}

// Active SO3-rotation counterclockwise around x-axis
template <typename T>
inline void RotateX(T &p, double angle) {
  const double costh = std::cos(angle);
  const double sinth = std::sin(angle);

  // Rotation matrix applied
  const double px = p.Px();
  const double py = costh * p.Py() - sinth * p.Pz();
  const double pz = sinth * p.Py() + costh * p.Pz();

  p = T(px, py, pz, p.E());
}

// Active SO3-rotation counterclockwise around y-axis
template <typename T>
inline void RotateY(T &p, double angle) {
  const double costh = std::cos(angle);
  const double sinth = std::sin(angle);

  // Rotation matrix applied
  const double px = costh * p.Px() + sinth * p.Pz();
  const double py = p.Py();
  const double pz = -sinth * p.Px() + costh * p.Pz();

  p = T(px, py, pz, p.E());
}

// Active SO3-rotation counterclockwise around z-axis
template <typename T>
inline void RotateZ(T &p, double angle) {
  const double costh = std::cos(angle);
  const double sinth = std::sin(angle);

  // Rotation matrix applied
  const double px = costh * p.Px() - sinth * p.Py();
  const double py = sinth * p.Px() + costh * p.Py();
  const double pz = p.Pz();

  p = T(px, py, pz, p.E());
}


// Transform off-shell to on-shell kinematics p1 + p2 -> {p} ...
// with optimization criteria being exact momentum conservation
// after the initial states have been put on-shell.
// 
// where
// p1, p2 are massless originally off-shell -> changed to on-shell
// {p} massive/massless final states
//
inline void OffShell2OnShell(M4Vec& p1, M4Vec& p2, std::vector<M4Vec>& p) {
  
  const int N = p.size();  
  const int MAXITER = 15;
  const double STOPEPS = 1e-10;

  auto psumfunc = [&] () {
    M4Vec sum(0,0,0,0);
    for (const auto& i : aux::indices(p)) { sum += p[i]; }
    return sum;
  };

  auto Esumfunc = [&] () {
    double sum = 0.0;
    for (const auto& i : aux::indices(p)) { sum += p[i].E(); }
    return sum;
  };

  // Set massless particles to on-shell (px,py,pz,E)
  p1.Set(p1.Px(), p1.Py(), p1.Pz(), p1.P3mod());
  p2.Set(p2.Px(), p2.Py(), p2.Pz(), p2.P3mod());

  // Get sum
  const M4Vec q = p1 + p2;
  
  // Normalize \vec{q} to unit length
  const std::vector<double> qvec = gra::matoper::Unit(q.P3());

  // Final state masses
  std::vector<double> m(p.size());
  for (const auto& i : aux::indices(p)) { m[i] = p[i].M(); }

  int iter = 0;
  
  //for (const auto& i : aux::indices(p)) { p[i].Print(); }
  //std::cout << "iter =" << iter << "  " << (q - psumfunc()).M2() << std::endl;
  
  while (true) {

    // --------------------------------------------------
    // 3-Momentum scaling and energy subtraction step

    // New energy difference
    const double dE = Esumfunc() - q.E();

    for (const auto& i : aux::indices(p)) {
      const double E = p[i].E() - dE / N; // Scaling distributed by 1/N
      const double a = gra::math::msqrt(E*E - m[i]*m[i]) / p[i].P3mod();  
      p[i].Set(a*p[i].Px(), a*p[i].Py(), a*p[i].Pz(), E);
    }

    // --------------------------------------------------
    // 3-Momentum subtraction step

    // New 3-momentum difference
    std::vector<double> D3 = (psumfunc() - q).P3();
    gra::matoper::ScaleVector(D3, 1.0/N); // Scaling distributed by 1/N

    for (const auto& i : aux::indices(p)) {
      p[i].SetP3(gra::matoper::Minus(p[i].P3(), D3));
      p[i].SetE(gra::math::msqrt(p[i].P3mod2() + m[i]*m[i]));
    }

    // --------------------------------------------------
    // 3-Momentum rotation step

    // Add and normalize to unit length, find rotation
    const std::vector<double> avec = gra::matoper::Unit(psumfunc().P3());
    const gra::MMatrix<double> R = gra::kinematics::RotOnetoOne(avec, qvec);

    // Rotate
    for (const auto& i : aux::indices(p)) { p[i].SetP3( R * p[i].P3()); }

    // --------------------------------------------------

    ++iter;
    const double IVM = std::abs( (q - (psumfunc())).M2() );
    if (iter > MAXITER || IVM < STOPEPS ) { break; }

    //std::cout << "iter =" << iter << "  " << IVM << std::endl;
  }

  //for (const auto& i : aux::indices(p)) { p[i].Print(); }
  //std::cout << " --------------------------------------------------- " << std::endl;
  
}


// A unified Lorentz Transform function to 'frametype' give below:
//
// "CS" : Collins-Soper
// "AH" : Anti-Helicity (Anti-CS)
// "HE" : Helicity
// "PG" : Pseudo-Gottfried-Jackson
//
// Collins-Soper: Quantization z-axis defined by the bi-sector vector between
// initial state
// p1 and (-p2) (NEGATIVE) directions in the (resonance) system rest frame,
// where p1 and p2 are the initial state proton 3-momentum.
//
// Anti-Helicity: Quantization z-axis defined by the bisector vector between
// initial state
// p1 and (p2) (POSITIVE) directions in the (resonance) system rest frame,
// where p1 and p2 are the initial state proton 3-momentum.
//
// Helicity: Quantization axis defined by the resonance system 3-momentum vector
// in
// the colliding beams frame (lab frame).
//
// Pseudo-Gottfried-Jackson: Quantization axis defined by the initial state
// proton
// p1 (or p2) 3-momentum vector in the (resonance) system rest frame.
//

template <typename T>
inline void LorentFramePrepare(const std::vector<T> &pf, const T &pbeam1, const T &pbeam2,
                               T &pb1boost, T &pb2boost, std::vector<T> &pfboost) {
  // 1. Central system 4-momentum as a sum
  T X(0, 0, 0, 0);
  for (const auto &k : gra::aux::indices(pf)) { X += pf[k]; }
  const double     M = X.M();

  // 2. Boost each final state to the system rest frame
  pfboost = pf;
  for (const auto &k : gra::aux::indices(pf)) {
    gra::kinematics::LorentzBoost(X, M, pfboost[k], -1);  // note minus sign
  }

  // 3. Boost initial state protons
  pb1boost = pbeam1;
  pb2boost = pbeam2;
  gra::kinematics::LorentzBoost(X, M, pb1boost, -1);  // note minus sign
  gra::kinematics::LorentzBoost(X, M, pb2boost, -1);  // note minus sign
}

template <typename T>
inline void LorentzFrame(std::vector<T> &pfout, const T &pb1boost, const T &pb2boost,
                         const std::vector<T> &pfboost, const std::string &frametype,
                         int direction) {
  const std::vector<double> pb1boost3 = {pb1boost.Px(), pb1boost.Py(), pb1boost.Pz()};
  const std::vector<double> pb2boost3 = {pb2boost.Px(), pb2boost.Py(), pb2boost.Pz()};

  // Frame rotation x-y-z-axes
  std::vector<double> zaxis;
  std::vector<double> yaxis;
  std::vector<double> xaxis;

  // @@ NON-ROTATED FRAME AXIS DEFINITION @@
  if (frametype == "SR") {
    zaxis = {0, 0, 1};
    yaxis = {0, 1, 0};
  }
  // @@ COLLINS-SOPER FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "CS") {
    zaxis = gra::matoper::Unit(
        gra::matoper::Minus(gra::matoper::Unit(pb1boost3), gra::matoper::Unit(pb2boost3)));
  }
  // @@ ANTI-HELICITY FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "AH") {
    zaxis = gra::matoper::Unit(
        gra::matoper::Plus(gra::matoper::Unit(pb1boost3), gra::matoper::Unit(pb2boost3)));
  }
  // @@ HELICITY FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "HE") {
    zaxis = gra::matoper::Unit(gra::matoper::Negat(gra::matoper::Plus(pb1boost3, pb2boost3)));
  }
  // @@ PSEUDO-GOTTFRIED-JACKSON AXIS DEFINITION: [1] or [2] @@
  else if (frametype == "PG") {
    if (direction == -1) {
      zaxis = gra::matoper::Unit(pb1boost3);
    } else if (direction == 1) {
      zaxis = gra::matoper::Unit(pb2boost3);
    } else {
      throw std::invalid_argument("gra::kinematics::LorentzFrame: Invalid direction <" +
                                  std::to_string(direction) + ">");
    }
  } else {
    throw std::invalid_argument("gra::kinematics::LorentzFrame: Unknown frame <" + frametype + ">");
  }

  // y-axis
  if (frametype != "SR") {
    yaxis = gra::matoper::Unit(
        gra::matoper::Cross(gra::matoper::Unit(pb2boost3), gra::matoper::Unit(pb1boost3)));
  }

  // x-axis
  xaxis = gra::matoper::Unit(gra::matoper::Cross(yaxis, zaxis));  // x = y [cross product] z

  // Create SO(3) rotation matrix for the new coordinate axes
  const MMatrix<double> R = {xaxis, yaxis, zaxis};  // Axes as rows

  // Rotate all vectors
  pfout = pfboost;
  for (const auto &k : gra::aux::indices(pfout)) {
    const std::vector<double> p3    = {pfout[k].Px(), pfout[k].Py(), pfout[k].Pz()};
    const std::vector<double> p3new = R * p3;                  // Spatial part rotation Rp -> p'
    pfout[k] = T(p3new[0], p3new[1], p3new[2], pfout[k].E());  // Full 4-momentum [px; py; pz; E]
  }
}

// ---------------------------------------------------------------------------
const bool   DEBUG   = false;
const double epsilon = 1e-8;

// From lab to Collins-Soper frame
// Quantization z-axis is defined as the bisector of two beams, in the rest
// frame of the resonance
template <typename T>
inline void CSframe(std::vector<T> &p) {
  // Sum to get the system 4-momentum
  T X(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { X += p[i]; }

  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::COLLINS-SOPER FRAME:: \n");
    printf("CSframe:: Pions in LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  const double Z_angle = -X.Phi();    // Note minus
  const double Y_angle = -X.Theta();  // Note minus

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::RotateZ(p[i], Z_angle);
    gra::kinematics::RotateY(p[i], -Y_angle);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("CSframe:: Pions in ROTATED LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  // Construct new central system vector
  T XNEW(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { XNEW += p[i]; }

  // Lorentz boost
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(XNEW, XNEW.M(), p[i], -1);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("CSframe:: Pions after boost in Collins-Soper FRAME: \n");
    p[0].Print();
    p[1].Print();

    printf("CSframe:: Central system in LAB FRAME: \n");
    X.Print();

    printf("CSframe:: Central system in ROTATED LAB FRAME: \n");
    XNEW.Print();

    printf("\n");
  }
  // ********************************************************************

  T sum;
  for (const auto &i : gra::aux::indices(p)) { sum += p[i]; }
  if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon ||
      std::abs(sum.Pz()) > epsilon) {
    printf("CSframe:: Not a rest frame!! \n");
  }
}

// From lab to Gottfried-Jackson frame
// Quantization z-axis spanned by the propagator (Pomeron, Gamma etc.) momentum
// in the rest frame of the central system
template <typename T>
inline void GJframe(std::vector<T> &p, const T &q) {
  // Take the lab frame propagator, we will boost and rotate this
  T propagator = q;

  // Sum to get the system 4-momentum
  T X(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { X += p[i]; }

  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::GOTTFRIED-JACKSON FRAME:: \n");

    printf("GJframe:: Pions in LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  // Boost particles to the central system rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, X.M(), p[i], -1);  // Note the minus sign
  }

  // Boost the propagator (Pomeron)
  gra::kinematics::LorentzBoost(X, X.M(), propagator, -1);  // Note the minus sign

  // Now get the rotation angle, note the minus
  const double Z_angle = -propagator.Phi();
  const double Y_angle = -propagator.Theta();

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::RotateZ(p[i], Z_angle);
    gra::kinematics::RotateY(p[i], Y_angle);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("GJframe:: Pions in Gottfried-Jackson FRAME: \n");
    p[0].Print();
    p[1].Print();

    // Rotate propagator
    gra::kinematics::RotateZ(propagator, Z_angle);
    gra::kinematics::RotateY(propagator, Y_angle);

    printf("GJframe:: Propagator in Gottfried-Jackson FRAME: \n");
    propagator.Print();
  }
  // ********************************************************************

  T sum;
  for (const auto &i : gra::aux::indices(p)) { sum += p[i]; }
  if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon ||
      std::abs(sum.Pz()) > epsilon) {
    printf("GJframe:: Not a rest frame!! \n");
  }
}

// From lab to the Helicity frame
// Quantization z-axis as the direction of the resonance in the lab frame
//
// Input is a vector of final state 4-momentum
template <typename T>
inline void HEframe(std::vector<T> &p) {
  // Sum to get the system 4-momentum
  T X(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { X += p[i]; }

  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::HELICITY FRAME:: \n");
    printf("HEframe:: Pions in LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  // ACTIVE (cf. PAsinthIVE) rotation of final states by z-y-z
  // (phi,theta,-phi)
  // (in the opposite direction -> minus signs) Euler angle sequence.
  const double Z_angle = -X.Phi();
  const double Y_angle = -X.Theta();

  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::RotateZ(p[i], Z_angle);
    gra::kinematics::RotateY(p[i], Y_angle);
    //    gra::kinematics::RotateZ(p[i], -Z_angle); // Comment this one
    //    out for tests
  }

  // ********************************************************************
  if (DEBUG) {
    printf("HEframe:: Pions in ROTATED LAB FRAME: \n");
    p[0].Print();
    p[1].Print();

    T ex(1, 0, 0, 0);
    T ey(0, 1, 0, 0);
    T ez(0, 0, 1, 0);

    // x -> x'
    gra::kinematics::RotateZ(ex, Z_angle);
    gra::kinematics::RotateY(ex, Y_angle);
    // gra::kinematics::RotateZ(ex, -Z_angle);  // Comment this one out
    // for tests
    // y -> y'
    gra::kinematics::RotateZ(ey, Z_angle);
    gra::kinematics::RotateY(ey, Y_angle);
    //   gra::kinematics::RotateZ(-Z_angle);  // Comment this one out
    //   for tests
    // z -> z'
    gra::kinematics::RotateZ(ez, Z_angle);
    gra::kinematics::RotateY(ez, Y_angle);
    //     gra::kinematics::RotateZ(ez, -Z_angle);  // Comment this one
    //     out for tests

    printf("HEframe:: AXIS vectors after rotation: \n");

    ex.Print();
    ey.Print();
    ez.Print();
  }
  // ********************************************************************

  // Construct the central system 4-momentum in ROTATED FRAME
  T XNEW(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { XNEW += p[i]; }

  // Boost particles to the central system rest frame
  // -> Helicity frame obtained
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(XNEW, XNEW.M(), p[i], -1);  // Note the minus sign
  }

  // ********************************************************************
  if (DEBUG) {
    printf("HEframe:: Pions after boost in HELICITY FRAME: \n");
    p[0].Print();
    p[1].Print();

    printf("HEframe:: Central system in LAB FRAME: \n");
    X.Print();

    printf("HEframe:: Central system in ROTATED LAB FRAME: \n");
    XNEW.Print();

    printf("\n");
  }
  // ********************************************************************

  T sum;
  for (const auto &i : gra::aux::indices(p)) { sum += p[i]; }
  if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon ||
      std::abs(sum.Pz()) > epsilon) {
    printf("HEframe:: Not a rest frame!! \n");
  }
}

// From lab to Pseudo-Gottfried-Jackson frame
// Quantization z-axis spanned by the beam proton +z (or -z) momentum
// in the rest frame of the central system
//
// Direction (1 = +z direction proton, -1 = -z direction proton)
template <typename T>
inline void PGframe(std::vector<T> &p, const int direction, const T &p_beam_plus,
                    const T &p_beam_minus) {
  // Sum to get the system 4-momentum
  T X(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { X += p[i]; }

  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::PSEUDO-GOTTFRIED-JACKSON FRAME:: \n");

    printf("PGframe:: Pions in LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  // Boost particles to the central system rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, X.M(), p[i], -1);  // Note the minus sign
  }

  // Boost initial state protons to the central system rest frame
  T proton_p = p_beam_plus;
  T proton_m = p_beam_minus;

  gra::kinematics::LorentzBoost(X, X.M(), proton_p, -1);  // Note the minus sign
  gra::kinematics::LorentzBoost(X, X.M(), proton_m, -1);  // Note the minus sign

  // Now get the rotation angle, note the minus
  double Z_angle = 0;
  double Y_angle = 0;

  // Choose which proton we use
  if (direction == 1) {
    Z_angle = -proton_p.Phi();
    Y_angle = -proton_p.Theta();
  } else if (direction == -1) {
    Z_angle = -proton_m.Phi();
    Y_angle = -proton_m.Theta();
  } else {
    printf("PGframe:: Input direction %d not valid (1,-1)!", direction);
    return;
  }

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::RotateZ(p[i], Z_angle);
    gra::kinematics::RotateY(p[i], Y_angle);
    //    gra::kinematics::RotateZ(p[i], -Z_angle);  // Comment this
    //    one out for tests
  }

  // ********************************************************************
  if (DEBUG) {
    printf("PGframe:: Pions in Pseudo-Gottfried-Jackson FRAME: \n");
    p[0].Print();
    p[1].Print();

    // Rotate proton
    gra::kinematics::RotateZ(proton_p, Z_angle);
    gra::kinematics::RotateY(proton_p, Y_angle);
    //   gra::kinematics::RotateZ(proton_p, -Z_angle);  // Comment this
    //   one out for tests

    printf("Protons are on (zx)-plane: \n");
    printf("USER chosen direction along %d beam direction \n", direction);
    printf("- +z Proton in Pseudo-Gottfried-Jackson FRAME: \n");
    proton_p.Print();

    // Rotate other proton
    gra::kinematics::RotateZ(proton_m, Z_angle);
    gra::kinematics::RotateY(proton_m, Y_angle);
    //   gra::kinematics::RotateZ(proton_m, -Z_angle);

    printf("PGframe:: -z Proton in Pseudo-Gottfried-Jackson FRAME: \n");
    proton_m.Print();

    // Mandelstam invariant
    printf("Mandelstam s = %0.5f \n", (proton_p + proton_m).M());
  }
  // ********************************************************************
  T sum;
  for (const auto &i : gra::aux::indices(p)) { sum += p[i]; }
  if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon ||
      std::abs(sum.Pz()) > epsilon) {
    printf("PGframe:: Not a rest frame!! \n");
  }
}

// "Rest frame" (no boost, beam axis as z-axis / spin quantization
// axis)
//
// Input is a vector of final state 4-momentum
template <typename T>
inline void SRframe(std::vector<T> &p) {
  // Sum to get the system 4-momentum
  T X(0, 0, 0, 0);
  for (const auto &i : gra::aux::indices(p)) { X += p[i]; }

  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::STANDARD REST FRAME:: \n");

    printf("- Pions in LAB FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************

  // Boost particles to the central system rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, X.M(), p[i], -1);  // Note the minus sign
  }

  // ********************************************************************
  if (DEBUG) {
    printf("- Pions in REST FRAME: \n");
    p[0].Print();
    p[1].Print();
  }
  // ********************************************************************
  T sum;
  for (const auto &i : gra::aux::indices(p)) { sum += p[i]; }
  if (std::abs(sum.Px()) > epsilon || std::abs(sum.Py()) > epsilon ||
      std::abs(sum.Pz()) > epsilon) {
    printf("SRframe:: Not a rest frame!! \n");
  }
}

}  // Namespace MKinematics

// Particle class
struct MParticle {
  // These values are read from the PDG-file
  std::string name;
  int         pdg      = 0;  // PDG code
  int         chargeX3 = 0;  // Q x 3
  int         spinX2   = 0;  // J x 2
  int         color    = 0;  // QCD color code

  double mass  = 0.0;
  double width = 0.0;
  double tau   = 0.0;  // hbar / width

  // J^PC
  int          P    = 1;      // Default even P-parity
  int          C    = 0;      // Default zero C-parity (e.g. fermions do not have)
  unsigned int L    = 0;      // For Mesons/Baryons
  bool         glue = false;  // Glueball state

  void setPCL(int _P, int _C, unsigned int _L) {
    P = _P;
    C = _C;
    L = _L;
  }

  // Width cut (in Breit-Wigner sampling)
  double wcut = 0.0;

  void print() const {
    std::cout << " NAME:   " << name << std::endl;
    std::cout << " ID:     " << pdg << std::endl;
    std::cout << " M:      " << mass << std::endl;
    std::cout << " W:      " << width << std::endl;
    std::cout << " J^PC:   " << aux::Spin2XtoString(spinX2) << "^" << aux::ParityToString(P)
              << aux::ParityToString(C) << std::endl
              << std::endl;
  }
};

// Recursive decay tree branch
struct MDecayBranch {
  // Offshell mass picked event by event
  double m_offshell = 0.0;

  MParticle                 p;               // PDG particle
  M4Vec                     p4;              // 4-momentum
  std::vector<MDecayBranch> legs;            // Daughters
  M4Vec                     decay_position;  // Decay 4-position

  // MC weight container
  gra::kinematics::MCW W;

  // Decay tree current level
  int depth = 0;
};

struct HELMatrix {
  void InitAlphaToZero() {
    const unsigned int N = 20;
    alpha                = MMatrix<std::complex<double>>(N, N, 0.0);
    alpha_set            = MMatrix<bool>(N, N, true);
  }

  // Lorentz frame of the spin distribution
  std::string FRAME = "";

  // Spin-density matrix (constant)
  MMatrix<std::complex<double>> rho;

  // Decay amplitude matrix (constant)
  MMatrix<std::complex<double>> T;

  // Helicity amplitude decay ls-couplings (constant)
  MMatrix<std::complex<double>> alpha;
  MMatrix<bool>                 alpha_set;

  // Parity conservation
  bool P_conservation = true;
};

// Resonance parameters
class PARAM_RES {
 public:
  PARAM_RES() {}
  void PrintParam(double sqrts) const {
    std::cout << rang::style::bold << "Custom resonance parameters:" << rang::style::reset
              << std::endl
              << std::endl;

    printf("- PDG ID:      %d \n", p.pdg);
    printf("- Mass M0:     %0.5f GeV \n", p.mass);
    printf("- Width W:     %0.5f GeV \n", p.width);
    printf("- Spin Jx2:    %d \n", p.spinX2);
    printf("- Parity P:    %d \n", p.P);
    std::cout << std::endl;

    printf("<Production> \n");
    printf("- Effective vertex constant g_i:  %0.1E x exp(i x %0.1f) \n", std::abs(g), std::arg(g));
    printf("- Form factor: %s \n", g_FF == true ? "true" : "false");

    std::cout << std::endl;
    printf("<Decay> \n");
    printf("- BR:          %0.3E\n", BR);
    printf("- Effective vertex constant g_f:  %0.3E", g_decay);
    std::cout << std::endl << std::endl;

    if (p.spinX2 != 0) {
      std::cout << "- Polarization Lorentz frame: " << hc.FRAME << std::endl;
      std::cout << std::endl;
      std::cout << "- Spin polarization density matrix [rho]:" << std::endl << std::endl;

      // Print elements
      for (std::size_t i = 0; i < hc.rho.size_row(); ++i) {
        for (std::size_t j = 0; j < hc.rho.size_col(); ++j) {
          std::string delim = (j < hc.rho.size_col() - 1) ? ", " : "";
          printf("%6.3f+i%6.3f%s ", std::real(hc.rho[i][j]), std::imag(hc.rho[i][j]),
                 delim.c_str());
        }
        std::cout << std::endl;
      }
      std::cout << std::endl << std::endl;
      // printf("Density matrix von Neumann entropy S = %0.3f
      // \n\n", MSpin::VonNeumannEntropy(rho));
    }
  }

  // Particle class
  MParticle p;

  // -------------------------------------------------------------------
  // Tensor Pomeron couplings
  std::vector<double> g_Tensor;
  // -------------------------------------------------------------------

  // (Complex) production coupling constant
  std::complex<double> g = 0.0;
  
  // Form factor
  bool g_FF = false;

  // Breit-Wigner type
  int BW = 0;

  // Branching Ratio to a particular decay
  double BR      = 1.0;
  double g_decay = 1.0;  // Equivalent decay coupling

  // --------------------------------------------------------------------
  // Helicity amplitude matrix
  HELMatrix hc;
  // --------------------------------------------------------------------
};

// Lorentz scalars and other common kinematic variables
class LORENTZSCALAR {

public:

  LORENTZSCALAR() {}
  ~LORENTZSCALAR() {
      delete GlobalSudakovPtr;
      delete GlobalPdfPtr;
  }
  
  // --------------------------------------------------------------------
  // PDF access
  // Sudakov/pdf routines
  MSudakov* GlobalSudakovPtr = nullptr;

  // Normal pdfs
  std::string  LHAPDFSET = "null";
  LHAPDF::PDF* GlobalPdfPtr = nullptr;
  int pdf_trials = 0;
  // --------------------------------------------------------------------


  // Resonances read from JSON input
  std::map<std::string, PARAM_RES> RESONANCES;

  // Helicity amplitudes returned by amplitude functions,
  // used in screening loop etc.
  std::vector<std::complex<double>> hamp;

  // Central system decaytree
  std::vector<MDecayBranch> decaytree;

  // Forward system decaytree
  MDecayBranch decayforward1;
  MDecayBranch decayforward2;

  // Integral weight container (central system phase space)
  gra::kinematics::MCW DW;

  // Sum containers
  gra::kinematics::MCWSUM DW_sum;
  gra::kinematics::MCWSUM DW_sum_exact;

  // Four momenta of initial and final states
  M4Vec              pbeam1;
  M4Vec              pbeam2;
  std::vector<M4Vec> pfinal;
  std::vector<M4Vec> pfinal_orig;

  // Basic Lorentz scalars
  double s      = 0.0;
  double t      = 0.0;
  double u      = 0.0;
  double sqrt_s = 0.0;

  // Sub-energies^2
  double s1 = 0.0;
  double s2 = 0.0;

  // 4-momentum transfer squared
  double t1 = 0.0;
  double t2 = 0.0;

  // Sub-Mandelstam variables
  double s_hat = 0.0;
  double t_hat = 0.0;
  double u_hat = 0.0;

  // Central system
  double m2 = 0.0;
  double Y  = 0.0;
  double Pt = 0.0;

  // Maximum 10 particles in the final state
  // (2 protons + 8 direct central system)
  double ss[10][10]    = {{0.0}};
  double tt_1[10]      = {0.0};
  double tt_2[10]      = {0.0};
  double tt_xy[10][10] = {{0.0}};

  // Longitudinal momentum losses
  double x1 = 0.0;
  double x2 = 0.0;

  // Bjorken-x
  double xbj1 = 0.0;
  double xbj2 = 0.0;

  // Forward proton excitation tag
  bool excite1 = false;
  bool excite2 = false;

  // Propagators from proton1 (up) and proton2 (down)
  M4Vec q1;
  M4Vec q2;

  // Propagator pt
  double qt1 = 0.0;
  double qt2 = 0.0;
};

// Generator cut (default parameters set here)
struct GENCUT {
  // Continuum phase space <C>
  double rap_min = -9.0;
  double rap_max = 9.0;
  double kt_min  = 0.0;
  double kt_max  = -1.0;  // Keep at -1 for user setup trigger

  // Factorized phase space <F>
  double Y_min = -9.0;
  double Y_max = 9.0;
  double M_min = 0.0;
  double M_max = 0.0;

  // Both <C> and <F> class forward legs
  double forward_pt_min = -1.0;  // Keep at -1 for user setup trigger
  double forward_pt_max = -1.0;  // Keep at -1 for user setup trigger

  // ---------------------------------------

  // Quasi-Elastic phase space <Q> or forward excitation
  double XI_min = 0.0;
  double XI_max = 1.0;
};

// Fiducial cuts (default parameters set here)
struct FIDCUT {
  bool active = false;

  // "Central particles"
  double eta_min = -30.0;
  double eta_max = 30.0;

  double rap_min = -30.0;
  double rap_max = 30.0;

  double pt_min = 0.0;
  double pt_max = 1000000.0;

  double Et_min = 0.0;
  double Et_max = 1000000.0;

  // "Central system"
  double M_min = 0.0;
  double M_max = 1000000.0;

  double Y_min = -30.0;
  double Y_max = 30.0;

  double Pt_min = 0.0;
  double Pt_max = 1000000.0;

  // "Forward system"
  double forward_t_min = 0.0;
  double forward_t_max = 1000000.0;

  double forward_M_min = 0.0;
  double forward_M_max = 1000000.0;
};

// For veto
struct VETODOMAIN {
  double eta_min = -30.0;
  double eta_max = 30.0;

  double pt_min = 0.0;
  double pt_max = 1000000.0;
};

struct VETOCUT {
  bool                    active = false;
  std::vector<VETODOMAIN> cuts;
};



}  // gra namespace

#endif