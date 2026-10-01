// Standard relativistic kinematics [HEADER ONLY file]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MKINEMATICS_H
#define MKINEMATICS_H

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <valarray>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Sampling/MCW.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace kinematics {

// Compute the inverse density of logarithmic transverse Helmert-ball sampling
inline double TransverseLogJacobian(double rho2, double radius2, double scale2, unsigned int multiplicity) {
  const double tolerance = 128.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, radius2);
  if (rho2 < 0.0 || rho2 > radius2 + tolerance) { return 0.0; }
  const unsigned int modes = multiplicity - 1;
  return multiplicity * std::pow(math::PI, modes) * std::pow(rho2, modes - 1) *
         std::log1p(radius2 / scale2) * (scale2 + rho2) / std::tgamma(static_cast<double>(modes));
}


constexpr double KINEMATICS_EPS = 1e-12;

// Compute the standard rejected phase-space weight
inline MCW InvalidPhaseSpacePoint() { return MCW(-1.0, 0.0, 0.0); }

// Compute whether M0 exceeds the sum of all daughter masses
// M0 > sum_i m_i
inline bool HasPhaseSpace(double M0, const std::vector<double> &m) {
  if (!std::isfinite(M0) || M0 <= 0.0) {
    return false;
  }

  double sum = 0.0;
  for (const auto &mass : m) {
    if (!std::isfinite(mass) || mass < 0.0) {
      return false;
    }
    sum += mass;
  }
  return (M0 - sum) > 0.0;
}

// Case m3 = m4
// Solve sqrt(s) = E3 + E4 + E5 and pz3 + pz4 + pz5 = 0 on shell
inline double SolvePz3_A(double m3, double pt3, double pt4, double pz5,
                         double E5, double s) {
  if (!(s > 0.0) || !std::isfinite(s)) {
    return -1.0;
  }
  const double t2 = E5 * E5;
  const double t3 = pt3 * pt3;
  const double t4 = pt4 * pt4;
  const double t5 = pz5 * pz5;
  const double t6 = m3 * m3;
  const double t7 = std::sqrt(s);
  const double t8 = s * t2 * 6.0;
  const double t9 = std::pow(s, 3.0 / 2.0);
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
  const double t22 = t8 + t10 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 +
                     t19 + t20 + t21 - E5 * t9 * 4.0 - s * t3 * 2.0 -
                     s * t4 * 2.0 - s * t5 * 2.0 - s * t6 * 4.0 -
                     t2 * t3 * 2.0 - t2 * t4 * 2.0 - t2 * t5 * 2.0 -
                     t3 * t4 * 2.0 - t2 * t6 * 4.0 - E5 * t2 * t7 * 4.0;
  if (t22 < -KINEMATICS_EPS) {
    return -1.0;
  }
  const double denom = s + t2 - t5 - E5 * t7 * 2.0;
  if (std::abs(denom) < KINEMATICS_EPS) {
    return -1.0;
  }
  const double t23 = gra::math::msqrt(t22);
  const double t0 =
      pz5 * (-1.0 / 2.0) - (E5 * t23 * (1.0 / 2.0) - t7 * t23 * (1.0 / 2.0) +
                            pz5 * (t3 * (1.0 / 2.0) - t4 * (1.0 / 2.0))) /
                               denom;

  return std::isfinite(t0) ? t0 : -1.0;
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
                         2.0 * pow2(pE) * pow2(pt2) - 2.0 * pow2(pt1) *
pow2(pt2) + pow4(pt2) - 2.0 * pow2(pE) * pow2(pz) + 2.0 * pow2(pt1) * pow2(pz)
+ 2.0 * pow2(pt2) * pow2(pz) + pow4(pz) - 4.0 * pow3(pE) * sqrt_s + 4.0 * pE *
pow2(pt1) * sqrt_s + 4.0 * pE * pow2(pt2) * sqrt_s + 4.0 * pE * pow2(pz) *
sqrt_s + 6.0 * pow2(pE) * s - 2.0 * pow2(pt1) * s - 2.0 * pow2(pt2) * s - 2.0 *
pow2(pz) * s - 4.0 * pE * std::pow(s, 1.5) + pow2(s) - 2.0 * pow2(m2) *
(pow2(pE) + pow2(pt1) - pow2(pt2) - pow2(pz) - 2.0 * pE * sqrt_s + s) - 2.0 *
pow2(m1) * (pow2(m2) + pow2(pE) - pow2(pt1) + pow2(pt2) - pow2(pz) - 2.0 * pE *
sqrt_s + s)))) / (2.0 * (pow4(pE) + pow2(pow2(pz) - s) - 2.0 * pow2(pE) *
(pow2(pz) + s)));

          return p1z;
}
*/

// Case m3, m4 same or different
// Solve sqrt(s) = E3 + E4 + E5 and pz3 + pz4 + pz5 = 0 on shell
inline double SolvePz3_B(double m3, double m4, double pt3, double pt4,
                         double pz5, double E5, double s) {
  if (!(s > 0.0) || !std::isfinite(s)) {
    return -1.0;
  }
  const double t2 = E5 * E5;
  const double t3 = m3 * m3;
  const double t4 = m4 * m4;
  const double t5 = pt3 * pt3;
  const double t6 = pt4 * pt4;
  const double t7 = pz5 * pz5;
  const double t8 = std::sqrt(s);
  const double t9 = s * t2 * 6.0;
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
  const double t29 =
      t9 + t11 + t12 + t13 + t14 + t15 + t16 + t17 + t18 + t19 + t20 + t21 +
      t22 + t23 + t24 + t25 + t26 + t27 + t28 - E5 * t10 * 4.0 - s * t3 * 2.0 -
      s * t4 * 2.0 - s * t5 * 2.0 - s * t6 * 2.0 - s * t7 * 2.0 -
      t2 * t3 * 2.0 - t2 * t4 * 2.0 - t2 * t5 * 2.0 - t3 * t4 * 2.0 -
      t2 * t6 * 2.0 - t2 * t7 * 2.0 - t3 * t6 * 2.0 - t4 * t5 * 2.0 -
      t5 * t6 * 2.0 - E5 * t2 * t8 * 4.0;
  if (t29 < -KINEMATICS_EPS) {
    return -1.0;
  }
  const double denom = s + t2 - t7 - E5 * t8 * 2.0;
  if (std::abs(denom) < KINEMATICS_EPS) {
    return -1.0;
  }
  const double t30 = gra::math::msqrt(t29);
  const double t0 =
      pz5 * (-1.0 / 2.0) - (E5 * t30 * (1.0 / 2.0) - t8 * t30 * (1.0 / 2.0) +
                            pz5 * (t3 * (1.0 / 2.0) - t4 * (1.0 / 2.0) +
                                   t5 * (1.0 / 2.0) - t6 * (1.0 / 2.0))) /
                               denom;

  return std::isfinite(t0) ? t0 : -1.0;
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

% Elastic case
sol_A = simplify(subs(S, [m4], [m3]), 'Steps', 25);
ccode(sol_A, 'File', 'temp.c'); readwritetext('temp.c','solution_A.c', str);

% Generic inelastic case
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
inline bool SolvePzInputHasPhaseSpace(double m3, double m4, double pt3,
                                      double pt4, double pz5, double E5,
                                      double s) {
  if (!(s > 0.0) || !std::isfinite(s) || !std::isfinite(m3) ||
      !std::isfinite(m4) || !std::isfinite(pt3) || !std::isfinite(pt4) ||
      !std::isfinite(pz5) || !std::isfinite(E5)) {
    return false;
  }
  if (m3 < 0.0 || m4 < 0.0 || pt3 < 0.0 || pt4 < 0.0 || E5 < 0.0) {
    return false;
  }

  const double sqrt_s = std::sqrt(s);
  const double remaining_E = sqrt_s - E5;
  if (!(remaining_E > 0.0)) {
    return false;
  }

  const double mt3 = std::sqrt(m3 * m3 + pt3 * pt3);
  const double mt4 = std::sqrt(m4 * m4 + pt4 * pt4);
  const double min_remaining_E =
      std::sqrt((mt3 + mt4) * (mt3 + mt4) + pz5 * pz5);
  const double scale =
      std::max(1.0, std::max(std::abs(sqrt_s), min_remaining_E));

  return remaining_E + 100.0 * KINEMATICS_EPS * scale >= min_remaining_E;
}

// Check sqrt(s) = E3 + E4 + E5 after imposing pz4 = -pz3-pz5
inline bool SolvePzSolutionConservesEnergy(double m3, double m4, double pt3,
                                           double pt4, double pz5, double E5,
                                           double s, double pz3) {
  if (!std::isfinite(pz3)) {
    return false;
  }

  const double pz4 = -pz3 - pz5;
  const double E3 = std::sqrt(m3 * m3 + pt3 * pt3 + pz3 * pz3);
  const double E4 = std::sqrt(m4 * m4 + pt4 * pt4 + pz4 * pz4);
  const double sqrt_s = std::sqrt(s);
  const double scale = std::max(1.0, std::abs(sqrt_s));

  return std::abs(E3 + E4 + E5 - sqrt_s) <= 1e-8 * scale;
}

// Compute the physical longitudinal solution or -1 for invalid kinematics
inline double SolvePz(double m3, double m4, double pt3, double pt4, double pz5,
                      double E5, double s) {
  constexpr double EPS = 1e-9;

  if (!SolvePzInputHasPhaseSpace(m3, m4, pt3, pt4, pz5, E5, s)) {
    return -1.0;
  }

  const double pz3 = (std::abs(m3 - m4) < EPS)
                         ? SolvePz3_A(m3, pt3, pt4, pz5, E5, s)
                         : SolvePz3_B(m3, m4, pt3, pt4, pz5, E5, s);

  return SolvePzSolutionConservesEnergy(m3, m4, pt3, pt4, pz5, E5, s, pz3)
             ? pz3
             : -1.0;
}

// Build an on-shell forward system of fixed mass from its large light-cone fraction
// p_large = fraction p_beam, p_small = (M_f^2+pT^2)/p_large
inline bool BuildForwardSystem(const M4Vec &beam, double fraction, double px,
                               double py, double final_mass2, bool plus_side,
                               M4Vec &forward) {
  if (!(fraction > 0.0) || fraction > 1.0) {
    return false;
  }

  if (!std::isfinite(final_mass2) || final_mass2 < 0.0) {
    return false;
  }
  const double pt2 = px * px + py * py;
  if (plus_side) {
    const double plus = fraction * beam.LightconePos();
    if (!(plus > 0.0)) {
      return false;
    }
    const double minus = (final_mass2 + pt2) / plus;
    forward = M4Vec(px, py, 0.5 * (plus - minus), 0.5 * (plus + minus));
  } else {
    const double minus = fraction * beam.LightconeNeg();
    if (!(minus > 0.0)) {
      return false;
    }
    const double plus = (final_mass2 + pt2) / minus;
    forward = M4Vec(px, py, 0.5 * (plus - minus), 0.5 * (plus + minus));
  }
  // The light-cone construction is algebraically on shell; avoid a large E2-pz2
  // cancellation here
  return std::isfinite(forward.E()) && std::isfinite(forward.Pz()) &&
         forward.E() > 0.0;
}

// Build an exactly on-shell elastic forward particle from its large light-cone fraction
inline bool BuildForwardParticle(const M4Vec &beam, double fraction, double px,
                                 double py, bool plus_side, M4Vec &forward) {
  return BuildForwardSystem(beam, fraction, px, py,
                            std::max(beam.M2(), 0.0), plus_side, forward);
}

// Build one EPA transfer at fixed xi, forward mass and transverse recoil
// t = -[qT^2+xi(M_f^2-M_i^2)+xi^2 M_i^2]/(1-xi)
inline bool BuildEPATransfer(const M4Vec &beam, double xi, double qx,
                            double qy, double initial_mass2,
                            double final_mass2, bool plus_side, M4Vec &forward,
                            M4Vec &transfer, double &t) {
  if (!(xi > 0.0) || xi >= 1.0 || !std::isfinite(qx) || !std::isfinite(qy) || !std::isfinite(initial_mass2) ||
      initial_mass2 < 0.0 || !std::isfinite(beam.E()) || !(beam.E() > 0.0) || std::fpclassify(beam.Pt2()) != FP_ZERO) {
    return false;
  }
  const auto        b               = beam.Contravariant<long double>();
  const long double beam_mass2      = beam.Invariant<long double>();
  const long double shell_scale     = std::max(static_cast<long double>(std::numeric_limits<double>::min()),
                                               SquaredNorm(b) + std::abs(static_cast<long double>(initial_mass2)));
  const long double shell_tolerance = 64.0L * std::numeric_limits<double>::epsilon() * shell_scale;
  if (!std::isfinite(beam_mass2) || std::abs(beam_mass2 - initial_mass2) > shell_tolerance) { return false; }
  const double retained = 1.0 - xi;
  if (!BuildForwardSystem(beam, retained, -qx, -qy, final_mass2,
                          plus_side, forward)) {
    return false;
  }

  const double qt2 = qx * qx + qy * qy;
  t = -(qt2 + xi * (final_mass2 - initial_mass2) +
        xi * xi * initial_mass2) /
      retained;
  if (!std::isfinite(t) || t > 0.0) {
    return false;
  }

  if (plus_side) {
    const double beam_plus = beam.LightconePos();
    if (!(beam_plus > 0.0)) { return false; }
    const double plus          = xi * beam_plus;
    const double forward_minus = (final_mass2 + qt2) / (retained * beam_plus);
    const double minus         = initial_mass2 / beam_plus - forward_minus;
    transfer = M4Vec(qx, qy, 0.5 * (plus - minus),
                     0.5 * (plus + minus));
  } else {
    const double beam_minus = beam.LightconeNeg();
    if (!(beam_minus > 0.0)) { return false; }
    const double minus        = xi * beam_minus;
    const double forward_plus = (final_mass2 + qt2) / (retained * beam_minus);
    const double plus         = initial_mass2 / beam_minus - forward_plus;
    transfer = M4Vec(qx, qy, 0.5 * (plus - minus),
                     0.5 * (plus + minus));
  }
  return std::isfinite(transfer.E()) && std::isfinite(transfer.Pz());
}

// Build an exactly on-shell forward particle from xi, t and azimuth
// pT^2 = -(1-xi)t - xi^2 m^2
inline bool BuildForwardParticleXiT(const M4Vec &beam, double xi, double t,
                                    double phi, bool plus_side,
                                    M4Vec &forward) {
  if (!(xi > 0.0) || xi >= 1.0 || t > 0.0) {
    return false;
  }

  const double fraction = 1.0 - xi;
  const double mass2 = std::max(beam.M2(), 0.0);
  double pt2 = -fraction * t - xi * xi * mass2;
  constexpr double tolerance = 1.0e-10;
  if (pt2 < -tolerance) {
    return false;
  }
  pt2 = std::max(pt2, 0.0);
  const double pt = math::msqrt(pt2);
  return BuildForwardParticle(beam, fraction, pt * std::cos(phi),
                              pt * std::sin(phi), plus_side, forward);
}

// Compute the exact large light-cone momentum loss of one forward beam particle
// xi = 1 - p_f^large/p_beam^large
// At high energy it is 1-pz_f/pz_b + O(mT^2/E^2)
inline double LongitudinalMomentumLoss(const M4Vec &beam, const M4Vec &forward,
                                       bool plus_side) {
  const double beam_component =
      plus_side ? beam.LightconePos() : beam.LightconeNeg();
  const double final_component =
      plus_side ? forward.LightconePos() : forward.LightconeNeg();
  if (!(beam_component > 0.0) || !std::isfinite(beam_component) ||
      !std::isfinite(final_component)) {
    throw std::invalid_argument("gra::kinematics::LongitudinalMomentumLoss: "
                                "invalid light-cone component");
  }
  return 1.0 - final_component / beam_component;
}

// Build continuum transverse momenta from adjacent differences and total recoil
// p_i - p_{i+1} = q_i and sum_i p_i = -p1f-p2f
inline bool BuildCentralTransverseMomenta(const M4Vec &p1f, const M4Vec &p2f,
                                          const std::vector<M4Vec> &q,
                                          std::vector<M4Vec> &p) {
  if (q.empty()) {
    p.clear();
    return false;
  }

  const std::size_t multiplicity = q.size() + 1;
  M4Vec weighted_difference_sum;
  for (const auto &i : aux::indices(q)) {
    weighted_difference_sum += q[i] * static_cast<double>(multiplicity - 1 - i);
  }

  p.assign(multiplicity, M4Vec());
  p[0] =
      (weighted_difference_sum - p1f - p2f) / static_cast<double>(multiplicity);
  for (std::size_t i = 1; i < multiplicity; ++i) {
    p[i] = p[i - 1] - q[i - 1];
  }
  return true;
}

// Center of mass momentum^2
// for e.g. 1+2 (p_{cm}^2) -> 3+4 (p_{cm}'^2)
//
// dsigma/dOmega (s,theta) = 1/(64*pi^2*s) p_{cm}'/p_{cm} |M|^2,
// where dOmega = dcostheta dphi, if no initial state polarization, no phi dep
//
// Then integrating phi: dcostheta dphi / (16*pi^2) = dcostheta/(8*pi), gives
// dsigma/dt (s,t) = 1/(64*pi*s*p_{cm}^2) |M(s,t)|^2
//
// p_cm^2 = [s-(m1+m2)^2][s-(m1-m2)^2]/(4s)
inline double pcm2(double s, double m1, double m2) {
  if (!std::isfinite(s) || !std::isfinite(m1) || !std::isfinite(m2) ||
      s <= 0.0 || m1 < 0.0 || m2 < 0.0 || std::sqrt(s) < m1 + m2) {
    return 0.0;
  }
  return 1 / (4.0 * s) * (s - gra::math::pow2(m1 + m2)) *
         (s - gra::math::pow2(m1 - m2));
}

// Compute the exact invariant incoming-state Moller flux
// F = 2 sqrt[lambda(s,m1^2,m2^2)]
inline double MollerFlux(const M4Vec &beam1, const M4Vec &beam2) {
  // (p1.p2)^2 - m1^2 m2^2 = |E1 p2 - E2 p1|^2 - |p1 cross p2|^2
  const std::array<long double, 3> p1 = {beam1.Px(), beam1.Py(), beam1.Pz()};
  const std::array<long double, 3> p2 = {beam2.Px(), beam2.Py(), beam2.Pz()};
  auto relative = p2;
  gra::Scale(relative, static_cast<long double>(beam1.E()));
  gra::AddScaled(relative, p1, -static_cast<long double>(beam2.E()));
  const long double flux2 = gra::SquaredNorm(relative) - gra::SquaredNorm(gra::CrossProduct(p1, p2));
  return static_cast<double>(4.0L * std::sqrt(std::max(0.0L, flux2)));
}

// Compute the collision CM rapidity in the laboratory frame
// y_CM = (1/2)log[(E1+E2+pz1+pz2)/(E1+E2-pz1-pz2)]
inline double CollisionRapidity(const M4Vec &beam1, const M4Vec &beam2) {
  const M4Vec total = beam1 + beam2;
  const double rapidity = total.Rap();
  if (!std::isfinite(rapidity) || !(total.M2() > 0.0)) {
    throw std::invalid_argument("CollisionRapidity: invalid incoming beams");
  }
  return rapidity;
}

// 2->2 scattering invariants (p1 + p2  -> p3 + p4)
//
// Compute t = (p1-p3)^2 including all external masses
inline double mandelstam_t(const M4Vec &p1, const M4Vec &p3) {
  return (p1 - p3).M2();
}

// Compute u = (p1-p4)^2 including all external masses
inline double mandelstam_u(const M4Vec &p1, const M4Vec &p4) {
  return (p1 - p4).M2();
}

// 1-body Lorentz invariant phase space volume integral for energy < MAX (GeV)
// \int dPS_1 = \int d^3p / (2E(2pi)^3)
//            = \int_0^{4\pi} d^2\Omega/(2pi)^3 \int_0^{\sqrt{MAX^2-m^2}} p^2
//            dp/(2E)
//            = 1/(8\pi^2) [p_max MAX - m^2 ln((p_max + MAX)/m)]
//
// with E^2 = p^2 + m^2
//
inline double dPhi1(double MAX, double m) {
  using gra::math::pow2;

  const double mass = std::abs(m);
  if (!(MAX > mass)) {
    return 0.0;
  }

  const double pmax = std::sqrt((MAX - mass) * (MAX + mass));
  if (math::IsZero(mass)) {
    return pow2(MAX) / (8.0 * pow2(gra::math::PI));
  }

  const double r = pmax / mass;
  if (r < 0.01) {
    // Integral of p^2/E through order (p/m)^9, avoiding threshold cancellation
    const double r2 = r * r;
    return pmax * pmax * r * (2.0 / 3.0 + r2 * (-1.0 / 5.0 + r2 * (3.0 / 28.0 - r2 * 5.0 / 72.0))) /
           (8.0 * pow2(gra::math::PI));
  }
  return (pmax * MAX - pow2(mass) * std::asinh(r)) /
         (8.0 * pow2(gra::math::PI));
}

// Two-body Lorentz invariant phase space volume, evaluated in the rest frame
// p_norm is the result given by DecayMomentum()
// Phi_2 = p_norm/(4pi E)
inline double dPhi2(double E, double p_norm) {
  if (E <= 0.0 || p_norm <= 0.0) {
    return 0.0;
  }
  // Sampling box volume cos(theta) ~ [-1, 1], phi ~ [0,2pi]
  const double V = 2.0 * 2.0 * gra::math::PI;
  return V * p_norm / (16.0 * gra::math::pow2(gra::math::PI) * E);
}

// Kallen triangle function (square root applied here)
// symmetric under interchange of x <-> y <-> z:
// lambda^(1/2) = ( (x - y - z)^2 - 4yz )^(1/2)
//
// 1 / lambda^(1/2)(s,m_a^2,m_b^2) -> 1/s when s >> m_a^2, m_b^2
//
// Note x,y,z are Lorentz scalars (ENERGY squared variables)
inline double SqrtKallenLambda(double x, double y, double z) {
  // Subtract the largest terms first to retain small phase space near threshold
  std::array<long double, 3> q = {x, y, z};
  std::sort(q.begin(), q.end());
  const long double difference = (q[2] - q[1]) - q[0];
  const long double lambda = std::fma(difference, difference, -4.0L * q[0] * q[1]);
  return static_cast<double>(std::sqrt(std::max(0.0L, lambda)));
}

// beta_12 = sqrt[lambda(s,m1^2,m2^2)]/s
inline double beta12(double s, double m1, double m2) {
  if (!std::isfinite(s) || !std::isfinite(m1) || !std::isfinite(m2) ||
      s <= 0.0 || m1 < 0.0 || m2 < 0.0 || std::sqrt(s) < m1 + m2) {
    return 0.0;
  }
  return SqrtKallenLambda(s, gra::math::pow2(m1), gra::math::pow2(m2)) / s;
}

// Kinematic (momentum) bound for phase space (decays)
// This is same as: SqrtKallenLambda(x^2, y^2, z^2) / (2x)
//
// Note x,y,z [mother, daughter1, daughter2] in units of ENERGY (not ENERGY^2)
inline double DecayMomentum(double x, double y, double z) {
  if (!std::isfinite(x) || !std::isfinite(y) || !std::isfinite(z) || x <= 0.0 ||
      y < 0.0 || z < 0.0 || x < y + z) {
    return 0.0;
  }

  const double radicand = (x - y - z) * (x + y + z) * (x + y - z) * (x - y + z);
  if (radicand <= 0.0) {
    return 0.0;
  }

  return std::sqrt(radicand) / (2.0 * x);
}

// Massive two body phase space integral volume PS^{2}(M^2; m_A^2, m_B^2)
// Phi_2 = sqrt[lambda(M^2,mA^2,mB^2)]/(8pi M^2)
inline double PS2Massive(double M2, double mA2, double mB2) {
  if (!std::isfinite(M2) || !std::isfinite(mA2) || !std::isfinite(mB2) ||
      M2 <= 0.0 || mA2 < 0.0 || mB2 < 0.0 ||
      gra::math::msqrt(M2) < gra::math::msqrt(mA2) + gra::math::msqrt(mB2)) {
    return 0.0;
  }
  return SqrtKallenLambda(M2, mA2, mB2) / (8.0 * gra::math::PI * M2);
}

// Massless phase space integral volume PS^{(n)} = \int dPS^{(n)}
// Phi_n(M^2) = (M^2)^(n-2)/[2(4pi)^(2n-3) Gamma(n) Gamma(n-1)]
// where M2 (GeV^2) is the total CMS energy^2 of the n-body system
// Use this e.g. as a (debug) reference of other routines
// Gives 1/(8*PI) for n = 2 case, regardless of M^2
//
// tgamma is C++/11 math Gamma function
inline double PSnMassless(double M2, int n) {
  if (M2 <= 0.0 || n < 2) {
    return 0.0;
  }
  // const double denom = tgamma(n) * tgamma(n - 1.0);
  const double denom =
      gra::math::factorial(n - 1) * gra::math::factorial(n - 2);
  const double E = gra::math::msqrt(M2);
  const double PS = (1.0 / (2.0 * std::pow(4.0 * gra::math::PI, 2 * n - 3))) *
                    (std::pow(E, 2 * n - 4) / denom);

  return PS;
}

// Exact 2-body partial decay width,
// is not dependent on momentum => phase space and matrix element factorize
// Gamma = lambda^(1/2)(M0^2,m1^2,m2^2) |M|^2 /(16 pi S M0^3)
//
// Input: Mother mass M0^2
//        Decay product 1 mass^2
//        Decay product 2 mass^2
//        Decay matrix element squared |M|^2 value
//        Final state symmetry factor S
//
// [REFERENCE: Alwall et al, https://arxiv.org/abs/1402.1178]
inline double PDW2body(double M0_2, double m1_2, double m2_2, double MatElem2,
                       double S_factor) {
  if (!std::isfinite(M0_2) || !std::isfinite(m1_2) || !std::isfinite(m2_2) ||
      !std::isfinite(MatElem2) || !std::isfinite(S_factor) || M0_2 <= 0.0 ||
      m1_2 < 0.0 || m2_2 < 0.0 || MatElem2 < 0.0 || S_factor <= 0.0 ||
      gra::math::msqrt(M0_2) <
          gra::math::msqrt(m1_2) + gra::math::msqrt(m2_2)) {
    return 0.0;
  }
  return SqrtKallenLambda(M0_2, m1_2, m2_2) * MatElem2 /
         (16.0 * gra::math::PI * S_factor * math::pow3(math::msqrt(M0_2)));
}

// Relativistic velocity of a particle with m^2
// in the system with invariant shat
// beta = sqrt(1 - 4m^2/shat)
inline double Beta(double m2, double shat) {
  if (!std::isfinite(m2) || !std::isfinite(shat) || m2 < 0.0 || shat <= 0.0 ||
      shat < 4.0 * m2) {
    return 0.0;
  }
  return gra::math::msqrt(1.0 - 4.0 * m2 / shat);
}

// 2->2 Scattering E1 + E2 -> E3 + E4 angle as a function of (s,t,m_i^2)
// cos(theta*) = [s^2+s(2t-sum_i m_i^2)+(m1^2-m2^2)(m3^2-m4^2)]/[sqrt(lambda12 lambda34)]
// [REFERENCE: Kallen, Elementary Particle Physics, Addison-Wesley, 1964]
inline double CosthetaStar(double s, double t, double E1, double E2, double E3,
                           double E4) {
  const double numer =
      s * s + s * (2.0 * t - E1 - E2 - E3 - E4) + (E1 - E2) * (E3 - E4);
  const double denom =
      SqrtKallenLambda(s, E1, E2) * SqrtKallenLambda(s, E3, E4);

  if (denom <= 0.0) {
    return 0.0;
  }
  return numer / denom;
}

// Two-body CM momentum norms and the forward Mandelstam transfer
struct ScatteringCM {
  long double pin = 0.0L;
  long double pout = 0.0L;
  long double tf = 0.0L;
};

// Compute physical two-body CM kinematics with a stable forward transfer
// [REFERENCE: PDG, Kinematics, https://pdg.lbl.gov/2025/reviews/rpp2025-rev-kinematics.pdf]
inline bool TwoBodyCM(double s, double a, double b, double c, double d, ScatteringCM &cm) {
  cm = {};
  if (!std::isfinite(s) || !(s > 0.0)) { return false; }
  for (const double mass2 : {a, b, c, d}) {
    if (!std::isfinite(mass2) || mass2 < 0.0) { return false; }
  }
  const long double root = std::sqrt(static_cast<long double>(s));
  const auto momentum = [root](double x, double y) {
    const long double mx = std::sqrt(static_cast<long double>(x));
    const long double my = std::sqrt(static_cast<long double>(y));
    if (!(root > mx + my)) { return 0.0L; }
    return std::sqrt((root - mx - my) * (root + mx + my) *
                     (root + mx - my) * (root - mx + my)) / (2 * root);
  };
  cm.pin = momentum(a, b);
  cm.pout = momentum(c, d);
  if (!(cm.pin > 0.0L) || !(cm.pout > 0.0L)) { return false; }
  const long double e1 = (s + static_cast<long double>(a) - b) / (2 * root);
  const long double e3 = (s + static_cast<long double>(c) - d) / (2 * root);
  const long double de = (static_cast<long double>(a) - b - c + d) / (2 * root);
  const long double dp = (de * (e1 + e3) - a + c) / (cm.pin + cm.pout);
  cm.tf = (de - dp) * (de + dp);
  return std::isfinite(cm.tf);
}

// Compute forward two-body scattering momenta without subtracting cos(theta) from one
inline bool ForwardScattering(double s, double a, double b, double c, double d,
                              double t, double &pz, double &pt) {
  pz = pt = 0.0;
  ScatteringCM cm;
  if (!std::isfinite(t) || !TwoBodyCM(s, a, b, c, d, cm)) { return false; }
  const long double loss = (cm.tf - t) / (2 * cm.pin);
  if (!std::isfinite(loss) || loss < 0.0L || loss > cm.pout) { return false; }
  pz = static_cast<double>(cm.pout - loss);
  pt = static_cast<double>(std::sqrt(loss * (2 * cm.pout - loss)));
  return std::isfinite(pz) && std::isfinite(pt);
}

// Compute the physical t interval or an empty interval below threshold
inline void Two2TwoLimit(double s, double a, double b, double c, double d,
                         double &tmin, double &tmax) {
  ScatteringCM cm;
  if (!TwoBodyCM(s, a, b, c, d, cm)) {
    tmin = 1.0;
    tmax = 0.0;
    return;
  }
  tmin = static_cast<double>(cm.tf - 4 * cm.pin * cm.pout);
  tmax = static_cast<double>(cm.tf);
}

// Standard conventions:
//
// E^2 = m^2 + |p|^2, E = gamma m, p = gamma m v in c = 1 units
// gamma = E/m = 1/sqrt(1-beta^2), beta = |v| = |p|/E

//
// d^3p / [2E(2pi)^3] = pT dpT dphi dy / [2(2pi)^3]
//

// Generic Lorentz boost into the direction given by 4-momentum of "boost"
// E' = gamma(E + beta.p), p' = p + [gamma E + gamma^2(beta.p)/(1+gamma)] beta
//
// A. sign = -1 for boost to the rest frame (out of lab, for example)
// B. sign =  1 out from the rest frame (back to lab, for example)
//
template <typename T>
inline void LorentzBoost(const T &boost, double M0, T &p, int sign) {
  if ((sign != -1 && sign != 1) || !std::isfinite(M0) || M0 <= 0.0 ||
      !std::isfinite(boost.Px()) || !std::isfinite(boost.Py()) ||
      !std::isfinite(boost.Pz()) || !std::isfinite(boost.E()) ||
      boost.E() <= 0.0 || !std::isfinite(p.Px()) || !std::isfinite(p.Py()) ||
      !std::isfinite(p.Pz()) || !std::isfinite(p.E())) {
    p = T(0.0, 0.0, 0.0, -1.0);
    return;
  }

  // Use gamma^2 = 1 + p^2/M^2 and gamma beta = p/M on the validated mass shell
  const long double mass         = M0;
  const long double bx           = boost.Px();
  const long double by           = boost.Py();
  const long double bz           = boost.Pz();
  const long double spatial2     = bx * bx + by * by + bz * bz;
  const long double invariant    = boost.template Invariant<long double>();
  const long double target_mass2 = mass * mass;
  const long double shell_scale =
      std::max(static_cast<long double>(std::numeric_limits<double>::min()),
               SquaredNorm(boost.template Contravariant<long double>()) + std::abs(target_mass2));
  const long double shell_tolerance = 128.0L * std::numeric_limits<double>::epsilon() * shell_scale;
  const long double gamma           = std::sqrt(1.0L + spatial2 / target_mass2);
  if (!std::isfinite(spatial2) || !std::isfinite(invariant) || !std::isfinite(gamma) ||
      !(invariant > 0.0L) || std::abs(invariant - target_mass2) > shell_tolerance) {
    p = T(0.0, 0.0, 0.0, -1.0);
    return;
  }

  const long double px                 = p.Px();
  const long double py                 = p.Py();
  const long double pz                 = p.Pz();
  const long double p_energy           = p.E();
  const long double dot                = bx * px + by * py + bz * pz;
  const long double parallel           = dot / (target_mass2 * (gamma + 1.0L));
  const long double signed_energy      = sign * p_energy / mass;
  const long double scale              = parallel + signed_energy;
  const long double transformed_energy = gamma * p_energy + sign * dot / mass;
  const long double transformed_px     = px + scale * bx;
  const long double transformed_py     = py + scale * by;
  const long double transformed_pz     = pz + scale * bz;
  const long double limit = std::numeric_limits<double>::max();
  if (!(std::abs(transformed_px) <= limit && std::abs(transformed_py) <= limit &&
        std::abs(transformed_pz) <= limit && std::abs(transformed_energy) <= limit)) {
    p = T(0.0, 0.0, 0.0, -1.0);
    return;
  }
  p = T(static_cast<double>(transformed_px), static_cast<double>(transformed_py), static_cast<double>(transformed_pz),
        static_cast<double>(transformed_energy));
}

// Validate the mother four-momentum against the supplied decay mass
template <typename T> inline bool ValidDecayMother(const T &mother, double mass) {
  T rest(0.0, 0.0, 0.0, mass);
  LorentzBoost(mother, mass, rest, 1);
  return rest.E() > 0.0;
}

// Build a soft remnant pair along its incoming axes and reject negative-energy production
inline bool BuildSoftRemnants(const std::array<M4Vec, 2> &incoming,
                              const std::array<double, 2> &x,
                              const std::array<double, 2> &mass2,
                              const std::array<M4Vec, 2> &transverse,
                              std::array<M4Vec, 2> &outgoing, M4Vec &system) {
  const M4Vec total = incoming[0] + incoming[1];
  const double mass = total.M();
  if (!ValidDecayMother(total, mass) || incoming[0].M2() < 0.0 || incoming[1].M2() < 0.0) { return false; }
  const double p = DecayMomentum(mass, incoming[0].M(), incoming[1].M());
  if (!(p > 0.0)) { return false; }
  M4Vec axis = incoming[0];
  LorentzBoost(total, mass, axis, -1);
  const auto rotation = M4Vec(0.0, 0.0, 1.0, 0.0).RotationTo(axis.P3());
  std::array<M4Vec, 2> remnant;
  for (const auto &i : aux::indices(remnant)) {
    if (!(x[i] > 0.0 && x[i] < 1.0) || !std::isfinite(mass2[i]) || mass2[i] < 0.0) { return false; }
    const double pz = (i == 0 ? 1.0 : -1.0) * (1.0 - x[i]) * p;
    remnant[i] = M4Vec(transverse[i].Px(), transverse[i].Py(), pz,
                       std::sqrt(mass2[i] + transverse[i].Pt2() + pz * pz));
    remnant[i].Rotate(rotation);
    LorentzBoost(total, mass, remnant[i], 1);
    if (!(remnant[i].E() > 0.0)) { return false; }
  }
  const M4Vec produced = total - remnant[0] - remnant[1];
  if (!(produced.E() > 0.0) || !std::isfinite(produced.M2()) || !(produced.M2() > 0.0)) { return false; }
  outgoing = remnant;
  system = produced;
  return true;
}

/*
template <typename T>
inline void LorentzBoostSpaceLike(const T &boost, T &p, int sign) {
  std::valarray<double> p3 = {p.Px(), p.Py(), p.Pz()};

  // Mother beta-vector
  const std::valarray<double> beta = {sign * boost.Px() / boost.E(), sign *
boost.Py() / boost.E(), sign * boost.Pz() / boost.E()};
  // Mother gamma-factor
  const double gamma = boost.E() / M0;

  // Apply transforms
  const double kappa1 = (beta * p3).sum();  // Inner product
  const double kappa2 = gamma * (gamma * kappa1 / (1.0 + gamma) + p.E());

  // Vector sum
  p3 += kappa2 * beta;

  p = T(p3[0], p3[1], p3[2], gamma * (p.E() + kappa1));
}
*/

// Flat variables in spherical coordinates
// cos(theta) ~ U[-1,1], phi ~ U[0,2pi]
template <typename T2>
inline void FlatIsotropic(double &costheta, double &sintheta, double &phi,
                          T2 &rng) {
  costheta = rng.U(-1.0, 1.0); // cos(theta) flat [-1,1]
  sintheta = gra::math::msqrt(1.0 - gra::math::pow2(costheta));
  phi = rng.U(0.0, 2.0 * gra::math::PI); // phi flat [0,2pi]
}

// Isotropic decay 0 -> 1 + 2 in spherical coordinates with flat cos(theta),phi
// p1 = (sqrt(m1^2+p^2),p n_hat), p2 = (sqrt(m2^2+p^2),-p n_hat)
template <typename T1, typename T2>
inline void Isotropic(double pnorm, T1 &p1, T1 &p2, double m1, double m2,
                      T2 &rng) {
  using gra::math::msqrt;
  using gra::math::pow2;

  double costheta, sintheta, phi;
  FlatIsotropic(costheta, sintheta, phi, rng);

  // Jacobian of spherical coordinates
  const std::valarray<double> k = {pnorm * sintheta * std::cos(phi),
                                   pnorm * sintheta * std::sin(phi),
                                   pnorm * costheta};

  // Energies by on-shell condition
  const std::valarray<double> e = {msqrt(pow2(m1) + pow2(pnorm)),
                                   msqrt(pow2(m2) + pow2(pnorm))};

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
inline MCW TwoBodyPhaseSpace(const T1 &mother, double M0,
                             const std::vector<double> &m, std::vector<T1> &p,
                             T2 &rng) {
  if (m.size() != 2 || !HasPhaseSpace(M0, m) || !ValidDecayMother(mother, M0)) {
    return InvalidPhaseSpacePoint();
  }

  // Two particles
  p.resize(2);

  // First isotropic decay
  const double pnorm = DecayMomentum(M0, m[0], m[1]);
  Isotropic(pnorm, p[0], p[1], m[0], m[1], rng);

  // Boost daughters to the lab frame
  const int sign = 1;
  for (const auto &i : {0, 1}) {
    LorentzBoost(mother, M0, p[i], sign);
    if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
  }

  // return phase space weight
  const double weight = dPhi2(M0, pnorm);
  return MCW(weight);
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
// Compute: phase space weight for a valid fragmentation and -1.0 for a
// kinematically impossible
// With unweight=true MCW stores every rejection trial, so its normalization is
// statistically valid only after MCW objects are accumulated across generated
// events
//
template <typename T1, typename T2>
inline MCW ThreeBodyPhaseSpace(const T1 &mother, double M0,
                               const std::vector<double> &m, std::vector<T1> &p,
                               bool unweight, T2 &rng) {
  if (m.size() != 3 || !HasPhaseSpace(M0, m) || !ValidDecayMother(mother, M0)) {
    return InvalidPhaseSpacePoint();
  }

  p.resize(3); // Three final states
  T1 p12;

  // Phase space boundaries [min,max]
  const std::array<double, 2> m12bound = {m[1] + m[2], M0 - m[0]};
  if (m12bound[1] - m12bound[0] <= 0.0) {
    return InvalidPhaseSpacePoint();
  }
  double w_max = DecayMomentum(M0, m[0], m12bound[0]) *
                 DecayMomentum(m12bound[1], m[1], m[2]);
  if (unweight == false) {
    w_max = 0;
  }

  // Accceptance-Rejection
  const unsigned int MAXTRIAL = 1e8;
  std::array<double, 2> pnorm = {0.0, 0.0};
  double m12 = 0;
  MCW x;
  do {
    m12 = rng.U(m12bound[0], m12bound[1]); // Flat mass (in GeV, not GeV^2)
    pnorm[0] = DecayMomentum(M0, m[0], m12);
    pnorm[1] = DecayMomentum(m12, m[1], m[2]);

    const double w =
        (m12 / gra::math::PI) * dPhi2(M0, pnorm[0]) * dPhi2(m12, pnorm[1]);
    x.Push(w);
    if (x.GetN() > MAXTRIAL) {
      return InvalidPhaseSpacePoint();
    }
  } while ((pnorm[0] * pnorm[1]) < rng.U(0.0, w_max));

  // pM -> p0 + p12 in the pM r.f. and p12 -> p1 + p2 in the p12 r.f
  Isotropic(pnorm[0], p[0], p12, m[0], m12, rng);
  Isotropic(pnorm[1], p[1], p[2], m[1], m[2], rng);

  // Boost p[1] and p[2] out from the p12 rest frame to the pM rest frame
  const int sign = 1; // Boost sign
  for (const auto &i : {1, 2}) {
    LorentzBoost(p12, m12, p[i], sign);
    if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
  }

  // Boost p[0], p[1] and p[2] to the lab frame
  for (const auto &i : {0, 1, 2}) {
    LorentzBoost(mother, M0, p[i], sign);
    if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
  }

  // Close the recursive decay against the mother four-vector after boosts
  p[0] = mother - p[1] - p[2];
  if (!(p[0].E() > 0.0)) { return InvalidPhaseSpacePoint(); }

  // Phasespace weight [note m versus m^2 jacobian, taken into account]
  // [m12 volume (GeV)] x [dm_12] x [dphi2] x [dphi2]
  const double volume = (m12bound[1] - m12bound[0]);
  x = x * volume; // operator overloaded

  // Compute phase space weight
  return x;
}

// Decay setup for NBodyPhaseSpace
// Uses the "sorting algorithm", see description in F. James, CERN/68
// bool unweighted for (un)weighted operation
//
// Compute: W for a valid and -1.0 for a kinematically impossible
//
template <typename T1, typename T2>
inline MCW NBodySetup(const T1 &mother, double M0, const std::vector<double> &m,
                      std::vector<double> &M_eff, std::vector<double> &pnorm,
                      bool unweight, T2 &rng) {
  // Decay multiplicity (1 -> N decay)
  const unsigned int N = m.size();
  if (N < 2 || !HasPhaseSpace(M0, m) || !ValidDecayMother(mother, M0)) {
    return InvalidPhaseSpacePoint();
  }
  M_eff.assign(N, 0.0);
  pnorm.assign(N, 0.0);

  std::vector<double> randvec(N - 2, 0.0);

  // Random variables functor
  auto fillrandom = [&](double &value) -> void { value = rng.U(0, 1); };

  // Effective (intermediate) masses
  M_eff[0] = m[0];   // First daughter
  M_eff[N - 1] = M0; // Mother

  // Cumulative sum of daughter masses
  std::vector<double> sumvec;
  gra::CumSum(m, sumvec);

  // Sampling mass interval
  const double DELTA = M_eff[N - 1] - sumvec[N - 1];
  if (DELTA <= 0.0) {
    return InvalidPhaseSpacePoint();
  }

  // Generate effective masses functor, where sumvec[i] defines the minimum
  auto calcmass = [&]() -> void {
    for (std::size_t i = 1; i < N - 1; ++i) {
      M_eff[i] = sumvec[i] + randvec[i - 1] * DELTA;
    }
  };

  // -------------------------------------------------------------------
  // Recursively find maximum weight for Acceptance-Rejection
  double w_max = 1.0;
  std::array<double, 2> m_minmax = {0.0, DELTA + m[0]};
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

  if (unweight == false) {
    w_max = 0;
  }

  const unsigned int MAXTRIAL = 1e8;
  MCW x;
  double w = 0.0;
  do {
    // 1. Ordered random numbers
    std::for_each(randvec.begin(), randvec.end(), fillrandom);
    std::sort(randvec.begin(), randvec.end()); // Ascending

    // 2. Calculate effective masses
    calcmass();

    // 3. Calculate momentum norm values
    w = calcmomentum();

    x.Push(w);
    if (x.GetN() > MAXTRIAL) {
      return InvalidPhaseSpacePoint();
    }
  } while (unweight && w < rng.U(0.0, w_max));

  // Normalize the phase space integral
  const double volume = 1.0 / (2.0 * std::pow(2.0 * gra::math::PI, 2 * N - 3)) *
                        std::pow(DELTA, N - 2) / gra::math::factorial(N - 2) /
                        M0;
  x = x * volume; // operator overloaded

  // Compute phase space weight
  return x;
}

// Flat N-body phase space by implementing a number of (N-2) recursive 2-body
// decays
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
// For (faster) alternatives, investigate:
// M.M. Block, Monte Carlo phase space evaluation, Comp. Phys. Commun
// 69, 459 (1992)
// S. Platzer, RAMBO on diet, https://arxiv.org/abs/1308.2922
//
// Compute: weight for a valid fragmentation and -1.0 for a kinematically
// impossible
//
// With unweight=true MCW stores every rejection trial, so its normalization is
// statistically valid only after MCW objects are accumulated across generated
// events
//
template <typename T1, typename T2>
inline MCW NBodyPhaseSpace(const T1 &mother, double M0,
                           const std::vector<double> &m, std::vector<T1> &p,
                           bool unweight, T2 &rng) {
  using gra::math::msqrt;
  using gra::math::pow2;

  const int N = m.size();
  if (N < 2) {
    return InvalidPhaseSpacePoint();
  }

  // Generate effective masses and decay momentum
  std::vector<double> M_eff(N, 0.0); // Effective masses
  std::vector<double> pnorm(N, 0.0); // Decay momentum weights
  const MCW x = NBodySetup(mother, M0, m, M_eff, pnorm, unweight, rng);
  if (x.GetW() < 0) {
    return x;
  } // Impossible kinematics

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
    const int sign = 1; // positive
    LorentzBoost(k, k.M(), p[i], sign);
    if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }

    // Subtract the momentum
    k -= p[i];
  }
  // Finally, we are left with only one daughter
  p[0] = k;
  if (!(p[0].E() > 0.0)) { return InvalidPhaseSpacePoint(); }

  return x; // Phase-space weight
}

// Supply unit integration coordinates to the shared phase-space sampler
class PhaseSpaceCoordinates {
 public:
  // Bind a fixed integration point without introducing a random generator
  explicit PhaseSpaceCoordinates(const std::vector<double> &coordinates) : point(coordinates) {}

  // Map the next unit coordinate to the requested interval
  double U(double lower, double upper) { return lower + (upper - lower) * point.at(index++); }

 private:
  const std::vector<double> &point;
  std::size_t index = 0;
};

// Integrate all 3N-4 coordinates using the same N-body phase-space measure
template <typename T1>
inline MCW NBodyPhaseSpace(const T1 &mother, double M0, const std::vector<double> &m,
                           std::vector<T1> &p, const std::vector<double> &coordinates) {
  PhaseSpaceCoordinates point(coordinates);
  return NBodyPhaseSpace(mother, M0, m, p, false, point);
}

// Pure RAMBO (massless) flat phase space
//
// [REFERENCE: Kleiss R., Stirling W.J., Ellis S.D, 1986]
//
// See: https://www.ictp-saifr.org/wp-content/uploads/2019/02/newmc.pdf
//
// Compute: weight for a valid and -1.0 for a kinematically impossible
//
template <typename T1, typename T2>
inline MCW RamboMassless(const T1 &mother, double M0, std::vector<T1> &p,
                         T2 &rng, bool boost = true) {
  const int N = p.size();
  if (N < 2 || !std::isfinite(M0) || M0 <= 0.0 || (boost && !ValidDecayMother(mother, M0))) {
    return InvalidPhaseSpacePoint();
  }

  // Generate energies q_j^0 according to: q_j^0 exp(-q_j^0)
  // Generate momenta \vec{q}_j isotropically with |\vec{q}_j| = q_j^0
  //
  std::valarray<T1> q(N);
  for (const auto &i : aux::indices(p)) {
    const double q0 = -std::log(rng.U(0, 1) * rng.U(0, 1));

    double costheta = 0;
    double sintheta = 0;
    double phi = 0;
    FlatIsotropic(costheta, sintheta, phi, rng);
    const double pt = q0 * sintheta;

    q[i] = T1(pt * std::cos(phi), pt * std::sin(phi), q0 * costheta, q0);
  }

  // Q^\mu = \sum_j q_j^\mu
  const T1 Q = q.sum();
  const double M_Q = Q.M();
  if (M_Q <= 0.0) {
    return InvalidPhaseSpacePoint();
  }

  M3Vec b = Q.P3();
  for (const auto &i : aux::indices(b)) {
    b[i] = -b[i] / M_Q;
  }

  const double gamma = Q.Gamma();
  const double a = 1.0 / (1.0 + gamma);
  const double x = M0 / M_Q;

  // Boost to the rest
  // p_j <-- \Lambda q_j / x
  T1 psum;
  for (const auto &i : aux::indices(p)) {
    const double qE = q[i].E();
    const M3Vec q3 = q[i].P3();
    const double BQ = b[0] * q3[0] + b[1] * q3[1] + b[2] * q3[2];

    // Conformal transform
    M3Vec p3{};
    for (const auto &k : aux::indices(p3)) {
      p3[k] = x * (q3[k] + b[k] * (qE + a * BQ));
    }
    const double p0 = x * (gamma * qE + BQ);

    p[i] = T1(p3[0], p3[1], p3[2], p0);
  }

  // ---------------------------------------------------------------------
  // Boost all particles to the lab frame

  if (boost) {
    for (const auto &i : aux::indices(p)) {
      // p[i].Print("particle " + std::to_string(i));
      LorentzBoost(mother, M0, p[i], 1); // note plus
      if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
    }
  }
  // ---------------------------------------------------------------------

  // Event weight
  const double weight = PSnMassless(M0 * M0, N);

  return MCW(weight);
}

// Solve the massive RAMBO kinetic energy constraint and rescale the momenta
// sum_i (xi q_i)^2/(sqrt(m_i^2+(xi q_i)^2)+m_i) = M0-sum_i m_i
template <typename T>
inline bool RamboScale(double M0, const std::vector<double> &m, std::vector<T> &p) {
  using gra::aux::indices;
  if (m.size() != p.size() || !HasPhaseSpace(M0, m)) { return false; }
  std::vector<long double> q(m.size());
  std::vector<double> energy(m.size());
  long double mass_sum = 0.0L;
  long double momentum_sum = 0.0L;
  for (const auto &i : indices(m)) {
    q[i] = p[i].P3mod();
    if (!std::isfinite(q[i]) || !(q[i] > 0.0L)) { return false; }
    mass_sum += static_cast<long double>(m[i]);
    momentum_sum += q[i];
  }
  const long double kinetic = static_cast<long double>(M0) - mass_sum;
  if (!(kinetic > 0.0L) || !std::isfinite(momentum_sum)) { return false; }
  const long double tolerance = 64.0L * std::numeric_limits<double>::epsilon() * kinetic;
  long double lower = 0.0L;
  long double upper = static_cast<long double>(M0) / momentum_sum;
  long double xi = upper * std::sqrt((kinetic / M0) * (1.0L + mass_sum / M0));

  for (unsigned int iteration = 0; iteration < 128; ++iteration) {
    long double kinetic_sum = 0.0L;
    long double derivative = 0.0L;
    for (const auto &i : indices(m)) {
      const long double momentum = xi * q[i];
      const long double e = std::hypot(static_cast<long double>(m[i]), momentum);
      if (!(e > 0.0L) || !std::isfinite(e)) { return false; }
      energy[i] = static_cast<double>(e);
      kinetic_sum += momentum * (momentum / (e + m[i]));
      derivative += q[i] * (momentum / e);
    }
    const long double residual = kinetic_sum - kinetic;
    if (std::abs(residual) <= tolerance) {
      for (const auto &i : indices(p)) {
        p[i] *= static_cast<double>(xi);
        p[i].SetE(energy[i]);
      }
      return true;
    }
    if (!std::isfinite(residual) || !(derivative > 0.0L) || !std::isfinite(derivative)) { return false; }
    if (residual > 0.0L) {
      upper = xi;
    } else {
      lower = xi;
    }
    const long double next = xi - residual / derivative;
    xi = next > lower && next < upper ? next : (lower + upper) / 2.0L;
  }
  return false;
}

// Compute the massive RAMBO weight from daughter momenta in the mother rest frame
// [REFERENCE: Kleiss, Stirling and Ellis, Comput. Phys. Commun. 40 (1986) 359]
template <typename T>
inline double RamboWeight(double M0, const std::vector<T> &p) {
  if (!(M0 > 0.0) || !std::isfinite(M0) || p.size() < 2) { return 0.0; }
  double prod = 1.0;
  double sumA = 0.0;
  double sumB = 0.0;
  for (const auto &momentum : p) {
    const double k_mag = momentum.P3mod();
    const double k0 = momentum.E();
    if (!(k_mag > 0.0) || !(k0 > 0.0) || !std::isfinite(k_mag) || !std::isfinite(k0)) { return 0.0; }
    sumA += k_mag / M0;
    sumB += math::pow2(k_mag) / (k0 * M0);
    prod *= k0 / k_mag;
  }
  if (!(sumB > 0.0) || !(sumA > 0.0) || !std::isfinite(sumB) || !std::isfinite(sumA)) { return 0.0; }
  const double weight = PSnMassless(M0 * M0, p.size()) * std::pow(sumA, 2.0 * p.size() - 3.0) / (prod * sumB);
  return std::isfinite(weight) && weight > 0.0 ? weight : 0.0;
}

// RAMBO algorithm for flat phase space with masses
//
// [REFERENCE: Kleiss R., Stirling W.J., Ellis S.D, 1986]
//
// Compute: weight for a valid and -1.0 for a kinematically impossible
//
template <typename T1, typename T2>
inline MCW RamboMassive(const T1 &mother, double M0,
                        const std::vector<double> &m, std::vector<T1> &p,
                        T2 &rng) {
  const int N = m.size();
  if (N < 2 || !HasPhaseSpace(M0, m) || !ValidDecayMother(mother, M0)) {
    return InvalidPhaseSpacePoint();
  }

  // Make sure it is of right size
  p.resize(N);

  // Generate massless, last parameter set to false, we do boost to the lab here
  const MCW W0 = RamboMassless(mother, M0, p, rng, false);
  if (W0.GetW() < 0.0) {
    return W0;
  }

  double XMT = 0.0;
  for (const auto &i : aux::indices(m)) {
    XMT += m[i];
  }

  if (math::IsZero(XMT)) { // Purely massless case
    for (const auto &i : aux::indices(p)) {
      LorentzBoost(mother, M0, p[i], 1); // note plus
      if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
    }
    return W0;
  }

  // =====================================================================
  // Generate masses while resolving the available kinetic energy
  if (!RamboScale(M0, m, p)) { return InvalidPhaseSpacePoint(); }

  // =====================================================================
  // Compute event weight

  const double weight = RamboWeight(M0, p);
  if (!(weight > 0.0)) { return InvalidPhaseSpacePoint(); }

  // ---------------------------------------------------------------------
  // Boost all particles to the lab frame
  for (const auto &i : aux::indices(p)) {
    LorentzBoost(mother, M0, p[i], 1); // note plus
    if (!(p[i].E() > 0.0)) { return InvalidPhaseSpacePoint(); }
  }
  // ---------------------------------------------------------------------

  return MCW(weight);
}

// Find the closest 4-vector on lightcone
//
// Solution by Lagrange multipliers \Nabla f = \lambda \Nabla g
//
// eq(1) = 2*(px-px0) == -2*lambda*px
// eq(2) = 2*(py-py0) == -2*lambda*py
// eq(3) = 2*(pz-pz0) == -2*lambda*pz
// eq(4) = 2*(E-E0)   ==  2*lambda*E
// eq(5) = E^2        == (px^2 + py^2 + pz^2)
//
// S = solve(eq, [E,px,py,pz,lambda])
//
inline M4Vec LagrangeLightCone(const M4Vec &p0) {
  const double p3mod2 = p0.P3mod2();
  const double p3mod = math::msqrt(p3mod2);
  const double energy_sign = p0.E() < 0.0 ? -1.0 : 1.0;
  const double radius = 0.5 * (p3mod + std::abs(p0.E()));
  if (math::IsZero(p3mod2)) {
    return M4Vec(0.0, 0.0, radius, energy_sign * radius);
  }

  const double alpha = radius / p3mod;
  const double energy = energy_sign * radius;
  if (!std::isfinite(energy) || !std::isfinite(alpha)) {
    return M4Vec(0.0, 0.0, 0.0, 0.0);
  }

  return M4Vec(alpha * p0.Px(), alpha * p0.Py(), alpha * p0.Pz(), energy);
}

// Compute one four-vector Euclidean norm squared in a timelike-system rest frame
// ||v*||_E^2 = 2(v.P)^2/P^2 - v^2
template <typename Vector, typename System>
inline auto RestFrameNorm2(const Vector &vector, const System &system,
                           const typename System::value_type system_mass2) {
  const auto rest_energy =
      MinkowskiProduct(vector, system) / std::sqrt(system_mass2);
  return 2 * rest_energy * rest_energy - MinkowskiProduct(vector, vector);
}

// Construct a covariant on-shell incoming pair while preserving the final state
//
// The axis is the component of q1 orthogonal to the hard-system momentum P
// and k1,2 = P/2 +/- sqrt(P^2) n/2, with n^2 = -1 and n.P = 0
//
inline bool CovariantOnShellInitialState(const M4Vec &q1, const M4Vec &q2,
                                         const std::vector<M4Vec> &final,
                                         M4Vec &k1, M4Vec &k2) {
  if (final.empty()) {
    return false;
  }

  using Four = std::array<long double, 4>;
  Four hard{};
  for (const auto &particle : final) {
    if (!std::isfinite(particle.E()) || !std::isfinite(particle.Px()) ||
        !std::isfinite(particle.Py()) || !std::isfinite(particle.Pz())) {
      return false;
    }
    AddScaled(hard, particle.Contravariant<long double>(), 1.0L);
  }

  const long double hard_mass2 = MinkowskiProduct(hard, hard);
  if (!std::isfinite(hard_mass2) || hard_mass2 <= KINEMATICS_EPS) {
    return false;
  }

  const Four q1_extended = q1.Contravariant<long double>();
  const Four q2_extended = q2.Contravariant<long double>();
  if (!AllFinite(q1_extended) || !AllFinite(q2_extended)) {
    return false;
  }
  Four input_delta = Add(q1_extended, q2_extended);
  AddScaled(input_delta, hard, -1.0L);
  const long double hard_mass = std::sqrt(hard_mass2);
  const long double input_delta_m2 = MinkowskiProduct(input_delta, input_delta);
  long double input_closure2 = RestFrameNorm2(input_delta, hard, hard_mass2);
  const long double delta_rest_e =
      MinkowskiProduct(input_delta, hard) / hard_mass;
  if (!std::isfinite(input_delta_m2) || !std::isfinite(input_closure2) ||
      !std::isfinite(delta_rest_e)) {
    return false;
  }
  const long double closure_roundoff =
      64.0L * std::numeric_limits<long double>::epsilon() *
      std::max(
          {1.0L, 2.0L * delta_rest_e * delta_rest_e, std::abs(input_delta_m2)});
  if (input_closure2 < -closure_roundoff) {
    return false;
  }
  input_closure2 = std::max(0.0L, input_closure2);
  if (std::sqrt(input_closure2) > 1.0e-9L * std::max(1.0L, hard_mass)) {
    return false;
  }

  Four axis = q1_extended;
  AddScaled(axis, hard, -MinkowskiProduct(q1_extended, hard) / hard_mass2);
  const long double axis_norm2 = -MinkowskiProduct(axis, axis);
  if (!std::isfinite(axis_norm2) || axis_norm2 <= KINEMATICS_EPS) {
    return false;
  }

  const long double half_mass = 0.5L * hard_mass;
  const long double axis_scale = half_mass / std::sqrt(axis_norm2);
  Four k1_extended = hard;
  Four k2_extended = hard;
  Scale(k1_extended, 0.5L);
  Scale(k2_extended, 0.5L);
  AddScaled(k1_extended, axis, axis_scale);
  AddScaled(k2_extended, axis, -axis_scale);
  if (!AllFinite(k1_extended) || !AllFinite(k2_extended)) {
    return false;
  }
  const long double double_max =
      static_cast<long double>(std::numeric_limits<double>::max());
  for (const auto component : k1_extended) {
    if (std::abs(component) > double_max) {
      return false;
    }
  }
  for (const auto component : k2_extended) {
    if (std::abs(component) > double_max) {
      return false;
    }
  }
  k1 = M4Vec(
      static_cast<double>(k1_extended[1]), static_cast<double>(k1_extended[2]),
      static_cast<double>(k1_extended[3]), static_cast<double>(k1_extended[0]));
  k2 = M4Vec(
      static_cast<double>(k2_extended[1]), static_cast<double>(k2_extended[2]),
      static_cast<double>(k2_extended[3]), static_cast<double>(k2_extended[0]));

  const Four projected_k1 = k1.Contravariant<long double>();
  const Four projected_k2 = k2.Contravariant<long double>();
  if (!AllFinite(projected_k1) || !AllFinite(projected_k2)) {
    return false;
  }
  Four projected_delta = Add(projected_k1, projected_k2);
  AddScaled(projected_delta, hard, -1.0L);
  long double projected_closure2 =
      RestFrameNorm2(projected_delta, hard, hard_mass2);
  const long double projected_k1_m2 =
      MinkowskiProduct(projected_k1, projected_k1);
  const long double projected_k2_m2 =
      MinkowskiProduct(projected_k2, projected_k2);
  if (!std::isfinite(projected_closure2) || !std::isfinite(projected_k1_m2) ||
      !std::isfinite(projected_k2_m2)) {
    return false;
  }
  if (projected_closure2 < -closure_roundoff) {
    return false;
  }
  projected_closure2 = std::max(0.0L, projected_closure2);
  const long double shell_tolerance = 1.0e-10L * std::max(1.0L, hard_mass2);
  return k1.E() > 0.0 && k2.E() > 0.0 &&
         std::abs(projected_k1_m2) <= shell_tolerance &&
         std::abs(projected_k2_m2) <= shell_tolerance &&
         std::sqrt(projected_closure2) <= 1.0e-10L * std::max(1.0L, hard_mass);
}

// Kinematic transform in process: p1 + p2 -> {p}, where
//
// p1, p2 are massless spacelike q^2 < 0 => transformed to lightlike q^2 = 0
// {p} massive/massless final states with q^2 = m^2
// Compute false when the on-shell momentum-conserving projection cannot be
// constructed
//
inline bool OffShell2LightCone(M4Vec &p1, M4Vec &p2, std::vector<M4Vec> &p) {
  const int N = p.size();
  const int MAXITER = 100;
  const double STOPEPS = 1e-10;
  if (N == 0) {
    return false;
  }

  const M4Vec input_p1 = p1;
  const M4Vec input_p2 = p2;
  const std::vector<M4Vec> input_p = p;
  auto fail = [&]() {
    p1 = input_p1;
    p2 = input_p2;
    p = input_p;
    return false;
  };

  // 4-momentum sum
  auto psumfunc = [&]() {
    M4Vec sum(0, 0, 0, 0);
    for (const auto &i : aux::indices(p)) {
      sum += p[i];
    }
    return sum;
  };

  // Energy sum
  auto Esumfunc = [&]() {
    double sum = 0.0;
    for (const auto &i : aux::indices(p)) {
      sum += p[i].E();
    }
    return sum;
  };

  /*
    // Fractions
    auto Efrac = [&] () {
      const double Esum = Esumfunc();
      std::vector<double> f(p.size());
      for (const auto& i : aux::indices(p)) { f[i] = p[i].E() / Esum; }
      return f;
    };
  */

  // -------------------------------------------------------------------------
  // Set initial state spacelike (q^2 < 0) particles to lightcone by
  // conserving 3-momentum and thus increasing energy
  p1.SetE(p1.P3mod());
  p2.SetE(p2.P3mod());
  // -------------------------------------------------------------------------

  // Get sum
  const M4Vec q = p1 + p2;

  // Final state masses
  std::vector<double> m(p.size());
  for (const auto &i : aux::indices(p)) {
    const double mass2 = p[i].M2();
    if (!std::isfinite(mass2)) {
      return fail();
    }
    m[i] = std::sqrt(std::max(0.0, mass2));
  }

  int iter = 0;


  while (true) {
    // --------------------------------------------------
    // 3-Momentum scaling and energy subtraction step

    // New energy difference
    const double dE = Esumfunc() - q.E();

    // Current energy fractions
    // std::vector<double> f = Efrac();

    for (const auto &i : aux::indices(p)) {
      const double E = p[i].E() - dE / N; // Scaling distributed by simple 1/N
      // const double E = p[i].E() - dE * f[i]; // Scaling distributed by
      // fractions

      const double p3mod = p[i].P3mod();
      const double rad = E * E - m[i] * m[i];

      if (!std::isfinite(E) || !std::isfinite(p3mod) || !std::isfinite(rad) ||
          p3mod <= KINEMATICS_EPS || rad <= KINEMATICS_EPS) {
        return fail();
      }
      const double a = std::sqrt(rad) / p3mod;
      p[i].Set(a * p[i].Px(), a * p[i].Py(), a * p[i].Pz(), E);
    }

    // --------------------------------------------------
    // 3-Momentum subtraction step

    // New 3-momentum difference
    const M3Vec D3 = (psumfunc() - q).P3();

    // Current energy fractions
    // f = Efrac();

    for (const auto &i : aux::indices(p)) {
      M3Vec D3_this = D3;
      gra::Scale(D3_this,
                 1.0 / N); // Scaling distributed by simple 1/N
      // gra::Scale(D3_this, f[i]); // Scaling distributed by fractions

      p[i].SetP3(gra::Subtract(p[i].P3(), D3_this));
      p[i].SetE(gra::math::msqrt(p[i].P3mod2() + m[i] * m[i]));
    }

    ++iter;
    const M4Vec delta = q - psumfunc();
    const double closure =
        std::max({std::abs(delta.Px()), std::abs(delta.Py()),
                  std::abs(delta.Pz()), std::abs(delta.E())});
    const double scale = std::max({1.0, std::abs(q.Px()), std::abs(q.Py()),
                                   std::abs(q.Pz()), std::abs(q.E())});
    if (closure < STOPEPS * scale) {
      break;
    }
    if (iter >= MAXITER) {
      return fail();
    }

  }

  return true;
}

// A unified Lorentz Transform function to 'frametype' give below:
//
// "CS" : Collins-Soper
// "AH" : Anti-Helicity (Anti-CS)
// "HE" : Helicity
// "PG" : Pseudo-Gottfried-Jackson
//
// Collins-Soper: Quantization z-axis defined by the bi-sector vector between
// initial state p1 and (-p2) (NEGATIVE) directions in the (resonance) system
// rest frame, where p1 and p2 are the initial state proton 3-momentum
//
// Anti-Helicity: Quantization z-axis defined by the bisector vector between
// initial state p1 and (p2) (POSITIVE) directions in the (resonance) system
// rest frame, where p1 and p2 are the initial state proton 3-momentum
//
// Helicity: Quantization axis defined by the resonance system
// 3-momentum vector in the colliding beams frame (lab frame). Use HXFrame() for
// generic helicity frame transforms
//
// Pseudo-Gottfried-Jackson: Quantization axis defined by the initial state
// proton p1 (or p2) 3-momentum vector in the (resonance) system rest frame
//
//
// N.B. For the helicity frame, this function is compatible with symmetric beam
// energies (LHC proton-proton type)
//

// Compute the mass of a validated positive-energy timelike rest system
// M = sqrt(P^2)
inline double ValidatedRestMass(const M4Vec &system,
                                const std::string &context) {
  const double mass2 = system.M2();
  if (!std::isfinite(system.E()) || !(system.E() > 0.0) ||
      !std::isfinite(mass2) || !(mass2 > 0.0)) {
    throw AmplitudeFailure(context +
                           ": rest system is not positive-energy timelike");
  }
  return std::sqrt(mass2);
}

template <typename T>
inline void LorentFramePrepare(const std::vector<T> &p, const T &X,
                               const T &pbeam1, const T &pbeam2, T &pb1boost,
                               T &pb2boost, std::vector<T> &pboost) {
  const double M = ValidatedRestMass(X, "LorentFramePrepare");

  // 2. Boost each particle to the system rest frame
  pboost = p;
  for (const auto &k : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, M, pboost[k], -1); // note minus sign
  }

  // 3. Boost initial state protons
  pb1boost = pbeam1;
  pb2boost = pbeam2;
  gra::kinematics::LorentzBoost(X, M, pb1boost, -1); // note minus sign
  gra::kinematics::LorentzBoost(X, M, pb2boost, -1); // note minus sign
}

// Normalize a frame axis with a deterministic nonzero fallback
// output = axis/|axis|
inline M3Vec NormalizeFrameAxis(const M3Vec &axis, const M3Vec &fallback) {
  const double norm2 = gra::SquaredNorm(axis);
  if (std::isfinite(norm2) && norm2 > KINEMATICS_EPS) {
    const double inverse_norm = 1.0 / math::msqrt(norm2);
    return {axis[0] * inverse_norm, axis[1] * inverse_norm,
            axis[2] * inverse_norm};
  }
  const double fallback_norm2 = gra::SquaredNorm(fallback);
  if (!std::isfinite(fallback_norm2) || fallback_norm2 <= KINEMATICS_EPS) {
    throw std::invalid_argument(
        "gra::kinematics::NormalizeFrameAxis: zero fallback axis");
  }
  const double inverse_norm = 1.0 / math::msqrt(fallback_norm2);
  return {fallback[0] * inverse_norm, fallback[1] * inverse_norm,
          fallback[2] * inverse_norm};
}

// Construct a deterministic unit vector perpendicular to a unit frame axis
// output = normalized(axis cross reference)
inline M3Vec PerpendicularFrameAxis(const M3Vec &axis) {
  const M3Vec zaxis = NormalizeFrameAxis(axis, {0.0, 0.0, 1.0});
  const M3Vec reference =
      std::abs(zaxis[0]) < 0.8 ? M3Vec{1.0, 0.0, 0.0} : M3Vec{0.0, 1.0, 0.0};
  return NormalizeFrameAxis(gra::CrossProduct(zaxis, reference),
                            {0.0, 1.0, 0.0});
}

// Construct a right-handed orthonormal frame including collinear beam limits
// x_hat = y_hat cross z_hat, y_hat = z_hat cross x_hat
inline MMatrix<double> FrameRotation(const M3Vec &beam1, const M3Vec &beam2,
                                     const M3Vec &preferred_z) {
  const M3Vec beam1_unit = NormalizeFrameAxis(beam1, {0.0, 0.0, 1.0});
  const M3Vec beam2_unit = NormalizeFrameAxis(beam2, {0.0, 0.0, -1.0});
  const M3Vec zaxis = NormalizeFrameAxis(preferred_z, beam1_unit);

  const M3Vec beam_normal = gra::CrossProduct(beam1_unit, beam2_unit);
  const M3Vec yaxis =
      NormalizeFrameAxis(beam_normal, PerpendicularFrameAxis(zaxis));
  const M3Vec xaxis = NormalizeFrameAxis(gra::CrossProduct(yaxis, zaxis),
                                         PerpendicularFrameAxis(zaxis));
  const M3Vec corrected_yaxis =
      NormalizeFrameAxis(gra::CrossProduct(zaxis, xaxis), yaxis);
  return {xaxis, corrected_yaxis, zaxis};
}

template <typename T>
inline void LorentzFrame(std::vector<T> &pfout, const T &pb1boost,
                         const T &pb2boost, const std::vector<T> &pfboost,
                         const std::string &frametype, int direction) {
  const M3Vec pb1boost3 = pb1boost.P3();
  const M3Vec pb2boost3 = pb2boost.P3();

  const M3Vec beam1_unit = NormalizeFrameAxis(pb1boost3, {0.0, 0.0, 1.0});
  const M3Vec beam2_unit = NormalizeFrameAxis(pb2boost3, {0.0, 0.0, -1.0});

  // Frame quantization axis before deterministic collinear completion
  M3Vec zaxis{};

  // @@ NON-ROTATED FRAME AXIS DEFINITION @@
  if (frametype == "CM") {
    zaxis = {0, 0, 1};
  }
  // @@ COLLINS-SOPER FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "CS" || frametype == "G1") {
    zaxis = gra::Subtract(beam1_unit, beam2_unit);
  }
  // @@ ANTI-HELICITY FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "AH" || frametype == "G2") {
    zaxis = gra::Add(beam1_unit, beam2_unit);
  }
  // @@ HELICITY FRAME POLARIZATION AXIS DEFINITION @@
  else if (frametype == "HX" || frametype == "G3") {
    zaxis = gra::Negated(gra::Add(pb1boost3, pb2boost3));
  }
  // @@ PSEUDO-GOTTFRIED-JACKSON AXIS DEFINITION: [1] or [2] @@
  else if (frametype == "PG" || frametype == "G4") {
    if (direction == -1) {
      zaxis = pb1boost3;
    } else if (direction == 1) {
      zaxis = pb2boost3;
    } else {
      throw std::invalid_argument(
          "gra::kinematics::LorentzFrame: Invalid direction <" +
          std::to_string(direction) + ">");
    }
  } else {
    throw std::invalid_argument(
        "gra::kinematics::LorentzFrame: Unknown frame <" + frametype + ">");
  }

  // Create SO(3) rotation matrix for the new coordinate axes
  const MMatrix<double> R = frametype == "CM"
                                ? MMatrix<double>(3, 3, "eye")
                                : FrameRotation(pb1boost3, pb2boost3, zaxis);

  // Rotate all vectors
  pfout = pfboost;
  for (const auto &k : gra::aux::indices(pfout)) {
    pfout[k].Rotate(R);
  }
}

// From lab to Anti-Helicity frame
// Quantization z-axis is defined as the bisector of two beams, in the rest
// frame of the system X
//
// Input:     p  =  Set of 4-momentum to be transformed
//            X  =  System 4-momentum
//            p1 =  Beam 1 4-momentum
//            p2 =  Beam 2 4-momentum
//
template <typename T>
inline void AHframe(std::vector<T> &p, const T &X, const T &p1, const T &p2,
                    bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::ANTI-HELICITY FRAME:: \n");
    printf("AHframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  const double rest_mass = ValidatedRestMass(X, "AHframe");

  // Boost particles to the system X rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  T p1b = p1;
  T p2b = p2;

  // Boost the beam particles
  gra::kinematics::LorentzBoost(X, rest_mass, p1b, -1); // Note the minus sign
  gra::kinematics::LorentzBoost(X, rest_mass, p2b, -1); // Note the minus sign

  // Now get the 3-momentum
  const M3Vec pb1boost3 = p1b.P3();
  const M3Vec pb2boost3 = p2b.P3();

  // Anti-Helicity positive bisector vector
  const M3Vec zaxis = gra::Add(NormalizeFrameAxis(pb1boost3, {0.0, 0.0, 1.0}),
                               NormalizeFrameAxis(pb2boost3, {0.0, 0.0, -1.0}));

  // Create SO(3) rotation matrix for the new coordinate axes
  const MMatrix<double> R = FrameRotation(pb1boost3, pb2boost3, zaxis);

  // Rotate all vectors
  for (const auto &k : gra::aux::indices(p)) {
    p[k].Rotate(R);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("AHframe:: Daughters in AH FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************
}

// From lab to Collins-Soper frame
// Quantization z-axis is defined as the bisector of two beams, in the rest
// frame of the system X
//
// Input:     p  =  Set of 4-momentum to be transformed
//            X  =  System 4-momentum
//            p1 =  Beam 1 4-momentum
//            p2 =  Beam 2 4-momentum
//
template <typename T>
inline void CSframe(std::vector<T> &p, const T &X, const T &p1, const T &p2,
                    bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::COLLINS-SOPER FRAME:: \n");
    printf("CSframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  const double rest_mass = ValidatedRestMass(X, "CSframe");

  // Boost particles to the system X rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  T p1b = p1;
  T p2b = p2;

  // Boost the beam particles
  gra::kinematics::LorentzBoost(X, rest_mass, p1b, -1); // Note the minus sign
  gra::kinematics::LorentzBoost(X, rest_mass, p2b, -1); // Note the minus sign

  // Now get the 3-momentum
  const M3Vec pb1boost3 = p1b.P3();
  const M3Vec pb2boost3 = p2b.P3();

  // Collins-Soper bisector vector
  const M3Vec zaxis =
      gra::Subtract(NormalizeFrameAxis(pb1boost3, {0.0, 0.0, 1.0}),
                    NormalizeFrameAxis(pb2boost3, {0.0, 0.0, -1.0}));

  /* ALTERNATIVE WAY, but gives random reflection (rotation around Z by PI)

  M4Vec bijector;
  bijector.SetP3(zaxis);

  // Now get the rotation angle, note the minus
  const double Z_angle = -bijector.Phi();
  const double Y_angle = -bijector.Theta();

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    p[i].RotateZ(Z_angle);
    p[i].RotateY(Y_angle);
    p[i].RotateZ(math::PI);
  }
  */

  // Create SO(3) rotation matrix for the new coordinate axes
  const MMatrix<double> R = FrameRotation(pb1boost3, pb2boost3, zaxis);

  // Rotate all vectors
  for (const auto &k : gra::aux::indices(p)) {
    p[k].Rotate(R);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("CSframe:: Daughters in CS FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************
}

// From lab to Gottfried-Jackson frame
// Quantization z-axis spanned by the propagator (momentum transfer vector)
// momentum in the rest frame of the system X
//
// Input:  p         = Set of 4-momentum to be transformed
//         X         = System 4-momentum
//         direction = -1 or 1
//         q1        = 4-momentum transfer vector 1
//         q2        = 4-momentum transfer vector 2
//
template <typename T>
inline void GJframe(std::vector<T> &p, const T &X, int direction, const T &q1,
                    const T &q2, bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::GOTTFRIED-JACKSON FRAME:: \n");

    printf("GJframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  const double rest_mass = ValidatedRestMass(X, "GJframe");

  // Boost particles to the system X rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  // Boost the propagators
  T q1boost = q1;
  T q2boost = q2;
  gra::kinematics::LorentzBoost(X, rest_mass, q1boost,
                                -1); // Note the minus sign
  gra::kinematics::LorentzBoost(X, rest_mass, q2boost,
                                -1); // Note the minus sign

  // Now get the rotation angle, note the minus
  double Z_angle = 0;
  double Y_angle = 0;
  if (direction == -1) {
    Z_angle = -q1boost.Phi();
    Y_angle = -q1boost.Theta();
  } else if (direction == 1) {
    Z_angle = -q2boost.Phi();
    Y_angle = -q2boost.Theta();
  } else {
    throw std::invalid_argument("GJframe: direction not -1 or 1");
  }

  const auto rotate = [&](T &p) {
    p.RotateZ(Z_angle);
    p.RotateY(Y_angle);
    // p.RotateZ(math::PI); // Reflection
  };

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    rotate(p[i]);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("GJframe:: Daughters in Gottfried-Jackson FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }

    // Rotate propagator
    rotate(q1boost);
    rotate(q2boost);

    printf("GJframe:: Propagators are along z-axis \n");
    printf("GJframe:: Propagator 1 in Gottfried-Jackson FRAME: \n");
    q1boost.Print();

    printf("GJframe:: Propagator 2 in Gottfried-Jackson FRAME: \n");
    q2boost.Print();
  }
  // ********************************************************************
}

// From lab to Pseudo-Gottfried-Jackson frame
// Quantization z-axis spanned by the beam proton +z (or -z) momentum
// in the rest frame of the system X
//
// Input:  p            = Set of 4-momentum to be transformed
//         X            = System 4-momentum
//         direction    = -1 or 1
//         p_beam_plus  = 4-momentum of beam1
//         p_beam_minus = 4-momentum of beam2
//
template <typename T>
inline void PGframe(std::vector<T> &p, const T &X, const int direction,
                    const T &p_beam_plus, const T &p_beam_minus,
                    bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::PSEUDO-GOTTFRIED-JACKSON FRAME:: \n");

    printf("PGframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  const double rest_mass = ValidatedRestMass(X, "PGframe");

  // Boost particles to the system X rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  // Boost initial state protons to the system X rest frame
  T proton_p = p_beam_plus;
  T proton_m = p_beam_minus;

  gra::kinematics::LorentzBoost(X, rest_mass, proton_p,
                                -1); // Note the minus sign
  gra::kinematics::LorentzBoost(X, rest_mass, proton_m,
                                -1); // Note the minus sign

  // Now get the rotation angles, note the minus
  double Z_angle = 0;
  double Y_angle = 0;

  if (direction == -1) {
    Z_angle = -proton_p.Phi();
    Y_angle = -proton_p.Theta();
  } else if (direction == 1) {
    Z_angle = -proton_m.Phi();
    Y_angle = -proton_m.Theta();
  } else {
    throw std::invalid_argument("PGframe: direction not -1 or 1");
  }

  // Rotation function
  const auto rotate = [&](T &p) {
    p.RotateZ(Z_angle);
    p.RotateY(Y_angle);
    // p.RotateZ(math::PI); // Reflection
  };

  // Rotate final states
  for (const auto &i : gra::aux::indices(p)) {
    rotate(p[i]);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("PGframe:: Daughters in Pseudo-Gottfried-Jackson FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }

    // Rotate proton
    rotate(proton_p);

    printf("Protons are on (zx)-plane: \n");
    printf("USER chosen direction along %d beam direction \n", direction);
    printf("  +z Proton in Pseudo-Gottfried-Jackson FRAME: \n");
    proton_p.Print();

    // Rotate other proton
    rotate(proton_m);

    printf("  -z Proton in Pseudo-Gottfried-Jackson FRAME: \n");
    proton_m.Print();

    // Mandelstam invariant
    printf("Mandelstam s = %0.5f \n", (proton_p + proton_m).M());
  }
  // ********************************************************************
}

// Rotate one four-vector so the selected direction becomes the helicity axis
// R = R_z(pi) R_y(-theta_X) R_z(-phi_X)
template <typename T> inline void RotateHelicityAxes(T &p, const T &X) {
  p.RotateZ(-X.Phi());
  p.RotateY(-X.Theta());
  p.RotateZ(math::PI);
}

// From lab to the "Helicity frame"
// Quantization z-axis as the direction spanned by X in the frame of reference
//
// Be careful with the frame definitions in multibody cascaded decays, i.e
// then the intermediate X is typically defined in _its_ own mother frame
//
//
// Input:  p = Set of 4-momentum in the lab to be transformed
//             (N.B. sum over [rotated] p is used to define the boost to their
//             rest frame)
//         X = Helicity direction 4-momentum
//             (N.B. this vector defines only the helicity direction but not the
//             boost here)
template <typename T>
inline void HXframe(std::vector<T> &p, const T &X, bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::HELICITY FRAME:: \n");
    printf("HXframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  for (const auto &i : gra::aux::indices(p)) {
    RotateHelicityAxes(p[i], X);
  }

  // ********************************************************************
  if (DEBUG) {
    printf("HXframe:: Daughters in ROTATED LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }

    T ex(1, 0, 0, 0);
    T ey(0, 1, 0, 0);
    T ez(0, 0, 1, 0);

    // x -> x', y -> y', z -> z'
    RotateHelicityAxes(ex, X);
    RotateHelicityAxes(ey, X);
    RotateHelicityAxes(ez, X);

    printf("HXframe:: AXIS vectors after rotation: \n");

    ex.Print();
    ey.Print();
    ez.Print();
  }
  // ********************************************************************

  // Boost direction defined as a sum over the rotated particles
  T B;
  for (const auto &i : gra::aux::indices(p)) {
    B += p[i];
  }

  const double rest_mass = ValidatedRestMass(B, "HXframe");

  // Boost particles to the system X rest frame
  // -> Helicity frame obtained
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(B, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  // ********************************************************************
  if (DEBUG) {
    printf("HXframe:: Daughters after boost in HELICITY FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }

    printf("HXframe:: Helicity rotation direction X: \n");
    X.Print();

    printf("HXframe:: Boost direction B: \n");
    B.Print();

    printf("\n");
  }
  // ********************************************************************
}

// "Rest frame" (no boost, beam axis as z-axis / spin quantization axis)
// of system X
//
// Input:  p = Set of 4-momentum to be transformed
//         X = System 4-momentum
//
template <typename T>
inline void CMframe(std::vector<T> &p, const T &X, bool DEBUG = false) {
  // ********************************************************************
  if (DEBUG) {
    printf("\n\n ::REST FRAME:: \n");

    printf("CMframe:: Daughters in LAB FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************

  const double rest_mass = ValidatedRestMass(X, "CMframe");

  // Boost particles to the system X rest frame
  for (const auto &i : gra::aux::indices(p)) {
    gra::kinematics::LorentzBoost(X, rest_mass, p[i],
                                  -1); // Note the minus sign
  }

  // ********************************************************************
  if (DEBUG) {
    printf("CMframe:: Daughters in REST FRAME: \n");
    for (const auto &i : gra::aux::indices(p)) {
      p[i].Print();
    }
  }
  // ********************************************************************
}

// Compute one four-vector boosted into a validated timelike rest frame
// Helicity angles require a physical parent rest system
inline M4Vec BoostToRestFrame(const M4Vec &p, const M4Vec &system,
                              const std::string &context) {
  const double mass = ValidatedRestMass(system, context);
  M4Vec out = p;
  LorentzBoost(system, mass, out, -1);
  return out;
}

// Store two momenta and their relative momentum in the pair rest frame
struct PairRestGeometry {
  M4Vec  total;
  M4Vec  u;
  M4Vec  relative_perp;
  M4Vec  first_in_rest;
  double momentum = 0.0;
};

// Compute whether all four generated momentum components are finite
inline bool FiniteFourVector(const M4Vec& momentum) {
  return std::isfinite(momentum.E()) && std::isfinite(momentum.Px()) && std::isfinite(momentum.Py()) &&
         std::isfinite(momentum.Pz());
}

// Build u, r_perp and the first-leg rest-frame momentum from two four-vectors
inline PairRestGeometry PairRest(const M4Vec& first, const M4Vec& second) {
  if (!FiniteFourVector(first) || !FiniteFourVector(second)) {
    throw AmplitudeFailure("PairRest: non-finite generated momentum");
  }
  PairRestGeometry geometry;
  geometry.total     = first + second;
  const double mass2 = geometry.total.M2();
  if (!(mass2 > 0.0) || !std::isfinite(mass2)) { throw AmplitudeFailure("PairRest: central momentum is not timelike"); }
  geometry.u             = geometry.total / std::sqrt(mass2);
  const M4Vec relative   = (first - second) * 0.5;
  geometry.relative_perp = relative - geometry.u * (geometry.u * relative);
  double momentum2       = -geometry.relative_perp.M2();
  if (momentum2 < 0.0 && momentum2 > -KINEMATICS_EPS * std::max(1.0, mass2)) { momentum2 = 0.0; }
  if (momentum2 < 0.0 || !std::isfinite(momentum2)) {
    throw AmplitudeFailure("PairRest: relative momentum is non-spacelike");
  }
  geometry.momentum      = std::sqrt(momentum2);
  geometry.first_in_rest = gra::kinematics::BoostToRestFrame(first, geometry.total, "PairRest");
  if (!FiniteFourVector(geometry.u) || !FiniteFourVector(geometry.relative_perp) ||
      !FiniteFourVector(geometry.first_in_rest)) {
    throw AmplitudeFailure("PairRest: generated geometry is non-finite");
  }
  return geometry;
}

// Boost one particle set into its validated total rest frame
// Two-body production angles are evaluated in the pair rest system
inline void BoostToRestFrame(std::vector<M4Vec> &particles,
                             const std::string &context) {
  M4Vec system;
  for (const auto &particle : particles) {
    system += particle;
  }
  const double mass = ValidatedRestMass(system, context);
  for (auto &particle : particles) {
    LorentzBoost(system, mass, particle, -1);
  }
}

} // namespace kinematics

} // namespace gra

#endif
