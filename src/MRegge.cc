// Regge Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MRegge.h"
#include "Graniitti/MSpin.h"
#include "Graniitti/MEikonal.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"


using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::zi;

using gra::aux::indices;

using namespace gra::form;

namespace gra {
// Generic Regge amplitude initialization
// PDG is the pion/kaon/proton
// TODO: we read couplings based on the first particle,
// generalize to pairs such as pi+pi-K+K- (4-body) (easy extension)
//
void MRegge::InitReggeAmplitude(int PDG, const std::string &MODELPARAM) {
  using json           = nlohmann::json;
  std::string fullpath = gra::aux::GetBasePath(2) + "/modeldata/" + MODELPARAM + "/GENERAL.json";

  // Read and parse
  try {
    std::string data = gra::aux::GetInputData(fullpath);
    json        j    = json::parse(data);

    // Setup the parameter string
    std::string str;
    if        (PDG == 211 || PDG == 111) {    // Charged or Neutral Pions
      str = "PARAM_PI";
    } else if (PDG == 321 || PDG == 311) {    // Charged or Neutral Kaons
      str = "PARAM_K";
    } else if (PDG == 2212 || PDG == 2112) {  // Protons or Neutrons
      str = "PARAM_P";
    } else if (PDG == 113 || PDG == 213) {    // Neutral or charged rho(770)
      str = "PARAM_RHO";
    } else if (PDG == 331) {                  // phi(1020)
      str = "PARAM_PHI"; 
    } else {  // The rest will use PARAM_X setup
      str = "PARAM_X";
    }
    
    // Reggeon couplings
    std::vector<double> c = j.at(str).at("c");
    std::vector<bool>   n = j.at(str).at("n");
    PARAM_REGGE::c        = c;
    PARAM_REGGE::n        = n;

    // PARAM_REGGE::PrintParam();
  } catch (...) {
    throw std::invalid_argument("MRegge::InitReggeAmplitude_: Parse error in " + fullpath);
  }
}

// Triple-Regge (Pomeron) limits
//
// discontinuity line (cut line)
//        |
//        |
//        |
//        .
//        .
// SD:
//
// a =gN= a a =gN= a
//     *      *
//  i   *    *  j
//       *3P*
//      k *  M_X^2
//        *
// a =====gN===== a
//
//
// DD:
//
// a =====gN===== a
//    k1  *
//        *  M_Y^2
//       *3P*
//   i  *    *  j
//      *    *
//       *3P*
//        *  M_X^2
//    k2  *
// a =====gN===== a
//
//
// CD:
//
// a =gN= a  a =gN= a
//     *      *
//  i1  *    *  j1
//       *3P*
//     k  *  M_X^2
//        *
//       *3P*
//  i2  *    *  j2
//     *      *
// a =gN= a  a =gN= a
//
//
//
// [REFERENCE: Gribov, A Reggeon Diagram Technique, Soviet JETP, 1968,
// http://jetp.ac.ru/cgi-bin/dn/e_026_02_0414.pdf]
// [REFERENCE: Gribov, Migdal, Strong coupling in the Pomeranchuk pole problem,
// Sov. Phys. JETP
// 28(4), 784-795 (1968)]
// [REFERENCE: Muller, 1972]
//
// For different forms of triple Pomeron coupling (weak/strong, scalar/vector),
// see:
// [REFERENCE: Luna, Khoze, Martin, Ryskin, https://arxiv.org/abs/1005.4864v1]
//
std::complex<double> MRegge::ME2(gra::LORENTZSCALAR &lts, int mode) const {
  const double s0 = 1.0;  // Energy scale GeV^2

  // Pomeron trajectory intercept + 1
  const double alpha_0 = 1.0 + PARAM_SOFT::DELTA_P;

  // Triple legs
  const double alpha_t_i = S3PomAlpha(lts.t);
  const double alpha_t_j = alpha_t_i;

  std::complex<double> A(0,0);

  if (mode == 1) {  // Elastic

    // Single Pomeron exchange
    // (this calls eikonal pomeron without eikonalization, just for test, use loop on)
    return MEikonal::SingleAmpElastic(lts.s, lts.t);
    
  } else if (mode == 2) {  // SD triple Pomeron

    // Which proton is excited
    const double MX_2 = (lts.ss[1][1] > 1.0) ? lts.ss[1][1] : lts.ss[2][2];

    const double amp2 = pow3(PARAM_SOFT::gN_P)  // Proton-Pomeron coupling ^ 3
                     * pow2(S3F(lts.t))      // Proton form factor ^ 2
                     * PARAM_SOFT::g3P       // Triple-Pomeron coupling

                     * std::pow(lts.s / MX_2,
                                alpha_t_i + alpha_t_j)  // LEFT + RIGHT LEG

                     * std::pow(MX_2 / s0, alpha_0);  // DISCONTINUITY PART

    // recover amplitude, we had the expression at cross section level
    A = gra::form::ReggeEta(alpha_t_i, 1) * msqrt(amp2);

  } else if (mode == 3) {  // DD triple Pomeron, note NO proton form
    // factors (both proton dissociated)

    const double MX_2 = lts.ss[1][1];
    const double MY_2 = lts.ss[2][2];

    const double amp2 = pow2(PARAM_SOFT::gN_P)   // Proton-Pomeron ^ 2
                     * pow2(PARAM_SOFT::g3P)  // Triple-Pomeron coupling ^ 2

                     * std::pow((lts.s * s0) / (MX_2 * MY_2),
                                alpha_t_i + alpha_t_j)  // LEFT + RIGHT LEG

                     * std::pow(MX_2 / s0, alpha_0)   // DISCONTINUITY PART
                     * std::pow(MY_2 / s0, alpha_0);  // DISCONTINUITY PART

    // recover amplitude, we had the expression at cross section level
    A = gra::form::ReggeEta(alpha_t_i, 1) * msqrt(amp2);
  } else {
    throw std::invalid_argument("MRegge::ME2: Unknown mode parameter: " + std::to_string(mode));
  }

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}

// ============================================================================
// Helicity matrix element for Pomeron-Pomeron resonances
//
// [THIS IS UNDER CONSTRUCTION!]
//
//  lambda1    lambda3
//   ==========>
//        *
//        *
//        x----> lambda5
//        *
//        *
//   ==========>
//  lambda2    lambda4
//
//
// [Check Witten-Sakai-Sugimoto model based holographic eta,eta':
// https://arxiv.org/pdf/1406.7010.pdf]
//
std::complex<double> MRegge::ME3HEL(gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const {
  const int J = resonance.p.spinX2 / 2.0;

  // --------------------------------------------------------------------------
  unsigned int       number = 0;
  const unsigned int N      = 16 * (2 * J + 1);

  // C++11 quarantees a local static variable initialization thread safety!
  static std::vector<std::vector<double>> lambda(N, std::vector<double>(5));

  for (double i = 0; i < 2; ++i) {                    // Proton 1
    for (double j = 0; j < 2; ++j) {                  // Proton 2
      for (double k = 0; k < 2; ++k) {                // Proton 3
        for (double l = 0; l < 2; ++l) {              // Proton 4
          for (double m = 0; m < (2 * J + 1); ++m) {  // Central Boson

            lambda[number] = {i - 0.5, j - 0.5, k - 0.5, l - 0.5, m - static_cast<double>(J)};
            ++number;
          }
        }
      }
    }
  }

  // --------------------------------------------------------------------------
  // Proton form factors
  const double FF_A =
      lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
  const double FF_B =
      lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

  // Forward proton deltaphi
  const double phi = lts.pfinal[1].DeltaPhi(lts.pfinal[2]);

  // Now loop over all helicity amplitudes
  // Common factor all amplitudes
  std::complex<double> common = PropOnly(lts.s1, lts.t1) * PARAM_SOFT::gN_P * FF_A *
                                CBW(lts, resonance) * resonance.g * PropOnly(lts.s2, lts.t2) *
                                PARAM_SOFT::gN_P * FF_B;

  lts.hamp.clear();
  unsigned int ok = 0;
  for (std::size_t i = 0; i < N; ++i) {
    /*
    // Apply upper vertex helicity conservation
    if (g_Vertex(lts.t1, lambda[i][0],  lambda[i][2])
                     != std::pow(-1, lambda[i][0] - lambda[i][2]) * g_Vertex(lts.t1,
    lambda[i][0], lambda[i][2]) ) {
                    continue;
    }
    // Apply lower vertex helicity conservation
    if (g_Vertex(lts.t2, lambda[i][1],  lambda[i][3])
                     != std::pow(-1, lambda[i][1] - lambda[i][3]) *
    g_Vertex(lts.t2,-lambda[i][1],-lambda[i][3])) {
                    continue;
    }
    */
    // printf("i = %d :: lambda0 = %0.1f, lambda2 = %0.1f, lambda1 = %0.1f,
    // lambda3 =
    // %0.1f \n",
    //		i, lambda[i][0], lambda[i][2], lambda[i][1], lambda[i][3]);

    // Calculate amplitude
    const std::complex<double> amp =
        g_Vertex(lts.t1, lambda[i][0], lambda[i][2]) *
        gik_Vertex(lts.t1, lts.t2, phi, lambda[i][4], resonance.p.spinX2 / 2.0, resonance.p.P) *
        common;
    g_Vertex(lts.t2, lambda[i][1], lambda[i][3]);

    // std::cout << amp << " :: " << gra::math::abs2(amp) << std::endl;

    // Is non-zero
    if (gra::math::abs2(amp) > 0) {
      lts.hamp.push_back(amp);
      ++ok;
    }
  }
  if (ok == 0) {
    throw std::invalid_argument("MRegge::ME3HEL: All " + std::to_string(N) +
                                " helicity amplitudes are zero! \n");
  }

  // 1/4 {|H1|^2 + |H2|^2 + ... + |HN|^2}
  double amp2 = 0.0;
  for (std::size_t i = 0; i < lts.hamp.size(); ++i) { amp2 += gra::math::abs2(lts.hamp[i]); }
  amp2 /= 4;  // Initial state average

  return msqrt(amp2); // we expect amplitude
}

// This function is used for parity conservation check:
//
// [ gamma_{m1,m2}^\lambda_h = (-1)^\lambda_h \xi_3 \gamma_{-m1-m2}^-\lambda ]
//
// Pomeron i  (Parity P, Signature sigma)
// Pomeron k  (Parity P, Signature sigma)
// Resonance  (Spin J, Parity P)
//
// "Naturality" \def= Reggeon Parity x Reggeon Regge Signature
//
int MRegge::xi3(int J, int P, int P_i, int sigma_i, int P_k, int sigma_k) const {
  return P * std::pow(-1, J) * P_i * sigma_i * P_k * sigma_k;
}

// --------------------------------------------------------------------------------------
// Proton-Pomeron-Proton vertex with helicity transition lambda_i-> lambda_f
//
// initial > final
//    -1/2 > -1/2
//    -1/2 >  1/2
//     1/2 > -1/2
//     1/2 >  1/2
//
// ~ works for for small |t|
//
// g_{\lambda_1,\lambda_3}(t) = (-1)^{\lambda_1 - \lambda_3} \xi_1
// g_{-\lambda_1-\lambda_3}(t)
// \xi_1 = \eta_1\eta_3 (-1)^{S_1 - S_3} P_i \sigma_i
//
double MRegge::g_Vertex(double t, double lambda_i, double lambda_f) const {
  return std::pow(-t, std::abs(lambda_i - lambda_f) / 2);
}

// Central vertex subfunction: For small |t1|, |t2|
//
// Selection rules:
// 1. m1 + m2 = lambda_h (angular momentum)
// 2. parity conservation
//
double MRegge::gammaLambda(double t1, double t2, double m1, double m2) const {
  return std::pow(-t1, std::abs(m1) / 2) * std::pow(-t2, std::abs(m2) / 2);
}

// Central vertex g_ik(t1,t2,\phi)^{\lambda_h} (central resonance J, P)
//
// Selection rules:
// 1. lambda_h = m1 + m2 (angular momentum)
// 2. parity conservation
//
// Naturality: intrinsic parity x (-1)^J (bosons)
//             intrinsic parity x (-1)^(J-1/2) (fermions)
//
//
std::complex<double> MRegge::gik_Vertex(double t1, double t2, double phi, int lambda_h, int J,
                                        int P) const {
  // Loop over Pomeron helicity m2, implicitly loop over m1 (see below)
  const int MAX_m = 1;

  // (we handle only Pomerons here for now, no secondary Reggeons - trivial
  // extension)
  const int pom1_P     = 1;  // Pomeron 1 parity +
  const int pom1_sigma = 1;  // Pomeron 1 signature +
  const int pom2_P     = 1;  // Pomeron 2 parity +
  const int pom2_sigma = 1;  // Pomeron 2 signature +

  const int            xi3val = xi3(J, P, pom1_P, pom1_sigma, pom2_P, pom2_sigma);
  std::complex<double> sum    = 0.0;

  for (int m2 = -MAX_m; m2 <= MAX_m; ++m2) {  // -\inf < m_2 < inf

    // Helicity relation lambda_h = m1 + m2
    const int m1 = lambda_h - m2;

    // Parity conservation
    const int lhs = 1;
    const int rhs = std::pow(-1, lambda_h) * xi3val;

    if (lhs != rhs) continue;

    // Analytical continuation
    std::complex<double> val = gra::math::sign(m1 * m2) *
                               std::exp(gra::math::zi * static_cast<double>(m2) * phi) *
                               gammaLambda(t1, t2, m1, m2);
    sum += val;
  }
  // std::cout << std::endl;
  return sum;
}

// ============================================================================
// Regge matrix element ansatz for 2->4 continuum spectrum
// with sign, one can choose the relative sign between sub-t and sub-u
// amplitudes -> e.g. baryons (fermions) (could) have different effective sign
// than mesons (bosons).
//
//  ======F======>
//        *
//        * P
//        *
//        ff ----> pi+/K+,p,rho0 ...
//        |
//        ff ----> pi-/K-,pbar,rho0 ...
//        *
//        * P
//        *
//  ======F======>
//
// For similar Regge constructions, see:
//
// [REFERENCE: Azimov, Khoze, Levin, Ryskin, Sov.J.Nucl.Phys.21, 215 (1975)]
// [REFERENCE: Pumplin, Henyey, Nucl.Phys. B117, 377 (1976)]
// [REFERENCE: Lebiedowicz, Szczurek, https://arxiv.org/abs/0912.0190]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, https://arxiv.org/abs/1312.4553]
//
std::complex<double> MRegge::ME4(gra::LORENTZSCALAR &lts, double sign) const {
  // Offshell propagator nominal mass^2 read from the outgoing particle PDG mass
  const double M2_ = pow2(lts.decaytree[0].p.mass);

  // Create function pointers
  double (*ff)(double, double);
  double (*prop)(double, double);

  if (sign > 0) {  // positive sign amplitudes

    ff   = &PARAM_REGGE::Meson_FF;
    prop = &PARAM_REGGE::Meson_prop;
  } else {  // negative sign (alternative model spin-statistics) amplitude

    ff   = &PARAM_REGGE::Baryon_FF;
    prop = &PARAM_REGGE::Baryon_prop;
  }

  const double FF_A =
      lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
  const double FF_B =
      lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

  // Particle-Particle-Pomeron coupling
  const double gpp_P = PARAM_REGGE::c[0] / PARAM_SOFT::gN_P;

  const std::complex<double> A_t = PropOnly(lts.ss[1][3], lts.t1) * PARAM_SOFT::gN_P * FF_A *
                                   (*ff)(lts.t_hat, M2_) * gpp_P * prop(lts.t_hat, M2_) *
                                   (*ff)(lts.t_hat, M2_) * gpp_P * PropOnly(lts.ss[2][4], lts.t2) *
                                   PARAM_SOFT::gN_P * FF_B;

  // sign applied here
  const std::complex<double> A_u = sign * PropOnly(lts.ss[1][4], lts.t1) * PARAM_SOFT::gN_P * FF_A *
                                   (*ff)(lts.u_hat, M2_) * gpp_P * prop(lts.u_hat, M2_) *
                                   (*ff)(lts.u_hat, M2_) * gpp_P * PropOnly(lts.ss[2][3], lts.t2) *
                                   PARAM_SOFT::gN_P * FF_B;

  // Total amplitude
  const std::complex<double> A = A_t + A_u;

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}

// ============================================================================
// Regge matrix element ansatz for 2->6 continuum spectrum
//
//
//  ======F======>
//        *
//        *
//        ff-----> a
//        |
//        ff-----> b
//        *
//        *
//        ff-----> c
//        |
//        ff-----> d
//        *
//        *
//  ======F======>
//
// For many similar amplitudes, see the literature of the era
// "pre-superstring/generalized Veneziano amplitude". For example:
//
// [REFERENCE: Bardakci, Ruegg,
// https://journals.aps.org/pr/pdf/10.1103/PhysRev.181.1884]
//
//
// Numerically, one may compare with:
//
// [REFERENCE: Kycia, Lebiedowicz, Szczurek, Turnau,
// https://arxiv.org/abs/1702.07572]
//
//
std::complex<double> MRegge::ME6(gra::LORENTZSCALAR &lts) const {

  // Amplitude_
  std::complex<double> A(0, 0);

  // Offshell propagator masses^2 [for mixed states, such as pi+ pi- K+ K-]
  const double M2_A = pow2(lts.decaytree[0].p.mass);
  const double M2_B = pow2(lts.decaytree[2].p.mass);

  // Create function pointers
  double (*ff)(double, double);
  double (*prop)(double, double);

  const double sign = 1;

  if (sign > 0) {  // positive sign amplitudes

    ff   = &PARAM_REGGE::Meson_FF;
    prop = &PARAM_REGGE::Meson_prop;
  } else {  // negative sign (alternative model spin-statistics) amplitude

    ff   = &PARAM_REGGE::Baryon_FF;
    prop = &PARAM_REGGE::Baryon_prop;
  }

  const double FF_A =
      lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
  const double FF_B =
      lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

  // Particle-Particle-Pomeron coupling
  const double gpp_P = PARAM_REGGE::c[0] / PARAM_SOFT::gN_P;

  // Loop over different final state permutations (max #16)
  for (const auto &i : indices(permutations4_)) {
    const unsigned int a = permutations4_[i][0];
    const unsigned int b = permutations4_[i][1];
    const unsigned int c = permutations4_[i][2];
    const unsigned int d = permutations4_[i][3];

    // Calculate t-type Lorentz scalars here [no need, done already]
    // const double tt_ab = (pbeam1_pfinal1 - pfinal[a]).M2();
    // const double tt_bc = (pbeam1_pfinal1 - pfinal[a] - pfinal[b]).M2(); // by
    // symmetry, same as (pbeam2 - pfinal[2] - pfinal[d] - pfinal[c]).M2()
    // const double tt_cd = (pbeam2_pfinal2 - pfinal[d]).M2();

    // Get Lorentz scalars
    const double tt_ab = lts.tt_1[a];
    const double tt_bc = lts.tt_xy[a][b];
    const double tt_cd = lts.tt_2[d];

    const std::complex<double> subamp =
        PropOnly(lts.ss[1][a], lts.t1) * PARAM_SOFT::gN_P * FF_A * (*ff)(tt_ab, M2_A) * gpp_P *
        prop(tt_ab, M2_A) * (*ff)(tt_bc, M2_A) * gpp_P * PropOnly(lts.ss[b][c], tt_bc) *
        (*ff)(tt_bc, M2_B) * gpp_P * prop(tt_cd, M2_B) * (*ff)(tt_cd, M2_B) * gpp_P *
        PropOnly(lts.ss[2][d], lts.t2) * PARAM_SOFT::gN_P * FF_B;

    A += subamp;
  }

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}

// ============================================================================
// Regge matrix element ansatz for 2->8 continuum spectrum (generalization of
// 2->6)
//
//  ======F======>
//        *
//        *
//        ff-----> a
//        |
//        ff-----> b
//        *
//        *
//        ff-----> c
//        |
//        ff-----> d
//        *
//        *
//        ff-----> e
//        |
//        ff-----> f
//        *
//        *
//  ======F======>
//
std::complex<double> MRegge::ME8(gra::LORENTZSCALAR &lts) const {
  const std::complex<double> nullc(0, 0);

  // Amplitude
  std::complex<double> A(0, 0);

  // Offshell propagator masses^2 [for mixed states, such as pi+ pi- K+ K- pi+
  // pi-]
  const double M2_A = pow2(lts.decaytree[0].p.mass);
  const double M2_B = pow2(lts.decaytree[2].p.mass);
  const double M2_C = pow2(lts.decaytree[4].p.mass);

  // Create function pointers
  double (*ff)(double, double);
  double (*prop)(double, double);

  const double sign = 1;

  if (sign > 0) {  // positive sign amplitudes

    ff   = &PARAM_REGGE::Meson_FF;
    prop = &PARAM_REGGE::Meson_prop;
  } else {  // negative sign (alternative model spin-statistics) amplitude

    ff   = &PARAM_REGGE::Baryon_FF;
    prop = &PARAM_REGGE::Baryon_prop;
  }
  const double FF_A =
      lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
  const double FF_B =
      lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

  // Particle-Particle-Pomeron coupling
  const double gpp_P = PARAM_REGGE::c[0] / PARAM_SOFT::gN_P;

  // Loop over different permutations (max #288)
  for (const auto &i : indices(permutations6_)) {
    const unsigned int a = permutations6_[i][0];
    const unsigned int b = permutations6_[i][1];
    const unsigned int c = permutations6_[i][2];
    const unsigned int d = permutations6_[i][3];
    const unsigned int e = permutations6_[i][4];
    const unsigned int f = permutations6_[i][5];

    // t-type Lorentz scalars [no need here, already calculated]
    // const double tt_ab = (pbeam1_pfinal1 - pfinal[a]).M2();
    // const double tt_bc = (pbeam1_pfinal1 - pfinal[a] - pfinal[b]).M2();
    // const double tt_cd = (pbeam1_pfinal1 - pfinal[b] - pfinal[c]).M2();
    // const double tt_de = (pbeam1_pfinal1 - pfinal[c] - pfinal[d]).M2();
    // const double tt_ef = (pbeam2_pfinal2 - pfinal[f]).M2();

    // Collect scalars
    const double tt_ab = lts.tt_1[a];
    const double tt_bc = lts.tt_xy[a][b];
    const double tt_cd = lts.tt_xy[b][c];
    const double tt_de = lts.tt_xy[c][d];
    const double tt_ef = lts.tt_2[f];

    const std::complex<double> subamp =
        PropOnly(lts.ss[1][a], lts.t1) * PARAM_SOFT::gN_P * FF_A * (*ff)(tt_ab, M2_A) * gpp_P *
        prop(tt_ab, M2_A) * (*ff)(tt_bc, M2_A) * gpp_P * PropOnly(lts.ss[b][c], tt_bc) *
        (*ff)(tt_bc, M2_B) * gpp_P * prop(tt_cd, M2_B) * (*ff)(tt_de, M2_B) * gpp_P *
        PropOnly(lts.ss[d][e], tt_de) * (*ff)(tt_de, M2_C) * gpp_P * prop(tt_ef, M2_C) *
        (*ff)(tt_ef, M2_C) * gpp_P * PropOnly(lts.ss[2][f], lts.t2) * PARAM_SOFT::gN_P * FF_B;

    A += subamp;
  }

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}

// ----------------------------------------------------------------------------
// Generic resonance (sub)-cross section:
//
//        \hat{\sigma} ~ \pi / \hat{s} [W_i(\hat{s} W_f(\hat{s})] /
//                                     [(\hat{s} - M_0)^2 + W_f(\hat{s})^2 ]
//
// where W_{i,f} denotes effective initial (final) state factors
// incorporating spin and color symmetry factors, couplings and the rest of the
// dynamics
//
//
// In general, the behavior can be different in domains:
//
// 1. \hat{s} << M0^2, where M0 is the on-shell mass
// 2. \hat{s} ~  M0^2
// 3. \hat{s} >> M0^2
//
// and isolated resonances might break unitarity when \hat{s} -> inf,
// but resonance + continuum amplitude does not.

// ============================================================================
// Simple matrix element ansatz for Pomeron-Pomeron resonances
//
// ===========
//      *
//      *
//      x---->
//      *
//      *
// ===========
//
std::complex<double> MRegge::ME3(gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const {
  // Proton form factors
  const double FF_A =
      lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
  const double FF_B =
      lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

  // s-channel
  // Factor 2 x from initial state (identical boson) statistics
  const std::complex<double> A_prod =
      2.0 * PropOnly(lts.s1, lts.t1) * FF_A * PARAM_SOFT::gN_P * CBW(lts, resonance) *
      PARAM_REGGE::ResonanceFormFactor(lts.m2, pow2(resonance.p.mass), resonance.g_FF) *
      PropOnly(lts.s2, lts.t2) * FF_B * PARAM_SOFT::gN_P;

  // Production and Decay amplitude
  const std::complex<double> A_spin = spin::ProdAmp(lts, resonance) * spin::SpinAmp(lts, resonance);

  // Flux
  // const double V = msqrt(lts.s / lts.s1 / lts.s2);
  const double V = msqrt(1.0 / lts.m2);

  // Full amplitude
  const std::complex<double> A = A_prod * A_spin * V;

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}

// ============================================================================
// Simple matrix element ansatz for Photoproduction of resonances
//
// ===========
//      $
//      $
//      x---->
//      *
//      *
// ===========
//
std::complex<double> MRegge::PhotoME3(gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const {
  // Check spin
  if (resonance.p.spinX2 != 2) {
    throw std::invalid_argument("MRegge::PhotoME3: Resonance spin J = " +
                                std::to_string(resonance.p.spinX2) + " (should be J = 1)!");
  }

  // Resonance part
  const std::complex<double> common =
      CBW(lts, resonance) *
      PARAM_REGGE::ResonanceFormFactor(lts.m2, pow2(resonance.p.mass), resonance.g_FF);

  double gammaflux1 = lts.excite1 ? IncohFlux(lts.x1, lts.t1, lts.qt1, lts.pfinal[1].M2())
                                  : CohFlux(lts.x1, lts.t1, lts.qt1);
  double gammaflux2 = lts.excite2 ? IncohFlux(lts.x2, lts.t2, lts.qt2, lts.pfinal[2].M2())
                                  : CohFlux(lts.x2, lts.t2, lts.qt2);

  // "To amplitude level"
  gammaflux1 = msqrt(gammaflux1);
  gammaflux2 = msqrt(gammaflux2);

  // Photon up (t1) x Pomeron down (t2)
  const std::complex<double> M1 =
      gammaflux1 / msqrt(lts.x1) * common *
      PhotoProp(lts.s2, lts.t2, pow2(resonance.p.mass), lts.excite2, lts.pfinal[2].M2());

  // Pomeron up (t1) x Photon down (t2)
  const std::complex<double> M2 =
      gammaflux2 / msqrt(lts.x2) * common *
      PhotoProp(lts.s1, lts.t1, pow2(resonance.p.mass), lts.excite1, lts.pfinal[1].M2());

  // Should sum here with negative sign if proton-antiproton initial state (anti-symmetric)
  const std::complex<double> A_prod = (lts.beam1.pdg == lts.beam2.pdg) ? (M1 + M2) : (M1 - M2);
  
  // Production and Decay amplitude
  const std::complex<double> A_spin = spin::ProdAmp(lts, resonance) * spin::SpinAmp(lts, resonance);

  // Full amplitude
  const std::complex<double> A = A_prod * A_spin;

  // --------------------------------------------------------------------
  // For screening loop
  lts.hamp = {A};
  // --------------------------------------------------------------------

  return A;
}


// ============================================================================
// Pomeron propagator x coupling to proton in photoproduction
//
// Essentially, parametrized Regge scaling of
// quasielastic photon + proton -> vectormeson + proton
//
// dsigma(vector meson)/dt ~ (W_yp^2/W0^2)^{2(alpha(t) - 1)} x exp(B_0 t),
//
// where the trajectory follows
// - soft pomeron for rho (delta ~ 0.08)
// - hard pomeron for j/psi (delta ~ 0.2) ... upsilon (delta ~ 0.3)
//
// HERA notation uses usually Q^2 for the lepton side momentum transfer
// and t for the proton (pomeron) side.
//
//
// [REFERENCE: Donnachie, Dosch, Landshoff, Nachtmann, Pomeron Physics and QCD, 2002]

std::complex<double> MRegge::PhotoProp(double s, double t, double m2, bool excite,
                                       double M2_forward) const {
  // Different vector meson slopes (HERA data),
  // dsigma/dt ~ exp(Bt),
  // where B = B0 + 4alpha'log(W/W0) = B0 + 2alpha'log(W^2/W0^2)
  // W0 reference scale usually 90...95 GeV at which B0 is fitted from
  // data
  // B0 ~ 11 for rho and 4.8 for J/psi, Upsilon
  //
  // A running slope fit to vector meson / photoproduction slope data
  const double m_freeze = 0.7;  // GeV, Minimum mass freezing of running
  const double B0 =
      std::pow(std::max(pow2(m_freeze), m2) / 2.0, -gra::math::PI / 2.0) + 3 * gra::math::PI / 2.0;

  // Normalization scale GeV^2
  // const double W02 = pow2(90.0); // Typical HERA data normalization
  const double W02 = pow2(1.0);  // We normalize to 1

  // N.B. The pomeron slope (a') also affects, so all parameters need
  // to be fitted in a proper way (full MC)
  const double alpha = PARAM_REGGE::a0[0] + t * PARAM_REGGE::ap[0];  // Pomeron trajectory
  const std::complex<double> eta = ReggeEta(alpha, 1);               // Pomeron signature

  // Proton form factor simply exponential here
  // Division by 2, because we are at amplitude level
  const double FF = excite ? gra::form::S3FINEL(t, M2_forward) : std::exp(B0 * t / 2.0);

  return eta * std::pow(s / W02, alpha) * FF * PARAM_SOFT::gN_P;
}

// ============================================================================
// Regge propagators in the form (s/s0)^alpha(t)
//
std::complex<double> MRegge::PropOnly(double s, double t) const {
  std::complex<double> M(0, 0);

  // Loop over Pomeron, Reggeon, Reggeon... exchanges
  for (std::size_t k = 0; k < PARAM_REGGE::a0.size(); ++k) {
    if (PARAM_REGGE::n[k] == true) {  // Trajectory active

      double alpha = 0.0;

      if (k == 0) {  // Pomeron
        // alpha = S3PomAlpha(t); // Non-Linear trajectory
        alpha = PARAM_REGGE::a0[k] + t * PARAM_REGGE::ap[k];
      } else {  // Reggeons
        alpha = PARAM_REGGE::a0[k] + t * PARAM_REGGE::ap[k];
      }

      // Complex Regge signature
      // (cannot be used with linear trajectory due to poles)
      // const std::complex<double> eta = ReggeEta(alpha, PARAM_REGGE::sgn[k]);

      // Evaluate at t = 0 to avoid poles (for linear trajectory)
      const std::complex<double> eta = ReggeEta(PARAM_REGGE::a0[k], PARAM_REGGE::sgn[k]);

      // Add to the sum
      M += eta * std::pow(s / PARAM_REGGE::s0, alpha);
    }
  }
  return M;
}

// ----------------------------------------------------------------------

// Regge amplitude parameters
namespace PARAM_REGGE {
bool initialized = false;

std::vector<double> a0;
std::vector<double> ap;
std::vector<double> sgn;

double s0 = 0.0;

int    offshellFF = 0;
double b_EXP      = 0.0;
double a_OREAR    = 0.0;
double b_OREAR    = 0.0;
double b_POW      = 0.0;

bool reggeize = false;

std::vector<double> c;  // coupling
std::vector<bool>   n;  // on/off

void PrintParam() {
  std::cout << "PARAM_REGGE:: Sub-amplitude parameters:" << std::endl << std::endl;
  for (unsigned int i = 0; i < a0.size(); ++i) {
    printf(
        "- Reggeon[%d]: trajectory alpha(t) = %0.3f + %0.3f "
        "[GeV^{-2}] t, sgn = %0.0f \n",
        i, a0[i], ap[i], sgn[i]);
  }
  printf("- offshellFF = %d \n", offshellFF);
  if (offshellFF == 0) { printf("  -- b_EXP = %0.3f [GeV^{-2}] (exponential) \n", b_EXP); }
  if (offshellFF == 1) {
    printf(
        "  -- a_OREAR = %0.3f [GeV^{-1}], b_OREAR = %0.3f "
        "[GeV^{-1}] (Orear)\n",
        a_OREAR, b_OREAR);
  }
  if (offshellFF == 2) { printf("  -- b_POW = %0.3f [GeV^{-2}] (powerlaw) \n", b_POW); }
  printf("- reggeize = %d \n", reggeize);
  printf("- s0 = %0.3f [GeV^2]", s0);
  std::cout << std::endl << std::endl;

  std::cout << "Couplings:" << std::endl;
  for (unsigned int i = 0; i < c.size(); ++i) {
    printf("- Reggeon[%d]: c = %0.3f [GeV^{-2}]", i, c[i]);
    if (n[i]) {
      std::cout << rang::fg::green << " [on]" << rang::fg::reset << std::endl;
    } else {
      std::cout << rang::fg::red << " [off]" << rang::fg::reset << std::endl;
    }
  }
  std::cout << std::endl;
}

// Non-Conserved Vector Current Pomeron [UNTESTED FUNCTION]
//
// [REFERENCE: Close, Schuler, https://arxiv.org/pdf/hep-ph/9902243.pdf]
// [REFERENCE: Close, Schuler, https://arxiv.org/pdf/hep-ph/9905305.pdf]
//
// see also:
// [REFERENCE: Kaidalov, Khoze, Martin, Ryskin,
// https://arxiv.org/pdf/hep-ph/0307064.pdf]
//
std::complex<double> JPC_CS_coupling(const gra::LORENTZSCALAR &lts,
                                     const gra::PARAM_RES &    resonance) {
  const unsigned int PP = 0;  // +1 +1
  const unsigned int PM = 1;  // +1 -1
  const unsigned int LL = 2;  //  0  0

  const unsigned int PL = 3;  // +1  0
  const unsigned int LP = 4;  //  0 +1

  std::vector<double> A(5, 0.0);

  // -------------------------------------------------------------------

  const int    J   = resonance.p.spinX2 / 2.0;
  const int    P   = resonance.p.P;
  const double phi = lts.pfinal[1].DeltaPhi(lts.pfinal[2]);

  // Dynamic structure (approx)
  const double Q1_2  = std::abs(lts.t1);
  const double Q2_2  = std::abs(lts.t2);
  const double D     = (Q1_2 - Q2_2) / lts.m2;
  const double delta = msqrt(Q1_2) * msqrt(Q2_2) / lts.m2;

  // Naturality
  const double eta1 = 1;                    // Pomeron 1 naturality (parity x signature)
  const double eta2 = 1;                    // Pomeron 2 naturality (parity x signature)
  const double etaM = P * std::pow(-1, J);  // Central boson naturality
  const double eta  = eta1 * eta2 * etaM;

  // -------------------------------------------------------------------

  // Coupling structure (capture local variables with &)
  auto Q = [&]() { return lts.t1 * lts.t2 * pow2(A[PP]) * pow2(std::sin(phi)); };
  auto W = [&]() {
    const double xi1 = gra::math::csgn(A[PP]) * gra::math::csgn(A[LL]);
    return pow2(msqrt(lts.t1 * lts.t2) * A[LL] -
                xi1 * msqrt(lts.t1 * lts.t2) * A[PP] * std::cos(phi));
  };
  auto E = [&]() {
    const double xi2   = gra::math::csgn(A[PL]) * gra::math::csgn(A[LP]);
    const double alpha = lts.qt2 * A[PL];
    const double beta  = eta * xi2 * lts.qt1 * A[LP];
    return (lts.q1 * alpha - lts.q2 * beta).Pt2();
  };
  auto R = [&]() { return lts.t1 * lts.t2 * 0.5 * pow2(A[PM]); };

  // -------------------------------------------------------------------

  const double nabla = 0;

  // Selection rules
  if (J == 0 && P == -1) {
    A[LL] = 0;
    A[PP] = 1;
    A[PM] = 0;

    return msqrt(Q());
  }
  if (J == 0 && P == 1) {
    A[LL] = nabla * delta;
    A[PP] = 1;
    A[PM] = 0;
    return msqrt(W());
  }
  if (J == 1 && P == -1) {
    A[LL] = nabla * D * delta;
    A[PP] = D;
    A[PM] = 0;

    return msqrt(Q() + W());
  }
  if (J == 1 && P == 1) {
    A[LL] = 0;
    A[PP] = D;
    A[PM] = 0;

    return msqrt(Q() + W());
  }
  if (J == 2 && P == -1) {
    A[LL] = 0;
    A[PP] = 1;
    A[PM] = D;
    A[PL] = msqrt(Q2_2 / lts.m2);
    A[LP] = msqrt(Q1_2 / lts.m2);

    return msqrt(Q() + E() + R());
  }
  if (J == 2 && P == 1) {
    A[LL] = nabla * delta;
    A[PP] = 1;
    A[PM] = 1;
    A[PL] = msqrt(Q2_2 / lts.m2);
    A[LP] = msqrt(Q1_2 / lts.m2);

    return msqrt(W() + E() + R());
  }

  return 1.0;
}

// For data, see:
//
// [REFERENCE: Kirk for WA102, https://arxiv.org/abs/hep-ph/9908253v1]

// Off-Shell meson form factor
double Meson_FF(double t_hat, double M2) {
  const double tprime = t_hat - M2;

  switch (PARAM_REGGE::offshellFF) {
    // Exponential
    case 0:
      return std::exp(-PARAM_REGGE::b_EXP * std::abs(tprime));
    // Orear
    case 1:
      return std::exp(-PARAM_REGGE::b_OREAR * msqrt(std::abs(tprime) + pow2(PARAM_REGGE::a_OREAR)) +
                      PARAM_REGGE::a_OREAR * PARAM_REGGE::b_OREAR);
    // Powerlaw
    case 2:
      return 1.0 / (1.0 - tprime / PARAM_REGGE::b_POW);

    default:
      throw std::invalid_argument("Meson_FF: Unknown parameter: " +
                                  std::to_string(PARAM_REGGE::offshellFF));
  }
}

// Off-Shell baryon form factor
double Baryon_FF(double t_hat, double M2) {
  const double tprime = t_hat - M2;

  switch (PARAM_REGGE::offshellFF) {
    // Exponential
    case 0:
      return std::exp(-PARAM_REGGE::b_EXP * std::abs(tprime));
    // Orear
    case 1:
      return std::exp(-PARAM_REGGE::b_OREAR * msqrt(std::abs(tprime) + pow2(PARAM_REGGE::a_OREAR)) +
                      PARAM_REGGE::a_OREAR * PARAM_REGGE::b_OREAR);
    // Powerlaw
    case 2:
      return 1.0 / (1.0 - tprime / PARAM_REGGE::b_POW);

    default:
      throw std::invalid_argument("Baryon_FF: Unknown parameter: " +
                                  std::to_string(PARAM_REGGE::offshellFF));
  }
}

// Off-shell meson M* propagator
double Meson_prop(double t_hat, double M2) {
  double val = 1.0 / (t_hat - M2);
  if (PARAM_REGGE::reggeize) {
    //"Reggeization" ~ s^J(t) (NOT IMPLEMENTED, return without)
    return val;
  } else {
    return val;
  }
}

// Off-shell baryon B* propagator
double Baryon_prop(double t_hat, double M2) {
  double val = 1.0 / (t_hat - M2);
  if (PARAM_REGGE::reggeize) {
    //"Reggeization" ~ s^J(t) (NOT IMPLEMENTED, return without)
    return val;
  } else {
    return val;
  }
}

// Final State Interaction (FSI) kt^2 loop integral
// NOT IMPLEMENTED
std::complex<double> FSI_prop(double t_hat, double M2) {
  std::complex<double> value(1, 0);
  return value;
}

// Resonance form factor
// (in addition to Breit-Wigner propagator)
//
// In general, without any interfering continuum process
// (which would provide the same effect), this may be needed?
//
double ResonanceFormFactor(double s_hat, double M2, bool active) {
  if (!active) {
    return 1.0;
  } else {
    const double scale4 = 1.0;  // GeV^4
    return std::exp(-pow2(s_hat - M2) / scale4);
  }
}

}  // PARAM_REGGE namespace ends
}  // gra namespace ends
