// Tensor Pomeron amplitudes
//
// See:
// [REFERENCE: Ewerz, Maniatis, Nachtmann, https://arxiv.org/pdf/1309.3478.pdf]
// [REFERENCE: Lebiodowicz, Nachtmann, Szczurek, arxiv.org/pdf/1601.04537.pdf]
// [REFERENCE: LNS, arxiv.org/pdf/1801.03902.pdf]
// [REFERENCE: LNS, arxiv.org/pdf/1901.11490.pdf]
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <iostream>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/M4Vec.h"
#include "Graniitti/MDirac.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MRegge.h"
#include "Graniitti/MTensorPomeron.h"

// FTensor algebra
#include "FTensor.hpp"

// LOOP MACROS
#define FOR_EACH_2(X)           \
  for (const auto &u : X) {     \
    for (const auto &v : X)     {
#define FOR_EACH_2_END \
  }                    \
  }

#define FOR_EACH_3(X)           \
  for (const auto &u : X) {     \
    for (const auto &v : X) {   \
      for (const auto &k : X)   {

#define FOR_EACH_3_END \
  }                    \
  }                    \
  }

#define FOR_EACH_4(X)           \
  for (const auto &u : X) {     \
    for (const auto &v : X) {   \
      for (const auto &k : X) { \
        for (const auto &l : X) {
#define FOR_EACH_4_END \
  }                    \
  }                    \
  }                    \
  }

#define FOR_EACH_5(X)             \
  for (const auto &u : X) {       \
    for (const auto &v : X) {     \
      for (const auto &k : X) {   \
        for (const auto &l : X) { \
          for (const auto &r : X) {
#define FOR_EACH_5_END \
  }                    \
  }                    \
  }                    \
  }                    \
  }

#define FOR_EACH_6(X)               \
  for (const auto &u : X) {         \
    for (const auto &v : X) {       \
      for (const auto &k : X) {     \
        for (const auto &l : X) {   \
          for (const auto &r : X) { \
            for (const auto &s : X) {
#define FOR_EACH_6_END \
  }                    \
  }                    \
  }                    \
  }                    \
  }                    \
  }

#define FOR_PP_HELICITY               \
  for (const auto &ha : {0, 1}) {     \
    for (const auto &hb : {0, 1}) {   \
      for (const auto &h1 : {0, 1}) { \
        for (const auto &h2 : {0, 1}) {
#define FOR_PP_HELICITY_END \
  }                         \
  }                         \
  }                         \
  }


using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;
using gra::math::zi;
using gra::math::PI;

using FTensor::Tensor0;
using FTensor::Tensor1;
using FTensor::Tensor2;
using FTensor::Tensor3;
using FTensor::Tensor4;

namespace gra {

// 2 -> 3 amplitudes
//
// return value: matrix element squared with helicities summed over
// lts.hamp:     individual helicity amplitudes
//
double MTensorPomeron::ME3(gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const {
  // Kinematics
  const M4Vec pa = lts.pbeam1;
  const M4Vec pb = lts.pbeam2;

  const M4Vec p1 = lts.pfinal[1];
  const M4Vec p2 = lts.pfinal[2];
  const M4Vec p3 = lts.decaytree[0].p4;
  const M4Vec p4 = lts.decaytree[1].p4;

  // ------------------------------------------------------------------

  // Spinors (2 helicities)
  const std::array<std::vector<std::complex<double>>, 2> u_a = SpinorStates(pa, "u");
  const std::array<std::vector<std::complex<double>>, 2> u_b = SpinorStates(pb, "u");

  const std::array<std::vector<std::complex<double>>, 2> ubar_1 = SpinorStates(p1, "ubar");
  const std::array<std::vector<std::complex<double>>, 2> ubar_2 = SpinorStates(p2, "ubar");

  // ------------------------------------------------------------------

  // t-channel Pomeron propagators
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_1 = iD_P(lts.s1, lts.t1);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_2 = iD_P(lts.s2, lts.t2);

  // ------------------------------------------------------------------

  // Resonance parameters
  const double M0    = resonance.p.mass;
  const double Gamma = resonance.p.width;
  const int    J     = resonance.p.spinX2 / 2.0;
  const int    P     = resonance.p.P;

  // ------------------------------------------------------------------

  // Reset helicity amplitude container
  lts.hamp.clear();

  // Free Lorentz indices [second parameter denotes the range of index]
  FTensor::Index<'a', 4> mu1;
  FTensor::Index<'b', 4> nu1;
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;
  FTensor::Index<'k', 4> mu2;
  FTensor::Index<'l', 4> nu2;

  // ====================================================================
  // Construct Central Vertex

  Tensor4<std::complex<double>, 4, 4, 4, 4> cvtx;

  // Pomeron-Pomeron-Scalar coupling
  if (J == 0 && P == 1) {

    cvtx = iG_PPS_total(lts.q1, lts.q2, M0, "scalar", resonance.g_Tensor);

    // Scalar BW-propagator
    std::complex<double> iD = iD_MES(lts.pfinal[0], M0, Gamma);

    // Decay coupling constant
    iD *= resonance.g_decay;

    FOR_EACH_4(LI);
    cvtx(u, v, k, l) *= iD;
    FOR_EACH_4_END;

  // Pomeron-Pomeron-Pseudoscalar coupling
  } else if (J == 0 && P == -1) {

    cvtx = iG_PPS_total(lts.q1, lts.q2, M0, "pseudoscalar", resonance.g_Tensor);

    // Scalar BW-propagator
    std::complex<double> iD = iD_MES(lts.pfinal[0], M0, Gamma);

    // Decay coupling constant
    iD *= resonance.g_decay;

    FOR_EACH_4(LI);
    cvtx(u, v, k, l) *= iD;
    FOR_EACH_4_END;

  // Pomeron-Pomeron-Tensor coupling
  } else if (J == 2) {

    const MTensor<std::complex<double>> iGPPf2 =
        iG_PPT_total(lts.q1, lts.q2, M0, resonance.g_Tensor);

    // f2 -> pipi decay structure
    const Tensor2<std::complex<double>, 4, 4> iGf2PIPI = iG_f2pipi(p3, p4, M0);
    
    // Tensor BW-propagator
    const bool INDEX_UP = true;
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iDf2 = iD_TMES(lts.pfinal[0], M0, Gamma, INDEX_UP);

    // Contract
    Tensor2<std::complex<double>, 4, 4> iD;
    iD(mu1, nu1) = iDf2(mu1, nu1, rho1, rho2) * iGf2PIPI(rho1, rho2);
 
    // Construct central vertex
    FOR_EACH_4(LI);
    cvtx(u, v, k, l) = 0.0;
    for (auto const &rho : LI) {
      for (auto const &sigma : LI) {
        cvtx(u, v, k, l) += iGPPf2({u, v, k, l, rho, sigma}) * iD(rho, sigma);
      }
    }
    FOR_EACH_4_END;
    
  } else {
    throw std::invalid_argument("MTensorPomeron::ME3: Unknown spin-parity input");
  }
  // ====================================================================

  // Two helicity states for incoming and outgoing protons
  FOR_PP_HELICITY;

  // Apply proton leg helicity conservation / No helicity flip
  // (high energy limit)
  // This gives at least 4 x speed improvement
  if (ha != h1 || hb != h2) { continue; }

  // Full proton-Pomeron-proton spinor structure (upper and
  // lower vertex)
  // const Tensor2<std::complex<double>,4,4> iG_1a =
  // iG_Ppp(p1, pa, ubar_1[h1], u_a[ha]);
  // const Tensor2<std::complex<double>,4,4> iG_2b =
  // iG_Ppp(p2, pb, ubar_2[h2], u_b[hb]);

  // High Energy Limit proton-Pomeron-proton spinor structure
  const Tensor2<std::complex<double>, 4, 4> iG_1a = iG_PppHE(p1, pa);
  const Tensor2<std::complex<double>, 4, 4> iG_2b = iG_PppHE(p2, pb);

  // s-channel amplitude
  const std::complex<double> A = (-zi) * iG_1a(mu1, nu1) * iDP_1(mu1, nu1, alpha1, beta1) *
                                 cvtx(alpha1, beta1, alpha2, beta2) * iDP_2(alpha2, beta2, mu2, nu2) *
                                 iG_2b(mu2, nu2);

  lts.hamp.push_back(A);

  FOR_PP_HELICITY_END;

  // Get total amplitude squared 1/4 \sum_h |A_h|^2
  double SumAmp2 = 0.0;
  for (const auto &i : indices(lts.hamp)) {
    SumAmp2 += gra::math::abs2(lts.hamp[i]);
  }
  SumAmp2 /= 4;  // Initial state helicity average

  return SumAmp2; // Amplitude squared
}

// 2 -> 4 amplitudes
//
// return value: matrix element squared with helicities summed over
// lts.hamp:     individual helicity amplitudes
//
double MTensorPomeron::ME4(gra::LORENTZSCALAR &lts) const {

  // Kinematics
  const M4Vec pa = lts.pbeam1;
  const M4Vec pb = lts.pbeam2;

  const M4Vec p1 = lts.pfinal[1];
  const M4Vec p2 = lts.pfinal[2];
  const M4Vec p3 = lts.decaytree[0].p4;
  const M4Vec p4 = lts.decaytree[1].p4;

  // Intermediate boson/fermion mass
  const double M_mes     = lts.decaytree[0].p.mass;
  //const double Gamma_mes = lts.decaytree[0].p.width;
  
  // Momentum convention of sub-diagrams
  //
  // ------<------ anti-particle p3
  // |
  // | \hat{t} (arrow down)
  // |
  // ------>------ particle p4
  //
  // ------<------ particle p4
  // |
  // | \hat{u} (arrow up)
  // |
  // ------>------ anti-particle p3
  //
  const M4Vec pt = pa - p1 - p3;  // => q1 = pt + p3, q2 = pt - p4
  const M4Vec pu = p4 - pa + p1;  // => q2 = pu + p3, q1 = pu - p3

  // ------------------------------------------------------------------

  // Spinors (2 helicities)
  const std::array<std::vector<std::complex<double>>, 2> u_a = SpinorStates(pa, "u");
  const std::array<std::vector<std::complex<double>>, 2> u_b = SpinorStates(pb, "u");

  const std::array<std::vector<std::complex<double>>, 2> ubar_1 = SpinorStates(p1, "ubar");
  const std::array<std::vector<std::complex<double>>, 2> ubar_2 = SpinorStates(p2, "ubar");

  // Massive polarization vectors (3 helicities)
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_3_conj = Spin1States(p3, "conj");
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_4_conj = Spin1States(p4, "conj");

  // ------------------------------------------------------------------

  // t-channel Pomeron propagators
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_13 = iD_P(lts.ss[1][3], lts.t1);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_24 = iD_P(lts.ss[2][4], lts.t2);

  // u-channel Pomeron propagators
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_14 = iD_P(lts.ss[1][4], lts.t1);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> iDP_23 = iD_P(lts.ss[2][3], lts.t2);

  // ------------------------------------------------------------------

  // Fermion propagator
  const double                        m     = lts.decaytree[0].p.mass;
  const MMatrix<std::complex<double>> iSF_t = iD_F(pt, m);
  const MMatrix<std::complex<double>> iSF_u = iD_F(pu, m);

  // Central spinors (2 helicities)
  const std::array<std::vector<std::complex<double>>, 2> v_3    = SpinorStates(p3, "v");
  const std::array<std::vector<std::complex<double>>, 2> ubar_4 = SpinorStates(p4, "ubar");

  // Reset
  lts.hamp.clear();

  // Free Lorentz indices [second parameter denotes the range of index]
  FTensor::Index<'a', 4> mu1;
  FTensor::Index<'b', 4> nu1;
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;
  FTensor::Index<'k', 4> mu2;
  FTensor::Index<'l', 4> nu2;

  // DEDUCE mode
  std::string SPINMODE;
  if        (lts.decaytree[0].p.spinX2 == 0 && lts.decaytree[1].p.spinX2 == 0) {
    SPINMODE = "2xS";
  } else if (lts.decaytree[0].p.spinX2 == 1 && lts.decaytree[1].p.spinX2 == 1) {
    SPINMODE = "2xF";
  } else if (lts.decaytree[0].p.spinX2 == 2 && lts.decaytree[1].p.spinX2 == 2) {
    SPINMODE = "2xV";
  } else {
    throw std::invalid_argument(
        "MTensorPomeron::ME4: Invalid daughter spin "
        "(J = 0,1/2,1 pairs supported)");
  }

  // Two helicity states for incoming and outgoing protons
  FOR_PP_HELICITY;
  // Apply proton leg helicity conservation / No helicity flip
  // (high energy limit)
  // This gives at least 4 x speed improvement
  if (ha != h1 || hb != h2) { continue; }

  // Full proton-Pomeron-proton spinor structure (upper and
  // lower vertex)
  // const Tensor2<std::complex<double>,4,4> iG_1a =
  // iG_Ppp(p1, pa, ubar_1[h1], u_a[ha]);
  // const Tensor2<std::complex<double>,4,4> iG_2b =
  // iG_Ppp(p2, pb, ubar_2[h2], u_b[hb]);

  // High Energy Limit proton-Pomeron-proton spinor structure
  const Tensor2<std::complex<double>, 4, 4> iG_1a = iG_PppHE(p1, pa);
  const Tensor2<std::complex<double>, 4, 4> iG_2b = iG_PppHE(p2, pb);

  // ==============================================================
  // 2 x pseudoscalar
  if (SPINMODE == "2xS") {

    // t-channel blocks
    const Tensor2<std::complex<double>, 4, 4> iG_ta = iG_Ppsps(pt, -p3);
    const std::complex<double> iDMES_t = iD_MES0(pt, M_mes);
    const Tensor2<std::complex<double>, 4, 4> iG_tb = iG_Ppsps(p4, pt);

    // u-channel blocks
    const Tensor2<std::complex<double>, 4, 4> iG_ua = iG_Ppsps(p4, pu);
    const std::complex<double> iDMES_u = iD_MES0(pu, M_mes);
    const Tensor2<std::complex<double>, 4, 4> iG_ub = iG_Ppsps(pu, -p3);

    std::complex<double> M_t;
    std::complex<double> M_u;
    // t-channel
    {
      // Upper block
      const std::complex<double> A = iG_1a(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) *
                                     iG_ta(alpha1, beta1) * iG_tb(alpha2, beta2) *
                                     iDP_24(alpha2, beta2, mu2, nu2) * iG_2b(mu2, nu2);
      // Apply off-shell form factors (two of them)
      M_t = A * iDMES_t * pow2(PARAM_REGGE::Meson_FF(pt.M2(), pow2(M_mes)));
    }

    // u-channel
    {
      // Upper block
      const std::complex<double> A = iG_1a(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) *
                                     iG_ua(alpha1, beta1) * iG_ub(alpha2, beta2) *
                                     iDP_23(alpha2, beta2, mu2, nu2) * iG_2b(mu2, nu2);
      // Apply off-shell form factors (two of them)
      M_u = A * iDMES_u * pow2(PARAM_REGGE::Meson_FF(pu.M2(), pow2(M_mes)));
    }

    // Full amplitude: iM = [ ... ]  <-> M = (-i)*[ ...
    // ]
    std::complex<double> M = (-zi) * (M_t + M_u);
    lts.hamp.push_back(M);
  }

  // ==============================================================
  // 2 x fermion
  if (SPINMODE == "2xF") {

    for (const auto &h3 : {0, 1}) {
      for (const auto &h4 : {0, 1}) {
        // Fermion propagator and connected
        // parts
        const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_t =
            iG_PppbarP(p4, ubar_4[h4], pt, iSF_t, v_3[h3], -p3);
        const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_u =
            iG_PppbarP(p4, ubar_4[h4], pu, iSF_u, v_3[h3], -p3);

        std::complex<double> M_t;
        std::complex<double> M_u;
        // t-channel
        {
          // Upper block
          const std::complex<double> A = iG_1a(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) *
                                         iG_t(alpha2, beta2, alpha1, beta1) *
                                         iDP_24(alpha2, beta2, mu2, nu2) * iG_2b(mu2, nu2);
          // Apply off-shell form
          // factors (two of them)
          M_t = A * pow2(PARAM_REGGE::Baryon_FF(pt.M2(), pow2(M_mes)));
        }

        // u-channel
        {
          // Upper block
          const std::complex<double> A = iG_1a(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) *
                                         iG_u(alpha1, beta1, alpha2, beta2) *
                                         iDP_23(alpha2, beta2, mu2, nu2) * iG_2b(mu2, nu2);
          // Apply off-shell form
          // factors (two of them)
          M_u = A * pow2(PARAM_REGGE::Baryon_FF(pu.M2(), pow2(M_mes)));
        }

        // Full amplitude: iM = [ ... ]  <->
        // M = (-i)*[ ... ]
        std::complex<double> M = (-zi) * (M_t + M_u);
        lts.hamp.push_back(M);
      }
    }
  }
  
  // ==============================================================
  // 2 x vector meson
  if (SPINMODE == "2xV") {

    // t-channel blocks
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_tA = iG_Pvv(pt, -p3);
    const Tensor2<std::complex<double>, 4, 4> iDMES_t = iD_V(pt, M_mes, lts.pfinal[0].M2());
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_tB = iG_Pvv(p4, pt);

    // u-channel blocks
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_uA = iG_Pvv(p4, pu);
    const Tensor2<std::complex<double>, 4, 4> iDMES_u = iD_V(pu, M_mes, lts.pfinal[0].M2());
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iG_uB = iG_Pvv(pu, -p3);

    FTensor::Index<'e', 4> rho3;
    FTensor::Index<'f', 4> rho4;

    // -------------------------------------------------------------------------
    // TENSOR LORENTZ INDEX CONTRACTION BLOCK
    // This is done in pieces, due to template <>
    // argument deduction constraints

    Tensor2<std::complex<double>, 4, 4> M_t;
    Tensor2<std::complex<double>, 4, 4> M_u;

    // t-channel
    {
      // Upper block
      Tensor2<std::complex<double>, 4, 4> A;
      A(rho3, rho1) =
          iG_1a(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) * iG_tA(rho3, rho1, alpha1, beta1);
      // Lower block
      Tensor2<std::complex<double>, 4, 4> B;
      B(rho2, rho4) =
          iG_2b(mu2, nu2) * iDP_24(alpha2, beta2, mu2, nu2) * iG_tB(rho2, rho4, alpha2, beta2);
      // Apply off-shell form factors (two of them)
      M_t(rho3, rho4) = A(rho3, rho1) * iDMES_t(rho1, rho2) * B(rho2, rho4) *
                        pow2(PARAM_REGGE::Meson_FF(pt.M2(), pow2(M_mes)));
    }

    // u-channel
    {
      // Upper block
      Tensor2<std::complex<double>, 4, 4> A;
      A(rho4, rho1) =
          iG_1a(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) * iG_uA(rho4, rho1, alpha1, beta1);
      // Lower block
      Tensor2<std::complex<double>, 4, 4> B;
      B(rho2, rho3) =
          iG_2b(mu2, nu2) * iDP_23(alpha2, beta2, mu2, nu2) * iG_uB(rho2, rho3, alpha2, beta2);
      // Apply off-shell form factors (two of them)
      M_u(rho3, rho4) = A(rho4, rho1) * iDMES_u(rho1, rho2) * B(rho2, rho3) *
                        pow2(PARAM_REGGE::Meson_FF(pu.M2(), pow2(M_mes)));
    }

    // Total amplitude: iM = [ ... ]  <-> M = (-i)*[ ...
    // ]
    Tensor2<std::complex<double>, 4, 4> M;
    FOR_EACH_2(LI);
      M(u,v) = (-zi) * (M_t(u,v) + M_u(u,v));
    FOR_EACH_2_END;

    const int OPTION = 1; // Use 1 for now

    // Polarization sum
    if (OPTION == 1) {
      double Amp2 = 0.0;

      FOR_EACH_4(LI);

        const std::complex<double> contract = std::conj(M(u, v)) * M(k, l) * g[u][k] * g[v][l];
        Amp2 += std::real(contract);  // real for casting to double

      FOR_EACH_4_END;

      lts.hamp.push_back(msqrt(Amp2)); // Phase information lost here
    }
    // --------------------------------------------------
    
    // Explicit polarization vectors [NEEDS TESTING]
    if (OPTION == 2) {
      // Loop over massive Spin-1 helicity states

      for (const auto &h3 : indices(eps_3_conj)) {
        for (const auto &h4 : indices(eps_4_conj)) {
          // Contract Lorentz indices
          // to get the full helicity
          // amplitude A(ha,hb,h1,h2,h3,h4)
          const std::complex<double> amp =
              eps_3_conj[h3](rho3) * eps_4_conj[h4](rho4) * M(rho3, rho4);

          lts.hamp.push_back(amp); // Phase information retained
        }
      }
    }
  }
  FOR_PP_HELICITY_END;

  // Get total amplitude squared 1/4 \sum_h |A_h|^2
  double SumAmp2 = 0.0;
  for (const auto &i : indices(lts.hamp)) {
    SumAmp2 += math::abs2(lts.hamp[i]);
  }
  SumAmp2 /= 4;  // Initial state helicity average

  return SumAmp2;  // Amplitude squared
}


// -------------------------------------------------------------------------
// Vertex functions

// Gamma-Lepton-Lepton vertex
// contracted in \bar{spinor} G_\mu \bar{spinor}
//
// iGamma_\mu (p',p)
//
// High Energy Limit:
// ubar(p',\lambda') \gamma_\mu u(p,\lambda ~= (p' + p)_\mu
// \delta_{\lambda',\lambda}
//
// Input as contravariant (upper index) 4-vectors
//
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_yee(
    const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
    const std::vector<std::complex<double>> &u) const {
  // const double q2 = (prime-p).M2();
  const double e = msqrt(alpha_QED * 4.0 * PI);  // ~ 0.3

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) {
    // \bar{spinor} [Gamma Matrix] \spinor product
    T(mu) = zi * e * gra::matoper::VecVecMultiply(ubar, gamma_lo[mu] * u);
  }
  return T;
}

// Gamma-Proton-Proton vertex: iGamma_\mu (p', p)
//
//
// contracted in \bar{spinor} iG_\mu \bar{spinor}
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_ypp(
    const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
    const std::vector<std::complex<double>> &u) const {
  const double t    = (prime - p).M2();
  const double e    = msqrt(alpha_QED * 4.0 * PI);  // ~ 0.3
  const M4Vec  psum = prime - p;

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) {
    MMatrix<std::complex<double>> SUM(4, 4, 0.0);
    for (const auto &                   nu : LI) { SUM += sigma_lo[mu][nu] * psum[nu]; }
    const MMatrix<std::complex<double>> A = gamma_lo[mu] * F1(t);
    const MMatrix<std::complex<double>> B = SUM * zi / (2 * mp) * F2(t);
    const MMatrix<std::complex<double>> G = (A + B) * ((-1.0) * zi * e);

    // \bar{spinor} [Gamma Matrix] \spinor product
    T(mu) = gra::matoper::VecVecMultiply(ubar, G * u);
  }
  return T;
}

// High-Energy limit proton-Pomeron-proton vertex times \delta^{lambda_prime,
// \lambda} (helicity conservation)
// (~1.5 faster evaluation than the exact spinor structure)
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_PppHE(const M4Vec &prime,
                                                             const M4Vec p) const {
  const M4Vec psum = prime + p;
  const std::complex<double> FACTOR = -zi * 3.0 * betaPNN * F1((prime - p).M2());

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u,v) = FACTOR * (psum % u) * (psum % v);
  FOR_EACH_2_END;

  return T;
}

// Pomeron-Proton-Proton [-Neutron-Neutron, -Antiproton-Antiproton] vertex
// function:
// contracted in \bar{spinor} G_{\mu\nu} \bar{spinor}
//
// i\Gamma_{\mu\nu} (p', p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_Ppp(
    const M4Vec &prime, const M4Vec &p, const std::vector<std::complex<double>> &ubar,
    const std::vector<std::complex<double>> &u) const {
  const double t    = (prime - p).M2();
  const M4Vec  psum = prime + p;

  Tensor2<std::complex<double>, 4, 4> T;
  const std::complex<double> FACTOR = -zi * 3.0 * betaPNN * F1(t);

  // Lorentz indices
  MMatrix<std::complex<double>> A;
  MMatrix<std::complex<double>> MAT;

  const MMatrix<std::complex<double>> slash = FSlash(psum);

  for (const auto &mu : LI) {
    for (const auto &nu : LI) {
      A = (gamma_lo[mu] * (psum % nu) + gamma_lo[nu] * (psum % mu)) * 0.5;
      if (mu == nu) {  // Speed it up
        MAT = (A - slash * 0.25 * g[mu][nu]) * FACTOR;
      } else {
        MAT = A * FACTOR;
      }

      // \bar{spinor} [Gamma Matrix] \spinor product
      T(mu, nu) = gra::matoper::VecVecMultiply(ubar, MAT * u);
    }
  }
  return T;
}

// Double vertex function in High Energy Limit
// Computationally very heavy, should optimize this [TBD]
//
// Pomeron - outgoing proton spinor ubar -
//           <antiproton/proton fermion propagator>
//         - outgoing antiproton spinor v - Pomeron
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PppbarP(
    const M4Vec &prime, const std::vector<std::complex<double>> &ubar, const M4Vec &pt,
    const MMatrix<std::complex<double>> &iSF, const std::vector<std::complex<double>> &v,
    const M4Vec &p) const {
  const std::complex<double> lhs_FACTOR = -zi * 3.0 * betaPNN * F1((prime - pt).M2());
  const std::complex<double> rhs_FACTOR = -zi * 3.0 * betaPNN * F1((pt - p).M2());

  const M4Vec lhs_psum = prime + pt;
  const M4Vec rhs_psum = pt + p;
  
  // Feynman slashes
  const MMatrix<std::complex<double>> lhs_slash = FSlash(lhs_psum);
  const MMatrix<std::complex<double>> rhs_slash = FSlash(rhs_psum);

  // Aux matrices
  MMatrix<std::complex<double>> lhs_A;
  MMatrix<std::complex<double>> lhs_MAT;

  MMatrix<std::complex<double>> rhs_A;
  MMatrix<std::complex<double>> rhs_MAT;

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  for (const auto &mu2 : LI) {
    for (const auto &nu2 : LI) {
      lhs_A = (gamma_lo[mu2] * (lhs_psum % nu2) + gamma_lo[nu2] * (lhs_psum % mu2)) * 0.5;
      if (mu2 == nu2) {  // Speed it up
        lhs_MAT = (lhs_A - lhs_slash * 0.25 * g[mu2][nu2]) * lhs_FACTOR;
      } else {
        lhs_MAT = lhs_A * lhs_FACTOR;
      }

      // Fermion propagator applied in the middle
      const std::vector<std::complex<double>> ubarM =
          gra::matoper::VecMatMultiply(ubar, lhs_MAT * iSF);

      for (const auto &mu1 : LI) {
        for (const auto &nu1 : LI) {
          rhs_A = (gamma_lo[mu1] * (rhs_psum % nu1) + gamma_lo[nu1] * (rhs_psum % mu1)) * 0.5;
          if (mu1 == nu1) {  // Speed it up
            rhs_MAT = (rhs_A - rhs_slash * 0.25 * g[mu1][nu1]) * rhs_FACTOR;
          } else {
            rhs_MAT = rhs_A * rhs_FACTOR;
          }

          // \bar{spinor} [Gamma Matrix] \spinor product
          T(mu2, nu2, mu1, nu1) = gra::matoper::VecVecMultiply(ubarM, rhs_MAT * v);
        }
      }
    }
  }
  return T;
}

// ======================================================================

// Pomeron (\mu\nu) - Pomeron (\kappa\lambda) - Scalar/Pseudoscalar
// vertex function: i\Gamma_{\mu\nu,\kappa\lambda,\rho\sigma}
//
// Input as contravariant (upper index) 4-vector, M0 the resonance peak mass
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_total(
    const M4Vec &q1, const M4Vec &q2, double M0, const std::string &mode,
    const std::vector<double> &g_PPS) const {
  auto FM_ = [](double t) -> double {
    const double LAMBDA = 1.0;  // GeV
    return 1.0 / (1.0 - t / pow2(LAMBDA));
  };
  auto F_ = [](double q2, double M0) -> double {
    const double LAMBDA = 1.0;  // GeV^2
    return std::exp(-pow2(q2 - pow2(M0)) / math::pow4(LAMBDA));
  };

  // Form FACTOR
  const double F_tilde = FM_(q1.M2()) * FM_(q2.M2()) * F_((q1 + q2).M2(), M0);

  // Tensor structures
  Tensor4<std::complex<double>, 4, 4, 4, 4> iG1;
  Tensor4<std::complex<double>, 4, 4, 4, 4> iG2;

  if        (mode == "scalar") {
    iG1 = MTensorPomeron::iG_PPS_0(g_PPS[0]);
    iG2 = MTensorPomeron::iG_PPS_1(q1, q2, g_PPS[1]);
  } else if (mode == "pseudoscalar") {
    iG1 = MTensorPomeron::iG_PPPS_0(q1, q2, g_PPS[0]);
    iG2 = MTensorPomeron::iG_PPPS_1(q1, q2, g_PPS[1]);
  } else {
    throw std::invalid_argument(
        "MTensorPomeron::iG_PPS_total: Unknown mode (should be scalar or "
        "pseudoscalar)");
  }

  // Construct vertex
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = (iG1(u, v, k, l) + iG2(u, v, k, l)) * F_tilde;
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Scalar coupling structure #0
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (0,0)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_0(double g_PPS) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = zi * g_PPS * S0;

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Scalar coupling structure #1
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (2,0)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_1(const M4Vec &q1, const M4Vec &q2,
                                                                   double g_PPS) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const double               q1q2   = q1 * q2;
  const std::complex<double> FACTOR = zi * g_PPS / (2 * S0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * ((q1 % (k)) * (q2 % (u)) * g[v][l] + (q1 % (k)) * (q2 % (v)) * g[u][l] +
                            (q1 % (l)) * (q2 % (u)) * g[v][k] + (q1 % (l)) * (q2 % (v)) * g[u][k] -
                            2.0 * q1q2 * (g[u][k] * g[v][l] + g[v][k] * g[u][l]));
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Pseudoscalar coupling structure #0
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (1,1)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPPS_0(const M4Vec &q1,
                                                                    const M4Vec &q2,
                                                                    double g_PPPS) const {
  const double S0 = 1.0;  // Mass scale (GeV)

  const M4Vec q1_q2 = q1 - q2;
  const M4Vec q     = q1 + q2;

  const std::complex<double> FACTOR = zi * g_PPPS / (2.0 * S0);
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;

  // Contract r and s indices
  for (const auto &r : LI) {
    for (const auto &s : LI) {
      // Note +=
      T(u, v, k, l) += FACTOR * (g[u][k] * eps_lo(v, l, r, s) + g[v][k] * eps_lo(u, l, r, s) +
                                 g[u][l] * eps_lo(v, k, r, s) + g[v][l] * eps_lo(u, k, r, s)) *
                       q1_q2[r] * q[s];
    }
  }
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Pseudoscalar coupling structure #1
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (3,3)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPPS_1(const M4Vec &q1,
                                                                    const M4Vec &q2,
                                                                    double       g_PPPS) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const M4Vec                q1_q2  = q1 - q2;
  const M4Vec                q      = q1 + q2;
  const double               q1q2   = q1 * q2;
  const std::complex<double> FACTOR = zi * g_PPPS / gra::math::pow3(S0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;

  // Contract r and s indices
  for (const auto &r : LI) {
    for (const auto &s : LI) {
      // Note +=
      T(u, v, k, l) += FACTOR * (eps_lo(v, l, r, s) * ((q1 % k) * (q2 % u) - q1q2 * g[u][k]) +
                                 eps_lo(u, l, r, s) * ((q1 % k) * (q2 % v) - q1q2 * g[v][k]) +
                                 eps_lo(v, k, r, s) * ((q1 % l) * (q2 % u) - q1q2 * g[u][l]) +
                                 eps_lo(u, k, r, s) * ((q1 % l) * (q2 % v) - q1q2 * g[v][l])) *
                       q1_q2[r] * q[s];
    }
  }
  FOR_EACH_4_END;

  return T;
}

// ======================================================================

// Pomeron (\mu\nu) - Pomeron (\kappa\lambda) - Tensor (\rho\sigma)
// vertex function: i\Gamma_{\mu\nu,\kappa\lambda,\rho\sigma}
//
// Input as contravariant (upper index) 4-vector, M0 the resonance peak mass
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_total(const M4Vec &q1, const M4Vec &q2,
                                                           double                     M0,
                                                           const std::vector<double> &g_PPT) const {
  auto FM_ = [](double t) -> double {
    const double LAMBDA = 1.0;  // GeV
    return 1.0 / (1.0 - t / pow2(LAMBDA));
  };
  auto F_ = [](double q2, double M0) -> double {
    const double LAMBDA = 1.0;  // GeV^2
    return std::exp(-pow2(q2 - pow2(M0)) / math::pow4(LAMBDA));
  };

  // Form FACTOR
  const double F_tilde = FM_(q1.M2()) * FM_(q2.M2()) * F_((q1 + q2).M2(), M0);

  // Kinematics independent tensor structure [0]
  static const MTensor<std::complex<double>> iG0 = iG_PPT_00();

  // Kinematics dependent tensor structures  [1 ... 6]
  std::vector<MTensor<std::complex<double>>> iGx;

  const double EPSILON = 1e-6;
  if (g_PPT[1] > EPSILON)
    iGx.push_back( iG_PPT_12(q1, q2, g_PPT[1], 1) );

  if (g_PPT[2] > EPSILON)
    iGx.push_back( iG_PPT_12(q1, q2, g_PPT[2], 2) );

  if (g_PPT[3] > EPSILON)
    iGx.push_back( iG_PPT_03(q1, q2, g_PPT[3]) );

  if (g_PPT[4] > EPSILON)
    iGx.push_back( iG_PPT_04(q1, q2, g_PPT[4]) );

  if (g_PPT[5] > EPSILON)
    iGx.push_back( iG_PPT_05(q1, q2, g_PPT[5]) );

  if (g_PPT[6] > EPSILON)
    iGx.push_back( iG_PPT_06(q1, q2, g_PPT[6]) );

  // Construct total vertex by summing up the tensors
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));

  FOR_EACH_6(LI);

  // Indexing
  const std::vector<size_t> ind = {u, v, k, l, r, s};

  // For number [0], we apply coupling here
  T(ind)                        = iG0(ind)*g_PPT[0];

  // For number [1-6], couplings have been already applied, just loop over
  for (const auto& i : indices(iGx)) {
    T(ind) += iGx[i](ind);
  }

  // Apply form factors
  T(ind) *= F_tilde;

  FOR_EACH_6_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #0
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (0,2)
//
// This is kinematics independent, can be pre-calculated.
//
// Remember to apply g_PPT (coupling constant) outside this
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_00() const {
  const double                  M0     = 1.0;
  const std::complex<double>    FACTOR = (2.0 * zi) * M0;

  // Init with zeros!
  MTensor<std::complex<double>> T      = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));

  FOR_EACH_6(LI);
  const std::vector<size_t> ind = {u, v, k, l, r, s};

  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &r1 : LI) {
          for (auto const &s1 : LI) {
            for (auto const &u1 : LI) {
              if (v1 != a1 || l1 != r1 || s1 != u1) { continue; }  // speed it up

              // Note +=
              T(ind) += FACTOR * R_DDDD(u, v, u1, v1) * R_DDDD(k, l, a1, l1) *
                                 R_DDDD(r, s, r1, s1) * g[v1][a1] * g[l1][r1] * g[s1][u1];
            }
          }
        }
      }
    }
  }
  FOR_EACH_6_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structures #1 and #2
// iG_{\mu \nu \kappa \lambda \rho \sigma}
// ~ (l,s) = (2,0) - (2,2)
// ~ (l,s) = (2,0) + (2,2)
// 
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_12(const M4Vec &q1, const M4Vec &q2,
                                                        double g_PPT, int mode) const {
  if (!(mode == 1 || mode == 2)) {
    throw std::invalid_argument("MTensorPomeron::iG_PPT_12: Error, mode should be 1 or 2");
  }
  const double                  M0     = 1.0;
  const double                  q1q2   = q1 * q2;
  const std::complex<double>    FACTOR = -(2.0 * zi) / M0 * g_PPT;

  // Coupling structure 1 or 2
  const double sign = (mode == 1) ? -1.0 : 1.0;

  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  const std::vector<size_t> ind = {u, v, k, l, r, s};

  // Pre-calculated tensor contraction structures T1,T2,T2 used here
  T(ind) += q1q2 * T1(ind);

  for (auto const &u1 : LI) {
    for (auto const &r1 : LI) {
      // Note +=
      T(ind) += sign * q2[u1] * (q1 % r1) * T2({u,v,k,l,r,s,u1,r1});
    }
  }
  for (auto const &u1 : LI) {
    for (auto const &s1 : LI) {
      // Note +=
      T(ind) += sign * q1[u1] * (q2 % s1) * T3({u,v,k,l,r,s,u1,s1});
    }
  }
  for (auto const &r1 : LI) {
    for (auto const &s1 : LI) {
      // Note +=
      T(ind) += (q1 % r1) * (q2 % s1) * R_DDDD(u, v, k, l) * R_DDUU(r, s, r1, s1);
    }
  }

  T(ind) *= FACTOR;
  FOR_EACH_6_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #3
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (2,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_03(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {

  const double                  M0     = 1.0;
  const std::complex<double>    FACTOR = - zi / M0 * g_PPT;

  Tensor3<double, 4,4,4> A;
  Tensor3<double, 4,4,4> B;
  Tensor3<double, 4,4,4> C;
  Tensor3<double, 4,4,4> D;  

  // Manual contractions
  FOR_EACH_3(LI);
  A(u,v,k) = 0.0;
  B(u,v,k) = 0.0;
  for (auto const &x : LI) {
    // Note +=
    A(u,v,k) += q2[x] * R_DDDD(u,v,x,k);
    B(u,v,k) += q1[x] * R_DDDD(u,v,x,k);
  }
  FOR_EACH_3_END;

  // Final rank-6 tensor
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  const std::vector<size_t> ind = {u, v, k, l, r, s};
  for (auto const &v1 : LI) {
    for (auto const &l1 : LI) {
      // Note +=
      T(ind) += (A(u,v,v1) * B(k,l,l1) + A(u,v,l1) * B(k,l,v1)) * R_UUDD(v1,l1,r,s);
    }
  }
  T(ind) *= FACTOR;
  /*

  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
          for (auto const &u1 : LI) {

            // Note +=
            T(ind) += FACTOR * (q2[u1] * R_DDDD(u,v,u1,v1) * q1[a1] * R_DDDD(k,l,a1,l1) +
                                q2[a1] * R_DDDD(u,v,a1,l1) * q1[u1] * R_DDDD(k,l,u1,v1)) * R_UUDD(v1,l1,r,s);
        }
      }
    }
  }
  */
  FOR_EACH_6_END;
  
  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #4
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (4,2)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_04(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  
  const double                  M0     = 1.0;
  const std::complex<double>    FACTOR = - 2.0 * zi / math::pow3(M0) * g_PPT;
  const double q1q2 = q1*q2;

  Tensor1<double, 4> q1_D = {q1 % 0, q1 % 1, q1 % 2, q1 % 3};
  Tensor1<double, 4> q2_D = {q2 % 0, q2 % 1, q2 % 2, q2 % 3};

  Tensor4<double, 4,4,4,4> A;
  Tensor4<double, 4,4,4,4> B;
  Tensor2<double, 4,4> C;

  // Manual contractions (not possible with autocontraction)
  FOR_EACH_4(LI);
  A(u,v,k,l) = 0.0;
  B(u,v,k,l) = 0.0;

  for (auto const &u1 : LI) {
    for (auto const &v1 : LI) {
      for (auto const &a1 : LI) {
        // Note +=
        A(u,v,k,l) += q2[v1] * R_DDDD(u,v,v1,a1) * q1[u1] * R_DDDU(k,l,u1,a1);
        B(u,v,k,l) += q2[u1] * R_DDDD(u,v,u1,a1) * q1[v1] * R_DDDU(k,l,v1,a1);
      }
    }
  }
  FOR_EACH_4_END;

  // Autocontractions
  {
    FTensor::Index<'a', 4> r;
    FTensor::Index<'b', 4> s;
    FTensor::Index<'c', 4> a1;
    FTensor::Index<'d', 4> l1;

    C(r,s) = q1_D(a1) * q2_D(l1) * R_UUDD(a1,l1,r,s);
  }

  // Final rank-6 expression
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  T({u,v,k,l,r,s}) = FACTOR * ( A(u,v,k,l) + B(u,v,k,l) - 2.0 * q1q2 * R_DDDD(u,v,k,l)) * C(r,s);

  /*
  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &u1 : LI) {

          // Note +=
          T(ind) += FACTOR * (q2[v1] * R_DDDD(u,v,v1,a1) * q1[u1] * R_DDDU(k,l,u1,a1)
                            + q2[u1] * R_DDDD(u,v,u1,a1) * q1[v1] * R_DDDU(k,l,v1,a1)
                            - 2.0 * q1q2 * R_DDDD(u,v,k,l)) * (q1 % a1) * (q2 % l1) * R_UUDD(a1,l1,r,s);
        }
      }
    }
  }
  */
  FOR_EACH_6_END;
  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #05
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (4,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_05(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  
  const double                  M0     = 1.0;
  const std::complex<double>    FACTOR = 2.0 * zi / math::pow3(M0) * g_PPT;

  Tensor1<double, 4> q1_U = {q1[0], q1[1], q1[2], q1[3]};
  Tensor1<double, 4> q2_U = {q2[0], q2[1], q2[2], q2[3]};
  Tensor1<double, 4> q1_D = {q1 % 0, q1 % 1, q1 % 2, q1 % 3};
  Tensor1<double, 4> q2_D = {q2 % 0, q2 % 1, q2 % 2, q2 % 3};

  Tensor4<double, 4,4,4,4> A;
  Tensor2<double, 4,4> B;
  Tensor4<double, 4,4,4,4> C;
  Tensor2<double, 4,4> D;

  // Manual contractions (not possible with autocontraction)
  MTensor<std::complex<double>> T      = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  A(u,v,r,s) = 0.0;
  C(k,l,r,s) = 0.0;

  for (auto const &u1 : LI) {
    for (auto const &v1 : LI) {
      for (auto const &r1 : LI) {
        // Note +=
        A(u,v,r,s) += q2_U(u1) * R_DDDD(u,v,u1,v1) * q2_D(r1) * R_UUDD(v1,r1,r,s);
        C(k,l,r,s) += q1_U(u1) * R_DDDD(k,l,u1,v1) * q1_D(r1) * R_UUDD(v1,r1,r,s);
      }
    }
  }
  FOR_EACH_6_END;

  // Autocontractions
  {
    FTensor::Index<'a', 4> u;
    FTensor::Index<'b', 4> v;
    FTensor::Index<'c', 4> k;
    FTensor::Index<'d', 4> l;
    FTensor::Index<'j', 4> a1;
    FTensor::Index<'k', 4> l1;

    B(k,l) = q1_U(a1) * q1_U(l1) * R_DDDD(k,l,a1,l1);
    D(u,v) = q2_U(a1) * q2_U(l1) * R_DDDD(u,v,a1,l1);
  }
  
  // Final rank-6 expression
  FOR_EACH_6(LI);
  T({u, v, k, l, r, s}) = FACTOR * (A(u,v,r,s)*B(k,l) + C(k,l,r,s)*D(u,v));

  /*
  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &u1 : LI) {
          for (auto const &r1 : LI) {

            // Note +=
            T(ind) += FACTOR * (q2[u1] * (q2 % r1) * R_DDDD(u,v,u1,v1) * R_UUDD(v1,r1,r,s) * q1[a1] * q1[l1] * R_DDDD(k,l,a1,l1) 
                             +  q1[u1] * (q1 % r1) * R_DDDD(k,l,u1,v1) * R_UUDD(v1,r1,r,s) * q2[a1] * q2[l1] * R_DDDD(u,v,a1,l1));
          }
        }
      }
    }
  }
  */
  FOR_EACH_6_END;
  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #6
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (6,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_06(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  
  const double               M0     = 1.0;
  const std::complex<double> FACTOR = - 2.0 * zi / math::pow5(M0) * g_PPT;

  FTensor::Index<'a', 4> a;
  FTensor::Index<'b', 4> b;
  FTensor::Index<'c', 4> c;
  FTensor::Index<'d', 4> d;

  Tensor1<double, 4> q1_U = {q1[0], q1[1], q1[2], q1[3]};
  Tensor1<double, 4> q2_U = {q2[0], q2[1], q2[2], q2[3]};

  Tensor2<double, 4,4> A;
  Tensor2<double, 4,4> B;
  Tensor2<double, 4,4> C;
  
  // Autocontractions
  A(a,b) = q2_U(c) * q2_U(d) * R_DDDD(a,b,c,d);
  B(a,b) = q1_U(c) * q1_U(d) * R_DDDD(a,b,c,d);
  C(a,b) = q1_U(c) * q2_U(d) * R_DDDD(a,b,c,d);

  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  T({u, v, k, l, r, s}) = FACTOR * A(u,v) * B(k,l) * C(r,s);
  FOR_EACH_6_END;

  return T;
}


// ======================================================================

// Rho-Pi-Pi vertex function: i\Gamma_\mu(k1,k2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_rhopipi(const M4Vec &k1,
                                                            const M4Vec &k2) const {
  const double               g_rhopipi = 11.51;
  const M4Vec                p         = k1 - k2;
  const std::complex<double> FACTOR    = -0.5 * zi * g_rhopipi;

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) { T(mu) = FACTOR * (p % mu); }
  return T;
}

// f2-Pi-Pi vertex function: i\Gamma_{\kappa\lambda}(k1,k2)
// 
// 
// Input as contravariant (upper index) 4-vectors, M0 is the f2 meson on-shell mass
// 
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_f2pipi(const M4Vec &k1,
                                                              const M4Vec &k2, double M0) const {
  // Form FACTOR
  const double lambda2  = 1.0; // GeV^2 (normalization scale)
  auto F_f2pipi                       = [&](double q2) { return std::exp(-(q2 - M0*M0) / lambda2); };

  const double               g_f2pipi = 9.26;  // Coupling
  const double               S0       = 1.0;   // Mass scale (GeV)
  const M4Vec                k1_k2    = k1 - k2;
  const std::complex<double> FACTOR   = -zi * g_f2pipi / (2 * S0) * F_f2pipi((k1+k2).M2());

  const double k1_k2_sq = k1_k2.M2();

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u,v) = FACTOR * ((k1_k2 % u) * (k1_k2 % v) - 0.25 * g[u][v] * k1_k2_sq);
  FOR_EACH_2_END;
  
  return T;
}

// f2-gamma-gamma vertex function: i\Gamma_{\mu\nu\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors, M0 is the f2 meson on-shell mass
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_f2yy(const M4Vec &k1,
                                                                  const M4Vec &k2, double M0) const {
  const double e = msqrt(alpha_QED * 4.0 * PI);  // ~ 0.3

  // Form FACTOR
  const double lambda2  = 1.0; // GeV^2 (normalization scale)
  auto F_f2yy = [&](double q2) { return std::exp(-(q2 - M0*M0) / lambda2); };

  // Couplings
  const double a_f2yy = pow2(e) / (4 * PI) * 1.45;  // GeV^{-3}
  const double b_f2yy = pow2(e) / (4 * PI) * 2.49;  // GeV^{-1}

  const Tensor4<std::complex<double>, 4, 4, 4, 4> G0 = Gamma0(k1, k2);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> G2 = Gamma2(k1, k2);

  const std::complex<double> FACTOR = zi * FM(k1.M2()) * FM(k2.M2()) * F_f2yy((k1 + k2).M2());

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (2.0 * a_f2yy * G0(u, v, k, l) - b_f2yy * G2(u, v, k, l));
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pseudoscalar-Pseudoscalar vertex function: i\Gamma_{\mu\nu}(p',p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_Ppsps(const M4Vec &prime,
                                                             const M4Vec &p) const {
  double beta_P = 0.0;
  if (prime.M() < 0.2) {  // Choose based on particle mass
    // Coupling (pion)
    beta_P = 1.76;  // GeV^{-1}
  } else {
    // Coupling (kaon)
    beta_P = 1.54;  // GeV^{-1}
  }

  const M4Vec  psum = prime + p;
  const double M2   = psum.M2();

  Tensor2<std::complex<double>, 4, 4> T;
  const std::complex<double> FACTOR = -zi * 2.0 * beta_P * FM((prime - p).M2());

  FOR_EACH_2(LI);
  T(u,v) = FACTOR * ((psum % u) * (psum % v) - 0.25 * g[u][v] * M2);
  FOR_EACH_2_END;
  
  return T;
}

// Pomeron-Vector(massive)-Vector(massive) vertex function:
// i\Gamma_{\alpha\beta\gamma\delta} (p',p)
//
// Input M0 is the vector meson on-shell mass
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_Pvv(const M4Vec &prime,
                                                                 const M4Vec &p) const {
  // Couplings (phi meson)
  const double a = 0.49;  // GeV^{-3}
  const double b = 4.27;  // GeV^{-1}

  const Tensor4<std::complex<double>, 4, 4, 4, 4> G0 = Gamma0(prime, -p);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> G2 = Gamma2(prime, -p);

  const std::complex<double> FACTOR = zi * FM((prime - p).M2());

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (2.0 * a * G0(u, v, k, l) - b * G2(u, v, k, l));
  FOR_EACH_4_END;

  return T;
}

// Gamma-Vector meson vertex
//
// type = "rho", "omega", "phi"
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::yV(std::string type) const {
  const double e      = msqrt(alpha_QED * 4.0 * PI);  // ~ 0.3
  double       gammaV = 0.0;
  double       mV     = 0.0;

  if (type == "rho") {
    gammaV = msqrt(4 * PI / 0.496);
    mV     = 0.770;
  }
  if (type == "omega") {
    gammaV = msqrt(4 * PI / 0.042);
    mV     = 0.785;
  }
  if (type == "phi") {
    gammaV = (-1.0) * msqrt(4 * PI / 0.0716);
    mV     = 1.020;
  }

  Tensor2<std::complex<double>, 4, 4> T;

  FOR_EACH_2(LI);
  T(u,v) = -zi * e * pow2(mV) / gammaV * g[u][v];
  FOR_EACH_2_END;

  return T;
}

// ----------------------------------------------------------------------
// Propagators

// Tensor Pomeron propagator: (covariant == contravariant)
//
// i\Delta_{\mu\nu,\kappa\lambda}(s,t) = i\Delta^{\mu\nu,\kappa\lambda}(s,t)
//
// Input as (sub)-Mandelstam invariants: s,t
//
//
// Symmetry relations:
// \Delta_{\mu\nu,\kappa\lambda} = \Delta_{\nu\mu,\kappa\lambda} =
// \Delta_{\mu\nu,\lambda\kappa} = \Delta_{\kappa\lambda,\mu\nu}
//
// Contraction with g^{\mu\nu} or g^{\kappa\lambda} gives 0
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_P(double s, double t) const {

  const std::complex<double> FACTOR = 1.0 / (4.0 * s) * std::pow(-zi * s * ap_P, alpha_P(t) - 1.0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Tensor Reggeon propagator iD (f2-a2-reggeon)
//
// Input as (sub)-Mandelstam invariants: s,t
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_2R(double s, double t) const {

  const std::complex<double> FACTOR =
      1.0 / (4.0 * s) * std::pow(-zi * s * ap_2R, alpha_2R(t) - 1.0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Vector Reggeon propagator iD (rho-reggeon, omega-reggeon)
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_1R(double s, double t) const {

  const std::complex<double> FACTOR =
      zi / pow2(M_1R) * std::pow(-zi * s * ap_1R, alpha_1R(t) - 1.0);

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR * g[u][v];
  FOR_EACH_2_END;

  return T;
}

// Odder propagator iD
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_O(double s, double t) const {

  const std::complex<double> FACTOR =
      -zi * eta_O / pow2(M_O) * std::pow(-zi * s * ap_O, alpha_O(t) - 1.0);

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR * g[u][v];
  FOR_EACH_2_END;

  return T;
}


// Vector propagator iD_{\mu \nu} == iD^{\mu \nu} with (naive) Reggeization
//
// Input as contravariant 4-vector and the particle peak mass M0
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_V(const M4Vec &p, double M0, double s34) const {

  const double m2       = p.M2();
  const double alpha    = 0.1 + 0.9 * m2;
  const double s_thresh = 4*pow2(M0);
  const double phi      = PI/2.0 * std::exp((s_thresh - s34) / s_thresh) - PI/2.0;

  // Reggeization
  std::complex<double> REGGEIZE = std::pow(std::exp(zi * phi) * s34 / s_thresh, alpha - 1.0);

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u,v) = -g[u][v] / (m2 - M0*M0) * REGGEIZE;
  FOR_EACH_2_END;

  return T;
}

// Scalar propagator iD with zero width
//
// Input as 4-vector and the particle peak mass M0
//
std::complex<double> MTensorPomeron::iD_MES0(const M4Vec &p, double M0) const {
  return zi / (p.M2() - pow2(M0));
}

// Scalar propagator iD with finite width
//
// Input as 4-vector and the particle peak mass M0 and full width
//
std::complex<double> MTensorPomeron::iD_MES(const M4Vec &p, double M0, double Gamma) const {
  const double m2 = p.M2();

  // Breit-Wigner
  const std::complex<double> Delta = 1.0 / (m2 - pow2(M0) + zi * M0 * Gamma);
  return zi * Delta;
}

// Massive vector propagator iD_{\mu \nu} or iD^{\mu \nu} with INDEX_UP == true
//
// Input as contravariant 4-vector and the particle peak mass M0 and full width Gamma
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_VMES(const M4Vec &p, double M0,
                                                            double Gamma, bool INDEX_UP) const {
  // Transverse part
  const double               m2      = p.M2();
  const std::complex<double> Delta_T = 1.0 / (m2 - pow2(M0) + zi * M0 * Gamma);

  // Longitudinal part (does not enter here)
  // const double Delta_L = 0.0;

  Tensor2<std::complex<double>, 4, 4> T;

  if (INDEX_UP) {
    FOR_EACH_2(LI);
    T(u,v) = zi * (-g[u][v] + (p [u]) * (p [v]) / m2) * Delta_T; // - zi*(p[u])*(p[v])/m2 * Delta_L; // this part does not enter
    FOR_EACH_2_END;
  } else {
    FOR_EACH_2(LI);
    T(u,v) = zi * (-g[u][v] + (p % u) * (p % v) / m2) * Delta_T; // - zi*(p%u)*(p%v)/m2 * Delta_L; // this part does not enter
    FOR_EACH_2_END;
  }

  return T;
}

// Massive tensor propagator iD_{\mu \nu \kappa \lambda} or iD^{} with INDEX_UP == true
//
// Input as contravariant (upper index) 4-vector and the particle peak mass M0
// and full width Gamma
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_TMES(const M4Vec &p, double M0,
                                                                  double Gamma, bool INDEX_UP) const {
  const double m2 = p.M2();
  
  // Inline auxialary function
  auto ghat = [&](int u, int v) {
    if (INDEX_UP) {
      return -g[u][v] + (p [u]) * (p [v]) / m2;
    } else {
      return -g[u][v] + (p % u) * (p % v) / m2;
    }
  };

  // Breit-Wigner
  const std::complex<double> Delta = 1.0 / (m2 - pow2(M0) + zi * M0 * Gamma);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = zi * Delta * (0.5 * (ghat(u, k) * ghat(v, l) + ghat(u, l) * ghat(v, k)) -
                                1.0 / 3.0 * ghat(u, v) * ghat(k, l));
  FOR_EACH_4_END;

  return T;
}

// ----------------------------------------------------------------------
// Trajectories (simple linear/affine here)

// Pomeron trajectory
double MTensorPomeron::alpha_P(double t) const { return 1.0 + delta_P + ap_P * t; }
// Odderon trajectory
double MTensorPomeron::alpha_O(double t) const { return 1.0 + delta_O + ap_O * t; }
// Reggeon (rho,omega) trajectory
double MTensorPomeron::alpha_1R(double t) const { return 1.0 + delta_1R + ap_1R * t; }
// Reggeon (f2,a2) trajectory
double MTensorPomeron::alpha_2R(double t) const { return 1.0 + delta_2R + ap_2R * t; }

// ----------------------------------------------------------------------
// Form factors

// Proton electromagnetic Dirac form FACTOR (electric)
double MTensorPomeron::F1(double t) const {
  return (1 - (t / (4 * mp * mp)) * mu_ratio) / (1 - t / (4 * mp * mp)) * GD(t);
}

// Proton electromagnetic Pauli form FACTOR (magnetic)
double MTensorPomeron::F2(double t) const {
  return (mu_ratio - 1) / (1 - t / (4 * mp * mp)) * GD(t);
}

// Dipole form FACTOR
double MTensorPomeron::GD(double t) const {
  const double m2D = 0.71;  // GeV^2
  return 1.0 / pow2(1 - t / m2D);
}

// Meson form FACTOR (pion electromagnetic type)
double MTensorPomeron::FM(double q2) const {
  const double A2 = 0.5;  // GeV^2
  return 1.0 / (1.0 - q2 / A2);
}

// ----------------------------------------------------------------------
// Tensor functions

// Output: \Gamma^{(0)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::Gamma0(const M4Vec &k1,
                                                                 const M4Vec &k2) const {
  const double k1k2 = k1 * k2;  // 4-dot product
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = (k1k2 * g[u][v] - (k2 % u) * (k1 % v)) *
                  ((k1 % k) * (k2 % l) + (k2 % k) * (k1 % l) - 0.5 * k1k2 * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Output: \Gamma^{(2)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::Gamma2(const M4Vec &k1,
                                                                 const M4Vec &k2) const {
  const double k1k2 = k1 * k2;  // 4-dot product
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = k1k2 * (g[u][k] * g[v][l] + g[u][l] * g[v][k]) +
                  g[u][v] * ((k1 % k) * (k2 % l) + (k2 % k) * (k1 % l)) -
                  (k1 % v) * (k2 % l) * g[u][k] - (k1 % v) * (k2 % k) * g[u][l] -
                  (k2 % u) * (k1 % l) * g[v][k] - (k2 % u) * (k1 % k) * g[v][l] -
                  (k1k2 * g[u][v] - (k2 % u) * (k1 % v)) * g[k][l];
  FOR_EACH_4_END;

  return T;
}

// Pre-calculate tensors
//
// 1. Minkowski metric tensor
// 2. Auxialary (helper) tensor
// 1/2 g_{\mu\kappa} g_{\nu\lambda} + 1/2 g_{\mu\lambda}g_{\nu\kappa} - 1/4
// g_{\mu\nu}g_{kappa\lambda}
//
void MTensorPomeron::CalcRTensor() {

  // Minkowski metric tensor (+1,-1,-1,-1)
  FOR_EACH_2(LI);
  gT(u,v) = g[u][v];
  FOR_EACH_2_END;

  // ------------------------------------------------------------------
  // Covariant [lower index] 4D-epsilon tensor
  
  const MTensor<int> etensor = math::EpsTensor(4);

  FOR_EACH_4(LI);
  std::vector<std::size_t> ind = {u, v, k, l};
  eps_lo(u, v, k, l) = static_cast<double>(etensor(ind));
  FOR_EACH_4_END;

  // Free Lorentz indices [second parameter denotes the range of index]
  FTensor::Index<'a', 4> a;
  FTensor::Index<'b', 4> b;
  FTensor::Index<'c', 4> c;
  FTensor::Index<'d', 4> d;
  FTensor::Index<'g', 4> alfa;
  FTensor::Index<'h', 4> beta;

  // Contravariant version (make it in two steps)
  eps_hi(a, b, c, d) = eps_lo(alfa, beta, c, d) * gT(a, alfa) * gT(b, beta);
  eps_hi(a, b, c, d) = eps_hi(a, b, alfa, beta) * gT(c, alfa) * gT(d, beta);
  
  // ------------------------------------------------------------------

  // Aux tensor R
  FTensor::Tensor4<double, 4, 4, 4, 4> R;

  R(a, b, c, d) =
      0.5 * gT(a, c) * gT(b, d) + 0.5 * gT(a, d) * gT(b, c) + 0.25 * gT(a, b) * gT(c, d);

  // Different mixed index covariant/contravariant versions
  // by contraction with g_{\mu\nu}
  R_DDDD = R;
  R_DDDU(a, b, c, d) = R(a, b, c, alfa) * gT(d, alfa);
  R_DDUU(a, b, c, d) = R(a, b, alfa, beta) * gT(c, alfa) * gT(d, beta);
  R_UUDD(a, b, c, d) = R(alfa, beta, c, d) * gT(a, alfa) * gT(b, beta);

  // -------------------------------------------------------------------
  // Pre-calculated tensor contractions for tensor resonances

  auto FIX1 = [&] () {
  MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4}, double(0.0));
  FOR_EACH_6(LI);
  const std::vector<size_t> ind = {u, v, k, l, r, s};
  for (auto const &a1 : LI) {
    for (auto const &r1 : LI) {
      for (auto const &s1 : LI) {
          A(ind) += R_DDDD(u, v, r1, a1) * R_DDDU(k, l, s1, a1) * R_DDUU(r, s, r1, s1);
      }
    }
  }
  FOR_EACH_6_END;
  return A;
  };

  auto FIX2 = [&] () {
  MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4, 4, 4}, double(0.0));
  FOR_EACH_6(LI);
  for (auto const &a1 : LI) {
    for (auto const &r1 : LI) {
      for (auto const &s1 : LI) {
        for (auto const &u1 : LI) {

          A({u,v,k,l,r,s,u1,r1}) += R_DDDD(u, v, u1, a1) * R_DDDU(k, l, s1, a1) * R_DDUU(r, s, r1, s1);
        }
      }
    }
  }
  FOR_EACH_6_END;
  return A;
  };

  auto FIX3 = [&] () {
  MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4, 4, 4}, double(0.0));
  FOR_EACH_6(LI);
  for (auto const &a1 : LI) {
    for (auto const &r1 : LI) {
      for (auto const &s1 : LI) {
        for (auto const &u1 : LI) {

          A({u,v,k,l,r,s,u1,s1}) += R_DDDD(u, v, r1, a1) * R_DDDU(k, l, u1, a1) * R_DDUU(r, s, r1, s1);
        }
      }
    }
  }
  FOR_EACH_6_END;
  return A;
  };

  T1 = FIX1();
  T2 = FIX2();
  T3 = FIX3();
}

}  // gra namespace ends
