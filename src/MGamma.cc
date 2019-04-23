// Gamma-Gamma Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MForm.h"
#include "Graniitti/MGamma.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"
#include "Graniitti/MPDG.h"
#include "Graniitti/MRegge.h"
#include "Graniitti/MSpin.h"

using gra::math::msqrt;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;
using gra::math::zi;

using namespace gra::form;

namespace gra {
// ============================================================================
// (yy -> fermion-antifermion pair)
// yy -> e+e-, mu+mu-, tau+tau-, qqbar or Monopole-Antimonopole (spin-1/2
// monopole)
// production
//
//
// This is the same amplitude as yy -> e+e- (two diagrams)
// obtained by crossing e+e- -> yy annihilation, see:
//
// http://theory.sinp.msu.ru/comphep_old/tutorial/QED/node4.html
//
// beta = v/c of the monopole (or antimonopole) in their CM system
//
//
// *************************************************************
// Dirac quantization condition:
//
// g = 2\pi \hbar n / (\mu_0 e),    where n = 1,2,3,...
// where alpha = e^2/(4*PI) is the running QED coupling
//
//
// Numerical values:
// alpha_g  = g^2 / (4*PI) ~ 34  (when n = 1)
// alpha_em = e^2 / (4*PI) ~ 1/137
// *************************************************************
//
// From lepton pair to monopole pair:
// replace e -> g*beta
//
// 16 * pow2(PI * alpha_EM) -> pow(g*beta, 4)
//
// *************************************************************
// Coupling schemes:
//
// Dirac:     alpha_g = g^2/(4pi)
// Beta-dirac alpha_g = (g*beta)^2 / (4pi)
//
// [REFERENCE: Rajantie,
// https://physicstoday.scitation.org/doi/pdf/10.1063/PT.3.3328]
// [REFERENCE: Dougall, Wick, https://arxiv.org/abs/0706.1042]
// [REFERENCE: Rels, Sauter, https://arxiv.org/abs/1707.04170v1]
//
std::complex<double> MGamma::yyffbar(gra::LORENTZSCALAR &lts) {
  // QED couplings
  double COUPL = 16.0 * pow2(gra::math::PI * gra::form::alpha_EM(0));  // = e^4
  const double mass = lts.decaytree[0].p4.M();  // lepton, quark (or monopole) mass
  const double mass2 = pow2(mass);
  const bool MONOPOLE_MODE = (lts.decaytree[0].p.pdg == PDG::PDG_monopole) ? true : false;

  if (PARAM_MONOPOLE::gn < 1) {
    throw std::invalid_argument("MGamma::yyffbar: Parameter Dirac n less than 1");
  }

  // Amplitude
  std::complex<double> A = 0.0;

  // Monopole-Antimonopole coupling
  if (MONOPOLE_MODE) {
    static const double g = 2.0 * math::PI * PARAM_MONOPOLE::gn / form::e_EM();

    if (PARAM_MONOPOLE::coupling == "beta-dirac") {
      // Calculate beta (velocity)
      // const M4Vec px = pfinal[3] + pfinal[4];
      // M4Vec p3 = pfinal[3];
      // LorentzBoost(px, px.M(), p3, -1); // // Boost to the
      // CM frame of M\bar{M} pair
      // double beta = std::sqrt( p3.Px()*p3.Px() +
      // p3.Py()*p3.Py() + p3.Pz()*p3.Pz()) / p3.E();

      // Faster way
      const double beta = msqrt(1.0 - 4.0 * pow2(PARAM_MONOPOLE::M0) / lts.s_hat);
      COUPL = pow4(g * beta);
    } else if (PARAM_MONOPOLE::coupling == "dirac") {
      COUPL = pow4(g);
    } else {
      throw std::invalid_argument("MGamma::yyffbar: Unknown PARAM_MONOPOLE::coupling " +
                                  PARAM_MONOPOLE::coupling);
    }

    // QED tree level amplitude squared |M|^2, spin averaged and
    // summed
    const double amp2 =
        2.0 * COUPL *
        ((lts.u_hat - mass2) / (lts.t_hat - mass2) + (lts.t_hat - mass2) / (lts.u_hat - mass2) +
         1.0 -
         pow2(1.0 + (2.0 * mass2) / (lts.t_hat - mass2) + (2.0 * mass2) / (lts.u_hat - mass2)));

    /*
      // FeynCalc result, less simplified, but exactly same result
      as above
            double amp2_ = 2 * COUPL  *
                      ( pow4(mass)*(3*pow2(lts.t_hat) +
      14*lts.t_hat*lts.u_hat + 3*pow2(lts.u_hat))
                            -mass2*(pow3(lts.t_hat) +
      7*pow2(lts.t_hat)*lts.u_hat + 7*pow2(lts.u_hat)*lts.t_hat +
      pow3(lts.u_hat) )
                            -6*pow8(mass) +
      lts.t_hat*lts.u_hat*(pow2(lts.t_hat) + pow2(lts.u_hat)) )/ (
      pow2(lts.t_hat - mass2) * pow2(lts.u_hat - mass2) );
    */

    // printf("%0.15f %0.15f \n", amp2, amp2_);

    A = msqrt(amp2);
  } else {  // MADGRAPH/HELAS, all helicity amplitudes individually -> for
    // the screening loop

    A = AmpMG5_yy_ll.CalcAmp(lts);
  }

  // ------------------------------------------------------------------
  // quark pair (charge 1/3 or 2/3), apply charge and color factors
  if (std::abs(lts.decaytree[0].p.pdg) <= 6) {  // we have a quark

    const double Q = lts.decaytree[0].p.chargeX3 / 3.0;
    const double NC = 3.0;  // quarks come in three colors

    const double factor = msqrt(pow4(Q) * NC);  // to amplitude level
    A *= factor;

    for (const auto &i : aux::indices(lts.hamp)) {  // helicity amplitudes
      lts.hamp[i] *= factor;
    }
  }
  // ------------------------------------------------------------------

  return A;
}

// --------------------------------------------------------------------------------------------
// A Monopolium (monopole-antimonopole) bound state process
//
// See e.g.
//
// [REFERENCE: Preskill,
// http://www.theory.caltech.edu/~preskill/pubs/preskill-1984-monopoles.pdf]
// [REFERENCE: Epele, Franchiotti, Garcia, Canal, Vento,
// https://arxiv.org/abs/hep-ph/0701133v2]
// [REFERENCE: Barrie, Sugamoto, Yamashita, https://arxiv.org/abs/1607.03987v3]
// [REFERENCE: Fanchiotti, Canal, Vento, https://arxiv.org/pdf/1703.06649.pdf]
// [REFERENCE: Reis, Sauter, https://arxiv.org/abs/1707.04170v1]
//
// *************************************************************
// Dirac quantization condition:
//
// g = 2\pi \hbar / (\mu_0 e) n,    where n = 1,2,3,...
// where alpha = e^2/(4*PI) is the running QED coupling
//
//
// Numerical values:
// alpha_g  = g^2 / (4*PI) ~ 34  (when n = 1)
// alpha_em = e^2 / (4*PI) ~ 1/137
// *************************************************************
//
std::complex<double> MGamma::yyMP(const gra::LORENTZSCALAR &lts) const {
  gra::g_mutex.lock();  // @@@
  PARAM_MONOPOLE::PrintParam(lts.sqrt_s);
  gra::g_mutex.unlock();  // @@@

  // Monopolium nominal mass and width parameters
  static const double M = 2.0 * PARAM_MONOPOLE::M0 + PARAM_MONOPOLE::EnergyMP(PARAM_MONOPOLE::En);
  static const double Gamma_M = PARAM_MONOPOLE::Gamma0;

  if (M < 0) {
    throw std::invalid_argument(
        "MGamma::yyMP: Increase ladder level parameter En. Monopolium "
        "nominal mass " +
        std::to_string(M) + " < 0!");
  }
  if (PARAM_MONOPOLE::gn < 1) {
    throw std::invalid_argument("MGamma::yyMP: Parameter Dirac n less than 1");
  }

  // Two coupling scenarios:
  static const double g = 2.0 * math::PI * PARAM_MONOPOLE::gn / form::e_EM();
  double beta = 0.0;

  if (PARAM_MONOPOLE::coupling == "beta-dirac") {
    beta = msqrt(1.0 - pow2(M) / lts.s_hat);
  } else if (PARAM_MONOPOLE::coupling == "dirac") {
    beta = 1.0;
  } else {
    throw std::invalid_argument("MGamma::yyMP: Unknown PARAM_MONOPOLE::coupling " +
                                PARAM_MONOPOLE::coupling);
  }

  // Magnetic coupling
  const double alpha_g = pow2(beta * g) / (4.0 * math::PI);

  // Running width
  const double Gamma_E = PARAM_MONOPOLE::GammaMP(PARAM_MONOPOLE::En, alpha_g);

  // printf("alpha_g = %0.3E, Gamma_E = %0.3E, Gamma_M = %0.3E, Psi_MP = %0.3E
  // \n",
  //	alpha_g, Gamma_E, Gamma_M, PARAM_MONOPOLE::PsiMP(PARAM_MONOPOLE::n));

  // Normalization
  double norm = 2.0 * math::PI * M * M;

  // Sub process
  double sigma_hat = norm * (Gamma_E * Gamma_M) / (pow2(lts.s_hat - M * M) + pow2(M * Gamma_M));

  return msqrt(sigma_hat);
}

// --------------------------------------------------------------------------------------------
// Gamma-Gamma to SM-Higgs 0++ helicity amplitudes
//
// Generic narrow width yy -> X cross section in terms of partial decay widths:
//
// \sigma(yy -> X) = 8\pi^2/M_X (2J+1) \Gamma(X -> yy) \delta(shat - M_X^2) (1 +
// h1h2)
//                 = (8 * \pi)  (2J+1) \Gamma(X -> yy) \Gamma_X (1 + h1h2) /
//                 ((shat - M_X^2)^2 +
//                 M_X^2\Gamma_X^2),
//
// where h1,h2 = +- gamma helicities
//
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0111078]
// [REFERENCE: Bernal, Lopez-Val, Sola, https://arxiv.org/pdf/0903.4978.pdf]
// [REFERENCE: Enterria, Lansberg,
// https://www.slac.stanford.edu/pubs/slacpubs/13750/slac-pub-13786.pdf]
//
std::complex<double> MGamma::yyHiggs(gra::LORENTZSCALAR &lts) const {
  lts.hamp.resize(4);

  const double M = 125.18;                  // Higgs mass measured (GeV)
  const double Gamma = 0.00415;             // Higgs total width calculated (GeV)
  const double Gamma_yy = 2.27E-3 * Gamma;  // Higgs to gamma-gamma calculated (GeV)

  // Normalization
  const double norm = 2.0 * math::PI * M * M;

  // Effective helicity amplitudes
  lts.hamp[0] =
      msqrt(norm * Gamma_yy * Gamma * (1 + 1) / (pow2(lts.s_hat - M * M) + pow2(M * Gamma)));  // --
  lts.hamp[1] = 0.0;                                                                           // -+
  lts.hamp[2] = 0.0;                                                                           // +-
  lts.hamp[3] = lts.hamp[0];                                                                   // ++

  // Sum over helicity amplitudes squared
  double sumA2 = 0.0;
  for (const auto &i : aux::indices(lts.hamp)) {
    sumA2 += gra::math::abs2(lts.hamp[i]);
  }
  sumA2 /= 4;  // Initial state polarization average

  return msqrt(sumA2);  // We take square later
}

// ============================================================================
// Simple matrix element ansatz for gamma-gamma resonances
//
// ===========
//      $
//      $
//      x---->
//      $
//      $
// ===========
//
// For Narrow-Width approximations, see:
//
// [REFERENCE: Uhlemann, Kauer, Narrow-width approximation accuracy,
// https://arxiv.org/pdf/0807.4112.pdf]
//
std::complex<double> MGamma::yyX(const gra::LORENTZSCALAR &lts, gra::PARAM_RES &resonance) const {
  // Factor of 2 x from (identical) initial state boson statistics
  std::complex<double> A_prod =
      2.0 * gra::form::CBW(lts, resonance) * resonance.g *
      PARAM_REGGE::ResonanceFormFactor(lts.m2, pow2(resonance.p.mass), resonance.g_FF) *
      PARAM_REGGE::JPCoupling(lts, resonance);

  // Spin and decay part
  const std::complex<double> A_decay = gra::spin::SpinAmp(lts, resonance);

  return A_prod * A_decay;
}

}  // gra namespace
