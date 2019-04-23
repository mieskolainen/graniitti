// 'Durham QCD' Processes and Amplitudes
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/MAux.h"
#include "Graniitti/MDurham.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MPDG.h"
#include "Graniitti/MSudakov.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

namespace gra {
using math::msqrt;
using math::pow2;
using math::zi;
using math::PI;
using PDG::mp;

// 2-body incoming or outgoing helicity combinations
// (keep it algebraic order to match with MadGraph)
constexpr int MM = 0;  // --
constexpr int MP = 1;  // -+
constexpr int PM = 2;  // +-
constexpr int PP = 3;  // ++

// 2-body final state helicity configurations
static const std::vector<int> fs2 = {MM, MP, PM, PP};

// Central system J^P = 0^+, 0^-, +2^+, -2^+
enum SPINPARITY { P0, M0, P2, M2 };  // Implicit conversion to int

// Durham loop integral discretization technical parameters
namespace Durham {
unsigned int N_qt = 0;   // Number of qt discretization intervals
unsigned int N_phi = 0;  // Number of phi discretization intervals

double qt2_MIN = 0;  // Loop momentum qt^2 minimum (GeV^2)
double qt2_MAX = 0;  // Loop momentum qt^2 maximum (GeV^2)

std::string PDF_scale = "MIN";  // Scheme
double alphas_scale = 4.0;      // PDF factorization scale
double MAXCOS = 0.9;            // Meson amplitude |cos(theta*)| < MAXCOS

void ReadParameters() {
  using json = nlohmann::json;

  const std::string inputfile =
      gra::aux::GetBasePath(2) + "/modeldata/" + gra::MODELPARAM + "/GENERAL.json";
  const std::string data = gra::aux::GetInputData(inputfile);
  json j;

  try {
    j = json::parse(data);

    // JSON block identifier
    const std::string XID = "PARAM_DURHAM_QCD";
    N_qt = j.at(XID).at("N_qt");
    N_phi = j.at(XID).at("N_phi");
    qt2_MIN = j.at(XID).at("qt2_MIN");
    qt2_MAX = j.at(XID).at("qt2_MAX");
    PDF_scale = j.at(XID).at("PDF_scale");
    alphas_scale = j.at(XID).at("alphas_scale");
    MAXCOS = j.at(XID).at("MAXCOS");

    // Now calculate rest
    qt_MIN = msqrt(Durham::qt2_MIN);
    qt_MAX = msqrt(Durham::qt2_MAX);

    qt_STEP = (Durham::qt_MIN - Durham::qt_MAX) / Durham::N_qt;
    phi_STEP = (2.0 * math::PI) / Durham::N_phi;

    initialized = true;
  } catch (...) {
    std::string str =
        "Durham::ReadParameters: Error parsing " + inputfile + " (Check for extra/missing commas)";
    throw std::invalid_argument(str);
  }
}

// THESE ARE CALCULATED FROM ABOVE
double qt_MIN = 0;
double qt_MAX = 0;
double qt_STEP = 0;
double phi_STEP = 0;
bool initialized = false;
}

// Durham QCD / KMR model
//
// [REFERENCE: Pumplin, Phys.Rev.D 52 (1995)]
// [REFERENCE: Khoze, Kaidalov, Martin, Ryskin, Stirling,
// https://arxiv.org/abs/hep-ph/0507040]
// [REFERENCE: Khoze, Martin, Ryskin, https://arxiv.org/abs/hep-ph/0605113]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, Stirling,
// https://arxiv.org/abs/1005.0695]
// [REFERENCE: Harland-Lang, Khoze, Ryskin, https://arxiv.org/abs/1409.4785]
//
std::complex<double> MDurham::DurhamQCD(gra::LORENTZSCALAR &lts, const std::string &process) {
  // First run, init parameters
  // @@ MULTITHREADING LOCK @@
  gra::g_mutex.lock();
  if (gra::GlobalSudakovPtr->initialized == false) {
    try {
      gra::GlobalSudakovPtr->Init(lts.sqrt_s, gra::LHAPDFSET, true);
    } catch (...) {
      gra::g_mutex.unlock();  // need to release here, otherwise get infinite lock
      throw;
    }
  }
  if (!Durham::initialized) {
    try {
      Durham::ReadParameters();
    } catch (...) {
      gra::g_mutex.unlock();  // need to release here, otherwise get infinite lock
      throw;
    }
  }
  gra::g_mutex.unlock();

  if (process == "gg") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    std::vector<std::vector<std::complex<double>>> Amp(4,
                                                       std::vector<std::complex<double>>(4, 0.0));

    // Madgraph
    gra::g_mutex.lock();
    const double alpha_s = gra::GlobalSudakovPtr->AlphaS_Q2(lts.s_hat);
    gra::g_mutex.unlock();
    AmpMG5_gg_gg.CalcAmp(lts, alpha_s);

    // Amplitude evaluated outside the Qt-loop (approximation)
    Dgg2gg(lts, Amp);

    // Run loop
    return DQtloop(lts, Amp);
  } else if (process == "qqbar") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    std::vector<std::vector<std::complex<double>>> Amp(4,
                                                       std::vector<std::complex<double>>(4, 0.0));

    // Madgraph
    gra::g_mutex.lock();
    const double alpha_s = gra::GlobalSudakovPtr->AlphaS_Q2(lts.s_hat);
    gra::g_mutex.unlock();
    AmpMG5_gg_qqbar.CalcAmp(lts, alpha_s);

    // Amplitude evaluated outside the Qt-loop (approximation)
    Dgg2qqbar(lts, Amp);

    // Run loop
    return DQtloop(lts, Amp);
  } else if (process == "MMbar") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    std::vector<std::vector<std::complex<double>>> Amp(1,
                                                       std::vector<std::complex<double>>(4, 0.0));

    // Amplitude evaluated outside the Qt-loop (approximation)
    Dgg2MMbar(lts, Amp);

    // Run loop
    return DQtloop(lts, Amp);
  } else if (process == "chic(0)") {
    // [Final state helicities/polarizations x 4 initial state helicities]
    std::vector<std::vector<std::complex<double>>> Amp(1,
                                                       std::vector<std::complex<double>>(4, 0.0));

    // Amplitude evaluated outside the Qt-loop (approximation)
    std::vector<double> null;
    Dgg2chic0(lts, Amp, null, null);

    // Run loop
    return DQtloop(lts, Amp);

    // ------------------------------------------------------------
    //
    // Implement more processes here ...
    //
    // ------------------------------------------------------------
  } else {
    std::string str = "MDurham::Durham: Unknown subprocess: " + process;
    throw std::invalid_argument(str);
  }
}

// Helicity basis decomposition
//
// In the forward limit q1_t = - q2_t = Q_t
// => gives gluon polarization vectors eps_1 = -eps_2 => central system J_z = 0
//
inline void MDurham::DHelicity(const std::vector<double> &q1, const std::vector<double> &q2,
                               std::vector<std::complex<double>> &JzP) const {
  const unsigned int X = 0;  // component for readability
  const unsigned int Y = 1;

  // 1/2  q1_t dot q2_t
  JzP[P0] = -0.5 * (q1[X] * q2[X] + q1[Y] * q2[Y]);

  // 1/2 i |q1_t x q2_t|
  JzP[M0] = -0.5 * math::zi * std::abs(q1[X] * q2[Y] - q1[Y] * q2[X]);

  const double Re = 0.5 * (q1[X] * q2[X] - q1[Y] * q2[Y]);
  const double Im = 0.5 * (q1[X] * q2[Y] + q1[Y] * q2[X]);

  // +1/2[ (xx - yy) + i(xy + yx) ]
  JzP[P2] = Re + math::zi * Im;

  // +1/2[ (xx - yy) - i(xy + yx) ]
  JzP[M2] = Re - math::zi * Im;
}

// M_{\lambda_1,\lambda_2}
// i.e. g(\lambda_1) g(\lambda_2) -> X (system) helicity amplitudes
//
// where lambda1, lambda2 are the initial state gluon helicities
//
// [REFERENCE: https://arxiv.org/abs/1405.0018v2, formula (20)]
//
inline std::complex<double> MDurham::DHelProj(const std::vector<std::complex<double>> &A,
                                              const std::vector<std::complex<double>> &JzP) const {
  // M_{++} + M_{--}
  // (this term gives J_z^PC = 0^++ selection rule in the forward pt->0 limit)
  const std::complex<double> aP0 = JzP[P0] * (A[PP] + A[MM]);
  // M_{++} - M_{--}
  const std::complex<double> aM0 = JzP[M0] * (A[PP] - A[MM]);
  // M_{-+}
  const std::complex<double> aP2 = JzP[P2] * A[MP];
  // M_{+-}
  const std::complex<double> aM2 = JzP[M2] * A[PM];

  /*
  // DEBUG
  std::cout << "0+ : " << aP0 << std::endl;
  std::cout << "0- : " << aM0 << std::endl;
  std::cout << "2+ : " << aP2 << std::endl;
  std::cout << "2- : " << aM2 << std::endl;
  std::cout << std::endl << std::endl;
  */

  return aP0 + aM0 + aP2 + aM2;
}

/*
// MADGRAPH HELICITY FORMAT [direct algebraic order]
const static int helicities[ncomb][nexternal] =
  {{-1, -1, -1, -1}, 0
   {-1, -1, -1,  1}, 1
   {-1, -1,  1, -1}, 2
   {-1, -1,  1,  1}, 3

   {-1,  1, -1, -1}, 4
   {-1,  1, -1,  1}, 5
   {-1,  1,  1, -1}, 6
   {-1,  1,  1,  1}, 7

   { 1, -1, -1, -1}, 8
   { 1, -1, -1,  1}, 9
   { 1, -1,  1, -1}, 10
   { 1, -1,  1,  1}, 11

   { 1,  1, -1, -1}, 12
   { 1,  1, -1,  1}, 13
   { 1,  1,  1, -1}, 14
   { 1,  1,  1,  1}}; 15
*/

// Alternative (semi-ad-hoc) scenarios for the scale choise
inline void MDurham::DScaleChoise(double qt2, double q1_2, double q2_2, double &Q1_2_scale,
                                  double &Q2_2_scale) const {
  if (Durham::PDF_scale == "MIN") {
    Q1_2_scale = std::min(qt2, q1_2);
    Q2_2_scale = std::min(qt2, q2_2);
  } else if (Durham::PDF_scale == "MAX") {
    Q1_2_scale = std::max(qt2, q1_2);
    Q2_2_scale = std::max(qt2, q2_2);
  } else if (Durham::PDF_scale == "IN") {
    Q1_2_scale = q1_2;
    Q2_2_scale = q2_2;
  } else if (Durham::PDF_scale == "EX") {
    Q1_2_scale = qt2;
    Q2_2_scale = qt2;
  } else if (Durham::PDF_scale == "AVG") {
    Q1_2_scale = (qt2 + q1_2) / 2.0;
    Q2_2_scale = (qt2 + q2_2) / 2.0;
  } else {
    throw std::invalid_argument("MDurham::DScaleChoise: Unknown 'Durham::PDF_scale' option!");
  }
}

// [REFERENCE: Khoze, Martin, Ryskin,
// https://journals.aps.org/prd/pdf/10.1103/PhysRevD.56.5867]
// [REFERENCE: Harland-Lang, Khoze, Martin, Ryskin, Stirling,
// https://arxiv.org/abs/1405.0018v2]
//
// See also: <https://arxiv.org/pdf/1608.03765.pdf> (LUND)
//
//
// Durham loop integral amplitude:
// A = pi^2 \int \frac{d^2 Q_t M(gg->X)}{Q_t^2(Q_t-{p_1t})^2(Q_t+p_{2t})^2}
//               x f_g(x_1,x_1',Q_1^2,\mu_2;t_1) x
//               f_g(x_2,x_2',Q_2^2,\mu_2;t_2)
//
std::complex<double> MDurham::DQtloop(gra::LORENTZSCALAR &lts,
                                      std::vector<std::vector<std::complex<double>>> Amp) {
  // Forward (proton) system pt-vectors
  const std::vector<double> pt1 = {lts.pfinal[1].Px(), lts.pfinal[1].Py()};
  const std::vector<double> pt2 = {lts.pfinal[2].Px(), lts.pfinal[2].Py()};

  // *************************************************************************
  // ** Process scale (GeV) / Sudakov suppression kt^2 integral upper bound **
  const double MU = msqrt(lts.s_hat / Durham::alphas_scale);
  // *************************************************************************

  // Init 2D-Simpson weight matrix (will be calculated only once, being
  // static), C++11 handles multithreaded static initialization
  const static MMatrix<double> WSimpson = math::Simpson38Weight2D(Durham::N_qt, Durham::N_phi);

  // NOTE N + 1, init with zero!
  std::vector<MMatrix<std::complex<double>>> f(
      Amp.size(), MMatrix<std::complex<double>>(Durham::N_qt + 1, Durham::N_phi + 1, 0.0));

  // Spin-Parity
  std::vector<std::complex<double>> JzP(4, 0.0);

  // 2D-loop integral
  //
  // \int d^2 \vec{qt} [...] = \int dphi \int dqt qt [...]
  //

  // Linearly discretized qt-loop, N+1!
  for (std::size_t i = 0; i < Durham::N_qt + 1; ++i) {
    const double qt = Durham::qt_MIN + i * Durham::qt_STEP;
    const double qt2 = pow2(qt);

    // Linearly discretized phi in [0,2pi), N+1!
    for (std::size_t j = 0; j < Durham::N_phi + 1; ++j) {
      const double qphi = j * Durham::phi_STEP;

      // --------------------------------------------------------------------------
      // Loop vector
      const std::vector<double> qt_ = {qt * std::cos(qphi), qt * std::sin(qphi)};

      // Fusing gluon pt-vectors
      const std::vector<double> q1 = {qt_[0] - pt1[0], qt_[1] - pt1[1]};
      const std::vector<double> q2 = {qt_[0] + pt2[0], qt_[1] + pt2[1]};

      const double q1_2 = math::vpow2(q1);
      const double q2_2 = math::vpow2(q2);

      // Get fusing gluon spin-parity (J_z^P) components
      // [q1,q2] -> [0^+,0^-,+2^+,-2^-]
      DHelicity(q1, q2, JzP);

      // ** Durham scale choise **
      double Q1_2_scale = 0.0;
      double Q2_2_scale = 0.0;
      DScaleChoise(qt2, q1_2, q2_2, Q1_2_scale, Q2_2_scale);

      // Minimum scale cutoff
      if (Q1_2_scale < Durham::qt2_MIN || Q2_2_scale < Durham::qt2_MIN) {
        continue;
      }

      // --------------------------------------------------------------------------
      // Get amplitude level pdfs
      const double fg_1 = gra::GlobalSudakovPtr->fg_xQ2M(lts.x1, Q1_2_scale, MU);
      const double fg_2 = gra::GlobalSudakovPtr->fg_xQ2M(lts.x2, Q2_2_scale, MU);

      // Amplitude weight:
      // * \pi^2 : see original KMR papers
      // *    2  : factor from initial state boson-statistics, check it:!
      // *    qt : jacobian of d^2qt -> dphi dqt qt
      std::complex<double> weight = fg_1 * fg_2 / (qt2 * q1_2 * q2_2);
      weight *= math::PIPI * 2.0 * qt;

      // Loop over (outgoing) helicity combinations.
      // Amp[h] contains initial state gluon helicity combinations --,-+,+-,++
      for (std::size_t h = 0; h < f.size(); ++h) {
        f[h][i][j] = weight * DHelProj(Amp[h], JzP);
      }

    }  // phi-loop
  }    // |q_t|-loop

  // Evaluate the total numerical integral for each helicity amplitude
  std::vector<std::complex<double>> sum(f.size(), 0.0);
  for (std::size_t h = 0; h < f.size(); ++h) {
    sum[h] = math::Simpson38Integral2D(f[h], WSimpson, Durham::qt_STEP, Durham::phi_STEP);
  }

  // *****Make sure it is of right size!*****
  lts.hamp.resize(f.size());

  // Final outcome as a coherent sum over helicity combinations
  // That is, sum at amplitude level in contrast to the usual \sum_{h \in
  // helicity set}
  // |A_h|^2
  std::complex<double> A(0, 0);
  for (std::size_t h = 0; h < f.size(); ++h) {
    lts.hamp[h] = sum[h];

    // Apply proton form factors
    lts.hamp[h] *=
        lts.excite1 ? gra::form::S3FINEL(lts.t1, lts.pfinal[1].M2()) : gra::form::S3F(lts.t1);
    lts.hamp[h] *=
        lts.excite2 ? gra::form::S3FINEL(lts.t2, lts.pfinal[2].M2()) : gra::form::S3F(lts.t2);

    // Apply phase space factors
    lts.hamp[h] *= msqrt(16.0 * math::PIPI);
    lts.hamp[h] *= msqrt(16.0 * math::PIPI);
    lts.hamp[h] *= msqrt(lts.s / lts.s_hat);

    // For the total complex amplitude (for no eikonal screening applied)
    A += lts.hamp[h];
  }
  return A;
}

// ======================================================================
// Exclusive \chi_c(0+) sub-amplitude
//
// [REFERENCE: Khoze, Martin, Ryskin, Stirling,
// https://arxiv.org/abs/hep-ph/0403218]
// [REFERENCE: Pasechnik, Szczurek, Teryaev, https://arxiv.org/abs/0709.0857v1]
//
void MDurham::Dgg2chic0(const gra::LORENTZSCALAR &lts,
                        std::vector<std::vector<std::complex<double>>> &Amp,
                        const std::vector<double> &qt1, const std::vector<double> &qt2) const {
  // @@ MUTEX LOCK @@
  gra::g_mutex.lock();
  const double alpha_s = gra::GlobalSudakovPtr->AlphaS_Q2(lts.s_hat / Durham::alphas_scale);
  gra::g_mutex.unlock();
  // @@ MUTEX LOCK @@

  const double gs2 = 4.0 * PI * alpha_s;  // coupling
  const double K_NLO = 1.68;              // NLO correction
  const double M0 = 3.41475;              // chi_c(0+) mass (GeV)
  const double W0 = 0.0108;               // chi_c(0+) width (GeV)
  const double NC = 3.0;                  // #colors

  // Gluonic width \Gamma(\chi_c(0+) -> gg), see references
  std::complex<double> A = K_NLO * 8.0 * math::zi * gs2 / M0 * msqrt(0.075) / msqrt(PI * M0 * NC);

  // Apply Breit-Wigner propagator shape delta-function
  A *= gra::form::deltaBWamp(lts.s_hat, M0, W0);

  // Initial state gluon helicity combinations
  // [RE-CHECK this point]
  std::vector<std::complex<double>> is(4, 0.0);
  is[MM] = A;
  is[MP] = A;
  is[PM] = A;
  is[PP] = A;

  // 0+ final state
  Amp[P0] = is;
}

// ======================================================================
// gg -> gg tree-level helicity amplitudes
//
//
// Basic result: d\hat{\sigma}/dt = 9/4 \pi \alpha_s^2 / E_T^4
//
// In pure gluon amplitudes: when gluon helicities are the same,
// or at most one is different from the rest, vanish for any n >= 4,
// where n is the total number of gluons (in+out)
//
// [REFERENCE: Dixon, https://arxiv.org/abs/1310.5353v1]
//
// Note that incoming gluons (a,b) are in a color singlet state,
// i.e., a = b or \delta^{ab} Tr [ .color algebra. ] has been applied
//
void MDurham::Dgg2gg(const gra::LORENTZSCALAR &lts,
                     std::vector<std::vector<std::complex<double>>> &Amp) {
  // @@ MUTEX LOCK @@
  gra::g_mutex.lock();
  const double alpha_s = gra::GlobalSudakovPtr->AlphaS_Q2(lts.s_hat / Durham::alphas_scale);
  gra::g_mutex.unlock();
  // @@ MUTEX LOCK @@

  // Vertex factor coupling, gs^2 = 4 pi alpha_s
  double norm = 4.0 * PI * alpha_s;

  // Color part
  const double NC = 3.0;
  norm *= NC / msqrt(NC * NC - 1);
  norm *= 2.0;

  // Helicity phase
  const double phi = lts.decaytree[0].p4.Phi();
  const std::complex<double> negphase = std::exp(-zi * phi);
  const std::complex<double> posphase = std::exp(zi * phi);

  const double aux01 = pow2(lts.s_hat) / (lts.u_hat * lts.t_hat);
  const double aux23 = lts.u_hat / lts.t_hat;

  // We use natural binary order / Madgraph order here
  std::vector<std::complex<double>> dualAmp(16, 0.0);

  // Loop over final state gluon pair helicity combinations
  for (std::size_t h = 0; h < fs2.size(); ++h) {
    std::vector<std::complex<double>> is(4, 0.0);
    if (fs2[h] == PP || fs2[h] == MM) {
      is[MM] = (fs2[h] == PP) ? 0 : norm * aux01;  // -- => ++ or --
      is[MP] = 0;                                  // -+ => ++ or --
      is[PM] = 0;                                  // +- => ++ or --
      is[PP] = (fs2[h] == PP) ? norm * aux01 : 0;  // ++ => ++ or --
    } else if (fs2[h] == PM || fs2[h] == MP) {
      is[MM] = 0;                                                       // -- => +- or -+
      is[MP] = negphase * norm * ((fs2[h] == PM) ? 1 / aux23 : aux23);  // -+ => +- or -+
      is[PM] = posphase * norm * ((fs2[h] == PM) ? aux23 : 1 / aux23);  // +- => +- or -+
      is[PP] = 0;                                                       // ++ => +- or -+
    }
    Amp[h] = is;

    for (std::size_t k = 0; k < 4; ++k) {
      dualAmp[4 * k + h] = is[k];
    }
  }

  // ------------------------------------------------------------------------
  const bool DEBUG = false;
  if (DEBUG) {
    // TEST amplitude squared WITH MADGRAPH (UNDER IMPLEMENTATION)
    double madsum = 0.0;
    double thissum = 0.0;
    for (unsigned int i = 0; i < 16; ++i) {
      printf("i=%2d : |mad|^2 = %0.3f , |amp|^2 = %0.3f \n", i, math::abs2(lts.hamp[i]),
             math::abs2(dualAmp[i]));
      madsum += math::abs2(lts.hamp[i]);
      thissum += math::abs2(dualAmp[i]);
    }
    Asum += thissum / madsum;
    Nsum += 1.0;
    printf("THIS/MADGRAPH = %0.10f \n\n", Asum / Nsum);
  }
  // ------------------------------------------------------------------------
}

// ======================================================================
// gg -> meson pair (gg -> q\bar{q} q\bar{q})
//
// [REFERENCE: Harland-Lang, Khoze, Ryskin, Stirling,
// https://arxiv.org/pdf/1105.1626.pdf]

// Meson wave-function (without scale Q⁻dependence)
// x the longitudinal momentum fraction of a parton within the meson
// fM the meson decay constant
//
// Normalization: \int_0^1 dx phi(x) \equiv fM/(2\sqrt{3})
//
double MDurham::phi_CZ(double x, double fM) const {
  return 5.0 * std::sqrt(3.0) * fM * x * (1.0 - x) * pow2(1.0 - 2.0 * x);
  // return 6.0*std::sqrt(3.0)*fM * x*(1.0 - x);
}

// For tabulating the meson wave function
//
// [REFERENCE: Takizawa,
// http://www-nuclth.phys.sci.osaka-u.ac.jp/jp-usw/talk/Takizawa.pdf]
//
// Here, one could perhaps update the phenomenology (check more recent papers).
//
// ----------------------------------------------------------------------
// Three charge neutral state are in the SU(3)_F quark model nonet
// (octet+singlet):
// pi0, eta8, eta0
//
// Rotation (mixing) to physical basis via:
//
// [eta > = [cos th  -sin th] [eta8>
// [eta'>   [sin th   cos th] [eta0>
//
// Decay constants:
//
// f_{8eta}  =  f_8 * cos th
// f_{8eta'} =  f_8 * sin th
//
// f_{0eta}  = -f_0 * sin th
// f_{0eta'} =  f_0 * cos th
//
std::vector<double> MDurham::EvalPhi(int N, int pdg) const {
  // Meson decay constants
  double fM = 0.0;
  static const std::vector<double> supported = {111, 211, 321, 311};  // pi0, pi+, K+, K0
  if (std::find(supported.begin(), supported.end(), pdg) != supported.end()) {
    fM = PDG::fM_meson.at(pdg);
  } else {
    fM = PDG::fM_meson.at(111);  // else, take pi0
  }
  const double fM_eta8 = fM * 1.0;  // just use pi0
  const double fM_eta0 = fM * 1.0;  // just use pi0

  // Mixing angle
  static const double THETA_eta = math::Deg2Rad(-15.4);
  const double costh = std::cos(THETA_eta);
  const double sinth = std::sin(THETA_eta);

  // Evaluate
  std::vector<double> f(N + 1);
  const double STEP = 1.0 / N;
  for (const auto &i : aux::indices(f)) {
    const double x = i * STEP;
    if (pdg == 221) {  // eta via mixing
      f[i] = costh * phi_CZ(x, fM_eta8) - sinth * phi_CZ(x, fM_eta0);
    } else if (pdg == 331) {  // etaprime via mixing
      f[i] = sinth * phi_CZ(x, fM_eta8) + costh * phi_CZ(x, fM_eta0);
    } else {
      f[i] = phi_CZ(x, fM);
    }
  }
  return f;
}

// Meson pair amplitude
//
void MDurham::Dgg2MMbar(const gra::LORENTZSCALAR &lts,
                        std::vector<std::vector<std::complex<double>>> &Amp) {
  // Wave function x-discretization
  // [Easy speed & accuracy improvement: Check lambda-functions below,
  // they could be pre-calculated and interpolated as a function of costheta]
  const int Nx = 96;
  const double STEPx = 1.0 / Nx;

  // Init 2D-Simpson weight matrix (will be calculated only once, being
  // static), C++11 handles multithreaded static initialization
  const static MMatrix<double> WSimpson = math::Simpson38Weight2D(Nx, Nx);

  const double delta_AB = 8;  // Sum over \delta^AB [gluons in with the same color]
  const double NC = 3;        // Three colors
  const double CF = (pow2(NC) - 1.0) / (NC * 2.0);  // SU(3) algebra

  // @@ MUTEX LOCK @@
  gra::g_mutex.lock();
  const double alpha_s = gra::GlobalSudakovPtr->AlphaS_Q2(lts.s_hat / Durham::alphas_scale);
  gra::g_mutex.unlock();
  // @@ MUTEX LOCK @@

  const double shat = lts.s_hat;
  const double m = lts.decaytree[0].p4.M();
  const double beta = kinematics::Beta(pow2(m), shat);
  const double costheta = (1.0 + 2.0 * lts.t_hat / lts.s_hat - 2.0 * pow2(m) / lts.s_hat) / beta;
  const double costheta2 = pow2(costheta);

  // ------------------------------------------------------------------
  // ** Hard angular cut-off **
  // some sub-amplitudes are singular when |costheta| -> 1

  if (std::abs(costheta) > Durham::MAXCOS) {
    Amp[0] = std::vector<std::complex<double>>(4, 0.0);  // Return zero
    return;
  }
  // ------------------------------------------------------------------

  // Helicity phase
  const double phi = lts.decaytree[0].p4.Phi();
  const std::complex<double> posphase = std::exp(2.0 * zi * phi);
  const std::complex<double> negphase = std::exp(-2.0 * zi * phi);

  // ------------------------------------------------------------------
  // Mesons scalar flavor octet (non-singlet) amplitude:
  // |\pi0\pi0>, |\pi+\pi->, |K+K->, |K0\bar{K0}>
  //
  // T_+- = T_-+
  auto T_SFO_PM = [&](double x, double y) {

    const double a = (1.0 - x) * (1.0 - y) + x * y;  // +
    const double b = (1.0 - x) * (1.0 - y) - x * y;  // -

    return 1.0 / (x * y * (1.0 - x) * (1.0 - y)) * (x * (1.0 - x) + y * (1.0 - y)) /
           (pow2(a) - pow2(b) * costheta2) * (NC / 2.0) * (costheta2 - 2.0 * CF / NC * a);
  };
  // ------------------------------------------------------------------

  // ------------------------------------------------------------------
  // SU(3)_F scalar flavor-singlet amplitude:
  // |\eta'\eta' >
  // |\eta\eta >, |\eta\eta' > (|\eta-\eta' > mixing)
  //
  // T_++ = T_--
  auto T_SFS_PP = [&](double x, double y) {

    return 1.0 / (x * y * (1.0 - x) * (1.0 - y)) * (1.0 + costheta2) / pow2(1.0 - costheta2);
  };

  // T_+- = T_-+
  auto T_SFS_PM = [&](double x, double y) {

    return 1.0 / (x * y * (1.0 - x) * (1.0 - y)) * (1.0 + 3.0 * costheta2) /
           (2.0 * pow2(1.0 - costheta2));
  };
  // ------------------------------------------------------------------

  // pi0, pi+, K+, K0
  static const std::vector<int> SFO_PDG = {111, 211, 321, 311};
  // eta, eta'
  static const std::vector<int> SFS_PDG = {221, 331};

  static const int pdg0 = std::abs(lts.decaytree[0].p.pdg);
  static const int pdg1 = std::abs(lts.decaytree[1].p.pdg);

  // Evaluate only once the meson wave functions
  static const std::vector<double> wfphi0 = EvalPhi(Nx, pdg0);
  static const std::vector<double> wfphi1 = EvalPhi(Nx, pdg1);

  // -------------------------------------------------------------------
  const double CUTOFF = 1e-15;  // To avoid singularity at x = 0, x = 1
  static const std::vector<double> xval = math::linspace<std::vector>(CUTOFF, 1.0 - CUTOFF, Nx + 1);

  // ------------------------------------------------------------------
  // Integral over meson wave functions:
  // M\int_0^1 dx dy \phi_M(x) \phi_\bar{M}(y) T_{\lambda\lambda'}
  // (x,y,\hat{s},\theta)

  // Normalization factor
  const double norm = (delta_AB / NC) * (64.0 * math::PIPI * pow2(alpha_s)) / shat;

  // Scalar flavor octet
  if ((std::find(SFO_PDG.begin(), SFO_PDG.end(), pdg0) != SFO_PDG.end())) {
    std::vector<std::complex<double>> is(4);

    MMatrix<double> f_PM(Nx + 1, Nx + 1, 0.0);

    for (std::size_t i = 0; i < Nx + 1; ++i) {
      for (std::size_t j = 0; j < Nx + 1; ++j) {
        f_PM[i][j] = wfphi0[i] * wfphi1[j] * T_SFO_PM(xval[i], xval[j]);
      }
    }

    is[PP] = 0.0;
    is[MM] = 0.0;
    is[PM] = norm * math::Simpson38Integral2D(f_PM, WSimpson, STEPx, STEPx) * posphase;
    is[MP] = is[PM] * negphase;

    Amp[0] = is;
  }

  // Scalar flavor singlet
  else if ((std::find(SFS_PDG.begin(), SFS_PDG.end(), pdg0) != SFS_PDG.end())) {
    std::vector<std::complex<double>> is(4);

    MMatrix<double> f_PP(Nx + 1, Nx + 1, 0.0);
    MMatrix<double> f_PM(Nx + 1, Nx + 1, 0.0);

    for (std::size_t i = 0; i < Nx + 1; ++i) {
      for (std::size_t j = 0; j < Nx + 1; ++j) {
        f_PP[i][j] = wfphi0[i] * wfphi1[j] * T_SFS_PP(xval[i], xval[j]);
        f_PM[i][j] = wfphi0[i] * wfphi1[j] * T_SFS_PM(xval[i], xval[j]);
      }
    }

    is[PP] = norm * math::Simpson38Integral2D(f_PP, WSimpson, STEPx, STEPx);
    is[MM] = is[PP];
    is[PM] = norm * math::Simpson38Integral2D(f_PM, WSimpson, STEPx, STEPx) * posphase;
    is[MP] = is[PM] * negphase;

    Amp[0] = is;
  } else {
    throw std::invalid_argument("MDurham::Dgg2MMbar: Unsupported meson: " +
                                std::to_string(lts.decaytree[0].p.pdg));
  }
}

// ======================================================================
// gg -> qqbar tree-level helicity amplitudes
//
//
void MDurham::Dgg2qqbar(const gra::LORENTZSCALAR &lts,
                        std::vector<std::vector<std::complex<double>>> &Amp) {
  throw std::invalid_argument("MDurham::DurhamQCD: qqbar amplitude in the next version");
}

}  // gra namespace ends
