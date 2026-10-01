// Exclusive heavy vector-meson photoproduction amplitude
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// [REFERENCE: S.P. Jones, A.D. Martin, M.G. Ryskin and T. Teubner, arXiv:1307.7099]
// [REFERENCE: https://arxiv.org/abs/2206.10161]
// [REFERENCE: https://arxiv.org/abs/1304.5162]
// [REFERENCE: https://arxiv.org/abs/2206.13343]
// [REFERENCE: https://arxiv.org/abs/0805.0717]
// [REFERENCE: https://arxiv.org/abs/hep-ph/0005250]

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <optional>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include <Eigen/SpecialFunctions>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MPhotoVM.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "json.hpp"

using gra::aux::indices;

using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

namespace {

// Reject nonfinite JMRT contributions before support cuts can turn them into zeros
template <typename T>
T JMRTFinite(T value) {
  if (!std::isfinite(std::real(value)) || !std::isfinite(std::imag(value))) {
    throw AmplitudeFailure("MPhotoVM: nonfinite JMRT contribution");
  }
  return value;
}

struct PhotoVMChannel {
  std::string label;
  int         vm_pdg    = 0;
  int         quark_pdg = 0;
};

struct PhotoVMFinalState {
  int fermion_index     = 0;
  int antifermion_index = 1;
  int fermion_pdg       = 0;
};

struct PhotoVMTargetChannels {
  std::array<std::complex<double>, 2>  amplitude{};
  std::optional<nuclear::PhotoCurrent> current;
  std::size_t                          count = 1;
};

struct PhotoVMProductionChannels {
  std::array<std::complex<double>, 4> amplitude{};
  std::size_t                         count = 1;
};

// Compute the supported heavy vector-meson channels
std::vector<PhotoVMChannel> PhotoVMChannels() {
  return {
      {"jpsi", 443, 4},           {"psi(2S)", 100443, 4},     {"Upsilon(1S)", 553, 5},
      {"Upsilon(2S)", 100553, 5}, {"Upsilon(3S)", 200553, 5},
  };
}

// Resolve one process channel label into vector meson parameters
PhotoVMChannel ResolvePhotoVMChannel(const std::string &channel) {
  for (const auto &entry : PhotoVMChannels()) {
    if (entry.label == channel) { return entry; }
  }
  throw std::invalid_argument("MPhotoVM::ResolvePhotoVMChannel: unsupported ygg[" + channel + "] heavy vector meson");
}

// Compute true for charged leptons supported by the direct PhotoVM syntax
bool IsPhotoVMLepton(int apdg) { return apdg == 11 || apdg == 13 || apdg == 15; }

// Compute true for one direct same-flavour PhotoVM charged-lepton pair
bool PhotoVMProcessAccepts(const std::vector<MDecayBranch> &tree) {
  if (tree.size() != 2 || !tree[0].legs.empty() || !tree[1].legs.empty()) { return false; }
  const int first_pdg = tree[0].p.pdg;
  return first_pdg == -tree[1].p.pdg && IsPhotoVMLepton(std::abs(first_pdg));
}

// Compute the fermion ordering of the initialized vector-meson final state
PhotoVMFinalState ResolvePhotoVMFinalState(const gra::LORENTZSCALAR &lts) {
  const int pdg0 = lts.decaytree[0].p.pdg;
  const int apdg = std::abs(pdg0);
  PhotoVMFinalState state;
  state.fermion_index     = (pdg0 > 0) ? 0 : 1;
  state.antifermion_index = 1 - state.fermion_index;
  state.fermion_pdg       = apdg;
  return state;
}

// Compute a vector-current coupling from the dilepton partial width
// g_Vll = sqrt(12 pi Gamma_ll/M_V)
double PhotoVMDecayCoupling(double mass, double gamma_ll) {
  if (!(mass > 0.0) || !(gamma_ll > 0.0)) { return 0.0; }
  return msqrt(12.0 * gra::math::PI * gamma_ll / mass);
}

// Compute the vector-meson timelike propagator
// D_V(q2) = 1/(q2-M_V^2+i M_V Gamma_V)
std::complex<double> PhotoVMPropagator(double q2, const gra::MParticle &vm) {
  return 1.0 / std::complex<double>(q2 - pow2(vm.mass), vm.mass * vm.width);
}

// Compute true when the amplitude uses one published JMRT gluon fit
bool JMRTUsesFittedGluon(const gra::MPhotoVMNumerics &num) {
  return num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_LO ||
         num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_NLO;
}

// Parse one explicit heavy-vector gluon-source selector
gra::MPhotoVMGluonSource ParsePhotoVMGluonSource(const std::string &source) {
  if (source == "LHAPDF_SHUVAEV") { return gra::MPhotoVMGluonSource::LHAPDF_SHUVAEV; }
  if (source == "JMRT_2013_LO") { return gra::MPhotoVMGluonSource::JMRT_2013_LO; }
  if (source == "JMRT_2013_NLO") { return gra::MPhotoVMGluonSource::JMRT_2013_NLO; }
  throw std::invalid_argument("MPhotoVM::ParsePhotoVMGluonSource: unknown gluon_source '" + source + "'");
}

// Compute the perturbative order of the selected JMRT running coupling
unsigned int JMRTRunningAlphaOrder(const gra::MPhotoVMNumerics &num) {
  return num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_NLO ? 2U : 1U;
}

// Compute the QCD beta-function coefficients for fixed active flavour count
// beta0 = 11-2nf/3, beta1 = 102-38nf/3
std::pair<double, double> JMRTBetaCoefficients(int nf) { return {11.0 - 2.0 * nf / 3.0, 102.0 - 38.0 * nf / 3.0}; }

// Evaluate the one-loop or two-loop coupling for fixed Lambda_QCD
// alpha_s = 4pi/(beta0 L)[1-beta1 ln(L)/(beta0^2 L)] at two loops
double JMRTAlphaSFromLambda(double q2, double lambda_qcd, int nf, unsigned int order) {
  if (!(q2 > pow2(lambda_qcd)) || !(lambda_qcd > 0.0) || (order != 1U && order != 2U)) { return 0.0; }
  const auto [beta0, beta1] = JMRTBetaCoefficients(nf);
  const double L            = std::log(q2 / pow2(lambda_qcd));
  double       alpha        = 4.0 * gra::math::PI / (beta0 * L);
  if (order == 2U) { alpha *= 1.0 - beta1 * std::log(L) / (pow2(beta0) * L); }
  return JMRTFinite(alpha) > 0.0 ? alpha : 0.0;
}

// Solve Lambda_QCD at fixed flavour number from one reference coupling
double JMRTLambdaForAlpha(double alpha, double q2, int nf, unsigned int order) {
  if (!(alpha > 0.0) || !(q2 > 0.0)) {
    throw std::invalid_argument("MPhotoVM::JMRTLambdaForAlpha: positive alpha and q2 are required");
  }
  double low  = 1.0e-6;
  double high = std::min(1.0, 0.5 * std::sqrt(q2));
  if (!(JMRTAlphaSFromLambda(q2, high, nf, order) > alpha)) {
    throw std::domain_error("MPhotoVM::JMRTLambdaForAlpha: failed to bracket Lambda_QCD");
  }
  for (unsigned int i = 0; i < 100; ++i) {
    const double mid = 0.5 * (low + high);
    if (JMRTAlphaSFromLambda(q2, mid, nf, order) < alpha) {
      low = mid;
    } else {
      high = mid;
    }
  }
  return 0.5 * (low + high);
}

// Build threshold-matched Lambda_QCD values for three to five flavours
std::map<int, double> JMRTAlphaLambdaThresholds(double alpha_mz, double mz, double charm_mass, double bottom_mass,
                                                unsigned int order) {
  std::map<int, double> out;
  out[5]                = JMRTLambdaForAlpha(alpha_mz, pow2(mz), 5, order);
  const double alpha_mb = JMRTAlphaSFromLambda(pow2(bottom_mass), out.at(5), 5, order);
  out[4]                = JMRTLambdaForAlpha(alpha_mb, pow2(bottom_mass), 4, order);
  const double alpha_mc = JMRTAlphaSFromLambda(pow2(charm_mass), out.at(4), 4, order);
  out[3]                = JMRTLambdaForAlpha(alpha_mc, pow2(charm_mass), 3, order);
  return out;
}

// Evaluate the threshold-matched running coupling used in the JMRT fit
double JMRTFittedAlphaS(const gra::MPhotoVMNumerics &num, double q2) {
  if (!(q2 > 0.0) || num.running_alpha_lambda_qcd.size() != 3) { return 0.0; }
  const double charm2  = pow2(num.HeavyQuarkMass(4));
  const double bottom2 = pow2(num.HeavyQuarkMass(5));
  const int    nf      = q2 >= bottom2 ? 5 : (q2 >= charm2 ? 4 : 3);
  return JMRTAlphaSFromLambda(q2, num.running_alpha_lambda_qcd.at(nf), nf, JMRTRunningAlphaOrder(num));
}

// Compute the selected published JMRT gluon-fit parameters
const gra::MPhotoVMGluonFitParameters &JMRTSelectedGluonFit(const gra::MPhotoVMNumerics &num) {
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_LO) { return num.jmrt_2013.lo; }
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_NLO) { return num.jmrt_2013.nlo; }
  throw std::invalid_argument("MPhotoVM::JMRTSelectedGluonFit: fitted gluon source is not active");
}

// Evaluate the fitted diagonal integrated gluon xg from JMRT Eqs. (16) and (17)
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTFittedXG(const gra::MPhotoVMNumerics &num, double x, double q2) {
  if (!(x > 0.0 && x < 1.0) || !(q2 > 0.0)) { return 0.0; }
  const auto  &fit  = JMRTSelectedGluonFit(num);
  const double logx = std::log(1.0 / x);
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_LO) {
    const double lambda = fit.a + fit.b * std::log(q2 / num.jmrt_2013.lo_scale2);
    const double out    = fit.normalization * std::exp(lambda * logx);
    return JMRTFinite(out);
  }

  const double lambda2 = pow2(num.jmrt_2013.nlo_lambda_qcd);
  const double q02     = num.jmrt_2013.nlo_q02;
  if (q2 < q02 || !(q02 > lambda2)) { return 0.0; }
  const double G          = std::log(q2 / lambda2) / std::log(q02 / lambda2);
  const double log_G      = std::max(0.0, std::log(G));
  const double double_log = std::sqrt(16.0 * 3.0 / 9.0 * logx * log_G);
  const double out        = fit.normalization * std::exp(fit.a * logx) * std::pow(q2, fit.b) * std::exp(double_log);
  return JMRTFinite(out);
}

// Compute the local small-x power of a fitted JMRT diagonal gluon
double JMRTFittedGluonLambda(const gra::MPhotoVMNumerics &num, double x, double q2) {
  if (!(x > 0.0 && x < 1.0) || !(q2 > 0.0)) { return 0.0; }
  const auto &fit = JMRTSelectedGluonFit(num);
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_LO) {
    return fit.a + fit.b * std::log(q2 / num.jmrt_2013.lo_scale2);
  }

  const double lambda2 = pow2(num.jmrt_2013.nlo_lambda_qcd);
  const double q02     = num.jmrt_2013.nlo_q02;
  if (q2 < q02 || !(q02 > lambda2)) { return 0.0; }
  const double G     = std::log(q2 / lambda2) / std::log(q02 / lambda2);
  const double log_G = std::max(0.0, std::log(G));
  const double logx  = std::log(1.0 / x);
  return fit.a + 0.5 * std::sqrt(16.0 * 3.0 / 9.0 * log_G / logx);
}

// Compute the analytic small-x skewness factor of JMRT Eq. (4)
// R_g = 2^(2lambda+3) Gamma(lambda+5/2)/[sqrt(pi) Gamma(lambda+4)]
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTSkewnessFactor(double lambda) {
  const double out = std::pow(2.0, 2.0 * lambda + 3.0) * std::tgamma(lambda + 2.5) /
                     (std::sqrt(gra::math::PI) * std::tgamma(lambda + 4.0));
  return JMRTFinite(out) > 0.0 ? out : 0.0;
}

// Compute the square-root Sudakov factor of JMRT Eq. (7)
// sqrt(T) = exp[-3 alpha_s ln^2(mu2/q2)/(8 pi)]
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTFittedSqrtSudakov(const gra::MPhotoVMNumerics &num, double q2, double mu2) {
  if (!(q2 > 0.0) || !(mu2 > 0.0)) { return 0.0; }
  if (q2 >= mu2) { return 1.0; }
  const double alpha     = JMRTFittedAlphaS(num, mu2);
  const double log_ratio = std::log(mu2 / q2);
  const double out       = std::exp(-3.0 * alpha * pow2(log_ratio) / (8.0 * gra::math::PI));
  return JMRTFinite(out);
}

// Compute the skewed integrated gluon differentiated in the last-step integral
// [REFERENCE: doi:10.1140/epjc/s10052-015-3832-8, Eq. (3); arXiv:1307.7099, Eq. (4)]
double JMRTFittedIntegratedFlux(const gra::MPhotoVMNumerics &num, double x, double q2, double mu2) {
  return JMRTSkewnessFactor(JMRTFittedGluonLambda(num, x, q2)) *
         JMRTFittedXG(num, x, q2) * JMRTFittedSqrtSudakov(num, q2, mu2);
}

// Compute d[R_g xg sqrt(T)]/d ln(q2), including the scale dependence of skewness
double JMRTFittedFluxDerivative(const gra::MPhotoVMNumerics &num, double x, double q2, double mu2) {
  const double L0 = std::log(num.jmrt_2013.nlo_q02 / pow2(num.jmrt_2013.nlo_lambda_qcd));
  const double l  = std::log(q2 / num.jmrt_2013.nlo_q02);
  const double derivative =
      num.jmrt_2013.nlo.b + std::sqrt(16.0 * 3.0 / 9.0 * std::log(1.0 / x) / std::log1p(l / L0)) / (2.0 * (L0 + l)) +
      3.0 * JMRTFittedAlphaS(num, mu2) * std::max(0.0, std::log(mu2 / q2)) / (4.0 * gra::math::PI);
  const double lambda = JMRTFittedGluonLambda(num, x, q2);
  const double lambda_derivative = (lambda - num.jmrt_2013.nlo.a) / (2.0 * (L0 + l) * std::log1p(l / L0));
  const double skew_derivative = 2.0 * std::log(2.0) + Eigen::numext::digamma(lambda + 2.5) -
                                 Eigen::numext::digamma(lambda + 4.0);
  return JMRTFittedIntegratedFlux(num, x, q2, mu2) * (derivative + skew_derivative * lambda_derivative);
}

// Compute the common hard mass squared for one quarkonium family
double JMRTHardMass2(const gra::MPhotoVMNumerics &num, int pdg) {
  return 4.0 * pow2(num.HeavyQuarkMass(num.Channel(pdg).quark_pdg));
}

// Compute the JMRT photoproduction scaling variable
// x = M_hard^2/W^2
double JMRTX(double hard_mass2, double w2) { return (w2 > 0.0) ? hard_mass2 / w2 : 0.0; }

// Compute the JMRT hard scale Qbar^2 for photoproduction
double JMRTQbar2(const gra::MPhotoVMNumerics &num, int pdg) {
  return pow2(num.HeavyQuarkMass(num.Channel(pdg).quark_pdg));
}

// Compute an integral over the positive logarithmic distance below the matching
// scale
double InfraredLogIntegral(unsigned int n, const std::function<double(double)> &f) {
  if (n == 0) { return 0.0; }
  const auto rule = gra::math::GaussLegendreRule(n, 0.0, 1.0);
  double     out  = 0.0;
  for (const auto &i : indices(rule.first)) {
    const double one_minus_t = 1.0 - rule.first[i];
    const double s           = rule.first[i] / one_minus_t;
    out += rule.second[i] * f(s) / pow2(one_minus_t);
  }
  return out;
}

// Compute the active upper k^2 bound of the JMRT last-step integral
double JMRTK2Max(const gra::MPhotoVMNumerics &num, double m2, double w2) {
  const double kinematic = 0.25 * std::max(0.0, w2 - m2);
  if (JMRTUsesFittedGluon(num)) { return kinematic; }
  return std::min(num.jmrt_k2_max, kinematic);
}

// Compute the JMRT scale mu for the last emitted parton
double JMRTMu(double k2, double qbar2) { return msqrt(std::max(k2, qbar2)); }

// Compute the coefficient-function coupling for the selected gluon source
double JMRTAlphaS(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double q2) {
  if (JMRTUsesFittedGluon(num)) { return JMRTFittedAlphaS(num, q2); }
  if (lts.GlobalSudakovPtr == nullptr) {
    throw std::invalid_argument("MPhotoVM::JMRTAlphaS: Sudakov tables are not initialized");
  }
  return lts.GlobalSudakovPtr->AlphaS_Q2(q2);
}

// Compute the perturbative JMRT last-step integral
double JMRTLastStepPerturbative(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2,
                                double k2_max) {
  if (lts.GlobalSudakovPtr == nullptr) {
    throw std::invalid_argument(
        "MPhotoVM::JMRTLastStepPerturbative: Sudakov "
        "tables are not initialized");
  }
  const double q02 = lts.GlobalSudakovPtr->GetQ2Min();
  if (!(x > 0.0 && x < 1.0) || !(qbar2 > 0.0) || !(k2_max > q02) || num.N_k == 0) { return 0.0; }

  return gra::math::LogMeasureGaussIntegral(num.N_k, q02, k2_max, [&](double k2) {
    const double mu = std::max(JMRTMu(k2, qbar2), lts.GlobalSudakovPtr->GetMuMin());
    try {
      const double alpha_fg = lts.GlobalSudakovPtr->AlphaSFlux_xQ2Mu(x, k2, mu, std::max(k2, qbar2));
      return JMRTFinite(alpha_fg / (qbar2 * (qbar2 + k2)));
    } catch (const std::exception &e) { throw AmplitudeFailure(e.what()); }
  });
}

// Compute the unresolved JMRT contribution from the selected matched infrared
// gluon
double JMRTLastStepInfrared(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2) {
  if (lts.GlobalSudakovPtr == nullptr) {
    throw std::invalid_argument("MPhotoVM::JMRTLastStepInfrared: Sudakov tables are not initialized");
  }
  const double q02 = lts.GlobalSudakovPtr->GetQ2Min();
  if (!(x > 0.0 && x < 1.0) || !(qbar2 > 0.0) || !(q02 > 0.0)) { return 0.0; }

  if (lts.GlobalSudakovPtr->GetIRMode() == gra::MSudakovIRMode::PERTURBATIVE_ONLY) { return 0.0; }
  const double mu_ir = std::max(JMRTMu(q02, qbar2), lts.GlobalSudakovPtr->GetMuMin());
  try {
    const double out = InfraredLogIntegral(num.N_k, [&](double s) {
      const double k2       = q02 * std::exp(-s);
      const double alpha_fg = lts.GlobalSudakovPtr->AlphaSFlux_xQ2Mu(x, k2, mu_ir, std::max(k2, qbar2));
      return JMRTFinite(alpha_fg / (qbar2 * (qbar2 + k2)));
    });
    return JMRTFinite(out);
  } catch (const std::exception &e) { throw AmplitudeFailure(e.what()); }
}

// Compute the fitted NLO last-step integral of JMRT Eq. (11)
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTFittedLastStepPerturbative(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x,
                                      double qbar2, double k2_max) {
  const double q02 = num.jmrt_2013.nlo_q02;
  if (!(x > 0.0 && x < 1.0) || !(qbar2 > 0.0) || !(k2_max > q02) || num.N_k == 0) { return 0.0; }
  // log(k2/Q0^2) = u^2 removes the integrable square-root singularity at Q0
  const auto [nodes, weights] = gra::math::GaussLegendreRule(num.N_k, 0.0, std::sqrt(std::log(k2_max / q02)));
  double integral             = 0.0;
  for (const auto &i : indices(nodes)) {
    const double k2  = q02 * std::exp(pow2(nodes[i]));
    const double mu2 = std::max(k2, qbar2);
    integral += weights[i] * 2.0 * nodes[i] * JMRTAlphaS(lts, num, mu2) * JMRTFittedFluxDerivative(num, x, k2, mu2) /
                (qbar2 * (qbar2 + k2));
  }
  return JMRTFinite(integral);
}

// Compute the strict linear infrared term of JMRT Eqs. (10) and (11)
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTFittedLastStepInfrared(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x,
                                  double qbar2) {
  const double q02 = num.jmrt_2013.nlo_q02;
  if (!(x > 0.0 && x < 1.0) || !(qbar2 > 0.0) || !(q02 > 0.0)) { return 0.0; }
  const double mu2      = std::max(q02, qbar2);
  const double boundary = JMRTFittedIntegratedFlux(num, x, q02, mu2);
  const double alpha    = JMRTAlphaS(lts, num, mu2);
  const double out      = std::log1p(q02 / qbar2) * alpha * boundary / (qbar2 * q02);
  return JMRTFinite(out);
}

// Compute the collinear LO fitted-gluon kernel of JMRT Eq. (1)
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTFittedLOKernel(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2) {
  const double xg    = JMRTFittedXG(num, x, qbar2);
  const double alpha = JMRTAlphaS(lts, num, qbar2);
  const double skew  = JMRTSkewnessFactor(JMRTFittedGluonLambda(num, x, qbar2));
  const double out   = skew * alpha * xg / pow2(qbar2);
  return JMRTFinite(out);
}

// Compute the fitted NLO kernel including skewness and the Eq. (11) infrared
// term
double JMRTFittedNLOKernel(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2,
                           double k2_max) {
  const double perturbative = JMRTFittedLastStepPerturbative(lts, num, x, qbar2, k2_max);
  const double infrared     = JMRTFittedLastStepInfrared(lts, num, x, qbar2);
  const double out          = perturbative + infrared;
  return JMRTFinite(out);
}

// Compute the full JMRT imaginary amplitude kernel
double JMRTImaginaryKernel(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2,
                           double k2_max) {
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_LO) { return JMRTFittedLOKernel(lts, num, x, qbar2); }
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_NLO) {
    return JMRTFittedNLOKernel(lts, num, x, qbar2, k2_max);
  }
  const double pert = JMRTLastStepPerturbative(lts, num, x, qbar2, k2_max);
  const double ir   = JMRTLastStepInfrared(lts, num, x, qbar2);
  const double out  = pert + ir;
  return JMRTFinite(out);
}

// Compute the local small-x power of the JMRT amplitude
double JMRTLambda(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, double x, double qbar2,
                  double k2_max) {
  if (!(x > 0.0 && x < 1.0)) { return 0.0; }

  const double h       = num.jmrt_lambda_step;
  const double x_plus  = x * std::exp(-h);
  const double x_minus = x * std::exp(h);
  if (!(x_plus > 0.0 && x_plus < 1.0 && x_minus > 0.0 && x_minus < 1.0)) { return 0.0; }

  const double a_plus  = std::abs(JMRTImaginaryKernel(lts, num, x_plus, qbar2, k2_max));
  const double a_minus = std::abs(JMRTImaginaryKernel(lts, num, x_minus, qbar2, k2_max));
  if (!(a_plus > 0.0 && a_minus > 0.0)) { return 0.0; }

  const double lambda = (std::log(a_plus) - std::log(a_minus)) / (2.0 * h);
  return JMRTFinite(lambda);
}

// Compute the JMRT derivative-dispersion real-to-imaginary ratio
// Re(A)/Im(A) = pi lambda/2
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099]
double JMRTRealPartRatio(double lambda) {
  const double rho = 0.5 * gra::math::PI * lambda;
  return JMRTFinite(rho);
}

// Compute the real-part correction as a complex amplitude factor
std::complex<double> JMRTRealPartFactor(double lambda) { return std::complex<double>(JMRTRealPartRatio(lambda), 1.0); }

// Compute the state-dependent diffractive momentum-transfer slope
// B(W) = B0+4 alpha' ln(W/W0)
double JMRTSlopeB(const gra::MPhotoVMNumerics &num, int pdg, double w2) {
  const double W          = msqrt(std::max(0.0, w2));
  const auto  &parameters = num.Channel(pdg);
  const double log_arg    = std::max(W / parameters.cross_section_slope_W0, 1e-12);
  return std::max(
      0.0, parameters.cross_section_slope_B0 + 4.0 * parameters.cross_section_slope_alpha_prime * std::log(log_arg));
}

// Derive the produced-vector nucleon profile from the forward gamma-p amplitude
// [REFERENCE: Bauer et al., Rev. Mod. Phys. 50 (1978) 261]
flux::PhotoTargetProfile VectorNucleonProfile(const MPhotoVMNumerics &num, int pdg, const MParticle &vm,
                                              const double w2, const std::complex<double> forward) {
  flux::PhotoTargetProfile profile;
  profile.slope                      = JMRTSlopeB(num, pdg, w2);
  const std::complex<double> reduced = forward / w2;
  if (std::abs(reduced.imag()) > 0.0) { profile.eta = reduced.real() / reduced.imag(); }
  const double width     = num.LeptonicWidthEE(vm.pdg);
  const double fv_over_e = std::sqrt(qed::alpha_QED() * vm.mass / (3.0 * width));
  profile.sigma_eff      = fv_over_e * std::abs(reduced.imag()) * PDG::GeV2mb;
  profile.x              = JMRTX(JMRTHardMass2(num, pdg), w2);
  profile.scale2         = JMRTQbar2(num, pdg);
  if (!std::isfinite(profile.sigma_eff) || profile.sigma_eff < 0.0) {
    throw AmplitudeFailure("MPhotoVM: invalid derived vector-nucleon profile");
  }
  return profile;
}

// Compute the signed forward JMRT amplitude with d sigma/dt = |A/W^2|^2/(16 pi)
// [REFERENCE: S.P. Jones et al., arXiv:1307.7099, Eqs. (1), (4), (5) and (11)]
std::complex<double> JMRTProductionAmplitude(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num, int pdg,
                                             const gra::MParticle &vm, double w2) {
  const double hard_mass2 = JMRTHardMass2(num, pdg);
  const double x          = JMRTX(hard_mass2, w2);
  const double qbar2      = JMRTQbar2(num, pdg);
  const double k2_max     = JMRTK2Max(num, hard_mass2, w2);
  if (!(x > 0.0 && x < 1.0)) { return 0.0; }
  if (num.gluon_source == gra::MPhotoVMGluonSource::LHAPDF_SHUVAEV && !(k2_max > lts.GlobalSudakovPtr->GetQ2Min())) {
    return 0.0;
  }
  if (num.gluon_source == gra::MPhotoVMGluonSource::JMRT_2013_NLO && !(k2_max > num.jmrt_2013.nlo_q02)) { return 0.0; }

  const double kernel = JMRTImaginaryKernel(lts, num, x, qbar2, k2_max);
  const double lambda = JMRTLambda(lts, num, x, qbar2, k2_max);
  const double norm   = std::sqrt(num.LeptonicWidthEE(pdg) * gra::math::pow3(vm.mass) * gra::math::pow4(gra::math::PI) /
                                  (3.0 * gra::qed::alpha_QED()));
  return JMRTFinite(w2 * norm * kernel * JMRTRealPartFactor(lambda) * num.Channel(pdg).wavefunction_correction);
}

// Compute the helicity-independent JMRT production amplitude for one photon
// direction
PhotoVMTargetChannels DirectionTargetAmplitudes(const gra::LORENTZSCALAR &lts, const gra::MPhotoVMNumerics &num,
                                                int pdg, const gra::MParticle &vm, bool photon_from_upper,
                                                const SoftModel &soft_model, const flux::PhotoDissParam *diss,
                                                const std::complex<double> diss_forward) {
  const gra::M4Vec     &q_photon = photon_from_upper ? lts.q1 : lts.q2;
  const ForwardLegState target_state =
      ResolveForwardLegState(lts, photon_from_upper ? ForwardBeamLeg::Lower : ForwardBeamLeg::Upper);
  PhotoVMTargetChannels out;

  if (!gra::flux::SupportsPhotoTarget(target_state)) { return out; }

  const gra::M4Vec p_target     = gra::flux::PhotoTargetMomentum(target_state);
  const double     w2           = (q_photon + p_target).M2();
  const double     target_mass2 = target_state.IsExcited() ? target_state.mass2 : p_target.M2();
  if (!(target_mass2 > 0.0) || !(w2 > pow2(msqrt(MPhotoQCD::CentralMass2(lts)) + msqrt(target_mass2)))) {
    auto target = gra::flux::ZeroPhotoTarget(target_state);
    out.count   = target.factor.size();
    out.current = std::move(target.current);
    return out;
  }

  const std::complex<double> forward = target_state.IsExcited() && diss != nullptr
      ? w2 * diss_forward : JMRTProductionAmplitude(lts, num, pdg, vm, w2);
  const double               elastic = gra::form::ExpSlopeAmplitude(JMRTSlopeB(num, pdg, w2), target_state.t);
  double                     hadron  = elastic;
  // Replace the elastic JMRT t slope by the inclusive transition on an excited
  // target leg
  if (target_state.IsExcited()) {
    hadron = diss != nullptr ? flux::PhotoDissFactor(*diss, w2, target_state.t, target_state.mass2)
        : soft_model.ForwardExcitationFactor(soft_model.ForwardExcitationExchange(), target_state.t, target_state.mass2);
  }
  auto target =
      gra::flux::ResolvePhotoTarget(target_state, VectorNucleonProfile(num, pdg, vm, w2, forward), elastic, hadron);
  out.current = std::move(target.current);
  out.count   = std::min(out.amplitude.size(), target.factor.size());
  for (std::size_t i = 0; i < out.count; ++i) {
    out.amplitude[i] = JMRTFinite(forward * target.factor[i]);
  }
  return out;
}

// Couple explicit emitter sectors to all target-sector amplitudes
PhotoVMProductionChannels DirectionProductionAmplitudes(const gra::LORENTZSCALAR    &lts,
                                                        const PhotoVMTargetChannels &target, const int photon_leg,
                                                        const int lambda) {
  const ForwardLegState photon_state =
      ResolveForwardLegState(lts, photon_leg == 1 ? ForwardBeamLeg::Upper : ForwardBeamLeg::Lower);
  const auto                source = gra::flux::PhotoSourceAmplitudes(lts, photon_state, lambda);
  PhotoVMProductionChannels out;
  out.count           = source.size() * target.count;
  std::size_t channel = 0;
  for (const std::complex<double> photon : source) {
    for (std::size_t transition = 0; transition < target.count; ++transition) {
      const std::complex<double> amplitude = photon * target.amplitude[transition];
      out.amplitude[channel++] = JMRTFinite(amplitude);
    }
  }
  return out;
}

// Compute all PhotoVM helicity amplitudes for one validated final state
std::vector<std::complex<double>> BuildPhotoVMHelicityAmplitudes(
    gra::LORENTZSCALAR &lts, const PhotoVMFinalState &state, int pdg, const gra::MPhotoVMNumerics &num,
    const gra::MParticle &vm, double q2, const gra::MDirac &dirac, const SoftModel &soft_model,
    const flux::PhotoDissParam *diss, const std::complex<double> diss_forward) {
  const gra::MDecayBranch   &fermion             = lts.decaytree[state.fermion_index];
  const gra::MDecayBranch   &anti                = lts.decaytree[state.antifermion_index];
  const gra::M4Vec           boson               = fermion.p4 + anti.p4;
  const ForwardLegState      upper_state         = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState      lower_state         = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const bool                 separate_directions = gra::flux::SplitPhotoDirections(lts, upper_state, lower_state);
  const std::complex<double> prop                = PhotoVMPropagator(q2, vm);
  const double               gll                 = PhotoVMDecayCoupling(vm.mass, num.LeptonicWidthEE(pdg));
  PhotoVMTargetChannels      upper_target        = DirectionTargetAmplitudes(lts, num, pdg, vm, true, soft_model, diss, diss_forward);
  PhotoVMTargetChannels      lower_target        = DirectionTargetAmplitudes(lts, num, pdg, vm, false, soft_model, diss, diss_forward);
  lts.screening.photo.current[1]                           = std::move(upper_target.current);
  lts.screening.photo.current[0]                           = std::move(lower_target.current);
  const std::size_t upper_count                  = gra::flux::PhotoSourceCount(upper_state) * upper_target.count;
  const std::size_t lower_count                  = gra::flux::PhotoSourceCount(lower_state) * lower_target.count;
  const std::size_t channel_count                = separate_directions ? upper_count + lower_count : 1;

  std::vector<std::complex<double>> amplitudes;
  amplitudes.assign(4 * channel_count, 0.0);
  const auto polarization = MPhotoQCD::TransverseStates(lts, boson);

  constexpr auto helicities = spin::BinaryHelicityLabelsX2();
  for (const int lambda : helicities) {
    const PhotoVMProductionChannels                 upper = DirectionProductionAmplitudes(lts, upper_target, 1, lambda);
    const PhotoVMProductionChannels                 lower = DirectionProductionAmplitudes(lts, lower_target, 2, lambda);
    const FTensor::Tensor1<std::complex<double>, 4> eps   = polarization[spin::BinaryHelicityIndexX2(lambda)];

    std::size_t row = 0;
    for (const int hf : helicities) {
      for (const int ha : helicities) {
        const MDirac::Current      current = gra::qed::FFVCurrent(dirac, fermion.p4, anti.p4, hf, ha, gll, gll);
        const std::complex<double> decay   = gra::MPhotoQCD::ContractCurrentPolarization(current, eps);
        if (!separate_directions) {
          amplitudes[row++] += (upper.amplitude[0] + lower.amplitude[0]) * prop * decay;
        } else {
          // Store emitter sectors outside target sectors in each direction
          for (std::size_t i = 0; i < upper.count; ++i) { amplitudes[row++] += upper.amplitude[i] * prop * decay; }
          for (std::size_t i = 0; i < lower.count; ++i) { amplitudes[row++] += lower.amplitude[i] * prop * decay; }
        }
      }
    }
  }
  return gra::flux::CompletePhotoInitialSpinStates(lts, amplitudes);
}

// Validate the published JMRT fit and running-coupling constants
void ValidateJMRT2013Parameters(const gra::MPhotoVMNumerics &num) {
  const auto valid_fit = [](const gra::MPhotoVMGluonFitParameters &fit) {
    return std::isfinite(fit.normalization) && fit.normalization > 0.0 && std::isfinite(fit.a) && std::isfinite(fit.b);
  };
  if (!valid_fit(num.jmrt_2013.lo) || !valid_fit(num.jmrt_2013.nlo)) {
    throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid jmrt_2013 gluon fit");
  }
  if (!std::isfinite(num.jmrt_2013.lo_scale2) || num.jmrt_2013.lo_scale2 <= 0.0 ||
      !std::isfinite(num.jmrt_2013.nlo_q02) || num.jmrt_2013.nlo_q02 <= 0.0 ||
      !std::isfinite(num.jmrt_2013.nlo_lambda_qcd) || num.jmrt_2013.nlo_lambda_qcd <= 0.0 ||
      num.jmrt_2013.nlo_q02 <= pow2(num.jmrt_2013.nlo_lambda_qcd) || !std::isfinite(num.jmrt_2013.alpha_s_mz) ||
      num.jmrt_2013.alpha_s_mz <= 0.0 || num.jmrt_2013.alpha_s_mz >= 1.0 || !std::isfinite(num.jmrt_2013.mz) ||
      num.jmrt_2013.mz <= 0.0) {
    throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid jmrt_2013 scale or coupling");
  }
  if (JMRTUsesFittedGluon(num)) {
    for (const int nf : {3, 4, 5}) {
      const auto it = num.running_alpha_lambda_qcd.find(nf);
      if (it == num.running_alpha_lambda_qcd.end() || !std::isfinite(it->second) || it->second <= 0.0) {
        throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid threshold-matched Lambda_QCD");
      }
    }
  }
}

}  // namespace

// Compute the heavy-quark mass for a vector-meson quark PDG id
double MPhotoVMNumerics::HeavyQuarkMass(int quark_pdg) const {
  const auto it = heavy_quark_mass.find(std::abs(quark_pdg));
  if (it == heavy_quark_mass.end()) {
    throw std::invalid_argument("MPhotoVMNumerics::HeavyQuarkMass: missing mass for heavy quark " +
                                std::to_string(quark_pdg));
  }
  return it->second;
}

// Compute the e+e- partial width used for vector-meson normalization
double MPhotoVMNumerics::LeptonicWidthEE(int vm_pdg) const { return Channel(vm_pdg).leptonic_width_ee; }

// Compute the configured physics parameters for one vector meson
const MPhotoVMChannelParameters &MPhotoVMNumerics::Channel(int vm_pdg) const {
  const auto it = channels.find(std::abs(vm_pdg));
  if (it == channels.end()) {
    throw std::invalid_argument("MPhotoVMNumerics::Channel: missing parameters for PDG " + std::to_string(vm_pdg));
  }
  return it->second;
}

// Validate vector-meson numerical steering
void MPhotoVMNumerics::Validate() const {
  if (!std::isfinite(jmrt_k2_max) || jmrt_k2_max <= 0.0) {
    throw std::invalid_argument("MPhotoVMNumerics::Validate: jmrt_k2_max must be positive");
  }
  if (!std::isfinite(jmrt_lambda_step) || jmrt_lambda_step <= 0.0 || jmrt_lambda_step >= 1.0) {
    throw std::invalid_argument("MPhotoVMNumerics::Validate: jmrt_lambda_step must be in (0,1)");
  }
  if (N_k < 4) { throw std::invalid_argument("MPhotoVMNumerics::Validate: N_k must be at least 4"); }
  ValidateJMRT2013Parameters(*this);
  if (heavy_quark_mass.empty()) {
    throw std::invalid_argument("MPhotoVMNumerics::Validate: heavy_quark_mass cannot be empty");
  }
  for (const auto &[quark_pdg, mass] : heavy_quark_mass) {
    if ((quark_pdg != 4 && quark_pdg != 5) || !std::isfinite(mass) || mass <= 0.0) {
      throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid heavy_quark_mass entry");
    }
  }
  if (channels.empty()) { throw std::invalid_argument("MPhotoVMNumerics::Validate: channels cannot be empty"); }
  for (const auto &[vm_pdg, parameters] : channels) {
    if (vm_pdg <= 0 || (std::abs(parameters.quark_pdg) != 4 && std::abs(parameters.quark_pdg) != 5) ||
        parameters.ground_state_pdg <= 0 || !std::isfinite(parameters.wavefunction_correction) ||
        parameters.wavefunction_correction <= 0.0 || !std::isfinite(parameters.leptonic_width_ee) ||
        parameters.leptonic_width_ee <= 0.0 || !std::isfinite(parameters.cross_section_slope_B0) ||
        parameters.cross_section_slope_B0 <= 0.0 || !std::isfinite(parameters.cross_section_slope_alpha_prime) ||
        parameters.cross_section_slope_alpha_prime < 0.0 || !std::isfinite(parameters.cross_section_slope_W0) ||
        parameters.cross_section_slope_W0 <= 0.0) {
      throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid channel entry for PDG " +
                                  std::to_string(vm_pdg));
    }
    (void)HeavyQuarkMass(parameters.quark_pdg);
  }
  for (const auto &entry : PhotoVMChannels()) {
    const auto &parameters = Channel(entry.vm_pdg);
    if (std::abs(parameters.quark_pdg) != std::abs(entry.quark_pdg)) {
      throw std::invalid_argument("MPhotoVMNumerics::Validate: channel quark mismatch for PDG " +
                                  std::to_string(entry.vm_pdg));
    }
    const auto &ground_state = Channel(parameters.ground_state_pdg);
    if (std::abs(ground_state.quark_pdg) != std::abs(parameters.quark_pdg) ||
        ground_state.ground_state_pdg != parameters.ground_state_pdg) {
      throw std::invalid_argument("MPhotoVMNumerics::Validate: invalid ground state for PDG " +
                                  std::to_string(entry.vm_pdg));
    }
  }
}

// Configure vector meson physics and numerics from immutable JSON text
void MPhotoVMNumerics::ConfigureFromJson(const std::string &general_file, const std::string &general_json,
                                         const std::string &numerics_file, const std::string &numerics_json) {
  using json = nlohmann::json;

  try {
    Configure(json::parse(general_json), general_file, json::parse(numerics_json), numerics_file);
  } catch (const json::exception &e) {
    std::string str = "MPhotoVMNumerics::ReadParameters: Error parsing " + general_file + " or " + numerics_file +
                      " (Check PARAM_PHOTOVM and NUMERICS_PHOTOVM): " + std::string(e.what());
    throw std::invalid_argument(str);
  }
}

// Configure vector meson physics and numerics from parsed card documents
void MPhotoVMNumerics::Configure(const nlohmann::json &general, const std::string &general_file,
                                 const nlohmann::json &numerics, const std::string &numerics_file) {
  using json = nlohmann::json;

  try {
    dissociation = flux::ReadPhotoDiss(general.at("PARAM_REGGE").at("photoprod_diss"));
    const auto &physics     = general.at("PARAM_PHOTOVM");
    const auto &integration = numerics.at("NUMERICS_PHOTOVM");
    jmrt_k2_max             = integration.at("jmrt_k2_max");
    jmrt_lambda_step        = integration.at("jmrt_lambda_step");
    const auto &count = integration.at("N_k");
    if (!count.is_number_integer() || count.get<long double>() < 4.0L ||
        count.get<long double>() > std::numeric_limits<unsigned int>::max()) {
      throw std::invalid_argument("MPhotoVMNumerics::Configure: N_k must be an integer in the unsigned range, at least four");
    }
    N_k = count.get<unsigned int>();
    gluon_source            = ParsePhotoVMGluonSource(physics.at("gluon_source").get<std::string>());

    const auto &jmrt     = physics.at("jmrt_2013");
    const auto  read_fit = [](const json &row, const std::string &name) {
      if (!row.is_array() || row.size() != 3) {
        throw std::invalid_argument("MPhotoVMNumerics::ReadParameters: " + name + " needs [N,a,b]");
      }
      return MPhotoVMGluonFitParameters{row.at(0).get<double>(), row.at(1).get<double>(), row.at(2).get<double>()};
    };
    jmrt_2013.lo             = read_fit(jmrt.at("lo_fit"), "jmrt_2013.lo_fit");
    jmrt_2013.nlo            = read_fit(jmrt.at("nlo_fit"), "jmrt_2013.nlo_fit");
    jmrt_2013.lo_scale2      = jmrt.at("lo_scale2");
    jmrt_2013.nlo_q02        = jmrt.at("nlo_q02");
    jmrt_2013.nlo_lambda_qcd = jmrt.at("nlo_lambda_qcd");
    jmrt_2013.alpha_s_mz     = jmrt.at("alpha_s_mz");
    jmrt_2013.mz             = jmrt.at("mz");

    heavy_quark_mass.clear();
    for (auto it = physics.at("heavy_quark_mass").begin(); it != physics.at("heavy_quark_mass").end(); ++it) {
      const int  quark_pdg = std::abs(std::stoi(it.key()));
      const auto inserted  = heavy_quark_mass.emplace(quark_pdg, it.value().get<double>()).second;
      if (!inserted) {
        throw std::invalid_argument("MPhotoVMNumerics::ReadParameters: duplicate heavy-quark PDG " +
                                    std::to_string(quark_pdg));
      }
    }
    running_alpha_lambda_qcd.clear();
    if (JMRTUsesFittedGluon(*this)) {
      running_alpha_lambda_qcd = JMRTAlphaLambdaThresholds(jmrt_2013.alpha_s_mz, jmrt_2013.mz, HeavyQuarkMass(4),
                                                           HeavyQuarkMass(5), JMRTRunningAlphaOrder(*this));
    }

    channels.clear();
    for (const auto &row : physics.at("channels")) {
      if (!row.is_array() || row.size() != 8) {
        throw std::invalid_argument(
            "MPhotoVMNumerics::ReadParameters: PARAM_PHOTOVM.channels rows "
            "need 8 values");
      }
      const int                       vm_pdg     = std::abs(row.at(0).get<int>());
      const MPhotoVMChannelParameters parameters = {
          row.at(1).get<int>(),    std::abs(row.at(2).get<int>()), row.at(3).get<double>(), row.at(4).get<double>(),
          row.at(5).get<double>(), row.at(6).get<double>(),        row.at(7).get<double>()};
      if (!channels.emplace(vm_pdg, parameters).second) {
        throw std::invalid_argument("MPhotoVMNumerics::ReadParameters: duplicate vector-meson PDG " +
                                    std::to_string(vm_pdg));
      }
    }
    Validate();
  } catch (const json::exception &e) {
    std::string str = "MPhotoVMNumerics::ReadParameters: Error parsing " + general_file + " or " + numerics_file +
                      " (Check PARAM_PHOTOVM and NUMERICS_PHOTOVM): " + std::string(e.what());
    throw std::invalid_argument(str);
  }

  initialized = true;
}

// Construct one immutable heavy-vector photoproduction block
MPhotoVMNumericsPtr ReadPhotoVMNumerics(const MModelTune &tune) {
  auto numerics = std::make_shared<MPhotoVMNumerics>();
  numerics->Configure(tune.General(), tune.GeneralFile(), tune.Numerics(), tune.NumericsFile());
  return numerics;
}

// Compute the run owned heavy-vector photoproduction block
MPhotoVMNumericsPtr GetPhotoVMNumerics(MModelCache &cache) {
  return cache.Get<MPhotoVMNumerics>("photovm", [&cache] { return ReadPhotoVMNumerics(cache.Tune()); });
}

// Build the immutable direct dilepton process for one vector-meson channel
std::shared_ptr<const amplitude::ProcessDefinition> MPhotoVM::ProcessDefinitionFor(const std::string &channel) {
  const PhotoVMChannel channel_info = ResolvePhotoVMChannel(channel);
  return std::make_shared<amplitude::AnalyticProcess>(
      "MPHOTOVM", "ygg_" + channel_info.label, "charged lepton pair from " + channel_info.label,
      DirectDecayStructure(), PhotoVMProcessAccepts, [](const LORENTZSCALAR &) { return DirectDecayStructure(); });
}

// Initialize vector-meson and mode-dependent Sudakov parameters before copies
void MPhotoVM::InitializeParameters(MProcessSetup &setup) {
  if (!PhotoVMProcessAccepts(setup.lts.decaytree)) {
    throw std::invalid_argument("MPhotoVM::InitializeParameters: unsupported direct fermion pair");
  }
  if (!setup.model_tune) { throw std::invalid_argument("MPhotoVM::InitializeParameters: missing model tune"); }
  MModelCache &cache    = RequireModelCache(setup.lts.model_cache, setup.model_tune, "MPhotoVM initialization");
  const auto   numerics = GetPhotoVMNumerics(cache);
  if (setup.excitation != 0 && setup.lts.process.PHOTO_DISSOCIATION == DissociationType::Hera &&
      !numerics->dissociation.contains(ResolvePhotoVMChannel(setup.channel).vm_pdg)) {
    throw std::invalid_argument("PARAM_NSTAR.MODEL: photoprod_diss has no row for ygg[" + setup.channel + "]");
  }
  if (JMRTUsesFittedGluon(*numerics)) { return; }
  MPhotoQCD::EnsureSudakov(setup.lts, setup.soft_model, "MPhotoVM::InitializeParameters");
}

// Construct the amplitude and bind its immutable channel process definition
MPhotoVM::MPhotoVM(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune_snapshot, const std::string &channel_in,
                   std::shared_ptr<const amplitude::ProcessDefinition> definition)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(RequireModelCache(lts.model_cache, model_tune_snapshot, "MPhotoVM").TunePtr()),
      soft_model(model_tune->Soft()),
      pdg(ResolvePhotoVMChannel(channel_in).vm_pdg),
      dirac("DIRAC") {
  if (!PhotoVMProcessAccepts(lts.decaytree)) {
    throw std::invalid_argument("MPhotoVM: unsupported direct fermion pair");
  }
  MModelCache &cache = *lts.model_cache;
  numerics           = GetPhotoVMNumerics(cache);
  EnsureSudakov(lts);
  const auto row = numerics->dissociation.find(pdg);
  if (GetNstarParam(cache)->Model("ygg").photo == DissociationType::Hera && row != numerics->dissociation.end()) {
    diss = &row->second;
    const double w02 = pow2(diss->W0);
    diss_forward = JMRTProductionAmplitude(lts, *numerics, pdg, lts.PDG.FindByPDG(pdg), w02) / w02;
  }
}

// Acquire shared read-only Sudakov and UGD tables
void MPhotoVM::EnsureSudakov(gra::LORENTZSCALAR &lts) const {
  if (JMRTUsesFittedGluon(*numerics)) { return; }
  gra::MPhotoQCD::EnsureSudakov(lts, soft_model, "MPhotoVM::EnsureSudakov");
}

// Evaluate gamma p to heavy vector meson to dilepton kinematics
double MPhotoVM::Amp2(gra::LORENTZSCALAR &lts) const {
  // Preserve both exact photon-emitter and target forward systems through LHE
  // conversion
  lts.exact_forward_photon_kinematics  = true;
  const PhotoVMFinalState state        = ResolvePhotoVMFinalState(lts);
  const gra::MParticle   &vm           = lts.PDG.FindByPDG(pdg);
  const double            q2           = gra::MPhotoQCD::CentralMass2(lts);
  if (!(q2 > 0.0) || !std::isfinite(q2)) {
    throw AmplitudeFailure("MPhotoVM::Amp2: central invariant mass must be positive");
  }

  lts.hard_color_flows.clear();
  lts.hamp = BuildPhotoVMHelicityAmplitudes(lts, state, pdg, *numerics, vm, q2, dirac, *soft_model, diss, diss_forward);
  gra::flux::ConfigurePhotoLayout(lts);
  return gra::flux::PhotoInitialSpinAverage(lts) * gra::SquaredNorm(lts.hamp);
}

// Compute the elastic photon-proton differential cross section in nb/GeV^2
// d sigma/dt = [d sigma/dt]_0 exp[B(W)t]
double MPhotoVM::GammaPDSigmaDt(gra::LORENTZSCALAR &lts, double W, double t) const {
  EnsureSudakov(lts);
  if (!(W > 0.0) || !std::isfinite(W)) { throw std::invalid_argument("MPhotoVM::GammaPDSigmaDt: W must be positive"); }
  if (!std::isfinite(t) || t > 0.0) {
    throw std::invalid_argument("MPhotoVM::GammaPDSigmaDt: t must be finite and non-positive");
  }

  const gra::MParticle &vm           = lts.PDG.FindByPDG(pdg);
  const double          w2           = pow2(W);
  if (!(W > vm.mass + gra::PDG::mp)) { return 0.0; }
  const double forward = std::norm(JMRTProductionAmplitude(lts, *numerics, pdg, vm, w2) / w2) / (16.0 * gra::math::PI);
  const double slope   = JMRTSlopeB(*numerics, pdg, w2);
  return forward * gra::form::ExpSlopeWeight(slope, t) * gra::PDG::GeV2barn * 1.0e9;
}

// Compute the elastic photon-proton cross section through the selected |t| limit
// sigma = [d sigma/dt]_0 [1-exp(-B|t|max)]/B
double MPhotoVM::GammaPCrossSection(gra::LORENTZSCALAR &lts, double W, double abs_t_max) const {
  EnsureSudakov(lts);
  if (!(W > 0.0) || !std::isfinite(W)) {
    throw std::invalid_argument("MPhotoVM::GammaPCrossSection: W must be positive");
  }
  if (!(abs_t_max > 0.0) || !std::isfinite(abs_t_max)) {
    throw std::invalid_argument("MPhotoVM::GammaPCrossSection: abs_t_max must be positive");
  }

  const gra::MParticle &vm           = lts.PDG.FindByPDG(pdg);
  const double          w2           = pow2(W);
  if (!(W > vm.mass + gra::PDG::mp)) { return 0.0; }
  const double forward = std::norm(JMRTProductionAmplitude(lts, *numerics, pdg, vm, w2) / w2) / (16.0 * gra::math::PI);
  const double slope   = JMRTSlopeB(*numerics, pdg, w2);
  if (!(slope > 0.0)) { return 0.0; }
  const double t_integral = -std::expm1(-slope * abs_t_max) / slope;
  return forward * t_integral * gra::PDG::GeV2barn * 1.0e9;
}

// Compute the channel t slope in GeV^-2
double MPhotoVM::GammaPTSlope(double W) const {
  if (!(W > 0.0) || !std::isfinite(W)) { throw std::invalid_argument("MPhotoVM::GammaPTSlope: W must be positive"); }
  return JMRTSlopeB(*numerics, pdg, pow2(W));
}

// Clear shower color flow for color-singlet vector-meson decays
void MPhotoVM::SampleColorFlow(gra::LORENTZSCALAR &lts) const {
  for (auto &branch : lts.decaytree) { branch.p.color_flow.clear(); }
}

}  // namespace gra
