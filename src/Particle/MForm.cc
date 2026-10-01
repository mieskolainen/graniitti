// Form factors and proton structure functions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <array>
#include <cmath>
#include <stdexcept>
#include <string>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"

// Libraries
#include "json.hpp"

using gra::math::pow2;
using gra::math::pow3;
using gra::math::pow4;

using gra::PDG::mp;

namespace gra {

namespace form {

// Read and validate one PARAM_STRUCTURE block
ParamStore ParamStore::Read(const std::string &source_file,
                            const std::string &json_text) {
  try {
    const nlohmann::json j = nlohmann::json::parse(json_text);
    const std::string XID = "PARAM_STRUCTURE";
    ParamStore param;
    param.F2 = j.at(XID).at("F2").get<std::string>();
    param.EM = j.at(XID).at("EM").get<std::string>();
    param.QED_alpha = j.at(XID).at("QED_alpha").get<std::string>();
    if (param.F2 != "DL" && param.F2 != "CKMT") {
      throw std::invalid_argument(
          "form::ParamStore::F2 has unsupported value " + param.F2);
    }
    if (param.EM != "DIPOLE" && param.EM != "KELLY") {
      throw std::invalid_argument(
          "form::ParamStore::EM has unsupported value " + param.EM);
    }
    if (param.QED_alpha != "ZERO" && param.QED_alpha != "MG" &&
        param.QED_alpha != "LL") {
      throw std::invalid_argument(
          "form::ParamStore::QED_alpha has unsupported value " +
          param.QED_alpha);
    }
    return param;
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "form::ParamStore::Read: Error reading " + source_file + ": " +
        error.what());
  }
}

// Read flat-amplitude parameters from one immutable GENERAL JSON text
FlatParam FlatParam::Read(const std::string &source_file, const std::string &json_text) {
  try {
    const nlohmann::json j = nlohmann::json::parse(json_text);

    const std::string XID = "PARAM_FLAT";

    if (!j.at(XID).at("B").is_number()) {
      throw std::invalid_argument("PARAM_FLAT::B must be numeric");
    }
    const double slope = j.at(XID).at("B");
    if (!std::isfinite(slope) || slope < 0.0) {
      throw std::invalid_argument(
          "PARAM_FLAT::B must be finite and nonnegative");
    }
    return {slope};
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "form::FlatParam::Read: Error reading " + source_file + ": " +
        error.what());
  }
}

// Evaluate an exponential cross-section slope as an amplitude factor
// A(t)/A(t0) = exp[B(t-t0)/2]
double ExpSlopeAmplitude(double B, double delta_t) {
  return std::exp(0.5 * B * delta_t);
}

// Evaluate an exponential cross-section slope as a weight factor
// |A(t)/A(t0)|^2 = exp[B(t-t0)]
double ExpSlopeWeight(double B, double delta_t) {
  return std::exp(B * delta_t);
}

// Proton inelastic structure function F2(x,Q^2) parametrization
//
// The basic idea is that at low-Q^2, a fully non-perturbative description
// (parametrization) is needed
// At high Q^2, DGLAP evolution could be done in log(Q^2) starting from the
// input description
//
// Now, some (very) classic ones have been implemented. Add new one here!
//
// [REFERENCE: Donnachie, Landshoff, arxiv.org/abs/hep-ph/9305319]
// [REFERENCE: Capella, Kaidalov, Merino, Tran Tranh Van, arxiv.org/abs/hep-ph/9405338v1]
//
double F2xQ2(double xbj, double Q2, const ParamStore &param) {
  if (!std::isfinite(xbj) || !std::isfinite(Q2) || !(xbj > 0.0) || xbj > 1.0 ||
      !(Q2 > 0.0)) {
    return 0.0;
  }
  if (param.F2 == "DL") {
    // Small-x fit of Eq. (4), without the paper's large-x extension
    constexpr double A = 0.324;
    constexpr double B = 0.098;

    constexpr double DELTA_P = 0.0808;
    constexpr double DELTA_R = 0.5475;

    constexpr double a = 0.561991692786383;
    constexpr double b = 0.011133;

    const double F2 =
        A * std::pow(xbj, -DELTA_P) * std::pow(Q2 / (Q2 + a), 1 + DELTA_P) +
        B * std::pow(xbj, 1 - DELTA_R) * std::pow(Q2 / (Q2 + b), DELTA_R);

    return F2;
  }

  else if (param.F2 == "CKMT") {
    constexpr double A = 0.1502;
    constexpr double B_u = 1.2064;
    constexpr double B_d = 0.1798;
    constexpr double alpha_R = 0.4150;
    constexpr double DELTA_0 = 0.07684;

    constexpr double a = 0.2631;
    constexpr double b = 0.6452;
    constexpr double c = 3.5489;
    constexpr double d = 1.1170;

    const double n_Q2 = (3.0 / 2.0) * (1 + Q2 / (Q2 + c));
    const double DELTA_Q2 = DELTA_0 * (1 + (2 * Q2) / (Q2 + d));

    const double C1 = std::pow(Q2 / (Q2 + a), 1.0 + DELTA_Q2);
    const double C2 = std::pow(Q2 / (Q2 + b), alpha_R);

    const double F2 =
        A * std::pow(xbj, -DELTA_Q2) * std::pow(1 - xbj, n_Q2 + 4.0) * C1 +
        std::pow(xbj, 1.0 - alpha_R) *
            (B_u * std::pow(1 - xbj, n_Q2) +
             B_d * std::pow(1 - xbj, n_Q2 + 1.0)) *
            C2;

    return F2;
  } else {
    throw std::invalid_argument(
        "gra::form::F2xQ2: unknown ParamStore::F2 = " + param.F2);
  }
}

// Compute the R1998 fit to R = sigma_L/sigma_T with its low-Q2 continuation
// R = (R_a+R_b+R_c)/3, scaled by Q2/0.2 below Q2 = 0.2 GeV2
// [REFERENCE: Abe et al., Phys. Lett. B452 (1999) 194, arXiv:hep-ex/9808028]
// [REFERENCE: Jefferson Lab Hall C resonance-data archive, r1998.f]
double RxQ2(double xbj, double Q2) {
  if (!std::isfinite(xbj) || !std::isfinite(Q2) || !(xbj > 0.0) || xbj > 1.0 ||
      !(Q2 > 0.0)) {
    return 0.0;
  }

  constexpr double q2_limit = 0.2;
  const double q2_eval = std::max(Q2, q2_limit);
  const double x2 = pow2(xbj);
  const double x3 = x2 * xbj;
  const double fac =
      1.0 + 12.0 * q2_eval / (q2_eval + 1.0) * pow2(0.125) / (pow2(0.125) + x2);
  const double rlog = fac / std::log(q2_eval / 0.04);

  constexpr std::array<double, 6> a = {4.8520e-2,  5.4704e-1, 2.0621,
                                       -3.8036e-1, 5.0896e-1, -2.8548e-2};
  constexpr std::array<double, 6> b = {4.8051e-2,  6.1130e-1, -3.5081e-1,
                                       -4.6076e-1, 7.1697e-1, -3.1726e-2};
  constexpr std::array<double, 6> c = {5.7654e-2, 4.6441e-1, 1.8288,
                                       1.2371e1,  -4.3104e1, 4.1741e1};

  const double qpa = a[1] / std::pow(pow4(q2_eval) + pow4(a[2]), 0.25) *
                     (1.0 + a[3] * xbj + a[4] * x2);
  const double ra = a[0] * rlog + qpa * std::pow(xbj, a[5]);

  const double qpb = (1.0 + b[3] * xbj + b[4] * x2) *
                     (b[1] / q2_eval + b[2] / (pow2(q2_eval) + 0.09));
  const double rb = b[0] * rlog + qpb * std::pow(xbj, b[5]);

  const double q2_threshold = c[3] * xbj + c[4] * x2 + c[5] * x3;
  const double qpc =
      c[1] / std::sqrt(pow2(q2_eval - q2_threshold) + pow2(c[2]));
  const double rc = c[0] * rlog + qpc;

  double ratio = (ra + rb + rc) / 3.0;
  if (Q2 < q2_limit) {
    ratio *= Q2 / q2_limit;
  }
  return std::max(ratio, 0.0);
}

// Compute F_L from F_2 and R including the exact target-mass coefficient
// F_L = (1+4x^2 mp^2/Q2) F_2 R/(1+R)
// [REFERENCE: Luszczak, Schafer and Szczurek, JHEP 05 (2018) 064, Eq. (2.14)]
double FLxQ2(double xbj, double Q2, const ParamStore &param) {
  if (!std::isfinite(xbj) || !std::isfinite(Q2) || !(xbj > 0.0) || xbj > 1.0 ||
      !(Q2 > 0.0)) {
    return 0.0;
  }
  const double rho = 4.0 * pow2(xbj * mp) / Q2;
  const double ratio = RxQ2(xbj, Q2);
  return (1.0 + rho) * F2xQ2(xbj, Q2, param) * ratio / (1.0 + ratio);
}

// Compute F_1 from F_2 and F_L including the exact target-mass coefficient
// F_1 = [(1+4x^2 mp^2/Q2)F_2-F_L]/(2x)
// [REFERENCE: Luszczak, Schafer and Szczurek, JHEP 05 (2018) 064, Eq. (2.14)]
double F1xQ2(double xbj, double Q2, const ParamStore &param) {
  if (!std::isfinite(xbj) || !std::isfinite(Q2) || !(xbj > 0.0) || xbj > 1.0 ||
      !(Q2 > 0.0)) {
    return 0.0;
  }
  const double rho = 4.0 * pow2(xbj * mp) / Q2;
  return ((1.0 + rho) * F2xQ2(xbj, Q2, param) -
          FLxQ2(xbj, Q2, param)) /
         (2.0 * xbj);
}

// Compute the proton magnetic moment in nuclear magneton units
// mu_p/mu_N = 2.792847337
double mu_ratio() { return 2.792847337; }

// ============================================================================
// Proton electromagnetic form factors, input Q^2 as positive

// kT unintegrated coherent EPA photon flux as in:
//
// [REFERENCE: Drees, Zeppenfeld, Phys. Rev. D39 (1989) 2536, doi:10.1103/PhysRevD.39.2536]
// [REFERENCE: Budnev, Ginzburg, Meledin, Serbo, Phys. Rept. 15 (1975) 181, doi:10.1016/0370-1573(75)90009-5]
// [REFERENCE: Luszczak, Schaefer, Szczurek, arxiv.org/abs/1802.03244]
//
// Form factors:
// [REFERENCE: Punjabi et al., arxiv.org/abs/1503.01452v4]
//
// Proton electromagnetic form factors: Basic notions, present
// achievements and future perspectives, Physics Reports, 2015
// <www.sciencedirect.com/science/article/pzi/S0370157314003184>
//
// [REFERENCE: Kelly, Phys. Rev. C70 (2004) 068202, doi:10.1103/PhysRevC.70.068202]
//
//
// Compute the proton Dirac electromagnetic form factor
// F_1 = (G_E + tau G_M)/(1+tau), tau = Q2/(4 mp^2)
double F1(double Q2, const ParamStore &param) {
  Q2 = std::abs(Q2);
  const double tau = Q2 / pow2(2 * mp);

  return 1.0 / (tau + 1) * G_E(Q2, param) +
         tau / (tau + 1) * G_M(Q2, param);
}

// Compute the proton Pauli electromagnetic form factor
// F_2 = (G_M - G_E)/(1+tau), tau = Q2/(4 mp^2)
double F2(double Q2, const ParamStore &param) {
  Q2 = std::abs(Q2);
  const double tau = Q2 / pow2(2 * mp);

  return -1.0 / (tau + 1) * G_E(Q2, param) +
         1.0 / (tau + 1) * G_M(Q2, param);
}

// Rosenbluth separation:
// low-Q^2 dominated by G_E, high-Q^2 dominated by G_M

// "Sachs Form Factor" goes as follows:
// G_E(0) = 1 for proton, 0 for neutron
// G_M(0) = mu_p for proton, mu_n for neutrons

// <http://www.scholarpedia.org/article/Nucleon_Form_factors>

// Compute the selected electric Sachs form factor G_E(Q2)
double G_E(double Q2, const ParamStore &param) {
  Q2 = std::abs(Q2); // For safety

  if (param.EM == "DIPOLE") {
    return G_E_DIPOLE(Q2);
  }
  return G_E_KELLY(Q2);
}

// Compute the selected magnetic Sachs form factor G_M(Q2)
double G_M(double Q2, const ParamStore &param) {
  Q2 = std::abs(Q2); // For safety

  if (param.EM == "DIPOLE") {
    return G_M_DIPOLE(Q2);
  }
  return G_M_KELLY(Q2);
}

// Evaluate the dipole parametrization of the electric Sachs form factor
// G_E(Q2) = [1+Q2/0.71 GeV2]^-2
double G_E_DIPOLE(double Q2) {
  Q2 = std::abs(Q2);
  constexpr double lambda2 = 0.71; // Dipole parameter GeV^2
  return 1.0 / pow2(1.0 + Q2 / lambda2);
}
// Compute G_M(Q2) = mu_p [1+Q2/0.71 GeV2]^-2
double G_M_DIPOLE(double Q2) {
  Q2 = std::abs(Q2);
  constexpr double lambda2 = 0.71; // Dipole parameter GeV^2
  return mu_ratio() / pow2(1.0 + Q2 / lambda2);
}

// Evaluate the Kelly parametrization of the electric Sachs form factor
// G_E = (1+a1 tau)/(1+b1 tau+b2 tau^2+b3 tau^3)
//
// [REFERENCE: Kelly, journals.aps.org/prc/pdf/10.1103/PhysRevC.70.068202]
double G_E_KELLY(double Q2) {
  Q2 = std::abs(Q2);
  constexpr std::array<double, 2> a = {1, -0.24};
  constexpr std::array<double, 3> b = {10.98, 12.82, 21.97};

  const double tau = Q2 / pow2(2 * mp);

  // Numerator
  double num = 0.0;  // 0
  num += a[0];       // a_0 tau^0
  num += a[1] * tau; // a_1 tau^1

  // Denominator
  double den = 1.0;        // 1.0
  den += b[0] * tau;       // b_1 tau^1
  den += b[1] * pow2(tau); // b_2 tau^2
  den += b[2] * pow3(tau); // b_3 tau^3

  return num / den;
}

// Evaluate the Kelly parametrization of the magnetic Sachs form factor
// G_M = mu_p(1+a1 tau)/(1+b1 tau+b2 tau^2+b3 tau^3)
double G_M_KELLY(double Q2) {
  Q2 = std::abs(Q2);
  constexpr std::array<double, 2> a = {1, 0.12};
  constexpr std::array<double, 3> b = {10.97, 18.86, 6.55};

  const double tau = Q2 / pow2(2 * mp);

  // Numerator
  double num = 0.0;  // 0
  num += a[0];       // a_0 tau^0
  num += a[1] * tau; // a_1 tau^1

  // Denominator
  double den = 1.0;        // 1.0
  den += b[0] * tau;       // b_1 tau^1
  den += b[1] * pow2(tau); // b_2 tau^2
  den += b[2] * pow3(tau); // b_3 tau^3

  return mu_ratio() * num / den;
}

} // namespace form
} // namespace gra
