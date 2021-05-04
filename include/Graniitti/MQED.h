// Running QED coupling [HEADER ONLY FILE]
//
// (c) 2017-2021 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.


#ifndef MQED_H
#define MQED_H

#include <vector>

// Own
#include "Graniitti/MMath.h"


namespace gra {
namespace qed {

constexpr double alpha_ref = 0.0072973525693;  // reference value of alpha at scale Q
constexpr double Q_ref     = 0.0005109989461;  // reference scale (GeV)
constexpr double alpha_0 =
    0.0072973525693;  // value of alpha at Q = 0 (fine structure constant)^{-1}

// -----------------------------------------------------------------------

// Get number of charged leptons at given scale
inline int get_N_leptons(double Q) {
  if (Q > 1.77686E+00) { return 3; }
  if (Q > 1.056583745E-01) { return 2; }
  return 1;
}

// The coefficient of one loop beta function
inline double beta0_QED(int N_lf) { return -1.0 / (3.0 * math::PI) * N_lf; }

// QED running coupling
inline double alpha_QED(double Q2, const std::string& order = "LL") {
  // Freeze the coupling at Q < Q_ref
  if (order == "ZERO" || Q2 < (Q_ref * Q_ref)) { return alpha_0; }

  // Input
  const double Q    = math::msqrt(Q2);
  const int    N_lf = get_N_leptons(Q);

  // Beta-function quantities
  const double beta0 = beta0_QED(N_lf);
  const double R     = std::log(Q2 / (Q_ref * Q_ref));

  if (order == "LL") {  // leading log [one loop geometric series resummation]
    return alpha_ref / (1.0 + alpha_ref * beta0 * R);
  } else {
    throw std::invalid_argument("alpha_QED:: Unknown order parameter: " + order);
  }
}

// QED coupling at Q = 0
inline double alpha_QED() { return alpha_QED(0, "ZERO"); }

// QED: Electric charge in natural units ~ 0.3
inline double e_QED(double Q2) { return math::msqrt(alpha_QED(Q2) * 4.0 * math::PI); }
inline double e_QED() { return math::msqrt(alpha_QED() * 4.0 * math::PI); }


}  // namespace qed
}  // namespace gra

#endif