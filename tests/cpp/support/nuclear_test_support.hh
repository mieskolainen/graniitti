// Shared deterministic nuclear physics controls for C++ tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef NUCLEAR_TEST_SUPPORT_HH
#define NUCLEAR_TEST_SUPPORT_HH

#include <cmath>
#include <string>

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Nuclear/MBreakup.h"
#include "Graniitti/Nuclear/MGlauber.h"
#include "Graniitti/Nuclear/MPhoto.h"
#include "Graniitti/Nuclear/MUPC.h"

namespace gra::test {

// Compute one normalized analytic inelastic nucleon profile for unit tests
inline nuclear::NNProfile NNProfile(const double sigma = 70.0,
                                    const double omega = 0.15) {
  nuclear::NNProfile profile;
  profile.sigma = sigma;
  profile.omega = omega;
  profile.fingerprint = "test-nn-profile-v1-" + std::to_string(sigma) + "-" +
                        std::to_string(omega);
  constexpr std::size_t intervals = 256;
  constexpr double b_max = 8.0;
  constexpr double central = 0.8;
  const double slope = 0.1 * sigma / (2.0 * math::PI * central);
  profile.b_node.resize(intervals + 1, 0.0);
  profile.inelastic.resize(intervals + 1, 0.0);
  for (std::size_t i = 0; i <= intervals; ++i) {
    const double b =
        b_max * static_cast<double>(i) / static_cast<double>(intervals);
    profile.b_node[i] = b;
    profile.inelastic[i] = central * std::exp(-b * b / (2.0 * slope));
  }
  const double area =
      math::LinearRadialIntegral(profile.b_node, profile.inelastic);
  const double normalization = 0.1 * sigma / (2.0 * math::PI * area);
  for (double &value : profile.inelastic) {
    value *= normalization;
  }
  return profile;
}

// Fill compact Glauber quadratures with one analytic elementary profile
inline void Glauber(nuclear::GlauberParam &param, const double omega = 0.15) {
  param.profile = NNProfile(70.0, omega);
  param.b_max = 20.0;
  param.q_max = 3.0;
  param.b_nodes = 32;
  param.q_nodes = 32;
  param.profile_nodes = 256;
  param.fluctuation = {omega, omega > 0.0 ? 8U : 1U, 2.0e-14, 10000, 100.0};
}

// Compute one channel-dependent elementary photonuclear profile for tests
inline nuclear::PhotoProfile PhotoProfile(const double sigma = 15.0,
                                          const double eta = 0.1,
                                          const double omega = 0.15) {
  return {sigma, 4.0, eta, omega, 1.0e-3, 3.0};
}

// Fill compact photonuclear geometry quadratures for tests
inline void Photo(nuclear::PhotoParam &param) {
  param.phase_scale = 1.5;
  param.b_nodes = 32;
  param.z_nodes = 32;
  param.table.qt_max = 0.5;
  param.table.qz_max = 0.25;
  param.table.qt_nodes = 129;
  param.table.qz_nodes = 65;
  param.table.series_terms = 40;
  param.table.series_abs_tol = 1.0e-11;
  param.fluctuation = {0.0, 8, 2.0e-14, 10000, 100.0};
  param.shadow.pdf_set         = "CT10nlo";
  param.shadow.pdf_member = 0;
  param.shadow.dpdf_set        = "GKG18_DPDF_FitB_NLO";
  param.shadow.dpdf_member = 0;
  param.shadow.alpha0 = 1.0988;
  param.shadow.alpha_prime = 0.0;
  param.shadow.flux_b = 7.0;
  param.shadow.dpdf_diss       = 1.21;
  param.shadow.b_diff = 6.0;
  param.shadow.x_min = 1.0e-6;
  param.shadow.x_max = 0.1;
  param.shadow.q2_min = 2.7225;
  param.shadow.q2_max          = 4.0;
  param.shadow.sigma3_anchor_x = 1.0e-4;
  param.shadow.sigma3_fade_x = 1.0e-2;
  param.shadow.sigma3_power = 0.06;
  param.shadow.x_nodes = 131;
  param.shadow.scale_nodes = 65;
  param.shadow.integral_nodes = 64;
}

// Fill compact sampled-profile convolution controls for tests
inline void Convolution(nuclear::ConvolutionParam &param) {
  param.smooth_b_nodes = 96;
  param.sample_b_nodes = 96;
  param.smooth_kt_nodes = 96;
  param.sample_kt_nodes = 96;
}

// Fill the proton form-factor quadrature used by EMD tests
inline void BreakupNumerics(nuclear::BreakupParam &param) {
  param.proton.q_max = 20.0;
  param.proton.q_nodes = 256;
  param.response.tail_rel_tol = 1.0e-3;
}

} // namespace gra::test

#endif // NUCLEAR_TEST_SUPPORT_HH
