// Unit tests for parton density classes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <catch.hpp>
#include <cmath>
#include <cstdint>
#include <limits>
#include <set>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/PDF/MLHAPDF.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Sampling/MRandom.h"

// LHAPDF
#include "LHAPDF/LHAPDF.h"

namespace {

// Integrate selected proton momentum fractions on a logarithmic x grid
double IntegratedMomentumFractionForTest(const LHAPDF::PDF &pdf,
                                         const std::vector<int> &flavours,
                                         double q2, unsigned int nodes) {
  const double x_min = pdf.info().get_entry_as<double>("XMin", 0.0);
  const double x_max = pdf.info().get_entry_as<double>("XMax", 1.0);
  const double log_x_min = std::log(x_min);
  const double log_x_max = std::log(x_max);
  const auto [coordinates, weights] =
      gra::math::GaussLegendreRule(nodes, log_x_min, log_x_max);

  double momentum_fraction = 0.0;
  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    const double x = std::exp(coordinates[i]);
    double xfx_sum = 0.0;
    for (const int flavour : flavours) {
      xfx_sum += pdf.xfxQ2(flavour, x, q2);
    }
    momentum_fraction += weights[i] * x * xfx_sum;
  }
  return momentum_fraction;
}

struct PartonLuminosityForTest {
  double quark_antiquark = 0.0;
  double antiquark_quark = 0.0;
};

// Integrate both charge-conjugate terms of one pp parton luminosity
PartonLuminosityForTest
IntegratedPartonLuminosityForTest(const LHAPDF::PDF &pdf, int quark_pdg,
                                  double tau, double q2, unsigned int nodes) {
  const auto [coordinates, weights] =
      gra::math::GaussLegendreRule(nodes, std::log(tau), 0.0);
  PartonLuminosityForTest luminosity;

  for (std::size_t i = 0; i < coordinates.size(); ++i) {
    const double x1 = std::exp(coordinates[i]);
    const double x2 = tau / x1;
    const double quark1 = pdf.xfxQ2(quark_pdg, x1, q2) / x1;
    const double antiquark1 = pdf.xfxQ2(-quark_pdg, x1, q2) / x1;
    const double quark2 = pdf.xfxQ2(quark_pdg, x2, q2) / x2;
    const double antiquark2 = pdf.xfxQ2(-quark_pdg, x2, q2) / x2;
    luminosity.quark_antiquark += weights[i] * quark1 * antiquark2;
    luminosity.antiquark_quark += weights[i] * antiquark1 * quark2;
  }
  return luminosity;
}

// Compute the massless tree-level QED annihilation cross section
double LeptonPairCrossSectionForTest(double s, double alpha) {
  return 4.0 * gra::math::PI * alpha * alpha / (3.0 * s);
}

} // namespace

TEST_CASE("LUXqed satisfies the proton momentum sum rule",
          "[LUXqed][LHAPDF][sum-rule]") {
  gra::MLHAPDFStore store;
  const auto pdf = store.GetPDF("LUXqed17_plus_PDF4LHC15_nnlo_100", 0);
  const std::vector<int> flavours = pdf->flavors();
  const std::set<int> flavour_set(flavours.begin(), flavours.end());
  const double x_min = pdf->info().get_entry_as<double>("XMin", 0.0);
  const double x_max = pdf->info().get_entry_as<double>("XMax", 1.0);
  REQUIRE(x_min > 0.0);
  REQUIRE(x_max > x_min);

  for (const int required : {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 22}) {
    REQUIRE(flavour_set.count(required) == 1);
  }

  for (const double q2 : {1e2, 1e4, 1e6, 1e8}) {
    CAPTURE(q2);
    REQUIRE(pdf->inRangeXQ2(std::sqrt(x_min * x_max), q2));
    const double coarse =
        IntegratedMomentumFractionForTest(*pdf, flavours, q2, 128);
    const double fine =
        IntegratedMomentumFractionForTest(*pdf, flavours, q2, 256);
    const double photon =
        IntegratedMomentumFractionForTest(*pdf, {22}, q2, 256);

    REQUIRE(std::abs(fine - coarse) < 2e-4);
    REQUIRE(fine == Approx(1.0).margin(5e-3));
    REQUIRE(photon > 0.0);
    REQUIRE(photon < fine);
  }
}

TEST_CASE(
    "Drell-Yan pp luminosity includes both charge-conjugate beam assignments",
    "[Drell-Yan][LHAPDF][quantum-numbers]") {
  gra::MLHAPDFStore store;
  const auto pdf = store.GetPDF("LUXqed17_plus_PDF4LHC15_nnlo_100", 0);
  const double collider_s = gra::math::pow2(13000.0);
  const double dilepton_mass = 90.0;
  const double q2 = dilepton_mass * dilepton_mass;
  const double tau = q2 / collider_s;

  for (const int quark_pdg : {1, 2, 3, 4, 5}) {
    CAPTURE(quark_pdg);
    const auto luminosity =
        IntegratedPartonLuminosityForTest(*pdf, quark_pdg, tau, q2, 192);
    REQUIRE(luminosity.quark_antiquark > 0.0);
    REQUIRE(luminosity.antiquark_quark > 0.0);
    REQUIRE(luminosity.quark_antiquark ==
            Approx(luminosity.antiquark_quark).epsilon(2e-12));
  }
}

TEST_CASE("TwoBodyPhaseSpace integrates the massless QED lepton-pair amplitude",
          "[M4Vec][TwoBodyPhaseSpace][QED]") {
  using gra::math::pow2;
  using gra::math::pow4;

  constexpr double alpha = 1.0 / 137.0;
  const double electric_charge = std::sqrt(4.0 * gra::math::PI * alpha);
  constexpr unsigned int samples = 20000;

  for (const double sqrt_s : {100.0, 1000.0}) {
    CAPTURE(sqrt_s);
    const double s = sqrt_s * sqrt_s;
    const gra::M4Vec beam1(0.0, 0.0, sqrt_s / 2.0, sqrt_s / 2.0);
    const gra::M4Vec beam2(0.0, 0.0, -sqrt_s / 2.0, sqrt_s / 2.0);
    const gra::M4Vec mother = beam1 + beam2;
    const std::vector<double> masses = {0.0, 0.0};
    std::vector<gra::M4Vec> final;
    gra::MRandom rng;
    rng.SetSeed(static_cast<std::uint64_t>(sqrt_s));
    gra::kinematics::MCW cross_section;

    for (unsigned int sample = 0; sample < samples; ++sample) {
      const auto phase_space = gra::kinematics::TwoBodyPhaseSpace(
          mother, mother.M(), masses, final, rng);
      REQUIRE(phase_space.GetN() == Approx(1.0));
      REQUIRE(phase_space.GetW() > 0.0);

      if (sample < 32) {
        REQUIRE(gra::math::CheckEMC(mother - final[0] - final[1], 1e-11));
        const double mass_shell_tolerance =
            32.0 * std::numeric_limits<double>::epsilon() *
            std::max({1.0, final[0].E() * final[0].E(),
                      final[1].E() * final[1].E()});
        REQUIRE(final[0].M2() == Approx(0.0).margin(mass_shell_tolerance));
        REQUIRE(final[1].M2() == Approx(0.0).margin(mass_shell_tolerance));
      }

      const double amplitude2 = 8.0 * pow4(electric_charge) / pow2(s) *
                                ((beam2 * final[0]) * (beam1 * final[1]) +
                                 (beam1 * final[0]) * (beam2 * final[1]));
      REQUIRE(amplitude2 >= 0.0);
      cross_section.Push(phase_space.GetW() * amplitude2 / (2.0 * s));
    }

    const double exact = LeptonPairCrossSectionForTest(s, alpha);
    const double roundoff =
        100.0 * std::numeric_limits<double>::epsilon() * exact;
    REQUIRE(std::abs(cross_section.Integral() - exact) <=
            6.0 * cross_section.IntegralError() + roundoff);
  }
}
