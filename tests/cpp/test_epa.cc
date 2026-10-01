// Unit tests for elastic and inelastic kT-EPA physics interfaces
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/HelAmps_sm_lepton_masses.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Tech/MException.h"

namespace {

using gra::aux::indices;

struct ElasticFluxComponents {
  double electric           = 0.0;
  double magnetic_per_state = 0.0;
};

struct InelasticFluxComponents {
  double electric           = 0.0;
  double magnetic_per_state = 0.0;
  double xbj                = 0.0;
  double Q2                 = 0.0;
};

struct CovariantMG5Point {
  gra::LORENTZSCALAR      lts;
  std::vector<gra::M4Vec> final;
};

// Build one collinear gamma-gamma point with exact massive proton beams
gra::LORENTZSCALAR MakeCollinearEPAEvent(double beam_pz, double x1, double x2) {
  gra::LORENTZSCALAR lts;
  const double       beam_energy = std::sqrt(beam_pz * beam_pz + gra::PDG::mp * gra::PDG::mp);
  lts.pbeam1                     = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  lts.pbeam2                     = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  lts.beam1.mass                 = gra::PDG::mp;
  lts.beam2.mass                 = gra::PDG::mp;
  lts.s                          = (lts.pbeam1 + lts.pbeam2).M2();
  lts.x1                         = x1;
  lts.x2                         = x2;
  lts.s_hat                      = x1 * x2 * lts.s;
  return lts;
}

// Build one exact elastic pp to pXp phase-space point
gra::LORENTZSCALAR MakeElasticEPAEvent(double beam_pz, double x1, double x2, double t1, double t2, double phi1,
                                       double phi2) {
  gra::LORENTZSCALAR lts;
  const double       beam_energy = std::sqrt(beam_pz * beam_pz + gra::PDG::mp * gra::PDG::mp);
  lts.pbeam1                     = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  lts.pbeam2                     = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  lts.beam1.pdg                  = gra::PDG::PDG_p;
  lts.beam1.chargeX3             = 3;
  lts.beam1.spinX2               = 1;
  lts.beam1.mass                 = gra::PDG::mp;
  lts.beam2                      = lts.beam1;
  if (!gra::kinematics::BuildForwardParticleXiT(lts.pbeam1, x1, t1, phi1, true, lts.pfinal[1]) ||
      !gra::kinematics::BuildForwardParticleXiT(lts.pbeam2, x2, t2, phi2, false, lts.pfinal[2])) {
    throw std::invalid_argument("MakeElasticEPAEvent: invalid forward point");
  }

  lts.q1            = lts.pbeam1 - lts.pfinal[1];
  lts.q2            = lts.pbeam2 - lts.pfinal[2];
  lts.pfinal[0]     = lts.q1 + lts.q2;
  lts.forward_mass2 = {gra::math::pow2(gra::PDG::mp), gra::math::pow2(gra::PDG::mp)};
  lts.s             = (lts.pbeam1 + lts.pbeam2).M2();
  lts.x1            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam1, lts.pfinal[1], true);
  lts.x2            = gra::kinematics::LongitudinalMomentumLoss(lts.pbeam2, lts.pfinal[2], false);
  lts.xi1           = lts.x1;
  lts.xi2           = lts.x2;
  lts.has_xi1       = true;
  lts.has_xi2       = true;
  lts.t1            = lts.q1.M2();
  lts.t2            = lts.q2.M2();
  lts.qt1           = lts.q1.Pt();
  lts.qt2           = lts.q2.Pt();
  lts.m2            = lts.pfinal[0].M2();
  lts.s_hat         = lts.m2;
  return lts;
}

// Compute central M2 and rapidity as functions of the exact emitter losses
std::array<double, 2> CentralMassRapidityFromX(const gra::LORENTZSCALAR &lts, double x1, double x2) {
  const double retained1     = 1.0 - x1;
  const double retained2     = 1.0 - x2;
  const double plus1         = retained1 * lts.pbeam1.LightconePos();
  const double minus2        = retained2 * lts.pbeam2.LightconeNeg();
  const double mt1sq         = lts.pfinal[1].M2() + lts.pfinal[1].Pt2();
  const double mt2sq         = lts.pfinal[2].M2() + lts.pfinal[2].Pt2();
  const double minus1        = mt1sq / plus1;
  const double plus2         = mt2sq / minus2;
  const double central_plus  = (lts.pbeam1 + lts.pbeam2).LightconePos() - plus1 - plus2;
  const double central_minus = (lts.pbeam1 + lts.pbeam2).LightconeNeg() - minus1 - minus2;
  const double central_pt2   = (lts.pfinal[1] + lts.pfinal[2]).Pt2();
  return {central_plus * central_minus - central_pt2, 0.5 * std::log(central_plus / central_minus)};
}

// Evaluate the exact dx1 dx2 over dM2 dY Jacobian by central differences
double NumericalEmitterJacobian(const gra::LORENTZSCALAR &lts) {
  constexpr double step     = 1.0e-6;
  const auto       x1_plus  = CentralMassRapidityFromX(lts, lts.x1 + step, lts.x2);
  const auto       x1_minus = CentralMassRapidityFromX(lts, lts.x1 - step, lts.x2);
  const auto       x2_plus  = CentralMassRapidityFromX(lts, lts.x1, lts.x2 + step);
  const auto       x2_minus = CentralMassRapidityFromX(lts, lts.x1, lts.x2 - step);
  const double     dm2_dx1  = (x1_plus[0] - x1_minus[0]) / (2.0 * step);
  const double     dy_dx1   = (x1_plus[1] - x1_minus[1]) / (2.0 * step);
  const double     dm2_dx2  = (x2_plus[0] - x2_minus[0]) / (2.0 * step);
  const double     dy_dx2   = (x2_plus[1] - x2_minus[1]) / (2.0 * step);
  return 1.0 / std::abs(dm2_dx1 * dy_dx2 - dm2_dx2 * dy_dx1);
}

// Split the current elastic scalar flux into electric and per-state magnetic
// terms
// [REFERENCE: Luszczak, Schafer and Szczurek, arXiv:1802.03244, Eq. (2.13)]
ElasticFluxComponents ElasticComponents(double x, double t, double pt) {
  const double pt2 = pt * pt;
  const double Q2  = std::abs(t);
  const double mp2 = gra::PDG::mp * gra::PDG::mp;
  const double electric_form =
      (4.0 * mp2 * std::pow(gra::form::G_E(Q2), 2) + Q2 * std::pow(gra::form::G_M(Q2), 2)) / (4.0 * mp2 + Q2);
  const double magnetic_form = std::pow(gra::form::G_M(Q2), 2);
  const double delta         = pt2 / (pt2 + x * x * mp2);
  const double common        = 16.0 * gra::math::PIPI * gra::qed::alpha_QED() / (gra::math::PI * x * pt2);

  return {common * (1.0 - x) * delta * delta * electric_form, common * x * x * delta * magnetic_form / 4.0};
}

// Split the current inelastic scalar flux into electric and per-state magnetic
// terms
// [REFERENCE: Luszczak, Schafer and Szczurek, arXiv:1802.03244, Eqs. (2.12) and (2.14)]
InelasticFluxComponents InelasticComponents(double x, double t, double pt, double remnant_mass2) {
  const double pt2         = pt * pt;
  const double Q2          = std::abs(t);
  const double mp2         = gra::PDG::mp * gra::PDG::mp;
  const double denom       = Q2 + remnant_mass2 - mp2;
  const double xbj         = Q2 / denom;
  const double delta_denom = pt2 + x * (remnant_mass2 - mp2) + x * x * mp2;
  const double delta       = pt2 / delta_denom;
  const double common      = 16.0 * gra::math::PIPI * gra::qed::alpha_QED() / (gra::math::PI * x * pt2);
  const double f2          = gra::form::F2xQ2(xbj, Q2);
  const double two_xf1     = 2.0 * xbj * gra::form::F1xQ2(xbj, Q2);

  return {common * (1.0 - x) * delta * delta * f2 / denom, common * x * x * delta * two_xf1 / (4.0 * xbj * xbj * denom),
          xbj, Q2};
}

// Integrate the defining DZ dipole spectrum in logarithmic Q2
// [REFERENCE: Drees and Zeppenfeld, Phys. Rev. D39 (1989) 2536]
double NumericalDZFlux(double x) {
  constexpr int    intervals = 20000;
  constexpr double log_span  = 30.0;
  const double     q2_min    = std::pow(gra::PDG::mp * x, 2) / (1.0 - x);
  const double     step      = log_span / static_cast<double>(intervals);
  const auto       integrand = [&](double log_ratio) {
    const double q2 = q2_min * std::exp(log_ratio);
    return std::pow(1.0 + q2 / 0.71, -4);
  };

  double integral = integrand(0.0) + integrand(log_span);
  for (int i = 1; i < intervals; ++i) {
    integral += (i % 2 == 0 ? 2.0 : 4.0) * integrand(static_cast<double>(i) * step);
  }
  integral *= step / 3.0;
  return gra::qed::alpha_QED() / (2.0 * gra::math::PI * x) * (1.0 + std::pow(1.0 - x, 2)) * integral;
}

// Construct one non-collinear closed hard point for the covariant MG5 embedding
CovariantMG5Point MakeCovariantMG5Point() {
  CovariantMG5Point point;
  point.lts.q1                              = gra::M4Vec(0.8, -0.5, 12.0, 10.0);
  point.lts.q2                              = gra::M4Vec(-0.8, 0.5, -12.0, 10.0);

  const double beam_pz     = 6500.0;
  const double beam_energy = std::sqrt(beam_pz * beam_pz + gra::PDG::mp * gra::PDG::mp);
  point.lts.pbeam1         = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  point.lts.pbeam2         = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  point.lts.pfinal[0]      = point.lts.q1 + point.lts.q2;
  point.lts.pfinal[1]      = point.lts.pbeam1 - point.lts.q1;
  point.lts.pfinal[2]      = point.lts.pbeam2 - point.lts.q2;
  point.lts.excite1        = true;
  point.lts.excite2        = true;
  point.lts.x1             = gra::kinematics::LongitudinalMomentumLoss(point.lts.pbeam1, point.lts.pfinal[1], true);
  point.lts.x2             = gra::kinematics::LongitudinalMomentumLoss(point.lts.pbeam2, point.lts.pfinal[2], false);
  point.lts.xi1            = point.lts.x1;
  point.lts.xi2            = point.lts.x2;
  point.lts.has_xi1        = true;
  point.lts.has_xi2        = true;
  point.lts.t1             = point.lts.q1.M2();
  point.lts.t2             = point.lts.q2.M2();
  point.lts.qt1            = point.lts.q1.Pt();
  point.lts.qt2            = point.lts.q2.Pt();
  point.lts.s              = (point.lts.pbeam1 + point.lts.pbeam2).M2();
  point.lts.m2             = point.lts.pfinal[0].M2();
  point.lts.s_hat          = point.lts.m2;

  const double final_pz = std::sqrt(71.0);
  point.final           = {gra::M4Vec(3.0, 4.0, final_pz, 10.0), gra::M4Vec(-3.0, -4.0, -final_pz, 10.0)};
  return point;
}

// Require two four-vectors to agree component by component
void RequireFourVectorNear(const gra::M4Vec &actual, const gra::M4Vec &expected, double tolerance) {
  REQUIRE(actual.Px() == Approx(expected.Px()).margin(tolerance));
  REQUIRE(actual.Py() == Approx(expected.Py()).margin(tolerance));
  REQUIRE(actual.Pz() == Approx(expected.Pz()).margin(tolerance));
  REQUIRE(actual.E() == Approx(expected.E()).margin(tolerance));
}

// Require two complex numbers to agree component by component
void RequireComplexNear(const std::complex<double> &actual, const std::complex<double> &expected, double tolerance) {
  REQUIRE(actual.real() == Approx(expected.real()).margin(tolerance));
  REQUIRE(actual.imag() == Approx(expected.imag()).margin(tolerance));
}

// Rotate one normalized real linear source in the MG5 helicity basis
std::array<std::complex<double>, 2> PerpendicularSource(const std::array<std::complex<double>, 2> &parallel) {
  return {-std::conj(parallel[1]), std::conj(parallel[0])};
}

// Project one source with the polarization vectors generated by MG5 HELAS
std::array<std::complex<double>, 2> HELASSourceCoefficients(const gra::M4Vec &source, const gra::M4Vec &k1,
                                                            const gra::M4Vec &k2) {
  const double     pair_product = k1 * k2;
  const gra::M4Vec transverse   = source - k1 * ((source * k2) / pair_product) - k2 * ((source * k1) / pair_product);
  const double     scale        = 1.0 / std::sqrt(-transverse.M2());
  double           momentum[4]  = {k1.E(), k1.Px(), k1.Py(), k1.Pz()};
  std::array<std::complex<double>, 2> coefficients = {};
  for (std::size_t index = 0; index < coefficients.size(); ++index) {
    const int            helicity        = index == 0 ? -1 : 1;
    std::complex<double> wavefunction[6] = {};
    MG5_sm_lepton_masses::vxxxxx(momentum, 0.0, helicity, -1, wavefunction);
    const std::complex<double> contraction =
        transverse.E() * std::conj(wavefunction[2]) - transverse.Px() * std::conj(wavefunction[3]) -
        transverse.Py() * std::conj(wavefunction[4]) - transverse.Pz() * std::conj(wavefunction[5]);
    coefficients[index] = -scale * contraction;
  }
  return coefficients;
}

}  // namespace

TEST_CASE("EPA loop transfers keep xi and both forward masses exact", "[EPA][kinematics][screening]") {
  constexpr double beam_pz = 6500.0;
  const auto       check   = [&](const double initial_mass, const double final_mass, const double xi, const double qx,
                         const double qy, const bool plus_side) {
    const double     initial_mass2 = gra::math::pow2(initial_mass);
    const double     final_mass2   = gra::math::pow2(final_mass);
    const double     energy        = std::sqrt(beam_pz * beam_pz + initial_mass2);
    const gra::M4Vec beam(0.0, 0.0, plus_side ? beam_pz : -beam_pz, energy);
    gra::M4Vec       forward;
    gra::M4Vec       transfer;
    double           t = 0.0;
    REQUIRE(gra::kinematics::BuildEPATransfer(beam, xi, qx, qy, initial_mass2, final_mass2, plus_side, forward,
                                                      transfer, t));
    const long double expected = -(static_cast<long double>(qx) * qx + static_cast<long double>(qy) * qy +
                                   static_cast<long double>(xi) * (final_mass2 - initial_mass2) +
                                   static_cast<long double>(xi) * xi * initial_mass2) /
                                 (1.0L - xi);
    REQUIRE(t == Approx(static_cast<double>(expected)).epsilon(5.0e-14));
    REQUIRE(forward.Px() == Approx(-qx).margin(1.0e-14));
    REQUIRE(forward.Py() == Approx(-qy).margin(1.0e-14));
    REQUIRE(gra::kinematics::LongitudinalMomentumLoss(beam, forward, plus_side) == Approx(xi).margin(2.0e-13));
    REQUIRE(forward.M2() == Approx(final_mass2).margin(2.0e-8));
    const gra::M4Vec difference = beam - forward;
    REQUIRE(transfer.Px() == Approx(difference.Px()).margin(2.0e-12));
    REQUIRE(transfer.Py() == Approx(difference.Py()).margin(2.0e-12));
    REQUIRE(transfer.Pz() == Approx(difference.Pz()).margin(2.0e-10));
    REQUIRE(transfer.E() == Approx(difference.E()).margin(2.0e-10));
    const auto        q  = transfer.Contravariant<long double>();
    const long double q2 = q[0] * q[0] - q[1] * q[1] - q[2] * q[2] - q[3] * q[3];
    const long double cancellation_scale =
        std::max({std::abs(expected), std::abs(q[0] * q[0]), std::abs(q[1] * q[1]), std::abs(q[2] * q[2]),
                  std::abs(q[3] * q[3]), std::numeric_limits<long double>::min()});
    REQUIRE(std::abs(q2 - expected) < 512.0L * std::numeric_limits<double>::epsilon() * cancellation_scale);
  };

  for (const bool plus_side : {false, true}) {
    for (const double xi : {1.0e-6, 0.037, 0.31}) {
      DYNAMIC_SECTION((plus_side ? "upper" : "lower") << ", xi=" << xi) {
        check(gra::PDG::mp, gra::PDG::mp, xi, 1.0e-5, -2.0e-5, plus_side);
        check(2.3, 3.1, xi, 0.43, -0.21, plus_side);
      }
    }
  }

  const double mp2         = gra::math::pow2(gra::PDG::mp);
  const double beam_energy = std::sqrt(beam_pz * beam_pz + mp2);
  gra::M4Vec   forward;
  gra::M4Vec   transfer;
  double       t = 0.0;
  CHECK_FALSE(gra::kinematics::BuildEPATransfer(gra::M4Vec(0.0, 0.0, beam_pz, beam_energy), 0.037, 0.43, -0.21, -mp2,
                                                mp2, true, forward, transfer, t));
  CHECK_FALSE(gra::kinematics::BuildEPATransfer(gra::M4Vec(0.0, 0.0, beam_pz, beam_energy), 0.037, 0.43, -0.21,
                                                std::numeric_limits<double>::quiet_NaN(), mp2, true, forward, transfer,
                                                t));
  CHECK_FALSE(gra::kinematics::BuildEPATransfer(gra::M4Vec(0.0, 0.0, beam_pz, beam_energy), 0.037, 0.43, -0.21,
                                                1.01 * mp2, mp2, true, forward, transfer, t));
}

TEST_CASE("Elastic photon currents are finite for exact forward legs", "[EPA][QED][current][validation]") {
  gra::LORENTZSCALAR valid = MakeElasticEPAEvent(6500.0, 1.0e-3, 2.0e-3, -0.011, -0.019, 0.37, -0.82);
  valid.model_cache =
      std::make_shared<gra::MModelCache>(gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
  const auto transitions = gra::spin::SpinHalfTransitions(false);

  const auto require_valid = [&](gra::LORENTZSCALAR &lts, const int leg) {
    const auto current = gra::qed::PhotonCurrentTransitions(lts, leg, transitions);
    REQUIRE(current.size() == transitions.size());
    for (const auto &row : current) { REQUIRE(gra::AllFinite(row)); }
  };
  for (const int leg : {1, 2}) {
    DYNAMIC_SECTION("beam leg " << leg) { require_valid(valid, leg); }
  }
}

TEST_CASE("kT-EPA uses the exact light-cone and hard-flux normalization", "[EPA][normalization]") {
  gra::LORENTZSCALAR lts = MakeElasticEPAEvent(100.0, 0.08, 0.13, -0.41, -0.76, 0.31, -1.14);
  lts.hamp               = {std::complex<double>(0.7, -0.2), std::complex<double>(-0.3, 0.5)};

  const double numerical_jacobian = NumericalEmitterJacobian(lts);
  const double velocity_difference =
      std::abs(lts.pfinal[1].Pz() / lts.pfinal[1].E() - lts.pfinal[2].Pz() / lts.pfinal[2].E());
  const double beta                  = gra::kinematics::beta12(lts.s, gra::PDG::mp, gra::PDG::mp);
  const double hard_mass2            = (lts.q1 + lts.q2).M2();
  const double reference_phase_space = 2.0 * lts.s * beta * lts.pfinal[1].E() * lts.pfinal[2].E() *
                                       velocity_difference * numerical_jacobian / hard_mass2;
  const double phase_space = gra::flux::ExactktEPAPhaseSpaceFactor(lts);
  REQUIRE(phase_space == Approx(reference_phase_space).epsilon(2.0e-8));
  REQUIRE(std::abs(hard_mass2 - lts.x1 * lts.x2 * lts.s) > 1.0e-6 * hard_mass2);

  const gra::LORENTZSCALAR high_energy = MakeElasticEPAEvent(6500.0, 0.001, 0.002, -0.01, -0.02, 0.4, -0.7);
  REQUIRE(gra::flux::ExactktEPAPhaseSpaceFactor(high_energy) ==
          Approx((1.0 - high_energy.x1) * (1.0 - high_energy.x2) / (high_energy.x1 * high_energy.x2)).epsilon(2.0e-4));

  const double amp2     = 2.3;
  const double flux1    = gra::flux::CohFlux(lts.x1, lts.t1, lts.qt1);
  const double flux2    = gra::flux::CohFlux(lts.x2, lts.t2, lts.qt2);
  const double factor   = flux1 * flux2 * phase_space;
  const auto   raw_hamp = lts.hamp;
  const double weighted = gra::flux::ApplyktEPAfluxes(amp2, lts).amp2;

  REQUIRE(weighted == Approx(factor * amp2).epsilon(1e-14));
  REQUIRE_FALSE(lts.proton_good_walker.has_value());
  REQUIRE(lts.pdf_xf1 == Approx(lts.xi1 * flux1).epsilon(1e-14));
  REQUIRE(lts.pdf_xf2 == Approx(lts.xi2 * flux2).epsilon(1e-14));
  REQUIRE(lts.exact_forward_photon_kinematics);
  REQUIRE(lts.muF == Approx(std::sqrt(lts.s_hat) / 2.0).epsilon(1e-14));
  REQUIRE(lts.muR == Approx(lts.muF));
  REQUIRE(lts.scalup == Approx(lts.muF));
  for (std::size_t i = 0; i < raw_hamp.size(); ++i) {
    RequireComplexNear(lts.hamp[i], raw_hamp[i] * std::sqrt(factor), 1e-14);
  }

  gra::LORENTZSCALAR fixed_lts   = MakeElasticEPAEvent(100.0, 0.08, 0.13, -0.41, -0.76, 0.31, -1.14);
  fixed_lts.hamp                 = raw_hamp;
  const double fixed_phase_space = 0.37 * phase_space;
  const double fixed_weighted    = gra::flux::ApplyktEPAcurrents(amp2, fixed_lts, fixed_phase_space).amp2;
  const double fixed_factor      = flux1 * flux2 * fixed_phase_space;
  REQUIRE(fixed_weighted == Approx(fixed_factor * amp2).epsilon(1e-14));
  for (const auto &i : indices(raw_hamp)) {
    RequireComplexNear(fixed_lts.hamp[i], raw_hamp[i] * std::sqrt(fixed_factor), 1e-14);
  }
}

// Check both recoil factors against independent proton spectra in dQ2 at finite photon fractions
// [REFERENCE: arXiv:1802.03244, Eq. (2.13)]
TEST_CASE("Two-photon EPA retains both finite-x recoil factors", "[EPA][normalization][recoil]") {
  const double x1 = GENERATE(0.03, 0.2, 0.7);
  const double x2 = GENERATE(0.03, 0.2, 0.7);
  const double phi = GENERATE(-0.9, 1.7);
  const double mass2 = gra::math::pow2(gra::PDG::mp);
  const double qmin1 = mass2 * x1 * x1 / (1.0 - x1);
  const double qmin2 = mass2 * x2 * x2 / (1.0 - x2);
  const double q1 = qmin1 + 0.4, q2 = qmin2 + 0.7;
  auto lts = MakeElasticEPAEvent(6500.0, x1, x2, -q1, -q2, phi, 0.3);
  const std::complex<double> hard(0.7, -0.2);
  lts.hamp = {hard};

  // Compute dn/(dx dQ2) without using the transverse EPA flux implementation
  const auto spectrum = [mass2](double x, double q, double qmin) {
    const double ge2 = gra::math::pow2(gra::form::G_E(q));
    const double gm2 = gra::math::pow2(gra::form::G_M(q));
    const double electric = (4.0 * mass2 * ge2 + q * gm2) / (4.0 * mass2 + q);
    return gra::qed::alpha_QED() / (gra::math::PI * x * q) *
           ((1.0 - x) * (1.0 - qmin / q) * electric + 0.5 * x * x * gm2);
  };
  // The high energy limit keeps finite x and neglects only masses and transfers over s
  const double expected = gra::math::pow2(16.0 * gra::math::PIPI) *
                          spectrum(x1, q1, qmin1) * spectrum(x2, q2, qmin2) / (x1 * x2);
  const auto weighted = gra::flux::ApplyktEPAfluxes(std::norm(hard), lts);
  CAPTURE(x1, x2, phi);
  REQUIRE(weighted.amp2 == Approx(expected * std::norm(hard)).epsilon(2.0e-4));
  REQUIRE(lts.hamp.size() == 1);
  RequireComplexNear(lts.hamp.front() / std::sqrt(expected), hard, 2.0e-4);
}

TEST_CASE("EPA sector weights preserve the hard-row ordering", "[EPA][nuclear][sector]") {
  const std::vector<std::complex<double>> raw        = {{0.7, -0.2}, {-0.3, 0.5}};
  std::vector<std::complex<double>>       amplitudes = raw;
  gra::flux::EPASectorWeights             weights;
  weights.amplitude_fraction = {0.2, 0.4, 0.6, 0.8};
  constexpr double common    = 1.7;

  REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(amplitudes, common, weights));
  REQUIRE(amplitudes.size() == raw.size() * weights.amplitude_fraction.size());
  for (std::size_t hard = 0; hard < raw.size(); ++hard) {
    for (std::size_t sector = 0; sector < weights.amplitude_fraction.size(); ++sector) {
      const std::size_t row = hard * weights.amplitude_fraction.size() + sector;
      RequireComplexNear(amplitudes[row], common * weights.amplitude_fraction[sector] * raw[hard], 1.0e-14);
    }
  }

  amplitudes                    = raw;
  weights.amplitude_fraction[2] = -0.1;
  REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(amplitudes, common, weights));
  REQUIRE(amplitudes.size() == raw.size() * weights.amplitude_fraction.size());
  RequireComplexNear(amplitudes[2], -0.1 * common * raw[0], 1.0e-14);
  weights.amplitude_fraction[2] = std::numeric_limits<double>::quiet_NaN();
  REQUIRE(gra::flux::ApplyEPAAmplitudeWeights(amplitudes, common, weights));
  REQUIRE_FALSE(std::isfinite(std::norm(amplitudes[2])));
}

TEST_CASE("kT-EPA throws AmplitudeFailure for absent generated forward states", "[EPA][inelastic][validation]") {
  gra::LORENTZSCALAR flux_lts = MakeElasticEPAEvent(100.0, 0.08, 0.13, -0.41, -0.76, 0.31, -1.14);
  flux_lts.excite1            = true;
  flux_lts.pfinal.resize(1);
  flux_lts.hamp = {std::complex<double>(0.7, -0.2)};
  REQUIRE_THROWS_AS(gra::flux::ApplyktEPAfluxes(1.0, flux_lts), gra::AmplitudeFailure);

  CovariantMG5Point point = MakeCovariantMG5Point();
  point.lts.excite1       = true;
  point.lts.pfinal.resize(1);
  const std::array<std::array<int, 2>, 4>       helicities = {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, 1},
                                                              std::array<int, 2>{1, -1}, std::array<int, 2>{1, 1}};
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (const auto &helicity : helicities) { components.push_back({helicity, {}, 0, 1.0, {}}); }

  gra::M4Vec k1;
  gra::M4Vec k2;
  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(point.lts, point.final, k1, k2));
  REQUIRE_THROWS_AS(gra::mg5helas::ContractEPAPhotonSources(point.lts, components, k1, k2), gra::AmplitudeFailure);

  CovariantMG5Point invalid_mass = MakeCovariantMG5Point();
  invalid_mass.lts.pfinal[1]     = gra::M4Vec();
  REQUIRE_THROWS_AS(gra::mg5helas::ContractEPAPhotonSources(invalid_mass.lts, components, k1, k2),
                    gra::AmplitudeFailure);
}

TEST_CASE("QED photon sources require complete spacelike forward kinematics", "[EPA][QED][Kinematics]") {
  const std::vector<std::pair<double, double>> transitions = {{-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<int>                       m_values    = {-1, 0, 1};

  SECTION("missing upper forward leg") {
    gra::LORENTZSCALAR lts = MakeElasticEPAEvent(6500.0, 0.02, 0.03, -0.2, -0.3, 0.4, -0.7);
    lts.pfinal.resize(1);
    CHECK_THROWS_AS(gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions, m_values, "QED"),
                    gra::AmplitudeFailure);
  }

  SECTION("missing lower forward leg") {
    gra::LORENTZSCALAR lts = MakeElasticEPAEvent(6500.0, 0.02, 0.03, -0.2, -0.3, 0.4, -0.7);
    lts.pfinal.resize(2);
    CHECK_THROWS_AS(gra::qed::PhotonSourceMatrixTransitions(lts, 2, transitions, m_values, "QED"),
                    gra::AmplitudeFailure);
  }

  SECTION("timelike upper transfer") {
    gra::LORENTZSCALAR lts = MakeElasticEPAEvent(6500.0, 0.02, 0.03, -0.2, -0.3, 0.4, -0.7);
    lts.model_cache =
        std::make_shared<gra::MModelCache>(gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
    lts.q1  = gra::M4Vec(0.0, 0.0, 0.0, 1.0);
    lts.t1  = lts.q1.M2();
    lts.qt1 = 0.0;
    REQUIRE(lts.q1.M2() > 0.0);
    CHECK_THROWS_AS(gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions, m_values, "QED"),
                    gra::AmplitudeFailure);
  }
}

TEST_CASE("DZ and LUX use collinear PDF and exact Moller-flux conventions", "[EPA][DZ][LUX][normalization]") {
  // Reference values divide out alpha/(2 pi) from the published DZ spectrum
  // [REFERENCE: Drees and Zeppenfeld, Phys. Rev. D39 (1989) 2536]
  const double inverse_prefactor = 2.0 * gra::math::PI / gra::qed::alpha_QED();
  REQUIRE(gra::flux::DZFlux(0.01) * inverse_prefactor == Approx(1416.2463108797074).epsilon(1e-13));
  REQUIRE(gra::flux::DZFlux(0.10) * inverse_prefactor == Approx(45.351233130774695).epsilon(1e-13));
  REQUIRE(gra::flux::DZFlux(0.50) * inverse_prefactor == Approx(0.18565645469573933).epsilon(1e-13));
  for (const double x : {0.01, 0.10, 0.50}) {
    REQUIRE(gra::flux::DZFlux(x) == Approx(NumericalDZFlux(x)).epsilon(2.0e-11));
  }
  const double boundary_x       = 1.0 - 1.0e-10;
  const double boundary_q2min   = std::pow(gra::PDG::mp * boundary_x, 2) / (1.0 - boundary_x);
  const double boundary_w       = 0.71 / (boundary_q2min + 0.71);
  const double boundary_leading = gra::qed::alpha_QED() / (2.0 * gra::math::PI * boundary_x) *
                                  (1.0 + std::pow(1.0 - boundary_x, 2)) * std::pow(boundary_w, 4) / 4.0;
  REQUIRE(gra::flux::DZFlux(boundary_x) > 0.0);
  REQUIRE(gra::flux::DZFlux(boundary_x) / boundary_leading == Approx(1.0).epsilon(2.0e-10));
  REQUIRE(gra::math::IsZero(gra::flux::DZFlux(0.0)));
  REQUIRE(gra::math::IsZero(gra::flux::DZFlux(1.0)));

  gra::LORENTZSCALAR dz     = MakeCollinearEPAEvent(3.0, 0.14, 0.23);
  dz.xi1                    = 0.01;
  dz.xi2                    = 0.02;
  dz.hamp                   = {std::complex<double>(0.8, -0.4)};
  const double beta         = gra::kinematics::beta12(dz.s, gra::PDG::mp, gra::PDG::mp);
  const double exact_factor = gra::flux::CollinearPhotonPhaseSpaceFactor(dz);
  REQUIRE(exact_factor == Approx(beta / (dz.x1 * dz.x2)).epsilon(1e-14));
  REQUIRE(std::abs(exact_factor - 1.0 / (dz.x1 * dz.x2)) > 1e-3);

  const double dz_flux1    = gra::flux::DZFlux(dz.x1);
  const double dz_flux2    = gra::flux::DZFlux(dz.x2);
  const auto   dz_hamp     = dz.hamp;
  const double dz_amp2     = 1.7;
  const double dz_weighted = gra::flux::ApplyDZfluxes(dz_amp2, dz);
  const double dz_weight   = dz_flux1 * dz_flux2 * exact_factor;
  REQUIRE(dz_weighted == Approx(dz_amp2 * dz_weight).epsilon(1e-14));
  REQUIRE(dz.id1 == gra::PDG::PDG_gamma);
  REQUIRE(dz.id2 == gra::PDG::PDG_gamma);
  REQUIRE(dz.pdf_xf1 == Approx(dz.x1 * dz_flux1).epsilon(1e-14));
  REQUIRE(dz.pdf_xf2 == Approx(dz.x2 * dz_flux2).epsilon(1e-14));
  REQUIRE(dz.muF == Approx(std::sqrt(dz.s_hat) / 2.0).epsilon(1e-14));
  REQUIRE(dz.muR == Approx(dz.muF));
  REQUIRE(dz.scalup == Approx(dz.muF));
  REQUIRE_FALSE(dz.exact_forward_photon_kinematics);
  RequireComplexNear(dz.hamp[0], dz_hamp[0] * std::sqrt(dz_weight), 1e-14);

  gra::LORENTZSCALAR invalid_dz = MakeCollinearEPAEvent(3.0, 0.14, 0.23);
  invalid_dz.hamp               = {1.0};
  REQUIRE_FALSE(std::isfinite(gra::flux::ApplyDZfluxes(std::numeric_limits<double>::quiet_NaN(), invalid_dz)));

  gra::LORENTZSCALAR lux     = MakeCollinearEPAEvent(6500.0, 0.08, 0.13);
  lux.xi1                    = 0.01;
  lux.xi2                    = 0.02;
  lux.LHAPDFSET              = "LUXqed17_plus_PDF4LHC15_nnlo_100";
  lux.hamp                   = {std::complex<double>(-0.3, 0.6)};
  const auto        lux_hamp = lux.hamp;
  gra::MLHAPDFStore pdf_store;
  const auto        pdf = pdf_store.GetPDF(lux.LHAPDFSET, 0);
  lux.GlobalPdfPtr      = pdf;
  const double Q2       = lux.s_hat / 4.0;
  lux.GlobalPdfPtr = pdf;
  REQUIRE(gra::flux::SetPhotonAlphaS(lux));
  REQUIRE(lux.alphaQCD == Approx(pdf->alphasQ2(Q2)).epsilon(1.0e-14));
  REQUIRE(lux.muF == Approx(std::sqrt(Q2)).epsilon(1.0e-14));
  REQUIRE(lux.muR == Approx(lux.muF));
  REQUIRE(lux.scalup == Approx(lux.muF));
  const double lux_flux1    = pdf->xfxQ2(gra::PDG::PDG_gamma, lux.x1, Q2) / lux.x1;
  const double lux_flux2    = pdf->xfxQ2(gra::PDG::PDG_gamma, lux.x2, Q2) / lux.x2;
  const double lux_factor   = gra::flux::CollinearPhotonPhaseSpaceFactor(lux);
  const double lux_amp2     = 2.1;
  const double lux_weighted = gra::flux::ApplyLUXfluxes(lux_amp2, lux);
  const double lux_weight   = lux_flux1 * lux_flux2 * lux_factor;
  REQUIRE(lux_weighted == Approx(lux_amp2 * lux_weight).epsilon(1e-13));
  REQUIRE(pdf->hasFlavor(gra::PDG::PDG_gamma));
  REQUIRE(lux.id1 == gra::PDG::PDG_gamma);
  REQUIRE(lux.id2 == gra::PDG::PDG_gamma);
  REQUIRE(lux.pdf_xf1 == Approx(lux.x1 * lux_flux1).epsilon(1e-13));
  REQUIRE(lux.pdf_xf2 == Approx(lux.x2 * lux_flux2).epsilon(1e-13));
  REQUIRE(lux.muF == Approx(std::sqrt(Q2)).epsilon(1e-14));
  REQUIRE(lux.muR == Approx(lux.muF));
  REQUIRE(lux.scalup == Approx(lux.muF));
  REQUIRE_FALSE(lux.exact_forward_photon_kinematics);
  RequireComplexNear(lux.hamp[0], lux_hamp[0] * std::sqrt(lux_weight), 1e-13);

  gra::LORENTZSCALAR below_grid = MakeCollinearEPAEvent(6500.0, 0.08, 0.13);
  below_grid.LHAPDFSET          = lux.LHAPDFSET;
  below_grid.GlobalPdfPtr       = pdf;
  const double q_min            = pdf->info().get_entry_as<double>("QMin");
  below_grid.s_hat              = q_min * q_min;
  below_grid.hamp               = {std::complex<double>(0.5, -0.2)};
  REQUIRE_FALSE(gra::flux::SetPhotonAlphaS(below_grid));
  REQUIRE(gra::math::IsZero(below_grid.alphaQCD));
  REQUIRE(gra::math::IsZero(gra::flux::ApplyLUXfluxes(1.0, below_grid)));
  RequireComplexNear(below_grid.hamp[0], 0.0, 0.0);


}

TEST_CASE("elastic scalar flux trace is distinct from one-direction EPA density", "[EPA][elastic][polarization]") {
  const double                          x                     = 0.17;
  const double                          pt                    = 0.63;
  const double                          t                     = -0.58;
  const ElasticFluxComponents           terms                 = ElasticComponents(x, t, pt);
  const gra::flux::TransversePhotonFlux density               = gra::flux::CohFluxTransverse(x, t, pt);
  const double                          scalar_trace          = terms.electric + 2.0 * terms.magnetic_per_state;
  const double                          parallel_density      = terms.electric + terms.magnetic_per_state;
  const double                          perpendicular_density = terms.magnetic_per_state;

  REQUIRE(terms.electric > 0.0);
  REQUIRE(terms.magnetic_per_state > 0.0);
  REQUIRE(gra::flux::CohFlux(x, t, pt) == Approx(scalar_trace).epsilon(1e-13));
  REQUIRE(density.parallel == Approx(parallel_density).epsilon(1e-13));
  REQUIRE(density.perpendicular == Approx(perpendicular_density).epsilon(1e-13));
  REQUIRE(density.Trace() == Approx(scalar_trace).epsilon(1e-13));
}

// Check the magnetic EPA density at the minimum photon virtuality
TEST_CASE("EPA magnetic densities retain their finite zero transverse momentum limit", "[EPA][flux][forward]") {
  const double xi = 0.17;
  const double mass2 = gra::math::pow2(gra::PDG::mp);
  const double elastic_q2 = xi * xi * mass2 / (1.0 - xi);
  const double elastic = 4.0 * gra::math::PI * gra::qed::alpha_QED() *
                         gra::math::pow2(gra::form::G_M(elastic_q2)) / (xi * mass2);
  const double remnant2 = 4.0;
  const double inelastic_q2 = (xi * (remnant2 - mass2) + xi * xi * mass2) / (1.0 - xi);
  const double xbj = inelastic_q2 / (inelastic_q2 + remnant2 - mass2);
  const double inelastic = 8.0 * gra::math::PI * gra::qed::alpha_QED() * xi *
                           gra::form::F1xQ2(xbj, inelastic_q2) / ((1.0 - xi) * inelastic_q2 * inelastic_q2);
  REQUIRE(elastic > 0.0);
  REQUIRE(inelastic > 0.0);
  for (const double pt : {0.0, 1.0e-170, 1.0e-155, 1.0e-12}) {
    CAPTURE(pt);
    const auto coherent = gra::flux::CohFluxTransverse(xi, -elastic_q2, pt);
    const auto incoherent = gra::flux::IncohFluxTransverse(xi, -inelastic_q2, pt, remnant2);
    CHECK(coherent.parallel == Approx(elastic).epsilon(1.0e-12));
    CHECK(coherent.perpendicular == Approx(elastic).epsilon(1.0e-12));
    CHECK(incoherent.parallel == Approx(inelastic).epsilon(1.0e-12));
    CHECK(incoherent.perpendicular == Approx(inelastic).epsilon(1.0e-12));
  }
}

// Compare F2 with the published CKMT fit rather than rounded intercept inputs
TEST_CASE("CKMT structure function reproduces the published low-Q2 fit", "[EPA][inelastic][structure-function]") {
  // [REFERENCE: Capella et al., arXiv:hep-ph/9405338, Eqs. (1), (3), (4), (6) and Fig. 2 caption]
  const std::array<std::array<double, 3>, 3> reference = {{{1.0e-4, 2.0, 0.6589329233481268},
                                                           {0.01, 1.0, 0.2902622242989286},
                                                           {0.4, 5.0, 0.2236194446201534}}};
  gra::form::ParamStore structure;
  structure.F2 = "CKMT";
  for (const auto &point : reference) {
    CAPTURE(point);
    CHECK(gra::form::F2xQ2(point[0], point[1], structure) == Approx(point[2]).epsilon(1.0e-12));
  }
  // At fixed photon energy F2 vanishes linearly in Q2 by electromagnetic gauge invariance
  const double two_m_nu = 100.0;
  const double q2 = 1.0e-7;
  const double first = gra::form::F2xQ2(q2 / two_m_nu, q2, structure) / q2;
  const double second = gra::form::F2xQ2(2.0 * q2 / two_m_nu, 2.0 * q2, structure) / (2.0 * q2);
  CHECK(first > 0.0);
  CHECK(first == Approx(second).epsilon(1.0e-5));
}

// Check separated structure functions and their inelastic photon flux
TEST_CASE("inelastic flux uses finite-Q2 F1 F2 and FL consistently", "[EPA][inelastic][structure-function]") {
  // Independent reference values from the official R1998 coefficient set
  // [REFERENCE: Abe et al., Phys. Lett. B452 (1999) 194, arXiv:hep-ex/9808028]
  REQUIRE(gra::form::RxQ2(0.12, 2.0) == Approx(0.29298796187101456).epsilon(1e-13));
  REQUIRE(gra::form::RxQ2(0.03, 1.5) == Approx(0.36770725679674027).epsilon(1e-13));
  REQUIRE(gra::form::RxQ2(0.40, 5.0) == Approx(0.12251399095833271).epsilon(1e-13));
  REQUIRE(gra::form::RxQ2(0.12, 0.2) == Approx(0.35594243692472544).epsilon(1e-13));
  REQUIRE(gra::form::RxQ2(0.12, 0.1) == Approx(0.17797121846236272).epsilon(1e-13));
  REQUIRE(gra::math::IsZero(gra::form::F2xQ2(0.0, 1.0)));
  REQUIRE(gra::math::IsZero(gra::form::F2xQ2(0.2, 0.0)));
  REQUIRE(gra::math::IsZero(gra::form::FLxQ2(1.1, 1.0)));
  REQUIRE(gra::math::IsZero(gra::form::F1xQ2(0.2, -1.0)));
  const double infinity = std::numeric_limits<double>::infinity();
  REQUIRE(gra::math::IsZero(gra::flux::CohFlux(0.1, -1.0, infinity)));
  REQUIRE(gra::math::IsZero(gra::flux::IncohFlux(0.1, -1.0, infinity, 4.0)));

  const double                          x             = 0.12;
  const double                          pt            = 0.65;
  const double                          remnant_mass2 = 2.4 * 2.4;
  const double                          mp2           = gra::PDG::mp * gra::PDG::mp;
  const double                          Q2            = (pt * pt + x * (remnant_mass2 - mp2) + x * x * mp2) / (1.0 - x);
  const double                          t             = -Q2;
  const InelasticFluxComponents         terms         = InelasticComponents(x, t, pt, remnant_mass2);
  const gra::flux::TransversePhotonFlux density       = gra::flux::IncohFluxTransverse(x, t, pt, remnant_mass2);
  const double                          scalar_trace  = terms.electric + 2.0 * terms.magnetic_per_state;

  REQUIRE(terms.xbj > 0.0);
  REQUIRE(terms.xbj < 1.0);
  REQUIRE(terms.Q2 == Approx(Q2).epsilon(1e-14));
  REQUIRE(gra::flux::IncohFlux(x, t, pt, remnant_mass2) == Approx(scalar_trace).epsilon(1e-13));
  REQUIRE(density.parallel == Approx(terms.electric + terms.magnetic_per_state).epsilon(1e-13));
  REQUIRE(density.perpendicular == Approx(terms.magnetic_per_state).epsilon(1e-13));

  const double f2          = gra::form::F2xQ2(terms.xbj, terms.Q2);
  const double fl          = gra::form::FLxQ2(terms.xbj, terms.Q2);
  const double code_2xf1   = 2.0 * terms.xbj * gra::form::F1xQ2(terms.xbj, terms.Q2);
  const double target_mass = 4.0 * terms.xbj * terms.xbj * mp2 / terms.Q2;
  const double ratio       = gra::form::RxQ2(terms.xbj, terms.Q2);

  REQUIRE(ratio > 0.0);
  REQUIRE(fl > 0.0);
  REQUIRE(code_2xf1 + fl == Approx((1.0 + target_mass) * f2).epsilon(1e-14));
  REQUIRE(fl / code_2xf1 == Approx(ratio).epsilon(1e-14));
  REQUIRE(gra::form::RxQ2(terms.xbj, 0.1) == Approx(0.5 * gra::form::RxQ2(terms.xbj, 0.2)).epsilon(1e-14));
}

TEST_CASE(
    "covariant MG5 embedding preserves closure and contracts normalized "
    "sources",
    "[EPA][MG2GRA][covariant]") {
  CovariantMG5Point point       = MakeCovariantMG5Point();
  const auto        input_final = point.final;
  gra::M4Vec        k1;
  gra::M4Vec        k2;

  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(point.lts, point.final, k1, k2));

  for (std::size_t i = 0; i < point.final.size(); ++i) { RequireFourVectorNear(point.final[i], input_final[i], 1e-14); }
  const gra::M4Vec hard = point.final[0] + point.final[1];
  RequireFourVectorNear(k1 + k2, hard, 1e-12);
  REQUIRE(std::abs(k1.M2()) < 1e-11);
  REQUIRE(std::abs(k2.M2()) < 1e-11);
  REQUIRE((k1 + k2).M2() == Approx(hard.M2()).margin(1e-11));

  const auto upper               = gra::mg5helas::TransverseHelicityCoefficients(point.lts.pbeam1, k1, k2, true);
  const auto lower               = gra::mg5helas::TransverseHelicityCoefficients(point.lts.pbeam2, k2, k1, true);
  const auto upper_perpendicular = PerpendicularSource(upper);
  const auto lower_perpendicular = PerpendicularSource(lower);
  const auto upper_helas         = HELASSourceCoefficients(point.lts.pbeam1, k1, k2);
  const auto lower_helas         = HELASSourceCoefficients(point.lts.pbeam2, k2, k1);
  REQUIRE(std::norm(upper[0]) + std::norm(upper[1]) == Approx(1.0).epsilon(1e-13));
  REQUIRE(std::norm(lower[0]) + std::norm(lower[1]) == Approx(1.0).epsilon(1e-13));
  REQUIRE(std::norm(upper_perpendicular[0]) + std::norm(upper_perpendicular[1]) == Approx(1.0).epsilon(1e-13));
  REQUIRE(std::norm(lower_perpendicular[0]) + std::norm(lower_perpendicular[1]) == Approx(1.0).epsilon(1e-13));
  RequireComplexNear(std::conj(upper[0]) * upper_perpendicular[0] + std::conj(upper[1]) * upper_perpendicular[1], 0.0,
                     1e-14);
  RequireComplexNear(std::conj(lower[0]) * lower_perpendicular[0] + std::conj(lower[1]) * lower_perpendicular[1], 0.0,
                     1e-14);
  for (std::size_t index = 0; index < 2; ++index) {
    RequireComplexNear(upper[index], upper_helas[index], 1e-14);
    RequireComplexNear(lower[index], lower_helas[index], 1e-14);
  }

  const std::array<std::complex<double>, 4> hard_amplitudes = {
      std::complex<double>(1.0, 0.2), std::complex<double>(-0.3, 0.7), std::complex<double>(0.4, -0.5),
      std::complex<double>(-0.8, 0.1)};
  const std::array<std::array<int, 2>, 4>       helicities = {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, 1},
                                                              std::array<int, 2>{1, -1}, std::array<int, 2>{1, 1}};
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (std::size_t i = 0; i < hard_amplitudes.size(); ++i) {
    components.push_back({helicities[i], {}, 0, hard_amplitudes[i], {}});
  }

  const gra::flux::TransversePhotonFlux density1 =
      gra::flux::ForwardPhotonFluxTransverse(gra::ResolveForwardLegState(point.lts, gra::ForwardBeamLeg::Upper));
  const gra::flux::TransversePhotonFlux density2 =
      gra::flux::ForwardPhotonFluxTransverse(gra::ResolveForwardLegState(point.lts, gra::ForwardBeamLeg::Lower));
  const std::array<std::array<std::complex<double>, 2>, 2> sources1               = {upper, upper_perpendicular};
  const std::array<std::array<std::complex<double>, 2>, 2> sources2               = {lower, lower_perpendicular};
  double                                                   fixed_helicity_average = 0.0;
  for (const auto &amplitude : hard_amplitudes) { fixed_helicity_average += std::norm(amplitude) / 4.0; }
  double rotated_unpolarized_average = 0.0;
  for (const auto &source1 : sources1) {
    for (const auto &source2 : sources2) {
      rotated_unpolarized_average +=
          std::norm(gra::mg5helas::ContractHelicitySources(hard_amplitudes, source1, source2)) / 4.0;
    }
  }
  REQUIRE(rotated_unpolarized_average == Approx(fixed_helicity_average).epsilon(1e-13));

  const std::array<double, 2> weights1   = {density1.parallel / density1.Trace(),
                                            density1.perpendicular / density1.Trace()};
  const std::array<double, 2> weights2   = {density2.parallel / density2.Trace(),
                                            density2.perpendicular / density2.Trace()};
  const auto                  contracted = gra::mg5helas::ContractEPAPhotonSources(point.lts, components, k1, k2);
  REQUIRE(contracted.size() == 4);
  std::size_t source_index        = 0;
  double      spin_averaged_norm  = 0.0;
  double      density_contraction = 0.0;
  for (std::size_t mode1 = 0; mode1 < 2; ++mode1) {
    for (std::size_t mode2 = 0; mode2 < 2; ++mode2) {
      const std::complex<double> hard_source =
          gra::mg5helas::ContractHelicitySources(hard_amplitudes, sources1[mode1], sources2[mode2]);
      const std::complex<double> expected = 2.0 * std::sqrt(weights1[mode1] * weights2[mode2]) * hard_source;
      RequireComplexNear(contracted[source_index], expected, 1e-13);
      spin_averaged_norm += std::norm(contracted[source_index]) / 4.0;
      density_contraction += weights1[mode1] * weights2[mode2] * std::norm(hard_source);
      ++source_index;
    }
  }
  REQUIRE(spin_averaged_norm == Approx(density_contraction).epsilon(1e-13));

}

// Check the Cartesian EPA tensor in the explicit HELAS beam convention
// [REFERENCE: Harland-Lang, Khoze and Ryskin, Eur. Phys. J. C 79 (2019) 39, arXiv:1810.06567, Eq. (15)]
TEST_CASE("EPA helicity sources reproduce the Cartesian photon tensor", "[EPA][MG2GRA][helicity][covariance]") {
  const gra::M4Vec                          k1(0.0, 0.0, 20.0, 20.0);
  const gra::M4Vec                          k2(0.0, 0.0, -20.0, 20.0);
  const gra::M4Vec                          q1(0.73, -0.41, 0.0, 0.0);
  const gra::M4Vec                          q2(-0.28, 0.64, 0.0, 0.0);
  const std::array<std::complex<double>, 4> hard = {std::complex<double>(0.37, -0.22),
                                                    std::complex<double>(-0.61, 0.19), std::complex<double>(0.48, 0.53),
                                                    std::complex<double>(-0.14, -0.72)};

  const double               dot        = q1.Px() * q2.Px() + q1.Py() * q2.Py();
  const double               cross      = q1.Px() * q2.Py() - q1.Py() * q2.Px();
  const double               quadrupole = q1.Px() * q2.Px() - q1.Py() * q2.Py();
  const double               mixed      = q1.Px() * q2.Py() + q1.Py() * q2.Px();
  const std::complex<double> expected   = -0.5 * dot * (hard[0] + hard[3]) +
                                        0.5 * gra::math::zi * cross * (hard[0] - hard[3]) +
                                        0.5 * std::complex<double>(quadrupole, mixed) * hard[1] +
                                        0.5 * std::complex<double>(quadrupole, -mixed) * hard[2];
  RequireComplexNear(gra::mg5helas::ContractTransverseSources(hard, q1, q2, k1, k2), expected, 2.0e-15);

  // A common transverse rotation is cancelled by the incoming hard-tensor phase
  const double angle  = 0.83;
  const auto   rotate = [angle](const gra::M4Vec &q) {
    return gra::M4Vec(std::cos(angle) * q.Px() - std::sin(angle) * q.Py(),
                        std::sin(angle) * q.Px() + std::cos(angle) * q.Py(), q.Pz(), q.E());
  };
  auto rotated_hard = hard;
  rotated_hard[1] *= std::polar(1.0, -2.0 * angle);
  rotated_hard[2] *= std::polar(1.0, 2.0 * angle);
  const auto rotated = gra::mg5helas::ContractTransverseSources(rotated_hard, rotate(q1), rotate(q2), k1, k2);
  RequireComplexNear(rotated, expected, 2.0e-15);
}

TEST_CASE("coherent MG5 EPA sources retain the emitter charge phase", "[EPA][MG2GRA][phase]") {
  CovariantMG5Point point  = MakeCovariantMG5Point();
  point.lts.beam1.pdg      = gra::PDG::PDG_p;
  point.lts.beam1.chargeX3 = 3;
  point.lts.beam2          = point.lts.beam1;

  gra::M4Vec k1;
  gra::M4Vec k2;
  REQUIRE(gra::mg5helas::PrepareOnShellKinematics(point.lts, point.final, k1, k2));

  const std::array<std::array<int, 2>, 4>   helicities = {std::array<int, 2>{-1, -1}, std::array<int, 2>{-1, 1},
                                                          std::array<int, 2>{1, -1}, std::array<int, 2>{1, 1}};
  const std::array<std::complex<double>, 4> hard = {std::complex<double>(0.7, -0.2), std::complex<double>(-0.4, 0.9),
                                                    std::complex<double>(0.3, 0.5), std::complex<double>(-0.8, -0.1)};
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (const auto &i : indices(helicities)) { components.push_back({helicities[i], {}, 0, hard[i], {}}); }

  const auto evaluate = [&](int upper_charge, int lower_charge) {
    auto lts           = point.lts;
    lts.beam1.pdg      = upper_charge > 0 ? gra::PDG::PDG_p : -gra::PDG::PDG_p;
    lts.beam1.chargeX3 = 3 * upper_charge;
    lts.beam2.pdg      = lower_charge > 0 ? gra::PDG::PDG_p : -gra::PDG::PDG_p;
    lts.beam2.chargeX3 = 3 * lower_charge;
    return gra::mg5helas::ContractEPAPhotonSources(lts, components, k1, k2);
  };

  const auto pp        = evaluate(1, 1);
  const auto pbar_p    = evaluate(-1, 1);
  const auto p_pbar    = evaluate(1, -1);
  const auto pbar_pbar = evaluate(-1, -1);
  REQUIRE_FALSE(pp.empty());
  REQUIRE(pbar_p.size() == pp.size());
  REQUIRE(p_pbar.size() == pp.size());
  REQUIRE(pbar_pbar.size() == pp.size());
  for (const auto &i : indices(pp)) {
    RequireComplexNear(pbar_p[i], -pp[i], 1.0e-13);
    RequireComplexNear(p_pbar[i], -pp[i], 1.0e-13);
    RequireComplexNear(pbar_pbar[i], pp[i], 1.0e-13);
  }
}
