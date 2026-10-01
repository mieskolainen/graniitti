// Physics regression tests for massive and unweighted phase space
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>

#include <algorithm>
#include <cmath>
#include <limits>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Sampling/MRandom.h"

namespace {

using gra::aux::indices;
using gra::kinematics::DecayMomentum;
using gra::kinematics::dPhi2;

// Keep null daughter masses stable under roundoff while rejecting spacelike momenta
TEST_CASE("physical invariant masses tolerate null roundoff without admitting spacelike states",
          "[Kinematics][Mass][Roundoff]") {
  for (const double energy : {1e-100, 1e-8, 1.0, 100.0, 1e6, 1e100}) {
    for (const double direction : {0.0, std::numeric_limits<double>::infinity()}) {
      const gra::M4Vec null(0.0, 0.0, std::nextafter(energy, direction), energy);
      double mass = -1.0;
      REQUIRE(null.PhysicalMass(mass));
      REQUIRE(mass == Approx(0.0).margin(1e-14));
      REQUIRE_FALSE(gra::M4Vec(0.0, 0.0, 1.01 * energy, energy).PhysicalMass(mass));
    }
  }
  for (const double scale : {1e-100, 1e-8, 1.0, 1e100}) {
    double mass = 0.0;
    REQUIRE(gra::M4Vec(0.0, 0.0, 0.0, scale).PhysicalMass(mass));
    REQUIRE(mass / scale == Approx(1.0).epsilon(1e-14));
    REQUIRE_FALSE(gra::M4Vec(scale, 0.0, 0.0, 0.0).PhysicalMass(mass));
  }
  gra::M4Vec massive(0.0, 0.0, 4.0, 5.0);
  massive.RotateX(0.43);
  massive.RotateY(-0.71);
  massive = massive.LorentzBoost({0.21, -0.13, 0.31});
  double mass = 0.0;
  REQUIRE(massive.PhysicalMass(mass));
  REQUIRE(mass == Approx(3.0).epsilon(1e-12));
  REQUIRE_FALSE(gra::M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN()).PhysicalMass(mass));
}

// Supply fixed unit coordinates through the phase-space random interface
struct PhaseSpaceRandom {
  std::vector<double> unit;
  std::size_t index = 0;

  // Map the next unit coordinate to the requested interval
  double U(double lower, double upper) { return std::lerp(lower, upper, unit.at(index++)); }
};

// Construct equal-energy massless tetrahedral directions with a common rotation
PhaseSpaceRandom TetrahedronRandom(double theta, double phi, int parity) {
  const double a = 1.0 / std::sqrt(3.0);
  std::vector<gra::M4Vec> direction = {{a, a, a, 1.0}, {a, -a, -a, 1.0},
                                      {-a, a, -a, 1.0}, {-a, -a, a, 1.0}};
  PhaseSpaceRandom random;
  for (auto &p : direction) {
    p.Rotate(theta, phi);
    p.SetPxPyPz(parity * p.Px(), parity * p.Py(), parity * p.Pz());
    const double azimuth = p.Phi() < 0.0 ? p.Phi() + 2.0 * gra::math::PI : p.Phi();
    random.unit.insert(random.unit.end(), {std::exp(-0.5), std::exp(-0.5),
                                          0.5 * (1.0 + p.Pz()), azimuth / (2.0 * gra::math::PI)});
  }
  return random;
}

// Solve the symmetric massive final state independently by energy bisection
long double SymmetricMomentum(double mass, const std::vector<double> &m) {
  long double lower = 0.0L;
  long double upper = static_cast<long double>(mass) / m.size();
  for (unsigned int i = 0; i < 128; ++i) {
    const long double p = (lower + upper) / 2.0L;
    long double energy = 0.0L;
    for (const double daughter : m) { energy += std::hypot(static_cast<long double>(daughter), p); }
    if (energy > mass) {
      upper = p;
    } else {
      lower = p;
    }
  }
  return (lower + upper) / 2.0L;
}

// Compute the RAMBO Jacobian for equal-energy massless seed momenta
// J = xi^(3n-5) M product_i(q_i/E_i) / sum_i(q_i^2/E_i)
double SymmetricRamboWeight(double mass, const std::vector<double> &m, long double p) {
  const long double q = static_cast<long double>(mass) / m.size();
  long double product = 1.0L;
  long double derivative = 0.0L;
  for (const double daughter : m) {
    const long double energy = std::hypot(static_cast<long double>(daughter), p);
    product *= q / energy;
    derivative += q * q / energy;
  }
  const long double jacobian = std::pow(p / q, 3 * static_cast<int>(m.size()) - 5) * mass * product / derivative;
  return static_cast<double>(jacobian * gra::kinematics::PSnMassless(mass * mass, m.size()));
}

// Integrate recursive two-body measures with an optional cut on the first pair mass
// Phi_4 = int dmu^2/(2pi) dnu^2/(2pi) Phi_2(M,mu,m3) Phi_2(mu,nu,m2) Phi_2(nu,m0,m1)
double FourBodyVolume(double mass, const std::vector<double> &m, double pair_max) {
  const auto [outer, outer_weight] = gra::math::GaussLegendreRule(96, m[0] + m[1] + m[2], mass - m[3]);
  double integral = 0.0;
  for (const auto &i : indices(outer)) {
    const double upper = std::min(pair_max, outer[i] - m[2]);
    if (!(upper > m[0] + m[1])) { continue; }
    const auto [inner, inner_weight] = gra::math::GaussLegendreRule(96, m[0] + m[1], upper);
    double partial = 0.0;
    for (const auto &j : indices(inner)) {
      partial += inner_weight[j] * inner[j] / gra::math::PI *
                 dPhi2(outer[i], DecayMomentum(outer[i], inner[j], m[2])) *
                 dPhi2(inner[j], DecayMomentum(inner[j], m[0], m[1]));
    }
    integral += outer_weight[i] * outer[i] / gra::math::PI *
                dPhi2(mass, DecayMomentum(mass, outer[i], m[3])) * partial;
  }
  return integral;
}

} // namespace

// Check rejection sampling and mass shells for every location of the heavy daughter
TEST_CASE("N-body rejection remains active for unequal daughter masses", "[PhaseSpace][NBody]") {
  const double mass = 1.6;
  const gra::M4Vec mother(0.0, 0.0, 0.0, mass);
  std::vector<double> m = {0.9, 0.1, 0.1, 0.1};
  for (std::size_t order = 0; order < m.size(); ++order) {
    CAPTURE(order);
    PhaseSpaceRandom random{{0.25, 0.75, 0.99, 0.25, 0.75, 1e-6, 0.2, 0.4, 0.7, 0.1, 0.6, 0.8}};
    std::vector<gra::M4Vec> p;
    const auto weight = gra::kinematics::NBodyPhaseSpace(mother, mass, m, p, true, random);
    REQUIRE(weight.GetN() == Approx(2.0));
    REQUIRE(weight.Integral() > 0.0);
    gra::M4Vec total;
    for (const auto &i : indices(p)) {
      CHECK(p[i].M2() == Approx(m[i] * m[i]).margin(1e-13));
      CHECK(p[i].E() > 0.0);
      total += p[i];
    }
    CHECK(gra::math::CheckEMC(total - mother, 1e-13));
    std::rotate(m.begin(), m.begin() + 1, m.end());
  }
}

// Compare unweighted pair-mass distributions and integrated volumes with recursive quadrature
TEST_CASE("Unequal-mass N-body phase space agrees with recursive integration", "[PhaseSpace][NBody]") {
  const double mass = 1.6;
  const gra::M4Vec mother(0.0, 0.0, 0.0, mass);
  for (const auto &m : {std::vector<double>{0.9, 0.1, 0.1, 0.1}, std::vector<double>{0.1, 0.1, 0.1, 0.9}}) {
    const double pair_max = 0.5 * (m[0] + m[1] + mass - m[2] - m[3]);
    const double volume = FourBodyVolume(mass, m, mass);
    const double probability = FourBodyVolume(mass, m, pair_max) / volume;
    gra::MRandom random;
    random.SetSeed(291847);
    gra::kinematics::MCW total;
    std::size_t below = 0;
    const std::size_t samples = 20000;
    std::vector<gra::M4Vec> p;
    for (std::size_t event = 0; event < samples; ++event) {
      const auto weight = gra::kinematics::NBodyPhaseSpace(mother, mass, m, p, true, random);
      REQUIRE(weight.Integral() > 0.0);
      total += weight;
      if ((p[0] + p[1]).M() < pair_max) { ++below; }
    }
    CHECK(total.GetN() > samples);
    CHECK(std::abs(total.Integral() / volume - 1.0) < 6.0 * total.IntegralError() / volume + 2e-4);
    const double error = std::sqrt(probability * (1.0 - probability) / samples);
    CHECK(std::abs(static_cast<double>(below) / samples - probability) < 6.0 * error + 2e-4);
  }
}

// Compare weighted RAMBO with the massive recursive measure close to threshold
TEST_CASE("Massive RAMBO threshold volume and pair masses match quadrature", "[PhaseSpace][Rambo]") {
  const std::vector<double> m = {0.1, 0.2, 0.3, 0.4};
  const double mass = 1.000001;
  const gra::M4Vec mother(0.0, 0.0, 0.0, mass);
  const double pair_max = 0.5 * (m[0] + m[1] + mass - m[2] - m[3]);
  const double volume = FourBodyVolume(mass, m, mass);
  const double probability = FourBodyVolume(mass, m, pair_max) / volume;
  gra::MRandom random;
  random.SetSeed(371629);
  gra::kinematics::MCW total;
  gra::kinematics::MCW below;
  std::vector<gra::M4Vec> p;
  for (std::size_t event = 0; event < 20000; ++event) {
    const auto weight = gra::kinematics::RamboMassive(mother, mass, m, p, random);
    REQUIRE(weight.Integral() > 0.0);
    total += weight;
    below.Push((p[0] + p[1]).M() < pair_max ? weight.Integral() : 0.0);
  }
  CHECK(std::abs(total.Integral() / volume - 1.0) < 6.0 * total.IntegralError() / volume + 2e-4);
  const double error = (below.IntegralError() + probability * total.IntegralError()) / total.Integral();
  CHECK(std::abs(below.Integral() / total.Integral() - probability) < 6.0 * error + 2e-4);
}

// Check threshold weights, mass shells, Lorentz covariance, rotations and parity
TEST_CASE("Massive RAMBO resolves kinetic energy near threshold", "[PhaseSpace][Rambo][Lorentz]") {
  const std::vector<double> m = {0.1, 0.2, 0.3, 0.4};
  for (const double excess : {1e-3, 1e-6, 1e-10, 1e-13}) {
    CAPTURE(excess);
    const double mass = 1.0 + excess;
    const gra::M4Vec mother(0.0, 0.0, 0.0, mass);
    const long double momentum = SymmetricMomentum(mass, m);
    const double expected = SymmetricRamboWeight(mass, m, momentum);
    auto random = TetrahedronRandom(0.0, 0.0, 1);
    std::vector<gra::M4Vec> p;
    const auto weight = gra::kinematics::RamboMassive(mother, mass, m, p, random);
    REQUIRE(weight.Integral() > 0.0);
    CHECK(weight.Integral() / expected == Approx(1.0).epsilon(2e-5));
    long double kinetic = 0.0L;
    long double available = mass;
    gra::M4Vec total;
    for (const auto &i : indices(p)) {
      CHECK(p[i].M2() == Approx(m[i] * m[i]).margin(1e-15));
      CHECK(p[i].P3mod() / static_cast<double>(momentum) == Approx(1.0).epsilon(3e-6));
      const long double k = p[i].P3mod();
      kinetic += k * (k / (static_cast<long double>(p[i].E()) + m[i]));
      available -= static_cast<long double>(m[i]);
      total += p[i];
    }
    CHECK(std::abs(kinetic / available - 1.0L) < 1e-10L);
    CHECK(gra::math::CheckEMC(total - mother, 2e-15));

    gra::M4Vec lab;
    lab.SetPxPyPzM(0.3, -0.2, 0.4, mass);
    for (const int parity : {-1, 1}) {
      auto transformed_random = TetrahedronRandom(0.47, -0.31, parity);
      std::vector<gra::M4Vec> transformed;
      const auto transformed_weight = gra::kinematics::RamboMassive(lab, mass, m, transformed, transformed_random);
      REQUIRE(transformed_weight.Integral() > 0.0);
      CHECK(transformed_weight.Integral() / weight.Integral() == Approx(1.0).epsilon(1e-11));
      for (const auto &i : indices(p)) {
        gra::M4Vec expected_p = p[i];
        expected_p.Rotate(0.47, -0.31);
        expected_p.SetPxPyPz(parity * expected_p.Px(), parity * expected_p.Py(), parity * expected_p.Pz());
        gra::kinematics::LorentzBoost(lab, mass, expected_p, 1);
        CHECK(gra::math::CheckEMC(transformed[i] - expected_p, 1e-13));
      }
    }
  }
}
