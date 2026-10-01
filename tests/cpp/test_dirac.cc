// Unit test for MDirac class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <catch.hpp>
#include <complex>
#include <random>
#include <utility>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Spin/MHelicity.h"

using namespace gra;

// Require a spinor to satisfy its Dirac equation and covariant normalization
void RequireDiracSpinorForTest(const gra::MDirac &dirac,
                               const gra::M4Vec &momentum,
                               const gra::MDirac::Spinor &spinor,
                               bool antiparticle) {
  const double mass = std::sqrt(std::max(0.0, momentum.M2()));
  const auto equation =
      dirac.FSlash(momentum) + dirac.I4 * (antiparticle ? mass : -mass);
  const auto residual = equation * spinor;
  const double scale = std::max(1.0, momentum.E());
  for (const auto &component : residual) {
    REQUIRE(std::abs(component) <= 2e-12 * scale);
  }

  const auto adjoint = dirac.Bar(spinor);
  const std::complex<double> norm = gra::BilinearProduct(adjoint, spinor);
  const double expected = antiparticle ? -2.0 * mass : 2.0 * mass;
  REQUIRE(norm.real() == Approx(expected).margin(2e-12 * scale));
  REQUIRE(norm.imag() == Approx(0.0).margin(2e-12 * scale));
}

// Require both raised and lowered gamma matrices to obey the Clifford algebra
void RequireCliffordAlgebra(const gra::MDirac &dirac, double tolerance) {
  for (const auto mu : dirac.LI) {
    for (const auto nu : dirac.LI) {
      const auto expected = dirac.I4 * (2.0 * dirac.g[mu][nu]);
      const auto lower = dirac.gamma_lo[mu] * dirac.gamma_lo[nu] +
                         dirac.gamma_lo[nu] * dirac.gamma_lo[mu];
      const auto upper = dirac.gamma_up[mu] * dirac.gamma_up[nu] +
                         dirac.gamma_up[nu] * dirac.gamma_up[mu];
      CAPTURE(mu, nu);
      REQUIRE((lower - expected).FrobNorm() < tolerance);
      REQUIRE((upper - expected).FrobNorm() < tolerance);
    }
  }
}

// Construct one physical spinor through the selected public representation
gra::MDirac::Spinor PhysicalSpinor(const gra::MDirac &dirac,
                                   const gra::M4Vec &momentum, int helicity,
                                   bool antiparticle, bool chiral,
                                   const std::string &mode) {
  if (mode == "gauge") {
    return antiparticle ? dirac.vGauge(momentum, helicity)
                        : dirac.uGauge(momentum, helicity);
  }
  if (mode == "spin") {
    return antiparticle ? dirac.vDirac(momentum, helicity)
                        : dirac.uDirac(momentum, helicity);
  }
  if (chiral) {
    return antiparticle ? dirac.vHelChiral(momentum, helicity)
                        : dirac.uHelChiral(momentum, helicity);
  }
  return antiparticle ? dirac.vHelDirac(momentum, helicity)
                      : dirac.uHelDirac(momentum, helicity);
}

// Require covariant normalization and the complete spinor projector identity
void RequireSpinorCompleteness(const gra::MDirac &dirac,
                               const gra::M4Vec &momentum,
                               bool antiparticle, bool chiral,
                               const std::string &mode, double tolerance) {
  gra::MMatrix<std::complex<double>> sum(4, 4, 0.0);
  const double sign = antiparticle ? -1.0 : 1.0;
  for (const int helicity : {-1, 1}) {
    const auto spinor =
        PhysicalSpinor(dirac, momentum, helicity, antiparticle, chiral, mode);
    const auto adjoint = dirac.Bar(spinor);
    const auto norm = gra::BilinearProduct(adjoint, spinor);
    const double expected = sign * 2.0 * momentum.M();
    CAPTURE(chiral, mode, antiparticle, helicity);
    REQUIRE(norm.real() == Approx(expected).epsilon(5e-3).margin(1e-10));
    REQUIRE(norm.imag() == Approx(0.0).margin(1e-10));
    sum.AddOuterProduct(spinor, adjoint, 1.0);
  }
  const auto projector = dirac.FSlash(momentum) + dirac.I4 * momentum.M() * sign;
  REQUIRE((sum - projector).FrobNorm() < tolerance);
}

// Require the massive spin-one polarization completeness identity
void RequireMassiveSpin1Completeness(const gra::MDirac &dirac,
                                     const gra::M4Vec &momentum,
                                     double tolerance) {
  REQUIRE(momentum.M2() > 0.0);
  for (const auto mu : dirac.LI) {
    for (const auto nu : dirac.LI) {
      std::complex<double> sum = 0.0;
      for (const int helicity : {-1, 0, 1}) {
        const auto polarization =
            dirac.EpsMassiveSpin1(momentum, helicity);
        sum += polarization(mu) * std::conj(polarization(nu));
      }
      const double expected =
          -dirac.g[mu][nu] + momentum[mu] * momentum[nu] / momentum.M2();
      CAPTURE(mu, nu, sum, expected);
      REQUIRE(std::abs(sum - expected) < tolerance);
    }
  }
}

// Dirac algebra tests
TEST_CASE("MDirac gamma matrices obey the Clifford algebra", "[MDirac]") {
  const double EPS = 1e-5;
  const MDirac dirac;
  const MDirac second;
  const MDirac chiral("CHIRAL");
  RequireCliffordAlgebra(dirac, EPS);
  REQUIRE(&dirac.gamma_up == &second.gamma_up);
  REQUIRE(&dirac.sigma_up == &second.sigma_up);
  REQUIRE(&dirac.PR() == &second.PR());
  REQUIRE(&dirac.gamma_up != &chiral.gamma_up);
  REQUIRE(&dirac.I4 == &chiral.I4);
  REQUIRE(&dirac.J_operator(1) == &chiral.J_operator(1));
}

// Check that the same physical helicity spinors are related by the basis map
TEST_CASE("MDirac Weyl and Dirac helicity spinors obey the basis transformation",
          "[MDirac][basis]") {
  const MDirac dirac("DIRAC");
  const MDirac chiral("CHIRAL");
  const std::array<M4Vec, 3> momenta = {
      M4Vec(0.0, 0.0, 0.0, 0.7),
      M4Vec(0.3, -0.4, 1.2, std::sqrt(0.7 * 0.7 + 1.69)),
      M4Vec(-1.1, 0.6, -0.2, std::sqrt(1.3 * 1.3 + 1.61))};

  for (const auto &momentum : momenta) {
    for (const int helicity : {-1, 1}) {
      const auto u_expected = dirac.uHelDirac(momentum, helicity);
      const auto v_expected = dirac.vHelDirac(momentum, helicity);
      const auto u_transformed =
          dirac.S_basis * chiral.uHelChiral(momentum, helicity);
      const auto v_transformed =
          dirac.S_basis * chiral.vHelChiral(momentum, helicity);
      CAPTURE(momentum, helicity);
      for (std::size_t component = 0; component < 4; ++component) {
        REQUIRE(std::abs(u_transformed[component] - u_expected[component]) <
                2e-12);
        REQUIRE(std::abs(v_transformed[component] - v_expected[component]) <
                2e-12);
      }
    }
  }
}

// Check the elastic proton current against the Gordon identity with q = p' - p
TEST_CASE("MDirac elastic proton current uses the outgoing-minus-incoming "
          "Pauli momentum",
          "[MDirac][photon]") {
  MDirac dirac("DIRAC");
  const double mp = gra::PDG::mp;
  const gra::M4Vec incoming(0.0, 0.0, 3.0, std::sqrt(mp * mp + 9.0));
  const gra::M4Vec outgoing(
      0.25, -0.15, 2.6,
      std::sqrt(mp * mp + 0.25 * 0.25 + 0.15 * 0.15 + 2.6 * 2.6));
  const gra::M4Vec vertex_q = outgoing - incoming;
  const gra::M4Vec momentum_sum = outgoing + incoming;
  const double f1 = gra::form::F1(vertex_q.M2());
  const double f2 = gra::form::F2(vertex_q.M2());

  for (const double lambda_in : {-0.5, 0.5}) {
    for (const double lambda_out : {-0.5, 0.5}) {
      const int h_in = static_cast<int>(2.0 * lambda_in);
      const int h_out = static_cast<int>(2.0 * lambda_out);
      const auto u = dirac.uHelDirac(incoming, h_in);
      const auto ubar = dirac.Bar(dirac.uHelDirac(outgoing, h_out));
      const auto current = dirac.DiracPauliCurrent(
          incoming, outgoing, vertex_q,
          gra::MDirac::FermionKind::Particle, lambda_in, lambda_out, mp, f1,
          f2);
      const std::complex<double> scalar = gra::BilinearProduct(ubar, u);

      for (const auto mu : dirac.LI) {
        const std::complex<double> gamma =
            dirac.gamma_lo[mu].BilinearForm(ubar, u);
        const std::complex<double> expected =
            (f1 + f2) * gamma - f2 * (momentum_sum % mu) / (2.0 * mp) * scalar;
        REQUIRE(std::abs(current[mu] - expected) == Approx(0.0).margin(2e-12));
      }
    }
  }
}

// Check the antiparticle current against its Gordon identity and Ward identity
TEST_CASE("MDirac elastic antiparticle current follows physical fermion flow",
          "[MDirac][photon][antiparticle]") {
  MDirac dirac("DIRAC");
  const double mass = gra::PDG::mp;
  const gra::M4Vec incoming(0.0, 0.0, 3.0, std::sqrt(mass * mass + 9.0));
  const gra::M4Vec outgoing(
      0.25, -0.15, 2.6,
      std::sqrt(mass * mass + 0.25 * 0.25 + 0.15 * 0.15 + 2.6 * 2.6));
  const gra::M4Vec vertex_q = outgoing - incoming;
  const gra::M4Vec momentum_sum = outgoing + incoming;
  constexpr double f1 = 0.83;
  constexpr double f2 = 1.47;

  for (const double lambda_in : {-0.5, 0.5}) {
    for (const double lambda_out : {-0.5, 0.5}) {
      const int h_in = static_cast<int>(2.0 * lambda_in);
      const int h_out = static_cast<int>(2.0 * lambda_out);
      const auto vbar = dirac.Bar(dirac.vHelDirac(incoming, h_in));
      const auto v = dirac.vHelDirac(outgoing, h_out);
      const auto current = dirac.DiracPauliCurrent(
          incoming, outgoing, vertex_q,
          gra::MDirac::FermionKind::Antiparticle, lambda_in, lambda_out, mass,
          f1, f2);
      const std::complex<double> scalar = gra::BilinearProduct(vbar, v);
      std::complex<double> ward = 0.0;

      for (const auto mu : dirac.LI) {
        const std::complex<double> gamma =
            dirac.gamma_lo[mu].BilinearForm(vbar, v);
        const std::complex<double> expected =
            (f1 + f2) * gamma +
            f2 * (momentum_sum % mu) / (2.0 * mass) * scalar;
        REQUIRE(std::abs(current[mu] - expected) == Approx(0.0).margin(2e-12));
        ward += vertex_q[mu] * current[mu];
      }
      CAPTURE(lambda_in, lambda_out, ward);
      REQUIRE(std::abs(ward) == Approx(0.0).margin(2e-12));
    }
  }
}

// Require conserved elastic currents to be independent of the gamma representation
TEST_CASE("MDirac elastic currents preserve the gamma basis and Ward identity",
          "[MDirac][photon][basis]") {
  using gra::aux::indices;
  const MDirac dirac("DIRAC");
  const MDirac chiral("CHIRAL");
  const double mass = gra::PDG::mp;
  const M4Vec incoming(0.31, -0.27, 3.0, std::sqrt(mass * mass + 9.169));
  const M4Vec outgoing(-0.25, 0.15, 2.6, std::sqrt(mass * mass + 6.845));
  const M4Vec q = outgoing - incoming;
  const std::array<std::array<double, 2>, 3> form_factors = {{{1.0, 0.0}, {0.0, 1.0}, {0.83, 1.47}}};

  for (const auto fermion : {MDirac::FermionKind::Particle, MDirac::FermionKind::Antiparticle}) {
    for (const auto &ff : form_factors) {
      const auto reference = dirac.DiracPauliCurrentBasis(incoming, outgoing, q, fermion, mass, ff[0], ff[1]);
      for (const MDirac *algebra : {&dirac, &chiral}) {
        const auto basis = algebra->DiracPauliCurrentBasis(incoming, outgoing, q, fermion, mass, ff[0], ff[1]);
        for (const auto &row : indices(basis)) {
          const double lambda_in = 0.5 * gra::spin::BinaryHelicityLabelX2(row % 2);
          const double lambda_out = 0.5 * gra::spin::BinaryHelicityLabelX2(row / 2);
          const auto current = algebra->DiracPauliCurrent(incoming, outgoing, q, fermion, lambda_in, lambda_out,
                                                         mass, ff[0], ff[1]);
          std::complex<double> ward = 0.0;
          std::complex<double> basis_ward = 0.0;
          CAPTURE(fermion, ff, row, algebra == &chiral);
          for (const auto &mu : indices(current)) {
            CHECK(std::abs(current[mu] - reference[row][mu]) < 2.0e-12);
            CHECK(std::abs(basis[row][mu] - reference[row][mu]) < 2.0e-12);
            ward += q[mu] * current[mu];
            basis_ward += q[mu] * basis[row][mu];
          }
          CHECK(std::abs(ward) < 2.0e-12);
          CHECK(std::abs(basis_ward) < 2.0e-12);
        }
      }
    }
  }
}

// Dirac algebra tests
TEST_CASE("MDirac multiple tests", "[MDirac]") {

  const double EPS = 1E-5;

  // First create some random particles
  MRandom rng;
  rng.SetSeed(123456);

  const double m0 = 100;
  const double mdaughter = 0.139;
  const int N = 2;

  M4Vec mother(0, 0, 0, m0);
  std::vector<double> m(N, mdaughter);

  // Create random 2-body decays
  for (std::size_t trial = 0; trial < 10; ++trial) {
    std::vector<M4Vec> p(N);
    gra::kinematics::TwoBodyPhaseSpace(mother, m0, m, p, rng);

    // Different basis
    std::vector<std::string> BASIS = {"DIRAC", "CHIRAL"};
    for (std::size_t k = 0; k < BASIS.size(); ++k) {

      MDirac dirac(BASIS[k]);

      // For each particle, do tests
      for (std::size_t i = 0; i < N; ++i) {

        std::vector<std::string> MODE = {"helicity", "gauge"};
        if (BASIS[k] == "DIRAC") {
          MODE.push_back("spin");
        }
        for (const auto &mode : MODE) {
          RequireSpinorCompleteness(dirac, p[i], false,
                                    BASIS[k] == "CHIRAL", mode, EPS);
          RequireSpinorCompleteness(dirac, p[i], true,
                                    BASIS[k] == "CHIRAL", mode, EPS);
        }

        RequireMassiveSpin1Completeness(dirac, p[i], EPS);
        REQUIRE((dirac.FSlash(p[i]) * dirac.FSlash(p[i]) -
                 dirac.I4 * p[i].M2())
                    .FrobNorm() < EPS);
      }
    }
  }
}

TEST_CASE("MDirac basis initialization preserves the Clifford algebra",
          "[MDirac][basis]") {
  const gra::M4Vec momentum(0.7, -1.1, 2.3, 3.4);

  for (const std::string basis : {"DIRAC", "CHIRAL"}) {
    CAPTURE(basis);
    const gra::MDirac dirac(basis);
    RequireCliffordAlgebra(dirac, 1e-12);
    const auto slash_square = dirac.FSlash(momentum) * dirac.FSlash(momentum);
    const auto expected = dirac.I4 * momentum.M2();
    for (std::size_t i = 0; i < 4; ++i) {
      for (std::size_t j = 0; j < 4; ++j) {
        REQUIRE(std::abs(slash_square(i, j) - expected(i, j)) < 2e-12);
      }
    }
  }
}

TEST_CASE("MDirac FSlash contracts a contravariant four-vector",
          "[MDirac][Clifford]") {
  gra::MDirac dirac;
  const gra::M4Vec vector(1.0, -2.0, 3.0, 5.0);
  const auto slash = dirac.FSlash(vector);

  gra::MMatrix<std::complex<double>> explicit_contraction(4, 4, 0.0);
  for (const auto mu : dirac.LI) {
    explicit_contraction += dirac.gamma_up[mu] * (vector % mu);
  }

  REQUIRE(slash.size_row() == 4);
  REQUIRE(slash.size_col() == 4);
  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      REQUIRE(std::abs(slash(i, j) - explicit_contraction(i, j)) ==
              Approx(0.0).margin(1e-15));
    }
  }
}

TEST_CASE("MDirac helicity spinors satisfy particle and antiparticle equations",
          "[MDirac][spinor][helicity]") {
  gra::MDirac dirac("DIRAC");
  gra::M4Vec momentum;
  momentum.SetPxPyPzM(1.0, -0.7, 2.1, 0.63);

  for (const int helicity : {-1, 1}) {
    const auto particle = dirac.uHelDirac(momentum, helicity);
    const auto antiparticle = dirac.vHelDirac(momentum, helicity);
    REQUIRE(particle.size() == 4);
    REQUIRE(antiparticle.size() == 4);
    RequireDiracSpinorForTest(dirac, momentum, particle, false);
    RequireDiracSpinorForTest(dirac, momentum, antiparticle, true);
  }
}

TEST_CASE("MDirac helicity spinors retain lower components near rest",
          "[MDirac][spinor][helicity][small-momentum]") {
  constexpr double mass = 1.0;
  constexpr double px = 1.0e-9;
  const gra::M4Vec momentum(px, 0.0, 0.0,
                            std::sqrt(mass * mass + px * px));
  const gra::MDirac dirac("DIRAC");

  for (const int helicity : {-1, 1}) {
    const auto particle = dirac.uHelDirac(momentum, helicity);
    const auto antiparticle = dirac.vHelDirac(momentum, helicity);
    const double particle_lower =
        std::sqrt(std::norm(particle[2]) + std::norm(particle[3]));
    const double antiparticle_upper =
        std::sqrt(std::norm(antiparticle[0]) + std::norm(antiparticle[1]));
    const double expected =
        std::sqrt(momentum.E() + mass) * px / (momentum.E() + mass);
    REQUIRE(particle_lower == Approx(expected).epsilon(2.0e-14));
    REQUIRE(antiparticle_upper == Approx(expected).epsilon(2.0e-14));
    REQUIRE(particle_lower > 0.0);
    REQUIRE(antiparticle_upper > 0.0);
    RequireDiracSpinorForTest(dirac, momentum, particle, false);
    RequireDiracSpinorForTest(dirac, momentum, antiparticle, true);
  }
}

// Check the distinct public JW and covariant spinor charts under z rotations
TEST_CASE("MDirac helicity spinor charts have explicit azimuth phases",
          "[MDirac][spinor][helicity][phase]") {
  const double delta = 0.67;
  gra::M4Vec momentum;
  momentum.SetPxPyPzM(0.8, -0.5, 1.7, 0.63);
  gra::M4Vec rotated = momentum;
  rotated.RotateZ(delta);

  const auto rotate_spinor =
      [delta](const gra::MDirac::Spinor &spinor,
              const std::complex<double> &section_phase) {
        gra::MDirac::Spinor out = spinor;
        const std::complex<double> up = std::exp(-0.5 * gra::math::zi * delta);
        const std::complex<double> down = std::exp(0.5 * gra::math::zi * delta);
        out[0] *= up;
        out[1] *= down;
        out[2] *= up;
        out[3] *= down;
        for (auto &component : out) {
          component *= section_phase;
        }
        return out;
      };
  const auto require_spinor = [](const gra::MDirac::Spinor &actual,
                                 const gra::MDirac::Spinor &expected) {
    for (std::size_t i = 0; i < actual.size(); ++i) {
      REQUIRE(std::abs(actual[i] - expected[i]) < 2.0e-12);
    }
  };

  for (const std::string basis : {"DIRAC", "CHIRAL"}) {
    gra::MDirac dirac(basis);
    for (const int helicity : {-1, 1}) {
      CAPTURE(basis, helicity);
      const std::complex<double> covariant_chart =
          std::exp(0.5 * gra::math::zi * delta);
      if (basis == "DIRAC") {
        require_spinor(dirac.uHelDirac(rotated, helicity),
                       rotate_spinor(dirac.uHelDirac(momentum, helicity),
                                     covariant_chart));
        require_spinor(dirac.vHelDirac(rotated, helicity),
                       rotate_spinor(dirac.vHelDirac(momentum, helicity),
                                     covariant_chart));
      } else {
        require_spinor(dirac.uHelChiral(rotated, helicity),
                       rotate_spinor(dirac.uHelChiral(momentum, helicity),
                                     covariant_chart));
        require_spinor(dirac.vHelChiral(rotated, helicity),
                       rotate_spinor(dirac.vHelChiral(momentum, helicity),
                                     covariant_chart));
      }
      require_spinor(
          dirac.uGauge(rotated, helicity),
          rotate_spinor(dirac.uGauge(momentum, helicity), covariant_chart));
      require_spinor(
          dirac.vGauge(rotated, helicity),
          rotate_spinor(dirac.vGauge(momentum, helicity), covariant_chart));
    }
  }

  gra::MDirac dirac("DIRAC");
  for (const int helicity : {-1, 1}) {
    const auto xi = dirac.XiSpinor(momentum, helicity);
    auto expected = xi;
    expected[0] *= std::exp(-0.5 * gra::math::zi * delta);
    expected[1] *= std::exp(0.5 * gra::math::zi * delta);
    const std::complex<double> jw_chart =
        std::exp(0.5 * gra::math::zi * static_cast<double>(helicity) * delta);
    for (auto &component : expected) {
      component *= jw_chart;
    }
    const auto actual = dirac.XiSpinor(rotated, helicity);
    REQUIRE(std::abs(actual[0] - expected[0]) < 2.0e-12);
    REQUIRE(std::abs(actual[1] - expected[1]) < 2.0e-12);
  }
}

// Check that polarization representatives rotate as Lorentz tensors without a
// JW phase
TEST_CASE("MDirac polarization charts are covariant under z rotations",
          "[MDirac][polarization][phase]") {
  const double delta = -0.73;
  gra::M4Vec massless(0.7, -0.4, 1.2,
                      std::sqrt(0.7 * 0.7 + 0.4 * 0.4 + 1.2 * 1.2));
  gra::M4Vec massive;
  massive.SetPxPyPzM(0.7, -0.4, 1.2, 0.9);
  gra::M4Vec massless_rotated = massless;
  gra::M4Vec massive_rotated = massive;
  massless_rotated.RotateZ(delta);
  massive_rotated.RotateZ(delta);
  const double c = std::cos(delta);
  const double s = std::sin(delta);
  const auto rotate_vector = [c, s](const auto &input) {
    auto out = input;
    out(1) = c * input(1) - s * input(2);
    out(2) = s * input(1) + c * input(2);
    return out;
  };
  gra::MDirac dirac("DIRAC");

  for (const int helicity : {-1, 1}) {
    const auto actual = dirac.EpsSpin1(massless_rotated, helicity);
    const auto expected = rotate_vector(dirac.EpsSpin1(massless, helicity));
    for (std::size_t mu = 0; mu < 4; ++mu) {
      REQUIRE(std::abs(actual(mu) - expected(mu)) < 2.0e-12);
    }
  }
  for (const int helicity : {-1, 0, 1}) {
    const auto actual = dirac.EpsMassiveSpin1(massive_rotated, helicity);
    const auto expected =
        rotate_vector(dirac.EpsMassiveSpin1(massive, helicity));
    for (std::size_t mu = 0; mu < 4; ++mu) {
      REQUIRE(std::abs(actual(mu) - expected(mu)) < 2.0e-12);
    }
  }

  const auto spin2 = dirac.EpsMassiveSpin2(massive, 1);
  const auto spin2_rotated = dirac.EpsMassiveSpin2(massive_rotated, 1);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      std::complex<double> expected = 0.0;
      for (std::size_t a = 0; a < 4; ++a) {
        const double r_mu_a = mu == 0 || mu == 3 ? (mu == a ? 1.0 : 0.0)
                              : mu == 1          ? (a == 1   ? c
                                                    : a == 2 ? -s
                                                             : 0.0)
                                                 : (a == 1   ? s
                                                    : a == 2 ? c
                                                             : 0.0);
        for (std::size_t b = 0; b < 4; ++b) {
          const double r_nu_b = nu == 0 || nu == 3 ? (nu == b ? 1.0 : 0.0)
                                : nu == 1          ? (b == 1   ? c
                                                      : b == 2 ? -s
                                                               : 0.0)
                                                   : (b == 1   ? s
                                                      : b == 2 ? c
                                                               : 0.0);
          expected += r_mu_a * r_nu_b * spin2(a, b);
        }
      }
      REQUIRE(std::abs(spin2_rotated(mu, nu) - expected) < 3.0e-12);
    }
  }
}

// Match the exact covariant FFV current to the Jacob-Wick decay matrix
TEST_CASE("MDirac FFV current matches the Jacob-Wick decay section",
          "[MDirac][current][gra::spin][phase]") {
  using Complex = std::complex<double>;
  gra::MDirac dirac("DIRAC");
  constexpr double mother_mass = 3.7;
  constexpr double fermion_mass = 0.61;
  const double momentum =
      std::sqrt(0.25 * mother_mass * mother_mass - fermion_mass * fermion_mass);
  const double energy = 0.5 * mother_mass;
  const gra::M4Vec parent(0.0, 0.0, 0.0, mother_mass);
  const Complex cL(0.73, -0.21);
  const Complex cR(-0.34, 0.18);
  constexpr std::array<int, 2> helicities = {-1, 1};

  const auto momenta = [momentum, energy](double theta, double phi) {
    const double st = std::sin(theta);
    const gra::M4Vec fermion(momentum * st * std::cos(phi),
                             momentum * st * std::sin(phi),
                             momentum * std::cos(theta), energy);
    return std::array<gra::M4Vec, 2>{
        fermion,
        gra::M4Vec(-fermion.Px(), -fermion.Py(), -fermion.Pz(), energy)};
  };
  const auto exact = [&](const std::array<gra::M4Vec, 2> &p, int hf, int ha,
                         int parent_spin) {
    const auto current =
        gra::qed::FFVCurrent(dirac, p[0], p[1], hf, ha, cL, cR);
    const auto eps = dirac.EpsMassiveSpin1(parent, parent_spin);
    Complex amplitude = 0.0;
    for (std::size_t mu = 0; mu < 4; ++mu) {
      amplitude += current[mu] * eps(mu);
    }
    return amplitude;
  };

  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, 1.0, 0.5, 0.5, {-1.0, 0.0, 1.0},
                              {-0.5, 0.5}, {-0.5, 0.5},
                              "Dirac Jacob-Wick current");
  hel.T = MMatrix<Complex>(2, 2, 0.0);
  const auto reference = momenta(0.0, 0.0);
  for (std::size_t i = 0; i < helicities.size(); ++i) {
    for (std::size_t j = 0; j < helicities.size(); ++j) {
      const int lambda = (helicities[i] - helicities[j]) / 2;
      // Continue the second-particle state through the south-pole JW section
      const double second_leg_phase =
          gra::spin::JacobWickSecondLegReversalPhase(
              0.5, 0.5 * static_cast<double>(helicities[j]));
      hel.T[i][j] = second_leg_phase *
                    exact(reference, helicities[i], helicities[j], lambda);
    }
  }

  for (const auto &[theta, phi] : std::array<std::pair<double, double>, 3>{
           std::pair{0.63, 0.0}, std::pair{0.63, 0.79},
           std::pair{2.14, -1.07}}) {
    CAPTURE(theta, phi);
    const auto p = momenta(theta, phi);
    const auto jw = gra::spin::fDecayMatrix(hel, theta, phi);
    for (std::size_t row = 0; row < hel.lambda_values.size_row(); ++row) {
      const int hf = helicities[row / 2];
      const int ha = helicities[row % 2];
      for (std::size_t col = 0; col < hel.Jz_values.size(); ++col) {
        const int parent_spin = static_cast<int>(hel.Jz_values[col]);
        const Complex covariant = exact(p, hf, ha, parent_spin);
        CAPTURE(row, col, covariant, jw[row][col]);
        REQUIRE(std::abs(covariant - jw[row][col]) <
                3.0e-11 * std::max(1.0, std::abs(covariant)));
      }
    }
  }
}

TEST_CASE("MDirac adjoint is conjugation followed by gamma zero",
          "[MDirac][adjoint]") {
  gra::MDirac dirac("DIRAC");
  const gra::MDirac::Spinor spinor = {
      std::complex<double>(1.0, 0.5), std::complex<double>(-0.2, 1.0),
      std::complex<double>(0.7, -0.4), std::complex<double>(-1.0, -0.3)};
  const gra::MDirac::Spinor expected = {
      std::conj(spinor[0]), std::conj(spinor[1]), -std::conj(spinor[2]),
      -std::conj(spinor[3])};

  const auto adjoint = dirac.Bar(spinor);
  REQUIRE(adjoint == expected);
  REQUIRE(dirac.Bar(adjoint) == spinor);
}

TEST_CASE("MDirac forward vector current has exact covariant normalization",
          "[MDirac][current]") {
  gra::MDirac dirac("DIRAC");
  gra::M4Vec momentum;
  momentum.SetPxPyPzM(0.7, -0.4, 2.3, 0.51);

  for (const int initial_helicity : {-1, 1}) {
    const auto particle = dirac.uHelDirac(momentum, initial_helicity);
    for (const int final_helicity : {-1, 1}) {
      const auto adjoint = dirac.Bar(dirac.uHelDirac(momentum, final_helicity));
      for (const auto mu : dirac.LI) {
        const std::complex<double> current =
            dirac.gamma_up[mu].BilinearForm(adjoint, particle);
        const double expected =
            initial_helicity == final_helicity ? 2.0 * (momentum ^ mu) : 0.0;
        REQUIRE(current.real() == Approx(expected).margin(2e-12));
        REQUIRE(current.imag() == Approx(0.0).margin(2e-12));
      }
    }
  }
}

// Test Spinor completeness relation
TEST_CASE("MDirac: Spinor completeness test") {
  gra::MDirac dirac;
  double m = 0.110;
  double p = 5.0;
  gra::M4Vec vec(p, 0.0, 0.0, std::sqrt(m * m + p * p));

  RequireSpinorCompleteness(dirac, vec, false, false, "helicity", 1e-8);
}

// Test polarization vectors for massive spin-1 particles
TEST_CASE("MDirac: Massive Spin-1 polarization test") {
  gra::MDirac dirac;
  double m = 0.110;
  double p = 5.0;
  gra::M4Vec vec(p, 0.0, 0.0, std::sqrt(m * m + p * p));

  RequireMassiveSpin1Completeness(dirac, vec, 1e-8);
}

TEST_CASE("MDirac: Basis validation, projectors and angular momentum operators",
          "[MDirac]") {
  using Complex = std::complex<double>;

  const gra::MDirac dirac;

  REQUIRE_THROWS_AS(gra::MDirac("BAD"), std::invalid_argument);

  const gra::MMatrix<Complex> PR = dirac.PR();
  const gra::MMatrix<Complex> PL = dirac.PL();
  const gra::MMatrix<Complex> I4 = dirac.I4;

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {
      const Complex expected_identity =
          (i == j) ? Complex(1.0, 0.0) : Complex(0.0, 0.0);
      REQUIRE(std::abs((PR + PL)(i, j) - expected_identity) ==
              Approx(0.0).margin(1e-12));
      REQUIRE(std::abs((PR * PR)(i, j) - PR(i, j)) ==
              Approx(0.0).margin(1e-12));
      REQUIRE(std::abs((PL * PL)(i, j) - PL(i, j)) ==
              Approx(0.0).margin(1e-12));
      REQUIRE(std::abs((PR * PL)(i, j)) == Approx(0.0).margin(1e-12));
      REQUIRE(std::abs(I4(i, j) - expected_identity) ==
              Approx(0.0).margin(1e-12));
    }
  }

  const gra::MMatrix<Complex> Jx = dirac.J_operator(1);
  const gra::MMatrix<Complex> Jy = dirac.J_operator(2);
  const gra::MMatrix<Complex> Jz = dirac.J_operator(3);

  REQUIRE(Jx(0, 1) == Complex(0.5, 0.0));
  REQUIRE(Jx(1, 0) == Complex(0.5, 0.0));
  REQUIRE(Jy(0, 1) == Complex(0.0, -0.5));
  REQUIRE(Jy(1, 0) == Complex(0.0, 0.5));
  REQUIRE(Jz(0, 0) == Complex(0.5, 0.0));
  REQUIRE(Jz(1, 1) == Complex(-0.5, 0.0));
  REQUIRE_THROWS_AS(dirac.J_operator(0), std::invalid_argument);
  REQUIRE_THROWS_AS(dirac.J_operator(4), std::invalid_argument);
}

TEST_CASE("MDirac: Charge conjugation and Xi spinor helpers", "[MDirac]") {
  using Complex = std::complex<double>;

  SECTION("Charge conjugation matrix is antisymmetric") {
    for (const auto &basis : {std::string("DIRAC"), std::string("CHIRAL")}) {
      gra::MDirac dirac(basis);
      const auto C = dirac.C_up();
      const auto CT = C.Transpose();
      for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
          REQUIRE(std::abs((CT + C)(i, j)) == Approx(0.0).margin(1e-12));
        }
      }
    }
  }

  SECTION("Xi spinors on the +z axis have the expected form") {
    gra::MDirac dirac;
    gra::M4Vec p(0.0, 0.0, 5.0, 5.0);

    const auto plus = dirac.XiSpinor(p, 1);
    const auto minus = dirac.XiSpinor(p, -1);

    REQUIRE(plus.size() == 2);
    REQUIRE(minus.size() == 2);
    REQUIRE(plus[0] == Complex(1.0, 0.0));
    REQUIRE(std::abs(plus[1]) == Approx(0.0).margin(1e-12));
    REQUIRE(std::abs(minus[0]) == Approx(0.0).margin(1e-12));
    REQUIRE(minus[1] == Complex(1.0, 0.0));
  }
}

TEST_CASE("MDirac: Propagators and polarization helpers", "[MDirac]") {
  using Complex = std::complex<double>;

  gra::MDirac dirac;
  const gra::M4Vec photon(0.0, 0.0, 3.0, 3.0);
  const gra::M4Vec massive(1.0, 2.0, 3.0, std::sqrt(30.0));

  SECTION("Photon propagator matches -i g_mu_nu / q2") {
    const double q2 = 7.5;
    const auto prop = dirac.iD_y(q2);
    for (const auto &mu : dirac.LI) {
      for (const auto &nu : dirac.LI) {
        REQUIRE(prop(mu, nu) == Complex(0.0, -dirac.g(mu, nu) / q2));
      }
    }
  }

  SECTION("Fermion propagator matches explicit formula") {
    const double m = 0.5;
    const auto lhs = dirac.iD_F(massive, m);
    const auto rhs = (dirac.FSlash(massive) + dirac.I4 * m) *
                     (gra::math::zi / (massive.M2() - m * m));
    for (std::size_t i = 0; i < 4; ++i) {
      for (std::size_t j = 0; j < 4; ++j) {
        REQUIRE(std::abs(lhs(i, j) - rhs(i, j)) == Approx(0.0).margin(1e-12));
      }
    }
  }

  SECTION("Massless polarization vectors are transverse and state helpers "
          "lower indices correctly") {
    for (const int helicity : {-1, 1}) {
      const auto eps = dirac.EpsSpin1(photon, helicity);
      Complex dot = 0.0;
      for (const auto &mu : dirac.LI) {
        dot += (photon % mu) * eps(mu);
      }
      REQUIRE(std::abs(dot) == Approx(0.0).margin(1e-12));
    }

    const auto up = dirac.MasslessSpin1States(photon, "none", true);
    const auto lo = dirac.MasslessSpin1States(photon, "none", false);
    const auto conj = dirac.MasslessSpin1States(photon, "conj", true);

    for (std::size_t s = 0; s < 2; ++s) {
      REQUIRE(lo[s](0) == up[s](0));
      for (std::size_t mu = 1; mu < 4; ++mu) {
        REQUIRE(lo[s](mu) == -up[s](mu));
        REQUIRE(conj[s](mu) == std::conj(up[s](mu)));
      }
    }

    REQUIRE_THROWS_AS(dirac.MasslessSpin1States(photon, "bad", true),
                      std::invalid_argument);
  }

  SECTION("Massive polarization state helpers and spin-2 tensor behave "
          "consistently") {
    const auto up = dirac.MassiveSpin1States(massive, "none", true);
    const auto lo = dirac.MassiveSpin1States(massive, "none", false);
    const auto conj = dirac.MassiveSpin1States(massive, "conj", true);

    for (std::size_t s = 0; s < 3; ++s) {
      REQUIRE(lo[s](0) == up[s](0));
      for (std::size_t mu = 1; mu < 4; ++mu) {
        REQUIRE(lo[s](mu) == -up[s](mu));
        REQUIRE(conj[s](mu) == std::conj(up[s](mu)));
      }
    }

    Complex completeness[4][4] = {};
    for (const int helicity : {-1, 0, 1}) {
      const auto eps = dirac.EpsMassiveSpin1(massive, helicity);
      Complex k_dot_eps = 0.0;
      Complex eps_dot_eps = 0.0;

      for (const auto &mu : dirac.LI) {
        const Complex eps_lower = (mu == 0) ? eps(mu) : -eps(mu);
        k_dot_eps += (massive % mu) * eps(mu);
        eps_dot_eps += std::conj(eps(mu)) * eps_lower;

        for (const auto &nu : dirac.LI) {
          completeness[mu][nu] += eps(mu) * std::conj(eps(nu));
        }
      }

      REQUIRE(std::abs(k_dot_eps) == Approx(0.0).margin(1e-12));
      REQUIRE(std::abs(eps_dot_eps + Complex(1.0, 0.0)) ==
              Approx(0.0).margin(1e-12));
    }

    const double mass2 = massive.M2();
    for (const auto &mu : dirac.LI) {
      for (const auto &nu : dirac.LI) {
        const double expected =
            -dirac.g(mu, nu) + (massive ^ mu) * (massive ^ nu) / mass2;
        REQUIRE(std::abs(completeness[mu][nu] - expected) ==
                Approx(0.0).margin(1e-12));
      }
    }

    Complex spin2_sum[4][4][4][4] = {};
    for (const int helicity : {-2, -1, 0, 1, 2}) {
      const auto eps2 = dirac.EpsMassiveSpin2(massive, helicity);
      Complex trace = 0.0;
      for (const auto &mu : dirac.LI) {
        Complex transverse = 0.0;
        for (const auto &nu : dirac.LI) {
          REQUIRE(std::abs(eps2(mu, nu) - eps2(nu, mu)) ==
                  Approx(0.0).margin(1e-12));
          trace += dirac.g(mu, nu) * eps2(mu, nu);
          transverse += (massive % nu) * eps2(nu, mu);
          for (const auto &rho : dirac.LI) {
            for (const auto &sigma : dirac.LI) {
              spin2_sum[mu][nu][rho][sigma] +=
                  eps2(mu, nu) * std::conj(eps2(rho, sigma));
            }
          }
        }
        REQUIRE(std::abs(transverse) == Approx(0.0).margin(2e-12));
      }
      REQUIRE(std::abs(trace) == Approx(0.0).margin(2e-12));
    }

    for (const auto &mu : dirac.LI) {
      for (const auto &nu : dirac.LI) {
        for (const auto &rho : dirac.LI) {
          for (const auto &sigma : dirac.LI) {
            const double p_mu_rho =
                -dirac.g(mu, rho) +
                (massive ^ mu) * (massive ^ rho) / mass2;
            const double p_nu_sigma =
                -dirac.g(nu, sigma) +
                (massive ^ nu) * (massive ^ sigma) / mass2;
            const double p_mu_sigma =
                -dirac.g(mu, sigma) +
                (massive ^ mu) * (massive ^ sigma) / mass2;
            const double p_nu_rho =
                -dirac.g(nu, rho) +
                (massive ^ nu) * (massive ^ rho) / mass2;
            const double p_mu_nu =
                -dirac.g(mu, nu) +
                (massive ^ mu) * (massive ^ nu) / mass2;
            const double p_rho_sigma =
                -dirac.g(rho, sigma) +
                (massive ^ rho) * (massive ^ sigma) / mass2;
            const double expected =
                0.5 * (p_mu_rho * p_nu_sigma + p_mu_sigma * p_nu_rho) -
                p_mu_nu * p_rho_sigma / 3.0;
            REQUIRE(std::abs(spin2_sum[mu][nu][rho][sigma] - expected) ==
                    Approx(0.0).margin(4e-12));
          }
        }
      }
    }

    REQUIRE_THROWS_AS(dirac.MassiveSpin1States(massive, "bad", true),
                      std::invalid_argument);
  }
}

TEST_CASE("MDirac: Spinor and gauge-state helper coverage", "[MDirac]") {
  SECTION("Dirac basis helper methods return valid four-component spinors") {
    gra::MDirac dirac("DIRAC");
    const gra::M4Vec p(0.3, 0.4, 1.2,
                       std::sqrt(1.0 + 0.3 * 0.3 + 0.4 * 0.4 + 1.2 * 1.2));

    for (const int helicity : {-1, 1}) {
      RequireDiracSpinorForTest(dirac, p, dirac.uHelDirac(p, helicity), false);
      RequireDiracSpinorForTest(dirac, p, dirac.vHelDirac(p, helicity), true);
      RequireDiracSpinorForTest(dirac, p, dirac.uDirac(p, helicity), false);
      RequireDiracSpinorForTest(dirac, p, dirac.vDirac(p, helicity), true);
      RequireDiracSpinorForTest(dirac, p, dirac.uGauge(p, helicity), false);
      RequireDiracSpinorForTest(dirac, p, dirac.vGauge(p, helicity), true);
    }

    const auto u = dirac.SpinorStates(p, "u");
    const auto ubar = dirac.SpinorStates(p, "ubar");
    const auto v = dirac.SpinorStates(p, "v");
    const auto vbar = dirac.SpinorStates(p, "vbar");
    for (std::size_t state = 0; state < u.size(); ++state) {
      REQUIRE(ubar[state] == dirac.Bar(u[state]));
      REQUIRE(vbar[state] == dirac.Bar(v[state]));
    }

    REQUIRE(std::abs(dirac.sProd(p, p, 1) - std::conj(dirac.sProd(p, p, -1))) ==
            Approx(0.0).margin(1e-12));
  }

  SECTION("Chiral basis helper methods return valid four-component spinors and "
          "reject wrong API usage") {
    gra::MDirac dirac("CHIRAL");
    const gra::M4Vec p(
        0.2, -0.3, 1.0,
        std::sqrt(0.5 * 0.5 + 0.2 * 0.2 + 0.3 * 0.3 + 1.0 * 1.0));

    for (const int helicity : {-1, 1}) {
      RequireDiracSpinorForTest(dirac, p, dirac.uHelChiral(p, helicity), false);
      RequireDiracSpinorForTest(dirac, p, dirac.vHelChiral(p, helicity), true);
      RequireDiracSpinorForTest(dirac, p, dirac.uGauge(p, helicity), false);
      RequireDiracSpinorForTest(dirac, p, dirac.vGauge(p, helicity), true);
    }

    const auto u = dirac.SpinorStates(p, "u");
    const auto vbar = dirac.SpinorStates(p, "vbar");
    for (std::size_t state = 0; state < u.size(); ++state) {
      REQUIRE(u[state].size() == 4);
      REQUIRE(vbar[state].size() == 4);
    }

  }

  SECTION("Generated spacelike momenta use the amplitude failure type") {
    gra::MDirac dirac("DIRAC");
    const gra::M4Vec spacelike(0.0, 0.0, 1.0, 0.5);

    REQUIRE_THROWS_AS(dirac.uHelDirac(spacelike, 1), gra::AmplitudeFailure);
  }
}
