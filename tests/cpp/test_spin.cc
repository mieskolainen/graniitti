// Spin and decay tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "support/models_test_support.hh"

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Spin/MSpin.h"
#include "Graniitti/Spin/MSpinDensity.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

// Check entanglement of physical scalar vector-pair LS states and their coherent phases
TEST_CASE("Scalar vector decays give normalized bipartite spin entanglement", "[gra::spin][entanglement]") {
  const auto scalar = ToyParticle(9000901, 0, 1, 1, "scalar");
  const auto vector = ToyParticle(113, 2, -1, -1, "vector");
  const auto parent = MMatrix<std::complex<double>>::IdentityMatrix(1);
  gra::HELMatrix hel;
  hel.alpha_ls.Set(0, 0, 1.0);
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "entangled vector pair", false, false);
  const auto rho = gra::spin::DecayDensity(hel, parent, 0.0, 0.0);
  const auto metrics = gra::spin::EntanglementMetrics(rho, 3, 3);
  std::vector<std::complex<double>> singlet(9, 0.0);
  singlet[0] = singlet[8] = 1.0 / std::sqrt(3.0);
  singlet[4] = -1.0 / std::sqrt(3.0);
  RequireMatrixNear(rho, gra::RankOneProjector(singlet), 1e-13);
  CHECK(metrics.pair.purity == Approx(1.0).margin(1e-13));
  CHECK(metrics.pair.entropy == Approx(0.0).margin(1e-13));
  CHECK(metrics.entropy1 == Approx(std::log2(3.0)).margin(1e-13));
  CHECK(metrics.entropy2 == Approx(std::log2(3.0)).margin(1e-13));
  CHECK(metrics.negativity == Approx(1.0).margin(1e-13));
  CHECK(metrics.log_negativity == Approx(std::log2(3.0)).margin(1e-13));
  RequireMatrixNear(gra::spin::DecayDensity(hel, parent, 0.81, -0.43), rho, 1e-13);

  // Local spin rotations and exchange of the two momentum modes preserve entanglement
  const auto rotation = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.7, -0.4), 1.0)
                            .Kronecker(gra::spin::SpinRotation(gra::spin::SpinHalfRotation(-0.3, 0.9), 1.0));
  const auto rotated = gra::spin::EntanglementMetrics(rotation * rho * rotation.Dagger(), 3, 3);
  CHECK(rotated.negativity == Approx(metrics.negativity).margin(1e-13));
  CHECK(rotated.entropy1 == Approx(metrics.entropy1).margin(1e-13));
  const std::vector<std::size_t> exchange = {0, 3, 6, 1, 4, 7, 2, 5, 8};
  RequireMatrixNear(rho.SelectRows(exchange).SelectColumns(exchange), rho, 1e-13);

  // Coherent S and D cancellation selects a separable longitudinal pair
  hel.alpha_ls.Clear();
  hel.alpha_ls.Set(0, 0, 1.0);
  hel.alpha_ls.Set(2, 4, -std::sqrt(2.0));
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "longitudinal vector pair", false, false);
  const auto longitudinal = gra::spin::EntanglementMetrics(gra::spin::DecayDensity(hel, parent, 0.0, 0.0), 3, 3);
  CHECK(longitudinal.pair.purity == Approx(1.0).margin(1e-13));
  CHECK(longitudinal.entropy1 == Approx(0.0).margin(1e-13));
  CHECK(longitudinal.negativity == Approx(0.0).margin(1e-13));

  // A pure D wave has Schmidt probabilities 1/6, 2/3, 1/6
  hel.alpha_ls.Clear();
  hel.alpha_ls.Set(2, 4, 1.0);
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "D-wave vector pair", false, false);
  const auto d_wave = gra::spin::EntanglementMetrics(gra::spin::DecayDensity(hel, parent, 0.0, 0.0), 3, 3);
  CHECK(d_wave.entropy1 == Approx(std::log2(6.0) / 3.0 - 2.0 * std::log2(2.0 / 3.0) / 3.0).margin(1e-13));
  CHECK(d_wave.negativity == Approx(5.0 / 6.0).margin(1e-13));

  // Real photons retain two physical polarizations in the embedded spin-one basis
  const auto photon = ToyParticle(22, 2, -1, -1, "photon");
  hel.alpha_ls.Clear();
  hel.alpha_ls.Set(0, 0, 1.0);
  gra::spin::InitTMatrix(hel, scalar, photon, photon, false, "entangled photon pair", false, false);
  const auto photons = gra::spin::EntanglementMetrics(gra::spin::DecayDensity(hel, parent, 0.0, 0.0), 3, 3);
  CHECK(photons.entropy1 == Approx(1.0).margin(1e-13));
  CHECK(photons.negativity == Approx(0.5).margin(1e-13));
}

// Distinguish classical spin correlations from negativity in mixed states
TEST_CASE("Mixed pair spin metrics distinguish entanglement from entropy", "[gra::spin][entanglement]") {
  using Complex = std::complex<double>;
  std::vector<Complex> state(9, 0.0);
  state[0] = state[4] = state[8] = 1.0 / std::sqrt(3.0);
  const auto pure = gra::RankOneProjector(state);
  const auto classical = MMatrix<Complex>::DiagonalMatrix(pure.GetDiag());
  const auto mixed = gra::spin::EntanglementMetrics(classical, 3, 3);
  CHECK(mixed.pair.purity == Approx(1.0 / 3.0).margin(1e-13));
  CHECK(mixed.pair.entropy == Approx(std::log2(3.0)).margin(1e-13));
  CHECK(mixed.entropy1 == Approx(std::log2(3.0)).margin(1e-13));
  CHECK(mixed.negativity == Approx(0.0).margin(1e-13));
  for (const double p : {0.0, 0.2, 0.25, 0.4, 1.0}) {
    const auto rho = pure * p + MMatrix<Complex>::IdentityMatrix(9) * ((1.0 - p) / 9.0);
    const auto metrics = gra::spin::EntanglementMetrics(rho, 3, 3);
    CHECK(metrics.pair.purity == Approx(p * p + (1.0 - p * p) / 9.0).margin(1e-13));
    CHECK(metrics.negativity == Approx(std::max(0.0, (4.0 * p - 1.0) / 3.0)).margin(1e-13));
  }
  CHECK(gra::spin::DensityMetrics(classical * 2.0).purity == Approx(mixed.pair.purity).margin(1e-13));
  CHECK(gra::spin::EntanglementMetrics(pure * 7.0, 3, 3).negativity == Approx(1.0).margin(1e-13));
  REQUIRE_THROWS(gra::spin::DensityMetrics(classical * 0.0));
  REQUIRE_THROWS(gra::spin::DensityMetrics(MMatrix<Complex>{{1.1, 0.0}, {0.0, -0.1}}));
  REQUIRE_THROWS(gra::spin::EntanglementMetrics(pure, 2, 3));
}

// Trace the parent spin before assigning entanglement to a nonzero-spin decay
TEST_CASE("Tensor decay pair density retains the parent spin mixture", "[gra::spin][entanglement]") {
  using Complex = std::complex<double>;
  const auto tensor = ToyParticle(9000902, 4, 1, 1, "tensor");
  const auto vector = ToyParticle(113, 2, -1, -1, "vector");
  gra::HELMatrix hel;
  hel.alpha_ls.Set(0, 4, 1.0);
  gra::spin::InitTMatrix(hel, tensor, vector, vector, false, "unpolarized tensor decay", false, false);
  const auto parent = MMatrix<Complex>::IdentityMatrix(5) / 5.0;
  const auto rho = gra::spin::DecayDensity(hel, parent, 0.0, 0.0);
  const auto mixed = gra::spin::EntanglementMetrics(rho, 3, 3);
  CHECK(mixed.pair.purity == Approx(0.2).margin(1e-13));
  CHECK(mixed.pair.entropy == Approx(std::log2(5.0)).margin(1e-13));
  CHECK(mixed.entropy1 == Approx(std::log2(3.0)).margin(1e-13));
  CHECK(mixed.negativity == Approx(0.0).margin(1e-13));
  RequireMatrixNear(gra::spin::DecayDensity(hel, parent, 0.91, -0.64), rho, 1e-13);

  MMatrix<Complex> aligned(5, 5, 0.0);
  aligned[2][2] = 1.0;
  const auto pure = gra::spin::EntanglementMetrics(gra::spin::DecayDensity(hel, aligned, 0.0, 0.0), 3, 3);
  CHECK(pure.pair.purity == Approx(1.0).margin(1e-13));
  CHECK(pure.negativity == Approx(5.0 / 6.0).margin(1e-13));
  hel.T *= 0.0;
  REQUIRE_THROWS(gra::spin::DecayDensity(hel, parent, 0.0, 0.0));
}

// Check physical transverse LS interference and continuity at the pole
TEST_CASE("Photon LS radial matrices retain their pole normalization and interference",
          "[gra::spin][radial][photon]") {
  const auto scalar = ToyParticle(9000901, 0, 1, 0, "scalar");
  const auto photon = ToyParticle(22, 2, -1, -1, "photon");
  for (const double alpha_d : {1.0, -0.5}) {
    CAPTURE(alpha_d);
    gra::HELMatrix hel;
    hel.alpha_ls.Set(0, 0, 1.0);
    hel.alpha_ls.Set(2, 4, alpha_d);
    gra::spin::InitTMatrix(hel, scalar, photon, photon, false, "radial photons", false, false);
    REQUIRE(hel.T.FrobNorm2() == Approx(1.0).epsilon(1e-12));
    MMatrix<std::complex<double>> cached(hel.T.size_row(), hel.T.size_col(), 0.0);
    for (const auto &component : hel.ls_components) { cached += component.matrix * component.alpha; }
    RequireMatrixNear(cached, hel.T, 1e-12);
    for (const double ratio : {0.0, 0.5, 1.0 - 1e-8, 1.0, 1.0 + 1e-8, std::sqrt(2.0), 2.0}) {
      CAPTURE(ratio);
      const auto amplitude = gra::spin::DecayLSHelicityMatrix(hel, ratio, true);
      const double expected = gra::math::pow2((1.0 + alpha_d * ratio * ratio) / (1.0 + alpha_d));
      CHECK(amplitude.FrobNorm2() == Approx(expected).margin(1e-12));
      CHECK(gra::spin::DecayLSIntensity(hel, ratio, true) == Approx(expected).margin(1e-12));
      CHECK(gra::spin::DecayLSIntensity(hel, ratio, false) == Approx(1.0));
    }
    REQUIRE_THROWS_AS(gra::spin::DecayLSHelicityMatrix(hel, -1.0, true), gra::AmplitudeFailure);
    REQUIRE_THROWS_AS(gra::spin::DecayLSIntensity(hel, std::numeric_limits<double>::infinity(), true),
                      gra::AmplitudeFailure);
    hel.T[0][0] = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_FALSE(gra::spin::DecayLSHelicityMatrix(hel, 1.0, true).IsFinite());
    REQUIRE_FALSE(gra::spin::DecayLSHelicityMatrix(hel, 0.5, false).IsFinite());
  }
}

// Check that a full massive spin sum preserves the orthogonality of LS partial waves
TEST_CASE("Massive LS radial intensities preserve the partial-wave norm", "[gra::spin][radial]") {
  const auto scalar = ToyParticle(9000920, 0, 1, 0, "X");
  auto vector = ToyParticle(9000921, 2, -1, 0, "V");
  vector.mass = 0.7;
  gra::HELMatrix hel;
  const std::complex<double> alpha_d(-0.3, 0.4);
  hel.alpha_ls.Set(0, 0, 1.0);
  hel.alpha_ls.Set(2, 4, alpha_d);
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "massive radial sum", false, false);
  for (const double ratio : {0.0, 0.5, 1.0, 1.7}) {
    const double expected = (1.0 + std::norm(alpha_d) * gra::math::pow4(ratio)) / (1.0 + std::norm(alpha_d));
    CHECK(gra::spin::DecayLSIntensity(hel, ratio, true) == Approx(expected).epsilon(1e-12));
    CHECK(gra::spin::DecayLSHelicityMatrix(hel, ratio, true).FrobNorm2() == Approx(expected).epsilon(1e-12));
  }
}

// Check finite high-spin rotations and Condon-Shortley singlet phases without factorial products
TEST_CASE("Wigner functions remain normalized for large supported spins", "[gra::spin][wigner]") {
  for (const double J : {20.0, 49.5, 50.0, 80.0}) {
    CAPTURE(J);
    CHECK(gra::wigner::d(0.7, J, J, J) == Approx(std::pow(std::cos(0.35), 2.0 * J)).epsilon(1e-10));
    CHECK(gra::wigner::W3j(J, J, 0.0, J, -J, 0.0) == Approx(1.0 / std::sqrt(2.0 * J + 1.0)).epsilon(1e-10));
    const auto projections = gra::spin::SpinRep::FromSpin(J, "Wigner test").Projections();
    const auto rotation = gra::wigner::dRows(0.7, projections, projections, J);
    const MMatrix<double> identity(projections.size(), projections.size(), "eye");
    CHECK((rotation * rotation.Transpose() - identity).FrobNorm() < 1e-10);
    CHECK((rotation * gra::wigner::dRows(-0.7, projections, projections, J) - identity).FrobNorm() < 1e-10);
  }
  CHECK(gra::wigner::d(0.7, 85.0, 85.0, 85.0) ==
        Approx(std::pow(std::cos(0.35), 170)).epsilon(1e-10));
  CHECK(gra::wigner::W3j(49.5, 49.5, 0.0, -49.5, 49.5, 0.0) == Approx(-0.1).epsilon(1e-10));
  // The stretched multiplet has positive CG coefficients even when its boundary values are tiny
  const double stretched = std::exp(gra::math::LogGamma(85.0) - 2.0 * gra::math::LogGamma(43.0) -
                                    0.5 * (gra::math::LogGamma(169.0) - 2.0 * gra::math::LogGamma(85.0)));
  CHECK(gra::wigner::CG(42.0, 42.0, 0.0, 0.0, 84.0, 0.0) == Approx(stretched).epsilon(1e-10));
}

// Check composition, spinor signs and inverse boosts in the common spin basis
TEST_CASE("Spin frame transport preserves SU2 phases and inverse boosts", "[gra::spin][cascade]") {
  const auto a = gra::spin::SpinHalfRotation(0.73, -0.41);
  const auto b = gra::spin::SpinHalfRotation(-0.26, 1.37);
  for (const double J : {0.0, 0.5, 1.0, 1.5, 2.0}) {
    CAPTURE(J);
    RequireMatrixNear(gra::spin::SpinRotation(a * b, J),
                       gra::spin::SpinRotation(a, J) * gra::spin::SpinRotation(b, J), 1e-12);
    const auto full_turn = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.0, 2.0 * gra::math::PI), J);
    const double phase = static_cast<int>(std::llround(2.0 * J)) % 2 == 0 ? 1.0 : -1.0;
    RequireMatrixNear(full_turn, MMatrix<std::complex<double>>(full_turn.size_row(), full_turn.size_col(), "eye") * phase, 1e-12);
  }
  const gra::M4Vec p(0.31, -0.27, 0.42, 1.3);
  const gra::M4Vec inverse(-p.Px(), -p.Py(), -p.Pz(), p.E());
  RequireMatrixNear(gra::spin::SpinHalfBoost(p) * gra::spin::SpinHalfBoost(inverse),
                     MMatrix<std::complex<double>>(2, 2, "eye"), 1e-12);
  for (const double scale : {1e-200, 1e200}) {
    RequireMatrixNear(gra::spin::SpinHalfBoost(p * scale), gra::spin::SpinHalfBoost(p), 1e-12);
  }
  const gra::M4Vec fast(0.3, 0.4, 1.0, std::sqrt(1.25 + 1e-8));
  const auto frame = gra::spin::SpinHalfBoost(gra::M4Vec(100.0, 20.0, 30.0, std::sqrt(11301.0)));
  const auto wigner = gra::spin::SpinHalfWigner(frame, fast);
  RequireMatrixNear(wigner.Dagger() * wigner, MMatrix<std::complex<double>>(2, 2, "eye"), 1e-12);
  const auto canonical_boost = frame * gra::spin::SpinHalfBoost(fast) * wigner.Dagger();
  CHECK((canonical_boost - canonical_boost.Dagger()).FrobNorm() < 1e-12 * canonical_boost.FrobNorm());
  CHECK(canonical_boost.SelfAdjointEigenvalues(1e-12).front() > 0.0);
  RequireMatrixNear(gra::spin::SpinHalfWigner(a * frame, inverse),
                     a * gra::spin::SpinHalfWigner(frame, inverse), 1e-12);
  REQUIRE_THROWS_AS(gra::spin::SpinHalfBoost(gra::M4Vec(0.0, 0.0, 1.0, 1.0)), gra::AmplitudeFailure);
}

// Check Cartesian frame transport against composed spin rotations including polar axes
TEST_CASE("Spin frames reproduce Cartesian rotations", "[gra::spin][frame]") {
  for (const double theta : {0.0, 0.73, gra::math::PI}) {
    for (const double phi : {0.0, -0.41, 2.63}) {
      for (const double gamma : {0.0, 1.17, -gra::math::PI}) {
        CAPTURE(theta, phi, gamma);
        gra::M4Vec x(1.0, 0.0, 0.0, 0.0), z(0.0, 0.0, 1.0, 0.0);
        for (auto *axis : {&x, &z}) {
          axis->RotateZ(gamma);
          axis->RotateY(theta);
          axis->RotateZ(phi);
        }
        const auto u = gra::spin::SpinHalfRotation(theta, phi) * gra::spin::SpinHalfRotation(0.0, gamma);
        for (const double J : {0.0, 1.0, 2.0, 3.0, 6.0}) {
          CAPTURE(J);
          RequireMatrixNear(gra::spin::SpinFrame(x, z, J), gra::spin::SpinRotation(u, J), 2.0e-12);
        }
      }
    }
  }
}

// Identical stable siblings are already antisymmetric through the LS selection rule
TEST_CASE("Scalar identical fermion decay survives coherent symmetrization",
          "[gra::spin][cascade][symmetry][fermion]") {
  gra::PARAM_RES res;
  res.p = ToyParticle(9000902, 0, 1, 0, "scalar");
  res.p.mass = 3.0;
  auto fermion = ToyParticle(9000903, 1, 1, 0, "fermion");
  fermion.mass = 1.0;
  res.hel_decay.alpha_ls.Set(0, 0, 1.0);
  gra::spin::InitTMatrix(res.hel_decay, res.p, fermion, fermion, false, "identical fermions", false, false);
  gra::LORENTZSCALAR lts;
  lts.process.SPINDEC = true;
  lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 3.0);
  gra::MDecayBranch a, b;
  a.p = b.p = fermion;
  a.p4 = gra::M4Vec(0.3, 0.4, 1.0, 1.5);
  b.p4 = gra::M4Vec(-0.3, -0.4, -1.0, 1.5);
  for (const double angle : {0.0, 0.7, 2.3}) {
    lts.decaytree = {a, b};
    for (auto &leaf : lts.decaytree) { leaf.p4.RotateY(angle); leaf.p4.RotateZ(0.43); }
    lts.amplitude.DECAY_SYM = false;
    const auto single = gra::spin::ResonanceDecayMatrix(lts, res, "CM");
    REQUIRE(single.FrobNorm2() == Approx(1.0).epsilon(1e-12));
    lts.amplitude.DECAY_SYM = true;
    RequireMatrixNear(gra::spin::ResonanceDecayMatrix(lts, res, "CM"), single, 1e-12);
    REQUIRE(lts.decay_symmetry_assignments.size() == 1);
  }
}

// Compare crossed massive and massless spin states after rotations and identical-particle exchange
TEST_CASE("Coherent spinful cascades preserve rotational covariance and exchange",
          "[gra::spin][cascade][symmetry][fermion][photon][rotation]") {
  const int spin_type = GENERATE(0, 1, 2);
  const int root_spin_x2 = GENERATE(0, 2);
  const bool photons = spin_type == 1;
  CAPTURE(spin_type, root_spin_x2);
  auto fermion = ToyParticle(photons ? 22 : 9000910, photons ? 2 : 1, photons ? -1 : 1, 0, "f");
  fermion.mass = spin_type == 0 ? 0.35 : 0.0;
  auto scalar = ToyParticle(9000911, 0, 1, 0, "s");
  scalar.mass = 0.2;
  gra::MDecayBranch f1, f2, s1, s2;
  f1.p = f2.p = fermion;
  s1.p = s2.p = scalar;
  s2.p.pdg = 9000912;
  f1.p4 = gra::M4Vec(0.31, -0.15, 0.42, std::sqrt(fermion.mass * fermion.mass + 0.31 * 0.31 + 0.15 * 0.15 + 0.42 * 0.42));
  f2.p4 = gra::M4Vec(-0.18, -0.23, -0.31, std::sqrt(fermion.mass * fermion.mass + 0.18 * 0.18 + 0.23 * 0.23 + 0.31 * 0.31));
  s1.p4 = gra::M4Vec(-0.21, 0.32, 0.13, std::sqrt(0.2 * 0.2 + 0.21 * 0.21 + 0.32 * 0.32 + 0.13 * 0.13));
  s2.p4 = gra::M4Vec(0.08, 0.06, -0.24, std::sqrt(0.2 * 0.2 + 0.08 * 0.08 + 0.06 * 0.06 + 0.24 * 0.24));
  gra::MDecayBranch left, right;
  left.p = ToyParticle(9000913, fermion.spinX2, fermion.P, 0, "F1");
  right.p = ToyParticle(9000914, fermion.spinX2, fermion.P, 0, "F2");
  left.legs = {f1, s1};
  right.legs = {f2, s2};
  for (auto *branch : {&left, &right}) {
    branch->p4 = branch->legs[0].p4 + branch->legs[1].p4;
    branch->p.mass = branch->p4.M();
    branch->p.width = 0.12;
    branch->hel.alpha_ls.Set(0, fermion.spinX2, 1.0);
    branch->hel.g_decay = 1.0;
    gra::spin::InitTMatrix(branch->hel, branch->p, branch->legs[0].p, branch->legs[1].p,
                           false, "fermion cascade", false, false);
  }
  gra::LORENTZSCALAR lts;
  lts.process.SPINDEC = true;
  lts.amplitude.DECAY_SYM = true;
  lts.decaytree = {left, right};
  lts.pfinal[0] = left.p4 + right.p4;
  gra::PARAM_RES res;
  res.p = ToyParticle(9000915, root_spin_x2, 1, 0, "X");
  res.p.mass = lts.pfinal[0].M();
  res.hel_decay.alpha_ls.Set(0, root_spin_x2, 1.0);
  gra::spin::InitTMatrix(res.hel_decay, res.p, left.p, right.p, false, "cascade root", false, false);
  const auto reference = gra::spin::ResonanceDecayMatrix(lts, res, "CM");
  const double norm = reference.FrobNorm2();
  const auto density = reference * reference.Dagger();
  REQUIRE(norm > 1e-10);
  REQUIRE(lts.decay_symmetry_assignments.size() == 2);
  auto reordered = lts;
  for (auto &branch : reordered.decaytree) {
    std::reverse(branch.hel.Jz_values.begin(), branch.hel.Jz_values.end());
    gra::spin::InitJWRotation(branch.hel);
  }
  RequireMatrixNear(gra::spin::ResonanceDecayMatrix(reordered, res, "CM"), reference, 1e-11);
  MMatrix<std::complex<double>> projected(reference.size_row(), reference.size_col(), 0.0);
  for (const double projection : lts.decaytree[0].hel.Jz_values) {
    auto selected = lts;
    selected.decaytree[0].hel.Jz_values = {projection};
    gra::spin::InitJWRotation(selected.decaytree[0].hel);
    projected += gra::spin::ResonanceDecayMatrix(selected, res, "CM");
  }
  RequireMatrixNear(projected, reference, 1e-11);
  for (const bool exchange : {false, true}) {
    for (const double angle : {0.0, 0.4, 1.2, -0.7}) {
      CAPTURE(exchange, angle);
      auto moved = lts;
      if (exchange) {
        REQUIRE(gra::spin::ApplyStableLeafSymmetryAssignment(moved.decaytree, lts.decay_symmetry_assignments[1]));
      }
      for (auto &branch : moved.decaytree) {
        for (auto &leaf : branch.legs) { leaf.p4.RotateY(angle); leaf.p4.RotateZ(0.3); }
        branch.p4 = branch.legs[0].p4 + branch.legs[1].p4;
      }
      moved.pfinal[0] = moved.decaytree[0].p4 + moved.decaytree[1].p4;
      const auto actual = gra::spin::ResonanceDecayMatrix(moved, res, "CM");
      CHECK(actual.FrobNorm2() ==
            Approx(norm).epsilon(1e-10));
      const auto rotation = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(angle, 0.3), 0.5 * root_spin_x2);
      const auto expected_density = rotation.Conj() * density * rotation.Transpose();
      const auto actual_density = actual * actual.Dagger();
      CHECK((actual_density - expected_density).FrobNorm() < 1e-10 * (1.0 + density.FrobNorm()));
    }
  }
}

// Close parity and Bose exchange for every representative of a four-state orbit
TEST_CASE("compact helicity input closes parity and identical vector exchange",
          "[gra::spin][helicity][parity][exchange]") {
  gra::MParticle mother;
  mother.pdg = 225;
  mother.spinX2 = 4;
  mother.P = mother.C = 1;
  mother.mass = 3.0;
  gra::MParticle vector;
  vector.pdg = 113;
  vector.spinX2 = 2;
  vector.P = vector.C = -1;
  vector.mass = 0.77;
  const std::complex<double> coupling = std::polar(0.8, 0.37);
  const std::vector<std::array<double, 2>> orbit = {{1, 0}, {-1, 0}, {0, 1}, {0, -1}};
  const auto complete = gra::BuildDirectHelicityCoupling(
      mother, {vector, vector}, orbit, std::vector<std::complex<double>>(orbit.size(), coupling),
      true, true, false, "complete vector orbit", false);
  for (const auto &helicity : orbit) {
    const auto compact = gra::BuildDirectHelicityCoupling(
        mother, {vector, vector}, {helicity}, {coupling}, true, true, false, "compact vector orbit", false);
    CHECK((compact.T - complete.T).FrobNorm2() == Approx(0.0).margin(1e-24));
    CHECK(compact.T.FrobNorm2() == Approx(4.0 * std::norm(coupling)).epsilon(1e-12));
    for (const auto &partner : orbit) {
      const auto row = gra::spin::SpinProjectionIndex(partner[0], 1.0, "vector 1");
      const auto col = gra::spin::SpinProjectionIndex(partner[1], 1.0, "vector 2");
      CHECK(std::abs(compact.T[row][col] - coupling) == Approx(0.0).margin(1e-12));
    }
  }
  REQUIRE_THROWS_AS(gra::BuildDirectHelicityCoupling(
      mother, {vector, vector}, {{1, 0}, {0, -1}}, {coupling, -coupling},
      true, true, false, "inconsistent vector orbit", false), std::invalid_argument);
}

// Reject nonphysical scalar spin projections before forming crossed photon keys
TEST_CASE("crossed photon helicity input validates original spin projections",
          "[gra::spin][helicity][photon][validation]") {
  gra::MParticle photon;
  photon.pdg = 22;
  photon.spinX2 = 2;
  photon.P = photon.C = -1;
  gra::MParticle pion;
  pion.pdg = 211;
  pion.spinX2 = 0;
  pion.P = -1;
  gra::MParticle antipion = pion;
  antipion.pdg = -211;
  gra::HelicityChannelMatch match;
  nlohmann::json block = {{"basis", "crossed_helicity"}, {"CP", {true, true}},
                          {"helicity", {{0.0, 0.0, 1.0, 0.0}}}};
  gra::HELMatrix valid;
  REQUIRE_NOTHROW(gra::ParseTwoBodyCouplings(
      valid, block, "scalar photon vertex", "CON_GP.json", photon, {pion, antipion}, match, true,
      gra::spin::VertexContext::SubTUChannelExchange, gra::regge::Signature::Negative, true, 2));
  CHECK(valid.T.FrobNorm2() == Approx(2.0).epsilon(1e-12));
  for (const double projection : {0.1, -0.1, 0.5, 1.0, 1e100,
                                   std::numeric_limits<double>::infinity(),
                                   std::numeric_limits<double>::quiet_NaN()}) {
    for (const std::size_t leg : {0U, 1U}) {
      CAPTURE(projection, leg);
      auto invalid = block;
      invalid["helicity"][0][leg] = projection;
      gra::HELMatrix hc;
      REQUIRE_THROWS_AS(gra::ParseTwoBodyCouplings(
          hc, invalid, "invalid scalar photon vertex", "CON_GP.json", photon, {pion, antipion}, match,
          true, gra::spin::VertexContext::SubTUChannelExchange, gra::regge::Signature::Negative, true, 2),
          std::invalid_argument);
    }
  }
}

// Validate physical branching fractions and finite decay phases at input time
TEST_CASE("decay input rejects unphysical branching ratios and non-finite couplings",
          "[gra::spin][decay][validation]") {
  nlohmann::json block = {{"BR", 0.37}, {"zeta", {{"MP", 0.2}, {"XP", -0.1}, {"GP", 0.4}, {"TP", -0.3}}}, {"g_decay_TP", {-0.3, 0.0, 0.5}}};
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    block["FF_decay"][model] = {{"type", "none"}};
  }
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    auto input = block;
    input["FF_decay"].erase(model);
    gra::HELMatrix hc;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "missing decay form factor", false, model),
                      std::invalid_argument);
  }

  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    gra::HELMatrix hc;
    gra::ParseDecayParameters(hc, block, "model decay phase", false, model);
    CHECK(hc.zeta == Approx(block.at("zeta").at(model).get<double>()));
  }
  for (const nlohmann::json &phases : {nlohmann::json(0.2), nlohmann::json{{"MP", 0.2}},
                                      nlohmann::json{{"MP", 0.2}, {"XP", 0.1}, {"GP", 0.4}, {"TP", "invalid"}}}) {
    auto input = block;
    input["zeta"] = phases;
    gra::HELMatrix hc;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "invalid phase table", false, "MP"), std::invalid_argument);
  }
  {
    auto input = block;
    input.erase("FF_decay");
    gra::HELMatrix hc;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "missing decay form factors", false, "GP"),
                      std::invalid_argument);
  }
  // Reusing a helicity matrix must not carry a prior channel's tensor vertex
  gra::HELMatrix reused;
  gra::ParseDecayParameters(reused, block, "tensor decay", false, "TP");
  REQUIRE(reused.g_decay_TP.size() == 3);
  auto without_tensor = block;
  without_tensor.erase("g_decay_TP");
  gra::ParseDecayParameters(reused, without_tensor, "new decay", false, "TP");
  CHECK(reused.g_decay_TP.empty());
  gra::ParseDecayParameters(reused, block, "tensor decay", false, "TP");
  gra::ParseDecayParameters(reused, {}, "production vertex", true, "TP");
  CHECK(reused.g_decay_TP.empty());
  CHECK_FALSE(reused.BR_set);
  for (const double br : {0.0, 0.37, 1.0}) {
    auto input = block;
    input["BR"] = br;
    gra::HELMatrix hc;
    REQUIRE_NOTHROW(gra::ParseDecayParameters(hc, input, "valid decay", false, "GP"));
    CHECK(hc.BR == Approx(br));
  }
  for (const double br : {-0.2, 1.2, std::numeric_limits<double>::infinity(),
                          std::numeric_limits<double>::quiet_NaN()}) {
    auto input = block;
    input["BR"] = br;
    gra::HELMatrix hc;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "invalid decay BR", false, "GP"),
                      std::invalid_argument);
  }
  for (const double value : {std::numeric_limits<double>::infinity(),
                             std::numeric_limits<double>::quiet_NaN()}) {
    gra::HELMatrix hc;
    auto input = block;
    input["zeta"]["GP"] = value;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "invalid decay phase", false, "GP"),
                      std::invalid_argument);
    input = block;
    input["g_decay_TP"][0] = value;
    REQUIRE_THROWS_AS(gra::ParseDecayParameters(hc, input, "invalid tensor decay", false, "TP"),
                      std::invalid_argument);
  }
}

namespace {

// Build one coherent photon fixture with the immutable TUNE0 cache bound
gra::LORENTZSCALAR MakeToyQEDLTS() {
  gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
  lts.model_cache =
      std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
  return lts;
}

// Convert leg-major negative-first rows to initial-final negative-first rows
MMatrix<std::complex<double>> CanonicalFullPairSpinRowsForTest(
    const MMatrix<std::complex<double>> &leg_major) {
  if (leg_major.size_row() != 16) {
    throw std::invalid_argument(
        "CanonicalFullPairSpinRowsForTest: expected sixteen rows");
  }

  MMatrix<std::complex<double>> canonical(leg_major.size_row(),
                                          leg_major.size_col(), 0.0);
  for (std::size_t i1 = 0; i1 < 2; ++i1) {
    for (std::size_t i2 = 0; i2 < 2; ++i2) {
      for (std::size_t f1 = 0; f1 < 2; ++f1) {
        for (std::size_t f2 = 0; f2 < 2; ++f2) {
          const std::size_t upper_row = 2 * i1 + f1;
          const std::size_t lower_row = 2 * i2 + f2;
          const std::size_t source = 4 * upper_row + lower_row;
          const std::size_t target = 8 * i1 + 4 * i2 + 2 * f1 + f2;
          for (std::size_t col = 0; col < leg_major.size_col(); ++col) {
            canonical[target][col] = leg_major[source][col];
          }
        }
      }
    }
  }
  return canonical;
}

// Convert compact leg-major no-flip rows to canonical negative-first rows
MMatrix<std::complex<double>> CanonicalCompactPairSpinRowsForTest(
    const MMatrix<std::complex<double>> &leg_major) {
  if (leg_major.size_row() != 4) {
    throw std::invalid_argument(
        "CanonicalCompactPairSpinRowsForTest: expected four rows");
  }

  MMatrix<std::complex<double>> canonical(leg_major.size_row(),
                                          leg_major.size_col(), 0.0);
  for (std::size_t i1 = 0; i1 < 2; ++i1) {
    for (std::size_t i2 = 0; i2 < 2; ++i2) {
      const std::size_t source = 2 * i1 + i2;
      const std::size_t target = 2 * i1 + i2;
      for (std::size_t col = 0; col < leg_major.size_col(); ++col) {
        canonical[target][col] = leg_major[source][col];
      }
    }
  }
  return canonical;
}

} // namespace

TEST_CASE("Common binary helicity bases define every pair row mapping",
          "[gra::spin][helicity-basis][pair-layout]") {
  constexpr std::array<int, 2> negative_labels = {-1, 1};
  constexpr std::array<std::array<int, 2>, 4> negative_pairs = {
      {{{-1, -1}}, {{-1, 1}}, {{1, -1}}, {{1, 1}}}};

  CHECK(gra::spin::BinaryHelicityLabelsX2() == negative_labels);
  for (std::size_t pair = 0; pair < 4; ++pair) {
    CHECK(gra::spin::BinaryPairHelicityLabelsX2(pair) == negative_pairs[pair]);
    CHECK(gra::spin::BinaryPairHelicityIndexX2(
              negative_pairs[pair][0], negative_pairs[pair][1]) == pair);
  }

  for (std::size_t initial = 0; initial < 4; ++initial) {
    for (std::size_t final = 0; final < 4; ++final) {
      CHECK(gra::spin::PairHelicityMatrixIndex(final, initial) ==
            4 * final + initial);
      CHECK(gra::spin::PairHelicityTransitionIndex(initial, final) ==
            4 * initial + final);
    }
  }

}

TEST_CASE("Independent pair-row references convert negative-first tokens",
          "[gra::spin][pair-layout][reference]") {
  MMatrix<std::complex<double>> full(16, 1, 0.0);
  for (std::size_t row = 0; row < full.size_row(); ++row) {
    full[row][0] = std::complex<double>(static_cast<double>(row + 1),
                                        -static_cast<double>(row + 1));
  }
  const auto canonical_full = CanonicalFullPairSpinRowsForTest(full);
  constexpr std::array<std::size_t, 16> full_source = {
      0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};
  for (std::size_t row = 0; row < canonical_full.size_row(); ++row) {
    REQUIRE(canonical_full[row][0] == full[full_source[row]][0]);
  }

  MMatrix<std::complex<double>> compact(4, 1, 0.0);
  for (std::size_t row = 0; row < compact.size_row(); ++row) {
    compact[row][0] = std::complex<double>(0.25 + row, 0.5 - row);
  }
  const auto canonical_compact = CanonicalCompactPairSpinRowsForTest(compact);
  constexpr std::array<std::size_t, 4> compact_source = {0, 1, 2, 3};
  for (std::size_t row = 0; row < canonical_compact.size_row(); ++row) {
    REQUIRE(canonical_compact[row][0] == compact[compact_source[row]][0]);
  }
}

TEST_CASE("gra::spin::SpinHalfTransitions returns the physical row order",
          "[gra::spin][transitions]") {
  const std::vector<std::pair<double, double>> expected_full = {
      {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<std::pair<double, double>> expected_no_flip = {{-0.5, -0.5},
                                                                   {0.5, 0.5}};

  REQUIRE(gra::spin::SpinHalfTransitions(false) == expected_full);
  REQUIRE(gra::spin::SpinHalfTransitions(true) == expected_no_flip);
}

// Check finite SU(2) representation dimensions, projections and Casimir sums
TEST_CASE("gra::spin finite representations obey SU(2) projection identities",
          "[gra::spin][representation]") {
  for (int spin_x2 = 0; spin_x2 <= 8; ++spin_x2) {
    const gra::spin::SpinRep rep(spin_x2, "test spin representation");
    const auto projections_x2 = rep.ProjectionsX2();
    const auto projections = rep.Projections();
    REQUIRE(rep.X2() == spin_x2);
    REQUIRE(rep.Dim() == static_cast<std::size_t>(spin_x2 + 1));
    REQUIRE(projections_x2.size() == rep.Dim());
    REQUIRE(projections.size() == rep.Dim());

    double sum = 0.0;
    double sum2 = 0.0;
    for (const auto &i : indices(projections)) {
      CAPTURE(spin_x2, i);
      REQUIRE(rep.IndexX2(projections_x2[i], "doubled projection") == i);
      REQUIRE(rep.Index(projections[i], "physical projection") == i);
      REQUIRE(projections[i] == Approx(0.5 * projections_x2[i]));
      sum += projections[i];
      sum2 += projections[i] * projections[i];
    }
    const double spin = rep.Spin();
    REQUIRE(sum == Approx(0.0).margin(1e-14));
    REQUIRE(sum2 ==
            Approx(spin * (spin + 1.0) * (2.0 * spin + 1.0) / 3.0)
                .margin(1e-13));
  }

  REQUIRE_THROWS_AS(gra::spin::SpinRep(-1, "negative spin"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(
      gra::spin::SpinRep::FromSpin(0.25, "quarter spin"),
      std::invalid_argument);
  const gra::spin::SpinRep spin_one(2, "spin one");
  REQUIRE_THROWS_AS(spin_one.IndexX2(1, "wrong projection parity"),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(spin_one.Index(2.0, "projection outside spin one"),
                    std::invalid_argument);
}

// Check the Jacob-Wick map is orthonormal on the coupled LS subspace
TEST_CASE("gra::spin Jacob-Wick operators form the coupled two-body basis",
          "[gra::spin][jacob-wick]") {
  const gra::spin::SpinRep vector(2, "vector");
  const gra::spin::TwoBodySpin spin(
      vector, gra::spin::SpinBasis(vector), gra::spin::SpinBasis(vector));
  std::vector<MMatrix<std::complex<double>>> operators;
  for (std::size_t l = 0; l <= 3; ++l) {
    for (int two_s = 0; two_s <= 4; ++two_s) {
      if (spin.AllowsLS(l, two_s)) {
        operators.push_back(
            spin.JacobWick(l, two_s, "Jacob-Wick identity"));
      }
    }
  }
  REQUIRE(operators.size() == 7);
  for (const auto &i : indices(operators)) {
    const auto left = operators[i].Flatten();
    for (const auto &j : indices(operators)) {
      const auto right = operators[j].Flatten();
      const auto overlap = gra::InnerProduct(left, right);
      CHECK(overlap.real() == Approx(i == j ? 1.0 : 0.0).margin(1e-13));
      CHECK(overlap.imag() == Approx(0.0).margin(1e-13));
    }
  }

  const gra::spin::TwoBodySpin transverse(
      gra::spin::SpinRep(0, "scalar"),
      gra::spin::SpinBasis(vector, {-2, 2}, "transverse leg 1"),
      gra::spin::SpinBasis(vector, {-2, 2}, "transverse leg 2"));
  const auto scalar = transverse.JacobWick(0, 0, "transverse scalar");
  CHECK(std::abs(scalar[1][0]) == Approx(0.0));
  CHECK(std::abs(scalar[0][1]) == Approx(0.0));
  CHECK(std::abs(scalar[1][1]) == Approx(0.0));
  CHECK(scalar.FrobNorm2() == Approx(2.0 / 3.0).margin(1e-13));
}

TEST_CASE("gra::wigner algebra functions satisfy reference identities",
          "[gra::wigner]") {
  constexpr double tol = 1e-12;

  SECTION("integer spin labels are classified correctly") {
    CHECK(gra::wigner::IsInt(0.0));
    CHECK(gra::wigner::IsInt(2.0));
    CHECK(gra::wigner::IsInt(-2.0));
    CHECK_FALSE(gra::wigner::IsInt(0.5));
    CHECK(gra::wigner::IsInt(1.0 + 1e-10));
    CHECK_FALSE(gra::wigner::IsInt(1.0 + 1e-6));
  }

  SECTION("spin and helicity projection bases are ordered from -s to +s") {
    const auto scalar = gra::spin::SpinProjections(0.0);
    REQUIRE(scalar.size() == 1);
    CHECK(scalar[0] == Approx(0.0));

    const auto spin_half = gra::spin::SpinProjections(0.5);
    REQUIRE(spin_half.size() == 2);
    CHECK(spin_half[0] == Approx(-0.5));
    CHECK(spin_half[1] == Approx(0.5));

    const auto spin_three_half = gra::spin::SpinProjections(1.5);
    REQUIRE(spin_three_half.size() == 4);
    CHECK(spin_three_half[0] == Approx(-1.5));
    CHECK(spin_three_half[1] == Approx(-0.5));
    CHECK(spin_three_half[2] == Approx(0.5));
    CHECK(spin_three_half[3] == Approx(1.5));
    CHECK_THROWS_AS(gra::spin::SpinProjections(0.25), std::invalid_argument);

    gra::HELMatrix basis;
    gra::spin::InitTwoBodyBasis(basis, 1.5, 0.5, 1.0,
                                gra::spin::SpinProjections(1.5), spin_half,
                                gra::spin::SpinProjections(1.0),
                                "mixed-spin basis test");
    REQUIRE(basis.lambda_values.size_row() == 6);
    REQUIRE(basis.lambda_values.size_col() == 2);
    REQUIRE(basis.lambda_idx.size_row() == 6);
    REQUIRE(basis.lambda_idx.size_col() == 2);
    CHECK(basis.lambda_values[0][0] == Approx(-0.5));
    CHECK(basis.lambda_values[0][1] == Approx(-1.0));
    CHECK(basis.lambda_values[2][0] == Approx(-0.5));
    CHECK(basis.lambda_values[2][1] == Approx(1.0));
    CHECK(basis.lambda_values[3][0] == Approx(0.5));
    CHECK(basis.lambda_values[3][1] == Approx(-1.0));
    CHECK(basis.lambda_values[5][0] == Approx(0.5));
    CHECK(basis.lambda_values[5][1] == Approx(1.0));
    CHECK(basis.lambda_idx[0][0] == 0);
    CHECK(basis.lambda_idx[0][1] == 0);
    CHECK(basis.lambda_idx[5][0] == 1);
    CHECK(basis.lambda_idx[5][1] == 2);

    CHECK(gra::spin::ColliderSpinHalfHelicityHarmonic(1, -1, 1, 1) == 1);
    CHECK(gra::spin::ColliderSpinHalfHelicityHarmonic(-1, 1, 1, -1) == -2);
    CHECK(gra::spin::ColliderSpinHalfHardHelicityHarmonic(-1, 1, -1, 1) == -1);
    CHECK(gra::spin::ColliderSpinHalfReciprocitySign(1, -1, 1, 1) ==
          Approx(-1.0));
    CHECK(gra::spin::ColliderSpinHalfReciprocitySign(-1, 1, 1, -1) ==
          Approx(1.0));
    CHECK_THROWS_AS(gra::spin::ColliderSpinHalfHelicityHarmonic(0, 1, 1, 1),
                    std::invalid_argument);
    CHECK_THROWS_AS(gra::spin::ColliderSpinHalfHardHelicityHarmonic(0, 1, 1, 1),
                    std::invalid_argument);
  }

  SECTION(
      "Clebsch-Gordan and Wigner-3j selection rules reject forbidden tuples") {
    CHECK(gra::math::IsZero(gra::wigner::CG(0.5, 0.5, 0.5, -0.5, 1.0, 1.0)));
    CHECK(gra::math::IsZero(gra::wigner::CG(1.0, 1.0, 0.0, 0.0, 3.0, 0.0)));
    CHECK(gra::math::IsZero(gra::wigner::CG(1.0, 1.0, 0.5, -0.5, 1.0, 0.0)));
    CHECK(gra::math::IsZero(gra::wigner::W3j(1.0, 1.0, 0.0, 1.0, 0.0, 0.0)));
    CHECK(gra::math::IsZero(gra::wigner::W3j(1.0, 1.0, 3.0, 0.0, 0.0, 0.0)));
    CHECK(gra::math::IsZero(gra::wigner::W3j(1.0, 1.0, 1.0, 0.5, -0.5, 0.0)));
  }

  SECTION("Wigner-3j symbols match closed-form values") {
    CHECK(gra::wigner::W3j(0.0, 0.0, 0.0, 0.0, 0.0, 0.0) ==
          Approx(1.0).epsilon(tol));
    CHECK(gra::wigner::W3j(1.0, 1.0, 0.0, 0.0, 0.0, 0.0) ==
          Approx(-1.0 / std::sqrt(3.0)).epsilon(tol));
    CHECK(gra::wigner::W3j(0.5, 0.5, 0.0, 0.5, -0.5, 0.0) ==
          Approx(1.0 / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::W3j(1.0, 1.0, 1.0, 1.0, -1.0, 0.0) ==
          Approx(1.0 / std::sqrt(6.0)).epsilon(tol));

    const double swapped = gra::wigner::W3j(1.0, 1.0, 1.0, -1.0, 1.0, 0.0);
    CHECK(
        swapped ==
        Approx(-gra::wigner::W3j(1.0, 1.0, 1.0, 1.0, -1.0, 0.0)).epsilon(tol));
    CHECK(gra::wigner::W3j(1.0, 1.0, 3.0, 0.0, 0.0, 0.0) ==
          Approx(0.0).margin(tol));
  }

  SECTION("Clebsch-Gordan coefficients match singlet and triplet spin-half "
          "states") {
    CHECK(gra::wigner::CG(0.5, 0.5, 0.5, 0.5, 1.0, 1.0) ==
          Approx(1.0).epsilon(tol));
    CHECK(gra::wigner::CG(0.5, 0.5, 0.5, -0.5, 1.0, 0.0) ==
          Approx(1.0 / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::CG(0.5, 0.5, -0.5, 0.5, 1.0, 0.0) ==
          Approx(1.0 / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::CG(0.5, 0.5, 0.5, -0.5, 0.0, 0.0) ==
          Approx(1.0 / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::CG(0.5, 0.5, -0.5, 0.5, 0.0, 0.0) ==
          Approx(-1.0 / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::CG(0.5, 0.5, 0.5, 0.5, 0.0, 0.0) ==
          Approx(0.0).margin(tol));
  }

  SECTION("Clebsch-Gordan coefficients form a unitary coupled basis") {
    const double j1 = 1.0;
    const double j2 = 0.5;
    const auto m1_values = gra::spin::SpinProjections(j1);
    const auto m2_values = gra::spin::SpinProjections(j2);

    for (const double m1 : m1_values) {
      for (const double m2 : m2_values) {
        for (const double m1_prime : m1_values) {
          for (const double m2_prime : m2_values) {
            double sum = 0.0;
            for (const double J : {0.5, 1.5}) {
              for (const double M : gra::spin::SpinProjections(J)) {
                sum += gra::wigner::CG(j1, j2, m1, m2, J, M) *
                       gra::wigner::CG(j1, j2, m1_prime, m2_prime, J, M);
              }
            }
            const double expected =
                (m1 == m1_prime && m2 == m2_prime) ? 1.0 : 0.0;
            CAPTURE(m1, m2, m1_prime, m2_prime);
            CHECK(sum == Approx(expected).margin(2e-13));
          }
        }
      }
    }
  }

  SECTION("Wigner-3j permutation and sign reversal symmetries hold") {
    const double base = gra::wigner::W3j(2.0, 1.5, 1.5, 1.0, -0.5, -0.5);
    const double phase = -1.0;
    CHECK(gra::wigner::W3j(1.5, 2.0, 1.5, -0.5, 1.0, -0.5) ==
          Approx(phase * base).margin(tol));
    CHECK(gra::wigner::W3j(2.0, 1.5, 1.5, -1.0, 0.5, 0.5) ==
          Approx(phase * base).margin(tol));
    CHECK(gra::wigner::W3j(1.5, 1.5, 2.0, -0.5, -0.5, 1.0) ==
          Approx(base).margin(tol));
  }

  SECTION("Wigner small-d matches spin-half and spin-one closed forms") {
    const double theta = 0.73;
    const double c = std::cos(theta / 2.0);
    const double s = std::sin(theta / 2.0);
    CHECK(gra::wigner::d(theta, 0.5, 0.5, 0.5) == Approx(c).epsilon(tol));
    CHECK(gra::wigner::d(theta, 0.5, -0.5, 0.5) == Approx(s).epsilon(tol));
    CHECK(gra::wigner::d(theta, -0.5, 0.5, 0.5) == Approx(-s).epsilon(tol));
    CHECK(gra::wigner::d(theta, -0.5, -0.5, 0.5) == Approx(c).epsilon(tol));

    CHECK(gra::wigner::d(theta, 1.0, 1.0, 1.0) == Approx(c * c).epsilon(tol));
    CHECK(gra::wigner::d(theta, 1.0, 0.0, 1.0) ==
          Approx(std::sin(theta) / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::d(theta, 0.0, 1.0, 1.0) ==
          Approx(-std::sin(theta) / std::sqrt(2.0)).epsilon(tol));
    CHECK(gra::wigner::d(theta, 0.0, 0.0, 1.0) ==
          Approx(std::cos(theta)).epsilon(tol));
    CHECK(gra::wigner::d(theta, 1.0, -1.0, 1.0) == Approx(s * s).epsilon(tol));
    CHECK(gra::wigner::d(theta, 2.0, 0.0, 1.0) == Approx(0.0).margin(tol));
    CHECK_THROWS_AS(gra::wigner::d(theta, 0.0, 0.0, 0.25),
                    std::invalid_argument);
    CHECK_THROWS_AS(
        gra::wigner::d(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 1.0),
        std::invalid_argument);
    CHECK_THROWS_AS(gra::wigner::D(theta,
                                   std::numeric_limits<double>::infinity(), 0.0,
                                   0.0, 1.0),
                    std::invalid_argument);
    CHECK_THROWS_AS(gra::wigner::DMatrix(
                        1.0, theta, std::numeric_limits<double>::quiet_NaN()),
                    std::invalid_argument);
    CHECK_THROWS_AS(gra::wigner::d(theta, 0.0, 0.0, 100.0),
                    std::invalid_argument);
  }

  SECTION("Regge functions reproduce ordinary integer spin values") {
    const auto w3j = gra::wigner::W3jRegge(2.0, 1.0, 2.0, 1.0, -1.0, 0.0);
    RequireComplexNear(w3j, gra::wigner::W3j(2.0, 1.0, 2.0, 1.0, -1.0, 0.0),
                       tol);

    const auto cg = gra::wigner::CGRegge(2.0, 1.0, 1.0, -1.0, 2.0, 0.0);
    RequireComplexNear(cg, gra::wigner::CG(2.0, 1.0, 1.0, -1.0, 2.0, 0.0), tol);

    const double theta = 0.73;
    const double phi = -0.41;
    RequireComplexNear(gra::wigner::dRegge(theta, 1.0, -1.0, 2.0),
                       gra::wigner::d(theta, 1.0, -1.0, 2.0), tol);
    RequireComplexNear(gra::wigner::DRegge(theta, phi, 1.0, -1.0, 2.0),
                       gra::wigner::D(theta, phi, 1.0, -1.0, 2.0), tol);
  }

  SECTION("Regge functions are continuous through integer trajectory spin") {
    constexpr double epsilon = 1.0e-6;
    const double theta = 0.73;
    const std::complex<double> cg_integer =
        gra::wigner::CG(2.0, 1.0, 1.0, -1.0, 2.0, 0.0);
    const auto cg_minus =
        gra::wigner::CGRegge(2.0 - epsilon, 1.0, 1.0, -1.0, 2.0, 0.0);
    const auto cg_plus =
        gra::wigner::CGRegge(2.0 + epsilon, 1.0, 1.0, -1.0, 2.0, 0.0);
    const std::complex<double> d_integer = gra::wigner::d(theta, 1.0, 0.0, 2.0);
    const auto d_minus = gra::wigner::dRegge(theta, 1.0, 0.0, 2.0 - epsilon);
    const auto d_plus = gra::wigner::dRegge(theta, 1.0, 0.0, 2.0 + epsilon);
    RequireComplexNear(cg_minus, cg_plus, 2.0e-5);
    RequireComplexNear(cg_minus, cg_integer, 2.0e-5);
    RequireComplexNear(cg_plus, cg_integer, 2.0e-5);
    RequireComplexNear(d_minus, d_plus, 2.0e-5);
    RequireComplexNear(d_minus, d_integer, 2.0e-5);
    RequireComplexNear(d_plus, d_integer, 2.0e-5);
  }

  SECTION("integer m rows stay analytic at half-integer trajectory spin") {
    constexpr double epsilon = 1.0e-7;
    const double theta = 0.73;
    const auto d_minus = gra::wigner::dRegge(theta, 1.0, 0.0, 1.5 - epsilon);
    const auto d_at = gra::wigner::dRegge(theta, 1.0, 0.0, 1.5);
    const auto d_plus = gra::wigner::dRegge(theta, 1.0, 0.0, 1.5 + epsilon);
    CHECK(std::isfinite(d_at.real()));
    CHECK(std::isfinite(d_at.imag()));
    RequireComplexNear(d_minus, d_at, 2.0e-5);
    RequireComplexNear(d_plus, d_at, 2.0e-5);

    const auto cg_minus =
        gra::wigner::CGRegge(2.0, 1.0, 0.0, 0.0, 1.5 - epsilon, 0.0);
    const auto cg_at = gra::wigner::CGRegge(2.0, 1.0, 0.0, 0.0, 1.5, 0.0);
    const auto cg_plus =
        gra::wigner::CGRegge(2.0, 1.0, 0.0, 0.0, 1.5 + epsilon, 0.0);
    CHECK(std::isfinite(cg_at.real()));
    CHECK(std::isfinite(cg_at.imag()));
    RequireComplexNear(cg_minus, cg_at, 2.0e-5);
    RequireComplexNear(cg_plus, cg_at, 2.0e-5);
  }

  SECTION("continued Wigner rotations approach every retained integer row") {
    constexpr double epsilon = 1.0e-10;
    for (const double theta : {0.31, 0.93, 2.40}) {
      for (int J = 0; J <= 4; ++J) {
        for (int m = -J; m <= J; ++m) {
          for (int mp = -J; mp <= J; ++mp) {
            CAPTURE(theta, J, m, mp);
            const double ordinary = gra::wigner::d(theta, m, mp, J);
            RequireComplexNear(gra::wigner::dRegge(theta, m, mp, J - epsilon),
                               ordinary, 2.0e-5);
            RequireComplexNear(gra::wigner::dRegge(theta, m, mp, J + epsilon),
                               ordinary, 2.0e-5);
          }
        }
      }
    }
  }

  SECTION("continued continuum Clebsch-Gordan rows approach integer spin") {
    constexpr double epsilon = 1.0e-12;
    for (int L = 0; L <= 4; ++L) {
      for (int S = 0; S <= 3; ++S) {
        for (int J = std::abs(L - S); J <= L + S; ++J) {
          for (int lambda = -S; lambda <= S; ++lambda) {
            CAPTURE(L, S, J, lambda);
            const double ordinary =
                gra::wigner::CG(L, S, 0.0, lambda, J, lambda);
            if (J > 0) {
              RequireComplexNear(
                  gra::wigner::CGRegge(L, S, 0.0, lambda, J - epsilon, lambda),
                  ordinary, 3.0e-5);
            }
            RequireComplexNear(
                gra::wigner::CGRegge(L, S, 0.0, lambda, J + epsilon, lambda),
                ordinary, 3.0e-5);
          }
        }
      }
    }
  }

  SECTION("Zero projection continuation agrees with the general Racah series") {
    for (const double j1 : {-0.3125, 0.4375, 1.0625, 1.5, 2.3125}) {
      for (const double j2 : {0.25, 0.8125, 1.09375, 2.125}) {
        for (const double j3 : {0.0, 0.625, 1.0, 2.0, 4.0}) {
          const long double triangle = std::tgamma(j1 + j2 - j3 + 1.0L) * std::tgamma(j1 - j2 + j3 + 1.0L) *
                                       std::tgamma(-j1 + j2 + j3 + 1.0L) / std::tgamma(j1 + j2 + j3 + 2.0L);
          const long double legs = std::tgamma(j1 + 1.0L) * std::tgamma(j2 + 1.0L) * std::tgamma(j3 + 1.0L);
          const long double denominator =
              std::tgamma(j1 + j2 - j3 + 1.0L) * std::tgamma(j1 + 1.0L) * std::tgamma(j2 + 1.0L);
          const double series =
              gra::math::RegularizedHyper3F2Unit(-j1 - j2 + j3, -j1, -j2, 1.0 - j2 + j3, 1.0 - j1 + j3);
          const auto expected = std::exp(gra::math::zi * gra::math::PI * (j1 - j2)) *
                                std::sqrt(std::complex<double>(static_cast<double>(triangle * legs * legs), 0.0)) *
                                static_cast<double>(series / denominator);
          CAPTURE(j1, j2, j3);
          RequireComplexNear(gra::wigner::W3jRegge(j1, j2, j3, 0.0, 0.0, 0.0), expected, 2.0e-12);
        }
      }
    }
  }

  SECTION("Regge selection zeros remain distinct from numerical failures") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    const double theta = 0.73;
    REQUIRE_THROWS_AS(gra::wigner::W3jRegge(nan, 1.0, 1.0, 0.0, 0.0, 0.0), gra::AmplitudeFailure);
    RequireComplexNear(gra::wigner::CGRegge(1.0, 1.0, 0.0, 1.0, 1.0, 0.0), 0.0,
                       0.0);
    RequireComplexNear(gra::wigner::dRegge(theta, 0.5, 0.0, 1.2), 0.0, 0.0);
    REQUIRE_FALSE(std::isfinite(std::norm(gra::wigner::DRegge(theta, nan, 0.0, 0.0, 1.2))));
  }

  SECTION("Regularized hypergeometric sums cross denominator Gamma poles") {
    CHECK(gra::math::RegularizedHyper2F1(-1.0, 2.0, 0.0, 0.25) ==
          Approx(-0.5).epsilon(tol));
    CHECK(gra::math::RegularizedHyper3F2Unit(-1.0, 1.0, 1.0, 0.0, 1.0) ==
          Approx(-1.0).epsilon(tol));
  }

  SECTION("Batched Wigner rows match scalar evaluation at edge angles") {
    const std::array<double, 7> angles = {0.0,
                                          1.0e-12,
                                          0.73,
                                          0.5 * gra::math::PI,
                                          gra::math::PI - 1.0e-12,
                                          gra::math::PI,
                                          -0.61};
    for (const double J : {0.0, 0.5, 1.0, 1.5, 2.0, 3.0}) {
      std::vector<double> m_values = gra::spin::SpinProjections(J);
      const std::vector<double> mp_values = gra::spin::SpinProjections(J);
      m_values.push_back(J + 1.0);
      const gra::wigner::Rotation rotation(m_values, mp_values, J);
      for (const double theta : angles) {
        CAPTURE(J, theta);
        const auto batched = rotation.Evaluate(theta);
        REQUIRE(batched.size_row() == m_values.size());
        REQUIRE(batched.size_col() == mp_values.size());
        for (std::size_t row = 0; row < m_values.size(); ++row) {
          for (std::size_t col = 0; col < mp_values.size(); ++col) {
            const double scalar =
                gra::wigner::d(theta, m_values[row], mp_values[col], J);
            CHECK(std::abs(batched[row][col] - scalar) < 2.0e-13);
          }
        }
      }
    }

    CHECK_FALSE(gra::wigner::Rotation({0.0}, {0.0}, 1.0)
                    .Evaluate(std::numeric_limits<double>::quiet_NaN()).IsFinite());
    CHECK_THROWS_AS(gra::wigner::dRows(0.3, {0.0}, {0.0}, 0.25),
                    std::invalid_argument);
  }

  SECTION("Wigner-D phase convention and DMatrix unitarity hold") {
    const double theta = 0.41;
    const double phi = -0.83;
    const auto D_element = gra::wigner::D(theta, phi, 1.0, -1.0, 1.0);
    const std::complex<double> expected =
        gra::wigner::d(theta, 1.0, -1.0, 1.0) *
        std::exp(std::complex<double>(0.0, -phi));
    RequireComplexNear(D_element, expected, tol);

    const auto D = gra::wigner::DMatrix(1.0, theta, phi);
    REQUIRE(D.size_row() == 3);
    REQUIRE(D.size_col() == 3);
    for (std::size_t i = 0; i < D.size_col(); ++i) {
      for (std::size_t j = 0; j < D.size_col(); ++j) {
        std::complex<double> overlap = 0.0;
        for (std::size_t k = 0; k < D.size_row(); ++k) {
          overlap += std::conj(D[k][i]) * D[k][j];
        }
        if (i == j) {
          RequireComplexNear(overlap, 1.0, 1e-11);
        } else {
          CHECK(std::abs(overlap) < 1e-11);
        }
      }
    }

    const auto D_half = gra::wigner::DMatrix(0.5, theta, phi);
    REQUIRE(D_half.size_row() == 2);
    REQUIRE(D_half.size_col() == 2);

    // Stored rows are the first API argument m and columns are mp, so this is
    // the transpose of the conventional D*_(mp,m) index layout
    const double ch = std::cos(theta / 2.0);
    const double sh = std::sin(theta / 2.0);
    const std::complex<double> phase_minus =
        std::exp(-0.5 * gra::math::zi * phi);
    const std::complex<double> phase_plus = std::exp(0.5 * gra::math::zi * phi);
    RequireComplexNear(D_half[0][0], ch * phase_minus, tol);
    RequireComplexNear(D_half[0][1], -sh * phase_plus, tol);
    RequireComplexNear(D_half[1][0], sh * phase_minus, tol);
    RequireComplexNear(D_half[1][1], ch * phase_plus, tol);

    const double phi1 = 0.37;
    const double phi2 = -0.81;
    const auto z1 = gra::wigner::DMatrix(0.5, 0.0, phi1);
    const auto z2 = gra::wigner::DMatrix(0.5, 0.0, phi2);
    const auto zsum = gra::wigner::DMatrix(0.5, 0.0, phi1 + phi2);
    RequireMatrixNear(z1 * z2, zsum, 1e-12);

    const auto full_turn = gra::wigner::DMatrix(0.5, 0.0, 2.0 * gra::math::PI);
    RequireComplexNear(full_turn[0][0], -1.0, tol);
    RequireComplexNear(full_turn[1][1], -1.0, tol);
    RequireComplexNear(full_turn[0][1], 0.0, tol);
    RequireComplexNear(full_turn[1][0], 0.0, tol);
  }

  SECTION("Wigner small-d matrices compose rotations about one axis") {
    const double theta1 = 0.37;
    const double theta2 = -0.81;
    for (const double J : {0.5, 1.0, 1.5, 2.0}) {
      const auto projections = gra::spin::SpinProjections(J);
      const auto d1 = gra::wigner::dRows(theta1, projections, projections, J);
      const auto d2 = gra::wigner::dRows(theta2, projections, projections, J);
      const auto dsum =
          gra::wigner::dRows(theta1 + theta2, projections, projections, J);
      CAPTURE(J);
      REQUIRE((d1 * d2).IsApprox(dsum, 2e-12));
    }
  }

  SECTION("fDecayMatrix is exactly the Jacob-Wick Wigner-D row times T") {
    gra::HELMatrix hel;
    gra::spin::InitTwoBodyBasis(hel, 1.0, 0.0, 0.0,
                                gra::spin::SpinProjections(1.0), {0.0}, {0.0},
                                "Jacob-Wick decay matrix test");
    hel.T =
        MMatrix<std::complex<double>>(1, 1, std::complex<double>(2.0, -0.5));

    const double theta = 0.52;
    const double phi = 0.31;
    const auto f = gra::spin::fDecayMatrix(hel, theta, phi);
    REQUIRE(hel.jw_rotation.ready);
    const auto cached = gra::spin::fDecayMatrix(hel, theta, phi);
    RequireMatrixNear(cached, f, tol);
    auto invalid_cache = hel;
    invalid_cache.lambda_idx = MMatrix<std::size_t>(1, 1, 0);
    REQUIRE_THROWS_AS(gra::spin::InitJWRotation(invalid_cache), std::invalid_argument);
    REQUIRE(f.size_row() == 1);
    REQUIRE(f.size_col() == 3);
    for (std::size_t j = 0; j < hel.Jz_values.size(); ++j) {
      const std::complex<double> expected =
          gra::wigner::D(theta, phi, 0.0, hel.Jz_values[j], hel.J) *
          hel.T[0][0];
      RequireComplexNear(f[0][j], expected, tol);
    }

    hel.Jz_values.pop_back();
    const auto selected = gra::spin::fDecayMatrix(hel, theta, phi);
    REQUIRE(selected.size_row() == 1);
    REQUIRE(selected.size_col() == 2);

    gra::HELMatrix scalar = hel;
    scalar.J = 0.0;
    scalar.Jz_values = {0.0};
    CHECK_THROWS_AS(gra::spin::fDecayMatrix(
                        scalar, std::numeric_limits<double>::quiet_NaN(), phi),
                    gra::AmplitudeFailure);
  }
}

TEST_CASE("gra::spin:: Contract emits canonical full proton-pair hard rows",
          "[gra::spin][pair-layout]") {
  const std::array<std::size_t, 16> expected_rows = {
      0, 1, 4, 5, 2, 3, 6, 7, 8, 9, 12, 13, 10, 11, 14, 15};
  const std::array<std::pair<int, int>, 4> negative_first_leg = {
      std::pair{-1, -1}, std::pair{-1, 1}, std::pair{1, -1}, std::pair{1, 1}};
  const MMatrix<std::complex<double>> central(1, 1, 1.0);
  const MMatrix<std::complex<double>> sub_t(1, 1,
                                            std::complex<double>(2.0, -0.5));
  const MMatrix<std::complex<double>> sub_u(1, 1,
                                            std::complex<double>(-0.3, 1.7));
  const auto destination_rows =
      gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(4, 4);
  REQUIRE(destination_rows.size() == 16);

  for (std::size_t upper_row = 0; upper_row < 4; ++upper_row) {
    for (std::size_t lower_row = 0; lower_row < 4; ++lower_row) {
      MMatrix<std::complex<double>> upper(4, 1, 0.0);
      MMatrix<std::complex<double>> lower(4, 1, 0.0);
      upper[upper_row][0] = 1.0;
      lower[lower_row][0] = 1.0;

      const std::size_t leg_major_row = 4 * upper_row + lower_row;
      const std::size_t expected_row = expected_rows[leg_major_row];
      REQUIRE(destination_rows[leg_major_row] == expected_row);
      const auto resonance = gra::spin::Contract(upper, lower, central);
      const auto continuum = std::pair{gra::spin::Contract(upper, lower, sub_t), gra::spin::Contract(upper, lower, sub_u)};

      const auto canonical_index = [](const int helicity) {
        return gra::spin::BinaryHelicityIndexX2(helicity);
      };
      const std::size_t i1 =
          canonical_index(negative_first_leg[upper_row].first);
      const std::size_t f1 =
          canonical_index(negative_first_leg[upper_row].second);
      const std::size_t i2 =
          canonical_index(negative_first_leg[lower_row].first);
      const std::size_t f2 =
          canonical_index(negative_first_leg[lower_row].second);
      REQUIRE(gra::spin::CanonicalProtonPairSpinLayout::HardRow(
                  i1, i2, f1, f2) == expected_row);
      REQUIRE(
          gra::spin::CanonicalProtonPairSpinLayout::FromNegativeFirstLegRows(
              upper_row, lower_row, 4, 4) == expected_row);

      for (std::size_t row = 0; row < 16; ++row) {
        const bool selected = row == expected_row;
        RequireComplexNear(resonance[row][0], selected ? 1.0 : 0.0);
        RequireComplexNear(continuum.first[row][0],
                           selected ? sub_t[0][0] : 0.0);
        RequireComplexNear(continuum.second[row][0],
                           selected ? sub_u[0][0] : 0.0);
      }
    }
  }
}

TEST_CASE("gra::spin:: Contract emits canonical compact no-flip pair rows",
          "[gra::spin][pair-layout]") {
  const std::array<std::size_t, 4> expected_rows = {0, 1, 2, 3};
  const MMatrix<std::complex<double>> central(1, 1,
                                              std::complex<double>(0.4, 0.8));
  const auto destination_rows =
      gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(2, 2);
  REQUIRE(destination_rows.size() == 4);

  for (std::size_t upper_row = 0; upper_row < 2; ++upper_row) {
    for (std::size_t lower_row = 0; lower_row < 2; ++lower_row) {
      MMatrix<std::complex<double>> upper(2, 1, 0.0);
      MMatrix<std::complex<double>> lower(2, 1, 0.0);
      upper[upper_row][0] = 1.0;
      lower[lower_row][0] = 1.0;

      const std::size_t expected_row = expected_rows[2 * upper_row + lower_row];
      REQUIRE(destination_rows[2 * upper_row + lower_row] == expected_row);
      const auto actual = gra::spin::Contract(upper, lower, central);
      REQUIRE(gra::spin::CanonicalProtonPairSpinLayout::CompactNoFlipHardRow(
                  upper_row, lower_row) == expected_row);
      REQUIRE(
          gra::spin::CanonicalProtonPairSpinLayout::FromNegativeFirstLegRows(
              upper_row, lower_row, 2, 2) == expected_row);

      for (std::size_t row = 0; row < 4; ++row) {
        RequireComplexNear(actual[row][0],
                           row == expected_row ? central[0][0] : 0.0);
      }
    }
  }
}

TEST_CASE("gra::spin:: Blind preserves forward transition sections",
          "[gra::spin][pair-layout][blind]") {
  using Complex = std::complex<double>;
  const std::vector<Complex> upper = {Complex(1.0, 0.0), Complex(0.2, -0.7),
                                      Complex(-0.4, 0.3), Complex(0.8, 0.1)};
  const std::vector<Complex> lower = {Complex(0.9, -0.2), Complex(-0.3, 0.6),
                                      Complex(0.5, 0.4), Complex(-0.7, -0.1)};
  const std::size_t spin_states = 3;
  const Complex scale(0.6, -0.2);
  const auto actual = gra::spin::Blind(upper, lower, spin_states, scale);
  const auto destination =
      gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(4, 4);
  const auto transition =
      gra::MappedKroneckerProduct(upper, lower, destination);

  REQUIRE(actual.size_row() == 16 * spin_states);
  REQUIRE(actual.size_col() == spin_states);
  for (const auto &row : gra::aux::indices(transition)) {
    for (std::size_t h = 0; h < spin_states; ++h) {
      for (std::size_t col = 0; col < spin_states; ++col) {
        const Complex expected =
            col == h ? scale * transition[row] / std::sqrt(3.0) : 0.0;
        RequireComplexNear(actual[row * spin_states + h][col], expected,
                           1.0e-13);
      }
    }
  }
}

// Check the fused contraction against an independent four index sum
TEST_CASE("gra::spin:: Contract preserves dense helicity phases and channels",
          "[gra::spin][pair-layout][contraction]") {
  using Complex = std::complex<double>;
  MMatrix<Complex> upper(4, 2, 0.0);
  MMatrix<Complex> lower(4, 3, 0.0);
  MMatrix<Complex> central(6, 3, 0.0);
  MMatrix<Complex> sub_t(6, 3, 0.0);
  MMatrix<Complex> sub_u(6, 3, 0.0);

  for (std::size_t row = 0; row < upper.size_row(); ++row) {
    for (std::size_t col = 0; col < upper.size_col(); ++col) {
      upper[row][col] = Complex(0.31 * (row + 1) + 0.07 * col,
                                -0.13 * row + 0.19 * (col + 1));
    }
  }
  for (std::size_t row = 0; row < lower.size_row(); ++row) {
    for (std::size_t col = 0; col < lower.size_col(); ++col) {
      lower[row][col] = Complex(-0.17 * (row + 1) + 0.11 * col,
                                0.23 * row + 0.05 * (col + 1));
    }
  }
  for (std::size_t row = 0; row < central.size_row(); ++row) {
    for (std::size_t col = 0; col < central.size_col(); ++col) {
      central[row][col] = Complex(0.09 * (row + 1) + 0.04 * col,
                                  -0.08 * row + 0.03 * (col + 1));
      sub_t[row][col] =
          central[row][col] * Complex(0.7 + 0.02 * row, -0.3 + 0.01 * col);
      sub_u[row][col] =
          central[row][col] * Complex(-0.2 + 0.03 * col, 0.8 + 0.01 * row);
    }
  }

  const Complex scale(0.63, -0.41);
  const auto resonance = gra::spin::Contract(upper, lower, central, scale);
  const auto continuum = std::pair{gra::spin::Contract(upper, lower, sub_t), gra::spin::Contract(upper, lower, sub_u)};
  MMatrix<Complex> reference_resonance(16, 3, 0.0);
  MMatrix<Complex> reference_t(16, 3, 0.0);
  MMatrix<Complex> reference_u(16, 3, 0.0);

  for (std::size_t upper_row = 0; upper_row < 4; ++upper_row) {
    for (std::size_t lower_row = 0; lower_row < 4; ++lower_row) {
      const std::size_t i1 = upper_row / 2;
      const std::size_t f1 = upper_row % 2;
      const std::size_t i2 = lower_row / 2;
      const std::size_t f2 = lower_row % 2;
      const std::size_t hard_row =
          gra::spin::CanonicalProtonPairSpinLayout::HardRow(i1, i2, f1, f2);
      for (std::size_t x = 0; x < central.size_col(); ++x) {
        for (std::size_t m = 0; m < upper.size_col(); ++m) {
          for (std::size_t n = 0; n < lower.size_col(); ++n) {
            const Complex forward = upper[upper_row][m] * lower[lower_row][n];
            const std::size_t central_row = m * lower.size_col() + n;
            reference_resonance[hard_row][x] +=
                scale * forward * central[central_row][x];
            reference_t[hard_row][x] += forward * sub_t[central_row][x];
            reference_u[hard_row][x] += forward * sub_u[central_row][x];
          }
        }
      }
    }
  }

  for (std::size_t row = 0; row < resonance.size_row(); ++row) {
    for (std::size_t col = 0; col < resonance.size_col(); ++col) {
      REQUIRE(std::abs(resonance[row][col]) > 0.0);
      REQUIRE(std::abs(continuum.first[row][col]) > 0.0);
      REQUIRE(std::abs(continuum.second[row][col]) > 0.0);
      RequireComplexNear(resonance[row][col], reference_resonance[row][col],
                         1e-12);
      RequireComplexNear(continuum.first[row][col], reference_t[row][col],
                         1e-12);
      RequireComplexNear(continuum.second[row][col], reference_u[row][col],
                         1e-12);
      CHECK(std::abs(continuum.first[row][col] - continuum.second[row][col]) >
            1e-12);
    }
  }
}

TEST_CASE("Canonical proton-pair layout rejects non-proton source rows",
          "[gra::spin][pair-layout][validation]") {
  using Layout = gra::spin::CanonicalProtonPairSpinLayout;
  CHECK_THROWS_AS(Layout::FromNegativeFirstLegRows(4, 0, 4, 4),
                  std::invalid_argument);
  CHECK_THROWS_AS(Layout::FromNegativeFirstLegRows(0, 2, 2, 2),
                  std::invalid_argument);
  CHECK_THROWS_AS(Layout::FromNegativeFirstLegRows(0, 0, 3, 3),
                  std::invalid_argument);
  CHECK_THROWS_AS(Layout::FromNegativeFirstLegRows(0, 0, 4, 2),
                  std::invalid_argument);
  CHECK_THROWS_AS(Layout::KroneckerDestinationRows(3, 3),
                  std::invalid_argument);
  CHECK_THROWS_AS(Layout::KroneckerDestinationRows(4, 2),
                  std::invalid_argument);
}

TEST_CASE(
    "Gamma-gamma MG5 W decay chains consume the GRANIITTI four-body cascade",
    "[MG2GRA][gamma-gamma][ww][cascade]") {
  const auto pdg_table = LoadedPDGTable();
  gra::AMP_MG5_yy_ww amplitude;
  const auto particles = amplitude.Particles();
  // Build one stable decay-tree leaf
  const auto stable_leaf = [&pdg_table, &particles](int pdg, const gra::M4Vec &p4) {
    gra::MDecayBranch leaf;
    leaf.p = pdg_table.FindByPDG(pdg);
    leaf.p.mass = particles.at(std::abs(pdg)).mass;
    leaf.p4 = p4;
    return leaf;
  };

  const double w_mass = particles.at(24).mass;
  const double w_energy = 200.0;
  const double w_momentum = std::sqrt(pow2(w_energy) - pow2(w_mass));
  const gra::M4Vec wplus_p4(w_momentum, 0.0, 0.0, w_energy);
  const gra::M4Vec wminus_p4(-w_momentum, 0.0, 0.0, w_energy);

  struct LeptonFlavor {
    int charged_abs;
    double mass;
  };
  const std::array<LeptonFlavor, 3> flavors = {LeptonFlavor{11, particles.at(11).mass},
                                               LeptonFlavor{13, particles.at(13).mass},
                                               LeptonFlavor{15, particles.at(15).mass}};

  // Build one physical flavour-diagonal leptonic W decay
  auto w_decay_branch = [&](bool positive_w, const LeptonFlavor &flavor,
                            const gra::M4Vec &w_p4, double theta, double phi) {
    const auto decay =
        TwoBodyRestKinematics(w_mass, flavor.mass, 0.0, theta, phi);
    gra::MDecayBranch branch;
    branch.p.pdg = positive_w ? 24 : -24;
    branch.p.mass = w_mass;
    branch.p.width = particles.at(24).width;
    branch.p4 = w_p4;
    const int charged_pdg =
        positive_w ? -flavor.charged_abs : flavor.charged_abs;
    const int neutrino_pdg =
        positive_w ? flavor.charged_abs + 1 : -(flavor.charged_abs + 1);
    branch.legs = {
        stable_leaf(charged_pdg, BoostFromRestFrame(decay[0], w_p4)),
        stable_leaf(neutrino_pdg, BoostFromRestFrame(decay[1], w_p4))};
    return branch;
  };

  // Build one full WW decay chain with independently chosen lepton flavours
  auto ww_decay_chain = [&](const LeptonFlavor &plus_flavor,
                            const LeptonFlavor &minus_flavor) {
    gra::LORENTZSCALAR lts;
    lts.model_cache =
        std::make_shared<gra::MModelCache>(gra::MModelTune::Load(modelfile));
    lts.q1 = gra::M4Vec(0.0, 0.0, w_energy, w_energy);
    lts.q2 = gra::M4Vec(0.0, 0.0, -w_energy, w_energy);
    lts.process.root_decay_mode = gra::RootDecayMode::Physical;
    lts.decaytree = {
        w_decay_branch(true, plus_flavor, wplus_p4, 0.73, 0.29),
        w_decay_branch(false, minus_flavor, wminus_p4, 1.18, -0.61)};
    return lts;
  };

  const gra::ProcessDescriptor description{"generated", "MG5", "yy", "", 1};
  gra::MGeneratedPhotonProc process("yy", "WW", description, "MG5_YY_WW");
  auto generated_chain = ww_decay_chain(flavors[0], flavors[0]);
  REQUIRE(process.MatchProcess(generated_chain.decaytree).has_value());
  const auto generated_result = amplitude.Evaluate(generated_chain, 0.0, false);
  REQUIRE(generated_result.Valid());
  REQUIRE(generated_result.amp2 > 0.0);
  REQUIRE(generated_chain.hamp.size() == 64);

  // Accept every generated physical charged-lepton flavour combination
  for (const auto &plus_flavor : flavors) {
    for (const auto &minus_flavor : flavors) {
      auto lts = ww_decay_chain(plus_flavor, minus_flavor);
      INFO("W+ charged lepton PDG " << -plus_flavor.charged_abs
                                    << ", W- charged lepton PDG "
                                    << minus_flavor.charged_abs);
      REQUIRE(process.MatchProcess(lts.decaytree).has_value());
      const auto result = amplitude.Evaluate(lts, 0.0, false);
      REQUIRE(result.Valid());
      REQUIRE(result.amp2 > 0.0);
    }
  }

  // Verify generated subprocesses and external states use the same model masses
  std::vector<gra::mg5::Subprocess<MG5_YY_WW::ProcessBase>> subprocesses;
  MG5_YY_WW::BuildSubprocesses(
      subprocesses,
      gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_WW/param_card.dat"));
  MG5_YY_WW::ProcessBase *tau_process = nullptr;
  const std::vector<int> tau_final = {-15, 16, 15, -16};
  for (auto &subprocess : subprocesses) {
    if (std::any_of(subprocess.channels.begin(), subprocess.channels.end(),
                    [&](const auto &channel) {
                      return channel.final == tau_final;
                    })) {
      tau_process = subprocess.process.get();
      break;
    }
  }
  REQUIRE(tau_process != nullptr);
  gra::mg5::SubprocessSum<MG5_YY_WW::ProcessBase> sum(
      std::move(subprocesses));
  auto tau_chain = ww_decay_chain(flavors[2], flavors[2]);
  tau_chain.id1 = gra::PDG::PDG_gamma;
  tau_chain.id2 = gra::PDG::PDG_gamma;
  tau_chain.alphaQCD = 0.0;
  REQUIRE(sum.Prepare(tau_chain, 0.0) == gra::mg5helas::EvaluationStatus::Success);
  double tau_amp2 = 0.0;
  REQUIRE(sum.CalcPreparedHelicityAmp2(tau_chain, tau_amp2) == gra::mg5helas::EvaluationStatus::Success);
  REQUIRE(tau_amp2 > 0.0);
  REQUIRE(tau_process->getMasses().size() == 6);
  REQUIRE(tau_process->getMasses()[2] == Approx(flavors[2].mass).margin(1e-10));
  REQUIRE(gra::math::IsZero(tau_process->getMasses()[3]));
  REQUIRE(tau_process->getMasses()[4] == Approx(flavors[2].mass).margin(1e-10));
  REQUIRE(gra::math::IsZero(tau_process->getMasses()[5]));

  const auto coherent_decay = process.DecayStructureFor(generated_chain);
  REQUIRE(coherent_decay.type == gra::DecayType::Full);
  REQUIRE((coherent_decay.type == gra::DecayType::Full));

  gra::MGeneratedPhotonProc dz_process("yy_DZ", "WW", description,
                                       "MG5_YY_WW");
  gra::MGeneratedPhotonProc lux_process("yy_LUX", "WW", description,
                                        "MG5_YY_WW");
  REQUIRE(dz_process.DecayStructureFor(generated_chain) == coherent_decay);
  REQUIRE(lux_process.DecayStructureFor(generated_chain) == coherent_decay);

  gra::LORENTZSCALAR direct_same_leaves = generated_chain;
  direct_same_leaves.decaytree.clear();
  for (const auto &branch : generated_chain.decaytree) {
    direct_same_leaves.decaytree.insert(direct_same_leaves.decaytree.end(),
                                        branch.legs.begin(), branch.legs.end());
  }
  REQUIRE(gra::PhotonStablePDGs(direct_same_leaves.decaytree) ==
          gra::PhotonStablePDGs(generated_chain.decaytree));
  REQUIRE_FALSE(process.MatchProcess(direct_same_leaves.decaytree).has_value());
  CHECK_THROWS_AS(process.DecayStructureFor(direct_same_leaves),
                  std::invalid_argument);

  auto mismatched_vertex = ww_decay_chain(flavors[1], flavors[0]);
  mismatched_vertex.decaytree[0].legs[1].p.pdg = 12;
  REQUIRE_FALSE(process.MatchProcess(mismatched_vertex.decaytree).has_value());

  auto isolated_lts = ww_decay_chain(flavors[0], flavors[0]);
  isolated_lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
  REQUIRE_THROWS(process.DecayStructureFor(isolated_lts));
  const auto isolated_result = amplitude.Evaluate(isolated_lts, 0.0, false);
  REQUIRE(isolated_result.status ==
          gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  REQUIRE(gra::math::IsZero(isolated_result.amp2));
}

TEST_CASE("Cascade decay syntax builds the physical Jacob-Wick vertex tree",
          "[gra::spin][cascade][syntax][helicity]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "CON",
                                "rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}");

  REQUIRE(process.state.lts.decaytree.size() == 2);
  for (const auto &branch : process.state.lts.decaytree) {
    REQUIRE(branch.p.pdg == 113);
    REQUIRE(branch.p.spinX2 == 2);
    REQUIRE(branch.p.P == -1);
    REQUIRE(branch.p.C == -1);
    REQUIRE(branch.legs.size() == 2);
    REQUIRE(branch.legs[0].p.pdg == 211);
    REQUIRE(branch.legs[1].p.pdg == -211);
    REQUIRE(branch.legs[0].p.spinX2 == 0);
    REQUIRE(branch.legs[1].p.spinX2 == 0);
    REQUIRE(branch.legs[0].legs.empty());
    REQUIRE(branch.legs[1].legs.empty());
  }

  REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
  for (const auto &branch : process.state.lts.decaytree) {
    const gra::HELMatrix &hel = branch.hel;
    REQUIRE(hel.J == Approx(1.0));
    REQUIRE(hel.s1 == Approx(0.0));
    REQUIRE(hel.s2 == Approx(0.0));
    REQUIRE(hel.Jz_values == std::vector<double>{-1.0, 0.0, 1.0});
    REQUIRE(hel.T.size_row() == 1);
    REQUIRE(hel.T.size_col() == 1);
    REQUIRE(std::norm(hel.T[0][0]) == Approx(1.0).margin(1e-14));
    REQUIRE(std::isfinite(std::real(hel.g_decay)));
    REQUIRE(std::isfinite(std::imag(hel.g_decay)));
    REQUIRE(std::abs(hel.g_decay) > 0.0);

    const auto decay = gra::spin::fDecayMatrix(hel, 0.83, -0.47);
    REQUIRE(decay.size_row() == 1);
    REQUIRE(decay.size_col() == 3);
    double helicity_probability = 0.0;
    for (std::size_t m = 0; m < decay.size_col(); ++m) {
      helicity_probability += std::norm(decay[0][m]);
    }
    REQUIRE(helicity_probability == Approx(1.0).margin(1e-13));

    const auto collinear = gra::spin::fDecayMatrix(hel, 0.0, 1.21);
    REQUIRE(std::abs(collinear[0][0]) == Approx(0.0).margin(1e-14));
    REQUIRE(std::abs(collinear[0][1]) == Approx(1.0).margin(1e-14));
    REQUIRE(std::abs(collinear[0][2]) == Approx(0.0).margin(1e-14));
  }

  SECTION("nested and asymmetric text preserves every ordered tree edge") {
    const std::string syntax = "rho(770)0>{pi0>{gamma gamma} pi0} "
                               "rho(770)0>{pi+ pi-}";
    ToyHelicityProcess nested;
    ConfigureToyProductionProcess(nested, "MP", "CON", syntax);
    REQUIRE(nested.state.lts.decaytree.size() == 2);

    const auto &left = nested.state.lts.decaytree[0];
    const auto &right = nested.state.lts.decaytree[1];
    REQUIRE(left.p.pdg == 113);
    REQUIRE(left.depth == 0);
    REQUIRE(left.legs.size() == 2);
    REQUIRE(left.legs[0].p.pdg == 111);
    REQUIRE(left.legs[0].depth == 1);
    REQUIRE(left.legs[0].legs.size() == 2);
    REQUIRE(left.legs[0].legs[0].p.pdg == 22);
    REQUIRE(left.legs[0].legs[1].p.pdg == 22);
    REQUIRE(left.legs[0].legs[0].depth == 2);
    REQUIRE(left.legs[0].legs[1].depth == 2);
    REQUIRE(left.legs[1].p.pdg == 111);
    REQUIRE(left.legs[1].depth == 1);
    REQUIRE(left.legs[1].legs.empty());

    REQUIRE(right.p.pdg == 113);
    REQUIRE(right.depth == 0);
    REQUIRE(right.legs.size() == 2);
    REQUIRE(right.legs[0].p.pdg == 211);
    REQUIRE(right.legs[1].p.pdg == -211);
    REQUIRE(right.legs[0].depth == 1);
    REQUIRE(right.legs[1].depth == 1);

    const auto check_names = [&](const auto &self,
                                 const gra::MDecayBranch &branch) -> void {
      REQUIRE(branch.name == std::to_string(branch.p.pdg) + "#" +
                                 std::to_string(branch.depth));
      for (const auto &daughter : branch.legs) {
        REQUIRE(daughter.depth == branch.depth + 1);
        self(self, daughter);
      }
    };
    check_names(check_names, left);
    check_names(check_names, right);
  }

  SECTION("whitespace does not alter the parsed daughter ordering") {
    ToyHelicityProcess compact;
    ToyHelicityProcess spaced;
    ConfigureToyProductionProcess(
        compact, "MP", "CON",
        "rho(770)0>{pi0>{gamma gamma} pi0}rho(770)0>{pi+ pi-}");
    ConfigureToyProductionProcess(
        spaced, "MP", "CON",
        " rho(770)0 > { pi0 > { gamma gamma } pi0 }   "
        "rho(770)0 > { pi+ pi- } ");
    const auto same_tree = [&](const auto &self, const gra::MDecayBranch &a,
                               const gra::MDecayBranch &b) -> void {
      REQUIRE(a.p.pdg == b.p.pdg);
      REQUIRE(a.depth == b.depth);
      REQUIRE(a.legs.size() == b.legs.size());
      for (std::size_t i = 0; i < a.legs.size(); ++i) {
        self(self, a.legs[i], b.legs[i]);
      }
    };
    REQUIRE(compact.state.lts.decaytree.size() ==
            spaced.state.lts.decaytree.size());
    for (std::size_t i = 0; i < compact.state.lts.decaytree.size(); ++i) {
      same_tree(same_tree, compact.state.lts.decaytree[i],
                spaced.state.lts.decaytree[i]);
    }
  }

  SECTION("unbalanced braces are rejected before particle lookup") {
    ToyHelicityProcess malformed;
    malformed.state.lts.PDG = LoadedPDGTable();
    REQUIRE_THROWS(malformed.SetDecayMode("rho(770)0 > {pi+ pi- rho(770)0"));
  }

  SECTION("every cascade arrow requires one daughter block") {
    ToyHelicityProcess malformed;
    malformed.state.lts.PDG = LoadedPDGTable();
    REQUIRE_THROWS(malformed.SetDecayMode("rho(770)0 > pi+ pi-"));
  }

  SECTION("empty and misplaced daughter blocks are rejected") {
    for (const std::string syntax :
         {"rho(770)0 > {} rho(770)0", "rho(770)0 > } pi+ pi- {",
          "rho(770)0 > {pi+ pi-} > {gamma gamma}",
          "rho(770)0 > {pi0 > {} pi0} rho(770)0"}) {
      CAPTURE(syntax);
      ToyHelicityProcess malformed;
      malformed.state.lts.PDG = LoadedPDGTable();
      REQUIRE_THROWS_AS(malformed.SetDecayMode(syntax), std::invalid_argument);
      REQUIRE(malformed.state.lts.decaytree.empty());
    }
  }

  SECTION("a parse failure leaves an existing tree unchanged") {
    std::vector<gra::MDecayBranch> tree;
    LoadedPDGTable().TokenizeProcess("pi+ pi-", 0, tree);
    REQUIRE(tree.size() == 2);
    const auto before = tree;
    REQUIRE_THROWS_AS(
        LoadedPDGTable().TokenizeProcess("rho(770)0>{pi0>{} pi0}", 0, tree),
        std::invalid_argument);
    REQUIRE(tree.size() == before.size());
    for (std::size_t i = 0; i < tree.size(); ++i) {
      REQUIRE(tree[i].p.pdg == before[i].p.pdg);
      REQUIRE(tree[i].depth == before[i].depth);
      REQUIRE(tree[i].legs.size() == before[i].legs.size());
    }
  }
}

TEST_CASE("Cascade DECAY_SYM coherently sums Bose assignments",
          "[gra::spin][cascade][symmetry][helicity]") {
  gra::LORENTZSCALAR unsymmetrized = ContinuumVectorCascadeLTSForTest(false);
  unsymmetrized.amplitude.DECAY_SYM = false;
  unsymmetrized.decay_symmetry_proposal_active = false;
  const auto amplitude_off =
      gra::spin::ContinuumDecayMatrix(unsymmetrized, "CM");

  gra::LORENTZSCALAR symmetrized = ContinuumVectorCascadeLTSForTest(false);
  symmetrized.amplitude.DECAY_SYM = true;
  symmetrized.decay_symmetry_proposal_active = false;
  const auto terms =
      gra::spin::StableLeafAmplitudeTrees(symmetrized, symmetrized.decaytree);
  REQUIRE(terms.size() == 2);
  for (const auto &term : terms) {
    REQUIRE(gra::math::IsExactEqual(term.statistics_sign, 1.0));
  }

  const auto amplitude_on = gra::spin::ContinuumDecayMatrix(symmetrized, "CM");
  MMatrix<std::complex<double>> coherent_reference;
  bool initialized = false;
  for (const auto &term : terms) {
    gra::LORENTZSCALAR assigned = symmetrized;
    assigned.amplitude.DECAY_SYM = false;
    assigned.decaytree = term.tree;
    const auto assignment_amplitude =
        gra::spin::ContinuumDecayMatrix(assigned, "CM");
    const auto contribution = assignment_amplitude * term.statistics_sign;
    if (!initialized) {
      coherent_reference = contribution;
      initialized = true;
    } else {
      coherent_reference = coherent_reference + contribution;
    }
  }
  RequireMatrixNear(amplitude_on, coherent_reference, 1e-11);
  REQUIRE(MatrixDiffNorm2(amplitude_on, amplitude_off) > 1e-10);

  gra::LORENTZSCALAR crossed_reference = symmetrized;
  crossed_reference.decaytree = terms[1].tree;
  const auto crossed_amplitude =
      gra::spin::ContinuumDecayMatrix(crossed_reference, "CM");
  RequireMatrixNear(crossed_amplitude, amplitude_on, 1e-11);
}

TEST_CASE("Stable-leaf assignments preserve Bose and Fermi statistics exactly",
          "[gra::spin][cascade][symmetry]") {
  const auto build_tree = [](int identical_spin_x2) {
    gra::MDecayBranch first;
    first.p = ToyParticle("x", 900301, identical_spin_x2, 0.10);
    first.p4 = gra::M4Vec(0.11, -0.03, 0.17, 0.31);
    gra::MDecayBranch second = first;
    second.p4 = gra::M4Vec(-0.07, 0.09, -0.13, 0.29);

    gra::MDecayBranch spectator_left;
    spectator_left.p = ToyParticle("a", 900302, 0, 0.12);
    spectator_left.p4 = gra::M4Vec(0.02, 0.06, 0.08, 0.22);
    gra::MDecayBranch spectator_right;
    spectator_right.p = ToyParticle("b", 900303, 0, 0.14);
    spectator_right.p4 = gra::M4Vec(-0.04, -0.05, 0.06, 0.24);

    gra::MDecayBranch left;
    left.p = ToyParticle("L", 900310, 0, 0.60);
    left.legs = {first, spectator_left};
    left.p4 = first.p4 + spectator_left.p4;
    gra::MDecayBranch right;
    right.p = ToyParticle("R", 900311, 0, 0.65);
    right.legs = {second, spectator_right};
    right.p4 = second.p4 + spectator_right.p4;
    return std::vector<gra::MDecayBranch>{left, right};
  };

  for (const int spin_x2 : {0, 1}) {
    CAPTURE(spin_x2);
    gra::LORENTZSCALAR lts;
    const auto reference = build_tree(spin_x2);
    gra::spin::PrepareStableLeafSymmetryAssignments(lts, reference);
    REQUIRE(lts.decay_symmetry_assignments.size() == 2);
    REQUIRE(lts.decay_symmetry_statistics_signs.size() == 2);
    REQUIRE(lts.decay_symmetry_assignments[0] ==
            std::vector<std::size_t>{0, 1, 2, 3});
    REQUIRE(lts.decay_symmetry_assignments[1] ==
            std::vector<std::size_t>{2, 1, 0, 3});
    REQUIRE(
        gra::math::IsExactEqual(lts.decay_symmetry_statistics_signs[0], 1.0));
    REQUIRE(gra::math::IsExactEqual(lts.decay_symmetry_statistics_signs[1],
                                    spin_x2 == 0 ? 1.0 : -1.0));

    const auto cached_assignments = lts.decay_symmetry_assignments;
    const auto cached_signs = lts.decay_symmetry_statistics_signs;
    const std::string cached_key = lts.decay_symmetry_topology_key;
    gra::spin::PrepareStableLeafSymmetryAssignments(lts, reference);
    REQUIRE(lts.decay_symmetry_assignments == cached_assignments);
    REQUIRE(lts.decay_symmetry_statistics_signs == cached_signs);
    REQUIRE(lts.decay_symmetry_topology_key == cached_key);

    auto crossed = reference;
    REQUIRE(gra::spin::ApplyStableLeafSymmetryAssignment(
        crossed, lts.decay_symmetry_assignments[1]));
    REQUIRE(crossed[0].legs[0].p4 == reference[1].legs[0].p4);
    REQUIRE(crossed[1].legs[0].p4 == reference[0].legs[0].p4);
    REQUIRE(crossed[0].p4 == crossed[0].legs[0].p4 + crossed[0].legs[1].p4);
    REQUIRE(crossed[1].p4 == crossed[1].legs[0].p4 + crossed[1].legs[1].p4);

    const auto unchanged = crossed;
    REQUIRE_FALSE(gra::spin::ApplyStableLeafSymmetryAssignment(
        crossed, std::vector<std::size_t>{0, 1, 1, 3}));
    REQUIRE(crossed[0].p4 == unchanged[0].p4);
    REQUIRE(crossed[1].p4 == unchanged[1].p4);
    REQUIRE_FALSE(gra::spin::ApplyStableLeafSymmetryAssignment(
        crossed, std::vector<std::size_t>{1, 0, 2, 3}));
    REQUIRE(crossed[0].p4 == unchanged[0].p4);
    REQUIRE(crossed[1].p4 == unchanged[1].p4);

    auto inconsistent = reference;
    inconsistent[1].legs[0].p.spinX2 = spin_x2 == 0 ? 1 : 0;
    gra::LORENTZSCALAR invalid_lts;
    REQUIRE_THROWS_AS(gra::spin::PrepareStableLeafSymmetryAssignments(
                          invalid_lts, inconsistent),
                      std::invalid_argument);
  }

  SECTION(
      "three-particle signs give the exact symmetric and alternating sums") {
    const auto build_three = [](int spin_x2) {
      std::vector<gra::MDecayBranch> tree;
      for (std::size_t i = 0; i < 3; ++i) {
        gra::MDecayBranch identical;
        identical.p = ToyParticle("x", 900320, spin_x2, 0.10);
        identical.p4 = gra::M4Vec(0.04 + 0.01 * i, -0.02 + 0.005 * i,
                                  0.03 - 0.007 * i, 0.20 + 0.02 * i);
        gra::MDecayBranch spectator;
        spectator.p =
            ToyParticle("s", 900330 + static_cast<int>(i), 0, 0.12 + 0.01 * i);
        spectator.p4 = gra::M4Vec(0.01 * i, 0.02, -0.01, 0.22 + 0.01 * i);
        gra::MDecayBranch parent;
        parent.p = ToyParticle("P", 900340 + static_cast<int>(i), 0, 0.7);
        parent.legs = {identical, spectator};
        parent.p4 = identical.p4 + spectator.p4;
        tree.push_back(parent);
      }
      return tree;
    };

    for (const int spin_x2 : {0, 1}) {
      CAPTURE(spin_x2);
      gra::LORENTZSCALAR lts;
      const auto reference = build_three(spin_x2);
      gra::spin::PrepareStableLeafSymmetryAssignments(lts, reference);
      REQUIRE(lts.decay_symmetry_assignments.size() == 6);
      double coherent_equal_kernel_sum = 0.0;
      for (std::size_t term = 0; term < lts.decay_symmetry_assignments.size();
           ++term) {
        const auto &assignment = lts.decay_symmetry_assignments[term];
        const std::array<std::size_t, 3> sources = {
            assignment[0] / 2, assignment[2] / 2, assignment[4] / 2};
        std::size_t inversions = 0;
        for (std::size_t i = 0; i < sources.size(); ++i) {
          for (std::size_t j = i + 1; j < sources.size(); ++j) {
            inversions += sources[i] > sources[j] ? 1 : 0;
          }
        }
        const double expected_sign =
            spin_x2 == 0 || inversions % 2 == 0 ? 1.0 : -1.0;
        REQUIRE(gra::math::IsExactEqual(
            lts.decay_symmetry_statistics_signs[term], expected_sign));
        coherent_equal_kernel_sum += expected_sign;

        auto assigned = reference;
        REQUIRE(
            gra::spin::ApplyStableLeafSymmetryAssignment(assigned, assignment));
        for (std::size_t branch = 0; branch < assigned.size(); ++branch) {
          REQUIRE(assigned[branch].legs[0].p4 ==
                  reference[sources[branch]].legs[0].p4);
          REQUIRE(assigned[branch].p4 ==
                  assigned[branch].legs[0].p4 + assigned[branch].legs[1].p4);
        }
      }
      REQUIRE(gra::math::IsExactEqual(coherent_equal_kernel_sum,
                                      spin_x2 == 0 ? 6.0 : 0.0));
    }
  }
}

TEST_CASE("gra::spin InitTMatrix enforces P and C selection rules",
          "[gra::spin]") {
  SECTION("half-integer helicity states have the physical projections") {
    const auto spin_half = gra::spin::SpinProjections(0.5);
    REQUIRE(spin_half.size() == 2);
    CHECK(spin_half[0] == Approx(-0.5));
    CHECK(spin_half[1] == Approx(0.5));
    const auto D = gra::wigner::DMatrix(0.5, 0.1, 0.2);
    CHECK(D.size_row() == 2);
    CHECK(D.size_col() == 2);
    CHECK_THROWS_AS(gra::spin::SpinProjections(0.25), std::invalid_argument);
  }

  SECTION("parity can allow or reject the same LS entry") {
    const auto mother_even = ToyParticle(9000001, 0, 1, 0, "0++");
    const auto mother_odd = ToyParticle(9000002, 0, -1, 0, "0-+");
    const auto pi_plus = ToyParticle(211, 0, -1, 0, "pi+");
    const auto pi_minus = ToyParticle(-211, 0, -1, 0, "pi-");

    auto allowed = ToyHelicityMatrix(0, 0, true, false);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(allowed, mother_even, pi_plus,
                                           pi_minus, false, "test P allowed",
                                           false, false));

    auto forbidden = ToyHelicityMatrix(0, 0, true, false);
    REQUIRE_THROWS_AS(gra::spin::InitTMatrix(forbidden, mother_odd, pi_plus,
                                             pi_minus, false,
                                             "test P forbidden", false, false),
                      std::invalid_argument);
  }

  SECTION("particle-antiparticle C uses (-1)^(L+S)") {
    const auto scalar_c_plus = ToyParticle(9000003, 0, 1, 1, "0++");
    const auto scalar_c_minus = ToyParticle(9000004, 0, 1, -1, "0+-");
    const auto vector_c_minus = ToyParticle(9000005, 2, -1, -1, "1--");
    const auto pi_plus = ToyParticle(211, 0, -1, 0, "pi+");
    const auto pi_minus = ToyParticle(-211, 0, -1, 0, "pi-");

    auto c_even = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(c_even, scalar_c_plus, pi_plus,
                                           pi_minus, false, "test C plus",
                                           false, false));

    auto c_even_forbidden = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_THROWS_AS(gra::spin::InitTMatrix(
                          c_even_forbidden, scalar_c_minus, pi_plus, pi_minus,
                          false, "test C minus forbidden", false, false),
                      std::invalid_argument);

    auto c_odd = ToyHelicityMatrix(1, 0, true, true);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(c_odd, vector_c_minus, pi_plus,
                                           pi_minus, false, "test C minus",
                                           false, false));
  }

  SECTION("C-eigen daughter pairs use C1*C2") {
    const auto scalar_c_plus = ToyParticle(9000006, 0, 1, 1, "0++");
    const auto scalar_c_minus = ToyParticle(9000007, 0, 1, -1, "0+-");
    const auto rho0 = ToyParticle(113, 2, -1, -1, "rho0");

    auto allowed = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(allowed, scalar_c_plus, rho0, rho0,
                                           false, "test C eigen", false,
                                           false));

    auto forbidden = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_THROWS_AS(
        gra::spin::InitTMatrix(forbidden, scalar_c_minus, rho0, rho0, false,
                               "test C eigen forbidden", false, false),
        std::invalid_argument);
  }

  SECTION("C_symmetry=false bypasses C filtering") {
    const auto scalar_c_minus = ToyParticle(9000008, 0, 1, -1, "0+-");
    const auto pi_plus = ToyParticle(211, 0, -1, 0, "pi+");
    const auto pi_minus = ToyParticle(-211, 0, -1, 0, "pi-");

    auto bypass = ToyHelicityMatrix(0, 0, true, false);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(bypass, scalar_c_minus, pi_plus,
                                           pi_minus, false, "test C bypass",
                                           false, false));
  }

  SECTION("allowed LS and direct-basis discovery share one row classifier") {
    const auto scalar_c_plus = ToyParticle(9000009, 0, 1, 1, "0++");
    const auto rho0 = ToyParticle(113, 2, -1, -1, "rho0");

    gra::HELMatrix prototype;
    prototype.P_symmetry = true;
    prototype.C_symmetry = true;
    prototype.alpha_ls.Clear();

    const auto allowed = gra::spin::AllowedLSCouplings(
        prototype, scalar_c_plus, rho0, rho0, false, "typed LS rejection",
        gra::spin::VertexContext::Auto);
    const auto raw_basis = gra::spin::BuildJWDirectHelicityRawBasis(
        prototype, scalar_c_plus, rho0, rho0, false, "typed LS rejection", "",
        gra::spin::VertexContext::Auto);

    REQUIRE_FALSE(allowed.empty());
    REQUIRE(raw_basis.size() == allowed.size());
  }

  SECTION("missing allowed rows are accepted while forbidden rows throw") {
    const auto scalar_c_plus = ToyParticle(9000009, 0, 1, 1, "0++");
    const auto vector_c_minus = ToyParticle(113, 2, -1, -1, "rho0");

    auto missing = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(missing, scalar_c_plus,
                                           vector_c_minus, vector_c_minus,
                                           false, "DECAYS.json", true, false));

    auto extra = ToyHelicityMatrix(1, 2, true, true);
    REQUIRE_THROWS(gra::spin::InitTMatrix(extra, scalar_c_plus, vector_c_minus,
                                          vector_c_minus, false, "DECAYS.json",
                                          true, false));
  }

  SECTION("direct helicity amplitudes are checked against the JW subspace") {
    const auto scalar_c_plus = ToyParticle(9000010, 0, 1, 1, "0++");
    const auto vector_c_minus = ToyParticle(113, 2, -1, -1, "rho0");

    auto reference = ToyHelicityMatrix(0, 0, true, true);
    REQUIRE_NOTHROW(gra::spin::InitTMatrix(reference, scalar_c_plus,
                                           vector_c_minus, vector_c_minus,
                                           false, "DECAYS.json", true, false));

    gra::HELMatrix direct;
    direct.BR = 1.0;
    direct.P_symmetry = true;
    direct.C_symmetry = true;
    direct.coupling_basis = gra::CouplingBasis::Helicity;
    direct.T = reference.T;
    direct.T_set =
        MMatrix<bool>(reference.T.size_row(), reference.T.size_col(), false);
    for (std::size_t i = 0; i < reference.T.size_row(); ++i) {
      for (std::size_t j = 0; j < reference.T.size_col(); ++j) {
        if (std::abs(reference.T[i][j]) > 1e-12) {
          direct.T_set[i][j] = true;
        }
      }
    }
    REQUIRE_NOTHROW(gra::spin::ValidateDirectTMatrix(
        direct, scalar_c_plus, vector_c_minus, vector_c_minus, false,
        "DECAYS.json", true, false));
    for (std::size_t i = 0; i < reference.T.size_row(); ++i) {
      for (std::size_t j = 0; j < reference.T.size_col(); ++j) {
        CHECK(std::abs(direct.T[i][j] - reference.T[i][j]) < 1e-12);
      }
    }

    auto forbidden = direct;
    forbidden.T_set[0][1] = true;
    REQUIRE_THROWS(gra::spin::ValidateDirectTMatrix(
        forbidden, scalar_c_plus, vector_c_minus, vector_c_minus, false,
        "DECAYS.json", true, false));

    auto symmetry_breaking = direct;
    symmetry_breaking.T = MMatrix<std::complex<double>>(3, 3, 0.0);
    symmetry_breaking.T[0][0] = 1.0;
    symmetry_breaking.T_set = MMatrix<bool>(3, 3, false);
    symmetry_breaking.T_set[0][0] = true;
    symmetry_breaking.T_set[1][1] = true;
    symmetry_breaking.T_set[2][2] = true;
    REQUIRE_THROWS(gra::spin::ValidateDirectTMatrix(
        symmetry_breaking, scalar_c_plus, vector_c_minus, vector_c_minus, false,
        "DECAYS.json", true, false));
  }

  SECTION("complex LS superposition equals the Jacob-Wick sum entrywise") {
    const auto mother = ToyParticle(9000020, 2, 1, 0, "J1");
    const auto first = ToyParticle(9000021, 2, -1, 0, "V1");
    const auto second = ToyParticle(9000022, 2, -1, 0, "V2");

    gra::HELMatrix hel;
    hel.BR = 1.0;
    hel.P_symmetry = false;
    hel.C_symmetry = false;
    hel.alpha_ls.Set(1, 0, {0.30, 0.40});
    hel.alpha_ls.Set(0, 2, {-0.20, 0.50});
    hel.alpha_ls.Set(2, 2, {0.60, -0.10});

    const std::array<gra::spin::LSCoupling, 3> rows = {
        gra::spin::LSCoupling{1, 0}, gra::spin::LSCoupling{0, 2},
        gra::spin::LSCoupling{2, 2}};
    const double alpha_norm2 = hel.alpha_ls.Norm2();

    MMatrix<std::complex<double>> expected(3, 3, 0.0);
    constexpr double J = 1.0;
    constexpr double s1 = 1.0;
    constexpr double s2 = 1.0;
    for (const auto &row : rows) {
      const double l = static_cast<double>(row.l);
      const double s = 0.5 * static_cast<double>(row.two_s);
      const std::complex<double> alpha =
          hel.alpha_ls.At(row.l, row.two_s) / std::sqrt(alpha_norm2);
      const double normalization = std::sqrt((2.0 * l + 1.0) / (2.0 * J + 1.0));
      for (std::size_t i = 0; i < 3; ++i) {
        const double lambda1 = -s1 + static_cast<double>(i);
        for (std::size_t j = 0; j < 3; ++j) {
          const double lambda2 = -s2 + static_cast<double>(j);
          const double lambda = lambda1 - lambda2;
          expected[i][j] +=
              alpha * normalization *
              gra::wigner::CG(l, s, 0.0, lambda, J, lambda) *
              gra::wigner::CG(s1, s2, lambda1, -lambda2, s, lambda);
        }
      }
    }

    gra::spin::InitTMatrix(hel, mother, first, second, false, "complex LS sum",
                           false, false);
    RequireMatrixNear(hel.T, expected, 1e-13);
    REQUIRE(hel.T.FrobNorm2() == Approx(1.0).margin(1e-13));
    for (const auto &row : rows) {
      REQUIRE(std::abs(hel.alpha_ls.At(row.l, row.two_s) -
                       (row.l == 1 && row.two_s == 0
                            ? std::complex<double>(0.30, 0.40)
                        : row.l == 0 ? std::complex<double>(-0.20, 0.50)
                                     : std::complex<double>(0.60, -0.10)) /
                           std::sqrt(alpha_norm2)) < 1e-14);
    }
  }

  SECTION("canonical minus-pi phase preserves the Jacob-Wick LS sign") {
    const auto mother = ToyParticle(9000030, 2, 1, 0, "J1");
    const auto first = ToyParticle(9000031, 2, -1, 0, "V1");
    const auto second = ToyParticle(9000032, 2, -1, 0, "V2");

    gra::HELMatrix phased;
    phased.BR = 1.0;
    phased.P_symmetry = false;
    phased.C_symmetry = false;
    phased.alpha_ls.Set(1, 0, 1.0);
    phased.alpha_ls.Set(0, 2, std::polar(0.4, -gra::math::PI));

    gra::HELMatrix signed_real = phased;
    signed_real.alpha_ls.Set(0, 2, {-0.4, 0.0});

    gra::spin::InitTMatrix(phased, mother, first, second, false,
                           "canonical phased LS sum", false, false);
    gra::spin::InitTMatrix(signed_real, mother, first, second, false,
                           "signed real LS sum", false, false);

    RequireMatrixNear(phased.T, signed_real.T, 1.0e-14);
    CHECK(phased.T.FrobNorm2() ==
          Approx(signed_real.T.FrobNorm2()).margin(1.0e-14));
  }
}

TEST_CASE("sparse LS coefficients stay unique and ordered", "[gra::spin][LS]") {
  gra::spin::LSCoefficients alpha;

  REQUIRE(alpha.Empty());
  REQUIRE(alpha.Insert(2, 4, {0.2, -0.1}));
  REQUIRE(alpha.Insert(0, 2, {0.5, 0.0}));
  REQUIRE(alpha.Insert(2, 0, {-0.3, 0.4}));
  REQUIRE_FALSE(alpha.Insert(2, 0, {9.0, 9.0}));
  REQUIRE(alpha.Size() == 3);

  auto term = alpha.cbegin();
  CHECK(term->l == 0);
  CHECK(term->two_s == 2);
  ++term;
  CHECK(term->l == 2);
  CHECK(term->two_s == 0);
  ++term;
  CHECK(term->l == 2);
  CHECK(term->two_s == 4);

  alpha.Set(2, 0, {0.6, 0.8});
  CHECK(std::abs(alpha.At(2, 0) - std::complex<double>(0.6, 0.8)) < 1e-14);
  CHECK(alpha.Norm2() == Approx(1.3));
  CHECK(alpha.IsFinite());
  CHECK_THROWS_AS(alpha.At(1, 1), std::out_of_range);

  alpha.Set(6, 2, 1.0e-7);
  alpha.RemoveBelow(1.0e-6);
  CHECK_FALSE(alpha.Contains(6, 2));
}

TEST_CASE("gra::spin excludes longitudinal real-vector final helicities",
          "[gra::spin][photon][physics]") {
  const auto scalar = ToyParticle(9000011, 0, 1, 1, "0++");
  const auto photon = ToyParticle(22, 2, -1, -1, "gamma");
  auto decay = ToyHelicityMatrix(0, 0, true, true);

  REQUIRE_NOTHROW(gra::spin::InitTMatrix(decay, scalar, photon, photon, false,
                                         "gamma gamma test", false, false));
  REQUIRE(gra::spin::FinalStateHelicityCount(photon, "gamma test") == 2);
  REQUIRE(decay.lambda_values.size_row() == 4);
  REQUIRE(decay.T.size_row() == 3);
  REQUIRE(decay.T.size_col() == 3);
  for (std::size_t row = 0; row < decay.lambda_values.size_row(); ++row) {
    CHECK(std::abs(decay.lambda_values[row][0]) == Approx(1.0));
    CHECK(std::abs(decay.lambda_values[row][1]) == Approx(1.0));
    CHECK(decay.lambda_idx[row][0] != 1);
    CHECK(decay.lambda_idx[row][1] != 1);
  }
  for (std::size_t i = 0; i < 3; ++i) {
    CHECK(std::abs(decay.T[1][i]) < 1e-14);
    CHECK(std::abs(decay.T[i][1]) < 1e-14);
  }
  CHECK(decay.T.FrobNorm2() == Approx(1.0).margin(1e-12));

  auto direct = decay;
  direct.coupling_basis = gra::CouplingBasis::Helicity;
  direct.T_set = MMatrix<bool>(3, 3, false);
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      if (std::abs(direct.T[i][j]) > 1e-12) {
        direct.T_set[i][j] = true;
      }
    }
  }
  direct.T[1][1] = 0.1;
  direct.T_set[1][1] = true;
  REQUIRE_THROWS(gra::spin::ValidateDirectTMatrix(
      direct, scalar, photon, photon, false, "gamma gamma direct test", false,
      false));

  const auto spin_one = ToyParticle(9000015, 2, -1, 1, "1-+");
  auto landau_yang = ToyHelicityMatrix(1, 2, true, true);
  REQUIRE_THROWS(gra::spin::InitTMatrix(landau_yang, spin_one, photon, photon,
                                        false, "Landau-Yang test", false,
                                        false));
}

TEST_CASE("gra::spin spin-blind decay keeps parent spin states incoherent",
          "[gra::spin][decay][physics]") {
  gra::LORENTZSCALAR lts;
  lts.process.SPINDEC = false;
  gra::MDecayBranch left;
  left.p = ToyParticle("left", 9000012, 0, 0.1);
  gra::MDecayBranch right;
  right.p = ToyParticle("right", 9000013, 0, 0.1);
  lts.decaytree = {left, right};

  gra::PARAM_RES resonance;
  resonance.production_model = gra::ReggeProductionModel::XP;
  resonance.p = ToyParticle("spin2", 9000014, 4, 1.0);
  gra::spin::DecayAmp(lts, resonance, "CM");

  REQUIRE(resonance.decay_f.size_row() == 5);
  REQUIRE(resonance.decay_f.size_col() == 5);
  const auto reduced = resonance.decay_f * resonance.decay_f.Dagger();
  for (std::size_t i = 0; i < 5; ++i) {
    for (std::size_t j = 0; j < 5; ++j) {
      const double expected = (i == j) ? 0.2 : 0.0;
      CHECK(std::real(reduced[i][j]) == Approx(expected).margin(1e-14));
      CHECK(std::imag(reduced[i][j]) == Approx(0.0).margin(1e-14));
    }
  }
  CHECK(std::real(reduced.Trace()) == Approx(1.0).margin(1e-14));

  // Isolated mode must ignore physical decay spin settings
  lts.process.SPINDEC = true;
  lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
  resonance.decay_f = MMatrix<std::complex<double>>();
  gra::spin::DecayAmp(lts, resonance, "CM");
  RequireMatrixNear(resonance.decay_f * resonance.decay_f.Dagger(), reduced,
                    1e-14);
}

TEST_CASE("gra::spin positivity requires the exact integer or half-integer "
          "spin dimension",
          "[gra::spin][density][physics]") {
  MMatrix<std::complex<double>> spin_half(2, 2, 0.0);
  spin_half[0][0] = 0.5;
  spin_half[1][1] = 0.5;
  REQUIRE(gra::spin::Positivity(spin_half, 0.5));
  REQUIRE_THROWS(gra::spin::Positivity(
      spin_half, std::numeric_limits<double>::quiet_NaN()));

  MMatrix<std::complex<double>> oversized(4, 4, 0.0);
  for (std::size_t i = 0; i < 4; ++i) {
    oversized[i][i] = 0.25;
  }
  REQUIRE_THROWS(gra::spin::Positivity(oversized, 1.0));

  MMatrix<std::complex<double>> rectangular(3, 4, 0.0);
  REQUIRE_THROWS(gra::spin::Positivity(rectangular, 1.0));

  MMatrix<std::complex<double>> unnormalized(2, 2, 0.0);
  unnormalized[0][0] = 0.5;
  unnormalized[1][1] = 0.495;
  REQUIRE_THROWS(gra::spin::Positivity(unnormalized, 0.5));
}

// Check stable public failure behavior for invalid density eigensystems
TEST_CASE("gra::spin density eigensystem failures are stable",
          "[gra::spin][density]") {
  using Complex = std::complex<double>;
  const MMatrix<Complex> nonhermitian{{Complex(0.5, 0.0), Complex(0.0, 0.2)},
                                      {Complex(0.0, 0.2), Complex(0.5, 0.0)}};

  REQUIRE(gra::spin::VonNeumannEntropy(nonhermitian) == Approx(-1.0));
  REQUIRE_THROWS(gra::spin::SpectralPureStates(nonhermitian));
}

TEST_CASE(
    "rho cascades retain assignment normalization when spin correlations or coherence are disabled",
    "[MProcess][symmetry][physics]") {
  ToyHelicityProcess proc;
  const auto tree = TensorRhoCascadeLTSForTest().decaytree;

  for (const std::string family : {"MP", "XP", "GP"}) {
    for (const std::string channel : {"RES", "CON", "RES+CON"}) {
      for (const bool spindec : {false, true}) {
        proc.state.lts.process.SPINDEC = spindec;
        for (const bool coherent : {false, true}) {
          CAPTURE(family, channel, spindec, coherent);
          CHECK(proc.DecaySymmetryCompensationFactorForTest(tree, coherent, family, channel) ==
                Approx(spindec && coherent ? 1.0 : 2.0));
          CHECK(proc.DecaySymmetryCompensationFactorForTest(tree, coherent, family, channel, true) == Approx(1.0));
        }
      }
    }
  }
  proc.SetFLATAMP(4);
  for (const bool coherent : {false, true}) {
    CHECK(proc.DecaySymmetryCompensationFactorForTest(tree, coherent) == Approx(1.0));
    CHECK_FALSE(proc.StableLeafProposalActiveForTest(tree, coherent, "MP", "RES"));
  }
  proc.SetFLATAMP(0);
  CHECK(proc.DecaySymmetryCompensationFactorForTest(tree, false, "TP", "RES") ==
        Approx(1.0));
  CHECK(proc.DecaySymmetryCompensationFactorForTest(tree, false, "TP", "CON") ==
        Approx(1.0));
  CHECK_THROWS_AS(
      proc.DecaySymmetryCompensationFactorForTest(tree, false, "TP", "RES+CON"),
      gra::AmplitudeFailure);
}

// Check the incoherent multi-body fallback uses its declaration for both sampling and normalization
TEST_CASE("multi-body resonance cascades preserve incoherent assignment normalization",
          "[MProcess][symmetry][cascade][physics]") {
  ToyHelicityProcess proc;
  auto physical = TensorRhoCascadeLTSForTest();
  gra::MDecayBranch spectator;
  spectator.p = ToyParticle("pi0", 111, 0, 0.135);
  spectator.p4 = gra::M4Vec(0.0, 0.0, 0.0, spectator.p.mass);
  physical.decaytree.push_back(spectator);
  physical.pfinal[0] += spectator.p4;
  gra::PARAM_RES res;
  res.p = ToyParticle("X", 900310, 0, physical.pfinal[0].M());
  physical.process.SPINDEC = false;
  const auto reference = gra::spin::ResonanceDecayMatrix(physical, res, "CM");
  REQUIRE(reference.FrobNorm2() > 0.0);

  for (const std::string family : {"MP", "XP", "GP"}) {
    for (const bool spindec : {false, true}) {
      physical.process.SPINDEC = spindec;
      proc.state.lts.process.SPINDEC = spindec;
      for (const bool coherent : {false, true}) {
        CAPTURE(family, spindec, coherent);
        physical.amplitude.DECAY_SYM = coherent;
        RequireMatrixNear(gra::spin::ResonanceDecayMatrix(physical, res, "CM"), reference, 1e-11);
        CHECK(proc.DecaySymmetryCompensationFactorForTest(physical.decaytree, coherent, family, "RES") ==
              Approx(2.0));
        CHECK(proc.state.lts.decay_structure.type == gra::DecayType::JacobWickIncoherent);
        CHECK_FALSE(proc.StableLeafProposalActiveForTest(physical.decaytree, coherent, family, "RES", spindec));
      }
    }
  }
}

TEST_CASE("gra::spin production helicity rows use strict C-aware matching",
          "[gra::spin]") {
  ToyHelicityProcess proc;

  const auto pomeron0 = ToyParticle(991, 0, 1, 1, "P0");
  const auto photon = ToyParticle(22, 2, -1, -1, "gamma");
  const auto proton = ToyParticle(2212, 1, 1, 0, "p");
  const auto pbar = ToyParticle(-2212, 1, -1, 0, "pbar");
  const auto f0 = ToyParticle(9000221, 0, 1, 1, "f0");
  const auto rho0 = ToyParticle(113, 2, -1, -1, "rho0");
  const auto phi0 = ToyParticle(333, 2, -1, -1, "phi0");

  SECTION("legacy C-even central rows are absent from the production table") {
    REQUIRE_THROWS(proc.ProcessHelicityStructure(f0, {pomeron0, pomeron0}, true,
                                                 true, "", false));
  }

  SECTION("legacy gamma-Pomeron central rows are absent from the production "
          "table") {
    REQUIRE_THROWS(proc.ProcessHelicityStructure(rho0, {photon, pomeron0}, true,
                                                 true, "", false));
  }

  SECTION("charge-conjugate reversed particle-antiparticle rows are accepted") {
    REQUIRE_NOTHROW(proc.ProcessHelicityStructure(pomeron0, {pbar, proton},
                                                  true, true, "", false));
  }

  SECTION("legacy reversed gamma-Pomeron central rows are also absent") {
    REQUIRE_THROWS(proc.ProcessHelicityStructure(rho0, {pomeron0, photon}, true,
                                                 true, "", false));
  }

  SECTION("legacy phi gamma-Pomeron central rows are absent") {
    REQUIRE_THROWS(proc.ProcessHelicityStructure(phi0, {pomeron0, photon}, true,
                                                 true, "", false));
  }
}

TEST_CASE("Branching setup rejects unsupported direct-continuum and Tensor "
          "topologies",
          "[gra::MProcess][topology]") {
  SECTION("Tensor axial resonance accepts its direct three-body decay") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "TP", "RES", "K+ K- pi0");
    MRandom rng;
    auto f1 = gra::resonance::Read("RES/f1_1420.json", rng, gra::ReggeProductionModel::TP);
    SetToyTensorChannel(f1, {1.0, 0.0});
    proc.SetResonances({{"f1_1420", f1}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto configured = proc.GetResonances().at("f1_1420");
    REQUIRE(configured.TP.channels.size() == 1);
    const auto &channel = ToyTensorChannel(configured);
    REQUIRE(channel.g_tensor.size() == 2);
    CHECK(std::all_of(
        channel.g_tensor.begin(), channel.g_tensor.end(),
        [](const double coupling) { return std::isfinite(coupling); }));
    CHECK(std::any_of(
        channel.g_tensor.begin(), channel.g_tensor.end(),
        [](const double coupling) { return !gra::math::IsZero(coupling); }));
    CHECK(std::all_of(channel.ff_transfer.param.begin(),
                      channel.ff_transfer.param.end(),
                      [](const double value) { return std::isfinite(value); }));
  }

  SECTION("Tensor axial binary cascades reject nested three-body helicity vertices") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "TP", "RES", "omega(782)0 > {pi+ pi- pi0} gamma");
    MRandom rng;
    auto f1 = gra::resonance::Read("RES/f1_1420.json", rng, gra::ReggeProductionModel::TP);
    SetToyTensorChannel(f1, {1.0, 0.0});
    proc.SetResonances({{"f1_1420", f1}});
    REQUIRE_THROWS_WITH(proc.InitializeProcessAmplitude(),
                        Catch::Contains("physical helicity decay vertices must be two-body"));
  }

  SECTION("Tensor axial resonance builds its two-body helicity decay") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "TP", "RES", "J/psi(1S)0 gamma");
    MRandom rng;
    auto chi_c1 = gra::resonance::Read("RES/chi_c1.json", rng, gra::ReggeProductionModel::TP);
    SetToyTensorChannel(chi_c1, {1.0, 0.0});
    proc.SetResonances({{"chi_c1", chi_c1}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto configured = proc.GetResonances().at("chi_c1");
    CHECK(configured.hel_decay.UsesHelicityCouplings());
    CHECK(configured.hel_decay.T.FrobNorm2() == Approx(1.0).epsilon(1e-12));
    REQUIRE_FALSE(configured.hel_decay.ls_components.empty());
    gra::MMatrix<std::complex<double>> reconstructed(
        configured.hel_decay.T.size_row(), configured.hel_decay.T.size_col(),
        0.0);
    for (const auto &component : configured.hel_decay.ls_components) {
      reconstructed += component.matrix * component.alpha;
    }
    RequireMatrixNear(reconstructed, configured.hel_decay.T, 1e-12);
    const auto pole =
        gra::spin::DecayLSHelicityMatrix(configured.hel_decay, 1.0, true);
    RequireMatrixNear(pole, configured.hel_decay.T, 1e-12);
    REQUIRE_NOTHROW(
        gra::spin::DecayLSHelicityMatrix(configured.hel_decay, 0.73, true));
    REQUIRE(gra::spin::DecayLSIntensity(configured.hel_decay, 0.73, true) >=
            0.0);
  }

  SECTION("MP direct 4-body matrix elements reject cascades") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "CON", "pi+ pi- K+ K-");
    proc.state.lts.decaytree[0].legs = {proc.state.lts.decaytree[0],
                                        proc.state.lts.decaytree[1]};
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("Tensor Pomeron initializes every direct channel with excitation") {
    for (const std::string channel : {"RES", "CON", "RES+CON"}) {
      for (const int nstars : {1, 2}) {
        CAPTURE(channel, nstars);
        ToyHelicityProcess proc;
        ConfigureToyProductionProcess(proc, "TP", channel, "pi+ pi-");
        if (channel != "CON") {
          MRandom rng;
          const auto f0 = gra::resonance::Read("RES/f0_980.json", rng, gra::ReggeProductionModel::TP);
          proc.SetResonances({{"f0_980", f0}});
        }
        proc.SetExcitation(nstars);
        REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
      }
    }
  }

  SECTION("Full QED initializes with single and double excitation") {
    for (const int nstars : {1, 2}) {
      CAPTURE(nstars);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "yy", "QED", "mu+ mu-");
      proc.SetExcitation(nstars);
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());
    }
  }

  SECTION("Tensor resonance rejects one-sided vector cascades") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "TP", "RES",
                                  "rho(770)0 > {pi+ pi-} rho(770)0");
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("Tensor resonance rejects deeper cascades") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(
        proc, "TP", "RES",
        "rho(770)0 > {pi0 > {gamma gamma} pi0} rho(770)0 > {pi+ pi-}");
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("combined Tensor resonance-continuum is direct two-body only") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(
        proc, "TP", "RES+CON", "rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}");
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("isolated continuum rejects its daughter-dependent matrix element") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "CON", "pi+ pi-");
    proc.SetRootDecayMode(gra::RootDecayMode::Isolated);
    proc.SetISOLATE(true);
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }

  SECTION("isolated non-axial Tensor resonance rejects intertwined decay "
          "vertices") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "TP", "RES", "pi+ pi-");
    proc.SetRootDecayMode(gra::RootDecayMode::Isolated);
    proc.SetISOLATE(true);
    MRandom rng;
    const auto f0 = gra::resonance::Read("RES/f0_980.json", rng, gra::ReggeProductionModel::TP);
    proc.SetResonances({{"f0_980", f0}});
    REQUIRE_THROWS(proc.InitializeProcessAmplitude());
  }
}

TEST_CASE("Isolated resonance setup retains a physical decay state",
          "[gra::MProcess][resonance][isolated]") {
  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
  proc.SetRootDecayMode(gra::RootDecayMode::Isolated);
  proc.SetISOLATE(true);

  MRandom rng;
  const auto f0 = gra::resonance::Read("RES/f0_980.json", rng, gra::ReggeProductionModel::MP);
  proc.SetResonances({{"f0_980", f0}});

  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

  const auto configured = proc.GetResonances().at("f0_980");
  CHECK(configured.hel_decay.BR_set);
  CHECK(configured.hel_decay.BR > 0.0);
  CHECK(configured.hel_decay.BR <= 1.0);
  CHECK(std::isfinite(std::real(configured.hel_decay.g_decay)));
  CHECK(std::isfinite(std::imag(configured.hel_decay.g_decay)));
  CHECK(std::abs(configured.hel_decay.g_decay) > 0.0);
}

// Check chi_cJ radiative cards produce finite normalized helicity matrices
TEST_CASE("chi_cJ radiative decays produce normalized helicity matrices",
          "[gra::spin][chic][decay]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  ToyHelicityProcess proc;
  proc.state.lts.PDG = LoadedPDGTable();
  const auto jpsi = proc.state.lts.PDG.FindByPDG(443);
  const auto gamma = proc.state.lts.PDG.FindByPDG(22);

  const auto chic0 = proc.state.lts.PDG.FindByPDG(10441);
  const auto hel0 = proc.ProcessHelicityStructure(
      chic0, {jpsi, gamma}, false, true, "chi_c0 E1 test", false);
  REQUIRE(hel0.UsesHelicityCouplings());
  REQUIRE(hel0.T.size_row() > 0);
  REQUIRE(hel0.T.size_col() > 0);
  REQUIRE(hel0.T.FrobNorm2() == Approx(1.0).margin(1e-12));

  const auto chic1 = proc.state.lts.PDG.FindByPDG(20443);
  const auto hel1 = proc.ProcessHelicityStructure(
      chic1, {jpsi, gamma}, false, true, "chi_c1 multipole test", false);
  REQUIRE(hel1.UsesHelicityCouplings());
  REQUIRE(hel1.T.size_row() > 0);
  REQUIRE(hel1.T.size_col() > 0);
  REQUIRE(hel1.T.FrobNorm2() == Approx(1.0).margin(1e-12));

  const auto chic2 = proc.state.lts.PDG.FindByPDG(445);
  const auto hel2 = proc.ProcessHelicityStructure(
      chic2, {jpsi, gamma}, false, true, "chi_c2 multipole test", false);
  REQUIRE(hel2.UsesHelicityCouplings());
  REQUIRE(hel2.T.size_row() > 0);
  REQUIRE(hel2.T.size_col() > 0);
  REQUIRE(hel2.T.FrobNorm2() == Approx(1.0).margin(1e-12));
}

TEST_CASE("Decay LS and helicity readers preserve common coupling phases",
          "[gra::spin][decay][phase]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  ToyHelicityProcess baseline_process;
  baseline_process.state.lts.PDG = LoadedPDGTable();
  const auto chic1 = baseline_process.state.lts.PDG.FindByPDG(20443);
  const auto jpsi = baseline_process.state.lts.PDG.FindByPDG(443);
  const auto gamma = baseline_process.state.lts.PDG.FindByPDG(22);
  const auto f0 = baseline_process.state.lts.PDG.FindByPDG(9010221);
  const auto kaon_plus = baseline_process.state.lts.PDG.FindByPDG(321);
  const auto kaon_minus = baseline_process.state.lts.PDG.FindByPDG(-321);
  const auto direct_baseline = baseline_process.ProcessHelicityStructure(
      chic1, {jpsi, gamma}, false, true, "decay phase baseline", false);
  const auto ls_baseline = baseline_process.ProcessHelicityStructure(
      f0, {kaon_plus, kaon_minus}, false, true, "decay LS phase baseline",
      false);

  constexpr double direct_phase = 0.43;
  constexpr double direct_zeta = -0.27;
  constexpr double ls_phase = -0.51;
  constexpr double ls_zeta = 0.18;
  const auto tune =
      WriteModifiedPhotoVMTune("decay_common_phase", [](auto &) {});
  const auto model = gra::MModelTune::Load(tune.second);
  const std::filesystem::path decays_path =
      std::filesystem::path(tune.first) / "DECAYS.json";
  auto decays =
      nlohmann::json::parse(gra::aux::GetInputData(decays_path.string()));
  auto &direct = decays["20443"]["[443,22]"];
  for (auto &row : direct["helicity"]) {
    row[3] = direct_phase;
  }
  direct["zeta"]["MP"] = direct_zeta;
  auto &ls = decays["9010221"]["[321,-321]"];
  ls["alpha_ls"][0][3] = ls_phase;
  ls["zeta"]["MP"] = ls_zeta;
  std::ofstream output(decays_path);
  REQUIRE(output.good());
  output << decays.dump(2);
  output.close();

  ToyHelicityProcess phased_process;
  phased_process.SetHelicityConfig(model);
  phased_process.SetModelTune(model);
  const auto direct_phased = phased_process.ProcessHelicityStructure(
      chic1, {jpsi, gamma}, false, true, "decay free helicity phase", false);
  const auto ls_phased = phased_process.ProcessHelicityStructure(
      f0, {kaon_plus, kaon_minus}, false, true, "decay free LS phase", false);

  auto direct_expected = direct_baseline.T;
  direct_expected *= std::polar(1.0, direct_phase);
  RequireMatrixNear(direct_phased.T, direct_expected, 1.0e-13);
  CHECK(direct_phased.zeta == Approx(direct_zeta).margin(1.0e-14));
  REQUIRE(ls_baseline.alpha_ls.Contains(0, 0));
  REQUIRE(ls_phased.alpha_ls.Contains(0, 0));
  RequireComplexNear(ls_phased.alpha_ls.At(0, 0),
                     ls_baseline.alpha_ls.At(0, 0) *
                         std::polar(1.0, ls_phase),
                     1.0e-13);
  CHECK(ls_phased.zeta == Approx(ls_zeta).margin(1.0e-14));
}

// Check transversality, parity and the diphoton width for the production amplitudes
TEST_CASE("Gamma-gamma resonance production is fixed by the diphoton width",
          "[gra::spin][resonance][gamma-gamma][normalization]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  MRandom rng;
  const auto nominal_f2 = gra::resonance::Read("RES/f2_1270_yy.json", rng, gra::ReggeProductionModel::GP);

  const double branching_ratio =
      gra::resonance::GammaGammaBranchingRatio(nominal_f2.p.pdg);
  const double partial_width = nominal_f2.p.width * branching_ratio;
  const double production_coupling =
      std::sqrt(32.0 * gra::math::PI * nominal_f2.p.mass *
                (nominal_f2.p.spinX2 + 1) * partial_width);
  CHECK(branching_ratio > 0.0);
  CHECK(branching_ratio < 1.0);
  CHECK(gra::resonance::GammaGammaPartialWidth(nominal_f2.p) ==
        Approx(partial_width));
  CHECK(gra::resonance::GammaGammaResonanceCoupling(nominal_f2.p) ==
        Approx(production_coupling));
  for (const auto *model :
       {static_cast<const gra::RES_PRODUCTION_MODEL *>(&nominal_f2.MP),
        static_cast<const gra::RES_PRODUCTION_MODEL *>(&nominal_f2.XP),
        static_cast<const gra::RES_PRODUCTION_MODEL *>(&nominal_f2.GP)}) {
    REQUIRE(model->channels.size() == 1);
    CHECK(model->channels.front().width_derived);
  }

  SECTION(
      "physical decay replaces the null card coupling by the width coupling") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetResonances({{"f2_1270_yy", nominal_f2}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = proc.GetResonances().at("f2_1270_yy");
    REQUIRE(configured.production.size() == 1);
    const auto pole = gra::spin::PoleLSReduced(
        configured.production.front().pole.value(), 0.5 * configured.p.mass);
    const double physical_norm2 = pole.FrobNorm2();
    CHECK(physical_norm2 == Approx(production_coupling * production_coupling));
    REQUIRE(pole.size_row() == 3);
    REQUIRE(pole.size_col() == 3);
    for (std::size_t i = 0; i < 3; ++i) {
      CHECK(std::abs(pole[1][i]) == Approx(0.0).margin(1e-14));
      CHECK(std::abs(pole[i][1]) == Approx(0.0).margin(1e-14));
    }
    CHECK(pole.FrobNorm2() ==
          Approx(2.0 * (std::norm(pole[0][0]) + std::norm(pole[0][2]))).epsilon(1e-12));
    RequireComplexNear(pole[0][0], pole[2][2], 1e-14);
    RequireComplexNear(pole[0][2], pole[2][0], 1e-14);
  }

  SECTION("derived magnitude retains the card phase") {
    std::vector<gra::HelAmp> reference;
    for (const double phase : {0.0, 0.37, -1.2, 2.4}) {
      CAPTURE(phase);
      ToyHelicityProcess proc;
      ConfigureToyProductionProcess(proc, "XP", "RES", "pi+ pi-");
      auto f2 = nominal_f2;
      const auto rotation = std::polar(1.0, phase);
      for (auto &term : f2.XP.channels.front().g_ls) { term.coefficient *= rotation; }
      for (auto &coupling : f2.XP.channels.front().g_helicity) { coupling *= rotation; }
      proc.SetResonances({{"f2_1270_yy", f2}});
      REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

      const auto &configured = proc.GetResonances().at("f2_1270_yy");
      REQUIRE(configured.production.size() == 1);
      const auto pole = gra::spin::PoleLSReduced(configured.production.front().pole.value(), 0.5 * configured.p.mass);
      CHECK(pole.FrobNorm2() == Approx(production_coupling * production_coupling));
      std::vector<gra::HelAmp> amplitudes;
      for (const double momentum : {0.2, 0.5 * configured.p.mass, 1.0, 1.7}) {
        amplitudes.push_back(gra::spin::PoleLSReduced(configured.production.front().pole.value(), momentum));
      }
      if (reference.empty()) { reference = amplitudes; }
      for (const auto &i : indices(amplitudes)) {
        CAPTURE(i);
        RequireMatrixNear(amplitudes[i], reference[i] * rotation, 1e-12);
      }
    }
  }

  SECTION(
      "isolated decay includes the normalized Breit-Wigner spectral residue") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "MP", "RES", "pi+ pi-");
    proc.SetRootDecayMode(gra::RootDecayMode::Isolated);
    proc.SetISOLATE(true);
    proc.SetResonances({{"f2_1270_yy", nominal_f2}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = proc.GetResonances().at("f2_1270_yy");
    const double expected =
        production_coupling *
        std::sqrt(2.0 * nominal_f2.p.mass * nominal_f2.p.width);
    REQUIRE(configured.production.size() == 1);
    const double physical_norm2 =
        gra::spin::PoleLSReduced(configured.production.front().pole.value(),
                                      0.5 * configured.p.mass)
            .FrobNorm2();
    CHECK(physical_norm2 == Approx(expected * expected));
    CHECK(configured.hel_decay.g_decay == std::complex<double>(1.0, 0.0));
  }

  SECTION("GP diphoton helicities retain the width normalization and parity") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    proc.SetResonances({{"f2_1270_yy", nominal_f2}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const gra::PARAM_RES configured = proc.GetResonances().at("f2_1270_yy");
    REQUIRE(configured.production.size() == 1);
    const auto &runtime = configured.production.front().hel;
    CHECK(runtime.UsesHelicityCouplings());
    CHECK(runtime.UsesReggeDomain());
    const auto negative = gra::gpom::AnalyticMIndex(-1, runtime.analytic_MMAX, "diphoton negative helicity");
    const auto positive = gra::gpom::AnalyticMIndex(1, runtime.analytic_MMAX, "diphoton positive helicity");
    CHECK(runtime.T.MaskedSquaredNorm(runtime.T_set) == Approx(production_coupling * production_coupling));
    CHECK(std::abs(runtime.T[negative][negative]) == Approx(0.0).margin(1e-14));
    CHECK(std::abs(runtime.T[positive][positive]) == Approx(0.0).margin(1e-14));
    CHECK(runtime.T.FrobNorm2() == Approx(2.0 * std::norm(runtime.T[negative][positive])).epsilon(1e-12));
    RequireComplexNear(runtime.T[negative][positive], runtime.T[positive][negative], 1e-14);
  }
}

TEST_CASE("DECAYS reader accepts reversed two-body particle-antiparticle order",
          "[gra::spin]") {
  ToyHelicityProcess proc;
  proc.state.lts.PDG = LoadedPDGTable();

  const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
  const auto pi_plus = proc.state.lts.PDG.FindByPDG(211);
  const auto pi_minus = proc.state.lts.PDG.FindByPDG(-211);

  const auto canonical = proc.ProcessHelicityStructure(
      rho0, {pi_plus, pi_minus}, false, true, "", false);
  const auto reversed = proc.ProcessHelicityStructure(rho0, {pi_minus, pi_plus},
                                                      false, true, "", false);

  CHECK(canonical.BR == Approx(1.0));
  CHECK(reversed.BR == Approx(1.0));
  CHECK(std::abs(canonical.alpha_ls.At(1, 0) + reversed.alpha_ls.At(1, 0)) <
        1e-12);
}

TEST_CASE("Decay LS construction drops inactive coefficients before caching",
          "[gra::spin][decay][LS][numerics]") {
  const auto tune = WriteModifiedPhotoVMTune("decay_ls_cutoff", [](auto &) {});
  const auto model = gra::MModelTune::Load(tune.second);
  const double cutoff = model->Global().coupling_min;

  const std::filesystem::path decays_path =
      std::filesystem::path(tune.first) / "DECAYS.json";
  auto decays =
      nlohmann::json::parse(gra::aux::GetInputData(decays_path.string()));
  decays["9080225"]["[113,113]"]["alpha_ls"] = {{0, 2, 2.0 * cutoff, 0.0},
                                                {2, 0, 0.5 * cutoff, 0.0}};
  std::ofstream output(decays_path);
  REQUIRE(output.good());
  output << decays.dump(2);
  output.close();

  ToyHelicityProcess proc;
  proc.SetHelicityConfig(model);
  proc.SetModelTune(model);
  const auto mother = proc.state.lts.PDG.FindByPDG(9080225);
  const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
  const auto hel = proc.ProcessHelicityStructure(mother, {rho0, rho0}, false,
                                                 true, "", false);

  REQUIRE(hel.alpha_ls.Size() == 1);
  CHECK(hel.alpha_ls.Contains(0, 4));
  CHECK_FALSE(hel.alpha_ls.Contains(2, 0));
}

TEST_CASE("Inactive coupling rows are validated before sparse reduction",
          "[gra::spin][GP][coupling][validation]") {
  ModelParamRestoreGuard restore;

  SECTION("a forbidden zero decay LS row is rejected") {
    const auto tune =
        WriteModifiedPhotoVMTune("decay_ls_zero_forbidden", [](auto &) {});
    const auto model = gra::MModelTune::Load(tune.second);
    const double cutoff = model->Global().coupling_min;
    const std::filesystem::path decays_path =
        std::filesystem::path(tune.first) / "DECAYS.json";
    auto decays =
        nlohmann::json::parse(gra::aux::GetInputData(decays_path.string()));
    decays["9080225"]["[113,113]"]["alpha_ls"] = {
        {0, 2, 2.0 * cutoff, 0.0}, {0, 3, 0.0, 0.0}};
    std::ofstream output(decays_path);
    REQUIRE(output.good());
    output << decays.dump(2);
    output.close();

    ToyHelicityProcess proc;
    proc.SetHelicityConfig(model);
    proc.SetModelTune(model);
    const auto mother = proc.state.lts.PDG.FindByPDG(9080225);
    const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
    bool forbidden = false;
    try {
      static_cast<void>(proc.ProcessHelicityStructure(
          mother, {rho0, rho0}, false, true, "", false));
    } catch (const std::invalid_argument &error) {
      forbidden = std::string(error.what()).find("forbidden LS row") !=
                  std::string::npos;
    }
    CHECK(forbidden);
  }

  SECTION("zero crossed helicity orbits are accepted then removed") {
    const std::string tune = WriteModifiedContinuumTune(
        "gp_crossed_zero_sparse", "GP", [](auto &card) {
          card.at("990").at("[113,113]").at("self").at("helicity") = {
              {-1, -1, 0, 1.0, 0.0}, {-1, 0, 0, 0.0, 0.0}};
        });
    const auto model = gra::MModelTune::Load(tune + "/GENERAL.json");
    gra::MODELPARAM = tune;
    ToyHelicityProcess proc;
    proc.SetProcessForTest("GP", "CON");
    proc.SetHelicityConfig(model);
    proc.SetModelTune(model);
    const auto trajectory = proc.state.lts.PDG.FindByPDG(990);
    const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
    const auto hel = proc.ProcessHelicityStructure(
        trajectory, {rho0, rho0}, true, true, "", false,
        gra::spin::VertexContext::SubTUChannelExchange);

    std::size_t active = 0;
    for (std::size_t row = 0; row < hel.T_set.size_row(); ++row) {
      for (std::size_t column = 0; column < hel.T_set.size_col(); ++column) {
        active += hel.T_set[row][column] ? 1U : 0U;
      }
    }
    CHECK(active == 2);
  }

  SECTION("an invalid zero crossed helicity row is still rejected") {
    const std::string tune = WriteModifiedContinuumTune(
        "gp_crossed_zero_forbidden", "GP", [](auto &card) {
          card.at("990").at("[113,113]").at("self").at("helicity") = {
              {-2, -2, 0, 0.0, 0.0}};
        });
    const auto model = gra::MModelTune::Load(tune + "/GENERAL.json");
    gra::MODELPARAM = tune;
    ToyHelicityProcess proc;
    proc.SetProcessForTest("GP", "CON");
    proc.SetHelicityConfig(model);
    proc.SetModelTune(model);
    const auto trajectory = proc.state.lts.PDG.FindByPDG(990);
    const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
    bool forbidden = false;
    try {
      static_cast<void>(proc.ProcessHelicityStructure(
          trajectory, {rho0, rho0}, true, true, "", false,
          gra::spin::VertexContext::SubTUChannelExchange));
    } catch (const std::invalid_argument &error) {
      forbidden = std::string(error.what()).find("outside spin basis") !=
                  std::string::npos;
    }
    CHECK(forbidden);
  }

  SECTION("zero GP fusion helicity orbits are accepted then removed") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto rho = gra::resonance::Read("RES/rho_770.json", proc.state.random, gra::ReggeProductionModel::GP);
    auto &channel = rho.GP.channels.front();
    channel.helicity = {{-1.0, 0.0}, {-1.0, -1.0}};
    channel.g_helicity = {1.0, 0.0};
    proc.SetResonances({{"rho_770", rho}});
    REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

    const auto &vertices = proc.GetResonances().at("rho_770").production;
    REQUIRE_FALSE(vertices.empty());
    for (const auto &production : vertices) {
      const auto &vertex = production.hel;
      std::size_t active = 0;
      for (std::size_t row = 0; row < vertex.T_set.size_row(); ++row) {
        for (std::size_t column = 0; column < vertex.T_set.size_col();
             ++column) {
          active += vertex.T_set[row][column] ? 1U : 0U;
        }
      }
      CHECK(active == 2);
    }
  }

  SECTION("an invalid all-zero GP fusion helicity table is still rejected") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto rho = gra::resonance::Read("RES/rho_770.json", proc.state.random, gra::ReggeProductionModel::GP);
    auto &channel = rho.GP.channels.front();
    channel.helicity = {{0.0, 0.0}};
    channel.g_helicity = {0.0};
    proc.SetResonances({{"rho_770", rho}});
    bool forbidden = false;
    try {
      proc.InitializeProcessAmplitude();
    } catch (const std::invalid_argument &error) {
      forbidden = std::string(error.what()).find(
                      "GP photon helicity must be transverse") !=
                  std::string::npos;
    }
    CHECK(forbidden);
  }

  SECTION("an invalid all-zero GP fusion LS table is still rejected") {
    ToyHelicityProcess proc;
    ConfigureToyProductionProcess(proc, "GP", "RES", "pi+ pi-");
    auto f2 = gra::resonance::Read("RES/f2_2150.json", proc.state.random, gra::ReggeProductionModel::GP);
    auto &coupling = f2.GP.channels.front().g_ls;
    coupling.Clear();
    coupling.Set(9, 0, 0.0);
    proc.SetResonances({{"f2_2150", f2}});
    CHECK_THROWS_AS(proc.InitializeProcessAmplitude(), std::invalid_argument);
  }
}

TEST_CASE("DECAYS reader requires explicit electromagnetic decay channels",
          "[gra::spin]") {
  ToyHelicityProcess proc;
  proc.state.lts.PDG = LoadedPDGTable();

  const auto rho0 = proc.state.lts.PDG.FindByPDG(113);
  const auto mu_minus = proc.state.lts.PDG.FindByPDG(13);
  const auto mu_plus = proc.state.lts.PDG.FindByPDG(-13);
  const auto e_minus = proc.state.lts.PDG.FindByPDG(11);
  const auto e_plus = proc.state.lts.PDG.FindByPDG(-11);

  const auto physical = proc.ProcessHelicityStructure(rho0, {mu_plus, mu_minus},
                                                      false, true, "", false);
  REQUIRE(physical.BR_set);
  REQUIRE(physical.BR > 0.0);
  REQUIRE(physical.BR < 1.0);
  REQUIRE(physical.UsesHelicityCouplings());
  CHECK(std::abs(physical.T[0][0]) < 1e-12);
  CHECK(std::abs(physical.T[0][1] - 1.0 / std::sqrt(2.0)) < 1e-12);
  CHECK(std::abs(physical.T[1][0] - 1.0 / std::sqrt(2.0)) < 1e-12);
  CHECK(std::abs(physical.T[1][1]) < 1e-12);
  RequireReducedHelicityNormalization(physical);

  // Physical lepton channels must not supply couplings for an unlisted species
  auto unlisted_minus = e_minus;
  auto unlisted_plus = e_plus;
  unlisted_minus.pdg = 9000991;
  unlisted_plus.pdg = -9000991;
  CHECK_THROWS_AS(proc.ProcessHelicityStructure(rho0, {unlisted_plus, unlisted_minus},
                                               false, false, "", false), gra::MissingHelicityData);
}

TEST_CASE(
    "BR-derived two-body decay coupling includes identical daughter factor",
    "[gra::spin][symmetry][physics]") {
  ToyHelicityProcess proc;
  proc.state.lts.PDG = LoadedPDGTable();

  auto f0 = proc.state.lts.PDG.FindByPDG(10331);
  f0.mass = 2.5;
  f0.width = 0.15;
  const auto rho0 = proc.state.lts.PDG.FindByPDG(113);

  const auto hc =
      proc.ProcessHelicityStructure(f0, {rho0, rho0}, false, true, "", false);
  const double identical_symmetry = 2.0;
  const double ps = gra::kinematics::PDW2body(
      gra::math::pow2(f0.mass), gra::math::pow2(rho0.mass),
      gra::math::pow2(rho0.mass), 1.0, identical_symmetry);
  const double expected = std::sqrt((f0.spinX2 + 1) * hc.BR * f0.width / ps);

  CHECK(hc.BR > 0.0);
  CHECK(hc.BR <= 1.0);
  CHECK(std::abs(hc.g_decay) == Approx(expected).epsilon(1e-12));
}

TEST_CASE(
    "gra::spin production vertices keep physical and auxiliary legs distinct",
    "[gra::spin]") {
  ToyHelicityProcess proc;
  ConfigureToyProductionProcess(proc, "MP", "CON", "p+ p-");

  const auto pomeron1 = proc.state.lts.PDG.FindByPDG(993);
  const auto odderon1 = proc.state.lts.PDG.FindByPDG(9993);
  const auto proton = proc.state.lts.PDG.FindByPDG(2212);
  const auto antiproton = proc.state.lts.PDG.FindByPDG(-2212);
  const auto neutron = proc.state.lts.PDG.FindByPDG(2112);
  const auto antineutron = proc.state.lts.PDG.FindByPDG(-2112);
  const auto pion = proc.state.lts.PDG.FindByPDG(211);
  const auto antipion = proc.state.lts.PDG.FindByPDG(-211);

  SECTION("forward proton source constructs auxiliary p pbar metadata") {
    const auto hel =
        gra::ForwardHadronSourceHelicityStructure(pomeron1, {proton, proton});

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.5));
    CHECK(hel.s2 == Approx(0.5));
    CHECK(hel.T.size_row() == 0);
  }

  SECTION("crossed antiproton leg uses the charge-conjugate pp row") {
    const auto hel = gra::ForwardHadronSourceHelicityStructure(
        pomeron1, {antiproton, proton});

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.5));
    CHECK(hel.s2 == Approx(0.5));
    CHECK(hel.T.size_row() == 0);
  }

  SECTION("C-odd crossed antiproton leg uses the charge-conjugate pp row") {
    const auto hel = gra::ForwardHadronSourceHelicityStructure(
        odderon1, {antiproton, proton});

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.5));
    CHECK(hel.s2 == Approx(0.5));
    CHECK(hel.T.size_row() == 0);
  }

  SECTION("crossed antineutron leg uses the charge-conjugate nn row") {
    const auto hel = gra::ForwardHadronSourceHelicityStructure(
        pomeron1, {antineutron, neutron});

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.5));
    CHECK(hel.s2 == Approx(0.5));
    CHECK(hel.T.size_row() == 0);
  }

  SECTION("physical pp subvertex selects same-sector LS data before crossing") {
    const auto hel = proc.ProcessHelicityStructure(
        pomeron1, {proton, proton}, true, true, "", false,
        gra::spin::VertexContext::SubTUChannelExchange);

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.5));
    CHECK(hel.s2 == Approx(0.5));
    CHECK(hel.T.FrobNorm2() > 0.0);
  }

  SECTION("production helicity tree constructs auxiliary p pbar metadata") {
    gra::MDecayBranch branch;
    branch.p = pomeron1;
    branch.legs.resize(2);
    branch.legs[0].p = proton;
    branch.legs[1].p = proton;

    REQUIRE_NOTHROW(proc.ProcessHelicityTree(branch, true, true));
    CHECK_FALSE(branch.hel.UsesHelicityCouplings());
    CHECK(branch.hel.s1 == Approx(0.5));
    CHECK(branch.hel.s2 == Approx(0.5));
    CHECK(branch.hel.T.size_row() == 0);
  }

  SECTION("physical pion antipion subvertex selects opposite odd-L data") {
    const auto hel = proc.ProcessHelicityStructure(
        pomeron1, {pion, antipion}, true, true, "", false,
        gra::spin::VertexContext::SubTUChannelExchange);

    CHECK_FALSE(hel.UsesHelicityCouplings());
    CHECK(hel.s1 == Approx(0.0));
    CHECK(hel.s2 == Approx(0.0));
    CHECK(hel.T.FrobNorm2() > 0.0);
  }
}

TEST_CASE(
    "gra::spin fixed-spin forward exchange residues use Regge-helicity m basis",
    "[gra::spin]") {
  gra::PARAM_RES res;
  const gra::LORENTZSCALAR lts = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  const double qt = (lts.pbeam1 - lts.pfinal[1]).Pt();
  const double s0 = 1.0;

  SECTION("J=0 Regge-helicity m basis") {
    const auto hel = RealisticProtonLegHelicityMatrix(991, 0, 1, 0, 0);
    const auto f   = gra::spin::Forward(lts, ForwardBranchForTest(hel), lts.pbeam1, lts.pfinal[1], false,
                                        gra::spin::Rows(ForwardBranchForTest(hel), false), lts.process.PHOTON_VERTEX,
                                        gra::spin::ForwardSpec{});

    REQUIRE(f.size_row() == 4);
    REQUIRE(f.size_col() == 1);
    CHECK(std::abs(f[0][0]) == Approx(1.0));
    CHECK(std::abs(f[3][0]) == Approx(1.0));
    CHECK(std::abs(f[1][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[2][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
  }

  SECTION("J=1 Regge-helicity m basis") {
    const auto hel = RealisticProtonLegHelicityMatrix(993, 2, -1, 1, 2);
    const auto f   = gra::spin::Forward(lts, ForwardBranchForTest(hel), lts.pbeam1, lts.pfinal[1], false,
                                        gra::spin::Rows(ForwardBranchForTest(hel), false), lts.process.PHOTON_VERTEX,
                                        gra::spin::ForwardSpec{});

    REQUIRE(f.size_row() == 4);
    REQUIRE(f.size_col() == 3);
    CHECK(std::abs(f[0][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[0][1]) == Approx(1.0));
    CHECK(std::abs(f[0][2]) ==
          Approx(ReggeResidueScaleForTest(hel, 2, qt, s0)));
    CHECK(std::abs(f[3][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[3][1]) == Approx(1.0));
    CHECK(std::abs(f[3][2]) ==
          Approx(ReggeResidueScaleForTest(hel, 2, qt, s0)));
    CHECK(std::abs(f[1][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[1][1]) ==
          Approx(ReggeResidueScaleForTest(hel, 1, qt, s0)));
    CHECK(std::abs(f[1][2]) ==
          Approx(ReggeResidueScaleForTest(hel, 2, qt, s0)));
    CHECK(std::abs(f[2][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[2][1]) ==
          Approx(ReggeResidueScaleForTest(hel, 1, qt, s0)));
    CHECK(std::abs(f[2][2]) ==
          Approx(ReggeResidueScaleForTest(hel, 2, qt, s0)));
  }

  SECTION("J=2 Regge-helicity m basis") {
    const auto hel = RealisticProtonLegHelicityMatrix(995, 4, 1, 2, 0);
    const auto f   = gra::spin::Forward(lts, ForwardBranchForTest(hel), lts.pbeam1, lts.pfinal[1], false,
                                        gra::spin::Rows(ForwardBranchForTest(hel), false), lts.process.PHOTON_VERTEX,
                                        gra::spin::ForwardSpec{});

    REQUIRE(f.size_row() == 4);
    REQUIRE(f.size_col() == 5);
    CHECK(std::abs(f[0][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[0][1]) ==
          Approx(ReggeResidueScaleForTest(hel, 1, qt, s0)));
    CHECK(std::abs(f[0][2]) == Approx(1.0));
    CHECK(std::abs(f[0][3]) ==
          Approx(ReggeResidueScaleForTest(hel, 3, qt, s0)));
    CHECK(std::abs(f[0][4]) ==
          Approx(ReggeResidueScaleForTest(hel, 4, qt, s0)));
    CHECK(std::abs(f[1][0]) ==
          Approx(ReggeResidueScaleForTest(hel, 0, qt, s0)));
    CHECK(std::abs(f[1][1]) ==
          Approx(ReggeResidueScaleForTest(hel, 1, qt, s0)));
    CHECK(std::abs(f[1][2]) ==
          Approx(ReggeResidueScaleForTest(hel, 2, qt, s0)));
    CHECK(std::abs(f[1][3]) ==
          Approx(ReggeResidueScaleForTest(hel, 3, qt, s0)));
    CHECK(std::abs(f[1][4]) ==
          Approx(ReggeResidueScaleForTest(hel, 4, qt, s0)));
  }
}

// Compare factorized and dense forward-source contractions in both pole bases
TEST_CASE("gra::spin forward factors preserve dense bilinear amplitudes",
          "[gra::spin][forward][factorized]") {
  using Complex = std::complex<double>;
  const gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.28, 4.22, -4.41);

  gra::MDecayBranch upper_branch;
  upper_branch.p.name = "test_tensor_exchange";
  upper_branch.p.pdg = 995;
  upper_branch.p.spinX2 = 4;
  upper_branch.hel = RealisticProtonLegHelicityMatrix(995, 4, 1, 2, 0);

  gra::MDecayBranch lower_branch;
  lower_branch.p.name = "test_vector_exchange";
  lower_branch.p.pdg = 993;
  lower_branch.p.spinX2 = 2;
  lower_branch.hel = RealisticProtonLegHelicityMatrix(993, 2, -1, 1, 2);

  gra::M4Vec lower_axis = lts.q2_in_X;
  lower_axis.Flip3();
  for (const auto mode : {gra::ForwardVertexMode::HelicityResidue,
                          gra::ForwardVertexMode::UnitResidue}) {
    for (const bool reduced_regge : {false, true}) {
      for (const bool no_flip : {true, false}) {
        CAPTURE(mode, reduced_regge, no_flip);
        const auto upper_rows = gra::spin::Rows(upper_branch, no_flip);
        const auto lower_rows = gra::spin::Rows(lower_branch, no_flip);
        const gra::spin::ForwardSpec forward{
            mode, 1.7,
            reduced_regge ? gra::ExchangeBasisType::ReducedRegge
                           : gra::ExchangeBasisType::HelicityTransport};
        const auto upper_factors = gra::spin::ForwardFactors(upper_branch, lts.pbeam1, lts.pfinal[1], false, upper_rows, forward);
        const auto lower_factors = gra::spin::ForwardFactors(lower_branch, lts.pbeam2, lts.pfinal[2], true, lower_rows, forward);
        REQUIRE(upper_factors.has_value());
        REQUIRE(lower_factors.has_value());

        const std::size_t upper_columns =
            reduced_regge ? 1 : upper_branch.hel.Jz_values.size();
        const std::size_t lower_columns =
            reduced_regge ? 1 : lower_branch.hel.Jz_values.size();
        REQUIRE(upper_factors->exchange_helicity.size() == upper_columns);
        REQUIRE(lower_factors->exchange_helicity.size() == lower_columns);

        const auto upper_dense =
            gra::spin::Forward(lts, upper_branch, lts.pbeam1, lts.pfinal[1], false, upper_rows, "EPA", forward);
        const auto lower_dense =
            gra::spin::Forward(lts, lower_branch, lts.pbeam2, lts.pfinal[2], true, lower_rows, "EPA", forward);
        RequireMatrixNear(gra::OuterProduct(upper_factors->beam_transition,
                                            upper_factors->exchange_helicity),
                          upper_dense, 1e-13);
        RequireMatrixNear(gra::OuterProduct(lower_factors->beam_transition,
                                            lower_factors->exchange_helicity),
                          lower_dense, 1e-13);

        MMatrix<Complex> central(
            upper_dense.size_col() * lower_dense.size_col(), 3, 0.0);
        for (std::size_t row = 0; row < central.size_row(); ++row) {
          for (std::size_t col = 0; col < central.size_col(); ++col) {
            if ((row + 2 * col) % 4 != 0) {
              central[row][col] =
                  Complex(0.07 * static_cast<double>(1 + row + col),
                          -0.03 * static_cast<double>(1 + 2 * row + col));
            }
          }
        }
        const MMatrix<Complex> sub_t = central * Complex(0.73, -0.21);
        const MMatrix<Complex> sub_u = central * Complex(-0.18, 0.64);
        const Complex scale(0.61, -0.37);

        RequireMatrixNear(
            gra::spin::Contract(*upper_factors, *lower_factors, central, scale),
            gra::spin::Contract(upper_dense, lower_dense, central, scale),
            2e-13);
        const auto factorized_continuum =
            std::pair{gra::spin::Contract(*upper_factors, *lower_factors, sub_t), gra::spin::Contract(*upper_factors, *lower_factors, sub_u)};
        const auto dense_continuum =
            std::pair{gra::spin::Contract(upper_dense, lower_dense, sub_t), gra::spin::Contract(upper_dense, lower_dense, sub_u)};
        RequireMatrixNear(factorized_continuum.first, dense_continuum.first,
                          2e-13);
        RequireMatrixNear(factorized_continuum.second, dense_continuum.second,
                          2e-13);
        const std::vector<Complex> projected_t =
            sub_t.Transpose() *
            gra::KroneckerProduct(upper_factors->exchange_helicity,
                                  lower_factors->exchange_helicity);
        const std::vector<Complex> projected_u =
            sub_u.Transpose() *
            gra::KroneckerProduct(upper_factors->exchange_helicity,
                                  lower_factors->exchange_helicity);
        const auto projected_continuum = gra::spin::Contract(
            *upper_factors, *lower_factors, projected_t, projected_u);
        RequireMatrixNear(projected_continuum.first, dense_continuum.first,
                          2e-13);
        RequireMatrixNear(projected_continuum.second, dense_continuum.second,
                          2e-13);
      }
    }
  }

  const auto rows = gra::spin::Rows(upper_branch, true);
  gra::MDecayBranch photon_branch = upper_branch;
  photon_branch.p.pdg = gra::PDG::PDG_gamma;
  CHECK_FALSE(gra::spin::ForwardFactors(photon_branch, lts.pbeam1, lts.pfinal[1], false, rows, {gra::ForwardVertexMode::HelicityResidue, 1.0})
                  .has_value());
}

// Check that both continuum contraction routes use the same bilinear source
TEST_CASE("gra::spin continuum sources use one bilinear lower residue",
          "[gra::spin][forward][factorized][continuum][regression]") {
  using Complex = std::complex<double>;
  const gra::spin::ForwardSourceFactors upper = {
      {Complex(0.7, -0.2), Complex(-0.4, 0.6), Complex(0.2, 0.5),
       Complex(-0.6, -0.3)},
      {Complex(0.8, 0.1), Complex(-0.3, -0.5)}};
  const gra::spin::ForwardSourceFactors lower = {
      {Complex(0.6, 0.3), Complex(-0.2, 0.7), Complex(0.5, -0.4),
       Complex(-0.8, -0.1)},
      {Complex(-0.7, 0.2), Complex(0.3, -0.6), Complex(0.4, 0.5)}};
  MMatrix<Complex> sub_t(6, 2, 0.0);
  MMatrix<Complex> sub_u(6, 2, 0.0);
  for (std::size_t row = 0; row < sub_t.size_row(); ++row) {
    for (std::size_t col = 0; col < sub_t.size_col(); ++col) {
      sub_t[row][col] = Complex(0.09 * static_cast<double>(1 + row + col),
                                -0.04 * static_cast<double>(1 + 2 * row + col));
      sub_u[row][col] =
          sub_t[row][col] * Complex(-0.31 + 0.02 * row, 0.57 - 0.03 * col);
    }
  }

  const auto destination =
      gra::spin::CanonicalProtonPairSpinLayout::KroneckerDestinationRows(
          upper.beam_transition.size(), lower.beam_transition.size());
  const auto actual = std::pair{gra::spin::Contract(upper, lower, sub_t), gra::spin::Contract(upper, lower, sub_u)};
  const auto expected_t = sub_t.FactorizedKroneckerMultiply(
      upper.beam_transition, upper.exchange_helicity, lower.beam_transition,
      lower.exchange_helicity, destination);
  const auto expected_u = sub_u.FactorizedKroneckerMultiply(
      upper.beam_transition, upper.exchange_helicity, lower.beam_transition,
      lower.exchange_helicity, destination);
  RequireMatrixNear(actual.first, expected_t, 2.0e-13);
  RequireMatrixNear(actual.second, expected_u, 2.0e-13);

  const auto hermitian = sub_t.FactorizedKroneckerMultiply(
      upper.beam_transition, upper.exchange_helicity,
      gra::Conjugated(lower.beam_transition),
      gra::Conjugated(lower.exchange_helicity), destination);
  REQUIRE(MatrixDiffNorm2(actual.first, hermitian) > 1.0e-4);

  const std::vector<Complex> projected_t =
      sub_t.Transpose() *
      gra::KroneckerProduct(upper.exchange_helicity,
                            lower.exchange_helicity);
  const std::vector<Complex> projected_u =
      sub_u.Transpose() *
      gra::KroneckerProduct(upper.exchange_helicity,
                            lower.exchange_helicity);
  const auto projected =
      gra::spin::Contract(upper, lower, projected_t, projected_u);
  const auto beam = gra::MappedKroneckerProduct(upper.beam_transition,
                                                lower.beam_transition,
                                                destination);
  RequireMatrixNear(projected.first, gra::OuterProduct(beam, projected_t),
                    2.0e-13);
  RequireMatrixNear(projected.second, gra::OuterProduct(beam, projected_u),
                    2.0e-13);
}

TEST_CASE("gra::spin decay matrices transpose without conjugating JW phases",
          "[gra::spin]") {
  gra::HELMatrix hel;
  hel.J = 1.0;
  hel.s1 = 0.0;
  hel.s2 = 0.0;
  hel.Jz_values = {-1.0, 0.0, 1.0};
  hel.lambda_values = MMatrix<double>(1, 2, 0.0);
  hel.lambda_idx = MMatrix<std::size_t>(1, 2, 0);
  hel.T = MMatrix<std::complex<double>>(1, 1, 0.0);
  hel.T[0][0] = std::complex<double>(1.0, 0.25);
  gra::spin::InitJWRotation(hel);

  const auto f = gra::spin::fDecayMatrix(hel, 0.9, 0.7);
  const auto decay = f.Transpose();

  REQUIRE(decay.size_row() == f.size_col());
  REQUIRE(decay.size_col() == f.size_row());

  bool saw_complex_phase = false;
  for (std::size_t j = 0; j < f.size_col(); ++j) {
    REQUIRE(decay[j][0].real() == Approx(f[0][j].real()).epsilon(1e-12));
    REQUIRE(decay[j][0].imag() == Approx(f[0][j].imag()).epsilon(1e-12));
    if (std::abs(f[0][j].imag()) > 1e-9) {
      saw_complex_phase = true;
      REQUIRE(decay[j][0].imag() != Approx(-f[0][j].imag()).margin(1e-9));
    }
  }
  REQUIRE(saw_complex_phase);
}

TEST_CASE("gra::spin forward sources obey the collider helicity section",
          "[gra::spin][forward][helicity][phase]") {
  constexpr std::array<int, 4> transition_in_x2 = {-1, -1, 1, 1};
  constexpr std::array<int, 4> transition_out_x2 = {-1, 1, -1, 1};
  const auto rotate = [](gra::M4Vec p, double angle) {
    p.RotateZ(angle);
    return p;
  };

  gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.27, 4.25, -4.38);
  const auto vector_hel = RealisticProtonLegHelicityMatrix(993, 2, -1, 1, 2);
  const double rotation = 0.58;

  for (const bool lower : {false, true}) {
    CAPTURE(lower);
    const gra::M4Vec incoming = lower ? lts.pbeam2 : lts.pbeam1;
    const gra::M4Vec outgoing = lower ? lts.pfinal[2] : lts.pfinal[1];
    const auto       reference =
        gra::spin::Forward(lts, ForwardBranchForTest(vector_hel), incoming, outgoing, lower,
                           gra::spin::Rows(ForwardBranchForTest(vector_hel), false), lts.process.PHOTON_VERTEX,
                           gra::spin::ForwardSpec{gra::ForwardVertexMode::HelicityResidue, 1.0});
    const auto rotated = gra::spin::Forward(
        lts, ForwardBranchForTest(vector_hel), rotate(incoming, rotation), rotate(outgoing, rotation), lower,
        gra::spin::Rows(ForwardBranchForTest(vector_hel), false), lts.process.PHOTON_VERTEX,
        gra::spin::ForwardSpec{gra::ForwardVertexMode::HelicityResidue, 1.0});

    for (std::size_t row = 0; row < reference.size_row(); ++row) {
      const int delta = (transition_in_x2[row] - transition_out_x2[row]) / 2;
      for (std::size_t col = 0; col < reference.size_col(); ++col) {
        const int m = static_cast<int>(vector_hel.Jz_values[col]);
        const int harmonic = (lower ? -1 : 1) * (delta + m);
        const std::complex<double> expected =
            reference[row][col] *
            std::exp(gra::math::zi * static_cast<double>(harmonic) * rotation);
        CAPTURE(row, col, harmonic, reference[row][col], rotated[row][col]);
        REQUIRE(std::abs(rotated[row][col] - expected) < 2.0e-11);
      }
    }
  }

  gra::MDecayBranch branch;
  branch.p.name = "test_vector_exchange";
  branch.p.pdg = 993;
  branch.p.spinX2 = 2;
  branch.p.P = -1;
  branch.hel = vector_hel;
  const auto rows = gra::spin::Rows(branch, false);
  const auto rotate_lts = [rotation](gra::LORENTZSCALAR event) {
    event.pbeam1.RotateZ(rotation);
    event.pbeam2.RotateZ(rotation);
    for (auto &p : event.pfinal) {
      p.RotateZ(rotation);
    }
    event.q1.RotateZ(rotation);
    event.q2.RotateZ(rotation);
    event.q1_in_X.RotateZ(rotation);
    event.q2_in_X.RotateZ(rotation);
    return event;
  };
  const gra::LORENTZSCALAR rotated_lts = rotate_lts(lts);

  for (const bool lower : {false, true}) {
    CAPTURE(lower);
    const gra::LORENTZSCALAR &event = lts;
    const gra::LORENTZSCALAR &rotated_event = rotated_lts;
    gra::M4Vec axis = lower ? event.q2_in_X : event.q1_in_X;
    gra::M4Vec rotated_axis =
        lower ? rotated_event.q2_in_X : rotated_event.q1_in_X;
    if (lower) {
      axis.Flip3();
      rotated_axis.Flip3();
    }
    const auto reference = gra::spin::Forward(event, branch, lower ? event.pbeam2 : event.pbeam1, lower ? event.pfinal[2] : event.pfinal[1], lower, rows, "EPA", gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});
    const auto rotated_source = gra::spin::Forward(rotated_event, branch, lower ? rotated_event.pbeam2 : rotated_event.pbeam1, lower ? rotated_event.pfinal[2] : rotated_event.pfinal[1], lower, rows, "EPA", gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});

    for (std::size_t row = 0; row < reference.size_row(); ++row) {
      const int delta = (transition_in_x2[row] - transition_out_x2[row]) / 2;
      for (std::size_t col = 0; col < reference.size_col(); ++col) {
        const int m = static_cast<int>(vector_hel.Jz_values[col]);
        const int harmonic = (lower ? -1 : 1) * (delta + m);
        const std::complex<double> expected =
            reference[row][col] *
            std::exp(gra::math::zi * static_cast<double>(harmonic) * rotation);
        CAPTURE(lower, row, col, harmonic);
        REQUIRE(std::abs(rotated_source[row][col] - expected) < 3.0e-10);
      }
    }
  }
}

TEST_CASE("gra::spin factorized scalar source has signed pair reciprocity",
          "[gra::spin][forward][helicity][reciprocity]") {
  const auto scalar_hel = RealisticProtonLegHelicityMatrix(991, 0, 1, 0, 0);
  gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.25, 4.2, -4.2);
  const double momentum = 5.0;
  const double transverse = 0.31;
  const double longitudinal =
      std::sqrt(momentum * momentum - transverse * transverse);
  const double energy =
      std::sqrt(momentum * momentum + gra::PDG::mp * gra::PDG::mp);
  const gra::M4Vec p1_in(0.0, 0.0, momentum, energy);
  const gra::M4Vec p2_in(0.0, 0.0, -momentum, energy);
  const gra::M4Vec p1_out(transverse, 0.0, longitudinal, energy);
  const gra::M4Vec p2_out(-transverse, 0.0, -longitudinal, energy);
  const std::vector<std::size_t> rows = {0, 1, 2, 3};
  const auto                     upper =
      gra::spin::Forward(lts, ForwardBranchForTest(scalar_hel), p1_in, p1_out, false,
                         gra::spin::Rows(ForwardBranchForTest(scalar_hel), false), lts.process.PHOTON_VERTEX,
                         gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});
  const auto lower =
      gra::spin::Forward(lts, ForwardBranchForTest(scalar_hel), p2_in, p2_out, true,
                         gra::spin::Rows(ForwardBranchForTest(scalar_hel), false), lts.process.PHOTON_VERTEX,
                         gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});
  REQUIRE(upper.size_col() == 1);
  REQUIRE(lower.size_col() == 1);
  const std::array<double, 4> upper_sign = {1.0, -1.0, 1.0, 1.0};
  const std::array<double, 4> lower_sign = {1.0, 1.0, -1.0, 1.0};
  for (std::size_t row : rows) {
    REQUIRE(std::real(upper[row][0] / std::abs(upper[row][0])) ==
            Approx(upper_sign[row]));
    REQUIRE(std::real(lower[row][0] / std::abs(lower[row][0])) ==
            Approx(lower_sign[row]));
  }

  const auto hard = gra::spin::Contract(
      upper, lower, MMatrix<std::complex<double>>(1, 1, 1.0));
  constexpr auto h = gra::spin::BinaryHelicityLabelsX2();
  for (std::size_t out = 0; out < 4; ++out) {
    for (std::size_t in = 0; in < 4; ++in) {
      const int h1_in = h[in / 2];
      const int h2_in = h[in % 2];
      const int h1_out = h[out / 2];
      const int h2_out = h[out % 2];
      const double sign = gra::spin::ColliderSpinHalfReciprocitySign(
          h1_in, h2_in, h1_out, h2_out);
      const std::complex<double> direct =
          hard[gra::spin::PairHelicityTransitionIndex(in, out)][0];
      const std::complex<double> reverse =
          hard[gra::spin::PairHelicityTransitionIndex(out, in)][0];
      CAPTURE(out, in, sign, direct, reverse);
      REQUIRE(std::abs(direct - sign * reverse) < 2.0e-12);
    }
  }
}

TEST_CASE("gra::spin stable-leaf BW mixture includes crossed decay phase-space "
          "density",
          "[gra::spin]") {
  const double mX = 2.4;
  const double mx = 0.14;
  const double ma = 0.20;
  const double mb = 0.30;
  const double mR1 = 0.82;
  const double mR2 = 0.93;

  const auto R_in_X = TwoBodyRestKinematics(mX, mR1, mR2, 0.62, -0.27);
  const auto x1_a_in_R1 = TwoBodyRestKinematics(mR1, mx, ma, 0.74, 0.51);
  const auto x2_b_in_R2 = TwoBodyRestKinematics(mR2, mx, mb, 1.38, -0.81);

  gra::MDecayBranch x1;
  x1.name = "x1";
  x1.p = ToyParticle("x", 111, 0, mx);
  x1.p4 = BoostFromRestFrame(x1_a_in_R1[0], R_in_X[0]);

  gra::MDecayBranch a;
  a.name = "a";
  a.p = ToyParticle("a", 211, 0, ma);
  a.p4 = BoostFromRestFrame(x1_a_in_R1[1], R_in_X[0]);

  gra::MDecayBranch x2;
  x2.name = "x2";
  x2.p = ToyParticle("x", 111, 0, mx);
  x2.p4 = BoostFromRestFrame(x2_b_in_R2[0], R_in_X[1]);

  gra::MDecayBranch b;
  b.name = "b";
  b.p = ToyParticle("b", 321, 0, mb);
  b.p4 = BoostFromRestFrame(x2_b_in_R2[1], R_in_X[1]);

  gra::MDecayBranch R1;
  R1.name = "R1";
  R1.p = ToyParticle("R1", 900001, 0, mR1);
  R1.p.width = 0.08;
  R1.mass_proposal_norm = 2.0;
  R1.mass_proposal_min2 = 0.0;
  R1.mass_proposal_max2 = 100.0;
  R1.mass_proposal = gra::MassProposal::BreitWigner;
  R1.legs = {x1, a};
  R1.p4 = x1.p4 + a.p4;

  gra::MDecayBranch R2;
  R2.name = "R2";
  R2.p = ToyParticle("R2", 900002, 0, mR2);
  R2.p.width = 0.11;
  R2.mass_proposal_norm = 3.0;
  R2.mass_proposal_min2 = 0.0;
  R2.mass_proposal_max2 = 100.0;
  R2.mass_proposal = gra::MassProposal::BreitWigner;
  R2.legs = {x2, b};
  R2.p4 = x2.p4 + b.p4;

  gra::LORENTZSCALAR lts;
  lts.amplitude.DECAY_SYM = true;
  lts.decay_symmetry_proposal_active = true;
  lts.pfinal.resize(1);
  lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, mX);

  const std::vector<gra::MDecayBranch> tree = {R1, R2};

  auto two_body_ps = [](const gra::MDecayBranch &branch) {
    const double M0 = branch.p4.M();
    const double m0 = branch.legs[0].p4.M();
    const double m1 = branch.legs[1].p4.M();
    return gra::kinematics::dPhi2(M0,
                                  gra::kinematics::DecayMomentum(M0, m0, m1));
  };
  auto bw_density = [](const gra::MDecayBranch &left,
                       const gra::MDecayBranch &right) {
    return gra::math::abs2(gra::resonance::FixedWidthLineShape(
                               left.p4.M2(), left.p.mass, left.p.width) *
                           gra::resonance::FixedWidthLineShape(
                               right.p4.M2(), right.p.mass, right.p.width)) /
           (left.mass_proposal_norm * right.mass_proposal_norm);
  };

  std::vector<gra::MDecayBranch> crossed = tree;
  crossed[0].legs[0].p4 = x2.p4;
  crossed[1].legs[0].p4 = x1.p4;
  crossed[0].p4 = crossed[0].legs[0].p4 + crossed[0].legs[1].p4;
  crossed[1].p4 = crossed[1].legs[0].p4 + crossed[1].legs[1].p4;

  const double reference_ps = two_body_ps(tree[0]) * two_body_ps(tree[1]);
  const double crossed_ps = two_body_ps(crossed[0]) * two_body_ps(crossed[1]);
  REQUIRE(crossed_ps > 0.0);
  REQUIRE(crossed_ps != Approx(reference_ps).epsilon(1e-3));

  const double expected =
      0.5 * (bw_density(tree[0], tree[1]) +
             bw_density(crossed[0], crossed[1]) * reference_ps / crossed_ps);
  const double old_density =
      0.5 * (bw_density(tree[0], tree[1]) + bw_density(crossed[0], crossed[1]));

  const double density = gra::decay::MixtureDensity(lts, tree);
  REQUIRE(density == Approx(expected).epsilon(1e-12));
  REQUIRE(density != Approx(old_density).epsilon(1e-6));

  auto root_ps = [&lts](const std::vector<gra::MDecayBranch> &branches) {
    const double M0 = lts.pfinal[0].M();
    const double m0 = branches[0].p4.M();
    const double m1 = branches[1].p4.M();
    return gra::kinematics::dPhi2(M0,
                                  gra::kinematics::DecayMomentum(M0, m0, m1));
  };

  lts.PS_active = true;
  lts.DW = gra::kinematics::MCW(1.0, 1.0, 1.0);
  const double reference_total_ps = root_ps(tree) * reference_ps;
  const double crossed_total_ps = root_ps(crossed) * crossed_ps;
  REQUIRE(crossed_total_ps > 0.0);

  const double expected_with_root =
      0.5 * (bw_density(tree[0], tree[1]) + bw_density(crossed[0], crossed[1]) *
                                                reference_total_ps /
                                                crossed_total_ps);
  const double density_with_root =
      gra::decay::MixtureDensity(lts, tree);
  REQUIRE(density_with_root == Approx(expected_with_root).epsilon(1e-12));

  lts.decay_symmetry_proposal_phase_space = crossed_total_ps;
  const double expected_with_selected_proposal =
      0.5 *
      (bw_density(tree[0], tree[1]) * crossed_total_ps / reference_total_ps +
       bw_density(crossed[0], crossed[1]));
  const double density_with_selected_proposal =
      gra::decay::MixtureDensity(lts, tree);
  REQUIRE(density_with_selected_proposal ==
          Approx(expected_with_selected_proposal).epsilon(1e-12));
  lts.decay_symmetry_proposal_phase_space = 0.0;

  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Factorized;
  lts.central_phase_space_mass_cut_min = 0.5;
  lts.central_phase_space_mass_max = 4.0;
  lts.central_phase_space_mass_margin = 1.0e-4;
  const auto central_mass_jacobian =
      [&lts](const std::vector<gra::MDecayBranch> &branches) {
        double branch_mass_sum = 0.0;
        for (const auto &branch : branches) {
          branch_mass_sum += branch.p4.M();
        }
        const double lower =
            std::max(lts.central_phase_space_mass_cut_min,
                     branch_mass_sum + lts.central_phase_space_mass_margin);
        return 2.0 * lts.pfinal[0].M2() *
               std::log(lts.central_phase_space_mass_max / lower);
      };
  const double reference_mass_jacobian = central_mass_jacobian(tree);
  const double crossed_mass_jacobian = central_mass_jacobian(crossed);
  REQUIRE(reference_mass_jacobian > 0.0);
  REQUIRE(crossed_mass_jacobian > 0.0);
  REQUIRE(reference_mass_jacobian !=
          Approx(crossed_mass_jacobian).epsilon(1e-6));
  lts.central_phase_space_generated_jacobian = reference_mass_jacobian;

  const double expected_factorized =
      0.5 *
      (bw_density(tree[0], tree[1]) +
       bw_density(crossed[0], crossed[1]) * reference_total_ps /
           crossed_total_ps * reference_mass_jacobian / crossed_mass_jacobian);
  const double factorized_density =
      gra::decay::MixtureDensity(lts, tree);
  REQUIRE(factorized_density == Approx(expected_factorized).epsilon(1e-12));
  REQUIRE(factorized_density != Approx(density_with_root).epsilon(1e-7));

  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Collinear;
  const double collinear_density = gra::decay::MixtureDensity(lts, tree);
  REQUIRE(collinear_density == Approx(expected_factorized).epsilon(1e-12));
  REQUIRE(collinear_density != Approx(density_with_root).epsilon(1e-7));
  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::HardDiffraction;
  REQUIRE(gra::decay::MixtureDensity(lts, tree) == Approx(expected_factorized).epsilon(1e-12));
  lts.central_phase_space_generated_jacobian = 0.0;
  REQUIRE_THROWS_AS(gra::decay::MixtureDensity(lts, tree), gra::PhaseSpaceFailure);
}

TEST_CASE("gra::spin stable-leaf mixture enforces exact mass-proposal support",
          "[gra::spin][proposal]") {
  gra::LORENTZSCALAR lts = ContinuumVectorCascadeLTSForTest(false);
  lts.amplitude.DECAY_SYM = true;
  lts.decay_symmetry_proposal_active = true;

  for (auto &branch : lts.decaytree) {
    const double center = branch.p4.M2();
    branch.mass_proposal = gra::MassProposal::BreitWigner;
    branch.mass_proposal_min2 = center - 1e-8;
    branch.mass_proposal_max2 = center + 1e-8;
  }

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() == 2);
  REQUIRE((terms[1].tree[0].p4.M2() !=
               Approx(lts.decaytree[0].p4.M2()).margin(1e-8) ||
           terms[1].tree[1].p4.M2() !=
               Approx(lts.decaytree[1].p4.M2()).margin(1e-8)));

  const double reference_density =
      gra::math::abs2(gra::spin::CascadeBWProduct(lts.decaytree)) /
      CascadeMassProposalNormForTest(lts.decaytree);
  const double density =
      gra::decay::MixtureDensity(lts, lts.decaytree);
  REQUIRE(density == Approx(0.5 * reference_density).epsilon(1e-12));

  gra::LORENTZSCALAR boosted = lts;
  const double collision_mass = (lts.pbeam1 + lts.pbeam2).M();
  constexpr double beam_rapidity = 1.1;
  const gra::M4Vec boost(0.0, 0.0, collision_mass * std::sinh(beam_rapidity),
                         collision_mass * std::cosh(beam_rapidity));
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pbeam1, 1);
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pbeam2, 1);
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pfinal[0], 1);
  for (auto &branch : boosted.decaytree) {
    BoostDecayBranchForTest(branch, boost, collision_mass);
  }

  REQUIRE(boosted.pbeam1.E() != Approx(boosted.pbeam2.E()).epsilon(1e-8));
  const double boosted_density =
      gra::decay::MixtureDensity(boosted, boosted.decaytree);
  REQUIRE(boosted_density == Approx(density).epsilon(1e-11));
}

// Compute the two-root transverse density directly from d^2q = pi d(q^2)
double CascadeTransverseJacobianForTest(const gra::LORENTZSCALAR &lts,
                                      const std::vector<gra::MDecayBranch> &tree) {
  REQUIRE(tree.size() == 2);
  const double q2 = (tree[0].p4 - tree[1].p4).Pt2();
  const double scale = std::max(tree[0].p4.M() + tree[1].p4.M(), 0.01 * lts.central_phase_space_mass_max);
  return gra::math::PI * (q2 + scale * scale) *
         std::log1p(2.0 * lts.central_phase_space_transverse_radius2 / (scale * scale));
}

// Prepare a transverse proposal with full support for the deterministic cascade fixture
void SetTransverseProposalForTest(gra::LORENTZSCALAR &lts) {
  lts.central_phase_space_mass_max = 10.0;
  lts.central_phase_space_transverse_radius2 = 50.0 + 0.5 * lts.pfinal[0].Pt2();
  lts.central_phase_space_generated_jacobian = CascadeTransverseJacobianForTest(lts, lts.decaytree);
}

TEST_CASE("central cascade mixtures include each history transverse density",
          "[gra::spin][continuum][proposal]") {
  gra::LORENTZSCALAR lts = ContinuumVectorCascadeLTSForTest(false);
  lts.amplitude.DECAY_SYM = true;
  lts.decay_symmetry_proposal_active = true;
  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Central;
  SetTransverseProposalForTest(lts);
  lts.central_phase_space_rap_min_cm = -9.0;
  lts.central_phase_space_rap_max_cm = 9.0;

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() == 2);
  const double phase = InternalCascadePhaseSpaceForTest(lts.decaytree);
  double expected = 0.0;
  for (const auto &term : terms) {
    expected += MixedMassDensityForTest(term.tree) * phase / InternalCascadePhaseSpaceForTest(term.tree) *
                lts.central_phase_space_generated_jacobian / CascadeTransverseJacobianForTest(lts, term.tree);
  }
  expected /= terms.size();
  const double density = gra::decay::MixtureDensity(lts, lts.decaytree);
  REQUIRE(density == Approx(expected).epsilon(1e-12));
  auto omitted = lts;
  omitted.central_phase_space_mode = gra::CentralPhaseSpaceMode::Unknown;
  REQUIRE(density != Approx(gra::decay::MixtureDensity(omitted, omitted.decaytree)).epsilon(1e-6));
}

TEST_CASE(
    "MProcess central phase space stable leaf proposal excludes the polar "
    "factor and root decay measure",
    "[MProcess][continuum][proposal]") {
  const gra::LORENTZSCALAR lts = ContinuumVectorCascadeLTSForTest(false);
  const double expected = lts.decaytree[0].W_event * lts.decaytree[1].W_event;

  ToyHelicityProcess process;
  const double direct = process.AppliedStableLeafProposalPhaseSpaceForTest(
      lts.decaytree, gra::CentralPhaseSpaceMode::Central, true);
  const double factorized = process.AppliedStableLeafProposalPhaseSpaceForTest(
      lts.decaytree, gra::CentralPhaseSpaceMode::Unknown, true);

  REQUIRE(direct == Approx(expected).epsilon(1e-12));
  REQUIRE(factorized == Approx(0.37 * expected).epsilon(1e-12));
}

TEST_CASE("gra::spin central phase space mixture enforces top level rapidity "
          "support",
          "[gra::spin][continuum][proposal]") {
  constexpr double mX = 2.4;
  constexpr double mx = 0.14;
  constexpr double ma = 0.20;
  constexpr double mb = 0.30;
  constexpr double mR1 = 0.82;
  constexpr double mR2 = 0.93;

  gra::LORENTZSCALAR lts = ContinuumVectorCascadeLTSForTest(false);
  const auto R_in_X = TwoBodyRestKinematics(mX, mR1, mR2, 0.74, 0.0);
  lts.decaytree = {ContinuumVectorBranchForTest(
                       "R1", 900001, R_in_X[0], 111, 211, mx, ma, 0.82, 0.0,
                       std::complex<double>(0.73, -0.18), 0.775, 0.08, 2.0),
                   ContinuumVectorBranchForTest(
                       "R2", 900002, R_in_X[1], 111, 321, mx, mb, 2.62, 0.0,
                       std::complex<double>(-0.45, 0.32), 0.890, 0.11, 3.0)};
  lts.amplitude.DECAY_SYM = true;
  lts.decay_symmetry_proposal_active = true;
  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Central;
  SetTransverseProposalForTest(lts);
  lts.central_phase_space_rap_min_cm = -0.6;
  lts.central_phase_space_rap_max_cm = 0.65;

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() == 2);
  REQUIRE(lts.decaytree[0].p4.Rap() > lts.central_phase_space_rap_min_cm);
  REQUIRE(lts.decaytree[0].p4.Rap() < lts.central_phase_space_rap_max_cm);
  REQUIRE(lts.decaytree[1].p4.Rap() > lts.central_phase_space_rap_min_cm);
  REQUIRE(lts.decaytree[1].p4.Rap() < lts.central_phase_space_rap_max_cm);
  REQUIRE(terms[1].tree[0].p4.Rap() < lts.central_phase_space_rap_min_cm);
  REQUIRE(terms[1].tree[1].p4.Rap() < lts.central_phase_space_rap_max_cm);

  const double reference_density =
      gra::math::abs2(gra::spin::CascadeBWProduct(lts.decaytree)) /
      CascadeMassProposalNormForTest(lts.decaytree);
  const double density =
      gra::decay::MixtureDensity(lts, lts.decaytree);
  REQUIRE(density == Approx(0.5 * reference_density).epsilon(1e-12));

  gra::LORENTZSCALAR boosted = lts;
  const double collision_mass = (lts.pbeam1 + lts.pbeam2).M();
  constexpr double beam_rapidity = 1.1;
  const gra::M4Vec boost(0.0, 0.0, collision_mass * std::sinh(beam_rapidity),
                         collision_mass * std::cosh(beam_rapidity));
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pbeam1, 1);
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pbeam2, 1);
  gra::kinematics::LorentzBoost(boost, collision_mass, boosted.pfinal[0], 1);
  for (auto &branch : boosted.decaytree) {
    BoostDecayBranchForTest(branch, boost, collision_mass);
  }

  REQUIRE(boosted.central_phase_space_mode ==
          gra::CentralPhaseSpaceMode::Central);
  REQUIRE(boosted.pbeam1.E() != Approx(boosted.pbeam2.E()).epsilon(1e-8));
  const double boosted_density =
      gra::decay::MixtureDensity(boosted, boosted.decaytree);
  REQUIRE(boosted_density == Approx(density).epsilon(1e-11));
}

TEST_CASE("gra::spin stable leaf mixture enforces fixed intermediate masses",
          "[gra::spin][continuum][proposal]") {
  gra::LORENTZSCALAR lts = ContinuumVectorCascadeLTSForTest(false);
  lts.amplitude.DECAY_SYM = true;
  lts.decay_symmetry_proposal_active = true;
  for (auto &branch : lts.decaytree) {
    branch.mass_proposal = gra::MassProposal::Fixed;
    branch.mass_proposal_min2 = branch.p4.M2();
    branch.mass_proposal_max2 = branch.p4.M2();
  }

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() == 2);
  REQUIRE(terms[1].tree[0].p4.M2() !=
          Approx(lts.decaytree[0].p4.M2()).margin(1e-8));

  const double density =
      gra::decay::MixtureDensity(lts, lts.decaytree);
  const double reference_density =
      gra::math::abs2(gra::spin::CascadeBWProduct(lts.decaytree)) / CascadeMassProposalNormForTest(lts.decaytree);
  REQUIRE(density == Approx(0.5 * reference_density).epsilon(1e-12));
}

// Check that a flat mass proposal leaves the physical cascade poles in the
// weight
TEST_CASE("MProcess full cascade with flat masses retains the physical BW "
          "density",
          "[MProcess][proposal][flat][decay]") {
  gra::LORENTZSCALAR lts = TensorRhoCascadeLTSForTest();
  lts.decay_symmetry_proposal_active = false;
  for (auto &branch : lts.decaytree) {
    branch.mass_proposal = gra::MassProposal::Uniform;
  }

  const double base_phase_space =
      InternalCascadePhaseSpaceForTest(lts.decaytree);
  const double flat_proposal_norm =
      CascadeMassProposalNormForTest(lts.decaytree);
  REQUIRE(flat_proposal_norm != Approx(1.0).epsilon(1e-12));

  ToyHelicityProcess sampler;
  const double cascade_phase_space = sampler.CascadePhaseSpaceForTest(
      lts, {gra::DecayType::Full});
  REQUIRE(cascade_phase_space ==
          Approx(base_phase_space * flat_proposal_norm).epsilon(1e-12));

  const std::complex<double> hard_amplitude(0.61, -0.19);
  const std::complex<double> bw = gra::spin::CascadeBWProduct(lts.decaytree);
  REQUIRE(gra::math::abs2(bw) != Approx(1.0).epsilon(1e-6));
  const std::complex<double> raw_amplitude = hard_amplitude * bw;
  const double sampled_weight =
      cascade_phase_space * gra::math::abs2(raw_amplitude);
  const double expected_weight = base_phase_space * flat_proposal_norm *
                                 gra::math::abs2(hard_amplitude) *
                                 gra::math::abs2(bw);
  REQUIRE(sampled_weight == Approx(expected_weight).epsilon(1e-12));
}

namespace {

// Reconstruct the conditional central mass Jacobian for a crossed history
double CascadeMassJacobianForTest(const gra::LORENTZSCALAR &lts, const std::vector<gra::MDecayBranch> &tree) {
  double threshold = lts.central_phase_space_mass_margin;
  for (const auto &branch : tree) { threshold += branch.p4.M(); }
  return 2.0 * lts.pfinal[0].M2() *
         std::log(lts.central_phase_space_mass_max / std::max(threshold, lts.central_phase_space_mass_cut_min));
}

// Reconstruct the RAMBO direct split from physical daughter momenta
double CascadeRamboWeightForTest(const gra::LORENTZSCALAR &lts, const std::vector<gra::MDecayBranch> &tree) {
  std::vector<gra::M4Vec> momenta;
  for (const auto &branch : tree) {
    momenta.push_back(gra::kinematics::BoostToRestFrame(branch.p4, lts.pfinal[0], "cascade RAMBO test"));
  }
  return gra::kinematics::RamboWeight(lts.pfinal[0].M(), momenta);
}

// Build a massive RAMBO root and construct both decaying daughters through MProcess
gra::LORENTZSCALAR RamboCascadeForTest() {
  ToyHelicityProcess sampler;
  sampler.state.random.SetSeed(80371);
  auto lts = ContinuumVectorCascadeLTSForTest(false);
  lts.PS_active = true;
  lts.pfinal.resize(11);
  for (const int pdg : {900341, 900342}) {
    gra::MDecayBranch spectator;
    spectator.p = ToyParticle("spectator", pdg, 0, 0.3);
    lts.decaytree.push_back(spectator);
  }
  lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, 6.0);
  std::vector<double> masses;
  for (const auto &branch : lts.decaytree) { masses.push_back(branch.p.mass); }
  std::vector<gra::M4Vec> momenta;
  lts.DW = gra::kinematics::RamboMassive(lts.pfinal[0], 6.0, masses, momenta, sampler.state.random);
  REQUIRE(lts.DW.Integral() > 0.0);
  sampler.state.lts.PS_active = true;
  for (const auto &i : indices(lts.decaytree)) {
    auto &branch = lts.decaytree[i];
    branch.p4 = momenta[i];
    REQUIRE(sampler.ConstructDecayKinematics(branch));
    if (branch.legs.empty()) { continue; }
    branch.mass_proposal = gra::MassProposal::Uniform;
    branch.mass_proposal_min2 = 0.0;
    branch.mass_proposal_max2 = 36.0;
    branch.mass_proposal_norm = 36.0;
  }
  lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Factorized;
  lts.central_phase_space_mass_cut_min = 0.1;
  lts.central_phase_space_mass_max = 10.0;
  lts.central_phase_space_mass_margin = 1e-4;
  lts.central_phase_space_generated_jacobian = CascadeMassJacobianForTest(lts, lts.decaytree);
  lts.decay_structure = {gra::DecayType::Full};
  lts.amplitude.DECAY_SYM = true;
  return lts;
}

// Build a third cascade vertex while retaining identical leaves on different branches
gra::LORENTZSCALAR NestedMassCascadeForTest(bool narrow_tail) {
  auto lts = ContinuumVectorCascadeLTSForTest(false);
  lts.PS_active = true;
  auto &branch = lts.decaytree[0].legs[1];
  branch.p.mass = narrow_tail ? 0.01 : 0.18;
  branch.p.width = narrow_tail ? 1e-20 : 0.03;
  const auto momenta = TwoBodyRestKinematics(branch.p4.M(), 0.02, 0.03, 0.49, -0.38);
  branch.legs.resize(2);
  for (const auto &i : indices(branch.legs)) {
    auto &leaf = branch.legs[i];
    leaf.p = ToyParticle("scalar", 900101 + static_cast<int>(i), 0, i == 0 ? 0.02 : 0.03);
    leaf.p.P = i == 0 ? 1 : -1;
    leaf.p4 = BoostFromRestFrame(momenta[i], branch.p4);
  }
  gra::spin::InitTwoBodyBasis(branch.hel, 0.0, 0.0, 0.0, {0.0}, {0.0}, {0.0}, "nested scalar decay");
  branch.hel.T = {{1.0}};
  branch.hel.g_decay = {0.67, -0.23};
  branch.W_event = BranchTwoBodyPhaseSpaceForTest(branch);
  branch.W = gra::kinematics::MCW(branch.W_event);
  return lts;
}

// Set the actual root measure and conditional support used by each central sampler
void SetCascadePhaseSpaceForTest(gra::LORENTZSCALAR &lts, gra::CentralPhaseSpaceMode mode) {
  lts.central_phase_space_mode = mode;
  lts.central_phase_space_rap_min_cm = -10.0;
  lts.central_phase_space_rap_max_cm = 10.0;
  lts.central_phase_space_mass_cut_min = 0.1;
  lts.central_phase_space_mass_max = 10.0;
  lts.central_phase_space_mass_margin = 1e-4;
  lts.central_phase_space_generated_jacobian = CascadeMassJacobianForTest(lts, lts.decaytree);
  if (mode == gra::CentralPhaseSpaceMode::Central) { SetTransverseProposalForTest(lts); }
  lts.DW = mode == gra::CentralPhaseSpaceMode::Factorized
               ? gra::kinematics::MCW(RootCascadePhaseSpaceForTest(lts, lts.decaytree))
               : gra::kinematics::MCW();
}

// Compute the full history density with the F conditional mass map or the C direct root measure
double CentralHistoryDensityForTest(gra::LORENTZSCALAR lts) {
  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  const double phase = StableLeafDensityPhaseSpaceForTest(lts, lts.decaytree);
  double density = 0.0;
  for (const auto &term : terms) {
    double ratio = phase / StableLeafDensityPhaseSpaceForTest(lts, term.tree);
    if (lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Factorized) {
      ratio *= lts.central_phase_space_generated_jacobian / CascadeMassJacobianForTest(lts, term.tree);
    } else if (lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Central) {
      ratio *= lts.central_phase_space_generated_jacobian / CascadeTransverseJacobianForTest(lts, term.tree);
    }
    density += MixedMassDensityForTest(term.tree) * ratio;
  }
  return density / static_cast<double>(terms.size());
}

}  // namespace

// Check the density of the actual F root generator and every inverse leaf assignment
TEST_CASE("four-branch cascades reverse the RAMBO root and every sampled history",
          "[gra::spin][MProcess][cascade][proposal][Rambo]") {
  const auto generated = RamboCascadeForTest();
  const double phase = generated.DW.Integral() * InternalCascadePhaseSpaceForTest(generated.decaytree);
  REQUIRE(CascadeRamboWeightForTest(generated, generated.decaytree) == Approx(generated.DW.Integral()).epsilon(1e-11));
  ToyHelicityProcess sampler;
  sampler.state.lts = generated;
  gra::spin::PrepareStableLeafSymmetryAssignments(sampler.state.lts, generated.decaytree);
  const auto assignments = sampler.state.lts.decay_symmetry_assignments;
  REQUIRE(assignments.size() == 2);
  for (const auto &selected : indices(assignments)) {
    CAPTURE(selected);
    sampler.state.lts = generated;
    sampler.state.lts.decay_symmetry_assignments = assignments;
    sampler.state.lts.decay_symmetry_proposal_index = selected;
    sampler.state.lts.decay_symmetry_proposal_active = true;
    REQUIRE(sampler.ApplyDecaySymmetryProposal());
    auto lts = sampler.state.lts;
    double expected = 0.0;
    const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
    for (const auto &term : terms) {
      const double term_phase = CascadeRamboWeightForTest(lts, term.tree) * InternalCascadePhaseSpaceForTest(term.tree);
      expected += MixedMassDensityForTest(term.tree) * phase / term_phase *
                  generated.central_phase_space_generated_jacobian / CascadeMassJacobianForTest(lts, term.tree);
    }
    expected /= terms.size();
    const double density = gra::decay::MixtureDensity(lts, lts.decaytree);
    REQUIRE(density == Approx(expected).epsilon(1e-10));
    REQUIRE(generated.DW.Integral() * sampler.CascadePhaseSpaceForTest(lts, generated.decay_structure) * density ==
            Approx(phase).epsilon(1e-10));

    auto restored = lts.decaytree;
    REQUIRE(gra::spin::ApplyStableLeafSymmetryAssignment(restored, assignments[selected]));
    const auto check_branch = [&](const auto &self, const gra::MDecayBranch &branch,
                                  const gra::MDecayBranch &original) -> void {
      REQUIRE(gra::math::CheckEMC(branch.p4 - original.p4, 1e-10));
      if (branch.legs.empty()) {
        REQUIRE(branch.p4.M2() == Approx(branch.p.mass * branch.p.mass).margin(1e-10));
        return;
      }
      gra::M4Vec sum;
      for (const auto &i : indices(branch.legs)) {
        sum += branch.legs[i].p4;
        self(self, branch.legs[i], original.legs[i]);
      }
      REQUIRE(gra::math::CheckEMC(branch.p4 - sum, 1e-10));
    };
    for (const auto &i : indices(restored)) { check_branch(check_branch, restored[i], generated.decaytree[i]); }

    auto boosted = lts;
    const gra::M4Vec boost(0.7, -0.4, 1.1, std::sqrt(1.0 + 0.49 + 0.16 + 1.21));
    gra::kinematics::LorentzBoost(boost, 1.0, boosted.pfinal[0], 1);
    for (auto &branch : boosted.decaytree) { BoostDecayBranchForTest(branch, boost, 1.0); }
    REQUIRE(gra::decay::MixtureDensity(boosted, boosted.decaytree) == Approx(density).epsilon(1e-10));
  }
}

// Check physical decay declarations select the coherent proposal
TEST_CASE("coherent cascade sampling follows the physical Jacob-Wick decay declarations",
          "[gra::spin][MProcess][cascade][proposal][declaration]") {
  ToyHelicityProcess sampler;
  auto physical = ContinuumVectorCascadeLTSForTest(true);
  physical.process.root_decay_mode = gra::RootDecayMode::Physical;
  physical.process.ROOT_RES_ACTIVE = true;
  for (const auto &[family, channel] : std::vector<std::pair<std::string, std::string>>{
           {"MP", "RES"}, {"XP", "RES"}, {"GP", "RES"}, {"MP", "CON"}, {"XP", "CON"}, {"GP", "CON"},
           {"MP", "RES+CON"}, {"XP", "RES+CON"}, {"GP", "RES+CON"}, {"gg", "chic(0)"}, {"gg", "chic(1)"}, {"gg", "chic(2)"},
           {"yy", "Higgs"}, {"yy", "monopolium(0)"}}) {
    CAPTURE(family, channel);
    for (const bool spindec : {false, true}) {
      for (const bool coherent : {false, true}) {
        sampler.state.lts = physical;
        sampler.state.lts.process.SPINDEC = spindec;
        sampler.state.lts.amplitude.DECAY_SYM = coherent;
        sampler.ProcPtr.Initialize(family, channel);
        sampler.PrepareDecaySymmetryProposal();
        const auto expected = spindec && coherent ? gra::DecayType::JacobWickCoherent : gra::DecayType::JacobWickIncoherent;
        REQUIRE(sampler.state.lts.decay_structure.type == expected);
        REQUIRE(sampler.state.lts.decay_symmetry_proposal_active == (spindec && coherent));
      }
    }
  }
  for (const std::string family : {"yy", "yy_DZ", "gg"}) {
    for (const bool coherent : {false, true}) {
      sampler.state.lts = physical;
      sampler.state.lts.amplitude.DECAY_SYM = coherent;
      sampler.ProcPtr.Initialize(family, "FLUX");
      sampler.PrepareDecaySymmetryProposal();
      REQUIRE(sampler.state.lts.decay_structure.type == gra::DecayType::None);
      REQUIRE_FALSE(sampler.state.lts.decay_symmetry_proposal_active);
      REQUIRE(sampler.DecaySymmetryCompensationFactor() == Approx(1.0));
    }
  }
  for (const std::string family : {"MP", "XP", "GP", "gg"}) {
    sampler.state.lts = physical;
    sampler.state.lts.decaytree.push_back(physical.decaytree.back().legs.back());
    sampler.state.lts.process.SPINDEC = true;
    sampler.state.lts.amplitude.DECAY_SYM = true;
    sampler.ProcPtr.Initialize(family, family == "gg" ? "chic(0)" : "RES");
    sampler.PrepareDecaySymmetryProposal();
    REQUIRE_FALSE((sampler.state.lts.decay_structure.type == gra::DecayType::JacobWickCoherent));
    REQUIRE_FALSE(sampler.state.lts.decay_symmetry_proposal_active);
  }
}

// Check proposal reversal point by point through three internal decay vertices
TEST_CASE("nested cascades give equal full and Jacob-Wick weights for mixed mass proposals",
          "[gra::spin][MProcess][cascade][proposal][physics]") {
  const auto mode = GENERATE(gra::CentralPhaseSpaceMode::Factorized, gra::CentralPhaseSpaceMode::Central);
  CAPTURE(static_cast<int>(mode));
  ToyHelicityProcess sampler;
  for (const bool narrow_tail : {false, true}) {
    MMatrix<std::complex<double>> reference;
    for (unsigned int selection = 0; selection < 8; ++selection) {
      CAPTURE(narrow_tail, selection);
      auto lts = NestedMassCascadeForTest(narrow_tail);
      SetCascadePhaseSpaceForTest(lts, mode);
      lts.amplitude.DECAY_SYM = false;
      unsigned int mask = selection;
      SetMixedMassProposalsForTest(lts.decaytree, mask);
      const auto jw = gra::spin::ContinuumDecayMatrix(lts, "CM");
      auto complete = lts;
      complete.decay_structure = {gra::DecayType::Full};
      const auto full = gra::spin::ContinuumDecayMatrix(complete, "CM");
      RequireMatrixNear(jw, full, 1e-11);
      if (selection == 0) { reference = full; }
      RequireMatrixNear(full, reference, 1e-10);

      const double density = MixedMassDensityForTest(lts.decaytree);
      const double phase = InternalCascadePhaseSpaceForTest(lts.decaytree);
      const double full_ps = sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full});
      REQUIRE(full_ps * density == Approx(phase).epsilon(1e-11));
      for (const std::string family : {"MP", "XP", "GP"}) {
        CAPTURE(family);
        sampler.ProcPtr.Initialize(family, "CON");
        const auto structure = sampler.ProcPtr.DecayStructureFor(lts);
        REQUIRE(structure.type == gra::DecayType::JacobWickIncoherent);
        const double jw_ps = sampler.CascadePhaseSpaceForTest(lts, structure);
        REQUIRE(jw_ps * MatrixNorm2(jw) == Approx(full_ps * MatrixNorm2(full)).epsilon(1e-11));
        REQUIRE(density * jw_ps * MatrixNorm2(jw) == Approx(phase * MatrixNorm2(reference)).epsilon(1e-11));
      }
    }
  }
}

// Check coherent history reversal in the process phase-space weight
TEST_CASE("nested coherent cascades reverse mixed proposals for either decay construction",
          "[gra::spin][MProcess][cascade][proposal][physics]") {
  const auto mode = GENERATE(gra::CentralPhaseSpaceMode::Factorized, gra::CentralPhaseSpaceMode::Central);
  CAPTURE(static_cast<int>(mode));
  ToyHelicityProcess sampler;
  for (unsigned int selection = 0; selection < 8; ++selection) {
    CAPTURE(selection);
    auto lts = NestedMassCascadeForTest(true);
    SetCascadePhaseSpaceForTest(lts, mode);
    unsigned int mask = selection;
    SetMixedMassProposalsForTest(lts.decaytree, mask);
    lts.amplitude.DECAY_SYM = true;
    lts.decay_symmetry_proposal_active = true;
    const double density = gra::decay::MixtureDensity(lts, lts.decaytree);
    REQUIRE(density == Approx(CentralHistoryDensityForTest(lts)).epsilon(1e-11));

    lts.decay_structure = {gra::DecayType::JacobWickCoherent};
    const auto jw = gra::spin::ContinuumDecayMatrix(lts, "CM");
    lts.decay_structure = {gra::DecayType::Full};
    const auto full = gra::spin::ContinuumDecayMatrix(lts, "CM");
    RequireMatrixNear(jw, full, 1e-10);
    const double full_ps = sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full});
    const double jw_ps = sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::JacobWickCoherent});
    REQUIRE(full_ps * density == Approx(InternalCascadePhaseSpaceForTest(lts.decaytree)).epsilon(1e-11));
    REQUIRE(jw_ps * MatrixNorm2(jw) == Approx(full_ps * MatrixNorm2(full)).epsilon(1e-11));

    auto rotated = lts;
    // Rotate every nested momentum without changing the central invariant mass
    const auto rotate_branch = [&](const auto &self, gra::MDecayBranch &branch) -> void {
      branch.p4.RotateZ(0.73);
      for (auto &leg : branch.legs) { self(self, leg); }
    };
    for (auto &branch : rotated.decaytree) { rotate_branch(rotate_branch, branch); }
    for (auto &p : rotated.pfinal) { p.RotateZ(0.73); }
    rotated.pbeam1.RotateZ(0.73);
    rotated.pbeam2.RotateZ(0.73);
    REQUIRE(gra::decay::MixtureDensity(rotated, rotated.decaytree) ==
            Approx(density).epsilon(1e-10));
    REQUIRE(MatrixNorm2(gra::spin::ContinuumDecayMatrix(rotated, "CM")) ==
            Approx(MatrixNorm2(full)).epsilon(1e-10));
  }
}

// Check that an explicit central root does not introduce a second mass integral
TEST_CASE("single central cascade root preserves the proposal density and phase-space normalization",
          "[gra::spin][MProcess][cascade][proposal][physics]") {
  ToyHelicityProcess sampler;
  for (unsigned int selection = 0; selection < 4; ++selection) {
    CAPTURE(selection);
    auto lts = ContinuumVectorCascadeLTSForTest(false);
    lts.PS_active = true;
    unsigned int mask = selection;
    SetMixedMassProposalsForTest(lts.decaytree, mask);
    const double daughters_ps = sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full});
    gra::MDecayBranch root;
    root.p = ToyParticle("X", 900310, 0, 2.3);
    root.p.width = 0.2;
    root.p4 = lts.pfinal[0];
    root.legs = lts.decaytree;
    root.W_event = BranchTwoBodyPhaseSpaceForTest(root);
    root.W = gra::kinematics::MCW(root.W_event);
    lts.decaytree = {root};
    lts.DW = gra::kinematics::MCW(1.0);
    lts.central_phase_space_mode = gra::CentralPhaseSpaceMode::Factorized;
    lts.central_phase_space_mass_cut_min = 1.8;
    lts.central_phase_space_mass_max = 3.0;
    lts.central_phase_space_mass_margin = 1e-4;
    lts.central_phase_space_generated_jacobian = 2.0 * root.p4.M2() * std::log(3.0 / 1.8);
    REQUIRE(sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full}) ==
            Approx(root.W_event * daughters_ps).epsilon(1e-11));

    lts.amplitude.DECAY_SYM = true;
    lts.decay_symmetry_proposal_active = true;
    const double density = gra::decay::MixtureDensity(lts, lts.decaytree);
    REQUIRE(density == Approx(MixedHistoryDensityForTest(lts)).epsilon(1e-11));
    const double phase = root.W_event * InternalCascadePhaseSpaceForTest(root.legs);
    REQUIRE(sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::Full}) * density ==
            Approx(phase).epsilon(1e-11));
    lts.decaytree[0].p.width *= 2.0;
    REQUIRE(gra::decay::MixtureDensity(lts, lts.decaytree) ==
            Approx(density).epsilon(1e-11));
  }
}

// Keep physical Jacob-Wick matrices independent of every proposal and angular switch
TEST_CASE("complete Jacob-Wick matrices retain all poles independently of the sampler",
          "[gra::spin][cascade][proposal][physics]") {
  for (const bool multibody : {false, true}) {
    for (const bool spin : {false, true}) {
      for (const bool coherent : {false, true}) {
        auto physical = NestedMassCascadeForTest(false);
        if (multibody) {
          gra::MDecayBranch spectator;
          spectator.p = ToyParticle("spectator", 900320, 0, 0.1);
          spectator.p4 = gra::M4Vec(0.0, 0.0, 0.0, spectator.p.mass);
          physical.decaytree.push_back(spectator);
          physical.pfinal[0] += spectator.p4;
        }
        physical.process.SPINDEC = spin;
        physical.amplitude.DECAY_SYM = coherent;
        physical.decay_structure = {gra::DecayType::Full};
        gra::PARAM_RES res;
        res.p = ToyParticle("X", 900310, 0, physical.pfinal[0].M());
        gra::spin::InitTwoBodyBasis(res.hel_decay, 0.0, 1.0, 1.0, {0.0}, {-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0}, "scalar vector decay");
        res.hel_decay.T = MMatrix<std::complex<double>>::IdentityMatrix(3) / std::sqrt(3.0);
        const auto reference = gra::spin::ResonanceDecayMatrix(physical, res, "CM");
        REQUIRE(reference.FrobNorm2() > 0.0);
        for (const bool mixture : {false, true}) {
          for (unsigned int selection = 0; selection < 8; ++selection) {
            CAPTURE(multibody, spin, coherent, mixture, selection);
            auto lts = physical;
            lts.decay_symmetry_proposal_active = mixture;
            unsigned int mask = selection;
            SetMixedMassProposalsForTest(lts.decaytree, mask);
            RequireMatrixNear(gra::spin::ResonanceDecayMatrix(lts, res, "CM"), reference, 1e-10);
          }
        }
      }
    }
  }
}

// Check that the multi-body root fallback retains descendant couplings and propagators
TEST_CASE("multi-body resonance roots preserve the spin-blind cascade normalization",
          "[gra::spin][resonance][cascade][physics]") {
  auto physical = NestedMassCascadeForTest(false);
  gra::MDecayBranch spectator;
  spectator.p = ToyParticle("spectator", 900320, 0, 0.1);
  spectator.p4 = gra::M4Vec(0.0, 0.0, 0.0, spectator.p.mass);
  physical.decaytree.push_back(spectator);
  physical.pfinal[0] += spectator.p4;
  unsigned int mask = 7;
  SetMixedMassProposalsForTest(physical.decaytree, mask);
  physical.amplitude.DECAY_SYM = false;
  gra::PARAM_RES res;
  res.p = ToyParticle("X", 900310, 0, physical.pfinal[0].M());
  const auto fallback = gra::spin::ResonanceDecayMatrix(physical, res, "CM");
  physical.process.SPINDEC = false;
  const auto spin_blind = gra::spin::ResonanceDecayMatrix(physical, res, "CM");
  RequireMatrixNear(fallback, spin_blind, 1e-11);
  REQUIRE(MatrixNorm2(fallback) > 0.0);
}

// Check complex Jacob-Wick amplitudes against an independent coherent contraction
TEST_CASE("physical cascade amplitudes do not depend on the mass proposal", "[gra::spin][cascade][proposal][physics]") {
  ToyHelicityProcess sampler;
  for (const bool coherent : {false, true}) {
    auto reference = ContinuumVectorCascadeLTSForTest(true);
    reference.amplitude.DECAY_SYM = coherent;
    const auto raw = coherent ? CoherentRawContinuumCascadeMatrixForTest(reference, reference.decaytree)
                              : RawContinuumCascadeMatrixForTest(reference, reference.decaytree);
    for (const bool mixture : {false, true}) {
      if (mixture && !coherent) { continue; }
      for (unsigned int selection = 0; selection < 4; ++selection) {
        auto lts = reference;
        unsigned int mask = selection;
        SetMixedMassProposalsForTest(lts.decaytree, mask);
        lts.decay_symmetry_proposal_active = mixture;
        const auto matrix = gra::spin::ContinuumDecayMatrix(lts, "CM");
        RequireMatrixNear(matrix, raw, 1e-11);
        const double density = mixture ? MixedHistoryDensityForTest(lts) : MixedMassDensityForTest(lts.decaytree);
        const double phase = sampler.CascadePhaseSpaceForTest(lts, {gra::DecayType::JacobWickCoherent});
        REQUIRE(phase * density * MatrixNorm2(matrix) ==
                Approx(InternalCascadePhaseSpaceForTest(lts.decaytree) * MatrixNorm2(raw)).epsilon(1e-11));
      }
    }
    // A virtuality supplied by the central mass map still has its physical pole
    for (auto &branch : reference.decaytree) { branch.mass_proposal = gra::MassProposal::None; }
    RequireMatrixNear(gra::spin::ContinuumDecayMatrix(reference, "CM"), raw, 1e-11);
  }
}

// Angular decorrelation retains the complex couplings and all internal propagators
TEST_CASE("spin-blind cascades retain physical couplings and phases", "[gra::spin][cascade][physics]") {
  auto physical = ContinuumVectorCascadeLTSForTest(false);
  physical.process.SPINDEC = false;
  auto unit = physical;
  for (auto &branch : unit.decaytree) { branch.hel.g_decay = 1.0; }
  const auto coupling = physical.decaytree[0].hel.g_decay * physical.decaytree[1].hel.g_decay;
  RequireMatrixNear(gra::spin::ContinuumDecayMatrix(physical, "CM"),
                    gra::spin::ContinuumDecayMatrix(unit, "CM") * coupling, 1e-11);
  gra::PARAM_RES res;
  res.p = ToyParticle("X", 900310, 0, physical.pfinal[0].M());
  const auto scalar = gra::spin::ResonanceDecayMatrix(physical, res, "CM");
  RequireMatrixNear(scalar, gra::spin::ResonanceDecayMatrix(unit, res, "CM") * coupling, 1e-11);
  std::complex<double> expected = coupling;
  for (const auto &branch : physical.decaytree) {
    expected /= std::complex<double>(branch.p4.M2() - gra::math::pow2(branch.p.mass), branch.p.mass * branch.p.width);
  }
  REQUIRE(scalar.size_row() == 1);
  REQUIRE(scalar.size_col() == 1);
  RequireComplexNear(scalar[0][0], expected, 1e-11);
}

TEST_CASE("Complex delta-BW follows the reduced amplitude phase convention",
          "[gra::form][physics]") {
  const double m0 = 3.41475;
  const double width = 0.0108;
  const std::array<double, 3> points = {gra::math::pow2(m0) - m0 * width,
                                        gra::math::pow2(m0),
                                        gra::math::pow2(m0) + m0 * width};

  for (const double m2 : points) {
    const std::complex<double> bw =
        gra::resonance::ComplexBreitWignerAmplitude(m2, m0, width);
    REQUIRE(gra::math::abs2(bw) ==
            Approx(gra::resonance::BreitWignerDensity(m2, m0, width))
                .epsilon(1e-12));
  }

  const std::complex<double> pole = gra::resonance::ComplexBreitWignerAmplitude(
      gra::math::pow2(m0), m0, width);
  REQUIRE(std::abs(pole.real()) <= 1e-12);
  REQUIRE(pole.imag() == Approx(-gra::resonance::BreitWignerAmplitude(
                             gra::math::pow2(m0), m0, width)));

  const std::complex<double> below =
      gra::resonance::ComplexBreitWignerAmplitude(
          gra::math::pow2(m0) - m0 * width, m0, width);
  const std::complex<double> above =
      gra::resonance::ComplexBreitWignerAmplitude(
          gra::math::pow2(m0) + m0 * width, m0, width);
  REQUIRE(below.real() < 0.0);
  REQUIRE(above.real() > 0.0);
  REQUIRE(below.imag() < 0.0);
  REQUIRE(above.imag() < 0.0);
  REQUIRE(below.real() == Approx(-above.real()).epsilon(1e-12));
  REQUIRE(below.imag() == Approx(above.imag()).epsilon(1e-12));

  // Require the common reduced pole phase for each spin prescription
  const auto require_reduced_pole = [m0,
                                     width](const std::complex<double> &line) {
    REQUIRE(std::abs(line.real()) <= 1e-12);
    REQUIRE(line.imag() == Approx(-1.0 / (m0 * width)).epsilon(1e-12));
  };
  require_reduced_pole(
      gra::resonance::FixedWidthLineShape(gra::math::pow2(m0), m0, width));
  require_reduced_pole(
      gra::resonance::KinematicWidthLineShape(gra::math::pow2(m0), m0, width));
  require_reduced_pole(gra::resonance::RunningWidthLineShape(
      gra::math::pow2(m0), m0, width, 1.0));

  // Verify all denominator formulas away from the pole
  const double off_pole = gra::math::pow2(m0) + 0.47;
  const double profile = 0.36;
  gra::PARAM_RES resonance;
  resonance.p.mass = m0;
  resonance.p.width = width;
  resonance.BW = gra::BreitWigner::FixedWidth;
  RequireComplexNear(
      gra::resonance::LineShape(off_pole, resonance),
      1.0 / (off_pole - gra::math::pow2(m0) + gra::math::zi * m0 * width));
  resonance.BW = gra::BreitWigner::KinematicWidth;
  RequireComplexNear(gra::resonance::LineShape(off_pole, resonance),
                     1.0 / (off_pole - gra::math::pow2(m0) +
                            gra::math::zi * std::sqrt(off_pole) * width));
  resonance.BW = gra::BreitWigner::RunningWidth;
  RequireComplexNear(gra::resonance::LineShape(off_pole, resonance, profile),
                     1.0 / (off_pole - gra::math::pow2(m0) +
                            gra::math::zi * m0 * width * profile));

  // A closed channel contributes no absorptive running width
  const auto closed =
      gra::resonance::RunningWidthLineShape(off_pole, m0, width, 0.0);
  REQUIRE(closed.imag() == Approx(0.0).margin(1e-15));
  REQUIRE(closed.real() ==
          Approx(1.0 / (off_pole - gra::math::pow2(m0))).margin(1e-15));
}

TEST_CASE("gra::spin root decay uses its explicit spin frame", "[gra::spin]") {
  for (const auto &frame : std::vector<std::string>{"CM", "HX", "CS"}) {
    CAPTURE(frame);

    gra::LORENTZSCALAR lts =
        MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    gra::MDecayBranch a;
    a.name = "a";
    a.p = ToyParticle("a", 900201, 0, 0.11);
    a.p4 = gra::M4Vec(0.13, -0.07, 0.19, 0.29);

    gra::MDecayBranch b;
    b.name = "b";
    b.p = ToyParticle("b", 900202, 0, 0.17);
    b.p4 = lts.pfinal[0] - a.p4;
    lts.decaytree = {a, b};

    gra::PARAM_RES res;
    res.production_model = gra::ReggeProductionModel::XP;
    lts.process.SPINDEC = true;
    res.p = ToyParticle("X", 900200, 2, lts.pfinal[0].M());
    res.hel_decay =
        SpinOneToScalarScalarHelicityMatrix(std::complex<double>(0.48, -0.27));

    const auto expected = RootDecayReferenceInFrame(lts, res, frame);
    gra::spin::DecayAmp(lts, res, frame);
    RequireMatrixNear(res.decay_f, expected, 1e-11);

    auto override_res = res;
    const auto cm_reference =
        RootDecayReferenceInFrame(lts, override_res, "CM");
    gra::spin::DecayAmp(lts, override_res, "CM");
    RequireMatrixNear(override_res.decay_f, cm_reference, 1e-11);

    if (frame == "HX" || frame == "CS") {
      const auto old_cm_reference = RootDecayReferenceInFrame(lts, res, "CM");
      REQUIRE(MatrixDiffNorm2(expected, old_cm_reference) > 1e-12);
    }
  }
}

TEST_CASE("gra::spin continuum decay uses an explicit root frame",
          "[gra::spin]") {
  const double mX = 2.4;
  const double mR = 1.15;
  const double ms = 0.25;
  const double ma = 0.14;
  const double mb = 0.20;

  gra::LORENTZSCALAR lts;
  lts.process.SPINDEC = true;
  lts.pfinal.resize(3);
  lts.pfinal[0] =
      gra::M4Vec(0.37, -0.22, 0.58,
                 std::sqrt(pow2(mX) + pow2(0.37) + pow2(-0.22) + pow2(0.58)));

  const auto R_s_in_X = TwoBodyRestKinematics(mX, mR, ms, 0.83, -0.41);
  const gra::M4Vec R_lab = BoostFromRestFrame(R_s_in_X[0], lts.pfinal[0]);
  const gra::M4Vec s_lab = BoostFromRestFrame(R_s_in_X[1], lts.pfinal[0]);

  const auto a_b_in_R = TwoBodyRestKinematics(mR, ma, mb, 1.21, 0.74);
  const gra::M4Vec a_lab = BoostFromRestFrame(a_b_in_R[0], R_lab);
  const gra::M4Vec b_lab = BoostFromRestFrame(a_b_in_R[1], R_lab);

  gra::MDecayBranch R;
  R.name = "R";
  R.p = ToyParticle("R", 900001, 2, mR);
  R.p4 = R_lab;
  R.hel =
      SpinOneToScalarScalarHelicityMatrix(std::complex<double>(0.72, -0.31));

  gra::MDecayBranch a;
  a.name = "a";
  a.p = ToyParticle("a", 900011, 0, ma);
  a.p4 = a_lab;

  gra::MDecayBranch b;
  b.name = "b";
  b.p = ToyParticle("b", 900012, 0, mb);
  b.p4 = b_lab;
  R.legs = {a, b};

  gra::MDecayBranch s;
  s.name = "s";
  s.p = ToyParticle("s", 900013, 0, ms);
  s.p4 = s_lab;
  lts.decaytree = {R, s};

  lts.amplitude.BeginCentral();
  const auto hx = gra::spin::ContinuumDecayMatrix(lts, "HX");
  const auto native = gra::spin::ContinuumDecayMatrix(lts, "CM");
  lts.screening.active = true;
  RequireMatrixNear(gra::spin::ContinuumDecayMatrix(lts, "HX"), hx, 1e-13);
  RequireMatrixNear(gra::spin::ContinuumDecayMatrix(lts, "CM"), native, 1e-13);
  lts.screening.active = false;
  lts.amplitude.AbortCentral();

  REQUIRE(MatrixDiffNorm2(hx, native) > 1e-12);
}

TEST_CASE("gra::spin cascaded decay uses consistent helicity frames",
          "[gra::spin]") {
  const double mX = 2.4;
  const double mR = 1.15;
  const double ms = 0.25;
  const double ma = 0.14;
  const double mb = 0.20;

  gra::LORENTZSCALAR lts;
  lts.pfinal.resize(3);
  lts.pfinal[0] =
      gra::M4Vec(0.37, -0.22, 0.58,
                 std::sqrt(pow2(mX) + pow2(0.37) + pow2(-0.22) + pow2(0.58)));

  const auto R_s_in_X = TwoBodyRestKinematics(mX, mR, ms, 0.83, -0.41);
  const gra::M4Vec R_lab = BoostFromRestFrame(R_s_in_X[0], lts.pfinal[0]);
  const gra::M4Vec s_lab = BoostFromRestFrame(R_s_in_X[1], lts.pfinal[0]);

  const auto a_b_in_R = TwoBodyRestKinematics(mR, ma, mb, 1.21, 0.74);
  const gra::M4Vec a_lab = BoostFromRestFrame(a_b_in_R[0], R_lab);
  const gra::M4Vec b_lab = BoostFromRestFrame(a_b_in_R[1], R_lab);

  gra::MDecayBranch R;
  R.name = "R";
  R.p = ToyParticle("R", 900001, 2, mR);
  R.p4 = R_lab;
  R.hel =
      SpinOneToScalarScalarHelicityMatrix(std::complex<double>(0.72, -0.31));

  gra::MDecayBranch a;
  a.name = "a";
  a.p = ToyParticle("a", 900011, 0, ma);
  a.p4 = a_lab;

  gra::MDecayBranch b;
  b.name = "b";
  b.p = ToyParticle("b", 900012, 0, mb);
  b.p4 = b_lab;
  R.legs = {a, b};

  gra::MDecayBranch s;
  s.name = "s";
  s.p = ToyParticle("s", 900013, 0, ms);
  s.p4 = s_lab;

  lts.decaytree = {R, s};
  lts.d0_in_X = BoostToRestFrame(R_lab, lts.pfinal[0]);
  lts.d1_in_X = BoostToRestFrame(s_lab, lts.pfinal[0]);

  gra::PARAM_RES res;
  res.production_model = gra::ReggeProductionModel::XP;
  lts.process.SPINDEC = true;
  res.p = ToyParticle("X", 900000, 2, mX);
  res.hel_decay = SpinOneToSpinOneScalarHelicityMatrix(
      {std::complex<double>(0.43, 0.17), std::complex<double>(-0.28, 0.39),
       std::complex<double>(0.61, -0.22)});

  const auto expected = SequentialOneLevelReference(lts, res);
  gra::spin::DecayAmp(lts, res, "CM");

  RequireMatrixNear(res.decay_f, expected, 1e-11);
}

TEST_CASE("gra::spin deep cascaded decay propagates helicity-frame transforms",
          "[gra::spin]") {
  const double mX = 3.1;
  const double mR = 1.85;
  const double mS = 0.92;
  const double mc = 0.24;
  const double md = 0.18;
  const double ma = 0.13;
  const double mb = 0.19;

  gra::LORENTZSCALAR lts;
  lts.pfinal.resize(3);
  lts.pfinal[0] =
      gra::M4Vec(-0.29, 0.31, 0.73,
                 std::sqrt(pow2(mX) + pow2(-0.29) + pow2(0.31) + pow2(0.73)));

  const auto R_c_in_X = TwoBodyRestKinematics(mX, mR, mc, 0.71, 0.57);
  const gra::M4Vec R_lab = BoostFromRestFrame(R_c_in_X[0], lts.pfinal[0]);
  const gra::M4Vec c_lab = BoostFromRestFrame(R_c_in_X[1], lts.pfinal[0]);

  const auto S_d_in_R = TwoBodyRestKinematics(mR, mS, md, 1.04, -0.88);
  const gra::M4Vec S_lab = BoostFromRestFrame(S_d_in_R[0], R_lab);
  const gra::M4Vec d_lab = BoostFromRestFrame(S_d_in_R[1], R_lab);

  const auto a_b_in_S = TwoBodyRestKinematics(mS, ma, mb, 0.93, 1.37);
  const gra::M4Vec a_lab = BoostFromRestFrame(a_b_in_S[0], S_lab);
  const gra::M4Vec b_lab = BoostFromRestFrame(a_b_in_S[1], S_lab);

  gra::MDecayBranch S;
  S.name = "S";
  S.p = ToyParticle("S", 900101, 2, mS);
  S.p4 = S_lab;
  S.hel =
      SpinOneToScalarScalarHelicityMatrix(std::complex<double>(-0.55, 0.21));

  gra::MDecayBranch a;
  a.name = "a";
  a.p = ToyParticle("a", 900111, 0, ma);
  a.p4 = a_lab;

  gra::MDecayBranch b;
  b.name = "b";
  b.p = ToyParticle("b", 900112, 0, mb);
  b.p4 = b_lab;
  S.legs = {a, b};

  gra::MDecayBranch d;
  d.name = "d";
  d.p = ToyParticle("d", 900113, 0, md);
  d.p4 = d_lab;

  gra::MDecayBranch R;
  R.name = "R";
  R.p = ToyParticle("R", 900102, 2, mR);
  R.p4 = R_lab;
  R.hel = SpinOneToSpinOneScalarHelicityMatrix(
      {std::complex<double>(0.32, -0.44), std::complex<double>(0.25, 0.36),
       std::complex<double>(-0.47, 0.12)});
  R.legs = {S, d};

  gra::MDecayBranch c;
  c.name = "c";
  c.p = ToyParticle("c", 900114, 0, mc);
  c.p4 = c_lab;

  lts.decaytree = {R, c};
  lts.d0_in_X = BoostToRestFrame(R_lab, lts.pfinal[0]);
  lts.d1_in_X = BoostToRestFrame(c_lab, lts.pfinal[0]);

  gra::PARAM_RES res;
  res.production_model = gra::ReggeProductionModel::XP;
  lts.process.SPINDEC = true;
  res.p = ToyParticle("X", 900100, 2, mX);
  res.hel_decay = SpinOneToSpinOneScalarHelicityMatrix(
      {std::complex<double>(0.18, 0.51), std::complex<double>(-0.37, -0.13),
       std::complex<double>(0.46, 0.29)});

  const auto expected = SequentialDeepReference(lts, res);
  gra::spin::DecayAmp(lts, res, "CM");

  RequireMatrixNear(res.decay_f, expected, 1e-11);
}

TEST_CASE("gra::spin balanced cascade preserves ordered unequal spin spaces",
          "[gra::spin][cascade][tensor]") {
  const double mX = 4.0;
  const double mR = 1.25;
  const double mS = 1.65;
  gra::LORENTZSCALAR lts;
  lts.process.SPINDEC = true;
  lts.pfinal.resize(1);
  lts.pfinal[0] =
      gra::M4Vec(0.23, -0.31, 0.47,
                 std::sqrt(pow2(mX) + pow2(0.23) + pow2(-0.31) + pow2(0.47)));

  const auto R_S_in_X = TwoBodyRestKinematics(mX, mR, mS, 0.79, -0.62);
  const gra::M4Vec R_lab = BoostFromRestFrame(R_S_in_X[0], lts.pfinal[0]);
  const gra::M4Vec S_lab = BoostFromRestFrame(R_S_in_X[1], lts.pfinal[0]);
  const auto a_b_in_R = TwoBodyRestKinematics(mR, 0.17, 0.21, 1.08, 0.43);
  const auto c_d_in_S = TwoBodyRestKinematics(mS, 0.19, 0.23, 0.66, -1.17);

  gra::MDecayBranch R;
  R.name = "R";
  R.p = ToyParticle("R", 900201, 2, mR);
  R.p4 = R_lab;
  R.hel = SpinOneToScalarScalarHelicityMatrix({0.49, -0.27});
  R.hel.g_decay = {-0.31, 0.22};
  gra::MDecayBranch a;
  a.name = "a";
  a.p = ToyParticle("a", 900211, 0, 0.17);
  a.p4 = BoostFromRestFrame(a_b_in_R[0], R_lab);
  gra::MDecayBranch b;
  b.name = "b";
  b.p = ToyParticle("b", 900212, 0, 0.21);
  b.p4 = BoostFromRestFrame(a_b_in_R[1], R_lab);
  R.legs = {a, b};

  gra::MDecayBranch S;
  S.name = "S";
  S.p = ToyParticle("S", 900202, 4, mS);
  S.p4 = S_lab;
  gra::spin::InitTwoBodyBasis(S.hel, 2.0, 0.0, 0.0,
                              {-2.0, -1.0, 0.0, 1.0, 2.0}, {0.0}, {0.0},
                              "balanced cascade S decay");
  S.hel.T = MMatrix<std::complex<double>>(1, 1, 0.0);
  S.hel.T[0][0] = {-0.38, 0.41};
  S.hel.g_decay = {0.24, 0.35};
  gra::MDecayBranch c;
  c.name = "c";
  c.p = ToyParticle("c", 900213, 0, 0.19);
  c.p4 = BoostFromRestFrame(c_d_in_S[0], S_lab);
  gra::MDecayBranch d;
  d.name = "d";
  d.p = ToyParticle("d", 900214, 0, 0.23);
  d.p4 = BoostFromRestFrame(c_d_in_S[1], S_lab);
  S.legs = {c, d};
  lts.decaytree = {R, S};

  std::vector<gra::M4Vec> R_daughters = {BoostToRestFrame(a.p4, lts.pfinal[0]),
                                         BoostToRestFrame(b.p4, lts.pfinal[0])};
  std::vector<gra::M4Vec> S_daughters = {BoostToRestFrame(c.p4, lts.pfinal[0]),
                                         BoostToRestFrame(d.p4, lts.pfinal[0])};
  const auto R_in_X = BoostToRestFrame(R.p4, lts.pfinal[0]);
  const auto S_in_X = BoostToRestFrame(S.p4, lts.pfinal[0]);
  gra::kinematics::HXframe(R_daughters, R_in_X);
  gra::kinematics::HXframe(S_daughters, S_in_X);
  const auto first_axis = gra::spin::SpinHalfRotation(R_in_X.Theta(), R_in_X.Phi());
  const MMatrix<std::complex<double>> flip = {{0.0, 1.0}, {1.0, 0.0}};
  const auto inverse_S = gra::spin::SpinHalfRotation(S_in_X.Theta(), S_in_X.Phi()) *
                          gra::spin::SpinHalfRotation(0.0, -gra::math::PI);
  const auto spin_R = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.0, gra::math::PI), 1.0);
  const auto spin_S = gra::spin::SpinRotation(inverse_S.Dagger() * first_axis * flip, 2.0);
  const auto fR = gra::spin::fDecayMatrix(R.hel, R_daughters[0].Theta(),
                                          R_daughters[0].Phi()) * spin_R *
                  R.hel.g_decay;
  const auto fS = gra::spin::fDecayMatrix(S.hel, S_daughters[0].Theta(),
                                          S_daughters[0].Phi()) * spin_S *
                  S.hel.g_decay;

  MMatrix<std::complex<double>> expected(15, 1, 0.0);
  MMatrix<std::complex<double>> reversed(15, 1, 0.0);
  for (std::size_t lambdaR = 0; lambdaR < 3; ++lambdaR) {
    for (std::size_t lambdaS = 0; lambdaS < 5; ++lambdaS) {
      // The ordered tensor maps do not commute: the right helicity is fastest
      expected[5 * lambdaR + lambdaS][0] = fR[0][lambdaR] * fS[0][lambdaS];
      reversed[3 * lambdaS + lambdaR][0] = fS[0][lambdaS] * fR[0][lambdaR];
    }
  }

  const auto actual = gra::spin::ContinuumDecayMatrix(lts, "CM");
  RequireMatrixNear(actual, expected, 1e-11);
  REQUIRE(MatrixDiffNorm2(actual, reversed) > 1e-8);
}

TEST_CASE("gra::spin:: HX central production distinguishes exchange and X axes",
          "[gra::spin]") {
  const gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  gra::PARAM_RES res = MakeToyResonance();

  std::vector<gra::M4Vec> from_q1 = {lts.q1, lts.q2}, from_X = from_q1;
  gra::spin::ProductionFrame(from_q1, lts, lts.process.MP_FRAME, lts.q1);
  gra::spin::ProductionFrame(from_X, lts, lts.process.MP_FRAME, lts.pfinal[0]);
  const auto f_x_from_q1 = gra::spin::fDecayMatrix(res.production[0].hel, from_q1[0].Theta(), from_q1[0].Phi());
  const auto f_x_from_X  = gra::spin::fDecayMatrix(res.production[0].hel, from_X[0].Theta(), from_X[0].Phi());

  double diff2 = 0.0;
  for (std::size_t i = 0; i < f_x_from_X.size_row(); ++i) {
    for (std::size_t j = 0; j < f_x_from_X.size_col(); ++j) {
      diff2 += std::norm(f_x_from_X[i][j] - f_x_from_q1[i][j]);
    }
  }
  REQUIRE(diff2 > 1e-12);
}

TEST_CASE("gra::spin:: FORWARD_VERTEX labels select distinct fixed-spin proton "
          "vertices",
          "[gra::spin]") {
  const gra::LORENTZSCALAR base_lts =
      MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
    CAPTURE(FORWARD_NOFLIP);
    gra::LORENTZSCALAR residue_lts = base_lts;
    residue_lts.process.FORWARD_VERTEX =
        gra::ForwardVertexMode::HelicityResidue;
    residue_lts.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
    gra::LORENTZSCALAR pole_lts = residue_lts;
    pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;
    gra::PARAM_RES residue_res = MakeToyMPResonance(false);
    gra::PARAM_RES pole_res = residue_res;

    const auto residue = gra::spin::Forward(
        residue_lts, residue_res.production[0].tree[0], residue_lts.pbeam1, residue_lts.pfinal[1], false,
        gra::spin::Rows(residue_res.production[0].tree[0], false), residue_lts.process.PHOTON_VERTEX,
        gra::spin::ForwardSpec{gra::ForwardVertexMode::HelicityResidue, 1.0});
    const auto pole =
        gra::spin::Forward(pole_lts, pole_res.production[0].tree[0], pole_lts.pbeam1, pole_lts.pfinal[1], false,
                           gra::spin::Rows(pole_res.production[0].tree[0], false), pole_lts.process.PHOTON_VERTEX,
                           gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});

    REQUIRE(residue.size_row() == pole.size_row());
    REQUIRE(residue.size_col() == pole.size_col());
    REQUIRE(MatrixDiffNorm2(pole, residue) > 1e-12);
    REQUIRE(MatrixNorm2(residue) > 1e-12);
    REQUIRE(MatrixNorm2(pole) > 1e-12);

    const auto residue_prod = gra::rspin::Resonance(residue_lts, residue_res, (residue_res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    const auto pole_prod = gra::rspin::Resonance(pole_lts, pole_res, (pole_res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    REQUIRE(residue_prod.size() == 1);
    REQUIRE(pole_prod.size() == 1);
    REQUIRE(residue_prod[0].size_row() == pole_prod[0].size_row());
    REQUIRE(residue_prod[0].size_col() == pole_prod[0].size_col());
    REQUIRE(residue_prod[0].size_row() == (FORWARD_NOFLIP ? 4 : 16));
    CHECK(MatrixDiffNorm2(residue_prod[0], pole_prod[0]) > 1e-12);
  }
}

TEST_CASE("gra::spin:: unit_residue removes exchange-helicity radial "
          "barrier",
          "[gra::spin]") {
  gra::LORENTZSCALAR residue_lts =
      MakeToyProductionLTSAsymmetric(0.22, 4.35, -4.25);
  residue_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
  gra::LORENTZSCALAR pole_lts = residue_lts;
  pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;
  gra::PARAM_RES res;

  const auto hel = RealisticProtonLegHelicityMatrix(995, 4, 1, 2, 0);
  const auto residue =
      gra::spin::Forward(residue_lts, ForwardBranchForTest(hel), residue_lts.pbeam1, residue_lts.pfinal[1], false,
                         gra::spin::Rows(ForwardBranchForTest(hel), false), residue_lts.process.PHOTON_VERTEX,
                         gra::spin::ForwardSpec{gra::ForwardVertexMode::HelicityResidue, 1.0});
  const auto pole =
      gra::spin::Forward(pole_lts, ForwardBranchForTest(hel), pole_lts.pbeam1, pole_lts.pfinal[1], false,
                         gra::spin::Rows(ForwardBranchForTest(hel), false), pole_lts.process.PHOTON_VERTEX,
                         gra::spin::ForwardSpec{gra::ForwardVertexMode::UnitResidue, 1.0});

  const gra::M4Vec q_up = residue_lts.pbeam1 - residue_lts.pfinal[1];
  const double qt = q_up.Pt();
  const double s0 = 1.0;

  REQUIRE(residue.size_row() == pole.size_row());
  REQUIRE(residue.size_col() == pole.size_col());
  REQUIRE(MatrixDiffNorm2(residue, pole) > 1e-12);

  bool saw_exchange_barrier = false;
  for (std::size_t row = 0; row < hel.lambda_values.size_row(); ++row) {
    for (std::size_t col = 0; col < hel.Jz_values.size(); ++col) {
      CAPTURE(row, col);
      CHECK(std::abs(residue[row][col]) ==
            Approx(ReggeResidueScaleForTest(hel, col, qt, s0)).epsilon(1e-12));
      CHECK(std::abs(pole[row][col]) ==
            Approx(ReggePoleResidueScaleForTest(hel, col)).epsilon(1e-12));
      if (std::abs(IntegerProjectionForTest(hel.Jz_values[col])) > 0) {
        saw_exchange_barrier =
            saw_exchange_barrier || std::abs(std::abs(residue[row][col]) -
                                             std::abs(pole[row][col])) > 1e-12;
      }
    }
  }
  REQUIRE(saw_exchange_barrier);
}

TEST_CASE("gra::spin:: FORWARD_VERTEX numerical split across exchange quantum "
          "numbers and kinematics",
          "[gra::spin]") {
  struct ExchangeCase {
    const char *label;
    int pdg;
    int spinX2;
    int parity;
    std::size_t proton_l;
    std::size_t proton_two_s;
  };

  const std::array<ExchangeCase, 4> exchanges = {
      {{"J=0 P=+", 991, 0, 1, 0, 0},
       {"J=0 P=-", 989, 0, -1, 1, 2},
       {"J=1 P=-", 993, 2, -1, 1, 2},
       {"J=2 P=+", 995, 4, 1, 2, 0}}};
  const std::array<gra::LORENTZSCALAR, 2> events = {
      {MakeToyProductionLTSAsymmetric(0.22, 4.35, -4.25),
       MakeToyProductionLTSDphi(0.85, 0.31)}};

  gra::PARAM_RES regge_res;

  double largest_relative_diff = 0.0;
  for (std::size_t event_index = 0; event_index < events.size();
       ++event_index) {
    const auto &event_lts = events[event_index];

    for (const auto &exchange : exchanges) {
      for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
        gra::MParticle exchange_particle;
        exchange_particle.name = exchange.label;
        exchange_particle.pdg = exchange.pdg;
        exchange_particle.spinX2 = exchange.spinX2;
        exchange_particle.P = exchange.parity;
        exchange_particle.C = 1;

        gra::PARAM_RES regge = regge_res;
        gra::LORENTZSCALAR regge_lts = event_lts;
        regge_lts.process.FORWARD_VERTEX =
            gra::ForwardVertexMode::HelicityResidue;
        regge_lts.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
        gra::LORENTZSCALAR pole_lts = regge_lts;
        pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;
        regge.production_model = gra::ReggeProductionModel::MP;
        regge.p = ToyParticle(9000000 + exchange.pdg, 0, 1, 1, "toy_scalar");
        gra::MDecayBranch up;
        up.p = exchange_particle;
        up.hel = RealisticProtonLegHelicityMatrix(
            exchange.pdg, exchange.spinX2, exchange.parity, exchange.proton_l,
            exchange.proton_two_s);
        gra::MDecayBranch dn = up;
        regge.production = {{{up, dn}}};

        gra::HELMatrix central;
        central.BR = 1.0;
        central.P_symmetry = true;
        central.C_symmetry = true;
        central.alpha_ls.Set(0, 0, 1.0);
        gra::spin::InitTMatrix(central, regge.p, exchange_particle,
                               exchange_particle, true,
                               "test forward vertex central", false, false);
        regge.production.front().hel = central;
        PrepareToyPoleOperators(regge, gra::ReggeProductionModel::MP,
                                {{{0, 0, 1.0}}});

        gra::PARAM_RES pole = regge;

        const auto regge_prod = gra::rspin::Resonance(regge_lts, regge, (regge.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
        const auto pole_prod =
            gra::rspin::Resonance(pole_lts, pole, (pole.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
        REQUIRE(regge_prod.size() == 1);
        REQUIRE(pole_prod.size() == 1);
        const double regge_norm2 = MatrixNorm2(regge_prod[0]);
        const double pole_norm2 = MatrixNorm2(pole_prod[0]);
        const double diff2 = MatrixDiffNorm2(regge_prod[0], pole_prod[0]);
        const double scale = std::max({regge_norm2, pole_norm2, 1e-300});
        const double relative_diff = diff2 / scale;
        largest_relative_diff = std::max(largest_relative_diff, relative_diff);

        REQUIRE(std::isfinite(regge_norm2));
        REQUIRE(std::isfinite(pole_norm2));
        REQUIRE(std::isfinite(diff2));
        REQUIRE(scale > 1e-16);
        REQUIRE(regge_prod[0].size_row() == (FORWARD_NOFLIP ? 4 : 16));
        if (exchange.spinX2 == 0) {
          CHECK(diff2 == Approx(0.0).margin(1e-12));
        } else {
          CHECK(diff2 > 1e-12);
        }
      }
    }
  }
  CHECK(largest_relative_diff > 1e-12);
}

TEST_CASE(
    "gra::spin:: FORWARD_VERTEX modes produce distinct full helicity matrices",
    "[gra::spin]") {
  gra::LORENTZSCALAR base_lts = MakeToyProductionLTSDphi(0.85, 0.31);
  gra::PARAM_RES regge = MakeRealisticTensorMPResonance();

  for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
    CAPTURE(FORWARD_NOFLIP);
    gra::LORENTZSCALAR regge_lts = base_lts;
    regge_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
    regge_lts.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
    gra::LORENTZSCALAR pole_lts = regge_lts;
    pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;
    gra::PARAM_RES pole = regge;

    const auto regge_channels =
        gra::rspin::Resonance(regge_lts, regge, (regge.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    const auto pole_channels =
        gra::rspin::Resonance(pole_lts, pole, (pole.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    REQUIRE(regge_channels.size() == 1);
    REQUIRE(pole_channels.size() == 1);
    const auto &r = regge_channels[0];
    const auto &c = pole_channels[0];
    REQUIRE(r.size_row() == c.size_row());
    REQUIRE(r.size_col() == c.size_col());

    const auto up_rows =
        FORWARD_NOFLIP
            ? HelicityConservingRows(regge.production[0].tree[0].hel)
            : AllHelicityRowsForTest(regge.production[0].tree[0].hel);
    const auto dn_rows =
        FORWARD_NOFLIP
            ? HelicityConservingRows(regge.production[0].tree[1].hel)
            : AllHelicityRowsForTest(regge.production[0].tree[1].hel);
    REQUIRE(r.size_row() == up_rows.size() * dn_rows.size());
    REQUIRE(r.size_col() == regge.production[0].hel.Jz_values.size());

    double max_relative_diff = 0.0;
    for (std::size_t row = 0; row < r.size_row(); ++row) {
      const std::size_t up_index = row / dn_rows.size();
      const std::size_t dn_index = row % dn_rows.size();
      REQUIRE(up_index < up_rows.size());
      REQUIRE(dn_index < dn_rows.size());

      for (std::size_t col = 0; col < r.size_col(); ++col) {
        const double abs_diff = std::abs(r[row][col] - c[row][col]);
        const double scale =
            std::max({std::abs(r[row][col]), std::abs(c[row][col]), 1e-300});
        const double relative_diff = abs_diff / scale;
        max_relative_diff = std::max(max_relative_diff, relative_diff);
      }
    }

    CHECK(max_relative_diff > 1e-12);
  }
}

TEST_CASE("gra::spin:: FORWARD_VERTEX modes produce distinct card-loaded "
          "helicity matrices",
          "[gra::spin]") {
  ToyHelicityProcess proc;
  proc.SetProcessForTest("XP", "RES");
  proc.state.lts.PDG = LoadedPDGTable();
  proc.state.lts.beam1 = proc.state.lts.PDG.FindByPDG(2212);
  proc.state.lts.beam2 = proc.state.lts.PDG.FindByPDG(2212);
  proc.SetDecayMode("pi+ pi-");

  gra::PARAM_RES f2 =
      gra::resonance::Read("RES/f2_1270.json", proc.state.random, gra::ReggeProductionModel::XP);
  proc.SetResonances({{"f2_1270", f2}});
  REQUIRE_NOTHROW(proc.InitializeProcessAmplitude());

  const auto resonances = proc.GetResonances();
  REQUIRE(resonances.count("f2_1270") == 1);
  gra::PARAM_RES regge = resonances.at("f2_1270");

  gra::LORENTZSCALAR lts = MakeToyProductionLTSDphi(0.85, 0.31);
  lts.process.MP_FRAME = "HX";
  lts.PDG = LoadedPDGTable();

  for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
    CAPTURE(FORWARD_NOFLIP);
    gra::LORENTZSCALAR regge_lts = lts;
    regge_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
    regge_lts.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
    gra::LORENTZSCALAR pole_lts = regge_lts;
    pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;
    gra::PARAM_RES pole = regge;

    const auto regge_channels =
        gra::rspin::Resonance(regge_lts, regge, (regge.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    const auto pole_channels =
        gra::rspin::Resonance(pole_lts, pole, (pole.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
    REQUIRE(regge_channels.size() == pole_channels.size());

    for (std::size_t channel = 0; channel < regge_channels.size(); ++channel) {
      CAPTURE(channel);
      const auto &r = regge_channels[channel];
      const auto &c = pole_channels[channel];
      REQUIRE(r.size_row() == c.size_row());
      REQUIRE(r.size_col() == c.size_col());
      REQUIRE(regge.production[channel].tree.size() == 2);

      const auto up_rows =
          FORWARD_NOFLIP
              ? HelicityConservingRows(regge.production[channel].tree[0].hel)
              : AllHelicityRowsForTest(regge.production[channel].tree[0].hel);
      const auto dn_rows =
          FORWARD_NOFLIP
              ? HelicityConservingRows(regge.production[channel].tree[1].hel)
              : AllHelicityRowsForTest(regge.production[channel].tree[1].hel);
      REQUIRE(r.size_row() == up_rows.size() * dn_rows.size());
      REQUIRE(r.size_col() == regge.production[channel].hel.Jz_values.size());

      double max_relative_diff = 0.0;
      for (std::size_t row = 0; row < r.size_row(); ++row) {
        const std::size_t up_index = row / dn_rows.size();
        const std::size_t dn_index = row % dn_rows.size();
        REQUIRE(up_index < up_rows.size());
        REQUIRE(dn_index < dn_rows.size());
        for (std::size_t col = 0; col < r.size_col(); ++col) {
          const double abs_diff = std::abs(r[row][col] - c[row][col]);
          const double scale =
              std::max({std::abs(r[row][col]), std::abs(c[row][col]), 1e-300});
          const double relative_diff = abs_diff / scale;
          max_relative_diff = std::max(max_relative_diff, relative_diff);
        }
      }

      CHECK(max_relative_diff > 1e-12);
    }
  }
}

TEST_CASE("gra::rspin:: MP Prod3 matches manual "
          "no-flip projection",
          "[gra::rspin]") {
  for (const auto &frame : std::vector<std::string>{"HX", "CS", "CM"}) {
    CAPTURE(frame);

    gra::LORENTZSCALAR lts =
        MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
    lts.process.MP_FRAME = frame;
    gra::PARAM_RES res = MakeToyMPResonance(false);

    const auto up_rows = HelicityConservingRows(res.production[0].tree[0].hel);
    const auto dn_rows = HelicityConservingRows(res.production[0].tree[1].hel);

    const auto f_up_first_full  = gra::spin::Forward(lts, res.production[0].tree[0], lts.pbeam1, lts.pfinal[1], false,
                                                     gra::spin::Rows(res.production[0].tree[0], false),
                                                     lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
    const auto f_dn_second_full = gra::spin::Forward(lts, res.production[0].tree[1], lts.pbeam2, lts.pfinal[2], true,
                                                     gra::spin::Rows(res.production[0].tree[1], false),
                                                     lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
    const auto f_x = gra::mpom::Fusion(lts, res.production[0].pole.value());

    const auto f_in = SelectRows(f_up_first_full, up_rows)
                          .Kronecker(SelectRows(f_dn_second_full, dn_rows));
    const auto expected = CanonicalCompactPairSpinRowsForTest(f_in * f_x);
    const auto actual =
        gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);

    REQUIRE(actual.size() == 1);
    REQUIRE(actual[0].size_row() == 4);
    RequireMatrixNear(actual[0], expected);
  }
}

TEST_CASE(
    "gra::spin:: coherent photon EPA flux follows the standard elastic formula",
    "[gra::spin][photon]") {
  const gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  const double reference1 =
      CoherentPhotonFluxReference(lts.x1, lts.t1, lts.qt1);
  const double reference2 =
      CoherentPhotonFluxReference(lts.x2, lts.t2, lts.qt2);

  REQUIRE(reference1 > 0.0);
  REQUIRE(reference2 > 0.0);
  REQUIRE(gra::flux::CohFlux(lts.x1, lts.t1, lts.qt1) ==
          Approx(reference1).epsilon(1e-12));
  REQUIRE(gra::flux::CohFlux(lts.x2, lts.t2, lts.qt2) ==
          Approx(reference2).epsilon(1e-12));
}

// Check explicit Dirac and chiral gamma representations and their identities
TEST_CASE("Dirac and chiral gamma matrices realize the same Clifford algebra",
          "[gra::spin][dirac][algebra]") {
  using Complex = std::complex<double>;
  using Matrix = gra::MMatrix<Complex>;
  const gra::MDirac dirac("DIRAC");
  const gra::MDirac chiral("CHIRAL");
  const Matrix identity(4, 4, "eye");
  const Matrix zero(4, 4, 0.0);
  const Matrix gamma0_dirac{{1.0, 0.0, 0.0, 0.0},
                            {0.0, 1.0, 0.0, 0.0},
                            {0.0, 0.0, -1.0, 0.0},
                            {0.0, 0.0, 0.0, -1.0}};
  const Matrix gamma0_chiral{{0.0, 0.0, 1.0, 0.0},
                             {0.0, 0.0, 0.0, 1.0},
                             {1.0, 0.0, 0.0, 0.0},
                             {0.0, 1.0, 0.0, 0.0}};
  const Matrix gamma1{{0.0, 0.0, 0.0, 1.0},
                      {0.0, 0.0, 1.0, 0.0},
                      {0.0, -1.0, 0.0, 0.0},
                      {-1.0, 0.0, 0.0, 0.0}};
  const Matrix gamma2{{0.0, 0.0, 0.0, -gra::math::zi},
                      {0.0, 0.0, gra::math::zi, 0.0},
                      {0.0, gra::math::zi, 0.0, 0.0},
                      {-gra::math::zi, 0.0, 0.0, 0.0}};
  const Matrix gamma3{{0.0, 0.0, 1.0, 0.0},
                      {0.0, 0.0, 0.0, -1.0},
                      {-1.0, 0.0, 0.0, 0.0},
                      {0.0, 1.0, 0.0, 0.0}};
  const Matrix gamma5_dirac = gamma0_chiral;
  const Matrix gamma5_chiral{{-1.0, 0.0, 0.0, 0.0},
                             {0.0, -1.0, 0.0, 0.0},
                             {0.0, 0.0, 1.0, 0.0},
                             {0.0, 0.0, 0.0, 1.0}};
  const std::array<Matrix, 5> expected_dirac = {
      gamma0_dirac, gamma1, gamma2, gamma3, gamma5_dirac};
  const std::array<Matrix, 5> expected_chiral = {
      gamma0_chiral, gamma1, gamma2, gamma3, gamma5_chiral};

  for (const auto &mu : gra::aux::indices(expected_dirac)) {
    REQUIRE(dirac.gamma_up[mu].IsApprox(expected_dirac[mu], 1.0e-14));
    REQUIRE(chiral.gamma_up[mu].IsApprox(expected_chiral[mu], 1.0e-14));
    REQUIRE((dirac.S_basis * chiral.gamma_up[mu] * dirac.S_basis.Dagger())
                .IsApprox(dirac.gamma_up[mu], 1.0e-14));
  }
  for (std::size_t mu = 0; mu < 4; ++mu) {
    for (std::size_t nu = 0; nu < 4; ++nu) {
      const Matrix target = identity * (2.0 * dirac.g[mu][nu]);
      REQUIRE((dirac.gamma_up[mu] * dirac.gamma_up[nu] +
               dirac.gamma_up[nu] * dirac.gamma_up[mu])
                  .IsApprox(target, 1.0e-14));
      REQUIRE((chiral.gamma_up[mu] * chiral.gamma_up[nu] +
               chiral.gamma_up[nu] * chiral.gamma_up[mu])
                  .IsApprox(target, 1.0e-14));
    }
  }
  for (const gra::MDirac *algebra : {&dirac, &chiral}) {
    const Matrix gamma5 = gra::math::zi * algebra->gamma_up[0] *
                          algebra->gamma_up[1] * algebra->gamma_up[2] *
                          algebra->gamma_up[3];
    REQUIRE(gamma5.IsApprox(algebra->gamma_up[4], 1.0e-14));
    REQUIRE((algebra->PR() + algebra->PL()).IsApprox(identity, 1.0e-14));
    REQUIRE((algebra->PR() * algebra->PL()).IsApprox(zero, 1.0e-14));
    REQUIRE((algebra->C_up() * algebra->C_up().Dagger())
                .IsApprox(identity, 1.0e-14));
    REQUIRE(algebra->C_up().Transpose().IsApprox(-algebra->C_up(), 1.0e-14));
    for (std::size_t mu = 0; mu < 4; ++mu) {
      REQUIRE((algebra->C_up() * algebra->gamma_lo[mu].Transpose() *
               algebra->C_up().Dagger())
                  .IsApprox(-algebra->gamma_lo[mu], 1.0e-14));
    }
  }
}

TEST_CASE("Dirac currents transform covariantly between screening frames",
          "[gra::spin][photon][covariance]") {
  const gra::MDirac::Current current = {
      std::complex<double>(1.2, -0.4), std::complex<double>(-0.7, 0.3),
      std::complex<double>(0.2, 0.9), std::complex<double>(1.1, -0.6)};
  const gra::M4Vec q(0.3, -0.8, 1.4, 2.6);
  const auto contract = [](const gra::MDirac::Current &j,
                           const gra::M4Vec &p) {
    return j[0] * p.E() + j[1] * p.Px() + j[2] * p.Py() +
           j[3] * p.Pz();
  };

  const double mass = 3.1;
  const gra::M4Vec boost(0.7, -0.4, 1.2,
                         std::sqrt(mass * mass + 2.09));
  gra::M4Vec boosted_q = q;
  gra::kinematics::LorentzBoost(boost, mass, boosted_q, -1);
  const auto boosted_current =
      gra::MDirac::BoostCurrent(current, boost, mass, -1);
  RequireComplexNear(contract(boosted_current, boosted_q),
                     contract(current, q), 2.0e-13);

  const auto rotation = boosted_q.RotationTo({0.0, 0.0, 1.0});
  gra::M4Vec rotated_q = boosted_q;
  rotated_q.Rotate(rotation);
  const auto rotated_current =
      gra::MDirac::RotateCurrent(boosted_current, rotation);
  RequireComplexNear(contract(rotated_current, rotated_q),
                     contract(current, q), 2.0e-13);
}

TEST_CASE("gra::spin:: QED photon current is conserved and EPA-normalized",
          "[gra::spin][photon]") {
  const gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  const std::vector<std::pair<double, double>> transitions = {
      {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};
  gra::MDirac dirac("DIRAC");

  for (const auto &[h_in, h_out] : transitions) {
    const auto current = dirac.DiracPauliCurrent(
        lts.pbeam1, lts.pfinal[1], -lts.q1,
        gra::MDirac::FermionKind::Particle, h_in, h_out, gra::PDG::mp,
        gra::form::F1(lts.t1), gra::form::F2(lts.t1));
    std::complex<double> ward = 0.0;
    for (std::size_t mu = 0; mu < 4; ++mu) {
      ward += lts.q1[mu] * current[mu];
    }
    CAPTURE(h_in, h_out, ward);
    REQUIRE(std::abs(ward) < 1e-7);
  }

  const auto qed = gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions,
                                                           m_values, "QED");
  const auto epa = gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions,
                                                           m_values, "EPA") *
                   gra::math::msqrt(gra::qed::ElasticPhotonDensity(lts, 1));
  const double target = gra::qed::ElasticPhotonDensity(lts, 1);

  REQUIRE(target > 0.0);
  REQUIRE(gra::spin::SourceSpinAveragedDensity(qed, 2, "test QED source") ==
          Approx(target).epsilon(1e-10));
  REQUIRE(gra::spin::SourceSpinAveragedDensity(epa, 2, "test EPA source") ==
          Approx(target).epsilon(1e-12));

  std::complex<double> overlap = 0.0;
  for (std::size_t row : {0U, 3U}) {
    for (const auto &col : gra::aux::indices(m_values)) {
      overlap += std::conj(epa[row][col]) * qed[row][col];
    }
  }
  CAPTURE(qed[0][0], qed[0][2], qed[3][0], qed[3][2], epa[0][0], epa[0][2],
          epa[3][0], epa[3][2], overlap);
  REQUIRE(overlap.real() > 0.0);
  REQUIRE(std::abs(overlap.imag()) <= 1.0e-10 * std::abs(overlap));

  const std::vector<std::pair<double, double>> noflip_transitions = {
      {-0.5, -0.5}, {0.5, 0.5}};
  const auto qed_noflip = gra::qed::PhotonSourceMatrixTransitions(
      lts, 1, noflip_transitions, m_values, "QED");
  for (std::size_t col = 0; col < m_values.size(); ++col) {
    REQUIRE(std::abs(qed_noflip[0][col] - qed[0][col]) < 1e-12);
    REQUIRE(std::abs(qed_noflip[1][col] - qed[3][col]) < 1e-12);
  }

  double qed_flip_norm = 0.0;
  double epa_flip_norm = 0.0;
  for (std::size_t row = 0; row < transitions.size(); ++row) {
    for (std::size_t col = 0; col < m_values.size(); ++col) {
      if (m_values[col] == 0) {
        REQUIRE(std::abs(qed[row][col]) < 1e-14);
        REQUIRE(std::abs(epa[row][col]) < 1e-14);
      }
      if (std::abs(transitions[row].first - transitions[row].second) > 1e-12) {
        qed_flip_norm += std::norm(qed[row][col]);
        epa_flip_norm += std::norm(epa[row][col]);
      }
    }
  }
  REQUIRE(qed_flip_norm > 1e-14);
  REQUIRE(epa_flip_norm < 1e-20);
}

TEST_CASE("gra::spin photon contractions retain transverse state averages",
          "[gra::spin][photon][normalization]") {
  const gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  const std::vector<std::pair<double, double>> transitions = {{-0.5, -0.5},
                                                              {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 1};
  const auto upper = gra::qed::PhotonSourceMatrixTransitions(
      lts, 1, transitions, m_values, "EPA");
  const auto lower = gra::qed::PhotonSourceMatrixTransitions(
      lts, 2, transitions, m_values, "EPA", true);
  CHECK(gra::spin::SourceSpinAveragedDensity(upper, 2,
                                             "test upper transverse average") ==
        Approx(1.0).epsilon(1.0e-13));
  CHECK(gra::spin::SourceSpinAveragedDensity(lower, 2,
                                             "test lower transverse average") ==
        Approx(1.0).epsilon(1.0e-13));

  MMatrix<std::complex<double>> scalar_source(2, 1, 1.0);
  MMatrix<std::complex<double>> gamma_scalar_central(2, 2, 0.0);
  for (std::size_t i = 0; i < 2; ++i) {
    gamma_scalar_central[i][i] = 1.0 / std::sqrt(2.0);
  }
  CHECK(gamma_scalar_central.FrobNorm2() == Approx(1.0).epsilon(1.0e-13));
  const auto gamma_scalar =
      gra::spin::Contract(upper, scalar_source, gamma_scalar_central);
  CHECK(gra::spin::SourceSpinAveragedDensity(gamma_scalar, 4,
                                             "test gamma scalar average") ==
        Approx(0.5).epsilon(1.0e-13));

  MMatrix<std::complex<double>> gamma_gamma_central(4, 4, 0.0);
  for (std::size_t i = 0; i < 4; ++i) {
    gamma_gamma_central[i][i] = 0.5;
  }
  CHECK(gamma_gamma_central.FrobNorm2() == Approx(1.0).epsilon(1.0e-13));
  const auto gamma_gamma =
      gra::spin::Contract(upper, lower, gamma_gamma_central);
  CHECK(gra::spin::SourceSpinAveragedDensity(gamma_gamma, 4,
                                             "test gamma gamma average") ==
        Approx(0.25).epsilon(1.0e-13));
}

TEST_CASE("gra::spin QED sources carry every external helicity phase",
          "[gra::spin][photon][helicity][phase]") {
  const std::vector<std::pair<double, double>> transitions = {
      {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};
  const double rotation = -0.64;
  const gra::LORENTZSCALAR proton_event = MakeToyQEDLTS();
  gra::MParticle antiproton = proton_event.PDG.FindByPDG(-gra::PDG::PDG_p);

  for (const int leg : {1, 2}) {
    for (const bool antiparticle : {false, true}) {
      CAPTURE(leg, antiparticle);
      gra::LORENTZSCALAR event = proton_event;
      if (antiparticle) {
        if (leg == 1) {
          event.beam1 = antiproton;
        } else {
          event.beam2 = antiproton;
        }
      }
      const auto reference = gra::qed::PhotonSourceMatrixTransitions(
          event, leg, transitions, m_values, "QED", leg == 2);
      const gra::LORENTZSCALAR rotated_event =
          RotateToyEventAroundZ(event, rotation);
      const auto rotated = gra::qed::PhotonSourceMatrixTransitions(
          rotated_event, leg, transitions, m_values, "QED", leg == 2);

      for (std::size_t row = 0; row < transitions.size(); ++row) {
        const int delta = static_cast<int>(
            std::llround(transitions[row].first - transitions[row].second));
        for (std::size_t col = 0; col < m_values.size(); ++col) {
          const int harmonic = (leg == 2 ? -1 : 1) * (delta + m_values[col]);
          const std::complex<double> expected =
              reference[row][col] *
              std::exp(gra::math::zi * static_cast<double>(harmonic) *
                       rotation);
          CAPTURE(row, col, harmonic, reference[row][col], rotated[row][col]);
          REQUIRE(std::abs(rotated[row][col] - expected) <
                  4.0e-10 * std::max(1.0, std::abs(expected)));
        }
      }
    }
  }
}

TEST_CASE("gra::spin:: source density uses the configured incoming spin count",
          "[gra::spin][normalization]") {
  MMatrix<std::complex<double>> source(3, 1, 1.0);
  REQUIRE(gra::spin::SourceSpinAveragedDensity(
              source, 3, "three-state source") == Approx(1.0).epsilon(1e-14));
  REQUIRE_THROWS_AS(
      gra::spin::SourceSpinAveragedDensity(source, 0, "invalid source"),
      std::invalid_argument);
}

TEST_CASE("gra::spin:: validated rest frames reject invalid configurations",
          "[gra::spin][frame][configuration]") {
  const gra::M4Vec particle(0.1, -0.2, 0.3, 1.0);
  const gra::M4Vec spacelike_system(1.0, 0.0, 0.0, 0.5);
  const gra::M4Vec negative_energy_system(0.0, 0.0, 0.0, -1.0);
  REQUIRE_THROWS_AS(gra::kinematics::BoostToRestFrame(
                        particle, spacelike_system, "spacelike test system"),
                    gra::AmplitudeFailure);
  REQUIRE_THROWS_AS(
      gra::kinematics::BoostToRestFrame(particle, negative_energy_system,
                                        "negative-energy test system"),
      gra::AmplitudeFailure);

  std::vector<gra::M4Vec> zero_total = {gra::M4Vec(1.0, 0.0, 0.0, 1.0),
                                        gra::M4Vec(-1.0, 0.0, 0.0, -1.0)};
  REQUIRE_THROWS_AS(
      gra::kinematics::BoostToRestFrame(zero_total, "zero-total test system"),
      gra::AmplitudeFailure);
}

TEST_CASE("gra::spin:: random density matrices support doubled spin labels",
          "[gra::spin][density]") {
  gra::MRandom rng;
  rng.SetSeed(1337);

  for (const int spinX2 : {0, 1, 2, 3, 4}) {
    const auto rho = gra::spin::RandomRho(spinX2, false, rng);
    CAPTURE(spinX2);
    REQUIRE(rho.size_row() == static_cast<std::size_t>(spinX2 + 1));
    REQUIRE(rho.size_col() == static_cast<std::size_t>(spinX2 + 1));
    REQUIRE(gra::spin::Positivity(rho, spinX2 / 2.0));
  }
  REQUIRE_THROWS_AS(gra::spin::RandomRho(-1, false, rng),
                    std::invalid_argument);

  constexpr int spinX2 = 3;
  const auto parity_rho = gra::spin::RandomRho(spinX2, true, rng);
  const std::size_t n = parity_rho.size_row();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      const auto parity_image = parity_rho[n - 1 - i][n - 1 - j];
      REQUIRE(std::abs(parity_rho[i][j] - parity_image) < 1e-12);
    }
  }
}

TEST_CASE("gra::spin:: elastic QED sources use beam particle metadata",
          "[gra::spin][photon][electron]") {
  gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  constexpr double electron_mass = 0.00051099895;
  lts.beam1.name = "e-";
  lts.beam1.pdg = 11;
  lts.beam1.chargeX3 = -3;
  lts.beam1.spinX2 = 1;
  lts.beam1.mass = electron_mass;
  lts.pbeam1.SetPxPyPzM(0.0, 0.0, 5.0, electron_mass);
  lts.pfinal[1].SetPxPyPzM(0.24, 0.11, 4.55, electron_mass);
  UpdateToyDerivedKinematics(lts);

  REQUIRE_NOTHROW(gra::qed::ValidateEmitter(lts, 1, "electron source test"));
  REQUIRE(gra::qed::Charge(lts, 1) == Approx(-1.0));
  const auto electron_emitter =
      gra::qed::Emitter(lts.beam1, lts.t1, "electron particle test");
  REQUIRE(electron_emitter.fermion_kind == gra::MDirac::FermionKind::Particle);

  gra::MDirac dirac("DIRAC");
  for (const double helicity : {-0.5, 0.5}) {
    const auto current = dirac.DiracPauliCurrent(
        lts.pbeam1, lts.pfinal[1], -lts.q1,
        gra::MDirac::FermionKind::Particle, helicity, helicity, electron_mass,
        1.0, 0.0);
    std::complex<double> ward = 0.0;
    for (std::size_t mu = 0; mu < 4; ++mu) {
      ward += lts.q1[mu] * current[mu];
    }
    CAPTURE(helicity, ward);
    REQUIRE(std::abs(ward) < 1e-7);
  }

  const std::vector<std::pair<double, double>> transitions = {
      {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};
  const auto electron_source = gra::qed::PhotonSourceMatrixTransitions(
      lts, 1, transitions, m_values, "QED");
  const double target = gra::qed::ElasticPhotonDensity(lts, 1);
  REQUIRE(target > 0.0);
  REQUIRE(gra::spin::SourceSpinAveragedDensity(electron_source, 2,
                                               "electron QED source") ==
          Approx(target).epsilon(1e-10));

  gra::LORENTZSCALAR positron = lts;
  positron.beam1.name = "e+";
  positron.beam1.pdg = -11;
  positron.beam1.chargeX3 = 3;
  const auto positron_emitter =
      gra::qed::Emitter(positron.beam1, positron.t1, "positron particle test");
  REQUIRE(positron_emitter.fermion_kind ==
          gra::MDirac::FermionKind::Antiparticle);
  REQUIRE(positron_emitter.charge == Approx(1.0));

  const auto positron_source = gra::qed::PhotonSourceMatrixTransitions(
      positron, 1, transitions, m_values, "QED");
  for (std::size_t row = 0; row < transitions.size(); ++row) {
    const double spin_phase =
        4.0 * transitions[row].first * transitions[row].second;
    for (std::size_t col = 0; col < m_values.size(); ++col) {
      const std::complex<double> expected =
          -spin_phase * electron_source[row][col];
      CAPTURE(row, col, electron_source[row][col], positron_source[row][col]);
      REQUIRE(std::abs(positron_source[row][col] - expected) < 2e-10);
    }
  }

  const auto proton = lts.PDG.FindByPDG(gra::PDG::PDG_p);
  const auto antiproton = lts.PDG.FindByPDG(-gra::PDG::PDG_p);
  const auto proton_emitter =
      gra::qed::Emitter(proton, -0.2, "proton particle test");
  const auto antiproton_emitter =
      gra::qed::Emitter(antiproton, -0.2, "antiproton particle test");
  REQUIRE(proton_emitter.fermion_kind == gra::MDirac::FermionKind::Particle);
  REQUIRE(antiproton_emitter.fermion_kind ==
          gra::MDirac::FermionKind::Antiparticle);
  REQUIRE(antiproton_emitter.charge == Approx(-1.0));
  REQUIRE(antiproton_emitter.F1 == Approx(proton_emitter.F1));
  REQUIRE(antiproton_emitter.F2 == Approx(proton_emitter.F2));
  REQUIRE(antiproton_emitter.charge * antiproton_emitter.F2 ==
          Approx(-proton_emitter.charge * proton_emitter.F2));

  gra::LORENTZSCALAR neutral = lts;
  neutral.beam1.pdg = gra::PDG::PDG_n;
  neutral.beam1.chargeX3 = 0;
  neutral.beam1.mass = gra::PDG::mp;
  REQUIRE_THROWS_AS(
      gra::qed::ValidateEmitter(neutral, 1, "neutral source test"),
      std::invalid_argument);

  gra::LORENTZSCALAR charged_scalar = lts;
  charged_scalar.beam1.pdg = gra::PDG::PDG_pip;
  charged_scalar.beam1.chargeX3 = 3;
  charged_scalar.beam1.spinX2 = 0;
  charged_scalar.beam1.mass = gra::PDG::mpi;
  REQUIRE_THROWS_AS(gra::qed::ValidateEmitter(charged_scalar, 1,
                                              "charged scalar source test"),
                    std::invalid_argument);
}

// Check all row-major current contractions against their explicit currents
TEST_CASE("gra::qed:: elastic photon exchange fills all pair helicities",
          "[gra::qed][photon][helicity]") {
  const auto beams = ProtonInitialState();
  const double momentum = 12.0;
  const double abs_t = 0.03;
  const double azimuth = 0.41;
  const double mass = beams[0].mass;
  const double energy = std::sqrt(momentum * momentum + mass * mass);
  const double cosine = 1.0 - abs_t / (2.0 * momentum * momentum);
  const double sine = std::sqrt(1.0 - cosine * cosine);
  const gra::M4Vec p1_in(0.0, 0.0, momentum, energy);
  const gra::M4Vec p2_in(0.0, 0.0, -momentum, energy);
  const gra::M4Vec p1_out(momentum * sine * std::cos(azimuth),
                          momentum * sine * std::sin(azimuth),
                          momentum * cosine, energy);
  const gra::M4Vec p2_out = p1_in + p2_in - p1_out;
  gra::MDirac dirac("DIRAC");
  const auto amplitude = gra::qed::ElasticSpinHalfPhotonExchange(
      dirac, beams[0], beams[1], p1_in, p2_in, p1_out, p2_out);
  const auto emitter1 =
      gra::qed::Emitter(beams[0], -abs_t, "pair-helicity test beam 1");
  const auto emitter2 =
      gra::qed::Emitter(beams[1], -abs_t, "pair-helicity test beam 2");
  constexpr auto helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  const double coefficient =
      4.0 * gra::math::PI * gra::qed::alpha_0 / (p1_out - p1_in).M2();

  for (std::size_t row = 0; row < 4; ++row) {
    const std::size_t out1 = row / 2;
    const std::size_t out2 = row % 2;
    for (std::size_t col = 0; col < 4; ++col) {
      const std::size_t in1 = col / 2;
      const std::size_t in2 = col % 2;
      const double lambda_in1 = helicity_x2[in1] / 2.0;
      const double lambda_in2 = helicity_x2[in2] / 2.0;
      const double lambda_out1 = helicity_x2[out1] / 2.0;
      const double lambda_out2 = helicity_x2[out2] / 2.0;
      const auto current1 = dirac.DiracPauliCurrent(
          p1_in, p1_out, p1_out - p1_in, emitter1.fermion_kind, lambda_in1,
          lambda_out1, emitter1.mass, emitter1.F1, emitter1.F2);
      const auto current2 = dirac.DiracPauliCurrent(
          p2_in, p2_out, p2_out - p2_in, emitter2.fermion_kind, lambda_in2,
          lambda_out2, emitter2.mass, emitter2.F1, emitter2.F2);
      const auto expected =
          gra::MDirac::ElasticSpinHalfColliderCurrentPhase(
              emitter1.fermion_kind, 1, lambda_in1, lambda_out1, azimuth) *
          gra::MDirac::ElasticSpinHalfColliderCurrentPhase(
              emitter2.fermion_kind, 2, lambda_in2, lambda_out2, azimuth) *
          coefficient * dirac.g.BilinearForm(current1, current2);
      const std::size_t index = gra::spin::PairHelicityMatrixIndex(row, col);
      CAPTURE(row, col, amplitude[index], expected);
      REQUIRE(std::abs(amplitude[index] - expected) <
              2.0e-13 * std::max(1.0, std::abs(expected)));
    }
  }
}

// Check the collider helicity section for every proton charge conjugation
TEST_CASE("gra::qed:: elastic photon collider section is covariant",
          "[gra::qed][photon][helicity][covariance]") {
  const auto protons = ProtonInitialState();
  gra::MParticle antiproton = protons[0];
  antiproton.name = "pbar";
  antiproton.pdg = -gra::PDG::PDG_p;
  antiproton.chargeX3 = -3;
  const std::vector<std::vector<gra::MParticle>> beam_states = {
      {protons[0], protons[1]},
      {protons[0], antiproton},
      {antiproton, protons[1]},
      {antiproton, antiproton}};
  constexpr auto helicity_x2 = gra::spin::BinaryHelicityLabelsX2();
  const double momentum = 18.0;
  const double abs_t = 0.07;
  const double energy =
      std::sqrt(momentum * momentum + gra::PDG::mp * gra::PDG::mp);
  const gra::MDirac dirac("DIRAC");

  // Construct one ordered elastic pair at the requested transfer azimuth
  const auto elastic_momenta = [momentum, energy](const double transfer,
                                                  const double azimuth) {
    const double cosine = 1.0 - transfer / (2.0 * momentum * momentum);
    const double sine = std::sqrt(1.0 - cosine * cosine);
    const gra::M4Vec p1_in(0.0, 0.0, momentum, energy);
    const gra::M4Vec p2_in(0.0, 0.0, -momentum, energy);
    const gra::M4Vec p1_out(momentum * sine * std::cos(azimuth),
                            momentum * sine * std::sin(azimuth),
                            momentum * cosine, energy);
    return std::array<gra::M4Vec, 4>{p1_in, p2_in, p1_out,
                                     p1_in + p2_in - p1_out};
  };

  // Evaluate the public exact current-current amplitude
  const auto photon_matrix = [&dirac](
                                 const std::vector<gra::MParticle> &state,
                                 const std::array<gra::M4Vec, 4> &momenta) {
    return gra::qed::ElasticSpinHalfPhotonExchange(dirac, state[0], state[1],
                                                   momenta[0], momenta[1],
                                                   momenta[2], momenta[3]);
  };

  for (const auto &state : beam_states) {
    CAPTURE(state[0].pdg, state[1].pdg);
    const auto reference = photon_matrix(state, elastic_momenta(abs_t, 0.0));
    for (const double azimuth : {-1.13, 0.41, 2.07}) {
      const auto rotated =
          photon_matrix(state, elastic_momenta(abs_t, azimuth));
      for (std::size_t row = 0; row < 4; ++row) {
        const std::size_t out1 = row / 2;
        const std::size_t out2 = row % 2;
        for (std::size_t col = 0; col < 4; ++col) {
          const std::size_t in1 = col / 2;
          const std::size_t in2 = col % 2;
          const int harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(
              helicity_x2[in1], helicity_x2[in2], helicity_x2[out1],
              helicity_x2[out2]);
          const std::size_t index =
              gra::spin::PairHelicityMatrixIndex(row, col);
          const std::complex<double> expected =
              reference[index] *
              std::exp(gra::math::zi * static_cast<double>(harmonic) * azimuth);
          CAPTURE(row, col, azimuth, harmonic, rotated[index], expected);
          REQUIRE(std::abs(rotated[index] - expected) <
                  2.0e-11 * std::max(1.0, std::abs(expected)));
        }
      }
    }

    for (std::size_t row = 0; row < 4; ++row) {
      const std::size_t out1 = row / 2;
      const std::size_t out2 = row % 2;
      for (std::size_t col = 0; col < 4; ++col) {
        const std::size_t in1 = col / 2;
        const std::size_t in2 = col % 2;
        const int harmonic = gra::spin::ColliderSpinHalfHelicityHarmonic(
            helicity_x2[in1], helicity_x2[in2], helicity_x2[out1],
            helicity_x2[out2]);
        const double reciprocity_sign =
            gra::spin::ColliderSpinHalfReciprocitySign(
                helicity_x2[in1], helicity_x2[in2], helicity_x2[out1],
                helicity_x2[out2]);
        const std::size_t index = gra::spin::PairHelicityMatrixIndex(row, col);
        const std::size_t transpose =
            gra::spin::PairHelicityMatrixIndex(col, row);
        const std::complex<double> expected =
            reciprocity_sign * reference[transpose];
        CAPTURE(row, col, harmonic, reference[index], expected);
        REQUIRE(std::abs(reference[index] - expected) <
                2.0e-11 * std::max(1.0, std::abs(expected)));
      }
    }

    const auto forward = photon_matrix(state, elastic_momenta(1.0e-8, 0.37));
    const auto forward_momenta = elastic_momenta(1.0e-8, 0.37);
    const double t = (forward_momenta[2] - forward_momenta[0]).M2();
    const double s = (forward_momenta[0] + forward_momenta[1]).M2();
    const double charge_product =
        static_cast<double>(state[0].chargeX3 * state[1].chargeX3) / 9.0;
    const double point_amplitude = 8.0 * gra::math::PI * gra::qed::alpha_0 *
                                   charge_product *
                                   (s - 2.0 * gra::PDG::mp * gra::PDG::mp) / t;
    for (std::size_t row = 0; row < 4; ++row) {
      const std::complex<double> diagonal = forward[5 * row] / point_amplitude;
      CAPTURE(row, diagonal);
      REQUIRE(std::abs(diagonal - 1.0) < 2.0e-6);
    }
  }

  const auto reference_momenta = elastic_momenta(abs_t, 0.0);
  const auto ppbar = photon_matrix(beam_states[1], reference_momenta);
  const auto pbarp = photon_matrix(beam_states[2], reference_momenta);
  for (std::size_t row = 0; row < 4; ++row) {
    const std::size_t out1 = row / 2;
    const std::size_t out2 = row % 2;
    const std::size_t swapped_row =
        gra::spin::BinaryPairHelicityIndex(out2, out1);
    for (std::size_t col = 0; col < 4; ++col) {
      const std::size_t in1 = col / 2;
      const std::size_t in2 = col % 2;
      const std::size_t swapped_col =
          gra::spin::BinaryPairHelicityIndex(in2, in1);
      const double sign = gra::spin::ColliderSpinHalfReciprocitySign(
          helicity_x2[in1], helicity_x2[in2], helicity_x2[out1],
          helicity_x2[out2]);
      const std::size_t index = gra::spin::PairHelicityMatrixIndex(row, col);
      const std::complex<double> expected =
          sign *
          pbarp[gra::spin::PairHelicityMatrixIndex(swapped_row, swapped_col)];
      CAPTURE(row, col, ppbar[index], expected);
      REQUIRE(std::abs(ppbar[index] - expected) <
              2.0e-11 * std::max(1.0, std::abs(expected)));
    }
  }

}

// Check photon emitter validation for direct and collinear production
TEST_CASE("Physical photon processes validate EPA and QED beam emitters",
          "[gra::spin][photon][initialization][validation]") {
  ModelParamRestoreGuard restore;
  for (const std::string mode : {"EPA", "QED"}) {
    CAPTURE(mode);
    const auto tune = WriteModifiedPhotoVMTune(
        "photon_emitter_" + mode, [&mode](auto &document) {
          for (auto &value : document.at("PARAM_REGGE").at("PHOTON_VERTEX")) { value = mode; }
        });

    const auto model = gra::MModelTune::Load(tune.second);
    const std::vector<std::array<std::string, 3>> direct_processes = {
        {"yy", "QED", "mu+ mu-"}, {"ygg", "Z", "mu+ mu-"}};

    // Configure complete gamma-p kinematics when the direct process uses a UGD
    const auto configure_direct_process =
        [](ToyHelicityProcess &process, const std::string &family,
           const std::string &channel, const std::string &decay) {
          if (family == "ygg") {
            process.state.lts = MakeToyPhotoZFFbar(13);
          }
          ConfigureToyProductionProcess(process, family, channel, decay);
        };

    for (const auto &[family, channel, decay] : direct_processes) {
      CAPTURE(family, channel);
      ToyHelicityProcess proton_process;
      configure_direct_process(proton_process, family, channel, decay);
      proton_process.SetModelTune(model);
      gra::MODELPARAM = tune.first;
      REQUIRE_NOTHROW(proton_process.InitializeProcessAmplitude());

      for (const int leg : {1, 2}) {
        CAPTURE(leg);
        ToyHelicityProcess neutral_process;
        configure_direct_process(neutral_process, family, channel, decay);
        const auto neutron =
            neutral_process.state.lts.PDG.FindByPDG(gra::PDG::PDG_n);
        if (leg == 1) {
          neutral_process.state.lts.beam1 = neutron;
        } else {
          neutral_process.state.lts.beam2 = neutron;
        }
        neutral_process.SetModelTune(model);
        gra::MODELPARAM = tune.first;
        REQUIRE_THROWS(neutral_process.InitializeProcessAmplitude());
      }
    }

    for (const std::string family : {"yy_DZ", "yy_LUX"}) {
      CAPTURE(family);
      ToyHelicityProcess collinear_process;
      ConfigureToyProductionProcess(collinear_process, family, "EPA",
                                    "mu+ mu-");
      if (family == "yy_LUX") {
        collinear_process.state.lts.LHAPDFSET =
            "LUXqed17_plus_PDF4LHC15_nnlo_100";
      }
      collinear_process.state.lts.beam1 =
          collinear_process.state.lts.PDG.FindByPDG(gra::PDG::PDG_n);
      collinear_process.SetModelTune(model);
      gra::MODELPARAM = tune.first;
      REQUIRE_NOTHROW(collinear_process.InitializeProcessAmplitude());
    }
  }
}

// Select independent spin conventions from one tune for each production model
TEST_CASE("Regge setup selects model specific spin steering", "[gra::spin][initialization][models]") {
  ModelParamRestoreGuard restore;
  const auto tune = WriteModifiedPhotoVMTune("regge_model_spin", [](auto &document) {
    auto &regge = document.at("PARAM_REGGE");
    regge.at("DECAY_BARRIERS") = {{"MP", false}, {"XP", true}, {"GP", false}};
    regge.at("FORWARD_VERTEX") = {{"MP", "unit_residue"}, {"XP", "helicity_residue"}, {"GP", "unit_residue"}};
    regge.at("PHOTON_VERTEX") = {{"MP", "EPA"}, {"XP", "QED"}, {"GP", "EPA"}};
    regge.at("TU_SIGN") = {{"MP", "positive"}, {"XP", "auto"}, {"GP", "positive"}};
  });
  const auto tune_model = gra::MModelTune::Load(tune.second);
  for (const std::string model : {"MP", "XP", "GP"}) {
    CAPTURE(model);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, model, "CON", "pi+ pi-");
    process.SetModelTune(tune_model);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto &state = process.state.lts.process;
    CHECK(state.DECAY_BARRIER == (model == "XP"));
    CHECK(state.FORWARD_VERTEX == (model == "XP" ? gra::ForwardVertexMode::HelicityResidue : gra::ForwardVertexMode::UnitResidue));
    CHECK(state.PHOTON_VERTEX == (model == "XP" ? "QED" : "EPA"));
    CHECK(state.TU_SIGN == (model == "XP" ? "auto" : "positive"));
  }
}

// Check that Regge spin steering is required only by MP, XP and GP
TEST_CASE("PARAM_REGGE spin steering is scoped to Regge production models",
          "[gra::spin][initialization][validation]") {
  ModelParamRestoreGuard restore;
  const auto tune =
      WriteModifiedPhotoVMTune("regge_spin_scope", [](auto &document) {
        for (const std::string key :
             {"DECAY_BARRIERS", "DERIVATIVE_FACTOR", "FORWARD_VERTEX",
              "PHOTON_VERTEX", "TU_SIGN"}) {
          document.at("PARAM_REGGE").erase(key);
        }
      });
  const auto model = gra::MModelTune::Load(tune.second);
  gra::MODELPARAM = tune.first;

  for (const auto &[family, channel, decay] :
       std::vector<std::array<std::string, 3>>{{"yy", "QED", "mu+ mu-"},
                                               {"TP", "CON", "pi+ pi-"}}) {
    CAPTURE(family, channel);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, channel, decay);
    process.SetModelTune(model);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
  }

  ToyHelicityProcess regge;
  ConfigureToyProductionProcess(regge, "MP", "CON", "pi+ pi-");
  regge.SetModelTune(model);
  REQUIRE_THROWS(regge.InitializeProcessAmplitude());
}

// Check that MP frame card and command steering are accepted only by MP
TEST_CASE("MP_FRAME steering is scoped to MP processes",
          "[gra::spin][initialization][validation]") {
  ModelParamRestoreGuard restore;
  const auto missing_tune =
      WriteModifiedPhotoVMTune("mp_frame_missing", [](auto &document) {
        document.at("PARAM_REGGE").erase("MP_FRAME");
      });
  const auto missing_model = gra::MModelTune::Load(missing_tune.second);
  gra::MODELPARAM = missing_tune.first;

  for (const auto &[family, channel, decay] :
       std::vector<std::array<std::string, 3>>{{"XP", "CON", "pi+ pi-"},
                                               {"GP", "CON", "pi+ pi-"},
                                               {"TP", "CON", "pi+ pi-"},
                                               {"yy", "QED", "mu+ mu-"}}) {
    CAPTURE(family, channel);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, channel, decay);
    process.SetModelTune(missing_model);
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
  }

  ToyHelicityProcess missing_mp;
  ConfigureToyProductionProcess(missing_mp, "MP", "CON", "pi+ pi-");
  missing_mp.SetModelTune(missing_model);
  REQUIRE_THROWS(missing_mp.InitializeProcessAmplitude());

  gra::MODELPARAM = modelfile;
  const auto model = gra::MModelTune::Load(modelfile);
  ToyHelicityProcess mp;
  ConfigureToyProductionProcess(mp, "MP", "CON", "pi+ pi-");
  mp.SetModelTune(model);
  mp.SetMPFrame("HX");
  REQUIRE_NOTHROW(mp.InitializeProcessAmplitude());
  REQUIRE(mp.state.lts.process.MP_FRAME == "HX");

  for (const auto &[family, channel, decay] :
       std::vector<std::array<std::string, 3>>{{"XP", "CON", "pi+ pi-"},
                                               {"GP", "CON", "pi+ pi-"},
                                               {"TP", "CON", "pi+ pi-"},
                                               {"yy", "QED", "mu+ mu-"}}) {
    CAPTURE(family, channel);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, channel, decay);
    process.SetModelTune(model);
    process.SetMPFrame("HX");
    REQUIRE_THROWS(process.InitializeProcessAmplitude());
  }

  ToyHelicityProcess invalid_mp;
  ConfigureToyProductionProcess(invalid_mp, "MP", "CON", "pi+ pi-");
  invalid_mp.SetModelTune(model);
  invalid_mp.SetMPFrame("GJ");
  REQUIRE_THROWS(invalid_mp.InitializeProcessAmplitude());
}

TEST_CASE("gra::spin:: QED and EPA photon columns rotate consistently",
          "[gra::spin][photon][covariance]") {
  const gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  constexpr double angle = 0.37;
  const gra::LORENTZSCALAR rotated = RotateToyEventAroundZ(lts, angle);
  const std::vector<std::pair<double, double>> transitions = {{-0.5, -0.5},
                                                              {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};

  for (const auto &[leg, second_exchange_daughter] :
       std::vector<std::pair<int, bool>>{{1, false}, {2, true}}) {
    const auto qed = gra::qed::PhotonSourceMatrixTransitions(
        lts, leg, transitions, m_values, "QED", second_exchange_daughter);
    const auto rotated_qed = gra::qed::PhotonSourceMatrixTransitions(
        rotated, leg, transitions, m_values, "QED", second_exchange_daughter);
    const auto epa = gra::qed::PhotonSourceMatrixTransitions(
        lts, leg, transitions, m_values, "EPA", second_exchange_daughter);
    const auto rotated_epa = gra::qed::PhotonSourceMatrixTransitions(
        rotated, leg, transitions, m_values, "EPA", second_exchange_daughter);

    for (std::size_t row = 0; row < transitions.size(); ++row) {
      CAPTURE(leg, row);
      REQUIRE(std::abs(qed[row][0]) > 1e-14);
      REQUIRE(std::abs(qed[row][2]) > 1e-14);
      const std::complex<double> qed_relative_rotation =
          rotated_qed[row][2] * qed[row][0] /
          (qed[row][2] * rotated_qed[row][0]);
      const std::complex<double> epa_relative_rotation =
          rotated_epa[row][2] * epa[row][0] /
          (epa[row][2] * rotated_epa[row][0]);
      REQUIRE(std::abs(qed_relative_rotation - epa_relative_rotation) < 1e-10);
    }
  }
}

// Check the electric photon phase in both oriented beam bases
TEST_CASE(
    "gra::spin:: EPA photon source carries exchange-helicity azimuth phases",
    "[gra::spin][photon]") {
  const gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  const std::vector<std::pair<double, double>> transitions = {{-0.5, -0.5},
                                                              {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};
  const auto upper = gra::qed::PhotonSourceMatrixTransitions(
      lts, 1, transitions, m_values, "EPA", false);
  const auto lower = gra::qed::PhotonSourceMatrixTransitions(
      lts, 2, transitions, m_values, "EPA", true);
  const std::complex<double> norm = 1.0 / std::sqrt(2.0);

  for (std::size_t row = 0; row < transitions.size(); ++row) {
    REQUIRE(std::abs(upper[row][1]) < 1e-14);
    REQUIRE(std::abs(lower[row][1]) < 1e-14);
    for (std::size_t col : {0U, 2U}) {
      const int m = m_values[col];
      const std::complex<double> expected_upper =
          -static_cast<double>(m) * norm * std::exp(gra::math::zi * static_cast<double>(m) * lts.q1.Phi());
      const std::complex<double> expected_lower =
          static_cast<double>(m) * norm * std::exp(-gra::math::zi * static_cast<double>(m) * lts.q2.Phi());
      REQUIRE(std::abs(upper[row][col] - expected_upper) < 1e-14);
      REQUIRE(std::abs(lower[row][col] - expected_lower) < 1e-14);
    }
  }
}

TEST_CASE("gra::spin:: inclusive QED photon source falls back to EPA rows",
          "[gra::spin][photon]") {
  gra::LORENTZSCALAR lts = MakeToyQEDLTS();
  lts.excite1 = true;
  const std::vector<std::pair<double, double>> transitions = {
      {-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}};
  const std::vector<int> m_values = {-1, 0, 1};

  const auto qed = gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions,
                                                           m_values, "QED");
  const auto epa = gra::qed::PhotonSourceMatrixTransitions(lts, 1, transitions,
                                                           m_values, "EPA");
  RequireMatrixNear(qed, epa);
}

TEST_CASE("gra::spin:: FORWARD_NOFLIP=false preserves forward proton flip rows "
          "in X-Pomeron 2to3 production",
          "[gra::spin]") {
  gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  gra::PARAM_RES res = MakeToyCovariantXPonance();
  lts.process.FORWARD_NOFLIP = false;

  const auto f_up_first_full  = gra::spin::Forward(lts, res.production[0].tree[0], lts.pbeam1, lts.pfinal[1], false,
                                                   gra::spin::Rows(res.production[0].tree[0], false),
                                                   lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
  const auto f_dn_second_full = gra::spin::Forward(lts, res.production[0].tree[1], lts.pbeam2, lts.pfinal[2], true,
                                                   gra::spin::Rows(res.production[0].tree[1], false),
                                                   lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{});
  const auto f_x              = gra::xpom::Fusion(lts, res.production[0].pole.value());

  const auto leg_major_expected =
      f_up_first_full.Kronecker(f_dn_second_full) * f_x;
  const auto expected = CanonicalFullPairSpinRowsForTest(leg_major_expected);
  const auto actual =
      gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);

  REQUIRE(actual.size() == 1);
  REQUIRE(actual[0].size_row() == 16);
  RequireMatrixNear(actual[0], expected);
}

// Check the beam exchange norm after summing the complete produced spin basis
TEST_CASE("gra::spin:: unpolarized MP HX production is beam-exchange invariant",
          "[gra::spin]") {
  gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  lts.process.MP_FRAME = "HX";
  const gra::PARAM_RES res = MakeRealisticTensorMPResonance();

  const gra::LORENTZSCALAR mirrored = BeamExchangeMirror(lts);
  const auto channels =
      gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  const auto mirrored_channels =
      gra::rspin::Resonance(mirrored, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  REQUIRE(channels.size() == 1);
  REQUIRE(mirrored_channels.size() == 1);
  REQUIRE(channels.front().size_row() == mirrored_channels.front().size_row());
  REQUIRE(channels.front().size_col() == mirrored_channels.front().size_col());

  const double amp2 = MatrixNorm2(channels.front()) /
                      static_cast<double>(channels.front().size_row());
  const double mirrored_amp2 =
      MatrixNorm2(mirrored_channels.front()) /
      static_cast<double>(mirrored_channels.front().size_row());
  CAPTURE(amp2, mirrored_amp2);
  REQUIRE(amp2 > 0.0);
  REQUIRE(mirrored_amp2 == Approx(amp2).epsilon(2.0e-12));
}

TEST_CASE("gra::spin:: fused subchannel contraction matches the indexed sum",
          "[gra::spin][subchannel][tensor]") {
  using Complex = std::complex<double>;
  gra::HELMatrix upper;
  gra::HELMatrix lower;
  upper.lambda_values =
      MMatrix<double>({{-0.5, -0.5}, {-0.5, 0.5}, {0.5, -0.5}, {0.5, 0.5}});
  lower.lambda_values = upper.lambda_values;

  MMatrix<Complex> upper_frame(4, 2, 0.0);
  MMatrix<Complex> lower_frame(4, 3, 0.0);
  for (std::size_t i = 0; i < upper_frame.size_row(); ++i) {
    for (std::size_t a = 0; a < upper_frame.size_col(); ++a) {
      upper_frame[i][a] = Complex(0.13 * static_cast<double>(1 + i + 2 * a),
                                  -0.07 * static_cast<double>(1 + 2 * i + a));
    }
    for (std::size_t b = 0; b < lower_frame.size_col(); ++b) {
      lower_frame[i][b] = Complex(-0.09 * static_cast<double>(1 + i + b),
                                  0.11 * static_cast<double>(1 + i + 3 * b));
    }
  }

  const auto left = ToyParticle("f1", 900400, 1, 0.2);
  const auto right = ToyParticle("f2", 900401, 1, 0.3);
  for (const bool analytic : {false, true}) {
    lower.domain = analytic ? gra::HelicityDomain::ReggeTrajectory
                            : gra::HelicityDomain::PhysicalPole;
    for (const bool swap : {false, true}) {
      CAPTURE(analytic, swap);
      const auto actual = gra::spin::Subchannel({upper, upper_frame}, {lower, lower_frame}, left, right, swap);
      MMatrix<Complex> expected(6, 4, 0.0);
      for (std::size_t i = 0; i < 4; ++i) {
        for (std::size_t j = 0; j < 4; ++j) {
          if (std::abs(upper.lambda_values[i][1] + lower.lambda_values[j][1]) >
              1e-14) {
            continue;
          }
          const std::size_t upper_final =
              static_cast<std::size_t>(upper.lambda_values[i][0] + 0.5);
          const std::size_t lower_final =
              static_cast<std::size_t>(lower.lambda_values[j][0] + 0.5);
          const std::size_t column = swap ? 2 * lower_final + upper_final
                                          : 2 * upper_final + lower_final;
          for (std::size_t a = 0; a < 2; ++a) {
            for (std::size_t b = 0; b < 3; ++b) {
              expected[3 * a + b][column] +=
                  upper_frame[i][a] * lower_frame[j][b];
            }
          }
        }
      }
      RequireMatrixNear(actual, expected, 1e-14);
    }
  }

  REQUIRE_THROWS_AS(gra::spin::Subchannel({upper, MMatrix<Complex>(3, 2, 0.0)}, {lower, lower_frame}, left, right, false),
                    std::invalid_argument);
}

TEST_CASE("gra::spin:: projected subchannels equal dense exchange contraction",
          "[gra::spin][subchannel][tensor]") {
  using Complex = std::complex<double>;
  gra::HELMatrix upper;
  upper.lambda_values =
      MMatrix<double>({{-0.5, -1.0}, {-0.5, 0.0}, {0.5, -1.0}, {0.5, 0.0}});
  gra::HELMatrix lower;
  lower.lambda_values =
      MMatrix<double>({{-0.5, 0.0}, {-0.5, 1.0}, {0.5, 0.0}, {0.5, 1.0}});

  MMatrix<Complex> upper_frame(4, 3, 0.0);
  MMatrix<Complex> lower_frame(4, 2, 0.0);
  for (std::size_t row = 0; row < upper_frame.size_row(); ++row) {
    for (std::size_t col = 0; col < upper_frame.size_col(); ++col) {
      upper_frame[row][col] =
          Complex(0.11 * static_cast<double>(1 + row + col),
                  -0.04 * static_cast<double>(1 + 2 * row + col));
    }
    for (std::size_t col = 0; col < lower_frame.size_col(); ++col) {
      lower_frame[row][col] =
          Complex(-0.08 * static_cast<double>(1 + row + 2 * col),
                  0.06 * static_cast<double>(1 + row + col));
    }
  }
  const std::vector<Complex> upper_exchange{
      Complex(0.7, -0.2), Complex(-0.3, 0.5), Complex(0.4, 0.1)};
  const std::vector<Complex> lower_exchange{Complex(-0.6, 0.3),
                                            Complex(0.2, -0.8)};
  const auto left = ToyParticle("f1", 900410, 1, 0.2);
  const auto right = ToyParticle("f2", 900411, 1, 0.3);

  for (const bool analytic : {false, true}) {
    lower.domain = analytic ? gra::HelicityDomain::ReggeTrajectory
                            : gra::HelicityDomain::PhysicalPole;
    for (const bool swap : {false, true}) {
      CAPTURE(analytic, swap);
      const auto dense = gra::spin::Subchannel({upper, upper_frame}, {lower, lower_frame}, left, right, swap);
      const auto projected = gra::spin::ProjectedSubchannel({upper, upper_frame}, {lower, lower_frame}, upper_exchange, lower_exchange, left, right, swap);
      REQUIRE(projected.size() == dense.size_col());
      for (std::size_t final = 0; final < projected.size(); ++final) {
        Complex expected = 0.0;
        for (std::size_t upper_col = 0; upper_col < upper_exchange.size();
             ++upper_col) {
          for (std::size_t lower_col = 0; lower_col < lower_exchange.size();
               ++lower_col) {
            const std::size_t exchange_row =
                upper_col * lower_exchange.size() + lower_col;
            expected += upper_exchange[upper_col] * lower_exchange[lower_col] *
                        dense[exchange_row][final];
          }
        }
        RequireComplexNear(projected[final], expected, 2e-13);
      }
    }
  }
}

TEST_CASE("gra::spin:: reduced continuum omits the physical exchange helicity "
          "barrier",
          "[gra::spin]") {
  for (const bool FORWARD_NOFLIP : std::vector<bool>{true, false}) {
    CAPTURE(FORWARD_NOFLIP);
    gra::LORENTZSCALAR residue_lts =
        MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
    residue_lts.process.MP_FRAME = "HX";
    residue_lts.process.FORWARD_VERTEX =
        gra::ForwardVertexMode::HelicityResidue;
    residue_lts.process.FORWARD_NOFLIP = FORWARD_NOFLIP;
    gra::LORENTZSCALAR pole_lts = residue_lts;
    pole_lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::UnitResidue;

    const auto residue = TestProd4Sum(residue_lts, 1.0);
    const auto pole = TestProd4Sum(pole_lts, 1.0);
    REQUIRE(residue.first.size_row() == (FORWARD_NOFLIP ? 4 : 16));
    REQUIRE(residue.second.size_row() == (FORWARD_NOFLIP ? 4 : 16));
    RequireMatrixNear(pole.first, residue.first, 1e-12);
    RequireMatrixNear(pole.second, residue.second, 1e-12);
    REQUIRE(MatrixNorm2(residue.first) > 1e-12);
    REQUIRE(MatrixNorm2(residue.second) > 1e-12);
    REQUIRE(MatrixNorm2(pole.first) > 1e-12);
    REQUIRE(MatrixNorm2(pole.second) > 1e-12);
  }
}

TEST_CASE("gra::rspin:: Prod4Sum keeps virtual subchannels frame independent",
          "[gra::rspin][subchannel]") {
  gra::LORENTZSCALAR reference_lts =
      MakeToyContinuumLTSAsymmetric(0.3, 4.2, -4.6);
  reference_lts.process.MP_FRAME = "HX";
  const auto reference = TestProd4Sum(reference_lts, 1.0);

  for (const auto &frame : std::vector<std::string>{"CS", "CM"}) {
    CAPTURE(frame);
    gra::LORENTZSCALAR lts = reference_lts;
    lts.process.MP_FRAME = frame;
    const auto actual = TestProd4Sum(lts, 1.0);
    RequireMatrixNear(actual.first, reference.first, 1e-12);
    RequireMatrixNear(actual.second, reference.second, 1e-12);
  }
}

TEST_CASE("gra::spin:: scalar MP proton legs keep compact eikonal rows",
          "[gra::spin]") {
  gra::LORENTZSCALAR lts =
      MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  gra::PARAM_RES res = MakeToyScalarMPResonance();

  const auto mats = gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  REQUIRE(mats.size() == 1);
  REQUIRE(mats[0].size_row() == 4);

  double compact_norm2 = 0.0;
  for (std::size_t i = 0; i < mats[0].size_row(); ++i) {
    for (std::size_t j = 0; j < mats[0].size_col(); ++j) {
      compact_norm2 += std::norm(mats[0][i][j]);
    }
  }
  REQUIRE(compact_norm2 > 1e-12);

  lts.process.FORWARD_NOFLIP = false;
  const auto full_mats =
      gra::rspin::Resonance(lts, res, (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  REQUIRE(full_mats.size() == 1);
  REQUIRE(full_mats[0].size_row() == 16);

  double full_norm2 = 0.0;
  for (std::size_t i = 0; i < full_mats[0].size_row(); ++i) {
    for (std::size_t j = 0; j < full_mats[0].size_col(); ++j) {
      full_norm2 += std::norm(full_mats[0][i][j]);
    }
  }
  REQUIRE(full_norm2 > 1e-12);
  REQUIRE(full_norm2 == Approx(4.0 * compact_norm2).epsilon(1e-12));
}

TEST_CASE("gra::spin:: supported exchange proton legs have unique compact rows",
          "[gra::spin]") {
  struct ProtonLegCase {
    int pdg;
    int spinX2;
    int parity;
    std::size_t l;
    std::size_t two_s;
  };

  const std::vector<ProtonLegCase> cases = {
      {22, 2, -1, 1, 2},  {991, 0, 1, 0, 0},  {993, 2, -1, 1, 2},
      {995, 4, 1, 2, 0},  {9991, 0, 1, 0, 0}, {9993, 2, -1, 1, 2},
  };

  for (const auto &tc : cases) {
    CAPTURE(tc.pdg, tc.spinX2, tc.parity, tc.l, tc.two_s);
    const auto hel = RealisticProtonLegHelicityMatrix(
        tc.pdg, tc.spinX2, tc.parity, tc.l, tc.two_s);
    const auto counts = ActiveProtonLegRowsPerIncoming(hel);

    REQUIRE(counts.size() == 2);
    REQUIRE(counts[0] == 1);
    REQUIRE(counts[1] == 1);
  }
}

TEST_CASE("gra::spin:: scalar MP continuum keeps compact eikonal rows",
          "[gra::spin]") {
  gra::LORENTZSCALAR lts =
      MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6);

  const auto mats = TestProd4Sum(lts, 1.0);
  REQUIRE(mats.first.size_row() == 4);
  REQUIRE(mats.second.size_row() == 4);

  double compact_norm2 = 0.0;
  for (std::size_t i = 0; i < mats.first.size_row(); ++i) {
    for (std::size_t j = 0; j < mats.first.size_col(); ++j) {
      compact_norm2 += std::norm(mats.first[i][j]);
      compact_norm2 += std::norm(mats.second[i][j]);
    }
  }
  REQUIRE(compact_norm2 > 1e-12);

  lts.process.FORWARD_NOFLIP = false;
  const auto full_mats = TestProd4Sum(lts, 1.0);
  REQUIRE(full_mats.first.size_row() == 16);
  REQUIRE(full_mats.second.size_row() == 16);

  double full_norm2 = 0.0;
  for (std::size_t i = 0; i < full_mats.first.size_row(); ++i) {
    for (std::size_t j = 0; j < full_mats.first.size_col(); ++j) {
      full_norm2 += std::norm(full_mats.first[i][j]);
      full_norm2 += std::norm(full_mats.second[i][j]);
    }
  }
  REQUIRE(full_norm2 > 1e-12);
  REQUIRE(full_norm2 == Approx(4.0 * compact_norm2).epsilon(1e-12));
}

// Stable final helicities are summed without an artificial parent spin average
TEST_CASE("gra::spin SPINDEC false preserves stable continuum helicities", "[gra::spin][continuum][physics]") {
  for (const int spin_x2 : {1, 2}) {
    gra::LORENTZSCALAR lts;
    gra::MDecayBranch left, right;
    left.p = ToyParticle("left", 9000012, spin_x2, 1.0);
    right.p = ToyParticle("right", 9000013, spin_x2, 1.0);
    lts.decaytree = {left, right};
    lts.process.SPINDEC = true;
    const auto physical = gra::spin::ContinuumDecayMatrix(lts, "CM");
    lts.process.SPINDEC = false;
    RequireMatrixNear(gra::spin::ContinuumDecayMatrix(lts, "CM"), physical, 1e-14);
  }

  auto lts = ContinuumVectorCascadeLTSForTest(false);
  lts.process.SPINDEC = false;
  lts.decaytree[0].legs.clear();
  const auto vector = gra::spin::ContinuumDecayMatrix(lts, "CM");
  lts.decaytree[0].p.spinX2 = 0;
  const auto scalar = gra::spin::ContinuumDecayMatrix(lts, "CM");
  const auto reference = MMatrix<std::complex<double>>(3, 3, "eye").Kronecker(scalar * scalar.Dagger());
  RequireMatrixNear(vector * vector.Dagger(), reference, 1e-13);
}
