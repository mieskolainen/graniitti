// Eikonal model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Nuclear/MSetup.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Eikonal/MProtonScreen.h"
#include "support/models_test_support.hh"

// Write a controlled Pomeron tune for screening and sampling regression tests
gra::MModelTunePtr EikonalRegressionTune(const std::string &name, const std::size_t n = 1,
                                        const double kappa = 0.0, const bool log_b = false,
                                        const unsigned int b_intervals = 32, const double g0 = 8.0) {
  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  auto numerics = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  auto &soft = general["PARAM_SOFT"];
  soft["active_model"] = "single";
  soft["EXCHANGE_DEF"]["P"]["trajectory_mode"] = "linear";
  auto &model = soft["MODEL"]["single"];
  model["GW"]["theta"] = std::vector<double>(n * (n - 1) / 2, 0.0);
  std::vector<double> resolved(n - 1, 0.0);
  if (!resolved.empty()) { resolved[0] = 1.0; }
  model["GW"]["a_c"] = resolved;
  model["EIKONAL"]["unitarization"] = "exp";
  model["EIKONAL"]["q"] = 1.0;
  model["EIKONAL"]["helicity"] = std::abs(kappa) > 0.0;
  for (auto &exchange : model["EXCHANGE"].items()) {
    exchange.value()["on"] = exchange.key() == "P";
    std::vector<std::vector<double>> coupling(n, std::vector<double>(n, 0.0));
    for (const auto &i : indices(coupling)) { coupling[i][i] = g0 + 2.0 * i; }
    exchange.value()["g"] = coupling;
  }
  model["EXCHANGE"]["P"]["alpha"] = {1.05, 0.2};
  model["EXCHANGE"]["P"]["eta_mode"] = "rotating";
  model["EXCHANGE"]["P"]["sign"] = 1;
  model["EXCHANGE"]["P"]["helicity"]["kappa"] = kappa;
  model["EXCHANGE"]["P"]["helicity"]["B_kappa"] = 0.0;
  for (auto &bank : model["FF"].items()) {
    bank.value()["type"] = "EXP";
    bank.value()["param"] = std::vector<std::vector<double>>(n, {4.0});
  }
  auto &num = numerics["NUMERICS_EIKONAL"];
  num["MinBT"] = 1.0e-5;
  num["MaxBT"] = 2.0;
  num["NumberBT"] = b_intervals;
  num["logBT"] = log_b;
  num["NumberKT2"] = 32;
  num["MinKT2"] = 1.0e-12;
  num["MaxKT2"] = 4.0;
  num["logKT2"] = false;
  num["FBIntegralMinKT"] = 0.0;
  num["FBIntegralMaxKT"] = 4.0;
  num["FBIntegralN"] = 512;
  num["LOOP_INTEGRAL"] = {{"kT_integrator", "GL"}, {"phi_integrator", "Trap"}, {"log_kT", false},
                           {"MinLoopKT", 0.0}, {"MaxLoopKT", 1.0}, {"NumberLoopKT", 4}, {"NumberLoopPHI", 4}};
  const auto directory = std::filesystem::path(gra::aux::GetBasePath(2)) / "tmp" / ("eikonal_" + name);
  std::filesystem::create_directories(directory);
  std::ofstream(directory / "GENERAL.json") << general.dump();
  std::ofstream(directory / "NUMERICS.json") << numerics.dump();
  for (const std::string filename : {"CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}) {
    std::ofstream(directory / filename) << gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", filename));
  }
  return gra::MModelTune::Load((directory / "GENERAL.json").string());
}

// Check the transparent strong limit without requiring a positive cut cross section
TEST_CASE("MEikonal admits a transparent Pomeron", "[gra::MEikonal][transparent][physics]") {
  const auto tune = EikonalRegressionTune("transparent", 1, 0.0, false, 32, 0.0);
  gra::MEikonal eikonal(tune);
  REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState()));
  REQUIRE(eikonal.IsInitialized());
  double total = 0.0, elastic = 0.0, inelastic = 0.0;
  eikonal.GetTotXS(total, elastic, inelastic);
  REQUIRE(total == Approx(0.0).margin(1.0e-14));
  REQUIRE(elastic == Approx(0.0).margin(1.0e-14));
  REQUIRE(inelastic == Approx(0.0).margin(1.0e-14));
  for (const double kt2 : {0.0, 0.02, 0.5, 2.0}) {
    REQUIRE(gra::SquaredNorm(eikonal.GetMatrixRuntime().HelicityMatrix(kt2, 0, 0)) < 1.0e-24);
    REQUIRE(std::abs(eikonal.S3PhysicalScreeningAmp(kt2)) < 1.0e-12);
  }
  std::mt19937 random(31415);
  unsigned int cuts = 0;
  double bt = 0.0;
  REQUIRE_THROWS_AS(eikonal.S3GetRandomCutsBt(cuts, bt, random), std::invalid_argument);
  REQUIRE_THROWS_AS(eikonal.S3GetRandomCuts(random), std::invalid_argument);
}

// Compare both Hankel transforms and integrated cross sections with a Gaussian Regge pole
TEST_CASE("MEikonal reproduces the analytic Gaussian eikonal", "[gra::MEikonal][gaussian][physics]") {
  const std::array<std::pair<std::string, double>, 3> cases = {{{"exp", 1.0}, {"q_exp", 0.63}, {"q_exp", 0.3}}};
  for (const auto &[unitarization, q] : cases) {
    const auto seed = EikonalRegressionTune("gaussian_" + std::to_string(q), 1, 0.0, false, 256);
    auto general = seed->General();
    auto numerics = seed->Numerics();
    auto &model = general["PARAM_SOFT"]["MODEL"]["single"];
    model["EXCHANGE"]["P"]["eta_mode"] = "rotating_t0";
    model["EIKONAL"]["unitarization"] = unitarization;
    model["EIKONAL"]["q"] = q;
    numerics["NUMERICS_EIKONAL"]["MaxBT"] = 5.0;
    numerics["NUMERICS_EIKONAL"]["FBIntegralN"] = 1024;
    std::ofstream(seed->GeneralFile()) << general.dump();
    std::ofstream(seed->NumericsFile()) << numerics.dump();
    const auto tune = gra::MModelTune::Load(seed->GeneralFile());
    const auto beams = ProtonInitialState();
    constexpr double s = 25.0;
    const double flux = s * gra::kinematics::beta12(s, beams[0].mass, beams[1].mass);
    const auto pomeron = tune->Soft()->ExchangeId("P");
    const auto &pole = tune->Soft()->Exchange(pomeron);
    const double slope = pole.FormFactorParameters()[0][0] + pole.AlphaPrime() * std::log(s);
    const auto born = gra::MEikonal::SingleAmpElastic(tune->Soft(), s, 0.0, pomeron, 0, 0);
    const auto chi0 = born / (8.0 * gra::math::PI * flux * slope);
    gra::MEikonal eikonal(tune);
    REQUIRE_NOTHROW(eikonal.S3Constructor(s, beams));
    const auto &runtime = eikonal.GetMatrixRuntime();
    for (const double b : runtime.ImpactParameterNodes()) {
      const auto expected = chi0 * std::exp(-b * b / (4.0 * slope));
      CAPTURE(q, b);
      REQUIRE(std::abs(eikonal.S3Density(b, 0, 0) - expected) < 2.0e-7 * std::abs(chi0));
    }
    // Integrate each Gaussian term of i[1-exp_q(i chi0 exp(-b^2/(4a)))] exactly
    const auto amplitude = [&](const double kt2) {
      const auto opacity = -gra::math::zi * chi0;
      std::complex<double> coefficient = opacity;
      std::complex<double> sum = 0.0;
      for (unsigned int n = 1; n <= 64; ++n) {
        sum += coefficient * std::exp(-slope * kt2 / n) / static_cast<double>(n);
        coefficient *= -opacity * (1.0 - n * (1.0 - q)) / static_cast<double>(n + 1);
      }
      return 8.0 * gra::math::PI * flux * slope * gra::math::zi * sum;
    };
    const double scale = std::abs(amplitude(0.0));
    for (const double kt2 : runtime.MomentumTransferNodes()) {
      CAPTURE(q, kt2);
      REQUIRE(std::abs(eikonal.S3ExclusiveAmp(kt2, 0, 0) - amplitude(kt2)) < 2.0e-6 * scale);
    }
    double total = 0.0, elastic = 0.0, inelastic = 0.0;
    eikonal.GetTotXS(total, elastic, inelastic);
    REQUIRE(total == Approx(std::imag(amplitude(0.0)) * gra::PDG::GeV2barn / flux).epsilon(2.0e-6));
    const auto [node, weight] = gra::math::GaussLegendreRule(128, 0.0, eikonal.Numerics.MaxKT2);
    double integrated = 0.0;
    for (const auto &i : indices(node)) { integrated += weight[i] * std::norm(amplitude(node[i])); }
    integrated *= gra::PDG::GeV2barn / (16.0 * gra::math::PI * flux * flux);
    REQUIRE(elastic == Approx(integrated).epsilon(2.0e-6));
  }
}

// Check peripheral coupled channel unitarity against a refined Born quadrature
TEST_CASE("MEikonal peripheral profiles converge with the Born quadrature", "[gra::MEikonal][peripheral][physics]") {
  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  general["PARAM_SOFT"]["active_model"] = "double";
  const auto model = gra::SoftModel::LoadFromJson(modelfile, general.dump());
  gra::MEikonalNumerics numerics;
  numerics.ReadParameters(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"));
  numerics.NumberBT = 128;
  numerics.NumberKT2 = 2;
  const auto coarse = gra::MEikonalMatrix::Build(model, 7000.0 * 7000.0, ProtonInitialState(), numerics, false);
  numerics.FBIntegralN *= 2;
  const auto fine = gra::MEikonalMatrix::Build(model, 7000.0 * 7000.0, ProtonInitialState(), numerics, false);
  const double tolerance = numerics.unitarity_tolerance;
  CAPTURE(coarse->MaxSingularValue(), fine->MaxSingularValue(), tolerance);
  REQUIRE(std::abs(coarse->MaxSingularValue() - fine->MaxSingularValue()) < tolerance);
  for (const double b : coarse->ImpactParameterNodes()) {
    for (std::size_t f1 = 0; f1 < model->GoodWalker().ChannelCount(); ++f1) {
      for (std::size_t f2 = 0; f2 < model->GoodWalker().ChannelCount(); ++f2) {
        const auto difference = coarse->PhysicalImpactSMatrix(b, f1, f2) - fine->PhysicalImpactSMatrix(b, f1, f2);
        CAPTURE(b, f1, f2);
        REQUIRE(difference.FrobNorm() < tolerance);
      }
    }
  }
}

// Compare the pair-space operator with its physical proton projection
void RequirePhysicalPairScreeningProjection(const gra::MEikonal &eikonal, const double kt2) {
  const auto        spin   = eikonal.PairScreeningSpinBank(kt2);
  const auto        named  = eikonal.PairScreeningHelicityBank(kt2);
  const auto       &mixing = eikonal.GetMixingMatrix();
  const std::size_t n      = eikonal.GetChannelCount();
  REQUIRE(mixing.size_row() == n);
  REQUIRE(mixing.size_col() == n);
  std::vector<double> proton_pair(n * n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) { proton_pair[i * n + j] = mixing(0, i) * mixing(0, j); }
  }
  const auto project = [&](const auto &matrix) {
    std::complex<double> value = 0.0;
    for (std::size_t i = 0; i < n * n; ++i) {
      for (std::size_t j = 0; j < n * n; ++j) { value += proton_pair[i] * matrix(i, j) * proton_pair[j]; }
    }
    return value;
  };
  const auto physical_spin = eikonal.GetMatrixRuntime().ScreeningHelicityMatrix(kt2, 0, 0);
  for (std::size_t entry = 0; entry < spin.size(); ++entry) {
    CAPTURE(n, entry, kt2);
    const auto actual = project(spin[entry]);
    REQUIRE(std::real(actual) == Approx(std::real(physical_spin[entry])).epsilon(1.0e-11).margin(1.0e-12));
    REQUIRE(std::imag(actual) == Approx(std::imag(physical_spin[entry])).epsilon(1.0e-11).margin(1.0e-12));
  }

  const auto physical = eikonal.GetMatrixRuntime().ScreeningHelicityAmplitudes(kt2, 0, 0);
  const std::array<std::complex<double>, 6>          expected  = {physical.phi1, physical.phi2, physical.phi3,
                                                                  physical.phi4, physical.phi5, physical.phi5_first};
  const std::array<gra::ElasticHelicityComponent, 6> component = {
      gra::ElasticHelicityComponent::Phi1, gra::ElasticHelicityComponent::Phi2,
      gra::ElasticHelicityComponent::Phi3, gra::ElasticHelicityComponent::Phi4,
      gra::ElasticHelicityComponent::Phi5, gra::ElasticHelicityComponent::Phi5First};
  for (std::size_t h = 0; h < expected.size(); ++h) {
    CAPTURE(n, h, kt2);
    const auto actual = project(named.Component(component[h]));
    REQUIRE(std::real(actual) == Approx(std::real(expected[h])).epsilon(1.0e-11).margin(1.0e-12));
    REQUIRE(std::imag(actual) == Approx(std::imag(expected[h])).epsilon(1.0e-11).margin(1.0e-12));
  }
}

// Construct one elastic proton pair in the center-of-mass frame
std::array<gra::M4Vec, 4> ElasticProtonMomenta(const double momentum, const double abs_t, const double azimuth) {
  const double mass   = ProtonInitialState()[0].mass;
  const double energy = std::sqrt(momentum * momentum + mass * mass);
  const double transverse = std::sqrt(abs_t * (1.0 - abs_t / (4.0 * momentum * momentum)));
  const double longitudinal = momentum - abs_t / (2.0 * momentum);
  return {
      gra::M4Vec(0.0, 0.0, momentum, energy), gra::M4Vec(0.0, 0.0, -momentum, energy),
      gra::M4Vec(transverse * std::cos(azimuth), transverse * std::sin(azimuth), longitudinal, energy),
      gra::M4Vec(-transverse * std::cos(azimuth), -transverse * std::sin(azimuth), -longitudinal, energy)};
}

// Check CNI mass shell roundoff, azimuthal covariance and rejection of physical mass shifts
TEST_CASE("Elastic CNI accepts on shell collider momenta", "[gra::MEikonal][CNI][validation]") {
  const auto tune = EikonalRegressionTune("cni_shell");
  const auto beams = ProtonInitialState();
  for (const double momentum : {6500.0, 7000.0, 50000.0}) {
    CAPTURE(momentum);
    const double s = 4.0 * (momentum * momentum + beams[0].mass * beams[0].mass);
    gra::MEikonal eikonal(tune);
    eikonal.S3Constructor(s, beams);
    eikonal.Numerics.CNI.NumberKT2 = 64;
    eikonal.Numerics.CNI.FBIntegralMaxKT = 2.0;
    eikonal.Numerics.CNI.FBIntegralN = 256;
    eikonal.Numerics.CNI.tail_max_qb = 128.0;
    eikonal.Numerics.CNI.TailIntegralN = 512;
    // This test isolates event validation from the separately tested transform convergence
    eikonal.Numerics.CNI.interp_rel_tol = 1.0e3;
    eikonal.InitializeElasticCNI(1.0e-4, 1.1);
    for (const double abs_t : {1.0e-4, 0.01, 0.1, 1.0}) {
      CAPTURE(abs_t);
      const auto p = ElasticProtonMomenta(momentum, abs_t, 0.0);
      const auto rotated = ElasticProtonMomenta(momentum, abs_t, 0.63);
      const auto amplitude = eikonal.PhysicalElasticHelicityMatrix(p[0], p[1], p[2], p[3]);
      const auto actual = eikonal.PhysicalElasticHelicityMatrix(rotated[0], rotated[1], rotated[2], rotated[3]);
      const auto expected = gra::RotateProtonHelicityMatrix(amplitude, 0.63);
      REQUIRE(gra::AllFinite(actual));
      REQUIRE(gra::SquaredNorm(actual) > 0.0);
      for (const auto &entry : indices(actual)) {
        REQUIRE(std::abs(actual[entry] - expected[entry]) <= 1.0e-8 * std::max(1.0, std::abs(expected[entry])));
      }
      auto off_shell = p;
      // Keep total four momentum fixed while shifting each outgoing mass
      off_shell[2].SetE(p[2].E() + 1.0e-3);
      off_shell[3].SetE(p[3].E() - 1.0e-3);
      REQUIRE_THROWS_AS(eikonal.PhysicalElasticHelicityMatrix(off_shell[0], off_shell[1], off_shell[2], off_shell[3]),
                        std::invalid_argument);
    }
  }
}

// Compare mixed scalar and spin dependent channels with the full physical matrix contraction
TEST_CASE("Good Walker screening handles scalar and helicity channels together", "[gra::MEikonal][helicity][GoodWalker]") {
  gra::MEikonal eikonal(EikonalRegressionTune("mixed_spin", 3, 0.01));
  eikonal.S3Constructor(25.0, ProtonInitialState());
  const auto &loop = eikonal.GetLoopConst(25.0);
  const auto type = gra::GWFinalClass::SingleDissociation1;
  const auto &channels = loop.good_walker_channels[static_cast<std::size_t>(type)];
  REQUIRE(channels.size() == 2);
  REQUIRE_FALSE(channels[0].spin_scalar);
  REQUIRE(channels[1].spin_scalar);
  REQUIRE(channels[1].screening_weight.empty());
  gra::ScreeningMetadata metadata(gra::ScreeningSpinBasis::ProtonHelicity,
                                   gra::ProtonScreeningMode::ForwardExcitation, 0.25,
                                   gra::ScreeningAmplitudeType::Physical, 16);
  metadata.PrepareSpinTransitions();
  std::vector<std::complex<double>> born(16);
  for (const auto &h : indices(born)) { born[h] = {0.1 * (h + 1), -0.03 * h}; }
  gra::eikonal::MProtonScreen screen({born}, metadata, loop, type);
  std::vector<std::complex<double>> expected(channels.size() * born.size());
  for (const auto &c : indices(channels)) {
    for (const auto &h : indices(born)) { expected[c * born.size() + h] = channels[c].born_coefficient * born[h]; }
  }
  const std::vector<std::span<const std::complex<double>>> shifted = {born};
  for (const auto &radial : indices(loop.kt2)) {
    for (std::size_t azimuth = 0; azimuth < loop.node_weight.size_col(); ++azimuth) {
      screen.Add(radial, azimuth, shifted);
      const double phi = std::atan2(loop.kt_y[radial][azimuth], loop.kt_x[radial][azimuth]);
      for (const auto &c : indices(channels)) {
        const auto f = channels[c].final_state;
        gra::ProtonHelicityMatrix soft{};
        for (const auto &channel : channels) {
          const auto i = channel.final_state;
          const auto matrix = eikonal.GetMatrixRuntime().ScreeningHelicityMatrix(loop.kt2[radial], f.f1, f.f2,
                                                                                i.f1, i.f2, phi);
          gra::AddScaled(soft, matrix, channel.born_coefficient);
        }
        for (std::size_t initial = 0; initial < 4; ++initial) {
          for (std::size_t final = 0; final < 4; ++final) {
            for (std::size_t intermediate = 0; intermediate < 4; ++intermediate) {
              expected[c * 16 + gra::spin::PairHelicityTransitionIndex(initial, final)] +=
                  loop.node_weight[radial][azimuth] * soft[gra::spin::PairHelicityMatrixIndex(final, intermediate)] *
                  born[gra::spin::PairHelicityTransitionIndex(initial, intermediate)];
            }
          }
        }
      }
    }
  }
  const auto actual = screen.Result().front();
  REQUIRE(actual.size() == expected.size());
  for (const auto &h : indices(actual)) {
    REQUIRE(std::abs(actual[h] - expected[h]) <= 1.0e-11 * std::max(1.0, std::abs(expected[h])));
  }
}

TEST_CASE("Elastic CNI uses a symmetric full-spin sandwich", "[gra::MEikonal][CNI][matrix]") {
  using Complex   = std::complex<double>;
  using Matrix    = gra::MElasticCNI::Matrix;
  Matrix strong_s = Matrix::IdentityMatrix(4);
  strong_s(0, 1)  = Complex(0.07, -0.02);
  strong_s(1, 0)  = Complex(-0.03, 0.04);
  strong_s(2, 3)  = Complex(0.02, 0.01);
  Matrix chi(4);
  chi(0, 0) = Complex(0.03, 0.01);
  chi(1, 2) = Complex(-0.02, 0.04);
  chi(2, 1) = Complex(0.01, -0.03);

  const Matrix half     = (chi * (0.5 * gra::math::zi)).Exp();
  const Matrix expected = half * strong_s * half;
  const Matrix actual   = gra::MElasticCNI::SymmetricShortRangeSMatrix(strong_s, chi);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) { REQUIRE(std::abs(actual(row, col) - expected(row, col)) < 2.0e-14); }
  }

  // Recover the commuting scalar all-orders Coulomb and nuclear S matrix
  // [REFERENCE: Kaspar, arXiv:2001.10227]
  const Complex scalar_strong(0.73, -0.11);
  const Complex scalar_chi(0.08, 0.03);
  const Matrix  scalar_strong_s = Matrix::IdentityMatrix(4) * scalar_strong;
  const Matrix  scalar_gamma    = Matrix::IdentityMatrix(4) * scalar_chi;
  const Matrix  scalar_combined = gra::MElasticCNI::SymmetricShortRangeSMatrix(scalar_strong_s, scalar_gamma);
  const Matrix  scalar_expected = Matrix::IdentityMatrix(4) * (std::exp(gra::math::zi * scalar_chi) * scalar_strong);
  REQUIRE((scalar_combined - scalar_expected).FrobNorm() < 2.0e-14);

  const Matrix zero(4);
  const Matrix pure_point = gra::MElasticCNI::FiniteRemainderProfile(Matrix::IdentityMatrix(4), zero, 0.37);
  REQUIRE(pure_point.FrobNorm() < 1.0e-14);

  const double epsilon    = 1.0e-7;
  const Matrix derivative = gra::MElasticCNI::FiniteRemainderProfile(strong_s, chi * epsilon, 0.0) / Complex(epsilon);
  const Matrix derivative_expected = (chi * strong_s + strong_s * chi) * 0.5 - chi;
  REQUIRE((derivative - derivative_expected).FrobNorm() < 2.0e-8);
}

TEST_CASE("Analytic point Coulomb factor is a pure phase", "[gra::MEikonal][CNI][point]") {
  const double eta   = gra::qed::alpha_0 * 0.997;
  const auto   phase = gra::MElasticCNI::RenormalizedPointCoulombPhase(eta, 1.0, 1.0e-3);
  REQUIRE(std::abs(phase) == Approx(1.0).margin(2.0e-15));
  REQUIRE(std::arg(phase) / eta == Approx(std::log(1.0e3)).epsilon(1.0e-4));
  REQUIRE(std::abs(gra::MElasticCNI::RenormalizedPointCoulombPhase(-eta, 1.0, 1.0e-3) - std::conj(phase)) < 2.0e-15);
  REQUIRE(gra::MElasticCNI::RenormalizedPointCoulombPhase(0.0, 1.0, 1.0e-3) == std::complex<double>(1.0, 0.0));
}

TEST_CASE("Exact elastic photon matrix has Rutherford and azimuth limits", "[gra::MEikonal][CNI][photon]") {
  const double      momentum  = 100.0;
  const double      abs_t     = 1.0e-4;
  const double      azimuth   = 0.63;
  const auto        reference = ElasticProtonMomenta(momentum, abs_t, 0.0);
  const auto        rotated   = ElasticProtonMomenta(momentum, abs_t, azimuth);
  const auto        pp_state  = ProtonInitialState();
  const gra::MDirac dirac("DIRAC");
  // Evaluate the reusable exact QED pair-helicity API for one beam state
  const auto photon_matrix = [&dirac](const std::vector<gra::MParticle> &state,
                                      const std::array<gra::M4Vec, 4>   &momenta) {
    return gra::qed::ElasticSpinHalfPhotonExchange(dirac, state[0], state[1], momenta[0], momenta[1], momenta[2],
                                                   momenta[3]);
  };
  const auto pp               = photon_matrix(pp_state, reference);
  const auto pp_rotated       = photon_matrix(pp_state, rotated);
  const auto expected_rotated = gra::RotateProtonHelicityMatrix(pp, azimuth);
  for (std::size_t entry = 0; entry < pp.size(); ++entry) {
    CAPTURE(entry);
    REQUIRE(std::abs(pp_rotated[entry] - expected_rotated[entry]) < 2.0e-9 * std::max(1.0, std::abs(pp[entry])));
  }

  const double energy = std::sqrt(momentum * momentum + gra::PDG::mp * gra::PDG::mp);
  const double s      = 4.0 * energy * energy;
  const double point  = -8.0 * gra::math::PI * gra::qed::alpha_0 * (s - 2.0 * gra::PDG::mp * gra::PDG::mp) / abs_t;
  const double tau    = abs_t / (4.0 * gra::PDG::mp * gra::PDG::mp);
  const double effective_form_factor2 =
      (gra::math::pow2(gra::form::G_E(abs_t)) + tau * gra::math::pow2(gra::form::G_M(abs_t))) / (1.0 + tau);
  const double spin_averaged = 0.25 * gra::SquaredNorm(pp);
  REQUIRE(std::real(pp[0] / point) == Approx(1.0).epsilon(2.0e-3));
  REQUIRE(spin_averaged == Approx(point * point * gra::math::pow2(effective_form_factor2)).epsilon(5.0e-4));

  const double beta       = 2.0 * momentum / std::sqrt(s);
  const double dsigma_dt  = spin_averaged / (16.0 * gra::math::PI * gra::math::pow2(s * beta));
  const double rutherford = 4.0 * gra::math::PI * gra::math::pow2(gra::qed::alpha_0) / gra::math::pow2(abs_t);
  REQUIRE(dsigma_dt == Approx(rutherford * gra::math::pow2(effective_form_factor2)).epsilon(5.0e-4));

  const double f1 = gra::form::F1(-abs_t);
  const double f2 = gra::form::F2(-abs_t);
  const double buttimore_phi2 =
      8.0 * gra::math::PI * gra::qed::alpha_0 * s * f2 * f2 / (4.0 * gra::PDG::mp * gra::PDG::mp);
  const double buttimore_phi5 =
      -8.0 * gra::math::PI * gra::qed::alpha_0 * s * f1 * f2 / (2.0 * gra::PDG::mp * std::sqrt(abs_t));
  REQUIRE(std::real(pp[gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi2)] / buttimore_phi2) ==
          Approx(1.0).epsilon(2.0e-2));
  REQUIRE(std::real(pp[gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi4)] / buttimore_phi2) ==
          Approx(-1.0).epsilon(2.0e-2));
  const auto phi5_reference = gra::ElasticHelicityReferenceTransition(gra::ElasticHelicityComponent::Phi5);
  REQUIRE(std::real(phi5_reference.sign * pp[gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi5)] /
                    buttimore_phi5) == Approx(1.0).epsilon(2.0e-2));

  auto ppbar_state        = pp_state;
  ppbar_state[1].name     = "pbar";
  ppbar_state[1].pdg      = -2212;
  ppbar_state[1].chargeX3 = -3;
  const auto ppbar        = photon_matrix(ppbar_state, reference);
  REQUIRE(std::real(ppbar[0] / pp[0]) == Approx(-1.0).epsilon(2.0e-3));
  REQUIRE(gra::SquaredNorm(ppbar) == Approx(gra::SquaredNorm(pp)).epsilon(2.0e-12));

  gra::ProtonHelicityMatrix strong{};
  for (std::size_t diagonal = 0; diagonal < 4; ++diagonal) { strong[5 * diagonal] = 0.1 * std::abs(point); }
  const auto interference = [&strong](const auto &photon) {
    auto coherent = strong;
    gra::AddScaled(coherent, photon, std::complex<double>(1.0));
    return gra::SquaredNorm(coherent) - gra::SquaredNorm(strong) - gra::SquaredNorm(photon);
  };
  REQUIRE(interference(pp) < 0.0);
  REQUIRE(interference(ppbar) > 0.0);

  auto pbarp_state = ppbar_state;
  std::swap(pbarp_state[0], pbarp_state[1]);
  const auto pbarp = photon_matrix(pbarp_state, reference);
  REQUIRE(gra::SquaredNorm(pbarp) == Approx(gra::SquaredNorm(pp)).epsilon(2.0e-12));

  auto pbarpbar_state = pp_state;
  pbarpbar_state[0]   = ppbar_state[1];
  pbarpbar_state[1]   = ppbar_state[1];
  const auto pbarpbar = photon_matrix(pbarpbar_state, reference);
  REQUIRE(gra::SquaredNorm(pbarpbar) == Approx(gra::SquaredNorm(pp)).epsilon(2.0e-12));

  auto rejected_state        = pp_state;
  rejected_state[1].pdg      = gra::PDG::PDG_n;
  rejected_state[1].chargeX3 = 0;
  // Reject unsupported emitters at initialization, before the photon amplitude
  gra::LORENTZSCALAR rejected_lts;
  rejected_lts.beam1 = rejected_state[0];
  rejected_lts.beam2 = rejected_state[1];
  REQUIRE_THROWS_AS(gra::qed::ValidateEmitter(rejected_lts, 2, "elastic photon test"), std::invalid_argument);
}

TEST_CASE("Elastic CNI converges under simultaneous transform refinement", "[gra::MEikonal][CNI][convergence]") {
  constexpr double s            = 25.0;
  constexpr double min_abs_t    = 0.01;
  constexpr double max_abs_t    = 0.2;
  const auto       initialstate = ProtonInitialState();
  const double     momentum = gra::kinematics::DecayMomentum(std::sqrt(s), initialstate[0].mass, initialstate[1].mass);
  const std::array<double, 3> transfer = {0.02, 0.07, 0.16};
  const double                beta     = 2.0 * momentum / std::sqrt(s);

  auto baseline                         = BuildTestEikonal({}, "cni_convergence_base", 128, 64);
  baseline.Numerics.CNI.NumberKT2       = 192;
  baseline.Numerics.CNI.logKT2          = true;
  baseline.Numerics.CNI.FBIntegralMaxKT = 2.0;
  baseline.Numerics.CNI.FBIntegralN     = 256;
  baseline.Numerics.CNI.tail_max_qb     = 128.0;
  baseline.Numerics.CNI.TailIntegralN   = 512;
  baseline.Numerics.CNI.interp_rel_tol  = 1.0e3;
  baseline.InitializeElasticCNI(min_abs_t, max_abs_t);

  const auto   cni   = gra::MElasticCNI::Build(baseline.GetMatrixRuntime(), baseline.SoftModelHandle(), s, initialstate,
                                               baseline.Numerics.CNI, baseline.Numerics.logBT, min_abs_t, max_abs_t,
                                               baseline.ModelTuneHandle()->Structure());
  const double delta = s - 2.0 * gra::math::pow2(initialstate[0].mass);
  const double eta   = gra::qed::alpha_0 * delta / (s * beta);
  const double point_coeff    = -8.0 * gra::math::PI * gra::qed::alpha_0 * delta;
  const double point_t        = transfer[1];
  const auto   point          = cni->PointHigherOrderHelicityMatrix(point_t);
  const auto   phase          = gra::MElasticCNI::RenormalizedPointCoulombPhase(eta, 1.0, point_t);
  const auto   expected_point = point_coeff * (phase - std::complex<double>(1.0, 0.0)) / point_t;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const auto expected = row == col ? expected_point : std::complex<double>(0.0, 0.0);
      REQUIRE(std::abs(point[gra::spin::PairHelicityMatrixIndex(row, col)] - expected) <
              2.0e-13 * std::max(1.0, std::abs(expected)));
    }
  }

  auto refined                         = BuildTestEikonal({}, "cni_convergence_refined", 256, 64);
  refined.Numerics.CNI.NumberKT2       = 384;
  refined.Numerics.CNI.logKT2          = true;
  refined.Numerics.CNI.FBIntegralMaxKT = 2.0;
  refined.Numerics.CNI.FBIntegralN     = 512;
  refined.Numerics.CNI.tail_max_qb     = 128.0;
  refined.Numerics.CNI.TailIntegralN   = 1024;
  refined.Numerics.CNI.interp_rel_tol  = 1.0e3;
  refined.InitializeElasticCNI(min_abs_t, max_abs_t);

  for (const double abs_t : transfer) {
    CAPTURE(abs_t);
    const auto   momenta = ElasticProtonMomenta(momentum, abs_t, 0.37);
    const auto   coarse  = baseline.PhysicalElasticHelicityMatrix(momenta[0], momenta[1], momenta[2], momenta[3]);
    const auto   fine    = refined.PhysicalElasticHelicityMatrix(momenta[0], momenta[1], momenta[2], momenta[3]);
    const double flux    = 16.0 * gra::math::PI * gra::math::pow2(s * beta);
    const double coarse_dsigma_dt = 0.25 * gra::SquaredNorm(coarse) / flux;
    const double fine_dsigma_dt   = 0.25 * gra::SquaredNorm(fine) / flux;
    REQUIRE(std::abs(coarse_dsigma_dt / fine_dsigma_dt - 1.0) < 5.0e-4);

    double amplitude_scale = 0.0;
    for (const auto &amplitude : fine) { amplitude_scale = std::max(amplitude_scale, std::abs(amplitude)); }
    REQUIRE(amplitude_scale > 0.0);
    for (std::size_t entry = 0; entry < fine.size(); ++entry) {
      CAPTURE(entry);
      const double physical_floor = 1.0e-12 * amplitude_scale;
      const double scale          = std::max(std::abs(fine[entry]), physical_floor);
      REQUIRE(std::abs(coarse[entry] - fine[entry]) / scale < 1.0e-3);
    }
  }
}

// Resolve the nonlinear Coulomb spin tail below and above qb=1 with the physical matrix API
TEST_CASE("Elastic CNI spin tails converge at very small transfer", "[gra::MEikonal][CNI][convergence]") {
  const auto    tune     = EikonalRegressionTune("cni_small_t", 1, 0.0, false, 128);
  const auto    beams    = ProtonInitialState();
  const double  momentum = 6500.0;
  const double  s        = 4.0 * (momentum * momentum + beams[0].mass * beams[0].mass);
  gra::MEikonal eikonal(tune);
  eikonal.S3Constructor(s, beams);
  auto numerics            = eikonal.Numerics.CNI;
  numerics.NumberKT2       = 192;
  numerics.FBIntegralMaxKT = 2.0;
  numerics.FBIntegralN     = 512;
  numerics.tail_max_qb     = 2048.0;
  numerics.TailIntegralN   = 8192;
  numerics.interp_rel_tol  = 5.0e-5;
  const auto coarse        = gra::MElasticCNI::Build(eikonal.GetMatrixRuntime(), tune->Soft(), s, beams, numerics,
                                                     eikonal.Numerics.logBT, 1.0e-10, 0.01, tune->Structure());
  numerics.TailIntegralN *= 4;
  numerics.tail_order *= 2;
  numerics.tail_split_qb *= 0.5;
  const auto fine = gra::MElasticCNI::Build(eikonal.GetMatrixRuntime(), tune->Soft(), s, beams, numerics,
                                            eikonal.Numerics.logBT, 1.0e-10, 0.01, tune->Structure());
  numerics.tail_max_qb *= 2;
  numerics.TailIntegralN *= 2;
  const auto extended = gra::MElasticCNI::Build(eikonal.GetMatrixRuntime(), tune->Soft(), s, beams, numerics,
                                                eikonal.Numerics.logBT, 1.0e-10, 0.01, tune->Structure());
  for (const double abs_t : {5.0e-11, 1.0e-12, 0.011}) {
    const auto p = ElasticProtonMomenta(momentum, abs_t, 0.37);
    CAPTURE(abs_t);
    REQUIRE_THROWS_AS(coarse->PhysicalHelicityMatrix(p[0], p[1], p[2], p[3]), std::out_of_range);
  }
  for (const double abs_t : {1.0e-10, 1.0e-8, 1.0e-6, 1.0e-4, 0.01}) {
    const auto p = ElasticProtonMomenta(momentum, abs_t, 0.37);
    const auto a = coarse->PhysicalHelicityMatrix(p[0], p[1], p[2], p[3]);
    const auto b = fine->PhysicalHelicityMatrix(p[0], p[1], p[2], p[3]);
    const auto c = extended->PhysicalHelicityMatrix(p[0], p[1], p[2], p[3]);
    for (const auto &entry : indices(a)) {
      CAPTURE(abs_t, entry, a[entry], b[entry], c[entry]);
      REQUIRE(std::abs(b[entry]) > 0.0);
      REQUIRE(std::abs(a[entry] - b[entry]) < 1.0e-8 * std::abs(b[entry]));
      REQUIRE(std::abs(c[entry] - b[entry]) < 1.0e-6 * std::abs(b[entry]));
    }
  }
}

TEST_CASE("Poisson veto uses the square root of the zero-secondary probability", "[gra::form][physics]") {
  gra::regge::VetoParam veto;
  veto.active = true;
  veto.M0     = 1.275;
  veto.c      = 0.7;

  CHECK(gra::regge::Veto(veto, 1.0) == Approx(1.0));
  CHECK(gra::regge::Veto(veto, veto.M0) == Approx(1.0));

  const double central_mass = 2.0;
  const double amplitude    = gra::regge::Veto(veto, central_mass);
  CHECK(amplitude == Approx(std::pow(veto.M0 / central_mass, veto.c)).margin(1e-15));
  CHECK(
      amplitude * amplitude ==
      Approx(std::exp(-veto.c * std::log(gra::math::pow2(central_mass / veto.M0)))).margin(1e-15));

  veto.active = false;
  CHECK(gra::regge::Veto(veto, central_mass) == Approx(1.0));
}

TEST_CASE("MRegge triple-Regge amplitudes use the triple-Pomeron eta mode", "[gra::MRegge][gra::MEikonal][physics]") {
  const auto rotating_tune  = WriteModifiedPhotoVMTune("triple_regge_rotating", [](auto &j) {
    const std::string model                                          = j.at("PARAM_SOFT").at("active_model");
    j.at("PARAM_SOFT").at("MODEL").at(model).at("3P").at("eta_mode") = "rotating";
  });
  const auto raw_tune       = WriteModifiedPhotoVMTune("triple_regge_raw", [](auto &j) {
    const std::string model                                          = j.at("PARAM_SOFT").at("active_model");
    j.at("PARAM_SOFT").at("MODEL").at(model).at("3P").at("eta_mode") = "raw";
  });
  const auto rotating_model = gra::MModelTune::Load(rotating_tune.second);
  const auto raw_model      = gra::MModelTune::Load(raw_tune.second);

  gra::LORENTZSCALAR reference_lts = MakeToyCoherentPhotonLTS();
  reference_lts.ss[1][1]           = 4.0;
  reference_lts.ss[2][2]           = 9.0;
  gra::LORENTZSCALAR raw_lts       = reference_lts;
  gra::MRegge        reference(reference_lts, rotating_model,
                               gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  gra::MRegge central_raw(raw_lts, raw_model, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));

  for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
    CAPTURE(mode);
    reference_lts.excite1 = true;
    reference_lts.excite2 = mode == gra::MReggeInclusive::DD;
    raw_lts.excite1       = true;
    raw_lts.excite2       = mode == gra::MReggeInclusive::DD;
    REQUIRE(std::abs(reference.ME2(reference_lts, mode)) > 0.0);
    REQUIRE(std::abs(central_raw.ME2(raw_lts, mode)) > 0.0);
    REQUIRE(reference_lts.proton_good_walker.has_value());
    REQUIRE(raw_lts.proton_good_walker.has_value());
    const auto &rotating_source = reference_lts.proton_good_walker->components.front().source;
    const auto &raw_source      = raw_lts.proton_good_walker->components.front().source;
    REQUIRE(rotating_source.size_row() == raw_source.size_row());
    REQUIRE(rotating_source.size_col() == raw_source.size_col());

    std::size_t row = 0;
    std::size_t col = 0;
    while (row < rotating_source.size_row() && std::abs(rotating_source(row, col)) <= 1.0e-15) {
      if (++col == rotating_source.size_col()) {
        col = 0;
        ++row;
      }
    }
    REQUIRE(row < rotating_source.size_row());
    const auto                 pomeron_id = rotating_model->Soft()->ExchangeId("P");
    const auto                &pomeron    = rotating_model->Soft()->Exchange(pomeron_id);
    const double               alpha      = rotating_model->Soft()->Alpha(pomeron_id, reference_lts.t);
    const std::complex<double> eta_ratio =
        gra::regge::EtaFactor(alpha, pomeron.Alpha0(), pomeron.Signature(), gra::EtaMode::Raw) /
        gra::regge::EtaFactor(alpha, pomeron.Alpha0(), pomeron.Signature(), gra::EtaMode::Rotating);
    const std::complex<double> expected = rotating_source(row, col) * eta_ratio;
    REQUIRE(std::real(raw_source(row, col)) == Approx(std::real(expected)).epsilon(1.0e-12).margin(1.0e-14));
    REQUIRE(std::imag(raw_source(row, col)) == Approx(std::imag(expected)).epsilon(1.0e-12).margin(1.0e-14));
  }
}

// Check inclusive amplitudes through the first unphysical spacelike signature pole
TEST_CASE("Inclusive diffraction stays smooth across zero Pomeron trajectory", "[gra::MRegge][physics][icepack]") {
  const auto tune = gra::MModelTune::Load(modelfile);
  const auto &soft = *tune->Soft();
  const auto pomeron = soft.ExchangeId("P");
  double lower = -50.0;
  double upper = 0.0;
  REQUIRE(soft.Alpha(pomeron, lower) < 0.0);
  REQUIRE(soft.Alpha(pomeron, upper) > 0.0);
  for (unsigned int i = 0; i < 60; ++i) {
    const double mid = (lower + upper) / 2.0;
    if (soft.Alpha(pomeron, mid) < 0.0) { lower = mid; } else { upper = mid; }
  }
  auto lts = MakeToyCoherentPhotonLTS();
  lts.ss[1][1] = 4.0;
  lts.ss[2][2] = 9.0;
  gra::MRegge regge(lts, tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  for (const auto mode : {gra::MReggeInclusive::SD, gra::MReggeInclusive::DD}) {
    lts.excite1 = true;
    lts.excite2 = mode == gra::MReggeInclusive::DD;
    lts.t = lts.t1 = (lower + upper) / 2.0 - 1.0e-4;
    const double reference = std::norm(regge.ME2(lts, mode));
    REQUIRE(reference > 0.0);
    for (const double offset : {-1.0e-6, 0.0, 1.0e-6, 1.0e-4}) {
      lts.t = lts.t1 = (lower + upper) / 2.0 + offset;
      const double amplitude2 = std::norm(regge.ME2(lts, mode));
      REQUIRE(std::isfinite(amplitude2));
      CHECK(amplitude2 / reference == Approx(1.0).epsilon(0.02));
    }
  }
}

// Reject cut sampling before the eikonal probability tables are initialized
TEST_CASE("MEikonal cut sampling rejects uninitialized probability tables", "[gra::MEikonal][validation]") {
  gra::MEikonal eikonal;
  REQUIRE(eikonal.SoftModelHandle() == nullptr);
  std::ranlux48 random(1234);
  unsigned int  cut              = 0;
  double        impact_parameter = 0.0;
  REQUIRE_THROWS_AS(eikonal.S3GetRandomCuts(random), std::invalid_argument);
  REQUIRE_THROWS_AS(eikonal.S3GetRandomCutsBt(cut, impact_parameter, random), std::invalid_argument);
  REQUIRE_THROWS(eikonal.S3Constructor(25.0, ProtonInitialState()));
  REQUIRE_THROWS(gra::MEikonal(gra::MModelTunePtr{}));
}

// Compare both cut samplers with integrated spectral Poisson moments on the same impact grid
TEST_CASE("MEikonal joint cut sampling reproduces its integrated measure", "[gra::MEikonal][CutSpectrum][sampling]") {
  for (const bool log_b : {false, true}) {
    CAPTURE(log_b);
    // Use a coarse grid to resolve the difference between node counting and the integrated measure
    gra::MEikonal eikonal(EikonalRegressionTune(log_b ? "cuts_log" : "cuts_linear", 1, 0.0, log_b, 2));
    eikonal.S3Constructor(25.0, ProtonInitialState());
    const auto &b = eikonal.GetMatrixRuntime().ImpactParameterNodes();
    std::array<std::vector<double>, 5> integrand;
    for (auto &row : integrand) { row.assign(b.size(), 0.0); }
    for (const auto &node : indices(b)) {
      const auto spectrum = eikonal.GetMatrixRuntime().CutSpectrum(b[node]);
      const double measure = 2.0 * gra::math::PI * (log_b ? b[node] * b[node] : b[node]);
      for (const auto &mode : indices(spectrum.eigenvalues)) {
        const double opacity = spectrum.eigenvalues[mode];
        double probability = std::exp(-opacity) * spectrum.incoming_weights[mode];
        for (unsigned int m = 1; m < 25; ++m) {
          probability *= opacity / m;
          const double weight = measure * probability;
          integrand[0][node] += weight;
          integrand[1][node] += weight * m;
          integrand[2][node] += weight * m * m;
          integrand[3][node] += weight * b[node];
          integrand[4][node] += weight * b[node] * b[node];
        }
      }
    }
    const double step = (log_b ? std::log(b.back() / b.front()) : b.back() - b.front()) / (b.size() - 1);
    const double norm = gra::math::CS13Integral(integrand[0], step);
    REQUIRE(norm > 0.0);
    const double mean_m = gra::math::CS13Integral(integrand[1], step) / norm;
    const double var_m = gra::math::CS13Integral(integrand[2], step) / norm - mean_m * mean_m;
    const double mean_b = gra::math::CS13Integral(integrand[3], step) / norm;
    const double var_b = gra::math::CS13Integral(integrand[4], step) / norm - mean_b * mean_b;
    std::mt19937 random(2026);
    constexpr unsigned int samples = 100000;
    double joint_m = 0.0;
    double marginal_m = 0.0;
    double joint_b = 0.0;
    for (unsigned int sample = 0; sample < samples; ++sample) {
      unsigned int m = 0;
      double bt = 0.0;
      eikonal.S3GetRandomCutsBt(m, bt, random);
      joint_m += m;
      joint_b += bt;
      marginal_m += eikonal.S3GetRandomCuts(random);
    }
    const double tolerance_m = 8.0 * std::sqrt(std::max(0.0, var_m) / samples) + 1.0 / samples;
    const double tolerance_b = 8.0 * std::sqrt(std::max(0.0, var_b) / samples) + 1.0 / samples;
    REQUIRE(std::abs(joint_m / samples - mean_m) < tolerance_m);
    REQUIRE(std::abs(marginal_m / samples - mean_m) < tolerance_m);
    REQUIRE(std::abs(joint_b / samples - mean_b) < tolerance_b);

    // Pending steering changes must not relabel the already constructed cut probability table
    auto reference_random = random;
    std::array<std::pair<unsigned int, double>, 16> reference{};
    for (auto &[m, bt] : reference) { eikonal.S3GetRandomCutsBt(m, bt, reference_random); }
    eikonal.Numerics.MinBT = 0.5;
    eikonal.Numerics.MaxBT = 3.0;
    eikonal.Numerics.NumberBT = 4;
    eikonal.Numerics.logBT = !log_b;
    for (const auto &[expected_m, expected_b] : reference) {
      unsigned int m = 0;
      double bt = 0.0;
      eikonal.S3GetRandomCutsBt(m, bt, random);
      REQUIRE(m == expected_m);
      REQUIRE(std::abs(bt - expected_b) <= std::numeric_limits<double>::epsilon() * std::abs(expected_b));
    }
    eikonal.S3Constructor(25.0, ProtonInitialState(), true);
    unsigned int m = 0;
    double bt = 0.0;
    REQUIRE_THROWS_AS(eikonal.S3GetRandomCutsBt(m, bt, random), std::invalid_argument);
    REQUIRE_THROWS_AS(eikonal.S3GetRandomCuts(random), std::invalid_argument);
  }
}

// Check the angular momentum zeros between momentum nodes against the impact profile
TEST_CASE("MEikonal spin interpolation preserves forward angular momentum", "[gra::MEikonal][helicity][physics]") {
  const auto tune = EikonalRegressionTune("spin_interpolation", 1, 0.05, false, 128);
  const auto beams = ProtonInitialState();
  constexpr double s = 25.0;
  const double flux = s * gra::kinematics::beta12(s, beams[0].mass, beams[1].mass);
  for (const double minimum : {0.0, 1.0e-12, 1.0e-8}) {
    gra::MEikonalNumerics numerics;
    numerics.ReadParameters(tune->NumericsFile(), tune->Numerics().dump());
    numerics.MinKT2 = minimum;
    numerics.MaxKT2 = 0.002;
    numerics.NumberKT2 = 2;
    numerics.logKT2 = minimum > 1.0e-9;
    const auto runtime = gra::MEikonalMatrix::Build(tune->Soft(), s, beams, numerics, true);
    const auto &b = runtime->ImpactParameterNodes();
    for (const double kt2 : {1.0e-5, 4.0e-5}) {
      const double q = std::sqrt(kt2);
      const auto amplitude = runtime->HelicityMatrix(kt2, 0, 0);
      const auto screening = runtime->ScreeningHelicityMatrix(kt2, 0, 0);
      const auto rotated = runtime->HelicityMatrix(kt2, 0, 0, 0, 0, 0.63);
      constexpr auto transitions = gra::CanonicalProtonHelicityTransitions();
      for (const auto &entry : indices(amplitude)) {
        const int harmonic = transitions[entry].azimuth_harmonic;
        if (harmonic == 0) { continue; }
        const int n = std::abs(harmonic);
        std::vector<std::complex<double>> integrand(b.size());
        for (const auto &i : indices(b)) {
          const auto profile = runtime->ImpactHelicityMatrix(b[i], 0, 0);
          integrand[i] = b[i] * gra::math::BesselJ012(b[i] * q)[n] * profile[entry];
        }
        const auto phase = n == 1 ? -gra::math::zi : std::complex<double>(-1.0, 0.0);
        const auto expected = 4.0 * gra::math::PI * flux * phase *
                              gra::math::CS13Integral(integrand, (b.back() - b.front()) / (b.size() - 1));
        const double scale = std::abs(expected);
        CAPTURE(minimum, numerics.logKT2, kt2, entry, harmonic, expected, amplitude[entry]);
        REQUIRE(scale > 1.0e-12);
        REQUIRE(std::abs(amplitude[entry] - expected) < 5.0e-4 * scale);
        REQUIRE(std::abs(screening[entry] - expected) < 5.0e-4 * scale);
        REQUIRE(std::abs(rotated[entry] - std::exp(gra::math::zi * (harmonic * 0.63)) * amplitude[entry]) <
                1.0e-12 * scale);
      }
    }
  }
}

// Distinguish rounding at the logarithmic table boundary from transfers beyond its support
TEST_CASE("MEikonal momentum bounds allow floating point roundoff", "[gra::MEikonal][validation][loop]") {
  const auto tune = EikonalRegressionTune("momentum_bounds");
  gra::MEikonalNumerics numerics;
  numerics.ReadParameters(tune->NumericsFile(), tune->Numerics().dump());
  numerics.MinKT2 = 0.001;
  numerics.MaxKT2 = 0.09;
  numerics.logKT2 = true;
  numerics.LOOP.r_max = 0.3;
  const auto runtime = gra::MEikonalMatrix::Build(tune->Soft(), 25.0, ProtonInitialState(), numerics, true);
  const double upper = runtime->MomentumTransferNodes().back();
  const auto expected = runtime->HelicityMatrix(upper, 0, 0);
  for (const double kt2 : {numerics.MaxKT2, std::nextafter(upper, std::numeric_limits<double>::infinity())}) {
    CAPTURE(kt2, upper);
    const auto actual = runtime->HelicityMatrix(kt2, 0, 0);
    REQUIRE(gra::AllFinite(actual));
    for (const auto &entry : indices(actual)) {
      REQUIRE(std::abs(actual[entry] - expected[entry]) <= 1.0e-12 * std::max(1.0, std::abs(expected[entry])));
    }
    REQUIRE_NOTHROW(runtime->PairScreeningSpinBank(kt2));
  }
  REQUIRE_THROWS_AS(runtime->HelicityMatrix(upper * (1.0 + 1.0e-6), 0, 0), std::out_of_range);
}

// Reject every changed quadrature definition before reusing cached screening operators
TEST_CASE("MEikonal rejects stale loop controls with unchanged node counts", "[gra::MEikonal][params][loop]") {
  gra::MEikonal eikonal(EikonalRegressionTune("loop_controls"));
  eikonal.S3Constructor(25.0, ProtonInitialState());
  auto &param = eikonal.Numerics.LOOP;
  param.radial_integrator = "1/3";
  param.radial_intervals = 12;
  param.r_min = 0.001;
  param.r_max = 0.2;
  eikonal.InitLoopWeightMatrix();
  REQUIRE_NOTHROW(eikonal.GetLoopConst(25.0));
  SECTION("lower bound") { param.r_min = 0.002; }
  SECTION("upper bound") { param.r_max = 0.3; }
  SECTION("radial map") { param.radial_map = gra::math::RadialMap::Log; }
  SECTION("radial rule") { param.radial_integrator = "3/8"; }
  REQUIRE_THROWS_AS(eikonal.GetLoopConst(25.0), std::invalid_argument);
  REQUIRE_THROWS_AS(eikonal.GetLoopWeightMatrix(), std::invalid_argument);
  eikonal.InitLoopWeightMatrix();
  REQUIRE_NOTHROW(eikonal.GetLoopConst(25.0));
  REQUIRE_NOTHROW(eikonal.GetLoopWeightMatrix());
  REQUIRE_THROWS_AS(eikonal.GetLoopConst(std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
}

// Check initialization against the NUMERICS text captured with the SOFT model
TEST_CASE("MEikonal uses its immutable numerical snapshot", "[gra::MEikonal][params][snapshot]") {
  ModelParamRestoreGuard restore;
  const std::string      model_file = WriteSoftModelFile({}, "immutable_numerics_snapshot");
  const auto             model      = gra::MModelTune::Load(model_file);

  auto changed                                                  = model->Numerics();
  changed["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"] = "BadKT";
  const std::filesystem::path numerics_file = std::filesystem::path(model_file).parent_path() / "NUMERICS.json";
  std::ofstream               out(numerics_file);
  REQUIRE(out.good());
  out << changed.dump(2);
  out.close();

  gra::MODELPARAM = "TUNE0";
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");
  gra::MEikonal eikonal(model);
  REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState(), false, 4, 4));
  REQUIRE(eikonal.IsInitialized());
}

TEST_CASE("MEikonal extends momentum support without changing loop controls", "[gra::MEikonal][params][loop]") {
  ModelParamRestoreGuard restore;
  const std::string      model_file = WriteSoftModelFile({}, "extended_momentum_support");
  const auto             model      = gra::MModelTune::Load(model_file);
  gra::MODELPARAM                   = "TUNE0";
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");

  gra::MEikonal baseline(model);
  REQUIRE_NOTHROW(baseline.S3Constructor(25.0, ProtonInitialState(), false, 4, 4, 1.0));
  const double configured_max_kt2 = model->Numerics().at("NUMERICS_EIKONAL").at("MaxKT2");
  CHECK(baseline.Numerics.MaxKT2 == Approx(configured_max_kt2));

  gra::MEikonal extended(model);
  const double  required_max_kt2 = configured_max_kt2 + 1.0;
  REQUIRE_NOTHROW(extended.S3Constructor(25.0, ProtonInitialState(), false, 4, 4, required_max_kt2));
  CHECK(extended.Numerics.MaxKT2 == Approx(required_max_kt2));
  CHECK(extended.GetMatrixRuntime().MomentumTransferNodes().back() == Approx(required_max_kt2));
  CHECK(extended.Numerics.LOOP.r_max == Approx(baseline.Numerics.LOOP.r_max));
  const auto &loop = extended.GetLoopConst(25.0);
  REQUIRE_FALSE(loop.kt2.empty());
  CHECK(loop.kt2.back() <= Approx(gra::math::pow2(extended.Numerics.LOOP.r_max)));
  REQUIRE(loop.pair_screening_cache_ready);
}

// Check odd grid rejection and cleanup after failed reinitialization
TEST_CASE("MEikonal rejects odd impact grids without stale runtime", "[gra::MEikonal][validation][lifecycle]") {
  ModelParamRestoreGuard restore;
  auto                   eikonal = BuildTestEikonal({}, "odd_bt_reinitialize");
  const auto             model   = eikonal.SoftModelHandle();
  auto                   invalid = eikonal.Numerics;
  invalid.NumberBT               = 5;

  REQUIRE_THROWS_AS(gra::MEikonalMatrix::Build(model, 25.0, ProtonInitialState(), invalid, false),
                    std::invalid_argument);
  REQUIRE(eikonal.IsInitialized());
  REQUIRE_NOTHROW(eikonal.GetMatrixRuntime());

  REQUIRE_THROWS_AS(eikonal.S3Constructor(25.0, ProtonInitialState(), false, 5, 4), std::invalid_argument);
  REQUIRE_FALSE(eikonal.IsInitialized());
  REQUIRE_THROWS_AS(eikonal.GetMatrixRuntime(), std::logic_error);
  REQUIRE(eikonal.GetChannelCount() == 0);
  REQUIRE(gra::math::IsZero(eikonal.InitializedMandelstamS()));
  REQUIRE(eikonal.InitialState().empty());
  REQUIRE(gra::math::IsZero(eikonal.MinOpacityEigenvalue()));
  REQUIRE(eikonal.GetExclusiveDiffXS().size_row() == 0);
  double total     = -1.0;
  double elastic   = -1.0;
  double inelastic = -1.0;
  eikonal.GetTotXS(total, elastic, inelastic);
  REQUIRE(gra::math::IsZero(total));
  REQUIRE(gra::math::IsZero(elastic));
  REQUIRE(gra::math::IsZero(inelastic));
  REQUIRE_THROWS(eikonal.GetLoopConst(25.0));
  std::ranlux48 random(1234);
  REQUIRE_THROWS(eikonal.S3GetRandomCuts(random));

  REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState(), false, 4, 4));
  REQUIRE(eikonal.IsInitialized());
}

TEST_CASE("Cylinder fragmentation rejects unsupported low multiplicity", "[gra::MFragment][validation]") {
  const auto       model_tune = gra::MModelTune::Load(modelfile);
  gra::MModelCache cache(model_tune);
  const auto       param = gra::GetNstarParam(cache);
  REQUIRE(param == gra::GetNstarParam(cache));
  gra::MRandom        random;
  std::vector<double> mass;
  std::vector<int>    pdg;
  REQUIRE_FALSE(gra::MFragment::PickParticles(2.0, 1, 0, 0, 0, mass, pdg, LoadedPDGTable(), random, *param));

  const gra::M4Vec        mother(0.0, 0.0, 0.0, 2.0);
  std::vector<gra::M4Vec> momenta;
  gra::MNstarParam        invalid_mode  = *param;
  invalid_mode.cylinder.pt_distribution = "EXP";
  REQUIRE_THROWS_AS(
      gra::MFragment::TubeFragment(mother, 2.0, {0.13957}, momenta, 1.0, 0.2, 0.5, 1.0, random, invalid_mode),
      std::invalid_argument);

  REQUIRE(gra::MFragment::TubeFragment(mother, 0.2, {0.13957, 0.13957}, momenta, 1.1, 0.2, 0.5, 1.0, random, *param) <
          0.0);

  std::vector<bool> stable;
  gra::MFragment::GetDecayStatus({111, 311, -311, 211}, stable);
  REQUIRE(stable == std::vector<bool>{false, false, false, true});

  std::vector<double> invalid_pt = {1.0};
  REQUIRE_FALSE(gra::MFragment::ExpPowRND(1.0, 0.2, 1.0, {gra::PDG::mpi}, invalid_pt, random, *param));
  REQUIRE(invalid_pt == std::vector<double>{0.0});

  gra::MNstarParam exhausted_pt   = *param;
  exhausted_pt.cylinder.pt_bins   = 1;
  exhausted_pt.cylinder.pt_trials = 1;
  REQUIRE_FALSE(gra::MFragment::ExpPowRND(1.1, 0.2, 1.0, {gra::PDG::mpi}, invalid_pt, random, exhausted_pt));

  std::vector<int> nstar_decay = {211};
  gra::MFragment::NstarDecayTable(1, gra::PDG::mp, nstar_decay, random, *param);
  REQUIRE(nstar_decay.empty());

  gra::MNstarParam antimatter       = *param;
  antimatter.cylinder.charged_prob  = 1.0e-12;
  antimatter.cylinder.particle_prob = 0.0;
  antimatter.cylinder.species_ratio = {0.0, 0.0, 1.0};
  REQUIRE(gra::MFragment::PickParticles(4.0, 2, -2, 0, 0, mass, pdg, LoadedPDGTable(), random, antimatter));
  REQUIRE(pdg == std::vector<int>{-gra::PDG::PDG_n, -gra::PDG::PDG_n});
}

// Reject empty hadronic phase space before drawing species, while retaining physical threshold states
TEST_CASE("Cylinder species sampling respects the baryon mass threshold", "[gra::MFragment][threshold]") {
  gra::MModelCache cache(gra::MModelTune::Load(modelfile));
  const auto param = gra::GetNstarParam(cache);
  const auto table = LoadedPDGTable();
  const double nucleon = table.FindByPDG(2212).mass;
  const double pion = table.FindByPDG(111).mass;
  gra::MRandom random;
  random.SetSeed(7813);
  std::vector<double> masses;
  std::vector<int> pdgs;
  for (const unsigned int multiplicity : {2, 4, 6}) {
    for (const int baryon : {-1, 0, 1}) {
      const double threshold = std::abs(baryon) * nucleon + (multiplicity - std::abs(baryon)) * pion;
      const auto before = random.rng;
      REQUIRE_FALSE(gra::MFragment::PickParticles(threshold - 1e-6, multiplicity, baryon, 0, baryon,
                                                  masses, pdgs, table, random, *param));
      REQUIRE(random.rng == before);
    }
  }
  REQUIRE(gra::MFragment::PickParticles(nucleon + pion + 1e-6, 2, 1, 0, 1,
                                       masses, pdgs, table, random, *param));
  std::sort(pdgs.begin(), pdgs.end());
  REQUIRE(pdgs == std::vector<int>{111, 2212});
  REQUIRE(gra::Sum(masses) == Approx(nucleon + pion).epsilon(1e-14));
}

TEST_CASE("MEikonal Pomeron rotating eta mode", "[gra::MEikonal]") {
  const auto  tune              = WriteModifiedPhotoVMTune("eikonal_rotating_eta", [](auto &j) {
    const std::string model                                      = j["PARAM_SOFT"]["active_model"];
    j["PARAM_SOFT"]["MODEL"][model]["EXCHANGE"]["P"]["eta_mode"] = "rotating";
  });
  const auto  model_tune        = gra::MModelTune::Load(tune.second);
  const auto  soft_model        = model_tune->Soft();
  const auto  pomeron_id        = soft_model->ExchangeId("P");
  const auto &pomeron_parameter = soft_model->Exchange(pomeron_id);
  REQUIRE(pomeron_parameter.Eta() == gra::EtaMode::Rotating);

  MEikonal eikonal(model_tune);
  REQUIRE(eikonal.SoftModelHandle() == soft_model);
  const double s     = 25.0;
  const double t     = -0.1;
  const double alpha = soft_model->Alpha(pomeron_id, t);
  const double F     = soft_model->FormFactor(pomeron_id, t, 0);

  const std::complex<double> expected = gra::regge::EtaPhase(alpha, gra::regge::Signature::Positive) *
                                        pomeron_parameter.Coupling(0, 0) * pomeron_parameter.Coupling(0, 0) * F * F *
                                        std::pow(s, alpha);
  const std::complex<double> actual = gra::MEikonal::SingleAmpElastic(soft_model, s, t, pomeron_id, 0, 0);

  REQUIRE(std::real(actual) == Approx(std::real(expected)).epsilon(1e-12));
  REQUIRE(std::imag(actual) == Approx(std::imag(expected)).epsilon(1e-12));
}

TEST_CASE("MEikonal Pomeron and Odderon eta modes are independent", "[gra::MEikonal][physics]") {
  const auto  tune              = WriteModifiedPhotoVMTune("eikonal_independent_eta", [](auto &j) {
    const std::string model    = j["PARAM_SOFT"]["active_model"];
    auto             &exchange = j["PARAM_SOFT"]["MODEL"][model]["EXCHANGE"];
    exchange["P"]["eta_mode"]  = "raw";
    exchange["O"]["eta_mode"]  = "rotating";
    exchange["O"]["g"][0][0]   = 1.7;
    exchange["O"]["sign"]      = -1;
  });
  const auto  model_tune        = gra::MModelTune::Load(tune.second);
  const auto  soft_model        = model_tune->Soft();
  const auto  pomeron_id        = soft_model->ExchangeId("P");
  const auto  odderon_id        = soft_model->ExchangeId("O");
  const auto &pomeron_parameter = soft_model->Exchange(pomeron_id);
  const auto &odderon_parameter = soft_model->Exchange(odderon_id);
  REQUIRE(pomeron_parameter.Eta() == gra::EtaMode::Raw);
  REQUIRE(odderon_parameter.Eta() == gra::EtaMode::Rotating);

  MEikonal eikonal(model_tune);
  REQUIRE(eikonal.SoftModelHandle() == soft_model);
  const double s       = 25.0;
  const double t       = -0.1;
  const double alpha_P = soft_model->Alpha(pomeron_id, t);
  const double alpha_O = soft_model->Alpha(odderon_id, t);
  const double F_P     = soft_model->FormFactor(pomeron_id, t, 0);
  const double F_O     = soft_model->FormFactor(odderon_id, t, 0);

  const std::complex<double> expected_P =
      gra::regge::EtaRaw(alpha_P, gra::regge::Signature::Positive, gra::regge::Rim::Lower) *
      pomeron_parameter.Coupling(0, 0) * pomeron_parameter.Coupling(0, 0) * F_P * F_P * std::pow(s, alpha_P);
  const std::complex<double> expected_O = static_cast<double>(odderon_parameter.ResidueSign()) *
                                          gra::regge::EtaPhase(alpha_O, gra::regge::Signature::Negative) *
                                          odderon_parameter.Coupling(0, 0) * odderon_parameter.Coupling(0, 0) * F_O *
                                          F_O * std::pow(s, alpha_O);
  const std::complex<double> actual_P = gra::MEikonal::SingleAmpElastic(soft_model, s, t, pomeron_id, 0, 0);
  const std::complex<double> actual_O = gra::MEikonal::SingleAmpElastic(soft_model, s, t, odderon_id, 0, 0);

  REQUIRE(std::real(actual_P) == Approx(std::real(expected_P)).epsilon(1e-12));
  REQUIRE(std::imag(actual_P) == Approx(std::imag(expected_P)).epsilon(1e-12));
  REQUIRE(std::real(actual_O) == Approx(std::real(expected_O)).epsilon(1e-12));
  REQUIRE(std::imag(actual_O) == Approx(std::imag(expected_O)).epsilon(1e-12));
}

TEST_CASE("MEikonal unitarization preserves the complex impact profile", "[gra::MEikonal][eikonal][phase][physics]") {
  struct UnitarizationCase {
    std::string name;
    double      q = 1.0;
  };
  const std::array<UnitarizationCase, 2> cases = {{{"exp", 1.0}, {"q_exp", 0.63}}};

  for (const auto &item : cases) {
    DYNAMIC_SECTION(item.name) {
      const auto  eikonal = BuildTestEikonal({}, "complex_profile_" + item.name, 32, 8, 0.0, item.name, item.q);
      const auto &runtime = eikonal.GetMatrixRuntime();
      const auto &b_node  = runtime.ImpactParameterNodes();
      REQUIRE(b_node.size() > 4);

      const std::array<std::size_t, 3> probe           = {b_node.size() / 4, b_node.size() / 2, 3 * b_node.size() / 4};
      double                           amplitude_scale = 0.0;
      double                           phase_distance  = 0.0;
      for (const std::size_t index : probe) {
        const std::complex<double> chi = eikonal.S3Density(b_node[index], 0, 0);
        const std::complex<double> smatrix =
            item.name == "exp"
                ? std::exp(gra::math::zi * chi)
                : std::pow(std::complex<double>(1.0, 0.0) + (1.0 - item.q) * gra::math::zi * chi, 1.0 / (1.0 - item.q));
        const std::complex<double> expected = gra::math::zi * (std::complex<double>(1.0, 0.0) - smatrix);
        const std::complex<double> actual   = runtime.ImpactHelicityMatrix(b_node[index], 0, 0)[0];
        const double               scale    = std::max(1.0, std::abs(expected));
        CAPTURE(item.name, index, b_node[index], chi, expected, actual);
        REQUIRE(std::abs(actual - expected) < 5.0e-12 * scale);
        amplitude_scale = std::max(amplitude_scale, std::abs(expected));
        phase_distance  = std::max(phase_distance, std::abs(expected - std::complex<double>(std::abs(expected), 0.0)));
      }
      REQUIRE(amplitude_scale > 0.0);
      REQUIRE(phase_distance > 1.0e-3 * amplitude_scale);
    }
  }
}

TEST_CASE("MEikonal Fourier-Bessel transform preserves complex helicity harmonics",
          "[gra::MEikonal][eikonal][fourier][phase][physics]") {
  const auto    tune  = WriteModifiedPhotoVMTune("eikonal_complex_hankel", [](auto &j) {
    const std::string model            = j.at("PARAM_SOFT").at("active_model");
    auto             &block            = j.at("PARAM_SOFT").at("MODEL").at(model);
    block.at("EIKONAL").at("helicity") = true;
    auto &exchange                     = block.at("EXCHANGE");
    for (auto it = exchange.begin(); it != exchange.end(); ++it) { it.value().at("on") = it.key() == "P"; }
    exchange.at("P").at("helicity").at("kappa") = 0.65;
  });
  const auto    model = gra::MModelTune::Load(tune.second);
  gra::MEikonal eikonal(model);
  REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState(), false, 256, 64));

  const auto &runtime = eikonal.GetMatrixRuntime();
  const auto &q2_node = runtime.MomentumTransferNodes();
  const auto &b_node  = runtime.ImpactParameterNodes();
  REQUIRE(q2_node.size() > 2);
  REQUIRE(b_node.size() == 257);
  const double q2   = q2_node[1];
  const double q    = std::sqrt(q2);
  const double beta = gra::kinematics::beta12(25.0, ProtonInitialState()[0].mass, ProtonInitialState()[1].mass);

  // Integrate the public impact profile with an independent Gaussian rule
  const auto transform = [&](const std::size_t order) {
    const auto [node, weight] = gra::math::GaussLegendreRule(order, b_node.front(), b_node.back());
    gra::ProtonHelicityMatrix integral{};
    for (const auto &i : indices(node)) {
      const auto profile = runtime.ImpactHelicityMatrix(node[i], 0, 0);
      const auto bessel  = gra::math::BesselJ012(node[i] * q);
      for (const auto &entry : indices(integral)) {
        const int         harmonic_signed = gra::CanonicalProtonHelicityTransitions()[entry].azimuth_harmonic;
        const std::size_t harmonic = static_cast<std::size_t>(harmonic_signed < 0 ? -harmonic_signed : harmonic_signed);
        integral[entry] += weight[i] * node[i] * bessel[harmonic] * profile[entry];
      }
    }
    for (const auto &entry : indices(integral)) {
      const int         harmonic_signed = gra::CanonicalProtonHelicityTransitions()[entry].azimuth_harmonic;
      const std::size_t harmonic = static_cast<std::size_t>(harmonic_signed < 0 ? -harmonic_signed : harmonic_signed);
      std::complex<double> phase = 1.0;
      for (std::size_t power = 0; power < harmonic; ++power) { phase *= -gra::math::zi; }
      integral[entry] *= 4.0 * gra::math::PI * 25.0 * beta * phase;
    }
    return integral;
  };

  const auto                       coarse    = transform(256);
  const auto                       fine      = transform(512);
  const auto                       full      = runtime.HelicityMatrix(q2, 0, 0);
  const auto                       screening = runtime.ScreeningHelicityMatrix(q2, 0, 0);
  const std::array<std::size_t, 3> entry     = {gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi1),
                                                gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi5),
                                                gra::ElasticHelicityReferenceIndex(gra::ElasticHelicityComponent::Phi4)};
  for (const std::size_t h : entry) {
    const int    harmonic_signed  = gra::CanonicalProtonHelicityTransitions()[h].azimuth_harmonic;
    const int    harmonic         = harmonic_signed < 0 ? -harmonic_signed : harmonic_signed;
    const double scale            = std::max({std::abs(fine[h]), std::abs(full[h]), std::abs(screening[h]), 1.0e-14});
    const double quadrature_error = std::abs(fine[h] - coarse[h]);
    const double tolerance        = 2.0e-3 * scale + 4.0 * quadrature_error + 1.0e-13;
    CAPTURE(h, harmonic, q2, coarse[h], fine[h], full[h], screening[h], quadrature_error, tolerance);
    REQUIRE(std::abs(fine[h]) > 1.0e-12);
    REQUIRE(quadrature_error < 5.0e-3 * scale);
    REQUIRE(std::abs(full[h] - fine[h]) < tolerance);
    REQUIRE(std::abs(screening[h] - fine[h]) < tolerance);
    REQUIRE(std::abs(full[h] - screening[h]) < 1.0e-11 * scale + 1.0e-13);
  }
}

TEST_CASE("PARAM_SOFT rejects invalid eta steering", "[gra::MEikonal][params][validation]") {
  SECTION("invalid exchange signature") {
    const auto tune = WriteModifiedPhotoVMTune(
        "soft_invalid_tau", [](auto &j) { j.at("PARAM_SOFT").at("EXCHANGE_DEF").at("O").at("tau") = 0; });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("unknown exchange name") {
    const auto tune = WriteModifiedPhotoVMTune("soft_unknown_exchange", [](auto &j) {
      const std::string model                                      = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE")["X"] = nlohmann::json::object();
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("missing eta array is rejected") {
    const auto tune = WriteModifiedPhotoVMTune("soft_missing_eta_mode", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE").at("P").erase("eta_mode");
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("unknown Pomeron mode") {
    const auto tune = WriteModifiedPhotoVMTune("soft_unknown_pomeron_eta_mode", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE").at("P").at("eta_mode") = "unknown";
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("removed full t0 mode") {
    const auto tune = WriteModifiedPhotoVMTune("soft_removed_full_t0", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE").at("P").at("eta_mode") = "full_t0";
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("unknown Odderon mode") {
    const auto tune = WriteModifiedPhotoVMTune("soft_unknown_odderon_eta_mode", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EXCHANGE").at("O").at("eta_mode") = "unknown";
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("missing triple-Pomeron mode") {
    const auto tune = WriteModifiedPhotoVMTune("soft_missing_triple_pomeron_eta_mode", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("3P").erase("eta_mode");
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("unknown triple-Pomeron mode") {
    const auto tune = WriteModifiedPhotoVMTune("soft_unknown_triple_pomeron_eta_mode", [](auto &j) {
      const std::string model                                          = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("3P").at("eta_mode") = "unknown";
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }
}

TEST_CASE("PARAM_SOFT named form factor banks support sharing and ragged rows", "[gra::MForm][gra::MEikonal][params]") {
  SECTION("f2 and rho share one bank") {
    const auto  model_tune = gra::MModelTune::Load(modelfile);
    const auto &soft_model = *model_tune->Soft();
    const auto &f2         = soft_model.Exchange(soft_model.ExchangeId("R_f2"));
    const auto &rho        = soft_model.Exchange(soft_model.ExchangeId("R_rho"));
    REQUIRE(f2.FormFactorName() == "R");
    REQUIRE(rho.FormFactorName() == "R");
    REQUIRE(f2.FormFactorType() == gra::SoftFormFactor::Exponential);
    REQUIRE(f2.FormFactorParameters() == rho.FormFactorParameters());
    REQUIRE(f2.FormFactorParameters().at(0).size() == 1);
    const double slope = f2.FormFactorParameters().at(0).at(0);
    const double t = -0.3;
    REQUIRE(soft_model.FormFactor(soft_model.ExchangeId("R_f2"), t, 0) ==
            Approx(std::exp(0.5 * slope * t)).epsilon(1e-12));
  }

  SECTION("Pomeron and Odderon select rows with different parameter counts") {
    auto              j             = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    const std::string model         = "double";
    j["PARAM_SOFT"]["active_model"] = model;
    auto &block                     = j["PARAM_SOFT"]["MODEL"][model];
    block["FF"]["O_GK"]             = {
                    {"type", "GKERNEL"},
                    {"param", {{1.0, 1.0, 0.0, 1.0, 1.2, 1.1, 0.0, 1.0}, {1.3, 1.2, 0.0, 1.0, 1.5, 1.3, 0.0, 1.0}}}};
    block["EXCHANGE"]["O"]["ff"] = "O_GK";

    const std::string path = WriteSoftModelFile({}, "ragged_form_factor_test");
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    const auto  model_tune = gra::MModelTune::Load(path);
    const auto &soft_model = *model_tune->Soft();
    const auto  pomeron_id = soft_model.ExchangeId("P");
    const auto  odderon_id = soft_model.ExchangeId("O");
    const auto &pomeron    = soft_model.Exchange(pomeron_id);
    const auto &odderon    = soft_model.Exchange(odderon_id);
    REQUIRE(pomeron.FormFactorType() == gra::SoftFormFactor::DualPower);
    REQUIRE(odderon.FormFactorType() == gra::SoftFormFactor::GeneralizedKernel);
    REQUIRE(pomeron.FormFactorParameters().at(0).size() == 3);
    REQUIRE(odderon.FormFactorParameters().at(0).size() == 8);
    REQUIRE(std::isfinite(soft_model.FormFactor(pomeron_id, -0.1, 0)));
    REQUIRE(std::isfinite(soft_model.FormFactor(odderon_id, -0.1, 0)));
  }
}

TEST_CASE("PARAM_SOFT supports generic even and odd soft Regge exchanges",
          "[gra::MForm][gra::MEikonal][params][physics]") {
  auto              j             = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const std::string model         = "single";
  auto             &soft          = j["PARAM_SOFT"];
  soft["active_model"]            = model;
  auto &definition                = soft["EXCHANGE_DEF"];
  definition["R_f2"]["tau"]       = 1;
  definition["R_f2"]["crossing"]  = 1;
  definition["R_rho"]["tau"]      = -1;
  definition["R_rho"]["crossing"] = -1;
  auto &exchange                  = soft["MODEL"][model]["EXCHANGE"];
  exchange["R_f2"]["alpha"]       = {0.55, 0.90};
  exchange["R_f2"]["g"]           = {{1.2}};
  exchange["R_f2"]["sign"]        = 1;
  exchange["R_rho"]["alpha"]      = {0.45, 0.90};
  exchange["R_rho"]["g"]          = {{0.8}};
  exchange["R_rho"]["sign"]       = -1;

  const std::string path = "tmp/graniitti_soft_generic_exchange_test.json";
  std::ofstream     out(path);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  const auto soft_model = gra::SoftModel::LoadFromJson(path, gra::aux::GetInputData(path));
  const auto f2         = soft_model->ExchangeId("R_f2");
  const auto pomeron    = soft_model->ExchangeId("P");
  const auto rho        = soft_model->ExchangeId("R_rho");
  REQUIRE(f2.Value() != pomeron.Value());
  REQUIRE(pomeron.Value() != rho.Value());
  REQUIRE(f2.Value() != rho.Value());

  const double t = -0.2;
  REQUIRE(soft_model->Alpha(f2, t) == Approx(0.55 + 0.90 * t));
  REQUIRE(soft_model->Alpha(rho, t) == Approx(0.45 + 0.90 * t));
  const auto &pomeron_parameter = soft_model->Exchange(pomeron);
  REQUIRE(std::abs(soft_model->Alpha(pomeron, t) - (pomeron_parameter.Alpha0() + pomeron_parameter.AlphaPrime() * t)) >
          1e-12);

  const double s = 25.0;
  for (const auto exchange : {f2, rho}) {
    const auto                &parameter   = soft_model->Exchange(exchange);
    const double               alpha       = soft_model->Alpha(exchange, t);
    const double               form_factor = soft_model->FormFactor(exchange, t, 0);
    const std::complex<double> expected =
        static_cast<double>(parameter.ResidueSign()) *
        gra::regge::EtaFactor(alpha, parameter.Alpha0(), parameter.Signature(), parameter.Eta()) *
        parameter.Coupling(0, 0) * parameter.Coupling(0, 0) * form_factor * form_factor * std::pow(s, alpha);
    const std::complex<double> actual = gra::MEikonal::SingleAmpElastic(soft_model, s, t, exchange, 0, 0);
    REQUIRE(std::real(actual) == Approx(std::real(expected)).epsilon(1e-12));
    REQUIRE(std::imag(actual) == Approx(std::imag(expected)).epsilon(1e-12));
  }
}

TEST_CASE("MMatrix real mixing preserves low-dimensional conventions", "[gra::MMatrix][gra::MEikonal]") {
  SECTION("one channel") {
    const auto actual = gra::MMatrix<double>::MixingReal({}, 1);
    REQUIRE(actual.size_row() == 1);
    REQUIRE(actual.size_col() == 1);
    REQUIRE(actual[0][0] == Approx(1.0));
  }

  SECTION("two channels") {
    const double               theta  = 0.31;
    const double               cosine = std::cos(theta);
    const double               sine   = std::sin(theta);
    const gra::MMatrix<double> expected{{cosine, sine}, {-sine, cosine}};
    const auto                 actual = gra::MMatrix<double>::MixingReal({theta}, 2);
    for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) { REQUIRE(actual[i][j] == Approx(expected[i][j]).epsilon(1e-14)); }
    }
  }

  SECTION("three channels") {
    const std::vector<double>  theta = {0.17, -0.11, 0.23};
    const double               c12   = std::cos(theta[0]);
    const double               s12   = std::sin(theta[0]);
    const double               c13   = std::cos(theta[1]);
    const double               s13   = std::sin(theta[1]);
    const double               c23   = std::cos(theta[2]);
    const double               s23   = std::sin(theta[2]);
    const gra::MMatrix<double> expected{{c12 * c13, s12 * c13, s13},
                                        {-s12 * c23 - c12 * s23 * s13, c12 * c23 - s12 * s23 * s13, s23 * c13},
                                        {s12 * s23 - c12 * c23 * s13, -c12 * s23 - s12 * c23 * s13, c23 * c13}};
    const auto                 actual = gra::MMatrix<double>::MixingReal(theta, 3);
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) { REQUIRE(actual[i][j] == Approx(expected[i][j]).epsilon(1e-14)); }
    }
  }

  SECTION("four channels are orthogonal") {
    const auto actual = gra::MMatrix<double>::MixingReal({0.11, -0.07, 0.19, 0.23, -0.13, 0.29}, 4);
    for (std::size_t i = 0; i < 4; ++i) {
      for (std::size_t j = 0; j < 4; ++j) {
        double inner = 0.0;
        for (std::size_t k = 0; k < 4; ++k) { inner += actual[i][k] * actual[j][k]; }
        REQUIRE(inner == Approx(i == j ? 1.0 : 0.0).margin(1e-14));
      }
    }
  }

  SECTION("invalid angle arrays are rejected") {
    REQUIRE_THROWS_AS(gra::MMatrix<double>::MixingReal({}, 0), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::MMatrix<double>::MixingReal({0.1, 0.2}, 3), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::MMatrix<double>::MixingReal({0.1, std::numeric_limits<double>::infinity(), 0.2}, 3),
                      std::invalid_argument);
  }
}

TEST_CASE("PARAM_SOFT validates the generic angle count", "[gra::PARAM_SOFT][gra::MEikonal]") {
  const std::string valid_path = WriteSoftModelFile({0.17, -0.11, 0.23}, "invalid_angle_count_source");
  auto              j          = nlohmann::json::parse(gra::aux::GetInputData(valid_path));
  j["PARAM_SOFT"]["MODEL"]["test_multichannel"]["GW"]["theta"] = {0.17, -0.11};

  const std::string invalid_path =
      (std::filesystem::path(valid_path).parent_path() / "GENERAL_invalid_angle_count.json").string();
  std::ofstream out(invalid_path);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  REQUIRE_THROWS(gra::MModelTune::Load(invalid_path));
}

TEST_CASE("PARAM_SOFT validates the resolved Good Walker direction", "[gra::PARAM_SOFT][gra::MEikonal]") {
  const std::string valid_path = WriteSoftModelFile({0.17, -0.11, 0.23}, "invalid_a_c_source");
  const auto        valid      = nlohmann::json::parse(gra::aux::GetInputData(valid_path));

  SECTION("channel count") {
    auto j                                                     = valid;
    j["PARAM_SOFT"]["MODEL"]["test_multichannel"]["GW"]["a_c"] = {1.0};
    const std::string path =
        (std::filesystem::path(valid_path).parent_path() / "GENERAL_invalid_a_c_count.json").string();
    std::ofstream out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    REQUIRE_THROWS(gra::MModelTune::Load(path));
  }

  SECTION("the resolved direction is scale invariant") {
    auto j                                                     = valid;
    j["PARAM_SOFT"]["MODEL"]["test_multichannel"]["GW"]["a_c"] = {0.5, 0.5};
    const std::string path = (std::filesystem::path(valid_path).parent_path() / "GENERAL_scaled_a_c.json").string();
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    const auto  scaled       = gra::MModelTune::Load(path);
    const auto &coefficients = scaled->Soft()->GoodWalker().ResolvedCoefficients();
    REQUIRE(coefficients.size() == 2);
    REQUIRE(coefficients[0] == Approx(1.0 / std::sqrt(2.0)));
    REQUIRE(coefficients[1] == Approx(1.0 / std::sqrt(2.0)));
  }

  SECTION("signed amplitude coordinates reach the screening cache") {
    auto j                                                     = valid;
    j["PARAM_SOFT"]["MODEL"]["test_multichannel"]["GW"]["a_c"] = {1.0, -2.0};
    const std::string path = (std::filesystem::path(valid_path).parent_path() / "GENERAL_signed_a_c.json").string();
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    const auto  signed_model = gra::MModelTune::Load(path);
    const auto &coefficient  = signed_model->Soft()->GoodWalker().ResolvedCoefficients();
    REQUIRE(coefficient.size() == 2);
    REQUIRE(coefficient[0] == Approx(1.0 / std::sqrt(5.0)).margin(1.0e-15));
    REQUIRE(coefficient[1] == Approx(-2.0 / std::sqrt(5.0)).margin(1.0e-15));

    ModelParamRestoreGuard restore;
    gra::MODELPARAM = "TUNE0";
    std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");
    gra::MEikonal eikonal(signed_model);
    REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState(), false, 4, 4));
    const auto &channel = eikonal.GetLoopConst(25.0)
                              .good_walker_channels[static_cast<std::size_t>(gra::GWFinalClass::SingleDissociation1)];
    REQUIRE(channel.size() == coefficient.size());
    for (const auto &i : indices(channel)) {
      REQUIRE(channel[i].born_coefficient == Approx(coefficient[i]).margin(1.0e-15));
    }
    REQUIRE(channel[0].born_coefficient > 0.0);
    REQUIRE(channel[1].born_coefficient < 0.0);
  }

  SECTION("zero resolved direction") {
    auto j                                                     = valid;
    j["PARAM_SOFT"]["MODEL"]["test_multichannel"]["GW"]["a_c"] = {0.0, 0.0};
    const std::string path = (std::filesystem::path(valid_path).parent_path() / "GENERAL_zero_a_c.json").string();
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    REQUIRE_THROWS(gra::MModelTune::Load(path));
  }
}

TEST_CASE("PARAM_SOFT validates screening exchange selection", "[gra::PARAM_SOFT][gra::MEikonal][validation]") {
  SECTION("missing selection") {
    const auto tune = WriteModifiedPhotoVMTune("soft_missing_screening_exchanges", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EIKONAL").erase("screening_exchanges");
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("unknown exchange") {
    const auto tune = WriteModifiedPhotoVMTune("soft_unknown_screening_exchange", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EIKONAL").at("screening_exchanges") = {"X"};
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }

  SECTION("duplicate exchange") {
    const auto tune = WriteModifiedPhotoVMTune("soft_duplicate_screening_exchange", [](auto &j) {
      const std::string model = j.at("PARAM_SOFT").at("active_model");
      j.at("PARAM_SOFT").at("MODEL").at(model).at("EIKONAL").at("screening_exchanges") = {"P", "P"};
    });
    REQUIRE_THROWS(gra::MModelTune::Load(tune.second));
  }
}

// Check selected noncommuting Pomeron, Reggeon and Odderon screening profiles
TEST_CASE("MEikonal screening exchange selectors preserve crossing parity",
          "[gra::MEikonal][screening_exchanges][matrix]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");

  // Build one deterministic two-channel model with the requested selection
  const auto build = [](const std::vector<std::string> &selection, const std::vector<gra::MParticle> &initial_state,
                        const std::string &suffix) {
    const std::string model_file = WriteSoftModelFile({0.31}, "screening_selector_" + suffix);
    auto              j          = nlohmann::json::parse(gra::aux::GetInputData(model_file));
    const std::string active     = j.at("PARAM_SOFT").at("active_model");
    auto             &model      = j.at("PARAM_SOFT").at("MODEL").at(active);
    auto             &exchange   = model.at("EXCHANGE");

    exchange.at("P").at("g")                      = {{2.20, 0.00}, {0.00, 1.30}};
    exchange.at("R_f2").at("g")                   = {{1.10, 0.28}, {0.28, 0.65}};
    exchange.at("O").at("g")                      = {{0.80, -0.22}, {-0.22, 0.45}};
    exchange.at("R_rho").at("on")                 = false;
    exchange.at("O3g").at("on")                   = false;
    model.at("EIKONAL").at("screening_exchanges") = selection;

    std::ofstream out(model_file);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    const auto    soft_model = gra::MModelTune::Load(model_file);
    gra::MEikonal eikonal(soft_model);
    eikonal.S3Constructor(25.0, initial_state, false, 4, 4);
    return eikonal;
  };

  const std::vector<std::string> p        = {"P"};
  const std::vector<std::string> pr       = {"P", "R_f2"};
  const std::vector<std::string> po       = {"P", "O"};
  const std::vector<std::string> all      = {"P", "R_f2", "O"};
  const std::vector<std::string> wildcard = {"*"};

  const auto p_pp        = build(p, ProtonInitialState(), "p_pp");
  const auto p_ppbar     = build(p, ProtonAntiprotonInitialState(), "p_ppbar");
  const auto pr_pp       = build(pr, ProtonInitialState(), "pr_pp");
  const auto pr_ppbar    = build(pr, ProtonAntiprotonInitialState(), "pr_ppbar");
  const auto po_pp       = build(po, ProtonInitialState(), "po_pp");
  const auto po_ppbar    = build(po, ProtonAntiprotonInitialState(), "po_ppbar");
  const auto all_pp      = build(all, ProtonInitialState(), "all_pp");
  const auto wildcard_pp = build(wildcard, ProtonInitialState(), "wildcard_pp");

  const auto  soft_model = p_pp.SoftModelHandle();
  const auto &p_matrix   = soft_model->Exchange(soft_model->ExchangeId("P")).CouplingMatrix();
  const auto &r_matrix   = soft_model->Exchange(soft_model->ExchangeId("R_f2")).CouplingMatrix();
  const auto &o_matrix   = soft_model->Exchange(soft_model->ExchangeId("O")).CouplingMatrix();
  REQUIRE((p_matrix * r_matrix - r_matrix * p_matrix).FrobNorm() > 0.1);
  REQUIRE((p_matrix * o_matrix - o_matrix * p_matrix).FrobNorm() > 0.1);
  REQUIRE((r_matrix * o_matrix - o_matrix * r_matrix).FrobNorm() > 0.1);

  const double kt2             = 0.0;
  const auto   p_pp_amp        = p_pp.S3PhysicalScreeningAmp(kt2);
  const auto   p_ppbar_amp     = p_ppbar.S3PhysicalScreeningAmp(kt2);
  const auto   pr_pp_amp       = pr_pp.S3PhysicalScreeningAmp(kt2);
  const auto   pr_ppbar_amp    = pr_ppbar.S3PhysicalScreeningAmp(kt2);
  const auto   po_pp_amp       = po_pp.S3PhysicalScreeningAmp(kt2);
  const auto   po_ppbar_amp    = po_ppbar.S3PhysicalScreeningAmp(kt2);
  const auto   all_pp_amp      = all_pp.S3PhysicalScreeningAmp(kt2);
  const auto   wildcard_pp_amp = wildcard_pp.S3PhysicalScreeningAmp(kt2);

  const auto scale = [](const std::complex<double> value) { return std::max(1.0, std::abs(value)); };
  REQUIRE(std::abs(p_pp_amp - p_ppbar_amp) < 1.0e-11 * scale(p_pp_amp));
  REQUIRE(std::abs(pr_pp_amp - pr_ppbar_amp) < 1.0e-11 * scale(pr_pp_amp));
  REQUIRE(std::abs(pr_pp_amp - p_pp_amp) > 1.0e-6 * scale(p_pp_amp));
  REQUIRE(std::abs(po_pp_amp - po_ppbar_amp) > 1.0e-6 * scale(po_pp_amp));
  REQUIRE(std::abs(all_pp_amp - wildcard_pp_amp) < 1.0e-11 * scale(all_pp_amp));
  REQUIRE(std::abs(all_pp_amp - pr_pp_amp) > 1.0e-6 * scale(all_pp_amp));
  REQUIRE(std::abs(all_pp_amp - po_pp_amp) > 1.0e-6 * scale(all_pp_amp));

  const auto all_bank      = all_pp.PairScreeningSpinBank(kt2);
  const auto wildcard_bank = wildcard_pp.PairScreeningSpinBank(kt2);
  for (const auto &entry : indices(all_bank)) {
    CAPTURE(entry);
    const double matrix_scale = std::max(1.0, all_bank[entry].FrobNorm());
    REQUIRE((all_bank[entry] - wildcard_bank[entry]).FrobNorm() < 1.0e-11 * matrix_scale);
  }
}

TEST_CASE("MEikonal exclusive amplitudes and cross sections", "[gra::MEikonal]") {
  SECTION("1-channel elastic only") {
    auto eikonal = BuildTestEikonal({}, "n1");
    REQUIRE(eikonal.SoftModelHandle() != nullptr);
    REQUIRE(std::isfinite(eikonal.MinOpacityEigenvalue()));
    REQUIRE(eikonal.MinOpacityEigenvalue() >= 0.0);
    REQUIRE(eikonal.GetMatrixRuntime().SoftModelHandle() == eikonal.SoftModelHandle());
    REQUIRE_THROWS(
        gra::MEikonalMatrix::Build(gra::SoftModelPtr{}, 25.0, ProtonInitialState(), eikonal.Numerics, false));
    const double kt2     = eikonal.Numerics.MinKT2 + 0.37 * (eikonal.Numerics.MaxKT2 - eikonal.Numerics.MinKT2);
    const double low_kt2 = 0.5 * eikonal.Numerics.MinKT2;

    const auto expected     = ExclusiveAmpReference(eikonal, kt2, 0, 0);
    const auto actual       = eikonal.S3ExclusiveAmp(kt2, 0, 0);
    const auto expected_t0  = ExclusiveAmpReference(eikonal, 0.0, 0, 0);
    const auto actual_t0    = eikonal.S3ExclusiveAmp(0.0, 0, 0);
    const auto expected_low = ExclusiveAmpReference(eikonal, low_kt2, 0, 0);
    const auto actual_low   = eikonal.S3ExclusiveAmp(low_kt2, 0, 0);

    REQUIRE(std::real(actual) == Approx(std::real(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
    REQUIRE(std::imag(actual) == Approx(std::imag(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
    REQUIRE(std::real(actual_t0) ==
            Approx(std::real(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
    REQUIRE(std::imag(actual_t0) ==
            Approx(std::imag(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
    REQUIRE(std::real(actual_low) ==
            Approx(std::real(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
    REQUIRE(std::imag(actual_low) ==
            Approx(std::imag(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::Elastic).size() == 1);
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation1).empty());
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation2).empty());
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::DoubleDissociation).empty());
    const auto &channel = eikonal.GetLoopConst(25.0).good_walker_channels;
    REQUIRE(channel[static_cast<std::size_t>(gra::GWFinalClass::Elastic)].size() == 1);
    REQUIRE(channel[static_cast<std::size_t>(gra::GWFinalClass::SingleDissociation1)].empty());
    REQUIRE(channel[static_cast<std::size_t>(gra::GWFinalClass::SingleDissociation2)].empty());
    REQUIRE(channel[static_cast<std::size_t>(gra::GWFinalClass::DoubleDissociation)].empty());
    REQUIRE(eikonal.S3DiffractiveAmpSquared(kt2, gra::GWFinalClass::Elastic) ==
            Approx(gra::math::abs2(actual)).epsilon(1e-12));
    RequirePhysicalPairScreeningProjection(eikonal, kt2);
  }

  SECTION("matrix Odderon projects pp and ppbar crossing eigenvalues") {
    auto pp    = BuildTestEikonalWithInitialState({}, "matrix_pp", ProtonInitialState(), 0.55);
    auto ppbar = BuildTestEikonalWithInitialState({}, "matrix_ppbar", ProtonAntiprotonInitialState(), 0.55);
    auto even  = BuildTestEikonalWithInitialState({}, "matrix_even", ProtonInitialState(), 0.0);
    auto pp_negative_residue =
        BuildTestEikonalWithInitialState({}, "matrix_pp_negative_residue", ProtonInitialState(), 0.55, -1);

    const double kt2 = 0.0;

    const auto pp_amp                  = pp.S3ExclusiveAmp(kt2, 0, 0);
    const auto ppbar_amp               = ppbar.S3ExclusiveAmp(kt2, 0, 0);
    const auto even_amp                = even.S3ExclusiveAmp(kt2, 0, 0);
    const auto pp_negative_residue_amp = pp_negative_residue.S3ExclusiveAmp(kt2, 0, 0);
    const auto average                 = 0.5 * (pp_amp + ppbar_amp);

    REQUIRE(std::abs(pp_amp - ppbar_amp) > 1e-10);
    REQUIRE(std::abs(average - even_amp) > 1e-10);
    REQUIRE(std::real(pp_negative_residue_amp) == Approx(std::real(ppbar_amp)).epsilon(1e-8));
    REQUIRE(std::imag(pp_negative_residue_amp) == Approx(std::imag(ppbar_amp)).epsilon(1e-8));
  }

  SECTION("screening exchange selection excludes the Odderon") {
    auto         pp    = BuildTestEikonalWithInitialState({}, "screening_pp", ProtonInitialState(), 0.55);
    auto         ppbar = BuildTestEikonalWithInitialState({}, "screening_ppbar", ProtonAntiprotonInitialState(), 0.55);
    auto         stronger = BuildTestEikonalWithInitialState({}, "screening_stronger_odd", ProtonInitialState(), 1.10);
    const double kt2      = 0.5 * (pp.Numerics.MinKT2 + pp.Numerics.MaxKT2);

    const auto pp_screening       = pp.S3PhysicalScreeningAmp(kt2);
    const auto ppbar_screening    = ppbar.S3PhysicalScreeningAmp(kt2);
    const auto stronger_screening = stronger.S3PhysicalScreeningAmp(kt2);

    REQUIRE(std::abs(pp.S3ExclusiveAmp(kt2, 0, 0) - ppbar.S3ExclusiveAmp(kt2, 0, 0)) > 1e-10);
    REQUIRE(pp_screening == pp.S3ExclusiveScreeningAmp(kt2, 0, 0));
    REQUIRE(std::real(pp_screening) == Approx(std::real(ppbar_screening)).epsilon(1e-10));
    REQUIRE(std::imag(pp_screening) == Approx(std::imag(ppbar_screening)).epsilon(1e-10));
    REQUIRE(std::real(pp_screening) == Approx(std::real(stronger_screening)).epsilon(1e-10));
    REQUIRE(std::imag(pp_screening) == Approx(std::imag(stronger_screening)).epsilon(1e-10));
  }

  SECTION("noncommuting Odderon transition reaches the elastic projection") {
    const auto diagonal = BuildTestEikonalWithInitialState({0.31}, "matrix_noncommuting_diagonal", ProtonInitialState(),
                                                           0.55, 1, 6, 6, 25.0, 0.0);
    const auto transition       = BuildTestEikonalWithInitialState({0.31}, "matrix_noncommuting_transition",
                                                                   ProtonInitialState(), 0.55, 1, 6, 6, 25.0, 0.15);
    const auto transition_ppbar = BuildTestEikonalWithInitialState(
        {0.31}, "matrix_noncommuting_ppbar", ProtonAntiprotonInitialState(), 0.55, 1, 6, 6, 25.0, 0.15);
    const double kt2 = 0.5 * (transition.Numerics.MinKT2 + transition.Numerics.MaxKT2);

    const auto helicity         = transition.GetMatrixRuntime().HelicityAmplitudes(kt2, 0, 0);
    const auto elastic          = transition.S3ExclusiveAmp(kt2, 0, 0);
    const auto diagonal_elastic = diagonal.S3ExclusiveAmp(kt2, 0, 0);
    const auto ppbar_elastic    = transition_ppbar.S3ExclusiveAmp(kt2, 0, 0);

    const auto helicity_average = 0.5 * (helicity.phi1 + helicity.phi3);
    REQUIRE(std::real(elastic) == Approx(std::real(helicity_average)));
    REQUIRE(std::imag(elastic) == Approx(std::imag(helicity_average)));
    REQUIRE(std::abs(elastic - diagonal_elastic) > 1e-10);
    REQUIRE(std::abs(elastic - ppbar_elastic) > 1e-10);
  }

  SECTION("2-channel exclusive amplitudes match the direct projection") {
    auto         eikonal = BuildTestEikonal({0.31}, "n2");
    const double kt2     = eikonal.Numerics.MinKT2 + 0.53 * (eikonal.Numerics.MaxKT2 - eikonal.Numerics.MinKT2);
    const double low_kt2 = 0.5 * eikonal.Numerics.MinKT2;

    for (std::size_t f1 = 0; f1 < 2; ++f1) {
      for (std::size_t f2 = 0; f2 < 2; ++f2) {
        const auto expected     = ExclusiveAmpReference(eikonal, kt2, f1, f2);
        const auto actual       = eikonal.S3ExclusiveAmp(kt2, f1, f2);
        const auto expected_t0  = ExclusiveAmpReference(eikonal, 0.0, f1, f2);
        const auto actual_t0    = eikonal.S3ExclusiveAmp(0.0, f1, f2);
        const auto expected_low = ExclusiveAmpReference(eikonal, low_kt2, f1, f2);
        const auto actual_low   = eikonal.S3ExclusiveAmp(low_kt2, f1, f2);

        REQUIRE(std::real(actual) == Approx(std::real(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
        REQUIRE(std::imag(actual) == Approx(std::imag(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
        REQUIRE(std::real(actual_t0) ==
                Approx(std::real(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
        REQUIRE(std::imag(actual_t0) ==
                Approx(std::imag(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
        REQUIRE(std::real(actual_low) ==
                Approx(std::real(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
        REQUIRE(std::imag(actual_low) ==
                Approx(std::imag(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
      }
    }

    const auto sd1 = eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation1);
    const auto sd2 = eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation2);
    const auto dd  = eikonal.ExclusiveFinalStates(gra::GWFinalClass::DoubleDissociation);
    REQUIRE(sd1.size() == 1);
    REQUIRE(sd1[0].f1 == 1);
    REQUIRE(sd1[0].f2 == 0);
    REQUIRE(sd2.size() == 1);
    REQUIRE(sd2[0].f1 == 0);
    REQUIRE(sd2[0].f2 == 1);
    REQUIRE(dd.size() == 1);
    REQUIRE(dd[0].f1 == 1);
    REQUIRE(dd[0].f2 == 1);
    RequirePhysicalPairScreeningProjection(eikonal, kt2);
  }

  SECTION("3-channel exclusive amplitudes drive inclusive cross sections") {
    auto         eikonal = BuildTestEikonal({0.17, -0.11, 0.23}, "n3");
    const double kt2     = eikonal.Numerics.MinKT2 + 0.61 * (eikonal.Numerics.MaxKT2 - eikonal.Numerics.MinKT2);
    const double low_kt2 = 0.5 * eikonal.Numerics.MinKT2;

    for (std::size_t f1 = 0; f1 < 3; ++f1) {
      for (std::size_t f2 = 0; f2 < 3; ++f2) {
        const auto expected     = ExclusiveAmpReference(eikonal, kt2, f1, f2);
        const auto actual       = eikonal.S3ExclusiveAmp(kt2, f1, f2);
        const auto expected_t0  = ExclusiveAmpReference(eikonal, 0.0, f1, f2);
        const auto actual_t0    = eikonal.S3ExclusiveAmp(0.0, f1, f2);
        const auto expected_low = ExclusiveAmpReference(eikonal, low_kt2, f1, f2);
        const auto actual_low   = eikonal.S3ExclusiveAmp(low_kt2, f1, f2);
        REQUIRE(std::real(actual) == Approx(std::real(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
        REQUIRE(std::imag(actual) == Approx(std::imag(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
        REQUIRE(std::real(actual_t0) ==
                Approx(std::real(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
        REQUIRE(std::imag(actual_t0) ==
                Approx(std::imag(expected_t0)).margin(ExclusiveProjectionMargin(eikonal, expected_t0)));
        REQUIRE(std::real(actual_low) ==
                Approx(std::real(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
        REQUIRE(std::imag(actual_low) ==
                Approx(std::imag(expected_low)).margin(ExclusiveProjectionMargin(eikonal, expected_low)));
      }
    }

    const auto sd1 = eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation1);
    const auto sd2 = eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation2);
    const auto dd  = eikonal.ExclusiveFinalStates(gra::GWFinalClass::DoubleDissociation);
    REQUIRE(sd1.size() == 2);
    REQUIRE(sd2.size() == 2);
    REQUIRE(dd.size() == 4);

    double expected_dd_intensity = 0.0;
    for (const auto &state : dd) {
      expected_dd_intensity += gra::math::abs2(eikonal.S3ExclusiveAmp(kt2, state.f1, state.f2));
    }
    REQUIRE(eikonal.S3DiffractiveAmpSquared(kt2, gra::GWFinalClass::DoubleDissociation) ==
            Approx(expected_dd_intensity).epsilon(1e-12));

    const auto &exclusive          = eikonal.GetExclusiveDiffXS();
    const auto &inclusive          = eikonal.GetInclusiveDiffXS();
    const auto  expected_inclusive = InclusiveFromExclusive(exclusive);
    for (std::size_t i = 0; i < 2; ++i) {
      for (std::size_t j = 0; j < 2; ++j) {
        REQUIRE(inclusive[i][j] == Approx(expected_inclusive[i][j]).epsilon(1e-10));
      }
    }

    double sigma_tot = 0.0;
    double sigma_el  = 0.0;
    double sigma_in  = 0.0;
    eikonal.GetTotXS(sigma_tot, sigma_el, sigma_in);
    REQUIRE(std::isfinite(sigma_tot));
    REQUIRE(std::isfinite(sigma_el));
    REQUIRE(std::isfinite(sigma_in));
    REQUIRE(sigma_tot != Approx(0.0).margin(1e-15));
    RequirePhysicalPairScreeningProjection(eikonal, kt2);
  }

  SECTION("4-channel amplitudes and classes are generic") {
    auto         eikonal = BuildTestEikonal({0.11, -0.07, 0.19, 0.23, -0.13, 0.29}, "n4");
    const double kt2     = eikonal.Numerics.MinKT2 + 0.43 * (eikonal.Numerics.MaxKT2 - eikonal.Numerics.MinKT2);

    REQUIRE(eikonal.GetChannelCount() == 4);
    for (std::size_t f1 = 0; f1 < 4; ++f1) {
      for (std::size_t f2 = 0; f2 < 4; ++f2) {
        const auto expected = ExclusiveAmpReference(eikonal, kt2, f1, f2);
        const auto actual   = eikonal.S3ExclusiveAmp(kt2, f1, f2);
        REQUIRE(std::real(actual) == Approx(std::real(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
        REQUIRE(std::imag(actual) == Approx(std::imag(expected)).margin(ExclusiveProjectionMargin(eikonal, expected)));
      }
    }

    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation1).size() == 3);
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation2).size() == 3);
    REQUIRE(eikonal.ExclusiveFinalStates(gra::GWFinalClass::DoubleDissociation).size() == 9);
  }

  SECTION("inclusive classes are invariant under excited-state rotations") {
    auto         first   = BuildTestEikonal({0.19, -0.07, 0.21}, "basis_first");
    auto         rotated = BuildTestEikonal({0.19, -0.07, -0.37}, "basis_rotated");
    const double kt2     = first.Numerics.MinKT2 + 0.47 * (first.Numerics.MaxKT2 - first.Numerics.MinKT2);

    for (const auto final_class : {gra::GWFinalClass::Elastic, gra::GWFinalClass::SingleDissociation1,
                                   gra::GWFinalClass::SingleDissociation2, gra::GWFinalClass::DoubleDissociation}) {
      REQUIRE(first.S3DiffractiveAmpSquared(kt2, final_class) ==
              Approx(rotated.S3DiffractiveAmpSquared(kt2, final_class)).epsilon(1e-10));
    }
  }
}

// Check scalar cut spectra against analytic exponential and q-eikonal maps
TEST_CASE("MEikonal scalar cut spectrum follows the selected unitarization", "[gra::MEikonal][CutSpectrum][physics]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");

  // Construct one scalar model with an optional second Pomeron
  const auto make_model = [](const bool auxiliary, const std::string &unitarization, const double q,
                             const std::string &suffix) {
    const std::string model_file =
        WriteSoftModelFile({}, "cut_spectrum_" + suffix, 0.25, 1, 0.0, 0.0, unitarization, q);
    auto              j                           = nlohmann::json::parse(gra::aux::GetInputData(model_file));
    auto             &soft                        = j.at("PARAM_SOFT");
    const std::string active                      = soft.at("active_model");
    auto             &model                       = soft.at("MODEL").at(active);
    soft.at("EXCHANGE_DEF")["P_aux"]              = soft.at("EXCHANGE_DEF").at("P");
    model.at("EXCHANGE")["P_aux"]                 = model.at("EXCHANGE").at("P");
    auto &p_aux                                   = model.at("EXCHANGE").at("P_aux");
    p_aux.at("on")                                = auxiliary;
    p_aux.at("g")[0][0]                           = 0.5 * model.at("EXCHANGE").at("P").at("g")[0][0].get<double>();
    model.at("EIKONAL").at("screening_exchanges") = {"O"};

    const std::filesystem::path numerics_file = std::filesystem::path(model_file).parent_path() / "NUMERICS.json";
    return gra::SoftModel::LoadFromJson(model_file, j.dump(), numerics_file.string(),
                                        gra::aux::GetInputData(numerics_file.string()));
  };

  // Build only the impact-parameter tables needed by the cut spectrum
  const auto build = [](const gra::SoftModelPtr &model) {
    gra::MEikonalNumerics numerics;
    numerics.ReadParameters(model->NumericsSourceFile(), model->NumericsSourceJson());
    numerics.NumberBT  = 8;
    numerics.NumberKT2 = 2;
    return gra::MEikonalMatrix::Build(model, 25.0, ProtonInitialState(), numerics, false);
  };

  const auto primary_model = make_model(false, "exp", 1.0, "primary");
  const auto exp_model     = make_model(true, "exp", 1.0, "exp");
  const auto q_model       = make_model(true, "q_exp", 0.63, "q_exp");
  const auto primary       = build(primary_model);
  const auto exponential   = build(exp_model);
  const auto q_eikonal     = build(q_model);

  // Count enabled Pomerons independently of the screening selector
  const auto pomeron_count = [](const gra::SoftModelPtr &model) {
    std::size_t count = 0;
    for (std::size_t exchange = 0; exchange < model->ExchangeCount(); ++exchange) {
      const auto &item = model->Exchange(gra::SoftExchangeId(exchange));
      if (item.Enabled() && item.Role() == gra::SoftExchangeRole::Pomeron) { ++count; }
    }
    return count;
  };
  REQUIRE(pomeron_count(primary_model) == 1);
  REQUIRE(pomeron_count(exp_model) == 2);
  REQUIRE(exp_model->Eikonal().ScreeningExchanges().size() == 1);
  REQUIRE(exp_model->Exchange(exp_model->Eikonal().ScreeningExchanges()[0]).Name() == "O");

  const double bt          = primary->ImpactParameterNodes()[1];
  const auto   primary_chi = primary->CutOpacity(bt);
  const auto   exp_chi     = exponential->CutOpacity(bt);
  const auto   q_chi       = q_eikonal->CutOpacity(bt);
  REQUIRE(primary_chi.size_row() == 1);
  REQUIRE(exp_chi.size_row() == 1);
  REQUIRE(std::abs(primary_chi(0, 0)) > 1.0e-12);
  // Independent transform builds agree within radial quadrature accuracy
  REQUIRE(std::abs(exp_chi(0, 0) - 1.25 * primary_chi(0, 0)) < 2.0e-7 * std::max(1.0, std::abs(exp_chi(0, 0))));
  REQUIRE(std::abs(q_chi(0, 0) - exp_chi(0, 0)) < 2.0e-7 * std::max(1.0, std::abs(exp_chi(0, 0))));

  const auto exp_spectrum = exponential->CutSpectrum(bt);
  REQUIRE(exp_spectrum.eigenvalues.size() == 1);
  REQUIRE(exp_spectrum.incoming_weights == std::vector<double>{1.0});
  const auto   exp_survival = std::exp(gra::math::zi * exp_chi(0, 0));
  const double exp_lambda   = 2.0 * std::imag(exp_chi(0, 0));
  REQUIRE(exp_spectrum.eigenvalues[0] == Approx(exp_lambda).epsilon(1.0e-11).margin(1.0e-13));
  REQUIRE(std::exp(-exp_spectrum.eigenvalues[0]) == Approx(std::norm(exp_survival)).epsilon(1.0e-11).margin(1.0e-13));

  const auto q_spectrum = q_eikonal->CutSpectrum(bt);
  REQUIRE(q_spectrum.eigenvalues.size() == 1);
  REQUIRE(q_spectrum.incoming_weights == std::vector<double>{1.0});
  const double q          = q_model->Eikonal().Q();
  const auto   q_base     = 1.0 + (1.0 - q) * gra::math::zi * q_chi(0, 0);
  const auto   q_survival = std::pow(q_base, 1.0 / (1.0 - q));
  const double q_lambda   = -std::log(std::norm(q_survival));
  REQUIRE(q_spectrum.eigenvalues[0] == Approx(q_lambda).epsilon(1.0e-11).margin(1.0e-13));
  REQUIRE(std::exp(-q_spectrum.eigenvalues[0]) == Approx(std::norm(q_survival)).epsilon(1.0e-11).margin(1.0e-13));
}

// Check finite and vanishing survival probabilities in the strongly absorptive limit
TEST_CASE("MEikonal cut spectrum reaches complete absorption", "[gra::MEikonal][CutSpectrum][physics]") {
  for (const double coupling : {48.0, 200.0, 400.0}) {
    CAPTURE(coupling);
    const auto tune = EikonalRegressionTune("black_disk_" + std::to_string(coupling), 1, 0.0, false, 64, coupling);
    gra::MEikonal eikonal(tune);
    REQUIRE_NOTHROW(eikonal.S3Constructor(25.0, ProtonInitialState()));
    const auto &runtime = eikonal.GetMatrixRuntime();
    const double b = runtime.ImpactParameterNodes().front();
    const double expected = 2.0 * std::imag(runtime.CutOpacity(b)(0, 0));
    const auto spectrum = runtime.CutSpectrum(b);
    REQUIRE(spectrum.eigenvalues.size() == 1);
    REQUIRE(spectrum.incoming_weights[0] == Approx(1.0).margin(1.0e-13));
    REQUIRE(expected > -std::log(1.0e-15));
    if (expected < 1400.0) {
      REQUIRE(spectrum.eigenvalues[0] == Approx(expected).epsilon(1.0e-10));
    } else {
      REQUIRE(std::isinf(spectrum.eigenvalues[0]));
      REQUIRE(spectrum.eigenvalues[0] > 0.0);
    }
    // Peripheral collisions still populate the retained finite multiplicities
    std::mt19937 random(31415);
    for (unsigned int sample = 0; sample < 32; ++sample) {
      unsigned int cuts = 0;
      double bt = 0.0;
      REQUIRE_NOTHROW(eikonal.S3GetRandomCutsBt(cuts, bt, random));
      REQUIRE(cuts > 0);
      REQUIRE(cuts < 25);
      REQUIRE(std::isfinite(bt));
    }
  }
}

// Check local and integrated unitarity identities for the active N=2 tune
TEST_CASE("Active two-channel MEikonal obeys unitarity identities",
          "[gra::MEikonal][unitarity][GoodWalker][CutSpectrum]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const std::string           model_file           = WriteSoftModelFile({0.31}, "strict_unitarity");
  const std::filesystem::path numerics_file        = std::filesystem::path(model_file).parent_path() / "NUMERICS.json";
  auto                        numerics             = nlohmann::json::parse(gra::aux::GetInputData(numerics_file));
  numerics["NUMERICS_EIKONAL"]["strict_unitarity"] = true;
  std::ofstream numerics_out(numerics_file);
  REQUIRE(numerics_out.good());
  numerics_out << numerics.dump(2);
  numerics_out.close();
  const auto model = gra::MModelTune::Load(model_file);
  REQUIRE(model->Soft()->GoodWalker().ChannelCount() == 2);

  gra::MEikonal    eikonal(model);
  constexpr double sqrts = 13000.0;
  REQUIRE_NOTHROW(eikonal.S3Constructor(sqrts * sqrts, ProtonInitialState(), false, 16, 4));
  REQUIRE(eikonal.Numerics.strict_unitarity);
  const double tolerance = eikonal.Numerics.unitarity_tolerance;
  const auto  &runtime   = eikonal.GetMatrixRuntime();
  REQUIRE(runtime.MaxSingularValue() <= 1.0 + tolerance);

  for (const double b : runtime.ImpactParameterNodes()) {
    gra::MEikonalMatrix::Matrix survival(4);
    for (std::size_t f1 = 0; f1 < 2; ++f1) {
      for (std::size_t f2 = 0; f2 < 2; ++f2) {
        const auto smatrix = runtime.PhysicalImpactSMatrix(b, f1, f2);
        survival += smatrix.Dagger() * smatrix;
      }
    }
    const auto residual = survival - survival.Dagger();
    CAPTURE(b, tolerance, residual.FrobNorm());
    REQUIRE(residual.FrobNorm() <= 1.0e-12 * std::max(1.0, survival.FrobNorm()));
    const auto hermitian   = (survival + survival.Dagger()) * std::complex<double>(0.5, 0.0);
    const auto eigenvalues = hermitian.SelfAdjointEigenvalues(1.0e-13);
    for (const double eigenvalue : eigenvalues) {
      CAPTURE(eigenvalue);
      REQUIRE(eigenvalue >= -tolerance);
      REQUIRE(eigenvalue <= 1.0 + tolerance);
    }
  }

  const double cut_b   = runtime.ImpactParameterNodes()[1];
  const auto   cut_chi = runtime.CutOpacity(cut_b);
  REQUIRE(cut_chi.size_row() == 4);
  REQUIRE(model->Soft()->Eikonal().Unitarization() == gra::SoftUnitarization::QExponential);
  const double q            = model->Soft()->Eikonal().Q();
  const auto   cut_base     = gra::MEikonalMatrix::Matrix::IdentityMatrix(4) + gra::math::zi * cut_chi * (1.0 - q);
  const auto   cut_smatrix  = cut_base.PrincipalPower(1.0 / (1.0 - q));
  const auto   survival_raw = cut_smatrix.Dagger() * cut_smatrix;
  const auto   survival     = (survival_raw + survival_raw.Dagger()) * std::complex<double>(0.5, 0.0);
  gra::MEikonalMatrix::Matrix eigenvectors;
  const auto                  survival_eigenvalues = survival.SelfAdjointEigenvalues(1.0e-13, &eigenvectors);
  const auto                  cut_spectrum         = runtime.CutSpectrum(cut_b);
  REQUIRE(cut_spectrum.eigenvalues.size() == survival_eigenvalues.size());
  REQUIRE(cut_spectrum.incoming_weights.size() == survival_eigenvalues.size());

  gra::MEikonalMatrix::Matrix reconstructed(4);
  double                      weight_sum        = 0.0;
  double                      spectral_survival = 0.0;
  for (const auto &mode : indices(survival_eigenvalues)) {
    CAPTURE(mode, survival_eigenvalues[mode]);
    REQUIRE(std::isfinite(cut_spectrum.eigenvalues[mode]));
    REQUIRE(cut_spectrum.eigenvalues[mode] >= 0.0);
    REQUIRE(std::isfinite(cut_spectrum.incoming_weights[mode]));
    REQUIRE(cut_spectrum.incoming_weights[mode] >= 0.0);
    const double clipped_survival = std::min(1.0, survival_eigenvalues[mode]);
    REQUIRE(std::exp(-cut_spectrum.eigenvalues[mode]) == Approx(clipped_survival).epsilon(1.0e-11).margin(1.0e-13));
    weight_sum += cut_spectrum.incoming_weights[mode];
    spectral_survival += cut_spectrum.incoming_weights[mode] * std::exp(-cut_spectrum.eigenvalues[mode]);

    const auto eigenvector = eigenvectors.Column(mode);
    auto       conjugate   = eigenvector;
    for (auto &value : conjugate) { value = std::conj(value); }
    reconstructed.AddOuterProduct(eigenvector, conjugate, survival_eigenvalues[mode]);
  }
  REQUIRE(weight_sum == Approx(1.0).margin(1.0e-13));
  REQUIRE((reconstructed - survival).FrobNorm() < 1.0e-11 * std::max(1.0, survival.FrobNorm()));

  const auto incoming_states = runtime.IncomingCutStates();
  REQUIRE_FALSE(incoming_states.empty());
  double direct_survival = 0.0;
  for (const auto &incoming : incoming_states) {
    const auto value = survival.MatrixElement(incoming, incoming);
    REQUIRE(std::abs(std::imag(value)) < 1.0e-12);
    direct_survival += std::real(value);
  }
  direct_survival /= static_cast<double>(incoming_states.size());
  REQUIRE(spectral_survival == Approx(direct_survival).epsilon(1.0e-11).margin(1.0e-13));

  double sigma_tot  = 0.0;
  double sigma_el   = 0.0;
  double sigma_inel = 0.0;
  eikonal.GetTotXS(sigma_tot, sigma_el, sigma_inel);
  REQUIRE(sigma_tot > 0.0);
  REQUIRE(sigma_el >= 0.0);
  REQUIRE(sigma_inel >= 0.0);
  REQUIRE(sigma_tot == Approx(sigma_el + sigma_inel).epsilon(1.0e-13).margin(1.0e-15));

  const auto &exclusive = eikonal.GetExclusiveDiffXS();
  const auto &inclusive = eikonal.GetInclusiveDiffXS();
  REQUIRE(exclusive.size_row() == 2);
  REQUIRE(exclusive.size_col() == 2);
  REQUIRE(inclusive.size_row() == 2);
  REQUIRE(inclusive.size_col() == 2);
  const std::array<gra::GWFinalClass, 4> classes = {gra::GWFinalClass::Elastic, gra::GWFinalClass::SingleDissociation1,
                                                    gra::GWFinalClass::SingleDissociation2,
                                                    gra::GWFinalClass::DoubleDissociation};
  const std::array<std::array<std::size_t, 2>, 4> inclusive_index = {
      std::array<std::size_t, 2>{0, 0}, {1, 0}, {0, 1}, {1, 1}};
  for (const auto &index : indices(classes)) {
    double exclusive_sum = 0.0;
    for (const auto &state : eikonal.ExclusiveFinalStates(classes[index])) {
      exclusive_sum += exclusive[state.f1][state.f2];
    }
    const auto row = inclusive_index[index][0];
    const auto col = inclusive_index[index][1];
    CAPTURE(index, row, col, exclusive_sum);
    REQUIRE(inclusive[row][col] == Approx(exclusive_sum).epsilon(1.0e-13).margin(1.0e-15));
  }
}

// Check the forward amplitude against the integrated total cross section
TEST_CASE("MEikonal optical theorem uses the initialized beam flux",
          "[gra::MEikonal][elastic][optical-theorem][physics]") {
  ModelParamRestoreGuard restore;
  auto                   initial_state = ProtonInitialState();
  initial_state[0].mass                = 1.31;
  initial_state[1].mass                = 1.47;
  constexpr double mandelstam_s        = 25.0;
  const auto eikonal = BuildTestEikonalWithInitialState({0.31}, "optical_theorem_beam_flux", initial_state, 0.25, 1, 16,
                                                        4, mandelstam_s);

  double sigma_tot  = 0.0;
  double sigma_el   = 0.0;
  double sigma_inel = 0.0;
  eikonal.GetTotXS(sigma_tot, sigma_el, sigma_inel);
  const double beta    = gra::kinematics::beta12(mandelstam_s, initial_state[0].mass, initial_state[1].mass);
  const double optical = std::imag(eikonal.S3ExclusiveAmp(0.0, 0, 0)) * gra::PDG::GeV2barn / (mandelstam_s * beta);

  REQUIRE(beta > 0.0);
  REQUIRE(sigma_tot == Approx(optical).epsilon(2.0e-12).margin(2.0e-14));
  REQUIRE(sigma_tot == Approx(sigma_el + sigma_inel).epsilon(1.0e-13).margin(1.0e-15));

  // The inverse Hankel profile scales with 1/sqrt(lambda) for fixed Born input
  gra::MEikonal light(eikonal.ModelTuneHandle());
  const auto    light_state = ProtonInitialState();
  REQUIRE_NOTHROW(light.S3Constructor(mandelstam_s, light_state, false, 16, 4));
  const double light_beta = gra::kinematics::beta12(mandelstam_s, light_state[0].mass, light_state[1].mass);
  REQUIRE(eikonal.GetMatrixRuntime().RuntimeFingerprint() != light.GetMatrixRuntime().RuntimeFingerprint());
  const auto  &impact    = eikonal.GetMatrixRuntime().ImpactParameterNodes();
  const double b         = impact[impact.size() / 3];
  const auto   heavy_chi = eikonal.S3Density(b, 0, 0);
  const auto   light_chi = light.S3Density(b, 0, 0);
  REQUIRE(std::abs(beta * heavy_chi - light_beta * light_chi) <
          2.0e-12 * std::max({1.0, std::abs(beta * heavy_chi), std::abs(light_beta * light_chi)}));
}

TEST_CASE("UPC nucleon profile uses the full eikonal inelastic cross section",
          "[gra::MEikonal][gra::nuclear][physics]") {
  ModelParamRestoreGuard restore;
  const auto             reference =
      BuildTestEikonalWithInitialState({0.31}, "upc_full_inelastic", ProtonInitialState(), 0.0, 1, 64, 4);
  gra::MEikonal density(reference.ModelTuneHandle());
  density.S3Constructor(25.0, ProtonInitialState(), true, 64, 2);

  double sigma_tot  = 0.0;
  double sigma_el   = 0.0;
  double sigma_inel = 0.0;
  density.GetTotXS(sigma_tot, sigma_el, sigma_inel);
  const auto   profile          = gra::nuclear::BuildNNProfile(density, std::nullopt);
  const double area             = gra::math::LinearRadialIntegral(profile.b_node, profile.inelastic);
  const auto   override_profile = gra::nuclear::BuildNNProfile(density, 70.0);
  const double override_area    = gra::math::LinearRadialIntegral(override_profile.b_node, override_profile.inelastic);

  REQUIRE(sigma_tot == Approx(sigma_el + sigma_inel).epsilon(1.0e-13).margin(1.0e-15));
  REQUIRE(profile.sigma == Approx(1.0e3 * sigma_inel).epsilon(1.0e-13));
  REQUIRE(10.0 * 2.0 * gra::math::PI * area == Approx(profile.sigma).epsilon(2.0e-13));
  REQUIRE(override_profile.sigma == Approx(70.0));
  REQUIRE(10.0 * 2.0 * gra::math::PI * override_area == Approx(override_profile.sigma).epsilon(2.0e-13));
  REQUIRE(override_profile.inelastic == profile.inelastic);
  const double radius_scale = std::sqrt(override_profile.sigma / profile.sigma);
  for (const auto &i : gra::aux::indices(profile.b_node)) {
    REQUIRE(override_profile.b_node[i] == Approx(radius_scale * profile.b_node[i]).epsilon(2.0e-14));
  }
  const auto  &runtime          = density.GetMatrixRuntime();
  const auto  &impact           = runtime.ImpactParameterNodes();
  const double tolerance        = runtime.UnitarityTolerance();
  const double unitarity_margin = tolerance * (2.0 + tolerance);
  REQUIRE(runtime.MaxSingularValue() <= 1.0 + tolerance);
  for (const auto &i : gra::aux::indices(profile.inelastic)) {
    const double probability = profile.inelastic[i];
    REQUIRE(probability >= 0.0);
    REQUIRE(probability <= 1.0);

    // Evaluate 1 - Tr(S_el^dagger S_el)/4 independently from the complete
    // complex impact-parameter helicity matrix, S_el = I + i A_el
    const auto amplitude        = runtime.ImpactHelicityMatrix(impact[i], 0, 0);
    double     elastic_survival = 0.0;
    for (std::size_t row = 0; row < 4; ++row) {
      for (std::size_t col = 0; col < 4; ++col) {
        const auto value = (row == col ? std::complex<double>(1.0, 0.0) : std::complex<double>(0.0, 0.0)) +
                           gra::math::zi * amplitude[gra::spin::PairHelicityMatrixIndex(row, col)];
        elastic_survival += std::norm(value) / 4.0;
      }
    }
    const double raw_probability = 1.0 - elastic_survival;
    CAPTURE(i, profile.b_node[i], raw_probability, unitarity_margin);
    REQUIRE(raw_probability >= -unitarity_margin);
    REQUIRE(raw_probability <= 1.0);
    REQUIRE(probability == Approx(std::clamp(raw_probability, 0.0, 1.0)).margin(2.0e-12));
  }
}

TEST_CASE("MEikonal matrix amplitudes reject momentum above the table", "[gra::MEikonal][validation]") {
  auto         eikonal = BuildTestEikonal({}, "matrix_momentum_range");
  const auto  &runtime = eikonal.GetMatrixRuntime();
  const double outside = runtime.MomentumTransferNodes().back() * 1.01 + 1.0e-12;

  REQUIRE_THROWS_AS(runtime.HelicityAmplitudes(outside, 0, 0), std::out_of_range);
  REQUIRE_THROWS_AS(runtime.CrossingHelicityAmplitudes(outside, 0, 0), std::out_of_range);
  REQUIRE_THROWS_AS(runtime.ScreeningHelicityAmplitudes(outside, 0, 0), std::out_of_range);
  REQUIRE_THROWS_AS(eikonal.S3ExclusiveAmp(outside, 0, 0), std::out_of_range);
  REQUIRE_THROWS_AS(eikonal.S3ExclusiveScreeningAmp(outside, 0, 0), std::out_of_range);
  REQUIRE_THROWS_AS(eikonal.S3PhysicalScreeningAmp(outside), std::out_of_range);
}

TEST_CASE("MEikonal caches the physical proton screening amplitude", "[gra::MEikonal]") {
  auto        eikonal    = BuildTestEikonal({0.17, -0.11, 0.23}, "exclusive_loop_cache");
  const auto &loop_const = eikonal.GetLoopConst(25.0);

  REQUIRE(loop_const.physical_screening_amplitude.size() == loop_const.kt2.size());
  REQUIRE(loop_const.pair_screening_spin.size() == loop_const.kt2.size());
  REQUIRE(loop_const.physical_screening_weight.size_row() == loop_const.kt2.size());
  REQUIRE(loop_const.physical_screening_weight.size_col() == loop_const.node_weight.size_col());
  REQUIRE(loop_const.physical_screening_helicity_weight.size() ==
          loop_const.kt2.size() * loop_const.node_weight.size_col());
  REQUIRE(loop_const.physical_screening_is_scalar);
  for (const std::size_t kt_index : {std::size_t{0}, loop_const.kt2.size() / 2, loop_const.kt2.size() - 1}) {
    const auto direct = eikonal.S3PhysicalScreeningAmp(loop_const.kt2[kt_index]);
    const auto stored = loop_const.physical_screening_amplitude[kt_index];
    REQUIRE(std::real(stored) == Approx(std::real(direct)).epsilon(1e-12));
    REQUIRE(std::imag(stored) == Approx(std::imag(direct)).epsilon(1e-12));
    for (std::size_t phi_index = 0; phi_index < loop_const.node_weight.size_col(); ++phi_index) {
      const auto expected = loop_const.node_weight[kt_index][phi_index] * direct;
      const auto weight   = loop_const.physical_screening_weight[kt_index][phi_index];
      REQUIRE(std::real(weight) == Approx(std::real(expected)).epsilon(1e-12));
      REQUIRE(std::imag(weight) == Approx(std::imag(expected)).epsilon(1e-12));

      const double azimuth = std::atan2(loop_const.kt_y[kt_index][phi_index], loop_const.kt_x[kt_index][phi_index]);
      auto         expected_matrix =
          eikonal.GetMatrixRuntime().ScreeningHelicityMatrix(loop_const.kt2[kt_index], 0, 0, 0, 0, azimuth);
      for (auto &component : expected_matrix) { component *= loop_const.node_weight[kt_index][phi_index]; }
      const auto &stored_matrix =
          loop_const.physical_screening_helicity_weight[kt_index * loop_const.node_weight.size_col() + phi_index];
      for (std::size_t component = 0; component < stored_matrix.size(); ++component) {
        REQUIRE(std::real(stored_matrix[component]) == Approx(std::real(expected_matrix[component])).epsilon(1e-12));
        REQUIRE(std::imag(stored_matrix[component]) == Approx(std::imag(expected_matrix[component])).epsilon(1e-12));
      }
    }
  }

  const auto &single_dissociation =
      loop_const.good_walker_channels[static_cast<std::size_t>(gra::GWFinalClass::SingleDissociation1)];
  const auto &resolved_coefficients = eikonal.SoftModelHandle()->GoodWalker().ResolvedCoefficients();
  REQUIRE(single_dissociation.size() == 2);
  for (std::size_t channel = 0; channel < single_dissociation.size(); ++channel) {
    REQUIRE(single_dissociation[channel].spin_scalar);
    REQUIRE(single_dissociation[channel].screening_weight.empty());
    REQUIRE(single_dissociation[channel].scalar_screening_weight.size() ==
            loop_const.kt2.size() * loop_const.node_weight.size_col());
    REQUIRE(single_dissociation[channel].final_state.f1 == channel + 1);
    REQUIRE(single_dissociation[channel].final_state.f2 == 0);
    REQUIRE(single_dissociation[channel].born_coefficient == Approx(resolved_coefficients[channel]).epsilon(1e-12));
  }

  const std::size_t         kt_index    = loop_const.kt2.size() / 2;
  const std::size_t         phi_index   = 0;
  const auto                final_state = single_dissociation.front().final_state;
  gra::ProtonHelicityMatrix effective{};
  for (std::size_t source = 1; source < eikonal.GetChannelCount(); ++source) {
    const auto   term = eikonal.GetMatrixRuntime().ScreeningHelicityMatrix(loop_const.kt2[kt_index], final_state.f1,
                                                                           final_state.f2, source, 0);
    const double coefficient = resolved_coefficients[source - 1];
    gra::AddScaled(effective, term, coefficient);
  }
  auto expected_transition = gra::RotateProtonHelicityMatrix(
      effective, std::atan2(loop_const.kt_y[kt_index][phi_index], loop_const.kt_x[kt_index][phi_index]));
  for (auto &component : expected_transition) { component *= loop_const.node_weight[kt_index][phi_index]; }
  const std::size_t node = kt_index * loop_const.node_weight.size_col() + phi_index;
  REQUIRE(std::real(single_dissociation.front().scalar_screening_weight[node]) ==
          Approx(std::real(expected_transition[0])).epsilon(1e-12));
  REQUIRE(std::imag(single_dissociation.front().scalar_screening_weight[node]) ==
          Approx(std::imag(expected_transition[0])).epsilon(1e-12));
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const auto        expected = row == col ? expected_transition[0] : 0.0;
      const std::size_t index    = gra::spin::PairHelicityMatrixIndex(row, col);
      REQUIRE(std::real(expected_transition[index]) == Approx(std::real(expected)).epsilon(1e-12));
      REQUIRE(std::imag(expected_transition[index]) == Approx(std::imag(expected)).epsilon(1e-12));
    }
  }
}

TEST_CASE("MEikonal helicity amplitude squared keeps leg-specific flips", "[gra::MEikonal][helicity]") {
  gra::ElasticHelicityAmplitudes amplitude;
  amplitude.phi5       = 2.0;
  amplitude.phi5_first = 3.0;
  REQUIRE(gra::UnpolarizedHelicityAmpSquared(amplitude) == Approx(13.0));
}

TEST_CASE("MEikonal constructs the complete pair-helicity matrix", "[gra::MEikonal][helicity]") {
  gra::ElasticHelicityAmplitudes amplitude;
  amplitude.phi1                           = {1.0, 0.2};
  amplitude.phi2                           = {2.0, 0.3};
  amplitude.phi3                           = {3.0, 0.4};
  amplitude.phi4                           = {4.0, 0.5};
  amplitude.phi5                           = {5.0, 0.6};
  amplitude.phi5_first                     = {6.0, 0.7};
  const auto                      matrix   = gra::BuildProtonHelicityMatrix(amplitude);
  const gra::ProtonHelicityMatrix expected = {
      amplitude.phi1, -amplitude.phi5,       amplitude.phi5_first,  amplitude.phi2, amplitude.phi5, amplitude.phi3,
      amplitude.phi4, amplitude.phi5_first,  -amplitude.phi5_first, amplitude.phi4, amplitude.phi3, -amplitude.phi5,
      amplitude.phi2, -amplitude.phi5_first, amplitude.phi5,        amplitude.phi1};
  REQUIRE(matrix == expected);

  constexpr std::array<double, 4> reciprocity = {1.0, -1.0, -1.0, 1.0};
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      REQUIRE(matrix[gra::spin::PairHelicityMatrixIndex(row, col)] ==
              reciprocity[row] * reciprocity[col] * matrix[gra::spin::PairHelicityMatrixIndex(col, row)]);
    }
  }

  const auto rotated             = gra::BuildProtonHelicityMatrix(amplitude, 0.73);
  double     matrix_amp_squared  = 0.0;
  double     rotated_amp_squared = 0.0;
  for (std::size_t component = 0; component < matrix.size(); ++component) {
    matrix_amp_squared += gra::math::abs2(matrix[component]);
    rotated_amp_squared += gra::math::abs2(rotated[component]);
  }
  REQUIRE(rotated_amp_squared == Approx(matrix_amp_squared).epsilon(1e-12));

  gra::ProtonHelicityMatrix arbitrary{};
  for (std::size_t entry = 0; entry < arbitrary.size(); ++entry) {
    arbitrary[entry] = {0.17 * static_cast<double>(entry + 1), -0.09 * static_cast<double>(entry + 2)};
  }
  const double   azimuth           = -0.81;
  const auto     arbitrary_rotated = gra::RotateProtonHelicityMatrix(arbitrary, azimuth);
  constexpr auto transition        = gra::CanonicalProtonHelicityTransitions();
  for (std::size_t entry = 0; entry < arbitrary.size(); ++entry) {
    const auto expected_phase =
        std::exp(gra::math::zi * static_cast<double>(transition[entry].azimuth_harmonic) * azimuth);
    REQUIRE(std::abs(arbitrary_rotated[entry] - arbitrary[entry] * expected_phase) < 1.0e-13);
  }
}

TEST_CASE("MEikonal retains oriented Good Walker spin blocks", "[gra::MEikonal][helicity][GoodWalker]") {
  auto eikonal = BuildTestEikonalWithInitialState({0.31}, "oriented_spin_blocks", ProtonInitialState(), 0.25, 1, 4, 4,
                                                  25.0, 0.18, 0.65);
  const double kt2  = eikonal.Numerics.MinKT2 + 0.43 * (eikonal.Numerics.MaxKT2 - eikonal.Numerics.MinKT2);
  const auto   spin = eikonal.PairScreeningSpinBank(kt2);
  constexpr std::array<double, 4> reciprocity    = {1.0, -1.0, -1.0, 1.0};
  const std::size_t               channel_count  = eikonal.GetMatrixRuntime().ChannelCount();
  const std::size_t               pair_dimension = channel_count * channel_count;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      for (std::size_t final = 0; final < pair_dimension; ++final) {
        for (std::size_t initial = 0; initial < pair_dimension; ++initial) {
          const auto actual = spin[gra::spin::PairHelicityMatrixIndex(row, col)](final, initial);
          const auto expected =
              reciprocity[row] * reciprocity[col] * spin[gra::spin::PairHelicityMatrixIndex(col, row)](initial, final);
          CAPTURE(row, col, final, initial, actual, expected);
          REQUIRE(std::abs(actual - expected) < 2.0e-10 * std::max(1.0, std::abs(expected)));
          // Parity reverses both helicities with the collider spin phases
          const auto parity = reciprocity[row] * reciprocity[col] *
                              spin[gra::spin::PairHelicityMatrixIndex(3 - row, 3 - col)](final, initial);
          REQUIRE(std::abs(actual - parity) < 2.0e-10 * std::max(1.0, std::abs(actual)));
          // Beam exchange reverses the transfer and interchanges both Good Walker indices
          const std::size_t swapped_row = 2 * (row % 2) + row / 2;
          const std::size_t swapped_col = 2 * (col % 2) + col / 2;
          const std::size_t swapped_final = channel_count * (final % channel_count) + final / channel_count;
          const std::size_t swapped_initial = channel_count * (initial % channel_count) + initial / channel_count;
          const auto exchanged = reciprocity[row] * reciprocity[col] *
                                 spin[gra::spin::PairHelicityMatrixIndex(swapped_row, swapped_col)](swapped_final,
                                                                                                  swapped_initial);
          REQUIRE(std::abs(actual - exchanged) < 2.0e-10 * std::max(1.0, std::abs(actual)));
        }
      }
    }
  }

  const auto elastic               = eikonal.GetMatrixRuntime().HelicityMatrix(kt2, 0, 0);
  const auto elastic_named         = eikonal.GetMatrixRuntime().HelicityAmplitudes(kt2, 0, 0);
  const auto elastic_reconstructed = gra::BuildProtonHelicityMatrix(elastic_named);
  for (std::size_t entry = 0; entry < elastic.size(); ++entry) {
    REQUIRE(std::abs(elastic[entry] - elastic_reconstructed[entry]) <
            2.0e-10 * std::max(1.0, std::abs(elastic[entry])));
  }

  const auto forward         = eikonal.GetMatrixRuntime().HelicityMatrix(kt2, 1, 0, 0, 0);
  const auto reverse         = eikonal.GetMatrixRuntime().HelicityMatrix(kt2, 0, 0, 1, 0);
  double     transition_norm = 0.0;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const auto actual   = forward[gra::spin::PairHelicityMatrixIndex(row, col)];
      const auto expected = reciprocity[row] * reciprocity[col] * reverse[gra::spin::PairHelicityMatrixIndex(col, row)];
      transition_norm += std::norm(actual);
      REQUIRE(std::abs(actual - expected) < 2.0e-10 * std::max(1.0, std::abs(expected)));
    }
  }
  REQUIRE(transition_norm > 1.0e-20);

  double expected_sd1 = 0.0;
  for (const auto &state : eikonal.ExclusiveFinalStates(gra::GWFinalClass::SingleDissociation1)) {
    expected_sd1 += 0.25 * gra::SquaredNorm(eikonal.GetMatrixRuntime().HelicityMatrix(kt2, state.f1, state.f2));
  }
  REQUIRE(eikonal.S3DiffractiveAmpSquared(kt2, gra::GWFinalClass::SingleDissociation1) ==
          Approx(expected_sd1).epsilon(1.0e-12));
}

TEST_CASE("MEikonal named pair-helicity bank follows the canonical transitions", "[gra::MEikonal][helicity]") {
  gra::ElasticHelicityAmplitudes amplitude;
  amplitude.phi1       = {1.0, 0.1};
  amplitude.phi2       = {2.0, 0.2};
  amplitude.phi3       = {3.0, 0.3};
  amplitude.phi4       = {4.0, 0.4};
  amplitude.phi5       = {5.0, 0.5};
  amplitude.phi5_first = {6.0, 0.6};

  gra::MEikonalMatrix::PairHelicityBank bank;
  bank.phi1       = {{amplitude.phi1}};
  bank.phi2       = {{amplitude.phi2}};
  bank.phi3       = {{amplitude.phi3}};
  bank.phi4       = {{amplitude.phi4}};
  bank.phi5       = {{amplitude.phi5}};
  bank.phi5_first = {{amplitude.phi5_first}};

  const double azimuth    = 0.37;
  const auto   expected   = gra::BuildProtonHelicityMatrix(amplitude, azimuth);
  const auto   transition = gra::CanonicalProtonHelicityTransitions();
  for (std::size_t index = 0; index < transition.size(); ++index) {
    const auto                &term = transition[index];
    const std::complex<double> actual =
        term.sign * bank.Component(term.component)(0, 0) *
        std::exp(gra::math::zi * (static_cast<double>(term.azimuth_harmonic) * azimuth));
    REQUIRE(std::real(actual) == Approx(std::real(expected[index])).epsilon(1e-12));
    REQUIRE(std::imag(actual) == Approx(std::imag(expected[index])).epsilon(1e-12));
  }
}

TEST_CASE("MEikonal screening loop quadrature coefficients", "[gra::MEikonal]") {
  gra::MEikonal eikonal;

  SECTION("1/3") {
    eikonal.Numerics.LOOP.radial_integrator = "1/3";
    REQUIRE(eikonal.GetLoopQuadratureCoefficient() == Approx(1.0 / 3.0));
  }

  SECTION("3/8") {
    eikonal.Numerics.LOOP.radial_integrator = "3/8";
    REQUIRE(eikonal.GetLoopQuadratureCoefficient() == Approx(3.0 / 8.0));
  }

  SECTION("Boole") {
    eikonal.Numerics.LOOP.radial_integrator = "Boole";
    REQUIRE(eikonal.GetLoopQuadratureCoefficient() == Approx(2.0 / 45.0));
  }

  SECTION("GL") {
    eikonal.Numerics.LOOP.radial_integrator = "GL";
    REQUIRE(eikonal.GetLoopQuadratureCoefficient() == Approx(1.0));
  }
}

TEST_CASE("MEikonal screening loop quadrature covers full 2D range", "[gra::MEikonal][MProcess]") {
  const double min_kt   = 0.2;
  const double max_kt   = 1.2;
  const double min_phi  = 0.0;
  const double max_phi  = 2.0 * gra::math::PI;
  const double expected = (max_kt - min_kt) * (max_phi - min_phi);

  SECTION("1/3") {
    const unsigned int    n_kt  = 2;
    const unsigned int    n_phi = 2;
    const MMatrix<double> f(n_kt + 1, n_phi + 1, 1.0);
    const MMatrix<double> w = gra::math::Simpson13Weight2D(n_kt, n_phi);
    const double value = gra::math::Simpson13Integral2D(f, w, (max_kt - min_kt) / n_kt, (max_phi - min_phi) / n_phi);

    REQUIRE(value == Approx(expected).epsilon(1e-14));
  }

  SECTION("3/8") {
    const unsigned int    n_kt  = 3;
    const unsigned int    n_phi = 3;
    const MMatrix<double> f(n_kt + 1, n_phi + 1, 1.0);
    const MMatrix<double> w = gra::math::Simpson38Weight2D(n_kt, n_phi);
    const double value = gra::math::Simpson38Integral2D(f, w, (max_kt - min_kt) / n_kt, (max_phi - min_phi) / n_phi);

    REQUIRE(value == Approx(expected).epsilon(1e-14));
  }

  SECTION("Boole") {
    const unsigned int    n_kt  = 4;
    const unsigned int    n_phi = 4;
    const MMatrix<double> f(n_kt + 1, n_phi + 1, 1.0);
    const MMatrix<double> w = gra::math::BooleWeight2D(n_kt, n_phi);
    const double value      = gra::math::BooleIntegral2D(f, w, (max_kt - min_kt) / n_kt, (max_phi - min_phi) / n_phi);

    REQUIRE(value == Approx(expected).epsilon(1e-14));
  }
}

TEST_CASE("MEikonal applies independent signed loop node offsets", "[gra::MEikonal][params]") {
  ModelParamRestoreGuard restore;

  auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = "GL";
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "Trap";
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopKT"]   = 5;
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopPHI"]  = 7;

  const std::filesystem::path dir = "tmp/graniitti_eikonal_gl_trap_counts";
  std::filesystem::create_directories(dir);

  const std::filesystem::path path = dir / "NUMERICS.json";
  std::ofstream               out(path);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  gra::MODELPARAM = dir.string();
  gra::MEikonalNumerics numerics;
  numerics.SetLoopDiscretization(2, -3);
  REQUIRE_NOTHROW(numerics.ReadParameters(path.string()));
  REQUIRE(numerics.LOOP.radial_integrator == "GL");
  REQUIRE(numerics.LOOP.azimuth_integrator == "Trap");
  REQUIRE(numerics.LOOP.radial_intervals == 7);
  REQUIRE(numerics.LOOP.azimuth_nodes == 4);
}

TEST_CASE("MEikonal validates decoupled loop integrator names and radial counts", "[gra::MEikonal][params]") {
  ModelParamRestoreGuard restore;

  auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  const std::filesystem::path dir = "tmp/graniitti_eikonal_decoupled_names";
  std::filesystem::create_directories(dir);
  const std::filesystem::path path = dir / "NUMERICS.json";

  auto write_card = [&]() {
    std::ofstream out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();
    gra::MODELPARAM = dir.string();
  };

  SECTION("unknown kT integrator") {
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = "BadKT";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "Trap";
    write_card();
    gra::MEikonalNumerics numerics;
    REQUIRE_THROWS(numerics.ReadParameters(path.string()));
  }

  SECTION("unknown phi integrator") {
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = "GL";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "BadPhi";
    write_card();
    gra::MEikonalNumerics numerics;
    REQUIRE_THROWS(numerics.ReadParameters(path.string()));
  }

  SECTION("incompatible radial count") {
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = "Boole";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "Trap";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopKT"]   = 10;
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopPHI"]  = 7;
    write_card();
    gra::MEikonalNumerics numerics;
    REQUIRE_THROWS(numerics.ReadParameters(path.string()));
  }

  SECTION("odd impact parameter interval count") {
    j["NUMERICS_EIKONAL"]["NumberBT"] = 5;
    write_card();
    gra::MEikonalNumerics numerics;
    REQUIRE_THROWS(numerics.ReadParameters(path.string()));
  }

  SECTION("arbitrary phi count for radial Newton-Cotes") {
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = "Boole";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "Trap";
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopKT"]   = 12;
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopPHI"]  = 7;
    write_card();
    gra::MEikonalNumerics numerics;
    REQUIRE_NOTHROW(numerics.ReadParameters(path.string()));
    REQUIRE(numerics.LOOP.azimuth_nodes == 7);
  }
}

TEST_CASE("MEikonal requires positive MinLoopKT for logarithmic loop grids", "[gra::MEikonal][params]") {
  ModelParamRestoreGuard restore;

  auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["log_kT"]        = true;
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["MinLoopKT"]     = 0.0;
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopKT"]  = 12;
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["NumberLoopPHI"] = 16;

  const std::filesystem::path dir = "tmp/graniitti_eikonal_log_loop_min";
  std::filesystem::create_directories(dir);

  const std::filesystem::path path = dir / "NUMERICS.json";

  for (const std::string kT_integrator : {"Boole", "GL"}) {
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["kT_integrator"]  = kT_integrator;
    j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["phi_integrator"] = "Trap";

    std::ofstream out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    gra::MODELPARAM = dir.string();
    gra::MEikonalNumerics numerics;
    REQUIRE_THROWS(numerics.ReadParameters(path.string()));
  }
}

TEST_CASE("MEikonal GL x Trap loop constants use node dimensions and periodic phi", "[gra::MEikonal]") {
  auto eikonal                             = BuildTestEikonal({}, "loop_gl_trap");
  eikonal.Numerics.LOOP.radial_integrator  = "GL";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
  eikonal.Numerics.LOOP.r_min              = 0.2;
  eikonal.Numerics.LOOP.r_max              = 1.2;
  eikonal.Numerics.LOOP.radial_intervals   = 5;
  eikonal.Numerics.LOOP.azimuth_nodes      = 7;
  eikonal.InitLoopWeightMatrix();

  const auto &loop_const = eikonal.GetLoopConst(25.0);
  REQUIRE(loop_const.kt2.size() == 5);
  REQUIRE(loop_const.W.size_row() == 5);
  REQUIRE(loop_const.W.size_col() == 7);
  REQUIRE(loop_const.measure_weight.size_row() == 5);
  REQUIRE(loop_const.measure_weight.size_col() == 7);
  REQUIRE(loop_const.kt_x.size_row() == 5);
  REQUIRE(loop_const.kt_x.size_col() == 7);

  REQUIRE(loop_const.StepPhi == Approx(2.0 * gra::math::PI / 7.0));
  REQUIRE(loop_const.kt_x[0][0] == Approx(loop_const.kt[0] * std::cos(gra::math::PI / 7.0)));
  REQUIRE(loop_const.kt_y[0][0] == Approx(loop_const.kt[0] * std::sin(gra::math::PI / 7.0)));
  REQUIRE(loop_const.kt_y[0][0] != Approx(loop_const.kt_y[0][6]).epsilon(1e-12));
  REQUIRE(loop_const.kt_y[0][6] != Approx(0.0).margin(1e-12));

  double weight_sum = 0.0;
  for (std::size_t i = 0; i < loop_const.measure_weight.size_row(); ++i) {
    for (std::size_t j = 0; j < loop_const.measure_weight.size_col(); ++j) {
      weight_sum += loop_const.measure_weight[i][j];
    }
  }
  const double expected = gra::math::PI * (pow2(1.2) - pow2(0.2));
  REQUIRE(weight_sum == Approx(expected).epsilon(1e-13));
}

TEST_CASE(
    "MEikonal Boole x Trap loop constants use radial intervals and "
    "periodic phi nodes",
    "[gra::MEikonal]") {
  auto eikonal                             = BuildTestEikonal({}, "loop_boole_trap");
  eikonal.Numerics.LOOP.radial_integrator  = "Boole";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Linear;
  eikonal.Numerics.LOOP.r_min              = 0.2;
  eikonal.Numerics.LOOP.r_max              = 1.2;
  eikonal.Numerics.LOOP.radial_intervals   = 4;
  eikonal.Numerics.LOOP.azimuth_nodes      = 7;
  eikonal.InitLoopWeightMatrix();

  const auto &loop_const = eikonal.GetLoopConst(25.0);
  REQUIRE(loop_const.kt2.size() == 5);
  REQUIRE(loop_const.W.size_row() == 5);
  REQUIRE(loop_const.W.size_col() == 7);
  REQUIRE(loop_const.kt_x.size_col() == 7);
  REQUIRE(loop_const.StepPhi == Approx(2.0 * gra::math::PI / 7.0));
  REQUIRE(loop_const.kt_x[1][0] == Approx(loop_const.kt[1] * std::cos(gra::math::PI / 7.0)));
  REQUIRE(loop_const.kt_y[1][0] == Approx(loop_const.kt[1] * std::sin(gra::math::PI / 7.0)));

  double weight_sum = 0.0;
  for (std::size_t i = 0; i < loop_const.measure_weight.size_row(); ++i) {
    for (std::size_t j = 0; j < loop_const.measure_weight.size_col(); ++j) {
      weight_sum += loop_const.measure_weight[i][j];
    }
  }
  const double expected = gra::math::PI * (pow2(1.2) - pow2(0.2));
  REQUIRE(weight_sum == Approx(expected).epsilon(1e-13));
}

TEST_CASE("MEikonal GL x Trap logarithmic kt integrates physical radial measure", "[gra::MEikonal]") {
  auto eikonal                             = BuildTestEikonal({}, "loop_gl_trap_log");
  eikonal.Numerics.LOOP.radial_integrator  = "GL";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.radial_map         = gra::math::RadialMap::Log;
  eikonal.Numerics.LOOP.r_min              = 0.2;
  eikonal.Numerics.LOOP.r_max              = 1.2;
  eikonal.Numerics.LOOP.radial_intervals   = 12;
  eikonal.Numerics.LOOP.azimuth_nodes      = 7;
  eikonal.InitLoopWeightMatrix();

  const auto &loop_const = eikonal.GetLoopConst(25.0);
  double      weight_sum = 0.0;
  for (std::size_t i = 0; i < loop_const.measure_weight.size_row(); ++i) {
    for (std::size_t j = 0; j < loop_const.measure_weight.size_col(); ++j) {
      weight_sum += loop_const.measure_weight[i][j];
    }
  }
  const double expected = gra::math::PI * (pow2(1.2) - pow2(0.2));
  REQUIRE(weight_sum == Approx(expected).epsilon(1e-12));
}

TEST_CASE(
    "MEikonal loop initializer requires positive MinLoopKT for "
    "logarithmic grids",
    "[gra::MEikonal]") {
  auto eikonal                           = BuildTestEikonal({}, "loop_log_min_guard");
  eikonal.Numerics.LOOP.radial_map       = gra::math::RadialMap::Log;
  eikonal.Numerics.LOOP.r_min            = 0.0;
  eikonal.Numerics.LOOP.r_max            = 1.2;
  eikonal.Numerics.LOOP.radial_intervals = 12;
  eikonal.Numerics.LOOP.azimuth_nodes    = 12;

  eikonal.Numerics.LOOP.radial_integrator  = "Boole";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  REQUIRE_THROWS(eikonal.InitLoopWeightMatrix());

  eikonal.Numerics.LOOP.radial_integrator  = "GL";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  REQUIRE_THROWS(eikonal.InitLoopWeightMatrix());
}

TEST_CASE("MEikonal loop initializer requires positive loop counts", "[gra::MEikonal]") {
  auto eikonal                             = BuildTestEikonal({}, "loop_positive_counts");
  eikonal.Numerics.LOOP.radial_integrator  = "GL";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.r_min              = 0.0;
  eikonal.Numerics.LOOP.r_max              = 1.2;
  eikonal.Numerics.LOOP.radial_intervals   = 0;
  eikonal.Numerics.LOOP.azimuth_nodes      = 7;
  REQUIRE_THROWS(eikonal.InitLoopWeightMatrix());

  eikonal.Numerics.LOOP.radial_intervals = 5;
  eikonal.Numerics.LOOP.azimuth_nodes    = 0;
  REQUIRE_THROWS(eikonal.InitLoopWeightMatrix());
}

TEST_CASE("MEikonal screening loop constants store grid values", "[gra::MEikonal]") {
  auto eikonal                             = BuildTestEikonal({}, "loop_const");
  eikonal.Numerics.LOOP.radial_integrator  = "1/3";
  eikonal.Numerics.LOOP.azimuth_integrator = "Trap";
  eikonal.Numerics.LOOP.r_min              = 0.0;
  eikonal.Numerics.LOOP.r_max              = 0.2;
  eikonal.Numerics.LOOP.radial_intervals   = 2;
  eikonal.Numerics.LOOP.azimuth_nodes      = 2;
  eikonal.InitLoopWeightMatrix();

  const auto                &loop_const = eikonal.GetLoopConst(25.0);
  const auto                 initial    = ProtonInitialState();
  const double               beta       = gra::kinematics::beta12(25.0, initial[0].mass, initial[1].mass);
  const std::complex<double> norm       = gra::math::zi / (8.0 * gra::math::PIPI * 25.0 * beta);

  REQUIRE(loop_const.kt2[0] == Approx(0.0));
  REQUIRE(loop_const.kt2[1] == Approx(0.01));
  REQUIRE(loop_const.kt2[2] == Approx(0.04));

  REQUIRE(loop_const.kt_x[1][0] == Approx(0.0).margin(1e-15));
  REQUIRE(loop_const.kt_y[1][0] == Approx(0.1));
  REQUIRE(loop_const.kt_x[1][1] == Approx(0.0).margin(1e-15));
  REQUIRE(loop_const.kt_y[1][1] == Approx(-0.1));

  REQUIRE(loop_const.node_weight[1][1].real() == Approx(0.0).margin(1e-15));
  REQUIRE(loop_const.W[1][1] == Approx((0.1 / 3.0) * 4.0 * gra::math::PI));
  REQUIRE(loop_const.node_weight[1][1].imag() == Approx((norm * loop_const.measure_weight[1][1]).imag()));
}

TEST_CASE("MEikonal rejects loop kt range outside interpolation domain", "[gra::MEikonal][params]") {
  ModelParamRestoreGuard restore;

  auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  j["NUMERICS_EIKONAL"]["MaxKT2"]                     = 1.0;
  j["NUMERICS_EIKONAL"]["LOOP_INTEGRAL"]["MaxLoopKT"] = 2.0;

  const std::filesystem::path dir = "tmp/graniitti_eikonal_invalid_loop";
  std::filesystem::create_directories(dir);

  const std::filesystem::path path = dir / "NUMERICS.json";
  std::ofstream               out(path);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  gra::MODELPARAM = dir.string();

  gra::MEikonalNumerics numerics;
  REQUIRE_THROWS(numerics.ReadParameters(path.string()));
}

// A constant hard source isolates the physical eikonal convolution and final-state projections
TEST_CASE("Constant hard sources obey the Good Walker screening sum", "[gra::MEikonal][screening][GoodWalker]") {
  for (const auto &angles : std::vector<std::vector<double>>{{}, {0.29}, {0.19, -0.07, 0.21},
                                                            {0.11, -0.07, 0.19, 0.23, -0.13, 0.29}}) {
    auto eikonal = BuildTestEikonal(angles, "screen_constant_" + std::to_string(angles.size()));
    const auto &loop = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
    const std::vector<std::complex<double>> born = {{1.25, -0.5}};
    gra::ScreeningMetadata metadata{gra::ScreeningSpinBasis::ProtonIdentity,
                                    gra::ProtonScreeningMode::ForwardExcitation, 1.0,
                                    gra::ScreeningAmplitudeType::Physical, 1, false};
    for (const auto final : {gra::GWFinalClass::Elastic, gra::GWFinalClass::SingleDissociation1,
                             gra::GWFinalClass::SingleDissociation2, gra::GWFinalClass::DoubleDissociation}) {
      CAPTURE(angles.size(), final);
      if (loop.good_walker_channels[static_cast<std::size_t>(final)].empty()) {
        REQUIRE(final != gra::GWFinalClass::Elastic);
        REQUIRE(eikonal.GetChannelCount() == 1);
        continue;
      }
      gra::eikonal::MProtonScreen screen({born}, metadata, loop, final);
      for (const auto &i : indices(loop.kt2)) {
        for (std::size_t j = 0; j < loop.node_weight.size_col(); ++j) { screen.Add(i, j, {born}); }
      }
      const auto result = screen.Result();
      REQUIRE(result.size() == 1);
      CHECK(gra::SquaredNorm(result.front()) == Approx(ExpectedToyGoodWalkerAmpSquared(eikonal, born[0], final)).epsilon(1e-10).margin(1e-24));
    }
  }
}

// Contract every proton helicity transition through the real physical kernel
TEST_CASE("Physical screening contracts every proton helicity transition", "[gra::MEikonal][screening][helicity]") {
  auto eikonal = BuildTestEikonal({}, "screen_full_helicity", 4, 4, 0.65);
  const auto &loop = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
  REQUIRE_FALSE(loop.physical_screening_is_scalar);
  std::vector<std::complex<double>> born(16);
  for (const auto &i : indices(born)) { born[i] = {0.11 * (i + 1), -0.04 * (i + 2)}; }
  gra::ScreeningMetadata metadata{gra::ScreeningSpinBasis::ProtonHelicity, gra::ProtonScreeningMode::Elastic,
                                  0.25, gra::ScreeningAmplitudeType::Physical, 16, false};
  metadata.PrepareSpinTransitions();
  gra::eikonal::MProtonScreen screen({born}, metadata, loop, gra::GWFinalClass::Elastic);
  for (const auto &i : indices(loop.kt2)) {
    for (std::size_t j = 0; j < loop.node_weight.size_col(); ++j) { screen.Add(i, j, {born}); }
  }
  const auto result = screen.Result();
  REQUIRE(result.size() == 1);
  CHECK(0.25 * gra::SquaredNorm(result.front()) == Approx(ExpectedToyHelicityScreenedAmpSquared(eikonal, born)).epsilon(1e-10));
}

// Exercise the actual Coulomb and nuclear process route with on-shell proton momenta
TEST_CASE("Elastic CNI keeps physical helicities and classifies invalid transfers", "[gra::MEikonal][MProcess][CNI]") {
  ProcessProbe process;

    process.eikonal = BuildTestEikonal({}, "screen_elastic_helicity", 4, 4, 0.65);
    process.ProcPtr.Initialize("X", "EL");
    gra::ScreeningMetadata metadata{gra::ScreeningSpinBasis::ProtonHelicity,
                                    gra::ProtonScreeningMode::Elastic,
                                    0.25,
                                    gra::ScreeningAmplitudeType::Physical,
                                    16,
                                    false};
    metadata.amplitude_type = gra::ScreeningAmplitudeType::ElasticCNI;
    process.state.lts.hamp.Configure(metadata);

    const auto   initialstate   = ProtonInitialState();
    const double momentum       = gra::kinematics::DecayMomentum(5.0, initialstate[0].mass, initialstate[1].mass);
    const auto   momenta        = ElasticProtonMomenta(momentum, 0.1, 0.0);
    process.state.lts.beam1     = initialstate[0];
    process.state.lts.beam2     = initialstate[1];
    process.state.lts.pbeam1    = momenta[0];
    process.state.lts.pbeam2    = momenta[1];
    process.state.lts.pfinal[1] = momenta[2];
    process.state.lts.pfinal[2] = momenta[3];
    process.state.lts.s         = 25.0;
    process.state.lts.t         = -0.1;

    double total_before     = 0.0;
    double elastic_before   = 0.0;
    double inelastic_before = 0.0;
    process.eikonal.GetTotXS(total_before, elastic_before, inelastic_before);
    process.eikonal.Numerics.CNI.NumberKT2       = 32;
    process.eikonal.Numerics.CNI.FBIntegralMaxKT = 2.0;
    process.eikonal.Numerics.CNI.FBIntegralN     = 64;
    process.eikonal.Numerics.CNI.tail_max_qb     = 64.0;
    process.eikonal.Numerics.CNI.TailIntegralN   = 128;
    process.eikonal.Numerics.CNI.interp_rel_tol  = 1.0e3;
    process.eikonal.InitializeElasticCNI(0.01, 0.2);

    const auto expected = process.eikonal.PhysicalElasticHelicityMatrix(momenta[0], momenta[1], momenta[2], momenta[3]);
    const double expected_amp2 = 0.25 * gra::SquaredNorm(expected);
    REQUIRE(process.ScreenedAmplitudeSquared(false) == Approx(expected_amp2).epsilon(1e-12));
    REQUIRE(process.ScreenedAmplitudeSquared(true) == Approx(expected_amp2).epsilon(1e-12));
    REQUIRE(process.state.lts.hamp.size() == 16);
    for (std::size_t entry = 0; entry < expected.size(); ++entry) {
      REQUIRE(process.state.lts.hamp[entry] == expected[entry]);
    }

    const double total_azimuth   = 0.47;
    const auto   rotated_momenta = ElasticProtonMomenta(momentum, 0.1, total_azimuth);
    const auto   rotated_total   = process.eikonal.PhysicalElasticHelicityMatrix(rotated_momenta[0], rotated_momenta[1],
                                                                                 rotated_momenta[2], rotated_momenta[3]);
    const auto   expected_rotated = gra::RotateProtonHelicityMatrix(expected, total_azimuth);
    for (std::size_t entry = 0; entry < rotated_total.size(); ++entry) {
      REQUIRE(std::abs(rotated_total[entry] - expected_rotated[entry]) <
              2.0e-11 * std::max(1.0, std::abs(expected_rotated[entry])));
    }

    const auto        limit_momenta = ElasticProtonMomenta(momentum, 0.01, 0.0);
    const auto        full_limit    = process.eikonal.PhysicalElasticHelicityMatrix(limit_momenta[0], limit_momenta[1],
                                                                                    limit_momenta[2], limit_momenta[3]);
    const auto        strong_limit  = process.eikonal.GetMatrixRuntime().HelicityMatrix(0.01, 0, 0);
    const gra::MDirac dirac("DIRAC");
    const auto        born_limit =
        gra::qed::ElasticSpinHalfPhotonExchange(dirac, initialstate[0], initialstate[1], limit_momenta[0],
                                                limit_momenta[1], limit_momenta[2], limit_momenta[3]);
    auto higher_order_limit = full_limit;
    gra::AddScaled(higher_order_limit, strong_limit, -1.0);
    gra::AddScaled(higher_order_limit, born_limit, -1.0);
    const double born_norm = std::sqrt(gra::SquaredNorm(born_limit));
    REQUIRE(born_norm > 0.0);
    REQUIRE(std::sqrt(gra::SquaredNorm(higher_order_limit)) / born_norm < 0.5);

    double total_after     = 0.0;
    double elastic_after   = 0.0;
    double inelastic_after = 0.0;
    process.eikonal.GetTotXS(total_after, elastic_after, inelastic_after);
    REQUIRE(gra::math::IsExactEqual(total_after, total_before));
    REQUIRE(gra::math::IsExactEqual(elastic_after, elastic_before));
    REQUIRE(gra::math::IsExactEqual(inelastic_after, inelastic_before));

    const auto overflow_momenta = ElasticProtonMomenta(momentum, 0.21, 0.0);
    process.state.lts.pfinal[1] = overflow_momenta[2];
    process.state.lts.pfinal[2] = overflow_momenta[3];
    process.state.lts.t         = -0.21;
    gra::MEventWeightState overflow_aux;
    double                 overflow_amp2 = -1.0;
    REQUIRE_NOTHROW(overflow_amp2 = process.EvaluateAmplitudeBoundary(true, overflow_aux));
    REQUIRE(overflow_amp2 == Approx(0.0));
    REQUIRE(overflow_aux.amplitude_failure);
    REQUIRE(overflow_aux.technical_failure);

    process.state.lts.hamp.Configure(metadata);
    process.state.lts.pfinal[1] = momenta[2];
    process.state.lts.pfinal[2] = momenta[3];
    process.state.lts.pfinal[1].SetPx(std::numeric_limits<double>::quiet_NaN());
    process.state.lts.t = std::numeric_limits<double>::quiet_NaN();
    gra::MEventWeightState nonfinite_aux;
    double                 nonfinite_amp2 = -1.0;
    REQUIRE_NOTHROW(nonfinite_amp2 = process.EvaluateAmplitudeBoundary(true, nonfinite_aux));
    REQUIRE(nonfinite_amp2 == Approx(0.0));
    REQUIRE(nonfinite_aux.amplitude_failure);
    REQUIRE(nonfinite_aux.technical_failure);

}

// Check generic validation of exceptional amplitude outputs and final state
// Validate actual amplitude banks directly instead of replacing a matrix element
TEST_CASE("The common amplitude validator rejects nonfinite scalar and helicity values", "[gra::MProcess][amplitude][exceptions]") {
  ProcessProbe process;
  process.state.lts.hamp = {1.0};
  REQUIRE_NOTHROW(process.ValidateAmplitude(1.0));
  for (const double invalid : {-1.0, std::numeric_limits<double>::infinity(),
                               std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(process.ValidateAmplitude(invalid), gra::AmplitudeFailure);
  }
  for (const double invalid : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    process.state.lts.hamp = {std::complex<double>(1.0, invalid)};
    CHECK_THROWS_AS(process.ValidateAmplitude(1.0), gra::AmplitudeFailure);
  }
  process.state.lts.hamp.clear();
  REQUIRE_NOTHROW(process.ValidateAmplitude(0.0));
}

// Reject incompatible source banks in the same kernel used for physical and color amplitudes
TEST_CASE("Proton screening rejects changed helicity and color bank dimensions", "[gra::MProcess][screening][failure]") {
  gra::MEikonal::LoopConst loop;
  loop.physical_screening_is_scalar = true;
  loop.physical_screening_weight = gra::MMatrix<std::complex<double>>(1, 1, std::complex<double>(0.2, -0.1));
  gra::ScreeningMetadata metadata{gra::ScreeningSpinBasis::ProtonIdentity, gra::ProtonScreeningMode::Elastic};
  const std::vector<std::complex<double>> source = {{1.0, 0.2}, {-0.3, 0.7}};
  gra::eikonal::MProtonScreen screen({source, source}, metadata, loop, gra::GWFinalClass::Elastic);
  REQUIRE_NOTHROW(screen.Add(0, 0, {source, source}));
  CHECK_THROWS_AS(screen.Add(0, 0, {source}), gra::AmplitudeFailure);
  const std::vector<std::complex<double>> changed = {{1.0, 0.2}};
  CHECK_THROWS_AS(screen.Add(0, 0, {changed, source}), gra::AmplitudeFailure);
  CHECK_THROWS_AS(screen.Add(0, 0, {source, changed}), gra::AmplitudeFailure);
}

// Validate the actual event containers before screening and color selection
TEST_CASE("Amplitude validation requires matching finite Good Walker and color sources", "[gra::MProcess][amplitude][failure]") {
  ProcessProbe process;
  auto &lts = process.state.lts;
  lts.hamp = {{1.0, 0.2}};
  REQUIRE_NOTHROW(process.ValidateAmplitude(1.04));
  lts.hamp.metadata.amplitude_type = gra::ScreeningAmplitudeType::GoodWalker;
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  lts.proton_good_walker.emplace();
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  lts.proton_good_walker->components.emplace_back();
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  auto &source = lts.proton_good_walker->components.front().source;
  source = gra::MMatrix<std::complex<double>>(1, 1, std::complex<double>(1.0, 0.2));
  REQUIRE_NOTHROW(process.ValidateAmplitude(1.04));
  source[0][0] = {std::numeric_limits<double>::quiet_NaN(), 0.0};
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  lts.hamp.metadata.amplitude_type = gra::ScreeningAmplitudeType::Physical;
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  lts.proton_good_walker.reset();
  lts.hard_color_flows.emplace_back();
  auto &flow = lts.hard_color_flows.front();
  flow.amplitudes = {{1.0, 0.2}};
  REQUIRE_NOTHROW(process.ValidateAmplitude(1.04));
  for (const double invalid : {-1.0, std::numeric_limits<double>::quiet_NaN(), std::numeric_limits<double>::infinity()}) {
    flow.screened_weight = invalid;
    CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
  }
  flow.screened_weight.reset();
  flow.amplitudes[0] = {0.0, std::numeric_limits<double>::infinity()};
  CHECK_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
}


TEST_CASE("Soft Pomeron residues use the physical multichannel proton",
          "[gra::MForm][gra::MEikonal][gra::MRegge][physics]") {
  const std::string data          = gra::aux::GetInputData(modelfile);
  auto              j             = nlohmann::json::parse(data);
  j["PARAM_SOFT"]["active_model"] = "triple";
  auto &soft                      = j["PARAM_SOFT"]["MODEL"]["triple"];
  soft["GW"]["theta"]             = {0.31, 0.47, 0.83};
  soft["EXCHANGE"]["P"]["g"]      = {{4.0, 0.0, 0.0}, {0.0, 7.0, 0.0}, {0.0, 0.0, 11.0}};
  soft["3P"]["g"]                 = 0.2;

  const std::string temp_modelfile = WriteSoftModelFile({}, "physical_pomeron_coupling_test");
  std::ofstream     out(temp_modelfile);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  const auto  model_tune        = gra::MModelTune::Load(temp_modelfile);
  const auto &soft_model        = *model_tune->Soft();
  const auto  pomeron_id        = soft_model.ExchangeId("P");
  const auto &pomeron           = soft_model.Exchange(pomeron_id);
  const auto &proton            = soft_model.GoodWalker().ProtonVector();
  double      expected_coupling = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) { expected_coupling += proton[i] * pomeron.Coupling(i, j) * proton[j]; }
  }
  REQUIRE(soft_model.PhysicalCoupling(pomeron_id) == Approx(expected_coupling).epsilon(1e-12));
  REQUIRE(soft_model.EffectiveTriplePomeronCoupling(pomeron_id) == Approx(0.2 * expected_coupling).epsilon(1e-12));

  const double t                = -0.19;
  const auto   residue          = soft_model.ResidueMatrix(pomeron_id, t);
  double       expected_residue = 0.0;
  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) { expected_residue += proton[i] * residue(i, j) * proton[j]; }
  }
  const double expected_form_factor = expected_residue / expected_coupling;
  REQUIRE(soft_model.PhysicalResidue(pomeron_id, t) == Approx(expected_residue).epsilon(1e-12));
  REQUIRE(soft_model.PhysicalResidue(pomeron_id, t) / expected_coupling == Approx(expected_form_factor).epsilon(1e-12));
  REQUIRE(soft_model.PhysicalResidue(pomeron_id, 0.0) / expected_coupling == Approx(1.0).epsilon(1e-12));

  const double beta_pi = (2.0 / 3.0) * expected_coupling;
  const double pion_coefficient =
      beta_pi * beta_pi * gra::PDG::mpi * gra::PDG::mpi / (32.0 * gra::math::pow3(gra::math::PI));
  const double pion_mass2       = gra::PDG::mpi * gra::PDG::mpi;
  const double loop_ratio       = 4.0 * pion_mass2 / std::abs(t);
  const double loop_root        = std::sqrt(1.0 + loop_ratio);
  const double pion_form_factor = 1.0 / (1.0 - t / soft_model.PionLoopScale2());
  const double pion_loop =
      (4.0 / loop_ratio) * pion_form_factor * pion_form_factor *
      (2.0 * loop_ratio - std::pow(1.0 + loop_ratio, 1.5) * std::log((loop_root + 1.0) / (loop_root - 1.0)) +
       std::log(1.0 / pion_mass2));
  const double expected_alpha = pomeron.Alpha0() + pomeron.AlphaPrime() * t - pion_coefficient * pion_loop;
  REQUIRE(soft_model.Alpha(pomeron_id, t) == Approx(expected_alpha).epsilon(1e-12));

  gra::LORENTZSCALAR lts = MakeToyCoherentPhotonLTS();
  lts.t                  = t;
  lts.t1                 = t;
  lts.ss[1][1]           = 4.0;
  lts.ss[2][2]           = 9.0;
  gra::MRegge regge(lts, model_tune, gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Generic, "test"));
  const auto  param = gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *regge.ModelTuneHandle());

  const gra::SoftExchangeId mapped_pomeron  = param.exchanges.at(param.pomeron_trajectory).soft_exchange;
  const double              expected_vertex = regge.SoftModelHandle()->PhysicalResidue(mapped_pomeron, lts.t1);
  REQUIRE(expected_vertex == Approx(expected_residue).epsilon(1e-12));

  const double alpha = soft_model.Alpha(mapped_pomeron, t);
  const double g3p   = soft_model.EffectiveTriplePomeronCoupling(mapped_pomeron);
  const double eta2 =
      std::norm(gra::regge::EtaFactor(alpha, pomeron.Alpha0(), pomeron.Signature(), soft_model.TriplePomeronEtaMode()));
  const double expected_sd2 = eta2 * expected_coupling * expected_residue * expected_residue * g3p *
                              std::pow(lts.s / 4.0, 2.0 * alpha) * std::pow(4.0, pomeron.Alpha0());
  const double expected_dd2 = eta2 * expected_coupling * expected_coupling * g3p * g3p *
                              std::pow(lts.s / 36.0, 2.0 * alpha) * std::pow(4.0, pomeron.Alpha0()) *
                              std::pow(9.0, pomeron.Alpha0());
  lts.excite1 = true;
  lts.excite2 = false;
  REQUIRE(std::norm(regge.ME2(lts, gra::MReggeInclusive::SD)) == Approx(expected_sd2).epsilon(1e-12));
  lts.excite2 = true;
  REQUIRE(std::norm(regge.ME2(lts, gra::MReggeInclusive::DD)) == Approx(expected_dd2).epsilon(1e-12));
}

// Check GGCF model isolation and one-sided forward diffraction
TEST_CASE("GGCF resolves a separate eikonal and forward variance", "[gra::nuclear][ggcf][physics]") {
  const bool helicity = GENERATE(false, true);
  auto tune = EikonalRegressionTune("ggcf", 2, helicity ? 0.1 : 0.0, false, 256);
  auto general = tune->General();
  general["PARAM_SOFT"]["MODEL"]["ggcf_test"] = general["PARAM_SOFT"]["MODEL"]["single"];
  general["PARAM_SOFT"]["MODEL"]["ggcf_test"]["GW"]["theta"] = {0.43};
  std::ofstream(tune->GeneralFile()) << general.dump();
  tune = gra::MModelTune::Load(tune->GeneralFile());
  gra::nuclear::GGCFParam param;
  param.eikonal = "ggcf_test";
  const auto selected = gra::nuclear::GGCFTune(tune, param);
  CHECK(tune->Soft()->ActiveModel() == "single");
  CHECK(selected->Soft()->ActiveModel() == param.eikonal);
  CHECK(selected->General("PARAM_NUCLEAR") == tune->General("PARAM_NUCLEAR"));
  const auto proton = ProtonInitialState().front();
  gra::MEikonal eikonal(selected);
  eikonal.S3Constructor(1600.0, {proton, proton}, false, 256, 32);
  const auto moments = gra::nuclear::ForwardGW(eikonal.GetMatrixRuntime());
  const auto exchanged = gra::nuclear::ForwardGW(eikonal.GetMatrixRuntime(), 1);
  const auto elastic = eikonal.S3ExclusiveAmp(0.0, 0, 0);
  const auto diss = eikonal.S3ExclusiveAmp(0.0, 1, 0);
  const auto reversed = eikonal.S3ExclusiveAmp(0.0, 0, 1);
  REQUIRE(moments.omega > 0.0);
  CHECK(std::abs(diss - reversed) < 1.0e-10 * std::abs(elastic));
  const double ratio = eikonal.S3DiffractiveAmpSquared(0.0, gra::GWFinalClass::SingleDissociation1) /
                       eikonal.S3DiffractiveAmpSquared(0.0, gra::GWFinalClass::Elastic);
  CHECK(moments.omega == Approx(ratio).epsilon(0.002));
  CHECK(moments.omega == Approx(exchanged.omega).epsilon(1.0e-10));
  double total = 0.0, el = 0.0, inel = 0.0;
  eikonal.GetTotXS(total, el, inel);
  CHECK(moments.sigma == Approx(1000.0 * total).epsilon(0.002));
  param.eikonal = "missing_ggcf_model";
  CHECK_THROWS_AS(gra::nuclear::GGCFTune(tune, param), std::invalid_argument);
  param.eikonal = "single";
  const auto single = EikonalRegressionTune("ggcf_single");
  CHECK_THROWS_AS(gra::nuclear::GGCFTune(single, param), std::invalid_argument);
  gra::MEikonal fixed(tune);
  fixed.S3Constructor(1600.0, {proton, proton}, true, 256, 2);
  CHECK(gra::nuclear::ForwardGW(fixed.GetMatrixRuntime()).omega == Approx(0.0).margin(1.0e-14));
}
