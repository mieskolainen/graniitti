// Dispersive vector propagation and coherent pion decay tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Tensor/MTensorPhoto.h"
#include "Graniitti/Tensor/MTensorVector.h"
#include "support/models_test_support.hh"

// Compare the analytic loop with a subtracted principal-value dispersion integral
TEST_CASE("Vector P-wave loop obeys its dispersion relation", "[tensor][vector][dispersion]") {
  const double mass = 0.23;
  const auto [node, weight] = gra::math::GaussLegendreRule(512, 0.0, 1.0);
  for (const double x : {-4.0, -0.2501, -0.2499, -1.0e-10, 0.0, 1.0e-10, 0.01, 0.2499, 0.2501,
                         0.75, 0.99, 1.01, 1.7, 30.0}) {
    CAPTURE(x);
    const double s = 4.0 * mass * mass * x;
    const double at_pole = x > 1.0 ? std::pow(1.0 - 1.0 / x, 1.5) : 0.0;
    long double integral = x > 1.0 ? -at_pole * std::log(x - 1.0) / x : 0.0;
    for (const auto &i : indices(node)) {
      integral += weight[i] * (std::pow(1.0 - node[i], 1.5) - at_pole) / (1.0 - x * node[i]);
    }
    const auto loop = gra::tensor::VectorLoop(s, mass);
    const double expected = s * x * static_cast<double>(integral) / (192.0 * pow2(gra::math::PI));
    CHECK(loop.real() == Approx(expected).epsilon(2.0e-8).margin(1.0e-24));
    CHECK(loop.imag() == Approx(s * at_pole / (192.0 * gra::math::PI)).margin(1.0e-18));
  }
}

// Check the real subtraction conditions and the complex matrix optical identity
TEST_CASE("Coupled vector propagators preserve subtractions and absorptive phases", "[tensor][vector][phase]") {
  const auto model = gra::MModelTune::Load(modelfile);
  const auto parameters = gra::ReadTensorPomeronParam(*model, LoadedPDGTable());
  const auto &mix = parameters->rho_omega;
  const auto zero = mix.Propagator(0.0);
  for (const auto &i : indices(mix.mass)) {
    CHECK(zero[i][i].real() == Approx(-1.0 / pow2(mix.mass[i])).epsilon(1.0e-13));
    CHECK(std::abs(zero[i][1-i]) < 1.0e-15);
    const auto inverse = mix.Inverse(pow2(mix.mass[i]));
    CHECK(std::abs(inverse[i][i].real()) < 1.0e-14);
  }
  CHECK(mix.Inverse(pow2(mix.mass[0]))[0][1].real() == Approx(pow2(mix.mass[0]) * mix.b).epsilon(1.0e-13));
  for (const double scale : {0.7, 1.0, 1.02, 1.6, 2.7}) {
    const double s = scale * pow2(mix.mass[0]);
    const auto inverse = mix.Inverse(s);
    const auto delta = mix.Propagator(s);
    const auto identity = inverse * delta;
    const auto absorptive = (inverse - inverse.Conj()) * std::complex<double>(0.0, -0.5);
    const auto cut = delta.Dagger() * absorptive * delta;
    for (const auto &i : indices(mix.mass)) {
      for (const auto &j : indices(mix.mass)) {
        RequireComplexNear(identity[i][j], i == j ? 1.0 : 0.0, 2.0e-12);
        RequireComplexNear(delta[i][j], delta[j][i], 2.0e-12);
        RequireComplexNear(cut[i][j], -delta[i][j].imag(), 2.0e-11);
      }
    }
  }
  for (const double s : {std::numeric_limits<double>::infinity(), std::numeric_limits<double>::quiet_NaN()}) {
    CHECK_THROWS_AS(mix.Propagator(s), gra::AmplitudeFailure);
  }
}

// Match the absorptive pion cut to the physical spin-averaged decay amplitude
TEST_CASE("Dispersive rho pion cut has the covariant decay normalization", "[tensor][vector][normalization]") {
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  const auto model = gra::MModelTune::Load(modelfile);
  gra::MTensorPomeron tensor(lts, model, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Generic));
  const auto parameters = gra::GetTensorParam(*lts.model_cache, lts.PDG, lts.process.RESONANCES);
  const auto &mix = parameters->rho_omega;
  const double mass = mix.mass[0];
  const double momentum = std::sqrt(pow2(mass) / 4.0 - pow2(mix.pion));
  const gra::M4Vec plus(momentum, 0.0, 0.0, mass / 2.0);
  const gra::M4Vec minus(-momentum, 0.0, 0.0, mass / 2.0);
  const auto vertex = tensor.iG_vpsps(plus, minus, mass, mix.g[0], gra::regge::FFParam{});
  const auto polarizations = tensor.MassiveSpin1States(plus + minus, "conj", true);
  double amplitude2 = 0.0;
  for (const auto &h : indices(polarizations)) {
    std::complex<double> amplitude = 0.0;
    for (const auto &mu : tensor.LI) { amplitude += polarizations[h](mu) * vertex(mu); }
    amplitude2 += std::norm(amplitude);
  }
  const double width = gra::kinematics::PDW2body(pow2(mass), pow2(mix.pion), pow2(mix.pion), amplitude2 / 3.0, 1.0);
  CHECK(mass * width == Approx(pow2(mix.g[0]) * gra::tensor::VectorLoop(pow2(mass), mix.pion).imag()).epsilon(2.0e-12));
  for (const int pdg : {113, 223}) {
    const auto row = mix.Index(pdg);
    const auto delta = mix.Propagator(pow2(mass));
    const auto current = tensor.VectorDecay(plus, minus, mix.mass[row], parameters->FindVector(pdg).width,
                                            pdg, mix.g[row], gra::regge::FFParam{});
    const auto spectral = delta[row][0] * mix.g[0] + delta[row][1] * mix.g[1];
    const auto conjugate = tensor.VectorDecay(minus, plus, mix.mass[row], parameters->FindVector(pdg).width,
                                              pdg, mix.g[row], gra::regge::FFParam{});
    for (const auto &mu : tensor.LI) {
      RequireComplexNear(current(mu), -0.5 * spectral * (plus - minus)[mu], 2.0e-12);
      RequireComplexNear(conjugate(mu), -current(mu), 2.0e-12);
    }
  }
}

// Check coherent vector sources in the real photon current and helicity basis
TEST_CASE("Rho and omega PHOTO amplitudes add coherently and rotate covariantly", "[tensor][vector][photo][phase]") {
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 0.91, -0.38);
  for (const int pdg : {113, 223}) {
    gra::PARAM_RES res;
    res.p = lts.PDG.FindByPDG(pdg);
    SetToyTensorChannel(res, {0.43, -0.18});
    res.hel_decay.g_decay_TP = {1.0};
    lts.process.RESONANCES.emplace(std::to_string(pdg), res);
  }
  const bool noflip = GENERATE(false, true);
  const auto tune = WriteModifiedPhotoVMTune("rho_omega_covariance_" + std::to_string(noflip), [&](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip;
  });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  const auto evaluate = [&](gra::LORENTZSCALAR event) {
    event.amplitude.BeginCentral();
    REQUIRE(tensor.MEPhoto(event) > 0.0);
    return std::vector<std::complex<double>>(event.hamp.begin(), event.hamp.end());
  };
  const auto total = evaluate(lts);
  auto rho = lts;
  rho.process.RESONANCES.erase("223");
  auto omega = lts;
  omega.process.RESONANCES.erase("113");
  auto continuum = lts;
  continuum.process.RESONANCES.clear();
  const auto r = evaluate(rho), w = evaluate(omega), ds = evaluate(continuum);
  for (const auto &i : indices(total)) { RequireComplexNear(total[i], r[i] + w[i] - ds[i], 3.0e-10); }
  std::swap(lts.decaytree[0], lts.decaytree[1]);
  RequireVectorNear(evaluate(lts), total, 3.0e-10);
  const double norm = gra::SquaredNorm(total);
  for (const auto &event : {RotateToyEventAroundZ(lts, 0.73), ReflectToyEventInXZ(lts), BeamExchangeMirrorWithDecay(lts)}) {
    CHECK(gra::SquaredNorm(evaluate(event)) == Approx(norm).epsilon(3.0e-9));
  }
  // A longitudinal boost mixes outgoing proton helicities through Wigner rotations
  if (!noflip) {
    CHECK(gra::SquaredNorm(evaluate(BoostToyEventAlongZ(lts, -0.31))) == Approx(norm).epsilon(3.0e-9));
  }
}

// Reject invalid coupled poles and ambiguous pion couplings during initialization
TEST_CASE("Coupled vector inputs are checked before sampling", "[tensor][vector][params]") {
  for (const int variation : {0, 1, 2, 3}) {
    const auto tune = WriteModifiedPhotoVMTune("vector_input_" + std::to_string(variation), [&](auto &j) {
      auto &tensor = j.at("PARAM_TENSORPOM");
      if (variation == 0) { tensor.at("RHO_OMEGA").at("g") = {0.7}; }
      if (variation == 1) { tensor.at("RHO_OMEGA").at("b") = "invalid"; }
      if (variation == 2) { tensor.at("VECTOR").at("Wmode").at(1) = "RHO_OMEGA"; }
      if (variation == 3) { tensor.at("VECTOR").at("dPDG").at(0) = 321; }
    });
    CHECK_THROWS(gra::ReadTensorPomeronParam(*gra::MModelTune::Load(tune.second), LoadedPDGTable()));
  }
}

// Bind the pion vertex to the spectral coupling while preserving generic BR normalization
TEST_CASE("Coupled pion vertices share the spectral parameters", "[tensor][vector][decay][normalization]") {
  const int pdg = GENERATE(113, 223);
  const auto evaluate = [&](bool varied) {
    const auto tune = WriteModifiedPhotoVMTune("vector_decay_" + std::to_string(pdg) + std::to_string(varied),
        [&](auto &j) {
          auto &coupling = j.at("PARAM_TENSORPOM").at("RHO_OMEGA").at("g").at(pdg == 113 ? 0 : 1);
          if (varied) { coupling = 1.3 * coupling.template get<double>(); }
        });
    ToyHelicityProcess process;
    process.SetTuneForTest(tune.first);
    process.SetProcessForTest("TP", "RES");
    process.state.lts.PDG = LoadedPDGTable();
    const auto &particles = process.state.lts.PDG;
    return process.ProcessHelicityStructure(particles.FindByPDG(pdg),
        {particles.FindByPDG(211), particles.FindByPDG(-211)}, false, true, "", false);
  };
  const auto nominal = evaluate(false);
  const auto varied = evaluate(true);
  REQUIRE(nominal.g_decay_TP.size() == 1);
  REQUIRE(varied.g_decay_TP.size() == 1);
  CHECK(varied.g_decay_TP[0] == Approx(1.3 * nominal.g_decay_TP[0]).epsilon(2.0e-13));
  RequireComplexNear(varied.g_decay, nominal.g_decay, 2.0e-13);
}
