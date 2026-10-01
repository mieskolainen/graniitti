// Covariant spin-three photoproduction and F-wave normalization tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Tensor/MTensorSpin3.h"
#include "support/models_test_support.hh"

// Test the irreducible spin-three projector, photon Ward identity and F-wave width
TEST_CASE("Spin-three current is transverse traceless and has the F-wave width", "[tensor][spin3][normalization]") {
  const double mass = 1.9, daughter = 0.2, width = 0.17, br = 0.23;
  const double momentum = std::sqrt(pow2(mass) / 4.0 - pow2(daughter));
  const double coupling = gra::MTensorPomeron::GDecay(3, mass, width, daughter, br);
  const gra::M4Vec q(0.23, -0.12, 1.3, 1.4);
  const gra::M4Vec zplus(0, 0, momentum, mass / 2), zminus(0, 0, -momentum, mass / 2);
  const auto axis = gra::tensor::Spin3Decay(zplus, zminus);
  const auto [nodes, weights] = gra::math::GaussLegendreRule(20, -1.0, 1.0);
  double average = 0.0;
  for (const auto &i : indices(nodes)) {
    const double z = nodes[i];
    const gra::M4Vec plus(momentum * std::sqrt(1 - z*z), 0, momentum*z, mass/2);
    const gra::M4Vec minus(-plus.Px(), 0, -plus.Pz(), mass/2);
    const auto p = plus + minus;
    const auto decay = gra::tensor::Spin3Decay(plus, minus);
    const auto current = gra::tensor::Spin3Current(q, p, decay);
    double amplitude = 0.0, norm = 0.0;
    for (int a = 0; a < 4; ++a) {
      double trace = 0.0;
      for (int b = 0; b < 4; ++b) {
        trace += (b == 0 ? 1 : -1) * decay(a,b,b);
        double ward = 0.0, transverse = 0.0;
        for (int c = 0; c < 4; ++c) {
          CHECK(decay(a,b,c) == Approx(decay(c,a,b)).margin(1e-13));
          ward += q[c] * current(c,a,b);
          transverse += p[c] * decay(c,a,b);
          amplitude += decay(a,b,c) * axis(a,b,c);
          norm += axis(a,b,c) * axis(a,b,c);
        }
        CHECK(std::abs(ward) < 1e-12);
        CHECK(std::abs(transverse) < 1e-12);
      }
      CHECK(std::abs(trace) < 1e-12);
    }
    amplitude *= coupling / std::sqrt(norm);
    const double expected = coupling * std::sqrt(2.0/5.0) * std::pow(2*momentum,3) * (5*z*z*z-3*z)/2;
    CHECK(amplitude == Approx(expected).margin(1e-12));
    average += weights[i] * amplitude * amplitude / 2;
  }
  CHECK(gra::kinematics::PDW2body(mass*mass, daughter*daughter, daughter*daughter, average, 1.0)
        == Approx(width*br).epsilon(2e-12));
}

// Test complex coherence and covariance of the full DS plus rho, omega and rho3 current
TEST_CASE("Spin-three PHOTO interferes coherently and respects beam and rotation symmetry", "[tensor][spin3][photo][phase]") {
  auto lts = AsymmetricCentralPairLTSForTest(211, -211, 1.69, -0.38);
  for (const int pdg : {113,223,117}) {
    gra::PARAM_RES res;
    res.p = lts.PDG.FindByPDG(pdg);
    SetToyTensorChannel(res, pdg == 117 ? std::vector<double>{0.43} : std::vector<double>{0.43,-0.18}).exchange = {22,995};
    res.hel_decay.g_decay_TP = {1.0};
    lts.process.RESONANCES.emplace(std::to_string(pdg), res);
  }
  const bool noflip = GENERATE(false,true);
  const auto tune = WriteModifiedPhotoVMTune("rho3_covariance_" + std::to_string(noflip), [&](auto &j) {
    j.at("PARAM_TENSORPOM").at("FORWARD_NOFLIP") = noflip;
  });
  gra::MTensorPomeron tensor(lts, gra::MModelTune::Load(tune.second),
      gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo), gra::MTensorPomeronMode::Photo);
  const auto evaluate = [&](gra::LORENTZSCALAR event) {
    event.screening.active = false;
    event.amplitude.BeginCentral();
    REQUIRE(tensor.MEPhoto(event) > 0);
    return std::vector<std::complex<double>>(event.hamp.begin(),event.hamp.end());
  };
  const auto total = evaluate(lts);
  // Keep the central decay cached while evaluating shifted screening transfers
  auto cached = lts;
  cached.amplitude.BeginCentral();
  REQUIRE(tensor.MEPhoto(cached) > 0);
  const gra::M4Vec shift(0.013, -0.009, 0, 0);
  cached.pfinal[1] -= shift;
  cached.pfinal[2] += shift;
  RefreshToyDerivedKinematicsPreserveDecay(cached);
  cached.screening.active = true;
  REQUIRE(tensor.MEPhoto(cached) > 0);
  const auto shifted = std::vector<std::complex<double>>(cached.hamp.begin(), cached.hamp.end());
  RequireVectorNear(shifted, evaluate(cached), 3e-10);
  REQUIRE(gra::SquaredNorm(shifted) != Approx(gra::SquaredNorm(total)).epsilon(1e-6));
  auto base = lts, spin3 = lts, continuum = lts;
  base.process.RESONANCES.erase("117");
  spin3.process.RESONANCES.erase("113");
  spin3.process.RESONANCES.erase("223");
  continuum.process.RESONANCES.clear();
  const auto b = evaluate(base), r = evaluate(spin3), ds = evaluate(continuum);
  for (const auto &i : indices(total)) { RequireComplexNear(total[i], b[i]+r[i]-ds[i],3e-10); }
  std::swap(lts.decaytree[0],lts.decaytree[1]);
  RequireVectorNear(evaluate(lts), total, 3e-10);
  const double norm = gra::SquaredNorm(total);
  for (const auto &event : {RotateToyEventAroundZ(lts,0.73), ReflectToyEventInXZ(lts), BeamExchangeMirrorWithDecay(lts)}) {
    CHECK(gra::SquaredNorm(evaluate(event)) == Approx(norm).epsilon(3e-8));
  }
  if (!noflip) {
    CHECK(gra::SquaredNorm(evaluate(BoostToyEventAlongZ(lts,-0.31))) == Approx(norm).epsilon(3e-8));
  }
}
