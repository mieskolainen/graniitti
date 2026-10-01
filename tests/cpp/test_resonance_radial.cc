// Resonance radial decay factor tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Own
#include "support/models_test_support.hh"

namespace {

// Evaluate one configured resonance at a fixed physical pion pair point
std::vector<std::complex<double>> EvaluateResonanceAmplitude(
    const ToyHelicityProcess &process, const gra::LORENTZSCALAR &point,
    gra::PARAM_RES resonance, const gra::ReggeProductionModel model,
    const std::string &family, const bool decay_barrier = true,
    const bool spin_decay = true,
    const gra::MModelTunePtr &model_tune = nullptr) {
  gra::LORENTZSCALAR lts = point;
  lts.process = process.state.lts.process;
  lts.process.DECAY_BARRIER = decay_barrier;
  lts.process.SPINDEC = spin_decay;
  gra::MRegge regge(lts, model_tune ? model_tune : process.GetModelTune(),
                    gra::MRegge::ProcessDefinitionFor(
                        gra::MReggeMode::Resonance, family + "_RES_RADIAL"));
  const double amp2 = TestReggeRes(regge, lts, resonance, model);
  REQUIRE(amp2 > 0.0);
  return {lts.hamp.begin(), lts.hamp.end()};
}

} // namespace

TEST_CASE("Decay LS barriers preserve amplitudes and the orthogonal LS trace",
          "[gra::spin][resonance][radial][physics]") {
  const std::array<double, 3> weights = {0.2, 0.3, 0.5};
  const std::array<double, 3> phases = {0.2, -0.4, 0.7};
  gra::HELMatrix hel;
  hel.T = gra::MMatrix<std::complex<double>>(3, 1, 0.0);

  for (std::size_t l = 0; l < weights.size(); ++l) {
    gra::MMatrix<std::complex<double>> basis(3, 1, 0.0);
    basis[l][0] = 1.0;
    const std::complex<double> alpha =
        std::polar(std::sqrt(weights[l]), phases[l]);
    hel.T += basis * alpha;
    hel.ls_components.push_back({l, 0, alpha, std::move(basis)});
  }

  constexpr double ratio = 0.42;
  const auto scaled = gra::spin::DecayLSHelicityMatrix(hel, ratio, true);
  for (const auto &l : gra::aux::indices(weights)) {
    RequireComplexNear(scaled[l][0], hel.T[l][0] * std::pow(ratio, l), 1.0e-14);
  }

  const double expected_intensity = weights[0] +
                                    weights[1] * gra::math::pow2(ratio) +
                                    weights[2] * gra::math::pow4(ratio);
  CHECK(gra::spin::DecayLSIntensity(hel, ratio, true) ==
        Approx(expected_intensity).margin(1.0e-14));
  CHECK(gra::spin::DecayLSIntensity(hel, ratio, false) == Approx(1.0));

  const auto pole = gra::spin::DecayLSHelicityMatrix(hel, 1.0, true);
  const auto disabled = gra::spin::DecayLSHelicityMatrix(hel, ratio, false);
  const auto threshold = gra::spin::DecayLSHelicityMatrix(hel, 0.0, true);
  for (const auto &l : gra::aux::indices(weights)) {
    RequireComplexNear(pole[l][0], hel.T[l][0], 1.0e-14);
    RequireComplexNear(disabled[l][0], hel.T[l][0], 1.0e-14);
    RequireComplexNear(threshold[l][0], l == 0 ? hel.T[l][0] : 0.0, 1.0e-14);
  }
  CHECK(gra::spin::DecayLSIntensity(hel, 0.0, true) ==
        Approx(weights[0]).margin(1.0e-14));
}

TEST_CASE("Generic Jacob Wick decay is independent of the MRegge barrier",
          "[gra::spin][resonance][radial][isolation]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  const auto input =
      gra::resonance::Read("RES/f2_1270.json", process.state.random, gra::ReggeProductionModel::MP);
  process.SetResonances({{"f2_1270", input}});
  REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
  const gra::PARAM_RES resonance = process.GetResonances().at("f2_1270");

  gra::LORENTZSCALAR enabled =
      ScalarPolePhasePointForTest(0.20, 1.1, -0.4, 1.55);
  enabled.process = process.state.lts.process;
  enabled.process.DECAY_BARRIER = true;
  gra::LORENTZSCALAR disabled = enabled;
  disabled.process.DECAY_BARRIER = false;
  const auto generic_enabled = gra::spin::ResonanceDecayMatrix(
      enabled, resonance, enabled.process.MP_FRAME);
  const auto generic_disabled = gra::spin::ResonanceDecayMatrix(
      disabled, resonance, disabled.process.MP_FRAME);
  RequireMatrixNear(generic_enabled, generic_disabled, 1e-12);

  const double q = gra::kinematics::DecayMomentum(std::sqrt(enabled.m2),
                                                  enabled.decaytree[0].p4.M(),
                                                  enabled.decaytree[1].p4.M());
  const double q0 = gra::kinematics::DecayMomentum(resonance.p.mass,
                                                   enabled.decaytree[0].p4.M(),
                                                   enabled.decaytree[1].p4.M());
  gra::HELMatrix scaled = resonance.hel_decay;
  scaled.T =
      gra::spin::DecayLSHelicityMatrix(resonance.hel_decay, q / q0, true);
  const auto mregge_scaled = gra::spin::ResonanceDecayMatrix(
      enabled, resonance, scaled, enabled.process.MP_FRAME);
  REQUIRE(MatrixDiffNorm2(mregge_scaled, generic_enabled) > 1e-12);
}

TEST_CASE("MP XP and GP amplitudes contain one scalar pair radial factor",
          "[gra::MRegge][resonance][radial][physics]") {
  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> cases = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};
  const double central_mass = 1.55;
  const gra::LORENTZSCALAR point =
      ScalarPolePhasePointForTest(0.20, 1.1, -0.4, central_mass);

  for (const auto &test : cases) {
    CAPTURE(test.family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES", "pi+ pi-");
    const auto input =
        gra::resonance::Read("RES/f2_1270.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    process.SetResonances({{"f2_1270", input}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::PARAM_RES physical = process.GetResonances().at("f2_1270");
    auto &physical_form = ReggeResonanceFormForTest(physical, test.model);
    physical_form.ff_prod = {};
    physical.hel_decay.ff_decay = {};
    physical.BW = gra::BreitWigner::FixedWidth;
    gra::PARAM_RES pole_at_event = physical;
    pole_at_event.p.mass = central_mass;

    const auto physical_amplitude = EvaluateResonanceAmplitude(
        process, point, physical, test.model, test.family);
    const auto pole_amplitude = EvaluateResonanceAmplitude(
        process, point, pole_at_event, test.model, test.family);
    REQUIRE(physical_amplitude.size() == pole_amplitude.size());

    const double mass_a = point.decaytree[0].p4.M();
    const double mass_b = point.decaytree[1].p4.M();
    const double q =
        gra::kinematics::DecayMomentum(central_mass, mass_a, mass_b);
    const double q0 =
        gra::kinematics::DecayMomentum(physical.p.mass, mass_a, mass_b);
    const double radial = gra::math::pow2(q / q0);
    const std::complex<double> line_ratio =
        gra::resonance::FixedWidthLineShape(gra::math::pow2(central_mass),
                                            physical.p.mass, physical.p.width) /
        gra::resonance::FixedWidthLineShape(gra::math::pow2(central_mass),
                                            pole_at_event.p.mass,
                                            pole_at_event.p.width);
    const std::complex<double> scale = radial * line_ratio;

    for (const auto &i : gra::aux::indices(physical_amplitude)) {
      RequireComplexNear(physical_amplitude[i], scale * pole_amplitude[i],
                         3.0e-10);
    }

    const auto barrier_disabled = EvaluateResonanceAmplitude(
        process, point, physical, test.model, test.family, false);
    for (const auto &i : gra::aux::indices(physical_amplitude)) {
      RequireComplexNear(physical_amplitude[i], radial * barrier_disabled[i],
                         3.0e-10);
    }

    gra::PARAM_RES running = physical;
    running.BW = gra::BreitWigner::RunningWidth;
    const auto running_enabled = EvaluateResonanceAmplitude(
        process, point, running, test.model, test.family, true);
    const auto running_disabled = EvaluateResonanceAmplitude(
        process, point, running, test.model, test.family, false);
    const double q_ratio = q / q0;
    const double profile_enabled =
        physical.p.mass / central_mass * q_ratio * gra::math::pow2(radial);
    const double profile_disabled = physical.p.mass / central_mass * q_ratio;
    const std::complex<double> running_scale =
        radial *
        gra::resonance::RunningWidthLineShape(
            point.m2, physical.p.mass, physical.p.width, profile_enabled) /
        gra::resonance::RunningWidthLineShape(
            point.m2, physical.p.mass, physical.p.width, profile_disabled);
    for (const auto &i : gra::aux::indices(running_enabled)) {
      RequireComplexNear(running_enabled[i],
                         running_scale * running_disabled[i], 3.0e-10);
    }

    const auto blind_enabled = EvaluateResonanceAmplitude(
        process, point, physical, test.model, test.family, true, false);
    const auto blind_disabled = EvaluateResonanceAmplitude(
        process, point, physical, test.model, test.family, false, false);
    for (const auto &i : gra::aux::indices(blind_enabled)) {
      RequireComplexNear(blind_enabled[i], radial * blind_disabled[i], 3.0e-10);
    }

    gra::PARAM_RES fallback_running = physical;
    fallback_running.BW = gra::BreitWigner::RunningWidth;
    fallback_running.hel_decay.ls_components.clear();
    gra::PARAM_RES fallback_fixed = fallback_running;
    fallback_fixed.BW = gra::BreitWigner::FixedWidth;
    const auto running_fallback = EvaluateResonanceAmplitude(
        process, point, fallback_running, test.model, test.family);
    const auto fixed_fallback = EvaluateResonanceAmplitude(
        process, point, fallback_fixed, test.model, test.family);
    for (const auto &i : gra::aux::indices(running_fallback)) {
      RequireComplexNear(running_fallback[i], fixed_fallback[i], 3.0e-10);
    }
  }
}

TEST_CASE("MP XP and GP apply the invariant mass factor at both vertices",
          "[gra::MRegge][resonance][form-factor][physics]") {
  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> cases = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};
  constexpr double first_cutoff = 0.80;
  constexpr double second_cutoff = 1.30;
  const gra::LORENTZSCALAR point =
      ScalarPolePhasePointForTest(0.20, 1.1, -0.4, 1.10);

  for (const auto &test : cases) {
    CAPTURE(test.family);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "RES", "pi+ pi-");
    const auto input =
        gra::resonance::Read("RES/f0_500.json", process.state.random, gra::ParseReggeProductionModel(test.family));
    process.SetResonances({{"f0_500", input}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    gra::PARAM_RES first = process.GetResonances().at("f0_500");
    gra::PARAM_RES second = first;
    first.BW = gra::BreitWigner::RunningWidth;
    second.BW = gra::BreitWigner::RunningWidth;
    auto &first_form = ReggeResonanceFormForTest(first, test.model);
    auto &second_form = ReggeResonanceFormForTest(second, test.model);
    first_form.ff_prod = {gra::regge::FFType::Gaussian,
                          gra::regge::FFNorm::Pole,
                          {first_cutoff * first_cutoff}};
    second_form.ff_prod = {gra::regge::FFType::Gaussian,
                           gra::regge::FFNorm::Pole,
                           {second_cutoff * second_cutoff}};
    first.hel_decay.ff_decay = first_form.ff_prod;
    second.hel_decay.ff_decay = second_form.ff_prod;
    const auto first_amplitude = EvaluateResonanceAmplitude(
        process, point, first, test.model, test.family);
    const auto second_amplitude = EvaluateResonanceAmplitude(
        process, point, second, test.model, test.family);
    REQUIRE(first_amplitude.size() == second_amplitude.size());

    const double vertex_ratio =
        gra::regge::MassFF(point.m2, gra::math::pow2(first.p.mass),
                           second_form.ff_prod) /
        gra::regge::MassFF(point.m2, gra::math::pow2(first.p.mass),
                           first_form.ff_prod);
    const double expected = gra::math::pow2(vertex_ratio);
    for (const auto &i : gra::aux::indices(first_amplitude)) {
      RequireComplexNear(second_amplitude[i], expected * first_amplitude[i],
                         3.0e-10);
    }
  }
}

TEST_CASE("Common resonance power factors act only in strong fusion",
          "[gra::MRegge][resonance][form-factor][QED][physics]") {
  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> models = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};
  struct TopologyCase {
    const char *card;
    double mass;
    std::size_t hadronic_legs;
  };
  const std::array<TopologyCase, 3> topologies = {{
      {"f0_500.json", 1.10, 2},
      {"rho_770.json", 0.90, 1},
      {"f2_1270_yy.json", 1.40, 0},
  }};
  constexpr double cutoff2 = 1.0;

  for (const auto &test : models) {
    for (const auto &topology : topologies) {
      CAPTURE(test.family, topology.card, topology.hadronic_legs);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, test.family, "RES", "pi+ pi-");
      const auto input = gra::resonance::Read(
          "RES/" + std::string(topology.card), process.state.random, gra::ParseReggeProductionModel(test.family));
      process.SetResonances({{"topology_resonance", input}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

      gra::PARAM_RES resonance =
          process.GetResonances().at("topology_resonance");
      // Probe the explicit resonance transfer factor
      // Avoid the exact two-source cancellation of symmetric vector photoproduction
      const gra::LORENTZSCALAR point = topology.hadronic_legs == 1
          ? AsymmetricF2PhasePointForTest(1.1, -0.4)
          : ScalarPolePhasePointForTest(0.20, 1.1, -0.4, topology.mass);
      auto &form = ReggeResonanceFormForTest(resonance, test.model);
      form.ff_transfer = gra::regge::ReadFF({{"type", "none"}}, "test transfer");
      const auto disabled_amplitude =
          EvaluateResonanceAmplitude(process, point, resonance, test.model,
                                     test.family, true, true);
      form.ff_transfer = {gra::regge::FFType::Power, gra::regge::FFNorm::Zero, {cutoff2, 1.0}};
      const auto enabled_amplitude =
          EvaluateResonanceAmplitude(process, point, resonance, test.model,
                                     test.family, true, true);
      REQUIRE(enabled_amplitude.size() == disabled_amplitude.size());

      double expected = 1.0;
      if (topology.hadronic_legs == 2) {
        const gra::regge::FFParam ff = {gra::regge::FFType::Power,
                                        gra::regge::FFNorm::Zero,
                                        {cutoff2, 1.0}};
        expected *= gra::regge::TransferFF(point.t1, ff) *
                    gra::regge::TransferFF(point.t2, ff);
      }
      for (const auto &i : gra::aux::indices(enabled_amplitude)) {
        RequireComplexNear(enabled_amplitude[i],
                           expected * disabled_amplitude[i], 4.0e-10);
      }
    }
  }
}

TEST_CASE("Resonance production and decay factors each enter once",
          "[gra::MRegge][resonance][form-factor][QED][physics]") {
  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> models = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};
  struct TopologyCase {
    const char *card;
    double mass;
    std::size_t hadronic_legs;
  };
  const std::array<TopologyCase, 3> topologies = {{
      {"f0_500.json", 1.10, 2},
      {"rho_770.json", 0.90, 1},
      {"f2_1270_yy.json", 1.40, 0},
  }};
  constexpr double first_cutoff = 0.80;
  constexpr double second_cutoff = 1.30;

  for (const auto &test : models) {
    for (const auto &topology : topologies) {
      CAPTURE(test.family, topology.card, topology.hadronic_legs);
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, test.family, "RES", "pi+ pi-");
      const auto input = gra::resonance::Read(
          "RES/" + std::string(topology.card), process.state.random, gra::ParseReggeProductionModel(test.family));
      process.SetResonances({{"topology_resonance", input}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

      gra::PARAM_RES first = process.GetResonances().at("topology_resonance");
      gra::PARAM_RES second = first;
      auto &first_form = ReggeResonanceFormForTest(first, test.model);
      auto &second_form = ReggeResonanceFormForTest(second, test.model);
      first_form.ff_prod = {gra::regge::FFType::Gaussian,
                            gra::regge::FFNorm::Pole,
                            {first_cutoff * first_cutoff}};
      second_form.ff_prod = {gra::regge::FFType::Gaussian,
                             gra::regge::FFNorm::Pole,
                             {second_cutoff * second_cutoff}};
      first.hel_decay.ff_decay = first_form.ff_prod;
      second.hel_decay.ff_decay = second_form.ff_prod;
      // Avoid the exact two-source cancellation of symmetric vector photoproduction
      const gra::LORENTZSCALAR point = topology.hadronic_legs == 1
          ? AsymmetricF2PhasePointForTest(1.1, -0.4)
          : ScalarPolePhasePointForTest(0.20, 1.1, -0.4, topology.mass);
      const auto first_amplitude = EvaluateResonanceAmplitude(
          process, point, first, test.model, test.family);
      const auto second_amplitude = EvaluateResonanceAmplitude(
          process, point, second, test.model, test.family);
      REQUIRE(first_amplitude.size() == second_amplitude.size());

      const double vertex_ratio =
          gra::regge::MassFF(point.m2, gra::math::pow2(first.p.mass),
                             second_form.ff_prod) /
          gra::regge::MassFF(point.m2, gra::math::pow2(first.p.mass),
                             first_form.ff_prod);
      const double expected = gra::math::pow2(vertex_ratio);
      for (const auto &i : gra::aux::indices(first_amplitude)) {
        RequireComplexNear(second_amplitude[i], expected * first_amplitude[i],
                           4.0e-10);
      }
    }
  }
}
