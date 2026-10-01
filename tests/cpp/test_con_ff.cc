// Continuum form factor tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Own
#include "support/models_test_support.hh"

namespace {

// Store one evaluated continuum amplitude and its helicity sum
struct ContinuumAmplitude {
  double amp2 = 0.0;
  std::vector<std::complex<double>> helicity;
};

// Compute one power transfer factor normalized at zero
gra::regge::FFParam TransferPower(const double scale2) {
  return {gra::regge::FFType::Power, gra::regge::FFNorm::Zero, {scale2, 1.0}};
}

// Evaluate one two body continuum point with an immutable model tune
ContinuumAmplitude EvaluateContinuumFF(const ToyHelicityProcess &process,
                                       const gra::LORENTZSCALAR &point,
                                       const gra::MModelTunePtr &model_tune,
                                       gra::ReggeProductionModel model,
                                       const std::string &family,
                                       const std::size_t hadronic_legs = 2) {
  gra::LORENTZSCALAR lts = point;
  lts.process = process.state.lts.process;
  if (hadronic_legs > 2) {
    throw std::invalid_argument(
        "EvaluateContinuumFF: hadronic leg count must not exceed two");
  }
  if (hadronic_legs < 2) {
    const int strong_pdg = lts.process.CONT_PRODUCTION.front()[1];
    const auto selected_production = lts.process.CONT_PRODUCTION.front();
    const auto selected_tree = lts.process.CONT_PRODUCTIONTREE.front();
    const double selected_tu_sign = lts.process.CONT_TU_SIGN.front();
    lts.process.CONT_PRODUCTION = {selected_production};
    lts.process.CONT_PRODUCTIONTREE = {selected_tree};
    lts.process.CONT_TU_SIGN = {selected_tu_sign};
    lts.process.CONTINUUM_POLE.clear();
    SetToyContinuumExchangePair(lts, gra::PDG::PDG_gamma,
                                hadronic_legs == 1 ? strong_pdg
                                                   : gra::PDG::PDG_gamma);
    if (model == gra::ReggeProductionModel::MP ||
        model == gra::ReggeProductionModel::XP) {
      PrepareToyContinuumOperators(lts, model);
    } else if (model == gra::ReggeProductionModel::GP) {
      UseToyGPSubchannelHelicity(lts);
    }
  }
  gra::MRegge regge(lts, model_tune,
                    gra::MRegge::ProcessDefinitionFor(
                        gra::MReggeMode::ContinuumTwoBody, family + "_CON_FM"));
  ContinuumAmplitude out;
  out.amp2 = TestReggeCon(regge, lts, model);
  out.helicity.assign(lts.hamp.begin(), lts.hamp.end());
  return out;
}

// Compute the analytic meson side vertex factor for a selected leg topology
double ExpectedTransferFF(const gra::LORENTZSCALAR &point,
                          const double Lambda_M2,
                          const std::size_t hadronic_legs) {
  double value = 1.0;
  const auto ff = TransferPower(Lambda_M2);
  if (hadronic_legs >= 1) {
    value *= gra::regge::TransferFF(point.t2, ff);
  }
  if (hadronic_legs == 2) {
    value *= gra::regge::TransferFF(point.t1, ff);
  }
  return value;
}

} // namespace

// Check the transfer option uses the same proton Dirac function as Tensor Pomeron
TEST_CASE("Dirac transfer uses the configured proton electromagnetic form", "[form-factor][dirac]") {
  for (const std::string em : {"DIPOLE", "KELLY"}) {
    const auto ff = gra::regge::ReadFF({{"type", "dirac"}, {"norm", "zero"}, {"EM", em}}, "test");
    gra::form::ParamStore structure;
    structure.EM = em;
    for (const double t : {0.0, -0.01, -0.1, -0.5, -2.0}) {
      CHECK(gra::regge::TransferFF(t, ff) == Approx(gra::form::F1(t, structure)).epsilon(1e-13));
    }
  }
  REQUIRE_THROWS_AS(gra::regge::ReadFF({{"type", "dirac"}, {"norm", "pole"}, {"EM", "KELLY"}}, "test"), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::regge::ReadFF({{"type", "dirac"}, {"norm", "zero"}, {"EM", "invalid"}}, "test"), std::invalid_argument);
}

// Check shared baryon transfer steering at complex amplitude level in all Regge models
TEST_CASE("MP XP and GP baryon continua use the shared Dirac transfer", "[form-factor][dirac][MP][XP][GP]") {
  const std::string em = GENERATE("DIPOLE", "KELLY");
  const auto general = [&](auto &j) { j.at("PARAM_STRUCTURE").at("EM") = em == "DIPOLE" ? "KELLY" : "DIPOLE"; };
  const auto inactive = WriteModifiedPhotoVMTune("baryon_transfer_none_" + em, general, {}, [](auto &card) {
    SetContinuumField(card, "[2212,2212]", "FF_transfer", {{"type", "none"}});
  });
  const auto active = WriteModifiedPhotoVMTune("baryon_transfer_dirac_" + em, general, {}, [&](auto &card) {
    SetContinuumField(card, "[2212,2212]", "FF_transfer", {{"type", "dirac"}, {"norm", "zero"}, {"EM", em}});
  });
  const auto model = gra::MModelTune::Load(active.second);
  const auto baseline = gra::MModelTune::Load(inactive.second);
  auto point = DirectCentralPairLTSForTest(2212, -2212);
  RefreshToyDerivedKinematicsPreserveDecay(point);
  gra::form::ParamStore structure;
  structure.EM = em;
  const double factor = gra::form::F1(point.t1, structure) * gra::form::F1(point.t2, structure);
  for (const auto &[family, type] : std::array{
           std::pair{"MP", gra::ReggeProductionModel::MP}, std::pair{"XP", gra::ReggeProductionModel::XP},
           std::pair{"GP", gra::ReggeProductionModel::GP}}) {
    CAPTURE(family, em);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, family, "CON", "p+ p-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto with_ff = EvaluateContinuumFF(process, point, model, type, family);
    const auto without_ff = EvaluateContinuumFF(process, point, baseline, type, family);
    REQUIRE(without_ff.amp2 > 0.0);
    CHECK(with_ff.amp2 / without_ff.amp2 == Approx(factor * factor).epsilon(2e-11));
    REQUIRE(with_ff.helicity.size() == without_ff.helicity.size());
    for (const auto &i : indices(with_ff.helicity)) {
      RequireComplexNear(with_ff.helicity[i], factor * without_ff.helicity[i], 2e-11);
    }
  }
}

// Reject missing continuum form factors before evaluating any amplitude
TEST_CASE("Continuum cards require explicit transfer and offshell form factors",
          "[gra::MRegge][MTensorPomeron][form-factor][validation]") {
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    for (const std::string field : {"FF_transfer", "FF_offshell"}) {
      CAPTURE(model, field);
      const std::string exchange = model == "GP" ? "990" : "995";
      const auto path = WriteModifiedContinuumTune("missing_" + model + "_" + field, model, [&](auto &card) {
        card.at(exchange).at("[211,211]").erase(field);
      });
      const auto tune = gra::MModelTune::Load(path + "/GENERAL.json");
      if (model == "TP") {
        REQUIRE_THROWS_AS(gra::ReadTensorPomeronParam(*tune, LoadedPDGTable()), std::invalid_argument);
      } else {
        REQUIRE_THROWS_AS(gra::regge::ReadParam({211, -211}, LoadedPDGTable(), *tune), std::invalid_argument);
      }
    }
  }
}

TEST_CASE("MP XP and GP apply an enabled meson side vertex factor per leg",
          "[gra::MRegge][continuum][form-factor][physics]") {
  constexpr double first_scale2 = 0.50;
  constexpr double second_scale2 = 1.30;
  const auto first_tune = WriteModifiedPhotoVMTune(
      "continuum_transfer_first", {}, {}, [first_scale2](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "power"}, {"norm", "zero"}, {"Lambda2", first_scale2}, {"n", 1.0}});
      });
  const auto second_tune = WriteModifiedPhotoVMTune(
      "continuum_transfer_second", {}, {}, [second_scale2](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "power"}, {"norm", "zero"}, {"Lambda2", second_scale2}, {"n", 1.0}});
      });
  const auto inactive_tune =
      WriteModifiedPhotoVMTune("continuum_transfer_inactive", {}, {}, [](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "none"}});
      });
  const auto first_model = gra::MModelTune::Load(first_tune.second);
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  const auto inactive_model = gra::MModelTune::Load(inactive_tune.second);
  const gra::LORENTZSCALAR point =
      ScalarPolePhasePointForTest(0.27, 1.18, -0.47, 1.35);
  const auto first_ff = TransferPower(first_scale2);
  const auto second_ff = TransferPower(second_scale2);
  const double expected = gra::regge::TransferFF(point.t1, second_ff) /
                          gra::regge::TransferFF(point.t1, first_ff) *
                          gra::regge::TransferFF(point.t2, second_ff) /
                          gra::regge::TransferFF(point.t2, first_ff);
  const double inactive_expected = 1.0 /
                                   gra::regge::TransferFF(point.t1, first_ff) /
                                   gra::regge::TransferFF(point.t2, first_ff);

  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> cases = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};

  for (const auto &test : cases) {
    CAPTURE(test.family, point.t1, point.t2, expected, inactive_expected);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());

    const auto first = EvaluateContinuumFF(process, point, first_model,
                                           test.model, test.family);
    const auto second = EvaluateContinuumFF(process, point, second_model,
                                            test.model, test.family);
    const auto inactive = EvaluateContinuumFF(process, point, inactive_model,
                                              test.model, test.family);
    REQUIRE(first.amp2 > 0.0);
    REQUIRE(first.helicity.size() == second.helicity.size());
    REQUIRE(first.helicity.size() == inactive.helicity.size());
    CHECK(second.amp2 / first.amp2 ==
          Approx(gra::math::pow2(expected)).epsilon(2.0e-11));
    for (const auto &i : gra::aux::indices(first.helicity)) {
      RequireComplexNear(second.helicity[i], expected * first.helicity[i],
                         2.0e-11);
      RequireComplexNear(inactive.helicity[i],
                         inactive_expected * first.helicity[i], 2.0e-11);
    }
  }
}

TEST_CASE("Enabled continuum power factors follow strong leg topology",
          "[gra::MRegge][continuum][form-factor][QED][physics]") {
  constexpr double cutoff2 = 0.50;
  const auto active_tune = WriteModifiedPhotoVMTune(
      "continuum_transfer_topology_active", {}, {}, [cutoff2](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "power"}, {"norm", "zero"}, {"Lambda2", cutoff2}, {"n", 1.0}});
        SetToyPhotonContinuum(card);
      });
  const auto inactive_tune = WriteModifiedPhotoVMTune(
      "continuum_transfer_topology_inactive", {}, {}, [](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_transfer", {{"type", "none"}});
        SetToyPhotonContinuum(card);
      });
  const auto active_model = gra::MModelTune::Load(active_tune.second);
  const auto inactive_model = gra::MModelTune::Load(inactive_tune.second);
  const gra::LORENTZSCALAR point =
      ScalarPolePhasePointForTest(0.27, 1.18, -0.47, 1.35);

  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> cases = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};

  for (const auto &test : cases) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    for (std::size_t hadronic_legs = 0; hadronic_legs <= 2; ++hadronic_legs) {
      CAPTURE(test.family, hadronic_legs);
      const auto active = EvaluateContinuumFF(
          process, point, active_model, test.model, test.family, hadronic_legs);
      const auto inactive =
          EvaluateContinuumFF(process, point, inactive_model, test.model,
                              test.family, hadronic_legs);
      REQUIRE(active.amp2 > 0.0);
      REQUIRE(active.helicity.size() == inactive.helicity.size());
      const double expected = ExpectedTransferFF(point, cutoff2, hadronic_legs);
      CHECK(active.amp2 / inactive.amp2 ==
            Approx(gra::math::pow2(expected)).epsilon(3.0e-10));
      for (const auto &i : gra::aux::indices(active.helicity)) {
        RequireComplexNear(active.helicity[i], expected * inactive.helicity[i],
                           3.0e-10);
      }
    }
  }
}

// Check both vertex factors through complex continuum helicity amplitudes
TEST_CASE("MP XP and GP apply the off-shell factor at both internal-line vertices",
          "[gra::MRegge][continuum][form-factor][offshell][physics]") {
  constexpr double first_scale2  = 0.70;
  constexpr double second_scale2 = 1.60;
  const std::string second_type = GENERATE("power", "logexp");
  const nlohmann::json second_card = second_type == "power"
      ? nlohmann::json{{"type", "power"}, {"norm", "pole"}, {"Lambda2", second_scale2}, {"n", 1.0}}
      : nlohmann::json{{"type", "logexp"}, {"norm", "pole"}, {"Lambda2", second_scale2}, {"b", 1.0}};
  const auto first_tune = WriteModifiedPhotoVMTune(
      "continuum_offshell_first", {}, {}, [first_scale2](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_offshell", {{"type", "power"}, {"norm", "pole"}, {"Lambda2", first_scale2}, {"n", 1.0}});
      });
  const auto second_tune = WriteModifiedPhotoVMTune(
      "continuum_offshell_second", {}, {}, [&second_card](auto &card) {
        SetContinuumField(card, "[211,211]", "FF_offshell", second_card);
      });
  const auto first_model  = gra::MModelTune::Load(first_tune.second);
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  const gra::LORENTZSCALAR point = ScalarPolePhasePointForTest(
      0.27, 0.5 * gra::math::PI, 0.5 * gra::math::PI, 1.35);
  REQUIRE(point.t_hat == Approx(point.u_hat).epsilon(2.0e-12));

  const double pion_m2 = gra::math::pow2(point.decaytree[0].p.mass);
  const gra::regge::FFParam first_ff = {
      gra::regge::FFType::Power, gra::regge::FFNorm::Pole, {first_scale2, 1.0}};
  const auto second_ff = gra::regge::ReadFF(second_card, "test");
  const double expected = gra::math::pow2(
      gra::regge::FormFactor(point.t_hat, pion_m2, second_ff) /
      gra::regge::FormFactor(point.t_hat, pion_m2, first_ff));

  struct ModelCase {
    const char *family;
    gra::ReggeProductionModel model;
  };
  const std::array<ModelCase, 3> cases = {{
      {"MP", gra::ReggeProductionModel::MP},
      {"XP", gra::ReggeProductionModel::XP},
      {"GP", gra::ReggeProductionModel::GP},
  }};
  for (const auto &test : cases) {
    CAPTURE(test.family, point.t_hat, expected);
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, test.family, "CON", "pi+ pi-");
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    const auto first = EvaluateContinuumFF(process, point, first_model, test.model, test.family);
    const auto second = EvaluateContinuumFF(process, point, second_model, test.model, test.family);
    REQUIRE(first.amp2 > 0.0);
    REQUIRE(first.helicity.size() == second.helicity.size());
    CHECK(second.amp2 / first.amp2 == Approx(gra::math::pow2(expected)).epsilon(3.0e-10));
    for (const auto &i : gra::aux::indices(first.helicity)) {
      RequireComplexNear(second.helicity[i], expected * first.helicity[i], 3.0e-10);
    }
  }
}
