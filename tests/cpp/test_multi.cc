// Four and six particle Regge amplitude tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <numeric>

#include "Graniitti/Regge/MReggeMulti.h"
#include "support/models_test_support.hh"

// Build a controlled two-channel production tune with an off-diagonal screening Odderon
std::string MultiReggeTune(const std::string &name, const std::vector<gra::regge::Topology> &topologies,
                          const double q, const bool exponential = false) {
  const auto directory = std::filesystem::path("tmp") / "multiregge_tests" / name;
  std::filesystem::create_directories(directory);
  std::filesystem::copy(std::filesystem::path(modelfile).parent_path(), directory,
                        std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
  auto card = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  card["PARAM_SOFT"]["active_model"] = "double";
  auto &soft = card["PARAM_SOFT"]["MODEL"]["double"];
  soft["EIKONAL"]["q"] = q;
  soft["EIKONAL"]["unitarization"] = exponential ? "exp" : "q_exp";
  soft["EIKONAL"]["screening_exchanges"] = {"P", "O"};
  soft["EIKONAL"]["helicity"] = true;
  soft["EXCHANGE"]["O"]["helicity"]["kappa"] = 0.1;
  card["PARAM_SPIN"]["FORWARD_NOFLIP"] = false;
  soft["EXCHANGE"]["P"]["g"] = {{2.0, 0.0}, {0.0, 3.0}};
  soft["EXCHANGE"]["O"]["g"] = {{0.4, 0.15}, {0.15, 0.6}};
  soft["FF"]["P"] = {{"type", "EXP"}, {"param", {{4.0}, {5.0}}}};
  soft["FF"]["O"] = {{"type", "EXP"}, {"param", {{3.0}, {4.0}}}};
  auto &continuum = card["PARAM_REGGE"]["PARAM_CON"];
  for (const std::string model : {"MP", "XP", "GP"}) {
    const int pdg = model == "GP" ? 990 : 995;
    continuum[model]["[211,-211]"] = {{pdg, pdg}};
  }
  int count = 0;
  for (const int size : topologies.front()) { count += size; }
  continuum["MULTI"]["partitions"][std::to_string(count)] = topologies;
  std::ofstream(directory / "GENERAL.json") << card.dump(2);

  auto gp = nlohmann::json::parse(gra::aux::GetInputData((directory / "CON_GP.json").string()));
  for (const std::string sector : {"same", "opposite"}) {
    auto &rows = gp.at("990").at("[211,211]").at(sector).at("helicity");
    rows.erase(std::remove_if(rows.begin(), rows.end(), [](const auto &row) { return row[2].template get<int>() != 0; }), rows.end());
  }
  std::ofstream(directory / "CON_GP.json") << gp.dump(2);
  auto numerics = nlohmann::json::parse(gra::aux::GetInputData((directory / "NUMERICS.json").string()));
  auto &loop = numerics["NUMERICS_REGGE"]["LOOP_INTEGRAL"];
  loop["NumberKT"] = 3;
  loop["NumberPHI"] = 8;
  loop["MaxKT"] = 0.6;
  loop["NumberBT"] = 32;
  loop["MaxBT"] = 3.0;
  auto &eikonal = numerics["NUMERICS_EIKONAL"];
  eikonal["NumberBT"] = 128;
  eikonal["MaxBT"] = 4.0;
  eikonal["NumberKT2"] = 64;
  eikonal["MaxKT2"] = 1.0;
  eikonal["FBIntegralN"] = 256;
  eikonal["FBIntegralMaxKT"] = 5.0;
  eikonal["LOOP_INTEGRAL"]["NumberLoopKT"] = 2;
  eikonal["LOOP_INTEGRAL"]["NumberLoopPHI"] = 8;
  eikonal["LOOP_INTEGRAL"]["MaxLoopKT"] = 0.08;
  std::ofstream(directory / "NUMERICS.json") << numerics.dump(2);
  return (directory / "GENERAL.json").string();
}

// Generate one reproducible physical F phase-space point through the real process API
std::unique_ptr<gra::MFactorized> MultiReggeEvent(const std::string &tune, const std::string &family,
                                               const int count, const bool screening, std::string final_state = "") {
  auto model = gra::MModelTune::Load(tune);
  auto process = std::make_unique<gra::MFactorized>(family + "[CON]<F>", std::vector<gra::aux::OneCMD>{}, model);
  process->SetHelicityConfig(model);
  process->SetInitialState({"p+", "p+"}, {100.0, 100.0});
  if (final_state.empty()) {
    for (int pair = 0; pair < count / 2; ++pair) { final_state += "pi+ pi- "; }
  }
  process->SetDecayMode(final_state);
  process->SetScreening(screening);
  auto &cuts = process->state.gcuts;
  cuts.forward_pt_min = 0.1;
  cuts.forward_pt_max = 0.4;
  cuts.M_min = 3.0;
  cuts.M_max = 5.0;
  cuts.Y_min = -0.2;
  cuts.Y_max = 0.2;
  cuts.XI_min = 0.0;
  cuts.XI_max = 1.0;
  process->PrepareRun();
  process->state.random.SetSeed(78231);
  gra::MEventWeightState status;
  status.adaptation_mode = true;
  status.include_screening = screening;
  const double weight = process->EventWeight({0.31, 0.59, 0.17, 0.68, 0.45, 0.64}, status);
  INFO(status.amplitude_failure_message);
  REQUIRE(status.Valid());
  REQUIRE(std::isfinite(weight));
  REQUIRE(weight > 0.0);
  REQUIRE(gra::SquaredNorm(process->state.lts.hamp) > 0.0);
  return process;
}

// Check complex scalar amplitudes under collider symmetries and identical pion exchange
void MultiReggeSymmetries(gra::MFactorized &process) {
  const auto base = process.state.lts;
  auto permuted = base;
  std::swap(permuted.decaytree[0].p4, permuted.decaytree[2].p4);
  const std::vector<gra::LORENTZSCALAR> transformed = {
      RotateToyEventAroundZ(base, gra::math::PI / 2.0), ReflectToyEventInXZ(base),
      BeamExchangeMirrorWithDecay(base), BoostToyEventAlongZ(base, 0.43), permuted};
  for (const auto &symmetry : indices(transformed)) {
    CAPTURE(symmetry);
    process.state.lts = transformed[symmetry];
    for (const auto &i : indices(process.state.lts.decaytree)) {
      process.state.lts.pfinal[i + 3] = process.state.lts.decaytree[i].p4;
    }
    REQUIRE(gra::kinematics::SetLorentzScalars(process.state, base.decaytree.size() + 2));
    process.ProcPtr.GetBareAmplitude2(process.state.lts);
    RequireVectorNear(process.state.lts.hamp, base.hamp, 1.0e-8);
  }
}

// Compare physical Born amplitudes while changing only the absorption functional
TEST_CASE("MultiRegge Born ladders are independent of q eikonal absorption",
          "[gra::MRegge][multiregge][physics][screening]") {
  ModelParamRestoreGuard restore;
  const std::vector<gra::regge::Topology> topologies = {{4}, {2, 2}, {6}, {4, 2}, {2, 2, 2}};
  for (const auto &topology : topologies) {
    const int count = std::accumulate(topology.begin(), topology.end(), 0);
    const std::string label = std::to_string(count) + "_" + std::to_string(topology.size());
    const auto first = MultiReggeTune(label + "_q06", {topology}, 0.6);
    const auto second = MultiReggeTune(label + "_q08", {topology}, 0.8);
    const auto exponential = MultiReggeTune(label + "_exp", {topology}, 1.0, true);
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model, topology);
      const auto reference = MultiReggeEvent(first, model, count, false);
      for (const auto &tune : {second, exponential}) {
        const auto changed = MultiReggeEvent(tune, model, count, false);
        const auto &a = reference->state.lts.hamp;
        const auto &b = changed->state.lts.hamp;
        REQUIRE(a.size() == b.size());
        double delta = 0.0;
        for (const auto &i : indices(a)) { delta += std::norm(a[i] - b[i]); }
        CHECK(delta < 1.0e-22 * gra::SquaredNorm(a));
      }
      MultiReggeSymmetries(*reference);
    }
  }
}

// Apply noncommuting Pomeron/Odderon screening to a coherent sum of real production topologies
TEST_CASE("MultiRegge coherent topologies share one matrix absorption operator",
          "[gra::MRegge][multiregge][physics][screening][GoodWalker]") {
  ModelParamRestoreGuard restore;
  std::map<std::string, std::vector<std::complex<double>>> screened;
  for (const int count : {4, 6}) {
    const std::vector<gra::regge::Topology> topologies = count == 4
        ? std::vector<gra::regge::Topology>{{4}, {2, 2}}
        : std::vector<gra::regge::Topology>{{6}, {4, 2}, {2, 2, 2}};
    for (const double q : {0.6, 0.8}) {
      const std::string label = "screen_" + std::to_string(count) + "_" + std::to_string(q);
      const auto tune = MultiReggeTune(label, topologies, q);
      for (const std::string model : {"MP", "XP", "GP"}) {
        CAPTURE(model, q, count);
        const auto combined = MultiReggeEvent(tune, model, count, true);
        const auto &sum = combined->state.lts.hamp;
        const auto &loop = combined->eikonal.GetLoopConst(combined->state.lts.s);
        REQUIRE_FALSE(loop.pair_screening_pair_diagonal);
        REQUIRE_FALSE(loop.pair_screening_spin_scalar);
        std::vector<std::complex<double>> expected(sum.size(), 0.0);
        for (const auto &sector : indices(topologies)) {
          const auto part_tune = MultiReggeTune(label + "_" + std::to_string(sector), {topologies[sector]}, q);
          const auto part = MultiReggeEvent(part_tune, model, count, true);
          REQUIRE(part->state.lts.hamp.size() == sum.size());
          for (const auto &i : indices(sum)) { expected[i] += part->state.lts.hamp[i]; }
        }
        RequireVectorNear(sum, expected, 1.0e-8);
        const auto [previous, inserted] = screened.try_emplace(model + std::to_string(count), sum);
        if (!inserted) {
          double difference = 0.0;
          for (const auto &i : indices(sum)) { difference += std::norm(sum[i] - previous->second[i]); }
          CHECK(difference > 1.0e-16 * gra::SquaredNorm(sum));
        }
        const auto born = MultiReggeEvent(tune, model, count, false);
        CHECK(std::abs(gra::SquaredNorm(sum) - gra::SquaredNorm(born->state.lts.hamp)) >
              1.0e-8 * gra::SquaredNorm(born->state.lts.hamp));
      }
    }
  }
}

// Check that independent parallel blocks can carry different fixed-spin Pomerons
TEST_CASE("MultiRegge parallel blocks do not require a serial exchange chain",
          "[gra::MRegge][multiregge][physics][exchange]") {
  ModelParamRestoreGuard restore;
  for (const auto &topology : std::vector<gra::regge::Topology>{{2, 2}, {4, 2}, {2, 2, 2}}) {
    const int count = std::accumulate(topology.begin(), topology.end(), 0);
    const auto tune = MultiReggeTune("mixed_" + std::to_string(count) + "_" + std::to_string(topology.size()), {topology}, 1.0, true);
    auto card = nlohmann::json::parse(gra::aux::GetInputData(tune));
    for (const std::string model : {"MP", "XP"}) {
      card["PARAM_REGGE"]["PARAM_CON"][model]["[321,-321]"] = {{991, 991}};
    }
    card["PARAM_REGGE"]["PARAM_CON"]["GP"]["[321,-321]"] = {{990, 990}};
    std::ofstream(tune) << card.dump(2);
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model, topology);
      const auto process = MultiReggeEvent(tune, model, count, false,
                                          count == 4 ? "pi+ pi- K+ K-" : "pi+ pi- pi+ pi- K+ K-");
      CHECK_FALSE(process->state.lts.process.CONT_LADDER_PERMUTATIONS.empty());
    }
  }
}

// Match every scalar continuum vertex across the LS and helicity inputs at the spin two pole
std::string ScalarMultiTune(const std::string &name, const gra::regge::Topology &topology,
                           const bool ls, const std::complex<double> coupling) {
  const auto tune = MultiReggeTune(name, {topology}, 1.0, true);
  const auto &pdg = LoadedPDGTable();
  const auto pole = gra::spin::PreparePoleLS(pdg.FindByPDG(995), pdg.FindByPDG(211), pdg.FindByPDG(-211),
                                           {{2, 0, coupling}}, 1.0, true, true, true,
                                           gra::spin::VertexContext::SubTUChannelExchange, 0.0, false);
  const auto coefficient = ls ? coupling : gra::spin::PoleLSReduced(pole, 1.0)(0, 0);
  auto general = nlohmann::json::parse(gra::aux::GetInputData(tune));
  general["PARAM_SOFT"]["EXCHANGE_DEF"]["P"]["trajectory_mode"] = "linear";
  general["PARAM_SOFT"]["MODEL"]["double"]["EXCHANGE"]["P"]["alpha"] = {2.0, 0.0};
  std::ofstream(tune) << general.dump(2);
  for (const std::string model : {"MP", "XP", "GP"}) {
    const auto path = std::filesystem::path(tune).parent_path() / ("CON_" + model + ".json");
    auto card = nlohmann::json::parse(gra::aux::GetInputData(path.string()));
    auto &pair = card[model == "GP" ? "990" : "995"]["[211,211]"];
    auto row = nlohmann::json::array({ls ? 2 : 0, 0});
    if (model == "GP") { row.push_back(0); }
    row.push_back(std::abs(coefficient));
    row.push_back(std::arg(coefficient));
    for (const std::string sector : {"same", "opposite"}) {
      pair[sector] = {{"basis", ls ? "crossed_ls" : "crossed_helicity"}, {"CP", {true, true}},
                      {ls ? "g_ls" : "helicity", {row}}};
    }
    pair["FF_transfer"] = {{"type", "none"}};
    pair["FF_offshell"] = {{"type", "none"}};
    for (auto &vertex : card) {
      if (!vertex.contains("[211,211]")) { continue; }
      vertex["[211,211]"]["reggeize"] = {{"active", false}, {"freeze_scale2", 1.0}};
      vertex["[211,211]"]["pveto"] = {{"active", false}, {"M0", 1.0}, {"c", 1.0}};
    }
    std::ofstream(path) << card.dump(2);
  }
  return tune;
}

// Compare complete complex amplitudes and one factor of the coupling per central particle
TEST_CASE("MultiRegge MP XP GP amplitudes share scalar pole normalization and phases",
          "[gra::MRegge][multiregge][physics][normalization]") {
  ModelParamRestoreGuard restore;
  const auto coupling = std::polar(0.73, -0.28);
  const auto scale = std::polar(1.17, 0.31);
  for (const auto &topology : std::vector<gra::regge::Topology>{{4}, {2, 2}, {6}, {4, 2}, {2, 2, 2}}) {
    const int count = std::accumulate(topology.begin(), topology.end(), 0);
    const std::string label = "pole_" + std::to_string(count) + "_" + std::to_string(topology.size());
    const auto helicity = ScalarMultiTune(label + "_hel", topology, false, coupling);
    const auto ls = ScalarMultiTune(label + "_ls", topology, true, coupling);
    const auto scaled = ScalarMultiTune(label + "_scaled", topology, false, coupling * scale);
    const auto reference = MultiReggeEvent(helicity, "MP", count, false);
    auto expected = reference->state.lts.hamp;
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model, topology);
      for (const auto &tune : {helicity, ls}) {
        CAPTURE(tune);
        const auto process = MultiReggeEvent(tune, model, count, false);
        RequireVectorNear(process->state.lts.hamp, expected, 2.0e-10);
      }
      const auto process = MultiReggeEvent(scaled, model, count, false);
      auto scaled_expected = expected;
      for (auto &value : scaled_expected) { value *= std::pow(scale, count); }
      RequireVectorNear(process->state.lts.hamp, scaled_expected, 2.0e-10);
    }
  }
}

// Resolve serial diagrams independently and reevaluate cached amplitudes at screening transfers
TEST_CASE("MultiRegge serial diagrams sum coherently and caches follow screening transfers",
          "[gra::MRegge][multiregge][physics][cache]") {
  ModelParamRestoreGuard restore;
  for (const int count : {4, 6}) {
    const auto tune = MultiReggeTune("serial_sum_" + std::to_string(count), {{count}}, 1.0, true);
    for (const std::string family : {"MP", "XP", "GP"}) {
      CAPTURE(family, count);
      const auto model = gra::ParseReggeProductionModel(family);
      const auto process = MultiReggeEvent(tune, family, count, false);
      const auto base = process->state.lts;
      const auto &permutations = base.process.CONT_LADDER_PERMUTATIONS;
      CHECK(permutations.size() == (count == 4 ? 16 : 288));
      CHECK(process->state.symmetry_factor == Approx(gra::math::pow2(gra::math::factorial(count / 2))));
      std::vector<std::complex<double>> sum(base.hamp.size(), 0.0);
      const auto definition = gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::ContinuumTwoFourSixBody, "serial_sum");
      for (const auto &permutation : permutations) {
        auto lts = base;
        lts.process.CONT_LADDER_PERMUTATIONS = {permutation};
        gra::MRegge diagram(lts, process->state.model_tune, definition);
        diagram.Amp2(lts, model, gra::MReggeMode::ContinuumTwoFourSixBody);
        REQUIRE(lts.hamp.size() == sum.size());
        for (const auto &i : indices(sum)) { sum[i] += lts.hamp[i]; }
      }
      RequireVectorNear(base.hamp, sum, 2.0e-10);
      for (const std::array<double, 2> shift : {std::array<double, 2>{0.07, -0.03}, {-0.04, 0.11}, {0.0, 0.0}}) {
        CAPTURE(shift);
        process->state.lts = base;
        process->state.lts.pfinal_orig = base.pfinal;
        REQUIRE(gra::kinematics::RebuildScreeningKinematics(
            process->state.lts, {base.pfinal[1].Px() - shift[0], base.pfinal[1].Py() - shift[1]},
            {base.pfinal[2].Px() + shift[0], base.pfinal[2].Py() + shift[1]}, true));
        REQUIRE(gra::kinematics::SetLorentzScalars(process->state, count + 2, true));
        auto fresh = process->state.lts;
        gra::MRegge amplitude(fresh, process->state.model_tune, definition);
        amplitude.Amp2(fresh, model, gra::MReggeMode::ContinuumTwoFourSixBody);
        process->ProcPtr.GetBareAmplitude2(process->state.lts);
        RequireVectorNear(process->state.lts.hamp, fresh.hamp, 2.0e-10);
      }
    }
  }
}

// Check rotations away from the quadrature axes after resolving the parallel loop integrals
TEST_CASE("MultiRegge amplitudes preserve arbitrary collider azimuths",
          "[gra::MRegge][multiregge][physics][rotation]") {
  ModelParamRestoreGuard restore;
  for (const auto &topology : std::vector<gra::regge::Topology>{{4}, {2, 2}, {6}, {4, 2}, {2, 2, 2}}) {
    const int count = std::accumulate(topology.begin(), topology.end(), 0);
    const auto tune = MultiReggeTune("rotation_" + std::to_string(count) + "_" + std::to_string(topology.size()), {topology}, 1.0, true);
    const auto path = std::filesystem::path(tune).parent_path() / "NUMERICS.json";
    auto numerics = nlohmann::json::parse(gra::aux::GetInputData(path.string()));
    auto &loop = numerics["NUMERICS_REGGE"]["LOOP_INTEGRAL"];
    loop["NumberKT"] = 32;
    loop["NumberPHI"] = 64;
    loop["MaxKT"] = 2.25;
    loop["NumberBT"] = 96;
    loop["MaxBT"] = 7.0;
    std::ofstream(path) << numerics.dump(2);
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model, topology);
      const auto process = MultiReggeEvent(tune, model, count, false);
      const auto base = process->state.lts;
      process->state.lts = RotateToyEventAroundZ(base, 0.37);
      REQUIRE(gra::kinematics::SetLorentzScalars(process->state, count + 2));
      process->ProcPtr.GetBareAmplitude2(process->state.lts);
      RequireVectorNear(process->state.lts.hamp, base.hamp, 2.0e-6);
    }
  }
}

// Include even and odd C secondary trajectories in every serial and parallel topology
TEST_CASE("MultiRegge secondary exchanges preserve collider symmetries",
          "[gra::MRegge][multiregge][physics][exchange][parity]") {
  ModelParamRestoreGuard restore;
  for (const auto &topology : std::vector<gra::regge::Topology>{{4}, {2, 2}, {6}, {4, 2}, {2, 2, 2}}) {
    const int count = std::accumulate(topology.begin(), topology.end(), 0);
    const auto tune = MultiReggeTune("secondary_" + std::to_string(count) + "_" + std::to_string(topology.size()), {topology}, 1.0, true);
    auto card = nlohmann::json::parse(gra::aux::GetInputData(tune));
    auto &continuum = card["PARAM_REGGE"]["PARAM_CON"];
    continuum["MULTI"]["secondary_exchanges"] = true;
    for (const std::string model : {"MP", "XP", "GP"}) {
      const int pomeron = model == "GP" ? 990 : 995;
      const int even = model == "GP" ? 9910 : 9915;
      const int odd = model == "GP" ? 9930 : 9933;
      continuum[model]["[211,-211]"] = {{pomeron, pomeron}, {pomeron, even}, {pomeron, odd},
                                        {even, even}, {even, odd}, {odd, odd}};
    }
    std::ofstream(tune) << card.dump(2);
    for (const std::string model : {"MP", "XP", "GP"}) {
      CAPTURE(model, topology);
      const auto process = MultiReggeEvent(tune, model, count, false);
      MultiReggeSymmetries(*process);
    }
  }
}
