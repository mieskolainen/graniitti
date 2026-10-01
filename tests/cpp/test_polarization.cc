// Direct production polarization and electric photon polarization tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Regge/MReggeMPInit.h"
#include "support/models_test_support.hh"

namespace {

using gra::aux::indices;

// Read a physical resonance with normalized complex spin weights
gra::PARAM_RES ReadSpinWeights(const std::string& name, bool jz0, bool random = false,
                               const std::string& dynamics = "auto_min_L") {
  const auto directory = std::filesystem::path("tmp") / "polarization";
  std::filesystem::create_directories(directory);
  auto card =
      nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "RES/" + name + ".json")));
  const int J = card["PARAM_RES"]["spinX2"].get<int>() / 2;
  // Remove common mass suppression when testing spin identities far from a pole
  card["PARAM_RES"]["MODELS"]["MP"]["FF_prod"] = {{"type", "none"}};
  for (auto& [key, block] : card["PARAM_RES"]["MODELS"]["MP"].items()) {
    if (key.front() != '[') { continue; }
    block.erase("g_ls");
    block.erase("helicity");
    block["g"] = {key == "[22,22]" ? nlohmann::json(nullptr) : nlohmann::json(1.0), 0.0};
    block["basis"] = dynamics;
    block["Lambda"] = 1.0;
    auto& polarization         = block["polarization"];
    polarization["mode"]       = random ? "rho" : "a_Jz";
    polarization["random_rho"] = random;
    polarization["a_Jz"]       = nlohmann::json::array();
    for (int m = -J; m <= 0; ++m) {
      const double magnitude = jz0 ? (m == 0 ? 1.0 : 0.0) : (m % 2 == 0 ? 1.0 : 0.0);
      polarization["a_Jz"].push_back({m, magnitude, 0.31 * m * m});
    }
  }
  const auto path = directory / (name + (random ? "_random.json" : jz0 ? "_jz0.json" : "_mixed.json"));
  {
    std::ofstream output(path);
    output << card;
  }
  gra::MRandom rng;
  return gra::resonance::Read(std::filesystem::relative(path, gra::ResolveModelTuneDir("TUNE0")).string(), rng, gra::ReggeProductionModel::MP,
                              "TUNE0");
}

// Build an elastic pp event in the small-x photon limit
gra::LORENTZSCALAR PhotonPoint(double xi, double phi) {
  auto lts = MakeToyCoherentPhotonLTS();
  gra::RequireModelCache(lts.model_cache, gra::MModelTune::Load(modelfile), "photon polarization");
  constexpr double pz     = 6500.0;
  const double     energy = std::hypot(pz, gra::PDG::mp);
  lts.pbeam1              = gra::M4Vec(0.0, 0.0, pz, energy);
  lts.pbeam2              = gra::M4Vec(0.0, 0.0, -pz, energy);
  REQUIRE(gra::kinematics::BuildForwardParticleXiT(lts.pbeam1, xi, -0.0001, phi, true, lts.pfinal[1]));
  REQUIRE(gra::kinematics::BuildForwardParticleXiT(lts.pbeam2, xi, -0.0002, phi + 0.71, false, lts.pfinal[2]));
  // Keep the central pair above threshold at the smallest photon fraction
  lts.decaytree[0].p = lts.PDG.FindByPDG(11);
  lts.decaytree[1].p = lts.PDG.FindByPDG(-11);
  UpdateToyDerivedKinematics(lts);
  return lts;
}

}  // namespace

// Check the physical density against S R S^dagger using the complete fusion amplitude
TEST_CASE("MP spin filters preserve full production dynamics", "[polarization][MP][density]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  process.SetResonances({{"f2_1270", ReadSpinWeights("f2_1270", true)}});
  process.InitializeProcessAmplitude();
  auto event = AsymmetricF2PhasePointForTest(0.91, -0.38);
  event.process = process.state.lts.process;
  auto res = event.process.RESONANCES.at("f2_1270");
  gra::HelVec state = {{0.2, 0.3}, 0.0, {0.4, -0.2}, 0.0, {0.2, 0.3}};
  gra::Scale(state, 1.0 / std::sqrt(gra::SquaredNorm(state)));
  for (const std::string frame : {"CM", "CS", "HX"}) {
    for (const bool noflip : {true, false}) {
      event.process.MP_FRAME = frame;
      event.process.FORWARD_NOFLIP = noflip;
      const auto fusion = gra::rspin::Resonance(event, res, gra::mpom::Fusion, 1.0).front();
      const auto density = fusion.Transpose() * fusion.Conj();
      const auto pure = gra::RankOneProjector(state);
      const auto mixed = pure * 0.6 + gra::HelAmp::IdentityMatrix(5) * 0.08;
      for (const auto& rho : {pure, mixed, gra::HelAmp::IdentityMatrix(5) / 5.0}) {
        const auto root = rho.PrincipalPositiveSemidefiniteSquareRoot();
        res.MP.filter = root.Transpose() * std::sqrt(5.0);
        const auto amplitude = gra::mpom::Resonance(event, res, 1.0).front();
        const auto actual = amplitude.Transpose() * amplitude.Conj();
        const auto expected = (root * density * root) * 5.0;
        CAPTURE(frame, noflip);
        RequireMatrixNear(actual, expected, 2e-11);
      }
      res.MP.filter = pure.Transpose() * std::sqrt(5.0);
      const auto amplitude = gra::mpom::Resonance(event, res, 1.0).front();
      const auto actual = amplitude.Transpose() * amplitude.Conj();
      RequireMatrixNear(actual / actual.Trace(), pure, 2e-12);
      res.MP.filter = gra::HelAmp::IdentityMatrix(5);
      RequireMatrixNear(gra::mpom::Resonance(event, res, 1.0).front(), fusion, 2e-12);
    }
  }
}

// Rotate the physical state with the decay axes and integrate the full decay solid angle
TEST_CASE("MP direct states obey frame covariance and decay normalization",
          "[polarization][MP][frame][normalization]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  process.SetResonances({{"f2_1270", ReadSpinWeights("f2_1270", true)}});
  process.InitializeProcessAmplitude();
  auto event                   = AsymmetricF2PhasePointForTest(0.91, -0.38);
  event.process                = process.state.lts.process;
  event.process.FORWARD_NOFLIP = true;
  auto res                    = event.process.RESONANCES.at("f2_1270");
  const auto  decay_cm         = gra::spin::ResonanceDecayMatrix(event, res, "CM");
  gra::HelVec state            = {{0.3, 0.1}, 0.0, {-0.2, 0.4}, 0.0, {0.3, 0.1}};
  gra::Scale(state, 1.0 / std::sqrt(gra::SquaredNorm(state)));
  event.process.MP_FRAME = "CM";
  res.MP.filter = gra::RankOneProjector(state).Transpose() * std::sqrt(5.0);
  const auto reference = gra::mpom::Resonance(event, res, 1.0).front() * decay_cm;
  for (const std::string frame : {"CM", "CS", "HX"}) {
    event.process.MP_FRAME = frame;
    const auto rotation    = gra::spin::ProductionRotation(event, frame, 2.0);
    const auto rotated     = rotation * state;
    res.MP.filter = gra::RankOneProjector(rotated).Transpose() * std::sqrt(5.0);
    const auto production = gra::mpom::Resonance(event, res, 1.0).front();
    RequireMatrixNear(production * rotation.Conj() * decay_cm, reference, 3e-11);
  }
  const auto [node, weight] = gra::math::GaussLegendreRule(12, -1.0, 1.0);
  for (const gra::HelVec& input : {state, gra::HelVec{0.0, 0.0, 1.0, 0.0, 0.0}, gra::HelVec{1.0, 0.0, 0.0, 0.0, 0.0}}) {
    res.MP.filter = gra::RankOneProjector(input).Transpose() * std::sqrt(5.0);
    const auto production = gra::mpom::Resonance(event, res, 1.0).front();
    double     integral   = 0.0;
    for (const auto i : indices(node)) {
      for (int phi = 0; phi < 16; ++phi) {
        const auto decay =
            gra::spin::fDecayMatrix(res.hel_decay, std::acos(node[i]), 2.0 * gra::math::PI * phi / 16).Transpose();
        integral += weight[i] * (production * decay).FrobNorm2() / 32.0;
      }
    }
    // Integral dOmega/(4 pi) B B^dagger = ||T||^2 I/(2J+1)
    CHECK(integral == Approx(production.FrobNorm2() * res.hel_decay.T.FrobNorm2() / 5.0).epsilon(2e-12));
  }
}

// Distinguish coherent preparation from a density without an overall quantum phase
TEST_CASE("MP spin filters preserve their interference policy",
          "[polarization][MP][interference]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES+CON", "pi+ pi-");
  process.SetResonances({{"f2_1270", ReadSpinWeights("f2_1270", true)}});
  process.InitializeProcessAmplitude();
  auto base    = AsymmetricF2PhasePointForTest(0.91, -0.38);
  base.process = process.state.lts.process;
  base.hamp.Configure(process.state.lts.hamp.metadata);
  const auto evaluate = [&](gra::LORENTZSCALAR event, gra::MReggeMode mode) {
    gra::MRegge regge(event, process.state.model_tune, gra::MRegge::ProcessDefinitionFor(mode, "direct_density"));
    return regge.Amp2(event, gra::ReggeProductionModel::MP, mode);
  };
  const auto   total     = gra::MReggeMode::ResonanceContinuumTwoBody;
  const auto   pole      = gra::MReggeMode::Resonance;
  const double continuum = evaluate(base, gra::MReggeMode::ContinuumTwoBody);
  const double coherent  = evaluate(base, pole);
  auto         shifted   = base;
  shifted.process.RESONANCES.at("f2_1270").MP.phi += gra::math::PI / 2.0;
  CHECK(evaluate(shifted, pole) == Approx(coherent).epsilon(2e-12));
  CHECK(std::abs(evaluate(shifted, total) - evaluate(base, total)) > 1e-6 * coherent);
  auto& res           = base.process.RESONANCES.at("f2_1270");
  res.spin_basis      = "rho";
  res.rho = gra::RankOneProjector(res.a_Jz);
  CHECK(evaluate(base, pole) == Approx(coherent).epsilon(2e-12));
  CHECK(evaluate(base, total) == Approx(coherent + continuum).epsilon(2e-12));
  const gra::HelVec first  = {{0.2, 0.1}, 0.0, {0.4, -0.2}, 0.0, {0.2, 0.1}};
  const gra::HelVec second = {0.0, {0.1, 0.2}, 0.0, {-0.1, -0.2}, 0.0};
  res.rho = gra::RankOneProjector(first) + gra::RankOneProjector(second);
  res.rho = res.rho / res.rho.Trace();
  gra::mpom::PreparePolarization(res, base);
  const double mixture = evaluate(base, total);
  res.MP.phi += 0.7;
  CHECK(evaluate(base, total) == Approx(mixture).epsilon(2e-12));
}

// Reject incompatible states and forbidden production channels during initialization
TEST_CASE("MP direct polarization validates physical production sectors at initialization",
          "[polarization][MP][validation]") {
  for (const std::string invalid : {"photon_C", "basis", "CS_odd", "CM_coherence"}) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
    auto res = ReadSpinWeights("f2_1270", true);
    if (invalid == "photon_C") { res.MP.channels.front().exchange[0] = 22; }
    if (invalid == "basis") { res.MP.channels.front().basis = static_cast<gra::ReggeVertexBasis>(999); }
    if (invalid == "CS_odd") { res.a_Jz = {0.0, 0.5, std::sqrt(0.5), -0.5, 0.0}; }
    if (invalid == "CM_coherence") {
      process.SetMPFrame("CM");
      res.a_Jz = {0.5, 0.0, std::sqrt(0.5), 0.0, 0.5};
    }
    process.SetResonances({{"f2_1270", res}});
    CAPTURE(invalid);
    CHECK_THROWS_AS(process.InitializeProcessAmplitude(), std::invalid_argument);
  }
}

// Sample only physical random densities and reconstruct them after frame preparation
TEST_CASE("MP random densities respect frame symmetries", "[polarization][MP][density][random]") {
  for (const std::string frame : {"CM", "CS", "HX"}) {
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
    process.SetMPFrame(frame);
    process.SetResonances({{"f2_1270", ReadSpinWeights("f2_1270", true, true)}});
    REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
    auto event      = AsymmetricF2PhasePointForTest(0.91, -0.38);
    event.process   = process.state.lts.process;
    const auto& res = event.process.RESONANCES.at("f2_1270");
    CHECK(gra::spin::Positivity(res.rho, 2.0));
    if (frame == "CM") { CHECK(res.rho.IsDiagonal()); }
    const auto filter = res.MP.filter.Transpose();
    RequireMatrixNear(filter * filter / 5.0, res.rho, 2e-12);
    CHECK(filter.IsHermitian(1e-12));
  }
}

// Compare the scalar direct vertex with its L=S=0 fixed-pole limit
TEST_CASE("MP direct scalar normalization matches fusion", "[polarization][MP][normalization]") {
  auto event                  = MakeToyProductionLTSAsymmetric(0.3, 4.2, -4.6);
  event.process.SPINGEN       = true;
  auto fusion                 = MakeToyScalarMPResonance();
  auto direct                 = fusion;
  direct.spin_basis           = "a_Jz";
  direct.a_Jz                 = {1.0};
  gra::mpom::PreparePolarization(direct, event);
  for (const bool noflip : {true, false}) {
    event.process.FORWARD_NOFLIP = noflip;
    const auto scalar            = gra::mpom::Resonance(event, direct, 1.0).front();
    const auto pole              = gra::mpom::Resonance(event, fusion, 1.0).front();
    RequireMatrixNear(scalar, pole, 2e-12);
  }
}

// Compare complex screened amplitudes to a linear combination of physical spin filters
TEST_CASE("MP spin filters commute with soft screening", "[polarization][MP][screening]") {
  const auto model = gra::ReggeProductionModel::MP;
  const auto mode = gra::MReggeMode::Resonance;
  auto eikonal = BuildTestEikonalWithInitialState({0.31}, "spin_filter", ProtonInitialState(), 0.25, 1, 4, 4,
                                                 MakeToyCoherentPhotonLTS().s, 0.13);
  ToyHelicityProcess setup;
  ConfigureToyProductionProcess(setup, "MP", "RES", "pi+ pi-");
  setup.SetResonances({{"f2_1270", ReadSpinWeights("f2_1270", true)}});
  setup.InitializeProcessAmplitude();
  for (const std::string frame : {"CM", "CS", "HX"}) {
    ToyReggePhaseScreeningProcess process(eikonal.ModelTuneHandle(), model, mode, frame);
    process.eikonal = eikonal;
    process.state.lts.process.RESONANCES = setup.state.lts.process.RESONANCES;
    auto& res = process.state.lts.process.RESONANCES.at("f2_1270");
    const auto pure = res.MP.filter;
    for (const bool screened : {false, true}) {
      const auto evaluate = [&](const gra::HelAmp& filter) {
        res.MP.filter = filter;
        const double norm = screened ? process.ScreenedAmp2() : process.BornAmp2();
        REQUIRE(norm > 0.0);
        if (screened) { REQUIRE(process.PairTrace().size() > 1); }
        return process.state.lts.hamp;
      };
      const auto unit = evaluate(gra::HelAmp::IdentityMatrix(5));
      const auto projected = evaluate(pure);
      const auto mixed = evaluate(gra::HelAmp::IdentityMatrix(5) * 0.4 + pure * 0.6);
      auto expected = unit;
      for (const auto i : indices(expected)) { expected[i] = 0.4 * unit[i] + 0.6 * projected[i]; }
      CAPTURE(frame, screened);
      RequireVectorNear(mixed, expected, 3e-10);
    }
  }
}

// Check physical direct states through the complete MP production and decay amplitude
TEST_CASE("MP direct polarization preserves parity and its rank-one density", "[polarization][MP][parity]") {
  struct Channel {
    const char* name;
    const char* decay;
    int         first;
    int         second;
  };
  const std::array<Channel, 11> channels = {{{"f0_980", "pi+ pi-", 211, -211},
                                             {"f2_1270", "pi+ pi-", 211, -211},
                                             {"f4_2300", "pi+ pi-", 211, -211},
                                             {"f6_2510", "pi+ pi-", 211, -211},
                                             {"eta", "gamma gamma", 22, 22},
                                             {"eta2_1645", "a(2)(1320)0 pi0", 115, 111},
                                             {"chi_c1", "J/psi(1S)0 gamma", 443, 22},
                                             {"rho_770", "pi+ pi-", 211, -211},
                                             {"rho_770_odd", "pi+ pi-", 211, -211},
                                             {"eta_odd", "gamma gamma", 22, 22},
                                             {"f2_1270_yy", "pi+ pi-", 211, -211}}};
  for (const auto& channel : channels) {
    for (const std::string frame : {"CS", "HX", "CM"}) {
      for (const int sector : {0, 1, 2}) {
        const bool jz0 = sector == 0;
        if (frame == "CM" && !jz0) { continue; }
        ToyHelicityProcess process;
        ConfigureToyProductionProcess(process, "MP", "RES", channel.decay);
        process.SetMPFrame(frame);
        auto input = ReadSpinWeights(channel.name, jz0);
        if (sector == 2) {
          const int J = input.p.spinX2 / 2;
          if (J == 0) { continue; }
          input.a_Jz.assign(2 * J + 1, 0.0);
          input.a_Jz[J - 1] = std::sqrt(0.5);
          input.a_Jz[J + 1] = -std::sqrt(0.5);
        }
        process.SetResonances({{channel.name, input}});
        REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
        auto event = DirectCentralPairLTSForTest(channel.first, channel.second);
        // Use unequal transfers to keep axial Jz=0 away from its physical zero
        auto& lower = event.pfinal[2];
        lower.SetPxPyPzM(1.37 * lower.Px(), 1.37 * lower.Py(), lower.Pz(), gra::PDG::mp);
        event.pfinal[0] = event.pbeam1 + event.pbeam2 - event.pfinal[1] - lower;
        const auto pair = TwoBodyRestKinematics(event.pfinal[0].M(), event.decaytree[0].p.mass,
                                                event.decaytree[1].p.mass, 0.91, -0.38);
        for (const auto& i : indices(event.decaytree)) {
          event.decaytree[i].p4 = BoostFromRestFrame(pair[i], event.pfinal[0]);
        }
        RefreshToyDerivedKinematicsPreserveDecay(event);
        REQUIRE(std::abs(event.t1 - event.t2) > 1.0e-3);
        event.process = process.state.lts.process;
        event.hamp.Configure(process.state.lts.hamp.metadata);
        const auto evaluate = [&](gra::LORENTZSCALAR point) {
          gra::MRegge amplitude(point, process.state.model_tune,
                                gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "spin_weights"));
          return amplitude.Amp2(point, gra::ReggeProductionModel::MP, gra::MReggeMode::Resonance);
        };
        const double coherent = evaluate(event);
        CAPTURE(channel.name, frame, sector, coherent);
        REQUIRE(coherent > 0.0);
        CHECK(evaluate(ReflectToyEventInXZ(event)) == Approx(coherent).epsilon(2.0e-9));
        CHECK(evaluate(RotateToyEventAroundZ(event, 0.63)) == Approx(coherent).epsilon(2.0e-9));
        CHECK(evaluate(BeamExchangeMirrorWithDecay(event)) == Approx(coherent).epsilon(2.0e-9));
        auto&      resonance      = event.process.RESONANCES.at(channel.name);
        const auto weights        = resonance.a_Jz;
        resonance.spin_basis      = "rho";
        resonance.rho             = gra::RankOneProjector(weights);
        gra::mpom::PreparePolarization(resonance, event);
        CHECK(evaluate(event) == Approx(coherent).epsilon(2.0e-9));
      }
    }
  }
}

// Reconstruct the electric field independently from the JW photon components
TEST_CASE("EPA photon helicities give a field parallel to the transfer", "[polarization][EPA][QED]") {
  const gra::MDirac dirac("DIRAC");
  const auto        transitions = gra::spin::SpinHalfTransitions(true);
  for (const double phi : {0.0, 0.37, 1.29, -2.41}) {
    const auto lts = PhotonPoint(0.0001, phi);
    for (const int leg : {1, 2}) {
      const auto&      q    = leg == 1 ? lts.q1 : lts.q2;
      const double     sign = leg == 1 ? 1.0 : -1.0;
      const gra::M4Vec axis(0.0, 0.0, sign, 1.0);
      const auto       epsilon = dirac.MasslessSpin1States(axis, "conj", true);
      const auto source = gra::qed::PhotonSourceMatrixTransitions(lts, leg, transitions, {-1, 1}, "EPA", leg == 2);
      for (const auto& row : indices(transitions)) {
        for (std::size_t mu = 0; mu < 4; ++mu) {
          const auto   field    = source[row][0] * epsilon[0](mu) + source[row][1] * epsilon[1](mu);
          const double expected = mu == 1 ? q.Px() / q.Pt() : mu == 2 ? q.Py() / q.Pt() : 0.0;
          CAPTURE(leg, phi, row, mu, field);
          CHECK(std::abs(field - expected) < 2.0e-12);
        }
      }
    }
  }
}

// Match the photon density matrix to the electric limit of the Dirac-Pauli current
TEST_CASE("QED and EPA photon polarizations converge at small x", "[polarization][EPA][QED]") {
  const auto transitions = gra::spin::SpinHalfTransitions(true);
  for (const double phi : {0.0, 0.37, 1.29, -2.41}) {
    for (const int leg : {1, 2}) {
      double previous = 1.0;
      for (const double xi : {0.001, 0.0001, 0.00001}) {
        const auto lts = PhotonPoint(xi, phi);
        const auto epa = gra::qed::PhotonSourceMatrixTransitions(lts, leg, transitions, {-1, 1}, "EPA", leg == 2);
        const auto qed = gra::qed::PhotonSourceMatrixTransitions(lts, leg, transitions, {-1, 1}, "QED", leg == 2);
        double     difference = 0.0;
        for (const auto& row : indices(transitions)) {
          const auto reference = gra::RankOneProjector(epa.Row(row));
          const auto actual    = gra::RankOneProjector(qed.Row(row)) / gra::SquaredNorm(qed.Row(row));
          difference           = std::max(difference, std::sqrt((actual - reference).FrobNorm2()));
        }
        CAPTURE(leg, phi, xi, difference, previous);
        CHECK(difference < 3.0 * xi);
        CHECK(difference < previous);
        previous = difference;
      }
    }
  }
}

// Require vector photoproduction to cancel at zero transverse momentum in symmetric pp collisions
TEST_CASE("Vector photoproduction has destructive interference at zero transverse momentum",
          "[polarization][MP][XP][GP][interference]") {
  for (const std::string model : {"MP", "XP", "GP"}) {
    for (const std::string source : {"EPA", "QED"}) {
      ToyHelicityProcess process;
      ConfigureToyProductionProcess(process, model, "RES", "pi+ pi-");
      gra::MRandom rng;
      const auto   resonance = gra::resonance::Read("RES/rho_770.json", rng, gra::ParseReggeProductionModel(model), "TUNE0");
      process.SetResonances({{"rho_770", resonance}});
      REQUIRE_NOTHROW(process.InitializeProcessAmplitude());
      auto event                  = ScalarPolePhasePointForTest(0.02, 0.91, -0.38, resonance.p.mass, 6500.0);
      event.process               = process.state.lts.process;
      event.process.PHOTON_VERTEX = source;
      event.hamp.Configure(process.state.lts.hamp.metadata);
      const auto production_model = model == "MP"   ? gra::ReggeProductionModel::MP
                                    : model == "XP" ? gra::ReggeProductionModel::XP
                                                    : gra::ReggeProductionModel::GP;
      const auto evaluate         = [&](gra::LORENTZSCALAR point) {
        gra::MRegge amplitude(point, process.state.model_tune,
                              gra::MRegge::ProcessDefinitionFor(gra::MReggeMode::Resonance, "interference"));
        return amplitude.Amp2(point, production_model, gra::MReggeMode::Resonance);
      };
      const double coherent   = evaluate(event);
      double       incoherent = 0.0;
      const auto   channels   = event.process.RESONANCES.at("rho_770").production;
      REQUIRE(channels.size() == 2);
      for (const auto& channel : channels) {
        auto single                                        = event;
        single.process.RESONANCES.at("rho_770").production = {channel};
        incoherent += evaluate(single);
      }
      CAPTURE(model, source, coherent, incoherent, coherent / incoherent);
      REQUIRE(incoherent > 0.0);
      CHECK(coherent < 1.0e-6 * incoherent);
      event.beam1 = event.PDG.FindByPDG(-gra::PDG::PDG_p);
      CHECK(evaluate(event) == Approx(2.0 * incoherent).epsilon(1.0e-6));
    }
  }
}

// Exercise both intrinsic parities at every integer spin through the actual production API
TEST_CASE("MP polarization supports either parity beyond even spin", "[polarization][MP][spin][frame]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  process.SetMPFrame("CS");
  auto setup            = process.CreateProcessSetup();
  auto event            = AsymmetricF2PhasePointForTest(0.91, -0.38);
  event.process         = process.state.lts.process;
  event.process.SPINGEN = true;
  for (const std::array<int, 2> exchanges : {std::array<int, 2>{995, 995}, std::array<int, 2>{995, 9993}}) {
    for (int J = 0; J <= 6; ++J) {
      for (const int parity : {-1, 1}) {
        auto res     = ReadSpinWeights("f2_1270", true);
        res.p.spinX2 = 2 * J;
        res.p.P      = parity;
        res.p.C      = setup.lts.PDG.FindByPDG(exchanges[0]).C * setup.lts.PDG.FindByPDG(exchanges[1]).C;
        res.MP.channels.front().exchange = exchanges;
        res.a_Jz.assign(2 * J + 1, 0.0);
        res.a_Jz[J]         = 1.0;
        gra::mpom::PreparePolarization(res, event);
        std::vector<gra::MParticle> legs;
        gra::RES_PRODUCTION         source;
        source.g = 1.0;
        for (const int pdg : res.MP.channels.front().exchange) {
          gra::MDecayBranch branch;
          branch.p = setup.lts.PDG.FindByPDG(pdg);
          branch.legs.resize(2);
          for (auto& daughter : branch.legs) { daughter.p = setup.lts.PDG.FindByPDG(2212); }
          setup.ProcessHelicityTree(branch, true, true);
          legs.push_back(branch.p);
          source.tree.push_back(branch);
        }
        source.pole = gra::mpom::PrepareResonance(setup, res, res.MP.channels.front(), legs);
        if (exchanges[0] != exchanges[1]) {
          const auto  reverse = gra::mpom::PrepareResonance(setup, res, res.MP.channels.front(), {legs[1], legs[0]});
          const auto& term    = source.pole->terms.front();
          const int   exponent =
              static_cast<int>(term.l) + (legs[0].spinX2 + legs[1].spinX2 - static_cast<int>(term.two_s)) / 2;
          RequireComplexNear(reverse.terms.front().coefficient, (exponent % 2 == 0 ? 1.0 : -1.0) * term.coefficient,
                             2e-12);
        }
        res.production         = {source};
        event.process.MP_FRAME = "CM";
        const auto reference   = gra::mpom::Resonance(event, res, 1.0).front();
        CAPTURE(J, parity, exchanges);
        REQUIRE(reference.IsFinite());
        REQUIRE(reference.FrobNorm2() > 0.0);
        for (const std::string frame : {"CM", "CS", "HX"}) {
          event.process.MP_FRAME = frame;
          const auto rotation    = gra::spin::ProductionRotation(event, frame, J);
          const auto state       = rotation * res.a_Jz;
          res.MP.filter = gra::RankOneProjector(state).Transpose() * std::sqrt(static_cast<double>(2 * J + 1));
          const auto amplitude = gra::mpom::Resonance(event, res, 1.0).front();
          RequireMatrixNear(amplitude * rotation.Conj(), reference, 3e-11);
          const auto density = amplitude.Transpose() * amplitude.Conj();
          RequireMatrixNear(density / density.Trace(), gra::RankOneProjector(state), 3e-11);
        }
      }
    }
  }
}

// Match unpolarized filters to fusion amplitudes for every dynamics and derivative choice
TEST_CASE("MP polarization dynamics preserve fusion normalization", "[polarization][MP][normalization][dynamics]") {
  ToyHelicityProcess process;
  ConfigureToyProductionProcess(process, "MP", "RES", "pi+ pi-");
  auto setup = process.CreateProcessSetup();
  auto event = AsymmetricF2PhasePointForTest(0.91, -0.38);
  event.process.SPINGEN = true;
  for (const std::string name : {"f0_980", "f2_1270", "eta", "rho_770", "f2_1270_yy"}) {
    for (const std::string dynamics : {"auto_min_L", "auto_min_S", "auto_equal_ls", "auto_equal_helicity"}) {
      for (const bool derivative : {false, true}) {
        auto res = ReadSpinWeights(name, true, false, dynamics);
        auto channel = res.MP.channels.front();
        setup.lts.process.DERIVATIVE_FACTOR = derivative;
        std::vector<gra::MParticle> legs;
        gra::RES_PRODUCTION production;
        for (const int pdg : channel.exchange) {
          gra::MDecayBranch branch;
          branch.p = setup.lts.PDG.FindByPDG(pdg);
          branch.legs.resize(2);
          for (auto& leg : branch.legs) { leg.p = setup.lts.PDG.FindByPDG(2212); }
          setup.ProcessHelicityTree(branch, true, true);
          legs.push_back(branch.p);
          production.tree.push_back(branch);
        }
        production.pole = gra::mpom::PrepareResonance(setup, res, channel, legs);
        CHECK(production.pole->derivative_factor == derivative);
        const auto reference = gra::mpom::Fusion(event, *production.pole);
        auto explicit_channel = channel;
        explicit_channel.basis = gra::ReggeVertexBasis::LS;
        explicit_channel.g_ls.Clear();
        for (const auto& row : gra::spin::CanonicalPoleOperators(res.p, legs[0], legs[1], true,
                                                               channel.C_symmetry, channel.P_symmetry)) {
          explicit_channel.g_ls.Set(row.coupling.l, row.coupling.two_s, 0.0);
        }
        for (const auto& term : production.pole->terms) {
          explicit_channel.g_ls.Set(term.l, term.two_s, term.coefficient);
        }
        const auto fusion = gra::mpom::PrepareResonance(setup, res, explicit_channel, legs);
        CAPTURE(name, dynamics, derivative);
        RequireMatrixNear(reference, gra::mpom::Fusion(event, fusion), 2e-12);
        const std::complex<double> scale(0.71, -0.23);
        channel.g *= scale;
        const auto scaled = gra::mpom::PrepareResonance(setup, res, channel, legs);
        RequireMatrixNear(gra::mpom::Fusion(event, scaled), reference * scale, 2e-12);
        res.production = {production};
        res.spin_basis = "rho";
        const auto n = static_cast<std::size_t>(res.p.spinX2 + 1);
        res.rho = gra::HelAmp::IdentityMatrix(n) / static_cast<double>(n);
        gra::mpom::PreparePolarization(res, event);
        const auto filtered = gra::mpom::Resonance(event, res, 1.0).front();
        res.spin_basis = "none";
        RequireMatrixNear(filtered, gra::mpom::Resonance(event, res, 1.0).front(), 2e-12);
      }
    }
  }
}
