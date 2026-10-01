// Forward excitation vertices and particle production tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "HepMC3/GenParticle.h"
#include "HepMC3/GenVertex.h"
#include "tests/cpp/support/models_test_support.hh"

namespace {

// Write a complete tune with short quadratures for event conservation checks
std::string ForwardTune(const std::string &model, double scale = 1.0) {
  const auto tune = WriteModifiedSudakovModelTune(
      "forward_" + model + "_" + std::to_string(scale),
      [&](auto &j) {
        j["PARAM_SOFT"]["active_model"] = model;
        j["PARAM_DURHAM"]["JET_pt_min"] = 0.0;
        auto &s0                        = j["PARAM_SOFT"]["FORWARD_EXCITATION"]["s0"];
        s0                              = scale * s0.template get<double>();
      },
      [](auto &j) {
        auto &eikonal                             = j["NUMERICS_EIKONAL"];
        eikonal["NumberBT"]                       = 128;
        eikonal["NumberKT2"]                      = 128;
        eikonal["FBIntegralN"]                    = 2048;
        eikonal["LOOP_INTEGRAL"]["NumberLoopKT"]  = 4;
        eikonal["LOOP_INTEGRAL"]["NumberLoopPHI"] = 4;
        j["NUMERICS_DURHAM"]["N_qt"]              = 6;
        j["NUMERICS_DURHAM"]["N_phi"]             = 4;
        for (const auto name : {"SHUV", "SUDA"}) {
          j["NUMERICS_SUDAKOV"][name]["N"]     = {32, 32};
          j["NUMERICS_SUDAKOV"][name]["N_int"] = 64;
        }
      });
  for (const auto file : {"DECAYS.json", "PDG_EXTRA.json"}) {
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", file), std::filesystem::path(tune) / file,
                               std::filesystem::copy_options::overwrite_existing);
  }
  return tune;
}

// Check charge, baryon number and color closure of a generated forward branch
void CheckForward(const gra::MDecayBranch &branch, const std::string &fragment) {
  if (fragment == "none") {
    CHECK(branch.legs.empty());
    return;
  }
  REQUIRE(branch.legs.size() >= 2);
  int        charge = 0, baryon3 = 0;
  gra::M4Vec sum;
  for (const auto &leg : branch.legs) {
    charge += leg.p.chargeX3;
    const int pdg = std::abs(leg.p.pdg), sign = gra::math::sign(leg.p.pdg);
    if (pdg == 2212 || pdg == 2112) { baryon3 += 3 * sign; }
    if (pdg == 1 || pdg == 2) { baryon3 += sign; }
    if (pdg == 2101 || pdg == 2203 || pdg == 1103) { baryon3 += 2 * sign; }
    sum += leg.p4;
    CHECK(leg.p4.M2() == Approx(leg.p.mass * leg.p.mass).margin(2.0e-9));
  }
  CHECK(charge == 3);
  CHECK(baryon3 == 3);
  CHECK(gra::math::CheckEMC(sum - branch.p4));
  if (fragment == "diquark") {
    REQUIRE(branch.legs.size() == 2);
    CHECK(branch.legs[0].p.color_flow.flow1 > 0);
    CHECK(branch.legs[0].p.color_flow.flow1 == branch.legs[1].p.color_flow.flow2);
  }
}

}  // namespace

// Rapidity sampling must preserve parity without assigning species to ordered slots
TEST_CASE("Cylinder daughters have no particle order or longitudinal parity bias", "[forward][fragmentation][parity]") {
  gra::MModelCache cache(gra::MModelTune::Load(modelfile));
  const auto       param = gra::GetNstarParam(cache);
  gra::MRandom     random;
  random.SetSeed(7919);
  for (const auto &masses : std::vector<std::vector<double>>{{gra::PDG::mp, gra::PDG::mpi0},
                                                             {gra::PDG::mpi0, gra::PDG::mp},
                                                             {gra::PDG::mp, gra::PDG::mpi, gra::PDG::mpi}}) {
    const gra::M4Vec          mother(0.0, 0.0, 0.0, 4.0);
    std::vector<unsigned int> positive(masses.size());
    unsigned int              accepted = 0;
    for (unsigned int trial = 0; trial < 512; ++trial) {
      std::vector<gra::M4Vec> p;
      const auto             &c = param->cylinder;
      if (gra::MFragment::TubeFragment(mother, mother.M(), masses, p, std::pow(masses.size(), c.q_power), c.T,
                                       c.exp_lambda, c.max_pt, random, *param) < 0.0) {
        continue;
      }
      ++accepted;
      gra::M4Vec sum;
      for (const auto &i : indices(p)) {
        positive[i] += p[i].Pz() > 0.0;
        sum += p[i];
        REQUIRE(p[i].M2() == Approx(masses[i] * masses[i]).margin(1.0e-10));
      }
      REQUIRE(gra::math::CheckEMC(sum - mother));
    }
    REQUIRE(accepted > 500);
    for (const auto count : positive) {
      CHECK(std::abs(static_cast<double>(count) - 0.5 * accepted) < 3.0 * std::sqrt(accepted));
    }
  }
}

// Reject power-law spectrum settings which have no sampleable transverse momentum
TEST_CASE("Cylinder input rejects empty power-law spectra", "[forward][input]") {
  const auto model                    = gra::MModelTune::Load(modelfile);
  auto       card                     = model->General("PARAM_NSTAR");
  card["CYLINDER"]["pt_distribution"] = "powexp";
  card["CYLINDER"]["q_power"]         = 0.04;
  card["CYLINDER"]["pt_bins"]         = 64;
  gra::MNstarParam param;
  REQUIRE_NOTHROW(param.Configure(card, "test"));
  for (const double power : {0.0, 1.0e-30, 1.0e30}) {
    auto invalid                   = card;
    invalid["CYLINDER"]["q_power"] = power;
    CHECK_THROWS_AS(param.Configure(invalid, "test"), std::invalid_argument);
  }
  card["CYLINDER"]["pt_bins"] = 1;
  CHECK_THROWS_AS(param.Configure(card, "test"), std::invalid_argument);
  card["CYLINDER"]["pt_distribution"] = "exp";
  card["CYLINDER"]["q_power"]         = 0.0;
  CHECK_NOTHROW(param.Configure(card, "test"));
}

// Exercise actual amplitudes, screening, fragmentation and HepMC3 vertices together
TEST_CASE("Forward particle production conserves both remnants across models", "[forward][process][screening][HepMC]") {
  const std::string fragment   = GENERATE("none", "fewbody", "cylinder", "diquark");
  const std::string model      = GENERATE("single", "double");
  const std::string channel    = GENERATE("MP[CON]<F> -> pi+ pi-", "XP[CON]<F> -> pi+ pi-", "GP[CON]<F> -> pi+ pi-", "TP[CON]<F> -> pi+ pi-",
                                          "yy[EPA]<F> -> mu+ mu-", "gg[QCD]<F> -> g g");
  const int         excitation = GENERATE(1, 2);
  CAPTURE(fragment, model, channel, excitation);
  ModelParamRestoreGuard restore;
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["MODELPARAM"]    = ForwardTune(model);
  card["GENERIC"]["HIST"]          = 0;
  card["GENERIC"]["CORES"]         = 1;
  card["GENERIC"]["RNDSEED"]       = 7919;
  card["GENERIC"]["INTEGRATOR"]    = "VEGAS";
  card["SCATTERING"]["PROCESS"]    = channel;
  card["SCATTERING"]["BEAMFRAG"]   = fragment;
  card["SCATTERING"]["RES"]        = nlohmann::json::array();
  card["SCATTERING"]["ENERGY"]     = {100.0, 100.0};
  card["SCATTERING"]["NSTARS"]     = excitation;
  card["SCATTERING"]["LOOPSCREEN"] = true;
  card["SCATTERING"]["LHAPDF"]     = "MMHT2014lo68cl";
  card["GENCUTS"]["<F>"] = {{"M", {4.0, 5.0}}, {"Rap", {-0.5, 0.5}}, {"Pt", {0.1, 0.5}}, {"Xi", {0.0001, 0.001}}};
  card["FIDCUTS"]        = {{"active", false}};
  card["VETOCUTS"]       = {{"active", false}};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto &process = *generator.proc;
  if (channel.starts_with("TP") && model == "double") {
    CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("multichannel Tensor forward excitation"));
    return;
  }
  process.PrepareRun();
  gra::MRandom random;
  random.SetSeed(3571);
  unsigned int accepted = 0;
  for (unsigned int trial = 0; trial < 100 && accepted < 2; ++trial) {
    std::vector<double> point(process.GetdLIPSDim());
    for (auto &unit : point) { unit = random.U(0.0, 1.0); }
    std::array<gra::M4Vec, 3> born;
    bool                      supported = true;
    for (const bool screening : {false, true}) {
      process.state.random.SetSeed(7919 + trial);
      gra::MEventWeightState weight;
      weight.include_screening = screening;
      const double value       = process.EventWeight(point, weight);
      REQUIRE_FALSE(weight.technical_failure);
      if (!(value > 0.0)) {
        supported = false;
        break;
      }
      const auto &lts = process.state.lts;
      REQUIRE(lts.excite1 + lts.excite2 == excitation);
      if (lts.excite1) { CheckForward(lts.decayforward1, fragment); }
      if (lts.excite2) { CheckForward(lts.decayforward2, fragment); }
      for (const auto &i : indices(born)) {
        if (!screening) { born[i] = lts.pfinal[i]; }
        CHECK(gra::math::CheckEMC(lts.pfinal[i] - born[i]));
      }
      REQUIRE_FALSE(lts.screening.active);
      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      REQUIRE(process.EventRecord(event));
      int remnants = 0;
      for (const auto &particle : event.particles()) {
        if (std::abs(particle->pid()) != gra::PDG::PDG_NSTAR) { continue; }
        ++remnants;
        if (fragment == "none") {
          CHECK_FALSE(particle->end_vertex());
        } else {
          REQUIRE(particle->end_vertex());
          CHECK(particle->end_vertex()->particles_out().size() >= 2);
        }
      }
      CHECK(remnants == excitation);
      for (const auto &vertex : event.vertices()) {
        gra::M4Vec balance;
        for (const auto &particle : vertex->particles_in()) { balance += gra::aux::HepMC2M4Vec(particle->momentum()); }
        for (const auto &particle : vertex->particles_out()) { balance -= gra::aux::HepMC2M4Vec(particle->momentum()); }
        CHECK(gra::math::CheckEMC(balance));
      }
      CHECK(process.state.lts.decayforward1.legs.empty());
      CHECK(process.state.lts.decayforward2.legs.empty());
    }
    accepted += supported;
  }
  REQUIRE(accepted == 2);
}

// Vary only the excitation scale to detect duplicate vertices and phase changes
TEST_CASE("Durham excitation scales each complex amplitude once per forward leg", "[forward][vertex][phase]") {
  ModelParamRestoreGuard restore;
  const auto             tune   = ForwardTune("double");
  const auto             scaled = ForwardTune("double", 4.0);
  const auto             first  = gra::MModelTune::Load(tune + "/GENERAL.json");
  const auto             second = gra::MModelTune::Load(scaled + "/GENERAL.json");
  const double           alpha  = first->Soft()->Exchange(first->Soft()->ForwardExcitationExchange()).Alpha0();
  for (const auto &[upper, lower] : std::array{std::pair{true, false}, std::pair{false, true}, std::pair{true, true}}) {
    auto lts    = MakeToyDurhamGG();
    lts.excite1 = upper;
    lts.excite2 = lower;
    for (const int side : {1, 2}) {
      if (!(side == 1 ? upper : lower)) { continue; }
      auto &p = lts.pfinal[side];
      p.SetPxPyPzM(p.Px(), p.Py(), p.Pz(), 1.4);
    }
    const auto central   = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
    const auto daughters = TwoBodyRestKinematics(central.M(), 0.0, 0.0, 0.31, 0.57);
    for (const auto &i : indices(lts.decaytree)) { lts.decaytree[i].p4 = BoostFromRestFrame(daughters[i], central); }
    lts.forward_mass2 = {lts.pfinal[1].M2(), lts.pfinal[2].M2()};
    UpdateToyDurhamDerivedKinematics(lts);
    lts.pfinal_orig   = lts.pfinal;
    auto other        = lts;
    other.model_cache = std::make_shared<gra::MModelCache>(second);
    gra::MRandom random;
    gra::MDurham a(lts, first, random, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham b(other, second, random, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham::DurhamProjectedAmp hard(1);
    hard[0] = {std::complex<double>(0.3, -0.2), {0.1, 0.4}, {-0.5, 0.2}, {0.2, -0.3}};
    for (const bool shifted : {false, true}) {
      if (shifted) {
        for (auto *point : {&lts, &other}) {
          REQUIRE(gra::kinematics::RebuildScreeningKinematics(
              *point, {point->pfinal[1].Px() - 0.08, point->pfinal[1].Py() + 0.03},
              {point->pfinal[2].Px() + 0.08, point->pfinal[2].Py() - 0.03}, true));
          UpdateToyDurhamDerivedKinematics(*point);
        }
      }
      const auto amplitude = a.DQtloopAmplitudes(lts, hard);
      const auto changed   = b.DQtloopAmplitudes(other, hard);
      REQUIRE(gra::SquaredNorm(amplitude) > 0.0);
      REQUIRE(changed.size() == amplitude.size());
      const double ratio = std::pow(4.0, 0.5 * alpha * (upper + lower));
      for (const auto &i : indices(amplitude)) { RequireComplexNear(changed[i], ratio * amplitude[i], 1.0e-11); }
    }
  }
}

// Resolve target rules separately from photon emission and reject ambiguous selectors
TEST_CASE("Dissociation steering separates photon emission and hadronic targets", "[forward][dissociation][input]") {
  auto card = gra::MModelTune::Load(modelfile)->General("PARAM_NSTAR");
  card["MODEL"] = {{"MP", {{"[*]", "soft"}, {"[22]", "structure"}, {"[22,P]", "hera"}}},
                   {"ygg", {{"[22]", "structure"}, {"[22,P]", "hera"}}},
                   {"yy", "structure"}, {"X", "triple_regge"}};
  gra::MNstarParam param;
  param.Configure(card, "dissociation test");
  CHECK(param.Model("MP").photo == gra::DissociationType::Hera);
  CHECK(param.Model("MP").hadron == gra::DissociationType::Soft);
  CHECK(param.Model("ygg").photo == gra::DissociationType::Hera);
  CHECK(param.Model("ygg").hadron == gra::DissociationType::None);
  CHECK(param.Model("yy").hadron == gra::DissociationType::Structure);
  CHECK(param.Model("X").hadron == gra::DissociationType::TripleRegge);
  CHECK(param.Model("gg").hadron == gra::DissociationType::None);
  for (const auto &invalid : {nlohmann::json("soft"), nlohmann::json{{"MP", "heraa"}}, nlohmann::json{{"MP", 1}},
      nlohmann::json{{"MP", {{"[*]", "soft"}, {"[22]", "hera"}}}},
      nlohmann::json{{"MP", {{"[*]", "soft"}, {"[22,22]", "hera"}}}},
      nlohmann::json{{"MP", {{"[22,P]", "hera"}}}},
      nlohmann::json{{"ygg", {{"[22]", "structure"}}}},
      nlohmann::json{{"ygg", {{"[22,P]", "heraa"}}}},
      nlohmann::json{{"gg", {{"[*]", "soft"}, {"[22,P]", "hera"}}}}}) {
    card["MODEL"] = invalid;
    CHECK_THROWS_AS(param.Configure(card, "invalid dissociation"), std::invalid_argument);
  }
}

// Reject invalid assignments in the common initialization before constructing amplitudes
TEST_CASE("Every excitation process validates its dissociation assignment", "[forward][dissociation][input]") {
  const gra::MSubProc registry;
  for (const auto &entry : registry.CreateAllProcesses()) {
    const auto models = entry->DissociationModels();
    if (models.empty()) { continue; }
    const std::string invalid = models.front() == gra::DissociationType::Structure ? "soft" : "structure";
    const auto tune = WriteModifiedPhotoVMTune("diss_invalid", [&](auto &card) {
      card["PARAM_NSTAR"]["MODEL"][entry->ISTATE] = invalid;
    });
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "yy", "EPA", "mu+ mu-");
    process.SetModelTune(gra::MModelTune::Load(tune.second));
    process.SetExcitation(1);
    CAPTURE(entry->Label());
    CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("PARAM_NSTAR.MODEL"));
  }
  const auto hera = WriteModifiedPhotoVMTune("diss_unsupported", [](auto &card) { card["PARAM_NSTAR"]["MODEL"]["MP"] = "hera"; });
  ToyHelicityProcess continuum;
  ConfigureToyProductionProcess(continuum, "MP", "CON", "pi+ pi-");
  continuum.SetModelTune(gra::MModelTune::Load(hera.second));
  continuum.SetExcitation(1);
  CHECK_THROWS_WITH(continuum.PrepareRun(), Catch::Contains("PARAM_NSTAR.MODEL"));
  for (const std::string key : {"MMP", "MP[RESS]", "yy_LUX", "IPp"}) {
    const auto tune = WriteModifiedPhotoVMTune("diss_unknown", [&](auto &card) { card["PARAM_NSTAR"]["MODEL"][key] = "soft"; });
    ToyHelicityProcess process;
    ConfigureToyProductionProcess(process, "yy", "EPA", "mu+ mu-");
    process.SetModelTune(gra::MModelTune::Load(tune.second));
    CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("PARAM_NSTAR.MODEL"));
  }
}

// Check proton-only NSTARS selection through generation and the complete HepMC record
TEST_CASE("Lepton proton NSTARS leaves the charged lepton intact", "[forward][lepton][process][HepMC]") {
  const std::string lepton = GENERATE("e-", "e+");
  const bool reverse = GENERATE(false, true);
  const int excitation = GENERATE(1, 2);
  CAPTURE(lepton, reverse, excitation);
  ModelParamRestoreGuard restore;
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["MODELPARAM"] = ForwardTune("single");
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["PROCESS"] = "yy[EPA]<F> -> mu+ mu-";
  card["SCATTERING"]["RES"] = nlohmann::json::array();
  card["SCATTERING"]["BEAM"] = reverse ? std::vector<std::string>{"p+", lepton} : std::vector<std::string>{lepton, "p+"};
  card["SCATTERING"]["ENERGY"] = reverse ? std::vector<double>{100.0, 25.0} : std::vector<double>{25.0, 100.0};
  card["SCATTERING"]["NSTARS"] = excitation;
  card["SCATTERING"]["BEAMFRAG"] = "fewbody";
  card["SCATTERING"]["LOOPSCREEN"] = false;
  card["GENCUTS"]["<F>"] = {{"M", {4.0, 5.0}}, {"Rap", {-0.5, 0.5}}, {"Pt", {0.1, 0.5}}, {"Xi", {0.0002, 0.001}}};
  card["FIDCUTS"] = {{"active", false}};
  card["VETOCUTS"] = {{"active", false}};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto &process = *generator.proc;
  if (excitation == 2) {
    CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("NSTARS"));
    return;
  }
  process.PrepareRun();
  gra::MRandom random;
  random.SetSeed(3571);
  unsigned int accepted = 0;
  for (unsigned int trial = 0; trial < 100 && accepted < 2; ++trial) {
    std::vector<double> point(process.GetdLIPSDim());
    for (auto &unit : point) { unit = random.U(0.0, 1.0); }
    gra::MEventWeightState weight;
    if (!(process.EventWeight(point, weight) > 0.0)) { continue; }
    const auto &lts = process.state.lts;
    REQUIRE(lts.excite1 == reverse);
    REQUIRE(lts.excite2 == !reverse);
    const auto &beam = reverse ? lts.beam2 : lts.beam1;
    const auto outgoing = lts.pfinal[reverse ? 2 : 1];
    const auto remnant = lts.pfinal[reverse ? 1 : 2];
    CHECK(outgoing.M2() == Approx(beam.mass * beam.mass).margin(1.0e-10));
    CHECK((reverse ? lts.decayforward2 : lts.decayforward1).legs.empty());
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    REQUIRE(process.EventRecord(event));
    int leptons = 0, remnants = 0;
    for (const auto &particle : event.particles()) {
      if (particle->status() == 1 && particle->pid() == beam.pdg) {
        ++leptons;
        CHECK_FALSE(particle->end_vertex());
        CHECK(gra::math::CheckEMC(gra::aux::HepMC2M4Vec(particle->momentum()) - outgoing));
      }
      if (std::abs(particle->pid()) == gra::PDG::PDG_NSTAR) {
        ++remnants;
        CHECK(gra::math::CheckEMC(gra::aux::HepMC2M4Vec(particle->momentum()) - remnant));
      }
    }
    CHECK(leptons == 1);
    CHECK(remnants == 1);
    ++accepted;
  }
  REQUIRE(accepted == 2);
}

// Reject absent and unsupported beam fragmentation modes during input parsing
TEST_CASE("Beam fragmentation is required process steering", "[forward][input][beamfrag]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["SCATTERING"]["PROCESS"] = "yy[EPA]<F> -> mu+ mu-";
  card["SCATTERING"]["RES"] = nlohmann::json::array();
  for (const auto &value : {nlohmann::json("string"), nlohmann::json("invalid"), nlohmann::json(1)}) {
    card["SCATTERING"]["BEAMFRAG"] = value;
    gra::MGraniitti generator;
    CHECK_THROWS(generator.ReadInput(card));
  }
  card["SCATTERING"].erase("BEAMFRAG");
  gra::MGraniitti generator;
  CHECK_THROWS(generator.ReadInput(card));
  auto nstar = gra::MModelTune::Load(modelfile)->General("PARAM_NSTAR");
  nstar["fragment"] = "none";
  gra::MNstarParam param;
  CHECK_THROWS_AS(param.Configure(nstar, "obsolete fragmentation selector"), std::invalid_argument);
}

// Check forward decay normalization constraints after selecting the process mode
TEST_CASE("Forward three body decays require unweighting", "[forward][input][beamfrag]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["SCATTERING"]["PROCESS"] = "yy[EPA]<F> -> mu+ mu-";
  card["SCATTERING"]["BEAMFRAG"] = "fewbody";
  card["SCATTERING"]["LOOPSCREEN"] = false;
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto &process = *generator.proc;
  auto param = std::make_shared<gra::MNstarParam>(*process.state.nstar_param);
  process.state.nstar_param = param;
  param->fewbody.body_br = {0.0, 1.0};
  param->fewbody.unweight = false;
  CHECK_THROWS_WITH(process.PrepareRun(), Catch::Contains("FEWBODY.unweight"));
  param->fewbody.unweight = true;
  CHECK_NOTHROW(process.PrepareRun());
  param->fewbody.body_br = {1.0, 0.0};
  param->fewbody.unweight = false;
  CHECK_NOTHROW(process.PrepareRun());
}
