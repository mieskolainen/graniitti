// Core model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <numeric>

#include "Graniitti/MGraniitti.h"
#include "Graniitti/MUserHistograms.h"
#include "Graniitti/MUserCuts.h"
#include "support/models_test_support.hh"

TEST_CASE("LORENTZSCALAR shared caches copy without owning-pointer duplication", "[gra::LORENTZSCALAR][threading]") {
  static_assert(std::is_same_v<decltype(std::declval<gra::LORENTZSCALAR>().GlobalSudakovPtr),
                               std::shared_ptr<const gra::MSudakov>>);
  static_assert(
      std::is_same_v<decltype(std::declval<gra::LORENTZSCALAR>().GlobalPdfPtr), std::shared_ptr<const LHAPDF::PDF>>);

  gra::LORENTZSCALAR first;
  first.GlobalSudakovPtr = std::make_shared<const gra::MSudakov>();
  first.model_cache      = std::make_shared<gra::MModelCache>();
  first.GlobalPdfPtr     = first.model_cache->pdf.GetPDF("CT10nlo", 0);

  gra::LORENTZSCALAR second = first;
  REQUIRE(first.GlobalSudakovPtr == second.GlobalSudakovPtr);
  REQUIRE(first.GlobalPdfPtr == second.GlobalPdfPtr);
  REQUIRE(first.model_cache == second.model_cache);

  second.GlobalSudakovPtr = nullptr;
  second.GlobalPdfPtr     = nullptr;
  REQUIRE(first.GlobalSudakovPtr != nullptr);
  REQUIRE(first.GlobalPdfPtr != nullptr);
  REQUIRE(second.GlobalSudakovPtr == nullptr);
  REQUIRE(second.GlobalPdfPtr == nullptr);
}

TEST_CASE("Generated MadGraph registry factories own independent amplitude state", "[MG2GRA][threading]") {
  static_assert(!std::is_copy_constructible_v<MG5_PP_Z::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_PP_ZJ::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_PP_JJ::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_PP_W::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_YY_JJ::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_YY_WW::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<MG5_YY_ZJJ::ProcessBase>);
  static_assert(!std::is_copy_constructible_v<gra::AMP_MG5_yy_ww>);
  static_assert(!std::is_copy_assignable_v<gra::AMP_MG5_yy_ww>);
  static_assert(!std::is_copy_constructible_v<gra::DurhamMG5Process>);
  static_assert(!std::is_copy_constructible_v<gra::PhotonMG5Process>);

  for (const auto &process : gra::amplitude::Processes("DURHAM")) {
    INFO(process.process_name);
    auto first  = gra::CreateDurhamMG5Process(process);
    auto second = gra::CreateDurhamMG5Process(process);
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    REQUIRE(first.get() != second.get());
    REQUIRE(first->Name() == process.process_name);
    REQUIRE(first->FinalPDGs() == second->FinalPDGs());
    REQUIRE(first->ColorCount() == second->ColorCount());
    REQUIRE(first->HelicityCount() == second->HelicityCount());
  }

  for (const auto &process : gra::amplitude::Processes("PHOTON")) {
    INFO(process.process_name);
    auto first  = gra::CreatePhotonMG5Process(process);
    auto second = gra::CreatePhotonMG5Process(process);
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    REQUIRE(first.get() != second.get());
  }
}

TEST_CASE(
    "MProcess final-state symmetry factor follows the integrated "
    "physical state",
    "[MProcess][symmetry][physics]") {
  ToyHelicityProcess proc;
  const auto         tree = TensorRhoCascadeLTSForTest().decaytree;

  CHECK(proc.SymmetryFactorForTest(tree, false) == Approx(4.0));
  CHECK(proc.SymmetryFactorForTest(tree, true) == Approx(4.0));
  CHECK(proc.SymmetryFactorForTest(tree, false, true) == Approx(1.0));
  CHECK(proc.SymmetryFactorForTest(tree, true, true) == Approx(1.0));
}

TEST_CASE("Fast angular histograms use one daughter for each angle pair", "[MUserHistograms][physics]") {
  gra::LORENTZSCALAR event;
  event.pbeam1 = gra::M4Vec(0.0, 0.0, 10.0, 10.1);
  event.pbeam2 = gra::M4Vec(0.0, 0.0, -10.0, 10.1);
  event.decaytree.resize(2);
  event.decaytree[0].p4 = gra::M4Vec(1.0, 0.0, 0.0, std::sqrt(2.0));
  event.decaytree[1].p4 = gra::M4Vec(-1.0, 0.0, 0.0, std::sqrt(2.0));
  event.pfinal          = {event.decaytree[0].p4 + event.decaytree[1].p4};
  event.q1              = gra::M4Vec(0.0, 0.0, 1.0, 0.0);
  event.q2              = gra::M4Vec(0.0, 0.0, -1.0, 0.0);

  gra::MUserHistograms histograms;
  histograms.InitHistograms();
  histograms.FillCosThetaPhi(1.0, event);

  int phi_bin = -1;
  histograms.h1.at("phi_CM").GetBinIdx(0.0, phi_bin);
  CHECK(histograms.h1.at("phi_CM").GetBinCount(phi_bin) == 1);

  int costheta_bin  = -1;
  int joint_phi_bin = -1;
  histograms.h2.at("costhetaphi_CM").GetBinIdx(0.0, 0.0, costheta_bin, joint_phi_bin);
  CHECK(histograms.h2.at("costhetaphi_CM").GetBinCount(costheta_bin, joint_phi_bin) == 1);
}

TEST_CASE("Fast histograms skip absent process observables safely", "[MUserHistograms][safety]") {
  gra::LORENTZSCALAR   event;
  gra::MUserHistograms histograms;
  histograms.SetHistograms(2);
  histograms.InitHistograms();

  CHECK_NOTHROW(histograms.FillHistograms(1.0, event));
  CHECK(histograms.h2.at("rap1rap2").FillCount() == 0);
  CHECK(histograms.h1.at("phi_CM").FillCount() == 0);
}

TEST_CASE("MKinematics uses the exact incoming-state Moller flux", "[gra::kinematics]") {
  ToyHelicityProcess proc;
  const double       p  = 3.0;
  const double       m1 = 2.0;
  const double       m2 = 0.7;
  proc.state.lts.pbeam1 = gra::M4Vec(0.0, 0.0, p, std::sqrt(pow2(p) + pow2(m1)));
  proc.state.lts.pbeam2 = gra::M4Vec(0.0, 0.0, -p, std::sqrt(pow2(p) + pow2(m2)));
  proc.state.lts.s      = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();

  const double lambda = pow2(proc.state.lts.s) + pow2(pow2(m1)) + pow2(pow2(m2)) -
                        2.0 * proc.state.lts.s * (pow2(m1) + pow2(m2)) - 2.0 * pow2(m1) * pow2(m2);
  REQUIRE(kinematics::MollerFlux(proc.state.lts.pbeam1, proc.state.lts.pbeam2) ==
          Approx(2.0 * std::sqrt(lambda)).epsilon(1e-13));
  REQUIRE(kinematics::MollerFlux(proc.state.lts.pbeam1, proc.state.lts.pbeam2) !=
          Approx(2.0 * proc.state.lts.s).epsilon(1e-6));

  proc.state.lts.pbeam1 = gra::M4Vec(0.0, 0.0, p, p);
  proc.state.lts.pbeam2 = gra::M4Vec(0.0, 0.0, -p, p);
  proc.state.lts.s      = (proc.state.lts.pbeam1 + proc.state.lts.pbeam2).M2();
  REQUIRE(kinematics::MollerFlux(proc.state.lts.pbeam1, proc.state.lts.pbeam2) ==
          Approx(2.0 * proc.state.lts.s).epsilon(1e-13));
}

// Check JZ steering preserves density populations without imposing coherent phases
TEST_CASE("Resonance JZ steering distinguishes coherent amplitudes and densities", "[gra::resonance][regression]") {
  const auto source = std::filesystem::path(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")).parent_path();
  const auto dir = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/test_jz_density"));
  std::filesystem::create_directories(dir);
  std::filesystem::copy(source, dir, std::filesystem::copy_options::recursive |
                                     std::filesystem::copy_options::overwrite_existing);
  const auto original = nlohmann::json::parse(gra::aux::GetInputData((source / "RES/eta2_1645.json").string()));
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["MODELPARAM"] = dir.string();
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["RES"] = {"eta2_1645"};
  for (const std::string mode : {"none", "rho", "a_Jz"}) {
    CAPTURE(mode);
    auto resonance = original;
    auto &pol = resonance["PARAM_RES"]["MODELS"]["MP"]["[995,995]"]["polarization"];
    pol["mode"] = mode;
    pol["a_Jz"] = {{-2, 0.0, 0.0}, {-1, std::sqrt(0.5), 0.0}, {0, 0.0, 0.0}};
    { std::ofstream output(dir / "RES/eta2_1645.json"); output << resonance.dump(); }
    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{JZ0:1}";
    gra::MGraniitti generator;
    if (mode == "none") {
      REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
    } else if (mode == "a_Jz") {
      // Intrinsic parity belongs to the production vertex and does not forbid a pure Jz=0 spin state
      REQUIRE_NOTHROW(generator.ReadInput(card));
      const auto &zero = generator.proc->state.lts.process.RESONANCES.at("eta2_1645");
      REQUIRE(zero.UsesCoherentSpinBasis());
      CHECK(std::norm(zero.a_Jz[2]) == Approx(1.0));
      CHECK(gra::SquaredNorm(zero.a_Jz) == Approx(1.0));
      card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{JZ1:1}";
      REQUIRE_NOTHROW(generator.ReadInput(card));
      const auto &res = generator.proc->state.lts.process.RESONANCES.at("eta2_1645");
      REQUIRE(res.UsesCoherentSpinBasis());
      CHECK(std::norm(res.a_Jz[2]) < 1e-12);
      CHECK(std::norm(res.a_Jz[1]) == Approx(0.5));
      CHECK(std::abs(res.a_Jz[3] + res.a_Jz[1]) < 1e-12);
    } else {
      REQUIRE_NOTHROW(generator.ReadInput(card));
      const auto &res = generator.proc->state.lts.process.RESONANCES.at("eta2_1645");
      REQUIRE(res.UsesDensitySpinBasis());
      CHECK(res.rho[2][2].real() == Approx(1.0));
      CHECK(res.rho.Trace().real() == Approx(1.0));
    }
    for (const std::string value : {"nan", "inf", "-1", "1junk", "1e999", "1e308"}) {
      CAPTURE(value);
      card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{JZ1:" + value + "}";
      REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
    }
  }

  // Mass, width and tensor overrides use the same strict numeric input rules
  gra::MGraniitti generator;
  for (const std::string key : {"M", "W", "g0"}) {
    for (const std::string value : {"nan", "inf", "1junk", "1e999"}) {
      CAPTURE(key, value);
      card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{" + key + ":" + value + "}";
      REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
    }
  }
  for (const std::string key : {"M", "W"}) {
    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{" + key + ":-1}";
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
  for (const std::string width : {"0.03", "0"}) {
    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> eta0 pi+ pi- @R[eta2_1645]{W:" + width + "}";
    REQUIRE_NOTHROW(generator.ReadInput(card));
    const auto &particle = generator.proc->state.lts.process.RESONANCES.at("eta2_1645").p;
    CHECK(particle.width == Approx(std::stod(width)));
    CHECK(particle.tau == Approx(particle.width > 0.0 ? gra::PDG::hbar / particle.width : 0.0));
  }
}

// Check short integrations preserve buffered weights and rejected trials in saved histograms
TEST_CASE("Fast histogram output includes unflushed automatic-range samples", "[MUserHistograms][regression]") {
  gra::MUserHistograms histograms;
  histograms.SetHistograms(1);
  histograms.InitHistograms();
  gra::LORENTZSCALAR event;
  event.m2 = 4.0;
  for (const double weight : {1.0, 0.0, 2.0}) { histograms.FillHistograms(weight, event); }
  const auto directory = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/test_histogram_buffer"));
  const auto filename = (directory / "nested/short.hfast").string();
  histograms.SaveHistograms(filename);
  const auto data = nlohmann::json::parse(gra::aux::GetInputData(filename)).at("h1").at("M");
  CHECK(data.at("fills").get<unsigned int>() == 3);
  const auto weights = data.at("weights").get<std::vector<double>>();
  const auto weights2 = data.at("weights2").get<std::vector<double>>();
  CHECK(std::accumulate(weights.begin(), weights.end(), 0.0) == Approx(3.0));
  CHECK(std::accumulate(weights2.begin(), weights2.end(), 0.0) == Approx(5.0));
}

// Check rejected trials contribute zero to angular Monte Carlo normalization
TEST_CASE("Fast histograms retain rejected trials without rest-frame boosts", "[MUserHistograms][regression]") {
  gra::LORENTZSCALAR event;
  event.decaytree.resize(2);
  event.pfinal.resize(11);
  gra::MUserHistograms histograms;
  histograms.SetHistograms(2);
  histograms.InitHistograms();
  REQUIRE_NOTHROW(histograms.FillHistograms(0.0, event));
  for (const auto &frame : {"CM", "CS", "HX", "AH", "GJ", "LA", "PG"}) {
    const std::string name(frame);
    CHECK(histograms.h1.at("phi_" + name).FillCount() == 1);
    CHECK(histograms.h2.at("costhetaphi_" + name).FillCount() == 1);
  }

  event.pbeam1 = gra::M4Vec(0.0, 0.0, 10.0, 10.1);
  event.pbeam2 = gra::M4Vec(0.0, 0.0, -10.0, 10.1);
  event.decaytree[0].p4 = gra::M4Vec(1.0, 0.0, 0.0, std::sqrt(2.0));
  event.decaytree[1].p4 = gra::M4Vec(-1.0, 0.0, 0.0, std::sqrt(2.0));
  event.pfinal[0] = event.decaytree[0].p4 + event.decaytree[1].p4;
  event.q1 = gra::M4Vec(0.0, 0.0, 1.0, 0.0);
  event.q2 = gra::M4Vec(0.0, 0.0, -1.0, 0.0);
  REQUIRE_NOTHROW(histograms.FillCosThetaPhi(2.0, event));
  CHECK(histograms.h1.at("phi_CM").FillCount() == 2);
  CHECK(histograms.h1.at("phi_CM").WeightMeanAndError().first == Approx(1.0));
  CHECK(histograms.h2.at("costhetaphi_CM").WeightMeanAndError().first == Approx(1.0));
}

// Check the real process rejection path with angular histograms enabled
TEST_CASE("Below-threshold sampling remains a counted kinematic failure", "[MProcess][MUserHistograms]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 2;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["PROCESS"] = "MP[CON]<F> -> pi+ pi-";
  gra::MGraniitti generator;
  generator.ReadInput(card);
  generator.proc->PrepareRun();
  auto cuts = generator.proc->state.gcuts;
  cuts.M_min = 0.001;
  cuts.M_max = 0.01;
  generator.proc->SetGenCuts(cuts);
  gra::MEventWeightState state;
  double weight = 1.0;
  REQUIRE_NOTHROW(weight = generator.proc->EventWeight(std::vector<double>(generator.proc->GetdLIPSDim(), 0.5), state));
  CHECK(weight == Approx(0.0));
  CHECK_FALSE(state.kinematics_ok);
}

// Check interleaved steering and worker creation retain the selected model cards
TEST_CASE("Generators isolate particle and numerical tune selection", "[MGraniitti][MModelTune][threading]") {
  const auto source = std::filesystem::path(gra::ResolveModelTuneDir("TUNE0"));
  const auto target = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/test_core_tune"));
  std::filesystem::create_directories(target);
  std::filesystem::copy(source, target, std::filesystem::copy_options::recursive |
                                        std::filesystem::copy_options::overwrite_existing);
  auto particles = nlohmann::json::parse(gra::aux::GetInputData((source / "PDG_EXTRA.json").string()));
  const double mass = particles.at("PARAM_PDG").at("monopolium").at("mass");
  particles["PARAM_PDG"]["monopolium"]["mass"] = mass * 1.1;
  { std::ofstream out(target / "PDG_EXTRA.json"); out << particles.dump(); }

  auto resonance = nlohmann::json::parse(gra::aux::GetInputData((source / "RES/f0_980.json").string()));
  const double pole_mass = resonance.at("PARAM_RES").at("MODELS").at("MP").at("mass");
  resonance["PARAM_RES"]["MODELS"]["MP"]["mass"] = pole_mass * 1.05;
  { std::ofstream out(target / "RES/f0_980.json"); out << resonance.dump(); }
  auto decays = nlohmann::json::parse(gra::aux::GetInputData((source / "DECAYS.json").string()));
  const double diphoton_br = decays.at("25").at("[22,22]").at("BR");
  decays["25"]["[22,22]"]["BR"] = diphoton_br * 0.5;
  { std::ofstream out(target / "DECAYS.json"); out << decays.dump(); }

  auto first_card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  first_card["GENERIC"]["CORES"] = 1;
  first_card["GENERIC"]["HIST"] = 0;
  first_card["GENERIC"]["MODELPARAM"] = "TUNE0";
  first_card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> pi+ pi-";
  first_card["SCATTERING"]["RES"] = {"f0_980"};
  auto second_card = first_card;
  second_card["GENERIC"]["MODELPARAM"] = target.string();
  const std::string standalone_tune = gra::MODELPARAM;
  gra::MGraniitti first;
  first.ReadInput(first_card);
  first.proc->PrepareRun();
  gra::MGraniitti second;
  REQUIRE_NOTHROW(first.InitMultiMemory());
  const auto first_tune = first.proc->GetModelTune();
  const auto first_seed = first.proc->state.random.GetSeed();
  second.ReadInput(second_card);
  first.SetCores(2);
  REQUIRE_NOTHROW(first.InitMultiMemory());
  CHECK(first.proc->state.random.GetSeed() == first_seed);
  REQUIRE_THROWS_AS(first.SetCores(-1), std::invalid_argument);
  CHECK(first.GetCores() == 2);
  CHECK(first.proc->GetModelTune() == first_tune);
  first.SetCores(1);
  REQUIRE_NOTHROW(first.InitMultiMemory());
  CHECK(first.proc->GetModelTune() == first_tune);
  CHECK(first.proc->state.lts.PDG.FindByPDG(881).mass == Approx(mass));
  CHECK(second.proc->state.lts.PDG.FindByPDG(881).mass == Approx(mass * 1.1));
  CHECK(first.proc->GetModelTune()->Directory() != second.proc->GetModelTune()->Directory());
  CHECK(gra::MODELPARAM == standalone_tune);
  const auto second_tune = second.proc->GetModelTune();
  { std::ofstream out(target / "NUMERICS.json"); out << "invalid after loading"; }
  REQUIRE_NOTHROW(second.InitMultiMemory());
  CHECK(second.proc->GetModelTune() == second_tune);
  CHECK(first.proc->state.lts.process.RESONANCES.at("f0_980").p.mass == Approx(pole_mass));
  CHECK(second.proc->state.lts.process.RESONANCES.at("f0_980").p.mass == Approx(pole_mass * 1.05));
  const auto higgs = first.proc->state.lts.PDG.FindByPDG(25);
  CHECK(gra::resonance::GammaGammaPartialWidth(higgs, first_tune->Directory()) == Approx(higgs.width * diphoton_br));
  CHECK(gra::resonance::GammaGammaPartialWidth(higgs, second.proc->GetModelTune()->Directory()) ==
        Approx(higgs.width * diphoton_br * 0.5));

  first.ReadGeneralParam(first_card);
  second.ReadGeneralParam(second_card);
  REQUIRE_NOTHROW(first.ReadProcessParam(first_card));
  CHECK(first.proc->state.lts.PDG.FindByPDG(881).mass == Approx(mass));
}

// Check non-finite central momenta cannot pass the STAR upper-pT conditions
TEST_CASE("User cuts reject non-finite central daughter momenta", "[MUserCuts][regression]") {
  gra::LORENTZSCALAR event;
  event.pfinal[1] = gra::M4Vec(0.0, 0.3, 10.0, 10.1);
  event.pfinal[2] = gra::M4Vec(0.0, -0.3, -10.0, 10.1);
  event.decaytree.resize(2);
  event.decaytree[0].p4 = gra::M4Vec(0.4, 0.0, 0.0, 1.0);
  event.decaytree[1].p4 = gra::M4Vec(-0.4, 0.0, 0.0, 1.0);
  for (const auto cut : {1792394010, 1792394020}) {
    REQUIRE(gra::UserCut(cut, event));
    for (const auto i : indices(event.decaytree)) {
      auto central = event.decaytree;
      central[i].p4 = gra::M4Vec(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 1.0);
      CHECK_FALSE(gra::UserCut(cut, event, central));
      central[i].p4 = gra::M4Vec(0.4, 0.0, 0.0, std::numeric_limits<double>::infinity());
      CHECK_FALSE(gra::UserCut(cut, event, central));
    }
  }
}

// Check malformed numerical process steering fails before any sampling
TEST_CASE("Process numerical overrides reject suffixes and invalid resonance widths", "[MGraniitti][input]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["SCATTERING"]["RES"] = {"f0_980"};
  gra::MGraniitti generator;
  for (const auto command : {"@PDG[211oops]{M:0.14}", "@PDG[211]{M:0.14GeV}",
                              "@FLATAMP:1oops", "@OFFSHELL:5oops", "@MMAX:2.5",
                              "@R[f0_980]{M:-1}", "@R[f0_980]{W:-0.1}", "@R[f0_980]{JZ0:nan}"}) {
    card["SCATTERING"]["PROCESS"] = std::string("MP[RES]<F> -> pi+ pi- ") + command;
    CHECK_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
  card["SCATTERING"]["PROCESS"] = "MP[RES]<F> -> pi+ pi- @R[f0_980]{W:0.1}";
  REQUIRE_NOTHROW(generator.ReadInput(card));
  const auto &particle = generator.proc->GetResonances().at("f0_980").p;
  CHECK(particle.width == Approx(0.1));
  CHECK(particle.tau == Approx(gra::PDG::hbar / particle.width));
}

// Check that process steering cannot contain ignored controls
TEST_CASE("Process steering rejects inapplicable controls", "[MGraniitti][input]") {
  const auto source  = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  const auto rejects = [&](nlohmann::json card) {
    gra::MGraniitti generator;
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  };

  auto card                     = source;
  card["SCATTERING"]["PROCESS"] = "GP[CON]<F> -> pi+ pi- @UNKNOWN:1";
  rejects(card);

  card                          = source;
  card["SCATTERING"]["PROCESS"] = "GP[CON]<F> -> pi+ pi- @MP_FRAME:CS";
  rejects(card);

  card                          = source;
  card["SCATTERING"]["PROCESS"] = "MP[CON]<F> -> pi+ pi- @MMAX:1";
  rejects(card);

  for (const std::string process : {"yy[EPA]<F> -> mu+ mu-", "yy[QED]<F> -> mu+ mu-",
                                    "yy_DZ[EPA]<P> -> mu+ mu-", "yy_LUX[EPA]<P> -> mu+ mu-",
                                    "ygg[jpsi]<F> -> mu+ mu-", "gg[chic(0)]<F> &> pi+ pi-",
                                    "MP[CON]<F> -> pi+ pi-", "XP[CON]<F> -> pi+ pi-",
                                    "GP[CON]<F> -> pi+ pi-", "TP[CON]<F> -> pi+ pi-"}) {
    for (const std::string command : {"@RES{f0_980:1}", "@R[f0_980]{W:0.1}"}) {
      CAPTURE(process, command);
      card = source;
      card["SCATTERING"]["PROCESS"] = process + " " + command;
      gra::MGraniitti generator;
      REQUIRE_THROWS_WITH(generator.ReadInput(card), Catch::Contains("requires a process with free resonance input"));
    }
  }

  card = source;
  card["SCATTERING"]["RES"] = nlohmann::json::array();
  card["SCATTERING"]["PROCESS"] = "TP[PHOTO]<F> -> pi+ pi- @R[rho_770]{W:0.1}";
  gra::MGraniitti generator;
  REQUIRE_THROWS_WITH(generator.ReadInput(card), Catch::Contains("invalid R[rho_770] not found"));
  card["SCATTERING"]["RES"] = {"rho_770"};
  REQUIRE_NOTHROW(generator.ReadInput(card));

  card                    = source;
  card["GENERIC"]["HIST"] = 3;
  rejects(card);
}

// Check charge conservation and strict physical decay input before sampling
TEST_CASE("Decay input rejects forbidden and unsupported cascades", "[MGraniitti][input][decay]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["CORES"] = 1;
  card["SCATTERING"]["LOOPSCREEN"] = false;
  for (const std::string arrow : {"->", "&>"}) {
    gra::MGraniitti generator;
    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> " + arrow + " pi0 22 @RES{rho_770:1}";
    REQUIRE_NOTHROW(generator.ReadInput(card));
    if (arrow == "->") {
      REQUIRE_THROWS_AS(generator.proc->InitializeProcessAmplitude(), gra::MissingHelicityData);
    } else {
      REQUIRE_NOTHROW(generator.proc->InitializeProcessAmplitude());
      const auto &decay = generator.proc->state.lts.process.RESONANCES.at("rho_770").hel_decay;
      // An isolated decay without data has unit normalization and no inferred LS couplings
      REQUIRE(decay.ls_components.empty());
      REQUIRE(std::abs(decay.g_decay - std::complex<double>(1.0, 0.0)) < 1e-12);
    }

    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> " + arrow +
        " rho(770)0 > {pi0 pi+} rho(770)0 > {pi+ pi-} @RES{f2_2150:1}";
    REQUIRE_THROWS_WITH(generator.ReadInput(card), Catch::Contains("violates charge conservation"));

    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> " + arrow +
        " rho(770)0 > {pi0 22} rho(770)0 > {pi+ pi-} @RES{f2_2150:1}";
    REQUIRE_NOTHROW(generator.ReadInput(card));
    // Isolating the central decay does not supply missing nested decay couplings
    REQUIRE_THROWS_AS(generator.proc->InitializeProcessAmplitude(), gra::MissingHelicityData);

    card["SCATTERING"]["PROCESS"] = "MP[RES]<F> " + arrow +
        " rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-} @RES{f2_2150:1}";
    REQUIRE_NOTHROW(generator.ReadInput(card));
    REQUIRE_NOTHROW(generator.proc->InitializeProcessAmplitude());
  }
}

// Check continuum ignores resonance cards while resonance processes require active states
TEST_CASE("Continuum steering ignores the configured resonance list", "[MGraniitti][input]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["SCATTERING"]["LOOPSCREEN"] = false;
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    CAPTURE(model);
    card["SCATTERING"]["PROCESS"] = model + "[CON]<F> -> pi+ pi-";
    for (const auto &names : {std::vector<std::string>{}, std::vector<std::string>{"unused_resonance"}}) {
      card["SCATTERING"]["RES"] = names;
      gra::MGraniitti generator;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      CHECK(generator.proc->GetResonances().empty());
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
    }
    for (const std::string channel : {"RES", "RES+CON"}) {
      CAPTURE(channel);
      card["SCATTERING"]["PROCESS"] = model + "[" + channel + "]<F> -> pi+ pi-";
      card["SCATTERING"]["RES"] = nlohmann::json::array();
      gra::MGraniitti generator;
      REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
    }
  }
}

// Check that a flat amplitude cannot silently bypass requested screening
TEST_CASE("Flat amplitudes reject Pomeron screening", "[MGraniitti][input]") {
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["SCATTERING"]["PROCESS"] = "GP[CON]<F> -> pi+ pi- @FLATAMP:1";
  card["SCATTERING"]["LOOPSCREEN"] = true;
  gra::MGraniitti generator;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
}

// Check every STAR 510 GeV Roman Pot branch and central pair pT condition
TEST_CASE("STAR 510 GeV fiducial cuts use the physical proton arm", "[MUserCuts][STAR510]") {
  gra::LORENTZSCALAR event;
  event.decaytree.resize(2);
  event.decaytree[0].p4 = gra::M4Vec(0.5, 0.0, 0.0, 1.2);
  event.decaytree[1].p4 = gra::M4Vec(-0.5, 0.0, 0.0, 1.2);
  // [REFERENCE: STAR Collaboration, arXiv:2510.27482, Eq. (4.1) and Table 1]
  const std::vector<std::pair<gra::M4Vec, bool>> points = {
      {{0.0, 0.6, -250.0, 251.0}, true},   {{0.0, 0.87, -250.0, 251.0}, false},
      {{0.0, -0.9, -250.0, 251.0}, true},  {{-0.24, -0.9, -250.0, 251.0}, false},
      {{0.0, 0.9, 250.0, 251.0}, true},   {{-0.20, 0.9, 250.0, 251.0}, false},
      {{0.0, -0.87, 250.0, 251.0}, true}, {{0.0, -0.89, 250.0, 251.0}, false},
      {{0.0, 0.44, -250.0, 251.0}, true}, {{0.0, 0.44, 250.0, 251.0}, false},
      {{0.5, 0.6, 250.0, 251.0}, false}, {{0.0, 0.0, 250.0, 251.0}, false}};
  for (const auto &[proton, passed] : points) {
    event.pfinal[1] = proton;
    event.pfinal[2] = gra::M4Vec(0.0, -0.6, -proton.Pz(), 251.0);
    for (const auto id : {3075716000LL, 3075716010LL, 3075716020LL}) {
      CHECK(gra::IsKnownUserCut(id));
      CHECK(gra::UserCut(id, event) == passed);
      std::swap(event.pfinal[1], event.pfinal[2]);
      CHECK(gra::UserCut(id, event) == passed);
      std::swap(event.pfinal[1], event.pfinal[2]);
    }
  }
  const std::array<std::array<double, 2>, 4> bounds = {{
      {-0.2300, 0.4200}, {-0.2500, 0.4800}, {-0.2100, 0.4600}, {-0.1900, 0.4600}}};
  for (const auto i : indices(bounds)) {
    const double z = i < 2 ? -250.0 : 250.0;
    const double sign = i % 2 == 0 ? 1.0 : -1.0;
    const auto &rp = bounds[i];
    event.pfinal[2] = gra::M4Vec(0.0, -0.6, -z, 251.0);
    event.pfinal[1] = gra::M4Vec(rp[0], sign * 0.6, z, 251.0);
    CHECK_FALSE(gra::UserCut(3075716000LL, event));
    event.pfinal[1] = gra::M4Vec(rp[0] + 1e-7, sign * 0.6, z, 251.0);
    CHECK(gra::UserCut(3075716000LL, event));
    event.pfinal[1] = gra::M4Vec(0.0, sign * rp[1], z, 251.0);
    CHECK_FALSE(gra::UserCut(3075716000LL, event));
    event.pfinal[1] = gra::M4Vec(0.0, sign * (rp[1] + 1e-7), z, 251.0);
    CHECK(gra::UserCut(3075716000LL, event));
  }
  // Probe both sides of each circular boundary with independent points on the published arcs
  const std::array<std::array<double, 5>, 6> arcs = {{
      {-250.0,  1.0,  0.40, 1.3600,  1.04}, {-250.0, -1.0,  0.35, 1.5000,  1.05},
      {-250.0, -1.0, -0.20, 0.9590, -0.45}, { 250.0,  1.0,  0.35, 1.3000,  0.95},
      { 250.0,  1.0, -0.15, 0.9460, -0.43}, { 250.0, -1.0,  0.35, 1.5000,  1.05}}};
  for (const auto &[z, sign, px, r2, dx] : arcs) {
    event.pfinal[2] = gra::M4Vec(0.0, -0.6, -z, 251.0);
    for (const double shift : {-1e-7, 1e-7}) {
      event.pfinal[1] = gra::M4Vec(px, sign * (std::sqrt(r2 - dx * dx) + shift), z, 251.0);
      CHECK(gra::UserCut(3075716000LL, event) == (shift < 0.0));
    }
  }
  for (const double z : {-250.0, 250.0}) {
    const double sign = z < 0.0 ? 1.0 : -1.0, py = z < 0.0 ? 0.8600 : 0.8800;
    event.pfinal[2] = gra::M4Vec(0.0, -0.6, -z, 251.0);
    event.pfinal[1] = gra::M4Vec(0.0, sign * py, z, 251.0);
    CHECK_FALSE(gra::UserCut(3075716000LL, event));
    event.pfinal[1] = gra::M4Vec(0.0, sign * (py - 1e-7), z, 251.0);
    CHECK(gra::UserCut(3075716000LL, event));
  }
  event.pfinal[1] = gra::M4Vec(0.0, 0.6, 250.0, 251.0);
  event.pfinal[2] = gra::M4Vec(0.0, -0.6, -250.0, 251.0);
  const std::array<double, 2> pt_max = {0.7000, 1.1000};
  for (const auto i : indices(pt_max)) {
    const auto id = i == 0 ? 3075716010LL : 3075716020LL;
    const double limit = pt_max[i];
    auto central = event.decaytree;
    central[0].p4 = gra::M4Vec(limit, 0.0, 0.0, 2.0);
    central[1].p4 = gra::M4Vec(-limit, 0.0, 0.0, 2.0);
    CHECK_FALSE(gra::UserCut(id, event, central));
    CHECK(gra::UserCut(3075716000LL, event, central));
    central[1].p4 = gra::M4Vec(-limit + 1e-7, 0.0, 0.0, 2.0);
    CHECK(gra::UserCut(id, event, central));
    std::swap(central[0], central[1]);
    CHECK(gra::UserCut(id, event, central));
    central[0].p4 = gra::M4Vec(std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0, 2.0);
    CHECK_FALSE(gra::UserCut(id, event, central));
  }
}


// Select independent pole parameters before decay normalization and line shape evaluation
TEST_CASE("Generator resonance poles follow the production model", "[MGraniitti][resonance][pole]") {
  const auto source = std::filesystem::path(gra::ResolveModelTuneDir("TUNE0"));
  const auto target = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/test_res_poles"));
  std::filesystem::create_directories(target);
  std::filesystem::copy(source, target, std::filesystem::copy_options::recursive |
                                        std::filesystem::copy_options::overwrite_existing);
  auto resonance = nlohmann::json::parse(gra::aux::GetInputData((source / "RES/f0_500.json").string()));
  auto &models = resonance["PARAM_RES"]["MODELS"];
  const std::array<std::string, 4> names = {"GP", "MP", "XP", "TP"};
  const std::array<std::string, 3> modes = {"fixed-width", "kinematic-width", "running-width"};
  const std::array<gra::BreitWigner, 4> bw = {gra::BreitWigner::FixedWidth, gra::BreitWigner::KinematicWidth,
                                           gra::BreitWigner::RunningWidth, gra::BreitWigner::FixedWidth};
  for (const auto &i : indices(names)) {
    auto &pole = models[names[i]];
    pole["mass"] = pole.at("mass").get<double>() * (1.0 + 0.1 * i);
    pole["width"] = pole.at("width").get<double>() * (1.0 + 0.2 * i);
    if (i < modes.size()) { pole["BW"] = modes[i]; }
  }
  { std::ofstream output(target / "RES/f0_500.json"); output << resonance; }
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["MODELPARAM"] = target.string();
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["LOOPSCREEN"] = false;
  for (const auto &i : indices(names)) {
    for (const bool override : {false, true}) {
      CAPTURE(names[i], override);
      const double mass = models[names[i]].at("mass").get<double>() * (override ? 1.05 : 1.0);
      const double width = models[names[i]].at("width").get<double>() * (override ? 1.1 : 1.0);
      card["SCATTERING"]["RES"] = {override ? "f0_980" : "f0_500"};
      card["SCATTERING"]["PROCESS"] = names[i] + "[RES]<F> -> pi+ pi-" +
          (override ? " @RES{f0_500:1,f0_980:0} @R[f0_500]{M:" + nlohmann::json(mass).dump() +
                      ",W:" + nlohmann::json(width).dump() + "}" : "");
      gra::MGraniitti generator;
      generator.ReadInput(card);
      generator.proc->PrepareRun();
      REQUIRE(generator.proc->GetResonances().size() == 1);
      const auto &res = generator.proc->GetResonances().at("f0_500");
      CHECK(res.p.mass == Approx(mass));
      CHECK(res.p.width == Approx(width));
      CHECK(res.p.tau == Approx(gra::PDG::hbar / width));
      CHECK(res.BW == bw[i]);
      if (names[i] == "TP") { continue; }
      const double mass2 = gra::math::pow2(1.2 * mass);
      const double imaginary = width * (i == 1 ? std::sqrt(mass2) : mass * (i == 2 ? 1.3 : 1.0));
      RequireComplexNear(gra::resonance::LineShape(mass2, res, 1.3),
                         1.0 / std::complex<double>(mass2 - mass * mass, imaginary), 1e-12);
    }
  }
}
