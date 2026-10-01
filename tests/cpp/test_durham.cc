// Durham model tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"
#include "Graniitti/Amplitude/Photon/AMP_gg_yy.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "support/models_test_support.hh"

namespace {

using gra::aux::indices;
using DurhamFrameTransform = std::function<void(gra::M4Vec &)>;

// Write a complete local tune for Durham covariance and rejection regressions
gra::MModelTunePtr DurhamTestTune(const std::string                           &name,
                                  const std::function<void(nlohmann::json &)> &general_edit  = {},
                                  const std::function<void(nlohmann::json &)> &numerics_edit = {}) {
  const auto dir = std::filesystem::path("tmp") / ("test_durham_" + name);
  std::filesystem::create_directories(dir);
  for (const std::string file :
       {"GENERAL.json", "NUMERICS.json", "CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}) {
    auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", file)));
    if (file == "GENERAL.json" && general_edit) { general_edit(card); }
    if (file == "NUMERICS.json" && numerics_edit) { numerics_edit(card); }
    std::ofstream out(dir / file);
    out << card.dump();
    if (!out.good()) { throw std::runtime_error("DurhamTestTune: failed to write " + file); }
  }
  return gra::MModelTune::Load((dir / "GENERAL.json").string());
}

// Transform one Durham decay branch and all stable daughters
void TransformDurhamBranch(gra::MDecayBranch &branch, const DurhamFrameTransform &transform) {
  transform(branch.p4);
  for (auto &daughter : branch.legs) { TransformDurhamBranch(daughter, transform); }
}

// Transform all physical momenta in one Durham hard point
gra::LORENTZSCALAR TransformDurhamPoint(gra::LORENTZSCALAR lts, const DurhamFrameTransform &transform) {
  transform(lts.pbeam1);
  transform(lts.pbeam2);
  transform(lts.q1);
  transform(lts.q2);
  for (auto &momentum : lts.pfinal) { transform(momentum); }
  for (auto &branch : lts.decaytree) { TransformDurhamBranch(branch, transform); }
  lts.screening.durham.hard_cached = false;
  lts.screening.durham.jet_cuts_cached        = false;
  return lts;
}

// Read the bottom mass used by the exact generated Durham process
double DurhamBottomMass(const gra::amplitude::Process &process) {
  const std::string card = gra::aux::ResolveProjectPath("MG5cards/Durham/" + process.process_name + "/param_card.dat");
  const SLHAReader  slha(card);
  const double      mass = slha.get_block_entry("mass", 5, std::numeric_limits<double>::quiet_NaN());
  if (!std::isfinite(mass) || mass <= 0.0) {
    throw std::invalid_argument("DurhamBottomMass: invalid generated bottom mass");
  }
  return mass;
}

// Build an above-threshold on-shell bottom final state with exact closure
gra::LORENTZSCALAR DurhamBottomPoint(const gra::amplitude::Process &process, bool radiative) {
  gra::LORENTZSCALAR lts = radiative ? MakeToyDurhamQQbarG(5) : MakeToyDurhamQQbar(5);

  const double beam_pz = 30.0;
  const double beam_e  = std::sqrt(beam_pz * beam_pz + gra::PDG::mp * gra::PDG::mp);
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, beam_pz, beam_e);
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -beam_pz, beam_e);

  const auto forward = [](double px, double py, double pz) {
    const double energy = std::sqrt(px * px + py * py + pz * pz + gra::PDG::mp * gra::PDG::mp);
    return gra::M4Vec(px, py, pz, energy);
  };
  lts.pfinal[1]                 = forward(0.18, -0.05, 18.0);
  lts.pfinal[2]                 = forward(-0.11, 0.04, -17.0);
  const gra::M4Vec  hard        = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  const double      hard_mass   = hard.M();
  const double      bottom_mass = DurhamBottomMass(process);
  const std::size_t final_count = radiative ? 3 : 2;
  if (lts.decaytree.size() != final_count || hard_mass <= 2.0 * bottom_mass) {
    throw std::invalid_argument("DurhamBottomPoint: invalid above-threshold hard system");
  }
  lts.decaytree[0].p.mass = bottom_mass;
  lts.decaytree[1].p.mass = bottom_mass;

  if (!radiative) {
    const double              energy   = 0.5 * hard_mass;
    const double              momentum = std::sqrt(energy * energy - bottom_mass * bottom_mass);
    std::array<gra::M4Vec, 2> rest = {gra::M4Vec(momentum, 0.0, 0.0, energy), gra::M4Vec(-momentum, 0.0, 0.0, energy)};
    for (const auto &i : indices(rest)) {
      gra::kinematics::LorentzBoost(hard, hard_mass, rest[i], +1);
      lts.decaytree[i].p4 = rest[i];
    }
  } else {
    const double              radicand      = hard_mass * hard_mass - 3.0 * bottom_mass * bottom_mass;
    const double              momentum      = (-hard_mass + 2.0 * std::sqrt(radicand)) / 3.0;
    const double              bottom_energy = std::sqrt(momentum * momentum + bottom_mass * bottom_mass);
    const double              sin120        = std::sqrt(3.0) / 2.0;
    std::array<gra::M4Vec, 3> rest          = {gra::M4Vec(momentum, 0.0, 0.0, bottom_energy),
                                               gra::M4Vec(-0.5 * momentum, sin120 * momentum, 0.0, bottom_energy),
                                               gra::M4Vec(-0.5 * momentum, -sin120 * momentum, 0.0, momentum)};
    for (const auto &i : indices(rest)) {
      gra::kinematics::LorentzBoost(hard, hard_mass, rest[i], +1);
      lts.decaytree[i].p4 = rest[i];
    }
  }
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

// Build one physical toy point for every registered Durham final state
gra::LORENTZSCALAR DurhamRegistryPoint(const gra::amplitude::Process &process) {
  const auto &pdgs = process.stable_pdgs;
  if (pdgs == std::vector<int>{22, 22}) {
    gra::LORENTZSCALAR lts = MakeToyDurhamGG();
    for (auto &branch : lts.decaytree) {
      branch.p.name  = "a";
      branch.p.pdg   = 22;
      branch.p.color = 1;
    }
    lts.q1              = gra::M4Vec(0.0, 0.0, 100.0, 100.0);
    lts.q2              = gra::M4Vec(0.0, 0.0, -100.0, 100.0);
    lts.decaytree[0].p4 = gra::M4Vec(80.0, 0.0, 60.0, 100.0);
    lts.decaytree[1].p4 = gra::M4Vec(-80.0, 0.0, -60.0, 100.0);
    lts.pfinal[0]       = lts.decaytree[0].p4 + lts.decaytree[1].p4;
    return lts;
  }
  if (pdgs == std::vector<int>{21, 21}) { return MakeToyDurhamGG(); }
  if (pdgs == std::vector<int>{21, 21, 21}) { return MakeToyDurhamGGG(); }
  if (pdgs == std::vector<int>{21, 21, 21, 21}) { return MakeToyDurhamGGGG(); }
  if (pdgs.size() == 2 && pdgs[0] > 0 && pdgs[1] == -pdgs[0]) {
    if (pdgs[0] == 5) { return DurhamBottomPoint(process, false); }
    return MakeToyDurhamQQbar(pdgs[0]);
  }
  if (pdgs.size() == 3 && pdgs[0] > 0 && pdgs[1] == -pdgs[0] && pdgs[2] == 21) {
    if (pdgs[0] == 5) { return DurhamBottomPoint(process, true); }
    return MakeToyDurhamQQbarG(pdgs[0]);
  }
  throw std::invalid_argument("DurhamRegistryPoint: unsupported registered final state");
}

// Compute the exact color projected hard norm after transverse source
// contraction
double ContractedDurhamNorm(gra::MDurham &durham, gra::DurhamMG5Process &matrix_element, gra::LORENTZSCALAR &lts,
                            const gra::M4Vec &source1, const gra::M4Vec &source2) {
  gra::MDurham::DurhamProjectedAmp hard;
  durham.Dgg2Generated(lts, matrix_element, hard);
  REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
  REQUIRE_FALSE(hard.empty());
  REQUIRE(lts.screening.durham.hard_cached);

  double norm = 0.0;
  for (const auto &channel : hard) {
    const std::complex<double> amplitude =
        gra::mg5helas::ContractTransverseSources(channel, source1, source2, lts.screening.durham.hard_k1, lts.screening.durham.hard_k2);
    REQUIRE(std::isfinite(amplitude.real()));
    REQUIRE(std::isfinite(amplitude.imag()));
    norm += std::norm(amplitude);
  }
  return norm;
}

}  // namespace

// Check public scale and helicity methods from an independent translation unit
TEST_CASE("Durham public projectors and PDF scales are linkable", "[gra::MDurham][helicity][params]") {
  ModelParamRestoreGuard restore;
  for (const auto &[scheme, expected] : std::map<std::string, std::array<double, 2>>{
           {"MIN", {1.0, 4.0}}, {"MAX", {4.0, 9.0}}, {"IN", {1.0, 9.0}}, {"EX", {4.0, 4.0}}, {"AVG", {2.5, 6.5}}}) {
    const auto tune = DurhamTestTune("public_" + scheme, [&](auto &j) { j["PARAM_DURHAM"]["PDF_scale"] = scheme; });
    auto       lts  = MakeToyDurhamGG();
    gra::MRandom rng;
    MDurham      durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
    double       first  = 0.0;
    double       second = 0.0;
    durham.DScaleChoise(4.0, 1.0, 9.0, first, second);
    REQUIRE(first == Approx(expected[0]));
    REQUIRE(second == Approx(expected[1]));

    const MDurham::DurhamTransverseMomentum q1 = {0.3, -0.7};
    const MDurham::DurhamTransverseMomentum q2 = {-0.2, 0.5};
    std::vector<std::complex<double>>       projector;
    durham.DHelicity(q1, q2, projector);
    REQUIRE(projector.size() == 4);
    const MDurham::DurhamInitialHelicity hard = {std::complex<double>(0.4, 0.2), {0.7, -0.3}, {-0.2, 0.1}, {0.5, 0.8}};
    const auto                           expected_amp = gra::mg5helas::ContractTransverseSources(
        hard, gra::M4Vec(q1[0], q1[1], 0.0, 0.0), gra::M4Vec(q2[0], q2[1], 0.0, 0.0), gra::M4Vec(0.0, 0.0, 10.0, 10.0),
        gra::M4Vec(0.0, 0.0, -10.0, 10.0));
    RequireComplexNear(durham.DHelProj(hard, projector), expected_amp);
    RequireComplexNear(durham.DHelProj(std::vector<std::complex<double>>(hard.begin(), hard.end()), projector),
                       expected_amp);
    for (const std::size_t size : {0U, 3U, 5U}) {
      const std::vector<std::complex<double>> invalid(size);
      REQUIRE_THROWS_AS(durham.DHelProj(invalid, projector), gra::AmplitudeFailure);
      REQUIRE_THROWS_AS(durham.DHelProj(hard, invalid), gra::AmplitudeFailure);
      REQUIRE_THROWS_AS(durham.DHelProj(std::vector<std::complex<double>>(hard.begin(), hard.end()), invalid),
                        gra::AmplitudeFailure);
    }
  }
}

// Check complete Durham evaluations classify resolved jet cuts as zero weight
TEST_CASE("Durham jet rejections preserve successful amplitude status", "[gra::MDurham][jets][color]") {
  ModelParamRestoreGuard restore;
  const std::string algorithm = GENERATE("anti-kt", "kt", "CA");
  for (const std::string cut : {"pt", "rap", "merge"}) {
    const auto tune = DurhamTestTune(algorithm + "_reject_" + cut, [&](auto &j) {
      j["PARAM_DURHAM"]["JET_ALGO"]    = algorithm;
      j["PARAM_DURHAM"]["JET_R"]       = cut == "merge" ? 10.0 : 0.4;
      j["PARAM_DURHAM"]["JET_pt_min"]  = cut == "pt" ? 100.0 : 0.0;
      j["PARAM_DURHAM"]["JET_rap_max"] = cut == "rap" ? 0.01 : 100.0;
    });
    CAPTURE(algorithm, cut);
    for (auto lts : {MakeToyDurhamGG(), MakeToyDurhamGGG(), MakeToyDurhamGGGG(), MakeToyDurhamQQbar(2),
                     MakeToyDurhamQQbarG(2)}) {
      gra::MRandom rng;
      MDurham      durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
      for (const bool screening : {false, true}) {
        lts.screening.active = screening;
        // Seed stale helicity amplitudes to check Born and cached screening rejection
        lts.hamp.assign(1, 1.0);
        REQUIRE(std::fpclassify(durham.DurhamQCD(lts, "MG5")) == FP_ZERO);
        REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
        REQUIRE(lts.screening.durham.jet_cuts_cached);
        REQUIRE_FALSE(lts.screening.durham.jet_cuts_pass);
        REQUIRE_FALSE(lts.screening.durham.hard_cached);
        REQUIRE_FALSE(lts.hamp.empty());
        REQUIRE(std::fpclassify(gra::SquaredNorm(lts.hamp)) == FP_ZERO);
        REQUIRE(lts.hard_color_flows.empty());
      }
    }
    auto photons = MakeToyDurhamGG();
    for (auto &branch : photons.decaytree) {
      branch.p.pdg   = 22;
      branch.p.color = 1;
      branch.p.name  = "a";
    }
    gra::MRandom rng;
    MDurham      durham(photons, tune, rng, gra::MDurham::ContinuumProcesses());
    const double photon_amp2 = durham.DurhamQCD(photons, "MG5");
    CAPTURE(static_cast<int>(durham.EvaluationStatus()));
    REQUIRE(photon_amp2 > 0.0);
    REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
    REQUIRE(photons.hard_color_flows.empty());
  }
}

// Check the sole color projection follows every shifted physical amplitude
TEST_CASE("Durham single color flow follows the screening momenta", "[gra::MDurham][color][screening]") {
  ModelParamRestoreGuard restore;
  const auto   tune = DurhamTestTune("screen_flow", [](auto &j) { j["PARAM_DURHAM"]["JET_pt_min"] = 0.0; });
  auto         lts  = MakeToyDurhamGG();
  gra::MRandom rng;
  MDurham      durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
  const double born = durham.DurhamQCD(lts, "MG5");
  REQUIRE(born > 0.0);
  REQUIRE(lts.hard_color_flows.size() == 1);
  const auto original = lts;
  for (const double shift : {0.15, 0.3, -0.1}) {
    lts                               = original;
    lts.pfinal_orig                   = lts.pfinal;
    lts.forward_mass2                 = {lts.pfinal[1].M2(), lts.pfinal[2].M2()};
    lts.screening.active              = true;
    const std::array<double, 2> upper = {lts.pfinal[1].Px() + shift, lts.pfinal[1].Py()};
    const std::array<double, 2> lower = {lts.pfinal[2].Px() - shift, lts.pfinal[2].Py()};
    REQUIRE(gra::kinematics::RebuildScreeningKinematics(lts, upper, lower, true));
    UpdateToyDurhamDerivedKinematics(lts);
    const double shifted = durham.DurhamQCD(lts, "MG5");
    REQUIRE(shifted > 0.0);
    REQUIRE(std::abs(shifted - born) > 1e-6 * born);
    REQUIRE(lts.hard_color_flows.front().amplitudes.size() == lts.hamp.size());
    for (const auto &h : indices(lts.hamp)) {
      RequireComplexNear(lts.hard_color_flows.front().amplitudes[h], lts.hamp[h]);
    }
  }
}

// Check rejected screening points cannot retain a previous color amplitude
TEST_CASE("Durham rejected screening points preserve zero color channels", "[gra::MDurham][color][screening]") {
  ModelParamRestoreGuard restore;
  const auto   tune = DurhamTestTune("rejected_screen", [](auto &j) { j["PARAM_DURHAM"]["JET_pt_min"] = 0.0; });
  auto         lts  = MakeToyDurhamGG();
  gra::MRandom rng;
  MDurham      durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
  REQUIRE(durham.DurhamQCD(lts, "MG5") > 0.0);
  REQUIRE_FALSE(lts.hard_color_flows.empty());
  const auto born = lts;
  for (const bool direct_loop : {false, true}) {
    lts                  = born;
    lts.screening.active = true;
    for (auto &flow : lts.hard_color_flows) { flow.screened_weight = 1.0; }
    lts.pfinal.clear();
    if (direct_loop) {
      MDurham::DurhamProjectedAmp hard(born.hamp.size());
      REQUIRE(durham.DQtloop(lts, hard) == Approx(0.0));
    } else {
      REQUIRE(durham.DurhamQCD(lts, "MG5") == Approx(0.0));
    }
    REQUIRE(lts.hamp.size() == born.hamp.size());
    REQUIRE(gra::SquaredNorm(lts.hamp) == Approx(0.0));
    REQUIRE(lts.hard_color_flows.size() == born.hard_color_flows.size());
    for (const auto &i : indices(lts.hard_color_flows)) {
      const auto &flow = lts.hard_color_flows[i];
      REQUIRE(flow.external.size() == born.hard_color_flows[i].external.size());
      for (const auto &leg : indices(flow.external)) {
        REQUIRE(flow.external[leg].color == born.hard_color_flows[i].external[leg].color);
        REQUIRE(flow.external[leg].anticolor == born.hard_color_flows[i].external[leg].anticolor);
      }
      REQUIRE(flow.amplitudes.size() == born.hard_color_flows[i].amplitudes.size());
      REQUIRE(gra::SquaredNorm(flow.amplitudes) == Approx(0.0));
      REQUIRE_FALSE(flow.screened_weight.has_value());
    }
  }
}

// Check direct generated calls reset status and cache validity on every point
TEST_CASE("Durham generated hard cache follows evaluation validity", "[gra::MDurham][MG2GRA][jets]") {
  ModelParamRestoreGuard restore;
  const auto             tune = DurhamTestTune("hard_cache", [](auto &j) {
    j["PARAM_DURHAM"]["JET_ALGO"]   = "anti-kt";
    j["PARAM_DURHAM"]["JET_pt_min"] = 0.5;
  });
  auto                   lts  = MakeToyDurhamGG();
  gra::MRandom           rng;
  MDurham                durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
  const auto             process = gra::amplitude::FindProcess("DURHAM", "gg_gg");
  REQUIRE(process.has_value());
  auto matrix_element = CreateDurhamMG5Process(*process);
  REQUIRE(matrix_element != nullptr);
  MDurham::DurhamProjectedAmp hard;
  durham.Dgg2Generated(lts, *matrix_element, hard);
  REQUIRE(lts.screening.durham.hard_cached);
  const auto original = lts;
  lts.decaytree.pop_back();
  durham.Dgg2Generated(lts, *matrix_element, hard);
  REQUIRE_FALSE(lts.screening.durham.hard_cached);
  REQUIRE_FALSE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
  lts = original;
  lts.decaytree.front().p4.SetPxPy(0.0, 0.0);
  durham.Dgg2Generated(lts, *matrix_element, hard);
  REQUIRE_FALSE(lts.screening.durham.hard_cached);
  REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
  for (const auto &channel : hard) { REQUIRE(gra::SquaredNorm(channel) == Approx(0.0)); }
  lts = original;
  durham.Dgg2Generated(lts, *matrix_element, hard);
  REQUIRE(lts.screening.durham.hard_cached);
  REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
}

// Check full loop amplitudes under common azimuthal rotations and beam exchange
TEST_CASE("Durham helicity prescriptions preserve azimuthal covariance", "[gra::MDurham][helicity][meson]") {
  ModelParamRestoreGuard restore;
  const bool log_qt = GENERATE(false, true);
  CAPTURE(log_qt);
  for (const std::string mode : {"transverse", "jzp"}) {
    const auto tune = DurhamTestTune(
        "rotation_" + mode, [](auto &j) { j["PARAM_DURHAM"]["JET_pt_min"] = 0.0; },
        [&](auto &j) {
          j["NUMERICS_DURHAM"]["HELICITY_PROJECTOR"] = mode;
          j["NUMERICS_DURHAM"]["qt_integrator"]      = "GL_IR";
          j["NUMERICS_DURHAM"]["log_qt"]             = log_qt;
          j["NUMERICS_DURHAM"]["N_phi"]              = 6;
          j["NUMERICS_DURHAM"]["N_qt"]               = 12;
        });
    for (const std::string process : {"MMbar", "MG5", "chic(0)", "chic(1)", "chic(2)"}) {
      CAPTURE(mode, process);
      auto         lts = process == "MMbar" ? MakeToyDurhamMesonPair(211, -211) : MakeToyDurhamGG();
      lts.PDG          = LoadedPDGTable();
      gra::MRandom rng;
      MDurham      durham(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
      const double born = durham.DurhamQCD(lts, process);
      REQUIRE(born > 0.0);
      if (process == "MMbar") {
        MDurham::DurhamProjectedAmp hard;
        REQUIRE_NOTHROW(durham.Dgg2MMbar(lts, hard));
        REQUIRE(hard.size() == 1);
        REQUIRE(gra::SquaredNorm(hard.front()) > 0.0);
      }
      const auto original = lts;
      // Include arbitrary rotations that do not permute a fixed azimuthal grid
      for (const double angle : {0.043, 0.137, gra::math::PI / 4.0, gra::math::PI / 2.0}) {
        lts = TransformDurhamPoint(original, [angle](auto &p) { p.RotateZ(angle); });
        UpdateToyDurhamDerivedKinematics(lts);
        REQUIRE(durham.DurhamQCD(lts, process) == Approx(born).epsilon(1e-9));
        if (process == "MMbar" || process == "chic(0)") { RequireComplexNear(lts.hamp[0], original.hamp[0], 1e-9); }
      }
      lts = TransformDurhamPoint(original, [](auto &p) { p.RotateY(gra::math::PI); });
      std::swap(lts.pbeam1, lts.pbeam2);
      std::swap(lts.pfinal[1], lts.pfinal[2]);
      UpdateToyDurhamDerivedKinematics(lts);
      REQUIRE(durham.DurhamQCD(lts, process) == Approx(born).epsilon(1e-9));
      lts = TransformDurhamPoint(original, [](auto &p) { p.SetPxPy(p.Px(), -p.Py()); });
      UpdateToyDurhamDerivedKinematics(lts);
      REQUIRE(durham.DurhamQCD(lts, process) == Approx(born).epsilon(1e-9));
    }
  }
}

// Check analytic meson helicities in the same physical basis as their sources
TEST_CASE("Durham meson source contractions are Lorentz covariant", "[gra::MDurham][meson][covariance]") {
  ModelParamRestoreGuard restore;
  const auto             tune =
      DurhamTestTune("meson_tensor", {}, [](auto &j) { j["NUMERICS_DURHAM"]["HELICITY_PROJECTOR"] = "transverse"; });
  const gra::M4Vec source1(0.37, -0.23, 0.0, 0.0);
  const gra::M4Vec source2(-0.19, 0.41, 0.0, 0.0);
  for (const auto &[first, second] : std::vector<std::pair<int, int>>{
           {211, -211}, {111, 111}, {321, -321}, {311, -311}, {221, 221}, {221, 331}, {331, 331}}) {
    CAPTURE(first, second);
    auto seed                            = MakeToyDurhamMesonPair(first, second);
    gra::MRandom rng;
    MDurham      durham(seed, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    // Contract the real public hard amplitude with arbitrary physical sources
    const auto contract = [&](const gra::LORENTZSCALAR &lts, const gra::M4Vec &q1, const gra::M4Vec &q2) {
      MDurham::DurhamProjectedAmp hard;
      durham.Dgg2MMbar(lts, hard);
      REQUIRE(hard.size() == 1);
      std::vector<gra::M4Vec> final = {lts.decaytree[0].p4, lts.decaytree[1].p4};
      gra::M4Vec              k1, k2;
      REQUIRE(gra::mg5helas::PrepareOnShellKinematics(lts, final, k1, k2));
      return gra::mg5helas::ContractTransverseSources(hard.front(), q1, q2, k1, k2);
    };
    const auto reference = contract(seed, source1, source2);
    REQUIRE(std::norm(reference) > 0.0);
    const std::vector<DurhamFrameTransform> transforms = {
        [](auto &p) { p = p.LorentzBoost({0.21, -0.13, 0.31}); },
        [](auto &p) { p = p.LorentzBoost({0.0, 0.0, std::tanh(2.0)}); },
        [](auto &p) {
          p.RotateX(0.41);
          p.RotateY(-0.72);
          p.RotateZ(1.19);
        },
        [](auto &p) { p = gra::M4Vec(-p.Px(), -p.Py(), -p.Pz(), p.E()); }};
    for (const auto &transform : transforms) {
      auto point = TransformDurhamPoint(seed, transform);
      UpdateToyDurhamDerivedKinematics(point);
      auto q1 = source1;
      auto q2 = source2;
      transform(q1);
      transform(q2);
      const auto value = contract(point, q1, q2);
      REQUIRE(std::abs(value - reference) < 2e-9 * std::abs(reference));
    }
    auto exchanged = seed;
    std::swap(exchanged.q1, exchanged.q2);
    const auto value = contract(exchanged, source2, source1);
    REQUIRE(std::abs(value - reference) < 2e-9 * std::abs(reference));
  }
}

// Distinguish rapidity acceptance from pseudorapidity for massive jets and cached screening
TEST_CASE("Durham massive jet acceptance uses rapidity", "[gra::MDurham][jets][flavour][screening]") {
  ModelParamRestoreGuard restore;
  gra::amplitude::ProcessRegistry registry(gra::amplitude::Processes("DURHAM"));
  const auto process = registry.MatchProcess(MakeToyDurhamQQbar(5).decaytree);
  REQUIRE(process.has_value());
  const double boost = GENERATE(-0.8, 0.8);
  const auto point = TransformDurhamPoint(DurhamBottomPoint(*process, false), [&](gra::M4Vec &p) {
    p = p.LorentzBoost({0.0, 0.0, std::tanh(boost)});
  });
  const double rap = std::max(std::abs(point.decaytree[0].p4.Rap()), std::abs(point.decaytree[1].p4.Rap()));
  const double eta = std::max(std::abs(point.decaytree[0].p4.Eta()), std::abs(point.decaytree[1].p4.Eta()));
  REQUIRE(rap > 0.0);
  REQUIRE(eta > rap + 1.0e-3);
  const std::string algorithm = GENERATE("anti-kt", "kt", "CA");
  for (const bool accepted : {true, false}) {
    CAPTURE(boost, algorithm, accepted, rap, eta);
    const auto tune = DurhamTestTune("jet_rap_" + algorithm, [&](auto &card) {
      card["PARAM_DURHAM"]["JET_ALGO"] = algorithm;
      card["PARAM_DURHAM"]["JET_R"] = 0.4;
      card["PARAM_DURHAM"]["JET_pt_min"] = 0.0;
      card["PARAM_DURHAM"]["JET_rap_max"] = accepted ? 0.5 * (rap + eta) : 0.5 * rap;
    });
    auto lts = point;
    lts.model_cache = std::make_shared<gra::MModelCache>(tune);
    gra::MRandom rng;
    MDurham durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
    for (const bool screening : {false, true}) {
      lts.screening.active = screening;
      const double amp2 = durham.DurhamQCD(lts, "MG5");
      REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
      CHECK(lts.screening.durham.jet_cuts_cached);
      CHECK(lts.screening.durham.jet_cuts_pass == accepted);
      if (accepted) { CHECK(amp2 > 0.0); }
      else { CHECK(amp2 == Approx(0.0)); }
    }
  }
}

// Check standard rapidity distances for a resolved bottom pair and gluon
TEST_CASE("Durham massive jet clustering uses rapidity distances", "[gra::MDurham][jets][flavour]") {
  ModelParamRestoreGuard          restore;
  gra::amplitude::ProcessRegistry registry(gra::amplitude::Processes("DURHAM"));
  auto                            seed = MakeToyDurhamQQbarG(5);
  const double                    mass = DurhamBottomMass(*registry.MatchProcess(seed.decaytree));
  for (const std::string algorithm : {"anti-kt", "kt", "CA"}) {
    const auto tune = DurhamTestTune("massive_jets_" + algorithm, [&](auto &j) {
      j["PARAM_DURHAM"]["JET_ALGO"]    = algorithm;
      j["PARAM_DURHAM"]["JET_R"]       = 0.4;
      j["PARAM_DURHAM"]["JET_pt_min"]  = 0.0;
      j["PARAM_DURHAM"]["JET_rap_max"] = 100.0;
    });
    for (const bool merge : {true, false}) {
      auto lts                = seed;
      // Each jet definition owns the cache for its own immutable tune
      lts.model_cache = std::make_shared<gra::MModelCache>(tune);
      lts.decaytree[0].p.mass = mass;
      lts.decaytree[1].p.mass = mass;
      lts.decaytree[0].p4.SetPxPyPzM(mass, 0.0, mass * std::sinh(1.0), mass);
      lts.decaytree[2].p4.SetPxPyPzM(mass, 0.0, mass * std::sinh(merge ? 0.5 : 0.0), 0.0);
      lts.decaytree[1].p4.SetPxPyPzM(-2.0 * mass, 0.0, -lts.decaytree[0].p4.Pz() - lts.decaytree[2].p4.Pz(), mass);
      const double energy = lts.decaytree[0].p4.E() + lts.decaytree[1].p4.E() + lts.decaytree[2].p4.E();
      lts.pbeam1.SetPxPyPzM(0.0, 0.0, 100.0, gra::PDG::mp);
      lts.pbeam2.SetPxPyPzM(0.0, 0.0, -100.0, gra::PDG::mp);
      const double proton_e = lts.pbeam1.E() - energy / 2.0;
      const double proton_z = std::sqrt(pow2(proton_e) - pow2(gra::PDG::mp) - 0.04);
      lts.pfinal[1]         = gra::M4Vec(0.2, 0.0, proton_z, proton_e);
      lts.pfinal[2]         = gra::M4Vec(-0.2, 0.0, -proton_z, proton_e);
      UpdateToyDurhamDerivedKinematics(lts);
      REQUIRE(std::abs(lts.decaytree[0].p4.Eta() - lts.decaytree[2].p4.Eta()) > 0.4);
      REQUIRE((std::abs(lts.decaytree[0].p4.Rap() - lts.decaytree[2].p4.Rap()) < 0.4) == merge);
      gra::MRandom rng;
      MDurham      durham(lts, tune, rng, gra::MDurham::ContinuumProcesses());
      const double amp2 = durham.DurhamQCD(lts, "MG5");
      REQUIRE(gra::mg5helas::EvaluationSucceeded(durham.EvaluationStatus()));
      if (merge) {
        REQUIRE(amp2 == Approx(0.0));
      } else {
        REQUIRE(amp2 > 0.0);
      }
    }
  }
}

TEST_CASE("Covariant MG5 transverse contraction reduces to the Durham projector", "[MG2GRA][gra::MDurham][helicity]") {
  const gra::M4Vec                          k1(0.0, 0.0, 50.0, 50.0);
  const gra::M4Vec                          k2(0.0, 0.0, -50.0, 50.0);
  const gra::M4Vec                          q1(0.7, -0.4, 0.0, 0.0);
  const gra::M4Vec                          q2(-0.2, 0.9, 0.0, 0.0);
  const std::array<std::complex<double>, 4> hard = {std::complex<double>(0.3, -0.2), std::complex<double>(-0.7, 0.4),
                                                    std::complex<double>(0.5, 0.8), std::complex<double>(-0.1, -0.6)};

  const double               dot   = q1.Px() * q2.Px() + q1.Py() * q2.Py();
  const double               cross = q1.Px() * q2.Py() - q1.Py() * q2.Px();
  const std::complex<double> p0    = -0.5 * dot;
  const std::complex<double> m0    = -0.5 * gra::math::zi * cross;
  const std::complex<double> p2 =
      0.5 * std::complex<double>(q1.Px() * q2.Px() - q1.Py() * q2.Py(), q1.Px() * q2.Py() + q1.Py() * q2.Px());
  const std::complex<double> m2 = std::conj(p2);
  const std::complex<double> expected =
      p0 * (hard[3] + hard[0]) + m0 * (hard[3] - hard[0]) + p2 * hard[1] + m2 * hard[2];

  const std::complex<double> projected = gra::mg5helas::ContractTransverseSources(hard, q1, q2, k1, k2);
  RequireComplexNear(projected, expected, 1e-12);
}

TEST_CASE(
    "All Durham hard kernels are Lorentz covariant after "
    "source contraction",
    "[gra::MDurham][MG5][registry][covariance]") {
  ModelParamRestoreGuard model_restore;
  const auto             tune = WriteModifiedDurhamTune("mg5_covariance", [](auto &general) {
    general["PARAM_DURHAM"]["JET_ALGO"]    = "none";
    general["PARAM_DURHAM"]["JET_pt_min"]  = 0.0;
    general["PARAM_DURHAM"]["JET_rap_max"] = 100.0;
  });
  gra::MODELPARAM             = tune.first;

  gra::LORENTZSCALAR config = MakeToyDurhamGG();
  gra::MRandom       rng;
  rng.SetSeed(23);
  gra::MDurham durham(config, gra::MModelTune::Load(tune.second), rng,
                      gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  std::vector<gra::amplitude::Process> processes = gra::amplitude::Processes("DURHAM");
  processes.push_back(gra::AMP_gg_yy::Definition());
  REQUIRE_FALSE(processes.empty());
  const gra::M4Vec source1(0.37, -0.23, 0.0, 0.0);
  const gra::M4Vec source2(-0.19, 0.41, 0.0, 0.0);

  for (const auto &process : processes) {
    CAPTURE(process.process_name);
    std::unique_ptr<gra::DurhamMG5Process> matrix_element;
    if (process.matrix_element_form == gra::amplitude::MatrixElementForm::Analytic) {
      REQUIRE_FALSE(gra::HasDurhamMG5Process(process));
      matrix_element = std::make_unique<gra::AMP_gg_yy>(config.model_cache->Tune().SM());
    } else {
      REQUIRE(gra::HasDurhamMG5Process(process));
      matrix_element = gra::CreateDurhamMG5Process(process);
    }
    REQUIRE(matrix_element != nullptr);
    gra::LORENTZSCALAR seed = DurhamRegistryPoint(process);
    seed.model_cache        = config.model_cache;
    seed.GlobalSudakovPtr   = config.GlobalSudakovPtr;

    gra::LORENTZSCALAR reference      = seed;
    const double       reference_norm = ContractedDurhamNorm(durham, *matrix_element, reference, source1, source2);
    REQUIRE(reference_norm > 0.0);

    // Compare one transformed physical source contraction with the reference
    const auto require_covariance = [&](const DurhamFrameTransform &transform, const double tolerance,
                                        const std::string &frame) {
      gra::LORENTZSCALAR transformed         = TransformDurhamPoint(seed, transform);
      gra::M4Vec         transformed_source1 = source1;
      gra::M4Vec         transformed_source2 = source2;
      transform(transformed_source1);
      transform(transformed_source2);
      const double transformed_norm =
          ContractedDurhamNorm(durham, *matrix_element, transformed, transformed_source1, transformed_source2);
      if (std::abs(process.stable_pdgs.front()) == 5) {
        constexpr std::array<std::size_t, 2> bottom = {0, 1};
        for (const auto &i : bottom) {
          CAPTURE(i);
          const double mass2 = seed.decaytree[i].p.mass * seed.decaytree[i].p.mass;
          REQUIRE(seed.decaytree[i].p4.M2() == Approx(mass2).epsilon(2.0e-12));
          REQUIRE(transformed.decaytree[i].p4.M2() == Approx(mass2).epsilon(2.0e-8));
        }
      }
      const double relative_difference = std::abs(transformed_norm - reference_norm) / reference_norm;
      CAPTURE(frame, reference_norm, transformed_norm, relative_difference);
      REQUIRE(transformed_norm == Approx(reference_norm).epsilon(tolerance));
    };

    // Apply one proper rotation which moves the incoming axis away from z
    const DurhamFrameTransform rotation = [](gra::M4Vec &momentum) {
      momentum.RotateX(0.41);
      momentum.RotateY(-0.72);
      momentum.RotateZ(1.19);
    };
    require_covariance(rotation, 2.0e-9, "rotation");

    // Apply one noncollinear boost to all hard and source momenta
    const DurhamFrameTransform noncollinear_boost = [](gra::M4Vec &momentum) {
      momentum = momentum.LorentzBoost({0.21, -0.13, 0.31});
    };
    require_covariance(noncollinear_boost, 5.0e-9, "noncollinear boost");

    for (const double rapidity : {-5.0, 5.0}) {
      // Apply a large longitudinal boost without changing any hard invariant
      const DurhamFrameTransform longitudinal_boost = [rapidity](gra::M4Vec &momentum) {
        momentum = momentum.LorentzBoost({0.0, 0.0, std::tanh(rapidity)});
      };
      require_covariance(longitudinal_boost, 2.0e-8, "rapidity " + std::to_string(rapidity));
    }
  }
}

TEST_CASE("Durham amplitudes are independent of the MRegge decay barrier", "[gra::MDurham][helicity][isolation]") {
  using gra::aux::indices;
  ModelParamRestoreGuard model_restore;
  const auto             relaxed_tune = WriteRelaxedDurhamTune("decay_barrier_isolation");
  gra::MODELPARAM                     = relaxed_tune.first;

  gra::LORENTZSCALAR enabled     = MakeToyDurhamGG();
  enabled.process.DECAY_BARRIER  = true;
  gra::LORENTZSCALAR disabled    = enabled;
  disabled.process.DECAY_BARRIER = false;
  gra::MRandom enabled_rng;
  gra::MRandom disabled_rng;
  enabled_rng.SetSeed(9);
  disabled_rng.SetSeed(9);
  MDurham enabled_durham(enabled, gra::MModelTune::Load(relaxed_tune.second), enabled_rng,
                         gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  MDurham disabled_durham(disabled, gra::MModelTune::Load(relaxed_tune.second), disabled_rng,
                          gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  MDurham::DurhamProjectedAmp enabled_amp(4);
  MDurham::DurhamProjectedAmp disabled_amp(4);
  for (auto &channel : enabled_amp) { channel.fill(0.0); }
  for (auto &channel : disabled_amp) { channel.fill(0.0); }
  EvaluateGeneratedDurham(enabled_durham, enabled, enabled_amp);
  EvaluateGeneratedDurham(disabled_durham, disabled, disabled_amp);
  REQUIRE(enabled_amp.size() == disabled_amp.size());
  for (const auto &c : indices(enabled_amp)) {
    for (const auto &h : indices(enabled_amp[c])) { RequireComplexNear(enabled_amp[c][h], disabled_amp[c][h], 1e-12); }
  }
}

TEST_CASE("MDurham reads model parameters and Durham loop numerics separately", "[gra::MDurham][params]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  gra::MDurhamParam param;
  REQUIRE_NOTHROW(param.ConfigureFromJson(modelfile, gra::aux::GetInputData(modelfile)));
  param.muF                       = 0.25;
  param.muR                       = 0.75;
  const gra::MDurhamScales scales = gra::DurhamCentralScales(20.0, param);
  REQUIRE(scales.muF == Approx(5.0));
  REQUIRE(scales.muR == Approx(15.0));

  gra::MDurhamNumerics numerics;
  const std::string    numerics_file = gra::ResolveModelDataFile("TUNE0", "NUMERICS.json");
  REQUIRE_NOTHROW(numerics.ConfigureFromJson(numerics_file, gra::aux::GetInputData(numerics_file)));
  REQUIRE(gra::math::PolarNodeCount(numerics.loop) == numerics.loop.radial_intervals);
  REQUIRE(numerics.loop.r_max == Approx(std::sqrt(numerics.qt2_MAX)).epsilon(1e-14));
}

// Check that the production rule enforces all three virtuality cuts
TEST_CASE("MDurham applies one inclusive cutoff to all loop virtualities", "[gra::MDurham][params][infrared]") {
  const auto   tune   = gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json"));
  const auto   config = gra::ReadDurhamConfig(*tune);
  const double cut    = config->param.loop_q2_cut;
  const std::vector<std::array<double, 2>> centres = {{0.3, 0.1}, {-0.4, 0.2}};
  const auto                               nodes   = config->loop_const->Nodes(centres, std::sqrt(cut));
  REQUIRE_FALSE(nodes.empty());
  for (const auto& node : nodes) {
    REQUIRE(node.x * node.x + node.y * node.y >= cut - 1e-12);
    for (const auto& centre : centres) { REQUIRE(pow2(node.x - centre[0]) + pow2(node.y - centre[1]) >= cut - 1e-12); }
  }
}

TEST_CASE("MDurham keeps the alpha_s cache separate from the factorization scale", "[gra::MDurham][params]") {
  ModelParamRestoreGuard restore;
  const auto             relaxed_tune = WriteModifiedDurhamTune("scale_cache", [](auto& j) {
    j["PARAM_DURHAM"]["muF"]      = 0.5;
    j["PARAM_DURHAM"]["JET_ALGO"] = "none";
  });
  gra::MODELPARAM                     = relaxed_tune.first;

  gra::LORENTZSCALAR seed = MakeToyDurhamGG();
  gra::MRandom       rng;
  rng.SetSeed(1);
  MDurham durham(seed, gra::MModelTune::Load(relaxed_tune.second), rng,
                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  gra::LORENTZSCALAR expected     = seed;
  gra::LORENTZSCALAR stale        = seed;
  const double       central_mass = std::sqrt(seed.s_hat);
  gra::MDurhamParam  param;
  param.ConfigureFromJson(relaxed_tune.second, gra::aux::GetInputData(relaxed_tune.second));
  const gra::MDurhamScales scales = gra::DurhamCentralScales(central_mass, param);

  expected.alphaQCD = expected.GlobalSudakovPtr->AlphaS_Q2(scales.muR * scales.muR);
  expected.muR      = scales.muR;
  expected.muF      = scales.muF;
  expected.scalup   = scales.muF;

  stale.alphaQCD = 1.0e-3;
  stale.muR      = scales.muF;
  stale.muF      = scales.muF;
  stale.scalup   = scales.muF;

  MDurham::DurhamProjectedAmp expected_amp;
  MDurham::DurhamProjectedAmp stale_amp;
  EvaluateGeneratedDurham(durham, expected, expected_amp);
  EvaluateGeneratedDurham(durham, stale, stale_amp);
  REQUIRE(expected_amp.size() == stale_amp.size());
  double norm = 0.0;
  for (std::size_t i = 0; i < expected_amp.size(); ++i) {
    for (std::size_t j = 0; j < expected_amp[i].size(); ++j) {
      norm += std::norm(expected_amp[i][j]);
      RequireComplexNear(stale_amp[i][j], expected_amp[i][j], 1e-12);
    }
  }
  REQUIRE(norm > 0.0);
}

TEST_CASE("MDurham rejects a perturbative loop cutoff below Q0 squared", "[gra::MDurham][params][infrared]") {
  ModelParamRestoreGuard restore;
  const std::string      tune = WriteModifiedSudakovModelTune(
           "perturbative_low_loop_cut",
           [](auto &j) {
        j["PARAM_SKEWED_UGD"]["mode"]        = "PERTURBATIVE_ONLY";
        j["PARAM_SKEWED_UGD"]["rkhs_radius"] = 0.0;
        j["PARAM_DURHAM"]["loop_q2_cut"]  = 0.4;
      },
           [](auto &) {});
  gra::MODELPARAM = tune;
  gra::LORENTZSCALAR lts;
  lts.sqrt_s    = 13000.0;
  lts.LHAPDFSET = "MMHT2014lo68cl";
  gra::MRandom rng;
  rng.SetSeed(1);
  REQUIRE_THROWS(MDurham(lts, gra::MModelTune::Load(gra::ResolveModelDataFile(tune, "GENERAL.json")), rng,
                         gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test")));
}

// Check the transverse integration measure and shared configuration for both radial maps
TEST_CASE("Durham configuration constructs safely from one tune across threads", "[gra::MDurham][params][threading]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  const bool log_qt = GENERATE(false, true);
  CAPTURE(log_qt);

  const auto model_tune = DurhamTestTune("radial_map", {}, [log_qt](auto &j) {
    j["NUMERICS_DURHAM"]["log_qt"] = log_qt;
  });
  gra::MModelCache cache(model_tune);
  const auto       first  = gra::GetDurhamConfig(cache);
  const auto       second = gra::GetDurhamConfig(cache);
  REQUIRE(first == second);
  REQUIRE(first->param.initialized);
  REQUIRE(first->numerics.initialized);
  REQUIRE(first->loop_const.has_value());
  double transverse_measure = 0.0;
  for (const auto& node : first->loop_const->Nodes({}, 0.0)) { transverse_measure += node.weight; }
  REQUIRE(transverse_measure ==
          Approx(gra::math::PI * (first->numerics.qt2_MAX - first->param.loop_q2_cut)).epsilon(1e-10));

  constexpr std::size_t                                          nthreads = 8;
  std::vector<std::shared_ptr<const gra::MDurham::DurhamConfig>> handles(nthreads);
  std::vector<std::thread>                                       workers;
  workers.reserve(nthreads);

  for (std::size_t i = 0; i < nthreads; ++i) {
    workers.emplace_back([i, &handles, &cache] { handles[i] = gra::GetDurhamConfig(cache); });
  }
  for (auto &worker : workers) { worker.join(); }

  for (const auto &handle : handles) { REQUIRE(handle == first); }
}

// Reject logarithmic integration without a positive nonempty momentum range at initialization
TEST_CASE("Durham logarithmic quadrature rejects invalid bounds", "[gra::MDurham][params][quadrature]") {
  ModelParamRestoreGuard restore;
  const auto tune = DurhamTestTune("log_zero_cut", [](auto &j) {
    j["PARAM_DURHAM"]["loop_q2_cut"] = 0.0;
  }, [](auto &j) { j["NUMERICS_DURHAM"]["log_qt"] = true; });
  REQUIRE_THROWS_AS(gra::ReadDurhamConfig(*tune), std::invalid_argument);

  const auto empty = DurhamTestTune("log_empty_range", {}, [](auto &j) {
    j["NUMERICS_DURHAM"]["log_qt"] = true;
    j["NUMERICS_DURHAM"]["qt2_MAX"] = 0.01;
  });
  REQUIRE_THROWS_AS(gra::ReadDurhamConfig(*empty), std::invalid_argument);
}

// Require Durham parameters to come from one complete tune block
TEST_CASE("Durham configuration rejects an incomplete tune block", "[gra::MDurham][params][snapshot]") {
  const auto tune =
      WriteModifiedDurhamTune("missing_parameters", [](auto &general) { general.erase("PARAM_DURHAM"); });
  REQUIRE_THROWS(gra::ReadDurhamConfig(*gra::MModelTune::Load(tune.second)));
}

TEST_CASE("MDurham rejects a Sudakov from a different SOFT snapshot", "[gra::MDurham][sudakov][snapshot]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto first_tune =
      WriteModifiedDurhamTune("first_soft", [](auto &general) { general["PARAM_REGGE"]["omega"]["MP"] = 0.81; });
  const auto second_tune =
      WriteModifiedDurhamTune("second_soft", [](auto &general) { general["PARAM_REGGE"]["omega"]["MP"] = 0.97; });
  const auto first_model  = gra::MModelTune::Load(first_tune.second);
  const auto second_model = gra::MModelTune::Load(second_tune.second);
  REQUIRE(first_model != second_model);

  gra::LORENTZSCALAR lts = MakeToyDurhamGG();
  lts.GlobalSudakovPtr   = lts.model_cache->sudakov.GetSudakov(lts.sqrt_s, lts.LHAPDFSET, first_model->Soft());
  REQUIRE(lts.GlobalSudakovPtr->SoftModelHandle() == first_model->Soft());

  gra::MRandom rng;
  rng.SetSeed(1);
  REQUIRE_THROWS(MDurham(lts, second_model, rng,
                         gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "snapshot_mismatch")));
}

TEST_CASE("MDurham rejects invalid Durham model parameters", "[gra::MDurham][params]") {
  auto require_invalid = [](const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate,
                            const std::string &message) {
    auto j = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    mutate(j);

    const std::string path = "tmp/graniitti_durham_invalid_model_" + suffix + ".json";
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    gra::MDurhamParam param;
    REQUIRE_THROWS(param.ConfigureFromJson(path, gra::aux::GetInputData(path)));
  };
  auto require_valid_algo = [](const std::string &algo) {
    auto j                           = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
    j["PARAM_DURHAM"]["JET_ALGO"] = algo;

    const std::string path = "tmp/graniitti_durham_valid_jet_algo_" + algo + ".json";
    std::ofstream     out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    gra::MDurhamParam param;
    REQUIRE_NOTHROW(param.ConfigureFromJson(path, gra::aux::GetInputData(path)));
    REQUIRE(param.JET_ALGO == algo);
  };

  require_invalid(
      "pdf_scale", [](auto &j) { j["PARAM_DURHAM"]["PDF_scale"] = "BAD"; }, "unknown PDF_scale");
  require_invalid(
      "muF", [](auto &j) { j["PARAM_DURHAM"]["muF"] = 0.0; }, "muF");
  require_invalid(
      "muR", [](auto &j) { j["PARAM_DURHAM"]["muR"] = 0.0; }, "muR");
  require_invalid(
      "loop_q2_cut", [](auto &j) { j["PARAM_DURHAM"]["loop_q2_cut"] = -0.1; }, "loop_q2_cut");
  require_invalid(
      "maxcos", [](auto &j) { j["PARAM_DURHAM"]["MAXCOS"] = 1.1; }, "MAXCOS");
  require_invalid(
      "maxcos_one", [](auto &j) { j["PARAM_DURHAM"]["MAXCOS"] = 1.0; }, "MAXCOS");
  require_invalid(
      "meson_pt", [](auto &j) { j["PARAM_DURHAM"]["MESON_pt_min"] = 0.0; }, "MESON_pt_min");
  require_invalid(
      "f_pi", [](auto &j) { j["PARAM_DURHAM"]["f_pi"] = 0.0; }, "f_pi");
  require_valid_algo("none");
  require_valid_algo("anti-kt");
  require_valid_algo("kt");
  require_valid_algo("CA");
  require_invalid(
      "jet_algo", [](auto &j) { j["PARAM_DURHAM"]["JET_ALGO"] = "antikt"; }, "JET_ALGO");
  require_invalid(
      "jet_r", [](auto &j) { j["PARAM_DURHAM"]["JET_R"] = 0.0; }, "JET_R");
  require_invalid(
      "jet_pt", [](auto &j) { j["PARAM_DURHAM"]["JET_pt_min"] = -1.0; }, "JET_pt_min");
  require_invalid(
      "jet_rap", [](auto &j) { j["PARAM_DURHAM"]["JET_rap_max"] = 0.0; }, "JET_rap_max");
  require_invalid(
      "eta_f8", [](auto &j) { j["PARAM_DURHAM"]["f_eta8_over_fpi"] = 0.0; }, "eta decay-constant");
  require_invalid(
      "eta_f0", [](auto &j) { j["PARAM_DURHAM"]["f_eta0_over_fpi"] = -1.0; }, "eta decay-constant");
  require_invalid(
      "eta_theta8", [](auto &j) { j["PARAM_DURHAM"]["eta_theta8_deg"] = "bad"; }, "type must be number");
  require_invalid(
      "eta_theta1", [](auto &j) { j["PARAM_DURHAM"]["eta_theta1_deg"] = "bad"; }, "type must be number");
  require_invalid(
      "chic_width_zero", [](auto &j) { j["PARAM_DURHAM"]["chic0_gg_width_fraction"] = 0.0; },
      "chic0_gg_width_fraction");
  require_invalid(
      "chic_width_large", [](auto &j) { j["PARAM_DURHAM"]["chic0_gg_width_fraction"] = 1.01; },
      "chic0_gg_width_fraction");
}

TEST_CASE("MDurham validates Durham loop numerics from NUMERICS.json", "[gra::MDurham][params]") {
  ModelParamRestoreGuard restore;

  auto require_invalid = [](const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate,
                            const std::string &message) {
    auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
    mutate(j);

    const std::filesystem::path dir = "tmp/graniitti_durham_invalid_numerics_" + suffix;
    std::filesystem::create_directories(dir);
    const std::filesystem::path path = dir / "NUMERICS.json";
    std::ofstream               out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    gra::MODELPARAM = dir.string();
    gra::MDurhamNumerics numerics;
    REQUIRE_THROWS(numerics.ConfigureFromJson(path.string(), gra::aux::GetInputData(path.string())));
  };

  require_invalid(
      "nqt", [](auto &j) { j["NUMERICS_DURHAM"]["N_qt"] = 0; }, "N_qt and N_phi");
  require_invalid(
      "nphi", [](auto &j) { j["NUMERICS_DURHAM"]["N_phi"] = 0; }, "N_qt and N_phi");
  for (const std::string key : {"N_qt", "N_phi", "N_x"}) {
    for (const auto& value : nlohmann::json::parse(R"([2.5,4294967304,true,"8"])")) {
      require_invalid("count_" + key, [&](auto& j) { j["NUMERICS_DURHAM"][key] = value; }, "integer");
    }
  }
  require_invalid(
      "bad_qt", [](auto &j) { j["NUMERICS_DURHAM"]["qt_integrator"] = "BadQT"; }, "Unknown qt_integrator");
  require_invalid(
      "bad_phi", [](auto &j) { j["NUMERICS_DURHAM"]["phi_integrator"] = "BadPhi"; }, "Unknown phi_integrator");
  require_invalid(
      "bad_helicity_projector", [](auto &j) { j["NUMERICS_DURHAM"]["HELICITY_PROJECTOR"] = "BAD"; },
      "Unknown HELICITY_PROJECTOR");
  require_invalid(
      "boole_count",
      [](auto &j) {
        j["NUMERICS_DURHAM"]["qt_integrator"] = "Boole";
        j["NUMERICS_DURHAM"]["N_qt"]          = 10;
      },
      "not compatible with qt_integrator Boole");
  require_invalid(
      "qt2_nonpositive", [](auto &j) { j["NUMERICS_DURHAM"]["qt2_MAX"] = 0.0; }, "qt2_MAX must be positive");
  require_invalid(
      "qt2_suda_max", [](auto &j) { j["NUMERICS_DURHAM"]["qt2_MAX"] = 1000.0; }, "must not exceed NUMERICS_SUDAKOV");
}

TEST_CASE("MDurham accepts infrared GL and Newton-Cotes loop node counts", "[gra::MDurham][params]") {
  ModelParamRestoreGuard restore;

  auto read_modified = [](const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate) {
    auto j = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
    mutate(j);

    const std::filesystem::path dir = "tmp/graniitti_durham_valid_numerics_" + suffix;
    std::filesystem::create_directories(dir);
    const std::filesystem::path path = dir / "NUMERICS.json";
    std::ofstream               out(path);
    REQUIRE(out.good());
    out << j.dump(2);
    out.close();

    gra::MODELPARAM = dir.string();
    gra::MDurhamNumerics numerics;
    REQUIRE_NOTHROW(numerics.ConfigureFromJson(path.string(), gra::aux::GetInputData(path.string())));
    return numerics;
  };

  const auto gl = read_modified("gl_trap", [](auto &j) {
    j["NUMERICS_DURHAM"]["qt_integrator"]  = "GL";
    j["NUMERICS_DURHAM"]["phi_integrator"] = "Trap";
    j["NUMERICS_DURHAM"]["N_qt"]           = 5;
    j["NUMERICS_DURHAM"]["N_phi"]          = 7;
  });
  REQUIRE(gra::math::PolarNodeCount(gl.loop) == 5);
  REQUIRE(gl.loop.azimuth_nodes == 7);
  REQUIRE(gl.helicity_projector == "transverse");

  const auto jzp = read_modified("jzp_projector", [](auto &j) { j["NUMERICS_DURHAM"]["HELICITY_PROJECTOR"] = "jzp"; });
  REQUIRE(jzp.helicity_projector == "jzp");

  const auto gl_ir = read_modified("gl_ir_trap", [](auto &j) {
    j["NUMERICS_DURHAM"]["qt_integrator"]  = "GL_IR";
    j["NUMERICS_DURHAM"]["phi_integrator"] = "Trap";
    j["NUMERICS_DURHAM"]["N_qt"]           = 5;
    j["NUMERICS_DURHAM"]["N_phi"]          = 7;
  });
  REQUIRE(gra::math::PolarNodeCount(gl_ir.loop) == 5);
  REQUIRE(gl_ir.loop.azimuth_nodes == 7);

  const auto boole = read_modified("boole_trap", [](auto &j) {
    j["NUMERICS_DURHAM"]["qt_integrator"]  = "Boole";
    j["NUMERICS_DURHAM"]["phi_integrator"] = "Trap";
    j["NUMERICS_DURHAM"]["N_qt"]           = 12;
    j["NUMERICS_DURHAM"]["N_phi"]          = 7;
  });
  REQUIRE(gra::math::PolarNodeCount(boole.loop) == 13);
  REQUIRE(boole.loop.azimuth_nodes == 7);

  const auto simpson13 = read_modified("simpson13_trap", [](auto &j) {
    j["NUMERICS_DURHAM"]["qt_integrator"]  = "1/3";
    j["NUMERICS_DURHAM"]["phi_integrator"] = "Trap";
    j["NUMERICS_DURHAM"]["N_qt"]           = 10;
    j["NUMERICS_DURHAM"]["N_phi"]          = 7;
  });
  REQUIRE(gra::math::PolarNodeCount(simpson13.loop) == 11);
  REQUIRE(simpson13.loop.azimuth_nodes == 7);

  const auto simpson38 = read_modified("simpson38_trap", [](auto &j) {
    j["NUMERICS_DURHAM"]["qt_integrator"]  = "3/8";
    j["NUMERICS_DURHAM"]["phi_integrator"] = "Trap";
    j["NUMERICS_DURHAM"]["N_qt"]           = 12;
    j["NUMERICS_DURHAM"]["N_phi"]          = 7;
  });
  REQUIRE(gra::math::PolarNodeCount(simpson38.loop) == 13);
  REQUIRE(simpson38.loop.azimuth_nodes == 7);
}

TEST_CASE("MSubProc binds Durham RNG lazily and resets it across copies", "[gra::MSubProc][gra::MDurham]") {
  ModelParamRestoreGuard model_restore;
  const auto             relaxed_tune = WriteRelaxedDurhamTune("subproc_rng");
  gra::MODELPARAM                     = relaxed_tune.first;

  ToyHelicityProcess master;
  gra::MODELPARAM = relaxed_tune.first;
  master.SetModelTune(gra::MModelTune::Load(relaxed_tune.second));
  master.state.lts = MakeToyDurhamGG();
  master.ProcPtr.Initialize("gg", "QCD");
  master.InitializeProcessAmplitude();

  gra::LORENTZSCALAR lts     = master.state.lts;
  gra::MSubProc      subproc = master.ProcPtr;
  REQUIRE_THROWS(subproc.GetBareAmplitude2(lts));

  gra::MRandom rng;
  rng.SetSeed(7);
  subproc.BindRandom(rng);
  const double amp2 = subproc.GetBareAmplitude2(lts);
  REQUIRE(std::isfinite(amp2));
  REQUIRE(amp2 >= 0.0);
  REQUIRE_FALSE(lts.hard_color_flows.empty());
  REQUIRE(lts.decaytree[0].p.color_flow.empty());
  subproc.SampleColorFlow(lts);
  REQUIRE_FALSE(lts.decaytree[0].p.color_flow.empty());

  gra::MSubProc copied = subproc;
  REQUIRE_THROWS(copied.GetBareAmplitude2(lts));
}

TEST_CASE("MDurham gg projector matches Durham-color-averaged SU(3) contractions", "[gra::MDurham][color]") {
  using Matrix3 = std::array<std::array<std::complex<double>, 3>, 3>;

  auto zero_matrix = []() {
    Matrix3 out{};
    for (auto &row : out) {
      for (auto &value : row) { value = 0.0; }
    }
    return out;
  };
  auto multiply = [](const Matrix3 &A, const Matrix3 &B) {
    Matrix3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        out[i][j] = 0.0;
        for (std::size_t k = 0; k < 3; ++k) { out[i][j] += A[i][k] * B[k][j]; }
      }
    }
    return out;
  };
  auto trace = [](const Matrix3 &A) { return A[0][0] + A[1][1] + A[2][2]; };

  std::array<Matrix3, 8> T{};
  T[0]       = zero_matrix();
  T[0][0][1] = 0.5;
  T[0][1][0] = 0.5;
  T[1]       = zero_matrix();
  T[1][0][1] = std::complex<double>(0.0, -0.5);
  T[1][1][0] = std::complex<double>(0.0, 0.5);
  T[2]       = zero_matrix();
  T[2][0][0] = 0.5;
  T[2][1][1] = -0.5;
  T[3]       = zero_matrix();
  T[3][0][2] = 0.5;
  T[3][2][0] = 0.5;
  T[4]       = zero_matrix();
  T[4][0][2] = std::complex<double>(0.0, -0.5);
  T[4][2][0] = std::complex<double>(0.0, 0.5);
  T[5]       = zero_matrix();
  T[5][1][2] = 0.5;
  T[5][2][1] = 0.5;
  T[6]       = zero_matrix();
  T[6][1][2] = std::complex<double>(0.0, -0.5);
  T[6][2][1] = std::complex<double>(0.0, 0.5);
  T[7]       = zero_matrix();
  T[7][0][0] = 1.0 / (2.0 * std::sqrt(3.0));
  T[7][1][1] = 1.0 / (2.0 * std::sqrt(3.0));
  T[7][2][2] = -1.0 / std::sqrt(3.0);

  const std::array<std::array<int, 4>, 6> perms = {
      {{{0, 1, 2, 3}}, {{0, 1, 3, 2}}, {{0, 2, 1, 3}}, {{0, 2, 3, 1}}, {{0, 3, 1, 2}}, {{0, 3, 2, 1}}}};

  // Incoming side is the Durham color average, final side is a normalized
  // singlet
  const double kmr_incoming_color_average = 1.0 / std::sqrt(8.0);
  const double gg_tensor_norm             = kmr_incoming_color_average / 8.0;

  std::vector<double> singlet_coeff(6, 0.0);
  for (std::size_t i = 0; i < perms.size(); ++i) {
    std::complex<double> coeff = 0.0;
    for (int a = 0; a < 8; ++a) {
      for (int c = 0; c < 8; ++c) {
        const auto ABCD = multiply(multiply(multiply(T[a], T[a]), T[c]), T[c]);
        const auto ACBD = multiply(multiply(multiply(T[a], T[c]), T[a]), T[c]);

        if (perms[i] == std::array<int, 4>{{0, 1, 2, 3}} || perms[i] == std::array<int, 4>{{0, 1, 3, 2}} ||
            perms[i] == std::array<int, 4>{{0, 2, 3, 1}} || perms[i] == std::array<int, 4>{{0, 3, 2, 1}}) {
          coeff += trace(ABCD) * gg_tensor_norm;
        } else {
          coeff += trace(ACBD) * gg_tensor_norm;
        }
      }
    }
    singlet_coeff[i] = coeff.real();
  }

  const std::vector<double> expected = {kmr_incoming_color_average * 2.0 / 3.0, kmr_incoming_color_average * 2.0 / 3.0,
                                        -kmr_incoming_color_average / 12.0,     kmr_incoming_color_average * 2.0 / 3.0,
                                        -kmr_incoming_color_average / 12.0,     kmr_incoming_color_average * 2.0 / 3.0};
  for (std::size_t i = 0; i < expected.size(); ++i) { REQUIRE(singlet_coeff[i] == Approx(expected[i]).epsilon(1e-12)); }
}

TEST_CASE(
    "MDurham qqbar and ggg projectors match Durham-color-averaged SU(3) "
    "contractions",
    "[gra::MDurham][color]") {
  using Matrix3 = std::array<std::array<std::complex<double>, 3>, 3>;

  auto zero_matrix = []() {
    Matrix3 out{};
    for (auto &row : out) {
      for (auto &value : row) { value = 0.0; }
    }
    return out;
  };
  auto add = [](const Matrix3 &A, const Matrix3 &B) {
    Matrix3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) { out[i][j] = A[i][j] + B[i][j]; }
    }
    return out;
  };
  auto subtract = [](const Matrix3 &A, const Matrix3 &B) {
    Matrix3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) { out[i][j] = A[i][j] - B[i][j]; }
    }
    return out;
  };
  auto multiply = [](const Matrix3 &A, const Matrix3 &B) {
    Matrix3 out{};
    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        out[i][j] = 0.0;
        for (std::size_t k = 0; k < 3; ++k) { out[i][j] += A[i][k] * B[k][j]; }
      }
    }
    return out;
  };
  auto trace = [](const Matrix3 &A) { return A[0][0] + A[1][1] + A[2][2]; };

  std::array<Matrix3, 8> T{};
  T[0]       = zero_matrix();
  T[0][0][1] = 0.5;
  T[0][1][0] = 0.5;
  T[1]       = zero_matrix();
  T[1][0][1] = std::complex<double>(0.0, -0.5);
  T[1][1][0] = std::complex<double>(0.0, 0.5);
  T[2]       = zero_matrix();
  T[2][0][0] = 0.5;
  T[2][1][1] = -0.5;
  T[3]       = zero_matrix();
  T[3][0][2] = 0.5;
  T[3][2][0] = 0.5;
  T[4]       = zero_matrix();
  T[4][0][2] = std::complex<double>(0.0, -0.5);
  T[4][2][0] = std::complex<double>(0.0, 0.5);
  T[5]       = zero_matrix();
  T[5][1][2] = 0.5;
  T[5][2][1] = 0.5;
  T[6]       = zero_matrix();
  T[6][1][2] = std::complex<double>(0.0, -0.5);
  T[6][2][1] = std::complex<double>(0.0, 0.5);
  T[7]       = zero_matrix();
  T[7][0][0] = 1.0 / (2.0 * std::sqrt(3.0));
  T[7][1][1] = 1.0 / (2.0 * std::sqrt(3.0));
  T[7][2][2] = -1.0 / std::sqrt(3.0);

  // Incoming side is the Durham color average, final channels are normalized
  // singlets
  const double kmr_incoming_color_average = 1.0 / std::sqrt(8.0);
  const double qqbar_tensor_norm          = 1.0 / (8.0 * std::sqrt(3.0));

  std::array<std::complex<double>, 2> qqbar_coeff{};
  for (std::size_t i = 0; i < qqbar_coeff.size(); ++i) {
    std::complex<double> coeff = 0.0;
    for (int a = 0; a < 8; ++a) { coeff += trace(multiply(T[a], T[a])) * qqbar_tensor_norm; }
    qqbar_coeff[i] = coeff;
  }
  REQUIRE(qqbar_coeff[0].real() == Approx(kmr_incoming_color_average * 2.0 / std::sqrt(6.0)).epsilon(1e-12));
  REQUIRE(qqbar_coeff[1].real() == Approx(kmr_incoming_color_average * 2.0 / std::sqrt(6.0)).epsilon(1e-12));
  REQUIRE(qqbar_coeff[0].imag() == Approx(0.0).margin(1e-12));
  REQUIRE(qqbar_coeff[1].imag() == Approx(0.0).margin(1e-12));

  auto fabc = [&](int a, int b, int c) {
    const auto comm = subtract(multiply(T[b], T[c]), multiply(T[c], T[b]));
    return ((2.0 / gra::math::zi) * trace(multiply(T[a], comm))).real();
  };
  auto dabc = [&](int a, int b, int c) {
    const auto anti = add(multiply(T[b], T[c]), multiply(T[c], T[b]));
    return (2.0 * trace(multiply(T[a], anti))).real();
  };

  const double              inv_sqrt3  = kmr_incoming_color_average / std::sqrt(3.0);
  const double              inv8_sqrt3 = inv_sqrt3 / 8.0;
  const double              d_leading  = kmr_incoming_color_average * std::sqrt(5.0 / 27.0);
  const double              d_sublead  = d_leading / 8.0;
  const std::vector<double> expected_f = {inv_sqrt3,   -inv_sqrt3,  -inv_sqrt3,  inv_sqrt3,  inv_sqrt3,   -inv_sqrt3,
                                          -inv8_sqrt3, inv8_sqrt3,  -inv8_sqrt3, inv_sqrt3,  inv8_sqrt3,  -inv_sqrt3,
                                          inv8_sqrt3,  -inv8_sqrt3, inv8_sqrt3,  -inv_sqrt3, -inv8_sqrt3, inv_sqrt3,
                                          -inv8_sqrt3, inv8_sqrt3,  -inv8_sqrt3, inv_sqrt3,  inv8_sqrt3,  -inv_sqrt3};
  const std::vector<double> expected_d = {d_leading,  d_leading,  d_leading,  d_leading, d_leading,  d_leading,
                                          -d_sublead, -d_sublead, -d_sublead, d_leading, -d_sublead, d_leading,
                                          -d_sublead, -d_sublead, -d_sublead, d_leading, -d_sublead, d_leading,
                                          -d_sublead, -d_sublead, -d_sublead, d_leading, -d_sublead, d_leading};

  std::array<int, 4> tail = {{1, 2, 3, 4}};
  std::size_t        idx  = 0;
  do {
    const std::array<int, 5> perm    = {{0, tail[0], tail[1], tail[2], tail[3]}};
    std::complex<double>     coeff_f = 0.0;
    std::complex<double>     coeff_d = 0.0;

    for (int a = 0; a < 8; ++a) {
      for (int b = 0; b < 8; ++b) {
        for (int c = 0; c < 8; ++c) {
          std::array<int, 5> color = {{a, a, b, c, 0}};
          for (int d = 0; d < 8; ++d) {
            color[4]     = d;
            Matrix3 prod = T[color[perm[0]]];
            for (std::size_t k = 1; k < perm.size(); ++k) { prod = multiply(prod, T[color[perm[k]]]); }
            coeff_f += trace(prod) * fabc(b, c, d) / (8.0 * std::sqrt(24.0));
            coeff_d += trace(prod) * dabc(b, c, d) / (8.0 * std::sqrt(40.0 / 3.0));
          }
        }
      }
    }

    coeff_f *= -gra::math::zi;  // Rotate the antisymmetric singlet to the real
                                // convention used in MDurham
    CAPTURE(idx, perm, coeff_f, coeff_d);
    REQUIRE(coeff_f.real() == Approx(expected_f[idx]).epsilon(1e-12));
    REQUIRE(coeff_f.imag() == Approx(0.0).margin(1e-12));
    REQUIRE(coeff_d.real() == Approx(expected_d[idx]).epsilon(1e-12));
    REQUIRE(coeff_d.imag() == Approx(0.0).margin(1e-12));
    ++idx;
  } while (std::next_permutation(tail.begin(), tail.end()));
}

TEST_CASE(
    "MDurham projected gg->gg helicity amplitudes use the HELAS matrix "
    "element",
    "[gra::MDurham][helicity]") {
  ModelParamRestoreGuard model_restore;
  const auto             relaxed_tune = WriteRelaxedDurhamTune("gg_helas");
  gra::MODELPARAM                     = relaxed_tune.first;

  gra::LORENTZSCALAR lts = MakeToyDurhamGG();
  REQUIRE(lts.pfinal[0].M() > 2.0);
  REQUIRE(lts.t_hat < 0.0);
  REQUIRE(lts.u_hat < 0.0);

  gra::MRandom rng;
  rng.SetSeed(2);
  MDurham durham(lts, gra::MModelTune::Load(relaxed_tune.second), rng,
                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  MDurham::DurhamProjectedAmp projected(4);
  for (auto &channel : projected) { channel.fill(0.0); }
  std::vector<MDurham::DurhamProjectedAmp> color_amp;
  EvaluateGeneratedDurham(durham, lts, projected, &color_amp);

  bool has_nonzero = false;
  for (std::size_t h = 0; h < projected.size(); ++h) {
    for (std::size_t i = 0; i < projected[h].size(); ++i) {
      CAPTURE(h, i, projected[h][i]);
      CHECK(std::isfinite(projected[h][i].real()));
      CHECK(std::isfinite(projected[h][i].imag()));
      has_nonzero = has_nonzero || std::abs(projected[h][i]) > 1e-12;
    }
  }
  REQUIRE(has_nonzero);

  // A common beam-axis rotation leaves raw HELAS hard components unchanged
  for (const double angle : {-1.17, 0.43, 2.06}) {
    gra::LORENTZSCALAR          rotated_lts = RotateToyEventAroundZ(lts, angle);
    MDurham::DurhamProjectedAmp rotated_generated(4);
    for (auto &channel : rotated_generated) { channel.fill(0.0); }
    EvaluateGeneratedDurham(durham, rotated_lts, rotated_generated);
    for (std::size_t h = 0; h < projected.size(); ++h) {
      for (std::size_t i = 0; i < projected[h].size(); ++i) {
        CAPTURE(angle, h, i, rotated_generated[h][i], projected[h][i]);
        RequireComplexNear(rotated_generated[h][i], projected[h][i], 1.0e-10);
      }
    }
  }

  const double projected_amp2 = durham.DQtloop(lts, projected);
  REQUIRE(std::isfinite(projected_amp2));
  REQUIRE(projected_amp2 > 0.0);
  REQUIRE_FALSE(lts.proton_good_walker.has_value());

  gra::LORENTZSCALAR excited_lts = lts;
  excited_lts.excite1            = true;
  const double excited_amp2      = durham.DQtloop(excited_lts, projected);
  const auto   model_tune        = gra::MModelTune::Load(relaxed_tune.second);
  const auto   pomeron           = model_tune->Soft()->ForwardExcitationExchange();
  const double elastic_factor    = model_tune->Soft()->NormalizedPhysicalResidue(pomeron, lts.t1);
  const double excitation_factor =
      model_tune->Soft()->ForwardExcitationFactor(pomeron, excited_lts.t1, excited_lts.pfinal[1].M2());
  REQUIRE(excited_amp2 ==
          Approx(projected_amp2 * gra::math::pow2(excitation_factor / elastic_factor)).epsilon(1.0e-10));
  REQUIRE_FALSE(excited_lts.proton_good_walker.has_value());
  REQUIRE(color_amp.size() == 1);
  REQUIRE(color_amp[0].size() == projected.size());
  for (std::size_t h = 0; h < projected.size(); ++h) {
    for (std::size_t i = 0; i < projected[h].size(); ++i) { CHECK(color_amp[0][h][i] == projected[h][i]); }
  }
}

TEST_CASE("MDurham elastic legs use the normalized physical beam residue", "[gra::MDurham][gra::SoftModel][residue]") {
  // Write one two-channel tune with a controlled off-diagonal bare coupling
  const auto write_tune = [](const std::string &suffix, const double off_diagonal) {
    return WriteModifiedSudakovModelTune(
        suffix,
        [off_diagonal](auto &card) {
          card["PARAM_SOFT"]["active_model"]      = "double";
          auto &model                             = card["PARAM_SOFT"]["MODEL"]["double"];
          model["EXCHANGE"]["P"]["g"]             = {{4.0, off_diagonal}, {off_diagonal, 2.0}};
          model["EXCHANGE"]["P"]["transition_ff"] = "diagonal";
          card["PARAM_DURHAM"]["JET_pt_min"]   = 0.0;
          card["PARAM_DURHAM"]["JET_rap_max"]  = 100.0;
        },
        [](auto &numerics) {
          numerics["NUMERICS_DURHAM"]["N_qt"]  = 6;
          numerics["NUMERICS_DURHAM"]["N_phi"] = 4;
        });
  };
  const std::string diagonal_tune        = write_tune("durham_diagonal_residue", 0.0);
  const std::string off_diagonal_tune    = write_tune("durham_offdiagonal_residue", 1.0);
  const auto        diagonal_model       = gra::MModelTune::Load(diagonal_tune + "/GENERAL.json");
  const auto        off_diagonal_model   = gra::MModelTune::Load(off_diagonal_tune + "/GENERAL.json");
  const auto        pomeron              = diagonal_model->Soft()->ForwardExcitationExchange();
  const auto        off_diagonal_pomeron = off_diagonal_model->Soft()->ForwardExcitationExchange();
  REQUIRE(std::abs(diagonal_model->Soft()->PhysicalCoupling(pomeron) -
                   off_diagonal_model->Soft()->PhysicalCoupling(off_diagonal_pomeron)) > 1.0e-3);
  CHECK(diagonal_model->Soft()->PhysicalResidue(pomeron, 0.0) ==
        Approx(off_diagonal_model->Soft()->PhysicalResidue(off_diagonal_pomeron, 0.0)).margin(1.0e-14));

  gra::LORENTZSCALAR diagonal_lts     = MakeToyDurhamGG();
  gra::LORENTZSCALAR off_diagonal_lts = MakeToyDurhamGG();
  gra::MRandom       diagonal_rng;
  gra::MRandom       off_diagonal_rng;
  diagonal_rng.SetSeed(81);
  off_diagonal_rng.SetSeed(81);
  MDurham diagonal(diagonal_lts, diagonal_model, diagonal_rng,
                   gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  MDurham off_diagonal(off_diagonal_lts, off_diagonal_model, off_diagonal_rng,
                       gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  MDurham::DurhamProjectedAmp hard(2);
  hard[0] = {std::complex<double>(0.3, -0.2), std::complex<double>(-0.1, 0.4), std::complex<double>(0.5, 0.1),
             std::complex<double>(-0.2, -0.3)};
  hard[1] = {std::complex<double>(-0.4, 0.2), std::complex<double>(0.1, 0.3), std::complex<double>(0.2, -0.5),
             std::complex<double>(0.6, 0.1)};
  const auto diagonal_amp     = diagonal.DQtloopAmplitudes(diagonal_lts, hard);
  const auto off_diagonal_amp = off_diagonal.DQtloopAmplitudes(off_diagonal_lts, hard);
  REQUIRE(diagonal_amp.size() == off_diagonal_amp.size());
  for (std::size_t row = 0; row < diagonal_amp.size(); ++row) {
    CAPTURE(row);
    RequireComplexNear(off_diagonal_amp[row], diagonal_amp[row], 2.0e-12);
  }
}

// Check the generated massless quark amplitude and Durham helicity ordering
// analytically
//
TEST_CASE("MDurham MG5 gg->qqbar matches the SuperChic 2 helicity amplitudes",
          "[gra::MDurham][helicity][MG2GRA][literature]") {
  constexpr double mass             = 100.0;
  constexpr double costheta         = 0.31;
  constexpr double phi              = 0.57;
  constexpr double alpha_s          = 0.118;
  const double     incoming_average = 1.0 / std::sqrt(8.0);
  const double     sintheta         = std::sqrt(1.0 - costheta * costheta);
  const double     momentum         = mass / 2.0;

  gra::LORENTZSCALAR lts;
  lts.q1                              = gra::M4Vec(0.0, 0.0, momentum, momentum);
  lts.q2                              = gra::M4Vec(0.0, 0.0, -momentum, momentum);
  lts.decaytree.resize(2);
  lts.decaytree[0].p.pdg = 2;
  lts.decaytree[1].p.pdg = -2;
  lts.decaytree[0].p4    = gra::M4Vec(momentum * sintheta * std::cos(phi), momentum * sintheta * std::sin(phi),
                                      momentum * costheta, momentum);
  lts.decaytree[1].p4    = -lts.decaytree[0].p4;
  lts.decaytree[1].p4.SetE(momentum);

  const auto process = gra::amplitude::FindProcess("DURHAM", "gg_uubar");
  REQUIRE(process.has_value());
  auto matrix_element = gra::CreateDurhamMG5Process(*process);
  REQUIRE(matrix_element != nullptr);
  const auto evaluation = matrix_element->Evaluate(lts, alpha_s);
  REQUIRE(evaluation.Valid());
  REQUIRE(matrix_element->ColorRank() == 1);
  REQUIRE(evaluation.projected.size() == matrix_element->HelicityCount());
  auto mg5 = evaluation.projected;
  for (auto &value : mg5) { value *= incoming_average; }

  // [REFERENCE: Harland-Lang et al., arXiv:1508.02718v2, Eqs. (A.1)--(A.2)]
  const double                                       norm        = 4.0 * gra::math::PI * alpha_s / std::sqrt(3.0);
  const auto                                         phase_plus  = std::polar(1.0, 2.0 * phi);
  const auto                                         phase_minus = std::conj(phase_plus);
  std::array<std::array<std::complex<double>, 4>, 4> analytic{};
  analytic[2][2] = -norm * (1.0 - costheta) / sintheta * phase_plus;
  analytic[2][3] = norm * (1.0 + costheta) / sintheta * phase_minus;
  analytic[3][2] = norm * (1.0 + costheta) / sintheta * phase_plus;
  analytic[3][3] = -norm * (1.0 - costheta) / sintheta * phase_minus;

  double mg5_norm      = 0.0;
  double analytic_norm = 0.0;
  for (const auto value : mg5) { mg5_norm += std::norm(value); }
  for (const auto &final : analytic) {
    for (const auto value : final) { analytic_norm += std::norm(value); }
  }
  REQUIRE(mg5_norm == Approx(analytic_norm).epsilon(1e-12));

  const gra::M4Vec           transverse_q1(0.7, -0.4, 0.0, 0.0);
  const gra::M4Vec           transverse_q2(-0.2, 0.9, 0.0, 0.0);
  const gra::M4Vec           hard_k1(0.0, 0.0, momentum, momentum);
  const gra::M4Vec           hard_k2(0.0, 0.0, -momentum, momentum);
  const double               dot   = transverse_q1.Px() * transverse_q2.Px() + transverse_q1.Py() * transverse_q2.Py();
  const double               cross = transverse_q1.Px() * transverse_q2.Py() - transverse_q1.Py() * transverse_q2.Px();
  const std::complex<double> pzero = -0.5 * dot;
  const std::complex<double> mzero = -0.5 * gra::math::zi * cross;
  const std::complex<double> ptwo =
      0.5 * std::complex<double>(transverse_q1.Px() * transverse_q2.Px() - transverse_q1.Py() * transverse_q2.Py(),
                                 transverse_q1.Px() * transverse_q2.Py() + transverse_q1.Py() * transverse_q2.Px());
  const std::complex<double> mtwo = std::conj(ptwo);

  double mg5_projected      = 0.0;
  double analytic_projected = 0.0;
  for (std::size_t final = 0; final < 4; ++final) {
    std::array<std::complex<double>, 4> hard{};
    for (std::size_t initial = 0; initial < 4; ++initial) {
      hard[initial] = mg5[gra::spin::PairHelicityTransitionIndex(initial, final)];
    }
    mg5_projected +=
        std::norm(gra::mg5helas::ContractTransverseSources(hard, transverse_q1, transverse_q2, hard_k1, hard_k2));
    const auto &reference = analytic[final];
    analytic_projected += std::norm(pzero * (reference[0] + reference[1]) + mzero * (reference[0] - reference[1]) +
                                    ptwo * reference[3] + mtwo * reference[2]);
  }
  REQUIRE(mg5_projected == Approx(analytic_projected).epsilon(1e-12));
}

TEST_CASE("Durham MG5 flow partitions reconstruct exact projections",
          "[gra::MDurham][color][MG2GRA][registry][projection]") {
  auto require_registry_projection = [](const std::string &process_name, const gra::LORENTZSCALAR &seed) {
    const auto process = gra::amplitude::FindProcess("DURHAM", process_name);
    REQUIRE(process.has_value());
    auto generated = gra::CreateDurhamMG5Process(*process);
    REQUIRE(generated != nullptr);

    gra::LORENTZSCALAR generated_lts = seed;
    const auto         evaluation    = generated->Evaluate(generated_lts, 0.118);
    REQUIRE(evaluation.Valid());
    REQUIRE(evaluation.projected.size() == generated->ColorRank() * generated->HelicityCount());
    REQUIRE(generated->ExactProjectors().size() == generated->ColorRank() * generated->ColorCount());
    REQUIRE(evaluation.flow_projected.size() == generated->FlowCandidates().size());

    for (std::size_t rank = 0; rank < generated->ColorRank(); ++rank) {
      for (std::size_t helicity = 0; helicity < generated->HelicityCount(); ++helicity) {
        const std::size_t    index    = rank * generated->HelicityCount() + helicity;
        std::complex<double> flow_sum = 0.0;
        for (const auto &flow : evaluation.flow_projected) {
          REQUIRE(flow.size() == evaluation.projected.size());
          flow_sum += flow[index];
        }
        CHECK(flow_sum.real() == Approx(evaluation.projected[index].real()).margin(5e-11));
        CHECK(flow_sum.imag() == Approx(evaluation.projected[index].imag()).margin(5e-11));
      }
    }
  };

  require_registry_projection("gg_gg", MakeToyDurhamGG());
  require_registry_projection("gg_uubar", MakeToyDurhamQQbar(2));
  require_registry_projection("gg_uubarg", MakeToyDurhamQQbarG(2));
  require_registry_projection("gg_ggg", MakeToyDurhamGGG());
}

TEST_CASE("Generated Durham registry evaluates the exact four-gluon color space", "[gra::MDurham][helicity][MG2GRA]") {
  ModelParamRestoreGuard model_restore;
  const std::vector<int> signature = {21, 21, 21, 21};
  const auto             process   = amplitude::FindProcess("DURHAM", "gg_gggg");
  REQUIRE(process.has_value());
  REQUIRE(HasDurhamMG5Process(*process));
  auto matrix_element = CreateDurhamMG5Process(*process);
  REQUIRE(matrix_element != nullptr);
  REQUIRE(matrix_element->Name() == "gg_gggg");
  REQUIRE(matrix_element->FinalPDGs() == signature);
  REQUIRE((matrix_element->FinalColorRepresentations() == std::vector<int>{8, 8, 8, 8}));
  REQUIRE(matrix_element->ColorCount() == 120);
  REQUIRE(matrix_element->ColorRank() == 8);
  REQUIRE(matrix_element->HelicityCount() == 64);
  REQUIRE(matrix_element->ExactProjectors().size() == 8 * 120);
  REQUIRE(matrix_element->FlowCandidates().size() == 9);

  gra::LORENTZSCALAR lts = MakeToyDurhamGGGG();
  gra::M4Vec         hard_k1;
  gra::M4Vec         hard_k2;
  const auto         evaluation = matrix_element->Evaluate(lts, 0.118, &hard_k1, &hard_k2);
  REQUIRE(evaluation.Valid());
  REQUIRE(evaluation.projected.size() == 8 * 64);
  REQUIRE(evaluation.flow_projected.size() == 9);
  REQUIRE(std::isfinite(hard_k1.E()));
  REQUIRE(std::isfinite(hard_k2.E()));

  std::size_t nonzero = 0;
  for (std::size_t i = 0; i < evaluation.projected.size(); ++i) {
    std::complex<double> flow_sum = 0.0;
    for (const auto &flow : evaluation.flow_projected) {
      REQUIRE(flow.size() == evaluation.projected.size());
      flow_sum += flow[i];
    }
    CAPTURE(i, evaluation.projected[i], flow_sum);
    CHECK(flow_sum.real() == Approx(evaluation.projected[i].real()).margin(5e-11));
    CHECK(flow_sum.imag() == Approx(evaluation.projected[i].imag()).margin(5e-11));
    if (std::abs(evaluation.projected[i]) > 1e-12) { ++nonzero; }
  }
  REQUIRE(nonzero > 0);

  const auto relaxed_tune = WriteRelaxedDurhamTune("gggg_registry");
  gra::MODELPARAM         = relaxed_tune.first;
  gra::MRandom rng;
  rng.SetSeed(31);
  MDurham                                       durham(lts, gra::MModelTune::Load(relaxed_tune.second), rng,
                                                       gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp              projected;
  std::vector<gra::MDurham::DurhamProjectedAmp> flow_projected;
  durham.Dgg2Generated(lts, *matrix_element, projected, &flow_projected);
  REQUIRE(projected.size() == 8 * 16);
  REQUIRE(flow_projected.size() == 9);
  for (const auto &flow : flow_projected) { REQUIRE(flow.size() == projected.size()); }
}

// Exact generated amplitudes must reject ungenerated resonance diagrams
TEST_CASE("Generated Durham amplitudes reject an unsupported nested resonance", "[gra::MDurham][MG2GRA][cascade]") {
  auto event = MakeToyDurhamQQbarG(2);
  const auto descriptor = gra::amplitude::FindProcess("DURHAM", "gg_uubarg");
  REQUIRE(descriptor.has_value());
  auto amplitude = gra::CreateDurhamMG5Process(*descriptor);
  REQUIRE(amplitude != nullptr);
  REQUIRE(amplitude->MatchProcess(event.decaytree).has_value());
  gra::MDecayBranch parent;
  parent.p = LoadedPDGTable().FindByPDG(25);
  parent.legs = {event.decaytree[0], event.decaytree[1]};
  parent.p4 = parent.legs[0].p4 + parent.legs[1].p4;
  event.decaytree = {parent, event.decaytree[2]};
  CHECK_FALSE(amplitude->MatchProcess(event.decaytree).has_value());
  const auto result = amplitude->Evaluate(event, 0.118);
  CHECK_FALSE(result.Valid());
}

TEST_CASE(
    "MDurham projected MadGraph matrix elements cover qqbar and ggg "
    "final states",
    "[gra::MDurham][helicity]") {
  ModelParamRestoreGuard model_restore;
  const auto             relaxed_tune = WriteRelaxedDurhamTune("matrix_element_coverage");
  gra::MODELPARAM                     = relaxed_tune.first;

  auto require_balanced_color_flow = [](const gra::LORENTZSCALAR &lts) {
    std::map<int, int> balance;
    for (const auto &branch : lts.decaytree) {
      REQUIRE_FALSE(branch.p.color_flow.empty());
      if (branch.p.color_flow.flow1 != 0) { ++balance[branch.p.color_flow.flow1]; }
      if (branch.p.color_flow.flow2 != 0) { --balance[branch.p.color_flow.flow2]; }
    }
    for (const auto &entry : balance) {
      CAPTURE(entry.first, entry.second);
      REQUIRE(entry.second == 0);
    }
  };

  auto finite_amp = [](const gra::MDurham::DurhamProjectedAmp &amp) {
    std::size_t nonzero = 0;
    for (const auto &channel : amp) {
      for (const auto &value : channel) {
        REQUIRE(std::isfinite(value.real()));
        REQUIRE(std::isfinite(value.imag()));
        if (std::abs(value) > 1e-12) { ++nonzero; }
      }
    }
    return nonzero;
  };

  gra::LORENTZSCALAR lts_gg = MakeToyDurhamGG();
  gra::MRandom       rng_gg;
  rng_gg.SetSeed(5);
  MDurham      durham_gg(lts_gg, gra::MModelTune::Load(relaxed_tune.second), rng_gg,
                         gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  const double gg_amp2 = durham_gg.DurhamQCD(lts_gg, "gg");
  REQUIRE(std::isfinite(gg_amp2));
  REQUIRE(gg_amp2 >= 0.0);
  REQUIRE_FALSE(lts_gg.proton_good_walker.has_value());
  REQUIRE(lts_gg.hamp.size() == 4);
  REQUIRE(lts_gg.hard_color_flows.size() == 1);
  REQUIRE(lts_gg.muF > 0.0);
  REQUIRE(lts_gg.muR > 0.0);
  REQUIRE(lts_gg.scalup == Approx(lts_gg.muF));
  durham_gg.SampleColorFlow(lts_gg);
  require_balanced_color_flow(lts_gg);

  gra::LORENTZSCALAR lts_qq = MakeToyDurhamQQbar(2);
  gra::MRandom       rng_qq;
  rng_qq.SetSeed(3);
  MDurham                          durham_qq(lts_qq, gra::MModelTune::Load(relaxed_tune.second), rng_qq,
                                             gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp qqbar_amp(4);
  for (auto &channel : qqbar_amp) { channel.fill(0.0); }
  EvaluateGeneratedDurham(durham_qq, lts_qq, qqbar_amp);
  REQUIRE(finite_amp(qqbar_amp) > 0);

  const double qqbar_amp2 = durham_qq.DurhamQCD(lts_qq, "MG5");
  REQUIRE(std::isfinite(qqbar_amp2));
  REQUIRE(qqbar_amp2 >= 0.0);
  REQUIRE_FALSE(lts_qq.proton_good_walker.has_value());
  REQUIRE(lts_qq.hamp.size() == 4);
  REQUIRE(lts_qq.hard_color_flows.size() == 1);
  REQUIRE(lts_qq.muF > 0.0);
  REQUIRE(lts_qq.muR > 0.0);
  REQUIRE(lts_qq.scalup == Approx(lts_qq.muF));
  durham_qq.SampleColorFlow(lts_qq);
  require_balanced_color_flow(lts_qq);

  gra::LORENTZSCALAR lts_ggg = MakeToyDurhamGGG();
  gra::MRandom       rng_ggg;
  rng_ggg.SetSeed(4);
  MDurham                          durham_ggg(lts_ggg, gra::MModelTune::Load(relaxed_tune.second), rng_ggg,
                                              gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp ggg_amp(16);
  for (auto &channel : ggg_amp) { channel.fill(0.0); }
  EvaluateGeneratedDurham(durham_ggg, lts_ggg, ggg_amp);
  REQUIRE(finite_amp(ggg_amp) > 0);

  const double ggg_amp2 = durham_ggg.DurhamQCD(lts_ggg, "MG5");
  REQUIRE(std::isfinite(ggg_amp2));
  REQUIRE(ggg_amp2 >= 0.0);
  REQUIRE_FALSE(lts_ggg.proton_good_walker.has_value());
  REQUIRE(lts_ggg.hamp.size() == 16);
  REQUIRE(lts_ggg.hard_color_flows.size() == 2);
  REQUIRE(lts_ggg.muF > 0.0);
  REQUIRE(lts_ggg.muR > 0.0);
  REQUIRE(lts_ggg.scalup == Approx(lts_ggg.muF));
  REQUIRE(lts_ggg.hard_color_flows[0].external.size() == 3);
  REQUIRE(lts_ggg.hard_color_flows[1].external.size() == 3);
  bool distinct_orientation = false;
  for (std::size_t i = 0; i < 3; ++i) {
    distinct_orientation =
        distinct_orientation ||
        lts_ggg.hard_color_flows[0].external[i].color != lts_ggg.hard_color_flows[1].external[i].color ||
        lts_ggg.hard_color_flows[0].external[i].anticolor != lts_ggg.hard_color_flows[1].external[i].anticolor;
  }
  REQUIRE(distinct_orientation);
  durham_ggg.SampleColorFlow(lts_ggg);
  require_balanced_color_flow(lts_ggg);

  const double valid_shat    = lts_gg.s_hat;
  double       rejected_amp2 = std::numeric_limits<double>::quiet_NaN();
  lts_gg.s_hat               = -1.0;
  REQUIRE_NOTHROW(rejected_amp2 = durham_gg.DurhamQCD(lts_gg, "gg"));
  REQUIRE(rejected_amp2 == Approx(0.0));
  REQUIRE(lts_gg.hamp.size() == 4);
  for (const auto &value : lts_gg.hamp) { REQUIRE(std::abs(value) == Approx(0.0)); }

  lts_gg.s_hat = valid_shat;
  lts_gg.GlobalSudakovPtr.reset();
  REQUIRE_NOTHROW(rejected_amp2 = durham_gg.DurhamQCD(lts_gg, "gg"));
  REQUIRE(rejected_amp2 == Approx(0.0));
}

TEST_CASE("MDurham supports explicit five-flavour qqbar and qqbar-g states", "[gra::MDurham][helicity][flavour]") {
  ModelParamRestoreGuard model_restore;
  const auto             relaxed_tune = WriteRelaxedDurhamTune("five_flavour");
  gra::MODELPARAM                     = relaxed_tune.first;

  gra::LORENTZSCALAR seed_lts = MakeToyDurhamGG();
  gra::MRandom       rng;
  rng.SetSeed(19);
  MDurham durham(seed_lts, gra::MModelTune::Load(relaxed_tune.second), rng,
                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  auto projected_norm = [](const gra::MDurham::DurhamProjectedAmp &amp) {
    double norm = 0.0;
    for (const auto &channel : amp) {
      for (const auto &value : channel) {
        REQUIRE(std::isfinite(value.real()));
        REQUIRE(std::isfinite(value.imag()));
        norm += gra::math::abs2(value);
      }
    }
    return norm;
  };

  std::vector<double> qqbar_norms;
  std::vector<double> qqbarg_norms;
  for (const int flavour : {1, 2, 3, 4, 5}) {
    gra::LORENTZSCALAR qqbar = flavour == 5
        ? DurhamBottomPoint(*gra::amplitude::FindProcess("DURHAM", "gg_bbbar"), false)
        : MakeToyDurhamQQbar(flavour);
    qqbar.GlobalSudakovPtr   = seed_lts.GlobalSudakovPtr;
    gra::MDurham::DurhamProjectedAmp qqbar_amp(4);
    EvaluateGeneratedDurham(durham, qqbar, qqbar_amp);
    CAPTURE(flavour);
    REQUIRE(qqbar_amp.size() == 4);
    qqbar_norms.push_back(projected_norm(qqbar_amp));
    REQUIRE(qqbar_norms.back() > 0.0);

    gra::LORENTZSCALAR qqbarg = flavour == 5
        ? DurhamBottomPoint(*gra::amplitude::FindProcess("DURHAM", "gg_bbbarg"), true)
        : MakeToyDurhamQQbarG(flavour);
    qqbarg.GlobalSudakovPtr   = seed_lts.GlobalSudakovPtr;
    gra::MDurham::DurhamProjectedAmp qqbarg_amp(8);
    EvaluateGeneratedDurham(durham, qqbarg, qqbarg_amp);
    REQUIRE(qqbarg_amp.size() == 8);
    qqbarg_norms.push_back(projected_norm(qqbarg_amp));
    REQUIRE(qqbarg_norms.back() > 0.0);
  }
  for (std::size_t flavour = 1; flavour < 4; ++flavour) {
    REQUIRE(qqbar_norms[flavour] == Approx(qqbar_norms[0]).epsilon(1e-12));
    REQUIRE(qqbarg_norms[flavour] == Approx(qqbarg_norms[0]).epsilon(1e-12));
  }

  auto require_balanced_color_flow = [](const gra::LORENTZSCALAR &lts) {
    std::map<int, int> balance;
    for (const auto &branch : lts.decaytree) {
      if (branch.p.color_flow.flow1 != 0) { ++balance[branch.p.color_flow.flow1]; }
      if (branch.p.color_flow.flow2 != 0) { --balance[branch.p.color_flow.flow2]; }
    }
    for (const auto &[tag, value] : balance) {
      CAPTURE(tag, value);
      REQUIRE(value == 0);
    }
  };

  gra::LORENTZSCALAR explicit_qqbar = MakeToyDurhamQQbar(4);
  explicit_qqbar.model_cache        = seed_lts.model_cache;
  explicit_qqbar.GlobalSudakovPtr   = seed_lts.GlobalSudakovPtr;
  REQUIRE(durham.DurhamQCD(explicit_qqbar, "MG5") > 0.0);
  REQUIRE(explicit_qqbar.hamp.size() == 4);
  durham.SampleColorFlow(explicit_qqbar);
  require_balanced_color_flow(explicit_qqbar);

  gra::LORENTZSCALAR explicit_qqbarg =
      DurhamBottomPoint(*gra::amplitude::FindProcess("DURHAM", "gg_bbbarg"), true);
  explicit_qqbarg.model_cache        = seed_lts.model_cache;
  gra::MDurham bottom(explicit_qqbarg, gra::MModelTune::Load(relaxed_tune.second), rng,
                      gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  REQUIRE(bottom.DurhamQCD(explicit_qqbarg, "MG5") > 0.0);
  REQUIRE(explicit_qqbarg.hamp.size() == 8);
  bottom.SampleColorFlow(explicit_qqbarg);
  REQUIRE(explicit_qqbarg.decaytree[0].p.color_flow.flow1 == 501);
  REQUIRE(explicit_qqbarg.decaytree[2].p.color_flow.flow2 == 501);
  REQUIRE(explicit_qqbarg.decaytree[2].p.color_flow.flow1 == 502);
  REQUIRE(explicit_qqbarg.decaytree[1].p.color_flow.flow2 == 502);
  require_balanced_color_flow(explicit_qqbarg);
}

TEST_CASE("MDurham applies resolved jet cuts before MadGraph amplitudes", "[gra::MDurham][jets]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  auto amp_norm = [](const gra::MDurham::DurhamProjectedAmp &amp) {
    double sum = 0.0;
    for (const auto &channel : amp) {
      for (const auto &value : channel) { sum += gra::math::abs2(value); }
    }
    return sum;
  };

  for (const std::string &algo : std::vector<std::string>{"anti-kt", "kt", "CA"}) {
    const auto         tune = WriteModifiedDurhamTune("jetcut_model_" + algo, [&algo](auto &general) {
      general["PARAM_DURHAM"]["JET_ALGO"]    = algo;
      general["PARAM_DURHAM"]["JET_R"]       = 10.0;
      general["PARAM_DURHAM"]["JET_pt_min"]  = 0.0;
      general["PARAM_DURHAM"]["JET_rap_max"] = 100.0;
    });
    const std::string &path = tune.second;

    gra::LORENTZSCALAR lts_gg = MakeToyDurhamGG();
    gra::MRandom       rng_gg;
    rng_gg.SetSeed(11);
    MDurham                          durham_gg(lts_gg, gra::MModelTune::Load(path), rng_gg,
                                               gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham::DurhamProjectedAmp gg_amp(4);
    for (auto &channel : gg_amp) { channel.fill(1.0); }
    std::vector<gra::MDurham::DurhamProjectedAmp> gg_color(1, gg_amp);
    EvaluateGeneratedDurham(durham_gg, lts_gg, gg_amp, &gg_color);
    REQUIRE(amp_norm(gg_amp) == Approx(0.0));
    REQUIRE(gg_color.empty());

    gra::LORENTZSCALAR lts_qq = MakeToyDurhamQQbar(2);
    gra::MRandom       rng_qq;
    rng_qq.SetSeed(12);
    MDurham                          durham_qq(lts_qq, gra::MModelTune::Load(path), rng_qq,
                                               gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham::DurhamProjectedAmp qq_amp(4);
    for (auto &channel : qq_amp) { channel.fill(1.0); }
    std::vector<gra::MDurham::DurhamProjectedAmp> qq_color(1, qq_amp);
    EvaluateGeneratedDurham(durham_qq, lts_qq, qq_amp, &qq_color);
    REQUIRE(amp_norm(qq_amp) == Approx(0.0));
    REQUIRE(qq_color.empty());

    gra::LORENTZSCALAR lts_qqg = MakeToyDurhamQQbarG(2);
    gra::MRandom       rng_qqg;
    rng_qqg.SetSeed(14);
    MDurham                          durham_qqg(lts_qqg, gra::MModelTune::Load(path), rng_qqg,
                                                gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham::DurhamProjectedAmp qqg_amp(8);
    for (auto &channel : qqg_amp) { channel.fill(1.0); }
    std::vector<gra::MDurham::DurhamProjectedAmp> qqg_color(1, qqg_amp);
    EvaluateGeneratedDurham(durham_qqg, lts_qqg, qqg_amp, &qqg_color);
    REQUIRE(amp_norm(qqg_amp) == Approx(0.0));
    REQUIRE(qqg_color.empty());

    gra::LORENTZSCALAR lts_ggg = MakeToyDurhamGGG();
    gra::MRandom       rng_ggg;
    rng_ggg.SetSeed(13);
    MDurham                          durham_ggg(lts_ggg, gra::MModelTune::Load(path), rng_ggg,
                                                gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    gra::MDurham::DurhamProjectedAmp ggg_amp(16);
    for (auto &channel : ggg_amp) { channel.fill(1.0); }
    std::vector<gra::MDurham::DurhamProjectedAmp> ggg_color(1, ggg_amp);
    EvaluateGeneratedDurham(durham_ggg, lts_ggg, ggg_amp, &ggg_color);
    REQUIRE(amp_norm(ggg_amp) == Approx(0.0));
    REQUIRE(ggg_color.empty());
  }

  const auto         none_tune = WriteModifiedDurhamTune("jetcut_model_none", [](auto &general) {
    general["PARAM_DURHAM"]["JET_ALGO"]    = "none";
    general["PARAM_DURHAM"]["JET_R"]       = 0.0;
    general["PARAM_DURHAM"]["JET_pt_min"]  = 1.0e9;
    general["PARAM_DURHAM"]["JET_rap_max"] = 0.0;
  });
  const std::string &path      = none_tune.second;

  gra::LORENTZSCALAR lts_gg = MakeToyDurhamGG();
  gra::MRandom       rng_gg;
  rng_gg.SetSeed(15);
  MDurham                                       durham_gg(lts_gg, gra::MModelTune::Load(path), rng_gg,
                                                          gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp              gg_amp(4);
  std::vector<gra::MDurham::DurhamProjectedAmp> gg_color;
  EvaluateGeneratedDurham(durham_gg, lts_gg, gg_amp, &gg_color);
  REQUIRE(amp_norm(gg_amp) > 0.0);
  REQUIRE_FALSE(gg_color.empty());

  // TopologyMode the SuperChic three-parton resolved anti-kT decision at R =
  // 0.6 [REFERENCE: SuperChic src/user/cuts.f, three-body jflag selection]
  const auto         anti_tune = WriteModifiedDurhamTune("jetcut_model_superchic_antikt", [](auto &general) {
    general["PARAM_DURHAM"]["JET_ALGO"]    = "anti-kt";
    general["PARAM_DURHAM"]["JET_R"]       = 0.6;
    general["PARAM_DURHAM"]["JET_pt_min"]  = 0.0;
    general["PARAM_DURHAM"]["JET_rap_max"] = 100.0;
  });
  const std::string &anti_path = anti_tune.second;

  gra::LORENTZSCALAR separated_lts = MakeToyDurhamGGGPair(0.9);
  gra::MRandom       separated_rng;
  separated_rng.SetSeed(16);
  MDurham                          separated_durham(separated_lts, gra::MModelTune::Load(anti_path), separated_rng,
                                                    gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp separated_amp(16);
  std::vector<gra::MDurham::DurhamProjectedAmp> separated_color;
  EvaluateGeneratedDurham(separated_durham, separated_lts, separated_amp, &separated_color);
  REQUIRE(amp_norm(separated_amp) > 0.0);
  REQUIRE_FALSE(separated_color.empty());

  gra::LORENTZSCALAR merged_lts = MakeToyDurhamGGGPair(0.3);
  gra::MRandom       merged_rng;
  merged_rng.SetSeed(17);
  MDurham                                       merged_durham(merged_lts, gra::MModelTune::Load(anti_path), merged_rng,
                                                              gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp              merged_amp(16);
  std::vector<gra::MDurham::DurhamProjectedAmp> merged_color(1, merged_amp);
  EvaluateGeneratedDurham(merged_durham, merged_lts, merged_amp, &merged_color);
  REQUIRE(amp_norm(merged_amp) == Approx(0.0));
  REQUIRE(merged_color.empty());
}

TEST_CASE("MDurham meson wave functions use explicit eta mixing parameters", "[gra::MDurham][meson]") {
  ModelParamRestoreGuard model_restore;
  gra::MODELPARAM = "TUNE0";

  gra::LORENTZSCALAR lts = MakeToyDurhamMesonPair(221, 331);
  gra::MRandom       rng;
  rng.SetSeed(11);
  MDurham durham(lts, gra::MModelTune::Load(modelfile), rng,
                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  const std::vector<double> xval     = {0.2, 0.4, 0.6};
  const auto                phi_pi   = durham.EvalPhi(xval, 111);
  const auto                phi_eta  = durham.EvalPhi(xval, 221);
  const auto                phi_etap = durham.EvalPhi(xval, 331);

  REQUIRE(phi_pi.size() == xval.size());
  REQUIRE(phi_eta.size() == xval.size());
  REQUIRE(phi_etap.size() == xval.size());

  const double theta8         = gra::math::Deg2Rad(-21.2);
  const double theta1         = gra::math::Deg2Rad(-9.2);
  const double fpi            = 0.1300;
  const auto [nodes, weights] = gra::math::GaussLegendreRule(48, 0.0, 1.0);
  double phi_pi_integral      = 0.0;
  for (std::size_t i = 0; i < nodes.size(); ++i) { phi_pi_integral += weights[i] * durham.phi_CZ(nodes[i], fpi); }
  // [REFERENCE: Harland-Lang et al., arXiv:1105.1626v2, Eqs. (3.7) and (3.9)]
  REQUIRE(phi_pi_integral == Approx(fpi / (2.0 * std::sqrt(3.0))).epsilon(1e-12));

  for (std::size_t i = 0; i < xval.size(); ++i) {
    const double phi8          = durham.phi_CZ(xval[i], fpi * 1.26);
    const double phi0          = durham.phi_CZ(xval[i], fpi * 1.17);
    const double eta_expected  = std::cos(theta8) * phi8 - std::sin(theta1) * phi0;
    const double etap_expected = std::sin(theta8) * phi8 + std::cos(theta1) * phi0;

    REQUIRE(phi_eta[i] == Approx(eta_expected).epsilon(1e-12));
    REQUIRE(phi_etap[i] == Approx(etap_expected).epsilon(1e-12));
  }
  REQUIRE_THROWS(durham.EvalPhi(xval, 999999));
}

TEST_CASE("MDurham meson-pair amplitude accepts mixed eta eta-prime singlets", "[gra::MDurham][meson]") {
  ModelParamRestoreGuard model_restore;
  gra::MODELPARAM = "TUNE0";

  auto finite_amp = [](const gra::MDurham::DurhamProjectedAmp &amp) {
    std::size_t nonzero = 0;
    for (const auto &channel : amp) {
      for (const auto &value : channel) {
        REQUIRE(std::isfinite(value.real()));
        REQUIRE(std::isfinite(value.imag()));
        if (std::abs(value) > 1e-12) { ++nonzero; }
      }
    }
    return nonzero;
  };

  for (const auto &pair : {std::pair<int, int>{221, 331}, std::pair<int, int>{331, 221}}) {
    gra::LORENTZSCALAR lts = MakeToyDurhamMesonPair(pair.first, pair.second);
    gra::MRandom       rng;
    rng.SetSeed(12);
    MDurham durham(lts, gra::MModelTune::Load(modelfile), rng,
                   gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

    gra::MDurham::DurhamProjectedAmp amp(1);
    amp[0].fill(0.0);
    durham.Dgg2MMbar(lts, amp);
    REQUIRE(finite_amp(amp) > 0);

    const double amp2 = durham.DurhamQCD(lts, "MMbar");
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 >= 0.0);
  }

  const auto definition = gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::MesonPair, "MMbar");
  const auto mixed_octet_singlet = MakeToyDurhamMesonPair(111, 221);
  REQUIRE_FALSE(definition->MatchProcess(mixed_octet_singlet.decaytree).has_value());
  auto cascade = MakeToyDurhamMesonPair(221, 221);
  cascade.decaytree[0].legs.emplace_back();
  REQUIRE_FALSE(definition->MatchProcess(cascade.decaytree).has_value());

  gra::LORENTZSCALAR low_pt = MakeToyDurhamMesonPair(221, 221, 0.8);
  gra::MRandom       low_pt_rng;
  low_pt_rng.SetSeed(15);
  MDurham                          low_pt_durham(low_pt, gra::MModelTune::Load(modelfile), low_pt_rng,
                                                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  gra::MDurham::DurhamProjectedAmp low_pt_amp(1);
  low_pt_amp[0].fill(1.0);
  low_pt_durham.Dgg2MMbar(low_pt, low_pt_amp);
  for (const auto &value : low_pt_amp[0]) { REQUIRE(std::abs(value) == Approx(0.0)); }
}

// Check charmonium production, invalid loop kinematics and physical decay normalization
TEST_CASE("MDurham charmonium amplitudes cover chi_c0 chi_c1 and chi_c2", "[gra::MDurham][chic]") {
  ModelParamRestoreGuard model_restore;
  gra::MODELPARAM = "TUNE0";

  auto require_finite_chic = [](const std::string &process, std::size_t channels) {
    gra::LORENTZSCALAR lts = MakeToyDurhamGG();
    lts.PDG                = LoadedPDGTable();
    gra::MRandom rng;
    rng.SetSeed(17);
    MDurham      durham(lts, gra::MModelTune::Load(modelfile), rng,
                        gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    const double amp2 = durham.DurhamQCD(lts, process);
    REQUIRE(std::isfinite(amp2));
    REQUIRE(amp2 >= 0.0);
    REQUIRE(lts.hamp.size() == channels);
    for (const auto &value : lts.hamp) {
      REQUIRE(std::isfinite(value.real()));
      REQUIRE(std::isfinite(value.imag()));
    }
  };

  require_finite_chic("chic(0)", 1);
  require_finite_chic("chic(1)", 3);
  require_finite_chic("chic(2)", 5);

  SECTION("invalid event and loop kinematics have zero weight") {
    gra::LORENTZSCALAR invalid_lts = MakeToyDurhamGG();
    invalid_lts.PDG                = LoadedPDGTable();
    gra::MRandom invalid_rng;
    invalid_rng.SetSeed(16);
    MDurham invalid_durham(invalid_lts, gra::MModelTune::Load(modelfile), invalid_rng,
                           gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

    // Prepare physical spin bases before testing invalid loop momenta
    const auto kernel0 = invalid_durham.Dgg2chic0(invalid_lts);
    const auto kernel1 = invalid_durham.Dgg2chic1(invalid_lts);
    const auto kernel2 = invalid_durham.Dgg2chic2(invalid_lts);
    invalid_lts.pfinal.clear();
    REQUIRE(gra::math::IsZero(invalid_durham.DurhamQCD(invalid_lts, "chic(2)")));
    REQUIRE(invalid_lts.hamp == std::vector<std::complex<double>>(5, 0.0));

    const double nan = std::numeric_limits<double>::quiet_NaN();
    REQUIRE(kernel0({nan, 0.0}, {0.1, -0.3}, {}) ==
            std::vector<std::complex<double>>(1, 0.0));
    REQUIRE(kernel1({0.2, 0.0}, {0.1, nan}, {}) ==
            std::vector<std::complex<double>>(3, 0.0));
    REQUIRE(kernel2({nan, 0.0}, {0.1, -0.3}, {}) ==
            std::vector<std::complex<double>>(5, 0.0));
  }

  gra::LORENTZSCALAR lts = MakeToyDurhamGG();
  lts.PDG                = LoadedPDGTable();
  gra::MRandom rng;
  rng.SetSeed(18);
  MDurham    durham(lts, gra::MModelTune::Load(modelfile), rng,
                    gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  const auto pdg_kernel = durham.Dgg2chic0(lts);
  const auto pdg_amp    = pdg_kernel({0.21, -0.13}, {-0.17, 0.09}, {});

  lts.process.ROOT_RES_ACTIVE = true;
  lts.process.ROOT_RES.p      = lts.PDG.FindByPDG(10441);
  lts.process.ROOT_RES.p.mass += 0.05;
  lts.process.ROOT_RES.p.width *= 2.0;
  const auto root_kernel = durham.Dgg2chic0(lts);
  const auto root_amp    = root_kernel({0.21, -0.13}, {-0.17, 0.09}, {});
  REQUIRE(std::abs(root_amp.at(0) - pdg_amp.at(0)) > 1e-12);

  lts.process.ROOT_RES.p = lts.PDG.FindByPDG(20443);
  REQUIRE_THROWS(durham.Dgg2chic0(lts));

  gra::LORENTZSCALAR decay_lts      = MakeToyDurhamGG();
  decay_lts.PDG                     = LoadedPDGTable();
  decay_lts.process.ROOT_RES_ACTIVE = true;
  decay_lts.process.ROOT_RES.p      = decay_lts.PDG.FindByPDG(10441);
  decay_lts.process.root_decay_mode = gra::RootDecayMode::None;
  gra::MRandom decay_rng;
  decay_rng.SetSeed(19);
  MDurham      decay_durham(decay_lts, gra::MModelTune::Load(modelfile), decay_rng,
                            gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  const double production_amp2 = decay_durham.DurhamQCD(decay_lts, "chic(0)");
  REQUIRE(production_amp2 > 0.0);

  decay_lts.amplitude.BeginCentral();
  MMatrix<std::complex<double>> cached_decay(1, 1, 0.0);
  cached_decay[0][0] = 1.0;
  decay_lts.amplitude.durham_decay.Store(decay_lts.amplitude.central, "root", std::move(cached_decay),
                                         "test Durham root decay");
  decay_lts.process.ROOT_RES.hel_decay.g_decay = std::complex<double>(0.3, -0.4);
  decay_lts.process.root_decay_mode            = gra::RootDecayMode::Physical;
  decay_lts.screening.active                 = true;
  const double decay_amp2                      = decay_durham.DurhamQCD(decay_lts, "chic(0)");
  REQUIRE(decay_lts.process.ROOT_RES.decay_f.size_row() == 0);
  REQUIRE(
      decay_lts.amplitude.durham_decay.Get(decay_lts.amplitude.central, "root", "test Durham root decay").size_row() ==
      1);
  const double expected_decay_factor = std::norm(decay_lts.process.ROOT_RES.hel_decay.g_decay) * gra::math::PI /
                                       (decay_lts.process.ROOT_RES.p.mass * decay_lts.process.ROOT_RES.p.width);
  REQUIRE(decay_amp2 / production_amp2 == Approx(expected_decay_factor).epsilon(1e-11));
}

TEST_CASE("MDurham chi_cJ vertices follow chi_c0 two-gluon width normalization", "[gra::MDurham][chic][params]") {
  ModelParamRestoreGuard                  restore;
  const MDurham::DurhamTransverseMomentum q1 = {0.21, -0.13};
  const MDurham::DurhamTransverseMomentum q2 = {-0.17, 0.09};

  gra::MODELPARAM                = "TUNE0";
  gra::LORENTZSCALAR nominal_lts = MakeToyDurhamGG();
  nominal_lts.PDG                = LoadedPDGTable();
  gra::MRandom nominal_rng;
  nominal_rng.SetSeed(20);
  MDurham nominal(nominal_lts, gra::MModelTune::Load(modelfile), nominal_rng,
                  gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  const auto        nominal_general        = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  const double      nominal_width_fraction = nominal_general.at("PARAM_DURHAM").at("chic0_gg_width_fraction");
  const std::string scaled_tune            = WriteModifiedSudakovModelTune(
                 "chic_width_fraction",
                 [nominal_width_fraction](auto &j) {
        j["PARAM_DURHAM"]["chic0_gg_width_fraction"] = nominal_width_fraction / 4.0;
      },
                 [](auto &) {});
  gra::MODELPARAM               = scaled_tune;
  gra::LORENTZSCALAR scaled_lts = MakeToyDurhamGG();
  scaled_lts.PDG                = LoadedPDGTable();
  gra::MRandom scaled_rng;
  scaled_rng.SetSeed(21);
  MDurham scaled(scaled_lts, gra::MModelTune::Load(scaled_tune + "/GENERAL.json"), scaled_rng,
                 gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));

  for (int spin = 0; spin <= 2; ++spin) {
    const auto nominal_kernel = spin == 0   ? nominal.Dgg2chic0(nominal_lts)
                                : spin == 1 ? nominal.Dgg2chic1(nominal_lts)
                                            : nominal.Dgg2chic2(nominal_lts);
    const auto scaled_kernel  = spin == 0   ? scaled.Dgg2chic0(scaled_lts)
                                : spin == 1 ? scaled.Dgg2chic1(scaled_lts)
                                            : scaled.Dgg2chic2(scaled_lts);
    const auto nominal_amp    = nominal_kernel(q1, q2, {});
    const auto scaled_amp     = scaled_kernel(q1, q2, {});

    REQUIRE(scaled_amp.size() == nominal_amp.size());
    std::size_t compared_channels = 0;
    for (std::size_t i = 0; i < nominal_amp.size(); ++i) {
      if (gra::math::IsZero(std::abs(nominal_amp[i]))) {
        REQUIRE(std::abs(scaled_amp[i]) == Approx(0.0));
        continue;
      }
      ++compared_channels;
      REQUIRE((scaled_amp[i] / nominal_amp[i]).real() == Approx(0.5).epsilon(1e-12));
      REQUIRE((scaled_amp[i] / nominal_amp[i]).imag() == Approx(0.0).margin(1e-12));
    }
    REQUIRE(compared_channels > 0);
  }
}

// Reject empty loop domains and scale choices with unsuppressed infrared poles
TEST_CASE("Durham cutoff validation covers every radial map and PDF scale", "[gra::MDurham][params][infrared]") {
  for (const bool logarithmic : {false, true}) {
    const auto empty = DurhamTestTune(
        "empty_cut_" + std::to_string(logarithmic), [](auto& j) { j["PARAM_DURHAM"]["loop_q2_cut"] = 4.0; },
        [=](auto& j) {
          j["NUMERICS_DURHAM"]["qt2_MAX"] = 4.0;
          j["NUMERICS_DURHAM"]["log_qt"]  = logarithmic;
        });
    REQUIRE_THROWS_AS(gra::ReadDurhamConfig(*empty), std::invalid_argument);
  }
  for (const std::string scale : {"MIN", "IN", "EX", "MAX", "AVG"}) {
    const auto tune = DurhamTestTune(
        "zero_cut_" + scale,
        [&](auto& j) {
          j["PARAM_DURHAM"]["PDF_scale"]   = scale;
          j["PARAM_DURHAM"]["loop_q2_cut"] = 0.0;
        },
        [](auto& j) { j["NUMERICS_DURHAM"]["log_qt"] = false; });
    if (scale == "MIN") {
      REQUIRE_NOTHROW(gra::ReadDurhamConfig(*tune));
    } else {
      REQUIRE_THROWS_AS(gra::ReadDurhamConfig(*tune), std::invalid_argument);
    }
  }
}

// Vary the physical factorization scale and test matching below Q0 without a mass veto
TEST_CASE("Durham physical scale variations reach the amplitude", "[gra::MDurham][scales][regression]") {
  std::vector<double> weights;
  for (const double scale : {0.25, 0.5, 1.0}) {
    const auto   tune = DurhamTestTune("factorization_" + std::to_string(scale), [=](auto& j) {
      j["PARAM_DURHAM"]["muF"]      = scale;
      j["PARAM_DURHAM"]["JET_ALGO"] = "none";
    });
    auto         lts  = MakeToyDurhamGG();
    gra::MRandom rng;
    gra::MDurham amp(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    const double value = amp.DurhamQCD(lts, "MG5");
    REQUIRE(value > 0.0);
    REQUIRE(lts.muF == Approx(scale * std::sqrt(lts.s_hat)));
    weights.push_back(value);
  }
  REQUIRE(weights[0] != Approx(weights[1]).epsilon(1e-3));
  REQUIRE(weights[1] != Approx(weights[2]).epsilon(1e-3));

  const auto tune = DurhamTestTune("small_muR", [](auto& j) {
    j["PARAM_DURHAM"]["muR"]      = 0.01;
    j["PARAM_DURHAM"]["JET_ALGO"] = "none";
  });
  auto       lts  = MakeToyDurhamGG();
  lts.PDG         = LoadedPDGTable();
  gra::MRandom rng;
  gra::MDurham amp(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  REQUIRE(amp.DurhamQCD(lts, "FLUX") > 0.0);
  for (const std::string process : {"chic(0)", "chic(1)", "chic(2)"}) { REQUIRE(amp.DurhamQCD(lts, process) > 0.0); }
  REQUIRE_THROWS_AS(amp.DurhamQCD(lts, "MG5"), gra::AmplitudeFailure);
}

// Keep a physical chi_c0 production point finite across the former matching mass cut
TEST_CASE("Durham low mass amplitudes remain finite through Q0", "[gra::MDurham][scales][chic]") {
  const auto tune = DurhamTestTune("matching_mass", [](auto& j) { j["PARAM_DURHAM"]["muF"] = 0.5; });
  auto       lts  = MakeToyDurhamGG();
  lts.PDG         = LoadedPDGTable();
  gra::MRandom rng;
  gra::MDurham amp(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
  const double boundary_mass = 2.0 * lts.GlobalSudakovPtr->GetQMin();
  for (const double ratio : {0.99, 0.9999, 1.0, 1.0001, 1.01}) {
    const double mass   = ratio * boundary_mass;
    const double energy = lts.pbeam1.E() - mass / 2.0;
    const double pz     = std::sqrt(energy * energy - 0.18 * 0.18 - 0.05 * 0.05 - gra::PDG::mp * gra::PDG::mp);
    lts.pfinal[1]       = gra::M4Vec(0.18, -0.05, pz, energy);
    lts.pfinal[2]       = gra::M4Vec(-0.18, 0.05, -pz, energy);
    lts.decaytree[0].p4 = gra::M4Vec(mass / 2.0, 0.0, 0.0, mass / 2.0);
    lts.decaytree[1].p4 = gra::M4Vec(-mass / 2.0, 0.0, 0.0, mass / 2.0);
    UpdateToyDurhamDerivedKinematics(lts);
    const double weight = amp.DurhamQCD(lts, "chic(0)");
    REQUIRE(std::isfinite(weight));
    REQUIRE(weight > 0.0);
  }
}

// Record invalid infrared members as amplitude failures while keeping the resolved flux usable
TEST_CASE("Durham infrared failures reach amplitude bookkeeping", "[gra::MDurham][infrared][failure]") {
  for (const double cut : {0.4, 4.0}) {
    const auto   tune = DurhamTestTune("infrared_failure_" + std::to_string(cut), [=](auto& j) {
      j["PARAM_SKEWED_UGD"]["mode"]        = "IR_RKHS";
      j["PARAM_SKEWED_UGD"]["rkhs_radius"] = 100.0;
      j["PARAM_DURHAM"]["loop_q2_cut"]     = cut;
    });
    auto         lts  = MakeToyDurhamGG();
    gra::MRandom rng;
    gra::MDurham amp(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
    if (cut < lts.GlobalSudakovPtr->GetQ2Min()) {
      REQUIRE_THROWS_AS(amp.DurhamQCD(lts, "FLUX"), gra::AmplitudeFailure);
      REQUIRE(amp.EvaluationStatus() == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
    } else {
      REQUIRE(amp.DurhamQCD(lts, "FLUX") > 0.0);
      REQUIRE(amp.EvaluationStatus() == gra::mg5helas::EvaluationStatus::Success);
    }
  }
}

// Compare coherent amplitudes against a denser loop with the default MIN prescription
TEST_CASE("Durham loop resolves suppressed helicities under refinement", "[gra::MDurham][convergence]") {
  const std::string scale = "MIN";
  std::map<std::string, MMatrix<std::complex<double>>> coarse;
  for (const bool refined : {false, true}) {
    const auto tune = DurhamTestTune(
        "convergence_" + scale + std::to_string(refined),
        [&](auto& j) {
          j["PARAM_DURHAM"]["PDF_scale"] = scale;
          j["PARAM_DURHAM"]["JET_ALGO"] = "none";
        },
        [=](auto& j) {
          if (refined) {
            j["NUMERICS_DURHAM"]["N_qt"] = 48;
            j["NUMERICS_DURHAM"]["N_phi"] = 24;
          }
        });
    const auto cache = std::make_shared<gra::MModelCache>(tune);
    for (const std::string process : {"MG5", "MMbar", "FLUX", "chic(0)", "chic(1)", "chic(2)"}) {
      auto lts = process == "MMbar" ? MakeToyDurhamMesonPair(211, -211) : MakeToyDurhamGG();
      lts.PDG = LoadedPDGTable();
      lts.model_cache = cache;
      if (process.starts_with("chic")) {
        const int pdg = process == "chic(0)" ? 10441 : process == "chic(1)" ? 20443 : 445;
        const double mass = lts.PDG.FindByPDG(pdg).mass;
        const double energy = lts.pbeam1.E() - mass / 2.0;
        const double pz = std::sqrt(energy * energy - 0.18 * 0.18 - 0.05 * 0.05 - gra::PDG::mp * gra::PDG::mp);
        lts.pfinal[1] = gra::M4Vec(0.18, -0.05, pz, energy);
        lts.pfinal[2] = gra::M4Vec(-0.18, 0.05, -pz, energy);
        lts.decaytree[0].p4 = gra::M4Vec(mass / 2.0, 0.0, 0.0, mass / 2.0);
        lts.decaytree[1].p4 = gra::M4Vec(-mass / 2.0, 0.0, 0.0, mass / 2.0);
        UpdateToyDurhamDerivedKinematics(lts);
      }
      gra::MRandom rng;
      gra::MDurham amp(lts, tune, rng, gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Generic, "test"));
      CAPTURE(scale, process, refined);
      REQUIRE(amp.DurhamQCD(lts, process) > 0.0);
      const MMatrix<std::complex<double>> current{lts.hamp};
      if (refined) {
        REQUIRE((current - coarse.at(process)).FrobNorm() / current.FrobNorm() < 0.01);
      } else {
        coarse.emplace(process, current);
      }
    }
  }
}
