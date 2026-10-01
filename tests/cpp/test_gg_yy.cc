// Analytic gluon fusion quark box physics tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <cmath>
#include <complex>
#include <filesystem>
#include <fstream>
#include <functional>
#include <future>
#include <catch.hpp>

#include "Graniitti/Amplitude/Photon/AMP_gg_yy.h"
#include "Graniitti/Amplitude/Photon/AMP_yy_yy.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Photon/MQED.h"

namespace {
using gra::aux::indices;

// Load a complete tune with modified Standard Model inputs through the real input reader
gra::MModelTunePtr BoxTune(const std::string &name, const std::function<void(nlohmann::json &)> &modify) {
  const std::filesystem::path source = gra::aux::ResolveProjectPath("modeldata/TUNE0");
  const std::filesystem::path target = source.parent_path().parent_path() / "tmp/param_sm_tests" / name;
  std::filesystem::create_directories(target);
  for (const auto &entry : std::filesystem::directory_iterator(source)) {
    if (entry.path().extension() == ".json") {
      std::filesystem::copy_file(entry.path(), target / entry.path().filename(),
                                 std::filesystem::copy_options::overwrite_existing);
    }
  }
  auto general = nlohmann::json::parse(gra::aux::GetInputData((source / "GENERAL.json").string()));
  modify(general);
  std::ofstream(target / "GENERAL.json") << general.dump(2);
  return gra::MModelTune::Load((target / "GENERAL.json").string());
}

// Construct an exactly closed massless partonic scattering point
// Apply the same rotation and boost to all four external particles
gra::LORENTZSCALAR BoxPoint(double energy, double cosine, double azimuth, bool transform) {
  gra::LORENTZSCALAR lts;
  const auto card = std::filesystem::path(__FILE__).parent_path().parent_path().parent_path() /
                    "modeldata/TUNE0/GENERAL.json";
  lts.model_cache = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(card.string()));
  lts.q1 = gra::M4Vec(0.0, 0.0, energy, energy);
  lts.q2 = gra::M4Vec(0.0, 0.0, -energy, energy);
  const double pt = energy * std::sqrt(1.0 - cosine * cosine);
  lts.decaytree.resize(2);
  lts.decaytree[0].p4 = gra::M4Vec(pt * std::cos(azimuth), pt * std::sin(azimuth), energy * cosine, energy);
  lts.decaytree[1].p4 = gra::M4Vec(-pt * std::cos(azimuth), -pt * std::sin(azimuth), -energy * cosine, energy);
  for (auto &branch : lts.decaytree) {
    branch.p.pdg = 22;
    branch.p.spinX2 = 2;
    branch.p.color = 1;
  }
  if (transform) {
    for (auto *p : {&lts.q1, &lts.q2, &lts.decaytree[0].p4, &lts.decaytree[1].p4}) {
      p->RotateY(0.37);
      p->RotateZ(-0.23);
      *p = p->LorentzBoost({0.13, -0.08, 0.21});
    }
  }
  lts.pfinal = {lts.decaytree[0].p4 + lts.decaytree[1].p4};
  lts.s_hat = lts.pfinal[0].M2();
  return lts;
}

}  // namespace

// Reject invalid Standard Model cards before any amplitude sampling
TEST_CASE("Standard Model loop inputs are validated during tune loading", "[gg_yy][input]") {
  REQUIRE_THROWS(BoxTune("missing_sm", [](auto &j) { j.erase("PARAM_SM"); }));
  REQUIRE_THROWS(BoxTune("missing_mass", [](auto &j) { j["PARAM_SM"]["mass"].erase("c"); }));
  REQUIRE_THROWS(BoxTune("massless", [](auto &j) { j["PARAM_SM"]["mass"]["u"] = 0.0; }));
  REQUIRE_THROWS(BoxTune("negative_mass", [](auto &j) { j["PARAM_SM"]["mass"]["w"] = -1.0; }));
  REQUIRE_THROWS(BoxTune("bad_alpha", [](auto &j) { j["PARAM_SM"]["alpha_em_inv"] = 0.0; }));
}

// Check shared mass steering reaches both loop kernels without mixing concurrent tunes
TEST_CASE("Both photon box amplitudes use independent Standard Model inputs", "[gg_yy][physics][thread]") {
  const auto nominal = BoxTune("nominal", [](auto &) {});
  const auto shifted = BoxTune("charm", [](auto &j) {
    auto &mass = j["PARAM_SM"]["mass"]["c"];
    mass = 2.0 * mass.template get<double>();
  });
  // Compute both helicity sums close to the charm threshold
  const auto evaluate = [](const gra::SMParam &sm) {
    const double s = 16.0;
    const double t = -0.35 * s;
    const double u = -s - t;
    gra::AMP_gg_yy gg(sm);
    const auto gluons = gg.Helicity(s, t, u, 0.118, 1.0 / sm.alpha_em_inv);
    std::array<std::complex<double>, 16> photons;
    if (!gra::lbyl::SMHelicityAmplitudes(s, t, u, 1.0 / sm.alpha_em_inv, sm, photons)) {
      throw std::runtime_error("massive photon box failed");
    }
    return std::array<double, 2>{gra::SquaredNorm(gluons), gra::SquaredNorm(photons)};
  };
  auto first = std::async(std::launch::async, evaluate, nominal->SM());
  auto second = std::async(std::launch::async, evaluate, shifted->SM());
  const auto base = first.get();
  const auto changed = second.get();
  const auto repeated = evaluate(nominal->SM());
  for (const auto &i : indices(base)) {
    REQUIRE(std::isfinite(base[i]));
    REQUIRE(base[i] > 0.0);
    REQUIRE(std::abs(changed[i] - base[i]) > 1.0e-4 * base[i]);
    REQUIRE(repeated[i] == Approx(base[i]).epsilon(1.0e-12));
  }
}

// Check ZERO, MG and LL coupling steering in both real photon box matrix elements
TEST_CASE("Photon boxes follow shared electromagnetic steering", "[gg_yy][physics]") {
  auto lts = BoxPoint(50.0, 0.31, 0.57, false);
  gra::AMP_gg_yy nominal(lts.model_cache->Tune().SM());
  const auto reference = nominal.Evaluate(lts, 0.118);
  REQUIRE(reference.Valid());
  const double norm = gra::SquaredNorm(reference.projected);
  gra::AMP_yy_yy photon_nominal(lts.model_cache->Tune().SM());
  const auto photon_reference = photon_nominal.Evaluate(lts, 0.0, false);
  REQUIRE(photon_reference.Valid());
  for (const std::string scheme : {"ZERO", "MG", "LL"}) {
    const auto tune = BoxTune("alpha_" + scheme, [&](auto &j) {
      j["PARAM_STRUCTURE"]["QED_alpha"] = scheme;
      auto &inverse = j["PARAM_SM"]["alpha_em_inv"];
      inverse = 2.0 * inverse.template get<double>();
    });
    lts.model_cache = std::make_shared<gra::MModelCache>(tune);
    gra::AMP_gg_yy amplitude(tune->SM());
    const auto result = amplitude.Evaluate(lts, 0.118);
    REQUIRE(result.Valid());
    const double alpha = gra::qed::alpha_QED((lts.q1 + lts.q2).M2(), scheme, 1.0 / tune->SM().alpha_em_inv);
    REQUIRE(gra::SquaredNorm(result.projected) ==
            Approx(norm * gra::math::pow2(alpha / gra::qed::alpha_QED())).epsilon(1.0e-11));
    gra::AMP_yy_yy photons(tune->SM());
    const auto photon_result = photons.Evaluate(lts, 0.0, false);
    REQUIRE(photon_result.Valid());
    REQUIRE(photon_result.amp2 ==
            Approx(photon_reference.amp2 * gra::math::pow4(alpha / gra::qed::alpha_QED())).epsilon(1.0e-11));
  }
}

// Test exact symmetries and numerical continuity across the quark thresholds
TEST_CASE("Analytic quark boxes remain covariant across mass thresholds", "[gg_yy][physics]") {
  const auto seed = BoxPoint(50.0, 0.3, 0.0, false);
  const auto &sm = seed.model_cache->Tune().SM();
  gra::AMP_gg_yy amplitude(sm);
  for (double energy : {0.2, 2.0, sm.c * 0.9999, sm.c * 1.0001, sm.b * 0.9999, sm.b * 1.0001,
                         sm.t * 0.9999, sm.t * 1.0001, 1000.0}) {
    auto lts = BoxPoint(energy, 0.31, 0.57, false);
    auto transformed = BoxPoint(energy, 0.31, 0.57, true);
    const auto direct = amplitude.Evaluate(lts, 0.118);
    const auto rotated = amplitude.Evaluate(transformed, 0.118);
    const auto scaled = amplitude.Evaluate(lts, 0.236);
    std::swap(lts.q1, lts.q2);
    const auto beams = amplitude.Evaluate(lts, 0.118);
    std::swap(lts.q1, lts.q2);
    std::swap(lts.decaytree[0], lts.decaytree[1]);
    const auto photons = amplitude.Evaluate(lts, 0.118);
    CAPTURE(energy);
    REQUIRE(direct.Valid());
    REQUIRE(rotated.Valid());
    REQUIRE(scaled.Valid());
    REQUIRE(beams.Valid());
    REQUIRE(photons.Valid());
    REQUIRE(gra::AllFinite(direct.projected));
    REQUIRE(gra::AllFinite(rotated.projected));
    const double norm = gra::SquaredNorm(direct.projected);
    REQUIRE(norm > 0.0);
    REQUIRE(gra::SquaredNorm(rotated.projected) == Approx(norm).epsilon(2.0e-9));
    REQUIRE(gra::SquaredNorm(scaled.projected) == Approx(4.0 * norm).epsilon(2.0e-12));
    for (const auto &row : indices(direct.projected)) {
      // Exchange identical gauge bosons and their ordered helicity labels
      const auto h1 = row / 8;
      const auto h2 = (row / 4) % 2;
      const auto h3 = (row / 2) % 2;
      const auto h4 = row % 2;
      const double bound = 2.0e-9 * std::abs(direct.projected[row]) + 1.0e-12;
      REQUIRE(std::abs(direct.projected[row] - beams.projected[8 * h2 + 4 * h1 + 2 * h3 + h4]) < bound);
      REQUIRE(std::abs(direct.projected[row] - photons.projected[8 * h1 + 4 * h2 + 2 * h4 + h3]) < bound);
      REQUIRE(std::norm(direct.projected[row]) == Approx(std::norm(direct.projected[15 - row])).epsilon(2.0e-9));
    }
  }
}
