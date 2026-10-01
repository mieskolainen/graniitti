// Physics tests for subnucleon incoherent photoproduction fluctuations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <json.hpp>
#include <limits>
#include <memory>
#include <stdexcept>
#include <vector>

#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MHotSpot.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MPhoto.h"
#include "Graniitti/Nuclear/MSteering.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "support/nuclear_test_support.hh"

using gra::aux::indices;

namespace {

// Construct a compact isoscalar target for hotspot physics tests
std::shared_ptr<const gra::nuclear::MNucleus> MakeOxygen() {
  gra::nuclear::NucleusParam param;
  param.pdg    = gra::nuclear::EncodeNuclearPDG(16, 8);
  param.mass   = 14.899;
  param.charge = {2.608, 0.513, 11.842, 64, 2048, 3.0, 4097, 1.0e-7};
  param.matter = param.charge;
  return std::make_shared<const gra::nuclear::MNucleus>(param);
}

// Construct one deterministic nested oxygen configuration bank
std::shared_ptr<const gra::nuclear::MConfigBank> MakeBank(const std::shared_ptr<const gra::nuclear::MNucleus> &nucleus,
                                                          const std::size_t count, const std::uint32_t seed = 71023) {
  gra::nuclear::ConfigParam param;
  param.count             = count;
  param.d_min             = 0.0;
  param.sweeps            = 0;
  param.max_trials        = 10000;
  param.neutron_nodes     = 4096;
  param.density_rel_tol   = 1.0e-6;
  param.negative_norm_tol = 1.0e-7;
  param.norm_tol          = 2.0e-5;
  const gra::nuclear::MConfigSampler sampler(*nucleus, param);
  gra::MRandom                       random;
  random.SetSeed(seed);
  return std::make_shared<const gra::nuclear::MConfigBank>(sampler, random);
}

// Compute compact optical quadratures for direct current tests
gra::nuclear::PhotoParam PhotoControls() {
  gra::nuclear::PhotoParam param;
  gra::test::Photo(param);
  param.b_nodes = 32;
  param.z_nodes = 32;
  return param;
}

// Compute the subnucleon parameters used for the ALICE incoherent model
gra::nuclear::HotSpotParam HotSpotControls() {
  gra::nuclear::HotSpotParam param;
  param.count          = 3;
  param.b_center       = 3.3;
  param.b_profile      = 0.7;
  param.strength_sigma = 0.5;
  return param;
}

// Sample one reproducible hotspot bank with an explicit test RNG
gra::nuclear::MHotSpot SampleHotSpot(const gra::nuclear::MConfigBank &bank, const gra::nuclear::HotSpotParam &param,
                                     const std::uint32_t seed) {
  gra::MRandom random;
  random.SetSeed(seed);
  return gra::nuclear::MHotSpot(bank, param, random);
}

// Store sampled complex-current moments
struct SampleMoment {
  std::complex<double> mean     = 0.0;
  double               second   = 0.0;
  double               variance = 0.0;
};

// Compute hotspot moments over every nucleon and nuclear configuration
SampleMoment HotSpotMoment(const gra::nuclear::MHotSpot &hotspot, const double qx, const double qy) {
  SampleMoment moment;
  const double count = static_cast<double>(hotspot.ConfigCount() * hotspot.NucleonCount());
  for (std::size_t sample = 0; sample < hotspot.ConfigCount(); ++sample) {
    for (std::size_t nucleon = 0; nucleon < hotspot.NucleonCount(); ++nucleon) {
      const auto value = hotspot.Factor(sample, nucleon, qx, qy);
      moment.mean += value;
      moment.second += std::norm(value);
    }
  }
  moment.mean /= count;
  moment.second /= count;
  moment.variance = std::max(0.0, moment.second - std::norm(moment.mean));
  return moment;
}

// Compute the analytic Gaussian-hotspot second moment
double HotSpotSecond(const gra::nuclear::HotSpotParam &param, const double qx, const double qy) {
  const double qt2      = qx * qx + qy * qy;
  const double diagonal = std::exp(param.strength_sigma * param.strength_sigma) / static_cast<double>(param.count);
  const double off_diagonal =
      (static_cast<double>(param.count) - 1.0) / static_cast<double>(param.count) * std::exp(-param.b_center * qt2);
  return std::exp(-param.b_profile * qt2) * (diagonal + off_diagonal);
}

// Compute the mean of one complex current bank
std::complex<double> CurrentMean(const std::vector<std::complex<double>> &current) {
  std::complex<double> mean = 0.0;
  for (const auto &value : current) { mean += value; }
  return mean / static_cast<double>(current.size());
}

}  // namespace

TEST_CASE("Gaussian hotspot currents converge to their analytic moments", "[gra::nuclear::MHotSpot][physics]") {
  const auto       nucleus = MakeOxygen();
  const auto       bank    = MakeBank(nucleus, 256);
  constexpr double qx      = 0.65;
  constexpr double qy      = -0.35;

  for (const std::uint32_t seed : {1907U, 2909U}) {
    const auto         param             = HotSpotControls();
    const auto         hotspot           = SampleHotSpot(*bank, param, seed);
    const SampleMoment moment            = HotSpotMoment(hotspot, qx, qy);
    const double       expected_mean     = hotspot.Mean(qx, qy);
    const double       expected_second   = HotSpotSecond(param, qx, qy);
    const double       expected_variance = expected_second - expected_mean * expected_mean;

    REQUIRE(moment.mean.real() == Approx(expected_mean).epsilon(0.06));
    REQUIRE(moment.mean.imag() == Approx(0.0).margin(0.02));
    REQUIRE(moment.second == Approx(expected_second).epsilon(0.06));
    REQUIRE(moment.variance == Approx(expected_variance).epsilon(0.08));
    REQUIRE(moment.variance > 0.0);
  }

  const auto         small_bank    = MakeBank(nucleus, 64);
  const auto         param         = HotSpotControls();
  const auto         small_hotspot = SampleHotSpot(*small_bank, param, 1907U);
  const SampleMoment small         = HotSpotMoment(small_hotspot, qx, qy);
  const double       expected_mean = small_hotspot.Mean(qx, qy);
  const double       expected      = HotSpotSecond(param, qx, qy) - expected_mean * expected_mean;
  REQUIRE(small.variance == Approx(expected).epsilon(0.18));
}

TEST_CASE("Hotspot target currents converge to the coherent Good-Walker moment",
          "[gra::nuclear::MPhoto][gra::nuclear::MHotSpot][physics]") {
  const auto nucleus       = MakeOxygen();
  const auto bank          = MakeBank(nucleus, 128);
  auto       base_param    = PhotoControls();
  auto       hotspot_param = base_param;
  hotspot_param.hotspot    = HotSpotControls();
  const gra::nuclear::MPhoto      base(nucleus, base_param);
  const gra::nuclear::MPhoto      hotspot(nucleus, hotspot_param);
  const auto                      hotspot_bank = SampleHotSpot(*bank, hotspot_param.hotspot, 92821U);
  const auto                      profile      = gra::test::PhotoProfile(12.0, 0.1);
  constexpr std::array<double, 3> q            = {0.43, -0.21, 0.08};

  const auto base_current =
      base.ShadowCurrents(profile, *bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ);
  const auto hotspot_current =
      hotspot.ShadowCurrents(profile, *bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ, &hotspot_bank);
  const auto   base_mean      = CurrentMean(base_current);
  const auto   hotspot_mean   = CurrentMean(hotspot_current);
  const double coherent_error = std::abs(hotspot_mean - base_mean) / std::max(1.0, std::abs(base_mean));
  REQUIRE(coherent_error < 0.15);

  double residual_second = 0.0;
  for (const auto &sample : indices(base_current)) {
    residual_second += std::norm(hotspot_current[sample] - base_current[sample]);
  }
  residual_second /= static_cast<double>(base_current.size());
  REQUIRE(residual_second > 0.0);

  const auto stat =
      hotspot.ShadowStat(profile, *bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ, &hotspot_bank);
  std::complex<double> direct_mean   = 0.0;
  double               direct_second = 0.0;
  for (const auto &value : hotspot_current) {
    direct_mean += value;
    direct_second += std::norm(value);
  }
  direct_mean /= static_cast<double>(hotspot_current.size());
  direct_second /= static_cast<double>(hotspot_current.size());
  const double direct_variance =
      gra::statistics::UnbiasedComplexVariance(direct_second, direct_mean, hotspot_current.size());
  REQUIRE(stat.mean.real() == Approx(base_mean.real()).margin(2.0e-13));
  REQUIRE(stat.mean.imag() == Approx(base_mean.imag()).margin(2.0e-13));
  REQUIRE(stat.second == Approx(std::norm(base_mean) + direct_variance).margin(2.0e-12));
  REQUIRE(stat.variance == Approx(direct_variance).margin(2.0e-12));
  REQUIRE(stat.variance >= 0.0);
  REQUIRE(stat.variance == Approx(stat.second - std::norm(stat.mean)).margin(1.0e-12));

  const auto moment_only =
      hotspot.Factors(profile, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ, bank.get(), &hotspot_bank);
  const auto stored_bank = hotspot.Factors(profile, q[0], q[1], q[2], gra::nuclear::PhotonDirection::NegativeZ,
                                           bank.get(), &hotspot_bank, true);
  REQUIRE(stored_bank.current.has_value());
  REQUIRE(moment_only.coherent.real() == Approx(stored_bank.coherent.real()).margin(2.0e-13));
  REQUIRE(moment_only.coherent.imag() == Approx(stored_bank.coherent.imag()).margin(2.0e-13));
  REQUIRE(moment_only.incoherent == Approx(stored_bank.incoherent).margin(2.0e-13));
}

TEST_CASE("Zero-width hotspots recover the nucleon-center target exactly",
          "[gra::nuclear::MPhoto][gra::nuclear::MHotSpot][limit]") {
  const auto nucleus        = MakeOxygen();
  const auto bank           = MakeBank(nucleus, 64);
  auto       base_param     = PhotoControls();
  auto       limit_param    = base_param;
  limit_param.hotspot.count = 1;
  const gra::nuclear::MPhoto      base(nucleus, base_param);
  const gra::nuclear::MPhoto      limit(nucleus, limit_param);
  const auto                      hotspot = SampleHotSpot(*bank, limit_param.hotspot, 331U);
  const auto                      profile = gra::test::PhotoProfile(12.0, 0.1);
  constexpr std::array<double, 3> q       = {0.37, 0.19, -0.06};
  const auto                      base_current =
      base.ShadowCurrents(profile, *bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::PositiveZ);
  const auto limit_current =
      limit.ShadowCurrents(profile, *bank, q[0], q[1], q[2], gra::nuclear::PhotonDirection::PositiveZ, &hotspot);
  for (const auto &sample : indices(base_current)) {
    REQUIRE(limit_current[sample].real() == Approx(base_current[sample].real()).margin(2.0e-13));
    REQUIRE(limit_current[sample].imag() == Approx(base_current[sample].imag()).margin(2.0e-13));
  }
}

TEST_CASE("Subnucleon geometry hardens the incoherent high-transfer tail",
          "[gra::nuclear::MPhoto][gra::nuclear::MHotSpot][physics]") {
  const auto nucleus                   = MakeOxygen();
  const auto bank                      = MakeBank(nucleus, 256);
  auto       base_param                = PhotoControls();
  auto       hotspot_param             = base_param;
  hotspot_param.hotspot                = HotSpotControls();
  hotspot_param.hotspot.strength_sigma = 0.0;
  const gra::nuclear::MPhoto base(nucleus, base_param);
  const gra::nuclear::MPhoto hotspot(nucleus, hotspot_param);
  const auto                 hotspot_bank = SampleHotSpot(*bank, hotspot_param.hotspot, 48163U);
  const auto                 profile      = gra::test::PhotoProfile(0.0, 0.1);

  const auto low_base = base.ShadowStat(profile, *bank, 0.08, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  const auto low_hotspot =
      hotspot.ShadowStat(profile, *bank, 0.08, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ, &hotspot_bank);
  const auto high_base = base.ShadowStat(profile, *bank, 0.80, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ);
  const auto high_hotspot =
      hotspot.ShadowStat(profile, *bank, 0.80, 0.0, 0.0, gra::nuclear::PhotonDirection::PositiveZ, &hotspot_bank);
  REQUIRE(low_base.variance > 0.0);
  REQUIRE(high_base.variance > 0.0);
  REQUIRE(high_hotspot.variance > high_base.variance);
  REQUIRE(high_hotspot.variance / high_base.variance > low_hotspot.variance / low_base.variance);
}

TEST_CASE("Lowercase hotspot steering reaches UPC target sectors",
          "[gra::nuclear::MSteering][gra::nuclear::MUPC][physics]") {
  nlohmann::json block = {{"emd", false},
                          {"fragmentation", nullptr},
                          {"emission", {"coherent", "coherent"}},
                          {"neutron_class", {"*", "*"}},
                          {"photoproduction", {{"target", {"coherent", "incoherent"}}, {"target_model", "impulse"}}},
                          {"sigma_NN", nullptr},
                          {"screening", "optical"},
                          {"structure", "hotspot"}};
  const auto     general =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  auto numerics =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json")));
  numerics["NUMERICS_NUCLEAR"]["SAMPLING"]["config_samples"]   = 64;
  numerics["NUMERICS_NUCLEAR"]["SAMPLING"]["current_samples"]  = 4;
  numerics["NUMERICS_NUCLEAR"]["SAMPLING"]["hard_core_sweeps"] = 0;
  const auto param =
      gra::nuclear::ReadUPCSteering(block, general.at("PARAM_NUCLEAR"), numerics.at("NUMERICS_NUCLEAR"), "TUNE0");
  REQUIRE(param.photo_model == gra::nuclear::PhotoModel::Impulse);
  REQUIRE(gra::nuclear::ParsePhotoModel("LTA") == gra::nuclear::PhotoModel::LTA);
  REQUIRE(gra::nuclear::PhotoModelName(gra::nuclear::PhotoModel::LTA) == "LTA");
  REQUIRE(param.structure == gra::nuclear::StructureType::Hotspot);
  REQUIRE(gra::nuclear::StructureName(param.structure) == "hotspot");
  REQUIRE(param.survival == gra::nuclear::SurvivalType::Optical);
  REQUIRE_FALSE(param.sigma_nn.has_value());
  REQUIRE_THROWS_AS(gra::nuclear::ParseSurvivalType("none"), std::invalid_argument);
  REQUIRE(gra::nuclear::ParseSurvivalType("mc_ggcf") == gra::nuclear::SurvivalType::MCGGCF);
  REQUIRE(gra::nuclear::SurvivalName(gra::nuclear::SurvivalType::MCGGCF) == "mc_ggcf");
  REQUIRE(param.current_count == 4);
  for (const auto& structure : {"nucleon", "hotspot"}) {
    auto configured         = block;
    configured["structure"] = structure;
    const auto parsed       = gra::nuclear::ReadUPCSteering(configured, general.at("PARAM_NUCLEAR"),
                                                            numerics.at("NUMERICS_NUCLEAR"), "TUNE0");
    REQUIRE(parsed.config.count == 64);
  }
  REQUIRE(param.photo[0].hotspot.count > 0);
  REQUIRE(param.photo[1].hotspot.count > 0);

  for (const auto& count :
       {nlohmann::json(-1), nlohmann::json(1.5), nlohmann::json("8"), nlohmann::json(4294967296ULL)}) {
    auto invalid                     = numerics.at("NUMERICS_NUCLEAR");
    invalid["GLAUBER"]["ggcf_nodes"] = count;
    REQUIRE_THROWS_AS(gra::nuclear::ReadUPCSteering(block, general.at("PARAM_NUCLEAR"), invalid, "TUNE0"),
                      std::invalid_argument);
  }
  auto invalid_samples                          = numerics.at("NUMERICS_NUCLEAR");
  invalid_samples["SAMPLING"]["config_samples"] = -1;
  REQUIRE_THROWS_AS(gra::nuclear::ReadUPCSteering(block, general.at("PARAM_NUCLEAR"), invalid_samples, "TUNE0"),
                    std::invalid_argument);

  auto reaction                                 = general.at("PARAM_NUCLEAR");
  auto invalid_decay                            = numerics.at("NUMERICS_NUCLEAR");
  invalid_decay["FRAGMENTATION"]["decay_steps"] = -1;
  REQUIRE_THROWS_AS(gra::nuclear::ReadUPCSteering(block, reaction, invalid_decay, "TUNE0"), std::invalid_argument);
  block["fragmentation"] = "internal";
  REQUIRE(gra::nuclear::ReadUPCSteering(block, reaction, numerics.at("NUMERICS_NUCLEAR"), "TUNE0").reaction ==
          gra::nuclear::ReactionType::Internal);
  block["fragmentation"] = "external";
  REQUIRE_THROWS_AS(gra::nuclear::ReadUPCSteering(block, reaction, numerics.at("NUMERICS_NUCLEAR"), "TUNE0"),
                    std::invalid_argument);
  reaction["FRAGMENTATION"]["external"]["library"] = "/test/libnuclear.so";
  REQUIRE(gra::nuclear::ReadUPCSteering(block, reaction, numerics.at("NUMERICS_NUCLEAR"), "TUNE0").reaction ==
          gra::nuclear::ReactionType::External);
  block["fragmentation"] = nullptr;
  REQUIRE_NOTHROW(gra::nuclear::ReadUPCSteering(block, reaction, numerics.at("NUMERICS_NUCLEAR"), "TUNE0"));
  reaction["EMD"]         = general.at("PARAM_NUCLEAR").at("EMD");
  auto obsolete           = block;
  obsolete["final_state"] = {{"library", "/test/libnuclear.so"}, {"config", nlohmann::json::object()}};
  REQUIRE_THROWS_AS(gra::nuclear::ReadUPCSteering(obsolete, reaction, numerics.at("NUMERICS_NUCLEAR"), "TUNE0"),
                    std::invalid_argument);

  auto override_block        = block;
  override_block["sigma_NN"] = 70.0;
  const auto override_param  = gra::nuclear::ReadUPCSteering(override_block, general.at("PARAM_NUCLEAR"),
                                                             numerics.at("NUMERICS_NUCLEAR"), "TUNE0");
  REQUIRE(override_param.sigma_nn.has_value());
  REQUIRE(*override_param.sigma_nn == Approx(70.0));

  const auto nucleus = MakeOxygen();
  REQUIRE_THROWS_AS(gra::nuclear::MUPC({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus},
                                       {nullptr, nucleus}, override_param),
                    std::invalid_argument);
  const gra::nuclear::MUPC model({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, nucleus},
                                 param);
  gra::MRandom             random;
  random.SetSeed(77123U);
  const auto upc = model.Sample(random);
  REQUIRE(upc->Photo(2)->Model() == gra::nuclear::PhotoModel::Impulse);
  const auto *bank = upc->Bank(2);
  REQUIRE(bank != nullptr);
  const gra::M3Vec q     = {0.52, -0.13, 0.07};
  const std::array ratio = {upc->TargetRatios(2, gra::nuclear::CoherenceType::Coherent, q),
                            upc->TargetRatios(2, gra::nuclear::CoherenceType::Incoherent, q)};
  REQUIRE(ratio[0].size() == bank->Size());
  REQUIRE(ratio[1].size() == bank->Size());
  REQUIRE(CurrentMean(ratio[0]).real() == Approx(1.0).margin(2.0e-12));
  REQUIRE(std::abs(CurrentMean(ratio[0]).imag()) < 2.0e-12);
  REQUIRE(std::abs(CurrentMean(ratio[1])) < 2.0e-12);
  double incoherent_second = 0.0;
  for (const auto &value : ratio[1]) { incoherent_second += std::norm(value); }
  incoherent_second /= static_cast<double>(ratio[1].size());
  REQUIRE(incoherent_second == Approx(1.0).margin(2.0e-12));

  const auto profile    = gra::test::PhotoProfile(12.0, 0.1);
  const auto transition = upc->Photo(2)->Factors(profile, q[0], q[1], q[2], gra::nuclear::TargetPhotonDirection(2),
                                                 bank, upc->HotSpot(2), true);
  REQUIRE(transition.current.has_value());
  const auto &cached_current = *transition.current;
  REQUIRE(cached_current.sample.size() == bank->Size());
  REQUIRE(transition.coherent.real() == Approx(cached_current.stat.mean.real()).margin(2.0e-12));
  REQUIRE(transition.coherent.imag() == Approx(cached_current.stat.mean.imag()).margin(2.0e-12));
  REQUIRE(transition.incoherent == Approx(std::sqrt(cached_current.stat.variance)).margin(2.0e-12));
  for (const auto target : {gra::nuclear::CoherenceType::Coherent, gra::nuclear::CoherenceType::Incoherent}) {
    const auto cached = upc->TargetCurrentRatios(2, target, cached_current);
    REQUIRE(cached.size() == bank->Size());
    REQUIRE(gra::AllFinite(cached));
  }

  gra::nuclear::ScreenLayout layout;
  layout.type  = gra::nuclear::ScreenType::Photo;
  layout.photo = gra::nuclear::PhotoChannels(*upc);
  gra::nuclear::ScreenPoint fallback_point;
  fallback_point.amplitude.assign(layout.photo.size(), {0.7, -0.2});
  fallback_point.transfer[1]    = q;
  auto cached_point             = fallback_point;
  cached_point.photo_current[1] = cached_current;
  const gra::nuclear::MUPCScreen fallback_screen(*upc, layout, fallback_point);
  const gra::nuclear::MUPCScreen cached_screen(*upc, layout, cached_point);
  const auto                     fallback_result = fallback_screen.Result();
  const auto                     cached_result   = cached_screen.Result();
  REQUIRE(cached_result.amplitude.size() == fallback_result.amplitude.size());
  REQUIRE(cached_result.helicity_norm.size() == fallback_result.helicity_norm.size());
  REQUIRE(gra::AllFinite(cached_result.amplitude));
  REQUIRE(gra::AllFinite(cached_result.helicity_norm));
  gra::nuclear::MUPCScreen shifted(*upc, layout, cached_point);
  gra::nuclear::LoopNode   node;
  node.weight = 0.2;
  REQUIRE_THROWS_AS(shifted.Add(node, fallback_point), gra::AmplitudeFailure);
  const auto unchanged = shifted.Result();
  for (const auto &h : indices(cached_result.helicity_norm)) {
    REQUIRE(unchanged.helicity_norm[h] == Approx(cached_result.helicity_norm[h]).margin(2.0e-13));
  }
  shifted.Add(node, cached_point);
  const auto accumulated = shifted.Result();
  for (const auto &h : indices(cached_result.helicity_norm)) {
    REQUIRE(accumulated.helicity_norm[h] == Approx(1.44 * cached_result.helicity_norm[h]).margin(2.0e-12));
  }


  auto invalid_current = cached_current;
  invalid_current.sample.pop_back();
  REQUIRE_THROWS_AS(upc->TargetCurrentRatios(2, gra::nuclear::CoherenceType::Incoherent, invalid_current),
                    std::invalid_argument);

  auto invalid        = block;
  invalid["screening"] = "config";
  REQUIRE_THROWS_AS(
      gra::nuclear::ReadUPCSteering(invalid, general.at("PARAM_NUCLEAR"), numerics.at("NUMERICS_NUCLEAR"), "TUNE0"),
      std::invalid_argument);

  invalid                                    = block;
  invalid["photoproduction"]["target_model"] = "unknown";
  REQUIRE_THROWS_AS(
      gra::nuclear::ReadUPCSteering(invalid, general.at("PARAM_NUCLEAR"), numerics.at("NUMERICS_NUCLEAR"), "TUNE0"),
      std::invalid_argument);

  for (const auto &value : {nlohmann::json(0.0), nlohmann::json(-1.0), nlohmann::json("70")}) {
    invalid             = block;
    invalid["sigma_NN"] = value;
    REQUIRE_THROWS_AS(
        gra::nuclear::ReadUPCSteering(invalid, general.at("PARAM_NUCLEAR"), numerics.at("NUMERICS_NUCLEAR"), "TUNE0"),
        std::invalid_argument);
  }
}

TEST_CASE("Photonuclear screening keeps distinct target profiles term local",
          "[gra::nuclear::MUPCScreen][photoproduction][screening][profile]") {
  nlohmann::json block = {{"emd", false},
                          {"fragmentation", nullptr},
                          {"emission", {"coherent", "coherent"}},
                          {"neutron_class", {"*", "*"}},
                          {"photoproduction", {{"target", {"coherent", "inclusive"}}, {"target_model", "glauber"}}},
                          {"sigma_NN", nullptr},
                          {"screening", "optical"},
                          {"structure", "nucleon"}};
  const auto     general =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  auto numerics =
      nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json")));
  numerics["NUMERICS_NUCLEAR"]["SAMPLING"]["current_samples"]  = 8;
  numerics["NUMERICS_NUCLEAR"]["SAMPLING"]["hard_core_sweeps"] = 0;
  const auto param =
      gra::nuclear::ReadUPCSteering(block, general.at("PARAM_NUCLEAR"), numerics.at("NUMERICS_NUCLEAR"), "TUNE0");
  const auto               nucleus = MakeOxygen();
  const gra::nuclear::MUPC model({gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Nucleus}, {nullptr, nucleus},
                                 param);
  gra::MRandom             random;
  random.SetSeed(88231U);
  const auto upc = model.Sample(random);
  REQUIRE(upc->HasSamples());
  REQUIRE(upc->Photo(2) != nullptr);
  REQUIRE(upc->Bank(2) != nullptr);

  gra::nuclear::ScreenLayout layout;
  layout.type  = gra::nuclear::ScreenType::Photo;
  layout.photo = gra::nuclear::PhotoChannels(*upc);
  REQUIRE(layout.photo.size() == 2);
  const auto selected = std::find_if(layout.photo.cbegin(), layout.photo.cend(), [](const auto &channel) {
    return channel.target == gra::nuclear::CoherenceType::Incoherent;
  });
  REQUIRE(selected != layout.photo.cend());
  const std::size_t channel = static_cast<std::size_t>(selected - layout.photo.cbegin());

  gra::nuclear::ScreenPoint point;
  point.transfer                                 = {gra::M3Vec{0.04, -0.03, 0.0}, gra::M3Vec{0.31, -0.12, 0.06}};
  const std::array<std::complex<double>, 2> hard = {std::complex<double>{1.0, 0.25}, std::complex<double>{-0.82, 0.11}};
  const std::array<gra::nuclear::PhotoProfile, 2> profile = {gra::test::PhotoProfile(5.0, 0.04),
                                                             gra::test::PhotoProfile(32.0, 0.24)};
  point.amplitude.assign(layout.photo.size(), 0.0);
  point.photo_terms.resize(hard.size());
  for (const auto &term : indices(hard)) {
    auto transition =
        upc->Photo(2)->Factors(profile[term], point.transfer[1][0], point.transfer[1][1], point.transfer[1][2],
                               gra::nuclear::TargetPhotonDirection(2), upc->Bank(2), upc->HotSpot(2), true);
    REQUIRE(transition.current.has_value());
    point.photo_terms[term].direction = gra::nuclear::PhotoDirection::Upper;
    point.photo_terms[term].amplitude.assign(layout.photo.size(), 0.0);
    point.photo_terms[term].amplitude[channel] = hard[term];
    point.photo_terms[term].photo_current[1]   = std::move(transition.current);
    point.amplitude[channel] += hard[term];
  }

  const auto emission = upc->EmissionRatios(1, layout.photo[channel].emission, point.transfer[0]);
  const auto first  = upc->TargetCurrentRatios(2, layout.photo[channel].target, *point.photo_terms[0].photo_current[1]);
  const auto second = upc->TargetCurrentRatios(2, layout.photo[channel].target, *point.photo_terms[1].photo_current[1]);
  REQUIRE(emission.size() == upc->SampleCount());
  REQUIRE(first.size() == upc->SampleCount());
  REQUIRE(second.size() == upc->SampleCount());

  std::vector<std::complex<double>> exact(upc->SampleCount(), 0.0);
  std::vector<std::complex<double>> overwritten(upc->SampleCount(), 0.0);
  for (const auto &sample : indices(exact)) {
    exact[sample] =
        std::complex<double>(1.0, 0.0) * emission[sample] * (hard[0] * first[sample] + hard[1] * second[sample]);
    overwritten[sample] = std::complex<double>(1.0, 0.0) * emission[sample] * (hard[0] + hard[1]) * second[sample];
  }
  const auto exact_projection       = gra::nuclear::NuclearGoodWalkerProject(exact, upc->SampleShape());
  const auto overwritten_projection = gra::nuclear::NuclearGoodWalkerProject(overwritten, upc->SampleShape());
  const gra::nuclear::MUPCScreen screen(*upc, layout, point);
  const auto                     result = screen.Result();
  const std::size_t              upper  = layout.photo[channel].Pair() / 2;
  const std::size_t              lower  = layout.photo[channel].Pair() % 2;
  REQUIRE(result.helicity_norm[channel] == Approx(exact_projection[upper][lower]).epsilon(2.0e-11));
  REQUIRE(std::abs(exact_projection[upper][lower] - overwritten_projection[upper][lower]) >
          1.0e-5 * std::max(1.0, exact_projection[upper][lower]));

  gra::nuclear::MUPCScreen accumulator(*upc, layout, point);
  gra::nuclear::LoopNode   node;
  node.weight = 0.2;
  for (const unsigned int failure : {0U, 1U, 2U, 3U, 4U}) {
    auto invalid = point;
    if (failure == 0) { invalid.photo_terms.clear(); }
    if (failure == 1) { invalid.photo_terms[0].amplitude.pop_back(); }
    if (failure == 2) { invalid.photo_terms[0].amplitude[channel] += 0.1; }
    if (failure == 3) { invalid.photo_terms[0].photo_current[1]->sample.pop_back(); }
    if (failure == 4) { invalid.photo_terms[0].photo_current[1]->sample[0] = std::numeric_limits<double>::quiet_NaN(); }
    REQUIRE_THROWS_AS(accumulator.Add(node, invalid), gra::AmplitudeFailure);
    const auto unchanged = accumulator.Result();
    REQUIRE(unchanged.helicity_norm[channel] == Approx(result.helicity_norm[channel]).margin(2.0e-13));
  }
  accumulator.Add(node, point);
  REQUIRE(accumulator.Result().helicity_norm[channel] == Approx(1.44 * result.helicity_norm[channel]).margin(2.0e-12));
}
