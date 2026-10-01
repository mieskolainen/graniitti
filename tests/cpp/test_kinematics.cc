// Unit tests for MKinematics class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <array>
#include <catch.hpp>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <limits>
#include <random>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Nuclear/MSteering.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonOverride.h"

using namespace gra::kinematics;
using namespace gra;

using gra::aux::indices;
using gra::math::pow2;

// Check the reconstructed RAMBO density against exact two-body and massless measures
TEST_CASE("RAMBO inverse density has the analytic phase-space normalization", "[kinematics][Rambo][cascade]") {
  gra::MRandom rng;
  rng.SetSeed(95173);
  const gra::M4Vec mother(0.0, 0.0, 0.0, 5.0);
  for (const auto &masses : std::vector<std::vector<double>>{{0.2, 0.9}, {1.7, 2.3},
                                                            {0.0, 0.0, 0.0, 0.0},
                                                            {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}}) {
    const double exact = masses.size() == 2 ? PS2Massive(25.0, pow2(masses[0]), pow2(masses[1]))
                                             : PSnMassless(25.0, masses.size());
    for (unsigned int event = 0; event < 16; ++event) {
      std::vector<gra::M4Vec> momenta;
      const auto weight = RamboMassive(mother, mother.M(), masses, momenta, rng);
      REQUIRE(weight.Integral() == Approx(exact).epsilon(1e-11));
      REQUIRE(RamboWeight(mother.M(), momenta) == Approx(exact).epsilon(1e-11));
    }
  }
}

struct JsonOverrideRegistryScope {
  // Reset override registry before the test section
  JsonOverrideRegistryScope() { gra::json_override::ClearCardOverrides(); }

  // Clear override registry after the test section
  ~JsonOverrideRegistryScope() { gra::json_override::ClearCardOverrides(); }
};

// Read nuclear mass tables unchanged while JSON model overrides are registered
TEST_CASE("Nuclear steering accepts model overrides with AME text tables", "[input][nuclear][icepack]") {
  JsonOverrideRegistryScope scope;
  const auto card = nlohmann::json::parse(aux::GetInputData(
      aux::ResolveProjectPath("icepack/UPC/EMD/ALICE_2149540/gencard.json")));
  const auto general = nlohmann::json::parse(aux::GetInputData(aux::ResolveProjectPath("modeldata/TUNE0/GENERAL.json")));
  const auto numerics = nlohmann::json::parse(aux::GetInputData(aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json")));
  REQUIRE_NOTHROW(nuclear::ReadUPCSteering(card.at("NUCLEAR"), general.at("PARAM_NUCLEAR"),
                                          numerics.at("NUMERICS_NUCLEAR"), "TUNE0"));
  json_override::RegisterCardOverrides({json_override::ParseSpec("NUMERICS.json:NUMERICS_MC.fast_adaptation=false")});
  (void)aux::GetInputData(aux::ResolveProjectPath("modeldata/TUNE0/NUMERICS.json"));
  json_override::RequireAllCardOverridesApplied();
  REQUIRE_NOTHROW(nuclear::ReadUPCSteering(card.at("NUCLEAR"), general.at("PARAM_NUCLEAR"),
                                          numerics.at("NUMERICS_NUCLEAR"), "TUNE0"));
}

// Sample the pure LS icepack amplitudes with their complete production rows
TEST_CASE("Pure LS icepacks produce finite nonzero event weights", "[input][spin][icepack]") {
  for (const std::string model : {"GP", "XP"}) {
    for (const std::string state : {"scalar", "pseudoscalar", "axialvector", "vector", "tensor", "pseudotensor"}) {
      const auto dataset = nlohmann::json::parse(aux::GetInputData(
          aux::ResolveProjectPath("icepack/SPIN/" + model + "/" + state + "_LS/dataset.json")));
      for (const auto &sample : dataset.at("samples")) {
        CAPTURE(model, state, sample.at("name"));
        JsonOverrideRegistryScope scope;
        auto card = nlohmann::json::parse(aux::GetInputData(aux::ResolveProjectPath(sample.at("gencard").get<std::string>())));
        std::vector<json_override::OverrideSpec> specs;
        for (const auto &[key, value] : sample.at("parameters").items()) {
          specs.push_back(json_override::ParseSpec(key + "=" + value.dump()));
        }
        json_override::ApplyInputOverrides(card, specs);
        json_override::RegisterCardOverrides(specs);
        card["GENERIC"]["CORES"] = 1;
        MGraniitti generator;
        generator.ReadInput(card);
        generator.proc->PrepareRun();
        MRandom random;
        random.SetSeed(48163);
        bool accepted = false;
        for (unsigned int i = 0; i < 10000 && !accepted; ++i) {
          std::vector<double> point(generator.proc->ProcPtr.LIPSDIM);
          for (double &x : point) { x = random.U(0.0, 1.0); }
          MEventWeightState state;
          const double weight = generator.proc->EventWeight(point, state);
          REQUIRE(std::isfinite(weight));
          accepted = state.Valid() && weight > 0.0;
        }
        REQUIRE(accepted);
      }
    }
  }
}

// Sum one particle set into its total four-momentum
gra::M4Vec psum(const std::vector<gra::M4Vec> &p) {
  gra::M4Vec ptot;
  for (const auto &i : indices(p)) { ptot += p[i]; }
  return ptot;
}

// Validate adaptive coordinates against analytic N-body volumes and momentum conservation
TEST_CASE("adaptive central coordinates preserve arbitrary multiplicity phase space", "[PhaseSpace][Adaptive]") {
  MRandom rng;
  rng.SetSeed(73521);
  const double mass = 7.0;
  const M4Vec mother(0.0, 0.0, 0.0, mass);
  for (const std::size_t n : {2, 3, 4, 5}) {
    const std::vector<double> masses(n, 0.0);
    std::vector<double> coordinates(3 * n - 4);
    std::vector<M4Vec> products;
    MCW integral;
    for (std::size_t event = 0; event < 20000; ++event) {
      for (double &unit : coordinates) { unit = rng.U(0.0, 1.0); }
      const auto weight = NBodyPhaseSpace(mother, mass, masses, products, coordinates);
      REQUIRE(weight.GetW() > 0.0);
      REQUIRE(gra::math::CheckEMC(mother - psum(products), 1e-9));
      integral += weight;
    }
    const double exact = PSnMassless(mass * mass, n);
    REQUIRE(std::abs(integral.Integral() - exact) < std::max(6.0 * integral.IntegralError(), 1e-10 * exact));
  }
}


// Integrate the one-particle phase-space measure numerically
double NumericalDPhi1(double max_energy, double mass) {
  const double        pmax = std::sqrt((max_energy - mass) * (max_energy + mass));
  const int           n    = 20000;
  const double        h    = pmax / n;
  std::vector<double> value(n + 1, 0.0);

  for (int i = 0; i <= n; ++i) {
    const double p = i * h;
    const double e = std::sqrt(p * p + mass * mass);
    value[i]       = p * p / e;
  }
  return gra::math::CS13Integral(value, h) / (4.0 * gra::math::pow2(gra::math::PI));
}

struct MCMeanForTest {
  double n    = 0.0;
  double sum  = 0.0;
  double sum2 = 0.0;

  // Add one Monte Carlo weight
  void Add(double weight) {
    ++n;
    sum += weight;
    sum2 += weight * weight;
  }

  // Compute the sample mean
  double Mean() const { return (n > 0.0) ? sum / n : 0.0; }

  // Compute the standard error of the sample mean
  double Error() const {
    if (n <= 1.0) { return 0.0; }
    const double mean = Mean();
    const double var  = std::max(0.0, sum2 / n - mean * mean);
    return std::sqrt(var / n);
  }
};

struct ContinuumPatchCutsForTest {
  double pt_min   = 0.0;
  double pt_max   = 1.5;
  double mass_min = 0.0;
  double mass_max = 8.0;
  double y_min    = -1.6;
  double y_max    = 1.6;
};

struct ContinuumObservablesForTest {
  double mass_volume_fraction   = 0.0;
  double radial_volume_fraction = 0.0;
  double rapidity_fraction      = 0.0;
};

// Expose the production central transverse map to the RAMBO comparison
class CentralMapProbeForTest : public gra::MCentral {
 public:
  // Compute the production Helmert ball radius
  static double Radius2(double mass_max, double recoil_pt2, unsigned int K) {
    return MassConditionedTransverseRadius2(mass_max, recoil_pt2, K);
  }

  // Apply the production logarithmic map and return its exact Jacobian
  static std::pair<std::vector<gra::M4Vec>, double> Map(const std::vector<double> &radial_units,
                                                        const std::vector<double> &angle_units,
                                                        const gra::M4Vec &central_transverse_momentum, double radius2,
                                                        double scale2) {
    double     jacobian = 0.0;
    const auto differences =
        MapTransverseLog(radial_units, angle_units, central_transverse_momentum, radius2, scale2, jacobian);
    return {differences, jacobian};
  }
};

// Compute the central phase space forward and rapidity integration volume
double ContinuumPatchVolumeForTest(unsigned int K, const ContinuumPatchCutsForTest &cuts) {
  return gra::math::pow2(cuts.pt_max - cuts.pt_min) * gra::math::pow2(2.0 * gra::math::PI) *
         std::pow(cuts.y_max - cuts.y_min, K);
}

// Compute the central transverse determinant squared for the central phase space
// map
double CentralTransverseJacobianForTest(unsigned int K) { return 1.0 / gra::math::pow2(static_cast<double>(K)); }

// Map unit coordinates through the production logarithmic Helmert map
std::pair<std::vector<gra::M4Vec>, double> MassConditionedDifferencesForTest(
    unsigned int K, const gra::M4Vec &central_transverse_momentum, double radius2, double scale2, gra::MRandom &rng) {
  const unsigned int  mode_count = K - 1;
  std::vector<double> radial_units(mode_count);
  std::vector<double> angle_units(mode_count);
  for (auto &unit : radial_units) { unit = rng.U(0.0, 1.0); }
  for (auto &unit : angle_units) { unit = rng.U(0.0, 1.0); }
  return CentralMapProbeForTest::Map(radial_units, angle_units, central_transverse_momentum, radius2, scale2);
}

// Build the central transverse momenta from difference vectors
std::vector<gra::M4Vec> BuildCentralTransverseMomentaForTest(const gra::M4Vec &p1f, const gra::M4Vec &p2f,
                                                             const std::vector<gra::M4Vec> &q) {
  std::vector<gra::M4Vec> p;
  if (!gra::kinematics::BuildCentralTransverseMomenta(p1f, p2f, q, p)) { return {}; }
  return p;
}

// Compute one central phase space patch Monte Carlo weight
double CentralPatchWeightForTest(unsigned int K, double sqrt_s, const ContinuumPatchCutsForTest &cuts,
                                 gra::MRandom &rng, ContinuumObservablesForTest &observables) {
  const double pt1  = rng.U(cuts.pt_min, cuts.pt_max);
  const double pt2  = rng.U(cuts.pt_min, cuts.pt_max);
  const double phi1 = rng.U(0.0, 2.0 * gra::math::PI);
  const double phi2 = rng.U(0.0, 2.0 * gra::math::PI);

  gra::M4Vec   p1(pt1 * std::cos(phi1), pt1 * std::sin(phi1), 0.0, 0.0);
  gra::M4Vec   p2(pt2 * std::cos(phi2), pt2 * std::sin(phi2), 0.0, 0.0);
  const double available_energy    = sqrt_s - pt1 - pt2;
  const double mass2_kinematic_max = gra::math::pow2(available_energy) - (p1 + p2).Pt2();
  if (!(mass2_kinematic_max > 0.0)) { return 0.0; }
  const double mass_max = std::min(cuts.mass_max, std::sqrt(mass2_kinematic_max));
  if (!(cuts.mass_min < mass_max)) { return 0.0; }
  const gra::M4Vec central_transverse_momentum = -(p1 + p2);
  const double     recoil_pt2                  = central_transverse_momentum.Pt2();
  const double     radius2                     = CentralMapProbeForTest::Radius2(mass_max, recoil_pt2, K);
  const double     scale                       = 0.01 * mass_max;
  const auto [q, transverse_map]               = MassConditionedDifferencesForTest(
                    K, central_transverse_momentum, radius2, gra::math::pow2(scale) / static_cast<double>(K), rng);

  std::vector<gra::M4Vec> central = BuildCentralTransverseMomentaForTest(p1, p2, q);
  gra::M4Vec              central_sum;
  for (std::size_t i = 0; i < K; ++i) {
    const double y  = rng.U(cuts.y_min, cuts.y_max);
    const double mt = central[i].Pt();
    central[i].SetPzE(mt * std::sinh(y), mt * std::cosh(y));
    central_sum += central[i];
  }
  if (central_sum.M() < cuts.mass_min || central_sum.M() > cuts.mass_max) { return 0.0; }

  const double s   = sqrt_s * sqrt_s;
  const double p1z = gra::kinematics::SolvePz(0.0, 0.0, pt1, pt2, central_sum.Pz(), central_sum.E(), s);
  if (!(p1z >= 0.0)) { return 0.0; }
  const double p2z = -central_sum.Pz() - p1z;
  if (!(p2z <= 0.0)) { return 0.0; }

  const double e1 = std::sqrt(pt1 * pt1 + p1z * p1z);
  const double e2 = std::sqrt(pt2 * pt2 + p2z * p2z);
  if (!(e1 > 0.0) || !(e2 > 0.0)) { return 0.0; }

  const double       jac_long = 1.0 / std::abs(p1z / e1 - p2z / e2);
  const unsigned int N        = K + 2;

  const double density = CentralTransverseJacobianForTest(K) * jac_long * pt1 * pt2 /
                         (std::pow(2.0, static_cast<double>(K + 2)) * e1 * e2 *
                          std::pow(2.0 * gra::math::PI, static_cast<double>(3 * N - 4)));

  double           rho2   = 0.0;
  const gra::M4Vec center = central_sum / static_cast<double>(K);
  for (const auto &momentum : central) { rho2 += (momentum - center).Pt2(); }
  const double radial_fraction       = std::clamp(rho2 / radius2, 0.0, 1.0);
  const double mass_fraction         = (central_sum.M() - cuts.mass_min) / (cuts.mass_max - cuts.mass_min);
  observables.mass_volume_fraction   = std::pow(std::clamp(mass_fraction, 0.0, 1.0), 2.0 * static_cast<double>(K - 1));
  observables.radial_volume_fraction = std::pow(radial_fraction, static_cast<double>(K - 1));
  observables.rapidity_fraction      = (central.front().Rap() - cuts.y_min) / (cuts.y_max - cuts.y_min);

  return density * transverse_map * ContinuumPatchVolumeForTest(K, cuts);
}

// Classify one RAMBO point in the same patch and return its observables
bool RamboPointInContinuumPatchForTest(const std::vector<gra::M4Vec> &p, double sqrt_s,
                                       const ContinuumPatchCutsForTest &cuts,
                                       ContinuumObservablesForTest     &observables) {
  if (p.size() < 4) { return false; }
  if (!(p[0].Pz() > 0.0) || !(p[1].Pz() < 0.0)) { return false; }
  if (p[0].Pt() < cuts.pt_min || p[0].Pt() > cuts.pt_max) { return false; }
  if (p[1].Pt() < cuts.pt_min || p[1].Pt() > cuts.pt_max) { return false; }

  for (std::size_t i = 2; i < p.size(); ++i) {
    const double y = p[i].Rap();
    if (!std::isfinite(y) || y < cuts.y_min || y > cuts.y_max) { return false; }
  }
  gra::M4Vec central_system;
  for (std::size_t i = 2; i < p.size(); ++i) { central_system += p[i]; }
  if (central_system.M() < cuts.mass_min || central_system.M() > cuts.mass_max) { return false; }

  const unsigned int K                   = p.size() - 2;
  const double       available_energy    = sqrt_s - p[0].Pt() - p[1].Pt();
  const double       mass2_kinematic_max = gra::math::pow2(available_energy) - (p[0] + p[1]).Pt2();
  if (!(mass2_kinematic_max > 0.0)) { return false; }
  const double     mass_max = std::min(cuts.mass_max, std::sqrt(mass2_kinematic_max));
  const double     radius2  = CentralMapProbeForTest::Radius2(mass_max, central_system.Pt2(), K);
  const gra::M4Vec center   = central_system / static_cast<double>(K);
  double           rho2     = 0.0;
  for (std::size_t i = 2; i < p.size(); ++i) { rho2 += (p[i] - center).Pt2(); }
  if (rho2 > radius2 * (1.0 + 1.0e-12)) { return false; }

  const double radial_fraction       = std::clamp(rho2 / radius2, 0.0, 1.0);
  const double mass_fraction         = (central_system.M() - cuts.mass_min) / (cuts.mass_max - cuts.mass_min);
  observables.mass_volume_fraction   = std::pow(std::clamp(mass_fraction, 0.0, 1.0), 2.0 * static_cast<double>(K - 1));
  observables.radial_volume_fraction = std::pow(radial_fraction, static_cast<double>(K - 1));
  observables.rapidity_fraction      = (p[2].Rap() - cuts.y_min) / (cuts.y_max - cuts.y_min);

  return true;
}

struct ContinuumEstimateForTest {
  // Initialize equally binned invariant distributions
  explicit ContinuumEstimateForTest(std::size_t bin_count) : mass(bin_count), radial(bin_count), rapidity(bin_count) {}

  MCMeanForTest              total;
  std::vector<MCMeanForTest> mass;
  std::vector<MCMeanForTest> radial;
  std::vector<MCMeanForTest> rapidity;
};

// Add one weighted unit-interval observable to every histogram estimator
void AddContinuumBinForTest(std::vector<MCMeanForTest> &histogram, double value, double weight) {
  std::size_t selected = histogram.size();
  if (weight > 0.0 && value >= 0.0 && value <= 1.0) {
    selected = std::min(static_cast<std::size_t>(value * static_cast<double>(histogram.size())), histogram.size() - 1);
  }
  for (const auto &i : indices(histogram)) { histogram[i].Add(i == selected ? weight : 0.0); }
}

// Add one weighted continuum point to the integral and differential estimates
void AddContinuumPointForTest(ContinuumEstimateForTest &estimate, const ContinuumObservablesForTest &observables,
                              double weight) {
  estimate.total.Add(weight);
  AddContinuumBinForTest(estimate.mass, observables.mass_volume_fraction, weight);
  AddContinuumBinForTest(estimate.radial, observables.radial_volume_fraction, weight);
  AddContinuumBinForTest(estimate.rapidity, observables.rapidity_fraction, weight);
}

// Require two independent weighted Monte Carlo estimates to be compatible
void RequireContinuumMeanForTest(const MCMeanForTest &direct, const MCMeanForTest &rambo) {
  const double combined_error = std::hypot(direct.Error(), rambo.Error());
  const double roundoff =
      100.0 * std::numeric_limits<double>::epsilon() *
      std::max({std::numeric_limits<double>::min(), std::abs(direct.Mean()), std::abs(rambo.Mean())});
  const double pull = combined_error > 0.0 ? std::abs(direct.Mean() - rambo.Mean()) / combined_error : 0.0;
  CAPTURE(direct.Mean(), direct.Error(), rambo.Mean(), rambo.Error(), combined_error, roundoff, pull);
  REQUIRE(std::abs(direct.Mean() - rambo.Mean()) <= 6.0 * combined_error + roundoff);
}

// Compare all bins of one direct and RAMBO differential estimate
void RequireContinuumHistogramForTest(const std::vector<MCMeanForTest> &direct, const std::vector<MCMeanForTest> &rambo,
                                      const char *observable) {
  REQUIRE(direct.size() == rambo.size());
  for (const auto &i : indices(direct)) {
    CAPTURE(observable, i);
    RequireContinuumMeanForTest(direct[i], rambo[i]);
  }
}

// Require the differential bins to contain the complete accepted integral
void RequireContinuumCoverageForTest(const MCMeanForTest &total, const std::vector<MCMeanForTest> &histogram,
                                     const char *observable) {
  double sum = 0.0;
  for (const auto &bin : histogram) { sum += bin.Mean(); }
  CAPTURE(observable, sum, total.Mean());
  REQUIRE(total.Mean() > 0.0);
  REQUIRE(sum == Approx(total.Mean()).epsilon(1.0e-12));
}

// Compare the central phase space patch against RAMBO for one central
// multiplicity
void ContinuumPatchRamboTest(unsigned int K, gra::MRandom &rng) {
  CAPTURE(K);
  const double                    sqrt_s    = 12.0;
  constexpr int                   samples   = 240000;
  constexpr std::size_t           bin_count = 6;
  const ContinuumPatchCutsForTest cuts;
  const gra::M4Vec                mother(0.0, 0.0, 0.0, sqrt_s);
  const unsigned int              N = K + 2;

  ContinuumEstimateForTest direct(bin_count);
  ContinuumEstimateForTest rambo(bin_count);
  std::vector<gra::M4Vec>  p(N);
  int                      invalid_rambo = 0;

  for (int i = 0; i < samples; ++i) {
    ContinuumObservablesForTest direct_observables;
    const double                direct_weight = CentralPatchWeightForTest(K, sqrt_s, cuts, rng, direct_observables);
    AddContinuumPointForTest(direct, direct_observables, direct_weight);

    const gra::kinematics::MCW rw = gra::kinematics::RamboMassless(mother, sqrt_s, p, rng, false);
    if (rw.GetW() > 0.0) {
      ContinuumObservablesForTest rambo_observables;
      const double                rambo_weight =
          RamboPointInContinuumPatchForTest(p, sqrt_s, cuts, rambo_observables) ? rw.GetW() : 0.0;
      AddContinuumPointForTest(rambo, rambo_observables, rambo_weight);
    } else {
      ++invalid_rambo;
    }
  }
  REQUIRE(invalid_rambo == 0);

  RequireContinuumCoverageForTest(direct.total, direct.mass, "direct mass");
  RequireContinuumCoverageForTest(direct.total, direct.radial, "direct radial");
  RequireContinuumCoverageForTest(direct.total, direct.rapidity, "direct rapidity");
  RequireContinuumCoverageForTest(rambo.total, rambo.mass, "RAMBO mass");
  RequireContinuumCoverageForTest(rambo.total, rambo.radial, "RAMBO radial");
  RequireContinuumCoverageForTest(rambo.total, rambo.rapidity, "RAMBO rapidity");
  RequireContinuumMeanForTest(direct.total, rambo.total);
  RequireContinuumHistogramForTest(direct.mass, rambo.mass, "central mass volume");
  RequireContinuumHistogramForTest(direct.radial, rambo.radial, "Helmert radial volume");
  RequireContinuumHistogramForTest(direct.rapidity, rambo.rapidity, "first central rapidity");
}

// MKinematics

// Test generic JSON override path parsing and input card mutation
TEST_CASE("JSON card overrides parse and apply nested input paths", "[json-override]") {
  using gra::json_override::ApplyInputOverrides;
  using gra::json_override::ParseSpecs;

  nlohmann::json card = {{"GENERIC", {{"NEVENTS", 100}, {"MODELPARAM", "TUNE0"}}},
                         {"FIDCUTS", {{"CENTRAL", {{"*", {{"Eta", {-2.5, 2.5}}}}}}}},
                         {"MATRIX", {{1, 2}, {3, 4}}},
                         {"weird.key", {{"a:b", 1}}}};

  const auto specs = ParseSpecs({"GENERIC.NEVENTS=10", "FIDCUTS.CENTRAL[\"*\"].Eta[0]=-3.0",
                                 "FIDCUTS.CENTRAL[\"*\"].Et=[0.0,100000.0]", "MATRIX[1,0]=5",
                                 "GENERIC.MODELPARAM=\"TUNE1\"", "[\"weird.key\"][\"a:b\"]=2"});

  ApplyInputOverrides(card, specs);

  REQUIRE(card["GENERIC"]["NEVENTS"].get<int>() == 10);
  REQUIRE(card["FIDCUTS"]["CENTRAL"]["*"]["Eta"][0].get<double>() == Approx(-3.0));
  REQUIRE(card["FIDCUTS"]["CENTRAL"]["*"]["Et"][1].get<double>() == Approx(100000.0));
  REQUIRE(card["MATRIX"][1][0].get<int>() == 5);
  REQUIRE(card["GENERIC"]["MODELPARAM"].get<std::string>() == "TUNE1");
  REQUIRE(card["weird.key"]["a:b"].get<int>() == 2);
}

// Test JSON object key renaming through the command line override syntax
TEST_CASE("JSON card overrides rename object keys", "[json-override]") {
  using gra::json_override::ApplyInputOverrides;
  using gra::json_override::ParseSpecs;

  nlohmann::json card = {
      {"MODELS", {{"MP", {{"[993,993]", {{"basis", "auto_min_L"}, {"polarization", {{"mode", "a_Jz"}}}}}}}}}};
  const auto specs =
      ParseSpecs({"MODELS.MP[\"[993,993]\"].@key=\"[991,991]\"", "MODELS.MP[\"[991,991]\"].polarization.mode=\"rho\""});

  ApplyInputOverrides(card, specs);

  REQUIRE_FALSE(card["MODELS"]["MP"].contains("[993,993]"));
  REQUIRE(card["MODELS"]["MP"]["[991,991]"]["polarization"]["mode"].get<std::string>() == "rho");

  const auto collision              = ParseSpecs({"MODELS.MP[\"[991,991]\"].@key=\"[995,995]\""});
  card["MODELS"]["MP"]["[995,995]"] = nlohmann::json::object();
  REQUIRE_THROWS_AS(ApplyInputOverrides(card, collision), std::invalid_argument);
}

// Test registered JSON override application through card file reads
TEST_CASE("JSON card overrides apply to selected model card reads", "[json-override]") {
  JsonOverrideRegistryScope registry_scope;

  using gra::json_override::ParseSpecs;
  using gra::json_override::RegisterCardOverrides;
  using gra::json_override::RequireAllCardOverridesApplied;
  using gra::json_override::ResolveCardOverrideTargets;
  using gra::json_override::ValidateCardOverrideSelectors;

  const auto specs = ParseSpecs({"GENERAL.json:PARAM_SOFT.MODEL.single.EXCHANGE.P.g[0,0]=8.4",
                                 "NUMERICS.json:NUMERICS_DURHAM.HELICITY_PROJECTOR=\"transverse\""});

  ValidateCardOverrideSelectors(specs, std::filesystem::path("modeldata/TUNE0"));
  const auto targets = ResolveCardOverrideTargets(specs, std::filesystem::path("modeldata/TUNE0"));
  REQUIRE(targets.size() == 2);
  RegisterCardOverrides(specs);

  for (const auto &target : targets) { (void)gra::aux::GetInputData(target); }
  REQUIRE_NOTHROW(RequireAllCardOverridesApplied());

  auto                       &registry = gra::json_override::OverrideRegistry();
  std::lock_guard<std::mutex> lock(registry.mutex);
  REQUIRE(registry.overrides.front().applied == 1);
  REQUIRE(registry.overrides.back().applied == 1);
  REQUIRE(registry.resolved_cards.size() == 2);
}

// Test invalid JSON override syntax and unresolved targets
TEST_CASE("JSON card overrides reject malformed or missing targets", "[json-override]") {
  JsonOverrideRegistryScope registry_scope;

  using gra::json_override::ApplyInputOverrides;
  using gra::json_override::ParseSpec;
  using gra::json_override::ParseSpecs;
  using gra::json_override::ValidateCardOverrideSelectors;

  REQUIRE_THROWS_AS(ParseSpec("GENERIC.NEVENTS"), std::invalid_argument);
  REQUIRE_THROWS_AS(ParseSpec("GENERIC.NEVENTS=abc"), std::invalid_argument);
  REQUIRE_THROWS_AS(ParseSpec("A.B[0,]=1"), std::invalid_argument);
  REQUIRE_THROWS_AS(ParseSpec(".A.B=1"), std::invalid_argument);
  REQUIRE_THROWS_AS(ParseSpec("A..B=1"), std::invalid_argument);

  nlohmann::json card              = {{"A", {{"B", {1}}}}};
  const auto     missing_key_specs = ParseSpecs({"A.C.D=1"});
  REQUIRE_THROWS_AS(ApplyInputOverrides(card, missing_key_specs), std::invalid_argument);

  const auto bad_selector_specs = ParseSpecs({"MISSING.json:A=1"});
  REQUIRE_THROWS_AS(ValidateCardOverrideSelectors(bad_selector_specs, std::filesystem::path("modeldata/TUNE0")),
                    std::invalid_argument);
}

TEST_CASE("Inline command parser accepts colon and equals singlet values", "[aux]") {
  const auto colon = gra::aux::SplitCommands("@MMAX:0 @SPINGEN:false");
  REQUIRE(colon.size() == 2);
  REQUIRE(colon[0].id == "MMAX");
  REQUIRE(colon[0].arg.at("_SINGLET_") == "0");
  REQUIRE(colon[1].id == "SPINGEN");
  REQUIRE(colon[1].arg.at("_SINGLET_") == "false");

  const auto equals = gra::aux::SplitCommands("@MMAX=0");
  REQUIRE(equals.size() == 1);
  REQUIRE(equals[0].id == "MMAX");
  REQUIRE(equals[0].arg.at("_SINGLET_") == "0");
}

TEST_CASE("Auxiliary numerical validators reject invalid domains", "[aux]") {
  CHECK(gra::aux::AssertRatio(2.0, 2.0, 0.0));
  CHECK_FALSE(gra::aux::AssertRatio(1.0, 0.0, 0.1));
  CHECK_FALSE(gra::aux::AssertRatio(1.0, 1.0, -0.1));
  CHECK_FALSE(gra::aux::AssertRange(1.0, std::vector<double>{2.0, 0.0}));
  CHECK_THROWS_AS(gra::aux::AssertCutRange(std::vector<double>{0.0, 1.0}, std::vector<double>{0.0}),
                  std::invalid_argument);
}

TEST_CASE("Forward leg state rejects invalid beam labels", "[Kinematics]") {
  gra::ForwardLegState state;
  state.leg = static_cast<gra::ForwardBeamLeg>(17);
  CHECK_THROWS_AS(state.Index(), std::invalid_argument);
  CHECK_THROWS_AS(state.FinalIndex(), std::invalid_argument);
}

TEST_CASE("Relativistic BW mass-squared proposal normalization matches sampler", "[MRandom][mass-proposal]") {
  const double m0    = 0.77526;
  const double width = 0.1491;
  const double limit = 5.0;
  const double mmin  = 2.0 * 0.13957;

  const double lower_mass = std::max(mmin, std::max(0.0, m0 - limit * width));
  const double lower2     = pow2(lower_mass);
  const double upper2     = pow2(m0 + limit * width);
  const double norm       = MRandom::RelativisticBWMass2Integral(m0, width, lower2, upper2);
  REQUIRE_THROWS_AS(MRandom::RelativisticBWMass2Integral(m0, 0.0, lower2, upper2), std::invalid_argument);
  REQUIRE_THROWS_AS(MRandom::RelativisticBWMass2Integral(m0, width, lower2, -1.0), std::invalid_argument);
  REQUIRE(MRandom::RelativisticBWMass2Integral(m0, width, upper2, lower2) == Approx(0.0));

  const double pole2    = pow2(m0);
  const double scale    = m0 * width;
  const double expected = (std::atan2(upper2 - pole2, scale) - std::atan2(lower2 - pole2, scale)) / scale;

  REQUIRE(norm == Approx(expected).epsilon(1e-14));
  REQUIRE(norm > 0.0);
  REQUIRE((upper2 - lower2) > 0.0);

  const double first_moment =
      pole2 * norm + 0.5 * std::log((pow2(upper2 - pole2) + pow2(scale)) / (pow2(lower2 - pole2) + pow2(scale)));
  const double expected_mean_s = first_moment / norm;

  MRandom rng;
  rng.SetSeed(12345);

  const unsigned int samples = 30000;
  double             mean_s  = 0.0;
  for (unsigned int i = 0; i < samples; ++i) {
    const double mass = rng.RelativisticBWRandom(m0, width, limit, mmin);
    const double s    = pow2(mass);
    REQUIRE(s >= lower2);
    REQUIRE(s <= upper2);
    mean_s += s;
  }
  mean_s /= samples;

  REQUIRE(mean_s == Approx(expected_mean_s).epsilon(0.02));

  // A broad state whose nominal lower window crosses zero must start at
  // threshold
  const double broad_m0    = 0.475;
  const double broad_width = 0.55;
  double       sampled_min = std::numeric_limits<double>::infinity();
  for (unsigned int i = 0; i < 1000; ++i) {
    const double mass = rng.RelativisticBWRandom(broad_m0, broad_width, limit, mmin);
    REQUIRE(mass >= mmin);
    sampled_min = std::min(sampled_min, mass);
  }
  REQUIRE(sampled_min < broad_m0);
}

// Resolve positive tail measures even when both ordinary arctangents round to pi/2
TEST_CASE("Breit-Wigner tail normalization and quantiles retain their mass support",
          "[MRandom][mass-proposal][tails]") {
  constexpr double width = 1e-20;
  for (const auto &[pole, lower, upper] : std::vector<std::array<double, 3>>{
           {1.0, 4.0, 16.0}, {4.0, 1.0, 4.0},
           {1.0, 4.0, std::nextafter(4.0, 5.0)}}) {
    // The narrow tail converges to the integral of 1/(s-M^2)^2
    const double expected = (upper - lower) / ((lower - pole * pole) * (upper - pole * pole));
    const double integral = MRandom::RelativisticBWMass2Integral(pole, width, lower, upper);
    REQUIRE(integral > 0.0);
    CHECK(integral == Approx(expected).epsilon(1e-13));
  }

  for (const bool relativistic : {false, true}) {
    MRandom random;
    random.SetSeed(91453);
    MRandom reference = random;
    for (unsigned int i = 0; i < 1000; ++i) {
      const double unit = reference.U(0.0, 1.0);
      const double lower = relativistic ? 4.0 : 2.0;
      const double upper = relativistic ? 16.0 : 4.0;
      const double inverse = (1.0 - unit) / (lower - 1.0) + unit / (upper - 1.0);
      const double quantile = 1.0 + 1.0 / inverse;
      const double mass = relativistic ? random.RelativisticBWRandom(1.0, width, 3e20, 2.0)
                                       : random.CauchyRandom(1.0, width, 3e20, 2.0);
      REQUIRE(mass >= 2.0);
      REQUIRE(mass <= 4.0);
      CHECK((relativistic ? mass * mass : mass) == Approx(quantile).epsilon(1e-12));
    }
  }
}

// Check that reseeding clears cached Gaussian draws and zero Poisson consumes no draws
TEST_CASE("Random reseeding and zero Poisson preserve the requested stream", "[MRandom]") {
  MRandom random;
  random.SetSeed(17);
  const double first = random.G(0.0, 1.0);
  random.SetSeed(17);
  CHECK(random.G(0.0, 1.0) == Approx(first).margin(0.0));

  MRandom copied = random;
  MRandom fresh;
  copied.SetSeed(29);
  fresh.SetSeed(29);
  for (int draw = 0; draw < 8; ++draw) {
    CHECK(copied.G(0.0, 1.0) == Approx(fresh.G(0.0, 1.0)).margin(0.0));
  }
  CHECK(random.PoissonRandom(0.0) == 0);
  random.SetSeed(53);
  fresh.SetSeed(53);
  CHECK(random.PoissonRandom(0.0) == 0);
  CHECK(random.U(0.0, 1.0) == Approx(fresh.U(0.0, 1.0)).margin(0.0));
}

// Check uniform support without subtracting bounds of opposite large magnitude
TEST_CASE("Uniform draws remain finite and exclude the upper bound", "[MRandom]") {
  MRandom random;
  const double lower = 1.0;
  const double upper = std::nextafter(lower, 2.0);
  for (int draw = 0; draw < 64; ++draw) {
    CHECK(random.U(lower, upper) < upper);
    const double value = random.U(-1e308, 1e308);
    CHECK(std::isfinite(value));
    CHECK(value >= -1e308);
    CHECK(value < 1e308);
  }
}

// Check the uniform limit of a bounded exponential at a tiny rate times range
TEST_CASE("Bounded exponential sampling resolves tiny dimensionless rates", "[MRandom]") {
  MRandom random;
  random.SetSeed(61);
  double mean = 0.0;
  for (int draw = 0; draw < 4096; ++draw) {
    const double x = random.ExpBoundedRandom(0.0, 1e-100, 1e-300);
    REQUIRE(std::isfinite(x));
    REQUIRE(x >= 0.0);
    REQUIRE(x <= 1e-100);
    CHECK(random.ExpBoundedPdf(x, 0.0, 1e-100, 1e-300) * 1e-100 == Approx(1.0));
    mean += x / 1e-100;
  }
  CHECK(mean / 4096 == Approx(0.5).margin(0.03));
  CHECK_THROWS_AS(random.ExpBoundedPdf(std::numeric_limits<double>::quiet_NaN(), 0.0, 1.0, 1.0),
                  std::invalid_argument);
}

// Check the logarithmic law in its small probability limit
TEST_CASE("Logarithmic multiplicities remain normalized for small probabilities", "[MRandom]") {
  MRandom random;
  for (const double p : {1e-20, 1e-200, std::numeric_limits<double>::denorm_min()}) {
    REQUIRE(random.LogPdf(1, p) == Approx(1.0));
    CHECK(random.LogRandom(p, 1) == 1);
  }
  double sum = 0.0;
  for (int n = 1; n <= 200; ++n) { sum += random.LogPdf(n, 0.5); }
  CHECK(sum == Approx(1.0).epsilon(1e-12));
}

// Check sparse Dirichlet draws against their analytic first and second moments
TEST_CASE("Sparse Dirichlet draws stay on the probability simplex", "[MRandom]") {
  MRandom random;
  random.SetSeed(0);
  constexpr int count = 8192;
  for (const double alpha : {0.001, 0.0001}) {
    double mean = 0.0;
    double second = 0.0;
    for (int draw = 0; draw < count; ++draw) {
      std::vector<double> values;
      random.DirRandom({alpha, alpha}, values);
      REQUIRE(values.size() == 2);
      REQUIRE(std::isfinite(values[0]));
      REQUIRE(std::isfinite(values[1]));
      REQUIRE(values[0] >= 0.0);
      REQUIRE(values[1] >= 0.0);
      CHECK(values[0] + values[1] == Approx(1.0).margin(1e-14));
      mean += values[0];
      second += values[0] * values[0];
    }
    CHECK(mean / count == Approx(0.5).margin(0.035));
    CHECK(second / count == Approx((alpha + 1.0) / (2.0 * (2.0 * alpha + 1.0))).margin(0.035));
  }
}

// Check truncated NBD sampling even when its unconditioned probability is tiny
TEST_CASE("Truncated NBD multiplicities follow their conditional law", "[MRandom]") {
  MRandom random;
  random.SetSeed(31);
  constexpr int count = 8192;
  std::vector<int> frequencies(11, 0);
  for (int draw = 0; draw < count; ++draw) {
    const int n = random.NBDRandom(1e9, 1.0, 10);
    REQUIRE(n >= 0);
    REQUIRE(n <= 10);
    ++frequencies[n];
  }
  // The truncated geometric law tends to a uniform distribution at large mean
  for (const int frequency : frequencies) {
    CHECK(frequency == Approx(count / 11.0).margin(140.0));
  }
  double mean = 0.0;
  for (int draw = 0; draw < count; ++draw) { mean += random.NBDRandom(4.0, 2.0, 8); }
  double norm = 0.0;
  double expected = 0.0;
  for (int n = 0; n <= 8; ++n) {
    const double probability = (n + 1.0) * std::pow(2.0 / 3.0, n) / 9.0;
    norm += probability;
    expected += n * probability;
  }
  CHECK(mean / count == Approx(expected / norm).margin(0.16));
  CHECK(random.NBDRandom(0.0, 2.0, 8) == 0);
  CHECK(random.NBDRandom(4.0, 2.0, 0) == 0);

  // At k=2 and very large mean the truncated weights are proportional to n+1
  constexpr int maximum = 1000000000;
  mean = 0.0;
  for (int draw = 0; draw < 2048; ++draw) {
    const int n = random.NBDRandom(1e300, 2.0, maximum);
    REQUIRE(n >= 0);
    REQUIRE(n <= maximum);
    mean += static_cast<double>(n) / maximum;
  }
  CHECK(mean / 2048 == Approx(2.0 / 3.0).margin(0.025));
}

// Check the NBD Poisson limit without subtracting nearly equal gamma functions
TEST_CASE("NBD probabilities and samples resolve the large shape limit", "[MRandom]") {
  MRandom random;
  random.SetSeed(71);
  for (const double shape : {1e20, 1e100, 1e300}) {
    double sum = 0.0;
    for (int n = 0; n <= 64; ++n) {
      const double poisson = std::exp(-10.0 + n * std::log(10.0) - gra::math::LogGamma(n + 1.0));
      CHECK(random.NBDPdf(n, 10.0, shape) == Approx(poisson).epsilon(1e-10));
      sum += random.NBDPdf(n, 10.0, shape);
    }
    CHECK(sum == Approx(1.0).epsilon(1e-10));
    double mean = 0.0;
    for (int draw = 0; draw < 4096; ++draw) { mean += random.NBDRandom(10.0, shape, 64); }
    CHECK(mean / 4096 == Approx(10.0).margin(0.3));
  }
  // Cover the finite shape transition to the asymptotic gamma ratio
  for (const double shape : {15.0, 16.0, 32.0}) {
    double sum = 0.0;
    for (int n = 0; n <= 300; ++n) { sum += random.NBDPdf(n, 30.0, shape); }
    CHECK(sum == Approx(1.0).epsilon(1e-10));
  }
}

// Check bounded power laws through their probability integral transform
TEST_CASE("Bounded power laws retain support and normalization at extreme slopes", "[MRandom]") {
  MRandom random;
  random.SetSeed(47);
  CHECK(random.PowerBoundedPdf(10.0, 10.0, 100.0, 400.0) == Approx(39.9));
  CHECK(random.PowerBoundedPdf(100.0, 10.0, 100.0, -400.0) == Approx(4.01));
  const double narrow_a = 1e300;
  const double narrow_b = std::nextafter(narrow_a, std::numeric_limits<double>::infinity());
  for (const double alpha : {-400.0, 1.0, 400.0}) {
    CHECK(random.PowerBoundedPdf(narrow_a, narrow_a, narrow_b, alpha) * (narrow_b - narrow_a) ==
          Approx(1.0).epsilon(1e-10));
  }
  const double log_span = 400.0 * std::log(10.0);
  CHECK(random.PowerBoundedPdf(1e-100, 1e-200, 1e200, 1.0 - 1e-12) * 1e-100 * log_span ==
        Approx(1.0).epsilon(1e-8));
  constexpr int count = 4096;
  for (const double alpha : {400.0, -400.0, 1.0, 1.0 - 1e-10, 1.0 + 1e-10}) {
    double mean_cdf = 0.0;
    for (int draw = 0; draw < count; ++draw) {
      const double x = random.PowerBoundedRandom(10.0, 100.0, alpha);
      REQUIRE(std::isfinite(x));
      REQUIRE(x >= 10.0);
      REQUIRE(x <= 100.0);
      CHECK(std::isfinite(random.PowerBoundedPdf(x, 10.0, 100.0, alpha)));
      double cdf = std::log(x / 10.0) / std::log(10.0);
      if (alpha > 2.0) { cdf = 1.0 - std::pow(10.0 / x, 399.0); }
      if (alpha < 0.0) { cdf = std::pow(x / 100.0, 401.0); }
      mean_cdf += cdf;
    }
    CHECK(mean_cdf / count == Approx(0.5).margin(0.03));
  }
}

TEST_CASE("Numerical utility preconditions reject invalid domains", "[MMath][MMatrix][MRandom]") {
  SECTION("empty combinatorics and zero-rank tensors") {
    std::vector<int> empty;
    const auto       permutations = gra::math::Permutations(empty);
    REQUIRE(permutations.size() == 1);
    REQUIRE(permutations.front().empty());
    REQUIRE_THROWS_AS(gra::math::EpsTensor(0), std::invalid_argument);

    const std::vector<int> repeated            = {2, 1, 1};
    const auto             unique_permutations = gra::math::Permutations(repeated);
    REQUIRE(unique_permutations.size() == 3);
    REQUIRE((unique_permutations.front() == std::vector<int>{1, 1, 2}));
    REQUIRE_THROWS_AS(gra::math::GetAmpPerm(3, 0), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::Ind2Vec(4, 2), std::overflow_error);
    REQUIRE(gra::math::Cbinom(2, 3) == 0);

    const auto increasing = gra::math::arange<std::vector>(0.0, 0.3, 1.0);
    REQUIRE(increasing.size() == 4);
    const auto decreasing = gra::math::arange<std::vector>(1.0, -0.4, -0.1);
    REQUIRE(decreasing.size() == 3);
    const auto decreasing_integer = gra::math::arange<std::vector>(3, -2, -4);
    REQUIRE(decreasing_integer == std::vector<int>{3, 1, -1, -3});
    const auto rounded = gra::math::arange<std::vector>(0.3, 0.1, 0.4);
    REQUIRE(rounded.size() == 1);
    REQUIRE(rounded.back() < 0.4);
    const auto reversed = gra::math::arange<std::vector>(-0.3, -0.1, -0.4);
    REQUIRE(reversed.size() == 1);
    REQUIRE(reversed.back() > -0.4);
    const auto adjacent = gra::math::arange<std::vector>(0.3, 0.1, std::nextafter(0.4, 1.0));
    REQUIRE(adjacent.size() == 2);
    REQUIRE(adjacent.back() < std::nextafter(0.4, 1.0));
    REQUIRE(gra::math::arange<std::vector>(0.0, 1.0, 1e-20).size() == 1);
    REQUIRE_THROWS_AS(gra::math::arange<std::vector>(0.0, 0.0, 1.0), std::invalid_argument);
    REQUIRE_FALSE(gra::math::CheckEMC(M4Vec(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN())));
  }

  SECTION("arithmetic progressions at integer and floating limits") {
    const auto maximum = std::numeric_limits<long long>::max();
    const auto minimum = std::numeric_limits<long long>::min();
    const auto increasing = gra::math::arange<std::vector>(0LL, maximum / 2, maximum);
    REQUIRE(increasing == std::vector<long long>{0LL, maximum / 2, maximum - 1});
    const auto decreasing = gra::math::arange<std::vector>(0LL, -(maximum / 2), -maximum);
    REQUIRE(decreasing == std::vector<long long>{0LL, -(maximum / 2), -(maximum - 1)});
    const auto crossing = gra::math::arange<std::vector>(maximum, minimum, minimum);
    REQUIRE(crossing == std::vector<long long>{maximum, -1LL});
    const auto unsigned_maximum = std::numeric_limits<unsigned long long>::max();
    const auto unsigned_range = gra::math::arange<std::vector>(0ULL, unsigned_maximum / 2, unsigned_maximum);
    REQUIRE(unsigned_range == std::vector<unsigned long long>{0ULL, unsigned_maximum / 2, unsigned_maximum - 1});
    const auto wide = gra::math::arange<std::vector>(-1e308, 1e308, 1.7e308);
    REQUIRE(wide.size() == 3);
    REQUIRE(wide.back() == Approx(1e308));
  }

  SECTION("unordered cut bounds") {
    const double nan = std::numeric_limits<double>::quiet_NaN();
    for (const auto &cut : {std::vector<double>{nan, 1.0}, std::vector<double>{0.0, nan}}) {
      REQUIRE_FALSE(gra::aux::AssertCut(cut));
      REQUIRE_THROWS_AS(gra::aux::AssertCut(cut, "test", true), std::invalid_argument);
      REQUIRE_FALSE(gra::aux::AssertCutRange(cut, std::vector<double>{0.0, 2.0}));
      REQUIRE_FALSE(gra::aux::AssertCutRange(std::vector<double>{0.2, 0.8}, cut));
      REQUIRE_THROWS_AS(gra::aux::AssertCutRange(std::vector<double>{0.2, 0.8}, cut, "test", true),
                        std::invalid_argument);
    }
    REQUIRE_FALSE(gra::aux::AssertCutRange(std::vector<double>{0.2, 0.8}, std::vector<double>{1.0, 0.0}));
  }

  SECTION("matrix row indexing") {
    gra::MMatrix<double> matrix(2, 2, 0.0);
    matrix[1][1] = 3.0;
    REQUIRE(matrix[1][1] == Approx(3.0));
    REQUIRE_THROWS_AS(matrix[2][0], std::out_of_range);
  }

  SECTION("random distribution parameters") {
    gra::MRandom random;
    REQUIRE(random.PoissonRandom(0.0) == 0);
    REQUIRE_THROWS_AS(random.PoissonRandom(-1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.PoissonRandom(std::numeric_limits<double>::infinity()), std::invalid_argument);
    REQUIRE_THROWS_AS(random.ExpRandom(0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.ExpPdf(1.0, -1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.ExpBoundedRandom(2.0, 1.0, 1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.ExpBoundedPdf(1.5, 1.0, 2.0, 0.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.U(2.0, 1.0), std::invalid_argument);
    REQUIRE_THROWS_AS(random.G(0.0, -1.0), std::invalid_argument);
    REQUIRE(random.NBDPdf(0, 2.0, 1.5) == Approx(std::pow(1.5 / 3.5, 1.5)));
    REQUIRE(random.NBDPdf(-1, 2.0, 1.5) == Approx(0.0));
    const double sample = random.ExpBoundedRandom(1.0, 2.0, 3.0);
    REQUIRE(sample >= 1.0);
    REQUIRE(sample <= 2.0);
  }
}

TEST_CASE("Discrete math algorithms preserve their defining maps", "[MMath][Combinatorics]") {
  for (unsigned int value = 0; value < 64; ++value) {
    REQUIRE(gra::math::Gray2Binary(gra::math::Binary2Gray(value)) == value);
    const std::vector<bool> bits = gra::math::Ind2Vec(value, 6);
    REQUIRE(gra::math::Vec2Ind(bits) == static_cast<int>(value));
  }

  REQUIRE(gra::math::LRsequence(3) == std::vector<unsigned int>{0, 4, 2, 6, 1, 5, 3, 7});
  REQUIRE(gra::math::BinaryMatrix(2) == std::vector<std::vector<int>>{{0, 0}, {0, 1}, {1, 0}, {1, 1}});

  const std::vector<std::vector<int>> axis = {{1, 2}, {3}, {4, 5}};
  std::vector<std::vector<int>>       product;
  gra::math::IndexComb(axis, 0, {}, product);
  REQUIRE(product == std::vector<std::vector<int>>{{1, 3, 4}, {1, 3, 5}, {2, 3, 4}, {2, 3, 5}});

  const auto epsilon = gra::math::EpsTensor(3);
  REQUIRE(epsilon(std::vector<std::size_t>{0, 1, 2}) == 1);
  REQUIRE(epsilon(std::vector<std::size_t>{1, 0, 2}) == -1);
  REQUIRE(epsilon(std::vector<std::size_t>{0, 0, 2}) == 0);

  REQUIRE(gra::math::Cbinom(8, 3) == 56);
  REQUIRE(gra::math::Cbinom(8, 5) == 56);
  REQUIRE_THROWS_AS(gra::math::Cbinom(34, 17), std::overflow_error);
  const std::vector<int> reference = {2, 4, 7, 9};
  REQUIRE(gra::math::PermutationSign(reference, reference) == 1);
  REQUIRE(gra::math::PermutationSign(reference, std::vector<int>{4, 2, 7, 9}) == -1);
  REQUIRE(gra::math::PermutationSign(reference, std::vector<int>{9, 7, 4, 2}) == 1);
  REQUIRE_THROWS_AS(gra::math::PermutationSign(reference, std::vector<int>{2, 4, 7, 7}), std::invalid_argument);
  REQUIRE(gra::math::GetAmpPerm(4, 0).size() == 16);
  REQUIRE(gra::math::GetAmpPerm(4, 1).size() == 24);
}

TEST_CASE("Generic construction algorithms cover serial and parallel limits", "[MMath][Algorithms]") {
  REQUIRE(gra::math::IntegerPower(2.0, 0) == Approx(1.0));
  REQUIRE(gra::math::IntegerPower(-2.0, 5) == Approx(-32.0));
  const auto                powers          = gra::math::IntegerPowers(-2.0, 4);
  const std::vector<double> expected_powers = {1.0, -2.0, 4.0, -8.0, 16.0};
  REQUIRE(powers.size() == expected_powers.size());
  for (const auto &i : indices(powers)) { REQUIRE(powers[i] == Approx(expected_powers[i])); }

  const auto singleton = gra::math::linspace(2.0, 5.0, 1);
  REQUIRE(singleton.size() == 1);
  REQUIRE(singleton.front() == Approx(2.0));
  const auto                grid          = gra::math::linspace(0.0, 1.0, 5);
  const std::vector<double> expected_grid = {0.0, 0.25, 0.5, 0.75, 1.0};
  REQUIRE(grid.size() == expected_grid.size());
  for (const auto &i : indices(grid)) { REQUIRE(grid[i] == Approx(expected_grid[i])); }
  REQUIRE_THROWS_AS(gra::math::linspace(0.0, 1.0, 0), std::invalid_argument);

  std::vector<unsigned int> visit(32, 0);
  gra::math::ParallelFor(visit.size(), [&](const std::size_t i) { ++visit[i]; });
  REQUIRE(std::all_of(visit.cbegin(), visit.cend(), [](const unsigned int count) { return count == 1; }));
  REQUIRE_THROWS_AS(gra::math::ParallelFor(8,
                                           [](const std::size_t i) {
                                             if (i == 3) { throw std::runtime_error("parallel failure"); }
                                           }),
                    std::runtime_error);
}

TEST_CASE("pcm2 obeys two-body center-of-mass kinematics", "[pcm2][Kinematics]") {
  const double sqrt_s = 10.0;
  const double s      = sqrt_s * sqrt_s;
  const double m1     = 1.0;
  const double m2     = 2.0;
  const double e1     = (s + m1 * m1 - m2 * m2) / (2.0 * sqrt_s);
  const double p2     = pcm2(s, m1, m2);

  REQUIRE(p2 == Approx(e1 * e1 - m1 * m1).margin(1e-14));
  REQUIRE(p2 == Approx(pow2(DecayMomentum(sqrt_s, m1, m2))).margin(1e-14));
  REQUIRE(pcm2(s, m2, m1) == Approx(p2).margin(1e-14));
  REQUIRE(pcm2(s, 0.0, 0.0) == Approx(s / 4.0).margin(1e-14));
  REQUIRE(pcm2(pow2(m1 + m2), m1, m2) == Approx(0.0));
  REQUIRE(pcm2(1.0, 3.0, 1.0) == Approx(0.0));
}

TEST_CASE("DecayMomentum constructs an on-shell two-body decay", "[DecayMomentum][Kinematics]") {
  const double mother_mass = 10.0;
  const double m1          = 1.0;
  const double m2          = 2.0;
  const double p           = DecayMomentum(mother_mass, m1, m2);
  const double e1          = std::sqrt(m1 * m1 + p * p);
  const double e2          = std::sqrt(m2 * m2 + p * p);
  const M4Vec  daughter1(0.0, 0.0, p, e1);
  const M4Vec  daughter2(0.0, 0.0, -p, e2);

  REQUIRE(daughter1.M2() == Approx(m1 * m1).margin(1e-14));
  REQUIRE(daughter2.M2() == Approx(m2 * m2).margin(1e-14));
  REQUIRE((daughter1 + daughter2).M() == Approx(mother_mass).margin(1e-14));
  REQUIRE(p * p == Approx(pcm2(mother_mass * mother_mass, m1, m2)).margin(1e-14));
  REQUIRE(DecayMomentum(m1 + m2, m1, m2) == Approx(0.0));
  REQUIRE(DecayMomentum(m1 + m2 - 1e-6, m1, m2) == Approx(0.0));
  REQUIRE(DecayMomentum(-mother_mass, m1, m2) == Approx(0.0));
  REQUIRE(DecayMomentum(1.0, 3.0, 1.0) == Approx(0.0));
  REQUIRE(gra::kinematics::PS2Massive(1.0, 9.0, 1.0) == Approx(0.0));
  REQUIRE(gra::kinematics::PDW2body(1.0, 9.0, 1.0, 1.0, 1.0) == Approx(0.0));
}

TEST_CASE("LorentzBoost interfaces implement the same Lorentz transformation", "[LorentzBoost][Kinematics]") {
  M4Vec boost;
  boost.SetPxPyPzM(1.7, -0.8, 2.3, 3.1);
  const M4Vec particle(0.75, -0.42, 0.31, 1.75);

  // Compare all components of two public boost interfaces
  auto require_same_vector = [](const M4Vec &lhs, const M4Vec &rhs) {
    REQUIRE(lhs.Px() == Approx(rhs.Px()).margin(1e-13));
    REQUIRE(lhs.Py() == Approx(rhs.Py()).margin(1e-13));
    REQUIRE(lhs.Pz() == Approx(rhs.Pz()).margin(1e-13));
    REQUIRE(lhs.E() == Approx(rhs.E()).margin(1e-13));
  };

  SECTION("forward transformation") {
    const M4Vec member_result = particle.LorentzBoost(boost.BetaVector(), 1);
    M4Vec       free_result   = particle;
    LorentzBoost(boost, boost.M(), free_result, 1);

    require_same_vector(free_result, member_result);
    REQUIRE(free_result.M2() == Approx(particle.M2()).margin(1e-12));
  }

  SECTION("inverse transformation") {
    const M4Vec member_result = particle.LorentzBoost(boost.BetaVector(), -1);
    M4Vec       free_result   = particle;
    LorentzBoost(boost, boost.M(), free_result, -1);

    require_same_vector(free_result, member_result);
    REQUIRE(free_result.M2() == Approx(particle.M2()).margin(1e-12));
  }

  SECTION("successive inverse boosts recover the four-vector") {
    M4Vec transformed = particle;
    LorentzBoost(boost, boost.M(), transformed, 1);
    LorentzBoost(boost, boost.M(), transformed, -1);

    require_same_vector(transformed, particle);
  }

  SECTION("boost four-momentum transforms to its rest frame") {
    M4Vec system_rest = boost;
    LorentzBoost(boost, boost.M(), system_rest, -1);

    REQUIRE(system_rest.P3mod() == Approx(0.0).margin(1e-13));
    REQUIRE(system_rest.E() == Approx(boost.M()).margin(1e-13));
  }

  SECTION("positive direction transforms the rest system to the lab") {
    M4Vec system_lab(0.0, 0.0, 0.0, boost.M());
    LorentzBoost(boost, boost.M(), system_lab, 1);

    require_same_vector(system_lab, boost);
  }

  SECTION("both directions preserve pairwise Minkowski products") {
    const M4Vec  second(-0.37, 0.91, 0.48, 2.30);
    const double product = particle.DotM(second);
    for (const int sign : {-1, 1}) {
      M4Vec first_transformed  = particle;
      M4Vec second_transformed = second;
      LorentzBoost(boost, boost.M(), first_transformed, sign);
      LorentzBoost(boost, boost.M(), second_transformed, sign);
      REQUIRE(first_transformed.DotM(second_transformed) == Approx(product).margin(2e-12));
    }
  }
}

TEST_CASE(
    "LorentzBoost preserves mass shell and handles invalid boosts "
    "without throwing",
    "[LorentzBoost][Kinematics][Hardening]") {
  M4Vec boost;
  boost.SetPxPyPzM(1.2, -0.7, 2.1, 3.5);

  const M4Vec  particle(0.25, -0.33, 0.72, 1.20);
  M4Vec        roundtrip = particle;
  const double m2_before = roundtrip.M2();

  LorentzBoost(boost, boost.M(), roundtrip, -1);
  REQUIRE(roundtrip.M2() == Approx(m2_before).epsilon(1e-12));

  LorentzBoost(boost, boost.M(), roundtrip, 1);
  REQUIRE(roundtrip.Px() == Approx(particle.Px()).margin(1e-11));
  REQUIRE(roundtrip.Py() == Approx(particle.Py()).margin(1e-11));
  REQUIRE(roundtrip.Pz() == Approx(particle.Pz()).margin(1e-11));
  REQUIRE(roundtrip.E() == Approx(particle.E()).margin(1e-11));

  M4Vec system_rest = boost;
  LorentzBoost(boost, boost.M(), system_rest, -1);
  REQUIRE(system_rest.Px() == Approx(0.0).margin(1e-11));
  REQUIRE(system_rest.Py() == Approx(0.0).margin(1e-11));
  REQUIRE(system_rest.Pz() == Approx(0.0).margin(1e-11));
  REQUIRE(system_rest.E() == Approx(boost.M()).margin(1e-11));

  M4Vec       invalid = particle;
  const M4Vec lightlike_boost(1.0, 0.0, 0.0, 1.0);
  REQUIRE_NOTHROW(LorentzBoost(lightlike_boost, lightlike_boost.M(), invalid, -1));
  REQUIRE(invalid.E() == Approx(-1.0));

  const M4Vec superluminal = particle.LorentzBoost({1.01, 0.0, 0.0}, 1);
  REQUIRE(superluminal.E() == Approx(-1.0));
  const M4Vec invalid_member = particle.LorentzBoost({0.1, 0.0, 0.0}, 0);
  REQUIRE(invalid_member.E() == Approx(-1.0));

  invalid = particle;
  REQUIRE_NOTHROW(LorentzBoost(boost, boost.M(), invalid, 0));
  REQUIRE(invalid.E() == Approx(-1.0));
}

// Check the velocity API against the analytic longitudinal Lorentz transformation
TEST_CASE("M4Vec boosts preserve the supplied velocity for beta^2 < 1", "[LorentzBoost][M4Vec][Kinematics]") {
  const M4Vec rest(0.0, 0.0, 0.0, 1.0);
  for (const double beta : {0.0, 1.0e-12, 0.6, std::nextafter(1.0, 0.0)}) {
    for (const int sign : {-1, 1}) {
      CAPTURE(beta, sign);
      const long double speed = beta;
      const long double gamma = 1.0L / std::sqrt((1.0L - speed) * (1.0L + speed));
      const M4Vec result = rest.LorentzBoost({beta, 0.0, 0.0}, sign);
      REQUIRE(result.E() == Approx(static_cast<double>(gamma)).epsilon(2.0e-15));
      REQUIRE(result.Px() == Approx(static_cast<double>(sign * gamma * speed)).margin(1.0e-25).epsilon(2.0e-15));
      REQUIRE(result.Py() == Approx(0.0).margin(1.0e-25));
      REQUIRE(result.Pz() == Approx(0.0).margin(1.0e-25));
    }
  }
  for (const M3Vec beta : {M3Vec{1.0, 0.0, 0.0}, M3Vec{std::nextafter(1.0, 2.0), 0.0, 0.0},
                          M3Vec{std::numeric_limits<double>::quiet_NaN(), 0.0, 0.0}}) {
    REQUIRE(rest.LorentzBoost(beta).E() < 0.0);
  }
  const M4Vec nonfinite(0.0, 0.0, 0.0, std::numeric_limits<double>::infinity());
  REQUIRE(nonfinite.LorentzBoost({0.0, 0.0, 0.0}).E() < 0.0);
  const M4Vec enormous(0.0, 0.0, 0.0, std::numeric_limits<double>::max());
  REQUIRE(enormous.LorentzBoost({0.9, 0.0, 0.0}).E() < 0.0);
}

// Preserve the identity boost when its supplied mass differs only by rounding
TEST_CASE("LorentzBoost accepts mass-shell roundoff at rest", "[LorentzBoost][Kinematics][Hardening]") {
  const M4Vec particle(0.25, -0.33, 0.72, 1.20);
  for (const double mass : {1.0, 240.0, 1.0e6}) {
    const M4Vec rest(0.0, 0.0, 0.0, mass);
    for (const double supplied : {std::nextafter(mass, 0.0), std::nextafter(mass, std::numeric_limits<double>::infinity())}) {
      for (const int sign : {-1, 1}) {
        CAPTURE(mass, supplied, sign);
        M4Vec transformed = particle;
        LorentzBoost(rest, supplied, transformed, sign);
        REQUIRE(transformed.E() > 0.0);
        const auto actual = transformed.Contravariant<double>();
        const auto expected = particle.Contravariant<double>();
        for (const auto &i : indices(actual)) {
          REQUIRE(actual[i] == Approx(expected[i]).margin(1.0e-14));
        }
      }
    }
  }
  M4Vec invalid = particle;
  LorentzBoost(M4Vec(0.0, 0.0, 0.0, 1.0), 1.001, invalid, 1);
  REQUIRE(invalid.E() < 0.0);
}

TEST_CASE("LorentzBoost remains covariant for LHC scale non-collinear boosts",
          "[LorentzBoost][Kinematics][high-gamma]") {
  constexpr double system_mass = 0.9382720813;
  constexpr double momentum    = 6500.0;
  const double     bx          = 1700.0;
  const double     by          = -900.0;
  const double     bz          = std::sqrt(momentum * momentum - bx * bx - by * by);
  const double     energy      = std::sqrt(momentum * momentum + system_mass * system_mass);
  const M4Vec      boost(bx, by, bz, energy);

  const double proton_mass2 = gra::math::pow2(gra::PDG::mp);
  const double beam_energy  = std::sqrt(momentum * momentum + proton_mass2);
  const M4Vec  incoming(0.0, 0.0, momentum, beam_energy);
  M4Vec        outgoing;
  M4Vec        transfer;
  double       t = 0.0;
  REQUIRE(gra::kinematics::BuildEPATransfer(incoming, 1.0e-6, 2.0e-5, -3.0e-5, proton_mass2, proton_mass2, true,
                                            outgoing, transfer, t));

  M4Vec transformed_in  = incoming;
  M4Vec transformed_out = outgoing;
  M4Vec transformed_q   = transfer;
  LorentzBoost(boost, system_mass, transformed_in, -1);
  LorentzBoost(boost, system_mass, transformed_out, -1);
  LorentzBoost(boost, system_mass, transformed_q, -1);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    const double scale =
        std::max({1.0, std::abs(transformed_in[mu]), std::abs(transformed_out[mu]), std::abs(transformed_q[mu])});
    REQUIRE(std::abs(transformed_q[mu] - (transformed_in[mu] - transformed_out[mu])) < 2.0e-11 * scale);
  }

  const auto require_invariant = [](const M4Vec &before, const M4Vec &after) {
    const long double scale = gra::SquaredNorm(after.Contravariant<long double>());
    const long double error = std::abs(after.Invariant<long double>() - before.Invariant<long double>());
    REQUIRE(error < 16.0L * std::numeric_limits<double>::epsilon() * scale);
  };
  require_invariant(incoming, transformed_in);
  require_invariant(outgoing, transformed_out);
  require_invariant(transfer, transformed_q);
  const auto        transformed_in_ld  = transformed_in.Contravariant<long double>();
  const auto        transformed_out_ld = transformed_out.Contravariant<long double>();
  const long double inner_before =
      gra::MinkowskiProduct(incoming.Contravariant<long double>(), outgoing.Contravariant<long double>());
  const long double inner_after     = gra::MinkowskiProduct(transformed_in_ld, transformed_out_ld);
  long double       inner_condition = 0.0L;
  for (std::size_t mu = 0; mu < 4; ++mu) {
    inner_condition += std::abs(transformed_in_ld[mu] * transformed_out_ld[mu]);
  }
  REQUIRE(std::abs(inner_after - inner_before) < 16.0L * std::numeric_limits<double>::epsilon() * inner_condition);

  M4Vec roundtrip = transformed_q;
  LorentzBoost(boost, system_mass, roundtrip, 1);
  for (std::size_t mu = 0; mu < 4; ++mu) {
    const double scale = std::max({1.0, std::abs(roundtrip[mu]), std::abs(transfer[mu])});
    REQUIRE(std::abs(roundtrip[mu] - transfer[mu]) < 2.0e-8 * scale);
  }

  M4Vec inconsistent = transfer;
  LorentzBoost(boost, 1.01 * system_mass, inconsistent, -1);
  REQUIRE(inconsistent.E() == Approx(-1.0));
}

TEST_CASE("dPhi1 matches direct one-particle phase-space integration", "[Kinematics][PhaseSpace]") {
  const double max_energy = 2.5;
  const double mass       = 0.35;

  REQUIRE(gra::kinematics::dPhi1(max_energy, mass) == Approx(NumericalDPhi1(max_energy, mass)).epsilon(1e-10));
  REQUIRE(gra::kinematics::dPhi1(max_energy, -mass) == Approx(gra::kinematics::dPhi1(max_energy, mass)).epsilon(1e-15));
  REQUIRE(gra::kinematics::dPhi1(max_energy, 0.0) ==
          Approx(max_energy * max_energy / (8.0 * gra::math::pow2(gra::math::PI))).epsilon(1e-15));
  REQUIRE(gra::kinematics::dPhi1(mass, mass) == Approx(0.0));
  REQUIRE(gra::kinematics::dPhi1(0.1, mass) == Approx(0.0));
}

TEST_CASE("massive Mandelstam t and u use the correct outgoing legs", "[Kinematics][Mandelstam]") {
  const double sqrt_s   = 10.0;
  const double m1       = 1.0;
  const double m2       = 2.0;
  const double m3       = 1.5;
  const double m4       = 0.5;
  const double costheta = 0.37;
  const double sintheta = std::sqrt(1.0 - costheta * costheta);
  const double pin      = gra::kinematics::DecayMomentum(sqrt_s, m1, m2);
  const double pout     = gra::kinematics::DecayMomentum(sqrt_s, m3, m4);
  const double e1       = std::sqrt(m1 * m1 + pin * pin);
  const double e2       = std::sqrt(m2 * m2 + pin * pin);
  const double e3       = std::sqrt(m3 * m3 + pout * pout);
  const double e4       = std::sqrt(m4 * m4 + pout * pout);

  const M4Vec  p1(0.0, 0.0, pin, e1);
  const M4Vec  p2(0.0, 0.0, -pin, e2);
  const M4Vec  p3(pout * sintheta, 0.0, pout * costheta, e3);
  const M4Vec  p4(-pout * sintheta, 0.0, -pout * costheta, e4);
  const double t_expected = m1 * m1 + m3 * m3 - 2.0 * (e1 * e3 - pin * pout * costheta);
  const double u_expected = m1 * m1 + m4 * m4 - 2.0 * (e1 * e4 + pin * pout * costheta);

  const double s = (p1 + p2).M2();
  const double t = gra::kinematics::mandelstam_t(p1, p3);
  const double u = gra::kinematics::mandelstam_u(p1, p4);
  REQUIRE(t == Approx(t_expected).margin(1e-12));
  REQUIRE(u == Approx(u_expected).margin(1e-12));
  REQUIRE(s + t + u == Approx(m1 * m1 + m2 * m2 + m3 * m3 + m4 * m4).margin(1e-11));
  REQUIRE(gra::math::CheckEMC(p1 + p2 - p3 - p4, 1e-12));
}

TEST_CASE("longitudinal momentum loss uses exact beam-side light-cone fractions", "[Kinematics][LightCone]") {
  const double mass = 0.938;
  M4Vec        plus_beam;
  M4Vec        minus_beam;
  plus_beam.SetPxPyPzM(0.0, 0.0, 2.0, mass);
  minus_beam.SetPxPyPzM(0.0, 0.0, -2.0, mass);
  M4Vec        plus_forward;
  M4Vec        minus_forward;
  const double xi_plus  = 0.17;
  const double xi_minus = 0.23;
  REQUIRE(gra::kinematics::BuildForwardParticle(plus_beam, 1.0 - xi_plus, 0.31, -0.22, true, plus_forward));
  REQUIRE(gra::kinematics::BuildForwardParticle(minus_beam, 1.0 - xi_minus, -0.27, 0.19, false, minus_forward));
  REQUIRE(gra::kinematics::LongitudinalMomentumLoss(plus_beam, plus_forward, true) == Approx(xi_plus).margin(2e-15));
  REQUIRE(gra::kinematics::LongitudinalMomentumLoss(minus_beam, minus_forward, false) ==
          Approx(xi_minus).margin(2e-15));
  REQUIRE(std::abs((1.0 - plus_forward.Pz() / plus_beam.Pz()) - xi_plus) > 1e-3);
  REQUIRE(std::abs((1.0 - minus_forward.Pz() / minus_beam.Pz()) - xi_minus) > 1e-3);

  M4Vec high_energy_beam;
  high_energy_beam.SetPxPyPzM(0.0, 0.0, 1.0e6, mass);
  M4Vec high_energy_forward;
  REQUIRE(
      gra::kinematics::BuildForwardParticle(high_energy_beam, 1.0 - xi_plus, 0.31, -0.22, true, high_energy_forward));
  const double high_energy_approximation = 1.0 - high_energy_forward.Pz() / high_energy_beam.Pz();
  REQUIRE(high_energy_approximation == Approx(xi_plus).margin(1e-11));
  REQUIRE_THROWS_AS(gra::kinematics::LongitudinalMomentumLoss(M4Vec(), plus_forward, true), std::invalid_argument);
}

TEST_CASE("SolvePz returns physical longitudinal momenta or an invalid sentinel", "[Kinematics][Hardening]") {
  auto check_solution = [](double m3, double m4, double pt3, double pt4, double pz3_seed, double pz4_seed, double E5) {
    const double pz5    = -pz3_seed - pz4_seed;
    const double E3seed = std::sqrt(m3 * m3 + pt3 * pt3 + pz3_seed * pz3_seed);
    const double E4seed = std::sqrt(m4 * m4 + pt4 * pt4 + pz4_seed * pz4_seed);
    const double sqrts  = E3seed + E4seed + E5;
    const double pz3    = gra::kinematics::SolvePz(m3, m4, pt3, pt4, pz5, E5, sqrts * sqrts);
    REQUIRE(pz3 != Approx(-1.0));
    REQUIRE(std::isfinite(pz3));

    const double pz4 = -pz3 - pz5;
    const double E3  = std::sqrt(m3 * m3 + pt3 * pt3 + pz3 * pz3);
    const double E4  = std::sqrt(m4 * m4 + pt4 * pt4 + pz4 * pz4);
    REQUIRE(E3 + E4 + E5 == Approx(sqrts).margin(2e-10));
    REQUIRE(std::sqrt(E3 * E3 - pz3 * pz3 - pt3 * pt3) == Approx(m3).margin(1e-10));
    REQUIRE(std::sqrt(E4 * E4 - pz4 * pz4 - pt4 * pt4) == Approx(m4).margin(1e-10));
  };

  check_solution(0.938, 0.938, 0.40, 0.70, 1.10, -0.35, 2.0);
  check_solution(0.938, 1.300, 0.40, 0.70, 1.10, -0.35, 2.0);

  REQUIRE(gra::kinematics::SolvePz(0.938, 0.938, 0.4, 0.7, 0.2, 2.0, -1.0) == Approx(-1.0));
  REQUIRE(gra::kinematics::SolvePz(0.938, 1.300, 0.4, 0.7, 0.2, 20.0, 1.0) == Approx(-1.0));
}

TEST_CASE(
    "MKinematics phase-space guards and longitudinal solution branches "
    "are physical",
    "[Kinematics][Coverage]") {
  const MCW invalid = gra::kinematics::InvalidPhaseSpacePoint();
  REQUIRE(invalid.GetW() == Approx(-1.0));
  REQUIRE(invalid.GetW2() == Approx(0.0));
  REQUIRE(invalid.GetN() == Approx(0.0));

  REQUIRE(gra::kinematics::HasPhaseSpace(4.0, {0.7, 1.1, 1.3}));
  REQUIRE_FALSE(gra::kinematics::HasPhaseSpace(3.1, {0.7, 1.1, 1.3}));
  REQUIRE_FALSE(gra::kinematics::HasPhaseSpace(4.0, {0.7, -0.1, 1.3}));
  REQUIRE_FALSE(gra::kinematics::HasPhaseSpace(std::numeric_limits<double>::infinity(), {0.7, 1.1}));

  const double pt3 = 0.43;
  const double pt4 = 0.68;
  const double pz5 = -0.37;
  const double E5  = 2.2;

  SECTION("equal-mass analytic branch") {
    const double m3       = 0.938;
    const double m4       = m3;
    const double pz3_seed = 1.14;
    const double pz4_seed = -pz3_seed - pz5;
    const double sqrt_s   = std::sqrt(m3 * m3 + pt3 * pt3 + pz3_seed * pz3_seed) +
                          std::sqrt(m4 * m4 + pt4 * pt4 + pz4_seed * pz4_seed) + E5;
    const double s        = sqrt_s * sqrt_s;
    const double solution = gra::kinematics::SolvePz3_A(m3, pt3, pt4, pz5, E5, s);

    REQUIRE(gra::kinematics::SolvePzInputHasPhaseSpace(m3, m4, pt3, pt4, pz5, E5, s));
    REQUIRE(gra::kinematics::SolvePzSolutionConservesEnergy(m3, m4, pt3, pt4, pz5, E5, s, solution));
    REQUIRE(gra::kinematics::SolvePz(m3, m4, pt3, pt4, pz5, E5, s) == Approx(solution).margin(1e-12));
  }

  SECTION("unequal-mass analytic branch") {
    const double m3       = 0.938;
    const double m4       = 1.42;
    const double pz3_seed = 1.14;
    const double pz4_seed = -pz3_seed - pz5;
    const double sqrt_s   = std::sqrt(m3 * m3 + pt3 * pt3 + pz3_seed * pz3_seed) +
                          std::sqrt(m4 * m4 + pt4 * pt4 + pz4_seed * pz4_seed) + E5;
    const double s        = sqrt_s * sqrt_s;
    const double solution = gra::kinematics::SolvePz3_B(m3, m4, pt3, pt4, pz5, E5, s);

    REQUIRE(gra::kinematics::SolvePzInputHasPhaseSpace(m3, m4, pt3, pt4, pz5, E5, s));
    REQUIRE(gra::kinematics::SolvePzSolutionConservesEnergy(m3, m4, pt3, pt4, pz5, E5, s, solution));
    REQUIRE(gra::kinematics::SolvePz(m3, m4, pt3, pt4, pz5, E5, s) == Approx(solution).margin(1e-12));
  }

  REQUIRE_FALSE(gra::kinematics::SolvePzInputHasPhaseSpace(0.938, 1.42, pt3, pt4, pz5, E5, 1.0));
  REQUIRE_FALSE(gra::kinematics::SolvePzSolutionConservesEnergy(0.938, 1.42, pt3, pt4, pz5, E5, 100.0,
                                                                std::numeric_limits<double>::quiet_NaN()));
}

TEST_CASE("forward xi-t construction closes exact light-cone and transfer kinematics",
          "[Kinematics][Forward][Coverage]") {
  const double mass = 0.938272;
  const double xi   = 0.19;
  const double t    = -0.73;
  const double phi  = 1.17;

  for (const bool plus_side : {true, false}) {
    M4Vec beam;
    beam.SetPxPyPzM(0.0, 0.0, plus_side ? 15.0 : -15.0, mass);
    M4Vec forward;
    REQUIRE(gra::kinematics::BuildForwardParticleXiT(beam, xi, t, phi, plus_side, forward));
    REQUIRE(forward.M2() == Approx(mass * mass).margin(1e-12));
    REQUIRE((beam - forward).M2() == Approx(t).margin(1e-12));
    REQUIRE(gra::kinematics::LongitudinalMomentumLoss(beam, forward, plus_side) == Approx(xi).margin(1e-14));
    REQUIRE(forward.Phi() == Approx(phi).margin(1e-14));
  }

  M4Vec beam;
  beam.SetPxPyPzM(0.0, 0.0, 15.0, mass);
  M4Vec forward;
  REQUIRE_FALSE(gra::kinematics::BuildForwardParticleXiT(beam, 0.0, t, phi, true, forward));
  REQUIRE_FALSE(gra::kinematics::BuildForwardParticleXiT(beam, xi, 0.1, phi, true, forward));
  REQUIRE_FALSE(gra::kinematics::BuildForwardParticleXiT(beam, xi, -0.01, phi, true, forward));
}

TEST_CASE("closed-form two-body helpers agree with explicit massive scattering", "[Kinematics][TwoBody][Coverage]") {
  const double sqrt_s   = 9.0;
  const double s        = sqrt_s * sqrt_s;
  const double m1       = 0.8;
  const double m2       = 1.3;
  const double m3       = 1.1;
  const double m4       = 0.6;
  const double costheta = -0.41;
  const double sintheta = std::sqrt(1.0 - costheta * costheta);
  const double pin      = gra::kinematics::DecayMomentum(sqrt_s, m1, m2);
  const double pout     = gra::kinematics::DecayMomentum(sqrt_s, m3, m4);
  const double e1       = std::sqrt(m1 * m1 + pin * pin);
  const double e3       = std::sqrt(m3 * m3 + pout * pout);
  const M4Vec  p1(0.0, 0.0, pin, e1);
  const M4Vec  p3(pout * sintheta, 0.0, pout * costheta, e3);
  const double t = (p1 - p3).M2();

  const double lambda_in = gra::kinematics::SqrtKallenLambda(s, m1 * m1, m2 * m2);
  REQUIRE(lambda_in == Approx(2.0 * sqrt_s * pin).margin(1e-12));
  REQUIRE(gra::kinematics::beta12(s, m1, m2) == Approx(2.0 * pin / sqrt_s).margin(1e-12));
  REQUIRE(gra::kinematics::dPhi2(sqrt_s, pout) == Approx(pout / (4.0 * gra::math::PI * sqrt_s)).margin(1e-15));
  REQUIRE(gra::kinematics::PS2Massive(s, m3 * m3, m4 * m4) ==
          Approx(pout / (4.0 * gra::math::PI * sqrt_s)).margin(1e-15));
  REQUIRE(gra::kinematics::CosthetaStar(s, t, m1 * m1, m2 * m2, m3 * m3, m4 * m4) == Approx(costheta).margin(1e-12));

  double tmin = 0.0;
  double tmax = 0.0;
  gra::kinematics::Two2TwoLimit(s, m1 * m1, m2 * m2, m3 * m3, m4 * m4, tmin, tmax);
  const double t_at_minus_one = m1 * m1 + m3 * m3 - 2.0 * e1 * e3 - 2.0 * pin * pout;
  const double t_at_plus_one  = m1 * m1 + m3 * m3 - 2.0 * e1 * e3 + 2.0 * pin * pout;
  REQUIRE(tmin == Approx(t_at_minus_one).margin(1e-12));
  REQUIRE(tmax == Approx(t_at_plus_one).margin(1e-12));

  const double daughter_m2 = 0.49;
  const double shat        = 16.0;
  REQUIRE(gra::kinematics::Beta(daughter_m2, shat) == Approx(std::sqrt(1.0 - 4.0 * daughter_m2 / shat)).margin(1e-15));

  const double matrix_element2 = 2.7;
  const double symmetry_factor = 2.0;
  const double expected_width  = gra::kinematics::SqrtKallenLambda(s, m3 * m3, m4 * m4) * matrix_element2 /
                                (16.0 * gra::math::PI * symmetry_factor * std::pow(sqrt_s, 3));
  REQUIRE(gra::kinematics::PDW2body(s, m3 * m3, m4 * m4, matrix_element2, symmetry_factor) ==
          Approx(expected_width).margin(1e-15));

  REQUIRE(gra::kinematics::dPhi2(0.0, pout) == Approx(0.0));
  REQUIRE(gra::kinematics::beta12(0.0, m1, m2) == Approx(0.0));
  REQUIRE(gra::kinematics::Beta(daughter_m2, 0.0) == Approx(0.0));
}

TEST_CASE("isotropic samplers produce non-trivial on-shell back-to-back momenta", "[Kinematics][Isotropic][Coverage]") {
  MRandom rng;
  rng.SetSeed(72531);

  double costheta = 0.0;
  double sintheta = 0.0;
  double phi      = 0.0;
  gra::kinematics::FlatIsotropic(costheta, sintheta, phi, rng);
  REQUIRE(costheta >= -1.0);
  REQUIRE(costheta <= 1.0);
  REQUIRE(phi >= 0.0);
  REQUIRE(phi <= 2.0 * gra::math::PI);
  REQUIRE(costheta * costheta + sintheta * sintheta == Approx(1.0).margin(1e-15));

  const double pnorm = 1.37;
  const double m1    = 0.31;
  const double m2    = 0.77;
  M4Vec        p1;
  M4Vec        p2;
  gra::kinematics::Isotropic(pnorm, p1, p2, m1, m2, rng);
  REQUIRE((p1 + p2).P3mod() == Approx(0.0).margin(1e-14));
  REQUIRE(p1.P3mod() == Approx(pnorm).margin(1e-14));
  REQUIRE(p2.P3mod() == Approx(pnorm).margin(1e-14));
  REQUIRE(p1.M2() == Approx(m1 * m1).margin(1e-14));
  REQUIRE(p2.M2() == Approx(m2 * m2).margin(1e-14));
  REQUIRE(std::abs(p1.Px()) + std::abs(p1.Py()) + std::abs(p1.Pz()) > 0.5);
}

TEST_CASE("NBodySetup generates ordered physical intermediate masses", "[Kinematics][NBody][Coverage]") {
  const double              mother_mass = 8.0;
  const M4Vec               mother(0.0, 0.0, 0.0, mother_mass);
  const std::vector<double> masses = {0.4, 0.7, 1.1, 0.9};
  std::vector<double>       effective_masses(masses.size(), 0.0);
  std::vector<double>       decay_momenta;
  effective_masses.clear();
  MRandom rng;
  rng.SetSeed(99173);

  const MCW weight =
      gra::kinematics::NBodySetup(mother, mother_mass, masses, effective_masses, decay_momenta, false, rng);
  REQUIRE(weight.GetW() > 0.0);
  REQUIRE(weight.GetN() == Approx(1.0));
  REQUIRE(effective_masses.front() == Approx(masses.front()).margin(1e-15));
  REQUIRE(effective_masses.back() == Approx(mother_mass).margin(1e-15));

  for (std::size_t i = 1; i < masses.size(); ++i) {
    REQUIRE(effective_masses[i] > effective_masses[i - 1] + masses[i]);
    REQUIRE(
        decay_momenta[i] ==
        Approx(gra::kinematics::DecayMomentum(effective_masses[i], effective_masses[i - 1], masses[i])).margin(1e-14));
  }
}

TEST_CASE("generic spatial rotation follows its stated theta-phi matrix", "[Kinematics][Rotation][Coverage]") {
  const double theta = 0.71;
  const double phi   = -0.46;
  const M4Vec  input(0.8, -1.1, 2.3, 3.4);
  M4Vec        rotated = input;
  rotated.Rotate(theta, phi);

  const double c1 = std::cos(theta);
  const double s1 = std::sin(theta);
  const double c2 = std::cos(phi);
  const double s2 = std::sin(phi);
  REQUIRE(rotated.Px() == Approx(c1 * c2 * input.Px() - s2 * input.Py() + s1 * c2 * input.Pz()).margin(1e-14));
  REQUIRE(rotated.Py() == Approx(c1 * s2 * input.Px() + c2 * input.Py() + s1 * s2 * input.Pz()).margin(1e-14));
  REQUIRE(rotated.Pz() == Approx(-s1 * input.Px() + c1 * input.Pz()).margin(1e-14));
  REQUIRE(rotated.E() == Approx(input.E()).margin(1e-15));
  REQUIRE(rotated.M2() == Approx(input.M2()).margin(1e-13));
}

TEST_CASE("M4Vec applies spatial rotation matrices without changing energy", "[Kinematics][Rotation][M4Vec]") {
  const MMatrix<double> rotation = {{0.0, -1.0, 0.0}, {1.0, 0.0, 0.0}, {0.0, 0.0, 1.0}};
  M4Vec                 vector(2.0, -3.0, 4.0, 9.0);

  vector.Rotate(rotation);

  REQUIRE(vector.Px() == Approx(3.0).margin(1e-15));
  REQUIRE(vector.Py() == Approx(2.0).margin(1e-15));
  REQUIRE(vector.Pz() == Approx(4.0).margin(1e-15));
  REQUIRE(vector.E() == Approx(9.0).margin(0.0));
}

TEST_CASE("axis rotations satisfy SO(3) matrices and inverse composition", "[RotateX][RotateY][RotateZ][Kinematics]") {
  const M4Vec  input(0.7, -1.2, 2.4, 3.2);
  const double angle  = 0.63;
  const double cosine = std::cos(angle);
  const double sine   = std::sin(angle);

  // Check four-vector equality after an inverse rotation
  auto require_original = [&](const M4Vec &result) {
    REQUIRE(result.Px() == Approx(input.Px()).margin(1e-14));
    REQUIRE(result.Py() == Approx(input.Py()).margin(1e-14));
    REQUIRE(result.Pz() == Approx(input.Pz()).margin(1e-14));
    REQUIRE(result.E() == Approx(input.E()).margin(1e-14));
  };

  SECTION("x-axis") {
    M4Vec rotated = input;
    rotated.RotateX(angle);
    REQUIRE(rotated.Px() == Approx(input.Px()).margin(1e-14));
    REQUIRE(rotated.Py() == Approx(cosine * input.Py() - sine * input.Pz()).margin(1e-14));
    REQUIRE(rotated.Pz() == Approx(sine * input.Py() + cosine * input.Pz()).margin(1e-14));
    REQUIRE(rotated.M2() == Approx(input.M2()).margin(1e-13));
    rotated.RotateX(-angle);
    require_original(rotated);
  }

  SECTION("y-axis") {
    M4Vec rotated = input;
    rotated.RotateY(angle);
    REQUIRE(rotated.Px() == Approx(cosine * input.Px() + sine * input.Pz()).margin(1e-14));
    REQUIRE(rotated.Py() == Approx(input.Py()).margin(1e-14));
    REQUIRE(rotated.Pz() == Approx(-sine * input.Px() + cosine * input.Pz()).margin(1e-14));
    REQUIRE(rotated.M2() == Approx(input.M2()).margin(1e-13));
    rotated.RotateY(-angle);
    require_original(rotated);
  }

  SECTION("z-axis") {
    M4Vec rotated = input;
    rotated.RotateZ(angle);
    REQUIRE(rotated.Px() == Approx(cosine * input.Px() - sine * input.Py()).margin(1e-14));
    REQUIRE(rotated.Py() == Approx(sine * input.Px() + cosine * input.Py()).margin(1e-14));
    REQUIRE(rotated.Pz() == Approx(input.Pz()).margin(1e-14));
    REQUIRE(rotated.M2() == Approx(input.M2()).margin(1e-13));
    rotated.RotateZ(-angle);
    require_original(rotated);
  }
}

struct PhaseSpaceEstimatesForTest {
  gra::kinematics::MCW dedicated;
  gra::kinematics::MCW nbody;
  gra::kinematics::MCW rambo;
};

// Sample phase space and verify every generated final state
PhaseSpaceEstimatesForTest EvaluatePhaseSpaceForTest(unsigned int multiplicity, const M4Vec &mother,
                                                     double daughter_mass, std::uint64_t seed, unsigned int samples) {
  MRandom rng;
  rng.SetSeed(seed);

  const double               mother_mass = mother.M();
  const std::vector<double>  masses(multiplicity, daughter_mass);
  std::vector<M4Vec>         momenta(multiplicity);
  PhaseSpaceEstimatesForTest estimates;

  // Check conservation, mass shells, and one valid Monte Carlo weight
  auto check_point = [&](const gra::kinematics::MCW &weight, const char *algorithm, unsigned int sample,
                         double conservation_tolerance) {
    CAPTURE(algorithm, sample, multiplicity, daughter_mass);
    REQUIRE(weight.GetN() == Approx(1.0));
    REQUIRE(std::isfinite(weight.GetW()));
    REQUIRE(weight.GetW() >= 0.0);
    const M4Vec residual = mother - psum(momenta);
    CAPTURE(residual.Px(), residual.Py(), residual.Pz(), residual.E());
    REQUIRE(math::CheckEMC(residual, conservation_tolerance));
    for (const auto &momentum : momenta) {
      REQUIRE(momentum.M2() == Approx(daughter_mass * daughter_mass).margin(1e-8));
      REQUIRE(momentum.E() > 0.0);
    }
  };

  for (unsigned int sample = 0; sample < samples; ++sample) {
    if (multiplicity == 2) {
      const auto weight = gra::kinematics::TwoBodyPhaseSpace(mother, mother_mass, masses, momenta, rng);
      check_point(weight, "TwoBodyPhaseSpace", sample, 1e-9);
      estimates.dedicated += weight;
    } else if (multiplicity == 3) {
      const auto weight = gra::kinematics::ThreeBodyPhaseSpace(mother, mother_mass, masses, momenta, false, rng);
      check_point(weight, "ThreeBodyPhaseSpace", sample, 1e-9);
      estimates.dedicated += weight;
    }

    const auto nbody_weight = gra::kinematics::NBodyPhaseSpace(mother, mother_mass, masses, momenta, false, rng);
    check_point(nbody_weight, "NBodyPhaseSpace", sample, 1e-9);
    estimates.nbody += nbody_weight;

    const auto   rambo_weight = gra::kinematics::RamboMassive(mother, mother_mass, masses, momenta, rng);
    const double rambo_tolerance = 1e-9 * std::max(1.0, std::abs(mother.E() / mother_mass));
    check_point(rambo_weight, "RamboMassive", sample, rambo_tolerance);
    estimates.rambo += rambo_weight;
  }
  return estimates;
}

// Require two independent Monte Carlo estimates to agree within uncertainty
void RequirePhaseSpaceCompatibilityForTest(const gra::kinematics::MCW &lhs, const gra::kinematics::MCW &rhs) {
  const double combined_error = std::hypot(lhs.IntegralError(), rhs.IntegralError());
  const double roundoff       = 100.0 * std::numeric_limits<double>::epsilon() *
                          std::max({1.0, std::abs(lhs.Integral()), std::abs(rhs.Integral())});
  REQUIRE(std::abs(lhs.Integral() - rhs.Integral()) <= 6.0 * combined_error + roundoff);
}

// Require one Monte Carlo estimate to agree with an analytic phase-space volume
void RequirePhaseSpaceNormalizationForTest(const gra::kinematics::MCW &estimate, double exact) {
  const double roundoff = 100.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(exact));
  REQUIRE(std::abs(estimate.Integral() - exact) <= 6.0 * estimate.IntegralError() + roundoff);
}

TEST_CASE(
    "phase-space generators conserve four-momentum and reproduce "
    "invariant volumes",
    "[PhaseSpace][Kinematics]") {
  constexpr double       mother_mass = 8.0;
  constexpr unsigned int samples     = 10000;
  M4Vec                  mother_rest(0.0, 0.0, 0.0, mother_mass);
  M4Vec                  mother_lab;
  mother_lab.SetPxPyPzM(1.7, -0.9, 2.4, mother_mass);

  for (const unsigned int multiplicity : {2U, 3U, 4U, 5U}) {
    for (const double daughter_mass : {0.0, 0.139}) {
      CAPTURE(multiplicity, daughter_mass);
      const std::uint64_t seed = 123456U + 100U * multiplicity + static_cast<std::uint64_t>(1000.0 * daughter_mass);
      const auto          rest = EvaluatePhaseSpaceForTest(multiplicity, mother_rest, daughter_mass, seed, samples);
      const auto          lab  = EvaluatePhaseSpaceForTest(multiplicity, mother_lab, daughter_mass, seed, samples);

      REQUIRE(lab.nbody.Integral() == Approx(rest.nbody.Integral()).margin(1e-13));
      REQUIRE(lab.rambo.Integral() == Approx(rest.rambo.Integral()).margin(1e-13));
      if (multiplicity <= 3) { REQUIRE(lab.dedicated.Integral() == Approx(rest.dedicated.Integral()).margin(1e-13)); }

      if (gra::math::IsZero(daughter_mass)) {
        const double exact = gra::kinematics::PSnMassless(mother_mass * mother_mass, multiplicity);
        RequirePhaseSpaceNormalizationForTest(rest.nbody, exact);
        RequirePhaseSpaceNormalizationForTest(rest.rambo, exact);
        if (multiplicity <= 3) { RequirePhaseSpaceNormalizationForTest(rest.dedicated, exact); }
      } else if (multiplicity == 2) {
        const double exact = gra::kinematics::PS2Massive(mother_mass * mother_mass, daughter_mass * daughter_mass,
                                                         daughter_mass * daughter_mass);
        RequirePhaseSpaceNormalizationForTest(rest.dedicated, exact);
        RequirePhaseSpaceNormalizationForTest(rest.nbody, exact);
        RequirePhaseSpaceNormalizationForTest(rest.rambo, exact);
      } else if (multiplicity == 3) {
        RequirePhaseSpaceCompatibilityForTest(rest.dedicated, rest.nbody);
        RequirePhaseSpaceCompatibilityForTest(rest.dedicated, rest.rambo);
      } else {
        RequirePhaseSpaceCompatibilityForTest(rest.nbody, rest.rambo);
      }
    }
  }
}

TEST_CASE("logarithmic central phase space reproduces flat RAMBO", "[PhaseSpace][Continuum]") {
  MRandom rng;
  rng.SetSeed(987654);

  ContinuumPatchRamboTest(2, rng);
  ContinuumPatchRamboTest(3, rng);
  ContinuumPatchRamboTest(4, rng);
  ContinuumPatchRamboTest(5, rng);
}

TEST_CASE("continuum transverse recurrence closes every supported topology", "[Kinematics][Continuum]") {
  const M4Vec p1f(0.31, -0.27, 0.0, 0.0);
  const M4Vec p2f(-0.19, 0.43, 0.0, 0.0);

  for (std::size_t multiplicity = 2; multiplicity <= 8; ++multiplicity) {
    std::vector<M4Vec> differences;
    for (std::size_t i = 0; i + 1 < multiplicity; ++i) {
      differences.emplace_back(0.07 * static_cast<double>(i + 1), -0.11 * static_cast<double>(i + 2), 0.0, 0.0);
    }

    std::vector<M4Vec> central;
    REQUIRE(gra::kinematics::BuildCentralTransverseMomenta(p1f, p2f, differences, central));
    REQUIRE(central.size() == multiplicity);
    for (std::size_t i = 0; i < differences.size(); ++i) {
      REQUIRE((central[i] - central[i + 1]).Px() == Approx(differences[i].Px()).margin(1e-14));
      REQUIRE((central[i] - central[i + 1]).Py() == Approx(differences[i].Py()).margin(1e-14));
    }
    M4Vec sum;
    for (const auto &momentum : central) { sum += momentum; }
    REQUIRE(sum.Px() == Approx(-(p1f + p2f).Px()).margin(1e-13));
    REQUIRE(sum.Py() == Approx(-(p1f + p2f).Py()).margin(1e-13));
  }

  std::vector<M4Vec> invalid;
  REQUIRE_FALSE(gra::kinematics::BuildCentralTransverseMomenta(p1f, p2f, {}, invalid));
}

TEST_CASE(
    "James massless phase-space volumes normalize weighted and "
    "unweighted generation",
    "[PhaseSpace][Literature]") {
  // F James, CERN 68-15, https://cds.cern.ch/record/275743/files/CERN-68-15.pdf
  // The tabulated result is Phi_n(s)=s^(n-2)/[2(4pi)^(2n-3)(n-1)!(n-2)!]
  const double mass = 7.0;
  const M4Vec  mother(0.0, 0.0, 0.0, mass);
  MRandom      rng;
  rng.SetSeed(24681357);

  for (const std::size_t multiplicity : {std::size_t{3}, std::size_t{4}}) {
    gra::kinematics::MCW      aggregate;
    std::vector<M4Vec>        daughters;
    const std::vector<double> masses(multiplicity, 0.0);
    for (std::size_t event = 0; event < 50000; ++event) {
      const auto weight = multiplicity == 3
                              ? gra::kinematics::ThreeBodyPhaseSpace(mother, mass, masses, daughters, false, rng)
                              : gra::kinematics::NBodyPhaseSpace(mother, mass, masses, daughters, false, rng);
      REQUIRE(weight.GetN() == Approx(1.0));
      const M4Vec residual = mother - psum(daughters);
      CAPTURE(event, multiplicity, residual.Px(), residual.Py(), residual.Pz(), residual.E());
      REQUIRE(gra::math::CheckEMC(residual, 1e-10));
      aggregate += weight;
    }
    const double s     = mass * mass;
    const double exact = multiplicity == 3 ? s / (256.0 * gra::math::pow3(gra::math::PI))
                                           : s * s / (24.0 * std::pow(4.0 * gra::math::PI, 5));
    REQUIRE(gra::kinematics::PSnMassless(s, multiplicity) == Approx(exact).epsilon(1e-14));
    REQUIRE(std::abs(aggregate.Integral() - exact) < std::max(6.0 * aggregate.IntegralError(), 0.005 * exact));
  }

  gra::kinematics::MCW unweighted_aggregate;
  std::vector<M4Vec>   daughters;
  for (std::size_t event = 0; event < 30000; ++event) {
    unweighted_aggregate += gra::kinematics::ThreeBodyPhaseSpace(mother, mass, {0.0, 0.0, 0.0}, daughters, true, rng);
  }
  const double exact = mass * mass / (256.0 * gra::math::pow3(gra::math::PI));
  REQUIRE(unweighted_aggregate.GetN() > 30000.0);
  REQUIRE(std::abs(unweighted_aggregate.Integral() - exact) <
          std::max(6.0 * unweighted_aggregate.IntegralError(), 0.007 * exact));
}

TEST_CASE("gra::kinematics phase space generators reject threshold violations", "[PhaseSpace][Hardening]") {
  MRandom rng;
  rng.SetSeed(123456);

  const M4Vec mother(0.0, 0.0, 0.0, 0.25);

  SECTION("Two-body threshold") {
    std::vector<M4Vec> p;
    const MCW          x = gra::kinematics::TwoBodyPhaseSpace(mother, mother.M(), {0.14, 0.14}, p, rng);
    REQUIRE(x.GetW() == Approx(-1.0));
    REQUIRE(x.GetN() == Approx(0.0));
  }

  SECTION("Three-body threshold") {
    std::vector<M4Vec> p;
    const MCW          x = gra::kinematics::ThreeBodyPhaseSpace(mother, mother.M(), {0.1, 0.1, 0.1}, p, false, rng);
    REQUIRE(x.GetW() == Approx(-1.0));
    REQUIRE(x.GetN() == Approx(0.0));
  }

  SECTION("N-body threshold") {
    std::vector<M4Vec> p;
    const MCW          x = gra::kinematics::NBodyPhaseSpace(mother, mother.M(), {0.1, 0.1, 0.1}, p, false, rng);
    REQUIRE(x.GetW() == Approx(-1.0));
    REQUIRE(x.GetN() == Approx(0.0));
  }

  SECTION("RamboMassive threshold") {
    std::vector<M4Vec> p;
    const MCW          x = gra::kinematics::RamboMassive(mother, mother.M(), {0.14, 0.14}, p, rng);
    REQUIRE(x.GetW() == Approx(-1.0));
    REQUIRE(x.GetN() == Approx(0.0));
  }
}

TEST_CASE("gra::kinematics singular helpers stay finite", "[Kinematics][Hardening]") {
  SECTION("RotationTo is an orthogonal 3D rotation for generic axes") {
    const gra::M3Vec a = {0.3, -0.4, 0.8660254037844386};
    const gra::M3Vec b = {-0.2, 0.9, 0.3872983346207417};
    const M4Vec      source(a[0], a[1], a[2], 7.0);

    const MMatrix<double> R  = source.RotationTo(b);
    const gra::M3Vec      ab = R * a;

    REQUIRE(ab[0] == Approx(b[0]).margin(1e-12));
    REQUIRE(ab[1] == Approx(b[1]).margin(1e-12));
    REQUIRE(ab[2] == Approx(b[2]).margin(1e-12));

    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) {
        double dot = 0.0;
        for (std::size_t k = 0; k < 3; ++k) { dot += R[k][i] * R[k][j]; }
        REQUIRE(dot == Approx(i == j ? 1.0 : 0.0).margin(1e-12));
      }
    }
  }

  SECTION("RotationTo handles anti-parallel vectors") {
    const gra::M3Vec a = {0.0, 0.0, 1.0};
    const gra::M3Vec b = {0.0, 0.0, -1.0};
    M4Vec            source(a[0], a[1], a[2], 7.0);

    const MMatrix<double> R  = source.RotationTo(b);
    const gra::M3Vec      ab = R * a;

    REQUIRE(ab[0] == Approx(b[0]).margin(1e-12));
    REQUIRE(ab[1] == Approx(b[1]).margin(1e-12));
    REQUIRE(ab[2] == Approx(b[2]).margin(1e-12));

    source.RotateTo(b);
    REQUIRE(source.Px() == Approx(b[0]).margin(1e-12));
    REQUIRE(source.Py() == Approx(b[1]).margin(1e-12));
    REQUIRE(source.Pz() == Approx(b[2]).margin(1e-12));
    REQUIRE(source.E() == Approx(7.0).margin(0.0));
  }

  SECTION("RotationTo keeps parallel directions unchanged") {
    M4Vec                 source(0.0, 0.0, 2.0, 7.0);
    const MMatrix<double> rotation = source.RotationTo({0.0, 0.0, 5.0});

    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) { REQUIRE(rotation[i][j] == Approx(i == j ? 1.0 : 0.0).margin(0.0)); }
    }

    source.Rotate(rotation);
    REQUIRE(source.Px() == Approx(0.0).margin(0.0));
    REQUIRE(source.Py() == Approx(0.0).margin(0.0));
    REQUIRE(source.Pz() == Approx(2.0).margin(0.0));
    REQUIRE(source.E() == Approx(7.0).margin(0.0));
  }

  SECTION("RotationTo keeps a zero spatial source unchanged") {
    M4Vec                 source(0.0, 0.0, 0.0, 7.0);
    const MMatrix<double> rotation = source.RotationTo({1.0, 0.0, 0.0});

    for (std::size_t i = 0; i < 3; ++i) {
      for (std::size_t j = 0; j < 3; ++j) { REQUIRE(rotation[i][j] == Approx(i == j ? 1.0 : 0.0).margin(0.0)); }
    }

    source.RotateTo({1.0, 0.0, 0.0});
    REQUIRE(source.P3mod2() == Approx(0.0).margin(0.0));
    REQUIRE(source.E() == Approx(7.0).margin(0.0));
  }

  SECTION(
      "LagrangeLightCone gives the Euclidean nearest point on either "
      "energy sheet") {
    const auto require_projection = [](const M4Vec &input, const M4Vec &expected) {
      const M4Vec q = gra::kinematics::LagrangeLightCone(input);
      REQUIRE(q.Px() == Approx(expected.Px()).margin(1e-13));
      REQUIRE(q.Py() == Approx(expected.Py()).margin(1e-13));
      REQUIRE(q.Pz() == Approx(expected.Pz()).margin(1e-13));
      REQUIRE(q.E() == Approx(expected.E()).margin(1e-13));
      REQUIRE(std::abs(q.M2()) < 1e-12);
      const double spatial_difference2 = gra::math::pow2(q.Px() - input.Px()) + gra::math::pow2(q.Py() - input.Py()) +
                                         gra::math::pow2(q.Pz() - input.Pz());
      const double distance2          = spatial_difference2 + gra::math::pow2(q.E() - input.E());
      const double expected_distance2 = 0.5 * gra::math::pow2(input.P3mod() - std::abs(input.E()));
      REQUIRE(distance2 == Approx(expected_distance2).margin(1e-12));
    };

    require_projection(M4Vec(3.0, 4.0, 0.0, 10.0), M4Vec(4.5, 6.0, 0.0, 7.5));
    require_projection(M4Vec(3.0, 4.0, 0.0, 2.0), M4Vec(2.1, 2.8, 0.0, 3.5));
    require_projection(M4Vec(3.0, 4.0, 0.0, -10.0), M4Vec(4.5, 6.0, 0.0, -7.5));
    require_projection(M4Vec(0.0, 0.0, 0.0, 2.0), M4Vec(0.0, 0.0, 1.0, 1.0));
    require_projection(M4Vec(0.0, 0.0, 0.0, -2.0), M4Vec(0.0, 0.0, 1.0, -1.0));
  }

  SECTION("CovariantOnShellInitialState preserves final-state four-vectors") {
    const M4Vec  q1(0.30, -0.20, 5.00, 4.95);
    const M4Vec  q2(-0.10, 0.25, -3.00, 3.10);
    const M4Vec  hard = q1 + q2;
    const double mass = hard.M();
    REQUIRE(mass > 0.0);

    std::vector<M4Vec> final = {M4Vec(0.7, -0.4, 0.2, 0.0), M4Vec(-0.5, 0.45, 1.8, 0.0)};
    for (auto &particle : final) { particle.SetE(particle.P3mod()); }
    const M4Vec current    = final[0] + final[1];
    const M4Vec correction = hard - current;
    final[1] += correction;
    REQUIRE((final[0] + final[1] - hard).P3mod() < 1e-12);

    const std::vector<M4Vec> original = final;
    M4Vec                    k1;
    M4Vec                    k2;
    REQUIRE(gra::kinematics::CovariantOnShellInitialState(q1, q2, final, k1, k2));
    REQUIRE(std::abs(k1.M2()) < 1e-10);
    REQUIRE(std::abs(k2.M2()) < 1e-10);
    const M4Vec delta = k1 + k2 - hard;
    REQUIRE(std::abs(delta.E()) < 1e-12);
    REQUIRE(delta.P3mod() < 1e-12);
    for (std::size_t i = 0; i < final.size(); ++i) {
      REQUIRE(gra::math::IsExactEqual(final[i].Px(), original[i].Px()));
      REQUIRE(gra::math::IsExactEqual(final[i].Py(), original[i].Py()));
      REQUIRE(gra::math::IsExactEqual(final[i].Pz(), original[i].Pz()));
      REQUIRE(gra::math::IsExactEqual(final[i].E(), original[i].E()));
    }
  }

  SECTION("CovariantOnShellInitialState preserves an on-shell incoming pair") {
    const M4Vec              q1(0.0, 0.0, 5.0, 5.0);
    const M4Vec              q2(0.0, 0.0, -5.0, 5.0);
    const std::vector<M4Vec> final = {M4Vec(4.0, 0.0, 3.0, 5.0), M4Vec(-4.0, 0.0, -3.0, 5.0)};
    M4Vec                    k1;
    M4Vec                    k2;
    REQUIRE(gra::kinematics::CovariantOnShellInitialState(q1, q2, final, k1, k2));
    REQUIRE((k1 - q1).P3mod() < 1e-12);
    REQUIRE(std::abs(k1.E() - q1.E()) < 1e-12);
    REQUIRE((k2 - q2).P3mod() < 1e-12);
    REQUIRE(std::abs(k2.E() - q2.E()) < 1e-12);
  }

  SECTION(
      "OffShell2LightCone reports an invalid phase-space point without "
      "throwing") {
    M4Vec              p1(0.0, 0.0, 1.0, 0.5);
    M4Vec              p2(0.0, 0.0, -1.0, 0.5);
    std::vector<M4Vec> final    = {M4Vec(0.0, 0.0, 0.0, 2.0), M4Vec(0.0, 0.0, 0.0, 2.0)};
    const M4Vec        input_p1 = p1;
    const M4Vec        input_p2 = p2;

    REQUIRE_FALSE(gra::kinematics::OffShell2LightCone(p1, p2, final));
    REQUIRE(p1.E() == Approx(input_p1.E()));
    REQUIRE(p1.Pz() == Approx(input_p1.Pz()));
    REQUIRE(p2.E() == Approx(input_p2.E()));
    REQUIRE(p2.Pz() == Approx(input_p2.Pz()));
    REQUIRE(final[0].E() == Approx(2.0));
    REQUIRE(final[1].E() == Approx(2.0));
  }

  SECTION("OffShell2LightCone returns a closed on-shell projection") {
    M4Vec              p1(0.0, 0.0, 1.0, 0.5);
    M4Vec              p2(0.0, 0.0, -1.0, 0.5);
    std::vector<M4Vec> final = {M4Vec(0.5, 0.0, 0.0, 0.5), M4Vec(-0.5, 0.0, 0.0, 0.5)};

    REQUIRE(gra::kinematics::OffShell2LightCone(p1, p2, final));
    const M4Vec delta = p1 + p2 - final[0] - final[1];
    REQUIRE(std::abs(p1.M2()) < 1e-12);
    REQUIRE(std::abs(p2.M2()) < 1e-12);
    REQUIRE(std::abs(final[0].M2()) < 1e-12);
    REQUIRE(std::abs(final[1].M2()) < 1e-12);
    REQUIRE(std::abs(delta.Px()) < 1e-12);
    REQUIRE(std::abs(delta.Py()) < 1e-12);
    REQUIRE(std::abs(delta.Pz()) < 1e-12);
    REQUIRE(std::abs(delta.E()) < 1e-12);
  }
}

TEST_CASE("Covariant incoming closure acceptance is Lorentz stable", "[Kinematics][Hardening][CovariantOnShell]") {
  constexpr double         mass = 200.0;
  const M4Vec              q1_base(0.17, -0.11, 100.0, 100.0);
  const M4Vec              q2_base(-0.17, 0.11, -100.0, 100.0);
  const std::vector<M4Vec> final_base = {M4Vec(80.0, 0.0, 60.0, 100.0), M4Vec(-80.0, 0.0, -60.0, 100.0)};

  // Test one rest-frame residual after a common longitudinal boost
  const auto accepts = [&](const double residual, const double rapidity) {
    const M3Vec beta = {0.0, 0.0, std::tanh(rapidity)};
    M4Vec       q1   = q1_base;
    q1.SetE(q1.E() + residual);
    q1                       = q1.LorentzBoost(beta);
    const M4Vec        q2    = q2_base.LorentzBoost(beta);
    std::vector<M4Vec> final = final_base;
    for (auto &particle : final) { particle = particle.LorentzBoost(beta); }
    M4Vec k1;
    M4Vec k2;
    return gra::kinematics::CovariantOnShellInitialState(q1, q2, final, k1, k2);
  };

  const double near_closure = 0.5e-9 * mass;
  const double bad_closure  = 2.0e-9 * mass;
  for (const double rapidity : {0.0, -5.0, 5.0}) {
    CAPTURE(rapidity);
    REQUIRE(accepts(near_closure, rapidity));
    REQUIRE_FALSE(accepts(bad_closure, rapidity));
  }

  // Reject non-finite incoming spatial components before projection
  M4Vec nonfinite = q1_base;
  nonfinite.SetPy(std::numeric_limits<double>::quiet_NaN());
  M4Vec k1;
  M4Vec k2;
  REQUIRE_FALSE(gra::kinematics::CovariantOnShellInitialState(nonfinite, q2_base, final_base, k1, k2));
}

TEST_CASE("frame-axis helpers construct deterministic orthonormal rotations", "[Kinematics][Frames][Coverage]") {
  const gra::M3Vec axis     = {1.2, -0.8, 2.1};
  const gra::M3Vec fallback = {0.0, 0.0, 3.0};
  const gra::M3Vec unit     = gra::kinematics::NormalizeFrameAxis(axis, fallback);
  REQUIRE(gra::SquaredNorm(unit) == Approx(1.0).margin(1e-14));

  const gra::M3Vec fallback_unit = gra::kinematics::NormalizeFrameAxis({0.0, 0.0, 0.0}, fallback);
  REQUIRE(fallback_unit[0] == Approx(0.0));
  REQUIRE(fallback_unit[1] == Approx(0.0));
  REQUIRE(fallback_unit[2] == Approx(1.0));

  const gra::M3Vec perpendicular = gra::kinematics::PerpendicularFrameAxis(unit);
  REQUIRE(gra::SquaredNorm(perpendicular) == Approx(1.0).margin(1e-14));
  REQUIRE(gra::InnerProduct(perpendicular, unit) == Approx(0.0).margin(1e-14));

  const gra::M3Vec      beam1        = {0.7, -0.2, 4.1};
  const gra::M3Vec      beam2        = {-0.4, 0.5, -3.3};
  const MMatrix<double> rotation     = gra::kinematics::FrameRotation(beam1, beam2, axis);
  const gra::M3Vec      rotated_axis = rotation * unit;
  REQUIRE(rotated_axis[0] == Approx(0.0).margin(1e-14));
  REQUIRE(rotated_axis[1] == Approx(0.0).margin(1e-14));
  REQUIRE(rotated_axis[2] == Approx(1.0).margin(1e-14));

  for (std::size_t i = 0; i < 3; ++i) {
    for (std::size_t j = 0; j < 3; ++j) {
      double dot = 0.0;
      for (std::size_t k = 0; k < 3; ++k) { dot += rotation[i][k] * rotation[j][k]; }
      REQUIRE(dot == Approx(i == j ? 1.0 : 0.0).margin(1e-14));
    }
  }

  REQUIRE_THROWS_AS(gra::kinematics::NormalizeFrameAxis({0.0, 0.0, 0.0}, {0.0, 0.0, 0.0}), std::invalid_argument);
}

TEST_CASE("all Lorentz-frame interfaces preserve rest-frame closure and mass shells",
          "[Kinematics][Frames][Coverage]") {
  const double system_mass = 4.0;
  M4Vec        system;
  system.SetPxPyPzM(1.1, -0.7, 2.3, system_mass);

  const double mass1    = 0.6;
  const double mass2    = 0.9;
  const double momentum = gra::kinematics::DecayMomentum(system_mass, mass1, mass2);
  M4Vec        daughter1_rest(0.37 * momentum, -0.48 * momentum, std::sqrt(1.0 - 0.37 * 0.37 - 0.48 * 0.48) * momentum,
                              std::sqrt(momentum * momentum + mass1 * mass1));
  M4Vec        daughter2_rest(-daughter1_rest.Px(), -daughter1_rest.Py(), -daughter1_rest.Pz(),
                              std::sqrt(momentum * momentum + mass2 * mass2));
  std::vector<M4Vec> daughters = {daughter1_rest, daughter2_rest};
  for (auto &daughter : daughters) { gra::kinematics::LorentzBoost(system, system_mass, daughter, 1); }
  REQUIRE((psum(daughters) - system).P3mod() < 1e-12);
  REQUIRE(std::abs(psum(daughters).E() - system.E()) < 1e-12);

  const double beam_mass = 0.938;
  M4Vec        beam1;
  M4Vec        beam2;
  beam1.SetPxPyPzM(0.0, 0.0, 7.0, beam_mass);
  beam2.SetPxPyPzM(0.0, 0.0, -7.0, beam_mass);

  M4Vec              beam1_boost;
  M4Vec              beam2_boost;
  std::vector<M4Vec> daughters_boost;
  gra::kinematics::LorentFramePrepare(daughters, system, beam1, beam2, beam1_boost, beam2_boost, daughters_boost);
  REQUIRE(psum(daughters_boost).P3mod() < 1e-12);
  REQUIRE(psum(daughters_boost).E() == Approx(system_mass).margin(1e-12));

  for (const std::string frame : {"CM", "CS", "AH", "HX"}) {
    std::vector<M4Vec> transformed;
    gra::kinematics::LorentzFrame(transformed, beam1_boost, beam2_boost, daughters_boost, frame, -1);
    REQUIRE(psum(transformed).P3mod() < 1e-11);
    REQUIRE(psum(transformed).E() == Approx(system_mass).margin(1e-12));
    REQUIRE(transformed[0].M2() == Approx(mass1 * mass1).margin(1e-12));
    REQUIRE(transformed[1].M2() == Approx(mass2 * mass2).margin(1e-12));
  }

  for (const int direction : {-1, 1}) {
    std::vector<M4Vec> transformed;
    gra::kinematics::LorentzFrame(transformed, beam1_boost, beam2_boost, daughters_boost, "PG", direction);
    REQUIRE(psum(transformed).P3mod() < 1e-11);
    REQUIRE(psum(transformed).E() == Approx(system_mass).margin(1e-12));
  }

  std::vector<M4Vec> generic_ah;
  gra::kinematics::LorentzFrame(generic_ah, beam1_boost, beam2_boost, daughters_boost, "AH", -1);
  std::vector<M4Vec> dedicated_ah = daughters;
  gra::kinematics::AHframe(dedicated_ah, system, beam1, beam2);

  std::vector<M4Vec> generic_cs;
  gra::kinematics::LorentzFrame(generic_cs, beam1_boost, beam2_boost, daughters_boost, "CS", -1);
  std::vector<M4Vec> dedicated_cs = daughters;
  gra::kinematics::CSframe(dedicated_cs, system, beam1, beam2);

  for (std::size_t i = 0; i < daughters.size(); ++i) {
    REQUIRE(dedicated_ah[i].Px() == Approx(generic_ah[i].Px()).margin(1e-12));
    REQUIRE(dedicated_ah[i].Py() == Approx(generic_ah[i].Py()).margin(1e-12));
    REQUIRE(dedicated_ah[i].Pz() == Approx(generic_ah[i].Pz()).margin(1e-12));
    REQUIRE(dedicated_cs[i].Px() == Approx(generic_cs[i].Px()).margin(1e-12));
    REQUIRE(dedicated_cs[i].Py() == Approx(generic_cs[i].Py()).margin(1e-12));
    REQUIRE(dedicated_cs[i].Pz() == Approx(generic_cs[i].Pz()).margin(1e-12));
  }

  std::vector<M4Vec> gj = {beam1};
  gra::kinematics::GJframe(gj, system, -1, beam1, beam2);
  REQUIRE(std::abs(gj[0].Px()) < 1e-12);
  REQUIRE(std::abs(gj[0].Py()) < 1e-12);

  std::vector<M4Vec> pg = {beam2};
  gra::kinematics::PGframe(pg, system, 1, beam1, beam2);
  REQUIRE(std::abs(pg[0].Px()) < 1e-12);
  REQUIRE(std::abs(pg[0].Py()) < 1e-12);

  std::vector<M4Vec> hx = daughters;
  gra::kinematics::HXframe(hx, system);
  REQUIRE(psum(hx).P3mod() < 1e-11);
  REQUIRE(psum(hx).E() == Approx(system_mass).margin(1e-12));

  std::vector<M4Vec> cm = daughters;
  gra::kinematics::CMframe(cm, system);
  REQUIRE(psum(cm).P3mod() < 1e-11);
  REQUIRE(psum(cm).E() == Approx(system_mass).margin(1e-12));

  std::vector<M4Vec> output;
  REQUIRE_THROWS_AS(gra::kinematics::LorentzFrame(output, beam1_boost, beam2_boost, daughters_boost, "PG", 0),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::kinematics::LorentzFrame(output, beam1_boost, beam2_boost, daughters_boost, "UNKNOWN", -1),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::kinematics::GJframe(output, system, 0, beam1, beam2), std::invalid_argument);
  REQUIRE_THROWS_AS(gra::kinematics::PGframe(output, system, 0, beam1, beam2), std::invalid_argument);
}

TEST_CASE("Legacy Lorentz frames reject nonphysical rest systems", "[Kinematics][LorentzFrame][failure]") {
  const M4Vec              spacelike_system(1.0, 0.0, 0.0, 0.5);
  const M4Vec              beam1(0.0, 0.0, 5.0, 5.1);
  const M4Vec              beam2(0.0, 0.0, -5.0, 5.1);
  const std::vector<M4Vec> particles = {M4Vec(0.7, 0.0, 0.0, 0.4)};

  SECTION("common frame preparation") {
    M4Vec              beam1_boost;
    M4Vec              beam2_boost;
    std::vector<M4Vec> boosted;
    REQUIRE_THROWS(gra::kinematics::LorentFramePrepare(particles, spacelike_system, beam1, beam2, beam1_boost,
                                                       beam2_boost, boosted));
  }

  SECTION("dedicated beam frames") {
    auto output = particles;
    REQUIRE_THROWS(gra::kinematics::AHframe(output, spacelike_system, beam1, beam2));
    output = particles;
    REQUIRE_THROWS(gra::kinematics::CSframe(output, spacelike_system, beam1, beam2));
    output = particles;
    REQUIRE_THROWS(gra::kinematics::GJframe(output, spacelike_system, -1, beam1, beam2));
    output = particles;
    REQUIRE_THROWS(gra::kinematics::PGframe(output, spacelike_system, 1, beam1, beam2));
  }

  SECTION("central rest and helicity frames") {
    auto output = particles;
    REQUIRE_THROWS(gra::kinematics::CMframe(output, spacelike_system));
    output = particles;
    REQUIRE_THROWS(gra::kinematics::HXframe(output, beam1));
  }
}

TEST_CASE("fiducial cut activity helpers detect each configured observable family",
          "[Kinematics][Fiducial][Coverage]") {
  FIDPDGCUT pdg_cut;
  REQUIRE_FALSE(pdg_cut.HasActiveRange());
  pdg_cut.Eta.active = true;
  REQUIRE(pdg_cut.HasActiveRange());

  FIDCUT cuts;
  REQUIRE_FALSE(cuts.HasCentralParticleCuts());
  REQUIRE_FALSE(cuts.HasCentralSystemCuts());
  REQUIRE_FALSE(cuts.HasForwardCuts());

  cuts.particle_pt_active = true;
  cuts.system_Rap_active  = true;
  cuts.forward_xi_active  = true;
  REQUIRE(cuts.HasCentralParticleCuts());
  REQUIRE(cuts.HasCentralSystemCuts());
  REQUIRE(cuts.HasForwardCuts());
}

// Factorials
//
//
TEST_CASE("gra::math::factorial: with values of 0,1,2,3,10", "[gra::math::factorial]") {
  REQUIRE(gra::math::IsExactEqual(gra::math::factorial(0), 1.0));
  REQUIRE(gra::math::IsExactEqual(gra::math::factorial(1), 1.0));
  REQUIRE(gra::math::IsExactEqual(gra::math::factorial(2), 2.0));
  REQUIRE(gra::math::IsExactEqual(gra::math::factorial(3), 6.0));
  REQUIRE(gra::math::IsExactEqual(gra::math::factorial(10), 3628800.0));
  REQUIRE(std::isfinite(gra::math::factorial(170)));
  REQUIRE(std::isinf(gra::math::factorial(171)));
}

TEST_CASE("special functions satisfy analytic series identities", "[MMath][SpecialFunctions]") {
  REQUIRE(gra::math::ReciprocalGamma(5.0) == Approx(1.0 / 24.0));
  REQUIRE(gra::math::ReciprocalGamma(0.0) == Approx(0.0));
  REQUIRE(gra::math::Hyper2F1(1.0, 2.0, 2.0, 0.25) == Approx(4.0 / 3.0).epsilon(2.0e-14));
  REQUIRE(gra::math::Hyper2F1(1.0, 1.0, 3.0, 1.0) == Approx(2.0).epsilon(2.0e-14));

  const double          magnitude = 0.7;
  constexpr std::size_t terms     = 8;
  double                partial   = 0.0;
  for (std::size_t n = 0; n < terms; ++n) {
    partial += std::pow(magnitude, static_cast<double>(n)) / gra::math::factorial(static_cast<int>(n));
  }
  REQUIRE(std::exp(magnitude) - partial <= gra::math::ExponentialTailBound(magnitude, terms));
  REQUIRE_THROWS_AS(gra::math::ExponentialTailBound(-1.0, terms), std::invalid_argument);

  // Check real to complex to real identity for every positive and negative m
  const std::array<double, 5> costheta = {-0.83, -0.37, 0.0, 0.41, 0.89};
  const std::array<double, 5> phi      = {-2.71, -0.81, 0.0, 1.13, 2.48};
  for (const double costheta_value : costheta) {
    for (const double phi_value : phi) {
      for (int l = 0; l <= 4; ++l) {
        for (int m = -l; m <= l; ++m) {
          CAPTURE(costheta_value, phi_value, l, m);
          const auto   complex  = gra::math::Y_complex_basis(costheta_value, phi_value, l, m);
          const auto   analytic = gra::math::Y_complex_basis_ref(costheta_value, phi_value, l, m);
          const double real     = gra::math::Y_real_basis(costheta_value, phi_value, l, m);
          REQUIRE(gra::math::NReY(complex, m) == Approx(real).margin(2.0e-14));
          REQUIRE(gra::math::NReY(analytic, m) == Approx(real).margin(2.0e-14));
        }
      }
    }
  }
}

// Legendre polynomials
//
//
TEST_CASE("gra::math::sf_legendre: numerically l = 0 ... 8, m = 0", "[gra::math::sf_legendre]") {
  const double EPS = 1e-5;

  std::vector<double> costheta = {-0.7, -0.3, 0.0, 0.3, 0.7};

  for (std::size_t i = 0; i < costheta.size(); ++i) {
    for (int l = 0; l < 8; ++l) {
      REQUIRE(gra::math::sf_legendre(l, 0, costheta[i]) == Approx(gra::math::LegendrePl(l, costheta[i])).epsilon(EPS));
    }
  }
}

// Complex spherical harmonics
//
//
TEST_CASE("gra::math::Y_complex_basis: numerically", "[gra::math::Y_complex_basis]") {
  const double EPS = 1e-5;

  std::vector<double> costheta = {-0.7, -0.3, 0.0, 0.3, 0.7};
  std::vector<double> phi      = {-2.5, 2.0, 1.0, 1.5, 2.5};

  for (std::size_t i = 0; i < costheta.size(); ++i) {
    for (std::size_t j = 0; j < phi.size(); ++j) {
      for (int l = 0; l <= 4; ++l) {
        for (int m = -l; m <= l; ++m) {
          SECTION("l = " + std::to_string(l) + " m = " + std::to_string(m)) {
            REQUIRE(std::real(gra::math::Y_complex_basis(costheta[i], phi[j], l, m)) ==
                    Approx(std::real(gra::math::Y_complex_basis_ref(costheta[i], phi[j], l, m))).epsilon(EPS));
            REQUIRE(std::imag(gra::math::Y_complex_basis(costheta[i], phi[j], l, m)) ==
                    Approx(std::imag(gra::math::Y_complex_basis_ref(costheta[i], phi[j], l, m))).epsilon(EPS));
          }
        }
      }
    }
  }
}

// Check every real spherical harmonic is normalized on the unit sphere
TEST_CASE("gra::math real spherical harmonics have unit norm",
          "[gra::math::Y_real_basis]") {
  constexpr std::size_t n_costheta = 640;
  constexpr std::size_t n_phi = 384;
  constexpr double volume = 4.0 * gra::math::PI;
  for (int l = 0; l <= 4; ++l) {
    for (int m = -l; m <= l; ++m) {
      double integral = 0.0;
      for (std::size_t i = 0; i < n_costheta; ++i) {
        const double costheta =
            -1.0 + 2.0 * (static_cast<double>(i) + 0.5) / n_costheta;
        for (std::size_t j = 0; j < n_phi; ++j) {
          const double phi = 2.0 * gra::math::PI *
                             (static_cast<double>(j) + 0.5) / n_phi;
          const double value =
              gra::math::Y_real_basis(costheta, phi, l, m);
          integral += value * value;
        }
      }
      CAPTURE(l, m);
      CHECK(volume * integral /
                static_cast<double>(n_costheta * n_phi) ==
            Approx(1.0).margin(1.0e-4));
    }
  }
}

// Reject boost failures before publishing positive decay phase-space weights
TEST_CASE("decay phase space rejects inconsistent mothers", "[Kinematics][PhaseSpace]") {
  gra::MRandom random;
  random.SetSeed(173);
  const auto mother = GENERATE(gra::M4Vec(0, 0, 0, 3), gra::M4Vec(0, 0, 0, -2),
                               gra::M4Vec(0, 0, 3, 2));
  std::vector<gra::M4Vec> products(4);
  REQUIRE(gra::kinematics::TwoBodyPhaseSpace(mother, 2.0, {0.1, 0.1}, products, random).GetW() < 0.0);
  REQUIRE(gra::kinematics::ThreeBodyPhaseSpace(mother, 2.0, {0.1, 0.1, 0.1}, products, false, random).GetW() < 0.0);
  REQUIRE(gra::kinematics::NBodyPhaseSpace(mother, 2.0, {0.1, 0.1, 0.1, 0.1}, products, false, random).GetW() < 0.0);
  REQUIRE(gra::kinematics::RamboMassless(mother, 2.0, products, random).GetW() < 0.0);
  REQUIRE(gra::kinematics::RamboMassive(mother, 2.0, {0.1, 0.1, 0.1, 0.1}, products, random).GetW() < 0.0);
}

// Preserve soft remnant directions, beam exchange and positive production energy
TEST_CASE("soft remnant recoil is physical", "[Kinematics][SoftChain]") {
  const double p = std::sqrt(10000.0 - gra::math::pow2(gra::PDG::mp));
  std::array<gra::M4Vec, 2> beams = {gra::M4Vec(0, 0, p, 100), gra::M4Vec(0, 0, -p, 100)};
  const std::array<gra::M4Vec, 2> transverse = {gra::M4Vec(0.4, 0, 0, 0), gra::M4Vec(-0.4, 0, 0, 0)};
  std::array<gra::M4Vec, 2> remnants;
  gra::M4Vec system;
  REQUIRE_FALSE(gra::kinematics::BuildSoftRemnants(beams, {0.001, 0.001}, {900.0, 900.0}, transverse, remnants, system));
  REQUIRE(gra::kinematics::BuildSoftRemnants(beams, {0.3, 0.2}, {100.0, 144.0}, transverse, remnants, system));
  REQUIRE(remnants[0].Pz() > 0.0);
  REQUIRE(remnants[1].Pz() < 0.0);
  REQUIRE(system.E() > 0.0);
  REQUIRE(remnants[0].M2() == Approx(100.0).margin(1e-9));
  REQUIRE(remnants[1].M2() == Approx(144.0).margin(1e-9));
  REQUIRE(gra::math::CheckEMC(beams[0] + beams[1] - remnants[0] - remnants[1] - system));
  const auto nominal = remnants;
  const auto nominal_system = system;
  std::array<gra::M4Vec, 2> exchanged;
  gra::M4Vec exchanged_system;
  REQUIRE(gra::kinematics::BuildSoftRemnants({beams[1], beams[0]}, {0.2, 0.3}, {144.0, 100.0},
                                             {transverse[1], transverse[0]}, exchanged, exchanged_system));
  REQUIRE(gra::math::CheckEMC(exchanged[0] - nominal[1]));
  REQUIRE(gra::math::CheckEMC(exchanged[1] - nominal[0]));
  REQUIRE(gra::math::CheckEMC(exchanged_system - nominal_system));
  for (auto &beam : beams) { beam.RotateY(0.7); }
  REQUIRE(gra::kinematics::BuildSoftRemnants(beams, {0.3, 0.2}, {100.0, 144.0}, transverse, remnants, system));
  for (const auto &i : gra::aux::indices(remnants)) {
    auto rotated = nominal[i];
    rotated.RotateY(0.7);
    REQUIRE(gra::math::CheckEMC(rotated - remnants[i]));
  }
  auto rotated_system = nominal_system;
  rotated_system.RotateY(0.7);
  REQUIRE(gra::math::CheckEMC(rotated_system - system));
}

// Resolve diffractive t bounds close to forward scattering and reject subthreshold states
TEST_CASE("two-body t limits preserve small forward transfers", "[Kinematics][Scattering]") {
  const double s = 13000.0 * 13000.0;
  const double m2 = gra::math::pow2(gra::PDG::mp);
  double tmin = 0.0;
  double tmax = 0.0;
  gra::kinematics::Two2TwoLimit(s, m2, m2, m2, m2, tmin, tmax);
  REQUIRE(tmax == Approx(0.0).margin(1e-30));
  REQUIRE(tmin == Approx(4 * m2 - s).epsilon(1e-14));
  gra::kinematics::Two2TwoLimit(s, m2, m2, 4.0, 9.0, tmin, tmax);
  REQUIRE(tmax < 0.0);
  REQUIRE(std::abs(tmax) < 1e-5);
  double pz = 0.0;
  double pt = 0.0;
  REQUIRE(gra::kinematics::ForwardScattering(s, m2, m2, 4.0, 9.0, 2 * tmax, pz, pt));
  REQUIRE_FALSE(gra::kinematics::ForwardScattering(s, m2, m2, 4.0, 9.0, 0.5 * tmax, pz, pt));
  gra::kinematics::Two2TwoLimit(4.0, m2, m2, 4.0, 9.0, tmin, tmax);
  REQUIRE(tmin > tmax);
  REQUIRE_FALSE(gra::kinematics::ForwardScattering(4.0, m2, m2, 4.0, 9.0, -0.1, pz, pt));
  REQUIRE_FALSE(gra::kinematics::ForwardScattering(s, -1.0, m2, m2, m2, -0.1, pz, pt));
}

// Retain finite angular variables and norms across changes of momentum units
TEST_CASE("M4Vec observables retain finite forward directions", "[M4Vec][Kinematics][Numerics]") {
  for (const double scale : {1e-200, 1.0, 1e200}) {
    const M4Vec p(3.0 * scale, 4.0 * scale, 12.0 * scale, 13.0 * scale);
    REQUIRE(p.Pt() / scale == Approx(5.0));
    REQUIRE(p.P3mod() / scale == Approx(13.0));
    REQUIRE(p.Eta() == Approx(std::asinh(12.0 / 5.0)));
    REQUIRE(p.Rap() == Approx(0.5 * std::log(25.0)));
  }
  REQUIRE(M4Vec(1.0, 0.0, 1e9, 1e9).Eta() == Approx(std::asinh(1e9)));
  REQUIRE(M4Vec(1e9, 0.0, 0.0, 1.0).Mt2() == Approx(1.0));
  for (const double sign : {-1.0, 1.0}) {
    REQUIRE(M4Vec(0.0, 0.0, sign * 1e308, sign * 1.5e308).Rap() == Approx(0.5 * std::log(5.0)));
    REQUIRE(M4Vec(0.0, 0.0, sign * 1e-20, 1.0).Rap() / (sign * 1e-20) == Approx(1.0));
  }
  REQUIRE(std::isnan(M4Vec().Eta()));
  REQUIRE(std::isnan(M4Vec().Rap()));
  REQUIRE(std::isnan(M4Vec(0.0, 0.0, 2.0, 1.0).Rap()));
}

// A rotation must align directions even near zero or pi and at any momentum scale
TEST_CASE("RotationTo resolves small angles independently of momentum scale", "[M4Vec][Kinematics][Numerics]") {
  for (const double scale : {1e-200, 1.0, 1e200}) {
    for (const double angle : {1e-7, 0.7, gra::math::PI - 1e-7}) {
      M4Vec source(0.0, 0.0, scale, 2.0 * scale);
      const M3Vec target = {scale * std::sin(angle), 0.0, scale * std::cos(angle)};
      const auto rotation = source.RotationTo(target);
      source.Rotate(rotation);
      REQUIRE(source.Px() / scale == Approx(std::sin(angle)).epsilon(1e-10));
      REQUIRE(source.Py() / scale == Approx(0.0).margin(1e-14));
      REQUIRE(source.Pz() / scale == Approx(std::cos(angle)).epsilon(1e-14));
      REQUIRE(source.P3mod() / scale == Approx(1.0).epsilon(1e-14));
      REQUIRE(source.E() / scale == Approx(2.0));
    }
  }
}

// Check alignment and orthogonality for general directions close to pi
TEST_CASE("RotationTo remains orthogonal near opposite general directions", "[M4Vec][Kinematics][Numerics]") {
  std::mt19937 engine(421);
  std::uniform_real_distribution<double> uniform(-1.0, 1.0);
  for (int trial = 0; trial < 100; ++trial) {
    const M3Vec direction = gra::NormalizedL2(M3Vec{uniform(engine), uniform(engine), uniform(engine)});
    const M3Vec reference = {uniform(engine), uniform(engine), uniform(engine)};
    const M3Vec axis = gra::NormalizedL2(gra::CrossProduct(direction, reference));
    for (const double angle : {0.0, 1e-16, 1e-10, 0.7, gra::math::PI - 1e-7, gra::math::PI - 1e-10, gra::math::PI}) {
      M3Vec target = direction;
      gra::Scale(target, std::cos(angle));
      gra::AddScaled(target, gra::CrossProduct(axis, direction), std::sin(angle));
      M4Vec p(direction[0], direction[1], direction[2], 2.0);
      const auto rotation = p.RotationTo(target);
      p.Rotate(rotation);
      const auto actual = p.P3();
      for (const auto &i : indices(target)) {
        REQUIRE(actual[i] == Approx(target[i]).epsilon(0.0).margin(5e-15));
      }
      const auto product = rotation * rotation.Transpose();
      for (std::size_t i = 0; i < 3; ++i) {
        for (std::size_t j = 0; j < 3; ++j) {
          REQUIRE(product(i, j) == Approx(i == j ? 1.0 : 0.0).epsilon(0.0).margin(5e-15));
        }
      }
    }
  }
}

// One massless daughter gives sqrt(lambda(x,0,z)) = abs(x-z) exactly
TEST_CASE("Kallen triangle retains phase space next to threshold", "[Kinematics][Numerics]") {
  for (const double scale : {1e-100, 1.0, 1e100}) {
    const double below = std::nextafter(scale, 0.0);
    std::array<double, 3> q = {0.0, below, scale};
    do {
      REQUIRE(SqrtKallenLambda(q[0], q[1], q[2]) / (scale - below) == Approx(1.0).epsilon(1e-14));
    } while (std::next_permutation(q.begin(), q.end()));
  }
  REQUIRE(SqrtKallenLambda(25.0, 1.0, 4.0) == Approx(std::sqrt(384.0)));
  REQUIRE(SqrtKallenLambda(1.0, 1.0, 1.0) == Approx(0.0));
  REQUIRE(SqrtKallenLambda(0.0, -1.0, 1.0) == Approx(2.0));
}

// Reject a boost whose extended precision result exceeds four-vector storage
TEST_CASE("LorentzBoost rejects overflow on conversion to double", "[LorentzBoost][Kinematics][Numerics]") {
  M4Vec p(0.0, 0.0, 0.0, std::numeric_limits<double>::max());
  LorentzBoost(M4Vec(0.0, 0.0, 3.0, 5.0), 4.0, p, 1);
  REQUIRE(p.E() < 0.0);
  REQUIRE(std::isfinite(p.Pz()));
}

// Compare the small phase-space volume with direct integration, including the first allowed double
TEST_CASE("dPhi1 remains positive and accurate near threshold", "[Kinematics][Numerics]") {
  for (const double energy : {std::nextafter(1.0, 2.0), 1.0 + 1e-10, 1.00001, 1.0001}) {
    const double reference = NumericalDPhi1(energy, 1.0);
    REQUIRE(reference > 0.0);
    REQUIRE(dPhi1(energy, 1.0) / reference == Approx(1.0).epsilon(1e-10));
  }
}

// For equal opposing momenta the invariant incoming flux is 8 E |p|
TEST_CASE("Moller flux retains the nonrelativistic limit", "[Kinematics][Numerics]") {
  for (const double momentum : {1e-12, 1e-8, 0.1, 10.0, 6500.0}) {
    const double energy = std::hypot(1.0, momentum);
    const M4Vec a(0.0, 0.0, momentum, energy), b(0.0, 0.0, -momentum, energy);
    const double reference = 8.0 * energy * momentum;
    REQUIRE(MollerFlux(a, b) / reference == Approx(1.0).epsilon(1e-14));
    REQUIRE(MollerFlux(b, a) / reference == Approx(1.0).epsilon(1e-14));
    REQUIRE(MollerFlux(a, a) == Approx(0.0).margin(1e-24));
  }
  const M4Vec a(1.0, 2.0, 3.0, 5.0), b(-2.0, 0.5, -1.0, 4.0);
  const double reference = 4.0 * std::sqrt(pow2(a.DotM(b)) - a.M2() * b.M2());
  REQUIRE(MollerFlux(a, b) == Approx(reference).epsilon(1e-14));
  REQUIRE(MollerFlux(a.LorentzBoost({0.2, -0.1, 0.3}), b.LorentzBoost({0.2, -0.1, 0.3})) ==
          Approx(reference).epsilon(1e-14));
}

// Construct finite on-shell energies without squaring beyond double range
TEST_CASE("M4Vec mass setters retain scale and invalid masses remain invalid", "[M4Vec][Kinematics][Numerics]") {
  for (const double scale : {1e-200, 1.0, 1e200}) {
    M4Vec p;
    p.SetPxPyPzM(3.0 * scale, 4.0 * scale, 0.0, 12.0 * scale);
    REQUIRE(p.E() / scale == Approx(13.0));
  }
  const M4Vec invalid(0.0, 0.0, 0.0, std::numeric_limits<double>::quiet_NaN());
  REQUIRE(std::isnan(invalid.M()));
  REQUIRE(std::isnan(invalid.Mt()));
}

// Check zero-truncated Poisson sampling against its normalized probability law
TEST_CASE("Positive Poisson sampling retains the conditional law", "[random][poisson]") {
  gra::MRandom random;
  random.SetSeed(41821);
  CHECK_THROWS_AS(random.PoissonRandom(0.0, true), std::invalid_argument);
  constexpr unsigned int trials = 100000;
  for (const double mean : {1.0e-12, 0.05, 0.8, 3.0, 25.0}) {
    double sum = 0.0, square = 0.0;
    unsigned int single = 0;
    for (unsigned int i = 0; i < trials; ++i) {
      const int count = random.PoissonRandom(mean, true);
      REQUIRE(count > 0);
      sum += count;
      square += static_cast<double>(count) * count;
      single += count == 1;
    }
    const double probability = -std::expm1(-mean), expected = mean / probability;
    const double variance = mean * (1.0 + mean) / probability - expected * expected;
    const double one = mean * std::exp(-mean) / probability;
    CAPTURE(mean, sum / trials, square / trials);
    CHECK(sum / trials == Approx(expected).margin(6.0 * std::sqrt(std::max(0.0, variance) / trials) + 1e-10));
    CHECK(static_cast<double>(single) / trials == Approx(one).margin(6.0 * std::sqrt(one * (1.0 - one) / trials) + 1e-10));
  }
}
