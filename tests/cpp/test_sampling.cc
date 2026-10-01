// Agreement of F and C event weights and cascaded Jacob-Wick sampling
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>

#include <array>
#include <cmath>
#include <complex>
#include <cstdio>
#include <vector>

#include "Graniitti/MGraniitti.h"
#include "Graniitti/Sampling/MVEGAS.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Spin/MSpin.h"

namespace {

using gra::aux::indices;
using Moments = std::array<gra::kinematics::MCW, 7>;

// Accumulate absolute rate, mass, rapidity and azimuth moments including rejected trials
void AddSamplingMoments(Moments &moments, double weight, const gra::LORENTZSCALAR &lts) {
  std::array<double, 7> values{};
  if (weight > 0.0) {
    const double phi = lts.pfinal[1].Phi() - lts.pfinal[2].Phi();
    const auto &pair = lts.decaytree.front();
    const double pair_mass2 = pair.legs.empty() ? (pair.p4 + lts.decaytree[1].p4).M2() : pair.p4.M2();
    values = {1.0, lts.pfinal[0].M2() / 10.0, gra::math::pow2(lts.pfinal[0].Rap()),
              gra::math::pow2(std::cos(phi)), pair_mass2 / 10.0, lts.pfinal[0].Rap() > 0.0 ? 1.0 : 0.0,
              lts.pfinal[0].Rap() <= 0.0 ? 1.0 : 0.0};
  }
  for (const auto &i : indices(moments)) { moments[i].Push(weight * values[i]); }
}

// Adapt VEGAS, freeze its proposal and integrate the real EventWeight API with independent random draws
Moments IntegrateSampling(const std::string &channel, const std::string &decay, const std::string &phase,
                          bool flat, unsigned int seed, unsigned int calls,
                          double energy1 = 80.0, double energy2 = 20.0,
                          double mass_min = 1.3, double mass_max = 3.2) {
  CAPTURE(channel, decay, phase, flat);
  auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
  card["GENERIC"]["CORES"] = 1;
  card["GENERIC"]["HIST"] = 0;
  card["GENERIC"]["RNDSEED"] = seed;
  card["GENERIC"]["INTEGRATOR"] = "VEGAS";
  card["SCATTERING"]["PROCESS"] = channel + "<" + phase + "> " + decay;
  card["SCATTERING"]["ENERGY"] = {energy1, energy2};
  card["SCATTERING"]["LOOPSCREEN"] = false;
  for (const std::string mode : {"<F>", "<C>"}) {
    card["GENCUTS"][mode]["M"] = {mass_min, mass_max};
    // Positive-energy daughter sums stay inside their common fiducial rapidity interval
    card["GENCUTS"][mode]["Rap"] = {-1.0, 1.0};
    card["GENCUTS"][mode]["Pt"] = {0.1, 0.6};
  }
  card["FIDCUTS"] = {{"active", true}, {"CENTRAL", {{"*", {{"Rap", {-1.0, 1.0}}}},
                     {"SYSTEM", {{"M", {mass_min, mass_max}}}}}}, {"USERCUTS", false}};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto &process = *generator.proc;
  if (phase == "F") {
    REQUIRE(dynamic_cast<gra::MFactorized *>(&process) != nullptr);
  } else {
    REQUIRE(dynamic_cast<gra::MCentral *>(&process) != nullptr);
  }
  process.SetFLATMASS2(flat);
  process.SetOFFSHELL(1.5);
  process.PrepareRun();
  gra::MRandom random;
  random.SetSeed(seed + 173);
  gra::VEGASPARAM param;
  param.BINS = 24;
  param.AUTOMATIC_CONVERGENCE = false;
  gra::MVEGASIntegrator vegas;
  vegas.SetDimension(process.GetdLIPSDim());
  vegas.Init(gra::VEGASStage::Adaptation, 4000, 6, param);
  const auto evaluate = [&](const gra::VEGASSample &sample, bool adapting) {
    gra::MEventWeightState aux;
    aux.include_screening = false;
    aux.adaptation_mode = adapting;
    aux.log_inverse_density = sample.log_inverse_density;
    const double weight = process.EventWeight(sample.point, aux) * std::exp(sample.log_inverse_density);
    if (aux.technical_failure || !std::isfinite(weight) || weight < 0.0) {
      CAPTURE(aux.technical_failure, weight);
      FAIL("Invalid sampling weight");
    }
    return weight;
  };
  while (!vegas.IsFrozen()) {
    const auto calls = vegas.NextAdaptationCalls();
    vegas.BeginAdaptationBatch();
    for (std::uint64_t i = 0; i < calls; ++i) {
      const auto sample = vegas.Sample([&]() { return random.U(0.0, 1.0); });
      vegas.AccumulateAdaptation(evaluate(sample, true), sample.indices);
    }
    const auto report = vegas.FinishAdaptationBatch(calls);
    REQUIRE(report.status != gra::VEGASAdaptationStatus::InsufficientSupport);
  }
  vegas.Init(gra::VEGASStage::Integration, calls, 1, param);
  Moments moments;
  for (unsigned int i = 0; i < calls; ++i) {
    const auto sample = vegas.Sample([&]() { return random.U(0.0, 1.0); });
    const double weight = evaluate(sample, false);
    AddSamplingMoments(moments, weight, process.state.lts);
  }
  for (const auto &i : indices(moments)) {
    std::printf("VEGAS %s <%s> %s moment %zu: %.4e +- %.4e\n", channel.c_str(), phase.c_str(),
                flat ? "uniform m2" : "Breit-Wigner", i, moments[i].Integral(), moments[i].IntegralError());
  }
  REQUIRE(moments.front().Integral() > 0.0);
  REQUIRE(moments.front().IntegralError() / moments.front().Integral() < 0.02);
  return moments;
}

// Compare independent integrals with a fixed precision requirement and combined statistical errors
void RequireSamplingAgreement(const Moments &reference, const Moments &candidate) {
  for (const auto &i : indices(reference)) {
    const double first = reference[i].Integral();
    const double second = candidate[i].Integral();
    const double error = std::hypot(reference[i].IntegralError(), candidate[i].IntegralError());
    CAPTURE(i, first, second, error);
    REQUIRE(first > 0.0);
    REQUIRE(second > 0.0);
    REQUIRE(error / first < 0.04);
    REQUIRE(std::abs(first - second) < 5.0 * error);
  }
}

}  // namespace

// Integrate the same physical final state and lab cuts with independent F and C maps
TEST_CASE("F and C agree on direct phase space rates and moments", "[Process][Sampling][VEGAS][FC]") {
  const auto multiplicity = GENERATE(2U, 4U);
  // Resolve the broader four pion weight distribution to the same precision
  const unsigned int calls = multiplicity == 2 ? 16000 : 160000;
  CAPTURE(multiplicity);
  const std::string channel = multiplicity == 2 ? "yy[EPA]" : "MP[CON]";
  const std::string decay = multiplicity == 2 ? "-> mu+ mu-" : "-> pi+ pi- pi+ pi-";
  const auto factorized = IntegrateSampling(channel, decay, "F", false, 27183, calls);
  const auto central = IntegrateSampling(channel, decay, "C", false, 81937, calls);
  RequireSamplingAgreement(factorized, central);
  if (multiplicity == 2) {
    const auto exchanged = IntegrateSampling(channel, decay, "C", false, 47149, calls, 20.0, 80.0);
    for (std::size_t i = 0; i < 5; ++i) {
      const double error = std::hypot(central[i].IntegralError(), exchanged[i].IntegralError());
      CHECK(std::abs(central[i].Integral() - exchanged[i].Integral()) < 5.0 * error);
    }
    for (const std::size_t i : {5U, 6U}) {
      const std::size_t opposite = 11U - i;
      const double error = std::hypot(central[i].IntegralError(), exchanged[opposite].IntegralError());
      CHECK(std::abs(central[i].Integral() - exchanged[opposite].Integral()) < 5.0 * error);
    }
  }
}

// Compare the real rho pair amplitudes and their pion decays with both mass proposals
TEST_CASE("F and C preserve rho cascade rates with Breit-Wigner and uniform mass squared proposals",
          "[Process][Sampling][VEGAS][FC][Cascade]") {
  constexpr unsigned int calls = 160000;
  const std::string decay = "-> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}";
  const auto factorized_bw = IntegrateSampling("MP[CON]", decay, "F", false, 73471, calls);
  const auto central_bw = IntegrateSampling("MP[CON]", decay, "C", false, 91573, calls);
  // Uniform mass squared proposals need more samples to resolve the Breit-Wigner peaks
  const auto factorized_flat = IntegrateSampling("MP[CON]", decay, "F", true, 51361, 2 * calls);
  const auto central_flat = IntegrateSampling("MP[CON]", decay, "C", true, 31991, 2 * calls);
  RequireSamplingAgreement(factorized_bw, central_bw);
  RequireSamplingAgreement(factorized_bw, factorized_flat);
  RequireSamplingAgreement(central_bw, central_flat);
  RequireSamplingAgreement(factorized_flat, central_flat);
}

// Integrate a generated amplitude with its own W propagators and no Jacob-Wick decay factors
TEST_CASE("MG5 WW cascades agree between F and C and both mass proposals", "[Process][Sampling][VEGAS][FC][Cascade][MG5]") {
  constexpr unsigned int calls = 80000;
  const std::string decay = "-> W+ > {e+ ve} W- > {mu- vm~}";
  const auto factorized_bw = IntegrateSampling("yy[WW]", decay, "F", false, 91373, calls, 650.0, 650.0, 155.0, 300.0);
  const auto central_bw = IntegrateSampling("yy[WW]", decay, "C", false, 73129, calls, 650.0, 650.0, 155.0, 300.0);
  const auto factorized_flat = IntegrateSampling("yy[WW]", decay, "F", true, 51913, 2 * calls, 650.0, 650.0, 155.0, 300.0);
  const auto central_flat = IntegrateSampling("yy[WW]", decay, "C", true, 32999, 2 * calls, 650.0, 650.0, 155.0, 300.0);
  RequireSamplingAgreement(factorized_bw, central_bw);
  RequireSamplingAgreement(factorized_bw, factorized_flat);
  RequireSamplingAgreement(central_bw, central_flat);
  RequireSamplingAgreement(factorized_flat, central_flat);
}
