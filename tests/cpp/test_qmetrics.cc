// Physical spin metrics, weighted ensembles and worker-local integration tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <array>
#include <catch.hpp>
#include <cmath>
#include <complex>
#include <limits>
#include <numeric>
#include <thread>
#include <vector>

#include "Graniitti/MGraniitti.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MSpin.h"
#include "Graniitti/Spin/MQMetrics.h"
#include "Graniitti/Tech/MException.h"
#include "support/models_test_support.hh"

namespace {

using gra::aux::indices;
using Density = gra::spin::MQMetrics::Density;
using Complex = std::complex<double>;

// Define a pair of physical spin spaces in negative-first helicity order
gra::spin::SpinSpec Pair(const int spin1, const int spin2) {
  gra::spin::SpinSpec spec;
  spec.name    = "final_pair";
  spec.type    = gra::spin::SpinType::FinalPair;
  spec.frame   = "CM";
  spec.measure = "analytic spin ensemble";
  for (int m = -spin1; m <= spin1; m += 2) { spec.helicity[0].push_back(m); }
  for (int m = -spin2; m <= spin2; m += 2) { spec.helicity[1].push_back(m); }
  return spec;
}

// Integrate an analytic amplitude without changing its physical phases
gra::spin::SpinResult Integrate(const gra::spin::SpinSpec &spec, const Density &amplitude) {
  gra::spin::MQMetrics metrics;
  metrics.Configure({spec}, 4, spec.Dimension());
  for (std::uint64_t i = 0; i < 16; ++i) {
    metrics.Begin(true);
    metrics.Add(0, amplitude, 1.0);
    metrics.Observe(1.0, 0.0, i);
  }
  return metrics.Results().at(0);
}

// Compare every reported metric and the density normalization
void Check(const gra::spin::SpinResult &result, const std::array<double, 6> &expected) {
  REQUIRE(result.status == gra::spin::SpinStatus::Ready);
  REQUIRE(result.rho.Trace().real() == Approx(1.0).margin(1e-12));
  REQUIRE((result.rho - result.rho.Dagger()).FrobNorm() < 1e-12);
  for (const auto i : indices(expected)) {
    CAPTURE(i);
    CHECK(result.value[i] == Approx(expected[i]).margin(1e-11));
  }
}

// Check that the proton extensions preserve the central ensemble and its integrated rate
void CheckProtons(const std::vector<gra::spin::SpinResult> &results) {
  const std::size_t count = results.size() / 3;
  REQUIRE(results.size() == 3 * count);
  for (std::size_t i = 0; i < count; ++i) {
    const auto &central = results[i];
    const auto &joint = results[count + 2 * i];
    const auto &pp = results[count + 2 * i + 1];
    CAPTURE(central.spec.name);
    REQUIRE(central.status == gra::spin::SpinStatus::Ready);
    REQUIRE(joint.status == gra::spin::SpinStatus::Ready);
    REQUIRE(pp.status == gra::spin::SpinStatus::Ready);
    const auto n = central.spec.Dimension();
    CHECK(joint.rho.size_row() == 4 * n);
    CHECK(pp.rho.size_row() == 4);
    CHECK((joint.rho.PartialTrace(4, n, gra::TensorFactor::First) - central.rho).FrobNorm() < 1e-11);
    CHECK((joint.rho.PartialTrace(4, n, gra::TensorFactor::Second) - pp.rho).FrobNorm() < 1e-11);
    CHECK(joint.integral == Approx(central.integral).epsilon(1e-11));
    CHECK(pp.integral == Approx(central.integral).epsilon(1e-11));
  }
}

// Integrate canonical proton transitions and the same central spin space
std::vector<gra::spin::SpinResult> ProtonMetrics(const Density &amplitude, const std::size_t rows) {
  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 16, "", {2212, 2212});
  for (std::uint64_t i = 0; i < 16; ++i) {
    metrics.Begin(true);
    metrics.Add(0, amplitude, 0.25, rows);
    metrics.Observe(1.0, 0.0, i);
  }
  const auto results = metrics.Results();
  CheckProtons(results);
  return results;
}

// Construct a particle with analytic quantum numbers independently of a model tune
gra::MParticle Particle(const int pdg, const int spin, const int parity) {
  gra::MParticle p;
  p.pdg    = pdg;
  p.spinX2 = spin;
  p.P      = parity;
  p.C      = parity;
  p.name   = "spin test";
  p.mass   = 1.0;
  return p;
}

// Compute Shannon entropy including exactly vanishing probabilities
double Entropy(const std::vector<double> &probability) {
  double value = 0.0;
  for (const double p : probability) {
    if (p > 0.0) { value -= p * std::log2(p); }
  }
  return value;
}

}  // namespace

// A product state and a complex Bell pair have different entanglement despite equal purity
TEST_CASE("Spin metrics recover product and Bell states with their phases", "[qmetrics]") {
  const auto   spec = Pair(1, 1);
  const double c    = 1.0 / std::sqrt(2.0);
  Density      product{{1.0}, {0.0}, {0.0}, {0.0}};
  Check(Integrate(spec, product), {1.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  Density    bell{{c}, {0.0}, {0.0}, {Complex(0.0, c)}};
  const auto result = Integrate(spec, bell);
  Check(result, {1.0, 0.0, 1.0, 1.0, 0.5, 1.0});
  CHECK(std::abs(result.rho[0][3] - Complex(0.0, -0.5)) < 1e-12);
  CHECK(result.errors);
  for (const double error : result.error) { CHECK(error < 1e-11); }

  const auto rotation = gra::spin::SpinHalfRotation(0.71, -0.32).Kronecker(gra::spin::SpinHalfRotation(-0.43, 0.29));
  Check(Integrate(spec, rotation * bell), result.value);
  const std::vector<std::size_t> exchange = {0, 2, 1, 3};
  Check(Integrate(spec, bell.SelectRows(exchange)), result.value);
}

// An unpolarized spin-one parent has three equal eigenvalues in any spin frame
TEST_CASE("Spin metrics recover an unpolarized production density", "[qmetrics]") {
  auto spec            = Pair(2, 0);
  spec.type            = gra::spin::SpinType::Production;
  const auto amplitude = Density::IdentityMatrix(3) / std::sqrt(3.0);
  Check(Integrate(spec, amplitude), {1.0 / 3.0, std::log2(3.0), 0.0, 0.0, 0.0, 0.0});
  const auto rotation = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.53, -0.87), 1.0);
  Check(Integrate(spec, rotation * amplitude), {1.0 / 3.0, std::log2(3.0), 0.0, 0.0, 0.0, 0.0});
}

// The two-qubit Werner family has an analytic negativity threshold at p=1/3
TEST_CASE("Spin metrics recover the Werner spectrum and separability threshold", "[qmetrics]") {
  for (const double p : {0.0, 0.2, 1.0 / 3.0, 0.5, 1.0}) {
    CAPTURE(p);
    Density amplitude(4, 5, 0.0);
    amplitude[1][0] = std::sqrt(p / 2.0);
    amplitude[2][0] = -std::sqrt(p / 2.0);
    for (std::size_t i = 0; i < 4; ++i) { amplitude[i][i + 1] = std::sqrt((1.0 - p) / 4.0); }
    const double n = std::max(0.0, (3.0 * p - 1.0) / 4.0);
    const double a = (1.0 + 3.0 * p) / 4.0, b = (1.0 - p) / 4.0;
    Check(Integrate(Pair(1, 1), amplitude),
          {(1.0 + 3.0 * p * p) / 4.0, Entropy({a, b, b, b}), 1.0, 1.0, n, std::log2(1.0 + 2.0 * n)});
  }
}

// Physical Jacob-Wick decay amplitudes reproduce scalar and tensor vector-pair densities
TEST_CASE("Spin metrics recover scalar vector singlets and unpolarized tensor decays", "[qmetrics]") {
  const auto vector = Particle(113, 2, -1);
  for (const int spin : {0, 4}) {
    const auto     mother = Particle(9000901 + spin, spin, 1);
    gra::HELMatrix hel;
    hel.alpha_ls.Set(0, spin, 1.0);
    gra::spin::InitTMatrix(hel, mother, vector, vector, false, "spin metrics", false, false);
    const auto amplitude = gra::spin::fDecayMatrix(hel, 0.63, -0.41) / std::sqrt(spin + 1.0);
    const auto result    = Integrate(Pair(2, 2), amplitude);
    if (spin == 0) {
      Check(result, {1.0, 0.0, std::log2(3.0), std::log2(3.0), 1.0, std::log2(3.0)});
      CHECK(result.rho[0][4].real() == Approx(-1.0 / 3.0).margin(1e-12));
    } else {
      Check(result, {0.2, std::log2(5.0), std::log2(3.0), std::log2(3.0), 0.0, 0.0});
    }
    const auto rotation = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.32, -0.14), 1.0);
    Check(Integrate(Pair(2, 2), rotation.Kronecker(rotation) * amplitude), result.value);
  }
}

// Coherent S and D waves select a separable longitudinal vector pair
TEST_CASE("Spin metrics retain LS interference and transverse photon states", "[qmetrics]") {
  const auto     scalar = Particle(9000901, 0, 1);
  const auto     vector = Particle(113, 2, -1);
  gra::HELMatrix hel;
  hel.alpha_ls.Set(0, 0, 1.0);
  hel.alpha_ls.Set(2, 4, -std::sqrt(2.0));
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "longitudinal pair", false, false);
  Check(Integrate(Pair(2, 2), gra::spin::fDecayMatrix(hel, 0.0, 0.0)), {1.0, 0.0, 0.0, 0.0, 0.0, 0.0});

  hel.alpha_ls.Clear();
  hel.alpha_ls.Set(2, 4, 1.0);
  gra::spin::InitTMatrix(hel, scalar, vector, vector, false, "D wave", false, false);
  const double s = Entropy({1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0});
  Check(Integrate(Pair(2, 2), gra::spin::fDecayMatrix(hel, 0.0, 0.0)),
        {1.0, 0.0, s, s, 5.0 / 6.0, std::log2(8.0 / 3.0)});

  auto photon = Particle(22, 2, -1);
  photon.mass = 0.0;
  photon.C    = -1;
  hel.alpha_ls.Clear();
  hel.alpha_ls.Set(0, 0, 1.0);
  gra::spin::InitTMatrix(hel, scalar, photon, photon, false, "photon pair", false, false);
  auto spec     = Pair(1, 1);
  spec.helicity = {std::vector<int>{-2, 2}, std::vector<int>{-2, 2}};
  Check(Integrate(spec, gra::spin::fDecayMatrix(hel, 0.0, 0.0)), {1.0, 0.0, 1.0, 1.0, 0.5, 1.0});
}

// Weighted density integration must precede nonlinear metrics and commute with worker merging
TEST_CASE("Spin ensembles retain importance weights across workers and serialization", "[qmetrics]") {
  gra::spin::MQMetrics serial;
  const bool protons = GENERATE(false, true);
  serial.Configure({Pair(1, 1)}, 4, 16, "", protons ? std::array<int, 2>{2212, 2212} : std::array<int, 2>{});
  std::array<gra::spin::MQMetrics, 3> workers = {serial, serial, serial};
  const auto                             sample  = [protons](gra::spin::MQMetrics &metrics, std::uint64_t i) {
    Density    amplitude(4, 1, 0.0);
    const bool first            = i % 3 == 0;
    amplitude[first ? 0 : 3][0] = first ? std::sqrt(2.0) : 1.0;
    metrics.Begin(true);
    if (protons) { amplitude = Density(4, 1, 0.5).Kronecker(amplitude); }
    metrics.Add(0, amplitude, 1.0, 4);
    metrics.Observe(3.0, -std::log(first ? 2.0 / 3.0 : 4.0 / 3.0), i);
  };
  std::vector<std::thread> threads;
  for (const auto w : indices(workers)) {
    threads.emplace_back([&, w] {
      for (std::uint64_t k = 81; k > 0; --k) { sample(workers[w], w + (k - 1) * workers.size()); }
    });
  }
  for (std::uint64_t i = 0; i < 243; ++i) { sample(serial, i); }
  for (auto &thread : threads) { thread.join(); }
  auto merged = workers[0];
  merged.Merge(workers[1]);
  merged.Merge(workers[2]);
  if (protons) { CheckProtons(merged.Results()); }
  const auto serial_results = serial.Results();
  const auto merged_results = merged.Results();
  for (const auto i : indices(serial_results)) {
    const auto &a = serial_results[i];
    const auto &b = merged_results[i];
    CHECK((a.rho - b.rho).FrobNorm() < 1e-12);
    Check(b, a.value);
    for (const auto j : indices(a.error)) { CHECK(b.error[j] == Approx(a.error[j]).margin(1e-12)); }
  }
  const auto expected = serial.Results().at(0);
  Check(expected, {5.0 / 9.0, Entropy({2.0 / 3.0, 1.0 / 3.0}), Entropy({2.0 / 3.0, 1.0 / 3.0}),
                   Entropy({2.0 / 3.0, 1.0 / 3.0}), 0.0, 0.0});
  CHECK(expected.integral == Approx(4.5).margin(1e-12));
  REQUIRE(expected.errors);
  CHECK(expected.error[0] == Approx(0.0).margin(1e-12));
  auto restored = serial;
  restored.Reset();
  restored.Deserialize(merged.Serialize());
  const auto restored_results = restored.Results();
  for (const auto i : indices(restored_results)) {
    CHECK((restored_results[i].rho - merged_results[i].rho).FrobNorm() < 1e-12);
  }
  for (const auto &metrics : {merged, restored}) {
    if (protons) { CheckProtons(metrics.Results()); }
    const auto result = metrics.Results().at(0);
    Check(result, expected.value);
    REQUIRE(result.samples == expected.samples);
    REQUIRE((result.rho - expected.rho).FrobNorm() < 1e-12);
    for (const auto i : indices(result.error)) { CHECK(result.error[i] == Approx(expected.error[i]).margin(1e-12)); }
  }
  auto invalid                          = merged.Serialize();
  invalid["specs"][0]["helicity_x2"][0] = {0, 2};
  REQUIRE_THROWS_AS(restored.Deserialize(invalid), std::invalid_argument);
}

// Unobserved row and column indices must both be traced, and nonlinear errors follow the mixture
TEST_CASE("Spin metrics trace spectators and recover analytic jackknife errors", "[qmetrics]") {
  Density amplitude(8, 2, 0.0);
  amplitude[0][0]      = std::sqrt(1.0 / 3.0);
  amplitude[7][1]      = std::sqrt(2.0 / 3.0);
  const double entropy = Entropy({1.0 / 3.0, 2.0 / 3.0});
  Check(Integrate(Pair(1, 1), amplitude), {5.0 / 9.0, entropy, entropy, entropy, 0.0, 0.0});

  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 4);
  for (std::uint64_t i = 0; i < 16; ++i) {
    Density state(4, 1, 0.0);
    state[i / 4 <= i % 4 ? 0 : 3][0] = 1.0;
    metrics.Begin(true);
    metrics.Add(0, state, 1.0);
    metrics.Observe(1.0, 0.0, i);
  }
  const auto   result = metrics.Results().at(0);
  const double s      = Entropy({5.0 / 8.0, 3.0 / 8.0});
  Check(result, {17.0 / 32.0, s, s, s, 0.0, 0.0});
  REQUIRE(result.errors);
  std::array<double, 4> purity;
  for (const auto g : indices(purity)) {
    const double p = (9.0 - g) / 12.0;
    purity[g]      = p * p + (1.0 - p) * (1.0 - p);
  }
  const double mean     = std::accumulate(purity.begin(), purity.end(), 0.0) / 4.0;
  double       variance = 0.0;
  for (const double p : purity) { variance += 0.75 * (p - mean) * (p - mean); }
  CHECK(result.error[0] == Approx(std::sqrt(variance)).margin(1e-12));
}

// Unequal spin spaces expose partial-trace ordering and complex rotation conventions
TEST_CASE("Spin metrics preserve unequal spin spaces and local covariance", "[qmetrics]") {
  Density amplitude(6, 3, 0.0);
  amplitude[0][0]     = std::sqrt(0.5);
  amplitude[1][1]     = 0.5;
  amplitude[5][2]     = 0.5;
  const auto   result = Integrate(Pair(1, 2), amplitude);
  const double s      = Entropy({0.75, 0.25});
  Check(result, {0.375, 1.5, s, 1.5, 0.0, 0.0});
  const auto rotation = gra::spin::SpinHalfRotation(0.71, -0.32)
                            .Kronecker(gra::spin::SpinRotation(gra::spin::SpinHalfRotation(-0.43, 0.29), 1.0));
  const auto rotated  = Integrate(Pair(1, 2), rotation * amplitude);
  Check(rotated, result.value);
  CHECK((rotated.rho - rotation * result.rho * rotation.Dagger()).FrobNorm() < 1e-12);
  Check(Integrate(Pair(2, 1), amplitude.SelectRows({0, 3, 1, 4, 2, 5})), {0.375, 1.5, 1.5, s, 0.0, 0.0});
}

// Opposite Bell phases cancel in the ensemble density but not in averaged entanglement
TEST_CASE("Spin metrics distinguish coherent blocks from mixed Bell pairs", "[qmetrics]") {
  const double  c = 1.0 / std::sqrt(2.0);
  const Density plus{{c}, {0.0}, {0.0}, {Complex(0.0, c)}};
  const Density minus = plus.Conj();
  Check(Integrate(Pair(1, 1), plus + minus), {1.0, 0.0, 0.0, 0.0, 0.0, 0.0});
  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 4);
  for (std::uint64_t i = 0; i < 16; ++i) {
    metrics.Begin(true);
    metrics.Add(0, plus, 0.5);
    metrics.Add(0, minus, 0.5);
    metrics.Observe(1.0, 0.0, i);
  }
  Check(metrics.Results().at(0), {0.5, 1.0, 1.0, 1.0, 0.0, 0.0});
}

// A dominant weight must not erase smaller groups when jackknife replicas are formed
TEST_CASE("Spin jackknife retains support across extreme weights", "[qmetrics]") {
  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 4);
  const double  c = 1.0 / std::sqrt(2.0);
  const Density bell{{c}, {0.0}, {0.0}, {Complex(0.0, c)}};
  for (std::uint64_t i = 0; i < 4; ++i) {
    Density state(4, 1, 0.0);
    state[i == 3 ? 3 : 0][0] = 1.0;
    metrics.Begin(true);
    metrics.Add(0, i == 0 ? bell : state, 1.0);
    metrics.Observe(i == 0 ? 1e30 : 1.0, 0.0, i);
  }
  const auto result = metrics.Results().at(0);
  Check(result, {1.0, 0.0, 1.0, 1.0, 0.5, 1.0});
  REQUIRE(result.errors);
  const double                s     = Entropy({2.0 / 3.0, 1.0 / 3.0});
  const std::array<double, 6> error = {1.0 / 3.0, 0.75 * s, 0.75 * (1.0 - s), 0.75 * (1.0 - s), 0.375, 0.75};
  for (const auto i : indices(error)) { CHECK(result.error[i] == Approx(error[i]).margin(1e-11)); }

  // Discard only the last extra trial when balancing unequal groups, even when it dominates
  metrics.Reset();
  for (std::uint64_t i = 0; i < 5; ++i) {
    metrics.Begin(true);
    metrics.Add(0, bell, 1.0);
    metrics.Observe(i == 4 ? 1e30 : 1.0, 0.0, i);
  }
  const auto balanced = metrics.Results().at(0);
  Check(balanced, result.value);
  REQUIRE(balanced.errors);
  for (const double error : balanced.error) { CHECK(error < 1e-11); }
}

// Reloading invalid counts or tiny unphysical densities must leave the existing state untouched
TEST_CASE("Spin metrics reject corrupt integration state atomically", "[qmetrics]") {
  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 4);
  const Density state{{1.0}, {0.0}, {0.0}, {0.0}};
  metrics.Begin(true);
  metrics.Add(0, state, 1.0);
  metrics.Observe(1.0, 0.0, 0);
  const auto saved = metrics.Serialize();
  for (int test = 0; test < 7; ++test) {
    CAPTURE(test);
    auto invalid = saved;
    if (test == 0) { invalid["failures"][0] = -1; }
    if (test == 1) { invalid["failures"][0] = std::uint64_t{2}; }
    if (test == 2) { invalid["sums"][0]["count"] = std::uint64_t{2}; }
    if (test == 3) { invalid["sums"][0]["last_sample"] = std::uint64_t{1}; }
    if (test == 4) { invalid["sums"][0]["last"][0][0] = {-1e-200, 0.0}; }
    if (test == 5) {
      invalid["sums"][0]["last"][0][1] = {1e-200, 1e-200};
      invalid["sums"][0]["last"][0][0] = {1e-200, 0.0};
    }
    if (test == 6) { invalid["sums"][0]["rest"] = invalid["sums"][0]["last"]; }
    REQUIRE_THROWS(metrics.Deserialize(invalid));
    CHECK(metrics.Serialize() == saved);
  }
}

// Disabled trials and rejected trials cannot reuse a previous amplitude
TEST_CASE("Spin metrics separate disabled zero and failed trials", "[qmetrics]") {
  gra::spin::MQMetrics metrics;
  metrics.Configure({Pair(1, 1)}, 4, 4);
  const Density state{{1.0}, {0.0}, {0.0}, {0.0}};
  metrics.Begin(false);
  metrics.Add(0, state, 1.0);
  metrics.Observe(1.0, 0.0, 0);
  CHECK(metrics.Results()[0].samples == 0);
  metrics.Begin(true);
  metrics.Add(0, state, 1.0);
  metrics.Observe(0.0, 0.0, 0);
  CHECK(metrics.Results()[0].status == gra::spin::SpinStatus::Empty);
  metrics.Begin(true);
  metrics.Observe(1.0, 0.0, 1);
  CHECK(metrics.Results()[0].status == gra::spin::SpinStatus::Incomplete);
  CHECK(metrics.Results()[0].failures == 1);
  metrics.Reset();
  Density invalid = state;
  invalid[0][0]   = std::numeric_limits<double>::quiet_NaN();
  metrics.Begin(true);
  REQUIRE_NOTHROW(metrics.Add(0, invalid, 1.0));
  REQUIRE_NOTHROW(metrics.Observe(1.0, 0.0, 0));
  CHECK(metrics.Results()[0].status == gra::spin::SpinStatus::Incomplete);
  CHECK(metrics.Results()[0].failures == 1);
  metrics.Reset();
  metrics.Begin(true);
  metrics.Add(0, {}, 1.0);
  metrics.Observe(1.0, 0.0, 0);
  CHECK(metrics.Results()[0].status == gra::spin::SpinStatus::Incomplete);
  metrics.Configure({Pair(0, 1)}, 4, 4);
  metrics.Begin(true);
  metrics.Add(0, Density::IdentityMatrix(2) * std::sqrt(0.75 * std::numeric_limits<double>::max()), 1.0);
  metrics.Observe(1.0, 0.0, 0);
  CHECK(metrics.Results()[0].status == gra::spin::SpinStatus::Incomplete);
  REQUIRE_THROWS_AS(metrics.Configure({Pair(1, 1)}, 1, 4), std::invalid_argument);
  REQUIRE_THROWS_AS(metrics.Configure({Pair(2, 2)}, 4, 4), std::invalid_argument);
}

// Diagnostic numerical failures must not change the acceptance of a finite physical amplitude
TEST_CASE("Spin diagnostic failures preserve the physical amplitude boundary", "[qmetrics]") {
  // Expose the actual protected validation method without overriding any physics
  struct Process : gra::MFactorized {
    using gra::MProcess::ValidateAmplitude;
  };
  Process process;
  auto   &lts                      = process.state.lts;
  lts.hamp                         = {{1.0, 0.2}};
  lts.hamp.metadata.amplitude_type = gra::ScreeningAmplitudeType::GoodWalker;
  lts.proton_good_walker.emplace();
  auto &components = lts.proton_good_walker->components;
  components.emplace_back();
  components.back().source = Density(1, 1, Complex(1.0, 0.2));
  components.emplace_back();
  components.back().spin_index = 0;
  components.back().source     = Density(1, 1, Complex(std::numeric_limits<double>::quiet_NaN(), 0.0));
  REQUIRE_NOTHROW(process.ValidateAmplitude(1.04));
  components.back().spin_index.reset();
  REQUIRE_THROWS_AS(process.ValidateAmplitude(1.04), gra::AmplitudeFailure);
}

// Proton-central coherence survives while distinct incoming helicities remain incoherent
TEST_CASE("Proton metrics retain complex entanglement and incoming spin mixtures", "[qmetrics][proton]") {
  Density amplitude(64, 1, 0.0);
  amplitude[0][0] = 1.0;
  amplitude[10][0] = Complex(0.0, 1.0);
  const auto pure = ProtonMetrics(amplitude, 16);
  Check(pure[1], {1.0, 0.0, 1.0, 1.0, 0.5, 1.0});
  Check(pure[2], {0.5, 1.0, 1.0, 0.0, 0.0, 0.0});
  CHECK(std::abs(pure[1].rho[0][10] - Complex(0.0, -0.5)) < 1e-12);

  // Opposite phases in a different incoming state must cancel without amplitude interference
  amplitude[16][0] = 1.0;
  amplitude[26][0] = Complex(0.0, -1.0);
  const auto mixed = ProtonMetrics(amplitude, 16);
  Check(mixed[1], {0.5, 1.0, 1.0, 1.0, 0.0, 0.0});
  CHECK(mixed[1].integral == Approx(2.0 * pure[1].integral));

  const auto rotation = gra::spin::SpinHalfRotation(0.71, -0.32)
                            .Kronecker(gra::spin::SpinHalfRotation(-0.43, 0.29))
                            .Kronecker(gra::spin::SpinHalfRotation(0.22, 0.61))
                            .Kronecker(gra::spin::SpinHalfRotation(-0.52, -0.11));
  const auto rotated = ProtonMetrics(Density::IdentityMatrix(4).Kronecker(rotation) * amplitude, 16);
  Check(rotated[1], mixed[1].value);
  Check(rotated[2], mixed[2].value);
  CHECK((rotated[1].rho - rotation * mixed[1].rho * rotation.Dagger()).FrobNorm() < 1e-12);

  const auto exchange = Density::IdentityMatrix(4).SelectRows({0, 2, 1, 3}).Kronecker(Density::IdentityMatrix(4));
  const auto swapped = ProtonMetrics(Density::IdentityMatrix(4).Kronecker(exchange) * amplitude, 16);
  Check(swapped[1], mixed[1].value);
  CHECK(swapped[2].value[2] == Approx(mixed[2].value[3]).margin(1e-12));
  CHECK(swapped[2].value[3] == Approx(mixed[2].value[2]).margin(1e-12));
}

// Compact no-flip storage and dense storage represent the same mixed proton state
TEST_CASE("Proton metrics expand no-flip rows and trace all spectators", "[qmetrics][proton]") {
  Density compact(32, 2, 0.0), dense(128, 2, 0.0);
  for (std::size_t i = 0; i < 4; ++i) {
    compact[i * 8][0] = dense[(i * 4 + i) * 8][0] = 1.0;
    compact[i * 8 + 3][0] = dense[(i * 4 + i) * 8 + 3][0] = Complex(0.0, 1.0);
    compact[i * 8 + 5][1] = dense[(i * 4 + i) * 8 + 5][1] = 0.5;
  }
  const auto a = ProtonMetrics(compact, 4);
  const auto b = ProtonMetrics(dense, 16);
  for (const auto i : indices(a)) {
    CHECK((a[i].rho - b[i].rho).FrobNorm() < 1e-12);
    CHECK(a[i].integral == Approx(b[i].integral));
  }
  Check(a[2], {0.25, 2.0, 1.0, 1.0, 0.0, 0.0});
  CHECK((a[1].rho - (Density::IdentityMatrix(4) / 4.0).Kronecker(a[0].rho)).FrobNorm() < 1e-12);
  CHECK(a[1].value[4] == Approx(0.0).margin(1e-12));
}

// Invalid proton dimensions remain diagnostic failures and cannot contaminate later events
TEST_CASE("Proton metrics bound storage and isolate invalid amplitudes", "[qmetrics][proton]") {
  gra::spin::MQMetrics metrics;
  REQUIRE_THROWS_AS(metrics.Configure({Pair(1, 1)}, 4, 15, "", {2212, 2212}), std::invalid_argument);
  metrics.Configure({Pair(1, 1)}, 4, 16, "", {2212, 2212});
  metrics.Begin(true);
  REQUIRE_NOTHROW(metrics.Add(0, Density(16, 1, 1.0), 0.25, 3));
  metrics.Observe(1.0, 0.0, 0);
  for (const auto &r : metrics.Results()) { CHECK(r.status == gra::spin::SpinStatus::Incomplete); }
  metrics.Reset();
  metrics.Begin(true);
  metrics.Add(0, Density(16, 1, 1.0), 0.25, 4);
  metrics.ClearEvent();
  metrics.Add(0, Density(16, 1, 2.0), 0.25, 4);
  metrics.Observe(1.0, 0.0, 0);
  CheckProtons(metrics.Results());
  CHECK(metrics.Results()[1].integral == Approx(16.0));
}

// The same scalar S-wave decay must remain a pure vector singlet in every production model
TEST_CASE("Regge metrics preserve amplitudes and recover scalar vector entanglement", "[qmetrics][regge]") {
  const std::string channel = GENERATE("RES", "CON", "RES+CON");
  const int decays = GENERATE(0, 1, 2);
  const std::string final = std::string("rho(770)0") + (decays > 0 ? " > {pi+ pi-}" : "") +
                            " rho(770)0" + (decays > 1 ? " > {pi+ pi-}" : "");
  CAPTURE(channel, decays);
  for (const std::string model : {"MP", "XP", "GP"}) {
    CAPTURE(model);
    auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
    card["GENERIC"]["HIST"]          = 0;
    card["GENERIC"]["INTEGRATOR"]    = "VEGAS";
    card["SCATTERING"]["PROCESS"]    = model + "[" + channel + "]<F> -> " + final + " @QMETRICS:true";
    card["SCATTERING"]["RES"]        = {"f0_1710"};
    card["SCATTERING"]["ENERGY"]     = {50.0, 50.0};
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["GENCUTS"]["<F>"]["M"]      = {2.0, 3.0};
    card["GENCUTS"]["<F>"]["Rap"]    = {-1.0, 1.0};
    card["GENCUTS"]["<F>"]["Pt"]     = {0.1, 0.5};
    card["FIDCUTS"]["active"]        = false;
    card["VETOCUTS"]["active"]       = false;
    gra::MGraniitti generator;
    generator.ReadInput(card);
    auto &process = dynamic_cast<gra::MFactorized &>(*generator.proc);
    process.PrepareRun();
    REQUIRE(process.state.lts.process.QMETRICS);
    if (channel != "CON") {
      auto &res = process.state.lts.process.RESONANCES.at("f0_1710");
      res.hel_decay.alpha_ls.Clear();
      res.hel_decay.alpha_ls.Set(0, 0, 1.0);
      gra::spin::InitTMatrix(res.hel_decay, res.p, process.state.lts.decaytree[0].p, process.state.lts.decaytree[1].p,
                             false, "scalar S wave", false, false);
    }
    card["SCATTERING"]["PROCESS"] = model + "[" + channel + "]<F> -> " + final + " @QMETRICS:false";
    gra::MGraniitti disabled;
    disabled.ReadInput(card);
    auto &reference = dynamic_cast<gra::MFactorized &>(*disabled.proc);
    reference.PrepareRun();
    if (channel != "CON") {
      reference.state.lts.process.RESONANCES.at("f0_1710").hel_decay =
          process.state.lts.process.RESONANCES.at("f0_1710").hel_decay;
    }
    REQUIRE_FALSE(reference.state.lts.process.QMETRICS);
    REQUIRE_FALSE(reference.state.lts.qmetrics.Configured());
    REQUIRE(reference.state.lts.qmetrics.Specs().empty());

    gra::MRandom random;
    random.SetSeed(82475);
    gra::VEGASPARAM param;
    param.BINS                  = 8;
    param.ROUNDS                = 1;
    param.AUTOMATIC_CONVERGENCE = false;
    gra::MVEGASIntegrator vegas;
    vegas.SetDimension(process.GetdLIPSDim());
    vegas.Init(gra::VEGASStage::Adaptation, 128, 1, param);
    const auto evaluate = [&](const gra::VEGASSample &sample, bool integration, bool screening, std::uint64_t index) {
      gra::MEventWeightState aux;
      aux.include_screening   = screening;
      aux.adaptation_mode     = !integration;
      aux.qmetrics        = integration;
      aux.sample_index        = index;
      aux.log_inverse_density = sample.log_inverse_density;
      auto off_aux            = aux;
      off_aux.qmetrics    = false;
      const double off        = reference.EventWeight(sample.point, off_aux);
      const double on         = process.EventWeight(sample.point, aux);
      REQUIRE_FALSE(aux.technical_failure);
      REQUIRE_FALSE(off_aux.technical_failure);
      CHECK(on == Approx(off).epsilon(1e-12).margin(1e-30));
      REQUIRE_FALSE(reference.state.lts.qmetrics.Active());
      REQUIRE(reference.state.lts.qmetrics.Results().empty());
      if (reference.state.lts.proton_good_walker) {
        for (const auto &component : reference.state.lts.proton_good_walker->components) {
          REQUIRE_FALSE(component.spin_index.has_value());
        }
      }
      const auto &a = process.state.lts.hamp;
      const auto &b = reference.state.lts.hamp;
      REQUIRE(a.size() == b.size());
      for (const auto h : indices(a)) {
        CHECK(std::abs(a[h] - b[h]) <= 1e-12 * std::max(std::abs(a[h]), std::abs(b[h])) + 1e-30);
      }
      return gra::statistics::ImportanceWeight(on, sample.log_inverse_density);
    };
    while (!vegas.IsFrozen()) {
      const auto calls = vegas.NextAdaptationCalls();
      vegas.BeginAdaptationBatch();
      for (std::uint64_t i = 0; i < calls; ++i) {
        const auto sample = vegas.Sample([&]() { return random.U(0.0, 1.0); });
        vegas.AccumulateAdaptation(evaluate(sample, false, false, i), sample.indices);
      }
      vegas.FinishAdaptationBatch(calls);
    }
    REQUIRE(process.state.lts.qmetrics.Results()[0].samples == 0);
    vegas.Init(gra::VEGASStage::Integration, 32, 1, param);
    for (const bool screening : {false, true}) {
      CAPTURE(screening);
      process.state.lts.qmetrics.Reset();
      double integral = 0.0;
      for (std::uint64_t i = 0; i < 32; ++i) {
        const auto sample = vegas.Sample([&]() { return random.U(0.0, 1.0); });
        integral += evaluate(sample, true, screening, i) / 32.0;
      }
      const auto results = process.state.lts.qmetrics.Results();
      REQUIRE(results.size() == (channel == "CON" ? 3 : 6));
      CheckProtons(results);
      REQUIRE(integral > 0.0);
      REQUIRE(results[0].status == gra::spin::SpinStatus::Ready);
      CHECK(results[0].rho.size_row() == 9);
      CHECK(results[0].rho.Trace().real() == Approx(1.0).margin(1e-12));
      if (channel == "RES") {
        Check(results[0], {1.0, 0.0, std::log2(3.0), std::log2(3.0), 1.0, std::log2(3.0)});
        CHECK(results[0].integral == Approx(results[1].integral).epsilon(1e-11));
      }
      if (channel != "CON") { Check(results[1], {1.0, 0.0, 0.0, 0.0, 0.0, 0.0}); }
      if (decays == 0) {
        CHECK(results[0].integral == Approx(integral).epsilon(1e-11));
      } else {
        CHECK(results[0].spec.type == gra::spin::SpinType::IntermediatePair);
        CHECK_FALSE(process.state.lts.qmetrics.FinalPair());
        auto restored = process.state.lts.qmetrics;
        restored.Deserialize(process.state.lts.qmetrics.Serialize());
        CHECK((restored.Results()[0].rho - results[0].rho).FrobNorm() < 1e-12);
      }
      CHECK(process.state.random.U(0.0, 1.0) == Approx(reference.state.random.U(0.0, 1.0)).margin(1e-15));
    }
  }
}

// Sewing a spin-two production density to its decay must reproduce the full pair density
TEST_CASE("Regge spin densities obey the production decay relation", "[qmetrics][regge]") {
  const auto tune = WriteModifiedPhotoVMTune("qmetrics_flip", [](auto &general) {
    general["PARAM_SPIN"]["FORWARD_NOFLIP"] = false;
  });
  const auto model_tune = gra::MModelTune::Load(tune.second);
  const bool cascade = GENERATE(false, true);
  const std::string final = cascade ? "rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}" : "rho(770)0 rho(770)0";
  CAPTURE(cascade);
  for (const std::string model : {"MP", "XP", "GP"}) {
    CAPTURE(model);
    auto card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
    card["GENERIC"]["HIST"]          = 0;
    card["GENERIC"]["INTEGRATOR"]    = "VEGAS";
    card["SCATTERING"]["PROCESS"]    = model + "[RES]<F> -> " + final + " @QMETRICS:true";
    card["SCATTERING"]["RES"]        = {"f2_2150"};
    card["SCATTERING"]["ENERGY"]     = {50.0, 50.0};
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["GENCUTS"]["<F>"]["M"]      = {2.0, 3.0};
    card["GENCUTS"]["<F>"]["Rap"]    = {-1.0, 1.0};
    card["GENCUTS"]["<F>"]["Pt"]     = {0.1, 0.5};
    card["FIDCUTS"]["active"]        = false;
    card["VETOCUTS"]["active"]       = false;
    gra::MGraniitti generator;
    generator.ReadInput(card);
    auto &process = *generator.proc;
    process.SetModelTune(model_tune);
    process.PrepareRun();
    auto &lts = process.state.lts;
    auto &res = lts.process.RESONANCES.at("f2_2150");
    res.hel_decay.alpha_ls.Clear();
    res.hel_decay.alpha_ls.Set(0, 4, 1.0);
    gra::spin::InitTMatrix(res.hel_decay, res.p, lts.decaytree[0].p, lts.decaytree[1].p, false, "tensor S wave", false,
                           false);
    gra::MRandom random;
    random.SetSeed(382541);
    for (const bool screening : {false, true}) {
      CAPTURE(screening);
      std::size_t accepted = 0;
      for (std::uint64_t i = 0; i < 100 && accepted < 3; ++i) {
        std::vector<double> point(process.GetdLIPSDim());
        for (auto &x : point) { x = random.U(0.0, 1.0); }
        lts.qmetrics.Reset();
        gra::MEventWeightState aux;
        aux.qmetrics      = true;
        aux.include_screening = screening;
        const double weight   = process.EventWeight(point, aux);
        REQUIRE_FALSE(aux.technical_failure);
        if (!aux.Valid()) { continue; }
        ++accepted;
        const auto results = lts.qmetrics.Results();
        REQUIRE(results.size() == 6);
        CheckProtons(results);
        const auto &pair   = results[0];
        const auto &parent = results[1];
        REQUIRE(pair.status == gra::spin::SpinStatus::Ready);
        REQUIRE(parent.status == gra::spin::SpinStatus::Ready);
        REQUIRE(parent.rho.size_row() == 5);
        auto decay = gra::spin::ResonanceDecayMatrix(lts, res, res.hel_decay, "CM", !cascade);
        if (cascade) {
          auto direct = lts;
          for (auto &branch : direct.decaytree) { branch.legs.clear(); }
          direct.amplitude.DECAY_SYM = false;
          const auto expected = gra::spin::ResonanceDecayMatrix(direct, res, "CM");
          CHECK((decay - expected).FrobNorm() < 1e-12);
          const double angle = 0.37;
          direct.pfinal[0].RotateZ(angle);
          for (auto &branch : direct.decaytree) { branch.p4.RotateZ(angle); }
          const auto rotated = gra::spin::ResonanceDecayMatrix(direct, res, res.hel_decay, "CM", false);
          const auto rotation = gra::spin::SpinRotation(gra::spin::SpinHalfRotation(0.0, angle), 2.0);
          CHECK((rotated - rotation.Conj() * decay).FrobNorm() < 1e-11);
        }
        if (model == "MP") { decay = gra::spin::ProductionRotation(lts, lts.process.MP_FRAME, 2.0).Conj() * decay; }
        const auto   rho  = decay.Transpose() * parent.rho * decay.Conj();
        const double norm = rho.Trace().real();
        REQUIRE(norm > 0.0);
        CHECK((pair.rho - rho / norm).FrobNorm() < 1e-10);
        CHECK(pair.integral == Approx(parent.integral * norm).epsilon(1e-10));
        const auto joint_decay = Density::IdentityMatrix(4).Kronecker(decay.Transpose());
        const auto joint_rho = joint_decay * results[4].rho * joint_decay.Dagger();
        CHECK((results[2].rho - joint_rho / joint_rho.Trace().real()).FrobNorm() < 1e-10);
        CHECK(results[2].integral == Approx(results[4].integral * joint_rho.Trace().real()).epsilon(1e-10));
        auto full = gra::spin::ResonanceDecayMatrix(lts, res, "CM");
        if (model == "MP") { full = gra::spin::ProductionRotation(lts, lts.process.MP_FRAME, 2.0).Conj() * full; }
        const double full_norm = (full.Transpose() * parent.rho * full.Conj()).Trace().real();
        CHECK(parent.integral * full_norm == Approx(weight).epsilon(1e-10));
      }
      REQUIRE(accepted == 3);
    }
  }
}

// Real VEGAS workers integrate the physical density once and preserve it through grid reuse
TEST_CASE("Spin metrics follow the multithreaded integration lifecycle", "[qmetrics][integration]") {
  for (const std::string phase : {"F", "C"}) {
    CAPTURE(phase);
    auto       card = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath("gencard/test.json")));
    const auto output                         = gra::aux::ResolveProjectPath("tmp/qmetrics_lifecycle_" + phase);
    card["GENERIC"]["OUTPUT"]                 = output;
    card["GENERIC"]["HIST"]                   = 0;
    card["GENERIC"]["CORES"]                  = 2;
    card["GENERIC"]["NEVENTS"]                = 0;
    card["GENERIC"]["WEIGHTED"]               = true;
    card["GENERIC"]["INTEGRATOR"]             = "VEGAS";
    card["SCATTERING"]["PROCESS"]             = "MP[CON]<" + phase + "> -> pi+ pi- @QMETRICS:true";
    card["GENCUTS"]["<" + phase + ">"]["M"]   = {0.5, 1.0};
    card["GENCUTS"]["<" + phase + ">"]["Rap"] = {-1.0, 1.0};
    card["GENCUTS"]["<" + phase + ">"]["Pt"]  = {0.1, 0.5};
    card["FIDCUTS"]["active"]                 = true;
    card["FIDCUTS"]["CENTRAL"]["*"]["Eta"]    = {-1.0, 1.0};
    card["VETOCUTS"]["active"]                = false;
    card["INTEGRATOR"]["min_samples"]         = 1024;
    card["INTEGRATOR"]["max_samples"]         = 1024;
    card["INTEGRATOR"]["precision"]           = 0.9;
    card["INTEGRATOR"]["VEGAS"]["bins"]       = 8;
    card["INTEGRATOR"]["VEGAS"]["rounds"]     = 1;
    card["INTEGRATOR"]["VEGAS"]["ncall"]      = 256;
    gra::MGraniitti generator;
    generator.ReadInput(card);
    generator.Initialize();
    const auto grid    = nlohmann::json::parse(gra::aux::GetInputData(output + ".vgrid"));
    auto       metrics = generator.proc->state.lts.qmetrics;
    metrics.Deserialize(grid.at("QMETRICS"));
    CheckProtons(metrics.Results());
    const auto result = metrics.Results().at(0);
    Check(result, {1.0, 0.0, 0.0, 0.0, 0.0, 0.0});
    REQUIRE(result.samples == 1024);
    CHECK(grid.at("STAT").at("fidcuts_ok").get<double>() < grid.at("STAT").at("kinematics_ok").get<double>());
    CHECK(result.samples == grid.at("STAT").at("integration_samples").get<std::uint64_t>());
    CHECK(result.integral == Approx(grid.at("STAT").at("sigma").get<double>()).epsilon(1e-11));

    card["GENERIC"]["OUTPUT"] = output + "_restored";
    gra::MGraniitti restored;
    restored.ReadInput(card);
    restored.SetVgridFile(output + ".vgrid");
    restored.Initialize();
    restored.SetNumberOfEvents(8);
    restored.Generate();
    const auto after = nlohmann::json::parse(gra::aux::GetInputData(output + "_restored.vgrid"));
    CHECK(after.at("QMETRICS") == grid.at("QMETRICS"));
  }
}
