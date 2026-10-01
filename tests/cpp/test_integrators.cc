// Unit tests for numerical integration methods
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>
#include <cmath>
#include <complex>
#include <cstdint>
#include <limits>
#include <random>
#include <set>
#include <vector>

#include "Graniitti/MGraniitti.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MPolarCut.h"
#include "Graniitti/Math/MTransport.h"
#include "Graniitti/Sampling/MVEGAS.h"
#include "Graniitti/Analysis/MSpherical.h"

using namespace gra;

using gra::math::pow2;

// Check all polar radial maps against the exact annular area
TEST_CASE("Polar quadrature maps retain the physical radial measure",
          "[integration][polar]") {
  gra::math::PolarParam param;
  param.radial_integrator = "GL";
  param.azimuth_integrator = "Trap";
  param.r_min = 0.2;
  param.r_max = 1.2;
  param.radial_intervals = 32;
  param.azimuth_nodes = 7;

  const auto azimuth = gra::math::PolarAzimuthRule(param);
  double azimuth_sum = 0.0;
  for (const double weight : azimuth.measure) {
    azimuth_sum += weight;
  }
  REQUIRE(azimuth_sum == Approx(2.0 * gra::math::PI).epsilon(1e-14));

  for (const auto map : {gra::math::RadialMap::Linear,
                         gra::math::RadialMap::Log,
                         gra::math::RadialMap::Square}) {
    param.radial_map = map;
    const auto radial = gra::math::PolarRadialRule(param);
    double radial_sum = 0.0;
    for (const double weight : radial.measure) {
      radial_sum += weight;
    }
    const double exact = 0.5 * (pow2(param.r_max) - pow2(param.r_min));
    REQUIRE(radial_sum == Approx(exact).epsilon(1e-13));
  }

  param.radial_integrator = "Boole";
  param.radial_map = gra::math::RadialMap::Linear;
  param.radial_intervals = 8;
  REQUIRE(gra::math::PolarIntervalMultiple(param) == 4);
  REQUIRE(gra::math::PolarNodeCount(param) == 9);
  const auto radial = gra::math::PolarRadialRule(param);
  double radial_sum = 0.0;
  for (const double weight : radial.measure) {
    radial_sum += weight;
  }
  REQUIRE(radial_sum ==
          Approx(0.5 * (pow2(param.r_max) - pow2(param.r_min)))
              .epsilon(1e-14));

  param.radial_map = gra::math::RadialMap::Square;
  REQUIRE_THROWS_AS(gra::math::PolarRadialRule(param), std::invalid_argument);
}

// Check deterministic worker seed cardinality and uniqueness
TEST_CASE("Worker seed generation preserves the fixed seed",
          "[integration][random]") {
  MGraniitti generator;
  constexpr uint32_t fixed_seed = 123456789U;

  CHECK(generator.GenerateUniqueSeeds(fixed_seed, 0).empty());

  const std::vector<uint32_t> one_seed =
      generator.GenerateUniqueSeeds(fixed_seed, 1);
  REQUIRE(one_seed == std::vector<uint32_t>{fixed_seed});

  const std::vector<uint32_t> worker_seeds =
      generator.GenerateUniqueSeeds(fixed_seed, 32);
  REQUIRE(worker_seeds.size() == 32);
  CHECK(worker_seeds.front() == fixed_seed);
  CHECK(std::set<uint32_t>(worker_seeds.begin(), worker_seeds.end()).size() ==
        worker_seeds.size());
  CHECK(worker_seeds == generator.GenerateUniqueSeeds(fixed_seed, 32));
}

// Check generation freezes its envelope and records overflow
TEST_CASE("Generation envelope and overflow state", "[integration]") {
  MAXWEIGHTSTATE maximum_state;
  maximum_state.envelope = 4.0;
  maximum_state.initialized = true;
  CHECK_THROWS_AS(maximum_state.GenerationEnvelope(), std::logic_error);
  maximum_state.ResetGeneration();
  CHECK(maximum_state.generation_envelope == Approx(4.0));
  CHECK(maximum_state.generation_active);
  CHECK(maximum_state.GenerationEnvelope() == Approx(4.0));
  CHECK_THROWS_AS(maximum_state.ResetGeneration(), std::logic_error);
  CHECK(maximum_state.ObserveGenerationOverflow(5.0));
  CHECK_FALSE(maximum_state.ObserveGenerationOverflow(6.0));
  CHECK(maximum_state.envelope == Approx(4.0));
  CHECK(maximum_state.generation_envelope == Approx(4.0));
  CHECK(maximum_state.largest_overflow == Approx(6.0));
  CHECK(maximum_state.generation_overflow);
  maximum_state.envelope = 7.0;
  CHECK_THROWS_AS(maximum_state.GenerationEnvelope(), std::logic_error);
  maximum_state.envelope = 4.0;
  maximum_state.EndGeneration();
  CHECK_FALSE(maximum_state.generation_active);
  CHECK_THROWS_AS(maximum_state.GenerationEnvelope(), std::logic_error);

}

// Check the standalone VEGAS sampler and adapted-grid retention policies
TEST_CASE("Standalone VEGAS owns sampling and grid retention",
          "[integration][VEGAS][adaptation]") {
  VEGASPARAM param;
  param.BINS = 4;
  param.MIN_SUPPORT = 2;
  param.AUTOMATIC_CONVERGENCE = false;

  SECTION("exact proposal sampling") {
    MVEGASIntegrator vegas;
    vegas.SetDimension(1);
    vegas.Init(VEGASStage::Adaptation, 2, 1, param);
    vegas.Data().xmat = {{0.1}, {0.4}, {0.8}, {1.0}};

    const std::vector<double> random_values = {0.5, 0.375};
    std::size_t draw = 0;
    const VEGASSample sample =
        vegas.Sample([&]() { return random_values.at(draw++); });
    REQUIRE(sample.indices == std::vector<std::size_t>{2});
    CHECK(sample.point[0] == Approx(0.25));

    const double grid_inverse_density = 4.0 * 0.3;
    const double mixed_density =
        (1.0 - param.UNIFORM_MIX) / grid_inverse_density + param.UNIFORM_MIX;
    CHECK(std::exp(-sample.log_inverse_density) == Approx(mixed_density));

    CHECK_THROWS_AS(vegas.Sample([]() { return 1.0; }), std::invalid_argument);
  }

  SECTION("last and best grid choices") {
    const auto adapt = [&param](VEGASGridChoice choice) {
      param.BEST_GRID_CHOICE = choice;
      MVEGASIntegrator vegas;
      vegas.SetDimension(1);
      vegas.Init(VEGASStage::Adaptation, 2, 2, param);
      const std::vector<std::vector<double>> initial = vegas.Data().xmat;

      vegas.BeginAdaptationBatch();
      vegas.AccumulateAdaptation(1.0, {1});
      vegas.AccumulateAdaptation(1.0, {2});
      CHECK(vegas.FinishAdaptationBatch(2).status ==
            VEGASAdaptationStatus::Continue);

      vegas.BeginAdaptationBatch();
      vegas.AccumulateAdaptation(1.0, {1});
      vegas.AccumulateAdaptation(9.0, {4});
      CHECK(vegas.FinishAdaptationBatch(2).status ==
            VEGASAdaptationStatus::RoundLimit);
      return std::make_pair(initial, vegas.Data().xmat);
    };

    const auto [best_initial, best_grid] = adapt(VEGASGridChoice::Best);
    const auto [last_initial, last_grid] = adapt(VEGASGridChoice::Last);
    CHECK(best_grid == best_initial);
    CHECK(last_grid != last_initial);
    CHECK(last_grid != best_grid);
  }

  SECTION("frozen proposal serialization") {
    MVEGASIntegrator vegas;
    vegas.SetDimension(1);
    vegas.Init(VEGASStage::Adaptation, 2, 1, param);
    vegas.BeginAdaptationBatch();
    vegas.AccumulateAdaptation(1.0, {1});
    vegas.AccumulateAdaptation(2.0, {4});
    REQUIRE(vegas.FinishAdaptationBatch(2).status ==
            VEGASAdaptationStatus::RoundLimit);

    const json serialized = vegas.Serialize();

    MVEGASIntegrator restored;
    restored.Deserialize(serialized, 1, param);
    CHECK(restored.IsFrozen());
    const auto before = restored.Serialize();
    VEGASPARAM invalid = param;
    invalid.BINS = 3;
    CHECK_THROWS_AS(restored.Init(VEGASStage::Generation, 2, 1, invalid), std::invalid_argument);
    CHECK(restored.Serialize() == before);
    json malformed = serialized;
    malformed["f2mat"] = json::array();
    CHECK_THROWS_AS(restored.Deserialize(malformed, 1, param), std::invalid_argument);
    malformed = serialized;
    malformed["fmat"][0] = json::array();
    CHECK_THROWS_AS(restored.Deserialize(malformed, 1, param), std::invalid_argument);
    invalid = param;
    invalid.ALPHA = -1.0;
    CHECK_THROWS_AS(restored.Deserialize(serialized, 1, invalid), std::invalid_argument);
    CHECK(restored.Serialize() == before);
    CHECK(restored.Dimension() == 1);
    CHECK(restored.Config().BEST_GRID_CHOICE == VEGASGridChoice::Last);

    json incompatible = serialized;
    incompatible["best_grid_choice"] = "best";
    CHECK_THROWS_AS(restored.Deserialize(incompatible, 1, param),
                    std::invalid_argument);
    CHECK(restored.IsFrozen());
    restored.SetDimension(2);
    CHECK_FALSE(restored.IsFrozen());
    CHECK_THROWS_AS(restored.Init(VEGASStage::Generation, 2, 1, param), std::logic_error);
    CHECK_THROWS_AS(restored.Sample([]() { return 0.5; }), std::logic_error);
    restored.Init(VEGASStage::Adaptation, 2, 1, param);
    CHECK(restored.Sample([]() { return 0.5; }).point.size() == 2);
  }

  CHECK(VegasGridChoiceName(ParseVegasGridChoice("last")) == "last");
  CHECK(VegasGridChoiceName(ParseVegasGridChoice("best")) == "best");
  CHECK_THROWS_AS(ParseVegasGridChoice("highest"), std::invalid_argument);
}

// Check the stationary grid when f equals the exact grid and uniform mixture
TEST_CASE("VEGAS preserves an optimal mixed proposal", "[integration][VEGAS][regression]") {
  VEGASPARAM param;
  param.BINS                  = 4;
  param.UNIFORM_MIX           = 0.5;
  param.ROUNDS                = 1;
  param.AUTOMATIC_CONVERGENCE = false;
  MVEGASIntegrator vegas;
  vegas.SetDimension(2);
  constexpr std::size_t calls = 8000;
  vegas.Init(VEGASStage::Adaptation, calls, 1, param);
  const std::vector<std::vector<double>> grid = {{0.1, 0.1}, {0.2, 0.2}, {0.3, 0.3}, {1.0, 1.0}};
  vegas.Data().xmat                           = grid;
  vegas.BeginAdaptationBatch();
  for (const auto &first : indices(grid)) {
    const double dx = grid[first][0] - (first > 0 ? grid[first - 1][0] : 0.0);
    for (const auto &second : indices(grid)) {
      const double dy = grid[second][1] - (second > 0 ? grid[second - 1][1] : 0.0);
      const auto   count =
          std::llround(calls * ((1.0 - param.UNIFORM_MIX) / (param.BINS * param.BINS) + param.UNIFORM_MIX * dx * dy));
      for (long long i = 0; i < count; ++i) { vegas.AccumulateAdaptation(1.0, {first + 1, second + 1}); }
    }
  }
  const auto report = vegas.FinishAdaptationBatch(calls);
  CHECK(report.ess_fraction == Approx(1.0));
  for (const auto &bin : indices(grid)) {
    for (const auto &dimension : indices(grid[bin])) {
      CHECK(vegas.Data().xmat[bin][dimension] == Approx(grid[bin][dimension]).margin(2e-11));
    }
  }
}

// Check normalization and correlated moments after VEGAS adaptation
TEST_CASE("VEGAS integrates separated correlated modes", "[integration][VEGAS][regression]") {
  VEGASPARAM param;
  param.BINS                  = 32;
  param.NCALL                 = 4096;
  param.ROUNDS                = 4;
  param.AUTOMATIC_CONVERGENCE = false;
  MVEGASIntegrator vegas;
  vegas.SetDimension(2);
  vegas.Init(VEGASStage::Adaptation, param.NCALL, param.ROUNDS, param);
  MRandom random;
  random.SetSeed(71);
  const auto target = [](const std::vector<double> &x) {
    return 36.0 * (0.4 * std::pow(x[0] * x[1], 5) + 0.6 * std::pow((1.0 - x[0]) * (1.0 - x[1]), 5));
  };
  while (!vegas.IsFrozen()) {
    const auto calls = vegas.NextAdaptationCalls();
    vegas.BeginAdaptationBatch();
    for (std::size_t i = 0; i < calls; ++i) {
      const auto sample = vegas.Sample([&random]() { return random.U(0.0, 1.0); });
      vegas.AccumulateAdaptation(target(sample.point) * std::exp(sample.log_inverse_density), sample.indices);
    }
    vegas.FinishAdaptationBatch(calls);
  }
  vegas.Init(VEGASStage::Integration, param.NCALL, 1, param);
  gra::statistics::RunningMoments integral, moment;
  constexpr std::size_t           count = 32768;
  for (std::size_t i = 0; i < count; ++i) {
    const auto   sample = vegas.Sample([&random]() { return random.U(0.0, 1.0); });
    const double weight = target(sample.point) * std::exp(sample.log_inverse_density);
    integral.Add(weight);
    moment.Add(sample.point[0] * sample.point[1] * weight);
  }
  for (const auto &[stats, expected] :
       std::vector<std::pair<gra::statistics::RunningMoments, double>>{{integral, 1.0}, {moment, 15.0 / 49.0}}) {
    CHECK(stats.IsFinite());
    CHECK(std::abs(stats.Mean() - expected) <= 5.0L * std::sqrt(stats.M2() / (count * (count - 1))));
  }
}

// Check that changing the bin count preserves the frozen quantile map and density
TEST_CASE("VEGAS resizes frozen grids in both directions", "[integration][VEGAS]") {
  VEGASPARAM param;
  param.BINS = 4;
  param.MIN_SUPPORT = 2;
  param.AUTOMATIC_CONVERGENCE = false;
  MVEGASIntegrator vegas;
  vegas.SetDimension(2);
  vegas.Init(VEGASStage::Adaptation, 2, 1, param);
  vegas.BeginAdaptationBatch();
  vegas.AccumulateAdaptation(1.0, {1, 1});
  vegas.AccumulateAdaptation(1.0, {2, 2});
  REQUIRE(vegas.FinishAdaptationBatch(2).status == VEGASAdaptationStatus::RoundLimit);
  const std::vector<std::vector<double>> original = {{0.1, 0.2}, {0.4, 0.3}, {0.8, 0.9}, {1.0, 1.0}};
  vegas.Data().xmat = original;
  for (const unsigned int bins : {8U, 4U}) {
    param.BINS = bins;
    vegas.Init(VEGASStage::Generation, 2, 1, param);
    const auto edges = vegas.QuantileEdges();
    REQUIRE(edges.size() == 2);
    REQUIRE(vegas.Data().fmat.size() == bins);
    REQUIRE(vegas.Data().f2mat.size() == bins);
    CHECK(vegas.Data().GridCalls() == 0);
    for (const auto &dimension : indices(edges)) {
      REQUIRE(edges[dimension].size() == bins + 1);
      for (std::size_t edge = 1; edge <= bins; ++edge) {
        CHECK(edges[dimension][edge] > edges[dimension][edge - 1]);
      }
      for (const auto &edge : indices(original)) {
        CHECK(edges[dimension][(edge + 1) * bins / 4] == Approx(original[edge][dimension]).margin(1e-10));
      }
    }
    const std::vector<double> draws = {0.5, 0.375, 0.625};
    std::size_t next = 0;
    const auto sample = vegas.Sample([&]() { return draws.at(next++); });
    std::vector<std::size_t> located;
    const double inverse = vegas.Data().GridLogInverseDensity(sample.point, bins, located);
    CHECK(sample.log_inverse_density == Approx(VEGASData::MixedLogInverseDensity(inverse, param.UNIFORM_MIX)));
    CHECK(sample.indices == located);
  }

  // A narrow bin must retain its sampled density when interpolation rounds up
  vegas.Data().xmat = {{0.5, 0.2}, {0.5000000000022204, 0.3}, {0.9, 0.9}, {1.0, 1.0}};
  const std::vector<double> draws = {0.5, 0.4999975, 0.625};
  std::size_t next = 0;
  const auto sample = vegas.Sample([&]() { return draws.at(next++); });
  REQUIRE(sample.indices[0] == 2);
  CHECK(sample.point[0] >= vegas.Data().xmat[0][0]);
  CHECK(sample.point[0] < vegas.Data().xmat[1][0]);
  std::vector<std::size_t> located;
  const double inverse = vegas.Data().GridLogInverseDensity(sample.point, param.BINS, located);
  CHECK(sample.indices == located);
  CHECK(sample.log_inverse_density == Approx(VEGASData::MixedLogInverseDensity(inverse, param.UNIFORM_MIX)));
}

// Check that fast adaptation only controls screening during adaptation
TEST_CASE("Fast adaptation screening controls are independent",
          "[integration][adaptation][screening]") {
  MEventWeightState production;
  CHECK_FALSE(production.adaptation_mode);
  CHECK(production.include_screening);

  MEventWeightState fast;
  fast.ConfigureAdaptation(true);
  CHECK(fast.adaptation_mode);
  CHECK_FALSE(fast.include_screening);

  MEventWeightState screened;
  screened.ConfigureAdaptation(false);
  CHECK(screened.adaptation_mode);
  CHECK(screened.include_screening);
}

// Check exact VEGAS grid and uniform mixture densities
TEST_CASE("VEGAS uniform mixture preserves support",
          "[integration][VEGAS][support]") {
  VEGASPARAM param;
  param.BINS = 4;

  VEGASData data;
  data.SetDimension(1);
  data.Init(VEGASStage::Adaptation, param);
  data.xmat = {{0.1}, {0.4}, {0.8}, {1.0}};

  std::vector<std::size_t> indices;
  const double grid_log_inverse =
      data.GridLogInverseDensity({0.2}, param.BINS, indices);
  REQUIRE(indices == std::vector<std::size_t>{2});
  CHECK(std::exp(grid_log_inverse) == Approx(1.2));

  const double mixed_log_inverse =
      VEGASData::MixedLogInverseDensity(grid_log_inverse, param.UNIFORM_MIX);
  const double grid_density = 1.0 / 1.2;
  const double mixed_density =
      (1.0 - param.UNIFORM_MIX) * grid_density + param.UNIFORM_MIX;
  CHECK(std::exp(-mixed_log_inverse) == Approx(mixed_density));
  CHECK(std::exp(-mixed_log_inverse) >= param.UNIFORM_MIX);

  data.ResetGrid();
  data.AccumulateGrid(1.0, {1});
  data.AccumulateGrid(3.0, {2});
  data.AccumulateGrid(0.0, {3});
  CHECK(data.GridCalls() == 3);
  CHECK(data.GridNonzeroCalls() == 2);
  CHECK(data.GridEffectiveFraction() == Approx(16.0 / 30.0));

  data.ResetGrid();
  CHECK(data.GridCalls() == 0);
  CHECK(data.GridNonzeroCalls() == 0);
  CHECK_THROWS_AS(data.AccumulateGrid(1.0, {0}), std::out_of_range);
  CHECK(data.GridCalls() == 0);
}

// Check rectangular transport kernels and short Sinkhorn runs
TEST_CASE("Optimal transport preserves rectangular matrix orientation",
          "[integration][transport]") {
  gra::MMatrix<double> kernel;
  gra::opt::ConvKernel(2, 3, 0.5, kernel);
  REQUIRE(kernel.size_row() == 2);
  REQUIRE(kernel.size_col() == 3);
  CHECK(kernel(0, 0) == Approx(1.0));
  CHECK(kernel(1, 2) == Approx(1.0));
  CHECK(kernel(0, 2) == Approx(std::exp(-2.0)));

  const gra::MMatrix<double> constant_kernel(2, 3, 1.0);
  const std::vector<double> p = {0.5, 0.5};
  const std::vector<double> q = {0.2, 0.3, 0.5};
  gra::MMatrix<double> coupling;
  const double source_residual =
      gra::opt::SinkHorn(coupling, constant_kernel, p, q, 5);
  CHECK(source_residual == Approx(0.0).margin(1.0e-14));
  REQUIRE(coupling.size_row() == 2);
  REQUIRE(coupling.size_col() == 3);
  CHECK(coupling(0, 0) == Approx(p[0] * q[0]));
  CHECK(coupling(1, 2) == Approx(p[1] * q[2]));

  gra::MMatrix<double> gibbs;
  CHECK_THROWS_AS(
      gra::opt::GibbsKernel(1.0, gra::MMatrix<double>{{0.0, -1.0}}, gibbs),
      std::invalid_argument);
}

// Check the lower bin count and lower-bounded VEGAS width transformation
TEST_CASE("VEGAS grid safeguards preserve positive bin widths",
          "[integration][VEGAS][grid]") {
  SECTION("at least four bins are required") {
    VEGASPARAM param;
    param.BINS = 2;

    VEGASData data;
    data.SetDimension(1);
    CHECK_THROWS_AS(data.Init(VEGASStage::Adaptation, param),
                    std::invalid_argument);
  }

  SECTION("the fixed normalized width floor is enforced") {
    VEGASPARAM param;
    param.BINS = 4;
    constexpr double minimum_width = 2.2204e-12;

    VEGASData data;
    data.SetDimension(1);
    data.Init(VEGASStage::Adaptation, param);
    data.f2mat[0][0] = 1e-300;
    data.f2mat[1][0] = 1e-300;
    data.f2mat[2][0] = 1e300;
    data.f2mat[3][0] = 1e-300;
    data.OptimizeGrid(param);

    double lower = 0.0;
    for (std::size_t bin = 0; bin < param.BINS; ++bin) {
      const double width = data.xmat[bin][0] - lower;
      CHECK(width + 1e-15 >= minimum_width);
      lower = data.xmat[bin][0];
    }
    CHECK(lower == Approx(1.0));
  }
}

// Check dimension-major VEGAS quantile export and malformed grid rejection
TEST_CASE("VEGAS quantile edges validate the adapted grid",
          "[integration][VEGAS][grid]") {
  VEGASData data;
  data.SetDimension(2);
  data.xmat = {{0.1, 0.25}, {0.4, 0.5}, {0.8, 0.75}, {1.0, 1.0}};

  const std::vector<std::vector<double>> edges = data.QuantileEdges(4);
  REQUIRE(edges.size() == 2);
  CHECK(edges[0] == std::vector<double>{0.0, 0.1, 0.4, 0.8, 1.0});
  CHECK(edges[1] == std::vector<double>{0.0, 0.25, 0.5, 0.75, 1.0});

  SECTION("bin count mismatch") {
    CHECK_THROWS_AS(data.QuantileEdges(3), std::invalid_argument);
  }
  SECTION("ragged grid") {
    data.xmat[1].pop_back();
    CHECK_THROWS_AS(data.QuantileEdges(4), std::invalid_argument);
  }
  SECTION("non-finite edge") {
    data.xmat[1][0] = std::numeric_limits<double>::infinity();
    CHECK_THROWS_AS(data.QuantileEdges(4), std::runtime_error);
  }
  SECTION("repeated edge") {
    data.xmat[1][0] = data.xmat[0][0];
    CHECK_THROWS_AS(data.QuantileEdges(4), std::runtime_error);
  }
  SECTION("missing unit boundary") {
    data.xmat.back()[0] = 0.9;
    CHECK_THROWS_AS(data.QuantileEdges(4), std::runtime_error);
  }
  SECTION("edge outside unit interval") {
    data.xmat.back()[0] = 1.1;
    CHECK_THROWS_AS(data.QuantileEdges(4), std::runtime_error);
  }
}

// Check that proposal adaptation leaves common statistics untouched
TEST_CASE("Adaptation samples do not accumulate integration statistics",
          "[integration][statistics][adaptation]") {
  Stats stats;
  MEventWeightState accepted;

  CHECK(stats.ObserveLogSample(accepted, 2.0, std::log(3.0),
                               SamplingStage::Adaptation) == Approx(6.0));
  CHECK(stats.evaluations == Approx(0.0));
  CHECK(stats.integration_samples == 0);
  CHECK(stats.kinematics_ok == Approx(0.0));
  CHECK(stats.fidcuts_ok == Approx(0.0));
  CHECK(stats.vetocuts_ok == Approx(0.0));
  CHECK(stats.amplitude_evaluations == Approx(0.0));
  CHECK(stats.amplitude_ok == Approx(0.0));
  CHECK(stats.amplitude_failures == Approx(0.0));
  CHECK(stats.technical_failures == Approx(0.0));
  CHECK(stats.all_ok == Approx(0.0));
  CHECK(stats.IntegrandWeightCount() == 0);
  CHECK(stats.ProposalDensityCount() == 0);
  CHECK(stats.SamplingWeightCount() == 0);
  CHECK(stats.SamplingWeightPositiveCount() == 0);
}

// Check common raw, proposal-density and importance-weight summaries
TEST_CASE("Frozen proposal samples populate common weight diagnostics",
          "[integration][statistics][weights]") {
  Stats stats;
  MEventWeightState accepted;
  const std::vector<double> raw_weights = {0.0, 2.0, 4.0};
  const std::vector<double> proposal_densities = {1.0, 2.0, 3.0};
  for (std::size_t i = 0; i < raw_weights.size(); ++i) {
    stats.ObserveLogSample(accepted, raw_weights[i],
                           -std::log(proposal_densities[i]),
                           SamplingStage::Integration);
  }

  CHECK(stats.IntegrandWeightCount() == 3);
  CHECK(stats.IntegrandPositiveCount() == 2);
  CHECK(static_cast<double>(stats.IntegrandWeightMean()) == Approx(2.0));
  CHECK(static_cast<double>(stats.IntegrandWeightRelativeStd()) == Approx(1.0));
  CHECK(static_cast<double>(stats.IntegrandWeightMinimumPositive()) ==
        Approx(2.0));
  CHECK(static_cast<double>(stats.IntegrandWeightMaximum()) == Approx(4.0));

  CHECK(stats.ProposalDensityCount() == 3);
  CHECK(stats.ProposalDensityPositiveCount() == 3);
  CHECK(static_cast<double>(stats.ProposalDensityMean()) == Approx(2.0));
  CHECK(static_cast<double>(stats.ProposalDensityRelativeStd()) == Approx(0.5));
  CHECK(static_cast<double>(stats.ProposalDensityMinimumPositive()) ==
        Approx(1.0));
  CHECK(static_cast<double>(stats.ProposalDensityMaximum()) == Approx(3.0));

  CHECK(stats.SamplingWeightCount() == 3);
  CHECK(stats.SamplingWeightPositiveCount() == 2);
  CHECK(static_cast<double>(stats.SamplingWeightMean()) == Approx(7.0 / 9.0));
  CHECK(static_cast<double>(stats.SamplingWeightRelativeStd()) ==
        Approx(std::sqrt(39.0) / 7.0));
  CHECK(static_cast<double>(stats.SamplingWeightMinimumPositive()) ==
        Approx(1.0));
  CHECK(static_cast<double>(stats.SamplingWeightMaximum()) ==
        Approx(4.0 / 3.0));
  CHECK(stats.SamplingEffectiveFraction() == Approx(49.0 / 75.0));
}

// Check logarithmic scaled moments across extreme and nearby weight scales
TEST_CASE("Scaled positive moments preserve logarithmic weight statistics",
          "[integration][statistics][range]") {
  statistics::ScaledPositiveMoments logarithmic;
  logarithmic.AddLogPositive(800.0L);
  logarithmic.AddLogPositive(801.0L);
  CHECK(logarithmic.Count() == 2);
  CHECK(logarithmic.PositiveCount() == 2);
  CHECK(logarithmic.LogScale() == Approx(801.0));
  CHECK(logarithmic.MinLogValue() == Approx(800.0));
  const double expected_logarithmic_relative_std =
      std::sqrt(2.0) * (std::exp(1.0) - 1.0) / (std::exp(1.0) + 1.0);
  CHECK(static_cast<double>(logarithmic.RelativeStandardDeviation()) ==
        Approx(expected_logarithmic_relative_std));

  statistics::ScaledPositiveMoments nearby;
  constexpr double DELTA = 1e-8;
  nearby.Add(1.0);
  nearby.Add(1.0 + DELTA);
  const double expected_nearby_relative_std =
      DELTA / std::sqrt(2.0) / (1.0 + DELTA / 2.0);
  CHECK(static_cast<double>(nearby.RelativeStandardDeviation()) ==
        Approx(expected_nearby_relative_std).epsilon(1e-8));
  CHECK(static_cast<double>(nearby.MinimumPositive()) == Approx(1.0));
  CHECK(static_cast<double>(nearby.Maximum()) == Approx(1.0 + DELTA));

  statistics::ScaledPositiveMoments invalid;
  CHECK_NOTHROW(invalid.Add(-1.0));
  CHECK_NOTHROW(
      invalid.AddLogPositive(std::numeric_limits<long double>::infinity()));
  CHECK_FALSE(invalid.IsFinite());

  statistics::ScaledPositiveMoments malformed;
  CHECK_THROWS_AS(malformed.Restore(2, 1, 10.0L, 0.0L, 0.0L, 0.0L, true),
                  std::invalid_argument);

  const long double below_half = std::nextafter(0.5L, 0.0L);
  const statistics::ExtendedFloatParts encoded =
      statistics::SplitExtendedFloat(below_half);
  CHECK(std::abs(encoded.mantissa) < 1.0);
  CHECK(static_cast<double>(statistics::JoinExtendedFloat(encoded)) ==
        Approx(static_cast<double>(below_half)));
}

// Check the reusable finite-sample and stable accumulation primitives
TEST_CASE("Reusable statistics reproduce analytic sample moments", "[integration][statistics]") {
  const std::vector<std::complex<double>> complex_sample = {{1.0, 1.0}, {3.0, -1.0}};
  const statistics::ComplexStat           complex_stat   = statistics::ComplexMoments(complex_sample);
  CHECK(complex_stat.mean.real() == Approx(2.0));
  CHECK(complex_stat.mean.imag() == Approx(0.0));
  CHECK(complex_stat.second == Approx(6.0));
  CHECK(complex_stat.variance == Approx(2.0));
  CHECK(statistics::UnbiasedComplexVariance(complex_stat.second, complex_stat.mean, complex_sample.size()) ==
        Approx(4.0));

  const std::vector<double> log_sample = {std::log(1.0), std::log(2.0), std::log(4.0)};
  CHECK(statistics::PoweredVariance(log_sample, 1.0) == Approx(2.0 / 7.0));
  CHECK(statistics::WeightedMeanVariance(2.0L, 2.0L, 2.0L) == Approx(1.0));

  statistics::CompensatedSum compensated;
  compensated.Add(1.0e20L);
  compensated.Add(1.0L);
  compensated.Add(-1.0e20L);
  CHECK(static_cast<double>(compensated.Value()) == Approx(1.0));

  const MMatrix<double> jackknife = statistics::JackknifeCovariance({{1.0, 2.0}, {2.0, 4.0}, {3.0, 6.0}});
  CHECK(jackknife[0][0] == Approx(4.0 / 3.0));
  CHECK(jackknife[0][1] == Approx(8.0 / 3.0));
  CHECK(jackknife[1][0] == Approx(8.0 / 3.0));
  CHECK(jackknife[1][1] == Approx(16.0 / 3.0));
  CHECK_THROWS_AS(statistics::JackknifeCovariance({{1.0}, {1.0, 2.0}}), std::invalid_argument);

  statistics::RunningMoments left;
  left.Add(1.0L);
  left.Add(3.0L);
  statistics::RunningMoments right;
  right.Add(5.0L);
  left.Merge(right);
  CHECK(left.Count() == 3);
  CHECK(static_cast<double>(left.Mean()) == Approx(3.0));
  CHECK(static_cast<double>(left.M2()) == Approx(8.0));
  left.Scale(2.0L);
  CHECK(static_cast<double>(left.Mean()) == Approx(6.0));
  CHECK(static_cast<double>(left.M2()) == Approx(32.0));

  statistics::WeightedEstimateMoments estimates;
  estimates.Add(2.0L, 1.0L);
  estimates.Add(4.0L, 3.0L);
  CHECK(static_cast<double>(estimates.Mean()) == Approx(2.5));
  CHECK(static_cast<double>(estimates.VarianceOfMean()) == Approx(0.75));
  CHECK(estimates.ReducedChi2() == Approx(1.0));
}

// Check direct-sampling mean errors and zero-support validity
TEST_CASE("Direct integration statistics use an unbiased sample variance",
          "[integration][statistics]") {
  Stats stats;
  MEventWeightState accepted;
  stats.ObserveLogSample(accepted, 1.0, 0.0, SamplingStage::Integration);
  stats.ObserveLogSample(accepted, 3.0, 0.0, SamplingStage::Integration);
  stats.CalculateCrossSection();

  REQUIRE(stats.ValidCrossSection());
  CHECK(stats.integration_samples == 2);
  CHECK(stats.sigma == Approx(2.0));
  CHECK(stats.sigma_err == Approx(1.0));
  CHECK(stats.RelativeError() == Approx(0.5));

  Stats zero;
  for (std::size_t i = 0; i < 100; ++i) {
    zero.ObserveLogSample(accepted, 0.0, 0.0, SamplingStage::Integration);
  }
  zero.CalculateCrossSection();
  CHECK_FALSE(zero.ValidCrossSection());
  CHECK(std::isinf(zero.RelativeError()));

  Stats nonfinite;
  nonfinite.ObserveLogSample(accepted, std::numeric_limits<double>::infinity(),
                             0.0, SamplingStage::Integration);
  nonfinite.ObserveLogSample(accepted, 1.0, 0.0, SamplingStage::Integration);
  nonfinite.CalculateCrossSection();
  CHECK_FALSE(nonfinite.ValidCrossSection());
  CHECK(std::isinf(nonfinite.RelativeError()));

  stats.trials = 20.0;
  stats.generated = 10;
  stats.N_overflow = 2;
  stats.ResetGeneration();
  CHECK(stats.trials == Approx(0.0));
  CHECK(stats.generated == 0);
  CHECK(stats.N_overflow == 0);
  CHECK(stats.evaluations == Approx(2.0));
  CHECK(stats.sigma == Approx(2.0));
}

// Check sequential cut-flow rates and sampled-envelope diagnostics
TEST_CASE("Integration statistics report conditional physics acceptance",
          "[integration][statistics]") {
  Stats stats;

  MEventWeightState accepted;
  CHECK(stats.ObserveLogSample(accepted, 2.0, std::log(3.0),
                               SamplingStage::Integration) == Approx(6.0));

  MEventWeightState amplitude_rejected;
  amplitude_rejected.amplitude_ok = false;
  stats.ObserveLogSample(amplitude_rejected, 0.0, 0.0,
                         SamplingStage::Integration);

  MEventWeightState veto_rejected;
  veto_rejected.vetocuts_ok = false;
  veto_rejected.amplitude_ok = false;
  stats.ObserveLogSample(veto_rejected, 0.0, 0.0, SamplingStage::Integration);

  MEventWeightState fiducial_rejected;
  fiducial_rejected.fidcuts_ok = false;
  fiducial_rejected.amplitude_ok = false;
  stats.ObserveLogSample(fiducial_rejected, 0.0, 0.0,
                         SamplingStage::Integration);

  MEventWeightState kinematics_rejected;
  kinematics_rejected.kinematics_ok = false;
  kinematics_rejected.amplitude_ok = false;
  stats.ObserveLogSample(kinematics_rejected, 0.0, 0.0,
                         SamplingStage::Integration);

  MEventWeightState technical_failure;
  technical_failure.technical_failure = true;
  technical_failure.amplitude_failure = true;
  technical_failure.amplitude_ok = false;
  technical_failure.forced_accept = true;
  CHECK_FALSE(technical_failure.Valid());
  stats.ObserveLogSample(technical_failure, 0.0, 0.0,
                         SamplingStage::Integration);

  CHECK(stats.evaluations == Approx(6.0));
  CHECK(stats.kinematics_ok == Approx(5.0));
  CHECK(stats.fidcuts_ok == Approx(4.0));
  CHECK(stats.vetocuts_ok == Approx(3.0));
  CHECK(stats.amplitude_evaluations == Approx(3.0));
  CHECK(stats.amplitude_ok == Approx(1.0));
  CHECK(stats.amplitude_failures == Approx(1.0));
  CHECK(stats.technical_failures == Approx(1.0));
  CHECK(stats.all_ok == Approx(1.0));
  // Include all five zero weights in the cross section mean and its uncertainty
  stats.CalculateCrossSection();
  CHECK(stats.integration_samples == 6);
  CHECK(stats.sigma == Approx(1.0));
  CHECK(stats.sigma_err == Approx(1.0));
  CHECK(static_cast<double>(stats.IntegrandWeightMaximum()) == Approx(2.0));
  CHECK(static_cast<double>(stats.SamplingWeightMaximum()) == Approx(6.0));
  CHECK(stats.SamplingWeightCount() == 6);
  CHECK(static_cast<double>(stats.SamplingWeightMean() *
                            stats.SamplingWeightCount()) == Approx(6.0));
  CHECK(stats.SamplingEffectiveFraction() == Approx(1.0 / 6.0));
  CHECK(
      Stats::ConditionalRate(stats.amplitude_ok, stats.amplitude_evaluations) ==
      Approx(1.0 / 3.0));
  CHECK(Stats::ConditionalRate(stats.amplitude_failures,
                               stats.amplitude_evaluations) ==
        Approx(1.0 / 3.0));
  CHECK(Stats::ConditionalRate(stats.technical_failures, stats.evaluations) ==
        Approx(1.0 / 6.0));

  const std::size_t integrand_count = stats.IntegrandWeightCount();
  const std::size_t proposal_count = stats.ProposalDensityCount();
  const std::size_t sampling_count = stats.SamplingWeightCount();
  const long double sampling_mean = stats.SamplingWeightMean();
  const double sampling_fraction = stats.SamplingEffectiveFraction();
  const double evaluations = stats.evaluations;
  const double accepted_events = stats.all_ok;
  stats.ObserveLogSample(accepted, 100.0, 0.0, SamplingStage::Generation);
  CHECK(static_cast<double>(stats.IntegrandWeightMaximum()) == Approx(2.0));
  CHECK(static_cast<double>(stats.SamplingWeightMaximum()) == Approx(6.0));
  stats.ObserveLogSample(accepted, 200.0, 0.0, SamplingStage::Adaptation);
  CHECK(static_cast<double>(stats.IntegrandWeightMaximum()) == Approx(2.0));
  CHECK(static_cast<double>(stats.SamplingWeightMaximum()) == Approx(6.0));
  CHECK(stats.IntegrandWeightCount() == integrand_count);
  CHECK(stats.ProposalDensityCount() == proposal_count);
  CHECK(stats.SamplingWeightCount() == sampling_count);
  CHECK(static_cast<double>(stats.SamplingWeightMean()) ==
        Approx(static_cast<double>(sampling_mean)));
  CHECK(stats.SamplingEffectiveFraction() == Approx(sampling_fraction));
  CHECK(stats.evaluations == Approx(evaluations));
  CHECK(stats.all_ok == Approx(accepted_events));
  stats.CalculateCrossSection();
  CHECK(stats.sigma == Approx(106.0 / 7.0));

  CHECK(stats.MaximumToMean(8.0) == Approx(8.0));
  CHECK(stats.EstimatedUnweightingEfficiency(8.0) == Approx(0.125));
}

// Check logarithmic proposal factors without premature range loss
TEST_CASE("Importance weights combine the integrand and density in log space",
          "[integration][statistics][range]") {
  const statistics::ImportanceSample direct =
      statistics::EvaluateImportanceSample(6.0, -std::log(2.0));
  REQUIRE(direct.HasMaterializedProposalDensity());
  CHECK(direct.proposal_density == Approx(2.0));
  CHECK(direct.weight == Approx(3.0));

  const double large_inverse = statistics::ImportanceWeight(1e-300, 710.0);
  REQUIRE(std::isfinite(large_inverse));
  CHECK(large_inverse ==
        Approx(std::exp(std::log(1e-300) + 710.0)).epsilon(2e-13));
  CHECK(std::isinf(1e-300 * std::exp(710.0)));
  const statistics::ImportanceSample extreme =
      statistics::EvaluateImportanceSample(1e-300, 710.0);
  CHECK_FALSE(extreme.HasMaterializedProposalDensity());
  CHECK(extreme.weight == Approx(large_inverse));

  const double small_inverse = statistics::ImportanceWeight(1e300, -800.0);
  REQUIRE(small_inverse > 0.0);
  CHECK(small_inverse ==
        Approx(std::exp(std::log(1e300) - 800.0)).epsilon(2e-13));
  CHECK(gra::math::IsZero(1e300 * std::exp(-800.0)));
}

// Check stable phase-space moments and logarithmic weighted combinations
TEST_CASE("MCW moments retain tiny scales and extreme logarithmic weights",
          "[integration][statistics][MCW]") {
  gra::kinematics::MCW tiny;
  tiny.Push(1e-200);
  tiny.Push(3e-200);
  const gra::kinematics::MCW restored_scale = tiny * 1e200;
  CHECK(restored_scale.Integral() == Approx(2.0));
  CHECK(restored_scale.IntegralError() == Approx(1.0));

  gra::kinematics::MCW first;
  first.Push(1.0);
  first.Push(3.0);
  gra::kinematics::MCW second;
  second.Push(2.0);
  second.Push(4.0);
  gra::kinematics::MCWSUM combined;
  combined.AddLogWeight(first, 800.0);
  combined.AddLogWeight(second, 799.0);
  const double ratio = std::exp(-1.0);
  CHECK(combined.Integral() == Approx((2.0 + 3.0 * ratio) / (1.0 + ratio)));
  CHECK(combined.IntegralError2() ==
        Approx((1.0 + ratio * ratio) / ((1.0 + ratio) * (1.0 + ratio))));
}

// Check that invalid observations remain invalid through estimates and merging
TEST_CASE("MCW propagates invalid observations", "[integration][statistics][MCW]") {
  for (const double invalid : {std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity()}) {
    gra::kinematics::MCW weight(2.0);
    CHECK_NOTHROW(weight.Push(invalid));
    CHECK(weight.GetN() == Approx(2.0));
    CHECK(std::isnan(weight.Integral()));
    CHECK(std::isnan(weight.IntegralError2()));
    CHECK(std::isnan(weight.IntegralError()));
    CHECK(std::isnan(weight.GetW()));
    CHECK(std::isnan(weight.GetW2()));
    CHECK(std::isnan((weight * 0.5).Integral()));
    CHECK(std::isnan((gra::kinematics::MCW(1.0) + weight).Integral()));
    CHECK(std::isnan(gra::kinematics::MCW(invalid).IntegralError2()));
    gra::kinematics::MCWSUM combined;
    CHECK_THROWS_AS(combined.Add(weight, 1.0), std::invalid_argument);
  }
}

// Check stable direct moments, batch chi2 and scale-normalized ESS
TEST_CASE("Direct integration statistics resist scale and cancellation",
          "[integration][statistics]") {
  MEventWeightState accepted;
  Stats flat;
  constexpr std::size_t CALLS = 1000;
  constexpr double DELTA = 1e-8;
  for (std::size_t i = 0; i < CALLS; ++i) {
    flat.ObserveLogSample(accepted, 1.0 + (i % 2 == 0 ? 0.0 : DELTA), 0.0,
                          SamplingStage::Integration);
  }
  flat.CalculateCrossSection();
  REQUIRE(flat.ValidCrossSection());
  const double expected_error =
      DELTA / (2.0 * std::sqrt(static_cast<double>(CALLS - 1)));
  CHECK(flat.sigma_err / expected_error == Approx(1.0).epsilon(1e-6));

  Stats scaled;
  scaled.ObserveLogSample(accepted, 1e-200, 0.0, SamplingStage::Integration);
  scaled.ObserveLogSample(accepted, 3e-200, 0.0, SamplingStage::Integration);
  scaled.CalculateCrossSection();
  REQUIRE(scaled.ValidCrossSection());
  CHECK(scaled.sigma / 2e-200 == Approx(1.0));
  CHECK(scaled.sigma_err / 1e-200 == Approx(1.0));
  CHECK(scaled.RelativeError() == Approx(0.5));

  Stats batches;
  batches.ResetIntegrationBatch();
  batches.ObserveLogSample(accepted, 1.0, 0.0, SamplingStage::Integration);
  batches.ObserveLogSample(accepted, 3.0, 0.0, SamplingStage::Integration);
  const IntegrationBatchSummary first_batch = batches.UpdateIntegrationChi2();
  REQUIRE(first_batch.valid);
  CHECK(first_batch.samples == 2);
  CHECK(first_batch.mean == Approx(2.0));
  CHECK(first_batch.variance_of_mean == Approx(1.0));
  CHECK(first_batch.ess_fraction == Approx(0.8));
  batches.ResetIntegrationBatch();
  batches.ObserveLogSample(accepted, 1.0 + DELTA, 0.0,
                           SamplingStage::Integration);
  batches.ObserveLogSample(accepted, 3.0 + DELTA, 0.0,
                           SamplingStage::Integration);
  batches.UpdateIntegrationChi2();
  const double expected_chi2 = DELTA * DELTA / 2.0;
  CHECK(batches.chi2 / expected_chi2 == Approx(1.0).epsilon(1e-7));

  WeightedEventStats weighted;
  weighted.Push(2e-200);
  weighted.Push(3e-200);
  CHECK(weighted.EffectiveSampleSize() == Approx(25.0 / 13.0));
  CHECK(weighted.EffectiveSampleFraction() == Approx(25.0 / 26.0));
}

// Check stable statistical state serialization and large-gradient clipping
TEST_CASE("Stable statistics round-trip and clip without overflow",
          "[integration][statistics]") {
  MEventWeightState accepted;
  Stats original;
  original.ResetIntegrationBatch();
  const std::vector<double> raw_weights = {1.0, 3.0, 2.0};
  const std::vector<double> proposal_densities = {1.0, 2.0, 4.0};
  for (std::size_t i = 0; i < raw_weights.size(); ++i) {
    original.ObserveLogSample(accepted, raw_weights[i],
                              -std::log(proposal_densities[i]),
                              SamplingStage::Integration);
  }
  MEventWeightState technical_failure;
  technical_failure.technical_failure = true;
  technical_failure.amplitude_failure = true;
  technical_failure.amplitude_ok = false;
  original.ObserveLogSample(technical_failure, 0.0, -std::log(3.0),
                            SamplingStage::Integration);
  original.ObserveLogSample(accepted, 5.0, 0.0, SamplingStage::Generation);
  original.UpdateIntegrationChi2();
  original.CalculateCrossSection();
  original.integration_runtime = 12.5;

  nlohmann::json payload;
  original.struct2json(payload);
  CHECK(payload.contains("importance_count"));
  CHECK(payload.at("cross_section_count").get<std::size_t>() == 5);
  CHECK_FALSE(payload.contains("sampling_count"));
  CHECK_FALSE(payload.contains("direct_count"));
  CHECK_FALSE(payload.contains("direct_mean"));
  CHECK_FALSE(payload.contains("direct_m2"));
  CHECK_FALSE(payload.contains("direct_finite"));
  CHECK_FALSE(payload.contains("maxW"));
  CHECK_FALSE(payload.contains("maxf"));
  Stats restored;
  restored.json2struct(payload);
  restored.CalculateCrossSection();
  CHECK(restored.sigma == Approx(original.sigma));
  CHECK(restored.sigma_err == Approx(original.sigma_err));
  CHECK(restored.SamplingEffectiveFraction() ==
        Approx(original.SamplingEffectiveFraction()));
  CHECK(restored.SamplingWeightCount() == original.SamplingWeightCount());
  CHECK(restored.SamplingWeightPositiveCount() ==
        original.SamplingWeightPositiveCount());
  CHECK(static_cast<double>(restored.SamplingWeightMean()) ==
        Approx(static_cast<double>(original.SamplingWeightMean())));
  CHECK(static_cast<double>(restored.SamplingWeightRelativeStd()) ==
        Approx(static_cast<double>(original.SamplingWeightRelativeStd())));
  CHECK(static_cast<double>(restored.SamplingWeightMinimumPositive()) ==
        Approx(static_cast<double>(original.SamplingWeightMinimumPositive())));
  CHECK(static_cast<double>(restored.SamplingWeightMaximum()) ==
        Approx(static_cast<double>(original.SamplingWeightMaximum())));
  CHECK(restored.integration_runtime == Approx(original.integration_runtime));
  CHECK(restored.IntegrandWeightCount() == original.IntegrandWeightCount());
  CHECK(restored.IntegrandPositiveCount() == original.IntegrandPositiveCount());
  CHECK(static_cast<double>(restored.IntegrandWeightMean()) ==
        Approx(static_cast<double>(original.IntegrandWeightMean())));
  CHECK(static_cast<double>(restored.IntegrandWeightRelativeStd()) ==
        Approx(static_cast<double>(original.IntegrandWeightRelativeStd())));
  CHECK(static_cast<double>(restored.IntegrandWeightMinimumPositive()) ==
        Approx(static_cast<double>(original.IntegrandWeightMinimumPositive())));
  CHECK(static_cast<double>(restored.IntegrandWeightMaximum()) ==
        Approx(static_cast<double>(original.IntegrandWeightMaximum())));
  CHECK(restored.ProposalDensityCount() == original.ProposalDensityCount());
  CHECK(restored.ProposalDensityPositiveCount() ==
        original.ProposalDensityPositiveCount());
  CHECK(static_cast<double>(restored.ProposalDensityMean()) ==
        Approx(static_cast<double>(original.ProposalDensityMean())));
  CHECK(static_cast<double>(restored.ProposalDensityRelativeStd()) ==
        Approx(static_cast<double>(original.ProposalDensityRelativeStd())));
  CHECK(static_cast<double>(restored.ProposalDensityMinimumPositive()) ==
        Approx(static_cast<double>(original.ProposalDensityMinimumPositive())));
  CHECK(static_cast<double>(restored.ProposalDensityMaximum()) ==
        Approx(static_cast<double>(original.ProposalDensityMaximum())));
  CHECK(restored.amplitude_evaluations ==
        Approx(original.amplitude_evaluations));
  CHECK(restored.amplitude_failures == Approx(original.amplitude_failures));
  CHECK(restored.amplitude_failures == Approx(1.0));
  CHECK(restored.technical_failures == Approx(original.technical_failures));
  CHECK(restored.technical_failures == Approx(1.0));

  Stats continued_original = original;
  Stats continued_restored = restored;
  continued_original.ObserveLogSample(accepted, 5.0, -std::log(6.0),
                                      SamplingStage::Integration);
  continued_restored.ObserveLogSample(accepted, 5.0, -std::log(6.0),
                                      SamplingStage::Integration);
  CHECK(static_cast<double>(continued_restored.IntegrandWeightMean()) ==
        Approx(static_cast<double>(continued_original.IntegrandWeightMean())));
  CHECK(static_cast<double>(continued_restored.ProposalDensityMean()) ==
        Approx(static_cast<double>(continued_original.ProposalDensityMean())));
  CHECK(static_cast<double>(continued_restored.SamplingWeightMean()) ==
        Approx(static_cast<double>(continued_original.SamplingWeightMean())));
  continued_original.CalculateCrossSection();
  continued_restored.CalculateCrossSection();
  CHECK(continued_restored.sigma == Approx(continued_original.sigma));
  CHECK(continued_restored.sigma_err == Approx(continued_original.sigma_err));

  nlohmann::json inconsistent = payload;
  inconsistent["proposal_count"] =
      inconsistent.at("proposal_count").get<std::size_t>() + 1;
  Stats rejected;
  CHECK_THROWS_AS(rejected.json2struct(inconsistent), std::invalid_argument);

  nlohmann::json missing_positive_proposal = payload;
  missing_positive_proposal["proposal_positive_count"] =
      missing_positive_proposal.at("proposal_positive_count")
          .get<std::size_t>() -
      1;
  Stats rejected_positive;
  CHECK_THROWS_AS(rejected_positive.json2struct(missing_positive_proposal),
                  std::invalid_argument);

  nlohmann::json inconsistent_importance = payload;
  inconsistent_importance["importance_count"] =
      inconsistent_importance.at("importance_count").get<std::size_t>() + 1;
  Stats rejected_importance;
  CHECK_THROWS_AS(rejected_importance.json2struct(inconsistent_importance),
                  std::invalid_argument);

  Stats tiny_stats;
  tiny_stats.ObserveLogSample(accepted, 1e-200, 0.0,
                              SamplingStage::Integration);
  tiny_stats.ObserveLogSample(accepted, 3e-200, 0.0,
                              SamplingStage::Integration);
  nlohmann::json tiny_payload;
  tiny_stats.struct2json(tiny_payload);
  Stats tiny_restored;
  tiny_restored.json2struct(tiny_payload);
  tiny_restored.CalculateCrossSection();
  CHECK(tiny_restored.sigma / 2e-200 == Approx(1.0));
  CHECK(tiny_restored.sigma_err / 1e-200 == Approx(1.0));

  statistics::ScaledWeightSums tiny;
  tiny.Add(1e-300);
  tiny.Add(2e-300);
  CHECK(tiny.HasSignificantSignedSum(std::numeric_limits<double>::epsilon()));
  CHECK(static_cast<double>(tiny.EffectiveSampleSize()) == Approx(9.0 / 5.0));

  statistics::ScaledWeightSums cancelled;
  cancelled.Add(1e-300);
  cancelled.Add(-1e-300);
  CHECK_FALSE(cancelled.HasSignificantSignedSum(
      std::numeric_limits<double>::epsilon()));

  std::vector<double> gradient = {std::numeric_limits<double>::max(),
                                  std::numeric_limits<double>::max()};
  statistics::ClipEuclideanNorm(gradient, 10.0);
  REQUIRE(gradient[0] > 0.0);
  REQUIRE(gradient[1] > 0.0);
  CHECK(std::hypot(gradient[0], gradient[1]) == Approx(10.0).epsilon(1e-12));
}

// Check VEGAS technical adaptation defaults without test-only source helpers
TEST_CASE("VEGAS adaptation support controls are explicit",
          "[integration][VEGAS][support]") {
  const VEGASPARAM param;
  CHECK(param.MIN_SUPPORT == 10);
  CHECK(std::size_t{1} < param.MIN_SUPPORT);
  CHECK(std::size_t{10} >= param.MIN_SUPPORT);
  CHECK(param.MAX_SUPPORT_CALLS == 1000000);
}

// Check rolling ordinary ESS convergence criteria
TEST_CASE("VEGAS adaptation convergence requires an ordinary ESS plateau",
          "[integration][VEGAS][convergence]") {
  VEGASPARAM param;
  param.CONVERGENCE_WINDOW = 3;
  param.ESS_REL_TOLERANCE = 0.10;

  CHECK(VegasAdaptationConverged({0.100, 0.105, 0.103}, param));
  CHECK_FALSE(VegasAdaptationConverged({0.100, 0.150, 0.200}, param));
  CHECK_FALSE(VegasAdaptationConverged({0.100, 0.101}, param));
}

// Check that integration retains the adapted grid and resets adaptation data
TEST_CASE("VEGAS integration transition resets adaptation data",
          "[integration][VEGAS]") {
  VEGASPARAM param;
  param.BINS = 4;

  VEGASData data;
  data.SetDimension(1);
  data.Init(VEGASStage::Adaptation, param);
  data.xmat = {{0.1}, {0.4}, {0.8}, {1.0}};
  data.AccumulateGrid(1.0, {2});

  const std::vector<std::vector<double>> adapted = data.xmat;
  data.Init(VEGASStage::Integration, param);

  CHECK(data.xmat == adapted);
  CHECK(gra::math::IsZero(data.GridEffectiveCalls()));
  for (std::size_t bin = 0; bin < param.BINS; ++bin) {
    CHECK(data.fmat[bin][0] == Approx(0.0));
    CHECK(data.f2mat[bin][0] == Approx(0.0));
  }
}

// Check that grid adaptation depends on shape rather than integrand scale
TEST_CASE("VEGAS grid adaptation is scale invariant", "[integration][VEGAS]") {
  VEGASPARAM param;
  param.BINS = 4;
  param.ALPHA = 1.5;

  VEGASData reference;
  reference.SetDimension(1);
  reference.Init(VEGASStage::Adaptation, param);
  const std::vector<double> scores = {0.0, 1.0, 4.0, 0.0};
  for (std::size_t i = 0; i < scores.size(); ++i) {
    reference.f2mat[i][0] = scores[i];
  }
  reference.OptimizeGrid(param);

  VEGASData scaled;
  scaled.SetDimension(1);
  scaled.Init(VEGASStage::Adaptation, param);
  for (std::size_t i = 0; i < scores.size(); ++i) {
    scaled.f2mat[i][0] = scores[i] * 1e-200;
  }
  scaled.OptimizeGrid(param);

  double previous = 0.0;
  for (std::size_t i = 0; i < scores.size(); ++i) {
    CHECK(scaled.xmat[i][0] == Approx(reference.xmat[i][0]));
    CHECK(scaled.xmat[i][0] > previous);
    previous = scaled.xmat[i][0];
  }

  const std::vector<std::vector<double>> adapted = scaled.xmat;
  for (std::size_t i = 0; i < scores.size(); ++i) {
    scaled.f2mat[i][0] = 0.0;
  }
  scaled.OptimizeGrid(param);
  CHECK(scaled.xmat == adapted);

  VEGASData accumulated_reference;
  accumulated_reference.SetDimension(1);
  accumulated_reference.Init(VEGASStage::Adaptation, param);
  accumulated_reference.AccumulateGrid(1.0, {1});
  accumulated_reference.AccumulateGrid(3.0, {2});

  VEGASData accumulated_tiny;
  accumulated_tiny.SetDimension(1);
  accumulated_tiny.Init(VEGASStage::Adaptation, param);
  accumulated_tiny.AccumulateGrid(1e-200, {1});
  accumulated_tiny.AccumulateGrid(3e-200, {2});
  CHECK(gra::math::IsZero(gra::math::pow2(1e-200)));
  CHECK(accumulated_tiny.f2mat[0][0] > 0.0);
  CHECK(accumulated_tiny.GridEffectiveCalls() == Approx(1.6).margin(2e-15));
  CHECK(accumulated_tiny.GridEffectiveCalls() ==
        Approx(accumulated_reference.GridEffectiveCalls()).margin(2e-15));
  for (std::size_t bin = 0; bin < param.BINS; ++bin) {
    CHECK(accumulated_tiny.fmat[bin][0] ==
          Approx(accumulated_reference.fmat[bin][0]).margin(2e-15));
    CHECK(accumulated_tiny.f2mat[bin][0] ==
          Approx(accumulated_reference.f2mat[bin][0]).margin(2e-15));
  }
  accumulated_reference.OptimizeGrid(param);
  accumulated_tiny.OptimizeGrid(param);
  for (std::size_t bin = 0; bin < param.BINS; ++bin) {
    CHECK(accumulated_tiny.xmat[bin][0] ==
          Approx(accumulated_reference.xmat[bin][0]).margin(2e-14));
  }
}

// Check common cross-section continuation for the VEGAS sampling path
TEST_CASE("VEGAS generation extends the common cross-section moments",
          "[integration][VEGAS]") {
  MEventWeightState accepted;
  Stats original;
  for (std::size_t i = 0; i < 1000; ++i) {
    original.ObserveLogSample(accepted, i % 2 == 0 ? 0.0 : 2.0, 0.0,
                              SamplingStage::Integration);
  }
  original.CalculateCrossSection();
  REQUIRE(original.ValidCrossSection());
  const double integration_error = original.sigma_err;

  for (std::size_t i = 0; i < 10; ++i) {
    original.ObserveLogSample(accepted, 1e-12, 0.0, SamplingStage::Generation);
  }
  original.CalculateCrossSection();
  CHECK(original.sigma == Approx((1000.0 + 1e-11) / 1010.0));
  CHECK(original.sigma_err < integration_error);
  CHECK(original.SamplingWeightCount() == 1000);

  nlohmann::json payload;
  original.struct2json(payload);
  CHECK(payload.at("cross_section_count").get<std::size_t>() == 1010);
  Stats restored;
  restored.json2struct(payload);
  restored.CalculateCrossSection();
  CHECK(restored.sigma == Approx(original.sigma));
  CHECK(restored.sigma_err == Approx(original.sigma_err));
}

// Sample a scalar function on a closed uniform grid
std::vector<double> generateFunctionValues(double (*func)(double), double a,
                                           double b, std::size_t N) {
  std::vector<double> values(N + 1);
  double hstep = (b - a) / N;
  for (std::size_t i = 0; i <= N; ++i) {
    values[i] = func(a + i * hstep);
  }
  return values;
}

// Sample a scalar function on a closed two-dimensional uniform grid
template <typename T>
MMatrix<T> generateFunctionValues2D(T (*func)(T, T), T a, T b, T c, T d,
                                    std::size_t M, std::size_t N) {
  MMatrix<T> values(M + 1, N + 1);
  double hstepM = (b - a) / M;
  double hstepN = (d - c) / N;
  for (std::size_t i = 0; i <= M; ++i) {
    for (std::size_t j = 0; j <= N; ++j) {
      values[i][j] = func(a + i * hstepM, c + j * hstepN);
    }
  }
  return values;
}

// Compute the degree-one integration test polynomial
double linearFunction(double x) { return x; }

// Compute the degree-two integration test polynomial
double quadraticFunction(double x) { return x * x; }

// Compute the degree-three integration test polynomial
double cubicFunction(double x) { return x * x * x; }

// Compute the degree-four integration test polynomial
double quarticFunction(double x) { return x * x * x * x; }

// Compute the separable degree-one two-dimensional polynomial
double linearFunction2D(double x, double y) { return x + y; }

// Compute the separable degree-two two-dimensional polynomial
double quadraticFunction2D(double x, double y) { return x * x + y * y; }

// Compute the separable degree-three two-dimensional polynomial
double cubicFunction2D(double x, double y) { return x * x * x + y * y * y; }

// Compute the separable degree-four two-dimensional polynomial
double quarticFunction2D(double x, double y) {
  return x * x * x * x + y * y * y * y;
}

// Tests
TEST_CASE("CSTrapzIntegral", "[integration]") {
  double a = 0.0;
  double b = 1.0;
  std::size_t N = 10;
  double hstep = (b - a) / N;

  SECTION("Linear function") {
    auto f = generateFunctionValues(linearFunction, a, b, N);
    double result = gra::math::CSTrapzIntegral(f, hstep);
    REQUIRE(result == Approx(0.5).margin(1e-15));
  }

  SECTION("Odd number of intervals") {
    const std::size_t odd_intervals = 11;
    const double odd_step = (b - a) / odd_intervals;
    const auto f = generateFunctionValues(linearFunction, a, b, odd_intervals);
    REQUIRE(gra::math::CSTrapzIntegral(f, odd_step) ==
            Approx(0.5).margin(1e-15));
  }

  SECTION("Invalid input") {
    REQUIRE_THROWS_AS(
        gra::math::CSTrapzIntegral(std::vector<double>{1.0}, hstep),
        std::invalid_argument);
    REQUIRE_THROWS_AS(
        gra::math::CSTrapzIntegral(std::vector<double>{0.0, 1.0},
                                   std::numeric_limits<double>::infinity()),
        std::invalid_argument);
  }
}

TEST_CASE("PeriodicTrapzIntegral", "[integration]") {
  const double a = 0.0;
  const double b = 2.0 * gra::math::PI;
  const unsigned int N = 7;
  const double hstep = (b - a) / N;

  SECTION("Constant function with arbitrary node count") {
    std::vector<double> f(N, 1.0);
    const double result = gra::math::PeriodicTrapzIntegral(f, hstep);
    REQUIRE(result == Approx(b - a).epsilon(1e-14));
  }

  SECTION("Midpoint periodic rule integrates first harmonic") {
    const auto [x, w] = gra::math::PeriodicTrapzRule(N, a, b);
    std::vector<double> f(N, 0.0);
    for (std::size_t i = 0; i < f.size(); ++i) {
      REQUIRE(w[i] == Approx(hstep).epsilon(1e-14));
      f[i] = std::cos(x[i]);
    }
    const double result = gra::math::PeriodicTrapzIntegral(f, hstep);
    REQUIRE(result == Approx(0.0).margin(1e-14));
    REQUIRE(x.front() == Approx(gra::math::PI / N).epsilon(1e-14));
    REQUIRE(x.back() < b);
  }

  SECTION("Invalid input") {
    REQUIRE_THROWS_AS(
        gra::math::PeriodicTrapzIntegral(std::vector<double>{}, hstep),
        std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::PeriodicTrapzRule(0, a, b),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::PeriodicTrapzRule(N, b, a),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(
        gra::math::PeriodicTrapzIntegral(std::vector<double>{1.0}, std::numeric_limits<double>::infinity()),
        std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::PeriodicTrapzRule(N, a, std::numeric_limits<double>::infinity()),
                      std::invalid_argument);
  }
}

TEST_CASE("CS13Integral", "[integration]") {
  double a = 0.0;
  double b = 1.0;
  std::size_t N = 10; // Even number of subintervals
  double hstep = (b - a) / N;

  SECTION("Cubic polynomial exactness") {
    auto f = generateFunctionValues(cubicFunction, a, b, N);
    double result = gra::math::CS13Integral(f, hstep);
    REQUIRE(result == Approx(1.0 / 4.0).margin(1e-15));
  }

  SECTION("Invalid number of points") {
    auto f = generateFunctionValues(quadraticFunction, a, b,
                                    N + 1); // Odd number of points
    REQUIRE_THROWS_AS(gra::math::CS13Integral(f, hstep), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::CS13Integral(generateFunctionValues(cubicFunction, a, b, N),
                                              std::numeric_limits<double>::infinity()),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::Simpson13Weight(0), std::invalid_argument);
  }
}

TEST_CASE("CS38Integral", "[integration]") {
  double a = 0.0;
  double b = 1.0;
  std::size_t N = 9; // Multiple of 3 number of subintervals
  double hstep = (b - a) / N;

  SECTION("Cubic function") {
    auto f = generateFunctionValues(cubicFunction, a, b, N);
    double result = gra::math::CS38Integral(f, hstep);
    REQUIRE(result == Approx(0.25).margin(1e-15));
  }

  SECTION("Invalid number of points") {
    auto f = generateFunctionValues(cubicFunction, a, b,
                                    N + 1); // Not a multiple of 3
    REQUIRE_THROWS_AS(gra::math::CS38Integral(f, hstep), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::CS38Integral(generateFunctionValues(cubicFunction, a, b, N),
                                              std::numeric_limits<double>::infinity()),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::Simpson38Weight(0), std::invalid_argument);
  }
}

TEST_CASE("CSBooleIntegral", "[integration]") {
  double a = 0.0;
  double b = 1.0;
  std::size_t N = 8; // Multiple of 4 number of subintervals
  double hstep = (b - a) / N;

  SECTION("Quartic function") {
    auto f = generateFunctionValues(quarticFunction, a, b, N);
    double result = gra::math::CSBooleIntegral(f, hstep);
    REQUIRE(result == Approx(0.2).margin(1e-15));
  }

  SECTION("Invalid number of points") {
    auto f = generateFunctionValues(quarticFunction, a, b,
                                    N + 1); // Not a multiple of 4
    REQUIRE_THROWS_AS(gra::math::CSBooleIntegral(f, hstep), std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::CSBooleIntegral(generateFunctionValues(quarticFunction, a, b, N),
                                                 std::numeric_limits<double>::infinity()),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::BooleWeight(0), std::invalid_argument);
  }
}

TEST_CASE("GaussLegendreRule", "[integration]") {
  SECTION("Constant function with one node") {
    const auto [x, w] = gra::math::GaussLegendreRule(1, 0.0, 2.0);
    REQUIRE(x.size() == 1);
    REQUIRE(w.size() == 1);
    REQUIRE(x[0] == Approx(1.0).epsilon(1e-14));
    REQUIRE(w[0] == Approx(2.0).epsilon(1e-14));
  }

  SECTION("Polynomial exactness on unit interval") {
    for (unsigned int n = 1; n <= 6; ++n) {
      const auto [x, w] = gra::math::GaussLegendreRule(n, 0.0, 1.0);
      for (unsigned int p = 0; p <= 2 * n - 1; ++p) {
        double integral = 0.0;
        for (std::size_t i = 0; i < x.size(); ++i) {
          integral += w[i] * std::pow(x[i], static_cast<double>(p));
        }
        REQUIRE(integral == Approx(1.0 / (p + 1.0)).epsilon(1e-12));
      }
    }
  }

  SECTION("Mapped interval") {
    const double a = 2.0;
    const double b = 5.0;
    const auto [x, w] = gra::math::GaussLegendreRule(3, a, b);
    double integral = 0.0;
    for (std::size_t i = 0; i < x.size(); ++i) {
      integral += w[i] * std::pow(x[i], 3.0);
    }
    REQUIRE(integral ==
            Approx((std::pow(b, 4.0) - std::pow(a, 4.0)) / 4.0).epsilon(1e-12));
  }

  SECTION("Invalid input") {
    REQUIRE_THROWS_AS(gra::math::GaussLegendreRule(0, 0.0, 1.0),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::GaussLegendreRule(2, 1.0, 1.0),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::GaussLegendreRule(2, 2.0, 1.0),
                      std::invalid_argument);
    REQUIRE_THROWS_AS(gra::math::GaussLegendreRule(2, 0.0, std::numeric_limits<double>::infinity()),
                      std::invalid_argument);
  }
}

TEST_CASE("Logarithmic Gauss-Legendre integration includes its Jacobian", "[MMath][integration]") {
  const double minimum = 0.1;
  const double maximum = 3.0;
  const double result =
      gra::math::LogGaussIntegral(12, minimum, maximum, [](const double value) { return value * value; });
  const double exact = (gra::math::pow3(maximum) - gra::math::pow3(minimum)) / 3.0;
  REQUIRE(result == Approx(exact).margin(2.0e-13));
  REQUIRE_THROWS_AS(gra::math::LogGaussIntegral(0, minimum, maximum, [](const double value) { return value; }),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::LogGaussIntegral(8, minimum, std::numeric_limits<double>::infinity(),
                                                [](const double value) { return value; }),
                    std::invalid_argument);
  REQUIRE_THROWS_AS(gra::math::CompositeWeight(0, 0), std::invalid_argument);

  const double log_measure =
      gra::math::LogMeasureGaussIntegral(12, minimum, maximum, [](const double value) { return value * value; });
  REQUIRE(log_measure == Approx((gra::math::pow2(maximum) - gra::math::pow2(minimum)) / 2.0).margin(2.0e-13));
  REQUIRE_THROWS_AS(gra::math::LogMeasureGaussIntegral(8, 0.0, maximum, [](const double value) { return value; }),
                    std::invalid_argument);
}

// ----------------------------------------------------------
// 2D integrators

TEST_CASE("Simpson13Integral2D", "[integration]") {
  double a = 0.0, b = 1.0, c = 0.0, d = 1.0;

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {

      std::size_t M =
          8 + i * 2; // Even number of subintervals in both directions
      std::size_t N = 8 + j * 2;
      double hstepM = (b - a) / M;
      double hstepN = (d - c) / N;

      // Basic
      {
        auto f = generateFunctionValues2D(cubicFunction2D, a, b, c, d, M, N);
        auto W = gra::math::Simpson13Weight2D(M, N);

        double result = gra::math::Simpson13Integral2D(f, W, hstepM, hstepN);
        REQUIRE(result == Approx(1.0 / 2.0).margin(2e-15));
      }

      // Input
      {
        auto f = generateFunctionValues2D(quadraticFunction2D, a, b, c, d,
                                          M + 1, N); // M is odd
        auto W = gra::math::Simpson13Weight2D(M, N);

        REQUIRE_THROWS_AS(gra::math::Simpson13Integral2D(f, W, hstepM, hstepN),
                          std::invalid_argument);
      }
    }
  }
}

TEST_CASE("Simpson38Integral2D", "[integration]") {
  double a = 0.0, b = 1.0, c = 0.0, d = 1.0;

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {

      std::size_t M =
          9 + i * 3; // Multiple of 3 number of subintervals in both directions
      std::size_t N = 9 + j * 3;
      double hstepM = (b - a) / M;
      double hstepN = (d - c) / N;

      // Basic
      {
        auto f = generateFunctionValues2D(cubicFunction2D, a, b, c, d, M, N);
        auto W = gra::math::Simpson38Weight2D(M, N);

        double result = gra::math::Simpson38Integral2D(f, W, hstepM, hstepN);
        REQUIRE(result == Approx(1.0 / 2.0).margin(2e-15));
      }

      // Input
      {
        auto f = generateFunctionValues2D(quadraticFunction2D, a, b, c, d,
                                          M + 1, N); // M is not multiple of 3
        auto W = gra::math::Simpson38Weight2D(M, N);

        REQUIRE_THROWS_AS(gra::math::Simpson38Integral2D(f, W, hstepM, hstepN),
                          std::invalid_argument);
      }
    }
  }
}

TEST_CASE("BooleIntegral2D uses tensor-product weights and quartic exactness",
          "[integration]") {
  double a = 0.0, b = 1.0, c = 0.0, d = 1.0;

  for (std::size_t i = 0; i < 4; ++i) {
    for (std::size_t j = 0; j < 4; ++j) {

      std::size_t M =
          8 + i * 4; // Multiple of 4 number of subintervals in both directions
      std::size_t N = 8 + j * 4;
      double hstepM = (b - a) / M;
      double hstepN = (d - c) / N;

      // Basic
      {
        const std::vector<double> row_weights = gra::math::BooleWeight(M);
        const std::vector<double> column_weights = gra::math::BooleWeight(N);
        auto f = generateFunctionValues2D(quarticFunction2D, a, b, c, d, M, N);
        auto W = gra::math::BooleWeight2D(M, N);

        REQUIRE(row_weights.size() == M + 1);
        REQUIRE(column_weights.size() == N + 1);
        REQUIRE(W.size_row() == M + 1);
        REQUIRE(W.size_col() == N + 1);
        for (std::size_t row = 0; row <= M; ++row) {
          for (std::size_t column = 0; column <= N; ++column) {
            REQUIRE(W[row][column] ==
                    Approx(row_weights[row] * column_weights[column]));
          }
        }

        double result = gra::math::BooleIntegral2D(f, W, hstepM, hstepN);
        REQUIRE(result == Approx(2.0 / 5.0).margin(2e-15));
      }

      // Input
      {
        auto f = generateFunctionValues2D(quadraticFunction2D, a, b, c, d,
                                          M + 1, N); // M is not multiple of 4
        auto W = gra::math::BooleWeight2D(M, N);

        REQUIRE_THROWS_AS(gra::math::BooleIntegral2D(f, W, hstepM, hstepN),
                          std::invalid_argument);
      }
    }
  }
}

// Check that tensor integration rejects nonfinite grid spacing
TEST_CASE("Tensor integration rules require finite steps", "[integration]") {
  const double          infinity = std::numeric_limits<double>::infinity();
  const MMatrix<double> simpson13(3, 3, 1.0);
  REQUIRE_THROWS_AS(gra::math::Simpson13Integral2D(simpson13, gra::math::Simpson13Weight2D(2, 2), infinity, 1.0),
                    std::invalid_argument);

  const MMatrix<double> simpson38(4, 4, 1.0);
  REQUIRE_THROWS_AS(gra::math::Simpson38Integral2D(simpson38, gra::math::Simpson38Weight2D(3, 3), 1.0, infinity),
                    std::invalid_argument);

  const MMatrix<double> boole(5, 5, 1.0);
  REQUIRE_THROWS_AS(gra::math::BooleIntegral2D(boole, gra::math::BooleWeight2D(4, 4), infinity, 1.0),
                    std::invalid_argument);
}

TEST_CASE("Spherical harmonic covariance propagation",
          "[spherical][covariance]") {
  gra::spherical::Omega first;
  first.costheta = 0.25;
  first.phi = 0.4;
  first.fiducial = true;
  first.selected = true;

  gra::spherical::Omega second = first;
  second.costheta = -0.35;
  second.phi = -0.7;
  const std::vector<gra::spherical::Omega> events = {first, second};
  const std::vector<std::size_t> indices = {0, 1};
  const auto estimate =
      gra::spherical::SphericalMoments(events, indices, 1, "det");

  REQUIRE(estimate.value.size() == 4);
  REQUIRE(estimate.covariance.size_row() == 4);
  REQUIRE(estimate.covariance.size_col() == 4);
  REQUIRE(std::abs(estimate.covariance[0][1]) > 0.0);
  for (std::size_t i = 0; i < 4; ++i) {
    REQUIRE(estimate.covariance[i][i] >= 0.0);
  }

  const gra::MMatrix<double> transform{{1.0, 2.0}, {0.0, 1.0}};
  const gra::MMatrix<double> covariance{{4.0, 1.0}, {1.0, 9.0}};
  const auto propagated = gra::spherical::CovarianceProp(transform, covariance);
  REQUIRE(propagated[0][0] == Approx(44.0));
  REQUIRE(propagated[0][1] == Approx(19.0));
  REQUIRE(propagated[1][0] == Approx(19.0));
  REQUIRE(propagated[1][1] == Approx(9.0));
}

TEST_CASE("Spherical response validation and jackknife covariance",
          "[spherical][validation]") {
  const std::vector<gra::spherical::Omega> empty_events;
  const std::vector<std::size_t> empty_indices;
  REQUIRE_THROWS_AS(
      gra::spherical::GetGMixing(empty_events, empty_indices, 1, "fla"),
      std::invalid_argument);
  REQUIRE_THROWS_AS(
      gra::spherical::GetELM(empty_events, empty_indices, 1, "fla"),
      std::invalid_argument);
  REQUIRE(gra::spherical::CalcError(1.0 - 1e-15, 1.0, 100.0) == Approx(0.0));
  REQUIRE_THROWS_AS(gra::spherical::CalcError(0.0, 1.0, 100.0),
                    std::domain_error);

  const std::vector<gra::MMatrix<double>> response = {
      gra::MMatrix<double>{{1.0}}, gra::MMatrix<double>{{2.0}}};
  const auto covariance = gra::spherical::ResponseJackknifeCovariance(
      response, {}, {2.0}, {true}, 0.0);
  REQUIRE(covariance[0][0] == Approx(0.25));
}

TEST_CASE("Spherical weighted moments retain yield covariance",
          "[spherical][weights]") {
  gra::spherical::Omega first;
  first.costheta = 0.2;
  first.phi = 0.3;
  first.weight = 2.0;
  first.fiducial = true;
  first.selected = true;

  gra::spherical::Omega second = first;
  second.costheta = -0.4;
  second.weight = 3.0;
  const std::vector<gra::spherical::Omega> events = {first, second};
  const auto estimate =
      gra::spherical::SphericalMoments(events, {0, 1}, 0, "det");

  REQUIRE(estimate.entries == 2);
  REQUIRE(estimate.sum_weight == Approx(5.0));
  REQUIRE(estimate.sum_weight2 == Approx(13.0));
  REQUIRE(estimate.value[0] == Approx(5.0));
  REQUIRE(estimate.covariance[0][0] == Approx(13.0));
  REQUIRE(gra::spherical::SummarizeWeights(events, {0, 1}, "det")
              .EffectiveEntries() == Approx(25.0 / 13.0));

  std::vector<gra::spherical::Omega> tiny_events = events;
  tiny_events[0].weight = 2e-200;
  tiny_events[1].weight = 3e-200;
  REQUIRE(gra::spherical::SummarizeWeights(tiny_events, {0, 1}, "det")
              .EffectiveEntries() == Approx(25.0 / 13.0));
}

// Check scale-invariant weighted spherical response uncertainties
TEST_CASE("Spherical weighted response errors ignore global weight scale",
          "[spherical][weights]") {
  std::vector<gra::spherical::Omega> reference(4);
  const std::vector<double> costheta = {-0.8, -0.2, 0.3, 0.9};
  const std::vector<double> phi = {-2.0, -0.4, 0.7, 2.4};
  const std::vector<double> weights = {1.0, 2.0, 3.0, 4.0};
  for (std::size_t i = 0; i < reference.size(); ++i) {
    reference[i].costheta = costheta[i];
    reference[i].phi = phi[i];
    reference[i].weight = weights[i];
    reference[i].fiducial = true;
    reference[i].selected = true;
  }
  std::vector<gra::spherical::Omega> scaled = reference;
  for (gra::spherical::Omega &event : scaled) {
    event.weight *= 1e-200;
  }
  const std::vector<std::size_t> indices = {0, 1, 2, 3};
  const auto expected = gra::spherical::GetELM(reference, indices, 1, "fla");
  const auto observed = gra::spherical::GetELM(scaled, indices, 1, "fla");
  REQUIRE(observed.first.size() == expected.first.size());
  REQUIRE(observed.second.size() == expected.second.size());
  for (std::size_t i = 0; i < expected.first.size(); ++i) {
    CHECK(observed.first[i] == Approx(expected.first[i]).epsilon(1e-12));
    CHECK(observed.second[i] == Approx(expected.second[i]).epsilon(1e-11));
  }
}

TEST_CASE("Spherical hypercell boundaries are disjoint",
          "[spherical][binning]") {
  gra::spherical::Omega lower;
  lower.M = 0.0;
  lower.Pt = 0.0;
  lower.Y = 0.0;
  gra::spherical::Omega boundary = lower;
  boundary.M = 1.0;
  gra::spherical::Omega upper = lower;
  upper.M = 2.0;
  const std::vector<gra::spherical::Omega> events = {lower, boundary, upper};

  const auto first =
      gra::spherical::GetIndices(events, {0.0, 1.0}, {-1.0, 1.0}, {-1.0, 1.0});
  const auto second = gra::spherical::GetIndices(
      events, {1.0, 2.0}, {-1.0, 1.0}, {-1.0, 1.0}, {true, false, false});
  REQUIRE(first == std::vector<std::size_t>{0});
  REQUIRE(second == std::vector<std::size_t>{1, 2});
}

TEST_CASE("Spherical active response inversion fixes inactive moments",
          "[spherical][response]") {
  const gra::MMatrix<double> response{{2.0, 10.0}, {0.0, 1.0}};
  const auto inverse =
      gra::spherical::ActiveResponsePseudoInverse(response, {true, false}, 0.0);
  const std::vector<double> unfolded = inverse * std::vector<double>{4.0, 7.0};
  REQUIRE(unfolded[0] == Approx(2.0));
  REQUIRE(unfolded[1] == Approx(0.0));

  const gra::MMatrix<double> covariance{{4.0, 1.0, 0.0, 0.0},
                                        {1.0, 9.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0},
                                        {0.0, 0.0, 0.0, 0.0}};
  REQUIRE(gra::spherical::HarmDotProdError({1.0, 2.0, 0.0, 0.0}, covariance,
                                           {true, true, false, false},
                                           1) == Approx(std::sqrt(44.0)));
  REQUIRE_THROWS_AS(
      gra::spherical::ResponseJackknifeCovariance(
          {gra::MMatrix<double>{{1.0}}, gra::MMatrix<double>{{1.0}}}, {},
          {1.0, 2.0}, {true, true}, 0.0),
      std::invalid_argument);
}

// Check circular cuts against their exact area, including disconnected radial intervals
TEST_CASE("Polar circular cuts preserve area and transverse moments", "[integration][polar][cuts]") {
  const int turns = GENERATE(-8, -2, 0, 2, 8);
  gra::math::PolarParam param;
  param.radial_integrator  = "GL";
  param.azimuth_integrator = "GL";
  param.radial_map         = gra::math::RadialMap::Log;
  param.r_min              = 0.2;
  param.r_max              = 2.0;
  param.radial_intervals   = 16;
  param.azimuth_nodes      = 24;
  const gra::math::PolarCutRule rule(param);
  for (const double angle : {0.0, 0.137, 1.83}) {
    CAPTURE(turns, angle);
    const std::array<double, 2> centre = {std::cos(angle), std::sin(angle)};
    const auto                  nodes  = rule.Nodes({centre}, 0.3, angle + 2.0 * gra::math::PI * turns);
    double                      area   = 0.0;
    double                      x      = 0.0;
    double                      y      = 0.0;
    for (const auto& node : nodes) {
      REQUIRE(std::hypot(node.x - centre[0], node.y - centre[1]) >= 0.3 * (1.0 - 1e-12));
      area += node.weight;
      x += node.weight * node.x;
      y += node.weight * node.y;
    }
    REQUIRE(area == Approx(gra::math::PI * (4.0 - 0.04 - 0.09)).epsilon(1e-5));
    REQUIRE(x == Approx(-gra::math::PI * 0.09 * centre[0]).epsilon(1e-4).margin(1e-8));
    REQUIRE(y == Approx(-gra::math::PI * 0.09 * centre[1]).epsilon(1e-4).margin(1e-8));
  }
}

// Integrate a continuous piecewise quadratic field across an independent matching circle
TEST_CASE("Polar matching circles resolve piecewise fields for every radial rule", "[integration][polar][matching]") {
  const int turns = GENERATE(-2, 0, 2);
  for (const std::string integrator : {"GL", "1/3", "3/8", "Boole"}) {
    for (const auto map : {gra::math::RadialMap::Linear, gra::math::RadialMap::Log, gra::math::RadialMap::Square}) {
      if (map == gra::math::RadialMap::Square && integrator != "GL") { continue; }
      gra::math::PolarParam param;
      param.radial_integrator  = integrator;
      param.azimuth_integrator = "GL";
      param.radial_map         = map;
      param.r_min             = 0.2;
      param.r_max             = 2.0;
      param.radial_intervals  = integrator == "GL" ? 12 : 48;
      param.azimuth_nodes     = 12;
      const gra::math::PolarCutRule rule(param);
      const auto nodes = rule.Nodes({{-1.0, 0.0}}, 0.3, 0.137 + 2.0 * gra::math::PI * turns, {{1.0, 0.0, 0.4}});
      double integral = 0.0;
      double moment   = 0.0;
      for (const auto& node : nodes) {
        const double field = std::max(0.0, 0.16 - gra::math::pow2(node.x - 1.0) - node.y * node.y);
        integral += node.weight * field;
        moment += node.weight * node.x * field;
      }
      CAPTURE(turns, integrator, static_cast<int>(map));
      REQUIRE(integral == Approx(gra::math::PI * 0.16 * 0.16 / 2.0).epsilon(5e-4));
      REQUIRE(moment == Approx(integral).epsilon(5e-4));
    }
  }
}
