// Spline flow neural importance sampler & integrator tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <catch.hpp>

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <random>
#include <string>
#include <thread>
#include <utility>
#include <vector>

#include "Graniitti/Sampling/MNeuroJac.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"

using gra::aux::indices;

namespace {

// Write one compact standalone NUMERICS_NEUROJAC card for optimizer tests
std::string WriteNeuroJacCard(
    const std::string &loss, const std::string &permutation = "random",
    double bias_decay = 0.0, const std::string &activation = "silu",
    const std::string &lr_schedule = "none",
    std::size_t validation_patience = 5, double validation_min_delta = 1e-3,
    std::size_t rounds = 1, double alpha_initial = 1.0,
    double alpha_final = 2.0, double annealing_fraction = 0.7,
    bool vegas_init = true, std::size_t flow_components = 1,
    bool vegas_spline_init = true) {
  const std::filesystem::path filename =
      "tmp/graniitti_neurojac_" + loss + "_" + permutation + "_" + activation +
      "_" + lr_schedule +
      (bias_decay > 0.0 ? "_bias_decay.json" : "_weight_decay.json");
  nlohmann::json block = {
      {"flow_components", flow_components},
      {"layers", 4},
      {"bins", 6},
      {"hidden", std::vector<std::size_t>{8}},
      {"activation", activation},
      {"permutation", permutation},
      {"buffer_size", 128},
      {"batch_size", 32},
      {"val_size", 32},
      {"val_patience", validation_patience},
      {"val_min_delta", validation_min_delta},
      {"replay_frac", 0.5},
      {"replay_cap", 256},
      {"rounds", rounds},
      {"epochs", 2},
      {"vegas_init", vegas_init},
      {"vegas_spline_init", vegas_init && vegas_spline_init},
      {"vegas_ncall", 128},
      {"vegas_rounds", 2},
      {"vegas_sampler_rounds", vegas_init ? 1 : 0},
      {"loss",
       {{"name", loss},
        {"alpha", {alpha_initial, alpha_final}},
        {"anneal_frac", annealing_fraction}}},
      {"learning_rate", 5e-3},
      {"final_learning_rate", 5e-4},
      {"lr_schedule", lr_schedule},
      {"uniform_mix", 0.05},
      {"min_bin", 1e-3},
      {"min_derivative", 1e-3},
      {"density_tol", 1e-6},
      {"gradient_clip", 10.0},
      {"adamw_beta1", 0.9},
      {"adamw_beta2", 0.999},
      {"adamw_epsilon", 1e-8},
      {"adamw_weight_decay", 1e-4},
      {"adamw_decay_biases", bias_decay}};
  nlohmann::json document = {{"NUMERICS_NEUROJAC", block}};

  std::ofstream output(filename);
  if (!output.good()) {
    throw std::runtime_error("WriteNeuroJacCard: cannot create test card");
  }
  output << document.dump(2);
  return filename.string();
}

// Compute a positive non-factorizing target with a known unit-square integral
double TrainingTarget(const std::vector<double> &point) {
  return 0.25 + point[0] + 2.0 * point[1] + point[0] * point[1];
}

// Copy one row-major proposal point for scalar reference checks
std::vector<double> BatchPoint(const gra::neurojac::MFlowSampleBatch &samples,
                               std::size_t event) {
  return std::vector<double>(
      samples.points.begin() +
          static_cast<std::ptrdiff_t>(event * samples.dimension),
      samples.points.begin() +
          static_cast<std::ptrdiff_t>((event + 1) * samples.dimension));
}

// Convert proposal samples into target events with exact collection densities
std::vector<gra::neurojac::MFlowTrainingEvent>
TrainingEvents(const gra::neurojac::MFlowSampleBatch &samples) {
  std::vector<gra::neurojac::MFlowTrainingEvent> events;
  events.reserve(samples.Size());
  for (std::size_t event = 0; event < samples.Size(); ++event) {
    std::vector<double> point = BatchPoint(samples, event);
    const double target = TrainingTarget(point);
    events.push_back({std::move(point), target, samples.log_densities[event]});
  }
  return events;
}

// Evaluate one reference mini-batch objective without automatic differentiation
double
FlowObjective(const gra::neurojac::MSplineFlow &flow,
              const std::vector<gra::neurojac::MFlowTrainingEvent> &events,
              const std::vector<double> &coefficients,
              gra::neurojac::FlowLoss loss, double uniform_mix) {
  const double coefficient_sum = coefficients[0] + coefficients[1];
  double objective = 0.0;
  for (std::size_t i = 0; i < events.size(); ++i) {
    const double flow_density = std::exp(flow.LogDensity(events[i].point));
    const double mixed_density =
        (1.0 - uniform_mix) * flow_density + uniform_mix;
    const double coefficient = coefficients[i] / coefficient_sum;
    objective += loss == gra::neurojac::FlowLoss::ForwardKL
                     ? -coefficient * std::log(mixed_density)
                     : coefficient / mixed_density;
  }
  return objective;
}

// Evaluate one normalized forward KL objective for a spline mixture
double
MixtureForwardKL(const gra::neurojac::MSplineFlowMixture &mixture,
                 const std::vector<gra::neurojac::MFlowTrainingEvent> &events,
                 const std::vector<double> &coefficients, double uniform_mix) {
  double objective = 0.0;
  for (std::size_t event = 0; event < events.size(); ++event) {
    const double mixed_density =
        (1.0 - uniform_mix) *
            std::exp(mixture.LogDensity(events[event].point)) +
        uniform_mix;
    objective -= coefficients[event] * std::log(mixed_density);
  }
  return objective;
}

// Evaluate one mixed flow and uniform proposal log density
double MixedLogDensityReference(const gra::neurojac::MSplineFlow &flow,
                                const std::vector<double> &point,
                                double uniform_mix) {
  const double flow_term = std::log1p(-uniform_mix) + flow.LogDensity(point);
  const double uniform_term = std::log(uniform_mix);
  const double maximum = std::max(flow_term, uniform_term);
  return maximum + std::log(std::exp(flow_term - maximum) +
                            std::exp(uniform_term - maximum));
}

// Build globally normalized empirical Renyi escort coefficients
std::vector<double> RenyiEscortCoefficients(
    const gra::neurojac::MSplineFlow &flow,
    const std::vector<gra::neurojac::MFlowTrainingEvent> &events, double alpha,
    double uniform_mix) {
  std::vector<double> log_coefficients(events.size());
  double maximum = -std::numeric_limits<double>::infinity();
  for (std::size_t event = 0; event < events.size(); ++event) {
    log_coefficients[event] =
        alpha * std::log(std::abs(events[event].target)) -
        events[event].collection_log_density +
        (1.0 - alpha) *
            MixedLogDensityReference(flow, events[event].point, uniform_mix);
    maximum = std::max(maximum, log_coefficients[event]);
  }
  std::vector<double> coefficients(events.size());
  double sum = 0.0;
  for (std::size_t event = 0; event < events.size(); ++event) {
    coefficients[event] = std::exp(log_coefficients[event] - maximum);
    sum += coefficients[event];
  }
  for (double &coefficient : coefficients) {
    coefficient /= sum;
  }
  return coefficients;
}

// Evaluate the empirical Renyi divergence up to constant target terms
double EmpiricalRenyiObjective(
    const gra::neurojac::MSplineFlow &flow,
    const std::vector<gra::neurojac::MFlowTrainingEvent> &events, double alpha,
    double uniform_mix) {
  std::vector<double> log_terms(events.size());
  double maximum = -std::numeric_limits<double>::infinity();
  for (std::size_t event = 0; event < events.size(); ++event) {
    log_terms[event] =
        alpha * std::log(std::abs(events[event].target)) -
        events[event].collection_log_density +
        (1.0 - alpha) *
            MixedLogDensityReference(flow, events[event].point, uniform_mix);
    maximum = std::max(maximum, log_terms[event]);
  }
  double scaled_sum = 0.0;
  for (double log_term : log_terms) {
    scaled_sum += std::exp(log_term - maximum);
  }
  return (maximum + std::log(scaled_sum)) / (alpha - 1.0);
}

} // namespace

TEST_CASE("MSpline1D has an analytic Jacobian-consistent inverse",
          "[neurojac][spline]") {
  const std::vector<double> widths = {0.15, 0.35, 0.50};
  const std::vector<double> heights = {0.25, 0.20, 0.55};
  const std::vector<double> derivatives = {0.8, 1.3, 0.6, 1.1};

  for (double x : std::vector<double>{0.0, 0.03, 0.15, 0.41, 0.82, 1.0}) {
    const gra::neurojac::MSplineResult forward =
        gra::neurojac::MSpline1D::Forward(x, widths, heights, derivatives);
    const gra::neurojac::MSplineResult inverse =
        gra::neurojac::MSpline1D::Inverse(forward.value, widths, heights,
                                          derivatives);
    CHECK(inverse.value == Approx(x).margin(2e-11));
    CHECK(inverse.log_jacobian + forward.log_jacobian ==
          Approx(0.0).margin(2e-10));
  }
}

// Check inverse branch selection under a common local-bin rescaling
TEST_CASE("MSpline1D inverse tolerances follow the coefficient scale",
          "[neurojac][spline]") {
  const std::vector<double> derivatives = {0.5, 3.0};
  const double coordinate = 0.37;
  const gra::neurojac::MSplineResult reference =
      gra::neurojac::MSpline1D::Forward(coordinate, {1.0}, {1.0}, derivatives);

  constexpr double scale = 1e-200;
  const std::vector<double> widths = {scale, 1.0 - scale};
  const std::vector<double> heights = {scale, 1.0 - scale};
  const std::vector<double> scaled_derivatives = {0.5, 3.0, 1.0};
  const gra::neurojac::MSplineResult forward =
      gra::neurojac::MSpline1D::Forward(scale * coordinate, widths, heights,
                                        scaled_derivatives);
  const gra::neurojac::MSplineResult inverse =
      gra::neurojac::MSpline1D::Inverse(forward.value, widths, heights,
                                        scaled_derivatives);

  CHECK(forward.value / scale == Approx(reference.value).epsilon(1e-12));
  CHECK(forward.log_jacobian == Approx(reference.log_jacobian).margin(1e-12));
  CHECK(inverse.value / scale == Approx(coordinate).epsilon(1e-12));
  CHECK(inverse.log_jacobian + forward.log_jacobian ==
        Approx(0.0).margin(1e-12));
}

// Check the inverse root and density when its quadratic coefficient b is negative
TEST_CASE("Steep splines retain their inverse and density", "[neurojac][spline][regression]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 1;
  config.bins   = 2;
  gra::neurojac::MSplineFlow flow;
  flow.Configure(1, config, 19);
  for (double slope : {1e6, 1e8}) {
    CAPTURE(slope);
    const std::vector<double> derivatives = {1.0, slope, 1.0};
    const auto                inverse = gra::neurojac::MSpline1D::Inverse(0.25, {0.5, 0.5}, {0.5, 0.5}, derivatives);
    const long double         s       = slope;
    const long double         t       = 2.0L / (s + 1.0L + std::sqrt((s - 1.0L) * (s - 1.0L) + 4.0L));
    CHECK(inverse.value == Approx(static_cast<double>((1.0L - t) / 2.0L)).margin(2e-15));
    const auto forward = gra::neurojac::MSpline1D::Forward(inverse.value, {0.5, 0.5}, {0.5, 0.5}, derivatives);
    CHECK(forward.value == Approx(0.25).margin(2e-8));
    CHECK(forward.log_jacobian + inverse.log_jacobian == Approx(0.0).margin(2e-8));

    auto parameters = flow.Parameters();
    parameters[5]   = slope - config.min_derivative;
    flow.SetParameters(parameters);
    CHECK(flow.LogDensity({0.25}) == Approx(inverse.log_jacobian).margin(2e-8));
    CHECK(flow.LogDensityBatch({0.25}).front() == Approx(inverse.log_jacobian).margin(2e-8));
    const auto          objective = gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::ForwardKL, 1.0);
    std::vector<double> gradient;
    CHECK(std::isfinite(flow.LossGradient({{{0.25}, 1.0, 0.0}}, {1.0}, {0}, 1, *objective, 0.05, gradient)));
    CHECK(gra::AllFinite(gradient));
    const double step = slope * 1e-3;
    parameters[5] += step;
    flow.SetParameters(parameters);
    const double upper = -std::log(0.95 * std::exp(flow.LogDensity({0.25})) + 0.05);
    parameters[5] -= 2.0 * step;
    flow.SetParameters(parameters);
    const double lower = -std::log(0.95 * std::exp(flow.LogDensity({0.25})) + 0.05);
    CHECK(gradient[5] == Approx((upper - lower) / (2.0 * step)).epsilon(2e-3).margin(1e-18));
  }
}

TEST_CASE("MSplineFlow starts as an exact unit-cube identity",
          "[neurojac][flow]") {
  for (gra::neurojac::FlowActivation activation :
       std::vector<gra::neurojac::FlowActivation>{
           gra::neurojac::FlowActivation::Relu,
           gra::neurojac::FlowActivation::Silu,
           gra::neurojac::FlowActivation::Tanh}) {
    for (gra::neurojac::FlowPermutation permutation :
         std::vector<gra::neurojac::FlowPermutation>{
             gra::neurojac::FlowPermutation::None,
             gra::neurojac::FlowPermutation::Interleave,
             gra::neurojac::FlowPermutation::Random}) {
      CAPTURE(gra::neurojac::FlowActivationName(activation));
      CAPTURE(gra::neurojac::FlowPermutationName(permutation));
      gra::neurojac::MNeuroJacConfig config;
      config.layers = 4;
      config.bins = 6;
      config.hidden = {8};
      config.activation = activation;
      config.permutation = permutation;

      gra::neurojac::MSplineFlow flow;
      flow.Configure(4, config, 17);
      const std::vector<double> base = {0.13, 0.37, 0.57, 0.91};
      const gra::neurojac::MFlowSample sample = flow.Forward(base);

      REQUIRE(sample.point.size() == base.size());
      for (std::size_t i = 0; i < base.size(); ++i) {
        CHECK(sample.point[i] == Approx(base[i]).margin(2e-12));
      }
      CHECK(sample.log_density == Approx(0.0).margin(2e-11));
      CHECK(flow.LogDensity(sample.point) ==
            Approx(sample.log_density).margin(2e-11));

      const std::vector<std::size_t> identity = {0, 1, 2, 3};
      std::vector<std::size_t> transformed_counts(base.size(), 0);
      bool has_layer_specific_permutation = false;
      for (std::size_t layer = 0; layer < flow.Couplings().size(); ++layer) {
        const auto &coupling = flow.Couplings()[layer];
        if (permutation == gra::neurojac::FlowPermutation::None) {
          CHECK(coupling.permutation_ == identity);
        } else if (layer > 0 && coupling.permutation_ !=
                                    flow.Couplings()[layer - 1].permutation_) {
          has_layer_specific_permutation = true;
        }
        for (std::size_t coordinate : coupling.transform_indices_) {
          ++transformed_counts[coupling.permutation_[coordinate]];
        }
      }
      const auto [minimum_count, maximum_count] = std::minmax_element(
          transformed_counts.begin(), transformed_counts.end());
      CHECK(*maximum_count - *minimum_count <= 1);
      if (permutation != gra::neurojac::FlowPermutation::None) {
        CHECK(has_layer_specific_permutation);
      }

      gra::neurojac::MSplineFlow repeated;
      repeated.Configure(4, config, 17);
      for (std::size_t layer = 0; layer < flow.Couplings().size(); ++layer) {
        CHECK(repeated.Couplings()[layer].permutation_ ==
              flow.Couplings()[layer].permutation_);
      }
    }
  }
}

// Check direct spline initialization from one frozen VEGAS inverse CDF
TEST_CASE("NEUROJAC initializes spline parameters from a frozen VEGAS grid",
          "[neurojac][flow][vegas-init]") {
  const std::vector<double> uniform_edges = {0.0, 0.2, 0.4, 0.6, 0.8, 1.0};
  const std::vector<std::vector<double>> identity_grid = {uniform_edges,
                                                          uniform_edges};
  const std::vector<std::vector<double>> skewed_grid = {
      {0.0, 0.02, 0.08, 0.24, 0.55, 1.0}, uniform_edges};
  const std::vector<double> spline_knots = {0.25, 0.5, 0.75};
  const std::vector<double> interpolated_edges = {0.035, 0.16, 0.4725};

  for (gra::neurojac::FlowPermutation permutation :
       {gra::neurojac::FlowPermutation::None,
        gra::neurojac::FlowPermutation::Interleave,
        gra::neurojac::FlowPermutation::Random}) {
    CAPTURE(permutation);
    gra::neurojac::MNeuroJacConfig config;
    config.layers = 4;
    config.bins = 4;
    config.hidden = {5};
    config.permutation = permutation;
    REQUIRE(skewed_grid.front().size() - 1 != config.bins);

    gra::neurojac::MSplineFlow identity;
    identity.Configure(2, config, 137);
    identity.InitializeVegasGrid(identity_grid);
    for (const std::vector<double> &base : std::vector<std::vector<double>>{
             {0.03, 0.17}, {0.51, 0.49}, {0.93, 0.81}}) {
      const gra::neurojac::MFlowSample sample = identity.Forward(base);
      REQUIRE(sample.point.size() == base.size());
      CHECK(sample.point[0] == Approx(base[0]).margin(2e-12));
      CHECK(sample.point[1] == Approx(base[1]).margin(2e-12));
      CHECK(sample.log_density == Approx(0.0).margin(2e-11));
    }

    gra::neurojac::MSplineFlow adapted;
    adapted.Configure(2, config, 137);
    adapted.InitializeVegasGrid(skewed_grid);
    for (std::size_t knot = 0; knot < spline_knots.size(); ++knot) {
      const std::vector<double> base = {spline_knots[knot], spline_knots[knot]};
      const gra::neurojac::MFlowSample sample = adapted.Forward(base);
      CHECK(sample.point[0] == Approx(interpolated_edges[knot]).margin(2e-11));
      CHECK(sample.point[1] == Approx(base[1]).margin(2e-11));
    }
    CHECK(adapted.LogDensity({0.01, 0.5}) > adapted.LogDensity({0.9, 0.5}));

    for (const std::vector<double> &base : std::vector<std::vector<double>>{
             {0.03, 0.17}, {0.51, 0.49}, {0.93, 0.81}}) {
      const gra::neurojac::MFlowSample sample = adapted.Forward(base);
      CHECK(std::isfinite(sample.log_density));
      CHECK(adapted.LogDensity(sample.point) ==
            Approx(sample.log_density).margin(2e-11));
    }
  }

  gra::neurojac::MNeuroJacConfig uncovered_config;
  uncovered_config.layers = 1;
  uncovered_config.bins = 4;
  uncovered_config.hidden = {5};
  uncovered_config.permutation = gra::neurojac::FlowPermutation::Interleave;
  gra::neurojac::MSplineFlow uncovered;
  uncovered.Configure(2, uncovered_config, 137);
  const std::vector<double> original_parameters = uncovered.Parameters();
  CHECK_THROWS_AS(uncovered.InitializeVegasGrid(skewed_grid), std::logic_error);
  CHECK(uncovered.Parameters() == original_parameters);
}

// Check that VEGAS steering cannot be silently ignored by direct training
TEST_CASE("NEUROJAC requires its VEGAS grid before the first training round",
          "[neurojac][flow][vegas-init]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("forward_kl"));
  integrator.Configure(2, 149);
  REQUIRE(integrator.Config().vegas_init);

  gra::MRandom random;
  random.SetSeed(151);
  integrator.BufferEvents(TrainingEvents(
      integrator.SampleBatch(random, integrator.FreshEventCount())));
  const std::vector<double> original_parameters =
      integrator.Flow().Parameters();
  CHECK_THROWS_AS(integrator.TrainRound(), std::logic_error);
  CHECK(integrator.Flow().Parameters() == original_parameters);
}

TEST_CASE("MSplineFlow batched kernels match scalar nonidentity evaluation",
          "[neurojac][batch]") {
  gra::neurojac::MNeuroJacConfig config;
  config.batch_size = 17;
  config.layers = 4;
  config.bins = 6;
  config.hidden = {8, 8};
  gra::neurojac::MSplineFlow flow;
  flow.Configure(4, config, 91);

  std::vector<double> parameters = flow.Parameters();
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    parameters[i] += 2e-3 * std::sin(static_cast<double>(i + 1));
  }
  flow.SetParameters(std::move(parameters));

  constexpr std::size_t dimension = 4;
  constexpr std::size_t events = 96;
  std::vector<double> bases(events * dimension);
  for (std::size_t event = 0; event < events; ++event) {
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      bases[event * dimension + coordinate] =
          std::fmod(0.137 * static_cast<double>(event + 1) +
                        0.223 * static_cast<double>(coordinate + 1),
                    1.0);
    }
  }

  const gra::neurojac::MFlowSampleBatch batched_forward =
      flow.ForwardBatch(bases);
  const std::vector<double> batched_density = flow.LogDensityBatch(bases);
  REQUIRE(batched_forward.Size() == events);
  REQUIRE(batched_density.size() == events);
  for (std::size_t event = 0; event < events; ++event) {
    const std::vector<double> base(
        bases.begin() + static_cast<std::ptrdiff_t>(event * dimension),
        bases.begin() + static_cast<std::ptrdiff_t>((event + 1) * dimension));
    const gra::neurojac::MFlowSample scalar_forward = flow.Forward(base);
    CHECK(batched_density[event] ==
          Approx(flow.LogDensity(base)).margin(3e-11));
    CHECK(batched_forward.log_densities[event] ==
          Approx(scalar_forward.log_density).margin(3e-11));
    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      CHECK(batched_forward.points[event * dimension + coordinate] ==
            Approx(scalar_forward.point[coordinate]).margin(3e-11));
    }
  }

  std::vector<double> checksums(4, 0.0);
  std::vector<std::thread> workers;
  for (std::size_t worker = 0; worker < checksums.size(); ++worker) {
    workers.emplace_back([&, worker] {
      for (std::size_t repeat = 0; repeat < 4; ++repeat) {
        const std::vector<double> densities = flow.LogDensityBatch(bases);
        checksums[worker] +=
            std::accumulate(densities.begin(), densities.end(), 0.0);
      }
    });
  }
  for (std::thread &worker : workers) {
    worker.join();
  }
  for (std::size_t worker = 1; worker < checksums.size(); ++worker) {
    CHECK(checksums[worker] == Approx(checksums.front()).margin(1e-12));
  }

  gra::neurojac::MSplineFlow one_dimensional;
  one_dimensional.Configure(1, config, 92);
  std::vector<double> one_dimensional_parameters = one_dimensional.Parameters();
  for (std::size_t i = 0; i < one_dimensional_parameters.size(); ++i) {
    one_dimensional_parameters[i] +=
        1e-3 * std::cos(static_cast<double>(i + 1));
  }
  one_dimensional.SetParameters(std::move(one_dimensional_parameters));
  const std::vector<double> one_dimensional_points = {0.03, 0.17, 0.51, 0.88};
  const gra::neurojac::MFlowSampleBatch one_dimensional_forward =
      one_dimensional.ForwardBatch(one_dimensional_points);
  const std::vector<double> one_dimensional_density =
      one_dimensional.LogDensityBatch(one_dimensional_points);
  for (std::size_t event = 0; event < one_dimensional_points.size(); ++event) {
    CHECK(one_dimensional_forward.points[event] ==
          Approx(
              one_dimensional.Forward({one_dimensional_points[event]}).point[0])
              .margin(3e-11));
    CHECK(one_dimensional_density[event] ==
          Approx(one_dimensional.LogDensity({one_dimensional_points[event]}))
              .margin(3e-11));
  }
}

TEST_CASE("Spline mixtures have exact densities and reverse weight gradients",
          "[neurojac][mixture][gradient]") {
  gra::neurojac::MNeuroJacConfig config;
  config.flow_components = 2;
  config.layers = 3;
  config.bins = 5;
  config.hidden = {7};
  config.permutation = gra::neurojac::FlowPermutation::Interleave;

  gra::neurojac::MSplineFlowMixture mixture;
  mixture.Configure(2, config, 2027);
  const std::size_t component_parameters =
      mixture.Component(0).Parameters().size();
  std::vector<double> parameters = mixture.Parameters();
  REQUIRE(parameters.size() == 2 * component_parameters + 1);
  for (std::size_t parameter = 0; parameter < component_parameters;
       ++parameter) {
    parameters[parameter] +=
        2e-3 * std::sin(static_cast<double>(parameter + 1));
    parameters[component_parameters + parameter] -=
        3e-3 * std::cos(static_cast<double>(parameter + 1));
  }
  parameters.back() = 0.6;
  mixture.SetParameters(parameters);

  const std::vector<double> weights = mixture.Weights();
  REQUIRE(weights.size() == 2);
  CHECK(weights[0] > 0.0);
  CHECK(weights[1] > 0.0);
  CHECK(weights[0] + weights[1] == Approx(1.0).margin(1e-15));
  for (const std::vector<double> &point :
       std::vector<std::vector<double>>{{0.13, 0.27}, {0.61, 0.83}}) {
    const double reference =
        std::log(weights[0] * std::exp(mixture.Component(0).LogDensity(point)) +
                 weights[1] * std::exp(mixture.Component(1).LogDensity(point)));
    CHECK(mixture.LogDensity(point) == Approx(reference).margin(2e-12));
  }

  const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
      {{0.17, 0.31}, 0.8, 0.0},
      {{0.43, 0.79}, 1.2, 0.0},
      {{0.88, 0.22}, 0.6, 0.0}};
  const std::vector<double> coefficients = {0.2, 0.3, 0.5};
  const std::vector<std::size_t> indices = {0, 1, 2};
  const auto objective = gra::neurojac::CreateFlowObjective(
      gra::neurojac::FlowLoss::ForwardKL, 1.0);
  std::vector<double> gradient;
  const double loss =
      mixture.LossGradient(events, coefficients, indices, indices.size(),
                           *objective, 0.05, gradient);
  CHECK(loss == Approx(MixtureForwardKL(mixture, events, coefficients, 0.05))
                    .margin(2e-12));
  REQUIRE(gradient.size() == parameters.size());

  std::vector<std::size_t> checks = {parameters.size() - 1};
  for (std::size_t component = 0; component < 2; ++component) {
    const auto begin = gradient.begin() + static_cast<std::ptrdiff_t>(
                                              component * component_parameters);
    const auto largest = std::max_element(
        begin, begin + static_cast<std::ptrdiff_t>(component_parameters),
        [](double left, double right) {
          return std::abs(left) < std::abs(right);
        });
    checks.push_back(static_cast<std::size_t>(largest - gradient.begin()));
  }
  constexpr double step = 1e-6;
  for (std::size_t parameter : checks) {
    std::vector<double> shifted = parameters;
    shifted[parameter] += step;
    mixture.SetParameters(shifted);
    const double upper = MixtureForwardKL(mixture, events, coefficients, 0.05);
    shifted[parameter] -= 2.0 * step;
    mixture.SetParameters(shifted);
    const double lower = MixtureForwardKL(mixture, events, coefficients, 0.05);
    const double numerical = (upper - lower) / (2.0 * step);
    CHECK(gradient[parameter] == Approx(numerical).epsilon(5e-5).margin(5e-7));
  }
  mixture.SetParameters(parameters);
}

// Check actual draws and analytic integrals for three distinct unequal-weight flows
TEST_CASE("Spline superposition samples its full probability density", "[neurojac][mixture][regression]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(
      WriteNeuroJacCard("renyi", "none", 0.0, "silu", "none", 5, 1e-3, 3, 1.0, 2.0, 0.7, false, 3, false));
  integrator.Configure(1, 59);
  integrator.Freeze();
  auto                      model       = nlohmann::json::parse(integrator.SerializeModel());
  const std::vector<double> weights     = {0.2, 0.3, 0.5};
  const std::vector<double> midpoints   = {0.1, 0.5, 0.9};
  const double              mix         = integrator.Config().uniform_mix;
  double                    probability = mix * 0.4;
  std::vector<double>       parameters;
  for (const auto &component : indices(weights)) {
    gra::neurojac::MSplineFlow flow;
    flow.Configure(1, integrator.Config(), 59);
    flow.InitializeVegasGrid({{0.0, midpoints[component], 1.0}});
    parameters.insert(parameters.end(), flow.Parameters().begin(), flow.Parameters().end());
    double lower = 0.0;
    double upper = 1.0;
    for (std::size_t step = 0; step < 60; ++step) {
      const double base = 0.5 * (lower + upper);
      if (flow.Forward({base}).point[0] < 0.4) {
        lower = base;
      } else {
        upper = base;
      }
    }
    probability += (1.0 - mix) * weights[component] * 0.5 * (lower + upper);
  }
  parameters.push_back(std::log(weights[0] / weights[2]));
  parameters.push_back(std::log(weights[1] / weights[2]));
  model["ACTIVE_COMPONENTS"] = weights.size();
  model["PARAMETERS"]        = parameters;
  model["PARAMETER_COUNT"]   = parameters.size();
  integrator.DeserializeModel(model.dump(), 1);
  for (bool batched : {false, true}) {
    CAPTURE(batched);
    gra::MRandom random;
    random.SetSeed(batched ? 61 : 67);
    constexpr std::size_t           count   = 8192;
    const auto                      samples = integrator.SampleBatch(random, batched ? count : 0);
    gra::statistics::RunningMoments volume, moment, fraction;
    double                          density_error = 0.0;
    for (std::size_t i = 0; i < count; ++i) {
      const auto   sample = batched ? gra::neurojac::MFlowSample{{samples.points[i]}, samples.log_densities[i]}
                                    : integrator.Sample(random);
      const double x      = sample.point[0];
      const double weight = std::exp(-sample.log_density);
      volume.Add(weight);
      moment.Add(x * x * weight);
      fraction.Add(x < 0.4 ? 1.0 : 0.0);
      density_error = std::max(density_error, std::abs(sample.log_density - integrator.LogDensity(sample.point)));
    }
    CHECK(density_error < 2e-9);
    for (const auto &[stats, expected] : std::vector<std::pair<gra::statistics::RunningMoments, double>>{
             {volume, 1.0}, {moment, 1.0 / 3.0}, {fraction, probability}}) {
      CHECK(stats.IsFinite());
      CHECK(std::abs(stats.Mean() - expected) <= 5.0L * std::sqrt(stats.M2() / (count * (count - 1))));
    }
  }
}

TEST_CASE("VEGAS presampling and spline initialization are independent",
          "[neurojac][flow][vegas-init]") {
  const std::vector<std::vector<double>> skewed_grid = {
      {0.0, 0.04, 0.16, 0.48, 1.0}, {0.0, 0.1, 0.3, 0.65, 1.0}};

  gra::neurojac::MFlowIntegrator uninitialized;
  uninitialized.ReadParameters(
      WriteNeuroJacCard("forward_kl", "interleave", 0.0, "silu", "none", 5,
                        1e-3, 1, 1.0, 2.0, 0.7, true, 1, false));
  uninitialized.Configure(2, 191);
  const std::vector<double> original = uninitialized.Mixture().Parameters();
  uninitialized.InitializeVegasGrid(skewed_grid);
  CHECK_FALSE(uninitialized.Config().vegas_spline_init);
  CHECK(uninitialized.Mixture().Parameters() == original);

  gra::neurojac::MFlowIntegrator initialized;
  initialized.ReadParameters(WriteNeuroJacCard("forward_kl", "interleave", 0.0,
                                               "silu", "none", 5, 1e-3, 1, 1.0,
                                               2.0, 0.7, true, 1, true));
  initialized.Configure(2, 191);
  initialized.InitializeVegasGrid(skewed_grid);
  CHECK(initialized.Config().vegas_spline_init);
  CHECK(initialized.Mixture().Parameters() != original);
}

TEST_CASE("NEUROJAC validates optimizer and VEGAS configuration",
          "[neurojac][config][validation]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("forward_kl", "none", 1e-4));
  integrator.Configure(2, 11);
  CHECK(integrator.Flow().Dimension() == 2);

  const std::vector<unsigned char> &decay_mask =
      integrator.Flow().WeightDecayMask();
  CHECK(std::count(decay_mask.begin(), decay_mask.end(), 0) > 0);
  CHECK(std::count(decay_mask.begin(), decay_mask.end(), 1) > 0);

  gra::neurojac::MFlowIntegrator all_decay;
  all_decay.ReadParameters(WriteNeuroJacCard("forward_kl", "none", 1e-4));
  all_decay.Configure(2, 11);
  CHECK(all_decay.Config().adamw_decay_biases == Approx(1e-4));
  CHECK(all_decay.Config().permutation == gra::neurojac::FlowPermutation::None);

  gra::neurojac::MNeuroJacConfig invalid_replay;
  for (const double minimum : {1.0, 2.0}) {
    gra::neurojac::MNeuroJacConfig invalid_slope;
    invalid_slope.min_derivative = minimum;
    CHECK_THROWS_AS(invalid_slope.Validate(2), std::invalid_argument);
  }
  invalid_replay.replay_fraction = 1.0;
  CHECK_THROWS_AS(invalid_replay.Validate(2), std::invalid_argument);
  invalid_replay.replay_fraction = 0.5;
  invalid_replay.replay_capacity = 0;
  CHECK_THROWS_AS(invalid_replay.Validate(2), std::invalid_argument);

  gra::neurojac::MNeuroJacConfig invalid_vegas;
  invalid_vegas.vegas_ncall = 0;
  CHECK_THROWS_AS(invalid_vegas.Validate(2), std::invalid_argument);
  invalid_vegas.vegas_ncall = 128;
  invalid_vegas.vegas_rounds = 0;
  CHECK_THROWS_AS(invalid_vegas.Validate(2), std::invalid_argument);
  invalid_vegas.vegas_rounds = 2;
  invalid_vegas.rounds = 2;
  invalid_vegas.vegas_sampler_rounds = 3;
  CHECK_THROWS_AS(invalid_vegas.Validate(2), std::invalid_argument);
  invalid_vegas.vegas_init = false;
  invalid_vegas.vegas_sampler_rounds = 1;
  CHECK_THROWS_AS(invalid_vegas.Validate(2), std::invalid_argument);

  gra::neurojac::MNeuroJacConfig invalid_mixture;
  invalid_mixture.flow_components = 0;
  gra::neurojac::MSplineFlowMixture mixture;
  CHECK_THROWS_AS(mixture.Configure(2, invalid_mixture, 17), std::invalid_argument);

  gra::neurojac::MFlowLossConfig invalid_loss;
  invalid_loss.alpha_initial = 0.9;
  CHECK_THROWS_AS(invalid_loss.Validate(), std::invalid_argument);
  invalid_loss.alpha_initial = 1.0;
  invalid_loss.alpha_final = 0.9;
  CHECK_THROWS_AS(invalid_loss.Validate(), std::invalid_argument);
  invalid_loss.alpha_final = 2.0;
  invalid_loss.annealing_fraction = 0.0;
  CHECK_THROWS_AS(invalid_loss.Validate(), std::invalid_argument);
}

TEST_CASE("Spline-flow reverse gradients agree with finite differences",
          "[neurojac][gradient]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 2;
  config.bins = 4;
  config.hidden = {5};
  gra::neurojac::MSplineFlow flow;
  flow.Configure(2, config, 29);
  std::vector<double> nonidentity_parameters = flow.Parameters();
  for (std::size_t i = 0; i < nonidentity_parameters.size(); ++i) {
    nonidentity_parameters[i] += 2e-3 * std::sin(static_cast<double>(i + 1));
  }
  flow.SetParameters(nonidentity_parameters);

  const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
      {{0.23, 0.71}, 1.0, 0.0}, {{0.82, 0.34}, 1.0, 0.0}};
  const std::vector<double> coefficients = {0.3, 0.7};
  const std::vector<std::size_t> indices = {0, 1};
  constexpr double uniform_mix = 0.05;

  for (gra::neurojac::FlowLoss loss : std::vector<gra::neurojac::FlowLoss>{
           gra::neurojac::FlowLoss::ForwardKL, gra::neurojac::FlowLoss::Chi2}) {
    const std::unique_ptr<gra::neurojac::MFlowObjective> loss_function =
        gra::neurojac::CreateFlowObjective(loss, 1.0);
    std::vector<double> gradient;
    const double objective =
        flow.LossGradient(events, coefficients, indices, events.size(),
                          *loss_function, uniform_mix, gradient);
    CHECK(objective ==
          Approx(FlowObjective(flow, events, coefficients, loss, uniform_mix))
              .margin(2e-12));

    double stochastic_objective = 0.0;
    std::vector<double> stochastic_gradient(gradient.size(), 0.0);
    for (std::size_t event = 0; event < events.size(); ++event) {
      std::vector<double> batch_gradient;
      const double batch_objective = flow.LossGradient(
          events, coefficients, std::vector<std::size_t>{event}, events.size(),
          *loss_function, uniform_mix, batch_gradient);
      stochastic_objective +=
          batch_objective / static_cast<double>(events.size());
      for (std::size_t parameter = 0; parameter < gradient.size();
           ++parameter) {
        stochastic_gradient[parameter] +=
            batch_gradient[parameter] / static_cast<double>(events.size());
      }
    }
    CHECK(stochastic_objective == Approx(objective).margin(2e-12));
    for (std::size_t parameter = 0; parameter < gradient.size(); ++parameter) {
      CHECK(stochastic_gradient[parameter] ==
            Approx(gradient[parameter]).margin(2e-11));
    }

    const auto largest = std::max_element(
        gradient.begin(), gradient.end(), [](double left, double right) {
          return std::abs(left) < std::abs(right);
        });
    REQUIRE(largest != gradient.end());
    REQUIRE(std::abs(*largest) > 1e-8);

    const std::vector<double> reference = flow.Parameters();
    const std::size_t checks = std::min<std::size_t>(8, gradient.size());
    for (std::size_t check = 0; check < checks; ++check) {
      const std::size_t parameter =
          check * (gradient.size() - 1) / std::max<std::size_t>(checks - 1, 1);
      std::vector<double> shifted = reference;
      constexpr double step = 1e-6;
      shifted[parameter] += step;
      flow.SetParameters(shifted);
      const double plus =
          FlowObjective(flow, events, coefficients, loss, uniform_mix);
      shifted[parameter] -= 2.0 * step;
      flow.SetParameters(shifted);
      const double minus =
          FlowObjective(flow, events, coefficients, loss, uniform_mix);

      const double numerical = (plus - minus) / (2.0 * step);
      CHECK(gradient[parameter] ==
            Approx(numerical).epsilon(3e-5).margin(3e-7));
    }
    flow.SetParameters(reference);
  }
}

// Check inverse derivatives at a vanishing quadratic or linear coefficient
TEST_CASE("Spline gradients remain correct at degenerate inverse coefficients", "[neurojac][gradient]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 1;
  config.bins = 2;
  const auto objective = gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::ForwardKL, 1.0);
  for (const std::vector<double> &values : std::vector<std::vector<double>>{
           {2.0, 3.0, 1.0 / 6.0}, {1.0, 3.0, 0.25}, {1.0, 1.0, 0.23}}) {
    gra::neurojac::MSplineFlow flow;
    flow.Configure(1, config, 17);
    auto parameters = flow.Parameters();
    parameters[4] = std::log(std::expm1(values[0] - config.min_derivative));
    parameters[5] = std::log(std::expm1(values[1] - config.min_derivative));
    flow.SetParameters(parameters);
    const std::vector<gra::neurojac::MFlowTrainingEvent> events = {{{values[2]}, 1.0, 0.0}};
    std::vector<double> gradient;
    flow.LossGradient(events, {1.0}, {0}, 1, *objective, 0.05, gradient);
    REQUIRE(gradient.size() == parameters.size());
    for (const auto &parameter : indices(parameters)) {
      auto shifted = parameters;
      shifted[parameter] += 1e-6;
      flow.SetParameters(shifted);
      const double plus = -std::log(0.05 + 0.95 * std::exp(flow.LogDensity(events[0].point)));
      shifted[parameter] -= 2e-6;
      flow.SetParameters(shifted);
      const double minus = -std::log(0.05 + 0.95 * std::exp(flow.LogDensity(events[0].point)));
      CAPTURE(values, parameter);
      CHECK(std::isfinite(gradient[parameter]));
      CHECK(gradient[parameter] == Approx((plus - minus) / 2e-6).epsilon(3e-5).margin(3e-7));
    }
  }
}

// Check every knot slope against scalar densities across all spline bins
TEST_CASE("Spline slope gradients remain local to the selected bin", "[neurojac][gradient]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 1;
  config.bins = 6;
  gra::neurojac::MSplineFlow flow;
  flow.Configure(1, config, 43);
  auto parameters = flow.Parameters();
  for (std::size_t knot = 0; knot <= config.bins; ++knot) { parameters[2 * config.bins + knot] += 0.1 * knot; }
  const auto objective = gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::ForwardKL, 1.0);
  for (std::size_t bin = 0; bin < config.bins; ++bin) {
    const std::vector<double> point = {(bin + 0.37) / config.bins};
    flow.SetParameters(parameters);
    std::vector<double> gradient;
    flow.LossGradient({{point, 1.0, 0.0}}, {1.0}, {0}, 1, *objective, config.uniform_mix, gradient);
    for (std::size_t knot = 0; knot <= config.bins; ++knot) {
      const std::size_t parameter = 2 * config.bins + knot;
      CAPTURE(bin, knot);
      if (knot != bin && knot != bin + 1) { CHECK(gradient[parameter] == Approx(0.0).epsilon(0.0).margin(1e-14)); }
      auto shifted = parameters;
      shifted[parameter] += 1e-6;
      flow.SetParameters(shifted);
      const double plus = -std::log(config.uniform_mix + (1.0 - config.uniform_mix) * std::exp(flow.LogDensity(point)));
      shifted[parameter] -= 2e-6;
      flow.SetParameters(shifted);
      const double minus = -std::log(config.uniform_mix + (1.0 - config.uniform_mix) * std::exp(flow.LogDensity(point)));
      CHECK(gradient[parameter] == Approx((plus - minus) / 2e-6).epsilon(3e-5).margin(3e-7));
    }
  }
}

TEST_CASE("Renyi order follows the configured cosine schedule",
          "[neurojac][loss][schedule]") {
  gra::neurojac::MFlowLossConfig config;
  config.kind = gra::neurojac::FlowLoss::Renyi;
  config.alpha_initial = 1.0;
  config.alpha_final = 2.0;
  config.annealing_fraction = 0.5;

  CHECK(gra::neurojac::ScheduledRenyiAlpha(config, 0, 9) == Approx(1.0));
  CHECK(gra::neurojac::ScheduledRenyiAlpha(config, 2, 9) == Approx(1.5));
  CHECK(gra::neurojac::ScheduledRenyiAlpha(config, 4, 9) == Approx(2.0));
  CHECK(gra::neurojac::ScheduledRenyiAlpha(config, 8, 9) == Approx(2.0));
  CHECK(gra::neurojac::ScheduledRenyiAlpha(config, 0, 1) == Approx(2.0));
  CHECK_THROWS_AS(gra::neurojac::ScheduledRenyiAlpha(config, 9, 9),
                  std::invalid_argument);
  CHECK(gra::neurojac::FlowLossName(gra::neurojac::FlowLoss::Renyi) == "renyi");
}

TEST_CASE("Renyi escort gradient matches the empirical divergence",
          "[neurojac][loss][gradient]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 2;
  config.bins = 4;
  config.hidden = {5};
  gra::neurojac::MSplineFlow flow;
  flow.Configure(2, config, 103);
  std::vector<double> parameters = flow.Parameters();
  for (std::size_t parameter = 0; parameter < parameters.size(); ++parameter) {
    parameters[parameter] +=
        3e-3 * std::cos(static_cast<double>(parameter + 1));
  }
  flow.SetParameters(parameters);

  const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
      {{0.17, 0.73}, 0.8, std::log(1.2)},
      {{0.48, 0.29}, 2.1, std::log(0.7)},
      {{0.86, 0.61}, 1.3, std::log(1.5)}};
  const std::vector<std::size_t> indices = {0, 1, 2};
  constexpr double alpha = 1.6;
  constexpr double uniform_mix = 0.05;
  const std::vector<double> coefficients =
      RenyiEscortCoefficients(flow, events, alpha, uniform_mix);
  const std::unique_ptr<gra::neurojac::MFlowObjective> loss_function =
      gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::Renyi, alpha);
  REQUIRE(loss_function->UsesProposalDensity());

  std::vector<double> gradient;
  const double surrogate =
      flow.LossGradient(events, coefficients, indices, events.size(),
                        *loss_function, uniform_mix, gradient);
  CHECK(std::isfinite(surrogate));

  const std::vector<double> reference = flow.Parameters();
  const std::size_t checks = std::min<std::size_t>(10, gradient.size());
  for (std::size_t check = 0; check < checks; ++check) {
    const std::size_t parameter =
        check * (gradient.size() - 1) / std::max<std::size_t>(checks - 1, 1);
    std::vector<double> shifted = reference;
    constexpr double step = 1e-6;
    shifted[parameter] += step;
    flow.SetParameters(shifted);
    const double plus =
        EmpiricalRenyiObjective(flow, events, alpha, uniform_mix);
    shifted[parameter] -= 2.0 * step;
    flow.SetParameters(shifted);
    const double minus =
        EmpiricalRenyiObjective(flow, events, alpha, uniform_mix);

    const double numerical = (plus - minus) / (2.0 * step);
    CHECK(gradient[parameter] == Approx(numerical).epsilon(4e-5).margin(4e-7));
  }
  flow.SetParameters(reference);
}

TEST_CASE("All conditioner activations have finite reverse gradients",
          "[neurojac][gradient][activation]") {
  for (gra::neurojac::FlowActivation activation :
       std::vector<gra::neurojac::FlowActivation>{
           gra::neurojac::FlowActivation::Relu,
           gra::neurojac::FlowActivation::Silu,
           gra::neurojac::FlowActivation::Tanh}) {
    CAPTURE(gra::neurojac::FlowActivationName(activation));
    gra::neurojac::MNeuroJacConfig config;
    config.layers = 3;
    config.bins = 4;
    config.hidden = {6, 6};
    config.activation = activation;
    gra::neurojac::MSplineFlow flow;
    flow.Configure(3, config, 71);

    std::vector<double> parameters = flow.Parameters();
    for (std::size_t i = 0; i < parameters.size(); ++i) {
      parameters[i] += 3e-3 * std::sin(static_cast<double>(i + 1));
    }
    flow.SetParameters(std::move(parameters));

    const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
        {{0.21, 0.47, 0.79}, 1.0, 0.0}, {{0.68, 0.31, 0.56}, 1.0, 0.0}};
    const std::vector<double> coefficients = {0.4, 0.6};
    const std::vector<std::size_t> indices = {0, 1};
    const std::unique_ptr<gra::neurojac::MFlowObjective> loss_function =
        gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::ForwardKL,
                                           1.0);
    std::vector<double> gradient;
    const double loss =
        flow.LossGradient(events, coefficients, indices, events.size(),
                          *loss_function, 0.05, gradient);
    CHECK(std::isfinite(loss));
    CHECK(std::all_of(gradient.begin(), gradient.end(),
                      [](double value) { return std::isfinite(value); }));
    CHECK(std::any_of(gradient.begin(), gradient.end(),
                      [](double value) { return std::abs(value) > 1e-10; }));
  }
}

TEST_CASE("Log-domain mixture loss stays finite above the exponential range",
          "[neurojac][gradient][stability]") {
  gra::neurojac::MNeuroJacConfig config;
  config.layers = 120;
  config.bins = 4;
  config.hidden = {};
  gra::neurojac::MSplineFlow flow;
  flow.Configure(4, config, 123);

  std::vector<double> parameters = flow.Parameters();
  const std::size_t stride = 3 * config.bins + 1;
  for (const gra::neurojac::MSplineCoupling &coupling : flow.Couplings()) {
    const auto &output = coupling.dense_shapes_.back();
    for (std::size_t coordinate = 0;
         coordinate < coupling.transform_indices_.size(); ++coordinate) {
      for (std::size_t knot = 0; knot <= config.bins; ++knot) {
        parameters[output.bias_offset + coordinate * stride + 2 * config.bins +
                   knot] = -1000.0;
      }
    }
  }
  flow.SetParameters(std::move(parameters));

  const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
      {{0.0, 0.0, 0.0, 0.0}, 1.0, 0.0}};
  const std::vector<double> coefficients = {1.0};
  const std::vector<std::size_t> indices = {0};
  REQUIRE(flow.LogDensity(events.front().point) >
          std::log(std::numeric_limits<double>::max()));
  for (gra::neurojac::FlowLoss loss : std::vector<gra::neurojac::FlowLoss>{
           gra::neurojac::FlowLoss::ForwardKL, gra::neurojac::FlowLoss::Chi2}) {
    const std::unique_ptr<gra::neurojac::MFlowObjective> loss_function =
        gra::neurojac::CreateFlowObjective(loss, 1.0);
    std::vector<double> gradient;
    const double objective =
        flow.LossGradient(events, coefficients, indices, events.size(),
                          *loss_function, 0.05, gradient);
    CHECK(std::isfinite(objective));
    CHECK(std::all_of(gradient.begin(), gradient.end(),
                      [](double value) { return std::isfinite(value); }));
  }
}

// Resolve a small mixture probability multiplied by a sharply concentrated density
TEST_CASE("Spline mixture densities retain small component probabilities",
          "[neurojac][mixture][gradient][stability]") {
  gra::neurojac::MNeuroJacConfig config;
  config.flow_components = 2;
  config.layers = 120;
  config.bins = 4;
  config.hidden = {};
  gra::neurojac::MSplineFlowMixture mixture;
  mixture.Configure(4, config, 123);
  auto parameters = mixture.Parameters();
  const std::size_t stride = 3 * config.bins + 1;
  for (const auto &coupling : mixture.Component(0).Couplings()) {
    const auto &output = coupling.dense_shapes_.back();
    for (const auto &coordinate : indices(coupling.transform_indices_)) {
      for (std::size_t knot = 0; knot <= config.bins; ++knot) {
        parameters[output.bias_offset + coordinate * stride +
                   2 * config.bins + knot] = -1000.0;
      }
    }
  }
  parameters.back() = -1000.0;
  mixture.SetParameters(parameters);

  const std::vector<double> point(4, 1e-200);
  const double expected =
      mixture.Component(0).LogDensity(point) + parameters.back();
  REQUIRE(expected > 100.0);
  CHECK(mixture.LogDensity(point) == Approx(expected).margin(1e-10));
  CHECK(mixture.LogDensityBatch(point).front() == Approx(expected).margin(1e-10));

  const auto objective = gra::neurojac::CreateFlowObjective(
      gra::neurojac::FlowLoss::ForwardKL, 1.0);
  std::vector<double> gradient;
  const double loss = mixture.LossGradient(
      {{point, 1.0, 0.0}}, {1.0}, {0}, 1, *objective, 0.05, gradient);
  CHECK(loss == Approx(-expected - std::log1p(-0.05)).margin(1e-10));
  CHECK(gradient.back() == Approx(-1.0).margin(1e-12));
  CHECK(gra::AllFinite(gradient));
}

// Reweight a dominant final component without subtracting probabilities near one
TEST_CASE("Spline mixture reweighting preserves suppressed component ratios",
          "[neurojac][mixture][stability]") {
  gra::neurojac::MNeuroJacConfig config;
  config.flow_components = 3;
  gra::neurojac::MSplineFlowMixture mixture;
  mixture.Configure(1, config, 17);
  auto parameters = mixture.Parameters();
  parameters[parameters.size() - 2] = -1000.0;
  parameters.back() = -1001.0;
  mixture.SetParameters(parameters);
  mixture.SetLastWeight(0.25);
  const auto weights = mixture.Weights();
  REQUIRE(gra::AllFinite(weights));
  CHECK(weights[0] == Approx(0.75 / (1.0 + std::exp(-1.0))).margin(1e-13));
  CHECK(weights[1] == Approx(0.75 / (1.0 + std::exp(1.0))).margin(1e-13));
  CHECK(weights[2] == Approx(0.25).margin(1e-13));
}

TEST_CASE("NUMERICS activation and learning-rate schedule options are applied",
          "[neurojac][config][schedule]") {
  const std::vector<std::pair<std::string, gra::neurojac::FlowActivation>>
      activations = {{"relu", gra::neurojac::FlowActivation::Relu},
                     {"silu", gra::neurojac::FlowActivation::Silu},
                     {"tanh", gra::neurojac::FlowActivation::Tanh}};
  for (const auto &[name, expected] : activations) {
    gra::neurojac::MFlowIntegrator integrator;
    integrator.ReadParameters(
        WriteNeuroJacCard("forward_kl", "random", 0.0, name));
    integrator.Configure(2, 73);
    CHECK(integrator.Config().activation == expected);
  }

  const std::vector<std::pair<std::string, gra::neurojac::FlowLRSchedule>>
      schedules = {{"none", gra::neurojac::FlowLRSchedule::None},
                   {"cosine", gra::neurojac::FlowLRSchedule::Cosine},
                   {"exponential", gra::neurojac::FlowLRSchedule::Exponential}};
  for (const auto &[name, expected] : schedules) {
    CAPTURE(name);
    gra::neurojac::MFlowIntegrator integrator;
    integrator.ReadParameters(WriteNeuroJacCard("forward_kl", "random", 0.0,
                                                "silu", name, 5, 1e-3, 1, 1.0,
                                                2.0, 0.7, false));
    integrator.Configure(2, 79);
    CHECK(integrator.Config().lr_schedule == expected);

    gra::MRandom random;
    random.SetSeed(83);
    integrator.BufferEvents(TrainingEvents(
        integrator.SampleBatch(random, integrator.Config().buffer_size)));
    const gra::neurojac::MFlowTrainingReport report = integrator.TrainRound();
    const double expected_last_rate =
        expected == gra::neurojac::FlowLRSchedule::None
            ? integrator.Config().learning_rate
            : integrator.Config().final_learning_rate;
    CHECK(report.learning_rate == Approx(expected_last_rate).epsilon(1e-12));
  }
}

TEST_CASE(
    "AdamW trains all NEUROJAC objectives and frozen sampling stays unbiased",
    "[neurojac][adamw]") {
  for (const std::string &loss :
       std::vector<std::string>{"forward_kl", "chi2", "renyi"}) {
    CAPTURE(loss);
    gra::neurojac::MFlowIntegrator integrator;
    integrator.ReadParameters(WriteNeuroJacCard(
        loss, "random", 0.0, "silu", "none", 5, 1e-3, 1, 1.0, 2.0, 0.7, false));
    integrator.Configure(2, 31415);
    const std::vector<double> initial_parameters =
        integrator.Flow().Parameters();

    gra::MRandom random;
    random.SetSeed(2718);
    std::vector<gra::neurojac::MFlowTrainingEvent> events;
    events.reserve(integrator.Config().buffer_size);
    for (std::size_t i = 0; i < integrator.Config().buffer_size; ++i) {
      const gra::neurojac::MFlowSample sample = integrator.Sample(random);
      events.push_back(
          {sample.point, TrainingTarget(sample.point), sample.log_density});
    }
    integrator.BufferEvents(events);
    const gra::neurojac::MFlowTrainingReport report = integrator.TrainRound();

    CHECK(std::isfinite(report.loss));
    CHECK(report.buffered_events == integrator.Config().buffer_size);
    CHECK(report.nonzero_events == report.buffered_events);
    CHECK(report.collection_ess_fraction > 0.0);

    double parameter_change = 0.0;
    for (std::size_t i = 0; i < initial_parameters.size(); ++i) {
      parameter_change +=
          std::abs(integrator.Flow().Parameters()[i] - initial_parameters[i]);
    }
    CHECK(parameter_change > 0.0);

    integrator.Freeze();
    CHECK(integrator.IsFrozen());
    CHECK_THROWS_AS(integrator.BufferEvents(events), std::logic_error);

    constexpr std::size_t samples = 6000;
    const gra::neurojac::MFlowSampleBatch generated =
        integrator.SampleBatch(random, samples);
    REQUIRE(generated.Size() == samples);
    double weighted_sum = 0.0;
    double weighted_sum2 = 0.0;
    for (std::size_t event = 0; event < generated.Size(); ++event) {
      const std::vector<double> point = BatchPoint(generated, event);
      CHECK(generated.log_densities[event] ==
            Approx(integrator.LogDensity(point)).margin(2e-9));
      const double weight =
          TrainingTarget(point) * std::exp(-generated.log_densities[event]);
      weighted_sum += weight;
      weighted_sum2 += weight * weight;
    }
    const double estimate = weighted_sum / static_cast<double>(samples);
    const double variance =
        (weighted_sum2 - static_cast<double>(samples) * estimate * estimate) /
        static_cast<double>(samples - 1);
    const double standard_error =
        std::sqrt(std::max(0.0, variance) / static_cast<double>(samples));
    CHECK(std::isfinite(standard_error));
    CHECK(standard_error > 0.0);
    CHECK(std::abs(estimate - 2.0) <= 5.0 * standard_error + 1e-12);
  }
}

TEST_CASE("Held-out validation restores the best flow checkpoint",
          "[neurojac][validation]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("forward_kl", "random", 0.0,
                                              "silu", "none", 1, 1.0, 2, 1.0,
                                              2.0, 0.7, false));
  integrator.Configure(2, 101);

  gra::MRandom random;
  random.SetSeed(103);
  std::vector<double> best_parameters;
  for (std::size_t round = 0; round < 2; ++round) {
    integrator.BufferEvents(TrainingEvents(
        integrator.SampleBatch(random, integrator.FreshEventCount())));
    integrator.BufferValidationEvents(TrainingEvents(
        integrator.SampleBatch(random, integrator.ValidationEventCount())));
    integrator.TrainRound();
    const gra::neurojac::MFlowValidationReport validation =
        integrator.ValidateRound();
    CHECK(validation.events == integrator.Config().validation_size);
    CHECK(validation.nonzero_events == validation.events);
    CHECK(validation.ess_fraction > 0.0);
    CHECK(validation.ess_fraction <= 1.0);
    CHECK(std::isfinite(validation.relative_variance));
    if (round == 0) {
      CHECK(validation.improved);
      CHECK(validation.best_round == 1);
      best_parameters = integrator.Flow().Parameters();
    } else {
      CHECK_FALSE(validation.improved);
      CHECK(validation.early_stop);
      CHECK(validation.best_round == 1);
    }
  }

  double distance_from_best = 0.0;
  for (std::size_t parameter = 0; parameter < best_parameters.size();
       ++parameter) {
    distance_from_best += std::abs(integrator.Flow().Parameters()[parameter] -
                                   best_parameters[parameter]);
  }
  CHECK(distance_from_best > 0.0);
  integrator.Freeze();
  CHECK(integrator.Flow().Parameters() == best_parameters);
  CHECK_THROWS_AS(integrator.BufferValidationEvents(
                      std::vector<gra::neurojac::MFlowTrainingEvent>{}),
                  std::logic_error);
}

// Check patience before the round limit and after adding the last flow
TEST_CASE("NEUROJAC stops stale training before the final round", "[neurojac][validation][regression]") {
  for (std::size_t components : {1U, 2U}) {
    CAPTURE(components);
    gra::neurojac::MFlowIntegrator integrator;
    integrator.ReadParameters(WriteNeuroJacCard("renyi", "interleave", 0.0, "silu", "none", 2, 1.0, 12, 1.0, 2.0, 0.7,
                                                false, components, false));
    integrator.Configure(2, 43);
    gra::MRandom random;
    random.SetSeed(47);
    const std::size_t   stop_round = components == 1 ? 3 : 8;
    std::vector<double> best;
    for (std::size_t round = 1; round <= stop_round; ++round) {
      integrator.BufferEvents(TrainingEvents(integrator.SampleBatch(random, integrator.FreshEventCount())));
      integrator.BufferValidationEvents(
          TrainingEvents(integrator.SampleBatch(random, integrator.ValidationEventCount())));
      integrator.TrainRound();
      const auto report = integrator.ValidateRound();
      CHECK(report.early_stop == (round == stop_round));
      CHECK(report.best_round == 1);
      if (round == 1) { best = integrator.Mixture().Parameters(); }
    }
    integrator.Freeze();
    CHECK(integrator.Mixture().ComponentCount() == 1);
    CHECK(integrator.Mixture().Parameters() == best);
  }
}

TEST_CASE("NEUROJAC replay mixes bounded history with fresh events",
          "[neurojac][replay]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("forward_kl"));
  integrator.Configure(2, 1618);
  integrator.InitializeVegasGrid(
      {{0.0, 0.25, 0.5, 0.75, 1.0}, {0.0, 0.25, 0.5, 0.75, 1.0}});

  gra::MRandom random;
  random.SetSeed(271828);
  integrator.BufferValidationEvents(TrainingEvents(
      integrator.SampleBatch(random, integrator.ValidationEventCount())));
  const gra::neurojac::MFlowValidationReport baseline =
      integrator.ValidateBaseline();
  CHECK(baseline.improved);
  CHECK(baseline.best_round == 0);

  for (std::size_t round = 0; round < 4; ++round) {
    const std::size_t fresh_events = integrator.FreshEventCount();
    CHECK(fresh_events == (round == 0 ? 128 : 64));
    const gra::neurojac::MFlowSampleBatch samples =
        integrator.SampleBatch(random, fresh_events);
    integrator.BufferEvents(TrainingEvents(samples));
    const gra::neurojac::MFlowTrainingReport report = integrator.TrainRound();

    CHECK(report.buffered_events == 128);
    CHECK(report.fresh_events == fresh_events);
    CHECK(report.replayed_events == (round == 0 ? 0 : 64));
    CHECK(report.vegas_initialized == (round == 0));
    CHECK(report.replay_pool_events ==
          std::min<std::size_t>(128 + round * 64, 256));
    CHECK(report.nonzero_events == report.buffered_events);
  }
}

TEST_CASE("NEUROJAC trains a sparse configured buffer without empty validation "
          "failure",
          "[neurojac][sparse]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("forward_kl", "random", 0.0,
                                              "silu", "cosine", 5, 1e-3, 1, 1.0,
                                              2.0, 0.7, false));
  integrator.Configure(2, 2026);

  gra::MRandom random;
  random.SetSeed(1792394);
  CHECK(integrator.FreshEventCount() == integrator.Config().buffer_size);
  CHECK(integrator.ValidationEventCount() ==
        integrator.Config().validation_size);
  const std::size_t configured_events = integrator.Config().buffer_size;
  const gra::neurojac::MFlowSampleBatch samples =
      integrator.SampleBatch(random, configured_events);
  const std::vector<double> supported_point =
      BatchPoint(samples, samples.Size() - 1);
  const double expected_integral =
      gra::statistics::ImportanceWeight(TrainingTarget(supported_point),
                                        -samples.log_densities.back()) /
      static_cast<double>(configured_events);
  std::vector<gra::neurojac::MFlowTrainingEvent> training;
  training.reserve(samples.Size());
  for (std::size_t i = 0; i < samples.Size(); ++i) {
    std::vector<double> point = BatchPoint(samples, i);
    const double target = i + 1 == samples.Size() ? TrainingTarget(point) : 0.0;
    training.push_back({std::move(point), target, samples.log_densities[i]});
  }
  integrator.BufferEvents(std::move(training));

  std::vector<gra::neurojac::MFlowTrainingEvent> validation;
  const gra::neurojac::MFlowSampleBatch validation_samples =
      integrator.SampleBatch(random, integrator.ValidationEventCount());
  validation.reserve(validation_samples.Size());
  for (std::size_t event = 0; event < validation_samples.Size(); ++event) {
    validation.push_back({BatchPoint(validation_samples, event), 0.0,
                          validation_samples.log_densities[event]});
  }
  integrator.BufferValidationEvents(std::move(validation));

  const gra::neurojac::MFlowTrainingReport training_report =
      integrator.TrainRound();
  CHECK(training_report.buffered_events == configured_events);
  CHECK(training_report.nonzero_events == 1);
  CHECK(training_report.integral == Approx(expected_integral).epsilon(1e-12));
  CHECK(training_report.relative_error == Approx(1.0).epsilon(1e-12));
  CHECK(std::isfinite(training_report.loss));
  CHECK(training_report.learning_rate ==
        Approx(integrator.Config().final_learning_rate).epsilon(1e-12));
  CHECK(training_report.collection_ess_fraction ==
        Approx(1.0 / static_cast<double>(configured_events)).epsilon(1e-12));

  const gra::neurojac::MFlowValidationReport validation_report =
      integrator.ValidateRound();
  CHECK(validation_report.events == integrator.Config().validation_size);
  CHECK(validation_report.nonzero_events == 0);
  CHECK(gra::math::IsZero(validation_report.ess_fraction));
  CHECK(std::isinf(validation_report.relative_variance));
  CHECK_FALSE(validation_report.improved);
  CHECK(validation_report.best_round == 0);
  CHECK(validation_report.stale_rounds == 0);
  CHECK_FALSE(validation_report.early_stop);

  integrator.BufferValidationEvents(TrainingEvents(
      integrator.SampleBatch(random, integrator.ValidationEventCount())));
  const gra::neurojac::MFlowValidationReport recovered_validation =
      integrator.ValidateRound();
  CHECK(recovered_validation.nonzero_events == recovered_validation.events);
  CHECK(recovered_validation.improved);
  CHECK(recovered_validation.best_round == 2);
  CHECK(integrator.FreshEventCount() == integrator.Config().buffer_size / 2);
  integrator.Freeze();
}

TEST_CASE("NEUROJAC adds residual flows stagewise and freezes earlier flows",
          "[neurojac][boosting]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard(
      "renyi", "interleave", 0.0, "silu", "none", 5, 1e-3, 4, 1.0,
      2.0, 0.7, false, 2, false));
  integrator.Configure(2, 31415);

  gra::MRandom random;
  random.SetSeed(92653);
  std::vector<double> frozen_first;
  std::vector<double> new_after_add;
  for (std::size_t round = 0; round < 4; ++round) {
    integrator.BufferEvents(TrainingEvents(
        integrator.SampleBatch(random, integrator.FreshEventCount())));
    const gra::neurojac::MFlowTrainingReport report = integrator.TrainRound();
    if (round < 2) {
      REQUIRE(report.component_weights.size() == 1);
      CHECK_FALSE(report.component_added);
      if (round == 1) {
        frozen_first = integrator.Mixture().Component(0).Parameters();
      }
    } else {
      REQUIRE(report.component_weights.size() == 2);
      CHECK(report.component_added == (round == 2));
      CHECK(report.component_weights[0] + report.component_weights[1] ==
            Approx(1.0).margin(1e-14));
      CHECK(std::abs(report.component_weights[0] -
                     report.component_weights[1]) > 1e-3);
      CHECK(integrator.Mixture().Component(0).Parameters() == frozen_first);
      CHECK(report.boost_weight > 0.0);
      CHECK(report.boost_weight < 1.0);
      if (round == 2) {
        CHECK(report.component_weights[1] >= 0.05 - 1e-14);
        new_after_add = integrator.Mixture().Component(1).Parameters();
      } else {
        CHECK(integrator.Mixture().Component(1).Parameters() != new_after_add);
      }
    }
  }

  integrator.Freeze();
  gra::neurojac::MFlowIntegrator restored;
  restored.DeserializeModel(integrator.SerializeModel(), 2);
  CHECK(restored.Mixture().ComponentCount() == 2);
  CHECK(restored.Mixture().Weights() == integrator.Mixture().Weights());
  const gra::neurojac::MFlowSampleBatch samples =
      integrator.SampleBatch(random, 1024);
  for (std::size_t event = 0; event < samples.Size(); ++event) {
    const std::vector<double> point = BatchPoint(samples, event);
    CHECK(samples.log_densities[event] ==
          Approx(integrator.LogDensity(point)).margin(2e-9));
    CHECK(restored.LogDensity(point) ==
          Approx(integrator.LogDensity(point)).margin(1e-13));
  }
}

// Check that a residual flow learns when the old density is large in a small volume
TEST_CASE("Superposition trains new flows on narrow support", "[neurojac][boosting][regression]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(
      WriteNeuroJacCard("renyi", "interleave", 0.0, "silu", "none", 5, 1e-3, 4, 1.0, 2.0, 0.7, true, 2, true));
  integrator.Configure(4, 79);
  integrator.InitializeVegasGrid(std::vector<std::vector<double>>(4, {0.0, 1e-4, 2e-4, 3e-4, 4e-4, 5e-4, 1.0}));
  gra::MRandom random;
  random.SetSeed(83);
  for (std::size_t round = 0; round < 3; ++round) {
    const auto samples = integrator.SampleBatch(random, integrator.FreshEventCount());
    std::vector<gra::neurojac::MFlowTrainingEvent> events;
    for (std::size_t i = 0; i < samples.Size(); ++i) {
      const auto   point  = BatchPoint(samples, i);
      const double target = std::all_of(point.begin(), point.end(), [](double x) { return x < 0.006; }) ? 1.0 : 0.0;
      events.push_back({point, target, samples.log_densities[i]});
    }
    integrator.BufferEvents(std::move(events));
    const auto report = integrator.TrainRound();
    REQUIRE(report.nonzero_events > 0);
  }
  REQUIRE(integrator.Mixture().ComponentCount() == 2);
  CHECK(integrator.Mixture().Component(1).LogDensity({0.001, 0.001, 0.001, 0.001}) > 0.01);
}

TEST_CASE("Frozen NEUROJAC proposals round-trip through the generic serializer",
          "[neurojac][serialization]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(WriteNeuroJacCard("renyi", "interleave", 0.0,
                                              "tanh", "cosine", 5, 1e-3, 1, 1.1,
                                              2.3, 0.6, false, 2, false));
  integrator.Configure(2, 31415);
  CHECK_FALSE(integrator.Config().vegas_init);
  CHECK_THROWS_AS(integrator.SerializeModel(), std::logic_error);

  gra::MRandom adaptation_random;
  adaptation_random.SetSeed(92653);
  integrator.BufferEvents(TrainingEvents(
      integrator.SampleBatch(adaptation_random, integrator.FreshEventCount())));
  integrator.BufferValidationEvents(TrainingEvents(integrator.SampleBatch(
      adaptation_random, integrator.ValidationEventCount())));
  const gra::neurojac::MFlowTrainingReport training_report =
      integrator.TrainRound();
  CHECK(training_report.loss_alpha == Approx(2.3));
  CHECK_FALSE(training_report.vegas_initialized);
  REQUIRE(training_report.component_weights.size() == 1);
  CHECK(training_report.component_weights[0] == Approx(1.0));
  CHECK_FALSE(training_report.component_added);
  integrator.ValidateRound();
  integrator.Freeze();

  const std::string serialized = integrator.SerializeModel();
  const nlohmann::json document = nlohmann::json::parse(serialized);
  CHECK(document.at("MODEL_FORMAT") == "GRANIITTI_NEUROJAC");
  CHECK(document.at("MODEL_VERSION") == 1);
  CHECK_FALSE(document.at("CONFIG").at("vegas_init").get<bool>());
  CHECK_FALSE(document.at("CONFIG").at("vegas_spline_init").get<bool>());
  CHECK(document.at("CONFIG").at("flow_components") == 2);
  CHECK(document.at("CONFIG").at("vegas_ncall") == 128);
  CHECK(document.at("CONFIG").at("vegas_rounds") == 2);
  CHECK(document.at("CONFIG").at("vegas_sampler_rounds") == 0);
  CHECK(document.at("ACTIVE_COMPONENTS") == 1);
  CHECK(document.at("DIMENSION") == 2);
  CHECK(document.at("PARAMETER_COUNT") ==
        integrator.Mixture().Parameters().size());
  CHECK(document.at("PARAMETERS").size() ==
        integrator.Mixture().Parameters().size());

  const std::string filename = "tmp/neurojac/nested/model.json";
  integrator.SaveModel(filename);
  gra::neurojac::MFlowIntegrator restored;
  restored.LoadModel(filename, 2);
  gra::neurojac::MFlowIntegrator adaptive;
  adaptive.Configure(2, 2718);
  CHECK_THROWS_AS(adaptive.SaveModel(filename), std::logic_error);
  gra::neurojac::MFlowIntegrator preserved;
  preserved.LoadModel(filename, 2);
  CHECK(preserved.SerializeModel() == serialized);
  CHECK(restored.IsConfigured());
  CHECK(restored.IsFrozen());
  CHECK(restored.Config().loss.kind == gra::neurojac::FlowLoss::Renyi);
  CHECK(restored.Config().loss.alpha_initial == Approx(1.1));
  CHECK(restored.Config().loss.alpha_final == Approx(2.3));
  CHECK(restored.Config().loss.annealing_fraction == Approx(0.6));
  CHECK(restored.Config().activation == gra::neurojac::FlowActivation::Tanh);
  CHECK(restored.Config().lr_schedule == gra::neurojac::FlowLRSchedule::Cosine);
  CHECK(restored.Config().permutation ==
        gra::neurojac::FlowPermutation::Interleave);
  CHECK_FALSE(restored.Config().vegas_init);
  CHECK_FALSE(restored.Config().vegas_spline_init);
  CHECK(restored.Config().flow_components == 2);
  CHECK(restored.Config().vegas_ncall == 128);
  CHECK(restored.Config().vegas_rounds == 2);
  CHECK(restored.Config().vegas_sampler_rounds == 0);
  CHECK(restored.Mixture().ComponentCount() == 1);
  CHECK(restored.Mixture().Parameters() == integrator.Mixture().Parameters());

  for (const std::vector<double> &point : std::vector<std::vector<double>>{
           {0.13, 0.27}, {0.51, 0.79}, {0.91, 0.07}}) {
    CHECK(restored.LogDensity(point) ==
          Approx(integrator.LogDensity(point)).margin(1e-13));
  }

  gra::MRandom original_random;
  gra::MRandom restored_random;
  original_random.SetSeed(58979);
  restored_random.SetSeed(58979);
  const auto original_samples = integrator.SampleBatch(original_random, 32);
  const auto restored_samples = restored.SampleBatch(restored_random, 32);
  REQUIRE(restored_samples.dimension == original_samples.dimension);
  CHECK(restored_samples.points == original_samples.points);
  REQUIRE(restored_samples.log_densities.size() ==
          original_samples.log_densities.size());
  for (std::size_t i = 0; i < original_samples.Size(); ++i) {
    CHECK(restored_samples.log_densities[i] ==
          Approx(original_samples.log_densities[i]).margin(1e-13));
  }

  gra::MRandom independent_node_random;
  independent_node_random.SetSeed(58980);
  const auto independent_samples =
      restored.SampleBatch(independent_node_random, 32);
  bool different_stream = false;
  for (std::size_t i = 0; i < independent_samples.Size(); ++i) {
    const std::vector<double> independent_point =
        BatchPoint(independent_samples, i);
    different_stream = different_stream ||
                       independent_point != BatchPoint(restored_samples, i);
    CHECK(independent_samples.log_densities[i] ==
          Approx(restored.LogDensity(independent_point)).margin(1e-13));
  }
  CHECK(different_stream);
  CHECK(restored.Mixture().Parameters() == integrator.Mixture().Parameters());

  gra::neurojac::MFlowIntegrator wrong_dimension;
  CHECK_THROWS_AS(wrong_dimension.LoadModel(filename, 3),
                  std::invalid_argument);
  nlohmann::json truncated = document;
  truncated.at("PARAMETERS").erase(truncated.at("PARAMETERS").end() - 1);
  gra::neurojac::MFlowIntegrator malformed;
  CHECK_THROWS_AS(malformed.DeserializeModel(truncated.dump(), 2),
                  std::invalid_argument);
}

// Reject invalid counts before unsigned conversion in cards and saved proposals
TEST_CASE("NEUROJAC rejects nonintegral and negative counts", "[neurojac][config]") {
  const auto filename = WriteNeuroJacCard("renyi");
  const auto document = nlohmann::json::parse(gra::aux::GetInputData(filename));
  for (const std::string key :
       {"flow_components", "layers", "bins", "buffer_size", "batch_size", "replay_cap", "val_size", "val_patience",
        "rounds", "epochs", "vegas_ncall", "vegas_rounds", "vegas_sampler_rounds", "hidden"}) {
    for (const nlohmann::json &value : {nlohmann::json(-1), nlohmann::json(1.5), nlohmann::json(true)}) {
      CAPTURE(key, value);
      auto invalid                      = document;
      invalid["NUMERICS_NEUROJAC"][key] = key == "hidden" ? nlohmann::json::array({value}) : value;
      std::ofstream(filename) << invalid;
      gra::neurojac::MFlowIntegrator integrator;
      CHECK_THROWS_AS(integrator.ReadParameters(filename), std::invalid_argument);
    }
  }
  gra::neurojac::MNeuroJacConfig config;
  config.rounds = std::numeric_limits<std::size_t>::max();
  CHECK_THROWS_AS(config.Validate(1), std::invalid_argument);

  gra::neurojac::MFlowIntegrator integrator;
  integrator.Configure(1, 13);
  integrator.Freeze();
  const auto model = nlohmann::json::parse(integrator.SerializeModel());
  for (const std::string key : {"MODEL_VERSION", "DIMENSION", "FLOW_SEED", "ACTIVE_COMPONENTS", "PARAMETER_COUNT"}) {
    auto invalid = model;
    invalid[key] = 1.5;
    CHECK_THROWS_AS(integrator.DeserializeModel(invalid.dump(), 1), std::invalid_argument);
  }
}

// Compare checkpoints on common data even when reported ESS is saturated
TEST_CASE("NEUROJAC replaces an optimistic validation baseline", "[neurojac][validation]") {
  const auto filename =
      WriteNeuroJacCard("forward_kl", "interleave", 0.0, "silu", "none", 5, 1e-3, 4, 1.0, 2.0, 0.7, true, 1, false);
  auto  card             = nlohmann::json::parse(gra::aux::GetInputData(filename));
  auto &block            = card["NUMERICS_NEUROJAC"];
  block["layers"]        = 2;
  block["bins"]          = 4;
  block["val_size"]      = 1;
  block["epochs"]        = 4;
  block["learning_rate"] = 0.03;
  std::ofstream(filename) << card;
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(filename);
  integrator.Configure(1, 1);
  integrator.InitializeVegasGrid({{0.0, 0.5, 1.0}});
  integrator.BufferValidationEvents({{{0.25}, std::exp(1.25), 0.0}});
  CHECK(integrator.ValidateBaseline().ess_fraction == Approx(1.0));
  bool improved = false;
  for (std::size_t round = 0; round < integrator.Config().rounds; ++round) {
    std::vector<gra::neurojac::MFlowTrainingEvent> events;
    const auto                                     count = integrator.FreshEventCount();
    for (std::size_t i = 0; i < count; ++i) {
      const double x = (i + 0.5) / count;
      events.push_back({{x}, std::exp(5.0 * x), 0.0});
    }
    integrator.BufferEvents(std::move(events));
    integrator.TrainRound();
    integrator.BufferValidationEvents({{{0.9}, std::exp(4.5), 0.0}});
    improved = integrator.ValidateRound().improved || improved;
  }
  CHECK(improved);
  integrator.Freeze();
  double                second_moment = 0.0;
  constexpr std::size_t points        = 2048;
  for (std::size_t i = 0; i < points; ++i) {
    const double x = (i + 0.5) / points;
    second_moment += std::exp(10.0 * x - integrator.LogDensity({x})) / points;
  }
  const double integral = std::expm1(5.0) / 5.0;
  CHECK(integral * integral / second_moment > 0.7);
}

// Check successive optimizer steps against freshly evaluated Renyi gradients
TEST_CASE("NEUROJAC refreshes Renyi coefficients between minibatches", "[neurojac][gradient]") {
  const auto filename    = WriteNeuroJacCard("renyi", "none", 0.0, "silu", "none", 5, 1e-3, 1, 2.0, 2.0, 1.0, false);
  auto       card        = nlohmann::json::parse(gra::aux::GetInputData(filename));
  auto      &block       = card["NUMERICS_NEUROJAC"];
  block["layers"]        = 1;
  block["bins"]          = 2;
  block["buffer_size"]   = 4;
  block["batch_size"]    = 2;
  block["epochs"]        = 1;
  block["learning_rate"] = 0.1;
  block["replay_frac"]   = 0.0;
  block["adamw_beta1"]   = 0.0;
  block["adamw_beta2"]   = 0.0;
  block["adamw_weight_decay"] = 0.0;
  std::ofstream(filename) << card;
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(filename);
  constexpr std::uint32_t seed = 19;
  integrator.Configure(1, seed);
  const auto                &config = integrator.Config();
  gra::neurojac::MSplineFlow reference;
  reference.Configure(1, config, seed);
  const std::vector<gra::neurojac::MFlowTrainingEvent> events = {
      {{0.13}, 1.0, 0.0}, {{0.39}, 2.0, 0.0}, {{0.61}, 3.0, 0.0}, {{0.87}, 8.0, 0.0}};
  std::vector<std::size_t> order = {0, 1, 2, 3};
  std::mt19937             random(seed ^ 0x9e3779b9U);
  std::shuffle(order.begin(), order.end(), random);
  const auto objective = gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::Renyi, 2.0);
  for (std::size_t batch = 0; batch < order.size(); batch += 2) {
    const auto          coefficients = RenyiEscortCoefficients(reference, events, 2.0, config.uniform_mix);
    std::vector<double> gradient;
    reference.LossGradient(events, coefficients, {order[batch], order[batch + 1]}, events.size(), *objective,
                           config.uniform_mix, gradient);
    gra::statistics::ClipEuclideanNorm(gradient, config.gradient_clip);
    for (double &value : gradient) { value = config.learning_rate * value / (std::abs(value) + config.adamw_epsilon); }
    reference.ApplyParameterStep(gradient, 0.0, 0.0);
  }
  integrator.BufferEvents(events);
  const auto report = integrator.TrainRound();
  REQUIRE(report.rejected_epochs == 0);
  for (const auto &i : indices(reference.Parameters())) {
    CHECK(integrator.Flow().Parameters()[i] == Approx(reference.Parameters()[i]).margin(1e-12));
  }
}

// Reject composed splines that lose numerical invertibility before production
TEST_CASE("NEUROJAC rejects inconsistent saved proposal densities", "[neurojac][stability]") {
  gra::neurojac::MFlowIntegrator integrator;
  integrator.Configure(1, 13);
  integrator.Freeze();
  auto model  = nlohmann::json::parse(integrator.SerializeModel());
  auto config = integrator.Config();
  config.bins = 2;
  gra::neurojac::MSplineFlow flow;
  flow.Configure(1, config, 13);
  auto parameters = flow.Parameters();
  for (const auto &coupling : flow.Couplings()) {
    const auto bias = coupling.dense_shapes_.back().bias_offset;
    for (std::size_t knot = 0; knot <= config.bins; ++knot) { parameters[bias + 2 * config.bins + knot] = 1e8; }
  }
  flow.SetParameters(parameters);
  const auto sample = flow.Forward({0.13});
  REQUIRE(std::abs(sample.log_density - flow.LogDensity(sample.point)) > config.density_tolerance);
  model["CONFIG"]["bins"]  = config.bins;
  model["PARAMETERS"]      = parameters;
  model["PARAMETER_COUNT"] = parameters.size();
  CHECK_THROWS_AS(integrator.DeserializeModel(model.dump(), 1), std::invalid_argument);
  CHECK(integrator.LogDensity({0.13}) == Approx(0.0).margin(1e-12));
}

// Restore the previous proposal when an optimizer epoch is numerically invalid
TEST_CASE("NEUROJAC rolls back unstable optimizer epochs", "[neurojac][stability]") {
  const auto filename = WriteNeuroJacCard("chi2", "none", 0.0, "silu", "none", 5, 1e-3, 1, 1.0, 2.0, 0.7, false);
  auto       card     = nlohmann::json::parse(gra::aux::GetInputData(filename));
  card["NUMERICS_NEUROJAC"]["learning_rate"] = std::numeric_limits<double>::max();
  std::ofstream(filename) << card;
  gra::neurojac::MFlowIntegrator integrator;
  integrator.ReadParameters(filename);
  integrator.Configure(2, 19);
  const auto   parameters = integrator.Mixture().Parameters();
  gra::MRandom random;
  random.SetSeed(23);
  integrator.BufferEvents(TrainingEvents(integrator.SampleBatch(random, integrator.FreshEventCount())));
  const auto report = integrator.TrainRound();
  CHECK(report.rejected_epochs == integrator.Config().epochs);
  CHECK(integrator.Mixture().Parameters() == parameters);
  CHECK_NOTHROW(integrator.Freeze());
}

// Preserve optimizer trajectories and sampled densities across training worker counts
TEST_CASE("NEUROJAC parallel density evaluation preserves training", "[neurojac][threading][gradient]") {
  const std::size_t threads = GENERATE(1U, 4U);
  for (const std::size_t components : {1U, 2U}) {
    const auto filename = WriteNeuroJacCard("renyi", "random", 0.0, "silu", "none", 5, 1e-3, 4,
                                          1.0, 2.0, 0.7, false, components, false);
    gra::neurojac::MFlowIntegrator serial, parallel;
    serial.ReadParameters(filename);
    parallel.ReadParameters(filename);
    serial.Configure(4, 29);
    parallel.Configure(4, 29);
    gra::MRandom random;
    random.SetSeed(31);
    for (std::size_t round = 0; round < serial.Config().rounds; ++round) {
      const auto events = TrainingEvents(serial.SampleBatch(random, serial.FreshEventCount()));
      serial.BufferEvents(events);
      parallel.BufferEvents(events);
      const auto reference = serial.TrainRound();
      CAPTURE(threads, components, round);
      const auto result = parallel.TrainRound(threads);
      CHECK(result.rejected_epochs == reference.rejected_epochs);
      CHECK(result.loss == Approx(reference.loss).epsilon(0.0).margin(1e-12));
      const auto expected = serial.Mixture().Parameters();
      const auto actual = parallel.Mixture().Parameters();
      REQUIRE(actual.size() == expected.size());
      for (const auto &i : indices(expected)) { CHECK(actual[i] == Approx(expected[i]).epsilon(0.0).margin(1e-12)); }
    }
    serial.Freeze();
    parallel.Freeze();
    gra::MRandom first, second;
    first.SetSeed(37);
    second.SetSeed(37);
    const auto expected = serial.SampleBatch(first, 128);
    const auto actual = parallel.SampleBatch(second, 128);
    for (const auto &i : indices(expected.points)) {
      CHECK(actual.points[i] == Approx(expected.points[i]).epsilon(0.0).margin(1e-12));
    }
    for (const auto &i : indices(expected.log_densities)) {
      CHECK(actual.log_densities[i] == Approx(expected.log_densities[i]).epsilon(0.0).margin(1e-11));
    }
  }
}
