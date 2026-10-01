// NEUROJAC scalar and batched kernel benchmark
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <chrono>
#include <cmath>
#include <cstddef>
#include <iomanip>
#include <iostream>
#include <memory>
#include <numeric>
#include <random>
#include <vector>

#include "Graniitti/Sampling/MNeuroJac.h"

namespace {

// Compute elapsed milliseconds for one benchmark callback
template <typename Function> double Milliseconds(Function &&function) {
  const auto begin = std::chrono::steady_clock::now();
  function();
  return std::chrono::duration<double, std::milli>(
             std::chrono::steady_clock::now() - begin)
      .count();
}

// Construct deterministic training points and normalized coefficients
void MakeTrainingBatch(std::size_t dimension, std::size_t batch_size,
                       std::vector<gra::neurojac::MFlowTrainingEvent> &events,
                       std::vector<double> &coefficients,
                       std::vector<std::size_t> &indices) {
  std::mt19937_64 generator(173);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  events.resize(batch_size);
  coefficients.resize(batch_size);
  indices.resize(batch_size);
  for (std::size_t event = 0; event < batch_size; ++event) {
    events[event].point.resize(dimension);
    double target = 0.2;
    for (double &coordinate : events[event].point) {
      coordinate = uniform(generator);
      target += coordinate * coordinate;
    }
    events[event].target = target;
    events[event].collection_log_density = 0.0;
    coefficients[event] = target;
    indices[event] = event;
  }
  const double sum =
      std::accumulate(coefficients.begin(), coefficients.end(), 0.0);
  for (double &coefficient : coefficients) {
    coefficient /= sum;
  }
}

} // namespace

// Benchmark representative NEUROJAC training and frozen density workloads
int main() {
  constexpr std::size_t dimension = 6;
  constexpr std::size_t batch_size = 512;
  constexpr std::size_t samples = 16384;

  gra::neurojac::MNeuroJacConfig config;
  config.layers = 4;
  config.bins = 8;
  config.hidden = {32, 32};

  gra::neurojac::MSplineFlow flow;
  flow.Configure(dimension, config, 12345);
  std::vector<double> parameters = flow.Parameters();
  for (std::size_t i = 0; i < parameters.size(); ++i) {
    parameters[i] += 2e-3 * std::sin(static_cast<double>(i + 1));
  }
  flow.SetParameters(std::move(parameters));

  std::vector<gra::neurojac::MFlowTrainingEvent> events;
  std::vector<double> coefficients;
  std::vector<std::size_t> indices;
  MakeTrainingBatch(dimension, batch_size, events, coefficients, indices);

  const std::unique_ptr<gra::neurojac::MFlowObjective> loss_function =
      gra::neurojac::CreateFlowObjective(gra::neurojac::FlowLoss::ForwardKL,
                                         1.0);
  std::vector<double> gradient;
  flow.LossGradient(events, coefficients, indices, events.size(),
                    *loss_function, config.uniform_mix, gradient);
  double loss = 0.0;
  const double training_ms = Milliseconds([&] {
    loss = flow.LossGradient(events, coefficients, indices, events.size(),
                             *loss_function, config.uniform_mix, gradient);
  });

  std::mt19937_64 generator(821);
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  std::vector<double> points(samples * dimension);
  for (double &coordinate : points) {
    coordinate = uniform(generator);
  }

  double density_checksum = 0.0;
  const double density_ms = Milliseconds([&] {
    const std::vector<double> densities = flow.LogDensityBatch(points);
    density_checksum = std::accumulate(densities.begin(), densities.end(), 0.0);
  });

  double forward_checksum = 0.0;
  const double forward_ms = Milliseconds([&] {
    const gra::neurojac::MFlowSampleBatch generated = flow.ForwardBatch(points);
    for (std::size_t event = 0; event < generated.Size(); ++event) {
      forward_checksum += generated.log_densities[event] +
                          generated.points[event * generated.dimension];
    }
  });

  config.flow_components = 2;
  gra::neurojac::MSplineFlowMixture mixture;
  mixture.Configure(dimension, config, 12345);
  std::vector<double> mixture_parameters = mixture.Parameters();
  for (std::size_t i = 0; i < mixture_parameters.size(); ++i) {
    mixture_parameters[i] += 2e-3 * std::sin(static_cast<double>(i + 1));
  }
  mixture.SetParameters(mixture_parameters);

  std::vector<double> mixture_gradient;
  mixture.LossGradient(events, coefficients, indices, events.size(),
                       *loss_function, config.uniform_mix, mixture_gradient);
  double mixture_loss = 0.0;
  const double mixture_training_ms = Milliseconds([&] {
    mixture_loss = mixture.LossGradient(events, coefficients, indices,
                                        events.size(), *loss_function,
                                        config.uniform_mix, mixture_gradient);
  });

  double mixture_density_checksum = 0.0;
  const double mixture_density_ms = Milliseconds([&] {
    const std::vector<double> densities = mixture.LogDensityBatch(points);
    mixture_density_checksum =
        std::accumulate(densities.begin(), densities.end(), 0.0);
  });

  std::cout << std::fixed << std::setprecision(3)
            << "training_batch_ms=" << training_ms << " loss=" << loss << '\n'
            << "density_16384_ms=" << density_ms
            << " checksum=" << density_checksum << '\n'
            << "forward_16384_ms=" << forward_ms
            << " checksum=" << forward_checksum << '\n'
            << "mixture2_training_batch_ms=" << mixture_training_ms
            << " ratio=" << mixture_training_ms / training_ms
            << " loss=" << mixture_loss << '\n'
            << "mixture2_density_16384_ms=" << mixture_density_ms
            << " ratio=" << mixture_density_ms / density_ms
            << " checksum=" << mixture_density_checksum << '\n';
}
