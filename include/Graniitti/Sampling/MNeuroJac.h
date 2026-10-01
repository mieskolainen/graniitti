// GRANIITTI spline-flow importance sampler
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>

#ifndef MNEUROJAC_H
#define MNEUROJAC_H

#include <cstddef>
#include <cstdint>
#include <limits>
#include <memory>
#include <random>
#include <string>
#include <vector>

#include "Graniitti/Sampling/MNeuroJacLoss.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {

namespace neurojac {

// Supported hidden-layer activation functions
enum class FlowActivation { Relu, Silu, Tanh };

// Supported learning-rate schedules
enum class FlowLRSchedule { None, Cosine, Exponential };

// Supported parameter-free coordinate mixing between coupling transforms
enum class FlowPermutation { None, Interleave, Random };

// Numerical steering for the NEUROJAC spline-flow integrator
struct MNeuroJacConfig {
  // Maximum number of stagewise spline flows
  std::size_t flow_components = 1;
  std::size_t layers = 4;
  std::size_t bins = 8;
  std::vector<std::size_t> hidden = {32, 32};
  std::size_t buffer_size = 4096;
  std::size_t batch_size = 128;
  std::size_t replay_capacity = 16384;
  std::size_t validation_size = 2048;
  std::size_t validation_patience = 5;
  std::size_t rounds = 4;
  std::size_t epochs = 4;
  std::size_t vegas_ncall = 100000;
  std::size_t vegas_rounds = 10;
  std::size_t vegas_sampler_rounds = 3;
  bool vegas_init = true;
  bool vegas_spline_init = false;
  double learning_rate = 1e-3;
  double final_learning_rate = 1e-4;
  double validation_min_delta = 1e-3;
  double uniform_mix = 0.05;
  double replay_fraction = 0.5;
  double min_bin = 1e-3;
  double min_derivative = 1e-3;
  double density_tolerance = 1e-6;
  double gradient_clip = 10.0;
  double adamw_beta1 = 0.9;
  double adamw_beta2 = 0.999;
  double adamw_epsilon = 1e-8;
  double adamw_weight_decay = 1e-4;
  double adamw_decay_biases = 0.0;
  MFlowLossConfig loss;
  FlowActivation activation = FlowActivation::Silu;
  FlowLRSchedule lr_schedule = FlowLRSchedule::Cosine;
  FlowPermutation permutation = FlowPermutation::Random;

  // Validate all numerical parameters for a phase-space dimension
  void Validate(std::size_t dimension) const;
};

// One scalar spline evaluation and its logarithmic Jacobian
struct MSplineResult {
  double value = 0.0;
  double log_jacobian = 0.0;
};

// One exact proposal draw from the unit hypercube
struct MFlowSample {
  std::vector<double> point;
  double log_density = 0.0;
};

// Contiguous row-major proposal draws and their exact log densities
struct MFlowSampleBatch {
  std::vector<double> points;
  std::vector<double> log_densities;
  std::size_t dimension = 0;

  // Compute the number of proposal draws in the batch
  std::size_t Size() const { return log_densities.size(); }
};

// One buffered target evaluation with its collection density
struct MFlowTrainingEvent {
  std::vector<double> point;
  double target = 0.0;
  double collection_log_density = 0.0;
};

// Diagnostics returned after one adaptive training round
struct MFlowTrainingReport {
  double integral = 0.0;
  double relative_error = std::numeric_limits<double>::infinity();
  double loss = 0.0;
  double learning_rate = 0.0;
  double loss_alpha = 0.0;
  double collection_ess_fraction = 0.0;
  double maximum_abs_weight = 0.0;
  double boost_weight = 0.0;
  std::size_t buffered_events = 0;
  std::size_t fresh_events = 0;
  std::size_t replayed_events = 0;
  std::size_t replay_pool_events = 0;
  std::size_t nonzero_events = 0;
  std::size_t rejected_epochs = 0;
  std::vector<double> component_weights;
  bool component_added = false;
  bool vegas_initialized = false;
};

// Diagnostics returned after one held-out validation round
struct MFlowValidationReport {
  double ess_fraction = 0.0;
  double relative_variance = 0.0;
  double best_ess_fraction = 0.0;
  std::size_t events = 0;
  std::size_t nonzero_events = 0;
  std::size_t best_round = 0;
  std::size_t stale_rounds = 0;
  bool improved = false;
  bool early_stop = false;
};

// Monotonic rational-quadratic spline on the closed unit interval
class MSpline1D {
public:
  // Evaluate the forward spline using normalized positive knot parameters
  static MSplineResult Forward(double x, const std::vector<double> &widths,
                               const std::vector<double> &heights,
                               const std::vector<double> &derivatives);

  // Evaluate the analytic inverse spline using normalized positive knot
  // parameters
  static MSplineResult Inverse(double y, const std::vector<double> &widths,
                               const std::vector<double> &heights,
                               const std::vector<double> &derivatives);
};

// One masked coupling transform with an MLP-conditioned spline per active
// coordinate
class MSplineCoupling {
public:
  // One flat dense-layer layout inside the shared flow parameter vector
  struct DenseShape {
    std::size_t inputs = 0;
    std::size_t outputs = 0;
    std::size_t weight_offset = 0;
    std::size_t bias_offset = 0;
  };

  // Construct an unconfigured coupling transform
  MSplineCoupling() = default;

  // Compute the number of owned flat network parameters
  std::size_t ParameterCount() const { return parameter_count_; }

  // Internal immutable layout populated only by MSplineFlow::Configure
  std::size_t dimension_ = 0;
  std::size_t bins_ = 0;
  std::size_t parameter_offset_ = 0;
  std::size_t parameter_count_ = 0;
  double min_bin_ = 0.0;
  double min_derivative_ = 0.0;
  FlowActivation activation_ = FlowActivation::Silu;
  std::vector<std::size_t> condition_indices_;
  std::vector<std::size_t> transform_indices_;
  std::vector<std::size_t> permutation_;
  std::vector<std::size_t> inverse_permutation_;
  std::vector<DenseShape> dense_shapes_;
};

// Composition of alternating bounded spline-coupling transforms
class MSplineFlow {
public:
  // Initialize an identity flow with deterministic conditioner parameters
  void Configure(std::size_t dimension, const MNeuroJacConfig &config,
                 std::uint32_t initialization_seed);

  // Transform one uniform base point and return its exact output density
  MFlowSample Forward(const std::vector<double> &base) const;

  // Transform row-major base points with one vectorized evaluation
  MFlowSampleBatch ForwardBatch(const std::vector<double> &bases) const;

  // Evaluate the exact flow density at one unit-hypercube point
  double LogDensity(const std::vector<double> &point) const;

  // Evaluate exact flow densities for row-major unit-hypercube points
  std::vector<double> LogDensityBatch(const std::vector<double> &points) const;

  // Build immutable backend-specific metadata for repeated evaluation
  void PrepareInference();

  // Compute the phase-space dimension
  std::size_t Dimension() const { return dimension_; }

  // Compute the flat trainable parameters
  const std::vector<double> &Parameters() const { return parameters_; }

  // Compute one for matrix weights and zero for bias-like parameters
  const std::vector<unsigned char> &WeightDecayMask() const {
    return weight_decay_mask_;
  }

  // Replace the flat trainable parameters after validation
  void SetParameters(std::vector<double> parameters);

  // Apply one validated additive optimizer step and decoupled decay in place
  void ApplyParameterStep(const std::vector<double> &steps,
                          double matrix_decay_factor, double bias_decay_factor);

  // Initialize one separable spline map from an adapted VEGAS grid
  void
  InitializeVegasGrid(const std::vector<std::vector<double>> &quantile_edges);

  // Evaluate one mini-batch objective and its reverse-mode gradient
  double LossGradient(const std::vector<MFlowTrainingEvent> &events,
                      const std::vector<double> &coefficients,
                      const std::vector<std::size_t> &indices,
                      std::size_t population_size,
                      const MFlowObjective &objective, double uniform_mix,
                      std::vector<double> &gradient) const;

  // Access the immutable sequence of coupling transforms
  const std::vector<MSplineCoupling> &Couplings() const { return couplings_; }

private:
  std::size_t dimension_ = 0;
  std::size_t batch_size_ = 0;
  std::vector<MSplineCoupling> couplings_;
  std::vector<double> parameters_;
  std::vector<unsigned char> weight_decay_mask_;
  std::shared_ptr<const void> inference_state_;
};

// Exact convex superposition of independently parameterized spline flows
class MSplineFlowMixture {
public:
  // Configure independent spline components with uniform initial weights
  void Configure(std::size_t dimension, const MNeuroJacConfig &config,
                 std::uint32_t initialization_seed);

  // Compute the common phase space dimension
  std::size_t Dimension() const;

  // Compute the number of spline components
  std::size_t ComponentCount() const { return flows_.size(); }

  // Compute one immutable spline component
  const MSplineFlow &Component(std::size_t index) const;

  // Append one independent spline component with a small initial weight
  void AddComponent(const MNeuroJacConfig &config,
                    std::uint32_t initialization_seed, double weight);

  // Replace the last component weight while preserving all earlier ratios
  void SetLastWeight(double weight);

  // Compute normalized non-negative component weights
  std::vector<double> Weights() const;

  // Compute normalized logarithmic component weights without underflow
  std::vector<double> LogWeights() const;

  // Compute all flow parameters followed by independent weight logits
  std::vector<double> Parameters() const;

  // Replace all flow parameters and independent weight logits
  void SetParameters(const std::vector<double> &parameters);

  // Apply one AdamW step to all spline parameters and weight logits
  void ApplyParameterStep(const std::vector<double> &steps,
                          double matrix_decay_factor, double bias_decay_factor);

  // Apply one optimizer step to a single spline component
  void ApplyComponentStep(std::size_t component,
                          const std::vector<double> &steps,
                          double matrix_decay_factor,
                          double bias_decay_factor);

  // Initialize every spline component from one adapted VEGAS grid
  void
  InitializeVegasGrid(const std::vector<std::vector<double>> &quantile_edges);

  // Build immutable inference data for every component
  void PrepareInference();

  // Evaluate the exact spline mixture density at one point
  double LogDensity(const std::vector<double> &point) const;

  // Evaluate exact spline mixture densities for row-major points
  std::vector<double> LogDensityBatch(const std::vector<double> &points) const;

  // Evaluate one mini-batch objective and its reverse-mode gradient
  double LossGradient(const std::vector<MFlowTrainingEvent> &events,
                      const std::vector<double> &coefficients,
                      const std::vector<std::size_t> &indices,
                      std::size_t population_size,
                      const MFlowObjective &objective, double uniform_mix,
                      std::vector<double> &gradient) const;

private:
  std::vector<MSplineFlow> flows_;
  std::vector<double> weight_logits_;
};

// Adaptive event buffer, optimizer and immutable production proposal
class MFlowIntegrator {
public:
  // Construct an empty adaptive spline-flow integrator
  MFlowIntegrator();

  // Destroy the source-owned AdamW optimizer state
  ~MFlowIntegrator();

  // Prevent copies of adaptive optimizer and random-engine state
  MFlowIntegrator(const MFlowIntegrator &) = delete;
  // Prevent copy assignment of adaptive optimizer and random-engine state
  MFlowIntegrator &operator=(const MFlowIntegrator &) = delete;
  // Prevent moves while worker methods may hold a stable object address
  MFlowIntegrator(MFlowIntegrator &&) = delete;
  // Prevent move assignment while worker methods may hold a stable object
  // address
  MFlowIntegrator &operator=(MFlowIntegrator &&) = delete;

  // Read NUMERICS_NEUROJAC from one model NUMERICS file
  void ReadParameters(const std::string &filename);

  // Configure a fresh spline flow for one phase-space dimension
  void Configure(std::size_t dimension, std::uint32_t process_seed);

  // Initialize the configured flow from one frozen VEGAS grid
  void
  InitializeVegasGrid(const std::vector<std::vector<double>> &quantile_edges);

  // Serialize one configured frozen proposal into a backend-neutral JSON
  // document
  std::string SerializeModel() const;

  // Restore one backend-neutral JSON proposal for an expected phase-space
  // dimension
  void DeserializeModel(const std::string &document,
                        std::size_t expected_dimension);

  // Save one configured frozen proposal directly to disk
  void SaveModel(const std::string &filename) const;

  // Load one frozen proposal directly from disk
  void LoadModel(const std::string &filename, std::size_t expected_dimension);

  // Draw one exact mixed flow/uniform proposal point
  MFlowSample Sample(MRandom &random) const;

  // Draw exact mixed proposal points into contiguous row-major storage
  MFlowSampleBatch SampleBatch(MRandom &random, std::size_t count) const;

  // Evaluate the exact mixed proposal density
  double LogDensity(const std::vector<double> &point) const;

  // Compute the fresh target evaluations needed for the next training round
  std::size_t FreshEventCount() const;

  // Compute the held-out target evaluations needed for one validation round
  std::size_t ValidationEventCount() const;

  // Clear all events from the current adaptive round
  void ClearBuffer();

  // Clear all held-out events from the current validation round
  void ClearValidationBuffer();

  // Append one collection buffer after worker synchronization
  void BufferEvents(std::vector<MFlowTrainingEvent> events);

  // Append held-out events that are never used for optimizer updates or replay
  void BufferValidationEvents(std::vector<MFlowTrainingEvent> events);

  // Optimize one adaptive round using the selected objective
  MFlowTrainingReport TrainRound(std::size_t threads = 1);

  // Score the VEGAS-mapped proposal as the round-zero checkpoint
  MFlowValidationReport ValidateBaseline();

  // Score the updated proposal and retain the best held-out checkpoint
  MFlowValidationReport ValidateRound();

  // Freeze trained parameters for integration and event generation
  void Freeze();

  // Check whether the proposal is immutable and ready for production
  bool IsFrozen() const { return frozen_; }

  // Check whether the flow dimension and optimizer state are configured
  bool IsConfigured() const { return configured_; }

  // Access immutable numerical steering
  const MNeuroJacConfig &Config() const { return config_; }

  // Compute the first immutable spline flow for numerical regression tests
  const MSplineFlow &Flow() const { return mixture_.Component(0); }

  // Access the immutable exact spline mixture
  const MSplineFlowMixture &Mixture() const { return mixture_; }

private:
  // Optimize one epoch with the current proposal coefficients
  double TrainEpoch(std::vector<std::size_t> &order,
                    const std::vector<double> &boost_coefficients,
                    std::size_t epoch, std::size_t threads, MFlowTrainingReport &report);

  // Append a random replay subset to the fresh current-round events
  std::size_t AppendReplayEvents();

  // Retain fresh current-round events in the bounded replay pool
  void RetainFreshEvents(std::size_t fresh_events);

  // Score one buffered proposal with optional neural-round advancement
  MFlowValidationReport ValidateBufferedProposal(bool advance_round);

  MNeuroJacConfig config_;
  MSplineFlowMixture mixture_;
  std::vector<MFlowTrainingEvent> buffer_;
  std::vector<MFlowTrainingEvent> replay_buffer_;
  std::vector<MFlowTrainingEvent> validation_buffer_;
  std::vector<double> best_parameters_;
  std::vector<double> adamw_first_moment_;
  std::vector<double> adamw_second_moment_;
  double best_validation_ess_fraction_ = 0.0;
  std::uint32_t flow_seed_ = 0;
  std::uint64_t adamw_step_ = 0;
  std::size_t training_round_ = 0;
  std::size_t validation_round_ = 0;
  std::size_t validation_stale_rounds_ = 0;
  std::size_t best_validation_round_ = 0;
  std::size_t best_component_count_ = 0;
  bool vegas_initialized_ = false;
  bool configured_ = false;
  bool frozen_ = false;
  std::mt19937 optimizer_random_;
  std::mt19937 replay_random_;
};

// Compute the stable steering label for one activation function
std::string FlowActivationName(FlowActivation activation);

// Compute the stable steering label for one learning-rate schedule
std::string FlowLRScheduleName(FlowLRSchedule schedule);

// Compute the stable steering label for one coordinate permutation mode
std::string FlowPermutationName(FlowPermutation permutation);

} // namespace neurojac
} // namespace gra

#endif
