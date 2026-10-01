// NEUROJAC: Spline flow neural importance sampler
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// For similar constructions, see:
// [REFERENCE: Durkan et al., arxiv.org/abs/1906.04032]
// [REFERENCE: Heimel et al., arxiv.org/abs/2212.06172]


#include "Graniitti/Sampling/MNeuroJac.h"

#include <algorithm>
#include <cmath>
#include <compare>
#include <filesystem>
#include <fstream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"

#ifdef GRANIITTI_USE_LIBTORCH
#include <torch/torch.h>
#else
#include <Eigen/Dense>
#endif

using gra::aux::indices;

namespace gra::neurojac {

namespace {

using json = nlohmann::json;

// Compute the stable identifier for the backend-neutral model format
std::string NeuroJacModelFormat() { return "GRANIITTI_NEUROJAC"; }

// Compute the current backend-neutral model schema version
std::uint32_t NeuroJacModelVersion() { return 1; }

// Parse one activation steering label
FlowActivation ParseFlowActivation(const std::string &name) {
  if (name == "relu") {
    return FlowActivation::Relu;
  }
  if (name == "silu") {
    return FlowActivation::Silu;
  }
  if (name == "tanh") {
    return FlowActivation::Tanh;
  }
  throw std::invalid_argument(
      "NUMERICS_NEUROJAC::activation must be relu, silu or tanh");
}

// Parse one learning-rate schedule steering label
FlowLRSchedule ParseFlowLRSchedule(const std::string &name) {
  if (name == "none") {
    return FlowLRSchedule::None;
  }
  if (name == "cosine") {
    return FlowLRSchedule::Cosine;
  }
  if (name == "exponential") {
    return FlowLRSchedule::Exponential;
  }
  throw std::invalid_argument(
      "NUMERICS_NEUROJAC::lr_schedule must be none, cosine or exponential");
}

// Parse one coordinate permutation steering label
FlowPermutation ParseFlowPermutation(const std::string &name) {
  if (name == "none") {
    return FlowPermutation::None;
  }
  if (name == "interleave") {
    return FlowPermutation::Interleave;
  }
  if (name == "random") {
    return FlowPermutation::Random;
  }
  throw std::invalid_argument(
      "NUMERICS_NEUROJAC::permutation must be none, interleave or random");
}

// Parse one adaptive objective steering label
FlowLoss ParseFlowLoss(const std::string &name) {
  if (name == "forward_kl") {
    return FlowLoss::ForwardKL;
  }
  if (name == "chi2") {
    return FlowLoss::Chi2;
  }
  if (name == "renyi") {
    return FlowLoss::Renyi;
  }
  throw std::invalid_argument(
      "NUMERICS_NEUROJAC::loss name must be forward_kl, chi2 or renyi");
}

// Parse one adaptive objective configuration
MFlowLossConfig ParseFlowLossConfig(const json &block) {
  MFlowLossConfig parsed;
  parsed.kind = ParseFlowLoss(block.at("name").get<std::string>());
  const std::vector<double> alpha =
      block.at("alpha").get<std::vector<double>>();
  if (alpha.size() != 2) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC::loss alpha must contain two values");
  }
  parsed.alpha_initial = alpha[0];
  parsed.alpha_final = alpha[1];
  parsed.annealing_fraction = block.at("anneal_frac");
  return parsed;
}

// Serialize one adaptive objective configuration
json SerializeFlowLossConfig(const MFlowLossConfig &config) {
  return {{"name", FlowLossName(config.kind)},
          {"alpha", {config.alpha_initial, config.alpha_final}},
          {"anneal_frac", config.annealing_fraction}};
}

// Parse a nonnegative integral count before unsigned conversion
std::size_t ParseFlowCount(const json &value, const std::string &name) {
  if (!value.is_number_integer() || (!value.is_number_unsigned() && value.get<std::int64_t>() < 0) ||
      value.get<std::uint64_t>() > std::numeric_limits<std::size_t>::max()) {
    throw std::invalid_argument("NEUROJAC: " + name + " must be a nonnegative integer count");
  }
  return value.get<std::size_t>();
}

// Parse the complete numerical flow configuration
MNeuroJacConfig ParseFlowConfig(const json &block) {
  MNeuroJacConfig parsed;
  parsed.flow_components = ParseFlowCount(block.at("flow_components"), "flow_components");
  parsed.layers          = ParseFlowCount(block.at("layers"), "layers");
  parsed.bins            = ParseFlowCount(block.at("bins"), "bins");
  const auto &hidden     = block.at("hidden");
  if (!hidden.is_array()) { throw std::invalid_argument("NUMERICS_NEUROJAC: hidden must be an array"); }
  parsed.hidden.clear();
  for (const auto &width : hidden) { parsed.hidden.push_back(ParseFlowCount(width, "hidden")); }
  parsed.activation =
      ParseFlowActivation(block.at("activation").get<std::string>());
  parsed.permutation =
      ParseFlowPermutation(block.at("permutation").get<std::string>());
  parsed.buffer_size          = ParseFlowCount(block.at("buffer_size"), "buffer_size");
  parsed.batch_size           = ParseFlowCount(block.at("batch_size"), "batch_size");
  parsed.replay_capacity      = ParseFlowCount(block.at("replay_cap"), "replay_cap");
  parsed.validation_size      = ParseFlowCount(block.at("val_size"), "val_size");
  parsed.validation_patience  = ParseFlowCount(block.at("val_patience"), "val_patience");
  parsed.rounds               = ParseFlowCount(block.at("rounds"), "rounds");
  parsed.epochs               = ParseFlowCount(block.at("epochs"), "epochs");
  parsed.vegas_init = block.at("vegas_init").get<bool>();
  parsed.vegas_spline_init = block.at("vegas_spline_init").get<bool>();
  parsed.vegas_ncall          = ParseFlowCount(block.at("vegas_ncall"), "vegas_ncall");
  parsed.vegas_rounds         = ParseFlowCount(block.at("vegas_rounds"), "vegas_rounds");
  parsed.vegas_sampler_rounds = ParseFlowCount(block.at("vegas_sampler_rounds"), "vegas_sampler_rounds");
  parsed.learning_rate = block.at("learning_rate");
  parsed.final_learning_rate = block.at("final_learning_rate");
  parsed.validation_min_delta = block.at("val_min_delta");
  parsed.lr_schedule =
      ParseFlowLRSchedule(block.at("lr_schedule").get<std::string>());
  parsed.uniform_mix = block.at("uniform_mix");
  parsed.replay_fraction = block.at("replay_frac");
  parsed.min_bin = block.at("min_bin");
  parsed.min_derivative = block.at("min_derivative");
  parsed.density_tolerance  = block.at("density_tol");
  parsed.gradient_clip = block.at("gradient_clip");
  parsed.adamw_beta1 = block.at("adamw_beta1");
  parsed.adamw_beta2 = block.at("adamw_beta2");
  parsed.adamw_epsilon = block.at("adamw_epsilon");
  parsed.adamw_weight_decay = block.at("adamw_weight_decay");
  parsed.adamw_decay_biases = block.at("adamw_decay_biases");
  parsed.loss = ParseFlowLossConfig(block.at("loss"));
  return parsed;
}

// Serialize the complete numerical flow configuration
json SerializeFlowConfig(const MNeuroJacConfig &config) {
  return {{"flow_components", config.flow_components},
          {"layers", config.layers},
          {"bins", config.bins},
          {"hidden", config.hidden},
          {"activation", FlowActivationName(config.activation)},
          {"permutation", FlowPermutationName(config.permutation)},
          {"buffer_size", config.buffer_size},
          {"batch_size", config.batch_size},
          {"replay_cap", config.replay_capacity},
          {"val_size", config.validation_size},
          {"val_patience", config.validation_patience},
          {"rounds", config.rounds},
          {"epochs", config.epochs},
          {"vegas_init", config.vegas_init},
          {"vegas_spline_init", config.vegas_spline_init},
          {"vegas_ncall", config.vegas_ncall},
          {"vegas_rounds", config.vegas_rounds},
          {"vegas_sampler_rounds", config.vegas_sampler_rounds},
          {"learning_rate", config.learning_rate},
          {"final_learning_rate", config.final_learning_rate},
          {"val_min_delta", config.validation_min_delta},
          {"lr_schedule", FlowLRScheduleName(config.lr_schedule)},
          {"uniform_mix", config.uniform_mix},
          {"replay_frac", config.replay_fraction},
          {"min_bin", config.min_bin},
          {"min_derivative", config.min_derivative},
          {"density_tol", config.density_tolerance},
          {"gradient_clip", config.gradient_clip},
          {"adamw_beta1", config.adamw_beta1},
          {"adamw_beta2", config.adamw_beta2},
          {"adamw_epsilon", config.adamw_epsilon},
          {"adamw_weight_decay", config.adamw_weight_decay},
          {"adamw_decay_biases", config.adamw_decay_biases},
          {"loss", SerializeFlowLossConfig(config.loss)}};
}

#ifndef GRANIITTI_USE_LIBTORCH

// Instance-owned reverse-mode tape with constant-time node insertion
class MReverseTape {
public:
  // Compute the sentinel used for nodes without a parent
  static constexpr std::size_t NoParent() {
    return std::numeric_limits<std::size_t>::max();
  }

  // One scalar operation and its local derivatives
  struct Node {
    double value = 0.0;
    double adjoint = 0.0;
    std::size_t parent_a = NoParent();
    std::size_t parent_b = NoParent();
    double derivative_a = 0.0;
    double derivative_b = 0.0;
  };

  // Insert an independent parameter or constant node
  std::size_t Leaf(double value) {
    nodes_.push_back({value, 0.0, NoParent(), NoParent(), 0.0, 0.0});
    return nodes_.size() - 1;
  }

  // Insert one unary operation with its evaluated local derivative
  std::size_t Unary(double value, std::size_t parent, double derivative) {
    nodes_.push_back({value, 0.0, parent, NoParent(), derivative, 0.0});
    return nodes_.size() - 1;
  }

  // Insert one binary operation with its evaluated local derivatives
  std::size_t Binary(double value, std::size_t parent_a, std::size_t parent_b,
                     double derivative_a, double derivative_b) {
    nodes_.push_back(
        {value, 0.0, parent_a, parent_b, derivative_a, derivative_b});
    return nodes_.size() - 1;
  }

  // Reuse the allocated node storage for one independent spline derivative
  void Reset() { nodes_.clear(); }

  // Reserve node storage before repeated small reverse sweeps
  void Reserve(std::size_t count) { nodes_.reserve(count); }

  // Propagate one scalar objective back through the tape
  void Backward(std::size_t root) {
    if (root >= nodes_.size()) {
      throw std::out_of_range("MReverseTape::Backward: invalid root");
    }
    for (Node &node : nodes_) {
      node.adjoint = 0.0;
    }
    nodes_[root].adjoint = 1.0;
    for (std::size_t i = root + 1; i-- > 0;) {
      const Node &node = nodes_[i];
      if (node.parent_a != NoParent()) {
        nodes_[node.parent_a].adjoint += node.adjoint * node.derivative_a;
      }
      if (node.parent_b != NoParent()) {
        nodes_[node.parent_b].adjoint += node.adjoint * node.derivative_b;
      }
    }
  }

  // Compute the forward value stored at one tape node
  double Value(std::size_t index) const { return nodes_[index].value; }

  // Compute the reverse derivative stored at one tape node
  double Adjoint(std::size_t index) const { return nodes_[index].adjoint; }

private:
  std::vector<Node> nodes_;
};

// Lightweight scalar handle into one instance-owned reverse tape
class MReverseScalar {
public:
  // Construct an empty handle for vector allocation before assignment
  MReverseScalar() = default;

  // Construct a scalar handle from one tape node
  MReverseScalar(MReverseTape &tape, std::size_t node)
      : tape_(&tape), node_(node) {}

  // Construct one independent trainable parameter
  static MReverseScalar Parameter(MReverseTape &tape, double value) {
    return MReverseScalar(tape, tape.Leaf(value));
  }

  // Construct one constant scalar on a tape
  static MReverseScalar Constant(MReverseTape &tape, double value) {
    return MReverseScalar(tape, tape.Leaf(value));
  }

  // Compute the owning tape after validating the handle
  MReverseTape &Tape() const {
    if (tape_ == nullptr) {
      throw std::logic_error("MReverseScalar::Tape: empty scalar handle");
    }
    return *tape_;
  }

  // Compute the node index on the owning tape
  std::size_t NodeIndex() const { return node_; }

  // Compute the scalar forward value
  double Value() const { return Tape().Value(node_); }

private:
  MReverseTape *tape_ = nullptr;
  std::size_t node_ = MReverseTape::NoParent();
};

// Validate that two differentiable operands belong to the same tape
MReverseTape &CommonTape(const MReverseScalar &left,
                         const MReverseScalar &right) {
  if (&left.Tape() != &right.Tape()) {
    throw std::invalid_argument(
        "MReverseScalar: operands belong to different tapes");
  }
  return left.Tape();
}

// Add two differentiable scalar values
MReverseScalar operator+(const MReverseScalar &left,
                         const MReverseScalar &right) {
  MReverseTape &tape = CommonTape(left, right);
  return MReverseScalar(tape, tape.Binary(left.Value() + right.Value(),
                                          left.NodeIndex(), right.NodeIndex(),
                                          1.0, 1.0));
}

// Subtract two differentiable scalar values
MReverseScalar operator-(const MReverseScalar &left,
                         const MReverseScalar &right) {
  MReverseTape &tape = CommonTape(left, right);
  return MReverseScalar(tape, tape.Binary(left.Value() - right.Value(),
                                          left.NodeIndex(), right.NodeIndex(),
                                          1.0, -1.0));
}

// Multiply two differentiable scalar values
MReverseScalar operator*(const MReverseScalar &left,
                         const MReverseScalar &right) {
  MReverseTape &tape = CommonTape(left, right);
  return MReverseScalar(tape, tape.Binary(left.Value() * right.Value(),
                                          left.NodeIndex(), right.NodeIndex(),
                                          right.Value(), left.Value()));
}

// Divide two differentiable scalar values
MReverseScalar operator/(const MReverseScalar &left,
                         const MReverseScalar &right) {
  MReverseTape &tape = CommonTape(left, right);
  const double inverse = 1.0 / right.Value();
  return MReverseScalar(tape,
                        tape.Binary(left.Value() * inverse, left.NodeIndex(),
                                    right.NodeIndex(), inverse,
                                    -left.Value() * inverse * inverse));
}

// Add a constant to a differentiable scalar
MReverseScalar operator+(const MReverseScalar &left, double right) {
  return MReverseScalar(left.Tape(), left.Tape().Unary(left.Value() + right,
                                                       left.NodeIndex(), 1.0));
}

// Add a differentiable scalar to a constant
MReverseScalar operator+(double left, const MReverseScalar &right) {
  return right + left;
}

// Subtract a constant from a differentiable scalar
MReverseScalar operator-(const MReverseScalar &left, double right) {
  return MReverseScalar(left.Tape(), left.Tape().Unary(left.Value() - right,
                                                       left.NodeIndex(), 1.0));
}

// Subtract a differentiable scalar from a constant
MReverseScalar operator-(double left, const MReverseScalar &right) {
  return MReverseScalar(
      right.Tape(),
      right.Tape().Unary(left - right.Value(), right.NodeIndex(), -1.0));
}

// Multiply a differentiable scalar by a constant
MReverseScalar operator*(const MReverseScalar &left, double right) {
  return MReverseScalar(
      left.Tape(),
      left.Tape().Unary(left.Value() * right, left.NodeIndex(), right));
}

// Multiply a constant by a differentiable scalar
MReverseScalar operator*(double left, const MReverseScalar &right) {
  return right * left;
}

// Negate one differentiable scalar
MReverseScalar operator-(const MReverseScalar &value) { return -1.0 * value; }

// Accumulate a differentiable scalar in place
MReverseScalar &operator+=(MReverseScalar &left, const MReverseScalar &right) {
  left = left + right;
  return left;
}

#endif

// Compute the fixed tolerance for closed-interval spline boundary checks
constexpr double SplineTolerance() { return 1e-10; }

// Compute the relative roundoff multiplier for the inverse quadratic
constexpr double SplineInverseRoundoffFactor() { return 64.0; }

// Compute a common scale for the inverse quadratic coefficients
double SplineCoefficientScale(double a, double b, double c) {
  const double scale = std::max({std::abs(a), std::abs(b), std::abs(c)});
  if (!(scale > 0.0) || !std::isfinite(scale)) {
    throw std::runtime_error(
        "SplineCoefficientScale: invalid inverse coefficients");
  }
  return scale;
}

// Compute whether the inverse quadratic is numerically linear
bool SplineUsesLinearRoot(double a, double b, double c) {
  const double scale = SplineCoefficientScale(a, b, c);
  const double tolerance = SplineInverseRoundoffFactor() *
                           std::numeric_limits<double>::epsilon() * scale;
  return std::abs(a) <= tolerance;
}

// Validate one discriminant formed from normalized inverse coefficients
double ValidatedSplineDiscriminant(double normalized_a, double normalized_b,
                                   double normalized_c,
                                   const std::string &context) {
  const double product = 4.0 * normalized_a * normalized_c;
  const double discriminant = normalized_b * normalized_b - product;
  const double scale = normalized_b * normalized_b + std::abs(product);
  const double tolerance = SplineInverseRoundoffFactor() *
                           std::numeric_limits<double>::epsilon() * scale;
  if (discriminant < -tolerance) {
    throw std::runtime_error(context + ": negative quadratic discriminant");
  }
  return std::max(discriminant, 0.0);
}

// Compute the numerical value used for piecewise spline branch selection
double ScalarValue(double value) { return value; }

#ifndef GRANIITTI_USE_LIBTORCH
// Compute the numerical value used for piecewise spline branch selection
double ScalarValue(const MReverseScalar &value) { return value.Value(); }
#endif

// Evaluate an exponential for a scalar double
double ScalarExp(double value) { return std::exp(value); }

#ifndef GRANIITTI_USE_LIBTORCH
// Evaluate an exponential while preserving the reverse-mode graph
MReverseScalar ScalarExp(const MReverseScalar &value) {
  const double result = std::exp(value.Value());
  return MReverseScalar(value.Tape(),
                        value.Tape().Unary(result, value.NodeIndex(), result));
}
#endif

// Evaluate a logarithm for a scalar double
double ScalarLog(double value) { return std::log(value); }

#ifndef GRANIITTI_USE_LIBTORCH
// Evaluate a logarithm while preserving the reverse-mode graph
MReverseScalar ScalarLog(const MReverseScalar &value) {
  return MReverseScalar(
      value.Tape(), value.Tape().Unary(std::log(value.Value()),
                                       value.NodeIndex(), 1.0 / value.Value()));
}
#endif

// Evaluate a square root for a scalar double
double ScalarSqrt(double value) { return std::sqrt(value); }

#ifndef GRANIITTI_USE_LIBTORCH
// Evaluate a square root while preserving the reverse-mode graph
MReverseScalar ScalarSqrt(const MReverseScalar &value) {
  const double result = std::sqrt(value.Value());
  return MReverseScalar(
      value.Tape(),
      value.Tape().Unary(result, value.NodeIndex(), 0.5 / result));
}
#endif

// Evaluate a stable hyperbolic tangent for a scalar double
double ScalarTanh(double value) { return std::tanh(value); }

// Evaluate a stable logistic sigmoid for a scalar double
double ScalarSigmoid(double value) {
  return value >= 0.0 ? 1.0 / (1.0 + std::exp(-value))
                      : std::exp(value) / (1.0 + std::exp(value));
}

// Apply one configured hidden-layer activation
double ScalarActivation(double value, FlowActivation activation) {
  switch (activation) {
  case FlowActivation::Relu:
    return std::max(value, 0.0);
  case FlowActivation::Silu:
    return value * ScalarSigmoid(value);
  case FlowActivation::Tanh:
    return ScalarTanh(value);
  }
  throw std::invalid_argument("ScalarActivation: unknown activation function");
}

// Evaluate a stable softplus for a scalar double
double ScalarSoftplus(double value) {
  return value > 0.0 ? value + std::log1p(std::exp(-value))
                     : std::log1p(std::exp(value));
}

#ifndef GRANIITTI_USE_LIBTORCH
// Evaluate a stable softplus while preserving the reverse-mode graph
MReverseScalar ScalarSoftplus(const MReverseScalar &value) {
  const double input = value.Value();
  const double result = ScalarSoftplus(input);
  const double derivative = input >= 0.0
                                ? 1.0 / (1.0 + std::exp(-input))
                                : std::exp(input) / (1.0 + std::exp(input));
  return MReverseScalar(
      value.Tape(), value.Tape().Unary(result, value.NodeIndex(), derivative));
}
#endif

// Construct a scalar constant matching one double reference
double ScalarConstant(double, double value) { return value; }

#ifndef GRANIITTI_USE_LIBTORCH
// Construct a scalar constant on the same reverse tape as a reference
MReverseScalar ScalarConstant(const MReverseScalar &reference, double value) {
  return MReverseScalar::Constant(reference.Tape(), value);
}

#endif

// Convert unconstrained logits into positive entries with a prescribed sum
template <typename Scalar>
std::vector<Scalar> NormalizeBins(const std::vector<Scalar> &raw,
                                  double minimum) {
  if (raw.empty()) {
    throw std::invalid_argument("NormalizeBins: empty bin vector");
  }
  const double available = 1.0 - minimum * static_cast<double>(raw.size());
  if (!(available > 0.0)) {
    throw std::invalid_argument(
        "NormalizeBins: minimum bin size leaves no available interval");
  }

  double maximum = ScalarValue(raw.front());
  for (const Scalar &entry : raw) {
    maximum = std::max(maximum, ScalarValue(entry));
  }

  std::vector<Scalar> probabilities(raw.size());
  Scalar denominator = ScalarConstant(raw.front(), 0.0);
  for (const auto &i : indices(raw)) {
    probabilities[i] = ScalarExp(raw[i] - maximum);
    denominator += probabilities[i];
  }
  for (Scalar &entry : probabilities) {
    entry = minimum + available * entry / denominator;
  }
  return probabilities;
}

// Convert unconstrained logits into normalized bins using caller-owned storage
template <typename Scalar>
void NormalizeBins(const std::vector<Scalar> &raw, double minimum,
                   std::vector<Scalar> &probabilities) {
  if (raw.empty()) {
    throw std::invalid_argument("NormalizeBins: empty bin vector");
  }
  const double available = 1.0 - minimum * static_cast<double>(raw.size());
  if (!(available > 0.0)) {
    throw std::invalid_argument(
        "NormalizeBins: minimum bin size leaves no available interval");
  }

  double maximum = ScalarValue(raw.front());
  for (const Scalar &entry : raw) {
    maximum = std::max(maximum, ScalarValue(entry));
  }
  probabilities.resize(raw.size());
  Scalar denominator = ScalarConstant(raw.front(), 0.0);
  for (const auto &i : indices(raw)) {
    probabilities[i] = ScalarExp(raw[i] - maximum);
    denominator += probabilities[i];
  }
  for (Scalar &entry : probabilities) {
    entry = minimum + available * entry / denominator;
  }
}

// Convert unconstrained derivative logits into strictly positive knot slopes
template <typename Scalar>
std::vector<Scalar> NormalizeDerivatives(const std::vector<Scalar> &raw,
                                         double minimum) {
  std::vector<Scalar> derivatives(raw.size());
  for (const auto &i : indices(raw)) {
    derivatives[i] = minimum + ScalarSoftplus(raw[i]);
  }
  return derivatives;
}

// Evaluate log((1-epsilon) q_flow + epsilon) without numerical overflow
template <typename Scalar>
Scalar MixedLogDensity(const Scalar &flow_log_density, double uniform_mix) {
  const Scalar flow_term = std::log1p(-uniform_mix) + flow_log_density;
  const double flat_term = std::log(uniform_mix);
  const double maximum = std::max(ScalarValue(flow_term), flat_term);
  return maximum + ScalarLog(ScalarExp(flow_term - maximum) +
                             std::exp(flat_term - maximum));
}

// Draw one spline component from normalized mixture weights
std::size_t DrawFlowComponent(const std::vector<double> &weights,
                              MRandom &random) {
  if (weights.size() == 1) {
    return 0;
  }
  const double draw = random.U(0.0, 1.0);
  double cumulative = 0.0;
  for (std::size_t component = 0; component + 1 < weights.size(); ++component) {
    cumulative += weights[component];
    if (draw < cumulative) {
      return component;
    }
  }
  return weights.size() - 1;
}

// Copy buffered phase space points into contiguous row major storage
std::vector<double>
FlattenFlowEventPoints(const std::vector<MFlowTrainingEvent> &events,
                       std::size_t dimension) {
  std::vector<double> points;
  points.reserve(events.size() * dimension);
  for (const MFlowTrainingEvent &event : events) {
    points.insert(points.end(), event.point.begin(), event.point.end());
  }
  return points;
}

// Evaluate independent density blocks without changing the event reduction order
template <typename Proposal>
std::vector<double> FlowLogDensities(const Proposal &flow, const std::vector<double> &points,
                                    std::size_t threads, std::size_t batch_size) {
  const std::size_t count = points.size() / flow.Dimension();
  if (count <= batch_size) { return flow.LogDensityBatch(points); }
  const std::size_t blocks = 1 + (count - 1) / batch_size;
  threads = std::clamp<std::size_t>(threads, 1, blocks);
  std::vector<double> densities(count);
  const auto evaluate = [&](std::size_t thread) {
    for (std::size_t block = thread; block < blocks; block += threads) {
      const std::size_t begin = block * batch_size;
      const std::size_t end = std::min(begin + batch_size, count);
      const std::vector<double> points_block(points.begin() + begin * flow.Dimension(),
                                              points.begin() + end * flow.Dimension());
      const auto values = flow.LogDensityBatch(points_block);
      std::copy(values.begin(), values.end(), densities.begin() + begin);
    }
  };
  if (threads == 1) { evaluate(0); } else { math::ParallelFor(threads, evaluate); }
  return densities;
}

// Evaluate mixed proposal log densities needed by one objective
template <typename Proposal>
std::vector<double> ObjectiveMixedLogDensities(
    const Proposal &flow, const std::vector<MFlowTrainingEvent> &events, const std::vector<double> &points,
    const MFlowObjective &objective, double uniform_mix, std::size_t threads, std::size_t batch_size) {
  std::vector<double> mixed_log_densities(events.size(), 0.0);
  if (!objective.UsesProposalDensity()) {
    return mixed_log_densities;
  }
  const std::vector<double> flow_log_densities = FlowLogDensities(flow, points, threads, batch_size);
  for (const auto &event : indices(events)) {
    mixed_log_densities[event] =
        MixedLogDensity(flow_log_densities[event], uniform_mix);
  }
  return mixed_log_densities;
}

// Build globally normalized event coefficients for one objective
template <typename Proposal>
std::vector<double> NormalizedObjectiveCoefficients(
    const Proposal &flow, const std::vector<MFlowTrainingEvent> &events, const std::vector<double> &points,
    const MFlowObjective &objective, double uniform_mix, std::size_t threads, std::size_t batch_size) {
  const std::vector<double> mixed_log_densities =
      ObjectiveMixedLogDensities(flow, events, points, objective, uniform_mix, threads, batch_size);
  std::vector<double> log_coefficients(
      events.size(), -std::numeric_limits<double>::infinity());
  double maximum_log_coefficient = -std::numeric_limits<double>::infinity();
  for (const auto &event : indices(events)) {
    const double absolute_target = std::abs(events[event].target);
    if (std::fpclassify(absolute_target) == FP_ZERO) {
      continue;
    }
    log_coefficients[event] = objective.LogCoefficient(
        std::log(absolute_target), events[event].collection_log_density,
        mixed_log_densities[event]);
    maximum_log_coefficient =
        std::max(maximum_log_coefficient, log_coefficients[event]);
  }
  if (!std::isfinite(maximum_log_coefficient)) {
    throw std::runtime_error(
        "NormalizedObjectiveCoefficients: all buffered targets are zero");
  }

  std::vector<double> coefficients(events.size(), 0.0);
  double coefficient_sum = 0.0;
  for (const auto &event : indices(events)) {
    if (!std::isfinite(log_coefficients[event])) {
      continue;
    }
    coefficients[event] =
        std::exp(log_coefficients[event] - maximum_log_coefficient);
    coefficient_sum += coefficients[event];
  }
  for (double &coefficient : coefficients) {
    coefficient /= coefficient_sum;
  }
  return coefficients;
}

// Evaluate the frozen old mixture and the newest component separately
std::pair<std::vector<double>, std::vector<double>> BoostLogDensities(
    const MSplineFlowMixture &mixture,
    const std::vector<MFlowTrainingEvent> &events) {
  const std::vector<double> points =
      FlattenFlowEventPoints(events, mixture.Dimension());
  const std::size_t last = mixture.ComponentCount() - 1;
  const std::vector<double> log_weights = mixture.LogWeights();
  const double maximum_weight =
      *std::max_element(log_weights.begin(), log_weights.end() - 1);
  double old_sum = 0.0;
  for (std::size_t component = 0; component < last; ++component) {
    old_sum += std::exp(log_weights[component] - maximum_weight);
  }
  const double log_sum = std::log(old_sum);
  std::vector<double> old_log(events.size(),
                              -std::numeric_limits<double>::infinity());
  for (std::size_t component = 0; component < last; ++component) {
    const std::vector<double> component_log =
        mixture.Component(component).LogDensityBatch(points);
    const double log_weight =
        (log_weights[component] - maximum_weight) - log_sum;
    for (const auto &event : indices(events)) {
      const double term = log_weight + component_log[event];
      const double maximum = std::max(old_log[event], term);
      old_log[event] =
          maximum + std::log(std::exp(old_log[event] - maximum) +
                             std::exp(term - maximum));
    }
  }
  return {std::move(old_log),
          mixture.Component(last).LogDensityBatch(points)};
}

// Compute normalized chi-squared mixture derivatives for the newest flow
std::vector<double> BoostChi2Coefficients(
    const MSplineFlowMixture &mixture,
    const std::vector<MFlowTrainingEvent> &events, const std::vector<double> &points,
    const std::vector<double> &coefficients, double uniform_mix, std::size_t threads, std::size_t batch_size) {
  const std::vector<double> mixture_log = FlowLogDensities(mixture, points, threads, batch_size);
  const std::size_t last = mixture.ComponentCount() - 1;
  const std::vector<double> component_log =
      FlowLogDensities(mixture.Component(last), points, threads, batch_size);
  const double log_weight = mixture.LogWeights().back();
  std::vector<double> output(events.size());
  for (const auto &event : indices(events)) {
    output[event] = std::log(coefficients[event]) + std::log1p(-uniform_mix) +
                    log_weight + component_log[event] -
                    2.0 * MixedLogDensity(mixture_log[event], uniform_mix);
  }
  // Preserve the gradient direction without falling below the AdamW floor
  const double maximum = *std::max_element(output.begin(), output.end());
  for (double &coefficient : output) {
    coefficient = std::exp(coefficient - maximum);
  }
  const double sum = gra::Sum(output);
  for (double &coefficient : output) {
    coefficient /= sum;
  }
  return output;
}

// Minimize the convex sampled second moment using its signed derivative
double OptimizeBoostWeight(const std::vector<MFlowTrainingEvent> &events,
                           const std::vector<double> &old_log,
                           const std::vector<double> &new_log,
                           double uniform_mix) {
  const auto derivative = [&](double weight) {
    double maximum = -std::numeric_limits<double>::infinity();
    long double sum = 0.0L;
    for (const auto &event : indices(events)) {
      const double target = std::abs(events[event].target);
      if (std::fpclassify(target) == FP_ZERO) { continue; }
      const double flow_log = std::max(std::log1p(-weight) + old_log[event],
                                       std::log(weight) + new_log[event]);
      const double flow_sum = std::exp(std::log1p(-weight) + old_log[event] - flow_log) +
                              std::exp(std::log(weight) + new_log[event] - flow_log);
      const double mixed_log = MixedLogDensity(flow_log + std::log(flow_sum), uniform_mix);
      const double term = 2.0 * std::log(target) - events[event].collection_log_density - mixed_log;
      const long double scale = static_cast<long double>(std::log1p(-uniform_mix)) - mixed_log;
      const long double slope = std::exp(scale + old_log[event]) - std::exp(scale + new_log[event]);
      if (term > maximum) {
        sum = sum * std::exp(static_cast<long double>(maximum) - term) + slope;
        maximum = term;
      } else {
        sum += std::exp(static_cast<long double>(term) - maximum) * slope;
      }
    }
    return sum;
  };

  // Keep both mixture fractions resolvable and stop at adjacent floating-point weights
  double lower = std::numeric_limits<double>::epsilon();
  double upper = 1.0 - lower;
  if (derivative(lower) >= 0.0L) { return lower; }
  if (derivative(upper) <= 0.0L) { return upper; }
  double middle = std::midpoint(lower, upper);
  while (lower < middle && middle < upper) {
    if (derivative(middle) > 0.0L) { upper = middle; } else { lower = middle; }
    middle = std::midpoint(lower, upper);
  }
  return middle;
}

// Check that one phase-space point belongs to the closed unit hypercube
void ValidateUnitPoint(const std::vector<double> &point, std::size_t dimension,
                       const std::string &context) {
  if (point.size() != dimension) {
    throw std::invalid_argument(context + ": point dimension mismatch");
  }
  for (double value : point) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
      throw std::invalid_argument(context +
                                  ": point is outside the unit hypercube");
    }
  }
}

// Validate flat row-major points and return their event count
std::size_t ValidateFlatUnitPoints(const std::vector<double> &points,
                                   std::size_t dimension,
                                   const std::string &context) {
  if (dimension == 0 || points.size() % dimension != 0) {
    throw std::invalid_argument(context + ": flat point dimension mismatch");
  }
  for (double value : points) {
    if (!std::isfinite(value) || value < 0.0 || value > 1.0) {
      throw std::invalid_argument(context +
                                  ": point is outside the unit hypercube");
    }
  }
  return points.size() / dimension;
}

// Validate target events and their collection proposal densities
void ValidateFlowEvents(const std::vector<MFlowTrainingEvent> &events,
                        std::size_t dimension, const std::string &context) {
  for (const MFlowTrainingEvent &event : events) {
    ValidateUnitPoint(event.point, dimension, context);
    if (!std::isfinite(event.target) ||
        !std::isfinite(event.collection_log_density)) {
      throw std::invalid_argument(context + ": non-finite target or density");
    }
  }
}

// Check normalized public spline knot parameters
void ValidateSplineParameters(const std::vector<double> &widths,
                              const std::vector<double> &heights,
                              const std::vector<double> &derivatives) {
  if (widths.empty() || heights.size() != widths.size() ||
      derivatives.size() != widths.size() + 1) {
    throw std::invalid_argument("MSpline1D: inconsistent knot parameter sizes");
  }
  const auto positive_finite = [](double value) {
    return std::isfinite(value) && value > 0.0;
  };
  if (!std::all_of(widths.begin(), widths.end(), positive_finite) ||
      !std::all_of(heights.begin(), heights.end(), positive_finite) ||
      !std::all_of(derivatives.begin(), derivatives.end(), positive_finite)) {
    throw std::invalid_argument(
        "MSpline1D: knot parameters must be finite and positive");
  }
  const double width_sum = gra::Sum(widths);
  const double height_sum = gra::Sum(heights);
  if (std::abs(width_sum - 1.0) > SplineTolerance() ||
      std::abs(height_sum - 1.0) > SplineTolerance()) {
    throw std::invalid_argument(
        "MSpline1D: widths and heights must each sum to one");
  }
}

// Locate a spline bin from a coordinate and differentiable cumulative knots
template <typename Scalar>
std::size_t FindBin(const Scalar &coordinate, const std::vector<Scalar> &bins) {
  double cumulative = 0.0;
  for (std::size_t i = 0; i + 1 < bins.size(); ++i) {
    cumulative += ScalarValue(bins[i]);
    if (ScalarValue(coordinate) < cumulative) {
      return i;
    }
  }
  return bins.size() - 1;
}

// Evaluate one normalized rational-quadratic spline in the forward direction
// y = y_k + h[delta xi^2+d_k xi(1-xi)]/[delta+(d_{k+1}+d_k-2delta)xi(1-xi)]
template <typename Scalar>
std::pair<Scalar, Scalar>
SplineForward(const Scalar &x, const std::vector<Scalar> &widths,
              const std::vector<Scalar> &heights,
              const std::vector<Scalar> &derivatives) {
  const double x_value = ScalarValue(x);
  if (x_value < -SplineTolerance() || x_value > 1.0 + SplineTolerance()) {
    throw std::invalid_argument("SplineForward: coordinate is outside [0,1]");
  }

  const std::size_t bin = FindBin(x, widths);
  Scalar x_k = ScalarConstant(x, 0.0);
  Scalar y_k = ScalarConstant(x, 0.0);
  for (std::size_t i = 0; i < bin; ++i) {
    x_k += widths[i];
    y_k += heights[i];
  }

  const Scalar width = widths[bin];
  const Scalar height = heights[bin];
  const Scalar delta = height / width;
  const Scalar xi = (x - x_k) / width;
  const Scalar one_minus_xi = 1.0 - xi;
  const Scalar xi_product = xi * one_minus_xi;
  const Scalar denominator =
      delta +
      (derivatives[bin + 1] + derivatives[bin] - 2.0 * delta) * xi_product;
  const Scalar numerator =
      height * (delta * xi * xi + derivatives[bin] * xi_product);
  const Scalar derivative_numerator =
      delta * delta *
      (derivatives[bin + 1] * xi * xi + 2.0 * delta * xi_product +
       derivatives[bin] * one_minus_xi * one_minus_xi);
  const Scalar jacobian = derivative_numerator / (denominator * denominator);
  return {y_k + numerator / denominator, ScalarLog(jacobian)};
}

// Evaluate one normalized rational-quadratic spline in the inverse direction
// Solve a xi^2+b xi+c = 0 in the selected monotonic bin
template <typename Scalar>
std::pair<Scalar, Scalar>
SplineInverse(const Scalar &y, const std::vector<Scalar> &widths,
              const std::vector<Scalar> &heights,
              const std::vector<Scalar> &derivatives) {
  const double y_value = ScalarValue(y);
  if (y_value < -SplineTolerance() || y_value > 1.0 + SplineTolerance()) {
    throw std::invalid_argument("SplineInverse: coordinate is outside [0,1]");
  }

  const std::size_t bin = FindBin(y, heights);
  Scalar x_k = ScalarConstant(y, 0.0);
  Scalar y_k = ScalarConstant(y, 0.0);
  for (std::size_t i = 0; i < bin; ++i) {
    x_k += widths[i];
    y_k += heights[i];
  }

  const Scalar width = widths[bin];
  const Scalar height = heights[bin];
  const Scalar delta = height / width;
  const Scalar y_delta = y - y_k;
  const Scalar slope_term =
      derivatives[bin] + derivatives[bin + 1] - 2.0 * delta;
  const Scalar a = y_delta * slope_term + height * (delta - derivatives[bin]);
  const Scalar b = height * derivatives[bin] - y_delta * slope_term;
  const Scalar c = -delta * y_delta;

  // The quadratic expression also preserves parameter derivatives at a = 0
  const double coefficient_scale =
      SplineCoefficientScale(ScalarValue(a), ScalarValue(b), ScalarValue(c));
  const Scalar scale = ScalarConstant(a, coefficient_scale);
  const Scalar normalized_a = a / scale;
  const Scalar normalized_b = b / scale;
  const Scalar normalized_c = c / scale;
  Scalar discriminant =
      normalized_b * normalized_b - 4.0 * normalized_a * normalized_c;
  const double validated = ValidatedSplineDiscriminant(
      ScalarValue(normalized_a), ScalarValue(normalized_b),
      ScalarValue(normalized_c), "SplineInverse");
  if (std::fpclassify(validated) == FP_ZERO) {
    discriminant = ScalarConstant(discriminant, 0.0);
  }
  const Scalar root = ScalarSqrt(discriminant);
  // Choose the equivalent root expression without subtractive cancellation
  Scalar xi = ScalarValue(normalized_b) < 0.0
                  ? (-normalized_b + root) / (2.0 * normalized_a)
                  : 2.0 * normalized_c / (-normalized_b - root);

  const double xi_value = ScalarValue(xi);
  if (xi_value < -1e-8 || xi_value > 1.0 + 1e-8) {
    throw std::runtime_error("SplineInverse: analytic root is outside its bin");
  }
  if (xi_value < 0.0) {
    xi = ScalarConstant(xi, 0.0);
  }
  if (xi_value > 1.0) {
    xi = ScalarConstant(xi, 1.0);
  }

  const Scalar one_minus_xi = 1.0 - xi;
  const Scalar xi_product = xi * one_minus_xi;
  const Scalar denominator =
      delta +
      (derivatives[bin + 1] + derivatives[bin] - 2.0 * delta) * xi_product;
  const Scalar derivative_numerator =
      delta * delta *
      (derivatives[bin + 1] * xi * xi + 2.0 * delta * xi_product +
       derivatives[bin] * one_minus_xi * one_minus_xi);
  const Scalar forward_jacobian =
      derivative_numerator / (denominator * denominator);
  return {x_k + xi * width, -ScalarLog(forward_jacobian)};
}

// Compute the inverse softplus used for identity derivative initialization
double InverseSoftplus(double value) {
  if (!(value > 0.0)) {
    throw std::invalid_argument("InverseSoftplus: value must be positive");
  }
  return value > 20.0 ? value : std::log(std::expm1(value));
}

// Resample equally spaced VEGAS quantiles into spline output heights
std::vector<double>
ResampleVegasQuantileHeights(const std::vector<double> &quantile_edges,
                             std::size_t spline_bins) {
  if (quantile_edges.size() < 2 ||
      !std::is_eq(quantile_edges.front() <=> 0.0) ||
      !std::is_eq(quantile_edges.back() <=> 1.0)) {
    throw std::invalid_argument(
        "ResampleVegasQuantileHeights: grid must include zero and one");
  }
  for (std::size_t edge = 1; edge < quantile_edges.size(); ++edge) {
    if (!std::isfinite(quantile_edges[edge]) ||
        !(quantile_edges[edge] > quantile_edges[edge - 1])) {
      throw std::invalid_argument(
          "ResampleVegasQuantileHeights: grid must be finite and increasing");
    }
  }

  const std::size_t source_bins = quantile_edges.size() - 1;
  std::vector<double> resampled_edges(spline_bins + 1, 0.0);
  resampled_edges.back() = 1.0;
  for (std::size_t edge = 1; edge < spline_bins; ++edge) {
    const double source_position = static_cast<double>(edge) *
                                   static_cast<double>(source_bins) /
                                   static_cast<double>(spline_bins);
    const std::size_t lower = static_cast<std::size_t>(source_position);
    const double fraction = source_position - static_cast<double>(lower);
    resampled_edges[edge] =
        quantile_edges[lower] +
        fraction * (quantile_edges[lower + 1] - quantile_edges[lower]);
  }

  std::vector<double> heights(spline_bins, 0.0);
  for (std::size_t bin = 0; bin < spline_bins; ++bin) {
    heights[bin] = resampled_edges[bin + 1] - resampled_edges[bin];
    if (!(heights[bin] > 0.0) || !std::isfinite(heights[bin])) {
      throw std::invalid_argument(
          "ResampleVegasQuantileHeights: resampled grid is invalid");
    }
  }
  return heights;
}

// Convert one positive grid into spline heights above the configured floor
std::vector<double> RegularizeVegasWidths(const std::vector<double> &widths,
                                          double minimum) {
  const double available = 1.0 - minimum * static_cast<double>(widths.size());
  const double total = gra::Sum(widths);
  if (!(available > 0.0) || !(total > 0.0)) {
    throw std::invalid_argument(
        "RegularizeVegasWidths: invalid grid or spline floor");
  }
  std::vector<double> regularized(widths.size(), 0.0);
  bool above_floor = true;
  for (const auto &bin : indices(widths)) {
    regularized[bin] = widths[bin] / total;
    above_floor = above_floor && regularized[bin] > minimum;
  }
  if (above_floor) {
    return regularized;
  }
  for (const auto &bin : indices(widths)) {
    regularized[bin] = minimum + available * widths[bin] / total;
  }
  return regularized;
}

// Construct smooth monotonic knot slopes for one spline initialization grid
std::vector<double> VegasKnotDerivatives(const std::vector<double> &heights,
                                         double minimum_derivative) {
  const double bins = static_cast<double>(heights.size());
  std::vector<double> secants(heights.size(), 0.0);
  for (const auto &bin : indices(heights)) {
    secants[bin] = bins * heights[bin];
  }

  std::vector<double> derivatives(heights.size() + 1, 0.0);
  derivatives.front() = secants.front();
  derivatives.back() = secants.back();
  for (std::size_t knot = 1; knot < heights.size(); ++knot) {
    derivatives[knot] = 2.0 * secants[knot - 1] * secants[knot] /
                        (secants[knot - 1] + secants[knot]);
  }
  const double offset = std::max(1e-12, minimum_derivative * 1e-9);
  for (double &derivative : derivatives) {
    derivative = std::max(derivative, minimum_derivative + offset);
  }
  return derivatives;
}

// Construct the identity coordinate permutation
std::vector<std::size_t> IdentityPermutation(std::size_t dimension) {
  std::vector<std::size_t> permutation(dimension);
  std::iota(permutation.begin(), permutation.end(), 0);
  return permutation;
}

// Construct the inverse of one validated coordinate permutation
std::vector<std::size_t>
InversePermutation(const std::vector<std::size_t> &permutation) {
  std::vector<std::size_t> inverse(permutation.size(), permutation.size());
  for (const auto &output : indices(permutation)) {
    const std::size_t input = permutation[output];
    if (input >= permutation.size() || inverse[input] != permutation.size()) {
      throw std::invalid_argument(
          "InversePermutation: invalid coordinate permutation");
    }
    inverse[input] = output;
  }
  return inverse;
}

// Compose restoration from one coupling with the next input permutation
std::vector<std::size_t>
CouplingTransition(const std::vector<std::size_t> &inverse_current,
                   const std::vector<std::size_t> &next_permutation) {
  if (inverse_current.size() != next_permutation.size()) {
    throw std::invalid_argument(
        "CouplingTransition: permutation size mismatch");
  }
  std::vector<std::size_t> transition(next_permutation.size());
  for (const auto &coordinate : indices(transition)) {
    transition[coordinate] = inverse_current[next_permutation[coordinate]];
  }
  return transition;
}

// Construct a deterministic even-odd interleaving permutation
std::vector<std::size_t> InterleavePermutation(std::size_t dimension) {
  std::vector<std::size_t> permutation;
  permutation.reserve(dimension);
  for (std::size_t coordinate = 0; coordinate < dimension; coordinate += 2) {
    permutation.push_back(coordinate);
  }
  for (std::size_t coordinate = 1; coordinate < dimension; coordinate += 2) {
    permutation.push_back(coordinate);
  }
  return permutation;
}

// Construct one balanced per-coupling coordinate assignment
std::vector<std::size_t> BalancedPermutation(
    std::size_t dimension, const std::vector<std::size_t> &transform_positions,
    FlowPermutation mode, std::size_t layer,
    std::vector<std::size_t> &transform_counts, std::mt19937 &generator) {
  if (mode == FlowPermutation::None) {
    return IdentityPermutation(dimension);
  }

  std::vector<std::size_t> order = mode == FlowPermutation::Interleave
                                       ? InterleavePermutation(dimension)
                                       : IdentityPermutation(dimension);
  if (mode == FlowPermutation::Interleave && dimension > 1) {
    std::rotate(order.begin(),
                order.begin() + static_cast<std::ptrdiff_t>(layer % dimension),
                order.end());
  } else if (mode == FlowPermutation::Random) {
    std::shuffle(order.begin(), order.end(), generator);
  }
  std::stable_sort(order.begin(), order.end(),
                   [&transform_counts](std::size_t left, std::size_t right) {
                     return transform_counts[left] < transform_counts[right];
                   });

  std::vector<unsigned char> is_transform(dimension, 0);
  for (std::size_t position : transform_positions) {
    is_transform[position] = 1;
  }
  std::vector<std::size_t> permutation(dimension);
  std::size_t transformed = 0;
  std::size_t conditioned = transform_positions.size();
  for (std::size_t position = 0; position < dimension; ++position) {
    permutation[position] = is_transform[position] != 0 ? order[transformed++]
                                                        : order[conditioned++];
  }
  for (std::size_t position : transform_positions) {
    ++transform_counts[permutation[position]];
  }
  return permutation;
}

// Compute the scheduled learning rate for one zero-based optimizer update
// eta_cos(t) = eta_f + (eta_0-eta_f)[1+cos(pi t)]/2
double ScheduledLearningRate(const MNeuroJacConfig &config,
                             std::uint64_t step) {
  if (config.lr_schedule == FlowLRSchedule::None) {
    return config.learning_rate;
  }
  const std::size_t   batches_per_epoch = 1 + (config.buffer_size - 1) / config.batch_size;
  const std::uint64_t total_updates =
      static_cast<std::uint64_t>(config.rounds) * config.epochs *
      batches_per_epoch;
  if (total_updates <= 1) {
    return config.learning_rate;
  }
  const double progress = std::min(
      static_cast<double>(step) / static_cast<double>(total_updates - 1), 1.0);
  if (config.lr_schedule == FlowLRSchedule::Cosine) {
    const double blend = 0.5 * (1.0 + std::cos(math::PI * progress));
    return config.final_learning_rate +
           (config.learning_rate - config.final_learning_rate) * blend;
  }
  if (config.lr_schedule == FlowLRSchedule::Exponential) {
    return config.learning_rate *
           std::pow(config.final_learning_rate / config.learning_rate,
                    progress);
  }
  throw std::invalid_argument(
      "ScheduledLearningRate: unknown learning-rate schedule");
}

// Apply one coordinate permutation with output[i] equal to
// input[permutation[i]]
template <typename Scalar>
std::vector<Scalar>
ApplyPermutation(const std::vector<Scalar> &input,
                 const std::vector<std::size_t> &permutation) {
  if (input.size() != permutation.size()) {
    throw std::invalid_argument("ApplyPermutation: coordinate size mismatch");
  }
  std::vector<Scalar> output(input.size());
  for (const auto &coordinate : indices(input)) {
    output[coordinate] = input[permutation[coordinate]];
  }
  return output;
}

} // namespace

// AdamW update kernel with decoupled weight decay and global gradient clipping
class MAdamW {
public:
  // Reset optimizer moments and hyperparameters for a new flow
  static void Configure(std::size_t parameter_count,
                        std::vector<double> &first_moment,
                        std::vector<double> &second_moment,
                        std::uint64_t &step) {
    first_moment.assign(parameter_count, 0.0);
    second_moment.assign(parameter_count, 0.0);
    step = 0;
  }

  // Convert one gradient into bias-corrected AdamW parameter steps
  // step_i = eta m_i/(1-beta1^t)/[sqrt(v_i/(1-beta2^t))+epsilon]
  static void Steps(std::vector<double> &gradient, double learning_rate,
                    double gradient_clip, const MNeuroJacConfig &config,
                    std::vector<double> &first_moment,
                    std::vector<double> &second_moment, std::uint64_t &step) {
    if (gradient.size() != first_moment.size() ||
        gradient.size() != second_moment.size()) {
      throw std::invalid_argument(
          "MAdamW::Steps: gradient and optimizer sizes differ");
    }

    for (double value : gradient) {
      if (!std::isfinite(value)) {
        throw std::runtime_error("MAdamW::Steps: non-finite gradient");
      }
    }
    statistics::ClipEuclideanNorm(gradient, gradient_clip);

    ++step;
    const double first_correction =
        1.0 - std::pow(config.adamw_beta1, static_cast<double>(step));
    const double second_correction =
        1.0 - std::pow(config.adamw_beta2, static_cast<double>(step));
    for (const auto &i : indices(gradient)) {
      const double clipped = gradient[i];
      first_moment[i] = config.adamw_beta1 * first_moment[i] +
                        (1.0 - config.adamw_beta1) * clipped;
      second_moment[i] = config.adamw_beta2 * second_moment[i] +
                         (1.0 - config.adamw_beta2) * clipped * clipped;
      const double first_unbiased = first_moment[i] / first_correction;
      const double second_unbiased = second_moment[i] / second_correction;
      gradient[i] = learning_rate * first_unbiased /
                    (std::sqrt(second_unbiased) + config.adamw_epsilon);
    }
  }
};

// Internal evaluator shared by the double and reverse-mode flow paths
class MSplineEvaluator {
public:
  // Evaluate a coupling conditioner MLP from the unchanged mask coordinates
  template <typename Scalar, typename Parameters>
  static std::vector<Scalar> Conditioner(const MSplineCoupling &coupling,
                                         const std::vector<Scalar> &point,
                                         const Parameters &parameters) {
    std::vector<Scalar> values;
    values.reserve(coupling.condition_indices_.size());
    for (std::size_t index : coupling.condition_indices_) {
      values.push_back(point[index]);
    }

    for (const auto &layer : indices(coupling.dense_shapes_)) {
      const MSplineCoupling::DenseShape &shape = coupling.dense_shapes_[layer];
      std::vector<Scalar> output(shape.outputs);
      for (std::size_t row = 0; row < shape.outputs; ++row) {
        output[row] = parameters[shape.bias_offset + row];
        for (std::size_t column = 0; column < shape.inputs; ++column) {
          output[row] +=
              parameters[shape.weight_offset + row * shape.inputs + column] *
              values[column];
        }
        if (layer + 1 < coupling.dense_shapes_.size()) {
          output[row] = ScalarActivation(output[row], coupling.activation_);
        }
      }
      values = std::move(output);
    }
    return values;
  }

  // Apply one coupling in the generative direction and accumulate log dx over
  // dz
  template <typename Scalar, typename Parameters>
  static std::pair<std::vector<Scalar>, Scalar>
  CouplingForward(const MSplineCoupling &coupling,
                  const std::vector<Scalar> &input,
                  const Parameters &parameters) {
    std::vector<Scalar> output = input;
    Scalar log_jacobian = ScalarConstant(input.front(), 0.0);
    const std::vector<Scalar> conditioned =
        Conditioner(coupling, input, parameters);
    const std::size_t stride = 3 * coupling.bins_ + 1;

    for (const auto &coordinate : indices(coupling.transform_indices_)) {
      const std::size_t offset = coordinate * stride;
      std::vector<Scalar> raw_widths(coupling.bins_);
      std::vector<Scalar> raw_heights(coupling.bins_);
      std::vector<Scalar> raw_derivatives(coupling.bins_ + 1);
      for (std::size_t bin = 0; bin < coupling.bins_; ++bin) {
        raw_widths[bin] = conditioned[offset + bin];
        raw_heights[bin] = conditioned[offset + coupling.bins_ + bin];
      }
      for (std::size_t knot = 0; knot <= coupling.bins_; ++knot) {
        raw_derivatives[knot] = conditioned[offset + 2 * coupling.bins_ + knot];
      }

      const std::vector<Scalar> widths =
          NormalizeBins(raw_widths, coupling.min_bin_);
      const std::vector<Scalar> heights =
          NormalizeBins(raw_heights, coupling.min_bin_);
      const std::vector<Scalar> derivatives =
          NormalizeDerivatives(raw_derivatives, coupling.min_derivative_);
      const auto transformed =
          SplineForward(input[coupling.transform_indices_[coordinate]], widths,
                        heights, derivatives);
      output[coupling.transform_indices_[coordinate]] = transformed.first;
      log_jacobian += transformed.second;
    }
    return {output, log_jacobian};
  }

  // Apply one coupling in the density direction and accumulate log dz over dx
  template <typename Scalar, typename Parameters>
  static std::pair<std::vector<Scalar>, Scalar>
  CouplingInverse(const MSplineCoupling &coupling,
                  const std::vector<Scalar> &input,
                  const Parameters &parameters) {
    std::vector<Scalar> output = input;
    Scalar log_jacobian = ScalarConstant(input.front(), 0.0);
    const std::vector<Scalar> conditioned =
        Conditioner(coupling, input, parameters);
    const std::size_t stride = 3 * coupling.bins_ + 1;

    for (const auto &coordinate : indices(coupling.transform_indices_)) {
      const std::size_t offset = coordinate * stride;
      std::vector<Scalar> raw_widths(coupling.bins_);
      std::vector<Scalar> raw_heights(coupling.bins_);
      std::vector<Scalar> raw_derivatives(coupling.bins_ + 1);
      for (std::size_t bin = 0; bin < coupling.bins_; ++bin) {
        raw_widths[bin] = conditioned[offset + bin];
        raw_heights[bin] = conditioned[offset + coupling.bins_ + bin];
      }
      for (std::size_t knot = 0; knot <= coupling.bins_; ++knot) {
        raw_derivatives[knot] = conditioned[offset + 2 * coupling.bins_ + knot];
      }

      const std::vector<Scalar> widths =
          NormalizeBins(raw_widths, coupling.min_bin_);
      const std::vector<Scalar> heights =
          NormalizeBins(raw_heights, coupling.min_bin_);
      const std::vector<Scalar> derivatives =
          NormalizeDerivatives(raw_derivatives, coupling.min_derivative_);
      const auto transformed =
          SplineInverse(input[coupling.transform_indices_[coordinate]], widths,
                        heights, derivatives);
      output[coupling.transform_indices_[coordinate]] = transformed.first;
      log_jacobian += transformed.second;
    }
    return {output, log_jacobian};
  }

  // Apply every coupling in the generative direction
  template <typename Scalar, typename Parameters>
  static std::pair<std::vector<Scalar>, Scalar>
  FlowForward(const MSplineFlow &flow, const std::vector<Scalar> &base,
              const Parameters &parameters) {
    std::vector<Scalar> point = base;
    Scalar log_jacobian = ScalarConstant(base.front(), 0.0);
    for (const MSplineCoupling &coupling : flow.Couplings()) {
      const std::vector<Scalar> permuted =
          ApplyPermutation(point, coupling.permutation_);
      auto transformed = CouplingForward(coupling, permuted, parameters);
      point =
          ApplyPermutation(transformed.first, coupling.inverse_permutation_);
      log_jacobian += transformed.second;
    }
    return {point, log_jacobian};
  }

  // Apply every coupling in reverse to evaluate the exact flow log density
  template <typename Scalar, typename Parameters>
  static Scalar FlowLogDensity(const MSplineFlow &flow,
                               const std::vector<Scalar> &point,
                               const Parameters &parameters) {
    std::vector<Scalar> base = point;
    Scalar log_density = ScalarConstant(point.front(), 0.0);
    for (auto coupling = flow.Couplings().rbegin();
         coupling != flow.Couplings().rend(); ++coupling) {
      const std::vector<Scalar> permuted =
          ApplyPermutation(base, coupling->permutation_);
      auto transformed = CouplingInverse(*coupling, permuted, parameters);
      base =
          ApplyPermutation(transformed.first, coupling->inverse_permutation_);
      log_density += transformed.second;
    }
    return log_density;
  }
};

#ifndef GRANIITTI_USE_LIBTORCH

using MNeuroJacBatch =
    Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>;
using MNeuroJacBatchVector = Eigen::VectorXd;

// Immutable standalone dense parameter views for one coupling
struct MStandaloneDenseInference {
  // Bind matrix and bias views to one owned flat parameter block
  MStandaloneDenseInference(const double *weight_data, const double *bias_data,
                            std::size_t inputs, std::size_t outputs)
      : weights(weight_data, static_cast<Eigen::Index>(outputs),
                static_cast<Eigen::Index>(inputs)),
        biases(bias_data, static_cast<Eigen::Index>(outputs)) {}

  Eigen::Map<const MNeuroJacBatch> weights;
  Eigen::Map<const MNeuroJacBatchVector> biases;
};

// Immutable standalone evaluation metadata for one coupling
struct MStandaloneCouplingInference {
  std::vector<MStandaloneDenseInference> dense_layers;
  std::vector<std::size_t> to_next;
  std::vector<std::size_t> to_previous;
};

// Immutable standalone evaluation state for one complete flow
struct MFlowInferenceState {
  std::vector<double> parameters;
  std::vector<MStandaloneCouplingInference> couplings;
};

// Build owned standalone parameters and persistent views for evaluation
std::shared_ptr<const void>
BuildInferenceState(const MSplineFlow &flow,
                    const std::vector<double> &parameters) {
  auto state = std::make_shared<MFlowInferenceState>();
  state->parameters = parameters;
  state->couplings.reserve(flow.Couplings().size());
  for (const auto &index : indices(flow.Couplings())) {
    const MSplineCoupling &coupling = flow.Couplings()[index];
    MStandaloneCouplingInference coupling_state;
    coupling_state.dense_layers.reserve(coupling.dense_shapes_.size());
    for (const MSplineCoupling::DenseShape &shape : coupling.dense_shapes_) {
      coupling_state.dense_layers.emplace_back(
          state->parameters.data() + shape.weight_offset,
          state->parameters.data() + shape.bias_offset, shape.inputs,
          shape.outputs);
    }
    coupling_state.to_next =
        index + 1 < flow.Couplings().size()
            ? CouplingTransition(coupling.inverse_permutation_,
                                 flow.Couplings()[index + 1].permutation_)
            : coupling.inverse_permutation_;
    coupling_state.to_previous =
        index > 0 ? CouplingTransition(coupling.inverse_permutation_,
                                       flow.Couplings()[index - 1].permutation_)
                  : coupling.inverse_permutation_;
    state->couplings.push_back(std::move(coupling_state));
  }
  return state;
}

// Cached inputs and activations for one standalone coupling reverse sweep
struct MStandaloneCouplingCache {
  std::size_t coupling_index = 0;
  MNeuroJacBatch input;
  std::vector<MNeuroJacBatch> activations;
  std::vector<MNeuroJacBatch> preactivations;
};

// Reusable scalar buffers for standalone spline batch evaluation
struct MStandaloneSplineWorkspace {
  // Allocate all scalar spline buffers for one configured bin count
  explicit MStandaloneSplineWorkspace(std::size_t bins)
      : widths(bins), heights(bins), derivatives(bins + 1) {}

  std::vector<double> widths;
  std::vector<double> heights;
  std::vector<double> derivatives;
};

// Reusable local reverse tape for one inverse spline and its log Jacobian
class MStandaloneSplineAdjoint {
public:
  // Allocate the differentiable spline buffers once per coupling sweep
  explicit MStandaloneSplineAdjoint(std::size_t bins)
      : raw_widths_(bins), raw_heights_(bins), raw_derivatives_(bins + 1),
        widths_(bins), heights_(bins), derivatives_(bins + 1) {
    tape_.Reserve(96 * bins + 128);
  }

  // Backpropagate weighted output and log-Jacobian seeds through one inverse
  // spline
  double Inverse(double coordinate, const double *raw, std::size_t bins,
                 double minimum_bin, double minimum_derivative,
                 double output_seed, double log_jacobian_seed,
                 double *raw_adjoint) {
    tape_.Reset();
    const MReverseScalar input = MReverseScalar::Parameter(tape_, coordinate);
    for (std::size_t bin = 0; bin < bins; ++bin) {
      raw_widths_[bin] = MReverseScalar::Parameter(tape_, raw[bin]);
      raw_heights_[bin] = MReverseScalar::Parameter(tape_, raw[bins + bin]);
    }
    for (std::size_t knot = 0; knot <= bins; ++knot) {
      raw_derivatives_[knot] =
          MReverseScalar::Parameter(tape_, raw[2 * bins + knot]);
    }
    NormalizeBins(raw_widths_, minimum_bin, widths_);
    NormalizeBins(raw_heights_, minimum_bin, heights_);
    const std::size_t bin = FindBin(input, heights_);
    for (std::size_t knot = bin; knot <= bin + 1; ++knot) {
      derivatives_[knot] = minimum_derivative + ScalarSoftplus(raw_derivatives_[knot]);
    }
    const auto transformed =
        SplineInverse(input, widths_, heights_, derivatives_);
    const MReverseScalar weighted = output_seed * transformed.first +
                                    log_jacobian_seed * transformed.second;
    tape_.Backward(weighted.NodeIndex());

    for (std::size_t bin = 0; bin < bins; ++bin) {
      raw_adjoint[bin] = tape_.Adjoint(raw_widths_[bin].NodeIndex());
      raw_adjoint[bins + bin] = tape_.Adjoint(raw_heights_[bin].NodeIndex());
    }
    for (std::size_t knot = 0; knot <= bins; ++knot) {
      raw_adjoint[2 * bins + knot] =
          tape_.Adjoint(raw_derivatives_[knot].NodeIndex());
    }
    return tape_.Adjoint(input.NodeIndex());
  }

private:
  MReverseTape tape_;
  std::vector<MReverseScalar> raw_widths_;
  std::vector<MReverseScalar> raw_heights_;
  std::vector<MReverseScalar> raw_derivatives_;
  std::vector<MReverseScalar> widths_;
  std::vector<MReverseScalar> heights_;
  std::vector<MReverseScalar> derivatives_;
};

// Convert unconstrained logits into normalized double-precision spline bins
void StandaloneNormalizeBins(const double *raw, std::size_t bins,
                             double minimum, std::vector<double> &normalized) {
  double maximum = raw[0];
  for (std::size_t bin = 1; bin < bins; ++bin) {
    maximum = std::max(maximum, raw[bin]);
  }
  double denominator = 0.0;
  for (std::size_t bin = 0; bin < bins; ++bin) {
    normalized[bin] = std::exp(raw[bin] - maximum);
    denominator += normalized[bin];
  }
  const double available = 1.0 - minimum * static_cast<double>(bins);
  const double scale = available / denominator;
  for (double &entry : normalized) {
    entry = minimum + scale * entry;
  }
}

// Evaluate one standalone spline directly from one conditioner output row
MSplineResult StandaloneSpline(double coordinate, const double *raw,
                               const MSplineCoupling &coupling, bool inverse,
                               MStandaloneSplineWorkspace &workspace) {
  StandaloneNormalizeBins(raw, coupling.bins_, coupling.min_bin_,
                          workspace.widths);
  StandaloneNormalizeBins(raw + coupling.bins_, coupling.bins_,
                          coupling.min_bin_, workspace.heights);
  const std::size_t bin = FindBin(coordinate, inverse ? workspace.heights : workspace.widths);
  for (std::size_t knot = bin; knot <= bin + 1; ++knot) {
    workspace.derivatives[knot] = coupling.min_derivative_ + ScalarSoftplus(raw[2 * coupling.bins_ + knot]);
  }
  const auto transformed =
      inverse ? SplineInverse(coordinate, workspace.widths, workspace.heights,
                              workspace.derivatives)
              : SplineForward(coordinate, workspace.widths, workspace.heights,
                              workspace.derivatives);
  return {transformed.first, transformed.second};
}

// Apply one coordinate permutation to every row of a standalone batch
MNeuroJacBatch
StandalonePermutation(const Eigen::Ref<const MNeuroJacBatch> &input,
                      const std::vector<std::size_t> &permutation) {
  MNeuroJacBatch output(input.rows(), input.cols());
  for (const auto &coordinate : indices(permutation)) {
    output.col(static_cast<Eigen::Index>(coordinate)) =
        input.col(static_cast<Eigen::Index>(permutation[coordinate]));
  }
  return output;
}

// Apply one hidden-layer activation to a standalone batch
MNeuroJacBatch StandaloneActivation(const MNeuroJacBatch &input,
                                    FlowActivation activation) {
  return input.unaryExpr([activation](double value) {
    return ScalarActivation(value, activation);
  });
}

// Compute one hidden-layer activation derivative for a standalone batch
MNeuroJacBatch StandaloneActivationDerivative(const MNeuroJacBatch &input,
                                              FlowActivation activation) {
  return input.unaryExpr([activation](double value) {
    if (activation == FlowActivation::Relu) {
      return value > 0.0 ? 1.0 : 0.0;
    }
    if (activation == FlowActivation::Silu) {
      const double sigmoid = ScalarSigmoid(value);
      return sigmoid * (1.0 + value * (1.0 - sigmoid));
    }
    if (activation == FlowActivation::Tanh) {
      const double result = std::tanh(value);
      return 1.0 - result * result;
    }
    throw std::invalid_argument(
        "StandaloneActivationDerivative: unknown activation function");
  });
}

// Evaluate one standalone conditioner with vectorized dense matrix kernels
MNeuroJacBatch
StandaloneConditioner(const MSplineCoupling &coupling,
                      const MNeuroJacBatch &points,
                      const std::vector<double> &parameters,
                      const MStandaloneCouplingInference *inference,
                      std::vector<MNeuroJacBatch> *activations,
                      std::vector<MNeuroJacBatch> *preactivations) {
  MNeuroJacBatch values(points.rows(), static_cast<Eigen::Index>(
                                           coupling.condition_indices_.size()));
  for (const auto &column : indices(coupling.condition_indices_)) {
    values.col(static_cast<Eigen::Index>(column)) = points.col(
        static_cast<Eigen::Index>(coupling.condition_indices_[column]));
  }
  if (activations != nullptr) {
    activations->clear();
    activations->reserve(coupling.dense_shapes_.size() + 1);
    activations->push_back(values);
    preactivations->clear();
    preactivations->reserve(coupling.dense_shapes_.size());
  }

  for (const auto &layer : indices(coupling.dense_shapes_)) {
    const MSplineCoupling::DenseShape &shape = coupling.dense_shapes_[layer];
    MNeuroJacBatch output;
    if (inference != nullptr) {
      const MStandaloneDenseInference &dense = inference->dense_layers[layer];
      output = values * dense.weights.transpose();
      output.rowwise() += dense.biases.transpose();
    } else {
      const Eigen::Map<const MNeuroJacBatch> weights(
          parameters.data() + shape.weight_offset,
          static_cast<Eigen::Index>(shape.outputs),
          static_cast<Eigen::Index>(shape.inputs));
      const Eigen::Map<const MNeuroJacBatchVector> biases(
          parameters.data() + shape.bias_offset,
          static_cast<Eigen::Index>(shape.outputs));
      output = values * weights.transpose();
      output.rowwise() += biases.transpose();
    }
    if (preactivations != nullptr) {
      preactivations->push_back(output);
    }
    if (layer + 1 < coupling.dense_shapes_.size()) {
      output = StandaloneActivation(output, coupling.activation_);
    }
    values.swap(output);
    if (activations != nullptr) {
      activations->push_back(values);
    }
  }
  return values;
}

// Apply one standalone spline coupling to a complete point batch
std::pair<MNeuroJacBatch, MNeuroJacBatchVector>
StandaloneCoupling(const MSplineCoupling &coupling, const MNeuroJacBatch &input,
                   const std::vector<double> &parameters, bool inverse,
                   std::vector<MNeuroJacBatch> *activations,
                   std::vector<MNeuroJacBatch> *preactivations) {
  MNeuroJacBatch output = input;
  MNeuroJacBatchVector log_jacobian = MNeuroJacBatchVector::Zero(input.rows());
  const MNeuroJacBatch conditioned = StandaloneConditioner(
      coupling, input, parameters, nullptr, activations, preactivations);
  const std::size_t stride = 3 * coupling.bins_ + 1;
  MStandaloneSplineWorkspace workspace(coupling.bins_);

  for (const auto &coordinate : indices(coupling.transform_indices_)) {
    const std::size_t input_coordinate =
        coupling.transform_indices_[coordinate];
    const std::size_t offset = coordinate * stride;
    for (Eigen::Index event = 0; event < input.rows(); ++event) {
      const double *raw = conditioned.data() + event * conditioned.cols() +
                          static_cast<Eigen::Index>(offset);
      const MSplineResult transformed = StandaloneSpline(
          input(event, static_cast<Eigen::Index>(input_coordinate)), raw,
          coupling, inverse, workspace);
      output(event, static_cast<Eigen::Index>(input_coordinate)) =
          transformed.value;
      log_jacobian[event] += transformed.log_jacobian;
    }
  }
  return {std::move(output), std::move(log_jacobian)};
}

// Normalize every row of one spline-logit block into positive bins
MNeuroJacBatch StandaloneNormalizedBinRows(const Eigen::Ref<const MNeuroJacBatch> &raw,
                                           Eigen::Index offset,
                                           Eigen::Index bins, double minimum) {
  MNeuroJacBatch normalized = raw.middleCols(offset, bins);
  const MNeuroJacBatchVector maxima = normalized.rowwise().maxCoeff();
  normalized.colwise() -= maxima;
  normalized = normalized.array().exp();
  const MNeuroJacBatchVector sums = normalized.rowwise().sum();
  normalized.array().colwise() /= sums.array();
  const double available = 1.0 - minimum * static_cast<double>(bins);
  return (minimum + available * normalized.array()).matrix();
}

// Selected spline-knot data for one flattened coordinate batch
struct MStandaloneSplineSelection {
  MNeuroJacBatchVector x_k;
  MNeuroJacBatchVector y_k;
  MNeuroJacBatchVector width;
  MNeuroJacBatchVector height;
  MNeuroJacBatchVector d0;
  MNeuroJacBatchVector d1;
};

// Select spline bins for all flattened transformed coordinates
MStandaloneSplineSelection
StandaloneSelectSplineBins(const Eigen::Ref<const MNeuroJacBatchVector> &coordinates,
                           const MNeuroJacBatch &widths,
                           const MNeuroJacBatch &heights,
                           const Eigen::Ref<const MNeuroJacBatch> &raw, double minimum_derivative, bool inverse) {
  const Eigen::Index rows = coordinates.size();
  MStandaloneSplineSelection selected{
      MNeuroJacBatchVector(rows), MNeuroJacBatchVector(rows),
      MNeuroJacBatchVector(rows), MNeuroJacBatchVector(rows),
      MNeuroJacBatchVector(rows), MNeuroJacBatchVector(rows)};
  for (Eigen::Index row = 0; row < rows; ++row) {
    double cumulative_width = 0.0;
    double cumulative_height = 0.0;
    Eigen::Index selected_bin = widths.cols() - 1;
    for (Eigen::Index bin = 0; bin + 1 < widths.cols(); ++bin) {
      const double boundary = inverse ? cumulative_height + heights(row, bin)
                                      : cumulative_width + widths(row, bin);
      if (coordinates[row] < boundary) {
        selected_bin = bin;
        break;
      }
      cumulative_width += widths(row, bin);
      cumulative_height += heights(row, bin);
    }
    selected.x_k[row] = cumulative_width;
    selected.y_k[row] = cumulative_height;
    selected.width[row] = widths(row, selected_bin);
    selected.height[row] = heights(row, selected_bin);
    selected.d0[row] = minimum_derivative + ScalarSoftplus(raw(row, 2 * widths.cols() + selected_bin));
    selected.d1[row] = minimum_derivative + ScalarSoftplus(raw(row, 2 * widths.cols() + selected_bin + 1));
  }
  return selected;
}

// Solve the analytic inverse coordinate for every selected spline bin
MNeuroJacBatchVector
StandaloneInverseSplineCoordinates(const Eigen::Ref<const MNeuroJacBatchVector> &coordinates,
                                   const MStandaloneSplineSelection &selected) {
  const MNeuroJacBatchVector delta =
      (selected.height.array() / selected.width.array()).matrix();
  const MNeuroJacBatchVector y_delta =
      (coordinates.array() - selected.y_k.array()).matrix();
  const MNeuroJacBatchVector slope_term =
      (selected.d0.array() + selected.d1.array() - 2.0 * delta.array())
          .matrix();
  const MNeuroJacBatchVector a =
      (y_delta.array() * slope_term.array() +
       selected.height.array() * (delta.array() - selected.d0.array()))
          .matrix();
  const MNeuroJacBatchVector b =
      (selected.height.array() * selected.d0.array() -
       y_delta.array() * slope_term.array())
          .matrix();
  const MNeuroJacBatchVector c = (-delta.array() * y_delta.array()).matrix();
  MNeuroJacBatchVector xi(coordinates.size());
  for (const auto &row : indices(coordinates)) {
    if (SplineUsesLinearRoot(a[row], b[row], c[row])) {
      xi[row] = -c[row] / b[row];
    } else {
      const double coefficient_scale =
          SplineCoefficientScale(a[row], b[row], c[row]);
      const double normalized_a = a[row] / coefficient_scale;
      const double normalized_b = b[row] / coefficient_scale;
      const double normalized_c = c[row] / coefficient_scale;
      const double discriminant =
          ValidatedSplineDiscriminant(normalized_a, normalized_b, normalized_c,
                                      "StandaloneInverseSplineCoordinates");
      const double root = std::sqrt(discriminant);
      xi[row] = normalized_b < 0.0
                    ? (-normalized_b + root) / (2.0 * normalized_a)
                    : 2.0 * normalized_c / (-normalized_b - root);
    }
    if (xi[row] < -1e-8 || xi[row] > 1.0 + 1e-8) {
      throw std::runtime_error(
          "StandaloneInverseSplineCoordinates: root outside its bin");
    }
    xi[row] = std::clamp(xi[row], 0.0, 1.0);
  }
  return xi;
}

// Evaluate all transformed coordinates of one spline coupling together
std::pair<MNeuroJacBatchVector, MNeuroJacBatchVector>
StandaloneSplineCoordinates(const Eigen::Ref<const MNeuroJacBatchVector> &coordinates,
                            const Eigen::Ref<const MNeuroJacBatch> &raw,
                            const MSplineCoupling &coupling, bool inverse) {
  const Eigen::Index bins = static_cast<Eigen::Index>(coupling.bins_);
  const MNeuroJacBatch widths =
      StandaloneNormalizedBinRows(raw, 0, bins, coupling.min_bin_);
  const MNeuroJacBatch heights =
      StandaloneNormalizedBinRows(raw, bins, bins, coupling.min_bin_);
  const MStandaloneSplineSelection selected = StandaloneSelectSplineBins(
      coordinates, widths, heights, raw, coupling.min_derivative_, inverse);
  const MNeuroJacBatchVector delta =
      (selected.height.array() / selected.width.array()).matrix();
  const MNeuroJacBatchVector xi =
      inverse ? StandaloneInverseSplineCoordinates(coordinates, selected)
              : ((coordinates.array() - selected.x_k.array()) /
                 selected.width.array())
                    .matrix();
  const MNeuroJacBatchVector one_minus_xi = (1.0 - xi.array()).matrix();
  const MNeuroJacBatchVector xi_product =
      (xi.array() * one_minus_xi.array()).matrix();
  const MNeuroJacBatchVector denominator =
      (delta.array() +
       (selected.d1.array() + selected.d0.array() - 2.0 * delta.array()) *
           xi_product.array())
          .matrix();
  const MNeuroJacBatchVector derivative_numerator =
      (delta.array().square() *
       (selected.d1.array() * xi.array().square() +
        2.0 * delta.array() * xi_product.array() +
        selected.d0.array() * one_minus_xi.array().square()))
          .matrix();
  const MNeuroJacBatchVector jacobian =
      (derivative_numerator.array() / denominator.array().square()).matrix();
  MNeuroJacBatchVector values(coordinates.size());
  MNeuroJacBatchVector log_jacobian(coordinates.size());
  if (inverse) {
    values =
        (selected.x_k.array() + xi.array() * selected.width.array()).matrix();
    log_jacobian = (-jacobian.array().log()).matrix();
  } else {
    values =
        (selected.y_k.array() + selected.height.array() *
                                    (delta.array() * xi.array().square() +
                                     selected.d0.array() * xi_product.array()) /
                                    denominator.array())
            .matrix();
    log_jacobian = jacobian.array().log().matrix();
  }
  return {std::move(values), std::move(log_jacobian)};
}

// Apply one standalone coupling with flattened coordinate vectorization
std::pair<MNeuroJacBatch, MNeuroJacBatchVector> StandaloneCouplingEvaluation(
    const MSplineCoupling &coupling, const MNeuroJacBatch &input,
    const std::vector<double> &parameters,
    const MStandaloneCouplingInference &inference, bool inverse) {
  const MNeuroJacBatch conditioned = StandaloneConditioner(
      coupling, input, parameters, &inference, nullptr, nullptr);
  const Eigen::Index transforms =
      static_cast<Eigen::Index>(coupling.transform_indices_.size());
  const Eigen::Index stride = static_cast<Eigen::Index>(3 * coupling.bins_ + 1);
  const Eigen::Map<const MNeuroJacBatch> raw(conditioned.data(),
                                             input.rows() * transforms, stride);
  MNeuroJacBatch coordinates(input.rows(), transforms);
  for (Eigen::Index transform = 0; transform < transforms; ++transform) {
    coordinates.col(transform) = input.col(static_cast<Eigen::Index>(
        coupling.transform_indices_[static_cast<std::size_t>(transform)]));
  }
  const Eigen::Map<const MNeuroJacBatchVector> flat_coordinates(
      coordinates.data(), coordinates.size());
  auto transformed =
      StandaloneSplineCoordinates(flat_coordinates, raw, coupling, inverse);
  const Eigen::Map<const MNeuroJacBatch> transformed_coordinates(
      transformed.first.data(), input.rows(), transforms);
  const Eigen::Map<const MNeuroJacBatch> coordinate_log_jacobians(
      transformed.second.data(), input.rows(), transforms);
  MNeuroJacBatch output = input;
  for (Eigen::Index transform = 0; transform < transforms; ++transform) {
    output.col(static_cast<Eigen::Index>(
        coupling.transform_indices_[static_cast<std::size_t>(transform)])) =
        transformed_coordinates.col(transform);
  }
  MNeuroJacBatchVector log_jacobian = coordinate_log_jacobians.rowwise().sum();
  return {std::move(output), std::move(log_jacobian)};
}

// Apply every standalone coupling in the generative direction
std::pair<MNeuroJacBatch, MNeuroJacBatchVector>
StandaloneFlowForward(const MSplineFlow &flow,
                      const Eigen::Ref<const MNeuroJacBatch> &base,
                      const std::vector<double> &parameters,
                      const MFlowInferenceState *inference) {
  MNeuroJacBatch points =
      StandalonePermutation(base, flow.Couplings().front().permutation_);
  MNeuroJacBatchVector log_jacobian = MNeuroJacBatchVector::Zero(base.rows());
  for (const auto &index : indices(flow.Couplings())) {
    const MSplineCoupling &coupling = flow.Couplings()[index];
    auto transformed = StandaloneCouplingEvaluation(
        coupling, points, parameters, inference->couplings[index], false);
    points = StandalonePermutation(transformed.first,
                                   inference->couplings[index].to_next);
    log_jacobian += transformed.second;
  }
  return {std::move(points), std::move(log_jacobian)};
}

// Apply every standalone coupling in the density direction
MNeuroJacBatchVector
StandaloneFlowLogDensity(const MSplineFlow &flow,
                         const Eigen::Ref<const MNeuroJacBatch> &points,
                         const std::vector<double> &parameters,
                         std::vector<MStandaloneCouplingCache> *caches,
                         const MFlowInferenceState *inference) {
  MNeuroJacBatchVector log_density = MNeuroJacBatchVector::Zero(points.rows());
  if (caches == nullptr) {
    MNeuroJacBatch base =
        StandalonePermutation(points, flow.Couplings().back().permutation_);
    for (std::size_t reverse = flow.Couplings().size(); reverse-- > 0;) {
      const MSplineCoupling &coupling = flow.Couplings()[reverse];
      auto transformed = StandaloneCouplingEvaluation(
          coupling, base, parameters, inference->couplings[reverse], true);
      base = StandalonePermutation(transformed.first,
                                   inference->couplings[reverse].to_previous);
      log_density += transformed.second;
    }
    return log_density;
  }

  MNeuroJacBatch base = points;
  caches->clear();
  caches->reserve(flow.Couplings().size());
  for (std::size_t reverse = flow.Couplings().size(); reverse-- > 0;) {
    const MSplineCoupling &coupling = flow.Couplings()[reverse];
    MNeuroJacBatch permuted =
        StandalonePermutation(base, coupling.permutation_);
    MStandaloneCouplingCache cache;
    cache.coupling_index = reverse;
    cache.input = permuted;
    auto transformed =
        StandaloneCoupling(coupling, permuted, parameters, true,
                           &cache.activations, &cache.preactivations);
    base =
        StandalonePermutation(transformed.first, coupling.inverse_permutation_);
    log_density += transformed.second;
    caches->push_back(std::move(cache));
  }
  return log_density;
}

// Backpropagate one conditioner batch and accumulate flat parameter gradients
MNeuroJacBatch StandaloneConditionerBackward(
    const MSplineCoupling &coupling,
    const std::vector<MNeuroJacBatch> &activations,
    const std::vector<MNeuroJacBatch> &preactivations, MNeuroJacBatch adjoint,
    const std::vector<double> &parameters, std::vector<double> &gradient) {
  for (std::size_t layer = coupling.dense_shapes_.size(); layer-- > 0;) {
    const MSplineCoupling::DenseShape &shape = coupling.dense_shapes_[layer];
    if (layer + 1 < coupling.dense_shapes_.size()) {
      adjoint.array() *= StandaloneActivationDerivative(preactivations[layer],
                                                        coupling.activation_)
                             .array();
    }
    Eigen::Map<MNeuroJacBatch> weight_gradient(
        gradient.data() + shape.weight_offset,
        static_cast<Eigen::Index>(shape.outputs),
        static_cast<Eigen::Index>(shape.inputs));
    Eigen::Map<MNeuroJacBatchVector> bias_gradient(
        gradient.data() + shape.bias_offset,
        static_cast<Eigen::Index>(shape.outputs));
    weight_gradient.noalias() += adjoint.transpose() * activations[layer];
    bias_gradient += adjoint.colwise().sum().transpose();

    const Eigen::Map<const MNeuroJacBatch> weights(
        parameters.data() + shape.weight_offset,
        static_cast<Eigen::Index>(shape.outputs),
        static_cast<Eigen::Index>(shape.inputs));
    MNeuroJacBatch previous_adjoint = adjoint * weights;
    adjoint.swap(previous_adjoint);
  }
  return adjoint;
}

// Backpropagate through one cached inverse coupling and its permutation
MNeuroJacBatch StandaloneCouplingBackward(
    const MSplineCoupling &coupling, const MStandaloneCouplingCache &cache,
    const MNeuroJacBatch &output_adjoint,
    const MNeuroJacBatchVector &log_jacobian_adjoint,
    const std::vector<double> &parameters, std::vector<double> &gradient) {
  MNeuroJacBatch transformed_adjoint =
      StandalonePermutation(output_adjoint, coupling.permutation_);
  MNeuroJacBatch input_adjoint = transformed_adjoint;
  MNeuroJacBatch conditioned_adjoint =
      MNeuroJacBatch::Zero(cache.input.rows(), cache.activations.back().cols());
  const MNeuroJacBatch &conditioned = cache.activations.back();
  const std::size_t stride = 3 * coupling.bins_ + 1;
  MStandaloneSplineAdjoint spline_adjoint(coupling.bins_);

  for (const auto &coordinate : indices(coupling.transform_indices_)) {
    const std::size_t input_coordinate =
        coupling.transform_indices_[coordinate];
    const std::size_t offset = coordinate * stride;
    for (Eigen::Index event = 0; event < cache.input.rows(); ++event) {
      const double *raw = conditioned.data() + event * conditioned.cols() +
                          static_cast<Eigen::Index>(offset);
      double *raw_adjoint = conditioned_adjoint.data() +
                            event * conditioned_adjoint.cols() +
                            static_cast<Eigen::Index>(offset);
      input_adjoint(event, static_cast<Eigen::Index>(input_coordinate)) =
          spline_adjoint.Inverse(
              cache.input(event, static_cast<Eigen::Index>(input_coordinate)),
              raw, coupling.bins_, coupling.min_bin_, coupling.min_derivative_,
              transformed_adjoint(event,
                                  static_cast<Eigen::Index>(input_coordinate)),
              log_jacobian_adjoint[event], raw_adjoint);
    }
  }

  const MNeuroJacBatch condition_adjoint = StandaloneConditionerBackward(
      coupling, cache.activations, cache.preactivations,
      std::move(conditioned_adjoint), parameters, gradient);
  for (const auto &column : indices(coupling.condition_indices_)) {
    input_adjoint.col(
        static_cast<Eigen::Index>(coupling.condition_indices_[column])) +=
        condition_adjoint.col(static_cast<Eigen::Index>(column));
  }
  return StandalonePermutation(input_adjoint, coupling.inverse_permutation_);
}

// Backpropagate a standalone flow-density batch into the flat parameter vector
void StandaloneFlowBackward(const MSplineFlow &flow,
                            const std::vector<MStandaloneCouplingCache> &caches,
                            const MNeuroJacBatchVector &log_density_adjoint,
                            const std::vector<double> &parameters,
                            std::vector<double> &gradient) {
  MNeuroJacBatch point_adjoint = MNeuroJacBatch::Zero(
      log_density_adjoint.size(), static_cast<Eigen::Index>(flow.Dimension()));
  for (std::size_t application = caches.size(); application-- > 0;) {
    const MStandaloneCouplingCache &cache = caches[application];
    point_adjoint = StandaloneCouplingBackward(
        flow.Couplings()[cache.coupling_index], cache, point_adjoint,
        log_density_adjoint, parameters, gradient);
  }
}

#endif

#ifdef GRANIITTI_USE_LIBTORCH

// Compute the common CPU double-precision tensor options for NEUROJAC
torch::TensorOptions NeuroJacTensorOptions() {
  return torch::TensorOptions().dtype(torch::kFloat64).device(torch::kCPU);
}

// Construct one owned tensor from the flat flow parameter vector
torch::Tensor TensorParameters(const std::vector<double> &parameters,
                               bool requires_gradient) {
  torch::Tensor tensor = torch::empty(
      {static_cast<std::int64_t>(parameters.size())}, NeuroJacTensorOptions());
  std::copy(parameters.begin(), parameters.end(), tensor.data_ptr<double>());
  tensor.set_requires_grad(requires_gradient);
  return tensor;
}

// Construct one immutable tensor view over flat row-major unit points
torch::Tensor TensorFlatPoints(const std::vector<double> &points,
                               std::size_t events, std::size_t dimension) {
  return torch::from_blob(
      const_cast<double *>(points.data()),
      {static_cast<std::int64_t>(events), static_cast<std::int64_t>(dimension)},
      NeuroJacTensorOptions());
}

// Construct one CPU index tensor from coordinate indices
torch::Tensor TensorIndices(const std::vector<std::size_t> &indices) {
  std::vector<std::int64_t> converted(indices.begin(), indices.end());
  return torch::tensor(converted, torch::TensorOptions().dtype(torch::kInt64));
}

// Immutable tensor views and coordinate indices for one dense layer
struct MTensorDenseInference {
  torch::Tensor weights;
  torch::Tensor biases;
};

// Immutable tensor evaluation metadata for one coupling
struct MTensorCouplingInference {
  torch::Tensor condition_indices;
  torch::Tensor transform_indices;
  torch::Tensor permutation;
  torch::Tensor to_next;
  torch::Tensor to_previous;
  std::vector<MTensorDenseInference> dense_layers;
};

// Immutable tensor evaluation state for one complete flow
struct MFlowInferenceState {
  torch::Tensor parameters;
  std::vector<MTensorCouplingInference> couplings;
};

// Build owned tensor parameters and persistent tensor views for evaluation
std::shared_ptr<const void>
BuildInferenceState(const MSplineFlow &flow,
                    const std::vector<double> &parameters) {
  auto state = std::make_shared<MFlowInferenceState>();
  state->parameters = TensorParameters(parameters, false);
  state->couplings.reserve(flow.Couplings().size());
  for (const auto &index : indices(flow.Couplings())) {
    const MSplineCoupling &coupling = flow.Couplings()[index];
    MTensorCouplingInference coupling_state;
    coupling_state.condition_indices =
        TensorIndices(coupling.condition_indices_);
    coupling_state.transform_indices =
        TensorIndices(coupling.transform_indices_);
    coupling_state.permutation = TensorIndices(coupling.permutation_);
    coupling_state.to_next = TensorIndices(
        index + 1 < flow.Couplings().size()
            ? CouplingTransition(coupling.inverse_permutation_,
                                 flow.Couplings()[index + 1].permutation_)
            : coupling.inverse_permutation_);
    coupling_state.to_previous = TensorIndices(
        index > 0 ? CouplingTransition(coupling.inverse_permutation_,
                                       flow.Couplings()[index - 1].permutation_)
                  : coupling.inverse_permutation_);
    coupling_state.dense_layers.reserve(coupling.dense_shapes_.size());
    for (const MSplineCoupling::DenseShape &shape : coupling.dense_shapes_) {
      const torch::Tensor weights =
          state->parameters
              .slice(0, static_cast<std::int64_t>(shape.weight_offset),
                     static_cast<std::int64_t>(shape.weight_offset +
                                               shape.inputs * shape.outputs))
              .view({static_cast<std::int64_t>(shape.outputs),
                     static_cast<std::int64_t>(shape.inputs)});
      const torch::Tensor biases = state->parameters.slice(
          0, static_cast<std::int64_t>(shape.bias_offset),
          static_cast<std::int64_t>(shape.bias_offset + shape.outputs));
      coupling_state.dense_layers.push_back({weights, biases});
    }
    state->couplings.push_back(std::move(coupling_state));
  }
  return state;
}

// Apply one coordinate permutation to every row of a tensor batch
torch::Tensor TensorPermutation(const torch::Tensor &points,
                                const std::vector<std::size_t> &permutation) {
  return points.index_select(1, TensorIndices(permutation));
}

// Apply one cached coordinate permutation to every tensor row
torch::Tensor TensorPermutation(const torch::Tensor &points,
                                const torch::Tensor &permutation) {
  return points.index_select(1, permutation);
}

// Normalize batched spline logits into positive bins with a prescribed row sum
torch::Tensor TensorNormalizeBins(const torch::Tensor &raw, double minimum) {
  const double available = 1.0 - minimum * static_cast<double>(raw.size(1));
  return minimum + available * torch::softmax(raw, 1);
}

// Normalize batched derivative logits into positive spline slopes
torch::Tensor TensorNormalizeDerivatives(const torch::Tensor &raw,
                                         double minimum) {
  return minimum + torch::softplus(raw);
}

// Gather one entry per row from a rank-two tensor
torch::Tensor TensorGatherRows(const torch::Tensor &values,
                               const torch::Tensor &indices) {
  return values.gather(1, indices.unsqueeze(1)).squeeze(1);
}

// Compute spline bin indices and cumulative knots for one batched coordinate
std::tuple<torch::Tensor, torch::Tensor, torch::Tensor>
TensorSplineBins(const torch::Tensor &coordinate, const torch::Tensor &widths,
                 const torch::Tensor &heights) {
  const torch::Tensor zero =
      torch::zeros({coordinate.size(0), 1}, coordinate.options());
  const torch::Tensor cumulative_widths =
      torch::cat({zero, widths.cumsum(1)}, 1);
  const torch::Tensor cumulative_heights =
      torch::cat({zero, heights.cumsum(1)}, 1);
  const torch::Tensor boundaries =
      cumulative_widths.slice(1, 1, widths.size(1));
  const torch::Tensor bins =
      (coordinate.unsqueeze(1) >= boundaries).sum(1).to(torch::kInt64);
  return {bins, cumulative_widths, cumulative_heights};
}

// Evaluate a batch of normalized rational-quadratic splines forward
std::pair<torch::Tensor, torch::Tensor>
TensorSplineForward(const torch::Tensor &x, const torch::Tensor &widths,
                    const torch::Tensor &heights,
                    const torch::Tensor &derivatives) {
  auto [bin, cumulative_widths, cumulative_heights] =
      TensorSplineBins(x, widths, heights);
  const torch::Tensor x_k = TensorGatherRows(cumulative_widths, bin);
  const torch::Tensor y_k = TensorGatherRows(cumulative_heights, bin);
  const torch::Tensor width = TensorGatherRows(widths, bin);
  const torch::Tensor height = TensorGatherRows(heights, bin);
  const torch::Tensor d0 = TensorGatherRows(derivatives, bin);
  const torch::Tensor d1 = TensorGatherRows(derivatives, bin + 1);
  const torch::Tensor delta = height / width;
  const torch::Tensor xi = (x - x_k) / width;
  const torch::Tensor one_minus_xi = 1.0 - xi;
  const torch::Tensor xi_product = xi * one_minus_xi;
  const torch::Tensor denominator =
      delta + (d1 + d0 - 2.0 * delta) * xi_product;
  const torch::Tensor numerator =
      height * (delta * xi.square() + d0 * xi_product);
  const torch::Tensor derivative_numerator =
      delta.square() * (d1 * xi.square() + 2.0 * delta * xi_product +
                        d0 * one_minus_xi.square());
  const torch::Tensor jacobian = derivative_numerator / denominator.square();
  return {y_k + numerator / denominator, torch::log(jacobian)};
}

// Evaluate a batch of normalized rational-quadratic splines analytically
// backward
std::pair<torch::Tensor, torch::Tensor>
TensorSplineInverse(const torch::Tensor &y, const torch::Tensor &widths,
                    const torch::Tensor &heights,
                    const torch::Tensor &derivatives) {
  auto [bin, cumulative_heights, cumulative_widths] =
      TensorSplineBins(y, heights, widths);
  const torch::Tensor y_k = TensorGatherRows(cumulative_heights, bin);
  const torch::Tensor x_k = TensorGatherRows(cumulative_widths, bin);
  const torch::Tensor width = TensorGatherRows(widths, bin);
  const torch::Tensor height = TensorGatherRows(heights, bin);
  const torch::Tensor d0 = TensorGatherRows(derivatives, bin);
  const torch::Tensor d1 = TensorGatherRows(derivatives, bin + 1);
  const torch::Tensor delta = height / width;
  const torch::Tensor y_delta = y - y_k;
  const torch::Tensor slope_term = d0 + d1 - 2.0 * delta;
  const torch::Tensor a = y_delta * slope_term + height * (delta - d0);
  const torch::Tensor b = height * d0 - y_delta * slope_term;
  const torch::Tensor c = -delta * y_delta;
  const torch::Tensor coefficient_scale = torch::maximum(
      torch::abs(a), torch::maximum(torch::abs(b), torch::abs(c)));
  if (torch::any(torch::logical_or(
                     torch::logical_not(torch::isfinite(coefficient_scale)),
                     coefficient_scale <= 0.0))
          .item<bool>()) {
    throw std::runtime_error(
        "TensorSplineInverse: invalid inverse coefficients");
  }
  const torch::Tensor normalized_a = a / coefficient_scale;
  const torch::Tensor normalized_b = b / coefficient_scale;
  const torch::Tensor normalized_c = c / coefficient_scale;
  const torch::Tensor product = 4.0 * normalized_a * normalized_c;
  torch::Tensor discriminant = normalized_b.square() - product;
  const torch::Tensor discriminant_scale =
      normalized_b.square() + torch::abs(product);
  const double relative_tolerance =
      SplineInverseRoundoffFactor() * std::numeric_limits<double>::epsilon();
  if (torch::any(discriminant < -relative_tolerance * discriminant_scale)
          .item<bool>()) {
    throw std::runtime_error(
        "TensorSplineInverse: negative quadratic discriminant");
  }
  discriminant = torch::where(discriminant < 0.0,
                              torch::zeros_like(discriminant), discriminant);
  const torch::Tensor root = torch::sqrt(discriminant);
  const torch::Tensor negative_b = normalized_b < 0.0;
  const torch::Tensor quadratic_xi =
      torch::where(negative_b, -normalized_b + root, 2.0 * normalized_c) /
      torch::where(negative_b, 2.0 * normalized_a, -normalized_b - root);
  const torch::Tensor xi = torch::clamp(quadratic_xi, 0.0, 1.0);
  const torch::Tensor one_minus_xi = 1.0 - xi;
  const torch::Tensor xi_product = xi * one_minus_xi;
  const torch::Tensor denominator =
      delta + (d1 + d0 - 2.0 * delta) * xi_product;
  const torch::Tensor derivative_numerator =
      delta.square() * (d1 * xi.square() + 2.0 * delta * xi_product +
                        d0 * one_minus_xi.square());
  const torch::Tensor jacobian = derivative_numerator / denominator.square();
  return {x_k + xi * width, -torch::log(jacobian)};
}

// Evaluate one conditioner MLP for a complete tensor batch
torch::Tensor TensorConditioner(const MSplineCoupling &coupling,
                                const torch::Tensor &points,
                                const torch::Tensor &parameters,
                                const MTensorCouplingInference *inference) {
  torch::Tensor values =
      inference != nullptr
          ? points.index_select(1, inference->condition_indices)
          : points.index_select(1, TensorIndices(coupling.condition_indices_));
  for (const auto &layer : indices(coupling.dense_shapes_)) {
    const MSplineCoupling::DenseShape &shape = coupling.dense_shapes_[layer];
    const torch::Tensor weights =
        inference != nullptr
            ? inference->dense_layers[layer].weights
            : parameters
                  .slice(
                      0, static_cast<std::int64_t>(shape.weight_offset),
                      static_cast<std::int64_t>(shape.weight_offset +
                                                shape.inputs * shape.outputs))
                  .view({static_cast<std::int64_t>(shape.outputs),
                         static_cast<std::int64_t>(shape.inputs)});
    const torch::Tensor biases =
        inference != nullptr
            ? inference->dense_layers[layer].biases
            : parameters.slice(
                  0, static_cast<std::int64_t>(shape.bias_offset),
                  static_cast<std::int64_t>(shape.bias_offset + shape.outputs));
    values = torch::matmul(values, weights.transpose(0, 1)) + biases;
    if (layer + 1 < coupling.dense_shapes_.size()) {
      if (coupling.activation_ == FlowActivation::Relu) {
        values = torch::relu(values);
      } else if (coupling.activation_ == FlowActivation::Silu) {
        values = torch::silu(values);
      } else if (coupling.activation_ == FlowActivation::Tanh) {
        values = torch::tanh(values);
      } else {
        throw std::invalid_argument(
            "TensorConditioner: unknown activation function");
      }
    }
  }
  return values;
}

// Apply one tensor spline coupling with flattened coordinate vectorization
std::pair<torch::Tensor, torch::Tensor>
TensorCoupling(const MSplineCoupling &coupling, const torch::Tensor &input,
               const torch::Tensor &parameters,
               const MTensorCouplingInference *inference, bool inverse) {
  const torch::Tensor conditioned =
      TensorConditioner(coupling, input, parameters, inference);
  const std::int64_t transforms =
      static_cast<std::int64_t>(coupling.transform_indices_.size());
  const std::int64_t bins = static_cast<std::int64_t>(coupling.bins_);
  const std::int64_t stride = 3 * bins + 1;
  const torch::Tensor raw =
      conditioned.reshape({input.size(0) * transforms, stride});
  const torch::Tensor widths =
      TensorNormalizeBins(raw.slice(1, 0, bins), coupling.min_bin_);
  const torch::Tensor heights =
      TensorNormalizeBins(raw.slice(1, bins, 2 * bins), coupling.min_bin_);
  const torch::Tensor derivatives = TensorNormalizeDerivatives(
      raw.slice(1, 2 * bins, 3 * bins + 1), coupling.min_derivative_);
  const torch::Tensor transform_indices =
      inference != nullptr ? inference->transform_indices
                           : TensorIndices(coupling.transform_indices_);
  const torch::Tensor coordinates =
      input.index_select(1, transform_indices).reshape({-1});
  auto transformed =
      inverse ? TensorSplineInverse(coordinates, widths, heights, derivatives)
              : TensorSplineForward(coordinates, widths, heights, derivatives);
  const torch::Tensor values =
      transformed.first.reshape({input.size(0), transforms});
  const torch::Tensor log_jacobian =
      transformed.second.reshape({input.size(0), transforms}).sum(1);
  return {input.index_copy(1, transform_indices, values), log_jacobian};
}

// Apply every spline coupling forward to a complete tensor batch
std::pair<torch::Tensor, torch::Tensor>
TensorFlowForward(const MSplineFlow &flow, const torch::Tensor &base,
                  const torch::Tensor &parameters,
                  const MFlowInferenceState *inference) {
  torch::Tensor log_jacobian = torch::zeros({base.size(0)}, base.options());
  if (inference != nullptr) {
    torch::Tensor points =
        TensorPermutation(base, inference->couplings.front().permutation);
    for (const auto &index : indices(flow.Couplings())) {
      const MSplineCoupling &coupling = flow.Couplings()[index];
      const MTensorCouplingInference &coupling_state =
          inference->couplings[index];
      auto transformed =
          TensorCoupling(coupling, points, parameters, &coupling_state, false);
      points = TensorPermutation(transformed.first, coupling_state.to_next);
      log_jacobian = log_jacobian + transformed.second;
    }
    return {points, log_jacobian};
  }

  torch::Tensor points = base;
  for (const auto &index : indices(flow.Couplings())) {
    const MSplineCoupling &coupling = flow.Couplings()[index];
    const torch::Tensor permuted =
        TensorPermutation(points, coupling.permutation_);
    auto transformed =
        TensorCoupling(coupling, permuted, parameters, nullptr, false);
    points =
        TensorPermutation(transformed.first, coupling.inverse_permutation_);
    log_jacobian = log_jacobian + transformed.second;
  }
  return {points, log_jacobian};
}

// Apply every spline coupling backward to a complete tensor batch
torch::Tensor TensorFlowLogDensity(const MSplineFlow &flow,
                                   const torch::Tensor &points,
                                   const torch::Tensor &parameters,
                                   const MFlowInferenceState *inference) {
  torch::Tensor log_density = torch::zeros({points.size(0)}, points.options());
  if (inference != nullptr) {
    torch::Tensor base =
        TensorPermutation(points, inference->couplings.back().permutation);
    for (std::size_t reverse = flow.Couplings().size(); reverse-- > 0;) {
      const MSplineCoupling &coupling = flow.Couplings()[reverse];
      const MTensorCouplingInference &coupling_state =
          inference->couplings[reverse];
      auto transformed =
          TensorCoupling(coupling, base, parameters, &coupling_state, true);
      base = TensorPermutation(transformed.first, coupling_state.to_previous);
      log_density = log_density + transformed.second;
    }
    return log_density;
  }

  torch::Tensor base = points;
  for (std::size_t reverse = flow.Couplings().size(); reverse-- > 0;) {
    const MSplineCoupling &coupling = flow.Couplings()[reverse];
    const torch::Tensor permuted =
        TensorPermutation(base, coupling.permutation_);
    auto transformed =
        TensorCoupling(coupling, permuted, parameters, nullptr, true);
    base = TensorPermutation(transformed.first, coupling.inverse_permutation_);
    log_density = log_density + transformed.second;
  }
  return log_density;
}

// Copy one contiguous rank-one CPU tensor into a standard vector
std::vector<double> TensorVector(const torch::Tensor &input) {
  const torch::Tensor contiguous = input.contiguous();
  std::vector<double> output(static_cast<std::size_t>(contiguous.numel()));
  std::copy_n(contiguous.data_ptr<double>(), output.size(), output.begin());
  return output;
}

#endif

// Validate all numerical parameters for a phase-space dimension
void MNeuroJacConfig::Validate(std::size_t dimension) const {
  if (dimension == 0) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: dimension must be positive");
  }
  if (flow_components == 0 || layers == 0 || bins < 2 || buffer_size == 0 ||
      batch_size == 0 || rounds == 0 || validation_size == 0 || epochs == 0) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: counts must be positive and bins at least two");
  }
  const auto        limit   = std::numeric_limits<std::size_t>::max();
  const std::size_t batches = 1 + (buffer_size - 1) / batch_size;
  if (rounds > limit / epochs || rounds * epochs > limit / batches) {
    throw std::invalid_argument("NUMERICS_NEUROJAC: optimizer schedule overflows");
  }
  if (std::any_of(hidden.begin(), hidden.end(),
                  [](std::size_t width) { return width == 0; })) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: hidden layer widths must be positive");
  }
  for (double value : {learning_rate, final_learning_rate, validation_min_delta, uniform_mix, replay_fraction, min_bin,
                       min_derivative, density_tolerance, gradient_clip, adamw_beta1, adamw_beta2, adamw_epsilon,
                       adamw_weight_decay, adamw_decay_biases}) {
    if (!std::isfinite(value)) {
      throw std::invalid_argument(
          "NUMERICS_NEUROJAC: all floating-point parameters must be finite");
    }
  }
  if (!(learning_rate > 0.0) || !(final_learning_rate > 0.0) || !(validation_min_delta >= 0.0) ||
      !(uniform_mix > 0.0 && uniform_mix < 1.0) || !(replay_fraction >= 0.0 && replay_fraction < 1.0) ||
      !(min_bin > 0.0 && min_bin * static_cast<double>(bins) < 1.0) ||
      !(min_derivative > 0.0 && min_derivative < 1.0) || !(density_tolerance > 0.0) || !(gradient_clip > 0.0)) {
    throw std::invalid_argument("NUMERICS_NEUROJAC: invalid learning rate, "
                                "mixture, replay fraction or spline bound");
  }
  if (replay_fraction > 0.0 && replay_capacity == 0) {
    throw std::invalid_argument("NUMERICS_NEUROJAC: replay_cap must be "
                                "positive when replay is enabled");
  }
  if (!vegas_init && vegas_sampler_rounds != 0) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: vegas_sampler_rounds must be zero when VEGAS "
        "initialization is disabled");
  }
  if (vegas_spline_init && !vegas_init) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: vegas_spline_init requires vegas_init");
  }
  if (vegas_sampler_rounds > rounds) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: vegas_sampler_rounds cannot exceed rounds");
  }
  if (vegas_init) {
    const std::size_t vegas_count_limit =
        std::numeric_limits<unsigned int>::max();
    if (vegas_ncall == 0 || vegas_rounds == 0 ||
        vegas_ncall > vegas_count_limit || vegas_rounds > vegas_count_limit) {
      throw std::invalid_argument(
          "NUMERICS_NEUROJAC: vegas_ncall and vegas_rounds must be positive "
          "unsigned int values when VEGAS initialization is enabled");
    }
    if (vegas_ncall > std::numeric_limits<std::size_t>::max() / vegas_rounds) {
      throw std::invalid_argument(
          "NUMERICS_NEUROJAC: vegas_ncall times vegas_rounds overflows");
    }
  }
  if (!(adamw_beta1 >= 0.0 && adamw_beta1 < 1.0) ||
      !(adamw_beta2 >= 0.0 && adamw_beta2 < 1.0) || !(adamw_epsilon > 0.0) ||
      !(adamw_weight_decay >= 0.0) || !(adamw_decay_biases >= 0.0)) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: invalid AdamW hyperparameters");
  }
  switch (activation) {
  case FlowActivation::Relu:
  case FlowActivation::Silu:
  case FlowActivation::Tanh:
    break;
  default:
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: invalid activation function");
  }
  switch (lr_schedule) {
  case FlowLRSchedule::None:
  case FlowLRSchedule::Cosine:
  case FlowLRSchedule::Exponential:
    break;
  default:
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: invalid learning-rate schedule");
  }
  switch (permutation) {
  case FlowPermutation::None:
  case FlowPermutation::Interleave:
  case FlowPermutation::Random:
    break;
  default:
    throw std::invalid_argument("NUMERICS_NEUROJAC: invalid permutation mode");
  }
  loss.Validate();
}

// Evaluate the forward spline using normalized positive knot parameters
MSplineResult MSpline1D::Forward(double x, const std::vector<double> &widths,
                                 const std::vector<double> &heights,
                                 const std::vector<double> &derivatives) {
  ValidateSplineParameters(widths, heights, derivatives);
  const auto result = SplineForward(x, widths, heights, derivatives);
  return {result.first, result.second};
}

// Evaluate the analytic inverse spline using normalized positive knot
// parameters
MSplineResult MSpline1D::Inverse(double y, const std::vector<double> &widths,
                                 const std::vector<double> &heights,
                                 const std::vector<double> &derivatives) {
  ValidateSplineParameters(widths, heights, derivatives);
  const auto result = SplineInverse(y, widths, heights, derivatives);
  return {result.first, result.second};
}

// Initialize an identity flow with deterministic conditioner parameters
void MSplineFlow::Configure(std::size_t dimension,
                            const MNeuroJacConfig &config,
                            std::uint32_t initialization_seed) {
  config.Validate(dimension);
  dimension_ = dimension;
  batch_size_ = config.batch_size;
  couplings_.clear();
  parameters_.clear();
  weight_decay_mask_.clear();
  inference_state_.reset();

  std::mt19937 generator(initialization_seed);
  std::normal_distribution<double> normal(0.0, 1.0);
  const double derivative_bias = InverseSoftplus(1.0 - config.min_derivative);

  for (std::size_t layer = 0; layer < config.layers; ++layer) {
    MSplineCoupling coupling;
    coupling.dimension_ = dimension;
    coupling.bins_ = config.bins;
    coupling.min_bin_ = config.min_bin;
    coupling.min_derivative_ = config.min_derivative;
    coupling.activation_ = config.activation;
    coupling.parameter_offset_ = parameters_.size();

    for (std::size_t coordinate = 0; coordinate < dimension; ++coordinate) {
      if (dimension > 1 && (coordinate + layer) % 2 == 0) {
        coupling.condition_indices_.push_back(coordinate);
      } else {
        coupling.transform_indices_.push_back(coordinate);
      }
    }
    if (coupling.transform_indices_.empty()) {
      coupling.transform_indices_.push_back(coupling.condition_indices_.back());
      coupling.condition_indices_.pop_back();
    }

    const std::size_t output_size =
        coupling.transform_indices_.size() * (3 * config.bins + 1);
    std::vector<std::size_t> network_widths;
    network_widths.push_back(coupling.condition_indices_.size());
    if (!coupling.condition_indices_.empty()) {
      network_widths.insert(network_widths.end(), config.hidden.begin(),
                            config.hidden.end());
    }
    network_widths.push_back(output_size);

    for (std::size_t dense = 0; dense + 1 < network_widths.size(); ++dense) {
      MSplineCoupling::DenseShape shape;
      shape.inputs = network_widths[dense];
      shape.outputs = network_widths[dense + 1];
      shape.weight_offset = parameters_.size();
      parameters_.resize(parameters_.size() + shape.inputs * shape.outputs,
                         0.0);
      weight_decay_mask_.resize(parameters_.size(), 1);
      shape.bias_offset = parameters_.size();
      parameters_.resize(parameters_.size() + shape.outputs, 0.0);
      weight_decay_mask_.resize(parameters_.size(), 0);

      const bool output_layer = dense + 2 == network_widths.size();
      if (!output_layer) {
        const double scale =
            config.activation == FlowActivation::Tanh
                ? std::sqrt(2.0 /
                            static_cast<double>(shape.inputs + shape.outputs))
                : std::sqrt(2.0 / static_cast<double>(shape.inputs));
        for (std::size_t i = 0; i < shape.inputs * shape.outputs; ++i) {
          parameters_[shape.weight_offset + i] = scale * normal(generator);
        }
      }
      coupling.dense_shapes_.push_back(shape);
    }

    const MSplineCoupling::DenseShape &output_shape =
        coupling.dense_shapes_.back();
    const std::size_t stride = 3 * config.bins + 1;
    for (const auto &coordinate : indices(coupling.transform_indices_)) {
      for (std::size_t knot = 0; knot <= config.bins; ++knot) {
        parameters_[output_shape.bias_offset + coordinate * stride +
                    2 * config.bins + knot] = derivative_bias;
      }
    }

    coupling.parameter_count_ = parameters_.size() - coupling.parameter_offset_;
    coupling.permutation_ = IdentityPermutation(dimension);
    coupling.inverse_permutation_ = coupling.permutation_;
    couplings_.push_back(std::move(coupling));
  }

  std::mt19937 permutation_generator(initialization_seed ^ 0xa341316cU);
  std::vector<std::size_t> transform_counts(dimension, 0);
  for (const auto &layer : indices(couplings_)) {
    const std::vector<std::size_t> permutation = BalancedPermutation(
        dimension, couplings_[layer].transform_indices_, config.permutation,
        layer, transform_counts, permutation_generator);
    const std::vector<std::size_t> inverse = InversePermutation(permutation);
    couplings_[layer].permutation_ = permutation;
    couplings_[layer].inverse_permutation_ = inverse;
  }
  PrepareInference();
}

// Transform one uniform base point and return its exact output density
MFlowSample MSplineFlow::Forward(const std::vector<double> &base) const {
  ValidateUnitPoint(base, dimension_, "MSplineFlow::Forward");
  const auto transformed =
      MSplineEvaluator::FlowForward(*this, base, parameters_);
  return {transformed.first, -transformed.second};
}

// Transform row-major base points with one vectorized evaluation
MFlowSampleBatch
MSplineFlow::ForwardBatch(const std::vector<double> &bases) const {
  MFlowSampleBatch samples;
  samples.dimension = dimension_;
  if (bases.empty()) {
    return samples;
  }
  const std::size_t events =
      ValidateFlatUnitPoints(bases, dimension_, "MSplineFlow::ForwardBatch");
#ifdef GRANIITTI_USE_LIBTORCH
  c10::InferenceMode guard;
#endif
  std::shared_ptr<const void> inference = inference_state_;
  if (inference == nullptr) {
    inference = BuildInferenceState(*this, parameters_);
  }
  const auto *state = static_cast<const MFlowInferenceState *>(inference.get());
#ifdef GRANIITTI_USE_LIBTORCH
  const torch::Tensor base = TensorFlatPoints(bases, events, dimension_);
  auto transformed = TensorFlowForward(*this, base, state->parameters, state);
  samples.points = TensorVector(transformed.first);
  samples.log_densities = TensorVector(-transformed.second);
#else
  samples.points.resize(bases.size());
  samples.log_densities.resize(events);
  for (std::size_t begin = 0; begin < events; begin += batch_size_) {
    const std::size_t count = std::min(batch_size_, events - begin);
    const Eigen::Map<const MNeuroJacBatch> base(bases.data() + begin * dimension_,
                                               static_cast<Eigen::Index>(count), static_cast<Eigen::Index>(dimension_));
    const auto transformed = StandaloneFlowForward(*this, base, parameters_, state);
    std::copy_n(transformed.first.data(), count * dimension_, samples.points.data() + begin * dimension_);
    for (std::size_t event = 0; event < count; ++event) {
      samples.log_densities[begin + event] = -transformed.second[static_cast<Eigen::Index>(event)];
    }
  }
#endif
  return samples;
}

// Evaluate the exact flow density at one unit-hypercube point
double MSplineFlow::LogDensity(const std::vector<double> &point) const {
  ValidateUnitPoint(point, dimension_, "MSplineFlow::LogDensity");
  return MSplineEvaluator::FlowLogDensity(*this, point, parameters_);
}

// Evaluate exact flow densities for row-major unit-hypercube points
std::vector<double>
MSplineFlow::LogDensityBatch(const std::vector<double> &points) const {
  if (points.empty()) {
    return {};
  }
  const std::size_t events = ValidateFlatUnitPoints(
      points, dimension_, "MSplineFlow::LogDensityBatch");
#ifdef GRANIITTI_USE_LIBTORCH
  c10::InferenceMode guard;
#endif
  std::shared_ptr<const void> inference = inference_state_;
  if (inference == nullptr) {
    inference = BuildInferenceState(*this, parameters_);
  }
  const auto *state = static_cast<const MFlowInferenceState *>(inference.get());
#ifdef GRANIITTI_USE_LIBTORCH
  const torch::Tensor point_tensor =
      TensorFlatPoints(points, events, dimension_);
  return TensorVector(
      TensorFlowLogDensity(*this, point_tensor, state->parameters, state));
#else
  std::vector<double> log_densities(events);
  for (std::size_t begin = 0; begin < events; begin += batch_size_) {
    const std::size_t count = std::min(batch_size_, events - begin);
    const Eigen::Map<const MNeuroJacBatch> point_batch(points.data() + begin * dimension_,
                                                      static_cast<Eigen::Index>(count), static_cast<Eigen::Index>(dimension_));
    const auto batch_density = StandaloneFlowLogDensity(*this, point_batch, parameters_, nullptr, state);
    std::copy_n(batch_density.data(), count, log_densities.begin() + begin);
  }
  return log_densities;
#endif
}

// Build immutable backend-specific metadata for repeated evaluation
void MSplineFlow::PrepareInference() {
#ifdef GRANIITTI_USE_LIBTORCH
  c10::InferenceMode guard;
#endif
  inference_state_ = BuildInferenceState(*this, parameters_);
}

// Replace the flat trainable parameters after validation
void MSplineFlow::SetParameters(std::vector<double> parameters) {
  if (parameters.size() != parameters_.size()) {
    throw std::invalid_argument(
        "MSplineFlow::SetParameters: parameter size mismatch");
  }
  if (!gra::AllFinite(parameters)) {
    throw std::invalid_argument(
        "MSplineFlow::SetParameters: parameters must be finite");
  }
  parameters_ = std::move(parameters);
  PrepareInference();
}

// Apply one validated additive optimizer step and decoupled decay in place
void MSplineFlow::ApplyParameterStep(const std::vector<double> &steps,
                                     double matrix_decay_factor,
                                     double bias_decay_factor) {
  if (steps.size() != parameters_.size()) {
    throw std::invalid_argument(
        "MSplineFlow::ApplyParameterStep: parameter size mismatch");
  }
  if (!std::isfinite(matrix_decay_factor) ||
      !std::isfinite(bias_decay_factor) || matrix_decay_factor < 0.0 ||
      bias_decay_factor < 0.0) {
    throw std::invalid_argument(
        "MSplineFlow::ApplyParameterStep: invalid decay factor");
  }
  for (const auto &i : indices(parameters_)) {
    if (!std::isfinite(steps[i])) {
      throw std::runtime_error(
          "MSplineFlow::ApplyParameterStep: non-finite optimizer step");
    }
    const double decay =
        weight_decay_mask_[i] != 0 ? matrix_decay_factor : bias_decay_factor;
    if (!std::isfinite(parameters_[i] * (1.0 - decay) - steps[i])) {
      throw std::runtime_error(
          "MSplineFlow::ApplyParameterStep: non-finite parameter update");
    }
  }
  inference_state_.reset();
  for (const auto &i : indices(parameters_)) {
    const double decay =
        weight_decay_mask_[i] != 0 ? matrix_decay_factor : bias_decay_factor;
    parameters_[i] = parameters_[i] * (1.0 - decay) - steps[i];
  }
}

// Initialize one separable spline map from an adapted VEGAS grid
void MSplineFlow::InitializeVegasGrid(
    const std::vector<std::vector<double>> &quantile_edges) {
  if (dimension_ == 0 || couplings_.empty()) {
    throw std::logic_error(
        "MSplineFlow::InitializeVegasGrid: flow is not configured");
  }
  if (quantile_edges.size() != dimension_) {
    throw std::invalid_argument(
        "MSplineFlow::InitializeVegasGrid: grid dimension mismatch");
  }

  const std::size_t bins = couplings_.front().bins_;
  std::vector<std::vector<double>> resampled_heights(dimension_);
  for (std::size_t coordinate = 0; coordinate < dimension_; ++coordinate) {
    resampled_heights[coordinate] =
        ResampleVegasQuantileHeights(quantile_edges[coordinate], bins);
  }

  std::vector<std::size_t> assigned_coupling(dimension_, couplings_.size());
  std::vector<std::size_t> assigned_output(dimension_, 0);
  for (const auto &layer : indices(couplings_)) {
    const MSplineCoupling &coupling = couplings_[layer];
    for (const auto &output : indices(coupling.transform_indices_)) {
      const std::size_t coordinate =
          coupling.permutation_[coupling.transform_indices_[output]];
      if (assigned_coupling[coordinate] == couplings_.size()) {
        assigned_coupling[coordinate] = layer;
        assigned_output[coordinate] = output;
      }
    }
  }

  if (std::find(assigned_coupling.begin(), assigned_coupling.end(),
                couplings_.size()) != assigned_coupling.end()) {
    throw std::logic_error(
        "MSplineFlow::InitializeVegasGrid: flow does not transform every "
        "coordinate");
  }

  for (std::size_t coordinate = 0; coordinate < dimension_; ++coordinate) {
    const MSplineCoupling &coupling = couplings_[assigned_coupling[coordinate]];
    const MSplineCoupling::DenseShape &output_shape =
        coupling.dense_shapes_.back();
    const std::size_t stride = 3 * bins + 1;
    const std::size_t row_offset = assigned_output[coordinate] * stride;
    const std::vector<double> heights =
        RegularizeVegasWidths(resampled_heights[coordinate], coupling.min_bin_);
    const std::vector<double> derivatives =
        VegasKnotDerivatives(heights, coupling.min_derivative_);
    const double available =
        1.0 - coupling.min_bin_ * static_cast<double>(bins);

    for (std::size_t row = 0; row < stride; ++row) {
      for (std::size_t input = 0; input < output_shape.inputs; ++input) {
        parameters_[output_shape.weight_offset +
                    (row_offset + row) * output_shape.inputs + input] = 0.0;
      }
    }
    for (std::size_t bin = 0; bin < bins; ++bin) {
      parameters_[output_shape.bias_offset + row_offset + bin] = 0.0;
      parameters_[output_shape.bias_offset + row_offset + bins + bin] =
          std::log((heights[bin] - coupling.min_bin_) / available);
    }
    for (std::size_t knot = 0; knot <= bins; ++knot) {
      parameters_[output_shape.bias_offset + row_offset + 2 * bins + knot] =
          InverseSoftplus(derivatives[knot] - coupling.min_derivative_);
    }
  }
  PrepareInference();
}

// Evaluate one mini-batch objective and its reverse-mode gradient
double MSplineFlow::LossGradient(const std::vector<MFlowTrainingEvent> &events,
                                 const std::vector<double> &coefficients,
                                 const std::vector<std::size_t> &indices,
                                 std::size_t population_size,
                                 const MFlowObjective &loss_function,
                                 double uniform_mix,
                                 std::vector<double> &gradient) const {
  if (events.size() != coefficients.size() || indices.empty() ||
      population_size < indices.size() || population_size > events.size()) {
    throw std::invalid_argument(
        "MSplineFlow::LossGradient: invalid mini-batch input");
  }

  double coefficient_sum = 0.0;
  for (std::size_t index : indices) {
    coefficient_sum += coefficients.at(index);
  }
  if (!(coefficient_sum > 0.0)) {
    gradient.assign(parameters_.size(), 0.0);
    return 0.0;
  }
  const double batch_scale = static_cast<double>(population_size) /
                             static_cast<double>(indices.size());

#ifdef GRANIITTI_USE_LIBTORCH
  std::vector<double> batch_points;
  std::vector<double> batch_coefficients;
  batch_points.reserve(indices.size() * dimension_);
  batch_coefficients.reserve(indices.size());
  for (std::size_t index : indices) {
    const std::vector<double> &point = events.at(index).point;
    ValidateUnitPoint(point, dimension_, "MSplineFlow::LossGradient");
    batch_points.insert(batch_points.end(), point.begin(), point.end());
    batch_coefficients.push_back(coefficients.at(index) * batch_scale);
  }

  const torch::Tensor point_tensor =
      TensorFlatPoints(batch_points, indices.size(), dimension_);
  torch::Tensor parameter_tensor = TensorParameters(parameters_, true);
  const torch::Tensor flow_log_density =
      TensorFlowLogDensity(*this, point_tensor, parameter_tensor, nullptr);
  const std::vector<double> flow_log_densities = TensorVector(flow_log_density);
  std::vector<double> flow_log_density_derivatives(indices.size(), 0.0);
  double objective_value = 0.0;
  for (const auto &event : gra::aux::indices(indices)) {
    const double mixed_log_density =
        MixedLogDensity(flow_log_densities[event], uniform_mix);
    const double flow_log_fraction = std::log1p(-uniform_mix) +
                                     flow_log_densities[event] -
                                     mixed_log_density;
    const MFlowLossEvaluation evaluation =
        loss_function.Evaluate(batch_coefficients[event], mixed_log_density,
                               std::exp(flow_log_fraction));
    objective_value += evaluation.value;
    flow_log_density_derivatives[event] =
        evaluation.flow_log_density_derivative;
  }
  torch::Tensor derivative_tensor = torch::empty(
      {static_cast<std::int64_t>(indices.size())}, NeuroJacTensorOptions());
  std::copy(flow_log_density_derivatives.begin(),
            flow_log_density_derivatives.end(),
            derivative_tensor.data_ptr<double>());
  (flow_log_density * derivative_tensor).sum().backward();
  if (!parameter_tensor.grad().defined()) {
    throw std::runtime_error(
        "MSplineFlow::LossGradient: autograd produced no gradient");
  }
  gradient = TensorVector(parameter_tensor.grad());
  if (!std::isfinite(objective_value)) {
    throw std::runtime_error("MSplineFlow::LossGradient: non-finite objective");
  }
  return objective_value;
#else
  std::vector<double> batch_points;
  batch_points.reserve(indices.size() * dimension_);
  MNeuroJacBatchVector batch_coefficients(
      static_cast<Eigen::Index>(indices.size()));
  for (const auto &batch : gra::aux::indices(indices)) {
    const std::size_t index = indices[batch];
    const std::vector<double> &point = events.at(index).point;
    ValidateUnitPoint(point, dimension_, "MSplineFlow::LossGradient");
    batch_points.insert(batch_points.end(), point.begin(), point.end());
    batch_coefficients[static_cast<Eigen::Index>(batch)] =
        coefficients.at(index) * batch_scale;
  }

  const Eigen::Map<const MNeuroJacBatch> point_batch(
      batch_points.data(), static_cast<Eigen::Index>(indices.size()),
      static_cast<Eigen::Index>(dimension_));
  std::vector<MStandaloneCouplingCache> caches;
  const MNeuroJacBatchVector flow_log_density = StandaloneFlowLogDensity(
      *this, point_batch, parameters_, &caches, nullptr);
  MNeuroJacBatchVector log_density_adjoint(flow_log_density.size());
  double objective_value = 0.0;
  for (const auto &event : gra::aux::indices(flow_log_density)) {
    if (!std::isfinite(flow_log_density[event])) {
      throw std::runtime_error(
          "MSplineFlow::LossGradient: non-finite flow log density");
    }
    const double mixed_log_density =
        MixedLogDensity(flow_log_density[event], uniform_mix);
    const double flow_log_fraction =
        std::log1p(-uniform_mix) + flow_log_density[event] - mixed_log_density;
    const double coefficient = batch_coefficients[event];
    const double flow_fraction = std::exp(flow_log_fraction);
    const MFlowLossEvaluation evaluation =
        loss_function.Evaluate(coefficient, mixed_log_density, flow_fraction);
    objective_value += evaluation.value;
    log_density_adjoint[event] = evaluation.flow_log_density_derivative;
  }
  if (!std::isfinite(objective_value)) {
    throw std::runtime_error("MSplineFlow::LossGradient: non-finite objective");
  }
  gradient.assign(parameters_.size(), 0.0);
  StandaloneFlowBackward(*this, caches, log_density_adjoint, parameters_,
                         gradient);
  return objective_value;
#endif
}

// Configure independent spline components with uniform initial weights
void MSplineFlowMixture::Configure(std::size_t dimension,
                                   const MNeuroJacConfig &config,
                                   std::uint32_t initialization_seed) {
  config.Validate(dimension);
  flows_.clear();
  flows_.resize(config.flow_components);
  for (const auto &component : indices(flows_)) {
    const std::uint32_t component_seed =
        initialization_seed ^
        (0x9e3779b9U * static_cast<std::uint32_t>(component));
    flows_[component].Configure(dimension, config, component_seed);
  }
  weight_logits_.assign(flows_.size() - 1, 0.0);
}

// Compute the common phase space dimension
std::size_t MSplineFlowMixture::Dimension() const {
  return flows_.empty() ? 0 : flows_.front().Dimension();
}

// Compute one immutable spline component
const MSplineFlow &MSplineFlowMixture::Component(std::size_t index) const {
  return flows_.at(index);
}

// Append one independently initialized spline component
void MSplineFlowMixture::AddComponent(const MNeuroJacConfig &config,
                                      std::uint32_t initialization_seed,
                                      double weight) {
  if (flows_.empty() || flows_.size() >= config.flow_components) {
    throw std::logic_error(
        "MSplineFlowMixture::AddComponent: component cannot be appended");
  }
  MSplineFlow flow;
  const std::size_t component = flows_.size();
  const std::uint32_t component_seed =
      initialization_seed ^
      (0x9e3779b9U * static_cast<std::uint32_t>(component));
  flow.Configure(Dimension(), config, component_seed);
  flows_.push_back(std::move(flow));
  weight_logits_.resize(flows_.size() - 1);
  SetLastWeight(weight);
}

// Set the last component weight without changing earlier relative weights
void MSplineFlowMixture::SetLastWeight(double weight) {
  if (flows_.size() < 2 || !std::isfinite(weight) || !(weight > 0.0) ||
      !(weight < 1.0)) {
    throw std::invalid_argument(
        "MSplineFlowMixture::SetLastWeight: invalid component weight");
  }
  const double maximum =
      *std::max_element(weight_logits_.begin(), weight_logits_.end());
  double old_sum = 0.0;
  for (double logit : weight_logits_) {
    old_sum += std::exp(logit - maximum);
  }
  const double log_sum = std::log(old_sum);
  const double log_ratio = std::log1p(-weight) - std::log(weight);
  for (const auto &component : indices(weight_logits_)) {
    weight_logits_[component] =
        (weight_logits_[component] - maximum) - log_sum + log_ratio;
  }
}

// Compute normalized non-negative component weights
std::vector<double> MSplineFlowMixture::Weights() const {
  std::vector<double> weights = LogWeights();
  for (double &weight : weights) {
    weight = std::exp(weight);
  }
  return weights;
}

// Compute normalized logarithmic component weights without underflow
std::vector<double> MSplineFlowMixture::LogWeights() const {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::LogWeights: mixture is not configured");
  }
  if (flows_.size() == 1) {
    return {0.0};
  }

  double maximum = 0.0;
  for (double logit : weight_logits_) {
    maximum = std::max(maximum, logit);
  }
  double denominator = std::exp(-maximum);
  for (double logit : weight_logits_) {
    denominator += std::exp(logit - maximum);
  }
  const double log_sum = std::log(denominator);
  std::vector<double> log_weights(flows_.size(), -maximum - log_sum);
  for (const auto &component : indices(weight_logits_)) {
    log_weights[component] =
        (weight_logits_[component] - maximum) - log_sum;
  }
  return log_weights;
}

// Compute all flow parameters followed by independent weight logits
std::vector<double> MSplineFlowMixture::Parameters() const {
  if (flows_.empty()) {
    return {};
  }
  const std::size_t flow_parameter_count = flows_.front().Parameters().size();
  std::vector<double> parameters;
  parameters.reserve(flows_.size() * flow_parameter_count +
                     weight_logits_.size());
  for (const MSplineFlow &flow : flows_) {
    const std::vector<double> &component_parameters = flow.Parameters();
    parameters.insert(parameters.end(), component_parameters.begin(),
                      component_parameters.end());
  }
  parameters.insert(parameters.end(), weight_logits_.begin(),
                    weight_logits_.end());
  return parameters;
}

// Replace all flow parameters and independent weight logits
void MSplineFlowMixture::SetParameters(const std::vector<double> &parameters) {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::SetParameters: mixture is not configured");
  }
  const std::size_t flow_parameter_count = flows_.front().Parameters().size();
  const std::size_t expected =
      flows_.size() * flow_parameter_count + weight_logits_.size();
  if (parameters.size() != expected) {
    throw std::invalid_argument(
        "MSplineFlowMixture::SetParameters: parameter size mismatch");
  }
  if (!gra::AllFinite(parameters)) {
    throw std::invalid_argument(
        "MSplineFlowMixture::SetParameters: parameters must be finite");
  }
  for (const auto &component : indices(flows_)) {
    const auto begin =
        parameters.begin() +
        static_cast<std::ptrdiff_t>(component * flow_parameter_count);
    flows_[component].SetParameters(
        std::vector<double>(begin, begin + flow_parameter_count));
  }
  std::copy(parameters.end() -
                static_cast<std::ptrdiff_t>(weight_logits_.size()),
            parameters.end(), weight_logits_.begin());
}

// Apply one AdamW step to all spline parameters and weight logits
void MSplineFlowMixture::ApplyParameterStep(const std::vector<double> &steps,
                                            double matrix_decay_factor,
                                            double bias_decay_factor) {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::ApplyParameterStep: mixture is not configured");
  }
  const std::size_t flow_parameter_count = flows_.front().Parameters().size();
  const std::size_t expected =
      flows_.size() * flow_parameter_count + weight_logits_.size();
  if (steps.size() != expected) {
    throw std::invalid_argument(
        "MSplineFlowMixture::ApplyParameterStep: parameter size mismatch");
  }
  for (const auto &component : indices(flows_)) {
    const auto begin = steps.begin() + static_cast<std::ptrdiff_t>(
                                           component * flow_parameter_count);
    flows_[component].ApplyParameterStep(
        std::vector<double>(begin, begin + flow_parameter_count),
        matrix_decay_factor, bias_decay_factor);
  }
  const std::size_t logit_offset = flows_.size() * flow_parameter_count;
  for (const auto &logit : indices(weight_logits_)) {
    if (!std::isfinite(steps[logit_offset + logit]) ||
        !std::isfinite(weight_logits_[logit] - steps[logit_offset + logit])) {
      throw std::runtime_error(
          "MSplineFlowMixture::ApplyParameterStep: non-finite weight update");
    }
    weight_logits_[logit] -= steps[logit_offset + logit];
  }
}

// Apply one optimizer step to one independently parameterized component
void MSplineFlowMixture::ApplyComponentStep(
    std::size_t component, const std::vector<double> &steps,
    double matrix_decay_factor, double bias_decay_factor) {
  flows_.at(component).ApplyParameterStep(
      steps, matrix_decay_factor, bias_decay_factor);
}

// Initialize every spline component from one adapted VEGAS grid
void MSplineFlowMixture::InitializeVegasGrid(
    const std::vector<std::vector<double>> &quantile_edges) {
  for (MSplineFlow &flow : flows_) {
    flow.InitializeVegasGrid(quantile_edges);
  }
}

// Build immutable inference data for every component
void MSplineFlowMixture::PrepareInference() {
  for (MSplineFlow &flow : flows_) {
    flow.PrepareInference();
  }
}

// Evaluate the exact spline mixture density at one point
double MSplineFlowMixture::LogDensity(const std::vector<double> &point) const {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::LogDensity: mixture is not configured");
  }
  if (flows_.size() == 1) {
    return flows_.front().LogDensity(point);
  }
  const std::vector<double> log_weights = LogWeights();
  std::vector<double> terms(flows_.size());
  double maximum = -std::numeric_limits<double>::infinity();
  for (const auto &component : indices(flows_)) {
    terms[component] =
        log_weights[component] + flows_[component].LogDensity(point);
    maximum = std::max(maximum, terms[component]);
  }
  double sum = 0.0;
  for (double term : terms) {
    sum += std::exp(term - maximum);
  }
  return maximum + std::log(sum);
}

// Evaluate exact spline mixture densities for row-major points
std::vector<double>
MSplineFlowMixture::LogDensityBatch(const std::vector<double> &points) const {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::LogDensityBatch: mixture is not configured");
  }
  if (flows_.size() == 1) {
    return flows_.front().LogDensityBatch(points);
  }
  const std::vector<double> log_weights = LogWeights();
  std::vector<std::vector<double>> component_log_densities(flows_.size());
  for (const auto &component : indices(flows_)) {
    component_log_densities[component] =
        flows_[component].LogDensityBatch(points);
  }
  const std::size_t events = component_log_densities.front().size();
  std::vector<double> output(events);
  for (std::size_t event = 0; event < events; ++event) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto &component : indices(flows_)) {
      maximum =
          std::max(maximum, log_weights[component] +
                                component_log_densities[component][event]);
    }
    double sum = 0.0;
    for (const auto &component : indices(flows_)) {
      sum += std::exp(log_weights[component] +
                      component_log_densities[component][event] - maximum);
    }
    output[event] = maximum + std::log(sum);
  }
  return output;
}

// Evaluate one mini-batch objective and its reverse-mode gradient
double MSplineFlowMixture::LossGradient(
    const std::vector<MFlowTrainingEvent> &events,
    const std::vector<double> &coefficients,
    const std::vector<std::size_t> &indices, std::size_t population_size,
    const MFlowObjective &loss_function, double uniform_mix,
    std::vector<double> &gradient) const {
  if (flows_.empty()) {
    throw std::logic_error(
        "MSplineFlowMixture::LossGradient: mixture is not configured");
  }
  if (flows_.size() == 1) {
    return flows_.front().LossGradient(events, coefficients, indices,
                                       population_size, loss_function,
                                       uniform_mix, gradient);
  }
  if (events.size() != coefficients.size() || indices.empty() ||
      population_size < indices.size() || population_size > events.size()) {
    throw std::invalid_argument(
        "MSplineFlowMixture::LossGradient: invalid mini-batch input");
  }
  double coefficient_sum = 0.0;
  for (std::size_t index : indices) {
    coefficient_sum += coefficients.at(index);
  }
  if (!(coefficient_sum > 0.0)) {
    gradient.assign(Parameters().size(), 0.0);
    return 0.0;
  }

  const double batch_scale = static_cast<double>(population_size) /
                             static_cast<double>(indices.size());
  std::vector<double> batch_points;
  std::vector<double> batch_coefficients;
  batch_points.reserve(indices.size() * Dimension());
  batch_coefficients.reserve(indices.size());
  for (std::size_t index : indices) {
    const std::vector<double> &point = events.at(index).point;
    ValidateUnitPoint(point, Dimension(), "MSplineFlowMixture::LossGradient");
    batch_points.insert(batch_points.end(), point.begin(), point.end());
    batch_coefficients.push_back(coefficients.at(index) * batch_scale);
  }

#ifdef GRANIITTI_USE_LIBTORCH
  const torch::Tensor point_tensor =
      TensorFlatPoints(batch_points, indices.size(), Dimension());
  std::vector<torch::Tensor> parameter_tensors;
  std::vector<torch::Tensor> component_log_densities;
  parameter_tensors.reserve(flows_.size());
  component_log_densities.reserve(flows_.size());
  for (const MSplineFlow &flow : flows_) {
    parameter_tensors.push_back(TensorParameters(flow.Parameters(), true));
    component_log_densities.push_back(TensorFlowLogDensity(
        flow, point_tensor, parameter_tensors.back(), nullptr));
  }
  torch::Tensor logit_tensor = TensorParameters(weight_logits_, true);
  const torch::Tensor full_logits =
      torch::cat({logit_tensor, torch::zeros({1}, NeuroJacTensorOptions())}, 0);
  const torch::Tensor log_weights = torch::log_softmax(full_logits, 0);
  const torch::Tensor flow_log_density = torch::logsumexp(
      torch::stack(component_log_densities, 1) + log_weights, 1);
  const std::vector<double> flow_log_densities = TensorVector(flow_log_density);
  std::vector<double> derivatives(indices.size(), 0.0);
  double objective_value = 0.0;
  for (const auto &event : gra::aux::indices(indices)) {
    const double mixed_log_density =
        MixedLogDensity(flow_log_densities[event], uniform_mix);
    const double flow_fraction =
        std::exp(std::log1p(-uniform_mix) + flow_log_densities[event] -
                 mixed_log_density);
    const MFlowLossEvaluation evaluation = loss_function.Evaluate(
        batch_coefficients[event], mixed_log_density, flow_fraction);
    objective_value += evaluation.value;
    derivatives[event] = evaluation.flow_log_density_derivative;
  }
  torch::Tensor derivative_tensor = torch::empty(
      {static_cast<std::int64_t>(indices.size())}, NeuroJacTensorOptions());
  std::copy(derivatives.begin(), derivatives.end(),
            derivative_tensor.data_ptr<double>());
  (flow_log_density * derivative_tensor).sum().backward();
  gradient.clear();
  for (const torch::Tensor &parameters : parameter_tensors) {
    if (!parameters.grad().defined()) {
      throw std::runtime_error("MSplineFlowMixture::LossGradient: autograd "
                               "produced no flow gradient");
    }
    const std::vector<double> component_gradient =
        TensorVector(parameters.grad());
    gradient.insert(gradient.end(), component_gradient.begin(),
                    component_gradient.end());
  }
  if (!logit_tensor.grad().defined()) {
    throw std::runtime_error("MSplineFlowMixture::LossGradient: autograd "
                             "produced no weight gradient");
  }
  const std::vector<double> logit_gradient = TensorVector(logit_tensor.grad());
  gradient.insert(gradient.end(), logit_gradient.begin(), logit_gradient.end());
#else
  const Eigen::Map<const MNeuroJacBatch> point_batch(
      batch_points.data(), static_cast<Eigen::Index>(indices.size()),
      static_cast<Eigen::Index>(Dimension()));
  std::vector<std::vector<MStandaloneCouplingCache>> caches(flows_.size());
  std::vector<MNeuroJacBatchVector> component_log_densities;
  component_log_densities.reserve(flows_.size());
  for (const auto &component : gra::aux::indices(flows_)) {
    component_log_densities.push_back(StandaloneFlowLogDensity(
        flows_[component], point_batch, flows_[component].Parameters(),
        &caches[component], nullptr));
  }
  const std::vector<double> log_weights = LogWeights();
  const std::vector<double> weights = Weights();
  std::vector<MNeuroJacBatchVector> component_adjoints(flows_.size());
  for (const auto &component : gra::aux::indices(flows_)) {
    component_adjoints[component] =
        MNeuroJacBatchVector::Zero(static_cast<Eigen::Index>(indices.size()));
  }
  std::vector<double> logit_gradient(weight_logits_.size(), 0.0);
  double objective_value = 0.0;
  for (const auto &event : gra::aux::indices(indices)) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (const auto &component : gra::aux::indices(flows_)) {
      maximum =
          std::max(maximum, log_weights[component] +
                                component_log_densities[component][event]);
    }
    double sum = 0.0;
    std::vector<double> responsibilities(flows_.size());
    for (const auto &component : gra::aux::indices(flows_)) {
      responsibilities[component] =
          std::exp(log_weights[component] +
                   component_log_densities[component][event] - maximum);
      sum += responsibilities[component];
    }
    const double flow_log_density = maximum + std::log(sum);
    const double mixed_log_density =
        MixedLogDensity(flow_log_density, uniform_mix);
    const double flow_fraction = std::exp(std::log1p(-uniform_mix) +
                                          flow_log_density - mixed_log_density);
    const MFlowLossEvaluation evaluation = loss_function.Evaluate(
        batch_coefficients[event], mixed_log_density, flow_fraction);
    objective_value += evaluation.value;
    for (const auto &component : gra::aux::indices(flows_)) {
      responsibilities[component] /= sum;
      component_adjoints[component][event] =
          evaluation.flow_log_density_derivative * responsibilities[component];
    }
    for (const auto &logit : gra::aux::indices(weight_logits_)) {
      logit_gradient[logit] += evaluation.flow_log_density_derivative *
                               (responsibilities[logit] - weights[logit]);
    }
  }
  gradient.clear();
  for (const auto &component : gra::aux::indices(flows_)) {
    std::vector<double> component_gradient(
        flows_[component].Parameters().size(), 0.0);
    StandaloneFlowBackward(flows_[component], caches[component],
                           component_adjoints[component],
                           flows_[component].Parameters(), component_gradient);
    gradient.insert(gradient.end(), component_gradient.begin(),
                    component_gradient.end());
  }
  gradient.insert(gradient.end(), logit_gradient.begin(), logit_gradient.end());
#endif
  if (!std::isfinite(objective_value)) {
    throw std::runtime_error(
        "MSplineFlowMixture::LossGradient: non-finite objective");
  }
  return objective_value;
}

// Check numerical density agreement before accepting a proposal
bool StableFlowDensity(const MSplineFlowMixture &mixture, const MNeuroJacConfig &config, std::uint32_t seed) {
  std::mt19937        random(seed);
  std::vector<double> bases(config.validation_size * mixture.Dimension());
  for (double &x : bases) { x = std::generate_canonical<double, std::numeric_limits<double>::digits>(random); }
  try {
    for (std::size_t component = 0; component < mixture.ComponentCount(); ++component) {
      const auto &flow    = mixture.Component(component);
      const auto  samples = flow.ForwardBatch(bases);
      const auto  inverse = flow.LogDensityBatch(samples.points);
      for (const auto &i : indices(inverse)) {
        if (!std::isfinite(inverse[i]) || !std::isfinite(samples.log_densities[i]) ||
            std::abs(inverse[i] - samples.log_densities[i]) > config.density_tolerance) {
          return false;
        }
      }
    }
  } catch (const std::runtime_error &) { return false; } catch (const std::invalid_argument &) {
    return false;
  }
  return true;
}

// Construct an empty adaptive spline-flow integrator
MFlowIntegrator::MFlowIntegrator() = default;

// Destroy the source-owned AdamW optimizer state
MFlowIntegrator::~MFlowIntegrator() = default;

// Read NUMERICS_NEUROJAC from one model NUMERICS file
void MFlowIntegrator::ReadParameters(const std::string &filename) {
  if (frozen_) {
    throw std::logic_error(
        "MFlowIntegrator::ReadParameters: proposal is frozen");
  }

  json document;
  try {
    document = json::parse(gra::aux::GetInputData(filename));
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MFlowIntegrator::ReadParameters: failed to parse '" + filename +
        "': " + error.what());
  }
  const json &block = document.at("NUMERICS_NEUROJAC");

  MNeuroJacConfig parsed = ParseFlowConfig(block);

  config_ = std::move(parsed);
  configured_ = false;
  buffer_.clear();
  replay_buffer_.clear();
  validation_buffer_.clear();
  best_parameters_.clear();
  adamw_first_moment_.clear();
  adamw_second_moment_.clear();
  best_validation_ess_fraction_ = 0.0;
  flow_seed_ = 0;
  adamw_step_ = 0;
  training_round_ = 0;
  validation_round_ = 0;
  validation_stale_rounds_ = 0;
  best_validation_round_ = 0;
  best_component_count_ = 0;
  vegas_initialized_ = false;
}

// Configure a fresh spline flow for one phase-space dimension
void MFlowIntegrator::Configure(std::size_t dimension,
                                std::uint32_t process_seed) {
  config_.Validate(dimension);
  const std::uint32_t seed = process_seed;
  flow_seed_ = seed;
  MNeuroJacConfig initial_config = config_;
  initial_config.flow_components = 1;
  mixture_.Configure(dimension, initial_config, seed);
  MAdamW::Configure(mixture_.Parameters().size(), adamw_first_moment_,
                    adamw_second_moment_, adamw_step_);
  optimizer_random_.seed(seed ^ 0x9e3779b9U);
  replay_random_.seed(seed ^ 0x85ebca6bU);
  buffer_.clear();
  buffer_.reserve(config_.buffer_size);
  replay_buffer_.clear();
  if (config_.replay_fraction > 0.0) {
    replay_buffer_.reserve(config_.replay_capacity);
  }
  validation_buffer_.clear();
  validation_buffer_.reserve(config_.validation_size);
  best_parameters_.clear();
  best_validation_ess_fraction_ = 0.0;
  training_round_ = 0;
  validation_round_ = 0;
  validation_stale_rounds_ = 0;
  best_validation_round_ = 0;
  best_component_count_ = 0;
  vegas_initialized_ = false;
  configured_ = true;
  frozen_ = false;
}

// Initialize the configured flow from one frozen VEGAS grid
void MFlowIntegrator::InitializeVegasGrid(
    const std::vector<std::vector<double>> &quantile_edges) {
  if (!configured_ || frozen_) {
    throw std::logic_error(
        "MFlowIntegrator::InitializeVegasGrid: adaptive flow is not available");
  }
  if (!config_.vegas_init) {
    throw std::logic_error("MFlowIntegrator::InitializeVegasGrid: VEGAS "
                           "initialization is disabled");
  }
  if (vegas_initialized_) {
    throw std::logic_error("MFlowIntegrator::InitializeVegasGrid: VEGAS grid "
                           "is already initialized");
  }
  if (training_round_ != 0 || validation_round_ != 0 || !buffer_.empty() ||
      !replay_buffer_.empty() || !validation_buffer_.empty() ||
      !best_parameters_.empty()) {
    throw std::logic_error(
        "MFlowIntegrator::InitializeVegasGrid: adaptation has already started");
  }
  if (config_.vegas_spline_init) {
    mixture_.InitializeVegasGrid(quantile_edges);
  }
  if (!StableFlowDensity(mixture_, config_, flow_seed_)) {
    throw std::invalid_argument("MFlowIntegrator::InitializeVegasGrid: inconsistent forward/inverse proposal density");
  }
  vegas_initialized_ = true;
}

// Serialize one configured frozen proposal into a backend-neutral JSON document
std::string MFlowIntegrator::SerializeModel() const {
  if (!configured_ || !frozen_) {
    throw std::logic_error("MFlowIntegrator::SerializeModel: proposal must be "
                           "configured and frozen");
  }
  json document;
  document["MODEL_FORMAT"] = NeuroJacModelFormat();
  document["MODEL_VERSION"] = NeuroJacModelVersion();
  document["DIMENSION"] = mixture_.Dimension();
  document["FLOW_SEED"] = flow_seed_;
  document["CONFIG"] = SerializeFlowConfig(config_);
  document["ACTIVE_COMPONENTS"] = mixture_.ComponentCount();
  document["PARAMETER_COUNT"] = mixture_.Parameters().size();
  document["PARAMETERS"] = mixture_.Parameters();
  return document.dump(2);
}

// Restore one backend-neutral JSON proposal for an expected phase-space
// dimension
void MFlowIntegrator::DeserializeModel(const std::string &document,
                                       std::size_t expected_dimension) {
  json model;
  try {
    model = json::parse(document);
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MFlowIntegrator::DeserializeModel: invalid JSON: " +
        std::string(error.what()));
  }

  if (model.at("MODEL_FORMAT").get<std::string>() != NeuroJacModelFormat()) {
    throw std::invalid_argument(
        "MFlowIntegrator::DeserializeModel: unknown model format");
  }
  const std::size_t version = ParseFlowCount(model.at("MODEL_VERSION"), "MODEL_VERSION");
  if (version != NeuroJacModelVersion()) {
    throw std::invalid_argument(
        "MFlowIntegrator::DeserializeModel: unsupported model version " +
        std::to_string(version));
  }
  const std::size_t dimension = ParseFlowCount(model.at("DIMENSION"), "DIMENSION");
  if (dimension != expected_dimension) {
    throw std::invalid_argument(
        "MFlowIntegrator::DeserializeModel: model phase-space dimension = " +
        std::to_string(dimension) +
        " but process dimension = " + std::to_string(expected_dimension));
  }

  MNeuroJacConfig restored_config = ParseFlowConfig(model.at("CONFIG"));
  restored_config.Validate(dimension);
  const std::size_t seed = ParseFlowCount(model.at("FLOW_SEED"), "FLOW_SEED");
  if (seed > std::numeric_limits<std::uint32_t>::max()) {
    throw std::invalid_argument("MFlowIntegrator::DeserializeModel: invalid flow seed");
  }
  const auto        restored_seed     = static_cast<std::uint32_t>(seed);
  const std::size_t active_components = ParseFlowCount(model.at("ACTIVE_COMPONENTS"), "ACTIVE_COMPONENTS");
  if (active_components == 0 ||
      active_components > restored_config.flow_components) {
    throw std::invalid_argument(
        "MFlowIntegrator::DeserializeModel: invalid active component count");
  }
  std::vector<double> restored_parameters =
      model.at("PARAMETERS").get<std::vector<double>>();
  const std::size_t parameter_count = ParseFlowCount(model.at("PARAMETER_COUNT"), "PARAMETER_COUNT");
  if (restored_parameters.size() != parameter_count) {
    throw std::invalid_argument("MFlowIntegrator::DeserializeModel: serialized "
                                "parameter count mismatch");
  }

  MSplineFlowMixture restored_mixture;
  MNeuroJacConfig active_config = restored_config;
  active_config.flow_components = active_components;
  restored_mixture.Configure(dimension, active_config, restored_seed);
  restored_mixture.SetParameters(restored_parameters);
  if (!StableFlowDensity(restored_mixture, restored_config, restored_seed)) {
    throw std::invalid_argument("MFlowIntegrator::DeserializeModel: inconsistent forward/inverse proposal density");
  }

  config_ = std::move(restored_config);
  mixture_ = std::move(restored_mixture);
  flow_seed_ = restored_seed;
  MAdamW::Configure(mixture_.Parameters().size(), adamw_first_moment_,
                    adamw_second_moment_, adamw_step_);
  optimizer_random_.seed(flow_seed_ ^ 0x9e3779b9U);
  replay_random_.seed(flow_seed_ ^ 0x85ebca6bU);
  buffer_.clear();
  buffer_.reserve(config_.buffer_size);
  replay_buffer_.clear();
  validation_buffer_.clear();
  validation_buffer_.reserve(config_.validation_size);
  best_parameters_.clear();
  best_validation_ess_fraction_ = 0.0;
  training_round_ = 0;
  validation_round_ = 0;
  validation_stale_rounds_ = 0;
  best_validation_round_ = 0;
  best_component_count_ = 0;
  vegas_initialized_ = false;
  configured_ = true;
  frozen_ = true;
}

// Save one configured frozen proposal directly to disk
void MFlowIntegrator::SaveModel(const std::string &filename) const {
  const std::string document = SerializeModel();
  gra::aux::CreateDirectory(std::filesystem::path(filename).parent_path().string());
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("MFlowIntegrator::SaveModel: cannot open '" +
                             filename + "'");
  }
  file << document << '\n';
  if (!file.good()) {
    throw std::runtime_error("MFlowIntegrator::SaveModel: failed to write '" +
                             filename + "'");
  }
}

// Load one frozen proposal directly from disk
void MFlowIntegrator::LoadModel(const std::string &filename,
                                std::size_t expected_dimension) {
  DeserializeModel(gra::aux::GetInputData(filename), expected_dimension);
}

// Draw one exact mixed flow/uniform proposal point
MFlowSample MFlowIntegrator::Sample(MRandom &random) const {
  if (!configured_) {
    throw std::logic_error("MFlowIntegrator::Sample: flow is not configured");
  }

  std::vector<double> base(mixture_.Dimension());
  for (double &coordinate : base) {
    coordinate = random.U(0.0, 1.0);
  }

  MFlowSample sample;
  if (random.U(0.0, 1.0) < config_.uniform_mix) {
    sample.point = std::move(base);
    sample.log_density =
        MixedLogDensity(mixture_.LogDensity(sample.point), config_.uniform_mix);
  } else {
    const std::size_t component = DrawFlowComponent(mixture_.Weights(), random);
    sample = mixture_.Component(component).Forward(base);
    const double flow_log_density = mixture_.ComponentCount() == 1
                                        ? sample.log_density
                                        : mixture_.LogDensity(sample.point);
    sample.log_density = MixedLogDensity(flow_log_density, config_.uniform_mix);
  }
  return sample;
}

// Draw exact mixed proposal points into contiguous row-major storage
MFlowSampleBatch MFlowIntegrator::SampleBatch(MRandom &random,
                                              std::size_t count) const {
  if (!configured_) {
    throw std::logic_error(
        "MFlowIntegrator::SampleBatch: flow is not configured");
  }
  MFlowSampleBatch samples;
  samples.dimension = mixture_.Dimension();
  if (count == 0) {
    return samples;
  }

  samples.points.resize(count * samples.dimension);
  samples.log_densities.resize(count);
  const std::size_t components = mixture_.ComponentCount();
  const std::vector<double> mixture_weights = mixture_.Weights();
  std::vector<std::vector<double>> component_bases(components);
  std::vector<std::vector<std::size_t>> component_indices(components);
  std::vector<double> uniform_points;
  std::vector<std::size_t> uniform_indices;
  uniform_points.reserve(count * samples.dimension);
  uniform_indices.reserve(count);
  for (std::size_t event = 0; event < count; ++event) {
    double *base = samples.points.data() + event * samples.dimension;
    for (std::size_t coordinate = 0; coordinate < samples.dimension;
         ++coordinate) {
      base[coordinate] = random.U(0.0, 1.0);
    }
    if (random.U(0.0, 1.0) < config_.uniform_mix) {
      uniform_indices.push_back(event);
      uniform_points.insert(uniform_points.end(), base,
                            base + samples.dimension);
    } else {
      const std::size_t component = DrawFlowComponent(mixture_weights, random);
      component_indices[component].push_back(event);
      component_bases[component].insert(component_bases[component].end(), base,
                                        base + samples.dimension);
    }
  }

  std::vector<std::size_t> source_component(count, components);
  std::vector<double> source_log_density(count, 0.0);
  for (std::size_t component = 0; component < components; ++component) {
    if (component_bases[component].empty()) {
      continue;
    }
    const MFlowSampleBatch flow_samples =
        mixture_.Component(component).ForwardBatch(component_bases[component]);
    for (const auto &event : indices(component_indices[component])) {
      const std::size_t output_event = component_indices[component][event];
      std::copy_n(flow_samples.points.data() + event * samples.dimension,
                  samples.dimension,
                  samples.points.data() + output_event * samples.dimension);
      source_component[output_event] = component;
      source_log_density[output_event] = flow_samples.log_densities[event];
    }
  }

  if (components == 1) {
    if (!uniform_points.empty()) {
      const std::vector<double> uniform_log_densities =
          mixture_.Component(0).LogDensityBatch(uniform_points);
      for (const auto &event : indices(uniform_indices)) {
        source_log_density[uniform_indices[event]] =
            uniform_log_densities[event];
      }
    }
    for (std::size_t event = 0; event < count; ++event) {
      samples.log_densities[event] =
          MixedLogDensity(source_log_density[event], config_.uniform_mix);
    }
    return samples;
  }

  std::vector<std::vector<double>> component_log_densities(
      components, std::vector<double>(count));
  const std::vector<double> log_weights = mixture_.LogWeights();
  for (std::size_t component = 0; component < components; ++component) {
    std::vector<double> evaluation_points;
    std::vector<std::size_t> evaluation_indices;
    evaluation_points.reserve(count * samples.dimension);
    evaluation_indices.reserve(count);
    for (std::size_t event = 0; event < count; ++event) {
      if (source_component[event] == component) {
        component_log_densities[component][event] = source_log_density[event];
        continue;
      }
      evaluation_indices.push_back(event);
      const double *point = samples.points.data() + event * samples.dimension;
      evaluation_points.insert(evaluation_points.end(), point,
                               point + samples.dimension);
    }
    if (!evaluation_points.empty()) {
      const std::vector<double> evaluated =
          mixture_.Component(component).LogDensityBatch(evaluation_points);
      for (const auto &event : indices(evaluation_indices)) {
        component_log_densities[component][evaluation_indices[event]] =
            evaluated[event];
      }
    }
  }
  for (std::size_t event = 0; event < count; ++event) {
    double maximum = -std::numeric_limits<double>::infinity();
    for (std::size_t component = 0; component < components; ++component) {
      maximum =
          std::max(maximum, log_weights[component] +
                                component_log_densities[component][event]);
    }
    double density_sum = 0.0;
    for (std::size_t component = 0; component < components; ++component) {
      density_sum +=
          std::exp(log_weights[component] +
                   component_log_densities[component][event] - maximum);
    }
    samples.log_densities[event] =
        MixedLogDensity(maximum + std::log(density_sum), config_.uniform_mix);
  }
  return samples;
}

// Evaluate the exact mixed proposal density
double MFlowIntegrator::LogDensity(const std::vector<double> &point) const {
  if (!configured_) {
    throw std::logic_error(
        "MFlowIntegrator::LogDensity: flow is not configured");
  }
  return MixedLogDensity(mixture_.LogDensity(point), config_.uniform_mix);
}

// Compute the fresh target evaluations needed for the next training round
std::size_t MFlowIntegrator::FreshEventCount() const {
  if (!configured_ || frozen_) {
    throw std::logic_error(
        "MFlowIntegrator::FreshEventCount: adaptive flow is not available");
  }
  const std::size_t target_replay = static_cast<std::size_t>(
      config_.replay_fraction * static_cast<double>(config_.buffer_size));
  const std::size_t replay_events =
      std::min(target_replay, replay_buffer_.size());
  return config_.buffer_size - replay_events;
}

// Compute the held-out target evaluations needed for one validation round
std::size_t MFlowIntegrator::ValidationEventCount() const {
  if (!configured_ || frozen_) {
    throw std::logic_error("MFlowIntegrator::ValidationEventCount: adaptive "
                           "flow is not available");
  }
  return config_.validation_size;
}

// Clear all events from the current adaptive round
void MFlowIntegrator::ClearBuffer() {
  if (frozen_) {
    throw std::logic_error("MFlowIntegrator::ClearBuffer: proposal is frozen");
  }
  buffer_.clear();
}

// Clear all held-out events from the current validation round
void MFlowIntegrator::ClearValidationBuffer() {
  if (frozen_) {
    throw std::logic_error(
        "MFlowIntegrator::ClearValidationBuffer: proposal is frozen");
  }
  validation_buffer_.clear();
}

// Append one collection buffer after worker synchronization
void MFlowIntegrator::BufferEvents(std::vector<MFlowTrainingEvent> events) {
  if (!configured_ || frozen_) {
    throw std::logic_error(
        "MFlowIntegrator::BufferEvents: proposal cannot accept training data");
  }
  if (buffer_.size() > config_.buffer_size ||
      events.size() > config_.buffer_size - buffer_.size()) {
    throw std::length_error(
        "MFlowIntegrator::BufferEvents: active buffer capacity exceeded");
  }
  ValidateFlowEvents(events, mixture_.Dimension(),
                     "MFlowIntegrator::BufferEvents");
  buffer_.insert(buffer_.end(), std::make_move_iterator(events.begin()),
                 std::make_move_iterator(events.end()));
}

// Append held-out events that are never used for optimizer updates or replay
void MFlowIntegrator::BufferValidationEvents(
    std::vector<MFlowTrainingEvent> events) {
  if (!configured_ || frozen_) {
    throw std::logic_error("MFlowIntegrator::BufferValidationEvents: proposal "
                           "cannot accept validation data");
  }
  if (validation_buffer_.size() > config_.validation_size ||
      events.size() > config_.validation_size - validation_buffer_.size()) {
    throw std::length_error("MFlowIntegrator::BufferValidationEvents: "
                            "configured validation capacity exceeded");
  }
  ValidateFlowEvents(events, mixture_.Dimension(),
                     "MFlowIntegrator::BufferValidationEvents");
  validation_buffer_.insert(validation_buffer_.end(),
                            std::make_move_iterator(events.begin()),
                            std::make_move_iterator(events.end()));
}

// Append a random replay subset to the fresh current-round events
// [REFERENCE: Heimel et al., arxiv.org/abs/2212.06172]
std::size_t MFlowIntegrator::AppendReplayEvents() {
  const std::size_t target_replay = static_cast<std::size_t>(
      config_.replay_fraction * static_cast<double>(config_.buffer_size));
  const std::size_t available_slots = buffer_.size() < config_.buffer_size
                                          ? config_.buffer_size - buffer_.size()
                                          : 0;
  const std::size_t replay_events =
      std::min({target_replay, replay_buffer_.size(), available_slots});
  if (replay_events == 0) {
    return 0;
  }

  std::vector<std::size_t> order(replay_buffer_.size());
  std::iota(order.begin(), order.end(), 0);
  std::shuffle(order.begin(), order.end(), replay_random_);
  for (std::size_t i = 0; i < replay_events; ++i) {
    buffer_.push_back(replay_buffer_[order[i]]);
  }
  return replay_events;
}

// Retain fresh current-round events in the bounded replay pool
void MFlowIntegrator::RetainFreshEvents(std::size_t fresh_events) {
  if (std::fpclassify(config_.replay_fraction) == FP_ZERO ||
      config_.replay_capacity == 0 ||
      fresh_events == 0) {
    return;
  }
  if (fresh_events >= config_.replay_capacity) {
    replay_buffer_.assign(
        buffer_.begin() +
            static_cast<std::ptrdiff_t>(fresh_events - config_.replay_capacity),
        buffer_.begin() + static_cast<std::ptrdiff_t>(fresh_events));
    return;
  }

  const std::size_t retained_capacity = config_.replay_capacity - fresh_events;
  if (replay_buffer_.size() > retained_capacity) {
    const std::size_t expired = replay_buffer_.size() - retained_capacity;
    replay_buffer_.erase(replay_buffer_.begin(),
                         replay_buffer_.begin() +
                             static_cast<std::ptrdiff_t>(expired));
  }
  replay_buffer_.insert(replay_buffer_.end(), buffer_.begin(),
                        buffer_.begin() +
                            static_cast<std::ptrdiff_t>(fresh_events));
}

// Optimize one epoch with coefficients evaluated at the current parameters
double MFlowIntegrator::TrainEpoch(std::vector<std::size_t> &order, const std::vector<double> &boost_coefficients,
                                   std::size_t epoch, std::size_t threads, MFlowTrainingReport &report) {
  const bool          boosting            = mixture_.ComponentCount() > 1;
  const std::size_t   active_batches      = 1 + (order.size() - 1) / config_.batch_size;
  const std::size_t   scheduled_batches   = 1 + (config_.buffer_size - 1) / config_.batch_size;
  const std::size_t   total_epochs        = config_.rounds * config_.epochs;
  const std::size_t   global_epoch        = training_round_ * config_.epochs + epoch;
  const double        alpha               = config_.loss.kind == FlowLoss::Renyi
                                                ? ScheduledRenyiAlpha(config_.loss, global_epoch, total_epochs)
                                                : (config_.loss.kind == FlowLoss::Chi2 ? 2.0 : 1.0);
  const auto          objective           = CreateFlowObjective(config_.loss.kind, alpha);
  const auto          component_objective = boosting ? CreateFlowObjective(FlowLoss::ForwardKL, 1.0) : nullptr;
  const std::vector<double> points = FlattenFlowEventPoints(buffer_, mixture_.Dimension());
  std::vector<double> coefficients;
  report.loss_alpha = boosting ? 2.0 : alpha;
  std::shuffle(order.begin(), order.end(), optimizer_random_);
  double      loss  = 0.0;
  std::size_t batch = 0;
  for (std::size_t begin = 0; begin < order.size(); begin += config_.batch_size) {
    if (coefficients.empty() || boosting || objective->UsesProposalDensity()) {
      if (boosting || objective->UsesProposalDensity()) { mixture_.PrepareInference(); }
      coefficients = boosting ? BoostChi2Coefficients(mixture_, buffer_, points, boost_coefficients, config_.uniform_mix, threads, config_.batch_size)
                              : NormalizedObjectiveCoefficients(mixture_, buffer_, points, *objective, config_.uniform_mix, threads, config_.batch_size);
    }
    const std::size_t              end = std::min(begin + config_.batch_size, order.size());
    const std::vector<std::size_t> batch_indices(order.begin() + begin, order.begin() + end);
    std::vector<double>            gradient;
    const std::size_t              component = mixture_.ComponentCount() - 1;
    loss += boosting ? mixture_.Component(component).LossGradient(buffer_, coefficients, batch_indices, order.size(),
                                                                  *component_objective, 0.0, gradient)
                     : mixture_.LossGradient(buffer_, coefficients, batch_indices, order.size(), *objective,
                                             config_.uniform_mix, gradient);
    const std::uint64_t schedule_step =
        global_epoch * scheduled_batches + (++batch * scheduled_batches / active_batches) - 1;
    report.learning_rate = ScheduledLearningRate(config_, schedule_step);
    MAdamW::Steps(gradient, report.learning_rate, config_.gradient_clip, config_, adamw_first_moment_,
                  adamw_second_moment_, adamw_step_);
    mixture_.ApplyComponentStep(component, gradient, report.learning_rate * config_.adamw_weight_decay,
                                report.learning_rate * config_.adamw_decay_biases);
  }
  return loss / active_batches;
}

// Optimize one adaptive round using the selected objective
MFlowTrainingReport MFlowIntegrator::TrainRound(std::size_t threads) {
  if (!configured_ || frozen_ || buffer_.empty()) {
    throw std::logic_error(
        "MFlowIntegrator::TrainRound: no adaptive buffer is available");
  }
  if (training_round_ == 0 && config_.vegas_init && !vegas_initialized_) {
    throw std::logic_error(
        "MFlowIntegrator::TrainRound: VEGAS grid is not initialized");
  }

  MFlowTrainingReport report;
  report.fresh_events = buffer_.size();
  report.replayed_events = AppendReplayEvents();
  report.buffered_events = buffer_.size();
  std::vector<double> log_weights(buffer_.size(),
                                  -std::numeric_limits<double>::infinity());
  double maximum_log_weight = -std::numeric_limits<double>::infinity();
  std::vector<std::size_t> order;
  order.reserve(buffer_.size());

  for (const auto &i : indices(buffer_)) {
    const double absolute_target = std::abs(buffer_[i].target);
    if (std::fpclassify(absolute_target) == FP_ZERO) {
      continue;
    }
    ++report.nonzero_events;
    log_weights[i] =
        std::log(absolute_target) - buffer_[i].collection_log_density;
    maximum_log_weight = std::max(maximum_log_weight, log_weights[i]);
    order.push_back(i);
  }
  if (report.nonzero_events == 0) {
    throw std::runtime_error(
        "MFlowIntegrator::TrainRound: all buffered targets are zero");
  }
  report.vegas_initialized = training_round_ == 0 && vegas_initialized_;

  double scaled_weight_sum = 0.0;
  double scaled_weight2_sum = 0.0;
  long double scaled_signed_weight_sum = 0.0L;
  for (const auto &i : indices(buffer_)) {
    if (!std::isfinite(log_weights[i])) {
      continue;
    }
    const double scaled_weight = std::exp(log_weights[i] - maximum_log_weight);
    scaled_weight_sum += scaled_weight;
    scaled_weight2_sum += scaled_weight * scaled_weight;
    scaled_signed_weight_sum +=
        buffer_[i].target > 0.0 ? scaled_weight : -scaled_weight;
  }

  const long double sample_count = static_cast<long double>(buffer_.size());
  const long double scaled_mean = scaled_signed_weight_sum / sample_count;
  if (!(std::fpclassify(scaled_mean) == FP_ZERO)) {
    const long double log_absolute_integral =
        static_cast<long double>(maximum_log_weight) +
        std::log(std::abs(scaled_mean));
    if (log_absolute_integral <= std::log(static_cast<long double>(
                                     std::numeric_limits<double>::max()))) {
      report.integral =
          std::copysign(static_cast<double>(std::exp(log_absolute_integral)),
                        static_cast<double>(scaled_mean));
    } else {
      report.integral = std::copysign(std::numeric_limits<double>::infinity(),
                                      static_cast<double>(scaled_mean));
    }
    if (buffer_.size() > 1) {
      const long double relative_second_moment =
          sample_count * static_cast<long double>(scaled_weight2_sum) /
          (scaled_signed_weight_sum * scaled_signed_weight_sum);
      const long double relative_variance = std::max(
          (relative_second_moment - 1.0L) / (sample_count - 1.0L), 0.0L);
      report.relative_error = static_cast<double>(std::sqrt(relative_variance));
    }
  }
  report.collection_ess_fraction =
      scaled_weight_sum * scaled_weight_sum /
      (scaled_weight2_sum * static_cast<double>(buffer_.size()));
  report.maximum_abs_weight =
      maximum_log_weight < std::log(std::numeric_limits<double>::max())
          ? std::exp(maximum_log_weight)
          : std::numeric_limits<double>::infinity();

  const std::size_t stage_rounds =
      std::max<std::size_t>(1, config_.rounds / config_.flow_components);
  if (training_round_ > 0 && training_round_ % stage_rounds == 0 &&
      mixture_.ComponentCount() < config_.flow_components) {
    mixture_.AddComponent(config_, flow_seed_, 0.05);
    MAdamW::Configure(mixture_.Component(mixture_.ComponentCount() - 1)
                          .Parameters()
                          .size(),
                      adamw_first_moment_, adamw_second_moment_, adamw_step_);
    validation_stale_rounds_ = 0;
    report.component_added = true;
  }
  const bool boosting = mixture_.ComponentCount() > 1;
  std::vector<double> boost_coefficients;
  if (boosting) {
    const std::unique_ptr<MFlowObjective> chi2 =
        CreateFlowObjective(FlowLoss::Chi2, 2.0);
    boost_coefficients = NormalizedObjectiveCoefficients(
        mixture_, buffer_, FlattenFlowEventPoints(buffer_, mixture_.Dimension()), *chi2,
        config_.uniform_mix, threads, config_.batch_size);
  }

  double accumulated_loss = 0.0;
  std::size_t accepted_epochs  = 0;
  for (std::size_t epoch = 0; epoch < config_.epochs; ++epoch) {
    const MSplineFlowMixture previous = mixture_;
    bool                     stable   = false;
    double                   loss     = 0.0;
    try {
      loss = TrainEpoch(order, boost_coefficients, epoch, threads, report);
      mixture_.PrepareInference();
      stable = StableFlowDensity(mixture_, config_, flow_seed_);
    } catch (const std::runtime_error &) {
      // Reject numerical optimizer failures before collecting further samples
    }
    if (stable) {
      accumulated_loss += loss;
      ++accepted_epochs;
    } else {
      mixture_ = previous;
      MAdamW::Configure(mixture_.Component(mixture_.ComponentCount() - 1).Parameters().size(), adamw_first_moment_,
                        adamw_second_moment_, adamw_step_);
      ++report.rejected_epochs;
    }
  }
  report.loss = accepted_epochs > 0 ? accumulated_loss / accepted_epochs : 0.0;
  const bool stage_complete =
      (training_round_ + 1) % stage_rounds == 0 ||
      training_round_ + 1 == config_.rounds;
  if (boosting) {
    const auto [old_log, new_log] = BoostLogDensities(mixture_, buffer_);
    report.boost_weight = OptimizeBoostWeight(
        buffer_, old_log, new_log, config_.uniform_mix);
    mixture_.SetLastWeight(
        stage_complete ? report.boost_weight
                       : std::max(0.05, report.boost_weight));
  }
  ++training_round_;
  RetainFreshEvents(report.fresh_events);
  report.replay_pool_events = replay_buffer_.size();
  buffer_.clear();
  mixture_.PrepareInference();
  report.component_weights = mixture_.Weights();
  return report;
}

// Score the VEGAS-mapped proposal as the round-zero checkpoint
MFlowValidationReport MFlowIntegrator::ValidateBaseline() {
  if (!vegas_initialized_ || training_round_ != 0 || validation_round_ != 0 ||
      !best_parameters_.empty()) {
    throw std::logic_error(
        "MFlowIntegrator::ValidateBaseline: VEGAS baseline is not available");
  }
  return ValidateBufferedProposal(false);
}

// Score the updated proposal and retain the best held-out checkpoint
MFlowValidationReport MFlowIntegrator::ValidateRound() {
  return ValidateBufferedProposal(true);
}

// Score one proposal on a common held-out sample without clipping its ranking
std::pair<MFlowValidationReport, double> ScoreFlow(const MSplineFlowMixture              &mixture,
                                                   const std::vector<MFlowTrainingEvent> &events, double uniform_mix) {
  std::vector<double> points;
  points.reserve(events.size() * mixture.Dimension());
  for (const MFlowTrainingEvent &event : events) {
    points.insert(points.end(), event.point.begin(), event.point.end());
  }
  const std::vector<double> flow_log_densities = mixture.LogDensityBatch(points);

  MFlowValidationReport report;
  double                log_variance_ratio = std::numeric_limits<double>::infinity();
  report.events                            = events.size();
  double maximum_log_absolute_integral =
      -std::numeric_limits<double>::infinity();
  double maximum_log_second_moment = -std::numeric_limits<double>::infinity();
  std::vector<double> log_absolute_integrands(events.size(), -std::numeric_limits<double>::infinity());
  std::vector<double> log_second_integrands(events.size(), -std::numeric_limits<double>::infinity());
  for (const auto &event : indices(events)) {
    const double absolute_target = std::abs(events[event].target);
    if (std::fpclassify(absolute_target) == FP_ZERO) {
      continue;
    }
    ++report.nonzero_events;
    const double log_absolute_target = std::log(absolute_target);
    const double current_log_density = MixedLogDensity(flow_log_densities[event], uniform_mix);
    log_absolute_integrands[event]   = log_absolute_target - events[event].collection_log_density;
    log_second_integrands[event] =
        2.0 * log_absolute_target - events[event].collection_log_density - current_log_density;
    maximum_log_absolute_integral =
        std::max(maximum_log_absolute_integral, log_absolute_integrands[event]);
    maximum_log_second_moment =
        std::max(maximum_log_second_moment, log_second_integrands[event]);
  }
  if (report.nonzero_events == 0) {
    report.relative_variance = std::numeric_limits<double>::infinity();
  } else {
    double absolute_integral_sum = 0.0;
    double second_moment_sum = 0.0;
    for (const auto &event : indices(events)) {
      if (!std::isfinite(log_absolute_integrands[event])) {
        continue;
      }
      absolute_integral_sum += std::exp(log_absolute_integrands[event] -
                                        maximum_log_absolute_integral);
      second_moment_sum +=
          std::exp(log_second_integrands[event] - maximum_log_second_moment);
    }
    const double log_event_count       = std::log(static_cast<double>(events.size()));
    const double log_absolute_integral = maximum_log_absolute_integral +
                                         std::log(absolute_integral_sum) -
                                         log_event_count;
    const double log_second_moment = maximum_log_second_moment +
                                     std::log(second_moment_sum) -
                                     log_event_count;
    log_variance_ratio = log_second_moment - 2.0 * log_absolute_integral;
    report.ess_fraction =
        std::min(1.0, std::exp(std::min(-log_variance_ratio, 0.0)));
    report.relative_variance =
        log_variance_ratio < std::log(std::numeric_limits<double>::max())
            ? std::max(std::exp(log_variance_ratio) - 1.0, 0.0)
            : std::numeric_limits<double>::infinity();
  }

  return {report, log_variance_ratio};
}

// Score one buffered proposal with optional neural-round advancement
MFlowValidationReport MFlowIntegrator::ValidateBufferedProposal(bool advance_round) {
  if (!configured_ || frozen_ || validation_buffer_.empty()) {
    throw std::logic_error(
        "MFlowIntegrator::ValidateRound: no held-out "
        "validation buffer is available");
  }

  auto [report, score] = ScoreFlow(mixture_, validation_buffer_, config_.uniform_mix);
  double best_score    = std::numeric_limits<double>::infinity();
  if (!best_parameters_.empty() && report.nonzero_events > 0) {
    MSplineFlowMixture best;
    MNeuroJacConfig    best_config = config_;
    best_config.flow_components    = best_component_count_;
    best.Configure(mixture_.Dimension(), best_config, flow_seed_);
    best.SetParameters(best_parameters_);
    const auto previous           = ScoreFlow(best, validation_buffer_, config_.uniform_mix);
    best_validation_ess_fraction_ = previous.first.ess_fraction;
    best_score                    = previous.second;
  }

  if (advance_round) {
    ++validation_round_;
  }
  if (report.nonzero_events > 0) {
    // Compare the absolute gain in unclipped ESS on these same events
    report.improved = best_parameters_.empty() ||
                      (score < best_score &&
                       -score + std::log(-std::expm1(score - best_score)) > std::log(config_.validation_min_delta));
    if (report.improved) {
      best_parameters_ = mixture_.Parameters();
      best_component_count_ = mixture_.ComponentCount();
      best_validation_ess_fraction_ = report.ess_fraction;
      best_validation_round_ = validation_round_;
      validation_stale_rounds_ = 0;
    } else {
      ++validation_stale_rounds_;
    }
  }
  report.best_ess_fraction = best_validation_ess_fraction_;
  report.best_round = best_validation_round_;
  report.stale_rounds = validation_stale_rounds_;
  report.early_stop = mixture_.ComponentCount() == config_.flow_components &&
                      config_.validation_patience > 0 &&
                      validation_stale_rounds_ >= config_.validation_patience;
  validation_buffer_.clear();
  return report;
}

// Freeze trained parameters for integration and event generation
void MFlowIntegrator::Freeze() {
  if (!configured_) {
    throw std::logic_error("MFlowIntegrator::Freeze: flow is not configured");
  }
  if (!buffer_.empty()) {
    throw std::logic_error(
        "MFlowIntegrator::Freeze: buffered events have not been trained");
  }
  if (!validation_buffer_.empty()) {
    throw std::logic_error("MFlowIntegrator::Freeze: buffered validation "
                           "events have not been scored");
  }
  if (!best_parameters_.empty()) {
    if (mixture_.ComponentCount() != best_component_count_) {
      MNeuroJacConfig best_config = config_;
      best_config.flow_components = best_component_count_;
      mixture_.Configure(mixture_.Dimension(), best_config, flow_seed_);
    }
    mixture_.SetParameters(best_parameters_);
  }
  mixture_.PrepareInference();
  if (!StableFlowDensity(mixture_, config_, flow_seed_)) {
    throw std::runtime_error("MFlowIntegrator::Freeze: inconsistent forward/inverse proposal density");
  }
  replay_buffer_.clear();
  vegas_initialized_ = false;
  frozen_ = true;
}

// Compute the stable steering label for one activation function
std::string FlowActivationName(FlowActivation activation) {
  switch (activation) {
  case FlowActivation::Relu:
    return "relu";
  case FlowActivation::Silu:
    return "silu";
  case FlowActivation::Tanh:
    return "tanh";
  }
  throw std::invalid_argument(
      "FlowActivationName: unknown activation function");
}

// Compute the stable steering label for one learning-rate schedule
std::string FlowLRScheduleName(FlowLRSchedule schedule) {
  switch (schedule) {
  case FlowLRSchedule::None:
    return "none";
  case FlowLRSchedule::Cosine:
    return "cosine";
  case FlowLRSchedule::Exponential:
    return "exponential";
  }
  throw std::invalid_argument(
      "FlowLRScheduleName: unknown learning-rate schedule");
}

// Compute the stable steering label for one coordinate permutation mode
std::string FlowPermutationName(FlowPermutation permutation) {
  switch (permutation) {
  case FlowPermutation::None:
    return "none";
  case FlowPermutation::Interleave:
    return "interleave";
  case FlowPermutation::Random:
    return "random";
  }
  throw std::invalid_argument("FlowPermutationName: unknown permutation mode");
}

} // namespace gra::neurojac
