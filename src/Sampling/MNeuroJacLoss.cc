// NEUROJAC adaptive normalizing flow objectives
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <algorithm>
#include <cmath>
#include <memory>
#include <stdexcept>
#include <string>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Sampling/MNeuroJacLoss.h"

namespace gra {

namespace neurojac {

namespace {

// Shared normalized weighted log density objective
class MWeightedLogObjective : public MFlowObjective {
public:
  // Evaluate one weighted negative log density contribution
  // L_i = -c_i log q_mix(x_i), dL_i/dlog q_flow = -c_i q_flow/q_mix
  MFlowLossEvaluation Evaluate(double coefficient, double mixed_log_density,
                               double flow_fraction) const final {
    return {-coefficient * mixed_log_density, -coefficient * flow_fraction};
  }
};

// Mass covering forward Kullback Leibler objective
class MForwardKLObjective final : public MWeightedLogObjective {
public:
  // Forward KL coefficients are independent of the proposal
  bool UsesProposalDensity() const final { return false; }

  // Compute the off proposal forward KL event coefficient
  // log c_i = log|f_i| - log q_collection(x_i)
  double LogCoefficient(double log_absolute_target,
                        double collection_log_density,
                        double mixed_log_density) const final {
    static_cast<void>(mixed_log_density);
    return log_absolute_target - collection_log_density;
  }
};

// Direct importance sampling second moment objective
class MChi2Objective final : public MFlowObjective {
public:
  // Chi squared coefficients are independent of the proposal
  bool UsesProposalDensity() const final { return false; }

  // Compute the off proposal chi squared event coefficient
  // log c_i = 2log|f_i| - log q_collection(x_i)
  double LogCoefficient(double log_absolute_target,
                        double collection_log_density,
                        double mixed_log_density) const final {
    static_cast<void>(mixed_log_density);
    return 2.0 * log_absolute_target - collection_log_density;
  }

  // Evaluate one weighted inverse density contribution
  // L_i = c_i/q_mix(x_i)
  MFlowLossEvaluation Evaluate(double coefficient, double mixed_log_density,
                               double flow_fraction) const final {
    const double value = coefficient * std::exp(-mixed_log_density);
    return {value, -value * flow_fraction};
  }
};

// Annealed Renyi objective represented by global escort weights
class MRenyiObjective final : public MWeightedLogObjective {
public:
  // Construct one fixed order Renyi objective
  explicit MRenyiObjective(double alpha) : alpha_(alpha) {
    if (!std::isfinite(alpha_) || alpha_ < 1.0) {
      throw std::invalid_argument(
          "MRenyiObjective: alpha must be finite and at least one");
    }
  }

  // Renyi escort weights depend on the current proposal
  bool UsesProposalDensity() const final { return true; }

  // Compute one logarithmic Renyi escort weight
  // log c_i = alpha log|f_i| - log q_collection + (1-alpha)log q_mix
  double LogCoefficient(double log_absolute_target,
                        double collection_log_density,
                        double mixed_log_density) const final {
    return alpha_ * log_absolute_target - collection_log_density +
           (1.0 - alpha_) * mixed_log_density;
  }

private:
  double alpha_ = 1.0;
};

} // namespace

// Validate one adaptive objective configuration
void MFlowLossConfig::Validate() const {
  if (!std::isfinite(alpha_initial) || !std::isfinite(alpha_final) ||
      !std::isfinite(annealing_fraction) || alpha_initial < 1.0 ||
      alpha_final < alpha_initial || !(annealing_fraction > 0.0) ||
      annealing_fraction > 1.0) {
    throw std::invalid_argument(
        "NUMERICS_NEUROJAC: invalid loss alpha interval or anneal_frac");
  }
  switch (kind) {
  case FlowLoss::ForwardKL:
  case FlowLoss::Chi2:
  case FlowLoss::Renyi:
    break;
  default:
    throw std::invalid_argument("NUMERICS_NEUROJAC: invalid loss mode");
  }
}

// Construct one immutable adaptive objective
std::unique_ptr<MFlowObjective> CreateFlowObjective(FlowLoss loss,
                                                    double renyi_alpha) {
  switch (loss) {
  case FlowLoss::ForwardKL:
    return std::make_unique<MForwardKLObjective>();
  case FlowLoss::Chi2:
    return std::make_unique<MChi2Objective>();
  case FlowLoss::Renyi:
    return std::make_unique<MRenyiObjective>(renyi_alpha);
  }
  throw std::invalid_argument("CreateFlowObjective: unknown loss mode");
}

// Compute the scheduled Renyi order for one training epoch
// alpha(t) = alpha_0 + (alpha_1-alpha_0)[1-cos(pi min(t/f,1))]/2
double ScheduledRenyiAlpha(const MFlowLossConfig &config, std::size_t epoch,
                           std::size_t total_epochs) {
  config.Validate();
  if (total_epochs == 0 || epoch >= total_epochs) {
    throw std::invalid_argument(
        "ScheduledRenyiAlpha: epoch is outside the training schedule");
  }
  if (total_epochs == 1) {
    return config.alpha_final;
  }
  const double progress =
      static_cast<double>(epoch) / static_cast<double>(total_epochs - 1);
  const double annealed_progress =
      std::min(progress / config.annealing_fraction, 1.0);
  const double blend = 0.5 * (1.0 - std::cos(math::PI * annealed_progress));
  return config.alpha_initial +
         (config.alpha_final - config.alpha_initial) * blend;
}

// Compute the stable steering label for one loss mode
std::string FlowLossName(FlowLoss loss) {
  switch (loss) {
  case FlowLoss::ForwardKL:
    return "forward_kl";
  case FlowLoss::Chi2:
    return "chi2";
  case FlowLoss::Renyi:
    return "renyi";
  }
  throw std::invalid_argument("FlowLossName: unknown loss mode");
}

} // namespace neurojac

} // namespace gra
