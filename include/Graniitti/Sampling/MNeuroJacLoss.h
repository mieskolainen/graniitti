// NEUROJAC adaptive normalizing flow objectives
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNEUROJACLOSS_H
#define MNEUROJACLOSS_H

#include <cstddef>
#include <memory>
#include <string>

namespace gra {
namespace neurojac {

// Supported adaptive normalizing flow objectives
enum class FlowLoss { ForwardKL, Chi2, Renyi };

// Steering for one adaptive normalizing flow objective
struct MFlowLossConfig {
  FlowLoss kind = FlowLoss::Renyi;
  double alpha_initial = 1.0;
  double alpha_final = 2.0;
  double annealing_fraction = 0.7;

  // Validate the objective parameters
  void Validate() const;
};

// One scalar objective value and its flow log density derivative
struct MFlowLossEvaluation {
  double value = 0.0;
  double flow_log_density_derivative = 0.0;
};

// Immutable adaptive objective interface shared by all flow backends
class MFlowObjective {
public:
  // Destroy one objective through its interface
  virtual ~MFlowObjective() = default;

  // Determine whether coefficient construction needs the current proposal
  virtual bool UsesProposalDensity() const = 0;

  // Compute one unnormalized logarithmic event coefficient
  virtual double LogCoefficient(double log_absolute_target,
                                double collection_log_density,
                                double mixed_log_density) const = 0;

  // Evaluate one normalized event contribution and reverse derivative
  virtual MFlowLossEvaluation Evaluate(double coefficient,
                                      double mixed_log_density,
                                      double flow_fraction) const = 0;
};

// Compute the stable steering label for one loss mode
std::string FlowLossName(FlowLoss loss);

// Construct one immutable adaptive objective
std::unique_ptr<MFlowObjective> CreateFlowObjective(FlowLoss loss,
                                                  double renyi_alpha);

// Compute the scheduled Renyi order for one training epoch
double ScheduledRenyiAlpha(const MFlowLossConfig &config, std::size_t epoch,
                          std::size_t total_epochs);

} // namespace neurojac
} // namespace gra

#endif
