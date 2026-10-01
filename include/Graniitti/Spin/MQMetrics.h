// Integrated spin densities and quantum correlation metrics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MQMETRICS_H
#define MQMETRICS_H

#include <array>
#include <cstdint>
#include <iosfwd>
#include <numeric>
#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Spin/MSpinDensity.h"
#include "json.hpp"

namespace gra::spin {

enum class SpinType { Production, FinalPair, IntermediatePair, ProtonCentral, ProtonPair };
enum class SpinStatus { Ready, Empty, Incomplete };

// Define the ordered physical helicities and the ensemble whose density is integrated
struct SpinSpec {
  std::string        name;
  SpinType           type = SpinType::Production;
  std::vector<int>   pdg  = {0, 0};
  // Helicity labels are twice the physical helicity
  std::vector<std::vector<int>> helicity = {{}, {}};
  std::string                     frame;
  std::string                     measure;

  // Compute the dimension of the retained spin space
  std::size_t Dimension() const {
    return std::accumulate(helicity.begin(), helicity.end(), std::size_t{1},
                           [](std::size_t n, const auto &h) { return n * h.size(); });
  }

  // Compare serialized physical spin definitions
  bool operator==(const SpinSpec &) const = default;
};

// Store the density, its integral and metrics in the order purity, S, S1, S2, N, log2(1+2N)
struct SpinResult {
  SpinSpec                      spec;
  SpinStatus                    status = SpinStatus::Empty;
  MMatrix<std::complex<double>> rho;
  double                        integral = 0.0;
  std::array<double, 6>         value    = {};
  std::array<double, 6>         error    = {};
  bool                          errors   = false;
  std::uint64_t                 samples  = 0;
  std::uint64_t                 failures = 0;
};

// Own independent worker sums and event-local densities without owning amplitude physics
class MQMetrics {
 public:
  using Density = MMatrix<std::complex<double>>;

  // Validate all spin spaces and allocate bounded integration storage
  void Configure(std::vector<SpinSpec> specs, std::size_t groups, std::size_t max_dimension, std::string note = "",
                  std::array<int, 2> protons = {});
  // Compute whether this process configured a spin metrics report
  bool Configured() const { return groups_ != 0; }
  // Compute whether the current integration trial requests amplitude densities
  bool Active() const { return active_ && !specs_.empty(); }
  // Compute the immutable spin definitions
  const std::vector<SpinSpec> &Specs() const { return specs_; }
  // Compute the physical final-pair density index when available
  std::optional<std::size_t> FinalPair() const { return pair_; }
  // Find a named spin space prepared during amplitude initialization
  std::optional<std::size_t> Find(const std::string &name) const;
  // Start one trial without retaining densities from any previous event
  void Begin(bool active);
  // Clear Born densities before projecting a completed screened amplitude
  void ClearEvent();
  // Add an orthogonal block ordered as (proton transition, spectator, central spin) x unobserved
  void Add(std::size_t index, const Density &amplitude, double normalization, std::size_t proton_rows = 0);
  // Commit one trial with its amplitude-free measure and globally unique sample index
  void Observe(double measure, double log_inverse_density, std::uint64_t sample);
  // Clear integration sums while retaining the configured physical spaces
  void Reset();
  // Combine independent worker sums with identical physical definitions
  void Merge(const MQMetrics &other);
  // Compute normalized densities and equal-group jackknife metric uncertainties
  std::vector<SpinResult> Results() const;
  // Save physical definitions, matrix sums and jackknife groups
  nlohmann::json Serialize() const;
  // Restore a compatible integration state after validating the input
  void Deserialize(const nlohmann::json &input);
  // Print integrated metrics and the definition of each ensemble
  void Print(std::ostream &output) const;

 private:
  // Keep the highest-index trial separate to balance groups without subtracting large weights
  struct Group {
    std::uint64_t        count       = 0;
    std::uint64_t        last_sample = 0;
    std::vector<Density> rest;
    std::vector<Density> last;
  };
  std::vector<SpinSpec>      specs_;
  std::vector<std::optional<std::size_t>> proton_;
  std::optional<std::size_t> pair_;
  std::size_t                groups_        = 0;
  std::size_t                max_dimension_ = 0;
  std::string                note_;
  std::vector<Group>         sums_;
  std::vector<Density>       event_;
  std::vector<bool>          ready_;
  std::vector<bool>          event_failed_;
  std::vector<std::uint64_t> failures_;
  bool                       active_ = false;
};

}  // namespace gra::spin

#endif
