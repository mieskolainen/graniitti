// Standalone VEGAS importance sampler and integration routines
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MVEGAS_H
#define MVEGAS_H

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"
using gra::aux::indices;

using json = nlohmann::json;

namespace gra {

// VEGAS execution stage controlling grid and integral state transitions
enum class VEGASStage { Adaptation, Integration, Generation };

// One VEGAS proposal point with its exact density and grid coordinates
struct VEGASSample {
  std::vector<double> point;
  std::vector<std::size_t> indices;
  double log_inverse_density = 0.0;
};

// Result of one proposal adaptation round or low-support retry
enum class VEGASAdaptationStatus {
  InsufficientSupport,
  Continue,
  Converged,
  RoundLimit
};

// Numerical diagnostics returned after one VEGAS adaptation batch
struct VEGASAdaptationReport {
  VEGASAdaptationStatus status = VEGASAdaptationStatus::Continue;
  std::size_t round = 0;
  std::size_t round_limit = 0;
  std::uint64_t evaluations = 0;
  std::uint64_t nonzero_calls = 0;
  std::size_t support = 0;
  std::size_t calls = 0;
  double effective_calls = 0.0;
  double ess_fraction = 0.0;
  double best_ess_fraction = 0.0;
  std::size_t best_round = 0;
};

// Adapted-grid retention policy after the final proposal round
enum class VEGASGridChoice { Last, Best };

// Compute the steering label for one adapted-grid retention policy
inline std::string VegasGridChoiceName(VEGASGridChoice choice) {
  if (choice == VEGASGridChoice::Last) {
    return "last";
  }
  if (choice == VEGASGridChoice::Best) {
    return "best";
  }
  throw std::invalid_argument("VegasGridChoiceName: unknown grid choice");
}

// Parse one adapted-grid retention policy from its steering label
inline VEGASGridChoice ParseVegasGridChoice(const std::string &choice) {
  if (choice == "last") {
    return VEGASGridChoice::Last;
  }
  if (choice == "best") {
    return VEGASGridChoice::Best;
  }
  throw std::invalid_argument(
      "ParseVegasGridChoice: choice must be 'last' or 'best'");
}

// Vegas MC default parameters
struct VEGASPARAM {
  unsigned int BINS =
      128;            // Maximum number of bins per dimension (EVEN NUMBER!)
  double ALPHA = 1.5; // Grid adaptation damping exponent
  double UNIFORM_MIX = 0.05;    // Uniform proposal fraction
  std::size_t MIN_SUPPORT = 10; // Minimum nonzero calls per adaptation round
  std::uint64_t MAX_SUPPORT_CALLS = 1000000; // Maximum low-support call target
  unsigned int MAX_ROUNDS = 24;              // Maximum adaptation rounds
  std::size_t CONVERGENCE_WINDOW = 3;        // ESS plateau window length
  double ESS_REL_TOLERANCE = 0.10;   // Relative ESS improvement tolerance
  VEGASGridChoice BEST_GRID_CHOICE =
      VEGASGridChoice::Last;          // Final adapted-grid retention policy
  bool AUTOMATIC_CONVERGENCE = true; // Extend adaptation until convergence

  // Proposal adaptation
  unsigned int NCALL = 20000; // Number of calls per adaptation round
  unsigned int ROUNDS = 5;    // Number of adaptation rounds
  int DEBUG = -1;             // Debug mode
};

// Test rolling VEGAS ESS improvement for convergence
// relative gain = [max(ESS_window)-ESS_first]/ESS_first
inline bool VegasAdaptationConverged(const std::vector<double> &ess_fractions,
                                     const VEGASPARAM &param) {
  const std::size_t window = param.CONVERGENCE_WINDOW;
  if (window < 2 || ess_fractions.size() < window) {
    return false;
  }

  const std::size_t begin = ess_fractions.size() - window;
  const double initial_ess = ess_fractions[begin];
  const double best_ess =
      *std::max_element(ess_fractions.begin() + begin, ess_fractions.end());
  if (!(initial_ess > 0.0) || !std::isfinite(initial_ess) ||
      !std::isfinite(best_ess)) {
    return false;
  }
  const double relative_gain =
      std::max(0.0, best_ess - initial_ess) / initial_ess;
  return relative_gain <= param.ESS_REL_TOLERANCE;
}

// Vegas MC adaptation data
struct VEGASData {

  // VEGAS initialization function
  void Init(VEGASStage stage, const VEGASPARAM &param) {
    ValidateParameters(param);
    if (FDIM == 0) {
      throw std::invalid_argument(
          "VEGASData::Init: phase-space dimension must be positive");
    }

    // First initialization: Create the grid and initial data
    if (stage == VEGASStage::Adaptation) {
      ClearAll(param);
      InitGridDependent(param);
    }
    // Start production integration from the adapted grid with reset integral
    // data
    else if (stage == VEGASStage::Integration) {
      ResetGrid();
      InitGridDependent(param);
    }
    // Generate events with the frozen grid, with or without integral data
    else if (stage == VEGASStage::Generation) {
      InitGridDependent(param);
    } else {
      throw std::invalid_argument("VEGASData::Init: unknown execution stage");
    }
  }

  // Create number of calls per thread, they need to sum to calls
  std::vector<unsigned int> GetLocalCalls(unsigned int calls,
                                          int N_threads) const {
    if (N_threads <= 0) {
      throw std::invalid_argument(
          "VEGASData::GetLocalCalls: N_threads must be positive");
    }
    const unsigned int thread_count = static_cast<unsigned int>(N_threads);
    std::vector<unsigned int> localcalls(thread_count, calls / thread_count);
    localcalls[0] += calls % thread_count;
    return localcalls;
  }

  // Initialize
  void InitGridDependent(const VEGASPARAM &param) {
    ValidateParameters(param);

    // If binning parameter changed from previous call
    if (param.BINS != BINS_prev) {
      if (BINS_prev == 0 || xmat.size() < BINS_prev ||
          !std::all_of(xmat.begin(), xmat.end(), [this](const auto &row) {
            return row.size() == FDIM;
          })) {
        throw std::logic_error("VEGASData::InitGridDependent: invalid source grid");
      }
      std::vector<std::vector<double>> grid(param.BINS, std::vector<double>(FDIM));
      const long double floor = 2.2204e-12L;
      const long double scale = 1.0L - param.BINS * floor;
      for (const auto &edge : indices(grid)) {
        const long double position =
            static_cast<long double>(edge + 1) * BINS_prev / param.BINS;
        const std::size_t source = std::min(
            static_cast<std::size_t>(position), static_cast<std::size_t>(BINS_prev - 1));
        const long double fraction = position - source;
        for (const auto &dimension : indices(grid[edge])) {
          const long double lower = source == 0 ? 0.0L : xmat[source - 1][dimension];
          const long double upper = xmat[source][dimension];
          grid[edge][dimension] = edge + 1 == param.BINS
              ? 1.0 : static_cast<double>((edge + 1) * floor + scale * (lower + fraction * (upper - lower)));
        }
      }
      xmat = std::move(grid);
      fmat.assign(param.BINS, std::vector<double>(FDIM, 0.0));
      f2mat.assign(param.BINS, std::vector<double>(FDIM, 0.0));
      BINS_prev = param.BINS;
      ResetGrid();
    }
  }

  // VEGAS grid optimizing function (algorithm adapted from Numerical Recipes)
  void OptimizeGrid(const VEGASPARAM &param) {
    ValidateParameters(param);

    std::vector<double> scores(param.BINS, 0.0);
    for (std::size_t j = 0; j < FDIM; ++j) {
      const double total = SmoothGridDimension(j, param.BINS);
      if (!(total > 0.0) || !std::isfinite(total)) {
        continue;
      }

      double score_sum = 0.0;
      for (std::size_t i = 0; i < param.BINS; ++i) {
        const double fraction = std::min(f2mat[i][j] / total, 1.0);
        scores[i] = GridAdaptationScore(fraction, param.ALPHA);
        score_sum += scores[i];
      }
      if (!(score_sum > 0.0) || !std::isfinite(score_sum)) {
        continue;
      }
      Rebin(score_sum / param.BINS, j, param.BINS, scores);
    }
  }

  // Set the phase-space dimension of the unit hypercube
  void SetDimension(unsigned int fdim) {
    if (fdim == 0) {
      throw std::invalid_argument(
          "VEGASData::SetDimension: dimension must be positive");
    }
    FDIM = fdim;
  }

  // Reset the scale-normalized grid moments for one VEGAS iteration
  void ResetGrid() {
    for (const auto &bin : indices(fmat)) {
      std::fill(fmat[bin].begin(), fmat[bin].end(), 0.0);
      std::fill(f2mat[bin].begin(), f2mat[bin].end(), 0.0);
    }
    grid_weight_scale = 0.0;
    grid_abs_weight_sum = 0.0L;
    grid_square_weight_sum = 0.0L;
    grid_calls = 0;
    grid_nonzero_calls = 0;
  }

  // Add one weight to each selected grid bin after common rescaling
  // f_{bd} += w/s, f2_{bd} += r_grid (w/s)^2 with s = max |w|
  void AccumulateGrid(double weight, const std::vector<std::size_t> &indices,
                      double grid_fraction = 1.0) {
    if (indices.size() != FDIM || !std::isfinite(weight) ||
        !(grid_fraction >= 0.0 && grid_fraction <= 1.0)) {
      throw std::invalid_argument(
          "VEGASData::AccumulateGrid: invalid weight, grid fraction or bin vector");
    }
    const bool valid_shape =
        fmat.size() == f2mat.size() && !fmat.empty() &&
        std::all_of(fmat.begin(), fmat.end(), [this](const auto &row) {
          return row.size() == FDIM;
        }) &&
        std::all_of(f2mat.begin(), f2mat.end(), [this](const auto &row) {
          return row.size() == FDIM;
        });
    if (!valid_shape) {
      throw std::logic_error(
          "VEGASData::AccumulateGrid: grid has invalid shape");
    }
    for (std::size_t dimension = 0; dimension < FDIM; ++dimension) {
      if (indices[dimension] == 0 || indices[dimension] > fmat.size()) {
        throw std::out_of_range(
            "VEGASData::AccumulateGrid: bin index outside the grid");
      }
    }
    if (grid_calls == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error(
          "VEGASData::AccumulateGrid: call counter overflow");
    }
    ++grid_calls;
    const double magnitude = std::abs(weight);
    if (math::IsZero(magnitude)) {
      return;
    }
    if (grid_nonzero_calls == std::numeric_limits<std::size_t>::max()) {
      throw std::overflow_error(
          "VEGASData::AccumulateGrid: support counter overflow");
    }
    ++grid_nonzero_calls;
    if (magnitude > grid_weight_scale) {
      const double ratio =
          grid_weight_scale > 0.0 ? grid_weight_scale / magnitude : 0.0;
      const long double extended_ratio = ratio;
      for (const auto &bin : gra::aux::indices(fmat)) {
        for (std::size_t dimension = 0; dimension < FDIM; ++dimension) {
          fmat[bin][dimension] *= ratio;
          f2mat[bin][dimension] *= ratio * ratio;
        }
      }
      grid_abs_weight_sum *= extended_ratio;
      grid_square_weight_sum *= extended_ratio * extended_ratio;
      grid_weight_scale = magnitude;
    }

    const double normalized = weight / grid_weight_scale;
    const long double normalized_magnitude = std::abs(normalized);
    const long double normalized_square =
        normalized_magnitude * normalized_magnitude;
    grid_abs_weight_sum += normalized_magnitude;
    grid_square_weight_sum += normalized_square;
    for (std::size_t dimension = 0; dimension < FDIM; ++dimension) {
      const std::size_t bin = indices[dimension] - 1;
      fmat[bin][dimension] += normalized;
      f2mat[bin][dimension] += grid_fraction * normalized * normalized;
    }
  }

  // Compute the effective sample fraction of the current grid round
  // ESS/N = (sum_i |w_i|)^2/[N sum_i w_i^2]
  double GridEffectiveFraction() const {
    return grid_calls > 0
               ? GridEffectiveCalls() / static_cast<double>(grid_calls)
               : 0.0;
  }

  // Compute the absolute-weight effective sample count of the current grid round
  // ESS = (sum_i |w_i|)^2/sum_i w_i^2
  double GridEffectiveCalls() const {
    if (!(grid_square_weight_sum > 0.0L)) {
      return 0.0;
    }
    const long double effective_calls =
        grid_abs_weight_sum * grid_abs_weight_sum / grid_square_weight_sum;
    return static_cast<double>(std::clamp(
        effective_calls, 0.0L, static_cast<long double>(grid_nonzero_calls)));
  }

  // Compute the calls accumulated in the current grid round
  std::size_t GridCalls() const { return grid_calls; }

  // Compute the nonzero calls accumulated in the current grid round
  std::size_t GridNonzeroCalls() const { return grid_nonzero_calls; }

  // Compute the validated VEGAS quantile edges in dimension-major order
  std::vector<std::vector<double>> QuantileEdges(std::size_t bins) const {
    if (FDIM == 0 || bins < 2 || xmat.size() != bins) {
      throw std::invalid_argument(
          "VEGASData::QuantileEdges: invalid dimension or bin count");
    }
    const bool valid_shape =
        std::all_of(xmat.begin(), xmat.end(),
                    [this](const auto &row) { return row.size() == FDIM; });
    if (!valid_shape) {
      throw std::invalid_argument(
          "VEGASData::QuantileEdges: invalid grid shape");
    }

    std::vector<std::vector<double>> edges(FDIM,
                                           std::vector<double>(bins + 1, 0.0));
    for (std::size_t dimension = 0; dimension < FDIM; ++dimension) {
      double previous = 0.0;
      for (std::size_t edge = 0; edge < bins; ++edge) {
        const double value = xmat[edge][dimension];
        if (!std::isfinite(value) || !(value > previous) || value > 1.0) {
          throw std::runtime_error(
              "VEGASData::QuantileEdges: grid edges must be finite and "
              "strictly increasing");
        }
        edges[dimension][edge + 1] = value;
        previous = value;
      }
      if (!math::IsExactEqual(previous, 1.0)) {
        throw std::runtime_error(
            "VEGASData::QuantileEdges: grid must end at one");
      }
    }
    return edges;
  }

  // Locate one point in the grid and return its log inverse density
  // log(1/q_grid) = sum_d log(B Delta x_{b_d,d})
  double GridLogInverseDensity(const std::vector<double> &point,
                               std::size_t bins,
                               std::vector<std::size_t> &indices) const {
    if (point.size() != FDIM || xmat.size() != bins || bins < 2) {
      throw std::invalid_argument(
          "VEGASData::GridLogInverseDensity: invalid point or grid");
    }
    const bool valid_shape =
        std::all_of(xmat.begin(), xmat.end(),
                    [this](const auto &row) { return row.size() >= FDIM; });
    if (!valid_shape) {
      throw std::invalid_argument(
          "VEGASData::GridLogInverseDensity: invalid grid shape");
    }
    indices.assign(FDIM, 0);
    double log_inverse_density = 0.0;
    for (std::size_t dimension = 0; dimension < FDIM; ++dimension) {
      if (!(point[dimension] >= 0.0 && point[dimension] <= 1.0) ||
          !std::isfinite(point[dimension])) {
        throw std::invalid_argument(
            "VEGASData::GridLogInverseDensity: point outside unit cube");
      }
      std::size_t bin = 0;
      while (bin + 1 < bins && point[dimension] >= xmat[bin][dimension]) {
        ++bin;
      }
      const double lower = bin == 0 ? 0.0 : xmat[bin - 1][dimension];
      const double upper = xmat[bin][dimension];
      const double width = upper - lower;
      if (!(width > 0.0) || !std::isfinite(width)) {
        throw std::runtime_error(
            "VEGASData::GridLogInverseDensity: non-positive grid width");
      }
      indices[dimension] = bin + 1;
      log_inverse_density += std::log(width * bins);
    }
    return log_inverse_density;
  }

  // Compute the exact log inverse density of a grid and uniform mixture
  // q_mix = (1-epsilon)q_grid + epsilon
  static double MixedLogInverseDensity(double grid_log_inverse_density,
                                       double uniform_mix) {
    if (!(uniform_mix > 0.0 && uniform_mix < 1.0) ||
        !std::isfinite(grid_log_inverse_density)) {
      throw std::invalid_argument(
          "VEGASData::MixedLogInverseDensity: invalid mixture input");
    }
    const double grid_term =
        std::log1p(-uniform_mix) - grid_log_inverse_density;
    const double uniform_term = std::log(uniform_mix);
    const double maximum = std::max(grid_term, uniform_term);
    const double log_density =
        maximum + std::log(std::exp(grid_term - maximum) +
                           std::exp(uniform_term - maximum));
    return -log_density;
  }

  // Full init
  void ClearAll(const VEGASPARAM &param) {
    // Matrices [BINS x FDIM]
    fmat = std::vector<std::vector<double>>(param.BINS,
                                            std::vector<double>(FDIM, 0.0));
    f2mat = std::vector<std::vector<double>>(param.BINS,
                                             std::vector<double>(FDIM, 0.0));
    xmat = std::vector<std::vector<double>>(param.BINS,
                                            std::vector<double>(FDIM, 0.0));

    // Init with 1!
    for (std::size_t j = 0; j < FDIM; ++j) {
      xmat[0][j] = 1.0;
    }

    ResetGrid();

    // Previous binning
    BINS_prev = 1;
  }

  // Save everything to JSON
  void struct2json(json &j) const {
    j["BINS_prev"] = BINS_prev;

    j["FDIM"] = FDIM;

    j["fmat"] = fmat;
    j["f2mat"] = f2mat;
    j["xmat"] = xmat;
  }

  // Read everything from JSON
  void json2struct(const json &jf) {
    jf.at("BINS_prev").get_to(BINS_prev);

    jf.at("FDIM").get_to(FDIM);

    jf.at("fmat").get_to(fmat);
    jf.at("f2mat").get_to(f2mat);
    jf.at("xmat").get_to(xmat);
  }
  // VEGAS scalars
  unsigned int BINS_prev = 0;
  unsigned int FDIM = 0;

  double grid_weight_scale = 0.0;
  long double grid_abs_weight_sum = 0.0L;
  long double grid_square_weight_sum = 0.0L;
  std::size_t grid_calls = 0;
  std::size_t grid_nonzero_calls = 0;

  // Matrices
  std::vector<std::vector<double>> fmat;
  std::vector<std::vector<double>> f2mat;
  std::vector<std::vector<double>> xmat;

private:
  // Validate technical controls before allocating or adapting the VEGAS map
  static void ValidateParameters(const VEGASPARAM &param) {
    if (param.BINS < 4 || (param.BINS % 2) != 0) {
      throw std::invalid_argument(
          "VEGASData: the bin count must be even and at least four");
    }
    if (param.MIN_SUPPORT == 0 || param.MAX_SUPPORT_CALLS < param.MIN_SUPPORT) {
      throw std::invalid_argument(
          "VEGASData: invalid adaptation support limits");
    }
    if (param.CONVERGENCE_WINDOW < 2 || param.MAX_ROUNDS < param.ROUNDS ||
        param.MAX_ROUNDS < param.CONVERGENCE_WINDOW) {
      throw std::invalid_argument(
          "VEGASData: invalid adaptation convergence round limits");
    }
    if (!(param.ESS_REL_TOLERANCE > 0.0) || !(param.ESS_REL_TOLERANCE <= 1.0) ||
        !std::isfinite(param.ESS_REL_TOLERANCE)) {
      throw std::invalid_argument(
          "VEGASData: ESS convergence tolerance must be positive and finite");
    }
    if (!std::isfinite(param.ALPHA) || param.ALPHA < 0.0 ||
        !std::isfinite(param.UNIFORM_MIX) || param.UNIFORM_MIX <= 0.0 ||
        param.UNIFORM_MIX >= 1.0 || param.NCALL == 0 || param.ROUNDS == 0) {
      throw std::invalid_argument(
          "VEGASData: invalid adaptation exponent, mixture or call controls");
    }
  }

  // Redistribute one grid dimension according to nonnegative bin scores
  void Rebin(double target, std::size_t dimension, std::size_t bins,
             const std::vector<double> &scores) {
    if (!(target > 0.0) || bins < 2 || scores.size() < bins) {
      throw std::invalid_argument("VEGASData::Rebin: invalid rebin input");
    }

    std::vector<double> new_grid(bins, 0.0);
    std::size_t source = 0;
    double accumulated = 0.0;
    double lower = 0.0;
    for (std::size_t edge = 0; edge + 1 < bins; ++edge) {
      while (target > accumulated && source < bins) {
        accumulated += scores[source];
        ++source;
      }
      if (target > accumulated || source == 0 || !(scores[source - 1] > 0.0)) {
        throw std::runtime_error(
            "VEGASData::Rebin: bin scores do not cover the target");
      }

      if (source > 1) {
        lower = xmat[source - 2][dimension];
      }
      const double upper = xmat[source - 1][dimension];
      accumulated -= target;
      new_grid[edge] =
          upper - (upper - lower) * accumulated / scores[source - 1];
    }
    new_grid[bins - 1] = 1.0;

    // Map raw widths onto the lower-bounded unit simplex
    const long double floor = 2.2204e-12L;
    const long double scale = 1.0L - static_cast<long double>(bins) * floor;
    long double raw_lower = 0.0L;
    long double protected_upper = 0.0L;
    for (std::size_t edge = 0; edge < bins; ++edge) {
      const long double raw_upper = new_grid[edge];
      const long double raw_width = raw_upper - raw_lower;
      if (!(raw_width >= 0.0L) || !std::isfinite(raw_width)) {
        throw std::runtime_error(
            "VEGASData::Rebin: generated an invalid raw bin width");
      }
      protected_upper += floor + scale * raw_width;
      xmat[edge][dimension] =
          edge + 1 == bins ? 1.0 : static_cast<double>(protected_upper);
      raw_lower = raw_upper;
    }
  }

  // Smooth one dimension of the squared-weight grid and return its sum
  // f2'_i = (f2_{i-1}+f2_i+f2_{i+1})/3 with half-weight boundaries
  double SmoothGridDimension(std::size_t dimension, std::size_t bins) {
    for (std::size_t i = 0; i < bins; ++i) {
      const double value = f2mat[i][dimension];
      if (!(value >= 0.0) || !std::isfinite(value)) {
        return std::numeric_limits<double>::quiet_NaN();
      }
    }

    std::vector<double> smoothed(bins, 0.0);
    smoothed[0] = 0.5 * f2mat[0][dimension] + 0.5 * f2mat[1][dimension];
    for (std::size_t i = 1; i + 1 < bins; ++i) {
      smoothed[i] = f2mat[i - 1][dimension] / 3.0 + f2mat[i][dimension] / 3.0 +
                    f2mat[i + 1][dimension] / 3.0;
    }
    smoothed[bins - 1] =
        0.5 * f2mat[bins - 2][dimension] + 0.5 * f2mat[bins - 1][dimension];

    long double total = 0.0L;
    for (std::size_t i = 0; i < bins; ++i) {
      f2mat[i][dimension] = smoothed[i];
      total += smoothed[i];
    }
    return static_cast<double>(total);
  }

  // Compute the damped logarithmic-mean score for one normalized grid bin
  // score(r) = [(1-r)/(-log r)]^alpha
  static double GridAdaptationScore(double fraction, double alpha) {
    if (math::IsZero(alpha)) {
      return 1.0;
    }
    if (!(fraction > 0.0)) {
      return 0.0;
    }
    if (fraction >= 1.0) {
      return 1.0;
    }

    const double complement = 1.0 - fraction;
    const double logarithm = -std::log(fraction);
    const double score = complement / logarithm;
    const double damped = std::pow(score, alpha);
    return (damped >= 0.0 && std::isfinite(damped)) ? damped : 0.0;
  }

}; // struct VEGASData

// Standalone adaptive VEGAS importance sampler
class MVEGASIntegrator {
public:
  // Set the unit-hypercube dimension before initializing a stage
  void SetDimension(std::size_t dimension) {
    if (dimension == 0 ||
        dimension > std::numeric_limits<unsigned int>::max()) {
      throw std::invalid_argument(
          "MVEGASIntegrator::SetDimension: invalid dimension");
    }
    if (dimension != data_.FDIM) {
      data_ = VEGASData{};
      proposal_frozen_ = false;
    }
    data_.SetDimension(static_cast<unsigned int>(dimension));
  }

  // Initialize one adaptation, integration or generation stage
  void Init(VEGASStage stage, std::uint64_t calls, unsigned int rounds,
            const VEGASPARAM &param) {
    if (calls == 0) {
      throw std::invalid_argument(
          "MVEGASIntegrator::Init: calls must be positive");
    }
    if (stage == VEGASStage::Adaptation && rounds == 0) {
      throw std::invalid_argument(
          "MVEGASIntegrator::Init: adaptation rounds must be positive");
    }
    if (stage != VEGASStage::Adaptation && !proposal_frozen_) {
      throw std::logic_error(
          "MVEGASIntegrator::Init: production requires an adapted proposal");
    }
    VEGASData initialized = data_;
    initialized.Init(stage, param);
    data_ = std::move(initialized);
    stage_ = stage;
    param_ = param;
    initial_calls_ = calls;
    minimum_rounds_ = rounds;
    if (stage_ == VEGASStage::Adaptation) {
      proposal_frozen_ = false;
      ResetAdaptation();
    }
  }

  // Compute the call count still needed by the current adaptation target
  std::uint64_t NextAdaptationCalls() const {
    RequireStage(VEGASStage::Adaptation);
    const std::uint64_t accumulated_calls =
        reset_adaptation_grid_ ? 0 : data_.GridCalls();
    if (accumulated_calls >= adaptation_calls_) {
      throw std::logic_error(
          "MVEGASIntegrator::NextAdaptationCalls: accumulated calls exceed "
          "the adaptation target");
    }
    return adaptation_calls_ - accumulated_calls;
  }

  // Reset grid accumulators before one proposal-adaptation batch
  void BeginAdaptationBatch() {
    RequireStage(VEGASStage::Adaptation);
    if (reset_adaptation_grid_) {
      data_.ResetGrid();
    }
  }

  // Draw one exact VEGAS and uniform mixture proposal point
  template <typename UniformRandom>
  VEGASSample Sample(UniformRandom &&uniform_random) const {
    if (data_.FDIM == 0 || data_.xmat.size() != param_.BINS) {
      throw std::logic_error(
          "MVEGASIntegrator::Sample: proposal is not initialized");
    }

    auto &&uniform_source = uniform_random;
    auto draw_uniform = [&uniform_source]() {
      const double value = uniform_source();
      if (!std::isfinite(value) || value < 0.0 || value >= 1.0) {
        throw std::invalid_argument(
            "MVEGASIntegrator::Sample: random draw outside [0,1)");
      }
      return value;
    };
    VEGASSample sample;
    sample.point.assign(data_.FDIM, 0.0);
    sample.indices.assign(data_.FDIM, 0);
    double grid_log_inverse_density = 0.0;

    if (draw_uniform() < param_.UNIFORM_MIX) {
      for (double &coordinate : sample.point) {
        coordinate = draw_uniform();
      }
      grid_log_inverse_density = data_.GridLogInverseDensity(
          sample.point, param_.BINS, sample.indices);
    } else {
      grid_log_inverse_density =
          DrawGridPoint(draw_uniform, sample.point, sample.indices);
    }
    sample.log_inverse_density = VEGASData::MixedLogInverseDensity(
        grid_log_inverse_density, param_.UNIFORM_MIX);
    return sample;
  }

  // Accumulate the variance gradient with the grid probability in the mixture
  // r_grid = (1-epsilon) q_grid/q_mix, d< w^2 >/d log q_grid = -< r_grid w^2 >
  void AccumulateAdaptation(double weight,
                            const std::vector<std::size_t> &bins) {
    RequireStage(VEGASStage::Adaptation);
    double log_inverse = 0.0;
    for (const auto &dimension : indices(bins)) {
      const std::size_t bin = bins[dimension] - 1;
      const double lower = bin == 0 ? 0.0 : data_.xmat.at(bin - 1).at(dimension);
      const double width = data_.xmat.at(bin).at(dimension) - lower;
      log_inverse += std::log(width * param_.BINS);
    }
    const double grid_fraction = std::exp(
        std::log1p(-param_.UNIFORM_MIX) - log_inverse +
        VEGASData::MixedLogInverseDensity(log_inverse, param_.UNIFORM_MIX));
    data_.AccumulateGrid(weight, bins, grid_fraction);
  }

  // Complete one adaptation batch and update the proposal grid
  VEGASAdaptationReport
  FinishAdaptationBatch(std::uint64_t batch_calls) {
    RequireStage(VEGASStage::Adaptation);
    if (data_.GridCalls() != adaptation_calls_) {
      throw std::logic_error(
          "MVEGASIntegrator::FinishAdaptationBatch: sampled call count does "
          "not match the adaptation target");
    }
    if (batch_calls >
        std::numeric_limits<std::uint64_t>::max() - adaptation_evaluations_) {
      throw std::overflow_error(
          "MVEGASIntegrator::FinishAdaptationBatch: evaluation counter overflow");
    }
    adaptation_evaluations_ += batch_calls;

    VEGASAdaptationReport report = CurrentAdaptationReport();
    if (data_.GridNonzeroCalls() < param_.MIN_SUPPORT) {
      report.status = VEGASAdaptationStatus::InsufficientSupport;
      if (data_.GridNonzeroCalls() > 0) {
        data_.OptimizeGrid(param_);
      }
      IncreaseAdaptationCalls(report);
      return report;
    }

    AcceptAdaptationRound();
    report = CurrentAdaptationReport();
    report.round = completed_adaptation_rounds_;
    report.status = AdaptationStatus();
    if (param_.BEST_GRID_CHOICE == VEGASGridChoice::Best &&
        (report.status == VEGASAdaptationStatus::Converged ||
         report.status == VEGASAdaptationStatus::RoundLimit)) {
      data_.xmat = best_adaptation_grid_;
    }
    if (report.status == VEGASAdaptationStatus::Converged ||
        report.status == VEGASAdaptationStatus::RoundLimit) {
      proposal_frozen_ = true;
    }
    return report;
  }

  // Compute whether the latest accepted adaptation round completes the stage
  bool AdaptationComplete() const {
    const VEGASAdaptationStatus status = AdaptationStatus();
    return status == VEGASAdaptationStatus::Converged ||
           status == VEGASAdaptationStatus::RoundLimit;
  }

  // Compute the accepted nonzero calls accumulated during adaptation
  std::uint64_t AdaptationNonzeroCalls() const {
    return adaptation_nonzero_calls_;
  }

  // Compute the maximum number of accepted adaptation rounds
  std::size_t AdaptationRoundLimit() const {
    return param_.AUTOMATIC_CONVERGENCE ? param_.MAX_ROUNDS : minimum_rounds_;
  }

  // Access the configured unit-hypercube dimension
  std::size_t Dimension() const { return data_.FDIM; }

  // Split one batch exactly across a positive worker count
  std::vector<unsigned int> LocalCalls(unsigned int calls,
                                       int workers) const {
    return data_.GetLocalCalls(calls, workers);
  }

  // Compute dimension-major quantile edges from the adapted proposal
  std::vector<std::vector<double>> QuantileEdges() const {
    return data_.QuantileEdges(param_.BINS);
  }

  // Check whether adaptation has produced a frozen proposal
  bool IsFrozen() const { return proposal_frozen_; }

  // Access immutable VEGAS numerical steering
  const VEGASPARAM &Config() const { return param_; }

  // Serialize the frozen proposal and its numerical steering
  json Serialize() const {
    if (!IsFrozen()) {
      throw std::logic_error(
          "MVEGASIntegrator::Serialize: proposal is not frozen");
    }
    json document;
    data_.struct2json(document);
    document["uniform_mix"] = param_.UNIFORM_MIX;
    document["min_support"] = param_.MIN_SUPPORT;
    document["max_support_calls"] = param_.MAX_SUPPORT_CALLS;
    document["max_rounds"] = param_.MAX_ROUNDS;
    document["convergence_window"] = param_.CONVERGENCE_WINDOW;
    document["ess_rel_tolerance"] = param_.ESS_REL_TOLERANCE;
    document["best_grid_choice"] =
        VegasGridChoiceName(param_.BEST_GRID_CHOICE);
    document["automatic_convergence"] = param_.AUTOMATIC_CONVERGENCE;
    return document;
  }

  // Restore one frozen proposal after validating its dimension and steering
  void Deserialize(const json &document, std::size_t expected_dimension,
                   const VEGASPARAM &param) {
    VEGASData restored;
    restored.json2struct(document);
    if (restored.BINS_prev != param.BINS) {
      throw std::invalid_argument(
          "MVEGASIntegrator::Deserialize: bin count does not match");
    }
    if (restored.FDIM != expected_dimension) {
      throw std::invalid_argument(
          "MVEGASIntegrator::Deserialize: dimension does not match");
    }
    if (!SerializedParametersMatch(document, param)) {
      throw std::invalid_argument(
          "MVEGASIntegrator::Deserialize: numerical steering does not match");
    }
    restored.QuantileEdges(param.BINS);
    for (const auto *matrix : {&restored.fmat, &restored.f2mat}) {
      if (matrix->size() != param.BINS ||
          !std::all_of(matrix->begin(), matrix->end(), [&restored](const auto &row) {
            return row.size() == restored.FDIM &&
                   std::all_of(row.begin(), row.end(), [](double value) { return std::isfinite(value); });
          })) {
        throw std::invalid_argument("MVEGASIntegrator::Deserialize: invalid grid moments");
      }
    }
    restored.Init(VEGASStage::Generation, param);
    data_ = std::move(restored);
    param_ = param;
    stage_ = VEGASStage::Generation;
    proposal_frozen_ = true;
  }

  // Compute mutable numerical data for serialization and diagnostics
  VEGASData &Data() { return data_; }

  // Access immutable numerical data for sampling diagnostics
  const VEGASData &Data() const { return data_; }

private:
  // Require one active VEGAS execution stage
  void RequireStage(VEGASStage expected) const {
    if (stage_ != expected) {
      throw std::logic_error("MVEGASIntegrator: unexpected execution stage");
    }
  }

  // Reset all proposal adaptation bookkeeping
  void ResetAdaptation() {
    adaptation_evaluations_ = 0;
    adaptation_nonzero_calls_ = 0;
    adaptation_calls_ = initial_calls_;
    completed_adaptation_rounds_ = 0;
    reset_adaptation_grid_ = true;
    adaptation_ess_fractions_.clear();
    best_adaptation_grid_.clear();
    best_adaptation_ess_ = -1.0;
    best_adaptation_round_ = 0;
  }

  // Compute diagnostics for the current adaptation grid
  VEGASAdaptationReport CurrentAdaptationReport() const {
    VEGASAdaptationReport report;
    report.round = completed_adaptation_rounds_ + 1;
    report.round_limit = AdaptationRoundLimit();
    report.evaluations = adaptation_evaluations_;
    report.nonzero_calls = adaptation_nonzero_calls_;
    report.support = data_.GridNonzeroCalls();
    report.calls = data_.GridCalls();
    report.effective_calls = data_.GridEffectiveCalls();
    report.ess_fraction = data_.GridEffectiveFraction();
    report.best_ess_fraction = best_adaptation_ess_;
    report.best_round = best_adaptation_round_;
    return report;
  }

  // Increase the cumulative call target after a low-support batch
  void IncreaseAdaptationCalls(const VEGASAdaptationReport &report) {
    if (adaptation_calls_ >= param_.MAX_SUPPORT_CALLS) {
      throw std::runtime_error(
          "MVEGASIntegrator: insufficient adaptation support after " +
          std::to_string(adaptation_calls_) + " calls: support " +
          std::to_string(report.support) + "/" +
          std::to_string(report.calls) + ", weight ESS calls " +
          std::to_string(report.effective_calls) + ", required support " +
          std::to_string(param_.MIN_SUPPORT) +
          ". Check generation and fiducial cuts.");
    }
    const std::uint64_t doubled =
        adaptation_calls_ > std::numeric_limits<std::uint64_t>::max() / 2
            ? std::numeric_limits<std::uint64_t>::max()
            : adaptation_calls_ * 2;
    adaptation_calls_ = std::min(doubled, param_.MAX_SUPPORT_CALLS);
    reset_adaptation_grid_ = data_.GridNonzeroCalls() > 0;
  }

  // Accept one supported adaptation round and optimize the grid
  void AcceptAdaptationRound() {
    const double ess_fraction = data_.GridEffectiveFraction();
    ++completed_adaptation_rounds_;
    if (ess_fraction > best_adaptation_ess_) {
      best_adaptation_ess_ = ess_fraction;
      best_adaptation_grid_ = data_.xmat;
      best_adaptation_round_ = completed_adaptation_rounds_;
    }
    if (data_.GridNonzeroCalls() >
        std::numeric_limits<std::uint64_t>::max() -
            adaptation_nonzero_calls_) {
      throw std::overflow_error(
          "MVEGASIntegrator::AcceptAdaptationRound: support counter overflow");
    }
    adaptation_nonzero_calls_ += data_.GridNonzeroCalls();
    data_.OptimizeGrid(param_);
    adaptation_ess_fractions_.push_back(ess_fraction);
    reset_adaptation_grid_ = true;
  }

  // Compute the convergence state after the latest supported round
  VEGASAdaptationStatus AdaptationStatus() const {
    const bool minimum_complete =
        completed_adaptation_rounds_ >= minimum_rounds_;
    if (param_.AUTOMATIC_CONVERGENCE && minimum_complete &&
        VegasAdaptationConverged(adaptation_ess_fractions_, param_)) {
      return VEGASAdaptationStatus::Converged;
    }
    if (completed_adaptation_rounds_ >= AdaptationRoundLimit()) {
      return VEGASAdaptationStatus::RoundLimit;
    }
    return VEGASAdaptationStatus::Continue;
  }

  // Compute whether serialized steering matches the active VEGAS controls
  static bool SerializedParametersMatch(const json &document,
                                        const VEGASPARAM &param) {
    return math::IsExactEqual(document.at("uniform_mix").get<double>(),
                              param.UNIFORM_MIX) &&
           document.at("min_support").get<std::size_t>() ==
               param.MIN_SUPPORT &&
           document.at("max_support_calls").get<std::uint64_t>() ==
               param.MAX_SUPPORT_CALLS &&
           document.at("max_rounds").get<unsigned int>() ==
               param.MAX_ROUNDS &&
           document.at("convergence_window").get<std::size_t>() ==
               param.CONVERGENCE_WINDOW &&
           math::IsExactEqual(
               document.at("ess_rel_tolerance").get<double>(),
               param.ESS_REL_TOLERANCE) &&
           document.at("best_grid_choice").get<std::string>() ==
               VegasGridChoiceName(param.BEST_GRID_CHOICE) &&
           document.at("automatic_convergence").get<bool>() ==
               param.AUTOMATIC_CONVERGENCE;
  }

  // Draw one point from the separable adapted grid
  template <typename UniformRandom>
  double DrawGridPoint(UniformRandom &uniform_random,
                       std::vector<double> &point,
                       std::vector<std::size_t> &indices) const {
    double log_inverse_density = 0.0;
    for (std::size_t dimension = 0; dimension < data_.FDIM; ++dimension) {
      const double bin_coordinate = uniform_random() * param_.BINS + 1.0;
      indices[dimension] = std::max(
          1U, std::min(static_cast<unsigned int>(bin_coordinate), param_.BINS));
      const std::size_t bin = indices[dimension] - 1;
      const double lower = bin == 0 ? 0.0 : data_.xmat[bin - 1][dimension];
      const double upper = data_.xmat[bin][dimension];
      const double width = upper - lower;
      if (!(width > 0.0) || !std::isfinite(width)) {
        throw std::runtime_error(
            "MVEGASIntegrator::DrawGridPoint: non-positive grid width");
      }
      // Keep rounding inside the bin whose density is used for this draw
      point[dimension] = std::min(
          lower + (bin_coordinate - indices[dimension]) * width,
          std::nextafter(upper, lower));
      log_inverse_density += std::log(width * param_.BINS);
    }
    return log_inverse_density;
  }

  VEGASData data_;
  VEGASPARAM param_;
  VEGASStage stage_ = VEGASStage::Adaptation;

  std::uint64_t initial_calls_ = 0;
  std::uint64_t adaptation_calls_ = 0;
  std::uint64_t adaptation_evaluations_ = 0;
  std::uint64_t adaptation_nonzero_calls_ = 0;
  std::size_t minimum_rounds_ = 0;
  std::size_t completed_adaptation_rounds_ = 0;
  bool reset_adaptation_grid_ = true;
  bool proposal_frozen_ = false;
  std::vector<double> adaptation_ess_fractions_;
  std::vector<std::vector<double>> best_adaptation_grid_;
  double best_adaptation_ess_ = -1.0;
  std::size_t best_adaptation_round_ = 0;
};

} // namespace gra

#endif
