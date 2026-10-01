// GRANIITTI Monte Carlo main class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGRANIITTI_H
#define MGRANIITTI_H

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <exception>
#include <limits>
#include <mutex>
#include <random>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

// HepMC3
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAscii.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Kinematics/MHardDiffraction.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Tech/MLHE.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Sampling/MNeuroJac.h"
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Tech/MTimer.h"
#include "Graniitti/Sampling/MVEGAS.h"

// Other
#include "cxxopts.hpp"
#include "json.hpp"

namespace gra {

// Common frozen-proposal integration controls
struct INTEGRATIONPARAM {
  std::uint64_t MIN_SAMPLES = 1000000;
  std::uint64_t MAX_SAMPLES = 100000000;
  double PRECISION = 0.01;
};

// Maximum-weight calibration modes
enum class MaxWeightMode { Estimate, Fixed };

// Generation actions after an overweight trial
enum class OverflowAction { KeepAsWeighted, Break };

// Common maximum-weight calibration controls
struct MAXWEIGHTPARAM {
  MaxWeightMode MODE = MaxWeightMode::Estimate;
  double MAX_W_VALUE = 0.0;
  double SAFETY_FACTOR = 1.5;
  double MAX_OVERFLOW_PROB = 1e-5;
  double CONFIDENCE_LEVEL = 0.999;
  OverflowAction OVERFLOW_ACTION = OverflowAction::KeepAsWeighted;
};

// Persisted maximum-weight calibration state
struct MAXWEIGHTSTATE {
  double envelope = 0.0;
  double largest_overflow = 0.0;
  double generation_envelope = 0.0;
  std::uint64_t validation_samples = 0;
  std::uint64_t envelope_updates = 0;
  bool initialized = false;
  bool generation_active = false;
  bool generation_overflow = false;

  // Reset transient state for one generation attempt
  void ResetGeneration() {
    if (generation_active) {
      throw std::logic_error(
          "MAX_WEIGHTSTATE: event generation is already active");
    }
    largest_overflow = 0.0;
    generation_envelope = envelope;
    generation_active = true;
    generation_overflow = false;
  }

  // End the interval in which the generation envelope is immutable
  void EndGeneration() { generation_active = false; }

  // Compute the envelope only while its generation snapshot is unchanged
  double GenerationEnvelope() const {
    if (!generation_active) {
      throw std::logic_error(
          "MAX_WEIGHT: generation envelope requested outside generation");
    }
    if (!initialized || !(generation_envelope > 0.0) ||
        !std::isfinite(generation_envelope)) {
      throw std::logic_error(
          "MAX_WEIGHT: generation requires an initialized finite envelope");
    }
    if (!math::IsExactEqual(envelope, generation_envelope)) {
      throw std::logic_error(
          "MAX_WEIGHT: envelope changed after event generation started");
    }
    return generation_envelope;
  }

  // Record one generation overflow without changing the frozen envelope
  bool ObserveGenerationOverflow(double weight) {
    const bool first_overflow = !generation_overflow;
    largest_overflow = std::max(largest_overflow, weight);
    generation_overflow = true;
    return first_overflow;
  }
};

// Sampling stage controlling integral and maximum weight observations
enum class SamplingStage { Adaptation, Integration, Generation };

// Summary of one independent integration batch
struct IntegrationBatchSummary {
  long double mean = 0.0L;
  long double variance_of_mean = 0.0L;
  double ess_fraction = 0.0;
  std::size_t samples = 0;
  bool valid = false;
};

// Integration statistics
class Stats {
public:
  // Reset counters that belong to one event-output run
  void ResetGeneration() {
    trials = 0.0;
    generated = 0;
    N_overflow = 0;
  }

private:
  // Accumulate one sequential cut-flow observation
  void AccumulateCutFlow(const gra::MEventWeightState &aux) {
    evaluations += 1.0;

    const bool kinematics_pass = aux.kinematics_ok;
    const bool fidcuts_pass = kinematics_pass && aux.fidcuts_ok;
    const bool vetocuts_pass = fidcuts_pass && aux.vetocuts_ok;
    const bool amplitude_evaluated =
        vetocuts_pass && (!aux.technical_failure || aux.amplitude_failure);
    const bool amplitude_pass = amplitude_evaluated && aux.amplitude_ok &&
                                !aux.amplitude_failure &&
                                !aux.technical_failure;

    kinematics_ok += kinematics_pass ? 1.0 : 0.0;
    fidcuts_ok += fidcuts_pass ? 1.0 : 0.0;
    vetocuts_ok += vetocuts_pass ? 1.0 : 0.0;
    amplitude_evaluations += amplitude_evaluated ? 1.0 : 0.0;
    amplitude_ok += amplitude_pass ? 1.0 : 0.0;
    amplitude_failures += aux.amplitude_failure ? 1.0 : 0.0;
    technical_failures += aux.technical_failure ? 1.0 : 0.0;
    all_ok += amplitude_pass ? 1.0 : 0.0;
  }

public:
  // Observe one phase-space sample through the common integrator interface
  // w = f(x)/p(x)
  double ObserveLogSample(const gra::MEventWeightState &aux, double raw_weight,
                          double log_inverse_proposal_density,
                          SamplingStage stage) {
    if (stage != SamplingStage::Integration) {
      const double weight = statistics::ImportanceWeight(
          raw_weight, log_inverse_proposal_density);
      if (stage == SamplingStage::Generation) {
        ContinueCrossSection(weight);
      }
      return weight;
    }

    const statistics::ImportanceSample importance_sample =
        statistics::EvaluateImportanceSample(raw_weight,
                                             log_inverse_proposal_density);
    AccumulateCutFlow(aux);
    ++integration_samples;
    integrand_weight_moments.Add(raw_weight);
    if (importance_sample.HasMaterializedProposalDensity()) {
      proposal_density_moments.Add(importance_sample.proposal_density);
    } else {
      proposal_density_moments.AddLogPositive(
          -static_cast<long double>(log_inverse_proposal_density));
    }
    importance_weight_moments.Add(importance_sample.weight);

    integration_batch_moments.Add(importance_sample.weight);
    if (continued_weight_moments.Count() > 0) {
      continued_weight_moments.Add(importance_sample.weight);
    }
    return importance_sample.weight;
  }

  // Compute a finite conditional cut-flow rate
  static double ConditionalRate(double passed, double tested) {
    return tested > 0.0 ? passed / tested : 0.0;
  }

  // Compute the maximum relative to the frozen-sample importance-weight mean
  // ratio = w_max/<w>
  double MaximumToMean(double maximum) const {
    const long double mean = SamplingWeightMean();
    return mean > 0.0L && maximum > 0.0 && std::isfinite(mean)
               ? static_cast<double>(static_cast<long double>(maximum) / mean)
               : std::numeric_limits<double>::infinity();
  }

  // Estimate the rejection-unweighting efficiency from the sampled envelope
  // epsilon = min(<w>/w_max,1)
  double EstimatedUnweightingEfficiency(double maximum) const {
    const double ratio = MaximumToMean(maximum);
    return std::isfinite(ratio) && ratio > 0.0 ? std::min(1.0 / ratio, 1.0)
                                               : 0.0;
  }

  // Compute the effective fraction of frozen-proposal importance weights
  double SamplingEffectiveFraction() const {
    return importance_weight_moments.EffectiveSampleFraction();
  }

  // Compute the arithmetic mean of frozen-proposal importance weights
  long double SamplingWeightMean() const {
    return importance_weight_moments.Mean();
  }

  // Compute the largest observed frozen-proposal importance weight
  long double SamplingWeightMaximum() const {
    return importance_weight_moments.Maximum();
  }

  // Compute the number of frozen-proposal importance weights
  std::size_t SamplingWeightCount() const {
    return importance_weight_moments.Count();
  }

  // Compute the number of positive frozen-proposal importance weights
  std::size_t SamplingWeightPositiveCount() const {
    return importance_weight_moments.PositiveCount();
  }

  // Compute the importance-weight sample standard deviation divided by its mean
  long double SamplingWeightRelativeStd() const {
    return importance_weight_moments.RelativeStandardDeviation();
  }

  // Compute the smallest positive frozen-proposal importance weight
  long double SamplingWeightMinimumPositive() const {
    return importance_weight_moments.MinimumPositive();
  }

  // Compute the number of frozen-proposal integrand observations
  std::size_t IntegrandWeightCount() const {
    return integrand_weight_moments.Count();
  }

  // Compute the number of positive frozen-proposal integrand observations
  std::size_t IntegrandPositiveCount() const {
    return integrand_weight_moments.PositiveCount();
  }

  // Compute the mean raw integrand weight
  long double IntegrandWeightMean() const {
    return integrand_weight_moments.Mean();
  }

  // Compute the raw integrand sample standard deviation divided by its mean
  long double IntegrandWeightRelativeStd() const {
    return integrand_weight_moments.RelativeStandardDeviation();
  }

  // Compute the smallest positive raw integrand weight
  long double IntegrandWeightMinimumPositive() const {
    return integrand_weight_moments.MinimumPositive();
  }

  // Compute the largest raw integrand weight
  long double IntegrandWeightMaximum() const {
    return integrand_weight_moments.Maximum();
  }

  // Compute the number of frozen proposal-density observations
  std::size_t ProposalDensityCount() const {
    return proposal_density_moments.Count();
  }

  // Compute the number of positive frozen proposal-density observations
  std::size_t ProposalDensityPositiveCount() const {
    return proposal_density_moments.PositiveCount();
  }

  // Compute the mean frozen proposal density
  long double ProposalDensityMean() const {
    return proposal_density_moments.Mean();
  }

  // Compute the proposal-density sample standard deviation divided by its mean
  long double ProposalDensityRelativeStd() const {
    return proposal_density_moments.RelativeStandardDeviation();
  }

  // Compute the smallest positive frozen proposal density
  long double ProposalDensityMinimumPositive() const {
    return proposal_density_moments.MinimumPositive();
  }

  // Compute the largest frozen proposal density
  long double ProposalDensityMaximum() const {
    return proposal_density_moments.Maximum();
  }

  // Calculate the pooled cross section and standard error for every integrator
  // sigma = <w>, delta sigma = s_w/sqrt(N)
  void CalculateCrossSection() {
    const statistics::ScaledPositiveMoments &moments = CrossSectionMoments();
    const std::size_t calls = moments.Count();
    const long double mean = moments.Mean();
    if (calls < 2) {
      sigma = calls > 0 ? static_cast<double>(mean) : 0.0;
      sigma_err = std::numeric_limits<double>::infinity();
      return;
    }
    const long double relative_std = moments.RelativeStandardDeviation();
    if (!std::isfinite(mean) || !std::isfinite(relative_std)) {
      sigma = std::numeric_limits<double>::infinity();
      sigma_err = std::numeric_limits<double>::infinity();
      return;
    }

    const long double count = static_cast<long double>(calls);
    const long double sample_std = std::abs(mean) * relative_std;
    sigma = static_cast<double>(mean);
    sigma_err = static_cast<double>(sample_std / std::sqrt(count));
  }

  // Compute whether pooled sampling has a finite positive integral estimate
  bool ValidCrossSection() const {
    return CrossSectionMoments().Count() >= 2 && sigma > 0.0 &&
           sigma_err >= 0.0 && std::isfinite(sigma) && std::isfinite(sigma_err);
  }

  // Compute the relative standard error or infinity for invalid support
  double RelativeError() const {
    return ValidCrossSection() ? std::abs(sigma_err / sigma)
                               : std::numeric_limits<double>::infinity();
  }

  // Reset the stable moments for one integration batch
  void ResetIntegrationBatch() { integration_batch_moments.Reset(); }

  // Update reduced chi2 and return the current integration-batch diagnostics
  // ESS/N = (sum w)^2/[N sum w^2]
  IntegrationBatchSummary UpdateIntegrationChi2() {
    IntegrationBatchSummary summary;
    summary.samples = integration_batch_moments.Count();
    if (summary.samples == 0 || !integration_batch_moments.IsFinite()) {
      return summary;
    }

    const long double count = static_cast<long double>(summary.samples);
    const long double mean = integration_batch_moments.Mean();
    const long double mean_squared = mean * mean;
    const long double square_sum =
        integration_batch_moments.M2() + count * mean_squared;
    if (square_sum > 0.0L && std::isfinite(square_sum)) {
      const long double fraction = count * mean_squared / square_sum;
      summary.ess_fraction =
          static_cast<double>(std::clamp(fraction, 0.0L, 1.0L));
    }
    if (summary.samples < 2 || !(mean > 0.0L) ||
        !(integration_batch_moments.M2() >= 0.0L)) {
      return summary;
    }

    const long double variance =
        integration_batch_moments.M2() / (count * (count - 1.0L));
    const long double epsilon = std::numeric_limits<double>::epsilon();
    const long double variance_floor = mean_squared * epsilon * epsilon;
    summary.mean = mean;
    summary.variance_of_mean = std::max(variance, variance_floor);
    summary.valid = std::isfinite(summary.variance_of_mean) &&
                    summary.variance_of_mean > 0.0L;
    if (summary.valid) {
      chi2_batch_estimates.Add(summary.mean, summary.variance_of_mean);
      chi2 = chi2_batch_estimates.ReducedChi2();
    }
    return summary;
  }

  // Serialize the complete integration state
  void struct2json(json &j) const {
    j["amplitude_ok"] = amplitude_ok;
    j["amplitude_evaluations"] = amplitude_evaluations;
    j["amplitude_failures"] = amplitude_failures;
    j["kinematics_ok"] = kinematics_ok;
    j["fidcuts_ok"] = fidcuts_ok;
    j["vetocuts_ok"] = vetocuts_ok;
    j["technical_failures"] = technical_failures;
    j["all_ok"] = all_ok;

    j["evaluations"] = evaluations;
    j["integration_samples"] = integration_samples;
    j["integration_runtime"] = integration_runtime;
    j["trials"] = trials;

    j["generated"] = generated;
    j["N_overflow"] = N_overflow;

    j["sigma"] = sigma;
    j["sigma_err"] = sigma_err;

    j["chi2"] = chi2;

    StorePositiveMoments(j, "integrand", integrand_weight_moments);
    StorePositiveMoments(j, "proposal", proposal_density_moments);
    StorePositiveMoments(j, "importance", importance_weight_moments);
    StorePositiveMoments(j, "cross_section", continued_weight_moments);

    j["integration_batch_count"] = integration_batch_moments.Count();
    StoreExtended(j, "integration_batch_mean",
                  integration_batch_moments.Mean());
    StoreExtended(j, "integration_batch_m2", integration_batch_moments.M2());
    j["integration_batch_finite"] = integration_batch_moments.IsFinite();

    j["chi2_estimate_count"] = chi2_batch_estimates.Samples();
    StoreExtended(j, "chi2_estimate_mean", chi2_batch_estimates.Mean());
    StoreExtended(j, "chi2_estimate_scaled_m2",
                  chi2_batch_estimates.ScaledM2());
    StoreExtended(j, "chi2_estimate_scaled_weight",
                  chi2_batch_estimates.ScaledWeight());
    StoreExtended(j, "chi2_estimate_max_log_weight",
                  chi2_batch_estimates.MaxLogWeight());
  }

  // Restore the complete integration state
  void json2struct(json &j) {
    j.at("amplitude_ok").get_to(amplitude_ok);
    j.at("amplitude_evaluations").get_to(amplitude_evaluations);
    j.at("amplitude_failures").get_to(amplitude_failures);
    j.at("kinematics_ok").get_to(kinematics_ok);
    j.at("fidcuts_ok").get_to(fidcuts_ok);
    j.at("vetocuts_ok").get_to(vetocuts_ok);
    j.at("technical_failures").get_to(technical_failures);
    j.at("all_ok").get_to(all_ok);

    j.at("evaluations").get_to(evaluations);
    j.at("integration_samples").get_to(integration_samples);
    j.at("integration_runtime").get_to(integration_runtime);
    if (!(integration_runtime >= 0.0) || !std::isfinite(integration_runtime)) {
      throw std::invalid_argument(
          "Stats::json2struct: invalid integration runtime");
    }
    j.at("trials").get_to(trials);

    j.at("generated").get_to(generated);
    j.at("N_overflow").get_to(N_overflow);

    j.at("sigma").get_to(sigma);
    j.at("sigma_err").get_to(sigma_err);

    j.at("chi2").get_to(chi2);

    LoadPositiveMoments(j, "integrand", integrand_weight_moments);
    LoadPositiveMoments(j, "proposal", proposal_density_moments);
    LoadPositiveMoments(j, "importance", importance_weight_moments);
    LoadPositiveMoments(j, "cross_section", continued_weight_moments);
    ValidateFrozenSampleCounts();
    integration_batch_moments.Restore(
        j.at("integration_batch_count").get<std::size_t>(),
        LoadExtended(j, "integration_batch_mean"),
        LoadExtended(j, "integration_batch_m2"),
        j.at("integration_batch_finite").get<bool>());
    chi2_batch_estimates.Restore(
        j.at("chi2_estimate_count").get<std::size_t>(),
        LoadExtended(j, "chi2_estimate_mean"),
        LoadExtended(j, "chi2_estimate_scaled_m2"),
        LoadExtended(j, "chi2_estimate_scaled_weight"),
        LoadExtended(j, "chi2_estimate_max_log_weight"));
  }

  // Sequential physics acceptance and independent technical failures
  double amplitude_ok = 0.0;
  double amplitude_evaluations = 0.0;
  double amplitude_failures = 0.0;
  double kinematics_ok = 0.0;
  double fidcuts_ok = 0.0;
  double vetocuts_ok = 0.0;
  double technical_failures = 0.0;
  double all_ok = 0.0;

  // Keep as double to avoid overflow of range
  double evaluations = 0.0;              // Integrand evaluations
  std::uint64_t integration_samples = 0; // Frozen-proposal samples
  double integration_runtime = 0.0; // Complete integration runtime in seconds
  double trials = 0.0;              // Event generation trials

  unsigned int generated = 0.0; // Event generation
  std::uint64_t N_overflow = 0; // Weight overflows

  // Cross section and its error
  double sigma = 0.0;
  double sigma_err = 0.0;

  // Reduced chi2 of independent integration estimates
  double chi2 = 0.0;

private:
  // Store one long double as a JSON mantissa and binary exponent
  static void StoreExtended(json &document, const std::string &key,
                            long double value) {
    const statistics::ExtendedFloatParts parts =
        statistics::SplitExtendedFloat(value);
    document[key] = {{"mantissa", parts.mantissa},
                     {"exponent", parts.exponent}};
  }

  // Load one long double from a JSON mantissa and binary exponent
  static long double LoadExtended(const json &document,
                                  const std::string &key) {
    const json &encoded = document.at(key);
    return statistics::JoinExtendedFloat({encoded.at("mantissa").get<double>(),
                                          encoded.at("exponent").get<int>()});
  }

  // Store one scale-normalized positive-moment state with a common key prefix
  static void
  StorePositiveMoments(json &document, const std::string &prefix,
                       const statistics::ScaledPositiveMoments &moments) {
    document[prefix + "_count"] = moments.Count();
    document[prefix + "_positive_count"] = moments.PositiveCount();
    StoreExtended(document, prefix + "_scaled_mean", moments.ScaledMean());
    StoreExtended(document, prefix + "_scaled_m2", moments.ScaledM2());
    StoreExtended(document, prefix + "_log_scale", moments.LogScale());
    StoreExtended(document, prefix + "_min_log_value", moments.MinLogValue());
    document[prefix + "_finite"] = moments.IsFinite();
  }

  // Restore one scale-normalized positive-moment state with a common key prefix
  static void LoadPositiveMoments(const json &document,
                                  const std::string &prefix,
                                  statistics::ScaledPositiveMoments &moments) {
    moments.Restore(document.at(prefix + "_count").get<std::size_t>(),
                    document.at(prefix + "_positive_count").get<std::size_t>(),
                    LoadExtended(document, prefix + "_scaled_mean"),
                    LoadExtended(document, prefix + "_scaled_m2"),
                    LoadExtended(document, prefix + "_log_scale"),
                    LoadExtended(document, prefix + "_min_log_value"),
                    document.at(prefix + "_finite").get<bool>());
  }

  // Validate that every frozen-proposal summary describes the same samples
  void ValidateFrozenSampleCounts() const {
    const std::size_t expected = static_cast<std::size_t>(integration_samples);
    if (static_cast<std::uint64_t>(expected) != integration_samples ||
        integrand_weight_moments.Count() != expected ||
        proposal_density_moments.Count() != expected ||
        importance_weight_moments.Count() != expected ||
        (proposal_density_moments.IsFinite() &&
         proposal_density_moments.PositiveCount() != expected)) {
      throw std::invalid_argument(
          "Stats::json2struct: inconsistent frozen-proposal sample counts");
    }
    if (continued_weight_moments.Count() != 0 &&
        continued_weight_moments.Count() < expected) {
      throw std::invalid_argument(
          "Stats::json2struct: inconsistent cross-section continuation");
    }
  }

  // Compute the moments representing the current pooled cross-section estimate
  const statistics::ScaledPositiveMoments &CrossSectionMoments() const {
    return continued_weight_moments.Count() > 0 ? continued_weight_moments
                                                : importance_weight_moments;
  }

  // Start or extend the pooled cross-section estimate with a generation trial
  void ContinueCrossSection(double weight) {
    if (continued_weight_moments.Count() == 0) {
      continued_weight_moments = importance_weight_moments;
    }
    continued_weight_moments.Add(weight);
  }

  statistics::ScaledPositiveMoments integrand_weight_moments;
  statistics::ScaledPositiveMoments proposal_density_moments;
  statistics::ScaledPositiveMoments importance_weight_moments;
  statistics::ScaledPositiveMoments continued_weight_moments;
  statistics::RunningMoments integration_batch_moments;
  statistics::WeightedEstimateMoments chi2_batch_estimates;
};

// Generated weighted-event statistics for effective sample size
class WeightedEventStats {
public:
  // Reset accumulated generated-event weight sums
  void Reset() { weights.Reset(); }

  // Add one accepted generated-event weight
  void Push(double weight) { weights.Add(weight); }

  // Compute the weighted effective sample size
  // N_eff = (sum_i w_i)^2/sum_i w_i^2
  double EffectiveSampleSize() const {
    return static_cast<double>(weights.EffectiveSampleSize());
  }

  // Compute the effective sample size fraction
  double EffectiveSampleFraction() const {
    return weights.EffectiveSampleFraction();
  }

private:
  statistics::ScaledWeightSums weights;
};

class MGraniitti {
public:
  // Destructor & Constructor
  MGraniitti();
  ~MGraniitti();

  // -------------------------------------------------------------------
  // Copy, assignment and move disabled
  MGraniitti(const MGraniitti &) = delete;
  MGraniitti &operator=(const MGraniitti &) = delete;
  MGraniitti(MGraniitti &&) = delete;
  MGraniitti &operator=(MGraniitti &&) = delete;
  // -------------------------------------------------------------------

  // For cross-section for HepMC3 output
  void ForceXS(double xs) { xsforced = xs; }

  // Read parameters
  void ReadInput(const json &j);

  // Initialize memory
  void InitProcessMemory(std::string process, unsigned int RNDSEED);
  void InitMultiMemory();
  std::vector<uint32_t> GenerateUniqueSeeds(uint32_t fixed_seed,
                                            std::size_t num_seeds);

  // Set common frozen-proposal integration parameters
  void SetIntegrationParam(const INTEGRATIONPARAM &in);

  // Set common maximum-weight calibration parameters
  void SetMaxWeightParam(const MAXWEIGHTPARAM &in);

  // Set VEGAS parameters
  void SetVegasParam(const VEGASPARAM &in);

  // Set external file handle
  void SetHepMC2Output(std::shared_ptr<HepMC3::WriterAsciiHepMC2> &hepmc,
                       const std::string &OUTPUTNAME) {
    OUTPUT = OUTPUTNAME;
    FORMAT = "hepmc2";
    outputHepMC2 = hepmc;
  }
  void SetHepMC3Output(std::shared_ptr<HepMC3::WriterAscii> &hepmc,
                       const std::string &OUTPUTNAME) {
    OUTPUT = OUTPUTNAME;
    FORMAT = "hepmc3";
    outputHepMC3 = hepmc;
  }

  // Maximum weight set/get
  void SetMaxweight(double w);
  double GetMaxweight() const;

  // Get number of CPU cores
  void SetCores(int N) {
    if (N < 0) { throw std::invalid_argument("MGraniitti::SetCores: number of cores must be nonnegative"); }
    CORES = N;
    if (CORES == 0) {
      // SETUP number of threads automatically
      // x additional factor for dead-time compensation
      CORES = std::round(std::thread::hardware_concurrency() * 1.25);

      // If autodetection fails, set 1
      if (CORES < 1) {
        CORES = 1;
      }
    }
  }
  int GetCores() const { return CORES; }
  void SetIntegrator(const std::string &integrator) { INTEGRATOR = integrator; }
  void SetWeighted(bool weighted) { WEIGHTED = weighted; }
  // Output file name
  void SetOutput(const std::string &output) { OUTPUT = output; }
  // Output file format
  void SetFormat(const std::string &format) {
    if (format == "hepmc3" || format == "hepmc2" || format == "hepevt" ||
        format == "lhe") {
      FORMAT = format;
    } else {
      throw std::invalid_argument(
          "MGraniitti::SetFormat: Unknown output format: " + format +
          " (valid: hepmc3, hepmc2, hepevt, lhe)");
    }
  }

  // Special method for combining histograms from each thread
  void HistogramFusion();

  // Get cross section and error
  void GetXS(double &xs, double &xs_err) const {
    xs = stat.sigma;
    xs_err = stat.sigma_err;
  }

  // Set number of (un)weighted events to be generated
  void SetNumberOfEvents(int n) {
    if (n < -1) {
      throw std::invalid_argument(
          "MGraniitti::SetNumberOfEvents: number of events must be at least "
          "-1");
    }
    NEVENTS = n;
  }
  int GetNumberOfEvents() const { return NEVENTS; }

  // Compute whether only the adaptive integration proposal is requested
  bool IntegrationProposalOnly() const { return NEVENTS == -1; }

  // Compute processes
  std::vector<std::string> GetProcessNumbers(const std::string &tune = gra::MODELPARAM) const;
  // Compute native and generated final-state rows for the process help table
  std::vector<std::vector<std::string>> GetFinalStateRows() const;

  // Initialize generator
  void Initialize();
  void Initialize(const MEikonal &eikonal_in);

  // Process object pointer, public so methods can be accessed
  MProcess *proc = nullptr;

  void SaveVGRID() const;
  void ReadVGRID(const std::string &inputfile);

  void PrintHistograms();
  void ReadGeneralParam(const json &j);
  void ReadProcessParam(const json &j);
  // Read QED ISR and FSR switches with tune numerical controls
  void ReadRadiativeParam(const json &j);
  // Read nuclear steering after the immutable model tune is bound
  void ReadNuclearParam(const json &j);
  void ReadIntegralParam(const json &j);
  // Read integration parameters from the model numerics card
  void ReadIntegratorNumerics(const std::string &filename);
  void ReadGenCuts(const json &j);
  void ReadFidCuts(const json &j);
  void ReadVetoCuts(const json &j);
  void ReadModelParam(const std::string &tune);

  void Generate();

  std::string FULL_PROCESS_STR = "null"; // Full process string
  std::string PROCESS = "null";          // Physics process identifier

  bool WEIGHTED = false;           // Unweighted or weighted event generation
  int NEVENTS = 0;                 // Number of events to be generated
  int CORES = 0;                   // Number of CPU cores (threads) in use
  std::string INTEGRATOR = "null"; // Integrator (VEGAS, FLAT, ...)

  // SILENT OUTPUT
  bool HILJAA = false;

  void SetVgridFile(const std::string &inputfile) {
    MC_VGRID_INPUT = inputfile;
  }

  // Terminal input
  void ConstructTerminal(cxxopts::Options &options) const;
  void ProcessTerminal(json &j, cxxopts::ParseResult const &r) const;

private:
  // Model cards selected by this generator
  std::string modelparam = "TUNE0";

  const double OUTPUT_XS_SCALE = 1E12; // Output in picobarns

  // Frozen output cross section and LHE weight conversion in picobarns
  double output_xs_pb = 0.0;
  double output_xs_err_pb = 0.0;
  double output_lhe_weight_scale = 1.0;

  // MC initialized from a file
  std::string MC_VGRID_INPUT = "null";

  // Integration statistics
  Stats stat;

  // Generated weighted-event statistics
  WeightedEventStats weighted_event_stats;

  // Histogram fusion
  bool hist_fusion_done = false;

  // Synchronize state shared only by this generator instance
  mutable std::mutex state_mutex;

  // Preserve the first exception raised by this instance's worker threads
  std::exception_ptr worker_exception;

  // Preserve the first rejected amplitude diagnostic across all workers
  std::string first_amplitude_failure;

  // HepMC outputfile
  std::string FULL_OUTPUT_STR = "null";
  std::string OUTPUT = "null";
  std::string FORMAT = "null"; // hepmc3 or hepmc2 or hepevt or lhe

  std::shared_ptr<HepMC3::GenRunInfo> runinfo = nullptr;
  std::shared_ptr<HepMC3::WriterAscii> outputHepMC3 = nullptr;
  std::shared_ptr<HepMC3::WriterAsciiHepMC2> outputHepMC2 = nullptr;
  std::shared_ptr<HepMC3::WriterHEPEVT> outputHEPEVT = nullptr;
  std::shared_ptr<gra::MLHEWriter> outputLHE = nullptr;
  bool owns_file_output = false;

  // VEGAS creates copies here
  std::vector<MProcess *> pvec;
  MCentral proc_C;
  MFactorized proc_F;
  MQuasiElastic proc_Q;
  MCollinear proc_P;
  MHardDiffraction proc_D;
  neurojac::MFlowIntegrator neurojac;

  // Forced cross-section for HepMC3 output
  double xsforced = -1;

  // Global timing
  MTimer global_tictoc;
  MTimer local_tictoc;
  MTimer atime;
  double time_t0 = 0.0;

  // -----------------------------------------------
  // Process MC integration generic variables

  // Generation mode (integration = 0, event generation = 1)
  unsigned int GMODE = 0;

  // -----------------------------------------------
  // Common frozen-proposal integration and envelope parameters
  INTEGRATIONPARAM integration_param;
  MAXWEIGHTPARAM max_weight_param;
  MAXWEIGHTSTATE max_weight_state;

  // Skip the event-by-event screening loop during proposal adaptation
  bool fast_adaptation = true;

  // -----------------------------------------------
  // VEGAS MC

  // Parameters
  VEGASPARAM vparam;

  // Standalone VEGAS proposal and batch diagnostics
  MVEGASIntegrator vegas;

  // Frozen VEGAS parameters used during the first NEUROJAC training rounds
  VEGASPARAM neurojac_vegas_param;

  std::uint64_t VEGAS(VEGASStage stage, std::uint64_t calls,
                      unsigned int rounds, unsigned int N,
                      const VEGASPARAM &param);
  void VEGASMultiThread(unsigned int N, unsigned int tid, VEGASStage stage,
                        unsigned int LOCALcalls, std::uint64_t first_sample);

  // -----------------------------------------------

  void UnifyHistogramBounds();

  // Calculate cross section
  void CalculateCrossSection();

  // Vegas wrapper
  double VegasWrapper(std::vector<double> &randvec, double wgt);

  // Event sampling/generation
  void CallIntegrator(unsigned int N);
  void SampleVegas(unsigned int N);
  void SampleFlat(unsigned int N);
  void SampleNeuroJac(unsigned int N);
  // Available proposals for collecting NEUROJAC adaptation events
  enum class NeuroJacSampler { Flow, FrozenVegas };
  // Collect one thread-local NEUROJAC adaptation buffer
  void NeuroJacCollect(std::size_t thread_id, std::size_t calls,
                       NeuroJacSampler sampler,
                       std::vector<neurojac::MFlowTrainingEvent> &events);
  // Collect one multithreaded NEUROJAC adaptation batch
  std::vector<neurojac::MFlowTrainingEvent>
  NeuroJacCollectBatch(std::size_t calls, NeuroJacSampler sampler);
  // Retry fixed-size NEUROJAC bootstrap batches until training has support
  std::vector<neurojac::MFlowTrainingEvent>
  NeuroJacCollectTraining(std::size_t &evaluations, NeuroJacSampler sampler);
  // Evaluate frozen NEUROJAC samples for integration or event generation
  void NeuroJacProduction(std::size_t thread_id, std::size_t calls,
                          unsigned int requested_events, SamplingStage stage, std::uint64_t first_sample);

  // Compute a merged copy of worker spin sums without changing any worker
  spin::MQMetrics QMetrics() const;

  // Compute an automatic frozen-proposal batch size
  std::size_t AutomaticIntegrationBatch() const;
  // Compute whether the serialized state has no dedicated integration samples
  bool ProposalOnlyState() const;
  // Compute whether weighted generation samples the cross section
  bool SamplesCrossSectionDuringGeneration() const;
  // Compute the immutable maximum weight captured at generation start
  double GetGenerationMaxweight() const;
  // Initialize or update maximum-weight calibration after one sample
  void ObserveMaximumWeight(double sampling_weight, SamplingStage stage);
  // Observe one sample and report its first amplitude failure centrally
  double ObserveSample(const gra::MEventWeightState &aux, double raw_weight,
                       double log_inverse_density, SamplingStage stage);
  // Initialize the estimated maximum-weight envelope after discovery
  void InitializeMaximumWeightEnvelope();
  // Compute the required independent envelope validation count
  std::uint64_t RequiredMaximumWeightValidation() const;
  // Compute whether common frozen-proposal integration may stop
  bool IntegrationConverged() const;
  // Reject an integration run that reaches its hard sample limit
  void CheckIntegrationLimit() const;
  // Print the common integration and maximum-weight controls
  void PrintIntegrationSetup() const;
  // Print the controls which affect VEGAS proposal adaptation
  void PrintVegasAdaptationSetup() const;
  // Print the controls which affect frozen VEGAS integration
  void PrintVegasIntegrationSetup() const;
  // Print NEUROJAC diagnostics and controls using their JSON keys
  void PrintNeuroJacParameters() const;
  // Check event writers for failed output streams
  void CheckFileOutput() const;
  // Close every internally owned event writer
  void CloseFileOutput();
  // Rename one invalid partial generation output for recovery
  void PreserveFailedOutput();

  // Helper functions
  void PrintInit() const;
  int SaveEvent(MProcess *pr, double W, double MAXW,
                const gra::MEventWeightState &aux);
  // Print one NEUROJAC training update using the common status layout
  void
  PrintNeuroJacTrainingStatus(const neurojac::MFlowTrainingReport &training,
                              const neurojac::MFlowValidationReport &validation,
                              std::size_t evaluations);
  // Print one VEGAS adaptation update without cross section statistics
  void PrintVegasAdaptationStatus(const VEGASAdaptationReport &report);
  // Print frozen-proposal performance, weight and unweighting diagnostics
  void PrintIntegrationDiagnostics(double runtime) const;
  // Print sequential event flow and independent failure diagnostics
  void PrintEventFlowStatistics() const;
  void PrintStatus(unsigned int events, unsigned int N, MTimer &tictoc,
                   double timercut);
  void PrintStatistics(unsigned int N);
  gra::PARAM_RES ReadFactorized(const std::string &resparam_str);
  void InitFileOutput();

  // Interpreter commands
  std::vector<aux::OneCMD> syntax;
};

// Launcher functions
void MLaunch(MGraniitti gen, int randomseed, int tid, int events);
void MThreader(MGraniitti &gen);

} // namespace gra

#endif
