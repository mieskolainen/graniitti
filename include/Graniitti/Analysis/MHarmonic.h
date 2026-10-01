// Detector aware spherical harmonics measurement model
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef ANALYSIS_MHARMONIC_H
#define ANALYSIS_MHARMONIC_H

// C++
#include <cstddef>
#include <cstdint>
#include <map>
#include <optional>
#include <string>
#include <vector>

// Own
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Kinematics/M4Vec.h"

namespace gra {
namespace harmonic {

// Define which measured degrees of freedom condition the angular expansion
enum class MeasurementMode { Central, Tagged };

// Identify the role of each input sample
enum class SampleType { Data, ResponseMC, ClosureMC };

// Identify the three explicitly separated analysis spaces
enum class PhaseSpaceLevel { Flat, Fiducial, Detector };

// Identify the statistical estimator applied to the shared detector response
enum class HarmonicEstimator { Algebraic, EML };

// Identify one coordinate of the conditional kinematic vector
enum class Coordinate { Mass, Momentum, Rapidity, AbsT1, AbsT2, DeltaPhiPP };

// Compute the canonical coordinate ordering for one measurement mode
std::vector<Coordinate> Coordinates(MeasurementMode mode);

// Compute a stable printable measurement-mode name
std::string ToString(MeasurementMode mode);

// Parse a measurement mode without accepting implicit aliases
MeasurementMode ParseMeasurementMode(const std::string &value);

// Compute a stable printable sample-type name
std::string ToString(SampleType type);

// Parse an input sample type without accepting implicit aliases
SampleType ParseSampleType(const std::string &value);

// Compute a stable printable phase-space-level name
std::string ToString(PhaseSpaceLevel level);

// Compute a stable printable estimator name
std::string ToString(HarmonicEstimator estimator);

// Parse an estimator without accepting implicit aliases
HarmonicEstimator ParseHarmonicEstimator(const std::string &value);

// Compute a stable printable coordinate name
std::string ToString(Coordinate coordinate);

// Parse a conditional coordinate without accepting implicit aliases
Coordinate ParseCoordinate(const std::string &value);

// Describe one regular conditional-coordinate axis
struct Axis {
  Coordinate coordinate = Coordinate::Mass;
  std::size_t bins = 0;
  double min = 0.0;
  double max = 0.0;

  // Validate the finite non-empty axis definition
  void Validate() const;
};

using Cell = std::vector<std::size_t>;

// Hold one angular observation and its conditional kinematics
struct Observation {
  double costheta = 0.0;
  double phi = 0.0;
  std::vector<double> z;
  double weight = 1.0;

  // Validate all angular, kinematic and weight components
  void Validate(std::size_t dimensions, const std::string &context) const;
};

// Keep central and forward particle-level fiducial decisions separate
struct FiducialDecision {
  bool central = false;
  bool forward = false;

  // Apply only the forward requirement present in the tagged mode
  bool Pass(MeasurementMode mode) const;
};

// Pair one generated event with its optional detector-level observation
struct ResponseEvent {
  std::uint64_t event_key = 0;
  std::string source;
  Observation truth;
  FiducialDecision fiducial;
  std::optional<Observation> reco;
  bool reco_selected = false;

  // Test reconstructed detector selection independently of truth fiducial cuts
  bool DetectorAccepted() const;
};

// Hold one selected detector-level event from data or closure MC
struct DataEvent {
  std::uint64_t event_key = 0;
  std::string source;
  Observation reco;
};

// Hold the particle momenta needed for physical detector reconstruction
struct EventKinematics {
  M4Vec pip;
  M4Vec pim;
  M4Vec beam_plus;
  M4Vec beam_minus;
  M4Vec proton_plus;
  M4Vec proton_minus;
  MeasurementMode mode = MeasurementMode::Central;
  bool has_forward_protons = false;
  double weight = 1.0;
};

// Describe a thread-safe detector response attachment point
class DetectorResponseModel {
public:
  // Destroy one response provider through its stable interface
  virtual ~DetectorResponseModel() = default;

  // Compute reconstruction efficiency before cuts on reconstructed particles
  virtual double Efficiency(const EventKinematics &truth) const = 0;

  // Sample reconstructed particle momenta before applying detector cuts
  virtual EventKinematics Reconstruct(const EventKinematics &truth,
                                      std::uint64_t event_key) const = 0;
};

// Hold a conditional angular coefficient prediction and its uncertainty
struct CoefficientPrediction {
  std::vector<double> coefficients;
  MMatrix<double> covariance;
};

// ** A neural network ready continuous coefficient-field attachment point **
class ConditionalCoefficientField {
public:
  // Destroy one coefficient provider through its stable interface
  virtual ~ConditionalCoefficientField() = default;

  // Predict angular coefficients and covariance at one kinematic point
  virtual CoefficientPrediction
  Evaluate(const std::vector<double> &z) const = 0;
};

// Configure the explicit closure only analytic detector model
struct ToyResponseConfig {
  double pion_efficiency = 1.0;
  double proton_efficiency = 1.0;
  // Central tracking resolutions in log(pT), pseudorapidity and azimuth
  double pion_logpt_sigma = 0.0;
  double pion_eta_sigma = 0.0;
  double pion_phi_sigma = 0.0;
  // Forward transverse resolution in GeV and longitudinal log(|pz|) resolution
  double proton_pt_sigma = 0.0;
  double proton_logpz_sigma = 0.0;
  std::uint64_t seed = 1;

  // Validate particle efficiencies and momentum resolutions
  void Validate() const;
};

// Provide an identity detector model for algebra and closure tests
class IdentityResponse final : public DetectorResponseModel {
public:
  // Compute unit efficiency
  double Efficiency(const EventKinematics &truth) const override;

  // Compute the unchanged generated particle momenta
  EventKinematics Reconstruct(const EventKinematics &truth,
                              std::uint64_t event_key) const override;
};

// Provide a deterministic analytic detector model for closure tests
class ToyResponse final : public DetectorResponseModel {
public:
  // Construct a validated analytic response
  explicit ToyResponse(const ToyResponseConfig &config);

  // Compute the efficiency to reconstruct all required charged particles
  double Efficiency(const EventKinematics &truth) const override;

  // Smear measured momenta while preserving particle masses
  EventKinematics Reconstruct(const EventKinematics &truth,
                              std::uint64_t event_key) const override;

private:
  ToyResponseConfig config_;
};

// Apply detector efficiency and response without shared random state
std::optional<EventKinematics> SimulateDetector(const DetectorResponseModel &model,
                                                const EventKinematics &truth,
                                                std::uint64_t event_key,
                                                std::uint64_t seed);

// Map conditional coordinates to occupied regular hypercells
class PhaseSpaceGrid {
public:
  // Construct and validate a mode-specific conditional grid
  PhaseSpaceGrid(MeasurementMode mode, const std::vector<Axis> &axes);

  // Locate one point and reject coordinates outside the analysis range
  std::optional<Cell> Locate(const std::vector<double> &z) const;

  // Compute the measurement mode
  MeasurementMode Mode() const;

  // Compute the ordered axes
  const std::vector<Axis> &Axes() const;

private:
  MeasurementMode mode_;
  std::vector<Axis> axes_;
};

// Configure the sparse-cell global angular response fit
struct FitConfig {
  HarmonicEstimator estimator = HarmonicEstimator::Algebraic;
  int lmax = 2;
  // Restrict production coefficients while retaining all measured moments
  bool remove_odd = false;
  bool remove_negative_m = false;
  double svd_relative_cut = 1e-3;
  std::size_t response_jackknife_bins = 16;
  std::size_t max_parameters = 1500;
  std::size_t eml_max_calls = 100000;
  std::size_t eml_positivity_costheta = 24;
  std::size_t eml_positivity_phi = 48;
  double eml_min_intensity = 1e-10;

  // Validate truncation, regularization and memory bounds
  void Validate() const;
};

// Hold one full-basis moment estimate in one conditional cell
struct CellEstimate {
  std::vector<double> value;
  MMatrix<double> covariance;
  double sum_weight = 0.0;
  double sum_weight2 = 0.0;
  std::size_t entries = 0;
  bool valid = false;
};

// Hold the three-level result and the full inter-cell covariances
struct MeasurementResult {
  MeasurementMode mode = MeasurementMode::Central;
  HarmonicEstimator estimator = HarmonicEstimator::Algebraic;
  std::vector<Cell> truth_cells;
  std::vector<Cell> detector_cells;
  std::map<Cell, CellEstimate> flat;
  std::map<Cell, CellEstimate> fiducial;
  std::map<Cell, CellEstimate> detector;
  MMatrix<double> flat_covariance;
  MMatrix<double> fiducial_covariance;
  MMatrix<double> detector_covariance;
  // Order flat-space covariances by the fitted production coefficients
  std::vector<std::size_t> active_indices;
  // Order fiducial and detector covariances by all retained angular moments
  std::vector<std::size_t> moment_indices;
  std::size_t active_coefficients = 0;
  std::size_t response_rows = 0;
  std::size_t response_columns = 0;
  std::size_t response_rank = 0;
  std::size_t retained_singular_values = 0;
  double response_condition_number = 0.0;
  std::size_t response_events = 0;
  std::size_t response_jackknife_replicas = 0;
  std::size_t data_events = 0;
  double data_sum_weight = 0.0;
  double data_sum_weight2 = 0.0;
  double objective = 0.0;
  double minimum_flat_intensity = 0.0;
  double minimum_detector_intensity = 0.0;
};

// Convert event yields to a normalized measurement with correlated uncertainty
void ApplyNormalization(MeasurementResult &result, double divisor,
                        double relative_uncertainty);

// Fit a migration-aware response across occupied conditional hypercells
class HarmonicMeasurement {
public:
  // Construct an immutable measurement definition
  HarmonicMeasurement(const PhaseSpaceGrid &grid, const FitConfig &config);

  // Invert detector moments to angular-flat space and project to fiducial space
  MeasurementResult Fit(const std::vector<ResponseEvent> &response,
                        const std::vector<DataEvent> &data) const;

private:
  PhaseSpaceGrid grid_;
  FitConfig config_;
};

} // namespace harmonic
} // namespace gra

#endif
