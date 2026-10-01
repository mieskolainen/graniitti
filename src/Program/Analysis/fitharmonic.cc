// Detector-aware exclusive two-pion spherical harmonic measurement
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <iterator>
#include <limits>
#include <map>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Libraries
#include "cxxopts.hpp"
#include "json.hpp"
#include "rang.hpp"

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/Program/Analysis/fitharmonic.h"
#include "Graniitti/Analysis/MHarmonicHepMC.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Analysis/MSpherical.h"

using gra::aux::indices;
using gra::program::Count;
using gra::program::WriteOutput;

namespace {

using gra::harmonic::Axis;
using gra::harmonic::Cell;
using gra::harmonic::CellEstimate;
using gra::harmonic::DataEvent;
using gra::harmonic::DetectorResponseModel;
using gra::harmonic::FiducialCuts;
using gra::harmonic::FiducialRange;
using gra::harmonic::FitConfig;
using gra::harmonic::HarmonicEstimator;
using gra::harmonic::HarmonicMeasurement;
using gra::harmonic::HepMCAnalysisReader;
using gra::harmonic::HepMCReadConfig;
using gra::harmonic::IdentityResponse;
using gra::harmonic::MeasurementMode;
using gra::harmonic::MeasurementResult;
using gra::harmonic::PhaseSpaceGrid;
using gra::harmonic::ResponseEvent;
using gra::harmonic::SampleType;
using gra::harmonic::ToyResponse;
using gra::harmonic::ToyResponseConfig;
using json = nlohmann::json;

// Read and parse one analysis card
json ReadCard(const std::string &path) {
  return json::parse(gra::aux::GetInputData(path));
}

// Parse one two-component finite interval
FiducialRange ParseRange(const json &values, const std::string &name) {
  if (!values.is_array() || values.size() != 2) {
    throw std::invalid_argument("fitharmonic:: " + name +
                                " must contain [min,max]");
  }
  const FiducialRange range{values.at(0).get<double>(),
                            values.at(1).get<double>()};
  range.Validate(name);
  return range;
}

// Parse particle-level central and forward fiducial cuts
FiducialCuts ParseFiducialCuts(const json &input, MeasurementMode mode) {
  const json &central = input.at("central");
  FiducialCuts cuts;
  cuts.pion_eta = ParseRange(central.at("pion_eta"), "pion_eta");
  cuts.pion_pt = ParseRange(central.at("pion_pt"), "pion_pt");
  cuts.proton_xi = FiducialRange{0.0, 1.0};
  cuts.proton_abs_t = FiducialRange{0.0, 1.0};
  if (mode == MeasurementMode::Tagged) {
    const json &forward = input.at("forward");
    cuts.proton_xi = ParseRange(forward.at("proton_xi"), "proton_xi");
    cuts.proton_abs_t = ParseRange(forward.at("proton_abs_t"), "proton_abs_t");
  }
  cuts.Validate(mode);
  return cuts;
}

// Parse the canonical conditional-coordinate grid
PhaseSpaceGrid ParseGrid(const json &input, MeasurementMode mode) {
  if (!input.is_array()) {
    throw std::invalid_argument("fitharmonic:: axes must be an array");
  }
  std::vector<Axis> axes;
  axes.reserve(input.size());
  for (const json &entry : input) {
    Axis axis;
    axis.coordinate = gra::harmonic::ParseCoordinate(
        entry.at("coordinate").get<std::string>());
    axis.bins = Count<std::size_t>(entry.at("bins"), "bins");
    axis.min = entry.at("min").get<double>();
    axis.max = entry.at("max").get<double>();
    axes.push_back(axis);
  }
  return PhaseSpaceGrid(mode, axes);
}

// Parse harmonic truncation and response-inversion controls
FitConfig ParseFitConfig(const json &input) {
  FitConfig config;
  config.estimator = gra::harmonic::ParseHarmonicEstimator(
      input.at("estimator").get<std::string>());
  config.lmax = Count<int>(input.at("lmax"), "lmax");
  config.remove_odd = input.value("remove_odd", false);
  config.remove_negative_m = input.value("remove_negative_m", false);
  config.svd_relative_cut = input.at("svd_relative_cut").get<double>();
  config.response_jackknife_bins =
      Count<std::size_t>(input.value("response_jackknife_bins", json(16)), "response_jackknife_bins");
  config.max_parameters = Count<std::size_t>(input.value("max_parameters", json(1500)), "max_parameters");
  if (config.estimator == HarmonicEstimator::EML) {
    const json &eml = input.at("eml");
    config.eml_max_calls = Count<std::size_t>(eml.value("max_calls", json(100000)), "max_calls");
    config.eml_positivity_costheta =
        Count<std::size_t>(eml.value("positivity_costheta", json(24)), "positivity_costheta");
    config.eml_positivity_phi = Count<std::size_t>(eml.value("positivity_phi", json(48)), "positivity_phi");
    config.eml_min_intensity = eml.value("minimum_intensity", 1e-10);
  }
  config.Validate();
  return config;
}

// Parse the common event-projection configuration
HepMCReadConfig ParseReadConfig(const json &card, MeasurementMode mode,
                                const FiducialCuts &cuts) {
  HepMCReadConfig config;
  config.mode = mode;
  config.frame = gra::harmonic::ParseAngularFrame(card.value("frame", "CS"));
  config.cuts = cuts;
  if (mode == MeasurementMode::Tagged) {
    const json &selection = card.at("detector_selection");
    config.selection.pt_balance_max =
        selection.at("pt_balance_max").get<double>();
    config.selection.mass_match_relative_max =
        selection.at("mass_match_relative_max").get<double>();
    config.selection.rapidity_match_max =
        selection.at("rapidity_match_max").get<double>();
  }
  config.sqrt_s = card.at("sqrt_s").get<double>();
  config.max_records =
      Count<std::uint64_t>(card.value("maximum", json(std::numeric_limits<std::uint64_t>::max())), "maximum");
  config.seed = Count<std::uint64_t>(card.value("seed", json(1)), "seed");
  config.weight_index = Count<std::size_t>(card.value("event_weight_index", json(0)), "event_weight_index");
  config.Validate();
  return config;
}

// Parse the closure-only analytic detector response
ToyResponseConfig ParseToyResponse(const json &input, std::uint64_t seed) {
  ToyResponseConfig config;
  config.pion_efficiency = input.at("pion_efficiency").get<double>();
  config.proton_efficiency = input.at("proton_efficiency").get<double>();
  config.pion_logpt_sigma = input.at("pion_logpt_sigma").get<double>();
  config.pion_eta_sigma = input.at("pion_eta_sigma").get<double>();
  config.pion_phi_sigma = input.at("pion_phi_sigma").get<double>();
  config.proton_pt_sigma = input.at("proton_pt_sigma").get<double>();
  config.proton_logpz_sigma = input.at("proton_logpz_sigma").get<double>();
  config.seed = Count<std::uint64_t>(input.value("seed", json(seed)), "seed");
  config.Validate();
  return config;
}

// Apply a positive sampling or importance-weight scale
void ScaleResponse(std::vector<ResponseEvent> &events, double scale) {
  if (!std::isfinite(scale) || scale <= 0.0) {
    throw std::invalid_argument(
        "fitharmonic:: response sample scale must be positive");
  }
  for (ResponseEvent &event : events) {
    event.truth.weight *= scale;
    if (event.reco.has_value()) {
      event.reco->weight *= scale;
    }
  }
}

// Append one event vector without retaining an intermediate copy
template <typename T>
void Append(std::vector<T> &target, std::vector<T> source) {
  target.reserve(target.size() + source.size());
  std::move(source.begin(), source.end(), std::back_inserter(target));
}

// Read one named detector-response group from one or more MC samples
std::pair<std::string, std::vector<ResponseEvent>>
ReadResponseGroup(const json &definition, const HepMCAnalysisReader &reader,
                  std::uint64_t seed) {
  const std::string name = definition.at("name").get<std::string>();
  const std::string mode = definition.at("mode").get<std::string>();
  const json &samples = definition.at("samples");
  if (name.empty() || !samples.is_array() || samples.empty()) {
    throw std::invalid_argument(
        "fitharmonic:: each response needs a name and samples");
  }

  std::unique_ptr<DetectorResponseModel> model;
  if (mode == "IDENTITY") {
    model = std::make_unique<IdentityResponse>();
  } else if (mode == "TOY") {
    model = std::make_unique<ToyResponse>(
        ParseToyResponse(definition.at("toy"), seed));
  } else if (mode != "PAIRED") {
    throw std::invalid_argument(
        "fitharmonic:: response mode must be IDENTITY, TOY or PAIRED");
  }

  std::vector<ResponseEvent> events;
  for (const json &sample : samples) {
    if (sample.at("reference").get<std::string>() != "ANGULAR_FLAT") {
      throw std::invalid_argument(
          "fitharmonic:: response samples must use ANGULAR_FLAT reference");
    }
    const std::string truth = sample.at("truth").get<std::string>();
    std::vector<ResponseEvent> current;
    if (mode == "PAIRED") {
      current = reader.ReadPairedResponse(truth,
                                          sample.at("reco").get<std::string>());
    } else {
      current = reader.ReadModelResponse(truth, *model);
    }
    ScaleResponse(current, sample.at("scale").get<double>());
    Append(events, std::move(current));
  }
  return {name, std::move(events)};
}

// Read all selected detector-level files belonging to one data sample
std::vector<DataEvent> ReadDataSample(const json &definition,
                                      const HepMCAnalysisReader &reader) {
  const json &paths = definition.at("paths");
  if (!paths.is_array() || paths.empty()) {
    throw std::invalid_argument(
        "fitharmonic:: each measurement sample needs paths");
  }
  std::vector<DataEvent> events;
  for (const json &path : paths) {
    Append(events, reader.ReadData(path.get<std::string>()));
  }
  return events;
}

// Convert one project matrix to a rectangular JSON array
json MatrixJSON(const gra::MMatrix<double> &matrix) {
  json output = json::array();
  for (std::size_t row = 0; row < matrix.size_row(); ++row) {
    json values = json::array();
    for (std::size_t col = 0; col < matrix.size_col(); ++col) {
      values.push_back(matrix[row][col]);
    }
    output.push_back(std::move(values));
  }
  return output;
}

// Compute physical bounds for one regular-grid cell
json CellBounds(const Cell &cell, const PhaseSpaceGrid &grid) {
  json output = json::array();
  for (const auto &axis : indices(cell)) {
    const Axis &definition = grid.Axes()[axis];
    const double width = (definition.max - definition.min) /
                         static_cast<double>(definition.bins);
    const double lower =
        definition.min + width * static_cast<double>(cell[axis]);
    const double upper =
        cell[axis] + 1 == definition.bins ? definition.max :
        definition.min + width * static_cast<double>(cell[axis] + 1);
    output.push_back(
        {{"coordinate", gra::harmonic::ToString(definition.coordinate)},
         {"min", lower},
         {"max", upper}});
  }
  return output;
}

// Convert one full-basis cell estimate to JSON
json EstimateJSON(const CellEstimate &estimate) {
  return {{"valid", estimate.valid},
          {"coefficients", estimate.value},
          {"covariance", MatrixJSON(estimate.covariance)},
          {"entries", estimate.entries},
          {"sum_weight", estimate.sum_weight},
          {"sum_weight2", estimate.sum_weight2}};
}

// Convert one phase-space level to ordered JSON cells
json LevelJSON(const std::vector<Cell> &cells,
               const std::map<Cell, CellEstimate> &estimates,
               const PhaseSpaceGrid &grid) {
  json output = json::array();
  for (const Cell &cell : cells) {
    const auto estimate = estimates.find(cell);
    if (estimate == estimates.end()) {
      continue;
    }
    output.push_back({{"index", cell},
                      {"bounds", CellBounds(cell, grid)},
                      {"estimate", EstimateJSON(estimate->second)}});
  }
  return output;
}

// Compute the linear-index convention for all stored coefficients
json CoefficientConvention(int lmax, const std::vector<std::size_t> &active) {
  json output = json::array();
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      const std::size_t index =
          static_cast<std::size_t>(gra::spherical::LinearInd(l, m));
      output.push_back({{"index", index},
                        {"l", l},
                        {"m", m},
                        {"production_active", std::find(active.begin(), active.end(),
                                             index) != active.end()}});
    }
  }
  return output;
}

// Compute the cell-major ordering of one global covariance matrix
json GlobalCovarianceOrder(const std::vector<Cell> &cells,
                           const std::vector<std::size_t> &active) {
  json output = json::array();
  for (const Cell &cell : cells) {
    for (const std::size_t coefficient : active) {
      output.push_back({{"cell", cell}, {"coefficient_index", coefficient}});
    }
  }
  return output;
}

// Convert a complete three-level measurement to JSON
json ResultJSON(const MeasurementResult &result, const PhaseSpaceGrid &grid,
                const FitConfig &config) {
  json estimator_diagnostics = json::object();
  if (result.estimator == HarmonicEstimator::EML) {
    estimator_diagnostics = {
        {"objective", result.objective},
        {"minimum_flat_intensity", result.minimum_flat_intensity},
        {"minimum_detector_intensity", result.minimum_detector_intensity}};
  }
  return {
      {"estimator", gra::harmonic::ToString(result.estimator)},
      {"active_coefficients", result.active_coefficients},
      {"coefficient_convention",
       CoefficientConvention(config.lmax, result.active_indices)},
      {"response_events", result.response_events},
      {"response_rows", result.response_rows},
      {"response_columns", result.response_columns},
      {"response_rank", result.response_rank},
      {"retained_singular_values", result.retained_singular_values},
      {"response_condition_number", result.response_condition_number},
      {"response_jackknife_replicas", result.response_jackknife_replicas},
      {"data_events", result.data_events},
      {"data_sum_weight", result.data_sum_weight},
      {"data_sum_weight2", result.data_sum_weight2},
      {"estimator_diagnostics", estimator_diagnostics},
      {"angular_flat",
       {{"cells", LevelJSON(result.truth_cells, result.flat, grid)},
        {"global_covariance_order",
         GlobalCovarianceOrder(result.truth_cells, result.active_indices)},
        {"global_covariance", MatrixJSON(result.flat_covariance)}}},
      {"fiducial",
       {{"cells", LevelJSON(result.truth_cells, result.fiducial, grid)},
        {"global_covariance_order",
         GlobalCovarianceOrder(result.truth_cells, result.moment_indices)},
        {"global_covariance", MatrixJSON(result.fiducial_covariance)}}},
      {"detector",
       {{"cells", LevelJSON(result.detector_cells, result.detector, grid)},
        {"global_covariance_order",
         GlobalCovarianceOrder(result.detector_cells, result.moment_indices)},
        {"global_covariance", MatrixJSON(result.detector_covariance)}}}};
}

// Normalize one result by integrated luminosity when requested
std::string ApplySampleNormalization(const json &definition,
                                     MeasurementResult &result) {
  const std::string normalization =
      definition.at("normalization").get<std::string>();
  if (normalization == "EVENT_YIELD") {
    return "events";
  }
  if (normalization != "CROSS_SECTION_PB") {
    throw std::invalid_argument(
        "fitharmonic:: normalization must be EVENT_YIELD or CROSS_SECTION_PB");
  }
  const double luminosity =
      definition.at("integrated_luminosity_pb").get<double>();
  const double relative_uncertainty =
      definition.value("luminosity_relative_uncertainty", 0.0);
  gra::harmonic::ApplyNormalization(result, luminosity, relative_uncertainty);
  return "pb";
}

} // namespace

// Run one card-defined detector-aware harmonic measurement
int main(int argc, char *argv[]) {
  gra::aux::PrintArgv(argc, argv);
  gra::aux::PrintFlashScreen(rang::fg::blue);
  std::cout << rang::style::bold
            << "GRANIITTI - Detector-aware spherical harmonic measurement"
            << rang::style::reset << std::endl
            << std::endl;
  gra::aux::PrintVersion();

  try {
    cxxopts::Options options(argv[0], "");
    options.add_options()("c,card", "Analysis card <path>",
                          cxxopts::value<std::string>())(
        "X,maximum", "Maximum event records per input <value>",
        cxxopts::value<std::uint64_t>())("H,help", "Help");
    const auto parsed = options.parse(argc, argv);
    if (parsed.count("help") || !parsed.count("card")) {
      std::cout << options.help({""}) << std::endl;
      return parsed.count("help") ? EXIT_SUCCESS : EXIT_FAILURE;
    }

    json card = ReadCard(parsed["card"].as<std::string>());
    if (parsed.count("maximum")) {
      const std::uint64_t maximum = parsed["maximum"].as<std::uint64_t>();
      if (maximum == 0) {
        throw std::invalid_argument("fitharmonic:: maximum must be positive");
      }
      card["maximum"] = maximum;
    }
    if (card.at("schema").get<std::string>() !=
        "GRANIITTI_HARMONIC_ANALYSIS_V1") {
      throw std::invalid_argument(
          "fitharmonic:: unsupported analysis-card schema");
    }
    const MeasurementMode mode = gra::harmonic::ParseMeasurementMode(
        card.at("measurement").get<std::string>());
    const FiducialCuts cuts = ParseFiducialCuts(card.at("fiducial"), mode);
    const PhaseSpaceGrid grid = ParseGrid(card.at("axes"), mode);
    const FitConfig fit_config = ParseFitConfig(card.at("fit"));
    const HepMCReadConfig read_config = ParseReadConfig(card, mode, cuts);
    const HepMCAnalysisReader reader(read_config);

    std::map<std::string, std::vector<ResponseEvent>> responses;
    const json &response_definitions = card.at("responses");
    if (!response_definitions.is_array() || response_definitions.empty()) {
      throw std::invalid_argument(
          "fitharmonic:: responses must be a non-empty array");
    }
    for (const json &definition : response_definitions) {
      auto [name, events] = ReadResponseGroup(
          definition, reader, read_config.seed);
      if (!responses.emplace(name, std::move(events)).second) {
        throw std::invalid_argument(
            "fitharmonic:: response names must be unique");
      }
    }

    json output;
    output["schema"] = "GRANIITTI_HARMONIC_MEASUREMENT_V1";
    output["configuration"] = card;
    output["measurement"] = gra::harmonic::ToString(mode);
    output["frame"] = gra::harmonic::ToString(read_config.frame);
    output["angular_daughter"] = "pi+";
    output["coordinate_convention"] = {
        {"ABST1", "positive z proton arm"},
        {"ABST2", "negative z proton arm"},
        {"DPHI_PP", "wrap(phi_positive_z - phi_negative_z) in (-pi,pi]"}};
    output["samples"] = json::array();

    const json &samples = card.at("samples");
    if (!samples.is_array() || samples.empty()) {
      throw std::invalid_argument(
          "fitharmonic:: samples must be a non-empty array");
    }
    std::set<std::string> sample_names;
    for (const json &definition : samples) {
      const std::string name = definition.at("name").get<std::string>();
      if (name.empty() || !sample_names.insert(name).second) {
        throw std::invalid_argument(
            "fitharmonic:: sample names must be non-empty and unique");
      }
      const SampleType type = gra::harmonic::ParseSampleType(
          definition.at("type").get<std::string>());
      if (type == SampleType::ResponseMC) {
        throw std::invalid_argument(
            "fitharmonic:: RESPONSE_MC belongs under responses");
      }
      const std::string response_name =
          definition.at("response").get<std::string>();
      const auto response = responses.find(response_name);
      if (response == responses.end()) {
        throw std::invalid_argument(
            "fitharmonic:: sample references an unknown response");
      }

      const std::vector<DataEvent> data = ReadDataSample(definition, reader);
      const HarmonicMeasurement measurement(grid, fit_config);
      MeasurementResult result = measurement.Fit(response->second, data);
      const std::string unit = ApplySampleNormalization(definition, result);
      output["samples"].push_back(
          {{"name", name},
           {"type", gra::harmonic::ToString(type)},
           {"response", response_name},
           {"coefficient_unit", unit},
           {"bin_normalization", "BIN_INTEGRAL"},
           {"result", ResultJSON(result, grid, fit_config)}});
    }

    WriteOutput(card.at("output").get<std::string>(), output);
    std::cout << "Wrote " << card.at("output").get<std::string>() << std::endl;
  } catch (const cxxopts::OptionException &error) {
    std::cerr << "fitharmonic:: command-line error: " << error.what()
              << std::endl;
    return EXIT_FAILURE;
  } catch (const json::exception &error) {
    std::cerr << "fitharmonic:: JSON error: " << error.what() << std::endl;
    return EXIT_FAILURE;
  } catch (const std::exception &error) {
    std::cerr << "fitharmonic:: error: " << error.what() << std::endl;
    return EXIT_FAILURE;
  }

  std::cout << "[fitharmonic: done]" << std::endl;
  return EXIT_SUCCESS;
}
