// Detector aware spherical harmonics measurement model
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <memory>
#include <numeric>
#include <random>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Eigen
#include <Eigen/Dense>

// ROOT
#include "Math/Factory.h"
#include "Math/Functor.h"
#include "Math/Minimizer.h"

// Own
#include "Graniitti/Analysis/MHarmonic.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Analysis/MSpherical.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra {
namespace harmonic {

namespace {

using BlockKey = std::pair<Cell, Cell>;

struct ResponseSums {
  std::map<Cell, double> generated;
  std::map<Cell, Eigen::MatrixXd> fiducial;
  std::map<BlockKey, Eigen::MatrixXd> detector;
};

struct LinearResponse {
  Eigen::MatrixXd detector;
  Eigen::MatrixXd fiducial;
};

struct DataSums {
  Eigen::VectorXd values;
  Eigen::MatrixXd covariance;
  std::map<Cell, CellEstimate> cells;
  std::size_t accepted = 0;
  double sum_weight = 0.0;
  double sum_weight2 = 0.0;
  double effective_entries = 0.0;
};

struct EMLCellData {
  std::size_t detector_cell = 0;
  Eigen::MatrixXd basis;
  Eigen::VectorXd weight;
};

struct EMLModel {
  Eigen::MatrixXd response;
  Eigen::MatrixXd angular_basis;
  Eigen::MatrixXd detector_basis;
  std::vector<EMLCellData> data;
  std::size_t nactive = 0;
  std::size_t nmoments = 0;
  std::size_t zero_index = 0;
};

struct EstimatorResult {
  Eigen::VectorXd flat;
  Eigen::MatrixXd flat_covariance;
  double objective = 0.0;
  double minimum_flat_intensity = 0.0;
  double minimum_detector_intensity = 0.0;
};

struct HarmonicPseudoInverse {
  Eigen::MatrixXd value;
  std::size_t numerical_rank = 0;
  std::size_t retained_rank = 0;
  double condition_number = 0.0;
};

// Build one harmonic fit pseudoinverse through the shared matrix algebra
HarmonicPseudoInverse PseudoInverse(const Eigen::MatrixXd &matrix,
                                    const double relative_cut) {
  PseudoInverseDiagnostics diagnostics;
  const MMatrix<double> inverse =
      MMatrix<double>::FromEigen(matrix).PseudoInverse(relative_cut,
                                                       &diagnostics);
  return {inverse.ToEigen(), diagnostics.numerical_rank,
          diagnostics.retained_rank, diagnostics.condition_number};
}

// Mix an event key into a deterministic random seed
std::uint64_t MixSeed(std::uint64_t seed, std::uint64_t event_key) {
  std::uint64_t value = seed + 0x9e3779b97f4a7c15ULL + event_key;
  value = (value ^ (value >> 30U)) * 0xbf58476d1ce4e5b9ULL;
  value = (value ^ (value >> 27U)) * 0x94d049bb133111ebULL;
  return value ^ (value >> 31U);
}

// Compute active full-basis coefficient indices
std::vector<std::size_t> ActiveIndices(const FitConfig &config) {
  std::vector<std::size_t> indices;
  for (int l = 0; l <= config.lmax; ++l) {
    if (config.remove_odd && l % 2 != 0) {
      continue;
    }
    for (int m = -l; m <= l; ++m) {
      if (config.remove_negative_m && m < 0) {
        continue;
      }
      indices.push_back(static_cast<std::size_t>(spherical::LinearInd(l, m)));
    }
  }
  return indices;
}

// Evaluate the normalized real spherical basis used by the response integral
// B_lm = sqrt(4pi) Y_lm^R(cos(theta),phi)
Eigen::VectorXd Basis(const Observation &event,
                      const std::vector<std::size_t> &active, int lmax) {
  const std::size_t ncoef = static_cast<std::size_t>((lmax + 1) * (lmax + 1));
  std::vector<double> full(ncoef, 0.0);
  const double normalization = std::sqrt(4.0 * gra::math::PI);
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      const std::size_t index =
          static_cast<std::size_t>(spherical::LinearInd(l, m));
      full[index] = normalization * gra::math::Y_real_basis(event.costheta, event.phi, l, m);
    }
  }
  Eigen::VectorXd output(static_cast<Eigen::Index>(active.size()));
  for (const auto &index : indices(active)) {
    output[static_cast<Eigen::Index>(index)] = full[active[index]];
  }
  return output;
}

// Add one weighted outer product to a sparse response block
// R += w b_left b_right^T
void AddOuter(std::map<Cell, Eigen::MatrixXd> &blocks, const Cell &key,
              const Eigen::VectorXd &left, const Eigen::VectorXd &right,
              double weight) {
  auto iterator =
      blocks.try_emplace(key, Eigen::MatrixXd::Zero(left.size(), right.size()))
          .first;
  iterator->second.noalias() += weight * left * right.transpose();
}

// Add one weighted outer product to a migration response block
// R += w b_left b_right^T
void AddOuter(std::map<BlockKey, Eigen::MatrixXd> &blocks, const BlockKey &key,
              const Eigen::VectorXd &left, const Eigen::VectorXd &right,
              double weight) {
  auto iterator =
      blocks.try_emplace(key, Eigen::MatrixXd::Zero(left.size(), right.size()))
          .first;
  iterator->second.noalias() += weight * left * right.transpose();
}

// Validate the response observations and compute a separate weight scale per truth cell
std::map<Cell, double> ResponseScales(const std::vector<ResponseEvent> &response,
                                      const PhaseSpaceGrid &grid) {
  std::map<Cell, statistics::ScaledWeightSums> weights;
  for (const auto &event : response) {
    event.truth.Validate(grid.Axes().size(), "harmonic::HarmonicMeasurement response truth");
    if (std::fpclassify(event.truth.weight) == FP_ZERO) { continue; }
    const auto cell = grid.Locate(event.truth.z);
    if (event.DetectorAccepted()) {
      event.reco->Validate(grid.Axes().size(), "harmonic::HarmonicMeasurement response reco");
      if (!cell && grid.Locate(event.reco->z)) {
        throw std::invalid_argument(
            "harmonic::HarmonicMeasurement has detector feed-in from outside "
            "the truth grid, add truth guard bins");
      }
    }
    if (cell) { weights[*cell].Add(event.truth.weight); }
  }
  std::map<Cell, double> scales;
  for (const auto &[cell, sum] : weights) {
    if (!(sum.ScaledSum() > 0.0L) || !sum.HasSignificantSignedSum(std::numeric_limits<double>::epsilon())) {
      throw std::invalid_argument("harmonic::HarmonicMeasurement response truth cell has non-positive or cancelling weights");
    }
    scales[cell] = static_cast<double>(sum.Scale());
  }
  return scales;
}

// Fold production coefficients into all fiducial and detector moments
// Acceptance can break production symmetries even within the same l truncation
// [REFERENCE: https://arxiv.org/abs/1503.04100, Sec. III B]
void AddResponseEvent(ResponseSums &sums, const ResponseEvent &event,
                      const PhaseSpaceGrid &grid,
                      const std::vector<std::size_t> &active,
                      const std::vector<std::size_t> &moments, int lmax,
                      const std::map<Cell, double> &scales) {
  if (std::fpclassify(event.truth.weight) == FP_ZERO) { return; }
  const std::optional<Cell> truth_cell = grid.Locate(event.truth.z);
  if (!truth_cell) { return; }
  const double weight = event.truth.weight / scales.at(*truth_cell);
  sums.generated[*truth_cell] += weight;
  const Eigen::VectorXd truth_basis = Basis(event.truth, active, lmax);

  if (event.fiducial.Pass(grid.Mode())) {
    const Eigen::VectorXd fiducial_basis = Basis(event.truth, moments, lmax);
    AddOuter(sums.fiducial, *truth_cell, fiducial_basis, truth_basis, weight);
  }
  if (!event.DetectorAccepted()) {
    return;
  }
  const std::optional<Cell> reco_cell = grid.Locate(event.reco->z);
  if (!reco_cell.has_value()) {
    return;
  }
  const Eigen::VectorXd reco_basis = Basis(*event.reco, moments, lmax);
  AddOuter(sums.detector, {*reco_cell, *truth_cell}, reco_basis, truth_basis,
           weight);
}

// Compute a stable ordered set of cells from map keys
template <typename T>
std::vector<Cell> MapCells(const std::map<Cell, T> &input) {
  std::vector<Cell> cells;
  cells.reserve(input.size());
  for (const auto &[cell, value] : input) {
    static_cast<void>(value);
    cells.push_back(cell);
  }
  return cells;
}

// Compute a cell to contiguous-block lookup
std::map<Cell, std::size_t> CellOffsets(const std::vector<Cell> &cells) {
  std::map<Cell, std::size_t> offsets;
  for (const auto &index : indices(cells)) {
    offsets.emplace(cells[index], index);
  }
  return offsets;
}

// Assemble normalized global detector and truth-fiducial response matrices
// R_ab = sum_i w_i B_ia^reco B_ib^truth / sum_i w_i by truth cell
LinearResponse AssembleResponse(const ResponseSums &sums,
                                const std::vector<Cell> &truth_cells,
                                const std::vector<Cell> &detector_cells,
                                std::size_t nactive, std::size_t nmoments) {
  const std::map<Cell, std::size_t> truth_offset = CellOffsets(truth_cells);
  const std::map<Cell, std::size_t> detector_offset =
      CellOffsets(detector_cells);
  LinearResponse response{
      Eigen::MatrixXd::Zero(
          static_cast<Eigen::Index>(detector_cells.size() * nmoments),
          static_cast<Eigen::Index>(truth_cells.size() * nactive)),
      Eigen::MatrixXd::Zero(
          static_cast<Eigen::Index>(truth_cells.size() * nmoments),
          static_cast<Eigen::Index>(truth_cells.size() * nactive))};

  for (const auto &[cell, index] : truth_offset) {
    const auto weight = sums.generated.find(cell);
    if (weight == sums.generated.end() || !(weight->second > 0.0) ||
        !std::isfinite(weight->second)) {
      throw std::invalid_argument(
          "harmonic::HarmonicMeasurement has a non-positive response truth "
          "cell weight");
    }
    const auto block = sums.fiducial.find(cell);
    if (block != sums.fiducial.end()) {
      response.fiducial.block(static_cast<Eigen::Index>(index * nmoments),
                              static_cast<Eigen::Index>(index * nactive),
                              static_cast<Eigen::Index>(nmoments),
                              static_cast<Eigen::Index>(nactive)) =
          block->second / weight->second;
    }
  }

  for (const auto &[cells, block] : sums.detector) {
    const auto reco = detector_offset.find(cells.first);
    const auto truth = truth_offset.find(cells.second);
    if (reco == detector_offset.end() || truth == truth_offset.end()) {
      continue;
    }
    const double weight = sums.generated.at(cells.second);
    response.detector.block(static_cast<Eigen::Index>(reco->second * nmoments),
                            static_cast<Eigen::Index>(truth->second * nactive),
                            static_cast<Eigen::Index>(nmoments),
                            static_cast<Eigen::Index>(nactive)) =
        block / weight;
  }
  return response;
}

// Require a positive response truth weight in every delete-group replica
void ValidateResponseJackknife(const ResponseSums &total,
                               const std::vector<ResponseSums> &deleted) {
  for (const auto &group : indices(deleted)) {
    for (const auto &[cell, weight] : deleted[group].generated) {
      const double retained = total.generated.at(cell) - weight;
      if (!(retained > 0.0) || !std::isfinite(retained)) {
        throw std::invalid_argument(
            "harmonic::HarmonicMeasurement response jackknife group " +
            std::to_string(group) + " leaves truth cell " + gra::aux::dvec2str(cell) +
            " without a positive finite weight, increase response statistics "
            "or use coarser truth bins");
      }
    }
  }
}

// Build detector-level moments and their exact weighted-event covariance
// m = sum_i w_i B_i, Cov(m) = sum_i w_i^2 B_i B_i^T
DataSums BuildDataSums(const std::vector<DataEvent> &data,
                       const PhaseSpaceGrid &grid,
                       const std::vector<Cell> &detector_cells,
                       const std::vector<std::size_t> &active, int lmax) {
  const std::size_t nactive = active.size();
  const std::map<Cell, std::size_t> offsets = CellOffsets(detector_cells);
  DataSums output{
      Eigen::VectorXd::Zero(
          static_cast<Eigen::Index>(detector_cells.size() * nactive)),
      Eigen::MatrixXd::Zero(
          static_cast<Eigen::Index>(detector_cells.size() * nactive),
          static_cast<Eigen::Index>(detector_cells.size() * nactive)),
      {},
      0,
      0.0,
      0.0};
  std::map<Cell, statistics::WeightedVectorSums> sums;
  statistics::ScaledWeightSums weights;
  for (const DataEvent &event : data) {
    event.reco.Validate(grid.Axes().size(), "harmonic::HarmonicMeasurement data reco");
    const auto cell = grid.Locate(event.reco.z);
    if (!cell || std::fpclassify(event.reco.weight) == FP_ZERO) { continue; }
    if (!offsets.contains(*cell)) {
      throw std::invalid_argument("harmonic::HarmonicMeasurement data occupy a cell without response");
    }
    sums.try_emplace(*cell, nactive).first->second.Add(Basis(event.reco, active, lmax), event.reco.weight);
    weights.Add(event.reco.weight);
  }
  for (const auto &[cell, sum] : sums) {
    const Eigen::Index start = static_cast<Eigen::Index>(offsets.at(cell) * nactive);
    const auto values = sum.Sum();
    output.values.segment(start, nactive) = Eigen::Map<const Eigen::VectorXd>(values.data(), nactive);
    output.covariance.block(start, start, nactive, nactive) = sum.Covariance().ToEigen();
    CellEstimate &estimate = output.cells[cell];
    estimate.sum_weight = static_cast<double>(sum.Weights().Sum());
    estimate.sum_weight2 = static_cast<double>(sum.Weights().SquareSum());
    estimate.entries = sum.Weights().Count();
  }
  output.accepted = weights.Count();
  output.sum_weight = static_cast<double>(weights.Sum());
  output.sum_weight2 = static_cast<double>(weights.SquareSum());
  output.effective_entries = static_cast<double>(weights.EffectiveSampleSize());
  return output;
}

// Expand an active cell vector and covariance into the full harmonic basis
CellEstimate ExpandCell(const Eigen::VectorXd &values,
                        const Eigen::MatrixXd &covariance,
                        const std::vector<std::size_t> &active,
                        std::size_t ncoef) {
  CellEstimate output;
  output.value.assign(ncoef, 0.0);
  output.covariance = MMatrix<double>(ncoef, ncoef, 0.0);
  for (const auto &row : indices(active)) {
    output.value[active[row]] = values[static_cast<Eigen::Index>(row)];
    for (const auto &col : indices(active)) {
      output.covariance[active[row]][active[col]] = covariance(
          static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col));
    }
  }
  output.valid = true;
  return output;
}

// Add delete-group jackknife covariance around the replica mean
Eigen::MatrixXd
JackknifeCovariance(const std::vector<Eigen::VectorXd> &replicas) {
  if (replicas.empty()) { return {}; }
  std::vector<std::vector<double>> sample;
  sample.reserve(replicas.size());
  for (const Eigen::VectorXd &replica : replicas) {
    std::vector<double> value(static_cast<std::size_t>(replica.size()), 0.0);
    for (std::size_t i = 0; i < value.size(); ++i) { value[i] = replica[static_cast<Eigen::Index>(i)]; }
    sample.push_back(std::move(value));
  }
  const MMatrix<double> covariance = statistics::JackknifeCovariance(sample);
  Eigen::MatrixXd       output(static_cast<Eigen::Index>(covariance.size_row()),
                               static_cast<Eigen::Index>(covariance.size_col()));
  for (std::size_t row = 0; row < covariance.size_row(); ++row) {
    for (std::size_t col = 0; col < covariance.size_col(); ++col) {
      output(static_cast<Eigen::Index>(row), static_cast<Eigen::Index>(col)) = covariance[row][col];
    }
  }
  return output;
}

// Build the angular grid used to reject non-physical fitted intensities
Eigen::MatrixXd PositivityBasis(const FitConfig &config,
                                const std::vector<std::size_t> &active) {
  const std::size_t rows =
      config.eml_positivity_costheta * config.eml_positivity_phi;
  Eigen::MatrixXd output(static_cast<Eigen::Index>(rows),
                         static_cast<Eigen::Index>(active.size()));
  std::size_t row = 0;
  for (std::size_t i = 0; i < config.eml_positivity_costheta; ++i) {
    const double costheta =
        -1.0 + 2.0 * static_cast<double>(i) /
                   static_cast<double>(config.eml_positivity_costheta - 1);
    for (std::size_t j = 0; j < config.eml_positivity_phi; ++j) {
      const double phi = -gra::math::PI +
                         2.0 * gra::math::PI * (static_cast<double>(j) + 0.5) /
                             static_cast<double>(config.eml_positivity_phi);
      const Observation point{costheta, phi, {}, 1.0};
      output.row(static_cast<Eigen::Index>(row)) =
          Basis(point, active, config.lmax).transpose();
      ++row;
    }
  }
  return output;
}

// Compute the active-vector location of the yield coefficient
std::size_t YieldIndex(const std::vector<std::size_t> &active) {
  const auto iterator = std::find(active.begin(), active.end(), 0);
  if (iterator == active.end()) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement requires the l=0,m=0 coefficient");
  }
  return static_cast<std::size_t>(std::distance(active.begin(), iterator));
}

// Group selected data by reconstructed conditional cell for the EML objective
EMLModel BuildEMLModel(const std::vector<DataEvent> &data,
                       const PhaseSpaceGrid &grid,
                       const std::vector<Cell> &detector_cells,
                       const std::vector<std::size_t> &active,
                       const std::vector<std::size_t> &moments,
                       const LinearResponse &response,
                       const FitConfig &config) {
  const std::map<Cell, std::size_t> offsets = CellOffsets(detector_cells);
  std::vector<std::vector<const DataEvent *>> grouped(detector_cells.size());
  for (const DataEvent &event : data) {
    event.reco.Validate(grid.Axes().size(),
                        "harmonic::HarmonicMeasurement EML data");
    const std::optional<Cell> cell = grid.Locate(event.reco.z);
    if (!cell.has_value() || std::fpclassify(event.reco.weight) == FP_ZERO) {
      continue;
    }
    if (event.reco.weight < 0.0) {
      throw std::invalid_argument(
          "harmonic::HarmonicMeasurement EML does not accept signed data "
          "weights, use ALGEBRAIC");
    }
    const auto offset = offsets.find(*cell);
    if (offset == offsets.end()) {
      throw std::invalid_argument(
          "harmonic::HarmonicMeasurement data occupy a cell without response");
    }
    grouped[offset->second].push_back(&event);
  }

  EMLModel model;
  model.response = response.detector;
  model.angular_basis = PositivityBasis(config, active);
  model.detector_basis = PositivityBasis(config, moments);
  model.nactive = active.size();
  model.nmoments = moments.size();
  model.zero_index = YieldIndex(active);
  for (const auto &cell : indices(grouped)) {
    if (grouped[cell].empty()) {
      continue;
    }
    EMLCellData selected;
    selected.detector_cell = cell;
    selected.basis =
        Eigen::MatrixXd(static_cast<Eigen::Index>(grouped[cell].size()),
                        static_cast<Eigen::Index>(moments.size()));
    selected.weight =
        Eigen::VectorXd(static_cast<Eigen::Index>(grouped[cell].size()));
    for (const auto &row : indices(grouped[cell])) {
      selected.basis.row(static_cast<Eigen::Index>(row)) =
          Basis(grouped[cell][row]->reco, moments, config.lmax).transpose();
      selected.weight[static_cast<Eigen::Index>(row)] =
          grouped[cell][row]->reco.weight;
    }
    model.data.push_back(std::move(selected));
  }
  return model;
}

// Compute the lowest angular-flat and reconstructed intensities on the grid
std::pair<double, double> MinimumIntensities(const EMLModel &model,
                                             const Eigen::VectorXd &flat) {
  double minimum_flat = std::numeric_limits<double>::infinity();
  for (Eigen::Index start = 0; start < flat.size();
       start += static_cast<Eigen::Index>(model.nactive)) {
    minimum_flat =
        std::min(minimum_flat,
                 (model.angular_basis *
                  flat.segment(start, static_cast<Eigen::Index>(model.nactive)))
                     .minCoeff());
  }
  const Eigen::VectorXd detector = model.response * flat;
  double minimum_detector = std::numeric_limits<double>::infinity();
  for (Eigen::Index start = 0; start < detector.size();
       start += static_cast<Eigen::Index>(model.nmoments)) {
    minimum_detector = std::min(
        minimum_detector,
        (model.detector_basis *
         detector.segment(start, static_cast<Eigen::Index>(model.nmoments)))
            .minCoeff());
  }
  return {minimum_flat, minimum_detector};
}

// Evaluate the weighted detector-space extended negative log likelihood
// NLL = sum_c mu_c-sum_i w_i log I_c(Omega_i)
double EMLObjective(const EMLModel &model, const Eigen::VectorXd &flat,
                    double minimum_intensity) {
  if (!flat.allFinite()) {
    return 1e100;
  }
  const auto minimum = MinimumIntensities(model, flat);
  if (minimum.first < -minimum_intensity ||
      minimum.second < -minimum_intensity) {
    return 1e100;
  }
  const Eigen::VectorXd detector = model.response * flat;
  double objective = 0.0;
  for (Eigen::Index start = 0; start < detector.size();
       start += static_cast<Eigen::Index>(model.nmoments)) {
    objective += detector[start + static_cast<Eigen::Index>(model.zero_index)];
  }
  if (!(objective > 0.0) || !std::isfinite(objective)) {
    return 1e100;
  }
  for (const EMLCellData &cell : model.data) {
    const Eigen::Index start =
        static_cast<Eigen::Index>(cell.detector_cell * model.nmoments);
    const Eigen::VectorXd intensity =
        cell.basis *
        detector.segment(start, static_cast<Eigen::Index>(model.nmoments));
    if ((intensity.array() <= minimum_intensity).any() ||
        !intensity.allFinite()) {
      return 1e100;
    }
    objective -= (cell.weight.array() * intensity.array().log()).matrix().sum();
  }
  return std::isfinite(objective) ? objective : 1e100;
}

// Evaluate the analytic EML gradient in the physical parameter region
// grad NLL = R^T[e_00-sum_i w_i B_i/I_i]
Eigen::VectorXd EMLGradient(const EMLModel &model, const Eigen::VectorXd &flat,
                            double minimum_intensity) {
  Eigen::VectorXd gradient = Eigen::VectorXd::Zero(flat.size());
  if (EMLObjective(model, flat, minimum_intensity) >= 1e99) {
    return gradient;
  }
  const Eigen::VectorXd detector = model.response * flat;
  Eigen::VectorXd detector_gradient = Eigen::VectorXd::Zero(detector.size());
  for (Eigen::Index start = 0; start < detector.size();
       start += static_cast<Eigen::Index>(model.nmoments)) {
    detector_gradient[start + static_cast<Eigen::Index>(model.zero_index)] =
        1.0;
  }
  for (const EMLCellData &cell : model.data) {
    const Eigen::Index start =
        static_cast<Eigen::Index>(cell.detector_cell * model.nmoments);
    const Eigen::VectorXd intensity =
        cell.basis *
        detector.segment(start, static_cast<Eigen::Index>(model.nmoments));
    detector_gradient.segment(start, static_cast<Eigen::Index>(model.nmoments))
        .noalias() -= cell.basis.transpose() *
                      (cell.weight.array() / intensity.array()).matrix();
  }
  gradient.noalias() = model.response.transpose() * detector_gradient;
  return gradient;
}

// Dampen algebraic anisotropies until they define a physical EML start
Eigen::VectorXd FeasibleStart(const EMLModel &model,
                              const Eigen::VectorXd &algebraic,
                              double data_weight, double minimum_intensity) {
  Eigen::VectorXd start = algebraic;
  const std::size_t truth_cells =
      static_cast<std::size_t>(start.size()) / model.nactive;
  const double average_yield =
      std::max(1.0, data_weight / static_cast<double>(truth_cells));
  for (std::size_t cell = 0; cell < truth_cells; ++cell) {
    const Eigen::Index zero =
        static_cast<Eigen::Index>(cell * model.nactive + model.zero_index);
    start[zero] = std::max(
        {minimum_intensity * 100.0, average_yield, std::abs(start[zero])});
  }
  for (std::size_t iteration = 0; iteration < 64; ++iteration) {
    if (EMLObjective(model, start, minimum_intensity) < 1e99) {
      return start;
    }
    for (std::size_t cell = 0; cell < truth_cells; ++cell) {
      for (std::size_t coefficient = 0; coefficient < model.nactive;
           ++coefficient) {
        if (coefficient != model.zero_index) {
          start[static_cast<Eigen::Index>(cell * model.nactive +
                                          coefficient)] *= 0.5;
        }
      }
    }
  }
  throw std::invalid_argument(
      "harmonic::HarmonicMeasurement cannot construct a positive EML model "
      "from this truncated response");
}

// Calculate the weighted sandwich covariance of one converged EML fit
// Cov = H^+ J H^{+T}
Eigen::MatrixXd EMLCovariance(const EMLModel &model,
                              const Eigen::VectorXd &flat) {
  const Eigen::Index parameters = flat.size();
  Eigen::MatrixXd hessian = Eigen::MatrixXd::Zero(parameters, parameters);
  Eigen::MatrixXd score = Eigen::MatrixXd::Zero(parameters, parameters);
  const Eigen::VectorXd detector = model.response * flat;
  for (const EMLCellData &cell : model.data) {
    const Eigen::Index start =
        static_cast<Eigen::Index>(cell.detector_cell * model.nmoments);
    const Eigen::MatrixXd block = model.response.block(
        start, 0, static_cast<Eigen::Index>(model.nmoments), parameters);
    const Eigen::VectorXd intensity =
        cell.basis *
        detector.segment(start, static_cast<Eigen::Index>(model.nmoments));
    const Eigen::ArrayXd inverse_square = intensity.array().square().inverse();
    const Eigen::MatrixXd local_hessian =
        cell.basis.transpose() *
        (cell.weight.array() * inverse_square).matrix().asDiagonal() *
        cell.basis;
    const Eigen::MatrixXd local_score =
        cell.basis.transpose() *
        (cell.weight.array().square() * inverse_square).matrix().asDiagonal() *
        cell.basis;
    hessian.noalias() += block.transpose() * local_hessian * block;
    score.noalias() += block.transpose() * local_score * block;
  }
  const HarmonicPseudoInverse inverse = PseudoInverse(hessian, 0.0);
  if (inverse.numerical_rank < static_cast<std::size_t>(parameters)) {
    throw std::runtime_error(
        "harmonic::HarmonicMeasurement EML covariance is not identifiable");
  }
  return inverse.value * score * inverse.value.transpose();
}

// Fit the global flat-space coefficients with the detector response folded in
EstimatorResult FitEML(const EMLModel &model,
                       const Eigen::VectorXd &algebraic_start,
                       double data_weight, const FitConfig &config,
                       bool calculate_covariance) {
  const Eigen::VectorXd start = FeasibleStart(
      model, algebraic_start, data_weight, config.eml_min_intensity);
  std::unique_ptr<ROOT::Math::Minimizer> minimizer(
      ROOT::Math::Factory::CreateMinimizer("Minuit2", "Migrad"));
  if (!minimizer) {
    throw std::runtime_error(
        "harmonic::HarmonicMeasurement failed to create Minuit2");
  }
  const std::function<double(const double *)> value =
      [&model, &config](const double *parameters) {
        const Eigen::Map<const Eigen::VectorXd> values(parameters,
                                                       model.response.cols());
        return EMLObjective(model, values, config.eml_min_intensity);
      };
  const std::function<void(const double *, double *)> gradient =
      [&model, &config](const double *parameters, double *output) {
        const Eigen::Map<const Eigen::VectorXd> values(parameters,
                                                       model.response.cols());
        Eigen::Map<Eigen::VectorXd> mapped(output, model.response.cols());
        mapped = EMLGradient(model, values, config.eml_min_intensity);
      };
  ROOT::Math::GradFunctor objective(
      value, static_cast<unsigned int>(start.size()), gradient);
  minimizer->SetFunction(objective);
  minimizer->SetMaxFunctionCalls(
      static_cast<unsigned int>(config.eml_max_calls));
  minimizer->SetMaxIterations(static_cast<unsigned int>(config.eml_max_calls));
  minimizer->SetTolerance(1e-6);
  minimizer->SetPrintLevel(-1);
  minimizer->SetErrorDef(0.5);

  for (const auto &index : indices(start)) {
    const auto zero = static_cast<Eigen::Index>(
        static_cast<std::size_t>(index) / model.nactive * model.nactive + model.zero_index);
    const double step = std::max(1e-4, 1e-3 * start[zero]);
    const auto parameter = static_cast<unsigned int>(index);
    const std::string name = "t_" + std::to_string(index);
    const bool configured = index == zero
        ? minimizer->SetLowerLimitedVariable(parameter, name, start[index], step,
                                             config.eml_min_intensity)
        : minimizer->SetVariable(parameter, name, start[index], step);
    if (!configured) {
      throw std::runtime_error(
          "harmonic::HarmonicMeasurement failed to configure EML parameter");
    }
  }
  if (!minimizer->Minimize() || minimizer->Status() != 0 ||
      !std::isfinite(minimizer->MinValue())) {
    throw std::runtime_error(
        "harmonic::HarmonicMeasurement EML did not converge, status " +
        std::to_string(minimizer->Status()));
  }

  EstimatorResult output;
  output.flat = Eigen::Map<const Eigen::VectorXd>(minimizer->X(), start.size());
  output.objective = minimizer->MinValue();
  const auto minimum = MinimumIntensities(model, output.flat);
  output.minimum_flat_intensity = minimum.first;
  output.minimum_detector_intensity = minimum.second;
  if (minimum.first < -config.eml_min_intensity ||
      minimum.second < -config.eml_min_intensity) {
    throw std::runtime_error(
        "harmonic::HarmonicMeasurement EML violated positivity");
  }
  if (calculate_covariance) {
    output.flat_covariance = EMLCovariance(model, output.flat);
  }
  return output;
}

// Require every fitted truth coefficient to be identifiable in the response
void ValidateResponseInverse(const HarmonicPseudoInverse &inverse,
                             Eigen::Index rows, Eigen::Index columns) {
  if (rows < columns ||
      inverse.numerical_rank < static_cast<std::size_t>(columns)) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement response is not identifiable");
  }
  if (inverse.retained_rank == 0) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement regularization removed all modes");
  }
}

// Solve the weighted moment equations with the shared pseudoinverse
EstimatorResult FitAlgebraic(const HarmonicPseudoInverse &inverse,
                             const DataSums &observed) {
  EstimatorResult output;
  output.flat = inverse.value * observed.values;
  output.flat_covariance =
      inverse.value * observed.covariance * inverse.value.transpose();
  return output;
}

// Insert one active global level into per-cell full-basis estimates
void FillCellEstimates(const std::vector<Cell> &cells,
                       const Eigen::VectorXd &values,
                       const Eigen::MatrixXd &covariance,
                       const std::vector<std::size_t> &active,
                       std::size_t ncoef,
                       std::map<Cell, CellEstimate> &output) {
  const std::size_t nactive = active.size();
  for (const auto &cell : indices(cells)) {
    const Eigen::Index start = static_cast<Eigen::Index>(cell * nactive);
    output[cells[cell]] = ExpandCell(
        values.segment(start, static_cast<Eigen::Index>(nactive)),
        covariance.block(start, start, static_cast<Eigen::Index>(nactive),
                         static_cast<Eigen::Index>(nactive)),
        active, ncoef);
  }
}

// Compute active global values in the stored cell order
Eigen::VectorXd GlobalValues(const std::vector<Cell> &cells,
                             const std::map<Cell, CellEstimate> &estimates,
                             const std::vector<std::size_t> &active) {
  Eigen::VectorXd output = Eigen::VectorXd::Zero(
      static_cast<Eigen::Index>(cells.size() * active.size()));
  for (const auto &cell : indices(cells)) {
    const CellEstimate &estimate = estimates.at(cells[cell]);
    for (const auto &index : indices(active)) {
      output[static_cast<Eigen::Index>(cell * active.size() + index)] =
          estimate.value[active[index]];
    }
  }
  return output;
}

// Scale one per-cell level and add its normalization covariance
void NormalizeCells(std::map<Cell, CellEstimate> &estimates, double divisor,
                    double relative_uncertainty) {
  for (auto &[cell, estimate] : estimates) {
    static_cast<void>(cell);
    for (double &value : estimate.value) { value /= divisor; }
    estimate.covariance = estimate.covariance / divisor / divisor;
    estimate.covariance.AddOuterProduct(
        estimate.value, estimate.value,
        relative_uncertainty * relative_uncertainty);
  }
}

// Scale one global covariance and add its normalization covariance
void NormalizeGlobal(MMatrix<double> &covariance,
                     const Eigen::VectorXd &scaled_values, double divisor,
                     double relative_uncertainty) {
  covariance = covariance / divisor / divisor;
  covariance.AddOuterProduct(scaled_values, scaled_values,
                             relative_uncertainty * relative_uncertainty);
}

} // namespace

// Compute the canonical coordinate ordering for one measurement mode
std::vector<Coordinate> Coordinates(MeasurementMode mode) {
  if (mode == MeasurementMode::Central) {
    return {Coordinate::Mass, Coordinate::Momentum, Coordinate::Rapidity};
  }
  return {Coordinate::Mass, Coordinate::Rapidity, Coordinate::AbsT1,
          Coordinate::AbsT2, Coordinate::DeltaPhiPP};
}

// Compute a stable printable measurement-mode name
std::string ToString(MeasurementMode mode) {
  return mode == MeasurementMode::Central ? "CENTRAL" : "TAGGED";
}

// Parse a measurement mode without accepting implicit aliases
MeasurementMode ParseMeasurementMode(const std::string &value) {
  if (value == "CENTRAL") {
    return MeasurementMode::Central;
  }
  if (value == "TAGGED") {
    return MeasurementMode::Tagged;
  }
  throw std::invalid_argument(
      "harmonic::ParseMeasurementMode expects CENTRAL or TAGGED");
}

// Compute a stable printable sample-type name
std::string ToString(SampleType type) {
  if (type == SampleType::Data) {
    return "DATA";
  }
  if (type == SampleType::ResponseMC) {
    return "RESPONSE_MC";
  }
  return "CLOSURE_MC";
}

// Parse an input sample type without accepting implicit aliases
SampleType ParseSampleType(const std::string &value) {
  if (value == "DATA") {
    return SampleType::Data;
  }
  if (value == "RESPONSE_MC") {
    return SampleType::ResponseMC;
  }
  if (value == "CLOSURE_MC") {
    return SampleType::ClosureMC;
  }
  throw std::invalid_argument(
      "harmonic::ParseSampleType expects DATA, RESPONSE_MC or CLOSURE_MC");
}

// Compute a stable printable phase-space-level name
std::string ToString(PhaseSpaceLevel level) {
  if (level == PhaseSpaceLevel::Flat) {
    return "ANGULAR_FLAT";
  }
  if (level == PhaseSpaceLevel::Fiducial) {
    return "FIDUCIAL";
  }
  return "DETECTOR";
}

// Compute a stable printable estimator name
std::string ToString(HarmonicEstimator estimator) {
  return estimator == HarmonicEstimator::Algebraic ? "ALGEBRAIC" : "EML";
}

// Parse an estimator without accepting implicit aliases
HarmonicEstimator ParseHarmonicEstimator(const std::string &value) {
  if (value == "ALGEBRAIC") {
    return HarmonicEstimator::Algebraic;
  }
  if (value == "EML") {
    return HarmonicEstimator::EML;
  }
  throw std::invalid_argument(
      "harmonic::ParseHarmonicEstimator expects ALGEBRAIC or EML");
}

// Compute a stable printable coordinate name
std::string ToString(Coordinate coordinate) {
  if (coordinate == Coordinate::Mass) {
    return "M";
  }
  if (coordinate == Coordinate::Momentum) {
    return "PT";
  }
  if (coordinate == Coordinate::Rapidity) {
    return "Y";
  }
  if (coordinate == Coordinate::AbsT1) {
    return "ABST1";
  }
  if (coordinate == Coordinate::AbsT2) {
    return "ABST2";
  }
  return "DPHI_PP";
}

// Parse a conditional coordinate without accepting implicit aliases
Coordinate ParseCoordinate(const std::string &value) {
  if (value == "M") {
    return Coordinate::Mass;
  }
  if (value == "PT") {
    return Coordinate::Momentum;
  }
  if (value == "Y") {
    return Coordinate::Rapidity;
  }
  if (value == "ABST1") {
    return Coordinate::AbsT1;
  }
  if (value == "ABST2") {
    return Coordinate::AbsT2;
  }
  if (value == "DPHI_PP") {
    return Coordinate::DeltaPhiPP;
  }
  throw std::invalid_argument(
      "harmonic::ParseCoordinate received an unknown coordinate");
}

// Validate the finite non-empty axis definition
void Axis::Validate() const {
  const double width = (max - min) / static_cast<double>(bins);
  if (bins == 0 || !std::isfinite(width) || !(min + width > min) || !(max - width < max) ||
      !std::isfinite(min) || !std::isfinite(max) || !(min < max)) {
    throw std::invalid_argument("harmonic::Axis has invalid binning");
  }
  if (coordinate == Coordinate::DeltaPhiPP &&
      (min < -gra::math::PI || max > gra::math::PI)) {
    throw std::invalid_argument(
        "harmonic::Axis DPHI_PP must use the signed [-pi,pi] convention");
  }
}

// Validate all angular, kinematic and weight components
void Observation::Validate(std::size_t dimensions,
                           const std::string &context) const {
  if (z.size() != dimensions || !std::isfinite(costheta) || costheta < -1.0 ||
      costheta > 1.0 || !std::isfinite(phi) || phi < -gra::math::PI ||
      phi > gra::math::PI || !std::isfinite(weight)) {
    throw std::invalid_argument(context + " has invalid dimensions or angles");
  }
  if (!gra::AllFinite(z)) {
    throw std::invalid_argument(context + " has non-finite coordinates");
  }
}

// Apply only the forward requirement present in the tagged mode
bool FiducialDecision::Pass(MeasurementMode mode) const {
  return central && (mode == MeasurementMode::Central || forward);
}

// Test reconstructed detector selection independently of truth fiducial cuts
bool ResponseEvent::DetectorAccepted() const {
  return reco.has_value() && reco_selected;
}

// Validate particle efficiencies and momentum resolutions
void ToyResponseConfig::Validate() const {
  for (const double efficiency : {pion_efficiency, proton_efficiency}) {
    if (!std::isfinite(efficiency) || efficiency < 0.0 || efficiency > 1.0) {
      throw std::invalid_argument("harmonic::ToyResponseConfig has invalid particle efficiency");
    }
  }
  for (const double sigma : {pion_logpt_sigma, pion_eta_sigma, pion_phi_sigma,
                              proton_pt_sigma, proton_logpz_sigma}) {
    if (!std::isfinite(sigma) || sigma < 0.0) {
      throw std::invalid_argument("harmonic::ToyResponseConfig has invalid momentum resolution");
    }
  }
}

// Compute unit efficiency
double IdentityResponse::Efficiency(const EventKinematics &truth) const {
  static_cast<void>(truth);
  return 1.0;
}

// Compute the unchanged generated particle momenta
EventKinematics IdentityResponse::Reconstruct(const EventKinematics &truth,
                                               std::uint64_t event_key) const {
  static_cast<void>(event_key);
  return truth;
}

// Construct a validated analytic response
ToyResponse::ToyResponse(const ToyResponseConfig &config) : config_(config) {
  config_.Validate();
}

// Compute the efficiency to reconstruct all required charged particles
double ToyResponse::Efficiency(const EventKinematics &truth) const {
  const double central = config_.pion_efficiency * config_.pion_efficiency;
  return truth.mode == MeasurementMode::Tagged
      ? central * config_.proton_efficiency * config_.proton_efficiency
      : central;
}

namespace {

// Recompute a measured energy from finite spatial momentum and the input mass
M4Vec MeasuredMomentum(const M4Vec &truth, double px, double py, double pz) {
  if (!std::isfinite(truth.E()) || !(truth.E() > 0.0) ||
      !std::isfinite(truth.M2()) || truth.M2() < 0.0 ||
      !std::isfinite(px) || !std::isfinite(py) || !std::isfinite(pz)) {
    throw std::invalid_argument("harmonic::ToyResponse has invalid particle momentum");
  }
  M4Vec reco;
  reco.SetPxPyPzM(px, py, pz, truth.M());
  if (!std::isfinite(reco.E())) {
    throw std::invalid_argument("harmonic::ToyResponse momentum resolution overflow");
  }
  return reco;
}

// Smear central tracks in log(pT), eta and phi without angular clipping
M4Vec SmearPion(const M4Vec &truth, const ToyResponseConfig &config,
                std::mt19937_64 &generator, std::normal_distribution<double> &normal) {
  const double sigma = config.pion_logpt_sigma;
  const double scale = std::exp(sigma * normal(generator) - 0.5 * sigma * sigma);
  const double eta = config.pion_eta_sigma * normal(generator);
  const double phi = config.pion_phi_sigma * normal(generator);
  M4Vec rotated = truth;
  rotated.RotateZ(phi);
  // Use pT sinh(eta + delta_eta) in a form also defined on the beam axis
  const double pz = scale * (truth.Pz() * std::cosh(eta) + truth.P3mod() * std::sinh(eta));
  return MeasuredMomentum(truth, scale * rotated.Px(), scale * rotated.Py(), pz);
}

// Smear forward transverse momentum and longitudinal magnitude in a fixed arm
M4Vec SmearProton(const M4Vec &truth, const ToyResponseConfig &config,
                  std::mt19937_64 &generator, std::normal_distribution<double> &normal) {
  const double px = truth.Px() + config.proton_pt_sigma * normal(generator);
  const double py = truth.Py() + config.proton_pt_sigma * normal(generator);
  const double sigma = config.proton_logpz_sigma;
  const double pz = truth.Pz() * std::exp(sigma * normal(generator) - 0.5 * sigma * sigma);
  return MeasuredMomentum(truth, px, py, pz);
}

} // namespace

// Smear measured momenta while preserving particle masses
EventKinematics ToyResponse::Reconstruct(const EventKinematics &truth,
                                           std::uint64_t event_key) const {
  // Keep the smearing stream independent of the acceptance decision
  std::mt19937_64 generator(MixSeed(MixSeed(config_.seed, 0x736d656172ULL), event_key));
  std::normal_distribution<double> normal(0.0, 1.0);
  EventKinematics reco = truth;
  reco.pip = SmearPion(truth.pip, config_, generator, normal);
  reco.pim = SmearPion(truth.pim, config_, generator, normal);
  if (truth.mode == MeasurementMode::Tagged && truth.has_forward_protons) {
    reco.proton_plus = SmearProton(truth.proton_plus, config_, generator, normal);
    reco.proton_minus = SmearProton(truth.proton_minus, config_, generator, normal);
  }
  return reco;
}

// Apply detector efficiency and response without shared random state
std::optional<EventKinematics> SimulateDetector(const DetectorResponseModel &model,
                                                const EventKinematics &truth,
                                                std::uint64_t event_key,
                                                std::uint64_t seed) {
  const double efficiency = model.Efficiency(truth);
  if (!std::isfinite(efficiency) || efficiency < 0.0 || efficiency > 1.0) {
    throw std::invalid_argument(
        "harmonic::SimulateDetector received an invalid efficiency");
  }
  // Use a distinct stream even when the response uses the same seed
  std::mt19937_64 generator(MixSeed(MixSeed(seed, 0x616363657074ULL), event_key));
  std::uniform_real_distribution<double> uniform(0.0, 1.0);
  if (uniform(generator) >= efficiency) {
    return std::nullopt;
  }
  return model.Reconstruct(truth, event_key);
}

// Construct and validate a mode-specific conditional grid
PhaseSpaceGrid::PhaseSpaceGrid(MeasurementMode mode,
                               const std::vector<Axis> &axes)
    : mode_(mode), axes_(axes) {
  const std::vector<Coordinate> expected = Coordinates(mode_);
  if (axes_.size() != expected.size()) {
    throw std::invalid_argument(
        "harmonic::PhaseSpaceGrid has incompatible dimensionality");
  }
  for (const auto &index : indices(axes_)) {
    axes_[index].Validate();
    if (axes_[index].coordinate != expected[index]) {
      throw std::invalid_argument(
          "harmonic::PhaseSpaceGrid axes are not in canonical order");
    }
  }
}

// Locate one point and reject coordinates outside the analysis range
std::optional<Cell> PhaseSpaceGrid::Locate(const std::vector<double> &z) const {
  if (z.size() != axes_.size()) {
    throw std::invalid_argument(
        "harmonic::PhaseSpaceGrid received incompatible coordinates");
  }
  Cell cell(axes_.size(), 0);
  for (const auto &index : indices(axes_)) {
    if (!std::isfinite(z[index]) || z[index] < axes_[index].min ||
        z[index] > axes_[index].max) {
      return std::nullopt;
    }
    if (std::is_eq(z[index] <=> axes_[index].max)) {
      cell[index] = axes_[index].bins - 1;
      continue;
    }
    const double fraction =
        (z[index] - axes_[index].min) / (axes_[index].max - axes_[index].min);
    cell[index] =
        std::min(axes_[index].bins - 1,
                 static_cast<std::size_t>(
                     fraction * static_cast<double>(axes_[index].bins)));
    // Resolve rounding against the same bin boundaries used in the measurement output
    const double width = (axes_[index].max - axes_[index].min) / static_cast<double>(axes_[index].bins);
    if (cell[index] + 1 < axes_[index].bins &&
        z[index] >= axes_[index].min + width * static_cast<double>(cell[index] + 1)) { ++cell[index]; }
    if (cell[index] > 0 &&
        z[index] < axes_[index].min + width * static_cast<double>(cell[index])) { --cell[index]; }
  }
  return cell;
}

// Compute the measurement mode
MeasurementMode PhaseSpaceGrid::Mode() const { return mode_; }

// Compute the ordered axes
const std::vector<Axis> &PhaseSpaceGrid::Axes() const { return axes_; }

// Validate truncation, regularization and memory bounds
void FitConfig::Validate() const {
  if (lmax < 0 || !std::isfinite(svd_relative_cut) ||
      svd_relative_cut < 0.0 || response_jackknife_bins == 1 ||
      max_parameters == 0 || eml_max_calls == 0 ||
      eml_positivity_costheta < 2 || eml_positivity_phi < 2 ||
      !std::isfinite(eml_min_intensity) || eml_min_intensity <= 0.0) {
    throw std::invalid_argument("harmonic::FitConfig has invalid parameters");
  }
  const std::size_t side = static_cast<std::size_t>(lmax) + 1;
  if (side > static_cast<std::size_t>(std::numeric_limits<int>::max()) / side) {
    throw std::invalid_argument("harmonic::FitConfig harmonic basis exceeds integer indexing");
  }
  // Bound the full measured basis even when production symmetries reduce the fit
  if (side * side > max_parameters) {
    throw std::invalid_argument("harmonic::FitConfig one angular basis exceeds max_parameters");
  }
  if (eml_max_calls > std::numeric_limits<unsigned int>::max() ||
      eml_positivity_costheta > static_cast<std::size_t>(std::numeric_limits<Eigen::Index>::max()) /
                                    eml_positivity_phi) {
    throw std::invalid_argument("harmonic::FitConfig EML dimensions exceed integer indexing");
  }
}

// Convert event yields to a normalized measurement with correlated uncertainty
void ApplyNormalization(MeasurementResult &result, double divisor,
                        double relative_uncertainty) {
  if (!std::isfinite(divisor) || divisor <= 0.0 ||
      !std::isfinite(relative_uncertainty) || relative_uncertainty < 0.0) {
    throw std::invalid_argument(
        "harmonic::ApplyNormalization has invalid normalization");
  }
  const Eigen::VectorXd flat_values =
      GlobalValues(result.truth_cells, result.flat, result.active_indices) /
      divisor;
  const Eigen::VectorXd fiducial_values =
      GlobalValues(result.truth_cells, result.fiducial, result.moment_indices) /
      divisor;
  const Eigen::VectorXd detector_values =
      GlobalValues(result.detector_cells, result.detector,
                   result.moment_indices) /
      divisor;
  NormalizeGlobal(result.flat_covariance, flat_values, divisor,
                  relative_uncertainty);
  NormalizeGlobal(result.fiducial_covariance, fiducial_values, divisor,
                  relative_uncertainty);
  NormalizeGlobal(result.detector_covariance, detector_values, divisor,
                  relative_uncertainty);
  NormalizeCells(result.flat, divisor, relative_uncertainty);
  NormalizeCells(result.fiducial, divisor, relative_uncertainty);
  NormalizeCells(result.detector, divisor, relative_uncertainty);
}

// Construct an immutable measurement definition
HarmonicMeasurement::HarmonicMeasurement(const PhaseSpaceGrid &grid,
                                         const FitConfig &config)
    : grid_(grid), config_(config) {
  config_.Validate();
}

// Unfold detector moments to angular-flat space and project to fiducial space
MeasurementResult
HarmonicMeasurement::Fit(const std::vector<ResponseEvent> &response,
                         const std::vector<DataEvent> &data) const {
  if (response.empty() || data.empty()) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement requires response and data events");
  }
  const std::vector<std::size_t> active = ActiveIndices(config_);
  const std::size_t nactive = active.size();
  const std::size_t ncoef =
      static_cast<std::size_t>((config_.lmax + 1) * (config_.lmax + 1));

  std::vector<std::size_t> moments(ncoef);
  std::iota(moments.begin(), moments.end(), 0);

  const auto scales = ResponseScales(response, grid_);
  ResponseSums total;
  std::vector<ResponseSums> deleted(config_.response_jackknife_bins);
  for (const ResponseEvent &event : response) {
    AddResponseEvent(total, event, grid_, active, moments, config_.lmax, scales);
    if (!deleted.empty()) {
      AddResponseEvent(deleted[event.event_key % deleted.size()], event, grid_,
                       active, moments, config_.lmax, scales);
    }
  }
  if (!deleted.empty() &&
      std::any_of(deleted.begin(), deleted.end(),
                  [](const ResponseSums &block) {
                    return block.generated.empty();
                  })) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement has an empty response jackknife "
        "group, reduce response_jackknife_bins");
  }
  const std::vector<Cell> truth_cells = MapCells(total.generated);
  if (truth_cells.empty()) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement response does not cover the grid");
  }

  std::set<Cell> detector_cell_set;
  for (const auto &[cells, block] : total.detector) {
    static_cast<void>(block);
    detector_cell_set.insert(cells.first);
  }
  const std::vector<Cell> detector_cells(detector_cell_set.begin(),
                                         detector_cell_set.end());
  if (detector_cells.empty()) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement has no selected detector cells");
  }
  if (truth_cells.size() > config_.max_parameters / ncoef ||
      detector_cells.size() > config_.max_parameters / ncoef) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement exceeds the configured fit size");
  }

  const LinearResponse nominal =
      AssembleResponse(total, truth_cells, detector_cells, nactive, ncoef);
  ValidateResponseJackknife(total, deleted);
  const DataSums observed =
      BuildDataSums(data, grid_, detector_cells, moments, config_.lmax);
  if (observed.accepted == 0) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement has no data events inside the grid");
  }
  const HarmonicPseudoInverse inverse =
      PseudoInverse(nominal.detector, config_.svd_relative_cut);
  ValidateResponseInverse(inverse, nominal.detector.rows(),
                          nominal.detector.cols());
  EstimatorResult estimate = FitAlgebraic(inverse, observed);
  if (config_.estimator == HarmonicEstimator::EML) {
    if (!(observed.sum_weight > 0.0) ||
        observed.effective_entries <= static_cast<double>(nominal.detector.cols())) {
      throw std::invalid_argument(
          "harmonic::HarmonicMeasurement EML requires positive total weight "
          "and more effective events than fitted coefficients");
    }
    const EMLModel model =
        BuildEMLModel(data, grid_, detector_cells, active, moments, nominal, config_);
    estimate = FitEML(model, estimate.flat, observed.sum_weight, config_, true);
  }
  const Eigen::VectorXd flat_values = estimate.flat;
  const Eigen::MatrixXd flat_data_covariance = estimate.flat_covariance;
  const Eigen::VectorXd fiducial_values = nominal.fiducial * flat_values;
  const Eigen::MatrixXd fiducial_data_covariance =
      nominal.fiducial * flat_data_covariance * nominal.fiducial.transpose();
  const Eigen::VectorXd detector_values =
      config_.estimator == HarmonicEstimator::EML
          ? nominal.detector * flat_values
          : observed.values;
  const Eigen::MatrixXd detector_data_covariance =
      config_.estimator == HarmonicEstimator::EML
          ? nominal.detector * flat_data_covariance *
                nominal.detector.transpose()
          : observed.covariance;

  std::vector<Eigen::VectorXd> flat_replicas;
  std::vector<Eigen::VectorXd> fiducial_replicas;
  std::vector<Eigen::VectorXd> detector_replicas;
  flat_replicas.reserve(deleted.size());
  fiducial_replicas.reserve(deleted.size());
  detector_replicas.reserve(deleted.size());
  for (const ResponseSums &block : deleted) {
    if (block.generated.empty()) {
      continue;
    }
    ResponseSums retained = total;
    for (const auto &[cell, weight] : block.generated) {
      retained.generated[cell] -= weight;
    }
    for (const auto &[cell, matrix] : block.fiducial) {
      retained.fiducial[cell] -= matrix;
    }
    for (const auto &[cells, matrix] : block.detector) {
      retained.detector[cells] -= matrix;
    }
    const LinearResponse replica =
        AssembleResponse(retained, truth_cells, detector_cells, nactive, ncoef);
    const HarmonicPseudoInverse replica_inverse =
        PseudoInverse(replica.detector, config_.svd_relative_cut);
    ValidateResponseInverse(replica_inverse, replica.detector.rows(),
                            replica.detector.cols());
    EstimatorResult replica_estimate = FitAlgebraic(replica_inverse, observed);
    if (config_.estimator == HarmonicEstimator::EML) {
      const EMLModel replica_model =
          BuildEMLModel(data, grid_, detector_cells, active, moments, replica, config_);
      replica_estimate = FitEML(replica_model, flat_values, observed.sum_weight,
                                config_, false);
    }
    flat_replicas.push_back(replica_estimate.flat);
    fiducial_replicas.push_back(replica.fiducial * replica_estimate.flat);
    if (config_.estimator == HarmonicEstimator::EML) {
      detector_replicas.push_back(replica.detector * replica_estimate.flat);
    }
  }
  if (config_.response_jackknife_bins >= 2 && flat_replicas.size() < 2) {
    throw std::invalid_argument(
        "harmonic::HarmonicMeasurement has fewer than two usable response "
        "jackknife replicas");
  }

  Eigen::MatrixXd flat_covariance = flat_data_covariance;
  Eigen::MatrixXd fiducial_covariance = fiducial_data_covariance;
  Eigen::MatrixXd detector_covariance = detector_data_covariance;
  const Eigen::MatrixXd flat_response_covariance =
      JackknifeCovariance(flat_replicas);
  const Eigen::MatrixXd fiducial_response_covariance =
      JackknifeCovariance(fiducial_replicas);
  const Eigen::MatrixXd detector_response_covariance =
      JackknifeCovariance(detector_replicas);
  if (flat_response_covariance.size() != 0) {
    flat_covariance += flat_response_covariance;
  }
  if (fiducial_response_covariance.size() != 0) {
    fiducial_covariance += fiducial_response_covariance;
  }
  if (detector_response_covariance.size() != 0) {
    detector_covariance += detector_response_covariance;
  }

  MeasurementResult result;
  result.mode = grid_.Mode();
  result.estimator = config_.estimator;
  result.truth_cells = truth_cells;
  result.detector_cells = detector_cells;
  result.flat_covariance = MMatrix<double>::FromEigen(flat_covariance);
  result.fiducial_covariance = MMatrix<double>::FromEigen(fiducial_covariance);
  result.detector_covariance =
      MMatrix<double>::FromEigen(detector_covariance);
  result.active_indices = active;
  result.moment_indices = moments;
  result.active_coefficients = nactive;
  result.response_rows = static_cast<std::size_t>(nominal.detector.rows());
  result.response_columns = static_cast<std::size_t>(nominal.detector.cols());
  result.response_rank = inverse.numerical_rank;
  result.retained_singular_values = inverse.retained_rank;
  result.response_condition_number = inverse.condition_number;
  result.response_events = response.size();
  result.response_jackknife_replicas = flat_replicas.size();
  result.data_events = observed.accepted;
  result.data_sum_weight = observed.sum_weight;
  result.data_sum_weight2 = observed.sum_weight2;
  result.objective = estimate.objective;
  result.minimum_flat_intensity = estimate.minimum_flat_intensity;
  result.minimum_detector_intensity = estimate.minimum_detector_intensity;
  FillCellEstimates(truth_cells, flat_values, flat_covariance, active, ncoef,
                    result.flat);
  FillCellEstimates(truth_cells, fiducial_values, fiducial_covariance, moments,
                    ncoef, result.fiducial);

  for (const auto &cell : indices(detector_cells)) {
    const Eigen::Index start = static_cast<Eigen::Index>(cell * ncoef);
    CellEstimate estimate = ExpandCell(
        detector_values.segment(start, static_cast<Eigen::Index>(ncoef)),
        detector_covariance.block(start, start,
                                  static_cast<Eigen::Index>(ncoef),
                                  static_cast<Eigen::Index>(ncoef)),
        moments, ncoef);
    const auto summary = observed.cells.find(detector_cells[cell]);
    if (summary != observed.cells.end()) {
      estimate.sum_weight = summary->second.sum_weight;
      estimate.sum_weight2 = summary->second.sum_weight2;
      estimate.entries = summary->second.entries;
    }
    result.detector[detector_cells[cell]] = std::move(estimate);
  }
  return result;
}

} // namespace harmonic
} // namespace gra
