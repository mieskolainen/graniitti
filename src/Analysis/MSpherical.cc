// Functional methods for Spherical Harmonic Expansions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>
#include <limits>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MStatistics.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Analysis/MSpherical.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::PI;

namespace gra {
namespace spherical {

namespace {

// Validate the full spherical basis before integer indexing or matrix allocation
int BasisSize(int lmax) {
  const long long side = static_cast<long long>(lmax) + 1;
  if (side <= 0 || side > std::numeric_limits<int>::max() / side) {
    throw std::invalid_argument("spherical::BasisSize: invalid angular truncation");
  }
  return static_cast<int>(side * side);
}

// Validate a nominal event weight before statistical accumulation
double EventWeight(const Omega &event, const std::string &context) {
  if (!std::isfinite(event.weight)) {
    throw std::invalid_argument(context + ": non-finite event weight");
  }
  return event.weight;
}

// Validate the angular coordinates used by the spherical basis
void ValidateAngles(const Omega &event, const std::string &context) {
  if (!std::isfinite(event.costheta) || !std::isfinite(event.phi) ||
      event.costheta < -1.0 || event.costheta > 1.0) {
    throw std::invalid_argument(context + ": invalid angular coordinates");
  }
}

// Evaluate each normalized real spherical harmonic once per angular observation
std::vector<double> Basis(const Omega &event, int lmax) {
  ValidateAngles(event, "spherical::Basis");
  std::vector<double> basis(BasisSize(lmax));
  for (int l = 0; l <= lmax; ++l) {
    for (int m = -l; m <= l; ++m) {
      basis[LinearInd(l, m)] = math::Y_real_basis(event.costheta, event.phi, l, m);
    }
  }
  return basis;
}

// Test whether an event belongs to the requested analysis level
bool AcceptEvent(const Omega &event, const std::string &mode) {
  if (mode == "fla") {
    return true;
  }
  if (mode == "fid") {
    return event.fiducial;
  }
  if (mode == "det") {
    return event.fiducial && event.selected;
  }
  throw std::invalid_argument("spherical::AcceptEvent: unknown mode " + mode);
}

// Validate scale-normalized reference weights used by spherical responses
void ValidateReferenceWeights(const statistics::ScaledWeightSums &weights,
                              const std::string &context) {
  if (!weights.IsFinite() ||
      !weights.HasSignificantSignedSum(
          std::numeric_limits<double>::epsilon()) ||
      !(weights.Sum() > 0.0L)) {
    throw std::invalid_argument(
        context + ": generated reference weight must be finite and positive");
  }
}

} // namespace

// Monte Carlo integral I ~= V 1/N \sum_{i=1}^N f(\vec{x}_i) = V <f(x)>
// True integral being  I = \int_\Omega f(\vec{x}) d\vec{x}
//
// This inner-product (overlap integral) matrix is identity if the phase space
// is flat (uniform). In a geometrically restricted phase space, these basis
// functions get mixed => lost orthogonality

MixingEstimate GetGMixing(const std::vector<Omega> &events,
                          const std::vector<std::size_t> &ind, int LMAX,
                          const std::string &mode) {
  if (mode != "fla" && mode != "fid" && mode != "det") {
    throw std::invalid_argument("spherical::GetGMixing: Unknown mode " + mode);
  }
  if (ind.size() < 2) {
    throw std::invalid_argument(
        "spherical::GetGMixing: at least two generated events are required");
  }
  for (const std::size_t k : ind) {
    if (k >= events.size()) {
      throw std::out_of_range(
          "spherical::GetGMixing: event index outside input array");
    }
  }
  const int NCOEF = BasisSize(LMAX);

  std::cout << "GetGMixing: mode = " << mode << std::endl;
  std::cout << "Generated flat MC phase space events = " << ind.size()
            << std::endl;

  // Construct the efficiency coefficients EPSILON_LM with linear
  // indexing
  MMatrix<double> E(NCOEF, NCOEF, 0.0);
  MMatrix<statistics::CompensatedSum> weighted_sum(NCOEF, NCOEF);
  const std::size_t jackknife_blocks = std::min<std::size_t>(32, ind.size());
  std::vector<MMatrix<statistics::CompensatedSum>> block_sum(
      jackknife_blocks, MMatrix<statistics::CompensatedSum>(NCOEF, NCOEF));
  std::vector<statistics::CompensatedSum> block_weight(jackknife_blocks);

  // Loop over GENERATED MC events and do the integral in effect via
  // uniform MC sampling
  // We evaluate the integral:
  //
  // eta_LM = \int eta(Omega) Re Y_LM(Omega) dOmega, Omega = (costheta,phi)
  std::size_t fiducial = 0;
  std::size_t selected = 0;
  statistics::ScaledWeightSums reference_weights;
  for (const std::size_t k : ind) {
    reference_weights.Add(EventWeight(events[k], "spherical::GetGMixing"));
  }
  ValidateReferenceWeights(reference_weights, "spherical::GetGMixing");
  const long double weight_scale = reference_weights.Scale();
  const long double generated_weight_scaled = reference_weights.ScaledSum();

  const double VOL = 4.0 * PI; // [costheta] x [phi] plane area

  for (const auto &sample : indices(ind)) {
    const std::size_t k = ind[sample];
    const std::size_t block = sample % jackknife_blocks;
    const double weight = EventWeight(events[k], "spherical::GetGMixing");
    const long double scaled_weight =
        static_cast<long double>(weight) / weight_scale;
    block_weight[block].Add(scaled_weight);
    const bool fid = events[k].fiducial;
    const bool sel = events[k].selected;

    // Flat phase space
    if (mode == "fla") {
      // all fine
    }
    // Geometric acceptance ("Fiducial phase space")
    else if (mode == "fid") {
      if (!fid) {
        continue;
      } else {
        ++fiducial;
      }
    }
    // Geometric x Efficiency ("Detector level")
    else if (mode == "det") {
      if (!fid) {
        continue;
      } else {
        ++fiducial;
      }
      if (!sel) {
        continue;
      } else {
        ++selected;
      }
    }

    const auto basis = Basis(events[k], LMAX);
    for (const auto &row : indices(basis)) {
      for (const auto &col : indices(basis)) {
        const long double f = static_cast<long double>(VOL) * basis[row] * basis[col];
        weighted_sum[row][col].Add(scaled_weight * f);
        block_sum[block][row][col].Add(scaled_weight * f);
      }
    }
  }

  if (mode == "fid" || mode == "det") {
    printf("Fiducial flat MC phase space events = %d (acceptance %0.3f "
           "percent) \n",
           static_cast<int>(fiducial),
           fiducial / static_cast<double>(ind.size()) * 100);
  }
  if ((mode == "fid" || mode == "det") && fiducial == 0) {
    throw std::invalid_argument(
        "spherical::GetGMixing: response has zero fiducial events");
  }
  if (mode == "det") {
    if (selected == 0) {
      throw std::invalid_argument(
          "spherical::GetGMixing: response has zero selected events");
    }
    printf("Selected flat MC phase space events = %d (efficiency %0.3f "
           "percent) \n",
           static_cast<int>(selected),
           selected / static_cast<double>(fiducial) * 100);
  }
  if (mode == "fid" || mode == "det") {
    if (fiducial < 50 || (mode == "det" && selected < 50)) {
      std::cout << rang::fg::red
                << "GetGMixing:: Too low MC event count in the mass bin!"
                << rang::fg::reset << std::endl;
    }
  }
  std::cout << std::endl;

  // Normalize with scale-free reference weights
  std::vector<MMatrix<double>> jackknife;
  jackknife.reserve(jackknife_blocks);
  for (std::size_t block = 0; block < jackknife_blocks; ++block) {
    const long double retained_weight =
        generated_weight_scaled - block_weight[block].Value();
    if (retained_weight <= 0.0L) {
      throw std::invalid_argument("spherical::GetGMixing: jackknife replica has no positive reference weight");
    }
    MMatrix<double> replica(NCOEF, NCOEF, 0.0);
    for (int row = 0; row < NCOEF; ++row) {
      for (int col = 0; col < NCOEF; ++col) {
        replica[row][col] =
            static_cast<double>((weighted_sum[row][col].Value() -
                                 block_sum[block][row][col].Value()) /
                                retained_weight);
      }
    }
    jackknife.push_back(std::move(replica));
  }

  for (int row = 0; row < NCOEF; ++row) {
    for (int col = 0; col < NCOEF; ++col) {
      E[row][col] = static_cast<double>(weighted_sum[row][col].Value() /
                                        generated_weight_scaled);
    }
  }

  MMatrix<long double> centered_sum(NCOEF, NCOEF, 0.0L);
  for (const std::size_t k : ind) {
    const double weight = EventWeight(events[k], "spherical::GetGMixing");
    const long double scaled_weight =
        static_cast<long double>(weight) / weight_scale;
    std::vector<double> basis(static_cast<std::size_t>(NCOEF), 0.0);
    if (AcceptEvent(events[k], mode)) {
      basis = Basis(events[k], LMAX);
    }
    for (int row = 0; row < NCOEF; ++row) {
      for (int col = 0; col < NCOEF; ++col) {
        const long double f =
            static_cast<long double>(VOL) * basis[row] * basis[col];
        const long double mean =
            weighted_sum[row][col].Value() / generated_weight_scaled;
        const long double residual = f - mean;
        const long double weighted_residual = scaled_weight * residual;
        centered_sum[row][col] += weighted_residual * weighted_residual;
      }
    }
  }
  std::cout << "Symmetric mixing matrix integral coefficients G_{ll'}^{mm'}:"
            << std::endl;
  std::cout << "         " << std::endl;

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      printf("G_%d%d  \t", l, m);
    }
  }
  std::cout << std::endl;

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = LinearInd(l, m);

      printf("G_%d%d : ", l, m);
      double colsum = 0;

      for (int lprime = 0; lprime <= LMAX; ++lprime) {
        for (int mprime = -lprime; mprime <= lprime; ++mprime) {
          const int indexprime = LinearInd(lprime, mprime);

          colsum += E[index][indexprime];

          const double error = std::sqrt(statistics::WeightedMeanVariance(
              generated_weight_scaled, reference_weights.ScaledSquareSum(), centered_sum[index][indexprime]));
          printf("%6.3f +- %6.3f \t", E[index][indexprime], error);
        }
      }
      printf(" | %6.3f \n", colsum);
    }
  }

  gra::aux::PrintBar("-");
  std::cout << "       " << std::endl;
  for (std::size_t j = 0; j < E.size_row(); ++j) {
    double rowsum = 0;
    for (std::size_t i = 0; i < E.size_col(); ++i) {
      rowsum += E[i][j];
    }
    printf("%6.3f\t", rowsum);
  }
  printf("\n\n\n\n");

  // Calculate matrix condition number via SVD
  const std::vector<double> singular_values = E.SingularValues();

  std::cout << "SVD singular values of the matrix:" << std::endl;
  for (const double value : singular_values) {
    printf("%0.4f ", value);
  }
  std::cout << std::endl;
  const double conditionnumber =
      singular_values.front() / singular_values.back();
  if (conditionnumber < 10) {
    std::cout << rang::fg::green;
  } else {
    std::cout << rang::fg::red;
  }

  printf("Condition number: sMax/sMin = %0.4f \n", conditionnumber);
  std::cout << std::endl << std::endl << std::endl;
  std::cout << rang::fg::reset;
  return {E, jackknife, static_cast<double>(reference_weights.Sum())};
}

// Acceptance expansion coefficients using MC events
std::pair<std::vector<double>, std::vector<double>>
GetELM(const std::vector<Omega> &MC, const std::vector<std::size_t> &ind,
       int LMAX, const std::string &mode) {
  if (mode != "fla" && mode != "fid" && mode != "det") {
    throw std::invalid_argument("spherical::GetELM: Unknown mode " + mode);
  }
  if (ind.empty()) {
    throw std::invalid_argument(
        "spherical::GetELM: generated event cell is empty");
  }
  for (const std::size_t k : ind) {
    if (k >= MC.size()) {
      throw std::out_of_range(
          "spherical::GetELM: event index outside input array");
    }
  }
  const int NCOEF = BasisSize(LMAX);

  std::cout << "GetELM: mode = " << mode << std::endl;
  std::cout << "Generated flat MC phase space events = " << ind.size()
            << std::endl;

  // Construct the efficiency coefficients E_LM with linear indexing
  std::vector<double> E(NCOEF, 0.0);
  std::vector<statistics::CompensatedSum> weighted_sum(
      static_cast<std::size_t>(NCOEF));

  // Loop over GENERATED MC events and do the integral in effect via
  // uniform MC sampling
  // We evaluate the integral:
  //
  // E_LM = \int E(Omega) Re Y_LM(Omega) dOmega, Omega = (costheta,phi)

  std::size_t fiducial = 0;
  std::size_t selected = 0;
  statistics::ScaledWeightSums reference_weights;
  for (const std::size_t k : ind) {
    reference_weights.Add(EventWeight(MC[k], "spherical::GetELM"));
  }
  ValidateReferenceWeights(reference_weights, "spherical::GetELM");
  const long double weight_scale = reference_weights.Scale();
  const long double generated_weight_scaled = reference_weights.ScaledSum();
  const double V = sqrt(4.0 * PI); // Normalization volume

  for (const auto &k : ind) {
    const double weight = EventWeight(MC[k], "spherical::GetELM");
    const long double scaled_weight =
        static_cast<long double>(weight) / weight_scale;
    bool fid = MC[k].fiducial;
    bool sel = MC[k].selected;

    // Flat flat phase space
    if (mode == "fla") {
      // all fine
    }
    // Geometric acceptance ("Fiducial phase space")
    else if (mode == "fid") {
      if (!fid) {
        continue;
      } else {
        ++fiducial;
      }
    }
    // Geometric x Efficiency ("Detector level")
    else if (mode == "det") {
      if (!fid) {
        continue;
      } else {
        ++fiducial;
      }
      if (!sel) {
        continue;
      } else {
        ++selected;
      }
    }

    const auto basis = Basis(MC[k], LMAX);
    for (const auto &index : indices(basis)) {
      weighted_sum[index].Add(scaled_weight * V * basis[index]);
    }
  }
  if (mode == "fid" || mode == "det") {
    printf("Fiducial flat MC phase space events = %d (geometric-kinematic "
           "acceptance "
           "%0.3f percent) \n\n",
           static_cast<int>(fiducial),
           fiducial / static_cast<double>(ind.size()) * 100);
  }
  if ((mode == "fid" || mode == "det") && fiducial == 0) {
    throw std::invalid_argument(
        "spherical::GetELM: response has zero fiducial events");
  }
  if (mode == "det") {
    if (selected == 0) {
      throw std::invalid_argument(
          "spherical::GetELM: response has zero selected events");
    }
    printf("Selected flat MC phase space events = %d (fiducial efficiency "
           "%0.3f percent) "
           "\n\n",
           static_cast<int>(selected),
           selected / static_cast<double>(fiducial) * 100);
  }
  if (mode == "fid") {
    if (fiducial < 1000) {
      std::cout << "GetELM:: Very low [fiducial] MC event count = " << fiducial
                << " in the mass bin!" << std::endl;
    }
  }
  if (mode == "det") {
    if (selected < 1000) {
      std::cout << "GetELM:: Very low [selected] MC event count = " << selected
                << " in the mass bin!" << std::endl;
    }
  }

  // Normalize with scale-free reference weights
  std::vector<double> E_error(E.size(), 0.0);
  for (int index = 0; index < NCOEF; ++index) {
    E[index] = static_cast<double>(weighted_sum[index].Value() /
                                   generated_weight_scaled);
  }

  std::vector<long double> centered_sum(E.size(), 0.0L);
  for (const std::size_t k : ind) {
    const double weight = EventWeight(MC[k], "spherical::GetELM");
    const long double scaled_weight =
        static_cast<long double>(weight) / weight_scale;
    std::vector<double> basis(static_cast<std::size_t>(NCOEF), 0.0);
    if (AcceptEvent(MC[k], mode)) {
      basis = Basis(MC[k], LMAX);
      gra::Scale(basis, V);
    }
    for (int index = 0; index < NCOEF; ++index) {
      const long double mean =
          weighted_sum[index].Value() / generated_weight_scaled;
      const long double residual =
          static_cast<long double>(basis[index]) - mean;
      const long double weighted_residual = scaled_weight * residual;
      centered_sum[index] += weighted_residual * weighted_residual;
    }
  }

  std::cout << "Acceptance decomposition coefficients:" << std::endl;
  double sum = 0;
  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = LinearInd(l, m);

      // 1 sigma MC integration uncertainty
      const double error = std::sqrt(statistics::WeightedMeanVariance(
          generated_weight_scaled, reference_weights.ScaledSquareSum(), centered_sum[index]));
      E_error[index] = error;
      sum += E[index];

      printf("E_%d%d \t= %12.8f +- %12.8f   (rel.error %9.3f percent) \n", l, m,
             E[index], error, std::abs(error / E[index]) * 100);
    }
  }
  printf("SUM_lm =   %0.8f \n", sum);
  std::cout << "Uncertanties by Monte Carlo errors = sqrt{(<f^2> - <f>^2) / n}"
            << std::endl;
  std::cout << std::endl;

  return {E, E_error};
}

// Calculate expansion coefficients directly via algebraic expansion
//
MomentEstimate SphericalMoments(const std::vector<Omega> &input,
                                const std::vector<std::size_t> &ind, int LMAX,
                                const std::string &mode) {
  if (mode != "fla" && mode != "fid" && mode != "det") {
    throw std::invalid_argument("spherical::SphericalMoments: Unknown mode " +
                                mode);
  }
  const double V = sqrt(4.0 * PI); // Normalization volume
  const std::size_t NCOEF = static_cast<std::size_t>(BasisSize(LMAX));
  statistics::WeightedVectorSums sums(NCOEF);
  for (const std::size_t k : ind) {
    if (k >= input.size()) {
      throw std::out_of_range("spherical::SphericalMoments: event index outside input array");
    }
    if (!AcceptEvent(input[k], mode)) { continue; }
    auto basis = Basis(input[k], LMAX);
    gra::Scale(basis, V);
    sums.Add(basis, EventWeight(input[k], "spherical::SphericalMoments"));
  }
  return {sums.Sum(), sums.Covariance(), static_cast<double>(sums.Weights().Sum()),
          static_cast<double>(sums.Weights().SquareSum()), sums.Weights().Count()};
}

// Calculate normalized real spherical harmonics for each event
// Y_(i,lm) = Y_lm(Omega_i)
MMatrix<double> YLM(const std::vector<Omega> &events, int LMAX) {
  MMatrix<double> Y_lm(events.size(), BasisSize(LMAX), 0.0);
  for (const auto &k : indices(events)) {
    const auto basis = Basis(events[k], LMAX);
    std::copy(basis.begin(), basis.end(), Y_lm.Row(k).begin());
  }
  return Y_lm;
}

// Spherical harmonic dot product between coefficients gives the integral
// This is the property of orthonormality
// \int G(Omega) x(omega) dOmega <=> \sum_i G_i x_i
double HarmDotProd(const std::vector<double> &G, const std::vector<double> &x,
                   const std::vector<bool> &ACTIVE, int LMAX) {
  const std::size_t ncoef = static_cast<std::size_t>(BasisSize(LMAX));
  if (G.size() != ncoef || x.size() != ncoef || ACTIVE.size() != ncoef) {
    throw std::invalid_argument(
        "spherical::HarmDotProd: incompatible dimensions");
  }
  return gra::MaskedBilinearProduct(G, x, ACTIVE);
}

// Propagate a full moment covariance through a harmonic dot product
// sigma = sqrt(G_active^T Cov G_active)
double HarmDotProdError(const std::vector<double> &G,
                        const MMatrix<double> &covariance,
                        const std::vector<bool> &ACTIVE, int LMAX) {
  const std::size_t ncoef = static_cast<std::size_t>(BasisSize(LMAX));
  if (G.size() != ncoef || ACTIVE.size() != ncoef ||
      covariance.size_row() != ncoef || covariance.size_col() != ncoef) {
    throw std::invalid_argument(
        "spherical::HarmDotProdError: incompatible dimensions");
  }
  const auto matrix = covariance.Transform([](double value) { return static_cast<long double>(value); });
  const std::vector<long double> g(G.begin(), G.end());
  auto magnitude = g;
  std::transform(g.begin(), g.end(), magnitude.begin(), [](long double value) { return std::abs(value); });
  const long double variance = matrix.MaskedBilinearForm(g, g, ACTIVE, ACTIVE);
  const long double scale = matrix.Transform([](long double value) { return std::abs(value); })
                               .MaskedBilinearForm(magnitude, magnitude, ACTIVE, ACTIVE);
  const long double tolerance = 64.0L * std::numeric_limits<double>::epsilon() * scale;
  if (!std::isfinite(variance) || variance < -tolerance) {
    throw std::domain_error(
        "spherical::HarmDotProdError: invalid propagated variance");
  }
  return static_cast<double>(std::sqrt(std::max(0.0L, variance)));
}

// Print active and inactive moments with their standard errors
void PrintOutMoments(const std::vector<double> &x,
                     const std::vector<double> &x_error,
                     const std::vector<bool> &ACTIVE, int LMAX) {
  const auto size = static_cast<std::size_t>(BasisSize(LMAX));
  if (x.size() != size || x_error.size() != size || ACTIVE.size() != size) {
    throw std::invalid_argument("spherical::PrintOutMoments: incompatible dimensions");
  }
  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const unsigned int index = LinearInd(l, m);
      printf("t_%d%d \t= %8.1f +- %5.1f   ", l, m, x[index], x_error[index]);

      if (ACTIVE[index]) {
        std::cout << "[" << rang::fg::green << "active" << rang::fg::reset
                  << "]" << std::endl;
      } else {
        std::cout << "[" << rang::fg::red << "inactive" << rang::fg::reset
                  << "]" << std::endl;
      }
    }
  }
  std::cout << std::endl;
}

// Calculate indices for this interval
//
std::vector<std::size_t> GetIndices(const std::vector<Omega> &events,
                                    const std::vector<double> &M,
                                    const std::vector<double> &Pt,
                                    const std::vector<double> &Y,
                                    const std::vector<bool> &include_upper) {
  if (M.size() != 2 || Pt.size() != 2 || Y.size() != 2 ||
      include_upper.size() != 3 || !(M[0] < M[1]) || !(Pt[0] < Pt[1]) ||
      !(Y[0] < Y[1])) {
    throw std::invalid_argument(
        "spherical::GetIndices: invalid interval definition");
  }
  std::vector<std::size_t> ind;

  for (const auto &i : indices(events)) {
    const bool in_mass =
        events[i].M >= M[0] &&
        (events[i].M < M[1] || (include_upper[0] && events[i].M <= M[1]));
    const bool in_pt =
        events[i].Pt >= Pt[0] &&
        (events[i].Pt < Pt[1] || (include_upper[1] && events[i].Pt <= Pt[1]));
    const bool in_y =
        events[i].Y >= Y[0] &&
        (events[i].Y < Y[1] || (include_upper[2] && events[i].Y <= Y[1]));
    if (in_mass && in_pt && in_y) {
      ind.push_back(i);
    }
  }

  gra::aux::PrintBar("-");
  std::cout << rang::fg::green;
  printf("MASS RANGE: [%0.3f, %0.3f] GeV, PT RANGE: [%0.3f, %0.3f] GeV, Y "
         "RANGE: [%0.3f, %0.3f] "
         ": Events in this hyperbin %lu/%lu \n\n",
         M[0], M[1], Pt[0], Pt[1], Y[0], Y[1], ind.size(), events.size());
  std::cout << rang::fg::reset;

  return ind;
}

// Summarize nominal weights after one analysis-level selection
WeightSummary SummarizeWeights(const std::vector<Omega> &events,
                               const std::vector<std::size_t> &ind,
                               const std::string &mode) {
  if (mode != "fla" && mode != "fid" && mode != "det") {
    throw std::invalid_argument("spherical::SummarizeWeights: unknown mode " + mode);
  }
  WeightSummary summary;
  for (const std::size_t index : ind) {
    if (index >= events.size()) {
      throw std::out_of_range(
          "spherical::SummarizeWeights: event index outside input array");
    }
    if (!AcceptEvent(events[index], mode)) {
      continue;
    }
    const double weight =
        EventWeight(events[index], "spherical::SummarizeWeights");
    summary.Add(weight);
  }
  return summary;
}

// Build a full-size pseudoinverse with inactive coefficient rows fixed to zero
MMatrix<double> ActiveResponsePseudoInverse(const MMatrix<double> &response,
                                            const std::vector<bool> &active,
                                            double svd_regularization) {
  if (response.size_row() != response.size_col() ||
      response.size_col() != active.size()) {
    throw std::invalid_argument(
        "spherical::ActiveResponsePseudoInverse: incompatible dimensions");
  }
  std::vector<std::size_t> active_index;
  for (const auto &index : indices(active)) {
    if (active[index]) {
      active_index.push_back(index);
    }
  }
  if (active_index.empty()) {
    throw std::invalid_argument(
        "spherical::ActiveResponsePseudoInverse: no active coefficients");
  }

  const MMatrix<double> reduced = response.SelectColumns(active_index);
  const double relative_cut =
      svd_regularization *
      static_cast<double>(std::max(reduced.size_row(), reduced.size_col()));
  const MMatrix<double> reduced_inverse = reduced.PseudoInverse(relative_cut);
  MMatrix<double> full_inverse(response.size_col(), response.size_row(), 0.0);
  full_inverse.SetRows(active_index, reduced_inverse);
  return full_inverse;
}

// Compute the linear index of one spherical-harmonic mode
// i = l(l+1)+m
int LinearInd(int l, int m) { return l * (l + 1) + m; }

// Compute the standard error from first and second moments
// sigma_mean = sqrt((<f^2>-<f>^2)/N)
double CalcError(double f2, double f, double N) {
  if (!std::isfinite(f2) || !std::isfinite(f) || !std::isfinite(N) ||
      N <= 0.0) {
    throw std::invalid_argument(
        "spherical::CalcError: finite moments and positive count required");
  }
  const long double square = static_cast<long double>(f) * f;
  const long double variance = static_cast<long double>(f2) - square;
  const long double tolerance = 64.0L * std::numeric_limits<double>::epsilon() *
                                std::max(std::abs(static_cast<long double>(f2)), square);
  if (variance < -tolerance) {
    throw std::domain_error(
        "spherical::CalcError: negative variance beyond roundoff");
  }
  return static_cast<double>(std::sqrt(std::max(0.0L, variance) / N));
}

// Print matrix to a file
void PrintMatrix(FILE *fp, const std::vector<std::vector<double>> &A) {
  // Print out coefficients
  for (const auto &i : indices(A)) {
    for (const auto &j : indices(A[i])) {
      fprintf(fp, "%0.1f ", A[i][j]);
    }
    fprintf(fp, "\n");
  }
  fprintf(fp, "\n");
}

// Synthesize distributions with real SH-basis
//
// f(theta,phi) = \sum_{l=0}\sum_{m=-l}^l c_lm x Y_R^{lm}(theta,phi)
//
MMatrix<double> Y_real_synthesize(const std::vector<double> &c_lm,
                                  const std::vector<bool> &ACTIVE,
                                  std::size_t N, std::vector<double> &costheta,
                                  std::vector<double> &phi, bool normalized) {
  if (N < 2 || ACTIVE.empty() || c_lm.size() != ACTIVE.size() || !gra::AllFinite(c_lm)) {
    throw std::invalid_argument(
        "spherical::Y_real_synthesize: invalid dimensions");
  }
  const int LMAX = msqrt(ACTIVE.size()) - 1;
  if (static_cast<std::size_t>(BasisSize(LMAX)) != ACTIVE.size()) {
    throw std::invalid_argument("spherical::Y_real_synthesize: coefficient "
                                "count is not a square basis");
  }

  // Cos(theta) and phi
  costheta = math::linspace(-1.0, 1.0, N);
  phi = math::linspace(-math::PI, math::PI, N);

  // Do the expansion
  MMatrix<double> Z(N, N, 0.0);

  for (int l = 0; l <= LMAX; ++l) {
    for (int m = -l; m <= l; ++m) {
      const int index = gra::spherical::LinearInd(l, m);
      if (!ACTIVE[index]) {
        continue;
      } // Not active

      for (const std::size_t &i : indices(costheta)) {
        for (const std::size_t &j : indices(phi)) {
          // Add value
          Z[i][j] +=
              c_lm[index] * math::Y_real_basis(costheta[i], phi[j], l, m);
        }
      }
    }
  }

  if (normalized) {
    double max_abs = 0.0;
    for (std::size_t i = 0; i < Z.size_row(); ++i) {
      for (std::size_t j = 0; j < Z.size_col(); ++j) {
        max_abs = std::max(max_abs, std::abs(Z[i][j]));
      }
    }
    if (max_abs > 0.0) {
      Z = Z / max_abs;
    }
  }

  return Z;
}

// Propagate independent input errors through one linear transformation
// sigma_y,m^2 = sum_n A_mn^2 sigma_x,n^2
std::vector<double> ErrorProp(const MMatrix<double> &A,
                              const std::vector<double> &x) {
  if (A.size_col() != x.size()) {
    throw std::invalid_argument(
        "spherical::ErrorProp: incompatible dimensions");
  }
  std::vector<double> y(A.size_row(), 0.0);

  for (std::size_t m = 0; m < A.size_row(); ++m) {
    y[m] = std::inner_product(A.Row(m).begin(), A.Row(m).end(), x.begin(), 0.0,
                             [](double a, double b) { return std::hypot(a, b); }, std::multiplies<double>{});
  }
  return y;
}

// Propagate a full covariance matrix through one linear transformation
// Cov_y = A Cov_x A^T
MMatrix<double> CovarianceProp(const MMatrix<double> &A,
                               const MMatrix<double> &covariance) {
  if (covariance.size_row() != covariance.size_col() ||
      A.size_col() != covariance.size_row()) {
    throw std::invalid_argument(
        "spherical::CovarianceProp: incompatible dimensions");
  }
  return A * covariance * A.Transpose();
}

// Compute standard deviations from a validated covariance diagonal
// sigma_i = sqrt(Cov_ii)
std::vector<double> CovarianceErrors(const MMatrix<double> &covariance) {
  if (covariance.size_row() != covariance.size_col()) {
    throw std::invalid_argument(
        "spherical::CovarianceErrors: covariance must be square");
  }
  const auto diagonal = covariance.GetDiag();
  const double scale = diagonal.empty() ? 0.0 : std::abs(*std::max_element(
      diagonal.begin(), diagonal.end(), [](double a, double b) { return std::abs(a) < std::abs(b); }));
  const double tolerance = 64.0 * std::numeric_limits<double>::epsilon() * scale;
  std::vector<double> errors(diagonal.size(), 0.0);
  for (const auto &i : indices(errors)) {
    if (!std::isfinite(covariance[i][i]) || covariance[i][i] < -tolerance) {
      throw std::domain_error(
          "spherical::CovarianceErrors: invalid covariance diagonal");
    }
    errors[i] = std::sqrt(std::max(0.0, covariance[i][i]));
  }
  return errors;
}

// Estimate finite-reference response covariance with delete-group jackknife
// replicas
MMatrix<double> ResponseJackknifeCovariance(
    const std::vector<MMatrix<double>> &inverse_response,
    const std::vector<MMatrix<double>> &forward_response,
    const std::vector<double> &observed, const std::vector<bool> &active,
    double svd_regularization) {
  if (inverse_response.size() < 2) {
    throw std::invalid_argument("spherical::ResponseJackknifeCovariance: at "
                                "least two replicas required");
  }
  if (!forward_response.empty() &&
      forward_response.size() != inverse_response.size()) {
    throw std::invalid_argument(
        "spherical::ResponseJackknifeCovariance: replica counts differ");
  }
  const std::size_t dimension = observed.size();
  if (dimension == 0 || active.size() != dimension) {
    throw std::invalid_argument("spherical::ResponseJackknifeCovariance: "
                                "incompatible vector dimensions");
  }
  for (const auto &replica : indices(inverse_response)) {
    if (inverse_response[replica].size_row() != dimension ||
        inverse_response[replica].size_col() != dimension ||
        (!forward_response.empty() &&
         (forward_response[replica].size_row() != dimension ||
          forward_response[replica].size_col() != dimension))) {
      throw std::invalid_argument("spherical::ResponseJackknifeCovariance: "
                                  "incompatible matrix dimensions");
    }
  }

  std::vector<std::vector<double>> replicas;
  replicas.reserve(inverse_response.size());
  for (const auto &replica : indices(inverse_response)) {
    std::vector<double> value =
        ActiveResponsePseudoInverse(inverse_response[replica], active,
                                    svd_regularization) *
        observed;
    if (!forward_response.empty()) {
      value = forward_response[replica] * value;
    }
    replicas.push_back(std::move(value));
  }

  return statistics::JackknifeCovariance(replicas);
}

} // namespace spherical
} // namespace gra
