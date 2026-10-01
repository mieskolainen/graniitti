// Elastic Coulomb and nuclear interference runtime
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <complex>
#include <cstddef>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <numbers>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Eikonal/MElasticCNI.h"
#include "Graniitti/Eikonal/MEikonalCache.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "json.hpp"

using gra::aux::indices;

namespace gra {

namespace {

using Complex = std::complex<double>;
using Matrix = MElasticCNI::Matrix;
using math::CompositeWeight;
using math::ParallelFor;
constexpr std::size_t kHelicityDimension = 4;
constexpr std::size_t kHelicityEntries = 16;
constexpr std::size_t kCacheVersion = 1;
constexpr std::size_t kMaximumTailSeriesOrder = 32;
constexpr double kPointIRScale = 1.0;

// Retain the largest normalized interpolation error and its coordinates
struct WorstInterpolation {
  double ratio = 0.0;
  double error = 0.0;
  double allowance = 0.0;
  double q2 = 0.0;
  std::size_t entry = 0;

  // Update the retained error when one validation point is worse
  void Update(const double new_error, const double new_allowance,
              const double new_q2, const std::size_t new_entry = 0) {
    const double new_ratio = new_error / new_allowance;
    if (!std::isfinite(new_ratio) || new_ratio > ratio) {
      ratio = new_ratio;
      error = new_error;
      allowance = new_allowance;
      q2 = new_q2;
      entry = new_entry;
    }
  }
};

// Build the absolute harmonics of all row-major proton spin transitions
constexpr std::array<int, kHelicityEntries> BuildAbsoluteHarmonicTable() {
  std::array<int, kHelicityEntries> harmonic{};
  constexpr auto transition = CanonicalProtonHelicityTransitions();
  for (const auto &entry : indices(harmonic)) {
    const int signed_harmonic = transition[entry].azimuth_harmonic;
    harmonic[entry] = signed_harmonic < 0 ? -signed_harmonic : signed_harmonic;
  }
  return harmonic;
}

constexpr auto kAbsoluteHarmonic = BuildAbsoluteHarmonicTable();

// Append one floating-point value in an exact cache-key representation
void AppendCacheDouble(std::ostringstream &stream, const double value) {
  stream << std::hexfloat << value << ';';
}

// Compute the absolute transfer harmonic for one row-major matrix entry
int AbsoluteHarmonic(const std::size_t index) {
  return kAbsoluteHarmonic.at(index);
}

// Convert one row-major helicity array to the reusable matrix container
Matrix HelicityArrayToMatrix(const ProtonHelicityMatrix &array) {
  Matrix matrix(kHelicityDimension);
  for (std::size_t row = 0; row < kHelicityDimension; ++row) {
    for (std::size_t col = 0; col < kHelicityDimension; ++col) {
      matrix(row, col) = array[kHelicityDimension * row + col];
    }
  }
  return matrix;
}

// Convert one four-dimensional matrix to a row-major helicity array
ProtonHelicityMatrix MatrixToHelicityArray(const Matrix &matrix) {
  if (matrix.size_row() != kHelicityDimension ||
      matrix.size_col() != kHelicityDimension) {
    throw std::invalid_argument(
        "MatrixToHelicityArray: expected a four-dimensional matrix");
  }
  ProtonHelicityMatrix array{};
  for (std::size_t row = 0; row < kHelicityDimension; ++row) {
    for (std::size_t col = 0; col < kHelicityDimension; ++col) {
      array[kHelicityDimension * row + col] = matrix(row, col);
    }
  }
  return array;
}

// Compute the exact two-body center-of-mass momentum
// p_cm = sqrt([s-(m1+m2)^2][s-(m1-m2)^2])/(2 sqrt(s))
double CenterOfMassMomentum(const double s, const double mass1,
                            const double mass2) {
  const double root_s = std::sqrt(s);
  const double first = s - math::pow2(mass1 + mass2);
  const double second = s - math::pow2(mass1 - mass2);
  if (!(first > 0.0) || !(second > 0.0)) {
    throw std::invalid_argument(
        "MElasticCNI: center-of-mass energy is below threshold");
  }
  return std::sqrt(first * second) / (2.0 * root_s);
}

// Construct exact elastic center-of-mass momenta in the reference plane
std::array<M4Vec, 4> ReferenceMomenta(const double s, const double mass1,
                                      const double mass2, const double abs_t) {
  const double root_s = std::sqrt(s);
  const double momentum = CenterOfMassMomentum(s, mass1, mass2);
  const double maximum_abs_t = 4.0 * momentum * momentum;
  if (!std::isfinite(abs_t) || abs_t < 0.0 ||
      abs_t > maximum_abs_t * (1.0 + 1.0e-12)) {
    throw std::invalid_argument(
        "MElasticCNI: momentum transfer is outside two-body phase space");
  }
  const double energy1 = (s + mass1 * mass1 - mass2 * mass2) / (2.0 * root_s);
  const double energy2 = (s + mass2 * mass2 - mass1 * mass1) / (2.0 * root_s);
  const double q = std::sqrt(abs_t);
  const double transverse =
      q * std::sqrt(std::max(0.0, 1.0 - abs_t / (4.0 * momentum * momentum)));
  const double longitudinal = momentum - abs_t / (2.0 * momentum);
  return {M4Vec(0.0, 0.0, momentum, energy1),
          M4Vec(0.0, 0.0, -momentum, energy2),
          M4Vec(transverse, 0.0, longitudinal, energy1),
          M4Vec(-transverse, 0.0, -longitudinal, energy2)};
}

// Hold the short-range Born table and signed single-flip pole coefficients
struct ResidualBornData {
  std::vector<ProtonHelicityMatrix> value;
  ProtonHelicityMatrix single_flip_coefficient{};
  ProtonHelicityMatrix double_flip_coefficient{};
  Matrix single_flip_eikonal_coefficient;
  Matrix double_flip_eikonal_coefficient;
};

// Extract signed small-q flip coefficients from the exact current
std::pair<ProtonHelicityMatrix, ProtonHelicityMatrix>
ExtractFlipCoefficients(const MDirac &dirac,
                        const std::vector<MParticle> &initialstate,
                        const double s, const double maximum_q,
                        const form::ParamStore &structure) {
  const double q_reference = std::min(1.0e-3, maximum_q * 1.0e-3);
  if (!(q_reference > 0.0)) {
    throw std::invalid_argument(
        "ExtractFlipCoefficients: invalid momentum range");
  }
  const std::array<double, 3> q = {q_reference, 0.5 * q_reference,
                                   0.25 * q_reference};
  std::array<ProtonHelicityMatrix, 3> amplitude;
  for (const auto &sample : indices(q)) {
    const auto momentum = ReferenceMomenta(
        s, initialstate[0].mass, initialstate[1].mass, q[sample] * q[sample]);
    amplitude[sample] = qed::ElasticSpinHalfPhotonExchange(
        dirac, initialstate[0], initialstate[1], momentum[0], momentum[1],
        momentum[2], momentum[3], structure);
  }

  ProtonHelicityMatrix coarse{};
  ProtonHelicityMatrix single_flip{};
  ProtonHelicityMatrix double_flip{};
  double coefficient_scale = 0.0;
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    const int harmonic = AbsoluteHarmonic(entry);
    if (harmonic != 1 && harmonic != 2) {
      continue;
    }
    const double first_scale = harmonic == 1 ? q[0] : 1.0;
    const double second_scale = harmonic == 1 ? q[1] : 1.0;
    const double third_scale = harmonic == 1 ? q[2] : 1.0;
    coarse[entry] = (4.0 * second_scale * amplitude[1][entry] -
                     first_scale * amplitude[0][entry]) /
                    3.0;
    const Complex coefficient = (4.0 * third_scale * amplitude[2][entry] -
                                 second_scale * amplitude[1][entry]) /
                                3.0;
    if (harmonic == 1) {
      single_flip[entry] = coefficient;
    } else {
      double_flip[entry] = coefficient;
    }
    coefficient_scale = std::max(coefficient_scale, std::abs(coefficient));
  }
  const double absolute_floor = std::max(1.0, coefficient_scale) * 64.0 *
                                std::numeric_limits<double>::epsilon();
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    const int harmonic = AbsoluteHarmonic(entry);
    if (harmonic != 1 && harmonic != 2) {
      continue;
    }
    const Complex coefficient =
        harmonic == 1 ? single_flip[entry] : double_flip[entry];
    const double relative_error =
        std::abs(coefficient - coarse[entry]) /
        std::max({std::abs(coefficient), absolute_floor});
    if (!std::isfinite(relative_error) || relative_error > 1.0e-6) {
      throw std::runtime_error(
          "MElasticCNI: unstable exact flip-tail extraction at entry " +
          std::to_string(entry));
    }
  }
  return {single_flip, double_flip};
}

// Compute the point-subtracted electromagnetic Born table and its flip pole
ResidualBornData BuildResidualBornTable(
    const MDirac &dirac, const std::vector<MParticle> &initialstate,
    const double s, const double sqrt_lambda, const double maximum_q,
    const std::size_t intervals, const double point_born_coefficient,
    const form::ParamStore &structure) {
  const double q_step = maximum_q / intervals;
  ResidualBornData residual;
  residual.value.resize(intervals + 1);
  std::tie(residual.single_flip_coefficient, residual.double_flip_coefficient) =
      ExtractFlipCoefficients(dirac, initialstate, s, maximum_q, structure);
  ProtonHelicityMatrix single_flip_eikonal{};
  ProtonHelicityMatrix double_flip_eikonal{};
  // A_1 -> D_1/q gives chi_1 -> i D_1/[4 pi sqrt(lambda) b]
  // A_2 -> D_2 gives chi_2 -> -D_2/[2 pi sqrt(lambda) b^2]
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    if (AbsoluteHarmonic(entry) == 1) {
      single_flip_eikonal[entry] = math::zi *
                                   residual.single_flip_coefficient[entry] /
                                   (4.0 * math::PI * sqrt_lambda);
    } else if (AbsoluteHarmonic(entry) == 2) {
      double_flip_eikonal[entry] =
          -residual.double_flip_coefficient[entry] / (2.0 * math::PI * sqrt_lambda);
    }
  }
  residual.single_flip_eikonal_coefficient =
      HelicityArrayToMatrix(single_flip_eikonal);
  residual.double_flip_eikonal_coefficient =
      HelicityArrayToMatrix(double_flip_eikonal);
  ParallelFor(
      intervals + 1,
      [&](const std::size_t index) {
        const double q = index * q_step;
        const double q2 = q * q;
        if (index == 0) {
          residual.value[index] = ProtonHelicityMatrix{};
          return;
        }
        const auto momentum =
            ReferenceMomenta(s, initialstate[0].mass, initialstate[1].mass, q2);
        residual.value[index] = qed::ElasticSpinHalfPhotonExchange(
            dirac, initialstate[0], initialstate[1], momentum[0], momentum[1],
            momentum[2], momentum[3], structure);
        const double point_born = point_born_coefficient / q2;
        for (std::size_t diagonal = 0; diagonal < kHelicityDimension;
             ++diagonal) {
          residual.value[index][kHelicityDimension * diagonal + diagonal] -=
              point_born;
        }
      },
      "Electromagnetic Born table");
  return residual;
}

// Preweight the residual Born table by its radial measure and inverse phases
std::vector<ProtonHelicityMatrix>
WeightedResidualBornTable(const std::vector<ProtonHelicityMatrix> &residual,
                          const double maximum_q, const std::size_t intervals,
                          const double sqrt_lambda) {
  if (residual.size() != intervals + 1 || !std::isfinite(maximum_q) ||
      maximum_q <= 0.0 || intervals < 2 || !std::isfinite(sqrt_lambda) || sqrt_lambda <= 0.0) {
    throw std::invalid_argument(
        "WeightedResidualBornTable: invalid residual transform input");
  }
  const double q_step = maximum_q / intervals;
  const std::array<Complex, 3> hankel_phase = {Complex(1.0), math::zi,
                                               Complex(-1.0)};
  std::vector<ProtonHelicityMatrix> weighted(residual.size());
  for (std::size_t index = 0; index <= intervals; ++index) {
    if (!gra::AllFinite(residual[index])) {
      throw std::runtime_error(
          "WeightedResidualBornTable: non-finite residual Born matrix");
    }
    const double q = index * q_step;
    const double measure = CompositeWeight(index, intervals) * q_step * q /
                           (4.0 * math::PI * sqrt_lambda);
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      weighted[index][entry] = residual[index][entry] * measure *
                               hankel_phase[AbsoluteHarmonic(entry)];
    }
  }
  return weighted;
}

// Inverse transform one preweighted residual electromagnetic Born matrix
Matrix
ResidualEikonalAtB(const std::vector<ProtonHelicityMatrix> &weighted_residual,
                   const double maximum_q, const std::size_t intervals,
                   const double b) {
  if (weighted_residual.size() != intervals + 1 || !std::isfinite(b) ||
      b < 0.0) {
    throw std::invalid_argument(
        "ResidualEikonalAtB: invalid weighted transform input");
  }
  ProtonHelicityMatrix integral{};
  const double q_step = maximum_q / intervals;
  for (std::size_t index = 0; index <= intervals; ++index) {
    const double q = index * q_step;
    const std::array<double, 3> bessel = math::BesselJ012(b * q);
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      integral[entry] +=
          weighted_residual[index][entry] * bessel[kAbsoluteHarmonic[entry]];
    }
  }
  return HelicityArrayToMatrix(integral);
}

// Compute the leading spin-flip residual eikonal tail
// chi_flip(b) = C1/b + C2/b^2
Matrix LinearFlipEikonalTail(const ResidualBornData &residual, const double b) {
  return residual.single_flip_eikonal_coefficient * (1.0 / b) +
         residual.double_flip_eikonal_coefficient * (1.0 / (b * b));
}

// Build the ordered Laurent series of the nonlinear asymptotic eikonal
std::vector<ProtonHelicityMatrix>
BuildNonlinearTailPolynomial(const ResidualBornData &residual,
                             const double boundary_b) {
  if (!std::isfinite(boundary_b) || boundary_b <= 0.0) {
    throw std::invalid_argument(
        "BuildNonlinearTailPolynomial: invalid impact-parameter boundary");
  }
  const Matrix boundary_chi = LinearFlipEikonalTail(residual, boundary_b);
  if (!std::isfinite(boundary_chi.FrobNorm()) ||
      boundary_chi.FrobNorm() > 0.05) {
    throw std::runtime_error(
        "MElasticCNI: asymptotic spin eikonal is not perturbative at the "
        "hadronic impact-parameter boundary");
  }

  std::vector<Matrix> base(3, Matrix(kHelicityDimension));
  base[1] = residual.single_flip_eikonal_coefficient;
  base[2] = residual.double_flip_eikonal_coefficient;
  std::vector<Matrix> power = base;
  std::vector<Matrix> coefficient(3, Matrix(kHelicityDimension));
  Matrix boundary_sum(kHelicityDimension);
  Complex series_coefficient = 0.5 * math::zi;
  bool converged = false;
  for (std::size_t order = 2; order <= kMaximumTailSeriesOrder; ++order) {
    std::vector<Matrix> next(power.size() + 2, Matrix(kHelicityDimension));
    for (const auto &exponent : indices(power)) {
      next[exponent + 1] +=
          power[exponent] * residual.single_flip_eikonal_coefficient;
      next[exponent + 2] +=
          power[exponent] * residual.double_flip_eikonal_coefficient;
    }
    power = std::move(next);
    coefficient.resize(power.size(), Matrix(kHelicityDimension));
    Matrix boundary_term(kHelicityDimension);
    for (std::size_t exponent = order; exponent <= 2 * order; ++exponent) {
      coefficient[exponent].AddScaled(power[exponent], series_coefficient);
      boundary_term.AddScaled(
          power[exponent],
          series_coefficient *
              std::pow(boundary_b, -static_cast<int>(exponent)));
    }
    boundary_sum += boundary_term;
    const double convergence_scale = std::max(1.0, boundary_sum.FrobNorm());
    if (boundary_term.FrobNorm() <=
        16.0 * std::numeric_limits<double>::epsilon() * convergence_scale) {
      converged = true;
      break;
    }
    series_coefficient *= math::zi / static_cast<double>(order + 1);
  }
  if (!converged) {
    throw std::runtime_error(
        "MElasticCNI: nonlinear asymptotic eikonal series did not converge");
  }

  std::vector<ProtonHelicityMatrix> polynomial(coefficient.size());
  for (std::size_t exponent = 2; exponent < coefficient.size(); ++exponent) {
    polynomial[exponent] = MatrixToHelicityArray(coefficient[exponent]);
  }
  return polynomial;
}

// Build the finite point-subtracted impact-parameter remainder profile
std::vector<ProtonHelicityMatrix> BuildFiniteRemainderProfile(
    const MEikonalMatrix &strong_runtime, const ResidualBornData &residual,
    const double maximum_q, const std::size_t intervals, const double sqrt_lambda,
    const double point_ir_scale, const double eta) {
  const auto &b_node = strong_runtime.ImpactParameterNodes();
  const auto weighted_residual = WeightedResidualBornTable(residual.value, maximum_q, intervals, sqrt_lambda);
  std::vector<ProtonHelicityMatrix> profile(b_node.size());
  ParallelFor(
      b_node.size(),
      [&](const std::size_t index) {
        const double b = b_node[index];
        const Matrix chi_residual =
            ResidualEikonalAtB(weighted_residual, maximum_q, intervals, b);
        const Matrix chi_asymptotic = LinearFlipEikonalTail(residual, b);
        const Matrix strong_s = strong_runtime.PhysicalImpactSMatrix(b);
        const double chi_point = 2.0 * eta *
                                 (std::log(point_ir_scale * b / 2.0) +
                                  std::numbers::egamma_v<double>);
        Matrix finite = MElasticCNI::FiniteRemainderProfile(
            strong_s, chi_residual, chi_point);
        // Remove only the linear asymptotic profile before its exact all-b
        // addback
        finite.AddScaled(chi_asymptotic, 1.0 - std::exp(math::zi * chi_point));
        profile[index] = MatrixToHelicityArray(finite);
      },
      "CNI impact-parameter profile");
  return profile;
}

// Forward transform one finite impact-parameter profile at one transfer
ProtonHelicityMatrix
ForwardTransformProfile(const std::vector<ProtonHelicityMatrix> &profile,
                        const std::vector<double> &b_node,
                        const bool logarithmic_b, const double sqrt_lambda,
                        const double q) {
  if (profile.size() != b_node.size() || b_node.size() < 3) {
    throw std::invalid_argument(
        "ForwardTransformProfile: inconsistent impact-parameter grid");
  }
  const std::size_t intervals = b_node.size() - 1;
  const double coordinate_min =
      logarithmic_b ? std::log(b_node.front()) : b_node.front();
  const double coordinate_max =
      logarithmic_b ? std::log(b_node.back()) : b_node.back();
  const double step = (coordinate_max - coordinate_min) / intervals;
  ProtonHelicityMatrix integral{};
  for (std::size_t index = 0; index <= intervals; ++index) {
    const double b = b_node[index];
    const double radial_measure = logarithmic_b ? b * b : b;
    const double measure =
        CompositeWeight(index, intervals) * step * radial_measure;
    const std::array<double, 3> bessel = math::BesselJ012(b * q);
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      integral[entry] +=
          profile[index][entry] * measure * bessel[AbsoluteHarmonic(entry)];
    }
  }
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    const int harmonic = AbsoluteHarmonic(entry);
    integral[entry] *= 4.0 * math::PI * sqrt_lambda * std::pow(-math::zi, harmonic);
  }
  return integral;
}

// Transform the nonlinear asymptotic profile beyond the hadronic b range
ProtonHelicityMatrix
NonlinearAsymptoticTailAtQ(const std::vector<ProtonHelicityMatrix> &polynomial,
                           const double boundary_b, const double sqrt_lambda,
                           const double eta, const double point_ir_scale,
                           const double q, const ElasticCNINumerics &numerics) {
  if (polynomial.size() <= 2) {
    return {};
  }
  const double minimum_qb = q * boundary_b;
  if (!(minimum_qb < numerics.tail_max_qb)) {
    throw std::invalid_argument(
        "MElasticCNI: CNI.tail_max_qb is below the output range");
  }
  std::vector<std::array<Complex, 3>> integral(polynomial.size());
  // Gauss Legendre panels retain the node budget while resolving the smooth tail accurately
  const auto [node, weight] = math::GaussLegendreRule(numerics.tail_order, 0.0, 1.0);
  // Below the configured x=qb split, dlog(x) resolves the leading J0(x)/x tail
  const auto integrate = [&](const double lower, const double upper, const std::size_t intervals,
                             const bool logarithmic) {
    const std::size_t panels = (intervals + node.size() - 1) / node.size();
    const double      span   = logarithmic ? std::log1p((upper - lower) / lower) : upper - lower;
    const double      step   = span / static_cast<double>(panels);
    for (std::size_t panel = 0; panel < panels; ++panel) {
      for (const auto &index : indices(node)) {
        const double  coordinate = (panel + node[index]) * step;
        const double  x          = logarithmic ? lower * std::exp(coordinate) : lower + coordinate;
        const double  b          = x / q;
        const double  chi_point  = 2.0 * eta * (std::log(point_ir_scale * b / 2.0) + std::numbers::egamma_v<double>);
        const Complex point_s    = std::exp(math::zi * chi_point);
        // Integrate exp(i chi_point) [i(I-exp(i chi_asym))-chi_asym] in x=q b
        const double                measure         = weight[index] * step * x * (logarithmic ? x : 1.0) / (q * q);
        const std::array<double, 3> bessel          = math::BesselJ012(x);
        double                      inverse_b_power = 1.0 / b;
        for (std::size_t exponent = 2; exponent < polynomial.size(); ++exponent) {
          inverse_b_power /= b;
          const Complex radial = point_s * measure * inverse_b_power;
          for (const auto &harmonic : indices(bessel)) { integral[exponent][harmonic] += radial * bessel[harmonic]; }
        }
      }
    }
  };
  const double split = std::min(numerics.tail_split_qb, numerics.tail_max_qb);
  if (minimum_qb < split) {
    const double span = std::log1p((split - minimum_qb) / minimum_qb);
    const std::size_t intervals =
        std::max<std::size_t>(1, static_cast<std::size_t>(std::ceil(numerics.TailIntegralN * span /
                                                                  (span + numerics.tail_max_qb - split))));
    integrate(minimum_qb, split, intervals, true);
  }
  if (numerics.tail_max_qb > std::max(minimum_qb, split)) {
    integrate(std::max(minimum_qb, split), numerics.tail_max_qb, numerics.TailIntegralN, false);
  }

  ProtonHelicityMatrix tail{};
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    const int harmonic = AbsoluteHarmonic(entry);
    for (std::size_t exponent = 2; exponent < polynomial.size(); ++exponent) {
      tail[entry] += polynomial[exponent][entry] * integral[exponent][harmonic];
    }
    tail[entry] *= 4.0 * math::PI * sqrt_lambda * std::pow(-math::zi, harmonic);
  }
  return tail;
}

// Compute the finite correction including the analytic all-b flip tails
ProtonHelicityMatrix FiniteCorrectionAtQ(
    const std::vector<ProtonHelicityMatrix> &profile,
    const std::vector<double> &b_node, const bool logarithmic_b, const double sqrt_lambda,
    const ResidualBornData &residual,
    const std::vector<ProtonHelicityMatrix> &nonlinear_polynomial,
    const double eta, const double point_ir_scale, const double q,
    const ElasticCNINumerics &numerics) {
  ProtonHelicityMatrix correction =
      ForwardTransformProfile(profile, b_node, logarithmic_b, sqrt_lambda, q);
  const Complex point_phase =
      MElasticCNI::RenormalizedPointCoulombPhase(eta, point_ir_scale, q * q);
  // The exact n=2 Mellin phase obeys P_2=P_1/(1-i eta)
  const Complex double_flip_phase =
      point_phase / (Complex(1.0) - math::zi * eta);
  const ProtonHelicityMatrix nonlinear = NonlinearAsymptoticTailAtQ(
      nonlinear_polynomial, b_node.back(), sqrt_lambda, eta, point_ir_scale, q, numerics);
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    if (AbsoluteHarmonic(entry) == 1) {
      correction[entry] += residual.single_flip_coefficient[entry] *
                           (point_phase - Complex(1.0)) / q;
    } else if (AbsoluteHarmonic(entry) == 2) {
      correction[entry] += residual.double_flip_coefficient[entry] *
                           (double_flip_phase - Complex(1.0));
    }
    correction[entry] += nonlinear[entry];
  }
  return correction;
}

// Build the logarithmic or linear finite correction momentum table
void BuildFiniteCorrectionTable(
    const std::vector<ProtonHelicityMatrix> &profile,
    const std::vector<double> &b_node, const bool logarithmic_b, const double sqrt_lambda,
    const ResidualBornData &residual,
    const std::vector<ProtonHelicityMatrix> &nonlinear_polynomial,
    const double eta, const double point_ir_scale,
    const ElasticCNINumerics &numerics, const double min_abs_t,
    const double max_abs_t, std::vector<double> &q2_node,
    std::vector<ProtonHelicityMatrix> &correction) {
  const std::size_t points = numerics.NumberKT2 + 1;
  const double coordinate_min =
      numerics.logKT2 ? std::log(min_abs_t) : min_abs_t;
  const double coordinate_max =
      numerics.logKT2 ? std::log(max_abs_t) : max_abs_t;
  const double step = (coordinate_max - coordinate_min) / numerics.NumberKT2;
  q2_node.resize(points);
  correction.resize(points);
  ParallelFor(
      points,
      [&](const std::size_t index) {
        const double coordinate = coordinate_min + index * step;
        const double q2 = numerics.logKT2 ? std::exp(coordinate) : coordinate;
        q2_node[index] = q2;
        correction[index] = FiniteCorrectionAtQ(
            profile, b_node, logarithmic_b, sqrt_lambda, residual,
            nonlinear_polynomial, eta, point_ir_scale, std::sqrt(q2), numerics);
      },
      "CNI momentum correction table");
  q2_node.front() = min_abs_t;
  q2_node.back() = max_abs_t;
}

// Interpolate one helicity matrix table with a local cubic stencil
ProtonHelicityMatrix InterpolateHelicityMatrixTable(
    const std::vector<double> &q2_node,
    const std::vector<ProtonHelicityMatrix> &matrix_table, const double q2) {
  if (!std::isfinite(q2) || q2_node.size() < 4 ||
      matrix_table.size() != q2_node.size()) {
    throw std::invalid_argument(
        "InterpolateHelicityMatrixTable: invalid interpolation table");
  }
  if (q2 <= q2_node.front()) {
    return matrix_table.front();
  }
  if (q2 >= q2_node.back()) {
    return matrix_table.back();
  }
  const auto upper = std::upper_bound(q2_node.cbegin(), q2_node.cend(), q2);
  const std::size_t hi = std::distance(q2_node.cbegin(), upper);
  const std::size_t lo = hi - 1;
  const std::size_t first = lo == 0                    ? 0
                            : hi + 1 >= q2_node.size() ? q2_node.size() - 4
                                                       : lo - 1;
  std::array<double, 4> node{};
  for (const auto &stencil : indices(node)) {
    node[stencil] = q2_node[first + stencil];
  }
  const std::array<double, 4> weight = math::CubicLagrangeWeights(node, q2);
  ProtonHelicityMatrix interpolated{};
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    std::array<Complex, 4> value{};
    for (const auto &stencil : indices(value)) {
      value[stencil] = matrix_table[first + stencil][entry];
    }
    interpolated[entry] = math::CubicLagrangeWeightedSum(value, weight);
  }
  return interpolated;
}

// Compute powers which remove the Coulomb poles of harmonics zero to two
// (q^2,q,1) makes q^(2-|n|) A_n(q) finite for |n| = 0,1,2
std::array<double, 3> PoleRegularizationFactors(const double q) {
  if (!std::isfinite(q) || q <= 0.0) {
    throw std::invalid_argument(
        "PoleRegularizationFactors: q must be finite and positive");
  }
  return {q * q, q, 1.0};
}

// Construct one physical electromagnetic plus higher-order CNI amplitude
ProtonHelicityMatrix
PhysicalCNIAmplitudeAtQ(const MDirac &dirac,
                        const std::vector<MParticle> &initialstate,
                        const double s, const double q2,
                        const ProtonHelicityMatrix &finite_correction,
                        const double point_born_coefficient, const double eta,
                        const double point_ir_scale,
                        const form::ParamStore &structure) {
  if (initialstate.size() != 2 || !std::isfinite(q2) || q2 <= 0.0) {
    throw std::invalid_argument("PhysicalCNIAmplitudeAtQ: invalid state or q2");
  }
  const auto reference =
      ReferenceMomenta(s, initialstate[0].mass, initialstate[1].mass, q2);
  ProtonHelicityMatrix amplitude = qed::ElasticSpinHalfPhotonExchange(
      dirac, initialstate[0], initialstate[1], reference[0], reference[1],
      reference[2], reference[3], structure);
  const Complex point =
      point_born_coefficient / q2 *
      (MElasticCNI::RenormalizedPointCoulombPhase(eta, point_ir_scale, q2) -
       Complex(1.0));
  for (std::size_t diagonal = 0; diagonal < kHelicityDimension; ++diagonal) {
    amplitude[kHelicityDimension * diagonal + diagonal] += point;
  }
  gra::AddScaled(amplitude, finite_correction, Complex(1.0));
  return amplitude;
}

// Remove the Coulomb poles from one complete physical helicity matrix
ProtonHelicityMatrix
RegularizePhysicalHelicityMatrix(const ProtonHelicityMatrix &amplitude,
                                 const double q) {
  const auto pole_factor = PoleRegularizationFactors(q);
  ProtonHelicityMatrix residue = amplitude;
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    residue[entry] *= pole_factor[AbsoluteHarmonic(entry)];
  }
  return residue;
}

// Build the pole-regular total physical strong plus CNI table
void BuildPhysicalTotalResidueTable(
    const MEikonalMatrix &strong_runtime, const MDirac &dirac,
    const std::vector<MParticle> &initialstate, const double s,
    const std::vector<double> &q2_node,
    const std::vector<ProtonHelicityMatrix> &finite_correction,
    const double point_born_coefficient, const double eta,
    const double point_ir_scale, const form::ParamStore &structure,
    std::vector<ProtonHelicityMatrix> &physical_total_residue) {
  if (initialstate.size() != 2 || q2_node.size() < 4 ||
      finite_correction.size() != q2_node.size()) {
    throw std::invalid_argument(
        "BuildPhysicalTotalResidueTable: invalid state or grid");
  }
  physical_total_residue.resize(q2_node.size());
  ParallelFor(q2_node.size(), [&](const std::size_t index) {
    const double q2 = q2_node[index];
    ProtonHelicityMatrix total = strong_runtime.HelicityMatrix(q2, 0, 0);
    const ProtonHelicityMatrix cni = PhysicalCNIAmplitudeAtQ(
        dirac, initialstate, s, q2_node[index], finite_correction[index],
        point_born_coefficient, eta, point_ir_scale, structure);
    gra::AddScaled(total, cni, Complex(1.0));
    physical_total_residue[index] =
        RegularizePhysicalHelicityMatrix(total, std::sqrt(q2));
    if (!gra::AllFinite(physical_total_residue[index])) {
      throw std::runtime_error(
          "BuildPhysicalTotalResidueTable: non-finite physical residue");
    }
  });
}

// Select uniformly spaced and high-curvature amplitude intervals for validation
std::vector<std::size_t> ResidueValidationIntervals(
    const std::vector<double> &q2_node,
    const std::vector<ProtonHelicityMatrix> &physical_total_residue,
    const std::array<double, kHelicityEntries> &entry_scale,
    const double tolerance, const bool logarithmic_q2) {
  const std::size_t intervals = q2_node.size() - 1;
  const std::size_t budget = std::min<std::size_t>(64, intervals);
  const std::size_t uniform_budget = std::min<std::size_t>(32, budget);
  std::vector<bool> selected(intervals, false);
  std::vector<std::size_t> index;
  index.reserve(budget);
  const auto add = [&](const std::size_t value) {
    if (!selected[value]) {
      selected[value] = true;
      index.push_back(value);
    }
  };
  for (std::size_t check = 0; check < uniform_budget; ++check) {
    add(uniform_budget == 1 ? 0
                            : check * (intervals - 1) / (uniform_budget - 1));
  }

  std::vector<std::pair<double, std::size_t>> curvature;
  curvature.reserve(intervals);
  for (std::size_t interval = 0; interval < intervals; ++interval) {
    const double q2 = logarithmic_q2
                          ? std::sqrt(q2_node[interval] * q2_node[interval + 1])
                          : 0.5 * (q2_node[interval] + q2_node[interval + 1]);
    const double fraction =
        (q2 - q2_node[interval]) / (q2_node[interval + 1] - q2_node[interval]);
    const ProtonHelicityMatrix cubic =
        InterpolateHelicityMatrixTable(q2_node, physical_total_residue, q2);
    double score = 0.0;
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      const Complex linear =
          physical_total_residue[interval][entry] * (1.0 - fraction) +
          physical_total_residue[interval + 1][entry] * fraction;
      const double scale =
          std::max({std::abs(cubic[entry]), std::abs(linear),
                    std::sqrt(tolerance) * entry_scale[entry], 1.0e-300});
      score = std::max(score, std::abs(cubic[entry] - linear) / scale);
    }
    curvature.push_back({score, interval});
  }
  std::sort(curvature.begin(), curvature.end(),
            [](const auto &first, const auto &second) {
              return first.first > second.first;
            });
  for (const auto &[score, interval] : curvature) {
    static_cast<void>(score);
    if (index.size() == budget) {
      break;
    }
    add(interval);
  }
  std::sort(index.begin(), index.end());
  return index;
}

// Restore the physical pole behavior of one regularized total amplitude
ProtonHelicityMatrix
DeRegularizePhysicalResidue(const ProtonHelicityMatrix &residue,
                            const double q) {
  const auto pole_factor = PoleRegularizationFactors(q);
  ProtonHelicityMatrix amplitude = residue;
  for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
    amplitude[entry] /= pole_factor[AbsoluteHarmonic(entry)];
  }
  return amplitude;
}

// Validate pole-regular total entries and the physical matrix at midpoints
void ValidatePhysicalTotalResidueInterpolation(
    const MEikonalMatrix &strong_runtime, const MDirac &dirac,
    const std::vector<MParticle> &initialstate,
    const std::vector<ProtonHelicityMatrix> &profile,
    const std::vector<double> &b_node, const bool logarithmic_b, const double s,
    const double sqrt_lambda,
    const ResidualBornData &residual,
    const std::vector<ProtonHelicityMatrix> &nonlinear_polynomial,
    const double point_born_coefficient, const double eta,
    const double point_ir_scale, const ElasticCNINumerics &numerics,
    const std::vector<double> &q2_node,
    const std::vector<ProtonHelicityMatrix> &physical_total_residue,
    const form::ParamStore &structure) {
  if (q2_node.size() < 4 || physical_total_residue.size() != q2_node.size()) {
    throw std::invalid_argument(
        "ValidatePhysicalTotalResidueInterpolation: invalid residue table");
  }
  std::array<double, kHelicityEntries> entry_scale{};
  for (const auto &matrix : physical_total_residue) {
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      entry_scale[entry] =
          std::max(entry_scale[entry], std::abs(matrix[entry]));
    }
  }
  const double tolerance = numerics.interp_rel_tol;
  const double epsilon = std::numeric_limits<double>::epsilon();
  WorstInterpolation worst_residue;
  WorstInterpolation worst_matrix;

  const auto validation_interval = ResidueValidationIntervals(
      q2_node, physical_total_residue, entry_scale, tolerance, numerics.logKT2);
  for (const std::size_t index : validation_interval) {
    const double q2 = numerics.logKT2
                          ? std::sqrt(q2_node[index] * q2_node[index + 1])
                          : 0.5 * (q2_node[index] + q2_node[index + 1]);
    const double q = std::sqrt(q2);
    const ProtonHelicityMatrix direct_finite = FiniteCorrectionAtQ(
        profile, b_node, logarithmic_b, sqrt_lambda, residual,
        nonlinear_polynomial, eta, point_ir_scale, q, numerics);
    const ProtonHelicityMatrix direct_cni =
        PhysicalCNIAmplitudeAtQ(dirac, initialstate, s, q2, direct_finite,
                                point_born_coefficient, eta, point_ir_scale,
                                structure);
    const ProtonHelicityMatrix strong = strong_runtime.HelicityMatrix(q2, 0, 0);
    ProtonHelicityMatrix direct_physical = strong;
    gra::AddScaled(direct_physical, direct_cni, Complex(1.0));
    const ProtonHelicityMatrix direct_residue =
        RegularizePhysicalHelicityMatrix(direct_physical, q);
    const ProtonHelicityMatrix interpolated_residue =
        InterpolateHelicityMatrixTable(q2_node, physical_total_residue, q2);
    for (std::size_t entry = 0; entry < kHelicityEntries; ++entry) {
      const double error =
          std::abs(direct_residue[entry] - interpolated_residue[entry]);
      const double scale =
          std::max({std::abs(direct_residue[entry]),
                    std::abs(interpolated_residue[entry]),
                    std::sqrt(tolerance) * entry_scale[entry]});
      const double allowance =
          tolerance * scale +
          256.0 * epsilon * std::max(1.0, entry_scale[entry]);
      worst_residue.Update(error, allowance, q2, entry);
    }

    const ProtonHelicityMatrix interpolated_physical =
        DeRegularizePhysicalResidue(interpolated_residue, q);
    ProtonHelicityMatrix difference = direct_physical;
    gra::AddScaled(difference, interpolated_physical, Complex(-1.0));
    const double error = std::sqrt(gra::SquaredNorm(difference));
    const double component_scale = std::sqrt(gra::SquaredNorm(strong)) +
                                   std::sqrt(gra::SquaredNorm(direct_cni));
    const double scale =
        std::max({std::sqrt(gra::SquaredNorm(direct_physical)),
                  std::sqrt(gra::SquaredNorm(interpolated_physical)),
                  std::sqrt(tolerance) * component_scale});
    const double allowance =
        tolerance * scale + 512.0 * epsilon * std::max(1.0, component_scale);
    worst_matrix.Update(error, allowance, q2);
  }
  if (!std::isfinite(worst_residue.ratio) ||
      !std::isfinite(worst_matrix.ratio) || worst_residue.ratio > 1.0 ||
      worst_matrix.ratio > 1.0) {
    std::ostringstream message;
    message << std::scientific << std::setprecision(6)
            << "MElasticCNI: physical total interpolation validation failed; "
            << "entry ratio=" << worst_residue.ratio
            << " at |t|=" << worst_residue.q2
            << ", helicity entry=" << worst_residue.entry
            << ", error=" << worst_residue.error
            << ", allowance=" << worst_residue.allowance
            << "; full-matrix ratio=" << worst_matrix.ratio
            << " at |t|=" << worst_matrix.q2
            << ", error=" << worst_matrix.error
            << ", allowance=" << worst_matrix.allowance
            << "; increase CNI.NumberKT2";
    throw std::runtime_error(message.str());
  }
}

// Construct the complete physics and numerical CNI cache key
std::string CNIKey(const SoftModel &model, const double s,
                   const std::vector<MParticle> &initialstate,
                   const ElasticCNINumerics &numerics,
                   const std::vector<double> &b_node,
                   const std::string &strong_fingerprint,
                   const bool strong_log_b, const double min_abs_t,
                   const double max_abs_t, const double maximum_q,
                   const double point_ir_scale,
                   const form::ParamStore &structure) {
  std::ostringstream stream;
  stream << kCacheVersion << ';' << model.Fingerprint() << ';'
         << structure.EM << ';' << strong_fingerprint << ';';
  AppendCacheDouble(stream, qed::alpha_0);
  AppendCacheDouble(stream, s);
  for (const auto &particle : initialstate) {
    stream << particle.pdg << ';' << particle.chargeX3 << ';' << particle.spinX2
           << ';';
    AppendCacheDouble(stream, particle.mass);
  }
  AppendCacheDouble(stream, min_abs_t);
  AppendCacheDouble(stream, max_abs_t);
  AppendCacheDouble(stream, point_ir_scale);
  AppendCacheDouble(stream, maximum_q);
  stream << numerics.NumberKT2 << ';' << numerics.logKT2 << ';';
  AppendCacheDouble(stream, numerics.FBIntegralMaxKT);
  stream << numerics.FBIntegralN << ';';
  AppendCacheDouble(stream, numerics.tail_max_qb);
  stream << numerics.TailIntegralN << ';';
  stream << numerics.tail_order << ';';
  AppendCacheDouble(stream, numerics.tail_split_qb);
  AppendCacheDouble(stream, numerics.interp_rel_tol);
  stream << strong_log_b << ';' << b_node.size() << ';';
  for (const double b : b_node) {
    AppendCacheDouble(stream, b);
  }
  return stream.str();
}

// Compute the persistent filename for one full CNI cache key
std::string CNIFilename(const std::string &key,
                        const std::vector<MParticle> &initialstate) {
  std::ostringstream name;
  name << aux::GetBasePath(2) << "/eikonal/ELASTIC_CNI";
  for (const auto &particle : initialstate) {
    name << '_' << particle.pdg;
  }
  name << '_' << aux::djb2hash(key) << ".json";
  return name.str();
}

// Flatten one complex helicity table for compact JSON storage
std::vector<double>
FlattenHelicityTable(const std::vector<ProtonHelicityMatrix> &matrix_table) {
  std::vector<double> data;
  data.reserve(2 * kHelicityEntries * matrix_table.size());
  for (const auto &matrix : matrix_table) {
    for (const Complex value : matrix) {
      if (!std::isfinite(std::real(value)) ||
          !std::isfinite(std::imag(value))) {
        throw std::runtime_error(
            "MElasticCNI cache contains a non-finite amplitude");
      }
      data.push_back(std::real(value));
      data.push_back(std::imag(value));
    }
  }
  return data;
}

// Restore and validate one pole-regular physical total cache
bool ReadCNICache(const std::string &filename, const std::string &key,
                  const std::size_t points, std::vector<double> &q2_node,
                  std::vector<ProtonHelicityMatrix> &physical_total_residue) {
  try {
    std::ifstream input(filename);
    if (!input.good()) {
      return false;
    }
    nlohmann::json cache;
    input >> cache;
    if (cache.at("version").get<std::size_t>() != kCacheVersion ||
        cache.at("key").get<std::string>() != key) {
      return false;
    }
    q2_node = cache.at("q2_node").get<std::vector<double>>();
    const auto data =
        cache.at("physical_total_residue").get<std::vector<double>>();
    if (q2_node.size() != points ||
        data.size() != 2 * kHelicityEntries * points) {
      return false;
    }
    math::ValidateInterpolationGrid(q2_node);
    physical_total_residue.assign(points, ProtonHelicityMatrix{});
    std::size_t data_index = 0;
    for (auto &matrix : physical_total_residue) {
      for (Complex &value : matrix) {
        value = Complex(data[data_index], data[data_index + 1]);
        data_index += 2;
        if (!std::isfinite(std::real(value)) ||
            !std::isfinite(std::imag(value))) {
          return false;
        }
      }
    }
    return true;
  } catch (const std::exception &) {
    return false;
  }
}

// Write one complete CNI cache and preserve any previous file
void WriteCNICache(
    const std::string &filename, const std::string &key,
    const std::vector<double> &q2_node,
    const std::vector<ProtonHelicityMatrix> &physical_total_residue) {
  nlohmann::json cache;
  cache["version"] = kCacheVersion;
  cache["key"] = key;
  cache["q2_node"] = q2_node;
  cache["physical_total_residue"] =
      FlattenHelicityTable(physical_total_residue);

  eikonal::PublishCache(filename, cache, "MElasticCNI cache");
}

// Validate one supplied elastic event and return its positive transfer
double ValidateElasticEvent(const double expected_s,
                            const std::vector<MParticle> &initialstate,
                            const M4Vec &p1_in, const M4Vec &p2_in,
                            const M4Vec &p1_out, const M4Vec &p2_out) {
  if (initialstate.size() != 2) {
    throw std::logic_error(
        "MElasticCNI::PhysicalHelicityMatrix: invalid initial state");
  }
  const auto on_shell = [](const M4Vec &p, const double mass) {
    const double mass2 = mass * mass;
    // Allow rounding of the external four momentum before the invariant subtraction
    const double tolerance = 1.0e-8 * std::max(1.0, mass2) +
                             16.0 * std::numeric_limits<double>::epsilon() *
                                 (math::pow2(p.E()) + p.P3mod2());
    return p.E() > 0.0 && std::isfinite(p.E()) && std::isfinite(p.M2()) &&
           std::isfinite(tolerance) && std::abs(p.M2() - mass2) <= tolerance;
  };
  if (!on_shell(p1_in, initialstate[0].mass) || !on_shell(p2_in, initialstate[1].mass) ||
      !on_shell(p1_out, initialstate[0].mass) || !on_shell(p2_out, initialstate[1].mass)) {
    throw std::invalid_argument("MElasticCNI::PhysicalHelicityMatrix: off-shell external momentum");
  }
  if (p1_in.Pt() > 1.0e-10 || p2_in.Pt() > 1.0e-10 || !(p1_in.Pz() > 0.0) || !(p2_in.Pz() < 0.0)) {
    throw std::invalid_argument(
        "MElasticCNI::PhysicalHelicityMatrix: ordered collinear beams are "
        "required");
  }
  if (!math::CheckEMC(p1_in + p2_in - p1_out - p2_out)) {
    throw std::invalid_argument("MElasticCNI::PhysicalHelicityMatrix: four-momentum is not conserved");
  }
  const M4Vec  incoming    = p1_in + p2_in;
  const double event_s     = incoming.M2();
  const double s_tolerance = 1.0e-10 * std::max({1.0, std::abs(event_s), std::abs(expected_s)});
  if (!std::isfinite(event_s) || std::abs(event_s - expected_s) > s_tolerance) {
    throw std::invalid_argument(
        "MElasticCNI::PhysicalHelicityMatrix: event energy disagrees with "
        "the initialized eikonal");
  }
  const double t1 = (p1_out - p1_in).M2();
  const double t2 = (p2_out - p2_in).M2();
  if (!std::isfinite(t1) || !std::isfinite(t2) || !(t1 < 0.0) || !(t2 < 0.0) ||
      std::abs(t1 - t2) > 1.0e-8 * std::max({1.0, std::abs(t1), std::abs(t2)})) {
    throw std::invalid_argument(
        "MElasticCNI::PhysicalHelicityMatrix: inconsistent elastic transfer");
  }
  return -0.5 * (t1 + t2);
}

// Compute the renormalized all-orders point-Coulomb phase ratio
Complex PointCoulombPhaseRatioImpl(const double eta,
                                   const double point_ir_scale,
                                   const double abs_t) {
  // [REFERENCE: arXiv:2001.10227]
  // In exp(2 i gamma_E eta) Gamma(1+i eta)/Gamma(1-i eta), the
  // Euler-constant term cancels exactly. The remaining convergent log-Gamma
  // series starts at eta^3. Truncation after eta^19 is below 1e-40 here
  Complex logarithm =
      math::zi * eta * std::log(point_ir_scale * point_ir_scale / abs_t);
  const Complex z = math::zi * eta;
  for (int order = 3; order <= 19; order += 2) {
    logarithm -= 2.0 * std::riemann_zeta(static_cast<double>(order)) *
                 std::pow(z, order) / static_cast<double>(order);
  }
  return std::exp(logarithm);
}

} // namespace

// Construct an empty mutable runtime before immutable publication
MElasticCNI::MElasticCNI(const double s, std::vector<MParticle> initialstate,
                         const double min_abs_t, const double max_abs_t,
                         const double point_ir_scale, const double sqrt_lambda,
                         const double eta)
    : s_(s), initialstate_(std::move(initialstate)), min_abs_t_(min_abs_t),
      max_abs_t_(max_abs_t),
      point_ir_scale_(point_ir_scale), eta_(eta),
      point_born_coefficient_(-8.0 * math::PI * sqrt_lambda * eta) {}

// Compute the symmetric residual electromagnetic sandwich of the strong S
MElasticCNI::Matrix
MElasticCNI::SymmetricShortRangeSMatrix(const Matrix &strong_s,
                                        const Matrix &chi_residual) {
  if (strong_s.size_row() == 0 || strong_s.size_row() != strong_s.size_col() ||
      chi_residual.size_row() != strong_s.size_row() ||
      chi_residual.size_col() != strong_s.size_col()) {
    throw std::invalid_argument(
        "MElasticCNI::SymmetricShortRangeSMatrix: matrix dimensions disagree");
  }
  const Matrix half = (chi_residual * (0.5 * math::zi)).Exp();
  return half * strong_s * half;
}

// Compute the finite impact-parameter remainder after both Born subtractions
MElasticCNI::Matrix
MElasticCNI::FiniteRemainderProfile(const Matrix &strong_s,
                                    const Matrix &chi_residual,
                                    const double chi_point) {
  if (!std::isfinite(chi_point)) {
    throw std::invalid_argument(
        "MElasticCNI::FiniteRemainderProfile: non-finite point eikonal");
  }
  const Matrix identity = Matrix::IdentityMatrix(strong_s.size_row());
  const Matrix strong_t = math::zi * (identity - strong_s);
  const Matrix short_t =
      math::zi *
      (identity - SymmetricShortRangeSMatrix(strong_s, chi_residual));
  return short_t * std::exp(math::zi * chi_point) - chi_residual - strong_t;
}

// Compute the IR-renormalized analytic point-Coulomb phase ratio
MElasticCNI::Complex MElasticCNI::RenormalizedPointCoulombPhase(
    const double eta, const double point_ir_scale, const double abs_t) {
  if (!std::isfinite(eta) || !std::isfinite(point_ir_scale) ||
      !std::isfinite(abs_t) || point_ir_scale <= 0.0 || abs_t <= 0.0) {
    throw std::invalid_argument(
        "MElasticCNI::RenormalizedPointCoulombPhase: invalid input");
  }
  return PointCoulombPhaseRatioImpl(eta, point_ir_scale, abs_t);
}

// Construct or load the pole-regular total physical helicity table
std::shared_ptr<const MElasticCNI>
MElasticCNI::Build(const MEikonalMatrix &strong_runtime,
                   const SoftModelPtr &model, const double s,
                   const std::vector<MParticle> &initialstate,
                   const ElasticCNINumerics &numerics, const bool strong_log_b,
                   const double min_abs_t, const double max_abs_t,
                   const form::ParamStore &structure) {
  if (model == nullptr || initialstate.size() != 2 || !std::isfinite(s) ||
      s <= 0.0 || !std::isfinite(min_abs_t) || !std::isfinite(max_abs_t) ||
      min_abs_t <= 0.0 || max_abs_t <= min_abs_t || numerics.NumberKT2 < 3 ||
      !std::isfinite(numerics.FBIntegralMaxKT) ||
      numerics.FBIntegralMaxKT <= 0.0 || numerics.FBIntegralN < 2 ||
      !std::isfinite(numerics.tail_max_qb) || numerics.tail_max_qb <= 0.0 ||
      numerics.TailIntegralN < 2 ||
      numerics.TailIntegralN % 2 != 0 ||
      numerics.tail_order < 1 || static_cast<unsigned int>(numerics.tail_order) > numerics.TailIntegralN ||
      !std::isfinite(numerics.tail_split_qb) || numerics.tail_split_qb <= 0.0 ||
      !std::isfinite(numerics.interp_rel_tol) ||
      numerics.interp_rel_tol <= 0.0) {
    throw std::invalid_argument("MElasticCNI::Build: invalid configuration");
  }
  if (std::abs(initialstate[0].pdg) != 2212 ||
      std::abs(initialstate[1].pdg) != 2212) {
    throw std::invalid_argument(
        "MElasticCNI::Build: elastic CNI currently supports protons and "
        "antiprotons only");
  }
  const auto emitter1 =
      qed::Emitter(initialstate[0], -min_abs_t, "MElasticCNI::Build beam 1",
                   structure);
  const auto emitter2 =
      qed::Emitter(initialstate[1], -min_abs_t, "MElasticCNI::Build beam 2",
                   structure);
  const double momentum = CenterOfMassMomentum(s, emitter1.mass, emitter2.mass);
  const double phase_space_max_t = 4.0 * momentum * momentum;
  if (max_abs_t > phase_space_max_t * (1.0 + 1.0e-12)) {
    throw std::invalid_argument(
        "MElasticCNI::Build: max_abs_t exceeds two-body phase space");
  }
  const double maximum_q = std::min(numerics.FBIntegralMaxKT, 2.0 * momentum);
  if (maximum_q < std::sqrt(max_abs_t) * (1.0 - 1.0e-12)) {
    throw std::invalid_argument(
        "MElasticCNI::Build: CNI.FBIntegralMaxKT is below the output range");
  }
  const auto &b_node = strong_runtime.ImpactParameterNodes();
  if (b_node.size() < 3 || !(b_node.front() > 0.0) ||
      !(b_node.back() > b_node.front())) {
    throw std::invalid_argument(
        "MElasticCNI::Build: invalid hadronic impact-parameter grid");
  }
  if (numerics.tail_max_qb <= std::sqrt(max_abs_t) * b_node.back()) {
    throw std::invalid_argument(
        "MElasticCNI::Build: CNI.tail_max_qb is below the output "
        "range at the hadronic impact-parameter boundary");
  }
  const double delta =
      s - emitter1.mass * emitter1.mass - emitter2.mass * emitter2.mass;
  const double beta = kinematics::beta12(s, emitter1.mass, emitter2.mass);
  if (!std::isfinite(beta) || !(beta > 0.0)) {
    throw std::invalid_argument("MElasticCNI::Build: initial state is below threshold");
  }
  const double sqrt_lambda = s * beta;
  // The relativistic Coulomb parameter is alpha z1 z2 (s-m1^2-m2^2)/sqrt(lambda)
  const double eta = qed::alpha_0 * emitter1.charge * emitter2.charge * delta / sqrt_lambda;
  auto runtime = std::shared_ptr<MElasticCNI>(new MElasticCNI(
      s, initialstate, min_abs_t, max_abs_t, kPointIRScale, sqrt_lambda, eta));

  const std::string cache_key =
      CNIKey(*model, s, initialstate, numerics, b_node,
             strong_runtime.RuntimeFingerprint(), strong_log_b, min_abs_t,
             max_abs_t, maximum_q, kPointIRScale, structure);
  const std::string cache_filename = CNIFilename(cache_key, initialstate);
  std::error_code filesystem_error;
  aux::CreateDirectory(std::filesystem::path(cache_filename).parent_path().string(), &filesystem_error);
  std::unique_ptr<eikonal::MCacheLock> cache_lock;
  if (!filesystem_error &&
      ReadCNICache(cache_filename, cache_key, numerics.NumberKT2 + 1,
                   runtime->q2_node_, runtime->physical_total_residue_)) {
    std::cout << "Loaded elastic CNI cache: " << cache_filename << std::endl;
    return runtime;
  }
  if (!filesystem_error) {
    cache_lock =
        std::make_unique<eikonal::MCacheLock>(cache_filename + ".lock");
    if (!cache_lock->Acquired()) {
      cache_lock.reset();
    } else if (ReadCNICache(cache_filename, cache_key, numerics.NumberKT2 + 1,
                            runtime->q2_node_,
                            runtime->physical_total_residue_)) {
      std::cout << "Loaded elastic CNI cache: " << cache_filename << std::endl;
      return runtime;
    }
  }

  std::cout << "Building full-spin elastic CNI table" << std::endl;
  const MDirac dirac("DIRAC");
  const auto residual = BuildResidualBornTable(
      dirac, initialstate, s, sqrt_lambda, maximum_q, numerics.FBIntegralN,
      runtime->point_born_coefficient_, structure);
  const auto nonlinear_polynomial =
      BuildNonlinearTailPolynomial(residual, b_node.back());
  const auto profile =
      BuildFiniteRemainderProfile(strong_runtime, residual, maximum_q,
                                  numerics.FBIntegralN, sqrt_lambda, kPointIRScale, eta);
  std::vector<ProtonHelicityMatrix> finite_correction;
  BuildFiniteCorrectionTable(profile, b_node, strong_log_b, sqrt_lambda,
                             residual, nonlinear_polynomial, eta, kPointIRScale,
                             numerics, min_abs_t, max_abs_t, runtime->q2_node_,
                             finite_correction);
  math::ValidateInterpolationGrid(runtime->q2_node_);
  BuildPhysicalTotalResidueTable(
      strong_runtime, dirac, initialstate, s, runtime->q2_node_,
      finite_correction, runtime->point_born_coefficient_, eta, kPointIRScale,
      structure, runtime->physical_total_residue_);
  ValidatePhysicalTotalResidueInterpolation(
      strong_runtime, dirac, initialstate, profile, b_node, strong_log_b, s,
      sqrt_lambda, residual, nonlinear_polynomial,
      runtime->point_born_coefficient_, eta, kPointIRScale, numerics, runtime->q2_node_,
      runtime->physical_total_residue_, structure);

  if (cache_lock) {
    try {
      WriteCNICache(cache_filename, cache_key, runtime->q2_node_,
                    runtime->physical_total_residue_);
      std::cout << "Saved elastic CNI cache: " << cache_filename << std::endl;
    } catch (const std::exception &error) {
      std::cerr << "WARNING: elastic CNI cache not saved: " << error.what()
                << std::endl;
    }
  }
  return runtime;
}

// Compute the analytic pure point-Coulomb correction beyond one photon
ProtonHelicityMatrix
MElasticCNI::PointHigherOrderHelicityMatrix(const double abs_t) const {
  if (!std::isfinite(abs_t) || abs_t <= 0.0) {
    throw std::invalid_argument(
        "MElasticCNI::PointHigherOrderHelicityMatrix: abs_t must be positive");
  }
  const Complex correction =
      point_born_coefficient_ / abs_t *
      (RenormalizedPointCoulombPhase(eta_, point_ir_scale_, abs_t) -
       Complex(1.0));
  ProtonHelicityMatrix matrix{};
  for (std::size_t diagonal = 0; diagonal < kHelicityDimension; ++diagonal) {
    matrix[kHelicityDimension * diagonal + diagonal] = correction;
  }
  return matrix;
}

// Compute the complete physical helicity amplitude for one elastic event
ProtonHelicityMatrix
MElasticCNI::PhysicalHelicityMatrix(const M4Vec &p1_in, const M4Vec &p2_in,
                                    const M4Vec &p1_out,
                                    const M4Vec &p2_out) const {
  double abs_t =
      ValidateElasticEvent(s_, initialstate_, p1_in, p2_in, p1_out, p2_out);
  // Scale roundoff tolerance to each bound to preserve the small-t Coulomb pole
  const double tolerance = 1.0e-10;
  if (abs_t < min_abs_t_ * (1.0 - tolerance) ||
      abs_t > max_abs_t_ * (1.0 + tolerance)) {
    throw std::out_of_range(
        "MElasticCNI::PhysicalHelicityMatrix: transfer is outside the "
        "initialized CNI range");
  }
  abs_t = std::clamp(abs_t, min_abs_t_, max_abs_t_);
  const double azimuth = (p1_out - p1_in).Phi();
  if (q2_node_.size() < 4 ||
      physical_total_residue_.size() != q2_node_.size()) {
    throw std::logic_error(
        "MElasticCNI::PhysicalHelicityMatrix: invalid physical total table");
  }
  const ProtonHelicityMatrix residue =
      InterpolateHelicityMatrixTable(q2_node_, physical_total_residue_, abs_t);
  return RotateProtonHelicityMatrix(
      DeRegularizePhysicalResidue(residue, std::sqrt(abs_t)), azimuth);
}

} // namespace gra
