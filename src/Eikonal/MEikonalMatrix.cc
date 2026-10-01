// Coupled-channel matrix eikonals and Good Walker geometry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <tuple>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Eikonal/MEikonalCache.h"
#include "Graniitti/Eikonal/MEikonalHelicity.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MAlgorithms.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MInterpolation.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MJsonZip.h"
#include "Graniitti/Tech/MTimer.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

namespace gra {

namespace {

using Complex = std::complex<double>;
using Matrix = MEikonalMatrix::Matrix;
using PairSpinBank = MEikonalMatrix::PairSpinBank;
using PairSpinTable = std::array<std::vector<Matrix>, 16>;
using math::CompositeWeight;
using math::ParallelFor;
using math::zi;
constexpr std::size_t kNHelicity = 5;
constexpr std::size_t kNSpinEntries = 16;
enum HelicityIndex : std::size_t { kPhi1, kPhi2, kPhi3, kPhi4, kPhi5 };
constexpr std::array<int, kNHelicity> kHarmonic = {0, 0, 0, 2, 1};
constexpr std::size_t kCacheVersion = 1;

// Check the configured contractivity bound for both built and cached S matrices
void CheckContractivity(const double sigma_max, const MEikonalNumerics &numerics) {
  if (sigma_max <= 1.0 + numerics.unitarity_tolerance) {
    return;
  }
  std::ostringstream message;
  message << std::setprecision(12)
          << "S-matrix contractivity violation: max singular value |S| = "
          << sigma_max << ", tolerance = " << numerics.unitarity_tolerance;
  if (numerics.strict_unitarity) {
    throw std::runtime_error(message.str());
  }
  std::cout << rang::fg::yellow << "WARNING: " << message.str()
            << rang::fg::reset << std::endl;
}

// Compute the nonnegative harmonic of one oriented spin transition
constexpr std::size_t SpinHarmonic(const std::size_t entry) {
  const int harmonic =
      CanonicalProtonHelicityTransitions()[entry].azimuth_harmonic;
  return static_cast<std::size_t>(harmonic < 0 ? -harmonic : harmonic);
}

// Append one floating-point value in an exact cache-key representation
void AppendCacheDouble(std::ostringstream &stream, const double value) {
  stream << std::hexfloat << value << ';';
}

// Construct the full physics and numerical key for one matrix table cache
std::string MatrixCacheKey(const double s,
                           const std::vector<MParticle> &initialstate,
                           const MEikonalNumerics &numerics,
                           const SoftModel &model) {
  std::ostringstream stream;
  AppendCacheDouble(stream, s);
  for (const auto &particle : initialstate) {
    stream << particle.pdg << ';';
    AppendCacheDouble(stream, particle.mass);
  }
  AppendCacheDouble(stream, numerics.MinKT2);
  AppendCacheDouble(stream, numerics.MaxKT2);
  stream << numerics.NumberKT2 << ';' << numerics.logKT2 << ';';
  AppendCacheDouble(stream, numerics.MinBT);
  AppendCacheDouble(stream, numerics.MaxBT);
  stream << numerics.NumberBT << ';' << numerics.logBT << ';';
  AppendCacheDouble(stream, numerics.FBIntegralMinKT);
  AppendCacheDouble(stream, numerics.FBIntegralMaxKT);
  stream << numerics.FBIntegralN << ';';
  stream << numerics.strict_unitarity << ';';
  AppendCacheDouble(stream, numerics.unitarity_tolerance);
  const auto &mixing = model.GoodWalker().MixingMatrix();
  for (std::size_t row = 0; row < mixing.size_row(); ++row) {
    for (std::size_t col = 0; col < mixing.size_col(); ++col) {
      AppendCacheDouble(stream, mixing(row, col));
    }
  }
  stream << model.Fingerprint();
  return stream.str();
}

// Compute the persistent cache filename for one full matrix-table key
std::string MatrixCacheFilename(const std::string &key,
                                const std::vector<MParticle> &initialstate,
                                const std::size_t channels) {
  const unsigned long hash = aux::djb2hash(key);
  std::ostringstream name;
  name << aux::GetBasePath(2) << "/eikonal/MATRIX_N" << channels;
  for (const auto &particle : initialstate) {
    name << '_' << particle.pdg;
  }
  name << '_' << hash << ".json";
  return name.str();
}

// Compress one fixed spin array of matrix tables
template <typename Bank>
nlohmann::json CompressMatrixTables(const Bank &source) {
  nlohmann::json output = nlohmann::json::array();
  for (const auto &table : source) {
    output.push_back(MJsonZip::CompressMatrixBank(table));
  }
  return output;
}

// Compress one fixed spin array of zero-transfer matrices
template <typename Bank>
nlohmann::json CompressMatrices(const Bank &source) {
  nlohmann::json output = nlohmann::json::array();
  for (const auto &matrix : source) {
    output.push_back(MJsonZip::CompressMatrix(matrix));
  }
  return output;
}

// Restore one fixed spin array of matrix tables
template <typename Bank>
bool ReadMatrixTables(const nlohmann::json &source, const std::size_t points,
                      const std::size_t dimension, Bank &target) {
  if (!source.is_array() || source.size() != target.size()) {
    return false;
  }
  for (const auto &entry : indices(target)) {
    if (!MJsonZip::DecompressMatrixBank(source[entry], target[entry], points,
                                        dimension, dimension)) {
      return false;
    }
  }
  return true;
}

// Restore one fixed spin array of zero-transfer matrices
template <typename Bank>
bool ReadMatrices(const nlohmann::json &source, const std::size_t dimension,
                  Bank &target) {
  if (!source.is_array() || source.size() != target.size()) {
    return false;
  }
  for (const auto &entry : indices(target)) {
    if (!MJsonZip::DecompressMatrix(source[entry], target[entry], dimension,
                                    dimension)) {
      return false;
    }
  }
  return true;
}

// Write all immutable coupled-channel tables to one JSON cache file
template <typename Tables>
void WriteMatrixJsonCache(const std::string &filename, const std::string &key,
                          const Tables &runtime, const bool has_amplitude) {
  nlohmann::json cache = {
      {"version", kCacheVersion},
      {"key", key},
      {"channels", runtime.N},
      {"dimension", runtime.D},
      {"has_amplitude", has_amplitude},
      {"max_singular_value", runtime.max_singular_value},
      {"U", MJsonZip::CompressVector(runtime.U)},
      {"b_node", MJsonZip::CompressVector(runtime.b_node)},
      {"chi_even_b", MJsonZip::CompressMatrixBank(runtime.chi_even_b)},
      {"chi_odd_b", MJsonZip::CompressMatrixBank(runtime.chi_odd_b)},
      {"chi_cut_b", MJsonZip::CompressMatrixBank(runtime.chi_cut_b)}};
  cache["amplitude_spin_b"] = CompressMatrixTables(runtime.amplitude_spin_b);
  cache["screening_spin_b"] = CompressMatrixTables(runtime.screening_spin_b);
  if (has_amplitude) {
    cache["q2_node"] = MJsonZip::CompressVector(runtime.q2_node);
    cache["amplitude_spin_q"] = CompressMatrixTables(runtime.amplitude_spin_q);
    cache["amplitude_spin_zero"] = CompressMatrices(runtime.amplitude_spin_zero);
    cache["crossing_even_spin_q"] =
        CompressMatrixTables(runtime.crossing_even_spin_q);
    cache["crossing_even_spin_zero"] =
        CompressMatrices(runtime.crossing_even_spin_zero);
    cache["crossing_odd_spin_q"] =
        CompressMatrixTables(runtime.crossing_odd_spin_q);
    cache["crossing_odd_spin_zero"] =
        CompressMatrices(runtime.crossing_odd_spin_zero);
    cache["screening_spin_q"] = CompressMatrixTables(runtime.screening_spin_q);
    cache["screening_spin_zero"] = CompressMatrices(runtime.screening_spin_zero);
  }

  eikonal::PublishCache(filename, cache, "MEikonalMatrix cache");
}

// Validate immutable impact-parameter and momentum interpolation grids once
template <typename Tables>
void ValidateInterpolationNodes(const Tables &runtime,
                                const bool require_amplitude) {
  math::ValidateInterpolationGrid(runtime.b_node);
  if (require_amplitude) {
    math::ValidateInterpolationGrid(runtime.q2_node);
  }
}

// Load one validated coupled-channel table cache from its JSON representation
template <typename Tables>
bool ReadMatrixJsonCache(const std::string &filename, const std::string &key,
                         const std::size_t channels,
                         const std::size_t dimension,
                         const std::size_t b_points, const std::size_t q_points,
                         const std::size_t cut_dimension,
                         const bool require_amplitude, Tables &runtime) {
  try {
    std::ifstream input(filename);
    if (!input.is_open()) {
      return false;
    }
    nlohmann::json cache;
    input >> cache;
    if (!cache.is_object() || cache.value("version", 0U) != kCacheVersion ||
        cache.value("key", std::string()) != key ||
        cache.at("channels").get<std::size_t>() != channels ||
        cache.at("dimension").get<std::size_t>() != dimension ||
        !cache.at("has_amplitude").is_boolean() ||
        (require_amplitude && !cache.at("has_amplitude").get<bool>()) ||
        !cache.at("max_singular_value").is_number()) {
      return false;
    }

    runtime.N = channels;
    runtime.D = dimension;
    runtime.max_singular_value = cache.at("max_singular_value").get<double>();
    if (!std::isfinite(runtime.max_singular_value) ||
        !MJsonZip::DecompressVector(cache.at("U"), runtime.U, dimension) ||
        !MJsonZip::DecompressVector(cache.at("b_node"), runtime.b_node,
                                    b_points) ||
        !MJsonZip::DecompressMatrixBank(cache.at("chi_even_b"),
                                        runtime.chi_even_b, b_points, dimension,
                                        dimension) ||
        !MJsonZip::DecompressMatrixBank(cache.at("chi_odd_b"),
                                        runtime.chi_odd_b, b_points, dimension,
                                        dimension) ||
        !MJsonZip::DecompressMatrixBank(cache.at("chi_cut_b"),
                                        runtime.chi_cut_b, b_points,
                                        cut_dimension, cut_dimension) ||
        !ReadMatrixTables(cache.at("amplitude_spin_b"), b_points, dimension,
                          runtime.amplitude_spin_b) ||
        !ReadMatrixTables(cache.at("screening_spin_b"), b_points, dimension,
                          runtime.screening_spin_b)) {
      return false;
    }

    if (!cache.at("has_amplitude").get<bool>()) {
      ValidateInterpolationNodes(runtime, false);
      return true;
    }
    if (!MJsonZip::DecompressVector(cache.at("q2_node"), runtime.q2_node,
                                    q_points) ||
        !ReadMatrixTables(cache.at("amplitude_spin_q"), q_points, dimension,
                          runtime.amplitude_spin_q) ||
        !ReadMatrices(cache.at("amplitude_spin_zero"), dimension,
                      runtime.amplitude_spin_zero) ||
        !ReadMatrixTables(cache.at("crossing_even_spin_q"), q_points, dimension,
                          runtime.crossing_even_spin_q) ||
        !ReadMatrices(cache.at("crossing_even_spin_zero"), dimension,
                      runtime.crossing_even_spin_zero) ||
        !ReadMatrixTables(cache.at("crossing_odd_spin_q"), q_points, dimension,
                          runtime.crossing_odd_spin_q) ||
        !ReadMatrices(cache.at("crossing_odd_spin_zero"), dimension,
                      runtime.crossing_odd_spin_zero) ||
        !ReadMatrixTables(cache.at("screening_spin_q"), q_points, dimension,
                          runtime.screening_spin_q) ||
        !ReadMatrices(cache.at("screening_spin_zero"), dimension,
                      runtime.screening_spin_zero)) {
      return false;
    }
    ValidateInterpolationNodes(runtime, true);
    return true;
  } catch (const std::exception &) {
    return false;
  }
}

// Compute the azimuthal harmonic phase with the selected Fourier sign
Complex HarmonicPhase(const int n, const int sign) {
  if (n == 0) {
    return 1.0;
  }
  if (n == 1) {
    return static_cast<double>(sign) * zi;
  }
  if (n == 2) {
    return -1.0;
  }
  return std::pow(static_cast<double>(sign) * zi, n);
}

// Interpolate one oriented spin block with its small momentum limit
template <typename Tables>
Matrix InterpolateSpinBlock(const Tables &runtime, const PairSpinTable &bank,
                            const PairSpinBank &zero, const std::size_t entry,
                            const double kt2) {
  if (entry >= kNSpinEntries || !std::isfinite(kt2) || kt2 < 0.0 ||
      runtime.q2_node.empty() || bank[entry].size() != runtime.q2_node.size()) {
    throw std::invalid_argument("InterpolateSpinBlock: invalid input");
  }
  if (std::fpclassify(kt2) == FP_ZERO) {
    return zero[entry];
  }
  if (kt2 <= runtime.q2_node.front()) {
    const std::size_t harmonic = SpinHarmonic(entry);
    const double power = harmonic == 1 ? 0.5 : 1.0;
    const double f = std::pow(kt2 / runtime.q2_node.front(), power);
    return zero[entry] * (1.0 - f) + bank[entry].front() * f;
  }
  if (kt2 >= runtime.q2_node.back()) {
    // Accept the rounding of the same upper bound through logarithms or squared loop momenta
    const double tolerance = 8.0 * std::numeric_limits<double>::epsilon() * kt2;
    if (kt2 - runtime.q2_node.back() <= tolerance) { return bank[entry].back(); }
    throw std::out_of_range(
        "InterpolateSpinBlock: kt2 exceeds NUMERICS_EIKONAL.MaxKT2");
  }
  const std::size_t harmonic = SpinHarmonic(entry);
  if (harmonic == 0) {
    return math::LinearInterpolateValidatedGrid(runtime.q2_node, bank[entry], kt2);
  }

  // Interpolate A_n(q)/q^n, which is regular in q^2 at the forward limit
  // [REFERENCE: Buttimore et al., hep-ph/9901339, Eq. (21)]
  const auto upper = std::lower_bound(runtime.q2_node.cbegin(), runtime.q2_node.cend(), kt2);
  std::size_t hi = static_cast<std::size_t>(upper - runtime.q2_node.cbegin());
  std::size_t lo = hi - 1;
  const double power = 0.5 * static_cast<double>(harmonic);
  if (!(runtime.q2_node[lo] > 0.0)) {
    // Extrapolate the regularized amplitude from the first two nonzero momentum nodes
    ++lo;
    ++hi;
  }
  const double fraction = (kt2 - runtime.q2_node[lo]) / (runtime.q2_node[hi] - runtime.q2_node[lo]);
  Matrix amplitude = bank[entry][lo] * ((1.0 - fraction) * std::pow(kt2 / runtime.q2_node[lo], power));
  amplitude.AddScaled(bank[entry][hi], fraction * std::pow(kt2 / runtime.q2_node[hi], power));
  return amplitude;
}

// Fourier-Bessel transform oriented impact-parameter spin-profile banks
// A_n(q) = 4 pi sqrt(lambda) (-i)^n integral b db J_n(bq) Gamma_n(b)
// The linear transform retains the complex profile sign and phase
template <typename Tables, std::size_t ProfileCount>
std::array<PairSpinBank, ProfileCount> MomentumTransformSpinBanks(
    const Tables &runtime,
    const std::array<const PairSpinTable *, ProfileCount> &profile,
    const MEikonalNumerics &numerics, const double bstep, const double sqrt_lambda,
    const double q) {
  for (const PairSpinTable *bank : profile) {
    if (bank == nullptr) {
      throw std::invalid_argument(
          "MomentumTransformSpinBanks: null profile bank");
    }
    for (const auto &entry : *bank) {
      if (entry.size() != runtime.b_node.size()) {
        throw std::invalid_argument(
            "MomentumTransformSpinBanks: inconsistent profile grid");
      }
    }
  }
  std::array<PairSpinBank, ProfileCount> integral;
  for (PairSpinBank &bank : integral) {
    for (Matrix &entry : bank) {
      entry = Matrix(runtime.D);
    }
  }
  for (const auto &ib : indices(runtime.b_node)) {
    const double b = runtime.b_node[ib];
    const double measure = numerics.logBT ? b * b : b;
    const double common =
        CompositeWeight(ib, numerics.NumberBT) * bstep * measure;
    const auto bessel = math::BesselJ012(b * q);
    for (const auto &family : indices(integral)) {
      for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
        integral[family][entry].AddScaled((*profile[family])[entry][ib],
                                          common * bessel[SpinHarmonic(entry)]);
      }
    }
  }
  for (PairSpinBank &bank : integral) {
    for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
      bank[entry] *=
          4.0 * math::PI * sqrt_lambda *
          HarmonicPhase(static_cast<int>(SpinHarmonic(entry)), -1);
    }
  }
  return integral;
}

// Build one physical two-proton state in the Good-Walker pair basis
template <typename Tables>
std::vector<double>
PairStateVector(const Tables &runtime, const GoodWalkerSpace &good_walker,
                const std::size_t f1, const std::size_t f2) {
  if (good_walker.ChannelCount() != runtime.N ||
      good_walker.PairDimension() != runtime.D) {
    throw std::invalid_argument(
        "PairStateVector: Good-Walker dimensions disagree");
  }
  if (f1 >= runtime.N || f2 >= runtime.N) {
    throw std::out_of_range(
        "PairStateVector: physical-state index is out of range");
  }
  const auto &mixing = good_walker.MixingMatrix();
  return gra::KroneckerProduct(mixing.Row(f1), mixing.Row(f2));
}

// Project one pair-basis amplitude between physical Good-Walker states
// output = <f1 f2|A|i1 i2> after the U x U rotation
template <typename Tables>
Complex ProjectAmplitude(const Tables &runtime, const Matrix &amplitude,
                         const GoodWalkerSpace &good_walker,
                         const std::size_t f1, const std::size_t f2,
                         const std::size_t i1 = 0, const std::size_t i2 = 0) {
  const auto initial = PairStateVector(runtime, good_walker, i1, i2);
  const auto final = PairStateVector(runtime, good_walker, f1, f2);
  return amplitude.MatrixElement(final, initial);
}

// Compute the six named definitions from one complete reference-plane matrix
ElasticHelicityAmplitudes
DefiningHelicityAmplitudes(const ProtonHelicityMatrix &matrix) {
  const auto value = [&matrix](const ElasticHelicityComponent component) {
    const auto transition = ElasticHelicityReferenceTransition(component);
    return transition.sign * matrix[ElasticHelicityReferenceIndex(component)];
  };
  ElasticHelicityAmplitudes out;
  out.phi1 = value(ElasticHelicityComponent::Phi1);
  out.phi2 = value(ElasticHelicityComponent::Phi2);
  out.phi3 = value(ElasticHelicityComponent::Phi3);
  out.phi4 = value(ElasticHelicityComponent::Phi4);
  out.phi5 = value(ElasticHelicityComponent::Phi5);
  out.phi5_first = value(ElasticHelicityComponent::Phi5First);
  return out;
}

// Compute the six named pair-space definitions from all oriented blocks
MEikonalMatrix::PairHelicityBank
DefiningPairHelicityBank(const PairSpinBank &spin) {
  const auto value = [&spin](const ElasticHelicityComponent component) {
    const auto transition = ElasticHelicityReferenceTransition(component);
    return spin[ElasticHelicityReferenceIndex(component)] *
           Complex(transition.sign);
  };
  MEikonalMatrix::PairHelicityBank out;
  out.phi1 = value(ElasticHelicityComponent::Phi1);
  out.phi2 = value(ElasticHelicityComponent::Phi2);
  out.phi3 = value(ElasticHelicityComponent::Phi3);
  out.phi4 = value(ElasticHelicityComponent::Phi4);
  out.phi5 = value(ElasticHelicityComponent::Phi5);
  out.phi5_first = value(ElasticHelicityComponent::Phi5First);
  return out;
}

// Interpolate all oriented spin blocks at one momentum transfer
template <typename Tables>
PairSpinBank InterpolateSpinBank(const Tables &runtime,
                                 const PairSpinTable &bank,
                                 const PairSpinBank &zero, const double kt2) {
  PairSpinBank out;
  for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
    out[entry] = InterpolateSpinBlock(runtime, bank, zero, entry, kt2);
  }
  return out;
}

// Project every oriented spin block between two physical Good Walker pairs
template <typename Tables>
ProtonHelicityMatrix
ProjectSpinBank(const Tables &runtime, const PairSpinBank &bank,
                const GoodWalkerSpace &good_walker, const std::size_t f1,
                const std::size_t f2, const std::size_t i1,
                const std::size_t i2, const double azimuth) {
  const auto initial = PairStateVector(runtime, good_walker, i1, i2);
  const auto final = PairStateVector(runtime, good_walker, f1, f2);
  ProtonHelicityMatrix reference{};
  for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
    reference[entry] = bank[entry].MatrixElement(final, initial);
  }
  return RotateProtonHelicityMatrix(reference, azimuth);
}

// Build one nonflip or Pauli flip proton Regge transition vertex
Matrix BuildVertexMatrix(const SoftModel &model, const SoftExchangeId exchange,
                         const double t, const std::size_t n,
                         const bool helicity_flip = false) {
  const auto residue = helicity_flip
                           ? model.HelicityFlipResidueMatrix(exchange, t)
                           : model.ResidueMatrix(exchange, t);
  if (residue.size_row() != n || residue.size_col() != n) {
    throw std::logic_error(
        "BuildVertexMatrix: exchange residue dimension mismatch");
  }
  return residue.Transform(
      [](const double value) { return Complex(value, 0.0); });
}

// Build crossing-separated Born helicity tables on the momentum grid
struct BornTables {
  std::vector<double> q;
  std::array<std::vector<Matrix>, kNHelicity> even;
  std::array<std::vector<Matrix>, kNHelicity> odd;
  std::vector<Matrix> even_phi5_first;
  std::vector<Matrix> odd_phi5_first;
  std::array<std::vector<Matrix>, kNHelicity> cut;
  std::vector<Matrix> cut_phi5_first;
  std::array<std::vector<Matrix>, kNHelicity> screening;
  std::vector<Matrix> screening_phi5_first;
};

// Construct crossing-separated Born helicity tables
BornTables BuildBornTables(const SoftModel &model, const double s,
                           const MEikonalNumerics &numerics,
                           const std::size_t n,
                           const double screening_odd_sign) {
  const std::size_t points = numerics.FBIntegralN + 1;
  const std::size_t d = model.GoodWalker().PairDimension();
  if (model.GoodWalker().ChannelCount() != n) {
    throw std::invalid_argument(
        "BuildBornTables: Good-Walker channel count disagrees");
  }
  BornTables table;
  table.q.resize(points);
  for (auto &entry : table.even) {
    entry.assign(points, Matrix(d));
  }
  for (auto &entry : table.odd) {
    entry.assign(points, Matrix(d));
  }
  table.even_phi5_first.assign(points, Matrix(d));
  table.odd_phi5_first.assign(points, Matrix(d));
  for (auto &entry : table.cut) {
    entry.assign(points, Matrix(d));
  }
  table.cut_phi5_first.assign(points, Matrix(d));
  for (auto &entry : table.screening) {
    entry.assign(points, Matrix(d));
  }
  table.screening_phi5_first.assign(points, Matrix(d));
  const double step = (numerics.FBIntegralMaxKT - numerics.FBIntegralMinKT) /
                      static_cast<double>(numerics.FBIntegralN);
  const auto &settings = model.Eikonal();
  const double mass = settings.HelicityMassScale();

  ParallelFor(
      points,
      [&](const std::size_t iq) {
        const double q = numerics.FBIntegralMinKT + iq * step;
        const double t = -q * q;
        table.q[iq] = q;
        for (std::size_t index = 0; index < model.ExchangeCount(); ++index) {
          const SoftExchangeId exchange(index);
          const auto &parameter = model.Exchange(exchange);
          if (!parameter.Enabled()) {
            continue;
          }
          const double alpha = model.Alpha(exchange, t);
          const Complex eta =
              regge::EtaFactor(alpha, parameter.Alpha0(), parameter.Signature(),
                               parameter.Eta());
          const Complex common = static_cast<double>(parameter.ResidueSign()) *
                                 eta * std::pow(s, alpha);
          const Matrix v0 = BuildVertexMatrix(model, exchange, t, n);
          const Matrix central = v0.Kronecker(v0) * common;
          std::array<Matrix, kNHelicity> term = {central, Matrix(d), central,
                                                 Matrix(d), Matrix(d)};
          Matrix phi5_first(d);
          if (settings.HelicityEnabled()) {
            const Matrix vf = BuildVertexMatrix(model, exchange, t, n, true);
            const Matrix double_flip =
                vf.Kronecker(vf) * (q * q * common / (4.0 * mass * mass));
            term[kPhi2] = -double_flip;
            term[kPhi4] = double_flip;
            phi5_first = vf.Kronecker(v0) * (q * common / (2.0 * mass));
            term[kPhi5] = v0.Kronecker(vf) * (q * common / (2.0 * mass));
          }
          auto &target =
              parameter.CrossingParity() == 1 ? table.even : table.odd;
          for (std::size_t h = 0; h < kNHelicity; ++h) {
            target[h][iq] += term[h];
          }
          auto &target_phi5_first = parameter.CrossingParity() == 1
                                        ? table.even_phi5_first
                                        : table.odd_phi5_first;
          target_phi5_first[iq] += phi5_first;
          if (parameter.Role() == SoftExchangeRole::Pomeron) {
            for (std::size_t h = 0; h < kNHelicity; ++h) {
              table.cut[h][iq] += term[h];
            }
            table.cut_phi5_first[iq] += phi5_first;
          }
          const auto &screening_exchanges = settings.ScreeningExchanges();
          if (std::find(screening_exchanges.begin(), screening_exchanges.end(),
                        exchange) != screening_exchanges.end()) {
            const double screening_sign =
                parameter.CrossingParity() == 1 ? 1.0 : screening_odd_sign;
            for (std::size_t h = 0; h < kNHelicity; ++h) {
              table.screening[h][iq].AddScaled(term[h], screening_sign);
            }
            table.screening_phi5_first[iq].AddScaled(phi5_first,
                                                     screening_sign);
          }
        }
      },
      "Born momentum table");
  return table;
}

// Compute the pair-space matrix for one named elastic helicity component
const Matrix &SpinComponentMatrix(const std::array<Matrix, kNHelicity> &chi,
                                  const Matrix &phi5_first,
                                  const Matrix &phi5_second,
                                  const ElasticHelicityComponent component) {
  switch (component) {
  case ElasticHelicityComponent::Phi1:
    return chi[kPhi1];
  case ElasticHelicityComponent::Phi2:
    return chi[kPhi2];
  case ElasticHelicityComponent::Phi3:
    return chi[kPhi3];
  case ElasticHelicityComponent::Phi4:
    return chi[kPhi4];
  case ElasticHelicityComponent::Phi5:
    return phi5_second;
  case ElasticHelicityComponent::Phi5First:
    return phi5_first;
  }
  throw std::invalid_argument("SpinComponentMatrix: unknown component");
}

// Embed helicity matrices with signed collider helicity reciprocity
Matrix BuildSpinMatrix(const std::array<Matrix, kNHelicity> &chi,
                       const Matrix &phi5_first, const Matrix &phi5_second) {
  const std::size_t d = chi[kPhi1].size_row();
  if (d > std::numeric_limits<std::size_t>::max() / 4) {
    throw std::length_error("BuildSpinMatrix: spin dimension overflow");
  }
  Matrix spin_matrix(4 * d);
  constexpr auto transitions = CanonicalProtonHelicityTransitions();
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const auto &transition =
          transitions[gra::spin::PairHelicityMatrixIndex(row, col)];
      spin_matrix.SetBlock(row, col,
                           SpinComponentMatrix(chi, phi5_first, phi5_second,
                                               transition.component),
                           transition.sign);
    }
  }
  return spin_matrix;
}

// Compute the selected matrix eikonal S matrix from one chi operator
//
// exp:   S = exp(i chi)
// q_exp: S = [I + (1-q) i chi]^(1/(1-q)), principal matrix branch
Matrix MatrixEikonalSMatrix(const Matrix &chi,
                            const SoftEikonalSettings &settings) {
  const Matrix argument = zi * chi;

  if (settings.Unitarization() == SoftUnitarization::Exponential) {
    return argument.Exp();
  }

  if (settings.Unitarization() == SoftUnitarization::QExponential) {
    const double q = settings.Q();
    if (!std::isfinite(q) || q <= 0.0) {
      throw std::invalid_argument(
          "MatrixEikonalSMatrix: q must be finite and positive");
    }
    if (std::abs(q - 1.0) <= 1.0e-12) {
      return argument.Exp();
    }

    const Matrix base =
        Matrix::IdentityMatrix(chi.size_row()) + argument * (1.0 - q);
    return base.PrincipalPower(1.0 / (1.0 - q));
  }

  throw std::logic_error("MatrixEikonalSMatrix: unknown unitarization");
}

// Split one coupled spin and Good Walker operator into oriented spin blocks
PairSpinBank SplitSpinBlocks(const Matrix &spin, const std::size_t d) {
  if (spin.size_row() != 4 * d || spin.size_col() != 4 * d) {
    throw std::invalid_argument("SplitSpinBlocks: spin dimensions disagree");
  }
  PairSpinBank out;
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      out[gra::spin::PairHelicityMatrixIndex(row, col)] =
          spin.Submatrix(row * d, col * d, d, d);
    }
  }
  return out;
}

// Form the complex Gamma = i(I-S) amplitude before probability projection
std::pair<PairSpinBank, double>
UnitarizeHelicityProfile(const std::array<Matrix, kNHelicity> &chi,
                         const Matrix &phi5_first,
                         const SoftEikonalSettings &settings) {
  const std::size_t d = chi[kPhi1].size_row();
  Matrix smatrix;
  PairSpinBank blocks;
  for (Matrix &block : blocks) {
    block = Matrix(d);
  }
  if (settings.HelicityEnabled()) {
    smatrix = MatrixEikonalSMatrix(BuildSpinMatrix(chi, phi5_first, chi[kPhi5]),
                                   settings);
    blocks = SplitSpinBlocks(zi * (Matrix::IdentityMatrix(4 * d) - smatrix), d);
  } else {
    smatrix = MatrixEikonalSMatrix(chi[kPhi1], settings);
    const Matrix amplitude = zi * (Matrix::IdentityMatrix(d) - smatrix);
    for (const std::size_t diagonal : {0U, 5U, 10U, 15U}) {
      blocks[diagonal] = amplitude;
    }
  }
  return {blocks, smatrix.MaxSingularValue()};
}

// Compute the C-odd sign selected by the incoming beam pair
double BeamOddSign(const std::vector<MParticle> &initialstate) {
  if (initialstate.size() != 2) {
    throw std::invalid_argument("BeamOddSign: exactly two beams are required");
  }
  if (initialstate[0].pdg == initialstate[1].pdg) {
    return 1.0;
  }
  if (initialstate[0].pdg == -initialstate[1].pdg) {
    return -1.0;
  }
  throw std::invalid_argument("BeamOddSign: unsupported incoming beam pair");
}

// Hold profile matrices separated by crossing parity for momentum transforms
struct ImpactParameterProfiles {
  PairSpinTable profile;
  PairSpinTable even;
  PairSpinTable odd;
  PairSpinTable screening;
};

// Fourier-Bessel transform and eikonalize the Born matrices on the b grid
template <typename Tables>
ImpactParameterProfiles BuildImpactParameterTables(
    Tables &runtime, const BornTables &born, const MEikonalNumerics &numerics,
    const std::vector<MParticle> &initialstate, const double sqrt_lambda,
    const SoftEikonalSettings &settings) {
  const std::size_t b_points = numerics.NumberBT + 1;
  const std::size_t d = runtime.D;
  const std::size_t cut_dimension = settings.HelicityEnabled() ? 4 * d : d;
  const double coordinate_min =
      numerics.logBT ? std::log(numerics.MinBT) : numerics.MinBT;
  const double coordinate_max =
      numerics.logBT ? std::log(numerics.MaxBT) : numerics.MaxBT;
  const double bstep = (coordinate_max - coordinate_min) / numerics.NumberBT;
  const double qstep = (numerics.FBIntegralMaxKT - numerics.FBIntegralMinKT) /
                       numerics.FBIntegralN;
  const double inverse_norm = 1.0 / (4.0 * math::PI * sqrt_lambda);
  const double beam_odd_sign = BeamOddSign(initialstate);

  runtime.b_node.resize(b_points);
  runtime.chi_even_b.assign(b_points, Matrix(d));
  runtime.chi_odd_b.assign(b_points, Matrix(d));
  runtime.chi_cut_b.assign(b_points, Matrix(cut_dimension));

  ImpactParameterProfiles tables;
  for (auto *bank : {&tables.profile, &tables.even, &tables.odd}) {
    for (auto &entry : *bank) {
      entry.assign(b_points, Matrix(d));
    }
  }
  for (auto &entry : tables.screening) {
    entry.assign(b_points, Matrix(d));
  }
  std::vector<double> singular_value(b_points);

  ParallelFor(
      b_points,
      [&](const std::size_t ib) {
        const double coordinate = coordinate_min + ib * bstep;
        const double b = numerics.logBT ? std::exp(coordinate) : coordinate;
        runtime.b_node[ib] = b;
        std::array<Matrix, kNHelicity> chi_even = {
            Matrix(d), Matrix(d), Matrix(d), Matrix(d), Matrix(d)};
        std::array<Matrix, kNHelicity> chi_odd = chi_even;
        std::array<Matrix, kNHelicity> chi_cut = chi_even;
        std::array<Matrix, kNHelicity> chi_screening = chi_even;
        Matrix chi_even_phi5_first(d);
        Matrix chi_odd_phi5_first(d);
        Matrix chi_cut_phi5_first(d);
        Matrix chi_screening_phi5_first(d);

        for (const auto &iq : indices(born.q)) {
          const double weight =
              CompositeWeight(iq, numerics.FBIntegralN) * qstep * born.q[iq];
          const auto bessel = math::BesselJ012(b * born.q[iq]);
          for (std::size_t helicity = 0; helicity < kNHelicity; ++helicity) {
            const Complex kernel = weight * bessel[kHarmonic[helicity]] *
                                   HarmonicPhase(kHarmonic[helicity], 1) *
                                   inverse_norm;
            chi_even[helicity].AddScaled(born.even[helicity][iq], kernel);
            chi_odd[helicity].AddScaled(born.odd[helicity][iq], kernel);
            chi_cut[helicity].AddScaled(born.cut[helicity][iq], kernel);
            chi_screening[helicity].AddScaled(born.screening[helicity][iq],
                                              kernel);
          }
          const Complex phi5_kernel = weight * bessel[kHarmonic[kPhi5]] *
                                      HarmonicPhase(kHarmonic[kPhi5], 1) *
                                      inverse_norm;
          chi_even_phi5_first.AddScaled(born.even_phi5_first[iq], phi5_kernel);
          chi_odd_phi5_first.AddScaled(born.odd_phi5_first[iq], phi5_kernel);
          chi_cut_phi5_first.AddScaled(born.cut_phi5_first[iq], phi5_kernel);
          chi_screening_phi5_first.AddScaled(born.screening_phi5_first[iq],
                                             phi5_kernel);
        }
        runtime.chi_even_b[ib] = chi_even[kPhi1];
        runtime.chi_odd_b[ib] = chi_odd[kPhi1];
        runtime.chi_cut_b[ib] =
            settings.HelicityEnabled()
                ? BuildSpinMatrix(chi_cut, chi_cut_phi5_first, chi_cut[kPhi5])
                : chi_cut[kPhi1];

        const auto unitarize = [&](const double sign) {
          std::array<Matrix, kNHelicity> chi = chi_even;
          for (std::size_t helicity = 0; helicity < kNHelicity; ++helicity) {
            chi[helicity].AddScaled(chi_odd[helicity], sign);
          }
          const Matrix chi_phi5_first =
              chi_even_phi5_first + chi_odd_phi5_first * sign;
          return UnitarizeHelicityProfile(chi, chi_phi5_first, settings);
        };

        const auto pp = unitarize(1.0);
        const auto pbar = unitarize(-1.0);
        for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
          tables.even[entry][ib] = (pp.first[entry] + pbar.first[entry]) * 0.5;
          tables.odd[entry][ib] = (pp.first[entry] - pbar.first[entry]) * 0.5;
          tables.profile[entry][ib] =
              tables.even[entry][ib] + tables.odd[entry][ib] * beam_odd_sign;
        }
        const auto screening = UnitarizeHelicityProfile(
            chi_screening, chi_screening_phi5_first, settings);
        for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
          tables.screening[entry][ib] = screening.first[entry];
        }
        singular_value[ib] =
            std::max({pp.second, pbar.second, screening.second});
      },
      "Impact-parameter transform and unitarization");
  runtime.max_singular_value =
      *std::max_element(singular_value.begin(), singular_value.end());
  return tables;
}

// Transform eikonalized impact-parameter profiles to momentum tables
template <typename Tables>
void BuildMomentumTables(Tables &runtime,
                         const ImpactParameterProfiles &profile,
                         const MEikonalNumerics &numerics, const double sqrt_lambda) {
  const std::size_t q_points = numerics.NumberKT2 + 1;
  const std::size_t d = runtime.D;
  const double coordinate_min =
      numerics.logKT2 ? std::log(numerics.MinKT2) : numerics.MinKT2;
  const double coordinate_max =
      numerics.logKT2 ? std::log(numerics.MaxKT2) : numerics.MaxKT2;
  const double q2step = (coordinate_max - coordinate_min) / numerics.NumberKT2;
  const double bstep =
      ((numerics.logBT ? std::log(numerics.MaxBT) : numerics.MaxBT) -
       (numerics.logBT ? std::log(numerics.MinBT) : numerics.MinBT)) /
      numerics.NumberBT;

  runtime.q2_node.resize(q_points);
  for (auto *bank : {&runtime.amplitude_spin_q, &runtime.crossing_even_spin_q,
                     &runtime.crossing_odd_spin_q, &runtime.screening_spin_q}) {
    for (auto &entry : *bank) {
      entry.assign(q_points, Matrix(d));
    }
  }
  const std::array<const PairSpinTable *, 4> profiles = {
      &profile.profile, &profile.even, &profile.odd, &profile.screening};
  const auto zero = MomentumTransformSpinBanks(runtime, profiles, numerics, bstep, sqrt_lambda, 0.0);
  runtime.amplitude_spin_zero = zero[0];
  runtime.crossing_even_spin_zero = zero[1];
  runtime.crossing_odd_spin_zero = zero[2];
  runtime.screening_spin_zero = zero[3];
  ParallelFor(
      q_points,
      [&](const std::size_t iq) {
        const double q2 = numerics.logKT2
                              ? std::exp(coordinate_min + iq * q2step)
                              : coordinate_min + iq * q2step;
        // Clamp only the validated q2 = 0 table boundary against roundoff
        const double q = std::sqrt(std::max(0.0, q2));
        runtime.q2_node[iq] = q2;
        const auto transformed = MomentumTransformSpinBanks(
            runtime, profiles, numerics, bstep, sqrt_lambda, q);
        for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
          runtime.amplitude_spin_q[entry][iq] = transformed[0][entry];
          runtime.crossing_even_spin_q[entry][iq] = transformed[1][entry];
          runtime.crossing_odd_spin_q[entry][iq] = transformed[2][entry];
          runtime.screening_spin_q[entry][iq] = transformed[3][entry];
        }
      },
      "Momentum-space transform");
}

// Select one named pair helicity matrix for const or mutable access
template <typename Bank>
decltype(auto) PairHelicityComponent(
    Bank &bank, const ElasticHelicityComponent component) {
  switch (component) {
  case ElasticHelicityComponent::Phi1:
    return (bank.phi1);
  case ElasticHelicityComponent::Phi2:
    return (bank.phi2);
  case ElasticHelicityComponent::Phi3:
    return (bank.phi3);
  case ElasticHelicityComponent::Phi4:
    return (bank.phi4);
  case ElasticHelicityComponent::Phi5:
    return (bank.phi5);
  case ElasticHelicityComponent::Phi5First:
    return (bank.phi5_first);
  }
  throw std::logic_error("PairHelicityComponent: unknown component");
}

} // namespace

// Compute one immutable named pair-space helicity matrix
const MEikonalMatrix::Matrix &MEikonalMatrix::PairHelicityBank::Component(
    const ElasticHelicityComponent component) const {
  return PairHelicityComponent(*this, component);
}

// Compute one mutable named pair-space helicity matrix
MEikonalMatrix::Matrix &MEikonalMatrix::PairHelicityBank::Component(
    const ElasticHelicityComponent component) {
  return PairHelicityComponent(*this, component);
}

// Construct the coupled Good-Walker and helicity matrix eikonal tables
std::shared_ptr<const MEikonalMatrix>
MEikonalMatrix::Build(const SoftModelPtr &model, const double s,
                      const std::vector<MParticle> &initialstate,
                      const MEikonalNumerics &numerics,
                      const bool build_amplitude) {
  if (model == nullptr) {
    throw std::invalid_argument("MEikonalMatrix::Build: null SOFT model");
  }
  const auto &good_walker = model->GoodWalker();
  const auto &model_mixing = good_walker.MixingMatrix();
  const std::size_t n = good_walker.ChannelCount();
  const std::size_t d = good_walker.PairDimension();

  if (!std::isfinite(s) || s <= 0.0 || initialstate.size() != 2 || n == 0 ||
      model_mixing.size_row() != n || model_mixing.size_col() != n ||
      numerics.NumberBT < 2 || numerics.NumberBT % 2 != 0 ||
      numerics.NumberKT2 < 2 || numerics.FBIntegralN < 2 ||
      !std::isfinite(numerics.MinBT) ||
      !std::isfinite(numerics.unitarity_tolerance) ||
      numerics.unitarity_tolerance <= 0.0 || !std::isfinite(numerics.MaxBT) ||
      numerics.MinBT < 0.0 || numerics.MaxBT <= numerics.MinBT ||
      (numerics.logBT && numerics.MinBT <= 0.0) ||
      !std::isfinite(numerics.MinKT2) || !std::isfinite(numerics.MaxKT2) ||
      numerics.MinKT2 < 0.0 || numerics.MaxKT2 <= numerics.MinKT2 ||
      (numerics.logKT2 && numerics.MinKT2 <= 0.0) ||
      !std::isfinite(numerics.FBIntegralMinKT) ||
      !std::isfinite(numerics.FBIntegralMaxKT) ||
      numerics.FBIntegralMinKT < 0.0 ||
      numerics.FBIntegralMaxKT <= numerics.FBIntegralMinKT) {
    throw std::invalid_argument(
        "MEikonalMatrix::Build: invalid matrix eikonal configuration");
  }
  const double beta = kinematics::beta12(s, initialstate[0].mass, initialstate[1].mass);
  if (!std::isfinite(beta) || !(beta > 0.0)) {
    throw std::invalid_argument("MEikonalMatrix::Build: initial state is below threshold");
  }
  // sqrt(lambda(s,m1^2,m2^2)) is the exact two-body amplitude flux scale
  const double sqrt_lambda = s * beta;
  const auto &settings = model->Eikonal();
  const std::size_t bp = numerics.NumberBT + 1;
  const std::size_t qp = numerics.NumberKT2 + 1;
  const std::size_t cut_dimension = settings.HelicityEnabled() ? 4 * d : d;
  const std::string cache_key =
      MatrixCacheKey(s, initialstate, numerics, *model);
  const std::string runtime_fingerprint =
      std::to_string(kCacheVersion) + ';' + cache_key;
  std::unique_ptr<eikonal::MCacheLock> cache_lock;
  std::string cache_filename;

  // Load one cache candidate into mutable storage before const publication
  const auto read_cache = [&](const std::string &filename,
                              const std::string &key) {
    auto runtime = std::shared_ptr<MEikonalMatrix>(new MEikonalMatrix(model));
    runtime->unitarity_tolerance = numerics.unitarity_tolerance;
    if (!ReadMatrixJsonCache(filename, key, n, d, bp, qp, cut_dimension,
                             build_amplitude, runtime->tables)) {
      runtime.reset();
      return runtime;
    }
    if (runtime->tables.U != model_mixing.Flatten()) {
      runtime.reset();
      return runtime;
    }
    CheckContractivity(runtime->tables.max_singular_value, numerics);
    runtime->runtime_fingerprint = runtime_fingerprint;
    return runtime;
  };

  try {
    cache_filename = MatrixCacheFilename(cache_key, initialstate, n);
    aux::CreateDirectory(std::filesystem::path(cache_filename).parent_path().string());
    if (const auto cached = read_cache(cache_filename, cache_key)) {
      std::cout << "Loaded matrix eikonal cache: " << cache_filename
                << std::endl;
      return cached;
    }
    cache_lock =
        std::make_unique<eikonal::MCacheLock>(cache_filename + ".lock");
    if (!cache_lock->Acquired()) {
      std::cerr << "WARNING: MEikonalMatrix cache lock unavailable, rebuilding "
                << "without disk cache" << std::endl;
      cache_lock.reset();
    } else if (const auto cached = read_cache(cache_filename, cache_key)) {
      std::cout << "Loaded matrix eikonal cache: " << cache_filename
                << std::endl;
      return cached;
    }
  } catch (const std::exception &error) {
    std::cerr << "WARNING: MEikonalMatrix cache unavailable: " << error.what()
              << std::endl;
    cache_lock.reset();
    cache_filename.clear();
  }

  // Publish only complete tables while the cache lock is held
  const auto publish_cache =
      [&](const std::shared_ptr<MEikonalMatrix> &runtime) {
        if (!cache_lock) {
          return;
        }
        try {
          WriteMatrixJsonCache(cache_filename, cache_key, runtime->tables,
                               build_amplitude);
          std::cout << "Saved matrix eikonal cache: " << cache_filename
                    << std::endl;
        } catch (const std::exception &error) {
          std::cerr << "WARNING: MEikonalMatrix cache not saved: "
                    << error.what() << std::endl;
        }
      };

  // Store the physical proton rotation in the Good-Walker pair basis
  auto out = std::shared_ptr<MEikonalMatrix>(new MEikonalMatrix(model));
  out->unitarity_tolerance = numerics.unitarity_tolerance;
  out->runtime_fingerprint = runtime_fingerprint;
  out->tables.N = n;
  out->tables.D = d;
  out->tables.U = model_mixing.Flatten();
  std::cout << "Building full " << d << " x " << d << " pair eikonal Matrix"
            << std::endl;
  MTimer timer(true);

  // Construct crossing separated Regge Born matrices on the transverse grid
  const BornTables born =
      BuildBornTables(*model, s, numerics, n, BeamOddSign(initialstate));

  // Build b-space profiles, eikonalize the spin matrix, and retain cut opacity
  const ImpactParameterProfiles profile = BuildImpactParameterTables(
      out->tables, born, numerics, initialstate, sqrt_lambda, settings);
  ValidateInterpolationNodes(out->tables, false);

  // Enforce the configured S-matrix contractivity tolerance
  CheckContractivity(out->tables.max_singular_value, numerics);

  // Retain impact-parameter amplitudes for unitarity-consistent cross sections
  out->tables.amplitude_spin_b = profile.profile;
  out->tables.screening_spin_b = profile.screening;
  if (!build_amplitude) {
    publish_cache(out);
    return out;
  }

  // Transform physical and crossing-separated helicity profiles to q^2 tables
  BuildMomentumTables(out->tables, profile, numerics, sqrt_lambda);
  ValidateInterpolationNodes(out->tables, true);
  std::cout << "Matrix/helicity Eikonal tables built in " << timer.ElapsedSec()
            << " s" << std::endl;
  publish_cache(out);
  return out;
}

// Compute separated crossing-even and crossing-odd diagonal profiles
std::pair<Complex, Complex>
MEikonalMatrix::DensityProfiles(const double bt, const std::size_t i,
                                const std::size_t k) const {
  if (i >= tables.N || k >= tables.N)
    throw std::out_of_range(
        "MEikonalMatrix::DensityProfiles: channel index out of range");
  const std::size_t pair = soft_model->GoodWalker().PairIndex(i, k);
  return {math::LinearInterpolateValidatedGrid(tables.b_node, tables.chi_even_b,
                                               bt)(pair, pair),
          math::LinearInterpolateValidatedGrid(tables.b_node, tables.chi_odd_b,
                                               bt)(pair, pair)};
}

// Compute the spin-averaged diagonal eigenchannel amplitude
Complex MEikonalMatrix::ScreeningAmplitude(const double kt2,
                                           const std::size_t i,
                                           const std::size_t k) const {
  if (i >= tables.N || k >= tables.N)
    throw std::out_of_range(
        "MEikonalMatrix::ScreeningAmplitude: channel index out of range");
  const std::size_t pair = soft_model->GoodWalker().PairIndex(i, k);
  Complex amplitude = 0.0;
  for (const std::size_t diagonal : {0U, 5U, 10U, 15U}) {
    amplitude += InterpolateSpinBlock(tables, tables.amplitude_spin_q,
                                      tables.amplitude_spin_zero, diagonal,
                                      kt2)(pair, pair);
  }
  return 0.25 * amplitude;
}

// Compute physical helicity amplitudes for one Good-Walker final state
ElasticHelicityAmplitudes
MEikonalMatrix::HelicityAmplitudes(const double kt2, const std::size_t f1,
                                   const std::size_t f2) const {
  return DefiningHelicityAmplitudes(HelicityMatrix(kt2, f1, f2));
}

// Compute the exact physical spin matrix for one Good Walker transition
ProtonHelicityMatrix MEikonalMatrix::HelicityMatrix(
    const double kt2, const std::size_t f1, const std::size_t f2,
    const std::size_t i1, const std::size_t i2, const double azimuth) const {
  return ProjectSpinBank(tables,
                         InterpolateSpinBank(tables, tables.amplitude_spin_q,
                                             tables.amplitude_spin_zero, kt2),
                         soft_model->GoodWalker(), f1, f2, i1, i2, azimuth);
}

// Compute exclusive amplitudes built from the selected screening exchanges
ElasticHelicityAmplitudes MEikonalMatrix::ScreeningHelicityAmplitudes(
    const double kt2, const std::size_t f1, const std::size_t f2,
    const std::size_t i1, const std::size_t i2) const {
  return DefiningHelicityAmplitudes(
      ScreeningHelicityMatrix(kt2, f1, f2, i1, i2));
}

// Compute the exact selected-exchange screening spin matrix
ProtonHelicityMatrix MEikonalMatrix::ScreeningHelicityMatrix(
    const double kt2, const std::size_t f1, const std::size_t f2,
    const std::size_t i1, const std::size_t i2, const double azimuth) const {
  return ProjectSpinBank(tables,
                         InterpolateSpinBank(tables, tables.screening_spin_q,
                                             tables.screening_spin_zero, kt2),
                         soft_model->GoodWalker(), f1, f2, i1, i2, azimuth);
}

// Compute the six named defining entries of the radial screening operator
MEikonalMatrix::PairHelicityBank
MEikonalMatrix::PairScreeningHelicityBank(const double kt2) const {
  return DefiningPairHelicityBank(PairScreeningSpinBank(kt2));
}

// Compute every radial pair-space screening spin transition
MEikonalMatrix::PairSpinBank
MEikonalMatrix::PairScreeningSpinBank(const double kt2) const {
  return InterpolateSpinBank(tables, tables.screening_spin_q,
                             tables.screening_spin_zero, kt2);
}

// Compute physical impact-parameter helicity amplitudes for one final state
ElasticHelicityAmplitudes
MEikonalMatrix::ImpactHelicityAmplitudes(const double bt, const std::size_t f1,
                                         const std::size_t f2) const {
  return DefiningHelicityAmplitudes(ImpactHelicityMatrix(bt, f1, f2));
}

// Compute the exact impact-parameter spin matrix for one physical transition
ProtonHelicityMatrix MEikonalMatrix::ImpactHelicityMatrix(
    const double bt, const std::size_t f1, const std::size_t f2,
    const std::size_t i1, const std::size_t i2, const double azimuth) const {
  if (!std::isfinite(bt) || bt < 0.0) {
    throw std::invalid_argument(
        "MEikonalMatrix::ImpactHelicityMatrix: negative impact parameter");
  }
  PairSpinBank bank;
  for (std::size_t entry = 0; entry < kNSpinEntries; ++entry) {
    bank[entry] = math::LinearInterpolateValidatedGrid(
        tables.b_node, tables.amplitude_spin_b[entry], bt);
  }
  return ProjectSpinBank(tables, bank, soft_model->GoodWalker(), f1, f2, i1, i2,
                         azimuth);
}

// Compute one physical impact-parameter helicity S-matrix transition
MEikonalMatrix::Matrix
MEikonalMatrix::PhysicalImpactSMatrix(const double bt, const std::size_t f1,
                                      const std::size_t f2,
                                      const double azimuth) const {
  if (tables.b_node.empty() || bt < tables.b_node.front() ||
      bt > tables.b_node.back()) {
    throw std::out_of_range(
        "MEikonalMatrix::PhysicalImpactSMatrix: impact parameter is outside "
        "the hadronic grid");
  }
  const ProtonHelicityMatrix amplitude =
      ImpactHelicityMatrix(bt, f1, f2, 0, 0, azimuth);
  Matrix smatrix(4);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      smatrix(row, col) =
          math::zi * amplitude[gra::spin::PairHelicityMatrixIndex(row, col)];
    }
  }
  if (f1 == 0 && f2 == 0) {
    smatrix += Matrix::IdentityMatrix(4);
  }
  return smatrix;
}

// Compute separately screened crossing helicity amplitudes
ElasticCrossingHelicityAmplitudes MEikonalMatrix::CrossingHelicityAmplitudes(
    const double kt2, const std::size_t f1, const std::size_t f2) const {
  ElasticCrossingHelicityAmplitudes out;
  const ProtonHelicityMatrix even =
      ProjectSpinBank(tables,
                      InterpolateSpinBank(tables, tables.crossing_even_spin_q,
                                          tables.crossing_even_spin_zero, kt2),
                      soft_model->GoodWalker(), f1, f2, 0, 0, 0.0);
  const ProtonHelicityMatrix odd =
      ProjectSpinBank(tables,
                      InterpolateSpinBank(tables, tables.crossing_odd_spin_q,
                                          tables.crossing_odd_spin_zero, kt2),
                      soft_model->GoodWalker(), f1, f2, 0, 0, 0.0);
  out.even = DefiningHelicityAmplitudes(even);
  out.odd = DefiningHelicityAmplitudes(odd);
  return out;
}

// Compute the Pomeron-only pre-unitarization chi operator at one impact
// parameter
Matrix MEikonalMatrix::CutOpacity(const double bt) const {
  if (!std::isfinite(bt) || bt < 0.0) {
    throw std::invalid_argument(
        "MEikonalMatrix::CutOpacity: invalid impact parameter");
  }
  return math::LinearInterpolateValidatedGrid(tables.b_node, tables.chi_cut_b,
                                              bt);
}

// Construct the effective Pomeron cut-opacity spectrum from the actual
// unitarized S matrix.  For a general non-normal chi, -i(chi-chi^dagger) is
// not equal to -Log(S^dagger S).  Defining
//
//   S_P^dagger S_P |a> = exp(-lambda_a) |a>
//
// makes the no-cut probability exact for the chosen q-eikonal.  The subsequent
// Poisson distribution in lambda_a is the explicit effective cut model
MEikonalMatrix::CutOpacitySpectrum
MEikonalMatrix::CutSpectrum(const double bt) const {
  const Matrix chi = CutOpacity(bt);
  const auto &settings = soft_model->Eikonal();
  const Matrix smatrix = MatrixEikonalSMatrix(chi, settings);
  // Right singular vectors diagonalize S^dagger S without squaring small survival amplitudes
  Matrix eigenvectors;
  std::vector<double> singular;
  try {
    singular = smatrix.SingularValues(&eigenvectors);
  } catch (const std::runtime_error &) {
    throw std::runtime_error(
        "MEikonalMatrix::CutSpectrum: survival eigensystem failed");
  }

  const auto incoming_states = IncomingCutStates();
  if (incoming_states.empty()) {
    throw std::runtime_error("MEikonalMatrix::CutSpectrum: no incoming states");
  }

  CutOpacitySpectrum spectrum;
  const std::size_t dimension = smatrix.size_row();
  spectrum.eigenvalues.resize(dimension, 0.0);
  spectrum.incoming_weights.resize(dimension, 0.0);

  const double contractivity_tolerance = unitarity_tolerance;

  double sigma_max = 0.0;
  for (const double sigma : singular) {
    if (!std::isfinite(sigma) || sigma < 0.0) {
      throw std::runtime_error(
          "MEikonalMatrix::CutSpectrum: invalid survival singular value");
    }
    sigma_max = std::max(sigma_max, sigma);
  }

  // A probabilistic cut interpretation requires contractivity. Tiny
  // roundoff overshoots are clipped, while a genuine violation is rejected
  if (sigma_max > 1.0 + contractivity_tolerance) {
    std::ostringstream message;
    message << std::scientific << std::setprecision(8)
            << "MEikonalMatrix::CutSpectrum: Pomeron S matrix is not "
               "contractive (b="
            << bt << ", sigma_max=" << sigma_max
            << ", tolerance=" << contractivity_tolerance << ')';
    throw std::runtime_error(message.str());
  }

  for (std::size_t mode = 0; mode < dimension; ++mode) {
    // Retain ascending survival order and allow the completely absorptive limit
    const std::size_t index = dimension - 1 - mode;
    const double sigma = std::min(1.0, singular[index]);
    spectrum.eigenvalues[mode] = sigma > 0.0 ? -2.0 * std::log(sigma) : std::numeric_limits<double>::infinity();

    const std::vector<Complex> eigenvector = eigenvectors.Column(index);
    double weight = 0.0;
    for (const auto &incoming : incoming_states) {
      if (incoming.size() != dimension) {
        throw std::runtime_error(
            "MEikonalMatrix::CutSpectrum: incoming-state dimension mismatch");
      }
      const Complex overlap = gra::InnerProduct(
          eigenvector.cbegin(), eigenvector.cend(), incoming.cbegin());
      weight += std::norm(overlap);
    }
    spectrum.incoming_weights[mode] =
        weight / static_cast<double>(incoming_states.size());
  }

  const double weight_sum = gra::Sum(spectrum.incoming_weights);
  if (!(weight_sum > 0.0) || !std::isfinite(weight_sum)) {
    throw std::runtime_error(
        "MEikonalMatrix::CutSpectrum: invalid incoming spectral normalization");
  }
  gra::Scale(spectrum.incoming_weights, 1.0 / weight_sum);

  return spectrum;
}

// Compute unpolarized incoming physical proton states in the cut-opacity basis
std::vector<std::vector<double>> MEikonalMatrix::IncomingCutStates() const {
  const std::vector<double> pair =
      PairStateVector(tables, soft_model->GoodWalker(), 0, 0);
  if (!soft_model->Eikonal().HelicityEnabled()) {
    return {pair};
  }
  std::vector<std::vector<double>> states(
      4, std::vector<double>(4 * tables.D, 0.0));
  for (const auto &spin : indices(states)) {
    std::copy(pair.cbegin(), pair.cend(),
              states[spin].begin() + spin * tables.D);
  }
  return states;
}

} // namespace gra
