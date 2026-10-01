// Eikonal density and screening class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <atomic>
#include <bit>
#include <cmath>
#include <cstdint>
#include <exception>
#include <filesystem>
#include <iostream>
#include <limits>
#include <memory>
#include <mutex>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Eikonal/MEikonalHelicity.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "json.hpp"
using gra::aux::indices;

using gra::math::zi;

namespace gra {

// Coupled-channel matrix helpers moved to MEikonalMatrix.cc

using Complex = std::complex<double>;
using Matrix = MEikonalMatrix::Matrix;

// Construct an empty eikonal helper. Heavy initialization is explicit via
// S3Constructor
MEikonal::MEikonal() {}

// Construct an eikonal helper bound to one immutable model tune
MEikonal::MEikonal(MModelTunePtr tune)
    : model_tune(std::move(tune)),
      soft_model(model_tune != nullptr ? model_tune->Soft() : nullptr) {
  if (model_tune == nullptr || soft_model == nullptr) {
    throw std::invalid_argument("MEikonal::MEikonal: null model tune");
  }
}

// Eikonal helper destructor
MEikonal::~MEikonal() = default;

// Clear every runtime value owned by one eikonal initialization
void MEikonal::ResetRuntime() {
  S3INIT = false;
  screening_amplitude_initialized = false;
  matrix_runtime.reset();
  elastic_cni_runtime.reset();
  U = MMatrix<double>();
  s = 0.0;
  INITIALSTATE.clear();
  sigma_tot = 0.0;
  sigma_el = 0.0;
  sigma_inel = 0.0;
  P_cut.clear();
  P_array.clear();
  cut_b.clear();
  MAX_P_array = 0.0;
  MAX_P_cut = 0.0;
  min_opacity_eigenvalue = 0.0;
  sigma_diff = MMatrix<double>(2, 2, 0.0);
  sigma_diff_exclusive = MMatrix<double>();
  loop_const = LoopConst{};
}

// Compute total, elastic, inelastic cross sections
void MEikonal::GetTotXS(double &tot, double &el, double &in) const {
  tot = sigma_tot;
  el = sigma_el;
  in = sigma_inel;
}

// Compute number of radial screening loop nodes for the current radial rule
unsigned int MEikonal::GetLoopKTNodeCount() const {
  return math::PolarNodeCount(Numerics.LOOP);
}

// Compute the radial quadrature coefficient for the absorptive loop integral
double MEikonal::GetLoopQuadratureCoefficient() const {
  return math::PolarCoefficient(Numerics.LOOP);
}

// Compare quadrature controls exactly as configuration data, including floating point bit patterns
bool MEikonal::LoopNumericsMatch() const {
  const auto &cached = loop_const.quadrature;
  const auto &current = Numerics.LOOP;
  return loop_const.initialized && cached.radial_integrator == current.radial_integrator &&
         cached.azimuth_integrator == current.azimuth_integrator && cached.radial_map == current.radial_map &&
         cached.radial_intervals == current.radial_intervals && cached.azimuth_nodes == current.azimuth_nodes &&
         std::bit_cast<std::uint64_t>(cached.r_min) == std::bit_cast<std::uint64_t>(current.r_min) &&
         std::bit_cast<std::uint64_t>(cached.r_max) == std::bit_cast<std::uint64_t>(current.r_max);
}

// Initialize the polar quadrature for the transverse screening convolution
// The node measure contains kT dkT dphi for both linear and logarithmic maps
void MEikonal::InitLoopWeightMatrix() {
  if (s <= 0.0) {
    throw std::invalid_argument("MEikonal::InitLoopWeightMatrix: Mandelstam s "
                                "must be initialized first");
  }

  const auto &loop = Numerics.LOOP;

  if (loop.radial_intervals == 0 || loop.azimuth_nodes == 0) {
    throw std::invalid_argument(
        "MEikonal::InitLoopWeightMatrix: LOOP.NumberLoopKT and "
        "LOOP.NumberLoopPHI must be positive");
  }
  if (loop.r_min < 0.0) {
    throw std::invalid_argument(
        "MEikonal::InitLoopWeightMatrix: LOOP.MinLoopKT must be >= 0");
  }
  if (loop.r_max <= loop.r_min) {
    throw std::invalid_argument("MEikonal::InitLoopWeightMatrix: "
                                "LOOP.MinLoopKT < LOOP.MaxLoopKT required");
  }
  if (loop.radial_map == math::RadialMap::Log && loop.r_min <= 0.0) {
    throw std::invalid_argument("MEikonal::InitLoopWeightMatrix: LOOP.log_kT = "
                                "true requires LOOP.MinLoopKT > 0");
  }

  const math::PolarRule1D kt_rule = math::PolarRadialRule(loop);
  const math::PolarRule1D phi_rule = math::PolarAzimuthRule(loop);

  if (INITIALSTATE.size() != 2) {
    throw std::invalid_argument("MEikonal::InitLoopWeightMatrix: exactly two beams are required");
  }
  const double beta = kinematics::beta12(s, INITIALSTATE[0].mass, INITIALSTATE[1].mass);
  if (!std::isfinite(beta) || !(beta > 0.0)) {
    throw std::invalid_argument("MEikonal::InitLoopWeightMatrix: initial state is below threshold");
  }
  // The absorptive correction carries i/[8 pi^2 sqrt(lambda)]
  const std::complex<double> norm = zi / (8.0 * gra::math::PIPI * s * beta);
  const std::size_t n_kt = kt_rule.node.size();
  const std::size_t n_phi = phi_rule.node.size();

  LoopConst loop_const;
  loop_const.quadrature = loop;
  loop_const.initialized = true;
  loop_const.log_kT = loop.radial_map == math::RadialMap::Log;
  loop_const.StepKT = kt_rule.step;
  loop_const.MinLogKT = kt_rule.map_min;
  loop_const.StepPhi = phi_rule.step;
  loop_const.quadrature_scale = 1.0;
  loop_const.s = s;
  loop_const.norm = norm;
  loop_const.kt = kt_rule.node;
  loop_const.kt2.resize(n_kt);
  std::transform(kt_rule.node.cbegin(), kt_rule.node.cend(),
                 loop_const.kt2.begin(),
                 [](const double kt) { return math::pow2(kt); });
  loop_const.jac = kt_rule.jac;
  loop_const.W = gra::OuterProduct(kt_rule.weight, phi_rule.weight);
  loop_const.measure_weight =
      gra::OuterProduct(kt_rule.measure, phi_rule.measure);

  std::vector<double> cos_phi(n_phi);
  std::vector<double> sin_phi(n_phi);
  std::transform(phi_rule.node.cbegin(), phi_rule.node.cend(), cos_phi.begin(),
                 [](const double phi) { return std::cos(phi); });
  std::transform(phi_rule.node.cbegin(), phi_rule.node.cend(), sin_phi.begin(),
                 [](const double phi) { return std::sin(phi); });
  loop_const.kt_x = gra::OuterProduct(kt_rule.node, cos_phi);
  loop_const.kt_y = gra::OuterProduct(kt_rule.node, sin_phi);
  loop_const.node_weight = loop_const.measure_weight * norm;

  this->loop_const = std::move(loop_const);
  if (screening_amplitude_initialized) {
    InitLoopPhysicalScreeningCache();
  }
}

// Compute event-by-event screening loop quadrature weights matching current
// numerics
const MMatrix<double> &MEikonal::GetLoopWeightMatrix() const {
  const unsigned int n_kt = GetLoopKTNodeCount();
  const unsigned int n_phi = Numerics.LOOP.azimuth_nodes;

  if (!LoopNumericsMatch() || loop_const.W.size_row() != n_kt ||
      loop_const.W.size_col() != n_phi) {
    throw std::invalid_argument(
        "MEikonal::GetLoopWeightMatrix: loop weights are not initialized for "
        "current loop numerics");
  }
  return loop_const.W;
}

// Compute event-by-event screening loop quadrature constants matching current
// numerics
const MEikonal::LoopConst &MEikonal::GetLoopConst(double current_s) const {
  const unsigned int n_kt = GetLoopKTNodeCount();
  const unsigned int n_phi = Numerics.LOOP.azimuth_nodes;

  if (!LoopNumericsMatch() || loop_const.kt2.size() != n_kt ||
      loop_const.kt_x.size_row() != n_kt ||
      loop_const.kt_x.size_col() != n_phi ||
      loop_const.kt_y.size_row() != n_kt ||
      loop_const.kt_y.size_col() != n_phi ||
      loop_const.measure_weight.size_row() != n_kt ||
      loop_const.measure_weight.size_col() != n_phi ||
      loop_const.node_weight.size_row() != n_kt ||
      loop_const.node_weight.size_col() != n_phi ||
      loop_const.physical_screening_amplitude.size() != n_kt ||
      loop_const.pair_screening_spin.size() != n_kt ||
      loop_const.physical_screening_weight.size_row() != n_kt ||
      loop_const.physical_screening_weight.size_col() != n_phi ||
      loop_const.physical_screening_helicity_weight.size() != n_kt * n_phi) {
    throw std::invalid_argument("MEikonal::GetLoopConst: loop constants are "
                                "not initialized for current loop numerics");
  }
  const double scale = (std::abs(current_s) > 1.0) ? std::abs(current_s) : 1.0;
  if (!std::isfinite(current_s) || std::abs(loop_const.s - current_s) > 1e-10 * scale) {
    throw std::invalid_argument("MEikonal::GetLoopConst: stored Mandelstam s "
                                "does not match current event");
  }

  return loop_const;
}

// Build the coherent Good Walker source for one diffractive class
// State zero is the proton and the resolved direction spans its excitations
static std::vector<std::pair<GWExclusiveState, double>>
GoodWalkerSource(const GoodWalkerSpace &good_walker,
                 const GWFinalClass final_class) {
  const std::size_t channel_count = good_walker.ChannelCount();
  if (channel_count == 1) {
    if (final_class == GWFinalClass::Elastic) {
      return {{{0, 0}, 1.0}};
    }
    return {};
  }
  const auto &resolved = good_walker.ResolvedCoefficients();
  if (resolved.size() + 1 != channel_count) {
    throw std::logic_error("GoodWalkerSource: resolved direction coefficient "
                           "count is inconsistent");
  }
  std::vector<std::pair<GWExclusiveState, double>> source;
  switch (final_class) {
  case GWFinalClass::Elastic:
    source.push_back({{0, 0}, 1.0});
    break;
  case GWFinalClass::SingleDissociation1:
    for (std::size_t c = 1; c < channel_count; ++c) {
      source.push_back({{c, 0}, resolved[c - 1]});
    }
    break;
  case GWFinalClass::SingleDissociation2:
    for (std::size_t c = 1; c < channel_count; ++c) {
      source.push_back({{0, c}, resolved[c - 1]});
    }
    break;
  case GWFinalClass::DoubleDissociation:
    for (std::size_t c = 1; c < channel_count; ++c) {
      for (std::size_t d = 1; d < channel_count; ++d) {
        source.push_back({{c, d}, resolved[c - 1] * resolved[d - 1]});
      }
    }
    break;
  }
  return source;
}

// Test whether a proton pair operator is proportional to the spin identity
static bool IsScalarProtonSpinOperator(const ProtonHelicityMatrix &matrix) {
  const std::complex<double> scalar = matrix[0];
  const double tolerance = 1.0e-12 * (1.0 + std::abs(scalar));
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const std::complex<double> expected =
          row == col ? scalar : std::complex<double>(0.0, 0.0);
      if (std::abs(matrix[gra::spin::PairHelicityMatrixIndex(row, col)] -
                   expected) > tolerance) {
        return false;
      }
    }
  }
  return true;
}

// Test whether the full Good Walker operator is scalar in proton spin space
static bool IsScalarPairSpinOperator(const MEikonalMatrix::PairSpinBank &bank) {
  const std::size_t rows = bank[0].size_row();
  const std::size_t cols = bank[0].size_col();
  if (rows == 0 || rows != cols) {
    return false;
  }
  const MMatrix<std::complex<double>> zero(rows, cols, 0.0);
  for (std::size_t row = 0; row < 4; ++row) {
    for (std::size_t col = 0; col < 4; ++col) {
      const auto &expected = row == col ? bank[0] : zero;
      if (!bank[gra::spin::PairHelicityMatrixIndex(row, col)].IsApprox(expected,
                                                                       0.0)) {
        return false;
      }
    }
  }
  return true;
}

// Validate one immutable radial pair-spin operator before event-loop reuse
static void
ValidatePairSpinOperatorCache(const MEikonalMatrix::PairSpinBank &bank,
                              const std::size_t pair_dimension) {
  for (const auto &matrix : bank) {
    if (matrix.size_row() != pair_dimension ||
        matrix.size_col() != pair_dimension || !matrix.IsFinite()) {
      throw std::logic_error(
          "MEikonal: invalid cached pair-space screening operator");
    }
  }
}

// Compute w_n = w exp(i n phi) for helicity harmonics n=-2,...,+2
static std::array<std::complex<double>, 5>
ScreeningHarmonicWeight(const double kt_x, const double kt_y,
                        const std::complex<double> node_weight) {
  const double azimuth = std::atan2(kt_y, kt_x);
  const std::complex<double> z = std::exp(std::complex<double>(0.0, azimuth));
  const std::complex<double> z2 = z * z;
  return {node_weight * std::conj(z2), node_weight * std::conj(z), node_weight,
          node_weight * z, node_weight * z2};
}

// Cache radial screening operators and their helicity azimuthal phases
// Good Walker transition amplitudes remain coherent until the final state sum
void MEikonal::InitLoopPhysicalScreeningCache() {
  const std::size_t n_kt = loop_const.kt2.size();
  const std::size_t n_phi = loop_const.node_weight.size_col();
  if (!loop_const.initialized || GetChannelCount() == 0 || n_kt == 0) {
    throw std::logic_error("MEikonal::InitLoopPhysicalScreeningCache: eikonal "
                           "loop is not initialized");
  }
  const std::size_t pair_dimension = soft_model->GoodWalker().PairDimension();

  loop_const.physical_screening_amplitude.assign(n_kt, 0.0);
  loop_const.pair_screening_spin.resize(n_kt);
  loop_const.physical_screening_weight =
      MMatrix<std::complex<double>>(n_kt, n_phi, 0.0);
  loop_const.physical_screening_helicity_weight.assign(n_kt * n_phi, {});
  loop_const.physical_screening_is_scalar = true;
  loop_const.pair_screening_cache_ready = false;
  loop_const.pair_screening_spin_scalar = true;
  loop_const.pair_screening_pair_diagonal = true;
  loop_const.pair_screening_diagonal.assign(n_kt, {});
  loop_const.pair_screening_harmonic_weight.assign(n_kt * n_phi, {});
  for (std::size_t kt_index = 0; kt_index < n_kt; ++kt_index) {
    loop_const.pair_screening_spin[kt_index] =
        PairScreeningSpinBank(loop_const.kt2[kt_index]);
    const auto &pair_bank = loop_const.pair_screening_spin[kt_index];
    ValidatePairSpinOperatorCache(pair_bank, pair_dimension);
    const auto classify = [&](const auto &bank, bool &all_scalar,
                              bool &all_diagonal, auto &diagonal) {
      const bool scalar = IsScalarPairSpinOperator(bank);
      const bool pair_diagonal = scalar && bank[0].IsDiagonal();
      all_scalar = all_scalar && scalar;
      all_diagonal = all_diagonal && pair_diagonal;
      if (pair_diagonal) {
        diagonal = bank[0].GetDiag();
      }
    };
    classify(pair_bank, loop_const.pair_screening_spin_scalar,
             loop_const.pair_screening_pair_diagonal,
             loop_const.pair_screening_diagonal[kt_index]);
    loop_const.physical_screening_amplitude[kt_index] =
        S3PhysicalScreeningAmp(loop_const.kt2[kt_index]);
    const ProtonHelicityMatrix reference =
        GetMatrixRuntime().ScreeningHelicityMatrix(loop_const.kt2[kt_index], 0,
                                                   0);
    for (std::size_t phi_index = 0; phi_index < n_phi; ++phi_index) {
      loop_const.pair_screening_harmonic_weight[kt_index * n_phi + phi_index] =
          ScreeningHarmonicWeight(loop_const.kt_x[kt_index][phi_index],
                                  loop_const.kt_y[kt_index][phi_index],
                                  loop_const.node_weight[kt_index][phi_index]);
      loop_const.physical_screening_weight[kt_index][phi_index] =
          loop_const.node_weight[kt_index][phi_index] *
          loop_const.physical_screening_amplitude[kt_index];
      auto matrix = RotateProtonHelicityMatrix(
          reference, std::atan2(loop_const.kt_y[kt_index][phi_index],
                                loop_const.kt_x[kt_index][phi_index]));
      gra::Scale(matrix, loop_const.node_weight[kt_index][phi_index]);
      const std::complex<double> scalar_weight =
          loop_const.physical_screening_weight[kt_index][phi_index];
      const double tolerance = 1.0e-12 * (1.0 + std::abs(scalar_weight));
      loop_const.physical_screening_is_scalar =
          loop_const.physical_screening_is_scalar &&
          IsScalarProtonSpinOperator(matrix) &&
          std::abs(matrix[0] - scalar_weight) <= tolerance;
      loop_const
          .physical_screening_helicity_weight[kt_index * n_phi + phi_index] =
          matrix;
    }
  }

  const std::array<GWFinalClass, 4> final_classes = {
      GWFinalClass::Elastic, GWFinalClass::SingleDissociation1,
      GWFinalClass::SingleDissociation2, GWFinalClass::DoubleDissociation};
  for (const auto &class_index : indices(final_classes)) {
    auto &channels = loop_const.good_walker_channels[class_index];
    channels.clear();
    const auto source =
        GoodWalkerSource(soft_model->GoodWalker(), final_classes[class_index]);
    for (const auto &[final_state, born_coefficient] : source) {
      LoopGoodWalkerChannel channel;
      channel.final_state = final_state;
      channel.born_coefficient = born_coefficient;
      channel.scalar_screening_weight.assign(n_kt * n_phi, 0.0);
      channel.screening_weight.assign(n_kt * n_phi, {});
      for (std::size_t kt_index = 0; kt_index < n_kt; ++kt_index) {
        ProtonHelicityMatrix effective{};
        for (const auto &[initial_state, coefficient] : source) {
          const auto transition = GetMatrixRuntime().ScreeningHelicityMatrix(
              loop_const.kt2[kt_index], final_state.f1, final_state.f2,
              initial_state.f1, initial_state.f2);
          gra::AddScaled(effective, transition, coefficient);
        }
        for (std::size_t phi_index = 0; phi_index < n_phi; ++phi_index) {
          auto matrix = RotateProtonHelicityMatrix(
              effective, std::atan2(loop_const.kt_y[kt_index][phi_index],
                                    loop_const.kt_x[kt_index][phi_index]));
          gra::Scale(matrix, loop_const.node_weight[kt_index][phi_index]);
          channel.spin_scalar =
              channel.spin_scalar && IsScalarProtonSpinOperator(matrix);
          channel.scalar_screening_weight[kt_index * n_phi + phi_index] =
              matrix[0];
          channel.screening_weight[kt_index * n_phi + phi_index] = matrix;
        }
      }
      if (channel.spin_scalar) {
        std::vector<ProtonHelicityMatrix>().swap(channel.screening_weight);
      } else {
        std::vector<std::complex<double>>().swap(
            channel.scalar_screening_weight);
      }
      channels.push_back(std::move(channel));
    }
  }
  if (!loop_const.pair_screening_pair_diagonal) {
    std::vector<std::vector<std::complex<double>>>().swap(
        loop_const.pair_screening_diagonal);
  }
  loop_const.pair_screening_cache_ready = true;
}

// Build the coupled Good Walker eikonal at fixed Mandelstam s
// Born exchanges are transformed to b space, unitarized as matrices and then
// transformed to qT space when momentum space amplitudes are requested
void MEikonal::S3Constructor(double s_in,
                             const std::vector<gra::MParticle> &initialstate_in,
                             bool onlydensity, int NumberBT, int NumberKT2,
                             double max_kt2) {
  std::cout << "MEikonal::S3Constructor:" << std::endl;

  if (soft_model == nullptr) {
    throw std::invalid_argument(
        "MEikonal::S3Constructor: immutable SOFT model is not set");
  }
  if (!std::isfinite(max_kt2) || max_kt2 < 0.0) {
    throw std::invalid_argument(
        "MEikonal::S3Constructor: max_kt2 must be finite and nonnegative");
  }

  const int extra_loop_kt = Numerics.extra_NumberLoopKT;
  const int extra_loop_phi = Numerics.extra_NumberLoopPHI;
  ResetRuntime();

  try {
    MEikonalNumerics numerics;
    numerics.SetLoopDiscretization(extra_loop_kt, extra_loop_phi);

    if (model_tune == nullptr) {
      throw std::logic_error(
          "MEikonal::S3Constructor: immutable model tune is not set");
    }
    numerics.ReadParameters(model_tune->NumericsFile(),
                            model_tune->Numerics().dump());

    // Apply explicit lookup-table overrides
    if (NumberBT > 0) {
      numerics.NumberBT = NumberBT;
    }
    if (NumberKT2 > 0) {
      numerics.NumberKT2 = NumberKT2;
    }
    numerics.MaxKT2 = std::max(numerics.MaxKT2, max_kt2);
    if (numerics.NumberBT < 2 || numerics.NumberBT % 2 != 0) {
      throw std::invalid_argument(
          "MEikonal::S3Constructor: NumberBT must be even and at least two");
    }
    Numerics = std::move(numerics);

    // Store Mandelstam s and the ordered initial state
    s = s_in;
    INITIALSTATE = initialstate_in;
    InitLoopWeightMatrix();

    // Construct the rotation from diffractive eigenstates to physical states
    if (soft_model->ExchangeCount() == 0) {
      throw std::invalid_argument(
          "MEikonal::S3Constructor: No soft exchanges configured");
    }
    const std::size_t N = soft_model->GoodWalker().ChannelCount();
    U = soft_model->GoodWalker().MixingMatrix();

    std::cout << "Using " << N << "-channel Eikonal model" << std::endl
              << std::endl;
    std::cout << "Eigenbasis orthogonal mixing matrix U:" << std::endl;
    std::cout << U << std::endl;

    const bool build_amplitude = !onlydensity;
    matrix_runtime = MEikonalMatrix::Build(soft_model, s, INITIALSTATE,
                                           Numerics, build_amplitude);
    S3CalcEikonalXS(build_amplitude);
    if (build_amplitude) {
      S3InitCutPomerons();
      screening_amplitude_initialized = true;
      InitLoopPhysicalScreeningCache();
    }
    S3INIT = true;
  } catch (...) {
    ResetRuntime();
    Numerics = MEikonalNumerics{};
    Numerics.SetLoopDiscretization(extra_loop_kt, extra_loop_phi);
    throw;
  }
}

// Compute one Regge pole Born amplitude between eigenstates i and k
// A_ik(s,t) = epsilon eta(alpha) beta_ii(t) beta_kk(t) (s/s0)^alpha(t)
// eta fixes the signature phase and epsilon fixes the fitted coupling sign
std::complex<double> MEikonal::SingleAmpElastic(const SoftModelPtr &model,
                                                const double s, const double t,
                                                const SoftExchangeId exchange,
                                                const std::size_t _i,
                                                const std::size_t _k) {
  if (model == nullptr) {
    throw std::invalid_argument("MEikonal::SingleAmpElastic: null SOFT model");
  }
  const auto &parameter = model->Exchange(exchange);
  const auto residue = model->ResidueMatrix(exchange, t);
  if (_i >= residue.size_row() || _k >= residue.size_row()) {
    throw std::out_of_range(
        "MEikonal::SingleAmpElastic: Good-Walker channel is out of range");
  }

  const double alpha = model->Alpha(exchange, t);
  const std::complex<double> eta =
      gra::regge::EtaFactor(alpha, parameter.Alpha0(), parameter.Signature(),
                            parameter.Eta());
  const double residue_sign = static_cast<double>(parameter.ResidueSign());
  constexpr double s0 = 1.0;

  return residue_sign * eta * residue(_i, _i) * residue(_k, _k) *
         std::pow(s / s0, alpha);
}

// Compute the C-odd crossing sign, +1 for pp and -1 for proton antiproton
double MEikonal::OddProfileBeamSign() const {
  if (INITIALSTATE.size() != 2) {
    throw std::invalid_argument(
        "MEikonal::OddProfileBeamSign: exactly two beams are required");
  }
  if (INITIALSTATE[0].pdg == INITIALSTATE[1].pdg) {
    return 1.0;
  }
  if (INITIALSTATE[0].pdg == -INITIALSTATE[1].pdg) {
    return -1.0;
  }
  throw std::invalid_argument(
      "MEikonal::OddProfileBeamSign: unsupported incoming beam pair");
}

// Interpolate the C-even and C-odd eigenchannel profiles Xi_ik(s,b)
// Keeping both crossing sectors separate lets the beam state set the odd sign
// [REFERENCE: Desgrolard, Giffon, Martynov, Predazzi, arxiv.org/abs/hep-ph/9907451v2]
// [REFERENCE: Desgrolard, Giffon, Martynov, Predazzi, arxiv.org/abs/hep-ph/0001149]
// [REFERENCE: Ewerz, Maniatis, Nachtmann, arxiv.org/abs/1309.3478]
std::pair<std::complex<double>, std::complex<double>>
MEikonal::S3DensityProfiles(double bt, std::size_t _i, std::size_t _k) const {
  if (bt < 0.0) {
    throw std::invalid_argument(
        "MEikonal::S3DensityProfiles: negative impact parameter");
  }
  return GetMatrixRuntime().DensityProfiles(bt, _i, _k);
}

// Form the beam eikonal for one Good Walker eigenstate pair
// Xi_beam(s,b) = Xi_even(s,b) + sigma_beam Xi_odd(s,b)
std::complex<double> MEikonal::S3Density(double bt, std::size_t _i,
                                         std::size_t _k) const {

  // pp and pbar p differ by the C-odd exchange sign
  const double sign = OddProfileBeamSign();
  const auto profiles = S3DensityProfiles(bt, _i, _k);
  return profiles.first + sign * profiles.second;
}

// Transform one unitarized eigenchannel profile from b to qT space
// A_ik(s,qT) = 4 pi sqrt(lambda) integral b db J0(b qT) Gamma_ik(s,b)
// Gamma=i(I-S) and sqrt(lambda)=s beta12 for the two incoming particles
std::complex<double> MEikonal::S3Screening(double kt2, std::size_t _i,
                                           std::size_t _k) const {
  if (kt2 < 0.0) {
    throw std::invalid_argument("MEikonal::S3Screening: negative kt2");
  }
  return GetMatrixRuntime().ScreeningAmplitude(kt2, _i, _k);
}

// Project the screening operator onto elastic physical proton states
// This is the spin averaged <pp|A_screen|pp> amplitude
std::complex<double> MEikonal::S3PhysicalScreeningAmp(const double kt2) const {
  if (kt2 < 0.0) {
    throw std::invalid_argument(
        "MEikonal::S3PhysicalScreeningAmp: negative kt2");
  }
  return S3ExclusiveScreeningAmp(kt2, 0, 0);
}

// Project the full eikonal onto an exclusive Good Walker final state
// State zero is the proton and higher indices are orthogonal excitations
// The returned scalar is one quarter of the pair helicity trace
std::complex<double> MEikonal::S3ExclusiveAmp(double kt2, std::size_t f1,
                                              std::size_t f2) const {
  if (kt2 < 0.0) {
    throw std::invalid_argument("MEikonal::S3ExclusiveAmp: negative kt2");
  }
  const auto amplitude = GetMatrixRuntime().HelicityMatrix(kt2, f1, f2);
  return 0.25 * (amplitude[0] + amplitude[5] + amplitude[10] + amplitude[15]);
}

// Project the selected screening exchanges between Good Walker pair states
// The returned scalar is one quarter of the pair helicity trace
std::complex<double>
MEikonal::S3ExclusiveScreeningAmp(const double kt2, const std::size_t f1,
                                  const std::size_t f2, const std::size_t i1,
                                  const std::size_t i2) const {
  if (kt2 < 0.0) {
    throw std::invalid_argument(
        "MEikonal::S3ExclusiveScreeningAmp: negative kt2");
  }
  const std::size_t channel_count = U.size_row();
  if (f1 >= channel_count || f2 >= channel_count || i1 >= channel_count ||
      i2 >= channel_count) {
    throw std::invalid_argument(
        "MEikonal::S3ExclusiveScreeningAmp: state index out of range");
  }
  const auto amplitude =
      GetMatrixRuntime().ScreeningHelicityMatrix(kt2, f1, f2, i1, i2);
  return 0.25 * (amplitude[0] + amplitude[5] + amplitude[10] + amplitude[15]);
}

// Enumerate the orthogonal final states in one Good Walker diffraction class
std::vector<GWExclusiveState>
MEikonal::ExclusiveFinalStates(GWFinalClass final_class) const {
  const std::size_t N = U.size_row();
  if (N == 0) {
    throw std::invalid_argument(
        "MEikonal::ExclusiveFinalStates: Eikonal basis not initialized");
  }

  std::vector<GWExclusiveState> states;
  if (final_class == GWFinalClass::Elastic) {
    states.push_back({0, 0});
  } else if (final_class == GWFinalClass::SingleDissociation1) {
    for (std::size_t f1 = 1; f1 < N; ++f1) {
      states.push_back({f1, 0});
    }
  } else if (final_class == GWFinalClass::SingleDissociation2) {
    for (std::size_t f2 = 1; f2 < N; ++f2) {
      states.push_back({0, f2});
    }
  } else if (final_class == GWFinalClass::DoubleDissociation) {
    for (std::size_t f1 = 1; f1 < N; ++f1) {
      for (std::size_t f2 = 1; f2 < N; ++f2) {
        states.push_back({f1, f2});
      }
    }
  } else {
    throw std::invalid_argument(
        "MEikonal::ExclusiveFinalStates: Unknown final-state class");
  }
  return states;
}

// Sum incoherently over orthogonal Good Walker final states
// <|A|^2> = 1/4 sum_(f1,f2,lambda_i,lambda_f) |A_f1f2|^2
double MEikonal::S3DiffractiveAmpSquared(double kt2,
                                         GWFinalClass final_class) const {
  double amp_squared = 0.0;
  for (const auto &state : ExclusiveFinalStates(final_class)) {
    amp_squared += 0.25 * gra::SquaredNorm(GetMatrixRuntime().HelicityMatrix(
                              kt2, state.f1, state.f2));
  }
  return amp_squared;
}

// Integrate total and exclusive Good Walker cross sections in b space
// sigma_tot = 4 pi Im integral b db Tr_spin[Gamma_00]/4
// sigma_f1f2 = 2 pi integral b db sum_spin |Gamma_f1f2|^2/4
// sigma_inel = sigma_tot-sigma_el includes diffractive and nondiffractive cuts
void MEikonal::S3CalcEikonalXS(const bool amplitude_diagnostics) {
  std::cout << "MEikonal::S3CalcEikonalXS:" << std::endl << std::endl;

  const std::size_t nch = U.size_row();
  if (nch == 0) {
    throw std::logic_error("MEikonal::S3CalcEikonalXS: eikonal basis is empty");
  }

  sigma_diff_exclusive = MMatrix<double>(nch, nch, 0.0);
  sigma_diff = MMatrix<double>(2, 2, 0.0);
  // Gamma=i(I-S) obeys the amplitude normalization used by the optical theorem
  const auto &b_node = matrix_runtime->ImpactParameterNodes();
  if (b_node.size() < 2) {
    throw std::logic_error(
        "MEikonal::S3CalcEikonalXS: incomplete impact-parameter table");
  }
  const double bstep =
      (Numerics.logBT ? std::log(Numerics.MaxBT) - std::log(Numerics.MinBT)
                      : Numerics.MaxBT - Numerics.MinBT) /
      Numerics.NumberBT;
  std::vector<std::vector<std::vector<std::complex<double>>>> integrand(
      nch, std::vector<std::vector<std::complex<double>>>(
               nch, std::vector<std::complex<double>>(b_node.size(), 0.0)));
  std::vector<std::complex<double>> optical_integrand(b_node.size(), 0.0);
  for (const auto &point : indices(b_node)) {
    const double b = b_node[point];
    const double measure = Numerics.logBT ? b * b : b;
    for (std::size_t f1 = 0; f1 < nch; ++f1) {
      for (std::size_t f2 = 0; f2 < nch; ++f2) {
        integrand[f1][f2][point] =
            measure * 0.25 *
            gra::SquaredNorm(matrix_runtime->ImpactHelicityMatrix(b, f1, f2));
      }
    }
    const ProtonHelicityMatrix amplitude =
        matrix_runtime->ImpactHelicityMatrix(b, 0, 0);
    optical_integrand[point] =
        measure * 0.25 *
        (amplitude[0] + amplitude[5] + amplitude[10] + amplitude[15]);
  }
  for (std::size_t f1 = 0; f1 < nch; ++f1) {
    for (std::size_t f2 = 0; f2 < nch; ++f2) {
      sigma_diff_exclusive[f1][f2] =
          2.0 * gra::math::PI *
          std::real(gra::math::CS13Integral(integrand[f1][f2], bstep)) *
          PDG::GeV2barn;
    }
  }
  sigma_tot = 4.0 * gra::math::PI *
              std::imag(gra::math::CS13Integral(optical_integrand, bstep)) *
              PDG::GeV2barn;

  sigma_diff[0][0] = sigma_diff_exclusive[0][0];
  for (std::size_t i = 1; i < nch; ++i) {
    sigma_diff[1][0] += sigma_diff_exclusive[i][0];
    sigma_diff[0][1] += sigma_diff_exclusive[0][i];
  }
  for (std::size_t i = 1; i < nch; ++i) {
    for (std::size_t j = 1; j < nch; ++j) {
      sigma_diff[1][1] += sigma_diff_exclusive[i][j];
    }
  }

  sigma_el = sigma_diff[0][0];
  sigma_inel = sigma_tot - sigma_el;

  if (!std::isfinite(sigma_tot) || !std::isfinite(sigma_el) ||
      !std::isfinite(sigma_inel)) {
    throw std::runtime_error(
        "MEikonal::S3CalcEikonalXS: non-finite cross section");
  }
  if (sigma_inel < -1.0e-10 * std::max(1.0, std::abs(sigma_tot))) {
    throw std::runtime_error(
        "MEikonal::S3CalcEikonalXS: negative inelastic "
        "cross section; check the S-matrix unitarity bound");
  }
  sigma_inel = std::max(0.0, sigma_inel);

  printf(" Total xs:      %0.3f mb \n", sigma_tot * 1E3);
  printf(" Inelastic xs:  %0.3f mb \n\n", sigma_inel * 1E3);
  printf(" Elastic xs:    %0.3f mb \n\n", sigma_el * 1E3);

  if (!amplitude_diagnostics) {
    return;
  }

  if (nch > 1) {
    printf(" Exclusive Good-Walker final-state matrix [mb]:\n");
    for (std::size_t f1 = 0; f1 < nch; ++f1) {
      printf("  f=%zu :", f1);
      for (std::size_t f2 = 0; f2 < nch; ++f2) {
        printf(" %0.6f", sigma_diff_exclusive[f1][f2] * 1E3);
      }
      printf("\n");
    }

    double min_im = std::numeric_limits<double>::max();
    double max_im = -std::numeric_limits<double>::max();
    double mean_im = 0.0;
    double mean2_im = 0.0;
    std::size_t count = 0;
    printf("\n Optical-point diagonal pair-state Im A_ik(0):\n");
    for (std::size_t i = 0; i < nch; ++i) {
      printf("  i=%zu :", i);
      for (std::size_t k = 0; k < nch; ++k) {
        const double value =
            std::imag(GetMatrixRuntime().ScreeningAmplitude(0.0, i, k));
        printf(" %0.3e", value);
        min_im = std::min(min_im, value);
        max_im = std::max(max_im, value);
        mean_im += value;
        mean2_im += value * value;
        ++count;
      }
      printf("\n");
    }
    mean_im /= static_cast<double>(count);
    mean2_im /= static_cast<double>(count);
    const double rms = std::sqrt(std::max(0.0, mean2_im - mean_im * mean_im));
    const double relative_rms =
        (std::abs(mean_im) > 1.0e-15) ? rms / std::abs(mean_im) : 0.0;
    printf(
        " Dispersion: min=%0.3e max=%0.3e mean=%0.3e rms=%0.3e rel.rms=%0.3e\n",
        min_im, max_im, mean_im, rms, relative_rms);
  }

  const ElasticHelicityAmplitudes optical =
      GetMatrixRuntime().HelicityAmplitudes(0.0, 0, 0);
  const double spin_flip_fraction =
      (UnpolarizedHelicityAmpSquared(optical) > 0.0)
          ? (gra::math::abs2(optical.phi5) +
             gra::math::abs2(optical.phi5_first)) /
                UnpolarizedHelicityAmpSquared(optical)
          : 0.0;
  printf(" Matrix diagnostic: max sigma(S)=%0.12f, forward phi5 "
         "fraction=%0.3e\n",
         GetMatrixRuntime().MaxSingularValue(), spin_flip_fraction);
  std::cout << std::endl;
}

// Construct the positive cut Pomeron multiplicity distribution
// Effective opacities diagonalize the Pomeron survival operator as
// S_P^dagger S_P |a> = exp(-lambda_a) |a>
// Each eigenmode follows P_a(m)=exp(-lambda_a) lambda_a^m/m!
// The incoming Good Walker weights give P(m,b)=sum_a w_a P_a(m,b)
// This preserves the exact no cut probability of the matrix unitarization
void MEikonal::S3InitCutPomerons() {
  std::cout << "MEikonal::S3InitCutPomerons:" << std::endl;

  if (!matrix_runtime) {
    throw std::logic_error(
        "MEikonal::S3InitCutPomerons: eikonal is not initialized");
  }

  const double coordinate_min =
      Numerics.logBT ? std::log(Numerics.MinBT) : Numerics.MinBT;
  const double coordinate_max =
      Numerics.logBT ? std::log(Numerics.MaxBT) : Numerics.MaxBT;
  const double step = (coordinate_max - coordinate_min) /
                      static_cast<double>(Numerics.NumberBT);

  cut_b = matrix_runtime->ImpactParameterNodes();
  P_array.assign(MCUT, std::vector<double>(cut_b.size(), 0.0));
  std::vector<double> exact_positive_array(cut_b.size(), 0.0);
  std::vector<double> omitted_tail_array(cut_b.size(), 0.0);

  MAX_P_array = 0.0;
  MAX_P_cut = 0.0;
  double minimum_opacity = std::numeric_limits<double>::infinity();

  for (const auto &node : indices(cut_b)) {
    const double b = cut_b[node];

    // Use d^2b=2 pi b db and db=b dx when x=log b
    const double jacobian = 2.0 * gra::math::PI * (Numerics.logBT ? b * b : b);

    const auto spectrum = matrix_runtime->CutSpectrum(b);
    const auto &opacity = spectrum.eigenvalues;
    const auto &state_weight = spectrum.incoming_weights;

    if (opacity.empty() || opacity.size() != state_weight.size()) {
      throw std::runtime_error(
          "MEikonal::S3InitCutPomerons: invalid opacity spectrum");
    }

    for (const auto &mode : indices(opacity)) {
      if (std::isnan(opacity[mode]) || opacity[mode] < 0.0 ||
          !std::isfinite(state_weight[mode]) || state_weight[mode] < 0.0) {
        throw std::runtime_error(
            "MEikonal::S3InitCutPomerons: invalid spectral opacity/weight");
      }
      minimum_opacity = std::min(minimum_opacity, opacity[mode]);
    }
    const double weight_sum = gra::Sum(state_weight);
    if (!(weight_sum > 0.0) || !std::isfinite(weight_sum) ||
        std::abs(weight_sum - 1.0) > 1.0e-8) {
      throw std::runtime_error(
          "MEikonal::S3InitCutPomerons: spectral weights do not sum to one");
    }

    double captured_positive = 0.0;
    for (std::size_t cuts = 1; cuts < MCUT; ++cuts) {
      double probability = 0.0;
      for (const auto &mode : indices(opacity)) {
        const double lambda = opacity[mode];
        // A completely absorbed mode has no finite-cut probability
        if (lambda > 0.0 && std::isfinite(lambda)) {
          probability += state_weight[mode] *
                         std::exp(static_cast<double>(cuts) * std::log(lambda) - lambda -
                                  std::lgamma(static_cast<double>(cuts) + 1.0));
        }
      }
      captured_positive += probability;
      P_array[cuts][node] = jacobian * probability;
    }

    double exact_positive = 0.0;
    for (const auto &mode : indices(opacity)) {
      exact_positive += state_weight[mode] * -std::expm1(-opacity[mode]);
    }
    exact_positive_array[node] = jacobian * exact_positive;
    omitted_tail_array[node] =
        jacobian * std::max(0.0, exact_positive - captured_positive);
  }

  min_opacity_eigenvalue = minimum_opacity;

  P_cut.assign(MCUT, 0.0);
  double captured_normalization = 0.0;
  for (std::size_t cuts = 1; cuts < MCUT; ++cuts) {
    P_cut[cuts] = gra::math::CS13Integral(P_array[cuts], step);
    captured_normalization += P_cut[cuts];
  }

  const double exact_positive =
      gra::math::CS13Integral(exact_positive_array, step);
  const double omitted_tail = gra::math::CS13Integral(omitted_tail_array, step);

  // A transparent Pomeron has no positive cuts to normalize or sample
  if (std::fpclassify(exact_positive) == FP_ZERO) {
    return;
  }

  if (!(captured_normalization > 0.0) ||
      !std::isfinite(captured_normalization) || !(exact_positive > 0.0) ||
      !std::isfinite(exact_positive)) {
    throw std::runtime_error(
        "MEikonal::S3InitCutPomerons: invalid cut normalization");
  }

  MAX_P_array = 0.0;
  for (std::size_t cuts = 1; cuts < MCUT; ++cuts) {
    P_cut[cuts] /= captured_normalization;
    MAX_P_cut = std::max(MAX_P_cut, P_cut[cuts]);
    // Sample the same discrete impact parameter measure used by the marginal integral
    for (const auto &node : indices(P_array[cuts])) {
      P_array[cuts][node] *= math::CompositeWeight(node, Numerics.NumberBT);
      MAX_P_array = std::max(MAX_P_array, P_array[cuts][node]);
    }
    printf("P_cut[m=%2zu] = %0.7f\n", cuts, P_cut[cuts]);
  }

  const double sum = gra::Sum(P_cut);
  double average = 0.0;
  for (std::size_t cuts = 1; cuts < MCUT; ++cuts) {
    average += static_cast<double>(cuts) * P_cut[cuts];
  }

  printf("--------------------------\n");
  printf("P_cut[SUM ] = %0.7f\n", sum);
  printf("Mean number of cut Pomerons: <m | m>0> = %0.3f\n", average);
  printf("Minimum effective opacity = %0.6e\n", minimum_opacity);
  printf("Omitted m >= %u tail fraction = %0.6e\n\n", MCUT,
         omitted_tail / exact_positive);
}

// Project the eikonal onto the six independent proton helicity amplitudes
// The Good Walker indices select the exclusive physical final state
ElasticHelicityAmplitudes GetElasticHelicityAmplitudes(const MEikonal &eikonal,
                                                       const double kt2,
                                                       const std::size_t f1,
                                                       const std::size_t f2) {
  if (kt2 < 0.0) {
    throw std::invalid_argument(
        "GetElasticHelicityAmplitudes: kt2 must be non-negative");
  }
  return eikonal.GetMatrixRuntime().HelicityAmplitudes(kt2, f1, f2);
}

// Separate the physical helicity amplitudes into crossing even and odd parts
// Separation is performed after the nonlinear coupled channel unitarization
ElasticCrossingHelicityAmplitudes
GetElasticCrossingHelicityAmplitudes(const MEikonal &eikonal, const double kt2,
                                     const std::size_t f1,
                                     const std::size_t f2) {
  if (kt2 < 0.0) {
    throw std::invalid_argument(
        "GetElasticCrossingHelicityAmplitudes: kt2 must be non-negative");
  }
  return eikonal.GetMatrixRuntime().CrossingHelicityAmplitudes(kt2, f1, f2);
}

// Select one Jacob Wick proton elastic helicity amplitude
static std::complex<double>
ElasticHelicityComponentValue(const ElasticHelicityAmplitudes &amplitude,
                              const ElasticHelicityComponent component) {
  switch (component) {
  case ElasticHelicityComponent::Phi1:
    return amplitude.phi1;
  case ElasticHelicityComponent::Phi2:
    return amplitude.phi2;
  case ElasticHelicityComponent::Phi3:
    return amplitude.phi3;
  case ElasticHelicityComponent::Phi4:
    return amplitude.phi4;
  case ElasticHelicityComponent::Phi5:
    return amplitude.phi5;
  case ElasticHelicityComponent::Phi5First:
    return amplitude.phi5_first;
  }
  throw std::logic_error(
      "ElasticHelicityComponentValue: unknown helicity component");
}

// Reconstruct the 4 by 4 proton pair helicity matrix at transfer azimuth phi
// Jacob Wick phases exp(i DeltaLambda phi) restore rotational covariance
ProtonHelicityMatrix
BuildProtonHelicityMatrix(const ElasticHelicityAmplitudes &amplitude,
                          const double azimuth) {
  if (!std::isfinite(azimuth)) {
    throw std::invalid_argument(
        "BuildProtonHelicityMatrix: azimuth must be finite");
  }
  ProtonHelicityMatrix matrix;
  const auto transition = CanonicalProtonHelicityTransitions();
  const std::complex<double> z = std::exp(std::complex<double>(0.0, azimuth));
  const std::complex<double> z2 = z * z;
  const std::array<std::complex<double>, 5> harmonic_phase = {
      std::conj(z2), std::conj(z), std::complex<double>(1.0, 0.0), z, z2};
  for (const auto &index : indices(transition)) {
    const auto &term = transition[index];
    matrix[index] = term.sign *
                    ElasticHelicityComponentValue(amplitude, term.component) *
                    harmonic_phase[term.azimuth_harmonic + 2];
  }
  return matrix;
}

// Rotate a reference plane pair helicity matrix by Jacob Wick phases
ProtonHelicityMatrix
RotateProtonHelicityMatrix(const ProtonHelicityMatrix &reference,
                           const double azimuth) {
  if (!std::isfinite(azimuth)) {
    throw std::invalid_argument(
        "RotateProtonHelicityMatrix: azimuth must be finite");
  }
  ProtonHelicityMatrix matrix;
  const auto transition = CanonicalProtonHelicityTransitions();
  const std::complex<double> z = std::exp(std::complex<double>(0.0, azimuth));
  const std::complex<double> z2 = z * z;
  const std::array<std::complex<double>, 5> harmonic_phase = {
      std::conj(z2), std::conj(z), std::complex<double>(1.0, 0.0), z, z2};
  for (const auto &index : indices(transition)) {
    matrix[index] = reference[index] *
                    harmonic_phase[transition[index].azimuth_harmonic + 2];
  }
  return matrix;
}

// Build coherent Coulomb nuclear interference on the physical elastic channel
void MEikonal::InitializeElasticCNI(const double min_abs_t,
                                    const double max_abs_t) {
  if (!S3INIT || !matrix_runtime) {
    throw std::logic_error(
        "MEikonal::InitializeElasticCNI requires initialized matrix eikonal "
        "mode");
  }
  if (!std::isfinite(min_abs_t) || !std::isfinite(max_abs_t) ||
      min_abs_t <= 0.0 || max_abs_t <= min_abs_t) {
    throw std::invalid_argument(
        "MEikonal::InitializeElasticCNI requires 0 < min_abs_t < max_abs_t");
  }
  const auto &strong_q2 = matrix_runtime->MomentumTransferNodes();
  if (strong_q2.empty() || max_abs_t > strong_q2.back() * (1.0 + 1.0e-12)) {
    throw std::invalid_argument(
        "MEikonal::InitializeElasticCNI: max_abs_t exceeds the hadronic "
        "momentum table");
  }
  if (INITIALSTATE.size() != 2) {
    throw std::logic_error(
        "MEikonal::InitializeElasticCNI: exactly two beams are required");
  }
  elastic_cni_runtime =
      MElasticCNI::Build(*matrix_runtime, soft_model, s, INITIALSTATE,
                         Numerics.CNI, Numerics.logBT, min_abs_t, max_abs_t,
                         model_tune->Structure());
}

// Compute the coherent strong plus electromagnetic proton helicity matrix
ProtonHelicityMatrix
MEikonal::PhysicalElasticHelicityMatrix(const M4Vec &p1_in, const M4Vec &p2_in,
                                        const M4Vec &p1_out,
                                        const M4Vec &p2_out) const {
  if (!elastic_cni_runtime) {
    throw std::logic_error(
        "MEikonal::PhysicalElasticHelicityMatrix: CNI is not initialized");
  }
  return elastic_cni_runtime->PhysicalHelicityMatrix(p1_in, p2_in, p1_out,
                                                     p2_out);
}

// Extract the six independent proton helicity screening operators
MEikonalMatrix::PairHelicityBank
MEikonal::PairScreeningHelicityBank(const double kt2) const {
  const auto spin = PairScreeningSpinBank(kt2);
  const auto value = [&spin](const ElasticHelicityComponent component) {
    const auto transition = ElasticHelicityReferenceTransition(component);
    return spin[ElasticHelicityReferenceIndex(component)] *
           std::complex<double>(transition.sign, 0.0);
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

// Compute all spin transitions of the radial Good Walker screening operator
MEikonalMatrix::PairSpinBank
MEikonal::PairScreeningSpinBank(const double kt2) const {
  return GetMatrixRuntime().PairScreeningSpinBank(kt2);
}

// Compute the unpolarized norm of the pair-helicity amplitude matrix
// <|A|^2> = (1/4) sum_(lambda_i,lambda_f) |A_fi|^2
double
UnpolarizedHelicityAmpSquared(const ElasticHelicityAmplitudes &amplitude) {
  const double value =
      0.25 * gra::SquaredNorm(BuildProtonHelicityMatrix(amplitude));
  if (!std::isfinite(value) || value < 0.0) {
    throw std::runtime_error(
        "UnpolarizedHelicityAmpSquared: non-finite amplitude squared");
  }
  return value;
}

} // namespace gra
