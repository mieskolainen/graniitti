// Tensor Pomeron amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// [REFERENCE: Ewerz, Maniatis, Nachtmann, arxiv.org/abs/1309.3478]
// [REFERENCE: Bolz et al., arxiv.org/abs/1409.8483]
// [REFERENCE: Lebiodowicz, Nachtmann, Szczurek, arxiv.org/abs/1601.04537]
// [REFERENCE: LNS, arxiv.org/abs/1606.05126]
// [REFERENCE: LNS, arxiv.org/abs/1801.03902]
// [REFERENCE: LNS, arxiv.org/abs/1901.11490]
// [REFERENCE: Lebiedowicz et al., arxiv.org/abs/2008.07452]

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tensor/MTensorPhoto.h"
#include "Graniitti/Tensor/MTensorPomeron.h"
#include "Graniitti/Tensor/MTensorSpin3.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Spin/MHelicityBasis.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

// FTensor algebra
#include "FTensor.hpp"

// LOOP MACROS
#define FOR_EACH_2(X)       \
  for (const auto &u : X) { \
    for (const auto &v : X) {
#define FOR_EACH_2_END \
  }                    \
  }

#define FOR_EACH_3(X)         \
  for (const auto &u : X) {   \
    for (const auto &v : X) { \
      for (const auto &k : X) {
#define FOR_EACH_3_END \
  }                    \
  }                    \
  }

#define FOR_EACH_4(X)           \
  for (const auto &u : X) {     \
    for (const auto &v : X) {   \
      for (const auto &k : X) { \
        for (const auto &l : X) {
#define FOR_EACH_4_END \
  }                    \
  }                    \
  }                    \
  }

#define FOR_EACH_5(X)             \
  for (const auto &u : X) {       \
    for (const auto &v : X) {     \
      for (const auto &k : X) {   \
        for (const auto &l : X) { \
          for (const auto &r : X) {
#define FOR_EACH_5_END \
  }                    \
  }                    \
  }                    \
  }                    \
  }

#define FOR_EACH_6(X)               \
  for (const auto &u : X) {         \
    for (const auto &v : X) {       \
      for (const auto &k : X) {     \
        for (const auto &l : X) {   \
          for (const auto &r : X) { \
            for (const auto &s : X) {
#define FOR_EACH_6_END \
  }                    \
  }                    \
  }                    \
  }                    \
  }                    \
  }

// SpinorStates and stored proton pair rows both place negative helicity first
#define FOR_PP_HELICITY                                          \
  for (const auto ha : gra::spin::BinaryHelicityIndices()) {     \
    for (const auto hb : gra::spin::BinaryHelicityIndices()) {   \
      for (const auto h1 : gra::spin::BinaryHelicityIndices()) { \
        for (const auto h2 : gra::spin::BinaryHelicityIndices()) {
#define FOR_PP_HELICITY_END \
  }                         \
  }                         \
  }                         \
  }

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::zi;

using FTensor::Tensor0;
using FTensor::Tensor1;
using FTensor::Tensor2;
using FTensor::Tensor3;
using FTensor::Tensor4;

namespace gra {

namespace {

// Evaluate the shared continuum off-shell form factor for one Tensor vertex
double TensorOffshell(const TensorVertexParam &vertex, const double q2, const double mass) {
  return regge::FormFactor(q2, pow2(mass), vertex.ff_offshell);
}

// Require one immutable SOFT model before Tensor parameter initialization
SoftModelPtr RequireTensorSoftModel(SoftModelPtr model) {
  if (model == nullptr) { throw std::invalid_argument("MTensorPomeron: null SOFT model"); }
  return model;
}

// Validate configured beam particles for one Tensor amplitude mode
void ValidateTensorBeamParticles(const LORENTZSCALAR &lts, const MTensorPomeronMode mode) {
  std::size_t configured = 0;
  std::size_t targets = 0;
  const bool photo = mode == MTensorPomeronMode::Photo || mode == MTensorPomeronMode::Resonance;
  for (const MParticle *beam : {&lts.beam1, &lts.beam2}) {
    if (beam->pdg == 0) { continue; }
    ++configured;
    const bool proton  = beam->pdg == PDG::PDG_p;
    const bool ion = photo && nuclear::IsNuclearPDG(beam->pdg);
    const bool lepton = photo && nuclear::IsChargedLepton(beam->pdg);
    if (!proton && !ion && !lepton) { throw std::invalid_argument("MTensorPomeron: unsupported beam for the selected mode"); }
    if (proton || ion) { ++targets; }
  }
  if (photo && configured == 2 && targets == 0) {
    throw std::invalid_argument("MTensorPomeron: photoproduction requires a proton or nuclear target");
  }
}

// Convert one stored Tensor fermion spin index into physical helicity
double TensorFermionHelicity(const std::size_t index) {
  return static_cast<double>(spin::BinaryHelicityLabelX2(index)) / 2.0;
}

// Compute the collider-section phase of one exact Tensor forward current
std::complex<double> TensorForwardColliderCurrentPhase(const ForwardLegState &state, const std::size_t initial_helicity,
                                                       const std::size_t final_helicity) {
  const bool  upper             = state.leg == ForwardBeamLeg::Upper;
  const M4Vec collider_transfer = upper ? state.outgoing - state.incoming : state.incoming - state.outgoing;
  const MDirac::FermionKind fermion_kind =
      state.emitter.pdg < 0 ? MDirac::FermionKind::Antiparticle : MDirac::FermionKind::Particle;
  return MDirac::ElasticSpinHalfColliderCurrentPhase(fermion_kind, state.Index(),
                                                     TensorFermionHelicity(initial_helicity),
                                                     TensorFermionHelicity(final_helicity), collider_transfer.Phi());
}

// Require the V to two-pseudoscalar current to be conserved
void RequireConservedVectorDecay(const M4Vec &k1, const M4Vec &k2, const std::string &context) {
  const double difference     = k1.M2() - k2.M2();
  const double mass_scale     = std::max({1.0, std::abs(k1.M2()), std::abs(k2.M2())});
  const auto   k1_components  = k1.Contravariant();
  const auto   k2_components  = k2.Contravariant();
  double       momentum_scale = 1.0;
  for (const auto &mu : indices(k1_components)) {
    momentum_scale = std::max({momentum_scale, math::pow2(k1_components[mu]), math::pow2(k2_components[mu])});
  }
  const double tolerance =
      std::max(1.0e-10 * mass_scale, 512.0 * std::numeric_limits<double>::epsilon() * momentum_scale);
  if (!std::isfinite(difference) || !std::isfinite(momentum_scale) || std::abs(difference) > tolerance) {
    throw AmplitudeFailure(context +
                                ": V -> PS PS requires equal daughter masses for the "
                                "conserved-current propagator of arXiv:1309.3478");
  }
}

// Contract one direct vector exchange after both Pomeron legs are reduced
Tensor2<std::complex<double>, 4, 4> TensorVectorDirectExchange(
    const Tensor2<std::complex<double>, 4, 4> &upper, const Tensor2<std::complex<double>, 4, 4> &lower,
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &upper_vertex,
    const Tensor2<std::complex<double>, 4, 4>       &vector_propagator,
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &lower_vertex, const std::complex<double> scale) {
  constexpr std::size_t               dimension = 4;
  Tensor2<std::complex<double>, 4, 4> upper_block;
  Tensor2<std::complex<double>, 4, 4> lower_block;
  Tensor2<std::complex<double>, 4, 4> propagated;
  Tensor2<std::complex<double>, 4, 4> out;

  for (std::size_t rho1 = 0; rho1 < dimension; ++rho1) {
    for (std::size_t rho3 = 0; rho3 < dimension; ++rho3) {
      upper_block(rho1, rho3) = 0.0;
      for (std::size_t alpha = 0; alpha < dimension; ++alpha) {
        for (std::size_t beta = 0; beta < dimension; ++beta) {
          upper_block(rho1, rho3) += upper(alpha, beta) * upper_vertex(rho1, rho3, alpha, beta);
        }
      }
    }
  }
  for (std::size_t rho4 = 0; rho4 < dimension; ++rho4) {
    for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
      lower_block(rho4, rho2) = 0.0;
      for (std::size_t alpha = 0; alpha < dimension; ++alpha) {
        for (std::size_t beta = 0; beta < dimension; ++beta) {
          lower_block(rho4, rho2) += lower(alpha, beta) * lower_vertex(rho4, rho2, alpha, beta);
        }
      }
    }
  }
  for (std::size_t rho3 = 0; rho3 < dimension; ++rho3) {
    for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
      propagated(rho3, rho2) = 0.0;
      for (std::size_t rho1 = 0; rho1 < dimension; ++rho1) {
        propagated(rho3, rho2) += upper_block(rho1, rho3) * vector_propagator(rho1, rho2);
      }
    }
  }
  for (std::size_t rho3 = 0; rho3 < dimension; ++rho3) {
    for (std::size_t rho4 = 0; rho4 < dimension; ++rho4) {
      out(rho3, rho4) = 0.0;
      for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
        out(rho3, rho4) += propagated(rho3, rho2) * lower_block(rho4, rho2);
      }
      out(rho3, rho4) *= scale;
    }
  }
  return out;
}

// Contract one crossed vector exchange after both Pomeron legs are reduced
Tensor2<std::complex<double>, 4, 4> TensorVectorCrossedExchange(
    const Tensor2<std::complex<double>, 4, 4> &upper, const Tensor2<std::complex<double>, 4, 4> &lower,
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &upper_vertex,
    const Tensor2<std::complex<double>, 4, 4>       &vector_propagator,
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &lower_vertex, const std::complex<double> scale) {
  constexpr std::size_t               dimension = 4;
  Tensor2<std::complex<double>, 4, 4> upper_block;
  Tensor2<std::complex<double>, 4, 4> lower_block;
  Tensor2<std::complex<double>, 4, 4> propagated;
  Tensor2<std::complex<double>, 4, 4> out;

  for (std::size_t rho4 = 0; rho4 < dimension; ++rho4) {
    for (std::size_t rho1 = 0; rho1 < dimension; ++rho1) {
      upper_block(rho4, rho1) = 0.0;
      for (std::size_t alpha = 0; alpha < dimension; ++alpha) {
        for (std::size_t beta = 0; beta < dimension; ++beta) {
          upper_block(rho4, rho1) += upper(alpha, beta) * upper_vertex(rho4, rho1, alpha, beta);
        }
      }
    }
  }
  for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
    for (std::size_t rho3 = 0; rho3 < dimension; ++rho3) {
      lower_block(rho2, rho3) = 0.0;
      for (std::size_t alpha = 0; alpha < dimension; ++alpha) {
        for (std::size_t beta = 0; beta < dimension; ++beta) {
          lower_block(rho2, rho3) += lower(alpha, beta) * lower_vertex(rho2, rho3, alpha, beta);
        }
      }
    }
  }
  for (std::size_t rho4 = 0; rho4 < dimension; ++rho4) {
    for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
      propagated(rho4, rho2) = 0.0;
      for (std::size_t rho1 = 0; rho1 < dimension; ++rho1) {
        propagated(rho4, rho2) += upper_block(rho4, rho1) * vector_propagator(rho1, rho2);
      }
    }
  }
  for (std::size_t rho3 = 0; rho3 < dimension; ++rho3) {
    for (std::size_t rho4 = 0; rho4 < dimension; ++rho4) {
      out(rho3, rho4) = 0.0;
      for (std::size_t rho2 = 0; rho2 < dimension; ++rho2) {
        out(rho3, rho4) += propagated(rho4, rho2) * lower_block(rho2, rho3);
      }
      out(rho3, rho4) *= scale;
    }
  }
  return out;
}

using TensorLegBank      = std::map<int, std::array<Tensor2<std::complex<double>, 4, 4>, 4>>;
using VectorLegBank      = std::map<int, std::array<Tensor1<std::complex<double>, 4>, 4>>;
using TensorExchangePair = std::array<int, 2>;

// Collect physical Tensor amplitudes resolved by ordered exchange pair
class TensorPairBank {
 public:
  // Append one block while retaining zeros for absent exchange pairs
  void Append(const std::map<TensorExchangePair, std::vector<std::complex<double>>> &block, const std::size_t count) {
    for (const auto &[pair, amplitude] : block) {
      (void)pair;
      if (amplitude.size() != count) { throw AmplitudeFailure("TensorPairBank::Append: exchange block size changed"); }
    }
    for (auto &[pair, amplitude] : channel) {
      (void)pair;
      amplitude.resize(size + count, 0.0);
    }
    for (const auto &[pair, amplitude] : block) {
      auto [entry, inserted] = channel.try_emplace(pair, std::vector<std::complex<double>>(size + count, 0.0));
      (void)inserted;
      for (const auto &i : indices(amplitude)) { entry->second[size + i] += amplitude[i]; }
    }
    size += count;
  }

  // Append one scalar amplitude for every contributing exchange pair
  void Append(const std::map<TensorExchangePair, std::complex<double>> &block) {
    std::map<TensorExchangePair, std::vector<std::complex<double>>> vector;
    for (const auto &[pair, amplitude] : block) { vector.emplace(pair, std::vector<std::complex<double>>{amplitude}); }
    Append(vector, 1);
  }

  // Add one coherent bank with an optional common complex scale
  void Add(const TensorPairBank &other, const std::complex<double> scale = 1.0) {
    if (size == 0) { size = other.size; }
    if (other.size != size) { throw AmplitudeFailure("TensorPairBank::Add: amplitude size changed"); }
    for (const auto &[pair, amplitude] : other.channel) {
      auto [entry, inserted] = channel.try_emplace(pair, std::vector<std::complex<double>>(size, 0.0));
      (void)inserted;
      for (const auto &i : indices(amplitude)) { entry->second[i] += scale * amplitude[i]; }
    }
  }

  // Repeat the complete first hard-spin block a fixed number of times
  void Repeat(const std::size_t copies) {
    const std::size_t block_size = size;
    for (auto &[pair, amplitude] : channel) {
      (void)pair;
      const auto block = amplitude;
      amplitude.reserve((copies + 1) * block_size);
      for (std::size_t copy = 0; copy < copies; ++copy) {
        amplitude.insert(amplitude.end(), block.cbegin(), block.cend());
      }
    }
    size *= copies + 1;
  }

  // Compute the common physical amplitude dimension
  std::size_t Size() const noexcept { return size; }

  // Compute every ordered exchange contribution
  const auto &Channels() const noexcept { return channel; }

 private:
  std::size_t                                                     size = 0;
  std::map<TensorExchangePair, std::vector<std::complex<double>>> channel;
};

// Store one summed resonance block and its ordered exchange decomposition
struct TensorChannelSum {
  std::vector<std::complex<double>>                               total;
  std::map<TensorExchangePair, std::vector<std::complex<double>>> channel;
};

// Store one projected rank-one or rank-two soft exchange leg
struct TensorSoftLeg {
  int                                 rank = 0;
  Tensor2<std::complex<double>, 4, 4> tensor;
  Tensor1<std::complex<double>, 4>    vector;
};

// Store one ordered pseudoscalar continuum exchange block
struct TensorScalarContinuumBlock {
  std::array<int, 2>                  exchange = {0, 0};
  Tensor2<std::complex<double>, 4, 4> t_upper;
  Tensor2<std::complex<double>, 4, 4> t_lower;
  Tensor2<std::complex<double>, 4, 4> u_upper;
  Tensor2<std::complex<double>, 4, 4> u_lower;
  Tensor1<std::complex<double>, 4>    t_upper_vector;
  Tensor1<std::complex<double>, 4>    t_lower_vector;
  Tensor1<std::complex<double>, 4>    u_upper_vector;
  Tensor1<std::complex<double>, 4>    u_lower_vector;
  std::complex<double>                t_scale = 0.0;
  std::complex<double>                u_scale = 0.0;
};

// Store one ordered outer pair with one internal vector transfer
struct TensorVectorBlock {
  std::array<int, 2>                        exchange = {0, 0};
  Tensor4<std::complex<double>, 4, 4, 4, 4> t_upper;
  Tensor2<std::complex<double>, 4, 4>       t_propagator;
  Tensor4<std::complex<double>, 4, 4, 4, 4> t_lower;
  Tensor4<std::complex<double>, 4, 4, 4, 4> u_upper;
  Tensor2<std::complex<double>, 4, 4>       u_propagator;
  Tensor4<std::complex<double>, 4, 4, 4, 4> u_lower;
  double                                    t_scale = 1.0;
  double                                    u_scale = 1.0;
};

// Compute the published phi-pair threshold suppression for Odderon exchange
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 99, 094034 (2019), Eq. (3.51)]
double TensorVectorThreshold(const double s34, const double vector_mass) {
  const double threshold = 4.0 * pow2(vector_mass);
  if (!(s34 > threshold)) { return 0.0; }
  return 1.0 - std::exp((threshold - s34) / threshold);
}

// Build all configured outer-pair and internal vector transfer blocks
std::vector<TensorVectorBlock> TensorVectorBlocks(
    const MTensorPomeron &tp, const MTensorPomeronParam &tensor_param, const int vector_pdg,
    const M4Vec &p3, const M4Vec &p4, const M4Vec &pt, const M4Vec &pu, const double s34) {
  std::vector<TensorVectorBlock> out;
  const double vector_mass = tensor_param.FindVector(vector_pdg).mass;
  const auto &pairs = tensor_param.exchange.FindContinuumPairs(vector_pdg, vector_pdg);
  for (const auto &pair : pairs) {
    for (const auto &exchange : MTensorExchangeModel::OrderedPairs(pair)) {
      if (tensor_param.exchange.FindExchange(exchange[0]).rank != 2 ||
          tensor_param.exchange.FindExchange(exchange[1]).rank != 2) {
        continue;
      }
      const auto transfers =
          tensor_param.exchange.FindActiveTransfers(exchange[0], exchange[1], vector_pdg, vector_pdg);
      for (const int transfer_pdg : transfers) {
        TensorVectorBlock block;
        block.exchange = exchange;
        if (std::abs(transfer_pdg) == std::abs(vector_pdg)) {
          const auto &upper  = tensor_param.exchange.FindVertex(exchange[0], vector_pdg);
          const auto &lower  = tensor_param.exchange.FindVertex(exchange[1], vector_pdg);
          block.t_upper      = tp.iG_Tvv(pt, -p3, exchange[0], vector_pdg);
          block.t_propagator = tp.iD_V(pt, vector_mass, s34, vector_pdg);
          block.t_lower      = tp.iG_Tvv(p4, pt, exchange[1], vector_pdg);
          block.u_upper      = tp.iG_Tvv(p4, pu, exchange[0], vector_pdg);
          block.u_propagator = tp.iD_V(pu, vector_mass, s34, vector_pdg);
          block.u_lower      = tp.iG_Tvv(pu, -p3, exchange[1], vector_pdg);
          block.t_scale =
              TensorOffshell(upper, pt.M2(), vector_mass) * TensorOffshell(lower, pt.M2(), vector_mass);
          block.u_scale =
              TensorOffshell(upper, pu.M2(), vector_mass) * TensorOffshell(lower, pu.M2(), vector_mass);
        } else {
          if (tensor_param.exchange.FindExchange(transfer_pdg).rank != 1) { continue; }
          const auto &upper  = tensor_param.exchange.FindVertex(exchange[0], vector_pdg, transfer_pdg);
          const auto &lower  = tensor_param.exchange.FindVertex(exchange[1], vector_pdg, transfer_pdg);
          block.t_upper      = tp.iG_Pvv(pt, -p3, upper.g_tensor[0], upper.g_tensor[1], upper.ff_transfer);
          block.t_propagator = tp.iD_VExchange(transfer_pdg, s34, pt.M2());
          block.t_lower      = tp.iG_Pvv(p4, pt, lower.g_tensor[0], lower.g_tensor[1], lower.ff_transfer);
          block.u_upper      = tp.iG_Pvv(p4, pu, upper.g_tensor[0], upper.g_tensor[1], upper.ff_transfer);
          block.u_propagator = tp.iD_VExchange(transfer_pdg, s34, pu.M2());
          block.u_lower      = tp.iG_Pvv(pu, -p3, lower.g_tensor[0], lower.g_tensor[1], lower.ff_transfer);
          block.t_scale =
              regge::FormFactor(pt.M2(), 0.0, upper.ff_offshell) * regge::FormFactor(pt.M2(), 0.0, lower.ff_offshell);
          block.u_scale =
              regge::FormFactor(pu.M2(), 0.0, upper.ff_offshell) * regge::FormFactor(pu.M2(), 0.0, lower.ff_offshell);
          if (upper.threshold || lower.threshold) {
            const double threshold = TensorVectorThreshold(s34, vector_mass);
            block.t_scale *= threshold;
            block.u_scale *= threshold;
          }
        }
        out.push_back(block);
      }
    }
  }
  return out;
}

// Sum all configured strong rank-two resonance channels coherently
template <typename Evaluation>
TensorChannelSum SumTensorResonanceChannels(const RES_TENSOR_MODEL &model, const TensorLegBank &upper,
                                            const TensorLegBank &lower, const std::size_t upper_helicity,
                                            const std::size_t lower_helicity, Evaluation &&evaluate) {
  TensorChannelSum out;
  bool             found = false;
  for (const auto &channel : model.channels) {
    if (channel.exchange[0] == 22 || channel.exchange[1] == 22) { continue; }
    if (channel.active_ready && channel.active_g_tensor.empty()) { continue; }
    for (const auto &ordered : MTensorExchangeModel::OrderedPairs(channel.exchange)) {
      const std::vector<std::complex<double>> amplitude =
          evaluate(upper.at(ordered[0])[upper_helicity], lower.at(ordered[1])[lower_helicity], channel);
      if (out.total.empty()) {
        out.total.resize(amplitude.size(), 0.0);
      } else if (out.total.size() != amplitude.size()) {
        throw AmplitudeFailure("Tensor resonance channels have incompatible helicity bases");
      }
      auto [entry, inserted] =
          out.channel.try_emplace(ordered, std::vector<std::complex<double>>(amplitude.size(), 0.0));
      (void)inserted;
      for (const auto &i : indices(out.total)) {
        out.total[i] += amplitude[i];
        entry->second[i] += amplitude[i];
      }
      found = true;
    }
  }
  if (!found) { throw AmplitudeFailure("Tensor resonance has no strong exchange channel"); }
  return out;
}

// Compute one normalized elastic Good Walker source for a Tensor exchange leg
std::vector<std::complex<double>> TensorElasticSource(const SoftModel            &model,
                                                      const MTensorExchangeModel &exchange_model,
                                                      const int exchange_pdg, const double t) {
  const auto                       &good_walker = model.GoodWalker();
  std::vector<std::complex<double>> proton(good_walker.ProtonVector().begin(), good_walker.ProtonVector().end());
  if (exchange_pdg == PDG::PDG_gamma) { return proton; }
  const SoftExchangeId exchange = exchange_model.SoftId(exchange_pdg, model);
  auto                 source   = model.ResidueMatrix(exchange, t) * proton;
  std::complex<double> physical = 0.0;
  for (const auto &i : indices(proton)) { physical += std::conj(proton[i]) * source[i]; }
  if (!std::isfinite(physical.real()) || !std::isfinite(physical.imag()) || std::abs(physical) <= 1.0e-14) {
    throw AmplitudeFailure("TensorElasticSource: vanishing physical SOFT residue");
  }
  gra::Scale(source, 1.0 / physical);
  return source;
}

// Build the exact elastic pair source while preserving the Tensor Born limit
PairSource TensorElasticPair(const LORENTZSCALAR &lts, const SoftModel &model,
                             const MTensorExchangeModel &exchange_model, const TensorExchangePair &exchange) {
  const ScreeningMetadata layout = ProtonSpinLayout(lts);
  const auto              upper  = TensorElasticSource(model, exchange_model, exchange[0], lts.t1);
  const auto              lower  = TensorElasticSource(model, exchange_model, exchange[1], lts.t2);
  const std::size_t       n      = model.GoodWalker().ChannelCount();
  PairSource              pair;
  pair.amplitude = MMatrix<std::complex<double>>(layout.spin_rows, n * n, 0.0);
  for (std::size_t row = 0; row < layout.spin_rows; ++row) {
    for (const auto &i : indices(upper)) {
      for (const auto &j : indices(lower)) { pair.amplitude[row][i * n + j] = upper[i] * lower[j]; }
    }
  }
  return pair;
}

// Attach channel-resolved Tensor Born terms to the proton Good Walker amplitude
void StoreTensorPair(LORENTZSCALAR &lts, const SoftModelPtr &model, const MTensorExchangeModel &exchange_model,
                     const TensorPairBank &bank, const std::string &context) {
  if (lts.hamp.metadata.amplitude_type != ScreeningAmplitudeType::GoodWalker) { return; }
  if (model == nullptr || lts.excite1 || lts.excite2 || bank.Channels().empty() || bank.Size() != lts.hamp.size()) {
    throw AmplitudeFailure(context + ": invalid Tensor pair source");
  }
  lts.proton_good_walker.reset();
  for (const auto &[exchange, amplitude] : bank.Channels()) {
    const PairSource              pair = TensorElasticPair(lts, *model, exchange_model, exchange);
    MMatrix<std::complex<double>> central(amplitude.size(), 1, 0.0);
    for (const auto &row : indices(amplitude)) { central[row][0] = amplitude[row]; }
    AddGoodWalker(lts, {pair}, central, 0, model);
  }
  const auto projected = ProjectGoodWalker(*lts.proton_good_walker);
  if (projected.size() != lts.hamp.size()) {
    throw AmplitudeFailure(context + ": Tensor Born projection size changed");
  }
  for (const auto &i : indices(projected)) {
    const double scale = std::max({1.0, std::abs(projected[i]), std::abs(lts.hamp[i])});
    if (std::abs(projected[i] - lts.hamp[i]) > 2.0e-10 * scale) {
      throw AmplitudeFailure(context + ": Tensor Born projection changed");
    }
  }
}

}  // namespace

// Build one immutable Tensor Pomeron process definition for a selected mode
std::shared_ptr<const amplitude::ProcessDefinition> MTensorPomeron::ProcessDefinitionFor(MTensorPomeronMode mode) {
  return std::make_shared<amplitude::AnalyticProcess>(
      "MTENSORPOMERON", tensor::ProcessName(mode), tensor::ProcessPattern(mode), MTensorPomeron::DirectDecayStructure(),
      [mode](const std::vector<MDecayBranch> &tree) { return tensor::ProcessAccepts(tree, mode); },
      [mode](const LORENTZSCALAR &lts) { return tensor::ProcessDecayStructure(mode, lts); });
}

// Initialize immutable Tensor Pomeron parameters before worker process copies
void MTensorPomeron::InitializeParameters(const MProcessSetup &setup) {
  if (setup.model_tune == nullptr) {
    throw std::invalid_argument("MTensorPomeron::InitializeParameters: missing model tune");
  }
  MModelCache &cache = RequireModelCache(setup.lts.model_cache, setup.model_tune, "MTensorPomeron initialization");
  (void)GetTensorParam(cache, setup.lts.PDG, setup.lts.process.RESONANCES);
}

// Constructor
MTensorPomeron::MTensorPomeron(gra::LORENTZSCALAR &lts, MModelTunePtr tune,
                               std::shared_ptr<const amplitude::ProcessDefinition> definition,
                               const MTensorPomeronMode                            mode)
    : amplitude::ProcessFamily(std::move(definition)),
      model_tune(RequireModelCache(lts.model_cache, tune, "MTensorPomeron").TunePtr()),
      soft_model(RequireTensorSoftModel(model_tune->Soft())),
      tensor_param_handle(GetTensorParam(*lts.model_cache, lts.PDG, lts.process.RESONANCES)),
      tensor_param(*tensor_param_handle) {
  ValidateTensorBeamParticles(lts, mode);
  CalcRTensor();  // Pre-Calculate tensors
}

// Compute the forward proton spin layout owned by one Tensor amplitude mode
bool MTensorPomeron::ForwardNoFlip(const MTensorPomeronMode mode, const gra::LORENTZSCALAR &lts) const {
  switch (mode) {
    case MTensorPomeronMode::QED:
      return lts.process.FORWARD_NOFLIP;
    case MTensorPomeronMode::Generic:
    case MTensorPomeronMode::Resonance:
    case MTensorPomeronMode::Continuum:
    case MTensorPomeronMode::ResonanceContinuum:
      return tensor_param.FORWARD_NOFLIP;
    case MTensorPomeronMode::Photo:
      return tensor_param.FORWARD_NOFLIP;
  }
  throw std::invalid_argument("MTensorPomeron::ForwardNoFlip: unknown mode");
}

namespace {

// Transform a lower-index complex axial current into the central rest frame
std::array<std::complex<double>, 4> CentralRestFrameCurrent(const LORENTZSCALAR                    &lts,
                                                            const Tensor1<std::complex<double>, 4> &current_lower) {
  M4Vec real_upper(-std::real(current_lower(1)), -std::real(current_lower(2)), -std::real(current_lower(3)),
                   std::real(current_lower(0)));
  M4Vec imag_upper(-std::imag(current_lower(1)), -std::imag(current_lower(2)), -std::imag(current_lower(3)),
                   std::imag(current_lower(0)));
  real_upper = kinematics::BoostToRestFrame(real_upper, lts.pfinal[0], "MTensorPomeron::ME3 axial real current");
  imag_upper = kinematics::BoostToRestFrame(imag_upper, lts.pfinal[0], "MTensorPomeron::ME3 axial imaginary current");

  std::array<std::complex<double>, 4> out = {};
  out[0]                                  = {real_upper.E(), imag_upper.E()};
  out[1]                                  = {-real_upper.Px(), -imag_upper.Px()};
  out[2]                                  = {-real_upper.Py(), -imag_upper.Py()};
  out[3]                                  = {-real_upper.Pz(), -imag_upper.Pz()};
  return out;
}

// Compute whether both top-level vector branches carry explicit two-body decays
bool HasVectorCascadeDecay(const std::vector<MDecayBranch> &tree) {
  return tree.size() == 2 && tree[0].p.spinX2 == 2 && tree[1].p.spinX2 == 2 && !tree[0].legs.empty() &&
         !tree[1].legs.empty();
}

// Build the covariant vector propagator and decay tensor block for V V -> 4PS
Tensor2<std::complex<double>, 4, 4> TensorCascadeDecayBlock(const MTensorPomeron            &tp,
                                                            const std::vector<MDecayBranch> &tree);

// Compute coherent scalar-resonance V V cascade decay contraction
std::complex<double> TensorCascadeScalarDecayAmplitude(const MTensorPomeron &tp, gra::LORENTZSCALAR &lts, double M0,
                                                       const std::vector<double> &couplings,
                                                       const regge::FFParam      &ff_decay);

// Compute coherent tensor-resonance V V cascade decay tensor
Tensor2<std::complex<double>, 4, 4> TensorCascadeSpin2DecayTensor(const MTensorPomeron &tp, gra::LORENTZSCALAR &lts,
                                                                  double M0, const std::vector<double> &couplings,
                                                                  const regge::FFParam &ff_decay);

// Build the scalar-resonance decay block from central event kinematics
TensorDecayState BuildScalarDecay(const MTensorPomeron &tp, LORENTZSCALAR &lts, const PARAM_RES &res,
                                  const regge::FFParam &ff_decay) {
  TensorDecayState out;
  out.type        = TensorResonanceType::Scalar;
  out.propagator  = tp.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);
  const M4Vec &p3 = lts.decaytree[0].p4;
  const M4Vec &p4 = lts.decaytree[1].p4;
  if (lts.decaytree[0].p.spinX2 == 0 && lts.decaytree[1].p.spinX2 == 0) {
    out.scalar.push_back(tp.iG_f0ss(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], ff_decay));
    return out;
  }
  if (HasVectorCascadeDecay(lts.decaytree)) {
    out.scalar.push_back(TensorCascadeScalarDecayAmplitude(tp, lts, res.p.mass, res.hel_decay.g_decay_TP, ff_decay));
    return out;
  }
  const auto vertex =
      tp.iG_f0vv(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.g_decay_TP[1], ff_decay);
  out.scalar = tp.MassiveSpin1PolSum(vertex, p3, p4);
  return out;
}

// Build the vector-resonance decay current from central event kinematics
TensorDecayState BuildVectorDecay(const MTensorPomeron &tp, const LORENTZSCALAR &lts, const PARAM_RES &res,
                                  const regge::FFParam &ff_decay) {
  TensorDecayState out;
  out.type        = TensorResonanceType::Vector;
  // Fix the particle-antiparticle order independently of the decay-tree order
  const bool reversed = lts.decaytree[0].p.pdg < 0 &&
                        lts.decaytree[0].p.pdg == -lts.decaytree[1].p.pdg;
  const M4Vec &p3 = lts.decaytree[reversed ? 1 : 0].p4;
  const M4Vec &p4 = lts.decaytree[reversed ? 0 : 1].p4;
  out.vector = tp.VectorDecay(p3, p4, res.p.mass, res.p.width, res.p.pdg,
                               res.hel_decay.g_decay_TP[0], ff_decay);
  return out;
}

// Build the axial-resonance decay matrix from central event kinematics
TensorDecayState BuildAxialDecay(const MTensorPomeron &tp, LORENTZSCALAR &lts, const PARAM_RES &res,
                                 const regge::FFParam &ff_decay) {
  TensorDecayState out;
  out.type            = TensorResonanceType::AxialVector;
  out.propagator      = tp.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);
  PARAM_RES axial_res = res;
  spin::DecayAmp(lts, axial_res, "CM");
  axial_res.decay_f = axial_res.decay_f * regge::MassFF(lts.m2, pow2(res.p.mass), ff_decay);
  out.axial = std::move(axial_res.decay_f);
  return out;
}

// Build the pseudoscalar-resonance decay block from central event kinematics
TensorDecayState BuildPseudoscalarDecay(const MTensorPomeron &tp, const LORENTZSCALAR &lts, const PARAM_RES &res,
                                        const regge::FFParam &ff_decay) {
  TensorDecayState out;
  out.type       = TensorResonanceType::Pseudoscalar;
  out.propagator = tp.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);
  const M4Vec &p3     = lts.decaytree[0].p4;
  const M4Vec &p4     = lts.decaytree[1].p4;
  const auto   vertex = tp.iG_psvv(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], ff_decay);
  out.scalar          = tp.MasslessSpin1PolSum(vertex, p3, p4);
  return out;
}

// Build the tensor-resonance decay tensors from central event kinematics
TensorDecayState BuildSpin2Decay(const MTensorPomeron &tp, LORENTZSCALAR &lts, const PARAM_RES &res,
                                 const regge::FFParam &ff_decay) {
  TensorDecayState out;
  out.type                          = TensorResonanceType::Tensor;
  const M4Vec           &p3         = lts.decaytree[0].p4;
  const M4Vec           &p4         = lts.decaytree[1].p4;
  const auto             propagator = tp.iD_TMES(lts.pfinal[0], res.p.mass, res.p.width, true);
  FTensor::Index<'a', 4> mu;
  FTensor::Index<'b', 4> nu;
  FTensor::Index<'c', 4> rho;
  FTensor::Index<'d', 4> sigma;

  if (lts.decaytree[0].p.spinX2 == 0 && lts.decaytree[1].p.spinX2 == 0) {
    const auto vertex = tp.iG_f2psps(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], ff_decay);
    Tensor2<std::complex<double>, 4, 4> block;
    block(mu, nu) = propagator(mu, nu, rho, sigma) * vertex(rho, sigma);
    out.tensor.push_back(block);
    return out;
  }

  const bool massive_vectors = lts.decaytree[0].p.spinX2 == 2 && lts.decaytree[1].p.spinX2 == 2 &&
                               lts.decaytree[0].p.pdg != PDG::PDG_gamma && lts.decaytree[1].p.pdg != PDG::PDG_gamma;
  if (massive_vectors) {
    if (HasVectorCascadeDecay(lts.decaytree)) {
      const auto vertex = TensorCascadeSpin2DecayTensor(tp, lts, res.p.mass, res.hel_decay.g_decay_TP, ff_decay);
      Tensor2<std::complex<double>, 4, 4> block;
      block(mu, nu) = propagator(mu, nu, rho, sigma) * vertex(rho, sigma);
      out.tensor.push_back(block);
      return out;
    }
    const auto vertex =
        tp.iG_f2vv(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.g_decay_TP[1], ff_decay);
    for (const auto &decay : tp.MassiveSpin1PolSum(vertex, p3, p4)) {
      Tensor2<std::complex<double>, 4, 4> block;
      block(mu, nu) = propagator(mu, nu, rho, sigma) * decay(rho, sigma);
      out.tensor.push_back(block);
    }
    return out;
  }

  const auto vertex =
      tp.iG_f2yy(p3, p4, res.p.mass, res.hel_decay.g_decay_TP[0], res.hel_decay.g_decay_TP[1], ff_decay);
  for (const auto &decay : tp.MasslessSpin1PolSum(vertex, p3, p4)) {
    Tensor2<std::complex<double>, 4, 4> block;
    block(mu, nu) = propagator(mu, nu, rho, sigma) * decay(rho, sigma);
    out.tensor.push_back(block);
  }
  return out;
}

// Build one resonance decay representation selected by its quantum numbers
TensorDecayState BuildTensorDecay(const MTensorPomeron &tp, LORENTZSCALAR &lts, const PARAM_RES &res,
                                  const TensorResonanceType type) {
  const regge::FFParam &ff_decay = res.hel_decay.ff_decay;
  switch (type) {
    case TensorResonanceType::Scalar:
      return BuildScalarDecay(tp, lts, res, ff_decay);
    case TensorResonanceType::Pseudoscalar:
      return BuildPseudoscalarDecay(tp, lts, res, ff_decay);
    case TensorResonanceType::Vector:
      return BuildVectorDecay(tp, lts, res, ff_decay);
    case TensorResonanceType::AxialVector:
      return BuildAxialDecay(tp, lts, res, ff_decay);
    case TensorResonanceType::Tensor:
      return BuildSpin2Decay(tp, lts, res, ff_decay);
    case TensorResonanceType::Spin3: {
      const auto positive = lts.decaytree[0].p.chargeX3 > 0 ? 0U : 1U;
      TensorDecayState out;
      out.type = type;
      out.spin3 = tensor::Spin3Decay(lts.decaytree[positive].p4, lts.decaytree[1 - positive].p4);
      out.propagator = res.hel_decay.g_decay_TP[0] * regge::MassFF(lts.m2, pow2(res.p.mass), ff_decay) /
                       std::complex<double>(lts.m2 - pow2(res.p.mass), res.p.mass * res.p.width);
      return out;
    }
  }
  throw std::invalid_argument("MTensorPomeron::ME3: unknown resonance type");
}

// Compute a decay block bound to this Born event or build an uncached fallback
const TensorDecayState &TensorDecayForEvent(const MTensorPomeron &tp, LORENTZSCALAR &lts, const std::string &name,
                                            const PARAM_RES &res, const TensorResonanceType type,
                                            TensorDecayState &local) {
  auto                   &cache = lts.amplitude.tensor_decay;
  const TensorDecayState *out   = nullptr;
  if (lts.screening.active && lts.amplitude.central.Active()) {
    out = &cache.Get(lts.amplitude.central, name, "MTensorPomeron::ME3");
  } else {
    local = BuildTensorDecay(tp, lts, res, type);
    if (lts.amplitude.central.Active()) {
      out = &cache.Store(lts.amplitude.central, name, std::move(local), "MTensorPomeron::ME3");
    } else {
      out = &local;
    }
  }
  if (out->type != type) { throw std::logic_error("MTensorPomeron::ME3: cached resonance decay type mismatch"); }
  return *out;
}

}  // namespace

// Compute decay coupling constant for resonance (M,Gamma) with decay daughter
// mass mf BR being the branching ratio BR = Width_partial / Width_total
//
// Decay: Mother (spin = 0/1/2/3) -> scalar or pseudoscalar daughters
//
double MTensorPomeron::GDecay(int J, double M, double Gamma, double mf, double BR, double symmetry) {
  const double S0           = 1.0;  // Should be set as the same scale as in decay amplitudes iG[]
  const double partialWidth = Gamma * BR * symmetry;

  const double P = sqrt(1 - 4 * pow2(mf / M));

  if (J == 0) {
    return sqrt(partialWidth / (1.0 / (16 * PI * M) * pow2(S0) * P));
  } else if (J == 1) {
    return sqrt(partialWidth / (M / (192 * PI) * math::pow3(P)));
  } else if (J == 2) {
    return sqrt(partialWidth / (M / (480 * PI) * pow2(M / S0) * math::pow5(P)));
  } else if (J == 3) {
    // The STF rank-three vertex has spin-averaged norm (2/35) M^6 beta^6
    return sqrt(partialWidth * 280.0 * PI / (math::pow5(M) * std::pow(P, 7)));
  } else {
    throw std::invalid_argument("MTensorPomeron::GDecay: Unknown input spin J = " + std::to_string(J));
  }
}

// Compute the coupling for i g/(2 S0) epsilon(mu,nu,rho,sigma) p1^rho p2^sigma
// including the identical-photon phase-space factor
//
double MTensorPomeron::GDecayPseudoscalarGammaGamma(double M, double Gamma, double BR) {
  const double S0 = 1.0;
  if (!(M > 0.0) || !(Gamma >= 0.0) || !(BR >= 0.0)) {
    throw std::invalid_argument(
        "MTensorPomeron::GDecayPseudoscalarGammaGamma: "
        "invalid mass, width or BR");
  }
  return std::sqrt(Gamma * BR * 256.0 * PI * pow2(S0) / math::pow3(M));
}

// Build the shared covariant vector production kernels including the physical decay
MTensorPomeron::VectorKernels MTensorPomeron::VectorProduction(
    const LORENTZSCALAR &lts, const PARAM_RES &res,
    const Tensor1<std::complex<double>, 4> &decay_current) const {
  const double M0 = res.p.mass;
  const double Gamma = res.p.width;
  const auto phase = std::polar(1.0, res.TP.phi);
  const auto decay_phase = tensor_param.use_zeta ? std::polar(1.0, res.hel_decay.zeta) : std::complex<double>(1.0, 0.0);

  // Gamma-Vector coupling
  const Tensor2<std::complex<double>, 4, 4> iGyV = iG_yV(0, res.p.pdg);

  // Vector-meson propagators
  const bool                                INDEX_UP          = true;
  const bool                                CONSERVED_CURRENT = true;
  const Tensor2<std::complex<double>, 4, 4> iDV_1 =
      iD_VMES(lts.q1, M0, Gamma, res.p.pdg, INDEX_UP, CONSERVED_CURRENT);
  const Tensor2<std::complex<double>, 4, 4> iDV_2 =
      iD_VMES(lts.q2, M0, Gamma, res.p.pdg, INDEX_UP, CONSERVED_CURRENT);

  // Build the effective photoproduction kernels
  //
  // K(mu,alpha,beta) =
  //   Gamma_gammaV(mu,nu) D_V(nu,rho)
  //   D_V(rho',nu') Gamma_VPSPS(nu')
  //   Gamma_PVV(rho',rho,alpha,beta)
  //
  // Direct and VMD-mixing paths are added coherently into the same kernel
  // before the proton helicity loop
  VectorKernels out;
  auto &kernel_yT = out.yT;
  auto &kernel_Ty = out.Ty;
  auto &kernel_VT = out.VT;
  auto &kernel_TV = out.TV;

  auto AddPhotoproductionKernel =
      [&](Tensor3<std::complex<double>, 4, 4, 4> &kernel, const Tensor2<std::complex<double>, 4, 4> &iGyV_path,
          const Tensor2<std::complex<double>, 4, 4>       &iDV_path,
          const Tensor4<std::complex<double>, 4, 4, 4, 4> &iGPvv_path, const double vector_ff) {
        Tensor2<std::complex<double>, 4, 4> photon_transfer;
        for (const auto &mu : LI) {
          for (const auto &rho : LI) {
            photon_transfer(mu, rho) = 0.0;
            for (const auto &nu : LI) {
              photon_transfer(mu, rho) += vector_ff * iGyV_path(mu, nu) * iDV_path(nu, rho);
            }
          }
        }

        Tensor3<std::complex<double>, 4, 4, 4> central;
        for (const auto &rho : LI) {
          for (const auto &alpha : LI) {
            for (const auto &beta : LI) {
              central(rho, alpha, beta) = 0.0;
              for (const auto &rho2 : LI) {
                central(rho, alpha, beta) +=
                    phase * decay_phase * decay_current(rho2) * iGPvv_path(rho2, rho, alpha, beta);
              }
            }
          }
        }

        for (const auto &mu : LI) {
          for (const auto &alpha : LI) {
            for (const auto &beta : LI) {
              for (const auto &rho : LI) {
                kernel(mu, alpha, beta) += photon_transfer(mu, rho) * central(rho, alpha, beta);
              }
            }
          }
        }
      };

  // Add one direct vector-exchange and tensor-exchange production kernel
  const auto AddStrongVectorKernel = [&](Tensor3<std::complex<double>, 4, 4, 4> &kernel, const Tensor4<std::complex<double>, 4, 4, 4, 4> &vertex, const double production_ff) {
    for (const auto &rho : LI) {
      for (const auto &alpha : LI) {
        for (const auto &beta : LI) {
          kernel(rho, alpha, beta) = 0.0;
          for (const auto &rho2 : LI) {
            kernel(rho, alpha, beta) +=
                phase * production_ff * decay_phase * decay_current(rho2) * vertex(rho2, rho, alpha, beta);
          }
        }
      }
    }
  };

  // Build each configured photon-tensor or vector-tensor channel
  for (const auto &channel : res.TP.channels) {
    if (channel.active_ready && channel.active_g_tensor.empty() && channel.active_vmd.empty()) { continue; }
    const regge::FFParam &channel_transfer = channel.ff_transfer;
    const regge::FFParam &channel_prod     = channel.ff_prod;
    // Dress both vector legs independently of the decay vertex
    // [REFERENCE: Lebiedowicz et al., arXiv:2508.06334v2, Eqs. (2.26)-(2.29)]
    const double production_ff = regge::MassFF(lts.m2, pow2(M0), channel_prod);
    const bool            first_photon     = channel.exchange[0] == 22;
    const bool            second_photon    = channel.exchange[1] == 22;
    if (first_photon && second_photon) { throw AmplitudeFailure("Vector resonance does not support a two-photon TP channel"); }
    if (!first_photon && !second_photon) {
      const auto &first  = tensor_param.exchange.FindExchange(channel.exchange[0]);
      const auto &second = tensor_param.exchange.FindExchange(channel.exchange[1]);
      if (first.rank == second.rank) { throw AmplitudeFailure("Vector resonance requires vector-tensor strong fusion"); }
      const int vector_pdg = first.rank == 1 ? channel.exchange[0] : channel.exchange[1];
      const int tensor_pdg = first.rank == 2 ? channel.exchange[0] : channel.exchange[1];
      if (!channel.VMD_MIXING.empty()) { throw AmplitudeFailure("VMD_MIXING is only valid for photon-tensor fusion"); }
      auto        &upper         = kernel_VT[{vector_pdg, tensor_pdg}];
      auto        &lower         = kernel_TV[{tensor_pdg, vector_pdg}];
      auto         upper_vertex  = iG_Pvv(lts.pfinal[0], lts.q1, channel.g_tensor[0], channel.g_tensor[1], channel_transfer);
      auto         lower_vertex  = iG_Pvv(lts.pfinal[0], lts.q2, channel.g_tensor[0], channel.g_tensor[1], channel_transfer);
      const double upper_factor  = regge::TransferFF(lts.q1.M2(), channel_transfer);
      const double lower_factor  = regge::TransferFF(lts.q2.M2(), channel_transfer);
      FOR_EACH_4(LI);
      upper_vertex(u, v, k, l) *= upper_factor;
      lower_vertex(u, v, k, l) *= lower_factor;
      FOR_EACH_4_END;
      AddStrongVectorKernel(upper, upper_vertex, production_ff);
      AddStrongVectorKernel(lower, lower_vertex, production_ff);
      continue;
    }
    const int exchange_pdg = first_photon ? channel.exchange[1] : channel.exchange[0];
    tensor_param.exchange.FindExchange(exchange_pdg);
    const bool hera = lts.process.PHOTO_DISSOCIATION == DissociationType::Hera &&
                      tensor_param.exchange.FindExchange(exchange_pdg).type == TensorExchangeType::Pomeron;
    const bool lower_diss = hera && ResolveForwardLegState(lts, ForwardBeamLeg::Lower).IsExcited();
    const bool upper_diss = hera && ResolveForwardLegState(lts, ForwardBeamLeg::Upper).IsExcited();

    auto [upper_kernel, upper_inserted] = kernel_yT.try_emplace(exchange_pdg);
    auto [lower_kernel, lower_inserted] = kernel_Ty.try_emplace(exchange_pdg);
    if (upper_inserted) {
      FOR_EACH_3(LI);
      upper_kernel->second(u, v, k) = 0.0;
      FOR_EACH_3_END;
    }
    if (lower_inserted) {
      FOR_EACH_3(LI);
      lower_kernel->second(u, v, k) = 0.0;
      FOR_EACH_3_END;
    }

    if (!channel.active_ready || !channel.active_g_tensor.empty()) {
      const auto vertex_1 =
          iG_Pvv(lts.pfinal[0], lts.q1, channel.g_tensor[0], channel.g_tensor[1], channel_transfer, lower_diss);
      const auto vertex_2 =
          iG_Pvv(lts.pfinal[0], lts.q2, channel.g_tensor[0], channel.g_tensor[1], channel_transfer, upper_diss);
      const double vector_ff_1 = regge::MassFF(lts.q1.M2(), pow2(M0), channel_prod);
      const double vector_ff_2 = regge::MassFF(lts.q2.M2(), pow2(M0), channel_prod);
      AddPhotoproductionKernel(upper_kernel->second, iGyV, iDV_1, vertex_1, production_ff * vector_ff_1);
      AddPhotoproductionKernel(lower_kernel->second, iGyV, iDV_2, vertex_2, production_ff * vector_ff_2);
    }

    // Add channel-local VMD mixing paths coherently
    // [REFERENCE: Lebiedowicz et al., arXiv:1911.01909, Appendix B]
    const auto AddVMDMixing = [&](const RES_TENSOR_VMD_MIXING &mixing) {
      const auto mixing_iGyV  = iG_yV(0, mixing.pdg);
      const auto mixing_iDV_1 = iD_VMD(lts.q1, mixing.pdg);
      const auto mixing_iDV_2 = iD_VMD(lts.q2, mixing.pdg);
      const auto mixing_vertex_1 =
          iG_Pvv(lts.pfinal[0], lts.q1, mixing.coupling[0], mixing.coupling[1], channel_transfer, lower_diss);
      const auto mixing_vertex_2 =
          iG_Pvv(lts.pfinal[0], lts.q2, mixing.coupling[0], mixing.coupling[1], channel_transfer, upper_diss);
      const double mixing_mass = tensor_param.FindVMD(mixing.pdg).mass;
      const double mixing_ff_1 = regge::MassFF(lts.q1.M2(), pow2(mixing_mass), channel_prod);
      const double mixing_ff_2 = regge::MassFF(lts.q2.M2(), pow2(mixing_mass), channel_prod);
      AddPhotoproductionKernel(upper_kernel->second, mixing_iGyV, mixing_iDV_1, mixing_vertex_1, production_ff * mixing_ff_1);
      AddPhotoproductionKernel(lower_kernel->second, mixing_iGyV, mixing_iDV_2, mixing_vertex_2, production_ff * mixing_ff_2);
    };
    if (channel.active_ready) {
      for (const std::size_t i : channel.active_vmd) { AddVMDMixing(channel.VMD_MIXING[i]); }
    } else {
      for (const auto &mixing : channel.VMD_MIXING) {
        if (!ActiveProductionCoupling(mixing.coupling[0]) && !ActiveProductionCoupling(mixing.coupling[1])) {
          continue;
        }
        AddVMDMixing(mixing);
      }
    }
  }

  return out;
}

// 2 -> 3 amplitudes
//
// return value: matrix element squared with helicities summed and averaged over
// lts.hamp:     individual complex helicity amplitudes
//
double MTensorPomeron::ME3(gra::LORENTZSCALAR &lts, const MTensorPomeronMode mode) const {
  lts.proton_good_walker.reset();
  // Tensor amplitudes use their model-local forward-spin approximation switch
  const bool forward_noflip = ForwardNoFlip(mode, lts);

  // Free Lorentz indices [second parameter denotes the range of index]
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;

  // Kinematics
  const M4Vec pa = lts.pbeam1;
  const M4Vec pb = lts.pbeam2;

  const M4Vec                   p1          = lts.pfinal[1];
  const M4Vec                   p2          = lts.pfinal[2];
  // [REFERENCE: Lebiedowicz, Nachtmann and Szczurek, arXiv:2506.04846, Eq. (2.11)]
  const double                  nu1x2       = (pa + p1) * lts.pfinal[0];
  const double                  nu2x2       = (pb + p2) * lts.pfinal[0];
  if (!(nu1x2 > 0.0) || !(nu2x2 > 0.0) || !std::isfinite(nu1x2) || !std::isfinite(nu2x2)) {
    throw AmplitudeFailure("MTensorPomeron::ME3: invalid crossing-symmetric Regge energy");
  }
  const ForwardLegState         upper_state = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState         lower_state = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const TensorForwardSourceBank source_bank = TensorForwardSourceBank::Hadronic(upper_state, lower_state);
  const bool upper_strong = mode != MTensorPomeronMode::Photo || std::abs(upper_state.emitter.pdg) == PDG::PDG_p;
  const bool lower_strong = mode != MTensorPomeronMode::Photo || std::abs(lower_state.emitter.pdg) == PDG::PDG_p;

  // ------------------------------------------------------------------

  // Spinors (2 helicities)
  bool exact_spinors_required = !forward_noflip;
  if (!exact_spinors_required) {
    for (const auto &x : lts.process.RESONANCES) {
      const auto &p = x.second.p;
      if (p.chargeX3 == 0 && (p.spinX2 == 2 || p.spinX2 == 6) && p.P == -1 && p.C == -1) {
        exact_spinors_required = true;
        break;
      }
    }
  }

  std::array<MDirac::Spinor, 2> u_a{};
  std::array<MDirac::Spinor, 2> u_b{};
  std::array<MDirac::Spinor, 2> ubar_1{};
  std::array<MDirac::Spinor, 2> ubar_2{};
  if (exact_spinors_required) {
    u_a = SpinorStates(pa, "u");
    u_b = SpinorStates(pb, "u");
    if (!upper_state.IsExcited()) { ubar_1 = SpinorStates(p1, "ubar"); }
    if (!lower_state.IsExcited()) { ubar_2 = SpinorStates(p2, "ubar"); }
  }

  // ------------------------------------------------------------------

  // Cache each projected rank-two leg across all resonance channels
  std::set<int> tensor_pdgs;
  std::set<int> vector_pdgs;
  for (const auto &[name, resonance] : lts.process.RESONANCES) {
    (void)name;
    for (const auto &channel : resonance.TP.channels) {
      if (channel.active_ready && channel.active_g_tensor.empty() && channel.active_vmd.empty()) { continue; }
      for (const int pdg : channel.exchange) {
        if (pdg != 22) {
          const auto &exchange = tensor_param.exchange.FindExchange(pdg);
          if (exchange.rank == 2) {
            tensor_pdgs.insert(pdg);
          } else {
            vector_pdgs.insert(pdg);
          }
        }
      }
    }
  }

  TensorLegBank tensor_leg_1;
  TensorLegBank tensor_leg_2;
  for (const int exchange_pdg : tensor_pdgs) {
    const std::complex<double> factor_1 = TensorPropagatorFactor(exchange_pdg, nu1x2, lts.t1);
    const std::complex<double> factor_2 = TensorPropagatorFactor(exchange_pdg, nu2x2, lts.t2);
    auto                      &upper    = tensor_leg_1[exchange_pdg];
    auto                      &lower    = tensor_leg_2[exchange_pdg];
    if (forward_noflip) {
      if (upper_strong) {
        const auto upper_leg = PomeronPropagatorCurrent(iG_TForwardHE(upper_state, exchange_pdg), factor_1);
        upper[0]             = upper_leg;
        upper[3]             = upper_leg;
      }
      if (lower_strong) {
        const auto lower_leg = PomeronPropagatorCurrent(iG_TForwardHE(lower_state, exchange_pdg), factor_2);
        lower[0]             = lower_leg;
        lower[3]             = lower_leg;
      }
      continue;
    }
    if (upper_strong) {
      for (const auto ha : spin::BinaryHelicityIndices()) {
        for (const auto h1 : spin::BinaryHelicityIndices()) {
          const std::size_t index = spin::BinaryPairHelicityIndex(ha, h1);
          upper[index] =
              PomeronPropagatorCurrent(iG_TForward(upper_state, exchange_pdg, ubar_1[h1], u_a[ha], ha, h1), factor_1);
        }
      }
    }
    if (lower_strong) {
      for (const auto hb : spin::BinaryHelicityIndices()) {
        for (const auto h2 : spin::BinaryHelicityIndices()) {
          const std::size_t index = spin::BinaryPairHelicityIndex(hb, h2);
          lower[index] =
              PomeronPropagatorCurrent(iG_TForward(lower_state, exchange_pdg, ubar_2[h2], u_b[hb], hb, h2), factor_2);
        }
      }
    }
  }

  VectorLegBank vector_leg_1;
  VectorLegBank vector_leg_2;
  for (const int exchange_pdg : vector_pdgs) {
    auto &upper = vector_leg_1[exchange_pdg];
    auto &lower = vector_leg_2[exchange_pdg];
    if (forward_noflip) {
      if (upper_strong) {
        const auto upper_leg =
            VectorPropagatorCurrent(iG_VForwardHE(upper_state, exchange_pdg), exchange_pdg, nu1x2, lts.t1);
        upper[0] = upper_leg;
        upper[3] = upper_leg;
      }
      if (lower_strong) {
        const auto lower_leg =
            VectorPropagatorCurrent(iG_VForwardHE(lower_state, exchange_pdg), exchange_pdg, nu2x2, lts.t2);
        lower[0] = lower_leg;
        lower[3] = lower_leg;
      }
      continue;
    }
    if (upper_strong) {
      for (const auto ha : spin::BinaryHelicityIndices()) {
        for (const auto h1 : spin::BinaryHelicityIndices()) {
          const std::size_t index = spin::BinaryPairHelicityIndex(ha, h1);
          upper[index] = VectorPropagatorCurrent(iG_VForward(upper_state, exchange_pdg, ubar_1[h1], u_a[ha], ha, h1),
                                                 exchange_pdg, nu1x2, lts.t1);
        }
      }
    }
    if (lower_strong) {
      for (const auto hb : spin::BinaryHelicityIndices()) {
        for (const auto h2 : spin::BinaryHelicityIndices()) {
          const std::size_t index = spin::BinaryPairHelicityIndex(hb, h2);
          lower[index] = VectorPropagatorCurrent(iG_VForward(lower_state, exchange_pdg, ubar_2[h2], u_b[hb], hb, h2),
                                                 exchange_pdg, nu2x2, lts.t2);
        }
      }
    }
  }
  // ------------------------------------------------------------------

  // Build each resonance in its exact helicity basis before the coherent sum
  lts.hamp.clear();
  TensorPairBank pair_bank;

  // Store one strong central-helicity block in the shared Tensor source basis
  const auto AppendPomeronPomeron = [&](std::vector<std::complex<double>> &target, TensorPairBank &target_bank,
                                        const TensorChannelSum &amplitudes) {
    const auto &component = source_bank.Components(TensorForwardMechanism::PomeronPomeron);
    if (component.size() != 1) { throw AmplitudeFailure("MTensorPomeron::ME3: invalid PP forward source basis"); }
    std::map<TensorExchangePair, std::vector<std::complex<double>>> block;
    for (const auto &[exchange, channel] : amplitudes.channel) {
      block.emplace(exchange, std::vector<std::complex<double>>(amplitudes.total.size() * source_bank.Size(), 0.0));
      for (const auto &i : indices(channel)) { block.at(exchange)[i * source_bank.Size() + component[0]] = channel[i]; }
    }
    for (const auto &amplitude : amplitudes.total) {
      std::vector<std::complex<double>> source(source_bank.Size(), 0.0);
      source[component[0]] = amplitude;
      target.insert(target.end(), source.begin(), source.end());
    }
    target_bank.Append(block, amplitudes.total.size() * source_bank.Size());
  };

  // >>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
  // 2. Coherent sum of Resonances (loop over)
  for (auto &x : lts.process.RESONANCES) {
    const PARAM_RES                  &res = x.second;
    std::vector<std::complex<double>> resonance_hamp;
    TensorPairBank                    resonance_bank;

    const TensorResonanceType resonance_type = tensor::ClassifyTensorResonance(res.p);

    // Resonance parameters
    const double               M0    = res.p.mass;
    const std::complex<double> phase = std::polar(1.0, res.TP.phi);
    const std::complex<double> decay_phase =
        tensor_param.use_zeta ? std::polar(1.0, res.hel_decay.zeta) : std::complex<double>(1.0, 0.0);

    TensorDecayState        local_decay;
    const TensorDecayState &decay =
        TensorDecayForEvent(*this, lts, x.first, res, resonance_type, local_decay);

    // Axial production is covariant while its generic decay tree uses g_decay

    // ------------------------------------------------------------------

    // =====================================================================
    // Pomeron-Pomeron-Scalar structure
    //
    //
    if (resonance_type == TensorResonanceType::Scalar) {
      const auto                &eps3eps4 = decay.scalar;
      const std::complex<double> iD       = decay.propagator;

      auto ScalarProductionAmplitudes = [&](const std::size_t upper_helicity, const std::size_t lower_helicity) {
        return SumTensorResonanceChannels(
            res.TP, tensor_leg_1, tensor_leg_2, upper_helicity, lower_helicity,
            [&](const Tensor2<std::complex<double>, 4, 4> &left, const Tensor2<std::complex<double>, 4, 4> &right,
                const RES_TENSOR_CHANNEL &channel) {
              auto vertex = iG_PPS_total(lts.q1, lts.q2, M0, resonance_type, channel.g_tensor, channel);
              FOR_EACH_4(LI);
              vertex(u, v, k, l) *= phase;
              FOR_EACH_4_END;
              Tensor2<std::complex<double>, 4, 4> central;
              central(alpha1, beta1) = vertex(alpha1, beta1, alpha2, beta2) * right(alpha2, beta2);
              const std::complex<double> block = (-zi) * iD * left(alpha1, beta1) * central(alpha1, beta1);

              std::vector<std::complex<double>> amplitudes;
              amplitudes.reserve(eps3eps4.size());
              for (const auto &decay : eps3eps4) { amplitudes.push_back(block * decay_phase * decay); }
              return amplitudes;
            });
      };

      TensorChannelSum noflip_hamp;
      if (forward_noflip) { noflip_hamp = ScalarProductionAmplitudes(0, 0); }

      // Two helicity states for incoming and outgoing protons
      FOR_PP_HELICITY;

      // Apply proton leg helicity conservation / No helicity flip (high energy
      // limit)
      if (forward_noflip && (ha != h1 || hb != h2)) { continue; }

      if (forward_noflip) {
        AppendPomeronPomeron(resonance_hamp, resonance_bank, noflip_hamp);
      } else {
        const std::size_t      index_1    = spin::BinaryPairHelicityIndex(ha, h1);
        const std::size_t      index_2    = spin::BinaryPairHelicityIndex(hb, h2);
        const TensorChannelSum amplitudes = ScalarProductionAmplitudes(index_1, index_2);
        AppendPomeronPomeron(resonance_hamp, resonance_bank, amplitudes);
      }

      FOR_PP_HELICITY_END;

      // =====================================================================
      // Pomeron-Gamma-Vector (Photoproduction of rho, phi ...) structure
      //
      // Should add vector meson mass dependent running of t-form factors
      //
      // p --x--------------- p
      //      *
      //     y *   rho0
      //        *x=====x===== rho0
      //               |
      //               | P
      //               |
      // p ------------x----- p
      //
    } else if (resonance_type == TensorResonanceType::Vector || resonance_type == TensorResonanceType::Spin3) {
      const auto kernels = resonance_type == TensorResonanceType::Spin3 ? Spin3Production(lts, res, decay)
                                                                      : VectorProduction(lts, res, decay.vector);
      const auto &kernel_yT = kernels.yT;
      const auto &kernel_Ty = kernels.Ty;
      const auto &kernel_VT = kernels.VT;
      const auto &kernel_TV = kernels.TV;
      // Keep the strong-exchange currents intact when photo and hadronic channels coexist
      const auto photo_legs = [&](const TensorLegBank &bank, const ForwardLegState &target,
                                  const auto &kernels, double nu, const M4Vec &photon) {
        auto out = bank;
        if (target.IsExcited() && lts.process.PHOTO_DISSOCIATION == DissociationType::Hera) {
          for (const auto &[exchange, kernel] : kernels) {
            (void)kernel;
            if (tensor_param.exchange.FindExchange(exchange).type != TensorExchangeType::Pomeron) { continue; }
            const auto current = PhotoDissCurrent(target, exchange, res.p.pdg, nu, (photon + target.incoming).M2());
            out.at(exchange)[0] = current;
            out.at(exchange)[3] = current;
          }
        }
        return out;
      };
      const auto photo_leg_1 = photo_legs(tensor_leg_1, upper_state, kernel_Ty, nu1x2, lts.q2);
      const auto photo_leg_2 = photo_legs(tensor_leg_2, lower_state, kernel_yT, nu2x2, lts.q1);

      // Final reduced contraction for one proton helicity combination
      auto PhotoproductionAmplitude = [&](const Tensor1<std::complex<double>, 4>       &gamma_current,
                                          const Tensor3<std::complex<double>, 4, 4, 4> &kernel,
                                          const Tensor2<std::complex<double>, 4, 4>    &tensor_leg) {
        std::complex<double> out = 0.0;
        for (const auto &mu : LI) {
          for (const auto &alpha : LI) {
            for (const auto &beta : LI) {
              out += gamma_current(mu) * kernel(mu, alpha, beta) * tensor_leg(alpha, beta);
            }
          }
        }
        return out;
      };

      // Proton electromagnetic currents and reduced Pomeron legs contain only
      // one proton-helicity pair each. Cache them instead of rebuilding them
      // inside the four-helicity Cartesian product
      std::array<std::vector<Tensor1<std::complex<double>, 4>>, 4> gamma_sources_1;
      std::array<std::vector<Tensor1<std::complex<double>, 4>>, 4> gamma_sources_2;
      const TensorPhotonVertex                                     photon_vertex =
          mode == MTensorPomeronMode::Photo && !tensor_param.photo.proton_pauli ? TensorPhotonVertex::Dirac
                                                                                                                         : TensorPhotonVertex::DiracPauli;
      const double photo_q2_1  = -lts.q1.M2();
      const double photo_q2_2  = -lts.q2.M2();
      const bool   use_gamma_1 = mode != MTensorPomeronMode::Photo || (std::isfinite(photo_q2_1) && photo_q2_1 >= 0.0 &&
                                                                     photo_q2_1 <= tensor_param.photo.q2_max);
      const bool   use_gamma_2 = mode != MTensorPomeronMode::Photo || (std::isfinite(photo_q2_2) && photo_q2_2 >= 0.0 &&
                                                                     photo_q2_2 <= tensor_param.photo.q2_max);

      if (forward_noflip) {
        if (use_gamma_1) {
          gamma_sources_1[0] = iG_yForwardSources(lts, upper_state, ubar_1[0], u_a[0], 0, 0, photon_vertex);
          gamma_sources_1[3] = iG_yForwardSources(lts, upper_state, ubar_1[1], u_a[1], 1, 1, photon_vertex);
        }
        if (use_gamma_2) {
          gamma_sources_2[0] = iG_yForwardSources(lts, lower_state, ubar_2[0], u_b[0], 0, 0, photon_vertex);
          gamma_sources_2[3] = iG_yForwardSources(lts, lower_state, ubar_2[1], u_b[1], 1, 1, photon_vertex);
        }
      } else {
        if (use_gamma_1) {
          for (const auto ha : spin::BinaryHelicityIndices()) {
            for (const auto h1 : spin::BinaryHelicityIndices()) {
              const std::size_t index = spin::BinaryPairHelicityIndex(ha, h1);
              gamma_sources_1[index] = iG_yForwardSources(lts, upper_state, ubar_1[h1], u_a[ha], ha, h1, photon_vertex);
            }
          }
        }
        if (use_gamma_2) {
          for (const auto hb : spin::BinaryHelicityIndices()) {
            for (const auto h2 : spin::BinaryHelicityIndices()) {
              const std::size_t index = spin::BinaryPairHelicityIndex(hb, h2);
              gamma_sources_2[index] = iG_yForwardSources(lts, lower_state, ubar_2[h2], u_b[hb], hb, h2, photon_vertex);
            }
          }
        }
      }

      // Two helicity states for incoming and outgoing protons
      FOR_PP_HELICITY;

      // Apply proton leg helicity conservation / No helicity flip (high energy
      // limit)
      if (forward_noflip && (ha != h1 || hb != h2)) { continue; }

      const std::size_t index_1 = spin::BinaryPairHelicityIndex(ha, h1);
      const std::size_t index_2 = spin::BinaryPairHelicityIndex(hb, h2);

      std::vector<std::complex<double>>                               source_amplitude(source_bank.Size(), 0.0);
      std::map<TensorExchangePair, std::vector<std::complex<double>>> source_channel;
      const auto                                                      AccumulateOrdering =
          [&](const TensorForwardMechanism mechanism, const std::vector<Tensor1<std::complex<double>, 4>> &sources,
              const Tensor3<std::complex<double>, 4, 4, 4> &kernel, const Tensor2<std::complex<double>, 4, 4> &pomeron,
              const TensorExchangePair &exchange) {
            if (sources.empty()) { return; }
            const auto &components = source_bank.Components(mechanism);
            if (components.size() != sources.size()) {
              throw AmplitudeFailure("MTensorPomeron::ME3: photon source-bank mismatch");
            }
            auto [entry, inserted] =
                source_channel.try_emplace(exchange, std::vector<std::complex<double>>(source_bank.Size(), 0.0));
            (void)inserted;
            for (const auto &source : indices(sources)) {
              const std::complex<double> amplitude = (-zi) * PhotoproductionAmplitude(sources[source], kernel, pomeron);
              source_amplitude[components[source]] += amplitude;
              entry->second[components[source]] += amplitude;
            }
          };

      // Excited photon and Pomeron transitions have no modeled mixed density,
      // so crossed orderings occupy separate source components when required
      for (const auto &[exchange_pdg, kernel] : kernel_yT) {
        AccumulateOrdering(TensorForwardMechanism::GammaPomeron, gamma_sources_1[index_1], kernel,
                           photo_leg_2.at(exchange_pdg)[index_2], {PDG::PDG_gamma, exchange_pdg});
      }
      for (const auto &[exchange_pdg, kernel] : kernel_Ty) {
        AccumulateOrdering(TensorForwardMechanism::PomeronGamma, gamma_sources_2[index_2], kernel,
                           photo_leg_1.at(exchange_pdg)[index_1], {exchange_pdg, PDG::PDG_gamma});
      }
      const auto &strong_component = source_bank.Components(TensorForwardMechanism::PomeronPomeron);
      if (strong_component.size() != 1) {
        throw AmplitudeFailure("MTensorPomeron::ME3: invalid strong forward source basis");
      }
      for (const auto &[exchange, kernel] : kernel_VT) {
        const std::complex<double> amplitude =
            (-zi) * PhotoproductionAmplitude(vector_leg_1.at(exchange.first)[index_1], kernel,
                                             tensor_leg_2.at(exchange.second)[index_2]);
        source_amplitude[strong_component[0]] += amplitude;
        auto &channel = source_channel[{exchange.first, exchange.second}];
        channel.resize(source_bank.Size(), 0.0);
        channel[strong_component[0]] += amplitude;
      }
      for (const auto &[exchange, kernel] : kernel_TV) {
        const std::complex<double> amplitude =
            (-zi) * PhotoproductionAmplitude(vector_leg_2.at(exchange.second)[index_2], kernel,
                                             tensor_leg_1.at(exchange.first)[index_1]);
        source_amplitude[strong_component[0]] += amplitude;
        auto &channel = source_channel[{exchange.first, exchange.second}];
        channel.resize(source_bank.Size(), 0.0);
        channel[strong_component[0]] += amplitude;
      }
      resonance_hamp.insert(resonance_hamp.end(), source_amplitude.begin(), source_amplitude.end());
      resonance_bank.Append(source_channel, source_bank.Size());

      FOR_PP_HELICITY_END;

      // =====================================================================
      // Pomeron-Pomeron-Axial-vector structure
      //
      //
    } else if (resonance_type == TensorResonanceType::AxialVector) {
      // The two couplings are the (l,S)=(2,2) and (4,4) PPf1 structures
      // [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020)]
      // The reference supplies the production vertex but no differential f1
      // decay vertex
      const std::complex<double> iD = decay.propagator;

      const M4Vec rest_momentum(0.0, 0.0, 0.0, lts.pfinal[0].M());
      const auto  eps_conj = MassiveSpin1States(rest_momentum, "conj", true);

      auto AxialProductionAmplitudes = [&](const std::size_t upper_helicity, const std::size_t lower_helicity) {
        return SumTensorResonanceChannels(
            res.TP, tensor_leg_1, tensor_leg_2, upper_helicity, lower_helicity,
            [&](const Tensor2<std::complex<double>, 4, 4> &left, const Tensor2<std::complex<double>, 4, 4> &right,
                const RES_TENSOR_CHANNEL &channel) {
              const MTensor<std::complex<double>>    vertex = iG_PPA_total(lts.q1, lts.q2, M0, channel);
              Tensor3<std::complex<double>, 4, 4, 4> central;
              for (const auto upper_mu : LI) {
                for (const auto upper_nu : LI) {
                  for (const auto alpha : LI) {
                    central(upper_mu, upper_nu, alpha) = 0.0;
                    for (const auto lower_mu : LI) {
                      for (const auto lower_nu : LI) {
                        central(upper_mu, upper_nu, alpha) +=
                            vertex({upper_mu, upper_nu, lower_mu, lower_nu, alpha}) * right(lower_mu, lower_nu);
                      }
                    }
                  }
                }
              }

              Tensor1<std::complex<double>, 4> current;
              for (const auto alpha : LI) {
                current(alpha) = 0.0;
                for (const auto upper_mu : LI) {
                  for (const auto upper_nu : LI) {
                    current(alpha) += left(upper_mu, upper_nu) * central(upper_mu, upper_nu, alpha);
                  }
                }
              }

              const auto                    rest_current = CentralRestFrameCurrent(lts, current);
              MMatrix<std::complex<double>> production(1, 3, 0.0);
              for (const auto &h : indices(eps_conj)) {
                for (const auto alpha : LI) { production[0][h] += rest_current[alpha] * eps_conj[h](alpha); }
                production[0][h] *= -zi * phase;
              }

              const std::complex<double> decay_coupling =
                  tensor_param.use_zeta ? res.hel_decay.g_decay
                                        : std::complex<double>(std::abs(res.hel_decay.g_decay), 0.0);
              const MMatrix<std::complex<double>> amplitude = (production * decay.axial) * (iD * decay_coupling);
              std::vector<std::complex<double>>   amplitudes;
              amplitudes.reserve(amplitude.size_col());
              for (std::size_t column = 0; column < amplitude.size_col(); ++column) {
                amplitudes.push_back(amplitude[0][column]);
              }
              return amplitudes;
            });
      };

      TensorChannelSum noflip_hamp;
      if (forward_noflip) { noflip_hamp = AxialProductionAmplitudes(0, 0); }

      FOR_PP_HELICITY;

      if (forward_noflip && (ha != h1 || hb != h2)) { continue; }

      if (forward_noflip) {
        AppendPomeronPomeron(resonance_hamp, resonance_bank, noflip_hamp);
      } else {
        const std::size_t      index_1    = spin::BinaryPairHelicityIndex(ha, h1);
        const std::size_t      index_2    = spin::BinaryPairHelicityIndex(hb, h2);
        const TensorChannelSum amplitudes = AxialProductionAmplitudes(index_1, index_2);
        AppendPomeronPomeron(resonance_hamp, resonance_bank, amplitudes);
      }

      FOR_PP_HELICITY_END;

      // =====================================================================
      // Pomeron-Pomeron-Pseudoscalar structure
      //
      //
    } else if (resonance_type == TensorResonanceType::Pseudoscalar) {
      const auto                &eps3eps4 = decay.scalar;
      const std::complex<double> iD       = decay.propagator;

      auto PseudoscalarProductionAmplitudes = [&](const std::size_t upper_helicity, const std::size_t lower_helicity) {
        return SumTensorResonanceChannels(
            res.TP, tensor_leg_1, tensor_leg_2, upper_helicity, lower_helicity,
            [&](const Tensor2<std::complex<double>, 4, 4> &left, const Tensor2<std::complex<double>, 4, 4> &right,
                const RES_TENSOR_CHANNEL &channel) {
              auto vertex = iG_PPS_total(lts.q1, lts.q2, M0, resonance_type, channel.g_tensor, channel);
              FOR_EACH_4(LI);
              vertex(u, v, k, l) *= phase;
              FOR_EACH_4_END;
              Tensor2<std::complex<double>, 4, 4> central;
              central(alpha1, beta1) = vertex(alpha1, beta1, alpha2, beta2) * right(alpha2, beta2);
              const std::complex<double> block = (-zi) * iD * left(alpha1, beta1) * central(alpha1, beta1);

              std::vector<std::complex<double>> amplitudes;
              amplitudes.reserve(eps3eps4.size());
              for (const auto &decay : eps3eps4) { amplitudes.push_back(block * decay_phase * decay); }
              return amplitudes;
            });
      };

      TensorChannelSum noflip_hamp;
      if (forward_noflip) { noflip_hamp = PseudoscalarProductionAmplitudes(0, 0); }

      // Two helicity states for incoming and outgoing protons
      FOR_PP_HELICITY;

      // Apply proton leg helicity conservation / No helicity flip (high energy
      // limit)
      if (forward_noflip && (ha != h1 || hb != h2)) { continue; }

      if (forward_noflip) {
        AppendPomeronPomeron(resonance_hamp, resonance_bank, noflip_hamp);
      } else {
        const std::size_t      index_1    = spin::BinaryPairHelicityIndex(ha, h1);
        const std::size_t      index_2    = spin::BinaryPairHelicityIndex(hb, h2);
        const TensorChannelSum amplitudes = PseudoscalarProductionAmplitudes(index_1, index_2);
        AppendPomeronPomeron(resonance_hamp, resonance_bank, amplitudes);
      }

      FOR_PP_HELICITY_END;

      // =====================================================================
      // Pomeron-Pomeron-Tensor structure
      //
      //
    } else if (resonance_type == TensorResonanceType::Tensor) {
      // Tensor production structure is contracted directly below
      const auto &iD = decay.tensor;

      // Contract the production side once before the central helicity loop
      auto TensorProductionAmplitudes = [&](const std::size_t upper_helicity, const std::size_t lower_helicity) {
        return SumTensorResonanceChannels(
            res.TP, tensor_leg_1, tensor_leg_2, upper_helicity, lower_helicity,
            [&](const Tensor2<std::complex<double>, 4, 4> &left, const Tensor2<std::complex<double>, 4, 4> &right,
                const RES_TENSOR_CHANNEL &channel) {
              const Tensor2<std::complex<double>, 4, 4> production =
                  iG_PPT_contract(left, right, lts.q1, lts.q2, M0, channel.g_tensor, channel);

              std::vector<std::complex<double>> amplitudes;
              amplitudes.reserve(iD.size());
              for (const auto &decay : iD) {
                amplitudes.push_back((-zi) * phase * production(rho1, rho2) * decay_phase * decay(rho1, rho2));
              }
              return amplitudes;
            });
      };

      // In the high-energy no-flip limit all four surviving proton helicity
      // combinations have the same covariant production block
      TensorChannelSum noflip_hamp;
      if (forward_noflip) { noflip_hamp = TensorProductionAmplitudes(0, 0); }

      // Two helicity states for incoming and outgoing protons
      FOR_PP_HELICITY;

      // Apply proton leg helicity conservation / No helicity flip (high energy
      // limit)
      if (forward_noflip && (ha != h1 || hb != h2)) { continue; }

      if (forward_noflip) {
        AppendPomeronPomeron(resonance_hamp, resonance_bank, noflip_hamp);
      } else {
        const std::size_t      index_1    = spin::BinaryPairHelicityIndex(ha, h1);
        const std::size_t      index_2    = spin::BinaryPairHelicityIndex(hb, h2);
        const TensorChannelSum amplitudes = TensorProductionAmplitudes(index_1, index_2);
        AppendPomeronPomeron(resonance_hamp, resonance_bank, amplitudes);
      }

      FOR_PP_HELICITY_END;

    } else {
      throw std::invalid_argument("MTensorPomeron::ME3: Unknown spin-parity input");
    }
    // ====================================================================

    // Coherent resonances must expose exactly the same forward and decay
    // helicity basis
    if (lts.hamp.empty()) {
      lts.hamp = resonance_hamp;
    } else {
      if (lts.hamp.size() != resonance_hamp.size()) {
        throw AmplitudeFailure(
            "MTensorPomeron::ME3: coherent resonances "
            "have incompatible helicity bases");
      }
      for (const auto &i : indices(lts.hamp)) { lts.hamp[i] += resonance_hamp[i]; }
    }
    pair_bank.Add(resonance_bank);

  }  // Loop over resonances
  // <<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<

  StoreTensorPair(lts, soft_model, tensor_param.exchange, pair_bank, "MTensorPomeron::ME3");

  // Get total amplitude squared 1/4 \sum_h |A_h|^2
  double SumAmp2 = 0.0;
  for (const auto &i : indices(lts.hamp)) { SumAmp2 += gra::math::abs2(lts.hamp[i]); }
  SumAmp2 /= 4;  // Initial state helicity average

  return SumAmp2;  // Amplitude squared
}

// Evaluate the coherent gauge-restored meson continuum and vector resonances
double MTensorPomeron::MEPhoto(gra::LORENTZSCALAR &lts) const {
  return PhotoProduction(lts, true);
}

// Compute the HERA target current with the elastic forward normalization at W0
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::PhotoDissCurrent(
    const ForwardLegState &state, int exchange_pdg, int vector_pdg, double nu, double w2) const {
  const auto &diss = tensor_param.photo_diss.at(vector_pdg);
  auto forward = state;
  forward.final_state = ForwardFinalState::Elastic;
  forward.t = 0.0;
  // Retain the tensor propagator denominator but replace its Regge power by the W0 anchor
  if (!(nu > 0.0) || !std::isfinite(nu)) { throw AmplitudeFailure("Tensor HERA current requires positive finite subenergy"); }
  const auto factor = TensorPropagatorFactor(exchange_pdg, pow2(diss.W0), 0.0) * (pow2(diss.W0) / nu) *
                      flux::PhotoDissFactor(diss, w2, state.t, state.mass2);
  return PomeronPropagatorCurrent(iG_TForwardHE(forward, exchange_pdg), factor);
}

// Derive the diagonal transverse vector-nucleon optical profile from Tensor couplings
// VMD_MIXING enters production only, nuclear propagation omits rho <-> omega conversion
// [REFERENCE: Ewerz, Maniatis and Nachtmann, arXiv:1309.3478, Eqs. (7.19)-(7.21)]
flux::PhotoTargetProfile MTensorPomeron::PhotoProfile(const PARAM_RES &res, double w2) const {
  const auto amplitude = [&](double t) {
    std::complex<double> out = 0.0;
    for (const auto &channel : res.TP.channels) {
      if (channel.active_g_tensor.empty()) { continue; }
      const int pdg = channel.exchange[channel.exchange[0] == PDG::PDG_gamma ? 1 : 0];
      const auto &exchange = tensor_param.exchange.FindExchange(pdg);
      const auto &proton = tensor_param.exchange.FindVertex(pdg, PDG::PDG_p);
      const double beta = tensor_param.exchange.BaryonCoupling(pdg, proton.g_tensor.front());
      const double residue = 3.0 * beta * (2.0 * pow2(res.p.mass) * channel.g_tensor[0] + channel.g_tensor[1]);
      const auto &ff = channel.ff_transfer;
      out += zi * residue * regge::TransferFF(t, proton.ff_transfer) * regge::TransferFF(t, ff) *
             std::pow(-zi * w2 * exchange.ap, exchange.delta + exchange.ap * t);
    }
    return out;
  };
  const double dt = tensor_param.photo_dt;
  const auto forward = amplitude(0.0);
  const auto nearby = amplitude(-dt);
  flux::PhotoTargetProfile out;
  out.sigma_eff = forward.imag() * PDG::GeV2mb;
  out.eta = forward.real() / forward.imag();
  out.slope = std::log(std::norm(forward) / std::norm(nearby)) / dt;
  out.x = pow2(res.p.mass) / w2;
  out.scale2 = pow2(res.p.mass) / 4.0;
  if (!(out.sigma_eff > 0.0) || !(out.slope > 0.0) || !std::isfinite(out.eta) || !std::isfinite(out.slope)) {
    throw AmplitudeFailure("Tensor photoproduction has an invalid vector-nucleon profile");
  }
  return out;
}

// Build the local (d^nu R^{abc}) F_{a nu} P_{bc} effective vertex with g in GeV^-1
MTensorPomeron::VectorKernels MTensorPomeron::Spin3Production(
    const LORENTZSCALAR &lts, const PARAM_RES &res, const TensorDecayState &decay) const {
  VectorKernels out;
  const auto phase = std::polar(1.0, res.TP.phi + (tensor_param.use_zeta ? res.hel_decay.zeta : 0.0));
  for (const auto &channel : res.TP.channels) {
    if (channel.active_ready && channel.active_g_tensor.empty()) { continue; }
    const int exchange = channel.exchange[channel.exchange[0] == 22 ? 1 : 0];
    for (const bool upper : {true, false}) {
      const auto &q = upper ? lts.q1 : lts.q2;
      const auto current = tensor::Spin3Current(q, lts.pfinal[0], decay.spin3);
      const auto factor = zi * qed::e_QED() * channel.g_tensor[0] * phase * decay.propagator *
          regge::MassFF(lts.m2, pow2(res.p.mass), channel.ff_prod) *
          regge::TransferFF((upper ? lts.q2 : lts.q1).M2(), channel.ff_transfer);
      auto [entry, inserted] = (upper ? out.yT : out.Ty).try_emplace(exchange);
      FOR_EACH_3(LI);
      if (inserted) { entry->second(u, v, k) = 0.0; }
      entry->second(u, v, k) += factor * current(u, v, k);
      FOR_EACH_3_END;
    }
  }
  return out;
}

// Evaluate Tensor photoproduction with physical lepton and nuclear beams
double MTensorPomeron::PhotoProduction(LORENTZSCALAR &lts, bool continuum) const {
  lts.proton_good_walker.reset();
  std::vector<nuclear::PhotoTerm> terms;
  for (const auto &[name, res] : lts.process.RESONANCES) {
    TensorDecayState local;
    const auto type = tensor::ClassifyTensorResonance(res.p);
    const auto &decay = TensorDecayForEvent(*this, lts, name, res, type, local);
    const auto kernels = type == TensorResonanceType::Spin3 ? Spin3Production(lts, res, decay)
                                                          : VectorProduction(lts, res, decay.vector);
    for (const auto target_leg : {ForwardBeamLeg::Lower, ForwardBeamLeg::Upper}) {
      const auto target = ResolveForwardLegState(lts, target_leg);
      if (!flux::SupportsPhotoTarget(target)) { continue; }
      const auto elementary = flux::ElementaryPhotoTarget(lts, target);
      const bool upper_gamma = target_leg == ForwardBeamLeg::Lower;
      const auto &q = upper_gamma ? lts.q1 : lts.q2;
      const auto &kernel = upper_gamma ? kernels.yT : kernels.Ty;
      const double nu = (elementary.incoming + elementary.outgoing) * lts.pfinal[0];
      auto profile = target.IsNuclear() && target.upc_model->Param().photo_model != nuclear::PhotoModel::Impulse
          ? PhotoProfile(res, (q + elementary.incoming).M2()) : flux::PhotoTargetProfile{};
      // The neutral isovector exchange changes sign between proton and neutron
      // [REFERENCE: arXiv:1309.3478, Eqs. (3.49)-(3.52)]
      for (int isospin = 0; isospin <= (target.IsNuclear() ? 2 : 0); isospin += 2) {
        if (target.IsNuclear() && std::none_of(kernel.begin(), kernel.end(), [&](const auto& entry) {
              return lts.PDG.FindByPDG(entry.first).isospinX2 == isospin;
            })) { continue; }
        TensorPhotoCurrent current;
        const bool noflip = tensor_param.FORWARD_NOFLIP || target.IsNuclear() || target.IsExcited();
        std::array<MDirac::Spinor, 2> incoming{}, outgoing{};
        if (!noflip) {
          incoming = SpinorStates(elementary.incoming, "u");
          outgoing = SpinorStates(elementary.outgoing, "ubar");
        }
        for (const auto h : spin::BinaryHelicityIndices()) {
          for (const auto prime : spin::BinaryHelicityIndices()) {
            auto &row = current[spin::BinaryPairHelicityIndex(h, prime)];
            for (const auto &mu : LI) { row(mu) = 0.0; }
            if (noflip && h != prime) { continue; }
            if (continuum && (-q.M2() < 0.0 || -q.M2() > tensor_param.photo.q2_max)) { continue; }
            for (const auto &[exchange, central] : kernel) {
              if (target.IsNuclear() && lts.PDG.FindByPDG(exchange).isospinX2 != isospin) { continue; }
              const auto vertex = noflip ? iG_TForwardHE(elementary, exchange) :
                  iG_TForward(elementary, exchange, outgoing[prime], incoming[h], h, prime);
              const auto leg = target.IsExcited() && lts.process.PHOTO_DISSOCIATION == DissociationType::Hera &&
                               tensor_param.exchange.FindExchange(exchange).type == TensorExchangeType::Pomeron
                  ? PhotoDissCurrent(target, exchange, res.p.pdg, nu, (q + target.incoming).M2())
                  : PomeronPropagatorCurrent(vertex, TensorPropagatorFactor(exchange, nu, target.t));
              for (const auto &mu : LI) {
                for (const auto &alpha : LI) {
                  for (const auto &beta : LI) { row(mu) += (-zi) * central(mu, alpha, beta) * leg(alpha, beta); }
                }
              }
            }
          }
        }
        profile.isospin = isospin == 0 ? nuclear::PhotoIsospin::Isoscalar : nuclear::PhotoIsospin::Isovector;
        terms.push_back(TensorPhotoTerm(*this, tensor_param, lts, target, current, profile, continuum));
      }
    }
  }
  if (continuum) {
    MTensorPhoto photo(*this, tensor_param_handle);
    auto additional = photo.Terms(lts);
    terms.insert(terms.end(), std::make_move_iterator(additional.begin()), std::make_move_iterator(additional.end()));
  }
  flux::SumPhotoTerms(lts, std::move(terms));
  return gra::SquaredNorm(lts.hamp) / 4.0;
}

// 2 -> 4 amplitudes
//
// return value: matrix element squared with helicities summed over
// lts.hamp:     individual helicity amplitudes
//
double MTensorPomeron::ME4(gra::LORENTZSCALAR &lts, TensorContinuumMode mode) const {
  lts.proton_good_walker.reset();
  const int spin_left = lts.decaytree[0].p.spinX2;
  const int  pdg_left          = lts.decaytree[0].p.pdg;
  const int  pdg_right         = lts.decaytree[1].p.pdg;
  const int  abs_pdg           = std::abs(pdg_left);
  const bool conjugate_pair    = pdg_left == -pdg_right;
  const bool is_quark          = abs_pdg >= 1 && abs_pdg <= 6;

  std::vector<std::array<int, 2>> ordered_pairs;
  if (mode == TensorContinuumMode::TensorPomeron) {
    for (const auto &pair : tensor_param.exchange.FindContinuumPairs(pdg_left, pdg_right)) {
      const auto ordered = MTensorExchangeModel::OrderedPairs(pair);
      for (const auto &entry : ordered) {
        if (spin_left == 2) {
          if (!tensor_param.exchange.FindActiveTransfers(entry[0], entry[1], pdg_left, pdg_right).empty()) {
            ordered_pairs.push_back(entry);
          }
          continue;
        }
        const auto &upper = tensor_param.exchange.FindVertex(entry[0], pdg_left);
        const auto &lower = tensor_param.exchange.FindVertex(entry[1], pdg_right);
        if (!upper.active_g_tensor.empty() && !lower.active_g_tensor.empty()) { ordered_pairs.push_back(entry); }
      }
    }
  }

  // Tensor continuum uses tensor-local steering while QED keeps generic spin
  // steering
  const bool forward_noflip =
      ForwardNoFlip(mode == TensorContinuumMode::QED ? MTensorPomeronMode::QED : MTensorPomeronMode::Continuum, lts);

  // Kinematics
  const M4Vec pa = lts.pbeam1;
  const M4Vec pb = lts.pbeam2;

  const M4Vec                   p1             = lts.pfinal[1];
  const M4Vec                   p2             = lts.pfinal[2];
  const ForwardLegState         upper_state    = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState         lower_state    = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const TensorForwardSourceBank source_bank    = mode == TensorContinuumMode::QED
                                                     ? TensorForwardSourceBank::PhotonFusion(upper_state, lower_state)
                                                     : TensorForwardSourceBank::Hadronic(upper_state, lower_state);
  std::size_t                   anti_index     = 0;
  std::size_t                   particle_index = 1;
  if ((spin_left == 0 || spin_left == 1) && conjugate_pair && pdg_left > 0) {
    anti_index     = 1;
    particle_index = 0;
  }
  M4Vec p3 = lts.decaytree[anti_index].p4;
  M4Vec p4 = lts.decaytree[particle_index].p4;

  // Intermediate boson/fermion mass
  const double M_ = lts.decaytree[particle_index].p.mass;

  // Momentum convention of sub-diagrams
  //
  // ------<------ anti-particle p3
  // |
  // | \hat{t} (arrow down)
  // |
  // ------>------ particle p4
  //
  // ------<------ particle p4
  // |
  // | \hat{u} (arrow up)
  // |
  // ------>------ anti-particle p3
  //
  const M4Vec pt = pa - p1 - p3;  // => q1 = pt + p3, q2 = pt - p4
  const M4Vec pu = p4 - pa + p1;  // => q2 = pu + p3, q1 = pu - p3

  // ------------------------------------------------------------------

  // Incoming and elastic outgoing proton spinors (2 helicities)
  std::array<MDirac::Spinor, 2> u_a{};
  std::array<MDirac::Spinor, 2> u_b{};
  std::array<MDirac::Spinor, 2> ubar_1{};
  std::array<MDirac::Spinor, 2> ubar_2{};
  const bool                    spinors_required = mode == TensorContinuumMode::QED || !forward_noflip;
  if (spinors_required) {
    u_a = SpinorStates(pa, "u");
    u_b = SpinorStates(pb, "u");
    if (!upper_state.IsExcited()) { ubar_1 = SpinorStates(p1, "ubar"); }
    if (!lower_state.IsExcited()) { ubar_2 = SpinorStates(p2, "ubar"); }
  }

  // ------------------------------------------------------------------

  // Reset
  lts.hamp.clear();
  TensorPairBank pair_bank;

  // Deduce the validated final-state spin basis
  std::string SPINMODE;

  if (lts.decaytree[0].p.spinX2 == 0 && lts.decaytree[1].p.spinX2 == 0) {
    SPINMODE = "2xS";
  } else if (lts.decaytree[0].p.spinX2 == 1 && lts.decaytree[1].p.spinX2 == 1) {
    SPINMODE = "2xF";
  } else if (lts.decaytree[0].p.spinX2 == 2 && lts.decaytree[1].p.spinX2 == 2) {
    SPINMODE = "2xV";
  }

  using ComplexTensor1 = Tensor1<std::complex<double>, 4>;
  using ComplexTensor2 = Tensor2<std::complex<double>, 4, 4>;

  // Build central QED blocks once because they do not depend on proton spin
  std::array<std::array<ComplexTensor2, 2>, 2> qed_t;
  std::array<std::array<ComplexTensor2, 2>, 2> qed_u;
  double                                       qed_factor = 1.0;

  // Build pseudoscalar continuum blocks once per ordered exchange channel
  std::vector<TensorScalarContinuumBlock> scalar_blocks;

  // Build baryon continuum blocks once per event
  MMatrix<std::complex<double>> baryon_iSF_t;
  MMatrix<std::complex<double>> baryon_iSF_u;
  std::array<MDirac::Spinor, 2> baryon_v;
  std::array<MDirac::Spinor, 2> baryon_ubar;

  // Build vector continuum blocks and external polarizations once per channel
  std::vector<TensorVectorBlock> vector_blocks;
  std::array<ComplexTensor1, 3>           vector_eps3;
  std::array<ComplexTensor1, 3>           vector_eps4;

  if (mode == TensorContinuumMode::QED) {
    const MMatrix<std::complex<double>> iSF_t  = iD_F(pt, M_);
    const MMatrix<std::complex<double>> iSF_u  = iD_F(pu, M_);
    const std::array<MDirac::Spinor, 2> v_3    = SpinorStates(p3, "v");
    const std::array<MDirac::Spinor, 2> ubar_4 = SpinorStates(p4, "ubar");
    for (const auto &h3 : indices(v_3)) {
      for (const auto &h4 : indices(ubar_4)) {
        qed_t[h3][h4] = iG_yeebary(ubar_4[h4], iSF_t, v_3[h3]);
        qed_u[h3][h4] = iG_yeebary(ubar_4[h4], iSF_u, v_3[h3]);
      }
    }
    if (is_quark) {
      const double charge = lts.decaytree[particle_index].p.chargeX3 / 3.0;
      qed_factor          = msqrt(math::pow4(charge) * 3.0);
    }
  } else if (SPINMODE == "2xS") {
    for (const auto &pair : ordered_pairs) {
      const auto                &upper_vertex = tensor_param.exchange.FindVertex(pair[0], abs_pdg);
      const auto                &lower_vertex = tensor_param.exchange.FindVertex(pair[1], abs_pdg);
      TensorScalarContinuumBlock block;
      block.exchange = pair;
      if (tensor_param.exchange.FindExchange(pair[0]).rank == 2) {
        block.t_upper = iG_Tpsps(pt, -p3, pair[0], abs_pdg);
        block.u_upper = iG_Tpsps(p4, pu, pair[0], abs_pdg);
      } else {
        // Both crossed vertices follow the same charged-pion line
        block.t_upper_vector = iG_Vpsps(pt, -p3, pair[0], abs_pdg);
        block.u_upper_vector = iG_Vpsps(p4, pu, pair[0], abs_pdg);
      }
      if (tensor_param.exchange.FindExchange(pair[1]).rank == 2) {
        block.t_lower = iG_Tpsps(p4, pt, pair[1], abs_pdg);
        block.u_lower = iG_Tpsps(pu, -p3, pair[1], abs_pdg);
      } else {
        block.t_lower_vector = iG_Vpsps(p4, pt, pair[1], abs_pdg);
        block.u_lower_vector = iG_Vpsps(pu, -p3, pair[1], abs_pdg);
      }
      block.t_scale =
          iD_MES0(pt, M_) * TensorOffshell(upper_vertex, pt.M2(), M_) * TensorOffshell(lower_vertex, pt.M2(), M_);
      block.u_scale =
          iD_MES0(pu, M_) * TensorOffshell(upper_vertex, pu.M2(), M_) * TensorOffshell(lower_vertex, pu.M2(), M_);
      scalar_blocks.push_back(block);
    }
  } else if (SPINMODE == "2xF") {
    baryon_iSF_t = iD_F(pt, M_);
    baryon_iSF_u = iD_F(pu, M_);
    baryon_v     = SpinorStates(p3, "v");
    baryon_ubar  = SpinorStates(p4, "ubar");
  } else if (SPINMODE == "2xV") {
    const int pdg = lts.decaytree[0].p.pdg;
    vector_blocks = TensorVectorBlocks(*this, tensor_param, pdg, p3, p4, pt, pu, lts.pfinal[0].M2());
    vector_eps3  = MassiveSpin1States(p3, "conj", true);
    vector_eps4  = MassiveSpin1States(p4, "conj", true);
  }

  const std::size_t central_states = SPINMODE == "2xS" ? 1U : SPINMODE == "2xF" ? 4U : 9U;
  const std::size_t forward_states = forward_noflip ? 4U : 16U;
  lts.hamp.reserve(forward_states * central_states * source_bank.Size());

  // Two helicity states for incoming and outgoing protons
  FOR_PP_HELICITY;
  // Apply proton leg helicity conservation / No helicity flip
  // (high energy limit)
  // This gives at least 4 x speed improvement
  if (forward_noflip && (ha != h1 || hb != h2)) { continue; }
  if (mode == TensorContinuumMode::TensorPomeron && forward_noflip && (ha != 0 || hb != 0)) { continue; }

  // ==============================================================
  // Lepton or quark pair via two photon fusion
  if (mode == TensorContinuumMode::QED) {
    // Elastic currents and inclusive density eigenvectors include photon
    // transfer
    const auto  photon_sources_1  = iG_yForwardSources(lts, upper_state, ubar_1[h1], u_a[ha], ha, h1);
    const auto  photon_sources_2  = iG_yForwardSources(lts, lower_state, ubar_2[h2], u_b[hb], hb, h2);
    const auto &source_components = source_bank.Components(TensorForwardMechanism::GammaGamma);
    if (source_components.size() != photon_sources_1.size() * photon_sources_2.size()) {
      throw AmplitudeFailure("MTensorPomeron::ME4: gamma-gamma source-bank mismatch");
    }

    std::vector<ComplexTensor1>       right(photon_sources_2.size());
    std::vector<std::complex<double>> source_amplitude(source_bank.Size());
    for (const auto h3 : spin::BinaryHelicityIndices()) {
      for (const auto h4 : spin::BinaryHelicityIndices()) {
        for (const auto &source_2 : indices(photon_sources_2)) {
          for (const auto &nu1 : LI) {
            right[source_2](nu1) = 0.0;
            for (const auto &nu2 : LI) {
              right[source_2](nu1) +=
                  (qed_t[h3][h4](nu2, nu1) + qed_u[h3][h4](nu1, nu2)) * photon_sources_2[source_2](nu2);
            }
          }
        }

        std::fill(source_amplitude.begin(), source_amplitude.end(), 0.0);
        std::size_t source_component = 0;
        for (const auto &source_1 : photon_sources_1) {
          for (const auto &source_2 : indices(photon_sources_2)) {
            std::complex<double> amplitude = 0.0;
            for (const auto &nu1 : LI) { amplitude += source_1(nu1) * right[source_2](nu1); }
            amplitude *= (-zi) * qed_factor;
            source_amplitude[source_components[source_component++]] += amplitude;
          }
        }
        lts.hamp.insert(lts.hamp.end(), source_amplitude.begin(), source_amplitude.end());
      }
    }

    continue;  // skip parts below, only for Pomeron amplitudes
  }

  // ------------------------------------------------------------------
  // ------------------------------------------------------------------

  const auto &pomeron_component = source_bank.Components(TensorForwardMechanism::PomeronPomeron);
  if (pomeron_component.size() != 1) { throw AmplitudeFailure("MTensorPomeron::ME4: invalid PP forward source basis"); }
  // Append one strong amplitude in the common forward source basis
  const auto AppendPomeronPomeron = [&](const std::map<TensorExchangePair, std::complex<double>> &channel) {
    std::complex<double> amplitude = 0.0;
    for (const auto &[exchange, value] : channel) {
      (void)exchange;
      amplitude += value;
    }
    const std::size_t offset = lts.hamp.size();
    lts.hamp.resize(offset + source_bank.Size(), 0.0);
    lts.hamp[offset + pomeron_component[0]] = amplitude;
    std::map<TensorExchangePair, std::vector<std::complex<double>>> block;
    for (const auto &[exchange, value] : channel) {
      auto &source = block[exchange];
      source.resize(source_bank.Size(), 0.0);
      source[pomeron_component[0]] = value;
    }
    pair_bank.Append(block, source_bank.Size());
  };

  // Project each selected exchange once for both continuum topologies
  const auto ProjectLeg = [&](const ForwardLegState &state, const int exchange_pdg, const MDirac::Spinor &outgoing,
                              const MDirac::Spinor &incoming, const std::size_t initial_helicity,
                              const std::size_t final_helicity, const double s, const double t) {
    TensorSoftLeg leg;
    leg.rank = tensor_param.exchange.FindExchange(exchange_pdg).rank;
    if (leg.rank == 2) {
      const auto current = forward_noflip
                               ? iG_TForwardHE(state, exchange_pdg)
                               : iG_TForward(state, exchange_pdg, outgoing, incoming, initial_helicity, final_helicity);
      leg.tensor         = PomeronPropagatorCurrent(current, TensorPropagatorFactor(exchange_pdg, s, t));
    } else {
      const auto current = forward_noflip
                               ? iG_VForwardHE(state, exchange_pdg)
                               : iG_VForward(state, exchange_pdg, outgoing, incoming, initial_helicity, final_helicity);
      leg.vector         = VectorPropagatorCurrent(current, exchange_pdg, s, t);
    }
    return leg;
  };

  std::set<int> active_exchanges;
  for (const auto &pair : ordered_pairs) {
    active_exchanges.insert(pair[0]);
    active_exchanges.insert(pair[1]);
  }
  std::map<int, TensorSoftLeg> upper_t;
  std::map<int, TensorSoftLeg> lower_t;
  std::map<int, TensorSoftLeg> upper_u;
  std::map<int, TensorSoftLeg> lower_u;
  for (const int exchange_pdg : active_exchanges) {
    upper_t.emplace(exchange_pdg,
                    ProjectLeg(upper_state, exchange_pdg, ubar_1[h1], u_a[ha], ha, h1, (p1 + p3).M2(), lts.t1));
    lower_t.emplace(exchange_pdg,
                    ProjectLeg(lower_state, exchange_pdg, ubar_2[h2], u_b[hb], hb, h2, (p2 + p4).M2(), lts.t2));
    upper_u.emplace(exchange_pdg,
                    ProjectLeg(upper_state, exchange_pdg, ubar_1[h1], u_a[ha], ha, h1, (p1 + p4).M2(), lts.t1));
    lower_u.emplace(exchange_pdg,
                    ProjectLeg(lower_state, exchange_pdg, ubar_2[h2], u_b[hb], hb, h2, (p2 + p3).M2(), lts.t2));
  }

  // ==============================================================
  // 2 x pseudoscalar (pion pair, kaon pair ...)
  if (SPINMODE == "2xS") {
    FTensor::Index<'c', 4>                             alpha;
    FTensor::Index<'d', 4>                             beta;
    FTensor::Index<'e', 4>                             mu;
    std::map<TensorExchangePair, std::complex<double>> amplitude;
    for (const auto &block : scalar_blocks) {
      const auto                &ut = upper_t.at(block.exchange[0]);
      const auto                &lt = lower_t.at(block.exchange[1]);
      const auto                &uu = upper_u.at(block.exchange[0]);
      const auto                &lu = lower_u.at(block.exchange[1]);
      const std::complex<double> upper_t_vertex =
          ut.rank == 2 ? ut.tensor(alpha, beta) * block.t_upper(alpha, beta) : ut.vector(mu) * block.t_upper_vector(mu);
      const std::complex<double> lower_t_vertex =
          lt.rank == 2 ? lt.tensor(alpha, beta) * block.t_lower(alpha, beta) : lt.vector(mu) * block.t_lower_vector(mu);
      const std::complex<double> upper_u_vertex =
          uu.rank == 2 ? uu.tensor(alpha, beta) * block.u_upper(alpha, beta) : uu.vector(mu) * block.u_upper_vector(mu);
      const std::complex<double> lower_u_vertex =
          lu.rank == 2 ? lu.tensor(alpha, beta) * block.u_lower(alpha, beta) : lu.vector(mu) * block.u_lower_vector(mu);
      amplitude[block.exchange] +=
          (-zi) * (upper_t_vertex * lower_t_vertex * block.t_scale + upper_u_vertex * lower_u_vertex * block.u_scale);
    }

    // Full amplitude: iM = [ ... ]  <-> M = (-i)*[ ... ]
    AppendPomeronPomeron(amplitude);
  }

  // ==============================================================
  // 2 x fermion (proton-antiproton pair, lambda pair ...)
  else if (SPINMODE == "2xF") {
    // Central fermion helicities
    for (const auto &h3 : indices(baryon_v)) {
      for (const auto &h4 : indices(baryon_ubar)) {
        std::map<TensorExchangePair, std::complex<double>> amplitude;
        for (const auto &pair : ordered_pairs) {
          const auto BaryonCurrent = [&](const TensorSoftLeg &leg, const M4Vec &prime, const M4Vec &momentum,
                                         const int exchange_pdg, const int baryon_pdg) {
            return leg.rank == 2 ? TensorBaryonCurrent(leg.tensor, prime, momentum, exchange_pdg, baryon_pdg)
                                 : VectorBaryonCurrent(leg.vector, prime, momentum, exchange_pdg, baryon_pdg);
          };
          const auto   t_left       = BaryonCurrent(lower_t.at(pair[1]), p4, pt, pair[1], abs_pdg);
          const auto   t_right      = BaryonCurrent(upper_t.at(pair[0]), pt, -p3, pair[0], abs_pdg);
          const auto   u_left       = BaryonCurrent(upper_u.at(pair[0]), p4, pu, pair[0], abs_pdg);
          const auto   u_right      = BaryonCurrent(lower_u.at(pair[1]), pu, -p3, pair[1], abs_pdg);
          const auto   t_propagated = baryon_iSF_t.LeftMultiply(t_left.LeftMultiply(baryon_ubar[h4]));
          const auto   u_propagated = baryon_iSF_u.LeftMultiply(u_left.LeftMultiply(baryon_ubar[h4]));
          const auto  &upper_vertex = tensor_param.exchange.FindVertex(pair[0], abs_pdg);
          const auto  &lower_vertex = tensor_param.exchange.FindVertex(pair[1], abs_pdg);
          const double t_scale = TensorOffshell(upper_vertex, pt.M2(), M_) * TensorOffshell(lower_vertex, pt.M2(), M_);
          const double u_scale = TensorOffshell(upper_vertex, pu.M2(), M_) * TensorOffshell(lower_vertex, pu.M2(), M_);
          amplitude[pair] += (-zi) * (t_right.BilinearForm(t_propagated, baryon_v[h3]) * t_scale +
                                      u_right.BilinearForm(u_propagated, baryon_v[h3]) * u_scale);
        }

        // Full amplitude: iM = [ ... ]  <-> M = (-i)*[ ... ]
        AppendPomeronPomeron(amplitude);
      }
    }
  }

  // ==============================================================
  // 2 x vector meson (rho pair, phi pair ...)
  else if (SPINMODE == "2xV") {
    std::map<TensorExchangePair, ComplexTensor2> amplitude_tensor;
    for (const auto &block : vector_blocks) {
      const ComplexTensor2 direct =
          TensorVectorDirectExchange(upper_t.at(block.exchange[0]).tensor, lower_t.at(block.exchange[1]).tensor,
                                     block.t_upper, block.t_propagator, block.t_lower, block.t_scale);
      const ComplexTensor2 crossed =
          TensorVectorCrossedExchange(upper_u.at(block.exchange[0]).tensor, lower_u.at(block.exchange[1]).tensor,
                                      block.u_upper, block.u_propagator, block.u_lower, block.u_scale);
      auto [entry, inserted] = amplitude_tensor.try_emplace(block.exchange);
      auto &amplitude = entry->second;
      FOR_EACH_2(LI);
      if (inserted) { amplitude(u, v) = 0.0; }
      amplitude(u, v) += direct(u, v) + crossed(u, v);
      FOR_EACH_2_END;
    }

    // Project the reduced covariant tensor on the cached vector helicities
    for (const auto &h3 : indices(vector_eps3)) {
      for (const auto &h4 : indices(vector_eps4)) {
        std::map<TensorExchangePair, std::complex<double>> amplitude;
        for (const auto &[exchange, tensor] : amplitude_tensor) {
          for (const auto &rho3 : LI) {
            for (const auto &rho4 : LI) {
              amplitude[exchange] += (-zi) * vector_eps3[h3](rho3) * tensor(rho3, rho4) * vector_eps4[h4](rho4);
            }
          }
        }
        AppendPomeronPomeron(amplitude);
      }
    }
  }  // 2xV

  FOR_PP_HELICITY_END;

  // The high-energy Tensor currents are identical for all four diagonal rows
  if (mode == TensorContinuumMode::TensorPomeron && forward_noflip) {
    const std::vector<std::complex<double>> diagonal_block = lts.hamp;
    for (std::size_t row = 1; row < 4; ++row) {
      lts.hamp.insert(lts.hamp.end(), diagonal_block.begin(), diagonal_block.end());
    }
    pair_bank.Repeat(3);
  }

  if (mode == TensorContinuumMode::TensorPomeron) {
    StoreTensorPair(lts, soft_model, tensor_param.exchange, pair_bank, "MTensorPomeron::ME4");
  }

  // Get total amplitude squared 1/4 \sum_h |A_h|^2
  double SumAmp2 = 0.0;
  for (const auto &i : indices(lts.hamp)) { SumAmp2 += math::abs2(lts.hamp[i]); }
  SumAmp2 /= 4;  // Initial state helicity average

  return SumAmp2;  // Amplitude squared
}

namespace {

// Build the covariant vector propagator and decay tensor block for V V -> 4PS
Tensor2<std::complex<double>, 4, 4> TensorCascadeDecayBlock(const MTensorPomeron            &tp,
                                                            const std::vector<MDecayBranch> &tree) {
  // Compute each propagated decay in fixed charge order
  const auto decay_current = [&](const MDecayBranch &branch) {
    const bool reversed = branch.legs[0].p.pdg < 0 &&
                          branch.legs[0].p.pdg == -branch.legs[1].p.pdg;
    return tp.VectorDecay(branch.legs[reversed ? 1 : 0].p4, branch.legs[reversed ? 0 : 1].p4,
        branch.p.mass, branch.p.width, branch.p.pdg, branch.hel.g_decay_TP[0], branch.hel.ff_decay);
  };
  const auto left = decay_current(tree[0]);
  const auto right = decay_current(tree[1]);

  Tensor2<std::complex<double>, 4, 4> block;
  for (const auto &rho3 : tp.LI) {
    for (const auto &rho4 : tp.LI) { block(rho3, rho4) = left(rho3) * right(rho4); }
  }
  return block;
}

// Add a scaled helicity vector into the output accumulator
void AddScaledHelicityAmplitudes(std::vector<std::complex<double>> &out, const std::vector<std::complex<double>> &in,
                                 std::complex<double> scale) {
  if (out.empty()) { out.assign(in.size(), 0.0); }
  if (out.size() != in.size()) { throw AmplitudeFailure("MTensorPomeron::ME6: helicity vector dimension mismatch"); }
  for (const auto &i : indices(in)) { out[i] += scale * in[i]; }
}

// Compute coherent scalar-resonance V V cascade decay contraction
std::complex<double> TensorCascadeScalarDecayAmplitude(const MTensorPomeron &tp, gra::LORENTZSCALAR &lts, double M0,
                                                       const std::vector<double> &couplings,
                                                       const regge::FFParam      &ff_decay) {
  FTensor::Index<'a', 4> rho1;
  FTensor::Index<'b', 4> rho2;

  std::complex<double> out   = 0.0;
  const auto           terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  for (const auto &term : terms) {
    const Tensor2<std::complex<double>, 4, 4> iGf0vv =
        tp.iG_f0vv(term.tree[0].p4, term.tree[1].p4, M0, couplings[0], couplings[1], ff_decay);
    const Tensor2<std::complex<double>, 4, 4> decay = TensorCascadeDecayBlock(tp, term.tree);
    out += term.statistics_sign * iGf0vv(rho1, rho2) * decay(rho1, rho2);
  }

  return out;
}

// Compute coherent tensor-resonance V V cascade decay tensor
Tensor2<std::complex<double>, 4, 4> TensorCascadeSpin2DecayTensor(const MTensorPomeron &tp, gra::LORENTZSCALAR &lts,
                                                                  double M0, const std::vector<double> &couplings,
                                                                  const regge::FFParam &ff_decay) {
  FTensor::Index<'a', 4> rho1;
  FTensor::Index<'b', 4> rho2;
  FTensor::Index<'c', 4> alpha;
  FTensor::Index<'d', 4> beta;

  Tensor2<std::complex<double>, 4, 4> out;
  out(alpha, beta) = 0.0;

  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  for (const auto &term : terms) {
    const Tensor4<std::complex<double>, 4, 4, 4, 4> iGf2vv =
        tp.iG_f2vv(term.tree[0].p4, term.tree[1].p4, M0, couplings[0], couplings[1], ff_decay);
    const Tensor2<std::complex<double>, 4, 4> decay = TensorCascadeDecayBlock(tp, term.tree);
    Tensor2<std::complex<double>, 4, 4>       block;
    block(alpha, beta) = iGf2vv(rho1, rho2, alpha, beta) * decay(rho1, rho2);
    out(alpha, beta)   = out(alpha, beta) + term.statistics_sign * block(alpha, beta);
  }

  return out;
}

// Store one cascade amplitude and its ordered exchange decomposition
struct TensorCascadeAmplitude {
  std::vector<std::complex<double>> total;
  TensorPairBank                    channel;
};

// Evaluate raw covariant cascade amplitudes before proposal compensation
TensorCascadeAmplitude TensorCascadeRawHelicityAmplitudes(const MTensorPomeron &tp, const gra::LORENTZSCALAR &lts,
                                                          const std::vector<MDecayBranch> &tree,
                                                          const MTensorPomeronParam       &tensor_param) {
  // Kinematics
  const M4Vec pa = lts.pbeam1;
  const M4Vec pb = lts.pbeam2;

  const M4Vec                   p1                = lts.pfinal[1];
  const M4Vec                   p2                = lts.pfinal[2];
  const ForwardLegState         upper_state       = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState         lower_state       = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);
  const TensorForwardSourceBank source_bank       = TensorForwardSourceBank::Hadronic(upper_state, lower_state);
  const auto                   &pomeron_component = source_bank.Components(TensorForwardMechanism::PomeronPomeron);
  if (pomeron_component.size() != 1) {
    throw AmplitudeFailure("MTensorPomeron::ME6: invalid hadronic forward source basis");
  }

  const bool forward_noflip = tp.ForwardNoFlip(MTensorPomeronMode::Continuum, lts);

  // Build exact proton spinors for the optional forward-flip rows
  std::array<MDirac::Spinor, 2> u_a{};
  std::array<MDirac::Spinor, 2> u_b{};
  std::array<MDirac::Spinor, 2> ubar_1{};
  std::array<MDirac::Spinor, 2> ubar_2{};
  if (!forward_noflip) {
    u_a = tp.SpinorStates(pa, "u");
    u_b = tp.SpinorStates(pb, "u");
    if (!upper_state.IsExcited()) { ubar_1 = tp.SpinorStates(p1, "ubar"); }
    if (!lower_state.IsExcited()) { ubar_2 = tp.SpinorStates(p2, "ubar"); }
  }

  const M4Vec p3 = tree[0].p4;
  const M4Vec p4 = tree[1].p4;

  // Momentum convention of sub-diagrams
  const M4Vec pt = pa - p1 - p3;
  const M4Vec pu = p4 - pa + p1;

  const double s13 = (lts.pfinal[1] + p3).M2();
  const double s24 = (lts.pfinal[2] + p4).M2();
  const double s14 = (lts.pfinal[1] + p4).M2();
  const double s23 = (lts.pfinal[2] + p3).M2();

  const Tensor2<std::complex<double>, 4, 4> decay = TensorCascadeDecayBlock(tp, tree);
  const auto                                transfer_blocks =
      TensorVectorBlocks(tp, tensor_param, tree[0].p.pdg, p3, p4, pt, pu, lts.pfinal[0].M2());
  if (transfer_blocks.empty()) { throw AmplitudeFailure("MTensorPomeron::ME6: no configured transfer path"); }

  std::vector<std::complex<double>> hamp;
  hamp.reserve((forward_noflip ? 4U : 16U) * source_bank.Size());
  TensorPairBank pair_bank;

  // Two helicity states for incoming and outgoing protons
  FOR_PP_HELICITY;

  // Apply proton leg helicity conservation / No helicity flip (high energy
  // limit) This gives at least 4 x speed improvement
  if (forward_noflip && (ha != h1 || hb != h2)) { continue; }
  if (forward_noflip && (ha != 0 || hb != 0)) { continue; }

  std::map<TensorExchangePair, Tensor2<std::complex<double>, 4, 4>> exchange_tensor;
  for (const auto &block : transfer_blocks) {
    const auto upper_current = forward_noflip
                                   ? tp.iG_TForwardHE(upper_state, block.exchange[0])
                                   : tp.iG_TForward(upper_state, block.exchange[0], ubar_1[h1], u_a[ha], ha, h1);
    const auto lower_current = forward_noflip
                                   ? tp.iG_TForwardHE(lower_state, block.exchange[1])
                                   : tp.iG_TForward(lower_state, block.exchange[1], ubar_2[h2], u_b[hb], hb, h2);
    const auto upper_t =
        tp.PomeronPropagatorCurrent(upper_current, tp.TensorPropagatorFactor(block.exchange[0], s13, lts.t1));
    const auto lower_t =
        tp.PomeronPropagatorCurrent(lower_current, tp.TensorPropagatorFactor(block.exchange[1], s24, lts.t2));
    const auto upper_u =
        tp.PomeronPropagatorCurrent(upper_current, tp.TensorPropagatorFactor(block.exchange[0], s14, lts.t1));
    const auto lower_u =
        tp.PomeronPropagatorCurrent(lower_current, tp.TensorPropagatorFactor(block.exchange[1], s23, lts.t2));
    const auto direct =
        TensorVectorDirectExchange(upper_t, lower_t, block.t_upper, block.t_propagator, block.t_lower, block.t_scale);
    const auto crossed =
        TensorVectorCrossedExchange(upper_u, lower_u, block.u_upper, block.u_propagator, block.u_lower, block.u_scale);
    auto [entry, inserted] = exchange_tensor.try_emplace(block.exchange);
    for (const auto &rho3 : tp.LI) {
      for (const auto &rho4 : tp.LI) {
        if (inserted) { entry->second(rho3, rho4) = 0.0; }
        entry->second(rho3, rho4) += direct(rho3, rho4) + crossed(rho3, rho4);
      }
    }
  }

  std::map<TensorExchangePair, std::complex<double>> exchange_amplitude;
  for (const auto &[exchange, tensor] : exchange_tensor) {
    for (const auto &rho3 : tp.LI) {
      for (const auto &rho4 : tp.LI) { exchange_amplitude[exchange] += (-zi) * tensor(rho3, rho4) * decay(rho3, rho4); }
    }
  }
  std::complex<double>                                            amp = 0.0;
  std::map<TensorExchangePair, std::vector<std::complex<double>>> block;
  for (const auto &[exchange, amplitude] : exchange_amplitude) {
    amp += amplitude;
    auto &source = block[exchange];
    source.resize(source_bank.Size(), 0.0);
    source[pomeron_component[0]] = amplitude;
  }
  const std::size_t offset = hamp.size();
  hamp.resize(offset + source_bank.Size(), 0.0);
  hamp[offset + pomeron_component[0]] = amp;
  pair_bank.Append(block, source_bank.Size());
  FOR_PP_HELICITY_END;

  // The high-energy Tensor currents are identical for all four diagonal rows
  if (forward_noflip) {
    const std::vector<std::complex<double>> diagonal_block = hamp;
    for (std::size_t row = 1; row < 4; ++row) { hamp.insert(hamp.end(), diagonal_block.begin(), diagonal_block.end()); }
    pair_bank.Repeat(3);
  }

  return {std::move(hamp), std::move(pair_bank)};
}

// Evaluate tensor cascade amplitudes with coherent stable-leaf assignments
TensorCascadeAmplitude TensorCascadeHelicityAmplitudes(const MTensorPomeron &tp, gra::LORENTZSCALAR &lts,
                                                       const MTensorPomeronParam &tensor_param) {
  const std::vector<gra::spin::StableLeafAmplitudeTerm> terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  if (terms.size() > 1) {
    TensorCascadeAmplitude out;
    for (const auto &term : terms) {
      const auto raw = TensorCascadeRawHelicityAmplitudes(tp, lts, term.tree, tensor_param);
      AddScaledHelicityAmplitudes(out.total, raw.total, term.statistics_sign);
      out.channel.Add(raw.channel, term.statistics_sign);
    }
    return out;
  }

  return TensorCascadeRawHelicityAmplitudes(tp, lts, lts.decaytree, tensor_param);
}

}  // namespace

// 2 -> 6 amplitudes through sequential vector-pair cascade proposal
//
// return value: matrix element squared with helicities summed over
// lts.hamp:     raw complete helicity amplitudes, with proposal compensation in
// MProcess
//
double MTensorPomeron::ME6(gra::LORENTZSCALAR &lts) const {
  lts.proton_good_walker.reset();
  if (!lts.process.SPINGEN || !lts.process.SPINDEC) {
    throw std::invalid_argument(
        "MTensorPomeron::ME6: SPINGEN and SPINDEC must both be enabled for "
        "covariant cascade amplitudes");
  }
  TensorCascadeAmplitude amplitude = TensorCascadeHelicityAmplitudes(*this, lts, tensor_param);
  lts.hamp                         = std::move(amplitude.total);
  StoreTensorPair(lts, soft_model, tensor_param.exchange, amplitude.channel, "MTensorPomeron::ME6");

  // Get total amplitude squared 1/4 \sum_h |A_h|^2
  double SumAmp2 = 0.0;
  for (const auto &i : indices(lts.hamp)) { SumAmp2 += math::abs2(lts.hamp[i]); }
  SumAmp2 /= 4;  // Initial state helicity average

  return SumAmp2;  // Amplitude squared
}

// External massless spin-1 polarization sum for outgoing states
//
// Input tensor M_{\mu \nu \kappa \rho} (indices down)
//
std::vector<Tensor2<std::complex<double>, 4, 4>> MTensorPomeron::MasslessSpin1PolSum(
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3, const M4Vec &p4) const {
  std::vector<Tensor2<std::complex<double>, 4, 4>> hamp;
  Tensor2<std::complex<double>, 4, 4>              temp;

  // Massless polarization vectors (2 helicities)
  const bool                                            INDEX_UP   = true;
  const std::array<Tensor1<std::complex<double>, 4>, 2> eps_3_conj = MasslessSpin1States(p3, "conj", INDEX_UP);
  const std::array<Tensor1<std::complex<double>, 4>, 2> eps_4_conj = MasslessSpin1States(p4, "conj", INDEX_UP);

  // Loop over states
  for (const auto &h3 : indices(eps_3_conj)) {
    Tensor3<std::complex<double>, 4, 4, 4> first;
    for (const auto &rho4 : LI) {
      for (const auto &mu3 : LI) {
        for (const auto &mu4 : LI) {
          first(rho4, mu3, mu4) = 0.0;
          for (const auto &rho3 : LI) { first(rho4, mu3, mu4) += eps_3_conj[h3](rho3) * M(rho3, rho4, mu3, mu4); }
        }
      }
    }
    for (const auto &h4 : indices(eps_4_conj)) {
      for (const auto &mu3 : LI) {
        for (const auto &mu4 : LI) {
          temp(mu3, mu4) = 0.0;
          for (const auto &rho4 : LI) { temp(mu3, mu4) += eps_4_conj[h4](rho4) * first(rho4, mu3, mu4); }
        }
      }
      hamp.push_back(temp);
    }
  }
  return hamp;
}

// External massless spin-1 polarization sum for outgoing states
//
// Input tensor M_{\mu \nu} (indices down)
//
std::vector<std::complex<double>> MTensorPomeron::MasslessSpin1PolSum(const Tensor2<std::complex<double>, 4, 4> &M,
                                                                      const M4Vec &p3, const M4Vec &p4) const {
  std::vector<std::complex<double>> hamp;

  // Massless polarization vectors (2 helicities)
  const bool                                            INDEX_UP   = true;
  const std::array<Tensor1<std::complex<double>, 4>, 2> eps_3_conj = MasslessSpin1States(p3, "conj", INDEX_UP);
  const std::array<Tensor1<std::complex<double>, 4>, 2> eps_4_conj = MasslessSpin1States(p4, "conj", INDEX_UP);

  // Loop over massive Spin-1 helicity states
  for (const auto &h3 : indices(eps_3_conj)) {
    Tensor1<std::complex<double>, 4> first;
    for (const auto &rho4 : LI) {
      first(rho4) = 0.0;
      for (const auto &rho3 : LI) { first(rho4) += eps_3_conj[h3](rho3) * M(rho3, rho4); }
    }
    for (const auto &h4 : indices(eps_4_conj)) {
      std::complex<double> amp = 0.0;
      for (const auto &rho4 : LI) { amp += eps_4_conj[h4](rho4) * first(rho4); }
      hamp.push_back(amp);
    }
  }
  return hamp;
}

// External massive spin-1 (vector) polarization sum for outgoing states
//
// Input tensor M_{\mu \nu \kappa \rho} (indices down)
//
std::vector<Tensor2<std::complex<double>, 4, 4>> MTensorPomeron::MassiveSpin1PolSum(
    const Tensor4<std::complex<double>, 4, 4, 4, 4> &M, const M4Vec &p3, const M4Vec &p4) const {
  std::vector<Tensor2<std::complex<double>, 4, 4>> hamp;
  Tensor2<std::complex<double>, 4, 4>              temp;

  // Massive polarization vectors (3 helicities)
  const bool                                            INDEX_UP   = true;
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_3_conj = MassiveSpin1States(p3, "conj", INDEX_UP);
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_4_conj = MassiveSpin1States(p4, "conj", INDEX_UP);

  // Loop over states
  for (const auto &h3 : indices(eps_3_conj)) {
    Tensor3<std::complex<double>, 4, 4, 4> first;
    for (const auto &rho4 : LI) {
      for (const auto &mu3 : LI) {
        for (const auto &mu4 : LI) {
          first(rho4, mu3, mu4) = 0.0;
          for (const auto &rho3 : LI) { first(rho4, mu3, mu4) += eps_3_conj[h3](rho3) * M(rho3, rho4, mu3, mu4); }
        }
      }
    }
    for (const auto &h4 : indices(eps_4_conj)) {
      for (const auto &mu3 : LI) {
        for (const auto &mu4 : LI) {
          temp(mu3, mu4) = 0.0;
          for (const auto &rho4 : LI) { temp(mu3, mu4) += eps_4_conj[h4](rho4) * first(rho4, mu3, mu4); }
        }
      }
      hamp.push_back(temp);
    }
  }
  return hamp;
}

// External massive spin-1 (vector) polarization sum for outgoing states
//
// Input tensor M_{\mu \nu \kappa \rho} (indices down)
//
std::vector<std::complex<double>> MTensorPomeron::MassiveSpin1PolSum(const Tensor2<std::complex<double>, 4, 4> &M,
                                                                     const M4Vec &p3, const M4Vec &p4) const {
  std::vector<std::complex<double>> hamp;

  // Massive polarization vectors (3 helicities)
  const bool                                            INDEX_UP   = true;
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_3_conj = MassiveSpin1States(p3, "conj", INDEX_UP);
  const std::array<Tensor1<std::complex<double>, 4>, 3> eps_4_conj = MassiveSpin1States(p4, "conj", INDEX_UP);

  // Loop over massive Spin-1 helicity states
  for (const auto &h3 : indices(eps_3_conj)) {
    Tensor1<std::complex<double>, 4> first;
    for (const auto &rho4 : LI) {
      first(rho4) = 0.0;
      for (const auto &rho3 : LI) { first(rho4) += eps_3_conj[h3](rho3) * M(rho3, rho4); }
    }
    for (const auto &h4 : indices(eps_4_conj)) {
      std::complex<double> amp = 0.0;
      for (const auto &rho4 : LI) { amp += eps_4_conj[h4](rho4) * first(rho4); }
      hamp.push_back(amp);
    }
  }

  /*
  // Polarization sum (completeness relation applied)
  // already simplified expression => second terms vanishes and leaves only
  g[u][k] * g[v][l] double Amp2 = 0.0;

  FOR_EACH_4(LI);
  const std::complex<double> contraction = std::conj(M(u, v)) * M(k, l) *
  g[u][k] * g[v][l]; Amp2 += std::real(contraction);     // real for casting to
  double FOR_EACH_4_END;

  hamp.push_back(msqrt(Amp2)); // Phase information lost here
  */

  return hamp;
}

// 2xvector -> 2xpseudoscalar decay block
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_vv2psps(const std::vector<M4Vec> &p, int PDG,
                                                               double decay_coupling, const regge::FFParam &ff_decay) const {
  const auto &vector = tensor_param.FindVector(PDG);
  const auto left = VectorDecay(p[0], p[1], vector.mass, vector.width, PDG, decay_coupling, ff_decay);
  const auto right = VectorDecay(p[2], p[3], vector.mass, vector.width, PDG, decay_coupling, ff_decay);

  Tensor2<std::complex<double>, 4, 4> BLOCK;
  for (const auto &rho3 : LI) {
    for (const auto &rho4 : LI) { BLOCK(rho3, rho4) = left(rho3) * right(rho4); }
  }

  return BLOCK;
}

// -------------------------------------------------------------------------
// Vertex functions

// Gamma-Lepton-Lepton vertex
// contracted in \bar{spinor} G_\mu \bar{spinor}
//
// iGamma_\mu (p',p)
//
// High Energy Limit:
// ubar(p',\lambda') \gamma_\mu u(p,\lambda ~= (p' + p)_\mu
// \delta_{\lambda',\lambda}
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_yee(const M4Vec &prime, const M4Vec &p, const MDirac::Spinor &ubar,
                                                        const MDirac::Spinor &u) const {
  // const double q2 = (prime-p).M2();
  const double e = msqrt(qed::alpha_QED() * 4.0 * PI);  // ~ 0.3, no running

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) {
    // \bar{spinor} [Gamma Matrix] \spinor product
    T(mu) = zi * e * gamma_lo[mu].BilinearForm(ubar, u);
  }
  return T;
}

// Gamma-Proton-Proton vertex iGamma_\mu (p', p)
// contracted in \bar{spinor} iG_\mu \bar{spinor}
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_ypp(const M4Vec &prime, const M4Vec &p, const MDirac::Spinor &ubar,
                                                        const MDirac::Spinor &u) const {
  const double t    = (prime - p).M2();
  const double e    = msqrt(qed::alpha_QED() * 4.0 * PI);  // ~ 0.3, no running
  const M4Vec  psum = prime - p;
  const double F1_  = form::F1(t, model_tune->Structure());
  const double F2_  = form::F2(t, model_tune->Structure());

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) {
    // \bar{spinor} [Gamma Matrix] \spinor product
    const std::complex<double> dirac = gamma_lo[mu].BilinearForm(ubar, u);

    std::complex<double> pauli = 0.0;
    for (const auto &nu : LI) { pauli += psum[nu] * sigma_lo[mu][nu].BilinearForm(ubar, u); }

    T(mu) = (-zi * e) * (F1_ * dirac + zi / (2 * PDG::mp) * F2_ * pauli);
  }

  return T;
}

// gamma - electron - fermion propagator - positron - gamma vertex
// iGamma_{\mu \nu}
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_yeebary(const MDirac::Spinor                &ubar,
                                                               const MMatrix<std::complex<double>> &iSF,
                                                               const MDirac::Spinor                &v) const {
  const double               e      = msqrt(qed::alpha_QED() * 4.0 * PI);  // ~ 0.3, no running
  const std::complex<double> vertex = zi * e;

  Tensor2<std::complex<double>, 4, 4> T;
  for (const auto &mu : LI) {
    const MDirac::Spinor lhs = (gamma_lo[mu] * iSF).LeftMultiply(ubar);

    for (const auto &nu : LI) { T(mu, nu) = vertex * gamma_lo[nu].BilinearForm(lhs, v) * vertex; }
  }
  return T;
}

// High-Energy limit proton-Pomeron-proton vertex times \delta^{lambda_prime,
// \lambda} (helicity conservation)
// (~x 2 faster evaluation than the exact spinor structure)
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_PppHE(const M4Vec &prime, const M4Vec p) const {
  const M4Vec                psum   = prime + p;
  const auto                &vertex = tensor_param.exchange.FindVertex(995, PDG::PDG_p);
  const std::complex<double> FACTOR =
      -zi * 3.0 * vertex.g_tensor.front() * regge::TransferFF((prime - p).M2(), vertex.ff_transfer);

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR * (psum % u) * (psum % v);
  FOR_EACH_2_END;

  return T;
}

// Compute the elastic or inclusive high-energy Tensor Pomeron beam current
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_PForwardHE(const ForwardLegState &state) const {
  return iG_TForwardHE(state, 995);
}

// Compute the high-energy beam current for one configured rank-two exchange
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.43), (3.49)]
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_TForwardHE(const ForwardLegState &state,
                                                                  const int              exchange_pdg) const {
  const M4Vec  psum     = state.outgoing + state.incoming;
  const auto  &vertex   = tensor_param.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  const double coupling = tensor_param.exchange.BaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  double       profile  = regge::TransferFF(state.t, vertex.ff_transfer);
  if (state.IsExcited()) {
    const SoftExchangeId exchange = tensor_param.exchange.SoftId(exchange_pdg, *soft_model);
    if (soft_model->Exchange(exchange).Role() != SoftExchangeRole::Pomeron) {
      Tensor2<std::complex<double>, 4, 4> out;
      FOR_EACH_2(LI);
      out(u, v) = 0.0;
      FOR_EACH_2_END;
      return out;
    }
    profile = soft_model->ForwardExcitationFactor(exchange, state.t, state.mass2);
  }
  if (!std::isfinite(profile) || profile < 0.0) {
    throw AmplitudeFailure("MTensorPomeron::iG_PForwardHE: invalid forward profile");
  }

  // The inclusive ansatz retains the elastic high-energy tensor shape
  const std::complex<double>          factor = -zi * 3.0 * coupling * profile;
  Tensor2<std::complex<double>, 4, 4> out;
  FOR_EACH_2(LI);
  out(u, v) = factor * (psum % u) * (psum % v);
  FOR_EACH_2_END;
  return out;
}

// Compute one elastic exact or inclusive identity Tensor Pomeron beam current
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_PForward(const ForwardLegState &state,
                                                                const MDirac::Spinor &ubar, const MDirac::Spinor &u,
                                                                const std::size_t initial_helicity,
                                                                const std::size_t final_helicity) const {
  return iG_TForward(state, 995, ubar, u, initial_helicity, final_helicity);
}

// Compute one exact elastic or inclusive current for a rank-two exchange
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.43), (3.49)]
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_TForward(const ForwardLegState &state, const int exchange_pdg,
                                                                const MDirac::Spinor &ubar, const MDirac::Spinor &u,
                                                                const std::size_t initial_helicity,
                                                                const std::size_t final_helicity) const {
  if (state.IsExcited()) {
    if (initial_helicity != final_helicity) {
      Tensor2<std::complex<double>, 4, 4> out;
      FOR_EACH_2(LI);
      out(u, v) = 0.0;
      FOR_EACH_2_END;
      return out;
    }
    return iG_TForwardHE(state, exchange_pdg);
  }
  const auto  &vertex   = tensor_param.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  const double coupling = tensor_param.exchange.BaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  const M4Vec  psum     = state.outgoing + state.incoming;
  const MMatrix<std::complex<double>> slash  = FSlash(psum);
  const std::complex<double>          factor = -zi * 3.0 * coupling * regge::TransferFF(state.t, vertex.ff_transfer);
  Tensor2<std::complex<double>, 4, 4> out;
  for (const auto &mu : LI) {
    for (const auto &nu : LI) {
      MMatrix<std::complex<double>> matrix = (gamma_lo[mu] * (psum % nu) + gamma_lo[nu] * (psum % mu)) * 0.5;
      if (mu == nu) { matrix = matrix - slash * 0.25 * g[mu][nu]; }
      out(mu, nu) = factor * matrix.BilinearForm(ubar, u);
    }
  }
  const std::complex<double> phase = TensorForwardColliderCurrentPhase(state, initial_helicity, final_helicity);
  FOR_EACH_2(LI);
  out(u, v) = phase * out(u, v);
  FOR_EACH_2_END;
  return out;
}

// Compute the high-energy beam current for one configured vector exchange
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.59), (3.68)]
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_VForwardHE(const ForwardLegState &state,
                                                               const int              exchange_pdg) const {
  const auto  &vertex   = tensor_param.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  const double coupling = tensor_param.exchange.VectorBaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  double       profile  = regge::TransferFF(state.t, vertex.ff_transfer);
  if (state.IsExcited()) {
    const SoftExchangeId exchange = tensor_param.exchange.SoftId(exchange_pdg, *soft_model);
    if (soft_model->Exchange(exchange).Role() != SoftExchangeRole::Pomeron) {
      Tensor1<std::complex<double>, 4> out;
      for (const auto &mu : LI) { out(mu) = 0.0; }
      return out;
    }
    profile = soft_model->ForwardExcitationFactor(exchange, state.t, state.mass2);
  }
  if (!std::isfinite(profile) || profile < 0.0) {
    throw AmplitudeFailure("MTensorPomeron::iG_VForwardHE: invalid forward profile");
  }
  const M4Vec                      psum = state.outgoing + state.incoming;
  Tensor1<std::complex<double>, 4> out;
  for (const auto &mu : LI) { out(mu) = -zi * coupling * profile * (psum % mu); }
  return out;
}

// Compute one exact elastic or inclusive current for a vector exchange
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.59), (3.68)]
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_VForward(const ForwardLegState &state, const int exchange_pdg,
                                                             const MDirac::Spinor &ubar, const MDirac::Spinor &u,
                                                             const std::size_t initial_helicity,
                                                             const std::size_t final_helicity) const {
  if (state.IsExcited()) {
    if (initial_helicity != final_helicity) {
      Tensor1<std::complex<double>, 4> out;
      for (const auto &mu : LI) { out(mu) = 0.0; }
      return out;
    }
    return iG_VForwardHE(state, exchange_pdg);
  }
  const auto  &vertex               = tensor_param.exchange.FindVertex(exchange_pdg, PDG::PDG_p);
  const double coupling             = tensor_param.exchange.VectorBaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  const std::complex<double> factor = -zi * coupling * regge::TransferFF(state.t, vertex.ff_transfer) *
                                      TensorForwardColliderCurrentPhase(state, initial_helicity, final_helicity);
  Tensor1<std::complex<double>, 4> out;
  for (const auto &mu : LI) { out(mu) = factor * gamma_lo[mu].BilinearForm(ubar, u); }
  return out;
}

// Compute orthogonal elastic or inclusive EPA sources after photon transfer
std::vector<Tensor1<std::complex<double>, 4>> MTensorPomeron::iG_yForwardSources(
    const LORENTZSCALAR &lts, const ForwardLegState &state, const MDirac::Spinor &ubar, const MDirac::Spinor &u,
    const std::size_t initial_helicity, const std::size_t final_helicity, const TensorPhotonVertex vertex) const {
  if (state.IsNuclear() || state.IsExcited()) {
    const auto sectors = flux::ForwardPhotonFluxSectors(state);
    std::vector<Tensor1<std::complex<double>, 4>> out(2 * sectors.size());
    // Define transverse EPA states in the beam CM and boost their currents to the event frame
    const M4Vec beam = lts.pbeam1 + lts.pbeam2;
    auto photon = state.transfer;
    kinematics::LorentzBoost(beam, beam.M(), photon, -1);
    auto polarization = MasslessSpin1States(photon, "none", true);
    for (auto &eps : polarization) {
      const auto current = MDirac::BoostCurrent({eps(0), -eps(1), -eps(2), -eps(3)}, beam, beam.M(), 1);
      for (const auto &mu : LI) { eps(mu) = (mu == 0 ? 1.0 : -1.0) * current[mu]; }
    }
    const double charge = qed::Charge(lts, state.Index()) < 0.0 ? -1.0 : 1.0;
    for (const auto &i : indices(sectors)) {
      for (const auto &mu : LI) {
        out[i](mu) = 0.0;
        out[sectors.size() + i](mu) = 0.0;
        if (initial_helicity != final_helicity || !(state.xi > 0.0)) { continue; }
        const double norm = charge / std::sqrt(2.0 * state.xi);
        out[i](mu) = norm * sectors[i].source.parallel * (polarization[0](mu) - polarization[1](mu));
        out[sectors.size() + i](mu) = norm * sectors[i].source.perpendicular * (polarization[0](mu) + polarization[1](mu));
      }
    }
    return out;
  }
  {
    Tensor1<std::complex<double>, 4> out;
    for (const auto &nu : LI) { out(nu) = 0.0; }
    Tensor1<std::complex<double>, 4> current;
    const int                        abs_pdg = std::abs(state.emitter.pdg);
    if (abs_pdg == 11 || abs_pdg == 13 || abs_pdg == 15) {
      const auto emitter = qed::Emitter(state.emitter, state.t, "MTensorPomeron::iG_yForwardSources lepton emitter");
      const auto elastic = DiracPauliCurrent(state.incoming, state.outgoing, -state.transfer, emitter.fermion_kind,
                                             TensorFermionHelicity(initial_helicity),
                                             TensorFermionHelicity(final_helicity), emitter.mass, 1.0, 0.0);
      for (const auto &mu : LI) { current(mu) = math::zi * qed::e_QED() * elastic[mu]; }
    } else if (abs_pdg == PDG::PDG_p && vertex == TensorPhotonVertex::Dirac) {
      const double factor = -qed::e_QED() * form::F1(state.t, model_tune->Structure());
      for (const auto &mu : LI) { current(mu) = zi * factor * gamma_lo[mu].BilinearForm(ubar, u); }
    } else if (abs_pdg == PDG::PDG_p) {
      current = iG_ypp(state.outgoing, state.incoming, ubar, u);
    } else {
      throw AmplitudeFailure("MTensorPomeron::iG_yForwardSources: unsupported emitter");
    }
    const Tensor2<std::complex<double>, 4, 4> propagator = iD_y(state.t);
    for (const auto &nu : LI) {
      for (const auto &mu : LI) { out(nu) += current(mu) * propagator(mu, nu); }
    }
    const std::complex<double> phase = TensorForwardColliderCurrentPhase(state, initial_helicity, final_helicity);
    for (const auto &nu : LI) { out(nu) *= phase; }
    return {out};
  }

}

// Pomeron-Proton-Proton [-Neutron-Neutron, -Antiproton-Antiproton] vertex
// contracted in \bar{spinor} G_{\mu\nu} \bar{spinor}
//
// i\Gamma_{\mu\nu} (p', p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_Ppp(const M4Vec &prime, const M4Vec &p,
                                                           const MDirac::Spinor &ubar, const MDirac::Spinor &u) const {
  const double t      = (prime - p).M2();
  const M4Vec  psum   = prime + p;
  const auto  &vertex = tensor_param.exchange.FindVertex(995, PDG::PDG_p);

  Tensor2<std::complex<double>, 4, 4> T;
  const std::complex<double>          FACTOR = -zi * 3.0 * vertex.g_tensor.front() * regge::TransferFF(t, vertex.ff_transfer);

  // Feynman slash
  const MMatrix<std::complex<double>> slash = FSlash(psum);

  // Aux matrix
  MMatrix<std::complex<double>> A;

  for (const auto &mu : LI) {
    for (const auto &nu : LI) {
      A = (gamma_lo[mu] * (psum % nu) + gamma_lo[nu] * (psum % mu)) * 0.5;
      if (mu == nu) { A = A - slash * 0.25 * g[mu][nu]; }

      // \bar{spinor} [Gamma Matrix] \spinor product
      T(mu, nu) = FACTOR * A.BilinearForm(ubar, u);
    }
  }
  return T;
}

// Pomeron x outgoing proton spinor ubar x
//           <antiproton/proton fermion propagator>
//         x outgoing antiproton spinor v x Pomeron vertex function
//
// iG_{\mu_2\nu_2\mu_1\nu_1}
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PppbarP(const M4Vec &prime, const MDirac::Spinor &ubar,
                                                                     const M4Vec                         &pt,
                                                                     const MMatrix<std::complex<double>> &iSF,
                                                                     const MDirac::Spinor &v, const M4Vec &p,
                                                                     double gPBB, const regge::FFParam &ff_transfer) const {
  const std::complex<double> lhs_FACTOR = -zi * 3.0 * gPBB * regge::TransferFF((prime - pt).M2(), ff_transfer);
  const std::complex<double> rhs_FACTOR = -zi * 3.0 * gPBB * regge::TransferFF((pt - p).M2(), ff_transfer);

  const M4Vec lhs_psum = prime + pt;
  const M4Vec rhs_psum = pt + p;

  // Feynman slashes
  const MMatrix<std::complex<double>> lhs_slash = FSlash(lhs_psum);
  const MMatrix<std::complex<double>> rhs_slash = FSlash(rhs_psum);

  // Aux matrices
  MMatrix<std::complex<double>> lhs_A;
  MMatrix<std::complex<double>> rhs_A;

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  for (const auto &mu2 : LI) {
    for (const auto &nu2 : LI) {
      lhs_A = (gamma_lo[mu2] * (lhs_psum % nu2) + gamma_lo[nu2] * (lhs_psum % mu2)) * 0.5;
      if (mu2 == nu2) { lhs_A = lhs_A - lhs_slash * 0.25 * g[mu2][nu2]; }

      // Fermion propagator applied in the middle
      const MDirac::Spinor ubarM = (lhs_A * iSF).LeftMultiply(ubar);

      for (const auto &mu1 : LI) {
        for (const auto &nu1 : LI) {
          rhs_A = (gamma_lo[mu1] * (rhs_psum % nu1) + gamma_lo[nu1] * (rhs_psum % mu1)) * 0.5;
          if (mu1 == nu1) { rhs_A = rhs_A - rhs_slash * 0.25 * g[mu1][nu1]; }

          // \bar{spinor} [Gamma Matrix] \spinor product
          T(mu2, nu2, mu1, nu1) = lhs_FACTOR * rhs_A.BilinearForm(ubarM, v) * rhs_FACTOR;
        }
      }
    }
  }
  return T;
}

// ======================================================================

// Compute whether one production coupling is large enough to evaluate
bool MTensorPomeron::ActiveProductionCoupling(double coupling) const {
  return std::abs(coupling) > model_tune->Global().coupling_min;
}

// Compute the common factorized PP-resonance form factor
double MTensorPomeron::PPResonanceFormFactor(const M4Vec &q1, const M4Vec &q2, double M0,
                                             const RES_TENSOR_CHANNEL &channel) const {
  const auto &transfer = channel.ff_transfer;
  const regge::FFParam &prod =
      channel.ff_prod;
  return regge::TransferFF(q1.M2(), transfer) * regge::TransferFF(q2.M2(), transfer) *
         regge::MassFF((q1 + q2).M2(), pow2(M0), prod);
}

// Pomeron (\mu\nu) - Pomeron (\kappa\lambda) - Scalar/Pseudoscalar resonance
// vertex function: i\Gamma_{\mu\nu,\kappa\lambda,\rho\sigma}
//
// Input as contravariant (upper index) 4-vector, M0 the resonance peak mass
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                                       TensorResonanceType        type,
                                                                       const std::vector<double> &g_PPS,
                                                                       const RES_TENSOR_CHANNEL &channel) const {


  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;
  FOR_EACH_4_END;

  if (type == TensorResonanceType::Scalar) {
    if (ActiveProductionCoupling(g_PPS[0])) {
      const auto iG0 = iG_PPS_0();
      FOR_EACH_4(LI);
      T(u, v, k, l) += g_PPS[0] * iG0(u, v, k, l);
      FOR_EACH_4_END;
    }
    if (ActiveProductionCoupling(g_PPS[1])) {
      const auto iG1 = iG_PPS_1(q1, q2, g_PPS[1]);
      FOR_EACH_4(LI);
      T(u, v, k, l) += iG1(u, v, k, l);
      FOR_EACH_4_END;
    }
  } else if (type == TensorResonanceType::Pseudoscalar) {
    if (ActiveProductionCoupling(g_PPS[0])) {
      const auto iG0 = iG_PPPS_0(q1, q2, g_PPS[0]);
      FOR_EACH_4(LI);
      T(u, v, k, l) += iG0(u, v, k, l);
      FOR_EACH_4_END;
    }
    if (ActiveProductionCoupling(g_PPS[1])) {
      const auto iG1 = iG_PPPS_1(q1, q2, g_PPS[1]);
      FOR_EACH_4(LI);
      T(u, v, k, l) += iG1(u, v, k, l);
      FOR_EACH_4_END;
    }
  } else {
    throw std::invalid_argument("MTensorPomeron::iG_PPS_total: type must be scalar or pseudoscalar");
  }

  const double form_factor = PPResonanceFormFactor(q1, q2, M0, channel);
  FOR_EACH_4(LI);
  T(u, v, k, l) *= form_factor;
  FOR_EACH_4_END;
  return T;
}

// Pomeron-Pomeron-Scalar coupling structure #0
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (0,0)
//
// Input as contravariant (upper index) 4-vectors
//
// Apply coupling g_PPS[0] outside this!
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_0() const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = zi * S0;

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Scalar coupling structure #1
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (2,2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPS_1(const M4Vec &q1, const M4Vec &q2,
                                                                   double g_PPS) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const double               q1q2   = q1 * q2;
  const std::complex<double> FACTOR = zi * g_PPS / (2 * S0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) =
      FACTOR * ((q1 % k) * (q2 % u) * g[v][l] + (q1 % k) * (q2 % v) * g[u][l] + (q1 % l) * (q2 % u) * g[v][k] +
                (q1 % l) * (q2 % v) * g[u][k] - 2.0 * q1q2 * (g[u][k] * g[v][l] + g[v][k] * g[u][l]));
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Pseudoscalar coupling structure #0
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (1,1)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPPS_0(const M4Vec &q1, const M4Vec &q2,
                                                                    double g_PPPS) const {
  const double S0 = 1.0;  // Mass scale (GeV)

  const M4Vec q1_q2 = q1 - q2;
  const M4Vec q     = q1 + q2;

  const std::complex<double>                FACTOR = zi * g_PPPS / (2.0 * S0);
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;

  // Contract r and s indices
  for (const auto &r : LI) {
    for (const auto &s : LI) {
      // Note +=
      T(u, v, k, l) += FACTOR *
                       (g[u][k] * eps_lo(v, l, r, s) + g[v][k] * eps_lo(u, l, r, s) + g[u][l] * eps_lo(v, k, r, s) +
                        g[v][l] * eps_lo(u, k, r, s)) *
                       q1_q2[r] * q[s];
    }
  }
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pomeron-Pseudoscalar coupling structure #1
// iG_{\mu \nu \kappa \lambda} ~ (l,s) = (3,3)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_PPPS_1(const M4Vec &q1, const M4Vec &q2,
                                                                    double g_PPPS) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const M4Vec                q1_q2  = q1 - q2;
  const M4Vec                q      = q1 + q2;
  const double               q1q2   = q1 * q2;
  const std::complex<double> FACTOR = zi * g_PPPS / gra::math::pow3(S0);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;

  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;

  // Contract r and s indices
  for (const auto &r : LI) {
    for (const auto &s : LI) {
      // Note +=
      T(u, v, k, l) += FACTOR *
                       (eps_lo(v, l, r, s) * ((q1 % k) * (q2 % u) - q1q2 * g[u][k]) +
                        eps_lo(u, l, r, s) * ((q1 % k) * (q2 % v) - q1q2 * g[v][k]) +
                        eps_lo(v, k, r, s) * ((q1 % l) * (q2 % u) - q1q2 * g[u][l]) +
                        eps_lo(u, k, r, s) * ((q1 % l) * (q2 % v) - q1q2 * g[v][l])) *
                       q1_q2[r] * q[s];
    }
  }
  FOR_EACH_4_END;

  return T;
}

// ======================================================================

// Compute the unsymmetrized rank-8 tensor entering the axial (2,2) coupling
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020), Eq. (A2)]
double MTensorPomeron::GammaPPA22Base(std::size_t kappa, std::size_t lambda, std::size_t rho, std::size_t sigma,
                                      std::size_t mu, std::size_t nu, std::size_t alpha, std::size_t beta) const {
  return g[kappa][rho] * g[mu][sigma] * eps_lo(lambda, nu, alpha, beta) +
         g[lambda][rho] * g[mu][sigma] * eps_lo(kappa, nu, alpha, beta) +
         g[kappa][sigma] * g[mu][rho] * eps_lo(lambda, nu, alpha, beta) +
         g[lambda][sigma] * g[mu][rho] * eps_lo(kappa, nu, alpha, beta) +
         g[kappa][rho] * g[mu][lambda] * eps_lo(sigma, nu, alpha, beta) +
         g[sigma][kappa] * g[mu][lambda] * eps_lo(rho, nu, alpha, beta) +
         g[rho][lambda] * g[mu][kappa] * eps_lo(sigma, nu, alpha, beta) +
         g[sigma][lambda] * g[mu][kappa] * eps_lo(rho, nu, alpha, beta) -
         g[kappa][lambda] * g[mu][rho] * eps_lo(sigma, nu, alpha, beta) -
         g[kappa][lambda] * g[mu][sigma] * eps_lo(rho, nu, alpha, beta) -
         g[kappa][mu] * g[rho][sigma] * eps_lo(lambda, nu, alpha, beta) -
         g[lambda][mu] * g[rho][sigma] * eps_lo(kappa, nu, alpha, beta);
}

// Construct the tensor-Pomeron axial-vector (l,S)=(2,2) vertex
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020), Eq. (2.11)]
MTensor<std::complex<double>> MTensorPomeron::iG_PPA_22(const M4Vec &q1, const M4Vec &q2, double g_PPA) const {
  const double                  S0     = 1.0;
  const double                  factor = -g_PPA / (8.0 * pow2(S0));
  const M4Vec                   d      = q1 - q2;
  const M4Vec                   p      = q1 + q2;
  MTensor<std::complex<double>> out({4, 4, 4, 4, 4}, 0.0);

  for (const auto kappa : LI) {
    for (const auto lambda : LI) {
      for (const auto rho : LI) {
        for (const auto sigma : LI) {
          for (const auto alpha : LI) {
            double contraction = 0.0;
            for (const auto mu : LI) {
              for (const auto nu : LI) {
                for (const auto beta : LI) {
                  contraction += d[mu] * d[nu] * p[beta] * PPA_22({kappa, lambda, rho, sigma, mu, nu, alpha, beta});
                }
              }
            }
            out(kappa, lambda, rho, sigma, alpha) = factor * contraction;
          }
        }
      }
    }
  }
  return out;
}

// Construct the tensor-Pomeron axial-vector (l,S)=(4,4) vertex
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020), Eq. (2.12)]
MTensor<std::complex<double>> MTensorPomeron::iG_PPA_44(const M4Vec &q1, const M4Vec &q2, double g_PPA) const {
  const double                  S0     = 1.0;
  const double                  factor = g_PPA / (4.0 * pow2(pow2(S0)));
  const M4Vec                   d      = q1 - q2;
  const M4Vec                   p      = q1 + q2;
  const double                  d2     = d.M2();
  MTensor<std::complex<double>> out({4, 4, 4, 4, 4}, 0.0);

  for (const auto kappa : LI) {
    for (const auto lambda : LI) {
      const double A_kl = (d % kappa) * (d % lambda) - 0.25 * g[kappa][lambda] * d2;
      for (const auto rho : LI) {
        for (const auto sigma : LI) {
          const double A_rs = (d % rho) * (d % sigma) - 0.25 * g[rho][sigma] * d2;
          for (const auto alpha : LI) {
            double contraction = 0.0;
            for (const auto beta : LI) {
              double B_rs = 0.0;
              double B_kl = 0.0;
              for (const auto mu : LI) {
                B_rs +=
                    d[mu] * ((d % rho) * eps_lo(sigma, mu, alpha, beta) + (d % sigma) * eps_lo(rho, mu, alpha, beta));
                B_kl += d[mu] *
                        ((d % kappa) * eps_lo(lambda, mu, alpha, beta) + (d % lambda) * eps_lo(kappa, mu, alpha, beta));
              }
              contraction += p[beta] * (A_kl * B_rs + A_rs * B_kl);
            }
            out(kappa, lambda, rho, sigma, alpha) = factor * contraction;
          }
        }
      }
    }
  }
  return out;
}

// Sum the two axial-vector tensor-Pomeron couplings and apply form factors
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 102, 114003 (2020)]
MTensor<std::complex<double>> MTensorPomeron::iG_PPA_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                           const RES_TENSOR_CHANNEL &channel) const {
  const auto &g_PPA = channel.g_tensor;

  const double                  form_factor = PPResonanceFormFactor(q1, q2, M0, channel);
  MTensor<std::complex<double>> out({4, 4, 4, 4, 4}, 0.0);

  if (ActiveProductionCoupling(g_PPA[0])) {
    const MTensor<std::complex<double>> structure = iG_PPA_22(q1, q2, g_PPA[0]);
    FOR_EACH_5(LI);
    out(u, v, k, l, r) += form_factor * structure(u, v, k, l, r);
    FOR_EACH_5_END;
  }
  if (ActiveProductionCoupling(g_PPA[1])) {
    const MTensor<std::complex<double>> structure = iG_PPA_44(q1, q2, g_PPA[1]);
    FOR_EACH_5(LI);
    out(u, v, k, l, r) += form_factor * structure(u, v, k, l, r);
    FOR_EACH_5_END;
  }

  return out;
}

// ======================================================================

// Pomeron (\mu\nu) - Pomeron (\kappa\lambda) - Tensor resonance (\rho\sigma)
// vertex function: i\Gamma_{\mu\nu,\kappa\lambda,\rho\sigma}
//
// Input as contravariant (upper index) 4-vector, M0 the resonance peak mass
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_total(const M4Vec &q1, const M4Vec &q2, double M0,
                                                           const std::vector<double> &g_PPT,
                                                           const RES_TENSOR_CHANNEL &channel) const {

  const double form_factor = PPResonanceFormFactor(q1, q2, M0, channel);

  // Construct total vertex by summing up the tensors, init with zeros!
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));

  // The (l,S) = (0,2) tensor is completely kinematics independent
  if (ActiveProductionCoupling(g_PPT[0])) {
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += g_PPT[0] * PPT_00(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[1])) {
    const MTensor<std::complex<double>> structure = iG_PPT_12(q1, q2, g_PPT[1], 1);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[2])) {
    const MTensor<std::complex<double>> structure = iG_PPT_12(q1, q2, g_PPT[2], 2);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[3])) {
    const MTensor<std::complex<double>> structure = iG_PPT_03(q1, q2, g_PPT[3]);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[4])) {
    const MTensor<std::complex<double>> structure = iG_PPT_04(q1, q2, g_PPT[4]);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[5])) {
    const MTensor<std::complex<double>> structure = iG_PPT_05(q1, q2, g_PPT[5]);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }
  if (ActiveProductionCoupling(g_PPT[6])) {
    const MTensor<std::complex<double>> structure = iG_PPT_06(q1, q2, g_PPT[6]);
    FOR_EACH_6(LI);
    T(u, v, k, l, r, s) += structure(u, v, k, l, r, s);
    FOR_EACH_6_END;
  }

  FOR_EACH_6(LI);
  T(u, v, k, l, r, s) *= form_factor;
  FOR_EACH_6_END;

  return T;
}

// Contract the PP-tensor-resonance vertex directly with both Pomeron currents
//
FTensor::Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_PPT_contract(
    const FTensor::Tensor2<std::complex<double>, 4, 4> &left, const FTensor::Tensor2<std::complex<double>, 4, 4> &right,
    const M4Vec &q1, const M4Vec &q2, double M0, const std::vector<double> &g_PPT,
    const RES_TENSOR_CHANNEL &channel) const {


  const double        form_factor = PPResonanceFormFactor(q1, q2, M0, channel);
  const double        S0          = 1.0;  // Mass scale (GeV)
  const double        q1q2        = q1 * q2;
  std::array<bool, 7> active{};
  for (const auto &i : indices(active)) { active[i] = ActiveProductionCoupling(g_PPT[i]); }

  const bool need_left_R       = active[0] || active[1] || active[2];
  const bool need_right_R      = active[0];
  const bool need_right_RU     = active[1] || active[2];
  const bool need_left_q2      = active[1] || active[2] || active[3] || active[4] || active[5];
  const bool need_right_q1     = active[3] || active[5];
  const bool need_right_q1U    = active[1] || active[2] || active[4];
  const bool need_left_right_R = active[1] || active[2] || active[4];
  const bool need_q1q2_DDUU    = active[1] || active[2];
  const bool need_q1q2_UUDD    = active[4];
  const bool need_q1q2_DDDD    = active[6];
  const bool need_q2q2_q1q1    = active[5] || active[6];

  Tensor1<double, 4> q1_U = {q1[0], q1[1], q1[2], q1[3]};
  Tensor1<double, 4> q2_U = {q2[0], q2[1], q2[2], q2[3]};
  Tensor1<double, 4> q1_D = {q1 % 0, q1 % 1, q1 % 2, q1 % 3};
  Tensor1<double, 4> q2_D = {q2 % 0, q2 % 1, q2 % 2, q2 % 3};

  // Contractions common to several tensor structures
  Tensor2<std::complex<double>, 4, 4> left_R;
  Tensor2<std::complex<double>, 4, 4> right_R;
  Tensor2<std::complex<double>, 4, 4> right_RU;

  FOR_EACH_2(LI);
  if (need_left_R) { left_R(u, v) = 0.0; }
  if (need_right_R) { right_R(u, v) = 0.0; }
  if (need_right_RU) { right_RU(u, v) = 0.0; }
  for (auto const &a : LI) {
    for (auto const &b : LI) {
      if (need_left_R) { left_R(u, v) += left(a, b) * R_DDDD(a, b, u, v); }
      if (need_right_R) { right_R(u, v) += right(a, b) * R_DDDD(a, b, u, v); }
      if (need_right_RU) { right_RU(u, v) += right(a, b) * R_DDDU(a, b, u, v); }
    }
  }
  FOR_EACH_2_END;

  Tensor1<std::complex<double>, 4> left_q2;
  Tensor1<std::complex<double>, 4> right_q1;
  Tensor1<std::complex<double>, 4> right_q1U;

  for (auto const &a1 : LI) {
    if (need_left_q2) { left_q2(a1) = 0.0; }
    if (need_right_q1) { right_q1(a1) = 0.0; }
    if (need_right_q1U) { right_q1U(a1) = 0.0; }
    for (auto const &u1 : LI) {
      for (auto const &a : LI) {
        for (auto const &b : LI) {
          if (need_left_q2) { left_q2(a1) += left(a, b) * q2_U(u1) * R_DDDD(a, b, u1, a1); }
          if (need_right_q1) { right_q1(a1) += right(a, b) * q1_U(u1) * R_DDDD(a, b, u1, a1); }
          if (need_right_q1U) { right_q1U(a1) += right(a, b) * q1_U(u1) * R_DDDU(a, b, u1, a1); }
        }
      }
    }
  }

  std::complex<double> left_right_R = 0.0;
  std::complex<double> left_q2q2    = 0.0;
  std::complex<double> right_q1q1   = 0.0;
  for (auto const &a : LI) {
    for (auto const &b : LI) {
      if (need_left_right_R) {
        for (auto const &k : LI) {
          for (auto const &l : LI) { left_right_R += left(a, b) * right(k, l) * R_DDDD(a, b, k, l); }
        }
      }
      if (need_q2q2_q1q1) {
        for (auto const &u1 : LI) {
          for (auto const &v1 : LI) {
            left_q2q2 += left(a, b) * q2_U(u1) * q2_U(v1) * R_DDDD(a, b, u1, v1);
            right_q1q1 += right(a, b) * q1_U(u1) * q1_U(v1) * R_DDDD(a, b, u1, v1);
          }
        }
      }
    }
  }

  Tensor2<double, 4, 4> q1q2_DDUU;
  Tensor2<double, 4, 4> q1q2_UUDD;
  Tensor2<double, 4, 4> q1q2_DDDD;
  FOR_EACH_2(LI);
  if (need_q1q2_DDUU) { q1q2_DDUU(u, v) = 0.0; }
  if (need_q1q2_UUDD) { q1q2_UUDD(u, v) = 0.0; }
  if (need_q1q2_DDDD) { q1q2_DDDD(u, v) = 0.0; }
  for (auto const &a : LI) {
    for (auto const &b : LI) {
      if (need_q1q2_DDUU) { q1q2_DDUU(u, v) += q1_D(a) * q2_D(b) * R_DDUU(u, v, a, b); }
      if (need_q1q2_UUDD) { q1q2_UUDD(u, v) += q1_D(a) * q2_D(b) * R_UUDD(a, b, u, v); }
      if (need_q1q2_DDDD) { q1q2_DDDD(u, v) += q1_U(a) * q2_U(b) * R_DDDD(u, v, a, b); }
    }
  }
  FOR_EACH_2_END;

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = 0.0;
  FOR_EACH_2_END;

  // Coupling structure #0, (l,s) = (0,2)
  if (active[0]) {
    const std::complex<double> FACTOR = 2.0 * zi * S0 * g_PPT[0];
    FOR_EACH_2(LI);
    std::complex<double> contraction = 0.0;
    for (auto const &v1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &s1 : LI) {
          contraction += left_R(s1, v1) * right_R(v1, l1) * R_DDDD(u, v, l1, s1) * g[v1][v1] * g[l1][l1] * g[s1][s1];
        }
      }
    }
    T(u, v) += FACTOR * contraction;
    FOR_EACH_2_END;
  }

  // Coupling structures #1 and #2, (l,s) = (2,0) -/+ (2,2)
  auto AddPPT12 = [&](const std::size_t index, double sign) {
    if (!active[index]) { return; }
    const std::complex<double> FACTOR = -2.0 * zi / S0 * g_PPT[index];

    FOR_EACH_2(LI);
    std::complex<double> T1_ = 0.0;
    std::complex<double> T2_ = 0.0;
    std::complex<double> T3_ = 0.0;
    for (auto const &a1 : LI) {
      for (auto const &r1 : LI) {
        for (auto const &s1 : LI) {
          T1_ += left_R(r1, a1) * right_RU(s1, a1) * R_DDUU(u, v, r1, s1);
          T2_ += left_q2(a1) * right_RU(s1, a1) * q1_D(r1) * R_DDUU(u, v, r1, s1);
          T3_ += left_R(r1, a1) * right_q1U(a1) * q2_D(s1) * R_DDUU(u, v, r1, s1);
        }
      }
    }
    T(u, v) += FACTOR * (q1q2 * T1_ + sign * T2_ + sign * T3_ + left_right_R * q1q2_DDUU(u, v));
    FOR_EACH_2_END;
  };

  AddPPT12(1, -1.0);
  AddPPT12(2, 1.0);

  // Coupling structure #3, (l,s) = (2,4)
  if (active[3]) {
    const std::complex<double> FACTOR = -zi / S0 * g_PPT[3];
    FOR_EACH_2(LI);
    std::complex<double> contraction = 0.0;
    for (auto const &v1 : LI) {
      for (auto const &l1 : LI) {
        contraction += (left_q2(v1) * right_q1(l1) + left_q2(l1) * right_q1(v1)) * R_UUDD(v1, l1, u, v);
      }
    }
    T(u, v) += FACTOR * contraction;
    FOR_EACH_2_END;
  }

  // Coupling structure #4, (l,s) = (4,2)
  if (active[4]) {
    const std::complex<double> FACTOR      = -2.0 * zi / math::pow3(S0) * g_PPT[4];
    std::complex<double>       contraction = -2.0 * q1q2 * left_right_R;
    for (auto const &a1 : LI) {
      contraction += left_q2(a1) * right_q1U(a1);
      contraction += left_q2(a1) * right_q1U(a1);
    }
    FOR_EACH_2(LI);
    T(u, v) += FACTOR * contraction * q1q2_UUDD(u, v);
    FOR_EACH_2_END;
  }

  // Coupling structure #5, (l,s) = (4,4)
  if (active[5]) {
    const std::complex<double> FACTOR = zi / math::pow3(S0) * g_PPT[5];
    FOR_EACH_2(LI);
    std::complex<double> A = 0.0;
    std::complex<double> C = 0.0;
    for (auto const &v1 : LI) {
      for (auto const &r1 : LI) {
        A += left_q2(v1) * q2_D(r1) * R_UUDD(v1, r1, u, v);
        C += right_q1(v1) * q1_D(r1) * R_UUDD(v1, r1, u, v);
      }
    }
    T(u, v) += FACTOR * (A * right_q1q1 + C * left_q2q2);
    FOR_EACH_2_END;
  }

  // Coupling structure #6, (l,s) = (6,4)
  if (active[6]) {
    const std::complex<double> FACTOR = -2.0 * zi / math::pow5(S0) * g_PPT[6];
    FOR_EACH_2(LI);
    T(u, v) += FACTOR * left_q2q2 * right_q1q1 * q1q2_DDDD(u, v);
    FOR_EACH_2_END;
  }

  FOR_EACH_2(LI);
  T(u, v) *= form_factor;
  FOR_EACH_2_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structures #1 and #2
// iG_{\mu \nu \kappa \lambda \rho \sigma}
// ~ (l,s) = (2,0) - (2,2)
// ~ (l,s) = (2,0) + (2,2)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_12(const M4Vec &q1, const M4Vec &q2, double g_PPT,
                                                        int mode) const {
  if (!(mode == 1 || mode == 2)) {
    throw std::invalid_argument("MTensorPomeron::iG_PPT_12: Error, mode should be 1 or 2");
  }
  const double               S0     = 1.0;  // Mass scale (GeV)
  const double               q1q2   = q1 * q2;
  const std::complex<double> FACTOR = -2.0 * zi / S0 * g_PPT;

  // Coupling structure 1 or 2
  const double sign = (mode == 1) ? -1.0 : 1.0;

  // Init with zeros!
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);

  // Pre-calculated tensor contraction structures T1,T2,T2 used here
  T(u, v, k, l, r, s) += q1q2 * T1(u, v, k, l, r, s);

  for (auto const &u1 : LI) {
    for (auto const &r1 : LI) {
      // Note +=
      T(u, v, k, l, r, s) += sign * q2[u1] * (q1 % r1) * T2({u, v, k, l, r, s, u1, r1});
    }
  }
  for (auto const &u1 : LI) {
    for (auto const &s1 : LI) {
      // Note +=
      T(u, v, k, l, r, s) += sign * q1[u1] * (q2 % s1) * T3({u, v, k, l, r, s, u1, s1});
    }
  }
  for (auto const &r1 : LI) {
    for (auto const &s1 : LI) {
      // Note +=
      T(u, v, k, l, r, s) += (q1 % r1) * (q2 % s1) * R_DDDD(u, v, k, l) * R_DDUU(r, s, r1, s1);
    }
  }

  T(u, v, k, l, r, s) *= FACTOR;
  FOR_EACH_6_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #3
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (2,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_03(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = -zi / S0 * g_PPT;

  Tensor3<double, 4, 4, 4> A;
  Tensor3<double, 4, 4, 4> B;
  Tensor3<double, 4, 4, 4> C;
  Tensor3<double, 4, 4, 4> D;

  // Manual contractions
  FOR_EACH_3(LI);
  A(u, v, k) = 0.0;
  B(u, v, k) = 0.0;
  for (auto const &x : LI) {
    // Note +=
    A(u, v, k) += q2[x] * R_DDDD(u, v, x, k);
    B(u, v, k) += q1[x] * R_DDDD(u, v, x, k);
  }
  FOR_EACH_3_END;

  // Final rank-6 tensor
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  for (auto const &v1 : LI) {
    for (auto const &l1 : LI) {
      // Note +=
      T(u, v, k, l, r, s) += (A(u, v, v1) * B(k, l, l1) + A(u, v, l1) * B(k, l, v1)) * R_UUDD(v1, l1, r, s);
    }
  }
  T(u, v, k, l, r, s) *= FACTOR;
  /*

  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
          for (auto const &u1 : LI) {

            // Note +=
            T(ind) += FACTOR * (q2[u1] * R_DDDD(u,v,u1,v1) * q1[a1] *
  R_DDDD(k,l,a1,l1) + q2[a1] * R_DDDD(u,v,a1,l1) * q1[u1] * R_DDDD(k,l,u1,v1)) *
  R_UUDD(v1,l1,r,s);
        }
      }
    }
  }
  */
  FOR_EACH_6_END;

  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #4
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (4,2)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_04(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = -2.0 * zi / math::pow3(S0) * g_PPT;
  const double               q1q2   = q1 * q2;

  Tensor1<double, 4> q1_D = {q1 % 0, q1 % 1, q1 % 2, q1 % 3};
  Tensor1<double, 4> q2_D = {q2 % 0, q2 % 1, q2 % 2, q2 % 3};

  Tensor4<double, 4, 4, 4, 4> A;
  Tensor4<double, 4, 4, 4, 4> B;
  Tensor2<double, 4, 4>       C;

  // Manual contractions (not possible with autocontraction)
  FOR_EACH_4(LI);
  A(u, v, k, l) = 0.0;
  B(u, v, k, l) = 0.0;

  for (auto const &u1 : LI) {
    for (auto const &v1 : LI) {
      for (auto const &a1 : LI) {
        // Note +=
        A(u, v, k, l) += q2[v1] * R_DDDD(u, v, v1, a1) * q1[u1] * R_DDDU(k, l, u1, a1);
        B(u, v, k, l) += q2[u1] * R_DDDD(u, v, u1, a1) * q1[v1] * R_DDDU(k, l, v1, a1);
      }
    }
  }
  FOR_EACH_4_END;

  // Autocontractions
  {
    FTensor::Index<'a', 4> r;
    FTensor::Index<'b', 4> s;
    FTensor::Index<'c', 4> a1;
    FTensor::Index<'d', 4> l1;

    C(r, s) = q1_D(a1) * q2_D(l1) * R_UUDD(a1, l1, r, s);
  }

  // Final rank-6 expression
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  T(u, v, k, l, r, s) = FACTOR * (A(u, v, k, l) + B(u, v, k, l) - 2.0 * q1q2 * R_DDDD(u, v, k, l)) * C(r, s);

  /*
  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &u1 : LI) {

          // Note +=
          T(ind) += FACTOR * (q2[v1] * R_DDDD(u,v,v1,a1) * q1[u1] *
  R_DDDU(k,l,u1,a1)
                            + q2[u1] * R_DDDD(u,v,u1,a1) * q1[v1] *
  R_DDDU(k,l,v1,a1)
                            - 2.0 * q1q2 * R_DDDD(u,v,k,l)) * (q1 % a1) * (q2 %
  l1) * R_UUDD(a1,l1,r,s);
        }
      }
    }
  }
  */
  FOR_EACH_6_END;
  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #05
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (4,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_05(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = zi / math::pow3(S0) * g_PPT;

  Tensor1<double, 4> q1_U = {q1[0], q1[1], q1[2], q1[3]};
  Tensor1<double, 4> q2_U = {q2[0], q2[1], q2[2], q2[3]};
  Tensor1<double, 4> q1_D = {q1 % 0, q1 % 1, q1 % 2, q1 % 3};
  Tensor1<double, 4> q2_D = {q2 % 0, q2 % 1, q2 % 2, q2 % 3};

  Tensor4<double, 4, 4, 4, 4> A;
  Tensor2<double, 4, 4>       B;
  Tensor4<double, 4, 4, 4, 4> C;
  Tensor2<double, 4, 4>       D;

  // Manual contractions (not possible with autocontraction)
  FOR_EACH_4(LI);
  A(u, v, k, l) = 0.0;
  C(u, v, k, l) = 0.0;

  for (auto const &u1 : LI) {
    for (auto const &v1 : LI) {
      for (auto const &r1 : LI) {
        // Note +=
        A(u, v, k, l) += q2_U(u1) * R_DDDD(u, v, u1, v1) * q2_D(r1) * R_UUDD(v1, r1, k, l);
        C(u, v, k, l) += q1_U(u1) * R_DDDD(u, v, u1, v1) * q1_D(r1) * R_UUDD(v1, r1, k, l);
      }
    }
  }
  FOR_EACH_4_END;

  // Autocontractions
  {
    FTensor::Index<'a', 4> u;
    FTensor::Index<'b', 4> v;
    FTensor::Index<'c', 4> k;
    FTensor::Index<'d', 4> l;
    FTensor::Index<'j', 4> a1;
    FTensor::Index<'k', 4> l1;

    B(k, l) = q1_U(a1) * q1_U(l1) * R_DDDD(k, l, a1, l1);
    D(u, v) = q2_U(a1) * q2_U(l1) * R_DDDD(u, v, a1, l1);
  }

  // Final rank-6 expression
  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  T({u, v, k, l, r, s}) = FACTOR * (A(u, v, r, s) * B(k, l) + C(k, l, r, s) * D(u, v));

  /*
  // NAIVE LOOP
  for (auto const &v1 : LI) {
    for (auto const &a1 : LI) {
      for (auto const &l1 : LI) {
        for (auto const &u1 : LI) {
          for (auto const &r1 : LI) {

            // Note +=
            T(ind) += FACTOR * (q2[u1] * (q2 % r1) * R_DDDD(u,v,u1,v1) *
  R_UUDD(v1,r1,r,s) * q1[a1]
  * q1[l1] * R_DDDD(k,l,a1,l1)
                             +  q1[u1] * (q1 % r1) * R_DDDD(k,l,u1,v1) *
  R_UUDD(v1,r1,r,s) * q2[a1]
  * q2[l1] * R_DDDD(u,v,a1,l1));
          }
        }
      }
    }
  }
  */
  FOR_EACH_6_END;
  return T;
}

// Pomeron-Pomeron-Tensor coupling structure #6
// iG_{\mu \nu \kappa \lambda \rho \sigma} ~ (l,s) = (6,4)
//
MTensor<std::complex<double>> MTensorPomeron::iG_PPT_06(const M4Vec &q1, const M4Vec &q2, double g_PPT) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const std::complex<double> FACTOR = -2.0 * zi / math::pow5(S0) * g_PPT;

  FTensor::Index<'a', 4> a;
  FTensor::Index<'b', 4> b;
  FTensor::Index<'c', 4> c;
  FTensor::Index<'d', 4> d;

  Tensor1<double, 4> q1_U = {q1[0], q1[1], q1[2], q1[3]};
  Tensor1<double, 4> q2_U = {q2[0], q2[1], q2[2], q2[3]};

  Tensor2<double, 4, 4> A;
  Tensor2<double, 4, 4> B;
  Tensor2<double, 4, 4> C;

  // Autocontractions
  A(a, b) = q2_U(c) * q2_U(d) * R_DDDD(a, b, c, d);
  B(a, b) = q1_U(c) * q1_U(d) * R_DDDD(a, b, c, d);
  C(a, b) = q1_U(c) * q2_U(d) * R_DDDD(a, b, c, d);

  MTensor<std::complex<double>> T = MTensor({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));
  FOR_EACH_6(LI);
  T(u, v, k, l, r, s) = FACTOR * A(u, v) * B(k, l) * C(r, s);
  FOR_EACH_6_END;

  return T;
}

// ======================================================================

// f0 (Scalar) - (Pseudo)Scalar - (Pseudo)Scalar vertex function
// i\Gamma
//
// Input as contravariant (upper index) 4-vectors, meson on-shell mass M0
//
// Relations: p3_\mu iG^{\mu \nu} = 0, p4_\nu iG^{\mu \nu} = 0
//
std::complex<double> MTensorPomeron::iG_f0ss(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                             const regge::FFParam &ff_decay) const {
  const double          S0   = 1.0;  // Mass scale (GeV)
  const regge::FFParam &decay = ff_decay;
  return zi * g1 * S0 * regge::MassFF((p3 + p4).M2(), pow2(M0), decay);
}

// f0 (Scalar) - Vector (Massive) - Vector (Massive) vertex function
// i\Gamma_{\mu \nu}
//
// Input as contravariant (upper index) 4-vectors, meson on-shell mass M0
//
// Relations: p3_\mu iG^{\mu \nu} = 0, p4_\nu iG^{\mu \nu} = 0
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_f0vv(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                                            double g2, const regge::FFParam &ff_decay) const {
  const double S0 = 1.0;  // Mass scale (GeV)

  const double p3p4 = p3 * p4;
  const double p3sq = p3.M2();
  const double p4sq = p4.M2();

  const regge::FFParam      &decay   = ff_decay;
  const double               F_      = regge::MassFF((p3 + p4).M2(), pow2(M0), decay);
  const std::complex<double> FACTOR1 = zi * g1 * 2.0 / math::pow3(S0) * F_;
  const std::complex<double> FACTOR2 = zi * g2 * 2.0 / S0 * F_;

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR1 * (p3sq * p4sq * g[u][v] - p4sq * (p3 % u) * (p3 % v) - p3sq * (p4 % u) * (p4 % v) +
                       p3p4 * (p3 % u) * (p4 % v)) +
            FACTOR2 * ((p4 % u) * (p3 % v) - p3p4 * g[u][v]);
  FOR_EACH_2_END;

  return T;
}

// Vector (Massive) - Pseudoscalar - Pseudoscalar vertex function
// i\Gamma_\mu(k1,k2)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_vpsps(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                          const regge::FFParam &ff_decay) const {
  const M4Vec                p      = k1 - k2;
  const regge::FFParam      &decay  = ff_decay;
  const double               F_     = regge::MassFF((k1 + k2).M2(), pow2(M0), decay);
  const std::complex<double> FACTOR = -0.5 * zi * g1 * F_;

  Tensor1<std::complex<double>, 4> T;
  for (const auto &mu : LI) { T(mu) = FACTOR * (p % mu); }
  return T;
}

// pseudoscalar - vector - vector vertex function
// i\Gamma_{\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors, M0 is the pseudoscalar
// on-shell mass
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_psvv(const M4Vec &p3, const M4Vec &p4, double M0, double g1,
                                                            const regge::FFParam &ff_decay) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const regge::FFParam      &decay  = ff_decay;
  const std::complex<double> FACTOR = zi * g1 / (2 * S0) * regge::MassFF((p3 + p4).M2(), pow2(M0), decay);

  Tensor1<std::complex<double>, 4> p3_ = {p3[0], p3[1], p3[2], p3[3]};
  Tensor1<std::complex<double>, 4> p4_ = {p4[0], p4[1], p4[2], p4[3]};

  FTensor::Index<'a', 4> mu;
  FTensor::Index<'b', 4> nu;
  FTensor::Index<'c', 4> rho;
  FTensor::Index<'d', 4> sigma;

  // Contract
  Tensor2<std::complex<double>, 4, 4> T;
  T(mu, nu) = FACTOR * p3_(rho) * p4_(sigma) * eps_lo(mu, nu, rho, sigma);

  return T;
}

// f2 - pseudoscalar - pseudoscalar vertex function
// i\Gamma_{\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors, M0 is the f2 meson on-shell
// mass
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_f2psps(const M4Vec &k1, const M4Vec &k2, double M0, double g1,
                                                              const regge::FFParam &ff_decay) const {
  const double               S0     = 1.0;  // Mass scale (GeV)
  const M4Vec                k1_k2  = k1 - k2;
  const regge::FFParam      &decay  = ff_decay;
  const std::complex<double> FACTOR = -zi * g1 / (2 * S0) * regge::MassFF((k1 + k2).M2(), pow2(M0), decay);

  const double k1_k2_sq = k1_k2.M2();

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR * ((k1_k2 % u) * (k1_k2 % v) - 0.25 * g[u][v] * k1_k2_sq);
  FOR_EACH_2_END;

  return T;
}

// f2 - Vector (massive) - Vector (Massive) vertex function
// i\Gamma_{\mu\nu\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors, M0 is the f2 meson on-shell
// mass
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_f2vv(const M4Vec &k1, const M4Vec &k2, double M0,
                                                                  double g1, double g2,
                                                                  const regge::FFParam &ff_decay) const {
  const double S0 = 1.0;  // Mass scale (GeV)

  const Tensor4<std::complex<double>, 4, 4, 4, 4> G0 = Gamma0(k1, k2);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> G2 = Gamma2(k1, k2);

  const regge::FFParam      &decay   = ff_decay;
  const std::complex<double> F_      = regge::MassFF((k1 + k2).M2(), pow2(M0), decay);
  const std::complex<double> FACTOR1 = zi * 2.0 / math::pow3(S0) * g1 * F_;
  const std::complex<double> FACTOR2 = -zi * 1.0 / S0 * g2 * F_;

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR1 * G0(u, v, k, l) + FACTOR2 * G2(u, v, k, l);
  FOR_EACH_4_END;

  return T;
}

// f2 - gamma - gamma vertex function
// i\Gamma_{\mu\nu\kappa\lambda}(k1,k2)
//
// Input as contravariant (upper index) 4-vectors, M0 is the f2 meson on-shell
// mass Couplings g1,g2
//
// Example couplings:
//
// const double e = msqrt(qed::alpha_QED() * 4.0 * PI);     // ~ 0.3, no running
// const double a_f2yy = pow2(e) / (4 * PI) * 1.45;  // GeV^{-3}
// const double b_f2yy = pow2(e) / (4 * PI) * 2.49;  // GeV^{-1}

Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_f2yy(const M4Vec &k1, const M4Vec &k2, double M0,
                                                                  double g1, double g2,
                                                                  const regge::FFParam &ff_decay) const {
  const Tensor4<std::complex<double>, 4, 4, 4, 4> G0 = Gamma0(k1, k2);
  const Tensor4<std::complex<double>, 4, 4, 4, 4> G2 = Gamma2(k1, k2);

  const regge::FFParam      &decay  = ff_decay;
  const std::complex<double> FACTOR = zi * regge::TransferFF(k1.M2(), tensor_param.meson_ff_transfer) *
                                      regge::TransferFF(k2.M2(), tensor_param.meson_ff_transfer) *
                                      regge::MassFF((k1 + k2).M2(), pow2(M0), decay);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (2.0 * g1 * G0(u, v, k, l) - g2 * G2(u, v, k, l));
  FOR_EACH_4_END;

  return T;
}

// Pomeron-Pseudoscalar-Pseudoscalar vertex function
// i\Gamma_{\mu\nu}(p',p)
//
// Input as contravariant (upper index) 4-vectors
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_Ppsps(const M4Vec &prime, const M4Vec &p, double g1) const {
  const M4Vec                psum = prime + p;
  const double               M2   = psum.M2();
  const std::complex<double> FACTOR =
      -zi * 2.0 * g1 * regge::TransferFF((prime - p).M2(), tensor_param.meson_ff_transfer);

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = FACTOR * ((psum % u) * (psum % v) - 0.25 * g[u][v] * M2);
  FOR_EACH_2_END;

  return T;
}

// Compute one exchange-specific pseudoscalar continuum vertex
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.53), (3.56)]
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_Tpsps(const M4Vec &prime, const M4Vec &p, const int exchange_pdg,
                                                             const int hadron_pdg) const {
  const auto  &vertex               = tensor_param.exchange.FindVertex(exchange_pdg, hadron_pdg);
  const double coupling             = tensor_param.exchange.PseudoscalarCoupling(exchange_pdg, vertex.g_tensor.front());
  const M4Vec  psum                 = prime + p;
  const double mass2                = psum.M2();
  const std::complex<double> factor = -zi * 2.0 * coupling * regge::TransferFF((prime - p).M2(), vertex.ff_transfer);

  Tensor2<std::complex<double>, 4, 4> out;
  FOR_EACH_2(LI);
  out(u, v) = factor * ((psum % u) * (psum % v) - 0.25 * g[u][v] * mass2);
  FOR_EACH_2_END;
  return out;
}

// Compute one vector-Reggeon pseudoscalar continuum vertex
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eq. (3.63)]
Tensor1<std::complex<double>, 4> MTensorPomeron::iG_Vpsps(const M4Vec &prime, const M4Vec &p, const int exchange_pdg,
                                                          const int hadron_pdg) const {
  const auto  &vertex        = tensor_param.exchange.FindVertex(exchange_pdg, hadron_pdg);
  const double coupling      = tensor_param.exchange.VectorPseudoscalarCoupling(exchange_pdg, vertex.g_tensor.front());
  const double particle_sign = hadron_pdg > 0 ? 1.0 : -1.0;
  const M4Vec  psum          = prime + p;
  const std::complex<double> factor =
      -zi * particle_sign * coupling * regge::TransferFF((prime - p).M2(), vertex.ff_transfer);
  Tensor1<std::complex<double>, 4> out;
  for (const auto &mu : LI) { out(mu) = factor * (psum % mu); }
  return out;
}

// Pomeron-Vector(massive)-Vector(massive) vertex function
// i\Gamma_{\alpha\beta\gamma\delta}(p',p)
//
// Input M0 is the vector meson on-shell mass
//
// Input as contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_Pvv(const M4Vec &prime, const M4Vec &p, double g1,
                                                                 double g2, const regge::FFParam &ff_transfer, bool forward) const {
  const std::complex<double> FACTOR = zi * regge::TransferFF(forward ? 0.0 : (prime - p).M2(), ff_transfer);
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = 0.0;
  FOR_EACH_4_END;

  if (ActiveProductionCoupling(g1)) {
    const auto G0 = Gamma0(prime, -p);
    FOR_EACH_4(LI);
    T(u, v, k, l) += FACTOR * 2.0 * g1 * G0(u, v, k, l);
    FOR_EACH_4_END;
  }
  if (ActiveProductionCoupling(g2)) {
    const auto G2 = Gamma2(prime, -p);
    FOR_EACH_4(LI);
    T(u, v, k, l) -= FACTOR * g2 * G2(u, v, k, l);
    FOR_EACH_4_END;
  }
  return T;
}

// Compute one exchange-specific vector continuum vertex
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iG_Tvv(const M4Vec &prime, const M4Vec &p,
                                                                 const int exchange_pdg, const int hadron_pdg) const {
  const auto &vertex = tensor_param.exchange.FindVertex(exchange_pdg, hadron_pdg);
  return iG_Pvv(prime, p, vertex.g_tensor[0], vertex.g_tensor[1], vertex.ff_transfer);
}

// Gamma-Vector meson transition vertex
// iGamma_{\mu \nu}
//
// VMD parameters are read from PARAM_TENSORPOM.VMD
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iG_yV(double q2, int pdg) const {
  const double e      = msqrt(qed::alpha_QED() * 4.0 * PI);  // ~ 0.3, no running
  const auto  &vmd    = tensor_param.FindVMD(pdg);
  const double gammaV = vmd.gammaV_sign * msqrt(4.0 * PI / vmd.gammaV2);
  const double mV     = vmd.mass;

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = -zi * e * pow2(mV) / gammaV * g[u][v];
  FOR_EACH_2_END;

  return T;
}

// ----------------------------------------------------------------------
// Propagators

// Tensor Pomeron propagator: (covariant == contravariant)
//
// i\Delta_{\mu\nu,\kappa\lambda}(s,t) = i\Delta^{\mu\nu,\kappa\lambda}(s,t)
//
// Input as (sub)-Mandelstam invariants: s,t
//
// Symmetry relations:
// \Delta_{\mu\nu,\kappa\lambda} = \Delta_{\nu\mu,\kappa\lambda} =
// \Delta_{\mu\nu,\lambda\kappa} = \Delta_{\kappa\lambda,\mu\nu}
//
// Contraction with g^{\mu\nu} or g^{\kappa\lambda} gives 0
//
// Compute the scalar Regge factor of the spin-2 Pomeron propagator
std::complex<double> MTensorPomeron::PomeronPropagatorFactor(const double s, const double t) const {
  return TensorPropagatorFactor(995, s, t);
}

// Compute the scalar Regge factor of one configured rank-two exchange
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.10), (3.12)]
std::complex<double> MTensorPomeron::TensorPropagatorFactor(const int exchange_pdg, const double s,
                                                            const double t) const {
  return tensor_param.exchange.PropagatorFactor(exchange_pdg, s, t);
}

// Contract one vector current with a configured vector propagator
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.14), (3.16)]
Tensor1<std::complex<double>, 4> MTensorPomeron::VectorPropagatorCurrent(
    const Tensor1<std::complex<double>, 4> &current, const int exchange_pdg, const double s, const double t) const {
  const auto &exchange = tensor_param.exchange.FindExchange(exchange_pdg);
  if (exchange.rank != 1) { throw std::invalid_argument("VectorPropagatorCurrent requires a rank-one exchange"); }
  const std::complex<double>       factor = TensorPropagatorFactor(exchange_pdg, s, t);
  Tensor1<std::complex<double>, 4> out;
  for (const auto &nu : LI) { out(nu) = factor * g[nu][nu] * current(nu); }
  return out;
}

// Contract one forward tensor current with the spin-2 Pomeron propagator
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::PomeronPropagatorCurrent(
    const Tensor2<std::complex<double>, 4, 4> &current, const double s, const double t) const {
  return PomeronPropagatorCurrent(current, PomeronPropagatorFactor(s, t));
}

// Project one tensor current with a precomputed Pomeron Regge factor
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::PomeronPropagatorCurrent(
    const Tensor2<std::complex<double>, 4, 4> &current, const std::complex<double> &factor) const {
  std::complex<double> trace = 0.0;
  for (const auto &mu : LI) { trace += g[mu][mu] * current(mu, mu); }

  Tensor2<std::complex<double>, 4, 4> out;
  for (const auto &alpha : LI) {
    for (const auto &beta : LI) {
      out(alpha, beta) = factor * (g[alpha][alpha] * g[beta][beta] * (current(alpha, beta) + current(beta, alpha)) -
                                   0.5 * g[alpha][beta] * trace);
    }
  }
  return out;
}

// Contract one Pomeron current into a baryon vertex Dirac matrix
MMatrix<std::complex<double>> MTensorPomeron::PomeronBaryonCurrent(const Tensor2<std::complex<double>, 4, 4> &current,
                                                                   const M4Vec &prime, const M4Vec &p,
                                                                   const double gPBB, const regge::FFParam &ff_transfer) const {
  const M4Vec                         psum  = prime + p;
  const MMatrix<std::complex<double>> slash = FSlash(psum);
  MMatrix<std::complex<double>>       out(4, 4, 0.0);

  for (const auto &mu : LI) {
    for (const auto &nu : LI) {
      const std::complex<double> weight = current(mu, nu);
      for (const auto &row : indices(gamma_lo[mu].Row(0))) {
        for (const auto &column : indices(out.Row(row))) {
          out[row][column] +=
              weight * (0.5 * (gamma_lo[mu][row][column] * (psum % nu) + gamma_lo[nu][row][column] * (psum % mu)) -
                        0.25 * g[mu][nu] * slash[row][column]);
        }
      }
    }
  }
  out *= -zi * 3.0 * gPBB * regge::TransferFF((prime - p).M2(), ff_transfer);
  return out;
}

// Contract one rank-two exchange current into a continuum baryon vertex
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.43), (3.49)]
MMatrix<std::complex<double>> MTensorPomeron::TensorBaryonCurrent(const Tensor2<std::complex<double>, 4, 4> &current,
                                                                  const M4Vec &prime, const M4Vec &p,
                                                                  const int exchange_pdg, const int baryon_pdg) const {
  const auto  &vertex   = tensor_param.exchange.FindVertex(exchange_pdg, baryon_pdg);
  const double coupling = tensor_param.exchange.BaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  return PomeronBaryonCurrent(current, prime, p, coupling, vertex.ff_transfer);
}

// Contract one vector exchange current into a continuum baryon vertex
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.59), (3.68)]
// [REFERENCE: Lebiedowicz et al., Phys. Rev. D 97, 094027 (2018), Eqs. (2.10) to (2.12)]
MMatrix<std::complex<double>> MTensorPomeron::VectorBaryonCurrent(const Tensor1<std::complex<double>, 4> &current,
                                                                  const M4Vec &prime, const M4Vec &p,
                                                                  const int exchange_pdg, const int baryon_pdg) const {
  const auto  &vertex   = tensor_param.exchange.FindVertex(exchange_pdg, baryon_pdg);
  const double coupling = tensor_param.exchange.VectorBaryonCoupling(exchange_pdg, vertex.g_tensor.front());
  MMatrix<std::complex<double>> out(4, 4, 0.0);
  for (const auto &mu : LI) { out += gamma_lo[mu] * current(mu); }
  // Both crossed vertices belong to one fermion line, so the Dirac chain
  // supplies its charge conjugation sign without an external antibaryon sign
  out *= -zi * coupling * regge::TransferFF((prime - p).M2(), vertex.ff_transfer);
  return out;
}

// Build the full spin-2 Pomeron propagator tensor
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_P(double s, double t) const {
  const std::complex<double> FACTOR = PomeronPropagatorFactor(s, t);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Tensor Reggeon propagator iD (f2-a2-reggeon)
//
// Input as (sub)-Mandelstam invariants: s,t
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_2R(double s, double t) const {
  const std::complex<double> FACTOR = TensorPropagatorFactor(9915, s, t);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = FACTOR * (g[u][k] * g[v][l] + g[u][l] * g[v][k] - 0.5 * g[u][v] * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Build one configured vector Odderon or vector Reggeon propagator
// [REFERENCE: Ewerz et al., Annals Phys. 342 (2014) 31, Eqs. (3.14), (3.16)]
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_VExchange(const int exchange_pdg, const double s,
                                                                 const double t) const {
  if (tensor_param.exchange.FindExchange(exchange_pdg).rank != 1) {
    throw std::invalid_argument("iD_VExchange requires a configured rank-one exchange");
  }
  const std::complex<double>          factor = TensorPropagatorFactor(exchange_pdg, s, t);
  Tensor2<std::complex<double>, 4, 4> out;
  FOR_EACH_2(LI);
  out(u, v) = factor * g[u][v];
  FOR_EACH_2_END;
  return out;
}

// Vector Reggeon propagator iD (rho-reggeon, omega-reggeon)
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_1R(double s, double t) const { return iD_VExchange(9933, s, t); }

// Odderon propagator iD
//
// Input as (sub)-Mandelstam invariants s,t
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_O(double s, double t) const { return iD_VExchange(9993, s, t); }

// Vector propagator iD_{\mu \nu} == iD^{\mu \nu} with (naive) Reggeization
//
// Input as contravariant 4-vector and the particle peak mass M0
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_V(const M4Vec &p, double M0, double s34, int pdg) const {
  const auto  &vector_param = tensor_param.FindVector(pdg);
  const double m2           = p.M2();
  const double alpha        = vector_param.trajectory_intercept + vector_param.trajectory_slope * m2;
  const double s_thresh     = 4 * pow2(M0);

  if (!std::isfinite(m2) || !std::isfinite(s34) || !(s34 > 0.0)) {
    throw AmplitudeFailure("MTensorPomeron::iD_V: invalid vector transfer or pair invariant mass");
  }

  // Apply the published species-specific Regge prescription above pair
  // threshold
  std::complex<double> reggeize = 1.0;
  if (s34 > s_thresh) {
    const double ratio = s34 / s_thresh;
    if (vector_param.regge_phase == TensorVectorReggePhase::None) {
      reggeize = std::pow(ratio, alpha - 1.0);
    } else {
      const double phase = PI / 2.0 * std::exp((s_thresh - s34) / s_thresh) - PI / 2.0;
      reggeize           = std::pow(std::exp(zi * phase) * ratio, alpha - 1.0);
    }
  }

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = -zi * g[u][v] / (m2 - pow2(M0)) * reggeize;
  FOR_EACH_2_END;

  return T;
}

// Stable transverse vector propagator for an off-shell VMD transition state
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_VMD(const M4Vec &p, int pdg) const {
  const double mass = tensor_param.FindVMD(pdg).mass;

  Tensor2<std::complex<double>, 4, 4> T;
  FOR_EACH_2(LI);
  T(u, v) = -zi * g[u][v] / (p.M2() - pow2(mass));
  FOR_EACH_2_END;
  return T;
}

// Scalar propagator iD with zero width
//
// Input as 4-vector and the particle peak mass M0
//
std::complex<double> MTensorPomeron::iD_MES0(const M4Vec &p, double M0) const { return zi / (p.M2() - pow2(M0)); }

// Scalar propagator iD with finite Breit-Wigner width
//
// Input as 4-vector and the particle peak mass M0 and full width
//
std::complex<double> MTensorPomeron::iD_MES(const M4Vec &p, double M0, double Gamma) const {
  const std::complex<double> Delta = 1.0 / (p.M2() - pow2(M0) + zi * M0 * Gamma);

  return zi * Delta;
}

// Compute the selected transverse spectral function, retaining bare initial VMD poles
std::complex<double> MTensorPomeron::VectorSpectrum(const double s, const double mass, const double width,
                                                    const int pdg) const {
  const auto &vector = tensor_param.FindVector(pdg);
  if (!std::isfinite(s)) { throw AmplitudeFailure("MTensorPomeron::VectorSpectrum: non-finite invariant mass"); }
  if (vector.width_model == TensorVectorWidthModel::RhoOmega) {
    if (s > 0.0) {
      const auto row = tensor_param.rho_omega.Index(pdg);
      return tensor_param.rho_omega.Propagator(s)[row][row];
    }
    // [REFERENCE: Lebiedowicz et al., arXiv:2508.06334v2, Eq. (2.30)]
    return 1.0 / (s - pow2(mass));
  }
  double imaginary = 0.0;
  if (vector.width_model == TensorVectorWidthModel::PWaveTwoBody) {
    const double threshold = 4.0 * pow2(vector.decay_daughter_mass);
    // [REFERENCE: Bolz et al., arXiv:1409.8483, Eqs. (B.20)-(B.21), rho-prime prescription]
    if (s > threshold) {
      imaginary = width * pow2(mass) / std::sqrt(s) * std::pow((s - threshold) / (pow2(mass) - threshold), 1.5);
    }
  } else if (s > 0.0) {
    imaginary = mass * width;
  }
  return 1.0 / (s - pow2(mass) + zi * imaginary);
}

// Contract the transverse propagator with the conserved pion or kaon decay current
Tensor1<std::complex<double>, 4> MTensorPomeron::VectorDecay(const M4Vec &positive, const M4Vec &negative,
    const double mass, const double width, const int pdg, const double coupling, const regge::FFParam &form) const {
  RequireConservedVectorDecay(positive, negative, "MTensorPomeron::VectorDecay");
  const double s = (positive + negative).M2();
  std::complex<double> spectral;
  if (tensor_param.FindVector(pdg).width_model == TensorVectorWidthModel::RhoOmega) {
    // [REFERENCE: Lebiedowicz et al., arXiv:2508.06334v2, Eqs. (2.25)-(2.28)]
    const auto &mixing = tensor_param.rho_omega;
    const auto delta = mixing.Propagator(s);
    const auto row = mixing.Index(pdg);
    spectral = delta[row][0] * mixing.g[0] + delta[row][1] * mixing.g[1];
  } else {
    spectral = VectorSpectrum(s, mass, width, pdg) * coupling;
  }
  const auto factor = -0.5 * spectral * regge::MassFF(s, pow2(mass), form);
  const M4Vec difference = positive - negative;
  Tensor1<std::complex<double>, 4> out;
  for (const auto &mu : LI) { out(mu) = factor * difference[mu]; }
  return out;
}

// Massive vector propagator iD_{\mu \nu} or iD^{\mu \nu} with INDEX_UP == true
//
// Input as contravariant 4-vector and the particle peak mass M0 and full width
// Gamma
//
Tensor2<std::complex<double>, 4, 4> MTensorPomeron::iD_VMES(const M4Vec &p, double M0, double Gamma, int pdg,
                                                            bool INDEX_UP, bool CONSERVED_CURRENT) const {
  const double m2 = p.M2();
  const auto Delta_T = VectorSpectrum(m2, M0, Gamma, pdg);

  // Longitudinal part (does not enter here)
  const double Delta_L = 0.0;

  Tensor2<std::complex<double>, 4, 4> T;
  if (CONSERVED_CURRENT) {
    FOR_EACH_2(LI);
    T(u, v) = -zi * g[u][v] * Delta_T;
    FOR_EACH_2_END;
    return T;
  }

  const double component_scale2 = std::max({1.0, pow2(p.E()), pow2(p.P3mod())});
  if (std::abs(m2) <= std::numeric_limits<double>::epsilon() * component_scale2) {
    throw AmplitudeFailure(
        "MTensorPomeron::iD_VMES: null momentum "
        "requires conserved-current contraction");
  }

  if (INDEX_UP) {
    FOR_EACH_2(LI);
    T(u, v) = zi * (-g[u][v] + (p[u]) * (p[v]) / m2) * Delta_T - zi * (p[u]) * (p[v]) / m2 * Delta_L;
    FOR_EACH_2_END;
  } else {
    FOR_EACH_2(LI);
    T(u, v) = zi * (-g[u][v] + (p % u) * (p % v) / m2) * Delta_T - zi * (p % u) * (p % v) / m2 * Delta_L;
    FOR_EACH_2_END;
  }
  return T;
}

// Massive tensor propagator iD_{\mu \nu \kappa \lambda} or iD^{} with INDEX_UP
// == true
//
// Input as contravariant (upper index) 4-vector and the particle peak mass M0
// and full width Gamma
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::iD_TMES(const M4Vec &p, double M0, double Gamma,
                                                                  bool INDEX_UP) const {
  const double m2 = p.M2();

  const double scale = std::max({1.0, pow2(p.E()), pow2(p.P3mod())});
  if (!std::isfinite(m2) || !std::isfinite(scale) ||
      std::abs(m2) <= std::numeric_limits<double>::epsilon() * scale) {
    throw AmplitudeFailure("MTensorPomeron::iD_TMES: null or non-finite tensor invariant mass");
  }

  // Pre-calculated auxialary tensor
  Tensor2<double, 4, 4> ghat;
  FOR_EACH_2(LI);
  if (INDEX_UP) {
    ghat(u, v) = -g[u][v] + (p[u]) * (p[v]) / m2;
  } else {
    ghat(u, v) = -g[u][v] + (p % u) * (p % v) / m2;
  }
  FOR_EACH_2_END;

  const std::complex<double> Delta = 1.0 / (m2 - pow2(M0) + zi * M0 * Gamma);

  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) =
      zi * Delta * (0.5 * (ghat(u, k) * ghat(v, l) + ghat(u, l) * ghat(v, k)) - 1.0 / 3.0 * ghat(u, v) * ghat(k, l));
  FOR_EACH_4_END;

  return T;
}

// ----------------------------------------------------------------------
// Trajectories (simple linear/affine here)

// Pomeron trajectory
double MTensorPomeron::alpha_P(double t) const {
  const auto &exchange = tensor_param.exchange.FindExchange(995);
  return 1.0 + exchange.delta + exchange.ap * t;
}
// Odderon trajectory
double MTensorPomeron::alpha_O(double t) const {
  const auto &exchange = tensor_param.exchange.FindExchange(9993);
  return 1.0 + exchange.delta + exchange.ap * t;
}
// Reggeon (rho,omega) trajectory
double MTensorPomeron::alpha_1R(double t) const {
  const auto &exchange = tensor_param.exchange.FindExchange(9933);
  return 1.0 + exchange.delta + exchange.ap * t;
}
// Reggeon (f2,a2) trajectory
double MTensorPomeron::alpha_2R(double t) const {
  const auto &exchange = tensor_param.exchange.FindExchange(9915);
  return 1.0 + exchange.delta + exchange.ap * t;
}

// ----------------------------------------------------------------------
// Form factors

// Proton electromagnetic Dirac form FACTOR
double MTensorPomeron::F1_(double t) const { return form::F1(t, model_tune->Structure()); }

// Proton electromagnetic Pauli form FACTOR
double MTensorPomeron::F2_(double t) const { return form::F2(t, model_tune->Structure()); }

// Dipole form FACTOR
double MTensorPomeron::GD(double t) const {
  const double m2D = 0.71;  // GeV^2
  return 1.0 / pow2(1 - t / m2D);
}

// ----------------------------------------------------------------------
// Tensor functions

// Output: \Gamma^{(0)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::Gamma0(const M4Vec &k1, const M4Vec &k2) const {
  const double                              k1k2 = k1 * k2;  // 4-dot product
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) =
      (k1k2 * g[u][v] - (k2 % u) * (k1 % v)) * ((k1 % k) * (k2 % l) + (k2 % k) * (k1 % l) - 0.5 * k1k2 * g[k][l]);
  FOR_EACH_4_END;

  return T;
}

// Output: \Gamma^{(2)}_{\mu\nu\kappa\lambda}
//
// Input must be contravariant (upper index) 4-vectors
//
Tensor4<std::complex<double>, 4, 4, 4, 4> MTensorPomeron::Gamma2(const M4Vec &k1, const M4Vec &k2) const {
  const double                              k1k2 = k1 * k2;  // 4-dot product
  Tensor4<std::complex<double>, 4, 4, 4, 4> T;
  FOR_EACH_4(LI);
  T(u, v, k, l) = k1k2 * (g[u][k] * g[v][l] + g[u][l] * g[v][k]) +
                  g[u][v] * ((k1 % k) * (k2 % l) + (k2 % k) * (k1 % l)) - (k1 % v) * (k2 % l) * g[u][k] -
                  (k1 % v) * (k2 % k) * g[u][l] - (k2 % u) * (k1 % l) * g[v][k] - (k2 % u) * (k1 % k) * g[v][l] -
                  (k1k2 * g[u][v] - (k2 % u) * (k1 % v)) * g[k][l];
  FOR_EACH_4_END;

  return T;
}

// Pre-calculate tensors for speed
//
// 1. Minkowski metric tensor
// 2. Auxialary (helper) tensor
// 1/2 g_{\mu\kappa} g_{\nu\lambda} + 1/2 g_{\mu\lambda}g_{\nu\kappa} - 1/4
// g_{\mu\nu}g_{kappa\lambda}
//
void MTensorPomeron::CalcRTensor() {
  // Minkowski metric tensor (+1,-1,-1,-1)
  FOR_EACH_2(LI);
  gT(u, v) = g[u][v];
  FOR_EACH_2_END;

  // ------------------------------------------------------------------
  // Covariant [lower index] 4D-epsilon tensor

  const MTensor<int> etensor = math::EpsTensor(4);

  FOR_EACH_4(LI);
  eps_lo(u, v, k, l) = static_cast<double>(etensor(u, v, k, l));
  FOR_EACH_4_END;

  // Free Lorentz indices [second parameter denotes the range of index]
  FTensor::Index<'a', 4> a;
  FTensor::Index<'b', 4> b;
  FTensor::Index<'c', 4> c;
  FTensor::Index<'d', 4> d;
  FTensor::Index<'g', 4> alfa;
  FTensor::Index<'h', 4> beta;

  // Contravariant version (make it in two steps)
  eps_hi(a, b, c, d) = eps_lo(alfa, beta, c, d) * gT(a, alfa) * gT(b, beta);
  eps_hi(a, b, c, d) = eps_hi(a, b, alfa, beta) * gT(c, alfa) * gT(d, beta);

  // ------------------------------------------------------------------

  // Aux tensor R
  FTensor::Tensor4<double, 4, 4, 4, 4> R;

  R(a, b, c, d) = 0.5 * gT(a, c) * gT(b, d) + 0.5 * gT(a, d) * gT(b, c) - 0.25 * gT(a, b) * gT(c, d);

  // Different mixed index covariant/contravariant versions
  // by contraction with g_{\mu\nu}
  R_DDDD             = R;
  R_DDDU(a, b, c, d) = R(a, b, c, alfa) * gT(d, alfa);
  R_DDUU(a, b, c, d) = R(a, b, alfa, beta) * gT(c, alfa) * gT(d, beta);
  R_UUDD(a, b, c, d) = R(alfa, beta, c, d) * gT(a, alfa) * gT(b, beta);

  // -------------------------------------------------------------------
  // Pre-calculated kinematics-independent resonance tensors

  {
    const double               S0     = 1.0;  // Mass scale (GeV)
    const std::complex<double> FACTOR = 2.0 * zi * S0;

    // Init with zeros!
    PPT_00 = MTensor<std::complex<double>>({4, 4, 4, 4, 4, 4}, std::complex<double>(0.0));

    FOR_EACH_6(LI);

    for (auto const &v1 : LI) {
      for (auto const &a1 : LI) {
        for (auto const &l1 : LI) {
          for (auto const &r1 : LI) {
            for (auto const &s1 : LI) {
              for (auto const &u1 : LI) {
                if (v1 != a1 || l1 != r1 || s1 != u1) { continue; }  // speed it up

                // Note +=
                PPT_00(u, v, k, l, r, s) += FACTOR * R_DDDD(u, v, u1, v1) * R_DDDD(k, l, a1, l1) *
                                            R_DDDD(r, s, r1, s1) * g[v1][a1] * g[l1][r1] * g[s1][u1];
              }
            }
          }
        }
      }
    }
    FOR_EACH_6_END;
  }

  // Symmetrized rank-8 tensor entering the axial (2,2) coupling
  PPA_22 = MTensor<double>({4, 4, 4, 4, 4, 4, 4, 4}, double(0.0));
  for (const auto kappa : LI) {
    for (const auto lambda : LI) {
      for (const auto rho : LI) {
        for (const auto sigma : LI) {
          for (const auto mu : LI) {
            for (const auto nu : LI) {
              for (const auto alpha : LI) {
                for (const auto beta : LI) {
                  PPA_22({kappa, lambda, rho, sigma, mu, nu, alpha, beta}) =
                      GammaPPA22Base(kappa, lambda, rho, sigma, mu, nu, alpha, beta) +
                      GammaPPA22Base(kappa, lambda, rho, sigma, nu, mu, alpha, beta);
                }
              }
            }
          }
        }
      }
    }
  }

  // -------------------------------------------------------------------
  // Pre-calculated tensor contractions for tensor resonances

  auto FIX1 = [&]() {
    MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4}, double(0.0));  // Init with zeros!
    FOR_EACH_6(LI);
    for (auto const &a1 : LI) {
      for (auto const &r1 : LI) {
        for (auto const &s1 : LI) {
          A(u, v, k, l, r, s) += R_DDDD(u, v, r1, a1) * R_DDDU(k, l, s1, a1) * R_DDUU(r, s, r1, s1);
        }
      }
    }
    FOR_EACH_6_END;
    return A;
  };

  auto FIX2 = [&]() {
    MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4, 4, 4}, double(0.0));  // Init with zeros!
    FOR_EACH_6(LI);
    for (auto const &a1 : LI) {
      for (auto const &r1 : LI) {
        for (auto const &s1 : LI) {
          for (auto const &u1 : LI) {
            A(u, v, k, l, r, s, u1, r1) += R_DDDD(u, v, u1, a1) * R_DDDU(k, l, s1, a1) * R_DDUU(r, s, r1, s1);
          }
        }
      }
    }
    FOR_EACH_6_END;
    return A;
  };

  auto FIX3 = [&]() {
    MTensor<double> A = MTensor({4, 4, 4, 4, 4, 4, 4, 4}, double(0.0));  // Init with zeros!
    FOR_EACH_6(LI);
    for (auto const &a1 : LI) {
      for (auto const &r1 : LI) {
        for (auto const &s1 : LI) {
          for (auto const &u1 : LI) {
            A(u, v, k, l, r, s, u1, s1) += R_DDDD(u, v, r1, a1) * R_DDDU(k, l, u1, a1) * R_DDUU(r, s, r1, s1);
          }
        }
      }
    }
    FOR_EACH_6_END;
    return A;
  };

  T1 = FIX1();
  T2 = FIX2();
  T3 = FIX3();
}

}  // namespace gra
