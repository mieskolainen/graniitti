// Serial four-particle and six-particle Regge ladders
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeMulti.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <optional>
#include <span>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {

namespace ladder {

// Compute one prepared local continuum pole by its exact structural key
const ReggeContinuumPole &Pole(const LORENTZSCALAR &lts, const ReggeContinuumPoleKey &key, const std::string &context) {
  const auto found = lts.process.CONT_LADDER_POLE.find(key);
  if (found == lts.process.CONT_LADDER_POLE.end()) { throw std::logic_error(context + ": missing prepared continuum spin key"); }
  return found->second;
}

// Compute the model-specific local continuum contraction basis
std::vector<double> Basis(const ReggeContinuumPole &pole, const ReggeProductionModel model, const std::size_t leg) {
  if (model == ReggeProductionModel::GP) { return gpom::LocalBasis(pole, leg); }
  if (model == ReggeProductionModel::MP || model == ReggeProductionModel::XP) { return (model == ReggeProductionModel::MP ? mpom::LocalBasis : xpom::LocalBasis)(pole, leg); }
  throw std::invalid_argument("MRegge ladder basis requires MP, XP or GP");
}

// Evaluate one model-specific local final-pair kernel
MMatrix<std::complex<double>> Kernel(const ReggeContinuumPole &pole, const regge::Param &param, const ReggeProductionModel model, const MDecayBranch &first, const MDecayBranch &second, const M4Vec &upper_exchange, const M4Vec &lower_exchange,
                                     const double alpha_upper, const double alpha_lower) {
  if (model == ReggeProductionModel::GP) { return gpom::PairKernel(pole, param, first, second, alpha_upper, alpha_lower); }
  if (model == ReggeProductionModel::MP || model == ReggeProductionModel::XP) { return (model == ReggeProductionModel::MP ? mpom::PairKernel : xpom::PairKernel)(pole, first, second, upper_exchange, lower_exchange); }
  throw std::invalid_argument("MRegge ladder kernel requires MP, XP or GP");
}

// Validate that adjacent kernels use the same retained spin basis
void CheckBasis(const std::vector<double> &upper, const std::vector<double> &lower, const std::string &context) {
  if (upper.size() != lower.size()) { throw std::invalid_argument(context + ": connected spin dimensions differ"); }
  for (const auto &i : indices(upper)) {
    if (std::abs(upper[i] - lower[i]) > 1.0e-9) { throw std::invalid_argument(context + ": connected helicity bases differ"); }
  }
}

// Contract a serial kernel chain between its two boundary vertices
std::complex<double> Contract(const std::vector<std::complex<double>> &upper, const std::span<const MMatrix<std::complex<double>>> kernels, const std::span<const MMatrix<std::complex<double>>> metrics,
                              const std::vector<std::complex<double>> &lower) {
  MMatrix<std::complex<double>> chain = kernels.front();
  for (const auto &i : indices(metrics)) { chain = chain * metrics[i] * kernels[i + 1]; }
  return chain.BilinearForm(upper, lower);
}

}  // namespace ladder

namespace regge {

// Build one immutable serial continuum graph and exact cache-slot plan
template <std::size_t N>
ContinuumPlan BuildContinuumPlan(const LORENTZSCALAR &lts, const Param &param, const ReggeProductionModel model) {
  static_assert(N == 4 || N == 6);
  const std::string context = N == 4 ? "MRegge::EvalContinuum4" : "MRegge::EvalContinuum6";
  ContinuumPlan out;
  out.central_count = N;
  out.permutations  = lts.process.CONT_LADDER_PERMUTATIONS.empty() ? LadderPermutations(lts, param, N, model) : lts.process.CONT_LADDER_PERMUTATIONS;
  const auto configured = param.con.multiregge_topologies.find(static_cast<int>(N));
  if (configured == param.con.multiregge_topologies.end()) { throw std::invalid_argument(context + ": missing continuum topology bank"); }
  out.topologies = lts.process.MULTIREGGE_TOPOLOGIES.empty() ? configured->second : lts.process.MULTIREGGE_TOPOLOGIES;
  CheckLadderPermutations(out.permutations, N, context);
  CheckMultiReggeTopologies(out.topologies, N, context);

  using PairKey = std::array<int, 2>;
  std::map<PairKey, std::size_t> pair_index;
  const auto pair_for = [&](const int first, const int second) {
    const PairKey key = {first, second};
    if (const auto found = pair_index.find(key); found != pair_index.end()) { return found->second; }
    ContinuumPairPlan pair;
    pair.first       = DecayIndex(first, lts, context);
    pair.second      = DecayIndex(second, lts, context);
    pair.param       = FindPair(param, {lts.decaytree[pair.first].p.pdg, lts.decaytree[pair.second].p.pdg}, model);
    pair.first_mass2 = pow2(lts.decaytree[pair.first].p.mass);
    if (pair.param == nullptr) { throw AmplitudeFailure(context + ": resolved continuum ladder lost its " + ReggeProductionModelName(model) + " entry"); }
    for (const auto &row : pair.param->channels) {
      if (!VertexAllowed(param, row)) { continue; }
      ContinuumVertexPlan vertex;
      vertex.param            = &row;
      vertex.upper_alias      = row.first;
      vertex.lower_alias      = row.second;
      vertex.upper_trajectory = TrajectoryIndex(param, row.first);
      vertex.lower_trajectory = TrajectoryIndex(param, row.second);
      const ReggeContinuumPoleKey pole_key = {lts.decaytree[pair.first].p.pdg, lts.decaytree[pair.second].p.pdg, row.first, row.second};
      vertex.pole = ladder::Pole(lts, pole_key, context);
      for (const auto &leg : indices(vertex.basis)) { vertex.basis[leg] = ladder::Basis(vertex.pole, model, leg); }
      pair.vertex.push_back(std::move(vertex));
    }
    const std::size_t index = out.pair.size();
    out.pair.push_back(std::move(pair));
    pair_index.emplace(key, index);
    return index;
  };
  const auto connected = [model](const ContinuumVertexPlan &upper, const ContinuumVertexPlan &lower) { return model == ReggeProductionModel::GP ? upper.lower_trajectory == lower.upper_trajectory : upper.lower_alias == lower.upper_alias; };
  const auto slot = [](auto &bank, const auto &key) { return bank.try_emplace(key, bank.size()).first->second; };
  std::array<std::map<PairKey, std::size_t>, 2> outer_meson;
  std::map<std::tuple<std::array<int, 3>, int, int>, std::size_t> middle_meson;
  std::map<std::tuple<int, int, int, int>, std::size_t> upper_internal;
  std::map<std::tuple<std::array<int, 4>, int, int, int>, std::size_t> lower_internal;
  std::map<std::tuple<int, int, int, int>, std::size_t> boundary_index;
  const auto boundary_for = [&](const int upper, const int first, const int lower, const int last) {
    const auto [found, inserted] = boundary_index.try_emplace(std::make_tuple(upper, first, lower, last), out.boundary.size());
    if (inserted) { out.boundary.push_back({upper, first, lower, last}); }
    return found->second;
  };

  for (const auto &permutation : out.permutations) {
    ContinuumDiagramPlan diagram;
    std::copy_n(permutation.cbegin(), N, diagram.order.begin());
    constexpr std::size_t pairs = N / 2;
    for (std::size_t i = 0; i < pairs; ++i) { diagram.pair[i] = pair_for(diagram.order[2 * i], diagram.order[2 * i + 1]); }
    diagram.meson_slot[0]         = slot(outer_meson[0], PairKey{diagram.order[0], diagram.order[1]});
    diagram.meson_slot[pairs - 1] = slot(outer_meson[1], PairKey{diagram.order[N - 2], diagram.order[N - 1]});
    if constexpr (N == 6) {
      std::array<int, 3> prefix = {diagram.order[0], diagram.order[1], diagram.order[2]};
      std::sort(prefix.begin(), prefix.end());
      diagram.meson_slot[1] = slot(middle_meson, std::make_tuple(prefix, diagram.order[2], diagram.order[3]));
    }
    const auto &first = out.pair[diagram.pair[0]].vertex;
    const auto &second = out.pair[diagram.pair[1]].vertex;
    for (const auto &i : indices(first)) {
      for (const auto &j : indices(second)) {
        if (!connected(first[i], second[j])) { continue; }
        ladder::CheckBasis(first[i].basis[1], second[j].basis[0], context);
        const std::size_t upper_slot = slot(upper_internal, std::make_tuple(diagram.order[0], diagram.order[1], diagram.order[2], first[i].lower_alias));
        if constexpr (N == 4) {
          diagram.chain.push_back({{i, j, 0}, {upper_slot, 0}, boundary_for(first[i].upper_alias, diagram.order[0], second[j].lower_alias, diagram.order[3])});
        } else {
          const auto &third = out.pair[diagram.pair[2]].vertex;
          std::array<int, 4> prefix = {diagram.order[0], diagram.order[1], diagram.order[2], diagram.order[3]};
          std::sort(prefix.begin(), prefix.end());
          for (const auto &k : indices(third)) {
            if (!connected(second[j], third[k])) { continue; }
            ladder::CheckBasis(second[j].basis[1], third[k].basis[0], context);
            const std::size_t lower_slot = slot(lower_internal, std::make_tuple(prefix, diagram.order[3], diagram.order[4], second[j].lower_alias));
            diagram.chain.push_back({{i, j, k}, {upper_slot, lower_slot}, boundary_for(first[i].upper_alias, diagram.order[0], third[k].lower_alias, diagram.order[5])});
          }
        }
      }
    }
    out.diagram.push_back(std::move(diagram));
  }
  out.meson_slots    = {outer_meson[0].size(), middle_meson.size(), outer_meson[1].size()};
  out.internal_slots = {upper_internal.size(), lower_internal.size()};
  return out;
}

// Build the immutable four-body continuum graph and cache-slot plan
ContinuumPlan BuildContinuum4Plan(const LORENTZSCALAR &lts, const Param &param, const ReggeProductionModel model) { return BuildContinuumPlan<4>(lts, param, model); }

// Build the immutable six-body continuum graph and cache-slot plan
ContinuumPlan BuildContinuum6Plan(const LORENTZSCALAR &lts, const Param &param, const ReggeProductionModel model) { return BuildContinuumPlan<6>(lts, param, model); }

}  // namespace regge

namespace {

// Complete one cached meson line with its ordered local vertex factor
std::complex<double> SerialPairFactor(const gra::LORENTZSCALAR &lts, const regge::Param &param, const ReggeProductionModel model, const regge::ContinuumPairPlan &pair, const regge::ContinuumVertexPlan &vertex, const std::complex<double> meson_exchange,
                                      const double meson_q2, const M4Vec &upper_exchange, const M4Vec &lower_exchange, const bool upper_ext, const bool lower_ext) {
  (void)model;
  const bool upper_ff = upper_ext ? param.con.multiregge_transfer_ext : param.con.multiregge_transfer_int;
  const bool lower_ff = lower_ext ? param.con.multiregge_transfer_ext : param.con.multiregge_transfer_int;
  return meson_exchange * regge::MesonVertexFactor(meson_q2, pair.first_mass2, *vertex.param) *
         regge::PairVertexFactor(param, *vertex.param, lts.decaytree[pair.first], lts.decaytree[pair.second], upper_exchange.M2(), lower_exchange.M2(), upper_ff, lower_ff);
}

// Contract one complete serial ladder through local kernels and spin metrics
std::complex<double> ContractLadderSpin(const gra::LORENTZSCALAR &lts, const regge::Param &param, const ReggeProductionModel model, const std::span<const regge::ContinuumPairPlan *const> pairs,
                                        const std::span<const regge::ContinuumVertexPlan *const> vertices, const std::span<const std::complex<double>> pair_scales, const std::span<const std::complex<double>> internal_propagators,
                                        const std::span<const M4Vec> exchange_momenta, const M4Vec &lower_transfer) {
  const std::size_t count = pairs.size();
  std::vector<MMatrix<std::complex<double>>> kernels;
  kernels.reserve(count);
  for (std::size_t index = 0; index < count; ++index) {
    const auto &cache = vertices[index]->pole;
    auto kernel = ladder::Kernel(cache, param, model, lts.decaytree[pairs[index]->first], lts.decaytree[pairs[index]->second],
                                 exchange_momenta[index], exchange_momenta[index + 1],
                                 regge::Alpha(param, vertices[index]->upper_alias, exchange_momenta[index].M2()),
                                 regge::Alpha(param, vertices[index]->lower_alias, exchange_momenta[index + 1].M2()));
    kernel *= pair_scales[index];
    kernels.push_back(std::move(kernel));
  }

  std::vector<MMatrix<std::complex<double>>> metrics;
  metrics.reserve(count - 1);
  for (std::size_t index = 0; index + 1 < count; ++index) {
    const auto &left_basis = vertices[index]->basis[1];
    auto        metric     = spin::ExchangeMetric(ExchangeBasisType::ReducedRegge, left_basis, 0.0);
    metric *= internal_propagators[index];
    metrics.push_back(std::move(metric));
  }

  const spin::ForwardSpec forward     = {lts.process.FORWARD_VERTEX, param.s0, ExchangeBasisType::ReducedRegge};
  const auto             &upper_basis = vertices.front()->basis[0];
  const auto             &lower_basis = vertices.back()->basis[1];
  const auto              upper       = spin::ExchangeBoundary(upper_basis, exchange_momenta.front(), false, forward);
  const auto              lower       = spin::ExchangeBoundary(lower_basis, lower_transfer, true, forward);
  return ladder::Contract(upper, kernels, metrics, lower);
}

// Compute one event-local value in its preassigned exact cache slot
template <typename Evaluator>
const std::complex<double> &CachedPlanValue(std::vector<std::optional<std::complex<double>>> &cache, const std::size_t slot, Evaluator &&evaluator) {
  if (!cache.at(slot)) { cache[slot] = evaluator(); }
  return *cache[slot];
}

// Evaluate the serial four- or six-body Feynman graphs in prepared order
std::vector<std::complex<double>> EvalSerialPlan(const MRegge &amplitude, const LORENTZSCALAR &lts, const regge::Param &param, const ReggeProductionModel model, const regge::ContinuumPlan &plan,
                                                 const ForwardLegState &upper_state, const ForwardLegState &lower_state) {
  const std::size_t pair_count = plan.central_count / 2;
  std::array<std::vector<std::optional<std::complex<double>>>, 3> meson_cache = {std::vector<std::optional<std::complex<double>>>(plan.meson_slots[0]), std::vector<std::optional<std::complex<double>>>(plan.meson_slots[1]),
                                                                                 std::vector<std::optional<std::complex<double>>>(plan.meson_slots[2])};
  std::array<std::vector<std::optional<std::complex<double>>>, 2> internal_cache = {std::vector<std::optional<std::complex<double>>>(plan.internal_slots[0]),
                                                                                    std::vector<std::optional<std::complex<double>>>(plan.internal_slots[1])};
  std::vector<std::complex<double>> boundary(plan.boundary.size(), 0.0);
  const M4Vec lower_transfer = lower_state.incoming - lower_state.outgoing;
  for (const auto &diagram : plan.diagram) {
    std::array<const regge::ContinuumPairPlan *, 3> pairs = {};
    for (std::size_t i = 0; i < pair_count; ++i) { pairs[i] = &plan.pair[diagram.pair[i]]; }
    std::array<M4Vec, 4> exchange;
    exchange[0] = upper_state.incoming - upper_state.outgoing;
    for (std::size_t i = 0; i < pair_count; ++i) { exchange[i + 1] = exchange[i] - lts.decaytree[pairs[i]->first].p4 - lts.decaytree[pairs[i]->second].p4; }

    std::array<double, 3> meson_q2 = {};
    meson_q2[0] = lts.tt_1[diagram.order[0]];
    meson_q2[pair_count - 1] = lts.tt_2[diagram.order[plan.central_count - 1]];
    std::array<std::complex<double>, 3> meson = {};
    meson[0] = CachedPlanValue(meson_cache[0], diagram.meson_slot[0], [&] { return regge::MesonPropagator(param, meson_q2[0], pairs[0]->first_mass2, *pairs[0]->param, lts.decaytree[pairs[0]->first], lts.decaytree[pairs[0]->second]); });
    meson[pair_count - 1] = CachedPlanValue(meson_cache[2], diagram.meson_slot[pair_count - 1], [&] {
      return regge::MesonPropagator(param, meson_q2[pair_count - 1], pairs[pair_count - 1]->first_mass2, *pairs[pair_count - 1]->param, lts.decaytree[pairs[pair_count - 1]->first], lts.decaytree[pairs[pair_count - 1]->second]);
    });
    if (pair_count == 3) {
      const M4Vec transfer = lts.q1 - lts.decaytree[pairs[0]->first].p4 - lts.decaytree[pairs[0]->second].p4 - lts.decaytree[pairs[1]->first].p4;
      meson_q2[1] = transfer.M2();
      meson[1] = CachedPlanValue(meson_cache[1], diagram.meson_slot[1], [&] { return regge::MesonPropagator(param, meson_q2[1], pairs[1]->first_mass2, *pairs[1]->param, lts.decaytree[pairs[1]->first], lts.decaytree[pairs[1]->second]); });
    }

    for (const auto &chain : diagram.chain) {
      std::array<const regge::ContinuumVertexPlan *, 3> vertices = {};
      std::array<std::complex<double>, 3> pair_scale = {};
      std::array<std::complex<double>, 2> internal = {};
      for (std::size_t i = 0; i < pair_count; ++i) {
        vertices[i]   = &pairs[i]->vertex[chain.vertex[i]];
        pair_scale[i] = SerialPairFactor(lts, param, model, *pairs[i], *vertices[i], meson[i], meson_q2[i], exchange[i], exchange[i + 1], i == 0, i + 1 == pair_count);
      }
      for (std::size_t i = 0; i + 1 < pair_count; ++i) {
        internal[i] = CachedPlanValue(internal_cache[i], chain.internal_slot[i], [&] {
          const double transfer = i == 0 ? lts.tt_xy[diagram.order[0]][diagram.order[1]] : exchange[i + 1].M2();
          return amplitude.ExchangeKernel(lts.ss[diagram.order[2 * i + 1]][diagram.order[2 * i + 2]], transfer, vertices[i]->lower_alias);
        });
      }
      boundary[chain.boundary] += ContractLadderSpin(lts, param, model, {pairs.data(), pair_count}, {vertices.data(), pair_count}, {pair_scale.data(), pair_count}, {internal.data(), pair_count - 1}, {exchange.data(), pair_count + 1}, lower_transfer);
    }
  }
  return boundary;
}

}  // namespace

// Evaluate one complete four- or six-particle central Regge amplitude
std::complex<double> MRegge::EvalContinuumLadder(gra::LORENTZSCALAR &lts, const ReggeProductionModel model, const std::size_t central_count) {
  const std::string context = central_count == 4 ? "MRegge::EvalContinuum4" : "MRegge::EvalContinuum6";
  lts.proton_good_walker.reset();
  const ForwardLegState upper_state = ResolveForwardLegState(lts, ForwardBeamLeg::Upper);
  const ForwardLegState lower_state = ResolveForwardLegState(lts, ForwardBeamLeg::Lower);

  for (const auto &topology : continuum_plan.topologies) {
    if (topology.size() > 1) { AccumulateParallelTopology(lts, model, topology, continuum_plan.permutations); }
  }

  const regge::Topology serial_topology = {static_cast<int>(central_count)};
  if (std::find(continuum_plan.topologies.cbegin(), continuum_plan.topologies.cend(), serial_topology) != continuum_plan.topologies.cend()) {
    const auto boundary = EvalSerialPlan(*this, lts, param, model, continuum_plan, upper_state, lower_state);
    for (const auto &i : indices(boundary)) {
      const auto &source = continuum_plan.boundary[i];
      AddLadderGoodWalker(lts, upper_state, lower_state, source.upper_exchange, lts.ss[1][source.first_particle], source.lower_exchange, lts.ss[2][source.last_particle], boundary[i]);
    }
  }

  if (!lts.proton_good_walker) { throw AmplitudeFailure(context + ": empty Good Walker source"); }
  ScaleGoodWalker(lts, SewingSign(central_count / 2), context);
  if (GoodWalkerOnly(lts, context)) { return 0.0; }
  ProjectGoodWalkerBorn(lts);
  return msqrt(BornNorm(lts));
}

// ============================================================================
// Regge matrix element ansatz for 2->6 with four central particles
//
// [4], one serial four-particle ladder:
//
//  ======F======>
//        *
//        *
//        ff-----> a
//        |
//        ff-----> b
//        *
//        *
//        ff-----> c
//        |
//        ff-----> d
//        *
//        *
//  ======F======>
//
// Charged-pion combinatorics with charge-balanced permutations:
// 16 ordered particle diagrams
//
// The upper and lower two-exchange proton vertices are combined through the
// order two coherent multi Reggeon production convolution
//
// A_[2,2](Q_T) = c_2 / (8*pi^2*s)
//   * sum_partitions sum_(sigma in S_2) integral d^2 k
//   * <pp| B_(sigma(1))(k) B_(sigma(2))(Q_T-k) |pp>
//
// with c_2 = i/2. The two attachment orders are summed coherently and
// each local upper plus lower transfer equals the momentum of its produced
// pair
//
// [2,2], two simultaneous two-particle ladders:
// Separate columns meet only through the upper and lower multi-Reggeon proton
// vertices
//
//  ======F==========F======>
//        *          *
//        *          *
//        ff-----> a ff-----> c
//        |          |
//        ff-----> b ff-----> d
//        *          *
//        *          *
//  ======F==========F======>
//
// Charged-pion combinatorics with charge-balanced permutations:
// 8 canonical particle partitions * 2 attachment orders = 16 diagram terms
//
// The configured coherent amplitude is
//
// A_4 = A_[4] + A_[2,2]
//
// Total: 16 + 16 = 32 diagram terms before trajectory-channel sums
// Exact Cartesian orientation grouping reduces arithmetic without changing
// these amplitude-term counts
//
// For many similar amplitudes, see the literature of the era
// "pre-superstring/generalized Veneziano amplitude". For example:
//
// [REFERENCE: Bardakci, Ruegg, journals.aps.org/pr/pdf/10.1103/PhysRev.181.1884]
// [REFERENCE: Kycia, Lebiedowicz, Szczurek, Turnau, arxiv.org/abs/1702.07572]
//
std::complex<double> MRegge::EvalContinuum4(gra::LORENTZSCALAR &lts, const ReggeProductionModel model) { return EvalContinuumLadder(lts, model, 4); }

// ============================================================================
// Regge matrix element ansatz for 2->8 with six central particles
//
// [6], one serial six-particle ladder:
//
//  ======F======>
//        *
//        *
//        ff-----> a
//        |
//        ff-----> b
//        *
//        *
//        ff-----> c
//        |
//        ff-----> d
//        *
//        *
//        ff-----> e
//        |
//        ff-----> f
//        *
//        *
//  ======F======>
//
// Charged-pion combinatorics with charge-balanced permutations:
// 288 ordered particle diagrams
//
// The upper and lower two-exchange proton vertices are combined through the
// order two coherent multi Reggeon production convolution
//
// A_[4,2](Q_T) = c_2 / (8*pi^2*s)
//   * sum_partitions sum_(sigma in S_2) integral d^2 k
//   * <pp| B_(sigma(1))(k) B_(sigma(2))(Q_T-k) |pp>
//
// with c_2 = i/2. Here the two B operators are the complete serial
// four-particle and two-particle subladders
//
// [4,2], simultaneous four-particle and two-particle ladders:
// Separate columns meet only through the upper and lower multi-Reggeon proton
// vertices
//
//  ======F==========F======>
//        *          *
//        *          *
//        ff-----> a ff-----> e
//        |          |
//        ff-----> b ff-----> f
//        *          *
//        *          *
//        ff-----> c *
//        |          *
//        ff-----> d *
//        *          *
//        *          *
//  ======F==========F======>
//
// Charged-pion combinatorics with charge-balanced permutations:
// 288 particle partitions * 2 attachment orders = 576 diagram terms
//
// The upper and lower three-exchange proton vertices are combined through the
// order three coherent multi Reggeon production convolution
//
// A_[2,2,2](Q_T) = c_3 / (8*pi^2*s)^2
//   * sum_partitions sum_(sigma in S_3) integral d^2 k d^2 l
//   * <pp| B_(sigma(1))(k) B_(sigma(2))(l)
//          B_(sigma(3))(Q_T-k-l) |pp>
//
// with c_3 = -1/6. The implementation Fourier transforms
// every B_r to impact space, forms 3!*B_1(b)*B_2(b)*B_3(b), and transforms
// the product back at Q_T. This is the same convolution and attachment sum
// for the supported diagonal Good Walker operators
//
// [2,2,2], three simultaneous two-particle ladders:
//
//  ======F==========F==========F======>
//        *          *          *
//        *          *          *
//        ff-----> a ff-----> c ff-----> e
//        |          |          |
//        ff-----> b ff-----> d ff-----> f
//        *          *          *
//        *          *          *
//  ======F==========F==========F======>
//
// Charged-pion combinatorics with charge-balanced permutations:
// 48 canonical particle partitions * 6 attachment orders = 288 diagram terms
//
// The configured coherent amplitude is
//
// A_6 = A_[6] + A_[4,2] + A_[2,2,2]
// Total: 288 + 576 + 288 = 1152 diagram terms before trajectory-channel sums
// Exact Cartesian orientation grouping reduces arithmetic without changing
// these amplitude-term counts
//
std::complex<double> MRegge::EvalContinuum6(gra::LORENTZSCALAR &lts, const ReggeProductionModel model) { return EvalContinuumLadder(lts, model, 6); }

}  // namespace gra
