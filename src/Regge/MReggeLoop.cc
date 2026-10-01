// Parallel multi-Regge production convolution loop integrals
//
// Continuum4: Serial [4] and Parallel [2,2]
// Continuum6: Parallel [4,2] and Parallel [2,2,2]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeLoop.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <map>
#include <numeric>
#include <set>
#include <span>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeMulti.h"
#include "Graniitti/Regge/MReggeSource.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;
using gra::math::pow2;

namespace gra {
namespace {

// Route one parallel block transfer with exact four-momentum conservation
std::pair<M4Vec, M4Vec> ParallelBlockTransfers(const M4Vec &block, const double upper_pz, const double upper_energy, const double upper_x, const double upper_y) {
  const M4Vec upper(upper_x, upper_y, upper_pz, upper_energy);
  const M4Vec lower = block - upper;
  return {upper, lower};
}

// Split one complete permutation into the configured serial block sizes
std::vector<std::vector<int>> SplitParallelBlocks(const std::vector<int> &permutation, const regge::Topology &topology) {
  std::vector<std::vector<int>> blocks;
  blocks.reserve(topology.size());
  std::size_t offset = 0;
  for (const int size : topology) {
    blocks.emplace_back(permutation.cbegin() + offset, permutation.cbegin() + offset + size);
    offset += static_cast<std::size_t>(size);
  }
  return blocks;
}

// Compute a canonical key under exchange of equal-size parallel blocks
std::vector<int> ParallelPartitionKey(const std::vector<std::vector<int>> &blocks, const regge::Topology &topology) {
  std::vector<int> key;
  for (std::size_t begin = 0; begin < blocks.size();) {
    std::size_t end = begin + 1;
    while (end < blocks.size() && topology[end] == topology[begin]) { ++end; }
    std::vector<std::vector<int>> equal(blocks.cbegin() + begin, blocks.cbegin() + end);
    std::sort(equal.begin(), equal.end());
    for (const auto &block : equal) { key.insert(key.end(), block.cbegin(), block.cend()); }
    begin = end;
  }
  return key;
}

// Store one connected internal trajectory contribution inside a subladder
struct PreparedParallelChannel {
  int            upper_alias          = 0;
  int            upper_internal_alias = 0;
  int            lower_internal_alias = 0;
  int            lower_alias          = 0;
  SoftExchangeId upper_exchange{0};
  SoftExchangeId internal_exchange{0};
  SoftExchangeId lower_exchange{0};
  bool           has_internal_exchange = false;
  double         upper_sign            = 1.0;
  double         lower_sign            = 1.0;
  double         multiplicity          = 0.0;
  const regge::VertexParam *first_vertex  = nullptr;
  const regge::VertexParam *second_vertex = nullptr;
};

// Store the event-local data defining one serial production subladder
struct PreparedParallelBlock {
  std::vector<int>                     indices;
  std::vector<int>                     orientation_key;
  M4Vec                                momentum;
  M4Vec                                first_momentum;
  M4Vec                                first_pair_momentum;
  M4Vec                                third_momentum;
  const regge::PairParam              *first_entry  = nullptr;
  const regge::PairParam              *second_entry = nullptr;
  std::size_t                          first        = 0;
  std::size_t                          second       = 0;
  std::size_t                          third        = 0;
  std::size_t                          fourth       = 0;
  double                               first_mass2  = 0.0;
  double                               second_mass2 = 0.0;
  double                               gap_s        = 0.0;
  double                               upper_s      = 0.0;
  double                               lower_s      = 0.0;
  double                               upper_pz     = 0.0;
  double                               upper_energy = 0.0;
  std::vector<PreparedParallelChannel> channels;
};

// Store one exact final-state partition in terms of prepared subladders
struct RawParallelPartition {
  std::vector<std::size_t> block;
};

// Store independent orientation variants for each physical subladder
struct PreparedParallelPartition {
  std::vector<std::vector<std::size_t>> block;
};

// Store all unique subladders and final-state partitions for one topology
struct PreparedParallelBasis {
  std::vector<PreparedParallelBlock>     block;
  std::vector<PreparedParallelPartition> partition;
};

// Store diagonal Good Walker subladders in one contiguous attachment table
struct DiagonalParallelBank {
  DiagonalParallelBank(const std::size_t attachments, const std::size_t blocks, const std::size_t pair_dimension) : block_count(blocks), pair_dimension(pair_dimension), amplitude(attachments * blocks, pair_dimension, 0.0) {}

  // Compute one mutable subladder row without copying
  std::span<std::complex<double>> Block(const std::size_t attachment, const std::size_t block) { return amplitude.Row(attachment * block_count + block); }

  // Compute one immutable subladder row without copying
  std::span<const std::complex<double>> Block(const std::size_t attachment, const std::size_t block) const { return amplitude.Row(attachment * block_count + block); }

  std::size_t                   block_count    = 0;
  std::size_t                   pair_dimension = 0;
  MMatrix<std::complex<double>> amplitude;
};

// Identify one transverse attachment inside a diagonal subladder bank
struct ParallelBlockSlot {
  const DiagonalParallelBank *bank       = nullptr;
  std::size_t                 attachment = 0;

  // Compute one immutable subladder row from this attachment
  std::span<const std::complex<double>> Block(const std::size_t block) const { return bank->Block(attachment, block); }
};

// Reuse diagonal vertex work space across transverse attachment evaluations
struct ParallelResidueWorkspace {
  explicit ParallelResidueWorkspace(const std::size_t channels) : upper(channels, 0.0), lower(channels, 0.0) {}

  std::vector<double> upper;
  std::vector<double> lower;
};

// Build the structural key which differs only under pair orientation changes
std::vector<int> ParallelBlockOrientationKey(const std::vector<int> &block_indices) {
  std::vector<int> key = block_indices;
  for (std::size_t pair = 0; pair < key.size(); pair += 2) { std::sort(key.begin() + static_cast<std::ptrdiff_t>(pair), key.begin() + static_cast<std::ptrdiff_t>(pair + 2)); }
  return key;
}

// Build the connected trajectory channels of one prepared parallel subladder
void PrepareParallelChannels(PreparedParallelBlock &out, const LORENTZSCALAR &lts, const regge::Param &param, const SoftModel &soft, const ReggeProductionModel model) {
  const auto                                                                          soft_exchange  = [&](const int alias) { return param.exchanges.at(regge::TrajectoryIndex(param, alias)).soft_exchange; };
  const auto                                                                          is_pomeron     = [&](const int alias) { return soft.Exchange(soft_exchange(alias)).Role() == SoftExchangeRole::Pomeron; };
  const int                                                                           upper_beam_pdg = ReggeBeamParticleForLeg(lts, 1).pdg;
  const int                                                                           lower_beam_pdg = ReggeBeamParticleForLeg(lts, 2).pdg;
  std::map<std::tuple<int, int, int, int, bool, bool, bool>, PreparedParallelChannel> combined;

  for (const auto &first : out.first_entry->channels) {
    if (!regge::VertexAllowed(param, first) || !is_pomeron(first.first)) { continue; }
    if (out.second_entry == nullptr) {
      if (!is_pomeron(first.second)) { continue; }
      const double upper_sign = regge::AntiparticleSign(param, first.first, upper_beam_pdg);
      const double lower_sign = regge::AntiparticleSign(param, first.second, lower_beam_pdg);
      const auto   key        = std::make_tuple(first.first, 0, 0, first.second, false, upper_sign < 0.0, lower_sign < 0.0);
      auto [found, inserted]  = combined.try_emplace(key, PreparedParallelChannel{first.first, 0, 0, first.second, soft_exchange(first.first), SoftExchangeId{0}, soft_exchange(first.second), false, upper_sign, lower_sign, 0.0, &first, nullptr});
      (void)inserted;
      found->second.multiplicity += 1.0;
      continue;
    }

    const std::size_t internal = regge::TrajectoryIndex(param, first.second);
    for (const auto &second : out.second_entry->channels) {
      if (!regge::VertexAllowed(param, second)) { continue; }
      const bool connected = model == ReggeProductionModel::GP ? regge::TrajectoryIndex(param, second.first) == internal : second.first == first.second;
      if (!connected || !is_pomeron(second.second)) { continue; }
      const double upper_sign = regge::AntiparticleSign(param, first.first, upper_beam_pdg);
      const double lower_sign = regge::AntiparticleSign(param, second.second, lower_beam_pdg);
      const auto   key        = std::make_tuple(first.first, first.second, second.first, second.second, true, upper_sign < 0.0, lower_sign < 0.0);
      auto [found, inserted] =
          combined.try_emplace(key, PreparedParallelChannel{first.first, first.second, second.first, second.second, soft_exchange(first.first), soft_exchange(first.second), soft_exchange(second.second), true, upper_sign, lower_sign, 0.0, &first, &second});
      (void)inserted;
      found->second.multiplicity += 1.0;
    }
  }

  out.channels.reserve(combined.size());
  for (const auto &[key, channel] : combined) {
    (void)key;
    if (channel.multiplicity > 0.0) { out.channels.push_back(channel); }
  }
}

// Prepare one connected Pomeron-to-Pomeron serial subladder
PreparedParallelBlock PrepareParallelBlock(const gra::LORENTZSCALAR &lts, const regge::Param &param, const gra::SoftModel &soft, const ReggeProductionModel model, const std::vector<int> &block_indices) {
  if (block_indices.size() != 2 && block_indices.size() != 4) { throw std::invalid_argument("MRegge parallel topology supports two- and four-body blocks"); }

  PreparedParallelBlock out;
  out.indices         = block_indices;
  out.orientation_key = ParallelBlockOrientationKey(block_indices);
  for (const int amplitude_index : block_indices) { out.momentum += lts.decaytree[regge::DecayIndex(amplitude_index, lts, "MRegge parallel topology")].p4; }
  const M4Vec  total       = lts.q1 + lts.q2;
  const double total_plus  = total.LightconePos();
  const double total_minus = total.LightconeNeg();
  if (!(total_plus > 0.0) || !(total_minus > 0.0) || !std::isfinite(total_plus) || !std::isfinite(total_minus)) { throw AmplitudeFailure("MRegge parallel topology has invalid light-cone transfer"); }
  const double upper_plus  = out.momentum.LightconePos() * lts.q1.LightconePos() / total_plus;
  const double upper_minus = out.momentum.LightconeNeg() * lts.q1.LightconeNeg() / total_minus;
  out.upper_pz             = 0.5 * (upper_plus - upper_minus);
  out.upper_energy         = 0.5 * (upper_plus + upper_minus);
  out.first                = regge::DecayIndex(block_indices[0], lts, "MRegge parallel topology");
  out.second               = regge::DecayIndex(block_indices[1], lts, "MRegge parallel topology");
  out.first_entry          = regge::FindPair(param, {lts.decaytree[out.first].p.pdg, lts.decaytree[out.second].p.pdg}, model);
  if (out.first_entry == nullptr) { return out; }
  out.first_mass2         = pow2(lts.decaytree[out.first].p.mass);
  out.first_momentum      = lts.decaytree[out.first].p4;
  out.first_pair_momentum = out.first_momentum + lts.decaytree[out.second].p4;
  if (block_indices.size() == 4) {
    out.third        = regge::DecayIndex(block_indices[2], lts, "MRegge parallel topology");
    out.fourth       = regge::DecayIndex(block_indices[3], lts, "MRegge parallel topology");
    out.second_entry = regge::FindPair(param, {lts.decaytree[out.third].p.pdg, lts.decaytree[out.fourth].p.pdg}, model);
    if (out.second_entry == nullptr) { return out; }
    out.second_mass2   = pow2(lts.decaytree[out.third].p.mass);
    out.third_momentum = lts.decaytree[out.third].p4;
    out.gap_s          = (lts.decaytree[out.second].p4 + lts.decaytree[out.third].p4).M2();
  }

  PrepareParallelChannels(out, lts, param, soft, model);

  out.upper_s                      = (lts.pfinal[1] + lts.decaytree[out.first].p4).M2();
  const std::size_t lower_particle = (out.second_entry == nullptr) ? out.second : out.fourth;
  out.lower_s                      = (lts.pfinal[2] + lts.decaytree[lower_particle].p4).M2();
  if (!std::isfinite(out.upper_s) || !(out.upper_s > 0.0) || !std::isfinite(out.lower_s) || !(out.lower_s > 0.0) || (out.second_entry != nullptr && (!std::isfinite(out.gap_s) || !(out.gap_s > 0.0)))) {
    throw AmplitudeFailure("MRegge parallel topology has invalid subenergy");
  }

  return out;
}

// Contract complete Cartesian pair orientations before attachment products
std::vector<PreparedParallelPartition> ContractParallelOrientations(const std::vector<PreparedParallelBlock> &blocks, const regge::Topology &topology, const std::vector<RawParallelPartition> &raw_partitions) {
  using StructureKey = std::vector<std::vector<int>>;
  std::map<StructureKey, std::set<std::vector<std::size_t>>> groups;
  for (const auto &raw : raw_partitions) {
    if (raw.block.size() != topology.size()) { throw std::logic_error("MRegge parallel orientation partition dimension disagrees"); }
    std::vector<std::size_t> ordered = raw.block;
    for (std::size_t begin = 0; begin < topology.size();) {
      std::size_t end = begin + 1;
      while (end < topology.size() && topology[end] == topology[begin]) { ++end; }
      std::sort(ordered.begin() + static_cast<std::ptrdiff_t>(begin), ordered.begin() + static_cast<std::ptrdiff_t>(end),
                [&](const std::size_t first, const std::size_t second) { return std::tie(blocks[first].orientation_key, first) < std::tie(blocks[second].orientation_key, second); });
      begin = end;
    }
    StructureKey structure;
    structure.reserve(ordered.size());
    for (const std::size_t block : ordered) { structure.push_back(blocks.at(block).orientation_key); }
    groups[structure].insert(std::move(ordered));
  }

  std::vector<PreparedParallelPartition> output;
  for (const auto &[structure, assignments] : groups) {
    (void)structure;
    PreparedParallelPartition contracted;
    contracted.block.resize(topology.size());
    for (const auto &assignment : assignments) {
      for (const auto &position : indices(assignment)) { contracted.block[position].push_back(assignment[position]); }
    }
    std::size_t cartesian_size = 1;
    for (auto &variants : contracted.block) {
      std::sort(variants.begin(), variants.end());
      variants.erase(std::unique(variants.begin(), variants.end()), variants.end());
      if (variants.empty() || cartesian_size > std::numeric_limits<std::size_t>::max() / variants.size()) { throw std::logic_error("MRegge parallel orientation contraction overflow"); }
      cartesian_size *= variants.size();
    }
    if (cartesian_size == assignments.size()) {
      output.push_back(std::move(contracted));
      continue;
    }
    for (const auto &assignment : assignments) {
      PreparedParallelPartition exact;
      exact.block.resize(assignment.size());
      for (const auto &position : indices(assignment)) { exact.block[position].push_back(assignment[position]); }
      output.push_back(std::move(exact));
    }
  }
  return output;
}

// Build every unique subladder and canonical particle partition
PreparedParallelBasis PrepareParallelBasis(const gra::LORENTZSCALAR &lts, const regge::Param &param, const gra::SoftModel &soft, const ReggeProductionModel model, const regge::Topology &topology,
                                           const std::vector<std::vector<int>> &permutations) {
  PreparedParallelBasis                   out;
  std::vector<RawParallelPartition>       raw_partitions;
  std::set<std::vector<int>>              seen;
  std::map<std::vector<int>, std::size_t> block_id;
  for (const auto &permutation : permutations) {
    const std::size_t expected = std::accumulate(topology.cbegin(), topology.cend(), std::size_t{0});
    if (permutation.size() != expected) { throw std::invalid_argument("MRegge parallel topology and permutation sizes disagree"); }
    const auto blocks = SplitParallelBlocks(permutation, topology);
    if (!seen.insert(ParallelPartitionKey(blocks, topology)).second) { continue; }
    RawParallelPartition partition;
    for (const auto &block : blocks) {
      auto found = block_id.find(block);
      if (found == block_id.end()) {
        auto prepared = PrepareParallelBlock(lts, param, soft, model, block);
        if (prepared.channels.empty()) {
          partition.block.clear();
          break;
        }
        const std::size_t id = out.block.size();
        out.block.push_back(std::move(prepared));
        found = block_id.emplace(block, id).first;
      }
      partition.block.push_back(found->second);
    }
    if (!partition.block.empty()) { raw_partitions.push_back(std::move(partition)); }
  }
  out.partition = ContractParallelOrientations(out.block, topology, raw_partitions);
  return out;
}

// Contract one parallel subladder through the same local spin kernels
std::complex<double> ParallelBlockSpinFactor(const gra::LORENTZSCALAR &lts, const regge::Param &param, const ReggeProductionModel model, const PreparedParallelBlock &block, const PreparedParallelChannel &channel,
                                             const std::pair<M4Vec, M4Vec> &transfer, const std::complex<double> first_meson, const std::complex<double> second_meson, const std::complex<double> internal_propagator) {
  if (channel.first_vertex == nullptr || (channel.has_internal_exchange && channel.second_vertex == nullptr)) { throw std::logic_error("MRegge parallel continuum channel has no resolved vertex parameters"); }
  const ReggeContinuumPoleKey first_key      = {lts.decaytree[block.first].p.pdg, lts.decaytree[block.second].p.pdg, channel.upper_alias, channel.has_internal_exchange ? channel.upper_internal_alias : channel.lower_alias};
  const auto                 *first_cache    = &ladder::Pole(lts, first_key, "MRegge::ParallelBlockSpinFactor first");
  const M4Vec                 lower_outgoing = -transfer.second;
  const M4Vec                 internal       = transfer.first - block.first_pair_momentum;
  const double first_meson_q2 = (transfer.first - block.first_momentum).M2();
  const double first_vertex_factor = regge::MesonVertexFactor(first_meson_q2, block.first_mass2, *channel.first_vertex) *
                                     regge::PairVertexFactor(param, *channel.first_vertex, lts.decaytree[block.first], lts.decaytree[block.second], transfer.first.M2(), channel.has_internal_exchange ? internal.M2() : lower_outgoing.M2(),
                                                             param.con.multiregge_transfer_ext, channel.has_internal_exchange ? param.con.multiregge_transfer_int : param.con.multiregge_transfer_ext);
  auto first_kernel = ladder::Kernel(
      *first_cache, param, model, lts.decaytree[block.first], lts.decaytree[block.second], transfer.first,
      channel.has_internal_exchange ? internal : lower_outgoing,
      regge::Alpha(param, channel.upper_alias, transfer.first.M2()),
      regge::Alpha(param, channel.has_internal_exchange ? channel.upper_internal_alias : channel.lower_alias,
                   channel.has_internal_exchange ? internal.M2() : lower_outgoing.M2()));
  first_kernel *= first_meson * first_vertex_factor;

  MMatrix<std::complex<double>> second_kernel;
  MMatrix<std::complex<double>> metric;
  const ReggeContinuumPole     *last_cache = first_cache;
  if (channel.has_internal_exchange) {
    const ReggeContinuumPoleKey second_key   = {lts.decaytree[block.third].p.pdg, lts.decaytree[block.fourth].p.pdg, channel.lower_internal_alias, channel.lower_alias};
    const auto                 *second_cache = &ladder::Pole(lts, second_key, "MRegge::ParallelBlockSpinFactor second");
    const auto                  first_basis  = ladder::Basis(*first_cache, model, 1);
    const auto                  second_basis = ladder::Basis(*second_cache, model, 0);
    ladder::CheckBasis(first_basis, second_basis, "MRegge::ParallelBlockSpinFactor");
    metric = spin::ExchangeMetric(ExchangeBasisType::ReducedRegge, first_basis, 0.0);
    metric *= internal_propagator;
    second_kernel =
        ladder::Kernel(*second_cache, param, model, lts.decaytree[block.third], lts.decaytree[block.fourth], internal,
                       lower_outgoing, regge::Alpha(param, channel.lower_internal_alias, internal.M2()),
                       regge::Alpha(param, channel.lower_alias, lower_outgoing.M2()));
    const double second_meson_q2 = (internal - block.third_momentum).M2();
    const double second_vertex_factor = regge::MesonVertexFactor(second_meson_q2, block.second_mass2, *channel.second_vertex) *
                                        regge::PairVertexFactor(param, *channel.second_vertex, lts.decaytree[block.third], lts.decaytree[block.fourth], internal.M2(), lower_outgoing.M2(), param.con.multiregge_transfer_int,
                                                                param.con.multiregge_transfer_ext);
    second_kernel *= second_meson * second_vertex_factor;
    last_cache = second_cache;
  }

  const spin::ForwardSpec forward     = {lts.process.FORWARD_VERTEX, param.s0, ExchangeBasisType::ReducedRegge};
  const auto              upper_basis = ladder::Basis(*first_cache, model, 0);
  const auto              lower_basis = ladder::Basis(*last_cache, model, 1);
  const auto              upper       = spin::ExchangeBoundary(upper_basis, transfer.first, false, forward);
  const auto              lower       = spin::ExchangeBoundary(lower_basis, transfer.second, true, forward);
  if (!channel.has_internal_exchange) {
    const std::array<MMatrix<std::complex<double>>, 1> kernels = {std::move(first_kernel)};
    return channel.multiplicity * ladder::Contract(upper, kernels, {}, lower);
  }
  const std::array<MMatrix<std::complex<double>>, 2> kernels = {std::move(first_kernel), std::move(second_kernel)};
  const std::array<MMatrix<std::complex<double>>, 1> metrics = {std::move(metric)};
  return channel.multiplicity * ladder::Contract(upper, kernels, metrics, lower);
}

// Build one diagonal Pomeron subladder at a local transverse transfer
template <typename Propagator>
void EvaluateParallelBlock(const gra::LORENTZSCALAR &lts, const gra::MRegge &regge, const regge::Param &param, const Propagator &propagator, const ReggeProductionModel model, const PreparedParallelBlock &block, const double upper_x,
                           const double upper_y, std::span<std::complex<double>> output, std::span<double> upper_residue, std::span<double> lower_residue) {
  const auto       &soft           = *regge.SoftModelHandle();
  const std::size_t channel_count  = soft.GoodWalker().ChannelCount();
  const std::size_t pair_dimension = soft.GoodWalker().PairDimension();
  if (output.size() != pair_dimension) { throw std::logic_error("MRegge parallel topology has an invalid output dimension"); }
  if (upper_residue.size() != channel_count || lower_residue.size() != channel_count) { throw std::logic_error("MRegge parallel topology has an invalid residue buffer"); }
  const auto                 transfer             = ParallelBlockTransfers(block.momentum, block.upper_pz, block.upper_energy, upper_x, upper_y);
  const double               upper_t              = transfer.first.M2();
  const double               lower_t              = transfer.second.M2();
  const M4Vec                first_meson_transfer = transfer.first - block.first_momentum;
  const std::complex<double> first_meson          = regge::MesonPropagator(param, first_meson_transfer.M2(), block.first_mass2, *block.first_entry, lts.decaytree[block.first], lts.decaytree[block.second]);
  const M4Vec                internal             = transfer.first - block.first_pair_momentum;
  std::complex<double>       second_meson         = 1.0;
  if (block.second_entry != nullptr) {
    const M4Vec second_meson_transfer = internal - block.third_momentum;
    second_meson                      = regge::MesonPropagator(param, second_meson_transfer.M2(), block.second_mass2, *block.second_entry, lts.decaytree[block.third], lts.decaytree[block.fourth]);
  }

  std::fill(output.begin(), output.end(), 0.0);
  for (const auto &channel : block.channels) {
    const std::complex<double> internal_propagator = !channel.has_internal_exchange ? std::complex<double>(1.0, 0.0) : propagator(block.gap_s, internal.M2(), channel.internal_exchange);
    const std::complex<double> central             = ParallelBlockSpinFactor(lts, param, model, block, channel, transfer, first_meson, second_meson, internal_propagator);
    const std::complex<double> kernel              = propagator(block.upper_s, upper_t, channel.upper_exchange) * propagator(block.lower_s, lower_t, channel.lower_exchange) * central * channel.upper_sign * channel.lower_sign;
    soft.DiagonalResidues(channel.upper_exchange, upper_t, upper_residue);
    soft.DiagonalResidues(channel.lower_exchange, lower_t, lower_residue);
    for (std::size_t upper = 0; upper < channel_count; ++upper) {
      for (std::size_t lower = 0; lower < channel_count; ++lower) {
        const std::size_t pair = upper * channel_count + lower;
        output[pair] += kernel * upper_residue[upper] * lower_residue[lower];
      }
    }
  }
}

// Evaluate every unique subladder at one transverse attachment
template <typename Propagator>
void EvaluateParallelBlockBank(const gra::LORENTZSCALAR &lts, const gra::MRegge &regge, const regge::Param &param, const Propagator &propagator, const ReggeProductionModel model, const std::vector<PreparedParallelBlock> &blocks,
                               const double x, const double y, DiagonalParallelBank &out, const std::size_t attachment, ParallelResidueWorkspace &workspace) {
  if (out.block_count != blocks.size() || out.pair_dimension != regge.SoftModelHandle()->GoodWalker().PairDimension() || attachment * out.block_count + out.block_count > out.amplitude.size_row()) {
    throw std::logic_error("MRegge parallel topology has an invalid subladder bank");
  }
  for (const auto &block : indices(blocks)) {
    auto output = out.Block(attachment, block);
    // Center each transverse transfer between the beams so the finite cutoff respects beam exchange
    EvaluateParallelBlock(lts, regge, param, propagator, model, blocks[block],
                          x + 0.5 * blocks[block].momentum.Px(), y + 0.5 * blocks[block].momentum.Py(),
                          output, workspace.upper, workspace.lower);
  }
}

// Add the symmetrized production vertex attachment for one particle partition
void AccumulateParallelInsertion(std::vector<std::complex<double>> &target, const std::vector<double> &proton, const std::span<const ParallelBlockSlot> slot, const PreparedParallelPartition &partition, const std::complex<double> scale) {
  if (slot.size() != partition.block.size() || (slot.size() != 2 && slot.size() != 3)) { throw std::logic_error("MRegge parallel insertion has an invalid block count"); }
  const auto block_sum = [&](const std::size_t attachment, const std::size_t position, const std::size_t pair) {
    std::complex<double> value = 0.0;
    for (const std::size_t block : partition.block[position]) { value += slot[attachment].Block(block)[pair]; }
    return value;
  };
  if (slot.size() == 2) {
    for (const auto &i : indices(target)) {
      const std::complex<double> a0 = block_sum(0, 0, i);
      const std::complex<double> b0 = block_sum(0, 1, i);
      const std::complex<double> a1 = block_sum(1, 0, i);
      const std::complex<double> b1 = block_sum(1, 1, i);
      target[i] += scale * proton[i] * (a0 * b1 + b0 * a1);
    }
    return;
  }

  for (const auto &i : indices(target)) {
    const std::complex<double> a0        = block_sum(0, 0, i);
    const std::complex<double> b0        = block_sum(0, 1, i);
    const std::complex<double> c0        = block_sum(0, 2, i);
    const std::complex<double> a1        = block_sum(1, 0, i);
    const std::complex<double> b1        = block_sum(1, 1, i);
    const std::complex<double> c1        = block_sum(1, 2, i);
    const std::complex<double> a2        = block_sum(2, 0, i);
    const std::complex<double> b2        = block_sum(2, 1, i);
    const std::complex<double> c2        = block_sum(2, 2, i);
    const std::complex<double> insertion = a0 * b1 * c2 + a0 * c1 * b2 + b0 * a1 * c2 + b0 * c1 * a2 + c0 * a1 * b2 + c0 * b1 * a2;
    target[i] += scale * proton[i] * insertion;
  }
}

// Store unique orientation sums and one group index for every partition slot
struct ParallelOrientationGroups {
  std::vector<std::vector<std::size_t>>   block;
  std::vector<std::array<std::size_t, 3>> partition;
};

// Intern the exact orientation sums used by the three ladder partitions
ParallelOrientationGroups PrepareParallelOrientationGroups(const PreparedParallelBasis &basis) {
  ParallelOrientationGroups                       output;
  std::map<std::vector<std::size_t>, std::size_t> group_index;
  output.partition.reserve(basis.partition.size());
  for (const auto &partition : basis.partition) {
    if (partition.block.size() != 3) { throw std::logic_error("MRegge parallel Fourier topology is not three-body"); }
    std::array<std::size_t, 3> mapped{};
    for (const auto &position : indices(mapped)) {
      std::vector<std::size_t> key = partition.block[position];
      if (key.empty()) { throw std::logic_error("MRegge parallel Fourier orientation bank is empty"); }
      std::sort(key.begin(), key.end());
      const auto [found, inserted] = group_index.try_emplace(key, output.block.size());
      if (inserted) { output.block.push_back(std::move(key)); }
      mapped[position] = found->second;
    }
    output.partition.push_back(mapped);
  }
  return output;
}

// Transform every distinct orientation sum into impact parameter space
std::vector<std::vector<math::PolarHarmonicField>> TransformParallelOrientationBank(const DiagonalParallelBank &bank, const ParallelOrientationGroups &orientation, const math::MPolarFourier &transform) {
  const std::size_t radial_count  = transform.MomentumRule().node.size();
  const std::size_t azimuth_count = transform.AzimuthNodes().size();
  if (bank.amplitude.size_row() != radial_count * azimuth_count * bank.block_count) { throw std::logic_error("MRegge parallel Fourier bank dimensions disagree"); }

  std::vector<std::vector<math::PolarHarmonicField>> output(orientation.block.size(), std::vector<math::PolarHarmonicField>(bank.pair_dimension));
  MMatrix<std::complex<double>>                      momentum_field(radial_count, azimuth_count, 0.0);
  for (const auto &group : indices(orientation.block)) {
    for (std::size_t pair = 0; pair < bank.pair_dimension; ++pair) {
      for (const auto &radial : indices(transform.MomentumRule().node)) {
        for (const auto &azimuth : indices(transform.AzimuthNodes())) {
          const std::size_t    attachment = radial * azimuth_count + azimuth;
          std::complex<double> sum        = 0.0;
          for (const std::size_t block : orientation.block[group]) { sum += bank.Block(attachment, block)[pair]; }
          momentum_field(radial, azimuth) = sum;
        }
      }
      output[group][pair] = transform.Forward(momentum_field);
    }
  }
  return output;
}

// Select the partition position with the fewest distinct orientation groups
std::size_t ParallelFourierAnchor(const ParallelOrientationGroups &orientation) {
  std::size_t best_position = 0;
  std::size_t best_count    = std::numeric_limits<std::size_t>::max();
  for (std::size_t position = 0; position < 3; ++position) {
    std::set<std::size_t> unique;
    for (const auto &partition : orientation.partition) { unique.insert(partition[position]); }
    if (unique.size() < best_count) {
      best_position = position;
      best_count    = unique.size();
    }
  }
  return best_position;
}

// Evaluate the three-subladder convolution as a product in impact space
void AccumulateParallelFourierTriple(std::vector<std::complex<double>> &target, const std::vector<double> &proton, const PreparedParallelBasis &basis, const DiagonalParallelBank &momentum_bank, const math::MPolarFourier &transform,
                                     const double qx, const double qy, const std::complex<double> normalization) {
  if (target.size() != proton.size() || target.size() != momentum_bank.pair_dimension) { throw std::logic_error("MRegge parallel Fourier pair dimensions disagree"); }
  const auto       orientation        = PrepareParallelOrientationGroups(basis);
  const auto       field              = TransformParallelOrientationBank(momentum_bank, orientation, transform);
  const double     q                  = std::hypot(qx, qy);
  const double     q_azimuth          = std::fpclassify(q) == FP_ZERO ? 0.0 : std::atan2(qy, qx);
  constexpr double linked_assignments = 6.0;

  std::vector<std::size_t>              active_pair;
  std::vector<math::PolarHarmonicField> topology_field;
  active_pair.reserve(target.size());
  topology_field.reserve(target.size());
  const std::size_t anchor       = ParallelFourierAnchor(orientation);
  const std::size_t first_other  = (anchor + 1) % 3;
  const std::size_t second_other = (anchor + 2) % 3;

  for (const auto &pair : indices(target)) {
    math::PolarHarmonicField                        pair_field;
    bool                                            initialized = false;
    std::map<std::size_t, math::PolarHarmonicField> inner_sum;
    for (const auto &partition : orientation.partition) {
      auto product                 = transform.Multiply(field[partition[first_other]][pair], field[partition[second_other]][pair]);
      const auto [found, inserted] = inner_sum.try_emplace(partition[anchor], std::move(product));
      if (!inserted) { found->second.coefficient += product.coefficient; }
    }
    for (auto &[anchor_group, inner] : inner_sum) {
      auto product = transform.Multiply(field[anchor_group][pair], inner);
      if (!initialized) {
        pair_field  = std::move(product);
        initialized = true;
      } else {
        pair_field.coefficient += product.coefficient;
      }
    }
    if (!initialized) { continue; }
    active_pair.push_back(pair);
    topology_field.push_back(std::move(pair_field));
  }
  if (topology_field.empty()) { return; }

  const auto inverse_kernel = transform.PrepareInverseKernel(q, q_azimuth, topology_field.front().max_harmonic);
  const auto inverse        = transform.InverseEvaluateBatch(topology_field, inverse_kernel);
  for (const auto &i : indices(active_pair)) {
    const std::size_t          pair        = active_pair[i];
    const std::complex<double> convolution = 4.0 * math::PIPI * inverse[i];
    target[pair] += linked_assignments * normalization * proton[pair] * convolution;
  }
}

// Store one completed parallel topology in the proton Good Walker source
void StoreParallelSource(LORENTZSCALAR &lts, const SoftModelPtr &soft, const std::vector<std::complex<double>> &pair_source) {
  const std::size_t rows = lts.hamp.metadata.spin_rows;
  if (rows != 1 && rows != 4 && rows != 16) { throw std::invalid_argument("MRegge parallel topology has an invalid proton spin row count"); }
  MMatrix<std::complex<double>> source(rows, pair_source.size(), 0.0);
  if (rows == 16) {
    for (std::size_t initial = 0; initial < 4; ++initial) {
      for (const auto &col : indices(pair_source)) { source(gra::spin::PairHelicityTransitionIndex(initial, initial), col) = pair_source[col]; }
    }
  } else {
    for (std::size_t row = 0; row < rows; ++row) {
      for (const auto &col : indices(pair_source)) { source(row, col) = pair_source[col]; }
    }
  }

  PrepareGoodWalker(lts, soft);
  auto      &components = lts.proton_good_walker->components;
  const auto match      = std::find_if(components.begin(), components.end(), [&](const auto &item) {
    return item.coherence_group == 0 && item.upper_sector == ProtonGoodWalkerSector::Elastic && item.lower_sector == ProtonGoodWalkerSector::Elastic;
  });
  if (match == components.end()) {
    components.push_back({0, ProtonGoodWalkerSector::Elastic, ProtonGoodWalkerSector::Elastic, std::move(source)});
  } else {
    match->source += source;
  }
}

}  // namespace

// For T=[n1,...,nm], accumulate the coherent multi Reggeon production amplitude
//
// A_T = c_m / (8*pi^2*s)^(m-1)
//       * sum_partitions sum_sigma integral prod_(r=1)^(m-1) d^2 k_r
//       * <pp| prod_(r=1)^m B_(sigma(r))(k_r + P_(sigma(r))T/2) |pp>
//
// Every model-supported final-state permutation is split into the configured
// blocks. Exchanges of equal-size blocks are stored once in sum_partitions,
// while sum_sigma restores all m! transverse attachment assignments without
// double counting them
//
// where Q_T = (q1T - q2T)/2 and k_m = Q_T - sum_(r=1)^(m-1) k_r
// The order-two convolution is integrated directly. For order three, each
// shifted B_r is Fourier transformed and their product is inverted at Q_T
// The factor 3! restores the linked attachment sum already divided by 3! in c_3
// Accumulate one bare coherent multi Reggeon production sector
void MRegge::AccumulateParallelTopology(gra::LORENTZSCALAR &lts, const ReggeProductionModel model, const regge::Topology &topology, const std::vector<std::vector<int>> &permutations) const {
  const auto basis = PrepareParallelBasis(lts, param, *soft_model, model, topology, permutations);
  if (basis.partition.empty()) { return; }
  const auto &rule = regge_numerics->ParallelRule();
  const std::size_t azimuth_count    = rule.kt_x.size_col();
  const std::size_t attachment_count = rule.kt.size() * azimuth_count;
  const std::size_t d                = soft_model->GoodWalker().PairDimension();
  const auto        proton           = gra::KroneckerProduct(soft_model->GoodWalker().ProtonVector(), soft_model->GoodWalker().ProtonVector());
  if (proton.size() != d) { throw std::logic_error("MRegge parallel topology has an invalid proton pair source"); }

  // Average the labelled attachments and include one longitudinal loop phase per sewing
  std::complex<double> normalization = 1.0 / math::factorial(topology.size());
  for (std::size_t block = 1; block < topology.size(); ++block) { normalization *= math::zi / (8.0 * math::PIPI * lts.s); }
  DiagonalParallelBank     free_bank(attachment_count, basis.block.size(), d);
  ParallelResidueWorkspace residue_workspace(soft_model->GoodWalker().ChannelCount());
  // Keep resolved propagator access inside the owning amplitude class
  const auto prepared_propagator = [this](const double subenergy, const double transfer, const SoftExchangeId exchange) { return PreparedPropagator(subenergy, transfer, exchange); };
  for (const auto &radial : indices(rule.kt)) {
    for (std::size_t azimuth = 0; azimuth < azimuth_count; ++azimuth) {
      const std::size_t attachment = radial * azimuth_count + azimuth;
      EvaluateParallelBlockBank(lts, *this, param, prepared_propagator, model, basis.block, rule.kt_x[radial][azimuth], rule.kt_y[radial][azimuth], free_bank, attachment, residue_workspace);
    }
  }

  // The shifted loop momenta sum to (q1T - q2T)/2
  const M4Vec transfer = (lts.q1 - lts.q2) * 0.5;
  std::vector<std::complex<double>> pair_source(d, 0.0);
  if (topology.size() == 2) {
    DiagonalParallelBank residual_bank(1, basis.block.size(), d);
    for (const auto &radial : indices(rule.kt)) {
      for (std::size_t azimuth = 0; azimuth < azimuth_count; ++azimuth) {
        const std::size_t attachment = radial * azimuth_count + azimuth;
        EvaluateParallelBlockBank(lts, *this, param, prepared_propagator, model, basis.block, transfer.Px() - rule.kt_x[radial][azimuth], transfer.Py() - rule.kt_y[radial][azimuth], residual_bank, 0, residue_workspace);
        const std::array<ParallelBlockSlot, 2> slot = {ParallelBlockSlot{&free_bank, attachment}, {&residual_bank, 0}};
        for (const auto &partition : basis.partition) { AccumulateParallelInsertion(pair_source, proton, slot, partition, normalization * rule.measure_weight[radial][azimuth]); }
      }
    }
  } else {
    if (parallel_fourier == nullptr) { throw std::logic_error("MRegge parallel Fourier transform is not initialized"); }
    AccumulateParallelFourierTriple(pair_source, proton, basis, free_bank, *parallel_fourier, transfer.Px(), transfer.Py(), normalization);
  }

  StoreParallelSource(lts, soft_model, pair_source);
}

}  // namespace gra
