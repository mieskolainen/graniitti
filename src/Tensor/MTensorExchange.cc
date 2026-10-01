// Tensor Pomeron and tensor Reggeon exchange steering
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorExchange.h"

// C++
#include <algorithm>
#include <cmath>
#include <set>
#include <stdexcept>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra {
namespace {

using gra::aux::indices;

// Parse one integer sign without narrowing oversized JSON values
int ParseUnitSign(const nlohmann::json &value, const std::string &context) {
  if (!value.is_number_integer()) {
    throw std::invalid_argument(context + " must be -1 or 1");
  }
  if (value.is_number_unsigned()) {
    if (value.get<nlohmann::json::number_unsigned_t>() == 1U) {
      return 1;
    }
  } else {
    const auto sign = value.get<nlohmann::json::number_integer_t>();
    if (sign == -1 || sign == 1) {
      return static_cast<int>(sign);
    }
  }
  throw std::invalid_argument(context + " must be -1 or 1");
}

// Reject unknown object fields in strict steering blocks
void ValidateFields(const nlohmann::json &block,
                    const std::set<std::string> &allowed,
                    const std::string &context) {
  if (!block.is_object()) {
    throw std::invalid_argument(context + " must be an object");
  }
  for (const auto &[field, value] : block.items()) {
    (void)value;
    if (!allowed.contains(field)) {
      throw std::invalid_argument(context + " has unknown field " + field);
    }
  }
}

// Parse one decimal PDG object key without accepting trailing text
int ParsePdgKey(const std::string &key, const std::string &context) {
  std::size_t end = 0;
  int pdg = 0;
  try {
    pdg = std::stoi(key, &end);
  } catch (const std::exception &) {
    throw std::invalid_argument(context + " has invalid PDG key " + key);
  }
  if (end != key.size() || pdg == 0) {
    throw std::invalid_argument(context + " has invalid PDG key " + key);
  }
  return pdg;
}

// Parse one compact two-PDG JSON object key
std::array<int, 2> ParsePairKey(const std::string &key,
                                const std::string &context) {
  nlohmann::json pair;
  try {
    pair = nlohmann::json::parse(key);
  } catch (const std::exception &) {
    throw std::invalid_argument(context + " has invalid pair key " + key);
  }
  if (!pair.is_array() || pair.size() != 2 || !pair[0].is_number_integer() ||
      !pair[1].is_number_integer()) {
    throw std::invalid_argument(context +
                                " pair key must contain two integers");
  }
  return {pair[0].get<int>(), pair[1].get<int>()};
}

// Compute an order-independent key for one signed PDG pair
std::array<int, 2> CanonicalPair(std::array<int, 2> pair) {
  if (pair[1] < pair[0]) {
    std::swap(pair[0], pair[1]);
  }
  return pair;
}

// Parse one published rank-two exchange normalization identifier
TensorExchangeType ParseExchangeType(const std::string &type,
                                     const std::string &context) {
  if (type == "P") {
    return TensorExchangeType::Pomeron;
  }
  if (type == "R2") {
    return TensorExchangeType::TensorReggeon;
  }
  if (type == "O") {
    return TensorExchangeType::Odderon;
  }
  if (type == "R1") {
    return TensorExchangeType::VectorReggeon;
  }
  throw std::invalid_argument(context + " type must be P, R2, O or R1");
}

// Compute the effective Lorentz rank selected by one exchange type
int ExchangeRank(const TensorExchangeType type) {
  return type == TensorExchangeType::Pomeron ||
                 type == TensorExchangeType::TensorReggeon
             ? 2
             : 1;
}

// Compute the required number of covariant continuum couplings by hadron spin
std::size_t CouplingCount(const int spinX2, const std::string &context) {
  if (spinX2 == 0 || spinX2 == 1) {
    return 1;
  }
  if (spinX2 == 2) {
    return 2;
  }
  throw std::invalid_argument(context + " supports hadron spinX2 = 0, 1 or 2");
}

// Validate the implemented currents and the physical crossed vertex before pruning
void ValidateVertex(const TensorExchangeParam &exchange, const MParticle &mother,
                    const MParticle &first, const MParticle &second,
                    const std::string &context) {
  if (exchange.rank == 1 && first.spinX2 == 2) {
    throw std::invalid_argument(context + " vector pairs require a rank-two exchange");
  }
  if (exchange.type == TensorExchangeType::Odderon && first.spinX2 == 0) {
    throw std::invalid_argument(context + " has no implemented Odderon pseudoscalar current");
  }
  if (mother.chargeX3 != first.chargeX3 + second.chargeX3) {
    throw std::invalid_argument(context + " violates charge conservation");
  }
  // [REFERENCE: Ewerz, Maniatis and Nachtmann, arXiv:1309.3478, Table 1 and Section 2]
  if (mother.G != 0 && first.G != 0 && second.G != 0 && mother.G != first.G * second.G) {
    throw std::invalid_argument(context + " violates G parity");
  }
  HELMatrix helicity;
  helicity.C_symmetry = true;
  helicity.P_symmetry = true;
  if (spin::AllowedLSCouplings(helicity, mother, first, second, false, context, spin::VertexContext::Auto).empty()) {
    throw std::invalid_argument(context + " has no allowed JPC and Bose or Fermi coupling");
  }
}

} // namespace

// Read Tensor exchange and continuum steering from one immutable tune
void MTensorExchangeModel::Configure(const nlohmann::json &general,
                                     const nlohmann::json &vertex_table,
                                     const gra::MPDG &pdg_table,
                                     const double coupling_min) {
  exchanges.clear();
  vertices.clear();
  continuum.clear();
  hadron_spins.clear();
  initialized = false;

  const auto &exchange_table = general.at("PARAM_TENSORPOM").at("EXCHANGES");
  if (!exchange_table.is_object() || exchange_table.empty()) {
    throw std::invalid_argument(
        "PARAM_TENSORPOM.EXCHANGES must be a nonempty object");
  }
  for (const auto &[key, block] : exchange_table.items()) {
    const std::string context = "PARAM_TENSORPOM.EXCHANGES." + key;
    ValidateFields(block, {"type", "delta", "ap", "eta", "scale"}, context);
    const int pdg = ParsePdgKey(key, context);
    const auto &particle = pdg_table.FindByPDG(pdg);
    const TensorExchangeType type =
        ParseExchangeType(block.at("type").get<std::string>(), context);
    const int rank = ExchangeRank(type);
    if ((rank == 2 &&
         (particle.spinX2 != 4 || particle.P != 1 || particle.C != 1)) ||
        (rank == 1 &&
         (particle.spinX2 != 2 || particle.P != -1 || particle.C != -1))) {
      throw std::invalid_argument(
          context + " fixed-spin PDG quantum numbers do not match type");
    }
    const double delta = block.at("delta").get<double>();
    const double ap = block.at("ap").get<double>();
    const int eta = block.contains("eta")
                        ? ParseUnitSign(block.at("eta"), context + " eta")
                        : 1;
    const double scale = block.value("scale", 1.0);
    if (!std::isfinite(delta) || !std::isfinite(ap) || ap <= 0.0 ||
        !std::isfinite(scale) || scale <= 0.0) {
      throw std::invalid_argument(
          context + " trajectory values must be finite and scales positive");
    }
    if (rank == 2 && (block.contains("eta") || block.contains("scale"))) {
      throw std::invalid_argument(
          context + " rank-two exchange must not set eta or scale");
    }
    if (type == TensorExchangeType::VectorReggeon && block.contains("eta")) {
      throw std::invalid_argument(context + " vector Reggeon must not set eta");
    }
    exchanges.emplace(
        pdg, TensorExchangeParam{pdg, type, rank, delta, ap, eta, scale});
  }

  const auto &continuum_param = general.at("PARAM_TENSORPOM").at("PARAM_CON");
  if (!continuum_param.is_object() || continuum_param.size() != 1 || !continuum_param.contains("TP")) {
    throw std::invalid_argument("PARAM_TENSORPOM.PARAM_CON must contain only TP");
  }
  const auto &selection = continuum_param.at("TP");
  if (!selection.is_object() || !selection.contains("[*]")) {
    throw std::invalid_argument(
        "PARAM_TENSORPOM.PARAM_CON.TP must be an object containing [*]");
  }
  std::set<std::array<int, 2>> final_seen;
  for (const auto &[key, pairs] : selection.items()) {
    TensorContinuumChannel channel;
    channel.fallback = key == "[*]";
    if (!channel.fallback) {
      channel.final_pdgs = CanonicalPair(ParsePairKey(key, "PARAM_TENSORPOM.PARAM_CON.TP"));
      if (!final_seen.insert(channel.final_pdgs).second) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.PARAM_CON.TP contains an order-equivalent final-state duplicate");
      }
    }
    if (!pairs.is_array() || pairs.empty()) {
      throw std::invalid_argument("PARAM_TENSORPOM.PARAM_CON.TP." + key +
                                  " must contain exchange pairs");
    }
    std::set<std::array<int, 2>> pair_seen;
    for (const auto &i : indices(pairs)) {
      const auto pair = ParsePairKey(pairs[i].dump(), "PARAM_TENSORPOM.PARAM_CON.TP." + key);
      FindExchange(pair[0]);
      FindExchange(pair[1]);
      const auto canonical = CanonicalPair(pair);
      if (!pair_seen.insert(canonical).second) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.PARAM_CON.TP." + key +
            " contains an unordered exchange duplicate");
      }
      channel.exchange_pairs.push_back(pair);
    }
    continuum.push_back(std::move(channel));
  }

  if (!vertex_table.is_object() || vertex_table.empty()) {
    throw std::invalid_argument("CON_TP.json must be a nonempty object");
  }
  for (const auto &[exchange_key, hadrons] : vertex_table.items()) {
    const int exchange_pdg = ParsePdgKey(exchange_key, "CON_TP.json");
    FindExchange(exchange_pdg);
    if (!hadrons.is_object() || hadrons.empty()) {
      throw std::invalid_argument("CON_TP.json." + exchange_key +
                                  " must be a nonempty object");
    }
    for (const auto &[pair_key, block] : hadrons.items()) {
      const std::string context =
          "CON_TP.json." + exchange_key + "." + pair_key;
      ValidateFields(block,
                     {"g_tensor", "FF_transfer", "FF_offshell", "threshold"},
                     context);
      const auto pair = ParsePairKey(pair_key, context);
      const auto canonical = CanonicalPair(pair);
      if (pair != canonical || pair[0] <= 0) {
        throw std::invalid_argument(
            context + " must identify one sorted positive-PDG pair");
      }
      const auto &first = pdg_table.FindByPDG(pair[0]);
      const auto &second = pdg_table.FindByPDG(pair[1]);
      const bool mixed = pair[0] != pair[1];
      const auto first_field = exchanges.find(pair[0]);
      const auto second_field = exchanges.find(pair[1]);
      const bool has_rank_one =
          (first_field != exchanges.end() && first_field->second.rank == 1) ||
          (second_field != exchanges.end() && second_field->second.rank == 1);
      if (mixed && (FindExchange(exchange_pdg).rank != 2 || first.spinX2 != 2 ||
                    second.spinX2 != 2 || !has_rank_one)) {
        throw std::invalid_argument(
            context +
            " mixed vertex requires a rank-two field and one rank-one field");
      }
      const auto coupling = block.at("g_tensor").get<std::vector<double>>();
      const std::size_t coupling_count =
          mixed ? 2 : CouplingCount(first.spinX2, context);
      if (coupling.size() != coupling_count ||
          !std::all_of(
              coupling.cbegin(), coupling.cend(),
              [](const double value) { return std::isfinite(value); })) {
        throw std::invalid_argument(
            context + " has invalid g_tensor cardinality or value");
      }
      const auto &crossed = !mixed && second.C == 0 ? pdg_table.FindByPDG(-second.pdg) : second;
      ValidateVertex(FindExchange(exchange_pdg), pdg_table.FindByPDG(exchange_pdg), first, crossed, context);
      const auto ff_transfer = regge::ReadFF(block.at("FF_transfer"), context + ".FF_transfer");
      if (ff_transfer.type != regge::FFType::None && ff_transfer.norm != regge::FFNorm::Zero) {
        throw std::invalid_argument(context + ".FF_transfer must use norm zero when active");
      }
      if (block.contains("threshold") && !block.at("threshold").is_boolean()) {
        throw std::invalid_argument(context + " threshold must be boolean");
      }
      const bool threshold = block.value("threshold", false);
      if (!mixed && block.contains("threshold")) {
        throw std::invalid_argument(context + " has invalid threshold field");
      }
      const auto ff_offshell = regge::ReadFF(block.at("FF_offshell"), context + ".FF_offshell");
      if (!mixed && ff_offshell.type != regge::FFType::None && ff_offshell.norm != regge::FFNorm::Pole) {
        throw std::invalid_argument(context + ".FF_offshell must use norm pole when active");
      }
      const std::array<int, 3> vertex_key = {exchange_pdg, pair[0], pair[1]};
      if (vertices.contains(vertex_key)) {
        throw std::invalid_argument(context + " duplicates a continuum vertex");
      }
      vertices.emplace(vertex_key,
                       TensorVertexParam{exchange_pdg, pair, coupling, {},
                                         ff_transfer, ff_offshell, threshold});
      if (!mixed) {
        hadron_spins.emplace(pair[0], first.spinX2);
      }
    }
  }

  for (const auto &channel : continuum) {
    if (channel.fallback) {
      continue;
    }
    const int hadron = std::abs(channel.final_pdgs[0]);
    for (const auto &pair : channel.exchange_pairs) {
      FindVertex(pair[0], hadron);
      FindVertex(pair[1], hadron);
    }
  }
  for (auto &[key, vertex] : vertices) {
    (void)key;
    vertex.active_g_tensor.clear();
    for (const auto &i : indices(vertex.g_tensor)) {
      if (std::abs(vertex.g_tensor[i]) > coupling_min) {
        vertex.active_g_tensor.push_back(i);
      }
    }
  }
  initialized = true;
}

// Compute one configured exchange by fixed-spin PDG alias
const TensorExchangeParam &
MTensorExchangeModel::FindExchange(const int pdg) const {
  const auto found = exchanges.find(pdg);
  if (found == exchanges.end()) {
    throw std::invalid_argument("Unknown tensor exchange PDG = " +
                                std::to_string(pdg));
  }
  return found->second;
}

// Compute the validated SOFT exchange matched by Lorentz and crossing type
SoftExchangeId MTensorExchangeModel::SoftId(const int pdg,
                                            const SoftModel &model) const {
  const TensorExchangeType type = FindExchange(pdg).type;
  std::string name;
  SoftExchangeRole role = SoftExchangeRole::Pomeron;
  int crossing = 1;
  switch (type) {
  case TensorExchangeType::Pomeron:
    name = "P";
    role = SoftExchangeRole::Pomeron;
    crossing = 1;
    break;
  case TensorExchangeType::TensorReggeon:
    name = "R_f2";
    role = SoftExchangeRole::Reggeon;
    crossing = 1;
    break;
  case TensorExchangeType::Odderon:
    name = "O";
    role = SoftExchangeRole::Odderon;
    crossing = -1;
    break;
  case TensorExchangeType::VectorReggeon:
    name = "R_rho";
    role = SoftExchangeRole::Reggeon;
    crossing = -1;
    break;
  }
  const SoftExchangeId exchange = model.ExchangeId(name);
  const auto &soft = model.Exchange(exchange);
  if (!soft.Enabled() || soft.Role() != role ||
      soft.CrossingParity() != crossing) {
    throw std::invalid_argument(
        "MTensorExchangeModel::SoftId: incompatible SOFT exchange " + name);
  }
  return exchange;
}

// Compute one configured exchange-to-hadron continuum vertex
const TensorVertexParam &
MTensorExchangeModel::FindVertex(const int exchange_pdg,
                                 const int hadron_pdg) const {
  return FindVertex(exchange_pdg, std::abs(hadron_pdg), std::abs(hadron_pdg));
}

// Compute one configured exchange-to-particle-pair vertex
const TensorVertexParam &
MTensorExchangeModel::FindVertex(const int exchange_pdg, const int first_pdg,
                                 const int second_pdg) const {
  const auto pair = CanonicalPair({std::abs(first_pdg), std::abs(second_pdg)});
  const auto found = vertices.find({exchange_pdg, pair[0], pair[1]});
  if (found == vertices.end()) {
    throw std::invalid_argument(
        "Missing tensor continuum vertex for exchange PDG = " +
        std::to_string(exchange_pdg) + " and particle PDGs = [" +
        std::to_string(first_pdg) + "," + std::to_string(second_pdg) + "]");
  }
  return found->second;
}

// Compute the selected unordered exchange pairs for one final-state pair
const std::vector<std::array<int, 2>> &
MTensorExchangeModel::FindContinuumPairs(const int first_pdg,
                                         const int second_pdg) const {
  const auto key = CanonicalPair({first_pdg, second_pdg});
  const TensorContinuumChannel *fallback = nullptr;
  for (const auto &channel : continuum) {
    if (channel.fallback) {
      fallback = &channel;
    } else if (channel.final_pdgs == key) {
      return channel.exchange_pairs;
    }
  }
  if (fallback == nullptr) {
    throw std::invalid_argument("PARAM_TENSORPOM.PARAM_CON.TP has no [*] channel");
  }
  return fallback->exchange_pairs;
}

// Compute common internal lines implied by two ordinary pair vertices
std::vector<int> MTensorExchangeModel::FindTransfers(
    const int first_exchange, const int second_exchange, const int first_hadron,
    const int second_hadron) const {
  const int first = std::abs(first_hadron);
  const int second = std::abs(second_hadron);
  std::set<int> transfers;
  for (const auto &[key, vertex] : vertices) {
    if (key[0] != first_exchange) {
      continue;
    }
    int transfer = 0;
    if (vertex.pair_pdgs[0] == first) {
      transfer = vertex.pair_pdgs[1];
    } else if (vertex.pair_pdgs[1] == first) {
      transfer = vertex.pair_pdgs[0];
    } else {
      continue;
    }
    const auto pair = CanonicalPair({second, transfer});
    if (vertices.contains({second_exchange, pair[0], pair[1]})) {
      transfers.insert(transfer);
    }
  }
  return {transfers.cbegin(), transfers.cend()};
}

// Compute active common internal lines implied by two pair vertices
std::vector<int> MTensorExchangeModel::FindActiveTransfers(
    const int first_exchange, const int second_exchange, const int first_hadron,
    const int second_hadron) const {
  std::vector<int> out;
  for (const int transfer : FindTransfers(first_exchange, second_exchange,
                                          first_hadron, second_hadron)) {
    const auto &first = FindVertex(first_exchange, first_hadron, transfer);
    const auto &second = FindVertex(second_exchange, second_hadron, transfer);
    if (!first.active_g_tensor.empty() && !second.active_g_tensor.empty()) {
      out.push_back(transfer);
    }
  }
  return out;
}

// Expand one unordered exchange pair into its physical beam orderings
std::vector<std::array<int, 2>>
MTensorExchangeModel::OrderedPairs(const std::array<int, 2> &pair) {
  if (pair[0] == pair[1]) {
    return {pair};
  }
  return {pair, {pair[1], pair[0]}};
}

// Compute the scalar Regge factor of one soft exchange propagator
std::complex<double>
MTensorExchangeModel::PropagatorFactor(const int exchange_pdg, const double s,
                                       const double t) const {
  const auto &exchange = FindExchange(exchange_pdg);
  if (!std::isfinite(s) || !std::isfinite(t) || s <= 0.0) {
    throw AmplitudeFailure(
        "Tensor exchange propagator requires finite s > 0 and finite t");
  }
  const double alpha = 1.0 + exchange.delta + exchange.ap * t;
  const std::complex<double> power =
      std::pow(-math::zi * s * exchange.ap, alpha - 1.0);
  if (exchange.rank == 2) {
    return power / (4.0 * s);
  }
  if (exchange.type == TensorExchangeType::Odderon) {
    return -math::zi * static_cast<double>(exchange.eta) * power /
           std::pow(exchange.scale, 2);
  }
  return math::zi * power / std::pow(exchange.scale, 2);
}

// Convert a card coupling to the common baryon tensor-vertex convention
double MTensorExchangeModel::BaryonCoupling(const int exchange_pdg,
                                            const double coupling) const {
  const auto &exchange = FindExchange(exchange_pdg);
  if (exchange.rank != 2) {
    throw std::invalid_argument(
        "Tensor baryon coupling requires a rank-two exchange");
  }
  return exchange.type == TensorExchangeType::Pomeron ? coupling
                                                      : coupling / 3.0;
}

// Convert a card coupling to the common pseudoscalar tensor vertex
double MTensorExchangeModel::PseudoscalarCoupling(const int exchange_pdg,
                                                  const double coupling) const {
  const auto &exchange = FindExchange(exchange_pdg);
  if (exchange.rank != 2) {
    throw std::invalid_argument(
        "Tensor pseudoscalar coupling requires a rank-two exchange");
  }
  return exchange.type == TensorExchangeType::Pomeron ? coupling
                                                      : coupling / 4.0;
}

// Convert a card coupling to the common vector baryon vertex
double MTensorExchangeModel::VectorBaryonCoupling(const int exchange_pdg,
                                                  const double coupling) const {
  const auto &exchange = FindExchange(exchange_pdg);
  if (exchange.rank != 1) {
    throw std::invalid_argument(
        "Vector baryon coupling requires a rank-one exchange");
  }
  return exchange.type == TensorExchangeType::Odderon
             ? 3.0 * exchange.scale * coupling
             : coupling;
}

// Convert a card coupling to the common vector pseudoscalar vertex
double
MTensorExchangeModel::VectorPseudoscalarCoupling(const int exchange_pdg,
                                                 const double coupling) const {
  const auto &exchange = FindExchange(exchange_pdg);
  if (exchange.type != TensorExchangeType::VectorReggeon) {
    throw std::invalid_argument(
        "Vector pseudoscalar coupling requires a vector Reggeon");
  }
  return coupling / 2.0;
}

// Compute all continuum hadrons with the requested spinX2
std::vector<int> MTensorExchangeModel::HadronPdgs(const int spinX2) const {
  std::vector<int> out;
  for (const auto &[pdg, spin] : hadron_spins) {
    if (spin == spinX2) {
      out.push_back(pdg);
    }
  }
  return out;
}

} // namespace gra
