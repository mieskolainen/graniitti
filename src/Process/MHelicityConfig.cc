// Helicity steering-card loading and vertex construction
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <complex>
#include <cstdio>
#include <filesystem>
#include <iostream>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Tensor/MTensorPomeron.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MSpinDensity.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra {
namespace {

// Parse one integer signature without narrowing oversized JSON values
regge::Signature ParseTrajectorySignature(const nlohmann::json &value, const std::string &path) {
  if (!value.is_number_integer()) { throw std::invalid_argument(path + " must be -1 or 1"); }
  if (value.is_number_unsigned()) {
    if (value.get<nlohmann::json::number_unsigned_t>() == 1U) { return regge::Signature::Positive; }
  } else {
    const auto tau = value.get<nlohmann::json::number_integer_t>();
    if (tau == 1) { return regge::Signature::Positive; }
    if (tau == -1) { return regge::Signature::Negative; }
  }
  throw std::invalid_argument(path + " must be -1 or 1");
}

}  // namespace

// Load immutable helicity steering cards for one model tune
MHelicityConfig::MHelicityConfig(MModelTunePtr tune) { Load(std::move(tune)); }

// Load immutable helicity steering cards for one model tune
void MHelicityConfig::Load(MModelTunePtr model_tune) {
  if (model_tune == nullptr) { throw std::invalid_argument("MHelicityConfig::Load: null model tune"); }
  tune = std::move(model_tune);
  continuum.clear();
  for (const std::string model : {"MP", "XP", "GP"}) { continuum[model] = &tune->Continuum(model); }
  decays_path = (std::filesystem::path(tune->GeneralFile()).parent_path() / "DECAYS.json").string();
  try {
    decays = std::make_shared<const nlohmann::json>(nlohmann::json::parse(gra::aux::GetInputData(decays_path)));
  } catch (const std::exception &error) {
    throw std::invalid_argument("MHelicityConfig::Load: Error reading " + decays_path + ": " + error.what());
  }
  trajectory_signature.clear();
  trajectory_pole_spinX2.clear();
  const auto &general = tune->General();
  const auto &soft    = general.at("PARAM_SOFT").at("EXCHANGE_DEF");
  for (const auto &row : general.at("PARAM_REGGE").at("EXCHANGES")) {
    const std::string      exchange  = row.at("soft_exchange").get<std::string>();
    const auto            &tau       = soft.at(exchange).at("tau");
    const regge::Signature signature = ParseTrajectorySignature(tau, "MHelicityConfig::Load: trajectory signature");
    if (!row.at("pole_spin").is_number_integer()) {
      throw std::invalid_argument("MHelicityConfig::Load: pole_spin must be an integer");
    }
    const int pole_spin = row.at("pole_spin").get<int>();
    if (pole_spin <= 0 || pole_spin > std::numeric_limits<int>::max() / 2) {
      throw std::invalid_argument("MHelicityConfig::Load: pole_spin must be positive");
    }
    for (const int pdg : row.at("pdg").get<std::vector<int>>()) {
      const auto [position, inserted] = trajectory_signature.emplace(pdg, signature);
      if (!inserted && position->second != signature) {
        throw std::invalid_argument("MHelicityConfig::Load: exchange PDG has conflicting signatures");
      }
      const auto [pole_position, pole_inserted] = trajectory_pole_spinX2.emplace(pdg, 2 * pole_spin);
      if (!pole_inserted && pole_position->second != 2 * pole_spin) {
        throw std::invalid_argument("MHelicityConfig::Load: exchange PDG has conflicting pole spins");
      }
    }
  }
  ValidateContinuum();
}

// Compute whether the helicity steering cards have been loaded
bool MHelicityConfig::IsLoaded() const noexcept {
  return continuum.size() == 3 && continuum.count("MP") != 0 && continuum.count("XP") != 0 &&
         continuum.count("GP") != 0 && continuum.at("MP") != nullptr && continuum.at("XP") != nullptr &&
         continuum.at("GP") != nullptr && decays != nullptr && !trajectory_signature.empty() &&
         !trajectory_pole_spinX2.empty();
}

// Read decay-only parameters in a two-body channel block
//
void ParseDecayParameters(HELMatrix &hc, const json &channel_block, const std::string &context, bool production_mode,
                          const std::string &model) {
  hc.BR                = 1.0;
  hc.BR_set            = false;
  hc.zeta              = 0.0;
  hc.g_decay_TP.clear();

  if (production_mode) { return; }

  if (!channel_block.contains("BR") || !channel_block["BR"].is_number()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " is missing numeric BR");
  }
  if (!channel_block.contains("zeta") || !channel_block.at("zeta").is_object() || channel_block.at("zeta").size() != 4) {
    throw std::invalid_argument(context + ": zeta must contain MP, XP, GP and TP");
  }
  const auto &phases = channel_block.at("zeta");
  for (const std::string key : {"MP", "XP", "GP", "TP"}) {
    if (!phases.contains(key) || !phases.at(key).is_number() || !std::isfinite(phases.at(key).get<double>())) {
      throw std::invalid_argument(context + ": zeta." + key + " must be finite and numeric");
    }
  }

  hc.BR     = channel_block["BR"].get<double>();
  hc.BR_set = true;
  // Other production processes have no Pomeron model decay phase
  if (phases.contains(model)) { hc.zeta = phases.at(model).get<double>(); }

  if (!std::isfinite(hc.BR) || hc.BR < 0.0 || hc.BR > 1.0) {
    throw std::invalid_argument(context + ": BR must be finite and in [0,1]");
  }

  if (!channel_block.contains("FF_decay")) {
    throw std::invalid_argument(context + ": missing FF_decay for MP, XP, GP and TP");
  }
  const auto &forms = channel_block.at("FF_decay");
  if (!forms.is_object() || forms.size() != 4) {
    throw std::invalid_argument(context + ": FF_decay must contain MP, XP, GP and TP");
  }
  for (const std::string key : {"MP", "XP", "GP", "TP"}) {
    if (!forms.contains(key)) { throw std::invalid_argument(context + ": missing FF_decay." + key); }
    const auto form = regge::ReadFF(forms.at(key), context + ".FF_decay." + key);
    if (form.type != regge::FFType::None &&
        ((form.type != regge::FFType::Gaussian && form.type != regge::FFType::Vector) ||
         form.norm != regge::FFNorm::Pole)) {
      throw std::invalid_argument(context + ".FF_decay." + key + " must be gaussian or vector with norm pole");
    }
    if (key == model) { hc.ff_decay = form; }
  }


  if (!channel_block.contains("g_decay_TP")) { return; }
  if (!channel_block["g_decay_TP"].is_array()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " has non-array g_decay_TP");
  }

  for (const auto &value : channel_block["g_decay_TP"]) {
    if (!value.is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " has non-numeric g_decay_TP value");
    }
    const double coupling = value.get<double>();
    if (!std::isfinite(coupling)) {
      throw std::invalid_argument(context + ": g_decay_TP must be finite");
    }
    hc.g_decay_TP.push_back(coupling);
  }
}

// Read LS-basis rows alpha_{lS} into sparse storage
//
void ParseAlphaLSRows(HELMatrix &hc, const json &channel_block, const std::string &context,
                      const std::string &table_name, bool production_mode) {
  const std::string field = production_mode ? "g_ls" : "alpha_ls";
  if (!channel_block.contains(field) || !channel_block[field].is_array()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " is missing array " + field);
  }
  if (channel_block[field].empty()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " " + field +
                                " array is empty");
  }

  hc.alpha_ls.Clear();
  for (const auto &a : indices(channel_block[field])) {
    const auto &row = channel_block[field][a];
    if (!row.is_array() || row.size() != 4) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " has malformed " + field +
                                  " row " + std::to_string(a) + " (expected [l, s, mag, phase])");
    }
    if (!row[0].is_number_integer() || !row[1].is_number() || !row[2].is_number() || !row[3].is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " has non-numeric " +
                                  field + " row " + std::to_string(a));
    }

    const long long l_raw = row[0].get<long long>();
    const double    s     = row[1].get<double>();
    const double    mag   = row[2].get<double>();
    const double    phase = row[3].get<double>();

    if (l_raw < 0 || l_raw > std::numeric_limits<int>::max() / 2) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " has invalid l in " +
                                  field + " row " + std::to_string(a));
    }

    const std::size_t l        = static_cast<std::size_t>(l_raw);
    const double      _2s_real = 2.0 * s;
    if (!std::isfinite(s) || !std::isfinite(mag) || !std::isfinite(phase) || s < 0.0 || mag < 0.0 ||
        _2s_real > static_cast<double>(std::numeric_limits<int>::max()) ||
        std::abs(_2s_real - std::round(_2s_real)) > 1e-9) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " has invalid " + field +
                                  " row " + std::to_string(a));
    }

    const std::size_t          _2s         = static_cast<std::size_t>(std::llround(_2s_real));
    const std::complex<double> coefficient = std::polar(mag, phase);
    if (!hc.alpha_ls.Insert(l, _2s, coefficient)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " duplicates " + field +
                                  " entry [l=" + std::to_string(l) + ", s=" + std::to_string(s) + "]");
    }
    std::cout << rang::fg::green << "LS couplings: " << table_name << ": "
              << "[l=" << l << ", s=" << s << ", mag=" << std::abs(coefficient) << ", phase=" << std::arg(coefficient)
              << "]" << rang::fg::reset << std::endl;
  }
}

// Validate every explicit finite-spin LS row before numerical pruning
void ValidateFiniteLSRows(const HELMatrix &hc, const MParticle &mother, const std::vector<MParticle> &legs,
                          bool production_mode, const std::string &source_name, const std::string &context,
                          gra::spin::VertexContext vertex_context) {
  if (legs.size() != 2) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " LS validation requires two legs");
  }
  const std::vector<spin::LSCoupling> allowed =
      spin::AllowedLSCouplings(hc, mother, legs[0], legs[1], production_mode, source_name, vertex_context);
  for (const spin::LSTerm &term : hc.alpha_ls) {
    const bool found = std::any_of(allowed.cbegin(), allowed.cend(), [&term](const spin::LSCoupling &row) {
      return row.l == term.l && row.two_s == term.two_s;
    });
    if (!found) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " has forbidden LS row [l=" + std::to_string(term.l) +
                                  ", s=" + gra::aux::Spin2XtoString(static_cast<int>(term.two_s)) + "]");
    }
  }
}

// Validate continued continuum LS rows using physical crossed-leg metadata
void ValidateCrossedLSRows(const HELMatrix &hc, const MParticle &mother, const std::vector<MParticle> &legs,
                           bool production_mode, const std::string &context, gra::spin::VertexContext vertex_context) {
  if (legs.size() != 2 || legs[0].spinX2 < 0 || legs[1].spinX2 < 0) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_ls requires two fixed-spin physical legs");
  }
  const auto validate = [&](const gra::spin::LSTerm &term) {
    const std::size_t l     = term.l;
    const std::size_t two_s = term.two_s;
    if (two_s % 2 != 0 || two_s < static_cast<std::size_t>(std::abs(legs[0].spinX2 - legs[1].spinX2)) ||
        two_s > static_cast<std::size_t>(legs[0].spinX2 + legs[1].spinX2) ||
        (legs[0].spinX2 + legs[1].spinX2 + static_cast<int>(two_s)) % 2 != 0) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls row has forbidden physical total spin S");
    }
    const int S = static_cast<int>(two_s / 2);
    if (hc.P_symmetry && mother.P != legs[0].P * legs[1].P * ((l % 2 == 0) ? 1 : -1)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls row violates parity");
    }
    if (!gra::spin::TwoBodyCParityAllowed(mother, legs[0], legs[1], l, static_cast<double>(S), hc.C_symmetry)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls row violates C parity");
    }
    if (!gra::spin::PhysicalIdenticalPair(legs[0], legs[1], production_mode, vertex_context)) { return; }
    if (legs[0].spinX2 % 2 == 0 && !gra::spin::BoseSymmetry(l, S)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls row violates Bose symmetry");
    }
    if (legs[0].spinX2 % 2 != 0 && !gra::spin::FermiSymmetry(l, S)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls row violates Fermi symmetry");
    }
  };

  if (hc.exchange_basis == ExchangeBasisType::ReggeHelicity) {
    const std::size_t nm = static_cast<std::size_t>(2 * hc.analytic_MMAX + 1);
    if (hc.m_ls.size() != nm) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls per-m basis has invalid dimensions");
    }
    for (const spin::LSCoefficients &terms : hc.m_ls) {
      for (const gra::spin::LSTerm &term : terms) { validate(term); }
    }
    return;
  }
  for (const gra::spin::LSTerm &term : hc.alpha_ls) { validate(term); }
}

// Collect particle PDG identifiers
//
std::vector<int> ParticlePDGList(const std::vector<MParticle> &legs) {
  std::vector<int> out;
  out.reserve(legs.size());
  for (const auto &leg : legs) { out.push_back(leg.pdg); }
  return out;
}

// Format PDG list
//
std::string FormatPDGList(const std::vector<int> &values) {
  std::string out = "[";
  for (const auto &i : indices(values)) {
    out += std::to_string(values[i]);
    if (i + 1 < values.size()) { out += " "; }
  }
  out += "]";
  return out;
}

// Print helper
//
std::string ChannelContext(const std::string &table_name, const std::string &pdg_str, const std::string &channel_id) {
  return table_name + " entry for PDG = " + pdg_str + ", channel[" + channel_id + "]";
}

// Check helper
//
bool IsParticleAntiparticleReverse(const std::vector<int> &requested, const std::vector<int> &stored) {
  return requested.size() == 2 && stored.size() == 2 && requested[0] == stored[1] && requested[1] == stored[0] &&
         stored[0] == -stored[1];
}

// Check if two PDG lists are charge conjugates in the same leg order
//
bool IsChargeConjugateSameOrder(const std::vector<int> &requested, const std::vector<int> &stored) {
  if (requested.size() != stored.size()) { return false; }
  for (const auto &i : indices(requested)) {
    if (requested[i] != -stored[i]) { return false; }
  }
  return true;
}

// Check whether two two-body vertices are the same physical vertex with legs
// exchanged
//
bool IsTwoBodyLegExchange(const std::vector<int> &requested, const std::vector<int> &stored) {
  return requested.size() == 2 && stored.size() == 2 && requested[0] == stored[1] && requested[1] == stored[0];
}

// Compute the canonical positive production pair used by continuum cards
//
std::vector<int> ProductionPairKey(const std::vector<int> &pdgs) {
  std::vector<int> key;
  key.reserve(pdgs.size());
  for (const int pdg : pdgs) { key.push_back(std::abs(pdg)); }
  std::sort(key.begin(), key.end());
  return key;
}

// Test whether a production leg is a self-conjugate neutral state
//
bool IsSelfConjugateProductionLeg(const MParticle &leg) {
  return leg.pdg == PDG::PDG_gamma || leg.pdg == PDG::PDG_gluon || (leg.chargeX3 == 0 && (leg.C == -1 || leg.C == 1));
}

// Compute the production sector selected by signed physical production legs
//
std::string ProductionSector(const std::vector<int> &requested, const std::vector<MParticle> &requested_legs) {
  if (requested.size() != 2 || requested_legs.size() != 2) {
    throw std::invalid_argument(
        "MHelicityConfig::ProcessHelicityStructure: "
        "central continuum card pair schema "
        "expects two legs");
  }
  if (IsSelfConjugateProductionLeg(requested_legs[0]) && IsSelfConjugateProductionLeg(requested_legs[1])) {
    return "self";
  }

  const bool first_positive  = requested[0] > 0;
  const bool second_positive = requested[1] > 0;
  return first_positive == second_positive ? "same" : "opposite";
}

// Compute whether a key names one supported production sector
//
bool IsProductionSectorName(const std::string &key) { return key == "same" || key == "opposite" || key == "self"; }

// Validate the two-body C and parity symmetry flag array
//
void ValidateTwoBodySymmetryFlags(const json &channel_block, const std::string &context) {
  if (!channel_block.contains("CP") || !channel_block["CP"].is_array() || channel_block["CP"].size() != 2 ||
      !channel_block["CP"][0].is_boolean() || !channel_block["CP"][1].is_boolean()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " requires CP = [C_symmetry, P_symmetry]");
  }
}

// Parse one compact PDG-array object key
//
std::vector<int> ParsePDGKey(const std::string &key, const std::string &context) {
  if (key.size() < 3 || key.front() != '[' || key.back() != ']') {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " key must be a compact PDG array");
  }

  std::vector<int> pdgs;
  std::size_t      begin = 1;
  while (begin < key.size() - 1) {
    const std::size_t comma = key.find(',', begin);
    const std::size_t end   = comma == std::string::npos ? key.size() - 1 : comma;
    if (end == begin) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " key contains an empty PDG value");
    }
    try {
      std::size_t       used  = 0;
      const std::string token = key.substr(begin, end - begin);
      const int         value = std::stoi(token, &used);
      if (used != token.size()) { throw std::invalid_argument("trailing characters"); }
      pdgs.push_back(value);
    } catch (const std::exception &) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " key contains a non-integer PDG value");
    }
    if (comma == std::string::npos) { break; }
    begin = comma + 1;
  }
  std::string canonical = "[";
  for (const auto &i : indices(pdgs)) {
    canonical += std::to_string(pdgs[i]);
    canonical += i + 1 < pdgs.size() ? "," : "]";
  }
  if (pdgs.empty() || canonical != key) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " key must use compact integer notation");
  }
  return pdgs;
}

// Compute whether one central vertex selects automatic spin coupling
//
bool IsAutomaticCoupling(const std::string &basis, ReggeVertexRole role) {
  return IsAutomaticReggeVertexBasis(ParseReggeVertexBasis(basis, role));
}

// Read one complete complex LS operator table
std::vector<spin::LSTerm> ReadCanonicalPoleTerms(const json &block, const std::string &field,
                                                 const std::string &context) {
  if (!block.contains(field) || !block.at(field).is_array() || block.at(field).empty()) {
    throw std::invalid_argument(context + " requires a non-empty " + field);
  }
  std::vector<spin::LSTerm> terms;
  terms.reserve(block.at(field).size());
  for (const auto &i : indices(block.at(field))) {
    const auto &row = block.at(field).at(i);
    if (!row.is_array() || (row.size() != 3 && row.size() != 4) || !row.at(0).is_number_integer() ||
        !row.at(1).is_number() || !row.at(2).is_number() || (row.size() == 4 && !row.at(3).is_number())) {
      throw std::invalid_argument(context + " has malformed " + field + " row " + std::to_string(i));
    }
    const long long l         = row.at(0).get<long long>();
    const double    s         = row.at(1).get<double>();
    const double    magnitude = row.at(2).get<double>();
    const double    phase     = row.size() == 4 ? row.at(3).get<double>() : 0.0;
    const double    two_s     = 2.0 * s;
    if (l < 0 || l > std::numeric_limits<int>::max() / 2 || !std::isfinite(s) || s < 0.0 || !std::isfinite(magnitude) ||
        magnitude < 0.0 || !std::isfinite(phase) || two_s > static_cast<double>(std::numeric_limits<int>::max()) ||
        std::abs(two_s - std::round(two_s)) > 1.0e-9) {
      throw std::invalid_argument(context + " has invalid " + field + " row " + std::to_string(i));
    }
    terms.push_back(
        {static_cast<std::size_t>(l), static_cast<std::size_t>(std::llround(two_s)), std::polar(magnitude, phase)});
  }
  return terms;
}

// Apply the canonical LS phase induced by exchanging the two physical legs
void ApplyCanonicalPoleLegExchange(std::vector<spin::LSTerm> &terms, const MParticle &leg1, const MParticle &leg2,
                                   const bool exchanged, const std::string &context) {
  if (!exchanged) { return; }
  for (auto &term : terms) {
    const long long twice_exponent =
        2 * static_cast<long long>(term.l) + leg1.spinX2 + leg2.spinX2 - static_cast<long long>(term.two_s);
    if (twice_exponent % 2 != 0) { throw std::invalid_argument(context + " has non-integer leg-exchange phase"); }
    if ((twice_exponent / 2) % 2 != 0) { term.coefficient *= -1.0; }
  }
}

// Validate one selected sector block before normal coupling parsing
//
void ValidateProductionSectorBlock(const json &sector_block, const std::string &context, bool analytic_production) {
  if (!sector_block.is_object()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production sector is not a JSON object");
  }
  if (!sector_block.contains("basis") || !sector_block["basis"].is_string()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production sector is missing string basis");
  }
  ValidateTwoBodySymmetryFlags(sector_block, context);

  const std::string      basis_name = sector_block["basis"].get<std::string>();
  const ReggeVertexBasis basis      = ParseReggeVertexBasis(basis_name, ReggeVertexRole::Continuum);
  if (analytic_production && (!sector_block["CP"][0].get<bool>() || !sector_block["CP"][1].get<bool>())) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " analytic Regge production requires CP = [true, true]");
  }
  std::set<std::string> allowed_fields = {"basis", "CP"};
  if (IsAutomaticReggeVertexBasis(basis)) {
    if (!sector_block.contains("g") || !sector_block.at("g").is_array() || sector_block.at("g").size() != 2 ||
        !sector_block.at("g").at(0).is_number() || !sector_block.at("g").at(1).is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " automatic sector requires g = [magnitude, phase]");
    }
    allowed_fields.insert("g");
    for (const auto &[field, value] : sector_block.items()) {
      (void)value;
      if (!allowed_fields.contains(field)) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " automatic sector has inactive field " + field);
      }
    }
    return;
  }
  if (basis == ReggeVertexBasis::LS) {
    if (!sector_block.contains("g_ls") || !sector_block["g_ls"].is_array() || sector_block["g_ls"].empty()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " production sector is missing non-empty g_ls");
    }
    allowed_fields.insert("g_ls");
    for (const auto &[field, value] : sector_block.items()) {
      (void)value;
      if (!allowed_fields.contains(field)) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " g_ls sector has inactive field " + field);
      }
    }
    return;
  }
  if (basis != ReggeVertexBasis::Helicity) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production sector has an unsupported basis");
  }
  const std::string field = "helicity";
  if (!sector_block.contains(field) || !sector_block[field].is_array() || sector_block[field].empty()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production sector is missing non-empty " + field);
  }
  allowed_fields.insert(field);
  for (const auto &[field, value] : sector_block.items()) {
    (void)value;
    if (!allowed_fields.contains(field)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " sector basis " +
                                  basis_name + " has inactive field " + field);
    }
  }
}

// Validate one central continuum card pair block
//
std::vector<int> ValidateProductionPairBlock(const std::string &pair_key, const json &pair_block,
                                             const std::string &context, bool analytic_production) {
  if (!pair_block.is_object()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production pair is not a JSON object");
  }
  if (pair_block.contains("basis")) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " contains 'basis'; expected pair sectors");
  }

  std::vector<int> pair = ParsePDGKey(pair_key, context);
  if (pair.size() != 2) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production pair key must contain two PDGs");
  }
  for (const int pdg : pair) {
    if (pdg <= 0) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " production pair PDGs must be positive");
    }
  }
  if (pair[0] > pair[1]) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " production pair PDGs must use ascending order");
  }

  regge::ValidateContinuumControls(pair_block, context);
  bool        found_sector = false;
  std::size_t sector_count = 0;
  bool        has_self     = false;
  for (auto it = pair_block.begin(); it != pair_block.end(); ++it) {
    const std::string key = it.key();
    if (key == "reggeize" || key == "pveto") { continue; }
    if (key == "FF_transfer") {
      const auto ff = regge::ReadFF(it.value(), context + ".FF_transfer");
      if (ff.type != regge::FFType::None && ff.norm != regge::FFNorm::Zero) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    ".FF_transfer must use norm zero when active");
      }
      continue;
    }
    if (key == "FF_offshell") {
      const auto ff = regge::ReadFF(it.value(), context + ".FF_offshell");
      if (ff.type != regge::FFType::None && ff.norm != regge::FFNorm::Pole) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    ".FF_offshell must use norm pole when active");
      }
      continue;
    }
    if (!IsProductionSectorName(key)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " has unsupported production sector '" + key + "'");
    }
    found_sector = true;
    ++sector_count;
    has_self = has_self || key == "self";
    ValidateProductionSectorBlock(it.value(), context + "/" + key, analytic_production);
  }
  if (!found_sector) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " has no same/opposite/self sector");
  }
  if (has_self && sector_count != 1) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " self sector cannot coexist with same or opposite sectors");
  }
  return pair;
}

// Validate every MP, XP and GP continuum card before channel lookup
void MHelicityConfig::ValidateContinuum() const {
  for (const std::string model : {"MP", "XP", "GP"}) {
    const json &card = *continuum.at(model);
    if (!card.is_object() || card.empty()) { throw std::invalid_argument("CON_" + model + ".json must be nonempty"); }
    for (const auto &[exchange, block] : card.items()) {
      std::size_t used = 0;
      int         pdg  = 0;
      try {
        pdg = std::stoi(exchange, &used);
      } catch (const std::exception &) {
        throw std::invalid_argument("CON_" + model + ".json invalid exchange " + exchange);
      }
      if (pdg <= 0 || used != exchange.size() || !block.is_object() || block.empty()) {
        throw std::invalid_argument("CON_" + model + ".json invalid exchange " + exchange);
      }
      for (const auto &[pair, sectors] : block.items()) {
        if (pair == "default_g") { throw std::invalid_argument("CON_" + model + ".json fallbacks are not supported"); }
        ValidateProductionPairBlock(pair, sectors, "CON_" + model + ".json." + exchange, model == "GP");
      }
    }
  }
}

// Compute true if two PDG vectors have identical values
//
bool SamePDGVector(const std::vector<int> &lhs, const std::vector<int> &rhs) {
  return lhs.size() == rhs.size() && std::equal(lhs.begin(), lhs.end(), rhs.begin());
}

// Test if a two-body vertex has exchangeable physical legs
//
bool AllowTwoBodyLegExchange(bool production_mode, gra::spin::VertexContext context, const std::vector<int> &refdecay) {
  if (refdecay.size() != 2) { return false; }
  if (!production_mode) { return true; }
  return context == gra::spin::VertexContext::Auto;
}

// Test if the auxiliary Jacob-Wick production vertex crosses its second leg
//
bool UsesCrossedSecondLeg(bool production_mode, gra::spin::VertexContext context) {
  return production_mode && context == gra::spin::VertexContext::CrossedBeamLeg;
}

// Compute antiparticle data for an outgoing leg crossed into the auxiliary JW
// vertex
//
MParticle CrossedLegMetadata(const MParticle &leg, const MPDG &pdg_table) {
  const auto anti = pdg_table.PDG_table.find(-leg.pdg);
  if (anti != pdg_table.PDG_table.end()) { return anti->second; }

  if (leg.chargeX3 == 0 && (leg.C == -1 || leg.C == 1 || leg.pdg == PDG::PDG_gamma || leg.pdg == PDG::PDG_gluon)) {
    return leg;
  }

  throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: cannot cross PDG " + std::to_string(leg.pdg) +
                              " because antiparticle PDG metadata is missing");
}

// Build auxiliary JW metadata, crossing the second production leg as an
// antiparticle
//
std::vector<MParticle> JacobWickVertexLegs(const std::vector<MParticle> &legs, bool production_mode,
                                           gra::spin::VertexContext context, const MPDG &pdg_table) {
  std::vector<MParticle> vertex_legs = legs;
  if (vertex_legs.size() == 2 && UsesCrossedSecondLeg(production_mode, context)) {
    vertex_legs[1] = CrossedLegMetadata(vertex_legs[1], pdg_table);
  }
  return vertex_legs;
}

// Build internal beam-side photon source metadata without central continuum
// card rows
//
HELMatrix ForwardPhotonSourceHelicityStructure(const MParticle &p, const std::vector<MParticle> &legs) {
  if (p.pdg != PDG::PDG_gamma || legs.size() != 2) {
    throw std::invalid_argument(
        "MHelicityConfig::ForwardPhotonSourceHelicityStructure: "
        "expected gamma <- beam beam");
  }

  // Trace the nuclear spectator spin analytically, retaining physical particle metadata
  const int spin1_x2 = nuclear::IsNuclearPDG(legs[0].pdg) ? 0 : legs[0].spinX2;
  const int spin2_x2 = nuclear::IsNuclearPDG(legs[1].pdg) ? 0 : legs[1].spinX2;
  if (!((spin1_x2 == 1 && spin2_x2 == 1) || (spin1_x2 == 0 && spin2_x2 == 0))) {
    throw std::invalid_argument(
        "MHelicityConfig::ForwardPhotonSourceHelicityStructure: "
        "EPA photon source expects spin-zero or spin-half beams");
  }
  const double s1 = spin1_x2 / 2.0;
  const double s2 = spin2_x2 / 2.0;

  HELMatrix hc;
  hc.BR         = 1.0;
  hc.zeta       = 0.0;
  hc.P_symmetry = false;
  hc.C_symmetry = false;
  hc.alpha_ls.Clear();
  // A real photon source contains only the two transverse pole helicities
  gra::spin::InitTwoBodyBasis(hc, 1.0, s1, s2, {-1.0, 1.0},
                              gra::spin::SpinProjections(s1),
                              gra::spin::SpinProjections(s2),
                              "forward photon source");
  return hc;
}

// Build beam-side fixed-spin exchange metadata without a continuum card row
//
HELMatrix ForwardHadronSourceHelicityStructure(const MParticle &p, const std::vector<MParticle> &legs) {
  if (legs.size() != 2) {
    throw std::invalid_argument(
        "MHelicityConfig::ForwardHadronSourceHelicityStructure: "
        "expected exchange <- beam beam");
  }
  if (p.spinX2 < 0 || p.spinX2 % 2 != 0) {
    throw std::invalid_argument(
        "MHelicityConfig::ForwardHadronSourceHelicityStructure: "
        "fixed Regge exchange spin must be a non-negative integer");
  }

  const double s1 = legs[0].spinX2 / 2.0;
  const double s2 = legs[1].spinX2 / 2.0;
  if (std::abs(s1 - 0.5) > 1e-12 || std::abs(s2 - 0.5) > 1e-12) {
    throw std::invalid_argument(
        "MHelicityConfig::ForwardHadronSourceHelicityStructure: "
        "fixed Regge source expects spin-half beam legs");
  }

  HELMatrix hc;
  hc.BR         = 1.0;
  hc.zeta       = 0.0;
  hc.P_symmetry = false;
  hc.C_symmetry = false;
  hc.alpha_ls.Clear();
  const double J = p.spinX2 / 2.0;
  gra::spin::InitTwoBodyBasis(hc, J, s1, s2,
                              gra::spin::SpinProjections(J),
                              gra::spin::SpinProjections(s1),
                              gra::spin::SpinProjections(s2),
                              "forward hadron source");
  return hc;
}

HelicityChannelMatch FindHelicityChannelMatch(const json &resonance_block, const std::vector<int> &refdecay,
                                              const std::string &table_name, const std::string &pdg_str,
                                              bool allow_charge_conjugate_match, bool allow_leg_exchange_match,
                                              bool production_pair_schema, bool analytic_production_schema,
                                              const std::vector<MParticle> *requested_legs) {
  HelicityChannelMatch match;
  HelicityChannelMatch exchanged_match;

  if (production_pair_schema) {
    if (requested_legs == nullptr) {
      throw std::invalid_argument(
          "MHelicityConfig::ProcessHelicityStructure: "
          "central continuum card pair schema "
          "needs leg metadata");
    }

    const std::vector<int> requested_pair   = ProductionPairKey(refdecay);
    const std::string      requested_sector = ProductionSector(refdecay, *requested_legs);

    std::vector<std::pair<std::vector<int>, std::string>> seen_sectors;
    for (auto it = resonance_block.begin(); it != resonance_block.end(); ++it) {
      const std::string      pair_id      = it.key();
      const auto            &pair_block   = it.value();
      const std::string      pair_context = ChannelContext(table_name, pdg_str, pair_id);
      const std::vector<int> pair =
          ValidateProductionPairBlock(pair_id, pair_block, pair_context, analytic_production_schema);

      for (auto sector_it = pair_block.begin(); sector_it != pair_block.end(); ++sector_it) {
        const std::string sector_id = sector_it.key();
        if (!IsProductionSectorName(sector_id)) { continue; }

        for (const auto &seen : seen_sectors) {
          if (seen.second == sector_id && SamePDGVector(seen.first, pair)) {
            throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: duplicate " + table_name +
                                        " pair/sector entries for PDG " + pdg_str + " pair " + FormatPDGList(pair) +
                                        " sector " + sector_id);
          }
        }
        seen_sectors.push_back({pair, sector_id});
      }

      const bool       exact_pair    = SamePDGVector(requested_pair, pair);
      std::vector<int> reversed_pair = pair;
      if (reversed_pair.size() == 2) { std::swap(reversed_pair[0], reversed_pair[1]); }
      const bool exchanged_pair = allow_leg_exchange_match && SamePDGVector(requested_pair, reversed_pair);
      if (!exact_pair && !exchanged_pair) { continue; }

      if (!pair_block.contains(requested_sector)) { continue; }
      if (match.found) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: ambiguous " + table_name +
                                    " pair entries for PDG " + pdg_str + " requested PDG " + FormatPDGList(refdecay) +
                                    " sector " + requested_sector);
      }

      const std::string channel_id = pair_id + "/" + requested_sector;
      ValidateProductionSectorBlock(pair_block[requested_sector], ChannelContext(table_name, pdg_str, channel_id),
                                    analytic_production_schema);

      match.found                         = true;
      match.channel_id                    = channel_id;
      match.matched_decay                 = pair;
      match.matched_by_leg_exchange       = exchanged_pair && !exact_pair;
      match.matched_by_charge_conjugation = false;
      match.channel_block                 = &pair_block[requested_sector];
    }

    return match;
  }

  for (auto it = resonance_block.begin(); it != resonance_block.end(); ++it) {
    const std::string channel_id    = it.key();
    const auto       &channel_block = it.value();

    const std::string context = ChannelContext(table_name, pdg_str, channel_id);
    if (!channel_block.is_object()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " decay channel is not a JSON object");
    }
    const std::vector<int> decay = ParsePDGKey(channel_id, context);

    if (decay.size() != refdecay.size()) { continue; }
    if (std::equal(refdecay.begin(), refdecay.end(), decay.begin())) {
      match.found                   = true;
      match.channel_id              = channel_id;
      match.matched_decay           = decay;
      match.matched_by_leg_exchange = false;
      match.channel_block           = &channel_block;
      return match;
    }

    const bool valid_charge_conjugate = allow_charge_conjugate_match && IsChargeConjugateSameOrder(refdecay, decay);
    const bool valid_exchange =
        !valid_charge_conjugate && allow_leg_exchange_match && IsTwoBodyLegExchange(refdecay, decay);
    const bool valid_charge_conjugate_reverse =
        !valid_charge_conjugate && allow_charge_conjugate_match && IsParticleAntiparticleReverse(refdecay, decay);
    if (valid_charge_conjugate || valid_exchange || valid_charge_conjugate_reverse) {
      if (exchanged_match.found) {
        throw std::invalid_argument(
            "MHelicityConfig::ProcessHelicityStructure:"
            " ambiguous exchanged two-body " +
            table_name + " entries for PDG " + pdg_str + " and requested PDG " + FormatPDGList(refdecay));
      }
      exchanged_match.found                         = true;
      exchanged_match.channel_id                    = channel_id;
      exchanged_match.matched_decay                 = decay;
      exchanged_match.matched_by_leg_exchange       = valid_exchange || valid_charge_conjugate_reverse;
      exchanged_match.matched_by_charge_conjugation = valid_charge_conjugate;
      exchanged_match.channel_block                 = &channel_block;
    }
  }

  if (exchanged_match.found) { return exchanged_match; }
  return match;
}

// Apply the two-body LS-basis leg-exchange phase to a canonical vertex
//
void ApplyTwoBodyLegExchangePhase(HELMatrix &hc, const MParticle &leg1, const MParticle &leg2,
                                  const std::vector<int> &requested_decay, const HelicityChannelMatch &match,
                                  const std::string &table_name) {
  const double s1 = leg1.spinX2 / 2.0;
  const double s2 = leg2.spinX2 / 2.0;

  std::cout << rang::fg::yellow
            << "MHelicityConfig::ApplyTwoBodyLegExchangePhase: exchanged legs "
               "requested PDG "
            << FormatPDGList(requested_decay) << " matched " << table_name << " channel " << match.channel_id << " PDG "
            << FormatPDGList(match.matched_decay) << rang::fg::reset << std::endl;

  for (gra::spin::LSTerm &term : hc.alpha_ls) {
    const double s           = 0.5 * static_cast<double>(term.two_s);
    const int    phase_power = static_cast<int>(std::llround(static_cast<double>(term.l) + s1 + s2 - s));
    const int    phase_sign  = phase_power % 2 == 0 ? 1 : -1;

    std::cout << rang::fg::yellow << "  alpha_ls[l=" << term.l << ", s=" << gra::spin::FormatSpinLabel(s)
              << "] leg-exchange phase: (-1)^(l + s1 + s2 - s) = (-1)^(" << term.l << " + "
              << gra::spin::FormatSpinLabel(s1) << " + " << gra::spin::FormatSpinLabel(s2) << " - "
              << gra::spin::FormatSpinLabel(s) << ") = " << phase_sign << rang::fg::reset << std::endl;

    if (phase_sign < 0) { term.coefficient *= -1.0; }
  }
}

// Apply the physical-leg exchange phase to every crossed per-m LS coefficient
void ApplyCrossedLSLegExchangePhase(HELMatrix &hc, const MParticle &leg1, const MParticle &leg2,
                                    const std::string &context) {
  if (hc.exchange_basis != ExchangeBasisType::ReggeHelicity ||
      hc.m_ls.size() != static_cast<std::size_t>(2 * hc.analytic_MMAX + 1)) {
    throw std::invalid_argument(context + ": invalid crossed per-m LS basis");
  }
  const double s1 = leg1.spinX2 / 2.0;
  const double s2 = leg2.spinX2 / 2.0;
  for (const auto &column : indices(hc.m_ls)) {
    const int m = static_cast<int>(column) - hc.analytic_MMAX;
    for (gra::spin::LSTerm &term : hc.m_ls[column]) {
      const double S        = 0.5 * static_cast<double>(term.two_s);
      const int    exponent = static_cast<int>(std::llround(static_cast<double>(term.l) + s1 + s2 - S));
      if (exponent % 2 != 0) { term.coefficient *= -1.0; }
      std::cout << rang::fg::yellow << "  g_ls[L=" << term.l << ", S=" << gra::spin::FormatSpinLabel(S) << ", m=" << m
                << "] leg-exchange phase = " << (exponent % 2 == 0 ? 1 : -1) << rang::fg::reset << std::endl;
    }
  }
}

// Compute the active two-body coupling basis requested by the card
//
std::string ParseTwoBodyCouplingBasis(const json &channel_block, const std::string &context, bool production_mode,
                                      gra::spin::VertexContext vertex_context) {
  if (!channel_block.contains("basis") || !channel_block["basis"].is_string()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " is missing string basis");
  }
  const std::string basis = channel_block["basis"].get<std::string>();
  if (!production_mode) {
    if (basis == "alpha_ls" || basis == "helicity") { return basis; }
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " has an unsupported basis");
  }
  const ReggeVertexRole role = vertex_context == gra::spin::VertexContext::SubTUChannelExchange
                                   ? ReggeVertexRole::Continuum
                                   : ReggeVertexRole::Resonance;
  (void)ParseReggeVertexBasis(basis, role);
  return basis;
}

// Compute the crossed Jacob-Wick parity phase at signature allowed spins
int CrossedParityPhase(const MParticle &mother, const std::vector<MParticle> &legs, const regge::Signature signature,
                       const std::string &context) {
  if (legs.size() != 2 || mother.P * mother.P != 1 || legs[0].P * legs[0].P != 1 || legs[1].P * legs[1].P != 1) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity parity requires P = +/-1 and tau = +/-1");
  }
  const int spin_sum_x2 = legs[0].spinX2 + legs[1].spinX2;
  if (spin_sum_x2 < 0 || spin_sum_x2 % 2 != 0) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity parity requires integer s1+s2");
  }
  const int spin_phase = (spin_sum_x2 / 2) % 2 == 0 ? 1 : -1;
  return mother.P * legs[0].P * legs[1].P * regge::Tau(signature) * spin_phase;
}

// Compute the crossed phase for exchanging the two physical legs
int CrossedLegExchangePhase(const std::vector<MParticle> &legs, const regge::Signature signature,
                            const std::string &context) {
  if (legs.size() != 2) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity leg exchange requires two legs and tau = +/-1");
  }
  const int spin_sum_x2 = legs[0].spinX2 + legs[1].spinX2;
  if (spin_sum_x2 < 0 || spin_sum_x2 % 2 != 0) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity leg exchange requires integer s1+s2");
  }
  const int spin_phase = (spin_sum_x2 / 2) % 2 == 0 ? 1 : -1;
  return regge::Tau(signature) * spin_phase;
}

// Read physical photon helicity couplings and complete their symmetry orbit
void ParsePhotonHelicityRows(HELMatrix &hc, const json &channel_block, const std::string &context,
                             const MParticle &mother, const std::vector<MParticle> &legs,
                             const HelicityChannelMatch &match, regge::Signature signature,
                             gra::spin::VertexContext vertex_context) {
  if (legs.size() != 2 || !channel_block.contains("helicity") || !channel_block.at("helicity").is_array() ||
      channel_block.at("helicity").empty()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity requires two legs and a non-empty helicity array");
  }

  const double s1 = legs[0].spinX2 / 2.0;
  const double s2 = legs[1].spinX2 / 2.0;
  gra::gpom::InitPhotonCrossed(hc, s1, s2, context);
  hc.coupling_basis = gra::CouplingBasis::Helicity;
  const int parity_phase = CrossedParityPhase(mother, legs, signature, context);
  const int exchange_phase = CrossedLegExchangePhase(legs, signature, context);
  const bool identical_legs = gra::spin::PhysicalIdenticalPair(
      legs[0], legs[1], true, gra::spin::VertexContext::SubTUChannelExchange);

  std::map<std::pair<long long, long long>, std::complex<double>> reduced;
  std::set<std::pair<long long, long long>>                       input_orbits;
  const auto add = [&](double lambda1, double lambda2, std::complex<double> value, const std::string &origin) {
    const auto key               = std::make_pair(std::llround(2.0 * lambda1), std::llround(2.0 * lambda2));
    const auto [found, inserted] = reduced.emplace(key, value);
    if (!inserted) {
      const double scale = std::max({1.0, std::abs(found->second), std::abs(value)});
      if (std::abs(found->second - value) > 1.0e-9 * scale) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " crossed_helicity symmetry orbit is inconsistent at " + origin);
      }
    }
  };

  for (const auto &i : indices(channel_block.at("helicity"))) {
    const auto &row = channel_block.at("helicity").at(i);
    if (!row.is_array() || row.size() != 4 || !row.at(0).is_number() || !row.at(1).is_number() ||
        !row.at(2).is_number() || !row.at(3).is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " helicity row should be [lambda1, lambda2, magnitude, phase]");
    }
    double       lambda1 = row.at(0).get<double>();
    double       lambda2 = row.at(1).get<double>();
    const double mag     = row.at(2).get<double>();
    const double phase   = row.at(3).get<double>();
    if (!std::isfinite(lambda1) || !std::isfinite(lambda2) || !std::isfinite(mag) || !std::isfinite(phase) ||
        mag < 0.0) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_helicity magnitude must be finite and nonnegative");
    }
    if (match.matched_by_leg_exchange) { std::swap(lambda1, lambda2); }
    // Validate the original projections before constructing integer orbit keys
    (void)gra::spin::SpinProjectionIndex(lambda1, s1, context + " lambda1");
    (void)gra::spin::SpinProjectionIndex(lambda2, s2, context + " lambda2");
    const auto key = std::make_pair(std::llround(2.0 * lambda1), std::llround(2.0 * lambda2));
    std::vector<std::pair<long long, long long>> orbit = {key, {-key.first, -key.second}};
    if (identical_legs) {
      orbit.push_back({key.second, key.first});
      orbit.push_back({-key.second, -key.first});
    }
    const auto canonical = *std::min_element(orbit.cbegin(), orbit.cend());
    if (!input_orbits.insert(canonical).second) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_helicity rows repeat one symmetry orbit");
    }
    if (vertex_context == gra::spin::VertexContext::SubTUChannelExchange && legs[0].pdg == PDG::PDG_gamma &&
        !math::IsExactEqual(std::abs(lambda1), 1.0)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " physical photon helicity must be transverse");
    }
    const std::complex<double> value = std::polar(mag, phase);
    add(lambda1, lambda2, value, "input row");
    add(-lambda1, -lambda2, static_cast<double>(parity_phase) * value, "parity partner");
    if (identical_legs) {
      add(lambda2, lambda1, static_cast<double>(exchange_phase) * value, "leg exchange partner");
      add(-lambda2, -lambda1, static_cast<double>(parity_phase * exchange_phase) * value,
          "parity and leg exchange partner");
    }
  }

  for (const auto &[key, value] : reduced) {
    const double      lambda1  = 0.5 * static_cast<double>(key.first);
    const double      lambda2  = 0.5 * static_cast<double>(key.second);
    const std::size_t physical = gra::spin::DirectHelicityCoordinateIndex(lambda1, lambda2, hc.s1, hc.s2, context);
    for (const int m : {-1, 1}) {
      const std::size_t analytic   = gra::gpom::AnalyticMIndex(m, hc.analytic_MMAX, context);
      hc.T[physical][analytic]     = value;
      hc.T_set[physical][analytic] = true;
    }
  }
}

// Read sparse angle-free Regge helicity couplings and complete their symmetries
void ParseCrossedReggeHelicityRows(HELMatrix &hc, const json &channel_block, const std::string &context,
                                   const MParticle &mother, const std::vector<MParticle> &legs,
                                   const HelicityChannelMatch &match, const regge::Signature signature,
                                   const int analytic_mmax) {
  constexpr const char *field = "helicity";
  if (legs.size() != 2 || mother.pdg == PDG::PDG_gamma || !channel_block.contains(field) ||
      !channel_block.at(field).is_array() || channel_block.at(field).empty()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_helicity requires a nonphoton trajectory, two legs and a "
                                "non-empty helicity array");
  }

  gra::gpom::InitCrossed(hc, legs[0].spinX2 / 2.0, legs[1].spinX2 / 2.0, analytic_mmax, context);
  hc.coupling_basis = gra::CouplingBasis::Helicity;
  const int  parity_phase   = CrossedParityPhase(mother, legs, signature, context);
  const int  exchange_phase = CrossedLegExchangePhase(legs, signature, context);
  const bool identical_legs = gra::spin::PhysicalIdenticalPair(
      legs[0], legs[1], true, gra::spin::VertexContext::SubTUChannelExchange);
  using Key                 = std::tuple<long long, long long, int>;
  std::map<Key, std::complex<double>> residue;
  std::set<Key>                       input_orbits;

  const auto add = [&](const double lambda1, const double lambda2, const int m, const std::complex<double> value,
                       const std::string &origin) {
    const Key key                = {std::llround(2.0 * lambda1), std::llround(2.0 * lambda2), m};
    const auto [found, inserted] = residue.emplace(key, value);
    if (!inserted) {
      const double scale = std::max({1.0, std::abs(found->second), std::abs(value)});
      if (std::abs(found->second - value) > 1.0e-9 * scale) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " crossed_helicity symmetry orbit is inconsistent at " + origin);
      }
    }
  };

  for (const auto &i : indices(channel_block.at(field))) {
    const auto &row = channel_block.at(field).at(i);
    if (!row.is_array() || row.size() != 5 || !row.at(0).is_number() || !row.at(1).is_number() ||
        !row.at(2).is_number_integer() || !row.at(3).is_number() || !row.at(4).is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " helicity row should be [lambda1, lambda2, m, magnitude, "
                                  "phase]");
    }
    double          lambda1 = row.at(0).get<double>();
    double          lambda2 = row.at(1).get<double>();
    const long long m_input = row.at(2).get<long long>();
    const double    mag     = row.at(3).get<double>();
    const double    phase   = row.at(4).get<double>();
    if (!std::isfinite(lambda1) || !std::isfinite(lambda2) || !std::isfinite(mag) || !std::isfinite(phase) ||
        mag < 0.0 || std::abs(lambda1) > 0.5 * std::numeric_limits<int>::max() ||
        std::abs(lambda2) > 0.5 * std::numeric_limits<int>::max() || m_input < -static_cast<long long>(analytic_mmax) ||
        m_input > static_cast<long long>(analytic_mmax)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " helicity row has invalid magnitude, phase or m outside "
                                  "MMAX");
    }
    const long long lambda1x2 = std::llround(2.0 * lambda1);
    const long long lambda2x2 = std::llround(2.0 * lambda2);
    if (std::abs(2.0 * lambda1 - lambda1x2) > 1.0e-9 || std::abs(2.0 * lambda2 - lambda2x2) > 1.0e-9) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed helicity is not integer or half-integer");
    }
    if (match.matched_by_leg_exchange) { std::swap(lambda1, lambda2); }
    const int        m     = static_cast<int>(m_input);
    const Key        key   = {std::llround(2.0 * lambda1), std::llround(2.0 * lambda2), m};
    std::vector<Key> orbit = {key, {-std::get<0>(key), -std::get<1>(key), -std::get<2>(key)}};
    if (identical_legs) {
      orbit.push_back({std::get<1>(key), std::get<0>(key), m});
      orbit.push_back({-std::get<1>(key), -std::get<0>(key), -std::get<2>(key)});
    }
    const Key canonical = *std::min_element(orbit.cbegin(), orbit.cend());
    if (!input_orbits.insert(canonical).second) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " helicity rows repeat one symmetry orbit");
    }

    const std::complex<double> value = std::polar(mag, phase);
    add(lambda1, lambda2, m, value, "input row");
    add(-lambda1, -lambda2, -m, static_cast<double>(parity_phase) * value, "parity partner");
    if (identical_legs) {
      add(lambda2, lambda1, m, static_cast<double>(exchange_phase) * value, "leg exchange partner");
      add(-lambda2, -lambda1, -m, static_cast<double>(parity_phase * exchange_phase) * value,
          "parity and leg exchange partner");
    }
  }

  for (const auto &[key, value] : residue) {
    const double      lambda1    = 0.5 * static_cast<double>(std::get<0>(key));
    const double      lambda2    = 0.5 * static_cast<double>(std::get<1>(key));
    const std::size_t physical   = gra::spin::DirectHelicityCoordinateIndex(lambda1, lambda2, hc.s1, hc.s2, context);
    const std::size_t analytic   = gra::gpom::AnalyticMIndex(std::get<2>(key), hc.analytic_MMAX, context);
    hc.T[physical][analytic]     = value;
    hc.T_set[physical][analytic] = true;
  }
}

// Read sparse per-m crossed LS couplings and complete their parity pairs
void ParseCrossedReggeLSRows(HELMatrix &hc, const json &channel_block, const std::string &context,
                             const std::string &table_name, const MParticle &mother, const std::vector<MParticle> &legs,
                             const int analytic_mmax) {
  constexpr const char *field = "g_ls";
  if (legs.size() != 2 || mother.pdg == PDG::PDG_gamma || !channel_block.contains(field) ||
      !channel_block.at(field).is_array() || channel_block.at(field).empty()) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " crossed_ls requires a nonphoton trajectory, two legs and a "
                                "non-empty g_ls array");
  }

  gra::gpom::InitCrossed(hc, legs[0].spinX2 / 2.0, legs[1].spinX2 / 2.0, analytic_mmax, context);
  hc.coupling_basis = gra::CouplingBasis::LS;
  const std::size_t nm = static_cast<std::size_t>(2 * analytic_mmax + 1);
  hc.m_ls.assign(nm, {});
  using Key = std::tuple<std::size_t, std::size_t, int>;
  std::set<Key> input_orbits;

  const auto add = [&](const std::size_t l, const std::size_t two_s, const int m,
                       const std::complex<double> coefficient) {
    const std::size_t column = gra::gpom::AnalyticMIndex(m, analytic_mmax, context);
    if (!hc.m_ls[column].Insert(l, two_s, coefficient)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " crossed_ls symmetry orbit has duplicate [L=" + std::to_string(l) + ", S=" +
                                  gra::aux::Spin2XtoString(static_cast<int>(two_s)) + ", m=" + std::to_string(m) + "]");
    }
  };

  for (const auto &i : indices(channel_block.at(field))) {
    const auto &row = channel_block.at(field).at(i);
    if (!row.is_array() || row.size() != 5 || !row.at(0).is_number_integer() || !row.at(1).is_number() ||
        !row.at(2).is_number_integer() || !row.at(3).is_number() || !row.at(4).is_number()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " g_ls row should be [L, S, m, magnitude, phase]");
    }
    const long long l_input     = row.at(0).get<long long>();
    const double    S           = row.at(1).get<double>();
    const long long m_input     = row.at(2).get<long long>();
    const double    mag         = row.at(3).get<double>();
    const double    phase       = row.at(4).get<double>();
    const double    two_s_value = 2.0 * S;
    if (l_input < 0 || l_input > std::numeric_limits<int>::max() / 2 || !std::isfinite(S) || S < 0.0 ||
        two_s_value > static_cast<double>(std::numeric_limits<int>::max()) ||
        std::abs(two_s_value - std::round(two_s_value)) > 1.0e-9 || m_input < -static_cast<long long>(analytic_mmax) ||
        m_input > static_cast<long long>(analytic_mmax) || !std::isfinite(mag) || mag < 0.0 || !std::isfinite(phase)) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " g_ls row has invalid L, S, m, magnitude or phase");
    }

    const std::size_t l     = static_cast<std::size_t>(l_input);
    const std::size_t two_s = static_cast<std::size_t>(std::llround(two_s_value));
    const int         m     = static_cast<int>(m_input);
    const Key         orbit = {l, two_s, std::abs(m)};
    if (!input_orbits.insert(orbit).second) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " g_ls rows repeat one m-reflection orbit");
    }
    const std::complex<double> coefficient = std::polar(mag, phase);
    add(l, two_s, m, coefficient);
    if (m != 0) { add(l, two_s, -m, coefficient); }
    std::cout << rang::fg::green << "LS couplings: " << table_name << ": "
              << "[L=" << l << ", S=" << S << ", m=" << m << ", mag=" << std::abs(coefficient)
              << ", phase=" << std::arg(coefficient) << "]" << rang::fg::reset << std::endl;
  }
}

// Check the generated crossed continuum parity partners
void ValidateCrossedHelicityParity(const HELMatrix &hc, const MParticle &mother, const std::vector<MParticle> &legs,
                                   const regge::Signature signature, const std::string &context) {
  if (!hc.P_symmetry) { return; }
  const int phase = CrossedParityPhase(mother, legs, signature, context);
  for (std::size_t row = 0; row < hc.lambda_values.size_row(); ++row) {
    const double      lambda1 = hc.lambda_values[row][0];
    const double      lambda2 = hc.lambda_values[row][1];
    const std::size_t partner_row =
        gra::spin::DirectHelicityCoordinateIndex(-lambda1, -lambda2, hc.s1, hc.s2, context + " parity");
    for (int m = -hc.analytic_MMAX; m <= hc.analytic_MMAX; ++m) {
      const std::size_t column = gra::gpom::AnalyticMIndex(m, hc.analytic_MMAX, context + " parity");
      if (!hc.T_set[row][column]) { continue; }
      const std::size_t partner_column = gra::gpom::AnalyticMIndex(-m, hc.analytic_MMAX, context + " parity");
      if (!hc.T_set[partner_row][partner_column]) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " crossed_helicity parity partner is missing");
      }
      const double scale = std::max({1.0, std::abs(hc.T[row][column]), std::abs(hc.T[partner_row][partner_column])});
      if (std::abs(hc.T[partner_row][partner_column] - static_cast<double>(phase) * hc.T[row][column]) >
          1.0e-9 * scale) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " crossed_helicity rows violate Jacob-Wick parity");
      }
    }
  }
}

// Read two-body P and C symmetry switches common to all coupling bases
//
void ParseTwoBodySymmetryFlags(HELMatrix &hc, const json &channel_block, const std::string &context) {
  ValidateTwoBodySymmetryFlags(channel_block, context);
  hc.C_symmetry = channel_block["CP"][0].get<bool>();
  hc.P_symmetry = channel_block["CP"][1].get<bool>();
}

// Read direct reduced helicity amplitudes H(lambda1,lambda2) into T
//
void ParseDirectHelicityRows(HELMatrix &hc, const json &channel_block, const std::string &context,
                             const std::string &table_name, const MParticle &mother,
                             const std::vector<MParticle> &vertex_legs, const HelicityChannelMatch &match,
                             bool production_mode, gra::spin::VertexContext vertex_context,
                             const std::vector<std::array<double, 2>> *typed_helicity = nullptr,
                             const std::vector<std::complex<double>> *typed_couplings = nullptr) {
  if (vertex_legs.size() != 2) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " direct helicity basis expects two legs");
  }
  const bool typed = typed_helicity != nullptr && typed_couplings != nullptr;
  if ((typed_helicity == nullptr) != (typed_couplings == nullptr)) {
    throw std::logic_error("MHelicityConfig::ProcessHelicityStructure: incomplete typed helicity input");
  }
  if (!typed && (!channel_block.contains("helicity") || !channel_block["helicity"].is_array())) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " is missing array helicity");
  }
  const std::size_t row_count = typed ? typed_helicity->size() : channel_block["helicity"].size();
  if (row_count == 0 || (typed && row_count != typed_couplings->size())) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " helicity array is empty");
  }

  hc.coupling_basis = gra::CouplingBasis::Helicity;
  hc.alpha_ls.Clear();
  const double J  = mother.spinX2 / 2.0;
  const double s1 = vertex_legs[0].spinX2 / 2.0;
  const double s2 = vertex_legs[1].spinX2 / 2.0;
  gra::spin::InitTwoBodyBasis(hc, J, s1, s2,
                              gra::spin::SpinProjections(J),
                              gra::spin::SpinProjections(s1),
                              gra::spin::SpinProjections(s2), context);

  const std::size_t n1 = static_cast<std::size_t>(vertex_legs[0].spinX2 + 1);
  const std::size_t n2 = static_cast<std::size_t>(vertex_legs[1].spinX2 + 1);
  hc.T     = MMatrix<std::complex<double>>(static_cast<unsigned int>(n1), static_cast<unsigned int>(n2), 0.0);
  hc.T_set = MMatrix<bool>(static_cast<unsigned int>(n1), static_cast<unsigned int>(n2), false);

  const double leg_exchange_phase = match.matched_by_leg_exchange ? gra::spin::DirectHelicityLegExchangePhase(
                                                                        mother, vertex_legs[0], vertex_legs[1], context)
                                                                  : 1.0;
  const bool   identical_legs =
      gra::spin::PhysicalIdenticalPair(vertex_legs[0], vertex_legs[1], production_mode, vertex_context);

  // Compute the flattened direct-helicity coordinate index
  const auto direct_index = [&](double lambda1, double lambda2) {
    return gra::spin::DirectHelicityCoordinateIndex(lambda1, lambda2, hc.s1, hc.s2, context + " compact");
  };

  const std::vector<std::vector<std::complex<double>>> raw_basis_vectors = gra::spin::BuildJWDirectHelicityRawBasis(
      hc, mother, vertex_legs[0], vertex_legs[1], production_mode, table_name, context, vertex_context);
  // Insert one direct-helicity amplitude after checking compact-basis conflicts
  const auto insert_helicity = [&](double lambda1, double lambda2, std::complex<double> value,
                                   const std::string &origin) {
    if (std::abs(lambda1 - lambda2) > hc.J + 1e-9) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context + " " + origin +
                                  " violates |lambda1-lambda2| <= J");
    }

    const std::size_t i1 = gra::spin::SpinProjectionIndex(lambda1, hc.s1, context + " " + origin);
    const std::size_t i2 = gra::spin::SpinProjectionIndex(lambda2, hc.s2, context + " " + origin);
    if (hc.T_set[i1][i2]) {
      const double scale = std::max({1.0, std::abs(hc.T[i1][i2]), std::abs(value)});
      if (std::abs(hc.T[i1][i2] - value) > 1e-9 * scale) {
        throw std::invalid_argument(
            "MHelicityConfig::ProcessHelicityStructure: " + context +
            " compact helicity expansion conflicts at [lambda1=" + gra::spin::FormatSpinLabel(lambda1) +
            ", lambda2=" + gra::spin::FormatSpinLabel(lambda2) + "]");
      }
      return;
    }

    hc.T[i1][i2]     = value;
    hc.T_set[i1][i2] = true;

    std::cout << rang::fg::green << "Helicity couplings: " << table_name << ": "
              << "[lambda1=" << gra::spin::FormatSpinLabel(lambda1)
              << ", lambda2=" << gra::spin::FormatSpinLabel(lambda2) << ", mag=" << std::abs(hc.T[i1][i2])
              << ", phase=" << std::arg(hc.T[i1][i2]) << (origin == "card" ? "" : ", generated=" + origin) << "]"
              << rang::fg::reset << std::endl;
  };

  for (std::size_t a = 0; a < row_count; ++a) {
    double               lambda1 = 0.0;
    double               lambda2 = 0.0;
    std::complex<double> coupling = 0.0;
    if (typed) {
      lambda1 = typed_helicity->at(a)[0];
      lambda2 = typed_helicity->at(a)[1];
      coupling = typed_couplings->at(a);
    } else {
      const auto &row = channel_block["helicity"][a];
      if (!row.is_array() || row.size() != 4 || !std::all_of(row.begin(), row.end(), [](const json &x) {
            return x.is_number();
          })) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " has malformed helicity row " + std::to_string(a) +
                                    " (expected numeric [lambda1, lambda2, mag, phase])");
      }
      lambda1          = row[0].get<double>();
      lambda2          = row[1].get<double>();
      const double mag = row[2].get<double>();
      const double phi = row[3].get<double>();
      if (!std::isfinite(mag) || !std::isfinite(phi) || mag < 0.0) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                    " has invalid helicity row " + std::to_string(a));
      }
      coupling = std::polar(mag, phi);
    }
    if (!std::isfinite(lambda1) || !std::isfinite(lambda2) || !std::isfinite(coupling.real()) ||
        !std::isfinite(coupling.imag()) || std::abs(lambda1) > 0.5 * std::numeric_limits<int>::max() ||
        std::abs(lambda2) > 0.5 * std::numeric_limits<int>::max()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " has invalid helicity row " + std::to_string(a));
    }

    const long long lambda1x2 = std::llround(2.0 * lambda1);
    const long long lambda2x2 = std::llround(2.0 * lambda2);
    if (std::abs(2.0 * lambda1 - lambda1x2) > 1e-9 || std::abs(2.0 * lambda2 - lambda2x2) > 1e-9) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " helicity is not integer or half-integer");
    }
    if (match.matched_by_leg_exchange) { std::swap(lambda1, lambda2); }
    std::vector<std::tuple<double, double, std::complex<double>, std::string>> orbit;
    orbit.push_back({lambda1, lambda2, leg_exchange_phase * coupling, "card"});

    // Visit appended partners until parity and leg exchange close the orbit
    for (std::size_t k = 0; k < orbit.size(); ++k) {
      const double               l1    = std::get<0>(orbit[k]);
      const double               l2    = std::get<1>(orbit[k]);
      const std::complex<double> value = std::get<2>(orbit[k]);

      const auto add_orbit = [&orbit](double new_l1, double new_l2, std::complex<double> new_value,
                                      const std::string &origin) {
        for (const auto &entry : orbit) {
          if (std::abs(std::get<0>(entry) - new_l1) < 1e-12 && std::abs(std::get<1>(entry) - new_l2) < 1e-12) {
            const double scale = std::max({1.0, std::abs(std::get<2>(entry)), std::abs(new_value)});
            if (std::abs(std::get<2>(entry) - new_value) > 1e-9 * scale) {
              throw std::invalid_argument(
                  "MHelicityConfig::ProcessHelicityStructure: compact helicity "
                  "orbit has "
                  "inconsistent " +
                  origin + " constraint");
            }
            return;
          }
        }
        orbit.push_back({new_l1, new_l2, new_value, origin});
      };

      if (hc.P_symmetry) {
        const gra::spin::DirectHelicityPhase parity_phase =
            gra::spin::DirectHelicityCoordinatePhase(raw_basis_vectors, direct_index(l1, l2), direct_index(-l1, -l2));
        if (parity_phase.found) { add_orbit(-l1, -l2, parity_phase.phase * value, "parity"); }
      }
      if (identical_legs) {
        const gra::spin::DirectHelicityPhase exchange_phase =
            gra::spin::DirectHelicityCoordinatePhase(raw_basis_vectors, direct_index(l1, l2), direct_index(l2, l1));
        if (exchange_phase.found) { add_orbit(l2, l1, exchange_phase.phase * value, "leg-exchange"); }
      }
    }

    for (const auto &entry : orbit) {
      insert_helicity(std::get<0>(entry), std::get<1>(entry), std::get<2>(entry), std::get<3>(entry));
    }
  }
}

// Build one direct two-body helicity matrix from typed couplings
HELMatrix BuildDirectHelicityCoupling(const MParticle &mother, const std::vector<MParticle> &legs,
                                      const std::vector<std::array<double, 2>> &helicity,
                                      const std::vector<std::complex<double>> &couplings, bool C_symmetry,
                                      bool P_symmetry, bool leg_exchange, const std::string &source_name,
                                      bool verbose_output,
                                      const std::string &verbose_label, gra::spin::VertexContext context) {
  if (helicity.size() != couplings.size()) {
    throw std::invalid_argument(source_name + " has mismatched helicity coordinates and couplings");
  }
  HELMatrix hc;
  hc.C_symmetry = C_symmetry;
  hc.P_symmetry = P_symmetry;
  HelicityChannelMatch match;
  match.matched_by_leg_exchange = leg_exchange;
  ParseDirectHelicityRows(hc, {}, source_name, source_name, mother, legs, match, true, context, &helicity, &couplings);
  gra::spin::ValidateDirectTMatrix(hc, mother, legs[0], legs[1], true, source_name, false, verbose_output,
                                   verbose_label, context);
  return hc;
}

// Convert one automatic coupling name to its spin-algebra selector
//
gra::spin::AutoCentralCouplingMode AutomaticCouplingMode(ReggeVertexBasis basis) {
  if (basis == ReggeVertexBasis::AutoMinL) { return gra::spin::AutoCentralCouplingMode::MinL; }
  if (basis == ReggeVertexBasis::AutoMinS) { return gra::spin::AutoCentralCouplingMode::MinS; }
  if (basis == ReggeVertexBasis::AutoEqualLS) { return gra::spin::AutoCentralCouplingMode::FlatLS; }
  if (basis == ReggeVertexBasis::AutoEqualHelicity) { return gra::spin::AutoCentralCouplingMode::FlatHelicity; }
  throw std::invalid_argument(
      "MHelicityConfig::ProcessHelicityStructure: automatic basis should "
      "select an automatic production mode");
}

// Read one two-body coupling block in LS or direct helicity basis
//
void ParseTwoBodyCouplings(HELMatrix &hc, const json &channel_block, const std::string &context,
                           const std::string &table_name, const MParticle &mother,
                           const std::vector<MParticle> &vertex_legs, const HelicityChannelMatch &match,
                           bool production_mode, gra::spin::VertexContext vertex_context,
                           regge::Signature trajectory_signature, bool analytic_production, int analytic_mmax) {
  ParseTwoBodySymmetryFlags(hc, channel_block, context);

  const std::string basis_name = ParseTwoBodyCouplingBasis(channel_block, context, production_mode, vertex_context);
  if (!production_mode) {
    if (basis_name == "helicity") {
      ParseDirectHelicityRows(hc, channel_block, context, table_name, mother, vertex_legs, match, false,
                              vertex_context);
      return;
    }
    hc.coupling_basis = gra::CouplingBasis::LS;
    ParseAlphaLSRows(hc, channel_block, context, table_name, false);
    return;
  }

  const ReggeVertexRole  role                = vertex_context == gra::spin::VertexContext::SubTUChannelExchange
                                                   ? ReggeVertexRole::Continuum
                                                   : ReggeVertexRole::Resonance;
  const ReggeVertexBasis basis               = ParseReggeVertexBasis(basis_name, role);
  const bool             analytic_trajectory = analytic_production && role == ReggeVertexRole::Continuum &&
                                   (mother.spinX2 == gra::aux::kNullSpinX2 || mother.pdg == PDG::PDG_gamma);
  if (analytic_trajectory) {
    if (!hc.C_symmetry || !hc.P_symmetry) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " analytic Regge production requires CP = [true, true]");
    }
    if (basis == ReggeVertexBasis::Helicity) {
      if (mother.pdg == PDG::PDG_gamma) {
        ParsePhotonHelicityRows(hc, channel_block, context, mother, vertex_legs, match, trajectory_signature,
                                vertex_context);
      } else {
        ParseCrossedReggeHelicityRows(hc, channel_block, context, mother, vertex_legs, match, trajectory_signature,
                                      analytic_mmax);
      }
      return;
    }
    if (basis != ReggeVertexBasis::LS) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                  " uses the wrong analytic production coupling basis");
    }
    if (mother.pdg == PDG::PDG_gamma) {
      gra::gpom::InitPhotonCrossed(hc, vertex_legs[0].spinX2 / 2.0, vertex_legs[1].spinX2 / 2.0, context);
      hc.coupling_basis = gra::CouplingBasis::LS;
      ParseAlphaLSRows(hc, channel_block, context, table_name, true);
    } else {
      ParseCrossedReggeLSRows(hc, channel_block, context, table_name, mother, vertex_legs, analytic_mmax);
    }
    ValidateCrossedLSRows(hc, mother, vertex_legs, production_mode, context, vertex_context);
    return;
  }
  if (mother.spinX2 == gra::aux::kNullSpinX2) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " direct production couplings require a physical pole spin in " + table_name);
  }
  if (basis == ReggeVertexBasis::Helicity) {
    ParseDirectHelicityRows(hc, channel_block, context, table_name, mother, vertex_legs, match, production_mode,
                            vertex_context);
    return;
  }

  hc.coupling_basis = gra::CouplingBasis::LS;
  if (basis != ReggeVertexBasis::LS) {
    throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                " uses the wrong production coupling basis");
  }
  ParseAlphaLSRows(hc, channel_block, context, table_name, production_mode);
}

// Remove validated helicity couplings at or below one magnitude threshold
void PruneHelicityCouplings(HELMatrix &hc, double coupling_min, const std::string &context) {
  if (!std::isfinite(coupling_min) || coupling_min < 0.0) {
    throw std::invalid_argument(context + ": coupling threshold must be nonnegative");
  }
  if (!hc.T.IsFinite() || hc.T_set.size_row() != hc.T.size_row() || hc.T_set.size_col() != hc.T.size_col()) {
    throw std::invalid_argument(context + ": helicity coupling matrix is invalid");
  }
  hc.T_active.clear();
  for (std::size_t row = 0; row < hc.T.size_row(); ++row) {
    for (std::size_t column = 0; column < hc.T.size_col(); ++column) {
      if (!hc.T_set[row][column]) { continue; }
      if (std::abs(hc.T[row][column]) <= coupling_min) {
        hc.T[row][column]     = 0.0;
        hc.T_set[row][column] = false;
        continue;
      }
      hc.T_active.emplace_back(row, column);
    }
  }
}

// Remove inactive per-m crossed LS couplings without changing their symmetry
void PruneCrossedLSCouplings(HELMatrix &hc, const double coupling_min, const std::string &context) {
  if (!std::isfinite(coupling_min) || coupling_min < 0.0 || hc.exchange_basis != ExchangeBasisType::ReggeHelicity ||
      hc.m_ls.size() != static_cast<std::size_t>(2 * hc.analytic_MMAX + 1)) {
    throw std::invalid_argument(context + ": invalid per-m LS basis");
  }
  bool active = false;
  for (spin::LSCoefficients &terms : hc.m_ls) {
    terms.RemoveBelow(coupling_min);
    active = active || !terms.Empty();
  }
  if (!active) { throw std::invalid_argument(context + ": no active g_ls coupling"); }
}

// Compute the two-body statistical factor for identical decay daughters
//
double TwoBodyDecaySymmetryFactor(const std::vector<MParticle> &legs) {
  if (legs.size() != 2) { return 1.0; }
  return (legs[0].pdg == legs[1].pdg) ? 2.0 : 1.0;
}

// Compute the effective decay coupling from branching ratio and two-body phase
// space
//
void FinalizeDecayCoupling(HELMatrix &hc, const MParticle &p, const std::vector<MParticle> &legs,
                           const std::string &decaymode, double zero_eps) {
  bool TP_computed     = false;
  bool g_decay_from_BR = false;

  if (legs.size() == 2) {
    const double amp2 = 1.0;
    const double sym  = TwoBodyDecaySymmetryFactor(legs);

    const bool scalar_pair = p.spinX2 <= 6 && legs[0].spinX2 == 0 && legs[1].spinX2 == 0 &&
                             std::abs(legs[0].mass - legs[1].mass) < zero_eps;
    double mass = p.mass;
    double PS = gra::kinematics::PDW2body(pow2(mass), pow2(legs[0].mass), pow2(legs[1].mass), amp2, sym);
    // Keep the BR normalization independent of an explicit Tensor decay coupling
    if (scalar_pair && std::abs(PS) < zero_eps) {
      for (std::size_t i = 1; i <= 3; ++i) {
        mass = p.mass + i * p.width;
        PS = gra::kinematics::PDW2body(pow2(mass), pow2(legs[0].mass), pow2(legs[1].mass), amp2, sym);
        if (PS > zero_eps) { break; }
      }
    }
    if (hc.g_decay_TP.empty()) {
      const bool gamma_gamma = p.spinX2 == 0 && p.P == -1 && legs[0].pdg == 22 && legs[1].pdg == 22;
      if (gamma_gamma) {
        hc.g_decay_TP = {MTensorPomeron::GDecayPseudoscalarGammaGamma(p.mass, p.width, hc.BR)};
        TP_computed = true;
      } else if (scalar_pair) {
        hc.g_decay_TP = {MTensorPomeron::GDecay(p.spinX2 / 2, mass, p.width, legs[0].mass, hc.BR, sym)};
        TP_computed = true;
      }
    }

    const double g = msqrt((p.spinX2 + 1) * hc.BR * p.width / PS);
    if (std::isnan(g) || std::isinf(g)) {
      throw std::invalid_argument(
          "MHelicityConfig::ProcessHelicityStructure: "
          "Kinematic coupling problem to '" +
          decaymode + "' with resonance PDG " + p.name + "(" + std::to_string(p.pdg) +
          ") (legs too heavy?) (check RES, DECAYMODE "
          "and DECAYS tables)");
    }

    hc.g_decay      = g;
    g_decay_from_BR = true;
  } else {
    hc.g_decay = 1.0;
  }

  hc.g_decay = hc.g_decay * std::exp(math::zi * hc.zeta);

  printf("(Mass, Full width):                    (%0.3E, %0.3E GeV) \n", p.mass, p.width);
  printf("Branching ratio (BR) || Partial width: %0.3E || %0.3E GeV \n", hc.BR, hc.BR * p.width);
  printf("Phase parameter (zeta):                %0.3E [rad] \n", hc.zeta);
  printf("=> Computed decay coupling             %0.3E x exp(i x %0.1f) \n", std::abs(hc.g_decay),
         std::arg(hc.g_decay));
  if (g_decay_from_BR) {
    printf(
        "=> Formula:                             g_decay = sqrt((p.spinX2 + "
        "1) * BR * width / "
        "PS2(S))\n");
  }
  if (!g_decay_from_BR) {
    printf(
        "=> Multi-body decay fallback:          g_decay set to 1 x exp(i x "
        "zeta); not computed "
        "from BR \n");
  }

  printf("=> Tensor Pomeron decay couplings:     [ ");
  for (const auto &i : indices(hc.g_decay_TP)) { printf("%0.3E ", hc.g_decay_TP[i]); }
  printf("] %s \n", TP_computed ? "(computed from J, width and BR)" : "");
}

// Process helicity amplitude by
// reading central continuum couplings or decay couplings from model cards
//
HELMatrix MHelicityConfig::ProcessHelicityStructure(const MParticle &p, const std::vector<MParticle> &legs,
                                                    const LORENTZSCALAR &lts, const MSubProc &subprocess,
                                                    const MModelTunePtr &model_tune, const std::string &decay_mode,
                                                    bool production_mode, bool strict_mode,
                                                    const std::string &verbose_label, bool verbose_output,
                                                    gra::spin::VertexContext context) const {
  if (!IsLoaded()) { throw std::logic_error("MHelicityConfig::ProcessHelicityStructure: cards are not loaded"); }
  HELMatrix  hc;
  const bool tensor_axial_decay =
      !production_mode && subprocess.ISTATE == "TP" && subprocess.CHANNEL == "RES" && !lts.process.RESONANCES.empty() &&
      std::all_of(lts.process.RESONANCES.begin(), lts.process.RESONANCES.end(), [](const auto &entry) {
        return entry.second.p.spinX2 == 2 && entry.second.p.P == 1 && entry.second.p.C == 1;
      });
  const bool                   uses_jw_helicity_algebra = subprocess.UsesJWHelicityAlgebra() || tensor_axial_decay;
  const bool                   verbose_jw               = uses_jw_helicity_algebra && verbose_output;
  const std::vector<MParticle> vertex_legs              = JacobWickVertexLegs(legs, production_mode, context, lts.PDG);
  bool                         production_matrix_ready  = false;

  // Select continuum sectors from physical legs before any JW crossing
  const std::vector<MParticle> &lookup_legs = production_mode ? legs : vertex_legs;
  const std::vector<int>        refdecay    = ParticlePDGList(lookup_legs);

  // Select the model-specific continuum table or the common decay table
  std::string table_name = "DECAYS.json";
  const json *table      = decays.get();
  if (production_mode) {
    const auto found_model = continuum.find(subprocess.ISTATE);
    if (found_model == continuum.end() || found_model->second == nullptr) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: model " + subprocess.ISTATE +
                                  " has no central continuum steering card");
    }
    table_name = "CON_" + subprocess.ISTATE + ".json";
    table      = found_model->second;
  }
  const json &j = *table;

  // --------------------------------------------------------------------

  const int pdg           = p.pdg;
  MParticle matrix_mother = p;

  // Find out if this production or decay object is present in the steering
  // table
  std::string PDG_STR = std::to_string(pdg);

  // Did we find the PDG block in the requested steering table
  bool found = j.count(PDG_STR);

  if (found == true) {
    if (!j[PDG_STR].is_object()) {
      throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + table_name + " entry for PDG " +
                                  PDG_STR + " is not a JSON object");
    }

    // Do we find it from PDG database
    std::string resonance_name = "not-found-from-PDG";
    try {
      MParticle part = lts.PDG.FindByPDG(pdg);
      resonance_name = part.name;
    } catch (...) {
      // continue;
    }
    if (verbose_jw) {
      std::cout << "MHelicityConfig::ProcessHelicityStructure:" << std::endl;
      std::cout << "Process: PDG = " + std::to_string(pdg) + " " + "[" + resonance_name << "]";

      if (production_mode) {
        std::cout << " <- ";
      } else {
        std::cout << " -> ";
      }

      for (const auto &i : indices(vertex_legs)) {
        std::cout << vertex_legs[i].pdg << " ";

        std::string resonance_name = "not-found-from-PDG";
        try {
          MParticle part = lts.PDG.FindByPDG(vertex_legs[i].pdg);
          resonance_name = part.name;
        } catch (...) {
          // continue;
        }
        std::cout << "[" << resonance_name << "]"
                  << " ";
      }
      std::cout << std::endl;
    }

    const bool           allow_charge_conjugate_match = production_mode && (p.C != 0);
    const bool           allow_leg_exchange_match     = AllowTwoBodyLegExchange(production_mode, context, refdecay);
    const bool           production_pair_schema       = production_mode;
    HelicityChannelMatch match                        = FindHelicityChannelMatch(
                               j[PDG_STR], refdecay, table_name, PDG_STR, allow_charge_conjugate_match, allow_leg_exchange_match,
                               production_pair_schema, subprocess.ISTATE == "GP", &lookup_legs);
    if (!match.found) {
      if (production_mode || strict_mode || subprocess.ISTATE == "MP" || subprocess.ISTATE == "XP" ||
          subprocess.ISTATE == "GP" || subprocess.ISTATE == "TP") {
        std::string str = "MHelicityConfig::ProcessHelicityStructure: Missing strict " + table_name + " entry for " +
                          std::to_string(pdg) + (production_mode ? " < " : " > ");
        for (const auto &i : indices(vertex_legs)) { str += std::to_string(vertex_legs[i].pdg) + " "; }
        throw MissingHelicityData(str);
      }
      gra::aux::PrintWarning();
      std::cout << rang::fg::red
                << "WARNING: " + table_name + " contains no information on this process: " + std::to_string(pdg) +
                       " > ";

      for (const auto &i : indices(vertex_legs)) { std::cout << vertex_legs[i].pdg << " "; }
      std::cout << std::endl;

      const bool select_lowest_decay_ls = uses_jw_helicity_algebra && vertex_legs.size() == 2;
      if (select_lowest_decay_ls) {
        std::cout << "(setting up BR = 1.0, lowest allowed alpha_ls = 1.0, "
                     "P_symmetry = true)";
      } else {
        std::cout << "(setting up BR = 1.0 and default decay-parameter fallback)";
      }
      std::cout << rang::fg::reset << std::endl;
      hc.BR         = 1.0;
      hc.zeta       = 0.0;
      hc.P_symmetry = true;
      hc.C_symmetry = false;
      if (select_lowest_decay_ls) {
        hc.alpha_ls.Clear();
        const auto allowed_ls = gra::spin::AllowedLSCouplings(hc, matrix_mother, vertex_legs[0], vertex_legs[1], false,
                                                              table_name + " missing-decay fallback", context);
        if (allowed_ls.empty()) {
          throw std::invalid_argument(
              "MHelicityConfig::ProcessHelicityStructure: No allowed fallback "
              "alpha_ls coupling for decay PDG " +
              std::to_string(pdg));
        }
        const auto lowest_ls = allowed_ls.front();
        hc.alpha_ls.Set(lowest_ls.l, lowest_ls.two_s, 1.0);
      } else if (!uses_jw_helicity_algebra) {
        hc.alpha_ls.Clear();
      } else {
        hc.alpha_ls.Clear();
      }
    } else {
      const gra::spin::VertexContext vertex_context = context;
      const std::string              context        = ChannelContext(table_name, PDG_STR, match.channel_id);
      MParticle                      vertex_mother  = p;

      ParseDecayParameters(hc, *match.channel_block, context, production_mode, subprocess.ISTATE);
      if (uses_jw_helicity_algebra) {
        const std::string basis = match.channel_block->contains("basis") && match.channel_block->at("basis").is_string()
                                      ? match.channel_block->at("basis").get<std::string>()
                                      : "";
        if (production_mode && IsAutomaticCoupling(basis, ReggeVertexRole::Continuum)) {
          if (!production_mode || table_name.rfind("CON_", 0) != 0) {
            throw std::invalid_argument(
                "MHelicityConfig::ProcessHelicityStructure: automatic "
                "couplings are supported only by continuum cards");
          }
          ParseTwoBodySymmetryFlags(hc, *match.channel_block, context);
          const ReggeVertexBasis parsed = ParseReggeVertexBasis(basis, ReggeVertexRole::Continuum);
          hc = gra::spin::BuildAutomaticCentralCoupling(vertex_mother, vertex_legs, AutomaticCouplingMode(parsed), true,
                                                        table_name + " " + context, verbose_label, verbose_jw,
                                                        hc.P_symmetry, hc.C_symmetry, vertex_context);
          const auto                &g        = match.channel_block->at("g");
          const std::complex<double> coupling = std::polar(g.at(0).get<double>(), g.at(1).get<double>());
          hc.alpha_ls.Scale(coupling);
          hc.T *= coupling;
          production_matrix_ready = true;
        } else {
          regge::Signature analytic_signature = regge::Signature::Negative;
          if (vertex_mother.spinX2 == gra::aux::kNullSpinX2 && vertex_mother.pdg != PDG::PDG_gamma) {
            const auto found = trajectory_signature.find(vertex_mother.pdg);
            if (found == trajectory_signature.end()) {
              throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                          " has no trajectory signature mapping");
            }
            analytic_signature = found->second;
          }
          ParseTwoBodyCouplings(hc, *match.channel_block, context, table_name, vertex_mother, vertex_legs, match,
                                production_mode, vertex_context, analytic_signature, subprocess.ISTATE == "GP",
                                lts.process.MMAX);
          if (hc.UsesReggeDomain() && hc.UsesLSCouplings()) {
            MParticle pole = vertex_mother;
            if (pole.pdg == PDG::PDG_gamma) {
              pole.spinX2 = 2;
            } else {
              const auto found = trajectory_pole_spinX2.find(pole.pdg);
              if (found == trajectory_pole_spinX2.end()) {
                throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                            " has no canonical trajectory pole spin");
              }
              pole.spinX2 = found->second;
            }
            const auto validate_terms = [&](const spin::LSCoefficients &ls) {
              const std::vector<spin::LSTerm> terms(ls.cbegin(), ls.cend());
              spin::ValidateCanonicalPoleTerms(pole, vertex_legs[0], vertex_legs[1], terms, production_mode,
                                               hc.C_symmetry, hc.P_symmetry, vertex_context);
            };
            if (hc.exchange_basis == ExchangeBasisType::ReggeHelicity) {
              for (const spin::LSCoefficients &terms : hc.m_ls) {
                if (!terms.Empty()) { validate_terms(terms); }
              }
            } else {
              validate_terms(hc.alpha_ls);
            }
          }
          if (hc.UsesReggeDomain() && hc.UsesHelicityCouplings()) {
            regge::Signature signature = regge::Signature::Negative;
            if (vertex_mother.pdg != PDG::PDG_gamma) {
              const auto found = trajectory_signature.find(vertex_mother.pdg);
              if (found == trajectory_signature.end()) {
                throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: " + context +
                                            " has no trajectory signature mapping");
              }
              signature = found->second;
            }
            if (match.matched_by_leg_exchange) {
              hc.T *= static_cast<double>(CrossedLegExchangePhase(vertex_legs, signature, context));
            }
            ValidateCrossedHelicityParity(hc, vertex_mother, vertex_legs, signature, context);
          }
        }
        if (match.matched_by_leg_exchange && !hc.UsesHelicityCouplings() &&
            (!hc.UsesReggeDomain() || hc.UsesLSCouplings()) && !production_matrix_ready) {
          if (hc.exchange_basis == ExchangeBasisType::ReggeHelicity) {
            ApplyCrossedLSLegExchangePhase(hc, vertex_legs[0], vertex_legs[1], context);
          } else {
            ApplyTwoBodyLegExchangePhase(hc, vertex_legs[0], vertex_legs[1], refdecay, match, table_name);
          }
        }
      } else {
        hc.alpha_ls.Clear();
        hc.m_ls.clear();
        hc.m_orbital.clear();
      }
      matrix_mother = vertex_mother;
    }

    // Drop inactive LS terms before constructing orbital matrices and caches
    if (!hc.UsesHelicityCouplings() && hc.UsesLSCouplings() && hc.exchange_basis == ExchangeBasisType::ReggeHelicity &&
        !production_matrix_ready) {
      if (model_tune == nullptr) {
        throw std::logic_error(
            "MHelicityConfig::ProcessHelicityStructure: model tune is not "
            "configured");
      }
      PruneCrossedLSCouplings(hc, model_tune->Global().coupling_min,
                              "MHelicityConfig::ProcessHelicityStructure: crossed_ls");
    }
    if (!hc.UsesHelicityCouplings() && !hc.alpha_ls.Empty() && !production_matrix_ready) {
      if (model_tune == nullptr) {
        throw std::logic_error(
            "MHelicityConfig::ProcessHelicityStructure: model tune is not "
            "configured");
      }
      if (!hc.UsesReggeDomain() && vertex_legs.size() == 2) {
        ValidateFiniteLSRows(hc, matrix_mother, vertex_legs, production_mode, table_name, "finite LS steering",
                             context);
      }
      hc.alpha_ls.RemoveBelow(model_tune->Global().coupling_min);
      if (hc.alpha_ls.Empty()) {
        throw std::invalid_argument("MHelicityConfig::ProcessHelicityStructure: no active " +
                                    std::string(production_mode ? "production g_ls" : "decay alpha_ls") + " coupling");
      }
    }
    if (hc.UsesHelicityCouplings() && hc.UsesReggeDomain() && !production_matrix_ready) {
      if (model_tune == nullptr) {
        throw std::logic_error(
            "MHelicityConfig::ProcessHelicityStructure: model tune is not "
            "configured");
      }
      PruneHelicityCouplings(hc, model_tune->Global().coupling_min,
                             "MHelicityConfig::ProcessHelicityStructure: crossed_helicity");
    }

    // ----------------------------------------------------------------------------
    // Init helicity structure matrix

    if (vertex_legs.size() == 2 && uses_jw_helicity_algebra) {
      MParticle named_p = matrix_mother;
      if (named_p.name.empty()) {
        try {
          named_p.name = lts.PDG.FindByPDG(named_p.pdg).name;
        } catch (...) {
          // Keep the resonance-card particle unchanged if no PDG-table name
          // exists
        }
      }
      try {
        if (production_mode && table_name.rfind("CON_", 0) == 0 && !hc.UsesHelicityCouplings() &&
            !hc.UsesReggeDomain() && !production_matrix_ready) {
          gra::spin::ApplyCanonicalPoleLSCoefficients(hc, named_p, vertex_legs[0], vertex_legs[1]);
          gra::spin::InitTMatrix(hc, named_p, vertex_legs[0], vertex_legs[1], true,
                                 table_name + " canonical pole operator", false, false, verbose_label, context);
          production_matrix_ready = true;
        }
        if (!production_matrix_ready && hc.UsesReggeDomain()) {
          if (!hc.UsesLSCouplings()) {
            gra::gpom::CheckHelicity(hc,
                                     "MHelicityConfig::ProcessHelicityStructure: analytic "
                                     "trajectory",
                                     verbose_jw);
            if (named_p.pdg == PDG::PDG_gamma) {
              if (hc.analytic_MMAX != 1) { throw std::invalid_argument("MHelicityConfig: photon residue requires MMAX=1"); }
              for (std::size_t row = 0; row < hc.T.size_row(); ++row) {
                if (!hc.T_set[row][0] || !hc.T_set[row][2]) { continue; }
                const double scale = std::max({1.0, std::abs(hc.T[row][0]), std::abs(hc.T[row][2])});
                if (std::abs(hc.T[row][0] - hc.T[row][2]) > 1.0e-10 * scale) {
                  throw std::invalid_argument("MHelicityConfig: photon residue depends on the external spin projection");
                }
              }
            }
          }
        } else if (!production_matrix_ready && hc.UsesHelicityCouplings()) {
          gra::spin::ValidateDirectTMatrix(hc, named_p, vertex_legs[0], vertex_legs[1], production_mode, table_name,
                                           true, verbose_jw, verbose_label, context);
        } else if (!production_matrix_ready) {
          gra::spin::InitTMatrix(hc, named_p, vertex_legs[0], vertex_legs[1], production_mode, table_name, match.found,
                                 verbose_jw, verbose_label, context);
        }
      } catch (std::invalid_argument &e) {
        std::string detail =
            production_mode ? "Problem with production vertex PDG = " : "Problem with resonance PDG = ";
        detail += std::to_string(pdg);
        if (match.found) {
          detail += " in " + table_name + " entry";
          if (!match.channel_id.empty()) { detail += " channel " + match.channel_id; }
          detail += " with requested PDG " + FormatPDGList(refdecay);
          detail += " and matched table PDG " + FormatPDGList(match.matched_decay);
          if (match.matched_by_leg_exchange) { detail += " (matched by physical leg exchange)"; }
          if (match.matched_by_charge_conjugation) { detail += " (matched by charge conjugation)"; }
        }
        throw std::invalid_argument(detail + " : " + e.what());
      }
    }
    // ----------------------------------------------------------------------------

    if (!production_mode) {
      constexpr double zero_eps = 1e-12;
      if (subprocess.ISTATE == "TP" && (p.pdg == 113 || p.pdg == 223) && legs.size() == 2 &&
          std::abs(legs[0].pdg) == 211 && legs[0].pdg == -legs[1].pdg) {
        const auto param = lts.model_cache ? GetTensorParam(*lts.model_cache, lts.PDG, lts.process.RESONANCES)
                                           : ReadTensorPomeronParam(*model_tune, lts.PDG);
        if (param->FindVector(p.pdg).width_model == TensorVectorWidthModel::RhoOmega) {
          if (!hc.g_decay_TP.empty()) {
            throw std::invalid_argument("RHO_OMEGA pion vertices must be set in PARAM_TENSORPOM.RHO_OMEGA.g");
          }
          hc.g_decay_TP = {param->rho_omega.g[param->rho_omega.Index(p.pdg)]};
        }
      }
      FinalizeDecayCoupling(hc, p, legs, decay_mode, zero_eps);
      if (verbose_jw && lts.process.QMETRICS && lts.process.SPINDEC && vertex_legs.size() == 2 &&
          (subprocess.ISTATE == "MP" || subprocess.ISTATE == "XP" || subprocess.ISTATE == "GP")) {
        spin::PrintDecayEntanglement(hc);
      }
    }

  } else {
    const std::string subject     = production_mode ? "production" : "branching ratio";
    const std::string object_type = production_mode ? "production vertex/exchange" : "resonance";
    std::string       str = "MHelicityConfig::ProcessHelicityStructure: Did not find any " + subject + " data in " +
                      table_name + " for " + object_type + " with PDG: " + std::to_string(p.pdg) +
                      " and requested legs " + FormatPDGList(refdecay);
    throw MissingHelicityData(str);
  }
  aux::PrintBar(".");

  return hc;
}

// Construct one canonical MP or XP pole operator from its continuum card
//
gra::spin::PoleLS MHelicityConfig::ProcessPoleOperatorStructure(ReggeProductionModel model, const MParticle &p,
                                                                const std::vector<MParticle> &legs,
                                                                const LORENTZSCALAR &lts, const MSubProc &subprocess,
                                                                const MModelTunePtr     &model_tune,
                                                                gra::spin::VertexContext context) const {
  if (!IsLoaded()) { throw std::logic_error("MHelicityConfig::ProcessPoleOperatorStructure: cards are not loaded"); }
  if ((model != ReggeProductionModel::MP && model != ReggeProductionModel::XP) ||
      subprocess.ISTATE != ReggeProductionModelName(model) ||
      context != gra::spin::VertexContext::SubTUChannelExchange) {
    throw std::invalid_argument(
        "MHelicityConfig::ProcessPoleOperatorStructure: requires the active "
        "MP or XP continuum sub-t/u model");
  }

  const std::string model_name = ReggeProductionModelName(model);
  const std::string table_name = "CON_" + model_name + ".json";
  const auto       &table      = *continuum.at(model_name);
  const std::string pdg_key    = std::to_string(p.pdg);
  if (!table.contains(pdg_key) || !table.at(pdg_key).is_object()) {
    throw MissingHelicityData("MHelicityConfig::ProcessPoleOperatorStructure: missing " + table_name +
                              " exchange PDG " + pdg_key);
  }

  const std::vector<MParticle> vertex_legs = JacobWickVertexLegs(legs, true, context, lts.PDG);
  const std::vector<int>       requested   = ParticlePDGList(vertex_legs);
  const HelicityChannelMatch   match =
      FindHelicityChannelMatch(table.at(pdg_key), requested, table_name, pdg_key, p.C != 0,
                               AllowTwoBodyLegExchange(true, context, requested), true, false, &vertex_legs);
  if (!match.found || match.channel_block == nullptr) {
    throw MissingHelicityData("MHelicityConfig::ProcessPoleOperatorStructure: missing " + table_name + " entry for " +
                              pdg_key + " < " + FormatPDGList(requested));
  }

  const json &block = *match.channel_block;
  ValidateProductionSectorBlock(block, table_name + " channel", false);
  HELMatrix symmetry;
  ParseTwoBodySymmetryFlags(symmetry, block, table_name + " channel");
  const ReggeVertexBasis basis =
      ParseReggeVertexBasis(block.at("basis").get<std::string>(), ReggeVertexRole::Continuum);
  double                         Lambda = 1.0;
  std::vector<gra::spin::LSTerm> terms;

  if (IsAutomaticReggeVertexBasis(basis)) {
    const auto mode = AutomaticCouplingMode(basis);
    const auto preset =
        gra::spin::BuildAutomaticCentralCoupling(p, vertex_legs, mode, true, table_name, "continuum subvertex", false,
                                                 symmetry.P_symmetry, symmetry.C_symmetry, context);
    if (!preset.UsesHelicityCouplings()) {
      terms.reserve(preset.alpha_ls.Size());
      for (const gra::spin::LSTerm &term : preset.alpha_ls) { terms.push_back(term); }
    } else {
      const auto expansion =
          gra::spin::DirectHelicityToLSCoefficients(preset, p, vertex_legs[0], vertex_legs[1], true, context);
      std::complex<double> scale = 0.0;
      for (const auto &i : indices(expansion.rows)) {
        const auto                &row = expansion.rows[i];
        const std::complex<double> coefficient =
            expansion.coefficients[i] / gra::spin::RawPoleLSNormalization(vertex_legs[0].spinX2, vertex_legs[1].spinX2,
                                                                          p.spinX2, row.l, static_cast<int>(row.two_s));
        if (std::fpclassify(std::abs(scale)) == FP_ZERO && std::abs(coefficient) > 1.0e-12) { scale = coefficient; }
        terms.push_back({row.l, row.two_s, coefficient});
      }
      if (std::abs(scale) <= 1.0e-12) {
        throw std::invalid_argument(
            "MHelicityConfig::ProcessPoleOperatorStructure: automatic "
            "helicity preset is zero");
      }
      for (auto &term : terms) { term.coefficient /= scale; }
    }
    const auto                &g        = block.at("g");
    const std::complex<double> coupling = std::polar(g.at(0).get<double>(), g.at(1).get<double>());
    for (auto &term : terms) { term.coefficient *= coupling; }
  } else if (basis == ReggeVertexBasis::LS) {
    terms = ReadCanonicalPoleTerms(block, "g_ls", table_name);
    gra::spin::ValidateCanonicalPoleTerms(p, vertex_legs[0], vertex_legs[1], terms, true, symmetry.C_symmetry,
                                          symmetry.P_symmetry, context);
    ApplyCanonicalPoleLegExchange(terms, vertex_legs[0], vertex_legs[1], match.matched_by_leg_exchange, table_name);
  } else if (basis == ReggeVertexBasis::Helicity) {
    const json &direct_block = block;
    HELMatrix   direct = symmetry;
    ParseDecayParameters(direct, direct_block, table_name, true, subprocess.ISTATE);
    ParseDirectHelicityRows(direct, direct_block, table_name, table_name, p, vertex_legs, match, true, context);
    gra::spin::ValidateDirectTMatrix(direct, p, vertex_legs[0], vertex_legs[1], true, table_name, true, false, "",
                                     context);
    const auto expansion =
        gra::spin::DirectHelicityToLSCoefficients(direct, p, vertex_legs[0], vertex_legs[1], true, context);
    if (expansion.residual_norm2 > 1.0e-18) {
      throw std::invalid_argument(
          "MHelicityConfig::ProcessPoleOperatorStructure: direct "
          "helicity tensor is outside its canonical pole basis");
    }
    for (const auto &i : indices(expansion.rows)) {
      const auto &row = expansion.rows[i];
      terms.push_back({row.l, row.two_s,
                       expansion.coefficients[i] /
                           gra::spin::RawPoleLSNormalization(vertex_legs[0].spinX2, vertex_legs[1].spinX2, p.spinX2,
                                                             row.l, static_cast<int>(row.two_s))});
    }
  } else {
    throw std::invalid_argument(
        "MHelicityConfig::ProcessPoleOperatorStructure:"
        " unsupported MP or XP basis");
  }

  if (model_tune == nullptr) {
    throw std::logic_error(
        "MHelicityConfig::ProcessPoleOperatorStructure: model tune is not "
        "configured");
  }
  return gra::spin::PreparePoleLS(p, vertex_legs[0], vertex_legs[1], terms, Lambda, true, symmetry.C_symmetry,
                                  symmetry.P_symmetry, context, model_tune->Global().coupling_min);
}

// Construct helicity vertices recursively through one decay branch
void MHelicityConfig::ProcessHelicityTree(MDecayBranch &branch, const LORENTZSCALAR &lts, const MSubProc &subprocess,
                                          const MModelTunePtr &model_tune, const std::string &decay_mode,
                                          bool production_mode, bool strict_mode) const {
  std::vector<MParticle> legs;
  legs.reserve(branch.legs.size());
  for (const auto &leg : branch.legs) { legs.push_back(leg.p); }

  if (production_mode) {
    branch.hel = branch.p.pdg == PDG::PDG_gamma ? ForwardPhotonSourceHelicityStructure(branch.p, legs)
                                                : ForwardHadronSourceHelicityStructure(branch.p, legs);
    if (subprocess.UsesJWHelicityAlgebra()) {
      std::cout << "MHelicityConfig::ProcessHelicityTree: forward source "
                   "metadata is constructed from the exchange and beam spins";
      if (branch.p.pdg == PDG::PDG_gamma) {
        std::cout << ", PARAM_REGGE.PHOTON_VERTEX = " << lts.process.PHOTON_VERTEX;
      }
      std::cout << std::endl;
    }
    return;
  }

  const gra::spin::VertexContext context =
      production_mode ? gra::spin::VertexContext::CrossedBeamLeg : gra::spin::VertexContext::Auto;
  branch.hel = ProcessHelicityStructure(branch.p, legs, lts, subprocess, model_tune, decay_mode, production_mode,
                                        strict_mode, "", true, context);

  for (auto &leg : branch.legs) {
    if (!leg.legs.empty()) {
      ProcessHelicityTree(leg, lts, subprocess, model_tune, decay_mode, production_mode, strict_mode);
    }
  }
}

}  // namespace gra
