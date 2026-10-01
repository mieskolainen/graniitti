// Regge amplitude parameter parsing and immutable storage
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MReggeParam.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <tuple>
#include <utility>

#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Tech/MAux.h"
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

namespace gra {
namespace {

// Build a charge-class order compatible with math::GetAmpPerm(type=0)
//
std::vector<int> ContinuumChargeClassOrder(const gra::LORENTZSCALAR &lts, std::size_t n_central, const std::string &context) {
  if (lts.decaytree.size() < n_central) { throw std::invalid_argument(context + ": decay tree has fewer central particles than requested"); }

  std::map<int, bool> has_antiparticle;
  for (std::size_t i = 0; i < n_central; ++i) {
    const int pdg = lts.decaytree[i].p.pdg;
    if (pdg < 0) { has_antiparticle[std::abs(pdg)] = true; }
  }

  std::map<int, std::size_t> neutral_seen;
  std::vector<int>           plus;
  std::vector<int>           minus;
  plus.reserve(n_central / 2);
  minus.reserve(n_central / 2);
  for (std::size_t i = 0; i < n_central; ++i) {
    const auto &particle = lts.decaytree[i].p;
    const int   pdg      = particle.pdg;
    const int   abs_pdg  = std::abs(pdg);
    int         sign     = 0;
    if (particle.chargeX3 > 0) {
      sign = 1;
    } else if (particle.chargeX3 < 0) {
      sign = -1;
    }
    if (has_antiparticle[abs_pdg]) { sign = (pdg >= 0) ? 1 : -1; }
    if (sign == 0) {
      const std::size_t seen = neutral_seen[pdg]++;
      sign                   = (seen % 2 == 0) ? 1 : -1;
    }

    if (sign > 0) {
      plus.push_back(static_cast<int>(i) + 3);
    } else {
      minus.push_back(static_cast<int>(i) + 3);
    }
  }

  if (plus.size() != minus.size()) {
    throw std::invalid_argument(context +
                                ": charged permutations require balanced "
                                "continuum charge classes");
  }

  std::vector<int> order;
  order.reserve(n_central);
  for (const auto &i : indices(plus)) {
    order.push_back(plus[i]);
    order.push_back(minus[i]);
  }
  return order;
}

// Generate continuum ladder permutations from the shared math utility
//
std::vector<std::vector<int>> ContinuumLadderPermutations(const gra::LORENTZSCALAR &lts, int n_central, const std::string &context, regge::PermType mode) {
  const int effective_mode = regge::PermCount(lts.decaytree, static_cast<std::size_t>(n_central), mode);

  auto permutations = gra::math::GetAmpPerm(n_central, effective_mode);
  if (effective_mode == 1) { return permutations; }

  constexpr int offset = 3;
  const auto    order  = ContinuumChargeClassOrder(lts, static_cast<std::size_t>(n_central), context);
  for (auto &permutation : permutations) {
    for (auto &index : permutation) {
      const std::size_t slot = static_cast<std::size_t>(index - offset);
      if (slot >= order.size()) {
        throw std::invalid_argument(context +
                                    ": library permutation index maps outside "
                                    "continuum charge-class order");
      }
      index = order[slot];
    }
  }
  return permutations;
}

}  // namespace

// Regge amplitude parameters
namespace regge {

// Collect the current decay-tree final-state PDG codes
std::vector<int> FinalPDGs(const LORENTZSCALAR &lts) {
  std::vector<int> pdg;
  pdg.reserve(lts.decaytree.size());
  for (const auto &branch : lts.decaytree) { pdg.push_back(branch.p.pdg); }
  return pdg;
}

// Convert one central amplitude index to a decay-tree slot
std::size_t DecayIndex(const int amplitude_index, const LORENTZSCALAR &lts, const std::string &context) {
  constexpr int offset = 3;
  if (amplitude_index < offset) { throw std::invalid_argument(context + ": invalid central index " + std::to_string(amplitude_index)); }
  const std::size_t index = static_cast<std::size_t>(amplitude_index - offset);
  if (index >= lts.decaytree.size()) { throw std::invalid_argument(context + ": central index outside decay tree"); }
  return index;
}

// Compute the C-induced sign of one parsed continuum vertex
double VertexSign(const VertexParam &row, const MDecayBranch &left, const MDecayBranch &right, const Param &param) { return AntiparticleSign(param, row.first, left.p.pdg) * AntiparticleSign(param, row.second, right.p.pdg); }

// Format a PDG vector for diagnostics
//
std::string FormatPDGVector(const std::vector<int> &pdg) {
  std::ostringstream ss;
  ss << "[";
  for (const auto &i : indices(pdg)) {
    ss << pdg[i];
    if (i + 1 < pdg.size()) { ss << " "; }
  }
  ss << "]";
  return ss.str();
}

// Parse one trajectory PDG alias array from an exchange mapping row
//
std::vector<int> ParsePDGAliasEntry(const nlohmann::json &entry, const std::string &context) {
  if (!entry.is_array() || entry.empty()) { throw std::invalid_argument(context + " must be a non-empty integer PDG array"); }
  std::vector<int> out;
  for (const auto &value : entry) {
    if (!value.is_number_integer()) { throw std::invalid_argument(context + " has a non-integer PDG alias"); }
    out.push_back(value.get<int>());
  }
  return out;
}

// Parse one exact lower-case trajectory role
SoftExchangeRole ParseExchangeParamRole(const nlohmann::json &value, const std::string &context) {
  if (!value.is_string()) { throw std::invalid_argument(context + " must be a string"); }
  const std::string role = value.get<std::string>();
  if (role == "pomeron") { return SoftExchangeRole::Pomeron; }
  if (role == "reggeon") { return SoftExchangeRole::Reggeon; }
  if (role == "odderon") { return SoftExchangeRole::Odderon; }
  throw std::invalid_argument(context + " must be pomeron, reggeon or odderon");
}

// Resolve and validate the C-parity of one exchange alias
//
int CParityFromPDG(const gra::MPDG &particle_data, int pdg, const std::string &context) {
  int c = 0;
  try {
    c = particle_data.FindByPDG(pdg).C;
  } catch (const std::exception &e) { throw std::invalid_argument(context + " cannot resolve exchange PDG " + std::to_string(pdg) + " from the PDG table: " + e.what()); }
  if (c != -1 && c != 1) { throw std::invalid_argument(context + " exchange PDG " + std::to_string(pdg) + " has undefined C-parity in the PDG table"); }
  return c;
}

// Register one alias row and return its common crossing parity
int RegisterExchangeAliases(Param &param, const ExchangeParam &mapping, const std::size_t trajectory, const gra::MPDG &particle_data, const std::string &context) {
  int trajectory_c = 0;
  for (const int alias : mapping.pdg) {
    if (param.PDG_TO_INDEX.count(alias)) { throw std::invalid_argument(context + " duplicates exchange PDG alias " + std::to_string(alias)); }
    const int alias_c = CParityFromPDG(particle_data, alias, context + ".pdg");
    if (trajectory_c == 0) {
      trajectory_c = alias_c;
    } else if (trajectory_c != alias_c) {
      throw std::invalid_argument(context + ".pdg aliases must have a common C-parity");
    }
    param.PDG_TO_INDEX.emplace(alias, trajectory);
    param.PDG_TO_C.emplace(alias, alias_c);
  }
  return trajectory_c;
}

// Compute whether two trajectory aliases carry the same internal quantum numbers
bool SamePoleQuantumNumbers(const MParticle &first, const MParticle &second) { return first.chargeX3 == second.chargeX3 && first.isospinX2 == second.isospinX2 && first.P == second.P && first.C == second.C && first.G == second.G; }

// Resolve the unique canonical pole matching one analytic trajectory alias
const MParticle &ResolvePole(const ExchangeParam &mapping, const MPDG &particle_data, const int pdg, const std::string &context) {
  const MParticle &analytic = particle_data.FindByPDG(pdg);
  if (analytic.spinX2 != gra::aux::kNullSpinX2) { throw std::invalid_argument(context + " exchange PDG " + std::to_string(pdg) + " is not an analytic trajectory alias"); }
  const MParticle *pole = nullptr;
  for (const int alias : mapping.pdg) {
    const MParticle &candidate = particle_data.FindByPDG(alias);
    if (candidate.spinX2 != mapping.pole_spinX2 || !SamePoleQuantumNumbers(analytic, candidate)) { continue; }
    if (pole != nullptr) { throw std::invalid_argument(context + " analytic exchange PDG " + std::to_string(pdg) + " has multiple canonical pole aliases"); }
    pole = &candidate;
  }
  if (pole == nullptr) { throw std::invalid_argument(context + " analytic exchange PDG " + std::to_string(pdg) + " has no canonical pole alias"); }
  return *pole;
}

// Validate the fixed spin representatives of one trajectory
void ValidateAliasSpins(const ExchangeParam &mapping, const gra::MPDG &particle_data, const regge::Signature signature, const std::string &context) {
  if (mapping.pole_spinX2 <= 0 || mapping.pole_spinX2 % 2 != 0) { throw std::invalid_argument(context + ".pole_spin must be a positive integer"); }
  const int pole_spin      = mapping.pole_spinX2 / 2;
  const int pole_signature = pole_spin % 2 == 0 ? 1 : -1;
  if (pole_signature != regge::Tau(signature)) { throw std::invalid_argument(context + ".pole_spin is incompatible with trajectory signature"); }

  bool has_analytic = false;
  for (const int alias : mapping.pdg) {
    const int spin_x2 = particle_data.FindByPDG(alias).spinX2;
    if (spin_x2 == gra::aux::kNullSpinX2) {
      has_analytic = true;
      (void)ResolvePole(mapping, particle_data, alias, context);
      continue;
    }
    if (spin_x2 < 0 || spin_x2 % 2 != 0) { throw std::invalid_argument(context + " has a non-bosonic fixed-spin alias"); }
  }
  if (!has_analytic) { throw std::invalid_argument(context + " has no analytic trajectory alias"); }
}

// Parse and validate every central trajectory to SOFT exchange mapping
void ReadExchangeParams(Param &param, const nlohmann::json &regge, const gra::MPDG &particle_data, const SoftModelPtr &soft_model) {
  if (!soft_model) { throw std::invalid_argument("PARAM_REGGE requires a SOFT model"); }
  const auto &rows = regge.at("EXCHANGES");
  if (!rows.is_array() || rows.empty()) { throw std::invalid_argument("PARAM_REGGE.EXCHANGES must be a non-empty row array"); }

  param.soft_model = soft_model;
  param.exchanges.clear();
  param.PDG_TO_INDEX.clear();
  param.PDG_TO_C.clear();
  std::size_t           pomeron_count = 0;
  std::set<std::size_t> mapped_soft_exchanges;

  for (const auto &k : indices(rows)) {
    const auto       &row     = rows.at(k);
    const std::string context = "PARAM_REGGE.EXCHANGES[" + std::to_string(k) + "]";
    if (!row.is_object()) { throw std::invalid_argument(context + " must be an object"); }
    if (row.size() != 4 || !row.contains("role") || !row.contains("pdg") || !row.contains("pole_spin") || !row.contains("soft_exchange")) {
      throw std::invalid_argument(context + " must contain only role, pdg, pole_spin and soft_exchange");
    }

    ExchangeParam mapping;
    mapping.role = ParseExchangeParamRole(row.at("role"), context + ".role");
    mapping.pdg  = ParsePDGAliasEntry(row.at("pdg"), context + ".pdg");
    if (!row.at("pole_spin").is_number_integer()) { throw std::invalid_argument(context + ".pole_spin must be an integer"); }
    const long long pole_spin = row.at("pole_spin").get<long long>();
    if (pole_spin <= 0 || pole_spin > std::numeric_limits<int>::max() / 2) { throw std::invalid_argument(context + ".pole_spin must be a supported positive integer"); }
    mapping.pole_spinX2 = 2 * static_cast<int>(pole_spin);
    if (!row.at("soft_exchange").is_string()) { throw std::invalid_argument(context + ".soft_exchange must be a string"); }
    mapping.soft_exchange_name = row.at("soft_exchange").get<std::string>();
    mapping.soft_exchange      = soft_model->ExchangeId(mapping.soft_exchange_name);
    if (!mapped_soft_exchanges.insert(mapping.soft_exchange.Value()).second) {
      throw std::invalid_argument(context +
                                  ".soft_exchange duplicates mapped SOFT "
                                  "exchange " +
                                  mapping.soft_exchange_name);
    }

    const auto &soft = soft_model->Exchange(mapping.soft_exchange);
    if (!soft.Enabled()) { throw std::invalid_argument(context + " maps to a disabled SOFT exchange"); }
    if (soft.Role() != mapping.role) { throw std::invalid_argument(context + " role does not match mapped SOFT exchange role"); }

    const int crossing = RegisterExchangeAliases(param, mapping, k, particle_data, context);
    if (soft.CrossingParity() != crossing || regge::Tau(soft.Signature()) != crossing) { throw std::invalid_argument(context + " crossing parity or signature does not match its aliases"); }
    ValidateAliasSpins(mapping, particle_data, soft.Signature(), context);
    const std::size_t n = soft_model->GoodWalker().ChannelCount();
    if (soft.CouplingMatrix().size_row() != n || soft.CouplingMatrix().size_col() != n || soft.FormFactorParameters().size() != n) { throw std::invalid_argument(context + " has inconsistent Good Walker dimensions"); }

    if (mapping.role == SoftExchangeRole::Pomeron) {
      param.pomeron_trajectory = k;
      ++pomeron_count;
    }
    param.exchanges.push_back(std::move(mapping));
  }

  if (pomeron_count != 1) { throw std::invalid_argument("PARAM_REGGE.EXCHANGES must contain exactly one pomeron row"); }
}


// Check the values of the central Regge normalization parameters
void ValidateReggeParameters(const Param &param) {
  if (!std::isfinite(param.s0) || param.s0 <= 0.0) { throw std::invalid_argument("PARAM_REGGE.s0 must be finite and positive"); }
}

// Compute true for the finite-pole secondary-Reggeon aliases
//
bool IsFiniteSecondaryReggeonPDG(int pdg) {
  const int id = std::abs(pdg);
  return id == 9915 || id == 9925 || id == 9933 || id == 9943;
}

// Compute whether one trajectory contains an f2, a2, rho or omega pole
//
bool IsSecondaryReggeonTrajectory(const Param &param, int pdg) {
  if (pdg == PDG::PDG_gamma) { return false; }
  const auto &aliases = param.exchanges.at(TrajectoryIndex(param, pdg)).pdg;
  return std::any_of(aliases.begin(), aliases.end(), IsFiniteSecondaryReggeonPDG);
}

// Compute whether one continuum vertex is enabled in direct multi-Regge ladders
bool VertexAllowed(const Param &param, const VertexParam &vertex) {
  if (vertex.first == PDG::PDG_gamma || vertex.second == PDG::PDG_gamma) { return false; }
  return param.con.multiregge_secondary_exchanges || (!IsSecondaryReggeonTrajectory(param, vertex.first) && !IsSecondaryReggeonTrajectory(param, vertex.second));
}

// Compute the PDG representative used for secondary Reggeon selection checks
const gra::MParticle &SecondaryReggeonRepresentative(const Param &param, const gra::MPDG &pdg_table, int exchange_pdg) {
  const auto &aliases = param.exchanges.at(TrajectoryIndex(param, exchange_pdg)).pdg;

  const gra::MParticle &requested = pdg_table.FindByPDG(exchange_pdg);
  if (requested.spinX2 >= 0 && IsFiniteSecondaryReggeonPDG(exchange_pdg)) { return requested; }

  for (const int alias : aliases) {
    const gra::MParticle &particle = pdg_table.FindByPDG(alias);
    if (particle.spinX2 >= 0 && IsFiniteSecondaryReggeonPDG(alias) && particle.chargeX3 == requested.chargeX3 && particle.isospinX2 == requested.isospinX2 && particle.P == requested.P && particle.C == requested.C &&
        particle.G == requested.G) {
      return particle;
    }
  }
  throw std::invalid_argument(
      "PARAM_REGGE: mapped secondary Reggeon trajectory has no fixed pole "
      "with the requested PDG quantum numbers");
}

// Read rows [PDG, W0, B_gammaPV_cross_section, alpha0, alpha_prime]
//
void ReadPhotoParams(Param &param, const nlohmann::json &rows) {
  if (!rows.is_array() || rows.empty()) { throw std::invalid_argument("PARAM_REGGE.photoprod must be a non-empty row array"); }

  param.photo_channels.clear();
  for (const auto &i : indices(rows)) {
    const auto       &row     = rows[i];
    const std::string context = "PARAM_REGGE.photoprod[" + std::to_string(i) + "]";
    if (!row.is_array() || row.size() != 5 || !row[0].is_number_integer()) {
      throw std::invalid_argument(context +
                                  " must be [PDG, W0, B_gammaPV_cross_section, "
                                  "alpha0, alpha_prime]");
    }
    for (std::size_t k = 1; k < row.size(); ++k) {
      if (!row[k].is_number()) { throw std::invalid_argument(context + " contains a non-numeric parameter"); }
    }

    const int        pdg = row[0].get<int>();
    const PhotoParam channel{row[1].get<double>(), row[2].get<double>(), row[3].get<double>(), row[4].get<double>()};
    if (pdg <= 0 || !(std::isfinite(channel.W0) && channel.W0 > 0.0) || !(std::isfinite(channel.B_gammaPV) && channel.B_gammaPV >= 0.0) || !(std::isfinite(channel.a0) && channel.a0 > 0.0) ||
        !(std::isfinite(channel.ap) && channel.ap >= 0.0)) {
      throw std::invalid_argument(context + " contains an invalid physical parameter");
    }
    const auto &pomeron = param.soft_model->Exchange(param.exchanges.at(param.pomeron_trajectory).soft_exchange);
    regge::CheckEta(channel.a0, channel.ap, pomeron.Signature(), param.photoprod_eta_mode, context);
    if (!param.photo_channels.emplace(pdg, channel).second) { throw std::invalid_argument(context + " duplicates absolute PDG " + std::to_string(pdg)); }
  }
}

// Read optional [PDG, x_scale_squared, double_log_coefficient] corrections to the energy law
void ReadPhotoDLog(Param &param, const nlohmann::json &rows, bool diss = false) {
  if (!rows.is_array()) { throw std::invalid_argument("Photoproduction double-log parameters must be a row array"); }
  std::set<std::pair<int, double>> seen;
  for (const auto &row : rows) {
    if (!row.is_array() || row.size() != 3 || !row[0].is_number_integer() || !row[1].is_number() || !row[2].is_number()) {
      throw std::invalid_argument("Photoproduction double-log rows require [PDG, scale2, c]");
    }
    const int pdg = row[0].get<int>();
    const auto found = param.photo_channels.find(pdg);
    if (found == param.photo_channels.end() || (diss && !found->second.diss)) {
      throw std::invalid_argument("Photoproduction double-log rows require a configured channel");
    }
    auto &channel = found->second;
    const double scale2 = row[1].get<double>(), c = row[2].get<double>();
    const double w0 = diss ? channel.diss->W0 : channel.W0;
    if (!std::isfinite(scale2) || !(scale2 > 0.0 && scale2 < math::pow2(w0)) || !std::isfinite(c) || (!diss && c < 0.0)) {
      throw std::invalid_argument("Photoproduction double-log rows require 0 < scale2 < W0^2 and finite c, elastic c >= 0");
    }
    if (!seen.emplace(pdg, diss ? scale2 : 0.0).second) {
      throw std::invalid_argument("Duplicate photoproduction double-log row");
    }
    if (std::abs(c) > 0.0) {
      const PhotoDLog term{scale2, c, std::sqrt(std::log(math::pow2(w0) / scale2))};
      if (diss) { channel.diss_dlog.push_back(term); }
      else { channel.dlog = term; }
    }
  }
}

// Read HERA dissociation rows and normalize the continuum in the reference mass interval
void ReadPhotoDissParams(Param &param, const nlohmann::json &rows) {
  for (const auto& [pdg, diss] : flux::ReadPhotoDiss(rows)) {
    const auto found = param.photo_channels.find(pdg);
    if (found == param.photo_channels.end()) {
      throw std::invalid_argument("PARAM_REGGE.photoprod_diss requires a configured photoproduction PDG");
    }
    const auto& pomeron = param.soft_model->Exchange(param.exchanges.at(param.pomeron_trajectory).soft_exchange);
    regge::CheckEta(1.0 + diss.delta / 4.0, 0.0, pomeron.Signature(), param.photoprod_eta_mode, "photoprod_diss");
    if (found->second.dlog && !(math::pow2(diss.W0) > found->second.dlog->scale2)) {
      throw std::invalid_argument("Dissociation reference W0 must exceed the elastic double-log scale");
    }
    found->second.diss = diss;
  }
}

// Read pole-anchored exchanged-meson trajectories [PDG, spin, alpha_prime]
//
void ReadMesonTrajectories(Param &param, const nlohmann::json &rows, const gra::MPDG &pdg_table) {
  if (!rows.is_array() || rows.empty()) { throw std::invalid_argument("PARAM_REGGE.meson_trajectories must be a non-empty row array"); }

  param.meson_trajectories.clear();
  for (const auto &i : indices(rows)) {
    const auto       &row     = rows[i];
    const std::string context = "PARAM_REGGE.meson_trajectories[" + std::to_string(i) + "]";
    if (!row.is_array() || row.size() != 3 || !row[0].is_number_integer() || !row[1].is_number() || !row[2].is_number()) { throw std::invalid_argument(context + " must be [PDG, spin, alpha_prime]"); }

    const int       pdg = row[0].get<int>();
    const MesonTraj trajectory{row[1].get<double>(), row[2].get<double>()};
    const double    two_spin = 2.0 * trajectory.spin;
    if (pdg <= 0 || !std::isfinite(trajectory.spin) || trajectory.spin < 0.0 || std::abs(two_spin - std::round(two_spin)) > 1e-12 || !(std::isfinite(trajectory.ap) && trajectory.ap > 0.0)) {
      throw std::invalid_argument(context + " contains an invalid physical parameter");
    }
    if (std::abs(two_spin - static_cast<double>(pdg_table.FindByPDG(pdg).spinX2)) > 1e-12) {
      throw std::invalid_argument(context + " spin does not match the physical pole PDG " + std::to_string(pdg));
    }
    if (!param.meson_trajectories.emplace(pdg, trajectory).second) { throw std::invalid_argument(context + " duplicates absolute PDG " + std::to_string(pdg)); }
  }
}

// Compute the trajectory index for an exchange PDG alias
//
std::size_t TrajectoryIndex(const Param &param, int pdg) {
  const auto it = param.PDG_TO_INDEX.find(pdg);
  if (it == param.PDG_TO_INDEX.end()) { throw std::invalid_argument("PARAM_REGGE: exchange PDG " + std::to_string(pdg) + " is not listed in PARAM_REGGE.EXCHANGES"); }
  return it->second;
}

// Compute the mapped trajectory spin at one generated transfer
double Alpha(const regge::Param &param, const int exchange_pdg, const double transfer) {
  const std::size_t trajectory = regge::TrajectoryIndex(param, exchange_pdg);
  return param.soft_model->Alpha(param.exchanges.at(trajectory).soft_exchange, transfer);
}

// Compute the canonical fixed-spin pole matching one trajectory alias
const MParticle &PoleRepresentative(const Param &param, const MPDG &pdg_table, const int pdg) {
  const ExchangeParam &mapping = param.exchanges.at(TrajectoryIndex(param, pdg));
  return ResolvePole(mapping, pdg_table, pdg, "PARAM_REGGE.PoleRepresentative");
}

// Compute the photoproduction parameters for one central vector-meson channel
//
const PhotoParam &Photo(const Param &param, int pdg) {
  const int  key = std::abs(pdg);
  const auto it  = param.photo_channels.find(key);
  if (it == param.photo_channels.end()) { throw std::invalid_argument("PARAM_REGGE.photoprod has no channel for absolute PDG " + std::to_string(key)); }
  return it->second;
}

// Compute the pole-anchored trajectory for one exchanged meson
//
const MesonTraj &MesonTrajectory(const Param &param, int pdg) {
  const int  key = std::abs(pdg);
  const auto it  = param.meson_trajectories.find(key);
  if (it == param.meson_trajectories.end()) { throw std::invalid_argument("PARAM_REGGE.meson_trajectories has no entry for absolute PDG " + std::to_string(key)); }
  return it->second;
}

// Resolve the selected continuum amplitude permutation construction
//
int PermCount(const std::vector<gra::MDecayBranch> &decaytree, std::size_t n_central, PermType mode) {
  if (decaytree.size() < n_central) {
    throw std::invalid_argument(
        "PARAM_CON.MULTI.permutations decay tree has fewer "
        "central particles than requested");
  }
  if (mode == PermType::All) { return 1; }
  if (mode == PermType::Charged) { return 0; }

  for (std::size_t i = 0; i < n_central; ++i) {
    if (decaytree[i].p.chargeX3 == 0) { return 1; }
  }
  return 0;
}

// Compute the C-parity assigned to an exchange PDG alias
//
int CParity(const Param &param, int pdg) {
  if (pdg == PDG::PDG_gamma) { return -1; }
  const auto it = param.PDG_TO_C.find(pdg);
  if (it == param.PDG_TO_C.end()) { throw std::invalid_argument("PARAM_REGGE: exchange PDG " + std::to_string(pdg) + " is not listed in PARAM_REGGE.EXCHANGES"); }
  return it->second;
}

// Compute the product C-parity of a two-exchange state
//
int PairCParity(const Param &param, int first_pdg, int second_pdg) { return CParity(param, first_pdg) * CParity(param, second_pdg); }

// Check a secondary-Reggeon vertex using the complete two-body JPC LS subspace
//
VertexCheck CheckVertex(const Param &param, const gra::MPDG &pdg_table, int exchange_pdg, const std::vector<int> &pair_pdgs) {
  if (pair_pdgs.size() != 2) {
    throw std::invalid_argument(
        "PARAM_CON: Reggeon vertex quantum-number check requires two "
        "final-state PDGs");
  }

  VertexCheck result;
  if (!IsSecondaryReggeonTrajectory(param, exchange_pdg)) { return result; }

  result.applies                         = true;
  const gra::MParticle &representative   = SecondaryReggeonRepresentative(param, pdg_table, exchange_pdg);
  result.representative_pdg              = representative.pdg;
  const std::vector<gra::MParticle> legs = {
      pdg_table.FindByPDG(pair_pdgs[0]),
      pdg_table.FindByPDG(pair_pdgs[1]),
  };
  if (representative.chargeX3 != legs[0].chargeX3 + legs[1].chargeX3) {
    result.allowed = false;
    result.reason  = "electric charge is not conserved";
    return result;
  }

  // [REFERENCE: Ewerz, Maniatis, Nachtmann, Annals Phys. 342 (2014) 31, arXiv:1309.3478]
  // [REFERENCE: Lebiedowicz, Nachtmann, Szczurek, Phys. Rev. D94 (2016) 034017, arXiv:1606.05126]
  try {
    (void)gra::spin::BuildAutomaticCentralCoupling(representative, legs, gra::spin::AutoCentralCouplingMode::MinL, false, "PARAM_CON Reggeon quantum-number check", "", false, true, true, gra::spin::VertexContext::Auto);
  } catch (const std::invalid_argument &e) {
    const std::string message = e.what();
    if (message.find("no allowed LS coupling") == std::string::npos) { throw; }
    result.allowed = false;
    result.reason =
        "no two-body LS state conserves J, P, C and "
        "identical-particle statistics";
  }
  return result;
}

// Compute the crossing sign picked up when an exchange couples to an
// antiparticle beam
//
// The proton-exchange vertex is the reference convention. Replacing the beam
// by an antiproton charge-conjugates the vertex, so C-even exchanges keep the
// same sign and C-odd exchanges flip sign. The phase is common to all Good
// Walker transitions and is not a Jacob-Wick helicity phase
double AntiparticleSign(const Param &param, int exchange_pdg, int particle_pdg) { return (particle_pdg < 0) ? static_cast<double>(CParity(param, exchange_pdg)) : 1.0; }

using ContinuumPairKey = std::pair<bool, std::pair<int, int>>;

// Parse one signed compact PDG pair key or wildcard entry
//
std::pair<ContinuumPairKey, std::vector<int>> ParseContinuumPairKey(const std::string &key, const std::string &context) {
  if (key == "[*]") { return {{true, {0, 0}}, {}}; }
  nlohmann::json pair;
  try {
    pair = nlohmann::json::parse(key);
  } catch (const nlohmann::json::parse_error &) { throw std::invalid_argument(context + " has invalid PDG pair key " + key); }
  if (!pair.is_array() || pair.size() != 2 || !pair[0].is_number_integer() || !pair[1].is_number_integer()) { throw std::invalid_argument(context + " PDG pair key must contain two integers: " + key); }
  const int  first   = pair[0].get<int>();
  const int  second  = pair[1].get<int>();
  const auto ordered = first <= second ? std::make_pair(first, second) : std::make_pair(second, first);
  return {{false, ordered}, {first, second}};
}

// Parse one decimal exchange PDG key
int ParseContinuumExchangeKey(const std::string &key, const std::string &context) {
  std::size_t end = 0;
  int         pdg = 0;
  try {
    pdg = std::stoi(key, &end);
  } catch (const std::exception &) { throw std::invalid_argument(context + " has invalid exchange key " + key); }
  if (pdg == 0 || end != key.size()) { throw std::invalid_argument(context + " has invalid exchange key " + key); }
  return pdg;
}

using VertexKey = std::tuple<ReggeProductionModel, int, int, int>;
struct VertexInput {
  VertexForm form;
  ReggeizeParam reggeize;
  VetoParam veto;
};
using VertexInputs = std::map<VertexKey, VertexInput>;

// Validate continuum propagation and zero-secondary-production settings
void ValidateContinuumControls(const nlohmann::json &block, const std::string &label) {
  const auto &reggeize = block.at("reggeize");
  if (!reggeize.is_object() || reggeize.size() != 2 || !reggeize.contains("active") ||
      !reggeize.at("active").is_boolean() || !reggeize.contains("freeze_scale2") || !reggeize.at("freeze_scale2").is_number()) {
    throw std::invalid_argument(label + " reggeize requires boolean active and numeric freeze_scale2");
  }
  const double scale2 = reggeize.at("freeze_scale2").get<double>();
  if (!std::isfinite(scale2) || scale2 <= 0.0) {
    throw std::invalid_argument(label + " reggeize requires finite freeze_scale2 > 0 in GeV^2");
  }
  const auto &veto = block.at("pveto");
  if (!veto.is_object() || veto.size() != 3 || !veto.contains("active") || !veto.at("active").is_boolean() ||
      !veto.contains("M0") || !veto.at("M0").is_number() || !veto.contains("c") || !veto.at("c").is_number()) {
    throw std::invalid_argument(label + " pveto must contain boolean active and numeric M0,c");
  }
  const double mass = veto.at("M0").get<double>();
  const double power = veto.at("c").get<double>();
  if (!std::isfinite(mass) || mass <= 0.0 || !std::isfinite(power) || power < 0.0) {
    throw std::invalid_argument(label + " pveto requires finite M0 > 0 and c >= 0");
  }
}

// Read explicit continuum vertex forms and propagation settings from the model cards
VertexInputs ReadContinuumVertices(const Param &param, const MModelTune &tune) {
  VertexInputs vertices;
  for (const auto model : {ReggeProductionModel::MP, ReggeProductionModel::XP, ReggeProductionModel::GP}) {
    const std::string name = ReggeProductionModelName(model);
    const auto &card = tune.Continuum(name);
    if (!card.is_object() || card.empty()) { throw std::invalid_argument("CON_" + name + ".json must be nonempty"); }
    for (const auto &[exchange_key, exchange] : card.items()) {
      const std::string context = "CON_" + name + ".json." + exchange_key;
      const int exchange_pdg = ParseContinuumExchangeKey(exchange_key, context);
      if (exchange_pdg != PDG::PDG_gamma) { TrajectoryIndex(param, exchange_pdg); }
      if (!exchange.is_object()) { throw std::invalid_argument(context + " must be an object"); }
      for (const auto &[pair_key, block] : exchange.items()) {
        const std::string label = context + "." + pair_key;
        const auto [parsed, pair] = ParseContinuumPairKey(pair_key, label);
        (void)parsed;
        if (pair.size() != 2 || pair[0] <= 0 || pair[1] <= 0 || pair[0] > pair[1] || pair_key != "[" + std::to_string(pair[0]) + "," + std::to_string(pair[1]) + "]") {
          throw std::invalid_argument(label + " must identify one compact sorted positive-PDG pair");
        }
        ValidateContinuumControls(block, label);
        VertexInput vertex;
        vertex.form.transfer = ReadFF(block.at("FF_transfer"), label + ".FF_transfer");
        vertex.form.offshell = ReadFF(block.at("FF_offshell"), label + ".FF_offshell");
        if (vertex.form.transfer.type != FFType::None && vertex.form.transfer.norm != FFNorm::Zero) {
          throw std::invalid_argument(label + ".FF_transfer must use norm zero when active");
        }
        if (vertex.form.offshell.type != FFType::None && vertex.form.offshell.norm != FFNorm::Pole) {
          throw std::invalid_argument(label + ".FF_offshell must use norm pole when active");
        }
        vertex.reggeize = {block.at("reggeize").at("active").get<bool>(), block.at("reggeize").at("freeze_scale2").get<double>()};
        const auto &veto = block.at("pveto");
        vertex.veto = {veto.at("active").get<bool>(), veto.at("M0").get<double>(), veto.at("c").get<double>()};
        if (!vertices.emplace(std::make_tuple(model, exchange_pdg, pair[0], pair[1]), vertex).second) {
          throw std::invalid_argument(label + " duplicates a continuum vertex");
        }
      }
    }
  }
  return vertices;
}

// Parse and expand the unordered trajectory pairs of one continuum row
void ReadContinuumChannels(PairParam &entry, const nlohmann::json &pairs, const std::string &context, const Param &param, const gra::MPDG &pdg_table, ReggeProductionModel model) {
  if (!pairs.is_array() || pairs.empty()) { throw std::invalid_argument(context + " exchange pair list must be a non-empty array"); }

  std::map<std::pair<int, int>, std::size_t> seen;
  for (const auto &index : indices(pairs)) {
    const auto &pair = pairs[index];
    if (!pair.is_array() || pair.size() != 2 || !pair[0].is_number_integer() || !pair[1].is_number_integer()) { throw std::invalid_argument(context + " exchange pair " + std::to_string(index) + " must be [exchange_top, exchange_bottom]"); }
    const int first  = pair[0].get<int>();
    const int second = pair[1].get<int>();
    for (const int exchange : {first, second}) {
      if (exchange == PDG::PDG_gamma) { continue; }
      TrajectoryIndex(param, exchange);
      const int spin_x2 = pdg_table.FindByPDG(exchange).spinX2;
      if (model == ReggeProductionModel::GP && spin_x2 != gra::aux::kNullSpinX2) { throw std::invalid_argument(context + " GP exchange PDG " + std::to_string(exchange) + " must be a null-spin analytic alias"); }
      if ((model == ReggeProductionModel::MP || model == ReggeProductionModel::XP) && spin_x2 == gra::aux::kNullSpinX2) {
        throw std::invalid_argument(context + " fixed-spin model exchange PDG " + std::to_string(exchange) + " must have a PDG spin");
      }
    }

    if (model != ReggeProductionModel::GP) {
      for (const int exchange : {first, second}) {
        const VertexCheck check = CheckVertex(param, pdg_table, exchange, entry.pdg);
        if (check.applies && !check.allowed) {
          throw std::invalid_argument(context + " exchange PDG " + std::to_string(exchange) + " cannot couple to final-state PDGs " + FormatPDGVector(entry.pdg) + " using finite Reggeon pole PDG " +
                                      std::to_string(check.representative_pdg) + ": " + check.reason);
        }
      }
    }

    const auto unordered = first <= second ? std::make_pair(first, second) : std::make_pair(second, first);
    const auto inserted  = seen.emplace(unordered, index);
    if (!inserted.second) {
      throw std::invalid_argument(context + " duplicate unordered exchange pair [" + std::to_string(first) + ", " + std::to_string(second) + "] duplicates exchange pair " + std::to_string(inserted.first->second) +
                                  "; list only one of [A, B] or [B, A]");
    }
    entry.channels.push_back({first, second});
    if (first != second) { entry.channels.push_back({second, first}); }
  }
}

// Parse one continuum exchange table from the steering card
//
std::vector<PairParam> ReadContinuumTable(const nlohmann::json &table, const std::string &label, const Param &param, const gra::MPDG &pdg_table, ReggeProductionModel model) {
  if (!table.is_object()) { throw std::invalid_argument(label + " must be an object"); }
  if (table.empty()) { throw std::invalid_argument(label + " must contain at least one explicit entry"); }

  std::vector<PairParam>                  out;
  std::map<ContinuumPairKey, std::size_t> seen_active;
  std::size_t                             row = 0;
  for (const auto &[key, pairs] : table.items()) {
    const std::string context = label + "." + key;
    PairParam         entry;
    const auto [active_key, pdgs] = ParseContinuumPairKey(key, label);
    if (active_key.first) { throw std::invalid_argument(label + " does not support a default entry"); }
    const auto inserted = seen_active.emplace(active_key, row);
    if (!inserted.second) { throw std::invalid_argument(context + " duplicates the order-equivalent final-state pair at row " + std::to_string(inserted.first->second)); }

    entry.pdg = pdgs;
    ReadContinuumChannels(entry, pairs, context, param, pdg_table, model);
    out.push_back(entry);
    ++row;
  }
  return out;
}

// Compute the requested continuum steering table by production model
//
const std::vector<PairParam> &ContinuumTableByModel(const Param &param, const ReggeProductionModel model) {
  if (model == ReggeProductionModel::MP) { return param.con.MP; }
  if (model == ReggeProductionModel::XP) { return param.con.XP; }
  if (model == ReggeProductionModel::GP) { return param.con.GP; }
  throw std::invalid_argument("PARAM_CON: unsupported continuum production model");
}

// Match one explicit unordered stable two-body continuum entry
const PairParam *ExplicitPair(const std::vector<PairParam> &table, const std::vector<int> &pdg) {
  for (const auto &entry : table) {
    if (entry.pdg.size() != pdg.size()) { continue; }
    if (entry.pdg == pdg || std::is_permutation(entry.pdg.begin(), entry.pdg.end(), pdg.begin(), pdg.end())) { return &entry; }
  }
  return nullptr;
}

// Find only an explicit continuum entry for an unordered stable two-body
// final-state pair
//
const PairParam *FindPair(const Param &param, const std::vector<int> &pdg, ReggeProductionModel model) { return ExplicitPair(ContinuumTableByModel(param, model), pdg); }

// Compute true when explicit entries form one model-defined exchange ladder
bool ContinuumLadderEntriesConnect(const Param &param, const std::vector<const PairParam *> &entries, const ReggeProductionModel model, const bool pomeron_only) {
  if (entries.empty() || entries.size() > 3) { throw std::invalid_argument("PARAM_CON ladder requires one to three pair entries"); }
  const auto pomeron = [&](const int pdg) { return param.exchanges.at(TrajectoryIndex(param, pdg)).role == SoftExchangeRole::Pomeron; };
  const auto connected = [&](const int upper, const int lower) {
    if (model == ReggeProductionModel::GP) { return TrajectoryIndex(param, upper) == TrajectoryIndex(param, lower); }
    return upper == lower;
  };
  const auto search = [&](const auto &self, const std::size_t index, const int previous) -> bool {
    for (const auto &vertex : entries[index]->channels) {
      if (!VertexAllowed(param, vertex) || (index > 0 && !connected(previous, vertex.first))) { continue; }
      if (pomeron_only && ((index == 0 && !pomeron(vertex.first)) || (index + 1 == entries.size() && !pomeron(vertex.second)))) { continue; }
      if (index + 1 == entries.size() || self(self, index + 1, vertex.second)) { return true; }
    }
    return false;
  };
  return search(search, 0, 0);
}

// Resolve direct multi-body continuum ladder permutations supported by a model
std::vector<std::vector<int>> LadderPermutations(const gra::LORENTZSCALAR &lts, const Param &param, std::size_t n_central, ReggeProductionModel model) {
  (void)ContinuumTableByModel(param, model);
  const std::string model_name = ReggeProductionModelName(model);
  const std::string context    = "PARAM_CON." + model_name + " ladder";
  if (n_central != 4 && n_central != 6) { throw std::invalid_argument(context + " requires four or six central particles"); }

  const auto                    candidates = ContinuumLadderPermutations(lts, static_cast<int>(n_central), context, param.con.permutations);
  std::vector<std::vector<int>> supported;
  for (const auto &permutation : candidates) {
    std::vector<const PairParam *> entries;
    entries.reserve(n_central / 2);
    for (std::size_t pair = 0; pair < n_central; pair += 2) {
      const std::size_t first  = DecayIndex(permutation[pair], lts, context);
      const std::size_t second = DecayIndex(permutation[pair + 1], lts, context);
      const auto       *entry  = FindPair(param, {lts.decaytree[first].p.pdg, lts.decaytree[second].p.pdg}, model);
      if (entry == nullptr) {
        entries.clear();
        break;
      }
      entries.push_back(entry);
    }
    if (entries.empty()) { continue; }
    const auto &topologies = param.con.multiregge_topologies.at(static_cast<int>(n_central));
    const bool allowed = std::any_of(topologies.begin(), topologies.end(), [&](const Topology &topology) {
      std::size_t begin = 0;
      for (const int count : topology) {
        const std::size_t end = begin + static_cast<std::size_t>(count) / 2;
        if (!ContinuumLadderEntriesConnect(param, {entries.begin() + begin, entries.begin() + end}, model, topology.size() > 1)) { return false; }
        begin = end;
      }
      return true;
    });
    if (allowed) { supported.push_back(permutation); }
  }
  return supported;
}

// Validate cached direct multi-body ladder permutations
void CheckLadderPermutations(const std::vector<std::vector<int>> &permutations, const std::size_t n_central, const std::string &context) {
  constexpr int offset = 3;
  if (permutations.empty()) { throw std::invalid_argument(context + " has no ladder permutations"); }
  std::set<std::vector<int>> unique;
  for (const auto &permutation : permutations) {
    if (permutation.size() != n_central) { throw std::invalid_argument(context + " has a malformed ladder permutation"); }
    std::vector<bool> seen(n_central, false);
    for (const int amplitude_index : permutation) {
      if (amplitude_index < offset) { throw std::invalid_argument(context + " has a ladder index below the central system"); }
      const std::size_t index = static_cast<std::size_t>(amplitude_index - offset);
      if (index >= n_central || seen[index]) { throw std::invalid_argument(context + " has a non-permutation ladder row"); }
      seen[index] = true;
    }
    if (!unique.insert(permutation).second) { throw std::invalid_argument(context + " has a duplicate ladder permutation"); }
  }
}

// Validate cached serial and parallel multi-Regge topology partitions
void CheckMultiReggeTopologies(const std::vector<Topology> &topologies, const std::size_t n_central, const std::string &context) {
  if (n_central != 4 && n_central != 6) { throw std::invalid_argument(context + " requires four or six central particles"); }
  if (topologies.empty()) { throw std::invalid_argument(context + " has no topology partitions"); }

  std::set<Topology> unique;
  for (const auto &topology : topologies) {
    if (topology.empty()) { throw std::invalid_argument(context + " has an empty topology"); }
    std::size_t total = 0;
    for (const auto &block : indices(topology)) {
      const int size = topology[block];
      if (size != 2 && size != 4 && size != 6) { throw std::invalid_argument(context + " has an unsupported ladder size"); }
      if (block > 0 && topology[block - 1] < size) { throw std::invalid_argument(context + " has a non-descending topology partition"); }
      total += static_cast<std::size_t>(size);
    }
    if (total != n_central) { throw std::invalid_argument(context + " has a topology with the wrong sum"); }
    if (!unique.insert(topology).second) { throw std::invalid_argument(context + " has a duplicate topology"); }
  }
}

// Find the continuum entry for an unordered stable two-body final-state pair
//
const PairParam &Pair(const Param &param, const std::vector<int> &pdg, ReggeProductionModel model) {
  if (const PairParam *entry = FindPair(param, pdg, model); entry != nullptr) { return *entry; }
  throw std::invalid_argument(std::string("PARAM_CON.") + ReggeProductionModelName(model) + ": no explicit channel entry for final-state PDGs " + FormatPDGVector(pdg));
}

// Attach explicit vertex forms and shared meson-line settings before sampling
void AttachContinuumVertices(std::vector<PairParam> &table, const VertexInputs &vertices,
                             const Param &param, ReggeProductionModel model) {
  for (auto &entry : table) {
    const int left = std::abs(entry.pdg.at(0));
    const int right = std::abs(entry.pdg.at(1));
    const auto pair = std::minmax(left, right);
    const std::string label = "CON_" + std::string(ReggeProductionModelName(model)) + " pair " + FormatPDGVector(entry.pdg);
    for (const auto &[key, vertex] : vertices) {
      const auto &[source_model, exchange, first, second] = key;
      if (source_model != model || first != pair.first || second != pair.second) { continue; }
      if (entry.forms.empty()) {
        entry.reggeize = vertex.reggeize;
        entry.veto = vertex.veto;
      } else if (entry.reggeize.active != vertex.reggeize.active ||
                 std::abs(entry.reggeize.freeze_scale2 - vertex.reggeize.freeze_scale2) > std::numeric_limits<double>::epsilon() * std::max(entry.reggeize.freeze_scale2, vertex.reggeize.freeze_scale2) || entry.veto.active != vertex.veto.active ||
                 std::abs(entry.veto.M0 - vertex.veto.M0) > std::numeric_limits<double>::epsilon() * std::max(entry.veto.M0, vertex.veto.M0) ||
                 std::abs(entry.veto.c - vertex.veto.c) > std::numeric_limits<double>::epsilon() * std::max(entry.veto.c, vertex.veto.c)) {
        throw std::invalid_argument(label + " exchanges must share reggeize and pveto for their common meson line");
      }
      entry.forms.emplace(exchange, vertex.form);
    }
    for (auto &channel : entry.channels) {
      for (const std::size_t leg : {std::size_t{0}, std::size_t{1}}) {
        const int exchange = leg == 0 ? channel.first : channel.second;
        const auto form = entry.forms.find(exchange);
        if (form == entry.forms.end()) { throw std::invalid_argument(label + " is missing vertex PDG " + std::to_string(exchange)); }
        channel.transfer[leg] = form->second.transfer;
        channel.offshell[leg] = form->second.offshell;
      }
    }
    if (entry.reggeize.active) {
      for (const int pdg : entry.pdg) { MesonTrajectory(param, pdg); }
    }
  }
}

namespace {

// Parse one descending multi-Regge ladder partition
Topology ReadMultiReggeTopology(const nlohmann::json &row, const std::string &path, const int multiplicity) {
  if (!row.is_array() || row.empty()) { throw std::invalid_argument(path + " must be a nonempty array"); }
  Topology topology;
  topology.reserve(row.size());
  int sum = 0;
  for (const auto &index : indices(row)) {
    const auto &component = row.at(index);
    if (!component.is_number_integer()) { throw std::invalid_argument(path + " must contain integers"); }
    const int size = component.get<int>();
    if (size != 2 && size != 4 && size != 6) { throw std::invalid_argument(path + " contains an unsupported ladder size " + std::to_string(size)); }
    if (!topology.empty() && topology.back() < size) { throw std::invalid_argument(path + " must list ladder sizes in descending order"); }
    topology.push_back(size);
    sum += size;
  }
  if (sum != multiplicity) { throw std::invalid_argument(path + " must sum to " + std::to_string(multiplicity)); }
  return topology;
}

// Read one coherent bank of serial and parallel multi-Regge topologies
std::map<int, std::vector<Topology>> ReadMultiReggeTopologies(const nlohmann::json &block) {
  std::map<int, std::vector<Topology>> output;
  const auto                          &partitions = block.at("partitions");
  if (!partitions.is_object()) { throw std::invalid_argument("PARAM_CON.MULTI.partitions must be an object"); }
  for (const int multiplicity : {4, 6}) {
    const std::string key  = std::to_string(multiplicity);
    const std::string path = "PARAM_CON.MULTI.partitions." + key;
    if (!partitions.contains(key) || !partitions.at(key).is_array() || partitions.at(key).empty()) { throw std::invalid_argument(path + " must be a nonempty array"); }

    std::set<Topology> unique;
    auto              &topologies = output[multiplicity];
    for (const auto &row_index : indices(partitions.at(key))) {
      const auto       &row      = partitions.at(key).at(row_index);
      const std::string row_path = path + "[" + std::to_string(row_index) + "]";
      Topology          topology = ReadMultiReggeTopology(row, row_path, multiplicity);
      if (!unique.insert(topology).second) { throw std::invalid_argument(path + " contains a duplicate topology"); }
      topologies.push_back(std::move(topology));
    }
  }
  return output;
}

// Read and validate the dedicated multi-Regge ladder parameters
void ReadMultiReggeParameters(Param &param, const nlohmann::json &block) {
  if (!block.is_object()) { throw std::invalid_argument("PARAM_CON.MULTI must be an object"); }
  const std::set<std::string> fields = {"FF_transfer_ext", "FF_transfer_int", "secondary_exchanges", "partitions", "permutations"};
  for (const auto &[field, value] : block.items()) {
    (void)value;
    if (!fields.contains(field)) { throw std::invalid_argument("PARAM_CON.MULTI has unknown field " + field); }
  }
  for (const auto &field : fields) {
    if (!block.contains(field)) { throw std::invalid_argument("PARAM_CON.MULTI is missing field " + field); }
  }
  if (!block.at("secondary_exchanges").is_boolean()) { throw std::invalid_argument("PARAM_CON.MULTI.secondary_exchanges must be boolean"); }
  if (!block.at("FF_transfer_ext").is_boolean() || !block.at("FF_transfer_int").is_boolean()) { throw std::invalid_argument("PARAM_CON.MULTI.FF_transfer_ext and FF_transfer_int must be boolean"); }
  if (!block.at("permutations").is_string()) { throw std::invalid_argument("PARAM_CON.MULTI.permutations must be auto, charged or all"); }

  param.con.multiregge_transfer_ext        = block.at("FF_transfer_ext").get<bool>();
  param.con.multiregge_transfer_int        = block.at("FF_transfer_int").get<bool>();
  param.con.multiregge_secondary_exchanges = block.at("secondary_exchanges").get<bool>();
  param.con.multiregge_topologies          = ReadMultiReggeTopologies(block);
  const std::string permutations           = block.at("permutations").get<std::string>();
  if (permutations == "auto") {
    param.con.permutations = PermType::Auto;
  } else if (permutations == "charged") {
    param.con.permutations = PermType::Charged;
  } else if (permutations == "all") {
    param.con.permutations = PermType::All;
  } else {
    throw std::invalid_argument("PARAM_CON.MULTI.permutations must be auto, charged or all");
  }
}

}  // namespace


// Read generic Regge amplitude parameters into an immutable block
//
Param ReadParam(const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table, const MModelTune &tune) {
  try {
    Param                 param;
    const nlohmann::json &j          = tune.General();
    const SoftModelPtr   &soft_model = tune.Soft();

    const std::string XID   = "PARAM_REGGE";
    const auto       &regge = j.at(XID);

    ReadExchangeParams(param, regge, pdg_table, soft_model);
    param.photoprod_eta_mode = regge::ParseEta(regge.at("photoprod_eta_mode").get<std::string>(), "PARAM_REGGE.photoprod_eta_mode");
    if (param.photoprod_eta_mode != EtaMode::RotatingT0 && param.photoprod_eta_mode != EtaMode::Rotating) { throw std::invalid_argument("PARAM_REGGE.photoprod_eta_mode must be rotating_t0 or rotating"); }

    if (!regge.at("s0").is_number()) { throw std::invalid_argument("PARAM_REGGE.s0 must be numeric"); }
    param.s0       = regge.at("s0").get<double>();
    param.gp_alpha_min = tune.Numerics().at("NUMERICS_REGGE").at("GP_alpha_min").get<double>();
    if (!std::isfinite(param.gp_alpha_min) || param.gp_alpha_min <= -0.5 || param.gp_alpha_min > 1.0) {
      throw std::invalid_argument("NUMERICS_REGGE.GP_alpha_min must satisfy -0.5 < GP_alpha_min <= 1");
    }
    for (const std::string key : {"use_zeta", "omega"}) {
      const auto &table = regge.at(key);
      if (!table.is_object() || table.size() != 3) { throw std::invalid_argument("PARAM_REGGE." + key + " must contain MP, XP and GP"); }
      for (const auto model : {ReggeProductionModel::MP, ReggeProductionModel::XP, ReggeProductionModel::GP}) {
        const std::string name = ReggeProductionModelName(model);
        if (!table.contains(name)) { throw std::invalid_argument("PARAM_REGGE." + key + " is missing " + name); }
        if (key == "use_zeta") {
          param.use_zeta[model] = table.at(name).get<bool>();
        } else {
          if (!table.at(name).is_number() || !std::isfinite(table.at(name).get<double>())) { throw std::invalid_argument("PARAM_REGGE.omega." + name + " must be finite and numeric"); }
          param.omega[model] = table.at(name).get<double>();
        }
      }
    }
    ValidateReggeParameters(param);
    ReadPhotoParams(param, regge.at("photoprod"));
    ReadPhotoDLog(param, regge.value("photoprod_dlog", nlohmann::json::array()));
    ReadPhotoDissParams(param, regge.at("photoprod_diss"));
    ReadPhotoDLog(param, regge.value("photoprod_diss_dlog", nlohmann::json::array()), true);
    ReadMesonTrajectories(param, regge.at("meson_trajectories"), pdg_table);
    const auto                 &continuum   = regge.at("PARAM_CON");
    const std::set<std::string> identifiers = {"MP", "XP", "GP", "MULTI"};
    if (!continuum.is_object()) { throw std::invalid_argument("PARAM_REGGE.PARAM_CON must be an object"); }
    for (const auto &[model, value] : continuum.items()) {
      (void)value;
      if (!identifiers.contains(model)) { throw std::invalid_argument("PARAM_REGGE.PARAM_CON has unknown identifier " + model); }
    }
    for (const auto &identifier : identifiers) {
      if (!continuum.contains(identifier)) { throw std::invalid_argument("PARAM_REGGE.PARAM_CON is missing identifier " + identifier); }
    }
    const auto vertices = ReadContinuumVertices(param, tune);
    ReadMultiReggeParameters(param, continuum.at("MULTI"));
    param.con.MP = ReadContinuumTable(continuum.at("MP"), "PARAM_REGGE.PARAM_CON.MP", param, pdg_table, ReggeProductionModel::MP);
    param.con.XP = ReadContinuumTable(continuum.at("XP"), "PARAM_REGGE.PARAM_CON.XP", param, pdg_table, ReggeProductionModel::XP);
    param.con.GP = ReadContinuumTable(continuum.at("GP"), "PARAM_REGGE.PARAM_CON.GP", param, pdg_table, ReggeProductionModel::GP);
    AttachContinuumVertices(param.con.MP, vertices, param, ReggeProductionModel::MP);
    AttachContinuumVertices(param.con.XP, vertices, param, ReggeProductionModel::XP);
    AttachContinuumVertices(param.con.GP, vertices, param, ReggeProductionModel::GP);

    if (final_pdgs.size() == 2) {
      param.con.form_pdgs = {final_pdgs[0], final_pdgs[1]};
    } else if (final_pdgs.size() == 1) {
      const int pdg       = final_pdgs[0];
      param.con.form_pdgs = {pdg, -pdg};
      if (pdg == 111 || pdg == 113 || pdg == 333) { param.con.form_pdgs = {pdg, pdg}; }
    } else {
      // Multi-body continuum vertices carry their own propagation settings
      param.con.form_pdgs = {0, 0};
    }

    std::cout << "MRegge::ReadParameters: [" << XID << "]" << std::endl;
    std::cout << regge << std::endl << std::endl;
    std::cout << "MRegge::ReadParameters: [PARAM_REGGE.PARAM_CON.MULTI]" << std::endl;
    std::cout << continuum.at("MULTI") << std::endl << std::endl;
    return param;
  } catch (const std::exception &e) { throw std::invalid_argument("MRegge::ReadParameters: Parse error in " + tune.GeneralFile() + ": " + e.what()); }
}

// Construct one shared immutable Regge parameter block
ParamPtr ReadParamPtr(const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table, const MModelTune &tune) { return std::make_shared<const Param>(ReadParam(final_pdgs, pdg_table, tune)); }

// Compute one run owned immutable Regge parameter block
ParamPtr GetParam(MModelCache &cache, const std::vector<int> &final_pdgs, const gra::MPDG &pdg_table) {
  const std::string key = "regge:" + FormatPDGVector(final_pdgs);
  return cache.Get<Param>(key, [&cache, &final_pdgs, &pdg_table] { return ReadParamPtr(final_pdgs, pdg_table, cache.Tune()); });
}

}  // namespace regge

}  // namespace gra
