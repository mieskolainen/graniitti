// Resonance card parsing and spin-state configuration
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <compare>
#include <complex>
#include <iostream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::zi;

namespace gra {
namespace resonance {

namespace {

// Read one discrete label before narrowing it to the supported integer range
int ReadInteger(const nlohmann::json &value, int minimum, int maximum, const std::string &context) {
  if (!value.is_number_integer()) { throw std::invalid_argument(context + " must be an integer"); }
  if (value.is_number_unsigned()) {
    const auto number = value.get<nlohmann::json::number_unsigned_t>();
    if (number > static_cast<nlohmann::json::number_unsigned_t>(maximum)) {
      throw std::invalid_argument(context + " is outside the supported range");
    }
  } else {
    const auto number = value.get<nlohmann::json::number_integer_t>();
    if (number < minimum || number > maximum) {
      throw std::invalid_argument(context + " is outside the supported range");
    }
  }
  const int number = value.get<int>();
  if (number < minimum) { throw std::invalid_argument(context + " is outside the supported range"); }
  return number;
}

// Read a finite nonnegative resonance mass or width
double ReadMassWidth(const nlohmann::json &value, const std::string &context) {
  if (!value.is_number()) { throw std::invalid_argument(context + " must be numeric"); }
  const double number = value.get<double>();
  if (!std::isfinite(number) || number < 0.0) {
    throw std::invalid_argument(context + " must be finite and nonnegative");
  }
  return number;
}

// Parse one integer parity without narrowing oversized JSON values
int ParseParity(const nlohmann::json &value, const std::string &context) {
  if (!value.is_number_integer()) { throw std::invalid_argument(context + " parity must be -1 or 1"); }
  if (value.is_number_unsigned()) {
    if (value.get<nlohmann::json::number_unsigned_t>() == 1U) { return 1; }
  } else {
    const auto parity = value.get<nlohmann::json::number_integer_t>();
    if (parity == -1 || parity == 1) { return static_cast<int>(parity); }
  }
  throw std::invalid_argument(context + " parity must be -1 or 1");
}

// Parse one strict resonance denominator mode
BreitWigner ParseBreitWigner(const nlohmann::json &value, const std::string &label) {
  if (!value.is_string()) { throw std::invalid_argument("gra::resonance::Read: <" + label + "> BW must be a string"); }
  const std::string mode = value.get<std::string>();
  if (mode == "fixed-width") { return BreitWigner::FixedWidth; }
  if (mode == "kinematic-width") { return BreitWigner::KinematicWidth; }
  if (mode == "running-width") { return BreitWigner::RunningWidth; }
  throw std::invalid_argument("gra::resonance::Read: <" + label +
                              "> BW must be 'fixed-width', 'kinematic-width' or 'running-width'");
}

// Compute the canonical inline production block key for one ordered channel
//
std::string ProductionChannelKey(const int pdg1, const int pdg2) {
  return "[" + std::to_string(pdg1) + "," + std::to_string(pdg2) + "]";
}

// Parse one canonical unordered production-channel object key
//
std::array<int, 2> ParseProductionChannelKey(const std::string &key, const std::string &label) {
  if (key.size() < 5 || key.front() != '[' || key.back() != ']') {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> invalid channel key " + key);
  }
  const std::size_t comma = key.find(',');
  if (comma == std::string::npos || key.find(',', comma + 1) != std::string::npos) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> invalid channel key " + key);
  }
  int first  = 0;
  int second = 0;
  try {
    std::size_t used_first  = 0;
    std::size_t used_second = 0;
    first                   = std::stoi(key.substr(1, comma - 1), &used_first);
    second                  = std::stoi(key.substr(comma + 1, key.size() - comma - 2), &used_second);
    if (used_first != comma - 1 || used_second != key.size() - comma - 2) {
      throw std::invalid_argument("trailing characters");
    }
  } catch (const std::exception &) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> invalid channel key " + key);
  }
  if (first > second || key != ProductionChannelKey(first, second)) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> channel key " + key +
                                " must use canonical ascending PDGs without whitespace");
  }
  return {first, second};
}

// Read and cache one canonical model level resonance phase
double ReadModelPhase(const nlohmann::json &model, const std::string &label) {
  if (!model.contains("phi") || !model.at("phi").is_number()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> requires numeric phi");
  }
  const double phi = model.at("phi").get<double>();
  if (!std::isfinite(phi) || phi < -PI || phi >= PI) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> phi must be finite in [-pi, pi)");
  }
  return phi;
}

// Read explicit model level resonance form factors
void ReadModelFormFactors(const nlohmann::json &model, RES_MODEL_FORM &out, const std::string &label) {
  out.ff_transfer = regge::ReadFF(model.at("FF_transfer"), "gra::resonance::Read: <" + label + "> FF_transfer");
  if (out.ff_transfer.type != regge::FFType::None && out.ff_transfer.norm != regge::FFNorm::Zero) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> FF_transfer must use norm zero");
  }
  out.ff_prod = regge::ReadFF(model.at("FF_prod"), "gra::resonance::Read: <" + label + "> FF_prod");
  if (out.ff_prod.type != regge::FFType::None &&
      ((out.ff_prod.type != regge::FFType::Gaussian && out.ff_prod.type != regge::FFType::Vector) ||
       out.ff_prod.norm != regge::FFNorm::Pole)) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> FF_prod must be gaussian or vector with norm pole");
  }
}

// Validate one MP polarization block while retaining both steering choices
void ValidateMPPolarization(const nlohmann::json &polarization, const std::string &label) {
  if (!polarization.is_object() || !polarization.contains("mode") || !polarization.at("mode").is_string()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> requires string mode");
  }
  const std::string mode = polarization.at("mode");
  if (mode != "none" && mode != "a_Jz" && mode != "rho") {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> mode should be 'none', 'a_Jz' or 'rho'");
  }
  const std::set<std::string> fields = {"mode", "a_Jz", "rho_mag", "rho_phase", "random_rho"};
  for (const auto &[field, value] : polarization.items()) {
    (void)value;
    if (!fields.contains(field)) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "> has unknown field " + field);
    }
  }

  std::set<std::string> required;
  if (mode == "a_Jz") {
    required = {"a_Jz"};
  } else if (mode == "rho") {
    required = {"rho_mag", "rho_phase", "random_rho"};
  }
  for (const auto &field : required) {
    if (!polarization.contains(field)) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "> is missing field " + field);
    }
  }

  for (const auto &field : {"a_Jz", "rho_mag", "rho_phase"}) {
    if (polarization.contains(field) && !polarization.at(field).is_array()) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "> field " + field + " should be an array");
    }
  }
  if (polarization.contains("random_rho") && !polarization.at("random_rho").is_boolean()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> field random_rho should be boolean");
  }
}

// Validate width-derived coupling notation
bool WidthDerived(const nlohmann::json &rows, ReggeVertexBasis basis, int pdg1, int pdg2, const std::string &label,
                  bool gp) {
  std::size_t null_count = 0;
  if (IsAutomaticReggeVertexBasis(basis)) {
    null_count = rows.at(0).is_null() ? 1 : 0;
  } else {
    for (const auto &row : rows) { null_count += row.at(2).is_null() ? 1 : 0; }
  }
  const bool diphoton = pdg1 == PDG::PDG_gamma && pdg2 == PDG::PDG_gamma;
  if ((diphoton && null_count != 1) || (!diphoton && null_count != 0)) {
    throw std::invalid_argument("gra::resonance::Read: <" + label +
                                "> requires exactly one null magnitude only for gamma-gamma");
  }
  const bool shaped = basis == ReggeVertexBasis::LS || (gp && basis == ReggeVertexBasis::Helicity);
  if (diphoton && !shaped && rows.size() != (IsAutomaticReggeVertexBasis(basis) ? 2 : 1)) {
    throw std::invalid_argument("gra::resonance::Read: <" + label +
                                "> multiple gamma-gamma width couplings require an LS shape or a GP "
                                "raw helicity shape");
  }
  return diphoton;
}

// Read one nonnegative magnitude and canonical phase
std::complex<double> ReadCoupling(const nlohmann::json &mag_value, const nlohmann::json &phase_value,
                                  const std::string &label) {
  if ((!mag_value.is_number() && !mag_value.is_null()) || !phase_value.is_number()) {
    throw std::invalid_argument(label + " coupling must be numeric");
  }
  const double mag   = mag_value.is_null() ? 1.0 : mag_value.get<double>();
  const double phase = phase_value.get<double>();
  if (!std::isfinite(mag) || mag < 0.0 || !std::isfinite(phase) || phase < -PI || phase >= PI ||
      (std::fpclassify(mag) == FP_ZERO && std::fpclassify(phase) != FP_ZERO)) {
    throw std::invalid_argument(label + " coupling magnitude or phase is invalid");
  }
  return std::polar(mag, phase);
}

// Read one sparse LS coupling table
spin::LSCoefficients ReadLS(const nlohmann::json &rows, const std::string &label) {
  if (!rows.is_array() || rows.empty()) { throw std::invalid_argument(label + " requires non-empty g_ls"); }
  spin::LSCoefficients out;
  for (const auto &i : indices(rows)) {
    const auto &row = rows[i];
    if (!row.is_array() || row.size() != 4 || !row[0].is_number_integer() || !row[1].is_number()) {
      throw std::invalid_argument(label + " malformed g_ls row");
    }
    const int       l     = ReadInteger(row[0], 0, std::numeric_limits<int>::max() / 2, label + " L");
    const double    s     = row[1];
    const double    two_s = 2.0 * s;
    if (l < 0 || l > std::numeric_limits<int>::max() / 2 || !std::isfinite(s) || s < 0.0 ||
        two_s > std::numeric_limits<int>::max() || std::abs(two_s - std::round(two_s)) > 1.0e-9) {
      throw std::invalid_argument(label + " invalid LS quantum numbers");
    }
    if (!out.Insert(static_cast<std::size_t>(l), static_cast<std::size_t>(std::llround(two_s)),
                    ReadCoupling(row[2], row[3], label))) {
      throw std::invalid_argument(label + " duplicate g_ls row");
    }
  }
  return out;
}

// Read one sparse direct helicity coupling table
void ReadHelicity(const nlohmann::json &rows, const std::string &label, RES_PRODUCTION_CHANNEL &out) {
  if (!rows.is_array() || rows.empty()) { throw std::invalid_argument(label + " requires non-empty helicity"); }
  std::set<std::array<long long, 2>> seen;
  for (const auto &row : rows) {
    if (!row.is_array() || row.size() != 4 || !row[0].is_number() || !row[1].is_number()) {
      throw std::invalid_argument(label + " malformed helicity row");
    }
    const std::array<double, 2>    index  = {row[0].get<double>(), row[1].get<double>()};
    if (!std::isfinite(index[0]) || !std::isfinite(index[1]) ||
        std::abs(index[0]) > std::numeric_limits<int>::max() / 2.0 ||
        std::abs(index[1]) > std::numeric_limits<int>::max() / 2.0) {
      throw std::invalid_argument(label + " helicity is outside the supported range");
    }
    const std::array<long long, 2> index2 = {std::llround(2.0 * index[0]), std::llround(2.0 * index[1])};
    if (std::abs(2.0 * index[0] - index2[0]) > 1.0e-9 ||
        std::abs(2.0 * index[1] - index2[1]) > 1.0e-9 || !seen.insert(index2).second) {
      throw std::invalid_argument(label + " invalid or duplicate helicity row");
    }
    out.helicity.push_back(index);
    out.g_helicity.push_back(ReadCoupling(row[2], row[3], label));
  }
}

// Read and validate one complete model-specific resonance channel
RES_PRODUCTION_CHANNEL ReadProductionChannel(const nlohmann::json &block, std::array<int, 2> exchange,
                                             ReggeProductionModel model, const std::string &label) {
  const std::set<std::string> allowed =
      model == ReggeProductionModel::MP
          ? std::set<std::string>{"basis", "Lambda", "CP", "g", "g_ls", "helicity", "polarization"}
          : std::set<std::string>{"basis", "Lambda", "CP", "g", "g_ls", "helicity"};
  if (!block.is_object()) { throw std::invalid_argument(label + " must be an object"); }
  for (const auto &[field, value] : block.items()) {
    (void)value;
    if (!allowed.contains(field)) { throw std::invalid_argument(label + " has unknown field " + field); }
  }
  if (!block.contains("basis") || !block["basis"].is_string() || !block.contains("CP") || !block["CP"].is_array() ||
      block["CP"].size() != 2 || !block["CP"][0].is_boolean() || !block["CP"][1].is_boolean()) {
    throw std::invalid_argument(label + " requires basis and boolean CP pair");
  }
  RES_PRODUCTION_CHANNEL out;
  out.exchange   = exchange;
  out.basis      = ParseReggeVertexBasis(block["basis"].get<std::string>(), ReggeVertexRole::Resonance);
  out.C_symmetry = block["CP"][0];
  out.P_symmetry = block["CP"][1];
  if (model == ReggeProductionModel::GP && (!out.C_symmetry || !out.P_symmetry)) {
    throw std::invalid_argument(label + " GP requires CP = [true,true]");
  }
  const bool automatic = IsAutomaticReggeVertexBasis(out.basis);
  if (model == ReggeProductionModel::GP && automatic) {
    throw std::invalid_argument(label + " GP does not support automatic basis");
  }
  const char *field = automatic ? "g" : out.basis == ReggeVertexBasis::LS ? "g_ls" : "helicity";
  for (const char *candidate : {"g", "g_ls", "helicity"}) {
    if (block.contains(candidate) != (std::string(candidate) == field)) {
      throw std::invalid_argument(label + " fields do not match basis");
    }
  }
  if (model == ReggeProductionModel::MP) {
    if (!block.contains("polarization")) { throw std::invalid_argument(label + " requires polarization"); }
    ValidateMPPolarization(block["polarization"], label + ".polarization");
  }
  if (model != ReggeProductionModel::GP || out.basis == ReggeVertexBasis::LS) {
    if (!block.contains("Lambda") || !block["Lambda"].is_number()) {
      throw std::invalid_argument(label + " requires numeric Lambda");
    }
    out.Lambda = block["Lambda"];
    if (!std::isfinite(out.Lambda) || out.Lambda <= 0.0) {
      throw std::invalid_argument(label + " Lambda must be finite and positive");
    }
  } else if (block.contains("Lambda")) {
    throw std::invalid_argument(label + " selected basis does not use Lambda");
  }
  const auto &rows  = block[field];
  out.width_derived =
      WidthDerived(rows, out.basis, exchange[0], exchange[1], label, model == ReggeProductionModel::GP);
  if (automatic) {
    if (!rows.is_array() || rows.size() != 2) { throw std::invalid_argument(label + " g must be [magnitude, phase]"); }
    out.g = ReadCoupling(rows[0], rows[1], label);
  } else if (out.basis == ReggeVertexBasis::LS) {
    out.g_ls = ReadLS(rows, label);
  } else {
    ReadHelicity(rows, label, out);
    const bool fractional = std::any_of(out.helicity.cbegin(), out.helicity.cend(), [](const auto &h) {
      return std::abs(h[0] - std::round(h[0])) > 1.0e-9 || std::abs(h[1] - std::round(h[1])) > 1.0e-9;
    });
    if (model == ReggeProductionModel::GP && fractional) {
      throw std::invalid_argument(label + " GP helicities must be integers");
    }
  }
  return out;
}

// Read one complete MP, XP or GP resonance model
RES_PRODUCTION_MODEL ReadProductionModel(const nlohmann::json &model_block, ReggeProductionModel model,
                                         const std::string &label) {
  if (!model_block.is_object()) { throw std::invalid_argument(label + " must be an object"); }
  RES_PRODUCTION_MODEL out;
  out.phi = ReadModelPhase(model_block, label);
  ReadModelFormFactors(model_block, out, label);
  for (const auto &[key, block] : model_block.items()) {
    if (key == "mass" || key == "width" || key == "BW" || key == "phi" || key == "FF_transfer" || key == "FF_prod") { continue; }
    const auto pair = ParseProductionChannelKey(key, label);
    out.channels.push_back(ReadProductionChannel(block, pair, model, label + "." + key));
  }
  if (out.channels.empty()) { throw std::invalid_argument(label + " requires at least one channel"); }
  std::sort(out.channels.begin(), out.channels.end(),
            [](const auto &a, const auto &b) { return a.exchange < b.exchange; });
  return out;
}

// Require the positive-reflectivity state a(+m) = (-1)^m a(-m) in production-plane axes
void RequireAJzParityPair(const std::vector<std::complex<double>>& a_Jz, const int abs_jz, const int J,
                          const std::string& resparam_str) {
  const int                  neg_idx  = -abs_jz + J;
  const int                  pos_idx  = abs_jz + J;
  const std::complex<double> expected = (abs_jz % 2 == 0 ? 1.0 : -1.0) * a_Jz[neg_idx];
  const double               scale    = std::max({1.0, std::abs(a_Jz[pos_idx]), std::abs(expected)});
  if (std::abs(a_Jz[pos_idx] - expected) > 1e-9 * scale) {
    throw std::invalid_argument("gra::resonance::Read: <" + resparam_str +
                                "> explicit a_Jz +/-Jz pair violates CP[1]=true "
                                "state reflection a(+m)=(-1)^m a(-m)");
  }
}

// Insert one coherent Jz amplitude after validating its index and magnitude
//
void InsertAJzAmplitude(std::vector<std::complex<double>> &a_Jz, std::vector<bool> &seen, const int spinX2,
                        const int Jz, const double mag, const double phase, const std::string &resparam_str) {
  if (Jz < -spinX2 / 2 || Jz > spinX2 / 2) {
    throw std::invalid_argument("gra::resonance::Read: <" + resparam_str + "> a_Jz contains Jz outside [-J, J]");
  }
  const int idx = Jz + spinX2 / 2;
  if (seen[idx]) {
    throw std::invalid_argument("gra::resonance::Read: <" + resparam_str + "> a_Jz contains duplicate Jz entries");
  }
  if (!std::isfinite(mag) || !std::isfinite(phase) || mag < 0.0) {
    throw std::invalid_argument("gra::resonance::Read: <" + resparam_str +
                                "> a_Jz must have finite mag/phase and mag >= 0");
  }
  seen[idx] = true;
  a_Jz[idx] = std::polar(mag, phase);
}

// Read coherent spin amplitudes and complete positive-reflectivity partners
std::vector<std::complex<double>> ReadCoherentAJz(const nlohmann::json& res_model, const int spinX2,
                                                  const bool P_symmetry, const std::string& resparam_str) {
  const int                         n = spinX2 + 1;  // n = 2J + 1
  std::vector<std::complex<double>> a_Jz(n, 0.0);

  if (!res_model.contains("a_Jz")) {
    throw std::invalid_argument("gra::resonance::Read: <" + resparam_str + "> Missing MP inline a_Jz");
  }

  std::vector<bool> seen(n, false);
  for (const auto &entry : res_model.at("a_Jz")) {
    if (!entry.is_array() || entry.size() != 3) {
      throw std::invalid_argument("gra::resonance::Read: <" + resparam_str +
                                  "> each a_Jz entry must have [Jz, magnitude, phase]");
    }
    if (!entry[0].is_number_integer() || !entry[1].is_number() || !entry[2].is_number()) {
      throw std::invalid_argument("gra::resonance::Read: <" + resparam_str + "> has a non-numeric a_Jz entry");
    }

    const int    Jz    = ReadInteger(entry.at(0), -spinX2 / 2, spinX2 / 2, resparam_str + " Jz");
    const double mag   = entry.at(1);
    const double phase = entry.at(2);
    InsertAJzAmplitude(a_Jz, seen, spinX2, Jz, mag, phase, resparam_str);
  }

  const int J = spinX2 / 2;
  if (P_symmetry) {
    for (int m = 1; m <= J; ++m) {
      const int neg_idx = -m + J;
      const int pos_idx = m + J;
      if (seen[neg_idx] && !seen[pos_idx]) {
        seen[pos_idx] = true;
        a_Jz[pos_idx] = (m % 2 == 0 ? 1.0 : -1.0) * a_Jz[neg_idx];
      } else if (!seen[neg_idx] && seen[pos_idx]) {
        seen[neg_idx] = true;
        a_Jz[neg_idx] = (m % 2 == 0 ? 1.0 : -1.0) * a_Jz[pos_idx];
      } else if (seen[neg_idx] && seen[pos_idx]) {
        RequireAJzParityPair(a_Jz, m, J, resparam_str);
      }
    }

    // Omitted projections remain zero, including both members of an absent pair
  }

  return a_Jz;
}

// Read the active resonance polarization mode from one MP block
//
std::string ReadMPPolarizationMode(const nlohmann::json &polarization, const std::string &label) {
  if (!polarization.contains("mode") || !polarization.at("mode").is_string()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> is missing string mode");
  }
  const std::string basis = polarization.at("mode").get<std::string>();
  if (basis != "none" && basis != "a_Jz" && basis != "rho") {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> mode should be 'none', 'a_Jz' or 'rho'");
  }
  return basis;
}

// Require the spin-density reflection in production-plane axes
// [REFERENCE: Mieskolainen, https://arxiv.org/abs/1910.06300, Eq. (103)]
void RequireRhoParity(const MMatrix<std::complex<double>> &rho, const std::string &label) {
  const std::size_t n = rho.size_row();
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      if (std::abs(rho[i][j] - ((i + j) % 2 == 0 ? 1.0 : -1.0) * rho[n - 1 - i][n - 1 - j]) > 1.0e-10) {
        throw std::invalid_argument("gra::resonance::Read: <" + label + "> rho violates CP[1]=true parity reflection");
      }
    }
  }
}

// Require an exact spin dimension before reading a density matrix
void RequireRhoShape(const nlohmann::json &rows, std::size_t n, const std::string &label) {
  if (!rows.is_array() || rows.size() != n ||
      !std::all_of(rows.begin(), rows.end(), [n](const auto &row) { return row.is_array() && row.size() == n; })) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> rho matrices must have dimension 2J+1");
  }
}

// Read the optional resonance specific tensor Pomeron transfer form factor
void ReadTensorFormFactor(const nlohmann::json &model, RES_TENSOR_CHANNEL &out, const std::string &label) {
  if (!model.contains("FF_transfer") && !model.contains("FF_prod")) { return; }
  if (model.contains("FF_transfer")) {
    out.ff_transfer = regge::ReadFF(model.at("FF_transfer"), "gra::resonance::Read: <" + label + "> TP.FF_transfer");
    if (out.ff_transfer.type != regge::FFType::None && out.ff_transfer.norm != regge::FFNorm::Zero) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "> TP.FF_transfer must use norm zero");
    }
  }
  if (model.contains("FF_prod")) {
    out.ff_prod          = regge::ReadFF(model.at("FF_prod"), "gra::resonance::Read: <" + label + "> TP.FF_prod");
    if (out.ff_prod.type != regge::FFType::None && ((out.ff_prod.type != regge::FFType::Gaussian && out.ff_prod.type != regge::FFType::Vector) || out.ff_prod.norm != regge::FFNorm::Pole)) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "> TP.FF_prod must be gaussian or vector with norm pole");
    }
  }
}

// Read pair-keyed covariant tensor-exchange resonance couplings
RES_TENSOR_MODEL ReadTensorModel(const nlohmann::json &model, const PARAM_RES &res, const std::string &label) {
  if (!model.is_object()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> requires at least one channel object");
  }

  RES_TENSOR_MODEL out;
  out.phi = ReadModelPhase(model, label);
  ReadModelFormFactors(model, out, label);
  const int spinX2 = res.p.spinX2;
  for (const auto &[key, block] : model.items()) {
    if (key == "mass" || key == "width" || key == "phi" || key == "FF_transfer" || key == "FF_prod") { continue; }
    const auto pair = ParseProductionChannelKey(key, label);
    if (!block.is_object()) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key + "> must be an object");
    }
    const std::set<std::string> fields = {"g_tensor", "VMD_MIXING", "FF_transfer", "FF_prod"};
    for (const auto &[field, value] : block.items()) {
      (void)value;
      if (!fields.contains(field)) {
        throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key + "> has unknown field " + field);
      }
    }

    RES_TENSOR_CHANNEL channel;
    channel.exchange           = pair;
    channel.g_tensor           = block.at("g_tensor").get<std::vector<double>>();
    const std::size_t expected = spinX2 == 6 ? 1U : (spinX2 == 4 ? 7U : 2U);
    if ((spinX2 == 0 || spinX2 == 2 || spinX2 == 4 || spinX2 == 6) && channel.g_tensor.size() != expected) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key +
                                  "> g_tensor has invalid cardinality for resonance spin");
    }
    if (!std::all_of(channel.g_tensor.cbegin(), channel.g_tensor.cend(),
                     [](const double value) { return std::isfinite(value); })) {
      throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key + "> g_tensor must be finite");
    }

    ReadTensorFormFactor(block, channel, label + "." + key);
    if (!block.contains("FF_transfer")) {
      channel.ff_transfer          = out.ff_transfer;
    }
    if (!block.contains("FF_prod")) {
      channel.ff_prod          = out.ff_prod;
    }

    if (block.contains("VMD_MIXING")) {
      if (spinX2 != 2 || res.p.P != -1 || res.p.C != -1 || (pair[0] != 22 && pair[1] != 22)) {
        throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key +
                                    "> VMD_MIXING requires a photon channel and J^PC = 1--");
      }
      if (!block.at("VMD_MIXING").is_array()) {
        throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key + "> VMD_MIXING must be an array");
      }
      std::set<int> mixing_pdg;
      for (const auto &entry : block.at("VMD_MIXING")) {
        if (!entry.is_object() || entry.size() != 2 || !entry.contains("PDG") || !entry.contains("g_tensor")) {
          throw std::invalid_argument("gra::resonance::Read: <" + label + "." + key +
                                      "> VMD_MIXING entries require only PDG and g_tensor");
        }
        const int  pdg      = ReadInteger(entry.at("PDG"), -std::numeric_limits<int>::max(),
                                         std::numeric_limits<int>::max(), label + " VMD_MIXING PDG");
        const auto coupling = entry.at("g_tensor").get<std::vector<double>>();
        if (pdg == 0 || std::abs(pdg) == std::abs(res.p.pdg) || coupling.size() != 2 || !std::isfinite(coupling[0]) ||
            !std::isfinite(coupling[1])) {
          throw std::invalid_argument("gra::resonance::Read: <" + label + "> invalid VMD_MIXING entry");
        }
        if (!mixing_pdg.insert(std::abs(pdg)).second) {
          throw std::invalid_argument("gra::resonance::Read: <" + label + "> duplicate VMD_MIXING PDG");
        }
        channel.VMD_MIXING.push_back({std::abs(pdg), {coupling[0], coupling[1]}});
      }
    }
    out.channels.push_back(std::move(channel));
  }
  if (out.channels.empty()) {
    throw std::invalid_argument("gra::resonance::Read: <" + label + "> requires at least one channel object");
  }
  std::sort(
      out.channels.begin(), out.channels.end(),
      [](const RES_TENSOR_CHANNEL &left, const RES_TENSOR_CHANNEL &right) { return left.exchange < right.exchange; });
  bool found_reference = false;
  for (const auto &channel : out.channels) {
    for (const double coupling : channel.g_tensor) {
      if (std::fpclassify(coupling) == FP_ZERO) { continue; }
      if (!found_reference && coupling < 0.0) {
        throw std::invalid_argument("gra::resonance::Read: <" + label + "> first nonzero g_tensor must be positive");
      }
      found_reference = true;
    }
  }
  return out;
}

}  // namespace

// Read resonance parameters
// Read a resonance card into particle data and production model settings
PARAM_RES Read(const std::string &resonance_file, MRandom &rng, ReggeProductionModel model,
               const std::string &modelparam) {
  // =====================================================================

  std::cout << "gra::resonance::Read: Reading " + resonance_file + " " << std::endl;

  // Create a JSON object from file
  std::string inputfile = gra::ResolveModelDataFile(modelparam, resonance_file);

  // Read and parse
  std::string    data;
  nlohmann::json j;

  try {
    data = gra::aux::GetInputData(inputfile);
    j    = nlohmann::json::parse(data);
  } catch (const std::exception &error) {
    throw std::invalid_argument("gra::resonance::Read: Error reading " + inputfile + ": " + error.what());
  }

  // Resonance parameters
  gra::PARAM_RES res;
  res.modelparam = gra::ResolveModelTuneDir(modelparam);
  const std::string active = ReggeProductionModelName(model == ReggeProductionModel::None ? ReggeProductionModel::GP : model);

  try {
    const auto                 &param_res         = j.at("PARAM_RES");
    const auto                 &models            = param_res.at("MODELS");
    const std::set<std::string> model_identifiers = {"MP", "XP", "GP", "TP"};
    if (!models.is_object()) {
      throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> MODELS must be an object");
    }
    for (const auto &[name, value] : models.items()) {
      (void)value;
      if (!model_identifiers.contains(name)) {
        throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> MODELS has unknown identifier " +
                                    name);
      }
    }
    for (const auto &name : model_identifiers) {
      if (!models.contains(name)) {
        throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> MODELS is missing identifier " +
                                    name);
      }
      const auto &block = models.at(name);
      const std::string label = resonance_file + "::MODELS." + name;
      const double mass = ReadMassWidth(block.at("mass"), label + " mass");
      const double width = ReadMassWidth(block.at("width"), label + " width");
      const BreitWigner bw = name == "TP" ? BreitWigner::FixedWidth : ParseBreitWigner(block.at("BW"), label);
      if (name == active) {
        res.p.mass = mass;
        res.p.width = width;
        res.BW = bw;
      }
    }
    res.p.tau = res.p.width > 0.0 ? PDG::hbar / res.p.width : 0.0;



    // Particle name
    if (!param_res.at("name").is_string()) {
      throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> name must be a string !");
    }
    res.p.name = param_res.at("name").get<std::string>();
    if (res.p.name.empty()) {
      throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> name must be non-empty !");
    }

    // PDG code
    if (!param_res.at("PDG").is_number_integer()) {
      throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> PDG must be an integer !");
    }
    res.p.pdg = ReadInteger(param_res.at("PDG"), -std::numeric_limits<int>::max(),
                            std::numeric_limits<int>::max(), resonance_file + " PDG");
    if (res.p.pdg == 0) { throw std::invalid_argument(resonance_file + " PDG must be nonzero"); }


    // Spin
    if (!param_res.at("spinX2").is_number_integer()) {
      throw std::invalid_argument("gra::resonance::Read:: <" + resonance_file + "> spinX2 must be an integer !");
    }
    res.p.spinX2 = ReadInteger(param_res.at("spinX2"), 0, std::numeric_limits<int>::max() - 1,
                               resonance_file + " spinX2");
    if (res.p.spinX2 % 2 != 0) {
      throw std::invalid_argument("gra::resonance::Read: <" + resonance_file + "> requires integer resonance spin");
    }

    // Parity
    res.p.P = ParseParity(param_res.at("P"), resonance_file);

    // C-parity
    res.p.C = ReadInteger(param_res.at("C"), -1, 1, resonance_file + " C");
    res.C_from_card = true;

    res.XP = ReadProductionModel(models.at("XP"), ReggeProductionModel::XP, resonance_file + "::MODELS.XP");
    res.GP = ReadProductionModel(models.at("GP"), ReggeProductionModel::GP, resonance_file + "::MODELS.GP");
    {
      const RES_PRODUCTION_MODEL mp =
          ReadProductionModel(models.at("MP"), ReggeProductionModel::MP, resonance_file + "::MODELS.MP");
      res.MP.phi                  = mp.phi;
      res.MP.ff_transfer          = mp.ff_transfer;
      res.MP.ff_prod              = mp.ff_prod;
      res.MP.channels             = mp.channels;
    }
    // ------------------------------------------------------------------
    // Pair-keyed covariant tensor Pomeron and f2-Reggeon couplings
    res.TP = ReadTensorModel(models.at("TP"), res, resonance_file + "::MODELS.TP");
    // ------------------------------------------------------------------

    const int                     n = res.p.spinX2 + 1;  // n = 2J + 1
    MMatrix<std::complex<double>> rho(n, n, 0.0);
    if (res.MP.channels.empty()) {
      throw std::invalid_argument("gra::resonance::Read: <" + resonance_file + "> Missing MP inline channel blocks");
    }
    const auto &mp_model = models.at("MP");
    const auto &first    = res.MP.channels.front();
    const auto &polarization =
        mp_model.at(ProductionChannelKey(first.exchange[0], first.exchange[1])).at("polarization");
    for (std::size_t k = 1; k < res.MP.channels.size(); ++k) {
      const auto &channel = res.MP.channels[k];
      const auto &block   = mp_model.at(ProductionChannelKey(channel.exchange[0], channel.exchange[1]));
      if (block.at("polarization") != polarization || channel.C_symmetry != first.C_symmetry ||
          channel.P_symmetry != first.P_symmetry) {
        throw std::invalid_argument("gra::resonance::Read: <" + resonance_file +
                                    "> MP channels must share spin steering");
      }
    }
    res.spin_basis        = ReadMPPolarizationMode(polarization, resonance_file + "::MODELS.MP.polarization");
    const bool P_symmetry = first.P_symmetry;
    if (res.UsesCoherentSpinBasis()) {
      std::vector<std::complex<double>> a_Jz  = ReadCoherentAJz(polarization, res.p.spinX2, P_symmetry, resonance_file);
      const double norm2 = gra::SquaredNorm(a_Jz);
      if (!std::isfinite(norm2) || norm2 <= 0.0) {
        throw std::invalid_argument("gra::resonance::Read: <" + resonance_file + "> a_Jz must have non-zero norm");
      }
      if (std::abs(norm2 - 1.0) > 1e-6) {
        gra::aux::PrintWarning();
        std::cout << rang::fg::yellow << "gra::resonance::Read: <" << resonance_file
                  << "> a_Jz is not normalized (sum |a_Jz|^2 = " << norm2 << "), normalizing it automatically"
                  << rang::fg::reset << std::endl;
      }
      gra::Scale(a_Jz, 1.0 / std::sqrt(norm2));
      // Validate physical amplitudes after removing the arbitrary input scale
      if (P_symmetry) {
        const int J = res.p.spinX2 / 2;
        for (int m = 1; m <= J; ++m) { RequireAJzParityPair(a_Jz, m, J, resonance_file); }
      }
      res.a_Jz            = a_Jz;
      rho                 = gra::RankOneProjector(a_Jz);
      res.MP.random_rho   = false;
      res.rho             = rho;
    } else if (res.UsesDensitySpinBasis()) {
      if (!polarization.at("random_rho").is_boolean()) {
        throw std::invalid_argument("gra::resonance::Read: <" + resonance_file + "> random_rho should be Boolean");
      }
      res.MP.random_rho = polarization.at("random_rho");
      if (res.MP.random_rho) {
        rho = gra::spin::RandomRho(res.p.spinX2, P_symmetry, rng);
        if (P_symmetry) {
          std::vector<std::complex<double>> phase(n, 1.0);
          for (int m = 1; m <= res.p.spinX2 / 2; ++m) { phase[res.p.spinX2 / 2 + m] = (m % 2 == 0 ? 1.0 : -1.0); }
          rho = rho.LeftDiagonalProduct(phase).RightDiagonalProduct(phase);
        }
      } else {
        const auto &rho_mag   = polarization.at("rho_mag");
        const auto &rho_phase = polarization.at("rho_phase");
        RequireRhoShape(rho_mag, rho.size_row(), resonance_file);
        RequireRhoShape(rho_phase, rho.size_row(), resonance_file);
        for (std::size_t a = 0; a < rho.size_row(); ++a) {
          for (std::size_t b = 0; b < rho.size_col(); ++b) {
            if (!rho_mag.at(a).at(b).is_number() || !rho_phase.at(a).at(b).is_number()) {
              throw std::invalid_argument("gra::resonance::Read: <" + resonance_file + "> rho entries must be numeric");
            }
            const double mag   = rho_mag.at(a).at(b);
            const double phase = rho_phase.at(a).at(b);
            if (!std::isfinite(mag) || !std::isfinite(phase) || mag < 0.0) {
              throw std::invalid_argument("gra::resonance::Read: <" + resonance_file +
                                          "> rho entries should be finite and "
                                          "rho_mag should be non-negative");
            }
            rho[a][b] = std::polar(mag, phase);
          }
        }
      }
      (void)gra::spin::Positivity(rho, res.p.spinX2 / 2.0);
      if (P_symmetry) { RequireRhoParity(rho, resonance_file); }
      res.rho              = rho;
      res.a_Jz.clear();
    } else {
      res.MP.random_rho = false;
      res.a_Jz.clear();
      res.rho = MMatrix<std::complex<double>>();
    }
    std::cout << param_res << std::endl;
    std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;

  } catch (nlohmann::json::exception &e) {
    throw std::invalid_argument("resonance::Read: Missing parameter in '" + resonance_file + "' : " + e.what());
  }

  return res;
}

// Build coherent spin amplitudes from diagonal Jz probabilities and card phases
// a_Jz = sqrt(p_Jz) exp(i phi_Jz) with positive-reflectivity partners
std::vector<std::complex<double>> CoherentAJzFromDiagonalWeights(const gra::PARAM_RES      &res,
                                                                 const std::vector<double> &probabilities) {
  if (res.p.spinX2 < 0 || res.p.spinX2 % 2 != 0 ||
      probabilities.size() != static_cast<std::size_t>(res.p.spinX2) + 1 ||
      probabilities.size() != res.a_Jz.size()) {
    throw std::invalid_argument(
        "gra::resonance::CoherentAJzFromDiagonalWeights: spin dimension "
        "mismatch");
  }
  if (res.MP.channels.empty()) {
    throw std::invalid_argument(
        "gra::resonance::CoherentAJzFromDiagonalWeights: "
        "missing MP channel metadata");
  }

  std::vector<std::complex<double>> a_Jz(probabilities.size(), 0.0);
  for (const auto &k : indices(probabilities)) {
    if (!std::isfinite(probabilities[k]) || probabilities[k] < 0.0) {
      throw std::invalid_argument(
          "gra::resonance::CoherentAJzFromDiagonalWeights: invalid Jz "
          "probability");
    }
    const std::complex<double> phase = std::abs(res.a_Jz[k]) > 1e-12 ? res.a_Jz[k] / std::abs(res.a_Jz[k]) : 1.0;
    a_Jz[k]                          = msqrt(probabilities[k]) * phase;
  }

  const auto &channel = res.MP.channels.front();
  if (!channel.P_symmetry) { return a_Jz; }

  const int J     = res.p.spinX2 / 2;
  for (int m = 1; m <= J; ++m) {
    const int    neg_idx = -m + J;
    const int    pos_idx = m + J;
    const double scale   = std::max({1.0, probabilities[neg_idx], probabilities[pos_idx]});
    if (std::abs(probabilities[neg_idx] - probabilities[pos_idx]) > 1e-12 * scale) {
      throw std::invalid_argument(
          "gra::resonance::CoherentAJzFromDiagonalWeights: "
          "parity partners need equal weights");
    }

    std::complex<double> neg_phase = 1.0;
    if (std::abs(res.a_Jz[neg_idx]) > 1e-12) {
      neg_phase = res.a_Jz[neg_idx] / std::abs(res.a_Jz[neg_idx]);
    } else if (std::abs(res.a_Jz[pos_idx]) > 1e-12) {
      neg_phase = (m % 2 == 0 ? 1.0 : -1.0) * res.a_Jz[pos_idx] / std::abs(res.a_Jz[pos_idx]);
    }
    a_Jz[neg_idx] = msqrt(probabilities[neg_idx]) * neg_phase;
    a_Jz[pos_idx] = (m % 2 == 0 ? 1.0 : -1.0) * a_Jz[neg_idx];
  }
  return a_Jz;
}

// Compute one branching ratio from the active decay table
double DecayBranchingRatio(const int resonance_pdg, const std::vector<int> &products, const std::string &context, const std::string &modelparam) {
  const std::string inputfile = gra::ResolveModelDataFile(modelparam, "DECAYS.json");
  nlohmann::json    table;
  try {
    table = nlohmann::json::parse(gra::aux::GetInputData(inputfile));
  } catch (const std::exception &error) {
    throw std::invalid_argument(context + ": Error reading " + inputfile + ": " + error.what());
  }

  const std::string pdg_key = std::to_string(resonance_pdg);
  if (!table.contains(pdg_key) || !table.at(pdg_key).is_object()) {
    throw std::invalid_argument(context + ": DECAYS.json has no resonance PDG " + pdg_key);
  }

  std::string channel_key = "[";
  for (const auto &i : indices(products)) {
    if (i > 0) { channel_key += ","; }
    channel_key += std::to_string(products[i]);
  }
  channel_key += "]";
  // A two-body partial width is independent of the daughter ordering
  if (!table.at(pdg_key).contains(channel_key) && products.size() == 2) {
    const std::string reversed = "[" + std::to_string(products[1]) + "," + std::to_string(products[0]) + "]";
    if (table.at(pdg_key).contains(reversed)) { channel_key = reversed; }
  }
  if (!table.at(pdg_key).contains(channel_key)) {
    throw std::invalid_argument(context + ": DECAYS.json has no " + channel_key + " entry for resonance PDG " +
                                pdg_key);
  }
  const auto &channel = table.at(pdg_key).at(channel_key);
  if (!channel.is_object() || !channel.contains("BR") || !channel.at("BR").is_number()) {
    throw std::invalid_argument(context + ": " + channel_key + " has no numeric BR for resonance PDG " + pdg_key);
  }
  const double branching_ratio = channel.at("BR").get<double>();
  if (!std::isfinite(branching_ratio) || branching_ratio <= 0.0 || branching_ratio > 1.0) {
    throw std::invalid_argument(context +
                                ": BR should be finite and in (0,1] for "
                                "resonance PDG " +
                                pdg_key);
  }
  return branching_ratio;
}

// Compute the diphoton branching ratio from the active decay table
double GammaGammaBranchingRatio(const int resonance_pdg, const std::string &modelparam) {
  return DecayBranchingRatio(resonance_pdg, {PDG::PDG_gamma, PDG::PDG_gamma}, "GammaGammaBranchingRatio", modelparam);
}

// Compute the electronic partial width from particle and decay tables
// Gamma_ee = Gamma_tot BR(X -> e+e-)
double ElectronicPartialWidth(const MParticle &resonance, const std::string &modelparam) {
  if (!std::isfinite(resonance.width) || resonance.width <= 0.0) {
    throw std::invalid_argument(
        "ElectronicPartialWidth: resonance width should be finite and "
        "positive for PDG " +
        std::to_string(resonance.pdg));
  }
  return resonance.width * DecayBranchingRatio(resonance.pdg, {11, -11}, "ElectronicPartialWidth", modelparam);
}

// Compute the on-shell diphoton partial width from particle and decay tables
// Gamma_gg = Gamma_tot BR(X -> gamma gamma)
double GammaGammaPartialWidth(const MParticle &resonance, const std::string &modelparam) {
  if (!std::isfinite(resonance.width) || resonance.width <= 0.0) {
    throw std::invalid_argument(
        "GammaGammaPartialWidth: resonance width "
        "should be finite and positive for PDG " +
        std::to_string(resonance.pdg));
  }
  return resonance.width * GammaGammaBranchingRatio(resonance.pdg, modelparam);
}

// Compute the width-normalized reduced diphoton resonance coupling
// g_Xgg = sqrt[32 pi M_X (2J+1) Gamma(X -> gamma gamma)]
double GammaGammaResonanceCoupling(const MParticle &resonance, const std::string &modelparam) {
  if (!std::isfinite(resonance.mass) || resonance.mass <= 0.0) {
    throw std::invalid_argument(
        "GammaGammaResonanceCoupling: resonance mass "
        "should be finite and positive for PDG " +
        std::to_string(resonance.pdg));
  }
  if (resonance.spinX2 < 0 || resonance.spinX2 % 2 != 0) {
    throw std::invalid_argument(
        "GammaGammaResonanceCoupling: gamma-gamma fusion requires integer "
        "resonance spin for PDG " +
        std::to_string(resonance.pdg));
  }
  if (resonance.spinX2 == 2) {
    throw std::invalid_argument(
        "GammaGammaResonanceCoupling: Landau-Yang forbids two real photons "
        "coupling to spin one for PDG " +
        std::to_string(resonance.pdg));
  }
  if (resonance.C != 1) {
    throw std::invalid_argument(
        "GammaGammaResonanceCoupling: two photons "
        "require resonance C = +1 for PDG " +
        std::to_string(resonance.pdg));
  }

  // [REFERENCE: Particle Data Group, Review of Particle Physics, Kinematics]
  // The reduced two-photon matrix has unit helicity norm and includes 1/2!
  // phase space
  const double spin_states = static_cast<double>(resonance.spinX2 + 1);
  return msqrt(32.0 * PI * resonance.mass * spin_states * GammaGammaPartialWidth(resonance, modelparam));
}

// Reduced Breit-Wigner line shapes and form factors
//
// Useful identity for normalization:
//
// \int dm^2 \frac{1}{(m^2 - M0^2)^2 + M0^2Gamma^2} \equiv \frac{\pi}{M0 Gamma}
//
// The full Feynman propagator is i times the reduced line shape used in the
// published invariant amplitude M. Vertex and propagator i factors belong to
// the complete iM graph and must not be attached to one selected topology
//
std::complex<double> LineShape(const double mass2, const PARAM_RES &resonance, const double running_profile) {
  switch (resonance.BW) {
    case BreitWigner::FixedWidth:
      return FixedWidthLineShape(mass2, resonance.p.mass, resonance.p.width);
    case BreitWigner::KinematicWidth:
      return KinematicWidthLineShape(mass2, resonance.p.mass, resonance.p.width);
    case BreitWigner::RunningWidth:
      return RunningWidthLineShape(mass2, resonance.p.mass, resonance.p.width, running_profile);
  }
  throw std::invalid_argument("resonance::LineShape: Unknown Breit Wigner mode");
}

// See e.g
// [REFERENCE: TASI Lectures on propagators, users.ictp.it/~smr2244/tait-supplemental.pdf]
// [REFERENCE: Cacciapaglia, Deandrea, Curtis, arxiv.org/abs/0906.3417v2]
// [REFERENCE: www.t2.ucsd.edu/twiki2/pub/UCSDTier2/Physics214Spring2015/ajw-breit-wigner-cbx99-55.pdf]

// ----------------------------------------------------------------------
// Delta function \delta(\hat{s} - M_0^2) replacement function:
//    \int d\hat{s} \delta(\hat{s} - M_0^2) -> \int d\hat{s}
//    deltaBW(\hat{s},M0,Gamma)
//
// Compute the Breit-Wigner density applied at cross-section level
// delta_BW(s) = M Gamma/[pi((s-M^2)^2+M^2 Gamma^2)]
double BreitWignerDensity(double shat, double M0, double Gamma) {
  return M0 * Gamma / PI / (pow2(shat - pow2(M0)) + pow2(M0 * Gamma));
}

// Compute the positive Breit-Wigner factor applied at amplitude level
// A_BW(s) = sqrt(delta_BW(s))
double BreitWignerAmplitude(double shat, double M0, double Gamma) {
  return std::sqrt(BreitWignerDensity(shat, M0, Gamma));
}

// Compute the complex amplitude whose norm squared is BreitWignerDensity
// A_BW(s) = sqrt(M Gamma/pi)/(s-M^2+iM Gamma)
std::complex<double> ComplexBreitWignerAmplitude(double shat, double M0, double Gamma) {
  return std::sqrt(M0 * Gamma / PI) * FixedWidthLineShape(shat, M0, Gamma);
}
// ----------------------------------------------------------------------

// ----------------------------------------------------------------------
// Breit-Wigner propagator parametrizations

// Compute the reduced fixed-width relativistic Breit-Wigner line shape
// BW(s) = 1/(s-M0^2+i M0 Gamma)
std::complex<double> FixedWidthLineShape(double m2, double M0, double Gamma) {
  return 1.0 / (m2 - M0 * M0 + zi * M0 * Gamma);
}

// Compute the reduced kinematic-width relativistic Breit-Wigner line shape
// BW(s) = 1/(s-M0^2+i sqrt(s) Gamma)
std::complex<double> KinematicWidthLineShape(double m2, double M0, double Gamma) {
  if (!(m2 >= 0.0) || !std::isfinite(m2)) {
    throw AmplitudeFailure(
        "KinematicWidthLineShape: mass squared must be "
        "finite and non-negative");
  }
  return 1.0 / (m2 - M0 * M0 + zi * std::sqrt(m2) * Gamma);
}

// Compute the reduced running-width relativistic Breit-Wigner line shape
// BW(s) = 1/(s-M0^2+i M0 Gamma R(s))
std::complex<double> RunningWidthLineShape(double m2, double M0, double Gamma, double running_profile) {
  if (!std::isfinite(running_profile) || running_profile < 0.0) {
    throw AmplitudeFailure("RunningWidthLineShape: profile must be finite and non-negative");
  }
  return 1.0 / (m2 - M0 * M0 + zi * M0 * Gamma * running_profile);
}

// Compute the MadWidth spin-dependent propagator estimate
// A_J(s) = N_J(s)/(s-M0^2+i M0 Gamma)
// This decay-width estimate does not replace a physical spin propagator
//
// Here m2 denotes E^2 in the MadWidth spin numerator
//
// [REFERENCE: Alwall et al., arXiv:1402.1178v2, Table 1]
//
std::complex<double> JacksonLineShape(double m2, double M0, double Gamma, double J) {
  const std::complex<double> denom = (m2 - M0 * M0 + zi * M0 * Gamma);
  if (!std::isfinite(J) || J < 0.0 || J > 2.0 || std::abs(2.0 * J - std::round(2.0 * J)) > 1.0e-12) {
    throw std::invalid_argument("JacksonLineShape: J must be an integer or half-integer in [0, 2]");
  }
  const int two_J = static_cast<int>(std::llround(2.0 * J));

  if (two_J == 0) {  // J = 0
    return 1.0 / denom;
  } else if (two_J == 1) {  // J = 1/2
    return std::sqrt(m2) / denom;
  } else if (two_J == 2) {  // J = 1
    return (1.0 - m2 / (M0 * M0)) / denom;
  } else if (two_J == 3) {  // J = 3/2
    return (2.0 / 3.0) * std::sqrt(m2) * (1.0 - m2 / (M0 * M0)) / denom;
  } else if (two_J == 4) {  // J = 2
    return (7.0 / 6.0 - (4.0 / 3.0) * (m2 / (M0 * M0)) + (2.0 / 3.0) * (m2 * m2) / (gra::math::pow4(M0))) / denom;
  }
  throw std::logic_error("JacksonLineShape: invalid spin dispatch");
}

}  // namespace resonance
}  // namespace gra
