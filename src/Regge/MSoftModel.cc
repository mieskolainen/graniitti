// Soft Regge model parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Regge/MSoftModel.h"

// C++
#include <algorithm>
#include <cmath>
#include <iomanip>
#include <limits>
#include <set>
#include <sstream>
#include <stdexcept>
#include <utility>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Regge/MFormFactor.h"
#include "Graniitti/Regge/MReggeSig.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

// Libraries
#include "json.hpp"

using gra::aux::indices;

namespace gra {
namespace {

using json = nlohmann::json;

// Exchange quantum numbers owned by PARAM_SOFT.EXCHANGE_DEF
struct ExchangeDefinition {
  std::string name;
  SoftExchangeRole role = SoftExchangeRole::Reggeon;
  SoftTrajectoryMode trajectory_mode = SoftTrajectoryMode::Linear;
  regge::Signature signature = regge::Signature::Positive;
  int crossing_parity = 1;
};

// Named form factor family and its eigenstate parameter rows
struct FormFactorBank {
  SoftFormFactor type = SoftFormFactor::DualPower;
  std::vector<std::vector<double>> parameters;
};

// Optional helicity controls shared by all exchange rows
struct HelicityDefinition {
  bool enabled = false;
  double mass_scale = PDG::mp;
};

// Parsed exchange values and their stable name lookup
struct ParsedExchanges {
  std::vector<SoftExchange> values;
  std::map<std::string, SoftExchangeId> ids;
};

// Triple Pomeron ratio and its Regge eta prescription
struct TriplePomeronDefinition {
  double ratio = 0.0;
  EtaMode eta_mode = EtaMode::Raw;
};

// Parse one integer sign without narrowing oversized JSON values
int ParseUnitSign(const json &value, const std::string &path) {
  if (!value.is_number_integer()) {
    throw std::invalid_argument(path + " must be -1 or 1");
  }
  if (value.is_number_unsigned()) {
    if (value.get<json::number_unsigned_t>() == 1U) {
      return 1;
    }
  } else {
    const auto sign = value.get<json::number_integer_t>();
    if (sign == -1 || sign == 1) {
      return static_cast<int>(sign);
    }
  }
  throw std::invalid_argument(path + " must be -1 or 1");
}

// Parse one exact soft exchange role
SoftExchangeRole ParseExchangeRole(const std::string &value,
                                   const std::string &path) {
  if (value == "pomeron") {
    return SoftExchangeRole::Pomeron;
  }
  if (value == "reggeon") {
    return SoftExchangeRole::Reggeon;
  }
  if (value == "odderon") {
    return SoftExchangeRole::Odderon;
  }
  throw std::invalid_argument(path + " has unknown role " + value);
}

// Parse one exact soft trajectory mode
SoftTrajectoryMode ParseTrajectoryMode(const std::string &value,
                                       const std::string &path) {
  if (value == "linear") {
    return SoftTrajectoryMode::Linear;
  }
  if (value == "pion_loop") {
    return SoftTrajectoryMode::PionLoop;
  }
  throw std::invalid_argument(path + " has unknown trajectory_mode " + value);
}

// Read one positive finite forward excitation profile scale
double ReadForwardExcitationScale(const json &profile,
                                  const std::string &name) {
  const std::string path = "PARAM_SOFT.FORWARD_EXCITATION." + name;
  if (!profile.contains(name) || !profile.at(name).is_number()) {
    throw std::invalid_argument(path + " must be a number");
  }
  const double value = profile.at(name).get<double>();
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::invalid_argument(path + " must be finite and positive");
  }
  return value;
}

// Parse one exact eigenstate transition form factor mode
SoftTransitionFormFactor ParseTransitionFormFactor(const std::string &value,
                                                   const std::string &path) {
  if (value == "diagonal") {
    return SoftTransitionFormFactor::Diagonal;
  }
  if (value == "arithmetic") {
    return SoftTransitionFormFactor::Arithmetic;
  }
  if (value == "geometric") {
    return SoftTransitionFormFactor::Geometric;
  }
  throw std::invalid_argument(path + " has unknown transition_ff " + value);
}

// Parse one exact eigenstate form factor family
SoftFormFactor ParseFormFactor(const std::string &value,
                               const std::string &path) {
  if (value == "EXP") {
    return SoftFormFactor::Exponential;
  }
  if (value == "EXPOW") {
    return SoftFormFactor::ExponentialPower;
  }
  if (value == "DPOW") {
    return SoftFormFactor::DualPower;
  }
  if (value == "MIXEXP") {
    return SoftFormFactor::MixedExponential;
  }
  if (value == "ODD3G_NODE") {
    return SoftFormFactor::OddThreeGluonNode;
  }
  if (value == "3G") {
    return SoftFormFactor::ThreeGluon;
  }
  if (value == "GKERNEL") {
    return SoftFormFactor::GeneralizedKernel;
  }
  throw std::invalid_argument(path + " has unknown form factor " + value);
}

// Validate one role and trajectory definition without relying on its name
void ValidateExchangeDefinition(const ExchangeDefinition &definition,
                                const std::string &path) {
  if (regge::Tau(definition.signature) != definition.crossing_parity) {
    throw std::invalid_argument(path +
                                " signature and crossing parity must agree");
  }
  if (definition.role == SoftExchangeRole::Pomeron) {
    if (definition.signature != regge::Signature::Positive) {
      throw std::invalid_argument(path +
                                  " pomeron requires positive signature");
    }
  } else if (definition.role == SoftExchangeRole::Odderon) {
    if (definition.signature != regge::Signature::Negative ||
        definition.trajectory_mode != SoftTrajectoryMode::Linear) {
      throw std::invalid_argument(
          path + " odderon requires negative signature and linear mode");
    }
  } else if (definition.trajectory_mode != SoftTrajectoryMode::Linear) {
    throw std::invalid_argument(path + " reggeon requires linear mode");
  }
}

// Convert and validate one finite square symmetric matrix
MMatrix<double> ReadSymmetricMatrix(const json &value, const std::size_t n,
                                    const std::string &path) {
  if (!value.is_array() || value.size() != n) {
    throw std::invalid_argument(path + " must be an explicit N x N matrix");
  }
  MMatrix<double> output(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    if (!value.at(i).is_array() || value.at(i).size() != n) {
      throw std::invalid_argument(path + " must be an explicit N x N matrix");
    }
    for (std::size_t j = 0; j < n; ++j) {
      if (!value.at(i).at(j).is_number()) {
        throw std::invalid_argument(path + " entries must be numeric");
      }
      output(i, j) = value.at(i).at(j).get<double>();
      if (!std::isfinite(output(i, j))) {
        throw std::invalid_argument(path + " entries must be finite");
      }
    }
  }
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = i + 1; j < n; ++j) {
      const double scale =
          std::max({1.0, std::abs(output(i, j)), std::abs(output(j, i))});
      if (std::abs(output(i, j) - output(j, i)) > 1.0e-12 * scale) {
        throw std::invalid_argument(path + " must be symmetric");
      }
    }
  }
  return output;
}

// Parse one scalar or explicit symmetric transition matrix
MMatrix<double> ReadTransitionMatrix(const json &value, const std::size_t n,
                                     const std::string &path,
                                     const bool require_nonnegative) {
  MMatrix<double> output;
  if (value.is_number()) {
    const double scalar = value.get<double>();
    if (!std::isfinite(scalar)) {
      throw std::invalid_argument(path + " scalar must be finite");
    }
    output = MMatrix<double>(n, n, scalar);
  } else if (value.is_array()) {
    output = ReadSymmetricMatrix(value, n, path);
  } else {
    throw std::invalid_argument(path +
                                " must be a scalar or explicit N x N matrix");
  }

  if (require_nonnegative) {
    for (std::size_t i = 0; i < n; ++i) {
      for (std::size_t j = 0; j < n; ++j) {
        if (output(i, j) < 0.0) {
          throw std::invalid_argument(path + " entries must be nonnegative");
        }
      }
    }
  }
  return output;
}

// Validate one form factor parameter row
void ValidateFormFactorRow(const SoftFormFactor type,
                           const std::vector<double> &row,
                           const std::string &path) {
  const bool size3 = type == SoftFormFactor::DualPower ||
                     type == SoftFormFactor::ExponentialPower ||
                     type == SoftFormFactor::OddThreeGluonNode;
  const bool size4 = type == SoftFormFactor::MixedExponential;
  const bool size1 =
      type == SoftFormFactor::Exponential || type == SoftFormFactor::ThreeGluon;
  const bool kernel_size = type == SoftFormFactor::GeneralizedKernel &&
                           ((!row.empty() && row.size() % 4 == 0) ||
                            (row.size() > 1 && (row.size() - 1) % 4 == 0));
  if (!((size3 && row.size() == 3) || (size4 && row.size() == 4) ||
        (size1 && row.size() == 1) || kernel_size)) {
    throw std::invalid_argument(path + " has an invalid parameter count");
  }
  for (const double parameter : row) {
    if (!std::isfinite(parameter)) {
      throw std::invalid_argument(path + " parameters must be finite");
    }
  }
  if ((type == SoftFormFactor::DualPower ||
       type == SoftFormFactor::ExponentialPower) &&
      (row[0] <= 0.0 || row[1] <= 0.0 || row[2] <= 0.0)) {
    throw std::invalid_argument(path + " requires three positive parameters");
  }
  if (type == SoftFormFactor::Exponential && row[0] <= 0.0) {
    throw std::invalid_argument(path + " requires a positive EXP slope");
  }
  if (type == SoftFormFactor::MixedExponential &&
      (row[0] < 0.0 || row[0] > 1.0 || row[1] < 0.0 || row[2] < 0.0 ||
       row[3] < 0.0)) {
    throw std::invalid_argument(path + " has invalid MIXEXP parameters");
  }
  if (type == SoftFormFactor::OddThreeGluonNode &&
      (row[0] < 0.0 || row[1] <= 0.0 || row[2] <= 0.0)) {
    throw std::invalid_argument(path + " has invalid ODD3G_NODE parameters");
  }
  if (type == SoftFormFactor::ThreeGluon && row[0] <= 0.0) {
    throw std::invalid_argument(path + " requires a positive 3G mass scale");
  }
  if (type == SoftFormFactor::GeneralizedKernel) {
    const std::size_t offset = row.size() % 4 == 1 ? 1 : 0;
    if (offset == 1 && row.front() < 0.0) {
      throw std::invalid_argument(path + " has invalid GKERNEL prefactor");
    }
    for (std::size_t i = offset; i < row.size(); i += 4) {
      if (row[i] < 0.0 || row[i + 1] <= 0.0 || row[i + 2] < 0.0 ||
          row[i + 3] < 0.0) {
        throw std::invalid_argument(path + " has invalid GKERNEL parameters");
      }
      if (math::IsZero(row[i + 3]) && row[i + 1] < 1.0) {
        throw std::invalid_argument(path +
                                    " has a nonfinite GKERNEL forward slope");
      }
    }
  }
}

// Evaluate the common generalized profile with its optional linear prefactor
double GeneralizedKernelProfile(double x, const std::vector<double> &params, const std::string &context) {
  const std::size_t offset    = params.size() % 4 == 1 ? 1 : 0;
  const double      prefactor = offset ? 1.0 + params.front() * x : 1.0;
  if (!std::isfinite(prefactor) || prefactor <= 0.0) { throw AmplitudeFailure(context + ": nonpositive prefactor"); }
  return regge::GKernel(x, std::span<const double>(params).subspan(offset), context, prefactor);
}

// Evaluate the pion loop trajectory function
double PionLoopFunction(const double ratio, const double t,
                        const double scale2) {
  const double mass2 = PDG::mpi * PDG::mpi;
  const double root = std::sqrt(1.0 + ratio);
  const double pion_form_factor = 1.0 / (1.0 - t / scale2);
  return (4.0 / ratio) * pion_form_factor * pion_form_factor *
         (2.0 * ratio -
          std::pow(1.0 + ratio, 1.5) * std::log((root + 1.0) / (root - 1.0)) +
          std::log(1.0 / mass2));
}

// Parse and validate every soft exchange definition
std::vector<ExchangeDefinition>
ReadExchangeDefinitions(const json &definition_block) {
  if (!definition_block.is_object() || definition_block.empty()) {
    throw std::invalid_argument(
        "PARAM_SOFT.EXCHANGE_DEF must be a nonempty object");
  }
  std::vector<ExchangeDefinition> definitions;
  definitions.reserve(definition_block.size());
  for (auto it = definition_block.begin(); it != definition_block.end(); ++it) {
    const std::string path = "PARAM_SOFT.EXCHANGE_DEF." + it.key();
    const json &entry = it.value();
    if (!entry.is_object()) {
      throw std::invalid_argument(path + " must be an object");
    }
    if (entry.size() != 4 || !entry.contains("role") ||
        !entry.contains("trajectory_mode") || !entry.contains("tau") ||
        !entry.contains("crossing")) {
      throw std::invalid_argument(
          path + " must contain only role, trajectory_mode, tau and crossing");
    }
    ExchangeDefinition definition;
    definition.name = it.key();
    definition.role =
        ParseExchangeRole(entry.at("role").get<std::string>(), path);
    definition.trajectory_mode = ParseTrajectoryMode(
        entry.at("trajectory_mode").get<std::string>(), path);
    if (!entry.contains("tau")) {
      throw std::invalid_argument(path + ".tau must be -1 or 1");
    }
    definition.signature = regge::ParseSignature(
        ParseUnitSign(entry.at("tau"), path + ".tau"), path + ".tau");
    definition.crossing_parity =
        ParseUnitSign(entry.at("crossing"), path + ".crossing");
    ValidateExchangeDefinition(definition, path);
    definitions.push_back(std::move(definition));
  }
  return definitions;
}

// Require exact agreement between exchange definitions and model rows
void ValidateExchangeRows(const std::vector<ExchangeDefinition> &definitions,
                          const json &definition_block,
                          const json &exchange_block,
                          const std::string &active_model) {
  if (!exchange_block.is_object()) {
    throw std::invalid_argument("PARAM_SOFT.MODEL." + active_model +
                                ".EXCHANGE must be an object");
  }
  for (const auto &definition : definitions) {
    if (!exchange_block.contains(definition.name)) {
      throw std::invalid_argument("PARAM_SOFT.MODEL." + active_model +
                                  ".EXCHANGE is missing " + definition.name);
    }
  }
  for (auto it = exchange_block.begin(); it != exchange_block.end(); ++it) {
    if (!definition_block.contains(it.key())) {
      throw std::invalid_argument("PARAM_SOFT.MODEL." + active_model +
                                  ".EXCHANGE has unknown exchange " + it.key());
    }
  }
}

// Determine the Good Walker dimension from the first Pomeron coupling matrix
std::size_t ReadChannelCount(const std::vector<ExchangeDefinition> &definitions,
                             const json &exchange_block) {
  for (const auto &definition : definitions) {
    if (definition.role != SoftExchangeRole::Pomeron) {
      continue;
    }
    const json &coupling = exchange_block.at(definition.name).at("g");
    if (!coupling.is_array() || coupling.empty()) {
      throw std::invalid_argument("Pomeron coupling matrix must not be empty");
    }
    return coupling.size();
  }
  throw std::invalid_argument(
      "PARAM_SOFT requires at least one pomeron exchange");
}

// Parse all named eigenstate form factor parameter banks
std::map<std::string, FormFactorBank>
ReadFormFactorBanks(const json &model, const std::string &active_model,
                    const std::size_t channel_count) {
  const json &form_factor_block = model.at("FF");
  if (!form_factor_block.is_object() || form_factor_block.empty()) {
    throw std::invalid_argument("PARAM_SOFT model FF must not be empty");
  }
  std::map<std::string, FormFactorBank> form_factors;
  for (auto it = form_factor_block.begin(); it != form_factor_block.end();
       ++it) {
    const std::string path =
        "PARAM_SOFT.MODEL." + active_model + ".FF." + it.key();
    FormFactorBank bank;
    bank.type = ParseFormFactor(it.value().at("type").get<std::string>(), path);
    bank.parameters =
        it.value().at("param").get<std::vector<std::vector<double>>>();
    if (bank.parameters.size() != channel_count) {
      throw std::invalid_argument(path + " must contain N parameter rows");
    }
    for (const auto &row : bank.parameters) {
      ValidateFormFactorRow(bank.type, row, path);
    }
    form_factors.emplace(it.key(), std::move(bank));
  }
  return form_factors;
}

// Parse the Good Walker mixing geometry and resolved direction
GoodWalkerSpace ReadGoodWalkerSpace(const json &model,
                                    const std::size_t channel_count) {
  const json &good_walker = model.at("GW");
  return GoodWalkerSpace(channel_count,
                         good_walker.at("theta").get<std::vector<double>>(),
                         good_walker.at("a_c").get<std::vector<double>>());
}

// Parse the proton helicity transition switch from the eikonal controls
HelicityDefinition ReadHelicityDefinition(const json &model) {
  HelicityDefinition definition;
  const json &eikonal = model.at("EIKONAL");
  if (!eikonal.contains("helicity") || !eikonal.at("helicity").is_boolean()) {
    throw std::invalid_argument(
        "PARAM_SOFT.EIKONAL.helicity must be a boolean");
  }
  definition.enabled = eikonal.at("helicity").get<bool>();
  definition.mass_scale = PDG::mp;
  return definition;
}

// Parse one complete set of named soft exchange parameter rows
ParsedExchanges
ReadSoftExchanges(const std::vector<ExchangeDefinition> &definitions,
                  const std::map<std::string, FormFactorBank> &form_factors,
                  const json &exchange_block, const std::string &active_model,
                  const std::size_t channel_count,
                  const bool helicity_enabled) {
  ParsedExchanges output;
  output.values.reserve(definitions.size());
  for (const auto &definition : definitions) {
    const json &entry = exchange_block.at(definition.name);
    const std::string path =
        "PARAM_SOFT.MODEL." + active_model + ".EXCHANGE." + definition.name;
    if (!entry.is_object()) {
      throw std::invalid_argument(path + " must be an object");
    }
    if (!entry.contains("on") || !entry.at("on").is_boolean()) {
      throw std::invalid_argument(path + ".on must be a boolean");
    }
    const bool enabled = entry.at("on").get<bool>();
    const auto alpha = entry.at("alpha").get<std::vector<double>>();
    if (alpha.size() != 2) {
      throw std::invalid_argument(path + ".alpha requires two parameters");
    }
    if (!entry.contains("eta_mode") || !entry.at("eta_mode").is_string()) {
      throw std::invalid_argument(path + ".eta_mode must be a string");
    }
    const EtaMode eta_mode = regge::ParseEta(
        entry.at("eta_mode").get<std::string>(), path + ".eta_mode");
    regge::CheckEta(alpha[0], alpha[1], definition.signature, eta_mode, path);
    const MMatrix<double> coupling =
        ReadSymmetricMatrix(entry.at("g"), channel_count, path + ".g");
    const int sign = ParseUnitSign(entry.at("sign"), path + ".sign");
    const std::string form_factor_name = entry.at("ff").get<std::string>();
    const auto form_factor = form_factors.find(form_factor_name);
    if (form_factor == form_factors.end()) {
      throw std::invalid_argument(path + ".ff selects an unknown bank");
    }
    if (definition.role == SoftExchangeRole::Pomeron &&
        (form_factor->second.type == SoftFormFactor::OddThreeGluonNode ||
         form_factor->second.type == SoftFormFactor::ThreeGluon)) {
      throw std::invalid_argument(
          path + " pomeron cannot use a three gluon form factor");
    }

    MMatrix<double> flip_coupling(channel_count, channel_count, 0.0);
    MMatrix<double> flip_slope(channel_count, channel_count, 0.0);
    if (helicity_enabled) {
      const json &helicity = entry.at("helicity");
      flip_coupling = ReadTransitionMatrix(helicity.at("kappa"), channel_count,
                                           path + ".helicity.kappa", false);
      flip_slope = ReadTransitionMatrix(helicity.at("B_kappa"), channel_count,
                                        path + ".helicity.B_kappa", true);
    }

    const SoftExchangeId id(output.values.size());
    output.ids.emplace(definition.name, id);
    output.values.emplace_back(
        definition.name, definition.role, definition.trajectory_mode, enabled,
        definition.crossing_parity, definition.signature, alpha[0], alpha[1],
        coupling, sign, eta_mode, form_factor_name, form_factor->second.type,
        form_factor->second.parameters,
        ParseTransitionFormFactor(entry.at("transition_ff").get<std::string>(),
                                  path),
        std::move(flip_coupling), std::move(flip_slope));
  }
  return output;
}

// Expand one explicit selection or wildcard over enabled exchanges
std::vector<SoftExchangeId>
ReadExchangeSelection(const json &selection, const std::string &path,
                      const std::vector<SoftExchange> &exchanges,
                      const std::map<std::string, SoftExchangeId> &exchange_ids,
                      const bool pomeron_only) {
  if (!selection.is_array() || selection.empty()) {
    throw std::invalid_argument(path + " must be a nonempty array");
  }
  if (selection.size() == 1 && selection.at(0).is_string() &&
      selection.at(0).get<std::string>() == "*") {
    std::vector<SoftExchangeId> selected;
    for (const auto &index : indices(exchanges)) {
      const auto &exchange = exchanges[index];
      if (exchange.Enabled() &&
          (!pomeron_only || exchange.Role() == SoftExchangeRole::Pomeron)) {
        selected.emplace_back(index);
      }
    }
    if (selected.empty()) {
      throw std::invalid_argument(path + " wildcard selects no exchanges");
    }
    return selected;
  }

  std::vector<SoftExchangeId> selected;
  std::set<std::size_t> selected_indices;
  for (const auto &value : selection) {
    if (!value.is_string()) {
      throw std::invalid_argument(path + " must contain exchange names");
    }
    const std::string name = value.get<std::string>();
    if (name == "*") {
      throw std::invalid_argument(path + " wildcard must be used alone");
    }
    const auto found = exchange_ids.find(name);
    if (found == exchange_ids.end()) {
      throw std::invalid_argument(path + " selects unknown exchange " + name);
    }
    const auto &exchange = exchanges.at(found->second.Value());
    if (!exchange.Enabled()) {
      throw std::invalid_argument(path + " selects disabled exchange " + name);
    }
    if (pomeron_only && exchange.Role() != SoftExchangeRole::Pomeron) {
      throw std::invalid_argument(path + " requires pomeron exchanges");
    }
    if (!selected_indices.insert(found->second.Value()).second) {
      throw std::invalid_argument(path + " duplicates exchange " + name);
    }
    selected.push_back(found->second);
  }
  return selected;
}

// Parse the matrix eikonal and screening exchange selection
SoftEikonalSettings
ReadEikonalSettings(const json &model, const HelicityDefinition &helicity,
                    const std::vector<SoftExchange> &exchanges,
                    const std::map<std::string, SoftExchangeId> &exchange_ids) {
  const json &eikonal = model.at("EIKONAL");
  const std::string unitarization_name =
      eikonal.at("unitarization").get<std::string>();
  const SoftUnitarization unitarization =
      unitarization_name == "exp" ? SoftUnitarization::Exponential
      : unitarization_name == "q_exp"
          ? SoftUnitarization::QExponential
          : throw std::invalid_argument(
                "PARAM_SOFT.EIKONAL has unknown unitarization");
  const std::string screening_path = "PARAM_SOFT.EIKONAL.screening_exchanges";
  if (!eikonal.contains("screening_exchanges")) {
    throw std::invalid_argument(screening_path + " must be a nonempty array");
  }
  const auto screening_exchanges =
      ReadExchangeSelection(eikonal.at("screening_exchanges"), screening_path,
                            exchanges, exchange_ids, false);
  return SoftEikonalSettings(unitarization, eikonal.at("q").get<double>(),
                             helicity.enabled, helicity.mass_scale,
                             screening_exchanges);
}

// Parse the single Pomeron selected by precontracted forward excitation
// amplitudes
SoftExchangeId ReadForwardExcitationExchange(
    const json &model, const std::string &active_model,
    const std::vector<SoftExchange> &exchanges,
    const std::map<std::string, SoftExchangeId> &exchange_ids) {
  const std::string path =
      "PARAM_SOFT.MODEL." + active_model + ".EIKONAL.excitation_exchanges";
  const json &eikonal = model.at("EIKONAL");
  if (!eikonal.contains("excitation_exchanges")) {
    throw std::invalid_argument(path + " is required");
  }
  const auto selected = ReadExchangeSelection(
      eikonal.at("excitation_exchanges"), path, exchanges, exchange_ids, true);
  if (selected.size() != 1) {
    throw std::invalid_argument(
        path + " must resolve to exactly one pomeron exchange");
  }
  return selected.front();
}

// Parse one positive pion loop scale
double ReadPionLoopScale2(const json &model) {
  const double value = model.at("pion_loop_scale2").get<double>();
  if (!std::isfinite(value) || value <= 0.0) {
    throw std::invalid_argument("PARAM_SOFT.pion_loop_scale2 must be positive");
  }
  return value;
}

// Parse the raw mapped triple Pomeron ratio and eta prescription
TriplePomeronDefinition ReadTriplePomeronDefinition(const json &model) {
  const json &triple = model.at("3P");
  TriplePomeronDefinition definition;
  definition.ratio = triple.at("g").get<double>();
  if (!std::isfinite(definition.ratio) || definition.ratio <= 0.0) {
    throw std::invalid_argument("PARAM_SOFT.3P.g ratio must be positive");
  }
  if (!triple.contains("eta_mode") || !triple.at("eta_mode").is_string()) {
    throw std::invalid_argument("PARAM_SOFT.3P.eta_mode must be a string");
  }
  definition.eta_mode = regge::ParseEta(
      triple.at("eta_mode").get<std::string>(), "PARAM_SOFT.3P.eta_mode");
  return definition;
}

// Build a deterministic diagnostic fingerprint from physics-owned blocks
std::string BuildFingerprint(const json &definition_block, const json &model,
                             const json &inelastic) {
  const std::string physics_json =
      definition_block.dump() + model.dump() + inelastic.dump();
  std::ostringstream fingerprint;
  fingerprint << std::hex << aux::djb2hash(physics_json);
  return fingerprint.str();
}

} // namespace

// Construct and validate one arbitrary-channel Good Walker space
GoodWalkerSpace::GoodWalkerSpace(const std::size_t channel_count,
                                 std::vector<double> mixing_angles,
                                 std::vector<double> resolved_coefficients)
    : channel_count_(channel_count), mixing_angles_(std::move(mixing_angles)),
      resolved_coefficients_(std::move(resolved_coefficients)) {
  if (channel_count_ == 0) {
    throw std::invalid_argument(
        "GoodWalkerSpace: channel count must be positive");
  }
  if (channel_count_ >
      std::numeric_limits<std::size_t>::max() / channel_count_) {
    throw std::invalid_argument(
        "GoodWalkerSpace: pair-space dimension overflows size_t");
  }
  pair_dimension_ = channel_count_ * channel_count_;
  if (resolved_coefficients_.size() + 1 != channel_count_) {
    throw std::invalid_argument(
        "GoodWalkerSpace: resolved coefficient count must equal N - 1");
  }

  if (!gra::AllFinite(resolved_coefficients_)) {
    throw std::invalid_argument(
        "GoodWalkerSpace: resolved coefficients must be finite");
  }
  // Normalize signed amplitude coordinates independently of their overall scale
  if (channel_count_ > 1) {
    resolved_coefficients_ = gra::NormalizedL2(resolved_coefficients_);
  }

  mixing_ = MMatrix<double>::MixingReal(mixing_angles_, channel_count_);
  proton_basis_ = mixing_.Submatrix(0, 0, 1, channel_count_).Transpose();
  excited_basis_ =
      mixing_.Submatrix(1, 0, channel_count_ - 1, channel_count_).Transpose();
  complete_basis_ = MMatrix<double>::IdentityMatrix(channel_count_);
  proton_ = proton_basis_.Flatten();
  resolved_ = excited_basis_ * resolved_coefficients_;

  const double proton_norm2 = gra::InnerProduct(proton_, proton_);
  const double resolved_norm2 = gra::InnerProduct(resolved_, resolved_);
  const double overlap = gra::InnerProduct(proton_, resolved_);
  if (std::abs(proton_norm2 - 1.0) > 1.0e-12 ||
      (channel_count_ > 1 && std::abs(resolved_norm2 - 1.0) > 1.0e-12) ||
      std::abs(overlap) > 1.0e-12) {
    throw std::logic_error(
        "GoodWalkerSpace: physical basis construction failed");
  }

  proton_projector_ = gra::RankOneProjector(proton_);
  resolved_projector_ = gra::RankOneProjector(resolved_);
  const MMatrix<double> identity =
      MMatrix<double>::IdentityMatrix(channel_count_);
  excited_projector_ = identity - proton_projector_;
  inclusive_projector_ = identity - resolved_projector_;
}

// Compute the row-major pair-space index for eigenstates i and k
std::size_t GoodWalkerSpace::PairIndex(const std::size_t i,
                                       const std::size_t k) const {
  if (i >= channel_count_ || k >= channel_count_) {
    throw std::out_of_range("GoodWalkerSpace: pair index is out of range");
  }
  return i * channel_count_ + k;
}

// Compute one cached orthonormal physical final-state basis
const MMatrix<double> &
GoodWalkerSpace::FinalBasis(const GoodWalkerFinalBasis basis) const {
  if (basis == GoodWalkerFinalBasis::Proton) {
    return proton_basis_;
  }
  if (basis == GoodWalkerFinalBasis::Excited) {
    return excited_basis_;
  }
  if (basis == GoodWalkerFinalBasis::Complete) {
    return complete_basis_;
  }
  throw std::invalid_argument("GoodWalkerSpace: unknown final basis");
}

// Project one pair-space source onto two cached physical final-state bases
std::vector<std::complex<double>>
GoodWalkerSpace::ProjectPair(const std::span<const std::complex<double>> source,
                             const GoodWalkerFinalBasis upper,
                             const GoodWalkerFinalBasis lower) const {
  return ProjectPair(source, FinalBasis(upper), FinalBasis(lower));
}

// Compute B_upper^T A B_lower for real Good Walker column bases
std::vector<std::complex<double>>
GoodWalkerSpace::ProjectPair(const std::span<const std::complex<double>> source,
                             const MMatrix<double> &upper_basis,
                             const MMatrix<double> &lower_basis) const {
  if (source.size() != PairDimension() ||
      upper_basis.size_row() != channel_count_ ||
      lower_basis.size_row() != channel_count_) {
    throw std::invalid_argument(
        "GoodWalkerSpace::ProjectPair: pair-space dimensions disagree");
  }
  return gra::TensorProductProjection(source, upper_basis, lower_basis);
}

// Compute the steering name of one soft exchange role
std::string SoftExchangeRoleName(const SoftExchangeRole role) {
  if (role == SoftExchangeRole::Pomeron) {
    return "pomeron";
  }
  if (role == SoftExchangeRole::Reggeon) {
    return "reggeon";
  }
  if (role == SoftExchangeRole::Odderon) {
    return "odderon";
  }
  throw std::invalid_argument("SoftExchangeRoleName: unknown role");
}

// Compute the steering name of one trajectory mode
std::string SoftTrajectoryModeName(const SoftTrajectoryMode mode) {
  if (mode == SoftTrajectoryMode::Linear) {
    return "linear";
  }
  if (mode == SoftTrajectoryMode::PionLoop) {
    return "pion_loop";
  }
  throw std::invalid_argument("SoftTrajectoryModeName: unknown mode");
}

// Construct one validated soft exchange value
SoftExchange::SoftExchange(
    std::string name, const SoftExchangeRole role,
    const SoftTrajectoryMode trajectory_mode, const bool enabled,
    const int crossing_parity, const regge::Signature signature,
    const double alpha0, const double alpha_prime, MMatrix<double> coupling,
    const int residue_sign, const EtaMode eta_mode,
    std::string form_factor_name, const SoftFormFactor form_factor,
    std::vector<std::vector<double>> form_factor_parameters,
    const SoftTransitionFormFactor transition_form_factor,
    MMatrix<double> helicity_flip_coupling, MMatrix<double> helicity_flip_slope)
    : name_(std::move(name)), role_(role), trajectory_mode_(trajectory_mode),
      enabled_(enabled), crossing_parity_(crossing_parity),
      signature_(signature), alpha0_(alpha0), alpha_prime_(alpha_prime),
      coupling_(std::move(coupling)), residue_sign_(residue_sign),
      eta_mode_(eta_mode), form_factor_name_(std::move(form_factor_name)),
      form_factor_(form_factor),
      form_factor_parameters_(std::move(form_factor_parameters)),
      transition_form_factor_(transition_form_factor),
      helicity_flip_coupling_(std::move(helicity_flip_coupling)),
      helicity_flip_slope_(std::move(helicity_flip_slope)) {
  const std::size_t n = coupling_.size_row();
  const bool invalid_signature = regge::Tau(signature_) != crossing_parity_;
  if (name_.empty() || n == 0 || coupling_.size_col() != n ||
      form_factor_parameters_.size() != n ||
      helicity_flip_coupling_.size_row() != n ||
      helicity_flip_coupling_.size_col() != n ||
      helicity_flip_slope_.size_row() != n ||
      helicity_flip_slope_.size_col() != n || invalid_signature ||
      (residue_sign_ != -1 && residue_sign_ != 1)) {
    throw std::invalid_argument("SoftExchange: inconsistent exchange data");
  }
}

// Compute one Good Walker transition coupling
double SoftExchange::Coupling(const std::size_t i, const std::size_t j) const {
  return coupling_(i, j);
}

// Compute one Pauli flip coupling
double SoftExchange::HelicityFlipCoupling(const std::size_t i,
                                          const std::size_t j) const {
  return helicity_flip_coupling_(i, j);
}

// Compute one Pauli flip slope
double SoftExchange::HelicityFlipSlope(const std::size_t i,
                                       const std::size_t j) const {
  return helicity_flip_slope_(i, j);
}

// Construct one validated set of matrix eikonal controls
SoftEikonalSettings::SoftEikonalSettings(
    const SoftUnitarization unitarization, const double q,
    const bool helicity_enabled, const double helicity_mass_scale,
    std::vector<SoftExchangeId> screening_exchanges)
    : unitarization_(unitarization), q_(q), helicity_enabled_(helicity_enabled),
      helicity_mass_scale_(helicity_mass_scale),
      screening_exchanges_(std::move(screening_exchanges)) {
  if (!std::isfinite(q_) || q_ <= 0.0 || q_ > 1.0 ||
      !std::isfinite(helicity_mass_scale_) || helicity_mass_scale_ <= 0.0 ||
      screening_exchanges_.empty()) {
    throw std::invalid_argument(
        "SoftEikonalSettings: invalid eikonal controls");
  }
  if (unitarization_ == SoftUnitarization::Exponential &&
      std::abs(q_ - 1.0) > 1.0e-12) {
    throw std::invalid_argument(
        "SoftEikonalSettings: exponential unitarization requires q = 1");
  }
}

// Construct one validated inelastic proton transition profile
SoftForwardExcitationProfile::SoftForwardExcitationProfile(const double s0,
                                                           const double a)
    : s0_(s0), a_(a) {
  if (!std::isfinite(s0_) || s0_ <= 0.0 || !std::isfinite(a_) || a_ <= 0.0) {
    throw std::invalid_argument(
        "SoftForwardExcitationProfile: s0 and a must be positive");
  }
}

// Construct one already parsed and validated immutable model snapshot
SoftModel::SoftModel(std::string source_file, std::string source_json,
                     std::string numerics_source_file,
                     std::string numerics_source_json, std::string active_model,
                     std::string fingerprint,
                     std::vector<SoftExchange> exchanges,
                     std::map<std::string, SoftExchangeId> exchange_ids,
                     GoodWalkerSpace good_walker, SoftEikonalSettings eikonal,
                     const SoftExchangeId forward_excitation_exchange,
                     const double pion_loop_scale2,
                     const double triple_pomeron_ratio,
                     const EtaMode triple_pomeron_eta_mode,
                     SoftForwardExcitationProfile forward_excitation)
    : source_file_(std::move(source_file)),
      source_json_(std::move(source_json)),
      numerics_source_file_(std::move(numerics_source_file)),
      numerics_source_json_(std::move(numerics_source_json)),
      active_model_(std::move(active_model)),
      fingerprint_(std::move(fingerprint)), exchanges_(std::move(exchanges)),
      exchange_ids_(std::move(exchange_ids)),
      good_walker_(std::move(good_walker)), eikonal_(std::move(eikonal)),
      forward_excitation_exchange_(forward_excitation_exchange),
      pion_loop_scale2_(pion_loop_scale2),
      triple_pomeron_ratio_(triple_pomeron_ratio),
      triple_pomeron_eta_mode_(triple_pomeron_eta_mode),
      forward_excitation_(forward_excitation) {}

// Parse one complete model snapshot from GENERAL JSON text
std::shared_ptr<const SoftModel>
SoftModel::LoadFromJson(const std::string &source_file,
                        const std::string &json_text) {
  return LoadFromJson(source_file, json_text, "", "");
}

// Parse one model snapshot with GENERAL and NUMERICS JSON text
std::shared_ptr<const SoftModel> SoftModel::LoadFromJson(
    const std::string &source_file, const std::string &json_text,
    const std::string &numerics_source_file, const std::string &numerics_json) {
  try {
    const json root = json::parse(json_text);
    const json &soft = root.at("PARAM_SOFT");
    const std::string active_model = soft.at("active_model").get<std::string>();
    const json &model = soft.at("MODEL").at(active_model);
    const json &definition_block = soft.at("EXCHANGE_DEF");
    const json &exchange_block = model.at("EXCHANGE");
    const std::vector<ExchangeDefinition> definitions =
        ReadExchangeDefinitions(definition_block);
    ValidateExchangeRows(definitions, definition_block, exchange_block,
                         active_model);
    const std::size_t channel_count =
        ReadChannelCount(definitions, exchange_block);
    const auto form_factors =
        ReadFormFactorBanks(model, active_model, channel_count);
    GoodWalkerSpace good_walker = ReadGoodWalkerSpace(model, channel_count);
    const HelicityDefinition helicity = ReadHelicityDefinition(model);
    ParsedExchanges exchanges =
        ReadSoftExchanges(definitions, form_factors, exchange_block,
                          active_model, channel_count, helicity.enabled);
    SoftEikonalSettings eikonal_settings =
        ReadEikonalSettings(model, helicity, exchanges.values, exchanges.ids);
    const SoftExchangeId forward_excitation_exchange =
        ReadForwardExcitationExchange(model, active_model, exchanges.values,
                                      exchanges.ids);
    const double pion_loop_scale2 = ReadPionLoopScale2(model);
    const TriplePomeronDefinition triple = ReadTriplePomeronDefinition(model);
    const json &inelastic = soft.at("FORWARD_EXCITATION");
    SoftForwardExcitationProfile forward_excitation(
        ReadForwardExcitationScale(inelastic, "s0"),
        ReadForwardExcitationScale(inelastic, "a"));

    auto snapshot = std::shared_ptr<const SoftModel>(new SoftModel(
        source_file, json_text, numerics_source_file, numerics_json,
        active_model, BuildFingerprint(definition_block, model, inelastic),
        std::move(exchanges.values), std::move(exchanges.ids),
        std::move(good_walker), std::move(eikonal_settings),
        forward_excitation_exchange, pion_loop_scale2, triple.ratio,
        triple.eta_mode, forward_excitation));

    return snapshot;
  } catch (const std::exception &error) {
    throw std::invalid_argument("SoftModel::LoadFromJson: Error reading " +
                                source_file + ": " + error.what());
  }
}

// Resolve one named exchange inside this snapshot
SoftExchangeId SoftModel::ExchangeId(const std::string &name) const {
  const auto found = exchange_ids_.find(name);
  if (found == exchange_ids_.end()) {
    throw std::invalid_argument("SoftModel exchange '" + name +
                                "' is not configured");
  }
  return found->second;
}

// Compute one exchange inside this snapshot
const SoftExchange &SoftModel::Exchange(const SoftExchangeId exchange) const {
  if (exchange.Value() >= exchanges_.size()) {
    throw std::out_of_range("SoftModel exchange index is out of range");
  }
  return exchanges_[exchange.Value()];
}

// Evaluate one mapped linear or pion loop trajectory
// alpha(t) = alpha0+alpha' t plus the optional pion loop correction
double SoftModel::Alpha(const SoftExchangeId exchange, const double t) const {
  if (!std::isfinite(t)) {
    throw AmplitudeFailure("SoftModel::Alpha: t must be finite");
  }
  const SoftExchange &parameters = Exchange(exchange);
  const double linear = parameters.Alpha0() + parameters.AlphaPrime() * t;
  if (parameters.TrajectoryMode() == SoftTrajectoryMode::Linear ||
      std::abs(t) <= 1.0e-12) {
    return linear;
  }
  const double beta_pi = (2.0 / 3.0) * PhysicalCoupling(exchange);
  const double coefficient =
      beta_pi * beta_pi * PDG::mpi * PDG::mpi / (32.0 * std::pow(math::PI, 3));
  const double ratio = 4.0 * PDG::mpi * PDG::mpi / std::abs(t);
  return linear - coefficient * PionLoopFunction(ratio, t, pion_loop_scale2_);
}

// Evaluate one diagonal eigenstate form factor
double SoftModel::FormFactor(const SoftExchangeId exchange, const double t,
                             const std::size_t i) const {
  if (!std::isfinite(t)) {
    throw AmplitudeFailure("SoftModel::FormFactor: t must be finite");
  }
  const SoftExchange &parameters = Exchange(exchange);
  const auto &rows = parameters.FormFactorParameters();
  if (i >= rows.size()) {
    throw std::out_of_range(
        "SoftModel::FormFactor: eigenstate index is out of range");
  }
  const auto &row = rows[i];
  if (parameters.FormFactorType() == SoftFormFactor::Exponential) {
    return std::exp(0.5 * row[0] * t);
  }
  if (parameters.FormFactorType() == SoftFormFactor::ExponentialPower) {
    return std::exp(-std::pow(row[0] * (row[1] - t), row[2]) +
                    std::pow(row[0] * row[1], row[2]));
  }
  if (parameters.FormFactorType() == SoftFormFactor::DualPower) {
    return std::pow(1.0 / ((1.0 - t / row[0]) * (1.0 - t / row[1])), row[2]);
  }
  if (parameters.FormFactorType() == SoftFormFactor::MixedExponential) {
    const double profile = (1.0 - row[0]) * std::exp(0.5 * row[1] * t) +
                           row[0] * std::exp(0.5 * row[2] * t);
    return row[3] > 0.0 ? (1.0 + t / row[3]) * profile : profile;
  }
  if (parameters.FormFactorType() == SoftFormFactor::OddThreeGluonNode) {
    return std::exp(0.5 * row[0] * t) * (1.0 + t / row[1]) /
           std::pow(1.0 - t / row[2], 3);
  }
  if (parameters.FormFactorType() == SoftFormFactor::ThreeGluon) {
    return 1.0 / std::pow(1.0 - t / row[0], 2);
  }
  return GeneralizedKernelProfile(
      -t, row, "SoftModel::FormFactor " + parameters.FormFactorName());
}

// Evaluate one symmetric eigenstate-transition form factor
// F_ij = delta_ij F_i, (F_i+F_j)/2, or sign(F_i)sqrt(F_i F_j)
double SoftModel::TransitionFormFactor(const SoftExchangeId exchange,
                                       const double t, const std::size_t i,
                                       const std::size_t j) const {
  const SoftExchange &parameters = Exchange(exchange);
  const double first = FormFactor(exchange, t, i);
  if (i == j) {
    return first;
  }
  if (parameters.TransitionFormFactor() == SoftTransitionFormFactor::Diagonal) {
    return 0.0;
  }
  const double second = FormFactor(exchange, t, j);
  if (parameters.TransitionFormFactor() ==
      SoftTransitionFormFactor::Arithmetic) {
    return 0.5 * (first + second);
  }
  const double product = first * second;
  if (product < -1.0e-14) {
    throw AmplitudeFailure(
        "SoftModel::TransitionFormFactor: geometric profiles have opposite "
        "signs");
  }
  // The geometric transition is real after the sign check above
  const double magnitude = std::sqrt(std::max(0.0, product));
  return (first < 0.0 && second < 0.0) ? -magnitude : magnitude;
}

// Compute one complete exchange-vertex matrix at fixed momentum transfer
// beta_ij(t) = g_ij F_ij(t)
MMatrix<double> SoftModel::ResidueMatrix(const SoftExchangeId exchange,
                                         const double t) const {
  const std::size_t n = good_walker_.ChannelCount();
  MMatrix<double> residue(n, n, 0.0);
  for (std::size_t i = 0; i < n; ++i) {
    for (std::size_t j = 0; j < n; ++j) {
      residue(i, j) = Exchange(exchange).Coupling(i, j) *
                      TransitionFormFactor(exchange, t, i, j);
    }
  }
  return residue;
}

// Compute one complete Pauli flip vertex matrix at fixed momentum transfer
MMatrix<double>
SoftModel::HelicityFlipResidueMatrix(const SoftExchangeId exchange,
                                     const double t) const {
  MMatrix<double> residue = ResidueMatrix(exchange, t);
  const SoftExchange &parameters = Exchange(exchange);
  for (std::size_t i = 0; i < residue.size_row(); ++i) {
    for (std::size_t j = 0; j < residue.size_col(); ++j) {
      residue[i][j] *= parameters.HelicityFlipCoupling(i, j) *
                       std::exp(0.5 * parameters.HelicityFlipSlope(i, j) * t);
    }
  }
  return residue;
}

// Evaluate only the diagonal Good Walker vertices into caller storage
void SoftModel::DiagonalResidues(const SoftExchangeId exchange, const double t,
                                 const std::span<double> output) const {
  const std::size_t n = good_walker_.ChannelCount();
  if (output.size() != n) {
    throw std::invalid_argument(
        "SoftModel::DiagonalResidues: output dimension disagrees");
  }
  const SoftExchange &parameters = Exchange(exchange);
  for (std::size_t i = 0; i < n; ++i) {
    output[i] = parameters.Coupling(i, i) * FormFactor(exchange, t, i);
  }
}

// Compute one exchange coupling projected onto the physical proton
double SoftModel::PhysicalCoupling(const SoftExchangeId exchange) const {
  const auto &proton = good_walker_.ProtonVector();
  const SoftExchange &parameters = Exchange(exchange);
  return parameters.CouplingMatrix().MatrixElement(proton, proton);
}

// Compute one exchange vertex projected onto the physical proton
double SoftModel::PhysicalResidue(const SoftExchangeId exchange,
                                  const double t) const {
  const auto &proton = good_walker_.ProtonVector();
  const MMatrix<double> residue = ResidueMatrix(exchange, t);
  return residue.MatrixElement(proton, proton);
}

// Compute the physical exchange form factor normalized at zero transfer
double SoftModel::NormalizedPhysicalResidue(const SoftExchangeId exchange,
                                            const double t) const {
  const double reference = PhysicalResidue(exchange, 0.0);
  const double residue = PhysicalResidue(exchange, t);
  if (!std::isfinite(reference) || math::IsZero(reference) ||
      !std::isfinite(residue)) {
    throw std::invalid_argument(
        "SoftModel normalized physical residue requires a finite nonzero "
        "zero-transfer residue");
  }
  const double value = residue / reference;
  if (!std::isfinite(value)) {
    throw std::invalid_argument(
        "SoftModel normalized physical residue is not finite");
  }
  return value;
}

// Compute the effective mapped triple-Pomeron coupling
// g_3P = r_3P beta_p(0)
double
SoftModel::EffectiveTriplePomeronCoupling(const SoftExchangeId exchange) const {
  const SoftExchange &parameters = Exchange(exchange);
  if (!parameters.Enabled() || parameters.Role() != SoftExchangeRole::Pomeron) {
    throw std::invalid_argument(
        "SoftModel triple Pomeron exchange must be an enabled pomeron");
  }
  const double value = triple_pomeron_ratio_ * PhysicalCoupling(exchange);
  if (!std::isfinite(value) || value < 0.0) {
    throw std::invalid_argument(
        "SoftModel effective triple Pomeron coupling must be nonnegative");
  }
  return value;
}

// Compute the real symmetric principal triple Pomeron matrix square root
MMatrix<double>
SoftModel::TriplePomeronCouplingRoot(const SoftExchangeId exchange) const {
  const SoftExchange &parameters = Exchange(exchange);
  const double coupling = EffectiveTriplePomeronCoupling(exchange);
  regge::CheckEta(parameters.Alpha0(), parameters.AlphaPrime(),
                  parameters.Signature(), triple_pomeron_eta_mode_,
                  "PARAM_SOFT.3P exchange " + parameters.Name());
  const MMatrix<double> target = parameters.CouplingMatrix() * coupling;
  try {
    return target.PrincipalPositiveSemidefiniteSquareRoot();
  } catch (const std::domain_error &) {
    throw std::invalid_argument(
        "SoftModel triple Pomeron coupling matrix is not positive "
        "semidefinite");
  } catch (const std::runtime_error &) {
    throw std::runtime_error(
        "SoftModel triple Pomeron matrix square root reconstruction failed");
  }
}

// Evaluate the full-range mapped inelastic proton-transition profile
// factor = [s0 |t|/(M2(|t|+a))]^(alpha0/2)
double SoftModel::ForwardExcitationFactor(const SoftExchangeId exchange,
                                          const double t,
                                          const double mass2) const {
  const SoftExchange &parameters = Exchange(exchange);
  if (!parameters.Enabled() || parameters.Role() != SoftExchangeRole::Pomeron) {
    throw std::invalid_argument(
        "SoftModel forward excitation requires an enabled pomeron exchange");
  }
  if (!std::isfinite(t) || !std::isfinite(mass2) || mass2 <= 0.0) {
    throw AmplitudeFailure(
        "SoftModel forward excitation requires finite t and positive mass2");
  }
  const double abs_t = std::abs(t);
  const double base = (forward_excitation_.S0() * abs_t) /
                      (mass2 * (abs_t + forward_excitation_.A()));
  return std::pow(base, 0.5 * parameters.Alpha0());
}

} // namespace gra
