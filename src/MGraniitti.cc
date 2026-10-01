// GRANIITTI Monte Carlo main class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <compare>
#include <complex>
#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <initializer_list>
#include <iostream>
#include <iterator>
#include <limits>
#include <memory>
#include <mutex>
#include <random>
#include <regex>
#include <set>
#include <stdexcept>
#include <string>
#include <string_view>
#include <thread>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Process/MProcessTable.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MSteering.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/Sampling/MNeuroJac.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MLHE.h"
#include "Graniitti/Tech/MTimer.h"

// HepMC3
#include "HepMC3/Attribute.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/WriterAsciiHepMC2.h"
#include "HepMC3/WriterHEPEVT.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

namespace gra {

using json = nlohmann::json;

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::PI;
using gra::math::pow2;
using gra::math::pow3;
using gra::math::zi;

namespace {

// Resolve an output stem and create its parent directory
std::string OutputPath(const std::string &output, const std::string &folder, const std::string &extension) {
  std::filesystem::path path(output);
  if (!path.has_parent_path()) { path = std::filesystem::path(aux::GetBasePath(2)) / folder / path; }
  path += "." + extension;
  aux::CreateDirectory(path.parent_path().string());
  return path.string();
}

// Validate steering against the capabilities of the selected process
void ValidateProcessCommands(const std::vector<aux::OneCMD> &commands, const ProcessInfo &process) {
  const std::array<std::string_view, 12> allowed = {
      "R", "RES", "PDG", "j", "FLATAMP", "FLATMASS2", "OFFSHELL", "SPINGEN", "SPINDEC", "QMETRICS", "MP_FRAME", "MMAX"};
  for (const auto &command : commands) {
    if (std::find(allowed.begin(), allowed.end(), command.id) == allowed.end()) {
      throw std::invalid_argument("MGraniitti::ReadProcessParam: unknown process command @" + command.id);
    }
    if (command.id == "MP_FRAME" && process.model != ReggeProductionModel::MP) {
      throw std::invalid_argument("MGraniitti::ReadProcessParam: @MP_FRAME applies only to MP processes");
    }
    if (command.id == "MMAX" && process.model != ReggeProductionModel::GP) {
      throw std::invalid_argument("MGraniitti::ReadProcessParam: @MMAX applies only to GP processes");
    }
    if ((command.id == "R" || command.id == "RES") && process.resonance == ResonanceType::None) {
      throw std::invalid_argument("MGraniitti::ReadProcessParam: @" + command.id +
                                  " requires a process with free resonance input");
    }
  }
}

// Read one complete finite resonance override, with nonnegative masses and weights
double ReadResValue(const std::string &text, const std::string &key, bool nonnegative) {
  try {
    std::size_t pos = 0;
    const double value = std::stod(text, &pos);
    if (pos != text.size() || !std::isfinite(value) || (nonnegative && value < 0.0)) {
      throw std::invalid_argument("invalid numeric value");
    }
    return value;
  } catch (const std::exception &error) {
    throw std::invalid_argument("MGraniitti::ReadProcessParam: @R[] " + key + ": " + error.what());
  }
}

// VGRID serialization shared with the icetune graniitti driver
constexpr unsigned int kVGridSchemaVersion = 1;

struct FiducialPDGSelectorSpec {
  std::vector<int> pdg;
  std::vector<bool> pdg_abs;
};

// Infer the maximum-weight mode from one nullable steering value
MaxWeightMode ParseMaximumWeightValue(const json &value, double &max_w_value) {
  max_w_value = 0.0;
  if (value.is_null()) {
    return MaxWeightMode::Estimate;
  }
  if (!value.is_number()) {
    throw std::invalid_argument(
        "INTEGRATOR::max_w_value must be null or a finite positive number");
  }
  max_w_value = value.get<double>();
  if (!(max_w_value > 0.0) || !std::isfinite(max_w_value)) {
    throw std::invalid_argument(
        "INTEGRATOR::max_w_value must be null or a finite positive number");
  }
  return MaxWeightMode::Fixed;
}

// Round a positive maximum upward at the displayed scientific precision
long double ConservativeMaximumForPrint(double value) {
  if (!(value > 0.0) || !std::isfinite(value)) {
    return static_cast<long double>(value);
  }
  constexpr int DECIMAL_PLACES = 3;
  const int exponent = static_cast<int>(std::floor(std::log10(value)));
  long double quantum =
      std::pow(10.0L, static_cast<long double>(exponent - DECIMAL_PLACES));
  if (!(quantum > 0.0L)) {
    quantum =
        static_cast<long double>(std::numeric_limits<double>::denorm_min());
  }
  return std::ceil(static_cast<long double>(value) / quantum) * quantum;
}

// Parse one generation overflow action
OverflowAction ParseOverflowAction(const std::string &action) {
  if (action == "keep_as_weighted") {
    return OverflowAction::KeepAsWeighted;
  }
  if (action == "break") {
    return OverflowAction::Break;
  }
  throw std::invalid_argument(
      "NUMERICS_MC::MAX_WEIGHT::overflow_action must be keep_as_weighted, "
      "or break");
}

// Compute the stable generation overflow action label
std::string OverflowActionName(OverflowAction action) {
  if (action == OverflowAction::KeepAsWeighted) {
    return "keep_as_weighted";
  }
  return "break";
}

// Parse a disabled boolean or one signed 64-bit custom-cut identifier
std::int64_t ParseUserCutID(const json &value) {
  if (value.is_boolean()) {
    if (!value.get<bool>()) {
      return 0;
    }
    throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::USERCUTS "
                                "must be false or a signed 64-bit integer");
  }
  if (value.is_number_unsigned()) {
    const std::uint64_t id = value.get<std::uint64_t>();
    if (id >
        static_cast<std::uint64_t>(std::numeric_limits<std::int64_t>::max())) {
      throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::USERCUTS "
                                  "exceeds the signed 64-bit range");
    }
    return static_cast<std::int64_t>(id);
  }
  if (!value.is_number_integer()) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::USERCUTS "
                                "must be false or a signed 64-bit integer");
  }
  return value.get<std::int64_t>();
}

// Compute a copy of one string without leading or trailing whitespace
std::string TrimSelectorString(const std::string &value) {
  const std::size_t first = value.find_first_not_of(" \t\n\r");
  if (first == std::string::npos) {
    return "";
  }
  const std::size_t last = value.find_last_not_of(" \t\n\r");
  return value.substr(first, last - first + 1);
}

// Split a bracketed PDG selector body into comma separated item tokens
std::vector<std::string> SplitPDGSelectorItems(const std::string &body,
                                               const std::string &key) {
  std::vector<std::string> items;
  std::size_t start = 0;
  while (start <= body.size()) {
    const std::size_t comma = body.find(',', start);
    const std::size_t stop = (comma == std::string::npos) ? body.size() : comma;
    const std::string item =
        TrimSelectorString(body.substr(start, stop - start));
    if (item.empty()) {
      throw std::invalid_argument(
          "MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL key '" + key +
          "' has an empty PDG selector entry");
    }
    items.push_back(item);
    if (comma == std::string::npos) {
      break;
    }
    start = comma + 1;
  }
  return items;
}

// Parse one signed integer PDG selector token without accepting floats
int ParsePDGIntegerToken(const std::string &token, const std::string &key) {
  try {
    std::size_t pos = 0;
    const long value = std::stol(token, &pos);
    if (pos != token.size()) {
      throw std::invalid_argument("non-integer suffix");
    }
    if (value < std::numeric_limits<int>::min() ||
        value > std::numeric_limits<int>::max()) {
      throw std::out_of_range("PDG id outside int range");
    }
    return static_cast<int>(value);
  } catch (const std::exception &e) {
    throw std::invalid_argument(
        "MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL key '" + key +
        "' contains invalid integer PDG selector '" + token + "': " + e.what());
  }
}

// Parse one PDG selector item, optionally wrapped as ABS(pdg)
void ParsePDGSelectorItem(const std::string &item, const std::string &key,
                          FiducialPDGSelectorSpec &selector) {
  const std::string abs_prefix = "ABS(";
  bool abs_match = false;
  std::string token = item;

  if (item.rfind(abs_prefix, 0) == 0) {
    if (item.size() <= abs_prefix.size() + 1 || item.back() != ')') {
      throw std::invalid_argument(
          "MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL key '" + key +
          "' has malformed ABS selector '" + item + "'");
    }
    abs_match = true;
    token = TrimSelectorString(
        item.substr(abs_prefix.size(), item.size() - abs_prefix.size() - 1));
  }

  int pdg = token == "j" ? PDG::PDG_hard_jet : ParsePDGIntegerToken(token, key);
  if (abs_match) {
    if (pdg == std::numeric_limits<int>::min()) {
      throw std::invalid_argument(
          "MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL key '" + key +
          "' has ABS selector outside the positive int range");
    }
    pdg = std::abs(pdg);
  }
  selector.pdg.push_back(pdg);
  selector.pdg_abs.push_back(abs_match);
}

// Parse one FIDCUTS::CENTRAL PDG selector key such as "22", "[13,-13]" or
// "ABS(13)"
FiducialPDGSelectorSpec ParseFiducialPDGSelector(const std::string &key) {
  const std::string selector_key = TrimSelectorString(key);
  if (selector_key.empty()) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL has "
                                "an empty PDG selector key");
  }

  std::vector<std::string> items;
  if (selector_key.front() == '[' || selector_key.back() == ']') {
    if (selector_key.size() < 2 || selector_key.front() != '[' ||
        selector_key.back() != ']') {
      throw std::invalid_argument(
          "MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL key '" + key +
          "' has malformed bracket PDG selector");
    }
    items = SplitPDGSelectorItems(
        selector_key.substr(1, selector_key.size() - 2), key);
  } else {
    items = {selector_key};
  }

  FiducialPDGSelectorSpec out;
  out.pdg.reserve(items.size());
  out.pdg_abs.reserve(items.size());
  for (const auto &item : items) {
    ParsePDGSelectorItem(item, key, out);
  }
  return out;
}

// Validate one finite fiducial range against optional physical bounds
void AssertFiducialRangeValues(const std::vector<double> &values,
                               const std::string &path,
                               const std::vector<double> &bounds = {}) {
  if (values.size() != 2) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: " + path +
                                " range must contain two values");
  }
  if (!std::isfinite(values[0]) || !std::isfinite(values[1])) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: " + path +
                                " range values must be finite");
  }
  if (bounds.empty()) {
    gra::aux::AssertCut(values, path, true);
  } else {
    gra::aux::AssertCutRange(values, bounds, path, true);
  }
}

// Reject unsupported keys in one fiducial steering block
void AssertFiducialBlockKeys(const json &block, const std::string &path,
                             std::initializer_list<const char *> allowed) {
  for (auto it = block.begin(); it != block.end(); ++it) {
    const bool supported =
        std::any_of(allowed.begin(), allowed.end(),
                    [&](const char *key) { return it.key() == key; });
    if (!supported) {
      throw std::invalid_argument("MGraniitti::ReadFidCuts: " + path +
                                  " has unsupported key '" + it.key() + "'");
    }
  }
}

// Read one optional observable range from a PDG fiducial cut block
void ReadOptionalPDGFiducialRange(const json &block, const std::string &pdg_key,
                                  const std::string &name,
                                  gra::FIDCUTRANGE &range) {
  if (!block.contains(name)) {
    return;
  }

  const std::vector<double> values = block.at(name).get<std::vector<double>>();
  const std::string path = "FIDCUTS::CENTRAL::" + pdg_key + "::" + name;
  if (name == "M" || name == "Pt" || name == "Et") {
    AssertFiducialRangeValues(values, path,
                              {0.0, std::numeric_limits<double>::max()});
  } else {
    AssertFiducialRangeValues(values, path);
  }
  range.active = true;
  range.min = values[0];
  range.max = values[1];
}

// Compute one optional FIDCUTS block after validating its shape
const json *OptionalFiducialBlock(const json &fidcuts_json,
                                  const std::string &block_name) {
  if (!fidcuts_json.contains(block_name)) {
    return nullptr;
  }
  const json &block = fidcuts_json.at(block_name);
  if (!block.is_object()) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::" +
                                block_name + " must be an object");
  }
  return &block;
}

// Read one optional fiducial observable range and enforce any physical bounds
void ReadOptionalFiducialRange(const json *block, const std::string &path,
                               const std::string &name, bool &active,
                               double &min, double &max,
                               const std::vector<double> &bounds = {}) {
  if (block == nullptr || !block->contains(name)) {
    return;
  }

  const std::vector<double> values = block->at(name).get<std::vector<double>>();
  AssertFiducialRangeValues(values, path + "::" + name, bounds);
  active = true;
  min = values[0];
  max = values[1];
}

// Read one PDG selector entry from FIDCUTS::CENTRAL
void ReadCentralPDGFiducialCut(const std::string &pdg_key, const json &block,
                               gra::FIDCUT &fcuts) {
  if (!block.is_object()) {
    throw std::invalid_argument("MGraniitti::ReadFidCuts: FIDCUTS::CENTRAL::" +
                                pdg_key + " must be an object");
  }
  AssertFiducialBlockKeys(block, "FIDCUTS::CENTRAL::" + pdg_key,
                          {"M", "Rap", "Eta", "Pt", "Et"});

  gra::FIDPDGCUT cut;
  const FiducialPDGSelectorSpec selector = ParseFiducialPDGSelector(pdg_key);
  cut.pdg = selector.pdg;
  cut.pdg_abs = selector.pdg_abs;
  ReadOptionalPDGFiducialRange(block, pdg_key, "M", cut.M);
  ReadOptionalPDGFiducialRange(block, pdg_key, "Rap", cut.Rap);
  ReadOptionalPDGFiducialRange(block, pdg_key, "Eta", cut.Eta);
  ReadOptionalPDGFiducialRange(block, pdg_key, "Pt", cut.Pt);
  ReadOptionalPDGFiducialRange(block, pdg_key, "Et", cut.Et);

  if (cut.HasActiveRange()) {
    fcuts.pdg_cuts.push_back(cut);
  }
}

// Read optional FIDCUTS::CENTRAL final-state, system and PDG-selected cuts
void ReadCentralFiducialCuts(const json &fidcuts_json, gra::FIDCUT &fcuts) {
  const json *central = OptionalFiducialBlock(fidcuts_json, "CENTRAL");
  if (central == nullptr) {
    return;
  }

  const json *particle = OptionalFiducialBlock(*central, "*");
  if (particle != nullptr) {
    AssertFiducialBlockKeys(*particle, "FIDCUTS::CENTRAL::*",
                            {"Eta", "Rap", "Pt", "Et"});
  }
  ReadOptionalFiducialRange(particle, "FIDCUTS::CENTRAL::*", "Eta",
                            fcuts.particle_eta_active, fcuts.eta_min,
                            fcuts.eta_max);
  ReadOptionalFiducialRange(particle, "FIDCUTS::CENTRAL::*", "Rap",
                            fcuts.particle_rap_active, fcuts.rap_min,
                            fcuts.rap_max);
  ReadOptionalFiducialRange(
      particle, "FIDCUTS::CENTRAL::*", "Pt", fcuts.particle_pt_active,
      fcuts.pt_min, fcuts.pt_max, {0.0, std::numeric_limits<double>::max()});
  ReadOptionalFiducialRange(
      particle, "FIDCUTS::CENTRAL::*", "Et", fcuts.particle_Et_active,
      fcuts.Et_min, fcuts.Et_max, {0.0, std::numeric_limits<double>::max()});

  const json *system = OptionalFiducialBlock(*central, "SYSTEM");
  if (system != nullptr) {
    AssertFiducialBlockKeys(*system, "FIDCUTS::CENTRAL::SYSTEM",
                            {"M", "Rap", "Pt"});
  }
  ReadOptionalFiducialRange(system, "FIDCUTS::CENTRAL::SYSTEM", "M",
                            fcuts.system_M_active, fcuts.M_min, fcuts.M_max,
                            {0.0, std::numeric_limits<double>::max()});
  ReadOptionalFiducialRange(system, "FIDCUTS::CENTRAL::SYSTEM", "Rap",
                            fcuts.system_Rap_active, fcuts.Y_min, fcuts.Y_max);
  ReadOptionalFiducialRange(system, "FIDCUTS::CENTRAL::SYSTEM", "Pt",
                            fcuts.system_Pt_active, fcuts.Pt_min, fcuts.Pt_max,
                            {0.0, std::numeric_limits<double>::max()});

  for (auto it = central->begin(); it != central->end(); ++it) {
    if (it.key() == "*" || it.key() == "SYSTEM") {
      continue;
    }
    ReadCentralPDGFiducialCut(it.key(), it.value(), fcuts);
  }
}

// Read system-level GENCUTS for one phase-space class key
void ReadSystemGenCutsBlock(const json &j, const std::string &class_key,
                            const std::string &class_name, gra::GENCUT &gcuts) {
  const std::string XID = "GENCUTS";

  const json *block = nullptr;
  try {
    block = &j.at(XID).at(class_key);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MGraniitti::ReadGenCuts: " + class_key + " " +
                                class_name +
                                " phase space class requires from user: "
                                "\"GENCUTS\" : { \"" +
                                class_key +
                                "\" : { \"Rap\" : [min, max], \"M\" : [min, "
                                "max] }}; underlying error: " +
                                error.what());
  }
  if (!block->is_object()) {
    throw std::invalid_argument("MGraniitti::ReadGenCuts: GENCUTS::" +
                                class_key + " must be an object");
  }

  std::vector<double> Y;
  try {
    Y = block->at("Rap").get<std::vector<double>>();
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MGraniitti::ReadGenCuts: " + class_key + " " + class_name +
        " phase space class requires from user: "
        "\"GENCUTS\" : { \"" +
        class_key +
        "\" : { \"Rap\" : [min, max] }}; underlying error: " + error.what());
  }
  gra::aux::AssertCut(Y, "GENCUTS::" + class_key + "::Rap", true);
  gcuts.Y_min = Y[0];
  gcuts.Y_max = Y[1];

  std::vector<double> M;
  try {
    M = block->at("M").get<std::vector<double>>();
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MGraniitti::ReadGenCuts: " + class_key + " " + class_name +
        " phase space class requires from user: "
        "\"GENCUTS\" : { \"" +
        class_key +
        "\" : { \"M\" : [min, max] }}; underlying error: " + error.what());
  }
  gra::aux::AssertCutRange(M, {0.0, 1e32}, "GENCUTS::" + class_key + "::M",
                           true);
  gcuts.M_min = M[0];
  gcuts.M_max = M[1];

  if (block->contains("Pt")) {
    const std::vector<double> pt = block->at("Pt").get<std::vector<double>>();
    gra::aux::AssertCutRange(pt, {0.0, 1e32}, "GENCUTS::" + class_key + "::Pt",
                             true);
    gcuts.forward_pt_min = pt[0];
    gcuts.forward_pt_max = pt[1];
  }

  if (block->contains("Xi")) {
    const std::vector<double> Xi = block->at("Xi").get<std::vector<double>>();
    gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::" + class_key + "::Xi",
                             true);
    gcuts.XI_min = Xi[0];
    gcuts.XI_max = Xi[1];
  }
}

// Compute true when two process descriptors have identical help metadata
bool SameProcessDescriptor(const ProcessDescriptor &first,
                           const ProcessDescriptor &second) {
  return first.process == second.process && first.model == second.model &&
         first.channels == second.channels &&
         first.comments == second.comments &&
         first.display_order == second.display_order;
}

// Replace one phase-space token in a process command or full final-state syntax
std::string ReplaceProcessMode(const std::string &command,
                               const std::string &current_mode,
                               const std::string &new_mode) {
  const std::string token = "<" + current_mode + ">";
  const std::size_t position = command.find(token);
  if (position == std::string::npos) {
    throw std::logic_error("Missing phase-space token " + token +
                           " in process command " + command);
  }
  const std::size_t token_end = position + token.size();
  if (token_end != command.size() &&
      command.compare(token_end, 4, " -> ") != 0) {
    throw std::logic_error(
        "Phase-space token is not terminal in process command " + command);
  }
  return command.substr(0, position) + "<" + new_mode + ">" +
         command.substr(token_end);
}

// Compute whether a C selector needs the factorized luminosity phase space
bool FactorizedCentralAlias(const std::string &process) {
  constexpr std::string_view suffix = "[FLUX]<C>";
  return process.size() >= suffix.size() &&
         process.compare(process.size() - suffix.size(), suffix.size(),
                         suffix) == 0;
}

// Build one shared process table after checking the F and C registries agree
std::vector<std::vector<std::string>>
SharedCentralProcessRows(const MSubProc &factorized,
                         const MSubProc &continuum, const MProcessTable &table) {
  if (factorized.ProcessRegistry.size() != continuum.ProcessRegistry.size()) {
    throw std::logic_error(
        "The <F> and <C> process registries have different sizes");
  }

  for (const auto &entry : factorized.ProcessRegistry) {
    const std::string continuum_command =
        ReplaceProcessMode(entry.first, "F", "C");
    const auto match = continuum.ProcessRegistry.find(continuum_command);
    if (match == continuum.ProcessRegistry.end() ||
        !SameProcessDescriptor(entry.second, match->second)) {
      throw std::logic_error("The <F> and <C> process registries differ at " +
                             entry.first);
    }
  }

  auto rows = table.Rows(factorized);
  for (auto &row : rows) {
    if (!row.empty() && !row[0].empty()) { row[0] = ReplaceProcessMode(row[0], "F", "F|C"); }
  }
  return rows;
}

// Compute process commands in the same order as the process help table
std::vector<std::string> OrderedProcessCommands(const MSubProc &subprocess) {
  std::vector<std::pair<std::string, ProcessDescriptor>> entries(
      subprocess.ProcessRegistry.begin(), subprocess.ProcessRegistry.end());
  std::sort(entries.begin(), entries.end(),
            [](const auto &first, const auto &second) {
              if (first.second.display_order != second.second.display_order) {
                return first.second.display_order < second.second.display_order;
              }
              return first.first < second.first;
            });

  std::vector<std::string> commands;
  commands.reserve(entries.size());
  for (const auto &entry : entries) {
    commands.push_back(entry.first);
  }
  return commands;
}

// Replace exact phase-space modes in generated final-state table rows
std::vector<std::vector<std::string>>
ReplaceFinalStateRowMode(std::vector<std::vector<std::string>> rows,
                         const std::string &current_mode,
                         const std::string &new_mode) {
  for (auto &row : rows) {
    if (row.size() != 4) {
      throw std::logic_error(
          "Generated final-state process rows must have four columns");
    }
    row[0] = ReplaceProcessMode(row[0], current_mode, new_mode);
  }
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

} // namespace

// Constructor
MGraniitti::MGraniitti() {
  // Print general layout
  PrintInit();
}

// Destructor
MGraniitti::~MGraniitti() {
  // Destroy processes
  for (const auto &i : indices(pvec)) {
    delete pvec[i];
  }

  std::cout << "~MGraniitti [DONE]" << std::endl;
}

void MGraniitti::Initialize() {
  proc->PrepareRun();

  // Integrate
  CallIntegrator(0);
}

// Initialize with external Eikonal
void MGraniitti::Initialize(const MEikonal &eikonal_in) {
  proc->PrepareRun(eikonal_in);

  // Integrate
  CallIntegrator(0);
}

// Save the configured integration proposal and statistics to a file
void MGraniitti::SaveVGRID() const {
  // MC stats
  json jstat;
  stat.struct2json(jstat);

  // Integrator-specific proposal
  json j;
  j["INTEGRATOR"] = INTEGRATOR;
  if (INTEGRATOR == "VEGAS") {
    j["VEGAS"] = vegas.Serialize();
  } else if (INTEGRATOR == "NEUROJAC") {
    j["NEUROJAC"] = json::parse(neurojac.SerializeModel());
  } else {
    throw std::logic_error("MGraniitti::SaveVGRID: integrator does not support "
                           "serialized proposals");
  }

  // Shared process metadata
  j["SCHEMA_VERSION"] = kVGridSchemaVersion;
  j["PROPOSAL_ONLY"] = ProposalOnlyState();
  j["BEAM_ENERGY"] = {nuclear::SteeringBeamEnergy(proc->state.lts.beam1.pdg,
                                                  proc->state.lts.pbeam1.E()),
                      nuclear::SteeringBeamEnergy(proc->state.lts.beam2.pdg,
                                                  proc->state.lts.pbeam2.E())};
  j["PROCESS"] = FULL_PROCESS_STR;
  j["STAT"] = jstat;
  if (proc->state.lts.process.QMETRICS) { j["QMETRICS"] = QMetrics().Serialize(); }
  j["INTEGRATION"] = {{"min_samples", integration_param.MIN_SAMPLES},
                      {"max_samples", integration_param.MAX_SAMPLES},
                      {"precision", integration_param.PRECISION}};
  j["MAX_WEIGHT"] = {
      {"max_w_value", max_weight_param.MODE == MaxWeightMode::Fixed
                          ? json(max_weight_param.MAX_W_VALUE)
                          : json(nullptr)},
      {"safety_factor", max_weight_param.SAFETY_FACTOR},
      {"max_overflow_prob", max_weight_param.MAX_OVERFLOW_PROB},
      {"confidence_level", max_weight_param.CONFIDENCE_LEVEL},
      {"overflow_action", OverflowActionName(max_weight_param.OVERFLOW_ACTION)},
      {"envelope", max_weight_state.envelope},
      {"validation_samples", max_weight_state.validation_samples},
      {"envelope_updates", max_weight_state.envelope_updates},
      {"initialized", max_weight_state.initialized}};
  j["LOOPSCREEN"] = proc->GetScreening();
  const double a1 =
      static_cast<double>(nuclear::BeamMassNumber(proc->state.lts.beam1.pdg));
  const double a2 =
      static_cast<double>(nuclear::BeamMassNumber(proc->state.lts.beam2.pdg));
  j["SQRTS"] = (proc->state.lts.pbeam1 / a1 + proc->state.lts.pbeam2 / a2).M();

  const std::string filename = OutputPath(OUTPUT, "vgrid", "vgrid");
  std::ofstream file(filename);
  if (!file.is_open()) {
    throw std::runtime_error("MGraniitti::SaveVGRID: cannot open " + filename);
  }
  file << j << std::endl;
  if (!file.good()) {
    throw std::runtime_error("MGraniitti::SaveVGRID: failed to write " +
                             filename);
  }

  std::cout << "MGraniitti::SaveVGRID: Saved MC integrator state to " + filename
            << std::endl;
}

// Read one serialized integration proposal and its statistics
void MGraniitti::ReadVGRID(const std::string &inputfile) {
  std::lock_guard<std::mutex> lifecycle_lock(state_mutex);
  if (max_weight_state.generation_active) {
    throw std::logic_error(
        "MGraniitti::ReadVGRID: cannot replace integration state during "
        "event generation");
  }
  std::cout << "MGraniitti::ReadVGRID: Reading pre-computed MC integration "
               "state from " +
                   inputfile
            << std::endl;
  std::cout << "** Re-cycling this proposal requires that the process energy, "
               "parameters and cuts "
               "match! **"
            << std::endl;

  // Read and parse
  const std::string data = gra::aux::GetInputData(inputfile);
  json j;
  try {
    j = json::parse(data);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MGraniitti::ReadVGRID: Error parsing " +
                                inputfile + ": " + error.what());
  }

  if (!j.contains("SCHEMA_VERSION") ||
      j.at("SCHEMA_VERSION").get<unsigned int>() != kVGridSchemaVersion) {
    throw std::invalid_argument("MGraniitti::ReadVGRID: incompatible VGRID "
                                "schema, rebuild the integration state");
  }

  const json &beam_energy = j.at("BEAM_ENERGY");
  if (!beam_energy.is_array() || beam_energy.size() != 2) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: invalid BEAM_ENERGY metadata");
  }
  const std::array<double, 2> current_beam_energy = {
      nuclear::SteeringBeamEnergy(proc->state.lts.beam1.pdg,
                                  proc->state.lts.pbeam1.E()),
      nuclear::SteeringBeamEnergy(proc->state.lts.beam2.pdg,
                                  proc->state.lts.pbeam2.E())};
  for (const auto &i : indices(current_beam_energy)) {
    if (!aux::AssertRatio(beam_energy.at(i).get<double>(),
                          current_beam_energy[i], 1e-12)) {
      throw std::invalid_argument("MGraniitti::ReadVGRID: beam energy does not "
                                  "match the current process");
    }
  }

  if (j.contains("INTEGRATOR") &&
      j.at("INTEGRATOR").get<std::string>() != INTEGRATOR) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: serialized integrator = " +
        j.at("INTEGRATOR").get<std::string>() +
        " but current integrator = " + INTEGRATOR);
  }

  // Restore the selected proposal
  if (INTEGRATOR == "VEGAS") {
    vegas.Deserialize(j.at("VEGAS"), proc->GetdLIPSDim(), vparam);
  } else if (INTEGRATOR == "NEUROJAC") {
    neurojac.DeserializeModel(j.at("NEUROJAC").dump(), proc->GetdLIPSDim());
  } else {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: current integrator has no serialized proposal");
  }

  const bool proposal_only = j.at("PROPOSAL_ONLY").get<bool>();

  // Require matching screening for state containing physical integration data
  const bool loopscreen = j.at("LOOPSCREEN");
  if (loopscreen != proc->GetScreening()) {
    if (!proposal_only) {
      throw std::invalid_argument(
          "MGraniitti::ReadVGRID: Input .vgrid with LOOPSCREEN = " +
          aux::bool_cast(loopscreen) + " but process with LOOPSCREEN = " +
          aux::bool_cast(proc->GetScreening()));
    }
    std::cout << rang::fg::red << "WARNING: proposal-only VGRID LOOPSCREEN = "
              << aux::bool_cast(loopscreen)
              << " differs from the current process LOOPSCREEN = "
              << aux::bool_cast(proc->GetScreening())
              << ". The proposal remains valid, but sampling efficiency may "
                 "not be optimal"
              << rang::fg::reset << std::endl;
  }

  // Check that the nucleon-pair CMS energy matches
  const double sqrts = j.at("SQRTS");
  const double a1 =
      static_cast<double>(nuclear::BeamMassNumber(proc->state.lts.beam1.pdg));
  const double a2 =
      static_cast<double>(nuclear::BeamMassNumber(proc->state.lts.beam2.pdg));
  const double current_sqrts =
      (proc->state.lts.pbeam1 / a1 + proc->state.lts.pbeam2 / a2).M();
  if (!aux::AssertRatio(sqrts, current_sqrts, 0.01)) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: Input .vgrid with SQRTS = " +
        std::to_string(sqrts) +
        " but process with SQRTS = " + std::to_string(current_sqrts));
  }

  // Read in MC stats
  stat.json2struct(j["STAT"]);
  if (proc->state.lts.process.QMETRICS) {
    for (auto *worker : pvec) { worker->state.lts.qmetrics.Reset(); }
  }
  if (proc->state.lts.process.QMETRICS && j.contains("QMETRICS")) {
    proc->state.lts.qmetrics.Deserialize(j.at("QMETRICS"));
  }
  if (proposal_only != ProposalOnlyState()) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: PROPOSAL_ONLY metadata does not match the "
        "serialized integration state");
  }

  // Require identical frozen-proposal integration settings
  const json &stored_integration = j.at("INTEGRATION");
  if (stored_integration.at("min_samples").get<std::uint64_t>() !=
          integration_param.MIN_SAMPLES ||
      stored_integration.at("max_samples").get<std::uint64_t>() !=
          integration_param.MAX_SAMPLES ||
      !std::is_eq(stored_integration.at("precision").get<double>() <=>
                  integration_param.PRECISION)) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: integration policy does not match the "
        "steering card");
  }

  // Restore maximum-weight calibration only under an identical policy
  const json &stored_maximum = j.at("MAX_WEIGHT");
  double stored_max_w_value = 0.0;
  const MaxWeightMode stored_mode = ParseMaximumWeightValue(
      stored_maximum.at("max_w_value"), stored_max_w_value);
  if (stored_mode != max_weight_param.MODE ||
      !std::is_eq(stored_maximum.at("safety_factor").get<double>() <=>
                  max_weight_param.SAFETY_FACTOR) ||
      !std::is_eq(stored_maximum.at("max_overflow_prob").get<double>() <=>
                  max_weight_param.MAX_OVERFLOW_PROB) ||
      !std::is_eq(stored_maximum.at("confidence_level").get<double>() <=>
                  max_weight_param.CONFIDENCE_LEVEL) ||
      stored_maximum.at("overflow_action").get<std::string>() !=
          OverflowActionName(max_weight_param.OVERFLOW_ACTION) ||
      (max_weight_param.MODE == MaxWeightMode::Fixed &&
       !std::is_eq(stored_max_w_value <=> max_weight_param.MAX_W_VALUE))) {
    throw std::invalid_argument(
        "MGraniitti::ReadVGRID: maximum-weight policy does not match the "
        "current configuration");
  }
  stored_maximum.at("envelope").get_to(max_weight_state.envelope);
  stored_maximum.at("validation_samples")
      .get_to(max_weight_state.validation_samples);
  stored_maximum.at("envelope_updates")
      .get_to(max_weight_state.envelope_updates);
  stored_maximum.at("initialized").get_to(max_weight_state.initialized);
}

// Print fast histograms
void MGraniitti::PrintHistograms() {
  if (proc->HIST != 0) {
    HistogramFusion();

    // First print to screen (flushes buffers)
    proc->PrintHistograms();

    // Then save
    const std::string filename = OutputPath(OUTPUT, "output", "hfast");
    proc->SaveHistograms(filename);
  }
}

// Unify histogram bounds across threads, so the histogram fusion is possible
void MGraniitti::UnifyHistogramBounds() {
  // Loop over all 1D histograms
  for (auto const &xpoint : proc->h1) {
    int xbins = 0;
    double xmin = 0;
    double xmax = 0;

    // Fuse buffer values [start index 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      proc->h1[xpoint.first].FuseBuffer(pvec[p]->h1[xpoint.first]);
    }

    proc->h1[xpoint.first].FlushBuffer();
    proc->h1[xpoint.first].GetBounds(xbins, xmin, xmax);
    proc->h1[xpoint.first].ResetBounds(xbins, xmin, xmax); // New start

    // Loop over processes and set histogram bounds [start index = 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      pvec[p]->h1[xpoint.first].ResetBounds(xbins, xmin, xmax);
    }
  }
  // Loop over all 2D histograms
  for (auto const &xpoint : proc->h2) {
    int xbins = 0;
    double xmin = 0;
    double xmax = 0;

    int ybins = 0;
    double ymin = 0;
    double ymax = 0;

    // Fuse buffer values [start index 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      proc->h2[xpoint.first].FuseBuffer(pvec[p]->h2[xpoint.first]);
    }

    proc->h2[xpoint.first].FlushBuffer();
    proc->h2[xpoint.first].GetBounds(xbins, xmin, xmax, ybins, ymin, ymax);
    proc->h2[xpoint.first].ResetBounds(xbins, xmin, xmax, ybins, ymin,
                                       ymax); // New start

    // Loop over processes and set histogram bounds [start index = 1]
    for (std::size_t p = 1; p < pvec.size(); ++p) {
      pvec[p]->h2[xpoint.first].ResetBounds(xbins, xmin, xmax, ybins, ymin,
                                            ymax);
    }
  }
}

// Fuse histograms for N-fold statistics
void MGraniitti::HistogramFusion() {
  if (hist_fusion_done == false) {
    // START with process index 1, because 0 is the base
    for (std::size_t i = 1; i < pvec.size(); ++i) {
      // Loop over all 1D histograms
      for (auto const &xpoint : proc->h1) {
        proc->h1[xpoint.first] =
            proc->h1[xpoint.first] + pvec[i]->h1[xpoint.first];
      }
      // Loop over all 2D histograms
      for (auto const &xpoint : proc->h2) {
        proc->h2[xpoint.first] =
            proc->h2[xpoint.first] + pvec[i]->h2[xpoint.first];
      }
    }
    hist_fusion_done = true;
  } else {
    std::cout
        << "MGraniitti::HistogramFusion: Multithreaded histograms have been "
           "fused already"
        << std::endl;
  }
}

// Generate events
//
void MGraniitti::Generate() {
  if (NEVENTS > 0) { // Generate events
    {
      std::lock_guard<std::mutex> lock(state_mutex);
      max_weight_state.ResetGeneration();
      stat.ResetGeneration();
      weighted_event_stats.Reset();
    }

    try {
      // Validate the frozen envelope before creating any output
      static_cast<void>(GetGenerationMaxweight());

      // Open output only after successful integration initialization
      InitFileOutput();
      CallIntegrator(NEVENTS);

      bool overflow_break = false;
      {
        std::lock_guard<std::mutex> lock(state_mutex);
        overflow_break =
            !WEIGHTED && max_weight_state.generation_overflow &&
            max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted;
      }
      if (overflow_break) {
        CloseFileOutput();
        PreserveFailedOutput();
        throw std::runtime_error(
            "MAX_WEIGHT: overflow_action=break stopped unweighted generation "
            "after an overweight event");
      }

      if (owns_file_output) { CloseFileOutput(); }
      std::lock_guard<std::mutex> lock(state_mutex);
      max_weight_state.EndGeneration();
      return;
    } catch (...) {
      bool preserve_lhe_overflow = false;
      {
        std::lock_guard<std::mutex> lock(state_mutex);
        preserve_lhe_overflow =
            FORMAT == "lhe" && !WEIGHTED &&
            max_weight_state.generation_overflow;
        max_weight_state.EndGeneration();
      }
      if (preserve_lhe_overflow && owns_file_output) {
        CloseFileOutput();
        PreserveFailedOutput();
      }
      throw;
    }
  }
}

// Compute concrete native and generated final states for all valid phase-space process classes
//
std::vector<std::vector<std::string>> MGraniitti::GetFinalStateRows() const {
  const auto factorized_rows = proc_F.ProcPtr.SupportedFinalStateRows();
  const auto continuum_rows = proc_C.ProcPtr.SupportedFinalStateRows();
  const auto factorized_comparison =
      ReplaceFinalStateRowMode(factorized_rows, "F", "MODE");
  const auto continuum_comparison =
      ReplaceFinalStateRowMode(continuum_rows, "C", "MODE");
  if (factorized_comparison != continuum_comparison) {
    throw std::logic_error(
        "The <F> and <C> final-state registries differ");
  }

  std::vector<std::vector<std::string>> rows =
      ReplaceFinalStateRowMode(factorized_rows, "F", "F|C");
  const auto append = [&rows](const MSubProc &subprocess) {
    const auto additional = subprocess.SupportedFinalStateRows();
    rows.insert(rows.end(), additional.begin(), additional.end());
  };
  append(proc_Q.ProcPtr);
  append(proc_P.ProcPtr);
  auto hard_rows = proc_D.ProcPtr.SupportedFinalStateRows();
  for (auto &row : hard_rows) {
    row[0] = ReplaceProcessMode(row[0], row[0].find("<C>") != std::string::npos ? "C" : "F", "F|C");
  }
  rows.insert(rows.end(), hard_rows.begin(), hard_rows.end());

  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

// Compute process numbers available
//
std::vector<std::string> MGraniitti::GetProcessNumbers(const std::string &tune) const {
  const std::string model = tune.empty() ? modelparam : tune;
  const MProcessTable table(model);
  std::cout << rang::style::bold;
  std::cout << "Available processes: ";
  std::cout << rang::style::reset << std::endl << std::endl;
  std::cout << std::endl;

  const auto central_rows =
      SharedCentralProcessRows(proc_F.ProcPtr, proc_C.ProcPtr, table);
  std::cout << rang::style::bold
            << "<F|C> central production, with <F> factorized phase space, <C> direct phase space:"
            << rang::style::reset << std::endl;
  aux::PrintTable({"Command", "Beams", "Model", "Exchanges", "Supported final states", "Example <F> final state"},
                  central_rows);
  std::cout << std::endl;

  std::vector<std::string> list1 = OrderedProcessCommands(proc_F.ProcPtr);
  std::vector<std::string> list2 = OrderedProcessCommands(proc_C.ProcPtr);

  std::cout << rang::style::bold << "Quasielastic:" << rang::style::reset
            << std::endl;
  aux::PrintTable({"Command", "Beams", "Model", "Exchanges", "Supported final states", "Example final state"}, table.Rows(proc_Q.ProcPtr));
  std::vector<std::string> list3 = OrderedProcessCommands(proc_Q.ProcPtr);
  std::cout << std::endl;

  std::cout << rang::style::bold
            << "2->1 x (1->M) collinear:" << rang::style::reset << std::endl;
  aux::PrintTable({"Command", "Beams", "Model", "Exchanges", "Supported final states", "Example final state"}, table.Rows(proc_P.ProcPtr));
  std::vector<std::string> list4 = OrderedProcessCommands(proc_P.ProcPtr);
  std::cout << std::endl;

  std::cout << rang::style::bold << "<F|C> hard diffraction, with <F> stochastic central phase space, <C> adaptive central phase space:" << rang::style::reset
            << std::endl;
  aux::PrintTable({"Command", "Beams", "Model", "Exchanges", "Supported final states", "Example final state"}, table.Rows(proc_D.ProcPtr));
  std::vector<std::string> list5 = OrderedProcessCommands(proc_D.ProcPtr);
  std::cout << std::endl;

  auto final_state_rows = GetFinalStateRows();
  // Separate registered final states by command model
  for (std::size_t i = final_state_rows.size(); i > 1; --i) {
    const auto &previous = final_state_rows[i - 2][0];
    const auto &current  = final_state_rows[i - 1][0];
    if (previous.substr(0, previous.find('[')) != current.substr(0, current.find('['))) {
      final_state_rows.emplace(final_state_rows.begin() + i - 1);
    }
  }
  if (!final_state_rows.empty()) {
    std::cout << rang::style::bold
              << "Registered final states:" << rang::style::reset
              << std::endl;
    aux::PrintTable({"Command", "Beams", "Final state", "Amplitude syntax"},
                    final_state_rows);
    std::cout << std::endl;
  }

  std::cout << "card = exchanges from tune " << model << "\n"
            << "Examples use production and decay cards from tune " << model << "\n"
               "Regge resonances with lepton or nuclear beams require photoproduction\n"
               "&> uses isolated kinematic decays with <F> or <P>, apply the decay branching fraction separately\n"
               "Registered final states list explicit amplitude topologies\n"
               "MG2GRA regeneration updates registered states and examples\n\n";

  // Concatenate all
  list1.insert(list1.end(), list2.begin(), list2.end());
  list1.insert(list1.end(), list3.begin(), list3.end());
  list1.insert(list1.end(), list4.begin(), list4.end());
  list1.insert(list1.end(), list5.begin(), list5.end());

  return list1;
}

// (Re-)assign the pointers to the local memory space
void MGraniitti::InitProcessMemory(std::string process, unsigned int seed) {
  // These must be here!
  PROCESS = process;

  const auto tune = MModelTune::Load(gra::ResolveModelDataFile(modelparam, "GENERAL.json"));

  // <Q> process
  if (proc_Q.ProcPtr.ProcessExist(process)) {
    proc_Q = MQuasiElastic(process, syntax, tune);
    proc = &proc_Q;

    // <F> processes
  } else if (proc_F.ProcPtr.ProcessExist(process)) {
    proc_F = MFactorized(process, syntax, tune);
    proc = &proc_F;

    // <C> processes
  } else if (proc_C.ProcPtr.ProcessExist(process)) {
    if (FactorizedCentralAlias(process)) {
      proc_F = MFactorized(process, syntax, "C", tune);
      proc = &proc_F;
    } else {
      proc_C = MCentral(process, syntax, tune);
      proc = &proc_C;
    }

    // <P> processes
  } else if (proc_P.ProcPtr.ProcessExist(process)) {
    proc_P = MCollinear(process, syntax, tune);
    proc = &proc_P;

    // Hard diffraction with factorized or adaptive central production
  } else if (proc_D.ProcPtr.ProcessExist(process)) {
    proc_D = MHardDiffraction(process, syntax, tune);
    proc = &proc_D;
  } else {
    std::string str =
        "MGraniitti::InitProcessMemory: Unknown PROCESS: " + process;
    // GetProcessNumbers(); // This is done by main program
    throw std::invalid_argument(str);
  }

  // Set random seed last!
  proc->state.random.SetSeed(seed);
}

// Generate the requested number of distinct deterministic worker seeds
std::vector<uint32_t> MGraniitti::GenerateUniqueSeeds(uint32_t fixed_seed,
                                                      std::size_t num_seeds) {
  if (num_seeds == 0) {
    return {};
  }

  // Seed the random number generator with the fixed seed
  std::mt19937 generator(fixed_seed);

  // Use a uniform distribution to generate seeds within the range
  std::uniform_int_distribution<uint32_t> distribution(
      0, std::numeric_limits<uint32_t>::max());

  // Set to store unique seeds
  std::set<uint32_t> unique_seeds;

  // Insert the fixed seed as the first element
  std::vector<uint32_t> seeds;
  seeds.push_back(fixed_seed);
  unique_seeds.insert(fixed_seed);

  // Generate additional unique seeds
  while (seeds.size() < num_seeds) {
    uint32_t seed = distribution(generator);

    // Keep each random seed only once
    if (unique_seeds.insert(seed).second) {
      seeds.push_back(seed);
    }
  }

  return seeds;
}

// Initialize one independent process object for each worker thread
void MGraniitti::InitMultiMemory() {
  if (CORES < 1 || proc == nullptr || proc->GetModelTune() == nullptr) {
    throw std::logic_error("MGraniitti::InitMultiMemory: require a configured process and at least one worker");
  }

  // Copy the current process before replacing any workers that own it
  const auto copy = [&]() -> std::unique_ptr<MProcess> {
    if (const auto *source = dynamic_cast<const MQuasiElastic *>(proc)) {
      return std::make_unique<MQuasiElastic>(*source);
    }
    if (const auto *source = dynamic_cast<const MFactorized *>(proc)) {
      return std::make_unique<MFactorized>(*source);
    }
    if (const auto *source = dynamic_cast<const MCentral *>(proc)) {
      return std::make_unique<MCentral>(*source);
    }
    if (const auto *source = dynamic_cast<const MCollinear *>(proc)) {
      return std::make_unique<MCollinear>(*source);
    }
    if (const auto *source = dynamic_cast<const MHardDiffraction *>(proc)) {
      return std::make_unique<MHardDiffraction>(*source);
    }
    throw std::logic_error("MGraniitti::InitMultiMemory: process has no worker implementation");
  };
  std::vector<std::unique_ptr<MProcess>> workers;
  workers.reserve(static_cast<std::size_t>(CORES));
  for (int i = 0; i < CORES; ++i) { workers.push_back(copy()); }

  // Use the same immutable numerical tune as the amplitudes
  const auto &numerics_mc = proc->GetModelTune()->Numerics("NUMERICS_MC");
  const std::string method = numerics_mc.at("thread_seeds");
  const int increment_factor = numerics_mc.at("increment_factor");
  const double width_min = numerics_mc.at("width_min");
  const double offshell_max = numerics_mc.at("offshell_max");
  // --------------------------------------------------

  for (const auto &i : indices(workers)) {
    if (workers[i] != nullptr) {
      workers[i]->SetWIDTHMIN(width_min);
      workers[i]->SetDefaultOFFSHELL(offshell_max);
    }
  }

  std::vector<uint32_t> thread_seed(workers.size(), 0);

  if (method == "random") {
    thread_seed =
        GenerateUniqueSeeds(workers[0]->state.random.GetSeed(), workers.size());
  } else if (method == "increment") {
    for (const auto &i : indices(workers)) {
      thread_seed[i] = static_cast<uint32_t>(workers[i]->state.random.GetSeed() +
                                             i * increment_factor + i);
    }
  } else {
    throw std::invalid_argument("InitMultiMemory: Unknown seeding method = " +
                                method);
  }

  std::cout << "InitMultiMemory: Generating unique seeds with 'method' = '"
            << method << "'" << std::endl;
  if (method == "increment") {
    std::cout << "InitMultiMemory: Note: for grid/HPC usage, use sufficiently "
                 "large spacing "
                 "between job (node) main seeds so incremented thread seeds do "
                 "not overlap between jobs."
              << std::endl;
  }
  std::cout << std::endl;

  for (const auto &i : indices(workers)) {
    if (workers[i]->state.lts.process.QMETRICS) { workers[i]->state.lts.qmetrics.Reset(); }
    workers[i]->state.random.SetSeed(thread_seed[i]);
    printf("Thread [%3lu] : RNG seed = %10u \n", i,
           workers[i]->state.random.GetSeed());
  }
  std::cout << std::endl;

  workers.front()->PrintInit(HILJAA);

  // Replace the process copies only after initialization has succeeded
  std::vector<MProcess *> next;
  next.reserve(workers.size());
  for (const auto &worker : workers) { next.push_back(worker.get()); }
  pvec.swap(next);
  proc = pvec.front();
  for (auto &worker : workers) { worker.release(); }
  for (auto *worker : next) { delete worker; }
}

// Combine independent worker densities only after worker sampling has stopped
spin::MQMetrics MGraniitti::QMetrics() const {
  auto result = proc->state.lts.qmetrics;
  for (const auto *worker : pvec) {
    if (worker != proc) { result.Merge(worker->state.lts.qmetrics); }
  }
  return result;
}

// Set common frozen-proposal integration parameters
void MGraniitti::SetIntegrationParam(const INTEGRATIONPARAM &in) {
  integration_param = in;
}

// Set common maximum-weight calibration parameters
void MGraniitti::SetMaxWeightParam(const MAXWEIGHTPARAM &in) {
  std::lock_guard<std::mutex> lock(state_mutex);
  if (max_weight_state.generation_active) {
    throw std::logic_error(
        "MGraniitti::SetMaxWeightParam: maximum-weight policy is frozen "
        "during event generation");
  }
  max_weight_param = in;
}

// Compute an automatic frozen-proposal batch size
std::size_t MGraniitti::AutomaticIntegrationBatch() const {
  const std::uint64_t target =
      std::clamp(integration_param.MIN_SAMPLES / 20, std::uint64_t{1024},
                 std::uint64_t{65536});
  const std::uint64_t threads = static_cast<std::uint64_t>(std::max(CORES, 1));
  const std::uint64_t aligned = ((target + threads - 1) / threads) * threads;
  return static_cast<std::size_t>(aligned);
}

// Compute whether the current state has no dedicated integration samples
bool MGraniitti::ProposalOnlyState() const {
  if (INTEGRATOR == "VEGAS" || INTEGRATOR == "NEUROJAC") {
    return stat.integration_samples == 0;
  }
  return false;
}

// Compute whether weighted generation samples the cross section
// (used e.g. with icetune)
bool MGraniitti::SamplesCrossSectionDuringGeneration() const {
  return NEVENTS > 0 && WEIGHTED && ProposalOnlyState();
}

// Initialize or update maximum-weight calibration after one sample
// envelope = safety_factor w when w exceeds the current maximum
void MGraniitti::ObserveMaximumWeight(double sampling_weight,
                                      SamplingStage stage) {
  if (WEIGHTED || stage != SamplingStage::Integration) {
    return;
  }
  if (!(sampling_weight >= 0.0) || !std::isfinite(sampling_weight)) {
    throw std::invalid_argument(
        "MAX_WEIGHT: unweighted generation requires finite nonnegative "
        "sampling weights");
  }
  if (!max_weight_state.initialized) {
    return;
  }
  if (sampling_weight > max_weight_state.envelope) {
    if (max_weight_param.MODE == MaxWeightMode::Fixed) {
      throw std::invalid_argument(
          "MAX_WEIGHT: sampled weight exceeds the fixed envelope");
    }
    max_weight_state.envelope =
        max_weight_param.SAFETY_FACTOR * sampling_weight;
    max_weight_state.validation_samples = 0;
    ++max_weight_state.envelope_updates;
    return;
  }
  ++max_weight_state.validation_samples;
}

// Observe one sample through the common statistics and failure-reporting path
double MGraniitti::ObserveSample(const gra::MEventWeightState &aux,
                                 const double raw_weight,
                                 const double log_inverse_density,
                                 const SamplingStage stage) {
  if (first_amplitude_failure.empty() && aux.amplitude_failure) {
    first_amplitude_failure = aux.amplitude_failure_message.empty()
                                  ? "unspecified amplitude failure"
                                  : aux.amplitude_failure_message;
    if (!HILJAA) {
      gra::aux::ClearProgress();
      std::cerr << rang::fg::red << "Amplitude failure (first occurrence, use --debug cli flag): "
                << first_amplitude_failure << rang::fg::reset << std::endl;
    }
  }
  return stat.ObserveLogSample(aux, raw_weight, log_inverse_density, stage);
}

// Initialize the estimated maximum-weight envelope after discovery
void MGraniitti::InitializeMaximumWeightEnvelope() {
  if (WEIGHTED || max_weight_state.initialized) {
    return;
  }
  if (max_weight_param.MODE == MaxWeightMode::Fixed) {
    max_weight_state.envelope = max_weight_param.MAX_W_VALUE;
    max_weight_state.initialized = true;
    return;
  }
  if (stat.integration_samples < integration_param.MIN_SAMPLES) {
    return;
  }
  const long double observed = INTEGRATOR == "FLAT"
                                   ? stat.IntegrandWeightMaximum()
                                   : stat.SamplingWeightMaximum();
  if (!(observed > 0.0) || !std::isfinite(observed)) {
    throw std::invalid_argument(
        "MAX_WEIGHT: no finite positive sampling weight was found during "
        "envelope discovery");
  }
  max_weight_state.envelope =
      max_weight_param.SAFETY_FACTOR * static_cast<double>(observed);
  max_weight_state.validation_samples = 0;
  max_weight_state.initialized = true;
}

// Compute the required independent envelope validation count
// N = ceil[log(1-C)/log(1-p_overflow)]
std::uint64_t MGraniitti::RequiredMaximumWeightValidation() const {
  const double numerator = std::log1p(-max_weight_param.CONFIDENCE_LEVEL);
  const double denominator = std::log1p(-max_weight_param.MAX_OVERFLOW_PROB);
  const double required = std::ceil(numerator / denominator);
  if (!(required > 0.0) ||
      required >
          static_cast<double>(std::numeric_limits<std::uint64_t>::max())) {
    throw std::overflow_error(
        "MAX_WEIGHT: requested validation count is outside uint64 range");
  }
  return static_cast<std::uint64_t>(required);
}

// Compute whether common frozen-proposal integration may stop
bool MGraniitti::IntegrationConverged() const {
  if (stat.integration_samples < integration_param.MIN_SAMPLES ||
      !(stat.sigma > 0.0) || !std::isfinite(stat.sigma) ||
      !std::isfinite(stat.sigma_err) ||
      std::abs(stat.sigma_err / stat.sigma) > integration_param.PRECISION) {
    return false;
  }
  if (WEIGHTED) {
    return true;
  }
  if (!max_weight_state.initialized) {
    return false;
  }
  return max_weight_param.MODE == MaxWeightMode::Fixed ||
         max_weight_state.validation_samples >=
             RequiredMaximumWeightValidation();
}

// Reject an integration run that reaches its hard sample limit
void MGraniitti::CheckIntegrationLimit() const {
  if (stat.integration_samples >= integration_param.MAX_SAMPLES &&
      !IntegrationConverged()) {
    throw std::runtime_error(
        "INTEGRATION: max_samples reached before all convergence criteria");
  }
}

// Print the common integration and maximum-weight controls
void MGraniitti::PrintIntegrationSetup() const {
  if (HILJAA) {
    return;
  }
  std::cout << "- min_samples = " << integration_param.MIN_SAMPLES << std::endl;
  std::cout << "- max_samples = " << integration_param.MAX_SAMPLES << std::endl;
  std::cout << "- precision = " << integration_param.PRECISION << std::endl;
  if (WEIGHTED) {
    return;
  }
  if (max_weight_param.MODE == MaxWeightMode::Estimate) {
    std::cout << "- max_w_value = null (estimated)" << std::endl;
    std::cout << "- safety_factor = " << max_weight_param.SAFETY_FACTOR
              << std::endl;
    std::cout << "- max_overflow_prob = " << max_weight_param.MAX_OVERFLOW_PROB
              << std::endl;
    std::cout << "- confidence_level = " << max_weight_param.CONFIDENCE_LEVEL
              << std::endl;
  } else {
    std::cout << "- max_w_value = " << max_weight_param.MAX_W_VALUE
              << std::endl;
  }
  std::cout << "- overflow_action = "
            << OverflowActionName(max_weight_param.OVERFLOW_ACTION)
            << std::endl;
}

// Print the controls which affect VEGAS proposal adaptation
void MGraniitti::PrintVegasAdaptationSetup() const {
  printf("- dimension = %zu\n", vegas.Dimension());
  printf("- bins = %u\n", vparam.BINS);
  printf("- alpha = %g\n", vparam.ALPHA);
  printf("- ncall = %u\n", vparam.NCALL);
  printf("- rounds = %u\n", vparam.ROUNDS);
  printf("- debug = %d\n", vparam.DEBUG);
  printf("- uniform_mix = %g\n", vparam.UNIFORM_MIX);
  printf("- min_support = %zu\n", vparam.MIN_SUPPORT);
  printf("- max_support_calls = %llu\n",
         static_cast<unsigned long long>(vparam.MAX_SUPPORT_CALLS));
  printf("- max_rounds = %u\n", vparam.MAX_ROUNDS);
  printf("- convergence_window = %zu\n", vparam.CONVERGENCE_WINDOW);
  printf("- ess_rel_tolerance = %.4f\n", vparam.ESS_REL_TOLERANCE);
  printf("- best_grid_choice = %s\n",
         VegasGridChoiceName(vparam.BEST_GRID_CHOICE).c_str());
  printf("- automatic_convergence = %s\n",
         aux::bool_cast(vparam.AUTOMATIC_CONVERGENCE).c_str());
}

// Print the controls which affect frozen VEGAS integration
void MGraniitti::PrintVegasIntegrationSetup() const {
  printf("- dimension = %zu\n", vegas.Dimension());
  printf("- bins = %u\n", vparam.BINS);
  printf("- debug = %d\n", vparam.DEBUG);
  printf("- uniform_mix = %g\n", vparam.UNIFORM_MIX);
  PrintIntegrationSetup();
}

// Print NEUROJAC diagnostics and controls using their JSON keys
void MGraniitti::PrintNeuroJacParameters() const {
#ifdef GRANIITTI_USE_LIBTORCH
  std::cout << "- backend = LIBTORCH" << std::endl;
#else
  std::cout << "- backend = INTERNAL" << std::endl;
#endif
  std::cout << "- dimension = " << neurojac.Flow().Dimension() << std::endl;

  const neurojac::MNeuroJacConfig &config = neurojac.Config();
  std::cout << "- uniform_mix = " << config.uniform_mix << std::endl;
  std::cout << "- flow_components = " << config.flow_components << std::endl;
  std::cout << "- vegas_init = " << gra::aux::bool_cast(config.vegas_init)
            << std::endl;
  std::cout << "- vegas_spline_init = "
            << gra::aux::bool_cast(config.vegas_spline_init) << std::endl;
  std::cout << "- vegas_ncall = " << config.vegas_ncall << std::endl;
  std::cout << "- vegas_rounds = " << config.vegas_rounds << std::endl;
  std::cout << "- vegas_sampler_rounds = " << config.vegas_sampler_rounds
            << std::endl;
  std::cout << "- layers = " << config.layers << std::endl;
  std::cout << "- bins = " << config.bins << std::endl;
  std::cout << "- hidden = [";
  for (const auto &i : indices(config.hidden)) {
    std::cout << (i == 0 ? "" : ", ") << config.hidden[i];
  }
  std::cout << "]" << std::endl;
  std::cout << "- activation = "
            << neurojac::FlowActivationName(config.activation) << std::endl;
  std::cout << "- min_bin = " << config.min_bin << std::endl;
  std::cout << "- min_derivative = " << config.min_derivative << std::endl;
  std::cout << "- permutation = "
            << neurojac::FlowPermutationName(config.permutation) << std::endl;
  std::cout << "- buffer_size = " << config.buffer_size << std::endl;
  std::cout << "- batch_size = " << config.batch_size << std::endl;
  std::cout << "- val_size = " << config.validation_size << std::endl;
  std::cout << "- val_patience = " << config.validation_patience << std::endl;
  std::cout << "- val_min_delta = " << config.validation_min_delta << std::endl;
  std::cout << "- replay_frac = " << config.replay_fraction << std::endl;
  std::cout << "- replay_cap = " << config.replay_capacity << std::endl;
  std::cout << "- rounds = " << config.rounds << std::endl;
  std::cout << "- epochs = " << config.epochs << std::endl;
  std::cout << "- loss.name = " << neurojac::FlowLossName(config.loss.kind)
            << std::endl;
  std::cout << "- loss.alpha = [" << config.loss.alpha_initial << ", "
            << config.loss.alpha_final << "]" << std::endl;
  std::cout << "- loss.anneal_frac = " << config.loss.annealing_fraction
            << std::endl;
  std::cout << "- learning_rate = " << config.learning_rate << std::endl;
  std::cout << "- final_learning_rate = " << config.final_learning_rate
            << std::endl;
  std::cout << "- lr_schedule = "
            << neurojac::FlowLRScheduleName(config.lr_schedule) << std::endl;
  std::cout << "- gradient_clip = " << config.gradient_clip << std::endl;
  std::cout << "- adamw_beta1 = " << config.adamw_beta1 << std::endl;
  std::cout << "- adamw_beta2 = " << config.adamw_beta2 << std::endl;
  std::cout << "- adamw_epsilon = " << config.adamw_epsilon << std::endl;
  std::cout << "- adamw_weight_decay = " << config.adamw_weight_decay
            << std::endl;
  std::cout << "- adamw_decay_biases = " << config.adamw_decay_biases
            << std::endl;
}

// Set VEGAS parameters
void MGraniitti::SetVegasParam(const VEGASPARAM &in) { vparam = in; }

// Read parameters from a single JSON file
void MGraniitti::ReadInput(const json &j) {
  ReadGeneralParam(j);
  ReadProcessParam(j);
  ReadModelParam(modelparam);
  ReadRadiativeParam(j);
  ReadNuclearParam(j);
}

// General parameter initialization
void MGraniitti::ReadGeneralParam(const json &j) {
  // JSON block identifier
  const std::string XID = "GENERIC";

  // Setup parameters (order is important)
  SetNumberOfEvents(j.at(XID).at("NEVENTS"));
  SetOutput(j.at(XID).at("OUTPUT"));
  SetFormat(j.at(XID).at("FORMAT"));
  SetWeighted(j.at(XID).at("WEIGHTED"));
  SetIntegrator(j.at(XID).at("INTEGRATOR"));
  SetCores(j.at(XID).at("CORES"));

  // Save for later use
  modelparam = j.at(XID).at("MODELPARAM");
}

// General model parameters initialized from .json file
//
void MGraniitti::ReadModelParam(const std::string &tune) {
  const std::string fullpath = gra::ResolveModelDataFile(tune, "GENERAL.json");

  // Read generic blocks
  modelparam = tune;
  const MModelTunePtr model_tune = MModelTune::Load(fullpath);
  proc->SetModelTune(model_tune);
  proc->SetHelicityConfig(model_tune);

  // The rest are handled by spesific amplitude classes
}

// Read radiative switches after the immutable numerical tune is bound
void MGraniitti::ReadRadiativeParam(const json &j) {
  const json modes =
      j.contains("RADIATIVE") ? j.at("RADIATIVE") : json::object();
  proc->SetRadiative(radiative::ReadConfig(
      modes, proc->GetModelTune()->Numerics("NUMERICS_RADIATIVE")));
}

// Process parameter initialization, Call proc->post_Constructor() after this
void MGraniitti::ReadProcessParam(const json &j) {
  const std::string XID = "SCATTERING";

  // ----------------------------------------------------------------
  // Initialize process

  std::string fullstring = j.at(XID).at("PROCESS");
  FULL_PROCESS_STR       = fullstring;  // Save it for later

  // ----------------------------------------------------------------
  // First separate possible extra arguments by @... ...
  syntax.clear();
  std::vector<std::size_t> markerpos = aux::FindOccurance(fullstring, "@");

  if (markerpos.size() != 0) {
    syntax     = aux::SplitCommands(fullstring.substr(markerpos[0]));
    fullstring = fullstring.substr(0, markerpos[0] - 1);
  }
  // ----------------------------------------------------------------

  // Now separate process and decay parts by "->"
  std::string PROCESS_PART = "";
  std::string DECAY_PART   = "";
  std::size_t pos          = 0;
  std::size_t pos1         = fullstring.find("->");
  std::size_t pos2         = fullstring.find("&>");  // ISOLATEd Phase-Space

  pos = (pos1 != std::string::npos) ? pos1 : std::string::npos;  // Try to find ->
  if (pos == std::string::npos) {
    pos = (pos2 != std::string::npos) ? pos2 : std::string::npos;  // Try to find &>
  }

  if (pos != std::string::npos) {
    PROCESS_PART = fullstring.substr(0, pos - 1);  // beginning to the pos-1
    DECAY_PART   = fullstring.substr(pos + 2);     // from pos+2 to the end
  } else {
    PROCESS_PART = fullstring;  // No decay defined
  }

  // Trim extra spaces away
  gra::aux::TrimExtraSpace(PROCESS_PART);
  gra::aux::TrimExtraSpace(DECAY_PART);

  InitProcessMemory(PROCESS_PART, j.at("GENERIC").at("RNDSEED"));
  const auto process = proc->ProcPtr.SelectedProcess()->Info();
  ValidateProcessCommands(syntax, process);

  // ----------------------------------------------------------------
  // SETUP process: process memory needs to be initialized before!

  proc->SetLHAPDF(j.at(XID).at("LHAPDF"));
  proc->SetScreening(j.at(XID).at("LOOPSCREEN"));
  proc->SetExcitation(j.at(XID).at("NSTARS"));
  proc->SetBeamFrag(j.at(XID).at("BEAMFRAG"));
  proc->SetHistograms(j.at("GENERIC").at("HIST"));

  if (pos2 != std::string::npos) {
    if (PROCESS_PART.find("<F>") == std::string::npos &&
        PROCESS_PART.find("<P>") == std::string::npos) {
      throw std::invalid_argument("MGraniitti::ReadProcessParam: Phase space "
                                  "isolation arrow '&>' to be "
                                  "used only with <F> or <P> class!");
    }
    proc->SetISOLATE(true);
  }

  if (pos != std::string::npos) {
    proc->SetRootDecayMode(pos2 != std::string::npos
                               ? gra::RootDecayMode::Isolated
                               : gra::RootDecayMode::Physical);
  } else {
    proc->SetRootDecayMode(gra::RootDecayMode::None);
  }

  // ----------------------------------------------------------------
  // Decaymode setup
  proc->SetDecayMode(DECAY_PART);

  // Read free resonance parameters for resonance channels or explicit steering
  std::map<std::string, gra::PARAM_RES> RESONANCES;
  const bool load_resonances = process.resonance == ResonanceType::Required ||
      std::any_of(syntax.begin(), syntax.end(), [](const auto &entry) { return entry.id == "RES" || entry.id == "R"; });
  if (load_resonances) {
    // From .json input
    std::vector<std::string> RES = j.at(XID).at("RES").get<std::vector<std::string>>();

    // @ Command syntax override
    for (const auto &i : indices(syntax)) {
      if (syntax[i].id == "RES") {
        std::vector<std::string> temp;
        for (const auto &x : syntax[i].arg) {
          if (aux::ParseBool(x.second, "@RES{" + x.first + "}")) {
            temp.push_back(x.first); // f0_980, ...
          }
        }
        RES = temp;
      }
    }

    if (RES.empty() && process.resonance == ResonanceType::Required) {
      throw std::invalid_argument(
          "MGraniitti::ReadProcessParam: RES and RES+CON processes require an active resonance");
    }

    // Read resonance data
    for (const auto &i : indices(RES)) {
      const std::string str = "RES/" + RES[i] + ".json";
      RESONANCES[RES[i]] = gra::resonance::Read(str, proc->state.random, process.model, modelparam);
    }

    // @ Command syntax override "on-the-flight parameters"
    for (const auto &i : indices(syntax)) {
      if (syntax[i].id == "R") {
        // Take target string
        std::string RESNAME;
        if (syntax[i].target.size() == 1) {
          RESNAME = syntax[i].target[0];
        } else if (syntax[i].target.size() == 0) {
          throw std::invalid_argument(
              "@Syntax error: invalid R[] without any target []");
        } else if (syntax[i].target.size() > 1) {
          throw std::invalid_argument(
              "@Syntax error: invalid R[] with multiple targets inside []");
        }

        // Do we find target
        if (RESONANCES.find(RESNAME) == RESONANCES.end()) {
          throw std::invalid_argument("@Syntax error: invalid R[" + RESNAME +
                                      "] not found");
        }
        auto &res = RESONANCES[RESNAME];
        bool couplings_touched = false;
        bool density_touched = false;

        const int N = res.p.spinX2 + 1;
        MMatrix<std::complex<double>> newrho(N, N, 0.0);

        // Loop over key:val arguments
        for (const auto &x : syntax[i].arg) {
          if (x.first == "M") {
            res.p.mass = ReadResValue(x.second, x.first, true);
            std::cout << rang::fg::green << "@R[" << RESNAME
                      << "] new mass: " << res.p.mass << rang::fg::reset
                      << std::endl;
          }
          if (x.first == "W") {
            res.p.width = ReadResValue(x.second, x.first, true);
            res.p.tau = res.p.width > 0.0 ? PDG::hbar / res.p.width : 0.0;
            std::cout << rang::fg::green << "@R[" << RESNAME
                      << "] new width: " << res.p.width << rang::fg::reset
                      << std::endl;
          }

          // Set couplings
          for (auto &channel : res.TP.channels) {
            for (const auto &n : indices(channel.g_tensor)) {
              if (x.first == ("g" + std::to_string(n))) {
                const double value = ReadResValue(x.second, x.first, false);
                channel.g_tensor[n] = value;
                couplings_touched = true;
              }
            }
          }

          // Set diagonal spin density elements
          for (std::size_t n = 0; n <= (newrho.size_row() - 1) / 2; ++n) {
            const int J = (newrho.size_row() - 1) / 2;
            if (x.first == ("JZ" + std::to_string(n))) {
              const double value = ReadResValue(x.second, x.first, true);

              newrho[newrho.size_row() - 1 - J - n]
                    [newrho.size_row() - 1 - J - n] = value;

              newrho[newrho.size_row() - 1 - J + n]
                    [newrho.size_row() - 1 - J + n] = value;
              density_touched = true;
            }
          }
        }

        // Normalize the trace and save it
        if (density_touched) {
          if (res.UsesUnrestrictedSpinBasis()) {
            throw std::invalid_argument("@R[] JZ populations require MP polarization.mode a_Jz or rho");
          }
          const std::complex<double> trace = newrho.Trace();
          if (std::isfinite(trace.real()) && std::isfinite(trace.imag()) && trace.real() > 0.0) {
            newrho = newrho * (1.0 / trace);
          } else {
            throw std::invalid_argument("MGraniitti::ReadProcessParam: @R[] "
                                        "Spin density requires a finite positive trace");
          }
          res.rho = newrho;

          std::vector<double> probabilities(newrho.size_row(), 0.0);

          for (std::size_t k = 0; k < newrho.size_row(); ++k) {
            const double prob = std::real(newrho[k][k]);
            if (prob < -1e-12) {
              throw std::invalid_argument(
                  "MGraniitti::ReadProcessParam: @R[] Negative diagonal "
                  "probability after normalization");
            }
            probabilities[k] = std::max(0.0, prob);

          }
          std::cout << rang::fg::green << "@R[" << RESNAME
                    << "] minimal-Pomeron spin steering updated from diagonal "
                       "JZ weights"
                    << std::endl;
          if (res.UsesCoherentSpinBasis()) {
            res.a_Jz = gra::resonance::CoherentAJzFromDiagonalWeights(res, probabilities);
            res.rho             = gra::RankOneProjector(res.a_Jz);
            const auto &rho_coh = res.rho;
            const int J = static_cast<int>((res.a_Jz.size() - 1) / 2);
            for (const auto &k : indices(res.a_Jz)) {
              const int Jz = static_cast<int>(k) - J;
              printf("Jz=%2d : %0.6f x exp(i x %0.1f)\n", Jz,
                     std::abs(res.a_Jz[k]), std::arg(res.a_Jz[k]));
            }
            std::cout << "Coherent rho = |a><a|:" << std::endl;
            rho_coh.Print();
          } else {
            // Diagonal density populations do not obey coherent amplitude phases
            res.spin_basis = "rho";
            res.a_Jz.clear();
            res.MP.random_rho = false;
          }
          std::cout << "MP rho steering:" << std::endl;
          res.rho.Print();
          std::cout << rang::fg::reset << std::endl;
        }

        // Print new coupling array
        if (couplings_touched) {
          std::cout << rang::fg::green << "@R[" << RESNAME
                    << "] new production couplings set";
          for (const auto &channel : res.TP.channels) {
            std::cout << ": [" << channel.exchange[0] << ","
                      << channel.exchange[1] << "] g_tensor = [";
            for (const auto &k : indices(channel.g_tensor)) {
              printf("%0.3E", channel.g_tensor[k]);
              if (k < channel.g_tensor.size() - 1) {
                std::cout << ", ";
              }
            }
            std::cout << "]";
          }
          std::cout << rang::fg::reset << std::endl;
        }
      }
    }
    proc->SetResonances(RESONANCES);
  }
  // ----------------------------------------------------------------

  // ----------------------------------------------------------------
  // Command syntax parameters

  for (const auto &i : indices(syntax)) {
    if (syntax[i].id == "FLATAMP") {
      proc->SetFLATAMP(aux::ParseInt(syntax[i].arg["_SINGLET_"], "@FLATAMP"));
    }
    if (syntax[i].id == "FLATMASS2") {
      proc->SetFLATMASS2(
          aux::ParseBool(syntax[i].arg["_SINGLET_"], "@FLATMASS2"));
    }
    if (syntax[i].id == "OFFSHELL") {
      proc->SetOFFSHELL(aux::ParseDouble(syntax[i].arg["_SINGLET_"], "@OFFSHELL"));
    }

    if (syntax[i].id == "SPINGEN") {
      proc->SetSPINGEN(aux::ParseBool(syntax[i].arg["_SINGLET_"], "@SPINGEN"));
    }
    if (syntax[i].id == "SPINDEC") {
      proc->SetSPINDEC(aux::ParseBool(syntax[i].arg["_SINGLET_"], "@SPINDEC"));
    }

    if (syntax[i].id == "QMETRICS") {
      proc->SetQMetrics(aux::ParseBool(syntax[i].arg["_SINGLET_"], "@QMETRICS"));
    }
    if (syntax[i].id == "MP_FRAME") {
      proc->SetMPFrame(syntax[i].arg["_SINGLET_"]);
    }
    if (syntax[i].id == "MMAX") {
      proc->SetMMAX(aux::ParseInt(syntax[i].arg["_SINGLET_"], "@MMAX"));
    }
  }
  // ----------------------------------------------------------------

  // Load after the PDG tables
  const std::vector<std::string> beam =
      j.at(XID).at("BEAM").get<std::vector<std::string>>();
  const std::vector<double> energy =
      j.at(XID).at("ENERGY").get<std::vector<double>>();
  proc->SetInitialState(beam, energy);

  // Now rest of the parameters
  ReadIntegralParam(j);
  ReadGenCuts(j);
  ReadFidCuts(j);
  ReadVetoCuts(j);
}

// Read heavy-ion UPC steering after beam and model-tune initialization
void MGraniitti::ReadNuclearParam(const json &j) {
  const auto initial = proc->GetInitialState();
  const bool ion1 = nuclear::IsNuclearPDG(initial[0].pdg);
  const bool ion2 = nuclear::IsNuclearPDG(initial[1].pdg);
  if ((ion1 || ion2) && !j.contains("NUCLEAR")) {
    throw std::invalid_argument(
        "MGraniitti::ReadNuclearParam: nuclear beams require NUCLEAR "
        "steering");
  }
  if (!j.contains("NUCLEAR")) { return; }
  if (!ion1 && !ion2) {
    throw std::invalid_argument("MGraniitti::ReadNuclearParam: NUCLEAR steering requires a nuclear beam");
  }
  const auto general =
      nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile(modelparam, "GENERAL.json")));
  const auto numerics = nlohmann::json::parse(gra::aux::GetInputData(
      gra::ResolveModelDataFile(modelparam, "NUMERICS.json")));
  auto upc = nuclear::ReadUPCSteering(
      j.at("NUCLEAR"), general.at("PARAM_NUCLEAR"),
      numerics.at("NUMERICS_NUCLEAR"), modelparam);
  proc->SetNuclear(upc);
}

// Read shared integration parameters from the model numerics card
void MGraniitti::ReadIntegratorNumerics(const std::string &filename) {
  const std::string data = gra::aux::GetInputData(filename);
  json document;
  try {
    document = json::parse(data);
  } catch (const std::exception &error) {
    throw std::invalid_argument(
        "MGraniitti::ReadIntegratorNumerics: invalid NUMERICS file: " +
        std::string(error.what()));
  }

  const json &numerics_mc = document.at("NUMERICS_MC");
  fast_adaptation = numerics_mc.at("fast_adaptation").get<bool>();

  const json &maximum = numerics_mc.at("MAX_WEIGHT");
  MAXWEIGHTPARAM maximum_config = max_weight_param;
  maximum_config.SAFETY_FACTOR = maximum.at("safety_factor");
  maximum_config.MAX_OVERFLOW_PROB = maximum.at("max_overflow_prob");
  maximum_config.CONFIDENCE_LEVEL = maximum.at("confidence_level");
  maximum_config.OVERFLOW_ACTION =
      ParseOverflowAction(maximum.at("overflow_action").get<std::string>());
  gra::aux::AssertRange(maximum_config.SAFETY_FACTOR, {1.0, 100.0},
                        "NUMERICS_MC::MAX_WEIGHT::safety_factor", true);
  gra::aux::AssertRange(maximum_config.MAX_OVERFLOW_PROB, {0.0, 1.0},
                        "NUMERICS_MC::MAX_WEIGHT::max_overflow_prob", true);
  gra::aux::AssertRange(maximum_config.CONFIDENCE_LEVEL, {0.0, 1.0},
                        "NUMERICS_MC::MAX_WEIGHT::confidence_level", true);
  if (!(maximum_config.MAX_OVERFLOW_PROB > 0.0) ||
      !(maximum_config.MAX_OVERFLOW_PROB < 1.0) ||
      !(maximum_config.CONFIDENCE_LEVEL > 0.0) ||
      !(maximum_config.CONFIDENCE_LEVEL < 1.0)) {
    throw std::invalid_argument(
        "NUMERICS_MC::MAX_WEIGHT probabilities must be strictly between "
        "zero and one");
  }
  SetMaxWeightParam(maximum_config);

  const json &vegas_numerics = document.at("NUMERICS_VEGAS");
  const double uniform_mix = vegas_numerics.at("uniform_mix");
  if (!(uniform_mix > 0.0 && uniform_mix < 1.0) ||
      !std::isfinite(uniform_mix)) {
    throw std::invalid_argument(
        "NUMERICS_VEGAS::uniform_mix must be finite and strictly between "
        "zero and one");
  }
  vparam.UNIFORM_MIX = uniform_mix;
  vparam.MIN_SUPPORT = vegas_numerics.at("min_support");
  vparam.MAX_SUPPORT_CALLS = vegas_numerics.at("max_support_calls");
  vparam.MAX_ROUNDS = vegas_numerics.at("max_rounds");
  vparam.CONVERGENCE_WINDOW = vegas_numerics.at("convergence_window");
  vparam.ESS_REL_TOLERANCE = vegas_numerics.at("ess_rel_tolerance");
  vparam.BEST_GRID_CHOICE = ParseVegasGridChoice(
      vegas_numerics.at("best_grid_choice").get<std::string>());
  vparam.AUTOMATIC_CONVERGENCE = vegas_numerics.at("automatic_convergence");
  if (vparam.MIN_SUPPORT == 0) {
    throw std::invalid_argument("NUMERICS_VEGAS::min_support must be positive");
  }
  if (vparam.MAX_SUPPORT_CALLS < vparam.MIN_SUPPORT) {
    throw std::invalid_argument(
        "NUMERICS_VEGAS::max_support_calls must not be smaller than "
        "min_support");
  }
  if (vparam.CONVERGENCE_WINDOW < 2 ||
      vparam.MAX_ROUNDS < vparam.CONVERGENCE_WINDOW) {
    throw std::invalid_argument(
        "NUMERICS_VEGAS convergence rounds are inconsistent");
  }
  if (!(vparam.ESS_REL_TOLERANCE > 0.0) || !(vparam.ESS_REL_TOLERANCE <= 1.0) ||
      !std::isfinite(vparam.ESS_REL_TOLERANCE)) {
    throw std::invalid_argument(
        "NUMERICS_VEGAS::ess_rel_tolerance must be positive and finite");
  }
}

// MC integrator parameter initialization
void MGraniitti::ReadIntegralParam(const json &j) {
  const std::string XID = "INTEGRATOR";

  // Numerical loop integration
  const int extra_kt = j.at(XID).at("LOOPSCREEN").at("extra_kt");
  const int extra_phi = j.at(XID).at("LOOPSCREEN").at("extra_phi");
  gra::aux::AssertRange(extra_kt, {-10, 100}, "LOOPSCREEN::extra_kt", true);
  gra::aux::AssertRange(extra_phi, {-10, 100}, "LOOPSCREEN::extra_phi", true);
  proc->eikonal.Numerics.SetLoopDiscretization(extra_kt, extra_phi);

  // Common frozen-proposal integration parameters
  const json &integration = j.at(XID);
  INTEGRATIONPARAM integration_config;
  integration_config.MIN_SAMPLES = integration.at("min_samples");
  integration_config.MAX_SAMPLES = integration.at("max_samples");
  integration_config.PRECISION = integration.at("precision");
  if (integration_config.MIN_SAMPLES < 10) {
    throw std::invalid_argument("INTEGRATOR::min_samples must be at least 10");
  }
  if (integration_config.MAX_SAMPLES < integration_config.MIN_SAMPLES) {
    throw std::invalid_argument(
        "INTEGRATOR::max_samples must not be smaller than min_samples");
  }
  gra::aux::AssertRange(integration_config.PRECISION, {0.0, 1.0},
                        "INTEGRATOR::precision", true);
  SetIntegrationParam(integration_config);

  // Process-specific maximum-weight override
  MAXWEIGHTPARAM maximum_config;
  maximum_config.MODE = ParseMaximumWeightValue(integration.at("max_w_value"),
                                                maximum_config.MAX_W_VALUE);
  SetMaxWeightParam(maximum_config);

  ReadIntegratorNumerics(
      gra::ResolveModelDataFile(modelparam, "NUMERICS.json"));

  // Keep the tailored one-dimensional elastic proposal at the requested rounds
  if (proc->ProcPtr.ISTATE == "X" && proc->ProcPtr.CHANNEL == "EL") {
    vparam.AUTOMATIC_CONVERGENCE = false;
  }

  if (!WEIGHTED && max_weight_param.MODE == MaxWeightMode::Estimate &&
      integration_config.MAX_SAMPLES - integration_config.MIN_SAMPLES <
          RequiredMaximumWeightValidation()) {
    throw std::invalid_argument(
        "INTEGRATOR::max_samples leaves too few samples for maximum weight "
        "validation");
  }

  if (!j.at(XID).contains("VEGAS")) {
    std::cout
        << "MGraniitti::ReadIntegralParam: Did not find VEGAS parameter block "
           "from json input, using default"
        << std::endl;
    return; // Did not find VEGAS parameter block at all, use default
  }

  // VEGAS parameters
  VEGASPARAM vpam = vparam;
  vpam.BINS = j.at(XID).at("VEGAS").at("bins");
  gra::aux::AssertRange(vpam.BINS, {4U, (unsigned int)1e9}, "VEGAS::bins",
                        true);
  if ((vpam.BINS % 2) != 0) {
    throw std::invalid_argument("VEGAS::bins = " + std::to_string(vpam.BINS) +
                                " should be even number");
  }

  vpam.ALPHA = j.at(XID).at("VEGAS").at("alpha");
  gra::aux::AssertRange(vpam.ALPHA, {0.0, 2.0}, "VEGAS::alpha", true);
  vpam.NCALL = j.at(XID).at("VEGAS").at("ncall");
  gra::aux::AssertRange(vpam.NCALL, {50, 10000000U}, "VEGAS::ncall", true);
  vpam.ROUNDS = j.at(XID).at("VEGAS").at("rounds");
  gra::aux::AssertRange(vpam.ROUNDS, {1, (unsigned int)1e9}, "VEGAS::rounds",
                        true);
  if (vpam.ROUNDS > vpam.MAX_ROUNDS) {
    throw std::invalid_argument(
        "VEGAS::rounds must not exceed NUMERICS_VEGAS::max_rounds");
  }
  vpam.DEBUG = j.at(XID).at("VEGAS").at("debug");
  gra::aux::AssertSet(vpam.DEBUG, {-1, 0, 1}, "VEGAS::debug", true);

  SetVegasParam(vpam);
}

// Generator cuts
void MGraniitti::ReadGenCuts(const json &j) {
  gra::GENCUT gcuts;
  const std::string XID = "GENCUTS";

  // Continuum phase space class
  const bool factorized_c = PROCESS.find("<C>") != std::string::npos &&
                            (proc->state.lts.central_phase_space_mode == CentralPhaseSpaceMode::Factorized ||
                             proc->state.lts.central_phase_space_mode == CentralPhaseSpaceMode::HardDiffraction);
  if (factorized_c) {
    ReadSystemGenCutsBlock(j, "<C>", "factorized", gcuts);
  } else if (PROCESS.find("<C>") != std::string::npos) {
    // Daughter rapidity and central system mass
    std::vector<double> rap;
    std::vector<double> mass;
    try {
      rap = j.at(XID).at("<C>").at("Rap").get<std::vector<double>>();
      mass = j.at(XID).at("<C>").at("M").get<std::vector<double>>();
    } catch (const std::exception &error) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: <C> phase space class requires from user: "
          "\"GENCUTS\" : { \"<C>\" : { \"Rap\" : [min, max], "
          "\"M\" : [min, max] }}; underlying "
          "error: " +
          std::string(error.what()));
    }
    gra::aux::AssertCut(rap, "GENCUTS::<C>::Rap", true);
    gra::aux::AssertCutRange(mass, {0.0, 1e32}, "GENCUTS::<C>::M", true);
    gcuts.rap_min = rap[0];
    gcuts.rap_max = rap[1];
    gcuts.M_min = mass[0];
    gcuts.M_max = mass[1];

    const auto &continuum_cuts = j.at(XID).at("<C>");

    // This is optional, forward leg pt
    if (continuum_cuts.contains("Pt")) {
      const std::vector<double> pt =
          continuum_cuts.at("Pt").get<std::vector<double>>();
      gra::aux::AssertCutRange(pt, {0.0, 1e32}, "GENCUTS::<C>::Pt", true);
      gcuts.forward_pt_min = pt[0];
      gcuts.forward_pt_max = pt[1];
    }

    // This is optional, forward leg excitation
    if (continuum_cuts.contains("Xi")) {
      const std::vector<double> Xi =
          continuum_cuts.at("Xi").get<std::vector<double>>();
      gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::<C>::Xi", true);
      gcuts.XI_min = Xi[0];
      gcuts.XI_max = Xi[1];
    }
  }

  // Factorized phase space class
  if (PROCESS.find("<F>") != std::string::npos) {
    ReadSystemGenCutsBlock(j, "<F>", "factorized", gcuts);
  }

  // Collinear parton phase space class
  if (PROCESS.find("<P>") != std::string::npos) {
    ReadSystemGenCutsBlock(j, "<P>", "collinear parton", gcuts);
  }

  // Quasielastic phase space class
  if (PROCESS.find("<Q>") != std::string::npos) {
    const auto &quasielastic_cuts = j.at(XID).at("<Q>");
    const bool elastic_cni =
        proc->ProcPtr.ISTATE == "X" && proc->ProcPtr.CHANNEL == "EL";
    const bool triple_regge =
        proc->ProcPtr.ISTATE == "X" &&
        (proc->ProcPtr.CHANNEL == "SD" || proc->ProcPtr.CHANNEL == "DD");
    if (elastic_cni && !proc->GetScreening()) {
      throw std::invalid_argument("MGraniitti::ReadGenCuts: X[EL]<Q> requires "
                                  "SCATTERING::LOOPSCREEN = true");
    }
    if (!elastic_cni) {
      std::vector<double> Xi;
      try {
        Xi = quasielastic_cuts.at("Xi").get<std::vector<double>>();
      } catch (const std::exception &error) {
        throw std::invalid_argument(
            "MGraniitti::ReadGenCuts: this <Q> process requires "
            "GENCUTS::<Q>::Xi = [min, max]; underlying error: " +
            std::string(error.what()));
      }
      gra::aux::AssertCutRange(Xi, {0.0, 1.0}, "GENCUTS::<Q>::Xi", true);
      gcuts.XI_min = Xi[0];
      gcuts.XI_max = Xi[1];
    }

    // Read the absolute momentum-transfer range in GeV^2
    if ((elastic_cni || triple_regge) &&
        (!quasielastic_cuts.contains("t") ||
         quasielastic_cuts.at("t").is_null())) {
      throw std::invalid_argument(
          "MGraniitti::ReadGenCuts: X[EL/SD/DD]<Q> requires an explicit "
          "finite GENCUTS::<Q>::t [min, max] range");
    }
    if (quasielastic_cuts.contains("t") &&
        !quasielastic_cuts.at("t").is_null()) {
      std::vector<double> abs_t;
      try {
        abs_t = quasielastic_cuts.at("t").get<std::vector<double>>();
      } catch (const std::exception &error) {
        throw std::invalid_argument("MGraniitti::ReadGenCuts: GENCUTS::<Q>::t "
                                    "must be a finite [min, max] range; "
                                    "underlying error: " +
                                    std::string(error.what()));
      }
      const bool valid_range =
          abs_t.size() == 2 && std::isfinite(abs_t[0]) &&
          std::isfinite(abs_t[1]) &&
          (elastic_cni ? abs_t[0] > 0.0 : abs_t[0] >= 0.0) &&
          abs_t[1] > abs_t[0];
      if (!valid_range) {
        const std::string range_requirement =
            elastic_cni ? "with 0 < min < max" : "with 0 <= min < max";
        throw std::invalid_argument(
            "MGraniitti::ReadGenCuts: GENCUTS::<Q>::t must be a finite "
            "[min, max] range " +
            range_requirement);
      }
      gcuts.q_t_abs_min = abs_t[0];
      gcuts.q_t_abs_max = abs_t[1];
    }
  }

  proc->SetGenCuts(gcuts);
}

// Fiducial cuts
void MGraniitti::ReadFidCuts(const json &j) {
  // Fiducial cuts
  gra::FIDCUT fcuts;
  const std::string XID = "FIDCUTS";

  if (!j.contains(XID)) {
    return;
  }
  const json &fidcuts_block = j.at(XID);
  if (!fidcuts_block.is_object()) {
    throw std::invalid_argument(
        "MGraniitti::ReadFidCuts: FIDCUTS must be an object");
  }
  AssertFiducialBlockKeys(fidcuts_block, "FIDCUTS",
                          {"active", "CENTRAL", "FORWARD", "USERCUTS"});
  fcuts.active = fidcuts_block.at("active");

  // Central final-state, system and PDG-selected cuts
  ReadCentralFiducialCuts(j.at(XID), fcuts);

  // Forward
  {
    const json *forward = OptionalFiducialBlock(fidcuts_block, "FORWARD");
    if (forward != nullptr) {
      AssertFiducialBlockKeys(*forward, "FIDCUTS::FORWARD",
                              {"t", "t1", "t2", "M", "Xi", "dPhi"});
    }
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "t",
                              fcuts.forward_t_active, fcuts.forward_t_min,
                              fcuts.forward_t_max,
                              {0.0, std::numeric_limits<double>::max()});
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "t1",
                              fcuts.forward_t1.active, fcuts.forward_t1.min, fcuts.forward_t1.max,
                              {0.0, std::numeric_limits<double>::max()});
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "t2",
                              fcuts.forward_t2.active, fcuts.forward_t2.min, fcuts.forward_t2.max,
                              {0.0, std::numeric_limits<double>::max()});
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "M",
                              fcuts.forward_M_active, fcuts.forward_M_min,
                              fcuts.forward_M_max,
                              {0.0, std::numeric_limits<double>::max()});
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "Xi",
                              fcuts.forward_xi_active, fcuts.forward_xi_min,
                              fcuts.forward_xi_max, {0.0, 1.0});
    ReadOptionalFiducialRange(forward, "FIDCUTS::FORWARD", "dPhi",
                              fcuts.forward_dPhi_active, fcuts.forward_dPhi_min,
                              fcuts.forward_dPhi_max, {0.0, 180.0});
  }

  // Read and validate the complete fiducial configuration before storing it
  const std::int64_t usercuts =
      fidcuts_block.contains("USERCUTS")
          ? ParseUserCutID(fidcuts_block.at("USERCUTS"))
          : 0;
  proc->ValidateFiducialCuts(fcuts, usercuts);
  proc->SetFidCuts(fcuts);
  proc->SetUserCuts(usercuts);
}

// Read veto-cut domains and their particle-source selectors
void MGraniitti::ReadVetoCuts(const json &j) {
  // Veto cuts
  gra::VETOCUT veto;
  const std::string XID = "VETOCUTS";

  if (!j.contains(XID)) {
    return;
  }
  veto.active = j.at(XID).at("active");
  const auto &veto_block = j.at(XID);

  // Find domains
  for (std::size_t i = 0; i < 100; ++i) {
    const std::string NUMBER = std::to_string(i);

    gra::VETODOMAIN domain;
    std::vector<double> eta;
    std::vector<double> pt;

    // Stop at the first domain index that is not configured
    if (!veto_block.contains(NUMBER)) {
      break;
    }
    const auto &domain_block = veto_block.at(NUMBER);
    for (auto entry = domain_block.begin(); entry != domain_block.end();
         ++entry) {
      if (entry.key() != "Eta" && entry.key() != "Pt" &&
          entry.key() != "sources" && entry.key() != "charge") {
        throw std::invalid_argument("VETOCUTS::" + NUMBER +
                                    " has unknown key " + entry.key());
      }
    }

    eta = domain_block.at("Eta").get<std::vector<double>>();
    pt = domain_block.at("Pt").get<std::vector<double>>();

    gra::aux::AssertCut(eta, "VETOCUTS::" + NUMBER + "::Eta", true);
    domain.eta_min = eta[0];
    domain.eta_max = eta[1];

    gra::aux::AssertCut(pt, "VETOCUTS::" + NUMBER + "::Pt", true);
    domain.pt_min = pt[0];
    domain.pt_max = pt[1];

    // Restrict an optional domain to named scattering-system sources
    if (domain_block.contains("sources")) {
      domain.source_forward = false;
      domain.source_central = false;
      const auto sources =
          domain_block.at("sources").get<std::vector<std::string>>();
      if (sources.empty()) {
        throw std::invalid_argument(
            "VETOCUTS::" + NUMBER +
            "::sources must contain at least one source");
      }
      for (const auto &source : sources) {
        if (source == "forward") {
          if (domain.source_forward) {
            throw std::invalid_argument("VETOCUTS::" + NUMBER +
                                        "::sources contains duplicate forward");
          }
          domain.source_forward = true;
        } else if (source == "central") {
          if (domain.source_central) {
            throw std::invalid_argument("VETOCUTS::" + NUMBER +
                                        "::sources contains duplicate central");
          }
          domain.source_central = true;
        } else {
          throw std::invalid_argument("VETOCUTS::" + NUMBER +
                                      "::sources supports forward and central");
        }
      }
    }

    // Select the particle charge category requested by the domain
    if (domain_block.contains("charge")) {
      const auto charge = domain_block.at("charge").get<std::string>();
      if (charge == "any") {
        domain.charge = gra::VetoCharge::Any;
      } else if (charge == "charged") {
        domain.charge = gra::VetoCharge::Charged;
      } else if (charge == "neutral") {
        domain.charge = gra::VetoCharge::Neutral;
      } else {
        throw std::invalid_argument(
            "VETOCUTS::" + NUMBER +
            "::charge supports any, charged and neutral");
      }
    }

    veto.cuts.push_back(domain);
  }

  // Set fiducial cuts
  proc->SetVetoCuts(veto);
}

// Get maximum integration weight
double MGraniitti::GetMaxweight() const {
  std::lock_guard<std::mutex> lock(state_mutex);
  if (max_weight_state.initialized) {
    return max_weight_state.envelope;
  }
  if (INTEGRATOR == "VEGAS" || INTEGRATOR == "NEUROJAC") {
    return static_cast<double>(stat.SamplingWeightMaximum());
  } else {
    return static_cast<double>(stat.IntegrandWeightMaximum());
  }
}

// Compute the immutable maximum weight captured at generation start
double MGraniitti::GetGenerationMaxweight() const {
  if (WEIGHTED) {
    return GetMaxweight();
  }
  std::lock_guard<std::mutex> lock(state_mutex);
  return max_weight_state.GenerationEnvelope();
}

// Set maximum weight for the integration process
void MGraniitti::SetMaxweight(double weight) {
  std::lock_guard<std::mutex> lock(state_mutex);
  if (max_weight_state.generation_active) {
    throw std::logic_error(
        "MGraniitti::SetMaxweight: maximum weight is frozen during event "
        "generation");
  }
  if (weight > 0 && std::isfinite(weight)) {
    max_weight_state.envelope = weight;
    max_weight_state.initialized = true;
  } else {
    std::string str =
        "MGraniitti::SetMaxweight: Maximum weight: " + std::to_string(weight) +
        " not valid!";
    throw std::invalid_argument(str);
  }
}

// Print the generator banner and runtime information
void MGraniitti::PrintInit() const {
  if (!HILJAA) {
    gra::aux::PrintFlashScreen(rang::fg::magenta);
    std::cout << rang::style::bold
              << "GRANIITTI - Monte Carlo event generator for "
                 "high energy diffraction"
              << rang::style::reset << std::endl
              << std::endl;
    gra::aux::PrintVersion();
    gra::aux::PrintBar("-");

    const double GB = pow3(1024.0);
    printf("Running on %s (%d CORE / %0.2f GB RAM) at %s \n",
           gra::aux::HostName().c_str(), std::thread::hardware_concurrency(),
           gra::aux::TotalSystemMemory() / GB, gra::aux::DateTime().c_str());
    int64_t size = 0;
    int64_t free = 0;
    int64_t used = 0;
    gra::aux::GetDiskUsage("/", size, free, used);
    printf("Path '/': size | used | free = %0.1f | %0.1f | %0.1f GB \n",
           size / GB, used / GB, free / GB);
    std::cout << "Program path: " << gra::aux::GetBasePath(2) << std::endl;
    std::cout << gra::aux::SystemName() << std::endl;
    gra::aux::PrintBar("-");
  }
}

// Main integrator steering
//
void MGraniitti::CallIntegrator(unsigned int N) {

  if (IntegrationProposalOnly() && INTEGRATOR == "FLAT") {
    throw std::invalid_argument(
        "MGraniitti::CallIntegrator: -n -1 requires VEGAS or NEUROJAC");
  }
  if (IntegrationProposalOnly() && MC_VGRID_INPUT != "null") {
    throw std::invalid_argument(
        "MGraniitti::CallIntegrator: -n -1 cannot be combined with a "
        "pre-computed integration proposal");
  }

  std::cout << std::endl;
  const std::string setup_title =
      N > 0 ? "Event generation:" : "Generic setup:";
  std::cout << rang::style::bold << setup_title << rang::style::reset
            << std::endl
            << std::endl;

  std::cout << "Output file:            " << OUTPUT << std::endl;
  std::cout << "Output format:          " << FORMAT << std::endl;
  std::cout << "Multithreading:         " << CORES << std::endl;
  std::cout << "Integrator:             " << INTEGRATOR << std::endl;
  std::cout << "Number of events:       " << NEVENTS << std::endl;
  std::cout << "Parameter setup:        " << modelparam << std::endl;

  std::string str = (WEIGHTED == true) ? "weighted" : "unweighted";
  std::cout << rang::fg::green << "Generation mode:        " << str
            << rang::fg::reset << std::endl;
  std::cout << std::endl;

  // =====================================================================

  // Initialize global clock
  if (N == 0) {
    global_tictoc = MTimer(true);
    first_amplitude_failure.clear();
  }

  // Sample the phase space
  if (INTEGRATOR == "VEGAS") {
    SampleVegas(N);
  } else if (INTEGRATOR == "FLAT") {
    SampleFlat(N);
  } else if (INTEGRATOR == "NEUROJAC") {
    SampleNeuroJac(N);
  } else {
    throw std::invalid_argument(
        "MGraniitti::CallIntegrator: Unknown INTEGRATOR = " + INTEGRATOR +
        " (use VEGAS, NEUROJAC, FLAT)");
  }
}

// Initialize and generate events using VEGAS MC
//
void MGraniitti::SampleVegas(unsigned int N) {
  if (N == 0) {
    InitMultiMemory();
    GMODE = 0; // Pure integration

    // ******************************************************
    // Set the unit-hypercube phase-space dimension
    vegas.SetDimension(proc->GetdLIPSDim());
    // ******************************************************
  }
  if (N > 0) {
    GMODE = 1; // Event generation
  }

  // A. Initalize from a pre-computed file
  if (MC_VGRID_INPUT != "null" && GMODE == 0) {
    ReadVGRID(MC_VGRID_INPUT);

    // Reuse full integration state or complete a proposal-only grid as needed
    if (stat.integration_samples > 0) {
      PrintStatistics(N);
      return;
    }
    if (SamplesCrossSectionDuringGeneration()) {
      return;
    }
    VEGAS(VEGASStage::Integration, AutomaticIntegrationBatch(), 0, N, vparam);
    SaveVGRID();
    return;
  }

  // B. On-the-flight initialization
  else if (GMODE == 0) {
    // Adapt the proposal with the configured screening mode
    VEGAS(VEGASStage::Adaptation, vparam.NCALL, vparam.ROUNDS, N, vparam);

    // Save the adapted proposal without production integration
    if (IntegrationProposalOnly()) {
      SaveVGRID();
      return;
    }

    // Reset adaptation statistics before frozen-proposal integration
    VEGAS(VEGASStage::Integration, AutomaticIntegrationBatch(), 0, N, vparam);

    // Save the frozen grid after integration
    SaveVGRID();
  }

  // Generate events with the frozen grid from integration
  else if (GMODE == 1) {
    if (!WEIGHTED && !max_weight_state.initialized) {
      throw std::logic_error(
          "MGraniitti::SampleVegas: unweighted generation requires a "
          "calibrated maximum weight");
    }
    VEGAS(VEGASStage::Generation, AutomaticIntegrationBatch(), 0, N, vparam);

    // Persist the integral improved by all frozen-proposal generation trials
    if (!max_weight_state.generation_overflow ||
        max_weight_param.OVERFLOW_ACTION == OverflowAction::KeepAsWeighted) {
      SaveVGRID();
    }
  }

  else {
    throw std::invalid_argument(
        "MGraniitti::SampleVegas: Unknown operation mode!");
  }
}

// Multithreaded VEGAS integrator
// [close to optimal importance sampling iff
//  integrand factorizes dimension by dimension]
//
// Original algorithm from:
// [REFERENCE: Lepage, G.P. Journal of Computational Physics, 1978]
// en.wikipedia.org/wiki/VEGAS_algorithm
//
std::uint64_t MGraniitti::VEGAS(VEGASStage stage, std::uint64_t calls,
                                unsigned int rounds, unsigned int N,
                                const VEGASPARAM &param) {
  vegas.Init(stage, calls, rounds, param);
  if (stage == VEGASStage::Integration &&
      max_weight_param.MODE == MaxWeightMode::Fixed) {
    InitializeMaximumWeightEnvelope();
  }

  if (stage == VEGASStage::Adaptation && !HILJAA) {
    gra::aux::ClearProgress();
    std::cout << rang::style::bold;
    if (param.AUTOMATIC_CONVERGENCE) {
      printf("VEGAS adaptation (minimum %u, maximum %u rounds): \n\n", rounds,
             param.MAX_ROUNDS);
    } else {
      printf("VEGAS adaptation (%u rounds): \n\n", rounds);
    }
    std::cout << rang::style::reset;
    PrintVegasAdaptationSetup();
    std::cout << std::endl;
  }

  if (stage == VEGASStage::Integration && !HILJAA) {
    gra::aux::ClearProgress();
    gra::aux::PrintBar("-");
    std::cout << rang::style::bold << "VEGAS integration:" << rang::style::reset
              << std::endl
              << std::endl;
    PrintVegasIntegrationSetup();
    gra::aux::PrintBar("-");
  }

  // Reset local timers
  local_tictoc = MTimer(true);

  MTimer stime = MTimer(true); // For statusprint
  atime = MTimer(true);        // For progressbar

  // VEGAS grid iterations
  for (std::size_t iter = 0;;) {
    std::uint64_t batch_calls =
        stage == VEGASStage::Adaptation ? vegas.NextAdaptationCalls() : calls;
    if (stage == VEGASStage::Integration) {
      if (stat.integration_samples >= integration_param.MAX_SAMPLES) {
        CheckIntegrationLimit();
        break;
      }
      const std::uint64_t remaining =
          integration_param.MAX_SAMPLES - stat.integration_samples;
      batch_calls = std::min(batch_calls, remaining);
    }
    if (batch_calls == 0 ||
        batch_calls > std::numeric_limits<unsigned int>::max()) {
      CheckIntegrationLimit();
      break;
    }
    const std::vector<unsigned int> local_calls =
        vegas.LocalCalls(static_cast<unsigned int>(batch_calls), CORES);

    if (stage == VEGASStage::Adaptation) {
      vegas.BeginAdaptationBatch();
    } else if (stage == VEGASStage::Integration) {
      stat.ResetIntegrationBatch();
    }
    {
      std::lock_guard<std::mutex> lock(state_mutex);
      worker_exception = nullptr;
    }

    std::vector<std::thread> threads;
    std::uint64_t first_sample = stat.integration_samples;
    for (int tid = 0; tid < CORES; ++tid) {
      threads.push_back(std::thread(
          [=, this] { VEGASMultiThread(N, tid, stage, local_calls[tid], first_sample); }));
      first_sample += local_calls[tid];
    }
    for (auto &t : threads) {
      t.join();
    }
    if (worker_exception) {
      std::rethrow_exception(worker_exception);
    }
    if (stage == VEGASStage::Generation &&
        max_weight_state.generation_overflow &&
        max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted) {
      return vegas.AdaptationNonzeroCalls();
    }

    // Report every adaptation round without integration statistics
    if (stage == VEGASStage::Adaptation) {
      const VEGASAdaptationReport report =
          vegas.FinishAdaptationBatch(batch_calls);

      if (!HILJAA) {
        PrintVegasAdaptationStatus(report);
        const double completed =
            report.status == VEGASAdaptationStatus::InsufficientSupport
                ? static_cast<double>(report.round - 1)
                : static_cast<double>(report.round);
        gra::aux::PrintProgress(
            std::min(completed / static_cast<double>(rounds), 1.0));
      }

      if (report.status == VEGASAdaptationStatus::InsufficientSupport) {
        continue;
      }
      if (report.status == VEGASAdaptationStatus::Converged) {
        if (!HILJAA) {
          gra::aux::ClearProgress();
          if (param.BEST_GRID_CHOICE == VEGASGridChoice::Best) {
            printf("VEGAS:: adaptation converged after %zu rounds, retaining "
                   "best ESS/N = %.4f grid from round %zu\n",
                   report.round, report.best_ess_fraction, report.best_round);
          } else {
            printf("VEGAS:: adaptation converged after %zu rounds, retaining "
                   "last adapted grid\n",
                   report.round);
          }
        }
        break;
      }
      if (report.status == VEGASAdaptationStatus::RoundLimit) {
        if (!HILJAA && param.AUTOMATIC_CONVERGENCE) {
          gra::aux::ClearProgress();
          if (param.BEST_GRID_CHOICE == VEGASGridChoice::Best) {
            printf("VEGAS:: adaptation reached max_rounds = %u without "
                   "convergence, retaining best ESS/N = %.4f grid from round "
                   "%zu\n",
                   param.MAX_ROUNDS, report.best_ess_fraction,
                   report.best_round);
          } else {
            printf("VEGAS:: adaptation reached max_rounds = %u without "
                   "convergence, retaining last adapted grid\n",
                   param.MAX_ROUNDS);
          }
        }
        break;
      }
      ++iter;
      continue;
    }

    stat.CalculateCrossSection();
    if (stage == VEGASStage::Integration) {
      const IntegrationBatchSummary batch = stat.UpdateIntegrationChi2();

      // Reject invalid frozen-proposal integration statistics
      const double relative_error =
          stat.sigma > 0.0 ? std::abs(stat.sigma_err / stat.sigma)
                           : std::numeric_limits<double>::infinity();
      if (!batch.valid || !(stat.sigma > 0.0) || !std::isfinite(stat.sigma) ||
          !std::isfinite(stat.sigma_err) || !std::isfinite(stat.chi2) ||
          !std::isfinite(relative_error)) {
        gra::aux::ClearProgress();
        throw std::invalid_argument(
            "VEGAS:: Invalid zero-support or non-finite integral estimate: "
            "check the process energy, phase space and cuts or use <F> phase "
            "space");
      }
      InitializeMaximumWeightEnvelope();

      if (param.DEBUG >= 0) {
        gra::aux::ClearProgress();
        printf("VEGAS:: local iter = %4lu integral = %0.5LE +- std = "
               "%0.5LE \t [global integral = %0.5E +- std = %0.5E] \t "
               "ESS/N = %0.3E chi2this = %0.2f \n",
               iter + 1, batch.mean, std::sqrt(batch.variance_of_mean),
               stat.sigma, stat.sigma_err, batch.ess_fraction, stat.chi2);
      }
    }
    if (param.DEBUG > 0) {
      const VEGASData &data = vegas.Data();
      for (std::size_t dimension = 0; dimension < data.FDIM; ++dimension) {
        printf("VEGAS:: data for dimension j = %lu (FDIM = %u) \n", dimension,
               data.FDIM);
        for (std::size_t bin = 0; bin < param.BINS; ++bin) {
          printf("xmat[%3lu][j] = %0.5E, fmat[%3lu][j] = %0.5E \n", bin,
                 data.xmat[bin][dimension], bin, data.fmat[bin][dimension]);
        }
      }
    }

    // Report frozen-proposal integration progress
    if (stage == VEGASStage::Integration) {
      if (stime.ElapsedSec() > 2.0) {
        PrintStatus(static_cast<unsigned int>(std::min<std::uint64_t>(
                        stat.integration_samples,
                        std::numeric_limits<unsigned int>::max())),
                    N, local_tictoc, -1.0);
        stime.Reset();
      }
      if (atime.ElapsedSec() > 0.01) {
        const double progress =
            stat.integration_samples /
            static_cast<double>(integration_param.MIN_SAMPLES);
        gra::aux::PrintProgress(std::min(progress, 1.0));
        atime.Reset();
      }
    }

    // Derive common histogram bounds from the first integration batch
    if (stage == VEGASStage::Integration && iter == 0) {
      UnifyHistogramBounds();
    }

    if (stage == VEGASStage::Integration && IntegrationConverged()) {
      if (!HILJAA) {
        printf("\nVEGAS:: integration criteria reached after %llu samples\n",
               static_cast<unsigned long long>(stat.integration_samples));
      }
      break;
    }
    if (stage == VEGASStage::Integration) {
      CheckIntegrationLimit();
    }

    if (stage == VEGASStage::Generation && stat.generated >= N) {
      break;
    }
    ++iter;
  }

  if (stage != VEGASStage::Adaptation) {
    PrintStatistics(N);
  }
  return vegas.AdaptationNonzeroCalls();
}

// This is called once for every VEGAS grid iteration
void MGraniitti::VEGASMultiThread(unsigned int N, unsigned int THREAD_ID,
                                  VEGASStage stage, unsigned int LOCALcalls, const std::uint64_t first_sample) {
  try {
    MProcess *worker = pvec.at(THREAD_ID);
    const bool qmetrics = stage == VEGASStage::Integration && worker->state.lts.process.QMETRICS;
    for (std::size_t k = 0; k < LOCALcalls;
         ++k) { // ** LOCALcalls ~= calls / CORES **
      if (stage == VEGASStage::Generation) {
        std::lock_guard<std::mutex> lock(state_mutex);
        if (stat.generated >= static_cast<unsigned int>(GetNumberOfEvents()) ||
            (max_weight_state.generation_overflow &&
             max_weight_param.OVERFLOW_ACTION !=
                 OverflowAction::KeepAsWeighted)) {
          break;
        }
      }

      VEGASSample sample =
          vegas.Sample([worker]() { return worker->state.random.U(0, 1); });

      // *******************************************************************
      // ****** Call the process under integration to get the weight *******

      gra::MEventWeightState aux;
      aux.log_inverse_density = sample.log_inverse_density;
      aux.qmetrics = qmetrics;
      if (qmetrics) { aux.sample_index = first_sample + k; }
      if (stage == VEGASStage::Adaptation) {
        aux.ConfigureAdaptation(fast_adaptation);
      }

      const double W = worker->EventWeight(sample.point, aux);

      // *******************************************************************

      const SamplingStage sampling_stage =
          stage == VEGASStage::Adaptation
              ? SamplingStage::Adaptation
              : (stage == VEGASStage::Integration ? SamplingStage::Integration
                                                  : SamplingStage::Generation);
      double f = 0.0;

      // Accumulate shared integration state under an exception-safe lock
      {
        std::lock_guard<std::mutex> lock(state_mutex);
        f = ObserveSample(aux, W, sample.log_inverse_density, sampling_stage);
        if (stage == VEGASStage::Adaptation) {
          vegas.AccumulateAdaptation(f, sample.indices);
        }
        ObserveMaximumWeight(f, sampling_stage);
      }

      // ----------------------------------------------------------
      // Event generation mode
      if (stage == VEGASStage::Generation) {
        // Enough events (>= instead of == for safety)
        unsigned int generated = 0;
        {
          std::lock_guard<std::mutex> lock(state_mutex);
          generated = stat.generated;
        }
        if (generated >= static_cast<unsigned int>(GetNumberOfEvents())) {
          break;
        }

        // Event trial
        SaveEvent(worker, f, GetGenerationMaxweight(), aux);

        if (THREAD_ID == 0 && atime.ElapsedSec() > 0.5) {
          {
            std::lock_guard<std::mutex> lock(state_mutex);
            generated = stat.generated;
          }
          PrintStatus(generated, N, local_tictoc, 10.0);
          gra::aux::PrintProgress(generated / static_cast<double>(N));
          atime.Reset();
        }
      }
    } // calls loop
  } catch (...) {
    // Preserve the first worker exception without racing other workers
    std::lock_guard<std::mutex> lock(state_mutex);
    if (!worker_exception) {
      worker_exception = std::current_exception();
    }
  }
}

// Generate events using plain simple MC (for reference/DEBUG purposes)
//
void MGraniitti::SampleFlat(unsigned int N) {
  // Integration mode
  if (N == 0) {
    GMODE = 0;
    proc->PrintInit(HILJAA);
    InitializeMaximumWeightEnvelope();
    stat.ResetIntegrationBatch();
    if (!HILJAA) {
      std::cout << rang::style::bold
                << "FLAT integration:" << rang::style::reset << std::endl
                << std::endl;
      PrintIntegrationSetup();
      std::cout << std::endl;
    }
  }
  // Event generation mode
  if (N > 0) {
    GMODE = 1;
  }

  // Get dimension of the phase space
  const unsigned int dim = proc->GetdLIPSDim();
  std::vector<double> randvec(dim, 0.0);

  // Reset local timer
  local_tictoc = MTimer(true);

  // Progressbar
  atime = MTimer(true);

  // Event loop
  while (true) {
    // ** Generate new random numbers **
    for (const auto &i : indices(randvec)) {
      randvec[i] = proc->state.random.U(0, 1);
    }

    // Generate event
    // Used for in-out control of the process
    gra::MEventWeightState aux;
    aux.log_inverse_density = 0.0;
    aux.qmetrics = GMODE == 0 && proc->state.lts.process.QMETRICS;
    aux.sample_index = stat.integration_samples;
    aux.adaptation_mode = false;
    const double W = proc->EventWeight(randvec, aux);

    const SamplingStage sampling_stage =
        GMODE == 0 ? SamplingStage::Integration : SamplingStage::Generation;
    const double sampling_weight = ObserveSample(aux, W, 0.0, sampling_stage);
    ObserveMaximumWeight(sampling_weight, sampling_stage);

    // Initialization
    if (GMODE == 0) {
      stat.CalculateCrossSection();
      if (stat.integration_samples % AutomaticIntegrationBatch() == 0) {
        stat.UpdateIntegrationChi2();
        stat.ResetIntegrationBatch();
      }
      PrintStatus(static_cast<unsigned int>(std::min<std::uint64_t>(
                      stat.integration_samples,
                      std::numeric_limits<unsigned int>::max())),
                  N, local_tictoc, 10.0);

      if (stat.integration_samples >= integration_param.MIN_SAMPLES) {
        if (!stat.ValidCrossSection()) {
          throw std::invalid_argument(
              "FLAT:: Invalid zero-support or non-finite integral estimate: "
              "check the process energy, phase space and cuts");
        }
        InitializeMaximumWeightEnvelope();
        if (IntegrationConverged()) {
          break;
        }
      }
      CheckIntegrationLimit();

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(std::min(
            1.0, stat.integration_samples /
                     static_cast<double>(integration_param.MIN_SAMPLES)));
        atime.Reset();
      }
    }

    // Event generation mode
    if (GMODE == 1) {
      SaveEvent(proc, sampling_weight, GetGenerationMaxweight(), aux);
      if (max_weight_state.generation_overflow &&
          max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted) {
        return;
      }
      PrintStatus(stat.generated, N, local_tictoc, 10.0);
      if (stat.generated >= N) {
        break;
      }

      // Progressbar
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
        atime.Reset();
      }
    }
  }
  PrintStatus(stat.generated, N, local_tictoc, -1.0);
  PrintStatistics(N);
}

// Collect one thread-local NEUROJAC adaptation buffer
//
void MGraniitti::NeuroJacCollect(
    std::size_t thread_id, std::size_t calls, NeuroJacSampler sampler,
    std::vector<neurojac::MFlowTrainingEvent> &events) {
  events.clear();
  events.reserve(calls);
  MProcess *worker = pvec.at(thread_id);
  if (sampler == NeuroJacSampler::Flow) {
    const neurojac::MFlowSampleBatch samples =
        neurojac.SampleBatch(worker->state.random, calls);
    for (std::size_t event = 0; event < samples.Size(); ++event) {
      std::vector<double> point(
          samples.points.begin() +
              static_cast<std::ptrdiff_t>(event * samples.dimension),
          samples.points.begin() +
              static_cast<std::ptrdiff_t>((event + 1) * samples.dimension));
      const double log_density = samples.log_densities[event];
      gra::MEventWeightState aux;
      aux.log_inverse_density = -log_density;
      aux.ConfigureAdaptation(fast_adaptation);

      // Evaluate the configured adaptation amplitude
      const double target = worker->EventWeight(point, aux);
      events.push_back({std::move(point), target, log_density});
    }
    return;
  }

  for (std::size_t event = 0; event < calls; ++event) {
    VEGASSample sample =
        vegas.Sample([worker]() { return worker->state.random.U(0, 1); });
    gra::MEventWeightState aux;
    aux.log_inverse_density = sample.log_inverse_density;
    aux.ConfigureAdaptation(fast_adaptation);

    // Evaluate the configured adaptation amplitude
    const double target = worker->EventWeight(sample.point, aux);
    events.push_back(
        {std::move(sample.point), target, -sample.log_inverse_density});
  }
}

// Collect one multithreaded NEUROJAC adaptation batch
std::vector<neurojac::MFlowTrainingEvent>
MGraniitti::NeuroJacCollectBatch(std::size_t calls, NeuroJacSampler sampler) {
  std::vector<std::size_t> worker_calls(
      static_cast<std::size_t>(CORES), calls / static_cast<std::size_t>(CORES));
  for (std::size_t i = 0; i < calls % static_cast<std::size_t>(CORES); ++i) {
    ++worker_calls[i];
  }
  std::vector<std::vector<neurojac::MFlowTrainingEvent>> worker_events(
      static_cast<std::size_t>(CORES));

  {
    std::lock_guard<std::mutex> lock(state_mutex);
    worker_exception = nullptr;
  }
  std::vector<std::thread> threads;
  for (int thread_id = 0; thread_id < CORES; ++thread_id) {
    threads.emplace_back([&, thread_id] {
      try {
        NeuroJacCollect(static_cast<std::size_t>(thread_id),
                        worker_calls[thread_id], sampler,
                        worker_events[thread_id]);
      } catch (...) {
        std::lock_guard<std::mutex> lock(state_mutex);
        if (!worker_exception) {
          worker_exception = std::current_exception();
        }
      }
    });
  }
  for (std::thread &thread : threads) {
    thread.join();
  }
  if (worker_exception) {
    std::rethrow_exception(worker_exception);
  }

  std::vector<neurojac::MFlowTrainingEvent> events;
  events.reserve(calls);
  for (auto &local_events : worker_events) {
    events.insert(events.end(), std::make_move_iterator(local_events.begin()),
                  std::make_move_iterator(local_events.end()));
  }
  return events;
}

// Retry fixed-size NEUROJAC bootstrap batches until training has support
//
std::vector<neurojac::MFlowTrainingEvent>
MGraniitti::NeuroJacCollectTraining(std::size_t &evaluations,
                                    NeuroJacSampler sampler) {
  constexpr std::size_t MAX_FLOW_EMPTY_BATCHES = 1024;
  constexpr std::size_t MAX_VEGAS_EMPTY_BATCHES = 8;
  const std::size_t max_empty_batches = sampler == NeuroJacSampler::FrozenVegas
                                            ? MAX_VEGAS_EMPTY_BATCHES
                                            : MAX_FLOW_EMPTY_BATCHES;
  const std::size_t batch_calls = neurojac.FreshEventCount();

  for (std::size_t batch = 0; batch < max_empty_batches; ++batch) {
    std::vector<neurojac::MFlowTrainingEvent> events =
        NeuroJacCollectBatch(batch_calls, sampler);
    evaluations += batch_calls;
    const bool has_support =
        std::any_of(events.begin(), events.end(),
                    [](const neurojac::MFlowTrainingEvent &event) {
                      return !(std::fpclassify(event.target) == FP_ZERO);
                    });
    if (has_support) {
      return events;
    }
    if (batch + 1 == max_empty_batches) {
      if (sampler == NeuroJacSampler::FrozenVegas) {
        throw std::runtime_error(
            "MGraniitti::NeuroJacCollectTraining: frozen VEGAS proposal found "
            "no nonzero target in " +
            std::to_string(max_empty_batches) + " consecutive batches of " +
            std::to_string(batch_calls) +
            " samples, increase NUMERICS_NEUROJAC::vegas_ncall or "
            "vegas_rounds. Check generation and fiducial cuts.");
      }
      throw std::runtime_error(
          "MGraniitti::NeuroJacCollectTraining: all targets are zero in " +
          std::to_string(max_empty_batches) + " consecutive batches of " +
          std::to_string(batch_calls) +
          " samples. Check generation and fiducial cuts.");
    }
    const std::size_t empty_batches = batch + 1;
    const bool report_retry = (empty_batches & (empty_batches - 1)) == 0;
    if (!HILJAA && report_retry) {
      gra::aux::ClearProgress();
      std::cout << "NEUROJAC: "
                << (sampler == NeuroJacSampler::FrozenVegas
                        ? "frozen VEGAS found "
                        : "")
                << "no nonzero training target in " << empty_batches
                << " consecutive batch" << (empty_batches == 1 ? "" : "es")
                << " of " << batch_calls
                << " samples, retrying the same batch size" << std::endl;
    }
  }
  throw std::logic_error(
      "MGraniitti::NeuroJacCollectTraining: unreachable bootstrap state");
}

// Evaluate frozen NEUROJAC samples for integration or event generation
//
void MGraniitti::NeuroJacProduction(std::size_t thread_id, std::size_t calls,
                                    unsigned int requested_events,
                                    SamplingStage stage, const std::uint64_t first_sample) {
  MProcess *worker = pvec.at(thread_id);
  const neurojac::MFlowSampleBatch samples =
      neurojac.SampleBatch(worker->state.random, calls);
  std::vector<double> point(samples.dimension);
  for (std::size_t event = 0; event < samples.Size(); ++event) {
    if (stage == SamplingStage::Generation) {
      std::lock_guard<std::mutex> lock(state_mutex);
      if (stat.generated >= requested_events ||
          (max_weight_state.generation_overflow &&
           max_weight_param.OVERFLOW_ACTION !=
               OverflowAction::KeepAsWeighted)) {
        break;
      }
    }
    std::copy_n(samples.points.data() + event * samples.dimension,
                samples.dimension, point.data());
    const double log_inverse_density = -samples.log_densities[event];
    gra::MEventWeightState aux;
    aux.log_inverse_density = log_inverse_density;
    aux.qmetrics = stage == SamplingStage::Integration && worker->state.lts.process.QMETRICS;
    if (aux.qmetrics) { aux.sample_index = first_sample + event; }
    aux.adaptation_mode = false;
    const double raw_weight = worker->EventWeight(point, aux);
    double weight = 0.0;

    {
      std::lock_guard<std::mutex> lock(state_mutex);
      weight = ObserveSample(aux, raw_weight, log_inverse_density, stage);
      ObserveMaximumWeight(weight, stage);
    }

    if (stage == SamplingStage::Generation) {
      SaveEvent(worker, weight, GetGenerationMaxweight(), aux);
    }
  }
}

// Train, freeze and sample the bounded NEUROJAC importance proposal
//
void MGraniitti::SampleNeuroJac(unsigned int N) {
  local_tictoc = MTimer(true);
  atime = MTimer(true);

  if (N == 0) {
    GMODE = 0;
    InitMultiMemory();
    neurojac.ReadParameters(
        gra::ResolveModelDataFile(modelparam, "NUMERICS.json"));
    if (MC_VGRID_INPUT != "null") {
      ReadVGRID(MC_VGRID_INPUT);
      if (!neurojac.IsFrozen()) {
        throw std::logic_error(
            "MGraniitti::SampleNeuroJac: reloaded proposal is not frozen");
      }
      if (stat.integration_samples > 0) {
        PrintStatistics(N);
        return;
      }
      if (SamplesCrossSectionDuringGeneration()) {
        return;
      }
    } else {
      neurojac.Configure(proc->GetdLIPSDim(), proc->state.random.GetSeed());

      if (!HILJAA) {
        gra::aux::PrintBar("-");
        std::cout << rang::style::bold
                  << "NEUROJAC spline-flow integrator:" << rang::style::reset
                  << std::endl
                  << std::endl;
        PrintNeuroJacParameters();
        gra::aux::PrintBar("-");
      }

      std::size_t adaptation_evaluations = 0;
      if (neurojac.Config().vegas_init) {
        neurojac_vegas_param = vparam;
        neurojac_vegas_param.NCALL =
            static_cast<unsigned int>(neurojac.Config().vegas_ncall);
        neurojac_vegas_param.ROUNDS =
            static_cast<unsigned int>(neurojac.Config().vegas_rounds);
        neurojac_vegas_param.UNIFORM_MIX = neurojac.Config().uniform_mix;
        neurojac_vegas_param.AUTOMATIC_CONVERGENCE = false;
        vegas.SetDimension(proc->GetdLIPSDim());

        if (!HILJAA) {
          std::cout << std::endl
                    << rang::style::bold
                    << "NEUROJAC VEGAS presampling:" << rang::style::reset
                    << std::endl
                    << std::endl;
        }
        const std::uint64_t pilot_nonzero_calls =
            VEGAS(VEGASStage::Adaptation, neurojac_vegas_param.NCALL,
                  neurojac_vegas_param.ROUNDS, N, neurojac_vegas_param);
        adaptation_evaluations =
            neurojac.Config().vegas_ncall * neurojac.Config().vegas_rounds;
        if (pilot_nonzero_calls == 0) {
          throw std::runtime_error(
              "MGraniitti::SampleNeuroJac: VEGAS initialization found no "
              "nonzero target in any pilot round, increase "
              "NUMERICS_NEUROJAC::vegas_ncall or vegas_rounds. Check "
              "generation and fiducial cuts.");
        }
        neurojac.InitializeVegasGrid(vegas.QuantileEdges());

        std::vector<neurojac::MFlowTrainingEvent> initial_validation =
            NeuroJacCollectBatch(neurojac.ValidationEventCount(),
                                 NeuroJacSampler::FrozenVegas);
        adaptation_evaluations += initial_validation.size();
        neurojac.BufferValidationEvents(std::move(initial_validation));
        const neurojac::MFlowValidationReport baseline =
            neurojac.ValidateBaseline();
        if (!HILJAA) {
          std::cout << std::endl;
          std::cout << "NEUROJAC: "
                    << (neurojac.Config().vegas_spline_init
                            ? "VEGAS initialized spline baseline"
                            : "initial spline mixture baseline")
                    << " [val ESS/N = " << baseline.ess_fraction << ", support "
                    << baseline.nonzero_events << "/" << baseline.events << "]"
                    << std::endl;
        }
      }

      if (!HILJAA) {
        std::cout << std::endl
                  << rang::style::bold
                  << "NEUROJAC adaptation:" << rang::style::reset << std::endl
                  << std::endl;
      }

      for (std::size_t round = 0; round < neurojac.Config().rounds; ++round) {
        const NeuroJacSampler sampler =
            neurojac.Config().vegas_init &&
                    round < neurojac.Config().vegas_sampler_rounds
                ? NeuroJacSampler::FrozenVegas
                : NeuroJacSampler::Flow;
        if (!HILJAA &&
            (round == 0 || (neurojac.Config().vegas_init &&
                            round == neurojac.Config().vegas_sampler_rounds))) {
          if (round > 0) {
            std::cout << std::endl;
          }
          std::cout << "NEUROJAC: training sampler = "
                    << (sampler == NeuroJacSampler::FrozenVegas ? "frozen VEGAS"
                                                                : "flow")
                    << std::endl;
        }
        std::vector<neurojac::MFlowTrainingEvent> training =
            NeuroJacCollectTraining(adaptation_evaluations, sampler);
        const std::size_t validation_events = neurojac.ValidationEventCount();
        std::vector<neurojac::MFlowTrainingEvent> validation_buffer =
            NeuroJacCollectBatch(validation_events, sampler);
        adaptation_evaluations += validation_events;

        neurojac.ClearBuffer();
        neurojac.ClearValidationBuffer();
        neurojac.BufferEvents(std::move(training));
        neurojac.BufferValidationEvents(std::move(validation_buffer));
        const neurojac::MFlowTrainingReport report = neurojac.TrainRound(static_cast<std::size_t>(CORES));
        const neurojac::MFlowValidationReport validation =
            neurojac.ValidateRound();
        if (!HILJAA) {
          PrintNeuroJacTrainingStatus(report, validation,
                                      adaptation_evaluations);
          gra::aux::PrintProgress(
              (round + 1) / static_cast<double>(neurojac.Config().rounds));
        }
        if (validation.early_stop) {
          if (!HILJAA) {
            gra::aux::ClearProgress();
            std::cout
                << "NEUROJAC: validation patience reached, restoring round "
                << validation.best_round << std::endl;
          }
          break;
        }
      }
      if (!HILJAA) {
        gra::aux::ClearProgress();
        std::cout << std::endl;
      }

      neurojac.Freeze();

      // Save the trained flow without production integration
      if (IntegrationProposalOnly()) {
        SaveVGRID();
        return;
      }
    }

    if (!HILJAA) {
      std::cout << rang::style::bold
                << "NEUROJAC integration:" << rang::style::reset << std::endl
                << std::endl;
      PrintIntegrationSetup();
      std::cout << std::endl;
    }
    InitializeMaximumWeightEnvelope();
  } else {
    GMODE = 1;
    if (!neurojac.IsFrozen()) {
      throw std::logic_error("MGraniitti::SampleNeuroJac: integration must "
                             "train and freeze NEUROJAC before generation");
    }
  }

  const SamplingStage production_stage =
      N == 0 ? SamplingStage::Integration : SamplingStage::Generation;
  bool histogram_bounds_ready = production_stage == SamplingStage::Generation;
  while (true) {
    if (GMODE == 0) {
      stat.ResetIntegrationBatch();
    }
    std::size_t batch_calls = AutomaticIntegrationBatch();
    if (GMODE == 0) {
      if (stat.integration_samples >= integration_param.MAX_SAMPLES) {
        CheckIntegrationLimit();
        break;
      }
      const std::uint64_t remaining =
          integration_param.MAX_SAMPLES - stat.integration_samples;
      batch_calls = static_cast<std::size_t>(
          std::min<std::uint64_t>(batch_calls, remaining));
    }
    if (batch_calls == 0) {
      CheckIntegrationLimit();
      break;
    }
    std::vector<std::size_t> calls(static_cast<std::size_t>(CORES),
                                   batch_calls /
                                       static_cast<std::size_t>(CORES));
    for (std::size_t i = 0; i < batch_calls % static_cast<std::size_t>(CORES);
         ++i) {
      ++calls[i];
    }

    {
      std::lock_guard<std::mutex> lock(state_mutex);
      worker_exception = nullptr;
    }
    std::vector<std::thread> threads;
    std::uint64_t first_sample = stat.integration_samples;
    for (int thread_id = 0; thread_id < CORES; ++thread_id) {
      const auto begin = first_sample;
      first_sample += calls[thread_id];
      threads.emplace_back([&, thread_id, begin] {
        try {
          NeuroJacProduction(static_cast<std::size_t>(thread_id),
                             calls[thread_id], N, production_stage, begin);
        } catch (...) {
          std::lock_guard<std::mutex> lock(state_mutex);
          if (!worker_exception) {
            worker_exception = std::current_exception();
          }
        }
      });
    }
    for (std::thread &thread : threads) {
      thread.join();
    }
    if (worker_exception) {
      std::rethrow_exception(worker_exception);
    }
    if (GMODE == 1 && max_weight_state.generation_overflow &&
        max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted) {
      return;
    }

    // Derive common histogram bounds from the first integration batch
    if (!histogram_bounds_ready) {
      UnifyHistogramBounds();
      histogram_bounds_ready = true;
    }

    if (GMODE == 0) {
      stat.UpdateIntegrationChi2();
    }
    stat.CalculateCrossSection();
    if (GMODE == 0) {
      PrintStatus(static_cast<unsigned int>(std::min<std::uint64_t>(
                      stat.integration_samples,
                      std::numeric_limits<unsigned int>::max())),
                  N, local_tictoc, 10.0);
      if (stat.integration_samples >= integration_param.MIN_SAMPLES) {
        if (!stat.ValidCrossSection()) {
          throw std::invalid_argument(
              "NEUROJAC:: Invalid zero-support or non-finite integral "
              "estimate. Check generation and fiducial cuts.");
        }
        InitializeMaximumWeightEnvelope();
        if (IntegrationConverged()) {
          break;
        }
      }
      CheckIntegrationLimit();
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(std::min(
            1.0, stat.integration_samples /
                     static_cast<double>(integration_param.MIN_SAMPLES)));
        atime.Reset();
      }
    } else {
      PrintStatus(stat.generated, N, local_tictoc, 10.0);
      if (stat.generated >= N) {
        break;
      }
      if (atime.ElapsedSec() > 0.1) {
        gra::aux::PrintProgress(stat.generated / static_cast<double>(N));
        atime.Reset();
      }
    }
  }

  PrintStatus(stat.generated, N, local_tictoc, -1.0);

  // Capture integration performance before serializing the frozen state
  if (GMODE == 0 && std::fpclassify(stat.integration_runtime) == FP_ZERO) {
    stat.integration_runtime = global_tictoc.ElapsedSec();
  }

  // Persist both the frozen flow and all subsequent production samples
  if (!max_weight_state.generation_overflow ||
      max_weight_param.OVERFLOW_ACTION == OverflowAction::KeepAsWeighted) {
    SaveVGRID();
  }
  PrintStatistics(N);
}

// Save unweighted or weighted event
//
int MGraniitti::SaveEvent(MProcess *pr, double weight, double MAXWEIGHT,
                          const gra::MEventWeightState &aux) {
  if (!WEIGHTED && (!(weight >= 0.0) || !std::isfinite(weight))) {
    throw std::invalid_argument(
        "MAX_WEIGHT: unweighted generation requires finite nonnegative "
        "sampling weights");
  }
  if (!WEIGHTED && (!(MAXWEIGHT > 0.0) || !std::isfinite(MAXWEIGHT))) {
    throw std::logic_error(
        "MAX_WEIGHT: generation maxweight envelope must be finite and positive");
  }
  const bool overweight = !WEIGHTED && weight > MAXWEIGHT;
  const bool kept_overweight = overweight && max_weight_param.OVERFLOW_ACTION ==
                                                 OverflowAction::KeepAsWeighted;
  {
    std::lock_guard<std::mutex> lock(state_mutex);

    // Serialize the generation boundary before accounting for this trial
    if (stat.generated >= static_cast<unsigned int>(GetNumberOfEvents()) ||
        (max_weight_state.generation_overflow &&
         max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted)) {
      return 1;
    }
    stat.trials += 1; // This is one trial more

    // Account for every overweight trial without affecting weighted generation
    if (overweight) {
      ++stat.N_overflow;
      const bool first_overflow =
          max_weight_state.ObserveGenerationOverflow(weight);
      if (first_overflow &&
          max_weight_param.OVERFLOW_ACTION == OverflowAction::KeepAsWeighted &&
          !HILJAA) {
        std::cout << std::endl;
        gra::aux::PrintWarning(true);
        printf("MAX_WEIGHT: accepting events with weight > maxweight (increase INTEGRATOR.min_samples), "
               "all overflow events are stored as weighted "
               "(w/wmax = %0.3f)\n",
               weight / MAXWEIGHT);
      }
    }
  }

  if (overweight &&
      max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted) {
    return 3;
  }

  bool hit_in = false;
  if (!WEIGHTED) {
    hit_in = kept_overweight ||
             pr->state.random.U(0, 1) < std::max(weight / MAXWEIGHT, 0.0);
  }

  // All acceptance paths require complete physics and technical validity
  if (aux.Valid() && (hit_in || WEIGHTED || aux.forced_accept)) {
    // Create HepMC3 event (do not lock yet for speed)
    HepMC3::GenEvent evt(HepMC3::Units::GEV, HepMC3::Units::MM);

    // ** Construct event record **
    if (!pr->EventRecord(evt)) { // Event not ok!
      return 2;
    }

    // Add QED radiation only after a complete accepted Born record exists
    if (!pr->ApplyRadiation(evt)) {
      return 2;
    }

    // Event ok, continue >>

    // Serialize event numbering, shared statistics and output writers
    std::lock_guard<std::mutex> lock(state_mutex);

    // ** This is a multithreading race-condition treatment **
    if (stat.generated >= static_cast<unsigned int>(GetNumberOfEvents()) ||
        (max_weight_state.generation_overflow &&
         max_weight_param.OVERFLOW_ACTION != OverflowAction::KeepAsWeighted)) {
      return 1;
    }

    stat.CalculateCrossSection();

    // Set event number
    evt.set_event_number(stat.generated);
    const unsigned int accepted_events = stat.generated + 1;

    // Save ordinary weighted, unit weight or explicit overweight events
    const double HepMC3_weight =
        WEIGHTED ? weight : (kept_overweight ? weight / MAXWEIGHT : 1.0);
    evt.weights().push_back(
        HepMC3_weight); // add more weight attributes with .push_back()

    if (kept_overweight && FORMAT == "hepmc3") {
      evt.add_attribute(
          "maximum_weight_overflow",
          std::make_shared<HepMC3::DoubleAttribute>(weight / MAXWEIGHT));
    }

    if (WEIGHTED) {
      weighted_event_stats.Push(HepMC3_weight);
      if (FORMAT == "hepmc3") {
        evt.add_attribute("weighted_ess",
                          std::make_shared<HepMC3::DoubleAttribute>(
                              weighted_event_stats.EffectiveSampleSize()));
        evt.add_attribute("weighted_ess_fraction",
                          std::make_shared<HepMC3::DoubleAttribute>(
                              weighted_event_stats.EffectiveSampleFraction()));
      }
    }

    // ** Save cross section information (HepMC3 wants event by event)
    std::shared_ptr<HepMC3::GenCrossSection> xsobj =
        std::make_shared<HepMC3::GenCrossSection>();
    evt.add_attribute("GenCrossSection", xsobj);

    // Now add the value in picobarns [HepMC3 convention]
    if (FORMAT == "lhe") {
      xsobj->set_cross_section(output_xs_pb, output_xs_err_pb);
    } else if (xsforced > 0) {
      xsobj->set_cross_section(xsforced * OUTPUT_XS_SCALE,
                               0); // external fixed one
    } else {
      xsobj->set_cross_section(stat.sigma * OUTPUT_XS_SCALE,
                               stat.sigma_err * OUTPUT_XS_SCALE);
    }

    // Add the number of generated and attempted
    xsobj->set_accepted_events(accepted_events);
    xsobj->set_attempted_events(
        stat.trials); // => allows computing xs = sum(weights) / attempted

    if (FORMAT == "hepmc3") {
      outputHepMC3->write_event(evt);
    } else if (FORMAT == "hepmc2") {
      outputHepMC2->write_event(evt);
    } else if (FORMAT == "hepevt") {
      outputHEPEVT->write_event(evt);
    } else if (FORMAT == "lhe") {
      outputLHE->WriteEvent(evt);
    } else {
      throw std::invalid_argument(
          "MGraniitti::SaveEvent: Unknown output FORMAT " + FORMAT);
    }

    CheckFileOutput();

    // LAST STEP
    stat.generated = accepted_events; // +1 event generated

    // Cleanup up the event from memory
    evt.clear();

    return 0;
  } else {
    return 1;
  }
}

// Report event stream failures before accepting an output event
void MGraniitti::CheckFileOutput() const {
  if ((outputHepMC3 != nullptr && outputHepMC3->failed()) ||
      (outputHepMC2 != nullptr && outputHepMC2->failed()) ||
      (outputHEPEVT != nullptr && outputHEPEVT->failed())) {
    throw std::runtime_error("MGraniitti: failed event output " + FULL_OUTPUT_STR);
  }
}

// Close every internally owned event writer
void MGraniitti::CloseFileOutput() {
  if (!owns_file_output) {
    throw std::runtime_error(
        "MAX_WEIGHT: overflow handling requires an internally owned "
        "output writer");
  }
  if (outputHepMC3 != nullptr) {
    outputHepMC3->close();
  }
  if (outputHepMC2 != nullptr) {
    outputHepMC2->close();
  }
  if (outputHEPEVT != nullptr) {
    outputHEPEVT->close();
  }
  CheckFileOutput();
  if (outputLHE != nullptr) { outputLHE->Close(); }
  outputHepMC3.reset();
  outputHepMC2.reset();
  outputHEPEVT.reset();
  outputLHE.reset();
  runinfo.reset();
  owns_file_output = false;
}

// Rename one invalid partial generation output for recovery
void MGraniitti::PreserveFailedOutput() {
  const std::filesystem::path source(FULL_OUTPUT_STR);
  if (!std::filesystem::exists(source)) {
    return;
  }
  std::filesystem::path destination(FULL_OUTPUT_STR + "._overflow");
  std::size_t suffix = 1;
  while (std::filesystem::exists(destination)) {
    destination = FULL_OUTPUT_STR + "._overflow." + std::to_string(suffix);
    ++suffix;
  }
  std::filesystem::rename(source, destination);
  if (!HILJAA) {
    std::cout << "MAX_WEIGHT: preserved invalid partial output as "
              << destination.string() << std::endl;
  }
}

// Initialize immediately before event generation
//
void MGraniitti::InitFileOutput() {
  if (NEVENTS > 0) {
    if (OUTPUT == "") { // OUTPUT must be set
      throw std::invalid_argument(
          "MGraniitti::InitFileOutput: OUTPUT filename not set!");
    }

    FULL_OUTPUT_STR = OutputPath(OUTPUT, "output", FORMAT);

    const bool generation_xs = SamplesCrossSectionDuringGeneration();
    const bool integrated_xs =
        stat.integration_samples > 0 && stat.ValidCrossSection();
    if (FORMAT == "lhe" && WEIGHTED && !integrated_xs) {
      throw std::invalid_argument(
          "MGraniitti::InitFileOutput: weighted LHE output requires a "
          "finite integrated cross section");
    }

    // Resolve fixed, sampled or generation-time cross section normalization
    if (xsforced > 0.0 || stat.ValidCrossSection()) {
      output_xs_pb =
          (xsforced > 0.0 ? xsforced : stat.sigma) * OUTPUT_XS_SCALE;
      output_xs_err_pb =
          xsforced > 0.0 ? 0.0 : stat.sigma_err * OUTPUT_XS_SCALE;
      if (!(output_xs_pb > 0.0) || !std::isfinite(output_xs_pb) ||
          !(output_xs_err_pb >= 0.0) || !std::isfinite(output_xs_err_pb)) {
        throw std::invalid_argument(
            "MGraniitti::InitFileOutput: output cross section is invalid");
      }
    } else if (generation_xs) {
      output_xs_pb = 0.0;
      output_xs_err_pb = 0.0;
    } else {
      throw std::invalid_argument(
          "MGraniitti::InitFileOutput: output cross section is invalid");
    }

    output_lhe_weight_scale = 1.0;
    if (FORMAT == "lhe" && WEIGHTED) {
      const double valid_fraction =
          stat.all_ok / static_cast<double>(stat.integration_samples);
      if (!(valid_fraction > 0.0) || valid_fraction > 1.0 ||
          !std::isfinite(valid_fraction)) {
        throw std::invalid_argument(
            "MGraniitti::InitFileOutput: weighted LHE output has no valid "
            "integration support");
      }
      output_lhe_weight_scale =
          output_xs_pb * valid_fraction / stat.sigma;
      if (!(output_lhe_weight_scale > 0.0) ||
          !std::isfinite(output_lhe_weight_scale)) {
        throw std::invalid_argument(
            "MGraniitti::InitFileOutput: weighted LHE scale is invalid");
      }
    }

    // --------------------------------------------------------------
    // Generator info
    runinfo = std::make_shared<HepMC3::GenRunInfo>();

    struct HepMC3::GenRunInfo::ToolInfo generator = {
        std::string("GRANIITTI (" + modelparam + ")"),
        std::to_string(aux::GetVersion()).substr(0, 5),
        std::string("Generator")};
    runinfo->tools().push_back(generator);

    struct HepMC3::GenRunInfo::ToolInfo config = {FULL_PROCESS_STR, "1.0",
                                                  std::string("Steering card")};
    runinfo->tools().push_back(config);
    if (proc->state.nuclear_final.has_value()) {
      runinfo->tools().push_back(proc->state.nuclear_final->Tool());
    }

    // ** Cross section is added also event by event later **
    // Now add the value in picobarns [HepMC3 convention]
    runinfo->add_attribute(
        "xs", std::make_shared<HepMC3::FloatAttribute>(output_xs_pb));
    runinfo->add_attribute(
        "xs_err", std::make_shared<HepMC3::FloatAttribute>(output_xs_err_pb));

    // --------------------------------------------------------------

    if (FORMAT == "hepmc3" && outputHepMC3 == nullptr) {
      outputHepMC3 =
          std::make_shared<HepMC3::WriterAscii>(FULL_OUTPUT_STR, runinfo);
      owns_file_output = true;
    } else if (FORMAT == "hepmc2" && outputHepMC2 == nullptr) {
      outputHepMC2 =
          std::make_shared<HepMC3::WriterAsciiHepMC2>(FULL_OUTPUT_STR, runinfo);
      owns_file_output = true;
    } else if (FORMAT == "hepevt" && outputHEPEVT == nullptr) {
      outputHEPEVT = std::make_shared<HepMC3::WriterHEPEVT>(FULL_OUTPUT_STR);
      owns_file_output = true;
    } else if (FORMAT == "lhe" && outputLHE == nullptr) {
      gra::MLHERunConfig config;
      if (WEIGHTED) {
        config.weight_type = gra::MLHEWeightType::Weighted;
        config.weight_scale = output_lhe_weight_scale;
      }
      outputLHE =
          std::make_shared<gra::MLHEWriter>(FULL_OUTPUT_STR, config);
      owns_file_output = true;
    }
    CheckFileOutput();
  }
}

// Print one NEUROJAC training update without cross section statistics
//
void MGraniitti::PrintNeuroJacTrainingStatus(
    const neurojac::MFlowTrainingReport &training,
    const neurojac::MFlowValidationReport &validation,
    std::size_t evaluations) {
  gra::aux::ClearProgress();

  double peak_use = 0.0;
  double resident_use = 0.0;
  constexpr double MB = 1024.0 * 1024.0;
  gra::aux::GetProcessMemory(peak_use, resident_use);
  resident_use /= MB;

  const double global_lap = global_tictoc.ElapsedSec();
  const double frequency = global_lap > 0.0 ? evaluations / global_lap : 0.0;
  printf("[%0.1f MB] loss: %9.3E, lr: %9.3E "
         "[train ESS/N = %6.4f, val ESS/N = %6.4f, "
         "saved ESS/N = %6.4f @ %zu, support %zu/%zu], %4.1f min ~ %0.1E Hz \n",
         resident_use, training.loss, training.learning_rate,
         training.collection_ess_fraction, validation.ess_fraction,
         validation.best_ess_fraction, validation.best_round,
         training.nonzero_events, training.buffered_events, global_lap / 60.0,
         frequency);
  if (training.rejected_epochs > 0) {
    printf("NEUROJAC: restored the previous proposal after %zu unstable optimizer epochs\n",
           training.rejected_epochs);
  }
  if (training.component_weights.size() > 1) {
    printf("flow weights = [");
    for (const auto &component : indices(training.component_weights)) {
      printf("%s%.4f", component == 0 ? "" : ", ",
             training.component_weights[component]);
    }
    printf("] (unconstrained newest = %.4f)\n", training.boost_weight);
  }
}

// Print one VEGAS adaptation update without cross section statistics
//
void MGraniitti::PrintVegasAdaptationStatus(
    const VEGASAdaptationReport &report) {
  gra::aux::ClearProgress();

  double peak_use = 0.0;
  double resident_use = 0.0;
  constexpr double MB = 1024.0 * 1024.0;
  gra::aux::GetProcessMemory(peak_use, resident_use);
  resident_use /= MB;

  const double lap = local_tictoc.ElapsedSec();
  const double frequency = lap > 0.0 ? report.evaluations / lap : 0.0;
  printf("[%0.1f MB] round: %zu/%zu [ESS/N = %6.4f, support %zu/%zu], "
         "%4.1f min ~ %0.1E Hz \n",
         resident_use, report.round, report.round_limit, report.ess_fraction,
         report.support, report.calls, lap / 60.0, frequency);
}

// Intermediate statistics
//
void MGraniitti::PrintStatus(unsigned int events, unsigned int N,
                             MTimer &tictoc, double timercut) {
  if (tictoc.ElapsedSec() > timercut) {
    tictoc.Reset();
    gra::aux::ClearProgress();

    double peak_use = 0.0;
    double resident_use = 0.0;
    const double MB = 1024 * 1024;
    const double GB = MB * 1024;
    gra::aux::GetProcessMemory(peak_use, resident_use);
    peak_use /= MB;
    resident_use /= MB;

    double sigma = 0.0;
    double sigma_err = 0.0;
    double chi2 = 0.0;
    double ess_fraction = 0.0;
    {
      std::lock_guard<std::mutex> lock(state_mutex);
      sigma = stat.sigma;
      sigma_err = stat.sigma_err;
      chi2 = stat.chi2;
      ess_fraction = stat.SamplingEffectiveFraction();
    }
    const double relative_error = !(std::fpclassify(sigma) == FP_ZERO)
                                      ? std::abs(sigma_err / sigma)
                                      : std::numeric_limits<double>::infinity();

    if (GMODE == 0) {
      const double global_lap = global_tictoc.ElapsedSec();
      printf("[%0.1f MB] xs: %9.3E, err: %7.5f "
             "[chi2/dof = %5.2f, ESS/N = %6.4f], %4.1f "
             "min ~ %0.1E Hz \n",
             resident_use, sigma, relative_error, chi2, ess_fraction,
             global_lap / 60.0, events / global_lap);
    }
    if (GMODE == 1) {
      const double global_lap = global_tictoc.ElapsedSec() - time_t0;
      double outputfilesize = gra::aux::GetFileSize(FULL_OUTPUT_STR) / GB;

      printf("[%0.1f MB/%0.2f GB] E: %9d, xs: %9.3E, err: %7.5f, %0.1f/%0.1f "
             "min ~ "
             "%0.1E Hz \n",
             resident_use, outputfilesize, events, sigma, relative_error,
             global_lap / 60.0,
             (N - events) * global_lap / (double)events / 60.0,
             events / global_lap);
    }
  }
}

// Print frozen-proposal performance, weight and unweighting diagnostics
void MGraniitti::PrintIntegrationDiagnostics(double runtime) const {
  const double frequency =
      runtime > 0.0 ? stat.integration_samples / runtime : 0.0;

  printf("%-29s  %10s  %s\n", "Performance", "Value", "Unit");
  printf("-----------------------------  ----------  ----\n");
  printf("%-29s  %10llu\n", "Integration samples",
         static_cast<unsigned long long>(stat.integration_samples));
  printf("%-29s  %10.2f  %s\n", "Integration runtime", runtime, "s");
  printf("%-29s  %10.2E  %s\n", "Integrand sampling frequency", frequency,
         "Hz");

  printf("\n");
  printf("%-30s  %10s  %10s  %12s  %10s  %10s\n", "Sampling", "Mean",
         "Std / mean", "Min positive", "Maximum", "ESS / N");
  printf("------------------------------  ----------  ----------  ------------"
         "  ----------  ----------\n");
  printf("%-30s  %10.3LE  %10.3LE  %12.3LE  %10.3LE  %10s\n", "Integrand f(x)",
         stat.IntegrandWeightMean(), stat.IntegrandWeightRelativeStd(),
         stat.IntegrandWeightMinimumPositive(), stat.IntegrandWeightMaximum(),
         "-");
  printf("%-30s  %10.3LE  %10.3LE  %12.3LE  %10.3LE  %10s\n",
         "Proposal density q(x)", stat.ProposalDensityMean(),
         stat.ProposalDensityRelativeStd(),
         stat.ProposalDensityMinimumPositive(), stat.ProposalDensityMaximum(),
         "-");
  printf("%-30s  %10.3LE  %10.3LE  %12.3LE  %10.3LE  %10.3E\n",
         "Importance weight f(x) / q(x)", stat.SamplingWeightMean(),
         stat.SamplingWeightRelativeStd(), stat.SamplingWeightMinimumPositive(),
         stat.SamplingWeightMaximum(), stat.SamplingEffectiveFraction());

  if (WEIGHTED || !max_weight_state.initialized) {
    return;
  }

  const double maximum = GetMaxweight();
  const long double card_maximum = ConservativeMaximumForPrint(maximum);
  printf("\n");
  printf("%-22s    %13s  %s\n", "Unweighting", "Value", "Definition");
  printf(
      "----------------------    -------------  --------------------------\n");
  printf("%-22s    %13.3LE  %s\n", "Maximum weight", card_maximum,
         "INTEGRATOR::max_w_value");
  printf("%-22s    %13.3E  %s\n", "Estimated efficiency",
         stat.EstimatedUnweightingEfficiency(static_cast<double>(card_maximum)),
         "mean / maximum weight");
  if (max_weight_param.MODE == MaxWeightMode::Estimate) {
    const std::string validation =
        std::to_string(max_weight_state.validation_samples) + " / " +
        std::to_string(RequiredMaximumWeightValidation());
    printf("%-22s    %13s  %s\n", "Validation samples", validation.c_str(),
           "validated / required");
  }
}

// Print sequential event flow and independent failure diagnostics
void MGraniitti::PrintEventFlowStatistics() const {
  printf("%-26s  %12s  %s\n", "Event flow", "Passing rate", "Denominator");
  printf("--------------------------  ------------  ---------------------\n");
  printf("%-26s  %12.3E  %s\n", "Kinematics",
         Stats::ConditionalRate(stat.kinematics_ok, stat.evaluations),
         "all evaluations");
  printf("%-26s  %12.3E  %s\n", "Fiducial cuts",
         Stats::ConditionalRate(stat.fidcuts_ok, stat.kinematics_ok),
         "passed kinematics");
  printf("%-26s  %12.3E  %s\n", "Veto cuts",
         Stats::ConditionalRate(stat.vetocuts_ok, stat.fidcuts_ok),
         "passed fiducial cuts");
  printf("%-26s  %12.3E  %s\n", "Amplitude",
         Stats::ConditionalRate(stat.amplitude_ok, stat.amplitude_evaluations),
         "amplitude evaluations");
  printf("--------------------------  ------------  ---------------------\n");
  printf("%-26s  %12.3E  %s\n", "Total",
         Stats::ConditionalRate(stat.all_ok, stat.evaluations),
         "all evaluations");

  printf("\n");
  printf("%-26s  %12s  %s\n", "Failures", "Failure rate", "Denominator");
  printf("--------------------------  ------------  ---------------------\n");
  printf("%-26s  %12.3E  %s\n", "Amplitude",
         Stats::ConditionalRate(stat.amplitude_failures,
                                stat.amplitude_evaluations),
         "amplitude evaluations");
  printf("%-26s  %12.3E  %s\n", "Technical",
         Stats::ConditionalRate(stat.technical_failures, stat.evaluations),
         "all evaluations");
}

// Final statistics
//
void MGraniitti::PrintStatistics(unsigned int N) {
  gra::aux::ClearProgress(); // Clear progressbar

  if (GMODE == 0) {
    time_t0 = global_tictoc.ElapsedSec();
    if (std::fpclassify(stat.integration_runtime) == FP_ZERO) {
      stat.integration_runtime = time_t0;
    }
    gra::aux::PrintBar("=");
    std::cout << rang::style::bold
              << "MC cross section summary:" << rang::style::reset << std::endl
              << std::endl;

    if (proc->GetISOLATE()) {
      std::cout
          << rang::fg::red
          << "NOTE: Central leg phase space isolation tag &> in use!"
          << std::endl
          << "      The decay tree is used only as a kinematic proposal for "
             "acceptance:"
          << " phase space with Breit-Wigner mass sampling" << std::endl
          << "      The recursive decay phase-space volume or couplings are not"
          << " included in the reported cross section, neither decay matrix"
          << " element / decay spin algebra; apply BRs instead" << std::endl
          << "      Resonance hel_decay.g_decay is fixed to exactly 1 + 0i in "
             "&> mode"
          << std::endl
          << "      Generated decay products still enter fiducial and veto "
             "cuts,"
          << " which impacts the printed cross-section";
      std::cout << rang::fg::reset;
      std::cout << std::endl << std::endl;
    }

    // Check if we have cascaded phase space turned on in the x-section
    // calculation
    const decay::GeneratedPhaseSpaceSummary phase_space =
        decay::GeneratedPhaseSpace(
            proc->state.lts.decaytree,
            proc->state.lts.decay_symmetry_proposal_active);
    const int Nf = phase_space.stable_leaves;
    unsigned int N_leg = std::max((int)proc->state.lts.decaytree.size(), Nf) +
                         2; // +2 forward legs

    // Special cases
    if (proc->GetISOLATE()) {
      N_leg = 3;
    } // ISOLATEd 2->3 process with <F> phase space
    if (proc->GetdLIPSDim() == 2) {
      N_leg = 2;
    } // <P> and <Q> class

    std::cout << rang::fg::yellow;
    printf("Fiducial cross section:    [%0.3E +- %0.3E] barn\n", stat.sigma,
           stat.sigma_err);
    std::cout << rang::fg::reset;
    const decay::BranchingRatioSummary branching =
        decay::IsolatedBranchingRatioProduct(proc->state.lts,
                                             proc->GetISOLATE());
    if (proc->GetISOLATE() && branching.factors > 0) {
      std::cout << rang::fg::yellow;
      printf("Fiducial cross section x BRs: [%0.3E +- %0.3E] barn\n",
             stat.sigma * branching.product,
             stat.sigma_err * branching.product);
      std::cout << rang::fg::reset;
      printf("Decay BR product (%u factors):     %0.3E\n", branching.factors,
             branching.product);
    }
    printf("Sampling uncertainty:     %0.3f %%\n",
           100 * stat.sigma_err / stat.sigma);
    printf("Reduced integration chi2: %0.3f\n", stat.chi2);
    std::cout << std::endl;

    if (proc->state.lts.DW_sum.Integral() > 0) { // Recursive phase space on

      std::cout << std::endl;
      const double PSvolume = proc->state.lts.DW_sum.Integral();
      // const double PSvolume_error       =
      // proc->state.lts.DW_sum.IntegralError();

      const double PSvolume_exact = proc->state.lts.DW_sum_exact.Integral();
      const double PSvolume_exact_error =
          proc->state.lts.DW_sum_exact.IntegralError();

      // Get sum of decay daughter masses
      double MSUM = 0.0;
      for (const auto &i : indices(proc->state.lts.decaytree)) {
        MSUM += proc->state.lts.decaytree[i].m_offshell;
      }
      // only 2-body exact or massless case
      if (proc->state.lts.decaytree.size() == 2 || MSUM < 1e-6) {
        printf("Analytic phase space volume:      %0.3E +- %0.3E \n",
               PSvolume_exact, PSvolume_exact_error);
        printf("RATIO: MC/analytic:               %0.6f \n",
               PSvolume / PSvolume_exact);
      }
      std::cout << std::endl;

      // Print out phase space weight
      printf("{1->%lu LIPS}:                      %0.3E +- %0.3E  ",
             proc->state.lts.decaytree.size(),
             proc->state.lts.DW_sum.Integral(),
             proc->state.lts.DW_sum.IntegralError());

      if (proc->GetISOLATE()) {
        std::cout << "[" << rang::fg::red << "INACTIVE " << rang::fg::reset
                  << "part of cross-section integral]" << std::endl;
      } else {
        std::cout << "[" << rang::fg::green << "ACTIVE " << rang::fg::reset
                  << "part of integral / (2PI)]" << std::endl;
      }
      // Recursion relation based (phase space factorization):
      // d^N PS(s; p_1, p2, ...p_N)
      // = 1/(2*PI) * d^3 PS(s; p1,p2,pX) d^{N-2} PS(M^2; pX,p3,p4,..,pN) dM^2
      //
      // Decaywidth = 1/(2M S) \int dPS |M_decay|^2, where M
      // = mother mass, S = final state symmetry factor

      if (!proc->GetISOLATE() && proc->GetdLIPSDim() != 2 &&
          proc->ProcPtr.UsesJWHelicityAlgebra()) { // Special cases
        printf("\n\n{2->3 cross section ~=~ [2->%u / (1->%lu LIPS) x 2PI]}:    "
               "   %0.3E \n",
               N_leg, proc->state.lts.decaytree.size(),
               stat.sigma / proc->state.lts.DW_sum.Integral() * (2 * PI));
        std::cout << std::endl;
        printf("** Remember to use &> operator instead of -> for phase space "
               "isolation ** \n\n");
      }
    }

    // Print recursively
    for (const auto &i : indices(proc->state.lts.decaytree)) {
      decay::PrintIntegratedPhaseSpace(proc->state.lts.decaytree[i],
                                       proc->state.lts.PS_active);
      if (proc->state.lts.decaytree[i].legs.size() != 0) {
        std::cout << std::endl;
      }
    }
    gra::aux::PrintBar("-");

    if (proc->state.lts.process.QMETRICS) {
      QMetrics().Print(std::cout);
      gra::aux::PrintBar("-");
    }

    std::cout << std::endl
              << rang::style::bold
              << "Integration statistics:" << rang::style::reset << std::endl
              << std::endl;
    PrintIntegrationDiagnostics(stat.integration_runtime);
    printf("\n");
    PrintEventFlowStatistics();
    printf("\n");

    std::cout << std::endl;
    printf("** All values include phase space generation and "
           "fiducial (+ veto) cuts ** \n");
    gra::aux::PrintBar("=");
  }
  if (GMODE == 1) {
    const double lap = global_tictoc.ElapsedSec() - time_t0;

    gra::aux::PrintBar("=");
    if (WEIGHTED) {
      gra::aux::PrintNotice();
      std::cout << rang::fg::red
                << "You did WEIGHTED event generation:" << std::endl
                << std::endl
                << std::endl;
    } else {
      std::cout << rang::fg::green
                << "You did UNWEIGHTED (acceptance-rejection) event generation:"
                << std::endl
                << std::endl
                << std::endl;
    }
    std::cout << rang::style::reset;

    if (WEIGHTED) {
      const double ESS = weighted_event_stats.EffectiveSampleSize();
      printf("Effective Sample Size:    %0.4g [(sum_i w_i)^2 / sum_i w_i^2] \n",
             ESS);
      printf("Relative ESS:             %0.4g (%0.4g / %d) \n",
             weighted_event_stats.EffectiveSampleFraction(), ESS, N);
      std::cout << std::endl;
    }

    printf("Generation efficiency:    %0.3g (%d / %0.0f) \n",
           N / (double)stat.trials, N, stat.trials);
    printf("Weight overflow:          %0.3E (%llu / %0.0f) \n",
           stat.N_overflow / static_cast<double>(stat.trials),
           static_cast<unsigned long long>(stat.N_overflow), stat.trials);
    if (!WEIGHTED && stat.N_overflow > 0 &&
        max_weight_param.OVERFLOW_ACTION == OverflowAction::KeepAsWeighted) {
      std::cout << rang::fg::yellow
                << "WARNING: output contains explicitly weighted maximum "
                   "weight overflow events"
                << rang::fg::reset << std::endl;
      printf("Largest overflow weight ratio: %0.3E \n",
             max_weight_state.largest_overflow / max_weight_state.envelope);
    }
    printf("Generation runtime:       %0.2f sec \n", lap);
    printf("Generation frequency:     %0.2E Hz \n", N / lap);

    double outputfilesize =
        gra::aux::GetFileSize(FULL_OUTPUT_STR) / (1024.0 * 1024.0 * 1024.0);
    std::cout << std::endl;

    printf("Outputfile size:          %0.3f GB [%s] \n", outputfilesize,
           FULL_OUTPUT_STR.c_str());
    gra::aux::PrintBar("=");
    std::cout << std::endl;
  }
}

// Construct terminal input parameters
//
void MGraniitti::ConstructTerminal(cxxopts::Options &options) const {
  options.add_options("GENERIC")("debug", "Print every nonfatal sampling failure");
  options.add_options("")("set", "Override JSON card entry <[card:]path=json>",
                          cxxopts::value<std::string>());

  options.add_options("GENERIC")("o,OUTPUT", "Output name or path      <string>",
                                 cxxopts::value<std::string>())(
      "f,FORMAT", "Output format           <hepmc3|hepmc2|hepevt|lhe>",
      cxxopts::value<std::string>())("c,CORES",
                                     "Number of CPU threads   <integer>",
                                     cxxopts::value<unsigned int>())(
      "n,NEVENTS", "Events (-1: proposal only) <integer32>",
      cxxopts::value<int>())("g,INTEGRATOR",
                             "Integrator              <VEGAS|FLAT|NEUROJAC>",
                             cxxopts::value<std::string>())(
      "w,WEIGHTED", "Weighted events         <true|false|1|0>",
      cxxopts::value<std::string>())("m,MODELPARAM",
                                     "Model tune              <string>",
                                     cxxopts::value<std::string>())(
      "h,HIST", "Histogramming           <0|1|2>",
      cxxopts::value<unsigned int>())("r,RNDSEED",
                                      "Random seed             <integer32>",
                                      cxxopts::value<unsigned int>());

  options.add_options("SCATTERING")("p,PROCESS",
                                    "Process                 <string>",
                                    cxxopts::value<std::string>())(
      "e,ENERGY", "CMS energy per nucleon pair for ions <double>",
      cxxopts::value<double>())(
      "l,LOOPSCREEN", "Soft survival screening <true|false|1|0>",
      cxxopts::value<std::string>())("s,NSTARS",
                                     "Excite protons          <0|1|2>",
                                     cxxopts::value<unsigned int>())(
      "q,LHAPDF", "Set LHAPDF              <string>",
      cxxopts::value<std::string>());
}

// Override json object parameters from the command line object r
//
void MGraniitti::ProcessTerminal(json &j, cxxopts::ParseResult const &r) const {
  // GENERIC

  if (r.count("n")) {
    j.at("GENERIC").at("NEVENTS") = r["n"].as<int>();
  }
  if (r.count("o")) {
    j.at("GENERIC").at("OUTPUT") = r["o"].as<std::string>();
  }
  if (r.count("f")) {
    j.at("GENERIC").at("FORMAT") = r["f"].as<std::string>();
  }
  if (r.count("m")) {
    j.at("GENERIC").at("MODELPARAM") = r["m"].as<std::string>();
  }
  if (r.count("g")) {
    j.at("GENERIC").at("INTEGRATOR") = r["g"].as<std::string>();
  }
  if (r.count("w")) {
    j.at("GENERIC").at("WEIGHTED") =
        aux::ParseBool(r["w"].as<std::string>(), "--WEIGHTED");
  }
  if (r.count("c")) {
    j.at("GENERIC").at("CORES") = r["c"].as<unsigned int>();
  }
  if (r.count("h")) {
    j.at("GENERIC").at("HIST") = r["h"].as<unsigned int>();
  }
  if (r.count("r")) {
    j.at("GENERIC").at("RNDSEED") = r["r"].as<unsigned int>();
  }

  // SCATTERING

  if (r.count("p")) {
    j.at("SCATTERING").at("PROCESS") = r["p"].as<std::string>();
  }
  if (r.count("l")) {
    j.at("SCATTERING").at("LOOPSCREEN") =
        aux::ParseBool(r["l"].as<std::string>(), "--LOOPSCREEN");
  }
  if (r.count("e")) {
    j.at("SCATTERING").at("ENERGY") =
        std::vector<double>(2, r["e"].as<double>() / 2.0);
  }
  if (r.count("s")) {
    j.at("SCATTERING").at("NSTARS") = r["s"].as<unsigned int>();
  }
  if (r.count("q")) {
    j.at("SCATTERING").at("LHAPDF") = r["q"].as<std::string>();
  }
}

} // namespace gra
