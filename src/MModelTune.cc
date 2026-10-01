// GRANIITTI model tune (immutable) container
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/MModelTune.h"

// C++
#include <array>
#include <cmath>
#include <limits>
#include <filesystem>
#include <stdexcept>
#include <utility>

// Own
#include "Graniitti/Tech/MAux.h"

namespace gra {

namespace {

// Compute one canonical model card path
std::string CanonicalCardPath(const std::string &path) {
  try {
    return std::filesystem::weakly_canonical(path).string();
  } catch (const std::filesystem::filesystem_error &error) {
    throw std::invalid_argument("MModelTune: cannot resolve " + path + ": " + error.what());
  }
}

// Parse one immutable model card document
nlohmann::json ReadCard(const std::string &path) {
  try {
    const nlohmann::json document = nlohmann::json::parse(aux::GetInputData(path));
    if (!document.is_object()) { throw std::invalid_argument("root must be an object"); }
    return document;
  } catch (const std::exception &error) {
    throw std::invalid_argument("MModelTune: error reading " + path + ": " + error.what());
  }
}

// Read and validate the program wide numerical block
MGlobalNumerics ReadGlobalNumerics(const nlohmann::json &document) {
  MGlobalNumerics global{};
  const auto     &block = document.at("NUMERICS_GLOBAL");
  global.coupling_min   = block.at("coupling_min").get<double>();
  if (!std::isfinite(global.coupling_min) || global.coupling_min < 0.0) {
    throw std::invalid_argument("NUMERICS_GLOBAL.coupling_min must be finite and nonnegative");
  }
  return global;
}

// Read and validate charged Standard Model inputs before sampling
SMParam ReadSM(const nlohmann::json &document) {
  const auto &block = document.at("PARAM_SM");
  const auto &mass = block.at("mass");
  // Compute one positive finite physical input
  const auto positive = [](const nlohmann::json &values, const std::string &key) {
    const double value = values.at(key).get<double>();
    if (!std::isfinite(value) || value <= 0.0) {
      throw std::invalid_argument("PARAM_SM." + key + " must be finite and positive");
    }
    return value;
  };
  return {positive(mass, "e"), positive(mass, "mu"), positive(mass, "tau"),
          positive(mass, "u"), positive(mass, "d"), positive(mass, "s"),
          positive(mass, "c"), positive(mass, "b"), positive(mass, "t"),
          positive(mass, "w"), positive(block, "alpha_em_inv")};
}

}  // namespace

// Construct one fully parsed immutable tune
MModelTune::MModelTune(std::string general_file, std::string numerics_file, nlohmann::json general,
                       nlohmann::json numerics, SoftModelPtr soft, MGlobalNumerics global, SMParam sm, form::ParamStore structure,
                       form::FlatParam flat,
                       std::map<std::string, nlohmann::json> continuum)
    : general_file_(std::move(general_file)),
      numerics_file_(std::move(numerics_file)),
      general_(std::move(general)),
      numerics_(std::move(numerics)),
      soft_(std::move(soft)),
      global_(global),
      sm_(sm),
      structure_(std::move(structure)),
      flat_(flat),
      continuum_(std::move(continuum)) {}

// Select another soft model without changing the production tune
std::shared_ptr<const MModelTune> MModelTune::WithSoft(const std::string &name) const {
  auto tune = std::shared_ptr<MModelTune>(new MModelTune(*this));
  tune->general_["PARAM_SOFT"]["active_model"] = name;
  tune->soft_ = SoftModel::LoadFromJson(general_file_, tune->general_.dump(), numerics_file_, numerics_.dump());
  return tune;
}

// Load and validate one complete model tune from its GENERAL file
std::shared_ptr<const MModelTune> MModelTune::Load(const std::string &general_file) {
  const std::string           canonical_general = CanonicalCardPath(general_file);
  const std::filesystem::path general_path(canonical_general);
  const std::string numerics_file = CanonicalCardPath((general_path.parent_path() / "NUMERICS.json").string());
  nlohmann::json    general       = ReadCard(canonical_general);
  nlohmann::json    numerics      = ReadCard(numerics_file);
  for (const auto &[name, minimum] : std::array<std::pair<const char *, int>, 2>{{{"metric_groups", 2}, {"metric_max_dimension", 1}}}) {
    const auto &value = numerics.at("NUMERICS_SPIN").at(name);
    if (!value.is_number_integer() || value.get<long double>() < minimum ||
        value.get<long double>() > std::numeric_limits<int>::max()) {
      throw std::invalid_argument(std::string("NUMERICS_SPIN.") + name + " must be a supported positive integer");
    }
  }
  SoftModelPtr      soft   = SoftModel::LoadFromJson(canonical_general, general.dump(), numerics_file, numerics.dump());
  MGlobalNumerics   global = ReadGlobalNumerics(numerics);
  const SMParam    sm = ReadSM(general);
  form::ParamStore  structure = form::ParamStore::Read(canonical_general, general.dump());
  const form::FlatParam flat = form::FlatParam::Read(canonical_general, general.dump());
  std::map<std::string, nlohmann::json> continuum;
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    continuum.emplace(model, ReadCard((general_path.parent_path() / ("CON_" + model + ".json")).string()));
  }
  return std::shared_ptr<const MModelTune>(new MModelTune(canonical_general, numerics_file, std::move(general),
                                                          std::move(numerics), std::move(soft), global,
                                                          sm, std::move(structure), flat, std::move(continuum)));
}

// Compute one immutable GENERAL block
const nlohmann::json &MModelTune::General(const std::string &name) const {
  try {
    return general_.at(name);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MModelTune: missing GENERAL block " + name + ": " + error.what());
  }
}

// Compute one immutable NUMERICS block
const nlohmann::json &MModelTune::Numerics(const std::string &name) const {
  try {
    return numerics_.at(name);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MModelTune: missing NUMERICS block " + name + ": " + error.what());
  }
}

// Compute one immutable MP, XP, GP or TP continuum steering card
const nlohmann::json &MModelTune::Continuum(const std::string &model) const {
  try {
    return continuum_.at(model);
  } catch (const std::exception &error) {
    throw std::invalid_argument("MModelTune: missing continuum model " + model + ": " + error.what());
  }
}

}  // namespace gra
