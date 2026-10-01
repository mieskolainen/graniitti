// Central Regge production model identifiers
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEPRODUCTIONMODEL_H
#define MREGGEPRODUCTIONMODEL_H

// C++
#include <stdexcept>
#include <string>
#include <string_view>

namespace gra {

// Select the central production dynamics without changing the common vertices
enum class ReggeProductionModel { None, MP, XP, GP, TP };

// Select the physical role of one central production vertex
enum class ReggeVertexRole { Resonance, Continuum };

// Select one typed central production operator basis
enum class ReggeVertexBasis { AutoMinL, AutoMinS, AutoEqualLS, AutoEqualHelicity, Helicity, LS };

// Compute whether the basis generates its central coupling automatically
inline bool IsAutomaticReggeVertexBasis(ReggeVertexBasis basis) { return basis == ReggeVertexBasis::AutoMinL || basis == ReggeVertexBasis::AutoMinS || basis == ReggeVertexBasis::AutoEqualLS || basis == ReggeVertexBasis::AutoEqualHelicity; }

// Compute the steering identifier of one central production model
inline const char *ReggeProductionModelName(ReggeProductionModel model) {
  switch (model) {
    case ReggeProductionModel::MP:
      return "MP";
    case ReggeProductionModel::XP:
      return "XP";
    case ReggeProductionModel::GP:
      return "GP";
    case ReggeProductionModel::TP:
      return "TP";
    case ReggeProductionModel::None:
      return "None";
  }
  throw std::invalid_argument("Unknown Regge production model enum");
}

// Parse one central production model steering identifier
inline ReggeProductionModel ParseReggeProductionModel(const std::string_view identifier) {
  if (identifier == "MP") { return ReggeProductionModel::MP; }
  if (identifier == "XP") { return ReggeProductionModel::XP; }
  if (identifier == "GP") { return ReggeProductionModel::GP; }
  if (identifier == "TP") { return ReggeProductionModel::TP; }
  throw std::invalid_argument("Unknown Regge production model '" + std::string(identifier) + "'");
}

// Compute the steering prefix selected by one production vertex role
inline const char *ReggeVertexRoleName(ReggeVertexRole role) {
  switch (role) {
    case ReggeVertexRole::Resonance:
      return "fusion";
    case ReggeVertexRole::Continuum:
      return "crossed";
  }
  throw std::invalid_argument("Unknown Regge vertex role enum");
}

// Compute the role-independent name of one central production basis
inline const char *ReggeVertexBasisStem(ReggeVertexBasis basis) {
  switch (basis) {
    case ReggeVertexBasis::AutoMinL:
      return "auto_min_L";
    case ReggeVertexBasis::AutoMinS:
      return "auto_min_S";
    case ReggeVertexBasis::AutoEqualLS:
      return "auto_equal_ls";
    case ReggeVertexBasis::AutoEqualHelicity:
      return "auto_equal_helicity";
    case ReggeVertexBasis::Helicity:
      return "helicity";
    case ReggeVertexBasis::LS:
      return "ls";
  }
  throw std::invalid_argument("Unknown Regge vertex basis enum");
}

// Compute the role-qualified steering name of one production basis
inline std::string ReggeVertexBasisName(ReggeVertexBasis basis, ReggeVertexRole role) {
  if (role == ReggeVertexRole::Resonance && basis == ReggeVertexBasis::LS) { return "g_ls"; }
  if (role == ReggeVertexRole::Resonance) { return ReggeVertexBasisStem(basis); }
  return std::string(ReggeVertexRoleName(role)) + "_" + ReggeVertexBasisStem(basis);
}

// Parse one role-qualified central production basis identifier
inline ReggeVertexBasis ParseReggeVertexBasis(const std::string_view basis, ReggeVertexRole role) {
  if (basis == "g_ls" && role == ReggeVertexRole::Resonance) { return ReggeVertexBasis::LS; }
  std::string_view stem = basis;
  if (role == ReggeVertexRole::Continuum) {
    const std::string prefix = std::string(ReggeVertexRoleName(role)) + "_";
    if (!basis.starts_with(prefix)) { throw std::invalid_argument("Production basis '" + std::string(basis) + "' does not match its crossed role"); }
    stem = basis.substr(prefix.size());
  }
  if (stem == "auto_min_L") { return ReggeVertexBasis::AutoMinL; }
  if (stem == "auto_min_S") { return ReggeVertexBasis::AutoMinS; }
  if (stem == "auto_equal_ls") { return ReggeVertexBasis::AutoEqualLS; }
  if (stem == "auto_equal_helicity") { return ReggeVertexBasis::AutoEqualHelicity; }
  if (stem == "helicity") { return ReggeVertexBasis::Helicity; }
  if (stem == "ls" && role == ReggeVertexRole::Continuum) { return ReggeVertexBasis::LS; }
  throw std::invalid_argument("Unknown Regge vertex basis '" + std::string(basis) + "'");
}

}  // namespace gra

#endif
