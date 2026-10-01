// Nuclear masses from evaluated tables and explicit theoretical completion
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MMass.h"

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <stdexcept>
#include <utility>

#include "Graniitti/Tech/MException.h"

namespace gra::nuclear {
namespace {

// Convert atomic mass excess in GeV to a bare nuclear mass, including electron binding
// [REFERENCE: Wang et al., Chinese Phys. C 45 (2021) 030003, Eqs. (1)-(2)]
double BareMass(unsigned int a, unsigned int z, double excess) {
  const double electrons = 1.0e-9 * (14.4381 * std::pow(z, 2.39) + 1.55468e-6 * std::pow(z, 5.35));
  return a * 0.93149410242 + excess - z * 0.00051099895000 + electrons;
}

// Read original fixed-column AME2020 or IAEA FRDM2012 data without changing their values
std::map<std::pair<unsigned int, unsigned int>, double> ReadMassTable(const MassTable& table) {
  std::ifstream input(table.file);
  if (!input) { throw std::invalid_argument("MMass: cannot read mass table " + table.file); }
  std::map<std::pair<unsigned int, unsigned int>, double> mass;
  std::string line;
  for (unsigned int row = 0; std::getline(input, line); ++row) {
    if (line.empty() || line.front() == '#' || (table.type == MassType::AME && row < 36)) { continue; }
    // The FRDM file also contains light AME-only rows without a theoretical mass column
    if (table.type == MassType::FRDM && (line.size() < 43 || line.substr(33, 10).find_first_not_of(' ') == std::string::npos)) { continue; }
    try {
      const auto a = static_cast<unsigned int>(std::stoul(line.substr(table.type == MassType::AME ? 14 : 4, table.type == MassType::AME ? 5 : 4)));
      const auto z = static_cast<unsigned int>(std::stoul(line.substr(table.type == MassType::AME ? 9 : 0, table.type == MassType::AME ? 5 : 4)));
      std::string excess = table.type == MassType::AME ? line.substr(28, 14) : line.substr(33, 10);
      std::replace(excess.begin(), excess.end(), '#', '.');
      const double value = BareMass(a, z, std::stod(excess) * (table.type == MassType::AME ? 1.0e-6 : 1.0e-3));
      if (a == 0 || a > 999 || z > a || !std::isfinite(value) || !(value > 0.0) ||
          !mass.emplace(std::make_pair(a, z), value).second) {
        throw std::invalid_argument("invalid or duplicate isotope");
      }
    } catch (const std::exception& error) {
      throw std::invalid_argument("MMass: " + table.file + ":" + std::to_string(row + 1) + ": " + error.what());
    }
  }
  if (mass.empty()) { throw std::invalid_argument("MMass: empty mass table " + table.file); }
  return mass;
}

// Compute the BWM binding energy in GeV only outside the configured mass tables
// [REFERENCE: arXiv:nucl-th/0405080, Eqs. (2)-(5), without the isotonic shift]
double Binding(unsigned int a, unsigned int z, const MassParam& param) {
  const auto& c = param.bwm;
  const double root = std::cbrt(a), imbalance = static_cast<double>(a) - 2.0 * z;
  const int pairing = a % 2 ? 0 : (z % 2 ? -1 : 1);
  return 1.0e-3 * (c[0] * a - c[1] * root * root - c[2] * z * (z - 1.0) / root -
                  c[3] * imbalance * imbalance / (a * (1.0 + std::exp(-static_cast<double>(a) / param.asymmetry))) +
                  c[4] * pairing * -std::expm1(-static_cast<double>(a) / param.pairing) / std::sqrt(a));
}

}  // namespace

// Merge evaluated and microscopic masses, then complete the finite daughter domain once
MMass::MMass(const MassParam& param) {
  if (param.tables.empty() || !std::isfinite(param.asymmetry) || !(param.asymmetry > 0.0) ||
      !std::isfinite(param.pairing) || !(param.pairing > 0.0) ||
      std::any_of(param.bwm.begin(), param.bwm.end(), [](double value) { return !std::isfinite(value) || !(value > 0.0); })) {
    throw std::invalid_argument("MMass: invalid mass tables or BWM coefficients");
  }
  std::map<std::pair<unsigned int, unsigned int>, double> tabulated;
  for (const auto& table : param.tables) {
    const auto values = ReadMassTable(table);
    tabulated.insert(values.begin(), values.end());
  }
  if (!tabulated.count({1, 0}) || !tabulated.count({1, 1})) {
    throw std::invalid_argument("MMass: neutron and proton masses are required");
  }
  const unsigned int maximum = tabulated.rbegin()->first.first;
  const double neutron = tabulated.at({1, 0}), proton = tabulated.at({1, 1});
  mass_ = MMatrix<double>(maximum + 1, maximum + 1, 0.0);
  for (unsigned int a = 1; a <= maximum; ++a) {
    for (unsigned int z = 0; z <= a; ++z) {
      const auto found = tabulated.find({a, z});
      // Free neutron or proton systems have no nuclear binding, and are not evaporation remnants
      const double value = found != tabulated.end() ? found->second :
          (a - z) * neutron + z * proton - (z > 0 && z < a ? Binding(a, z, param) : 0.0);
      if (!std::isfinite(value) || !(value > 0.0)) { throw std::invalid_argument("MMass: nonphysical completed mass"); }
      mass_[a][z] = value;
    }
  }
}

// A cascade only decreases A, so initialization covers every subsequent mass lookup
// Model masses complete coverage without asserting experimental accuracy outside the tables
double MMass::Mass(unsigned int a, unsigned int z) const {
  if (a == 0 || a >= mass_.size_row() || z > a) { throw PhaseSpaceFailure("MMass: isotope outside initialized mass domain"); }
  return mass_[a][z];
}

}  // namespace gra::nuclear
