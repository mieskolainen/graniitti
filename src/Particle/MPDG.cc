// PDG class
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <cctype>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <functional>
#include <iterator>
#include <limits>
#include <map>
#include <set>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Particle/MPDG.h"

// Libraries
#include "json.hpp"
#include "rang.hpp"

using gra::aux::indices;

namespace gra {

namespace {

// Compute true when the doubled spin label is a physical half-integer
//
bool IsFermionSpinX2(int spinX2) {
  return spinX2 >= 0 && spinX2 % 2 != 0;
}

// Compute true for a neutral meson with different valence flavors
//
bool IsOpenFlavorNeutralMeson(const gra::MParticle &particle) {
  const int pdg = std::abs(particle.pdg);
  if (pdg == 130 || pdg == 310) { return false; }
  if (particle.chargeX3 != 0 || particle.spinX2 < 0 ||
      particle.spinX2 % 2 != 0 || pdg < 100) {
    return false;
  }
  const int first_flavor = (pdg / 100) % 10;
  const int second_flavor = (pdg / 10) % 10;
  return first_flavor >= 1 && first_flavor <= 6 && second_flavor >= 1 &&
         second_flavor <= 6 && first_flavor != second_flavor;
}

// Convert a particle width to its proper lifetime
// tau = hbar/Gamma for Gamma > 0
double ParticleTau(double width) { return width > 0.0 ? PDG::hbar / width : 0.0; }

// Decode the LS assignment only for meson codes carrying quark-model labels
// [REFERENCE: PDG 2025, Monte Carlo Particle Numbering Scheme, Table 45.1 and section 5(i)]
void SetMesonPCL(gra::MParticle &p) {
  if (p.pdg / 1000000 == 9) {
    // These enumeration codes specify J^PC without a unique quark-model L
    // [REFERENCE: PDG 2025, Mesons Summary Tables]
    switch (p.pdg) {
      case 9000221: case 9010221: case 9030221: case 9050221:
      case 9000111: case 9000211:
      case 9000311: case 9000321: case 9020311: case 9020321:
        p.setPCL(1, 1, 1); return;
      case 9020221: case 9010111: case 9010211:
        p.setPCL(-1, 1, 0); return;
      case 9010225: case 9050225: case 9060225: case 9070225: case 9080225: case 9090225:
      case 9000115: case 9000215: case 9010315: case 9010325:
        p.setPCL(1, 1, 1); return;
      case 9010113: case 9010213:
        p.setPCL(-1, 1, 0); return;
      case 9020113: case 9020213: case 9000313: case 9000323:
        p.setPCL(1, 1, 1); return;
      case 9000443: case 9010443: case 9020443: case 9000553: case 9010553:
        p.setPCL(-1, -1, 0); return;
      default:
        throw std::invalid_argument("MPDG::ReadParticleData: missing J^PC assignment for meson PDG " +
                                    std::to_string(p.pdg));
    }
  }

  const int J = p.spinX2 / 2;
  const int nL = (p.pdg / 10000) % 10;
  int L = 0;
  int S = 0;
  if (J == 0 && (nL == 0 || nL == 1)) {
    L = S = nL;
  } else if (J > 0 && nL >= 0 && nL <= 3) {
    L = J + (nL == 0 ? -1 : nL == 3 ? 1 : 0);
    S = nL == 1 ? 0 : 1;
  } else {
    throw std::invalid_argument("MPDG::ReadParticleData: invalid meson LS code " + std::to_string(p.pdg));
  }
  p.setPCL(L % 2 == 0 ? -1 : 1, (L + S) % 2 == 0 ? 1 : -1, L);
}

// Decode meson isospin and G parity before clearing the charged-state C label
// [REFERENCE: PDG 2025, Monte Carlo Particle Numbering Scheme, section 5(b)]
void SetMesonIG(gra::MParticle &p) {
  const int first = (p.pdg / 100) % 10;
  const int second = (p.pdg / 10) % 10;
  if ((first == 1 && second == 1) || (first == 2 && second == 1)) {
    p.isospinX2 = 2;
    p.G = -p.C;
  } else if (first == second && first >= 2 && first <= 6) {
    p.isospinX2 = 0;
    p.G = p.C;
  } else if (first >= 3 && first <= 6 && second >= 1 && second <= 2) {
    p.isospinX2 = 1;
  }
}

// Assign baryon parity using the excitation band and explicit table exceptions
// [REFERENCE: PDG 2025, Monte Carlo Particle Numbering Scheme, section 6 and Table 45.2]
void SetBaryonPCL(gra::MParticle &p) {
  // These mass-table IDs do not encode the tabulated quark-model spin assignments
  // [REFERENCE: PDG 2025, Charmed Baryons Summary Tables, Lambda_c(2625), Xi_c(2790) and Xi_c(2815)]
  switch (p.pdg) {
    case 104122: case 104312: case 104322: p.spinX2 = 3; break;
    case 104314: case 104324: p.spinX2 = 1; break;
    default: break;
  }
  const int nr = (p.pdg / 100000) % 10;
  p.setPCL(nr % 2 == 0 ? 1 : -1, 0, nr % 2);

  // The mass table retains excited baryon IDs predating the band encoding
  // [REFERENCE: PDG 2025, Baryons Summary Tables, pdg.lbl.gov/2025/tables/rpp2025-qtab-baryons.pdf]
  switch (p.pdg) {
    case 1214: case 2124: case 22112: case 22212: case 32112: case 32212:
    case 2116: case 2216: case 21214: case 22124: case 1218: case 2128:
    case 1112: case 1212: case 2122: case 2222:
    case 11114: case 12114: case 12214: case 12224:
    case 11112: case 11212: case 12122: case 12222:
    case 11116: case 11216: case 12126: case 12226:
    case 13122: case 3124: case 33122: case 13124: case 43122: case 13126: case 3128:
    case 13114: case 13214: case 13224: case 23112: case 23212: case 23222:
    case 23114: case 23214: case 23224:
    case 3116: case 3216: case 3226: case 13314: case 13324: case 14122:
      p.setPCL(-1, 0, 1); break;
    default: break;
  }
  // Use the minimum constituent orbital compatible with J and intrinsic parity
  int L = std::max(0, (p.spinX2 - 3) / 2);
  if ((L % 2 == 0 ? 1 : -1) != p.P) { ++L; }
  p.L = static_cast<unsigned int>(L);
}

// Read one complete numeric field without accepting a partial conversion
double ReadPDGNumber(const std::string &field, const std::string &context, bool optional = false) {
  if (optional && field.find_first_not_of(" \t\r") == std::string::npos) { return 0.0; }
  std::istringstream input(field);
  double value = 0.0;
  if (!(input >> value) || !(input >> std::ws).eof() || !std::isfinite(value)) {
    throw std::invalid_argument("MPDG::ReadParticleData: invalid numeric " + context);
  }
  return value;
}

// Store the parallel charge states of one fixed-column PDG row
struct PDGRow {
  std::vector<int> id;
  std::vector<int> chargeX3;
  double mass = 0.0;
  double width = 0.0;
  std::string name;
};

// Parse the documented PDG fixed columns, including absent width fields
PDGRow ReadPDGRow(const std::string &line, const std::string &context) {
  if (line.size() < 108) { throw std::invalid_argument("MPDG::ReadParticleData: short row in " + context); }
  PDGRow row;
  for (std::size_t offset = 0; offset < 32; offset += 8) {
    const std::string field = line.substr(offset, 8);
    if (field.find_first_not_of(' ') == std::string::npos) { continue; }
    std::istringstream input(field);
    int pdg = 0;
    if (!(input >> pdg) || !(input >> std::ws).eof() || pdg <= 0) {
      throw std::invalid_argument("MPDG::ReadParticleData: invalid PDG field in " + context);
    }
    row.id.push_back(pdg);
  }
  row.mass = ReadPDGNumber(line.substr(33, 18), "mass in " + context);
  row.width = ReadPDGNumber(line.substr(70, 18), "width in " + context, true);
  if (row.mass < 0.0 || row.width < 0.0) {
    throw std::invalid_argument("MPDG::ReadParticleData: negative mass or width in " + context);
  }
  for (const std::size_t offset : {52U, 61U, 89U, 98U}) {
    (void)ReadPDGNumber(line.substr(offset, 8), "uncertainty in " + context, true);
  }
  std::istringstream names(line.substr(107));
  names >> row.name;
  if (line.at(line.find_last_not_of(" \t\r")) == ',') {
    throw std::invalid_argument("MPDG::ReadParticleData: missing charge after comma in " + context);
  }
  const std::map<std::string, int> charges = {{"--", -6}, {"-", -3}, {"-1/3", -1},
                                               {"0", 0}, {"+2/3", 2}, {"+", 3}, {"++", 6}};
  std::string field;
  while (std::getline(names, field, ',')) {
    field.erase(std::remove_if(field.begin(), field.end(), [](unsigned char c) { return std::isspace(c); }), field.end());
    const auto found = charges.find(field);
    if (found == charges.end()) {
      throw std::invalid_argument("MPDG::ReadParticleData: invalid charge in " + context);
    }
    row.chargeX3.push_back(found->second);
  }
  if (row.id.empty() || row.name.empty() || row.id.size() != row.chargeX3.size()) {
    throw std::invalid_argument("MPDG::ReadParticleData: malformed particle row in " + context);
  }
  return row;
}

// Insert one unique PDG particle and process name
void InsertPDGParticle(std::map<int, gra::MParticle> &table, const gra::MParticle &particle,
                       const std::string &source) {
  if (particle.pdg == 0) {
    throw std::invalid_argument("MPDG::ReadParticleData: PDG id 0 in " + source);
  }
  if (table.count(particle.pdg)) {
    throw std::invalid_argument("MPDG::ReadParticleData: duplicate PDG id " +
                                std::to_string(particle.pdg) + " in " + source);
  }
  for (const auto &[pdg, existing] : table) {
    if (existing.name == particle.name) {
      throw std::invalid_argument("MPDG::ReadParticleData: duplicate particle name " + particle.name +
                                  " for PDGs " + std::to_string(pdg) + " and " + std::to_string(particle.pdg) +
                                  " in " + source);
    }
  }
  table.insert(std::make_pair(particle.pdg, particle));
}

// Insert a particle and its distinct charge conjugate when required
void InsertPDGParticleWithAnti(std::map<int, gra::MParticle> &table, gra::MParticle particle,
                               const std::string &anti_name_seed, bool auto_antiparticle,
                               const std::string &source) {
  particle.tau = ParticleTau(particle.width);
  InsertPDGParticle(table, particle, source);

  if (!auto_antiparticle ||
      !(IsFermionSpinX2(particle.spinX2) || std::abs(particle.chargeX3) > 0 ||
        IsOpenFlavorNeutralMeson(particle))) {
    return;
  }

  gra::MParticle antip = particle;

  if ((1 <= particle.pdg && particle.pdg <= 6) || particle.pdg == PDG::PDG_hard_jet ||
      (particle.pdg == 12 || particle.pdg == 14 || particle.pdg == 16)) {
    antip.name = anti_name_seed + "~";
  } else if (particle.chargeX3 != 0) {
    const std::string signstr = particle.chargeX3 > 0
                                    ? std::string(std::abs(particle.chargeX3 / 3), '-')
                                    : std::string(std::abs(particle.chargeX3 / 3), '+');
    antip.name = anti_name_seed + signstr;
    // Charged baryon multiplets need an explicit antibaryon marker
    if (IsFermionSpinX2(particle.spinX2) && std::abs(particle.pdg) >= 1000 &&
        std::abs(particle.pdg) != PDG::PDG_p) {
      antip.name += "~";
    }
  } else {
    antip.name = anti_name_seed + "0~";
  }

  antip.pdg      = -particle.pdg;
  antip.chargeX3 = -particle.chargeX3;
  if (std::abs(particle.color) == 3 || std::abs(particle.color) == 6) { antip.color = -particle.color; }

  if (IsFermionSpinX2(particle.spinX2)) { antip.P = -particle.P; }

  InsertPDGParticle(table, antip, source);
}

// Read a discrete particle label without narrowing out-of-range JSON integers
int RequiredJsonInt(const nlohmann::json &entry, const std::string &key,
                    const std::string &context) {
  if (!entry.contains(key)) {
    throw std::invalid_argument("MPDG::ReadParticleData: missing '" + key + "' in " + context);
  }
  if (!entry.at(key).is_number_integer()) {
    throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' must be an integer in " +
                                context);
  }
  const auto &value = entry.at(key);
  if (value.is_number_unsigned()) {
    if (value.get<nlohmann::json::number_unsigned_t>() > static_cast<unsigned int>(std::numeric_limits<int>::max())) {
      throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' exceeds int range in " + context);
    }
  } else {
    const auto number = value.get<nlohmann::json::number_integer_t>();
    if (number < -std::numeric_limits<int>::max() || number > std::numeric_limits<int>::max()) {
      throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' exceeds signed particle range in " + context);
    }
  }
  return value.get<int>();
}

// Read one finite particle mass or width
double RequiredJsonDouble(const nlohmann::json &entry, const std::string &key,
                          const std::string &context) {
  if (!entry.contains(key)) {
    throw std::invalid_argument("MPDG::ReadParticleData: missing '" + key + "' in " + context);
  }
  if (!entry.at(key).is_number()) {
    throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' must be numeric in " +
                                context);
  }
  const double value = entry.at(key).get<double>();
  if (!std::isfinite(value)) {
    throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' must be finite in " +
                                context);
  }
  return value;
}

// Read one particle name
std::string RequiredJsonString(const nlohmann::json &entry, const std::string &key,
                               const std::string &context) {
  if (!entry.contains(key)) {
    throw std::invalid_argument("MPDG::ReadParticleData: missing '" + key + "' in " + context);
  }
  if (!entry.at(key).is_string()) {
    throw std::invalid_argument("MPDG::ReadParticleData: '" + key + "' must be a string in " +
                                context);
  }
  return entry.at(key).get<std::string>();
}

// Reject unsupported keys in one extra-particle entry
//
void ValidateExtraParticleKeys(const nlohmann::json &entry, const std::string &context) {
  const std::set<std::string> allowed = {"PDG",      "name",       "mass", "width",
                                                "chargeX3", "spinX2",     "P",    "C",
                                                "isospinX2", "G",         "L",    "color"};
  for (auto it = entry.begin(); it != entry.end(); ++it) {
    if (!allowed.count(it.key())) {
      throw std::invalid_argument("MPDG::ReadParticleData: unknown key '" + it.key() + "' in " +
                                  context);
    }
  }
}

// Validate discrete quantum numbers of one extra particle
//
void ValidateExtraParticleQuantumNumbers(const gra::MParticle &particle, int L,
                                         const std::string &context) {
  if (L < 0) {
    throw std::invalid_argument(
        "MPDG::ReadParticleData: L must be nonnegative in " + context);
  }
  if (!(particle.P == -1 || particle.P == 0 || particle.P == 1) ||
      !(particle.C == -1 || particle.C == 0 || particle.C == 1)) {
    throw std::invalid_argument(
        "MPDG::ReadParticleData: P and C must be -1, 0 or 1 in " + context);
  }
  if (particle.chargeX3 != 0 && particle.C != 0) {
    throw std::invalid_argument(
        "MPDG::ReadParticleData: charged particles cannot have defined C-parity in " + context);
  }
  if (particle.isospinX2 < gra::aux::kNullIsospinX2) {
    throw std::invalid_argument(
        "MPDG::ReadParticleData: isospinX2 must be -1 or nonnegative in " +
        context);
  }
  if (!(particle.G == -1 || particle.G == 0 || particle.G == 1)) {
    throw std::invalid_argument(
        "MPDG::ReadParticleData: G must be -1, 0 or 1 in " + context);
  }
  if (particle.G != 0 && (particle.isospinX2 == gra::aux::kNullIsospinX2 ||
                          particle.isospinX2 % 2 != 0)) {
    throw std::invalid_argument("MPDG::ReadParticleData: defined G-parity "
                                "requires nonnegative integer isospin in " +
                                context);
  }
  if (particle.G != 0 && particle.chargeX3 == 0 && particle.C != 0) {
    const int expected_G =
        ((particle.isospinX2 / 2) % 2 == 0) ? particle.C : -particle.C;
    if (particle.G != expected_G) {
      throw std::invalid_argument("MPDG::ReadParticleData: neutral particle "
                                  "G-parity must equal C*(-1)^I in " +
                                  context);
    }
  }
}

// Read the tune-specific extra-particle table
//
void ReadExtraParticleData(const std::string &modelparam, std::map<int, gra::MParticle> &table) {
  const std::filesystem::path extra_path =
      (modelparam.empty() || modelparam == "null")
          ? std::filesystem::path()
          : std::filesystem::path(gra::ResolveModelDataFile(modelparam, "PDG_EXTRA.json"));
  if (!std::filesystem::exists(extra_path)) { return; }

  nlohmann::json j;
  try {
    j = nlohmann::json::parse(gra::aux::GetInputData(extra_path.string()));
  } catch (const nlohmann::json::exception &e) {
    throw std::invalid_argument("MPDG::ReadParticleData: Error parsing " + extra_path.string() +
                                ": " + e.what());
  }

  if (!j.contains("PARAM_PDG") || !j.at("PARAM_PDG").is_object()) {
    throw std::invalid_argument("MPDG::ReadParticleData: " + extra_path.string() +
                                " must contain a 'PARAM_PDG' object");
  }

  for (auto it = j.at("PARAM_PDG").begin(); it != j.at("PARAM_PDG").end(); ++it) {
    const auto &entry       = it.value();
    const std::string label = extra_path.string() + ":PARAM_PDG[\"" + it.key() + "\"]";
    if (!entry.is_object()) {
      throw std::invalid_argument("MPDG::ReadParticleData: " + label + " must be an object");
    }
    ValidateExtraParticleKeys(entry, label);

    gra::MParticle p;
    p.pdg      = RequiredJsonInt(entry, "PDG", label);
    p.name     = RequiredJsonString(entry, "name", label);
    p.mass     = RequiredJsonDouble(entry, "mass", label);
    p.width    = RequiredJsonDouble(entry, "width", label);
    p.chargeX3 = RequiredJsonInt(entry, "chargeX3", label);
    p.spinX2   = RequiredJsonInt(entry, "spinX2", label);
    p.P        = RequiredJsonInt(entry, "P", label);
    p.C        = RequiredJsonInt(entry, "C", label);
    p.isospinX2 = RequiredJsonInt(entry, "isospinX2", label);
    p.G         = RequiredJsonInt(entry, "G", label);
    const int L = RequiredJsonInt(entry, "L", label);
    p.L         = static_cast<unsigned int>(L);
    p.color     = RequiredJsonInt(entry, "color", label);

    if (p.mass < 0.0 || p.width < 0.0) {
      throw std::invalid_argument("MPDG::ReadParticleData: mass and width must be nonnegative in " +
                                  label);
    }
    if (p.spinX2 < 0 && p.spinX2 != gra::aux::kNullSpinX2) {
      throw std::invalid_argument(
          "MPDG::ReadParticleData: spinX2 must be nonnegative or the null "
          "trajectory marker in " +
          label);
    }
    ValidateExtraParticleQuantumNumbers(p, L, label);

    InsertPDGParticleWithAnti(table, p, p.name, true, label);
  }
}

}  // namespace

// Compute the default PDG mass and width table path
std::string MPDG::DataFile() {
  return gra::aux::ResolveProjectPath("modeldata/mass_width_2026.mcd");
}

// Read the default PDG mass and width table
void MPDG::ReadParticleData() { ReadParticleData(DataFile()); }

// Read PDG particle data in
// Input as the full path to PDG .mcd file
void MPDG::ReadParticleData(const std::string &filepath, const std::string &modelparam) {
  std::cout << "MPDG::ReadParticleData: Reading in PDG tables: ";
  std::ifstream infile(filepath);

  if (!infile.is_open()) {
    throw std::invalid_argument("MPDG::ReadParticleData: Error: Cannot open inputfile " + filepath);
  }

  std::string line;

  // Keep the previous table intact until all input has been validated
  std::map<int, gra::MParticle> table;

  while (std::getline(infile, line)) {
    std::istringstream iss(line);

    // Check * (comment) lines
    std::string first;
    iss >> first;
    if (first.empty() || first.front() == '*') {
      continue;
    }

    const PDGRow row = ReadPDGRow(line, filepath);
    const auto &id = row.id;
    const auto &name = row.name;

    // Loop over different charge assignments
    for (const auto &k : indices(id)) {

      // New particle
      gra::MParticle p;

      // Add properties
      p.name = name;

      p.pdg = id[k];
      p.mass = row.mass;
      p.width = row.width;
      p.chargeX3 = row.chargeX3[k];

      int lastdigit = p.pdg % 10;       // Get last digit
      p.spinX2      = (lastdigit - 1);  // Get 2J
      p.wcut        = 0;                // Off-shell mass (width) cut

      // Neutral mesons/baryons get 0 for their name
      if (std::abs(p.chargeX3) == 0 && p.pdg > 100) {
        p.name = p.name + "0";
      } else if ((p.chargeX3 != 0) &&
                 !((p.pdg >= 1 && p.pdg <= 6) ||
                   (p.pdg == 12 || p.pdg == 14 || p.pdg == 16))) {  // quarks & neutrinos

        std::string signstr = p.chargeX3 > 0 ? std::string(std::abs(p.chargeX3 / 3), '+')
                                             : std::string(std::abs(p.chargeX3 / 3), '-');

        p.name = name + signstr;
      }

      // -------------------------------------------------
      // SPECIAL CASES (not following spin numbering)

      // SM bosons
      if (p.pdg == 21) { // gluon
        p.spinX2 = 2;
        p.color  = 8;
        p.P      = -1;
        p.C      = 0;    // Not defined
      }
      if (p.pdg == 22) { // gamma
        p.spinX2 = 2;
        p.P      = -1;
        p.C      = -1;
      }
      if (std::abs(p.pdg) == 24) {  // W+-
        p.spinX2 = 2;
        p.P      = 0;   // Not defined
        p.C      = 0;   // Not defined
      }
      if (p.pdg == 23) {  // Z
        p.spinX2 = 2;
        p.P      = 0;   // Not defined
        p.C      = 0;   // Not defined
      }
      if (p.pdg == 25) {  // H
        p.spinX2 = 0;
        p.P      = 1;   // Scalar SM higgs
        p.C      = 1;   // Scalar SM higgs
      }

      // SM fermions
      if (p.pdg >= 11 && p.pdg <= 16) {  // leptons and neutrinos
        p.spinX2 = 1;
        p.P      = 1;
        p.C      = 0;   // Not defined
      }
      if (p.pdg >= 1 && p.pdg <= 6) {  // quarks
        p.spinX2 = 1;
        p.color  = 3;
        p.P      = 1;
        p.C      = 0;   // Not defined
      }

      // K_S and K_L are spin-zero pseudoscalars with special PDG codes
      if (p.pdg == 130 || p.pdg == 310) {
        p.spinX2 = 0;
        p.setPCL(-1, 0, 0);
      } else if (p.pdg >= 100 && (p.pdg / 1000) % 10 == 0 && p.spinX2 % 2 == 0) {
        SetMesonPCL(p);
        SetMesonIG(p);
      } else if (p.pdg >= 1000 && IsFermionSpinX2(p.spinX2)) {
        SetBaryonPCL(p);
      }

      // Charged and open-flavor neutral particles are not C-parity eigenstates
      if (p.chargeX3 != 0 || IsOpenFlavorNeutralMeson(p)) { p.C = 0; }

      InsertPDGParticleWithAnti(table, p, name, IsFermionSpinX2(p.spinX2) || id.size() < 3, filepath);
    }
  }  // PDG file line while loop

  if (infile.bad()) {
    throw std::invalid_argument("MPDG::ReadParticleData: Error reading inputfile " + filepath);
  }
  ReadExtraParticleData(modelparam, table);
  PDG_table.swap(table);

  // Process now here fast variable setup >>
  std::cout << rang::fg::green << "[DONE]" << rang::fg::reset << std::endl;
}

// Recursive function to read the decay process string
void MPDG::TokenizeProcess(const std::string &str, int depth,
                           std::vector<gra::MDecayBranch> &branches) const {
  std::vector<std::string> tokens;
  for (std::size_t i = 0; i < str.size();) {
    if (std::isspace(static_cast<unsigned char>(str[i]))) {
      ++i;
      continue;
    }
    if (str[i] == '{' || str[i] == '}' || str[i] == '>') {
      tokens.emplace_back(1, str[i++]);
      continue;
    }
    const std::size_t begin = i;
    while (i < str.size() &&
           !std::isspace(static_cast<unsigned char>(str[i])) && str[i] != '{' &&
           str[i] != '}' && str[i] != '>') {
      ++i;
    }
    tokens.push_back(str.substr(begin, i - begin));
  }

  std::size_t position = 0;
  std::vector<gra::MDecayBranch> parsed;
  std::function<void(int, std::vector<gra::MDecayBranch> &, bool)> parse;
  parse = [&](int current_depth, std::vector<gra::MDecayBranch> &current,
              bool requires_closing_brace) {
    bool has_particle = false;
    while (position < tokens.size()) {
      const std::string &token = tokens[position];
      if (token == "}") {
        if (!requires_closing_brace) {
          throw std::invalid_argument(
              "MPDG::TokenizeProcess: closing brace without a decay block");
        }
        if (!has_particle) {
          throw std::invalid_argument(
              "MPDG::TokenizeProcess: decay daughter block is empty");
        }
        ++position;
        return;
      }
      if (token == "{" || token == ">") {
        throw std::invalid_argument(
            "MPDG::TokenizeProcess: expected a particle before '" + token +
            "'");
      }

      MDecayBranch branch;
      branch.p = FindByPDGName(token);
      branch.depth = current_depth;
      branch.name =
          std::to_string(branch.p.pdg) + "#" + std::to_string(current_depth);
      current.push_back(std::move(branch));
      has_particle = true;
      ++position;

      if (position < tokens.size() && tokens[position] == ">") {
        ++position;
        if (position >= tokens.size() || tokens[position] != "{") {
          throw std::invalid_argument(
              "MPDG::TokenizeProcess: decay arrow must be followed by a "
              "daughter block");
        }
        ++position;
        parse(current_depth + 1, current.back().legs, true);
      }
    }
    if (requires_closing_brace) {
      throw std::invalid_argument(
          "MPDG::TokenizeProcess: decay daughter block is not closed");
    }
  };

  parse(depth, parsed, false);
  branches.insert(branches.end(), std::make_move_iterator(parsed.begin()),
                  std::make_move_iterator(parsed.end()));
}

// Print out PDG table
void MPDG::PrintPDGTable() const {
  std::cout << "MPDG::PrintPDGTable:" << std::endl << std::endl;
  printf("\t\tpdg\tmass\t\twidth\t\tctau\t\tcharge\tJ^PC\tname\n");

  std::map<int, gra::MParticle>::const_iterator it = PDG_table.begin();

  unsigned int counter = 0;
  while (it != PDG_table.end()) {
    const gra::MParticle p = it->second;
    
    printf("%d\t%10d\t%0.6f\t%0.3E\t%0.1E\t\t%2s\t%s%s%s\t%s \n", ++counter, p.pdg, p.mass,
           p.width, PDG::c * p.tau,
           gra::aux::Charge3XtoString(p.chargeX3).c_str(),
           gra::aux::NullableSpin2XtoString(p.spinX2).c_str(), gra::aux::ParityToString(p.P).c_str(),
           gra::aux::ParityToString(p.C).c_str(), p.name.c_str());
    ++it;
  }
}

// Find particle by PDG ID
const gra::MParticle &MPDG::FindByPDG(int pdgcode) const {
  std::map<int, gra::MParticle>::const_iterator it = PDG_table.find(pdgcode);
  if (it != PDG_table.end()) { return it->second; }

  // Throw a fatal error, did not found the particle
  PrintPDGTable();
  std::string str = "MProcess::FindByPDG: Unknown PDG ID: " + std::to_string(pdgcode);
  throw std::invalid_argument(str);
}

// Find particle by PDG name
const gra::MParticle &MPDG::FindByPDGName(const std::string &pdgname) const {
  std::map<int, gra::MParticle>::const_iterator it = PDG_table.begin();
  while (it != PDG_table.end()) {
    if (it->second.name.compare(pdgname) == 0) { return it->second; }
    ++it;
  }

  // Accept the standard MadGraph neutrino names at the process interface
  if (pdgname == "ve") { return FindByPDG(12); }
  if (pdgname == "ve~") { return FindByPDG(-12); }
  if (pdgname == "vm") { return FindByPDG(14); }
  if (pdgname == "vm~") { return FindByPDG(-14); }
  if (pdgname == "vt") { return FindByPDG(16); }
  if (pdgname == "vt~") { return FindByPDG(-16); }
  if (pdgname == "~j") { return FindByPDG(-PDG::PDG_hard_jet); }

  // Check if we have PDG number input (we allow that too)
  if (gra::aux::IsIntegerDigits(pdgname)) {
    const int pdgcode = std::stoi(pdgname);
    return FindByPDG(pdgcode);
  }

  // Throw a fatal error, did not found the particle
  PrintPDGTable();
  std::string str = "MProcess::FindByPDGName: Unknown PDG name: " + pdgname;
  throw std::invalid_argument(str);
}

}  // namespace gra
