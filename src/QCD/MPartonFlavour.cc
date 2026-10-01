// Process-local final-state parton flavour selection
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <charconv>
#include <cstdlib>
#include <stdexcept>
#include <string>
#include <system_error>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/QCD/MPartonFlavour.h"
#include "Graniitti/Particle/MPDG.h"

namespace gra {

// Construct the default five-flavour quark and gluon selection
MPartonFlavour::MPartonFlavour()
    : flavours({2, 1, 3, 4, 5, PDG::PDG_gluon}) {}

// Parse one symbolic or numeric physical parton species
int MPartonFlavour::ParseToken(const std::string &input) const {
  std::string token = input;
  aux::TrimExtraSpace(token);
  if (token == "u") { return 2; }
  if (token == "d") { return 1; }
  if (token == "s") { return 3; }
  if (token == "c") { return 4; }
  if (token == "b") { return 5; }
  if (token == "g") { return PDG::PDG_gluon; }
  if (token.empty()) {
    throw std::invalid_argument("@j: empty final-state parton flavour");
  }

  int value = 0;
  const auto result = std::from_chars(token.data(), token.data() + token.size(), value);
  if (result.ec != std::errc() || result.ptr != token.data() + token.size()) {
    throw std::invalid_argument("@j: unsupported final-state parton flavour '" + token + "'");
  }
  if (!((value >= 1 && value <= 5) || value == PDG::PDG_gluon)) {
    throw std::invalid_argument(
        "@j: PDG id " + std::to_string(value) +
        " is not a supported final-state parton species, use 1,2,3,4,5 or 21");
  }
  return value;
}

// Replace the default selection using one validated @j value list
void MPartonFlavour::Configure(const std::vector<std::string> &tokens) {
  if (tokens.empty()) {
    throw std::invalid_argument("@j: final-state parton flavour list is empty");
  }

  std::vector<int> parsed;
  parsed.reserve(tokens.size());
  for (const auto &token : tokens) {
    const int pdg = ParseToken(token);
    if (std::find(parsed.begin(), parsed.end(), pdg) != parsed.end()) {
      throw std::invalid_argument("@j: duplicate final-state parton flavour PDG " +
                                  std::to_string(pdg));
    }
    parsed.push_back(pdg);
  }
  flavours = std::move(parsed);
  configured = true;
}

// Compute the selected quark species without antiparticle signs
std::vector<int> MPartonFlavour::QuarkFlavours() const {
  std::vector<int> quarks;
  for (const int pdg : flavours) {
    if (pdg >= 1 && pdg <= 5) { quarks.push_back(pdg); }
  }
  return quarks;
}

// Compute true when the physical parton species is selected
bool MPartonFlavour::Accepts(int pdg) const {
  int species = 0;
  if (pdg == PDG::PDG_gluon) {
    species = pdg;
  } else if ((pdg >= 1 && pdg <= 5) || (pdg <= -1 && pdg >= -5)) {
    species = std::abs(pdg);
  } else {
    return false;
  }
  return std::find(flavours.begin(), flavours.end(), species) != flavours.end();
}

}  // namespace gra
