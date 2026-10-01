// Process help derived from amplitude definitions and model cards
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPROCESSTABLE_H
#define MPROCESSTABLE_H

#include <map>
#include <string>
#include <vector>

#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Process/MSubProc.h"

namespace gra {

// Read one selected tune for algorithmic process examples
class MProcessTable {
 public:
  // Read particle names, production channels and decays once for the help tables
  explicit MProcessTable(const std::string& modelparam);

  // Compute final-state examples from process definitions and active model cards
  std::vector<std::string> Examples(const MProc& process,
                                   nuclear::CollisionType collision = nuclear::CollisionType::PP) const;

  // Compute compact overview rows in process registration order
  std::vector<std::vector<std::string>> Rows(const MSubProc& subprocess) const;

 private:
  // Construct a physical decay tree from a signed PDG channel key
  std::vector<MDecayBranch> Tree(const std::string& key) const;

  // Expand unstable vector daughters using their largest configured branching ratios
  bool Expand(std::vector<MDecayBranch>& tree) const;

  // Compute steering syntax using the same particle names as the input parser
  std::string Syntax(const std::vector<MDecayBranch>& tree) const;

  MPDG                                  pdg;
  nlohmann::json                        general;
  nlohmann::json                        decays;
  std::map<std::string, nlohmann::json> resonances;
  double                                coupling_min = 0.0;
};

}  // namespace gra

#endif
