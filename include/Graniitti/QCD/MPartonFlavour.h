// Process-local final-state parton flavour selection
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPARTONFLAVOUR_H
#define MPARTONFLAVOUR_H

// C++
#include <string>
#include <vector>

namespace gra {

// Validated process-local selection for generic final-state parton aliases
class MPartonFlavour {
 public:
  // Construct the default five-flavour quark and gluon selection
  MPartonFlavour();

  // Replace the default selection using one validated @j value list
  void Configure(const std::vector<std::string> &tokens);

  // Compute all selected physical parton PDG ids
  const std::vector<int> &Flavours() const { return flavours; }

  // Compute the selected quark species without antiparticle signs
  std::vector<int> QuarkFlavours() const;

  // Compute true when the physical parton species is selected
  bool Accepts(int pdg) const;

  // Check whether @j replaced the default parton selection
  bool IsConfigured() const { return configured; }

 private:
  // Parse one symbolic or numeric physical parton species
  int ParseToken(const std::string &token) const;

  std::vector<int> flavours;
  bool configured = false;
};

}  // namespace gra

#endif
