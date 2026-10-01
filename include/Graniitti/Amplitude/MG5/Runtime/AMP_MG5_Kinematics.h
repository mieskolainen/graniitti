// Shared kinematics for generated MG5 processes
//
// (c) 2017-2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_KINEMATICS_H
#define AMP_MG5_KINEMATICS_H

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"

namespace gra {
namespace mg5 {

// Collect the stable leaves of a decay branch recursively
inline void CollectStableDecayLeaves(const MDecayBranch &branch,
                                     std::vector<const MDecayBranch *> &leaves) {
  if (branch.legs.empty()) {
    leaves.push_back(&branch);
    return;
  }
  for (const auto &leg : branch.legs) {
    CollectStableDecayLeaves(leg, leaves);
  }
}

// Compute the stable leaves of a decay tree in tree order
inline std::vector<const MDecayBranch *>
StableDecayLeaves(const std::vector<MDecayBranch> &tree) {
  std::vector<const MDecayBranch *> leaves;
  for (const auto &branch : tree) {
    CollectStableDecayLeaves(branch, leaves);
  }
  return leaves;
}

// Extract stable-leaf four-momenta
inline std::vector<M4Vec>
StableLeafMomenta(const std::vector<const MDecayBranch *> &leaves) {
  std::vector<M4Vec> momenta;
  momenta.reserve(leaves.size());
  for (const auto *leaf : leaves) {
    momenta.push_back(leaf->p4);
  }
  return momenta;
}

// Check external momenta against fixed model masses without changing parameters
inline bool OnShellFinal(const std::vector<M4Vec> &momenta, const std::vector<double> &masses) {
  if (masses.size() != momenta.size() + 2) { return false; }
  for (std::size_t i = 0; i < momenta.size(); ++i) {
    const auto &p = momenta[i];
    const double mass2 = masses[i + 2] * masses[i + 2];
    const double spatial2 = p.P3mod2();
    const double residual = p.E() * p.E() - spatial2 - mass2;
    const double scale = std::max({1.0, p.E() * p.E(), spatial2, mass2});
    if (!(p.E() > 0.0) || !std::isfinite(residual) ||
        std::abs(residual) > 64.0 * std::numeric_limits<double>::epsilon() * scale) { return false; }
  }
  return true;
}

// Build MadGraph momentum pointers with [E, px, py, pz] convention
inline void BuildMG5Momenta(const M4Vec &p1, const M4Vec &p2,
                            const std::vector<M4Vec> &pf,
                            std::vector<std::array<double, 4>> &storage,
                            std::vector<double *> &momenta) {
  storage.clear();
  storage.reserve(2 + pf.size());
  storage.push_back({p1.E(), p1.Px(), p1.Py(), p1.Pz()});
  storage.push_back({p2.E(), p2.Px(), p2.Py(), p2.Pz()});
  for (const auto &p : pf) {
    storage.push_back({p.E(), p.Px(), p.Py(), p.Pz()});
  }

  momenta.clear();
  momenta.reserve(storage.size());
  for (auto &p : storage) {
    momenta.push_back(p.data());
  }
}

// Compute true when GRANIITTI samples the complete generated cascade phase space
inline bool UsesFullDecayChainMode(const LORENTZSCALAR &lts) {
  return lts.process.root_decay_mode != RootDecayMode::Isolated;
}

}  // namespace mg5
}  // namespace gra

#endif
