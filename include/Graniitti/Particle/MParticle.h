// Particle and branch objects [HEADER ONLY file]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPARTICLE_H
#define MPARTICLE_H

// C++
#include <algorithm>
#include <iostream>
#include <random>
#include <string>
#include <valarray>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Sampling/MCW.h"
#include "Graniitti/Spin/MHELMatrix.h"

namespace gra {

// Event color-flow tags for HepMC/LHE style color lines
struct MColorFlow {
  int flow1 = 0;
  int flow2 = 0;

  // Clear event color-flow tags
  void clear() {
    flow1 = 0;
    flow2 = 0;
  }

  // Compute true if no event color-flow tags are assigned
  bool empty() const { return flow1 == 0 && flow2 == 0; }
};

// Particle class
struct MParticle {
  // These values are read from the PDG-file
  std::string name;
  int         pdg      = 0;  // PDG code
  int         chargeX3 = 0;  // Q x 3
  int         spinX2   = 0;  // J x 2
  int         color    = 0;  // QCD color code
  MColorFlow  color_flow;    // Event color-flow tags

  double mass  = 0.0;
  double width = 0.0;
  double tau   = 0.0;  // hbar / width in seconds, zero marks a stable particle

  // J^PC and I^G quantum numbers
  int          P           = 1;                       // Default even P-parity
  int          C           = 0;                       // Zero when C-parity is undefined
  int          isospinX2   = aux::kNullIsospinX2;    // I x 2, -1 when undefined
  int          G           = 0;                       // Zero when G-parity is undefined
  unsigned int L           = 0;                       // For Mesons/Baryons
  bool         glue        = false;                   // Glueball state

  // Set parity, charge conjugation and orbital angular momentum
  void setPCL(int _P, int _C, unsigned int _L) {
    P = _P;
    C = _C;
    L = _L;
  }

  // Width cut (in Breit-Wigner sampling)
  double wcut = 0.0;

  // Print the particle quantum numbers used by the process setup
  void print() const {
    std::cout << " NAME:   " << name << std::endl;
    std::cout << " ID:     " << pdg << std::endl;
    std::cout << " M:      " << mass << std::endl;
    std::cout << " W:      " << width << std::endl;
    std::cout << " J^PC:   " << aux::NullableSpin2XtoString(spinX2) << "^" << aux::ParityToString(P)
              << aux::ParityToString(C) << std::endl
              << std::endl;
  }
};


// Select the invariant-mass measure for one decay branch
enum class MassProposal { None, BreitWigner, Uniform, Fixed };

// Recursive decay tree branch
struct MDecayBranch {
  // Construct an empty branch with a neutral identifier
  MDecayBranch() = default;

  // Offshell mass picked event by event
  double m_offshell = 0.0;

  // Integral normalization of the event-by-event intermediate mass proposal
  double mass_proposal_norm = 1.0;

  // Exact generated support in invariant mass squared for crossed proposal terms
  double mass_proposal_min2 = 0.0;
  double mass_proposal_max2 = 0.0;

  // None denotes a virtuality already integrated by the central mass map
  MassProposal mass_proposal = MassProposal::None;

  MParticle                 p;               // PDG particle
  M4Vec                     p4;              // 4-momentum
  std::vector<MDecayBranch> legs;            // Daughters
  M4Vec                     decay_position;  // Decay 4-position
  gra::HELMatrix            hel;             // Decay helicity information

  std::string name = "null"; // Stable branch identifier

  // MC weight container
  gra::kinematics::MCW W;

  // Current event decay phase-space weight
  double W_event = 0.0;

  // Decay tree current level
  int depth = 0;
};

// Clear event color-flow tags from one complete decay branch
inline void ClearDecayBranchColorFlow(MDecayBranch &branch) {
  branch.p.color_flow.clear();
  for (auto &leg : branch.legs) {
    ClearDecayBranchColorFlow(leg);
  }
}

// Compute true when one colored particle is an intermediate decay node
inline bool DecayBranchHasColoredIntermediate(const MDecayBranch &branch) {
  if (!branch.legs.empty() && branch.p.color != 0) {
    return true;
  }
  return std::any_of(
      branch.legs.begin(), branch.legs.end(),
      [](const auto &leg) { return DecayBranchHasColoredIntermediate(leg); });
}

// Compute true when an ordered decay tree contains a colored intermediate node
inline bool
DecayTreeHasColoredIntermediate(const std::vector<MDecayBranch> &tree) {
  return std::any_of(tree.begin(), tree.end(), [](const auto &branch) {
    return DecayBranchHasColoredIntermediate(branch);
  });
}

}  // namespace gra

#endif
