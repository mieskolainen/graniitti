// Shared decay-tree topology and normalization operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDECAYTREE_H
#define MDECAYTREE_H

// C++
#include <vector>

// Own
#include "Graniitti/Kinematics/MDecaySampling.h"

namespace gra::decay {

// Branching-ratio product and the number of contributing decay vertices
struct BranchingRatioSummary {
  double product = 1.0;
  unsigned int factors = 0;
};

// Generated cascade phase space factors for one sampled event
struct GeneratedPhaseSpaceSummary {
  double weight = 1.0;
  double two_pi = 1.0;
  double proposal_volume = 1.0;
  int stable_leaves = 0;
};

// Collect stable terminal PDG ids from one decay branch recursively
void CollectStableLeafPDGs(const MDecayBranch &branch, std::vector<int> &pdgs);

// Compute whether the central final state contains explicit cascade daughters
bool HasCascadedDecay(const std::vector<MDecayBranch> &tree);

// Resolve a physical Jacob-Wick decay from the initialized spin and topology settings
DecayStructure JacobWickStructure(const LORENTZSCALAR &lts);

// Compute the generated cascade phase space factors for one event
GeneratedPhaseSpaceSummary
GeneratedPhaseSpace(const std::vector<MDecayBranch> &tree,
                    bool symmetry_proposal_active);

// Compute the QFT statistical factor for identical stable final states
double FinalStateSymmetryFactor(const std::vector<MDecayBranch> &tree,
                                bool isolated);

// Compute the branching-ratio product for one decay branch
BranchingRatioSummary BranchingRatioProduct(const MDecayBranch &branch);

// Compute the branching-ratio product omitted from an isolated cross section
BranchingRatioSummary IsolatedBranchingRatioProduct(const LORENTZSCALAR &lts,
                                                    bool isolated);

// Print one decay tree recursively
void PrintTree(const MDecayBranch &branch);

// Print integrated phase space factors for one branch recursively
void PrintIntegratedPhaseSpace(const MDecayBranch &branch, bool active);

} // namespace gra::decay

#endif
