// Cascaded spin decays and coherent final-state symmetrization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef SPIN_MHELICITYDECAY_H
#define SPIN_MHELICITYDECAY_H

// C++ standard
#include <complex>
#include <string>
#include <vector>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {
namespace spin {

struct StableLeafAmplitudeTerm {
  std::vector<MDecayBranch> tree;
  double statistics_sign = 1.0;
};

// Construct an event-dependent resonance decay operator
MMatrix<std::complex<double>> ResonanceDecayMatrix(gra::LORENTZSCALAR &lts,
                                                   const gra::PARAM_RES &res,
                                                   const std::string &frame);

// Construct a root helicity decay operator with optional subsequent cascade decays
MMatrix<std::complex<double>>
ResonanceDecayMatrix(gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res,
                     const gra::HELMatrix &root_hel, const std::string &frame, bool cascade = true);

// Store the event-dependent decay operator in a mutable resonance instance
void DecayAmp(gra::LORENTZSCALAR &lts, gra::PARAM_RES &res,
              const std::string &frame);

// Build or return the cached continuum decay matrix for the current event
MMatrix<std::complex<double>> ContinuumDecayMatrix(gra::LORENTZSCALAR &lts,
                                                   const std::string &frame);

// Cache the unique Bose and Fermi stable-leaf assignments for one topology
void PrepareStableLeafSymmetryAssignments(
    gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &reference_tree);

// Apply one cached stable-leaf assignment and reconstruct internal momenta
bool ApplyStableLeafSymmetryAssignment(
    std::vector<MDecayBranch> &tree,
    const std::vector<std::size_t> &assignment);

// Compute the coherent stable-leaf terms for external amplitude builders
std::vector<StableLeafAmplitudeTerm>
StableLeafAmplitudeTrees(gra::LORENTZSCALAR &lts,
                         const std::vector<MDecayBranch> &reference_tree);

// Compute the fixed-width Breit-Wigner product of generated internal branches
std::complex<double> CascadeBWProduct(const MDecayBranch &branch);

// Compute the fixed-width Breit-Wigner product of a generated decay tree
std::complex<double> CascadeBWProduct(const std::vector<MDecayBranch> &tree);

} // namespace spin
} // namespace gra

#endif
