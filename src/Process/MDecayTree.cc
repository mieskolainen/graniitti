// Shared decay-tree topology and normalization operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <map>
#include <string>
#include <vector>

// Own
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MResonance.h"

// Libraries
#include "rang.hpp"

using gra::aux::indices;
using gra::math::msqrt;
using gra::math::pow2;

namespace gra::decay {

// Collect stable terminal PDG ids from one decay branch recursively
void CollectStableLeafPDGs(const MDecayBranch &branch, std::vector<int> &pdgs) {
  if (branch.legs.empty()) {
    pdgs.push_back(branch.p.pdg);
    return;
  }
  for (const auto &leg : branch.legs) {
    CollectStableLeafPDGs(leg, pdgs);
  }
}

// Compute whether the central final state contains explicit cascade daughters
bool HasCascadedDecay(const std::vector<MDecayBranch> &tree) {
  return std::any_of(tree.begin(), tree.end(), [](const MDecayBranch &branch) {
    return !branch.legs.empty();
  });
}

// Resolve a physical Jacob-Wick decay from the initialized spin and topology settings
DecayStructure JacobWickStructure(const LORENTZSCALAR &lts) {
  if (lts.process.root_decay_mode == RootDecayMode::Isolated) { return {DecayType::None}; }
  const bool coherent = lts.decaytree.size() == 2 && lts.process.SPINDEC && lts.amplitude.DECAY_SYM;
  return {coherent ? DecayType::JacobWickCoherent : DecayType::JacobWickIncoherent};
}

namespace {

// Compute true when one internal branch has a generated mass proposal
bool GeneratedMassProposalActive(const MDecayBranch &branch) {
  return !branch.legs.empty() && branch.mass_proposal != MassProposal::None;
}

// Accumulate generated phase space factors from one decay branch
void AccumulateGeneratedPhaseSpace(const MDecayBranch &branch,
                                   const bool symmetry_proposal_active,
                                   GeneratedPhaseSpaceSummary &summary) {
  if (!branch.legs.empty()) {
    if (!std::isfinite(branch.W_event) || branch.W_event <= 0.0) {
      throw PhaseSpaceFailure(
          "decay::GeneratedPhaseSpace: invalid decay phase space weight");
    }
    summary.weight *= branch.W_event;
    // The central mass map already includes the sole production root virtuality
    if (GeneratedMassProposalActive(branch)) { summary.two_pi *= 2.0 * math::PI; }

    if (!symmetry_proposal_active && GeneratedMassProposalActive(branch)) {
      summary.proposal_volume /= MassDensity(branch);
    }

    for (const auto &leg : branch.legs) {
      AccumulateGeneratedPhaseSpace(leg, symmetry_proposal_active,
                                    summary);
    }
  }
  if (branch.legs.empty()) {
    ++summary.stable_leaves;
  }
}

} // namespace

// Compute the generated cascade phase space factors for one event
GeneratedPhaseSpaceSummary
GeneratedPhaseSpace(const std::vector<MDecayBranch> &tree,
                    const bool symmetry_proposal_active) {
  GeneratedPhaseSpaceSummary summary;
  for (const auto &branch : tree) {
    AccumulateGeneratedPhaseSpace(branch, symmetry_proposal_active, summary);
  }
  return summary;
}

// Compute the QFT statistical factor for identical stable final states
// S = product_a n_a!
double FinalStateSymmetryFactor(const std::vector<MDecayBranch> &tree,
                                bool isolated) {
  if (isolated) {
    return 1.0;
  }

  std::vector<int> pdgs;
  for (const auto &branch : tree) {
    CollectStableLeafPDGs(branch, pdgs);
  }

  std::map<int, unsigned int> multiplicity;
  for (const int pdg : pdgs) {
    ++multiplicity[pdg];
  }

  double factor = 1.0;
  for (const auto &[pdg, count] : multiplicity) {
    (void)pdg;
    factor *= gra::math::factorial(count);
  }
  return factor;
}

// Compute the branching-ratio product for one decay branch
// BR = product_v BR_v
BranchingRatioSummary BranchingRatioProduct(const MDecayBranch &branch) {
  BranchingRatioSummary summary;
  if (!branch.legs.empty() && branch.hel.BR_set) {
    summary.product *= branch.hel.BR;
    ++summary.factors;
  }

  for (const auto &leg : branch.legs) {
    const BranchingRatioSummary child = BranchingRatioProduct(leg);
    summary.product *= child.product;
    summary.factors += child.factors;
  }
  return summary;
}

// Compute the branching-ratio product omitted from an isolated cross section
BranchingRatioSummary IsolatedBranchingRatioProduct(const LORENTZSCALAR &lts,
                                                    bool isolated) {
  BranchingRatioSummary summary;
  if (!isolated) {
    return summary;
  }

  if (lts.process.ROOT_RES_ACTIVE && lts.process.ROOT_RES.hel_decay.BR_set) {
    summary.product *= lts.process.ROOT_RES.hel_decay.BR;
    ++summary.factors;
  }
  if (lts.process.RESONANCES.size() == 1) {
    const auto &resonance = lts.process.RESONANCES.begin()->second;
    if (resonance.hel_decay.BR_set) {
      summary.product *= resonance.hel_decay.BR;
      ++summary.factors;
    }
  }

  for (const auto &branch : lts.decaytree) {
    const BranchingRatioSummary child = BranchingRatioProduct(branch);
    summary.product *= child.product;
    summary.factors += child.factors;
  }
  return summary;
}

// Print one decay tree recursively
void PrintTree(const MDecayBranch &branch) {
  const std::string spaces(branch.depth * 2, ' ');
  const std::string lines = spaces + "\\-->  ";

  if ((branch.depth + 1) % 3 == 0) {
    std::cout << rang::fg::blue << lines << rang::fg::reset;
  } else if ((branch.depth + 1) % 2 == 0) {
    std::cout << rang::fg::yellow << lines << rang::fg::reset;
  } else {
    std::cout << rang::fg::green << lines << rang::fg::reset;
  }

  std::printf("M = %0.2E, W = %0.2E (GeV) [<tau> = %0.1E s] PDG = [%s%d | %s] "
              "[Q=%s, J=%s]",
              branch.p.mass, branch.p.width, branch.p.tau,
              branch.p.pdg > 0 ? " " : "", branch.p.pdg, branch.p.name.c_str(),
              aux::Charge3XtoString(branch.p.chargeX3).c_str(),
              aux::NullableSpin2XtoString(branch.p.spinX2).c_str());
  if (!branch.legs.empty()) {
    std::printf(" [BR=%0.3E]", branch.hel.BR);
  }
  std::printf(" \n");

  for (const auto &leg : branch.legs) {
    PrintTree(leg);
  }
}

// Print integrated phase space factors for one branch recursively
void PrintIntegratedPhaseSpace(const MDecayBranch &branch, const bool active) {
  const std::string spaces(branch.depth * 2, ' ');
  const std::string lines = spaces + "|--";

  if (branch.W.Integral() > 0.0) {
    std::printf("%s> \t {1->%lu LIPS}:  %0.3E +- %0.3E ", lines.c_str(),
                branch.legs.size(), branch.W.Integral(),
                branch.W.IntegralError());
    if (!branch.legs.empty()) {
      std::printf("[BR=%0.3E] ", branch.hel.BR);
    }

    if (!active || branch.legs.empty()) {
      std::cout << "[" << rang::fg::red << "INACTIVE " << rang::fg::reset
                << "part of cross-section integral]" << std::endl;
    } else {
      std::cout << "[" << rang::fg::green << "ACTIVE   " << rang::fg::reset
                << "part of integral / (2PI)]" << std::endl;
    }

    for (const auto &leg : branch.legs) {
      PrintIntegratedPhaseSpace(leg, active);
    }
  }
}

} // namespace gra::decay

namespace gra {

namespace {

// Check whether a branch width is large enough for numerical mass sampling
bool HasResolvableWidth(const gra::MDecayBranch &branch, double width_min) {
  return std::isfinite(branch.p.width) && branch.p.width > 0.0 && branch.p.width >= width_min;
}

// Store one pole mass as a discrete fixed-mass proposal
bool SetFixedMassProposal(gra::MDecayBranch &branch, double minimum_mass, double &mass) {
  if (!std::isfinite(branch.p.mass) || branch.p.mass < minimum_mass) { return false; }
  mass                          = branch.p.mass;
  // Integrate the unresolved pole in the narrow-width approximation
  // [REFERENCE: PDG, Review of Particle Physics, Kinematics, narrow-width approximation]
  branch.mass_proposal_norm = !branch.legs.empty() && branch.p.width > 0.0 ? math::PI / (mass * branch.p.width) : 1.0;
  branch.mass_proposal_min2     = gra::math::pow2(mass);
  branch.mass_proposal_max2     = gra::math::pow2(mass);
  branch.mass_proposal = gra::MassProposal::Fixed;
  return true;
}

// Clear any event-local mass proposal stored on one decay branch
void ResetMassProposal(gra::MDecayBranch &branch) {
  branch.mass_proposal = MassProposal::None;
  branch.mass_proposal_norm     = 1.0;
  branch.mass_proposal_min2     = 0.0;
  branch.mass_proposal_max2     = 0.0;
}

// Compute the minimum kinematically allowed mass for one decay branch
double MinimumBranchMass(const gra::MDecayBranch &branch, double offshell, double width_min) {
  const double pole_floor =
      !HasResolvableWidth(branch, width_min) ? branch.p.mass : std::max(0.0, branch.p.mass - offshell * branch.p.width);
  if (branch.legs.empty()) { return pole_floor; }

  double threshold = 0.0;
  for (const auto &leg : branch.legs) { threshold += MinimumBranchMass(leg, offshell, width_min); }
  threshold += 1e-5;
  return std::max(pole_floor, threshold);
}

// Compute the minimum daughter-mass threshold for one decay branch
double MinimumDaughterMassSum(const gra::MDecayBranch &branch, double offshell, double width_min) {
  if (branch.legs.empty()) { return 0.0; }

  double threshold = 0.0;
  for (const auto &leg : branch.legs) { threshold += MinimumBranchMass(leg, offshell, width_min); }
  return threshold + 1e-5;
}

// Compute the continuous virtuality window for one decaying production root
bool SingleCentralRootMassWindow(const gra::MDecayBranch &root, double offshell, double width_min, double &lower,
                                 double &upper) {
  lower = 0.0;
  upper = 0.0;
  if (root.legs.empty() || !(offshell > 0.0) || !HasResolvableWidth(root, width_min) || !std::isfinite(root.p.mass) ||
      !(root.p.mass > 0.0)) {
    return false;
  }

  const double daughter_masses = MinimumDaughterMassSum(root, offshell, width_min);
  lower                        = std::max(daughter_masses, std::max(0.0, root.p.mass - offshell * root.p.width));
  upper                        = root.p.mass + offshell * root.p.width;
  return std::isfinite(lower) && std::isfinite(upper) && upper > lower;
}

// Compute the generated proposal phase-space product for one branch
double GeneratedDecayPhaseSpaceProduct(const gra::MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double out = branch.W_event;
  for (const auto &leg : branch.legs) { out *= GeneratedDecayPhaseSpaceProduct(leg); }
  return out;
}

// Compute the generated proposal phase-space product for the full cascade
double GeneratedStableLeafProposalPhaseSpace(const gra::LORENTZSCALAR &lts) {
  double out = 1.0;
  if (lts.central_phase_space_mode != gra::CentralPhaseSpaceMode::Central && lts.PS_active && lts.DW.GetN() > 0.0 &&
      !lts.decaytree.empty()) {
    out *= lts.DW.Integral();
  }
  for (const auto &branch : lts.decaytree) { out *= GeneratedDecayPhaseSpaceProduct(branch); }
  return out;
}

}  // namespace

// Compute the cascaded phase-space proposal and Monte Carlo factor
// factor = V_proposal W_PS/(2pi)^N
double MProcess::CascadePS() {
  if (!state.lts.PS_active) { return 1.0; }

  const decay::GeneratedPhaseSpaceSummary phase_space = decay::GeneratedPhaseSpace(
      state.lts.decaytree, state.lts.decay_symmetry_proposal_active);
  double proposal_volume = phase_space.proposal_volume;
  if (state.lts.decay_symmetry_proposal_active) {
    const double mixture_density = decay::MixtureDensity(state.lts, state.lts.decaytree);
    if (!std::isfinite(mixture_density) || mixture_density <= 0.0) {
      throw PhaseSpaceFailure("MProcess::CascadePS: invalid stable-leaf mixture density");
    }
    proposal_volume /= mixture_density;
  }
  return proposal_volume * phase_space.weight / phase_space.two_pi;
}

// Cache stable-leaf assignment plan for symmetrized or diagnostic decay
// handling
void MProcess::CacheStableLeafSymmetryAssignments() {
  spin::PrepareStableLeafSymmetryAssignments(state.lts, state.lts.decaytree);
}

// Compute missing stable-leaf histories for an incoherent physical cascade
double MProcess::DecaySymmetryCompensationFactor() {
  try {
    ProcPtr.DecayStructureFor(state.lts);
  } catch (const AmplitudeFailure &) { throw; } catch (const std::exception &error) {
    throw AmplitudeFailure(std::string("MProcess::DecaySymmetryCompensationFactor: ") + error.what());
  }
  if (state.flat_amplitude != 0 || GetISOLATE() ||
      state.lts.decay_structure.type != DecayType::JacobWickIncoherent ||
      state.lts.decaytree.empty()) { return 1.0; }

  CacheStableLeafSymmetryAssignments();
  return static_cast<double>(std::max<std::size_t>(state.lts.decay_symmetry_assignments.size(), 1));
}

// Select one stable-leaf proposal channel for symmetrized cascade decays
void MProcess::PrepareDecaySymmetryProposal() {
  state.lts.decay_symmetry_proposal_active      = false;
  state.lts.decay_symmetry_proposal_index       = 0;
  state.lts.decay_symmetry_proposal_phase_space = 0.0;

  ProcPtr.DecayStructureFor(state.lts);
  if (state.flat_amplitude != 0 || GetISOLATE() || !state.lts.decay_structure.Coherent() ||
      state.lts.decaytree.empty()) { return; }

  CacheStableLeafSymmetryAssignments();

  if (state.lts.decay_symmetry_assignments.size() <= 1) { return; }

  const double u = state.random.U(0.0, static_cast<double>(state.lts.decay_symmetry_assignments.size()));
  state.lts.decay_symmetry_proposal_index =
      std::min(static_cast<std::size_t>(u), state.lts.decay_symmetry_assignments.size() - 1);
  state.lts.decay_symmetry_proposal_active = true;
}

// Apply the selected stable-leaf proposal channel to generated kinematics
bool MProcess::ApplyDecaySymmetryProposal() {
  if (!state.lts.decay_symmetry_proposal_active) { return true; }
  if (state.lts.decay_symmetry_proposal_index >= state.lts.decay_symmetry_assignments.size()) { return false; }

  state.lts.decay_symmetry_proposal_phase_space = GeneratedStableLeafProposalPhaseSpace(state.lts);
  if (!std::isfinite(state.lts.decay_symmetry_proposal_phase_space) ||
      state.lts.decay_symmetry_proposal_phase_space <= 0.0) {
    return false;
  }

  const auto &assignment = state.lts.decay_symmetry_assignments[state.lts.decay_symmetry_proposal_index];
  if (!spin::ApplyStableLeafSymmetryAssignment(state.lts.decaytree, assignment)) { return false; }

  const unsigned int offset = 3;
  for (const auto &i : indices(state.lts.decaytree)) { state.lts.pfinal[i + offset] = state.lts.decaytree[i].p4; }

  return true;
}


// Get intermediate off-shell mass
void MProcess::GetOffShellMass(gra::MDecayBranch &branch, double &mass) {
  mass = -1.0;
  ResetMassProposal(branch);

  const double daughter_masses = MinimumDaughterMassSum(branch, state.offshell_widths, state.width_min);

  // Fix the mass for a numerically narrow state or a zero-width sampling window
  if (!HasResolvableWidth(branch, state.width_min) || std::fpclassify(state.offshell_widths) == FP_ZERO) {
    SetFixedMassProposal(branch, daughter_masses, mass);
    return;
  }

  const unsigned int INNERMAXTRIAL = 1e4;

  const double M     = branch.p.mass;
  const double W     = branch.p.width;
  const double lower = std::max(daughter_masses, std::max(0.0, M - state.offshell_widths * W));
  const double upper = M + state.offshell_widths * W;

  if (!(upper > lower)) {
    if (std::is_eq(upper <=> lower)) { SetFixedMassProposal(branch, daughter_masses, mass); }
    return;
  }

  if (state.flat_mass2) {
    const double lower2 = pow2(lower);
    const double upper2 = std::min(state.lts.s, pow2(upper));
    if (!(upper2 > lower2)) { return; }
    branch.mass_proposal_norm    = upper2 - lower2;
    branch.mass_proposal_min2    = lower2;
    branch.mass_proposal_max2    = upper2;
    branch.mass_proposal = gra::MassProposal::Uniform;
    mass                         = msqrt(state.random.U(lower2, upper2));
    return;
  }

  const double lower2       = pow2(lower);
  const double upper2       = pow2(upper);
  branch.mass_proposal_norm = gra::MRandom::RelativisticBWMass2Integral(M, W, lower2, upper2);
  if (!std::isfinite(branch.mass_proposal_norm) || branch.mass_proposal_norm <= 0.0) { return; }
  branch.mass_proposal_min2  = lower2;
  branch.mass_proposal_max2  = upper2;
  branch.mass_proposal = gra::MassProposal::BreitWigner;

  unsigned int innertrials = 0;
  while (innertrials < INNERMAXTRIAL) {
    const double candidate = state.random.RelativisticBWRandom(M, W, state.offshell_widths, daughter_masses);
    if (candidate >= daughter_masses && candidate >= lower && candidate <= upper) {
      mass = candidate;
      return;
    }
    ++innertrials;
  }
}

// Prepare top-level branch masses without independently sampling a sole root
bool MProcess::PrepareCentralBranchMasses(double &mass_sum, double *mass_max) {
  mass_sum = 0.0;
  if (state.lts.decaytree.size() == 1) {
    MDecayBranch &root  = state.lts.decaytree.front();
    double        lower = 0.0;
    double        upper = 0.0;
    if (!SingleCentralRootMassWindow(root, state.offshell_widths, state.width_min, lower, upper)) { return false; }
    ResetMassProposal(root);
    root.m_offshell = lower;
    mass_sum        = lower;
    if (mass_max != nullptr) { *mass_max = std::min(*mass_max, upper); }
    return true;
  }

  for (auto &branch : state.lts.decaytree) {
    GetOffShellMass(branch, branch.m_offshell);
    if (branch.m_offshell < 0.0) { return false; }
    mass_sum += branch.m_offshell;
  }
  return true;
}

// Bind one sole production root to the complete central-system momentum
bool MProcess::SetSingleCentralRootKinematics() {
  if (state.lts.decaytree.size() != 1 || state.lts.pfinal.size() < 4) { return false; }
  const M4Vec &central = state.lts.pfinal[0];
  const double mass2   = central.M2();
  if (!std::isfinite(central.Px()) || !std::isfinite(central.Py()) || !std::isfinite(central.Pz()) ||
      !std::isfinite(central.E()) || !std::isfinite(mass2) || !(mass2 > 0.0)) {
    return false;
  }

  MDecayBranch &root  = state.lts.decaytree.front();
  double        lower = 0.0;
  double        upper = 0.0;
  if (!SingleCentralRootMassWindow(root, state.offshell_widths, state.width_min, lower, upper)) { return false; }
  const double mass = msqrt(mass2);
  const double tolerance =
      128.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(lower), std::abs(upper), mass});
  if (mass < lower - tolerance || mass > upper + tolerance) { return false; }

  root.p4         = central;
  root.m_offshell = mass;
  ResetMassProposal(root);
  state.lts.pfinal[3] = central;
  state.lts.DW        = kinematics::MCW(1.0);
  return true;
}

// Compute the direct central phase-space volume for the current root count
double MProcess::CentralDecayWidthPS() const {
  if (state.lts.decaytree.size() == 1) {
    double lower = 0.0;
    double upper = 0.0;
    if (!SingleCentralRootMassWindow(state.lts.decaytree.front(), state.offshell_widths, state.width_min, lower,
                                     upper)) {
      return 0.0;
    }
    const double mass = state.lts.decaytree.front().p4.M();
    if (!std::isfinite(mass) || mass < lower || mass > upper) { return 0.0; }
    return 1.0;
  }
  if (state.lts.decaytree.size() == 2) {
    return kinematics::PS2Massive(state.lts.m2, pow2(state.lts.decaytree[0].p4.M()),
                                  pow2(state.lts.decaytree[1].p4.M()));
  }
  if (state.lts.decaytree.size() > 2) { return kinematics::PSnMassless(state.lts.m2, state.lts.decaytree.size()); }
  return 0.0;
}


// Recursive decay tree kinematics (called event by event from inhereting
// classes)
bool MProcess::ConstructDecayKinematics(gra::MDecayBranch &branch) {
  branch.W_event = 0.0;

  // This leg has any daughters
  if (branch.legs.size() != 0) {
    // Generate decay product masses
    std::vector<double> m(branch.legs.size(), 0.0);

    for (const auto &i : indices(branch.legs)) {
      GetOffShellMass(branch.legs[i], branch.legs[i].m_offshell);
      if (branch.legs[i].m_offshell < 0.0) { return false; }
      m[i] = branch.legs[i].m_offshell;
    }
    if (branch.p4.M() < gra::Sum(m)) { return false; }

    // Active decay phase space must carry its local event weight into the
    // integrand
    const bool         UNWEIGHT = !state.lts.PS_active;
    std::vector<M4Vec> p;

    // 2-body
    gra::kinematics::MCW w;
    if (branch.legs.size() == 2) {
      w = gra::kinematics::TwoBodyPhaseSpace(branch.p4, branch.p4.M(), m, p, state.random);
      // 3-body
    } else if (branch.legs.size() == 3) {
      w = gra::kinematics::ThreeBodyPhaseSpace(branch.p4, branch.p4.M(), m, p, UNWEIGHT, state.random);
      // N-body
    } else {
      w = gra::kinematics::NBodyPhaseSpace(branch.p4, branch.p4.M(), m, p, UNWEIGHT, state.random);
    }
    if (w.GetW() < 0) {
      std::string str =
          "MProcess::ConstructDecayKinematics: Fatal error: "
          "Weight < 0 (Check your decay tree)";
      std::cout << str << std::endl;
      return false;
    }

    // Collect weight
    branch.W_event = w.Integral();
    branch.W += w;

    // Collect decay product 4-momenta
    for (const auto &i : indices(branch.legs)) { branch.legs[i].p4 = p[i]; }

    //  ** Now the INNER recursion **
    for (const auto &i : indices(branch.legs)) {
      if (!ConstructDecayKinematics(branch.legs[i])) { return false; }
    }
  }
  return true;
}


}  // namespace gra
