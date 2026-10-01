// Cascaded spin decay amplitudes and coherent final-state symmetrization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// Identical stable leaves are coherently reassigned across cascade branches
// Crossed terms carry relative fixed-width relativistic Breit-Wigner phases

#include "Graniitti/Spin/MHelicityDecay.h"

// C++ standard
#include <algorithm>
#include <cmath>
#include <complex>
#include <functional>
#include <limits>
#include <map>
#include <numeric>
#include <set>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Math/MCombinatorics.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MSpecialFunctions.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MSpin.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

using gra::aux::indices;

using gra::math::msqrt;

namespace gra {
namespace spin {

namespace {

using FrameTransform = std::function<gra::M4Vec(const gra::M4Vec &)>;

// Rotate one four-vector so the selected axis becomes the helicity z-axis
// The z-axis follows the Jacob-Wick mother direction
//
gra::M4Vec RotateToHelicityAxis(const gra::M4Vec &p, const gra::M4Vec &axis) {
  gra::M4Vec out = p;
  out.RotateZ(-axis.Phi());
  out.RotateY(-axis.Theta());
  out.RotateZ(gra::math::PI);
  return out;
}

// Apply one Jacob-Wick helicity-frame step in the current grandmother frame
// Daughter helicities are measured after the standard JW boost
//
gra::M4Vec ApplyHelicityStep(const gra::M4Vec &p, const gra::M4Vec &axis,
                             const gra::M4Vec &system) {
  gra::M4Vec out = RotateToHelicityAxis(p, axis);
  gra::M4Vec boost = RotateToHelicityAxis(system, axis);
  return gra::kinematics::BoostToRestFrame(out, boost,
                                           "MHelicityDecay::ApplyHelicityStep");
}

// Calculate a branch decay matrix from daughters already transformed to the mother frame
// The branch amplitude is D* times the reduced helicity coupling
//
MMatrix<std::complex<double>> CalculateFMatrixInFrame(
    const MDecayBranch &branch, const gra::M4Vec &branch_axis_in_mother,
    const gra::M4Vec &legA_in_mother, const gra::M4Vec &legB_in_mother) {
  const unsigned int A = 0; // left  '<'

  std::vector<M4Vec> daughters_in_branch = {legA_in_mother, legB_in_mother};
  gra::kinematics::HXframe(daughters_in_branch, branch_axis_in_mother);

  return fDecayMatrix(branch.hel, daughters_in_branch[A].Theta(),
                      daughters_in_branch[A].Phi()) *
         branch.hel.g_decay;
}

// Build the identity carried by one stable final-state helicity space
MMatrix<std::complex<double>>
StableHelicityIdentity(const MDecayBranch &branch, const std::string &context) {
  return MMatrix<std::complex<double>>::IdentityMatrix(
      FinalStateHelicityCount(branch.p, context));
}

// Test whether an ordered two-body vertex terminates in stable daughters
bool HasStableDaughterPair(const MDecayBranch &left,
                           const MDecayBranch &right) {
  return left.legs.empty() && right.legs.empty();
}

// Recursively contract one binary decay branch in its local helicity frames
// Daughter operators generally do not commute, so the left then right
// Kronecker order and the daughter before mother multiplication must not change
MMatrix<std::complex<double>>
CascadeBranchOperator(const MDecayBranch &branch,
                      const FrameTransform &to_mother_frame,
                      const MMatrix<std::complex<double>> &parent_axis) {
  if (branch.legs.empty()) {
    return StableHelicityIdentity(
        branch,
        "MHelicityDecay::CascadeBranchOperator stable branch " + branch.name);
  }

  const gra::M4Vec axis = to_mother_frame(branch.p4);
  const gra::M4Vec left_p4 = to_mother_frame(branch.legs[0].p4);
  const gra::M4Vec right_p4 = to_mother_frame(branch.legs[1].p4);
  const gra::M4Vec system = left_p4 + right_p4;
  // Match the produced helicity ket to the actual HX rest-frame spin axes
  const auto inverse_rotation = SpinHalfRotation(axis.Theta(), axis.Phi()) *
                                SpinHalfRotation(0.0, -gra::math::PI);
  const auto spin = SpinRep::FromSpin(0.5 * branch.p.spinX2, "CascadeBranchOperator");
  std::vector<std::size_t> rows;
  for (const double projection : branch.hel.Jz_values) {
    rows.push_back(spin.Index(projection, "CascadeBranchOperator decay spin"));
  }
  std::vector<std::size_t> columns;
  for (const double helicity : FinalStateHelicities(branch.p, "CascadeBranchOperator")) {
    columns.push_back(spin.Index(helicity, "CascadeBranchOperator produced helicity"));
  }
  const auto parent_rotation = SpinRotation(inverse_rotation.Dagger() * parent_axis, spin.Spin())
                                   .SelectRows(rows).SelectColumns(columns);
  const MMatrix<std::complex<double>> vertex =
      CalculateFMatrixInFrame(branch, axis, left_p4, right_p4) * parent_rotation;
  if (HasStableDaughterPair(branch.legs[0], branch.legs[1])) {
    return vertex;
  }
  const FrameTransform to_branch_frame = [to_mother_frame, axis,
                                          system](const gra::M4Vec &p) {
    return ApplyHelicityStep(to_mother_frame(p), axis, system);
  };
  const auto left_in_branch = to_branch_frame(branch.legs[0].p4);
  const auto daughter_axis = SpinHalfRotation(left_in_branch.Theta(), left_in_branch.Phi());
  const MMatrix<std::complex<double>> flip = {{0.0, 1.0}, {1.0, 0.0}};
  const MMatrix<std::complex<double>> left =
      CascadeBranchOperator(branch.legs[0], to_branch_frame, daughter_axis);
  const MMatrix<std::complex<double>> right =
      CascadeBranchOperator(branch.legs[1], to_branch_frame, daughter_axis * flip);
  return left.KroneckerMultiply(right, vertex);
}

// Create a transform from the lab frame into the rest frame of a system
// CM-frame decay angles use only the parent rest boost
//
FrameTransform RestFrameTransform(const gra::M4Vec &system) {
  return [system](const gra::M4Vec &p) {
    return gra::kinematics::BoostToRestFrame(
        p, system, "MHelicityDecay::RestFrameTransform");
  };
}

// Build the root decay transform matching the configured central spin basis
// CS, HX and CM choose the spin quantization axes of X
//
// [REFERENCE: Trueman, Phys. Rev. D 18 (1978) 3423]
// [REFERENCE: WA102 Collaboration, Phys. Lett. B 467 (1999) 165]
// [REFERENCE: WA76 Collaboration, CERN-EP-89-117]
// [REFERENCE: DM2 Collaboration, Nucl. Phys. B 292 (1987) 653]
FrameTransform RootDecayFrameTransform(const gra::LORENTZSCALAR &lts,
                                       const std::string &frame,
                                       const std::string &label) {
  if (frame == "CM") {
    return RestFrameTransform(lts.pfinal[0]);
  }

  if (frame == "HX") {
    return [system = lts.pfinal[0]](const gra::M4Vec &p) {
      return ApplyHelicityStep(p, system, system);
    };
  }

  if (frame == "CS") {
    return [system = lts.pfinal[0], beam1 = lts.pbeam1,
            beam2 = lts.pbeam2](const gra::M4Vec &p) {
      std::vector<gra::M4Vec> out = {p};
      gra::kinematics::CSframe(out, system, beam1, beam2);
      return out[0];
    };
  }

  throw std::invalid_argument(label + ": Unsupported spin decay frame '" +
                              frame + "'; allowed values are CS, HX, CM");
}

// Build an incoherent spin-independent decay operator with a purified parent-spin label
// The decay operator obeys D D^dagger = I/(2J+1) and hidden columns prevent coherent Jz interference
MMatrix<std::complex<double>> SpinBlindDecayMatrix(std::size_t rows,
                                                   std::size_t cols) {
  if (rows == 0 || cols == 0) {
    return MMatrix<std::complex<double>>();
  }

  MMatrix<std::complex<double>> out(rows, rows * cols, 0.0);
  const std::complex<double> weight =
      1.0 / msqrt(static_cast<double>(rows * cols));
  for (std::size_t spin = 0; spin < rows; ++spin) {
    for (std::size_t final = 0; final < cols; ++final) {
      out[spin][spin * cols + final] = weight;
    }
  }
  return out;
}

// Compute the product of terminal helicity dimensions below one branch
// Each stable leaf contributes its 2s+1 helicity states
//
std::size_t LeafHelicityDimension(const gra::MDecayBranch &branch) {
  if (branch.legs.empty()) {
    return FinalStateHelicityCount(branch.p,
                                   "LeafHelicityDimension terminal branch");
  }

  std::size_t dim = 1;
  for (const auto &child : branch.legs) {
    dim *= LeafHelicityDimension(child);
  }
  return dim;
}

// Preserve stable helicities and average only an unstable parent spin
MMatrix<std::complex<double>> SpinBlindBranch(const gra::MDecayBranch &branch) {
  if (branch.legs.empty()) { return StableHelicityIdentity(branch, "SpinBlindBranch"); }
  return SpinBlindDecayMatrix(FinalStateHelicityCount(branch.p, "SpinBlindBranch"), LeafHelicityDimension(branch));
}

// Compute the full terminal helicity dimension of a central decay tree
// Final-state spin space is the tensor product of leaf spaces
//
std::size_t FinalStateHelicityDimension(const gra::LORENTZSCALAR &lts) {
  std::size_t dim = 1;
  for (const auto &branch : lts.decaytree) {
    dim *= LeafHelicityDimension(branch);
  }
  return dim;
}

using DecayPath = std::vector<std::size_t>;

struct StableLeafSlot {
  DecayPath path;
  int pdg = 0;
  int spinX2 = 0;
  M4Vec p4;
};

struct StableLeafTerm {
  std::vector<MDecayBranch> tree;
  double statistics_sign = 1.0;
  std::vector<std::size_t> assignment;
};

struct StableLeafAssignment {
  std::vector<std::size_t> assignment;
  double statistics_sign = 1.0;
};

struct StableLeafTopologyPlan {
  std::vector<StableLeafAssignment> assignments;
};

// Collect stable leaf slots and their paths in a decay tree
// Identical-particle symmetrization acts only on stable leaves
//
void CollectStableLeafSlots(const MDecayBranch &branch, DecayPath path,
                            std::vector<StableLeafSlot> &slots) {
  if (branch.legs.empty()) {
    slots.push_back({path, branch.p.pdg, branch.p.spinX2, branch.p4});
    return;
  }
  for (const auto &i : indices(branch.legs)) {
    DecayPath child_path = path;
    child_path.push_back(i);
    CollectStableLeafSlots(branch.legs[i], child_path, slots);
  }
}

// Collect stable leaf slots from all top-level central branches
// Top-level branches together define the central final state
//
void CollectStableLeafSlots(const std::vector<MDecayBranch> &tree,
                            std::vector<StableLeafSlot> &slots) {
  for (const auto &i : indices(tree)) {
    CollectStableLeafSlots(tree[i], DecayPath{i}, slots);
  }
}

// Compute a mutable branch reference by its stable decay-tree path
// Leaf permutations replace momenta at fixed decay topology
//
MDecayBranch &BranchAtPath(std::vector<MDecayBranch> &tree,
                           const DecayPath &path) {
  MDecayBranch *branch = &tree.at(path.at(0));
  for (std::size_t i = 1; i < path.size(); ++i) {
    branch = &branch->legs.at(path[i]);
  }
  return *branch;
}

// Recompute one branch four-momentum from its current daughters
// Permuted leaf momenta determine internal off-shell masses
//
M4Vec RecomputeBranchMomentum(MDecayBranch &branch) {
  if (branch.legs.empty()) {
    if (!branch.p4.PhysicalMass(branch.m_offshell)) {
      throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid daughter momentum");
    }
    return branch.p4;
  }

  M4Vec sum(0.0, 0.0, 0.0, 0.0);
  for (auto &leg : branch.legs) {
    sum += RecomputeBranchMomentum(leg);
  }
  branch.p4 = sum;
  if (!branch.p4.PhysicalMass(branch.m_offshell)) {
    throw PhaseSpaceFailure("Stable-leaf symmetrization: invalid branch momentum");
  }
  return branch.p4;
}

// Recompute all internal four-momenta in a permuted decay tree
// Every coherent permutation gets its own internal kinematics
//
void RecomputeTreeMomenta(std::vector<MDecayBranch> &tree) {
  for (auto &branch : tree) {
    RecomputeBranchMomentum(branch);
  }
}

// Build a topology key that identifies isomorphic decay subtrees
// Symmetric boson subtrees can be canonically reordered
//
std::string DecayTopologyKey(const MDecayBranch &branch) {
  if (branch.legs.empty()) {
    return std::to_string(branch.p.pdg);
  }

  std::vector<std::string> child_keys;
  for (const auto &leg : branch.legs) {
    child_keys.push_back(DecayTopologyKey(leg));
  }
  std::sort(child_keys.begin(), child_keys.end());

  std::string out = std::to_string(branch.p.pdg) + "{";
  for (const auto &key : child_keys) {
    out += key + ",";
  }
  out += "}";
  return out;
}

// Check whether a doubled spin quantum number corresponds to a physical fermion
// Odd spinX2 means half-integer spin and Fermi statistics
//
bool IsFermionSpinX2(int spinX2) { return spinX2 >= 0 && (spinX2 % 2) == 1; }

// Count stable fermion leaves by PDG code under one subtree
// Repeated identical fermions determine antisymmetry signs
//
void CollectFermionLeafCounts(const MDecayBranch &branch,
                              std::map<int, std::size_t> &counts) {
  if (branch.legs.empty()) {
    if (IsFermionSpinX2(branch.p.spinX2)) {
      ++counts[branch.p.pdg];
    }
    return;
  }
  for (const auto &leg : branch.legs) {
    CollectFermionLeafCounts(leg, counts);
  }
}

// Test whether a subtree contains repeated identical stable fermion leaves
// Repeated fermions prevent blind canonical bosonic sorting
//
bool HasRepeatedFermionLeaf(const MDecayBranch &branch) {
  std::map<int, std::size_t> counts;
  CollectFermionLeafCounts(branch, counts);
  return std::any_of(counts.begin(), counts.end(),
                     [](const auto &entry) { return entry.second > 1; });
}

// Test whether a full decay tree contains repeated identical stable fermion leaves
// Full-tree fermion repeats control global assignment signs
//
bool HasRepeatedFermionLeaf(const std::vector<MDecayBranch> &tree) {
  std::map<int, std::size_t> counts;
  for (const auto &branch : tree) {
    CollectFermionLeafCounts(branch, counts);
  }
  return std::any_of(counts.begin(), counts.end(),
                     [](const auto &entry) { return entry.second > 1; });
}

// Build a canonical key for one stable-leaf assignment under a branch
// Equivalent Bose assignments are identified coherently
//
std::string
LeafAssignmentKey(const std::map<DecayPath, std::size_t> &slot_by_path,
                  const std::vector<std::size_t> &assignment,
                  const MDecayBranch &branch, const DecayPath &path) {
  if (branch.legs.empty()) {
    const auto it = slot_by_path.find(path);
    if (it == slot_by_path.end()) {
      throw std::invalid_argument(
          "Stable-leaf symmetrization: leaf path not found");
    }
    return std::to_string(branch.p.pdg) + "#" +
           std::to_string(assignment[it->second]);
  }

  std::vector<std::pair<std::string, std::string>> child_keys;
  for (const auto &i : indices(branch.legs)) {
    DecayPath child_path = path;
    child_path.push_back(i);
    child_keys.push_back({DecayTopologyKey(branch.legs[i]),
                          LeafAssignmentKey(slot_by_path, assignment,
                                            branch.legs[i], child_path)});
  }

  // Identical stable siblings already obey the two-body LS statistics rule
  const bool stable_pair = branch.legs.size() == 2 &&
      HasStableDaughterPair(branch.legs[0], branch.legs[1]);
  if (stable_pair || !HasRepeatedFermionLeaf(branch)) {
    std::stable_sort(child_keys.begin(), child_keys.end(),
                     [](const auto &a, const auto &b) {
                       return (a.first == b.first) ? (a.second < b.second)
                                                   : (a.first < b.first);
                     });
  }

  std::string out = std::to_string(branch.p.pdg) + "{";
  for (const auto &key : child_keys) {
    out += key.second + ",";
  }
  out += "}";
  return out;
}

// Build a canonical key for one stable-leaf assignment over the whole tree
// The full final-state assignment is keyed independent of branch
// labels
//
std::string
RootAssignmentKey(const std::vector<MDecayBranch> &tree,
                  const std::map<DecayPath, std::size_t> &slot_by_path,
                  const std::vector<std::size_t> &assignment) {
  std::vector<std::pair<std::string, std::string>> child_keys;
  for (const auto &i : indices(tree)) {
    child_keys.push_back(
        {DecayTopologyKey(tree[i]),
         LeafAssignmentKey(slot_by_path, assignment, tree[i], DecayPath{i})});
  }
  const bool stable_pair = tree.size() == 2 && HasStableDaughterPair(tree[0], tree[1]);
  if (stable_pair || !HasRepeatedFermionLeaf(tree)) {
    std::stable_sort(child_keys.begin(), child_keys.end(),
                     [](const auto &a, const auto &b) {
                       return (a.first == b.first) ? (a.second < b.second)
                                                   : (a.first < b.first);
                     });
  }

  std::string out = "ROOT{";
  for (const auto &key : child_keys) {
    out += key.second + ",";
  }
  out += "}";
  return out;
}

// Build a deterministic string key for one decay path
// Paths identify which leaf momentum is reassigned
std::string DecayPathKey(const DecayPath &path) {
  std::string out;
  for (const auto index : path) {
    out += std::to_string(index) + ".";
  }
  return out;
}

// Build a cache key for stable-leaf assignment topology
// Event-local symmetry caches are valid only for the same topology
std::string StableLeafTopologyCacheKey(const std::vector<MDecayBranch> &tree) {
  std::vector<StableLeafSlot> slots;
  CollectStableLeafSlots(tree, slots);

  std::string key = "ROOT{";
  for (const auto &branch : tree) {
    key += DecayTopologyKey(branch) + ",";
  }
  key += "}|";
  for (const auto &slot : slots) {
    key += DecayPathKey(slot.path) + ":" + std::to_string(slot.pdg) + ":" +
           std::to_string(slot.spinX2) + ";";
  }
  return key;
}

// Generate coherent stable-leaf assignments with Bose/Fermi statistics signs
// Amplitudes are summed over indistinguishable final-state assignments
//
void GenerateStableLeafAssignments(
    const std::vector<std::vector<std::size_t>> &groups,
    std::size_t group_index, const std::vector<bool> &fermion_group,
    std::vector<std::size_t> &assignment,
    const std::vector<MDecayBranch> &reference_tree,
    const std::map<DecayPath, std::size_t> &slot_by_path,
    std::set<std::string> &seen, double statistics_sign,
    std::vector<StableLeafAssignment> &out) {
  if (group_index == groups.size()) {
    const std::string key =
        RootAssignmentKey(reference_tree, slot_by_path, assignment);
    if (seen.insert(key).second) {
      out.push_back({assignment, statistics_sign});
    }
    return;
  }

  const auto &slots = groups[group_index];
  std::vector<std::size_t> perm = slots;
  do {
    for (const auto &i : indices(slots)) {
      assignment[slots[i]] = perm[i];
    }
    const double group_sign =
        fermion_group[group_index] ? static_cast<double>(math::PermutationSign(slots, perm)) : 1.0;
    GenerateStableLeafAssignments(groups, group_index + 1, fermion_group,
                                  assignment, reference_tree, slot_by_path,
                                  seen, statistics_sign * group_sign, out);
  } while (std::next_permutation(perm.begin(), perm.end()));
}

// Build the stable-leaf assignment topology and statistics signs
// Each assignment is one coherent Bose/Fermi amplitude term
//
StableLeafTopologyPlan
BuildStableLeafTopologyPlan(const std::vector<MDecayBranch> &reference_tree) {
  std::vector<StableLeafSlot> slots;
  CollectStableLeafSlots(reference_tree, slots);

  std::map<int, std::vector<std::size_t>> slots_by_pdg;
  std::map<DecayPath, std::size_t> slot_by_path;
  for (const auto &i : indices(slots)) {
    slots_by_pdg[slots[i].pdg].push_back(i);
    slot_by_path[slots[i].path] = i;
  }

  std::vector<std::vector<std::size_t>> groups;
  std::vector<bool> fermion_group;
  for (const auto &[pdg, group] : slots_by_pdg) {
    (void)pdg;
    if (group.size() > 1) {
      const int spinX2 = slots[group.front()].spinX2;
      if (spinX2 < 0 ||
          std::any_of(group.begin(), group.end(), [&](std::size_t slot) {
            return slots[slot].spinX2 != spinX2;
          })) {
        throw std::invalid_argument(
            "Stable-leaf symmetrization: identical PDG leaves have "
            "inconsistent spin metadata");
      }
      groups.push_back(group);
      fermion_group.push_back(IsFermionSpinX2(spinX2));
    }
  }

  std::vector<std::size_t> assignment(slots.size(), 0);
  std::iota(assignment.begin(), assignment.end(), 0);
  if (slots.size() <= 1 || groups.empty()) {
    return {{{assignment, 1.0}}};
  }

  std::set<std::string> seen;
  std::vector<StableLeafAssignment> assignments;
  GenerateStableLeafAssignments(groups, 0, fermion_group, assignment,
                                reference_tree, slot_by_path, seen, 1.0,
                                assignments);

  return {assignments};
}

// Cache the stable-leaf assignment topology in event-local storage
// The proposal and coherent amplitude use exactly the same assignments
void CacheStableLeafSymmetryAssignments(
    gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &reference_tree) {
  const std::string key = StableLeafTopologyCacheKey(reference_tree);
  if (lts.decay_symmetry_topology_key != key ||
      lts.decay_symmetry_assignments.empty() ||
      lts.decay_symmetry_assignments.size() !=
          lts.decay_symmetry_statistics_signs.size()) {
    const StableLeafTopologyPlan plan =
        BuildStableLeafTopologyPlan(reference_tree);
    lts.decay_symmetry_topology_key = key;
    lts.decay_symmetry_assignments.clear();
    lts.decay_symmetry_statistics_signs.clear();
    lts.decay_symmetry_assignments.reserve(plan.assignments.size());
    lts.decay_symmetry_statistics_signs.reserve(plan.assignments.size());
    for (const auto &assignment : plan.assignments) {
      lts.decay_symmetry_assignments.push_back(assignment.assignment);
      lts.decay_symmetry_statistics_signs.push_back(assignment.statistics_sign);
    }
  }
}

// Apply one leaf assignment using a previously collected ordered slot basis
// The assignment entry at each destination selects the source four-momentum
bool ApplyStableLeafAssignment(std::vector<MDecayBranch> &tree,
                               const std::vector<std::size_t> &assignment,
                               const std::vector<StableLeafSlot> &slots) {
  if (assignment.size() != slots.size()) {
    return false;
  }
  std::vector<bool> used(slots.size(), false);
  for (const auto &destination : indices(assignment)) {
    const std::size_t source = assignment[destination];
    if (source >= slots.size() || used[source] ||
        slots[source].pdg != slots[destination].pdg ||
        slots[source].spinX2 != slots[destination].spinX2) {
      return false;
    }
    used[source] = true;
  }
  for (const auto &destination : indices(slots)) {
    BranchAtPath(tree, slots[destination].path).p4 =
        slots[assignment[destination]].p4;
  }
  RecomputeTreeMomenta(tree);
  return true;
}

// Build all stable-leaf permuted decay trees needed for coherent symmetrization
// Each tree carries a statistics sign and recomputed off-shell masses
//
std::vector<StableLeafTerm>
StableLeafSymmetryTrees(gra::LORENTZSCALAR &lts,
                        const std::vector<MDecayBranch> &reference_tree) {
  std::vector<StableLeafSlot> slots;
  CollectStableLeafSlots(reference_tree, slots);
  if (slots.size() <= 1) {
    return {{reference_tree, 1.0}};
  }

  CacheStableLeafSymmetryAssignments(lts, reference_tree);
  if (lts.decay_symmetry_assignments.size() == 1) {
    return {{reference_tree, 1.0}};
  }

  std::vector<StableLeafTerm> terms;
  terms.reserve(lts.decay_symmetry_assignments.size());
  for (const auto &term : indices(lts.decay_symmetry_assignments)) {
    std::vector<MDecayBranch> tree = reference_tree;
    if (!ApplyStableLeafAssignment(tree, lts.decay_symmetry_assignments[term],
                                   slots)) {
      throw std::invalid_argument(
          "Stable-leaf symmetrization: cached assignment is invalid");
    }
    terms.push_back(
        {std::move(tree), lts.decay_symmetry_statistics_signs[term],
         lts.decay_symmetry_assignments[term]});
  }

  return terms;
}

// Compute the product of physical couplings for all explicitly generated
// cascade vertices
std::complex<double> InternalDecayCouplingProduct(const MDecayBranch &branch) {
  if (branch.legs.empty()) {
    return 1.0;
  }

  std::complex<double> out = branch.hel.g_decay;
  for (const auto &leg : branch.legs) {
    out *= InternalDecayCouplingProduct(leg);
  }
  return out;
}

// Compute the physical internal decay-coupling product for a full tree
std::complex<double>
InternalDecayCouplingProduct(const std::vector<MDecayBranch> &tree) {
  std::complex<double> out = 1.0;
  for (const auto &branch : tree) {
    out *= InternalDecayCouplingProduct(branch);
  }
  return out;
}

// Build the full resonance decay matrix for a two-body sequential cascade
// Sequential amplitudes are tensor-contracted from the leaves to X
//
MMatrix<std::complex<double>> ResonanceCascadeDecayMatrix(
    const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res,
    const gra::HELMatrix &root_hel, const std::vector<MDecayBranch> &tree,
    const std::string &frame, const bool cascade = true) {

  const FrameTransform to_X_frame =
      RootDecayFrameTransform(lts, frame, "gra::spin::DecayAmp");

  const gra::M4Vec left_in_X = to_X_frame(tree[0].p4);
  const MMatrix<std::complex<double>> root =
      fDecayMatrix(root_hel, left_in_X.Theta(), left_in_X.Phi());
  if (!cascade || HasStableDaughterPair(tree[0], tree[1])) {
    return root.Transpose();
  }
  const auto daughter_axis = SpinHalfRotation(left_in_X.Theta(), left_in_X.Phi());
  const MMatrix<std::complex<double>> flip = {{0.0, 1.0}, {1.0, 0.0}};
  const MMatrix<std::complex<double>> left =
      CascadeBranchOperator(tree[0], to_X_frame, daughter_axis);
  const MMatrix<std::complex<double>> right =
      CascadeBranchOperator(tree[1], to_X_frame, daughter_axis * flip);
  // Daughter decay operators are noncommuting maps in ordered helicity spaces
  return left.KroneckerMultiply(right, root).Transpose();
}

// Build the continuum decay matrix for a two-body sequential cascade
// Continuum root spin space is the direct final-pair helicity space
//
MMatrix<std::complex<double>>
ContinuumCascadeDecayMatrix(const gra::LORENTZSCALAR &lts,
                            const std::vector<MDecayBranch> &tree,
                            std::size_t n_top, const std::string &frame) {

  const FrameTransform to_X_frame =
      RootDecayFrameTransform(lts, frame, "gra::spin::ContinuumDecayMatrix");

  const MMatrix<std::complex<double>> root =
      MMatrix<std::complex<double>>::IdentityMatrix(n_top);
  if (HasStableDaughterPair(tree[0], tree[1])) {
    return root.Transpose();
  }
  const auto left_in_X = to_X_frame(tree[0].p4);
  const auto daughter_axis = SpinHalfRotation(left_in_X.Theta(), left_in_X.Phi());
  const MMatrix<std::complex<double>> flip = {{0.0, 1.0}, {1.0, 0.0}};
  const MMatrix<std::complex<double>> left =
      CascadeBranchOperator(tree[0], to_X_frame, daughter_axis);
  const MMatrix<std::complex<double>> right =
      CascadeBranchOperator(tree[1], to_X_frame, daughter_axis * flip);
  // Daughter decay operators are noncommuting maps in ordered helicity spaces
  return left.KroneckerMultiply(right, root).Transpose();
}

// Match one external spin to canonical massive or helicity massless states in the root frame
// [REFERENCE: Marangotto, Adv. High Energy Phys. 2020 (2020) 6674595, arXiv:1911.10025]
MMatrix<std::complex<double>> LeafSpinFrame(
    const MDecayBranch &leaf, const FrameTransform &to_local,
    const FrameTransform &to_root, const MMatrix<std::complex<double>> &from_local,
    const MMatrix<std::complex<double>> &axis) {
  const auto helicities = FinalStateHelicities(leaf.p, "LeafSpinFrame");
  if (leaf.p.spinX2 == 0) { return {{1.0}}; }
  const M4Vec local = to_local(leaf.p4);
  const M4Vec root = to_root(leaf.p4);
  if (std::fpclassify(leaf.p.mass) == FP_ZERO) {
    const auto little = SpinHalfRotation(root.Theta(), root.Phi()).Dagger() * from_local * axis;
    MMatrix<std::complex<double>> out(helicities.size(), helicities.size(), 0.0);
    for (const auto &i : indices(helicities)) {
      const std::size_t spinor = helicities[i] < 0.0 ? 0 : 1;
      const double norm = std::abs(little[spinor][spinor]);
      if (!(norm > 0.0) || !std::isfinite(norm)) {
        throw AmplitudeFailure("LeafSpinFrame: undefined massless helicity phase");
      }
      const int power = static_cast<int>(std::llround(2.0 * std::abs(helicities[i])));
      out[i][i] = std::pow(little[spinor][spinor] / norm, power);
    }
    return out;
  }
  const auto rotation = SpinHalfWigner(from_local, local) * axis;
  return SpinRotation(rotation, 0.5 * leaf.p.spinX2);
}

// Follow the same helicity boosts as the decay contraction and retain their spinor signs
void CollectLeafSpinFrames(
    const std::vector<MDecayBranch> &tree, const FrameTransform &to_local,
    const FrameTransform &to_root, const MMatrix<std::complex<double>> &from_local,
    std::vector<MMatrix<std::complex<double>>> &frames) {
  const M4Vec left = to_local(tree[0].p4);
  const auto first_axis = SpinHalfRotation(left.Theta(), left.Phi());
  // The second Jacob-Wick helicity label denotes the opposite spin projection
  const MMatrix<std::complex<double>> flip = {{0.0, 1.0}, {1.0, 0.0}};
  for (const auto &i : indices(tree)) {
    const auto &branch = tree[i];
    if (branch.legs.empty()) {
      frames.push_back(LeafSpinFrame(branch, to_local, to_root, from_local,
                                     i == 0 ? first_axis : first_axis * flip));
      continue;
    }
    const M4Vec axis = to_local(branch.p4);
    const M4Vec system = to_local(branch.legs[0].p4 + branch.legs[1].p4);
    const auto inverse_rotation = SpinHalfRotation(axis.Theta(), axis.Phi()) *
                                  SpinHalfRotation(0.0, -gra::math::PI);
    const auto from_branch = from_local * inverse_rotation *
                             SpinHalfBoost(RotateToHelicityAxis(system, axis));
    const FrameTransform to_branch = [to_local, axis, system](const M4Vec &p) {
      return ApplyHelicityStep(to_local(p), axis, system);
    };
    CollectLeafSpinFrames(branch.legs, to_branch, to_root, from_branch, frames);
  }
}

// Compute every external spin frame in the ordered stable-leaf tensor basis
std::vector<MMatrix<std::complex<double>>> LeafSpinFrames(
    const std::vector<MDecayBranch> &tree, const FrameTransform &to_root) {
  std::vector<MMatrix<std::complex<double>>> frames;
  CollectLeafSpinFrames(tree, to_root, to_root,
                        MMatrix<std::complex<double>>::IdentityMatrix(2), frames);
  return frames;
}

// Rotate each final spin into the reference chain and permute the external helicity slots
MMatrix<std::complex<double>> MatchLeafSpins(
    MMatrix<std::complex<double>> matrix, const StableLeafTerm &term,
    const FrameTransform &to_root,
    const std::vector<MMatrix<std::complex<double>>> &reference) {
  if (reference.empty()) { return matrix; }
  const auto frames = LeafSpinFrames(term.tree, to_root);
  std::size_t before = 1;
  std::vector<std::size_t> strides(frames.size(), 1);
  for (std::size_t i = frames.size(); i > 1; --i) {
    strides[i - 2] = strides[i - 1] * reference[i - 1].size_row();
  }
  for (const auto &i : indices(frames)) {
    const std::size_t dim = frames[i].size_row();
    if (dim > 1) {
      const auto rotation = reference[term.assignment[i]].Dagger() * frames[i];
      matrix = matrix.TransformColumnAxis(rotation, before, matrix.size_col() / before / dim);
    }
    before *= dim;
  }
  std::vector<std::size_t> columns(matrix.size_col(), 0);
  for (const auto &column : indices(columns)) {
    for (const auto &i : indices(frames)) {
      const std::size_t source = term.assignment[i];
      const std::size_t helicity = column / strides[source] % reference[source].size_row();
      columns[column] = columns[column] * frames[i].size_row() + helicity;
    }
  }
  return matrix.SelectColumns(columns);
}

// Sum decay matrices over coherent stable-leaf Bose/Fermi permutations
// Indistinguishable final states interfere at amplitude level
//
template <class MatrixBuilder>
MMatrix<std::complex<double>> StableLeafSymmetrizedDecayMatrix(
    gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &reference_tree,
    const std::string &frame, MatrixBuilder build_matrix) {
  const std::vector<StableLeafTerm> terms =
      StableLeafSymmetryTrees(lts, reference_tree);
  if (terms.size() == 1) {
    return build_matrix(terms[0].tree) *
           CascadeBWProduct(terms[0].tree);
  }

  const auto to_root = RootDecayFrameTransform(lts, frame, "StableLeafSymmetrizedDecayMatrix");
  const bool scalar_leaves = std::all_of(reference_tree.begin(), reference_tree.end(),
      [](const MDecayBranch &branch) { return LeafHelicityDimension(branch) == 1; });
  const auto reference_spins = scalar_leaves ? std::vector<MMatrix<std::complex<double>>>{}
                                             : LeafSpinFrames(reference_tree, to_root);

  MMatrix<std::complex<double>> out;
  bool initialized = false;
  for (const auto &term : terms) {
    const auto amplitude = MatchLeafSpins(build_matrix(term.tree), term, to_root, reference_spins) *
                           (term.statistics_sign * CascadeBWProduct(term.tree));
    if (!initialized) {
      out = amplitude;
      initialized = true;
    } else {
      out = out + amplitude;
    }
  }
  return out;
}

} // namespace

// Cache the unique Bose and Fermi stable-leaf assignments for one topology
void PrepareStableLeafSymmetryAssignments(
    gra::LORENTZSCALAR &lts, const std::vector<MDecayBranch> &reference_tree) {
  CacheStableLeafSymmetryAssignments(lts, reference_tree);
}

// Apply one cached assignment and reconstruct every internal four-momentum
bool ApplyStableLeafSymmetryAssignment(
    std::vector<MDecayBranch> &tree,
    const std::vector<std::size_t> &assignment) {
  std::vector<StableLeafSlot> slots;
  CollectStableLeafSlots(tree, slots);
  return ApplyStableLeafAssignment(tree, assignment, slots);
}

// Compute the product of internal fixed-width Breit-Wigner factors for one branch
// Coherent cascade permutations differ by internal BW phases
std::complex<double> CascadeBWProduct(const MDecayBranch &branch) {
  if (branch.legs.empty()) {
    return 1.0;
  }

  std::complex<double> out = 1.0;
  if (branch.p.width > 0.0) {
    out *= gra::resonance::FixedWidthLineShape(branch.p4.M2(), branch.p.mass,
                                               branch.p.width);
  }
  for (const auto &leg : branch.legs) { out *= CascadeBWProduct(leg); }
  return out;
}

// Compute the product of internal fixed-width Breit-Wigner factors for a tree
// The tree BW product weights one coherent decay history
std::complex<double> CascadeBWProduct(const std::vector<MDecayBranch> &tree) {
  std::complex<double> out = 1.0;
  for (const auto &branch : tree) { out *= CascadeBWProduct(branch); }
  return out;
}

// Compute stable-leaf terms for non-Jacob-Wick coherent amplitude builders
// External amplitude code receives the same Bose/Fermi assignments
std::vector<StableLeafAmplitudeTerm>
StableLeafAmplitudeTrees(gra::LORENTZSCALAR &lts,
                         const std::vector<MDecayBranch> &reference_tree) {
  const std::vector<StableLeafTerm> terms =
      StableLeafSymmetryTrees(lts, reference_tree);
  std::vector<StableLeafAmplitudeTerm> out;
  out.reserve(terms.size());
  for (const auto &term : terms) {
    out.push_back({term.tree, term.statistics_sign});
  }
  return out;
}

// Build a recursive two-body decay chain with spin correlations to arbitrary depth
// X decay spin is propagated through all sequential helicity frames
//
// X -> A B
// X -> A > {A1 A2} B > {B1 B2}
// ..
// X -> A > {A1 > {C1 C2} A2} B > ...
//
MMatrix<std::complex<double>>
ResonanceDecayMatrix(gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res,
                     const gra::HELMatrix &root_hel, const std::string &frame, const bool cascade) {

  // Retain the ordered root daughter helicities before applying subsequent decays
  if (!cascade) { return ResonanceCascadeDecayMatrix(lts, res, root_hel, lts.decaytree, frame, false); }

  // Isolated decays are kinematic proposals and must not contain decay spin
  // algebra
  if (lts.process.root_decay_mode == gra::RootDecayMode::Isolated) {
    const unsigned int rows = static_cast<unsigned int>(
        SpinStateCount(0.5 * static_cast<double>(res.p.spinX2),
                       "gra::spin::DecayAmp isolated resonance"));
    const unsigned int cols =
        static_cast<unsigned int>(FinalStateHelicityDimension(lts));
    return SpinBlindDecayMatrix(rows, cols);
  }

  if (lts.process.SPINDEC == false) {
    const unsigned int rows = static_cast<unsigned int>(
        SpinStateCount(0.5 * static_cast<double>(res.p.spinX2),
                       "gra::spin::DecayAmp resonance"));
    const unsigned int cols =
        static_cast<unsigned int>(FinalStateHelicityDimension(lts));
    // Turning off angular correlations preserves physical couplings and propagators
    return SpinBlindDecayMatrix(rows, cols) *
           InternalDecayCouplingProduct(lts.decaytree) *
           CascadeBWProduct(lts.decaytree);
  }

  // For direct multi-body root decays, fall back to a spin-blind final-state
  // operator instead of leaving decay_f empty
  if (lts.decaytree.size() != 2) {
    const unsigned int rows = static_cast<unsigned int>(
        SpinStateCount(0.5 * static_cast<double>(res.p.spinX2),
                       "gra::spin::DecayAmp resonance"));
    const unsigned int cols =
        static_cast<unsigned int>(FinalStateHelicityDimension(lts));
    return SpinBlindDecayMatrix(rows, cols) *
           InternalDecayCouplingProduct(lts.decaytree) *
           CascadeBWProduct(lts.decaytree);
  }

  auto build_matrix = [&lts, &res, &root_hel,
                       &frame](const std::vector<MDecayBranch> &tree) {
    return ResonanceCascadeDecayMatrix(lts, res, root_hel, tree, frame);
  };

  const bool decay_sym =
      lts.amplitude.DECAY_SYM &&
      lts.process.root_decay_mode != gra::RootDecayMode::Isolated;
  if (decay_sym) {
    return StableLeafSymmetrizedDecayMatrix(lts, lts.decaytree, frame, build_matrix);
  }
  return build_matrix(lts.decaytree) *
         CascadeBWProduct(lts.decaytree);
}

// Construct an event-dependent resonance decay operator
MMatrix<std::complex<double>> ResonanceDecayMatrix(gra::LORENTZSCALAR &lts,
                                                   const gra::PARAM_RES &res,
                                                   const std::string &frame) {
  return ResonanceDecayMatrix(lts, res, res.hel_decay, frame);
}

// Store the event-dependent decay operator in a mutable resonance instance
void DecayAmp(gra::LORENTZSCALAR &lts, gra::PARAM_RES &res,
              const std::string &frame) {
  res.decay_f = ResonanceDecayMatrix(lts, res, frame);
}

// Build or return the cached continuum decay matrix for the current event
// Continuum spin decay acts on the final-pair helicity basis
//
MMatrix<std::complex<double>> ContinuumDecayMatrix(gra::LORENTZSCALAR &lts,
                                                   const std::string &frame) {

  // Reuse only a decay map bound to the current Born event during screening
  if (lts.screening.active && lts.amplitude.central.Active()) {
    return lts.amplitude.continuum_decay.Get(lts.amplitude.central, frame,
                                             "ContinuumDecayMatrix");
  }

  const std::size_t n_left = FinalStateHelicityCount(
      lts.decaytree[0].p, "ContinuumDecayMatrix left branch");
  const std::size_t n_right = FinalStateHelicityCount(
      lts.decaytree[1].p, "ContinuumDecayMatrix right branch");
  const std::size_t n_top = n_left * n_right;

  if (!lts.process.SPINDEC) {
    // Keep physical cascade normalization while replacing angular matrices by a
    // flat map
    MMatrix<std::complex<double>> out =
        SpinBlindBranch(lts.decaytree[0]).Kronecker(SpinBlindBranch(lts.decaytree[1])) *
        InternalDecayCouplingProduct(lts.decaytree) *
        CascadeBWProduct(lts.decaytree);
    if (!lts.amplitude.central.Active()) {
      return out;
    }
    return lts.amplitude.continuum_decay.Store(
        lts.amplitude.central, frame, std::move(out), "ContinuumDecayMatrix");
  }

  auto build_matrix = [&lts, n_top,
                       &frame](const std::vector<MDecayBranch> &tree) {
    return ContinuumCascadeDecayMatrix(lts, tree, n_top, frame);
  };

  const bool decay_sym =
      lts.amplitude.DECAY_SYM &&
      lts.process.root_decay_mode != gra::RootDecayMode::Isolated;
  MMatrix<std::complex<double>> out =
      decay_sym
          ? StableLeafSymmetrizedDecayMatrix(lts, lts.decaytree, frame, build_matrix)
          : build_matrix(lts.decaytree) *
                CascadeBWProduct(lts.decaytree);

  if (!lts.amplitude.central.Active()) {
    return out;
  }
  return lts.amplitude.continuum_decay.Store(
      lts.amplitude.central, frame, std::move(out), "ContinuumDecayMatrix");
}

} // namespace spin
} // namespace gra
