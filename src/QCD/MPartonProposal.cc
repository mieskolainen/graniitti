// Generic final-state parton proposal operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <exception>
#include <limits>
#include <map>
#include <optional>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"

namespace gra::parton {

using gra::aux::indices;

namespace {

// Collect mutable stable decay leaves recursively
void CollectMutableStableLeaves(MDecayBranch &branch,
                                std::vector<MDecayBranch *> &leaves) {
  if (branch.legs.empty()) {
    leaves.push_back(&branch);
    return;
  }
  for (auto &leg : branch.legs) {
    CollectMutableStableLeaves(leg, leaves);
  }
}

// Collect mutable stable decay leaves from the central tree
std::vector<MDecayBranch *>
CollectMutableStableLeaves(std::vector<MDecayBranch> &tree) {
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : tree) {
    CollectMutableStableLeaves(branch, leaves);
  }
  return leaves;
}

// Stable leaf location with its immediate sibling group
struct StableLeafSlot {
  std::size_t leaf_index = 0;
  std::size_t sibling_group = 0;
  std::size_t sibling_index = 0;
  int pdg = 0;
};

// Collect stable leaf locations while preserving immediate parentage
void CollectStableLeafSlots(const std::vector<MDecayBranch> &siblings,
                            const std::size_t sibling_group,
                            std::size_t &next_sibling_group,
                            std::size_t &next_leaf_index,
                            std::vector<StableLeafSlot> &slots) {
  for (std::size_t i = 0; i < siblings.size(); ++i) {
    const auto &branch = siblings[i];
    if (branch.legs.empty()) {
      slots.push_back({next_leaf_index++, sibling_group, i, branch.p.pdg});
      continue;
    }
    const std::size_t child_group = next_sibling_group++;
    CollectStableLeafSlots(branch.legs, child_group, next_sibling_group,
                           next_leaf_index, slots);
  }
}

// Collect stable leaf locations from the full central decay tree
std::vector<StableLeafSlot>
CollectStableLeafSlots(const std::vector<MDecayBranch> &tree) {
  std::vector<StableLeafSlot> slots;
  std::size_t next_sibling_group = 1;
  std::size_t next_leaf_index = 0;
  CollectStableLeafSlots(tree, 0, next_sibling_group, next_leaf_index, slots);
  return slots;
}

// One correlated group of generic parton leaves and its physical choices
struct ProposalGroup {
  std::vector<std::size_t> leaf_indices;
  std::vector<std::vector<int>> choices;
};

// Compute one physical PDG id with the sign carried by a generic alias
int SignedPartonPDG(const int species, const int alias) {
  if (species == PDG::PDG_gluon) {
    return species;
  }
  return alias > 0 ? species : -species;
}

// Build independent or particle-antiparticle correlated alias groups
std::vector<ProposalGroup>
BuildProposalGroups(const std::vector<StableLeafSlot> &leaves,
                    const std::vector<int> &species) {
  std::vector<ProposalGroup> groups;
  std::size_t leaf_cursor = 0;
  while (leaf_cursor < leaves.size()) {
    const std::size_t i = leaf_cursor;
    const int alias = leaves[i].pdg;
    if (std::abs(alias) != PDG::PDG_hard_jet) {
      ++leaf_cursor;
      continue;
    }

    const bool correlated =
        i + 1 < leaves.size() &&
        leaves[i].sibling_group == leaves[i + 1].sibling_group &&
        leaves[i].sibling_index + 1 == leaves[i + 1].sibling_index &&
        std::abs(leaves[i + 1].pdg) == PDG::PDG_hard_jet &&
        alias == -leaves[i + 1].pdg;
    ProposalGroup group;
    group.leaf_indices =
        correlated ? std::vector<std::size_t>{leaves[i].leaf_index,
                                              leaves[i + 1].leaf_index}
                   : std::vector<std::size_t>{leaves[i].leaf_index};
    for (const int flavour : species) {
      if (correlated) {
        std::vector<int> choice = {SignedPartonPDG(flavour, alias)};
        choice.push_back(SignedPartonPDG(flavour, leaves[i + 1].pdg));
        group.choices.push_back(std::move(choice));
      } else if (flavour == PDG::PDG_gluon) {
        group.choices.push_back({flavour});
      } else {
        group.choices.push_back({flavour});
        group.choices.push_back({-flavour});
      }
    }
    groups.push_back(std::move(group));
    leaf_cursor += correlated ? 2 : 1;
  }
  return groups;
}

// Compute the Cartesian physical mode count with overflow checking
std::size_t ProposalModeCount(const std::vector<ProposalGroup> &groups) {
  std::size_t count = 1;
  for (const auto &group : groups) {
    if (group.choices.empty() ||
        count >
            std::numeric_limits<std::size_t>::max() / group.choices.size()) {
      throw PhaseSpaceFailure(
          "parton::PrepareFinalStateProposal: invalid proposal mode count");
    }
    count *= group.choices.size();
  }
  return count;
}

// Compute true when one generated node can represent the generic branch
bool ProcessNodeMatchesTemplate(const AmplitudeTopologyNode &node,
                                const MDecayBranch &branch) {
  if (node.daughters.size() != branch.legs.size()) {
    return false;
  }
  if (std::abs(branch.p.pdg) == PDG::PDG_hard_jet) {
    if (!branch.legs.empty()) {
      return false;
    }
    return std::any_of(node.allowed_pdgs.begin(), node.allowed_pdgs.end(),
                       [](const int pdg) {
                         return pdg == PDG::PDG_gluon ||
                                (std::abs(pdg) >= 1 && std::abs(pdg) <= 6) ||
                                std::abs(pdg) == PDG::PDG_hard_jet;
                       });
  }
  if (std::find(node.allowed_pdgs.begin(), node.allowed_pdgs.end(),
                branch.p.pdg) == node.allowed_pdgs.end()) {
    return false;
  }
  for (std::size_t i = 0; i < branch.legs.size(); ++i) {
    if (!ProcessNodeMatchesTemplate(node.daughters[i], branch.legs[i])) {
      return false;
    }
  }
  return true;
}

// Compute true when one generated topology can represent the generic tree
bool ProcessMatchesTemplate(const AmplitudeTopology &topology,
                            const std::vector<MDecayBranch> &template_tree) {
  if (topology.size() != template_tree.size()) {
    return false;
  }
  for (std::size_t i = 0; i < template_tree.size(); ++i) {
    if (!ProcessNodeMatchesTemplate(topology[i], template_tree[i])) {
      return false;
    }
  }
  return true;
}

// Append exact stable PDGs from one generated topology branch
void AppendExactStablePDGs(const AmplitudeTopologyNode &node,
                           std::vector<int> &pdgs) {
  if (node.daughters.empty()) {
    if (node.allowed_pdgs.size() != 1) {
      throw PhaseSpaceFailure(
          "parton::PrepareFinalStateProposal: generated topology is not exact");
    }
    pdgs.push_back(node.allowed_pdgs.front());
    return;
  }
  for (const auto &daughter : node.daughters) {
    AppendExactStablePDGs(daughter, pdgs);
  }
}

// Compute true when one exact stable mode obeys all alias choices
bool ExactModeAllowed(const std::vector<int> &mode,
                      const std::vector<ProposalGroup> &groups) {
  for (const auto &group : groups) {
    std::vector<int> choice;
    choice.reserve(group.leaf_indices.size());
    for (const std::size_t index : group.leaf_indices) {
      if (index >= mode.size()) {
        return false;
      }
      choice.push_back(mode[index]);
    }
    if (std::find(group.choices.begin(), group.choices.end(), choice) ==
        group.choices.end()) {
      return false;
    }
  }
  return true;
}

// Collect exact generated modes compatible with one generic alias template
std::pair<bool, std::vector<std::vector<int>>>
GeneratedProposalModes(const std::vector<amplitude::Process> &processes,
                       const std::vector<MDecayBranch> &template_tree,
                       const std::vector<ProposalGroup> &groups,
                       const std::size_t stable_leaf_count) {
  bool constrained = false;
  std::vector<std::vector<int>> modes;
  for (const auto &process_record : processes) {
    if (!ProcessMatchesTemplate(process_record.topology, template_tree)) {
      continue;
    }
    std::vector<AmplitudeTopology> exact_topologies =
        process_record.channel_topologies;
    if (exact_topologies.empty() &&
        process_record.matrix_element_form ==
            amplitude::MatrixElementForm::Generated &&
        process_record.topology_mode == amplitude::TopologyMode::Exact) {
      exact_topologies.push_back(process_record.topology);
    }
    if (exact_topologies.empty()) {
      continue;
    }
    constrained = true;
    for (const auto &topology : exact_topologies) {
      if (!ProcessMatchesTemplate(topology, template_tree)) {
        continue;
      }
      std::vector<int> mode;
      for (const auto &node : topology) {
        AppendExactStablePDGs(node, mode);
      }
      if (mode.size() != stable_leaf_count || !ExactModeAllowed(mode, groups)) {
        continue;
      }
      if (std::find(modes.begin(), modes.end(), mode) == modes.end()) {
        modes.push_back(std::move(mode));
      }
    }
  }
  return {constrained, modes};
}

// Collect Cartesian alias modes under an optional amplitude constraint
std::vector<std::vector<int>>
AmplitudeProposalModes(const std::vector<MDecayBranch> &template_tree,
                       const std::vector<StableLeafSlot> &template_leaves,
                       const std::vector<ProposalGroup> &groups,
                       MSubProc &subprocess) {
  const std::size_t count = ProposalModeCount(groups);
  bool constrained = false;
  std::vector<std::vector<int>> cartesian_modes;
  std::vector<std::vector<int>> modes;
  for (std::size_t index = 0; index < count; ++index) {
    std::vector<int> mode;
    mode.reserve(template_leaves.size());
    for (const auto &leaf : template_leaves) {
      mode.push_back(leaf.pdg);
    }

    std::size_t cursor = index;
    for (const auto &group : groups) {
      const auto &choice = group.choices[cursor % group.choices.size()];
      cursor /= group.choices.size();
      for (std::size_t i = 0; i < group.leaf_indices.size(); ++i) {
        mode[group.leaf_indices[i]] = choice[i];
      }
    }

    std::vector<MDecayBranch> candidate = template_tree;
    const auto leaves = CollectMutableStableLeaves(candidate);
    if (leaves.size() != mode.size()) {
      throw std::logic_error(
          "parton::AmplitudeProposalModes: stable topology changed");
    }
    for (std::size_t i = 0; i < leaves.size(); ++i) {
      leaves[i]->p.pdg = mode[i];
    }
    const auto accepted = subprocess.AcceptsPartonMode(candidate);
    cartesian_modes.push_back(mode);
    if (accepted.has_value()) {
      constrained = true;
      if (*accepted) { modes.push_back(std::move(mode)); }
    }
  }
  return constrained ? modes : cartesian_modes;
}

// Assign one physical parton and its process-local mass to a stable leaf
void AssignFinalStateParton(MDecayBranch &leaf, const int pdg,
                            const LORENTZSCALAR &lts) {
  MParticle particle;
  try {
    particle = lts.PDG.FindByPDG(pdg);
  } catch (const std::exception &) {
    throw PhaseSpaceFailure(
        "parton::PrepareFinalStateProposal: selected parton is absent from "
        "the process PDG table");
  }
  const auto mass = lts.final_state_parton_masses.find(std::abs(pdg));
  if (mass != lts.final_state_parton_masses.end()) {
    particle.mass = mass->second;
  }
  leaf.p = particle;
}

// Store one model and its outgoing quark species during initialization
struct MG5MassData {
  mg5::ParticleMap particles;
  std::set<int> flavours;
};

// Collect evaluated models and quark species from the allowed physical modes
std::map<std::string, MG5MassData> MG5MassCards(const LORENTZSCALAR &lts, MSubProc &subprocess) {
  LORENTZSCALAR proposal = lts;
  ConfigureFinalStateProposal(proposal, subprocess);
  if (proposal.final_state_parton_modes.empty()) {
    std::vector<int> mode;
    for (const auto &leaf : CollectStableLeafSlots(lts.decaytree)) { mode.push_back(leaf.pdg); }
    proposal.final_state_parton_modes.push_back(std::move(mode));
  }

  auto                                 tree   = lts.decaytree;
  const auto                           leaves = CollectMutableStableLeaves(tree);
  std::map<std::string, MG5MassData> cards;
  for (const auto &mode : proposal.final_state_parton_modes) {
    for (const auto &i : indices(leaves)) { leaves[i]->p.pdg = mode[i]; }
    const auto process = subprocess.MatchProcess(tree);
    if (!process.has_value()) { continue; }
    const auto card = amplitude::ParameterCard(process->process_family, process->process_name);
    if (!card.has_value()) { continue; }
    auto [entry, inserted] = cards.try_emplace(*card);
    if (inserted) { entry->second.particles = MG5Particles(*process); }
    auto &flavours = entry->second.flavours;
    for (const int pdg : mode) {
      const int flavour = std::abs(pdg);
      if (flavour >= 1 && flavour <= 5) { flavours.insert(flavour); }
    }
  }
  return cards;
}

// Compare generated mass parameters with a relative roundoff tolerance
bool SameMassParameter(double first, double second) {
  return std::abs(first - second) <=
         64.0 * std::numeric_limits<double>::epsilon() * std::max({1.0, std::abs(first), std::abs(second)});
}

// Require common pole parameters for internal branches shared by physical modes
void CheckDecayMasses(const std::vector<MDecayBranch> &first, const std::vector<MDecayBranch> &second) {
  for (const auto &i : indices(first)) {
    if (!SameMassParameter(first[i].p.mass, second[i].p.mass) ||
        !SameMassParameter(first[i].p.width, second[i].p.width)) {
      throw std::invalid_argument("parton::ConfigureFinalStateMasses: conflicting generated decay pole parameters");
    }
    CheckDecayMasses(first[i].legs, second[i].legs);
  }
}

} // namespace

// Configure generated decay parameters and outgoing quark masses
void ConfigureFinalStateMasses(LORENTZSCALAR &lts, MSubProc &subprocess) {
  lts.final_state_parton_masses.clear();
  bool decay_configured = false;
  for (const auto &[card, model] : MG5MassCards(lts, subprocess)) {
    for (const auto &[pdg, mass] : lts.particle_mass_overrides) {
      const auto found = model.particles.find(pdg);
      if (found != model.particles.end() && !SameMassParameter(mass, found->second.mass)) {
        throw std::invalid_argument("MG5: @PDG mass conflicts with " + card + " for PDG " + std::to_string(pdg));
      }
    }
    for (const auto &[pdg, width] : lts.particle_width_overrides) {
      const auto found = model.particles.find(pdg);
      if (found != model.particles.end() && !SameMassParameter(width, found->second.width)) {
        throw std::invalid_argument("MG5: @PDG width conflicts with " + card + " for PDG " + std::to_string(pdg));
      }
    }
    auto tree = lts.decaytree;
    SynchronizeMG5DecayParameters(tree, model.particles);
    if (decay_configured) { CheckDecayMasses(lts.decaytree, tree); }
    lts.decaytree    = std::move(tree);
    decay_configured = true;

    for (const int flavour : model.flavours) {
      const double mass = model.particles.at(flavour).mass;
      const auto [entry, inserted] = lts.final_state_parton_masses.emplace(flavour, mass);
      if (!inserted && !SameMassParameter(entry->second, mass)) {
        throw std::invalid_argument("parton::ConfigureFinalStateMasses: conflicting generated quark masses for PDG " +
                                    std::to_string(flavour));
      }
    }
  }
}

// Cache physical alias modes accepted by the active amplitude definition
void ConfigureFinalStateProposal(LORENTZSCALAR &lts, MSubProc &subprocess) {
  lts.final_state_parton_modes.clear();
  const auto alias_leaves = CollectStableLeafSlots(lts.decaytree);
  if (std::none_of(alias_leaves.begin(), alias_leaves.end(),
                   [](const StableLeafSlot &leaf) {
                     return std::abs(leaf.pdg) == PDG::PDG_hard_jet;
                   })) {
    return;
  }
  const std::vector<int> species = lts.final_state_partons.Flavours();
  if (species.empty()) {
    throw std::invalid_argument(
        "parton::ConfigureFinalStateProposal: @j selects no supported parton "
        "species");
  }
  const auto groups = BuildProposalGroups(alias_leaves, species);
  auto [generated_constraint, modes] = GeneratedProposalModes(
      subprocess.Processes(), lts.decaytree, groups, alias_leaves.size());
  if (!generated_constraint) {
    modes = AmplitudeProposalModes(lts.decaytree, alias_leaves, groups,
                                   subprocess);
  }
  if (modes.empty()) {
    throw std::invalid_argument(
        "parton::ConfigureFinalStateProposal: generic aliases select no "
        "amplitude final state");
  }
  lts.final_state_parton_modes = std::move(modes);
}

// Propose one physical on-shell mode for generic final-state parton aliases
// q_i = 1/N_mode for each accepted mode
bool PrepareFinalStateProposal(LORENTZSCALAR &lts, MRandom &random) {
  lts.parton_proposal_weight = 1.0;
  if (lts.generic_parton_decaytree.empty()) {
    return false;
  }

  const auto alias_leaves =
      CollectStableLeafSlots(lts.generic_parton_decaytree);
  std::vector<MDecayBranch *> leaves =
      CollectMutableStableLeaves(lts.decaytree);
  if (leaves.size() != alias_leaves.size()) {
    throw PhaseSpaceFailure(
        "parton::PrepareFinalStateProposal: decay tree topology changed after "
        "generic parton setup");
  }
  if (lts.final_state_parton_modes.empty()) {
    throw PhaseSpaceFailure(
        "parton::PrepareFinalStateProposal: no configured amplitude mode");
  }

  const std::size_t mode_count = lts.final_state_parton_modes.size();
  const double draw = random.U(0.0, static_cast<double>(mode_count));
  const std::size_t mode_index =
      std::min(static_cast<std::size_t>(draw), mode_count - 1);
  const auto &mode = lts.final_state_parton_modes[mode_index];
  if (mode.size() != leaves.size()) {
    throw PhaseSpaceFailure(
        "parton::PrepareFinalStateProposal: configured mode size changed");
  }
  for (std::size_t i = 0; i < alias_leaves.size(); ++i) {
    if (std::abs(alias_leaves[i].pdg) == PDG::PDG_hard_jet) {
      AssignFinalStateParton(*leaves[i], mode[i], lts);
    }
  }

  lts.parton_proposal_weight = static_cast<double>(mode_count);
  return true;
}

// Compute the inverse probability of the current physical parton proposal
// 1/q_i = N_mode
double ProposalWeight(const LORENTZSCALAR &lts) {
  const double weight = lts.parton_proposal_weight;
  if (!std::isfinite(weight) || weight < 1.0) {
    throw PhaseSpaceFailure(
        "parton::ProposalWeight: invalid inverse proposal probability");
  }
  return weight;
}

} // namespace gra::parton
