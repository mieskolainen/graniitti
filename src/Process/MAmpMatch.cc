// Amplitude process matching
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra {
namespace amplitude {
namespace {

// Compute one readable topology token for an event particle
std::string AmplitudeParticleToken(const MParticle &particle) {
  const int pdg = particle.pdg;
  switch (pdg) {
    case 1:
      return "d";
    case -1:
      return "d~";
    case 2:
      return "u";
    case -2:
      return "u~";
    case 3:
      return "s";
    case -3:
      return "s~";
    case 4:
      return "c";
    case -4:
      return "c~";
    case 5:
      return "b";
    case -5:
      return "b~";
    case 6:
      return "t";
    case -6:
      return "t~";
    case 11:
      return "e-";
    case -11:
      return "e+";
    case 12:
      return "ve";
    case -12:
      return "ve~";
    case 13:
      return "mu-";
    case -13:
      return "mu+";
    case 14:
      return "vm";
    case -14:
      return "vm~";
    case 15:
      return "ta-";
    case -15:
      return "ta+";
    case 16:
      return "vt";
    case -16:
      return "vt~";
    case 21:
      return "g";
    case 22:
      return "a";
    case 23:
      return "z";
    case 24:
      return "w+";
    case -24:
      return "w-";
    case 25:
      return "h";
    case PDG::PDG_hard_jet:
    case -PDG::PDG_hard_jet:
      return "j";
    default:
      return "pdg:" + std::to_string(pdg);
  }
}

// Render one event decay branch with exact daughter ordering
std::string EventTopologyBranch(const MDecayBranch &branch) {
  std::string signature = AmplitudeParticleToken(branch.p);
  if (branch.legs.empty()) { return signature; }

  signature += "(";
  for (const auto &i : indices(branch.legs)) {
    if (i != 0) { signature += ","; }
    signature += EventTopologyBranch(branch.legs[i]);
  }
  signature += ")";
  return signature;
}

// Build one exact numeric topology node from an event decay branch
AmplitudeTopologyNode EventTopologyNode(const MDecayBranch &branch) {
  AmplitudeTopologyNode node;
  node.allowed_pdgs = {branch.p.pdg};
  node.daughters.reserve(branch.legs.size());
  for (const auto &daughter : branch.legs) { node.daughters.push_back(EventTopologyNode(daughter)); }
  return node;
}

// Append stable PDG codes from one branch in event momentum order
void CollectStablePDGs(const MDecayBranch &branch, std::vector<int> &pdgs) {
  if (branch.legs.empty()) {
    pdgs.push_back(branch.p.pdg);
    return;
  }
  for (const auto &leg : branch.legs) { CollectStablePDGs(leg, pdgs); }
}

}  // namespace

// Render one readable ordered event decay tree
std::string AmplitudeTopologySignature(const std::vector<MDecayBranch> &tree) {
  std::string signature;
  for (const auto &i : indices(tree)) {
    if (i != 0) { signature += ","; }
    signature += EventTopologyBranch(tree[i]);
  }
  return signature;
}

// Compute the exact ordered numeric topology of one event decay tree
AmplitudeTopology AmplitudeTopologyFromDecayTree(const std::vector<MDecayBranch> &decaytree) {
  AmplitudeTopology topology;
  topology.reserve(decaytree.size());
  for (const auto &branch : decaytree) { topology.push_back(EventTopologyNode(branch)); }
  return topology;
}

// Compute stable final-state PDGs in event momentum order
std::vector<int> StableFinalStatePDGs(const std::vector<MDecayBranch> &decaytree) {
  std::vector<int> pdgs;
  for (const auto &branch : decaytree) { CollectStablePDGs(branch, pdgs); }
  return pdgs;
}

namespace {

// Compute true when one explicit particle set contains the event PDG
bool ParticleMatches(const std::vector<int> &allowed_pdgs, int pdg) {
  return std::find(allowed_pdgs.begin(), allowed_pdgs.end(), pdg) != allowed_pdgs.end();
}

// Compute true when one topology node matches one event branch
bool TopologyNodeMatches(const AmplitudeTopologyNode &topology, const MDecayBranch &branch) {
  if (!ParticleMatches(topology.allowed_pdgs, branch.p.pdg) || topology.daughters.size() != branch.legs.size()) {
    return false;
  }
  for (const auto &i : indices(topology.daughters)) {
    if (!TopologyNodeMatches(topology.daughters[i], branch.legs[i])) { return false; }
  }
  return true;
}

// Compute true when one ordered topology matches an event decay tree
bool TopologyMatches(const AmplitudeTopology &topology, const std::vector<MDecayBranch> &tree) {
  if (topology.size() != tree.size()) { return false; }
  for (const auto &i : indices(topology)) {
    if (!TopologyNodeMatches(topology[i], tree[i])) { return false; }
  }
  return true;
}

// Compute true when ordered stable PDGs agree with the topology leaves
bool StablePDGsMatch(const std::vector<MDecayBranch> &tree, const std::vector<int> &expected) {
  const std::vector<int> actual = StableFinalStatePDGs(tree);
  if (actual.size() != expected.size()) { return false; }
  for (const auto &i : indices(actual)) {
    if (expected[i] != 0 && expected[i] != actual[i]) { return false; }
  }
  return true;
}

// Compute true when two explicit PDG particle sets have a nonempty intersection
bool ParticleSetsOverlap(const std::vector<int> &first, const std::vector<int> &second) {
  return std::any_of(first.begin(), first.end(), [&](int pdg) { return ParticleMatches(second, pdg); });
}

// Compute true when two topology nodes can accept at least one common branch
bool TopologyNodesOverlap(const AmplitudeTopologyNode &first, const AmplitudeTopologyNode &second) {
  if (!ParticleSetsOverlap(first.allowed_pdgs, second.allowed_pdgs) ||
      first.daughters.size() != second.daughters.size()) {
    return false;
  }
  for (const auto &i : indices(first.daughters)) {
    if (!TopologyNodesOverlap(first.daughters[i], second.daughters[i])) { return false; }
  }
  return true;
}

// Compute true when two ordered topologies have a common decay tree
bool TopologiesOverlap(const AmplitudeTopology &first, const AmplitudeTopology &second) {
  if (first.size() != second.size()) { return false; }
  for (const auto &i : indices(first)) {
    if (!TopologyNodesOverlap(first[i], second[i])) { return false; }
  }
  return true;
}

// Compute true when every PDG accepted by the first particle set is in the
// second
bool ParticleSetIsSubset(const std::vector<int> &first, const std::vector<int> &second) {
  return std::all_of(first.begin(), first.end(), [&](int pdg) { return ParticleMatches(second, pdg); });
}

// Compute true when the first topology node is contained in the second
bool TopologyNodeIsSubset(const AmplitudeTopologyNode &first, const AmplitudeTopologyNode &second) {
  if (!ParticleSetIsSubset(first.allowed_pdgs, second.allowed_pdgs) ||
      first.daughters.size() != second.daughters.size()) {
    return false;
  }
  for (const auto &i : indices(first.daughters)) {
    if (!TopologyNodeIsSubset(first.daughters[i], second.daughters[i])) { return false; }
  }
  return true;
}

// Compute true when every decay tree accepted by the first is in the second
bool TopologyIsSubset(const AmplitudeTopology &first, const AmplitudeTopology &second) {
  if (first.size() != second.size()) { return false; }
  for (const auto &i : indices(first)) {
    if (!TopologyNodeIsSubset(first[i], second[i])) { return false; }
  }
  return true;
}

// Compute true when one exact topology contains an unresolved hard-jet alias
bool TopologyHasHardJetAlias(const AmplitudeTopologyNode &node) {
  if (std::find(node.allowed_pdgs.begin(), node.allowed_pdgs.end(), PDG::PDG_hard_jet) != node.allowed_pdgs.end() ||
      std::find(node.allowed_pdgs.begin(), node.allowed_pdgs.end(), -PDG::PDG_hard_jet) != node.allowed_pdgs.end()) {
    return true;
  }
  return std::any_of(node.daughters.begin(), node.daughters.end(),
                     [](const auto &daughter) { return TopologyHasHardJetAlias(daughter); });
}

// Compute true when one ordered topology contains an unresolved hard-jet alias
bool TopologyHasHardJetAlias(const AmplitudeTopology &topology) {
  return std::any_of(topology.begin(), topology.end(), [](const auto &node) { return TopologyHasHardJetAlias(node); });
}

// Compute true when one exact node is compatible with a partly resolved request
bool ExactTopologyNodeMatchesRequest(const AmplitudeTopologyNode &exact, const AmplitudeTopologyNode &request) {
  if (exact.allowed_pdgs.size() != 1 || request.allowed_pdgs.size() != 1 ||
      exact.daughters.size() != request.daughters.size()) {
    return false;
  }
  const int requested_pdg = request.allowed_pdgs.front();
  if (std::abs(requested_pdg) == PDG::PDG_hard_jet) {
    const int exact_pdg = exact.allowed_pdgs.front();
    return request.daughters.empty() && exact.daughters.empty() &&
           (exact_pdg == PDG::PDG_gluon || (std::abs(exact_pdg) >= 1 && std::abs(exact_pdg) <= 6));
  }
  if (exact.allowed_pdgs.front() != requested_pdg) { return false; }
  for (const auto &index : indices(exact.daughters)) {
    if (!ExactTopologyNodeMatchesRequest(exact.daughters[index], request.daughters[index])) { return false; }
  }
  return true;
}

// Compute true when one exact topology is compatible with a partial request
bool ExactTopologyMatchesRequest(const AmplitudeTopology &exact, const AmplitudeTopology &request) {
  if (exact.size() != request.size()) { return false; }
  for (const auto &index : indices(exact)) {
    if (!ExactTopologyNodeMatchesRequest(exact[index], request[index])) { return false; }
  }
  return true;
}

// Compute true when one constrained generated process accepts the event tree
bool GeneratedTopologyMatches(const Process &process, const std::vector<MDecayBranch> &tree) {
  if (process.channel_topologies.empty()) { return true; }
  const AmplitudeTopology event = AmplitudeTopologyFromDecayTree(tree);
  return std::any_of(process.channel_topologies.begin(), process.channel_topologies.end(),
                     [&](const auto &exact) { return ExactTopologyMatchesRequest(exact, event); });
}

// Validate one generated topology node and append its stable particle set form
void ValidateTopologyNode(const AmplitudeTopologyNode &node, std::vector<int> &stable_pdgs, bool &has_particle_set) {
  if (node.allowed_pdgs.empty() ||
      std::find(node.allowed_pdgs.begin(), node.allowed_pdgs.end(), 0) != node.allowed_pdgs.end() ||
      !std::is_sorted(node.allowed_pdgs.begin(), node.allowed_pdgs.end()) ||
      std::adjacent_find(node.allowed_pdgs.begin(), node.allowed_pdgs.end()) != node.allowed_pdgs.end()) {
    throw std::invalid_argument(
        "Generated amplitude topology particle sets must be nonempty, sorted, "
        "unique and nonzero");
  }
  has_particle_set = has_particle_set || node.allowed_pdgs.size() > 1;
  if (node.daughters.empty()) {
    stable_pdgs.push_back(node.allowed_pdgs.size() == 1 ? node.allowed_pdgs.front() : 0);
    return;
  }
  for (const auto &daughter : node.daughters) { ValidateTopologyNode(daughter, stable_pdgs, has_particle_set); }
}

// Compute true when one registered process is narrower than another
bool ProcessIsStrictSubset(const Process &first, const Process &second) {
  return first.topology != second.topology && TopologyIsSubset(first.topology, second.topology);
}

// Store one decay tree in canonical MG5 external-particle order
struct MG5Order {
  AmplitudeTopology         topology;
  std::vector<MDecayBranch> tree;
};

// Compute true when one branch contains an unresolved final-state parton alias
bool HasMG5PartonAlias(const MDecayBranch &branch) {
  return std::abs(branch.p.pdg) == PDG::PDG_hard_jet ||
         std::any_of(branch.legs.begin(), branch.legs.end(), [](const auto &leg) { return HasMG5PartonAlias(leg); });
}

// Compute true when one tree contains an unresolved final-state parton alias
bool HasMG5PartonAlias(const std::vector<MDecayBranch> &tree) {
  return std::any_of(tree.begin(), tree.end(), [](const auto &branch) { return HasMG5PartonAlias(branch); });
}

// Match one user particle to an MG5 external particle or parton slot
bool MG5ParticleMatches(const AmplitudeTopologyNode &node, const MDecayBranch &branch) {
  if (std::abs(branch.p.pdg) != PDG::PDG_hard_jet) { return ParticleMatches(node.allowed_pdgs, branch.p.pdg); }
  if (!branch.legs.empty()) { return false; }
  return std::any_of(node.allowed_pdgs.begin(), node.allowed_pdgs.end(), [&](const int pdg) {
    const bool parton =
        pdg == PDG::PDG_gluon || (std::abs(pdg) >= 1 && std::abs(pdg) <= 6) || std::abs(pdg) == PDG::PDG_hard_jet;
    return parton && (branch.p.pdg > 0 || pdg <= 0 || pdg == PDG::PDG_gluon);
  });
}

// Match one user branch against one MG5 node without changing particle data
bool MG5NodeMatches(const AmplitudeTopologyNode &node, const MDecayBranch &branch) {
  if (!MG5ParticleMatches(node, branch) || node.daughters.size() != branch.legs.size()) { return false; }
  for (const auto &i : indices(node.daughters)) {
    if (!MG5NodeMatches(node.daughters[i], branch.legs[i])) { return false; }
  }
  return true;
}

// Match one user tree already in MG5 external-particle order
bool MG5TreeMatches(const AmplitudeTopology &topology, const std::vector<MDecayBranch> &tree) {
  if (topology.size() != tree.size()) { return false; }
  for (const auto &i : indices(topology)) {
    if (!MG5NodeMatches(topology[i], tree[i])) { return false; }
  }
  return true;
}

// Compute one recursive bijection in MG5 order without changing particle data
std::optional<std::vector<MDecayBranch>> MG5OrderedTree(const AmplitudeTopology         &topology,
                                                        const std::vector<MDecayBranch> &tree) {
  std::function<bool(const AmplitudeTopologyNode &, const MDecayBranch &, MDecayBranch &)> order_node;
  std::function<bool(const AmplitudeTopology &, const std::vector<MDecayBranch> &, std::vector<MDecayBranch> &)>
      order_tree;

  order_node = [&](const AmplitudeTopologyNode &node, const MDecayBranch &branch, MDecayBranch &ordered) {
    if (!MG5ParticleMatches(node, branch)) { return false; }
    std::vector<MDecayBranch> daughters;
    if (!order_tree(node.daughters, branch.legs, daughters)) { return false; }
    ordered      = branch;
    ordered.legs = std::move(daughters);
    return true;
  };

  order_tree = [&](const AmplitudeTopology &nodes, const std::vector<MDecayBranch> &branches,
                   std::vector<MDecayBranch> &ordered) {
    if (nodes.size() != branches.size()) { return false; }
    ordered.resize(nodes.size());
    std::vector<bool>                      used(branches.size(), false);
    const std::function<bool(std::size_t)> match_permutation = [&](const std::size_t node_index) {
      if (node_index == nodes.size()) { return true; }
      for (const auto &branch_index : indices(branches)) {
        MDecayBranch candidate;
        if (used[branch_index] || !order_node(nodes[node_index], branches[branch_index], candidate)) { continue; }
        used[branch_index]  = true;
        ordered[node_index] = std::move(candidate);
        if (match_permutation(node_index + 1)) { return true; }
        used[branch_index] = false;
      }
      return false;
    };
    return match_permutation(0);
  };

  std::vector<MDecayBranch> ordered;
  if (!order_tree(topology, tree, ordered)) { return std::nullopt; }
  return ordered;
}

// Match one user decay tree against every exact channel of an MG5 process
bool MG5ProcessSyntaxMatches(const Process &process, const std::vector<MDecayBranch> &tree) {
  if (!IsGenerated(process.matrix_element_form)) { return false; }
  if (process.channel_topologies.empty()) { return MG5TreeMatches(process.topology, tree); }
  return std::any_of(process.channel_topologies.begin(), process.channel_topologies.end(),
                     [&](const AmplitudeTopology &topology) { return MG5TreeMatches(topology, tree); });
}

// Compute every canonical tree reachable through an exact MG5 channel
std::vector<MG5Order> MG5ProcessOrderings(const Process &process, const std::vector<MDecayBranch> &tree) {
  std::vector<MG5Order> orderings;
  const auto            append = [&](const AmplitudeTopology &topology) {
    auto candidate = MG5OrderedTree(topology, tree);
    if (candidate.has_value()) {
      MG5Order order;
      order.topology = AmplitudeTopologyFromDecayTree(*candidate);
      order.tree     = std::move(*candidate);
      orderings.push_back(std::move(order));
    }
  };
  if (process.channel_topologies.empty()) {
    append(process.topology);
  } else {
    for (const auto &topology : process.channel_topologies) { append(topology); }
  }
  return orderings;
}

// Append one distinct canonical MG5 external-particle ordering
void AppendMG5Order(MG5Order order, std::vector<MG5Order> &orderings) {
  const bool duplicate = std::any_of(orderings.begin(), orderings.end(),
                                     [&](const auto &current) { return current.topology == order.topology; });
  if (!duplicate) { orderings.push_back(std::move(order)); }
}

// Append one unique user-facing MG5 final-state syntax
void AppendMG5Syntax(const Process &process, std::vector<std::string> &syntax) {
  if (std::find(syntax.begin(), syntax.end(), process.final_state_syntax) == syntax.end()) {
    syntax.push_back(process.final_state_syntax);
  }
}

// Join accepted MG5 final-state syntax in one concise diagnostic
std::string JoinMG5Syntax(const std::vector<std::string> &syntax) {
  std::string joined;
  for (const auto &i : indices(syntax)) {
    if (i != 0) { joined += ", "; }
    joined += "'" + syntax[i] + "'";
  }
  return joined;
}

// Find the narrowest registered process matching one event topology
std::optional<Process> MatchTopology(const std::vector<Process> &processes, const std::vector<MDecayBranch> &tree) {
  std::optional<Process> matched_process;
  for (const auto &process : processes) {
    if (!TopologyMatches(process.topology, tree) || !StablePDGsMatch(tree, process.stable_pdgs) ||
        !GeneratedTopologyMatches(process, tree)) {
      continue;
    }
    if (!matched_process.has_value() || ProcessIsStrictSubset(process, *matched_process)) {
      matched_process = process;
      continue;
    }
    if (!ProcessIsStrictSubset(*matched_process, process)) {
      throw std::logic_error(
          "Amplitude process definition resolved incomparable event "
          "topologies");
    }
  }
  return matched_process;
}

// Reject a decay structure incompatible with the requested root decay mode
void ValidateDecayStructure(const DecayStructure &decay_structure, const LORENTZSCALAR &lts) {
  const bool cascade = std::any_of(lts.decaytree.begin(), lts.decaytree.end(),
                                    [](const MDecayBranch &branch) { return !branch.legs.empty(); });
  if (lts.process.root_decay_mode == RootDecayMode::Isolated &&
      ((cascade && decay_structure.type == DecayType::Full) || !decay_structure.allows_isolated_resonance)) {
    throw std::invalid_argument("Amplitude process does not support isolated decay sampling");
  }
}

// Validate definitions and combine their advertised processes in order
std::vector<Process> CombinedProcesses(const std::vector<std::shared_ptr<const ProcessDefinition>> &definitions) {
  if (definitions.empty()) {
    throw std::invalid_argument("Composite amplitude process definition sequence must be nonempty");
  }

  std::vector<Process> processes;
  for (const auto &definition : definitions) {
    if (definition == nullptr) {
      throw std::invalid_argument("Composite amplitude process definition contains a null definition");
    }
    const auto supported = definition->Processes();
    processes.insert(processes.end(), supported.begin(), supported.end());
  }
  if (processes.empty()) {
    throw std::invalid_argument("Composite amplitude process definition process data must be nonempty");
  }
  return processes;
}

}  // namespace

// Order user syntax one-to-one with a generated MG5 external state
void OrderMG5ProcessSyntax(const std::vector<Process> &processes, std::vector<MDecayBranch> &decaytree,
                           const std::string &syntax) {
  std::vector<std::string> generated_syntax;
  std::vector<std::string> reordered_syntax;
  std::vector<MG5Order>    orderings;
  const bool               has_generic_parton = HasMG5PartonAlias(decaytree);
  bool                     all_generated      = !processes.empty();
  for (const auto &process : processes) {
    if (!IsGenerated(process.matrix_element_form)) {
      all_generated = false;
      continue;
    }
    AppendMG5Syntax(process, generated_syntax);
    if (MG5ProcessSyntaxMatches(process, decaytree)) { return; }
    if (has_generic_parton) { continue; }
    for (auto &order : MG5ProcessOrderings(process, decaytree)) {
      AppendMG5Syntax(process, reordered_syntax);
      AppendMG5Order(std::move(order), orderings);
    }
  }

  if (orderings.size() == 1) {
    auto       ordered  = std::move(orderings.front().tree);
    const bool accepted = std::any_of(processes.begin(), processes.end(),
                                      [&](const auto &process) { return MG5ProcessSyntaxMatches(process, ordered); });
    if (!accepted) { throw std::logic_error("MG5 final-state ordering did not produce an exact channel"); }
    decaytree = std::move(ordered);
    return;
  }
  if (orderings.size() > 1) {
    throw std::invalid_argument("Ambiguous MG5 final-state ordering in '" + syntax +
                                "': no unique one-to-one PDG ordering among " + JoinMG5Syntax(reordered_syntax));
  }
  if (all_generated) {
    throw std::invalid_argument("Unsupported MG5 final state '" + syntax + "': expected " +
                                JoinMG5Syntax(generated_syntax));
  }
}

// Resolve one exact topology with its event-dependent decay structure
std::optional<Process> ProcessDefinition::ResolveProcess(const LORENTZSCALAR &lts) const {
  auto process = MatchProcess(lts.decaytree);
  if (!process.has_value()) { return std::nullopt; }
  process->decay_structure = DecayStructureFor(lts);
  if (process->matrix_element_form == MatrixElementForm::Generated) {
    process->topology_signature = AmplitudeTopologySignature(lts.decaytree);
    process->topology           = AmplitudeTopologyFromDecayTree(lts.decaytree);
    process->channel_topologies = {process->topology};
    process->stable_pdgs        = StableFinalStatePDGs(lts.decaytree);
    process->topology_mode      = TopologyMode::Event;
  }
  return process;
}

// Use unconstrained parton selection unless an amplitude restricts it
std::optional<bool> ProcessDefinition::AcceptsPartonMode(const std::vector<MDecayBranch> &) const {
  return std::nullopt;
}

// Store one nonnull immutable process definition
ProcessFamily::ProcessFamily(std::shared_ptr<const ProcessDefinition> definition) : definition_(std::move(definition)) {
  if (definition_ == nullptr) { throw std::invalid_argument("Amplitude family requires a nonnull process definition"); }
}

// Compute the processes in this amplitude family
std::vector<Process> ProcessFamily::Processes() const { return definition_->Processes(); }

// Match one event topology in this amplitude family
std::optional<Process> ProcessFamily::MatchProcess(const std::vector<MDecayBranch> &decaytree) const {
  return definition_->MatchProcess(decaytree);
}

// Compute the optional generic-parton constraint of this amplitude family
std::optional<bool> ProcessFamily::AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const {
  return definition_->AcceptsPartonMode(decaytree);
}

// Compute the matrix-element decay structure for one event
DecayStructure ProcessFamily::DecayStructureFor(const LORENTZSCALAR &lts) const {
  return definition_->DecayStructureFor(lts);
}

// Store one nonempty generated process family
ProcessRegistry::ProcessRegistry(std::vector<Process> processes) : processes_(std::move(processes)) {
  if (processes_.empty()) { throw std::invalid_argument("Amplitude process registry must be nonempty"); }

  const std::string &process_family = processes_.front().process_family;
  if (process_family.empty()) { throw std::invalid_argument("Amplitude process family must be nonempty"); }
  for (const auto &process : processes_) {
    if (process.process_family != process_family) {
      throw std::invalid_argument("Registered amplitude processes must share one process family");
    }
    if (process.process_name.empty() || process.stable_final_state.empty() || process.final_state_syntax.empty() ||
        process.topology_signature.empty() || process.process_syntax.empty()) {
      throw std::invalid_argument("Registered amplitude process fields must be nonempty");
    }
    if (process.topology.empty() || process.stable_pdgs.empty()) {
      throw std::invalid_argument(
          "Registered amplitude processes require explicit topology and "
          "stable-PDG sets");
    }
    std::vector<int> expected_stable_pdgs;
    bool             has_particle_set = false;
    for (const auto &node : process.topology) { ValidateTopologyNode(node, expected_stable_pdgs, has_particle_set); }
    if (IsGenerated(process.matrix_element_form)) {
      if (process.decay_structure.type != DecayType::Full) {
        throw std::invalid_argument("Generated matrix elements must supply the complete stable-particle amplitude");
      }
    }
    if (process.stable_pdgs != expected_stable_pdgs) {
      throw std::invalid_argument("Registered amplitude stable PDGs disagree with its topology");
    }
    const TopologyMode expected_mode = has_particle_set ? TopologyMode::FlavourSet : TopologyMode::Exact;
    if (process.topology_mode != expected_mode) {
      throw std::invalid_argument(
          "Registered amplitude topology mode "
          "disagrees with its particle sets");
    }
    for (const auto &exact_index : indices(process.channel_topologies)) {
      const auto &exact = process.channel_topologies[exact_index];
      if (exact.empty()) { throw std::invalid_argument("Registered amplitude generated topology must be nonempty"); }
      std::vector<int> exact_stable_pdgs;
      bool             exact_has_particle_set = false;
      for (const auto &node : exact) { ValidateTopologyNode(node, exact_stable_pdgs, exact_has_particle_set); }
      if (exact_has_particle_set || TopologyHasHardJetAlias(exact)) {
        throw std::invalid_argument(
            "Registered amplitude generated topology must be exact and "
            "fully specified");
      }
      if (!TopologyIsSubset(exact, process.topology)) {
        throw std::invalid_argument(
            "Registered amplitude generated topology is outside its "
            "registered particle sets");
      }
      for (std::size_t previous = 0; previous < exact_index; ++previous) {
        if (process.channel_topologies[previous] == exact) {
          throw std::invalid_argument("Registered amplitude generated topology is duplicated");
        }
      }
    }
  }

  for (const auto &first : indices(processes_)) {
    for (std::size_t second = first + 1; second < processes_.size(); ++second) {
      const auto &first_process  = processes_[first];
      const auto &second_process = processes_[second];
      if (TopologiesOverlap(first_process.topology, second_process.topology) &&
          !ProcessIsStrictSubset(first_process, second_process) &&
          !ProcessIsStrictSubset(second_process, first_process)) {
        throw std::invalid_argument("Amplitude process registry has ambiguous topologies");
      }
    }
  }
}

// Compute a copy of the registered processes
std::vector<Process> ProcessRegistry::Processes() const { return processes_; }

// Match one ordered event decay topology against generated particle sets
std::optional<Process> ProcessRegistry::MatchProcess(const std::vector<MDecayBranch> &decaytree) const {
  if (decaytree.empty()) { return std::nullopt; }
  return MatchTopology(processes_, decaytree);
}

// Resolve the exact registered topology or reject the event
DecayStructure ProcessRegistry::DecayStructureFor(const LORENTZSCALAR &lts) const {
  const auto process = MatchProcess(lts.decaytree);
  if (!process.has_value()) {
    throw std::invalid_argument(
        "Amplitude process registry has no decay structure for "
        "the requested decay topology");
  }
  ValidateDecayStructure(process->decay_structure, lts);
  return process->decay_structure;
}

// Store the analytic amplitude identity and topology functions by value
AnalyticProcess::AnalyticProcess(std::string process_family, std::string process_name,
                                 const std::string &final_state_set, DecayStructure default_decay_structure,
                                 TopologyCondition topology_condition, DecayStructureFunction decay_structure_function,
                                 TopologyCondition parton_mode_condition)
    : process_set_{std::move(process_family),
                   std::move(process_name),
                   final_state_set,
                   final_state_set,
                   "particles(" + final_state_set + ")",
                   "analytic",
                   {},
                   {},
                   {},
                   default_decay_structure,
                   MatrixElementForm::Analytic,
                   TopologyMode::FlavourSet},
      topology_condition_(std::move(topology_condition)),
      decay_structure_function_(std::move(decay_structure_function)),
      parton_mode_condition_(std::move(parton_mode_condition)) {
  if (process_set_.process_family.empty() || process_set_.process_name.empty() || final_state_set.empty()) {
    throw std::invalid_argument(
        "Analytic amplitude process identity and "
        "final-state particle set must be nonempty");
  }
  if (!topology_condition_ || !decay_structure_function_) {
    throw std::invalid_argument(
        "Analytic amplitude process requires "
        "topology and decay-structure functions");
  }
}

// Compute the nonempty human-readable analytic particle set
std::vector<Process> AnalyticProcess::Processes() const { return {process_set_}; }

// Build one exact process from accepted event decay data
Process AnalyticProcess::ExactProcess(const std::vector<MDecayBranch> &decaytree,
                                      const DecayStructure            &decay_structure) const {
  Process process            = process_set_;
  process.topology_signature = AmplitudeTopologySignature(decaytree);
  process.topology           = AmplitudeTopologyFromDecayTree(decaytree);
  process.stable_pdgs        = StableFinalStatePDGs(decaytree);
  process.decay_structure    = decay_structure;
  process.topology_mode      = TopologyMode::Event;
  return process;
}

// Compute an exact process with the default analytic decay structure
std::optional<Process> AnalyticProcess::MatchProcess(const std::vector<MDecayBranch> &decaytree) const {
  if (decaytree.empty() || !topology_condition_(decaytree)) { return std::nullopt; }
  return ExactProcess(decaytree, process_set_.decay_structure);
}

// Compute an optional analytic constraint for a generic-parton mode
std::optional<bool> AnalyticProcess::AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const {
  if (!parton_mode_condition_) { return std::nullopt; }
  return parton_mode_condition_(decaytree);
}

// Compute an accepted event decay structure or reject the topology
DecayStructure AnalyticProcess::DecayStructureFor(const LORENTZSCALAR &lts) const {
  if (!topology_condition_(lts.decaytree)) {
    throw std::invalid_argument("Analytic amplitude has no decay structure for the requested topology");
  }
  const DecayStructure decay_structure = decay_structure_function_(lts);
  ValidateDecayStructure(decay_structure, lts);
  return decay_structure;
}

// Store a nonempty ordered sequence of process definitions
ProcessAlternatives::ProcessAlternatives(std::vector<std::shared_ptr<const ProcessDefinition>> definitions)
    : definitions_(std::move(definitions)), processes_(CombinedProcesses(definitions_)) {}

// Compute the combined processes in matching priority order
std::vector<Process> ProcessAlternatives::Processes() const { return processes_; }

// Compute the first exact event match in physical priority order
std::optional<Process> ProcessAlternatives::MatchProcess(const std::vector<MDecayBranch> &decaytree) const {
  for (const auto &definition : definitions_) {
    auto process = definition->MatchProcess(decaytree);
    if (process.has_value()) { return process; }
  }
  return std::nullopt;
}

// Combine optional generic-parton constraints in physical priority order
std::optional<bool> ProcessAlternatives::AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const {
  bool constrained = false;
  for (const auto &definition : definitions_) {
    const auto accepted = definition->AcceptsPartonMode(decaytree);
    if (!accepted.has_value()) { continue; }
    constrained = true;
    if (*accepted) { return true; }
  }
  return constrained ? std::optional<bool>(false) : std::nullopt;
}

// Compute the first matching decay structure or reject the event
DecayStructure ProcessAlternatives::DecayStructureFor(const LORENTZSCALAR &lts) const {
  for (const auto &definition : definitions_) {
    if (definition->MatchProcess(lts.decaytree).has_value()) { return definition->DecayStructureFor(lts); }
  }
  throw std::invalid_argument(
      "Amplitude process alternatives have no decay structure for "
      "the requested decay topology");
}

}  // namespace amplitude
}  // namespace gra
