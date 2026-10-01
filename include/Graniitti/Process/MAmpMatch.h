// Amplitude process matching shared by GRANIITTI and MG2GRA generated registries
//
// MG2GRA generated registries include this header and use its matching API
// Keep develop/MG2GRA/modules/process_registry.py synchronized with this interface
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef GRANIITTI_AMP_MATCH_H
#define GRANIITTI_AMP_MATCH_H

#include <functional>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Process/MProcessState.h"

namespace gra {

// Describe one ordered decay-tree node with an explicit particle set
struct AmplitudeTopologyNode {
  std::vector<int>                   allowed_pdgs;
  std::vector<AmplitudeTopologyNode> daughters;

  // Compare the complete particle set and daughter tree
  bool operator==(const AmplitudeTopologyNode &other) const {
    return allowed_pdgs == other.allowed_pdgs && daughters == other.daughters;
  }

  // Compute true when any topology-node field differs
  bool operator!=(const AmplitudeTopologyNode &other) const { return !(*this == other); }
};

// Store one ordered production and decay topology
using AmplitudeTopology = std::vector<AmplitudeTopologyNode>;

namespace amplitude {

// Compute a readable ordered topology for one event
std::string AmplitudeTopologySignature(const std::vector<MDecayBranch> &decaytree);

// Compute the exact ordered numeric topology of one event decay tree
AmplitudeTopology AmplitudeTopologyFromDecayTree(const std::vector<MDecayBranch> &decaytree);

// Compute stable final-state PDGs in event momentum order
std::vector<int> StableFinalStatePDGs(const std::vector<MDecayBranch> &decaytree);

// Identify whether a process uses a generated or analytic matrix element
enum class MatrixElementForm { Generated, Analytic };

// Identify whether a topology is exact, a particle set or event specific
enum class TopologyMode { Exact, FlavourSet, Event };

// Compute true when a process uses a generated matrix element
constexpr bool IsGenerated(MatrixElementForm matrix_element_form) {
  return matrix_element_form == MatrixElementForm::Generated;
}

// Describe one advertised or event-resolved amplitude process
struct Process {
  std::string       process_family;
  std::string       process_name;
  std::string       stable_final_state;
  std::string       final_state_syntax;
  std::string       topology_signature;
  std::string       process_syntax;
  AmplitudeTopology topology;
  // Exact generated channels represented by one particle-set topology
  std::vector<AmplitudeTopology> channel_topologies;
  std::vector<int>               stable_pdgs;
  DecayStructure                 decay_structure;
  MatrixElementForm              matrix_element_form = MatrixElementForm::Generated;
  TopologyMode                   topology_mode       = TopologyMode::Exact;

  // Compare the complete process identity, topology and decay structure
  bool operator==(const Process &other) const {
    return process_family == other.process_family && process_name == other.process_name &&
           stable_final_state == other.stable_final_state && final_state_syntax == other.final_state_syntax &&
           topology_signature == other.topology_signature && process_syntax == other.process_syntax &&
           topology == other.topology && channel_topologies == other.channel_topologies &&
           stable_pdgs == other.stable_pdgs && decay_structure == other.decay_structure &&
           matrix_element_form == other.matrix_element_form && topology_mode == other.topology_mode;
  }

  // Compute true when any complete process field differs
  bool operator!=(const Process &other) const { return !(*this == other); }
};

// Order user syntax one-to-one with a generated MG5 external state
void OrderMG5ProcessSyntax(const std::vector<Process> &processes, std::vector<MDecayBranch> &decaytree,
                           const std::string &syntax);

// Define the processes and decay structures of one amplitude
class ProcessDefinition {
 public:
  // Destroy one process definition through its physical base class
  virtual ~ProcessDefinition() = default;

  // Compute supported processes without constructing an event
  virtual std::vector<Process> Processes() const = 0;

  // Match one decay tree against exact particles or particle sets
  virtual std::optional<Process> MatchProcess(const std::vector<MDecayBranch> &decaytree) const = 0;

  // Compute an optional constraint for one concrete generic-parton mode
  virtual std::optional<bool> AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const;

  // Compute the matrix-element decay structure for one event
  virtual DecayStructure DecayStructureFor(const LORENTZSCALAR &lts) const = 0;

  // Resolve one exact event topology and its decay structure
  std::optional<Process> ResolveProcess(const LORENTZSCALAR &lts) const;
};

// Represent one physical family of amplitudes
class ProcessFamily : public ProcessDefinition {
 public:
  // Store one nonnull immutable process definition
  explicit ProcessFamily(std::shared_ptr<const ProcessDefinition> definition);

  // Compute the processes in this amplitude family
  std::vector<Process> Processes() const final;

  // Match one event topology in this amplitude family
  std::optional<Process> MatchProcess(const std::vector<MDecayBranch> &decaytree) const final;

  // Compute the optional generic-parton constraint of this amplitude family
  std::optional<bool> AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const final;

  // Compute the decay structure for one event
  DecayStructure DecayStructureFor(const LORENTZSCALAR &lts) const final;

 protected:
  // Access the immutable process definition for amplitude composition
  std::shared_ptr<const ProcessDefinition> Family() const { return definition_; }

 private:
  const std::shared_ptr<const ProcessDefinition> definition_;
};

// Store implemented processes and match their physical topologies
class ProcessRegistry : public ProcessDefinition {
 public:
  // Store one nonempty generated process family
  explicit ProcessRegistry(std::vector<Process> processes);

  // Compute a copy of the registered processes
  std::vector<Process> Processes() const final;

  // Match one ordered topology and prefer the narrowest particle set
  std::optional<Process> MatchProcess(const std::vector<MDecayBranch> &decaytree) const final;

  // Resolve the exact registered topology or reject the event
  DecayStructure DecayStructureFor(const LORENTZSCALAR &lts) const final;

 private:
  const std::vector<Process> processes_;
};

// Select event decay trees supported by one analytic amplitude
using TopologyCondition = std::function<bool(const std::vector<MDecayBranch> &)>;

// Compute the decay structure of one analytic amplitude
using DecayStructureFunction = std::function<DecayStructure(const LORENTZSCALAR &)>;

// Resolve processes supported by one analytic amplitude
class AnalyticProcess final : public ProcessDefinition {
 public:
  // Store the analytic amplitude identity and topology functions by value
  AnalyticProcess(std::string process_family, std::string process_name, const std::string &final_state_set,
                  DecayStructure default_decay_structure, TopologyCondition topology_condition,
                  DecayStructureFunction decay_structure_function, TopologyCondition parton_mode_condition = {});

  // Compute the nonempty human-readable analytic particle set
  std::vector<Process> Processes() const final;

  // Compute an exact process with the default analytic decay structure
  std::optional<Process> MatchProcess(const std::vector<MDecayBranch> &decaytree) const final;

  // Compute an optional analytic constraint for a generic-parton mode
  std::optional<bool> AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const final;

  // Compute the accepted event decay structure or reject the topology
  DecayStructure DecayStructureFor(const LORENTZSCALAR &lts) const final;

 private:
  // Build one exact process from accepted event decay data
  Process ExactProcess(const std::vector<MDecayBranch> &decaytree, const DecayStructure &decay_structure) const;

  const Process                process_set_;
  const TopologyCondition      topology_condition_;
  const DecayStructureFunction decay_structure_function_;
  const TopologyCondition      parton_mode_condition_;
};

// Match immutable process definitions in physical priority order
class ProcessAlternatives final : public ProcessDefinition {
 public:
  // Store a nonempty ordered sequence of process definitions
  explicit ProcessAlternatives(std::vector<std::shared_ptr<const ProcessDefinition>> definitions);

  // Compute the combined processes in matching priority order
  std::vector<Process> Processes() const final;

  // Compute the first exact event match in physical priority order
  std::optional<Process> MatchProcess(const std::vector<MDecayBranch> &decaytree) const final;

  // Combine optional generic-parton constraints in physical priority order
  std::optional<bool> AcceptsPartonMode(const std::vector<MDecayBranch> &decaytree) const final;

  // Compute the first matching decay structure or reject the event
  DecayStructure DecayStructureFor(const LORENTZSCALAR &lts) const final;

 private:
  const std::vector<std::shared_ptr<const ProcessDefinition>> definitions_;
  const std::vector<Process>                                  processes_;
};

}  // namespace amplitude
}  // namespace gra

#endif
