// Matching for generated MG5 hard scattering subprocesses
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_SUBPROCESSSUM_H
#define AMP_MG5_SUBPROCESSSUM_H

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Tech/MAux.h"

using gra::aux::indices;

namespace gra {
namespace mg5 {

// One generated MG5 subprocess channel with exact incoming and final topology
struct Channel {
  std::array<int, 2>                       initial = {0, 0};
  std::vector<int>                         final;
  AmplitudeTopology                        topology;
  std::vector<int>                         external_color_representations;
  std::vector<mg5helas::ExternalColorFlow> external_color_flows;
};

// Compute whether one exact channel has a colored stable final-state leg
inline bool ChannelHasFinalStateColor(const Channel &channel) {
  if (channel.external_color_representations.size() != channel.final.size() + 2) { return false; }
  return std::any_of(channel.external_color_representations.begin() + 2, channel.external_color_representations.end(),
                     [](int representation) { return representation != 1; });
}

// One generated MG5 subprocess instance with all channels it represents
template <class ProcessBase>
struct Subprocess {
  std::unique_ptr<ProcessBase> process;
  std::vector<Channel>         channels;
};

// Store exact ordered incoming pairs accepted by one process family
using InitialStates = std::vector<std::array<int, 2>>;

// Compute true for one exact ordered stable final PDG sequence
inline bool FinalStateMatches(const std::vector<int> &channel, const std::vector<int> &request) {
  return channel == request;
}

// Compute true for one exact ordered numeric production and decay topology
inline bool TopologyMatches(const AmplitudeTopology &channel, const AmplitudeTopology &request) {
  return channel == request;
}

// Compute true when one channel matches the complete final state and topology
inline bool ChannelMatches(const Channel &channel, const AmplitudeTopology &topology, const std::vector<int> &final) {
  return TopologyMatches(channel.topology, topology) && FinalStateMatches(channel.final, final);
}

// Compute true when every external momentum component is unchanged
inline bool PreparedMomentumMatches(const M4Vec &first, const M4Vec &second) {
  return math::IsExactEqual(first.E(), second.E()) && math::IsExactEqual(first.Px(), second.Px()) &&
         math::IsExactEqual(first.Py(), second.Py()) && math::IsExactEqual(first.Pz(), second.Pz());
}

// Compute true for incoming particles supported by the massless MG5 matrix
// elements
inline bool IsSupportedInitialPDG(int pdg) {
  const int abs_pdg = std::abs(pdg);
  return (abs_pdg >= 1 && abs_pdg <= 5) || pdg == PDG::PDG_gluon || pdg == PDG::PDG_gamma ||
         (abs_pdg >= 11 && abs_pdg <= 14) || abs_pdg == 16;
}

// Reject one external color slot that disagrees with its signed representation
inline void ValidateColorFlowLeg(int representation, const mg5helas::ColorFlowLeg &leg) {
  if (std::abs(representation) == 6) {
    if ((representation == 6 && leg.color > 0 && leg.anticolor < 0) ||
        (representation == -6 && leg.color < 0 && leg.anticolor > 0)) { return; }
    throw std::invalid_argument("mg5::SubprocessSum sextet leg has invalid color-flow tags");
  }
  if (leg.color < 0 || leg.anticolor < 0) {
    throw std::invalid_argument("mg5::SubprocessSum channel has a negative color-flow tag");
  }
  if (representation == 1) {
    if (leg.color != 0 || leg.anticolor != 0) {
      throw std::invalid_argument("mg5::SubprocessSum singlet leg has color-flow tags");
    }
    return;
  }
  if (representation == 3) {
    if (leg.color == 0 || leg.anticolor != 0) {
      throw std::invalid_argument("mg5::SubprocessSum triplet leg has invalid color-flow tags");
    }
    return;
  }
  if (representation == -3) {
    if (leg.color != 0 || leg.anticolor == 0) {
      throw std::invalid_argument("mg5::SubprocessSum antitriplet leg has invalid color-flow tags");
    }
    return;
  }
  if (representation == 8) {
    if (leg.color == 0 || leg.anticolor == 0 || leg.color == leg.anticolor) {
      throw std::invalid_argument("mg5::SubprocessSum octet leg has invalid color-flow tags");
    }
    return;
  }
  throw std::invalid_argument("mg5::SubprocessSum channel has an unsupported color representation");
}

// Validate one raw MG5 color-flow row after crossing the incoming color slots
inline void ValidateExternalColorFlow(const std::vector<int>            &representations,
                                      const mg5helas::ExternalColorFlow &flow) {
  if (flow.size() != representations.size()) {
    throw std::invalid_argument(
        "mg5::SubprocessSum external color-flow leg count disagrees with "
        "the generated representations");
  }

  std::map<int, std::array<std::size_t, 2>> tag_counts;
  for (const auto &leg_index : indices(flow)) {
    ValidateColorFlowLeg(representations[leg_index], flow[leg_index]);
    int color     = flow[leg_index].color;
    int anticolor = flow[leg_index].anticolor;
    if (leg_index < 2) { std::swap(color, anticolor); }
    if (color != 0) { ++tag_counts[std::abs(color)][color < 0 ? 1 : 0]; }
    if (anticolor != 0) { ++tag_counts[std::abs(anticolor)][anticolor > 0 ? 1 : 0]; }
  }
  for (const auto &entry : tag_counts) {
    if (entry.second[0] != 1 || entry.second[1] != 1) {
      throw std::invalid_argument("mg5::SubprocessSum external color-flow tags do not close");
    }
  }
}

// Validate the converter-owned color color structure for one exact generated
// channel
inline void ValidateChannelColorStructure(const Channel &channel) {
  if (channel.external_color_representations.empty() && channel.external_color_flows.empty()) { return; }
  if (channel.external_color_representations.size() != 2 + channel.final.size()) {
    throw std::invalid_argument(
        "mg5::SubprocessSum channel color representations have the wrong "
        "external-leg count");
  }
  if (channel.external_color_flows.empty()) {
    throw std::invalid_argument(
        "mg5::SubprocessSum channel has color "
        "representations but no flows");
  }
  for (const auto &flow : channel.external_color_flows) {
    ValidateExternalColorFlow(channel.external_color_representations, flow);
  }
}

// Pair raw JAMP amplitude rows with the exact selected generated channel color
// structure
inline bool PairHardColorFlows(const Channel &channel, std::vector<std::vector<std::complex<double>>> amplitudes,
                               std::vector<mg5helas::HardColorFlow> &flows) {
  if (channel.external_color_flows.empty()) {
    if (amplitudes.size() > 1) { return false; }
  } else if (channel.external_color_flows.size() != amplitudes.size()) {
    return false;
  }

  flows.clear();
  flows.reserve(amplitudes.size());
  for (const auto &flow : indices(amplitudes)) {
    mg5helas::ExternalColorFlow external;
    if (!channel.external_color_flows.empty()) { external = channel.external_color_flows[flow]; }
    flows.push_back({std::move(amplitudes[flow]), std::move(external)});
  }
  return true;
}

// Build the complete ordered light-parton pair domain
inline InitialStates MasslessQCDInitialStates() {
  const std::vector<int> partons = {-5, -4, -3, -2, -1, 1, 2, 3, 4, 5, PDG::PDG_gluon};
  InitialStates          domain;
  domain.reserve(partons.size() * partons.size());
  for (const int first : partons) {
    for (const int second : partons) { domain.push_back({first, second}); }
  }
  return domain;
}

// Build the exact two-photon incoming domain
inline InitialStates PhotonInitialStates() { return {{{PDG::PDG_gamma, PDG::PDG_gamma}}}; }

// Reject generated channels outside the process-family incoming states
template <class ProcessBase>
void ValidateSubprocessInitialStates(const std::vector<Subprocess<ProcessBase>> &subprocesses,
                                     const InitialStates &domain, const std::string &context) {
  if (domain.empty()) { throw std::invalid_argument(context + ": incoming domain is empty"); }
  for (const auto &initial : domain) {
    if (!IsSupportedInitialPDG(initial[0]) || !IsSupportedInitialPDG(initial[1])) {
      throw std::invalid_argument(context + ": incoming domain contains an unsupported particle");
    }
  }
  for (const auto &subprocess : subprocesses) {
    for (const auto &channel : subprocess.channels) {
      if (std::find(domain.begin(), domain.end(), channel.initial) == domain.end()) {
        throw std::invalid_argument(context + ": generated channel is outside its incoming domain");
      }
    }
  }
}

// Resolve the event-scale coupling used by generated subprocesses
inline bool ResolveEffectiveAlphaQCD(const LORENTZSCALAR &lts, double alphas, double &effective_alpha_s) {
  if (!std::isfinite(lts.alphaQCD) || lts.alphaQCD < 0.0 || !std::isfinite(alphas) || alphas < 0.0) {
    effective_alpha_s = 0.0;
    return false;
  }
  effective_alpha_s = lts.alphaQCD > 0.0 ? lts.alphaQCD : alphas;
  return true;
}

// Validate one exact generated topology node and collect its stable leaves
inline void ValidateTopologyNode(const AmplitudeTopologyNode &node, std::vector<int> &stable_pdgs) {
  if (node.allowed_pdgs.size() != 1 || node.allowed_pdgs.front() == 0 ||
      std::abs(node.allowed_pdgs.front()) == PDG::PDG_hard_jet) {
    throw std::invalid_argument("mg5::SubprocessSum channel topology is not exact");
  }
  if (node.daughters.empty()) {
    stable_pdgs.push_back(node.allowed_pdgs.front());
    return;
  }
  for (const auto &daughter : node.daughters) { ValidateTopologyNode(daughter, stable_pdgs); }
}

// Reject unresolved outgoing hard-jet aliases before stable final matching
inline bool HasExplicitStableFinalState(const std::vector<int> &pdgs) {
  for (const int pdg : pdgs) {
    if (std::abs(pdg) == PDG::PDG_hard_jet) { return false; }
  }
  return true;
}

// Contract generated components with EPA photon sources when requested
inline std::vector<std::complex<double>> PhotonAmplitudes(LORENTZSCALAR                                  &lts,
                                                          const std::vector<mg5helas::HelicityComponent> &components,
                                                          const M4Vec &p1, const M4Vec &p2, bool coherent_epa,
                                                          const mg5helas::EPAHardFrame *frame = nullptr) {
  if (coherent_epa) {
    if (frame != nullptr) { return mg5helas::ContractEPAHardPhotonSources(lts, components, *frame); }
    return mg5helas::ContractEPAPhotonSources(lts, components, p1, p2);
  }
  std::vector<std::complex<double>> amplitudes;
  amplitudes.reserve(components.size());
  for (const auto &component : components) { amplitudes.push_back(component.value); }
  return amplitudes;
}

// Event local subprocess sum for generated MG5 hard process families
template <class ProcessBase>
class SubprocessSum {
 public:
  // Construct an empty subprocess sum
  SubprocessSum() = default;

  // Construct a subprocess sum with generated subprocesses
  explicit SubprocessSum(std::vector<Subprocess<ProcessBase>> subprocesses_in) {
    SetSubprocesses(std::move(subprocesses_in));
  }

  // Replace generated subprocesses
  void SetSubprocesses(std::vector<Subprocess<ProcessBase>> subprocesses_in) {
    InvalidatePreparedState();
    ValidateSubprocesses(subprocesses_in);
    subprocesses = std::move(subprocesses_in);
  }

  // Compute model parameters shared by every subprocess in this family
  ParticleMap Particles() const { return subprocesses.at(0).process->Particles(); }

  // Compute the evaluated UFO electromagnetic coupling
  double AlphaQED() const { return subprocesses.at(0).process->AlphaQED(); }

  // Initialize every subprocess from the same model parameter card
  void InitParameters(SLHAReader card) {
    InvalidatePreparedState();
    for (auto &subprocess : subprocesses) {
      subprocess.process->InitParameters(card);
      ValidateIncomingMasses(*subprocess.process);
    }
  }

  // Discard all event-local prepared kinematics and evaluation outputs
  void InvalidatePreparedState() {
    momenta.clear();
    storage.clear();
    topology.clear();
    stable_final_pdgs.clear();
    input_initial = {};
    input_final.clear();
    projected_initial         = {};
    alpha_qcd                 = 0.0;
    alpha_qed_zero            = false;
    restore_final_symmetry    = false;
    prepared                  = false;
    contributing_subprocesses = 0;
    ClearHardColorFlows();
  }

  // Compute the unique generated channel matching one complete event topology
  const Channel *MatchingChannel(const LORENTZSCALAR &lts, const std::array<int, 2> &initial) const {
    const std::vector<const MDecayBranch *> leaves = mg5::StableDecayLeaves(lts.decaytree);
    std::vector<int>                        final;
    final.reserve(leaves.size());
    for (const MDecayBranch *leaf : leaves) { final.push_back(leaf->p.pdg); }
    if (!HasExplicitStableFinalState(final)) { return nullptr; }
    const AmplitudeTopology event_topology = amplitude::AmplitudeTopologyFromDecayTree(lts.decaytree);

    const Channel *selected = nullptr;
    for (const auto &subprocess : subprocesses) {
      for (const auto &channel : subprocess.channels) {
        if (channel.initial != initial || !ChannelMatches(channel, event_topology, final)) { continue; }
        if (selected != nullptr) { return nullptr; }
        selected = &channel;
      }
    }
    if (selected == nullptr) { return nullptr; }
    return selected;
  }

  // Compute the unique channel matching the event incoming and final states
  const Channel *MatchingChannel(const LORENTZSCALAR &lts) const { return MatchingChannel(lts, {lts.id1, lts.id2}); }

  // Compute true when prepared final-state quantum numbers match one event
  bool PreparedFinalStateMatches(const LORENTZSCALAR &lts) const {
    if (!prepared || topology != amplitude::AmplitudeTopologyFromDecayTree(lts.decaytree)) { return false; }
    const auto leaves = mg5::StableDecayLeaves(lts.decaytree);
    if (leaves.size() != stable_final_pdgs.size()) { return false; }
    for (const auto &i : indices(leaves)) {
      if (leaves[i]->p.pdg != stable_final_pdgs[i]) { return false; }
    }
    return true;
  }

  // Compute true when prepared external kinematics match one event
  bool PreparedKinematicsMatches(const LORENTZSCALAR &lts) const {
    if (!prepared || !PreparedMomentumMatches(input_initial[0], lts.q1) ||
        !PreparedMomentumMatches(input_initial[1], lts.q2) ||
        restore_final_symmetry != (lts.process.root_decay_mode != RootDecayMode::Isolated)) {
      return false;
    }
    const auto leaves = mg5::StableDecayLeaves(lts.decaytree);
    if (leaves.size() != input_final.size()) { return false; }
    for (const auto &i : indices(leaves)) {
      if (!PreparedMomentumMatches(input_final[i], leaves[i]->p4)) { return false; }
    }
    return true;
  }

  // Compute true when the complete prepared event key matches one event
  bool PreparedStateMatches(const LORENTZSCALAR &lts, double alphas) const {
    double effective_alpha_s = 0.0;
    return PreparedFinalStateMatches(lts) && PreparedKinematicsMatches(lts) &&
           alpha_qed_zero == mg5helas::AlphaQEDAtZero(lts) &&
           ResolveEffectiveAlphaQCD(lts, alphas, effective_alpha_s) && math::IsExactEqual(alpha_qcd, effective_alpha_s);
  }

  // Reject an unregistered decay topology and invalidate prepared event state
  bool RequireGeneratedTopology(const amplitude::ProcessDefinition &definition, LORENTZSCALAR &lts) {
    // Cache only topology membership, independent of momenta and screening
    if (validated_definition != &definition) {
      validated_topologies.clear();
      validated_definition = &definition;
    }
    const auto requested = amplitude::AmplitudeTopologyFromDecayTree(lts.decaytree);
    const bool known = std::find(validated_topologies.begin(), validated_topologies.end(), requested) !=
                       validated_topologies.end();
    if (known || definition.MatchProcess(lts.decaytree).has_value()) {
      if (!known) { validated_topologies.push_back(requested); }
      if (prepared && !PreparedFinalStateMatches(lts)) { InvalidatePreparedState(); }
      return true;
    }
    InvalidatePreparedState();
    lts.hamp.clear();
    return false;
  }

  // Prepare momenta and stable final filters once for one phase-space point
  mg5helas::EvaluationStatus Prepare(LORENTZSCALAR &lts, double alphas, bool use_epa_hard = false) {
    InvalidatePreparedState();
    lts.hamp.clear();
    lts.epa_hard.Clear();
    alpha_qed_zero = mg5helas::AlphaQEDAtZero(lts);
    for (auto &subprocess : subprocesses) { subprocess.process->setAlphaQEDZero(alpha_qed_zero); }
    const std::vector<const MDecayBranch *> leaves = mg5::StableDecayLeaves(lts.decaytree);
    topology                                       = amplitude::AmplitudeTopologyFromDecayTree(lts.decaytree);
    stable_final_pdgs.reserve(leaves.size());
    for (const MDecayBranch *leaf : leaves) { stable_final_pdgs.push_back(leaf->p.pdg); }
    if (!HasExplicitStableFinalState(stable_final_pdgs)) { return mg5helas::EvaluationStatus::AmplitudeFailure; }
    const std::vector<double> *model_masses = nullptr;
    for (const auto &subprocess : subprocesses) {
      for (const auto &channel : subprocess.channels) {
        if (ChannelMatches(channel, topology, stable_final_pdgs)) {
          model_masses = &subprocess.process->getMasses();
          break;
        }
      }
      if (model_masses != nullptr) { break; }
    }
    if (model_masses == nullptr) {
      InvalidatePreparedState();
      lts.hamp.clear();
      return mg5helas::EvaluationStatus::AmplitudeFailure;
    }
    restore_final_symmetry = lts.process.root_decay_mode != RootDecayMode::Isolated;
    input_initial          = {lts.q1, lts.q2};
    input_final            = mg5::StableLeafMomenta(leaves);

    if (!mg5::OnShellFinal(input_final, *model_masses)) {
      InvalidatePreparedState();
      return mg5helas::EvaluationStatus::KinematicsFailure;
    }

    std::vector<M4Vec> pf = input_final;
    M4Vec              p1 = lts.q1;
    M4Vec              p2 = lts.q2;
    if (use_epa_hard) {
      if (!mg5helas::PrepareEPAHardFrame(lts, pf, lts.epa_hard.frame)) {
        return mg5helas::EvaluationStatus::KinematicsFailure;
      }
      p1 = lts.epa_hard.frame.incoming[0];
      p2 = lts.epa_hard.frame.incoming[1];
    } else if (!mg5helas::PrepareOnShellKinematics(lts, pf, p1, p2)) {
      return mg5helas::EvaluationStatus::KinematicsFailure;
    }
    projected_initial = {p1, p2};
    mg5::BuildMG5Momenta(p1, p2, pf, storage, momenta);

    if (!ResolveEffectiveAlphaQCD(lts, alphas, alpha_qcd)) { return mg5helas::EvaluationStatus::AmplitudeFailure; }
    if (!mg5::UsesFullDecayChainMode(lts)) { return mg5helas::EvaluationStatus::AmplitudeFailure; }

    prepared                  = true;
    contributing_subprocesses = 0;
    return mg5helas::EvaluationStatus::Success;
  }

  // Calculate exact complex helicity components for the current incoming
  // flavours
  mg5helas::EvaluationStatus CalcPreparedHelicityAmp2(LORENTZSCALAR &lts, double &amp2) {
    lts.hamp.clear();
    amp2 = 0.0;
    // Color-flow amplitudes belong to this prepared phase-space evaluation only
    ClearHardColorFlows();
    if (!prepared) {
      contributing_subprocesses = 0;
      return mg5helas::EvaluationStatus::KinematicsFailure;
    }

    const std::array<int, 2> initial = {lts.id1, lts.id2};
    contributing_subprocesses        = 0;

    for (auto &subprocess : subprocesses) {
      const double final_state_weight = SubprocessWeight(subprocess, initial);
      if (!(final_state_weight > 0.0)) { continue; }
      const Channel *selected_channel = MatchingChannel(lts, initial);
      if (selected_channel == nullptr) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      ++contributing_subprocesses;
      subprocess.process->setMomenta(momenta);
      subprocess.process->setInitial(lts.id1, lts.id2);
      subprocess.process->setAlphaS(alpha_qcd);

      const double scale = std::sqrt(final_state_weight);
      // Accumulate each MG5 flow over the same helicity rows as the matrix
      // element
      std::vector<std::vector<std::complex<double>>> channel_flows;
      bool                                           flow_dimensions_set = false;
      const auto                                     components          = subprocess.process->helicityAmplitudes();
      if (components.empty()) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      for (const auto &component : components) {
        // Hard parton components are screened without an explicit HELAS source
        const std::complex<double> transport = mg5helas::IncomingPairTransport(
            momenta[0], lts.id1, component.incoming[0], momenta[1], lts.id2, component.incoming[1]);
        const std::complex<double> value = scale * transport * component.value;
        if (!std::isfinite(value.real()) || !std::isfinite(value.imag())) {
          lts.hamp.clear();
          amp2 = 0.0;
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        lts.hamp.push_back(value);
        amp2 += std::norm(value);
        if (!std::isfinite(amp2)) {
          lts.hamp.clear();
          amp2 = 0.0;
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        if (component.flow_values.empty()) {
          if (component.color == 0) {
            lts.hamp.clear();
            ClearHardColorFlows();
            return mg5helas::EvaluationStatus::AmplitudeFailure;
          }
          continue;
        }
        if (component.color != 0) {
          lts.hamp.clear();
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        if (!flow_dimensions_set) {
          channel_flows.resize(component.flow_values.size());
          flow_dimensions_set = true;
        } else if (channel_flows.size() != component.flow_values.size()) {
          lts.hamp.clear();
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        for (const auto &flow : indices(component.flow_values)) {
          const std::complex<double> flow_value = scale * transport * component.flow_values[flow];
          if (!std::isfinite(flow_value.real()) || !std::isfinite(flow_value.imag())) {
            lts.hamp.clear();
            amp2 = 0.0;
            ClearHardColorFlows();
            return mg5helas::EvaluationStatus::AmplitudeFailure;
          }
          channel_flows[flow].push_back(flow_value);
        }
      }
      if (!StoreHardColorFlows(*selected_channel, std::move(channel_flows))) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
    }
    return mg5helas::EvaluationStatus::Success;
  }

  // Calculate photon helicity components with optional coherent EPA sources
  mg5helas::EvaluationStatus CalcPreparedPhotonHelicityAmp2(LORENTZSCALAR &lts, bool coherent_epa, double &amp2) {
    lts.hamp.clear();
    ClearHardColorFlows();
    lts.epa_hard.amplitude.clear();
    lts.epa_hard.color.clear();
    lts.epa_hard.normalization = 1.0;
    lts.hard_color_flows.clear();
    amp2 = 0.0;
    if (!prepared) {
      contributing_subprocesses = 0;
      return mg5helas::EvaluationStatus::KinematicsFailure;
    }
    const std::array<int, 2> initial        = {PDG::PDG_gamma, PDG::PDG_gamma};
    const bool               store_epa_hard = coherent_epa && lts.epa_hard.frame.valid;
    const auto              *hard_frame     = store_epa_hard ? &lts.epa_hard.frame : nullptr;
    contributing_subprocesses               = 0;
    for (auto &subprocess : subprocesses) {
      const double final_state_weight = SubprocessWeight(subprocess, initial);
      if (!(final_state_weight > 0.0)) { continue; }
      const Channel *selected_channel = MatchingChannel(lts, initial);
      if (selected_channel == nullptr) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      ++contributing_subprocesses;
      subprocess.process->setMomenta(momenta);
      subprocess.process->setInitial(initial[0], initial[1]);
      subprocess.process->setAlphaS(alpha_qcd);

      const double                                          scale = std::sqrt(final_state_weight);
      std::vector<mg5helas::HelicityComponent>              exact_components;
      std::vector<std::vector<mg5helas::HelicityComponent>> flow_components;
      bool                                                  flow_dimensions_set = false;
      const auto                                            components = subprocess.process->helicityAmplitudes();
      if (components.empty()) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      for (auto component : components) {
        // Keep uncontracted photon helicities in one screening section
        const std::complex<double> transport =
            coherent_epa ? std::complex<double>(1.0, 0.0) : mg5helas::IncomingPairTransport(momenta[0], initial[0], component.incoming[0], momenta[1],
                                                           initial[1], component.incoming[1]);
        const auto raw_flows = std::move(component.flow_values);
        component.flow_values.clear();
        component.value *= scale * transport;
        if (!std::isfinite(component.value.real()) || !std::isfinite(component.value.imag())) {
          lts.hamp.clear();
          amp2 = 0.0;
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        exact_components.push_back(component);
        if (raw_flows.empty()) {
          if (component.color == 0) {
            lts.hamp.clear();
            ClearHardColorFlows();
            return mg5helas::EvaluationStatus::AmplitudeFailure;
          }
          continue;
        }
        if (component.color != 0) {
          lts.hamp.clear();
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        if (!flow_dimensions_set) {
          flow_components.resize(raw_flows.size());
          flow_dimensions_set = true;
        } else if (raw_flows.size() != flow_components.size()) {
          lts.hamp.clear();
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        for (const auto &flow : indices(raw_flows)) {
          auto flow_component  = component;
          flow_component.value = scale * transport * raw_flows[flow];
          if (!std::isfinite(flow_component.value.real()) || !std::isfinite(flow_component.value.imag())) {
            lts.hamp.clear();
            amp2 = 0.0;
            ClearHardColorFlows();
            return mg5helas::EvaluationStatus::AmplitudeFailure;
          }
          flow_components[flow].push_back(std::move(flow_component));
        }
      }

      const auto amplitudes =
          PhotonAmplitudes(lts, exact_components, projected_initial[0], projected_initial[1], coherent_epa, hard_frame);
      if (amplitudes.empty() || !gra::AllFinite(amplitudes)) {
        lts.hamp.clear();
        amp2 = 0.0;
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      amp2 += gra::SquaredNorm(amplitudes);
      if (!std::isfinite(amp2)) {
        lts.hamp.clear();
        amp2 = 0.0;
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      for (const auto &amplitude : amplitudes) {
        // MG5 contains the 1/4 photon spin average, while screening applies
        // that average after its coherent loop
        lts.hamp.push_back(2.0 * amplitude);
      }
      std::vector<std::vector<std::complex<double>>> channel_flows;
      channel_flows.reserve(flow_components.size());
      for (const auto &components : flow_components) {
        auto flow_amplitudes =
            PhotonAmplitudes(lts, components, projected_initial[0], projected_initial[1], coherent_epa, hard_frame);
        if (flow_amplitudes.empty() || !gra::AllFinite(flow_amplitudes)) {
          lts.hamp.clear();
          amp2 = 0.0;
          ClearHardColorFlows();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        channel_flows.push_back(std::move(flow_amplitudes));
      }
      if (!StoreHardColorFlows(*selected_channel, std::move(channel_flows))) {
        lts.hamp.clear();
        ClearHardColorFlows();
        return mg5helas::EvaluationStatus::AmplitudeFailure;
      }
      if (store_epa_hard) {
        lts.epa_hard.amplitude.insert(lts.epa_hard.amplitude.end(), exact_components.begin(), exact_components.end());
        if (flow_components.size() != selected_channel->external_color_flows.size()) {
          lts.epa_hard.Clear();
          return mg5helas::EvaluationStatus::AmplitudeFailure;
        }
        for (const auto &flow : indices(flow_components)) {
          lts.epa_hard.color.push_back(
              {std::move(flow_components[flow]), selected_channel->external_color_flows[flow]});
        }
      }
    }
    if (store_epa_hard && !lts.epa_hard.amplitude.empty()) { lts.epa_hard.normalization = 2.0; }
    lts.hard_color_flows = last_hard_color_flows;
    return mg5helas::EvaluationStatus::Success;
  }

  // Compute the number of subprocesses used in the previous call
  std::size_t ContributingSubprocesses() const { return contributing_subprocesses; }

  // Compute the number of generated subprocesses
  std::size_t SubprocessCount() const { return subprocesses.size(); }

  // Compute the maximum generated QCD amplitude order
  int AlphaSPower() const {
    int power = 0;
    for (const auto &subprocess : subprocesses) { power = std::max(power, subprocess.process->AlphaSPower()); }
    return power;
  }

  // Compute the maximum generated QED amplitude order
  int AlphaQEDPower() const {
    int power = 0;
    for (const auto &subprocess : subprocesses) { power = std::max(power, subprocess.process->AlphaQEDPower()); }
    return power;
  }

  // Compute raw JAMP rows paired with the exact generated external color flows
  const std::vector<mg5helas::HardColorFlow> &LastHardColorFlows() const { return last_hard_color_flows; }

  // Check whether event kinematics are prepared
  bool IsPrepared() const { return prepared; }

 private:
  // Clear the event-local generated color-flow result and numeric projection
  void ClearHardColorFlows() {
    has_prepared_color_channel = false;
    last_hard_color_flows.clear();
  }

  // Store one selected channel's raw JAMP rows and external flows atomically
  bool StoreHardColorFlows(const Channel &channel, std::vector<std::vector<std::complex<double>>> amplitudes) {
    if (has_prepared_color_channel) { return false; }
    std::vector<mg5helas::HardColorFlow> flows;
    if (!PairHardColorFlows(channel, std::move(amplitudes), flows)) { return false; }
    last_hard_color_flows      = std::move(flows);
    has_prepared_color_channel = true;
    return true;
  }

  // Require model masses compatible with the massless incoming momentum construction
  static void ValidateIncomingMasses(const ProcessBase &process) {
    const auto &masses = process.getMasses();
    if (masses.size() < 2 || std::fpclassify(masses[0]) != FP_ZERO || std::fpclassify(masses[1]) != FP_ZERO) {
      throw std::invalid_argument("MG5: incoming model masses must be zero for the massless incoming kinematics");
    }
  }

  // Validate converter-generated subprocesses before storing them
  static void ValidateSubprocesses(const std::vector<Subprocess<ProcessBase>> &subprocesses_in) {
    if (subprocesses_in.empty()) {
      throw std::invalid_argument("mg5::SubprocessSum requires at least one generated subprocess");
    }
    for (const auto &subprocess_index : indices(subprocesses_in)) {
      const auto &subprocess = subprocesses_in[subprocess_index];
      if (subprocess.process == nullptr) {
        throw std::invalid_argument("mg5::SubprocessSum received a null generated process");
      }
      ValidateIncomingMasses(*subprocess.process);
      if (subprocess.channels.empty()) {
        throw std::invalid_argument("mg5::SubprocessSum subprocess has no generated channels");
      }
      for (const auto &channel_index : indices(subprocess.channels)) {
        const auto &channel = subprocess.channels[channel_index];
        if (!IsSupportedInitialPDG(channel.initial[0]) || !IsSupportedInitialPDG(channel.initial[1])) {
          throw std::invalid_argument(
              "mg5::SubprocessSum channel has "
              "unresolved incoming particles");
        }
        if (channel.topology.empty() || channel.final.empty()) {
          throw std::invalid_argument("mg5::SubprocessSum channel has an empty final topology");
        }
        std::vector<int> topology_final;
        for (const auto &node : channel.topology) { ValidateTopologyNode(node, topology_final); }
        if (topology_final != channel.final) {
          throw std::invalid_argument(
              "mg5::SubprocessSum channel topology "
              "and stable final disagree");
        }
        ValidateChannelColorStructure(channel);
        for (std::size_t previous_subprocess = 0; previous_subprocess <= subprocess_index; ++previous_subprocess) {
          const std::size_t previous_limit = previous_subprocess == subprocess_index
                                                 ? channel_index
                                                 : subprocesses_in[previous_subprocess].channels.size();
          for (std::size_t previous_channel = 0; previous_channel < previous_limit; ++previous_channel) {
            const auto &other = subprocesses_in[previous_subprocess].channels[previous_channel];
            if (other.initial == channel.initial && ChannelMatches(other, channel.topology, channel.final)) {
              throw std::invalid_argument("mg5::SubprocessSum has duplicate generated channels");
            }
          }
        }
      }
    }
  }

  // Compute the selected fraction with resolved final-state symmetry restored
  double SubprocessWeight(const Subprocess<ProcessBase> &subprocess, const std::array<int, 2> &initial) const {
    std::size_t available       = 0;
    double      selected_weight = 0.0;
    for (const auto &channel : subprocess.channels) {
      if (channel.initial != initial) { continue; }
      ++available;
      if (!ChannelMatches(channel, topology, stable_final_pdgs)) { continue; }
      if (!restore_final_symmetry) {
        selected_weight += 1.0;
        continue;
      }
      selected_weight += mg5helas::FinalStateSymmetryFactor(channel.final);
    }
    if (available == 0 || !(selected_weight > 0.0)) { return 0.0; }
    return selected_weight / static_cast<double>(available);
  }

  std::vector<Subprocess<ProcessBase>> subprocesses;
  std::vector<std::array<double, 4>>   storage;
  std::vector<double *>                momenta;
  std::array<M4Vec, 2>                 projected_initial;
  AmplitudeTopology                    topology;
  const amplitude::ProcessDefinition *validated_definition = nullptr;
  std::vector<AmplitudeTopology>       validated_topologies;
  std::vector<int>                     stable_final_pdgs;
  std::array<M4Vec, 2>                 input_initial;
  std::vector<M4Vec>                   input_final;
  double                               alpha_qcd                  = 0.0;
  bool                                 alpha_qed_zero             = false;
  bool                                 restore_final_symmetry     = false;
  bool                                 prepared                   = false;
  std::size_t                          contributing_subprocesses  = 0;
  bool                                 has_prepared_color_channel = false;
  // Prepared raw JAMP rows paired with exact generated external flow color
  // structure
  std::vector<mg5helas::HardColorFlow> last_hard_color_flows;
  // Numeric projection kept for current hard-process amplitudes
};

}  // namespace mg5
}  // namespace gra

#endif
