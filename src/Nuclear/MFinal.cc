// Nuclear reaction completion and conservation checks on HepMC3 graphs
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MFinal.h"

#include <dlfcn.h>

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <set>
#include <stdexcept>
#include <utility>
#include <vector>

#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenVertex.h"

using gra::aux::indices;

namespace gra::nuclear {
namespace {

// Compute the local floating-point tolerance in a Cartesian four-vector
double Roundoff(const HepMC3::FourVector& p) {
  return 512.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(p.e()) + p.p3mod());
}

// Check Cartesian conservation without subtracting boosted invariant masses
bool Close(const HepMC3::FourVector& a, const HepMC3::FourVector& b) {
  const double tolerance = Roundoff(a) + Roundoff(b);
  return std::abs(a.e() - b.e()) <= tolerance && std::abs(a.px() - b.px()) <= tolerance &&
         std::abs(a.py() - b.py()) <= tolerance && std::abs(a.pz() - b.pz()) <= tolerance;
}

// Compute baryon number and electric charge from nuclear IDs and the PDG table
std::array<int, 2> Quantum(const int id, const MPDG& pdg) {
  const int sign = id < 0 ? -1 : 1;
  if (IsNuclearPDG(id)) {
    const auto nucleus = DecodeNuclearPDG(id);
    return {sign * static_cast<int>(nucleus.a), sign * 3 * static_cast<int>(nucleus.z)};
  }
  const auto found = pdg.PDG_table.find(id);
  if (found == pdg.PDG_table.cend() || std::abs(id) < 11 || std::abs(id) == 21 ||
      (std::abs(id) >= 81 && std::abs(id) <= 100)) {
    throw PhaseSpaceFailure("ValidateNuclearDecay: unphysical or unknown particle PDG " + std::to_string(id));
  }
  const auto& particle = found->second;
  // The three valence-quark digits distinguish ordinary baryons from mesons
  const int code   = std::abs(id);
  const int baryon = code >= 1000 && code % 10000 >= 1000 && particle.spinX2 % 2 == 1 ? sign : 0;
  return {baryon, particle.chargeX3};
}

// Check a physical particle against its declared generated mass
void ValidateParticle(const HepMC3::ConstGenParticlePtr& particle) {
  const auto&  p    = particle->momentum();
  const double mass = particle->generated_mass();
  if (!gra::AllFinite(std::array<double, 5>{p.px(), p.py(), p.pz(), p.e(), mass}) || !(p.e() > 0.0) || mass < 0.0 ||
      std::abs(p.m2() - mass * mass) > Roundoff(p) * std::max(1.0, p.e() + p.p3mod())) {
    throw PhaseSpaceFailure("ValidateNuclearDecay: invalid physical particle momentum or mass");
  }
  const bool terminal = particle->end_vertex() == nullptr;
  if ((terminal && particle->status() != 1) || (!terminal && particle->status() != 2)) {
    throw PhaseSpaceFailure("ValidateNuclearDecay: require status 1 final particles and status 2 decaying particles");
  }
}

// Check one local nuclear reaction vertex and its causal displacement
void ValidateVertex(const HepMC3::ConstGenVertexPtr& vertex, const MPDG& pdg) {
  HepMC3::FourVector incoming, outgoing;
  std::array<int, 2> in{}, out{};
  const auto&        position = vertex->position();
  if (vertex->particles_in().empty() || vertex->particles_out().empty() ||
      !gra::AllFinite(std::array<double, 4>{position.x(), position.y(), position.z(), position.t()})) {
    throw PhaseSpaceFailure("ValidateNuclearDecay: empty vertex or non-finite position");
  }
  for (const auto& particle : vertex->particles_in()) {
    incoming           = incoming + particle->momentum();
    const auto quantum = Quantum(particle->pid(), pdg);
    for (const auto& i : indices(in)) { in[i] += quantum[i]; }
    const auto production = particle->production_vertex();
    if (production != nullptr) {
      const auto displacement = position - production->position();
      if (displacement.t() + Roundoff(position) < displacement.p3mod()) {
        throw PhaseSpaceFailure("ValidateNuclearDecay: acausal decay displacement");
      }
    }
  }
  for (const auto& particle : vertex->particles_out()) {
    outgoing           = outgoing + particle->momentum();
    const auto quantum = Quantum(particle->pid(), pdg);
    for (const auto& i : indices(out)) { out[i] += quantum[i]; }
  }
  if (in != out || !Close(incoming, outgoing)) {
    throw PhaseSpaceFailure("ValidateNuclearDecay: four-momentum, charge or baryon number is not conserved");
  }
}

// Load the normal particle data once and share it across cloned workers
std::shared_ptr<const MPDG> ParticleData() {
  auto data = std::make_shared<MPDG>();
  data->ReadParticleData();
  return data;
}

// Compare graph connections by the stable HepMC3 particle IDs
template <typename Particles>
std::set<int> ParticleIDs(const Particles& particles) {
  std::set<int> ids;
  for (const auto& particle : particles) { ids.insert(particle->id()); }
  return ids;
}

// Check that external completion only adds particles below the forward roots
void ValidateHard(const HepMC3::GenEvent& before, const HepMC3::GenEvent& after, const std::set<int>& roots) {
  if (before.momentum_unit() != after.momentum_unit() || before.length_unit() != after.length_unit() ||
      before.particles().size() > after.particles().size() || before.vertices().size() > after.vertices().size() ||
      before.weights().size() != after.weights().size()) {
    throw PhaseSpaceFailure("MFinal: backend changed hard event units, particles or weights");
  }
  for (const auto& particle : before.particles()) {
    const auto& current = after.particles().at(static_cast<std::size_t>(particle->id() - 1));
    if (current->id() != particle->id() || !Close(current->momentum(), particle->momentum()) ||
        std::abs(current->generated_mass() - particle->generated_mass()) >
            512.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, particle->generated_mass()) ||
        (!roots.count(particle->id()) &&
         (current->pid() != particle->pid() || current->status() != particle->status()))) {
      throw PhaseSpaceFailure("MFinal: backend changed an existing hard particle");
    }
    if (!roots.count(particle->id())) {
      const auto previous = particle->end_vertex();
      const auto next     = current->end_vertex();
      if ((previous == nullptr) != (next == nullptr) || (previous != nullptr && previous->id() != next->id())) {
        throw PhaseSpaceFailure("MFinal: backend changed hard event ancestry");
      }
    }
  }
  for (const auto& vertex : before.vertices()) {
    const auto& current = after.vertices().at(static_cast<std::size_t>(-vertex->id() - 1));
    if (current->id() != vertex->id() || current->status() != vertex->status() ||
        !Close(current->position(), vertex->position()) ||
        ParticleIDs(current->particles_in()) != ParticleIDs(vertex->particles_in()) ||
        ParticleIDs(current->particles_out()) != ParticleIDs(vertex->particles_out())) {
      throw PhaseSpaceFailure("MFinal: backend changed a hard vertex");
    }
  }
  for (const auto& i : indices(before.weights())) {
    if (!std::isfinite(after.weights()[i]) ||
        std::abs(before.weights()[i] - after.weights()[i]) >
            512.0 * std::numeric_limits<double>::epsilon() * std::max(1.0, std::abs(before.weights()[i]))) {
      throw PhaseSpaceFailure("MFinal: backend changed a hard event weight");
    }
  }
}

}  // namespace

// Validate a closed acyclic forward branch and collect its particle IDs
std::set<int> ValidateNuclearDecay(const HepMC3::GenEvent& event, const HepMC3::ConstGenParticlePtr& parent,
                                   const MPDG& pdg) {
  std::set<int>                                                 particles, vertices, active;
  const std::function<void(const HepMC3::ConstGenParticlePtr&)> visit = [&](const auto& particle) {
    if (particle == nullptr || particle->parent_event() != &event) {
      throw PhaseSpaceFailure("ValidateNuclearDecay: foreign particle");
    }
    if (active.count(particle->id())) { throw PhaseSpaceFailure("ValidateNuclearDecay: cyclic decay graph"); }
    if (!particles.insert(particle->id()).second) { return; }
    active.insert(particle->id());
    ValidateParticle(particle);
    (void)Quantum(particle->pid(), pdg);
    const auto vertex = particle->end_vertex();
    if (vertex != nullptr) {
      vertices.insert(vertex->id());
      for (const auto& daughter : vertex->particles_out()) { visit(daughter); }
    }
    active.erase(particle->id());
  };
  visit(parent);
  for (const auto id : vertices) {
    const auto vertex = event.vertices().at(static_cast<std::size_t>(-id - 1));
    for (const auto& particle : vertex->particles_in()) {
      if (!particles.count(particle->id())) { throw PhaseSpaceFailure("ValidateNuclearDecay: external decay input"); }
    }
    ValidateVertex(vertex, pdg);
  }
  return particles;
}

// Find the ordered forward roots using the existing production leg labels
std::array<HepMC3::GenParticlePtr, 2> ForwardParents(HepMC3::GenEvent& event) {
  std::array<HepMC3::GenParticlePtr, 2> roots{};
  for (const auto& particle : event.particles()) {
    const auto leg    = particle->attribute<HepMC3::IntAttribute>("graniitti_upc_leg");
    const auto vertex = particle->production_vertex();
    if (leg == nullptr || vertex == nullptr || vertex->particles_in().size() != 1 ||
        vertex->particles_in().front()->status() != 4) {
      continue;
    }
    if (leg->value() < 1 || leg->value() > 2 || roots[leg->value() - 1] != nullptr) {
      throw PhaseSpaceFailure("ForwardParents: invalid or repeated beam leg");
    }
    roots[leg->value() - 1] = particle;
  }
  if (roots[0] == nullptr || roots[1] == nullptr) {
    throw PhaseSpaceFailure("ForwardParents: both forward systems must be identified");
  }
  return roots;
}

// Own a reaction model with the normal particle data
MFinal::MFinal(std::unique_ptr<MReaction> model, std::array<NeutronSelection, 2> neutron)
    : model_(std::move(model)), pdg_(ParticleData()), neutron_(neutron) {
  if (model_ == nullptr || (model_->Tool().name.empty() || model_->Tool().version.empty())) {
    throw std::invalid_argument("MFinal: missing reaction model");
  }
}

// Load the optional model without introducing an external particle format
MFinal::MFinal(const std::string& library, const nlohmann::json& config, const UPCParam& upc)
    : pdg_(ParticleData()), neutron_(upc.neutron) {
  void* handle = dlopen(library.c_str(), RTLD_NOW | RTLD_LOCAL);
  if (handle == nullptr) { throw std::invalid_argument("MFinal: cannot load " + library + ": " + dlerror()); }
  library_     = std::shared_ptr<void>(handle, [](void* value) { dlclose(value); });
  auto factory = reinterpret_cast<ReactionFactory>(dlsym(handle, "graniitti_nuclear_v0"));
  if (factory == nullptr) { throw std::invalid_argument("MFinal: missing graniitti_nuclear_v0 factory"); }
  model_.reset(factory(config.dump().c_str(), &upc));
  if (model_ == nullptr || (model_->Tool().name.empty() || model_->Tool().version.empty())) {
    throw std::invalid_argument("MFinal: factory returned no model");
  }
}

// Clone only the model inputs, leaving event-local proposals empty
MFinal::MFinal(const MFinal& other)
    : library_(other.library_), model_(other.model_->Clone()), pdg_(other.pdg_), neutron_(other.neutron_) {
  if (model_ == nullptr) { throw std::invalid_argument("MFinal: model clone is null"); }
}

// Copy through an independent worker to preserve external-library lifetimes
MFinal& MFinal::operator=(const MFinal& other) {
  if (this != &other) { *this = MFinal(other); }
  return *this;
}

// Release the previous model before unloading its external library
MFinal& MFinal::operator=(MFinal&& other) noexcept {
  if (this != &other) {
    model_.reset();
    library_  = std::move(other.library_);
    model_    = std::move(other.model_);
    pdg_      = std::move(other.pdg_);
    mass_     = other.mass_;
    neutron_  = other.neutron_;
    prepared_ = other.prepared_;
    other.Reset();
  }
  return *this;
}

// Sample a fresh mass proposal without changing the beam momenta
RecoilMass MFinal::SampleMasses(const HepMC3::GenEvent& beams, MRandom& random) {
  Reset();
  try {
    const auto beam_particles = beams.beams();
    if (beams.momentum_unit() != HepMC3::Units::GEV || beams.length_unit() != HepMC3::Units::MM ||
        beam_particles.size() != 2) {
      throw PhaseSpaceFailure("MFinal: expected two HepMC3 beams in GeV and mm");
    }
    const auto recoil = model_->SampleMasses(beams, random);
    mass_             = recoil.mass;
    for (const auto& i : indices(mass_)) {
      const auto& beam = beam_particles[i];
      if (!std::isfinite(mass_[i]) || !(mass_[i] > 0.0) || !std::isfinite(recoil.emd[i]) || recoil.emd[i] < 0.0 ||
          mass_[i] - recoil.emd[i] + Roundoff(beam->momentum()) < beam->generated_mass()) {
        throw PhaseSpaceFailure("MFinal: invalid forward invariant mass");
      }
      if (!IsNuclearPDG(beam->pid()) && std::abs(mass_[i] - beam->generated_mass()) > Roundoff(beam->momentum())) {
        throw PhaseSpaceFailure("MFinal: non-nuclear beams must remain elastic");
      }
    }
    prepared_ = true;
    return recoil;
  } catch (const std::exception& error) {
    throw PhaseSpaceFailure(std::string("MFinal::SampleMasses: ") + error.what());
  } catch (...) { throw PhaseSpaceFailure("MFinal::SampleMasses: external model failure"); }
}

// Complete the sampled reaction at the fixed hard kinematics
double MFinal::Complete(HepMC3::GenEvent& event, MRandom& random) {
  try {
    if (!prepared_) { throw PhaseSpaceFailure("MFinal: no mass proposal for this event"); }
    prepared_        = false;
    const auto roots = ForwardParents(event);
    std::set<int> nuclear_roots;
    for (const auto& i : indices(roots)) {
      if (!IsNuclearPDG(roots[i]->production_vertex()->particles_in().front()->pid())) { continue; }
      roots[i]->set_generated_mass(mass_[i]);
      nuclear_roots.insert(roots[i]->id());
    }
    const HepMC3::GenEvent before(event);
    const double           weight = model_->Complete(event, roots, random);
    if (!std::isfinite(weight) || weight < 0.0) { throw PhaseSpaceFailure("MFinal: invalid nuclear response weight"); }
    if (!(weight > 0.0)) { return 0.0; }
    ValidateHard(before, event, nuclear_roots);
    std::set<int> completed;
    bool          selected = true;
    for (const auto& leg : indices(roots)) {
      const auto& root = roots[leg];
      const auto  beam = root->production_vertex()->particles_in().front();
      if (!IsNuclearPDG(beam->pid())) { continue; }
      if (Quantum(root->pid(), *pdg_) != Quantum(beam->pid(), *pdg_)) {
        throw PhaseSpaceFailure("MFinal: forward system changed the beam charge or baryon number");
      }
      const auto branch  = ValidateNuclearDecay(event, root, *pdg_);
      std::size_t neutron = 0;
      for (const auto id : branch) {
        const auto& particle = event.particles()[id - 1];
        neutron += particle->status() == 1 && std::abs(particle->pid()) == PDG::PDG_n;
        if (!completed.insert(id).second) { throw PhaseSpaceFailure("MFinal: overlapping forward systems"); }
      }
      selected &= neutron_[leg].Accept(neutron);
    }
    for (const auto& particle : event.particles()) {
      if (static_cast<std::size_t>(particle->id()) > before.particles().size() && !completed.count(particle->id())) {
        throw PhaseSpaceFailure("MFinal: particle outside the forward decay branches");
      }
    }
    for (const auto& vertex : event.vertices()) {
      if (static_cast<std::size_t>(-vertex->id()) <= before.vertices().size()) { continue; }
      if (vertex->particles_in().empty() || vertex->particles_out().empty()) {
        throw PhaseSpaceFailure("MFinal: disconnected reaction vertex");
      }
      for (const auto& particle : vertex->particles_in()) {
        if (!completed.count(particle->id())) {
          throw PhaseSpaceFailure("MFinal: vertex outside the nuclear branches");
        }
      }
    }
    return selected ? weight : 0.0;
  } catch (const std::exception& error) {
    throw PhaseSpaceFailure(std::string("MFinal::Complete: ") + error.what());
  } catch (...) { throw PhaseSpaceFailure("MFinal::Complete: external model failure"); }
}

// Clear the event-local proposal
void MFinal::Reset() noexcept {
  mass_     = {};
  prepared_ = false;
}

// Compute the external physics model description
HepMC3::GenRunInfo::ToolInfo MFinal::Tool() const { return model_->Tool(); }

}  // namespace gra::nuclear
