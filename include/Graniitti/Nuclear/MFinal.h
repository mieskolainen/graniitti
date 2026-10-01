// Nuclear reaction and de-excitation interfaces using HepMC3 records
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARFINAL_H
#define MNUCLEARFINAL_H

#include <array>
#include <memory>
#include <set>
#include <string>

#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Sampling/MRandom.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"
#include "HepMC3/GenRunInfo.h"
#include "json.hpp"

namespace gra::nuclear {

// Carry net recoil masses and independent EMD excitation in GeV
struct RecoilMass {
  std::array<double, 2> mass{};
  std::array<double, 2> emd{};
  std::shared_ptr<const ExcitationChannel> channel;
};

// Resolve the declared nuclear reaction with one worker-local model
class MReaction {
 public:
  // Release the external model in its own library
  virtual ~MReaction() = default;

  // Construct an independent worker with the same immutable physics inputs
  virtual std::unique_ptr<MReaction> Clone() const = 0;

  // Identify the physics model and its version
  virtual HepMC3::GenRunInfo::ToolInfo Tool() const = 0;

  // Propose both forward invariant masses in GeV before hard phase space
  // Keep the sampled nuclear state until Complete and replace it on every call
  virtual RecoilMass SampleMasses(const HepMC3::GenEvent& beams, MRandom& random) = 0;

  // Fill forward branches without changing the supplied hard particles
  // Compute the joint differential response divided by the mass/state proposal
  // Include all modeled internal phase space and any additional EMD exactly once
  // Match any mass-dependent production correction to the inclusive hard kernel
  virtual double Complete(HepMC3::GenEvent& event, const std::array<HepMC3::GenParticlePtr, 2>& forward,
                          MRandom& random) = 0;
};

// Find the two forward systems in a GRANIITTI central production record
std::array<HepMC3::GenParticlePtr, 2> ForwardParents(HepMC3::GenEvent& event);

// Check a complete nuclear branch, including momentum, charge and baryon number
std::set<int> ValidateNuclearDecay(const HepMC3::GenEvent& event, const HepMC3::ConstGenParticlePtr& parent, const MPDG& pdg);

// Own an optional external reaction model and its event-local proposal
class MFinal {
 public:
  // Own a directly supplied reaction model
  explicit MFinal(std::unique_ptr<MReaction> model, std::array<NeutronSelection, 2> neutron = {});

  // Load a reaction factory from a shared library built with this HepMC3 ABI
  MFinal(const std::string& library, const nlohmann::json& config, const UPCParam& upc);

  // Clone the physics model when a process worker is copied
  MFinal(const MFinal& other);

  // Replace a worker with an independent clone
  MFinal& operator=(const MFinal& other);

  // Move ownership without sharing mutable model state
  MFinal(MFinal&& other) noexcept = default;
  MFinal& operator=(MFinal&& other) noexcept;

  // Sample and check masses before constructing the hard momenta
  RecoilMass SampleMasses(const HepMC3::GenEvent& beams, MRandom& random);

  // Validate one backend completion before using its integration weight
  double Complete(HepMC3::GenEvent& event, MRandom& random);

  // Clear the proposal when a sampled point ends or fails
  void Reset() noexcept;

  // Compute the active model description
  HepMC3::GenRunInfo::ToolInfo Tool() const;

 private:
  std::shared_ptr<void>       library_;
  std::unique_ptr<MReaction>  model_;
  std::shared_ptr<const MPDG> pdg_;
  std::array<double, 2>       mass_{};
  std::array<NeutronSelection, 2>  neutron_{};
  bool                        prepared_ = false;
};

// External libraries export graniitti_nuclear_v0 with this factory signature
using ReactionFactory = MReaction* (*)(const char* config, const UPCParam* upc);

}  // namespace gra::nuclear

#endif
