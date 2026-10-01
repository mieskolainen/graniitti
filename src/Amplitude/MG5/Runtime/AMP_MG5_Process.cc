// Common generated MadGraph process operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"

#include <cmath>
#include <stdexcept>
#include <string>

#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_DurhamRegistry.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Particle/MPDG.h"

namespace gra {

// Compute the evaluated model parameters for one registered generated process
mg5::ParticleMap MG5Particles(const amplitude::Process &process) {
  std::unique_ptr<MG5Process> runtime;
  if (process.process_family == "DURHAM") {
    runtime = CreateDurhamMG5Process(process);
  } else if (process.process_family == "PHOTON") {
    runtime = CreatePhotonMG5Process(process);
  } else {
    runtime = CreatePhotonMG5Process(process.process_family);
    if (runtime == nullptr) { runtime = CreatePartonMG5Process(process.process_family); }
  }
  if (runtime == nullptr) {
    throw std::invalid_argument("MG5: no generated model for " + process.process_name);
  }
  return runtime->Particles();
}

// Synchronize all generated particle masses and decay poles before sampling
void SynchronizeMG5DecayParameters(std::vector<MDecayBranch> &tree, const mg5::ParticleMap &particles) {
  for (auto &branch : tree) {
    if (std::abs(branch.p.pdg) == PDG::PDG_hard_jet) { continue; }
    const auto found = particles.find(std::abs(branch.p.pdg));
    if (found == particles.end()) {
      throw std::invalid_argument("MG5: missing model parameters for PDG " + std::to_string(branch.p.pdg));
    }
    branch.p.mass = found->second.mass;
    if (!branch.legs.empty()) {
      branch.p.width = found->second.width;
      branch.p.tau = branch.p.width > 0.0 ? PDG::hbar / branch.p.width : 0.0;
      SynchronizeMG5DecayParameters(branch.legs, particles);
    }
  }
}

}  // namespace gra
