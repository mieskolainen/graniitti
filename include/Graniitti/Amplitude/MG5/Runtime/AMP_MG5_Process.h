// Common generated MadGraph process interface
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PROCESS_H
#define AMP_MG5_PROCESS_H

#include <cstddef>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"
#include "Graniitti/Particle/MParticle.h"

namespace gra {

// Manifest-defined public route for one generated matrix-element family
struct MG5ProcessInfo {
  std::string process_family;
  std::string channel;
  std::string parameter_card;
};

// Compute the evaluated model parameters for one registered generated process
mg5::ParticleMap MG5Particles(const amplitude::Process &process);

// Synchronize all generated particle masses and decay poles before sampling
void SynchronizeMG5DecayParameters(std::vector<MDecayBranch> &tree, const mg5::ParticleMap &particles);

// Common immutable process ownership shared by every generated runtime
class MG5Process : public amplitude::ProcessRegistry {
 public:
  // Store the exact processes represented by this generated runtime
  explicit MG5Process(std::vector<amplitude::Process> processes)
      : amplitude::ProcessRegistry(std::move(processes)) {}

  // Destroy one generated runtime through the common base class
  virtual ~MG5Process() = default;

  // Initialize independent parameters and all derived couplings before sampling
  virtual void InitParameters(SLHAReader) {
    throw std::logic_error("MG5: model parameter initialization is unavailable");
  }

  // Compute evaluated model pole parameters for input and phase space
  virtual mg5::ParticleMap Particles() const {
    throw std::logic_error("MG5: generated model particle parameters are unavailable");
  }

  // Compute the number of generated subprocess matrix elements
  virtual std::size_t SubprocessCount() const = 0;
};

}  // namespace gra

#endif
