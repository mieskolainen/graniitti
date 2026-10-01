// Generated Durham MadGraph process registry
//
// (c) 2017-2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_DURHAMREGISTRY_H
#define AMP_MG5_DURHAMREGISTRY_H

#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"

namespace gra {

// Amplitudes from one generated Durham hard-process evaluation
struct DurhamMG5Evaluation {
  mg5helas::EvaluationStatus status = mg5helas::EvaluationStatus::Success;
  std::vector<std::complex<double>> projected;
  std::vector<std::vector<std::complex<double>>> flow_projected;

  // Compute true when event kinematics and amplitudes were evaluated
  bool Valid() const { return mg5helas::EvaluationSucceeded(status); }
};

// Base class for one generated finite-Nc Durham process
class DurhamMG5Process : public MG5Process {
 public:
  // Store the exact process represented by this matrix element
  explicit DurhamMG5Process(
      std::vector<amplitude::Process> processes)
      : MG5Process(std::move(processes)) {}

  // Destroy one generated matrix element through its base class
  virtual ~DurhamMG5Process() = default;

  // Compute the stable MG2GRA process name
  virtual const std::string &Name() const = 0;
  // Compute final-state PDG codes in MadGraph momentum order
  virtual const std::vector<int> &FinalPDGs() const = 0;
  // Compute final-state SU(3) representations in MadGraph momentum order
  virtual const std::vector<int> &FinalColorRepresentations() const = 0;
  // Expose the event-dependent decay structure of the registered process
  using amplitude::ProcessRegistry::DecayStructureFor;
  // Compute the generated matrix-element decay structure
  DecayStructure DecayStructureFor() const {
    return Processes().at(0).decay_structure;
  }
  // Compute the number of MadGraph color-basis tensors
  virtual std::size_t ColorCount() const = 0;
  // Compute the rank of the incoming-singlet restricted color space
  virtual std::size_t ColorRank() const = 0;
  // Compute the complete MadGraph helicity count
  virtual std::size_t HelicityCount() const = 0;
  // Compute exact finite-Nc projector rows in MadGraph color-basis order
  virtual const std::vector<std::complex<double>> &ExactProjectors() const = 0;
  // Compute leading-color final-state candidates used only for shower tags
  virtual const std::vector<std::vector<MColorFlow>> &FlowCandidates() const = 0;
  // Evaluate exact and shower-partitioned projected helicity amplitudes
  virtual DurhamMG5Evaluation Evaluate(LORENTZSCALAR &lts, double alpha_s,
                                       M4Vec *hard_k1 = nullptr,
                                       M4Vec *hard_k2 = nullptr) = 0;
};

// Assign one generated shower-flow candidate to ordered stable decay leaves
bool AssignDurhamColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate);

// Compute true when an exact generated Durham matrix element matches this process
bool HasDurhamMG5Process(const amplitude::Process &process);

// Construct the generated Durham matrix element for this exact process
std::unique_ptr<DurhamMG5Process> CreateDurhamMG5Process(
    const amplitude::Process &process);

}  // namespace gra

#endif
