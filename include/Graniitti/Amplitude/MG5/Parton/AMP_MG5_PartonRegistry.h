// Generated parton MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PARTONREGISTRY_H
#define AMP_MG5_PARTONREGISTRY_H

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Kinematics/MKinematics.h"

namespace gra {

// Own all values produced by one generated parton-process evaluation
struct PartonMG5Evaluation {
  mg5helas::EvaluationStatus status = mg5helas::EvaluationStatus::Success;
  double                     amp2   = 0.0;
  std::size_t                contributing_subprocesses = 0;
  std::vector<mg5helas::HardColorFlow>           color_flows;

  // Compute true when event kinematics and amplitudes were evaluated
  bool Valid() const { return mg5helas::EvaluationSucceeded(status); }
};

// Base class for one generated parton subprocess family
class PartonMG5Process : public MG5Process {
 public:
  // Store the exact process family represented by this matrix element
  explicit PartonMG5Process(std::vector<amplitude::Process> processes)
      : MG5Process(std::move(processes)) {}

  // Destroy one generated parton process family through its base class
  virtual ~PartonMG5Process() = default;

  // Prepare event momenta and final-state channels
  virtual mg5helas::EvaluationStatus Prepare(LORENTZSCALAR &lts, double alpha_s) = 0;

  // Evaluate one prepared event and return all event-local results by value
  virtual PartonMG5Evaluation EvaluatePrepared(LORENTZSCALAR &lts, double alpha_s) = 0;

  // Compute the number of generated subprocess channels
  virtual std::size_t SubprocessCount() const = 0;
};

// Compute generated parton families and their public process channels
std::vector<MG5ProcessInfo> PartonMG5ProcessInfos();

// Construct one generated parton process family
std::unique_ptr<PartonMG5Process> CreatePartonMG5Process(
    const std::string &process_family);

}  // namespace gra

#endif
