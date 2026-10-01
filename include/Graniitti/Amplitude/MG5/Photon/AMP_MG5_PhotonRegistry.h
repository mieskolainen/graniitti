// Generated photon MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PHOTONREGISTRY_H
#define AMP_MG5_PHOTONREGISTRY_H

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {

namespace flux {
struct EPASectorWeights;
}

// Clear color tags from all stable photon-process final states
bool ClearHardColorFlow(LORENTZSCALAR &lts);

// Common runtime interface for every generated gamma-gamma amplitude
class PhotonMG5Process : public MG5Process {
 public:
  // Store the exact processes represented by this generated runtime
  PhotonMG5Process(std::vector<amplitude::Process> processes, double alpha_qed)
      : MG5Process(std::move(processes)), alpha_qed_(alpha_qed) {}

  // Destroy one generated photon runtime through its base class
  virtual ~PhotonMG5Process() = default;

  // Evaluate the generated photon helicity amplitudes
  virtual mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) = 0;

  // Compute the Born matrix-element power of alpha_s
  virtual int AlphaSPower() const noexcept = 0;

  // Compute the Born matrix-element power of alpha_QED
  virtual int AlphaQEDPower() const noexcept = 0;

  // Compute the reference electromagnetic coupling
  virtual double DefaultAlphaQED() const { return alpha_qed_; }

  // Sample one generated final-state color flow
  virtual bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) = 0;

  // Compute whether the selected event channel has colored stable final states
  virtual bool HasFinalStateColor() const = 0;

 private:
  double alpha_qed_ = 0.0;
};

// Assign one shower-flow candidate to stable final-state leaves
bool AssignHardColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate);

// Assign the unique color-singlet flow of one quark-antiquark pair
bool AssignPhotonSingletQuarkPairColorFlow(LORENTZSCALAR &lts);

// Sample one prepared generated photon color flow
bool SampleHardColorFlow(
    LORENTZSCALAR &lts, MRandom &random,
    const std::vector<mg5helas::HardColorFlow> &flows,
    std::size_t expected_flows);

// Compute true when a generated photon matrix element matches this exact process
bool HasPhotonMG5Process(const amplitude::Process &process);

// Construct the generated photon matrix element for this exact process
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const amplitude::Process &process);

// Compute generated photon families and their public process channels
std::vector<MG5ProcessInfo> PhotonMG5ProcessInfos();

// Construct one generated photon subprocess family
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const std::string &process_family);

}  // namespace gra

#endif
