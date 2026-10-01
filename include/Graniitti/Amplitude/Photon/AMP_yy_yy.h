// Standard Model gamma gamma to gamma gamma one loop amplitude
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_YY_YY_H
#define AMP_YY_YY_H

#include <array>
#include <complex>
#include <cstddef>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/MModelTune.h"

namespace gra::lbyl {

// Evaluate all helicity amplitudes for real photon scattering
bool SMHelicityAmplitudes(double s, double t, double u, double alpha, const SMParam &sm, std::array<std::complex<double>, 16> &amplitude);

}  // namespace gra::lbyl

namespace gra {

// Compute the exact analytic light-by-light process definition
std::vector<amplitude::Process> LightByLightProcesses();

// Evaluate the analytic charged Standard Model light-by-light amplitude
class AMP_yy_yy final : public PhotonMG5Process {
 public:
  // Construct the analytic amplitude with validated Standard Model inputs
  explicit AMP_yy_yy(const SMParam &sm);

  // Destroy one event-local analytic amplitude
  ~AMP_yy_yy() override = default;

  // Compute the number of subprocess matrix elements
  std::size_t SubprocessCount() const override;
  // Compute zero because the charged-particle loop has no QCD vertex
  int AlphaSPower() const noexcept override;

  // Compute the Born matrix-element power of alpha_QED
  int AlphaQEDPower() const noexcept override;
  // Evaluate the complete charged Standard Model loop amplitude
  mg5helas::MatrixElementEvaluation Evaluate(LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) override;
  // Clear color tags for the colorless final state
  bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) override;
  // Compute false because both final-state photons are color singlets
  bool HasFinalStateColor() const override;

 private:
  const SMParam sm_;
  std::vector<std::complex<double>>              amplitudes_;
  std::vector<M4Vec>                             final_buffer_;
};

}  // namespace gra

#endif
