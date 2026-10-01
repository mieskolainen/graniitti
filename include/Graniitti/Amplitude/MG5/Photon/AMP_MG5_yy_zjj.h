// Generated MG5 amplitude for gamma-gamma Zjj production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_YY_ZJJ_H
#define AMP_MG5_YY_ZJJ_H

#include <complex>
#include <cstddef>
#include <vector>

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/ProcessBase.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {

// Sum every generated subprocess for gamma-gamma Zjj production
class AMP_MG5_yy_zjj : public PhotonMG5Process {
 public:
  // Construct all generated subprocesses
  AMP_MG5_yy_zjj();

  // Destroy all generated subprocesses
  ~AMP_MG5_yy_zjj() override = default;

  AMP_MG5_yy_zjj(const AMP_MG5_yy_zjj &) = delete;
  AMP_MG5_yy_zjj &operator=(const AMP_MG5_yy_zjj &) = delete;

  // Evaluate the summed matrix element squared
  mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) override;

  // Compute the Born matrix-element power of alpha_s
  int AlphaSPower() const noexcept override;

  // Compute the Born matrix-element power of alpha_QED
  int AlphaQEDPower() const noexcept override;

  // Sample one generated final-state color flow
  bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) override;

  // Compute whether the selected event channel has colored stable final states
  bool HasFinalStateColor() const override;

  // Compute the number of generated subprocesses
  std::size_t SubprocessCount() const override;

  // Initialize all model parameters before sampling
  void InitParameters(SLHAReader card) override { subprocess_sum.InitParameters(std::move(card)); }

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override { return subprocess_sum.Particles(); }

  // Compute the evaluated UFO electromagnetic coupling
  double DefaultAlphaQED() const override { return subprocess_sum.AlphaQED(); }

 private:
  mg5::SubprocessSum<MG5_YY_ZJJ::ProcessBase> subprocess_sum;
  bool has_final_state_color_ = false;
};

}  // namespace gra

#endif
