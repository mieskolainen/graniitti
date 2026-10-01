// Generated MG5 amplitude for parton W production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PP_W_H
#define AMP_MG5_PP_W_H

#include <cstddef>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_SubprocessSum.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/ProcessBase.h"
#include "Graniitti/Kinematics/MKinematics.h"

namespace gra {

// Sum every generated subprocess for parton W production
class AMP_MG5_pp_w : public PartonMG5Process {
 public:
  // Construct all generated subprocesses
  AMP_MG5_pp_w();

  // Destroy all generated subprocesses
  ~AMP_MG5_pp_w() override = default;

  AMP_MG5_pp_w(const AMP_MG5_pp_w &) = delete;
  AMP_MG5_pp_w &operator=(const AMP_MG5_pp_w &) = delete;

  // Prepare event-local momenta and final-state channels
  mg5helas::EvaluationStatus Prepare(LORENTZSCALAR &lts,
                                     double alpha_s) override;

  // Evaluate one prepared event and return all event-local results
  PartonMG5Evaluation EvaluatePrepared(LORENTZSCALAR &lts,
                                     double alpha_s) override;

  // Compute the number of generated subprocesses
  std::size_t SubprocessCount() const override;

  // Initialize all model parameters before sampling
  void InitParameters(SLHAReader card) override { subprocess_sum.InitParameters(std::move(card)); }

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override { return subprocess_sum.Particles(); }

 private:
  mg5::SubprocessSum<MG5_PP_W::ProcessBase> subprocess_sum;
};

}  // namespace gra

#endif
