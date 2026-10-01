// Analytic quark box amplitudes for gg -> gamma gamma
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_GG_YY_H
#define AMP_GG_YY_H

#include <array>
#include <complex>
#include <cmath>
#include <string>
#include <vector>

#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_DurhamRegistry.h"
#include "Graniitti/MModelTune.h"

namespace gra {

class AMP_gg_yy final : public DurhamMG5Process {
 public:
  // Define the exact analytic two-photon final state
  static amplitude::Process Definition();

  // Read the quark masses and electromagnetic coupling before sampling
  explicit AMP_gg_yy(const SMParam &sm);

  // Compute the registered gluon fusion process name
  const std::string &Name() const override { return name_; }
  // Compute the two outgoing photon PDG codes
  const std::vector<int> &FinalPDGs() const override { return pdgs_; }
  // Compute the outgoing color singlet representations
  const std::vector<int> &FinalColorRepresentations() const override { return colors_; }
  // Compute the dimension of the delta_ab color basis
  std::size_t ColorCount() const override { return 1; }
  // Compute the rank after incoming gluon singlet projection
  std::size_t ColorRank() const override { return 1; }
  // Compute the number of physical helicity configurations
  std::size_t HelicityCount() const override { return 16; }
  // Compute the one quark box subprocess count
  std::size_t SubprocessCount() const override { return 1; }
  // Compute the normalized incoming singlet projection of delta_ab
  const std::vector<std::complex<double>> &ExactProjectors() const override { return projector_; }
  // Compute the empty shower color flow space for two photons
  const std::vector<std::vector<MColorFlow>> &FlowCandidates() const override { return flows_; }

  // Compute planar helicity amplitudes multiplying delta_ab
  std::array<std::complex<double>, 16> Helicity(double s, double t, double u, double alpha_s, double alpha) const;
  // Compute singlet projected amplitudes in the external HELAS polarization basis
  DurhamMG5Evaluation Evaluate(LORENTZSCALAR &lts, double alpha_s, M4Vec *hard_k1 = nullptr,
                              M4Vec *hard_k2 = nullptr) override;

 private:
  struct Quark {
    double mass;
    double charge2;
  };
  std::array<Quark, 6> quarks_;
  double alpha_ = 0.0;
  std::string name_ = "gg_yy";
  std::vector<int> pdgs_ = {22, 22};
  std::vector<int> colors_ = {1, 1};
  std::vector<std::complex<double>> projector_ = {std::sqrt(8.0)};
  std::vector<std::vector<MColorFlow>> flows_;
};

}  // namespace gra

#endif
