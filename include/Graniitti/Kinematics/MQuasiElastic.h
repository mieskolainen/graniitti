// Quasielastic EL, SD, DD and soft ND phase space generator <Q>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MQUASIELASTIC_H
#define MQUASIELASTIC_H

// C++
#include <complex>
#include <random>
#include <vector>

// HepMC33
#include "HepMC3/GenEvent.h"

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MQuasiElastic : public MProcess {
 public:
  MQuasiElastic();
  MQuasiElastic(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});
  virtual ~MQuasiElastic();

  void FinalizeProcessConfiguration() override;

  void PrintInit(bool silent) const;

 private:
  void Initialize();

 protected:
  // Compute one quasielastic process weight
  double ComputeEventWeight(const std::vector<double> &randvec, MEventWeightState &aux) override;

  // Construct one quasielastic event record
  bool BuildEventRecord(HepMC3::GenEvent &evt) override;

  bool LoopKinematics(const std::array<double, 2> &p1p, const std::array<double, 2> &p2p);

  // Refresh quasielastic invariants from the restored Born four-momenta
  bool RefreshBornKinematics() override;

  // Compute the eikonal kT2 table maximum required by elastic scattering
  double EikonalMaxKT2() const override;

  // Access the configured minimum sampled absolute momentum transfer
  double MinimumSampledAbsT() const;

  // Access the configured maximum sampled absolute momentum transfer
  double MaximumSampledAbsT() const;

  // Sample absolute momentum transfer from the elastic mixed proposal
  static double SampleElasticAbsT(double unit, double min_abs_t, double max_abs_t);

  // Compute the inverse density of the elastic mixed proposal
  static double ElasticAbsTJacobian(double abs_t, double min_abs_t, double max_abs_t);

  // 2/3-Dim phase space, 2->2 quasielastic
  bool B3BuildKin(double s3, double s4, double t);

  // Apply quasielastic forward fiducial cuts
  bool FiducialCuts() const;

 private:
  // Construct the sampled soft cut Pomeron chain
  bool BuildSoftChain();

  // Write the cut Pomeron chain and its two final remnant decays
  bool BuildSoftEventRecord(HepMC3::GenEvent &evt);

  // Fragment one soft system attached to its existing production particle
  bool WriteSoftDecay(const HepMC3::GenParticlePtr &particle, const M4Vec &momentum,
                      int baryon, int charge, HepMC3::GenEvent &evt);

  // 2/3-Dim phase space, 2->2 quasielastic
  bool   B3RandomKin(const std::vector<double> &randvec);
  bool   B3GetLorentzScalars();
  double B3IntegralVolume() const;
  double B3PhaseSpaceWeight() const;

  // Event by event integration boundaries
  double t_max               = 0.0;
  double t_min               = 0.0;
  double DD_M2_1_max         = 0.0;
  double DD_M2_max           = 0.0;
  double log_DD_M2_1_max     = 0.0;
  double log_DD_M2_max       = 0.0;
  double t_sampling_jacobian = 0.0;

  // Multipomeron weight table
  std::vector<double> MAXPOMW;
};

}  // namespace gra

#endif
