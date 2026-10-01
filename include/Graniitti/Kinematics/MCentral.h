// Central multiparticle phase space generator <C>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MCENTRAL_H
#define MCENTRAL_H

// C++
#include <complex>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Tech/MAux.h"

// HepMC33
#include "HepMC3/GenEvent.h"

namespace gra {
class MCentral : public MProcess {
public:
  MCentral();
  MCentral(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});
  virtual ~MCentral();

  void FinalizeProcessConfiguration() override;

  void PrintInit(bool silent) const;

protected:
  // Compute one central phase space process weight
  double ComputeEventWeight(const std::vector<double> &randvec,
                            MEventWeightState &aux) override;

  // Compute the central mass ceiling allowed by one sampled forward point
  static double CentralMassKinematicMaximum(double sqrt_s, double forward_mass1,
                                            double forward_mass2,
                                            const M4Vec &forward1,
                                            const M4Vec &forward2);

  // Compute the squared Helmert ball radius for one mass and recoil ceiling
  static double
  MassConditionedTransverseRadius2(double mass_max, double recoil_pt2,
                                   unsigned int central_multiplicity);

  // Map a logarithmic Helmert hyperradius onto transverse differences
  static std::vector<M4Vec>
  MapTransverseLog(const std::vector<double> &radial_units,
                   const std::vector<double> &angle_units,
                   const M4Vec &central_transverse_momentum, double radius2,
                   double scale2, double &jacobian);

private:
  void Initialize();

protected:
  bool LoopKinematics(const std::array<double, 2> &p1p,
                      const std::array<double, 2> &p2p);

  // Refresh central phase space invariants from the restored Born four-momenta
  bool RefreshBornKinematics() override;

  // 3*N-4 dimensional phase space, 2->N
  bool BNBuildKin(unsigned int Nf, double pt1, double pt2, double phi1,
                  double phi2, const std::vector<double> &kt,
                  const std::vector<double> &phi, const std::vector<double> &y,
                  double m1, double m2);

private:
  // 3*N-4 dimensional phase space, 2->N
  bool BNRandomKin(unsigned int Nf, const std::vector<double> &randvec);

  double BNIntegralVolume() const;
  double BNPhaseSpaceWeight() const;

  // Auxiliary transverse difference vectors
  std::vector<M4Vec> pkt_;
};

} // namespace gra

#endif
