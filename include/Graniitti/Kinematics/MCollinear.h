// Collinear 2->2 and 2->N phase space generator <P>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MCOLLINEAR_H
#define MCOLLINEAR_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Process/MProcess.h"

// HepMC33
#include "HepMC3/GenEvent.h"

namespace gra {
class MCollinear : public MProcess {
public:
  MCollinear();
  MCollinear(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});
  virtual ~MCollinear();

  void FinalizeProcessConfiguration() override;

  void PrintInit(bool silent) const;

protected:
  // Compute one parton process weight
  double ComputeEventWeight(const std::vector<double> &randvec,
                            MEventWeightState &aux) override;

  // Map one unit coordinate logarithmically onto the central mass squared
  static double SampleCentralMassSquared(double unit, double mass_min,
                                         double mass_max);

  // Compute the event local Jacobian of the logarithmic central mass map
  static double CentralMassSquaredJacobian(double mass_squared, double mass_min,
                                           double mass_max);

private:
  void Initialize();

protected:
  bool LoopKinematics(const std::array<double, 2> &p1p,
                      const std::array<double, 2> &p2p);

  // Build collinear 2->1 kinematics followed by the central decay
  bool B2BuildKin(double xbj1, double xbj2);

private:
  // 2->2 dim phase space
  bool B2RandomKin(const std::vector<double> &randvec);
  void B2RecordEvent(HepMC3::GenEvent &evt);

  double B2IntegralVolume() const;
  double B2PhaseSpaceWeight() const;

  void DecayWidthPS(double &exact) const;

  double collinear_integral_volume = 0.0;
};

} // namespace gra

#endif
