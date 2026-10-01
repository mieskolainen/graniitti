// Factorized 2->3 phase space generator <F>
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFACTORIZED_H
#define MFACTORIZED_H

// C++
#include <complex>
#include <random>
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
class MFactorized : public MProcess {
public:
  // Construct the default factorized F phase space
  MFactorized();

  // Construct the factorized F phase space for one process
  MFactorized(std::string process, const std::vector<aux::OneCMD> &syntax, MModelTunePtr tune = {});

  // Construct factorized kinematics behind a compatible F or C tag
  MFactorized(std::string process, const std::vector<aux::OneCMD> &syntax,
              const std::string &mode, MModelTunePtr tune = {});

  // Destroy one factorized process
  virtual ~MFactorized();

  void FinalizeProcessConfiguration() override;

  void PrintInit(bool silent) const;

protected:
  // Compute one factorized process weight
  double ComputeEventWeight(const std::vector<double> &randvec,
                            MEventWeightState &aux) override;

  // Map one unit coordinate logarithmically onto the central mass squared
  static double SampleCentralMassSquared(double unit, double mass_min,
                                         double mass_max);

  // Compute the event local Jacobian of the logarithmic central mass map
  static double CentralMassSquaredJacobian(double mass_squared, double mass_min,
                                           double mass_max);

private:
  // Initialize factorized kinematics behind one compatible selector tag
  void Initialize(const std::string &mode = "F");

protected:
  bool LoopKinematics(const std::array<double, 2> &p1p,
                      const std::array<double, 2> &p2p);

  // Refresh factorized invariants from the restored Born four-momenta
  bool RefreshBornKinematics() override;

  // 5+1-Dim phase space, 2->3
  bool B51BuildKin(double pt1, double pt2, double phi1, double phi2, double yX,
                   double m2X, double m1, double m2);

private:
  // 5+1-Dim phase space, 2->3
  bool B51RandomKin(const std::vector<double> &randvec);
  void B51RecordEvent(HepMC3::GenEvent &evt);

  double B51IntegralVolume() const;
  double B51PhaseSpaceWeight() const;

  void DecayWidthPS(double &exact) const;

  // Dynamic sampling boundaries based on resonance position and width
  double M_MIN = 0.0;
  double M_MAX = 0.0;
};

} // namespace gra

#endif
