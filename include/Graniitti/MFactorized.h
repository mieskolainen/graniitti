// Factorized 2->3 phase space class
//
// (c) 2017-2019 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFACTORIZED_H
#define MFACTORIZED_H

// C++
#include <complex>
#include <random>
#include <vector>

// Own
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"

// HepMC33
#include "HepMC3/GenEvent.h"

namespace gra {
class MFactorized : public MProcess {
 public:
  MFactorized();
  MFactorized(std::string process, const std::vector<aux::OneCMD> &syntax);
  virtual ~MFactorized();

  void post_Constructor();

  double operator()(const std::vector<double> &randvec, AuxIntData &aux) {
    return EventWeight(randvec, aux);
  }
  double EventWeight(const std::vector<double> &randvec, AuxIntData &aux);
  bool   EventRecord(HepMC3::GenEvent &evt);
  void   PrintInit(bool silent) const;

 private:
  void Initialize();
  bool LoopKinematics(const std::vector<double> &p1p, const std::vector<double> &p2p);
  bool FiducialCuts() const;

  // 5+1-Dim phase space, 2->3
  bool B51RandomKin(const std::vector<double> &randvec);
  bool B51BuildKin(double pt1, double pt2, double phi1, double phi2, double yX, double m2X,
                   double m1, double m2);
  void B51RecordEvent(HepMC3::GenEvent &evt);

  double B51IntegralVolume() const;
  double B51PhaseSpaceWeight() const;

  void DecayWidthPS(double &exact) const;

  // Dynamic sampling boundaries based on resonance position and width
  double M_MIN = 0.0;
  double M_MAX = 0.0;
};

}  // namespace gra

#endif
