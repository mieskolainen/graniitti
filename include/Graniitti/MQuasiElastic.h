// QuasiElastic (EL,SD,DD) and soft ND class
//
// (c) 2017-2019 Mikael Mieskolainen
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
#include "Graniitti/M4Vec.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MMatrix.h"
#include "Graniitti/MProcess.h"
#include "Graniitti/MSpin.h"

namespace gra {
// Matrix element dimension: " GeV^" << -(2*external_legs - 8)
class MQuasiElastic : public MProcess {
 public:
  MQuasiElastic();
  MQuasiElastic(std::string process, const std::vector<aux::OneCMD> &syntax);
  virtual ~MQuasiElastic();

  void post_Constructor();

  double operator()(const std::vector<double> &randvec, AuxIntData &aux) {
    return EventWeight(randvec, aux);
  }
  double EventWeight(const std::vector<double> &randvec, AuxIntData &aux);
  bool EventRecord(HepMC3::GenEvent &evt);
  void PrintInit(bool silent) const;

 private:
  bool LoopKinematics(const std::vector<double> &p1p, const std::vector<double> &p2p);
  bool FiducialCuts() const;
  void ConstructProcesses();

  // Non-Diffractive
  std::complex<double> PolySoft(const std::vector<double> &randvec);

  // 2/3-Dim phase space, 2->2 quasielastic
  bool B3RandomKin(const std::vector<double> &randvec);
  bool B3BuildKin(double s3, double s4, double t);
  bool   B3GetLorentzScalars();
  double B3IntegralVolume() const;
  double B3PhaseSpaceWeight() const;

  // Event by event integration boundaries
  double t_max     = 0.0;
  double t_min     = 0.0;
  double DD_M2_max = 0.0;

  // Multipomeron weight table
  std::vector<double> MAXPOMW;
};

}  // gra namespace ends

#endif
