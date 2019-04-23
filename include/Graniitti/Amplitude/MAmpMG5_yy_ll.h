//==========================================================================
// This file has been automatically generated for C++ Standalone
// MadGraph5_aMC@NLO v. 3.0.0, 2018-05-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef MAMPMG5_yy_ll_H
#define MAMPMG5_yy_ll_H

#include "Graniitti/Amplitude/Parameters_sm.h"
#include "Graniitti/MForm.h"

class MAmpMG5_yy_ll {
 public:
  MAmpMG5_yy_ll();
  ~MAmpMG5_yy_ll();
  std::complex<double> CalcAmp(gra::LORENTZSCALAR &lts);

 private:
  // Constants for array limits
  static const int ninitial   = 2;
  static const int nexternal  = 4;
  static const int nprocesses = 1;

  // Private functions to calculate the matrix element for all
  // subprocesses
  // Calculate wavefunctions
  void calculate_wavefunctions(const int perm[], const int hel[]);
  static const int     nwavefuncs = 6;
  std::complex<double> w[nwavefuncs][18];
  static const int     namplitudes = 2;
  std::complex<double> amp[namplitudes];
  double               matrix_1_aa_epem();

  // Store the matrix element value from sigmaKin
  double matrix_element[nprocesses];

  // Color flows, used when selecting color
  double *jamp2[nprocesses];

  // Pointer to the model parameters
  Parameters_sm *pars;

  // vector with external particle masses
  vector<double> mME;

  // vector with momenta (to be changed each event)
  vector<double *> p;

  // Initial particle ids
  int id1, id2;
};

#endif