//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.6.7, 2019-10-16
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//==========================================================================

#ifndef MG5_Sigma_sm_aa_epem_H
#define MG5_Sigma_sm_aa_epem_H

#include <complex>
#include <vector>

#include "Graniitti/Amplitude/Parameters_sm.h"
#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"

using namespace std;

//==========================================================================
// A class for calculating the matrix elements for
// Process: a a > e+ e- WEIGHTED<=4 @1
//--------------------------------------------------------------------------

class AMP_MG5_yy_ll {
 public:
  // Constructor.
  AMP_MG5_yy_ll() { initProc(gra::aux::GetBasePath(2) + "/MG5cards/param_card_yy_ll.dat"); }

  // Initialize process.
  void initProc(string param_card_name);

  // Calculate flavour-independent parts of cross section.
  double CalcAmp2(gra::LORENTZSCALAR &lts, double alphas);

  // Evaluate sigmaHat(sHat).
  double sigmaHat();

  // Info on the subprocess.
  string name() const { return "a a > e+ e- (sm)"; }

  int code() const { return 1; }

  const vector<double> &getMasses() const { return mME; }

  // Get and set momenta for matrix element evaluation
  vector<double *> getMomenta() { return p; }
  void             setMomenta(vector<double *> &momenta) { p = momenta; }
  void             setInitial(int inid1, int inid2) {
    id1 = inid1;
    id2 = inid2;
  }

  // Get matrix element vector
  const double *getMatrixElements() const { return matrix_element; }

  // Constants for array limits
  static const int ninitial   = 2;
  static const int nexternal  = 4;
  static const int nprocesses = 1;

 private:
  // Private functions to calculate the matrix element for all subprocesses
  // Calculate wavefunctions
  void                 calculate_wavefunctions(const int perm[], const int hel[]);
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
  Parameters_sm pars;  // GRANIITTI

  // vector with external particle masses
  vector<double> mME;

  // vector with momenta (to be changed each event)
  vector<double *> p;
  // Initial particle ids
  int id1, id2;
};


#endif  // MG5_Sigma_sm_aa_epem_H
