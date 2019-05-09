//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.0.0, 2018-05-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
// Modified from the MadGraph C++ Standalone output
//
// ** This class is UNDER CONSTRUCTION ** :
// colorflow and color stripped helicity amplitudes needs to be propagated out
//
// Mikael Mieskolainen, 2018

#include <cmath>
#include <complex>
#include <vector>

// Own
#include "Graniitti/Amplitude/HelAmps_sm_gg_qqbar.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_qqbar.h"
#include "Graniitti/Amplitude/Parameters_sm.h"

#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"

MAmpMG5_gg_qqbar::MAmpMG5_gg_qqbar() {
  std::string param_card_name = gra::aux::GetBasePath(2) + "/MG5cards/" + "gg_qqbar_param_card.dat";

  // Instantiate the model class and set parameters that stay fixed
  // during run
  pars = Parameters_sm();
  SLHAReader slha(param_card_name);
  pars.setIndependentParameters(slha);
  pars.setIndependentCouplings();
  // pars.printIndependentParameters();
  // pars.printIndependentCouplings();

  // Set external particle masses for this matrix element
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);
  mME.push_back(pars.ZERO);

  jamp2 = std::vector<std::vector<double>>(nprocesses, std::vector<double>(2));
}

MAmpMG5_gg_qqbar::~MAmpMG5_gg_qqbar() {}

// Get amplitude
std::complex<double> MAmpMG5_gg_qqbar::CalcAmp(gra::LORENTZSCALAR &lts, double alpS) {
  // *** Set the parameters which change event by event ***
  pars.setDependentParameters(alpS);  // alphaS
  pars.setDependentCouplings();
  // ******************************************************

  const double mgluon = 0;

  // *** Set masses for HELAS ***
  const std::vector<double> masses = {mgluon, mgluon, lts.decaytree[0].p4.M(),
                                      lts.decaytree[1].p4.M()};
  mME = masses;

  // *** Set particle 4-momentum: [E,px,py,pz] convention here! ***
  gra::M4Vec p1_ = lts.q1;
  gra::M4Vec p2_ = lts.q2;
  gra::M4Vec p3_ = lts.decaytree[0].p4;
  gra::M4Vec p4_ = lts.decaytree[1].p4;  
  
  // Do kinematic transform
  gra::kinematics::OffShell2OnShell(p1_,p2_,p3_,p4_);

  // Set components
  double p1[] = {p1_.E(), p1_.Px(), p1_.Py(), p1_.Pz()};
  double p2[] = {p2_.E(), p2_.Px(), p2_.Py(), p2_.Pz()};
  double p3[] = {p3_.E(), p3_.Px(), p3_.Py(), p3_.Pz()};
  double p4[] = {p4_.E(), p4_.Px(), p4_.Py(), p4_.Pz()};

  p.clear();
  p.push_back(&p1[0]);
  p.push_back(&p2[0]);
  p.push_back(&p3[0]);
  p.push_back(&p4[0]);
  // --------------------------------------------------------------------

  const int ncomb = 16;

  lts.hamp.resize(ncomb);  // init helicity amplitudes

  const static std::vector<int> nonzero = {0, 1, 2,  3,  4,  5,  6,  7,
                                           8, 9, 10, 11, 12, 13, 14, 15};  // All
  // const static std::vector<int> nonzero = {5,6,9,10}; // High energy
  // limit / Helicity conserving

  // Reset color flows
  for (int i = 0; i < 2; i++) jamp2[0][i] = 0.;

  // Local variables and constants
  static bool goodhel[ncomb] = {false};
  static int  ntry = 0, sum_hel = 0, ngood = 0;
  static int  igood[ncomb];
  static int  jhel;
  // std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1}, {-1, -1, -1, 1}, {-1, -1, 1, -1}, {-1, -1, 1, 1},
      {-1, 1, -1, -1},  {-1, 1, -1, 1},  {-1, 1, 1, -1},  {-1, 1, 1, 1},
      {1, -1, -1, -1},  {1, -1, -1, 1},  {1, -1, 1, -1},  {1, -1, 1, 1},
      {1, 1, -1, -1},   {1, 1, -1, 1},   {1, 1, 1, -1},   {1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {256};

  // Reset the matrix elements
  for (int i = 0; i < nprocesses; i++) { matrix_element[i] = 0.; }
  // Define permutation
  int perm[nexternal];
  for (int i = 0; i < nexternal; i++) { perm[i] = i; }

  if (sum_hel == 0 || ntry < 10) {
    // Calculate the matrix element for all helicities
    for (int ihel = 0; ihel < ncomb; ihel++) {
      if (goodhel[ihel] || ntry < 2) {
        calculate_wavefunctions(perm, helicities[ihel]);
        t[0]           = matrix_1_gg_uux();
        lts.hamp[ihel] = 0.0;  // ** SET HELICITY AMPLITUDE: not enough, need
        // color structure **

        double tsum = 0;
        for (int iproc = 0; iproc < nprocesses; iproc++) {
          matrix_element[iproc] += t[iproc];
          tsum += t[iproc];
        }
        // Store which helicities give non-zero result
        if (tsum != 0. && !goodhel[ihel]) {
          goodhel[ihel] = true;
          ngood++;
          igood[ngood] = ihel;
        }
      }
    }
    jhel    = 0;
    sum_hel = min(sum_hel, ngood);
  } else {
    // Only use the "good" helicities
    for (int j = 0; j < sum_hel; j++) {
      jhel++;
      if (jhel >= ngood) jhel = 0;
      double hwgt             = double(ngood) / double(sum_hel);
      int    ihel             = igood[jhel];
      calculate_wavefunctions(perm, helicities[ihel]);
      t[0]           = matrix_1_gg_uux();
      lts.hamp[ihel] = 0.0;  // ** SET HELICITY AMPLITUDE: not enough, need color structure **

      for (int iproc = 0; iproc < nprocesses; iproc++) { matrix_element[iproc] += t[iproc] * hwgt; }
    }
  }

  for (int i = 0; i < nprocesses; i++) matrix_element[i] /= denominators[i];

  return std::sqrt(matrix_element[0]);  // square root, we take square later
}

// ------------------------------------------------------------------------------------------------
// From MadGraph automatic output

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void MAmpMG5_gg_qqbar::calculate_wavefunctions(const int perm[], const int hel[]) {
  // Calculate wavefunctions for all processes
  // int i, j;

  // Calculate all wavefunctions
  MG5_sm_gg_qqbar::vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  MG5_sm_gg_qqbar::vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  MG5_sm_gg_qqbar::oxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  MG5_sm_gg_qqbar::ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  MG5_sm_gg_qqbar::VVV1P0_1(w[0], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[4]);
  MG5_sm_gg_qqbar::FFV1_1(w[2], w[0], pars.GC_11, pars.ZERO, pars.ZERO, w[5]);
  MG5_sm_gg_qqbar::FFV1_2(w[3], w[0], pars.GC_11, pars.ZERO, pars.ZERO, w[6]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  MG5_sm_gg_qqbar::FFV1_0(w[3], w[2], w[4], pars.GC_11, amp[0]);
  MG5_sm_gg_qqbar::FFV1_0(w[3], w[5], w[1], pars.GC_11, amp[1]);
  MG5_sm_gg_qqbar::FFV1_0(w[6], w[2], w[1], pars.GC_11, amp[2]);
}
double MAmpMG5_gg_qqbar::matrix_1_gg_uux() {
  int i, j;
  // Local variables
  // const int ngraphs = 3;
  const int            ncolor = 2;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor]      = {3, 3};
  static const double cf[ncolor][ncolor] = {{16, -2}, {-2, 16}};

  // Calculate color flows
  jamp[0] = +std::complex<double>(0, 1) * amp[0] - amp[1];
  jamp[1] = -std::complex<double>(0, 1) * amp[0] - amp[2];

  // Sum and square the color flows to get the matrix element
  // jamp^*_i CF^{ij} jamp_j
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++) ztemp = ztemp + cf[i][j] * jamp[j];
    matrix                             = matrix + real(ztemp * conj(jamp[i])) / denom[i];
  }

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++) jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

  return matrix;
}
