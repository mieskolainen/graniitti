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
#include "Graniitti/Amplitude/HelAmps_sm_gg_gg.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_gg.h"
#include "Graniitti/Amplitude/Parameters_sm.h"

#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"

MAmpMG5_gg_gg::MAmpMG5_gg_gg() {
  std::string param_card_name = gra::aux::GetBasePath(2) + "/MG5cards/" + "gg_gg_param_card.dat";

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

  jamp2 = std::vector<std::vector<double>>(nprocesses, std::vector<double>(6));
}

MAmpMG5_gg_gg::~MAmpMG5_gg_gg() {}

// Get amplitude
std::complex<double> MAmpMG5_gg_gg::CalcAmp(gra::LORENTZSCALAR &lts, double alpS) {
  // *** Set the parameters which change event by event ***
  pars.setDependentParameters(alpS);  // alphaS
  pars.setDependentCouplings();
  // ******************************************************

  const double mgluon = 0;

  // *** Set masses for HELAS ***
  const std::vector<double> masses = {mgluon, mgluon, mgluon, mgluon};
  mME                              = masses;

  // *** Set particle 4-momentum: [E,px,py,pz] convention here! ***
  double p1[] = {lts.q1.E(), lts.q1.Px(), lts.q1.Py(), lts.q1.Pz()};
  double p2[] = {lts.q2.E(), lts.q2.Px(), lts.q2.Py(), lts.q2.Pz()};
  double p3[] = {lts.decaytree[0].p4.E(), lts.decaytree[0].p4.Px(), lts.decaytree[0].p4.Py(),
                 lts.decaytree[0].p4.Pz()};
  double p4[] = {lts.decaytree[1].p4.E(), lts.decaytree[1].p4.Px(), lts.decaytree[1].p4.Py(),
                 lts.decaytree[1].p4.Pz()};

  /*
  // TEST INPUT
     double p1[] = {7.500000e+02,  0.000000e+00,  0.000000e+00,
  7.500000e+02};
     double p2[] = {7.500000e+02,  0.000000e+00,  0.000000e+00,
  -7.500000e+02};
     double p3[] = {7.500000e+02,  1.663864e+02,  6.672462e+02,
  -2.993294e+02};
     double p4[] = {7.500000e+02, -1.663864e+02, -6.672462e+02,
  2.993294e+02};
  */

  p.clear();
  p.push_back(&p1[0]);
  p.push_back(&p2[0]);
  p.push_back(&p3[0]);
  p.push_back(&p4[0]);

  const int ncomb = 16;

  lts.hamp.resize(ncomb);  // init helicity amplitudes

  const static std::vector<int> nonzero = {0, 1, 2,  3,  4,  5,  6,  7,
                                           8, 9, 10, 11, 12, 13, 14, 15};  // All
  // const static std::vector<int> nonzero = {5,6,9,10}; // High energy
  // limit / Helicity conserving

  // Reset color flows
  for (int i = 0; i < 6; i++) jamp2[0][i] = 0.;

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
  const int denominators[nprocesses] = {512};

  ntry = ntry + 1;

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
        t[0]           = matrix_1_gg_gg();
        lts.hamp[ihel] = t[0] / denominators[0];  // ** SET HELICITY AMPLITUDE **

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
      t[0]           = matrix_1_gg_gg();
      lts.hamp[ihel] = 0.0;  // ** SET HELICITY AMPLITUDE: not enough, need color structure **

      for (int iproc = 0; iproc < nprocesses; iproc++) { matrix_element[iproc] += t[iproc] * hwgt; }
    }
  }

  for (int i = 0; i < nprocesses; i++) matrix_element[i] /= denominators[i];

  return std::sqrt(matrix_element[0]);  // square root, we take square later
}

// ------------------------------------------------------------------------------------------------
// From MadGraph automatic output

// Calculate feynman diagrams amp[i] for the spesific helicity permutation
void MAmpMG5_gg_gg::calculate_wavefunctions(const int perm[], const int hel[]) {
  // std::cout << "GC_10 = " << pars.GC_10 << " GC_12 = " << pars.GC_12
  // << std::endl;

  // Calculate all wavefunctions
  MG5_sm_gg_gg::vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  MG5_sm_gg_gg::vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  MG5_sm_gg_gg::vxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  MG5_sm_gg_gg::vxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  MG5_sm_gg_gg::VVV1P0_1(w[0], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[4]);
  MG5_sm_gg_gg::VVV1P0_1(w[0], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[5]);
  MG5_sm_gg_gg::VVV1P0_1(w[0], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[6]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  MG5_sm_gg_gg::VVVV1_0(w[0], w[1], w[2], w[3], pars.GC_12, amp[0]);
  MG5_sm_gg_gg::VVVV3_0(w[0], w[1], w[2], w[3], pars.GC_12, amp[1]);
  MG5_sm_gg_gg::VVVV4_0(w[0], w[1], w[2], w[3], pars.GC_12, amp[2]);
  MG5_sm_gg_gg::VVV1_0(w[2], w[3], w[4], pars.GC_10, amp[3]);
  MG5_sm_gg_gg::VVV1_0(w[1], w[3], w[5], pars.GC_10, amp[4]);
  MG5_sm_gg_gg::VVV1_0(w[1], w[2], w[6], pars.GC_10, amp[5]);
}

/*

// Pythia gluon color coding

// Six different [color,anticolor] flows
static int colors[6][8] =
                {{4, 1, 1, 2, 3, 2, 4, 3},
                 {3, 1, 1, 2, 3, 4, 4, 2},
                 {4, 1, 3, 2, 3, 1, 4, 2},
                 {2, 1, 4, 2, 3, 1, 4, 3},
                 {3, 1, 4, 2, 3, 2, 4, 1},
                 {2, 1, 3, 2, 3, 4, 4, 1}};
*/

// Returns matrix element squared |M|^2 for particular helicity
double MAmpMG5_gg_gg::matrix_1_gg_gg() {
  int i, j;
  // Local variables
  // const int ngraphs = 6;
  const int            ncolor = 6;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];

  // The color matrix;
  static const double denom[ncolor]      = {6, 6, 6, 6, 6, 6};
  static const double cf[ncolor][ncolor] = {{19, -2, -2, -2, -2, 4}, {-2, 19, -2, 4, -2, -2},
                                            {-2, -2, 19, -2, 4, -2}, {-2, 4, -2, 19, -2, -2},
                                            {-2, -2, 4, -2, 19, -2}, {4, -2, -2, -2, -2, 19}};

  // Calculate color flows
  jamp[0] = +2. * (+amp[2] - amp[0] - amp[3] + amp[5]);
  jamp[1] = +2. * (+amp[0] + amp[1] + amp[3] + amp[4]);
  jamp[2] = +2. * (-amp[2] - amp[1] - amp[4] - amp[5]);
  jamp[3] = +2. * (+amp[0] + amp[1] + amp[3] + amp[4]);
  jamp[4] = +2. * (-amp[2] - amp[1] - amp[4] - amp[5]);
  jamp[5] = +2. * (+amp[2] - amp[0] - amp[3] + amp[5]);

  // Sum and square the color flows to get the matrix element
  // jamp^*_i CF^{ij} jamp_j
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {
    ztemp = 0.;
    for (j = 0; j < ncolor; j++) ztemp = ztemp + cf[i][j] * jamp[j];
    matrix                             = matrix + real(conj(jamp[i]) * ztemp) / denom[i];
  }

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++) {
                  jamp2[0][i] = real(jamp[i] * conj(jamp[i]));
  }

  return matrix;
}
