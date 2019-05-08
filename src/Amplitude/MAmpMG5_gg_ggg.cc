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
#include "Graniitti/Amplitude/HelAmps_sm_gg_ggg.h"
#include "Graniitti/Amplitude/MAmpMG5_gg_ggg.h"
#include "Graniitti/Amplitude/Parameters_sm.h"

#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"

MAmpMG5_gg_ggg::MAmpMG5_gg_ggg() {
  std::string param_card_name = gra::aux::GetBasePath(2) + "/MG5cards/" + "gg_ggg_param_card.dat";

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
  mME.push_back(pars.ZERO);

  jamp2 = std::vector<std::vector<double>>(nprocesses, std::vector<double>(24));
}

MAmpMG5_gg_ggg::~MAmpMG5_gg_ggg() {}

// Get amplitude
std::complex<double> MAmpMG5_gg_ggg::CalcAmp(gra::LORENTZSCALAR &lts, double alpS) {
  // *** Set the parameters which change event by event ***
  pars.setDependentParameters(alpS);  // alphaS
  pars.setDependentCouplings();
  // ******************************************************

  const double mgluon = 0;

  // *** Set masses for HELAS ***
  const std::vector<double> masses = {mgluon, mgluon, mgluon, mgluon, mgluon};
  mME                              = masses;

  // *** Set particle 4-momentum: [E,px,py,pz] convention here! ***
  double p1[] = {lts.q1.E(), lts.q1.Px(), lts.q1.Py(), lts.q1.Pz()};
  double p2[] = {lts.q2.E(), lts.q2.Px(), lts.q2.Py(), lts.q2.Pz()};
  double p3[] = {lts.decaytree[0].p4.E(), lts.decaytree[0].p4.Px(), lts.decaytree[0].p4.Py(),
                 lts.decaytree[0].p4.Pz()};
  double p4[] = {lts.decaytree[1].p4.E(), lts.decaytree[1].p4.Px(), lts.decaytree[1].p4.Py(),
                 lts.decaytree[1].p4.Pz()};
  double p5[] = {lts.decaytree[2].p4.E(), lts.decaytree[2].p4.Px(), lts.decaytree[2].p4.Py(),
                 lts.decaytree[2].p4.Pz()};

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
  p.push_back(&p5[0]);

  const int ncomb = 32;

  lts.hamp.resize(ncomb);  // init helicity amplitudes

  // const static std::vector<int> nonzero = {
  //    0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15}; // All
  // const static std::vector<int> nonzero = {5,6,9,10}; // High energy
  // limit / Helicity conserving

  // Reset color flows
  for (int i = 0; i < 24; i++) jamp2[0][i] = 0.;

  // Local variables and constants
  static bool goodhel[ncomb] = {false};
  static int  ntry = 0, sum_hel = 0, ngood = 0;
  static int  igood[ncomb];
  static int  jhel;
  // std::complex<double> * * wfs;
  double t[nprocesses];
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1, -1}, {-1, -1, -1, -1, 1}, {-1, -1, -1, 1, -1}, {-1, -1, -1, 1, 1},
      {-1, -1, 1, -1, -1},  {-1, -1, 1, -1, 1},  {-1, -1, 1, 1, -1},  {-1, -1, 1, 1, 1},
      {-1, 1, -1, -1, -1},  {-1, 1, -1, -1, 1},  {-1, 1, -1, 1, -1},  {-1, 1, -1, 1, 1},
      {-1, 1, 1, -1, -1},   {-1, 1, 1, -1, 1},   {-1, 1, 1, 1, -1},   {-1, 1, 1, 1, 1},
      {1, -1, -1, -1, -1},  {1, -1, -1, -1, 1},  {1, -1, -1, 1, -1},  {1, -1, -1, 1, 1},
      {1, -1, 1, -1, -1},   {1, -1, 1, -1, 1},   {1, -1, 1, 1, -1},   {1, -1, 1, 1, 1},
      {1, 1, -1, -1, -1},   {1, 1, -1, -1, 1},   {1, 1, -1, 1, -1},   {1, 1, -1, 1, 1},
      {1, 1, 1, -1, -1},    {1, 1, 1, -1, 1},    {1, 1, 1, 1, -1},    {1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {1536};

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
        t[0] = matrix_1_gg_ggg();

        lts.hamp[ihel] = gra::math::msqrt(t[0] / denominators[0]);  // ** SET HELICITY AMPLITUDE **

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
      t[0]           = matrix_1_gg_ggg();
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
void MAmpMG5_gg_ggg::calculate_wavefunctions(const int perm[], const int hel[]) {
  // Calculate wavefunctions for all processes
  // int i, j;

  // Calculate all wavefunctions
  MG5_sm_gg_ggg::vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  MG5_sm_gg_ggg::vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  MG5_sm_gg_ggg::vxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]);
  MG5_sm_gg_ggg::vxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]);
  MG5_sm_gg_ggg::vxxxxx(p[perm[4]], mME[4], hel[4], +1, w[4]);
  MG5_sm_gg_ggg::VVV1P0_1(w[0], w[1], pars.GC_10, pars.ZERO, pars.ZERO, w[5]);
  MG5_sm_gg_ggg::VVV1P0_1(w[2], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[6]);
  MG5_sm_gg_ggg::VVV1P0_1(w[2], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[7]);
  MG5_sm_gg_ggg::VVV1P0_1(w[3], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[8]);
  MG5_sm_gg_ggg::VVV1P0_1(w[0], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[9]);
  MG5_sm_gg_ggg::VVV1P0_1(w[1], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[10]);
  MG5_sm_gg_ggg::VVV1P0_1(w[1], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[11]);
  MG5_sm_gg_ggg::VVV1P0_1(w[0], w[3], pars.GC_10, pars.ZERO, pars.ZERO, w[12]);
  MG5_sm_gg_ggg::VVV1P0_1(w[1], w[2], pars.GC_10, pars.ZERO, pars.ZERO, w[13]);
  MG5_sm_gg_ggg::VVV1P0_1(w[0], w[4], pars.GC_10, pars.ZERO, pars.ZERO, w[14]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[15]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[16]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[1], w[2], pars.GC_12, pars.ZERO, pars.ZERO, w[17]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[18]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[19]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[1], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[20]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[21]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[22]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[1], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[23]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[24]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[25]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[2], w[3], pars.GC_12, pars.ZERO, pars.ZERO, w[26]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[27]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[28]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[2], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[29]);
  MG5_sm_gg_ggg::VVVV1P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[30]);
  MG5_sm_gg_ggg::VVVV3P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[31]);
  MG5_sm_gg_ggg::VVVV4P0_1(w[0], w[3], w[4], pars.GC_12, pars.ZERO, pars.ZERO, w[32]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  MG5_sm_gg_ggg::VVV1_0(w[5], w[6], w[4], pars.GC_10, amp[0]);
  MG5_sm_gg_ggg::VVV1_0(w[5], w[7], w[3], pars.GC_10, amp[1]);
  MG5_sm_gg_ggg::VVV1_0(w[5], w[2], w[8], pars.GC_10, amp[2]);
  MG5_sm_gg_ggg::VVVV1_0(w[2], w[3], w[4], w[5], pars.GC_12, amp[3]);
  MG5_sm_gg_ggg::VVVV3_0(w[2], w[3], w[4], w[5], pars.GC_12, amp[4]);
  MG5_sm_gg_ggg::VVVV4_0(w[2], w[3], w[4], w[5], pars.GC_12, amp[5]);
  MG5_sm_gg_ggg::VVV1_0(w[9], w[10], w[4], pars.GC_10, amp[6]);
  MG5_sm_gg_ggg::VVV1_0(w[9], w[11], w[3], pars.GC_10, amp[7]);
  MG5_sm_gg_ggg::VVV1_0(w[9], w[1], w[8], pars.GC_10, amp[8]);
  MG5_sm_gg_ggg::VVVV1_0(w[1], w[3], w[4], w[9], pars.GC_12, amp[9]);
  MG5_sm_gg_ggg::VVVV3_0(w[1], w[3], w[4], w[9], pars.GC_12, amp[10]);
  MG5_sm_gg_ggg::VVVV4_0(w[1], w[3], w[4], w[9], pars.GC_12, amp[11]);
  MG5_sm_gg_ggg::VVV1_0(w[12], w[13], w[4], pars.GC_10, amp[12]);
  MG5_sm_gg_ggg::VVV1_0(w[12], w[11], w[2], pars.GC_10, amp[13]);
  MG5_sm_gg_ggg::VVV1_0(w[12], w[1], w[7], pars.GC_10, amp[14]);
  MG5_sm_gg_ggg::VVVV1_0(w[1], w[2], w[4], w[12], pars.GC_12, amp[15]);
  MG5_sm_gg_ggg::VVVV3_0(w[1], w[2], w[4], w[12], pars.GC_12, amp[16]);
  MG5_sm_gg_ggg::VVVV4_0(w[1], w[2], w[4], w[12], pars.GC_12, amp[17]);
  MG5_sm_gg_ggg::VVV1_0(w[14], w[13], w[3], pars.GC_10, amp[18]);
  MG5_sm_gg_ggg::VVV1_0(w[14], w[10], w[2], pars.GC_10, amp[19]);
  MG5_sm_gg_ggg::VVV1_0(w[14], w[1], w[6], pars.GC_10, amp[20]);
  MG5_sm_gg_ggg::VVVV1_0(w[1], w[2], w[3], w[14], pars.GC_12, amp[21]);
  MG5_sm_gg_ggg::VVVV3_0(w[1], w[2], w[3], w[14], pars.GC_12, amp[22]);
  MG5_sm_gg_ggg::VVVV4_0(w[1], w[2], w[3], w[14], pars.GC_12, amp[23]);
  MG5_sm_gg_ggg::VVV1_0(w[0], w[13], w[8], pars.GC_10, amp[24]);
  MG5_sm_gg_ggg::VVV1_0(w[0], w[10], w[7], pars.GC_10, amp[25]);
  MG5_sm_gg_ggg::VVV1_0(w[0], w[11], w[6], pars.GC_10, amp[26]);
  MG5_sm_gg_ggg::VVV1_0(w[3], w[4], w[15], pars.GC_10, amp[27]);
  MG5_sm_gg_ggg::VVV1_0(w[3], w[4], w[16], pars.GC_10, amp[28]);
  MG5_sm_gg_ggg::VVV1_0(w[3], w[4], w[17], pars.GC_10, amp[29]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[4], w[18], pars.GC_10, amp[30]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[4], w[19], pars.GC_10, amp[31]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[4], w[20], pars.GC_10, amp[32]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[3], w[21], pars.GC_10, amp[33]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[3], w[22], pars.GC_10, amp[34]);
  MG5_sm_gg_ggg::VVV1_0(w[2], w[3], w[23], pars.GC_10, amp[35]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[4], w[24], pars.GC_10, amp[36]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[4], w[25], pars.GC_10, amp[37]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[4], w[26], pars.GC_10, amp[38]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[3], w[27], pars.GC_10, amp[39]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[3], w[28], pars.GC_10, amp[40]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[3], w[29], pars.GC_10, amp[41]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[2], w[30], pars.GC_10, amp[42]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[2], w[31], pars.GC_10, amp[43]);
  MG5_sm_gg_ggg::VVV1_0(w[1], w[2], w[32], pars.GC_10, amp[44]);
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
double MAmpMG5_gg_ggg::matrix_1_gg_ggg() {
  int i, j;
  // Local variables
  // const int ngraphs = 45;
  const int            ncolor = 24;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor] = {108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108,
                                       108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108, 108};
  static const double cf[ncolor][ncolor] = {{455, -58, -58, 14, 14, 68, -58, -4, 14, -58, 5,  -4,
                                             14,  5,   68,  -4, 14, 68, -58, -4, -4, 68,  68, -40},
                                            {-58, 455, 14, 68, -58, 14,  -4, -58, 5,  -4, 14, -58,
                                             -58, -4,  -4, 68, 68,  -40, 14, 5,   68, -4, 14, 68},
                                            {-58, 14, 455, -58, 68, 14, 14, 5,   68, -4,  14, 68,
                                             -58, -4, 14,  -58, 5,  -4, -4, -58, 68, -40, -4, 68},
                                            {14, 68,  -58, 455, 14, -58, -58, -4, -4, 68, 68, -40,
                                             -4, -58, 5,   -4,  14, -58, 5,   14, 14, 68, 68, -4},
                                            {14, -58, 68, 14,  455, -58, 5,   14, 14, 68,  68, -4,
                                             -4, -58, 68, -40, -4,  68,  -58, -4, 14, -58, 5,  -4},
                                            {68, 14, 14, -58, -58, 455, -4, -58, 68, -40, -4, 68,
                                             5,  14, 14, 68,  68,  -4,  -4, -58, 5,  -4,  14, -58},
                                            {-58, -4, 14, -58, 5,  -4, 455, -58, -58, 14, 14,  68,
                                             68,  -4, 14, 5,   68, 14, -4,  68,  -58, -4, -40, 68},
                                            {-4, -58, 5,   -4, 14,  -58, -58, 455, 14, 68, -58, 14,
                                             -4, 68,  -58, -4, -40, 68,  68,  -4,  14, 5,  68,  14},
                                            {14, 5,   68,  -4, 14, 68, -58, 14,  455, -58, 68, 14,
                                             14, -58, -58, -4, -4, 5,  68,  -40, -4,  -58, 68, -4},
                                            {-58, -4, -4, 68,  68,  -40, 14, 68, -58, 455, 14, -58,
                                             5,   -4, -4, -58, -58, 14,  14, 68, 5,   14,  -4, 68},
                                            {5,  14,  14, 68,  68, -4, 14, -58, 68,  14, 455, -58,
                                             68, -40, -4, -58, 68, -4, 14, -58, -58, -4, -4,  5},
                                            {-4, -58, 68, -40, -4, 68, 68, 14, 14, -58, -58, 455,
                                             14, 68,  5,  14,  -4, 68, 5,  -4, -4, -58, -58, 14},
                                            {14,  -58, -58, -4, -4, 5,  68, -4, 14,  5,  68,  14,
                                             455, -58, -58, 14, 14, 68, 68, -4, -40, 68, -58, -4},
                                            {5,   -4,  -4, -58, -58, 14, -4, 68, -58, -4, -40, 68,
                                             -58, 455, 14, 68,  -58, 14, -4, 68, 68,  14, 14,  5},
                                            {68,  -4, 14,  5,   68, 14, 14,  -58, -58, -4, -4, 5,
                                             -58, 14, 455, -58, 68, 14, -40, 68,  68,  -4, -4, -58},
                                            {-4, 68, -58, -4,  -40, 68,  5,  -4, -4, -58, -58, 14,
                                             14, 68, -58, 455, 14,  -58, 68, 14, -4, 68,  5,   14},
                                            {14, 68,  5,  14, -4,  68,  68,  -40, -4, -58, 68,  -4,
                                             14, -58, 68, 14, 455, -58, -58, 14,  -4, 5,   -58, -4},
                                            {68, -40, -4, -58, 68,  -4,  14, 68, 5,   14, -4, 68,
                                             68, 14,  14, -58, -58, 455, -4, 5,  -58, 14, -4, -58},
                                            {-58, 14, -4,  5,  -58, -4, -4,  68,  68,  14, 14, 5,
                                             68,  -4, -40, 68, -58, -4, 455, -58, -58, 14, 14, 68},
                                            {-4, 5,  -58, 14, -4, -58, 68,  -4,  -40, 68, -58, -4,
                                             -4, 68, 68,  14, 14, 5,   -58, 455, 14,  68, -58, 14},
                                            {-4,  68, 68, 14, 14, 5,   -58, 14, -4,  5,   -58, -4,
                                             -40, 68, 68, -4, -4, -58, -58, 14, 455, -58, 68,  14},
                                            {68, -4, -40, 68, -58, -4, -4, 5,  -58, 14,  -4, -58,
                                             68, 14, -4,  68, 5,   14, 14, 68, -58, 455, 14, -58},
                                            {68,  14, -4, 68, 5,   14, -40, 68,  68, -4, -4,  -58,
                                             -58, 14, -4, 5,  -58, -4, 14,  -58, 68, 14, 455, -58},
                                            {-40, 68, 68,  -4, -4, -58, 68, 14, -4, 68,  5,   14,
                                             -4,  5,  -58, 14, -4, -58, 68, 14, 14, -58, -58, 455}};

  // Calculate color flows
  jamp[0] = +2. * (+std::complex<double>(0, 1) * amp[0] + std::complex<double>(0, 1) * amp[2] +
                   std::complex<double>(0, 1) * amp[3] - std::complex<double>(0, 1) * amp[5] -
                   std::complex<double>(0, 1) * amp[18] - std::complex<double>(0, 1) * amp[20] -
                   std::complex<double>(0, 1) * amp[21] + std::complex<double>(0, 1) * amp[23] +
                   std::complex<double>(0, 1) * amp[24] - std::complex<double>(0, 1) * amp[29] +
                   std::complex<double>(0, 1) * amp[27] + std::complex<double>(0, 1) * amp[35] +
                   std::complex<double>(0, 1) * amp[34] - std::complex<double>(0, 1) * amp[43] -
                   std::complex<double>(0, 1) * amp[42]);
  jamp[1] = +2. * (+std::complex<double>(0, 1) * amp[1] - std::complex<double>(0, 1) * amp[2] +
                   std::complex<double>(0, 1) * amp[4] + std::complex<double>(0, 1) * amp[5] -
                   std::complex<double>(0, 1) * amp[12] - std::complex<double>(0, 1) * amp[14] -
                   std::complex<double>(0, 1) * amp[15] + std::complex<double>(0, 1) * amp[17] -
                   std::complex<double>(0, 1) * amp[24] + std::complex<double>(0, 1) * amp[29] -
                   std::complex<double>(0, 1) * amp[27] + std::complex<double>(0, 1) * amp[32] +
                   std::complex<double>(0, 1) * amp[31] - std::complex<double>(0, 1) * amp[44] +
                   std::complex<double>(0, 1) * amp[42]);
  jamp[2] = +2. * (-std::complex<double>(0, 1) * amp[0] - std::complex<double>(0, 1) * amp[1] -
                   std::complex<double>(0, 1) * amp[4] - std::complex<double>(0, 1) * amp[3] -
                   std::complex<double>(0, 1) * amp[19] + std::complex<double>(0, 1) * amp[20] -
                   std::complex<double>(0, 1) * amp[22] - std::complex<double>(0, 1) * amp[23] +
                   std::complex<double>(0, 1) * amp[25] - std::complex<double>(0, 1) * amp[32] +
                   std::complex<double>(0, 1) * amp[30] - std::complex<double>(0, 1) * amp[35] -
                   std::complex<double>(0, 1) * amp[34] - std::complex<double>(0, 1) * amp[40] -
                   std::complex<double>(0, 1) * amp[39]);
  jamp[3] = +2. * (+std::complex<double>(0, 1) * amp[1] - std::complex<double>(0, 1) * amp[2] +
                   std::complex<double>(0, 1) * amp[4] + std::complex<double>(0, 1) * amp[5] -
                   std::complex<double>(0, 1) * amp[6] - std::complex<double>(0, 1) * amp[8] -
                   std::complex<double>(0, 1) * amp[9] + std::complex<double>(0, 1) * amp[11] -
                   std::complex<double>(0, 1) * amp[25] + std::complex<double>(0, 1) * amp[29] +
                   std::complex<double>(0, 1) * amp[28] + std::complex<double>(0, 1) * amp[32] -
                   std::complex<double>(0, 1) * amp[30] - std::complex<double>(0, 1) * amp[41] +
                   std::complex<double>(0, 1) * amp[39]);
  jamp[4] = +2. * (-std::complex<double>(0, 1) * amp[0] - std::complex<double>(0, 1) * amp[1] -
                   std::complex<double>(0, 1) * amp[4] - std::complex<double>(0, 1) * amp[3] -
                   std::complex<double>(0, 1) * amp[13] + std::complex<double>(0, 1) * amp[14] -
                   std::complex<double>(0, 1) * amp[16] - std::complex<double>(0, 1) * amp[17] +
                   std::complex<double>(0, 1) * amp[26] - std::complex<double>(0, 1) * amp[32] -
                   std::complex<double>(0, 1) * amp[31] - std::complex<double>(0, 1) * amp[35] +
                   std::complex<double>(0, 1) * amp[33] - std::complex<double>(0, 1) * amp[37] -
                   std::complex<double>(0, 1) * amp[36]);
  jamp[5] = +2. * (+std::complex<double>(0, 1) * amp[0] + std::complex<double>(0, 1) * amp[2] +
                   std::complex<double>(0, 1) * amp[3] - std::complex<double>(0, 1) * amp[5] -
                   std::complex<double>(0, 1) * amp[7] + std::complex<double>(0, 1) * amp[8] -
                   std::complex<double>(0, 1) * amp[10] - std::complex<double>(0, 1) * amp[11] -
                   std::complex<double>(0, 1) * amp[26] - std::complex<double>(0, 1) * amp[29] -
                   std::complex<double>(0, 1) * amp[28] + std::complex<double>(0, 1) * amp[35] -
                   std::complex<double>(0, 1) * amp[33] - std::complex<double>(0, 1) * amp[38] +
                   std::complex<double>(0, 1) * amp[36]);
  jamp[6] = +2. * (+std::complex<double>(0, 1) * amp[6] + std::complex<double>(0, 1) * amp[8] +
                   std::complex<double>(0, 1) * amp[9] - std::complex<double>(0, 1) * amp[11] +
                   std::complex<double>(0, 1) * amp[18] + std::complex<double>(0, 1) * amp[19] +
                   std::complex<double>(0, 1) * amp[22] + std::complex<double>(0, 1) * amp[21] -
                   std::complex<double>(0, 1) * amp[24] - std::complex<double>(0, 1) * amp[28] -
                   std::complex<double>(0, 1) * amp[27] + std::complex<double>(0, 1) * amp[41] +
                   std::complex<double>(0, 1) * amp[40] + std::complex<double>(0, 1) * amp[43] +
                   std::complex<double>(0, 1) * amp[42]);
  jamp[7] = +2. * (+std::complex<double>(0, 1) * amp[7] - std::complex<double>(0, 1) * amp[8] +
                   std::complex<double>(0, 1) * amp[10] + std::complex<double>(0, 1) * amp[11] +
                   std::complex<double>(0, 1) * amp[12] + std::complex<double>(0, 1) * amp[13] +
                   std::complex<double>(0, 1) * amp[16] + std::complex<double>(0, 1) * amp[15] +
                   std::complex<double>(0, 1) * amp[24] + std::complex<double>(0, 1) * amp[28] +
                   std::complex<double>(0, 1) * amp[27] + std::complex<double>(0, 1) * amp[38] +
                   std::complex<double>(0, 1) * amp[37] + std::complex<double>(0, 1) * amp[44] -
                   std::complex<double>(0, 1) * amp[42]);
  jamp[8] = +2. * (-std::complex<double>(0, 1) * amp[6] - std::complex<double>(0, 1) * amp[7] -
                   std::complex<double>(0, 1) * amp[10] - std::complex<double>(0, 1) * amp[9] -
                   std::complex<double>(0, 1) * amp[19] + std::complex<double>(0, 1) * amp[20] -
                   std::complex<double>(0, 1) * amp[22] - std::complex<double>(0, 1) * amp[23] -
                   std::complex<double>(0, 1) * amp[26] - std::complex<double>(0, 1) * amp[34] -
                   std::complex<double>(0, 1) * amp[33] - std::complex<double>(0, 1) * amp[38] +
                   std::complex<double>(0, 1) * amp[36] - std::complex<double>(0, 1) * amp[41] -
                   std::complex<double>(0, 1) * amp[40]);
  jamp[9] = +2. * (-std::complex<double>(0, 1) * amp[0] - std::complex<double>(0, 1) * amp[2] -
                   std::complex<double>(0, 1) * amp[3] + std::complex<double>(0, 1) * amp[5] +
                   std::complex<double>(0, 1) * amp[7] - std::complex<double>(0, 1) * amp[8] +
                   std::complex<double>(0, 1) * amp[10] + std::complex<double>(0, 1) * amp[11] +
                   std::complex<double>(0, 1) * amp[26] + std::complex<double>(0, 1) * amp[29] +
                   std::complex<double>(0, 1) * amp[28] - std::complex<double>(0, 1) * amp[35] +
                   std::complex<double>(0, 1) * amp[33] + std::complex<double>(0, 1) * amp[38] -
                   std::complex<double>(0, 1) * amp[36]);
  jamp[10] = +2. * (-std::complex<double>(0, 1) * amp[6] - std::complex<double>(0, 1) * amp[7] -
                    std::complex<double>(0, 1) * amp[10] - std::complex<double>(0, 1) * amp[9] -
                    std::complex<double>(0, 1) * amp[13] + std::complex<double>(0, 1) * amp[14] -
                    std::complex<double>(0, 1) * amp[16] - std::complex<double>(0, 1) * amp[17] -
                    std::complex<double>(0, 1) * amp[25] - std::complex<double>(0, 1) * amp[31] -
                    std::complex<double>(0, 1) * amp[30] - std::complex<double>(0, 1) * amp[38] -
                    std::complex<double>(0, 1) * amp[37] - std::complex<double>(0, 1) * amp[41] +
                    std::complex<double>(0, 1) * amp[39]);
  jamp[11] = +2. * (-std::complex<double>(0, 1) * amp[1] + std::complex<double>(0, 1) * amp[2] -
                    std::complex<double>(0, 1) * amp[4] - std::complex<double>(0, 1) * amp[5] +
                    std::complex<double>(0, 1) * amp[6] + std::complex<double>(0, 1) * amp[8] +
                    std::complex<double>(0, 1) * amp[9] - std::complex<double>(0, 1) * amp[11] +
                    std::complex<double>(0, 1) * amp[25] - std::complex<double>(0, 1) * amp[29] -
                    std::complex<double>(0, 1) * amp[28] - std::complex<double>(0, 1) * amp[32] +
                    std::complex<double>(0, 1) * amp[30] + std::complex<double>(0, 1) * amp[41] -
                    std::complex<double>(0, 1) * amp[39]);
  jamp[12] = +2. * (+std::complex<double>(0, 1) * amp[12] + std::complex<double>(0, 1) * amp[14] +
                    std::complex<double>(0, 1) * amp[15] - std::complex<double>(0, 1) * amp[17] +
                    std::complex<double>(0, 1) * amp[18] + std::complex<double>(0, 1) * amp[19] +
                    std::complex<double>(0, 1) * amp[22] + std::complex<double>(0, 1) * amp[21] -
                    std::complex<double>(0, 1) * amp[25] - std::complex<double>(0, 1) * amp[31] -
                    std::complex<double>(0, 1) * amp[30] + std::complex<double>(0, 1) * amp[40] +
                    std::complex<double>(0, 1) * amp[39] + std::complex<double>(0, 1) * amp[44] +
                    std::complex<double>(0, 1) * amp[43]);
  jamp[13] = +2. * (+std::complex<double>(0, 1) * amp[6] + std::complex<double>(0, 1) * amp[7] +
                    std::complex<double>(0, 1) * amp[10] + std::complex<double>(0, 1) * amp[9] +
                    std::complex<double>(0, 1) * amp[13] - std::complex<double>(0, 1) * amp[14] +
                    std::complex<double>(0, 1) * amp[16] + std::complex<double>(0, 1) * amp[17] +
                    std::complex<double>(0, 1) * amp[25] + std::complex<double>(0, 1) * amp[31] +
                    std::complex<double>(0, 1) * amp[30] + std::complex<double>(0, 1) * amp[38] +
                    std::complex<double>(0, 1) * amp[37] + std::complex<double>(0, 1) * amp[41] -
                    std::complex<double>(0, 1) * amp[39]);
  jamp[14] = +2. * (-std::complex<double>(0, 1) * amp[12] - std::complex<double>(0, 1) * amp[13] -
                    std::complex<double>(0, 1) * amp[16] - std::complex<double>(0, 1) * amp[15] -
                    std::complex<double>(0, 1) * amp[18] - std::complex<double>(0, 1) * amp[20] -
                    std::complex<double>(0, 1) * amp[21] + std::complex<double>(0, 1) * amp[23] +
                    std::complex<double>(0, 1) * amp[26] + std::complex<double>(0, 1) * amp[34] +
                    std::complex<double>(0, 1) * amp[33] - std::complex<double>(0, 1) * amp[37] -
                    std::complex<double>(0, 1) * amp[36] - std::complex<double>(0, 1) * amp[44] -
                    std::complex<double>(0, 1) * amp[43]);
  jamp[15] = +2. * (+std::complex<double>(0, 1) * amp[0] + std::complex<double>(0, 1) * amp[1] +
                    std::complex<double>(0, 1) * amp[4] + std::complex<double>(0, 1) * amp[3] +
                    std::complex<double>(0, 1) * amp[13] - std::complex<double>(0, 1) * amp[14] +
                    std::complex<double>(0, 1) * amp[16] + std::complex<double>(0, 1) * amp[17] -
                    std::complex<double>(0, 1) * amp[26] + std::complex<double>(0, 1) * amp[32] +
                    std::complex<double>(0, 1) * amp[31] + std::complex<double>(0, 1) * amp[35] -
                    std::complex<double>(0, 1) * amp[33] + std::complex<double>(0, 1) * amp[37] +
                    std::complex<double>(0, 1) * amp[36]);
  jamp[16] = +2. * (-std::complex<double>(0, 1) * amp[7] + std::complex<double>(0, 1) * amp[8] -
                    std::complex<double>(0, 1) * amp[10] - std::complex<double>(0, 1) * amp[11] -
                    std::complex<double>(0, 1) * amp[12] - std::complex<double>(0, 1) * amp[13] -
                    std::complex<double>(0, 1) * amp[16] - std::complex<double>(0, 1) * amp[15] -
                    std::complex<double>(0, 1) * amp[24] - std::complex<double>(0, 1) * amp[28] -
                    std::complex<double>(0, 1) * amp[27] - std::complex<double>(0, 1) * amp[38] -
                    std::complex<double>(0, 1) * amp[37] - std::complex<double>(0, 1) * amp[44] +
                    std::complex<double>(0, 1) * amp[42]);
  jamp[17] = +2. * (-std::complex<double>(0, 1) * amp[1] + std::complex<double>(0, 1) * amp[2] -
                    std::complex<double>(0, 1) * amp[4] - std::complex<double>(0, 1) * amp[5] +
                    std::complex<double>(0, 1) * amp[12] + std::complex<double>(0, 1) * amp[14] +
                    std::complex<double>(0, 1) * amp[15] - std::complex<double>(0, 1) * amp[17] +
                    std::complex<double>(0, 1) * amp[24] - std::complex<double>(0, 1) * amp[29] +
                    std::complex<double>(0, 1) * amp[27] - std::complex<double>(0, 1) * amp[32] -
                    std::complex<double>(0, 1) * amp[31] + std::complex<double>(0, 1) * amp[44] -
                    std::complex<double>(0, 1) * amp[42]);
  jamp[18] = +2. * (+std::complex<double>(0, 1) * amp[12] + std::complex<double>(0, 1) * amp[13] +
                    std::complex<double>(0, 1) * amp[16] + std::complex<double>(0, 1) * amp[15] +
                    std::complex<double>(0, 1) * amp[18] + std::complex<double>(0, 1) * amp[20] +
                    std::complex<double>(0, 1) * amp[21] - std::complex<double>(0, 1) * amp[23] -
                    std::complex<double>(0, 1) * amp[26] - std::complex<double>(0, 1) * amp[34] -
                    std::complex<double>(0, 1) * amp[33] + std::complex<double>(0, 1) * amp[37] +
                    std::complex<double>(0, 1) * amp[36] + std::complex<double>(0, 1) * amp[44] +
                    std::complex<double>(0, 1) * amp[43]);
  jamp[19] = +2. * (+std::complex<double>(0, 1) * amp[6] + std::complex<double>(0, 1) * amp[7] +
                    std::complex<double>(0, 1) * amp[10] + std::complex<double>(0, 1) * amp[9] +
                    std::complex<double>(0, 1) * amp[19] - std::complex<double>(0, 1) * amp[20] +
                    std::complex<double>(0, 1) * amp[22] + std::complex<double>(0, 1) * amp[23] +
                    std::complex<double>(0, 1) * amp[26] + std::complex<double>(0, 1) * amp[34] +
                    std::complex<double>(0, 1) * amp[33] + std::complex<double>(0, 1) * amp[38] -
                    std::complex<double>(0, 1) * amp[36] + std::complex<double>(0, 1) * amp[41] +
                    std::complex<double>(0, 1) * amp[40]);
  jamp[20] = +2. * (-std::complex<double>(0, 1) * amp[12] - std::complex<double>(0, 1) * amp[14] -
                    std::complex<double>(0, 1) * amp[15] + std::complex<double>(0, 1) * amp[17] -
                    std::complex<double>(0, 1) * amp[18] - std::complex<double>(0, 1) * amp[19] -
                    std::complex<double>(0, 1) * amp[22] - std::complex<double>(0, 1) * amp[21] +
                    std::complex<double>(0, 1) * amp[25] + std::complex<double>(0, 1) * amp[31] +
                    std::complex<double>(0, 1) * amp[30] - std::complex<double>(0, 1) * amp[40] -
                    std::complex<double>(0, 1) * amp[39] - std::complex<double>(0, 1) * amp[44] -
                    std::complex<double>(0, 1) * amp[43]);
  jamp[21] = +2. * (+std::complex<double>(0, 1) * amp[0] + std::complex<double>(0, 1) * amp[1] +
                    std::complex<double>(0, 1) * amp[4] + std::complex<double>(0, 1) * amp[3] +
                    std::complex<double>(0, 1) * amp[19] - std::complex<double>(0, 1) * amp[20] +
                    std::complex<double>(0, 1) * amp[22] + std::complex<double>(0, 1) * amp[23] -
                    std::complex<double>(0, 1) * amp[25] + std::complex<double>(0, 1) * amp[32] -
                    std::complex<double>(0, 1) * amp[30] + std::complex<double>(0, 1) * amp[35] +
                    std::complex<double>(0, 1) * amp[34] + std::complex<double>(0, 1) * amp[40] +
                    std::complex<double>(0, 1) * amp[39]);
  jamp[22] = +2. * (-std::complex<double>(0, 1) * amp[6] - std::complex<double>(0, 1) * amp[8] -
                    std::complex<double>(0, 1) * amp[9] + std::complex<double>(0, 1) * amp[11] -
                    std::complex<double>(0, 1) * amp[18] - std::complex<double>(0, 1) * amp[19] -
                    std::complex<double>(0, 1) * amp[22] - std::complex<double>(0, 1) * amp[21] +
                    std::complex<double>(0, 1) * amp[24] + std::complex<double>(0, 1) * amp[28] +
                    std::complex<double>(0, 1) * amp[27] - std::complex<double>(0, 1) * amp[41] -
                    std::complex<double>(0, 1) * amp[40] - std::complex<double>(0, 1) * amp[43] -
                    std::complex<double>(0, 1) * amp[42]);
  jamp[23] = +2. * (-std::complex<double>(0, 1) * amp[0] - std::complex<double>(0, 1) * amp[2] -
                    std::complex<double>(0, 1) * amp[3] + std::complex<double>(0, 1) * amp[5] +
                    std::complex<double>(0, 1) * amp[18] + std::complex<double>(0, 1) * amp[20] +
                    std::complex<double>(0, 1) * amp[21] - std::complex<double>(0, 1) * amp[23] -
                    std::complex<double>(0, 1) * amp[24] + std::complex<double>(0, 1) * amp[29] -
                    std::complex<double>(0, 1) * amp[27] - std::complex<double>(0, 1) * amp[35] -
                    std::complex<double>(0, 1) * amp[34] + std::complex<double>(0, 1) * amp[43] +
                    std::complex<double>(0, 1) * amp[42]);

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
