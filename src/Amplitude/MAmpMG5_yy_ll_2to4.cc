//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 3.0.0, 2018-05-01
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================
// Modified from the MadGraph C++ Standalone output
// Mikael Mieskolainen, 2018

// [FEASABILITY TEST CLASS]

// Own
#include "Graniitti/Amplitude/MAmpMG5_yy_ll_2to4.h"
#include "Graniitti/Amplitude/HelAmps_sm_yy_ll_2to4.h"
#include "Graniitti/Amplitude/Parameters_sm.h"

#include "Graniitti/MAux.h"
#include "Graniitti/MForm.h"
#include "Graniitti/MKinematics.h"
#include "Graniitti/MMath.h"

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: mu+ mu+ > mu+ e+ e- mu+ WEIGHTED<=8 / g z mu @1

//--------------------------------------------------------------------------
// Initialize process.

using namespace MG5_sm_yy_ll_2to4;

MAmpMG5_yy_ll_2to4::MAmpMG5_yy_ll_2to4() {
  std::string param_card_name = gra::aux::GetBasePath(2) + "/MG5cards/" + "yy_ll_param_card.dat";

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
  mME.push_back(pars.ZERO);

  jamp2 = std::vector<std::vector<double>>(nprocesses, std::vector<double>(2));
}

MAmpMG5_yy_ll_2to4::~MAmpMG5_yy_ll_2to4() {}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.
//
//
// 2 u -----------> u  4
//          |
//          |-----> e+ 6
//          |
//          |-----> e- 5
//          |
// 1 u -----------> u  3

// Get amplitude squared
double MAmpMG5_yy_ll_2to4::CalcAmp2(gra::LORENTZSCALAR &lts) {
  // Set the parameters which change event by event
  // pars.setDependentParameters();
  // pars.setDependentCouplings();
  // static bool firsttime = true;
  // if (firsttime)
  //{
  //  pars.printDependentParameters();
  //  pars.printDependentCouplings();
  //  firsttime = false;
  //}

  // *** Set masses for HELAS ***
  const std::vector<double> masses = {lts.pbeam2.M(),          lts.pbeam1.M(),
                                      lts.pfinal[1].M(),       lts.pfinal[0].M(),
                                      lts.decaytree[0].p4.M(), lts.decaytree[1].p4.M()};
  mME = masses;

  // *** Set particle 4-momentum: [E,px,py,pz] convention here! ***
  double p1[] = {lts.pbeam2.E(), lts.pbeam2.Px(), lts.pbeam2.Py(), lts.pbeam2.Pz()};
  double p2[] = {lts.pbeam1.E(), lts.pbeam1.Px(), lts.pbeam1.Py(), lts.pbeam1.Pz()};
  double p3[] = {lts.pfinal[1].E(), lts.pfinal[1].Px(), lts.pfinal[1].Py(), lts.pfinal[1].Pz()};
  double p4[] = {lts.pfinal[0].E(), lts.pfinal[0].Px(), lts.pfinal[0].Py(), lts.pfinal[0].Pz()};
  double p6[] = {lts.decaytree[0].p4.E(), lts.decaytree[0].p4.Px(), lts.decaytree[0].p4.Py(),
                 lts.decaytree[0].p4.Pz()};
  double p5[] = {lts.decaytree[1].p4.E(), lts.decaytree[1].p4.Px(), lts.decaytree[1].p4.Py(),
                 lts.decaytree[1].p4.Pz()};

  p.clear();
  p.push_back(&p1[0]);
  p.push_back(&p2[0]);
  p.push_back(&p3[0]);
  p.push_back(&p4[0]);
  p.push_back(&p5[0]);
  p.push_back(&p6[0]);

  // Local variables and constants
  const int                         ncomb = 64;
  std::vector<std::complex<double>> temp(ncomb, 0.0);
  lts.hamp = temp;

  // static bool goodhel[ncomb] = {ncomb * false};
  // static int ntry = 0, sum_hel = 0, ngood = 0;
  // static int igood[ncomb];
  // static int jhel;
  // double t[nprocesses];

  static const int helicities[ncomb][nexternal] = {
      {-1, -1, -1, -1, -1, -1}, {-1, -1, -1, -1, -1, 1}, {-1, -1, -1, -1, 1, -1},
      {-1, -1, -1, -1, 1, 1},   {-1, -1, -1, 1, -1, -1}, {-1, -1, -1, 1, -1, 1},
      {-1, -1, -1, 1, 1, -1},   {-1, -1, -1, 1, 1, 1},   {-1, -1, 1, -1, -1, -1},
      {-1, -1, 1, -1, -1, 1},   {-1, -1, 1, -1, 1, -1},  {-1, -1, 1, -1, 1, 1},
      {-1, -1, 1, 1, -1, -1},   {-1, -1, 1, 1, -1, 1},   {-1, -1, 1, 1, 1, -1},
      {-1, -1, 1, 1, 1, 1},     {-1, 1, -1, -1, -1, -1}, {-1, 1, -1, -1, -1, 1},
      {-1, 1, -1, -1, 1, -1},   {-1, 1, -1, -1, 1, 1},   {-1, 1, -1, 1, -1, -1},
      {-1, 1, -1, 1, -1, 1},    {-1, 1, -1, 1, 1, -1},   {-1, 1, -1, 1, 1, 1},
      {-1, 1, 1, -1, -1, -1},   {-1, 1, 1, -1, -1, 1},   {-1, 1, 1, -1, 1, -1},
      {-1, 1, 1, -1, 1, 1},     {-1, 1, 1, 1, -1, -1},   {-1, 1, 1, 1, -1, 1},
      {-1, 1, 1, 1, 1, -1},     {-1, 1, 1, 1, 1, 1},     {1, -1, -1, -1, -1, -1},
      {1, -1, -1, -1, -1, 1},   {1, -1, -1, -1, 1, -1},  {1, -1, -1, -1, 1, 1},
      {1, -1, -1, 1, -1, -1},   {1, -1, -1, 1, -1, 1},   {1, -1, -1, 1, 1, -1},
      {1, -1, -1, 1, 1, 1},     {1, -1, 1, -1, -1, -1},  {1, -1, 1, -1, -1, 1},
      {1, -1, 1, -1, 1, -1},    {1, -1, 1, -1, 1, 1},    {1, -1, 1, 1, -1, -1},
      {1, -1, 1, 1, -1, 1},     {1, -1, 1, 1, 1, -1},    {1, -1, 1, 1, 1, 1},
      {1, 1, -1, -1, -1, -1},   {1, 1, -1, -1, -1, 1},   {1, 1, -1, -1, 1, -1},
      {1, 1, -1, -1, 1, 1},     {1, 1, -1, 1, -1, -1},   {1, 1, -1, 1, -1, 1},
      {1, 1, -1, 1, 1, -1},     {1, 1, -1, 1, 1, 1},     {1, 1, 1, -1, -1, -1},
      {1, 1, 1, -1, -1, 1},     {1, 1, 1, -1, 1, -1},    {1, 1, 1, -1, 1, 1},
      {1, 1, 1, 1, -1, -1},     {1, 1, 1, 1, -1, 1},     {1, 1, 1, 1, 1, -1},
      {1, 1, 1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {8};

  // ntry = ntry + 1;

  // Reset the matrix elements
  // for (int i = 0; i < nprocesses; ++i ) {
  //  matrix_element[i] = 0.;
  //}

  // Define permutation
  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  // const static std::vector<int> nonzero =
  // {0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15}; // All
  std::vector<int> nonzero;  // High energy limit / Helicity conserving
  if (nonzero.size() == 0) {
    for (unsigned int i = 0; i < ncomb; ++i) { nonzero.push_back(i); }
  }

  // Loop over helicity combinations
  for (auto ihel : nonzero) {
    if ((helicities[ihel][1 - 1] == helicities[ihel][3 - 1]) &&
        (helicities[ihel][2 - 1] == helicities[ihel][4 - 1])) {  // no flip
      calculate_wavefunctions(perm, helicities[ihel]);

      // Sum of subamplitudes (s,t,u,...)
      for (int k = 0; k < namplitudes; ++k) { lts.hamp[ihel] += amp[k]; }
    }
  }

  // Total amplitude squared over all helicity combinations individually
  double amp2 = 0.0;
  for (auto ihel : nonzero) {
    // printf("%4d : %0.5E \n", ihel, gra::math::abs2(lts.hamp[ihel]));
    amp2 += gra::math::abs2(lts.hamp[ihel]);
  }
  amp2 /= denominators[0];  // spin average matrix element squared

  // printf("\n");

  return amp2;  // amplitude squared

  /*
  if (sum_hel == 0 || ntry < 10) {

    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ ) {

          if (goodhel[ihel] || ntry < 2) {

            calculate_wavefunctions(perm, helicities[ihel]);
            t[0] = matrix_1();

            double tsum = 0;
            for(int iproc = 0; iproc < nprocesses; iproc++ )
            {
                  matrix_element[iproc] += t[iproc];
                  tsum += t[iproc];
            }
            // Store which helicities give non-zero result
            if (tsum != 0. && !goodhel[ihel])
            {
                  goodhel[ihel] = true;
                  ngood++;
                  igood[ngood] = ihel;
            }
          }
    }
    jhel = 0;
    sum_hel = min(sum_hel, ngood);
  }
  else
  {
    // Only use the "good" helicities
    for(int j = 0; j < sum_hel; j++ )
    {
          jhel++;
          if (jhel >= ngood)
            jhel = 0;
          double hwgt = double(ngood)/double(sum_hel);
          int ihel = igood[jhel];
          calculate_wavefunctions(perm, helicities[ihel]);
          t[0] = matrix_1();

          for(int iproc = 0; iproc < nprocesses; iproc++ )
          {
            matrix_element[iproc] += t[iproc] * hwgt;
          }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i];


  return matrix_element[0];
  */
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void MAmpMG5_yy_ll_2to4::calculate_wavefunctions(const int perm[], const int hel[]) {
  // *** MODIFIED: Constant QED coupling, instead of GC_3 ***
  const std::complex<double> GC_3(0, -gra::form::e_EM());

  // Calculate all wavefunctions
  oxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]);
  oxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]);
  ixxxxx(p[perm[2]], mME[2], hel[2], -1, w[2]);
  ixxxxx(p[perm[3]], mME[3], hel[3], -1, w[3]);
  ixxxxx(p[perm[4]], mME[4], hel[4], -1, w[4]);
  oxxxxx(p[perm[5]], mME[5], hel[5], +1, w[5]);
  FFV1P0_3(w[2], w[0], GC_3, pars.ZERO, pars.ZERO, w[6]);
  FFV1P0_3(w[3], w[1], GC_3, pars.ZERO, pars.ZERO, w[7]);
  FFV1_2(w[4], w[6], GC_3, pars.ZERO, pars.ZERO, w[8]);
  FFV1_1(w[5], w[6], GC_3, pars.ZERO, pars.ZERO, w[9]);
  FFV1P0_3(w[3], w[0], GC_3, pars.ZERO, pars.ZERO, w[10]);
  FFV1P0_3(w[2], w[1], GC_3, pars.ZERO, pars.ZERO, w[11]);
  FFV1_2(w[4], w[10], GC_3, pars.ZERO, pars.ZERO, w[12]);
  FFV1_1(w[5], w[10], GC_3, pars.ZERO, pars.ZERO, w[13]);

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  FFV1_0(w[8], w[5], w[7], GC_3, amp[0]);
  FFV1_0(w[4], w[9], w[7], GC_3, amp[1]);
  FFV1_0(w[12], w[5], w[11], GC_3, amp[2]);
  FFV1_0(w[4], w[13], w[11], GC_3, amp[3]);
}

double MAmpMG5_yy_ll_2to4::matrix_1() {
  int i, j;
  // Local variables
  // const int ngraphs = 4;
  const int            ncolor = 1;
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  // The color matrix;
  static const double denom[ncolor]      = {1};
  static const double cf[ncolor][ncolor] = {{1}};

  // Calculate color flows
  jamp[0] = -amp[0] - amp[1] + amp[2] + amp[3];

  // Sum and square the color flows to get the matrix element
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
