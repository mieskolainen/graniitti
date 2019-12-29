//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.6.7, 2019-10-16
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch 
// @@@@ MadGraph to GRANIITTI autoconversion done @@@@
//==========================================================================

#include "Graniitti/Amplitude/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/HelAmps_sm.h"

using namespace MG5_sm; 

//==========================================================================
// Class member functions for calculating the matrix elements for
// Process: a a > w+ w- WEIGHTED<=4 @1

//--------------------------------------------------------------------------
// Initialize process.

void AMP_MG5_yy_ww::initProc(string param_card_name) 
{
  // Instantiate the model class and set parameters that stay fixed during run
  pars = Parameters_sm(); // GRANIITTI 
  SLHAReader slha(param_card_name); 
  pars.setIndependentParameters(slha); 
  pars.setIndependentCouplings(); 
// pars.printIndependentParameters(); // GRANIITTI 
// pars.printIndependentCouplings(); // GRANIITTI 
  // Set external particle masses for this matrix element
  mME.push_back(pars.ZERO); 
  mME.push_back(pars.ZERO); 
  mME.push_back(pars.mdl_MW); 
  mME.push_back(pars.mdl_MW); 
  jamp2[0] = new double[1]; 
}

//--------------------------------------------------------------------------
// Evaluate |M|^2, part independent of incoming flavour.

double AMP_MG5_yy_ww::CalcAmp2(gra::LORENTZSCALAR &lts, double alphas) 
{
  // Set the parameters which change event by event
  pars.setDependentParameters(alphas); 
if (gra::PARAM_STRUCTURE::QED_alpha == "ZERO") 
pars.setAlphaQEDZero(); // GRANIITTI: Set at scale alpha_QED(Q2=0) 

  pars.setDependentCouplings(); 
  static bool firsttime = false; // GRANIITTI 
  if (firsttime)
  {
    pars.printDependentParameters(); 
    pars.printDependentCouplings(); 
    firsttime = false; 
  }

  // Reset color flows
  for(int i = 0; i < 1; i++ )
    jamp2[0][i] = 0.; 

  // Local variables and constants

	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	// *** MADGRAPH CONVENTION IS [E,px,py,pz] ! ***

	// *** Set masses for HELAS ***
	std::vector<double> masses = {0, 0}; // Massless two initial states
	std::vector<gra::M4Vec> pf;

	for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
		masses.push_back(lts.decaytree[i].p4.M());
		pf.push_back(lts.decaytree[i].p4);
	}
	mME = masses;
	
	gra::M4Vec p1_ = lts.q1;
	gra::M4Vec p2_ = lts.q2;
		
	
	// Do kinematic transform
	gra::kinematics::OffShell2LightCone(p1_, p2_, pf);

	// Set initial state 4-momentum
	p.clear();
	double p1[] = {p1_.E(), p1_.Px(), p1_.Py(), p1_.Pz()}; 
	p.push_back(&p1[0]);

	double p2[] = {p2_.E(), p2_.Px(), p2_.Py(), p2_.Pz()};
	p.push_back(&p2[0]);

	// Set final state 4-momentum
	double pthis[pf.size()][4];
	for (std::size_t i = 0; i < pf.size(); ++i) {
		pthis[i][0] = pf[i].E();
		pthis[i][1] = pf[i].Px();
		pthis[i][2] = pf[i].Py();
		pthis[i][3] = pf[i].Pz();
		p.push_back(&pthis[i][0]);
	}
  	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	  const int ncomb = 36; 
lts.hamp = std::vector<std::complex<double>>(ncomb); // GRANIITTI 

  static bool goodhel[ncomb] = {ncomb * false}; 
  static int ntry = 0, sum_hel = 0, ngood = 0; 
  static int igood[ncomb]; 
  static int jhel; 
  std::complex<double> * * wfs; 
  double t[nprocesses]; 
  // Helicities for the process
  static const int helicities[ncomb][nexternal] = {{-1, -1, -1, -1}, {-1, -1,
      -1, 0}, {-1, -1, -1, 1}, {-1, -1, 0, -1}, {-1, -1, 0, 0}, {-1, -1, 0, 1},
      {-1, -1, 1, -1}, {-1, -1, 1, 0}, {-1, -1, 1, 1}, {-1, 1, -1, -1}, {-1, 1,
      -1, 0}, {-1, 1, -1, 1}, {-1, 1, 0, -1}, {-1, 1, 0, 0}, {-1, 1, 0, 1},
      {-1, 1, 1, -1}, {-1, 1, 1, 0}, {-1, 1, 1, 1}, {1, -1, -1, -1}, {1, -1,
      -1, 0}, {1, -1, -1, 1}, {1, -1, 0, -1}, {1, -1, 0, 0}, {1, -1, 0, 1}, {1,
      -1, 1, -1}, {1, -1, 1, 0}, {1, -1, 1, 1}, {1, 1, -1, -1}, {1, 1, -1, 0},
      {1, 1, -1, 1}, {1, 1, 0, -1}, {1, 1, 0, 0}, {1, 1, 0, 1}, {1, 1, 1, -1},
      {1, 1, 1, 0}, {1, 1, 1, 1}};
  // Denominators: spins, colors and identical particles
  const int denominators[nprocesses] = {4}; 

  ntry = 1; // GRANIITTI 

  // Reset the matrix elements
  for(int i = 0; i < nprocesses; i++ )
  {
    matrix_element[i] = 0.; 
  }
  // Define permutation
  int perm[nexternal]; 
  for(int i = 0; i < nexternal; i++ )
  {
    perm[i] = i; 
  }

goto SKIPLABEL; // GRANIITTI: Skip this block 
  if (sum_hel == 0 || ntry < 10)
  {
    // Calculate the matrix element for all helicities
    for(int ihel = 0; ihel < ncomb; ihel++ )
    {
      if (goodhel[ihel] || ntry < 2)
      {
        calculate_wavefunctions(perm, helicities[ihel]); 
        t[0] = matrix_1_aa_wpwm(); 

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
      t[0] = matrix_1_aa_wpwm(); 

      for(int iproc = 0; iproc < nprocesses; iproc++ )
      {
        matrix_element[iproc] += t[iproc] * hwgt; 
      }
    }
  }

  for (int i = 0; i < nprocesses; i++ )
    matrix_element[i] /= denominators[i]; 

	SKIPLABEL:
	
	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
  	// Define permutation
  	for (int i = 0; i < nexternal; ++i) { perm[i] = i; }

  	// Loop over helicity combinations
  	for (int ihel = 0; ihel < ncomb; ++ihel) {
    		calculate_wavefunctions(perm, helicities[ihel]);

    		// Sum of subamplitudes (s,t,u,...)
    		for (int k = 0; k < namplitudes; ++k) { lts.hamp[ihel] += amp[k]; }
	}

	// Total amplitude squared over all helicity combinations individually
	double amp2 = 0.0;
	for (int ihel = 0; ihel < ncomb; ++ihel) {
	amp2 += gra::math::abs2(lts.hamp[ihel]);
	}
	amp2 /= denominators[0];  // spin average matrix element squared

	return amp2;  // amplitude squared
	// @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@
	


}

//--------------------------------------------------------------------------
// Evaluate |M|^2, including incoming flavour dependence.

double AMP_MG5_yy_ww::sigmaHat() 
{
  // Select between the different processes
  if(id1 == 22 && id2 == 22)
  {
    // Add matrix elements for processes with beams (22, 22)
    return matrix_element[0]; 
  }
  else
  {
    // Return 0 if not correct initial state assignment
    return 0.; 
  }
}

//==========================================================================
// Private class member functions

//--------------------------------------------------------------------------
// Evaluate |M|^2 for each subprocess

void AMP_MG5_yy_ww::calculate_wavefunctions(const int perm[], const int hel[])
{
  // Calculate wavefunctions for all processes
  int i, j; 

  // Calculate all wavefunctions
  vxxxxx(p[perm[0]], mME[0], hel[0], -1, w[0]); 
  vxxxxx(p[perm[1]], mME[1], hel[1], -1, w[1]); 
  vxxxxx(p[perm[2]], mME[2], hel[2], +1, w[2]); 
  vxxxxx(p[perm[3]], mME[3], hel[3], +1, w[3]); 
  VVV1_2(w[0], w[2], -pars.GC_3, pars.mdl_MW, pars.mdl_WW, w[4]); 
  VVV1_3(w[0], w[3], -pars.GC_3, pars.mdl_MW, pars.mdl_WW, w[5]); 

  // Calculate all amplitudes
  // Amplitude(s) for diagram number 0
  VVVV2_0(w[0], w[1], w[3], w[2], pars.GC_5, amp[0]); 
  VVV1_0(w[1], w[3], w[4], -pars.GC_3, amp[1]); 
  VVV1_0(w[1], w[5], w[2], -pars.GC_3, amp[2]); 

}
double AMP_MG5_yy_ww::matrix_1_aa_wpwm() 
{
  int i, j; 
  // Local variables
  const int ngraphs = 3; 
  const int ncolor = 1; 
  std::complex<double> ztemp; 
  std::complex<double> jamp[ncolor]; 
  // The color matrix;
  static const double denom[ncolor] = {1}; 
  static const double cf[ncolor][ncolor] = {{1}}; 

  // Calculate color flows
  jamp[0] = +amp[0] + amp[1] + amp[2]; 

  // Sum and square the color flows to get the matrix element
  double matrix = 0; 
  for(i = 0; i < ncolor; i++ )
  {
    ztemp = 0.; 
    for(j = 0; j < ncolor; j++ )
      ztemp = ztemp + cf[i][j] * jamp[j]; 
    matrix = matrix + real(ztemp * conj(jamp[i]))/denom[i]; 
  }

  // Store the leading color flows for choice of color
  for(i = 0; i < ncolor; i++ )
    jamp2[0][i] += real(jamp[i] * conj(jamp[i])); 

  return matrix; 
}



