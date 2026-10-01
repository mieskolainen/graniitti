//==========================================================================
// This file has been automatically generated for C++ by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#include <cmath>
#include <iostream>
#include <iomanip>
#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/Parameters_sm_lepton_masses.h"




void Parameters_sm_lepton_masses::setIndependentParameters(SLHAReader& slha)
{
  // Define "zero"
  zero = 0;
  ZERO = 0;
  // Prepare a vector for indices
  std::vector<int> indices(2, 0);
  mdl_WH = slha.get_block_entry("decay", 25, 6.382339e-03);
  mdl_WW = slha.get_block_entry("decay", 24, 2.047600e+00);
  mdl_WZ = slha.get_block_entry("decay", 23, 2.441404e+00);
  mdl_WTau = slha.get_block_entry("decay", 15, 2.270000e-12);
  mdl_WT = slha.get_block_entry("decay", 6, 1.491500e+00);
  mdl_ymtau = slha.get_block_entry("yukawa", 15, 1.777000e+00);
  mdl_ymm = slha.get_block_entry("yukawa", 13, 1.056600e-01);
  mdl_yme = slha.get_block_entry("yukawa", 11, 5.110000e-04);
  mdl_ymt = slha.get_block_entry("yukawa", 6, 1.730000e+02);
  mdl_ymb = slha.get_block_entry("yukawa", 5, 4.700000e+00);
  aS = slha.get_block_entry("sminputs", 3, 1.180000e-01);
  mdl_Gf = slha.get_block_entry("sminputs", 2, 1.166390e-05);
  aEWM1 = slha.get_block_entry("sminputs", 1, 1.325070e+02);
  mdl_MH = slha.get_block_entry("mass", 25, 1.250000e+02);
  mdl_MZ = slha.get_block_entry("mass", 23, 9.118800e+01);
  mdl_MTA = slha.get_block_entry("mass", 15, 1.777000e+00);
  mdl_MM = slha.get_block_entry("mass", 13, 1.056600e-01);
  mdl_Me = slha.get_block_entry("mass", 11, 5.110000e-04);
  mdl_MT = slha.get_block_entry("mass", 6, 1.730000e+02);
  mdl_MB = slha.get_block_entry("mass", 5, 4.700000e+00);
  mdl_conjg__CKM3x3 = 1.;
  mdl_CKM3x3 = 1.;
  mdl_conjg__CKM1x1 = 1.;
  mdl_complexi = std::complex<double> (0., 1.);
  mdl_MZ__exp__2 = ((mdl_MZ) * (mdl_MZ));
  mdl_MZ__exp__4 = ((mdl_MZ) * (mdl_MZ) * (mdl_MZ) * (mdl_MZ));
  mdl_sqrt__2 = std::sqrt(2.);
  mdl_MH__exp__2 = ((mdl_MH) * (mdl_MH));
  mdl_aEW = 1./aEWM1;
  mdl_MW = std::sqrt(mdl_MZ__exp__2/2. + std::sqrt(mdl_MZ__exp__4/4. - (mdl_aEW * M_PI *
      mdl_MZ__exp__2)/(mdl_Gf * mdl_sqrt__2)));
  mdl_sqrt__aEW = std::sqrt(mdl_aEW);
  mdl_ee = 2. * mdl_sqrt__aEW * std::sqrt(M_PI);
  mdl_MW__exp__2 = ((mdl_MW) * (mdl_MW));
  mdl_sw2 = 1. - mdl_MW__exp__2/mdl_MZ__exp__2;
  mdl_cw = std::sqrt(1. - mdl_sw2);
  mdl_sqrt__sw2 = std::sqrt(mdl_sw2);
  mdl_sw = mdl_sqrt__sw2;
  mdl_g1 = mdl_ee/mdl_cw;
  mdl_gw = mdl_ee/mdl_sw;
  mdl_vev = (2. * mdl_MW * mdl_sw)/mdl_ee;
  mdl_vev__exp__2 = ((mdl_vev) * (mdl_vev));
  mdl_lam = mdl_MH__exp__2/(2. * mdl_vev__exp__2);
  mdl_yb = (mdl_ymb * mdl_sqrt__2)/mdl_vev;
  mdl_ye = (mdl_yme * mdl_sqrt__2)/mdl_vev;
  mdl_ym = (mdl_ymm * mdl_sqrt__2)/mdl_vev;
  mdl_yt = (mdl_ymt * mdl_sqrt__2)/mdl_vev;
  mdl_ytau = (mdl_ymtau * mdl_sqrt__2)/mdl_vev;
  mdl_muH = std::sqrt(mdl_lam * mdl_vev__exp__2);
  mdl_I1x33 = mdl_yb * mdl_conjg__CKM3x3;
  mdl_I2x33 = mdl_yt * mdl_conjg__CKM3x3;
  mdl_I3x33 = mdl_CKM3x3 * mdl_yt;
  mdl_I4x33 = mdl_CKM3x3 * mdl_yb;
  mdl_ee__exp__2 = ((mdl_ee) * (mdl_ee));
  mdl_sw__exp__2 = ((mdl_sw) * (mdl_sw));
  mdl_cw__exp__2 = ((mdl_cw) * (mdl_cw));
}
void Parameters_sm_lepton_masses::setIndependentCouplings()
{
  GC_2 = (2. * mdl_ee * mdl_complexi)/3.;
  GC_3 = -(mdl_ee * mdl_complexi);
}
void Parameters_sm_lepton_masses::setDependentParameters(double alpS)
{
  aS = alpS;
  mdl_sqrt__aS = std::sqrt(aS);
  G = 2. * mdl_sqrt__aS * std::sqrt(M_PI);
  mdl_G__exp__2 = ((G) * (G));
}
void Parameters_sm_lepton_masses::setDependentCouplings()
{
  GC_10 = -G;
  GC_11 = mdl_complexi * G;
  GC_12 = mdl_complexi * mdl_G__exp__2;
}
void Parameters_sm_lepton_masses::setAlphaQEDZero() {
  const double mdl_ee_NEW = 2.0 * std::sqrt(1.0 / 137.03599908000) * std::sqrt(M_PI);

  GC_2 = (2. * mdl_ee_NEW * mdl_complexi)/3.;
  GC_3 = -(mdl_ee_NEW * mdl_complexi);
}

// Routines for printing out parameters
void Parameters_sm_lepton_masses::printIndependentParameters()
{
  std::cout <<  "sm_lepton_masses model parameters independent of event kinematics:"
      << std::endl;
  std::cout << std::setw(20) <<  "mdl_WH " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_WH << std::endl;
  std::cout << std::setw(20) <<  "mdl_WW " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_WW << std::endl;
  std::cout << std::setw(20) <<  "mdl_WZ " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_WZ << std::endl;
  std::cout << std::setw(20) <<  "mdl_WTau " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_WTau << std::endl;
  std::cout << std::setw(20) <<  "mdl_WT " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_WT << std::endl;
  std::cout << std::setw(20) <<  "mdl_ymtau " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ymtau << std::endl;
  std::cout << std::setw(20) <<  "mdl_ymm " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ymm << std::endl;
  std::cout << std::setw(20) <<  "mdl_yme " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_yme << std::endl;
  std::cout << std::setw(20) <<  "mdl_ymt " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ymt << std::endl;
  std::cout << std::setw(20) <<  "mdl_ymb " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ymb << std::endl;
  std::cout << std::setw(20) <<  "aS " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << aS << std::endl;
  std::cout << std::setw(20) <<  "mdl_Gf " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_Gf << std::endl;
  std::cout << std::setw(20) <<  "aEWM1 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << aEWM1 << std::endl;
  std::cout << std::setw(20) <<  "mdl_MH " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MH << std::endl;
  std::cout << std::setw(20) <<  "mdl_MZ " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MZ << std::endl;
  std::cout << std::setw(20) <<  "mdl_MTA " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MTA << std::endl;
  std::cout << std::setw(20) <<  "mdl_MM " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MM << std::endl;
  std::cout << std::setw(20) <<  "mdl_Me " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_Me << std::endl;
  std::cout << std::setw(20) <<  "mdl_MT " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MT << std::endl;
  std::cout << std::setw(20) <<  "mdl_MB " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MB << std::endl;
  std::cout << std::setw(20) <<  "mdl_conjg__CKM3x3 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_conjg__CKM3x3 << std::endl;
  std::cout << std::setw(20) <<  "mdl_CKM3x3 " <<  "= " << std::setiosflags(std::ios::scientific)
      << std::setw(10) << mdl_CKM3x3 << std::endl;
  std::cout << std::setw(20) <<  "mdl_conjg__CKM1x1 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_conjg__CKM1x1 << std::endl;
  std::cout << std::setw(20) <<  "mdl_complexi " <<  "= " << std::setiosflags(std::ios::scientific)
      << std::setw(10) << mdl_complexi << std::endl;
  std::cout << std::setw(20) <<  "mdl_MZ__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_MZ__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_MZ__exp__4 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_MZ__exp__4 << std::endl;
  std::cout << std::setw(20) <<  "mdl_sqrt__2 " <<  "= " << std::setiosflags(std::ios::scientific)
      << std::setw(10) << mdl_sqrt__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_MH__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_MH__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_aEW " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_aEW << std::endl;
  std::cout << std::setw(20) <<  "mdl_MW " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_MW << std::endl;
  std::cout << std::setw(20) <<  "mdl_sqrt__aEW " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_sqrt__aEW << std::endl;
  std::cout << std::setw(20) <<  "mdl_ee " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ee << std::endl;
  std::cout << std::setw(20) <<  "mdl_MW__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_MW__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_sw2 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_sw2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_cw " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_cw << std::endl;
  std::cout << std::setw(20) <<  "mdl_sqrt__sw2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_sqrt__sw2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_sw " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_sw << std::endl;
  std::cout << std::setw(20) <<  "mdl_g1 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_g1 << std::endl;
  std::cout << std::setw(20) <<  "mdl_gw " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_gw << std::endl;
  std::cout << std::setw(20) <<  "mdl_vev " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_vev << std::endl;
  std::cout << std::setw(20) <<  "mdl_vev__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_vev__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_lam " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_lam << std::endl;
  std::cout << std::setw(20) <<  "mdl_yb " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_yb << std::endl;
  std::cout << std::setw(20) <<  "mdl_ye " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ye << std::endl;
  std::cout << std::setw(20) <<  "mdl_ym " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ym << std::endl;
  std::cout << std::setw(20) <<  "mdl_yt " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_yt << std::endl;
  std::cout << std::setw(20) <<  "mdl_ytau " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_ytau << std::endl;
  std::cout << std::setw(20) <<  "mdl_muH " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_muH << std::endl;
  std::cout << std::setw(20) <<  "mdl_I1x33 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_I1x33 << std::endl;
  std::cout << std::setw(20) <<  "mdl_I2x33 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_I2x33 << std::endl;
  std::cout << std::setw(20) <<  "mdl_I3x33 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_I3x33 << std::endl;
  std::cout << std::setw(20) <<  "mdl_I4x33 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << mdl_I4x33 << std::endl;
  std::cout << std::setw(20) <<  "mdl_ee__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_ee__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_sw__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_sw__exp__2 << std::endl;
  std::cout << std::setw(20) <<  "mdl_cw__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_cw__exp__2 << std::endl;
}
void Parameters_sm_lepton_masses::printIndependentCouplings()
{
  std::cout <<  "sm_lepton_masses model couplings independent of event kinematics:"
      << std::endl;
  std::cout << std::setw(20) <<  "GC_2 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << GC_2 << std::endl;
  std::cout << std::setw(20) <<  "GC_3 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << GC_3 << std::endl;
}
void Parameters_sm_lepton_masses::printDependentParameters()
{
  std::cout <<  "sm_lepton_masses model parameters dependent on event kinematics:"
      << std::endl;
  std::cout << std::setw(20) <<  "mdl_sqrt__aS " <<  "= " << std::setiosflags(std::ios::scientific)
      << std::setw(10) << mdl_sqrt__aS << std::endl;
  std::cout << std::setw(20) <<  "G " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << G << std::endl;
  std::cout << std::setw(20) <<  "mdl_G__exp__2 " <<  "= " <<
      std::setiosflags(std::ios::scientific) << std::setw(10) << mdl_G__exp__2 << std::endl;
}
void Parameters_sm_lepton_masses::printDependentCouplings()
{
  std::cout <<  "sm_lepton_masses model couplings dependent on event kinematics:" <<
      std::endl;
  std::cout << std::setw(20) <<  "GC_10 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << GC_10 << std::endl;
  std::cout << std::setw(20) <<  "GC_11 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << GC_11 << std::endl;
  std::cout << std::setw(20) <<  "GC_12 " <<  "= " << std::setiosflags(std::ios::scientific) <<
      std::setw(10) << GC_12 << std::endl;
}
