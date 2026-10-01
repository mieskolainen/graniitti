//==========================================================================
// This file has been automatically generated for C++
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef GRANIITTI_AMPLITUDE_MG5_YY_ZJJ_PARAMETERS_SM_LEPTON_MASSES_H
#define GRANIITTI_AMPLITUDE_MG5_YY_ZJJ_PARAMETERS_SM_LEPTON_MASSES_H

#include <complex>
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"

#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"


namespace MG5_YY_ZJJ {

class Parameters_sm_lepton_masses
{
  public:
    // Compute the evaluated UFO electromagnetic coupling
    double AlphaQED() const { return std::norm(mdl_ee) / (4.0 * M_PI); }

    // Compute masses and widths from the evaluated UFO parameters
    gra::mg5::ParticleMap Particles() const {
      return {
      {22, {ZERO, ZERO}},
      {23, {mdl_MZ, mdl_WZ}},
      {24, {mdl_MW, mdl_WW}},
      {21, {ZERO, ZERO}},
      {12, {ZERO, ZERO}},
      {14, {ZERO, ZERO}},
      {16, {ZERO, ZERO}},
      {2, {ZERO, ZERO}},
      {4, {ZERO, ZERO}},
      {6, {mdl_MT, mdl_WT}},
      {1, {ZERO, ZERO}},
      {3, {ZERO, ZERO}},
      {5, {mdl_MB, ZERO}},
      {25, {mdl_MH, mdl_WH}},
      {11, {mdl_Me, ZERO}},
      {13, {mdl_MM, ZERO}},
      {15, {mdl_MTA, mdl_WTau}}
      };
    }


    // Define "zero"
    double zero, ZERO;
    // Model parameters independent of aS
    double mdl_WH, mdl_WW, mdl_WZ, mdl_WTau, mdl_WT, mdl_ymtau, mdl_ymm,
        mdl_yme, mdl_ymt, mdl_ymb, aS, mdl_Gf, aEWM1, mdl_MH, mdl_MZ, mdl_MTA,
        mdl_MM, mdl_Me, mdl_MT, mdl_MB, mdl_conjg__CKM3x3, mdl_CKM3x3,
        mdl_conjg__CKM1x1, mdl_MZ__exp__2, mdl_MZ__exp__4, mdl_sqrt__2,
        mdl_MH__exp__2, mdl_aEW, mdl_MW, mdl_sqrt__aEW, mdl_ee, mdl_MW__exp__2,
        mdl_sw2, mdl_cw, mdl_sqrt__sw2, mdl_sw, mdl_g1, mdl_gw, mdl_vev,
        mdl_vev__exp__2, mdl_lam, mdl_yb, mdl_ye, mdl_ym, mdl_yt, mdl_ytau,
        mdl_muH, mdl_ee__exp__2, mdl_sw__exp__2, mdl_cw__exp__2;
    std::complex<double> mdl_complexi, mdl_I1x33, mdl_I2x33, mdl_I3x33,
        mdl_I4x33;
    // Model parameters dependent on aS
    double mdl_sqrt__aS, G, mdl_G__exp__2;
    // Model couplings independent of aS
    std::complex<double> GC_1, GC_2, GC_50, GC_58, GC_59;
    // Model couplings dependent on aS


    // Set parameters that are unchanged during the run
    void setIndependentParameters(SLHAReader& slha);
    // Set couplings that are unchanged during the run
    void setIndependentCouplings();
    // Set parameters that are changed event by event
    void setDependentParameters(double alpS);
    // Set couplings that are changed event by event
    void setDependentCouplings();
    // Set electromagnetic couplings at Q2 = 0
    void setAlphaQEDZero();

    // Print parameters that are unchanged during the run
    void printIndependentParameters();
    // Print couplings that are unchanged during the run
    void printIndependentCouplings();
    // Print parameters that are changed event by event
    void printDependentParameters();
    // Print couplings that are changed event by event
    void printDependentCouplings();


  private:

};

}  // namespace MG5_YY_ZJJ

#endif  // Parameters_sm_lepton_masses_H
