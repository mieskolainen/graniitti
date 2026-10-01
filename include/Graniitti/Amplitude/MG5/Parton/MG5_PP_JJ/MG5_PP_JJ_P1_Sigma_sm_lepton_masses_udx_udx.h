//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef GRANIITTI_AMPLITUDE_MG5_PP_JJ_MG5_PP_JJ_P1_Sigma_sm_lepton_masses_udx_udx_H
#define GRANIITTI_AMPLITUDE_MG5_PP_JJ_MG5_PP_JJ_P1_Sigma_sm_lepton_masses_udx_udx_H

#include <complex>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_JJ/Parameters_sm_lepton_masses.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_JJ/ProcessBase.h"



//==========================================================================
// A class for calculating the matrix elements for
// Process: u d~ > u d~ QED=0 @1
// Process: u s~ > u s~ QED=0 @1
// Process: u c~ > u c~ QED=0 @1
// Process: d u~ > d u~ QED=0 @1
// Process: d s~ > d s~ QED=0 @1
// Process: d c~ > d c~ QED=0 @1
// Process: s u~ > s u~ QED=0 @1
// Process: s d~ > s d~ QED=0 @1
// Process: s c~ > s c~ QED=0 @1
// Process: c u~ > c u~ QED=0 @1
// Process: c d~ > c d~ QED=0 @1
// Process: c s~ > c s~ QED=0 @1
//--------------------------------------------------------------------------

namespace MG5_PP_JJ {

class MG5_PP_JJ_P1_Sigma_sm_lepton_masses_udx_udx : public ProcessBase
{
  public:
    // Compute the maximum generated amplitude order
    int AlphaQEDPower() const override { return 0; }
    // Compute the maximum generated amplitude order
    int AlphaSPower() const override { return 2; }

    // Constructor.
    MG5_PP_JJ_P1_Sigma_sm_lepton_masses_udx_udx() = default;
    ~MG5_PP_JJ_P1_Sigma_sm_lepton_masses_udx_udx() override = default;

    // Initialize process.
    virtual void initProc(std::string param_card_name);
    // Initialize every mass and coupling through the generated model
    void InitParameters(SLHAReader slha) override;
    // Compute evaluated model particle parameters
    gra::mg5::ParticleMap Particles() const override { return pars.Particles(); }
    // Compute the evaluated UFO electromagnetic coupling
    double AlphaQED() const override { return pars.AlphaQED(); }

    // Calculate flavour-independent parts of cross section.
    virtual void sigmaKin() override;

    // Evaluate sigmaHat(sHat).
    virtual double sigmaHat() override;

    // Compute fixed-basis complex helicity and color components
    std::vector<gra::mg5helas::HelicityComponent> helicityAmplitudes() override;

    // Info on the subprocess.
    virtual std::string name() const {return "u d~ > u d~ (sm_lepton_masses)";}

    virtual int code() const {return 1;}

    const std::vector<double> & getMasses() const override {return mME;}

    // Get and set momenta for matrix element evaluation
    std::vector < double * > getMomenta(){return p;}
    void setMomenta(std::vector < double * > & momenta) override {p = momenta;}
    void setInitial(int inid1, int inid2) override {id1 = inid1; id2 = inid2;}
    void setAlphaS(double in) override {alphaS = in;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2;
    static const int nexternal = 4;
    static const int nprocesses = 2;

  private:

    // Compute the generated subprocess selected by the incoming flavours
    int selectedProcess() const;

    // Compute the identical-flavour multiplicity of the selected subprocess
    double selectedProcessMultiplicity() const;

    // Compute orthogonal complex components of the generated color sum
    std::vector<std::complex<double>> colorAmplitudes() const;

    // Compute raw MG5 leading-color flow amplitudes
    std::vector<std::complex<double>> leadingColorAmplitudes() const;

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]);
    static const int nwavefuncs = 5;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 1;
    std::complex<double> amp[namplitudes];
    double matrix_1_udx_udx();

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses];

    // Color flows, used when selecting color
    std::vector<std::vector<double>> jamp2 =
        std::vector<std::vector<double>>(nprocesses);

    // Pointer to the model parameters
    Parameters_sm_lepton_masses  pars;  // GRANIITTI

    // vector with external particle masses
    std::vector<double> mME;

    // vector with momenta (to be changed each event)
    std::vector < double * > p;
    // Event-dependent strong coupling
    double alphaS = 0.118;

    // Initial particle ids
    int id1, id2;

};


}  // namespace MG5_PP_JJ

#endif  // MG5_Sigma_sm_lepton_masses_udx_udx_H
