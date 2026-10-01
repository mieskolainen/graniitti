//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
//==========================================================================

#ifndef GRANIITTI_AMPLITUDE_MG5_PP_ZJ_MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg_H
#define GRANIITTI_AMPLITUDE_MG5_PP_ZJ_MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg_H

#include <complex>
#include <vector>

#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/Parameters_sm_lepton_masses.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/ProcessBase.h"



//==========================================================================
// A class for calculating the matrix elements for
// Process: u u~ > z g WEIGHTED<=3 @3
// *   Decay: z > ta+ ta- WEIGHTED<=2
// Process: c c~ > z g WEIGHTED<=3 @3
// *   Decay: z > ta+ ta- WEIGHTED<=2
//--------------------------------------------------------------------------

namespace MG5_PP_ZJ {

class MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg : public ProcessBase
{
  public:
    // Compute the maximum generated amplitude order
    int AlphaQEDPower() const override { return 2; }
    // Compute the maximum generated amplitude order
    int AlphaSPower() const override { return 1; }

    // Constructor.
    MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg() = default;
    ~MG5_PP_ZJ_P3_Sigma_sm_lepton_masses_uux_taptamg() override = default;

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
    virtual std::string name() const {return "u u~ > ta+ ta- g (sm_lepton_masses)";}

    virtual int code() const {return 3;}

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
    static const int nexternal = 5;
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
    static const int nwavefuncs = 8;
    std::complex<double> w[nwavefuncs][18];
    static const int namplitudes = 2;
    std::complex<double> amp[namplitudes];
    double matrix_3_uux_zg_z_taptam();

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


}  // namespace MG5_PP_ZJ

#endif  // MG5_Sigma_sm_lepton_masses_uux_taptamg_H
