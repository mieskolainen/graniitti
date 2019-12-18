//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.6.7, 2019-10-16
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch 
// @@@@ MadGraph to GRANIITTI autoconversion done @@@@
//==========================================================================

#ifndef MG5_Sigma_sm_gg_ggg_H
#define MG5_Sigma_sm_gg_ggg_H

#include <complex> 
#include <vector> 

#include "Graniitti/Amplitude/Parameters_sm.h" 
#include "Graniitti/MForm.h" 
#include "Graniitti/MAux.h" 
#include "Graniitti/MKinematics.h" 
#include "Graniitti/MMath.h"

using namespace std; 

//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > g g g WEIGHTED<=3 @1
//--------------------------------------------------------------------------

class AMP_MG5_gg_ggg
{
  public:

    // Constructor.
    AMP_MG5_gg_ggg() {initProc(gra::aux::GetBasePath(2) + "/MG5cards/param_card_gg_ggg.dat"); }

    // Initialize process.
     void initProc(string param_card_name); 

    // Calculate flavour-independent parts of cross section.
     double CalcAmp2(gra::LORENTZSCALAR &lts, double alphas); 

    // Evaluate sigmaHat(sHat).
     double sigmaHat(); 

    // Info on the subprocess.
     string name() const {return "g g > g g g (sm)";}

     int code() const {return 1;}

    const vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    vector < double * > getMomenta(){return p;}
    void setMomenta(vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
    static const int ninitial = 2; 
    static const int nexternal = 5; 
    static const int nprocesses = 1; 

  private:

    // Private functions to calculate the matrix element for all subprocesses
    // Calculate wavefunctions
    void calculate_wavefunctions(const int perm[], const int hel[]); 
    static const int nwavefuncs = 33; 
    std::complex<double> w[nwavefuncs][18]; 
    static const int namplitudes = 45; 
    std::complex<double> amp[namplitudes]; 
    double matrix_1_gg_ggg(); 

    // Store the matrix element value from sigmaKin
    double matrix_element[nprocesses]; 

    // Color flows, used when selecting color
    double * jamp2[nprocesses]; 

    // Pointer to the model parameters
    Parameters_sm pars; // GRANIITTI 

    // vector with external particle masses
    vector<double> mME; 

    // vector with momenta (to be changed each event)
    vector < double * > p; 
    // Initial particle ids
    int id1, id2; 

}; 


#endif  // MG5_Sigma_sm_gg_ggg_H
