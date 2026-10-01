//==========================================================================
// This file has been automatically generated for C++ Standalone by
// MadGraph5_aMC@NLO v. 2.9.27, 2026-01-05
// By the MadGraph5_aMC@NLO Development Team
// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch
// @@@@ MadGraph to GRANIITTI conversion done @@@@
//==========================================================================

#ifndef MG5_Sigma_sm_lepton_masses_gg_uuxg_H
#define MG5_Sigma_sm_lepton_masses_gg_uuxg_H

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <stdexcept>
#include <vector>

#include "Graniitti/Amplitude/MG5/Runtime/Models/sm/Parameters_sm_lepton_masses.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Kinematics.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Math/MMath.h"



//==========================================================================
// A class for calculating the matrix elements for
// Process: g g > u u~ g QED=0 @1
//--------------------------------------------------------------------------

class AMP_MG5_gg_uubarg
    : public gra::amplitude::MG5ProcessRegistry_gg_uubarg {
 public:


  using ColorFlowVector          = std::vector<std::complex<double>>;
  using ColorFlowHelicityMatrix  = std::vector<ColorFlowVector>;
// Constructor.
    AMP_MG5_gg_uubarg() { initProc(gra::aux::ResolveProjectPath("MG5cards/Durham/gg_uubarg/param_card.dat")); }
  // Keep event-local MG5 work arrays isolated between process instances
  AMP_MG5_gg_uubarg(const AMP_MG5_gg_uubarg &) = delete;
  AMP_MG5_gg_uubarg &operator=(const AMP_MG5_gg_uubarg &) = delete;

    // Initialize process.
    void initProc(std::string param_card_name);
    // Initialize every mass and coupling through the generated model
    void InitParameters(SLHAReader slha);
    // Compute the evaluated particle parameters of the generated model
    gra::mg5::ParticleMap Particles() const { return pars.Particles(); }
    // Compute the evaluated UFO electromagnetic coupling
    double AlphaQED() const { return pars.AlphaQED(); }

    // Calculate flavour-independent parts of cross section.
    gra::mg5helas::MatrixElementEvaluation Evaluate(gra::LORENTZSCALAR &lts, double alphas);
  void   CalcColorFlowHelicity(gra::LORENTZSCALAR &lts, double alphas,
                               ColorFlowHelicityMatrix &jamp_matrix);
  // Contract arbitrary color projectors while the generated amplitudes are live
  gra::mg5helas::EvaluationStatus CalcColorProjectedHelicity(
      gra::LORENTZSCALAR &lts, double alphas,
      const std::complex<double> *color_projectors, int projector_count,
      std::complex<double> *projected, gra::M4Vec *hard_k1 = nullptr,
      gra::M4Vec *hard_k2 = nullptr);

    // Evaluate sigmaHat(sHat).
    double sigmaHat();

    // Info on the subprocess.
    std::string name() const {return "g g > u u~ g (sm_lepton_masses)";}

    int code() const {return 1;}

    const std::vector<double> & getMasses() const {return mME;}

    // Get and set momenta for matrix element evaluation
    std::vector < double * > getMomenta(){return p;}
    void setMomenta(std::vector < double * > & momenta){p = momenta;}
    void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}

    // Get matrix element vector
    const double * getMatrixElements() const {return matrix_element;}

    // Constants for array limits
  static const int ninitial   = 2;
  static const int nexternal  = 5;
  static const int ncolor     = 6;
  static constexpr int nhelicity = 32;
  static const int nprocesses = 1;

  std::vector<double>              ColorDenominators() const;
  std::vector<std::vector<double>> ColorMetric() const;
 private:
  // Private functions to calculate the matrix element for all subprocesses
  // Prepare one physical on-shell HELAS phase-space point
  bool setup_kinematics(gra::LORENTZSCALAR &lts);
  void calculate_color_flows(std::complex<double> jamp[ncolor]) const;
// Calculate wavefunctions
  void calculate_wavefunctions(const int perm[], const int hel[]);
  static const int nwavefuncs = 18;
  std::complex<double> w[nwavefuncs][18];
  static const int namplitudes = 18;
  std::complex<double> amp[namplitudes];
  double matrix_1_gg_uuxg();

    // Store the matrix element value from sigmaKin
  double matrix_element[nprocesses];

    // Color flows, used when selecting color
  std::array<std::array<double, ncolor>, nprocesses> jamp2 = {};
// Pointer to the model parameters
  Parameters_sm_lepton_masses pars;  // GRANIITTI

    // vector with external particle masses
  std::vector<double> mME;

    // vector with momenta (to be changed each event)
  std::vector<double *> p;
  std::vector<std::array<double, 4>> momentum_buffer;
  std::vector<gra::M4Vec> final_buffer;
  // Initial particle ids
  int id1, id2;

};


#endif  // MG5_Sigma_sm_lepton_masses_gg_uuxg_H
