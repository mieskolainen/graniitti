// Gamma amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MGAMMA_H
#define MGAMMA_H

// C++
#include <complex>
#include <memory>
#include <map>
#include <mutex>
#include <random>
#include <string>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MSoftModel.h"

namespace gra {

// Select the gamma-gamma amplitude family exposed through the process
// definition
enum class MGammaMode { Generic, Higgs, Monopolium, FermionPair, Flux };

struct MMonopoleParam {
  double      monopole_mass    = 0.0;
  double      monopolium_mass  = 0.0;
  double      monopolium_width = 0.0;
  std::string coupling         = "null";
  std::string wavefunction     = "null";
  int         gn               = 0;
  bool        initialized      = false;

  // Validate parameters shared by monopole and monopolium production
  void Validate() const;

  // Validate the physical monopolium bound-state pole
  void ValidateMonopolium() const;

  // Compute the monopolium binding energy relative to two free monopoles
  double BindingEnergy() const;

  // Compute the Coulombic monopolium wavefunction at the origin
  double PsiAtOrigin() const;

  // Compute the monopolium diphoton width for a magnetic coupling
  double GammaGamma(double alpha_g) const;
};

using MMonopoleParamPtr = std::shared_ptr<const MMonopoleParam>;

// Construct one immutable monopole parameter block
MMonopoleParamPtr ReadMonopoleParam(const MModelTune &tune, double monopole_mass, double monopolium_mass,
                                    double monopolium_width);

// Compute one run owned immutable monopole parameter block
MMonopoleParamPtr GetMonopoleParam(MModelCache &cache, double monopole_mass, double monopolium_mass,
                                   double monopolium_width);

class MGamma : public amplitude::ProcessFamily {
 public:
  // Construct the gamma-gamma amplitude with one immutable process definition
  MGamma(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune,
         std::shared_ptr<const amplitude::ProcessDefinition> definition);

  // Destroy the gamma-gamma amplitude
  ~MGamma() {}

  // Build one immutable gamma-gamma process definition for an analytic mode
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor(MGammaMode mode);

  // Initialize mode-dependent monopole parameters before worker process copies
  static void InitializeParameters(const MProcessSetup &setup, MGammaMode mode);

  // Evaluate yy->lepton pair, quark pair, or monopole antimonopole
  mg5helas::MatrixElementEvaluation yyffbar(gra::LORENTZSCALAR &lts, bool coherent_epa);

  // Assign the unique color-singlet flow of a direct quark pair
  bool SampleColorFlow(gra::LORENTZSCALAR &lts);

  // yy->SM Higgs
  double yyHiggs(gra::LORENTZSCALAR &lts);

  // yy->monopolium
  double yyMP(gra::LORENTZSCALAR &lts);

 protected:
  // Compute the physical coupling multiplier of one direct fermion pair
  bool PairScale(const gra::LORENTZSCALAR &lts, double &scale) const;

  // Build and contract one complex CP-even scalar diphoton hard tensor
  double ScalarYY(gra::LORENTZSCALAR &lts, double mass, double gamma_yy, double gamma_tot);

  // Compute immutable parameters required by a monopole amplitude
  const MMonopoleParam &MonopoleParameters(const gra::LORENTZSCALAR &lts) const;

  // Print monopole parameters once for this process helper instance
  void PrintMonopoleParameters(const MMonopoleParam &parameters, double sqrts) const;

  // Higgs pole data copied from the active process PDG and decay tables
  MModelTunePtr                     model_tune;
  MParticle                         higgs;
  double                            higgs_gamma_gamma_width = 0.0;
  mutable MMonopoleParamPtr         monopole_param_handle;
  mutable std::once_flag            monopole_parameter_load_once;
  mutable std::once_flag            monopole_parameters_print_once;
  // Own a separately initialized massive Dirac model for each allowed species
  std::map<int, std::unique_ptr<PhotonMG5Process>> pair_amplitudes;
  PhotonMG5Process *pair_amplitude = nullptr;
};

}  // namespace gra

#endif
