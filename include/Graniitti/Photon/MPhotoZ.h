// Exclusive Z photoproduction amplitude
//
// [REFERENCE: A. Cisek, W. Schafer and A. Szczurek, PRD 80 (2009) 074013, arXiv:0906.1739]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPHOTOZ_H
#define MPHOTOZ_H

// C++
#include <memory>
#include <string>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Photon/MPhotoQCD.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Spin/MDirac.h"

namespace gra {

// Dipole and kT-factorized gamma-Pomeron Z photoproduction amplitude
// Direct neutral-current photoproduction in the transverse high-energy SCHC
// approximation
class MPhotoZ : public amplitude::ProcessFamily {
 public:
  // Construct the photoproduction amplitude with one immutable process
  // definition
  MPhotoZ(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune,
          std::shared_ptr<const amplitude::ProcessDefinition> definition);

  // Build the immutable direct neutral-current process definition
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor();

  // Initialize photoproduction and Sudakov parameters before worker copies
  static void InitializeParameters(MProcessSetup &setup);

  // Compute the direct neutral-current final-state decay structure
  static constexpr DecayStructure DirectDecayStructure() { return {DecayType::Full}; }

  // Evaluate gamma p to neutral-current f fbar kinematics
  double Amp2(gra::LORENTZSCALAR &lts) const;

  // Assign shower-compatible color flow for quark-pair final states
  void SampleColorFlow(gra::LORENTZSCALAR &lts) const;

  // Compute the on-shell gamma p to Z p total cross section in nanobarns
  double GammaPTotalCrossSectionNb(gra::LORENTZSCALAR &lts, double W) const;

 private:
  // Acquire shared read-only Sudakov and UGD tables
  void EnsureSudakov(gra::LORENTZSCALAR &lts) const;

  MModelTunePtr        model_tune;
  SoftModelPtr         soft_model;
  MDirac               dirac;
  MPhotoQCDParamPtr    param;
  MPhotoQCDNumericsPtr numerics;
};

}  // namespace gra

#endif
