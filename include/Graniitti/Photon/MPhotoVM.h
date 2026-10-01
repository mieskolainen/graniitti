// Exclusive heavy vector-meson photoproduction amplitude
//
// [REFERENCE: S.P. Jones, A.D. Martin, M.G. Ryskin and T. Teubner, arXiv:1307.7099]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPHOTOVM_H
#define MPHOTOVM_H

// C++
#include <map>
#include <memory>
#include <string>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Photon/MPhotoQCD.h"
#include "Graniitti/Photon/MPhotoDiss.h"
#include "Graniitti/Process/MAmpMatch.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Spin/MDirac.h"

namespace gra {

// Available gluon inputs for heavy-vector photoproduction
enum class MPhotoVMGluonSource { LHAPDF_SHUVAEV, JMRT_2013_LO, JMRT_2013_NLO };

// Three-parameter integrated-gluon fit used by JMRT
struct MPhotoVMGluonFitParameters {
  double normalization = 0.0;
  double a             = 0.0;
  double b             = 0.0;
};

// Published constants and gluon fits from the 2013 JMRT analysis
struct MPhotoVMJMRT2013Parameters {
  MPhotoVMGluonFitParameters lo;
  MPhotoVMGluonFitParameters nlo;
  double                     lo_scale2      = 0.45;
  double                     nlo_q02        = 1.0;
  double                     nlo_lambda_qcd = 0.2;
  double                     alpha_s_mz     = 0.118;
  double                     mz             = 91.1876;
};

// Physics parameters for one heavy-vector photoproduction channel
struct MPhotoVMChannelParameters {
  int    quark_pdg                       = 0;
  int    ground_state_pdg                = 0;
  double wavefunction_correction         = 1.0;
  double leptonic_width_ee               = 0.0;
  double cross_section_slope_B0          = 0.0;
  double cross_section_slope_alpha_prime = 0.0;
  double cross_section_slope_W0          = 0.0;
};

// Numerical steering for heavy-vector photoproduction
struct MPhotoVMNumerics {
  double                                   jmrt_k2_max      = 125.0;
  double                                   jmrt_lambda_step = 0.05;
  unsigned int                             N_k              = 48;
  MPhotoVMGluonSource                      gluon_source     = MPhotoVMGluonSource::LHAPDF_SHUVAEV;
  MPhotoVMJMRT2013Parameters               jmrt_2013;
  std::map<int, double>                    running_alpha_lambda_qcd;
  std::map<int, double>                    heavy_quark_mass;
  std::map<int, MPhotoVMChannelParameters> channels;
  std::map<int, flux::PhotoDissParam> dissociation;

  bool initialized = false;

  // Compute the heavy-quark mass for a vector-meson quark PDG id
  double HeavyQuarkMass(int quark_pdg) const;

  // Compute the e+e- partial width used for vector-meson normalization
  double LeptonicWidthEE(int vm_pdg) const;

  // Compute the configured physics parameters for one vector meson
  const MPhotoVMChannelParameters &Channel(int vm_pdg) const;

  // Validate vector-meson numerical steering
  void Validate() const;

  // Configure vector meson physics and numerics from immutable JSON text
  void ConfigureFromJson(const std::string &general_file, const std::string &general_json,
                         const std::string &numerics_file, const std::string &numerics_json);

  // Configure vector meson physics and numerics from parsed card documents
  void Configure(const nlohmann::json &general, const std::string &general_file, const nlohmann::json &numerics,
                 const std::string &numerics_file);
};

using MPhotoVMNumericsPtr = std::shared_ptr<const MPhotoVMNumerics>;

// Construct one immutable heavy-vector photoproduction block
MPhotoVMNumericsPtr ReadPhotoVMNumerics(const MModelTune &tune);

// Compute the run owned heavy-vector photoproduction block
MPhotoVMNumericsPtr GetPhotoVMNumerics(MModelCache &cache);

// JMRT gamma-Pomeron vector-meson photoproduction amplitude
class MPhotoVM : public amplitude::ProcessFamily {
 public:
  // Construct the amplitude and bind its immutable channel process definition
  MPhotoVM(gra::LORENTZSCALAR &lts, MModelTunePtr model_tune, const std::string &channel,
           std::shared_ptr<const amplitude::ProcessDefinition> definition);

  // Build the immutable direct dilepton process for one vector-meson channel
  static std::shared_ptr<const amplitude::ProcessDefinition> ProcessDefinitionFor(const std::string &channel);

  // Initialize vector-meson and mode-dependent Sudakov parameters before copies
  static void InitializeParameters(MProcessSetup &setup);

  // Compute the direct dilepton final-state decay structure
  static constexpr DecayStructure DirectDecayStructure() { return {DecayType::Full}; }

  // Evaluate gamma p to heavy vector meson to dilepton kinematics
  double Amp2(gra::LORENTZSCALAR &lts) const;

  // Compute the elastic gamma-p differential cross section in nb/GeV^2
  double GammaPDSigmaDt(gra::LORENTZSCALAR &lts, double W, double t) const;

  // Compute the elastic gamma-p cross section in nb over zero to |t|max
  double GammaPCrossSection(gra::LORENTZSCALAR &lts, double W, double abs_t_max) const;

  // Compute the channel t slope in GeV^-2
  double GammaPTSlope(double W) const;

  // Clear shower color flow for color-singlet vector-meson decays
  void SampleColorFlow(gra::LORENTZSCALAR &lts) const;

 private:
  // Acquire shared read-only Sudakov and UGD tables
  void EnsureSudakov(gra::LORENTZSCALAR &lts) const;

  MModelTunePtr       model_tune;
  SoftModelPtr        soft_model;
  int                 pdg;
  MDirac              dirac;
  MPhotoVMNumericsPtr numerics;
  const flux::PhotoDissParam *diss = nullptr;
  std::complex<double> diss_forward = 0.0;
};

}  // namespace gra

#endif
