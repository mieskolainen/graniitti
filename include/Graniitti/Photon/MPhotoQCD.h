// QCD photoproduction utilities
//
// [REFERENCE: A. Cisek, W. Schafer and A. Szczurek, PRD 80 (2009) 074013, arXiv:0906.1739]
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPHOTOQCD_H
#define MPHOTOQCD_H

// C++
#include <array>
#include <complex>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Regge/MSoftModel.h"

namespace gra {

// Physical steering for CSS-like photoproduction impact factors
struct MPhotoQCDParam {
  double xg_coefficient = 1.0;
  double mu_MIN = 1.5;
  double mu_over_m = 0.5;
  double light_quark_mass_MIN = 1e-3;
  double real_part_delta = 0.20;
  double t_slope_B0 = 3.5;
  double t_slope_alpha_prime = 0.164;
  double t_slope_W0 = 95.0;

  bool initialized = false;

  // Read photoproduction physics steering from one immutable GENERAL JSON text
  void ConfigureFromJson(const std::string &source_file,
                         const std::string &json_text,
                         const std::string &block_name);

  // Read photoproduction physics steering from one parsed GENERAL document
  void Configure(const nlohmann::json &document,
                 const std::string &source_file,
                 const std::string &block_name);

  // Validate photoproduction physics steering
  void Validate(const std::string &block_name) const;
};

// Numerical steering for CSS-like photoproduction impact-factor integrals
struct MPhotoQCDNumerics {
  double kappa2_MIN = 2.0;
  double kappa2_MAX = 100.0;
  double k2_MIN = 1e-4;
  double k2_MAX = 100.0;
  bool k2_MAX_use_q2 = true;

  unsigned int N_kappa = 24;
  unsigned int N_k = 24;

  bool initialized = false;

  std::vector<double> kappa2_nodes;
  std::vector<double> kappa2_weights;

  // Compute the active k2 upper integration bound
  double K2Max(double q2) const;

  // Prepare fixed integration rules for repeated impact-factor calls
  void PrepareIntegrationRules();

  // Read numerical steering from one immutable NUMERICS JSON text
  void ConfigureFromJson(const std::string &source_file,
                         const std::string &json_text,
                         const std::string &block_name);

  // Read numerical steering from one parsed NUMERICS document
  void Configure(const nlohmann::json &document,
                 const std::string &source_file,
                 const std::string &block_name);

  // Validate photoproduction numerical steering against the Sudakov grid
  void Validate(double sudakov_q2_max, const std::string &block_name) const;
};

using MPhotoQCDParamPtr = std::shared_ptr<const MPhotoQCDParam>;
using MPhotoQCDNumericsPtr = std::shared_ptr<const MPhotoQCDNumerics>;

// Construct one immutable photoproduction physics block
MPhotoQCDParamPtr ReadPhotoQCDParam(const MModelTune &tune,
                                    const std::string &block_name);

// Construct one immutable photoproduction numerical block
MPhotoQCDNumericsPtr ReadPhotoQCDNumerics(const MModelTune &tune,
                                          const std::string &block_name);

// Compute one run owned photoproduction physics block
MPhotoQCDParamPtr GetPhotoQCDParam(MModelCache &cache,
                                   const std::string &block_name);

// Compute one run owned photoproduction numerical block
MPhotoQCDNumericsPtr GetPhotoQCDNumerics(MModelCache &cache,
                                         const std::string &block_name);

// Shared functions for gamma-Pomeron QCD photoproduction amplitudes
class MPhotoQCD {
public:
  // Acquire shared read-only Sudakov and UGD tables
  static void EnsureSudakov(gra::LORENTZSCALAR &lts,
                            const SoftModelPtr &soft_model,
                            const std::string &context);

  // Compute a positive invariant mass squared for the generated central system
  static double CentralMass2(const gra::LORENTZSCALAR &lts);

  // Compute a steerable real-part correction for high-energy vector production
  static std::complex<double> RealPartFactor(const MPhotoQCDParam &param);

  // Compute the elastic t-slope factor used around the UGD impact factor
  static double TSlopeFactor(const MPhotoQCDParam &param, double w2, double t);

  // Compute the elastic or dissociative target transition factor
  static double TargetTransitionFactor(const ForwardLegState &state,
                                       const MPhotoQCDParam &param, double w2,
                                       const SoftModel &soft_model);

  // Compute the cross-section t slope associated with the amplitude factor
  static double TSlope(const MPhotoQCDParam &param, double w2);

  // Expand spin-averaged source amplitudes over the four incoming proton spin
  // rows
  static std::vector<std::complex<double>>
  InitialProtonSpinCopies(const std::vector<std::complex<double>> &amplitudes);

  // Compute the small-x value used by the impact-factor model
  static double XGluon(const MPhotoQCDParam &param, double q2, double w2);

  // Compute the UGD hard scale used by the impact-factor model
  static double HardScale(const MPhotoQCDParam &param, double q2);

  // Contract a lower-index current with an upper-index polarization vector
  static std::complex<double> ContractCurrentPolarization(
      const MDirac::Current &current,
      const FTensor::Tensor1<std::complex<double>, 4> &eps);

  // Compute transported transverse spin states dual to exp(i m phi) photon sources
  static std::array<FTensor::Tensor1<std::complex<double>, 4>, 2>
  TransverseStates(const LORENTZSCALAR &lts, const M4Vec &boson);

  // Compute one flavour contribution to the CSS z and k2 integral
  static std::complex<double>
  CSSFlavorIntegral(const gra::LORENTZSCALAR &lts, const MPhotoQCDParam &param,
                    const MPhotoQCDNumerics &num, double xg,
                    double quark_mass, double q2, double mu);
};

} // namespace gra

#endif
