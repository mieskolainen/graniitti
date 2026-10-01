// Coherent and incoherent photonuclear target transitions
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARPHOTO_H
#define MNUCLEARPHOTO_H

#include <array>
#include <complex>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MFluctuation.h"
#include "Graniitti/Nuclear/MHotSpot.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Nuclear/MShadow.h"
#include "Graniitti/Nuclear/MTypes.h"

namespace gra::nuclear {

// Store the optical response interpolation and series controls
struct PhotoTableParam {
  double       qt_max         = 0.0;  // Transverse momentum table limit [GeV]
  double       qz_max         = 0.0;  // Longitudinal momentum table limit [GeV]
  unsigned int qt_nodes       = 0;    // Uniform transverse momentum nodes
  unsigned int qz_nodes       = 0;    // Uniform nonnegative longitudinal nodes
  unsigned int series_terms   = 0;    // Exponential moment terms
  double       series_abs_tol = 0.0;  // Maximum accepted Taylor tail
};

// Store photonuclear shadowing inputs and radial integration controls
struct PhotoParam {
  unsigned int     b_nodes     = 0;    // Radial quadrature nodes
  unsigned int     z_nodes     = 0;    // Longitudinal quadrature nodes
  double           phase_scale = 0.0;  // Gauss order relative to half the maximum Fourier phase
  PhotoTableParam  table;              // Optical response table controls
  FluctuationParam fluctuation;        // Numerical cross-section eigenstate rule
  HotSpotParam     hotspot;            // Subnucleon fluctuation controls
  ShadowParam      shadow;             // Leading twist gluon shadowing controls
};

// Select the relative neutron sign of one elementary production amplitude
enum class PhotoIsospin { Isoscalar, Isovector };

// Store one energy and channel dependent elementary photoproduction profile
struct PhotoProfile {
  double sigma_eff = 0.0;  // Effective produced-state cross section [mb]
  double slope     = 0.0;  // Elementary produced-state slope [GeV^-2]
  double eta       = 0.0;  // Real-to-imaginary amplitude ratio
  double omega     = 0.0;  // Relative cross-section variance
  double x         = 0.0;  // Hard gluon momentum fraction
  double scale2    = 0.0;  // Hard factorization scale squared [GeV^2]
  PhotoIsospin isospin = PhotoIsospin::Isoscalar;
};

// Select the photonuclear target model
enum class PhotoModel { Impulse, Glauber, LTA };

// Select the photon propagation direction in the nuclear rest frame
enum class PhotonDirection { PositiveZ, NegativeZ };

// Compute the incident photon direction for one target beam leg
PhotonDirection TargetPhotonDirection(int target_leg);

// Store fluctuations around the optical coherent mean and their Good-Walker moments
struct PhotoCurrent {
  CurrentStat                       stat;
  std::vector<std::complex<double>> sample;
};

// Store coherent and incoherent nuclear amplitude multipliers
struct PhotoTransition {
  std::complex<double>        coherent   = {0.0, 0.0};
  double                      incoherent = 0.0;
  std::optional<PhotoCurrent> current;
};

// Evaluate optical shadowing of one elementary photon-nucleon amplitude
class MPhoto {
 public:
  // Construct one immutable photonuclear target model
  MPhoto(const MNucleus& target, PhotoParam param, PhotoModel model = PhotoModel::Glauber);

  // Construct one target transition sharing an immutable nuclear model
  MPhoto(std::shared_ptr<const MNucleus> target, PhotoParam param, PhotoModel model = PhotoModel::Glauber);

  // Compute the validated photonuclear controls
  const PhotoParam& Param() const { return param_; }

  // Compute the selected photonuclear target model
  PhotoModel Model() const { return model_; }

  // Compute the optical coherent amplitude and sampled incoherent transitions
  // Stored currents retain centered fluctuations around the optical mean
  PhotoTransition Factors(const PhotoProfile& profile, double qx, double qy, double qz, PhotonDirection direction,
                          const MConfigBank* bank = nullptr, const MHotSpot* hotspot = nullptr,
                          bool with_current = false, CoherenceType sector = CoherenceType::Inclusive) const;

  // Compute one nucleon-center shadowed configuration current
  std::complex<double> ShadowCurrent(const PhotoProfile& profile, const MConfig& config, double qx, double qy,
                                     double qz, PhotonDirection direction) const;

  // Compute the mean-preserving hotspot currents of one fixed bank
  std::vector<std::complex<double>> ShadowCurrents(const PhotoProfile& profile, const MConfigBank& bank, double qx,
                                                   double qy, double qz, PhotonDirection direction,
                                                   const MHotSpot* hotspot = nullptr) const;

  // Compute the shadowed current moments of one fixed configuration bank
  CurrentStat ShadowStat(const PhotoProfile& profile, const MConfigBank& bank, double qx, double qy, double qz,
                         PhotonDirection direction, const MHotSpot* hotspot = nullptr) const;

 private:
  // Store one exact event-profile attenuation kernel
  struct AttenuationKernel {
    std::complex<double> strength = 0.0;
    std::vector<double>  inverse_width, production;
    PhotoIsospin         isospin   = PhotoIsospin::Isoscalar;
    double               direct    = 0.0;
    double               rescatter = 1.0;
  };

  std::shared_ptr<const MNucleus>   target_;
  PhotoParam                        param_;
  PhotoModel                        model_ = PhotoModel::Glauber;
  std::vector<double>               b_node_;
  std::vector<double>               b_weight_;
  std::vector<double>               z_node_;
  std::vector<double>               z_weight_;
  std::vector<double>               radial_weight_;
  std::vector<double>               density_weight_, charge_weight_;
  std::vector<double>               optical_depth_;
  std::vector<std::complex<double>> response_, charge_response_;
  std::vector<double>               longitudinal_, charge_longitudinal_;
  std::vector<double>               thickness_;
  std::vector<double>               sigma_scale_;
  std::unique_ptr<MShadow>          shadow_model_;
  double                            r_max_       = 0.0;
  double                            depth_max_   = 0.0;

  // Validate controls and construct the radial integration rule
  void Prepare();

  // Tabulate the common radial and matter-density integration weights
  void PrepareDensityKernel();

  // Integrate outgoing thickness directly on the optical quadrature
  void PrepareOpticalKernel();

  // Tabulate momentum transforms of the optical-depth moment basis
  void PrepareResponse();

  // Compute configured or event-specific cross-section eigenvalues
  const std::vector<double>& CrossSectionScales(const PhotoProfile& profile, std::vector<double>& scratch) const;

  // Compute optical current moments from one full momentum transfer
  CurrentStat OpticalStat(const PhotoProfile& profile, double qt, double qz, PhotonDirection direction) const;

  // Compute direct optical current moments from the original quadrature
  CurrentStat OpticalStatDirect(const PhotoProfile& profile, double qt, double qz) const;

  // Interpolate momentum transforms of the cached geometry moments
  std::vector<std::complex<double>> Moments(double qt, double qz, std::size_t terms,
                                            PhotoIsospin isospin = PhotoIsospin::Isoscalar) const;

  // Compute the production density relative to the matter density normalization
  double Source(double matter, double charge, PhotoIsospin isospin) const;

  // Interpolate the full spherical nuclear thickness
  double Thickness(double b) const;

  // Compute the long-coherence Glauber density variance
  double GlauberVariance(const PhotoProfile& profile, double qt, double qz) const;

  // Compute leading twist coherent and incoherent current moments
  CurrentStat LeadingStat(const PhotoProfile& profile, double qt, double qz) const;

  // Compute nucleon-center and hotspot shadowed currents in that order
  std::array<std::complex<double>, 2> ShadowPair(const MConfig& config, std::size_t sample, double qx, double qy,
                                                 double qz, PhotonDirection direction, const MHotSpot* hotspot,
                                                 const AttenuationKernel& kernel) const;

  // Compute one finite-range configuration attenuation eigen-amplitude
  std::complex<double> ConfigAttenuation(const MConfig& config, std::size_t produced, double inverse_width,
                                         PhotonDirection direction, const AttenuationKernel& kernel) const;

  // Factor the elementary profile into exact attenuation coefficients
  AttenuationKernel PrepareAttenuation(const PhotoProfile& profile) const;

  // Construct the exact nuclear-geometry cache key
  std::string CacheKey() const;

  // Load one exact kinematics-independent optical kernel
  bool ReadCache(const std::string& filename, const std::string& key);

  // Save one exact kinematics-independent optical kernel
  void WriteCache(const std::string& filename, const std::string& key) const;
};

}  // namespace gra::nuclear

#endif
