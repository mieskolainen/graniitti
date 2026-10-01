// Leading twist nuclear gluon shadowing
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARSHADOW_H
#define MNUCLEARSHADOW_H

#include <memory>
#include <string>
#include <vector>

#include "Graniitti/PDF/MLHAPDF.h"

namespace gra::nuclear {

// Store published PDF inputs and numerical table controls
struct ShadowParam {
  std::string  pdf_set;
  int          pdf_member = 0;
  std::string  dpdf_set;
  int          dpdf_member     = 0;
  double       alpha0          = 0.0;
  double       alpha_prime     = 0.0;
  double       flux_b          = 0.0;
  double       dpdf_diss       = 0.0;  // (Elastic + low-mass dissociation) / elastic
  double       b_diff          = 0.0;
  double       x_min           = 0.0;  // Lower validated gluon x edge
  double       x_max           = 0.0;  // Upper leading twist x edge
  double       q2_min          = 0.0;  // Lower scale squared [GeV2]
  double       q2_max          = 0.0;  // Upper scale squared [GeV2]
  double       sigma3_anchor_x = 0.0;  // FGS10_H lower bound anchor
  double       sigma3_fade_x   = 0.0;  // Start of the large x linear fade
  double       sigma3_power    = 0.0;  // Intermediate x continuation power
  unsigned int x_nodes         = 0;
  unsigned int scale_nodes     = 0;
  unsigned int integral_nodes  = 0;
};

// Store the two and three nucleon effective cross sections
struct ShadowXS {
  double sigma2    = 0.0;  // Two nucleon effective cross section [mb]
  double sigma3    = 0.0;  // Higher rescattering effective cross section [mb]
  double sigma3_in = 0.0;  // Inelastic higher rescattering cross section [mb]
  double ratio     = 0.0;  // sigma2 / sigma3
};

// Tabulate the leading twist shadowing input
class MShadow {
 public:
  // Construct immutable PDF and interpolation tables
  explicit MShadow(ShadowParam param);

  // Compute the validated shadowing controls
  const ShadowParam &Param() const { return param_; }

  // Compute effective cross sections at one gluon x and scale
  ShadowXS CrossSections(double x, double scale2) const;

 private:
  ShadowParam                        param_;
  MLHAPDFStore                       pdf_store_;
  std::shared_ptr<const LHAPDF::PDF> pdf_;
  std::shared_ptr<const LHAPDF::PDF> dpdf_;
  std::vector<double>                log_x_;
  std::vector<double>                log_scale2_;
  std::vector<double>                sigma2_;

  // Validate controls and prepare the persistent interpolation table
  void Prepare();

  // Evaluate the two nucleon effective cross section directly
  double Sigma2(double x, double scale2) const;

  // Interpolate the two nucleon effective cross section
  double Interpolate(double x, double scale2) const;

  // Construct the complete table identity
  std::string CacheKey() const;

  // Load one compressed leading twist table
  bool ReadCache(const std::string &filename, const std::string &key);

  // Save one compressed leading twist table
  void WriteCache(const std::string &filename, const std::string &key) const;
};

}  // namespace gra::nuclear

#endif
