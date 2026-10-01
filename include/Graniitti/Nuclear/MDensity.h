// Spherical nuclear charge and matter densities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MDENSITY_H
#define MDENSITY_H

#include <vector>

#include "Graniitti/Nuclear/MTypes.h"

namespace gra::nuclear {

// Own one normalized spherical two-parameter Fermi density
class MDensity {
 public:
  // Construct and normalize one density on a fixed radial quadrature
  explicit MDensity(DensityParam param);

  // Compute the validated density parameters
  const DensityParam &Param() const { return param_; }

  // Compute the normalized three-dimensional density in fm^-3
  double Rho(double r) const;

  // Compute the normalized spherical form factor for momentum in GeV
  double Form(double q) const;

  // Compute the normalized transverse thickness in fm^-2
  double Thick(double b) const;

  // Compute the charge fraction inside one transverse cylinder
  double Cylinder(double b) const;

  // Compute the root-mean-square radius in fm
  double Rms() const { return rms_; }

  // Map one uniform variate to the normalized radial density
  double Radius(double u) const;

 private:
  DensityParam        param_;
  double              norm_ = 0.0;
  double              rms_  = 0.0;
  std::vector<double> r_node_;
  std::vector<double> r_weight_;
  std::vector<double> cdf_node_;
  std::vector<double> cdf_;
  std::vector<double> form_r_node_;
  std::vector<double> form_r_weight_;
  std::vector<double> form_;
  double              form_step_ = 0.0;

  // Compute the unnormalized two-parameter Fermi profile
  double Profile(double r) const;

  // Validate parameters and select the finite tail limit
  void Validate();

  // Prepare normalization, moments and inverse radial sampling data
  void Prepare();

  // Compute the direct radial-quadrature form factor
  double FormDirect(double q) const;

  // Evaluate the form factor with one explicit radial rule
  double FormRule(double q, const std::vector<double> &node, const std::vector<double> &weight) const;

  // Compute a radial node count resolving one oscillatory form factor
  std::size_t FormNodeCount(double q, double scale) const;

  // Interpolate the cached form factor on its uniform momentum grid
  double FormCached(double q) const;

  // Prepare and validate the immutable form-factor table
  void PrepareForm();
};

}  // namespace gra::nuclear

#endif
