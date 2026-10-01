// Nuclear geometry and configuration value types
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARTYPES_H
#define MNUCLEARTYPES_H

#include <array>
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace gra::nuclear {

// Select a coherent, incoherent or inclusive nuclear sector
enum class CoherenceType { Coherent, Incoherent, Inclusive };

// Select the nucleon isospin state
enum class NucleonType { Proton, Neutron };

// Store one decoded 10LZZZAAAI nuclear PDG identity
struct NucleusID {
  int          pdg    = 0;
  unsigned int a      = 0;
  unsigned int z      = 0;
  unsigned int lambda = 0;
  unsigned int isomer = 0;
  bool         anti   = false;
};

// Store one spherical two-parameter Fermi density
struct DensityParam {
  double       radius        = 0.0;  // Half-density radius [fm]
  double       skin          = 0.0;  // Surface diffuseness [fm]
  double       r_max         = 0.0;  // Radial integration limit [fm]
  unsigned int nodes         = 0;    // Gauss-Legendre radial nodes
  std::size_t  cdf_nodes     = 0;    // Inverse-CDF intervals
  double       form_q_max    = 0.0;  // Cached form-factor momentum limit [GeV]
  std::size_t  form_nodes    = 0;    // Cached form-factor grid nodes
  double       form_abs_tol  = 0.0;  // Form-factor interpolation tolerance
};

// Store one nucleus and its separate charge and matter densities
struct NucleusParam {
  int          pdg  = 0;    // Nuclear PDG code
  double       mass = 0.0;  // Bare nuclear mass [GeV]
  DensityParam charge;      // Charge density controls
  DensityParam matter;      // Matter density controls
};

// Store one explicit charge and matter Fermi-density override
struct NucleusShape {
  unsigned int          a      = 0;           // Nuclear mass number
  unsigned int          z      = 0;           // Nuclear proton number
  std::array<double, 2> charge = {0.0, 0.0};  // Charge radius and skin [fm]
  std::array<double, 2> matter = {0.0, 0.0};  // Matter radius and skin [fm]
};

// Store the generic density rule and explicit nuclear overrides
struct GeometryParam {
  double                    radius_scale  = 0.0;  // Generic A^(1/3) radius coefficient [fm]
  double                    skin          = 0.0;  // Generic Fermi surface diffuseness [fm]
  unsigned int              nodes         = 0;    // Density Gauss-Legendre radial nodes
  double                    tail_skin     = 0.0;  // Retained density tail in skin lengths
  std::size_t               cdf_nodes     = 0;    // Density sampling intervals
  double                    form_q_max    = 0.0;  // Cached form-factor momentum limit [GeV]
  std::size_t               form_nodes    = 0;    // Cached form-factor grid nodes
  double                    form_abs_tol  = 0.0;  // Form-factor interpolation tolerance
  std::vector<NucleusShape> nucleus;              // Explicit isotope density inputs
};

// Store correlated configuration sampling controls
struct ConfigParam {
  std::size_t count             = 0;    // Configurations per nuclear leg
  double      d_min             = 0.0;  // Minimum nucleon-center distance [fm]
  std::size_t max_trials        = 0;    // Placement trials per nucleon
  std::size_t sweeps            = 0;    // Hard-core equilibration sweeps
  std::size_t neutron_nodes     = 0;    // Neutron-density intervals
  double      density_rel_tol   = 0.0;  // Relative negative-density tolerance
  double      negative_norm_tol = 0.0;  // Integrated negative-density tolerance
  double      norm_tol          = 0.0;  // Density normalization tolerance
};

// Store one proton or neutron center in the nuclear rest frame
struct Nucleon {
  std::array<double, 3> x    = {0.0, 0.0, 0.0};  // Position [fm]
  NucleonType           type = NucleonType::Neutron;
};

// Store the first two ensemble moments of one configuration current
struct CurrentStat {
  std::complex<double> mean     = {0.0, 0.0};
  double               second   = 0.0;
  double               variance = 0.0;
};

}  // namespace gra::nuclear

#endif
