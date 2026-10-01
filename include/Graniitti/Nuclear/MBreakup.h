// Electromagnetic absorption and deposited excitation for nuclear decay
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARBREAKUP_H
#define MNUCLEARBREAKUP_H

#include <array>
#include <memory>
#include <string>
#include <vector>

#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra::nuclear {

// Select the electromagnetic form factor of the photon emitter
enum class EmitterType { Point, Proton, Nuclear };

// Store generic SLO giant dipole systematics
struct GDRSystematics {
  double energy_a    = 0.0;  // [MeV] A^(-1/3) energy coefficient
  double energy_b    = 0.0;  // [MeV] A^(-1/6) energy coefficient
  double width_norm  = 0.0;  // Width normalization
  double width_power = 0.0;  // Width energy power
  double trk         = 0.0;  // [mb MeV] TRK sum coefficient
  double strength    = 0.0;  // Sum-rule strength
};

// Store one isotope-specific giant dipole response
struct GDRIsotope {
  unsigned int a        = 0;    // Mass number
  unsigned int z        = 0;    // Proton number
  double       energy   = 0.0;  // [MeV] Pole energy
  double       width    = 0.0;  // [MeV] Width
  double       strength = 0.0;  // Sum-rule strength
};

// Store the complete giant dipole response
struct GDRParam {
  GDRSystematics          systematics;
  std::vector<GDRIsotope> isotope;
};

// Validate shared E1 inputs before absorption or statistical decay
void ValidateGDR(const GDRParam &param);

// Store one resolved E1 absorption line in GeV and mb
struct Dipole {
  double energy = 0.0, width = 0.0, peak = 0.0;

  // Construct an inactive line before model initialization
  Dipole() = default;

  // Resolve isotope measurements or the common mass systematics
  Dipole(unsigned int a, unsigned int z, const GDRParam &param);

  // Compute the standard Lorentzian photoabsorption cross section [mb]
  double Sigma(double omega) const;
};

// Store the deuteron photoabsorption response
struct DeuteronParam {
  double threshold = 0.0;  // [MeV] Deuteron breakup threshold
  double norm      = 0.0;  // [mb MeV^(3/2)] Cross-section normalization
};

// Store the Pauli-blocking response
struct PauliParam {
  double                low_edge  = 0.0;  // [MeV] Low-energy boundary
  double                high_edge = 0.0;  // [MeV] High-energy boundary
  double                low_exp   = 0.0;  // [MeV] Low-energy exponential coefficient
  std::array<double, 5> poly{};           // Polynomial coefficients
  double                high_exp = 0.0;   // [MeV] High-energy exponential coefficient
};

// Store the quasi-deuteron response
struct QDParam {
  DeuteronParam deuteron;
  double        levinger = 0.0;
  PauliParam    pauli;
};

// Store one isotope-specific nucleon resonance response
struct ResonanceIsotope {
  unsigned int a      = 0;    // Mass number
  unsigned int z      = 0;    // Proton number
  double       energy = 0.0;  // [MeV] Gaussian center
  double       width  = 0.0;  // [MeV] Gaussian width
  double       area   = 0.0;  // [mb MeV] Integrated strength
};

// Store the high-energy photoabsorption response
struct ResonanceParam {
  double                        threshold = 0.0;  // [MeV] Response threshold
  std::vector<ResonanceIsotope> isotope;
};

// Store the isotope-specific high-energy photonuclear continuum response
struct ContinuumParam {
  unsigned int a              = 0;    // Mass number of the fitted target
  unsigned int z              = 0;    // Proton number of the fitted target
  double       threshold      = 0.0;  // [MeV] Start of the smooth continuum onset
  double       match          = 0.0;  // [MeV] End of the smooth continuum onset
  double       omega0         = 0.0;  // [MeV] Logarithmic reference energy
  double       constant       = 0.0;  // [mb] Pomeron contribution
  double       log2           = 0.0;  // [mb] Squared-log coefficient
  double       norm           = 0.0;  // [mb/MeV] Continuum normalization
  double       mean           = 0.0;  // [MeV] Continuum energy offset
  double       width          = 0.0;  // [MeV] Continuum decay scale
};

// Store the total photoabsorption response
struct PhotoAbsParam {
  double         energy_min = 0.0;  // [MeV] Photoexcitation lower limit
  GDRParam       gdr;
  QDParam        qd;
  ResonanceParam resonance;
  ContinuumParam continuum;
};

// Store excitation-energy transfer controls
struct TransferParam {
  double E0    = 0.0;  // [MeV] Two-component thermalization scale
  double match = 0.0;  // [MeV] High-energy probability matching scale
};

// Store photoabsorption quadrature controls
struct ResponseNumerics {
  unsigned int nodes        = 0;    // Nodes per smooth response interval
  double       rel_tol      = 0.0;  // Fine versus coarse response tolerance
  double       tail_rel_tol = 0.0;  // Photon-energy tail convergence tolerance
};

// Store impact-profile interpolation controls
struct ProfileNumerics {
  std::size_t nodes = 0;    // Logarithmic impact-parameter nodes
  double      b_min = 0.0;  // [fm] Minimum impact parameter
};

// Store proton charge-transform controls
struct ProtonEmitterNumerics {
  double       q_max   = 0.0;  // [GeV] Form-factor transform limit
  unsigned int q_nodes = 0;    // Form-factor transform nodes
};

// Store the photonuclear electromagnetic-dissociation controls
struct BreakupParam {
  unsigned int          a       = 0;                   // Target mass number
  unsigned int          z       = 0;                   // Target proton number
  double                z_emit  = 0.0;                 // Absolute emitter charge in units of e
  double                gamma   = 0.0;                 // Emitter-target relative Lorentz factor
  EmitterType           emitter = EmitterType::Point;  // Emitter charge form factor
  PhotoAbsParam         photo;
  TransferParam         transfer;
  ResponseNumerics      response;
  ProfileNumerics       profile;
  ProtonEmitterNumerics proton;
};

// Store shared E1 inputs and complete responses for each target isotope
struct EMDParam {
  GDRParam                  gdr;
  std::vector<BreakupParam> isotope;
};

// Sample electromagnetic excitation independently on each ion at common impact parameter
class MBreakup {
 public:
  // Construct one immutable electromagnetic breakup model
  explicit MBreakup(BreakupParam param, std::shared_ptr<const MNucleus> emitter = nullptr,
                    form::ParamStore structure = {});

  // Compute the validated breakup controls
  const BreakupParam &Param() const { return param_; }

  // Compute whether electromagnetic breakup is active
  bool Enabled() const;

  // Compute the total photoabsorption cross section [mb]
  double PhotoAbsorption(double omega) const;

  // Compute the mean number of absorbed photons
  double Mean(double b) const;

  // Compute the derived upper impact-parameter support in fm
  double ImpactMax() const;

  // Compute the normalized discrete absorption spectrum on the energy quadrature
  std::vector<double> Spectrum(double b) const;

  // Compute the absorption quadrature energies in GeV
  const std::vector<double>& Energies() const { return omega_node_; }

  // Sample one absorbed photon from the impact-dependent spectrum
  double SampleEnergy(double b, MRandom &random) const;

  // Sample one physical deposited energy from the continuous TCM distribution
  double SampleTransfer(double omega, MRandom &random) const;

  // Sample excitation, optionally conditional on absorption with probability 1-exp(-Mean(b))
  double SampleExcitation(double b, MRandom &random, bool absorbed = false) const;

  // Compute the finite-size impact-space EPA density [GeV^-1 fm^-2]
  double PhotonDensity(double omega, double b) const;

 private:
  BreakupParam                    param_;
  std::shared_ptr<const MNucleus> emitter_;
  form::ParamStore                structure_;
  Dipole                          gdr_;
  double                          omega_min_    = 0.0;  // [GeV] Photoexcitation lower limit
  double                          omega_max_    = 0.0;  // [GeV] Photoabsorption upper limit
  double                          kernel_limit_ = 0.0;  // Dimensionless EPA support
  std::vector<double>             omega_edge_;          // [GeV] Response boundaries
  std::vector<double>             omega_node_;          // [GeV] Photon energies
  std::vector<double>             omega_weight_;        // [GeV] Quadrature weights
  std::vector<double>             absorption_;          // [mb] Photoabsorption values
  std::vector<double>             log_b_node_;          // Log impact parameter in fm
  std::vector<double>             mean_;                // Mean absorbed-photon multiplicity
  std::vector<double>             charge_b_node_;       // [fm] Charge-transform radii
  std::vector<double>             charge_fraction_;     // Enclosed transverse charge fraction

  // Compute the emitter charge fraction inside one transverse cylinder
  double ChargeFraction(double b) const;

  // Compute the stable dimensionless Bessel-function EPA kernel
  double BesselKernel(double x) const;

  // Fold the photon density and photoabsorption cross section directly
  double FoldMean(double b) const;

  // Tabulate the projected charge of the shared proton form factor
  void PrepareCharge();

  // Prepare and validate the piecewise photon energy quadrature
  void PrepareEnergy();

  // Prepare the validated logarithmic impact-parameter interpolation table
  void PrepareMean();
};

}  // namespace gra::nuclear

#endif
