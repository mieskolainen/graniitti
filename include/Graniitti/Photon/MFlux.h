// Photon and other fluxes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MFLUX_H
#define MFLUX_H

// C++
#include <complex>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

// Own
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {
namespace flux {

// Store the two eigenvalues of one transverse EPA photon density
struct TransversePhotonFlux {
  double parallel      = 0.0;
  double perpendicular = 0.0;

  // Compute the scalar unpolarized flux as the density-matrix trace
  double Trace() const { return parallel + perpendicular; }
};

// Store ordered amplitude-level weights for unresolved nuclear EPA sectors
struct EPASectorWeights {
  std::vector<double> amplitude_fraction = {1.0};
};

// Store one applied EPA normalization and its resolved source weights
struct EPAWeight {
  double           amp2            = 0.0;
  double           amplitude_scale = 0.0;
  EPASectorWeights sector;
};

// Store real coefficients of a transverse photon density eigenbasis
struct TransversePhotonSource {
  double parallel      = 0.0;
  double perpendicular = 0.0;
};

// Store one resolved forward photon source, density and process-sector code
struct PhotonFluxSector {
  TransversePhotonFlux   density;
  TransversePhotonSource source;
  std::uint8_t           code = 0;
};

// Store elementary profile inputs needed by a resolved photon target
struct PhotoTargetProfile {
  double sigma_eff = 0.0;
  double slope     = 0.0;
  double eta       = 0.0;
  double x         = 0.0;
  double scale2    = 0.0;
  nuclear::PhotoIsospin isospin = nuclear::PhotoIsospin::Isoscalar;
};

// Store ordered target-transition factors for one photon direction
struct PhotoTargetFactors {
  std::vector<std::complex<double>>    factor;
  std::optional<nuclear::PhotoCurrent> current;
};

// Construct zero target amplitudes with the configured sectors and current bank
PhotoTargetFactors ZeroPhotoTarget(const ForwardLegState &state);

// Compute ordered forward photon densities without exposing their source model
std::vector<PhotonFluxSector> ForwardPhotonFluxSectors(const ForwardLegState &state);

// Validate one configured photon emitter through the common source layer
void ValidatePhotonEmitter(const LORENTZSCALAR &lts, int leg, const std::string &context);

// Compute whether one forward state supports the elementary photon target
bool SupportsPhotoTarget(const ForwardLegState &state);

// Compute whether one ordered photon beam direction has a physical target
bool SupportsPhotoDirection(const LORENTZSCALAR &lts, int photon_leg);

// Compute the elementary proton state for a nuclear photoproduction target
ForwardLegState ElementaryPhotoTarget(const LORENTZSCALAR &lts, const ForwardLegState &state);

// Sum physical photon directions and preserve each nuclear target current for screening
void SumPhotoTerms(LORENTZSCALAR &lts, std::vector<nuclear::PhotoTerm> terms);

// Compute the hadron or constituent-nucleon momentum seen by the hard amplitude
M4Vec PhotoTargetMomentum(const ForwardLegState &state);

// Resolve hadron and nuclear target factors outside the hard amplitude module
PhotoTargetFactors ResolvePhotoTarget(const ForwardLegState &state, const PhotoTargetProfile &profile,
                                      std::complex<double> elastic_factor, std::complex<double> hadron_factor);

// Complete the elementary photoproduction amplitudes over initial spin states
std::vector<std::complex<double>> CompletePhotoInitialSpinStates(const LORENTZSCALAR                     &lts,
                                                                 const std::vector<std::complex<double>> &amplitudes);

// Compute the initial-state spin average for one photoproduction process
double PhotoInitialSpinAverage(const LORENTZSCALAR &lts);

// Compute the number of resolved photon-source sectors for one forward leg
std::size_t PhotoSourceCount(const ForwardLegState &state);

// Compute all resolved photon-source amplitudes for one helicity
std::vector<std::complex<double>> PhotoSourceAmplitudes(const LORENTZSCALAR &lts, const ForwardLegState &state,
                                                        int helicity);

// Compute resolved scalar photon-source factors without helicity phases
std::vector<std::complex<double>> PhotoSourceScalarAmplitudes(const LORENTZSCALAR &lts, const ForwardLegState &state);

// Compute whether two photon directions need distinct amplitude channels
bool SplitPhotoDirections(const LORENTZSCALAR &lts, const ForwardLegState &upper, const ForwardLegState &lower);

// Store the canonical direction and nuclear-sector order for screening
void ConfigurePhotoLayout(LORENTZSCALAR &lts);

// Apply the common EPA normalization and ordered sector weights to amplitudes
bool ApplyEPAAmplitudeWeights(std::vector<std::complex<double>> &amplitudes, double common,
                              const EPASectorWeights &weights);

// Compute the elastic transverse EPA density of a charged spin-half emitter
TransversePhotonFlux ElasticSpinHalfFluxTransverse(double x, double t, double pt, double mass,
                                                   double electric_form_factor, double magnetic_form_factor,
                                                   double charge);

// Compute the elastic transverse EPA density of a proton
TransversePhotonFlux CohFluxTransverse(double x, double t, double pt, const form::ParamStore &structure = {});

// Compute the inelastic transverse EPA density of a dissociated proton
TransversePhotonFlux IncohFluxTransverse(double x, double t, double pt, double M2,
                                         const form::ParamStore &structure = {});

// Compute the scalar coherent EPA flux
double CohFlux(double x, double t, double pt, const form::ParamStore &structure = {});

// Compute the scalar inelastic EPA flux
double IncohFlux(double x, double t, double pt, double M2, const form::ParamStore &structure = {});

// Compute the t-integrated Drees-Zeppenfeld coherent photon flux
double DZFlux(double x);

// Compute the coherent or inclusive photon flux for one generated forward leg
double ForwardPhotonFlux(const ForwardLegState &state);

// Compute the transverse coherent or inclusive photon density of one leg
TransversePhotonFlux ForwardPhotonFluxTransverse(const ForwardLegState &state);

// Compute one selected nuclear photon-emission density of a nuclear leg
TransversePhotonFlux NuclearPhotonFluxTransverse(const ForwardLegState &state, nuclear::CoherenceType emission);

// Compute the phase-space conversion for exact non-collinear EPA kinematics
double ExactktEPAPhaseSpaceFactor(const LORENTZSCALAR &lts);

// Compute the incoming pp to collinear gamma-gamma flux conversion
double CollinearPhotonPhaseSpaceFactor(const LORENTZSCALAR &lts);

// Set the photon hard scales and alpha_s from the configured PDF
bool SetPhotonAlphaS(LORENTZSCALAR &lts);

// Apply exact non-collinear EPA fluxes at cross-section level
EPAWeight ApplyktEPAfluxes(double amp2, LORENTZSCALAR &lts);

// Apply transverse EPA currents with one fixed outer phase-space conversion
EPAWeight ApplyktEPAcurrents(double amp2, LORENTZSCALAR &lts, double phase_space);

// Apply collinear Drees-Zeppenfeld photon fluxes at cross-section level
double ApplyDZfluxes(double amp2, LORENTZSCALAR &lts);

// Apply collinear LUX photon densities at cross-section level
double ApplyLUXfluxes(double amp2, LORENTZSCALAR &lts);

}  // namespace flux
}  // namespace gra

#endif
