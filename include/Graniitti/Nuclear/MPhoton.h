// Coherent and incoherent nuclear transverse photon densities
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARPHOTON_H
#define MNUCLEARPHOTON_H

#include <complex>
#include <memory>

#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Nuclear/MConfig.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Particle/MForm.h"

namespace gra::nuclear {

// Store the two eigenvalues of one transverse photon density
struct PhotonDensity {
  double parallel      = 0.0;
  double perpendicular = 0.0;

  // Compute the unpolarized photon density
  double Trace() const;
};

// Compute one transverse nuclear photon-source helicity amplitude
std::complex<double> PhotonAmp(double source, double xi, double charge, double phi, int helicity);

// Compute the dimensionless point charge impact EPA kernel
double PhotonKernel(double x, double gamma);

// Compute the impact EPA density [GeV^-1 fm^-2] with an enclosed charge fraction
double ImpactPhotonDensity(double omega, double b, double gamma, double charge, double fraction);

// Evaluate coherent and Good-Walker incoherent nuclear photon densities
class MPhoton {
 public:
  // Construct one immutable nuclear photon source
  explicit MPhoton(const MNucleus &nucleus, form::ParamStore structure = {});

  // Construct one photon source sharing an immutable nuclear model
  explicit MPhoton(std::shared_ptr<const MNucleus> nucleus, form::ParamStore structure = {});

  // Compute the nuclear source model
  const MNucleus &Nucleus() const { return *nucleus_; }

  // Compute one coherent, incoherent or inclusive nuclear photon density
  PhotonDensity Density(CoherenceType type, double xi, double t, double pt) const;

  // Compute a configuration density using the full nuclear-rest-frame transfer
  PhotonDensity Density(CoherenceType type, double xi, double t, double pt, const M3Vec& q,
                        const MConfigBank* bank) const;

  // Compute the signed coherent transverse photon source
  double CoherentAmp(double xi, double t, double pt) const;

 private:
  std::shared_ptr<const MNucleus> nucleus_;
  form::ParamStore                structure_;

  // Compute the elastic bound-proton transverse density
  PhotonDensity Proton(double xi, double t, double pt) const;

  // Compute the charge-current variance in proton-count units
  double Variance(const M3Vec& q, const MConfigBank* bank) const;

  // Compute the incoherent bound-proton fluctuation density
  PhotonDensity Incoherent(double xi, double t, double pt, const M3Vec& q, const MConfigBank* bank) const;
};

}  // namespace gra::nuclear

#endif
