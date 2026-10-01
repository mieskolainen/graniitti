// Statistical nuclear particle and photon emission with exact recoil kinematics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEAREVAPORATION_H
#define MNUCLEAREVAPORATION_H

#include "Graniitti/Nuclear/MBreakup.h"
#include "Graniitti/Nuclear/MMass.h"
#include "HepMC3/GenEvent.h"
#include "HepMC3/GenParticle.h"

namespace gra::nuclear {

// Own completed nuclear masses and a spin averaged Weisskopf decay calculation
class MEvaporation {
 public:
  // Prepare nuclear masses and independent density and E1 inputs
  MEvaporation(const MassParam& mass, double radius, const GDRParam& gdr, unsigned int nodes, unsigned int steps);

  // Compute a bare ground state mass in GeV using the shared mass prescription
  double Mass(unsigned int a, unsigned int z) const;

  // Compute the Fermi momentum from the species density in a uniform sphere
  double Fermi(unsigned int a, unsigned int z, bool proton) const;

  // Append competing n, p, light ion and E1 gamma decays to a HepMC3 nucleus
  void Decay(HepMC3::GenEvent& event, HepMC3::GenParticlePtr parent, MRandom& random) const;

 private:
  MMass                                                   mass_;
  double                                                  radius_;
  GDRParam                                                gdr_;
  unsigned int                                            nodes_;
  unsigned int                                            steps_;

  // Append one isotropic two-body emission and compute the residual particle
  HepMC3::GenParticlePtr Emit(HepMC3::GenEvent& event, const HepMC3::GenParticlePtr& parent, int id1, double mass1,
                              int id2, double mass2, MRandom& random) const;

  // Compute the logarithmic continuum level density from Z(beta) = exp(a/beta)
  double LogDensity(unsigned int a, unsigned int z, double energy) const;
};

}  // namespace gra::nuclear

#endif
