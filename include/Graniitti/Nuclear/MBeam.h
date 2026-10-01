// Charged beam classification and process support
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARBEAM_H
#define MNUCLEARBEAM_H

#include <array>
#include <string>
#include <vector>

#include "Graniitti/Kinematics/M4Vec.h"

namespace gra::nuclear {

// Select one unordered physical beam pair
enum class CollisionType { Invalid, PP, EE, EP, PA, EA, AA };

// Compute whether one PDG code is a supported charged elementary lepton
bool IsChargedLepton(int pdg);

// Compute whether one PDG code is a proton or antiproton
bool IsProton(int pdg);

// Compute whether one PDG code is a charged nucleus or antinucleus
bool IsChargedNucleus(int pdg);

// Compute whether one PDG code can emit an EPA photon
bool IsEPAEmitter(int pdg);

// Compute one for an elementary beam or A for a nuclear beam
unsigned int BeamMassNumber(int pdg);

// Compute the full beam energy from the steering convention
double FullBeamEnergy(int pdg, double energy);

// Compute the steering energy from one full beam energy
double SteeringBeamEnergy(int pdg, double energy);

// Compute the charge encoded by one supported beam PDG code in units of e/3
int BeamChargeX3(int pdg);

// Validate one beam PDG code and its particle-table charge
void ValidateBeam(int pdg, int charge_x3);

// Classify one unordered pair of supported charged beams
CollisionType ClassifyCollision(int pdg1, int pdg2);

// Compute the compact process-table name of one collision class
std::string CollisionName(CollisionType type);

// Declare the collision classes implemented by a physical process
struct BeamSupport {
  std::vector<CollisionType> collisions = {CollisionType::PP};
  bool antiprotons = true;

  // Compute the charged collision classes for photon fusion
  static BeamSupport EPA();

  // Compute the charged collision classes for photoproduction
  static BeamSupport Photo(bool antiprotons = true);

  // Compute whether a collision class is supported
  bool Supports(CollisionType type) const;

  // Compute support for physical beam PDGs including charge conjugation
  bool Supports(int first, int second) const;

  // Compute the supported collision names for process help
  std::string Names() const;
};

// Compute the spatial momentum transfer in the incoming-particle rest frame
M3Vec RestTransfer(const M4Vec &beam, const M4Vec &transfer, double mass);

// Close nuclear transfers after subtracting high-energy beam momenta
bool CloseTransfers(const std::array<int, 2> &pdg, const std::array<M4Vec, 2> &beam,
                    const std::array<M4Vec, 2> &forward, const M4Vec &central, std::array<M4Vec, 2> &transfer);

}  // namespace gra::nuclear

#endif
