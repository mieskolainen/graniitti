// Charged beam classification and process support
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MBeam.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Nuclear/MNucleus.h"
#include "Graniitti/Particle/MPDG.h"

namespace gra::nuclear {
namespace {

// Compute the largest absolute four-vector component
double MaxComponent(const M4Vec &p) {
  return std::max({std::abs(p.E()), std::abs(p.Px()), std::abs(p.Py()), std::abs(p.Pz())});
}

// Compute true when every four-vector component is finite
bool Finite(const M4Vec &p) {
  return std::isfinite(p.E()) && std::isfinite(p.Px()) && std::isfinite(p.Py()) && std::isfinite(p.Pz());
}

}  // namespace

// Compute whether one PDG code is a supported charged elementary lepton
bool IsChargedLepton(const int pdg) {
  const int apdg = std::abs(pdg);
  return apdg == 11 || apdg == 13 || apdg == 15;
}

// Compute whether one PDG code is a proton or antiproton
bool IsProton(const int pdg) { return std::abs(pdg) == PDG::PDG_p; }

// Compute whether one PDG code is a charged nucleus or antinucleus
bool IsChargedNucleus(const int pdg) { return IsNuclearPDG(pdg) && DecodeNuclearPDG(pdg).z > 0; }

// Compute whether one PDG code can emit an EPA photon
bool IsEPAEmitter(const int pdg) { return IsChargedLepton(pdg) || IsProton(pdg) || IsChargedNucleus(pdg); }

// Compute one for an elementary beam or A for a nuclear beam
unsigned int BeamMassNumber(const int pdg) { return IsNuclearPDG(pdg) ? DecodeNuclearPDG(pdg).a : 1U; }

// Compute the full beam energy from the steering convention
// E_full = A E_nucleon
double FullBeamEnergy(const int pdg, const double energy) { return energy * static_cast<double>(BeamMassNumber(pdg)); }

// Compute the steering energy from one full beam energy
// E_nucleon = E_full / A
double SteeringBeamEnergy(const int pdg, const double energy) {
  return energy / static_cast<double>(BeamMassNumber(pdg));
}

// Compute the charge encoded by one supported beam PDG code in units of e/3
int BeamChargeX3(const int pdg) {
  if (IsChargedLepton(pdg)) { return pdg > 0 ? -3 : 3; }
  if (IsProton(pdg)) { return pdg > 0 ? 3 : -3; }
  if (IsChargedNucleus(pdg)) {
    const int charge = 3 * static_cast<int>(DecodeNuclearPDG(pdg).z);
    return pdg > 0 ? charge : -charge;
  }
  throw std::invalid_argument("BeamChargeX3: unsupported or neutral beam PDG code");
}

// Validate one beam PDG code and its particle-table charge
void ValidateBeam(const int pdg, const int charge_x3) {
  if (!IsEPAEmitter(pdg)) { throw std::invalid_argument("ValidateBeam: unsupported or neutral EPA beam"); }
  if (charge_x3 != BeamChargeX3(pdg)) {
    throw std::invalid_argument("ValidateBeam: particle-table beam charge disagrees with PDG code");
  }
}

// Classify one unordered pair of supported charged beams
CollisionType ClassifyCollision(const int pdg1, const int pdg2) {
  if (!IsEPAEmitter(pdg1) || !IsEPAEmitter(pdg2)) { return CollisionType::Invalid; }
  const bool lepton1  = IsChargedLepton(pdg1);
  const bool lepton2  = IsChargedLepton(pdg2);
  const bool proton1  = IsProton(pdg1);
  const bool proton2  = IsProton(pdg2);
  const bool nucleus1 = IsChargedNucleus(pdg1);
  const bool nucleus2 = IsChargedNucleus(pdg2);

  if (proton1 && proton2) { return CollisionType::PP; }
  if (lepton1 && lepton2) { return CollisionType::EE; }
  if ((lepton1 && proton2) || (proton1 && lepton2)) { return CollisionType::EP; }
  if ((proton1 && nucleus2) || (nucleus1 && proton2)) { return CollisionType::PA; }
  if ((lepton1 && nucleus2) || (nucleus1 && lepton2)) { return CollisionType::EA; }
  if (nucleus1 && nucleus2) { return CollisionType::AA; }
  return CollisionType::Invalid;
}

// Compute the compact process-table name of one collision class
std::string CollisionName(const CollisionType type) {
  switch (type) {
    case CollisionType::PP:
      return "pp";
    case CollisionType::EE:
      return "ee";
    case CollisionType::EP:
      return "ep";
    case CollisionType::PA:
      return "pA";
    case CollisionType::EA:
      return "eA";
    case CollisionType::AA:
      return "AA";
    case CollisionType::Invalid:
      return "invalid";
  }
  throw std::invalid_argument("CollisionName: unknown collision type");
}

// Compute the charged collision classes for photon fusion
BeamSupport BeamSupport::EPA() {
  using enum CollisionType;
  return {{PP, EE, EP, PA, EA, AA}};
}

// Compute the charged collision classes for photoproduction
BeamSupport BeamSupport::Photo(const bool antiprotons) {
  using enum CollisionType;
  return {{PP, EP, PA, EA, AA}, antiprotons};
}

// Compute whether a declared collision class is supported
bool BeamSupport::Supports(const CollisionType type) const {
  return type != CollisionType::Invalid && std::find(collisions.begin(), collisions.end(), type) != collisions.end();
}

// Compute support for the physical beam PDGs including charge conjugation
bool BeamSupport::Supports(const int first, const int second) const {
  return Supports(ClassifyCollision(first, second)) && (antiprotons || (first != -PDG::PDG_p && second != -PDG::PDG_p));
}

// Compute all supported collision classes for one process table row
std::string BeamSupport::Names() const {
  constexpr std::array<CollisionType, 6> types = {CollisionType::PP, CollisionType::EE, CollisionType::EP,
                                                  CollisionType::PA, CollisionType::EA, CollisionType::AA};
  std::string                            output;
  for (const auto type : types) {
    if (!Supports(type)) { continue; }
    if (!output.empty()) { output += ", "; }
    output += CollisionName(type);
  }
  return output;
}

// Compute the spatial momentum transfer in the incoming-particle rest frame
// q_rest = Lambda(-p_beam/E_beam) q
M3Vec RestTransfer(const M4Vec &beam, const M4Vec &transfer, const double mass) {
  const auto        p               = beam.Contravariant<long double>();
  const long double mass2           = static_cast<long double>(mass) * mass;
  const long double shell           = beam.Invariant<long double>();
  const long double shell_scale     = std::max(1.0L, SquaredNorm(p) + std::abs(mass2));
  const long double shell_tolerance = 128.0L * std::numeric_limits<double>::epsilon() * shell_scale;
  if (!Finite(beam) || !Finite(transfer) || !(beam.E() > 0.0) || !(mass > 0.0) || !std::isfinite(mass) ||
      !std::isfinite(shell) || !(shell > 0.0L) || std::abs(shell - mass2) > shell_tolerance) {
    throw std::invalid_argument("RestTransfer: beam and transfer must be finite and beam timelike");
  }
  M4Vec rest = transfer;
  kinematics::LorentzBoost(beam, mass, rest, -1);
  if (!Finite(rest)) { throw std::runtime_error("RestTransfer: Lorentz boost failed"); }
  return rest.P3();
}

// Validate nuclear transfers after subtracting high-energy beam momenta
// q_i = p_i - p'_i with q_1 + q_2 = p_X
bool CloseTransfers(const std::array<int, 2> &pdg, const std::array<M4Vec, 2> &beam,
                    const std::array<M4Vec, 2> &forward, const M4Vec &central, std::array<M4Vec, 2> &transfer) {
  transfer = {beam[0] - forward[0], beam[1] - forward[1]};
  if (!IsChargedNucleus(pdg[0]) && !IsChargedNucleus(pdg[1])) { return true; }
  if (!Finite(beam[0]) || !Finite(beam[1]) || !Finite(forward[0]) || !Finite(forward[1]) || !Finite(central) ||
      !Finite(transfer[0]) || !Finite(transfer[1])) {
    return false;
  }
  const M4Vec  residual  = central - transfer[0] - transfer[1];
  const double scale     = std::max({1.0, MaxComponent(beam[0]), MaxComponent(beam[1]), MaxComponent(forward[0]),
                                     MaxComponent(forward[1]), MaxComponent(central)});
  const double tolerance = std::max(1.0e-6, 128.0 * std::numeric_limits<double>::epsilon() * scale);
  if (!Finite(residual) || MaxComponent(residual) > tolerance) { return false; }

  // Preserve the directly resolved nonnuclear transfer and close on the nuclear leg
  if (IsChargedNucleus(pdg[0]) && !IsChargedNucleus(pdg[1])) {
    transfer[0] = central - transfer[1];
  } else {
    transfer[1] = central - transfer[0];
  }
  return Finite(transfer[0]) && Finite(transfer[1]);
}

}  // namespace gra::nuclear
