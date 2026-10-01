// Conditional nuclear impulse fragmentation through the common HepMC3 reaction API
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MNUCLEARINCOHERENT_H
#define MNUCLEARINCOHERENT_H

#include "Graniitti/Nuclear/MEMD.h"
#include "Graniitti/Nuclear/MEvaporation.h"
#include "Graniitti/Nuclear/MFinal.h"

namespace gra::nuclear {

// Resolve a quasifree nucleon and its evaporating spectator remnant
class MIncoherent final : public MReaction {
 public:
  // Use the same nuclear geometry and GDR inputs as the production calculation
  explicit MIncoherent(std::shared_ptr<const MUPC> upc);

  // Copy immutable physics inputs without sharing an event proposal
  std::unique_ptr<MReaction> Clone() const override;

  // Identify the impulse and electromagnetic excitation approximations
  HepMC3::GenRunInfo::ToolInfo Tool() const override;

  // Propose coherent or continuum masses before solving the hard kinematics
  RecoilMass SampleMasses(const HepMC3::GenEvent& beams, MRandom& random) override;

  // Condition the knockout on the actual transfer and decay the spectator
  double Complete(HepMC3::GenEvent& event, const std::array<HepMC3::GenParticlePtr, 2>& forward,
                  MRandom& random) override;

 private:
  std::array<std::shared_ptr<const MEvaporation>, 2>  decay_;
  std::array<bool, 2>                                 active_{}, charge_{};
  std::pair<std::vector<double>, std::vector<double>> rule_;
  std::shared_ptr<const MEMD>                         emd_;
  MEMD::State                                         excitation_;
  struct State {
    bool   excited = false, proton = false;
    double p = 0.0, mass = 0.0, residue = 0.0, pdf = 0.0;
  };
  std::array<State, 2> state_{};

  // Integrate the impulse strength over both occupied momentum spheres
  double Norm(std::size_t leg, const NucleusID& ion, double mass, double t) const;
};

}  // namespace gra::nuclear

#endif
