// Generated event invariants and forward-leg kinematics
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MEVENTKINEMATICS_H
#define MEVENTKINEMATICS_H

// C++
#include <array>
#include <cstddef>
#include <memory>
#include <stdexcept>

// Own
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {

namespace kinematics {

// Select the physical virtuality prescription for forward transfer momenta
enum class TransferVirtualityPolicy { Spacelike, CollinearLightlike, HardDiffractive };

// Reconstruct all Lorentz scalars from one generated event state
bool SetLorentzScalars(MProcessState &state, unsigned int final_count, bool loop_active = false,
                       TransferVirtualityPolicy transfer_policy = TransferVirtualityPolicy::Spacelike);

// Reconstruct exact on-shell forward momenta for one screening-loop transfer
bool RebuildScreeningKinematics(LORENTZSCALAR &lts, const std::array<double, 2> &upper_transverse,
                                const std::array<double, 2> &lower_transverse, bool central_system_active);

}  // namespace kinematics

// Select one generated beam side
enum class ForwardBeamLeg { Upper, Lower };

// Select the physical final state carried by one forward leg
enum class ForwardFinalState { Elastic, InclusiveExcitation };

// Resolve generated recoil kinematics and the corresponding hard nuclear current
struct ForwardLegState {
  ForwardBeamLeg                       leg         = ForwardBeamLeg::Upper;
  ForwardFinalState                    final_state = ForwardFinalState::Elastic;
  bool                                 has_xi      = false;
  double                               xi          = 0.0;
  double                               t           = 0.0;
  double                               qt          = 0.0;
  double                               mass2       = 0.0;
  MParticle                            emitter;
  form::ParamStore                     structure;
  bool                                 is_nuclear = false;
  nuclear::CoherenceType               emission   = nuclear::CoherenceType::Coherent;
  nuclear::CoherenceType               target     = nuclear::CoherenceType::Coherent;
  nuclear::NeutronSelection                 neutron    = {};
  std::shared_ptr<const nuclear::MUPC> upc_model;
  M4Vec                                incoming;
  M4Vec                                outgoing;
  M4Vec                                transfer;

  // Compute whether the generated forward state is inclusively excited
  bool IsExcited() const noexcept { return final_state == ForwardFinalState::InclusiveExcitation; }

  // Check whether the photon emitter is a nucleus
  bool IsNuclear() const noexcept { return is_nuclear; }

  // Compute whether the selected nuclear photon emission is incoherent
  bool IsNuclearIncoherent() const noexcept { return is_nuclear && emission == nuclear::CoherenceType::Incoherent; }

  // Compute whether the selected photonuclear target is incoherent
  bool IsTargetIncoherent() const noexcept { return is_nuclear && target == nuclear::CoherenceType::Incoherent; }

  // Compute the conventional one-based beam-leg index
  int Index() const {
    if (leg == ForwardBeamLeg::Upper) { return 1; }
    if (leg == ForwardBeamLeg::Lower) { return 2; }
    throw std::invalid_argument("ForwardLegState::Index: invalid beam leg");
  }

  // Compute the outgoing forward-system position in pfinal
  std::size_t FinalIndex() const { return static_cast<std::size_t>(Index()); }
};

// Resolve and validate one generated forward leg from event kinematics
ForwardLegState ResolveForwardLegState(const LORENTZSCALAR &lts, ForwardBeamLeg leg);

}  // namespace gra

#endif
