// Two-body Regge resonance and continuum amplitudes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGECEP_H
#define MREGGECEP_H

#include <array>
#include <cstddef>
#include <optional>
#include <string>

#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeSource.h"

namespace gra::reggeamp {

// Compute the root decay frame sewn to one Regge production basis
std::string RootFrame(const LORENTZSCALAR &lts, ReggeProductionModel model);

// Resolve one transfer-independent decay matrix for the active event
const MReggeDecayState *Decay(LORENTZSCALAR &lts, const std::string &name, const PARAM_RES &resonance, const std::string &root_frame, MReggeDecayState &local);

// Evaluate one selected model-specific resonance amplitude
void EvalResonance(const MRegge &regge, LORENTZSCALAR &lts, const regge::Param &param, const std::string &name, PARAM_RES &resonance, ReggeProductionModel model, gpom::AmpCache &gp_amp, ReggeSourceCache &source_cache,
                   std::size_t coherence_group, std::optional<std::size_t> channel_filter = std::nullopt, const MReggeDecayState *prepared_decay = nullptr, const std::array<ForwardLegState, 2> *forward_state = nullptr);

}  // namespace gra::reggeamp

#endif
