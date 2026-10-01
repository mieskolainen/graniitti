// Generic final-state parton proposal operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPARTONPROPOSAL_H
#define MPARTONPROPOSAL_H

// C++
#include <string>

// Own
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/Process/MSubProc.h"

namespace gra::parton {

// Configure generated decay parameters and outgoing quark masses
void ConfigureFinalStateMasses(LORENTZSCALAR &lts, MSubProc &subprocess);

// Cache physical alias modes accepted by the active amplitude definition
void ConfigureFinalStateProposal(LORENTZSCALAR &lts, MSubProc &subprocess);

// Propose one physical on-shell mode for generic final-state parton aliases
bool PrepareFinalStateProposal(LORENTZSCALAR &lts, MRandom &random);

// Compute the inverse probability of the current physical parton proposal
double ProposalWeight(const LORENTZSCALAR &lts);

} // namespace gra::parton

#endif
