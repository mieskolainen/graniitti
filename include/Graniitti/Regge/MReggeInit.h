// Regge process initialization and prepared amplitude data
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MREGGEINIT_H
#define MREGGEINIT_H

#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Process/MProcessSetup.h"
#include "Graniitti/Regge/MRegge.h"

namespace gra {

// Compute the diphoton scale from its physical transverse norm
double GammaGammaScale(const PARAM_RES &res, double norm2, bool isolated_decay);
// Compute the sign of an integral spin phase
int IntegerPhaseSign(double exponent, const std::string &context);

// Store the active branching modes selected during process initialization
struct BranchingProcessFlags {
  ReggeProductionModel model               = ReggeProductionModel::None;
  bool                 resonance           = false;
  bool                 continuum           = false;
  bool                 jw_helicity_algebra = false;

  // Compute whether resonance dynamics uses one production model
  bool ResonanceModel(ReggeProductionModel candidate) const { return resonance && model == candidate; }

  // Compute whether continuum dynamics uses one production model
  bool ContinuumModel(ReggeProductionModel candidate) const { return continuum && model == candidate; }
};

// Compute the validated immutable soft model of one process setup
const SoftModelPtr &RequireSetupSoftModel(const MProcessSetup &setup);

// Compute the support restriction for one active photoproduction channel
std::string PhotoProductionError(ReggeProductionModel model, nuclear::CollisionType collision,
                               const MParticle &resonance, const std::array<int, 2> &exchange,
                               const MPDG &pdg, const nlohmann::json &general);

// Classify one subprocess into its active branching modes
BranchingProcessFlags ClassifyBranchingProcess(const MProcessSetup &setup);

// Initialize the shared physical Regge or tensor process state
void InitializePhysicalProcessState(MProcessSetup &setup, const BranchingProcessFlags &flags);

// Prepare the shared physical spin spaces for integrated Regge metrics
void ConfigureReggeQMetrics(MProcessSetup &setup, ReggeProductionModel model);

// Build the immutable two-body continuum production plan
void BuildContinuum2Plan(MProcessSetup &setup, ReggeProductionModel model);

// Build the immutable four-body continuum production plan
void BuildContinuum4Plan(MProcessSetup &setup, ReggeProductionModel model);

// Build the immutable six-body continuum production plan
void BuildContinuum6Plan(MProcessSetup &setup, ReggeProductionModel model);

}  // namespace gra

#endif
