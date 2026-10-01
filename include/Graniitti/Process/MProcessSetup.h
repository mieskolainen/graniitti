// Process initialization context for physical amplitude classes
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef MPROCESSSETUP_H
#define MPROCESSSETUP_H

// C++
#include <functional>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Process/MProcessState.h"
#include "Graniitti/MModelTune.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Regge/MSoftModel.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MHELMatrix.h"
#include "Graniitti/Spin/MHelicity.h"

namespace gra {

// Select whether a physical process accepts free resonance cards
enum class ResonanceType { None, Optional, Required };

// Declare process capabilities shared by steering, initialization and help
struct ProcessInfo {
  ReggeProductionModel model = ReggeProductionModel::None;
  ResonanceType resonance = ResonanceType::None;
  bool continuum = false;
  bool jw_helicity_algebra = true;
  nuclear::BeamSupport beams;
  RootDecayMode example_decay = RootDecayMode::Physical;
};

// Expose process-independent initialization services to physical amplitudes
struct MProcessSetup {
  using HelicityStructureFunction = std::function<HELMatrix(
      const MParticle &, const std::vector<MParticle> &, bool, bool,
      const std::string &, bool, gra::spin::VertexContext)>;
  using HelicityTreeFunction = std::function<void(MDecayBranch &, bool, bool)>;
  using PoleOperatorStructureFunction =
      std::function<spin::PoleLS(
          ReggeProductionModel, const MParticle &,
          const std::vector<MParticle> &, gra::spin::VertexContext)>;

  LORENTZSCALAR &lts;
  MModelTunePtr model_tune;
  SoftModelPtr soft_model;
  std::string istate;
  std::string channel;
  const ProcessInfo info;
  bool flat_mass2 = false;
  bool isolate = false;
  int excitation = 0;
  bool spingen_user = false;
  bool spindec_user = false;
  bool mmax_user = false;
  double &symmetry_factor;
  HelicityStructureFunction helicity_structure;
  HelicityTreeFunction helicity_tree;
  PoleOperatorStructureFunction pole_operator_structure;

  // Construct one generic process initialization context
  MProcessSetup(LORENTZSCALAR &state, MModelTunePtr tune,
                std::string initial_state, std::string process_channel, ProcessInfo process_info,
                bool use_flat_mass2, bool isolate_decay, int forward_excitation,
                bool has_spingen_override, bool has_spindec_override,
                bool has_mmax_override, double &statistical_factor,
                HelicityStructureFunction structure_function,
                HelicityTreeFunction tree_function,
                PoleOperatorStructureFunction pole_operator_function)
      : lts(state), model_tune(std::move(tune)),
        soft_model(model_tune != nullptr ? model_tune->Soft() : nullptr),
        istate(std::move(initial_state)), channel(std::move(process_channel)),
        info(std::move(process_info)),
        flat_mass2(use_flat_mass2), isolate(isolate_decay),
        excitation(forward_excitation), spingen_user(has_spingen_override),
        spindec_user(has_spindec_override), mmax_user(has_mmax_override),
        symmetry_factor(statistical_factor),
        helicity_structure(std::move(structure_function)),
        helicity_tree(std::move(tree_function)),
        pole_operator_structure(std::move(pole_operator_function)) {}

  // Construct one generic helicity vertex
  HELMatrix ProcessHelicityStructure(
      const MParticle &particle, const std::vector<MParticle> &legs,
      bool production_mode = false, bool strict_mode = false,
      const std::string &verbose_label = "", bool verbose_output = true,
      gra::spin::VertexContext context = gra::spin::VertexContext::Auto) const {
    return helicity_structure(particle, legs, production_mode, strict_mode,
                              verbose_label, verbose_output, context);
  }

  // Construct generic helicity vertices recursively through one decay branch
  void ProcessHelicityTree(MDecayBranch &branch, bool production_mode = false,
                           bool strict_mode = false) const {
    helicity_tree(branch, production_mode, strict_mode);
  }

  // Construct one canonical MP or XP continuum pole operator
  spin::PoleLS
  ProcessPoleOperatorStructure(ReggeProductionModel model,
                               const MParticle &particle,
                               const std::vector<MParticle> &legs,
                               gra::spin::VertexContext context) const {
    return pole_operator_structure(model, particle, legs, context);
  }
};

// Initialize process-independent root decay and phase-space state
void InitializeGenericProcessState(MProcessSetup &setup);

} // namespace gra

#endif
