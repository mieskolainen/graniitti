// Nuclear UPC event-record metadata
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MRecord.h"

#include <cmath>
#include <memory>
#include <string>
#include <utility>

#include "Graniitti/Tech/MAux.h"
#include "HepMC3/Attribute.h"
#include "HepMC3/GenHeavyIon.h"

namespace gra::nuclear {

using gra::aux::indices;

namespace {

// Build one stable per-leg HepMC attribute name
std::string LegKey(const std::string &name, const int leg) {
  return "graniitti_upc_" + name + "_" + std::to_string(leg);
}

// Compute spectator neutron and proton counts for one UPC beam
std::pair<int, int> Spectators(const MUPC &upc, const int leg) {
  const MNucleus *nucleus = upc.Nucleus(leg);
  if (nucleus != nullptr) { return {static_cast<int>(nucleus->A() - nucleus->Z()), static_cast<int>(nucleus->Z())}; }
  if (upc.Type(leg) == BeamType::Proton) { return {0, 1}; }
  return {-1, -1};
}

// Construct the standard HepMC heavy-ion information for one UPC event
std::shared_ptr<HepMC3::GenHeavyIon> HeavyIonRecord(const MUPC &upc, const double impact_parameter) {
  auto record                          = std::make_shared<HepMC3::GenHeavyIon>();
  record->Ncoll_hard                   = 0;
  record->Ncoll                        = 0;
  record->N_Nwounded_collisions        = 0;
  record->Nwounded_N_collisions        = 0;
  record->Nwounded_Nwounded_collisions = 0;
  record->Npart_proj                   = upc.Type(1) == BeamType::Lepton ? -1 : 0;
  record->Npart_targ                   = upc.Type(2) == BeamType::Lepton ? -1 : 0;

  const auto projectile = Spectators(upc, 1);
  const auto target     = Spectators(upc, 2);
  record->Nspec_proj_n  = projectile.first;
  record->Nspec_proj_p  = projectile.second;
  record->Nspec_targ_n  = target.first;
  record->Nspec_targ_p  = target.second;
  if (std::isfinite(impact_parameter) && impact_parameter >= 0.0) { record->impact_parameter = impact_parameter; }
  const double sigma_inel = upc.Param().glauber.profile.sigma;
  if (upc.HasHadronicPair() && std::isfinite(sigma_inel) && sigma_inel > 0.0) { record->sigma_inel_NN = sigma_inel; }
  return record;
}

}  // namespace

// Attach complete UPC model metadata to one HepMC event
void AttachRecord(const MUPC &upc, const FinalState &final, HepMC3::GenEvent &event) {
  const UPCParam &param = upc.Param();
  event.set_heavy_ion(HeavyIonRecord(upc, -1.0));
  event.add_attribute("graniitti_upc_survival",
                      std::make_shared<HepMC3::StringAttribute>(SurvivalName(param.survival)));
  event.add_attribute("graniitti_upc_hadronic_screening",
                      std::make_shared<HepMC3::IntAttribute>(upc.HadronicConvolution()));
  event.add_attribute("graniitti_upc_structure",
                      std::make_shared<HepMC3::StringAttribute>(StructureName(param.structure)));
  event.add_attribute("graniitti_upc_structure_samples",
                      std::make_shared<HepMC3::IntAttribute>(static_cast<int>(param.config.count)));
  event.add_attribute("graniitti_upc_sample_count",
                      std::make_shared<HepMC3::IntAttribute>(static_cast<int>(upc.SampleCount())));
  event.add_attribute("graniitti_upc_ggcf_states",
                      std::make_shared<HepMC3::IntAttribute>(static_cast<int>(param.glauber.fluctuation.nodes)));
  event.add_attribute("graniitti_upc_nn_eikonal", std::make_shared<HepMC3::StringAttribute>(param.survival_eikonal));
  event.add_attribute("graniitti_upc_nn_omega", std::make_shared<HepMC3::DoubleAttribute>(param.glauber.profile.omega));
  event.add_attribute("graniitti_upc_nn_s", std::make_shared<HepMC3::DoubleAttribute>(param.s_nn));
  event.add_attribute("graniitti_upc_photo_model", std::make_shared<HepMC3::StringAttribute>(PhotoModelName(param.photo_model)));
  for (const auto &index : indices(param.emission)) {
    const int leg = static_cast<int>(index + 1);
    event.add_attribute(LegKey("beam", leg), std::make_shared<HepMC3::StringAttribute>(BeamName(upc.Type(leg))));
    event.add_attribute(LegKey("emission", leg),
                        std::make_shared<HepMC3::StringAttribute>(CoherenceName(param.emission[index])));
    event.add_attribute(LegKey("target", leg),
                        std::make_shared<HepMC3::StringAttribute>(CoherenceName(param.target[index])));
    event.add_attribute(LegKey("neutrons", leg),
                        std::make_shared<HepMC3::StringAttribute>(NeutronName(param.neutron[index])));
    event.add_attribute(LegKey("structure_samples", leg),
                        std::make_shared<HepMC3::IntAttribute>(static_cast<int>(upc.SampleShape()[index])));
    event.add_attribute(
        LegKey("final_sector", leg),
        std::make_shared<HepMC3::StringAttribute>(final.valid ? CoherenceName(final.leg[index]) : "unspecified"));
  }
}

}  // namespace gra::nuclear
