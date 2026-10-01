// Nuclear UPC steering parser
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Nuclear/MSteering.h"

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/MGlobals.h"
#include "Graniitti/Tech/MAux.h"

namespace gra::nuclear {
namespace {

using gra::aux::indices;

// Require one object to contain exactly the declared keys
void RequireKeys(const nlohmann::json &block, const std::set<std::string> &required, const std::string &context) {
  if (!block.is_object()) { throw std::invalid_argument("ReadUPCSteering: " + context + " must be an object"); }
  for (const auto &key : required) {
    if (!block.contains(key)) { throw std::invalid_argument("ReadUPCSteering: missing " + context + "::" + key); }
  }
  for (auto item = block.begin(); item != block.end(); ++item) {
    if (required.count(item.key()) == 0) {
      throw std::invalid_argument("ReadUPCSteering: unknown " + context + " key '" + item.key() + "'");
    }
  }
}

// Read ordered mass tables and the published formula used only outside their coverage
MassParam ReadMass(const nlohmann::json& block, const std::string& tune) {
  RequireKeys(block, {"tables", "bwm"}, "FRAGMENTATION::mass");
  RequireKeys(block.at("bwm"), {"coefficients", "asymmetry_scale", "pairing_scale"}, "FRAGMENTATION::mass::bwm");
  if (!block.at("tables").is_array() || block.at("tables").empty()) {
    throw std::invalid_argument("ReadUPCSteering: mass tables must be a nonempty array");
  }
  MassParam param;
  for (const auto& row : block.at("tables")) {
    RequireKeys(row, {"type", "file"}, "FRAGMENTATION::mass::tables[]");
    const auto type = row.at("type").get<std::string>();
    if (type != "ame2020" && type != "frdm2012") { throw std::invalid_argument("ReadUPCSteering: unknown mass table type"); }
    param.tables.push_back({type == "ame2020" ? MassType::AME : MassType::FRDM,
                            gra::ResolveModelDataFile(tune, row.at("file").get<std::string>())});
  }
  param.bwm = block.at("bwm").at("coefficients").get<std::array<double, 5>>();
  param.asymmetry = block.at("bwm").at("asymmetry_scale").get<double>();
  param.pairing = block.at("bwm").at("pairing_scale").get<double>();
  return param;
}

// Parse one integer count without truncating fractions or wrapping negative inputs
template <class T>
T ReadCount(const nlohmann::json& value) {
  if (!value.is_number_integer() || value < 0 || value.get<std::uint64_t>() > std::numeric_limits<T>::max()) {
    throw std::invalid_argument("ReadUPCSteering: count must be a nonnegative integer within range");
  }
  return value.get<T>();
}

// Read one required pair of lowercase nuclear leg selectors
std::vector<std::string> ReadLegNames(const nlohmann::json& block, const std::string& name) {
  const auto values = block.at(name).get<std::vector<std::string>>();
  if (values.size() != 2) {
    throw std::invalid_argument("ReadUPCSteering: NUCLEAR::" + name + " must contain two lowercase values");
  }
  return values;
}

// Parse one moment-matched Good-Walker cross-section rule
FluctuationParam ReadFluctuation(const unsigned int nodes, const nlohmann::json& numerics) {
  FluctuationParam param;
  param.nodes            = nodes;
  param.cdf_tol          = numerics.at("cdf_tol").get<double>();
  param.cdf_max_iter     = ReadCount<unsigned int>(numerics.at("cdf_max_iter"));
  param.normal_shape_min = numerics.at("normal_shape_min").get<double>();
  return param;
}

// Parse the generic density rule and explicit nuclear shapes
GeometryParam ReadGeometry(const nlohmann::json &physics, const nlohmann::json &numerics) {
  GeometryParam geometry;
  geometry.radius_scale = physics.at("radius_scale").get<double>();
  geometry.skin         = physics.at("skin").get<double>();
  geometry.nodes        = ReadCount<unsigned int>(numerics.at("density_nodes"));
  geometry.tail_skin    = numerics.at("tail_skin").get<double>();
  geometry.cdf_nodes    = ReadCount<std::size_t>(numerics.at("cdf_nodes"));
  geometry.form_q_max   = numerics.at("form_q_max").get<double>();
  geometry.form_nodes   = ReadCount<std::size_t>(numerics.at("form_nodes"));
  geometry.form_abs_tol = numerics.at("form_abs_tol").get<double>();
  if (!physics.at("nuclei").is_array()) {
    throw std::invalid_argument("ReadUPCSteering: PARAM_NUCLEAR::STRUCTURE::nuclei must be an array");
  }
  std::set<std::pair<unsigned int, unsigned int>> identity;
  for (const auto &item : physics.at("nuclei")) {
    RequireKeys(item, {"A", "Z", "charge", "matter"}, "PARAM_NUCLEAR::STRUCTURE::nuclei[]");
    NucleusShape shape;
    shape.a      = ReadCount<unsigned int>(item.at("A"));
    shape.z      = ReadCount<unsigned int>(item.at("Z"));
    shape.charge = item.at("charge").get<std::array<double, 2>>();
    shape.matter = item.at("matter").get<std::array<double, 2>>();
    if (shape.a == 0 || shape.z > shape.a || !identity.insert({shape.a, shape.z}).second) {
      throw std::invalid_argument("ReadUPCSteering: invalid or duplicate nuclear geometry");
    }
    geometry.nucleus.push_back(shape);
  }
  return geometry;
}

// Parse one isotope-specific giant dipole response
GDRIsotope ReadGDRIsotope(const nlohmann::json &item) {
  RequireKeys(item, {"A", "Z", "energy", "strength", "width"}, "PARAM_NUCLEAR::EMD::gdr::isotopes[]");
  return {ReadCount<unsigned int>(item.at("A")), ReadCount<unsigned int>(item.at("Z")), item.at("energy").get<double>(),
          item.at("width").get<double>(), item.at("strength").get<double>()};
}

// Parse one target isotope response and its numerical controls
BreakupParam ReadEMDIsotope(const nlohmann::json &physics, const nlohmann::json &numerics) {
  RequireKeys(physics, {"A", "Z", "photoabsorption", "energy_transfer"}, "PARAM_NUCLEAR::EMD::isotopes[]");
  const auto &photo = physics.at("photoabsorption");
  RequireKeys(photo, {"energy_min", "quasi_deuteron", "resonance", "continuum"}, "EMD::isotopes[]::photoabsorption");
  const auto &qd = photo.at("quasi_deuteron");
  RequireKeys(qd, {"deuteron", "levinger", "pauli"}, "EMD::quasi_deuteron");
  const auto &deuteron = qd.at("deuteron");
  RequireKeys(deuteron, {"norm", "threshold"}, "EMD::quasi_deuteron::deuteron");
  const auto &pauli = qd.at("pauli");
  RequireKeys(pauli, {"high_edge", "high_exp", "low_edge", "low_exp", "poly"}, "EMD::quasi_deuteron::pauli");
  const auto &resonance = photo.at("resonance");
  RequireKeys(resonance, {"threshold", "energy", "width", "area"}, "EMD::resonance");
  const auto &continuum = photo.at("continuum");
  RequireKeys(continuum, {"constant", "match", "log2", "mean", "norm", "omega0", "threshold", "width"}, "EMD::continuum");
  const auto &transfer = physics.at("energy_transfer");
  RequireKeys(transfer, {"E0", "match"}, "EMD::energy_transfer");
  BreakupParam emd;
  emd.a = ReadCount<unsigned int>(physics.at("A"));
  emd.z = ReadCount<unsigned int>(physics.at("Z"));
  emd.photo.energy_min = photo.at("energy_min").get<double>();
  emd.photo.qd.deuteron = {deuteron.at("threshold").get<double>(), deuteron.at("norm").get<double>()};
  emd.photo.qd.levinger = qd.at("levinger").get<double>();
  emd.photo.qd.pauli = {pauli.at("low_edge").get<double>(), pauli.at("high_edge").get<double>(),
                       pauli.at("low_exp").get<double>(), pauli.at("poly").get<std::array<double, 5>>(),
                       pauli.at("high_exp").get<double>()};
  emd.photo.resonance.threshold = resonance.at("threshold").get<double>();
  emd.photo.resonance.isotope = {{emd.a, emd.z, resonance.at("energy").get<double>(),
                                 resonance.at("width").get<double>(), resonance.at("area").get<double>()}};
  emd.photo.continuum = {emd.a, emd.z, continuum.at("threshold").get<double>(), continuum.at("match").get<double>(),
                         continuum.at("omega0").get<double>(), continuum.at("constant").get<double>(),
                         continuum.at("log2").get<double>(), continuum.at("norm").get<double>(),
                         continuum.at("mean").get<double>(), continuum.at("width").get<double>()};
  emd.transfer = {transfer.at("E0").get<double>(), transfer.at("match").get<double>()};
  emd.response = {ReadCount<unsigned int>(numerics.at("response").at("nodes")),
                  numerics.at("response").at("rel_tol").get<double>(),
                  numerics.at("response").at("tail_rel_tol").get<double>()};
  emd.profile = {ReadCount<std::size_t>(numerics.at("profile").at("nodes")),
                 numerics.at("profile").at("b_min").get<double>()};
  emd.proton = {numerics.at("proton_emitter").at("q_max").get<double>(),
                ReadCount<unsigned int>(numerics.at("proton_emitter").at("q_nodes"))};
  return emd;
}

// Parse all target responses and shared E1 systematics before model construction
EMDParam ReadEMD(const nlohmann::json &physics, const nlohmann::json &numerics) {
  RequireKeys(physics, {"gdr", "isotopes"}, "PARAM_NUCLEAR::EMD");
  RequireKeys(numerics, {"profile", "proton_emitter", "response", "sampling"}, "NUMERICS_NUCLEAR::EMD");
  RequireKeys(numerics.at("sampling"), {"focus_fraction", "nonzero_photons"}, "NUMERICS_NUCLEAR::EMD::sampling");
  RequireKeys(numerics.at("response"), {"nodes", "rel_tol", "tail_rel_tol"}, "NUMERICS_NUCLEAR::EMD::response");
  RequireKeys(numerics.at("profile"), {"b_min", "nodes"}, "NUMERICS_NUCLEAR::EMD::profile");
  RequireKeys(numerics.at("proton_emitter"), {"q_max", "q_nodes"}, "NUMERICS_NUCLEAR::EMD::proton_emitter");
  const auto &gdr = physics.at("gdr");
  RequireKeys(gdr, {"systematics", "isotopes"}, "PARAM_NUCLEAR::EMD::gdr");
  const auto &systematics = gdr.at("systematics");
  RequireKeys(systematics, {"energy_a", "energy_b", "width_norm", "width_power", "trk", "strength"}, "EMD::gdr::systematics");
  if (!gdr.at("isotopes").is_array() || !physics.at("isotopes").is_array()) {
    throw std::invalid_argument("ReadUPCSteering: EMD isotope responses must be arrays");
  }
  EMDParam emd;
  emd.gdr.systematics = {systematics.at("energy_a").get<double>(), systematics.at("energy_b").get<double>(),
                         systematics.at("width_norm").get<double>(), systematics.at("width_power").get<double>(),
                         systematics.at("trk").get<double>(), systematics.at("strength").get<double>()};
  for (const auto &item : gdr.at("isotopes")) { emd.gdr.isotope.push_back(ReadGDRIsotope(item)); }
  ValidateGDR(emd.gdr);
  std::set<std::pair<unsigned int, unsigned int>> identity;
  for (const auto &item : physics.at("isotopes")) {
    auto response = ReadEMDIsotope(item, numerics);
    if (response.a <= 1 || response.z == 0 || response.z >= response.a ||
        !identity.insert({response.a, response.z}).second) {
      throw std::invalid_argument("ReadUPCSteering: invalid or duplicate EMD isotope");
    }
    emd.isotope.push_back(std::move(response));
  }
  return emd;
}

}  // namespace

// Parse selectors and complete nuclear model-card blocks
UPCParam ReadUPCSteering(const nlohmann::json &block, const nlohmann::json &physics, const nlohmann::json &numerics,
                         const std::string &tune) {
  RequireKeys(
      block,
      {"emd", "fragmentation", "emission", "neutron_class", "photoproduction", "sigma_NN", "structure", "screening"},
      "NUCLEAR");
  RequireKeys(block.at("photoproduction"), {"target", "target_model"}, "NUCLEAR::photoproduction");
  RequireKeys(physics, {"EMD", "FRAGMENTATION", "GGCF", "LTA", "STRUCTURE"}, "PARAM_NUCLEAR");
  RequireKeys(physics.at("FRAGMENTATION"), {"mass", "external"}, "PARAM_NUCLEAR::FRAGMENTATION");
  RequireKeys(
      numerics,
      {"EMD", "FLUCTUATION", "FRAGMENTATION", "GLAUBER", "IMPACT_PROFILE", "PHOTOPROD", "SAMPLING", "STRUCTURE"},
      "NUMERICS_NUCLEAR");
  RequireKeys(physics.at("STRUCTURE"), {"hotspot", "nuclei", "nucleon_d_min", "radius_scale", "skin"},
              "PARAM_NUCLEAR::STRUCTURE");
  RequireKeys(physics.at("STRUCTURE").at("hotspot"), {"b_center", "b_profile", "count", "strength_sigma"},
              "PARAM_NUCLEAR::STRUCTURE::hotspot");
  RequireKeys(physics.at("LTA"),
              {"alpha0", "alpha_prime", "b_diff", "dpdf_diss", "dpdf_member", "dpdf_set", "flux_b", "pdf_member",
               "pdf_set", "q2_max", "q2_min", "sigma3_anchor_x", "sigma3_fade_x", "sigma3_power", "x_max", "x_min"},
              "PARAM_NUCLEAR::LTA");
  RequireKeys(
      numerics.at("STRUCTURE"),
      {"cdf_nodes", "density_nodes", "form_abs_tol", "form_nodes", "form_q_max", "tail_skin"},
      "NUMERICS_NUCLEAR::STRUCTURE");
  RequireKeys(numerics.at("SAMPLING"),
              {"config_samples", "current_samples", "hard_core_sweeps", "neutron_density_rel_tol", "neutron_nodes",
               "neutron_negative_norm_tol", "neutron_norm_tol", "placement_trials"},
              "NUMERICS_NUCLEAR::SAMPLING");
  RequireKeys(numerics.at("FLUCTUATION"), {"cdf_tol", "cdf_max_iter", "normal_shape_min"},
              "NUMERICS_NUCLEAR::FLUCTUATION");
  RequireKeys(numerics.at("GLAUBER"), {"b_max", "b_nodes", "ggcf_nodes", "profile_nodes", "q_max", "q_nodes"},
              "NUMERICS_NUCLEAR::GLAUBER");
  RequireKeys(numerics.at("PHOTOPROD"),
              {"cf_nodes", "lta_integral_nodes", "lta_q2_nodes", "lta_x_nodes", "phase_scale", "table_qt_max",
               "table_qt_nodes", "table_qz_max", "table_qz_nodes", "table_series_abs_tol", "table_series_terms"},
              "NUMERICS_NUCLEAR::PHOTOPROD");
  RequireKeys(numerics.at("IMPACT_PROFILE"),
              {"b_max", "b_phi_nodes", "kt_max", "kt_min", "kt_nodes", "log_kt", "phi_nodes", "sample_b_nodes",
               "sample_kt_nodes", "smooth_b_nodes", "smooth_kt_nodes"},
              "NUMERICS_NUCLEAR::IMPACT_PROFILE");

  UPCParam param;
  param.additional_emd      = block.at("emd").get<bool>();
  const auto &fragmentation = physics.at("FRAGMENTATION");
  const auto &external      = fragmentation.at("external");
  RequireKeys(external, {"library", "config"}, "FRAGMENTATION::external");
  RequireKeys(numerics.at("FRAGMENTATION"), {"decay_nodes", "decay_steps", "impulse_nodes"},
              "NUMERICS_NUCLEAR::FRAGMENTATION");
  if (!block.at("fragmentation").is_null()) {
    const auto model = block.at("fragmentation").get<std::string>();
    if (model == "internal") {
      param.reaction = ReactionType::Internal;
    } else if (model == "external") {
      param.reaction = ReactionType::External;
    } else {
      throw std::invalid_argument("ReadUPCSteering: NUCLEAR::fragmentation must be null, internal or external");
    }
  }
  param.final_library = external.at("library").get<std::string>();
  param.final_config  = external.at("config").dump();
  param.mass          = ReadMass(fragmentation.at("mass"), tune);
  param.decay_nodes   = ReadCount<unsigned int>(numerics.at("FRAGMENTATION").at("decay_nodes"));
  param.decay_steps   = ReadCount<unsigned int>(numerics.at("FRAGMENTATION").at("decay_steps"));
  param.impulse_nodes = ReadCount<unsigned int>(numerics.at("FRAGMENTATION").at("impulse_nodes"));
  param.emd_focus     = numerics.at("EMD").at("sampling").at("focus_fraction").get<double>();
  param.emd_condition = numerics.at("EMD").at("sampling").at("nonzero_photons").get<bool>();
  if (!external.at("config").is_object() || param.decay_nodes < 16 || param.decay_nodes > 4096 ||
      param.decay_steps == 0 || param.impulse_nodes < 16 || param.impulse_nodes > 4096 ||
      (param.reaction == ReactionType::External && param.final_library.empty())) {
    throw std::invalid_argument("ReadUPCSteering: invalid nuclear reaction inputs or decay quadrature");
  }
  const auto  emission        = ReadLegNames(block, "emission");
  const auto& photoproduction = block.at("photoproduction");
  const auto  target          = ReadLegNames(photoproduction, "target");
  const auto  neutron         = ReadLegNames(block, "neutron_class");
  for (const auto &leg : indices(param.emission)) {
    param.emission[leg] = ParseCoherenceType(emission[leg]);
    param.target[leg]   = ParseCoherenceType(target[leg]);
    param.neutron[leg]  = ParseNeutronSelection(neutron[leg]);
  }
  if (!param.reaction && (param.additional_emd || std::any_of(param.neutron.begin(), param.neutron.end(),
                                                              [](auto n) { return n.type != NeutronSelection::Type::Any; }))) {
    throw std::invalid_argument("ReadUPCSteering: EMD and neutron classes require internal or external fragmentation");
  }
  const auto &ggcf = physics.at("GGCF");
  RequireKeys(ggcf, {"eikonal"}, "PARAM_NUCLEAR::GGCF");
  param.ggcf.eikonal = ggcf.at("eikonal").get<std::string>();
  param.photo_model    = ParsePhotoModel(photoproduction.at("target_model").get<std::string>());
  param.structure      = ParseStructureType(block.at("structure").get<std::string>());
  param.survival       = ParseSurvivalType(block.at("screening").get<std::string>());
  const auto &sigma_nn = block.at("sigma_NN");
  if (!sigma_nn.is_null()) {
    if (!sigma_nn.is_number()) {
      throw std::invalid_argument("ReadUPCSteering: NUCLEAR::sigma_NN must be null or a positive number in mb");
    }
    const double value = sigma_nn.get<double>();
    if (!std::isfinite(value) || !(value > 0.0)) {
      throw std::invalid_argument("ReadUPCSteering: NUCLEAR::sigma_NN must be null or a positive number in mb");
    }
    param.sigma_nn = value;
  }

  const auto& structure     = physics.at("STRUCTURE");
  const auto& structure_num = numerics.at("STRUCTURE");
  const auto& sampling      = numerics.at("SAMPLING");
  param.config.count =
      param.structure == StructureType::Smooth ? 0 : ReadCount<std::size_t>(sampling.at("config_samples"));
  param.current_count =
      param.structure == StructureType::Smooth ? 0 : ReadCount<std::size_t>(sampling.at("current_samples"));
  param.config.d_min             = structure.at("nucleon_d_min").get<double>();
  param.config.sweeps            = ReadCount<std::size_t>(sampling.at("hard_core_sweeps"));
  param.config.max_trials        = ReadCount<std::size_t>(sampling.at("placement_trials"));
  param.config.neutron_nodes     = ReadCount<std::size_t>(sampling.at("neutron_nodes"));
  param.config.density_rel_tol   = sampling.at("neutron_density_rel_tol").get<double>();
  param.config.negative_norm_tol = sampling.at("neutron_negative_norm_tol").get<double>();
  param.config.norm_tol          = sampling.at("neutron_norm_tol").get<double>();
  param.geometry                 = ReadGeometry(structure, structure_num);

  const auto& glauber_num = numerics.at("GLAUBER");
  param.glauber.fluctuation =
      ReadFluctuation(ReadCount<unsigned int>(glauber_num.at("ggcf_nodes")), numerics.at("FLUCTUATION"));
  param.glauber.b_max         = glauber_num.at("b_max").get<double>();
  param.glauber.b_nodes       = ReadCount<unsigned int>(glauber_num.at("b_nodes"));
  param.glauber.q_max         = glauber_num.at("q_max").get<double>();
  param.glauber.q_nodes       = ReadCount<unsigned int>(glauber_num.at("q_nodes"));
  param.glauber.profile_nodes = ReadCount<unsigned int>(glauber_num.at("profile_nodes"));

  const auto& profile_num       = numerics.at("IMPACT_PROFILE");
  param.loop.radial_integrator  = "GL";
  param.loop.azimuth_integrator = "Trap";
  param.loop.radial_map         = profile_num.at("log_kt").get<bool>() ? math::RadialMap::Log : math::RadialMap::Square;
  param.loop.r_min              = profile_num.at("kt_min").get<double>();
  param.loop.r_max              = profile_num.at("kt_max").get<double>();
  param.loop.radial_intervals   = ReadCount<unsigned int>(profile_num.at("kt_nodes"));
  param.loop.azimuth_nodes      = ReadCount<unsigned int>(profile_num.at("phi_nodes"));
  param.convolution.b_max       = profile_num.at("b_max").get<double>();
  param.convolution.smooth_b_nodes  = ReadCount<unsigned int>(profile_num.at("smooth_b_nodes"));
  param.convolution.sample_b_nodes  = ReadCount<unsigned int>(profile_num.at("sample_b_nodes"));
  param.convolution.b_phi_nodes     = ReadCount<unsigned int>(profile_num.at("b_phi_nodes"));
  param.convolution.smooth_kt_nodes = ReadCount<unsigned int>(profile_num.at("smooth_kt_nodes"));
  param.convolution.sample_kt_nodes = ReadCount<unsigned int>(profile_num.at("sample_kt_nodes"));

  param.emd = ReadEMD(physics.at("EMD"), numerics.at("EMD"));

  const auto &photo_num = numerics.at("PHOTOPROD");
  const auto &lta       = physics.at("LTA");
  for (const auto &leg : indices(param.photo)) {
    auto &photo          = param.photo[leg];
    photo.fluctuation    = ReadFluctuation(ReadCount<unsigned int>(photo_num.at("cf_nodes")), numerics.at("FLUCTUATION"));
    photo.phase_scale    = photo_num.at("phase_scale").get<double>();
    photo.b_nodes        = param.geometry.nodes;
    photo.z_nodes        = param.geometry.nodes;
    photo.table.qt_max   = photo_num.at("table_qt_max").get<double>();
    photo.table.qz_max   = photo_num.at("table_qz_max").get<double>();
    photo.table.qt_nodes = ReadCount<unsigned int>(photo_num.at("table_qt_nodes"));
    photo.table.qz_nodes = ReadCount<unsigned int>(photo_num.at("table_qz_nodes"));
    photo.table.series_terms     = ReadCount<unsigned int>(photo_num.at("table_series_terms"));
    photo.table.series_abs_tol   = photo_num.at("table_series_abs_tol").get<double>();
    photo.shadow.pdf_set         = lta.at("pdf_set").get<std::string>();
    photo.shadow.pdf_member      = lta.at("pdf_member").get<int>();
    photo.shadow.dpdf_set        = lta.at("dpdf_set").get<std::string>();
    photo.shadow.dpdf_member     = lta.at("dpdf_member").get<int>();
    photo.shadow.alpha0          = lta.at("alpha0").get<double>();
    photo.shadow.alpha_prime     = lta.at("alpha_prime").get<double>();
    photo.shadow.flux_b          = lta.at("flux_b").get<double>();
    photo.shadow.dpdf_diss       = lta.at("dpdf_diss").get<double>();
    photo.shadow.b_diff          = lta.at("b_diff").get<double>();
    photo.shadow.x_min           = lta.at("x_min").get<double>();
    photo.shadow.x_max           = lta.at("x_max").get<double>();
    photo.shadow.q2_min          = lta.at("q2_min").get<double>();
    photo.shadow.q2_max          = lta.at("q2_max").get<double>();
    photo.shadow.sigma3_anchor_x = lta.at("sigma3_anchor_x").get<double>();
    photo.shadow.sigma3_fade_x   = lta.at("sigma3_fade_x").get<double>();
    photo.shadow.sigma3_power    = lta.at("sigma3_power").get<double>();
    photo.shadow.x_nodes         = ReadCount<unsigned int>(photo_num.at("lta_x_nodes"));
    photo.shadow.scale_nodes     = ReadCount<unsigned int>(photo_num.at("lta_q2_nodes"));
    photo.shadow.integral_nodes  = ReadCount<unsigned int>(photo_num.at("lta_integral_nodes"));
    if (param.structure == StructureType::Hotspot) {
      const auto& hotspot          = structure.at("hotspot");
      photo.hotspot.count          = ReadCount<unsigned int>(hotspot.at("count"));
      photo.hotspot.b_center       = hotspot.at("b_center").get<double>();
      photo.hotspot.b_profile      = hotspot.at("b_profile").get<double>();
      photo.hotspot.strength_sigma = hotspot.at("strength_sigma").get<double>();
      photo.hotspot.Validate();
      if (param.config.count == 0) {
        throw std::invalid_argument("ReadUPCSteering: hotspot configuration counts must be positive");
      }
    }
  }
  std::string mass_data;
  for (const auto& table : param.mass.tables) { mass_data += gra::aux::ReadFile(table.file); }
  param.table_fingerprint = std::to_string(gra::aux::djb2hash(block.dump() + physics.dump() + numerics.dump() + mass_data));
  return param;
}

}  // namespace gra::nuclear
