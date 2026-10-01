// Tensor Pomeron parameters
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorParam.h"

#include <cmath>
#include <iostream>
#include <set>
#include <stdexcept>
#include <utility>

#include "Graniitti/Tech/MAux.h"

#include "json.hpp"

using gra::aux::indices;

namespace gra {
namespace {

using json = nlohmann::json;

// Parse one integer sign without narrowing oversized JSON values
int ParseUnitSign(const json &value, const std::string &path) {
  if (!value.is_number_integer()) {
    throw std::invalid_argument(path + " must be -1 or 1");
  }
  if (value.is_number_unsigned()) {
    if (value.get<json::number_unsigned_t>() == 1U) {
      return 1;
    }
  } else {
    const auto sign = value.get<json::number_integer_t>();
    if (sign == -1 || sign == 1) {
      return static_cast<int>(sign);
    }
  }
  throw std::invalid_argument(path + " must be -1 or 1");
}

} // namespace

// Parse one exchanged-vector Regge phase prescription
TensorVectorReggePhase
MTensorPomeronParam::ParseVectorReggePhase(const std::string &mode) {
  if (mode == "NONE") {
    return TensorVectorReggePhase::None;
  }
  if (mode == "EXPONENTIAL") {
    return TensorVectorReggePhase::Exponential;
  }
  throw std::invalid_argument(
      "PARAM_TENSORPOM.VECTOR.phase contains unsupported mode '" + mode + "'");
}

// Parse one transverse vector-meson width prescription
TensorVectorWidthModel
MTensorPomeronParam::ParseVectorWidthModel(const std::string &mode) {
  if (mode == "P_WAVE_2BODY") {
    return TensorVectorWidthModel::PWaveTwoBody;
  }
  if (mode == "RHO_OMEGA") { return TensorVectorWidthModel::RhoOmega; }
  if (mode == "CONSTANT") {
    return TensorVectorWidthModel::Constant;
  }
  throw std::invalid_argument(
      "PARAM_TENSORPOM.VECTOR.Wmode contains unsupported mode '" + mode + "'");
}

// Find pseudoscalar parameters by absolute PDG id
const MTensorPomeronParam::PseudoscalarParam &
MTensorPomeronParam::FindPseudoscalar(const int pdg) const {
  for (const auto &entry : pseudoscalars) {
    if (entry.pdg == std::abs(pdg)) {
      return entry;
    }
  }
  throw std::invalid_argument(
      "MTensorPomeronParam::FindPseudoscalar: Unsupported pseudoscalar pdg = " +
      std::to_string(pdg));
}

// Find baryon parameters by absolute PDG id
const MTensorPomeronParam::BaryonParam &
MTensorPomeronParam::FindBaryon(const int pdg) const {
  for (const auto &entry : baryons) {
    if (entry.pdg == std::abs(pdg)) {
      return entry;
    }
  }
  throw std::invalid_argument(
      "MTensorPomeronParam::FindBaryon: Unsupported baryon pdg = " +
      std::to_string(pdg));
}

// Find vector parameters by absolute PDG id
const MTensorPomeronParam::VectorParam &
MTensorPomeronParam::FindVector(const int pdg) const {
  for (const auto &entry : vectors) {
    if (entry.pdg == std::abs(pdg)) {
      return entry;
    }
  }
  throw std::invalid_argument(
      "MTensorPomeronParam::FindVector: Unsupported vector meson pdg = " +
      std::to_string(pdg));
}

// Find VMD parameters by PDG id
const MTensorPomeronParam::VMDParam &
MTensorPomeronParam::FindVMD(const int pdg) const {
  for (const auto &vmd : VMD) {
    if (vmd.pdg == pdg) {
      return vmd;
    }
  }
  throw std::invalid_argument(
      "MTensorPomeronParam::FindVMD: Unknown vector meson pdg = " +
      std::to_string(pdg));
}

// Read parameters from immutable model cards and one PDG snapshot
void MTensorPomeronParam::Configure(const nlohmann::json &j,
                                    const nlohmann::json &continuum,
                                    const std::string &source,
                                    const MPDG &pdg_table,
                                    const double coupling_min) {
  using json = nlohmann::json;
  try {
    const std::string XID = "PARAM_TENSORPOM";

    FORWARD_NOFLIP = j.at(XID).at("FORWARD_NOFLIP").get<bool>();
    use_zeta       = j.at(XID).at("use_zeta").get<bool>();
    photo_diss     = flux::ReadPhotoDiss(j.at("PARAM_REGGE").at("photoprod_diss"));
    exchange.Configure(j, continuum, pdg_table, coupling_min);

    const json &photo_card = j.at(XID).at("PHOTO");
    const std::set<std::string> photo_fields = {"exchanges", "photon_exchange",
                                                "M0_2",      "Lambda", "offshell",
                                                "q2_max",    "proton_pauli"};
    if (!photo_card.is_object() || photo_card.size() != photo_fields.size()) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM.PHOTO has invalid fields");
    }
    for (const auto &[field, value] : photo_card.items()) {
      (void)value;
      if (!photo_fields.contains(field)) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.PHOTO has unknown field " + field);
      }
    }
    photo.proton_pauli = photo_card.at("proton_pauli").get<bool>();
    photo.photon_exchange = photo_card.at("photon_exchange").get<bool>();
    photo.m0_2 = photo_card.at("M0_2").get<double>();
    photo.q2_max = photo_card.at("q2_max").get<double>();
    photo.exchanges = photo_card.at("exchanges").get<std::vector<int>>();
    if (!(photo.m0_2 > 0.0) ||
        !(photo.q2_max >= 0.0) || photo.exchanges.empty() ||
        !std::isfinite(photo.m0_2) ||
        !std::isfinite(photo.q2_max)) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM.PHOTO has invalid numerical values");
    }
    photo.lambda.clear();
    photo.offshell.clear();
    const auto &offshell = photo_card.at("offshell");
    if (!offshell.is_object() || offshell.size() != photo_card.at("Lambda").size()) {
      throw std::invalid_argument("PARAM_TENSORPOM.PHOTO.offshell requires the charged-meson species in Lambda");
    }
    for (const auto &[key, value] : photo_card.at("Lambda").items()) {
      const int pdg = std::stoi(key);
      const double lambda = value.get<double>();
      if ((pdg != PDG::PDG_pip && pdg != PDG::PDG_Kp) || !(lambda > 0.0) || !std::isfinite(lambda)) {
        throw std::invalid_argument("PARAM_TENSORPOM.PHOTO.Lambda requires positive charged-meson scales");
      }
      photo.lambda.emplace(pdg, lambda);
      const auto &vertex = offshell.at(key);
      const auto ff = regge::ReadFF(vertex.at("FF_offshell"), "PARAM_TENSORPOM.PHOTO.offshell.FF_offshell");
      if (vertex.size() != 2 || (ff.type != regge::FFType::None &&
          (ff.norm != regge::FFNorm::Pole || (ff.type != regge::FFType::Exponential && ff.type != regge::FFType::Gaussian)))) {
        throw std::invalid_argument("PARAM_TENSORPOM.PHOTO.offshell requires FF_transfer and pole-normalized none, exp or gaussian FF_offshell");
      }
      photo.offshell.emplace(pdg, std::make_pair(regge::ReadFF(vertex.at("FF_transfer"),
          "PARAM_TENSORPOM.PHOTO.offshell.FF_transfer"), ff));
    }
    std::set<int> photo_seen;
    for (const int pdg : photo.exchanges) {
      const auto &photo_exchange = exchange.FindExchange(pdg);
      if (!photo_seen.insert(pdg).second ||
          (photo_exchange.type != TensorExchangeType::Pomeron &&
           photo_exchange.type != TensorExchangeType::TensorReggeon &&
           photo_exchange.type != TensorExchangeType::VectorReggeon)) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.PHOTO.exchanges has an invalid entry");
      }
      (void)exchange.FindVertex(pdg, PDG::PDG_p);
    }

    pseudoscalars.clear();
    for (const int pdg : exchange.HadronPdgs(0)) {
      const auto &vertex = exchange.FindVertex(995, pdg);
      pseudoscalars.push_back({pdg, vertex.g_tensor[0], vertex.ff_offshell});
    }
    if (pseudoscalars.empty()) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM requires a pseudoscalar continuum channel");
    }
    meson_ff_transfer =
        exchange.FindVertex(995, pseudoscalars.front().pdg).ff_transfer;

    baryons.clear();
    for (const int pdg : exchange.HadronPdgs(1)) {
      const auto &vertex = exchange.FindVertex(995, pdg);
      baryons.push_back({pdg, vertex.g_tensor[0], vertex.ff_offshell});
    }

    const auto &mixing = j.at(XID).at("RHO_OMEGA");
    rho_omega.mass = {pdg_table.FindByPDG(113).mass, pdg_table.FindByPDG(223).mass};
    rho_omega.width = pdg_table.FindByPDG(223).width;
    rho_omega.pion = pdg_table.FindByPDG(211).mass;
    rho_omega.kaon = pdg_table.FindByPDG(321).mass;
    if (!mixing.at("g").is_array() || mixing.at("g").size() != 2) {
      throw std::invalid_argument("PARAM_TENSORPOM.RHO_OMEGA.g requires rho and omega pion couplings");
    }
    rho_omega.g = mixing.at("g").get<std::array<double, 2>>();
    rho_omega.b = mixing.at("b").get<double>();
    if (mixing.size() != 2 || !std::isfinite(rho_omega.b) ||
        !std::isfinite(rho_omega.g[0]) || !std::isfinite(rho_omega.g[1]) ||
        !(rho_omega.pion > 0.0) || !(rho_omega.kaon > 0.0) ||
        !std::isfinite(rho_omega.pion) || !std::isfinite(rho_omega.kaon) ||
        !std::isfinite(rho_omega.width) || !(rho_omega.width > 0.0)) {
      throw std::invalid_argument("PARAM_TENSORPOM.RHO_OMEGA has invalid pion couplings or pole data");
    }
    for (const double mass : rho_omega.mass) {
      if (!std::isfinite(mass) || !(mass > 2.0 * rho_omega.pion)) {
        throw std::invalid_argument("PARAM_TENSORPOM.RHO_OMEGA requires poles above the two-pion threshold");
      }
    }

    const json &vector = j.at(XID).at("VECTOR");
    const auto vector_pdg = vector.at("PDG").get<std::vector<int>>();
    const auto vector_a0 = vector.at("a0").get<std::vector<double>>();
    const auto vector_ap = vector.at("ap").get<std::vector<double>>();
    const auto vector_regge_phase =
        vector.at("phase").get<std::vector<std::string>>();
    const auto vector_width_model =
        vector.at("Wmode").get<std::vector<std::string>>();
    const auto vector_daughter = vector.at("dPDG").get<std::vector<int>>();
    const std::size_t nvector = vector_pdg.size();
    if (nvector == 0 || vector_a0.size() != nvector ||
        vector_ap.size() != nvector || vector_regge_phase.size() != nvector ||
        vector_width_model.size() != nvector ||
        vector_daughter.size() != nvector) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM.VECTOR arrays must have the same nonzero size");
    }
    vectors.clear();
    std::set<int> vector_seen;
    for (const auto &i : indices(vector_pdg)) {
      const auto &particle = pdg_table.FindByPDG(vector_pdg[i]);
      if (!vector_seen.insert(vector_pdg[i]).second || particle.spinX2 != 2 ||
          particle.P != -1 || particle.C != -1 || particle.chargeX3 != 0) {
        throw std::invalid_argument("PARAM_TENSORPOM.VECTOR requires distinct neutral vector mesons");
      }
      const auto &daughter = pdg_table.FindByPDG(vector_daughter[i]);
      if (daughter.spinX2 != 0 || daughter.C != 0) {
        throw std::invalid_argument("PARAM_TENSORPOM.VECTOR.dPDG requires a spin-zero particle-antiparticle pair");
      }
      // Require a distinct antiparticle, excluding K_S and K_L
      (void)pdg_table.FindByPDG(-daughter.pdg);
      const double daughter_mass = daughter.mass;
      const TensorVectorReggePhase regge_phase =
          ParseVectorReggePhase(vector_regge_phase[i]);
      const TensorVectorWidthModel width_model =
          ParseVectorWidthModel(vector_width_model[i]);
      if (width_model == TensorVectorWidthModel::RhoOmega &&
          ((vector_pdg[i] != 113 && vector_pdg[i] != 223) || vector_daughter[i] != 211)) {
        throw std::invalid_argument("RHO_OMEGA propagation requires rho or omega with charged-pion decay");
      }
      const auto &vertex = exchange.FindVertex(995, vector_pdg[i]);
      if (vector_pdg[i] <= 0 || !std::isfinite(vector_a0[i]) ||
          !std::isfinite(vector_ap[i]) || !std::isfinite(daughter_mass) ||
          !std::isfinite(particle.mass) || !std::isfinite(particle.width) ||
          !(vector_ap[i] > 0.0) || vector_daughter[i] <= 0 ||
          !(daughter_mass > 0.0) || !(particle.mass > 0.0) ||
          !(particle.width > 0.0)) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.VECTOR contains invalid values");
      }
      if (width_model == TensorVectorWidthModel::PWaveTwoBody &&
          !(particle.mass > 2.0 * daughter_mass)) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.VECTOR P-wave pole must be above its two-body threshold");
      }
      vectors.push_back({vector_pdg[i],
                         {vertex.g_tensor[0], vertex.g_tensor[1]},
                         vertex.ff_transfer,
                         vector_a0[i],
                         vector_ap[i],
                         vertex.ff_offshell,
                         regge_phase,
                         width_model,
                         vector_daughter[i],
                         daughter_mass,
                         particle.mass,
                         particle.width});
    }

    const json &vmd = j.at(XID).at("VMD");
    const auto vmd_pdg = vmd.at("PDG").get<std::vector<int>>();
    const auto vmd_gammaV2 = vmd.at("gammaV2").get<std::vector<double>>();
    const json &vmd_gammaV_sign_values = vmd.at("gammaV_sign");
    if (!vmd_gammaV_sign_values.is_array()) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM.VMD.gammaV_sign must contain integers");
    }
    std::vector<int> vmd_gammaV_sign;
    vmd_gammaV_sign.reserve(vmd_gammaV_sign_values.size());
    for (const auto &value : vmd_gammaV_sign_values) {
      vmd_gammaV_sign.push_back(
          ParseUnitSign(value, "PARAM_TENSORPOM.VMD.gammaV_sign"));
    }
    if (vmd_pdg.size() != vmd_gammaV2.size() ||
        vmd_pdg.size() != vmd_gammaV_sign.size()) {
      throw std::invalid_argument(
          "PARAM_TENSORPOM.VMD arrays must have the same size");
    }

    VMD.clear();
    std::set<int> vmd_seen;
    for (const auto &i : indices(vmd_pdg)) {
      const auto &particle = pdg_table.FindByPDG(vmd_pdg[i]);
      if (!vmd_seen.insert(vmd_pdg[i]).second || particle.spinX2 != 2 ||
          particle.P != -1 || particle.C != -1 || particle.chargeX3 != 0) {
        throw std::invalid_argument("PARAM_TENSORPOM.VMD requires distinct neutral vector mesons");
      }
      const double mass = particle.mass;
      if (vmd_pdg[i] == 0 || !std::isfinite(vmd_gammaV2[i]) ||
          vmd_gammaV2[i] <= 0.0 || !std::isfinite(mass) || mass <= 0.0) {
        throw std::invalid_argument(
            "PARAM_TENSORPOM.VMD contains invalid values");
      }
      VMD.push_back({vmd_pdg[i], vmd_gammaV2[i], vmd_gammaV_sign[i], mass});
    }

    std::cout << "MTensorPomeron::ReadParameters: [PARAM_TENSORPOM]"
              << std::endl;
    std::cout << j.at(XID) << std::endl << std::endl;
    initialized = true;
  } catch (const std::exception &e) {
    throw std::invalid_argument(
        "MTensorPomeronParam::Configure: Error parsing " + source + " (" +
        e.what() + ")");
  }
}

// Read one immutable Tensor Pomeron parameter block from a tune
MTensorPomeronParamPtr ReadTensorPomeronParam(const MModelTune &tune,
                                              const MPDG &pdg_table) {
  auto param = std::make_shared<MTensorPomeronParam>();
  param->Configure(tune.General(), tune.Continuum("TP"), tune.GeneralFile(), pdg_table,
                   tune.Global().coupling_min);
  param->photo_dt = tune.Numerics("NUMERICS_REGGE").at("photo_dt").get<double>();
  if (!(param->photo_dt > 0.0) || !std::isfinite(param->photo_dt)) {
    throw std::invalid_argument("NUMERICS_REGGE.photo_dt must be finite and positive");
  }
  return param;
}

// Compute the run owned immutable Tensor Pomeron parameter block
MTensorPomeronParamPtr GetTensorParam(MModelCache &cache,
                                      const MPDG &pdg_table, const std::map<std::string, PARAM_RES> &resonances) {
  return cache.Get<MTensorPomeronParam>("tensor", [&cache, &pdg_table, &resonances] {
    auto poles = pdg_table;
    for (const auto &[name, res] : resonances) {
      (void)name;
      if (res.p.spinX2 != 2 || res.p.P != -1 || res.p.C != -1) { continue; }
      auto &pole = poles.PDG_table.at(res.p.pdg);
      pole.mass = res.p.mass;
      pole.width = res.p.width;
    }
    return ReadTensorPomeronParam(cache.Tune(), poles);
  });
}

} // namespace gra
