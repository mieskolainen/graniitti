// Tensor Pomeron process initialization
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Tensor/MTensorInit.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "Graniitti/MGlobals.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tensor/MTensorInit.h"
#include "Graniitti/Tensor/MTensorPomeron.h"

using gra::aux::indices;

namespace gra {

// Compute the beam restriction for the configured Tensor photo continuum
std::string TensorPhotoBeamError(nuclear::CollisionType collision, const nlohmann::json &general, const MPDG &pdg) {
  if (collision != nuclear::CollisionType::PA && collision != nuclear::CollisionType::EA &&
      collision != nuclear::CollisionType::AA) { return {}; }
  const auto &photo = general.at("PARAM_TENSORPOM").at("PHOTO");
  if (photo.at("photon_exchange").get<bool>()) {
    return "Nuclear TP[PHOTO] requires PARAM_TENSORPOM.PHOTO.photon_exchange=false";
  }
  for (const int exchange : photo.at("exchanges").get<std::vector<int>>()) {
    if (pdg.FindByPDG(exchange).isospinX2 != 0) {
      return "Nuclear TP[PHOTO] requires isoscalar exchanges in PARAM_TENSORPOM.PHOTO.exchanges";
    }
  }
  return {};
}

namespace {

// Validate one raw Tensor resonance production channel
void ValidateTensorChannel(const PARAM_RES &resonance,
                           const RES_TENSOR_CHANNEL &channel,
                           const MTensorPomeronParam &parameters,
                           const std::string &name) {
  const bool first_photon = channel.exchange[0] == PDG::PDG_gamma;
  const bool second_photon = channel.exchange[1] == PDG::PDG_gamma;
  const bool vector_resonance = resonance.p.spinX2 == 2 &&
                                resonance.p.P == -1 && resonance.p.C == -1 &&
                                resonance.p.chargeX3 == 0;
  if (!std::all_of(
          channel.g_tensor.cbegin(), channel.g_tensor.cend(),
          [](const double coupling) { return std::isfinite(coupling); })) {
    throw std::invalid_argument("Tensor resonance " + name +
                                " has a nonfinite production coupling");
  }
  for (const auto &mixing : channel.VMD_MIXING) {
    if (!std::all_of(
            mixing.coupling.cbegin(), mixing.coupling.cend(),
            [](const double coupling) { return std::isfinite(coupling); })) {
      throw std::invalid_argument("Tensor resonance " + name +
                                  " has a nonfinite VMD mixing coupling");
    }
    (void)parameters.FindVMD(mixing.pdg);
  }
  const bool spin3 = resonance.p.spinX2 == 6 && resonance.p.P == -1 && resonance.p.C == -1 && resonance.p.chargeX3 == 0;
  if (spin3 && ((!first_photon && !second_photon) || !channel.VMD_MIXING.empty())) {
    throw std::invalid_argument("Tensor spin-three production requires one photon without VMD mixing");
  }
  if (vector_resonance || spin3) {
    if (first_photon && second_photon) {
      throw std::invalid_argument("Tensor resonance " + name +
                                  " cannot use a two-photon TP channel");
    }
    if (first_photon || second_photon) {
      const int strong_pdg =
          first_photon ? channel.exchange[1] : channel.exchange[0];
      if (parameters.exchange.FindExchange(strong_pdg).rank != 2) {
        throw std::invalid_argument(
            "Tensor resonance " + name +
            " photon fusion requires one rank-two exchange");
      }
      return;
    }
    if (!channel.VMD_MIXING.empty()) {
      throw std::invalid_argument("Tensor resonance " + name +
                                  " has VMD mixing without a photon");
    }
    const int first_rank =
        parameters.exchange.FindExchange(channel.exchange[0]).rank;
    const int second_rank =
        parameters.exchange.FindExchange(channel.exchange[1]).rank;
    if (first_rank == second_rank) {
      throw std::invalid_argument(
          "Tensor resonance " + name +
          " strong vector production requires rank-one and rank-two "
          "exchange fusion");
    }
    return;
  }
  if (first_photon || second_photon || !channel.VMD_MIXING.empty() ||
      parameters.exchange.FindExchange(channel.exchange[0]).rank != 2 ||
      parameters.exchange.FindExchange(channel.exchange[1]).rank != 2) {
    throw std::invalid_argument("Tensor resonance " + name +
                                " requires two configured rank-two exchanges");
  }
}

// Build compact active Tensor resonance coupling indices
bool PrepareTensorChannel(RES_TENSOR_CHANNEL &channel,
                          const double coupling_min) {
  channel.active_g_tensor.clear();
  for (const auto &i : indices(channel.g_tensor)) {
    if (std::abs(channel.g_tensor[i]) > coupling_min) {
      channel.active_g_tensor.push_back(i);
    }
  }
  channel.active_vmd.clear();
  for (const auto &i : indices(channel.VMD_MIXING)) {
    auto &mixing = channel.VMD_MIXING[i];
    mixing.active_coupling.clear();
    for (const auto &j : indices(mixing.coupling)) {
      if (std::abs(mixing.coupling[j]) > coupling_min) {
        mixing.active_coupling.push_back(j);
      }
    }
    if (!mixing.active_coupling.empty()) {
      channel.active_vmd.push_back(i);
    }
  }
  channel.active_ready = true;
  return !channel.active_g_tensor.empty() || !channel.active_vmd.empty();
}

// Compute whether one direct continuum exchange ordering is active
bool ActiveDirectTensorPair(const MTensorExchangeModel &exchange,
                            const std::array<int, 2> &pair,
                            const int first_hadron, const int second_hadron) {
  const auto &first = exchange.FindVertex(pair[0], first_hadron);
  const auto &second = exchange.FindVertex(pair[1], second_hadron);
  return !first.active_g_tensor.empty() && !second.active_g_tensor.empty();
}

// Compute whether one vector exchange pair has an active internal line
bool ActiveVectorTensorPair(const MTensorExchangeModel &exchange, const std::array<int, 2> &pair,
                            const int vector_pdg) {
  for (const auto &ordered : MTensorExchangeModel::OrderedPairs(pair)) {
    if (!exchange
             .FindActiveTransfers(ordered[0], ordered[1], vector_pdg,
                                  vector_pdg)
             .empty()) {
      return true;
    }
  }
  return false;
}

// Validate the elementary amplitudes available for lepton and nuclear beams
void ValidateTensorPhoto(const MProcessSetup &setup, MTensorPomeronMode mode,
                               const MTensorPomeronParam &parameters) {
  const auto &lts = setup.lts;
  const bool ion = nuclear::IsNuclearPDG(lts.beam1.pdg) || nuclear::IsNuclearPDG(lts.beam2.pdg);
  const bool lepton = nuclear::IsChargedLepton(lts.beam1.pdg) || nuclear::IsChargedLepton(lts.beam2.pdg);
  if (!ion && !lepton && mode != MTensorPomeronMode::Photo) { return; }
  if (mode != MTensorPomeronMode::Resonance && mode != MTensorPomeronMode::Photo) {
    throw std::invalid_argument("Tensor photoproduction requires RES or PHOTO");
  }
  if (setup.excitation != 0 && (ion || mode == MTensorPomeronMode::Photo)) {
    throw std::invalid_argument("Tensor nuclear and continuum photoproduction require elastic forward particles");
  }
  for (const auto &[name, res] : lts.process.RESONANCES) {
    if (ion && lts.upc_model->Param().photo_model != nuclear::PhotoModel::Impulse &&
        std::none_of(res.TP.channels.begin(), res.TP.channels.end(), [](const auto& channel) {
          return !channel.active_g_tensor.empty();
        })) {
      throw std::invalid_argument("Tensor nuclear shadowing requires a diagonal vector-nucleon amplitude");
    }
    for (const auto &channel : res.TP.channels) {
      if (channel.active_g_tensor.empty() && channel.active_vmd.empty()) { continue; }
      // PHOTO uses the same elementary vector channels also in proton collisions
      const auto collision = nuclear::ClassifyCollision(lts.beam1.pdg, lts.beam2.pdg);
      const auto error = PhotoProductionError(ReggeProductionModel::TP, collision, res.p, channel.exchange,
                                            lts.PDG, setup.model_tune->General());
      if (!error.empty()) { throw std::invalid_argument("Tensor photoproduction: " + name + ": " + error); }
      const int exchange = channel.exchange[channel.exchange[0] == PDG::PDG_gamma ? 1 : 0];
      if (ion) {
        // Strong VMD transitions obey the neutral-component isospin selection rule
        const auto allowed = [&](int incoming) {
          return !math::IsZero(wigner::CG(0.5 * lts.PDG.FindByPDG(incoming).isospinX2,
              0.5 * lts.PDG.FindByPDG(exchange).isospinX2, 0.0, 0.0,
              0.5 * lts.PDG.FindByPDG(res.p.pdg).isospinX2, 0.0));
        };
        if ((!channel.active_g_tensor.empty() && !allowed(res.p.pdg)) ||
            std::any_of(channel.active_vmd.begin(), channel.active_vmd.end(), [&](std::size_t i) {
              return !allowed(channel.VMD_MIXING[i].pdg);
            })) {
          throw std::invalid_argument("Tensor nuclear VMD production violates strong isospin conservation");
        }
      }
      if (ion && lts.PDG.FindByPDG(exchange).isospinX2 == 2) {
        const auto& upc = lts.upc_model->Param();
        if (upc.photo_model == nuclear::PhotoModel::LTA) {
          throw std::invalid_argument("Isovector photoproduction requires impulse or Glauber, not gluon LTA shadowing");
        }
        for (const auto leg : indices(upc.target)) {
          if (lts.upc_model->Type(leg + 1) == nuclear::BeamType::Nucleus &&
              upc.target[leg] != nuclear::CoherenceType::Coherent && upc.structure == nuclear::StructureType::Smooth) {
            throw std::invalid_argument("Isovector incoherent photoproduction requires sampled nucleon or hotspot structure");
          }
        }
      }
    }
  }
  if (ion && mode == MTensorPomeronMode::Photo) {
    if (!lts.upc_model || lts.upc_model->Param().photo_model != nuclear::PhotoModel::Impulse) {
      throw std::invalid_argument("Nuclear TP[PHOTO] requires photoproduction.target_model=impulse");
    }
    const auto error = TensorPhotoBeamError(nuclear::ClassifyCollision(lts.beam1.pdg, lts.beam2.pdg),
                                            setup.model_tune->General(), lts.PDG);
    if (!error.empty()) { throw std::invalid_argument(error); }
  }
  if (mode != MTensorPomeronMode::Photo) { return; }
  const int meson_pdg = std::abs(lts.decaytree[0].p.pdg);
  if (!parameters.photo.lambda.contains(meson_pdg)) {
    throw std::invalid_argument("TP[PHOTO] requires a charged-meson scale under PARAM_TENSORPOM.PHOTO.Lambda");
  }
  // Resonance channels have already passed their active-coupling validation
  bool active = parameters.photo.photon_exchange || !setup.lts.process.RESONANCES.empty();
  for (const int pdg : parameters.photo.exchanges) {
    const auto &meson = parameters.exchange.FindVertex(pdg, meson_pdg);
    const auto &proton = parameters.exchange.FindVertex(pdg, PDG::PDG_p);
    if (!meson.active_g_tensor.empty() && !proton.active_g_tensor.empty()) {
      active = true;
      (void)parameters.exchange.SoftId(pdg, *setup.soft_model);
    }
  }
  if (!active) {
    throw std::invalid_argument(
        "TP[PHOTO] production model is entirely zero");
  }
}

} // namespace

// Resolve tune-specific Tensor Pomeron final-state process_record data once
void SetupTensorProcessModel(MProcessSetup &setup) {
  setup.lts.process.TENSOR_MODEL_READY = false;
  setup.lts.process.TENSOR_PSEUDOSCALAR_PDGS.clear();
  setup.lts.process.TENSOR_BARYON_PDGS.clear();
  setup.lts.process.TENSOR_VECTOR_DECAY_PDGS.clear();
  if (setup.istate != "TP") {
    return;
  }

  const auto param =
      GetTensorParam(*setup.lts.model_cache, setup.lts.PDG, setup.lts.process.RESONANCES);
  for (const auto &entry : param->pseudoscalars) {
    setup.lts.process.TENSOR_PSEUDOSCALAR_PDGS.push_back(std::abs(entry.pdg));
  }
  for (const auto &entry : param->baryons) {
    setup.lts.process.TENSOR_BARYON_PDGS.push_back(std::abs(entry.pdg));
  }
  for (const auto &entry : param->vectors) {
    setup.lts.process.TENSOR_VECTOR_DECAY_PDGS.emplace(
        std::abs(entry.pdg), std::abs(entry.decay_daughter_pdg));
  }
  setup.lts.process.TENSOR_MODEL_READY = true;
}

// Validate Tensor resonance exchange ranks before event sampling
void ValidateTensorResonanceChannels(MProcessSetup &setup) {
  const auto parameters =
      GetTensorParam(*setup.lts.model_cache, setup.lts.PDG, setup.lts.process.RESONANCES);
  const double coupling_min = setup.model_tune->Global().coupling_min;
  for (auto &[name, resonance] : setup.lts.process.RESONANCES) {
    for (const auto &channel : resonance.TP.channels) {
      ValidateTensorChannel(resonance, channel, *parameters, name);
    }
    bool active = false;
    for (auto &channel : resonance.TP.channels) {
      active = PrepareTensorChannel(channel, coupling_min) || active;
    }
    if (!active) {
      throw std::invalid_argument("Tensor resonance " + name +
                                  " production model is entirely zero");
    }
  }
}

// Validate every Tensor exchange used by a forward leg against the SOFT model
void ValidateTensorSoftSources(const MProcessSetup &setup,
                               const BranchingProcessFlags &flags) {
  const auto parameters =
      GetTensorParam(*setup.lts.model_cache, setup.lts.PDG, setup.lts.process.RESONANCES);
  const auto validate = [&](const int exchange_pdg) {
    if (exchange_pdg != PDG::PDG_gamma) {
      (void)parameters->exchange.SoftId(exchange_pdg, *setup.soft_model);
    }
  };

  if (flags.resonance) {
    for (const auto &[name, resonance] : setup.lts.process.RESONANCES) {
      (void)name;
      for (const auto &channel : resonance.TP.channels) {
        validate(channel.exchange[0]);
        validate(channel.exchange[1]);
      }
    }
  }
  if (flags.continuum) {
    if (setup.lts.decaytree.size() != 2) {
      throw std::invalid_argument(
          "Tensor continuum requires two top-level final states");
    }
    const auto &pairs = parameters->exchange.FindContinuumPairs(
        setup.lts.decaytree[0].p.pdg, setup.lts.decaytree[1].p.pdg);
    const bool vector_pair = setup.lts.decaytree[0].p.spinX2 == 2;
    bool active = false;
    for (const auto &pair : pairs) {
      validate(pair[0]);
      validate(pair[1]);
      if (vector_pair) {
        active = ActiveVectorTensorPair(parameters->exchange, pair, setup.lts.decaytree[0].p.pdg) || active;
        continue;
      }
      for (const auto &ordered : MTensorExchangeModel::OrderedPairs(pair)) {
        active = ActiveDirectTensorPair(parameters->exchange, ordered,
                                        setup.lts.decaytree[0].p.pdg,
                                        setup.lts.decaytree[1].p.pdg) ||
                 active;
      }
    }
    if (!active) {
      throw std::invalid_argument(
          "Tensor continuum production model is entirely zero");
    }
  }
}

// Validate HERA only for vector gamma-Pomeron target transitions with measured rows
void ValidateTensorDissociation(const MProcessSetup &setup, MTensorPomeronMode mode,
                               const MTensorPomeronParam &parameters) {
  if (setup.excitation == 0) { return; }
  const auto &state = setup.lts.process;
  if ((mode == MTensorPomeronMode::Continuum || mode == MTensorPomeronMode::ResonanceContinuum) &&
      state.DISSOCIATION == DissociationType::Hera) {
    throw std::invalid_argument("PARAM_NSTAR.MODEL: hera has no TP hadronic continuum profile");
  }
  if (mode == MTensorPomeronMode::Photo && state.PHOTO_DISSOCIATION == DissociationType::Hera) {
    throw std::invalid_argument("PARAM_NSTAR.MODEL: TP[PHOTO] continuum has no HERA vector-meson profile");
  }
  for (const auto &[name, res] : state.RESONANCES) {
    for (const auto &channel : res.TP.channels) {
      if (channel.active_g_tensor.empty() && channel.active_vmd.empty()) { continue; }
      const bool first = channel.exchange[0] == PDG::PDG_gamma;
      const bool second = channel.exchange[1] == PDG::PDG_gamma;
      const bool photo = first != second &&
          parameters.exchange.FindExchange(channel.exchange[first ? 1 : 0]).type == TensorExchangeType::Pomeron;
      const auto type = photo ? state.PHOTO_DISSOCIATION : state.DISSOCIATION;
      if (type != DissociationType::Hera) { continue; }
      if (res.p.spinX2 != 2 || first == second ||
          parameters.exchange.FindExchange(channel.exchange[first ? 1 : 0]).type != TensorExchangeType::Pomeron) {
        throw std::invalid_argument("PARAM_NSTAR.MODEL: TP hera requires vector gamma-Pomeron production: " + name);
      }
      if (!parameters.photo_diss.contains(res.p.pdg)) {
        throw std::invalid_argument("PARAM_NSTAR.MODEL: photoprod_diss has no TP row for " + name);
      }
    }
  }
}

// Initialize tensor Pomeron resonance and continuum branching structures
void MTensorPomeron::InitializeBranching(MProcessSetup &setup,
                                         MTensorPomeronMode mode) {
  const BranchingProcessFlags flags = ClassifyBranchingProcess(setup);
  if (flags.model != ReggeProductionModel::TP) {
    throw std::invalid_argument(
        "MTensorPomeron::InitializeBranching: production model mismatch");
  }
  InitializePhysicalProcessState(setup, flags);
  ValidateTensorResonanceChannels(setup);
  ValidateTensorSoftSources(setup, flags);
  const auto parameters = GetTensorParam(*setup.lts.model_cache, setup.lts.PDG, setup.lts.process.RESONANCES);
  ValidateTensorDissociation(setup, mode, *parameters);
  ValidateTensorPhoto(setup, mode, *parameters);
}

} // namespace gra
