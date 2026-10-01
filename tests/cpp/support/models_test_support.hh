// Shared model test fixtures and builders
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#pragma once

#if defined(__GNUC__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wunused-function"
#endif

#include <algorithm>
#include <array>
#include <atomic>
#include <catch.hpp>
#include <cctype>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <functional>
#include <iomanip>
#include <iostream>
#include <limits>
#include <map>
#include <memory>
#include <optional>
#include <random>
#include <sstream>
#include <thread>
#include <type_traits>
#include <utility>

#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_JJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_W/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_Z/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Parton/MG5_PP_ZJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_JJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_WW/Processes.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/ProcessBase.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Eikonal/MEikonal.h"
#include "Graniitti/Eikonal/MProtonScreen.h"
#include "Graniitti/Kinematics/M4Vec.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MUserCuts.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MForm.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MGamma.h"
#include "Graniitti/Photon/MPhotoQCD.h"
#include "Graniitti/Photon/MPhotoVM.h"
#include "Graniitti/Photon/MPhotoZ.h"
#include "Graniitti/Photon/MQED.h"
#include "Graniitti/Process/MHelicityConfig.h"
#include "Graniitti/Process/MProcess.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/QCD/MDurham.h"
#include "Graniitti/Regge/MFragment.h"
#include "Graniitti/Regge/MRegge.h"
#include "Graniitti/Regge/MReggeMPXP.h"
#include "Graniitti/Regge/MReggeForm.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Regge/MReggeGPInit.h"
#include "Graniitti/Regge/MReggeMP.h"
#include "Graniitti/Regge/MReggeXP.h"
#include "Graniitti/Sampling/MRandom.h"
#include "Graniitti/Spin/MDirac.h"
#include "Graniitti/Spin/MHelicity.h"
#include "Graniitti/Spin/MHelicityDecay.h"
#include "Graniitti/Spin/MHelicityNorm.h"
#include "Graniitti/Spin/MHelicityScatter.h"
#include "Graniitti/Spin/MPoleLS.h"
#include "Graniitti/Spin/MWigner.h"
#include "Graniitti/Tech/MLHE.h"
#include "Graniitti/Tensor/MTensorPomeron.h"
#include "HepMC3/Attribute.h"
#include "json.hpp"

using namespace gra;

using gra::aux::indices;
using gra::math::pow2;

const std::string modelfile = (std::filesystem::path(__FILE__).parent_path().parent_path().parent_path().parent_path() /
                               "modeldata/TUNE0/GENERAL.json")
                                  .string();
const std::string pdgfile = gra::MPDG::DataFile();

namespace {

// Set one toy Tensor channel using the production pair fixed by its JPC
gra::RES_TENSOR_CHANNEL &SetToyTensorChannel(gra::PARAM_RES &res, const std::vector<double> &couplings) {
  gra::RES_TENSOR_CHANNEL channel;
  const bool              vector = res.p.spinX2 == 2 && res.p.P == -1 && res.p.C == -1;
  channel.exchange               = vector ? std::array<int, 2>{22, 995} : std::array<int, 2>{995, 995};
  channel.g_tensor               = couplings;
  res.TP.channels                = {std::move(channel)};
  return res.TP.channels.front();
}

// Compute the sole Tensor channel used by one focused model test
const gra::RES_TENSOR_CHANNEL &ToyTensorChannel(const gra::PARAM_RES &res) {
  REQUIRE(res.TP.channels.size() == 1);
  return res.TP.channels.front();
}

// Evaluate one resonance through the public Regge amplitude interface
double TestReggeRes(gra::MRegge &regge, gra::LORENTZSCALAR &lts, gra::PARAM_RES &res, gra::ReggeProductionModel spin) {
  gra::PARAM_RES input = res;
  auto           saved = std::move(lts.process.RESONANCES);
  lts.process.RESONANCES.emplace("__test_resonance__", std::move(input));
  try {
    const double amp2      = regge.Amp2(lts, spin, gra::MReggeMode::Resonance);
    res                    = lts.process.RESONANCES.begin()->second;
    lts.process.RESONANCES = std::move(saved);
    return amp2;
  } catch (...) {
    res                    = lts.process.RESONANCES.begin()->second;
    lts.process.RESONANCES = std::move(saved);
    throw;
  }
}

// Evaluate one continuum through the public Regge amplitude interface
double TestReggeCon(gra::MRegge &regge, gra::LORENTZSCALAR &lts, gra::ReggeProductionModel spin) {
  if (lts.process.CONT_TU_SIGN.empty()) { gra::MRegge::InitializeContinuumInterference(lts); }
  return regge.Amp2(lts, spin, gra::MReggeMode::ContinuumTwoBody);
}

// Sum finite-spin resonance production channels for matrix-level tests
void TestProd3Sum(const gra::LORENTZSCALAR &lts, gra::PARAM_RES &res, double s0 = 1.0) {
  auto channels = gra::rspin::Resonance(
      lts, res,
      (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), s0);
  res.prod_f    = gra::MMatrix<std::complex<double>>();
  for (auto &channel : channels) {
    if (res.prod_f.size_row() == 0) {
      res.prod_f = std::move(channel);
    } else {
      res.prod_f += channel;
    }
  }
}

// Sum finite-spin continuum production channels for matrix-level tests
std::pair<gra::MMatrix<std::complex<double>>, gra::MMatrix<std::complex<double>>> TestProd4Sum(
    const gra::LORENTZSCALAR &lts, double s0 = 1.0) {
  auto                               channels = gra::mpom::Continuum(lts, s0);
  gra::MMatrix<std::complex<double>> sum_t;
  gra::MMatrix<std::complex<double>> sum_u;
  for (std::size_t i = 0; i < channels.size(); ++i) {
    if (i == 0) {
      sum_t = std::move(channels[i].first);
      sum_u = std::move(channels[i].second);
    } else {
      sum_t += channels[i].first;
      sum_u += channels[i].second;
    }
  }
  return {std::move(sum_t), std::move(sum_u)};
}

std::vector<gra::MParticle> ProtonInitialState();

// Load an independent particle table without changing another worker's tune
gra::MPDG LoadedPDGTable() {
  gra::MPDG pdg;
  pdg.ReadParticleData(pdgfile, "TUNE0");
  return pdg;
}

// Write one temporary Durham tune with complete GENERAL and NUMERICS snapshots
std::pair<std::string, std::string> WriteModifiedDurhamTune(
    const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate_general) {
  const std::filesystem::path dir = "tmp/graniitti_durham_" + suffix;
  std::filesystem::create_directories(dir);

  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  if (mutate_general) { mutate_general(general); }
  const std::filesystem::path general_path = dir / "GENERAL.json";
  std::ofstream               general_out(general_path);
  if (!general_out.good()) { throw std::runtime_error("WriteModifiedDurhamTune: failed to write GENERAL.json"); }
  general_out << general.dump(2);
  general_out.close();

  const std::filesystem::path numerics_path = dir / "NUMERICS.json";
  std::ofstream               numerics_out(numerics_path);
  if (!numerics_out.good()) { throw std::runtime_error("WriteModifiedDurhamTune: failed to write NUMERICS.json"); }
  numerics_out << gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"));
  numerics_out.close();

  // Preserve the complete continuum snapshot required by MModelTune
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    const std::string card = "CON_" + model + ".json";
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", card), dir / card,
                               std::filesystem::copy_options::overwrite_existing);
  }

  return {dir.string(), general_path.string()};
}

// Write one temporary Durham tune with relaxed parton level jet cuts
std::pair<std::string, std::string> WriteRelaxedDurhamTune(const std::string &suffix) {
  return WriteModifiedDurhamTune(suffix, [](auto &general) {
    general["PARAM_DURHAM"]["JET_pt_min"]  = 0.0;
    general["PARAM_DURHAM"]["JET_rap_max"] = 100.0;
  });
}

// Write one temporary PhotoZ tune with modified physical and numerical steering
std::pair<std::string, std::string> WriteModifiedPhotoZTune(
    const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate_general,
    const std::function<void(nlohmann::json &)> &mutate_numerics) {
  const std::filesystem::path dir = "tmp/graniitti_photoz_" + suffix;
  std::filesystem::create_directories(dir);

  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  if (mutate_general) { mutate_general(general); }
  const std::filesystem::path general_path = dir / "GENERAL.json";
  std::ofstream               general_out(general_path);
  if (!general_out.good()) { throw std::runtime_error("WriteModifiedPhotoZTune: failed to write GENERAL.json"); }
  general_out << general.dump(2);
  general_out.close();

  auto numerics = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  mutate_numerics(numerics);

  const std::filesystem::path numerics_path = dir / "NUMERICS.json";
  std::ofstream               numerics_out(numerics_path);
  if (!numerics_out.good()) { throw std::runtime_error("WriteModifiedPhotoZTune: failed to write NUMERICS.json"); }
  numerics_out << numerics.dump(2);
  numerics_out.close();

  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    const std::string filename = "CON_" + model + ".json";
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", filename), dir / filename,
                               std::filesystem::copy_options::overwrite_existing);
  }
  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "DECAYS.json"), dir / "DECAYS.json",
                             std::filesystem::copy_options::overwrite_existing);

  return {dir.string(), general_path.string()};
}

// Set one explicit continuum field for every exchange of a selected hadron pair
void SetContinuumField(nlohmann::json &card, const std::string &pair, const std::string &field, const nlohmann::json &value) {
  for (auto &[exchange, pairs] : card.items()) {
    (void)exchange;
    if (pairs.contains(pair)) { pairs.at(pair)[field] = value; }
  }
}

// Configure an explicit unit photon-pion vertex for synthetic continuum channels
void SetToyPhotonContinuum(nlohmann::json &card) {
  const auto &pion = card.begin().value().at("[211,211]");
  const nlohmann::json photon = {
      {"opposite", {{"basis", "crossed_ls"}, {"CP", {true, true}}, {"g_ls", {{1, 0, 1.0, 0.0}}}}},
      {"FF_transfer", {{"type", "none"}}},
      {"FF_offshell", {{"type", "none"}}},
      {"reggeize", pion.at("reggeize")},
      {"pveto", pion.at("pveto")}};
  card["22"]["[211,211]"] = photon;
}

// Write one temporary PhotoVM tune with modified physical steering
std::pair<std::string, std::string> WriteModifiedPhotoVMTune(
    const std::string &suffix, const std::function<void(nlohmann::json &)> &mutate_general,
    const std::function<void(nlohmann::json &)> &mutate_tensor    = {},
    const std::function<void(nlohmann::json &)> &mutate_continuum = {}) {
  const auto dir = std::filesystem::path(gra::aux::GetBasePath(2)) / "tmp" / ("graniitti_photovm_" + suffix);
  std::filesystem::create_directories(dir);

  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  if (mutate_general) { mutate_general(general); }
  const std::filesystem::path general_path = dir / "GENERAL.json";
  std::ofstream               general_out(general_path);
  if (!general_out.good()) { throw std::runtime_error("WriteModifiedPhotoVMTune: failed to write GENERAL.json"); }
  general_out << general.dump(2);
  general_out.close();

  const std::filesystem::path numerics_path = dir / "NUMERICS.json";
  std::ofstream               numerics_out(numerics_path);
  if (!numerics_out.good()) { throw std::runtime_error("WriteModifiedPhotoVMTune: failed to write NUMERICS.json"); }
  numerics_out << gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"));
  numerics_out.close();

  auto tensor = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "CON_TP.json")));
  if (mutate_tensor) { mutate_tensor(tensor); }
  const std::filesystem::path tensor_path = dir / "CON_TP.json";
  std::ofstream               tensor_out(tensor_path);
  if (!tensor_out.good()) { throw std::runtime_error("WriteModifiedPhotoVMTune: failed to write CON_TP.json"); }
  tensor_out << tensor.dump(2);
  tensor_out.close();

  for (const std::string model : {"MP", "XP", "GP"}) {
    const std::string filename = "CON_" + model + ".json";
    auto continuum = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", filename)));
    if (mutate_continuum) { mutate_continuum(continuum); }
    std::ofstream continuum_out(dir / filename);
    if (!continuum_out.good()) { throw std::runtime_error("WriteModifiedPhotoVMTune: failed to write " + filename); }
    continuum_out << continuum.dump(2);
    continuum_out.close();
  }
  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "DECAYS.json"), dir / "DECAYS.json",
                             std::filesystem::copy_options::overwrite_existing);

  return {dir.string(), general_path.string()};
}

// Write one temporary tune with one modified central continuum card
//
std::string WriteModifiedContinuumTune(const std::string &suffix, const std::string &model,
                                       const std::function<void(nlohmann::json &)> &mutate_continuum) {
  const std::filesystem::path dir = "tmp/graniitti_continuum_" + suffix;
  std::filesystem::create_directories(dir);

  for (const std::string card_model : {"MP", "XP", "GP", "TP"}) {
    const std::string filename = "CON_" + card_model + ".json";
    auto continuum = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", filename)));
    if (card_model == model) { mutate_continuum(continuum); }
    const std::filesystem::path continuum_path = dir / filename;
    std::ofstream               continuum_out(continuum_path);
    if (!continuum_out.good()) { throw std::runtime_error("WriteModifiedContinuumTune: failed to write " + filename); }
    continuum_out << continuum.dump(2);
    continuum_out.close();
  }

  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "DECAYS.json"), dir / "DECAYS.json",
                             std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "GENERAL.json"), dir / "GENERAL.json",
                             std::filesystem::copy_options::overwrite_existing);
  std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"), dir / "NUMERICS.json",
                             std::filesystem::copy_options::overwrite_existing);
  return dir.string();
}

// Write one temporary tune with independently modified skewed-UGD and numerical
// steering
std::string WriteModifiedSudakovModelTune(const std::string                           &suffix,
                                          const std::function<void(nlohmann::json &)> &mutate_general,
                                          const std::function<void(nlohmann::json &)> &mutate_numerics) {
  const std::filesystem::path dir = "tmp/graniitti_sudakov_model_" + suffix;
  std::filesystem::create_directories(dir);

  auto general = nlohmann::json::parse(gra::aux::GetInputData(modelfile));
  if (mutate_general) { mutate_general(general); }
  const std::filesystem::path general_path = dir / "GENERAL.json";
  std::ofstream               general_out(general_path);
  if (!general_out.good()) { throw std::runtime_error("WriteModifiedSudakovModelTune: failed to write GENERAL.json"); }
  general_out << general.dump(2);
  general_out.close();

  auto numerics = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  mutate_numerics(numerics);
  const std::filesystem::path numerics_path = dir / "NUMERICS.json";
  std::ofstream               numerics_out(numerics_path);
  if (!numerics_out.good()) {
    throw std::runtime_error("WriteModifiedSudakovModelTune: failed to write NUMERICS.json");
  }
  numerics_out << numerics.dump(2);
  numerics_out.close();

  // Preserve the complete continuum snapshot required by MModelTune
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    const std::string card = "CON_" + model + ".json";
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", card), dir / card,
                               std::filesystem::copy_options::overwrite_existing);
  }

  return dir.string();
}

// Construct a nonphoton forward branch from prepared helicity metadata
gra::MDecayBranch ForwardBranchForTest(const gra::HELMatrix &hel) {
  gra::MDecayBranch branch;
  branch.hel = hel;
  return branch;
}

// Boost the four-vector into the rest frame of the system
gra::M4Vec BoostToRestFrame(const gra::M4Vec &p, const gra::M4Vec &system) {
  gra::M4Vec out = p;
  REQUIRE(system.M2() > 0.0);
  gra::kinematics::LorentzBoost(system, system.M(), out, -1);
  return out;
}

// Compute the same energy with the spatial direction reversed
gra::M4Vec OppositeSpatialDirectionForTest(const gra::M4Vec &p) {
  gra::M4Vec out = p;
  out.Flip3();
  return out;
}

// Compute the Jacob-Wick second-leg crossing phase used by the production helper
double CrossedSecondLegPhaseForTest(double spin, double physical_lambda) {
  const double exponent = spin - physical_lambda;
  const double rounded  = std::round(exponent);
  REQUIRE(std::abs(exponent - rounded) < 1e-9);
  const long long n = std::llround(rounded);
  return (n % 2 == 0) ? 1.0 : -1.0;
}

// Find the explicit helicity row matching one lambda pair
std::size_t FindHelicityPairRowForTest(const gra::HELMatrix &hel, double lambda1, double lambda2) {
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) {
    if (std::abs(hel.lambda_values[i][0] - lambda1) < 1e-9 && std::abs(hel.lambda_values[i][1] - lambda2) < 1e-9) {
      return i;
    }
  }
  FAIL("Could not find helicity pair row");
  return 0;
}

// Convert an auxiliary crossed matrix into the physical second-leg helicity
// basis
MMatrix<std::complex<double>> CrossedSecondRowsToPhysicalForTest(const gra::HELMatrix                &hel,
                                                                 const MMatrix<std::complex<double>> &aux) {
  REQUIRE(aux.size_row() == hel.lambda_values.size_row());
  const double                  s2 = (hel.s2 >= 0.0) ? hel.s2 : 0.0;
  MMatrix<std::complex<double>> out(aux.size_row(), aux.size_col(), 0.0);

  for (std::size_t row = 0; row < hel.lambda_values.size_row(); ++row) {
    const double      lambda1 = hel.lambda_values[row][0];
    const double      lambda2 = hel.lambda_values[row][1];
    const std::size_t aux_row = FindHelicityPairRowForTest(hel, lambda1, -lambda2);
    const double      phase   = CrossedSecondLegPhaseForTest(s2, lambda2);

    for (std::size_t col = 0; col < aux.size_col(); ++col) { out[row][col] = phase * aux[aux_row][col]; }
  }
  return out;
}

// Apply the ordered upper/lower phase of one virtual Regge subchannel
void ApplyVirtualSubchannelPhaseForTest(MMatrix<std::complex<double>> &f, const gra::HELMatrix &hel,
                                        const gra::M4Vec &axis_in_X, bool second_exchange_daughter) {
  REQUIRE(f.size_col() == hel.Jz_values.size());
  const double phase_sign = second_exchange_daughter ? 1.0 : -1.0;
  for (std::size_t col = 0; col < f.size_col(); ++col) {
    const std::complex<double> phase = std::exp(gra::math::zi * phase_sign * hel.Jz_values[col] * axis_in_X.Phi());
    for (std::size_t row = 0; row < f.size_row(); ++row) { f[row][col] *= phase; }
  }
}

gra::HELMatrix SpinHalfLegHelicityMatrix() {
  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, 1.0, 0.5, 0.5, {-1.0, 0.0, 1.0}, {-0.5, 0.5}, {-0.5, 0.5},
                              "SpinHalfLegHelicityMatrix");
  hel.T       = MMatrix<std::complex<double>>(2, 2, 0.0);
  hel.T[1][0] = 1.0;
  hel.T[0][1] = std::complex<double>(0.75, -0.2);
  return hel;
}

gra::HELMatrix RealisticProtonLegHelicityMatrix(int exchange_pdg, int exchange_spinX2, int exchange_parity,
                                                std::size_t l, std::size_t two_s) {
  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = true;
  hel.alpha_ls.Set(l, two_s, 1.0);

  gra::MParticle exchange;
  exchange.name   = "toy_exchange";
  exchange.pdg    = exchange_pdg;
  exchange.spinX2 = exchange_spinX2;
  exchange.P      = exchange_parity;

  const auto beams = ProtonInitialState();
  gra::spin::InitTMatrix(hel, exchange, beams[0], beams[1], true, "test proton leg", false, false);
  return hel;
}

std::vector<std::size_t> ActiveProtonLegRowsPerIncoming(const gra::HELMatrix &hel) {
  std::vector<std::size_t> counts(hel.T.size_row(), 0);
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) {
    const std::size_t i1 = hel.lambda_idx[i][0];
    const std::size_t i2 = hel.lambda_idx[i][1];
    if (std::abs(hel.T[i1][i2]) > 1e-12) { ++counts[i1]; }
  }
  return counts;
}

gra::HELMatrix PseudoscalarFusionHelicityMatrix(int exchange_pdg, int exchange_spinX2, int exchange_parity,
                                                std::size_t l, std::size_t two_s) {
  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = true;
  hel.alpha_ls.Set(l, two_s, 1.0);

  gra::MParticle resonance;
  resonance.name   = "toy_pseudoscalar";
  resonance.pdg    = 9000221;
  resonance.spinX2 = 0;
  resonance.P      = -1;

  gra::MParticle exchange;
  exchange.name   = "toy_exchange";
  exchange.pdg    = exchange_pdg;
  exchange.spinX2 = exchange_spinX2;
  exchange.P      = exchange_parity;

  gra::spin::InitTMatrix(hel, resonance, exchange, exchange, true, "test pseudoscalar fusion", false, false);
  return hel;
}

gra::MParticle ToyParticle(int pdg, int spinX2, int P, int C, const std::string &name) {
  gra::MParticle p;
  p.pdg    = pdg;
  p.spinX2 = spinX2;
  p.P      = P;
  p.C      = C;
  p.name   = name;
  return p;
}

gra::HELMatrix ToyHelicityMatrix(std::size_t l, std::size_t two_s, bool P_symmetry = true, bool C_symmetry = true) {
  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = P_symmetry;
  hel.C_symmetry = C_symmetry;
  hel.alpha_ls.Set(l, two_s, 1.0);
  return hel;
}

std::vector<std::complex<double>> MHelicityScatterCentralAmplitudes(const gra::HELMatrix &hel, const gra::M4Vec &q1,
                                                                    const gra::M4Vec &q2) {
  gra::LORENTZSCALAR lts;
  lts.pfinal.resize(1);
  lts.pfinal[0] = q1 + q2;

  gra::PARAM_RES res;

  std::vector<gra::M4Vec> pair = {q1, q2};
  gra::spin::ProductionFrame(pair, lts, "CM", lts.pfinal[0]);
  const auto matrix = gra::spin::fDecayMatrix(hel, pair[0].Theta(), pair[0].Phi());
  REQUIRE(matrix.size_col() == 1);

  std::vector<std::complex<double>> out;
  out.reserve(matrix.size_row());
  for (std::size_t i = 0; i < matrix.size_row(); ++i) { out.push_back(matrix[i][0]); }
  return out;
}

std::vector<std::complex<double>> ProjectVectorPseudoscalarVertex(const gra::MTensorPomeron &tensor,
                                                                  const gra::M4Vec &q1, const gra::M4Vec &q2) {
  const auto vertex = tensor.iG_psvv(q1, q2, (q1 + q2).M(), 1.0, gra::regge::ReadFF({{"type", "none"}}, "test decay"));
  return tensor.MassiveSpin1PolSum(vertex, q1, q2);
}

std::vector<std::complex<double>> ProjectTensorPseudoscalarVertex(const gra::MTensorPomeron &tensor,
                                                                  const gra::M4Vec &q1, const gra::M4Vec &q2,
                                                                  int structure) {
  const auto vertex = (structure == 0) ? tensor.iG_PPPS_0(q1, q2, 1.0) : tensor.iG_PPPS_1(q1, q2, 1.0);

  std::vector<std::complex<double>> out;
  out.reserve(25);
  for (int h1 = -2; h1 <= 2; ++h1) {
    const auto eps1 = tensor.EpsMassiveSpin2(q1, h1);
    for (int h2 = -2; h2 <= 2; ++h2) {
      const auto           eps2 = tensor.EpsMassiveSpin2(q2, h2);
      std::complex<double> amp  = 0.0;
      for (std::size_t mu = 0; mu < 4; ++mu) {
        for (std::size_t nu = 0; nu < 4; ++nu) {
          for (std::size_t kappa = 0; kappa < 4; ++kappa) {
            for (std::size_t lambda = 0; lambda < 4; ++lambda) {
              amp += std::conj(eps1(mu, nu)) * std::conj(eps2(kappa, lambda)) * vertex(mu, nu, kappa, lambda);
            }
          }
        }
      }
      out.push_back(amp);
    }
  }
  return out;
}

// Project one covariant PPf1 vertex onto spin-2, spin-2 and spin-1 states
std::vector<std::complex<double>> ProjectTensorAxialVertex(const gra::MTensorPomeron &tensor, const gra::M4Vec &q1,
                                                           const gra::M4Vec &q2, int structure) {
  const auto vertex       = (structure == 0) ? tensor.iG_PPA_22(q1, q2, 1.0) : tensor.iG_PPA_44(q1, q2, 1.0);
  const auto axial_states = tensor.MassiveSpin1States(q1 + q2, "conj", true);

  std::vector<std::complex<double>> out;
  out.reserve(75);
  for (int h1 = -2; h1 <= 2; ++h1) {
    const auto eps1 = tensor.EpsMassiveSpin2(q1, h1);
    for (int h2 = -2; h2 <= 2; ++h2) {
      const auto eps2 = tensor.EpsMassiveSpin2(q2, h2);
      // TopologyMode the Appendix B chi2 phases to the MDirac helicity phases
      const double paper_second_leg_phase = std::abs(h2) == 1 ? -1.0 : 1.0;
      for (std::size_t hX = 0; hX < axial_states.size(); ++hX) {
        std::complex<double> amplitude = 0.0;
        for (std::size_t kappa = 0; kappa < 4; ++kappa) {
          for (std::size_t lambda = 0; lambda < 4; ++lambda) {
            for (std::size_t rho = 0; rho < 4; ++rho) {
              for (std::size_t sigma = 0; sigma < 4; ++sigma) {
                for (std::size_t alpha = 0; alpha < 4; ++alpha) {
                  amplitude += eps1(kappa, lambda) * (paper_second_leg_phase * eps2(rho, sigma)) *
                               vertex({kappa, lambda, rho, sigma, alpha}) * axial_states[hX](alpha);
                }
              }
            }
          }
        }
        out.push_back(amplitude);
      }
    }
  }
  return out;
}

// Compute the reduced PPf1 helicity table in arXiv:2008.07452 Appendix B
std::vector<std::complex<double>> PaperTensorAxialReducedTable(double exchange_mass, double momentum, int structure) {
  const double longitudinal = (pow2(exchange_mass) + 4.0 * pow2(momentum)) / (std::sqrt(6.0) * pow2(exchange_mass));
  std::vector<std::complex<double>> out(75, 0.0);
  auto                              set = [&](int h1, int h2, int hX, double value) {
    const std::size_t index =
        static_cast<std::size_t>(h1 + 2) * 15 + static_cast<std::size_t>(h2 + 2) * 3 + static_cast<std::size_t>(hX + 1);
    out[index] = value;
  };

  if (structure == 0) {
    set(2, 1, 1, -1.0);
    set(1, 0, 1, longitudinal);
    set(0, -1, 1, longitudinal);
    set(-1, -2, 1, -1.0);
    set(1, 2, -1, 1.0);
    set(0, 1, -1, -longitudinal);
    set(-1, 0, -1, -longitudinal);
    set(-2, -1, -1, 1.0);
  } else {
    set(1, 0, 1, 1.0);
    set(0, -1, 1, 1.0);
    set(0, 1, -1, -1.0);
    set(-1, 0, -1, -1.0);
  }
  return out;
}

// Check the covariant PPf1 tensor symmetries, traces and axial transversality
void RequireAxialVertexIdentities(const gra::MTensor<std::complex<double>> &vertex, const gra::M4Vec &central,
                                  double tolerance = 1e-11) {
  for (std::size_t kappa = 0; kappa < 4; ++kappa) {
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
      for (std::size_t rho = 0; rho < 4; ++rho) {
        for (std::size_t sigma = 0; sigma < 4; ++sigma) {
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            const std::vector<std::size_t> index       = {kappa, lambda, rho, sigma, alpha};
            const std::vector<std::size_t> swap_first  = {lambda, kappa, rho, sigma, alpha};
            const std::vector<std::size_t> swap_second = {kappa, lambda, sigma, rho, alpha};
            CHECK(std::abs(vertex(index) - vertex(swap_first)) < tolerance);
            CHECK(std::abs(vertex(index) - vertex(swap_second)) < tolerance);
          }

          std::complex<double> longitudinal = 0.0;
          for (std::size_t alpha = 0; alpha < 4; ++alpha) {
            longitudinal += central[alpha] * vertex({kappa, lambda, rho, sigma, alpha});
          }
          CHECK(std::abs(longitudinal) < tolerance);
        }
      }
    }
  }

  for (std::size_t rho = 0; rho < 4; ++rho) {
    for (std::size_t sigma = 0; sigma < 4; ++sigma) {
      for (std::size_t alpha = 0; alpha < 4; ++alpha) {
        std::complex<double> trace = 0.0;
        for (std::size_t kappa = 0; kappa < 4; ++kappa) {
          const double metric = kappa == 0 ? 1.0 : -1.0;
          trace += metric * vertex({kappa, kappa, rho, sigma, alpha});
        }
        CHECK(std::abs(trace) < tolerance);
      }
    }
  }
  for (std::size_t kappa = 0; kappa < 4; ++kappa) {
    for (std::size_t lambda = 0; lambda < 4; ++lambda) {
      for (std::size_t alpha = 0; alpha < 4; ++alpha) {
        std::complex<double> trace = 0.0;
        for (std::size_t rho = 0; rho < 4; ++rho) {
          const double metric = rho == 0 ? 1.0 : -1.0;
          trace += metric * vertex({kappa, lambda, rho, rho, alpha});
        }
        CHECK(std::abs(trace) < tolerance);
      }
    }
  }
}

void RequireProjectivelyEqual(const std::vector<std::complex<double>> &actual,
                              const std::vector<std::complex<double>> &expected, double tol = 1e-10) {
  REQUIRE(actual.size() == expected.size());

  double actual_norm2   = 0.0;
  double expected_norm2 = 0.0;
  for (std::size_t i = 0; i < actual.size(); ++i) {
    actual_norm2 += std::norm(actual[i]);
    expected_norm2 += std::norm(expected[i]);
  }
  REQUIRE(actual_norm2 > tol);
  REQUIRE(expected_norm2 > tol);

  std::size_t pivot = actual.size();
  for (std::size_t i = 0; i < expected.size(); ++i) {
    if (std::abs(expected[i]) > tol) {
      pivot = i;
      break;
    }
  }
  REQUIRE(pivot < expected.size());

  const std::complex<double> scale = actual[pivot] / expected[pivot];
  for (std::size_t i = 0; i < actual.size(); ++i) {
    CAPTURE(i, actual[i], expected[i], scale);
    const std::complex<double> ref       = scale * expected[i];
    const double               scale_abs = std::max({1.0, std::abs(actual[i]), std::abs(ref)});
    REQUIRE(std::abs(actual[i] - ref) <= tol * scale_abs);
  }
}

std::pair<std::complex<double>, std::complex<double>> TwoVectorSpanCoefficients(
    const std::vector<std::complex<double>> &target, const std::vector<std::complex<double>> &basis0,
    const std::vector<std::complex<double>> &basis1, double tol = 1e-10) {
  REQUIRE(target.size() == basis0.size());
  REQUIRE(target.size() == basis1.size());

  const std::complex<double> g00 = gra::InnerProduct(basis0, basis0);
  const std::complex<double> g01 = gra::InnerProduct(basis0, basis1);
  const std::complex<double> g10 = gra::InnerProduct(basis1, basis0);
  const std::complex<double> g11 = gra::InnerProduct(basis1, basis1);
  const std::complex<double> r0  = gra::InnerProduct(basis0, target);
  const std::complex<double> r1  = gra::InnerProduct(basis1, target);
  const std::complex<double> det = g00 * g11 - g01 * g10;

  REQUIRE(std::abs(det) > tol);

  const std::complex<double> c0 = (r0 * g11 - g01 * r1) / det;
  const std::complex<double> c1 = (g00 * r1 - r0 * g10) / det;
  return {c0, c1};
}

void RequireInTwoVectorSpan(const std::vector<std::complex<double>> &target,
                            const std::vector<std::complex<double>> &basis0,
                            const std::vector<std::complex<double>> &basis1, double tol = 1e-10) {
  const auto [c0, c1] = TwoVectorSpanCoefficients(target, basis0, basis1, tol);

  double residual2 = 0.0;
  double target2   = 0.0;
  for (std::size_t i = 0; i < target.size(); ++i) {
    const std::complex<double> residual = target[i] - c0 * basis0[i] - c1 * basis1[i];
    residual2 += std::norm(residual);
    target2 += std::norm(target[i]);
  }

  CAPTURE(c0, c1, residual2, target2);
  REQUIRE(target2 > tol);
  REQUIRE(residual2 <= tol * target2);
}

gra::HELMatrix SpinZeroCentralResonanceHelicityMatrix() {
  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, 0.0, 0.0, 0.0, {0.0}, {0.0}, {0.0}, "SpinZeroCentralResonanceHelicityMatrix");
  hel.T = MMatrix<std::complex<double>>(1, 1, 1.0);
  return hel;
}

gra::HELMatrix SpinOneCentralResonanceHelicityMatrix() {
  gra::HELMatrix            hel;
  const std::vector<double> spin_one = {-1.0, 0.0, 1.0};
  gra::spin::InitTwoBodyBasis(hel, 1.0, 1.0, 1.0, spin_one, spin_one, spin_one,
                              "SpinOneCentralResonanceHelicityMatrix");
  hel.T       = MMatrix<std::complex<double>>(3, 3, 0.0);
  hel.T[2][1] = 1.0;
  return hel;
}

gra::HELMatrix SpinZeroContinuumSubchannelHelicityMatrix() {
  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, 0.0, 0.0, 0.0, {0.0}, {0.0}, {0.0}, "SpinZeroContinuumSubchannelHelicityMatrix");
  hel.T = MMatrix<std::complex<double>>(1, 1, 1.0);
  return hel;
}

gra::HELMatrix SpinOneContinuumSubchannelHelicityMatrix(std::size_t i_final, std::size_t i_exchange) {
  gra::HELMatrix            hel;
  const std::vector<double> spin_one = {-1.0, 0.0, 1.0};
  gra::spin::InitTwoBodyBasis(hel, 1.0, 1.0, 1.0, spin_one, spin_one, spin_one,
                              "SpinOneContinuumSubchannelHelicityMatrix");
  hel.T                      = MMatrix<std::complex<double>>(3, 3, 0.0);
  hel.T[i_final][i_exchange] = 1.0;
  return hel;
}

std::vector<double> SpinProjectionValues(int spinX2) {
  std::vector<double> out;
  for (int i = 0; i <= spinX2; ++i) { out.push_back(-0.5 * spinX2 + i); }
  return out;
}

gra::HELMatrix ToyContinuumSubchannelHelicityMatrix(int final_spinX2, int other_spinX2, int exchange_spinX2, int seed) {
  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, exchange_spinX2 / 2.0, final_spinX2 / 2.0, other_spinX2 / 2.0,
                              SpinProjectionValues(exchange_spinX2), SpinProjectionValues(final_spinX2),
                              SpinProjectionValues(other_spinX2), "ToyContinuumSubchannelHelicityMatrix");

  const std::size_t n1 = static_cast<std::size_t>(final_spinX2 + 1);
  const std::size_t n2 = static_cast<std::size_t>(other_spinX2 + 1);
  hel.T                = MMatrix<std::complex<double>>(n1, n2, 0.0);
  for (std::size_t i = 0; i < n1; ++i) {
    for (std::size_t j = 0; j < n2; ++j) {
      hel.T[i][j] =
          std::complex<double>(0.17 * (seed + 1) + 0.11 * i - 0.07 * j, 0.05 * (seed + 2) + 0.03 * (i + 1) * (j + 1));
    }
  }
  return hel;
}

// Build a compact analytic trajectory subchannel matrix for GP tests
gra::HELMatrix ToyAnalyticTrajectorySubchannelHelicityMatrix(int final_spinX2, int other_spinX2, int MMAX, int seed) {
  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = false;
  hel.C_symmetry = false;
  gra::gpom::InitHelicity(hel, final_spinX2 / 2.0, other_spinX2 / 2.0, MMAX,
                          "ToyAnalyticTrajectorySubchannelHelicityMatrix");

  if (MMAX == 0 && final_spinX2 == 0 && other_spinX2 == 0) {
    hel.T[0][0]     = 1.0;
    hel.T_set[0][0] = true;
    return hel;
  }

  for (std::size_t row = 0; row < hel.T.size_row(); ++row) {
    for (std::size_t col = 0; col < hel.T.size_col(); ++col) {
      hel.T[row][col] =
          std::complex<double>(0.13 * (seed + 1) + 0.07 * static_cast<double>(row) - 0.03 * static_cast<double>(col),
                               0.02 * (seed + 2) + 0.01 * static_cast<double>((row + 1) * (col + 1)));
      hel.T_set[row][col] = true;
    }
  }
  gra::gpom::CheckHelicity(hel, "ToyAnalyticTrajectorySubchannelHelicityMatrix", false);
  return hel;
}

gra::HELMatrix ToyProtonLegHelicityMatrix(int exchange_spinX2) {
  gra::HELMatrix hel = SpinHalfLegHelicityMatrix();
  hel.J              = exchange_spinX2 / 2.0;
  hel.Jz_values      = SpinProjectionValues(exchange_spinX2);
  return hel;
}

// Generate an on-shell pair at a reproducible rest-frame orientation
void SetCentralPair(gra::LORENTZSCALAR &lts) {
  REQUIRE(lts.decaytree.size() == 2);
  gra::MRandom random;
  random.SetSeed(173);
  std::vector<gra::M4Vec> daughters;
  const auto weight = gra::kinematics::TwoBodyPhaseSpace(
      lts.pfinal[0], lts.pfinal[0].M(), {lts.decaytree[0].p.mass, lts.decaytree[1].p.mass}, daughters, random);
  REQUIRE(weight.Integral() > 0.0);
  for (const auto &i : indices(lts.decaytree)) { lts.decaytree[i].p4 = daughters[i]; }
  REQUIRE(gra::math::CheckEMC(lts.pfinal[0] - daughters[0] - daughters[1]));
  for (const auto &i : indices(daughters)) {
    REQUIRE(daughters[i].M2() == Approx(pow2(lts.decaytree[i].p.mass)).epsilon(1e-10).margin(1e-12));
  }
}

// Build proton production kinematics with asymmetric longitudinal recoil
gra::LORENTZSCALAR MakeToyProductionLTSAsymmetric(double proton_pt, double pz1, double pz2) {
  gra::LORENTZSCALAR lts;
  lts.PDG              = LoadedPDGTable();
  const auto beams     = ProtonInitialState();
  lts.beam1            = beams[0];
  lts.beam2            = beams[1];
  lts.process.MP_FRAME = "HX";
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, 5.0, std::hypot(5.0, lts.beam1.mass));
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -5.0, std::hypot(5.0, lts.beam2.mass));
  lts.pfinal.resize(3);
  lts.pfinal[1].SetPxPyPzM(proton_pt, 0.0, pz1, lts.beam1.mass);
  lts.pfinal[2].SetPxPyPzM(-proton_pt, 0.0, pz2, lts.beam2.mass);
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1        = lts.pbeam1 - lts.pfinal[1];
  lts.q2        = lts.pbeam2 - lts.pfinal[2];
  lts.t1        = lts.q1.M2();
  lts.t2        = lts.q2.M2();
  lts.q1_in_X   = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X   = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  return lts;
}

// Build a symmetric forward-proton event with a controlled transverse angle
gra::LORENTZSCALAR MakeToyProductionLTSDphi(double dphi, double phi0 = 0.0) {
  const double beam_energy   = 100.0;
  const double proton_pt     = 0.2;
  const double proton_pz     = 98.0;
  const double proton_energy = std::sqrt(proton_pz * proton_pz + proton_pt * proton_pt + gra::PDG::mp * gra::PDG::mp);

  gra::LORENTZSCALAR lts;
  lts.PDG              = LoadedPDGTable();
  const auto beams     = ProtonInitialState();
  lts.beam1            = beams[0];
  lts.beam2            = beams[1];
  lts.process.MP_FRAME = "HX";
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, beam_energy, std::hypot(beam_energy, lts.beam1.mass));
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -beam_energy, std::hypot(beam_energy, lts.beam2.mass));
  lts.pfinal.resize(3);
  lts.pfinal[1] = gra::M4Vec(proton_pt * std::cos(phi0), proton_pt * std::sin(phi0), proton_pz, proton_energy);
  lts.pfinal[2] =
      gra::M4Vec(proton_pt * std::cos(phi0 + dphi), proton_pt * std::sin(phi0 + dphi), -proton_pz, proton_energy);
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1        = lts.pbeam1 - lts.pfinal[1];
  lts.q2        = lts.pbeam2 - lts.pfinal[2];
  lts.t1        = lts.q1.M2();
  lts.t2        = lts.q2.M2();
  lts.q1_in_X   = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X   = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  lts.s         = (lts.pbeam1 + lts.pbeam2).M2();
  lts.sqrt_s    = std::sqrt(std::max(0.0, lts.s));
  lts.m2        = lts.pfinal[0].M2();
  lts.s1        = (lts.pfinal[0] + lts.pfinal[1]).M2();
  lts.s2        = (lts.pfinal[0] + lts.pfinal[2]).M2();
  lts.qt1       = lts.q1.Pt();
  lts.qt2       = lts.q2.Pt();
  return lts;
}

gra::PARAM_RES MakeToyResonance() {
  gra::PARAM_RES res;

  gra::MParticle exchange;
  exchange.name   = "toy_exchange";
  exchange.pdg    = 993;
  exchange.spinX2 = 2;
  exchange.P      = -1;
  exchange.C      = 1;

  gra::MParticle mother;
  mother.name   = "toy_resonance";
  mother.pdg    = 9000001;
  mother.spinX2 = 2;
  mother.P      = 1;
  res.p         = mother;

  gra::MDecayBranch up;
  up.p   = exchange;
  up.hel = SpinHalfLegHelicityMatrix();

  gra::MDecayBranch dn = up;
  res.production       = {{{up, dn}, SpinOneCentralResonanceHelicityMatrix()}};
  return res;
}

// Build canonical pole caches for one fixed-spin toy continuum channel
void PrepareToyContinuumOperators(gra::LORENTZSCALAR &lts, gra::ReggeProductionModel model) {
  REQUIRE((model == gra::ReggeProductionModel::MP || model == gra::ReggeProductionModel::XP));
  REQUIRE_FALSE(lts.process.CONT_PRODUCTIONTREE.empty());
  REQUIRE(lts.decaytree.size() == 2);
  const std::array<std::vector<gra::MParticle>, 4> legs = {
      std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p},
      std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
      std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
      std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p}};
  std::vector<std::vector<gra::spin::PoleResidue>> all_cache;
  all_cache.reserve(lts.process.CONT_PRODUCTIONTREE.size());
  for (const auto &tree : lts.process.CONT_PRODUCTIONTREE) {
    REQUIRE(tree.size() == 2);
    const std::array<const gra::MParticle *, 4> mothers = {&tree[0].p, &tree[1].p, &tree[0].p, &tree[1].p};
    std::vector<gra::spin::PoleResidue>         cache;
    cache.reserve(4);
    for (const auto &i : indices(legs)) {
      const auto operators = gra::spin::CanonicalPoleOperators(*mothers[i], legs[i][0], legs[i][1], true, false, false,
                                                               gra::spin::VertexContext::SubTUChannelExchange);
      REQUIRE_FALSE(operators.empty());
      const auto &coupling = operators.front().coupling;
      auto vertex = gra::spin::PreparePoleLS(*mothers[i], legs[i][0], legs[i][1], {{coupling.l, coupling.two_s, 1.0}},
                                             1.0, true, false, false, gra::spin::VertexContext::SubTUChannelExchange);
      cache.emplace_back(vertex);
    }
    all_cache.push_back(std::move(cache));
  }
  lts.process.CONTINUUM_POLE = std::move(all_cache);
}

gra::LORENTZSCALAR MakeToyContinuumLTSAsymmetric(double proton_pt, double pz1, double pz2) {
  gra::LORENTZSCALAR lts;
  lts.PDG              = LoadedPDGTable();
  const auto beams     = ProtonInitialState();
  lts.beam1            = beams[0];
  lts.beam2            = beams[1];
  lts.process.MP_FRAME = "HX";
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, 5.0, std::hypot(5.0, lts.beam1.mass));
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -5.0, std::hypot(5.0, lts.beam2.mass));
  lts.pfinal.resize(3);
  lts.pfinal[1].SetPxPyPzM(proton_pt, 0.1, pz1, lts.beam1.mass);
  lts.pfinal[2].SetPxPyPzM(-0.6 * proton_pt, 0.04, pz2, lts.beam2.mass);
  lts.pfinal[0]       = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1              = lts.pbeam1 - lts.pfinal[1];
  lts.q2              = lts.pbeam2 - lts.pfinal[2];
  lts.t1              = lts.q1.M2();
  lts.t2              = lts.q2.M2();
  lts.q1_in_X         = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X         = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  lts.s_hat           = lts.pfinal[0].M2();
  lts.process.SPINGEN = true;
  lts.decaytree.resize(2);

  gra::MParticle vec;
  vec.name   = "toy_vector";
  vec.pdg    = 113;
  vec.spinX2 = 2;
  vec.P      = -1;
  vec.C      = -1;
  vec.mass   = 0.1;

  lts.decaytree[0].p = vec;
  lts.decaytree[1].p = vec;

  SetCentralPair(lts);

  gra::MParticle exchange;
  exchange.name   = "toy_exchange";
  exchange.pdg    = 993;
  exchange.spinX2 = 2;
  exchange.P      = -1;
  exchange.C      = 1;

  gra::MDecayBranch up;
  up.p                            = exchange;
  up.hel                          = SpinHalfLegHelicityMatrix();
  gra::MDecayBranch dn            = up;
  lts.process.CONT_PRODUCTIONTREE = {{up, dn}};
  lts.process.CONT_PRODUCTION     = {{exchange.pdg, exchange.pdg}};
  PrepareToyContinuumOperators(lts, gra::ReggeProductionModel::MP);
  return lts;
}

gra::LORENTZSCALAR MakeToyScalarContinuumLTSAsymmetric(double proton_pt, double pz1, double pz2) {
  gra::LORENTZSCALAR lts;
  lts.PDG              = LoadedPDGTable();
  const auto beams     = ProtonInitialState();
  lts.beam1            = beams[0];
  lts.beam2            = beams[1];
  lts.process.MP_FRAME = "HX";
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, 5.0, std::hypot(5.0, lts.beam1.mass));
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -5.0, std::hypot(5.0, lts.beam2.mass));
  lts.pfinal.resize(3);
  lts.pfinal[1].SetPxPyPzM(proton_pt, 0.1, pz1, lts.beam1.mass);
  lts.pfinal[2].SetPxPyPzM(-0.6 * proton_pt, 0.04, pz2, lts.beam2.mass);
  lts.pfinal[0]       = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1              = lts.pbeam1 - lts.pfinal[1];
  lts.q2              = lts.pbeam2 - lts.pfinal[2];
  lts.t1              = lts.q1.M2();
  lts.t2              = lts.q2.M2();
  lts.q1_in_X         = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X         = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  lts.s_hat           = lts.pfinal[0].M2();
  lts.process.SPINGEN = true;
  lts.decaytree.resize(2);

  gra::MParticle pion;
  pion.name   = "pi";
  pion.pdg    = 211;
  pion.spinX2 = 0;
  pion.P      = -1;
  pion.mass   = 0.14;

  lts.decaytree[0].p     = pion;
  lts.decaytree[1].p     = pion;
  lts.decaytree[1].p.pdg = -211;

  SetCentralPair(lts);

  gra::MParticle exchange;
  exchange.name   = "toy_scalar_exchange";
  exchange.pdg    = 991;
  exchange.spinX2 = 0;
  exchange.P      = 1;

  gra::MDecayBranch up;
  up.p                            = exchange;
  up.hel                          = RealisticProtonLegHelicityMatrix(991, 0, 1, 0, 0);
  gra::MDecayBranch dn            = up;
  lts.process.CONT_PRODUCTIONTREE = {{up, dn}};
  lts.process.CONT_PRODUCTION     = {{exchange.pdg, exchange.pdg}};
  PrepareToyContinuumOperators(lts, gra::ReggeProductionModel::MP);
  return lts;
}

// Replace one toy continuum cache by four raw canonical XP operators
void PrepareToyXPContinuumOperators(gra::LORENTZSCALAR &lts, const std::vector<gra::spin::LSTerm> &terms) {
  REQUIRE(lts.process.CONT_PRODUCTIONTREE.size() == 1);
  REQUIRE(lts.process.CONT_PRODUCTIONTREE.front().size() == 2);
  REQUIRE(lts.decaytree.size() == 2);

  const auto                                      &tree = lts.process.CONT_PRODUCTIONTREE.front();
  const std::array<std::vector<gra::MParticle>, 4> legs = {
      std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p},
      std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
      std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
      std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p}};
  const std::array<const gra::MParticle *, 4> mothers = {&tree[0].p, &tree[1].p, &tree[0].p, &tree[1].p};

  std::vector<gra::spin::PoleResidue> cache;
  cache.reserve(4);
  for (const auto &i : indices(legs)) {
    auto vertex = gra::spin::PreparePoleLS(*mothers[i], legs[i][0], legs[i][1], terms, 1.0, true, false, true,
                                           gra::spin::VertexContext::SubTUChannelExchange);
    cache.emplace_back(vertex);
  }
  lts.process.CONTINUUM_POLE = {std::move(cache)};
}

// Build one virtual subchannel vertex without a crossed pair rest boost
MMatrix<std::complex<double>> ExplicitVirtualSubchannelVertexForTest(const gra::LORENTZSCALAR &lts,
                                                                     const gra::HELMatrix     &hel,
                                                                     const gra::M4Vec         &final_state,
                                                                     const gra::M4Vec         &axis_in_X,
                                                                     bool second_exchange_daughter) {
  if (hel.UsesReggeDomain()) {
    auto out = CrossedSecondRowsToPhysicalForTest(hel, gra::spin::DirectFrame(hel, "test virtual subchannel"));
    ApplyVirtualSubchannelPhaseForTest(out, hel, axis_in_X, second_exchange_daughter);
    return out;
  }

  gra::M4Vec final_in_X = BoostToRestFrame(final_state, lts.pfinal[0]);
  gra::kinematics::RotateHelicityAxes(final_in_X, axis_in_X);
  auto out =
      CrossedSecondRowsToPhysicalForTest(hel, gra::spin::fDecayMatrix(hel, final_in_X.Theta(), final_in_X.Phi()));
  ApplyVirtualSubchannelPhaseForTest(out, hel, axis_in_X, second_exchange_daughter);
  return out;
}

// Contract two explicit X rest virtual vertices over exchange helicity
MMatrix<std::complex<double>> ManualContinuumSubchannelMatrix(
    const gra::LORENTZSCALAR &lts, const gra::HELMatrix &upper, const gra::HELMatrix &lower,
    const gra::M4Vec &upper_final, const gra::M4Vec &lower_final, const gra::M4Vec &upper_parent_dir,
    const gra::M4Vec &lower_parent_dir, double s_left, double s_right, bool swap_final_order) {
  const auto f_up = ExplicitVirtualSubchannelVertexForTest(lts, upper, upper_final, upper_parent_dir, false);
  const auto f_dn = ExplicitVirtualSubchannelVertexForTest(lts, lower, lower_final, lower_parent_dir, true);

  const std::size_t n_left  = static_cast<std::size_t>(std::round(2.0 * s_left + 1.0));
  const std::size_t n_right = static_cast<std::size_t>(std::round(2.0 * s_right + 1.0));

  MMatrix<std::complex<double>> out(f_up.size_col() * f_dn.size_col(), n_left * n_right, 0.0);
  const double                  tol = 1e-9;

  for (std::size_t i = 0; i < upper.lambda_values.size_row(); ++i) {
    const double lambda_up    = upper.lambda_values[i][0];
    const double lambda_ex_up = upper.lambda_values[i][1];

    for (std::size_t j = 0; j < lower.lambda_values.size_row(); ++j) {
      const double lambda_dn    = lower.lambda_values[j][0];
      const double lambda_ex_dn = lower.lambda_values[j][1];

      if (std::abs(lambda_ex_up + lambda_ex_dn) > tol) { continue; }

      const double      lambda_left  = swap_final_order ? lambda_dn : lambda_up;
      const double      lambda_right = swap_final_order ? lambda_up : lambda_dn;
      const std::size_t i_left       = static_cast<std::size_t>(std::llround(lambda_left + s_left));
      const std::size_t i_right      = static_cast<std::size_t>(std::llround(lambda_right + s_right));
      const std::size_t final_col    = i_left * n_right + i_right;

      for (std::size_t a = 0; a < f_up.size_col(); ++a) {
        for (std::size_t b = 0; b < f_dn.size_col(); ++b) {
          const std::size_t init_row = a * f_dn.size_col() + b;
          out[init_row][final_col] += f_up[i][a] * f_dn[j][b];
        }
      }
    }
  }
  return out;
}

std::vector<std::size_t> HelicityConservingRows(const gra::HELMatrix &hel) {
  const double                          s1         = (hel.s1 >= 0.0) ? hel.s1 : 0.0;
  const std::size_t                     n_incoming = static_cast<std::size_t>(std::llround(2.0 * s1 + 1.0));
  std::vector<std::vector<std::size_t>> rows_by_incoming(n_incoming);

  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) {
    const double lambda_in  = hel.lambda_values[i][0];
    const double lambda_out = hel.lambda_values[i][1];
    if (std::abs(lambda_in - lambda_out) > 1e-9) { continue; }
    const long long incoming = std::llround(lambda_in + s1);
    REQUIRE(incoming >= 0);
    REQUIRE(static_cast<std::size_t>(incoming) < rows_by_incoming.size());
    rows_by_incoming[static_cast<std::size_t>(incoming)].push_back(i);
  }

  std::vector<std::size_t> rows;
  for (const auto &incoming_rows : rows_by_incoming) {
    REQUIRE(incoming_rows.size() == 1);
    rows.push_back(incoming_rows[0]);
  }
  return rows;
}

// Compute all explicit helicity rows
std::vector<std::size_t> AllHelicityRowsForTest(const gra::HELMatrix &hel) {
  std::vector<std::size_t> rows;
  rows.reserve(hel.lambda_values.size_row());
  for (std::size_t i = 0; i < hel.lambda_values.size_row(); ++i) { rows.push_back(i); }
  return rows;
}

MMatrix<std::complex<double>> SelectRows(const MMatrix<std::complex<double>> &input,
                                         const std::vector<std::size_t>      &rows) {
  MMatrix<std::complex<double>> out(rows.size(), input.size_col(), 0.0);
  for (std::size_t i = 0; i < rows.size(); ++i) {
    for (std::size_t j = 0; j < input.size_col(); ++j) { out[i][j] = input[rows[i]][j]; }
  }
  return out;
}

// Compute one selected non-photon source in the elastic matched convention
MMatrix<std::complex<double>> NormalizeForwardRowsForTest(const MMatrix<std::complex<double>> &input) { return input; }

// Build the analytic vertex R_(i,m) before outer product factorization
MMatrix<std::complex<double>> NestedAnalyticReggeResidueForTest(const std::vector<std::pair<double, double>> &rows,
                                                                const std::vector<int>                       &m_values,
                                                                const std::vector<double>                    &nonsense,
                                                                const gra::M4Vec &q, double s0,
                                                                bool second_exchange_daughter,
                                                                bool use_exchange_helicity_barrier) {
  REQUIRE(m_values.size() == nonsense.size());
  MMatrix<std::complex<double>> residue(rows.size(), m_values.size(), 0.0);
  for (std::size_t row = 0; row < rows.size(); ++row) {
    for (std::size_t col = 0; col < m_values.size(); ++col) {
      residue[row][col] = gra::gpom::Residue(m_values[col], 0.0, q.Pt(), q.Phi(), s0, second_exchange_daughter,
                                             use_exchange_helicity_barrier) *
                          nonsense[col];
    }
  }
  return residue;
}

// Build the transverse photon source used on photoproduction branches
MMatrix<std::complex<double>> PhotonLegMatrixForTest(const gra::LORENTZSCALAR &lts, const gra::HELMatrix &hel, int leg,
                                                     bool second_exchange_daughter) {
  std::vector<std::pair<double, double>> transitions;
  transitions.reserve(hel.lambda_values.size_row());
  for (std::size_t row = 0; row < hel.lambda_values.size_row(); ++row) {
    transitions.emplace_back(hel.lambda_values[row][0], hel.lambda_values[row][1]);
  }
  std::vector<int> m_values;
  m_values.reserve(hel.Jz_values.size());
  for (const double m : hel.Jz_values) { m_values.push_back(static_cast<int>(std::llround(m))); }
  return gra::qed::PhotonSourceMatrixTransitions(lts, leg, transitions, m_values, lts.process.PHOTON_VERTEX,
                                                 second_exchange_daughter);
}

std::pair<MMatrix<std::complex<double>>, MMatrix<std::complex<double>>> ExpectedIdenticalContinuumMatrices(
    const gra::LORENTZSCALAR &lts) {
  const auto q1_in_X        = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  const auto q2_in_X        = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  const auto q2_second_axis = OppositeSpatialDirectionForTest(q2_in_X);
  const auto up_rows        = HelicityConservingRows(lts.process.CONT_PRODUCTIONTREE[0][0].hel);
  const auto dn_rows        = HelicityConservingRows(lts.process.CONT_PRODUCTIONTREE[0][1].hel);
  const auto f_up_first     = NormalizeForwardRowsForTest(
          SelectRows(gra::spin::Forward(lts, lts.process.CONT_PRODUCTIONTREE[0][0], lts.pbeam1, lts.pfinal[1], false,
                                        gra::spin::Rows(lts.process.CONT_PRODUCTIONTREE[0][0], false),
                                        lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{}),
                     up_rows));
  const auto f_dn_second = NormalizeForwardRowsForTest(
      SelectRows(gra::spin::Forward(lts, lts.process.CONT_PRODUCTIONTREE[0][1], lts.pbeam2, lts.pfinal[2], true,
                                    gra::spin::Rows(lts.process.CONT_PRODUCTIONTREE[0][1], false),
                                    lts.process.PHOTON_VERTEX, gra::spin::ForwardSpec{}),
                 dn_rows));

  const double s3 = lts.decaytree[0].p.spinX2 / 2.0;
  const double s4 = lts.decaytree[1].p.spinX2 / 2.0;
  auto         t  = f_up_first.Kronecker(f_dn_second) *
           ManualContinuumSubchannelMatrix(lts, lts.process.CONTINUUM_POLE[0][0].Pole().helicity,
                                           lts.process.CONTINUUM_POLE[0][1].Pole().helicity, lts.decaytree[0].p4,
                                           lts.decaytree[1].p4, q1_in_X, q2_second_axis, s3, s4, false);
  auto u = f_up_first.Kronecker(f_dn_second) *
           ManualContinuumSubchannelMatrix(lts, lts.process.CONTINUUM_POLE[0][2].Pole().helicity,
                                           lts.process.CONTINUUM_POLE[0][3].Pole().helicity, lts.decaytree[1].p4,
                                           lts.decaytree[0].p4, q1_in_X, q2_second_axis, s3, s4, true);

  return {t, u};
}

// Compute the squared Frobenius distance between two amplitude matrices
double MatrixDiffNorm2(const MMatrix<std::complex<double>> &a, const MMatrix<std::complex<double>> &b) {
  REQUIRE(a.size_row() == b.size_row());
  REQUIRE(a.size_col() == b.size_col());

  return (a - b).FrobNorm2();
}

// Compute the squared Frobenius norm of one amplitude matrix
double MatrixNorm2(const MMatrix<std::complex<double>> &a) {
  return a.FrobNorm2();
}

// Compute the squared norm of one dynamically ranked tensor
double DynamicTensorNorm2(const gra::MTensor<std::complex<double>> &tensor) {
  if (tensor.empty()) { return 0.0; }
  std::vector<std::size_t> index(tensor.rank(), 0);
  double                   out = 0.0;
  while (true) {
    out += std::norm(tensor(index));
    std::size_t axis = tensor.rank();
    while (axis > 0) {
      --axis;
      if (++index[axis] < tensor.size(axis)) { break; }
      index[axis] = 0;
    }
    if (axis == 0 && index[0] == 0) { break; }
  }
  return out;
}

// Compute the squared norm of one matrix row
double MatrixRowNorm2(const MMatrix<std::complex<double>> &a, std::size_t row) {
  REQUIRE(row < a.size_row());
  return gra::SquaredNorm(a.Row(row));
}

// Convert an exchange spin projection label into an integer m value
int IntegerProjectionForTest(double value) {
  const double rounded = std::round(value);
  REQUIRE(std::abs(value - rounded) < 1e-9);
  return static_cast<int>(std::llround(rounded));
}

// Compute the exchange-helicity scale without modifying the proton vertex
double ReggeResidueScaleForTest(const gra::HELMatrix &hel, std::size_t col, double qt, double s0) {
  REQUIRE(col < hel.Jz_values.size());
  const int    m              = IntegerProjectionForTest(hel.Jz_values[col]);
  const double exchange_power = std::abs(static_cast<double>(m));
  return std::pow(qt / std::sqrt(s0), exchange_power);
}

// Compute the unit exchange-vertex scale without a helicity barrier
double ReggePoleResidueScaleForTest(const gra::HELMatrix &hel, std::size_t col) {
  REQUIRE(col < hel.Jz_values.size());
  return 1.0;
}

// Require finite amplitudes before comparing their real and imaginary components
void RequireComplexNear(const std::complex<double> &actual, const std::complex<double> &expected,
                        double epsilon = 1e-12) {
  const double actual_re = std::real(actual);
  const double actual_im = std::imag(actual);
  const double expect_re = std::real(expected);
  const double expect_im = std::imag(expected);

  REQUIRE(std::isfinite(actual_re));
  REQUIRE(std::isfinite(actual_im));
  REQUIRE(std::isfinite(expect_re));
  REQUIRE(std::isfinite(expect_im));
  REQUIRE(actual_re == Approx(expect_re).epsilon(epsilon).margin(epsilon));
  REQUIRE(actual_im == Approx(expect_im).epsilon(epsilon).margin(epsilon));
}

// Compare finite amplitude matrices with matching dimensions
void RequireMatrixNear(const MMatrix<std::complex<double>> &actual, const MMatrix<std::complex<double>> &expected,
                       double epsilon = 1e-12) {
  REQUIRE(actual.size_row() == expected.size_row());
  REQUIRE(actual.size_col() == expected.size_col());

  for (std::size_t i = 0; i < actual.size_row(); ++i) {
    for (std::size_t j = 0; j < actual.size_col(); ++j) { RequireComplexNear(actual[i][j], expected[i][j], epsilon); }
  }
}

// Compare finite helicity amplitudes in the same ordering
void RequireVectorNear(const std::vector<std::complex<double>> &actual,
                       const std::vector<std::complex<double>> &expected, double epsilon = 1e-12) {
  REQUIRE(actual.size() == expected.size());
  for (const auto &i : indices(actual)) { RequireComplexNear(actual[i], expected[i], epsilon); }
}

// Check the component helicity phases under one beam-axis rotation
void RequireAzimuthalCovariance(const std::vector<std::complex<double>> &rotated,
                                const std::vector<std::complex<double>> &reference, const std::vector<int> &harmonics,
                                double angle, double epsilon = 1e-12) {
  REQUIRE(rotated.size() == reference.size());
  REQUIRE(harmonics.size() == reference.size());
  for (std::size_t i = 0; i < reference.size(); ++i) {
    const std::complex<double> phase = std::exp(gra::math::zi * static_cast<double>(harmonics[i]) * angle);
    RequireComplexNear(rotated[i], phase * reference[i], epsilon);
  }
}

gra::M4Vec BoostFromRestFrame(const gra::M4Vec &p, const gra::M4Vec &system) {
  gra::M4Vec out = p;
  if (system.M() > 0.0) { gra::kinematics::LorentzBoost(system, system.M(), out, 1); }
  return out;
}

// Apply one common Lorentz boost recursively to a decay branch
void BoostDecayBranchForTest(gra::MDecayBranch &branch, const gra::M4Vec &boost, double boost_mass) {
  gra::kinematics::LorentzBoost(boost, boost_mass, branch.p4, 1);
  for (auto &leg : branch.legs) { BoostDecayBranchForTest(leg, boost, boost_mass); }
}

std::array<gra::M4Vec, 2> TwoBodyRestKinematics(double mother_mass, double m1, double m2, double theta, double phi) {
  const double pnorm = gra::kinematics::SqrtKallenLambda(pow2(mother_mass), pow2(m1), pow2(m2)) / (2.0 * mother_mass);
  const double px    = pnorm * std::sin(theta) * std::cos(phi);
  const double py    = pnorm * std::sin(theta) * std::sin(phi);
  const double pz    = pnorm * std::cos(theta);
  const double e1    = std::sqrt(pow2(m1) + pow2(pnorm));
  const double e2    = std::sqrt(pow2(m2) + pow2(pnorm));
  return {gra::M4Vec(px, py, pz, e1), gra::M4Vec(-px, -py, -pz, e2)};
}

gra::MParticle ToyParticle(const std::string &name, int pdg, int spinX2, double mass) {
  gra::MParticle p;
  p.name   = name;
  p.pdg    = pdg;
  p.spinX2 = spinX2;
  p.P      = 1;
  p.mass   = mass;
  return p;
}

gra::HELMatrix SpinOneToSpinOneScalarHelicityMatrix(const std::array<std::complex<double>, 3> &couplings) {
  gra::HELMatrix            hel;
  const std::vector<double> spin_one = {-1.0, 0.0, 1.0};
  gra::spin::InitTwoBodyBasis(hel, 1.0, 1.0, 0.0, spin_one, spin_one, {0.0}, "SpinOneToSpinOneScalarHelicityMatrix");
  hel.T = MMatrix<std::complex<double>>(3, 1, 0.0);

  for (std::size_t i = 0; i < 3; ++i) { hel.T[i][0] = couplings[i]; }
  return hel;
}

gra::HELMatrix SpinOneToScalarScalarHelicityMatrix(std::complex<double> coupling) {
  gra::HELMatrix hel;
  gra::spin::InitTwoBodyBasis(hel, 1.0, 0.0, 0.0, {-1.0, 0.0, 1.0}, {0.0}, {0.0},
                              "SpinOneToScalarScalarHelicityMatrix");
  hel.T       = MMatrix<std::complex<double>>(1, 1, 0.0);
  hel.T[0][0] = coupling;
  return hel;
}

gra::M4Vec ApplyHXStepInCurrentFrame(const gra::M4Vec &p, const gra::M4Vec &axis, const gra::M4Vec &system) {
  auto rotate = [&axis](gra::M4Vec q) {
    q.RotateZ(-axis.Phi());
    q.RotateY(-axis.Theta());
    q.RotateZ(gra::math::PI);
    return q;
  };

  gra::M4Vec out   = rotate(p);
  gra::M4Vec boost = rotate(system);
  if (boost.M() > 0.0) { gra::kinematics::LorentzBoost(boost, boost.M(), out, -1); }
  return out;
}

// Contract the explicit cascade with the first-daughter HX spin phase at each vertex
MMatrix<std::complex<double>> SequentialOneLevelReference(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res) {
  const auto &R = lts.decaytree[0];
  const auto &s = lts.decaytree[1];
  (void)s;

  const auto R_in_X = BoostToRestFrame(R.p4, lts.pfinal[0]);
  const auto fX     = gra::spin::fDecayMatrix(res.hel_decay, R_in_X.Theta(), R_in_X.Phi());

  std::vector<gra::M4Vec> R_daughters_in_X = {BoostToRestFrame(R.legs[0].p4, lts.pfinal[0]),
                                              BoostToRestFrame(R.legs[1].p4, lts.pfinal[0])};
  gra::kinematics::HXframe(R_daughters_in_X, R_in_X);
  const auto fR = gra::spin::fDecayMatrix(R.hel, R_daughters_in_X[0].Theta(), R_daughters_in_X[0].Phi());

  MMatrix<std::complex<double>> expected(3, 1, 0.0);
  for (std::size_t M = 0; M < 3; ++M) {
    for (std::size_t lambdaR = 0; lambdaR < 3; ++lambdaR) {
      const auto phase = std::exp(std::complex<double>(0.0, gra::math::PI * (static_cast<double>(lambdaR) - 1.0)));
      expected[M][0] += fX[lambdaR][M] * fR[0][lambdaR] * phase;
    }
  }
  return expected;
}

// Contract two intermediate spin sums with their explicit HX azimuthal phases
MMatrix<std::complex<double>> SequentialDeepReference(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res) {
  const auto &R = lts.decaytree[0];
  const auto &S = R.legs[0];

  const auto R_in_X = BoostToRestFrame(R.p4, lts.pfinal[0]);
  const auto fX     = gra::spin::fDecayMatrix(res.hel_decay, R_in_X.Theta(), R_in_X.Phi());

  const auto S_in_X        = BoostToRestFrame(S.p4, lts.pfinal[0]);
  const auto d_in_X        = BoostToRestFrame(R.legs[1].p4, lts.pfinal[0]);
  const auto R_system_in_X = S_in_X + d_in_X;

  std::vector<gra::M4Vec> R_daughters_in_X = {S_in_X, d_in_X};
  gra::kinematics::HXframe(R_daughters_in_X, R_in_X);
  const auto fR = gra::spin::fDecayMatrix(R.hel, R_daughters_in_X[0].Theta(), R_daughters_in_X[0].Phi());

  std::vector<gra::M4Vec> S_daughters_in_R_HX = {
      ApplyHXStepInCurrentFrame(BoostToRestFrame(S.legs[0].p4, lts.pfinal[0]), R_in_X, R_system_in_X),
      ApplyHXStepInCurrentFrame(BoostToRestFrame(S.legs[1].p4, lts.pfinal[0]), R_in_X, R_system_in_X)};
  gra::kinematics::HXframe(S_daughters_in_R_HX, R_daughters_in_X[0]);
  const auto fS = gra::spin::fDecayMatrix(S.hel, S_daughters_in_R_HX[0].Theta(), S_daughters_in_R_HX[0].Phi());

  MMatrix<std::complex<double>> expected(3, 1, 0.0);
  for (std::size_t M = 0; M < 3; ++M) {
    for (std::size_t lambdaR = 0; lambdaR < 3; ++lambdaR) {
      for (std::size_t lambdaS = 0; lambdaS < 3; ++lambdaS) {
        // Sequential helicity maps do not commute, so keep the X, R, S order
        const auto phase = std::exp(std::complex<double>(0.0, gra::math::PI *
            (static_cast<double>(lambdaR + lambdaS) - 2.0)));
        expected[M][0] += fX[lambdaR][M] * fR[lambdaS][lambdaR] * fS[0][lambdaS] * phase;
      }
    }
  }
  return expected;
}

MMatrix<std::complex<double>> RootDecayReferenceInFrame(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res,
                                                        const std::string &frame) {
  std::vector<gra::M4Vec> daughters = {lts.decaytree[0].p4, lts.decaytree[1].p4};

  if (frame == "HX") {
    gra::kinematics::HXframe(daughters, lts.pfinal[0]);
  } else if (frame == "CS") {
    gra::kinematics::CSframe(daughters, lts.pfinal[0], lts.pbeam1, lts.pbeam2);
  } else if (frame == "CM") {
    gra::kinematics::CMframe(daughters, lts.pfinal[0]);
  } else {
    throw std::invalid_argument("RootDecayReferenceInFrame: unsupported frame " + frame);
  }

  return gra::spin::fDecayMatrix(res.hel_decay, daughters[0].Theta(), daughters[0].Phi()).Transpose();
}

gra::M4Vec RotatePiAroundY(const gra::M4Vec &p) { return gra::M4Vec(-p.Px(), p.Py(), -p.Pz(), p.E()); }

gra::LORENTZSCALAR BeamExchangeMirror(const gra::LORENTZSCALAR &lts) {
  gra::LORENTZSCALAR out = lts;
  out.pbeam1             = lts.pbeam1;
  out.pbeam2             = lts.pbeam2;
  out.pfinal[1]          = RotatePiAroundY(lts.pfinal[2]);
  out.pfinal[2]          = RotatePiAroundY(lts.pfinal[1]);
  out.pfinal[0]          = RotatePiAroundY(lts.pfinal[0]);
  out.q1                 = out.pbeam1 - out.pfinal[1];
  out.q2                 = out.pbeam2 - out.pfinal[2];
  out.t1                 = out.q1.M2();
  out.t2                 = out.q2.M2();
  out.q1_in_X            = BoostToRestFrame(out.q1, out.pfinal[0]);
  out.q2_in_X            = BoostToRestFrame(out.q2, out.pfinal[0]);
  out.s1                 = (out.pbeam1 + out.q2).M2();
  out.s2                 = (out.pbeam2 + out.q1).M2();
  out.m2                 = out.pfinal[0].M2();
  return out;
}

gra::PARAM_RES MakeRealisticTensorMPResonance() {
  gra::PARAM_RES res;
  res.production_model = gra::ReggeProductionModel::MP;
  res.p                = ToyParticle("X2", 900225, 4, 1.2754);
  res.p.P              = 1;
  res.p.C              = 1;

  gra::MParticle exchange;
  exchange.name   = "pomeron(1)";
  exchange.pdg    = 993;
  exchange.spinX2 = 2;
  exchange.P      = -1;
  exchange.C      = 1;

  gra::MDecayBranch up;
  up.p                 = exchange;
  up.hel               = RealisticProtonLegHelicityMatrix(993, 2, -1, 1, 2);
  gra::MDecayBranch dn = up;
  res.production       = {{{up, dn}}};

  gra::HELMatrix central;
  central.BR         = 1.0;
  central.P_symmetry = true;
  central.alpha_ls.Set(0, 4, 1.0);
  gra::spin::InitTMatrix(central, res.p, exchange, exchange, true, "test MP tensor central", false, false);
  res.production.front().hel = central;
  res.a_Jz                   = {std::complex<double>(1.0 / std::sqrt(3.0), 0.0), std::complex<double>(0.0, 0.0),
                                std::complex<double>(1.0 / std::sqrt(3.0), 0.0), std::complex<double>(0.0, 0.0),
                                std::complex<double>(1.0 / std::sqrt(3.0), 0.0)};
  auto vertex = gra::spin::PreparePoleLS(res.p, exchange, exchange, {{0, 4, 1.0}}, 1.0, true, false, true);
  res.production.front().hel  = vertex.helicity;
  res.production.front().pole = std::move(vertex);
  return res;
}

double SteeredProductionAmp2(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res) {
  const auto channels = gra::rspin::Resonance(
      lts, res,
      (res.production_model == gra::ReggeProductionModel::MP ? gra::mpom::Fusion : gra::xpom::Fusion), 1.0);
  REQUIRE(channels.size() == 1);
  double out = 0.0;
  for (std::size_t row = 0; row < channels[0].size_row(); ++row) {
    std::complex<double> amp = 0.0;
    for (std::size_t jz = 0; jz < channels[0].size_col(); ++jz) { amp += channels[0][row][jz] * res.a_Jz[jz]; }
    out += std::norm(amp);
  }
  return out / static_cast<double>(channels[0].size_row());
}

std::string ComplexTableEntry(const std::complex<double> &z) {
  std::ostringstream out;
  out << std::scientific << std::setprecision(6) << "(" << z.real() << "," << z.imag() << ")";
  return out.str();
}

MMatrix<std::complex<double>> ReshapeFlatMatrix(const std::vector<std::complex<double>> &flat, std::size_t rows,
                                                std::size_t cols) {
  REQUIRE(flat.size() == rows * cols);
  MMatrix<std::complex<double>> out(rows, cols, 0.0);
  for (std::size_t i = 0; i < rows; ++i) {
    for (std::size_t j = 0; j < cols; ++j) { out[i][j] = flat[i * cols + j]; }
  }
  return out;
}

// Compute Regge parameters keyed by the exact event and model snapshots
std::shared_ptr<const gra::regge::Param> ReggeParametersForTest(const gra::MRegge        &regge,
                                                                const gra::LORENTZSCALAR &lts) {
  std::vector<int> final_pdgs;
  final_pdgs.reserve(lts.decaytree.size());
  for (const auto &branch : lts.decaytree) { final_pdgs.push_back(branch.p.pdg); }
  return gra::regge::ReadParamPtr(final_pdgs, lts.PDG, *regge.ModelTuneHandle());
}

// Compute the selected Regge resonance model form factor block
gra::RES_MODEL_FORM &ReggeResonanceFormForTest(gra::PARAM_RES &res, const gra::ReggeProductionModel model) {
  if (model == gra::ReggeProductionModel::MP) { return res.MP; }
  if (model == gra::ReggeProductionModel::XP) { return res.XP; }
  if (model == gra::ReggeProductionModel::GP) { return res.GP; }
  throw std::invalid_argument("ReggeResonanceFormForTest requires MP, XP or GP");
}

// Compute the selected Regge resonance model form factor block
const gra::RES_MODEL_FORM &ReggeResonanceFormForTest(const gra::PARAM_RES &res, const gra::ReggeProductionModel model) {
  if (model == gra::ReggeProductionModel::MP) { return res.MP; }
  if (model == gra::ReggeProductionModel::XP) { return res.XP; }
  if (model == gra::ReggeProductionModel::GP) { return res.GP; }
  throw std::invalid_argument("ReggeResonanceFormForTest requires MP, XP or GP");
}

// Compute the resonance form factors for one ordered production channel
double ResonanceFormFactorReference(const gra::LORENTZSCALAR &lts, const gra::PARAM_RES &res,
                                    const std::vector<gra::MDecayBranch> &tree) {
  REQUIRE(tree.size() == 2);
  const auto  &model        = ReggeResonanceFormForTest(res, res.production_model);
  const auto  &ff_prod      = model.ff_prod;
  const auto  &ff_decay     = res.hel_decay.ff_decay;
  const auto  &ff_transfer  = model.ff_transfer;
  double       value        = gra::regge::MassFF(lts.m2, gra::math::pow2(res.p.mass), ff_decay);
  value *= gra::regge::MassFF(lts.m2, gra::math::pow2(res.p.mass), ff_prod);
  const bool   upper_strong = tree[0].p.pdg != gra::PDG::PDG_gamma;
  const bool   lower_strong = tree[1].p.pdg != gra::PDG::PDG_gamma;
  if (upper_strong && lower_strong) {
    value *= gra::regge::TransferFF(lts.t1, ff_transfer);
    value *= gra::regge::TransferFF(lts.t2, ff_transfer);
  }
  return value;
}

// Compute the common scalar factor for one two-strong-leg MP resonance
std::complex<double> MPCommonFactor(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                    const gra::PARAM_RES &res) {
  const auto                 param       = ReggeParametersForTest(regge, lts);
  const gra::SoftExchangeId  pomeron     = param->exchanges.at(param->pomeron_trajectory).soft_exchange;
  const auto                &soft        = *regge.SoftModelHandle();
  const std::complex<double> A_prop_prop = soft.PhysicalResidue(pomeron, lts.t1) * regge.PomeronKernel(lts.s1, lts.t1) *
                                           regge.PomeronKernel(lts.s2, lts.t2) * soft.PhysicalResidue(pomeron, lts.t2);
  const std::complex<double> V = std::pow(param->s0 / lts.m2, param->omega.at(res.production_model));
  const std::complex<double> A_bw =
      gra::resonance::LineShape(lts.m2, res) * ResonanceFormFactorReference(lts, res, res.production[0].tree);

  return -res.hel_decay.g_decay * (A_prop_prop * V * A_bw);
}

// Compute the photon-emitter charge sign used at amplitude level
//
double PhotonBeamChargeSignReference(const gra::LORENTZSCALAR &lts, int leg) {
  const auto &beam = (leg == 1) ? lts.beam1 : lts.beam2;
  return (beam.chargeX3 < 0) ? -1.0 : 1.0;
}

double CoherentPhotonFluxReference(double x, double t, double pt) {
  const double pt2 = pow2(pt);
  if (!(x > 0.0 && x < 1.0) || !(pt2 > 0.0)) { return 0.0; }
  const double Q2 = std::abs(t);
  const double FE = (4.0 * pow2(gra::PDG::mp) * pow2(gra::form::G_E(Q2)) + Q2 * pow2(gra::form::G_M(Q2))) /
                    (4.0 * pow2(gra::PDG::mp) + Q2);
  const double FM    = pow2(gra::form::G_M(Q2));
  const double delta = pt2 / (pt2 + pow2(x * gra::PDG::mp));
  double flux = gra::qed::alpha_QED() / gra::math::PI * ((1.0 - x) * pow2(delta) * FE + (pow2(x) / 2.0) * delta * FM);
  flux /= x;
  flux /= pt2;
  flux *= 16.0 * gra::math::PIPI;
  return flux;
}

std::complex<double> PhotonAmplitudeFactorReference(const gra::LORENTZSCALAR &lts, int leg) {
  const bool   excite = (leg == 1) ? lts.excite1 : lts.excite2;
  const double x      = (leg == 1) ? lts.x1 : lts.x2;
  const double t      = (leg == 1) ? lts.t1 : lts.t2;
  const double qt     = (leg == 1) ? lts.qt1 : lts.qt2;
  const double m2     = (leg == 1) ? lts.pfinal[1].M2() : lts.pfinal[2].M2();
  if (lts.process.PHOTON_VERTEX == "QED" && !excite) { return std::complex<double>(1.0, 0.0); }
  REQUIRE(lts.model_cache != nullptr);
  const auto  &structure = lts.model_cache->Tune().Structure();
  const double flux = excite ? gra::flux::IncohFlux(x, t, qt, m2, structure) : gra::flux::CohFlux(x, t, qt, structure);
  return PhotonBeamChargeSignReference(lts, leg) * gra::math::msqrt(flux / x);
}

// Build one ordinary non-photon Regge forward-leg factor
//
std::complex<double> ReggeForwardLegFactorReference(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts, int leg,
                                                    int exchange_pdg) {
  REQUIRE((leg == 1 || leg == 2));
  const auto                param      = ReggeParametersForTest(regge, lts);
  const std::size_t         trajectory = gra::regge::TrajectoryIndex(*param, exchange_pdg);
  const gra::SoftExchangeId exchange   = param->exchanges.at(trajectory).soft_exchange;
  const double              s_forward  = (leg == 1) ? lts.s1 : lts.s2;
  const double              t_forward  = (leg == 1) ? lts.t1 : lts.t2;
  const int                 beam_pdg   = gra::ReggeBeamParticleForLeg(lts, leg).pdg;
  const double              beam_sign  = gra::regge::AntiparticleSign(*param, exchange_pdg, beam_pdg);
  return beam_sign * regge.SoftModelHandle()->PhysicalResidue(exchange, t_forward) *
         regge.ExchangeKernel(s_forward, t_forward, exchange_pdg);
}

// Check whether an exchange belongs to the configured Pomeron trajectory slot
//
bool IsPomeronTrajectoryReference(const gra::regge::Param &param, int exchange_pdg) {
  return gra::regge::TrajectoryIndex(param, exchange_pdg) == 0;
}

// Select the reference scalar prefactor for photon-containing production rows
//
std::complex<double> PhotoParamFactorReference(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                               const gra::PARAM_RES &res, const std::vector<gra::MDecayBranch> &tree) {
  const auto param = ReggeParametersForTest(regge, lts);
  REQUIRE_FALSE(lts.excite1);
  REQUIRE_FALSE(lts.excite2);
  const gra::SoftExchangeId pomeron   = param->exchanges.at(param->pomeron_trajectory).soft_exchange;
  const auto                photo_leg = [&regge, &lts, &res, pomeron](const int leg) {
    const double s = (leg == 1) ? lts.s1 : lts.s2;
    const double t = (leg == 1) ? lts.t1 : lts.t2;
    return regge.SoftModelHandle()->PhysicalResidue(pomeron, t) * regge.PhotoKernel(s, t, res.p.pdg);
  };
  const std::complex<double> A_1       = PhotonAmplitudeFactorReference(lts, 2) * photo_leg(1);
  const std::complex<double> A_2       = PhotonAmplitudeFactorReference(lts, 1) * photo_leg(2);
  const std::complex<double> A_default = A_1 + A_2;

  if (tree.size() != 2) { return A_default; }

  const bool upper_gamma = (tree[0].p.pdg == 22);
  const bool lower_gamma = (tree[1].p.pdg == 22);
  if (upper_gamma && lower_gamma) {
    return PhotonAmplitudeFactorReference(lts, 1) * PhotonAmplitudeFactorReference(lts, 2);
  }
  if (upper_gamma == lower_gamma) { return A_default; }

  const int regge_pdg = upper_gamma ? tree[1].p.pdg : tree[0].p.pdg;
  if (IsPomeronTrajectoryReference(*param, regge_pdg)) { return lower_gamma ? A_1 : A_2; }

  const int photon_leg = upper_gamma ? 1 : 2;
  const int regge_leg  = upper_gamma ? 2 : 1;
  return PhotonAmplitudeFactorReference(lts, photon_leg) *
         ReggeForwardLegFactorReference(regge, lts, regge_leg, regge_pdg);
}

std::complex<double> PhotoParamFactorReference(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                               const gra::PARAM_RES &res) {
  return PhotoParamFactorReference(regge, lts, res, res.production[0].tree);
}

std::complex<double> ReggeChannelFactorReference(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                                 const std::vector<gra::MDecayBranch> &tree) {
  REQUIRE(tree.size() == 2);
  const int up = tree[0].p.pdg;
  const int dn = tree[1].p.pdg;
  return ReggeForwardLegFactorReference(regge, lts, 1, up) * ReggeForwardLegFactorReference(regge, lts, 2, dn);
}

// Compute the scalar reference factor for one mixed production row
std::complex<double> MixedProductionCommonFactor(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                                 const gra::PARAM_RES                 &res,
                                                 const std::vector<gra::MDecayBranch> &tree) {
  const auto                 param = ReggeParametersForTest(regge, lts);
  const std::complex<double> A_bw =
      gra::resonance::LineShape(lts.m2, res) * ResonanceFormFactorReference(lts, res, tree);
  const bool                 has_photon = tree.size() == 2 && (tree[0].p.pdg == 22 || tree[1].p.pdg == 22);
  const std::complex<double> production =
      has_photon ? PhotoParamFactorReference(regge, lts, res, tree)
                 : ReggeChannelFactorReference(regge, lts, tree) * std::pow(param->s0 / lts.m2, param->omega.at(res.production_model));
  return -production * A_bw * res.hel_decay.g_decay;
}

// Compute the scalar reference factor for one photon production row
std::complex<double> PhotoCommonFactor(const gra::MRegge &regge, const gra::LORENTZSCALAR &lts,
                                       const gra::PARAM_RES &res) {
  const auto                 param = ReggeParametersForTest(regge, lts);
  const std::complex<double> A_bw =
      gra::resonance::LineShape(lts.m2, res) * ResonanceFormFactorReference(lts, res, res.production[0].tree);
  return -PhotoParamFactorReference(regge, lts, res) * A_bw * res.hel_decay.g_decay;
}

// Compute a finite helicity-averaged amplitude squared
//
double FiniteHelicityAmp2(const std::vector<std::complex<double>> &hamp) {
  const double amp2 = gra::SquaredNorm(hamp);
  REQUIRE(std::isfinite(amp2));
  return amp2;
}

// Update the derived toy kinematics and reset the simple two-body decay momenta
void UpdateToyDerivedKinematics(gra::LORENTZSCALAR &lts) {
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1        = lts.pbeam1 - lts.pfinal[1];
  lts.q2        = lts.pbeam2 - lts.pfinal[2];
  lts.t1        = lts.q1.M2();
  lts.t2        = lts.q2.M2();
  lts.s         = (lts.pbeam1 + lts.pbeam2).M2();
  lts.sqrt_s    = std::sqrt(std::max(0.0, lts.s));
  lts.m2        = lts.pfinal[0].M2();
  lts.s1        = (lts.pfinal[0] + lts.pfinal[1]).M2();
  lts.s2        = (lts.pfinal[0] + lts.pfinal[2]).M2();
  lts.x1        = 1.0 - lts.pfinal[1].Pz() / lts.pbeam1.Pz();
  lts.x2        = 1.0 - lts.pfinal[2].Pz() / lts.pbeam2.Pz();
  lts.xi1       = lts.x1;
  lts.xi2       = lts.x2;
  lts.has_xi1   = true;
  lts.has_xi2   = true;
  lts.qt1       = lts.q1.Pt();
  lts.qt2       = lts.q2.Pt();
  lts.q1_in_X   = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X   = BoostToRestFrame(lts.q2, lts.pfinal[0]);

  if (lts.decaytree.size() == 2) {
    SetCentralPair(lts);
    lts.ss[1][3] = lts.ss[3][1] = (lts.pfinal[1] + lts.decaytree[0].p4).M2();
    lts.ss[2][4] = lts.ss[4][2] = (lts.pfinal[2] + lts.decaytree[1].p4).M2();
    lts.ss[1][4] = lts.ss[4][1] = (lts.pfinal[1] + lts.decaytree[1].p4).M2();
    lts.ss[2][3] = lts.ss[3][2] = (lts.pfinal[2] + lts.decaytree[0].p4).M2();
    lts.s_hat                   = lts.pfinal[0].M2();
    lts.t_hat                   = (lts.q1 - lts.decaytree[0].p4).M2();
    lts.u_hat                   = (lts.q1 - lts.decaytree[1].p4).M2();
    lts.d0_in_X                 = BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0]);
    lts.d1_in_X                 = BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0]);
  }
}

// Refresh the derived toy invariants without overwriting the central decay
// momenta
void RefreshToyDerivedKinematicsPreserveDecay(gra::LORENTZSCALAR &lts) {
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1        = lts.pbeam1 - lts.pfinal[1];
  lts.q2        = lts.pbeam2 - lts.pfinal[2];
  lts.t1        = lts.q1.M2();
  lts.t2        = lts.q2.M2();
  lts.s         = (lts.pbeam1 + lts.pbeam2).M2();
  lts.sqrt_s    = std::sqrt(std::max(0.0, lts.s));
  lts.m2        = lts.pfinal[0].M2();
  lts.s1        = (lts.pfinal[0] + lts.pfinal[1]).M2();
  lts.s2        = (lts.pfinal[0] + lts.pfinal[2]).M2();
  lts.x1        = 1.0 - lts.pfinal[1].Pz() / lts.pbeam1.Pz();
  lts.x2        = 1.0 - lts.pfinal[2].Pz() / lts.pbeam2.Pz();
  lts.xi1       = lts.x1;
  lts.xi2       = lts.x2;
  lts.has_xi1   = true;
  lts.has_xi2   = true;
  lts.qt1       = lts.q1.Pt();
  lts.qt2       = lts.q2.Pt();
  lts.q1_in_X   = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X   = BoostToRestFrame(lts.q2, lts.pfinal[0]);

  if (lts.decaytree.size() == 2) {
    lts.ss[1][3] = lts.ss[3][1] = (lts.pfinal[1] + lts.decaytree[0].p4).M2();
    lts.ss[2][4] = lts.ss[4][2] = (lts.pfinal[2] + lts.decaytree[1].p4).M2();
    lts.ss[1][4] = lts.ss[4][1] = (lts.pfinal[1] + lts.decaytree[1].p4).M2();
    lts.ss[2][3] = lts.ss[3][2] = (lts.pfinal[2] + lts.decaytree[0].p4).M2();
    lts.s_hat                   = lts.pfinal[0].M2();
    lts.t_hat                   = (lts.q1 - lts.decaytree[0].p4).M2();
    lts.u_hat                   = (lts.q1 - lts.decaytree[1].p4).M2();
    lts.d0_in_X                 = BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0]);
    lts.d1_in_X                 = BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0]);
  }
}

// Rotate one four-vector around the beam axis for azimuthal covariance tests
gra::M4Vec RotateZForTest(gra::M4Vec p, double angle) {
  p.RotateZ(angle);
  return p;
}

// Rotate all event momenta around the beam axis while preserving the toy decay
gra::LORENTZSCALAR RotateToyEventAroundZ(gra::LORENTZSCALAR lts, double angle) {
  lts.pbeam1.RotateZ(angle);
  lts.pbeam2.RotateZ(angle);
  for (auto &p : lts.pfinal) { p.RotateZ(angle); }
  const auto rotate = [angle](const auto &self, gra::MDecayBranch &branch) -> void {
    branch.p4.RotateZ(angle);
    for (auto &leg : branch.legs) { self(self, leg); }
  };
  for (auto &branch : lts.decaytree) { rotate(rotate, branch); }
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Reflect all generated momenta in the collider xz plane
gra::LORENTZSCALAR ReflectToyEventInXZ(gra::LORENTZSCALAR lts) {
  const auto reflect = [](const gra::M4Vec &p) { return gra::M4Vec(p.Px(), -p.Py(), p.Pz(), p.E()); };
  lts.pbeam1         = reflect(lts.pbeam1);
  lts.pbeam2         = reflect(lts.pbeam2);
  for (auto &p : lts.pfinal) { p = reflect(p); }
  const auto reflect_branch = [&](const auto &self, gra::MDecayBranch &current) -> void {
    current.p4 = reflect(current.p4);
    for (auto &leg : current.legs) { self(self, leg); }
  };
  for (auto &branch : lts.decaytree) { reflect_branch(reflect_branch, branch); }
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Exchange the beam directions and rotate the complete central decay tree
gra::LORENTZSCALAR BeamExchangeMirrorWithDecay(gra::LORENTZSCALAR lts) {
  std::swap(lts.beam1, lts.beam2);
  std::swap(lts.pbeam1, lts.pbeam2);
  lts.pbeam1 = RotatePiAroundY(lts.pbeam1);
  lts.pbeam2 = RotatePiAroundY(lts.pbeam2);
  std::swap(lts.excite1, lts.excite2);
  const gra::M4Vec upper   = lts.pfinal[1];
  const gra::M4Vec lower   = lts.pfinal[2];
  lts.pfinal[1]            = RotatePiAroundY(lower);
  lts.pfinal[2]            = RotatePiAroundY(upper);
  lts.pfinal[0]            = RotatePiAroundY(lts.pfinal[0]);
  const auto rotate_branch = [&](const auto &self, gra::MDecayBranch &current) -> void {
    current.p4 = RotatePiAroundY(current.p4);
    for (auto &leg : current.legs) { self(self, leg); }
  };
  for (auto &branch : lts.decaytree) { rotate_branch(rotate_branch, branch); }
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Boost all toy event momenta along the beam axis and refresh invariants
gra::LORENTZSCALAR BoostToyEventAlongZ(gra::LORENTZSCALAR lts, const double rapidity) {
  const gra::M4Vec boost(0.0, 0.0, std::sinh(rapidity), std::cosh(rapidity));
  gra::kinematics::LorentzBoost(boost, 1.0, lts.pbeam1, +1);
  gra::kinematics::LorentzBoost(boost, 1.0, lts.pbeam2, +1);
  gra::kinematics::LorentzBoost(boost, 1.0, lts.pfinal[1], +1);
  gra::kinematics::LorentzBoost(boost, 1.0, lts.pfinal[2], +1);
  for (auto &branch : lts.decaytree) { BoostDecayBranchForTest(branch, boost, 1.0); }
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

gra::LORENTZSCALAR MakeToyReggeLTSAsymmetric(double proton_pt, double pz1, double pz2) {
  gra::LORENTZSCALAR lts = MakeToyContinuumLTSAsymmetric(proton_pt, pz1, pz2);
  UpdateToyDerivedKinematics(lts);
  return lts;
}

gra::LORENTZSCALAR MakeToyCoherentPhotonLTS() {
  gra::LORENTZSCALAR lts   = MakeToyContinuumLTSAsymmetric(0.25, 4.55, -4.50);
  const auto         beams = ProtonInitialState();
  lts.beam1                = beams[0];
  lts.beam2                = beams[1];
  const double pz_beam     = 5.0;
  const double e_beam      = std::sqrt(pow2(pz_beam) + pow2(gra::PDG::mp));
  lts.pbeam1               = gra::M4Vec(0.0, 0.0, pz_beam, e_beam);
  lts.pbeam2               = gra::M4Vec(0.0, 0.0, -pz_beam, e_beam);
  lts.pfinal[1].SetPxPyPzM(0.24, 0.11, 4.55, gra::PDG::mp);
  lts.pfinal[2].SetPxPyPzM(-0.18, 0.07, -4.50, gra::PDG::mp);
  UpdateToyDerivedKinematics(lts);
  return lts;
}

void SetToyContinuumExchangePair(gra::LORENTZSCALAR &lts, int upper_pdg, int lower_pdg) {
  lts.process.CONT_TU_SIGN.clear();
  REQUIRE(lts.process.CONT_PRODUCTIONTREE.size() == 1);
  REQUIRE(lts.process.CONT_PRODUCTIONTREE[0].size() == 2);
  lts.process.CONT_PRODUCTION             = {{upper_pdg, lower_pdg}};
  lts.process.CONT_PRODUCTIONTREE[0][0].p = lts.PDG.FindByPDG(upper_pdg);
  lts.process.CONT_PRODUCTIONTREE[0][1].p = lts.PDG.FindByPDG(lower_pdg);
  for (auto &branch : lts.process.CONT_PRODUCTIONTREE[0]) {
    if (branch.p.spinX2 >= 0) { branch.hel = ToyProtonLegHelicityMatrix(branch.p.spinX2); }
  }
  if (!lts.process.CONTINUUM_POLE.empty()) {
    if (lts.process.CONT_PRODUCTIONTREE[0][0].p.spinX2 < 0 || lts.process.CONT_PRODUCTIONTREE[0][1].p.spinX2 < 0) {
      lts.process.CONTINUUM_POLE.clear();
    } else {
      PrepareToyContinuumOperators(lts, gra::ReggeProductionModel::MP);
    }
  }
}

// Replace toy subchannel data by explicit analytic GP LS vertices
void UseToyGPSubchannelHelicity(gra::LORENTZSCALAR &lts) {
  REQUIRE(lts.decaytree.size() == 2);

  auto trajectory_spin = [](const gra::MParticle &particle) {
    if (particle.pdg == 22) { return 1U; }
    if (particle.pdg == 990 || particle.pdg == 9910 || particle.pdg == 9920) { return 2U; }
    if (particle.pdg == 9930 || particle.pdg == 9940 || particle.pdg == 9990) { return 1U; }
    throw std::invalid_argument("UseToyGPSubchannelHelicity: unsupported analytic exchange");
  };

  lts.process.CONTINUUM_POLE.clear();
  lts.process.CONTINUUM_GP.clear();
  for (const auto &channel : indices(lts.process.CONT_PRODUCTIONTREE)) {
    const auto &tree = lts.process.CONT_PRODUCTIONTREE[channel];
    REQUIRE(tree.size() == 2);
    const std::array<std::vector<gra::MParticle>, 4> legs = {
        std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p},
        std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
        std::vector<gra::MParticle>{lts.decaytree[1].p, lts.decaytree[0].p},
        std::vector<gra::MParticle>{lts.decaytree[0].p, lts.decaytree[1].p}};
    const std::array<std::size_t, 4> orbital = {trajectory_spin(tree[0].p), trajectory_spin(tree[1].p),
                                                trajectory_spin(tree[0].p), trajectory_spin(tree[1].p)};
    std::vector<gra::HELMatrix>      cache;
    cache.reserve(4);
    for (const auto &vertex : indices(legs)) {
      gra::HELMatrix hel;
      const bool     photon = tree[vertex % 2].p.pdg == gra::PDG::PDG_gamma;
      if (photon) {
        gra::gpom::InitPhotonCrossed(hel, legs[vertex][0].spinX2 / 2.0, legs[vertex][1].spinX2 / 2.0,
                                     "UseToyGPSubchannelHelicity");
      } else {
        gra::gpom::InitCrossed(hel, legs[vertex][0].spinX2 / 2.0, legs[vertex][1].spinX2 / 2.0, lts.process.MMAX,
                               "UseToyGPSubchannelHelicity");
      }
      hel.coupling_basis = gra::CouplingBasis::LS;
      if (photon) {
        hel.alpha_ls.Set(orbital[vertex], 0, 1.0);
      } else {
        hel.m_ls.resize(2 * lts.process.MMAX + 1);
        hel.m_ls[lts.process.MMAX].Set(orbital[vertex], 0, 1.0);
      }
      gra::gpom::InitCrossedLS(hel, static_cast<int>(trajectory_spin(tree[vertex % 2].p)));
      cache.push_back(std::move(hel));
    }
    lts.process.CONTINUUM_GP.push_back(std::move(cache));
  }
}

// Replace toy production matrices by canonical fixed-spin pole operators
void PrepareToyPoleOperators(gra::PARAM_RES &res, const gra::ReggeProductionModel model,
                             const std::vector<std::vector<gra::spin::LSTerm>> &terms, bool C_symmetry = false,
                             bool P_symmetry = true) {
  REQUIRE((model == gra::ReggeProductionModel::MP || model == gra::ReggeProductionModel::XP));
  REQUIRE(terms.size() == res.production.size());
  res.production_model = model;
  for (const auto &channel : indices(res.production)) {
    REQUIRE(res.production[channel].tree.size() == 2);
    auto vertex = gra::spin::PreparePoleLS(res.p, res.production[channel].tree[0].p, res.production[channel].tree[1].p,
                                           terms[channel], 1.0, true, C_symmetry, P_symmetry);
    res.production[channel].hel  = vertex.helicity;
    res.production[channel].pole = std::move(vertex);
  }
}

gra::PARAM_RES MakeToyMPResonance(bool xp_mode = false) {
  gra::PARAM_RES res    = MakeToyResonance();
  res.production_model  = xp_mode ? gra::ReggeProductionModel::XP : gra::ReggeProductionModel::MP;
  auto &form            = ReggeResonanceFormForTest(res, res.production_model);
  form.ff_prod          = {gra::regge::FFType::Gaussian, gra::regge::FFNorm::Pole, {0.64}};
  res.BW                = gra::BreitWigner::FixedWidth;
  res.p.mass            = 0.75;
  res.p.width           = 0.08;
  res.hel_decay         = SpinOneCentralResonanceHelicityMatrix();
  res.hel_decay.ff_decay          = form.ff_prod;
  res.hel_decay.g_decay = std::complex<double>(0.6, 0.1);
  res.a_Jz = {std::complex<double>(1.0, 0.0), std::complex<double>(0.3, -0.2), std::complex<double>(-0.25, 0.15)};
  if (!xp_mode) { PrepareToyPoleOperators(res, gra::ReggeProductionModel::MP, {{{2, 4, 1.0}}}); }
  return res;
}

// Replace toy production matrices by raw canonical XP operator caches
void PrepareToyXPOperators(gra::PARAM_RES &res, const std::vector<std::vector<gra::spin::LSTerm>> &terms,
                           bool C_symmetry = false, bool P_symmetry = true) {
  PrepareToyPoleOperators(res, gra::ReggeProductionModel::XP, terms, C_symmetry, P_symmetry);
}

// Build the same toy resonance with one absolute u-STF XP operator
gra::PARAM_RES MakeToyCovariantXPonance() {
  gra::PARAM_RES res = MakeToyMPResonance(true);
  PrepareToyXPOperators(res, {{{2, 4, 1.0}}});
  return res;
}

gra::PARAM_RES MakeToyScalarMPResonance() {
  gra::PARAM_RES res;
  res.production_model    = gra::ReggeProductionModel::MP;
  res.MP.ff_prod          = {gra::regge::FFType::Gaussian, gra::regge::FFNorm::Pole, {0.64}};
  res.BW                  = gra::BreitWigner::FixedWidth;

  gra::MParticle exchange;
  exchange.name   = "toy_scalar_exchange";
  exchange.pdg    = 991;
  exchange.spinX2 = 0;
  exchange.P      = 1;

  gra::MParticle mother;
  mother.name   = "toy_scalar_resonance";
  mother.pdg    = 9001001;
  mother.spinX2 = 0;
  mother.P      = 1;
  mother.mass   = 0.75;
  mother.width  = 0.08;
  res.p         = mother;

  gra::MDecayBranch up;
  up.p   = exchange;
  up.hel = RealisticProtonLegHelicityMatrix(991, 0, 1, 0, 0);

  gra::MDecayBranch dn  = up;
  res.production        = {{{up, dn}, SpinZeroCentralResonanceHelicityMatrix()}};
  res.hel_decay         = SpinZeroCentralResonanceHelicityMatrix();
  res.hel_decay.ff_decay          = res.MP.ff_prod;
  res.hel_decay.g_decay = std::complex<double>(1.0, 0.0);
  res.a_Jz              = {std::complex<double>(1.0, 0.0)};
  PrepareToyPoleOperators(res, gra::ReggeProductionModel::MP, {{{0, 0, 1.0}}});
  return res;
}

// Build the vector-vector pseudoscalar XP model used by the dphi regression
gra::PARAM_RES MakeToyPseudoscalarVectorXP() {
  gra::PARAM_RES res;
  res.production_model = gra::ReggeProductionModel::XP;

  gra::MParticle exchange;
  exchange.name   = "pomeron(1)";
  exchange.pdg    = 993;
  exchange.spinX2 = 2;
  exchange.P      = -1;
  exchange.C      = 1;

  gra::MParticle mother;
  mother.name   = "toy_eta";
  mother.pdg    = 9000221;
  mother.spinX2 = 0;
  mother.P      = -1;
  mother.C      = 1;
  mother.mass   = 0.957;
  mother.width  = 0.0002;
  res.p         = mother;

  gra::MDecayBranch up;
  up.p                 = exchange;
  up.hel               = RealisticProtonLegHelicityMatrix(993, 2, -1, 1, 2);
  gra::MDecayBranch dn = up;
  res.production       = {{{up, dn}}};

  res.hel_decay         = SpinZeroCentralResonanceHelicityMatrix();
  res.hel_decay.g_decay = std::complex<double>(1.0, 0.0);
  res.a_Jz              = {std::complex<double>(1.0, 0.0)};
  PrepareToyXPOperators(res, {{{1, 2, 1.0}}});
  return res;
}

// Build the vector-vector pseudoscalar MP model used by the dphi regression
gra::PARAM_RES MakeToyPseudoscalarVectorRes() {
  gra::PARAM_RES res = MakeToyPseudoscalarVectorXP();
  PrepareToyPoleOperators(res, gra::ReggeProductionModel::MP, {{{1, 2, 1.0}}});
  return res;
}

// Build the tensor-tensor pseudoscalar XP model used by the dphi regression
gra::PARAM_RES MakeToyPseudoscalarTensorXP() {
  gra::PARAM_RES res;
  res.production_model = gra::ReggeProductionModel::XP;

  gra::MParticle exchange;
  exchange.name   = "pomeron(2)";
  exchange.pdg    = 995;
  exchange.spinX2 = 4;
  exchange.P      = 1;
  exchange.C      = 1;

  gra::MParticle mother;
  mother.name   = "toy_eta_tensor";
  mother.pdg    = 9000222;
  mother.spinX2 = 0;
  mother.P      = -1;
  mother.C      = 1;
  mother.mass   = 0.957;
  mother.width  = 0.0002;
  res.p         = mother;

  gra::MDecayBranch up;
  up.p                 = exchange;
  up.hel               = RealisticProtonLegHelicityMatrix(995, 4, 1, 2, 0);
  gra::MDecayBranch dn = up;
  res.production       = {{{up, dn}}};

  res.hel_decay         = SpinZeroCentralResonanceHelicityMatrix();
  res.hel_decay.g_decay = std::complex<double>(1.0, 0.0);
  res.a_Jz              = {std::complex<double>(1.0, 0.0)};
  PrepareToyXPOperators(res, {{{1, 2, 1.0}}});
  return res;
}

gra::PARAM_RES MakeToyPhotoMPResonance(bool xp_mode = false) {
  gra::PARAM_RES res = MakeToyMPResonance(xp_mode);
  res.p.P            = -1;
  res.p.pdg          = 113;

  gra::MParticle photon;
  photon.name   = "gamma";
  photon.pdg    = 22;
  photon.spinX2 = 2;
  photon.P      = -1;

  gra::MParticle exchange;
  exchange.name   = "toy_exchange";
  exchange.pdg    = 993;
  exchange.spinX2 = 2;
  exchange.P      = 1;

  gra::MDecayBranch up;
  up.p   = photon;
  up.hel = SpinHalfLegHelicityMatrix();

  gra::MDecayBranch dn;
  dn.p   = exchange;
  dn.hel = SpinHalfLegHelicityMatrix();

  res.production = {{{up, dn}, SpinOneCentralResonanceHelicityMatrix()}};
  if (!xp_mode) { PrepareToyPoleOperators(res, gra::ReggeProductionModel::MP, {{{0, 2, 1.0}}}, false, true); }
  return res;
}

// Build the spin-one photon plus exchange toy with one XP operator
gra::PARAM_RES MakeToyCovariantPhotoXPonance() {
  gra::PARAM_RES res                                = MakeToyPhotoMPResonance(true);
  res.production.front().tree.front().hel.Jz_values = {-1.0, 1.0};
  PrepareToyXPOperators(res, {{{0, 2, 1.0}}});
  return res;
}

// Build the l=0, S=1 rho <- gamma + scalar-Pomeron central matrix
gra::HELMatrix RhoGammaScalarCentralHelicityMatrix(const gra::MParticle &first, const gra::MParticle &second) {
  gra::MParticle rho;
  rho.name   = "toy_rho";
  rho.pdg    = 9000113;
  rho.spinX2 = 2;
  rho.P      = -1;
  rho.C      = -1;

  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = true;
  hel.C_symmetry = true;
  hel.alpha_ls.Set(0, 2, 1.0);
  gra::spin::InitTMatrix(hel, rho, first, second, true, "test rho gamma scalar central", false, false);
  return hel;
}

// Build a rho-like gamma plus scalar-Pomeron XP model
gra::PARAM_RES MakeToyScalarPomeronPhotoXP() {
  gra::PARAM_RES res = MakeToyPhotoMPResonance(true);
  res.p.name         = "toy_rho";
  res.p.pdg          = 9000113;
  res.p.spinX2       = 2;
  res.p.P            = -1;
  res.p.C            = -1;

  gra::MParticle photon;
  photon.name   = "gamma";
  photon.pdg    = 22;
  photon.spinX2 = 2;
  photon.P      = -1;
  photon.C      = -1;

  gra::MParticle pomeron0;
  pomeron0.name   = "pomeron(0)";
  pomeron0.pdg    = 991;
  pomeron0.spinX2 = 0;
  pomeron0.P      = 1;
  pomeron0.C      = 1;

  gra::MDecayBranch up;
  up.p             = photon;
  up.hel           = RealisticProtonLegHelicityMatrix(22, 2, -1, 1, 2);
  up.hel.Jz_values = {-1.0, 1.0};

  gra::MDecayBranch dn;
  dn.p   = pomeron0;
  dn.hel = RealisticProtonLegHelicityMatrix(991, 0, 1, 0, 0);

  res.production = {{{up, dn}, RhoGammaScalarCentralHelicityMatrix(photon, pomeron0)}};
  res.hel_decay  = SpinOneToScalarScalarHelicityMatrix(std::complex<double>(1.0, 0.0));
  res.a_Jz       = {std::complex<double>(0.9, 0.1), std::complex<double>(0.2, -0.3), std::complex<double>(-0.4, 0.25)};
  return res;
}

// Build the photon plus scalar-Pomeron toy with an absolute XP operator
gra::PARAM_RES MakeToyCovariantScalarPomeronPhotoXP() {
  gra::PARAM_RES res                                = MakeToyScalarPomeronPhotoXP();
  res.p.pdg                                         = 113;
  res.production.front().tree.front().hel.Jz_values = {-1.0, 1.0};
  const auto &tree                                  = res.production.front().tree;
  auto        vertex = gra::spin::PreparePoleLS(res.p, tree[0].p, tree[1].p, {{0, 2, 1.0}}, 1.0, true, true, true);
  res.production.front().hel  = vertex.helicity;
  res.production.front().pole = std::move(vertex);
  return res;
}

// Build a photon plus analytic trajectory GP resonance model
gra::PARAM_RES MakeToyPhotoGPResonance(int MMAX) {
  gra::PARAM_RES res = MakeToyPhotoMPResonance(true);
  res.p.C            = -1;

  gra::MParticle exchange;
  exchange.name   = "pomeron(J=a(t))";
  exchange.pdg    = 990;
  exchange.spinX2 = gra::aux::kNullSpinX2;
  exchange.P      = 1;
  exchange.C      = 1;

  res.production[0].tree[1].p = exchange;

  gra::HELMatrix central;
  gra::gpom::InitHelicity(central, 0.0, 0.0, MMAX, "MakeToyPhotoGPResonance");
  central.BR             = 1.0;
  central.P_symmetry     = true;
  central.C_symmetry     = true;
  central.coupling_basis = gra::CouplingBasis::LS;
  central.alpha_ls.Set(0, 2, 1.0);
  gra::gpom::InitResonanceLS(central, res.p.spinX2 / 2, 2, 4);
  res.production.front().hel = central;
  return res;
}

// Build a hadronic analytic-trajectory GP model
//
gra::PARAM_RES MakeToyGPResonance(int MMAX) {
  gra::PARAM_RES res      = MakeToyMPResonance(true);
  res.production_model    = gra::ReggeProductionModel::GP;
  res.GP.ff_prod          = res.XP.ff_prod;

  gra::MParticle exchange;
  exchange.name   = "pomeron(J=a(t))";
  exchange.pdg    = 990;
  exchange.spinX2 = gra::aux::kNullSpinX2;
  exchange.P      = 1;
  exchange.C      = 1;

  res.production[0].tree[0].p = exchange;
  res.production[0].tree[1].p = exchange;

  gra::HELMatrix central;
  gra::gpom::InitHelicity(central, 0.0, 0.0, MMAX, "MakeToyGPResonance");
  central.BR             = 1.0;
  central.P_symmetry     = true;
  central.C_symmetry     = true;
  central.coupling_basis = gra::CouplingBasis::LS;
  central.alpha_ls.Set(0, 2, 1.0);
  gra::gpom::InitResonanceLS(central, res.p.spinX2 / 2, 4, 4);
  res.production.front().hel = central;
  return res;
}

// Build a scalar analytic-trajectory GP model with an S=0 central row
//
gra::PARAM_RES MakeToyScalarGPResonance(int MMAX) {
  gra::PARAM_RES res = MakeToyScalarMPResonance();

  gra::MParticle exchange;
  exchange.name   = "pomeron(J=a(t))";
  exchange.pdg    = 990;
  exchange.spinX2 = gra::aux::kNullSpinX2;
  exchange.P      = 1;
  exchange.C      = 1;

  res.production_model        = gra::ReggeProductionModel::GP;
  res.GP.ff_prod              = res.MP.ff_prod;
  res.production[0].tree[0].p = exchange;
  res.production[0].tree[1].p = exchange;

  gra::HELMatrix central;
  gra::gpom::InitHelicity(central, 0.0, 0.0, MMAX, "MakeToyScalarGPResonance");
  central.BR             = 1.0;
  central.P_symmetry     = true;
  central.C_symmetry     = true;
  central.coupling_basis = gra::CouplingBasis::LS;
  central.alpha_ls.Set(0, 0, 1.0);
  gra::gpom::InitResonanceLS(central, res.p.spinX2 / 2, 4, 4);
  res.production.front().hel = central;
  return res;
}

// Build a tensor analogue of the toy analytic GP model
//
gra::PARAM_RES MakeToyTensorGPResonance(int MMAX) {
  gra::PARAM_RES res = MakeToyGPResonance(MMAX);
  res.p.spinX2       = 4;

  gra::HELMatrix decay;
  gra::spin::InitTwoBodyBasis(decay, 2.0, 1.0, 1.0, {-2.0, -1.0, 0.0, 1.0, 2.0}, {-1.0, 0.0, 1.0}, {-1.0, 0.0, 1.0},
                              "MakeToyTensorGPResonance decay");
  decay.T       = MMatrix<std::complex<double>>(3, 3, 0.0);
  decay.T[1][1] = 1.0;
  decay.g_decay = std::complex<double>(0.6, 0.1);
  res.hel_decay = decay;

  gra::HELMatrix central;
  gra::gpom::InitHelicity(central, 2.0, 2.0, MMAX, "MakeToyTensorGPResonance");
  central.BR             = 1.0;
  central.P_symmetry     = true;
  central.C_symmetry     = true;
  central.coupling_basis = gra::CouplingBasis::LS;
  central.alpha_ls.Set(0, 4, 1.0);
  gra::gpom::InitResonanceLS(central, res.p.spinX2 / 2, 4, 4);
  res.production.front().hel = central;
  return res;
}

struct ModelParamRestoreGuard {
  std::string value = gra::MODELPARAM;
  ~ModelParamRestoreGuard() { gra::MODELPARAM = value; }
};

// Compute the canonical proton pair from the particle data used by processes
std::vector<gra::MParticle> ProtonInitialState() {
  const gra::MParticle proton = LoadedPDGTable().FindByPDG(gra::PDG::PDG_p);
  return {proton, proton};
}

// Build a proton-antiproton initial state for crossing-sign tests
std::vector<gra::MParticle> ProtonAntiprotonInitialState() {
  const gra::MParticle proton     = LoadedPDGTable().FindByPDG(gra::PDG::PDG_p);
  const gra::MParticle antiproton = LoadedPDGTable().FindByPDG(-gra::PDG::PDG_p);
  return {proton, antiproton};
}

// Write a temporary soft model file with controlled multichannel couplings
std::string WriteSoftModelFile(const std::vector<double> &theta, const std::string &suffix,
                               double odderon_coupling = 0.25, int odderon_sign = 1, double odderon_transition = 0.0,
                               double pomeron_kappa = 0.0, const std::string &unitarization = {},
                               double eikonal_q = std::numeric_limits<double>::quiet_NaN()) {
  const std::string data      = gra::aux::GetInputData(modelfile);
  auto              j         = nlohmann::json::parse(data);
  std::size_t       nchannels = 1;
  while (nchannels * (nchannels - 1) / 2 < theta.size()) { ++nchannels; }
  REQUIRE(nchannels * (nchannels - 1) / 2 == theta.size());

  const std::string model = "test_multichannel";
  auto             &soft  = j["PARAM_SOFT"];
  soft["MODEL"][model]    = soft["MODEL"]["single"];

  soft["active_model"]               = model;
  auto &model_block                  = soft["MODEL"][model];
  model_block["GW"]["theta"]         = theta;
  model_block["EIKONAL"]["helicity"] = !gra::math::IsZero(pomeron_kappa);
  if (!unitarization.empty()) { model_block["EIKONAL"]["unitarization"] = unitarization; }
  if (std::isfinite(eikonal_q)) { model_block["EIKONAL"]["q"] = eikonal_q; }
  model_block["EXCHANGE"]["P"]["helicity"]["kappa"] = pomeron_kappa;
  std::vector<double> resolved_direction(nchannels > 0 ? nchannels - 1 : 0, 0.0);
  double              direction_norm2 = 0.0;
  for (std::size_t c = 0; c < resolved_direction.size(); ++c) {
    resolved_direction[c] = static_cast<double>(c + 1);
    direction_norm2 += gra::math::pow2(resolved_direction[c]);
  }
  for (double &coefficient : resolved_direction) { coefficient /= std::sqrt(direction_norm2); }
  model_block["GW"]["a_c"] = resolved_direction;
  std::vector<double> pomeron_coupling(nchannels, 1.0);
  std::vector<double> odderon_couplings(nchannels, odderon_coupling);
  for (std::size_t i = 0; i < nchannels; ++i) {
    pomeron_coupling[i] *= 1.0 + 0.2 * i;
    odderon_couplings[i] *= 1.0 + 0.1 * i;
  }

  auto       &exchange_block = model_block["EXCHANGE"];
  const auto &definitions    = soft["EXCHANGE_DEF"];
  for (auto it = definitions.begin(); it != definitions.end(); ++it) {
    const std::string name     = it.key();
    const double      value    = (it.value().at("crossing").get<int>() == -1 && name != "O")
                                     ? 0.0
                                     : exchange_block.at(name).at("g").at(0).at(0).get<double>();
    auto             &coupling = exchange_block[name]["g"];
    coupling                   = nlohmann::json::array();
    for (std::size_t i = 0; i < nchannels; ++i) {
      coupling.push_back(nlohmann::json::array());
      for (std::size_t j = 0; j < nchannels; ++j) { coupling[i].push_back((i == j) ? value : 0.0); }
    }
  }
  for (std::size_t i = 0; i < nchannels; ++i) {
    exchange_block["P"]["g"][i][i] = pomeron_coupling[i];
    exchange_block["O"]["g"][i][i] = odderon_couplings[i];
  }
  for (std::size_t i = 0; i < nchannels; ++i) {
    for (std::size_t j = i + 1; j < nchannels; ++j) {
      exchange_block["O"]["g"][i][j] = odderon_transition;
      exchange_block["O"]["g"][j][i] = odderon_transition;
    }
  }
  REQUIRE((odderon_sign == -1 || odderon_sign == 1));
  exchange_block["O"]["sign"] = odderon_sign;
  for (auto it = model_block["FF"].begin(); it != model_block["FF"].end(); ++it) {
    auto      &parameters = it.value()["param"];
    const auto row        = parameters.at(0);
    parameters            = std::vector<nlohmann::json>(nchannels, row);
  }

  const std::filesystem::path directory = "tmp/graniitti_soft_multichannel_" + suffix;
  std::filesystem::create_directories(directory);
  const std::filesystem::path path = directory / "GENERAL.json";
  std::ofstream               out(path);
  REQUIRE(out.good());
  out << j.dump(2);
  out.close();

  std::ofstream numerics_out(directory / "NUMERICS.json");
  REQUIRE(numerics_out.good());
  numerics_out << gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json"));
  numerics_out.close();
  for (const std::string model : {"MP", "XP", "GP", "TP"}) {
    const std::string filename = "CON_" + model + ".json";
    std::filesystem::copy_file(gra::ResolveModelDataFile("TUNE0", filename), directory / filename,
                               std::filesystem::copy_options::overwrite_existing);
  }
  return path.string();
}

// Build a small eikonal object for a chosen initial state
gra::MEikonal BuildTestEikonalWithInitialState(const std::vector<double> &theta, const std::string &suffix,
                                               const std::vector<gra::MParticle> &initialstate,
                                               double odderon_coupling = 0.25, int odderon_sign = 1, int number_bt = 4,
                                               int number_kt2 = 4, double mandelstam_s = 25.0,
                                               double odderon_transition = 0.0, double pomeron_kappa = 0.0,
                                               int number_loop_kt = 0, int number_loop_phi = 0,
                                               const std::string &unitarization = {},
                                               double eikonal_q = std::numeric_limits<double>::quiet_NaN()) {
  const std::string path = WriteSoftModelFile(theta, suffix, odderon_coupling, odderon_sign, odderon_transition,
                                              pomeron_kappa, unitarization, eikonal_q);
  gra::MODELPARAM        = "TUNE0";
  std::filesystem::create_directories(gra::aux::GetBasePath(2) + "/eikonal");

  const auto    model_tune = gra::MModelTune::Load(path);
  gra::MEikonal eikonal(model_tune);
  if (number_loop_kt > 0 || number_loop_phi > 0) {
    const int extra_loop_kt  = number_loop_kt > 0 ? number_loop_kt - eikonal.Numerics.LOOP.radial_intervals : 0;
    const int extra_loop_phi = number_loop_phi > 0 ? number_loop_phi - eikonal.Numerics.LOOP.azimuth_nodes : 0;
    eikonal.Numerics.SetLoopDiscretization(extra_loop_kt, extra_loop_phi);
  }
  eikonal.S3Constructor(mandelstam_s, initialstate, false, number_bt, number_kt2);
  return eikonal;
}

// Build a small proton-proton eikonal object for tests
gra::MEikonal BuildTestEikonal(const std::vector<double> &theta, const std::string &suffix, int number_bt = 4,
                               int number_kt2 = 4, double pomeron_kappa = 0.0, const std::string &unitarization = {},
                               double eikonal_q = std::numeric_limits<double>::quiet_NaN()) {
  return BuildTestEikonalWithInitialState(theta, suffix, ProtonInitialState(), 0.25, 1, number_bt, number_kt2, 25.0,
                                          0.0, pomeron_kappa, 0, 0, unitarization, eikonal_q);
}

// Build the spin-averaged exclusive amplitude from the public helicity bank
std::complex<double> ExclusiveAmpReference(const gra::MEikonal &eikonal, double kt2, std::size_t f1, std::size_t f2) {
  const auto amplitude = eikonal.GetMatrixRuntime().HelicityAmplitudes(kt2, f1, f2);
  return 0.5 * (amplitude.phi1 + amplitude.phi3);
}

// Bound roundoff from the direct D by D matrix projection
double ExclusiveProjectionMargin(const gra::MEikonal &eikonal, std::complex<double> reference) {
  const double dimension =
      static_cast<double>(eikonal.GetMatrixRuntime().ChannelCount() * eikonal.GetMatrixRuntime().ChannelCount());
  return 32.0 * std::numeric_limits<double>::epsilon() * dimension * dimension * std::max(1.0, std::abs(reference));
}

MMatrix<double> InclusiveFromExclusive(const MMatrix<double> &exclusive) {
  MMatrix<double> inclusive(2, 2, 0.0);
  inclusive[0][0] = exclusive[0][0];
  for (std::size_t i = 1; i < exclusive.size_row(); ++i) { inclusive[1][0] += exclusive[i][0]; }
  for (std::size_t j = 1; j < exclusive.size_col(); ++j) { inclusive[0][1] += exclusive[0][j]; }
  for (std::size_t i = 1; i < exclusive.size_row(); ++i) {
    for (std::size_t j = 1; j < exclusive.size_col(); ++j) { inclusive[1][1] += exclusive[i][j]; }
  }
  return inclusive;
}

// Prepare an on-shell Born state through the production invariant calculation
void PrepareScreeningPoint(gra::MProcessState &state) {
  auto &lts = state.lts;
  lts.pfinal.resize(lts.decaytree.size() + 3);
  for (const auto &i : indices(lts.decaytree)) { lts.pfinal[i + 3] = lts.decaytree[i].p4; }
  REQUIRE(gra::math::CheckEMC(lts.pbeam1 + lts.pbeam2 - lts.pfinal[0] - lts.pfinal[1] - lts.pfinal[2]));
  REQUIRE(lts.pbeam1.M2() == Approx(pow2(lts.beam1.mass)).epsilon(1e-7).margin(1e-9));
  REQUIRE(lts.pbeam2.M2() == Approx(pow2(lts.beam2.mass)).epsilon(1e-7).margin(1e-9));
  REQUIRE(lts.pfinal[1].M2() == Approx(pow2(lts.beam1.mass)).epsilon(1e-7).margin(1e-9));
  REQUIRE(lts.pfinal[2].M2() == Approx(pow2(lts.beam2.mass)).epsilon(1e-7).margin(1e-9));
  REQUIRE(gra::kinematics::SetLorentzScalars(state, lts.decaytree.size() + 2));
}

// Expose production mass sampling and numerical validation without replacing physics
class ProcessProbe : public gra::MFactorized {
 public:
  using gra::MProcess::BookkeepAmplitudeWeight;
  using gra::MProcess::EvaluateAmplitudeBoundary;
  using gra::MProcess::ScreenedAmplitudeSquared;
  using gra::MProcess::ValidateAmplitude;

  // Sample the real process mass proposal
  double SampleBranchMass(gra::MDecayBranch &branch) {
    double mass = 0.0;
    GetOffShellMass(branch, mass);
    return mass;
  }
};

enum class ToyPhotonTopology { Upper, Lower, Coherent, Hadronic };

// Install an on-shell charged-pion pair with a nontrivial decay azimuth
void ConfigureToyPionPair(gra::LORENTZSCALAR &lts) {
  lts.decaytree[0].p          = lts.PDG.FindByPDG(211);
  lts.decaytree[1].p          = lts.PDG.FindByPDG(-211);
  const double     mass       = lts.decaytree[0].p.mass;
  const double     momentum   = gra::kinematics::DecayMomentum(lts.pfinal[0].M(), mass, mass);
  const double     theta      = 0.91;
  const double     phi        = 0.37;
  const double     transverse = momentum * std::sin(theta);
  const double     energy     = std::sqrt(momentum * momentum + mass * mass);
  const gra::M4Vec first_rest(transverse * std::cos(phi), transverse * std::sin(phi), momentum * std::cos(theta),
                              energy);
  const gra::M4Vec second_rest(-first_rest.Px(), -first_rest.Py(), -first_rest.Pz(), energy);
  lts.decaytree[0].p4 = BoostFromRestFrame(first_rest, lts.pfinal[0]);
  lts.decaytree[1].p4 = BoostFromRestFrame(second_rest, lts.pfinal[0]);
  RefreshToyDerivedKinematicsPreserveDecay(lts);
}

// Configure one photon-Pomeron resonance or Pomeron-pair continuum route
void ConfigureToyReggePhaseRoute(gra::LORENTZSCALAR &lts, gra::ReggeProductionModel spin, gra::MReggeMode mode,
                                 ToyPhotonTopology topology) {
  if (spin == gra::ReggeProductionModel::MP) { lts.process.MP_FRAME = "CM"; }
  lts.process.SPINGEN        = true;
  lts.process.FORWARD_NOFLIP = true;
  lts.process.FORWARD_VERTEX = gra::ForwardVertexMode::HelicityResidue;
  lts.process.PHOTON_VERTEX  = "EPA";
  lts.process.MMAX           = 1;
  lts.process.TU_SIGN        = "positive";
  lts.process.RESONANCES.clear();
  ConfigureToyPionPair(lts);

  if (mode == gra::MReggeMode::Resonance) {
    gra::PARAM_RES resonance = spin == gra::ReggeProductionModel::GP
                                   ? MakeToyPhotoGPResonance(lts.process.MMAX)
                                   : (spin == gra::ReggeProductionModel::XP ? MakeToyCovariantScalarPomeronPhotoXP()
                                                                            : MakeToyScalarPomeronPhotoXP());
    resonance.p              = lts.PDG.FindByPDG(113);
    auto &form               = ReggeResonanceFormForTest(resonance, spin);
    form.ff_prod             = {gra::regge::FFType::Gaussian, gra::regge::FFNorm::Pole, {0.64}};
    resonance.BW             = gra::BreitWigner::FixedWidth;
    resonance.hel_decay      = SpinOneToScalarScalarHelicityMatrix({0.73, -0.19});
    resonance.hel_decay.ff_decay          = form.ff_prod;
    if (spin == gra::ReggeProductionModel::MP) {
      PrepareToyPoleOperators(resonance, gra::ReggeProductionModel::MP, {{{0, 2, 1.0}}}, true, true);
      resonance.spin_basis = "none";
    }

    gra::RES_PRODUCTION lower = resonance.production.front();
    std::swap(lower.tree[0], lower.tree[1]);
    if (spin == gra::ReggeProductionModel::GP) {
      lower.hel.T     = lower.hel.T.Transpose();
      lower.hel.T_set = lower.hel.T_set.Transpose();
    } else if (spin == gra::ReggeProductionModel::XP) {
      lower.pole =
          gra::spin::PreparePoleLS(resonance.p, lower.tree[0].p, lower.tree[1].p, {{0, 2, 1.0}}, 1.0, true, true, true);
      lower.hel  = lower.pole->helicity;
    } else {
      lower.pole =
          gra::spin::PreparePoleLS(resonance.p, lower.tree[0].p, lower.tree[1].p, {{0, 2, 1.0}}, 1.0, true, true, true);
      lower.hel  = lower.pole->helicity;
    }

    if (topology == ToyPhotonTopology::Upper) {
      resonance.production.resize(1);
    } else if (topology == ToyPhotonTopology::Lower) {
      resonance.production = {std::move(lower)};
    } else if (topology == ToyPhotonTopology::Coherent) {
      resonance.production.push_back(std::move(lower));
    }
    lts.process.RESONANCES = {{"rho_phase_test", std::move(resonance)}};
    return;
  }

  REQUIRE(mode == gra::MReggeMode::ContinuumTwoBody);
  const bool        analytic    = spin == gra::ReggeProductionModel::GP;
  const int         pomeron_pdg = analytic ? 990 : 993;
  gra::MDecayBranch upper;
  upper.p                         = lts.PDG.FindByPDG(pomeron_pdg);
  upper.hel                       = RealisticProtonLegHelicityMatrix(pomeron_pdg, 2, -1, 1, 2);
  gra::MDecayBranch lower         = upper;
  lts.process.CONT_PRODUCTION     = {{pomeron_pdg, pomeron_pdg}};
  lts.process.CONT_PRODUCTIONTREE = {{upper, lower}};
  if (topology != ToyPhotonTopology::Hadronic) {
    gra::MDecayBranch photon;
    photon.p   = lts.PDG.FindByPDG(22);
    photon.hel = RealisticProtonLegHelicityMatrix(22, 2, -1, 1, 2);
    if (topology == ToyPhotonTopology::Upper) {
      lts.process.CONT_PRODUCTION     = {{22, pomeron_pdg}};
      lts.process.CONT_PRODUCTIONTREE = {{photon, lower}};
    } else if (topology == ToyPhotonTopology::Lower) {
      lts.process.CONT_PRODUCTION     = {{pomeron_pdg, 22}};
      lts.process.CONT_PRODUCTIONTREE = {{upper, photon}};
    } else {
      lts.process.CONT_PRODUCTION     = {{22, pomeron_pdg}, {pomeron_pdg, 22}};
      lts.process.CONT_PRODUCTIONTREE = {{photon, lower}, {upper, photon}};
    }
  }
  lts.process.CONT_TU_SIGN.clear();
  if (spin == gra::ReggeProductionModel::MP) {
    PrepareToyContinuumOperators(lts, gra::ReggeProductionModel::MP);
  } else if (spin == gra::ReggeProductionModel::XP) {
    PrepareToyContinuumOperators(lts, gra::ReggeProductionModel::XP);
  } else if (spin == gra::ReggeProductionModel::GP) {
    UseToyGPSubchannelHelicity(lts);
  }
}

// Run one real Regge route through the production and LOOPSCREEN interfaces
class ToyReggePhaseScreeningProcess : public gra::MFactorized {
 public:
  ToyReggePhaseScreeningProcess(const gra::MModelTunePtr &model_tune, gra::ReggeProductionModel spin,
                                gra::MReggeMode mode, const std::string &frame = "CM", double azimuth_rotation = 0.0,
                                ToyPhotonTopology                      topology      = ToyPhotonTopology::Coherent,
                                const std::vector<std::array<int, 2>> &hard_channels = {})
      : spin_(spin), mode_(mode) {
    state.screening = true;
    ProcPtr.ISTATE =
        spin == gra::ReggeProductionModel::MP ? "MP" : (spin == gra::ReggeProductionModel::XP ? "XP" : "GP");
    ProcPtr.CHANNEL = mode == gra::MReggeMode::Resonance ? "RES" : "CON";
    state.lts       = MakeToyCoherentPhotonLTS();
    ConfigureToyReggePhaseRoute(state.lts, spin_, mode_, topology);
    if (!hard_channels.empty()) {
      REQUIRE(spin_ == gra::ReggeProductionModel::MP);
      REQUIRE(mode_ == gra::MReggeMode::ContinuumTwoBody);
      REQUIRE(topology == ToyPhotonTopology::Hadronic);
      state.lts.process.CONT_PRODUCTION.clear();
      state.lts.process.CONT_PRODUCTIONTREE.clear();
      for (const auto &channel : hard_channels) {
        gra::MDecayBranch upper;
        upper.p   = state.lts.PDG.FindByPDG(channel[0]);
        upper.hel = ToyProtonLegHelicityMatrix(upper.p.spinX2);
        gra::MDecayBranch lower;
        lower.p   = state.lts.PDG.FindByPDG(channel[1]);
        lower.hel = ToyProtonLegHelicityMatrix(lower.p.spinX2);
        state.lts.process.CONT_PRODUCTION.push_back({channel[0], channel[1]});
        state.lts.process.CONT_PRODUCTIONTREE.push_back({std::move(upper), std::move(lower)});
      }
      PrepareToyContinuumOperators(state.lts, gra::ReggeProductionModel::MP);
    }
    if (spin_ == gra::ReggeProductionModel::MP) { state.lts.process.MP_FRAME = frame; }
    if (!gra::math::IsZero(azimuth_rotation)) { state.lts = RotateToyEventAroundZ(state.lts, azimuth_rotation); }

    gra::ScreeningMetadata metadata;
    metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonHelicity;
    metadata.proton_mode             = gra::ProtonScreeningMode::ForwardExcitation;
    metadata.amplitude_normalization = 0.25;
    metadata.spin_rows               = 4;
    metadata.forward_noflip          = true;
    metadata.amplitude_type          = gra::ScreeningAmplitudeType::GoodWalker;
    metadata.PrepareSpinTransitions();
    state.lts.hamp.Configure(metadata);
    SetModelTune(model_tune);
    PrepareScreeningPoint(state);
    regge_ = std::make_unique<gra::MRegge>(state.lts, model_tune,
                                           gra::MRegge::ProcessDefinitionFor(mode_, "phase_screening_test"));
  }

  // Evaluate only the hard Born amplitude
  double BornAmp2() {
    pair_trace_.clear();
    state.lts.proton_good_walker.reset();
    return ScreenedAmplitudeSquared(false);
  }

  // Evaluate the hard amplitude and its eikonal convolution
  double ScreenedAmp2() {
    pair_trace_.clear();
    state.lts.proton_good_walker.reset();
    return ScreenedAmplitudeSquared(true);
  }

  // Compute the Born followed by every shifted Regge-pair source
  const std::vector<gra::ProtonGoodWalkerAmplitude> &PairTrace() const { return pair_trace_; }

 protected:
  // Evaluate and retain one pair-space source for the independent sum
  double EvaluateBareAmplitude() override {
    const double amp2 = regge_->Amp2(state.lts, spin_, mode_);
    REQUIRE(state.lts.proton_good_walker.has_value());
    pair_trace_.push_back(*state.lts.proton_good_walker);
    return amp2;
  }

 private:
  gra::ReggeProductionModel                   spin_;
  gra::MReggeMode                             mode_;
  std::unique_ptr<gra::MRegge>                regge_;
  std::vector<gra::ProtonGoodWalkerAmplitude> pair_trace_;
};

// Hold the exact Born, interference and loop-squared decomposition
struct ToyReggeScreeningReference {
  std::vector<std::complex<double>> amplitude;
  double                            born         = 0.0;
  double                            interference = 0.0;
  double                            loop_squared = 0.0;
};

// Contract traced shifted amplitudes with the public coupled-channel matrix
ToyReggeScreeningReference ManualReggeScreeningFromTrace(const std::vector<gra::ProtonGoodWalkerAmplitude> &trace,
                                                         const gra::MEikonal                               &eikonal) {
  REQUIRE_FALSE(trace.empty());
  REQUIRE(trace.front().model == eikonal.SoftModelHandle());
  REQUIRE(trace.front().channel_count == eikonal.GetChannelCount());
  REQUIRE(trace.front().components.size() == 1);
  gra::ScreeningMetadata metadata;
  metadata.spin_basis              = gra::ScreeningSpinBasis::ProtonHelicity;
  metadata.spin_rows               = 4;
  metadata.forward_noflip          = true;
  metadata.amplitude_normalization = 0.25;
  metadata.PrepareSpinTransitions();
  const auto &component   = trace.front().components.front();
  const auto &born_source = component.source;
  REQUIRE(metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
  REQUIRE(metadata.spin_rows == 4);
  REQUIRE(metadata.forward_noflip);
  const std::size_t pair_dimension = trace.front().model->GoodWalker().PairDimension();
  REQUIRE(born_source.size_col() == pair_dimension);
  REQUIRE(born_source.size_row() % metadata.spin_rows == 0);
  const std::size_t spectators = born_source.size_row() / metadata.spin_rows;
  const auto       &loop       = eikonal.GetLoopConst(eikonal.InitializedMandelstamS());
  const std::size_t nodes      = loop.kt2.size() * loop.node_weight.size_col();
  REQUIRE(trace.size() == nodes + 1);

  std::vector<std::complex<double>> born(16 * spectators * pair_dimension, 0.0);
  std::vector<std::complex<double>> correction(16 * spectators * pair_dimension, 0.0);
  for (std::size_t index = 0; index < metadata.spin_transition_count; ++index) {
    const auto &transition = metadata.spin_transition[index];
    for (std::size_t spectator = 0; spectator < spectators; ++spectator) {
      const std::size_t source = transition.source_row * spectators + spectator;
      const std::size_t destination =
          (gra::spin::PairHelicityTransitionIndex(transition.initial, transition.intermediate) * spectators +
           spectator) *
          pair_dimension;
      for (std::size_t pair = 0; pair < pair_dimension; ++pair) {
        born[destination + pair] += born_source[source][pair];
      }
    }
  }

  constexpr auto    helicity_transition = gra::CanonicalProtonHelicityTransitions();
  const std::size_t azimuth_count       = loop.node_weight.size_col();
  for (std::size_t node = 0; node < nodes; ++node) {
    const std::size_t radial  = node / azimuth_count;
    const std::size_t azimuth = node % azimuth_count;
    const auto &soft = loop.pair_screening_spin.at(radial);
    REQUIRE(trace[node + 1].model == trace.front().model);
    REQUIRE(trace[node + 1].channel_count == trace.front().channel_count);
    REQUIRE(trace[node + 1].components.size() == 1);
    REQUIRE(trace[node + 1].components.front().coherence_group == component.coherence_group);
    REQUIRE(trace[node + 1].components.front().upper_sector == component.upper_sector);
    REQUIRE(trace[node + 1].components.front().lower_sector == component.lower_sector);
    const auto &hard = trace[node + 1].components.front().source;
    REQUIRE(hard.size_row() == born_source.size_row());
    REQUIRE(hard.size_col() == pair_dimension);
    for (std::size_t index = 0; index < metadata.spin_transition_count; ++index) {
      const auto &transition = metadata.spin_transition[index];
      for (std::size_t final = 0; final < 4; ++final) {
        const std::size_t spin_index = gra::spin::PairHelicityMatrixIndex(final, transition.intermediate);
        REQUIRE(soft[spin_index].size_row() == pair_dimension);
        REQUIRE(soft[spin_index].size_col() == pair_dimension);
        const double               phi      = std::atan2(loop.kt_y(radial, azimuth), loop.kt_x(radial, azimuth));
        const std::complex<double> rotation = std::polar(1.0, helicity_transition[spin_index].azimuth_harmonic * phi);
        const std::complex<double> weight   = loop.node_weight(radial, azimuth) * rotation;
        for (std::size_t spectator = 0; spectator < spectators; ++spectator) {
          const std::size_t source = transition.source_row * spectators + spectator;
          const std::size_t destination =
              (gra::spin::PairHelicityTransitionIndex(transition.initial, final) * spectators + spectator) *
              pair_dimension;
          for (std::size_t a = 0; a < pair_dimension; ++a) {
            for (std::size_t b = 0; b < pair_dimension; ++b) {
              correction[destination + a] += weight * soft[spin_index][a][b] * hard[source][b];
            }
          }
        }
      }
    }
  }

  const auto project = [&](const auto &pair_source) {
    std::vector<std::complex<double>> physical;
    const std::size_t                 blocks = pair_source.size() / pair_dimension;
    for (std::size_t block = 0; block < blocks; ++block) {
      const std::span<const std::complex<double>> source(pair_source.data() + block * pair_dimension, pair_dimension);
      const auto                                  projected = trace.front().model->GoodWalker().ProjectPair(
                                           source, gra::SectorFinalBasis(component.upper_sector), gra::SectorFinalBasis(component.lower_sector));
      physical.insert(physical.end(), projected.cbegin(), projected.cend());
    }
    return physical;
  };
  const auto born_physical       = project(born);
  const auto correction_physical = project(correction);
  REQUIRE(born_physical.size() == correction_physical.size());

  ToyReggeScreeningReference out;
  out.amplitude = born_physical;
  std::transform(out.amplitude.cbegin(), out.amplitude.cend(), correction_physical.cbegin(), out.amplitude.begin(),
                 std::plus<>{});
  for (const auto &i : indices(born_physical)) {
    out.born += std::norm(born_physical[i]);
    out.interference += 2.0 * std::real(std::conj(born_physical[i]) * correction_physical[i]);
    out.loop_squared += std::norm(correction_physical[i]);
  }
  out.born *= metadata.amplitude_normalization;
  out.interference *= metadata.amplitude_normalization;
  out.loop_squared *= metadata.amplitude_normalization;
  return out;
}

// Contract a full toy proton-helicity Born matrix with the cached loop weights
double ExpectedToyHelicityScreenedAmpSquared(const gra::MEikonal                     &eikonal,
                                             const std::vector<std::complex<double>> &born) {
  REQUIRE(born.size() == 16);
  const auto &loop_const  = eikonal.GetLoopConst(25.0);
  double      amp_squared = 0.0;
  for (std::size_t initial = 0; initial < 4; ++initial) {
    for (std::size_t final = 0; final < 4; ++final) {
      std::complex<double> amplitude = born[gra::spin::PairHelicityTransitionIndex(initial, final)];
      for (std::size_t node = 0; node < loop_const.physical_screening_helicity_weight.size(); ++node) {
        const auto &weight = loop_const.physical_screening_helicity_weight[node];
        for (std::size_t intermediate = 0; intermediate < 4; ++intermediate) {
          amplitude += weight[gra::spin::PairHelicityMatrixIndex(final, intermediate)] *
                       born[gra::spin::PairHelicityTransitionIndex(initial, intermediate)];
        }
      }
      amp_squared += gra::math::abs2(amplitude);
    }
  }
  return 0.25 * amp_squared;
}

// Contract one spin-independent toy Born amplitude with a cached GW class
double ExpectedToyGoodWalkerAmpSquared(const gra::MEikonal &eikonal, const std::complex<double> born,
                                       const gra::GWFinalClass final_class) {
  const auto &loop_const  = eikonal.GetLoopConst(25.0);
  const auto &channels    = loop_const.good_walker_channels[static_cast<std::size_t>(final_class)];
  double      amp_squared = 0.0;
  for (const auto &channel : channels) {
    REQUIRE(channel.spin_scalar);
    for (std::size_t initial = 0; initial < 4; ++initial) {
      for (std::size_t final = 0; final < 4; ++final) {
        std::complex<double> amplitude = final == initial ? channel.born_coefficient * born : 0.0;
        for (const auto weight : channel.scalar_screening_weight) {
          if (final == initial) { amplitude += weight * born; }
        }
        amp_squared += gra::math::abs2(amplitude);
      }
    }
  }
  return 0.25 * amp_squared;
}

void UpdateToyDurhamDerivedKinematics(gra::LORENTZSCALAR &lts) {
  lts.pfinal[0] = gra::M4Vec();
  for (const auto &particle : lts.decaytree) { lts.pfinal[0] += particle.p4; }
  lts.q1      = lts.pbeam1 - lts.pfinal[1];
  lts.q2      = lts.pbeam2 - lts.pfinal[2];
  lts.t1      = lts.q1.M2();
  lts.t2      = lts.q2.M2();
  lts.s       = (lts.pbeam1 + lts.pbeam2).M2();
  lts.sqrt_s  = std::sqrt(std::max(0.0, lts.s));
  lts.m2      = lts.pfinal[0].M2();
  lts.s_hat   = lts.m2;
  lts.Y       = lts.pfinal[0].Rap();
  lts.Pt      = lts.pfinal[0].Pt();
  lts.s1      = (lts.pfinal[0] + lts.pfinal[1]).M2();
  lts.s2      = (lts.pfinal[0] + lts.pfinal[2]).M2();
  lts.x1      = 1.0 - lts.pfinal[1].Pz() / lts.pbeam1.Pz();
  lts.x2      = 1.0 - lts.pfinal[2].Pz() / lts.pbeam2.Pz();
  lts.xi1     = lts.x1;
  lts.xi2     = lts.x2;
  lts.has_xi1 = true;
  lts.has_xi2 = true;
  lts.qt1     = lts.q1.Pt();
  lts.qt2     = lts.q2.Pt();
  lts.q1_in_X = BoostToRestFrame(lts.q1, lts.pfinal[0]);
  lts.q2_in_X = BoostToRestFrame(lts.q2, lts.pfinal[0]);
  if (lts.decaytree.size() >= 2) {
    lts.t_hat   = (lts.q1 - lts.decaytree[0].p4).M2();
    lts.u_hat   = (lts.q1 - lts.decaytree[1].p4).M2();
    lts.d0_in_X = BoostToRestFrame(lts.decaytree[0].p4, lts.pfinal[0]);
    lts.d1_in_X = BoostToRestFrame(lts.decaytree[1].p4, lts.pfinal[0]);
  }
}

gra::LORENTZSCALAR MakeToyDurhamGG() {
  gra::LORENTZSCALAR lts;
  lts.model_cache  = std::make_shared<gra::MModelCache>();
  lts.LHAPDFSET    = "MMHT2014lo68cl";
  const auto beams = ProtonInitialState();
  lts.beam1        = beams[0];
  lts.beam2        = beams[1];

  const double beam_pz = 10.0;
  const double beam_e  = std::sqrt(pow2(beam_pz) + pow2(gra::PDG::mp));
  lts.pbeam1           = gra::M4Vec(0.0, 0.0, beam_pz, beam_e);
  lts.pbeam2           = gra::M4Vec(0.0, 0.0, -beam_pz, beam_e);

  lts.pfinal.resize(3);
  const double p1x = 0.18;
  const double p1y = -0.05;
  const double p1z = 8.15;
  const double p2x = -0.11;
  const double p2y = 0.04;
  const double p2z = -7.55;
  lts.pfinal[1]    = gra::M4Vec(p1x, p1y, p1z, std::sqrt(pow2(p1x) + pow2(p1y) + pow2(p1z) + pow2(gra::PDG::mp)));
  lts.pfinal[2]    = gra::M4Vec(p2x, p2y, p2z, std::sqrt(pow2(p2x) + pow2(p2y) + pow2(p2z) + pow2(gra::PDG::mp)));
  lts.pfinal[0]    = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];

  gra::MParticle gluon;
  gluon.name   = "g";
  gluon.pdg    = 21;
  gluon.spinX2 = 2;
  gluon.color  = 8;
  gluon.mass   = 0.0;

  lts.decaytree.resize(2);
  lts.decaytree[0].p = gluon;
  lts.decaytree[1].p = gluon;

  const double mass     = lts.pfinal[0].M();
  const double costheta = 0.31;
  const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
  const double phi      = 0.57;
  const double momentum = mass / 2.0;

  gra::M4Vec p3(momentum * sintheta * std::cos(phi), momentum * sintheta * std::sin(phi), momentum * costheta,
                momentum);
  gra::M4Vec p4(-p3.Px(), -p3.Py(), -p3.Pz(), momentum);
  gra::kinematics::LorentzBoost(lts.pfinal[0], mass, p3, +1);
  gra::kinematics::LorentzBoost(lts.pfinal[0], mass, p4, +1);

  lts.decaytree[0].p4 = p3;
  lts.decaytree[1].p4 = p4;
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

gra::LORENTZSCALAR MakeToyDurhamQQbar(int quark_pdg) {
  auto lts = MakeToyDurhamGG();

  gra::MParticle quark;
  quark.pdg    = quark_pdg;
  quark.spinX2 = 1;
  quark.color  = 3;
  quark.mass   = 0.0;
  quark.name   = (std::abs(quark_pdg) == 2) ? "u" : "d";

  gra::MParticle antiquark = quark;
  antiquark.pdg            = -quark_pdg;
  antiquark.name += "~";

  lts.decaytree[0].p = quark;
  lts.decaytree[1].p = antiquark;
  return lts;
}

// Build a high-mass direct PhotoZ event with elastic proton legs
gra::LORENTZSCALAR MakeToyPhotoZFFbar(int fermion_pdg, double central_mass = 91.1876) {
  gra::LORENTZSCALAR lts;
  lts.model_cache             = std::make_shared<gra::MModelCache>();
  lts.PDG                     = LoadedPDGTable();
  lts.LHAPDFSET               = "MMHT2014lo68cl";
  lts.process.root_decay_mode = gra::RootDecayMode::None;
  lts.beam1                   = lts.PDG.FindByPDG(2212);
  lts.beam2                   = lts.PDG.FindByPDG(2212);

  const double beam_e    = 200.0;
  const double beam_pz   = std::sqrt(pow2(beam_e) - pow2(gra::PDG::mp));
  const double proton_pt = 0.18;
  const double final_e   = beam_e - 0.5 * central_mass;
  const double final_pz  = std::sqrt(std::max(0.0, pow2(final_e) - pow2(gra::PDG::mp) - pow2(proton_pt)));

  lts.pbeam1 = gra::M4Vec(0.0, 0.0, beam_pz, beam_e);
  lts.pbeam2 = gra::M4Vec(0.0, 0.0, -beam_pz, beam_e);
  lts.pfinal.resize(3);
  lts.pfinal[1] = gra::M4Vec(proton_pt, 0.0, final_pz, final_e);
  lts.pfinal[2] = gra::M4Vec(-proton_pt, 0.0, -final_pz, final_e);

  const gra::MParticle fermion     = lts.PDG.FindByPDG(fermion_pdg);
  const gra::MParticle antifermion = lts.PDG.FindByPDG(-fermion_pdg);
  const double         mf          = std::max(0.0, fermion.mass);
  const double         energy      = 0.5 * central_mass;
  const double         momentum    = std::sqrt(std::max(0.0, pow2(energy) - pow2(mf)));
  const double         costheta    = 0.23;
  const double         sintheta    = std::sqrt(std::max(0.0, 1.0 - pow2(costheta)));
  const double         phi         = 0.41;

  lts.decaytree.resize(2);
  lts.decaytree[0].p = fermion;
  lts.decaytree[1].p = antifermion;
  lts.decaytree[0].p4 =
      gra::M4Vec(momentum * sintheta * std::cos(phi), momentum * sintheta * std::sin(phi), momentum * costheta, energy);
  lts.decaytree[1].p4 =
      gra::M4Vec(-lts.decaytree[0].p4.Px(), -lts.decaytree[0].p4.Py(), -lts.decaytree[0].p4.Pz(), energy);

  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

// Replace one elastic toy proton by a forward system of fixed invariant mass
void SetToyPhotoForwardExcitation(gra::LORENTZSCALAR &lts, int leg, double mass) {
  gra::M4Vec  &forward = lts.pfinal.at(static_cast<std::size_t>(leg));
  const double pz_sign = forward.Pz() >= 0.0 ? 1.0 : -1.0;
  const double pz2     = pow2(forward.E()) - pow2(mass) - forward.Pt2();
  if (!(pz2 > 0.0)) { throw std::invalid_argument("SetToyPhotoForwardExcitation: mass exceeds forward energy"); }
  forward.SetPzE(pz_sign * std::sqrt(pz2), forward.E());
  if (leg == 1) {
    lts.excite1 = true;
  } else if (leg == 2) {
    lts.excite2 = true;
  } else {
    throw std::invalid_argument("SetToyPhotoForwardExcitation: leg should be 1 or 2");
  }
  UpdateToyDurhamDerivedKinematics(lts);
}

// Read one HepMC integer attribute with a zero fallback
int ReadHepMCIntAttributeForTest(const HepMC3::ConstGenParticlePtr &particle, const std::string &name) {
  const auto attribute = particle->attribute<HepMC3::IntAttribute>(name);
  return attribute ? attribute->value() : 0;
}

class ToyPhotoZRecordProcess : public gra::MFactorized {
 public:
  // Build the common HepMC record from prepared test kinematics
  bool Record(HepMC3::GenEvent &evt) { return MProcess::BuildEventRecord(evt); }

};

class ToyPhotoZScreeningProcess : public gra::MFactorized {
 public:
  // Construct a prepared PhotoZ process for screening-loop checks
  explicit ToyPhotoZScreeningProcess(gra::MModelTunePtr model_tune) {
    state.screening = true;
    SetModelTune(std::move(model_tune));
    ProcPtr.CHANNEL = "Z";
    ProcPtr.ISTATE  = "ygg";
    state.lts       = MakeToyPhotoZFFbar(13);
    PrepareScreeningPoint(state);
    photoz          = std::make_unique<gra::MPhotoZ>(state.lts, state.model_tune, gra::MPhotoZ::ProcessDefinitionFor());
  }

  // Evaluate the bare PhotoZ amplitude through the protected hook
  double BareAmp2() { return EvaluateBareAmplitude(); }

  // Evaluate the screened PhotoZ amplitude through the existing loop machinery
  double ScreenedAmp2() { return ScreenedAmplitudeSquared(true); }

 protected:
  // Evaluate the PhotoZ matrix element for the current kinematics
  double EvaluateBareAmplitude() override { return photoz->Amp2(state.lts); }

 private:
  std::unique_ptr<gra::MPhotoZ> photoz;
};

// Build a two-body Durham meson-pair toy event with physical meson masses
gra::LORENTZSCALAR MakeToyDurhamMesonPair(int pdg0, int pdg1, double costheta = 0.23) {
  auto lts = MakeToyDurhamGG();

  // Raise the dedicated meson-pair toy mass above the default hard-pT cutoff
  const double p1z = 7.9;
  const double p2z = -7.3;
  lts.pfinal[1] =
      gra::M4Vec(lts.pfinal[1].Px(), lts.pfinal[1].Py(), p1z,
                 std::sqrt(pow2(lts.pfinal[1].Px()) + pow2(lts.pfinal[1].Py()) + pow2(p1z) + pow2(gra::PDG::mp)));
  lts.pfinal[2] =
      gra::M4Vec(lts.pfinal[2].Px(), lts.pfinal[2].Py(), p2z,
                 std::sqrt(pow2(lts.pfinal[2].Px()) + pow2(lts.pfinal[2].Py()) + pow2(p2z) + pow2(gra::PDG::mp)));
  lts.pfinal[0] = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];

  const auto          &pdg    = LoadedPDGTable();
  const gra::MParticle meson0 = pdg.FindByPDG(pdg0);
  const gra::MParticle meson1 = pdg.FindByPDG(pdg1);

  lts.decaytree.resize(2);
  lts.decaytree[0].p = meson0;
  lts.decaytree[1].p = meson1;

  const gra::M4Vec system = lts.pfinal[0];
  const double     mass   = system.M();
  const double     m0     = meson0.mass;
  const double     m1     = meson1.mass;
  const double     lambda = (pow2(mass) - pow2(m0 + m1)) * (pow2(mass) - pow2(m0 - m1));
  if (lambda <= 0.0) {
    throw std::invalid_argument("MakeToyDurhamMesonPair: central system is below meson-pair threshold");
  }

  const double momentum = std::sqrt(lambda) / (2.0 * mass);
  const double sintheta = std::sqrt(std::max(0.0, 1.0 - costheta * costheta));
  const double phi      = 0.41;

  gra::M4Vec p3(momentum * sintheta * std::cos(phi), momentum * sintheta * std::sin(phi), momentum * costheta,
                std::sqrt(pow2(momentum) + pow2(m0)));
  gra::M4Vec p4(-p3.Px(), -p3.Py(), -p3.Pz(), std::sqrt(pow2(momentum) + pow2(m1)));

  gra::kinematics::LorentzBoost(system, mass, p3, +1);
  gra::kinematics::LorentzBoost(system, mass, p4, +1);

  lts.decaytree[0].p4 = p3;
  lts.decaytree[1].p4 = p4;
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

gra::LORENTZSCALAR MakeToyDurhamGGG() {
  auto lts = MakeToyDurhamGG();

  gra::MParticle gluon;
  gluon.name   = "g";
  gluon.pdg    = 21;
  gluon.spinX2 = 2;
  gluon.color  = 8;
  gluon.mass   = 0.0;

  lts.decaytree.resize(3);
  for (auto &particle : lts.decaytree) { particle.p = gluon; }

  const double mass     = lts.pfinal[0].M();
  const double momentum = mass / 3.0;
  const double cosphi   = -0.5;
  const double sinphi   = std::sqrt(3.0) / 2.0;

  gra::M4Vec p3(momentum, 0.0, 0.0, momentum);
  gra::M4Vec p4(momentum * cosphi, momentum * sinphi, 0.0, momentum);
  gra::M4Vec p5(momentum * cosphi, -momentum * sinphi, 0.0, momentum);

  gra::kinematics::LorentzBoost(lts.pfinal[0], mass, p3, +1);
  gra::kinematics::LorentzBoost(lts.pfinal[0], mass, p4, +1);
  gra::kinematics::LorentzBoost(lts.pfinal[0], mass, p5, +1);

  lts.decaytree[0].p4 = p3;
  lts.decaytree[1].p4 = p4;
  lts.decaytree[2].p4 = p5;
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

// Build a massless three-parton state with one controlled azimuthal separation
gra::LORENTZSCALAR MakeToyDurhamGGGPair(double separation) {
  auto             lts      = MakeToyDurhamGGG();
  const gra::M4Vec system   = lts.pfinal[0];
  const double     mass     = system.M();
  const double     half     = 0.5 * separation;
  const double     momentum = mass / (2.0 * (1.0 + std::cos(half)));

  std::array<gra::M4Vec, 3> partons = {
      gra::M4Vec(momentum * std::cos(half), -momentum * std::sin(half), 0.0, momentum),
      gra::M4Vec(momentum * std::cos(half), momentum * std::sin(half), 0.0, momentum),
      gra::M4Vec(-2.0 * momentum * std::cos(half), 0.0, 0.0, 2.0 * momentum * std::cos(half)),
  };
  for (std::size_t i = 0; i < partons.size(); ++i) {
    gra::kinematics::LorentzBoost(system, mass, partons[i], +1);
    lts.decaytree[i].p4 = partons[i];
  }
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

// Build a massless four-gluon Durham reference point
gra::LORENTZSCALAR MakeToyDurhamGGGG() {
  auto                 lts   = MakeToyDurhamGGG();
  const gra::MParticle gluon = lts.decaytree.front().p;
  lts.decaytree.resize(4);
  for (auto &particle : lts.decaytree) { particle.p = gluon; }

  const double                               mass       = lts.pfinal[0].M();
  const double                               momentum   = mass / 4.0;
  const double                               unit       = 1.0 / std::sqrt(3.0);
  const std::array<std::array<double, 3>, 4> directions = {
      {{{unit, unit, unit}}, {{unit, -unit, -unit}}, {{-unit, unit, -unit}}, {{-unit, -unit, unit}}}};

  for (std::size_t i = 0; i < directions.size(); ++i) {
    gra::M4Vec momentum_i(momentum * directions[i][0], momentum * directions[i][1], momentum * directions[i][2],
                          momentum);
    gra::kinematics::LorentzBoost(lts.pfinal[0], mass, momentum_i, +1);
    lts.decaytree[i].p4 = momentum_i;
  }
  UpdateToyDurhamDerivedKinematics(lts);
  return lts;
}

// Build a massless three-parton Durham q q~ g reference point
gra::LORENTZSCALAR MakeToyDurhamQQbarG(int quark_pdg) {
  auto lts = MakeToyDurhamGGG();

  gra::MParticle quark;
  quark.pdg    = quark_pdg;
  quark.spinX2 = 1;
  quark.color  = 3;
  quark.mass   = 0.0;
  quark.name   = (quark_pdg == gra::PDG::PDG_hard_jet) ? "j" : "q";

  gra::MParticle antiquark = quark;
  antiquark.pdg            = -quark_pdg;
  antiquark.name += "~";

  lts.decaytree[0].p = quark;
  lts.decaytree[1].p = antiquark;
  return lts;
}

// Evaluate one toy event through the automated Durham MG5 registry
void EvaluateGeneratedDurham(gra::MDurham &durham, gra::LORENTZSCALAR &lts, gra::MDurham::DurhamProjectedAmp &projected,
                             std::vector<gra::MDurham::DurhamProjectedAmp> *flow_projected = nullptr) {
  gra::amplitude::ProcessRegistry registry(gra::amplitude::Processes("DURHAM"));
  const auto                      process = registry.MatchProcess(lts.decaytree);
  if (!process.has_value()) {
    throw std::invalid_argument("EvaluateGeneratedDurham: no generated process for toy event");
  }
  auto matrix_element = gra::CreateDurhamMG5Process(*process);
  if (matrix_element == nullptr) {
    throw std::invalid_argument("EvaluateGeneratedDurham: generated process construction failed");
  }
  durham.Dgg2Generated(lts, *matrix_element, projected, flow_projected);
}

// Compute tensor-Pomeron production couplings from the active test tune
std::array<double, 2> TensorVectorCouplingsForTest(int pdg) {
  const auto  tune      = gra::MModelTune::Load(modelfile);
  const auto  params    = gra::ReadTensorPomeronParam(*tune, LoadedPDGTable());
  const auto &couplings = params->exchange.FindVertex(995, pdg).g_tensor;
  REQUIRE(couplings.size() == 2);
  return {couplings[0], couplings[1]};
}

// Compute the channel specific Pomeron vector vector transfer form factor
gra::regge::FFParam TensorVectorTransferForTest(int pdg) {
  const auto tune   = gra::MModelTune::Load(modelfile);
  const auto params = gra::ReadTensorPomeronParam(*tune, LoadedPDGTable());
  return params->FindVector(pdg).ff_transfer;
}

// Compute the channel specific exchanged vector subchannel form factor
double TensorVectorOffShellForTest(int pdg, double q2, double mass) {
  const auto  tune   = gra::MModelTune::Load(modelfile);
  const auto  params = gra::ReadTensorPomeronParam(*tune, LoadedPDGTable());
  const auto &vector = params->FindVector(pdg);
  return gra::regge::FormFactor(q2, gra::math::pow2(mass), vector.ff_offshell);
}

// Compute the default relativistic BW mass-squared proposal normalization
double TensorCascadeMassProposalNormForTest(const gra::MParticle &p, double daughter_mass_sum, double offshell = 5.0) {
  const double lower = std::max(daughter_mass_sum, std::max(0.0, p.mass - offshell * p.width));
  const double upper = p.mass + offshell * p.width;
  return gra::MRandom::RelativisticBWMass2Integral(p.mass, p.width, gra::math::pow2(lower), gra::math::pow2(upper));
}

// Build one vector branch with deterministic off-shell daughters
gra::MDecayBranch TensorVectorBranchForTest(int vector_pdg, double pole_mass, double width, const gra::M4Vec &vector_p4,
                                            int d1_pdg, int d2_pdg, double daughter_mass, double theta, double phi,
                                            double tensor_decay_coupling) {
  gra::MDecayBranch branch;
  branch.name                = (vector_pdg == 333) ? "phi" : "rho";
  branch.p                   = ToyParticle(branch.name, vector_pdg, 2, pole_mass);
  branch.p.P                 = -1;
  branch.p.C                 = -1;
  branch.p.width             = width;
  branch.p4                  = vector_p4;
  branch.mass_proposal_norm  = TensorCascadeMassProposalNormForTest(branch.p, 2.0 * daughter_mass);
  branch.mass_proposal_min2  = 0.0;
  branch.mass_proposal_max2  = 100.0;
  branch.mass_proposal = gra::MassProposal::BreitWigner;
  branch.hel.g_decay_TP  = {tensor_decay_coupling};

  const auto daughters_in_vector = TwoBodyRestKinematics(vector_p4.M(), daughter_mass, daughter_mass, theta, phi);

  gra::MDecayBranch d1;
  d1.name = (d1_pdg > 0) ? "d+" : "d-";
  d1.p    = ToyParticle(d1.name, d1_pdg, 0, daughter_mass);
  d1.p.P  = -1;
  d1.p4   = BoostFromRestFrame(daughters_in_vector[0], vector_p4);

  gra::MDecayBranch d2;
  d2.name = (d2_pdg > 0) ? "d+" : "d-";
  d2.p    = ToyParticle(d2.name, d2_pdg, 0, daughter_mass);
  d2.p.P  = -1;
  d2.p4   = BoostFromRestFrame(daughters_in_vector[1], vector_p4);

  branch.legs    = {d1, d2};
  branch.W_event = gra::kinematics::dPhi2(vector_p4.M(),
                                          gra::kinematics::DecayMomentum(vector_p4.M(), daughter_mass, daughter_mass));
  branch.W       = gra::kinematics::MCW(branch.W_event, gra::math::pow2(branch.W_event), 1.0);
  return branch;
}

// Build a physical tensor rho0rho0 cascade phase-space point
gra::LORENTZSCALAR TensorRhoCascadeLTSForTest() {
  gra::LORENTZSCALAR lts = MakeToyProductionLTSDphi(gra::math::PI);
  lts.PDG                = LoadedPDGTable();
  lts.process.SPINGEN    = true;
  lts.process.SPINDEC    = true;
  lts.PS_active          = true;
  lts.process.TENSOR_MODEL_READY = true;
  lts.process.TENSOR_VECTOR_DECAY_PDGS = {{113, 211}};

  const double     rho_pole     = 0.77526;
  const double     rho_width    = 0.1491;
  const double     m_pi         = 0.13957061;
  const double     m_rho_left   = 0.83;
  const double     m_rho_right  = 0.91;
  const auto       vectors_in_X = TwoBodyRestKinematics(lts.pfinal[0].M(), m_rho_left, m_rho_right, 0.91, -0.42);
  const gra::M4Vec rho_left     = BoostFromRestFrame(vectors_in_X[0], lts.pfinal[0]);
  const gra::M4Vec rho_right    = BoostFromRestFrame(vectors_in_X[1], lts.pfinal[0]);

  lts.decaytree = {TensorVectorBranchForTest(113, rho_pole, rho_width, rho_left, gra::PDG::PDG_pip, gra::PDG::PDG_pim,
                                             m_pi, 0.74, 0.51, 11.95),
                   TensorVectorBranchForTest(113, rho_pole, rho_width, rho_right, gra::PDG::PDG_pip, gra::PDG::PDG_pim,
                                             m_pi, 1.38, -0.81, 11.95)};

  const double root_ps = gra::kinematics::dPhi2(
      lts.pfinal[0].M(), gra::kinematics::DecayMomentum(lts.pfinal[0].M(), rho_left.M(), rho_right.M()));
  lts.DW = gra::kinematics::MCW(root_ps);
  return lts;
}

// Build a physical tensor phi phi cascade phase-space point
gra::LORENTZSCALAR TensorPhiCascadeLTSForTest() {
  gra::LORENTZSCALAR lts = MakeToyProductionLTSDphi(2.34, 0.17);
  lts.PDG                = LoadedPDGTable();
  lts.process.SPINGEN    = true;
  lts.process.SPINDEC    = true;
  lts.PS_active          = true;
  lts.process.TENSOR_MODEL_READY = true;
  lts.process.TENSOR_VECTOR_DECAY_PDGS = {{333, 321}};

  const double     phi_pole       = 1.019461;
  const double     phi_width      = 0.004249;
  const double     kaon_mass      = 0.493677;
  const double     phi_left_mass  = 1.035;
  const double     phi_right_mass = 1.045;
  const auto       vectors_in_X = TwoBodyRestKinematics(lts.pfinal[0].M(), phi_left_mass, phi_right_mass, 0.73, -0.29);
  const gra::M4Vec phi_left     = BoostFromRestFrame(vectors_in_X[0], lts.pfinal[0]);
  const gra::M4Vec phi_right    = BoostFromRestFrame(vectors_in_X[1], lts.pfinal[0]);

  lts.decaytree = {TensorVectorBranchForTest(333, phi_pole, phi_width, phi_left, gra::PDG::PDG_Kp, gra::PDG::PDG_Km,
                                             kaon_mass, 0.81, 0.37, 4.48),
                   TensorVectorBranchForTest(333, phi_pole, phi_width, phi_right, gra::PDG::PDG_Kp, gra::PDG::PDG_Km,
                                             kaon_mass, 1.12, -0.64, 4.48)};

  const double root_ps = gra::kinematics::dPhi2(
      lts.pfinal[0].M(), gra::kinematics::DecayMomentum(lts.pfinal[0].M(), phi_left.M(), phi_right.M()));
  lts.DW = gra::kinematics::MCW(root_ps);
  return lts;
}

// Build an exact two-body central state for Tensor-Pomeron or full-QED
// continuum tests
gra::LORENTZSCALAR DirectCentralPairLTSForTest(int left_pdg, int right_pdg) {
  gra::LORENTZSCALAR lts         = MakeToyProductionLTSDphi(0.73, 0.21);
  const double       beam_pz     = lts.pbeam1.Pz();
  const double       beam_energy = std::sqrt(pow2(beam_pz) + pow2(gra::PDG::mp));
  lts.pbeam1.SetPzE(beam_pz, beam_energy);
  lts.pbeam2.SetPzE(-beam_pz, beam_energy);
  lts.pfinal[0]       = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  lts.q1              = lts.pbeam1 - lts.pfinal[1];
  lts.q2              = lts.pbeam2 - lts.pfinal[2];
  lts.t1              = lts.q1.M2();
  lts.t2              = lts.q2.M2();
  lts.s               = (lts.pbeam1 + lts.pbeam2).M2();
  lts.sqrt_s          = std::sqrt(lts.s);
  lts.m2              = lts.pfinal[0].M2();
  lts.s1              = (lts.pfinal[0] + lts.pfinal[1]).M2();
  lts.s2              = (lts.pfinal[0] + lts.pfinal[2]).M2();
  lts.process.SPINGEN = true;
  lts.process.SPINDEC = true;

  gra::MDecayBranch left;
  left.p = lts.PDG.FindByPDG(left_pdg);
  gra::MDecayBranch right;
  right.p              = lts.PDG.FindByPDG(right_pdg);
  const auto pair_in_X = TwoBodyRestKinematics(lts.pfinal[0].M(), left.p.mass, right.p.mass, 0.91, -0.38);
  left.p4              = BoostFromRestFrame(pair_in_X[0], lts.pfinal[0]);
  right.p4             = BoostFromRestFrame(pair_in_X[1], lts.pfinal[0]);
  lts.decaytree        = {left, right};
  return lts;
}

// Build one physical pion-pair point at a selected pole and proton transfer
gra::LORENTZSCALAR ScalarPolePhasePointForTest(double proton_pt, double theta, double phi, double central_mass,
                                               double beam_pz = 5.0) {
  gra::LORENTZSCALAR lts           = MakeToyScalarContinuumLTSAsymmetric(0.2, 4.2, -4.2);
  const double       beam_energy   = std::sqrt(gra::math::pow2(beam_pz) + gra::math::pow2(gra::PDG::mp));
  const double       proton_energy = beam_energy - 0.5 * central_mass;
  const double       proton_pz =
      std::sqrt(gra::math::pow2(proton_energy) - gra::math::pow2(gra::PDG::mp) - gra::math::pow2(proton_pt));
  lts.pbeam1      = gra::M4Vec(0.0, 0.0, beam_pz, beam_energy);
  lts.pbeam2      = gra::M4Vec(0.0, 0.0, -beam_pz, beam_energy);
  lts.pfinal[1]   = gra::M4Vec(proton_pt, 0.0, proton_pz, proton_energy);
  lts.pfinal[2]   = gra::M4Vec(-proton_pt, 0.0, -proton_pz, proton_energy);
  lts.pfinal[0]   = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  const auto pair = TwoBodyRestKinematics(central_mass, lts.decaytree[0].p.mass, lts.decaytree[1].p.mass, theta, phi);
  lts.decaytree[0].p4        = pair[0];
  lts.decaytree[1].p4        = pair[1];
  lts.process.MP_FRAME       = "CM";
  lts.process.FORWARD_NOFLIP = true;
  lts.process.TU_SIGN        = "positive";
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Build one asymmetric pion-pair point near the f2 pole
gra::LORENTZSCALAR AsymmetricF2PhasePointForTest(double theta, double phi) {
  gra::LORENTZSCALAR lts    = ScalarPolePhasePointForTest(0.20, theta, phi, 1.30);
  const auto         proton = [](const double pt, const double azimuth, const double pz) {
    const double energy = std::sqrt(gra::math::pow2(pz) + gra::math::pow2(pt) + gra::math::pow2(gra::PDG::mp));
    return gra::M4Vec(pt * std::cos(azimuth), pt * std::sin(azimuth), pz, energy);
  };
  constexpr double upper_phi = 0.21;
  lts.pfinal[1]              = proton(0.16, upper_phi, 4.30);
  lts.pfinal[2]              = proton(0.27, upper_phi + 2.40, -4.35);
  lts.pfinal[0]              = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  const double central_mass  = lts.pfinal[0].M();
  const auto   pair = TwoBodyRestKinematics(central_mass, lts.decaytree[0].p.mass, lts.decaytree[1].p.mass, theta, phi);
  lts.decaytree[0].p4 = BoostFromRestFrame(pair[0], lts.pfinal[0]);
  lts.decaytree[1].p4 = BoostFromRestFrame(pair[1], lts.pfinal[0]);
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Build an asymmetric two-body central state with selected final particles
gra::LORENTZSCALAR AsymmetricCentralPairLTSForTest(int left_pdg, int right_pdg, double theta, double phi) {
  gra::LORENTZSCALAR       lts       = DirectCentralPairLTSForTest(left_pdg, right_pdg);
  const gra::LORENTZSCALAR reference = AsymmetricF2PhasePointForTest(theta, phi);
  lts.pbeam1                         = reference.pbeam1;
  lts.pbeam2                         = reference.pbeam2;
  lts.pfinal[1]                      = reference.pfinal[1];
  lts.pfinal[2]                      = reference.pfinal[2];
  lts.pfinal[0]                      = lts.pbeam1 + lts.pbeam2 - lts.pfinal[1] - lts.pfinal[2];
  const auto pair =
      TwoBodyRestKinematics(lts.pfinal[0].M(), lts.decaytree[0].p.mass, lts.decaytree[1].p.mass, theta, phi);
  lts.decaytree[0].p4 = BoostFromRestFrame(pair[0], lts.pfinal[0]);
  lts.decaytree[1].p4 = BoostFromRestFrame(pair[1], lts.pfinal[0]);
  RefreshToyDerivedKinematicsPreserveDecay(lts);
  return lts;
}

// Compute the complete reduced scalar Tensor Pomeron t and u graphs
std::array<std::complex<double>, 2> TensorPionContinuumTUForTest(const gra::MTensorPomeron      &tensor,
                                                                 const gra::MTensorPomeronParam &param,
                                                                 const gra::LORENTZSCALAR       &lts) {
  const auto       &meson          = param.FindPseudoscalar(211);
  const gra::M4Vec &pa             = lts.pbeam1;
  const gra::M4Vec &p1             = lts.pfinal[1];
  const gra::M4Vec &p2             = lts.pfinal[2];
  const gra::M4Vec &p3             = lts.decaytree[0].p4;
  const gra::M4Vec &p4             = lts.decaytree[1].p4;
  const gra::M4Vec  pt             = pa - p1 - p3;
  const gra::M4Vec  pu             = p4 - pa + p1;
  const auto        upper_source   = tensor.iG_PForwardHE(gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper));
  const auto        lower_source   = tensor.iG_PForwardHE(gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower));
  const auto        upper_t        = tensor.PomeronPropagatorCurrent(upper_source, (p1 + p3).M2(), lts.t1);
  const auto        lower_t        = tensor.PomeronPropagatorCurrent(lower_source, (p2 + p4).M2(), lts.t2);
  const auto        upper_u        = tensor.PomeronPropagatorCurrent(upper_source, (p1 + p4).M2(), lts.t1);
  const auto        lower_u        = tensor.PomeronPropagatorCurrent(lower_source, (p2 + p3).M2(), lts.t2);
  const auto        vertex_t_upper = tensor.iG_Ppsps(pt, -p3, meson.gPPS);
  const auto        vertex_t_lower = tensor.iG_Ppsps(p4, pt, meson.gPPS);
  const auto        vertex_u_upper = tensor.iG_Ppsps(p4, pu, meson.gPPS);
  const auto        vertex_u_lower = tensor.iG_Ppsps(pu, -p3, meson.gPPS);
  const std::complex<double> scale_t =
      tensor.iD_MES0(pt, lts.decaytree[0].p.mass) *
      gra::math::pow2(gra::regge::FormFactor(pt.M2(), gra::math::pow2(lts.decaytree[0].p.mass), meson.ff_offshell));
  const std::complex<double> scale_u =
      tensor.iD_MES0(pu, lts.decaytree[1].p.mass) *
      gra::math::pow2(gra::regge::FormFactor(pu.M2(), gra::math::pow2(lts.decaytree[1].p.mass), meson.ff_offshell));
  FTensor::Index<'c', 4>     alpha;
  FTensor::Index<'d', 4>     beta;
  const std::complex<double> t = (-gra::math::zi) * (upper_t(alpha, beta) * vertex_t_upper(alpha, beta)) *
                                 (lower_t(alpha, beta) * vertex_t_lower(alpha, beta)) * scale_t;
  const std::complex<double> u = (-gra::math::zi) * (upper_u(alpha, beta) * vertex_u_upper(alpha, beta)) *
                                 (lower_u(alpha, beta) * vertex_u_lower(alpha, beta)) * scale_u;
  return {t, u};
}

// Compute the helicity-summed complex overlap of two amplitude sections
std::complex<double> ComplexSectionOverlapForTest(const std::vector<std::complex<double>> &left,
                                                  const std::vector<std::complex<double>> &right) {
  REQUIRE(left.size() == right.size());
  std::complex<double> out = 0.0;
  for (const auto &i : indices(left)) { out += left[i] * std::conj(right[i]); }
  REQUIRE(std::abs(out) > 0.0);
  return out;
}

// Form one linear combination of two equal helicity sections
std::vector<std::complex<double>> LinearAmplitudeSectionForTest(const std::vector<std::complex<double>> &left,
                                                                const std::vector<std::complex<double>> &right,
                                                                double a, double b) {
  REQUIRE(left.size() == right.size());
  auto out = left;
  for (const auto &i : indices(out)) { out[i] = a * left[i] + b * right[i]; }
  return out;
}

// Compute the factorized direct 4-body phase-space density for the cascade point
double TensorDirectFourBodyPhaseSpaceForTest(const gra::LORENTZSCALAR &lts) {
  const double two_pi = 2.0 * gra::math::PI;
  return (lts.DW.Integral() / two_pi) * (lts.decaytree[0].W_event * lts.decaytree[1].W_event) / gra::math::pow2(two_pi);
}

// Compute the deterministic two-body decay phase-space of one branch
double BranchTwoBodyPhaseSpaceForTest(const gra::MDecayBranch &branch) {
  REQUIRE(branch.legs.size() == 2);
  double M0 = 0.0;
  double m0 = 0.0;
  double m1 = 0.0;
  REQUIRE(branch.p4.PhysicalMass(M0));
  REQUIRE(branch.legs[0].p4.PhysicalMass(m0));
  REQUIRE(branch.legs[1].p4.PhysicalMass(m1));
  return gra::kinematics::dPhi2(M0, gra::kinematics::DecayMomentum(M0, m0, m1));
}

// Compute the recursive internal phase-space factor of one branch
double BranchInternalCascadePhaseSpaceForTest(const gra::MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double product = BranchTwoBodyPhaseSpaceForTest(branch);
  if ((branch.mass_proposal == gra::MassProposal::BreitWigner) || (branch.mass_proposal == gra::MassProposal::Uniform) || (branch.mass_proposal == gra::MassProposal::Fixed)) {
    product /= 2.0 * gra::math::PI;
  }
  for (const auto &leg : branch.legs) { product *= BranchInternalCascadePhaseSpaceForTest(leg); }
  return product;
}

// Compute the internal recursive cascade phase-space factor used by <C>
double InternalCascadePhaseSpaceForTest(const std::vector<gra::MDecayBranch> &tree) {
  double product = 1.0;
  for (const auto &branch : tree) { product *= BranchInternalCascadePhaseSpaceForTest(branch); }
  return product;
}

// Compute the recursive product of BW proposal normalizations for one branch
double BranchMassProposalNormForTest(const gra::MDecayBranch &branch) {
  if (branch.legs.empty()) { return 1.0; }

  double product = (branch.p.width >= 1e-40) ? branch.mass_proposal_norm : 1.0;
  for (const auto &leg : branch.legs) { product *= BranchMassProposalNormForTest(leg); }
  return product;
}

// Compute the recursive product of intermediate BW proposal normalizations
double CascadeMassProposalNormForTest(const std::vector<gra::MDecayBranch> &tree) {
  double product = 1.0;
  for (const auto &branch : tree) { product *= BranchMassProposalNormForTest(branch); }
  return product;
}

// Compute the root two-body phase-space factor used by <F> before cascades
double RootCascadePhaseSpaceForTest(const gra::LORENTZSCALAR &lts, const std::vector<gra::MDecayBranch> &tree) {
  REQUIRE(tree.size() == 2);
  gra::MDecayBranch root;
  root.p4 = lts.pfinal[0];
  root.legs = tree;
  return BranchTwoBodyPhaseSpaceForTest(root);
}

// Compute the phase-space density used by stable-leaf mixture-density ratios
double StableLeafDensityPhaseSpaceForTest(const gra::LORENTZSCALAR &lts, const std::vector<gra::MDecayBranch> &tree) {
  double phase_space = InternalCascadePhaseSpaceForTest(tree);
  if (lts.PS_active && lts.DW.GetN() > 0.0 && tree.size() > 1) { phase_space *= RootCascadePhaseSpaceForTest(lts, tree); }
  return phase_space;
}

// Assign independently normalized BW and flat mass proposals at every cascade depth
void SetMixedMassProposalsForTest(std::vector<gra::MDecayBranch> &tree, unsigned int &mask) {
  for (auto &branch : tree) {
    if (branch.legs.empty()) { continue; }
    branch.mass_proposal = (mask & 1U) != 0 ? gra::MassProposal::Uniform : gra::MassProposal::BreitWigner;
    mask >>= 1U;
    branch.mass_proposal_min2 = 0.25 * branch.p4.M2();
    branch.mass_proposal_max2 = 4.0 * branch.p4.M2();
    branch.mass_proposal_norm = (branch.mass_proposal == gra::MassProposal::Uniform)
        ? branch.mass_proposal_max2 - branch.mass_proposal_min2
        : gra::MRandom::RelativisticBWMass2Integral(branch.p.mass, branch.p.width,
              branch.mass_proposal_min2, branch.mass_proposal_max2);
    REQUIRE(branch.mass_proposal_norm > 0.0);
    SetMixedMassProposalsForTest(branch.legs, mask);
  }
}

// Compute the normalized proposal directly from its real Lorentzian denominator
double MixedMassDensityForTest(const std::vector<gra::MDecayBranch> &tree) {
  double density = 1.0;
  for (const auto &branch : tree) {
    if (branch.legs.empty()) { continue; }
    if ((branch.mass_proposal == gra::MassProposal::BreitWigner) || (branch.mass_proposal == gra::MassProposal::Uniform)) {
      const double s = branch.p4.M2();
      if (s < branch.mass_proposal_min2 || s > branch.mass_proposal_max2) { return 0.0; }
      density /= branch.mass_proposal_norm;
      if ((branch.mass_proposal == gra::MassProposal::BreitWigner)) {
        density /= gra::math::pow2(s - gra::math::pow2(branch.p.mass)) +
                   gra::math::pow2(branch.p.mass * branch.p.width);
      }
    }
    density *= MixedMassDensityForTest(branch.legs);
  }
  return density;
}

// Compute the coherent mass density including the central decay phase-space ratio
double MixedHistoryDensityForTest(gra::LORENTZSCALAR lts) {
  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, lts.decaytree);
  REQUIRE(terms.size() > 1);
  const double phase = StableLeafDensityPhaseSpaceForTest(lts, lts.decaytree);
  double density = 0.0;
  for (const auto &term : terms) {
    density += MixedMassDensityForTest(term.tree) * phase / StableLeafDensityPhaseSpaceForTest(lts, term.tree);
  }
  return density / static_cast<double>(terms.size());
}

// Build one spin-1 continuum cascade branch with deterministic daughters
gra::MDecayBranch ContinuumVectorBranchForTest(const std::string &name, int pdg, const gra::M4Vec &vector_p4,
                                               int d1_pdg, int d2_pdg, double m1, double m2, double theta, double phi,
                                               std::complex<double> coupling, double pole_mass, double width,
                                               double mass_proposal_norm) {
  gra::MDecayBranch branch;
  branch.name                = name;
  branch.p                   = ToyParticle(name, pdg, 2, pole_mass);
  branch.p.P                 = -1;
  branch.p.C                 = -1;
  branch.p.width             = width;
  branch.p4                  = vector_p4;
  branch.mass_proposal_norm  = mass_proposal_norm;
  branch.mass_proposal_min2  = 0.0;
  branch.mass_proposal_max2  = 100.0;
  branch.mass_proposal = gra::MassProposal::BreitWigner;
  branch.hel                 = SpinOneToScalarScalarHelicityMatrix(coupling);

  const auto daughters_in_vector = TwoBodyRestKinematics(vector_p4.M(), m1, m2, theta, phi);

  gra::MDecayBranch d1;
  d1.name = "d1";
  d1.p    = ToyParticle("d1", d1_pdg, 0, m1);
  d1.p.P  = -1;
  d1.p4   = BoostFromRestFrame(daughters_in_vector[0], vector_p4);

  gra::MDecayBranch d2;
  d2.name = "d2";
  d2.p    = ToyParticle("d2", d2_pdg, 0, m2);
  d2.p.P  = -1;
  d2.p4   = BoostFromRestFrame(daughters_in_vector[1], vector_p4);

  branch.legs    = {d1, d2};
  branch.W_event = BranchTwoBodyPhaseSpaceForTest(branch);
  branch.W       = gra::kinematics::MCW(branch.W_event, gra::math::pow2(branch.W_event), 1.0);
  return branch;
}

// Build a continuum two-vector cascade point with repeated stable leaves
gra::LORENTZSCALAR ContinuumVectorCascadeLTSForTest(bool root_phase_space_active) {
  const double mX  = 2.4;
  const double mx  = 0.14;
  const double ma  = 0.20;
  const double mb  = 0.30;
  const double mR1 = 0.82;
  const double mR2 = 0.93;

  gra::LORENTZSCALAR lts;
  lts.PDG             = LoadedPDGTable();
  lts.process.SPINDEC = true;
  lts.process.SPINGEN = true;
  lts.PS_active       = root_phase_space_active;
  lts.pfinal.resize(3);
  lts.pfinal[0] = gra::M4Vec(0.0, 0.0, 0.0, mX);
  lts.pbeam1    = gra::M4Vec(0.0, 0.0, 5.0, std::hypot(5.0, lts.beam1.mass));
  lts.pbeam2    = gra::M4Vec(0.0, 0.0, -5.0, std::hypot(5.0, lts.beam2.mass));

  const auto R_in_X = TwoBodyRestKinematics(mX, mR1, mR2, 0.62, -0.27);

  lts.decaytree = {ContinuumVectorBranchForTest("R1", 900001, R_in_X[0], 111, 211, mx, ma, 0.74, 0.51,
                                                std::complex<double>(0.73, -0.18), 0.775, 0.08, 2.0),
                   ContinuumVectorBranchForTest("R2", 900002, R_in_X[1], 111, 321, mx, mb, 1.38, -0.81,
                                                std::complex<double>(-0.45, 0.32), 0.890, 0.11, 3.0)};

  if (root_phase_space_active) {
    const double root_ps = RootCascadePhaseSpaceForTest(lts, lts.decaytree);
    lts.DW               = gra::kinematics::MCW(root_ps);
  }
  return lts;
}

// Compute a raw continuum cascade matrix including all internal BW factors
MMatrix<std::complex<double>> RawContinuumCascadeMatrixForTest(const gra::LORENTZSCALAR             &lts,
                                                               const std::vector<gra::MDecayBranch> &tree) {
  gra::LORENTZSCALAR term_lts             = lts;
  term_lts.amplitude.DECAY_SYM            = false;
  term_lts.decay_symmetry_proposal_active = false;
  term_lts.decaytree                      = tree;
  term_lts.decay_structure                = {gra::DecayType::Full};
  return gra::spin::ContinuumDecayMatrix(term_lts, "CM");
}

// Compute the coherent raw stable-leaf symmetrized continuum matrix
MMatrix<std::complex<double>> CoherentRawContinuumCascadeMatrixForTest(gra::LORENTZSCALAR                    lts,
                                                                       const std::vector<gra::MDecayBranch> &tree) {
  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, tree);
  REQUIRE_FALSE(terms.empty());

  MMatrix<std::complex<double>> out;
  bool                          initialized = false;
  for (const auto &term : terms) {
    const auto raw = RawContinuumCascadeMatrixForTest(lts, term.tree) * std::complex<double>(term.statistics_sign, 0.0);
    if (!initialized) {
      out         = raw;
      initialized = true;
    } else {
      out = out + raw;
    }
  }
  return out;
}

// Compute the explicit stable-leaf mixture density used by the proposal
double ExplicitStableLeafMixtureDensityForTest(gra::LORENTZSCALAR                    lts,
                                               const std::vector<gra::MDecayBranch> &reference_tree) {
  const auto   terms                 = gra::spin::StableLeafAmplitudeTrees(lts, reference_tree);
  const double reference_phase_space = StableLeafDensityPhaseSpaceForTest(lts, reference_tree);
  REQUIRE(reference_phase_space > 0.0);

  double density = 0.0;
  for (const auto &term : terms) {
    const double norm             = CascadeMassProposalNormForTest(term.tree);
    const double term_phase_space = StableLeafDensityPhaseSpaceForTest(lts, term.tree);
    REQUIRE(norm > 0.0);
    REQUIRE(term_phase_space > 0.0);
    density +=
        gra::math::abs2(gra::spin::CascadeBWProduct(term.tree)) / norm * reference_phase_space / term_phase_space;
  }
  return density / static_cast<double>(terms.size());
}

// Express the mixed decay as an effective vertex for the diagonal Lorentz propagator
FTensor::Tensor1<std::complex<double>, 4> TensorVectorDecayVertexForTest(
    const gra::MTensorPomeron &tp, const gra::LORENTZSCALAR &lts, const gra::MDecayBranch &branch) {
  const auto parameters = gra::GetTensorParam(*lts.model_cache, lts.PDG, lts.process.RESONANCES);
  std::complex<double> coupling = branch.hel.g_decay_TP[0];
  if (parameters->FindVector(branch.p.pdg).width_model == gra::TensorVectorWidthModel::RhoOmega) {
    const auto &mix = parameters->rho_omega;
    const auto delta = mix.Propagator(branch.p4.M2());
    const auto row = mix.Index(branch.p.pdg);
    coupling = (delta[row][0] * mix.g[0] + delta[row][1] * mix.g[1]) / delta[row][row];
  }
  auto vertex = tp.iG_vpsps(branch.legs[0].p4, branch.legs[1].p4, branch.p.mass, 1.0, branch.hel.ff_decay);
  for (const auto &mu : tp.LI) { vertex(mu) *= coupling; }
  return vertex;
}

// Build the covariant raw tensor cascade amplitudes independently of ME6
std::vector<std::complex<double>> TensorRawCascadeAmplitudesForTest(const gra::MTensorPomeron            &tp,
                                                                    const gra::LORENTZSCALAR             &lts,
                                                                    const std::vector<gra::MDecayBranch> &tree) {
  const auto             g = TensorVectorCouplingsForTest(tree[0].p.pdg);
  FTensor::Index<'a', 4> mu1;
  FTensor::Index<'b', 4> nu1;
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'e', 4> rho3;
  FTensor::Index<'f', 4> rho4;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;
  FTensor::Index<'k', 4> mu2;
  FTensor::Index<'l', 4> nu2;
  FTensor::Index<'m', 4> kappa3;
  FTensor::Index<'n', 4> kappa4;

  const gra::M4Vec pa = lts.pbeam1;
  const gra::M4Vec pb = lts.pbeam2;
  const gra::M4Vec p1 = lts.pfinal[1];
  const gra::M4Vec p2 = lts.pfinal[2];
  const gra::M4Vec p3 = tree[0].p4;
  const gra::M4Vec p4 = tree[1].p4;
  const gra::M4Vec pt = pa - p1 - p3;
  const gra::M4Vec pu = p4 - pa + p1;

  const auto iG_1a  = tp.iG_PppHE(p1, pa);
  const auto iG_2b  = tp.iG_PppHE(p2, pb);
  const auto iDP_13 = tp.iD_P((p1 + p3).M2(), lts.t1);
  const auto iDP_24 = tp.iD_P((p2 + p4).M2(), lts.t2);
  const auto iDP_14 = tp.iD_P((p1 + p4).M2(), lts.t1);
  const auto iDP_23 = tp.iD_P((p2 + p3).M2(), lts.t2);

  const double vector_mass = tree[0].p.mass;
  const auto   ff_transfer = TensorVectorTransferForTest(tree[0].p.pdg);
  const auto   iG_tA       = tp.iG_Pvv(pt, -p3, g[0], g[1], ff_transfer);
  const auto   iDV_t       = tp.iD_V(pt, vector_mass, lts.pfinal[0].M2(), tree[0].p.pdg);
  const auto   iG_tB       = tp.iG_Pvv(p4, pt, g[0], g[1], ff_transfer);
  const auto   iG_uA       = tp.iG_Pvv(p4, pu, g[0], g[1], ff_transfer);
  const auto   iDV_u       = tp.iD_V(pu, vector_mass, lts.pfinal[0].M2(), tree[0].p.pdg);
  const auto   iG_uB       = tp.iG_Pvv(pu, -p3, g[0], g[1], ff_transfer);

  const bool index_up          = true;
  const bool conserved_current = true;
  const auto iD_V1 =
      tp.iD_VMES(tree[0].p4, tree[0].p.mass, tree[0].p.width, tree[0].p.pdg, index_up, conserved_current);
  const auto iD_V2 =
      tp.iD_VMES(tree[1].p4, tree[1].p.mass, tree[1].p.width, tree[1].p.pdg, index_up, conserved_current);
  const auto iG_D1 = TensorVectorDecayVertexForTest(tp, lts, tree[0]);
  const auto iG_D2 = TensorVectorDecayVertexForTest(tp, lts, tree[1]);

  FTensor::Tensor2<std::complex<double>, 4, 4> decay;
  decay(rho3, rho4) = iD_V1(rho3, kappa3) * iG_D1(kappa3) * iD_V2(rho4, kappa4) * iG_D2(kappa4);

  FTensor::Tensor2<std::complex<double>, 4, 4> M_t;
  {
    FTensor::Tensor2<std::complex<double>, 4, 4> A;
    A(rho1, rho3) = iG_1a(mu1, nu1) * iDP_13(mu1, nu1, alpha1, beta1) * iG_tA(rho1, rho3, alpha1, beta1);
    FTensor::Tensor2<std::complex<double>, 4, 4> B;
    B(rho4, rho2)   = iG_2b(mu2, nu2) * iDP_24(alpha2, beta2, mu2, nu2) * iG_tB(rho4, rho2, alpha2, beta2);
    M_t(rho3, rho4) = A(rho1, rho3) * iDV_t(rho1, rho2) * B(rho4, rho2) *
                      gra::math::pow2(TensorVectorOffShellForTest(tree[0].p.pdg, pt.M2(), vector_mass));
  }

  FTensor::Tensor2<std::complex<double>, 4, 4> M_u;
  {
    FTensor::Tensor2<std::complex<double>, 4, 4> A;
    A(rho4, rho1) = iG_1a(mu1, nu1) * iDP_14(mu1, nu1, alpha1, beta1) * iG_uA(rho4, rho1, alpha1, beta1);
    FTensor::Tensor2<std::complex<double>, 4, 4> B;
    B(rho2, rho3)   = iG_2b(mu2, nu2) * iDP_23(alpha2, beta2, mu2, nu2) * iG_uB(rho2, rho3, alpha2, beta2);
    M_u(rho3, rho4) = A(rho4, rho1) * iDV_u(rho1, rho2) * B(rho2, rho3) *
                      gra::math::pow2(TensorVectorOffShellForTest(tree[0].p.pdg, pu.M2(), vector_mass));
  }

  const std::complex<double> amp = (-gra::math::zi) * (M_t(rho3, rho4) + M_u(rho3, rho4)) * decay(rho3, rho4);

  std::vector<std::complex<double>> out;
  for (int ha = 0; ha < 2; ++ha) {
    for (int hb = 0; hb < 2; ++hb) {
      for (int h1 = 0; h1 < 2; ++h1) {
        for (int h2 = 0; h2 < 2; ++h2) {
          if (ha == h1 && hb == h2) { out.push_back(amp); }
        }
      }
    }
  }
  return out;
}

// Build the raw scalar-resonance vector-cascade amplitudes independently of ME3
std::vector<std::complex<double>> ScalarRawCascadeAmplitudesForTest(const gra::MTensorPomeron &tp,
                                                                    const gra::LORENTZSCALAR  &lts,
                                                                    const gra::PARAM_RES      &res) {
  FTensor::Index<'a', 4> mu1;
  FTensor::Index<'b', 4> nu1;
  FTensor::Index<'c', 4> rho1;
  FTensor::Index<'d', 4> rho2;
  FTensor::Index<'g', 4> alpha1;
  FTensor::Index<'h', 4> beta1;
  FTensor::Index<'i', 4> alpha2;
  FTensor::Index<'j', 4> beta2;
  FTensor::Index<'k', 4> mu2;
  FTensor::Index<'l', 4> nu2;
  FTensor::Index<'m', 4> kappa1;
  FTensor::Index<'n', 4> kappa2;

  const auto iG_1a = tp.iG_PppHE(lts.pfinal[1], lts.pbeam1);
  const auto iG_2b = tp.iG_PppHE(lts.pfinal[2], lts.pbeam2);
  const double nu1x2 = (lts.pbeam1 + lts.pfinal[1]) * lts.pfinal[0];
  const double nu2x2 = (lts.pbeam2 + lts.pfinal[2]) * lts.pfinal[0];
  const auto iDP_1 = tp.iD_P(nu1x2, lts.t1);
  const auto iDP_2 = tp.iD_P(nu2x2, lts.t2);
  const auto cvtx =
      tp.iG_PPS_total(lts.q1, lts.q2, res.p.mass, gra::TensorResonanceType::Scalar, ToyTensorChannel(res).g_tensor, ToyTensorChannel(res));
  const std::complex<double> iD = tp.iD_MES(lts.pfinal[0], res.p.mass, res.p.width);

  const auto iGf0vv = tp.iG_f0vv(lts.decaytree[0].p4, lts.decaytree[1].p4, res.p.mass, res.hel_decay.g_decay_TP[0],
                                 res.hel_decay.g_decay_TP[1], res.hel_decay.ff_decay);
  const bool index_up          = true;
  const bool conserved_current = true;
  const auto iD_V1             = tp.iD_VMES(lts.decaytree[0].p4, lts.decaytree[0].p.mass, lts.decaytree[0].p.width,
                                            lts.decaytree[0].p.pdg, index_up, conserved_current);
  const auto iD_V2             = tp.iD_VMES(lts.decaytree[1].p4, lts.decaytree[1].p.mass, lts.decaytree[1].p.width,
                                            lts.decaytree[1].p.pdg, index_up, conserved_current);
  const auto iG_D1 = TensorVectorDecayVertexForTest(tp, lts, lts.decaytree[0]);
  const auto iG_D2 = TensorVectorDecayVertexForTest(tp, lts, lts.decaytree[1]);

  FTensor::Tensor2<std::complex<double>, 4, 4> decay;
  decay(rho1, rho2) = iD_V1(rho1, kappa1) * iG_D1(kappa1) * iD_V2(rho2, kappa2) * iG_D2(kappa2);

  const std::complex<double> decay_scalar = iGf0vv(rho1, rho2) * decay(rho1, rho2);
  const std::complex<double> amp          = (-gra::math::zi) * iG_1a(mu1, nu1) * iDP_1(mu1, nu1, alpha1, beta1) *
                                   cvtx(alpha1, beta1, alpha2, beta2) * iD * iDP_2(alpha2, beta2, mu2, nu2) *
                                   iG_2b(mu2, nu2) * decay_scalar;

  return {amp, amp, amp, amp};
}

// Compute the coherent raw stable-leaf symmetrized tensor amplitudes
std::vector<std::complex<double>> CoherentRawTensorCascadeAmplitudesForTest(
    gra::LORENTZSCALAR lts, const gra::MTensorPomeron &tp, const std::vector<gra::MDecayBranch> &tree) {
  const auto terms = gra::spin::StableLeafAmplitudeTrees(lts, tree);
  REQUIRE_FALSE(terms.empty());

  std::vector<std::complex<double>> out;
  for (const auto &term : terms) {
    const auto raw = TensorRawCascadeAmplitudesForTest(tp, lts, term.tree);
    if (out.empty()) { out.assign(raw.size(), 0.0); }
    REQUIRE(out.size() == raw.size());
    for (std::size_t i = 0; i < raw.size(); ++i) { out[i] += term.statistics_sign * raw[i]; }
  }
  return out;
}

// Compute helicity averaged amplitude squared for a test helicity vector
double HelicityAverageAmp2ForTest(const std::vector<std::complex<double>> &hamp) {
  double out = 0.0;
  for (const auto &amp : hamp) { out += gra::math::abs2(amp); }
  return out / 4.0;
}

// Compute scalar screened toy color-flow amplitudes

}  // namespace

// Build one direct pion pair accepted by the toy amplitude process routes
std::vector<gra::MDecayBranch> ToyPionPairDecayTree() {
  gra::MDecayBranch positive;
  positive.p = ToyParticle("pi+", gra::PDG::PDG_pip, 0, 0.13957061);
  gra::MDecayBranch negative;
  negative.p = ToyParticle("pi-", gra::PDG::PDG_pim, 0, 0.13957061);
  return {positive, negative};
}

// Populate the minimal event-local data required by one real process route
void ConfigureToyAmplitude(gra::MSubProc &subprocess, gra::LORENTZSCALAR &lts, const std::string &family,
                           const std::string &channel) {
  subprocess.ISTATE  = family;
  subprocess.CHANNEL = channel;
  lts.process.RESONANCES.clear();
  lts.process.CONT_PRODUCTION.clear();
  lts.process.CONT_PRODUCTIONTREE.clear();
  lts.process.CONTINUUM_POLE.clear();
  lts.process.CONTINUUM_GP.clear();
  lts.process.TENSOR_MODEL_READY = false;
  lts.process.TENSOR_PSEUDOSCALAR_PDGS.clear();
  lts.process.TENSOR_BARYON_PDGS.clear();
  lts.process.TENSOR_VECTOR_DECAY_PDGS.clear();
  lts.process.SPINGEN = true;

  const bool regge     = family == "MP" || family == "XP" || family == "GP";
  const bool resonance = channel == "RES" || channel == "RES+CON";
  const bool continuum = channel == "CON" || channel == "RES+CON";
  if (regge && continuum) {
    lts.process.CONT_PRODUCTION     = {{991, 991}};
    lts.process.CONT_PRODUCTIONTREE = {std::vector<gra::MDecayBranch>(2)};
    lts.process.CONTINUUM_POLE      = MakeToyScalarContinuumLTSAsymmetric(0.3, 4.2, -4.6).process.CONTINUUM_POLE;
  }
  if (regge && resonance) { lts.process.RESONANCES = {{"toy_regge", gra::PARAM_RES{}}}; }
  if (family != "TP") { return; }

  if (resonance) {
    gra::PARAM_RES scalar;
    scalar.p                        = ToyParticle("f0", 9000021, 0, 1.0);
    scalar.p.C                      = 1;
    scalar.p.width                  = 0.1;
    scalar.hel_decay.g_decay_TP = {0.7, -0.2};
    lts.process.RESONANCES          = {{"toy_tensor", scalar}};
  }
  const bool vector_cascade =
      std::any_of(lts.decaytree.begin(), lts.decaytree.end(), [](const auto &branch) { return !branch.legs.empty(); });
  if (!continuum && !vector_cascade) { return; }

  lts.process.TENSOR_MODEL_READY = true;
  for (const auto &branch : lts.decaytree) {
    const int pdg = std::abs(branch.p.pdg);
    if (!branch.legs.empty()) {
      lts.process.TENSOR_VECTOR_DECAY_PDGS[pdg] = std::abs(branch.legs.front().p.pdg);
    } else if (branch.p.spinX2 == 0) {
      lts.process.TENSOR_PSEUDOSCALAR_PDGS.push_back(pdg);
    } else if (branch.p.spinX2 == 1) {
      lts.process.TENSOR_BARYON_PDGS.push_back(pdg);
    }
  }
}

class ToyHelicityProcess : public gra::MFactorized {
 public:
  // Exercise the real screening shift and Born restoration paths
  using gra::MFactorized::LoopKinematics;
  using gra::MProcess::RestoreBornKinematics;

  // Exercise the production tree construction and mixture selection directly
  using gra::MProcess::ApplyDecaySymmetryProposal;
  using gra::MProcess::ConstructDecayKinematics;
  using gra::MProcess::CreateProcessSetup;
  using gra::MProcess::DecaySymmetryCompensationFactor;
  using gra::MProcess::PrepareDecaySymmetryProposal;

  // Initialize the complete immutable TUNE0 model and particle snapshots
  ToyHelicityProcess() {
    gra::MODELPARAM = "TUNE0";
    state.lts.PDG   = LoadedPDGTable();
    SetModelTune(gra::MModelTune::Load(modelfile));
    SetHelicityConfig(GetModelTune());
    ProcPtr.ISTATE  = "MP";
    ProcPtr.CHANNEL = "RES";
  }

  // Load one temporary model tune for both process and helicity readers
  void SetTuneForTest(const std::string &path) {
    const auto tune = gra::MModelTune::Load((std::filesystem::path(path) / "GENERAL.json").string());
    SetModelTune(tune);
    SetHelicityConfig(tune);
  }

  // Set the subprocess family and channel used by setup-level reader tests
  void SetProcessForTest(const std::string &family, const std::string &channel) {
    ProcPtr.ISTATE  = family;
    ProcPtr.CHANNEL = channel;
  }

  // Access the configured signed 64-bit custom-cut identifier
  std::int64_t UserCutsForTest() const { return state.usercuts; }

  // Access the configured veto-cut domains
  const gra::VETOCUT &VetoCutsForTest() const { return state.vetocuts; }

  // Access the configured quasielastic sampling minimum
  double QuasiElasticAbsTMinForTest() const { return state.gcuts.q_t_abs_min; }

  // Access the configured quasielastic sampling maximum
  double QuasiElasticAbsTMaxForTest() const { return state.gcuts.q_t_abs_max; }

  // Configure phase-space metadata for fiducial validation tests
  void ConfigureFiducialValidationForTest(const std::string &phase_space, const std::string &channel, int excitation) {
    state.phase_space_class = phase_space;
    ProcPtr.CHANNEL         = channel;
    state.excitation        = excitation;
  }

  // Compute the process-level final-state symmetry factor for a test tree
  double SymmetryFactorForTest(const std::vector<gra::MDecayBranch> &tree, bool decay_sym, bool isolated = false) {
    SetISOLATE(isolated);
    state.lts.decaytree           = tree;
    state.lts.amplitude.DECAY_SYM = decay_sym;
    CalculateSymmetryFactor();
    return state.symmetry_factor;
  }

  // Compute the unsymmetrized diagnostic assignment compensation factor
  double DecaySymmetryCompensationFactorForTest(const std::vector<gra::MDecayBranch> &tree, bool decay_sym,
                                                const std::string &family = "MP", const std::string &channel = "RES",
                                                bool isolated = false) {
    SetISOLATE(isolated);
    ProcPtr.ISTATE                = family;
    ProcPtr.CHANNEL               = channel;
    state.lts.decaytree           = tree;
    state.lts.amplitude.DECAY_SYM = decay_sym;
    ConfigureToyAmplitude(ProcPtr, state.lts, family, channel);
    return DecaySymmetryCompensationFactor();
  }

  // Compute whether stable-leaf proposal sampling activates for a test tree
  bool StableLeafProposalActiveForTest(const std::vector<gra::MDecayBranch> &tree, bool decay_sym,
                                       const std::string &family, const std::string &channel, bool spindec = true) {
    ProcPtr.ISTATE                = family;
    ProcPtr.CHANNEL               = channel;
    state.lts.decaytree           = tree;
    state.lts.amplitude.DECAY_SYM = decay_sym;
    state.lts.process.SPINDEC     = spindec;
    ConfigureToyAmplitude(ProcPtr, state.lts, family, channel);
    PrepareDecaySymmetryProposal();
    return state.lts.decay_symmetry_proposal_active;
  }

  // Apply one prepared stable-leaf proposal and return its generated phase
  // space
  double AppliedStableLeafProposalPhaseSpaceForTest(const std::vector<gra::MDecayBranch> &tree,
                                                    gra::CentralPhaseSpaceMode mode, bool root_phase_space_active) {
    state.lts.decaytree                      = tree;
    state.lts.central_phase_space_mode       = mode;
    state.lts.PS_active                      = root_phase_space_active;
    state.lts.DW                             = gra::kinematics::MCW(0.37, gra::math::pow2(0.37), 1.0);
    state.lts.decay_symmetry_proposal_active = true;
    state.lts.decay_symmetry_proposal_index  = 0;
    state.lts.decay_symmetry_assignments     = {{0, 1, 2, 3}};
    state.lts.pfinal.assign(tree.size() + 3, gra::M4Vec(0.0, 0.0, 0.0, 0.0));
    if (!ApplyDecaySymmetryProposal()) { throw std::runtime_error("Failed to apply toy stable-leaf proposal"); }
    return state.lts.decay_symmetry_proposal_phase_space;
  }

  // Compute the cascade factor after common amplitude decay compensation
  double CascadePhaseSpaceForTest(const gra::LORENTZSCALAR &event, gra::DecayStructure decay_structure) {
    state.lts                 = event;
    state.lts.decay_structure = decay_structure;
    return CascadePS();
  }

  // Compute the common fiducial-cut decision for a prepared central tree
  bool CommonCutsForTest(const gra::FIDCUT &cuts, const std::vector<gra::MDecayBranch> &tree, double system_mass = 10.0,
                         double system_rapidity = 0.0, double system_pt = 0.0) {
    state.fcuts             = cuts;
    state.lts.decaytree     = tree;
    state.lts.m2            = system_mass * system_mass;
    state.lts.Y             = system_rapidity;
    state.lts.Pt            = system_pt;
    state.phase_space_class = "P";
    state.usercuts          = 0;
    return CommonCuts();
  }

  // Compute the common fiducial-cut decision for prepared forward-system
  // kinematics
  bool CommonForwardCutsForTest(const gra::FIDCUT &cuts, const gra::LORENTZSCALAR &event) {
    state.fcuts             = cuts;
    state.lts               = event;
    state.phase_space_class = "F";
    state.usercuts          = 0;
    return CommonCuts();
  }

};

// Testable quasielastic process exposing prepared fiducial kinematics
class ToyQuasiElasticProcess : public gra::MQuasiElastic {
 public:
  // Compute the quasielastic fiducial decision for one prepared event
  bool FiducialCutsForTest(const gra::FIDCUT &cuts, const gra::LORENTZSCALAR &event, const std::string &channel) {
    state.fcuts     = cuts;
    state.lts       = event;
    ProcPtr.CHANNEL = channel;
    state.usercuts  = 0;
    return FiducialCuts();
  }
};

// Configure a setup-level toy process with common pp beams and one decay mode
//
void ConfigureToyProductionProcess(ToyHelicityProcess &proc, const std::string &family, const std::string &channel,
                                   const std::string &decay_mode) {
  gra::MODELPARAM = "TUNE0";
  proc.SetProcessForTest(family, channel);
  proc.state.lts.PDG   = LoadedPDGTable();
  proc.state.lts.beam1 = proc.state.lts.PDG.FindByPDG(2212);
  proc.state.lts.beam2 = proc.state.lts.PDG.FindByPDG(2212);
  proc.SetDecayMode(decay_mode);
}

// Build one stable leaf with a toy particle and fixed four-momentum
gra::MDecayBranch FiducialLeafForTest(int pdg, const gra::M4Vec &p4) {
  gra::MDecayBranch branch;
  branch.p  = ToyParticle(pdg, 0, 1, 0, "fid_leaf");
  branch.p4 = p4;
  return branch;
}

// Build one massless stable leaf from transverse momentum and pseudorapidity
gra::MDecayBranch MasslessFiducialLeafForTest(int pdg, double pt, double eta) {
  gra::M4Vec p4;
  p4.SetPxPyPzM(pt, 0.0, pt * std::sinh(eta), 0.0);
  return FiducialLeafForTest(pdg, p4);
}

// Compute a broad default fiducial-cut configuration for focused PDG tests
gra::FIDCUT BroadFiducialCutsForTest() {
  gra::FIDCUT cuts;
  cuts.active              = true;
  cuts.particle_eta_active = true;
  cuts.eta_min             = -100.0;
  cuts.eta_max             = 100.0;
  cuts.particle_rap_active = true;
  cuts.rap_min             = -100.0;
  cuts.rap_max             = 100.0;
  cuts.particle_pt_active  = true;
  cuts.pt_min              = 0.0;
  cuts.pt_max              = 1000000.0;
  cuts.particle_Et_active  = true;
  cuts.Et_min              = 0.0;
  cuts.Et_max              = 1000000.0;
  cuts.system_M_active     = true;
  cuts.M_min               = 0.0;
  cuts.M_max               = 1000000.0;
  cuts.system_Rap_active   = true;
  cuts.Y_min               = -100.0;
  cuts.Y_max               = 100.0;
  cuts.system_Pt_active    = true;
  cuts.Pt_min              = 0.0;
  cuts.Pt_max              = 1000000.0;
  return cuts;
}

// Compute one active PDG fiducial range
gra::FIDCUTRANGE FiducialRangeForTest(double min, double max) {
  gra::FIDCUTRANGE range;
  range.active = true;
  range.min    = min;
  range.max    = max;
  return range;
}

// Compute one explicit PDG selector test cut
gra::FIDPDGCUT PDGFiducialCutForTest(std::vector<int> pdg, std::vector<bool> pdg_abs) {
  if (pdg.size() != pdg_abs.size()) {
    throw std::invalid_argument("PDGFiducialCutForTest: selector dimension mismatch");
  }
  gra::FIDPDGCUT cut;
  cut.pdg     = std::move(pdg);
  cut.pdg_abs = std::move(pdg_abs);
  return cut;
}

// Parse one FIDCUTS object through the real setup reader
gra::FIDCUT ParsedFiducialCutsForTest(const nlohmann::json &card) {
  gra::MGraniitti    generator;
  ToyHelicityProcess proc;
  generator.proc = &proc;
  generator.ReadFidCuts(card);
  return proc.state.fcuts;
}

// Parse one quasielastic GENCUTS object through the real setup reader
std::array<double, 2> ParsedQuasiElasticAbsTRangeForTest(const nlohmann::json &card, const std::string &family = "X",
                                                         const std::string &channel = "EL", bool screening = true) {
  auto generator     = std::make_unique<gra::MGraniitti>();
  auto proc          = std::make_unique<ToyHelicityProcess>();
  generator->PROCESS = "X[EL]<Q>";
  generator->proc    = proc.get();
  proc->SetProcessForTest(family, channel);
  proc->SetScreening(screening);
  generator->ReadGenCuts(card);
  return {proc->QuasiElasticAbsTMinForTest(), proc->QuasiElasticAbsTMaxForTest()};
}

// Parse one VETOCUTS object through the real setup reader
gra::VETOCUT ParsedVetoCutsForTest(const nlohmann::json &card) {
  gra::MGraniitti    generator;
  ToyHelicityProcess proc;
  generator.proc = &proc;
  generator.ReadVetoCuts(card);
  return proc.VetoCutsForTest();
}

// Parse and return one USERCUTS value through the real setup reader
std::int64_t ParsedUserCutIDForTest(const nlohmann::json &card) {
  gra::MGraniitti    generator;
  ToyHelicityProcess proc;
  generator.proc = &proc;
  generator.ReadFidCuts(card);
  return proc.UserCutsForTest();
}

// Require one finite reduced helicity tensor without imposing a norm
//
void RequireReducedHelicityNormalization(const gra::HELMatrix &hel) {
  REQUIRE(hel.T.IsFinite());
  REQUIRE(hel.T.FrobNorm2() > 0.0);
}

// Require one initialized canonical pole vertex without imposing a norm
void RequireCanonicalPoleVertex(const gra::spin::PoleLS &vertex) {
  REQUIRE(vertex.ready);
  REQUIRE_FALSE(vertex.terms.empty());
  REQUIRE(vertex.terms.size() == vertex.raw_normalization.size());
  REQUIRE(vertex.terms.size() == vertex.reduced_basis.size());
  REQUIRE(vertex.rows > 0);
  REQUIRE(vertex.cols > 0);
  REQUIRE(vertex.helicity.jw_rotation.ready);
  REQUIRE(vertex.helicity.jw_rotation.row.size() == vertex.helicity.lambda_values.size_row());
  REQUIRE(vertex.dual_spin.size() == vertex.helicity.Jz_values.size());
  REQUIRE(vertex.spin_metric.size() == vertex.helicity.Jz_values.size());
  for (const auto &basis : vertex.reduced_basis) {
    REQUIRE(basis.size_row() == vertex.rows);
    REQUIRE(basis.size_col() == vertex.cols);
  }
  const auto reduced = gra::spin::PoleLSReduced(vertex, vertex.Lambda);
  REQUIRE(reduced.IsFinite());
  REQUIRE(reduced.FrobNorm2() > 0.0);
}

// Compute the unscaled leading pole density of one fixed-spin resonance
double ConfiguredPoleReferenceDensity(const gra::PARAM_RES &res, const gra::ReggeProductionModel model,
                                      const std::size_t channel) {
  REQUIRE(channel < res.production.size());
  REQUIRE(res.production[channel].tree.size() == 2);

  if (model != gra::ReggeProductionModel::MP && model != gra::ReggeProductionModel::XP) {
    throw std::invalid_argument("ConfiguredPoleReferenceDensity: model must be MP or XP");
  }

  REQUIRE(channel < res.production.size());
  return gra::spin::LeadingPoleDensity(res.production[channel].pole.value(),
                                       res.production[channel].pole.value().Lambda);
}

// Compute the number of active alpha_ls rows
//
std::size_t ActiveAlphaRowCount(const gra::HELMatrix &hel) { return hel.alpha_ls.Size(); }

// Require that one normalized LS row is the only active row
//
void RequireSingleAlphaRow(const gra::HELMatrix &hel, std::size_t l, std::size_t two_s) {
  REQUIRE(ActiveAlphaRowCount(hel) == 1);
  REQUIRE(hel.alpha_ls.Contains(l, two_s));
  REQUIRE(std::abs(hel.alpha_ls.At(l, two_s)) == Approx(1.0).margin(1e-12));
}

// Require that two production LS rows have equal relative coefficients
//
void RequireEqualAlphaRows(const gra::HELMatrix &hel, std::size_t l1, std::size_t two_s1, std::size_t l2,
                           std::size_t two_s2) {
  REQUIRE(hel.alpha_ls.Contains(l1, two_s1));
  REQUIRE(hel.alpha_ls.Contains(l2, two_s2));
  REQUIRE(std::abs(hel.alpha_ls.At(l1, two_s1)) == Approx(1.0).margin(1e-12));
  REQUIRE(std::abs(hel.alpha_ls.At(l2, two_s2)) == Approx(1.0).margin(1e-12));
}

// Compute one normalized single-LS helicity matrix definition for direct
// validation
//
gra::HELMatrix SingleLSDefinition(std::size_t l, std::size_t two_s) {
  gra::HELMatrix hel;
  hel.BR         = 1.0;
  hel.P_symmetry = true;
  hel.C_symmetry = true;
  hel.alpha_ls.Set(l, two_s, 1.0);
  return hel;
}

// Convert a computed LS matrix into a direct-helicity input matrix
//
gra::HELMatrix DirectDefinitionFromComputedT(const gra::HELMatrix &computed) {
  gra::HELMatrix direct;
  direct.BR             = 1.0;
  direct.P_symmetry     = true;
  direct.C_symmetry     = true;
  direct.coupling_basis = gra::CouplingBasis::Helicity;
  direct.T              = computed.T;
  direct.T_set          = MMatrix<bool>(computed.T.size_row(), computed.T.size_col(), false);
  for (std::size_t i = 0; i < computed.T.size_row(); ++i) {
    for (std::size_t j = 0; j < computed.T.size_col(); ++j) {
      if (std::abs(computed.T[i][j]) > 1e-12) { direct.T_set[i][j] = true; }
    }
  }
  return direct;
}

// Compute a direct rho <- gamma + Pomeron helicity validation matrix
//
gra::HELMatrix DirectRhoGammaPomeronMatrix() {
  gra::HELMatrix direct;
  direct.BR             = 1.0;
  direct.P_symmetry     = true;
  direct.C_symmetry     = true;
  direct.coupling_basis = gra::CouplingBasis::Helicity;
  direct.T              = MMatrix<std::complex<double>>(3, 1, 0.0);
  direct.T_set          = MMatrix<bool>(3, 1, false);

  const double schc  = 1.0 / std::sqrt(2.0);
  direct.T[0][0]     = schc;
  direct.T[2][0]     = schc;
  direct.T_set[0][0] = true;
  direct.T_set[2][0] = true;
  return direct;
}

#if defined(__GNUC__)
#pragma GCC diagnostic pop
#endif
