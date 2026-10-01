// Nuclear EPA beam and phase-space integration tests
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdint>
#include <filesystem>
#include <fstream>
#include <limits>
#include <memory>
#include <numeric>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Kinematics/MEventKinematics.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/Math/MIntegration.h"
#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Nuclear/MFinal.h"
#include "Graniitti/Nuclear/MEMD.h"
#include "Graniitti/Nuclear/MIncoherent.h"
#include "Graniitti/Nuclear/MSteering.h"
#include "Graniitti/Nuclear/MUPC.h"
#include "Graniitti/Nuclear/MUPCScreen.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Regge/MReggeGP.h"
#include "Graniitti/Tech/MAux.h"
#include "Graniitti/Tech/MException.h"
#include "Graniitti/Tensor/MTensorPhoto.h"
#include "support/nuclear_test_support.hh"

// Libraries
#include <HepMC3/Attribute.h>
#include <HepMC3/GenEvent.h>
#include <HepMC3/GenHeavyIon.h>
#include <HepMC3/GenVertex.h>

#include <catch.hpp>
#include <json.hpp>

namespace {

using gra::aux::indices;

// Test a radiative ion channel through the production mass and event APIs
class RadiativeIon final : public gra::nuclear::MReaction {
 public:
  // Create one independent worker without a sampled state
  std::unique_ptr<gra::nuclear::MReaction> Clone() const override { return std::make_unique<RadiativeIon>(); }

  // Identify the analytic integration fixture
  HepMC3::GenRunInfo::ToolInfo Tool() const override {
    return {"radiative ion", "test", "Production mass integration fixture"};
  }

  // Propose a discrete excitation above the isotope ground state
  gra::nuclear::RecoilMass SampleMasses(const HepMC3::GenEvent &event, gra::MRandom &) override {
    std::array<double, 2> mass{};
    for (const auto &i : indices(mass)) {
      ground_[i] = event.beams()[i]->generated_mass();
      mass[i]    = ground_[i] + 0.02;
    }
    return {mass, {0.02, 0.02}};
  }

  // Attach a normalized isotropic photon decay without changing the hard recoil
  double Complete(HepMC3::GenEvent &event, const std::array<HepMC3::GenParticlePtr, 2> &parents,
                  gra::MRandom &random) override {
    for (const auto &i : indices(parents)) {
      const auto &parent = parents[i];
      const auto  beam   = parent->production_vertex()->particles_in().front();
      parent->set_pid(beam->pid());
      parent->set_status(2);
      const auto             &p = parent->momentum();
      std::vector<gra::M4Vec> daughters;
      const auto              phase = gra::kinematics::TwoBodyPhaseSpace(
                       gra::M4Vec(p.px(), p.py(), p.pz(), p.e()), parent->generated_mass(), {0.0, ground_[i]}, daughters, random);
      REQUIRE(phase.GetW() > 0.0);
      auto vertex = std::make_shared<HepMC3::GenVertex>();
      vertex->add_particle_in(parent);
      auto gamma = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(daughters[0]), 22, 1);
      auto ion   = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(daughters[1]), beam->pid(), 1);
      gamma->set_generated_mass(0.0);
      ion->set_generated_mass(ground_[i]);
      vertex->add_particle_out(gamma);
      vertex->add_particle_out(ion);
      event.add_vertex(vertex);
    }
    return 1.0;
  }

 private:
  std::array<double, 2> ground_{};
};

// Store one physical beam pair used by the supported-process matrix
struct BeamCase {
  std::string                           name;
  std::array<std::string, 2>            beam;
  std::array<double, 2>                 energy;
  std::array<gra::nuclear::BeamType, 2> type;
  bool                                  nuclear = false;
};

// Compute the supported pp, pA, AA, ep, eA and ee beam matrix
std::vector<BeamCase> SupportedBeams() {
  using gra::nuclear::BeamType;
  return {
      {"pp", {"p+", "p+"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"pA", {"p+", "Pb208"}, {6500.0, 2510.0}, {BeamType::Proton, BeamType::Nucleus}, true},
      {"AA", {"Pb208", "Pb208"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
      {"ep", {"e-", "p+"}, {60.0, 6500.0}, {BeamType::Lepton, BeamType::Proton}, false},
      {"eA", {"e-", "Pb208"}, {60.0, 2510.0}, {BeamType::Lepton, BeamType::Nucleus}, true},
      {"ee", {"e-", "e+"}, {100.0, 100.0}, {BeamType::Lepton, BeamType::Lepton}, false},
  };
}

// Compute every hard gamma-gamma amplitude using noncollinear EPA sources
std::vector<std::string> HardPhotonProcesses() {
  std::vector<std::string> processes = {
      "yy[EPA]<F> -> mu+ mu-",
      "yy[Higgs]<F> &> b b~",
      "yy[monopolium(0)]<F> &> gamma gamma "
      "@PDG[881]{M:527.57673,W:10}",
      "yy[FLUX]<F>",
      "yy[jj]<F> -> u u~",
      "yy[WW]<F> -> W+ > {mu+ vm} W- > {e- ve~}",
      "yy[Zjj]<F> -> Z > {mu+ mu-} u u~",
      "yy[yy]<F> -> gamma gamma",
  };
  return processes;
}

// Compute reversed and charge-conjugated representatives of every beam class
std::vector<BeamCase> OrderedChargeBeams() {
  using gra::nuclear::BeamType;
  return {
      {"p pbar", {"p+", "p-"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"pbar p", {"p-", "p+"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"e- e+", {"e-", "e+"}, {100.0, 100.0}, {BeamType::Lepton, BeamType::Lepton}, false},
      {"e+ e-", {"e+", "e-"}, {100.0, 100.0}, {BeamType::Lepton, BeamType::Lepton}, false},
      {"e+ pbar", {"e+", "p-"}, {60.0, 6500.0}, {BeamType::Lepton, BeamType::Proton}, false},
      {"pbar e+", {"p-", "e+"}, {6500.0, 60.0}, {BeamType::Proton, BeamType::Lepton}, false},
      {"pbar anti-A", {"p-", "-1000822080"}, {6500.0, 2510.0}, {BeamType::Proton, BeamType::Nucleus}, true},
      {"anti-A pbar", {"-1000822080", "p-"}, {2510.0, 6500.0}, {BeamType::Nucleus, BeamType::Proton}, true},
      {"e+ anti-A", {"e+", "-1000822080"}, {60.0, 2510.0}, {BeamType::Lepton, BeamType::Nucleus}, true},
      {"anti-A e+", {"-1000822080", "e+"}, {2510.0, 60.0}, {BeamType::Nucleus, BeamType::Lepton}, true},
      {"A anti-A", {"Pb208", "-1000822080"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
      {"anti-A A", {"-1000822080", "Pb208"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
  };
}

// Compute a bounded ordered and charge-conjugated photoproduction event matrix
std::vector<BeamCase> PhotoEventBeams() {
  using gra::nuclear::BeamType;
  return {
      {"pp", {"p+", "p+"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"ep", {"e-", "p+"}, {60.0, 6500.0}, {BeamType::Lepton, BeamType::Proton}, false},
      {"pbar e+", {"p-", "e+"}, {6500.0, 60.0}, {BeamType::Proton, BeamType::Lepton}, false},
      {"pA", {"p+", "Pb208"}, {6500.0, 2510.0}, {BeamType::Proton, BeamType::Nucleus}, true},
      {"anti-A pbar", {"-1000822080", "p-"}, {2510.0, 6500.0}, {BeamType::Nucleus, BeamType::Proton}, true},
      {"eA", {"e-", "Pb208"}, {60.0, 2510.0}, {BeamType::Lepton, BeamType::Nucleus}, true},
      {"anti-A e+", {"-1000822080", "e+"}, {2510.0, 60.0}, {BeamType::Nucleus, BeamType::Lepton}, true},
      {"AA", {"Pb208", "Pb208"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
      {"A anti-A", {"Pb208", "-1000822080"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
      {"anti-A A", {"-1000822080", "Pb208"}, {2510.0, 2510.0}, {BeamType::Nucleus, BeamType::Nucleus}, true},
  };
}

// Compute the complete ordered and charge-conjugated light-by-light beam matrix
std::vector<BeamCase> LightByLightEventBeams() {
  using gra::nuclear::BeamType;
  struct ChargedBeam {
    std::string name;
    std::string particle;
    double      energy  = 0.0;
    BeamType    type    = BeamType::Lepton;
    bool        nuclear = false;
  };
  const std::array<ChargedBeam, 2> proton   = {{
        {"p", "p+", 6500.0, BeamType::Proton, false},
        {"pbar", "p-", 6500.0, BeamType::Proton, false},
  }};
  const std::array<ChargedBeam, 2> electron = {{
      {"e-", "e-", 60.0, BeamType::Lepton, false},
      {"e+", "e+", 60.0, BeamType::Lepton, false},
  }};
  const std::array<ChargedBeam, 2> ion      = {{
           {"A", "Pb208", 2510.0, BeamType::Nucleus, true},
           {"anti-A", "-1000822080", 2510.0, BeamType::Nucleus, true},
  }};

  std::vector<BeamCase> output;
  const auto            append = [&output](const auto &left, const auto &right) {
    for (const auto &first : left) {
      for (const auto &second : right) {
        output.push_back({first.name + " " + second.name,
                          {first.particle, second.particle},
                          {first.energy, second.energy},
                          {first.type, second.type},
                          first.nuclear || second.nuclear});
      }
    }
  };
  append(proton, proton);
  append(electron, electron);
  append(electron, proton);
  append(proton, electron);
  append(proton, ion);
  append(ion, proton);
  append(electron, ion);
  append(ion, electron);
  append(ion, ion);
  return output;
}

// Load one complete checked-in UPC generator card through the production reader
nlohmann::json BaseCard() {
  const std::string path           = gra::aux::ResolveProjectPath("icepack/UPC/GAMMA/ATLAS_1832628/mumu/gencard.json");
  nlohmann::json    card           = nlohmann::json::parse(gra::aux::GetInputData(path));
  card["GENERIC"]["CORES"]         = 1;
  card["GENERIC"]["NEVENTS"]       = 0;
  card["SCATTERING"]["LOOPSCREEN"] = false;
  card["SCATTERING"]["NSTARS"]     = 0;
  card["RADIATIVE"]                = {{"ISR_QED", "none"}, {"FSR_QED", "none"}};
  return card;
}

// Write one temporary model tune with the requested adaptation mode
std::string FastAdaptationTune(bool fast = true) {
  const std::filesystem::path directory =
      gra::aux::ResolveProjectPath(std::string("tmp/nuclear_adaptation_") + (fast ? "fast" : "full"));
  std::filesystem::create_directories(directory);

  std::filesystem::copy(gra::ResolveModelTuneDir("TUNE0"), directory,
                        std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);

  auto numerics = nlohmann::json::parse(gra::aux::GetInputData(gra::ResolveModelDataFile("TUNE0", "NUMERICS.json")));
  numerics["NUMERICS_MC"]["fast_adaptation"] = fast;
  std::ofstream numerics_out(directory / "NUMERICS.json");
  if (!numerics_out.good()) { throw std::runtime_error("FastAdaptationTune: cannot write NUMERICS.json"); }
  numerics_out << numerics.dump(2);
  return directory.string();
}

// Configure low-cost but physical optical UPC quadratures for setup tests
void SetOpticalNuclearBlock(nlohmann::json &card) {
  card["NUCLEAR"] = {
      {"emd", false},
      {"fragmentation", nullptr},
      {"emission", {"coherent", "coherent"}},
      {"neutron_class", {"*", "*"}},
      {"photoproduction", {{"target", {"coherent", "coherent"}}, {"target_model", "glauber"}}},
      {"sigma_NN", nullptr},
      {"screening", "optical"},
      {"structure", "smooth"},
  };
}

// Configure one card for a supported beam pair and physical subprocess
nlohmann::json ProcessCard(const BeamCase &beam_case, const std::string &process) {
  nlohmann::json card           = BaseCard();
  card["SCATTERING"]["BEAM"]    = beam_case.beam;
  card["SCATTERING"]["ENERGY"]  = beam_case.energy;
  card["SCATTERING"]["PROCESS"] = process;
  card["GENCUTS"]["<C>"]["Rap"] = {-2.4, 2.4};
  card["GENCUTS"]["<F>"]["Rap"] = {-2.4, 2.4};
  if (beam_case.nuclear) {
    SetOpticalNuclearBlock(card);
  } else {
    card.erase("NUCLEAR");
  }
  return card;
}

// Read and prepare one process through the generator and process public APIs
void RequirePrepared(const BeamCase &beam_case, const std::string &process) {
  gra::MGraniitti      generator;
  const nlohmann::json card = ProcessCard(beam_case, process);
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE(generator.proc != nullptr);
  REQUIRE(generator.proc->HasNuclear() == beam_case.nuclear);
  REQUIRE_NOTHROW(generator.proc->PrepareRun());

  const auto initial = generator.proc->GetInitialState();
  REQUIRE(initial.size() == 2);
  const std::array<double, 2> full_energy = {generator.proc->state.lts.pbeam1.E(),
                                             generator.proc->state.lts.pbeam2.E()};
  for (const auto &i : indices(full_energy)) {
    double expected = beam_case.energy[i];
    if (beam_case.type[i] == gra::nuclear::BeamType::Nucleus) {
      expected *= static_cast<double>(gra::nuclear::DecodeNuclearPDG(initial[i].pdg).a);
    }
    CHECK(full_energy[i] == Approx(expected).epsilon(1.0e-12));
  }
  REQUIRE(generator.proc->state.lts.hamp.metadata.spin_basis != gra::ScreeningSpinBasis::Unset);
  REQUIRE(generator.proc->state.lts.hamp.metadata.spin_rows > 0);

  const auto &upc = generator.proc->state.lts.upc_model;
  if (beam_case.nuclear) {
    REQUIRE(upc != nullptr);
    CHECK(upc->Type(1) == beam_case.type[0]);
    CHECK(upc->Type(2) == beam_case.type[1]);
  } else {
    CHECK(upc == nullptr);
  }
}

// Replace the exact phase-space tag in one process selector
std::string ReplaceMode(const std::string &process, const std::string &from, const std::string &to) {
  const std::string token    = "<" + from + ">";
  const std::size_t position = process.find(token);
  if (position == std::string::npos) { throw std::invalid_argument("ReplaceMode: missing phase-space tag"); }
  return process.substr(0, position) + "<" + to + ">" + process.substr(position + token.size());
}

constexpr std::size_t EPA_EVENT_TRIALS = 64;

// Compute one component of a deterministic low-discrepancy sequence
double RadicalInverse(std::size_t index, const std::size_t base) {
  double value = 0.0;
  double scale = 1.0 / static_cast<double>(base);
  while (index > 0) {
    value += scale * static_cast<double>(index % base);
    index /= base;
    scale /= static_cast<double>(base);
  }
  return value;
}

// Compute one interior phase-space point for an event-level regression
std::vector<double> EPAEventPoint(const std::size_t dimension, const std::size_t trial, bool knockout = false) {
  const std::array<std::array<double, 8>, 6> seed = {{
      {0.82, 0.80, 0.13, 0.67, 0.50, 0.35, 0.42, 0.58},
      {0.76, 0.84, 0.21, 0.73, 0.45, 0.61, 0.38, 0.62},
      {0.88, 0.79, 0.31, 0.82, 0.55, 0.47, 0.46, 0.54},
      {0.72, 0.86, 0.17, 0.59, 0.40, 0.72, 0.34, 0.66},
      {0.85, 0.74, 0.39, 0.91, 0.60, 0.28, 0.44, 0.56},
      {0.79, 0.89, 0.07, 0.53, 0.48, 0.83, 0.36, 0.64},
  }};
  if (dimension > seed.front().size() || trial >= EPA_EVENT_TRIALS) {
    throw std::invalid_argument("EPAEventPoint: unsupported dimensions");
  }
  if (trial < seed.size()) { return {seed[trial].begin(), seed[trial].begin() + dimension}; }

  constexpr std::array<std::size_t, 8> base = {2, 3, 5, 7, 11, 13, 17, 19};
  std::vector<double>                  point(dimension);
  for (const auto &i : indices(point)) { point[i] = 0.02 + 0.96 * RadicalInverse(trial - seed.size() + 1, base[i]); }
  if (knockout) {
    // Cover nucleon recoil as well as the small coherent photon momentum on the logarithmic map
    const double target = 0.925 + 0.005 * (trial % 8), photon = 0.78 + 0.005 * ((trial / 2) % 5);
    point[0] = trial % 2 ? target : photon;
    point[1] = trial % 2 ? photon : target;
  }
  return point;
}

// Require one accepted event with a finite positive EPA weight
double RequirePositiveEPAEvent(gra::MProcess &process, const bool screening = false) {
  std::size_t kinematics_failure  = 0;
  std::size_t amplitude_failure   = 0;
  std::size_t technical_failure   = 0;
  std::size_t amplitude_rejection = 0;
  std::size_t fiducial_rejection  = 0;
  std::size_t veto_rejection      = 0;
  std::size_t valid_zero          = 0;
  std::string amplitude_message;
  const auto& upc = process.state.lts.upc_model;
  const bool knockout = process.state.nuclear_final.has_value() && upc &&
      std::any_of(upc->Param().target.begin(), upc->Param().target.end(), [](auto sector) {
        return sector != gra::nuclear::CoherenceType::Coherent;
      });
  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    gra::MEventWeightState event_state;
    event_state.include_screening = screening;
    const double weight           = process.EventWeight(EPAEventPoint(process.ProcPtr.LIPSDIM, trial, knockout), event_state);
    kinematics_failure += event_state.kinematics_ok ? 0 : 1;
    amplitude_failure += event_state.amplitude_failure ? 1 : 0;
    if (amplitude_message.empty()) { amplitude_message = event_state.amplitude_failure_message; }
    technical_failure += event_state.technical_failure ? 1 : 0;
    amplitude_rejection += event_state.amplitude_ok ? 0 : 1;
    fiducial_rejection += event_state.fidcuts_ok ? 0 : 1;
    veto_rejection += event_state.vetocuts_ok ? 0 : 1;
    valid_zero += event_state.Valid() && !(weight > 0.0) ? 1 : 0;
    if (event_state.Valid() && weight > 0.0) {
      REQUIRE(std::isfinite(weight));
      return weight;
    }
  }
  FAIL("No positive EPA event in deterministic interior phase-space points: "
       << "kinematics failures=" << kinematics_failure << ", amplitude failures=" << amplitude_failure
       << ", technical failures=" << technical_failure << ", amplitude rejections=" << amplitude_rejection
       << ", fiducial rejections=" << fiducial_rejection << ", veto rejections=" << veto_rejection
       << ", valid zero weights=" << valid_zero << ", first amplitude failure: " << amplitude_message);
  return 0.0;
}

// Check equal phase space normalization using the independently fitted photoproduction amplitudes
void RequirePhotoNormalization(const BeamCase &beam_case) {
  const std::array<std::string, 2> processes = {"ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}"};
  std::array<std::unique_ptr<gra::MGraniitti>, 2> generator;
  for (const auto &i : indices(generator)) {
    generator[i]                     = std::make_unique<gra::MGraniitti>();
    auto card                        = ProcessCard(beam_case, processes[i]);
    card["FIDCUTS"]["active"]        = false;
    card["SCATTERING"]["LOOPSCREEN"] = false;
    card["GENCUTS"]["<F>"]["M"]      = {3.0, 3.2};
    generator[i]->ReadInput(card);
    generator[i]->proc->PrepareRun();
  }
  REQUIRE(generator[0]->proc->ProcPtr.LIPSDIM == generator[1]->proc->ProcPtr.LIPSDIM);

  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    std::array<double, 2> weight{};
    bool                  valid = true;
    for (const auto &i : indices(generator)) {
      gra::MEventWeightState event_state;
      event_state.include_screening = false;
      weight[i] =
          generator[i]->proc->EventWeight(EPAEventPoint(generator[i]->proc->ProcPtr.LIPSDIM, trial), event_state);
      valid = valid && event_state.Valid() && weight[i] > 0.0 && std::isfinite(weight[i]);
    }
    if (valid) {
      const auto &direct = generator[0]->proc->state.lts;
      const auto &regge  = generator[1]->proc->state.lts;
      // M4Vec equality compares all components with an absolute tolerance of 1e-10 GeV
      REQUIRE(direct.pfinal == regge.pfinal);
      REQUIRE(direct.decaytree.size() == regge.decaytree.size());
      for (const auto &i : indices(direct.decaytree)) { REQUIRE(direct.decaytree[i].p4 == regge.decaytree[i].p4); }
      REQUIRE(direct.PS_active == regge.PS_active);
      REQUIRE(direct.DW.Integral() == Approx(regge.DW.Integral()).epsilon(1.0e-12));
      REQUIRE(direct.central_phase_space_generated_jacobian ==
              Approx(regge.central_phase_space_generated_jacobian).epsilon(1.0e-12));
      const double direct_norm = direct.hamp.metadata.amplitude_normalization * gra::SquaredNorm(direct.hamp);
      const double regge_norm  = regge.hamp.metadata.amplitude_normalization * gra::SquaredNorm(regge.hamp);
      const double ratio = weight[1] / weight[0], expected = regge_norm / direct_norm;
      CAPTURE(ratio, expected);
      CHECK(ratio == Approx(expected).epsilon(1.0e-12));
      return;
    }
  }
  throw std::runtime_error("RequirePhotoNormalization: no common positive phase-space point");
}

// Require one accepted photoproduction event after the requested screening
double RequirePositiveScreenedPhotoEvent(gra::MProcess &process) {
  const auto &upc = process.state.lts.upc_model;
  if (upc != nullptr) {
    REQUIRE(process.GetScreening());
    if (upc->HasHadronicPair()) { REQUIRE(upc->HadronicConvolution()); }
  } else {
    REQUIRE(process.GetScreening());
  }
  std::size_t kinematics_failure = 0;
  std::size_t amplitude_failure  = 0;
  std::size_t technical_failure  = 0;
  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    gra::MEventWeightState event_state;
    REQUIRE(event_state.include_screening);
    const double weight = process.EventWeight(EPAEventPoint(process.ProcPtr.LIPSDIM, trial), event_state);
    kinematics_failure += event_state.kinematics_ok ? 0 : 1;
    amplitude_failure += event_state.amplitude_failure ? 1 : 0;
    technical_failure += event_state.technical_failure ? 1 : 0;
    if (!event_state.Valid() || !(weight > 0.0)) { continue; }
    REQUIRE(std::isfinite(weight));
    REQUIRE_FALSE(process.state.lts.hamp.empty());
    REQUIRE(gra::AllFinite(process.state.lts.hamp));
    REQUIRE_FALSE(process.state.lts.screening.active);
    if (process.state.lts.upc_model != nullptr) { REQUIRE(process.state.upc_final.valid); }
    return weight;
  }
  FAIL(
      "No positive screened ygg[jpsi] event in deterministic phase-space "
      "points: "
      << "kinematics failures=" << kinematics_failure << ", amplitude failures=" << amplitude_failure
      << ", technical failures=" << technical_failure);
  return 0.0;
}

// Set a coherent sliding LS mixture in both ordered photoproduction vertices
void SetGPPhotoLS(gra::LORENTZSCALAR &lts) {
  REQUIRE_FALSE(lts.process.RESONANCES.empty());
  lts.process.DERIVATIVE_FACTOR = true;
  for (auto &[name, resonance] : lts.process.RESONANCES) {
    REQUIRE_FALSE(resonance.production.empty());
    for (auto &production : resonance.production) {
      auto &hel          = production.hel;
      hel.coupling_basis = gra::CouplingBasis::LS;
      hel.alpha_ls.Set(0, 2, {1.0, 0.2});
      hel.alpha_ls.Set(2, 2, {0.3, 0.1});
      const bool upper = production.tree[0].p.pdg == gra::PDG::PDG_gamma;
      gra::gpom::InitResonanceLS(hel, 1, upper ? 2 : 4, upper ? 4 : 2);
    }
  }
}

// Expose the real factorized kinematics and screening calls for boundary tests
class ScreenNodeProbe : public gra::MFactorized {
 public:
  // Copy one initialized physical process without replacing its amplitude
  explicit ScreenNodeProbe(const gra::MFactorized &source) : gra::MFactorized(source) {}
  using gra::MFactorized::B51BuildKin;
  using gra::MFactorized::LoopKinematics;
  using gra::MProcess::ComputeScreeningNode;
  using gra::MProcess::EvaluateBareAmplitude;
  using gra::MProcess::evaluation_status;
  using gra::MProcess::PrepareUPCEvent;
  using gra::MProcess::ScreenedAmplitudeSquared;
};

// Store one complete external photoproduction evaluation at a screening-loop point
struct PhotoNodeState {
  std::array<gra::M4Vec, 2>            q;
  std::vector<std::complex<double>>    amplitude;
  std::vector<gra::nuclear::PhotoTerm> photo_terms;
};

// Trace the real factorized photoproduction amplitude without replacing its physics API
class PhotoNodeProbe : public gra::MFactorized {
 public:
  // Copy one fully prepared factorized process and rebuild its worker runtime
  explicit PhotoNodeProbe(const gra::MFactorized &source) : gra::MFactorized(source) {}

  // Evaluate the complete nuclear screening convolution and retain every node
  double ScreenedAmp2(const bool screening = true) {
    trace_.clear();
    return ScreenedAmplitudeSquared(screening);
  }

  // Evaluate one independently shifted point through the same process runtime
  bool EvaluateNode(const std::array<double, 2> &upper, const std::array<double, 2> &lower) {
    trace_.clear();
    state.lts.pfinal_orig      = state.lts.pfinal;
    state.lts.screening.active = true;
    if (!LoopKinematics(upper, lower)) {
      state.lts.screening.active = false;
      return false;
    }
    try {
      EvaluateAmplitude();
      state.lts.screening.active = false;
      return true;
    } catch (...) {
      state.lts.screening.active = false;
      throw;
    }
  }

  // Compute the Born entry followed by every accepted shifted loop entry
  const std::vector<PhotoNodeState> &Trace() const { return trace_; }

 protected:
  // Evaluate the real photoproduction amplitude and copy its node-local source payload
  double EvaluateBareAmplitude() override {
    const double value = gra::MProcess::EvaluateBareAmplitude();
    trace_.push_back({{state.lts.q1, state.lts.q2}, state.lts.hamp, state.lts.screening.photo.term});
    return value;
  }

 private:
  std::vector<PhotoNodeState> trace_;
};

// Require one accepted event with finite nonzero weighted helicity amplitudes
double RequirePositiveLightByLightEvent(gra::MProcess &process) {
  std::size_t kinematics_failure = 0;
  std::size_t amplitude_failure  = 0;
  std::size_t technical_failure  = 0;
  std::size_t valid_zero         = 0;
  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    gra::MEventWeightState event_state;
    event_state.include_screening = false;
    const double weight           = process.EventWeight(EPAEventPoint(process.ProcPtr.LIPSDIM, trial), event_state);
    kinematics_failure += event_state.kinematics_ok ? 0 : 1;
    amplitude_failure += event_state.amplitude_failure ? 1 : 0;
    technical_failure += event_state.technical_failure ? 1 : 0;
    valid_zero += event_state.Valid() && !(weight > 0.0) ? 1 : 0;
    if (!event_state.Valid() || !(weight > 0.0)) { continue; }
    REQUIRE(std::isfinite(weight));
    REQUIRE_FALSE(process.state.lts.hamp.empty());
    REQUIRE(gra::AllFinite(process.state.lts.hamp));
    REQUIRE(gra::SquaredNorm(process.state.lts.hamp) > 0.0);
    REQUIRE(process.state.lts.exact_forward_photon_kinematics);
    return weight;
  }
  FAIL("No positive yy[yy] event in deterministic interior phase-space points: "
       << "kinematics failures=" << kinematics_failure << ", amplitude failures=" << amplitude_failure
       << ", technical failures=" << technical_failure << ", valid zero weights=" << valid_zero);
  return 0.0;
}

// Check the derived Glauber eigenstates against their first two moments
void CheckGlauberRule(const gra::nuclear::MUPC &upc) {
  const auto *glauber = upc.Glauber();
  REQUIRE(glauber != nullptr);
  const auto &param = upc.Param().glauber;
  REQUIRE(param.profile.sigma > 0.0);
  REQUIRE(glauber->SigmaCount() > 0);
  double mean   = 0.0;
  double second = 0.0;
  for (std::size_t i = 0; i < glauber->SigmaCount(); ++i) {
    for (std::size_t j = 0; j < glauber->SigmaCount(); ++j) {
      const double sigma = glauber->Sigma(i, j);
      mean += sigma;
      second += sigma * sigma;
    }
  }
  mean /= static_cast<double>(glauber->SigmaCount() * glauber->SigmaCount());
  second /= static_cast<double>(glauber->SigmaCount() * glauber->SigmaCount());
  const double omega = second / (mean * mean) - 1.0;
  CHECK(mean == Approx(param.profile.sigma).epsilon(2.0e-12));
  CHECK(omega == Approx(std::pow(1.0 + param.fluctuation.omega, 2) - 1.0).margin(2.0e-12));
  if (std::fpclassify(param.fluctuation.omega) == FP_ZERO) {
    CHECK(glauber->SigmaCount() == 1);
    CHECK(glauber->GGCFAmp(4.0) == Approx(glauber->OpticalAmp(4.0)).margin(2.0e-12));
  } else {
    CHECK(glauber->SigmaCount() > 1);
    CHECK(glauber->GGCFAmp(4.0) >= glauber->OpticalAmp(4.0));
  }
}

}  // namespace

// Preserve light central masses through the asymmetric beam-frame transformation
TEST_CASE("Asymmetric UPC phase spaces retain physical lepton mass shells", "[nuclear][EPA][kinematics][mass-shell]") {
  const auto input = nlohmann::json::parse(gra::aux::GetInputData(
      gra::aux::ResolveProjectPath("icepack/UPC/GAMMA/ALICE_2654315/mumu/gencard.json")));
  for (const bool reverse : {false, true}) {
    for (const std::string mode : {"F", "C"}) {
      DYNAMIC_SECTION(mode << " reversed=" << reverse) {
        auto card = input;
        card["GENERIC"]["CORES"] = 1;
        card["SCATTERING"]["PROCESS"] = "yy[EPA]<" + mode + "> -> mu+ mu-";
        card["SCATTERING"]["LOOPSCREEN"] = false;
        card["NUCLEAR"]["emd"] = false;
        card["RADIATIVE"]["FSR_QED"] = "none";
        const auto &pair = card["FIDCUTS"]["CENTRAL"]["[13,-13]"];
        card["GENCUTS"]["<F>"]["M"] = pair.at("M");
        card["GENCUTS"]["<C>"]["M"] = pair.at("M");
        card["GENCUTS"]["<F>"]["Rap"] = pair.at("Rap");
        if (reverse) {
          std::swap(card["SCATTERING"]["BEAM"][0], card["SCATTERING"]["BEAM"][1]);
          std::swap(card["SCATTERING"]["ENERGY"][0], card["SCATTERING"]["ENERGY"][1]);
          const auto range = card["GENCUTS"]["<F>"]["Rap"].get<std::array<double, 2>>();
          card["GENCUTS"]["<F>"]["Rap"] = {-range[1], -range[0]};
          card["FIDCUTS"]["CENTRAL"]["[13,-13]"]["Rap"] = {-range[1], -range[0]};
        }
        // These points have both daughters within the pair rapidity interval
        card["GENCUTS"]["<C>"]["Rap"] = card["GENCUTS"]["<F>"]["Rap"];
        gra::MGraniitti generator;
        generator.ReadInput(card);
        generator.proc->PrepareRun();
        for (const double rapidity : {0.2, 0.5, 0.8}) {
          for (unsigned int trial = 0; trial < 8; ++trial) {
            generator.proc->state.random.SetSeed(7103 + trial);
            const std::vector<double> point = mode == "F"
                ? std::vector<double>{0.75, 0.65, 0.2, 0.7, rapidity, 0.4}
                : std::vector<double>{0.75, 0.65, 0.2, 0.7, 0.8, 0.4, rapidity, rapidity};
            gra::MEventWeightState state;
            state.include_screening = false;
            const double weight = generator.proc->EventWeight(point, state);
            CAPTURE(rapidity, trial, state.kinematics_ok, state.technical_failure);
            REQUIRE(state.Valid());
            REQUIRE(weight > 0.0);
            const auto &lts = generator.proc->state.lts;
            REQUIRE(gra::math::CheckEMC(lts.decaytree[0].p4 + lts.decaytree[1].p4 - lts.pfinal[0]));
            for (const auto &branch : lts.decaytree) {
              const double scale = std::max(1.0, branch.p4.E() * branch.p4.E());
              REQUIRE(std::abs(branch.p4.M2() - branch.p.mass * branch.p.mass) <
                      64.0 * std::numeric_limits<double>::epsilon() * scale);
            }
          }
        }
      }
    }
  }
}

// Keep a rejected EPA phase-space node local to the convolution
TEST_CASE("Rejected screening nodes do not retain an event kinematics failure",
          "[EPA][screening][kinematics][screen-fix]") {
  for (const std::string amplitude : {"yy[EPA]<F> -> mu+ mu-", "yy[yy]<F> -> gamma gamma"}) {
    CAPTURE(amplitude);
    gra::MGraniitti generator;
    generator.ReadInput(ProcessCard(SupportedBeams()[0], amplitude));
    generator.proc->PrepareRun();
    ScreenNodeProbe process(*dynamic_cast<gra::MFactorized *>(generator.proc));
    auto           &lts = process.state.lts;
    for (auto &branch : lts.decaytree) { branch.m_offshell = branch.p.mass; }
    REQUIRE(process.B51BuildKin(0.1, 0.2, 0.35, 2.1, 0.0, 400.0, lts.beam1.mass, lts.beam2.mass));
    REQUIRE(process.ScreenedAmplitudeSquared(false) > 0.0);
    REQUIRE(lts.epa_hard.Ready());
    lts.pfinal_orig                   = lts.pfinal;
    lts.screening.active              = true;
    const std::array<double, 2> upper = {lts.pfinal[1].Px(), lts.pfinal[1].Py()};
    const std::array<double, 2> lower = {lts.pfinal[2].Px(), lts.pfinal[2].Py()};
    const double                s     = lts.s;
    for (const bool rejected_first : {false, true}) {
      for (const bool rejected : {rejected_first, !rejected_first}) {
        // Exercise the real EPA phase-space failure while preserving reconstructible four-momenta
        lts.s = rejected ? 0.0 : s;
        REQUIRE(process.ComputeScreeningNode(upper, lower) == !rejected);
        CHECK(process.ProcPtr.EvaluationStatus() == (rejected ? gra::mg5helas::EvaluationStatus::KinematicsFailure
                                                              : gra::mg5helas::EvaluationStatus::Success));
        CHECK(process.evaluation_status == gra::mg5helas::EvaluationStatus::Success);
      }
    }
    lts.s = s;
    REQUIRE(process.ComputeScreeningNode(upper, lower));
    const auto accepted = std::vector<std::complex<double>>(lts.hamp);
    REQUIRE_FALSE(process.ComputeScreeningNode({lts.sqrt_s, 0.0}, {-lts.sqrt_s, 0.0}));
    CHECK(process.evaluation_status == gra::mg5helas::EvaluationStatus::Success);
    REQUIRE(process.ComputeScreeningNode(upper, lower));
    for (const auto &h : indices(accepted)) {
      CHECK(std::abs(lts.hamp[h] - accepted[h]) <= 1.0e-12 * std::max(1.0, std::abs(accepted[h])));
    }
    lts.epa_hard.amplitude.front().value = std::numeric_limits<double>::quiet_NaN();
    REQUIRE_THROWS_AS(process.ComputeScreeningNode(upper, lower), gra::AmplitudeFailure);
  }
}

// Cross the elementary target threshold with conserved screening-loop kinematics
TEST_CASE("Nuclear photoproduction keeps zero target currents at shifted thresholds",
          "[nuclear][photoproduction][screening][kinematics][screen-fix]") {
  using namespace gra::nuclear;
  auto beams   = SupportedBeams()[2];
  beams.energy = {500.0, 500.0};
  for (const auto &channel : {std::string("jpsi"), std::string("Z")}) {
    CAPTURE(channel);
    gra::MGraniitti generator;
    auto            card                         = ProcessCard(beams, "ygg[" + channel + "]<F> -> mu+ mu-");
    card["FIDCUTS"]["active"]                    = false;
    card["NUCLEAR"]["structure"]                 = "nucleon";
    card["NUCLEAR"]["emission"]                  = {"inclusive", "inclusive"};
    card["NUCLEAR"]["photoproduction"]["target"] = {"inclusive", "inclusive"};
    generator.ReadInput(card);
    generator.proc->PrepareRun();
    ScreenNodeProbe process(*dynamic_cast<gra::MFactorized *>(generator.proc));
    auto           &lts = process.state.lts;
    for (auto &branch : lts.decaytree) { branch.m_offshell = branch.p.mass; }
    const double mass = lts.PDG.FindByPDG(channel == "jpsi" ? 443 : 23).mass;
    REQUIRE(process.B51BuildKin(0.18, 0.18, 0.0, gra::math::PI, 0.0, mass * mass, lts.beam1.mass, lts.beam2.mass));
    const auto   target = gra::flux::PhotoTargetMomentum(gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower));
    const double threshold = gra::math::pow2(mass + target.M());
    double       left      = channel == "jpsi" ? -7.0 : -4.0;
    double       right     = 0.0;
    for (unsigned int i = 0; i < 32; ++i) {
      const double rapidity = 0.5 * (left + right);
      REQUIRE(
          process.B51BuildKin(0.18, 0.18, 0.0, gra::math::PI, rapidity, mass * mass, lts.beam1.mass, lts.beam2.mass));
      if ((lts.q1 + target).M2() > threshold + 0.01) {
        right = rapidity;
      } else {
        left = rapidity;
      }
    }
    process.PrepareUPCEvent(false);
    REQUIRE(lts.upc_event->HasSamples());
    REQUIRE(process.EvaluateBareAmplitude() > 0.0);
    ScreenLayout layout;
    layout.type  = ScreenType::Photo;
    layout.photo = PhotoChannels(*lts.upc_event);
    // Capture the native amplitude and its event-local nuclear currents
    const auto point = [&]() {
      ScreenPoint output;
      output.amplitude     = lts.hamp;
      output.photo_current = lts.screening.photo.current;
      for (const auto &leg : indices(output.transfer)) {
        const auto state =
            gra::ResolveForwardLegState(lts, leg == 0 ? gra::ForwardBeamLeg::Upper : gra::ForwardBeamLeg::Lower);
        output.transfer[leg] = RestTransfer(state.incoming, state.transfer, state.emitter.mass);
      }
      return output;
    };
    const auto born = point();
    MUPCScreen screen(*lts.upc_event, layout, born);
    lts.pfinal_orig      = lts.pfinal;
    lts.screening.active = true;
    REQUIRE(process.ComputeScreeningNode({0.48, 0.0}, {-0.48, 0.0}));
    REQUIRE((lts.q1 + target).M2() < threshold);
    REQUIRE(lts.hamp.size() == born.amplitude.size());
    REQUIRE(lts.screening.photo.current[1].has_value());
    CHECK(lts.screening.photo.current[1]->sample.size() == lts.upc_event->Bank(2)->Size());
    CHECK(gra::SquaredNorm(lts.screening.photo.current[1]->sample) == Approx(0.0).margin(1.0e-30));
    for (const auto &h : indices(lts.hamp)) {
      if (layout.photo[h % layout.photo.size()].direction == PhotoDirection::Upper) {
        CHECK(std::abs(lts.hamp[h]) == Approx(0.0).margin(1.0e-30));
      }
    }
    REQUIRE(gra::SquaredNorm(lts.hamp) > 0.0);
    MUPCScreen closed(*lts.upc_event, layout, point());
    LoopNode   node;
    node.weight = -0.1;
    REQUIRE_NOTHROW(screen.Add(node, point()));
    CHECK(gra::Sum(screen.Result().helicity_norm) > 0.0);
    REQUIRE(process.ComputeScreeningNode({0.18, 0.0}, {-0.18, 0.0}));
    REQUIRE(lts.screening.photo.current[1].has_value());
    CHECK(gra::SquaredNorm(lts.screening.photo.current[1]->sample) > 0.0);
    REQUIRE_NOTHROW(screen.Add(node, point()));
    REQUIRE_NOTHROW(closed.Add(node, point()));
    CHECK(gra::Sum(closed.Result().helicity_norm) > 0.0);
  }
}

// Construct one explicit one-pair screening state for numerical checks
TEST_CASE("Nuclear screening accepts an explicit one-pair sample", "[nuclear][EPA][config][screening]") {
  gra::MGraniitti generator;
  auto            card             = BaseCard();
  card["SCATTERING"]["LOOPSCREEN"] = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  const auto &model = generator.proc->state.lts.upc_model;
  REQUIRE(model != nullptr);
  CHECK(model->Nucleus(1) == model->Nucleus(2));
  gra::MRandom random;
  random.SetSeed(314159);
  const auto event = model->Sample(random, true, 1);
  REQUIRE(event != nullptr);
  const std::array<std::size_t, 2> shape = {1, 1};
  CHECK(event->SampleShape() == shape);
  CHECK_FALSE(event->Nodes({}).empty());
}

// Keep coherent configuration work aligned with the common screening switch
TEST_CASE("Coherent nuclear photon fusion follows LOOPSCREEN steering", "[nuclear][EPA][config][convolution][event]") {
  for (const bool screening : {false, true}) {
    DYNAMIC_SECTION("LOOPSCREEN=" << screening) {
      gra::MGraniitti generator;
      auto            card             = BaseCard();
      card["FIDCUTS"]["active"]        = false;
      card["SCATTERING"]["LOOPSCREEN"] = screening;
      card["NUCLEAR"]["screening"]      = "mc_ggcf";
      card["NUCLEAR"]["structure"]     = "nucleon";
      // This test samples hadronic configurations without an external EMD response
      card["NUCLEAR"]["emd"] = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());

      bool        accepted           = false;
      std::size_t kinematics_failure = 0;
      std::size_t amplitude_failure  = 0;
      std::size_t technical_failure  = 0;
      std::size_t valid_zero         = 0;
      for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS && !accepted; ++trial) {
        gra::MEventWeightState state;
        const double weight = generator.proc->EventWeight(EPAEventPoint(generator.proc->ProcPtr.LIPSDIM, trial), state);
        kinematics_failure += state.kinematics_ok ? 0 : 1;
        amplitude_failure += state.amplitude_failure ? 1 : 0;
        technical_failure += state.technical_failure ? 1 : 0;
        valid_zero += state.Valid() && !(weight > 0.0) ? 1 : 0;
        accepted = state.Valid() && weight > 0.0;
      }
      INFO("kinematics failures=" << kinematics_failure << ", amplitude failures=" << amplitude_failure
                                  << ", technical failures=" << technical_failure << ", valid zero=" << valid_zero);
      REQUIRE(accepted);

      const auto &model = generator.proc->state.lts.upc_model;
      const auto &event = generator.proc->state.lts.upc_event;
      REQUIRE(model != nullptr);
      REQUIRE(event != nullptr);
      REQUIRE(model->Photon(1) != nullptr);
      REQUIRE(model->Photon(2) != nullptr);
      CHECK(model->Photo(1) == nullptr);
      CHECK(model->Photo(2) == nullptr);
      if (!screening) {
        CHECK(model->Glauber() == nullptr);
        CHECK_FALSE(model->HadronicConvolution());
        CHECK(event == model);
        CHECK_FALSE(event->HasSamples());
        CHECK(event->Nodes({}).empty());
      } else {
        CHECK(model->Glauber() != nullptr);
        CHECK(model->HadronicConvolution());
        CHECK(event != model);
        const std::array<std::size_t, 2> shape = {event->Param().config.count, event->Param().config.count};
        CHECK(event->SampleShape() == shape);
        CHECK_FALSE(event->Nodes({}).empty());
      }
    }
  }
}

// Apply LOOPSCREEN to hadronic nuclear beams and ignore hadronic survival for leptons
TEST_CASE("Nuclear beams use LOOPSCREEN steering", "[nuclear][EPA][convolution][steering]") {
  for (const auto &beam_case : SupportedBeams()) {
    if (!beam_case.nuclear) { continue; }
    DYNAMIC_SECTION(beam_case.name) {
      auto card = ProcessCard(beam_case, "yy[EPA]<F> -> mu+ mu-");

      gra::MGraniitti accepted;
      card["SCATTERING"]["LOOPSCREEN"] = false;
      REQUIRE_NOTHROW(accepted.ReadInput(card));
      REQUIRE(accepted.proc->HasNuclear());
      REQUIRE_FALSE(accepted.proc->state.lts.upc_model->Screening());
      REQUIRE_FALSE(accepted.proc->state.lts.upc_model->HadronicConvolution());

      gra::MGraniitti screened;
      card["SCATTERING"]["LOOPSCREEN"] = true;
      REQUIRE_NOTHROW(screened.ReadInput(card));
      const auto &upc = screened.proc->state.lts.upc_model;
      REQUIRE(upc != nullptr);
      CHECK(upc->Screening());
      CHECK(upc->HadronicConvolution() == upc->HasHadronicPair());
    }
  }
}

// Keep optical nuclear survival finite for analytic light by light scattering
TEST_CASE("Optical nuclear light by light retains positive event weights",
          "[nuclear][EPA][convolution][optical][event]") {
  gra::MGraniitti   generator;
  const std::string path     = gra::aux::ResolveProjectPath("icepack/UPC/GAMMA/ATLAS_1811464/gammagamma/gencard.json");
  auto              card     = nlohmann::json::parse(gra::aux::GetInputData(path));
  card["GENERIC"]["CORES"]   = 1;
  card["GENERIC"]["NEVENTS"] = 0;
  card["SCATTERING"]["LOOPSCREEN"] = true;
  card["NUCLEAR"]["screening"]      = "optical";
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());

  std::size_t kinematics_failure = 0;
  std::size_t amplitude_failure  = 0;
  std::size_t technical_failure  = 0;
  std::size_t valid_zero         = 0;
  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    gra::MEventWeightState event_state;
    const double           weight =
        generator.proc->EventWeight(EPAEventPoint(generator.proc->ProcPtr.LIPSDIM, trial), event_state);
    kinematics_failure += event_state.kinematics_ok ? 0 : 1;
    amplitude_failure += event_state.amplitude_failure ? 1 : 0;
    technical_failure += event_state.technical_failure ? 1 : 0;
    valid_zero += event_state.Valid() && !(weight > 0.0) ? 1 : 0;
    if (event_state.Valid() && weight > 0.0) {
      REQUIRE(std::isfinite(weight));
      return;
    }
  }
  FAIL("No positive optical yy[yy] event in deterministic phase-space points: "
       << "kinematics failures=" << kinematics_failure << ", amplitude failures=" << amplitude_failure
       << ", technical failures=" << technical_failure << ", valid zero weights=" << valid_zero);
}

// Prepare every hard photon amplitude with the full charged EPA class matrix
TEST_CASE("Hard photon amplitudes prepare pp, pA, AA, ep, eA and ee beams", "[nuclear][EPA][MG5][initial-state]") {
  for (const auto &process : HardPhotonProcesses()) {
    for (auto beam_case : SupportedBeams()) {
      // Use sufficient lepton collision energy for the heavy resonance channels
      if (beam_case.name == "ee") { beam_case.energy = {6500.0, 6500.0}; }
      DYNAMIC_SECTION(beam_case.name << " with " << process) { RequirePrepared(beam_case, process); }
    }
  }
}

// Reject unsupported QED beams before constructing an amplitude worker
TEST_CASE("Full QED requires proton beams during initialization", "[nuclear][QED][initial-state][validation]") {
  auto       beams     = SupportedBeams();
  const auto conjugate = OrderedChargeBeams();
  beams.insert(beams.end(), conjugate.begin(), conjugate.end());
  for (const std::string mode : {"F", "C"}) {
    for (const auto &beam : beams) {
      DYNAMIC_SECTION(mode << " with " << beam.name) {
        gra::MGraniitti generator;
        generator.ReadInput(ProcessCard(beam, "yy[QED]<" + mode + "> -> mu+ mu-"));
        if (beam.name == "pp") {
          REQUIRE_NOTHROW(generator.proc->PrepareRun());
          CHECK(generator.proc->state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonHelicity);
        } else {
          REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
        }
      }
    }
  }
}

// Reject every requested Pomeron fusion resonance on lepton and nuclear beams
TEST_CASE("Resonance beam restrictions fail before amplitude workers", "[nuclear][Regge][initial-state][validation]") {
  for (const std::string mode : {"F", "C"}) {
    for (const std::string model : {"MP", "XP", "GP", "TP"}) {
      for (const auto &beam : SupportedBeams()) {
        for (const std::string resonances : {"f0_980:1", "rho_770:1,f0_980:1"}) {
          DYNAMIC_SECTION(model << " " << mode << " " << beam.name << " " << resonances) {
            gra::MGraniitti generator;
            auto            card = ProcessCard(beam, model + "[RES]<" + mode + "> -> pi+ pi- @RES{" + resonances + "}");
            card["SCATTERING"]["RES"] = nlohmann::json::array();
            generator.ReadInput(card);
            if (beam.name == "pp") {
              REQUIRE_NOTHROW(generator.proc->PrepareRun());
            } else {
              REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
            }
          }
        }
      }
    }
  }
}

// Retain the implemented photo channels when applying resonance beam restrictions
TEST_CASE("Supported Regge photo channels still initialize", "[nuclear][Regge][initial-state][validation]") {
  for (const std::string mode : {"F", "C"}) {
    for (const std::string model : {"MP", "XP", "GP"}) {
      for (const std::string resonance : {"rho_770", "rho3_1690"}) {
        DYNAMIC_SECTION(model << " " << mode << " ep " << resonance) {
          RequirePrepared(SupportedBeams()[3], model + "[RES]<" + mode + "> -> pi+ pi- @RES{" + resonance + ":1}");
        }
      }
    }
    for (const std::string model : {"GP", "TP"}) {
      for (const auto &beam : {SupportedBeams()[1], SupportedBeams()[2], SupportedBeams()[4]}) {
        DYNAMIC_SECTION(model << " " << mode << " " << beam.name) {
          RequirePrepared(beam, model + "[RES]<" + mode + "> -> pi+ pi- @RES{rho_770:1}");
        }
      }
    }
  }
}

// Route every antiproton ordering through signed generic kT-EPA currents
TEST_CASE("Antiproton EPA continua use charged photon currents", "[nuclear][EPA][QED][event][charge]") {
  using gra::nuclear::BeamType;
  const std::array<BeamCase, 3> beams = {{
      {"p pbar", {"p+", "p-"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"pbar p", {"p-", "p+"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
      {"pbar pbar", {"p-", "p-"}, {6500.0, 6500.0}, {BeamType::Proton, BeamType::Proton}, false},
  }};

  for (const std::string mode : {"F", "C"}) {
    for (const auto &beam_case : beams) {
      DYNAMIC_SECTION(mode << " phase space with " << beam_case.name) {
        gra::MGraniitti generator;
        auto            card      = ProcessCard(beam_case, "yy[EPA]<" + mode + "> -> mu+ mu-");
        card["FIDCUTS"]["active"] = false;
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        CHECK(generator.proc->state.lts.hamp.metadata.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      }
    }
  }
}

// Evaluate real light-by-light events for all beam orders, charges and modes
TEST_CASE("Analytic light by light produces events for every charged EPA beam",
          "[nuclear][EPA][lbyl][event][initial-state]") {
  for (const std::string mode : {"F", "C"}) {
    for (const auto &beam_case : LightByLightEventBeams()) {
      DYNAMIC_SECTION(mode << " phase space with " << beam_case.name) {
        gra::MGraniitti generator;
        auto            card      = ProcessCard(beam_case, "yy[yy]<" + mode + "> -> gamma gamma");
        card["FIDCUTS"]["active"] = false;
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        REQUIRE(RequirePositiveLightByLightEvent(*generator.proc) > 0.0);

        const auto &lts = generator.proc->state.lts;
        REQUIRE(lts.id1 == gra::PDG::PDG_gamma);
        REQUIRE(lts.id2 == gra::PDG::PDG_gamma);
        REQUIRE(std::isfinite(lts.pdf_xf1));
        REQUIRE(std::isfinite(lts.pdf_xf2));
        REQUIRE(lts.pdf_xf1 > 0.0);
        REQUIRE(lts.pdf_xf2 > 0.0);

        HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
        REQUIRE(generator.proc->EventRecord(event));
        REQUIRE_FALSE(event.particles().empty());
      }
    }
  }
}

// Match inclusive event sectors against coherent and incoherent flux fractions
TEST_CASE("Inclusive nuclear light by light resolves physical EPA sectors", "[nuclear][EPA][lbyl][event][incoherent]") {
  gra::MGraniitti generator;
  auto            card        = ProcessCard(SupportedBeams()[2], "yy[yy]<F> -> gamma gamma");
  card["FIDCUTS"]["active"]   = false;
  card["NUCLEAR"]["emission"] = {"inclusive", "inclusive"};
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveLightByLightEvent(*generator.proc) > 0.0);

  const auto &lts      = generator.proc->state.lts;
  const auto &metadata = lts.hamp.layout;
  using gra::nuclear::CoherenceType;
  const auto coherent   = static_cast<std::uint8_t>(CoherenceType::Coherent);
  const auto incoherent = static_cast<std::uint8_t>(CoherenceType::Incoherent);
  REQUIRE(metadata.epa_sector_resolved);
  REQUIRE(metadata.epa_sector_count == std::array<std::uint8_t, 2>{2, 2});
  REQUIRE(metadata.epa_sector_type[0] == std::array<std::uint8_t, 2>{coherent, incoherent});
  REQUIRE(metadata.epa_sector_type[1] == std::array<std::uint8_t, 2>{coherent, incoherent});
  REQUIRE(lts.hamp.size() % 4 == 0);

  const auto                  upper      = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
  const auto                  lower      = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Lower);
  const std::array<double, 2> upper_flux = {
      gra::flux::NuclearPhotonFluxTransverse(upper, CoherenceType::Coherent).Trace(),
      gra::flux::NuclearPhotonFluxTransverse(upper, CoherenceType::Incoherent).Trace()};
  const std::array<double, 2> lower_flux = {
      gra::flux::NuclearPhotonFluxTransverse(lower, CoherenceType::Coherent).Trace(),
      gra::flux::NuclearPhotonFluxTransverse(lower, CoherenceType::Incoherent).Trace()};
  const double upper_total = upper_flux[0] + upper_flux[1];
  const double lower_total = lower_flux[0] + lower_flux[1];
  REQUIRE(upper_flux[0] > 0.0);
  REQUIRE(upper_flux[1] > 0.0);
  REQUIRE(lower_flux[0] > 0.0);
  REQUIRE(lower_flux[1] > 0.0);
  REQUIRE(gra::flux::ForwardPhotonFlux(upper) == Approx(upper_total).epsilon(2.0e-12));
  REQUIRE(gra::flux::ForwardPhotonFlux(lower) == Approx(lower_total).epsilon(2.0e-12));

  // The resolved source carries the signed Pb charge form factor through the
  // same public flux API used inside every shifted screening-loop evaluation
  auto diffraction_probe        = upper;
  diffraction_probe.emission    = CoherenceType::Coherent;
  diffraction_probe.t           = -0.18 * 0.18;
  diffraction_probe.qt          = 0.18;
  const auto diffraction_sector = gra::flux::ForwardPhotonFluxSectors(diffraction_probe);
  REQUIRE(diffraction_sector.size() == 1);
  REQUIRE(diffraction_sector[0].source.parallel < 0.0);
  REQUIRE(diffraction_sector[0].source.perpendicular == Approx(0.0));
  REQUIRE(diffraction_sector[0].density.parallel ==
          Approx(gra::math::pow2(diffraction_sector[0].source.parallel)).epsilon(2.0e-13));

  const std::size_t rows = metadata.epa_rows_per_sector;
  REQUIRE(rows == 2);
  const std::size_t lower_rows = rows * metadata.epa_sector_count[1];
  const std::size_t group_rows = rows * metadata.epa_sector_count[0] * rows * metadata.epa_sector_count[1];
  REQUIRE(lts.hamp.size() % group_rows == 0);

  std::array<double, 4> sector_norm = {};
  for (const auto &row : indices(lts.hamp)) {
    const std::size_t local        = row % group_rows;
    const std::size_t upper_sector = (local / lower_rows) / rows;
    const std::size_t lower_sector = (local % lower_rows) / rows;
    sector_norm[2 * upper_sector + lower_sector] += std::norm(lts.hamp[row]);
  }
  const double total_norm = gra::Sum(sector_norm);
  REQUIRE(total_norm > 0.0);
  CHECK(total_norm == Approx(gra::SquaredNorm(lts.hamp)).epsilon(2.0e-12));
  for (const auto &upper_sector : indices(upper_flux)) {
    for (const auto &lower_sector : indices(lower_flux)) {
      const std::size_t sector = 2 * upper_sector + lower_sector;
      REQUIRE(sector_norm[sector] > 0.0);
    }
  }

  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  const auto       selected = generator.proc->state.upc_final;
  REQUIRE(selected.valid);
  REQUIRE(generator.proc->EventRecord(event));
  const auto roots = gra::nuclear::ForwardParents(event);
  for (const auto &leg : indices(roots)) {
    CHECK(roots[leg]->pid() == 1000822080);
    CHECK(roots[leg]->status() == 3);
    CHECK(roots[leg]->end_vertex() == nullptr);
    CHECK(roots[leg]->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() ==
          gra::nuclear::CoherenceName(selected.leg[leg]));
  }
}

// Prepare charge-conjugated yy beams and every asymmetric beam order
TEST_CASE("Photon fusion prepares reversed and charge-conjugated beams", "[nuclear][EPA][charge][initial-state]") {
  for (const auto &beam_case : OrderedChargeBeams()) {
    DYNAMIC_SECTION(beam_case.name) { RequirePrepared(beam_case, "yy[EPA]<F> -> mu+ mu-"); }
  }
}

// Prepare charge-conjugated ygg beams and every physical asymmetric order
TEST_CASE("Photoproduction prepares reversed and charge-conjugated beams",
          "[nuclear][EPA][photoproduction][charge][initial-state]") {
  for (const auto &beam_case : OrderedChargeBeams()) {
    if (beam_case.type[0] == gra::nuclear::BeamType::Lepton && beam_case.type[1] == gra::nuclear::BeamType::Lepton) {
      continue;
    }
    DYNAMIC_SECTION(beam_case.name) { RequirePrepared(beam_case, "ygg[jpsi]<F> -> mu+ mu-"); }
  }
}

// Produce screened ygg events for every hadronic beam class and beam order
TEST_CASE("Photoproduction produces screened events for ordered charged beams",
          "[nuclear][EPA][photoproduction][charge][event][screening]") {
  for (const std::string mode : {"F", "C"}) {
    for (const auto &beam_case : PhotoEventBeams()) {
      DYNAMIC_SECTION(mode << " phase space with " << beam_case.name) {
        gra::MGraniitti generator;
        auto            card             = ProcessCard(beam_case, "ygg[jpsi]<" + mode + "> -> mu+ mu-");
        card["FIDCUTS"]["active"]        = false;
        card["SCATTERING"]["LOOPSCREEN"] = true;
        REQUIRE_NOTHROW(generator.ReadInput(card));
        if (beam_case.type[0] == gra::nuclear::BeamType::Lepton ||
            beam_case.type[1] == gra::nuclear::BeamType::Lepton) {
          REQUIRE_THROWS_WITH(generator.proc->PrepareRun(),
                              "MProcess::PrepareRun: LOOPSCREEN requires two hadronic beams");
          continue;
        }
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        REQUIRE(RequirePositiveScreenedPhotoEvent(*generator.proc) > 0.0);
        REQUIRE(gra::SquaredNorm(generator.proc->state.lts.hamp) > 0.0);

        HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
        REQUIRE(generator.proc->EventRecord(event));
        REQUIRE_FALSE(event.particles().empty());
      }
    }
  }
}

// Keep the charged EPA screening basis fixed across worker construction
TEST_CASE("Charged EPA workers preserve their screening layout", "[nuclear][EPA][screening][initial-state]") {
  for (const auto &beam_case : SupportedBeams()) {
    if (beam_case.name == "pp") { continue; }
    DYNAMIC_SECTION(beam_case.name) {
      gra::MGraniitti      generator;
      const nlohmann::json card = ProcessCard(beam_case, "yy[EPA]<F> -> mu+ mu-");
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      const auto before = generator.proc->state.lts.hamp.metadata;
      REQUIRE(before.spin_basis == gra::ScreeningSpinBasis::ProtonIdentity);
      REQUIRE(before.spin_rows == 1);

      REQUIRE_NOTHROW(generator.proc->ProcPtr.PrepareBareAmplitude(generator.proc->state.lts));
      const auto &after = generator.proc->state.lts.hamp.metadata;
      CHECK(after.spin_basis == before.spin_basis);
      CHECK(after.spin_rows == before.spin_rows);
      CHECK(after.forward_noflip == before.forward_noflip);
      CHECK(after.amplitude_normalization == Approx(before.amplitude_normalization).margin(2.0e-14));
    }
  }
}

// Prepare photoproduction whenever at least one physical hadron is available
TEST_CASE("Photoproduction prepares pp, pA, AA, ep and eA beams", "[nuclear][EPA][photoproduction][initial-state]") {
  const auto beam_cases = SupportedBeams();
  for (const auto &beam_case : beam_cases) {
    if (beam_case.name == "ee") { continue; }
    DYNAMIC_SECTION(beam_case.name) { RequirePrepared(beam_case, "ygg[jpsi]<F> -> mu+ mu-"); }
  }
}

// Produce direct and Regge photoproduction events in both ep beam orders
TEST_CASE("Direct and Regge photoproduction produce ordered ep events", "[EPA][photoproduction][HERA][event]") {
  using gra::nuclear::BeamType;
  const std::array<BeamCase, 2> beam_cases = {
      BeamCase{"ep", {"e-", "p+"}, {60.0, 6500.0}, {BeamType::Lepton, BeamType::Proton}, false},
      BeamCase{"pe", {"p+", "e-"}, {6500.0, 60.0}, {BeamType::Proton, BeamType::Lepton}, false}};
  const std::array<std::string, 2> processes = {"ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}"};

  for (const auto &process : processes) {
    for (const auto &beam_case : beam_cases) {
      DYNAMIC_SECTION(process << " with " << beam_case.name) {
        gra::MGraniitti generator;
        auto            card             = ProcessCard(beam_case, process);
        card["FIDCUTS"]["active"]        = false;
        card["SCATTERING"]["LOOPSCREEN"] = false;
        card["GENCUTS"]["<F>"]["M"]      = {3.0, 3.2};
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
        REQUIRE_FALSE(generator.proc->state.lts.hamp.empty());
        REQUIRE(gra::AllFinite(generator.proc->state.lts.hamp));
      }
    }
  }
}

// Compare generator weights and physical spin sums without equating independent model fits
TEST_CASE("Regge and direct photoproduction share phase space normalization",
          "[nuclear][EPA][photoproduction][normalization]") {
  for (const auto &beam_case : PhotoEventBeams()) {
    DYNAMIC_SECTION(beam_case.name) { RequirePhotoNormalization(beam_case); }
  }
}

// Check the absolute GP coupling through physical beam fluxes and the nuclear impulse limit
TEST_CASE("GP coupling scales physical pp ep pA eA AA amplitudes",
          "[nuclear][EPA][photoproduction][GP][normalization]") {
  for (const bool ls : {false, true}) {
    for (const auto &beam_case : PhotoEventBeams()) {
      DYNAMIC_SECTION(beam_case.name << " LS=" << ls) {
        gra::MGraniitti generator;
        auto            card             = ProcessCard(beam_case, "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
        card["FIDCUTS"]["active"]        = false;
        card["SCATTERING"]["LOOPSCREEN"] = false;
        card["GENCUTS"]["<F>"]["M"]      = {3.0, 3.2};
        if (beam_case.nuclear) { card["NUCLEAR"]["photoproduction"]["target_model"] = "impulse"; }
        generator.ReadInput(card);
        generator.proc->PrepareRun();
        if (ls) { SetGPPhotoLS(generator.proc->state.lts); }
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
        auto                      &lts       = generator.proc->state.lts;
        const double               original  = generator.proc->ProcPtr.GetBareAmplitude2(lts);
        const auto                 amplitude = lts.hamp;
        const std::complex<double> factor(1.6, 0.3);
        for (auto &[name, resonance] : lts.process.RESONANCES) {
          for (auto &production : resonance.production) {
            if (ls) {
              const auto terms = production.hel.alpha_ls;
              for (const auto &term : terms) {
                production.hel.alpha_ls.Set(term.l, term.two_s, term.coefficient * factor);
              }
            } else {
              production.hel.T *= factor;
            }
          }
        }
        const double scaled = generator.proc->ProcPtr.GetBareAmplitude2(lts);
        CHECK(scaled == Approx(std::norm(factor) * original).epsilon(1.0e-10));
        REQUIRE(lts.hamp.size() == amplitude.size());
        auto difference = lts.hamp;
        gra::AddScaled(difference, amplitude, -factor);
        CHECK(gra::SquaredNorm(difference) < 1.0e-20 * gra::SquaredNorm(amplitude));
      }
    }
  }
}

// Build only the active nuclear target service in lepton-ion photoproduction
TEST_CASE("Lepton-ion photoproduction requests only physical nuclear services",
          "[nuclear][EPA][photoproduction][service][event]") {
  using gra::nuclear::BeamType;
  const std::array<BeamCase, 2> beam_cases = {
      BeamCase{"eA", {"e-", "Pb208"}, {60.0, 2510.0}, {BeamType::Lepton, BeamType::Nucleus}, true},
      BeamCase{"Ae", {"Pb208", "e-"}, {2510.0, 60.0}, {BeamType::Nucleus, BeamType::Lepton}, true}};
  const std::array<std::string, 2> processes = {"ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}"};

  for (const auto &process : processes) {
    for (const auto &beam_case : beam_cases) {
      DYNAMIC_SECTION(process << " with " << beam_case.name) {
        gra::MGraniitti generator;
        auto            card             = ProcessCard(beam_case, process);
        card["FIDCUTS"]["active"]        = false;
        card["SCATTERING"]["LOOPSCREEN"] = false;
        card["GENCUTS"]["<F>"]["M"]      = {3.0, 3.2};
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        const auto &upc = generator.proc->state.lts.upc_model;
        REQUIRE(upc != nullptr);
        const int ion_leg = beam_case.type[0] == BeamType::Nucleus ? 1 : 2;
        CHECK(upc->Photon(ion_leg) == nullptr);
        REQUIRE(upc->Photo(ion_leg) != nullptr);
        CHECK(upc->Glauber() == nullptr);
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
        CHECK(generator.proc->state.lts.upc_event == upc);
      }
    }
  }
}

// Prepare eA and Ae inclusive ygg configuration-current ensembles
TEST_CASE("Lepton-ion ygg cards retain configuration currents in both orders",
          "[nuclear][EPA][photoproduction][good-walker][config]") {
  using gra::nuclear::BeamType;
  const std::array<BeamCase, 2> beam_cases = {
      BeamCase{"eA", {"e-", "Pb208"}, {60.0, 2510.0}, {BeamType::Lepton, BeamType::Nucleus}, true},
      BeamCase{"Ae", {"Pb208", "e-"}, {2510.0, 60.0}, {BeamType::Nucleus, BeamType::Lepton}, true}};

  for (const auto &beam_case : beam_cases) {
    DYNAMIC_SECTION(beam_case.name) {
      gra::MGraniitti   generator;
      auto              card                                = ProcessCard(beam_case, "ygg[jpsi]<F> -> mu+ mu-");
      const std::size_t ion_leg                             = beam_case.type[0] == BeamType::Nucleus ? 0U : 1U;
      card["SCATTERING"]["LOOPSCREEN"]                      = false;
      card["NUCLEAR"]["screening"]                           = "optical";
      card["NUCLEAR"]["structure"]                          = "nucleon";
      card["NUCLEAR"]["emission"]                           = {"coherent", "coherent"};
      card["NUCLEAR"]["photoproduction"]["target"]          = {"coherent", "coherent"};
      card["NUCLEAR"]["emission"][ion_leg]                  = "inclusive";
      card["NUCLEAR"]["photoproduction"]["target"][ion_leg] = "inclusive";

      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      const auto &upc = generator.proc->state.lts.upc_model;
      REQUIRE(upc != nullptr);
      REQUIRE_FALSE(upc->HasHadronicPair());
      REQUIRE(upc->Glauber() == nullptr);
      REQUIRE_FALSE(upc->HasSamples());
      gra::MRandom random;
      random.SetSeed(1618033);
      const auto event = upc->Sample(random);
      REQUIRE(event->HasSamples());
      REQUIRE(event->SampleCount() == event->Param().config.count);
      REQUIRE(event->Nodes({}).empty());

      const int        leg      = static_cast<int>(ion_leg + 1);
      const auto       emission = event->EmissionComponents(leg, {0.11, -0.04, 0.07});
      const std::array target   = {
            event->TargetRatios(leg, gra::nuclear::CoherenceType::Coherent, {0.11, -0.04, 0.07}),
            event->TargetRatios(leg, gra::nuclear::CoherenceType::Incoherent, {0.11, -0.04, 0.07})};
      for (const auto &sector : {emission, target}) {
        REQUIRE(sector[0].size() == event->SampleCount());
        REQUIRE(sector[1].size() == event->SampleCount());
        REQUIRE(gra::AllFinite(sector[0]));
        REQUIRE(gra::AllFinite(sector[1]));
      }
      REQUIRE_NOTHROW(generator.proc->ProcPtr.PrepareBareAmplitude(generator.proc->state.lts));
    }
  }
}

// Prepare the Tensor pion photoproduction process for pp and H1 e+p beams
TEST_CASE("Tensor pion photoproduction prepares pp and H1 e+p beams",
          "[nuclear][EPA][tensor][photoproduction][initial-state]") {
  const BeamCase pp = SupportedBeams()[0];
  const BeamCase h1 = {
      "H1 e+p", {"e+", "p+"}, {27.6, 920.0}, {gra::nuclear::BeamType::Lepton, gra::nuclear::BeamType::Proton}, false};
  for (const auto &beam_case : {pp, h1}) {
    DYNAMIC_SECTION(beam_case.name) { RequirePrepared(beam_case, "TP[PHOTO]<F> -> pi+ pi-"); }
  }
}

// Check Tensor resonance normalization and azimuthal covariance through physical beam sources
TEST_CASE("Tensor vector photoproduction preserves beam and coupling normalization",
          "[nuclear][tensor][photoproduction][normalization]") {
  for (const std::string resonance : {"rho_770", "phi_1020"}) {
    for (const auto &beam : SupportedBeams()) {
      if (beam.name == "ee") { continue; }
      DYNAMIC_SECTION(resonance << " " << beam.name) {
        const bool rho = resonance == "rho_770";
        auto       card =
            ProcessCard(beam, "TP[RES]<F> -> " + std::string(rho ? "pi+ pi-" : "K+ K-") + " @RES{" + resonance + ":1}");
        card["FIDCUTS"]["active"]   = false;
        card["GENCUTS"]["<F>"]["M"] = rho ? std::vector<double>{0.6, 0.9} : std::vector<double>{1.01, 1.03};
        if (beam.nuclear) { card["NUCLEAR"]["photoproduction"]["target_model"] = "impulse"; }
        gra::MGraniitti generator;
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
        auto        &lts       = generator.proc->state.lts;
        const double original  = generator.proc->ProcPtr.GetBareAmplitude2(lts);
        const auto   amplitude = lts.hamp;

        // Compare the common photon-current contraction with the pp Tensor source bank
        const auto          model = generator.proc->ProcPtr.GetModelTune();
        gra::MTensorPomeron tensor(lts, model,
                                   gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Resonance),
                                   gra::MTensorPomeronMode::Resonance);
        if (beam.name == "pp") {
          CHECK(tensor.PhotoProduction(lts, false) == Approx(original).epsilon(1.0e-10));
          REQUIRE(lts.hamp.size() == amplitude.size());
          auto difference = lts.hamp;
          gra::AddScaled(difference, amplitude, -1.0);
          CHECK(gra::SquaredNorm(difference) < 1.0e-20 * gra::SquaredNorm(amplitude));
        }

        auto rotated      = lts;
        rotated.amplitude = {};
        for (auto &p : rotated.pfinal) { p.RotateZ(0.73); }
        rotated.q1.RotateZ(0.73);
        rotated.q2.RotateZ(0.73);
        for (auto &branch : rotated.decaytree) { branch.p4.RotateZ(0.73); }
        CHECK(tensor.PhotoProduction(rotated, false) == Approx(original).epsilon(2.0e-7));

        // Use an exact binary rescaling in the boosted conserved-current contraction
        const double factor = 2.0;
        for (auto &[name, res] : lts.process.RESONANCES) {
          for (auto &channel : res.TP.channels) {
            for (auto &g : channel.g_tensor) { g *= factor; }
            for (auto& mixing : channel.VMD_MIXING) {
              for (auto& g : mixing.coupling) { g *= factor; }
            }
          }
        }
        CHECK(generator.proc->ProcPtr.GetBareAmplitude2(lts) == Approx(factor * factor * original).epsilon(1.0e-10));
      }
    }
  }
}

// Exercise the isoscalar nuclear meson continuum without assigning a vector absorption profile
TEST_CASE("Tensor nuclear meson currents preserve their Ward identity and vector interference",
          "[nuclear][tensor][photoproduction][continuum][gauge]") {
  const bool kaon      = GENERATE(false, true);
  const bool resonance = GENERATE(false, true);
  CAPTURE(kaon, resonance);
  const auto directory = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/tensor_nuclear_impulse"));
  std::filesystem::create_directories(directory);
  std::filesystem::copy(gra::ResolveModelTuneDir("TUNE0"), directory,
                        std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
  gra::MPDG pdg;
  pdg.ReadParticleData(gra::MPDG::DataFile(), directory.string());
  auto  general            = nlohmann::json::parse(gra::aux::GetInputData((directory / "GENERAL.json").string()));
  auto &photo              = general["PARAM_TENSORPOM"]["PHOTO"];
  photo["photon_exchange"] = false;
  auto exchanges           = photo["exchanges"].get<std::vector<int>>();
  exchanges.erase(
      std::remove_if(exchanges.begin(), exchanges.end(), [&](int code) { return pdg.FindByPDG(code).isospinX2 != 0; }),
      exchanges.end());
  photo["exchanges"] = exchanges;
  std::ofstream(directory / "GENERAL.json") << general.dump(2);
  for (const auto &beam : SupportedBeams()) {
    if (beam.name == "ee") { continue; }
    DYNAMIC_SECTION(beam.name) {
      auto card                     = ProcessCard(beam, "TP[PHOTO]<F> -> " + std::string(kaon ? "K+ K-" : "pi+ pi-") +
                                                            (resonance ? (kaon ? " @RES{phi_1020:1}" : " @RES{rho_770:1}") : ""));
      card["GENERIC"]["MODELPARAM"] = directory.string();
      card["FIDCUTS"]["active"]     = false;
      card["GENCUTS"]["<F>"]["M"]   = kaon ? std::vector<double>{1.01, 1.03} : std::vector<double>{0.6, 0.9};
      if (beam.nuclear) { card["NUCLEAR"]["photoproduction"]["target_model"] = "impulse"; }
      gra::MGraniitti generator;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      auto               &lts   = generator.proc->state.lts;
      const auto          model = generator.proc->ProcPtr.GetModelTune();
      gra::MTensorPomeron tensor(lts, model, gra::MTensorPomeron::ProcessDefinitionFor(gra::MTensorPomeronMode::Photo),
                                 gra::MTensorPomeronMode::Photo);
      gra::MTensorPhoto   pion(tensor, gra::ReadTensorPomeronParam(*model, lts.PDG));
      for (const bool upper : {false, true}) {
        if (!gra::flux::SupportsPhotoDirection(lts, upper ? 1 : 2)) { continue; }
        const auto           current = pion.DSCurrent(lts, upper, 0, 0);
        const auto          &q       = upper ? lts.q1 : lts.q2;
        std::complex<double> ward    = 0.0;
        double               scale   = 0.0;
        for (const auto &mu : indices(current)) {
          ward += q[mu] * current[mu];
          scale += std::abs(q[mu] * current[mu]);
        }
        CHECK(std::abs(ward) <= 1.0e-7 * scale);
      }
    }
  }
}

// Check that shared F and C selectors own identical physical definitions
TEST_CASE("F and C steering tags share photon-process definitions", "[nuclear][EPA][phase-space]") {
  gra::MFactorized factorized;
  gra::MCentral    central;

  for (const auto &selector : {std::array<std::string, 3>{"yy[EPA]<F>", "yy", "EPA"},
                               std::array<std::string, 3>{"ygg[jpsi]<F>", "ygg", "jpsi"}}) {
    const std::string central_selector = ReplaceMode(selector[0], "F", "C");
    INFO("Process selectors " << selector[0] << " and " << central_selector);
    REQUIRE(factorized.ProcPtr.ProcessExist(selector[0]));
    REQUIRE(central.ProcPtr.ProcessExist(central_selector));
    CHECK(factorized.ProcPtr.GetProcessDescriptor(selector[0]) ==
          central.ProcPtr.GetProcessDescriptor(central_selector));

    factorized.ProcPtr.Initialize(selector[1], selector[2]);
    central.ProcPtr.Initialize(selector[1], selector[2]);
    CHECK(factorized.ProcPtr.Processes() == central.ProcPtr.Processes());
  }
}

// Prepare actual nuclear cards through both shared central phase spaces
TEST_CASE("F and C steering tags both prepare nuclear photon processes", "[nuclear][EPA][phase-space]") {
  const BeamCase aa = SupportedBeams()[2];
  for (const std::string mode : {"F", "C"}) {
    DYNAMIC_SECTION("photon fusion " << mode) { RequirePrepared(aa, "yy[EPA]<" + mode + "> -> mu+ mu-"); }
    DYNAMIC_SECTION("photoproduction " << mode) { RequirePrepared(aa, "ygg[jpsi]<" + mode + "> -> mu+ mu-"); }
    DYNAMIC_SECTION("factorized flux " << mode) {
      gra::MGraniitti generator;
      auto            card      = ProcessCard(aa, "yy[FLUX]<" + mode + ">");
      card["FIDCUTS"]["active"] = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      CHECK(generator.proc->state.phase_space_class == mode);
      CHECK(generator.proc->state.lts.central_phase_space_mode == gra::CentralPhaseSpaceMode::Factorized);
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);

      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      REQUIRE(generator.proc->EventRecord(event));
      REQUIRE_FALSE(event.particles().empty());
    }
  }
}

// Record the sampled CI or IC target-incoherent photoproduction sector
TEST_CASE("Target-incoherent photoproduction records one unresolved ion",
          "[nuclear][EPA][photoproduction][event][record]") {
  gra::MGraniitti generator;
  auto            card                         = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
  card["FIDCUTS"]["active"]                    = false;
  card["NUCLEAR"]["emission"]                  = {"coherent", "coherent"};
  card["NUCLEAR"]["photoproduction"]["target"] = {"incoherent", "incoherent"};
  card["NUCLEAR"]["emd"]                       = false;
  card["NUCLEAR"]["neutron_class"]             = {"*", "*"};
  card["SCATTERING"]["LOOPSCREEN"]             = false;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);

  const auto selected = generator.proc->state.upc_final;
  REQUIRE(selected.valid);
  const int selected_incoherent =
      static_cast<int>(std::count(selected.leg.begin(), selected.leg.end(), gra::nuclear::CoherenceType::Incoherent));
  REQUIRE(selected_incoherent == 1);

  HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
  REQUIRE(generator.proc->EventRecord(event));
  const auto roots = gra::nuclear::ForwardParents(event);
  for (const auto &leg : indices(roots)) {
    CHECK(roots[leg]->status() == 3);
    CHECK(roots[leg]->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() ==
          gra::nuclear::CoherenceName(selected.leg[leg]));
  }
}

// Produce bare direct and Regge target-incoherent ion-ion events
TEST_CASE("Bare target-incoherent photoproduction produces ion-ion events",
          "[nuclear][EPA][photoproduction][event][incoherent]") {
  const BeamCase aa = SupportedBeams()[2];
  for (const std::string process : {"ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}"}) {
    DYNAMIC_SECTION(process) {
      gra::MGraniitti generator;
      auto            card                         = ProcessCard(aa, process);
      card["FIDCUTS"]["active"]                    = false;
      card["GENCUTS"]["<F>"]["M"]                  = {3.0, 3.2};
      card["NUCLEAR"]["structure"]                 = "nucleon";
      card["NUCLEAR"]["photoproduction"]["target"] = {"incoherent", "incoherent"};
      card["NUCLEAR"]["screening"]                  = "optical";
      card["SCATTERING"]["LOOPSCREEN"]             = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      const auto &event = generator.proc->state.lts.upc_event;
      REQUIRE(event != nullptr);
      for (int leg = 1; leg <= 2; ++leg) {
        REQUIRE(event->Bank(leg) != nullptr);
        CHECK(event->Bank(leg)->Size() == event->Param().current_count);
      }
      const std::array<std::size_t, 2> shape = {event->Bank(1)->Size(), event->Bank(2)->Size()};
      CHECK(event->SampleShape() == shape);
      CHECK(event->SampleCount() == shape[0] * shape[1]);
    }
  }
}

TEST_CASE("Screened optical photoproduction retains incoherent current banks",
          "[nuclear][EPA][photoproduction][event][incoherent][screening]") {
  gra::MGraniitti generator;
  auto            card                         = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
  card["FIDCUTS"]["active"]                    = false;
  card["GENCUTS"]["<F>"]["M"]                  = {3.0, 3.2};
  card["NUCLEAR"]["structure"]                 = "nucleon";
  card["NUCLEAR"]["photoproduction"]["target"] = {"incoherent", "incoherent"};
  card["NUCLEAR"]["screening"]                  = "optical";
  card["SCATTERING"]["LOOPSCREEN"]             = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveScreenedPhotoEvent(*generator.proc) > 0.0);
  REQUIRE(generator.proc->state.upc_final.valid);
  CHECK(generator.proc->state.upc_final.leg[0] != generator.proc->state.upc_final.leg[1]);
  double final_weight_sum = 0.0;
  for (const double weight : generator.proc->state.upc_final_weight) {
    REQUIRE(std::isfinite(weight));
    REQUIRE(weight >= 0.0);
    final_weight_sum += weight;
  }
  CHECK(final_weight_sum == Approx(1.0).margin(2.0e-14));
  const std::size_t selected =
      gra::nuclear::FinalPairIndex(generator.proc->state.upc_final.leg[0], generator.proc->state.upc_final.leg[1]);
  REQUIRE(generator.proc->state.upc_final_weight[selected] > 0.0);
  // Orthogonal incoherent sectors carry probability, not one coherent vector
  CHECK(std::fpclassify(gra::SquaredNorm(generator.proc->state.lts.hamp)) == FP_ZERO);

  const auto &event = generator.proc->state.lts.upc_event;
  REQUIRE(event != nullptr);
  REQUIRE(event->HasSamples());
  const std::array<std::size_t, 2> shape = {event->Param().current_count, event->Param().current_count};
  CHECK(event->SampleShape() == shape);
}

TEST_CASE("Nuclear photoproduction screening nodes equal independent shifted evaluations",
          "[nuclear][EPA][photoproduction][GP][tensor][screening][kinematics]") {
  const int  model = GENERATE(0, 1, 2);
  const bool ls    = model == 1;
  CAPTURE(ls);
  gra::MGraniitti generator;
  auto            card         = ProcessCard(SupportedBeams()[2],
                          model == 2 ? "TP[RES]<F> -> pi+ pi- @RES{rho_770:1}" : "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  card["FIDCUTS"]["active"]    = false;
  card["GENCUTS"]["<F>"]["M"]  = model == 2 ? std::vector<double>{0.6, 0.9} : std::vector<double>{3.0, 3.2};
  card["NUCLEAR"]["structure"] = "nucleon";
  card["NUCLEAR"]["photoproduction"]["target"] = {"inclusive", "inclusive"};
  card["NUCLEAR"]["screening"]                  = "optical";
  card["SCATTERING"]["LOOPSCREEN"]             = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  if (ls) { SetGPPhotoLS(generator.proc->state.lts); }
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);

  auto *master = dynamic_cast<gra::MFactorized *>(generator.proc);
  REQUIRE(master != nullptr);
  REQUIRE(master->state.lts.upc_event != nullptr);
  const std::array<double, 2>     upper    = {master->state.lts.pfinal[1].Px(), master->state.lts.pfinal[1].Py()};
  const std::array<double, 2>     lower    = {master->state.lts.pfinal[2].Px(), master->state.lts.pfinal[2].Py()};
  const std::array<gra::M3Vec, 2> transfer = {
      gra::nuclear::RestTransfer(master->state.lts.pbeam1, master->state.lts.q1, master->state.lts.beam1.mass),
      gra::nuclear::RestTransfer(master->state.lts.pbeam2, master->state.lts.q2, master->state.lts.beam2.mass)};
  const auto nodes = master->state.lts.upc_event->Nodes(transfer);
  REQUIRE(nodes.size() >= 3);

  PhotoNodeProbe screened(*master);
  if (model != 2) {
    REQUIRE(screened.state.lts.process.RESONANCES.at("jpsi").production.front().hel.UsesLSCouplings() == ls);
  }
  REQUIRE(screened.ScreenedAmp2() > 0.0);
  REQUIRE(screened.Trace().size() == nodes.size() + 1);
  const std::array<std::size_t, 3> selected          = {0, nodes.size() / 2, nodes.size() - 1};
  bool                             shifted_amplitude = false;
  bool                             shifted_current   = false;
  for (const std::size_t index : selected) {
    CAPTURE(index, nodes[index].kx, nodes[index].ky);
    PhotoNodeProbe direct(*master);
    REQUIRE(direct.EvaluateNode({upper[0] - nodes[index].kx, upper[1] - nodes[index].ky},
                                {lower[0] + nodes[index].kx, lower[1] + nodes[index].ky}));
    REQUIRE(direct.Trace().size() == 1);
    const auto &actual   = screened.Trace()[index + 1];
    const auto &expected = direct.Trace().front();
    for (std::size_t leg = 0; leg < 2; ++leg) {
      for (std::size_t mu = 0; mu < 4; ++mu) {
        const double scale = std::max({1.0, std::abs(actual.q[leg][mu]), std::abs(expected.q[leg][mu])});
        REQUIRE(std::abs(actual.q[leg][mu] - expected.q[leg][mu]) < 2.0e-12 * scale);
      }
    }
    REQUIRE(actual.amplitude.size() == expected.amplitude.size());
    for (const auto &h : indices(actual.amplitude)) {
      const double scale = std::max({1.0, std::abs(actual.amplitude[h]), std::abs(expected.amplitude[h])});
      REQUIRE(std::abs(actual.amplitude[h] - expected.amplitude[h]) < 3.0e-11 * scale);
      shifted_amplitude =
          shifted_amplitude || std::abs(actual.amplitude[h] - screened.Trace().front().amplitude[h]) > 1.0e-12 * scale;
    }
    REQUIRE(actual.photo_terms.size() == expected.photo_terms.size());
    for (const auto &term : indices(actual.photo_terms)) {
      REQUIRE(actual.photo_terms[term].direction == expected.photo_terms[term].direction);
      REQUIRE(actual.photo_terms[term].amplitude.size() == expected.photo_terms[term].amplitude.size());
      for (const auto &h : indices(actual.photo_terms[term].amplitude)) {
        const auto   value     = actual.photo_terms[term].amplitude[h];
        const auto   reference = expected.photo_terms[term].amplitude[h];
        const double scale     = std::max({1.0, std::abs(value), std::abs(reference)});
        REQUIRE(std::abs(value - reference) < 3.0e-11 * scale);
      }
      for (std::size_t leg = 0; leg < 2; ++leg) {
        const auto &value     = actual.photo_terms[term].photo_current[leg];
        const auto &reference = expected.photo_terms[term].photo_current[leg];
        REQUIRE(value.has_value() == reference.has_value());
        if (!value.has_value()) { continue; }
        REQUIRE(value->sample.size() == reference->sample.size());
        for (const auto &sample : indices(value->sample)) {
          const double scale = std::max({1.0, std::abs(value->sample[sample]), std::abs(reference->sample[sample])});
          REQUIRE(std::abs(value->sample[sample] - reference->sample[sample]) < 3.0e-11 * scale);
          const auto &born = screened.Trace().front().photo_terms[term].photo_current[leg];
          REQUIRE(born.has_value());
          REQUIRE(born->sample.size() == value->sample.size());
          const double born_scale = std::max({1.0, std::abs(value->sample[sample]), std::abs(born->sample[sample])});
          shifted_current =
              shifted_current || std::abs(value->sample[sample] - born->sample[sample]) > 1.0e-12 * born_scale;
        }
      }
    }
  }
  REQUIRE(shifted_amplitude);
  REQUIRE(shifted_current);
}

TEST_CASE("Screened nuclear GP keeps every resonance target current node local",
          "[nuclear][EPA][photoproduction][GP][screening][profile]") {
  gra::MGraniitti generator;
  auto            card         = ProcessCard(SupportedBeams()[2], "GP[RES]<F> -> mu+ mu- @RES{jpsi:1,psi2S:1}");
  card["FIDCUTS"]["active"]    = false;
  card["GENCUTS"]["<F>"]["M"]  = {3.0, 3.8};
  card["NUCLEAR"]["structure"] = "nucleon";
  card["NUCLEAR"]["photoproduction"]["target"] = {"inclusive", "inclusive"};
  card["NUCLEAR"]["screening"]                  = "optical";
  card["SCATTERING"]["LOOPSCREEN"]             = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);

  REQUIRE(generator.proc->ProcPtr.GetBareAmplitude2(generator.proc->state.lts) > 0.0);
  const auto &lts = generator.proc->state.lts;
  REQUIRE(lts.screening.photo.term.size() == 4);
  for (const auto &current : lts.screening.photo.current) { REQUIRE_FALSE(current.has_value()); }
  std::vector<std::complex<double>> sum(lts.hamp.size(), 0.0);
  for (const auto &term : lts.screening.photo.term) {
    REQUIRE(term.amplitude.size() == lts.hamp.size());
    REQUIRE(gra::AllFinite(term.amplitude));
    gra::AddScaled(sum, term.amplitude, 1.0);
    const std::size_t currents = static_cast<std::size_t>(term.photo_current[0].has_value()) +
                                 static_cast<std::size_t>(term.photo_current[1].has_value());
    REQUIRE(currents == 1);
  }
  for (const auto &h : indices(sum)) {
    const double scale = std::max({1.0, std::abs(sum[h]), std::abs(lts.hamp[h])});
    REQUIRE(std::abs(sum[h] - lts.hamp[h]) < 2.0e-12 * scale);
  }

  bool distinct_profile_current = false;
  for (std::size_t first = 0; first < lts.screening.photo.term.size(); ++first) {
    for (std::size_t second = first + 1; second < lts.screening.photo.term.size(); ++second) {
      for (std::size_t leg = 0; leg < 2; ++leg) {
        const auto &left  = lts.screening.photo.term[first].photo_current[leg];
        const auto &right = lts.screening.photo.term[second].photo_current[leg];
        if (!left.has_value() || !right.has_value() || left->sample.size() != right->sample.size()) { continue; }
        double difference = 0.0;
        for (const auto &sample : indices(left->sample)) {
          difference += std::abs(left->sample[sample] - right->sample[sample]);
        }
        distinct_profile_current = distinct_profile_current || difference > 1.0e-10;
      }
    }
  }
  REQUIRE(distinct_profile_current);
  auto *factorized = dynamic_cast<gra::MFactorized *>(generator.proc);
  REQUIRE(factorized != nullptr);
  PhotoNodeProbe screened(*factorized);
  REQUIRE(screened.ScreenedAmp2() > 0.0);
}

// Produce one bare event from each checked-in ALICE GP generator card
TEST_CASE("ALICE GP cards produce bare ion-ion events", "[nuclear][EPA][photoproduction][event][ALICE]") {
  for (const std::string path : {
           "icepack/UPC/PHOTOPROD/ALICE_1782227/rho_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_1840600/jpsi_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_2658375/jpsi_incoherent/gencard.json",
       }) {
    for (const std::string model : {"impulse", "glauber", "LTA"}) {
      DYNAMIC_SECTION(path << ", " << model) {
        gra::MGraniitti   generator;
        const std::string resolved                         = gra::aux::ResolveProjectPath(path);
        auto              card                             = nlohmann::json::parse(gra::aux::GetInputData(resolved));
        card["GENERIC"]["CORES"]                           = 1;
        card["GENERIC"]["NEVENTS"]                         = 0;
        card["SCATTERING"]["LOOPSCREEN"]                   = false;
        card["NUCLEAR"]["photoproduction"]["target_model"] = model;
        if (model == "LTA" && path.find("rho_coherent") != std::string::npos) {
          REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
          continue;
        }
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);

        auto *master = dynamic_cast<gra::MFactorized *>(generator.proc);
        REQUIRE(master != nullptr);
        gra::MFactorized worker(*master);
        worker.state.random.SetSeed(2658375);
        REQUIRE(RequirePositiveEPAEvent(worker) > 0.0);
      }
    }
  }
}

// Produce one bare event from each HEPData backed LHC UPC vector meson card
TEST_CASE("HEPData UPC vector meson cards retain positive event support",
          "[nuclear][EPA][photoproduction][event][HEPData]") {
  for (const std::string path : {
           "icepack/UPC/PHOTOPROD/ALICE_1782227/rho_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_1840601/jpsi_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_1840601/psi2s_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_2666011/jpsi_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/ATLAS_2966819/jpsi_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/CMS_2899343/jpsi_incoherent/gencard.json",
           "icepack/UPC/PHOTOPROD/CMS_2908607/phi_coherent/gencard.json",
           "icepack/UPC/PHOTOPROD/CMS_3140573/upsilon_coherent/gencard.json",
       }) {
    DYNAMIC_SECTION(path) {
      gra::MGraniitti   generator;
      const std::string resolved       = gra::aux::ResolveProjectPath(path);
      auto              card           = nlohmann::json::parse(gra::aux::GetInputData(resolved));
      card["GENERIC"]["CORES"]         = 1;
      card["GENERIC"]["NEVENTS"]       = 0;
      card["SCATTERING"]["LOOPSCREEN"] = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
    }
  }
}

// Prepare and evaluate the asymmetric pPb processes from both ALICE cards
TEST_CASE("ALICE pPb UPC cards retain positive event support", "[nuclear][EPA][pA][ALICE]") {
  for (const std::string path : {
           "icepack/UPC/GAMMA/ALICE_2654315/mumu/gencard.json",
           "icepack/UPC/PHOTOPROD/ALICE_2654315/jpsi/gencard.json",
       }) {
    DYNAMIC_SECTION(path) {
      gra::MGraniitti   generator;
      const std::string resolved       = gra::aux::ResolveProjectPath(path);
      auto              card           = nlohmann::json::parse(gra::aux::GetInputData(resolved));
      card["GENERIC"]["CORES"]         = 1;
      card["GENERIC"]["NEVENTS"]       = 0;
      card["SCATTERING"]["LOOPSCREEN"] = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(generator.proc->HasNuclear());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
    }
  }
}

// Preserve the fixed coherent final state under optical and GGCF survival
TEST_CASE("Coherent photoproduction records intact optical and GGCF ions",
          "[nuclear][EPA][photoproduction][event][record][ggcf]") {
  for (const std::string survival : {"optical", "optical_ggcf"}) {
    DYNAMIC_SECTION(survival) {
      gra::MGraniitti generator;
      auto            card             = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
      card["FIDCUTS"]["active"]        = false;
      card["NUCLEAR"]["screening"]      = survival;
      card["NUCLEAR"]["emd"]           = false;
      card["NUCLEAR"]["neutron_class"] = {"*", "*"};
      card["SCATTERING"]["LOOPSCREEN"] = true;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      const auto selected = generator.proc->state.upc_final;
      REQUIRE(selected.valid);
      CHECK(selected.leg[0] == gra::nuclear::CoherenceType::Coherent);
      CHECK(selected.leg[1] == gra::nuclear::CoherenceType::Coherent);

      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      REQUIRE(generator.proc->EventRecord(event));
      for (const auto &root : gra::nuclear::ForwardParents(event)) {
        CHECK(root->status() == 3);
        CHECK(root->end_vertex() == nullptr);
      }
    }
  }
}

// Construct tagged incoherent production with its calculated hard nuclear response
TEST_CASE("Neutron tags derive leg-specific EPA photoabsorption controls", "[nuclear][EPA][breakup][steering]") {
  const BeamCase  pa = SupportedBeams()[1];
  gra::MGraniitti generator;
  auto            card             = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"]           = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["neutron_class"] = {"*", "n == 0"};
  REQUIRE_NOTHROW(generator.ReadInput(card));
  const auto &upc = generator.proc->state.lts.upc_model;
  REQUIRE(upc != nullptr);
  const auto &breakup = upc->Param().breakup[1];
  CHECK(breakup.a == 208);
  CHECK(breakup.z == 82);
  CHECK(breakup.z_emit == Approx(1.0).margin(2.0e-14));
  CHECK(breakup.gamma > 1.0);
  CHECK(breakup.emitter == gra::nuclear::EmitterType::Proton);
  CHECK(breakup.proton.q_max > 0.0);
  CHECK(breakup.proton.q_nodes > 0);
  CHECK(breakup.photo.energy_min > 0.0);
  CHECK(breakup.photo.continuum.threshold > breakup.photo.energy_min);
  CHECK(breakup.transfer.match == Approx(200.0));
  CHECK(breakup.response.nodes == upc->Param().emd.isotope.front().response.nodes);
  CHECK(upc->Breakup(2)->Mean(20.0) > 0.0);
  CHECK(upc->Breakup(2)->Mean(20.0) > upc->Breakup(2)->Mean(100.0));
}

// Check pPb generation rapidity steering directly in the laboratory frame
TEST_CASE("Asymmetric pPb cards steer laboratory rapidity", "[nuclear][EPA][rapidity][pA]") {
  struct RapidityCase {
    std::string path;
    std::string process;
    bool        direct = false;
  };
  const std::array<RapidityCase, 3> cases = {{
      {"icepack/UPC/GAMMA/ALICE_2654315/mumu/gencard.json", "yy[EPA]<F> -> mu+ mu-", false},
      {"icepack/UPC/GAMMA/ALICE_2654315/mumu/gencard.json", "yy[EPA]<C> -> mu+ mu-", true},
      {"icepack/UPC/PHOTOPROD/ALICE_2654315/jpsi/gencard.json", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}", false},
  }};
  for (const auto &item : cases) {
    DYNAMIC_SECTION(item.process) {
      gra::MGraniitti generator;
      auto            card     = nlohmann::json::parse(gra::aux::GetInputData(gra::aux::ResolveProjectPath(item.path)));
      card["GENERIC"]["CORES"] = 1;
      card["GENERIC"]["NEVENTS"]    = 0;
      card["SCATTERING"]["PROCESS"] = item.process;
      card["NUCLEAR"]["emd"] = false;
      card["RADIATIVE"]["FSR_QED"] = "none";
      const auto &central = card["FIDCUTS"]["CENTRAL"];
      const auto &pair = central.at(central.contains("SYSTEM") ? "SYSTEM" : "[13,-13]");
      const auto range = pair.at("Rap").get<std::array<double, 2>>();
      for (const auto mode : {"<F>", "<C>"}) {
        card["GENCUTS"][mode]["M"] = pair.at("M");
        card["GENCUTS"][mode]["Rap"] = range;
      }
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      const double rapidity = generator.proc->state.lts.pfinal[0].Rap();
      CHECK(rapidity >= range[0]);
      CHECK(rapidity <= range[1]);
      if (item.direct) {
        REQUIRE(generator.proc->state.lts.decaytree.size() == 2);
        for (const auto &branch : generator.proc->state.lts.decaytree) {
          CHECK(branch.p4.Rap() >= range[0]);
          CHECK(branch.p4.Rap() <= range[1]);
        }
      }
    }
  }
}

// Parse analytic and configuration GGCF controls through production cards
TEST_CASE("Nuclear cards configure Gamma Glauber fluctuations", "[nuclear][EPA][glauber][steering]") {
  const BeamCase pa = SupportedBeams()[1];

  SECTION("analytic GGCF") {
    gra::MGraniitti generator;
    auto            card             = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["NUCLEAR"]["screening"]      = "optical_ggcf";
    card["NUCLEAR"]["sigma_NN"]      = 70.0;
    REQUIRE_NOTHROW(generator.ReadInput(card));
    const auto &upc = generator.proc->state.lts.upc_model;
    REQUIRE(upc != nullptr);
    REQUIRE(upc->Glauber() != nullptr);
    CHECK(upc->Param().survival == gra::nuclear::SurvivalType::OpticalGGCF);
    CHECK(upc->Param().glauber.profile.sigma == Approx(70.0));
    CHECK(upc->Param().survival_eikonal == upc->Param().ggcf.eikonal);
    CHECK(upc->Param().glauber.profile.omega > 0.0);
    CHECK(generator.proc->state.model_tune->Soft()->ActiveModel() ==
          generator.proc->state.model_tune->General("PARAM_SOFT").at("active_model").get<std::string>());
    CheckGlauberRule(*upc);
  }

  SECTION("configuration GGCF") {
    gra::MGraniitti generator;
    auto            card             = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["NUCLEAR"]["screening"]      = "mc_ggcf";
    card["NUCLEAR"]["structure"]     = "nucleon";
    card["NUCLEAR"]["sigma_NN"]      = 70.0;
    REQUIRE_NOTHROW(generator.ReadInput(card));
    const auto &upc = generator.proc->state.lts.upc_model;
    REQUIRE(upc != nullptr);
    REQUIRE(upc->Glauber() != nullptr);
    CHECK(upc->Param().survival == gra::nuclear::SurvivalType::MCGGCF);
    CHECK(upc->Param().config.count > 0);
    CHECK(upc->Param().config.sweeps > 0);
    CHECK(upc->Param().glauber.profile.sigma == Approx(70.0));
    CHECK(upc->Param().glauber.profile.b_node.size() > 2);
    CHECK(upc->Param().glauber.profile.b_node.size() == upc->Param().glauber.profile.inelastic.size());
    CHECK_FALSE(upc->Param().glauber.profile.fingerprint.empty());
    CHECK(upc->Param().survival_eikonal == upc->Param().ggcf.eikonal);
    CHECK(upc->Param().glauber.profile.omega > 0.0);
    CHECK(generator.proc->state.model_tune->Soft()->ActiveModel() ==
          generator.proc->state.model_tune->General("PARAM_SOFT").at("active_model").get<std::string>());
    CheckGlauberRule(*upc);
    CHECK_FALSE(upc->HasSamples());
    gra::MRandom random;
    random.SetSeed(314159);
    const auto event = upc->Sample(random);
    CHECK(event->SampleCount() == upc->Param().config.count);
  }

  SECTION("AA optical NN override and inactive rejection") {
    gra::MGraniitti generator;
    auto            card             = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["NUCLEAR"]["sigma_NN"]      = 70.0;
    REQUIRE_NOTHROW(generator.ReadInput(card));
    REQUIRE(generator.proc->state.lts.upc_model->Glauber() != nullptr);
    CHECK(generator.proc->state.lts.upc_model->Param().glauber.profile.sigma == Approx(70.0));

    gra::MGraniitti inactive;
    card                             = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["SCATTERING"]["LOOPSCREEN"] = false;
    card["NUCLEAR"]["sigma_NN"]      = 70.0;
    CHECK_THROWS_AS(inactive.ReadInput(card), std::invalid_argument);
  }

  SECTION("misplaced numerical configuration control") {
    gra::MGraniitti generator;
    auto            card               = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["screening"]        = "mc_ggcf";
    card["NUCLEAR"]["structure"]       = "nucleon";
    card["NUCLEAR"]["config_samples"]  = 1;
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("misplaced Glauber model control") {
    gra::MGraniitti generator;
    auto            card         = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["screening"]  = "mc_ggcf";
    card["NUCLEAR"]["structure"] = "nucleon";
    card["NUCLEAR"]["omega"]     = 0.2;
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
}

// Reject beam and subprocess combinations outside the implemented EPA domain
TEST_CASE("Invalid EPA initial states are rejected during setup", "[nuclear][EPA][initial-state]") {
  SECTION("pion emitter") {
    BeamCase        invalid{"pi p",
                     {"pi+", "p+"},
                     {100.0, 6500.0},
                     {gra::nuclear::BeamType::Proton, gra::nuclear::BeamType::Proton},
                     false};
    gra::MGraniitti generator;
    const auto      card = ProcessCard(invalid, "yy[EPA]<F> -> mu+ mu-");
    REQUIRE_NOTHROW(generator.ReadInput(card));
    REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
  }

  SECTION("photoproduction without a hadronic target") {
    const BeamCase  ee = SupportedBeams()[5];
    gra::MGraniitti generator;
    const auto      card = ProcessCard(ee, "ygg[jpsi]<F> -> mu+ mu-");
    REQUIRE_NOTHROW(generator.ReadInput(card));
    REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
  }

  SECTION("nonphoton nuclear subprocess") {
    const BeamCase  aa = SupportedBeams()[2];
    gra::MGraniitti generator;
    const auto      card = ProcessCard(aa, "MP[CON]<F> -> pi+ pi-");
    REQUIRE_NOTHROW(generator.ReadInput(card));
    REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
  }

  SECTION("missing nuclear steering") {
    const BeamCase  aa = SupportedBeams()[2];
    gra::MGraniitti generator;
    auto            card = ProcessCard(aa, "yy[EPA]<F> -> mu+ mu-");
    card.erase("NUCLEAR");
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("inactive nuclear steering without ions") {
    const BeamCase  pp = SupportedBeams()[0];
    gra::MGraniitti generator;
    auto            card = ProcessCard(pp, "yy[EPA]<F> -> mu+ mu-");
    SetOpticalNuclearBlock(card);
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("invalid nuclear steering") {
    const BeamCase  aa = SupportedBeams()[2];
    gra::MGraniitti generator;
    auto            card        = ProcessCard(aa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["screening"] = "unknown";
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("nonlowercase nuclear key") {
    const BeamCase  aa = SupportedBeams()[2];
    gra::MGraniitti generator;
    auto            card        = ProcessCard(aa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["SCREENING"] = "optical";
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("nonlowercase nuclear selector") {
    const BeamCase  aa = SupportedBeams()[2];
    gra::MGraniitti generator;
    auto            card        = ProcessCard(aa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["emission"] = {"Coherent", "coherent"};
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("incoherent emission from the proton leg") {
    const BeamCase  pa = SupportedBeams()[1];
    gra::MGraniitti generator;
    auto            card        = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["emission"] = {"incoherent", "coherent"};
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }

  SECTION("particle neutron selection does not activate bare hadronic screening") {
    const BeamCase  pa = SupportedBeams()[1];
    gra::MGraniitti generator;
    auto            card             = ProcessCard(pa, "yy[EPA]<F> -> mu+ mu-");
    card["NUCLEAR"]["emd"]           = true;
    card["NUCLEAR"]["fragmentation"] = "internal";
    card["NUCLEAR"]["neutron_class"] = {"*", "n == 0"};
    card["NUCLEAR"]["screening"]      = "optical";
    card["SCATTERING"]["LOOPSCREEN"] = false;
    REQUIRE_NOTHROW(generator.ReadInput(card));
    const auto &upc = generator.proc->state.lts.upc_model;
    REQUIRE(upc != nullptr);
    CHECK_FALSE(upc->HadronicConvolution());
    CHECK(upc->Glauber() == nullptr);
    CHECK(upc->Nodes({}).empty());
    REQUIRE_NOTHROW(generator.proc->PrepareRun());
  }

  SECTION("lepton-ion survival models are inactive") {
    const BeamCase ea = SupportedBeams()[4];
    for (const std::string survival : {"optical", "optical_ggcf", "mc_ggcf"}) {
      DYNAMIC_SECTION(survival) {
        gra::MGraniitti generator;
        auto            card        = ProcessCard(ea, "yy[EPA]<F> -> mu+ mu-");
        card["NUCLEAR"]["screening"] = survival;
        if (survival == "mc_ggcf") { card["NUCLEAR"]["structure"] = "nucleon"; }
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_FALSE(generator.proc->state.lts.upc_model->HadronicConvolution());
      }
    }
  }

  SECTION("EMD requires explicit particle fragmentation") {
    const BeamCase ea                = SupportedBeams()[4];
    auto           card              = ProcessCard(ea, "yy[EPA]<F> -> mu+ mu-");
    card["GENERIC"]["MODELPARAM"]    = FastAdaptationTune();
    card["NUCLEAR"]["emd"]           = true;
    card["NUCLEAR"]["neutron_class"] = {"*", "n > 0"};
    card["NUCLEAR"]["screening"]      = "optical";
    gra::MGraniitti generator;
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
}

// Check beam support through the registered physical process API
TEST_CASE("Registered processes declare physical beam support", "[nuclear][EPA][initial-state][validation]") {
  using gra::nuclear::CollisionType;
  gra::MSubProc registry;
  gra::MPDG     pdg;
  pdg.ReadParticleData(gra::MPDG::DataFile(), "TUNE0");
  const auto processes = registry.CreateAllProcesses();
  for (const auto &process : processes) {
    CAPTURE(process->Label());
    const auto &beams = process->Info().beams;
    CHECK(beams.Supports(CollisionType::PP));
    CHECK_FALSE(beams.Supports(CollisionType::Invalid));
    CHECK_FALSE(beams.Supports(2112, 2212));
    for (const auto &beam : SupportedBeams()) {
      const int first  = pdg.FindByPDGName(beam.beam[0]).pdg;
      const int second = pdg.FindByPDGName(beam.beam[1]).pdg;
      CHECK(beams.Supports(first, second) == beams.Supports(second, first));
    }
  }
  const auto qed = gra::PROC_004_QED_YY_QED().Info().beams;
  CHECK(qed.Supports(2212, 2212));
  CHECK_FALSE(qed.Supports(-2212, 2212));
  const auto epa = gra::PROC_003_QED_YY_EPA().Info().beams;
  for (const auto type : {CollisionType::PP, CollisionType::EE, CollisionType::EP, CollisionType::PA, CollisionType::EA,
                          CollisionType::AA}) {
    CHECK(epa.Supports(type));
  }
  CHECK(epa.Supports(-2212, 2212));

  CHECK(gra::nuclear::BeamChargeX3(11) == -3);
  CHECK(gra::nuclear::BeamChargeX3(-11) == 3);
  CHECK(gra::nuclear::BeamChargeX3(2212) == 3);
  CHECK(gra::nuclear::BeamChargeX3(-2212) == -3);
  CHECK(gra::nuclear::BeamChargeX3(1000822080) == 246);
  CHECK(gra::nuclear::BeamChargeX3(-1000822080) == -246);
  CHECK(gra::nuclear::BeamMassNumber(2212) == 1);
  CHECK(gra::nuclear::BeamMassNumber(1000822080) == 208);
  CHECK(gra::nuclear::BeamMassNumber(-1000822080) == 208);
  CHECK(gra::nuclear::FullBeamEnergy(1000822080, 2510.0) == Approx(522080.0));
  CHECK(gra::nuclear::FullBeamEnergy(-1000822080, 2510.0) == Approx(522080.0));
  CHECK(gra::nuclear::SteeringBeamEnergy(1000822080, 522080.0) == Approx(2510.0));
  CHECK(gra::nuclear::SteeringBeamEnergy(2212, 6500.0) == Approx(6500.0));
  CHECK(gra::nuclear::ClassifyCollision(11, 1000822080) == CollisionType::EA);
  CHECK(gra::nuclear::ClassifyCollision(1000822080, 11) == CollisionType::EA);
  CHECK(gra::nuclear::ClassifyCollision(-11, -1000822080) == CollisionType::EA);
  CHECK(gra::nuclear::ClassifyCollision(1000822080, 2212) == CollisionType::PA);
  CHECK(gra::nuclear::ClassifyCollision(-2212, -1000822080) == CollisionType::PA);
  CHECK(gra::nuclear::ClassifyCollision(-11, -2212) == CollisionType::EP);
  CHECK(gra::nuclear::ClassifyCollision(-2212, -11) == CollisionType::EP);
  CHECK(gra::nuclear::ClassifyCollision(-1000822080, 1000822080) == CollisionType::AA);

  CHECK(gra::nuclear::ClassifyCollision(211, 2212) == CollisionType::Invalid);
  CHECK(gra::nuclear::ClassifyCollision(2112, 2212) == CollisionType::Invalid);
  CHECK(gra::nuclear::ClassifyCollision(1000000010, 2212) == CollisionType::Invalid);
  CHECK_THROWS_AS(gra::nuclear::BeamChargeX3(2112), std::invalid_argument);
  CHECK_THROWS_AS(gra::nuclear::ValidateBeam(-11, -3), std::invalid_argument);
  CHECK_NOTHROW(gra::nuclear::ValidateBeam(-11, 3));
  CHECK_NOTHROW(gra::nuclear::ValidateBeam(-2212, -3));
  CHECK_NOTHROW(gra::nuclear::ValidateBeam(-1000822080, -246));
}

// Check that ISR sampling dimensions follow the charged lepton beam content
TEST_CASE("QED ISR adds one convolution coordinate per lepton beam", "[radiative][ISR][initial-state]") {
  {
    gra::MGraniitti generator;
    auto            card         = ProcessCard(SupportedBeams()[0], "yy[EPA]<F> -> mu+ mu-");
    card["RADIATIVE"]["ISR_QED"] = "YFS";
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
  for (const auto &[beam_index, expected] : std::array<std::pair<std::size_t, unsigned int>, 2>{
           std::pair<std::size_t, unsigned int>{3, 1}, std::pair<std::size_t, unsigned int>{5, 2}}) {
    gra::MGraniitti generator;
    auto            card         = ProcessCard(SupportedBeams()[beam_index], "yy[EPA]<F> -> mu+ mu-");
    card["RADIATIVE"]["ISR_QED"] = "YFS";
    REQUIRE_NOTHROW(generator.ReadInput(card));
    CHECK(generator.proc->GetdLIPSDim() == generator.proc->ProcPtr.LIPSDIM + expected);
  }
}

// Check that screening rejects beams without two hadronic overlap profiles
TEST_CASE("Pomeron screening rejects charged-lepton beams", "[radiative][screening][initial-state]") {
  gra::MGraniitti generator;
  auto            card             = ProcessCard(SupportedBeams()[3], "yy[EPA]<F> -> mu+ mu-");
  card["SCATTERING"]["LOOPSCREEN"] = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
}

// Compare bare GP event weights with the complete configuration projection
TEST_CASE("Bare nuclear GP uses the same current projection as screening",
          "[nuclear][EPA][photoproduction][GP][config]") {
  gra::MGraniitti generator;
  auto            card                         = ProcessCard(SupportedBeams()[2], "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  card["FIDCUTS"]["active"]                    = false;
  card["GENCUTS"]["<F>"]["M"]                  = {3.0, 3.2};
  card["NUCLEAR"]["structure"]                 = "nucleon";
  card["NUCLEAR"]["emission"]                  = {"inclusive", "inclusive"};
  card["NUCLEAR"]["photoproduction"]["target"] = {"inclusive", "inclusive"};
  card["NUCLEAR"]["screening"]                  = "optical";
  card["SCATTERING"]["LOOPSCREEN"]             = false;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
  auto *master = dynamic_cast<gra::MFactorized *>(generator.proc);
  REQUIRE(master != nullptr);
  PhotoNodeProbe process(*master);
  const double   actual = process.ScreenedAmp2(false);
  REQUIRE(actual > 0.0);
  REQUIRE(process.Trace().size() == 1);
  const auto &lts = process.state.lts;
  REQUIRE(lts.upc_event != nullptr);
  REQUIRE(lts.upc_event->HasSamples());
  const auto &trace = process.Trace().front();
  REQUIRE_FALSE(trace.photo_terms.empty());
  gra::nuclear::ScreenPoint born;
  born.amplitude   = trace.amplitude;
  born.photo_terms = trace.photo_terms;
  born.transfer    = {gra::nuclear::RestTransfer(lts.pbeam1, trace.q[0], lts.beam1.mass),
                      gra::nuclear::RestTransfer(lts.pbeam2, trace.q[1], lts.beam2.mass)};
  gra::nuclear::ScreenLayout layout;
  layout.type  = gra::nuclear::ScreenType::Photo;
  layout.photo = gra::nuclear::PhotoChannels(*lts.upc_event);
  const gra::nuclear::MUPCScreen projection(*lts.upc_event, layout, born);
  const auto expected = process.ProcPtr.NormalizeAmplitude(projection.Result().helicity_norm, lts.hamp.metadata);
  REQUIRE(actual == Approx(expected.total_weight).epsilon(2.0e-12));
  const auto   final = gra::nuclear::FinalSectorWeights(layout, expected.helicity_weight);
  const double total = std::accumulate(final.cbegin(), final.cend(), 0.0);
  REQUIRE(total > 0.0);
  for (const auto &i : indices(final)) {
    REQUIRE(process.state.upc_final_weight[i] == Approx(final[i] / total).epsilon(2.0e-12));
  }
}

// Exercise retained optical target terms through the real scalar GP convolution
TEST_CASE("Smooth nuclear GP screens resolved target terms", "[nuclear][EPA][photoproduction][GP][screening]") {
  gra::MGraniitti generator;
  auto            card                         = ProcessCard(SupportedBeams()[1], "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  card["FIDCUTS"]["active"]                    = false;
  card["GENCUTS"]["<F>"]["M"]                  = {3.0, 3.2};
  card["NUCLEAR"]["structure"]                 = "smooth";
  card["NUCLEAR"]["photoproduction"]["target"] = {"coherent", "inclusive"};
  card["NUCLEAR"]["screening"]                  = "optical";
  card["SCATTERING"]["LOOPSCREEN"]             = true;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
  auto *master = dynamic_cast<gra::MFactorized *>(generator.proc);
  REQUIRE(master != nullptr);
  PhotoNodeProbe process(*master);
  REQUIRE_FALSE(process.state.lts.upc_event->HasSamples());
  REQUIRE(process.state.lts.upc_event->HadronicConvolution());
  REQUIRE(process.ScreenedAmp2() > 0.0);
  REQUIRE(process.Trace().size() > 1);
  for (const auto &node : process.Trace()) { REQUIRE_FALSE(node.photo_terms.empty()); }
}

// Exercise mass proposals, completion and cached HepMC3 output in both exact phase spaces
TEST_CASE("Nuclear backend proposals enter production before forward kinematics", "[nuclear][final][event]") {
  for (const std::string mode : {"F", "C"}) {
    DYNAMIC_SECTION(mode) {
      gra::MGraniitti generator;
      auto            card             = ProcessCard(SupportedBeams()[2], "yy[EPA]<" + mode + "> -> mu+ mu-");
      card["SCATTERING"]["LOOPSCREEN"] = false;
      card["NUCLEAR"]["emd"]           = true;
      card["NUCLEAR"]["fragmentation"] = "internal";
      card["NUCLEAR"]["neutron_class"] = {"n == 0", "n == 0"};
      card["FIDCUTS"]["active"]        = false;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      generator.proc->SetNuclearFinalState(std::make_unique<RadiativeIon>());
      CHECK(generator.proc->state.lts.upc_model->Breakup(1)->Enabled());
      CHECK(generator.proc->state.upc_param.reaction == gra::nuclear::ReactionType::External);
      CHECK(generator.proc->state.lts.upc_model->Param().reaction == gra::nuclear::ReactionType::External);
      generator.proc->PrepareRun();
      REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
      const auto &state = generator.proc->state;
      REQUIRE(state.nuclear_event.has_value());
      REQUIRE(state.lts.forward_mass2[0] == Approx(gra::math::pow2(state.lts.beam1.mass + 0.02)));
      REQUIRE(state.lts.forward_mass2[1] == Approx(gra::math::pow2(state.lts.beam2.mass + 0.02)));
      HepMC3::GenEvent event;
      REQUIRE(generator.proc->EventRecord(event));
      REQUIRE_FALSE(state.nuclear_event.has_value());
      for (const auto &parent : gra::nuclear::ForwardParents(event)) {
        REQUIRE(parent->status() == 2);
        REQUIRE(parent->end_vertex()->particles_out().size() == 2);
        REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, state.lts.PDG));
      }
    }
  }
}

// Test AME thresholds and competing particle channels through the actual decay API
TEST_CASE("Internal nuclear decay conserves charge and recoil across thresholds", "[nuclear][final][evaporation]") {
  gra::MGraniitti generator;
  generator.ReadInput(ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-"));
  const auto                &param = generator.proc->state.upc_param;
  gra::nuclear::MEvaporation decay(param.mass, param.geometry.radius_scale, param.emd.gdr,
                                   param.decay_nodes, param.decay_steps);
  // Beam masses and decay thresholds must use the same bare-nucleus convention
  for (const std::string name : {"O16", "Au197", "Pb208"}) {
    const auto particle = generator.proc->state.lts.PDG.FindByPDGName(name);
    const auto isotope = gra::nuclear::DecodeNuclearPDG(particle.pdg);
    CHECK(particle.mass == Approx(decay.Mass(isotope.a, isotope.z)).margin(1.0e-10));
  }
  const double               ground = decay.Mass(16, 8);
  const double               sn     = decay.Mass(15, 8) + decay.Mass(1, 0) - ground;
  const double               sp     = decay.Mass(15, 7) + decay.Mass(1, 1) - ground;
  CHECK(sn == Approx(0.015664).margin(0.000002));
  CHECK(sp == Approx(0.012127).margin(0.000002));
  gra::MRandom random;
  random.SetSeed(84031);
  for (const double excitation : {0.001, 0.030}) {
    std::map<int, unsigned int> counts;
    for (unsigned int trial = 0; trial < 256; ++trial) {
      HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
      const double     mass = ground + excitation, pz = (trial % 2 ? -1.0 : 1.0) * 2700.0 * mass;
      auto             parent =
          std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, pz, std::hypot(mass, pz)), 1000080160, 3);
      parent->set_generated_mass(mass);
      event.add_particle(parent);
      REQUIRE_NOTHROW(decay.Decay(event, parent, random));
      REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
      for (const auto &particle : event.particles()) {
        if (particle->status() == 1) { ++counts[particle->pid()]; }
      }
    }
    CHECK(counts[22] > 0);
    if (excitation < std::min(sn, sp)) {
      CHECK(counts[2112] == 0);
      CHECK(counts[2212] == 0);
    } else {
      CHECK(counts[2112] > 0);
      CHECK(counts[2212] > 0);
    }
  }
  // Resolve the measured Be-8 ground state even though alpha emission is below the Coulomb barrier
  const double be8 = decay.Mass(8, 4), alpha = decay.Mass(4, 2);
  CHECK(be8 - 2.0 * alpha == Approx(0.00009184).margin(0.000001));
  for (const auto &ion : std::array<std::array<int, 2>, 3>{{{5, 2}, {5, 3}, {8, 4}}}) {
    for (const double excitation : {0.0, 0.001}) {
      for (const double boost : {0.0, -2700.0, 2700.0}) {
        CAPTURE(ion, excitation, boost);
        HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
        const double     mass = decay.Mass(ion[0], ion[1]) + excitation, pz = boost * mass;
        auto parent = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, pz, std::hypot(mass, pz)),
                                                            1000000000 + ion[1] * 10000 + ion[0] * 10, 1);
        parent->set_generated_mass(mass);
        event.add_particle(parent);
        REQUIRE_NOTHROW(decay.Decay(event, parent, random));
        REQUIRE(parent->end_vertex());
        REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
        unsigned int alphas = 0;
        for (const auto &particle : event.particles()) {
          if (particle->status() != 1) { continue; }
          CHECK((particle->pid() == 1000020040 || particle->pid() == 22 ||
                 particle->pid() == (ion[1] == 2 ? 2112 : 2212)));
          if (particle->pid() == 1000020040) {
            ++alphas;
            CHECK(particle->generated_mass() == Approx(alpha).epsilon(1.0e-14));
          }
        }
        CHECK(alphas == (ion[0] == 8 ? 2 : 1));
      }
    }
  }
}

// Construct a central system with fixed nuclear t and the proposed forward masses
HepMC3::GenEvent ImpulseEvent(const HepMC3::GenEvent &beams, const std::array<double, 2> &mass, double t) {
  HepMC3::GenEvent   event(beams.momentum_unit(), beams.length_unit());
  auto               central = std::make_shared<HepMC3::GenVertex>();
  HepMC3::FourVector total;
  for (const auto &i : indices(mass)) {
    const auto   beam = std::make_shared<HepMC3::GenParticle>(*beams.beams().at(i));
    const double ma = beam->generated_mass(), omega = ((mass[i] - ma) * (mass[i] + ma) - t) / (2.0 * ma);
    gra::M4Vec   recoil(0, 0, (i == 0 ? -1.0 : 1.0) * std::sqrt(omega * omega - t), ma + omega);
    const auto  &p = beam->momentum();
    gra::kinematics::LorentzBoost(gra::M4Vec(p.px(), p.py(), p.pz(), p.e()), ma, recoil, 1);
    const auto q       = p - gra::aux::M4Vec2HepMC3(recoil);
    auto       forward = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(recoil), beam->pid(), 3);
    forward->set_generated_mass(mass[i]);
    auto photon = std::make_shared<HepMC3::GenParticle>(q, 22, 3);
    auto vertex = std::make_shared<HepMC3::GenVertex>();
    vertex->add_particle_in(beam);
    vertex->add_particle_out(forward);
    vertex->add_particle_out(photon);
    event.add_vertex(vertex);
    forward->add_attribute("graniitti_upc_leg", std::make_shared<HepMC3::IntAttribute>(i + 1));
    event.add_attribute("graniitti_upc_final_sector_" + std::to_string(i + 1),
                        std::make_shared<HepMC3::StringAttribute>(i == 0 ? "incoherent" : "coherent"));
    central->add_particle_in(photon);
    total = total + q;
  }
  central->add_particle_out(std::make_shared<HepMC3::GenParticle>(total, 23, 1));
  event.add_vertex(central);
  return event;
}

// Check inclusive normalization and a resolved charge knockout at fixed hard transfer
TEST_CASE("Internal impulse response normalizes the joint recoil distribution", "[nuclear][final][impulse]") {
  gra::MGraniitti generator;
  auto            card         = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emission"]  = {"incoherent", "coherent"};
  card["NUCLEAR"]["structure"] = "nucleon";
  generator.ReadInput(card);
  const auto      &state = generator.proc->state;
  HepMC3::GenEvent beams(HepMC3::Units::GEV, HepMC3::Units::MM);
  auto first  = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam1), state.lts.beam1.pdg, 4);
  auto second = std::make_shared<HepMC3::GenParticle>(gra::aux::M4Vec2HepMC3(state.lts.pbeam2), state.lts.beam2.pdg, 4);
  first->set_generated_mass(state.lts.beam1.mass);
  second->set_generated_mass(state.lts.beam2.mass);
  beams.set_beam_particles(first, second);
  gra::nuclear::MFinal model(std::make_unique<gra::nuclear::MIncoherent>(state.lts.upc_model));
  const auto          &nucleus = *state.lts.upc_model->Nucleus(1);
  const auto          &param   = state.upc_param;
  const double         radius  = std::sqrt(5.0 / 3.0) * nucleus.MatterDensity().Rms() / std::cbrt(nucleus.A());
  const gra::nuclear::MEvaporation decay(param.mass, radius, param.emd.gdr, param.decay_nodes,
                                         param.decay_steps);
  const double                     ef = std::hypot(decay.Mass(1, 1), decay.Fermi(nucleus.A(), nucleus.Z(), true));
  gra::MRandom                     random;
  random.SetSeed(710204);
  for (const double t : {-0.04, -0.36}) {
    double                 sum = 0.0, sum2 = 0.0;
    double                 energy = 0.0, rapidity = 0.0, correlation = 0.0;
    constexpr unsigned int trials = 4096;
    for (unsigned int trial = 0; trial < trials; ++trial) {
      const auto   mass    = model.SampleMasses(beams, random);
      auto         event   = ImpulseEvent(beams, mass.mass, t);
      const auto   central = event.particles().back();
      const double weight  = model.Complete(event, random);
      sum += weight;
      sum2 += weight * weight;
      if (!(weight > 0.0)) { continue; }
      const auto roots = gra::nuclear::ForwardParents(event);
      REQUIRE(roots[0]->end_vertex() != nullptr);
      const auto proton = roots[0]->end_vertex()->particles_out().front();
      CHECK(proton->pid() == 2212);
      CHECK(gra::aux::HepMC2M4Vec(proton->momentum()).DotM(state.lts.pbeam1) / state.lts.beam1.mass >= ef - 1e-7);
      CHECK(roots[1]->end_vertex() == nullptr);
      const double y = central->momentum().rap(), e = proton->momentum().e();
      energy += weight * e;
      rapidity += weight * y;
      correlation += weight * e * y;
    }
    const double mean = sum / trials, error = std::sqrt(std::max(0.0, sum2 / trials - mean * mean) / trials);
    CHECK(mean == Approx(1.0).margin(6.0 * error + 0.002));
    // At fixed spacelike t the leading proton energy follows the central rapidity
    REQUIRE(sum > 0.0);
    CHECK(correlation / sum - energy * rapidity / (sum * sum) > 0.0);
  }
}

// Reject a photonuclear fit for another isotope before internal EMD sampling
TEST_CASE("Internal EMD requires a photonuclear fit for the beam isotope", "[nuclear][breakup][steering][isotope]") {
  gra::MGraniitti generator;
  auto            card = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  REQUIRE_NOTHROW(generator.ReadInput(card));

  auto param = generator.proc->state.upc_param;
  ++param.emd.isotope.front().a;
  REQUIRE_NOTHROW(generator.proc->SetNuclear(param));
  param.additional_emd = true;
  param.reaction       = gra::nuclear::ReactionType::Internal;
  REQUIRE_THROWS_WITH(generator.proc->SetNuclear(param),
                      Catch::Matchers::Contains("no EMD response for the target isotope"));
}

// Exercise the internal backend in the same production path as external reactions
TEST_CASE("Internal nuclear production retains completed forward particles", "[nuclear][final][event][impulse]") {
  gra::MGraniitti generator;
  auto            card                         = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
  card["NUCLEAR"]["emission"]                  = {"coherent", "coherent"};
  card["NUCLEAR"]["photoproduction"]["target"] = {"incoherent", "incoherent"};
  card["NUCLEAR"]["structure"]                 = "nucleon";
  card["FIDCUTS"]["active"]                    = false;
  card["GENCUTS"]["<F>"]["Pt"]                 = {0.2, 0.5};
  generator.ReadInput(card);
  auto param           = generator.proc->state.upc_param;
  param.additional_emd = true;
  param.reaction       = gra::nuclear::ReactionType::Internal;
  generator.proc->SetNuclear(param);
  generator.proc->PrepareRun();
  HepMC3::GenEvent event;
  bool             accepted = false;
  for (std::size_t trial = 0; trial < 4096 && !accepted; ++trial) {
    gra::MEventWeightState state;
    const double weight = generator.proc->EventWeight(EPAEventPoint(generator.proc->ProcPtr.LIPSDIM, trial), state);
    if (!state.Valid() || !(weight > 0.0)) { continue; }
    HepMC3::GenEvent candidate;
    REQUIRE(generator.proc->EventRecord(candidate));
    for (const auto &parent : gra::nuclear::ForwardParents(candidate)) {
      if (parent->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() == "coherent" &&
          parent->end_vertex() != nullptr) {
        accepted = true;
      }
    }
    if (accepted) { event = candidate; }
  }
  REQUIRE(accepted);
  const auto roots = gra::nuclear::ForwardParents(event);
  REQUIRE(event.heavy_ion() != nullptr);
  REQUIRE(event.heavy_ion()->impact_parameter < 0.0);
  unsigned int knockout = 0;
  for (const auto &parent : roots) {
    if (parent->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() == "incoherent") {
      REQUIRE(parent->end_vertex() != nullptr);
      const int pid = parent->end_vertex()->particles_out().front()->pid();
      REQUIRE((pid == 2212 || pid == 2112));
      ++knockout;
    }
    REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
  }
  CHECK(knockout == 1);
}


// Reject hard sources whose impact dependence is not supplied by coherent EPA
TEST_CASE("Internal EMD rejects unsupported hard source matching", "[nuclear][final][EMD][steering]") {
  for (const bool fusion : {false, true}) {
    gra::MGraniitti generator;
    auto card = ProcessCard(SupportedBeams()[2], fusion ? "yy[EPA]<F> -> mu+ mu-" : "ygg[jpsi]<F> -> mu+ mu-");
    card["NUCLEAR"]["emission"]  = {"incoherent", "coherent"};
    card["NUCLEAR"]["structure"] = "nucleon";
    generator.ReadInput(card);
    auto param           = generator.proc->state.upc_param;
    param.additional_emd = true;
    param.reaction       = gra::nuclear::ReactionType::Internal;
    REQUIRE_THROWS_AS(generator.proc->SetNuclear(param), std::invalid_argument);
  }
}

// Initialize EMD with screened and unscreened configuration survival
TEST_CASE("Internal EMD supports configuration survival",
          "[nuclear][final][EMD][mc_ggcf]") {
  auto beams        = std::vector<BeamCase>{SupportedBeams()[1], SupportedBeams()[2], SupportedBeams()[1]};
  beams.back().name = "Ap";
  std::swap(beams.back().beam[0], beams.back().beam[1]);
  std::swap(beams.back().energy[0], beams.back().energy[1]);
  std::swap(beams.back().type[0], beams.back().type[1]);
  for (const auto &beam : beams) {
    for (const bool screening : {false, true}) {
      DYNAMIC_SECTION(beam.name << " LOOPSCREEN=" << screening) {
        auto card                        = ProcessCard(beam, "yy[EPA]<F> -> mu+ mu-");
        card["SCATTERING"]["LOOPSCREEN"] = screening;
        card["NUCLEAR"]["structure"]     = "nucleon";
        card["NUCLEAR"]["screening"]      = "mc_ggcf";
        card["NUCLEAR"]["fragmentation"] = "internal";
        card["NUCLEAR"]["emd"]           = true;
        gra::MGraniitti generator;
        REQUIRE_NOTHROW(generator.ReadInput(card));
        REQUIRE_NOTHROW(generator.proc->PrepareRun());
        CHECK(generator.proc->state.lts.upc_model->HadronicConvolution() == screening);
        REQUIRE(RequirePositiveEPAEvent(*generator.proc, false) > 0.0);
        CHECK(generator.proc->state.lts.upc_event->Convolution());
        CHECK_FALSE(generator.proc->state.lts.upc_event->HadronicConvolution());
      }
    }
  }
}


// Reject unsupported scales and small-x support before either photon direction is sampled
TEST_CASE("LTA validates the production domain before sampling", "[nuclear][EPA][LTA][validation]") {
  const auto base = ProcessCard(SupportedBeams()[2], "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  for (const bool hard : {false, true}) {
    auto card                                          = base;
    card["NUCLEAR"]["photoproduction"]["target_model"] = "LTA";
    card["GENCUTS"]["<F>"]["M"]   = hard ? std::array<double, 2>{3.0, 20.0} : std::array<double, 2>{3.0, 3.2};
    card["GENCUTS"]["<F>"]["Rap"] = hard ? std::array<double, 2>{-0.8, 0.8} : std::array<double, 2>{-9.0, 9.0};
    gra::MGraniitti generator;
    REQUIRE_THROWS_AS(generator.ReadInput(card), std::invalid_argument);
  }
}




// Partition the actual decay outcomes once without analytic neutron weights
TEST_CASE("Coherent EMD particle neutron classes close in production", "[nuclear][final][EMD][event][steering]") {
  const std::string process = GENERATE("yy[EPA]<F> -> mu+ mu-", "ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> pi+ pi- @RES{rho_770:1}",
                                       "TP[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  INFO(process);
  std::array<gra::MGraniitti, 8>                  generator;
  const std::array<std::array<std::string, 2>, 8> tags = {
      {{"*", "*"}, {"n == 0", "n == 0"}, {"n == 0", "n > 0"}, {"n > 0", "n == 0"}, {"n > 0", "n > 0"}, {"n == 0", "*"}, {"n == 1", "*"}, {"n >= 2", "*"}}};
  for (const auto &i : indices(generator)) {
    auto card                        = ProcessCard(SupportedBeams()[2], process);
    card["NUCLEAR"]["emd"]           = true;
    card["NUCLEAR"]["fragmentation"] = "internal";
    card["NUCLEAR"]["neutron_class"] = tags[i];
    card["SCATTERING"]["LOOPSCREEN"] = true;
    card["FIDCUTS"]["active"]        = false;
    card["GENCUTS"]["<F>"]["M"] = process.starts_with("ygg") ? std::array<double, 2>{3.0, 3.2}
                                                               : std::array<double, 2>{0.6, 0.9};
    card["GENCUTS"]["<F>"]["Pt"]     = {0.001, 0.06};
    generator[i].ReadInput(card);
    // Use a common proposal for the exact pointwise partition of the neutron classes
    auto param = generator[i].proc->state.upc_param;
    param.emd_focus = 0.0;
    param.emd_condition = false;
    generator[i].proc->SetNuclear(param);
    generator[i].proc->PrepareRun();
  }
  std::array<double, 8> sum{};
  bool                  breakup = false, current_checked = false;
  unsigned int          kinematic = 0, technical = 0, amplitude = 0;
  for (unsigned int trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    std::array<double, 8> weight{};
    for (const auto &i : indices(generator)) {
      auto &process = *generator[i].proc;
      process.state.random.SetSeed(78501 + trial);
      gra::MEventWeightState state;
      weight[i] = process.EventWeight(EPAEventPoint(process.ProcPtr.LIPSDIM, trial), state);
      if (i == 0) {
        kinematic += !state.kinematics_ok;
        technical += state.technical_failure;
        amplitude += !state.amplitude_ok;
      }
      if (!(weight[i] > 0.0)) { continue; }
      REQUIRE(state.Valid());
      if (i == 0 && !current_checked) {
        // The coherent current remains elastic while the net recoil carries EMD excitation
        auto             boosted = process.state.lts;
        const gra::M4Vec boost(0.2, -0.1, 0.3, 1.5);
        for (auto *momentum :
             {&boosted.pbeam1, &boosted.pbeam2, &boosted.q1, &boosted.q2, &boosted.pfinal[1], &boosted.pfinal[2]}) {
          gra::kinematics::LorentzBoost(boost, boost.M(), *momentum, 1);
        }
        for (const auto leg : {gra::ForwardBeamLeg::Upper, gra::ForwardBeamLeg::Lower}) {
          const auto   hard  = gra::ResolveForwardLegState(process.state.lts, leg);
          const double mass2 = gra::math::pow2(hard.emitter.mass);
          CHECK(hard.outgoing.M2() == Approx(mass2).epsilon(1e-8));
          CHECK(hard.t == Approx(-(hard.qt * hard.qt + hard.xi * hard.xi * mass2) / (1.0 - hard.xi)).margin(1e-8));
          auto transfer = hard.transfer;
          gra::kinematics::LorentzBoost(boost, boost.M(), transfer, 1);
          const auto covariant = gra::ResolveForwardLegState(boosted, leg);
          for (unsigned int component = 0; component < 4; ++component) {
            CHECK(covariant.transfer[component] == Approx(transfer[component]).margin(1e-6));
          }
        }
        current_checked = true;
      }
      HepMC3::GenEvent event;
      REQUIRE(process.EventRecord(event));
      REQUIRE(event.heavy_ion() != nullptr);
      CHECK(event.heavy_ion()->impact_parameter < 0.0);
      const auto roots = gra::nuclear::ForwardParents(event);
      for (const auto &leg : indices(roots)) {
        CHECK(roots[leg]->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() == "coherent");
        const auto branch  = gra::nuclear::ValidateNuclearDecay(event, roots[leg], process.state.lts.PDG);
        std::size_t neutron = 0;
        for (const auto id : branch) {
          const auto &particle = event.particles()[id - 1];
          neutron += particle->status() == 1 && particle->pid() == 2112;
        }
        CHECK(gra::nuclear::ParseNeutronSelection(tags[i][leg]).Accept(neutron));
        breakup |= neutron > 0;
      }
      sum[i] += weight[i];
    }
    CHECK(weight[0] == Approx(weight[5] + weight[6] + weight[7]).epsilon(1e-10));
    CHECK(weight[0] == Approx(weight[1] + weight[2] + weight[3] + weight[4]).epsilon(1e-10));
  }
  CAPTURE(kinematic, technical, amplitude);
  CHECK(sum[0] > 0.0);
  CHECK(breakup);
  CHECK(current_checked);
}

// Exercise the projected coherent currents through the complete optical screening loop
TEST_CASE("Coherent EMD retains screened particle support", "[nuclear][final][EMD][screening]") {
  auto card                        = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"]           = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["SCATTERING"]["LOOPSCREEN"] = true;
  card["FIDCUTS"]["active"]        = false;
  card["GENCUTS"]["<F>"]["M"]      = {3.0, 3.2};
  card["GENCUTS"]["<F>"]["Pt"]     = {0.001, 0.06};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  generator.proc->PrepareRun();
  REQUIRE(generator.proc->state.lts.upc_model->HadronicConvolution());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc, true) > 0.0);
  HepMC3::GenEvent event;
  REQUIRE(generator.proc->EventRecord(event));
  CHECK(event.heavy_ion()->impact_parameter < 0.0);
  for (const auto &root : gra::nuclear::ForwardParents(event)) {
    CHECK_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, root, generator.proc->state.lts.PDG));
  }
}

// Switching EMD off keeps coherent ions intact even when a decay backend is selected
TEST_CASE("Fragmentation does not implicitly enable EMD", "[nuclear][final][EMD][event][steering]") {
  gra::MGraniitti generator;
  auto            card             = ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-");
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["neutron_class"] = {"n == 0", "n == 0"};
  generator.ReadInput(card);
  generator.proc->PrepareRun();
  REQUIRE_FALSE(generator.proc->state.lts.upc_model->Breakup(1)->Enabled());
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
  HepMC3::GenEvent event;
  REQUIRE(generator.proc->EventRecord(event));
  for (const auto &root : gra::nuclear::ForwardParents(event)) {
    CHECK(root->end_vertex() == nullptr);
    CHECK(root->status() == 1);
  }
  card["NUCLEAR"]["neutron_class"] = {"n > 0", "*"};
  gra::MGraniitti rejected;
  REQUIRE_THROWS_AS(rejected.ReadInput(card), std::invalid_argument);
}

// Check neutron thresholds with AME masses and the actual particle decay API
TEST_CASE("Particle evaporation respects the AME neutron threshold", "[nuclear][final][EMD][response]") {
  gra::MGraniitti generator;
  auto            card             = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"]           = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["neutron_class"] = {"n == 0", "n == 0"};

  generator.ReadInput(card);
  const auto  &param   = generator.proc->state.upc_param;
  const auto  &upc     = *generator.proc->state.lts.upc_model;
  const auto  &nucleus = *upc.Nucleus(1);
  const double radius  = std::sqrt(5.0 / 3.0) * nucleus.MatterDensity().Rms() / std::cbrt(nucleus.A());
  const gra::nuclear::MEvaporation decay(param.mass, radius, param.emd.gdr, param.decay_nodes,
                                         param.decay_steps);
  const double                     ground    = decay.Mass(nucleus.A(), nucleus.Z());
  const double                     threshold = decay.Mass(nucleus.A() - 1, nucleus.Z()) + decay.Mass(1, 0) - ground;
  gra::MRandom                     random;
  random.SetSeed(730124);
  for (const double energy : {0.5 * threshold, 1.5 * threshold, 3.0 * threshold}) {
    unsigned int zero = 0;
    for (unsigned int trial = 0; trial < 256; ++trial) {
      HepMC3::GenEvent event;
      auto parent = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 0, ground + energy), 1000822080, 3);
      parent->set_generated_mass(ground + energy);
      event.add_particle(parent);
      decay.Decay(event, parent, random);
      REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
      zero += std::none_of(event.particles().begin(), event.particles().end(),
                           [](const auto &p) { return p->status() == 1 && p->pid() == 2112; });
    }
    const double particle = static_cast<double>(zero) / 256;
    if (energy < threshold) {
      CHECK(particle == Approx(1.0));
    } else {
      CHECK(particle < 0.5);
    }
  }
}

// Validate shared decay inputs without EMD and retain finite high excitation cascades
TEST_CASE("Nuclear evaporation remains finite beyond exponential overflow",
          "[nuclear][final][evaporation][validation]") {
  gra::MGraniitti generator;
  auto            card             = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["emd"]           = false;
  generator.ReadInput(card);
  auto         param   = generator.proc->state.upc_param;
  const auto  &nucleus = *generator.proc->state.lts.upc_model->Nucleus(1);
  const double radius  = std::sqrt(5.0 / 3.0) * nucleus.MatterDensity().Rms() / std::cbrt(nucleus.A());
  const gra::nuclear::MEvaporation decay(param.mass, radius, param.emd.gdr, param.decay_nodes,
                                         param.decay_steps);
  gra::MRandom                     random;
  random.SetSeed(93173);
  for (const double excitation : {0.02, 5.0, 10.0, 20.0}) {
    HepMC3::GenEvent event;
    const double     mass   = decay.Mass(nucleus.A(), nucleus.Z()) + excitation;
    auto             parent = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, 0, mass), 1000822080, 1);
    parent->set_generated_mass(mass);
    event.add_particle(parent);
    REQUIRE_NOTHROW(decay.Decay(event, parent, random));
    CHECK_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
  }
  const auto original = param.emd.gdr;
  for (const unsigned int fault : {0U, 1U, 2U}) {
    param.emd.gdr = original;
    if (fault == 0) { param.emd.gdr.systematics.energy_b = -1.0; }
    if (fault == 1) { param.emd.gdr.isotope.front().width = -1.0; }
    if (fault == 2) { param.emd.gdr.isotope.push_back(param.emd.gdr.isotope.front()); }
    CHECK_THROWS_AS(generator.proc->SetNuclear(param), std::invalid_argument);
  }
}

// Compare identical physical recoils across coherent and inclusive target steering
TEST_CASE("Hard recoil projection removes only independent EMD excitation", "[nuclear][final][EMD][covariance]") {
  gra::MGraniitti generator;
  generator.ReadInput(ProcessCard(SupportedBeams()[2], "ygg[jpsi]<F> -> mu+ mu-"));
  auto param = generator.proc->state.upc_param;
  for (const double knockout : {0.0, 0.1}) {
    std::array<double, 2> transfer{};
    for (const auto &i : indices(transfer)) {
      param.target[0] = i == 0 ? gra::nuclear::CoherenceType::Coherent : gra::nuclear::CoherenceType::Inclusive;
      generator.proc->SetNuclear(param);
      auto         lts = generator.proc->state.lts;
      const double pt = 0.02, xi = 1.0e-5, emd = 0.02;
      const double mass = lts.beam1.mass + knockout + emd;
      const double plus = (1.0 - xi) * lts.pbeam1.LightconePos(), minus = (mass * mass + pt * pt) / plus;
      lts.pfinal.resize(3);
      lts.pfinal[1]      = gra::M4Vec(pt, 0, (plus - minus) / 2.0, (plus + minus) / 2.0);
      lts.q1             = lts.pbeam1 - lts.pfinal[1];
      lts.t1             = lts.q1.M2();
      lts.qt1            = pt;
      lts.forward_emd[0] = emd;
      const auto hard    = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
      transfer[i]        = hard.t;
      CHECK(hard.mass2 == Approx(gra::math::pow2(lts.beam1.mass + knockout)).epsilon(1e-8));
      lts.forward_emd[0] = 0.0;
      const auto net     = gra::ResolveForwardLegState(lts, gra::ForwardBeamLeg::Upper);
      CHECK(net.mass2 == Approx(mass * mass).epsilon(1e-8));
    }
    CHECK(transfer[0] == Approx(transfer[1]).margin(1e-12));
  }
}

// Retain the complex interference between diagonal VMD and neutral isovector conversion
TEST_CASE("Nuclear tensor VMD amplitudes preserve isospin interference", "[nuclear][tensor][photoproduction][isospin]") {
  const auto beam = SupportedBeams()[4];
  auto card = ProcessCard(beam, "TP[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  card["FIDCUTS"]["active"] = false;
  card["GENCUTS"]["<F>"]["M"] = {0.6, 0.9};
  card["NUCLEAR"]["photoproduction"]["target_model"] = "impulse";
  gra::MGraniitti generator;
  generator.ReadInput(card);
  generator.proc->PrepareRun();
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
  auto original = generator.proc->state.lts;
  generator.proc->ProcPtr.GetBareAmplitude2(original);
  const auto full = original.hamp;
  std::array<std::vector<std::complex<double>>, 2> part;
  for (const auto component : indices(part)) {
    auto lts = original;
    lts.amplitude = {};
    for (auto& [name, res] : lts.process.RESONANCES) {
      for (auto& channel : res.TP.channels) {
        const int exchange = channel.exchange[channel.exchange[0] == gra::PDG::PDG_gamma ? 1 : 0];
        if (lts.PDG.FindByPDG(exchange).isospinX2 != static_cast<int>(2 * component)) {
          channel.active_g_tensor.clear();
          channel.active_vmd.clear();
        }
      }
    }
    generator.proc->ProcPtr.GetBareAmplitude2(lts);
    part[component].assign(lts.hamp.begin(), lts.hamp.end());
    REQUIRE(gra::SquaredNorm(part[component]) > 0.0);
  }
  REQUIRE(part[0].size() == full.size());
  REQUIRE(part[1].size() == full.size());
  auto difference = part[0];
  gra::AddScaled(difference, part[1], 1.0);
  gra::AddScaled(difference, full, -1.0);
  CHECK(gra::SquaredNorm(difference) < 1.0e-22 * gra::SquaredNorm(full));

  auto reversed = original;
  reversed.amplitude = {};
  for (auto& [name, res] : reversed.process.RESONANCES) {
    for (auto& channel : res.TP.channels) {
      const int exchange = channel.exchange[channel.exchange[0] == gra::PDG::PDG_gamma ? 1 : 0];
      if (reversed.PDG.FindByPDG(exchange).isospinX2 != 2) { continue; }
      for (auto& g : channel.g_tensor) { g = -g; }
      for (auto& mixing : channel.VMD_MIXING) {
        for (auto& g : mixing.coupling) { g = -g; }
      }
    }
  }
  generator.proc->ProcPtr.GetBareAmplitude2(reversed);
  difference = part[0];
  gra::AddScaled(difference, part[1], -1.0);
  gra::AddScaled(difference, reversed.hamp, -1.0);
  CHECK(gra::SquaredNorm(difference) < 1.0e-22 * gra::SquaredNorm(full));
  const double diagonal = gra::SquaredNorm(part[0]) + gra::SquaredNorm(part[1]);
  CHECK(std::abs(gra::SquaredNorm(full) - diagonal) > 1.0e-8 * diagonal);
}

// Reject an unspecified outgoing absorption amplitude before event sampling
TEST_CASE("Nuclear VMD shadowing requires a diagonal amplitude", "[nuclear][tensor][isospin][initialization]") {
  const auto directory = std::filesystem::path(gra::aux::ResolveProjectPath("tmp/tensor_vmd_missing_diag"));
  std::filesystem::create_directories(directory);
  std::filesystem::copy(gra::ResolveModelTuneDir("TUNE0"), directory,
      std::filesystem::copy_options::recursive | std::filesystem::copy_options::overwrite_existing);
  const auto filename = directory / "RES/rho_770.json";
  auto resonance = nlohmann::json::parse(gra::aux::GetInputData(filename.string()));
  for (auto& [key, value] : resonance["PARAM_RES"]["MODELS"]["TP"].items()) {
    if (key.front() == '[') {
      for (auto& coupling : value["g_tensor"]) { coupling = 0.0; }
    }
  }
  std::ofstream(filename) << resonance.dump(2);
  auto card = ProcessCard(SupportedBeams()[4], "TP[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  card["GENERIC"]["MODELPARAM"] = directory.string();
  gra::MGraniitti generator;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_THROWS_WITH(generator.proc->PrepareRun(), Catch::Contains("diagonal vector-nucleon amplitude"));
}

// Preserve single proton excitation and its normalization in either pA beam order
TEST_CASE("pA photoproduction keeps proton dissociation through nuclear completion", "[nuclear][photoproduction][dissociation]") {
  const bool reverse = GENERATE(false, true);
  auto beam = SupportedBeams()[1];
  if (reverse) {
    std::swap(beam.beam[0], beam.beam[1]);
    std::swap(beam.energy[0], beam.energy[1]);
    std::swap(beam.type[0], beam.type[1]);
  }
  auto card = ProcessCard(beam, "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  card["SCATTERING"]["NSTARS"] = 1;
  card["GENCUTS"]["<F>"]["M"] = {3.0, 3.2};
  card["GENCUTS"]["<F>"]["Xi"] = {0.0, 1.0e-8};
  card["GENCUTS"]["<F>"]["Pt"] = {0.0, 3.0};
  card["FIDCUTS"]["active"] = false;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["emd"] = true;
  gra::MGraniitti generator;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  auto& process = *generator.proc;
  bool accepted = false;
  for (std::size_t trial = 0; trial < EPA_EVENT_TRIALS; ++trial) {
    auto point = EPAEventPoint(8, trial);
    point.resize(process.GetdLIPSDim(), 0.25);
    gra::MEventWeightState state;
    state.include_screening = false;
    const double weight = process.EventWeight(point, state);
    if (!state.Valid() || !(weight > 0.0)) { continue; }
    REQUIRE(std::isfinite(weight));
    CHECK(process.state.lts.excite1 == !reverse);
    CHECK(process.state.lts.excite2 == reverse);
    auto lts = process.state.lts;
    const double norm = process.ProcPtr.GetBareAmplitude2(lts);
    const auto amplitude = lts.hamp;
    REQUIRE(norm > 0.0);
    // Check rotational covariance and recover the complex amplitude after the inverse rotation
    for (const double angle : {0.73, -0.73}) {
      lts.amplitude = {};
      for (auto& momentum : lts.pfinal) { momentum.RotateZ(angle); }
      lts.q1.RotateZ(angle);
      lts.q2.RotateZ(angle);
      for (auto& branch : lts.decaytree) { branch.p4.RotateZ(angle); }
      CHECK(process.ProcPtr.GetBareAmplitude2(lts) == Approx(norm).epsilon(2.0e-7));
    }
    auto difference = lts.hamp;
    gra::AddScaled(difference, amplitude, -1.0);
    CHECK(gra::SquaredNorm(difference) < 1.0e-12 * gra::SquaredNorm(amplitude));
    HepMC3::GenEvent event(HepMC3::Units::GEV, HepMC3::Units::MM);
    REQUIRE(process.EventRecord(event));
    const auto roots = gra::nuclear::ForwardParents(event);
    const auto& proton = roots[reverse ? 1 : 0];
    CHECK(std::abs(proton->pid()) == gra::PDG::PDG_NSTAR);
    CHECK(proton->generated_mass() > process.GetInitialState()[reverse ? 1 : 0].mass);
    CHECK(proton->generated_mass() == Approx(process.state.lts.pfinal[reverse ? 2 : 1].M()).epsilon(1.0e-6));
    accepted = true;
    break;
  }
  REQUIRE(accepted);
}

// Keep ordinary photonuclear shadowing independent of GGCF ion survival
TEST_CASE("Photonuclear Glauber shadowing is independent of GGCF survival", "[nuclear][ggcf][steering][photoproduction]") {
  const auto process = GENERATE("GP[RES]<F> -> pi+ pi- @RES{rho_770:1}", "GP[RES]<F> -> mu+ mu- @RES{jpsi:1}");
  auto card = ProcessCard(SupportedBeams()[1], process);
  card["NUCLEAR"]["screening"] = "optical_ggcf";
  card["NUCLEAR"]["photoproduction"]["target_model"] = "glauber";
  std::complex<double> reference = 0.0;
  for (const bool screening : {false, true}) {
    gra::MGraniitti generator;
    card["SCATTERING"]["LOOPSCREEN"] = screening;
    REQUIRE_NOTHROW(generator.ReadInput(card));
    const auto &upc = generator.proc->state.lts.upc_model;
    REQUIRE(upc != nullptr);
    CHECK(upc->HadronicConvolution() == screening);
    const auto profile = gra::test::PhotoProfile(20.0, 0.1, 0.0);
    const int leg = upc->Photo(1) != nullptr ? 1 : 2;
    CHECK(upc->Param().photo[leg - 1].fluctuation.omega == Approx(0.0).margin(1.0e-14));
    const auto amplitude = upc->Photo(leg)->Factors(profile, 0.02, -0.03, 0.01,
                                                   gra::nuclear::PhotonDirection::PositiveZ).coherent;
    if (!screening) {
      reference = amplitude;
      CHECK(upc->Glauber() == nullptr);
    } else {
      CHECK(upc->Param().glauber.profile.omega > 0.0);
      CHECK(upc->Param().survival_eikonal == upc->Param().ggcf.eikonal);
      CHECK(std::abs(amplitude - reference) < 1.0e-10 * std::abs(reference));
    }
  }
  CHECK_THROWS_AS(gra::nuclear::ParsePhotoModel("glauber_ggcf"), std::invalid_argument);
}

// Reject ambiguous isotope tables through the real model-card reader
TEST_CASE("EMD isotope arrays reject duplicate and invalid identities", "[nuclear][breakup][steering][isotope]") {
  const auto tune = gra::ResolveModelTuneDir("TUNE0");
  auto physics = nlohmann::json::parse(gra::aux::GetInputData(tune + "/GENERAL.json")).at("PARAM_NUCLEAR");
  const auto numerics = nlohmann::json::parse(gra::aux::GetInputData(tune + "/NUMERICS.json")).at("NUMERICS_NUCLEAR");
  const auto card = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  REQUIRE_NOTHROW(gra::nuclear::ReadUPCSteering(card.at("NUCLEAR"), physics, numerics, tune));
  auto &rows = physics["EMD"]["isotopes"];
  SECTION("duplicate target") { rows.push_back(rows.front()); }
  SECTION("fractional mass number") { rows.front()["A"] = 16.5; }
  SECTION("invalid charge") { rows.front()["Z"] = rows.front()["A"]; }
  SECTION("object instead of array") { rows = rows.front(); }
  CHECK_THROWS_AS(gra::nuclear::ReadUPCSteering(card.at("NUCLEAR"), physics, numerics, tune), std::invalid_argument);
}

// Check that a response is selected by target identity rather than table position
TEST_CASE("EMD target selection is independent of isotope ordering", "[nuclear][breakup][steering][isotope]") {
  gra::MGraniitti generator;
  auto card = ProcessCard(SupportedBeams()[1], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  REQUIRE_NOTHROW(generator.ReadInput(card));
  const auto &upc = generator.proc->state.lts.upc_model;
  const double original = upc->Breakup(2)->PhotoAbsorption(1.0);
  auto param = generator.proc->state.upc_param;
  std::reverse(param.emd.isotope.begin(), param.emd.isotope.end());
  const auto a = upc->Nucleus(2)->A(), z = upc->Nucleus(2)->Z();
  for (auto &row : param.emd.isotope) {
    if (row.a != a || row.z != z) { row.photo.continuum.norm *= 2.0; }
  }
  REQUIRE_NOTHROW(generator.proc->SetNuclear(param));
  CHECK(generator.proc->state.lts.upc_model->Breakup(2)->PhotoAbsorption(1.0) == Approx(original).epsilon(1.0e-12));
  for (auto &row : param.emd.isotope) {
    if (row.a == a && row.z == z) { row.photo.continuum.norm *= 2.0; }
  }
  REQUIRE_NOTHROW(generator.proc->SetNuclear(param));
  CHECK(generator.proc->state.lts.upc_model->Breakup(2)->PhotoAbsorption(1.0) > original);
}

// Trace spectator spin without changing coherent amplitudes or their azimuthal phases
TEST_CASE("Gold beams use spin-averaged nuclear currents", "[nuclear][isotope][spin][gold]") {
  const std::string model = GENERATE("GP", "TP", "yy");
  const bool reverse = GENERATE(false, true);
  auto beam = SupportedBeams()[2];
  beam.beam = reverse ? std::array<std::string, 2>{"Pb208", "Au197"}
                      : std::array<std::string, 2>{"Au197", "Pb208"};
  beam.energy = {100.0, 100.0};
  auto card = ProcessCard(beam, model == "yy" ? "yy[EPA]<F> -> mu+ mu-"
      : model + "[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  card["GENCUTS"]["<F>"]["M"] = {0.6, 0.9};
  card["FIDCUTS"]["active"] = false;
  gra::MGraniitti generator;
  REQUIRE_NOTHROW(generator.ReadInput(card));
  REQUIRE_NOTHROW(generator.proc->PrepareRun());
  auto& process = *generator.proc;
  REQUIRE(RequirePositiveEPAEvent(process) > 0.0);
  const auto physical = process.state.lts;
  REQUIRE((reverse ? physical.beam2 : physical.beam1).spinX2 > 0);
  auto scalar = physical;
  scalar.beam1.spinX2 = scalar.beam2.spinX2 = 0;
  scalar.amplitude = {};
  auto lts = physical;
  const double norm = process.ProcPtr.GetBareAmplitude2(lts);
  const auto amplitude = lts.hamp;
  REQUIRE(norm > 0.0);
  CHECK(process.ProcPtr.GetBareAmplitude2(scalar) == Approx(norm).epsilon(1.0e-12));
  REQUIRE(scalar.hamp.size() == amplitude.size());
  auto difference = scalar.hamp;
  gra::AddScaled(difference, amplitude, -1.0);
  CHECK(gra::SquaredNorm(difference) < 1.0e-22 * gra::SquaredNorm(amplitude));
  for (const double angle : {0.73, -0.73}) {
    lts.amplitude = {};
    for (auto& momentum : lts.pfinal) { momentum.RotateZ(angle); }
    lts.q1.RotateZ(angle);
    lts.q2.RotateZ(angle);
    for (auto& branch : lts.decaytree) { branch.p4.RotateZ(angle); }
    CHECK(process.ProcPtr.GetBareAmplitude2(lts) == Approx(norm).epsilon(2.0e-7));
  }
  difference = lts.hamp;
  gra::AddScaled(difference, amplitude, -1.0);
  CHECK(gra::SquaredNorm(difference) < 1.0e-12 * gra::SquaredNorm(amplitude));
  // An unspecified Regge spin is invalid for a physical ion, even after tracing its current
  process.state.lts.beam1.spinX2 = gra::aux::kNullSpinX2;
  CHECK_THROWS_AS(process.SetNuclear(process.state.upc_param), std::invalid_argument);
}


// Keep zero-absorption support when the hard nuclear transition can supply the selected neutrons
TEST_CASE("EMD conditioning respects hard breakup and proposal support", "[nuclear][EMD][importance][steering]") {
  auto card = ProcessCard(SupportedBeams()[2], "GP[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["structure"] = "nucleon";
  card["NUCLEAR"]["neutron_class"] = {"n > 0", "*"};
  card["NUCLEAR"]["photoproduction"]["target"] = {"inclusive", "coherent"};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  const auto param = generator.proc->state.upc_param;
  const gra::nuclear::MEMD emd(generator.proc->state.lts.upc_model);
  gra::MRandom random;
  random.SetSeed(9274);
  unsigned int zero = 0;
  for (unsigned int trial = 0; trial < 128; ++trial) {
    const auto state = emd.Sample(random);
    zero += !(state.excitation[0] > 0.0);
  }
  CHECK(zero > 0);
  for (const double fraction : {-0.1, 1.0, std::numeric_limits<double>::infinity()}) {
    auto invalid = param;
    invalid.emd_focus = fraction;
    CHECK_THROWS_AS(generator.proc->SetNuclear(invalid), std::invalid_argument);
  }
}

// Cover the complete daughter domain and the Ta-155 cascade that previously lacked a neutron-channel mass
TEST_CASE("Nuclear mass completion preserves data and closes evaporation cascades", "[nuclear][evaporation][mass]") {
  gra::MGraniitti generator;
  generator.ReadInput(ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-"));
  const auto& param = generator.proc->state.upc_param;
  const gra::nuclear::MMass mass(param.mass);
  auto shifted = param.mass;
  shifted.bwm[0] *= 1.01;
  const gra::nuclear::MMass shifted_mass(shifted);
  const auto& nucleus = *generator.proc->state.lts.upc_model->Nucleus(1);
  for (unsigned int a = 1; a <= nucleus.A(); ++a) {
    for (unsigned int z = 0; z <= a; ++z) {
      const double value = mass.Mass(a, z);
      REQUIRE(std::isfinite(value));
      REQUIRE(value > 0.0);
    }
  }
  // Evaluated and FRDM masses are independent of the formula for untabulated isotopes
  for (const auto& ion : std::array<std::array<unsigned int, 2>, 2>{{{197, 79}, {154, 73}}}) {
    CHECK(shifted_mass.Mass(ion[0], ion[1]) == Approx(mass.Mass(ion[0], ion[1])).epsilon(1.0e-14));
  }
  const double change = -0.1 * (shifted.bwm[0] - param.mass.bwm[0]);
  CHECK(shifted_mass.Mass(100, 2) - mass.Mass(100, 2) == Approx(change).epsilon(1.0e-10));
  CHECK_THROWS_AS(mass.Mass(3, 4), gra::PhaseSpaceFailure);
  for (const unsigned int fault : {0U, 1U, 2U}) {
    auto invalid = param.mass;
    if (fault == 0) { invalid.tables.clear(); }
    if (fault == 1) { invalid.bwm[0] = std::numeric_limits<double>::quiet_NaN(); }
    if (fault == 2) { invalid.pairing = 0.0; }
    CHECK_THROWS_AS(gra::nuclear::MMass(invalid), std::invalid_argument);
  }
  const double radius = std::sqrt(5.0 / 3.0) * nucleus.MatterDensity().Rms() / std::cbrt(nucleus.A());
  const gra::nuclear::MEvaporation decay(param.mass, radius, param.emd.gdr, param.decay_nodes, param.decay_steps);
  const double ground = decay.Mass(155, 73), excitation = 0.001905;
  REQUIRE(decay.Mass(154, 73) + decay.Mass(1, 0) - ground > excitation);
  gra::MRandom random;
  random.SetSeed(17491);
  for (const double boost : {0.0, -100.0, 100.0}) {
    for (unsigned int trial = 0; trial < 32; ++trial) {
      HepMC3::GenEvent event;
      const double m = ground + excitation, pz = boost * m;
      auto parent = std::make_shared<HepMC3::GenParticle>(HepMC3::FourVector(0, 0, pz, std::hypot(m, pz)),
                                                         gra::nuclear::EncodeNuclearPDG(155, 73), 1);
      parent->set_generated_mass(m);
      event.add_particle(parent);
      REQUIRE_NOTHROW(decay.Decay(event, parent, random));
      REQUIRE_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, parent, generator.proc->state.lts.PDG));
      REQUIRE(parent->end_vertex());
      for (const auto& particle : event.particles()) {
        if (particle->status() == 1) { CHECK(particle->pid() != 2112); }
      }
    }
  }
}

// Trace orthogonal histories at fixed impact before any hard-amplitude approximation
TEST_CASE("Resolved EMD amplitudes close the compound Poisson probabilities", "[nuclear][EMD][amplitude][closure]") {
  auto card = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  const bool conditioned = GENERATE(false, true);
  card["NUCLEAR"]["neutron_class"] = conditioned ? std::array<std::string, 2>{"n > 0", "n > 0"}
                                                   : std::array<std::string, 2>{"*", "*"};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  if (!GENERATE(false, true)) {
    auto param = generator.proc->state.upc_param;
    param.emd_focus = 0.0;
    generator.proc->SetNuclear(param);
  }
  const auto upc = generator.proc->state.lts.upc_model;
  const gra::nuclear::MEMD model(upc);
  gra::MRandom random;
  random.SetSeed(921583);
  std::array<gra::kinematics::MCW, 3> total, vacuum;
  gra::kinematics::MCW interference;
  double overlap = 0.0;
  std::array<std::size_t, 3> selected{};
  std::array<double, 3> probability{}, inclusive{};
  for (unsigned int trial = 0; trial < 2048; ++trial) {
    const auto state = model.Sample(random);
    const auto& channel = *state.channel;
    if (trial == 0) {
      for (const auto& j : indices(selected)) {
        const double b = (2.0 + 2.0 * j) * upc->Nucleus(1)->MatterDensity().Rms();
        selected[j] = std::lower_bound(channel.impact.begin(), channel.impact.end(), b) - channel.impact.begin();
        const double radius = channel.impact[selected[j]];
        probability[j] = conditioned ? 0.0 : std::exp(-upc->Breakup(1)->Mean(radius) - upc->Breakup(2)->Mean(radius));
        inclusive[j] = conditioned ? std::expm1(-upc->Breakup(1)->Mean(radius)) *
                                    std::expm1(-upc->Breakup(2)->Mean(radius)) : 1.0;
      }
      double exponent = 0.0;
      for (const int leg : {1, 2}) {
        const auto& breakup = *upc->Breakup(leg);
        const double b = channel.impact[selected[0]], other = channel.impact[selected[1]];
        const auto first = breakup.Spectrum(b), second = breakup.Spectrum(other);
        double fidelity = 0.0;
        for (const auto& k : indices(first)) { fidelity += std::sqrt(first[k] * second[k]); }
        const double common = std::sqrt(breakup.Mean(b) * breakup.Mean(other)) * fidelity;
        exponent += -0.5 * (breakup.Mean(b) + breakup.Mean(other)) +
                    (conditioned ? std::log(std::expm1(common)) : common);
      }
      overlap = std::exp(exponent);
    }
    interference.Push(channel.amplitude[selected[0]] * channel.amplitude[selected[1]]);
    for (const auto& j : indices(selected)) {
      const double weight = gra::math::pow2(channel.amplitude[selected[j]]);
      total[j].Push(weight);
      vacuum[j].Push(channel.born > 0.0 ? weight : 0.0);
    }
    REQUIRE(gra::AllFinite(channel.bare));
    REQUIRE(gra::AllFinite(channel.screened));
    if (state.excitation[0] > 0.0 || state.excitation[1] > 0.0) { CHECK(channel.born == Approx(0.0)); }
  }
  CHECK(interference.Integral() == Approx(overlap).margin(6.0 * interference.IntegralError()));
  for (const auto& j : indices(total)) {
    CHECK(total[j].Integral() == Approx(inclusive[j]).margin(6.0 * total[j].IntegralError()));
    CHECK(vacuum[j].Integral() == Approx(probability[j]).margin(6.0 * vacuum[j].IntegralError()));
  }
}

// Exercise the same excitation operator for both hard processes and both pA beam orders
TEST_CASE("Resolved EMD preserves elastic proton recoils", "[nuclear][EMD][amplitude][event]") {
  auto beam = SupportedBeams()[1];
  const bool reverse = GENERATE(false, true);
  if (reverse) {
    std::swap(beam.beam[0], beam.beam[1]);
    std::swap(beam.energy[0], beam.energy[1]);
    std::swap(beam.type[0], beam.type[1]);
  }
  const std::string process = GENERATE("yy[EPA]<F> -> mu+ mu-", "ygg[jpsi]<F> -> mu+ mu-");
  auto card = ProcessCard(beam, process);
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["FIDCUTS"]["active"] = false;
  card["GENCUTS"]["<F>"]["M"] = {3.0, 3.2};
  card["GENCUTS"]["<F>"]["Pt"] = {0.001, 0.06};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  generator.proc->PrepareRun();
  REQUIRE(RequirePositiveEPAEvent(*generator.proc) > 0.0);
  const auto& lts = generator.proc->state.lts;
  REQUIRE(lts.upc_excitation != nullptr);
  REQUIRE(lts.upc_event->Convolution());
  REQUIRE_FALSE(lts.upc_event->HadronicConvolution());
  const std::size_t proton = reverse ? 1 : 0;
  CHECK(lts.forward_emd[proton] == Approx(0.0));
  HepMC3::GenEvent event;
  REQUIRE(generator.proc->EventRecord(event));
  const auto roots = gra::nuclear::ForwardParents(event);
  CHECK(roots[proton]->end_vertex() == nullptr);
  CHECK_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, roots[1 - proton], lts.PDG));
}

// Compare the physical EMD convolution with independent impact-space Gaussian quadrature
TEST_CASE("Resolved EMD convolution converges to impact space", "[nuclear][EMD][amplitude][convergence]") {
  auto card = ProcessCard(SupportedBeams()[2], "yy[EPA]<F> -> mu+ mu-");
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  const bool survival = GENERATE(false, true);
  const bool absorbed = GENERATE(false, true);
  card["SCATTERING"]["LOOPSCREEN"] = survival;
  card["NUCLEAR"]["neutron_class"] = absorbed ? std::array<std::string, 2>{"n > 0", "n > 0"}
                                               : std::array<std::string, 2>{"*", "*"};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  const double slope = 3000.0, hbarc = gra::PDG::GeV2fm;
  const std::complex<double> phase(0.6, 0.8);
  const auto [node, weight] = gra::math::GaussLegendreRule(12, 0.0, 1.0);
  std::array<std::complex<double>, 3> amplitude{};
  std::complex<double> expected;
  const auto controls = generator.proc->state.upc_param;
  for (const auto& level : indices(amplitude)) {
    auto param = controls;
    if (!absorbed) {
      for (auto& isotope : param.emd.isotope) {
        isotope.profile.nodes = std::max<std::size_t>(2, isotope.profile.nodes >> (2 - level));
      }
    }
    param.loop.radial_intervals = 12U << level;
    param.loop.azimuth_nodes = 12U << level;
    param.convolution.smooth_kt_nodes = 256U << level;
    generator.proc->SetNuclear(param);
    const auto upc = generator.proc->state.lts.upc_model;
    gra::MRandom random;
    random.SetSeed(87124);
    const gra::nuclear::MEMD model(upc);
    auto state = model.Sample(random);
    if (!absorbed) {
      while (!(state.channel->born > 0.0)) { state = model.Sample(random); }
    }
    const auto& channel = *state.channel;
    const auto event = upc->WithExcitation(state.channel, survival);
    gra::nuclear::ScreenPoint born;
    born.amplitude = {phase};
    gra::nuclear::MUPCScreen screen(*event, {}, born);
    for (const auto& k : event->Nodes(born.transfer)) {
      auto shifted = born;
      shifted.amplitude[0] *= std::exp(-2.0 * slope * k.kt * k.kt);
      screen.Add(k, shifted);
    }
    amplitude[level] = screen.Result().amplitude.front();
    // The first bin is constant and subsequent bins are linear in log(b)
    double integral = 0.0;
    for (const auto& i : indices(channel.impact)) {
      const double low = i == 0 ? 0.0 : channel.impact[i - 1], high = channel.impact[i];
      const double upper = upc->SurvivalAmp(high) * channel.amplitude[i];
      const double lower = i == 0 ? upper : upc->SurvivalAmp(low) * channel.amplitude[i - 1];
      for (const auto& j : indices(node)) {
        const double b = low + (high - low) * node[j];
        const double fraction = i == 0 ? 0.0 : std::log(b / low) / std::log(high / low);
        integral += weight[j] * (high - low) * b * (lower + fraction * (upper - lower)) *
                    std::exp(-b * b / (8.0 * slope * hbarc * hbarc));
      }
    }
    expected = phase * integral / (4.0 * slope * hbarc * hbarc);
    // The unique vacuum history permits profile refinement after removing its proposal normalization
    if (!absorbed) {
      amplitude[level] /= channel.born;
      expected /= channel.born;
    }
    CAPTURE(level, survival, absorbed, amplitude[level], expected);
    // Relative amplitude accuracy also controls the cross section to twice this order
    CHECK(std::abs(amplitude[level] - expected) < (level == 0 ? 0.03 : 0.01) * std::abs(expected));
  }
  CHECK(std::abs(amplitude.back() - amplitude[1]) < 0.005 * std::abs(expected));
}

// Combine target knockout and EMD before the remnant evaporation cascade
TEST_CASE("Incoherent photoproduction retains resolved EMD particles", "[nuclear][EMD][incoherent][event]") {
  const std::string process = GENERATE("ygg[jpsi]<F> -> mu+ mu-", "GP[RES]<F> -> pi+ pi- @RES{rho_770:1}");
  auto card = ProcessCard(SupportedBeams()[2], process);
  card["NUCLEAR"]["emd"] = true;
  card["NUCLEAR"]["fragmentation"] = "internal";
  card["NUCLEAR"]["structure"] = "nucleon";
  card["NUCLEAR"]["screening"] = GENERATE("optical", "mc_ggcf");
  card["NUCLEAR"]["photoproduction"]["target"] = {"incoherent", "incoherent"};
  card["SCATTERING"]["LOOPSCREEN"] = true;
  card["FIDCUTS"]["active"] = false;
  card["GENCUTS"]["<F>"]["Pt"] = {0.2, 0.5};
  card["GENCUTS"]["<F>"]["M"] = process.starts_with("ygg") ? std::array<double, 2>{3.0, 3.2}
                                                             : std::array<double, 2>{0.6, 0.9};
  gra::MGraniitti generator;
  generator.ReadInput(card);
  auto param = generator.proc->state.upc_param;
  param.config.count = param.current_count = 2;
  generator.proc->SetNuclear(param);
  generator.proc->PrepareRun();
  REQUIRE(RequirePositiveEPAEvent(*generator.proc, true) > 0.0);
  const auto& lts = generator.proc->state.lts;
  REQUIRE(lts.upc_excitation != nullptr);
  REQUIRE(lts.upc_event->HadronicConvolution());
  HepMC3::GenEvent event;
  REQUIRE(generator.proc->EventRecord(event));
  unsigned int incoherent = 0;
  for (const auto& root : gra::nuclear::ForwardParents(event)) {
    CHECK_NOTHROW(gra::nuclear::ValidateNuclearDecay(event, root, lts.PDG));
    if (root->attribute<HepMC3::StringAttribute>("graniitti_upc_final_sector")->value() != "incoherent") { continue; }
    ++incoherent;
    REQUIRE(root->end_vertex() != nullptr);
    const auto& particles = root->end_vertex()->particles_out();
    CHECK(std::any_of(particles.begin(), particles.end(), [](const auto& p) { return p->pid() == 2112 || p->pid() == 2212; }));
  }
  CHECK(incoherent == 1);
}
