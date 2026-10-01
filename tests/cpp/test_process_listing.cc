// Unit tests for algorithmic process registries
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iterator>
#include <limits>
#include <set>
#include <stdexcept>
#include <string>
#include <type_traits>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Durham/AMP_MG5_DurhamRegistry.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_PartonRegistry.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_jj.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_w.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_z.h"
#include "Graniitti/Amplitude/MG5/Parton/AMP_MG5_pp_zj.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_jj.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_zjj.h"
#include "Graniitti/Amplitude/MG5/Photon/MG5_YY_ZJJ/Processes.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Process.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_MadLoop.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"
#include "Graniitti/MGraniitti.h"
#include "Graniitti/MModelCache.h"
#include "Graniitti/Math/MFloat.h"
#include "Graniitti/Particle/MResonance.h"
#include "Graniitti/Kinematics/MCentral.h"
#include "Graniitti/Kinematics/MCollinear.h"
#include "Graniitti/Process/MDecayTree.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/Kinematics/MHardDiffraction.h"
#include "Graniitti/Kinematics/MQuasiElastic.h"
#include "Graniitti/Process/MSubProc.h"
#include "Graniitti/Process/MProcessTable.h"

// Libraries
#include <catch.hpp>

namespace {

// Bind the immutable TUNE0 cache used by direct generated amplitude fixtures
void BindProcessModelCache(gra::LORENTZSCALAR &lts) {
  lts.model_cache =
      std::make_shared<gra::MModelCache>(gra::MModelTune::Load(gra::ResolveModelDataFile("TUNE0", "GENERAL.json")));
}


// Compute true when one process row contains both requested fragments
bool ContainsProcessRow(const std::vector<std::vector<std::string>> &rows, const std::string &command,
                        const std::string &matrix_element) {
  return std::any_of(rows.begin(), rows.end(), [&](const auto &row) {
    return row.size() == 4 && row[0].find(command) != std::string::npos &&
           row[3].find(matrix_element) != std::string::npos;
  });
}

// Replace one exact phase-space token in a process command
std::string ReplaceProcessMode(const std::string &command, const std::string &current_mode,
                               const std::string &new_mode) {
  const std::string token    = "<" + current_mode + ">";
  const std::size_t position = command.find(token);
  if (position == std::string::npos) {
    throw std::logic_error("ReplaceProcessMode: missing token " + token + " in " + command);
  }
  return command.substr(0, position) + "<" + new_mode + ">" + command.substr(position + token.size());
}

// Replace phase-space tokens and sort one generated final-state row set
std::vector<std::vector<std::string>> ReplaceFinalStateRowModes(std::vector<std::vector<std::string>> rows,
                                                                const std::string                    &current_mode,
                                                                const std::string                    &new_mode) {
  for (auto &row : rows) {
    if (row.size() != 4) { throw std::logic_error("ReplaceFinalStateRowModes: expected four columns"); }
    row[0] = ReplaceProcessMode(row[0], current_mode, new_mode);
  }
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
  return rows;
}

// Append one generated final-state row set and restore sorted uniqueness
void AppendFinalStateRows(std::vector<std::vector<std::string>> &rows, const gra::MSubProc &subprocess) {
  const auto additional = subprocess.SupportedFinalStateRows();
  rows.insert(rows.end(), additional.begin(), additional.end());
  std::sort(rows.begin(), rows.end());
  rows.erase(std::unique(rows.begin(), rows.end()), rows.end());
}

// Compute sorted exact process commands from one phase-space registry
std::vector<std::string> ProcessCommandKeys(const gra::MSubProc &subprocess) {
  std::vector<std::string> commands;
  commands.reserve(subprocess.ProcessRegistry.size());
  for (const auto &entry : subprocess.ProcessRegistry) { commands.push_back(entry.first); }
  std::sort(commands.begin(), commands.end());
  return commands;
}

// Append exact process commands from one phase-space registry
void AppendProcessCommandKeys(std::vector<std::string> &commands, const gra::MSubProc &subprocess) {
  const auto additional = ProcessCommandKeys(subprocess);
  commands.insert(commands.end(), additional.begin(), additional.end());
  std::sort(commands.begin(), commands.end());
  commands.erase(std::unique(commands.begin(), commands.end()), commands.end());
}

// Build one stable branch for generated amplitude topology tests
gra::MDecayBranch StableBranch(int pdg) {
  gra::MDecayBranch branch;
  branch.p.pdg = pdg;
  return branch;
}

// Build one internal branch with ordered daughters for topology tests
gra::MDecayBranch CascadeBranch(int pdg, std::vector<gra::MDecayBranch> daughters) {
  gra::MDecayBranch branch;
  branch.p.pdg = pdg;
  branch.legs  = std::move(daughters);
  return branch;
}

// Build one stable branch with fixed benchmark momentum
gra::MDecayBranch MomentumBranch(int pdg, const gra::M4Vec &p4) {
  gra::MDecayBranch branch = StableBranch(pdg);
  branch.p4                = p4;
  return branch;
}

// Build one Z cascade with fixed ordered dimuon momenta
gra::MDecayBranch ZMomentumBranch(const gra::M4Vec &mup, const gra::M4Vec &mum) {
  gra::MDecayBranch branch = CascadeBranch(23, {MomentumBranch(-13, mup), MomentumBranch(13, mum)});
  branch.p4                = mup + mum;
  return branch;
}

// Build a fixed direct dimuon hard-process point
gra::LORENTZSCALAR DirectZKinematics() {
  gra::LORENTZSCALAR lts;
  const double mass = gra::AMP_MG5_pp_z().Particles().at(13).mass;
  const double momentum = std::sqrt(100.0 * 100.0 - mass * mass);
  lts.q1        = gra::M4Vec(0.0, 0.0, 100.0, 100.0);
  lts.q2        = gra::M4Vec(0.0, 0.0, -100.0, 100.0);
  lts.id1       = 2;
  lts.id2       = -2;
  lts.alphaQCD  = 0.118;
  lts.decaytree = {
      MomentumBranch(-13, gra::M4Vec(0.8 * momentum, 0.0, 0.6 * momentum, 100.0)),
      MomentumBranch(13, gra::M4Vec(-0.8 * momentum, 0.0, -0.6 * momentum, 100.0)),
  };
  return lts;
}

// Build a fixed Z plus one-parton hard-process point
gra::LORENTZSCALAR ZJetKinematics() {
  gra::LORENTZSCALAR lts;
  const double mass = gra::AMP_MG5_pp_zj().Particles().at(13).mass;
  const double transverse = std::sqrt(80.0 * 80.0 - mass * mass);
  lts.q1        = gra::M4Vec(0.0, 0.0, 160.0, 160.0);
  lts.q2        = gra::M4Vec(0.0, 0.0, -160.0, 160.0);
  lts.id1       = 21;
  lts.id2       = 2;
  lts.alphaQCD  = 0.118;
  lts.decaytree = {
      ZMomentumBranch(gra::M4Vec(transverse, 0.0, 60.0, 100.0), gra::M4Vec(-transverse, 0.0, 60.0, 100.0)),
      MomentumBranch(2, gra::M4Vec(0.0, 0.0, -120.0, 120.0)),
  };
  return lts;
}

// Build a fixed Z plus two-parton hard-process point
gra::LORENTZSCALAR ZDijetKinematics() {
  gra::LORENTZSCALAR lts;
  const double mass = gra::AMP_MG5_yy_zjj().Particles().at(13).mass;
  const double momentum = std::sqrt(100.0 * 100.0 - mass * mass);
  lts.q1        = gra::M4Vec(0.0, 0.0, 220.0, 220.0);
  lts.q2        = gra::M4Vec(0.0, 0.0, -220.0, 220.0);
  lts.id1       = 21;
  lts.id2       = 21;
  lts.alphaQCD  = 0.118;
  lts.decaytree = {
      ZMomentumBranch(gra::M4Vec(momentum, 0.0, 0.0, 100.0), gra::M4Vec(-momentum, 0.0, 0.0, 100.0)),
      MomentumBranch(2, gra::M4Vec(0.0, 120.0, 0.0, 120.0)),
      MomentumBranch(-2, gra::M4Vec(0.0, -120.0, 0.0, 120.0)),
  };
  return lts;
}

// Build a fixed direct two-parton hard-process point
gra::LORENTZSCALAR DijetKinematics() {
  gra::LORENTZSCALAR lts;
  lts.q1        = gra::M4Vec(0.0, 0.0, 160.0, 160.0);
  lts.q2        = gra::M4Vec(0.0, 0.0, -160.0, 160.0);
  lts.id1       = 21;
  lts.id2       = 21;
  lts.alphaQCD  = 0.118;
  lts.decaytree = {
      MomentumBranch(2, gra::M4Vec(160.0, 0.0, 0.0, 160.0)),
      MomentumBranch(-2, gra::M4Vec(-160.0, 0.0, 0.0, 160.0)),
  };
  return lts;
}

// Flatten the first cascade while preserving stable-leaf momenta and order
gra::LORENTZSCALAR FlattenFirstCascade(gra::LORENTZSCALAR lts) {
  if (lts.decaytree.empty() || lts.decaytree.front().legs.empty()) {
    throw std::invalid_argument("FlattenFirstCascade: first branch is not a cascade");
  }
  std::vector<gra::MDecayBranch> flattened = lts.decaytree.front().legs;
  flattened.insert(flattened.end(), lts.decaytree.begin() + 1, lts.decaytree.end());
  lts.decaytree = std::move(flattened);
  return lts;
}

// Check stale prepared hard-process state cannot bypass an exact topology guard
template <class Amplitude>
void CheckPreparedTopologyRejection(Amplitude &amplitude, gra::LORENTZSCALAR valid, gra::LORENTZSCALAR wrong) {
  BindProcessModelCache(valid);
  wrong.model_cache       = valid.model_cache;
  const auto valid_result = amplitude.EvaluatePrepared(valid, 0.118);
  REQUIRE(valid_result.Valid());
  REQUIRE(valid_result.amp2 > 0.0);
  REQUIRE(valid_result.contributing_subprocesses > 0);
  REQUIRE_FALSE(valid_result.color_flows.empty());

  wrong.hamp              = {1.0};
  const auto wrong_result = amplitude.EvaluatePrepared(wrong, 0.118);
  CHECK_FALSE(wrong_result.Valid());
  CHECK(wrong_result.status == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  CHECK(wrong.hamp.empty());
  CHECK(wrong_result.contributing_subprocesses == 0);
  CHECK(wrong_result.color_flows.empty());

  REQUIRE(amplitude.EvaluatePrepared(valid, 0.118).amp2 > 0.0);
  wrong.hamp = {1.0};
  CHECK(amplitude.Prepare(wrong, 0.118) == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  CHECK(wrong.hamp.empty());
  CHECK_FALSE(amplitude.EvaluatePrepared(wrong, 0.118).Valid());
}

// Insert one particle token into the local process parser table
void InsertProcessParticle(gra::MPDG &pdg, int code, const std::string &name) {
  gra::MParticle particle;
  particle.pdg  = code;
  particle.name = name;
  pdg.PDG_table.emplace(code, std::move(particle));
}

// Build a self-contained particle table for every registered final-state token
gra::MPDG ProcessPDGTable() {
  gra::MPDG                                      pdg;
  const std::vector<std::pair<int, std::string>> particles = {
      {1, "d"},   {-1, "d~"},    {2, "u"},     {-2, "u~"},    {3, "s"},    {-3, "s~"},
      {4, "c"},   {-4, "c~"},    {5, "b"},     {-5, "b~"},    {6, "t"},    {-6, "t~"},
      {11, "e-"}, {-11, "e+"},   {12, "ve"},   {-12, "ve~"},  {13, "mu-"}, {-13, "mu+"},
      {14, "vm"}, {-14, "vm~"},  {15, "tau-"}, {-15, "tau+"}, {16, "vt"},  {-16, "vt~"},
      {21, "g"},  {22, "gamma"}, {23, "Z"},    {24, "W+"},    {-24, "W-"}, {gra::PDG::PDG_hard_jet, "j"},
  };
  for (const auto &[code, name] : particles) { InsertProcessParticle(pdg, code, name); }
  return pdg;
}

// Parse one registry final state into an event decay tree
std::vector<gra::MDecayBranch> ProcessTree(const gra::MPDG &pdg, const gra::amplitude::Process &process) {
  std::vector<gra::MDecayBranch> tree;
  pdg.TokenizeProcess(process.final_state_syntax, 0, tree);
  return tree;
}

// Compute the first stable leaf of one mutable decay branch
gra::MDecayBranch *FirstStableLeaf(gra::MDecayBranch &branch) {
  if (branch.legs.empty()) { return &branch; }
  for (auto &daughter : branch.legs) {
    if (auto *leaf = FirstStableLeaf(daughter); leaf != nullptr) { return leaf; }
  }
  return nullptr;
}

// Minimal matrix element used to test subprocess sum event keys
struct SubprocessTestMatrixElement {
  const std::vector<double> masses = {0.0, 0.0, 0.0, 0.0};

  // Compute the fixed massless model masses used by these event-key fixtures
  const std::vector<double> &getMasses() const { return masses; }

  // Accept the process-local QED coupling selector used by generated families
  void setAlphaQEDZero(bool) noexcept {}
};

// Change one stable particle so an otherwise valid topology must not match
std::vector<gra::MDecayBranch> MismatchedTree(std::vector<gra::MDecayBranch> tree) {
  for (auto &branch : tree) {
    if (auto *leaf = FirstStableLeaf(branch); leaf != nullptr) {
      const int replacement = std::abs(leaf->p.pdg) == gra::PDG::PDG_gamma ? 13 : gra::PDG::PDG_gamma;
      leaf->p               = gra::MParticle{};
      leaf->p.pdg           = replacement;
      leaf->p.name          = replacement == gra::PDG::PDG_gamma ? "gamma" : "mu-";
      return tree;
    }
  }
  throw std::invalid_argument("MismatchedTree: decay tree has no stable leaf");
}

// Compute all converter-generated entries from one process sequence
std::vector<gra::amplitude::Process> GeneratedProcesses(const std::vector<gra::amplitude::Process> &processes) {
  std::vector<gra::amplitude::Process> generated;
  std::copy_if(processes.begin(), processes.end(), std::back_inserter(generated),
               [](const auto &process) { return gra::amplitude::IsGenerated(process.matrix_element_form); });
  return generated;
}

// Build one self-contained generated process for selector matching tests
gra::amplitude::Process GeneratedProcess(const std::string &name, const std::string &signature,
                                         gra::AmplitudeTopology topology, std::vector<int> stable_pdgs,
                                         gra::amplitude::TopologyMode topology_mode) {
  return {"TEST",
          name,
          signature,
          signature,
          signature,
          "test > " + signature,
          std::move(topology),
          {},
          std::move(stable_pdgs),
          {gra::DecayType::Full},
          gra::amplitude::MatrixElementForm::Generated,
          topology_mode};
}

// Check that one process route follows one complete process definition
void CheckProcessRoute(const gra::amplitude::ProcessDefinition &definition, gra::MProc &route, const gra::MPDG &pdg,
                       bool strict_process_set) {
  const auto processes = definition.Processes();
  REQUIRE_FALSE(processes.empty());
  const auto supported_processes = route.Processes();
  if (strict_process_set) {
    REQUIRE(supported_processes == processes);
  } else {
    REQUIRE(GeneratedProcesses(supported_processes) == processes);
    REQUIRE(supported_processes.size() > processes.size());
  }

  for (const auto &process : processes) {
    INFO("Process route " << route.ISTATE << "[" << route.CHANNEL << "] -> " << process.process_name);
    const auto tree = ProcessTree(pdg, process);
    REQUIRE_FALSE(tree.empty());
    REQUIRE(definition.MatchProcess(tree) == process);
    REQUIRE(route.MatchProcess(tree) == process);

    gra::LORENTZSCALAR lts;
    lts.decaytree = tree;
    REQUIRE(route.DecayStructureFor(lts) == process.decay_structure);

    const auto wrong_tree = MismatchedTree(tree);
    CHECK_FALSE(definition.MatchProcess(wrong_tree).has_value());
    CHECK_FALSE(route.MatchProcess(wrong_tree).has_value());
    lts.decaytree = wrong_tree;
    CHECK_THROWS_AS(definition.DecayStructureFor(lts), std::invalid_argument);
    CHECK_THROWS_AS(route.DecayStructureFor(lts), std::invalid_argument);
  }
}

}  // namespace

// Check finite-Nc color contraction and nonthrowing invalid-metric handling
TEST_CASE("MG5 color metric contraction fails through its return value", "[MG2GRA][color]") {
  const std::vector<std::complex<double>> jamp         = {{1.0, 1.0}, {2.0, -1.0}};
  const std::vector<double>               denominators = {1.0, 1.0};
  const std::vector<std::vector<double>>  metric       = {{2.0, 1.0}, {1.0, 2.0}};

  const auto amplitudes = gra::mg5helas::ColorMetricAmplitudes(jamp, denominators, metric);
  REQUIRE(amplitudes.size() == jamp.size());
  const std::vector<std::complex<double>> expected_amplitudes = {std::sqrt(2.0) * jamp[0] + jamp[1] / std::sqrt(2.0),
                                                                 std::sqrt(3.0 / 2.0) * jamp[1]};
  for (std::size_t component = 0; component < amplitudes.size(); ++component) {
    CHECK(amplitudes[component].real() == Approx(expected_amplitudes[component].real()).margin(1.0e-12));
    CHECK(amplitudes[component].imag() == Approx(expected_amplitudes[component].imag()).margin(1.0e-12));
  }
  double norm = 0.0;
  for (const auto &amplitude : amplitudes) { norm += std::norm(amplitude); }
  const double expected =
      std::real(std::conj(jamp[0]) * (2.0 * jamp[0] + jamp[1]) + std::conj(jamp[1]) * (jamp[0] + 2.0 * jamp[1]));
  CHECK(norm == Approx(expected).epsilon(1.0e-12));

  CHECK(gra::mg5helas::ColorMetricAmplitudes(jamp, {1.0}, metric).empty());
  CHECK(gra::mg5helas::ColorMetricAmplitudes(jamp, denominators, {{1.0, 2.0}, {2.0, 1.0}}).empty());
  CHECK(gra::mg5helas::ColorMetricAmplitudes(jamp, denominators, {{2.0, 1.0}, {1.25, 2.0}}).empty());
  auto nonfinite_jamp = jamp;
  nonfinite_jamp[0]   = {std::numeric_limits<double>::quiet_NaN(), 0.0};
  CHECK(gra::mg5helas::ColorMetricAmplitudes(nonfinite_jamp, denominators, metric).empty());
}

// Check typed, finite and fixed-rank SLHA parsing without model-card constants
TEST_CASE("MG5 SLHA reader handles general block ranks and atomic failure", "[MG2GRA][SLHA]") {
  SLHABlock block("tensor");
  CHECK(block.get_indices() == 0);
  block.set_entry({1, 2, 3}, 4.5);
  CHECK(block.get_indices() == 3);
  CHECK(block.get_entry({1, 2, 3}) == Approx(4.5));
  REQUIRE(block.find_entry({1, 2, 3}).has_value());
  CHECK(*block.find_entry({1, 2, 3}) == Approx(4.5));
  CHECK_FALSE(block.find_entry({3, 2, 1}).has_value());
  CHECK_THROWS_AS(block.set_entry({1, 2}, 4.5), std::invalid_argument);
  CHECK_THROWS_AS(block.set_entry({1, 2, 3}, std::numeric_limits<double>::quiet_NaN()), std::invalid_argument);
  SLHABlock scalar("alpha");
  scalar.set_entry({}, -0.0078);
  CHECK(scalar.get_indices() == 0);
  CHECK(scalar.get_entry({}) == Approx(-0.0078));

  const std::filesystem::path directory = std::filesystem::path("tmp") / "test_mg5_slha_reader";
  std::filesystem::create_directories(directory);
  const std::filesystem::path valid_path = directory / "valid.dat";
  {
    std::ofstream output(valid_path);
    REQUIRE(output.good());
    output << "BLOCK TENSOR\n"
           << std::string(300, ' ') << "1 2 3 4.5D+00\n"
           << "BLOCK ALPHA\n"
           << "-7.8D-03\n"
           << "DECAY 23 2.5D+00\n";
  }
  SLHAReader card(valid_path.string());
  CHECK(card.get_block_entry("tensor", {1, 2, 3}) == Approx(4.5));
  CHECK(card.get_block_entry("ALPHA", std::vector<int>{}) == Approx(-0.0078));
  CHECK(card.get_block_entry("decay", 23) == Approx(2.5));
  REQUIRE(card.find_block_entry("tensor", {1, 2, 3}).has_value());
  CHECK(*card.find_block_entry("tensor", {1, 2, 3}) == Approx(4.5));
  REQUIRE(card.find_block_entry("DECAY", 23).has_value());
  CHECK(*card.find_block_entry("DECAY", 23) == Approx(2.5));
  CHECK_FALSE(card.find_block_entry("mass", 23).has_value());

  const std::filesystem::path invalid_path = directory / "invalid.dat";
  {
    std::ofstream output(invalid_path);
    REQUIRE(output.good());
    output << "BLOCK BROKEN\n1 not-a-finite-value\n";
  }
  CHECK_THROWS_AS(card.read_slha_file(invalid_path.string()), std::runtime_error);
  CHECK(card.get_block_entry("tensor", {1, 2, 3}) == Approx(4.5));
}

// Check generated decay branches use their live family card parameters
TEST_CASE("Generated decay proposals use family parameter cards", "[process][MG5][SLHA][decay]") {
  const std::string z_card = gra::aux::ResolveProjectPath("MG5cards/Parton/MG5_PP_ZJ/param_card.dat");
  const SLHAReader  z_values(z_card);
  const auto        z_mass  = z_values.find_block_entry("mass", 23);
  const auto        z_width = z_values.find_block_entry("decay", 23);
  REQUIRE(z_mass.has_value());
  REQUIRE(z_width.has_value());

  auto electron_tree       = std::vector<gra::MDecayBranch>{CascadeBranch(23, {StableBranch(-11), StableBranch(11)}),
                                                            StableBranch(gra::PDG::PDG_hard_jet)};
  electron_tree[0].p.mass  = 17.0;
  electron_tree[0].p.width = 3.0;
  electron_tree[0].legs[0].p.mass = 0.123;
  electron_tree[0].legs[1].p.mass = 0.456;
  const auto z_parameters = gra::CreatePartonMG5Process("MG5_PP_ZJ")->Particles();
  gra::SynchronizeMG5DecayParameters(electron_tree, z_parameters);
  CHECK(electron_tree[0].p.mass == Approx(*z_mass));
  CHECK(electron_tree[0].p.width == Approx(*z_width));
  CHECK(electron_tree[0].legs[0].p.mass == Approx(z_parameters.at(11).mass));
  CHECK(electron_tree[0].legs[1].p.mass == Approx(z_parameters.at(11).mass));

  auto muon_tree              = std::vector<gra::MDecayBranch>{CascadeBranch(23, {StableBranch(-13), StableBranch(13)}),
                                                               StableBranch(gra::PDG::PDG_hard_jet)};
  muon_tree[0].legs[0].p.mass = 0.789;
  muon_tree[0].legs[1].p.mass = 0.987;
  gra::SynchronizeMG5DecayParameters(muon_tree, z_parameters);
  CHECK(muon_tree[0].p.mass == Approx(*z_mass));
  CHECK(muon_tree[0].p.width == Approx(*z_width));
  CHECK(muon_tree[0].legs[0].p.mass == Approx(z_parameters.at(13).mass));
  CHECK(muon_tree[0].legs[1].p.mass == Approx(z_parameters.at(13).mass));

  const gra::amplitude::ProcessRegistry zj(gra::amplitude::Processes("MG5_PP_ZJ"));
  REQUIRE(zj.MatchProcess(electron_tree).has_value());
  REQUIRE(zj.MatchProcess(muon_tree).has_value());

  const std::string ww_card = gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_WW/param_card.dat");
  const SLHAReader  ww_values(ww_card);
  const auto        w_mass  = ww_values.find_block_entry("mass", 24);
  const auto        w_width = ww_values.find_block_entry("decay", 24);
  REQUIRE(w_mass.has_value());
  REQUIRE(w_width.has_value());
  auto ww_tree = std::vector<gra::MDecayBranch>{CascadeBranch(24, {StableBranch(-11), StableBranch(12)}),
                                                CascadeBranch(-24, {StableBranch(13), StableBranch(-14)})};
  const auto w_parameters = gra::CreatePhotonMG5Process("MG5_YY_WW")->Particles();
  gra::SynchronizeMG5DecayParameters(ww_tree, w_parameters);
  for (const auto &branch : ww_tree) {
    CHECK(branch.p.mass == Approx(w_parameters.at(24).mass));
    CHECK(branch.p.width == Approx(*w_width));
  }

  auto resonance                = electron_tree[0];
  resonance.p4                  = gra::M4Vec(0.0, 0.0, 0.0, *z_mass + 0.5 * *z_width);
  resonance.W_event             = 1.0;
  resonance.mass_proposal_norm  = 2.75;
  resonance.mass_proposal = gra::MassProposal::BreitWigner;
  const auto   phase_space      = gra::decay::GeneratedPhaseSpace({resonance}, false);
  const double propagator_density =
      std::norm(gra::resonance::FixedWidthLineShape(resonance.p4.M2(), resonance.p.mass, resonance.p.width));
  REQUIRE(propagator_density > 0.0);
  CHECK(phase_space.proposal_volume * propagator_density == Approx(resonance.mass_proposal_norm).epsilon(1.0e-13));
}

// Check the signed event-count modes accepted by the generator interface
TEST_CASE("Generator event count distinguishes proposal-only integration", "[generator][integration]") {
  gra::MGraniitti generator;

  generator.SetNumberOfEvents(-1);
  CHECK(generator.IntegrationProposalOnly());
  CHECK(generator.GetNumberOfEvents() == -1);

  generator.SetNumberOfEvents(0);
  CHECK_FALSE(generator.IntegrationProposalOnly());
  CHECK(generator.GetNumberOfEvents() == 0);

  generator.SetNumberOfEvents(10);
  CHECK_FALSE(generator.IntegrationProposalOnly());
  CHECK(generator.GetNumberOfEvents() == 10);
  CHECK_THROWS_AS(generator.SetNumberOfEvents(-2), std::invalid_argument);
}

// Check generated rows and process command keys against every phase-space
// registry
TEST_CASE("Process registry covers generated MG5 processes algorithmically", "[process][MG5]") {
  const gra::MFactorized      factorized_phase_space;
  const gra::MCentral         central_phase_space;
  const gra::MQuasiElastic    quasielastic_phase_space;
  const gra::MCollinear       collinear_phase_space;
  const gra::MHardDiffraction hard_diffraction_phase_space;
  const gra::MSubProc        &factorized       = factorized_phase_space.ProcPtr;
  const gra::MSubProc        &continuum        = central_phase_space.ProcPtr;
  const gra::MSubProc        &quasielastic     = quasielastic_phase_space.ProcPtr;
  const gra::MSubProc        &collinear        = collinear_phase_space.ProcPtr;
  const gra::MSubProc        &hard_diffraction = hard_diffraction_phase_space.ProcPtr;

  gra::MGraniitti generator;
  const auto      rows = generator.GetFinalStateRows();
  REQUIRE_FALSE(rows.empty());

  for (const auto &row : rows) {
    REQUIRE(row.size() == 4);
    CHECK_FALSE(row[0].empty());
    CHECK_FALSE(row[1].empty());
    CHECK_FALSE(row[2].empty());
    CHECK_FALSE(row[3].empty());
    CHECK(row[0].find(" / ") == std::string::npos);
  }

  const auto factorized_rows = factorized.SupportedFinalStateRows();
  const auto continuum_rows  = continuum.SupportedFinalStateRows();
  REQUIRE(ReplaceFinalStateRowModes(factorized_rows, "F", "MODE") ==
          ReplaceFinalStateRowModes(continuum_rows, "C", "MODE"));

  auto expected_rows = ReplaceFinalStateRowModes(factorized_rows, "F", "F|C");
  AppendFinalStateRows(expected_rows, quasielastic);
  AppendFinalStateRows(expected_rows, collinear);
  const gra::MSubProc hard_f({"IPp", "IPIP"}, "F");
  const auto hard_rows = ReplaceFinalStateRowModes(hard_f.SupportedFinalStateRows(), "F", "F|C");
  expected_rows.insert(expected_rows.end(), hard_rows.begin(), hard_rows.end());
  std::sort(expected_rows.begin(), expected_rows.end());
  CHECK(rows == expected_rows);

  CHECK(ContainsProcessRow(rows, "yy[WW]<F|C> -> W+ > {e+ ve} W- > {e- ve~}", "a a > w+ w-, w+ > e+ ve, w- > e- ve~"));
  CHECK(ContainsProcessRow(rows, "yy_LUX[WW]<P> -> W+ > {tau+ vt} W- > {mu- vm~}",
                           "a a > w+ w-, w+ > ta+ vt, w- > mu- vm~"));
  CHECK(ContainsProcessRow(rows, "gg[QCD]<F|C> -> g g g g", "g g > g g g g QED=0"));
  CHECK(ContainsProcessRow(rows, "yy[Zjj]<F|C> -> Z > {mu+ mu-} u u~", "a a > z u u~, z > mu+ mu-"));
  CHECK(ContainsProcessRow(rows, "yy[jj]<F|C> -> j j", "a a > j j"));
  CHECK(ContainsProcessRow(rows, "IPp[Z]<F|C> -> mu+ mu-", "p p > mu+ mu-"));
  CHECK(ContainsProcessRow(rows, "IPIP[Z]<F|C> -> tau+ tau-", "p p > ta+ ta-"));
  CHECK(ContainsProcessRow(rows, "IPp[Zj]<F|C> -> Z > {mu+ mu-} j", "p p > z j, z > mu+ mu-"));
  CHECK(ContainsProcessRow(rows, "IPIP[jj]<F|C> -> j j", "p p > j j QED=0"));
  CHECK(ContainsProcessRow(rows, "IPp[W]<F|C> -> mu+ vm", "p p > mu+ vm"));
  CHECK(ContainsProcessRow(rows, "IPIP[W]<F|C> -> tau- vt~", "p p > ta- vt~"));

  std::vector<std::string> commands = generator.GetProcessNumbers();

  auto expected_commands = ProcessCommandKeys(factorized);
  AppendProcessCommandKeys(expected_commands, continuum);
  AppendProcessCommandKeys(expected_commands, quasielastic);
  AppendProcessCommandKeys(expected_commands, collinear);
  AppendProcessCommandKeys(expected_commands, hard_diffraction);
  REQUIRE(commands.size() == expected_commands.size());
  std::sort(commands.begin(), commands.end());
  commands.erase(std::unique(commands.begin(), commands.end()), commands.end());
  REQUIRE(commands == expected_commands);

  for (const auto &entry : hard_diffraction.ProcessRegistry) {
    const std::string removed_command = ReplaceProcessMode(entry.first, entry.first.find("<C>") != std::string::npos ? "C" : "F", "H");
    CHECK(std::binary_search(commands.begin(), commands.end(), entry.first));
    CHECK_FALSE(std::binary_search(commands.begin(), commands.end(), removed_command));
  }
}

// Check every process in the shared F and C registry owns a process definition
TEST_CASE("Every shared central process owns a process definition", "[process][process][central]") {
  const gra::MFactorized factorized_phase_space;
  const gra::MCentral    central_phase_space;
  const gra::MSubProc   &factorized = factorized_phase_space.ProcPtr;
  const gra::MSubProc   &continuum  = central_phase_space.ProcPtr;

  REQUIRE_FALSE(factorized.ProcessRegistry.empty());
  REQUIRE(factorized.ProcessRegistry.size() == continuum.ProcessRegistry.size());

  std::size_t advertised_processes = 0;
  for (const auto &process : factorized.CreateAllProcesses()) {
    const std::string factorized_command = process->ISTATE + "[" + process->CHANNEL + "]<F>";
    if (!factorized.ProcessExist(factorized_command)) { continue; }

    ++advertised_processes;
    const std::string continuum_command = process->ISTATE + "[" + process->CHANNEL + "]<C>";
    INFO("Shared central process " << factorized_command << " / " << continuum_command);
    REQUIRE(continuum.ProcessExist(continuum_command));

    const auto processes = process->Processes();
    REQUIRE_FALSE(processes.empty());
  }
  REQUIRE(advertised_processes == factorized.ProcessRegistry.size());
}

// Check the QED continuum decay structure through both shared phase spaces
TEST_CASE("QED direct pairs expose one strict shared F and C decay structure", "[process][process][central][QED]") {
  const gra::DecayStructure expected{gra::DecayType::Full};

  for (const std::string mode : {"F", "C"}) {
    INFO("Central phase-space mode " << mode);
    gra::MSubProc process({"yy"}, mode);
    process.Initialize("yy", "QED");

    gra::LORENTZSCALAR direct;
    direct.decaytree             = {StableBranch(-13), StableBranch(13)};
    direct.decaytree[0].p.spinX2 = 1;
    direct.decaytree[1].p.spinX2 = 1;
    REQUIRE(process.DecayStructureFor(direct) == expected);
    REQUIRE(direct.decay_structure == expected);

    gra::LORENTZSCALAR wrong_spin    = direct;
    wrong_spin.decaytree[0].p.spinX2 = 0;
    CHECK_THROWS_AS(process.DecayStructureFor(wrong_spin), std::invalid_argument);

    gra::LORENTZSCALAR malformed;
    malformed.decaytree = {CascadeBranch(23, {StableBranch(-13), StableBranch(13)})};
    CHECK_THROWS_AS(process.DecayStructureFor(malformed), std::invalid_argument);

    for (const auto &tree : {std::vector{StableBranch(2), StableBranch(-2), StableBranch(21)},
                             std::vector{StableBranch(2), StableBranch(-2), StableBranch(21), StableBranch(21)},
                             std::vector{StableBranch(882), StableBranch(-882)}}) {
      gra::LORENTZSCALAR unsupported;
      unsupported.decaytree = tree;
      CHECK_FALSE(process.MatchProcess(tree).has_value());
      CHECK_THROWS_AS(process.DecayStructureFor(unsupported), std::invalid_argument);
    }
  }
}

// Reject a C forbidden production channel even when another channel is allowed
TEST_CASE("Regge initialization rejects a C forbidden production channel", "[process][Regge][C-parity][validation]") {
  for (const std::string model : {"MP", "XP", "GP"}) {
    DYNAMIC_SECTION(model) {
      auto card = nlohmann::json::parse(gra::aux::GetInputData("gencard/test.json"));
      card["SCATTERING"]["PROCESS"] = model + "[RES]<F> -> pi+ pi- @RES{f0_980:1}";
      card["SCATTERING"]["LOOPSCREEN"] = false;
      card["GENERIC"]["HIST"] = 0;
      gra::MGraniitti generator;
      generator.ReadInput(card);
      auto &resonance = generator.proc->state.lts.process.RESONANCES.at("f0_980");
      auto &channels = model == "MP" ? resonance.MP.channels : model == "XP" ? resonance.XP.channels : resonance.GP.channels;
      auto forbidden = channels.front();
      forbidden.exchange[0] = gra::PDG::PDG_gamma;
      channels.push_back(std::move(forbidden));
      CHECK_THROWS_AS(generator.proc->PrepareRun(), std::invalid_argument);
    }
  }
}

// Check Durham continuum processes follow the matrix element and meson domain
TEST_CASE("Durham continuum processes match their amplitude selection", "[process][process][Durham]") {
  const std::vector<gra::MDecayBranch> pion_pair         = {StableBranch(211), StableBranch(-211)};
  const std::vector<gra::MDecayBranch> neutral_pions     = {StableBranch(111), StableBranch(111)};
  const std::vector<gra::MDecayBranch> eta_mixture       = {StableBranch(221), StableBranch(331)};
  const std::vector<gra::MDecayBranch> mixed_meson_pair  = {StableBranch(111), StableBranch(221)};
  const std::vector<gra::MDecayBranch> charged_same_sign = {StableBranch(211), StableBranch(211)};
  const std::vector<gra::MDecayBranch> dimuon            = {StableBranch(-13), StableBranch(13)};

  CHECK(gra::DurhamMesonPairSupported(pion_pair));
  CHECK(gra::DurhamMesonPairSupported(neutral_pions));
  CHECK(gra::DurhamMesonPairSupported(eta_mixture));
  CHECK_FALSE(gra::DurhamMesonPairSupported(mixed_meson_pair));
  CHECK_FALSE(gra::DurhamMesonPairSupported(charged_same_sign));
  CHECK_FALSE(gra::DurhamMesonPairSupported(dimuon));

  const std::vector<gra::MDecayBranch> gluon_pair         = {StableBranch(21), StableBranch(21)};
  const std::vector<gra::MDecayBranch> nested_gluon_roots = {CascadeBranch(21, {StableBranch(2), StableBranch(-2)}),
                                                             CascadeBranch(21, {StableBranch(21), StableBranch(21)})};
  const std::vector<gra::MDecayBranch> gluon_triplet      = {StableBranch(21), StableBranch(21), StableBranch(21)};
  const auto partons = gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::Parton, "gg_QCD");
  const auto mesons = gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::MesonPair, "gg_MM");
  const auto diphotons = gra::MDurham::ProcessDefinitionFor(gra::MDurhamMode::PhotonPair, "gg_yy");
  CHECK(partons->MatchProcess(gluon_pair).has_value());
  CHECK(partons->MatchProcess(gluon_triplet).has_value());
  CHECK_FALSE(partons->MatchProcess(nested_gluon_roots).has_value());
  CHECK_FALSE(partons->MatchProcess(pion_pair).has_value());
  CHECK_FALSE(partons->MatchProcess({StableBranch(22), StableBranch(22)}).has_value());
  CHECK(mesons->MatchProcess(pion_pair).has_value());
  CHECK_FALSE(mesons->MatchProcess(gluon_pair).has_value());
  CHECK(diphotons->MatchProcess({StableBranch(22), StableBranch(22)}).has_value());
  CHECK_FALSE(diphotons->MatchProcess(pion_pair).has_value());
  const auto continuum = gra::MDurham::ContinuumProcesses();

  const auto photons = continuum->MatchProcess({StableBranch(22), StableBranch(22)});
  REQUIRE(photons.has_value());
  CHECK(photons->process_name == "gg_yy");
  CHECK(photons->matrix_element_form == gra::amplitude::MatrixElementForm::Analytic);
  CHECK_FALSE(gra::HasDurhamMG5Process(*photons));
  CHECK(gra::CreateDurhamMG5Process(*photons) == nullptr);
  CHECK_FALSE(gra::amplitude::FindProcess("DURHAM", "gg_yy").has_value());

  const auto generated_pair = continuum->MatchProcess(gluon_pair);
  REQUIRE(generated_pair.has_value());
  CHECK(generated_pair->process_family == "DURHAM");
  CHECK(generated_pair->process_name == "gg_gg");
  CHECK(generated_pair->matrix_element_form == gra::amplitude::MatrixElementForm::Generated);
  CHECK(generated_pair->topology_mode == gra::amplitude::TopologyMode::Exact);

  const auto generated_triplet = continuum->MatchProcess(gluon_triplet);
  REQUIRE(generated_triplet.has_value());
  CHECK(generated_triplet->process_family == "DURHAM");
  CHECK(generated_triplet->process_name == "gg_ggg");
  CHECK(continuum->MatchProcess(pion_pair).has_value());
  CHECK_FALSE(continuum->MatchProcess(dimuon).has_value());
}

// Initialize algorithmically selected vector cascades through the real model readers
TEST_CASE("Model card examples initialize physical vector cascades", "[process][examples][cascade]") {
  const gra::MProcessTable table("TUNE0");
  const gra::MSubProc subprocess;
  for (const auto &process : subprocess.CreateAllProcesses()) {
    if (process->ISTATE != "MP" && process->ISTATE != "XP" && process->ISTATE != "GP") { continue; }
    if (process->CHANNEL != "CON" && process->CHANNEL != "RES") { continue; }
    const auto examples = table.Examples(*process);
    const auto nested = std::find_if(examples.begin(), examples.end(), [](const auto &value) {
      return value.find(" > {") != std::string::npos;
    });
    REQUIRE(nested != examples.end());
    DYNAMIC_SECTION(process->ISTATE << " " << process->CHANNEL << " " << *nested) {
      auto card = nlohmann::json::parse(gra::aux::GetInputData("gencard/test.json"));
      card["SCATTERING"]["PROCESS"] = process->ISTATE + "[" + process->CHANNEL + "]<F> -> " + *nested;
      card["SCATTERING"]["RES"] = nlohmann::json::array();
      card["SCATTERING"]["LOOPSCREEN"] = false;
      card["GENERIC"]["HIST"] = 0;
      gra::MGraniitti generator;
      REQUIRE_NOTHROW(generator.ReadInput(card));
      REQUIRE_NOTHROW(generator.proc->PrepareRun());
      CHECK(generator.proc->ProcPtr.MatchProcess(generator.proc->state.lts.decaytree).has_value());
    }
  }
}

// Distinguish generic scalar helicity decays from complete generated decay amplitudes
TEST_CASE("scalar photon processes declare physical Jacob-Wick cascades",
          "[process][photon][cascade][proposal]") {
  gra::MSubProc process;
  gra::LORENTZSCALAR lts;
  lts.decaytree = {CascadeBranch(23, {StableBranch(13), StableBranch(-13)}),
                   CascadeBranch(23, {StableBranch(11), StableBranch(-11)})};
  for (const std::string channel : {"Higgs", "monopolium(0)"}) {
    CAPTURE(channel);
    process.Initialize("yy", channel);
    lts.process.root_decay_mode = gra::RootDecayMode::Physical;
    for (const bool coherent : {false, true}) {
      lts.amplitude.DECAY_SYM = coherent;
      const auto physical = process.DecayStructureFor(lts);
      REQUIRE(physical.type == (coherent ? gra::DecayType::JacobWickCoherent : gra::DecayType::JacobWickIncoherent));
    }
    lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
    REQUIRE(process.DecayStructureFor(lts).type == gra::DecayType::None);
  }
}

// Apply the same decay ownership checks to every analytic process definition
TEST_CASE("analytic decay ownership enforces isolated production support", "[process][analytic][cascade][proposal]") {
  gra::LORENTZSCALAR lts;
  lts.decaytree = {CascadeBranch(23, {StableBranch(13), StableBranch(-13)}), StableBranch(22)};
  for (const auto type : {gra::DecayType::None, gra::DecayType::JacobWickIncoherent, gra::DecayType::JacobWickCoherent, gra::DecayType::Full}) {
    const gra::DecayStructure structure{type};
    gra::amplitude::AnalyticProcess process("TEST", "nested", "nested decay", structure,
        [](const auto &tree) { return !tree.empty(); }, [structure](const auto &) { return structure; });
    lts.process.root_decay_mode = gra::RootDecayMode::Physical;
    REQUIRE(process.DecayStructureFor(lts) == structure);
    lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
    if (type == gra::DecayType::Full) {
      REQUIRE_THROWS_AS(process.DecayStructureFor(lts), std::invalid_argument);
    } else {
      REQUIRE(process.DecayStructureFor(lts) == structure);
    }
  }
}

// Require generated matrix elements to declare their complete stable-particle amplitude
TEST_CASE("All generated processes require complete decay amplitudes", "[process][MG5][decay][proposal]") {
  const auto processes = gra::amplitude::AllProcesses();
  REQUIRE_FALSE(processes.empty());
  for (const auto &process : processes) {
    INFO(process.process_family << " " << process.process_name);
    REQUIRE(process.decay_structure.type == gra::DecayType::Full);
    REQUIRE_NOTHROW(gra::amplitude::ProcessRegistry({process}));
    for (const auto type : {gra::DecayType::None, gra::DecayType::JacobWickIncoherent, gra::DecayType::JacobWickCoherent}) {
      auto wrong = process;
      wrong.decay_structure.type = type;
      REQUIRE_THROWS_AS(gra::amplitude::ProcessRegistry({wrong}), std::invalid_argument);
    }
  }
}

// Check the same declaration for future generated chains at arbitrary depth
TEST_CASE("New generated cascade topologies inherit the full amplitude sampling requirement",
          "[process][MG5][decay][proposal]") {
  std::vector<gra::MDecayBranch> tree = {StableBranch(-13), StableBranch(13)};
  for (int depth = 1; depth <= 4; ++depth) {
    tree = {CascadeBranch(900000 + depth, tree), StableBranch(22)};
    const auto process = GeneratedProcess("nested", "nested", gra::amplitude::AmplitudeTopologyFromDecayTree(tree),
                                          gra::amplitude::StableFinalStatePDGs(tree), gra::amplitude::TopologyMode::Exact);
    const gra::amplitude::ProcessRegistry registry({process});
    gra::LORENTZSCALAR lts;
    lts.decaytree = tree;
    const auto decay = registry.DecayStructureFor(lts);
    REQUIRE(decay.type == gra::DecayType::Full);
    auto wrong = process;
    wrong.decay_structure.type = gra::DecayType::JacobWickIncoherent;
    REQUIRE_THROWS_AS(gra::amplitude::ProcessRegistry({wrong}), std::invalid_argument);
  }
}

// Check generated decay structures require the exact registered topology
TEST_CASE("Generated amplitudes require exact decay topology", "[process][MG5][decay]") {
  const gra::ProcessDescriptor hard_description{"generated", "MG5", "IPp", "", 1};
  gra::MGeneratedPartonProc    pp_z("IPp", "Z", hard_description, "MG5_PP_Z");
  gra::LORENTZSCALAR           direct_z;
  direct_z.decaytree      = {StableBranch(-13), StableBranch(13)};
  const auto direct_decay = pp_z.DecayStructureFor(direct_z);
  CHECK(direct_decay.type == gra::DecayType::Full);
  CHECK_FALSE(direct_decay.allows_isolated_resonance);
  direct_z.process.root_decay_mode = gra::RootDecayMode::Isolated;
  CHECK_THROWS_AS(pp_z.DecayStructureFor(direct_z), std::invalid_argument);
  direct_z.process.root_decay_mode = gra::RootDecayMode::None;

  gra::MGeneratedPartonProc pp_zj("IPp", "Zj", hard_description, "MG5_PP_ZJ");
  gra::LORENTZSCALAR        cascaded_zj;
  cascaded_zj.decaytree = {
      CascadeBranch(23, {StableBranch(-13), StableBranch(13)}),
      StableBranch(gra::PDG::PDG_hard_jet),
  };
  const auto cascaded_decay = pp_zj.DecayStructureFor(cascaded_zj);
  CHECK(cascaded_decay.type == gra::DecayType::Full);

  gra::LORENTZSCALAR wrong_zj;
  wrong_zj.decaytree = {
      StableBranch(-13),
      StableBranch(13),
      StableBranch(gra::PDG::PDG_hard_jet),
  };
  CHECK_THROWS_AS(pp_zj.DecayStructureFor(wrong_zj), std::invalid_argument);

  gra::MGeneratedPhotonProc photon("yy", "WW", {"generated", "MG5", "yy", "", 1}, "MG5_YY_WW");
  gra::LORENTZSCALAR        cascaded_ww;
  cascaded_ww.decaytree = {
      CascadeBranch(24, {StableBranch(-11), StableBranch(12)}),
      CascadeBranch(-24, {StableBranch(11), StableBranch(-12)}),
  };
  const auto ww_process = photon.MatchProcess(cascaded_ww.decaytree);
  REQUIRE(ww_process.has_value());
  CHECK(ww_process->process_family == "MG5_YY_WW");
  CHECK(ww_process->decay_structure.type == gra::DecayType::Full);

  gra::LORENTZSCALAR direct_ww_leaves;
  direct_ww_leaves.decaytree = {
      StableBranch(-11),
      StableBranch(12),
      StableBranch(11),
      StableBranch(-12),
  };
  CHECK_FALSE(photon.MatchProcess(direct_ww_leaves.decaytree).has_value());

  auto ww_matrix_element = gra::CreatePhotonMG5Process("MG5_YY_WW");
  REQUIRE(ww_matrix_element != nullptr);
  const auto wrong_evaluation = ww_matrix_element->Evaluate(direct_ww_leaves, 0.118, false);
  CHECK_FALSE(wrong_evaluation.Valid());
  CHECK(direct_ww_leaves.hamp.empty());
  CHECK(wrong_evaluation.amp2 == Approx(0.0));
}

// Check failed process matching clears the previous event decay structure
TEST_CASE("Failed amplitude resolution clears stale sampling state", "[process][decay]") {
  gra::MSubProc subprocess;
  subprocess.Initialize("IPp", "Z");
  gra::LORENTZSCALAR lts             = DirectZKinematics();
  const auto         decay_structure = subprocess.DecayStructureFor(lts);
  REQUIRE((decay_structure.type == gra::DecayType::Full));
  REQUIRE_FALSE(decay_structure.allows_isolated_resonance);

  lts.decaytree = {StableBranch(22)};
  CHECK_THROWS_AS(subprocess.DecayStructureFor(lts), std::invalid_argument);
  CHECK(lts.decay_structure == gra::DecayStructure{});
}

// Check every standalone factory owns one exact immutable process
TEST_CASE("Standalone MG5 matrix elements own exact generated processes", "[process][MG5][matrix-element]") {
  const gra::MPDG pdg = ProcessPDGTable();

  const auto photon_processes = gra::amplitude::Processes("PHOTON");
  REQUIRE_FALSE(photon_processes.empty());
  for (const auto &process : photon_processes) {
    INFO("Photon matrix element " << process.process_name);
    REQUIRE(gra::HasPhotonMG5Process(process));
    auto matrix_element = gra::CreatePhotonMG5Process(process);
    REQUIRE(matrix_element != nullptr);
    REQUIRE(matrix_element->Processes() == std::vector<gra::amplitude::Process>{process});
    REQUIRE(matrix_element->Processes().at(0).decay_structure == process.decay_structure);
    REQUIRE(matrix_element->SubprocessCount() == 1);

    const auto tree = ProcessTree(pdg, process);
    REQUIRE(matrix_element->MatchProcess(tree) == process);
    gra::LORENTZSCALAR lts;
    lts.decaytree = tree;
    REQUIRE(matrix_element->DecayStructureFor(lts) == process.decay_structure);

    lts.decaytree = MismatchedTree(tree);
    CHECK_THROWS_AS(matrix_element->DecayStructureFor(lts), std::invalid_argument);
    lts.hamp          = {1.0};
    const auto result = matrix_element->Evaluate(lts, 0.118, false);
    CHECK_FALSE(result.Valid());
    CHECK(lts.hamp.empty());
    CHECK(lts.hamp.empty());
    CHECK(lts.hard_color_flows.empty());

    auto forged = process;
    forged.topology_signature += ",forged";
    CHECK_FALSE(gra::HasPhotonMG5Process(forged));
    CHECK(gra::CreatePhotonMG5Process(forged) == nullptr);
  }

  const auto durham_processes = gra::amplitude::Processes("DURHAM");
  REQUIRE_FALSE(durham_processes.empty());
  for (const auto &process : durham_processes) {
    INFO("Durham matrix element " << process.process_name);
    REQUIRE(gra::HasDurhamMG5Process(process));
    auto matrix_element = gra::CreateDurhamMG5Process(process);
    REQUIRE(matrix_element != nullptr);
    REQUIRE(matrix_element->Name() == process.process_name);
    REQUIRE(matrix_element->FinalPDGs() == process.stable_pdgs);
    REQUIRE(matrix_element->Processes() == std::vector<gra::amplitude::Process>{process});
    REQUIRE(matrix_element->Processes().at(0).decay_structure == process.decay_structure);

    const auto tree = ProcessTree(pdg, process);
    REQUIRE(matrix_element->MatchProcess(tree) == process);
    gra::LORENTZSCALAR lts;
    lts.decaytree = tree;
    REQUIRE(matrix_element->DecayStructureFor(lts) == process.decay_structure);

    lts.decaytree = MismatchedTree(tree);
    CHECK_THROWS_AS(matrix_element->DecayStructureFor(lts), std::invalid_argument);
    const auto result = matrix_element->Evaluate(lts, 0.118);
    CHECK_FALSE(result.Valid());
    CHECK(result.projected.empty());
    CHECK(result.flow_projected.empty());

    auto forged = process;
    forged.process_syntax += " forged";
    CHECK_FALSE(gra::HasDurhamMG5Process(forged));
    CHECK(gra::CreateDurhamMG5Process(forged) == nullptr);
  }
}

// Check hard process families construct independent subprocess sums
TEST_CASE("Hard MG5 amplitudes are selected by process family", "[process][MG5][hard-router]") {
  const std::vector<std::string> expected = {"MG5_PP_Z", "MG5_PP_ZJ", "MG5_PP_JJ", "MG5_PP_W"};
  std::vector<std::string>       process_families;
  for (const auto &info : gra::PartonMG5ProcessInfos()) { process_families.push_back(info.process_family); }
  REQUIRE(process_families == expected);
  CHECK(std::set<std::string>(process_families.begin(), process_families.end()).size() == process_families.size());

  for (const auto &process_family : process_families) {
    INFO("Hard process family " << process_family);
    REQUIRE_FALSE(gra::amplitude::Processes(process_family).empty());
    auto first  = gra::CreatePartonMG5Process(process_family);
    auto second = gra::CreatePartonMG5Process(process_family);
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    CHECK(first.get() != second.get());
    CHECK(first->SubprocessCount() > 0);
    CHECK(second->SubprocessCount() == first->SubprocessCount());
  }

  CHECK(gra::CreatePartonMG5Process("MG5_YY_ZJJ") == nullptr);
  CHECK(gra::CreatePartonMG5Process("UNKNOWN") == nullptr);
}

// Check photon process families construct independent subprocess sums
TEST_CASE("Photon MG5 families are selected by process family", "[process][MG5][photon-router]") {
  const auto infos = gra::PhotonMG5ProcessInfos();
  REQUIRE(infos.size() == 3);
  for (const auto &info : infos) {
    INFO("Photon process family " << info.process_family);
    auto first  = gra::CreatePhotonMG5Process(info.process_family);
    auto second = gra::CreatePhotonMG5Process(info.process_family);
    REQUIRE(first != nullptr);
    REQUIRE(second != nullptr);
    CHECK(first.get() != second.get());
    CHECK(first->SubprocessCount() > 0);
    CHECK(second->SubprocessCount() == first->SubprocessCount());
  }
  CHECK(gra::CreatePhotonMG5Process("MG5_PP_JJ") == nullptr);
  CHECK(gra::CreatePhotonMG5Process("UNKNOWN") == nullptr);
}

// Check every GRANIITTI route follows its amplitude process definition
TEST_CASE("MG5 process routes cannot diverge from amplitude definitions", "[process][MG5][delegation]") {
  const gra::MPDG pdg = ProcessPDGTable();

  gra::amplitude::ProcessRegistry photon_registry(gra::amplitude::Processes("PHOTON"));
  gra::PROC_003_QED_YY_EPA        photon_kt;
  gra::PROC_030_QED_YY_DZ_EPA     photon_dz;
  gra::PROC_040_QED_YY_LUX_EPA    photon_lux;
  CheckProcessRoute(photon_registry, photon_kt, pdg, false);
  CheckProcessRoute(photon_registry, photon_dz, pdg, false);
  CheckProcessRoute(photon_registry, photon_lux, pdg, false);

  gra::amplitude::ProcessRegistry durham_registry(gra::amplitude::Processes("DURHAM"));
  gra::MDurhamContinuumProc durham("QCD", gra::MDurhamMode::Parton, 703);
  CheckProcessRoute(durham_registry, durham, pdg, true);

  const gra::ProcessDescriptor photon_description{"generated", "MG5", "yy", "", 1};
  gra::AMP_MG5_yy_zjj          yy_zjj_amplitude;
  gra::MGeneratedPhotonProc    yy_zjj_kt("yy", "Zjj", photon_description, "MG5_YY_ZJJ");
  gra::MGeneratedPhotonProc    yy_zjj_dz("yy_DZ", "Zjj", photon_description, "MG5_YY_ZJJ");
  gra::MGeneratedPhotonProc    yy_zjj_lux("yy_LUX", "Zjj", photon_description, "MG5_YY_ZJJ");
  CheckProcessRoute(yy_zjj_amplitude, yy_zjj_kt, pdg, true);
  CheckProcessRoute(yy_zjj_amplitude, yy_zjj_dz, pdg, true);
  CheckProcessRoute(yy_zjj_amplitude, yy_zjj_lux, pdg, true);

  gra::AMP_MG5_yy_jj yy_jj_amplitude;
  for (const std::string initial : {"yy", "yy_DZ", "yy_LUX"}) {
    gra::MGeneratedPhotonProc yy_jj(initial, "jj", photon_description, "MG5_YY_JJ");
    CheckProcessRoute(yy_jj_amplitude, yy_jj, pdg, true);
  }

  gra::AMP_MG5_yy_ww yy_ww_amplitude;
  for (const std::string initial : {"yy", "yy_DZ", "yy_LUX"}) {
    gra::MGeneratedPhotonProc yy_ww(initial, "WW", photon_description, "MG5_YY_WW");
    CheckProcessRoute(yy_ww_amplitude, yy_ww, pdg, true);
  }

  const gra::ProcessDescriptor hard_description{"generated", "MG5", "IPp", "", 2};
  gra::AMP_MG5_pp_z            pp_z_amplitude;
  gra::MGeneratedPartonProc    pp_z_ip("IPp", "Z", hard_description, "MG5_PP_Z");
  gra::MGeneratedPartonProc    pp_z_ipip("IPIP", "Z", hard_description, "MG5_PP_Z");
  CheckProcessRoute(pp_z_amplitude, pp_z_ip, pdg, true);
  CheckProcessRoute(pp_z_amplitude, pp_z_ipip, pdg, true);

  gra::AMP_MG5_pp_zj        pp_zj_amplitude;
  gra::MGeneratedPartonProc pp_zj("IPp", "Zj", hard_description, "MG5_PP_ZJ");
  CheckProcessRoute(pp_zj_amplitude, pp_zj, pdg, true);

  gra::AMP_MG5_pp_jj        pp_jj_amplitude;
  gra::MGeneratedPartonProc pp_jj("IPIP", "jj", hard_description, "MG5_PP_JJ");
  CheckProcessRoute(pp_jj_amplitude, pp_jj, pdg, true);

  gra::AMP_MG5_pp_w         pp_w_amplitude;
  gra::MGeneratedPartonProc pp_w("IPp", "W", hard_description, "MG5_PP_W");
  CheckProcessRoute(pp_w_amplitude, pp_w, pdg, true);

  const auto zj_process = pp_zj_amplitude.Processes().at(0);
  auto       top_tree   = ProcessTree(pdg, zj_process);
  REQUIRE(top_tree.size() == 2);
  top_tree[1].p.pdg  = 6;
  top_tree[1].p.name = "t";
  CHECK_FALSE(pp_zj_amplitude.MatchProcess(top_tree).has_value());
  CHECK_FALSE(pp_zj.MatchProcess(top_tree).has_value());

  auto forged_name_tree = ProcessTree(pdg, zj_process);
  REQUIRE_FALSE(forged_name_tree.at(0).legs.empty());
  forged_name_tree.at(0).legs.at(0).p.pdg  = 9000013;
  forged_name_tree.at(0).legs.at(0).p.name = "mu+";
  CHECK_FALSE(pp_zj_amplitude.MatchProcess(forged_name_tree).has_value());
  CHECK_FALSE(pp_zj.MatchProcess(forged_name_tree).has_value());

  auto reversed_decay = ProcessTree(pdg, zj_process);
  REQUIRE(reversed_decay.at(0).legs.size() == 2);
  std::swap(reversed_decay.at(0).legs.at(0), reversed_decay.at(0).legs.at(1));
  CHECK_FALSE(pp_zj_amplitude.MatchProcess(reversed_decay).has_value());
  CHECK_FALSE(pp_zj.MatchProcess(reversed_decay).has_value());
}

// Check every registered generated process has a matrix element in MSubProc
TEST_CASE("Registered MG5 processes have matrix elements", "[process][MG5][registry]") {
  using ProcessKey = std::pair<std::string, std::string>;
  std::set<ProcessKey> advertised;
  for (const auto &process : gra::amplitude::AllProcesses()) {
    advertised.emplace(process.process_family, process.process_name);
  }
  for (const auto &process : gra::LightByLightProcesses()) {
    advertised.emplace(process.process_family, process.process_name);
  }

  std::set<ProcessKey> routed;
  const gra::MSubProc  subprocess;
  for (const auto &process : subprocess.CreateAllProcesses()) {
    for (const auto &process : process->Processes()) {
      if (gra::amplitude::IsGenerated(process.matrix_element_form)) {
        routed.emplace(process.process_family, process.process_name);
      }
    }
  }
  REQUIRE(routed.size() == advertised.size());
  for (const auto &key : advertised) {
    INFO("Registered MG5 process " << key.first << ":" << key.second);
    CHECK(routed.count(key) == 1);
  }
}

// Check duplicate event topology claims are rejected independent of names
TEST_CASE("Amplitude process definitions reject ambiguous topologies", "[process][MG5][topology]") {
  const gra::amplitude::Process first =
      GeneratedProcess("first", "u,u~", {{{2}, {}}, {{-2}, {}}}, {2, -2}, gra::amplitude::TopologyMode::Exact);
  auto second                                          = first;
  second.process_name                                  = "second";
  second.process_syntax                                = "a a > u u~";
  const std::vector<gra::amplitude::Process> ambiguous = {first, second};
  CHECK_THROWS_AS(gra::amplitude::ProcessRegistry(ambiguous), std::invalid_argument);
}

// Check generated selectors match independently at every topology node
TEST_CASE("Amplitude process selectors are local to each topology node", "[process][MG5][topology]") {
  const std::vector<int> partons = {-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89};
  const auto             mixed =
      GeneratedProcess("mixed", "g,j", {{{21}, {}}, {partons, {}}}, {21, 0}, gra::amplitude::TopologyMode::FlavourSet);
  const auto nested = GeneratedProcess("nested", "z(mu+,j),g", {{{23}, {{{-13}, {}}, {partons, {}}}}, {{21}, {}}},
                                       {-13, 0, 21}, gra::amplitude::TopologyMode::FlavourSet);
  gra::amplitude::ProcessRegistry definition({mixed, nested});

  const auto mixed_match = definition.MatchProcess({StableBranch(21), StableBranch(2)});
  REQUIRE(mixed_match.has_value());
  CHECK(mixed_match->process_name == "mixed");
  CHECK_FALSE(definition.MatchProcess({StableBranch(2), StableBranch(21)}).has_value());

  const std::vector<gra::MDecayBranch> nested_tree = {
      CascadeBranch(23, {StableBranch(-13), StableBranch(-1)}),
      StableBranch(21),
  };
  const auto nested_match = definition.MatchProcess(nested_tree);
  REQUIRE(nested_match.has_value());
  CHECK(nested_match->process_name == "nested");

  gra::LORENTZSCALAR lts;
  lts.decaytree       = nested_tree;
  const auto resolved = definition.ResolveProcess(lts);
  REQUIRE(resolved.has_value());
  CHECK(resolved->matrix_element_form == gra::amplitude::MatrixElementForm::Generated);
  CHECK(resolved->topology_mode == gra::amplitude::TopologyMode::Event);
  CHECK(resolved->stable_pdgs == std::vector<int>{-13, -1, 21});
  CHECK(resolved->topology == gra::amplitude::AmplitudeTopologyFromDecayTree(nested_tree));
  REQUIRE(resolved->channel_topologies.size() == 1);
  CHECK(resolved->channel_topologies.front() == resolved->topology);
}

// Check generated channels constrain event states but accept aliases
TEST_CASE("Amplitude process patterns retain exact generated correlations", "[process][MG5][topology]") {
  const std::vector<int> partons     = {-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89};
  auto                   constrained = GeneratedProcess("constrained", "j,j", {{partons, {}}, {partons, {}}}, {0, 0},
                                                        gra::amplitude::TopologyMode::FlavourSet);
  constrained.channel_topologies     = {
          {{{21}, {}}, {{2}, {}}},
          {{{2}, {}}, {{2}, {}}},
  };
  gra::amplitude::ProcessRegistry definition({constrained});

  CHECK(definition.MatchProcess({StableBranch(gra::PDG::PDG_hard_jet), StableBranch(gra::PDG::PDG_hard_jet)})
            .has_value());
  CHECK(definition.MatchProcess({StableBranch(gra::PDG::PDG_hard_jet), StableBranch(2)}).has_value());
  CHECK_FALSE(definition.MatchProcess({StableBranch(gra::PDG::PDG_hard_jet), StableBranch(21)}).has_value());
  CHECK_FALSE(definition.MatchProcess({StableBranch(1), StableBranch(gra::PDG::PDG_hard_jet)}).has_value());
  CHECK(definition.MatchProcess({StableBranch(21), StableBranch(2)}).has_value());
  CHECK(definition.MatchProcess({StableBranch(2), StableBranch(2)}).has_value());
  CHECK_FALSE(definition.MatchProcess({StableBranch(2), StableBranch(21)}).has_value());

  gra::LORENTZSCALAR lts;
  lts.decaytree       = {StableBranch(21), StableBranch(2)};
  const auto resolved = definition.ResolveProcess(lts);
  REQUIRE(resolved.has_value());
  CHECK(resolved->channel_topologies == std::vector<gra::AmplitudeTopology>{resolved->topology});

  auto unresolved               = constrained;
  unresolved.channel_topologies = {{{{89}, {}}, {{2}, {}}}};
  CHECK_THROWS_AS(gra::amplitude::ProcessRegistry({unresolved}), std::invalid_argument);

  auto outside               = constrained;
  outside.channel_topologies = {{{{22}, {}}, {{2}, {}}}};
  CHECK_THROWS_AS(gra::amplitude::ProcessRegistry({outside}), std::invalid_argument);
}

// Check narrower generated selectors win and incomparable overlaps are invalid
TEST_CASE("Amplitude process selector overlap has deterministic semantics", "[process][MG5][topology]") {
  const std::vector<int> partons   = {-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89};
  const auto             inclusive = GeneratedProcess("inclusive", "j,j", {{partons, {}}, {partons, {}}}, {0, 0},
                                                      gra::amplitude::TopologyMode::FlavourSet);
  const auto             exact =
      GeneratedProcess("exact", "g,g", {{{21}, {}}, {{21}, {}}}, {21, 21}, gra::amplitude::TopologyMode::Exact);
  gra::amplitude::ProcessRegistry definition({inclusive, exact});

  const auto exact_match = definition.MatchProcess({StableBranch(21), StableBranch(21)});
  REQUIRE(exact_match.has_value());
  CHECK(exact_match->process_name == "exact");
  const auto inclusive_match = definition.MatchProcess({StableBranch(21), StableBranch(2)});
  REQUIRE(inclusive_match.has_value());
  CHECK(inclusive_match->process_name == "inclusive");

  const auto left =
      GeneratedProcess("left", "g,j", {{{21}, {}}, {partons, {}}}, {21, 0}, gra::amplitude::TopologyMode::FlavourSet);
  const auto right =
      GeneratedProcess("right", "j,g", {{partons, {}}, {{21}, {}}}, {0, 21}, gra::amplitude::TopologyMode::FlavourSet);
  CHECK_THROWS_AS(gra::amplitude::ProcessRegistry({left, right}), std::invalid_argument);
}

// Check numeric internal selectors do not depend on mutable particle names
TEST_CASE("Amplitude topology identity uses internal PDG selectors", "[process][MG5][topology]") {
  const auto first  = GeneratedProcess("first_bsm", "x(mu+,mu-)", {{{9000001}, {{{-13}, {}}, {{13}, {}}}}}, {-13, 13},
                                       gra::amplitude::TopologyMode::Exact);
  const auto second = GeneratedProcess("second_bsm", "x(mu+,mu-)", {{{9000002}, {{{-13}, {}}, {{13}, {}}}}}, {-13, 13},
                                       gra::amplitude::TopologyMode::Exact);
  gra::amplitude::ProcessRegistry definition({first, second});

  auto first_tree       = std::vector<gra::MDecayBranch>{CascadeBranch(9000001, {StableBranch(-13), StableBranch(13)})};
  auto second_tree      = first_tree;
  first_tree[0].p.name  = "same name";
  second_tree[0].p.pdg  = 9000002;
  second_tree[0].p.name = "same name";
  REQUIRE(definition.MatchProcess(first_tree).has_value());
  REQUIRE(definition.MatchProcess(second_tree).has_value());
  CHECK(definition.MatchProcess(first_tree)->process_name == "first_bsm");
  CHECK(definition.MatchProcess(second_tree)->process_name == "second_bsm");

  auto renamed_first_tree      = first_tree;
  renamed_first_tree[0].p.name = "different mutable name";
  CHECK(gra::amplitude::AmplitudeTopologySignature(renamed_first_tree) ==
        gra::amplitude::AmplitudeTopologySignature(first_tree));
  CHECK(gra::amplitude::AmplitudeTopologySignature(first_tree) == "pdg:9000001(mu+,mu-)");
}

// Check analytic process definitions reject empty and isolated trees
TEST_CASE("Analytic processes share topology invariants", "[process][topology]") {
  gra::amplitude::AnalyticProcess definition(
      "TEST", "analytic", "any nonempty tree", {gra::DecayType::Full},
      [](const std::vector<gra::MDecayBranch> &) { return true; },
      [](const gra::LORENTZSCALAR &) {
        return gra::DecayStructure{gra::DecayType::Full};
      });
  CHECK_FALSE(definition.MatchProcess({}).has_value());

  gra::LORENTZSCALAR lts;
  lts.decaytree               = {CascadeBranch(113, {StableBranch(211), StableBranch(-211)})};
  lts.process.root_decay_mode = gra::RootDecayMode::Isolated;
  CHECK_THROWS_AS(definition.DecayStructureFor(lts), std::invalid_argument);
}

// Check converter-generated yy Zjj channels match the registered topologies
TEST_CASE("Gamma-gamma Zjj subprocess channels follow generated processes", "[process][MG5][channel]") {
  const gra::MPDG pdg                 = ProcessPDGTable();
  const auto      processes           = gra::amplitude::Processes("MG5_YY_ZJJ");
  const auto      subprocess_channels = MG5_YY_ZJJ::SubprocessChannels();
  REQUIRE(subprocess_channels.size() == processes.size());
  gra::AMP_MG5_yy_zjj amplitude;
  CHECK(amplitude.Processes() == processes);

  const std::array<int, 2> photons       = {22, 22};
  std::size_t              channel_count = 0;
  for (const auto &channels : subprocess_channels) { channel_count += channels.size(); }
  REQUIRE(channel_count == processes.size());

  for (const auto &process : processes) {
    INFO("Generated yy Zjj channel " << process.process_name);
    const auto tree = ProcessTree(pdg, process);
    REQUIRE(tree.size() == 3);
    REQUIRE(tree[0].p.pdg == 23);
    REQUIRE(tree[0].legs.size() == 2);
    const std::vector<int> expected_final    = gra::amplitude::StableFinalStatePDGs(tree);
    const auto             expected_topology = gra::amplitude::AmplitudeTopologyFromDecayTree(tree);

    std::size_t matches = 0;
    for (const auto &channels : subprocess_channels) {
      matches += static_cast<std::size_t>(
          std::count_if(channels.begin(), channels.end(), [&](const gra::mg5::Channel &channel) {
            return channel.initial == photons && channel.final == expected_final &&
                   channel.topology == expected_topology;
          }));
    }
    CHECK(matches == 1);
  }
}

// Check subprocess channels distinguish resonance trees with equal stable
// leaves
TEST_CASE("MG5 subprocess channels match complete numeric topologies", "[process][MG5][topology]") {
  const auto z_tree = std::vector<gra::MDecayBranch>{CascadeBranch(23, {StableBranch(-13), StableBranch(13)})};
  const auto h_tree = std::vector<gra::MDecayBranch>{CascadeBranch(25, {StableBranch(-13), StableBranch(13)})};
  const std::vector<int> stable_pdgs = gra::amplitude::StableFinalStatePDGs(z_tree);
  REQUIRE(stable_pdgs == gra::amplitude::StableFinalStatePDGs(h_tree));

  const gra::mg5::Channel      z_channel{{21, 21}, stable_pdgs, gra::amplitude::AmplitudeTopologyFromDecayTree(z_tree)};
  const gra::mg5::Channel      h_channel{{21, 21}, stable_pdgs, gra::amplitude::AmplitudeTopologyFromDecayTree(h_tree)};
  const gra::AmplitudeTopology event_topology = gra::amplitude::AmplitudeTopologyFromDecayTree(z_tree);

  CHECK(gra::mg5::FinalStateMatches(z_channel.final, stable_pdgs));
  CHECK(gra::mg5::FinalStateMatches(h_channel.final, stable_pdgs));
  CHECK(gra::mg5::ChannelMatches(z_channel, event_topology, stable_pdgs));
  CHECK_FALSE(gra::mg5::ChannelMatches(h_channel, event_topology, stable_pdgs));
}

// Check the generated subprocess sum rejects unsupported incoming states
TEST_CASE("MG5 channel selection accepts only supported massless incoming particles", "[process][MG5][channel]") {
  for (const int pdg : {-16, -14, -13, -12, -11, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 11, 12, 13, 14, 16, 21, 22}) {
    CHECK(gra::mg5::IsSupportedInitialPDG(pdg));
  }
  for (const int pdg : {-24, -21, -17, -15, -10, -6, 0, 6, 10, 15, 17, 23, 24, 25, gra::PDG::PDG_hard_jet}) {
    CHECK_FALSE(gra::mg5::IsSupportedInitialPDG(pdg));
  }

  std::vector<gra::mg5::Subprocess<SubprocessTestMatrixElement>> photon_channels;
  photon_channels.push_back(
      {std::make_unique<SubprocessTestMatrixElement>(),
       {{{22, 22}, {-13, 13}, gra::amplitude::AmplitudeTopologyFromDecayTree({StableBranch(-13), StableBranch(13)})}}});
  CHECK_NOTHROW(gra::mg5::ValidateSubprocessInitialStates(photon_channels, gra::mg5::PhotonInitialStates(),
                                                          "photon domain test"));
  CHECK_THROWS_AS(gra::mg5::ValidateSubprocessInitialStates(photon_channels, gra::mg5::MasslessQCDInitialStates(),
                                                            "QCD domain test"),
                  std::invalid_argument);

  std::vector<gra::mg5::Subprocess<SubprocessTestMatrixElement>> qcd_channels;
  qcd_channels.push_back(
      {std::make_unique<SubprocessTestMatrixElement>(),
       {{{21, 2}, {-13, 13}, gra::amplitude::AmplitudeTopologyFromDecayTree({StableBranch(-13), StableBranch(13)})}}});
  CHECK_NOTHROW(
      gra::mg5::ValidateSubprocessInitialStates(qcd_channels, gra::mg5::MasslessQCDInitialStates(), "QCD domain test"));
  CHECK_THROWS_AS(
      gra::mg5::ValidateSubprocessInitialStates(qcd_channels, gra::mg5::PhotonInitialStates(), "photon domain test"),
      std::invalid_argument);

  std::vector<gra::mg5::Subprocess<SubprocessTestMatrixElement>> mixed_channels;
  mixed_channels.push_back(
      {std::make_unique<SubprocessTestMatrixElement>(),
       {{{22, 2}, {-13, 13}, gra::amplitude::AmplitudeTopologyFromDecayTree({StableBranch(-13), StableBranch(13)})}}});
  const gra::mg5::InitialStates mixed_domain = {{{22, 2}}};
  CHECK_NOTHROW(gra::mg5::ValidateSubprocessInitialStates(mixed_channels, mixed_domain, "mixed domain test"));
  CHECK_THROWS_AS(gra::mg5::ValidateSubprocessInitialStates(mixed_channels, gra::mg5::InitialStates{{{2, 22}}},
                                                            "reversed mixed domain test"),
                  std::invalid_argument);
}

// Check the generated subprocess sum rejects nonphysical alpha_s inputs
TEST_CASE("MG5 subprocess sum validates the effective alpha_s", "[process][MG5][coupling]") {
  gra::LORENTZSCALAR lts;
  double             effective_alpha_s = -1.0;
  CHECK(gra::mg5::ResolveEffectiveAlphaQCD(lts, 0.0, effective_alpha_s));
  CHECK(gra::math::IsZero(effective_alpha_s));
  CHECK(gra::mg5::ResolveEffectiveAlphaQCD(lts, 0.118, effective_alpha_s));
  CHECK(gra::math::IsExactEqual(effective_alpha_s, 0.118));
  lts.alphaQCD = 0.130;
  CHECK(gra::mg5::ResolveEffectiveAlphaQCD(lts, 0.118, effective_alpha_s));
  CHECK(gra::math::IsExactEqual(effective_alpha_s, 0.130));

  lts.alphaQCD = -0.1;
  CHECK_FALSE(gra::mg5::ResolveEffectiveAlphaQCD(lts, 0.118, effective_alpha_s));
  lts.alphaQCD = std::numeric_limits<double>::quiet_NaN();
  CHECK_FALSE(gra::mg5::ResolveEffectiveAlphaQCD(lts, 0.118, effective_alpha_s));
  lts.alphaQCD = 0.0;
  CHECK_FALSE(gra::mg5::ResolveEffectiveAlphaQCD(lts, -0.1, effective_alpha_s));
  CHECK_FALSE(gra::mg5::ResolveEffectiveAlphaQCD(lts, std::numeric_limits<double>::infinity(), effective_alpha_s));
}

// Check a valid topology change invalidates prepared event-local state
TEST_CASE("MG5 prepared event state includes the complete final topology", "[process][MG5][topology]") {
  const auto z_tree    = std::vector<gra::MDecayBranch>{CascadeBranch(23, {StableBranch(-13), StableBranch(13)})};
  const auto h_tree    = std::vector<gra::MDecayBranch>{CascadeBranch(25, {StableBranch(-13), StableBranch(13)})};
  const auto z_process = GeneratedProcess("z_mode", "z(mu+,mu-)", {{{23}, {{{-13}, {}}, {{13}, {}}}}}, {-13, 13},
                                          gra::amplitude::TopologyMode::Exact);
  const auto h_process = GeneratedProcess("h_mode", "h(mu+,mu-)", {{{25}, {{{-13}, {}}, {{13}, {}}}}}, {-13, 13},
                                          gra::amplitude::TopologyMode::Exact);
  gra::amplitude::ProcessRegistry definition({z_process, h_process});

  std::vector<gra::mg5::Subprocess<SubprocessTestMatrixElement>> channels;
  channels.push_back({std::make_unique<SubprocessTestMatrixElement>(),
                      {{{21, 21}, {-13, 13}, gra::amplitude::AmplitudeTopologyFromDecayTree(z_tree)}}});
  channels.push_back({std::make_unique<SubprocessTestMatrixElement>(),
                      {{{21, 21}, {-13, 13}, gra::amplitude::AmplitudeTopologyFromDecayTree(h_tree)}}});
  gra::mg5::SubprocessSum<SubprocessTestMatrixElement> subprocess_sum(std::move(channels));

  gra::LORENTZSCALAR z_event;
  z_event.q1        = gra::M4Vec(0.0, 0.0, 100.0, 100.0);
  z_event.q2        = gra::M4Vec(0.0, 0.0, -100.0, 100.0);
  z_event.decaytree = {ZMomentumBranch(gra::M4Vec(80.0, 0.0, 60.0, 100.0), gra::M4Vec(-80.0, 0.0, -60.0, 100.0))};
  BindProcessModelCache(z_event);
  gra::LORENTZSCALAR h_event = z_event;
  h_event.decaytree[0].p.pdg = 25;

  REQUIRE(subprocess_sum.Prepare(z_event, 0.118) == gra::mg5helas::EvaluationStatus::Success);
  REQUIRE(subprocess_sum.IsPrepared());
  REQUIRE(subprocess_sum.PreparedFinalStateMatches(z_event));
  REQUIRE(subprocess_sum.PreparedStateMatches(z_event, 0.118));
  CHECK_FALSE(subprocess_sum.PreparedStateMatches(z_event, 0.130));
  auto changed_initial = z_event;
  changed_initial.q1.SetPx(1.0e-12);
  CHECK_FALSE(subprocess_sum.PreparedStateMatches(changed_initial, 0.118));
  auto changed_final = z_event;
  changed_final.decaytree[0].legs[0].p4.SetPx(79.0);
  CHECK_FALSE(subprocess_sum.PreparedStateMatches(changed_final, 0.118));
  auto changed_decay_mode                    = z_event;
  changed_decay_mode.process.root_decay_mode = gra::RootDecayMode::Isolated;
  CHECK_FALSE(subprocess_sum.PreparedStateMatches(changed_decay_mode, 0.118));
  z_event.id1 = 1;
  z_event.id2 = -1;
  CHECK(subprocess_sum.PreparedStateMatches(z_event, 0.118));
  REQUIRE(definition.MatchProcess(h_event.decaytree).has_value());
  REQUIRE(subprocess_sum.RequireGeneratedTopology(definition, h_event));
  CHECK_FALSE(subprocess_sum.IsPrepared());
  CHECK_FALSE(subprocess_sum.PreparedFinalStateMatches(h_event));
}

// Check prepared wrappers refresh generated couplings when alpha_s changes
TEST_CASE("MG5 prepared event state includes the effective alpha_s", "[process][MG5][coupling]") {
  gra::LORENTZSCALAR first_event = ZJetKinematics();
  BindProcessModelCache(first_event);
  first_event.alphaQCD = 0.100;
  gra::AMP_MG5_pp_zj amplitude;
  REQUIRE(amplitude.Prepare(first_event, 0.100) == gra::mg5helas::EvaluationStatus::Success);
  const double first = amplitude.EvaluatePrepared(first_event, 0.100).amp2;

  gra::LORENTZSCALAR changed_event = first_event;
  changed_event.alphaQCD           = 0.200;
  const double       refreshed     = amplitude.EvaluatePrepared(changed_event, 0.200).amp2;
  gra::LORENTZSCALAR fresh_event   = changed_event;
  const double       fresh         = amplitude.EvaluatePrepared(fresh_event, 0.200).amp2;

  REQUIRE(first > 0.0);
  REQUIRE(refreshed > first);
  CHECK(refreshed == Approx(fresh).epsilon(1.0e-12));
}

// Check prepared hard amplitudes refresh the QED scheme at identical momenta
TEST_CASE("MG5 prepared event state includes the QED scheme", "[process][MG5][coupling]") {
  const auto directory = std::filesystem::path("tmp/test_mg5_prepared_qed");
  std::filesystem::create_directories(directory);
  const auto general = std::filesystem::path(gra::ResolveModelDataFile("TUNE0", "GENERAL.json"));
  for (const auto &filename : {"NUMERICS.json", "CON_MP.json", "CON_XP.json", "CON_GP.json", "CON_TP.json"}) {
    std::filesystem::copy_file(general.parent_path() / filename, directory / filename,
                              std::filesystem::copy_options::overwrite_existing);
  }
  gra::AMP_MG5_pp_z reused;
  auto event = DirectZKinematics();
  double previous = 0.0;
  for (const auto &scheme : {"ZERO", "MG", "ZERO"}) {
    auto card = nlohmann::json::parse(gra::aux::GetInputData(general.string()));
    card["PARAM_STRUCTURE"]["QED_alpha"] = scheme;
    const auto path = directory / "GENERAL.json";
    { std::ofstream output(path); output << card; }
    event.model_cache = std::make_shared<gra::MModelCache>(gra::MModelTune::Load(path.string()));
    gra::AMP_MG5_pp_z fresh;
    auto reference = event;
    const auto expected = fresh.EvaluatePrepared(reference, reference.alphaQCD);
    const auto result = reused.EvaluatePrepared(event, event.alphaQCD);
    REQUIRE(result.Valid());
    REQUIRE(expected.Valid());
    CHECK(result.amp2 == Approx(expected.amp2).epsilon(2e-12));
    REQUIRE(event.hamp.size() == reference.hamp.size());
    for (std::size_t i = 0; i < event.hamp.size(); ++i) {
      CHECK(std::abs(event.hamp[i] - reference.hamp[i]) <= 2e-12 * std::max(1.0, std::abs(reference.hamp[i])));
    }
    if (previous > 0.0) { CHECK(std::abs(result.amp2 - previous) > 1e-3 * previous); }
    previous = result.amp2;
  }
}

// Check physical Majorana poles preserve the signed SLHA convention
TEST_CASE("MG5 accepts signed Majorana poles", "[MG2GRA][SLHA]") {
  SLHAReader card;
  card.set_block_entry("mass", 1000022, -200.0);
  card.set_block_entry("decay", 1000022, -0.1);
  gra::mg5::ParticleMap poles{{1000022, {200.0, 0.1, true}}};
  REQUIRE_NOTHROW(gra::mg5::ValidateModel(card, poles));
  card.set_block_entry("mass", 1000022, 200.0);
  REQUIRE_THROWS_AS(gra::mg5::ValidateModel(card, poles), std::invalid_argument);
  card.set_block_entry("mass", 1000022, -200.0);
  poles.at(1000022).signed_mass = false;
  REQUIRE_THROWS_AS(gra::mg5::ValidateModel(card, poles), std::invalid_argument);
}

// Check sextet color lines close after crossing the incoming partons
TEST_CASE("MG5 validates sextet color flows", "[MG2GRA][color]") {
  using gra::mg5::ValidateExternalColorFlow;
  const std::vector<int> representations{3, 3, 6, 1};
  gra::mg5helas::ExternalColorFlow flow{{1, 0}, {2, 0}, {1, -2}, {0, 0}};
  REQUIRE_NOTHROW(ValidateExternalColorFlow(representations, flow));
  for (auto &leg : flow) { std::swap(leg.color, leg.anticolor); }
  REQUIRE_NOTHROW(ValidateExternalColorFlow({-3, -3, -6, 1}, flow));
  flow[2].color = -3;
  REQUIRE_THROWS_AS(ValidateExternalColorFlow({-3, -3, -6, 1}, flow), std::invalid_argument);
}

// Check model updates against a separate real gg to gamma gamma MadLoop runtime
TEST_CASE("MadLoop model parameters initialize physical poles", "[MG2GRA][MadLoop]") {
  const char *directory = std::getenv("GRANIITTI_TEST_MADLOOP");
  if (directory == nullptr) {
    WARN("Set GRANIITTI_TEST_MADLOOP to a generated gg to gamma gamma runtime");
    return;
  }
  gra::mg5::MadLoop loop(directory, "gg_aa", 4, 16, 1);
  INFO(loop.Error());
  REQUIRE(loop.Ready());
  const auto poles = loop.Particles();
  REQUIRE(poles.at(6).mass > 0.0);
  SLHAReader card(std::string(directory) + "/param_card.dat");
  const std::vector<double> momentum{250, 0, 0, 250, 250, 0, 0, -250,
                                      250, 200, 0, 150, 250, -200, 0, -150};
  gra::mg5::MadLoopResult before, changed, expected;
  REQUIRE(loop.Evaluate(momentum, before));
  card.set_block_entry("mass", 6, 1.2 * poles.at(6).mass);
  loop.InitParameters(card);
  CHECK(loop.Particles().at(6).mass == Approx(1.2 * poles.at(6).mass));
  REQUIRE(loop.Evaluate(momentum, changed));
  gra::mg5::MadLoop fresh(directory, "gg_aa", 4, 16, 1);
  fresh.InitParameters(card);
  REQUIRE(fresh.Evaluate(momentum, expected));
  double difference = 0.0;
  for (std::size_t i = 0; i < changed.amplitude.size(); ++i) {
    CHECK(std::abs(changed.amplitude[i] - expected.amplitude[i]) < 1e-10);
    difference += std::norm(changed.amplitude[i] - before.amplitude[i]);
  }
  CHECK(difference > 1e-15);
}

// Check direct family wrappers reject leaf-equivalent unregistered topologies
TEST_CASE("MG5 family wrappers reject direct topology bypasses", "[process][MG5][topology]") {
  SECTION("pp Z rejects an unregistered dimuon mother") {
    gra::LORENTZSCALAR valid = DirectZKinematics();
    gra::AMP_MG5_pp_z  amplitude;
    REQUIRE(amplitude.MatchProcess(valid.decaytree).has_value());

    gra::LORENTZSCALAR wrong     = valid;
    auto               daughters = std::move(wrong.decaytree);
    wrong.decaytree              = {CascadeBranch(25, std::move(daughters))};
    CHECK_FALSE(amplitude.MatchProcess(wrong.decaytree).has_value());
    CheckPreparedTopologyRejection(amplitude, valid, wrong);
  }

  SECTION("pp Z plus jet requires the registered Z cascade") {
    gra::LORENTZSCALAR valid = ZJetKinematics();
    gra::AMP_MG5_pp_zj amplitude;
    REQUIRE(amplitude.MatchProcess(valid.decaytree).has_value());

    gra::LORENTZSCALAR wrong = FlattenFirstCascade(valid);
    CHECK_FALSE(amplitude.MatchProcess(wrong.decaytree).has_value());
    CheckPreparedTopologyRejection(amplitude, valid, wrong);
  }

  SECTION("pp dijets require the registered two-parton topology") {
    gra::LORENTZSCALAR valid = DijetKinematics();
    gra::AMP_MG5_pp_jj amplitude;
    REQUIRE(amplitude.MatchProcess(valid.decaytree).has_value());

    gra::LORENTZSCALAR wrong = valid;
    wrong.decaytree.push_back(StableBranch(21));
    CHECK_FALSE(amplitude.MatchProcess(wrong.decaytree).has_value());
    CheckPreparedTopologyRejection(amplitude, valid, wrong);
  }

  SECTION("gamma-gamma Zjj rejects topology before HELAS and color flow") {
    gra::LORENTZSCALAR valid = ZDijetKinematics();
    valid.id1                = 22;
    valid.id2                = 22;
    gra::AMP_MG5_yy_zjj amplitude;
    REQUIRE(amplitude.MatchProcess(valid.decaytree).has_value());

    gra::LORENTZSCALAR wrong = FlattenFirstCascade(valid);
    REQUIRE(wrong.decaytree.size() == 4);
    CHECK_FALSE(amplitude.MatchProcess(wrong.decaytree).has_value());
    wrong.hamp                            = {1.0};
    wrong.decaytree[2].p.color_flow.flow1 = 701;
    wrong.decaytree[3].p.color_flow.flow2 = 702;
    const auto result                     = amplitude.Evaluate(wrong, 0.0, false);
    CHECK(result.status == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
    CHECK(gra::math::IsZero(result.amp2));
    CHECK(wrong.hamp.empty());

    gra::MRandom random;
    CHECK_FALSE(amplitude.SampleColorFlow(wrong, random));
    CHECK(wrong.decaytree[2].p.color_flow.flow1 == 701);
    CHECK(wrong.decaytree[3].p.color_flow.flow2 == 702);

    gra::LORENTZSCALAR reordered = valid;
    std::swap(reordered.decaytree[1], reordered.decaytree[2]);
    CHECK_FALSE(amplitude.MatchProcess(reordered.decaytree).has_value());
    CHECK(amplitude.Evaluate(reordered, 0.0, false).status == gra::mg5helas::EvaluationStatus::AmplitudeFailure);
  }
}
