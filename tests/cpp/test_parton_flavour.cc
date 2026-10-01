// Unit tests for generic final-state parton selection
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <iterator>
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

// Own
#include "Graniitti/Amplitude/MG5/Runtime/read_slha.h"
#include "Graniitti/MGlobals.h"
#include "Graniitti/Kinematics/MFactorized.h"
#include "Graniitti/QCD/MPartonFlavour.h"
#include "Graniitti/QCD/MPartonProposal.h"
#include "Graniitti/Tech/MAux.h"

// Libraries
#include <catch.hpp>

namespace {

// Restore the process tune after a focused selector test
struct ModelParamRestoreGuard {
  std::string value = gra::MODELPARAM;

  // Restore the saved model tune when the guard leaves scope
  ~ModelParamRestoreGuard() { gra::MODELPARAM = value; }
};

// Expose process setup and final-state proposal state for focused tests
class PartonProposalProbe : public gra::MFactorized {
public:
  // Construct a deterministic proposal probe
  PartonProposalProbe() { state.random.SetSeed(24680); }

  // Configure one real subprocess and its central decay tree
  void ConfigureProcess(const std::string &process, const std::string &decay,
                        const std::string &istate,
                        const std::string &phase_space,
                        const std::string &selector = "") {
    ProcPtr = gra::MSubProc({istate}, phase_space);
    std::vector<gra::aux::OneCMD> syntax;
    if (!selector.empty()) {
      syntax = gra::aux::SplitCommands(selector);
    }
    std::string mutable_process = process;
    SetProcess(mutable_process, syntax);
    SetDecayMode(decay);
  }

  // Draw one physical final-state parton proposal
  void Propose() {
    if (gra::parton::PrepareFinalStateProposal(state.lts, state.random)) {
      CalculateSymmetryFactor();
    }
  }

  // Compute stable central leaves in decay-tree order
  std::vector<const gra::MDecayBranch *> Leaves() const {
    std::vector<const gra::MDecayBranch *> leaves;
    for (const auto &branch : state.lts.decaytree) {
      CollectLeaves(branch, leaves);
    }
    return leaves;
  }

  // Compute the current inverse proposal probability
  double ProposalWeight() const {
    return gra::parton::ProposalWeight(state.lts);
  }

  // Compute the current physical identical-particle symmetry factor
  double SymmetryFactor() const { return state.symmetry_factor; }

  // Compute whether the canonical decay tree selects a generated matrix element
  bool MatchesMG5() {
    const auto process = ProcPtr.MatchProcess(state.lts.decaytree);
    return process.has_value() &&
           gra::amplitude::IsGenerated(process->matrix_element_form);
  }

private:
  // Collect stable leaves recursively
  static void CollectLeaves(const gra::MDecayBranch &branch,
                            std::vector<const gra::MDecayBranch *> &leaves) {
    if (branch.legs.empty()) {
      leaves.push_back(&branch);
      return;
    }
    for (const auto &leg : branch.legs) {
      CollectLeaves(leg, leaves);
    }
  }
};

// Compute stable-leaf PDG ids from one proposal probe
std::vector<int> LeafPDGs(const PartonProposalProbe &probe) {
  std::vector<int> pdgs;
  for (const auto *leaf : probe.Leaves()) {
    pdgs.push_back(leaf->p.pdg);
  }
  return pdgs;
}

// Append stable PDGs from one exact generated process branch
void AppendExactProcessPDGs(const gra::AmplitudeTopologyNode &node,
                            std::vector<int> &pdgs) {
  if (node.daughters.empty()) {
    REQUIRE(node.allowed_pdgs.size() == 1);
    pdgs.push_back(node.allowed_pdgs.front());
    return;
  }
  for (const auto &daughter : node.daughters) {
    AppendExactProcessPDGs(daughter, pdgs);
  }
}

// Compute deduplicated stable modes represented by exact generated topologies
std::set<std::vector<int>>
GeneratedStableModes(const gra::amplitude::Process &process) {
  std::set<std::vector<int>> modes;
  for (const auto &topology : process.channel_topologies) {
    std::vector<int> pdgs;
    for (const auto &node : topology) {
      AppendExactProcessPDGs(node, pdgs);
    }
    modes.insert(std::move(pdgs));
  }
  return modes;
}

// Build one exact MG5 decay branch from generated topology data
gra::MDecayBranch ExactMG5Branch(const gra::AmplitudeTopologyNode &node,
                                 const gra::MPDG &pdg) {
  if (node.allowed_pdgs.size() != 1) {
    throw std::logic_error("ExactMG5Branch: topology is not exact");
  }
  gra::MDecayBranch branch;
  branch.p = pdg.FindByPDG(node.allowed_pdgs.front());
  branch.name = std::to_string(branch.p.pdg);
  for (const auto &daughter : node.daughters) {
    branch.legs.push_back(ExactMG5Branch(daughter, pdg));
  }
  return branch;
}

// Build one exact MG5 decay tree from generated topology data
std::vector<gra::MDecayBranch>
ExactMG5Tree(const gra::AmplitudeTopology &topology, const gra::MPDG &pdg) {
  std::vector<gra::MDecayBranch> tree;
  tree.reserve(topology.size());
  for (const auto &node : topology) {
    tree.push_back(ExactMG5Branch(node, pdg));
  }
  return tree;
}

// Compute every unique recursive sibling permutation of one MG5 decay tree
std::vector<std::vector<gra::MDecayBranch>>
MG5TreePermutations(const std::vector<gra::MDecayBranch> &tree) {
  std::vector<std::vector<gra::MDecayBranch>> combinations(1);
  for (const auto &branch : tree) {
    const auto daughter_permutations = MG5TreePermutations(branch.legs);
    std::vector<std::vector<gra::MDecayBranch>> next;
    for (const auto &combination : combinations) {
      for (const auto &daughters : daughter_permutations) {
        auto expanded = combination;
        auto node = branch;
        node.legs = daughters;
        expanded.push_back(std::move(node));
        next.push_back(std::move(expanded));
      }
    }
    combinations = std::move(next);
  }

  std::vector<std::vector<gra::MDecayBranch>> permutations;
  std::vector<gra::AmplitudeTopology> seen;
  std::vector<std::size_t> order;
  for (const auto &i : gra::aux::indices(tree)) {
    order.push_back(i);
  }
  for (const auto &combination : combinations) {
    do {
      std::vector<gra::MDecayBranch> permuted;
      permuted.reserve(order.size());
      for (const auto &i : order) {
        permuted.push_back(combination[i]);
      }
      const auto topology =
          gra::amplitude::AmplitudeTopologyFromDecayTree(permuted);
      if (std::find(seen.begin(), seen.end(), topology) == seen.end()) {
        seen.push_back(topology);
        permutations.push_back(std::move(permuted));
      }
    } while (std::next_permutation(order.begin(), order.end()));
  }
  return permutations;
}

// Check reordered MG5 branches retain data belonging to their own PDG id
bool MG5BranchDataMatches(const gra::MDecayBranch &branch) {
  return branch.name == std::to_string(branch.p.pdg) &&
         std::all_of(branch.legs.begin(), branch.legs.end(),
                     [](const auto &daughter) {
                       return MG5BranchDataMatches(daughter);
                     });
}

// Check one stable leaf against its phase-space mass shell
void RequireConfiguredMassShell(const gra::MDecayBranch &leaf) {
  gra::M4Vec p4;
  p4.SetPxPyPzM(3.0, 4.0, 5.0, leaf.p.mass);
  REQUIRE(p4.M2() == Approx(leaf.p.mass * leaf.p.mass).margin(1e-12));
}

} // namespace

// Check symbolic and numeric final-state selector values
TEST_CASE("Final-state parton selector accepts physical species",
          "[parton-selector]") {
  const auto commands =
      gra::aux::SplitCommands("gg[QCD]<F> -> j ~j @j={u,1,s,4,b,21}");
  REQUIRE(commands.size() == 1);
  REQUIRE(commands.front().id == "j");
  REQUIRE((commands.front().values ==
           std::vector<std::string>{"u", "1", "s", "4", "b", "21"}));

  gra::MPartonFlavour selector;
  selector.Configure(commands.front().values);
  REQUIRE((selector.Flavours() == std::vector<int>{2, 1, 3, 4, 5, 21}));
  REQUIRE((selector.QuarkFlavours() == std::vector<int>{2, 1, 3, 4, 5}));
  REQUIRE(selector.Accepts(-2));
  REQUIRE(selector.Accepts(21));
  REQUIRE_FALSE(selector.Accepts(-21));
  REQUIRE(selector.IsConfigured());
}

// Check one allowed trailing comma and the ordinary no-command path
TEST_CASE(
    "Final-state parton selector parser has an unambiguous optional comma",
    "[parton-selector]") {
  const auto commands = gra::aux::SplitCommands("@j={u,d,s,}");
  REQUIRE(commands.size() == 1);
  REQUIRE((commands.front().values == std::vector<std::string>{"u", "d", "s"}));
  REQUIRE(gra::aux::SplitCommands("gg[QCD]<F>").empty());
}

// Check malformed selector syntax and unsupported species
TEST_CASE("Final-state parton selector rejects ambiguous input",
          "[parton-selector]") {
  REQUIRE_THROWS(gra::aux::SplitCommands("@j={u,d} @j={s}"));
  REQUIRE_THROWS(gra::aux::SplitCommands("@j={u,,d}"));
  REQUIRE_THROWS(gra::aux::SplitCommands("@j={u,d"));
  REQUIRE_THROWS(gra::aux::SplitCommands("@j={u,d} trailing"));
  REQUIRE_THROWS(gra::aux::SplitCommands("@q={u,d}"));

  gra::MPartonFlavour selector;
  REQUIRE_THROWS((selector.Configure({"u", "2"})));
  REQUIRE_THROWS(selector.Configure({"-1"}));
  REQUIRE_THROWS(selector.Configure({"6"}));
  REQUIRE_THROWS(selector.Configure({"22"}));
  REQUIRE_THROWS(selector.Configure({""}));
}

// Check lone-jet charge modes and correlated pair symmetry compensation
TEST_CASE(
    "Final-state parton proposal closes physical species and symmetry sums",
    "[parton-selector][proposal]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  PartonProposalProbe explicit_state;
  explicit_state.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} u", "IPp", "F");
  explicit_state.Propose();
  REQUIRE(explicit_state.ProposalWeight() == Approx(1.0));
  REQUIRE(explicit_state.Leaves().back()->p.pdg == 2);

  PartonProposalProbe default_lone;
  default_lone.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} j", "IPp", "F");
  default_lone.Propose();
  REQUIRE(default_lone.ProposalWeight() == Approx(11.0));

  PartonProposalProbe lone;
  lone.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} j", "IPp", "F",
                        "@j={u,d,g}");
  std::set<int> lone_modes;
  for (std::size_t draw = 0; draw < 400; ++draw) {
    lone.Propose();
    const auto leaves = lone.Leaves();
    REQUIRE(leaves.size() == 3);
    lone_modes.insert(leaves.back()->p.pdg);
    REQUIRE(lone.ProposalWeight() == Approx(5.0));
    REQUIRE(lone.SymmetryFactor() == Approx(1.0));
  }
  REQUIRE((lone_modes == std::set<int>{-2, -1, 1, 2, 21}));

  PartonProposalProbe pair;
  pair.ConfigureProcess("gg[QCD]<F>", "j ~j", "gg", "F", "@j={c,b,g}");
  std::set<std::vector<int>> pair_modes;
  for (std::size_t draw = 0; draw < 300; ++draw) {
    pair.Propose();
    const auto pdgs = LeafPDGs(pair);
    pair_modes.insert(pdgs);
    REQUIRE(pair.ProposalWeight() == Approx(3.0));
    if (pdgs == std::vector<int>{21, 21}) {
      REQUIRE(pair.SymmetryFactor() == Approx(2.0));
    } else {
      REQUIRE(pair.SymmetryFactor() == Approx(1.0));
    }
  }
  REQUIRE(
      (pair_modes == std::set<std::vector<int>>{{4, -4}, {5, -5}, {21, 21}}));

  PartonProposalProbe default_pair;
  default_pair.ConfigureProcess("gg[QCD]<F>", "j ~j", "gg", "F");
  default_pair.Propose();
  REQUIRE(default_pair.ProposalWeight() == Approx(6.0));
}

// Check generic hard dijets sample only exact correlated MG5 modes
TEST_CASE("Generated hard-dijet proposals follow exact family topologies",
          "[parton-selector][proposal][MG5]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const auto processes = gra::amplitude::Processes("MG5_PP_JJ");
  REQUIRE(processes.size() == 1);
  const auto exact_modes = GeneratedStableModes(processes.front());
  REQUIRE(exact_modes.size() == 66);

  PartonProposalProbe generic;
  generic.ConfigureProcess("IPIP[jj]<F>", "j j", "IPIP", "F");
  std::set<std::vector<int>> sampled;
  for (std::size_t draw = 0; draw < 1000; ++draw) {
    generic.Propose();
    const auto mode = LeafPDGs(generic);
    sampled.insert(mode);
    REQUIRE(exact_modes.count(mode) == 1);
    REQUIRE(generic.ProposalWeight() == Approx(66.0));
  }
  CHECK(sampled.size() > 50);

  std::set<std::vector<int>> selected_modes;
  std::copy_if(exact_modes.begin(), exact_modes.end(),
               std::inserter(selected_modes, selected_modes.end()),
               [](const auto &mode) {
                 const std::set<int> selected = {-2, 2, 21};
                 return mode.size() == 2 && selected.count(mode[0]) == 1 &&
                        selected.count(mode[1]) == 1;
               });
  REQUIRE(selected_modes.size() == 6);

  PartonProposalProbe selected;
  selected.ConfigureProcess("IPIP[jj]<F>", "j j", "IPIP", "F",
                            "@j={u,g}");
  std::set<std::vector<int>> selected_sample;
  for (std::size_t draw = 0; draw < 300; ++draw) {
    selected.Propose();
    const auto mode = LeafPDGs(selected);
    selected_sample.insert(mode);
    REQUIRE(selected_modes.count(mode) == 1);
    REQUIRE(selected.ProposalWeight() == Approx(6.0));
  }
  CHECK(selected_sample == selected_modes);

  PartonProposalProbe correlated;
  correlated.ConfigureProcess("IPIP[jj]<F>", "j ~j", "IPIP", "F");
  correlated.Propose();
  REQUIRE(exact_modes.count(LeafPDGs(correlated)) == 1);
  REQUIRE(correlated.ProposalWeight() == Approx(6.0));
}

// Check that flattened aliases under different parents remain independent
TEST_CASE("Final-state parton pairs correlate only between immediate siblings",
          "[parton-selector][proposal]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  PartonProposalProbe disjoint_parents;
  disjoint_parents.ConfigureProcess("MP[RES]<F>", "Z > {j} Z > {~j}", "MP", "F",
                                    "@j={u}");
  std::set<std::vector<int>> modes;
  for (std::size_t draw = 0; draw < 200; ++draw) {
    disjoint_parents.Propose();
    modes.insert(LeafPDGs(disjoint_parents));
    REQUIRE(disjoint_parents.ProposalWeight() == Approx(4.0));
  }
  REQUIRE((modes ==
           std::set<std::vector<int>>{{-2, -2}, {-2, 2}, {2, -2}, {2, 2}}));
}

// Check generated matrix elements follow each registered ordering
TEST_CASE(
    "Generated parton-pair matrix elements follow registered ordering",
    "[parton-selector][proposal]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  PartonProposalProbe photon;
  REQUIRE_THROWS(
      photon.ConfigureProcess("yy[EPA]<F>", "~j j", "yy", "F", "@j={u}"));

  PartonProposalProbe reversed_hard_dijet;
  REQUIRE_THROWS(reversed_hard_dijet.ConfigureProcess(
      "IPIP[jj]<F>", "~j j", "IPIP", "F", "@j={u}"));

  PartonProposalProbe hard_dijet;
  REQUIRE_NOTHROW(hard_dijet.ConfigureProcess("IPIP[jj]<F>", "j ~j", "IPIP",
                                              "F", "@j={u}"));
  hard_dijet.Propose();
  REQUIRE((LeafPDGs(hard_dijet) == std::vector<int>{2, -2}));
}

// Canonicalize every unique top-level and nested MG5 process permutation
TEST_CASE("MG5 process syntax follows the unique generated particle ordering",
          "[parton-selector][MG5][ordering]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  PartonProposalProbe direct_zj;
  REQUIRE_NOTHROW(
      direct_zj.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} u", "IPp", "F"));

  PartonProposalProbe reversed_top_level;
  REQUIRE_NOTHROW(reversed_top_level.ConfigureProcess(
      "IPp[Zj]<F>", "u Z > {mu+ mu-}", "IPp", "F"));
  REQUIRE((LeafPDGs(reversed_top_level) == std::vector<int>{-13, 13, 2}));
  REQUIRE(reversed_top_level.MatchesMG5());

  PartonProposalProbe reversed_daughters;
  REQUIRE_NOTHROW(reversed_daughters.ConfigureProcess(
      "IPp[Zj]<F>", "Z > {mu- mu+} u", "IPp", "F"));
  REQUIRE((LeafPDGs(reversed_daughters) == std::vector<int>{-13, 13, 2}));
  REQUIRE(reversed_daughters.MatchesMG5());

  PartonProposalProbe reversed_pair;
  REQUIRE_NOTHROW(
      reversed_pair.ConfigureProcess("IPIP[Z]<F>", "mu- mu+", "IPIP", "F"));
  REQUIRE((LeafPDGs(reversed_pair) == std::vector<int>{-13, 13}));
  REQUIRE(reversed_pair.MatchesMG5());

  PartonProposalProbe reversed_photon_pair;
  REQUIRE_NOTHROW(
      reversed_photon_pair.ConfigureProcess("yy[EPA]<F>", "e- e+", "yy", "F"));
  REQUIRE((LeafPDGs(reversed_photon_pair) == std::vector<int>{-11, 11}));
  REQUIRE(reversed_photon_pair.MatchesMG5());

  PartonProposalProbe canonical_three_body;
  REQUIRE_NOTHROW(
      canonical_three_body.ConfigureProcess("yy[EPA]<F>", "u u~ g", "yy", "F"));
  REQUIRE((LeafPDGs(canonical_three_body) == std::vector<int>{2, -2, 21}));
  REQUIRE(canonical_three_body.MatchesMG5());

  PartonProposalProbe permuted_three_body;
  REQUIRE_NOTHROW(
      permuted_three_body.ConfigureProcess("yy[EPA]<F>", "g u u~", "yy", "F"));
  REQUIRE((LeafPDGs(permuted_three_body) == std::vector<int>{2, -2, 21}));
  REQUIRE(permuted_three_body.MatchesMG5());

  PartonProposalProbe unsupported;
  REQUIRE_THROWS_AS(
      unsupported.ConfigureProcess("IPIP[Z]<F>", "pi+ pi-", "IPIP", "F"),
      std::invalid_argument);

  PartonProposalProbe missing_particle;
  REQUIRE_THROWS_AS(missing_particle.ConfigureProcess(
                        "IPp[Zj]<F>", "Z > {mu+ mu-}", "IPp", "F"),
                    std::invalid_argument);

  PartonProposalProbe duplicate_particle;
  REQUIRE_THROWS_AS(
      duplicate_particle.ConfigureProcess("yy[EPA]<F>", "u u g", "yy", "F"),
      std::invalid_argument);

  PartonProposalProbe wrong_decay_branch;
  REQUIRE_THROWS_AS(wrong_decay_branch.ConfigureProcess(
                        "IPp[Zj]<F>", "mu+ Z > {mu- u}", "IPp", "F"),
                    std::invalid_argument);
}

// Check a recursive permutation of every exact generated topology remains safe
TEST_CASE("Every MG5 final state canonicalizes to an exact generated channel",
          "[parton-selector][MG5][ordering][registry]") {
  gra::MPDG pdg;
  pdg.ReadParticleData();

  std::set<std::string> families;
  for (const auto &process : gra::amplitude::AllProcesses()) {
    families.insert(process.process_family);
  }

  for (const auto &family : families) {
    const auto processes = gra::amplitude::Processes(family);
    const gra::amplitude::ProcessRegistry registry(processes);
    for (const auto &process : processes) {
      const std::vector<gra::AmplitudeTopology> topologies =
          process.channel_topologies.empty()
              ? std::vector<gra::AmplitudeTopology>{process.topology}
              : process.channel_topologies;
      for (const auto &topology : topologies) {
        CAPTURE(family, process.process_name, process.final_state_syntax);
        const auto canonical = ExactMG5Tree(topology, pdg);
        for (auto tree : MG5TreePermutations(canonical)) {
          gra::amplitude::OrderMG5ProcessSyntax(processes, tree,
                                                process.final_state_syntax);
          REQUIRE(registry.MatchProcess(tree).has_value());
          REQUIRE(std::all_of(tree.begin(), tree.end(), [](const auto &branch) {
            return MG5BranchDataMatches(branch);
          }));
        }
      }
    }
  }
}

// Check generated external masses for explicit and generic final states
TEST_CASE("Final-state partons use process-local generated mass schemes",
          "[parton-selector][mass-shell]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  PartonProposalProbe analytic_epa;
  analytic_epa.ConfigureProcess("yy[EPA]<F>", "c c~", "yy", "F");
  REQUIRE(analytic_epa.Leaves().front()->p.mass > 0.0);

  PartonProposalProbe durham;
  durham.ConfigureProcess("gg[QCD]<F>", "c c~", "gg", "F");
  REQUIRE(durham.Leaves().front()->p.mass == Approx(0.0));
  RequireConfiguredMassShell(*durham.Leaves().front());
  durham.SetDecayMode("b b~");
  const double configured_bottom_mass = durham.Leaves().front()->p.mass;
  REQUIRE(configured_bottom_mass > 0.0);
  RequireConfiguredMassShell(*durham.Leaves().front());

  PartonProposalProbe photon;
  photon.ConfigureProcess("yy[Zjj]<F>", "Z > {mu+ mu-} c c~", "yy", "F");
  REQUIRE(photon.Leaves().at(2)->p.mass == Approx(0.0));
  RequireConfiguredMassShell(*photon.Leaves().at(2));
  photon.SetDecayMode("Z > {mu+ mu-} b b~");
  REQUIRE(photon.Leaves().at(2)->p.mass == Approx(configured_bottom_mass));
  RequireConfiguredMassShell(*photon.Leaves().at(2));

  PartonProposalProbe single_diffraction;
  single_diffraction.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} b", "IPp",
                                      "F");
  REQUIRE(single_diffraction.Leaves().back()->p.mass == Approx(0.0));
  RequireConfiguredMassShell(*single_diffraction.Leaves().back());

  PartonProposalProbe double_diffraction;
  double_diffraction.ConfigureProcess("IPIP[jj]<F>", "c c~", "IPIP", "F");
  REQUIRE(double_diffraction.Leaves().front()->p.mass == Approx(0.0));
  RequireConfiguredMassShell(*double_diffraction.Leaves().front());

  PartonProposalProbe selected;
  selected.ConfigureProcess("yy[Zjj]<F>", "Z > {mu+ mu-} j ~j", "yy", "F",
                            "@j={b}");
  selected.Propose();
  REQUIRE((LeafPDGs(selected) == std::vector<int>{-13, 13, 5, -5}));
  REQUIRE(selected.Leaves().at(2)->p.mass == Approx(configured_bottom_mass));
  RequireConfiguredMassShell(*selected.Leaves().at(2));
}

// Compare explicit and generic Durham masses before phase space construction
TEST_CASE("Generic Durham modes preserve the generated quark mass schemes", "[parton-selector][mass-shell][durham]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  for (const std::string suffix : {"", " g"}) {
    std::map<int, double> masses;
    for (const auto &[pdg, name] : std::map<int, std::string>{{1, "d"}, {2, "u"}, {3, "s"}, {4, "c"}, {5, "b"}}) {
      PartonProposalProbe explicit_mode;
      explicit_mode.ConfigureProcess("gg[QCD]<F>", name + " " + name + "~" + suffix, "gg", "F");
      masses[pdg] = explicit_mode.Leaves().front()->p.mass;

      PartonProposalProbe generic;
      generic.ConfigureProcess("gg[QCD]<F>", "j ~j" + suffix, "gg", "F", "@j={" + name + "}");
      generic.Propose();
      REQUIRE(generic.MatchesMG5());
      REQUIRE(generic.Leaves().front()->p.pdg == pdg);
      for (const auto *leaf : generic.Leaves()) {
        if (std::abs(leaf->p.pdg) != pdg) { continue; }
        REQUIRE(leaf->p.mass == Approx(masses.at(pdg)).margin(1e-12));
        RequireConfiguredMassShell(*leaf);
      }
    }

    PartonProposalProbe mixed;
    mixed.ConfigureProcess("gg[QCD]<F>", "j ~j" + suffix, "gg", "F", "@j={u,d,s,c,b}");
    std::set<int> seen;
    for (int draw = 0; draw < 100; ++draw) {
      mixed.Propose();
      const auto leaves = mixed.Leaves();
      const int  pdg    = leaves.front()->p.pdg;
      seen.insert(pdg);
      REQUIRE(mixed.MatchesMG5());
      REQUIRE(leaves[0]->p.mass == Approx(masses.at(pdg)).margin(1e-12));
      REQUIRE(leaves[1]->p.mass == Approx(masses.at(pdg)).margin(1e-12));
      REQUIRE(mixed.ProposalWeight() == Approx(5.0));
    }
    REQUIRE(seen.size() == masses.size());
  }
}

// Check process setup synchronizes internal resonances from generated cards
TEST_CASE("Generated process setup synchronizes resonance parameters",
          "[parton-selector][MG5][decay]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";

  const SLHAReader z_card(
      gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_ZJJ/param_card.dat"));
  const auto z_mass = z_card.find_block_entry("mass", 23);
  const auto z_width = z_card.find_block_entry("decay", 23);
  REQUIRE(z_mass.has_value());
  REQUIRE(z_width.has_value());

  PartonProposalProbe zj;
  zj.ConfigureProcess("yy[Zjj]<F>", "Z > {mu+ mu-} c c~", "yy", "F");
  REQUIRE(zj.state.lts.decaytree.size() == 3);
  CHECK(zj.state.lts.decaytree.front().p.mass == Approx(*z_mass));
  CHECK(zj.state.lts.decaytree.front().p.width == Approx(*z_width));
  REQUIRE(zj.state.lts.final_state_parton_masses.size() == 1);
  CHECK(zj.state.lts.final_state_parton_masses.contains(4));

  const SLHAReader w_card(
      gra::aux::ResolveProjectPath("MG5cards/Photon/MG5_YY_WW/param_card.dat"));
  const auto w_mass = w_card.find_block_entry("mass", 24);
  const auto w_width = w_card.find_block_entry("decay", 24);
  REQUIRE(w_mass.has_value());
  REQUIRE(w_width.has_value());

  PartonProposalProbe ww;
  ww.ConfigureProcess("yy[WW]<F>", "W+ > {e+ ve} W- > {mu- vm~}", "yy", "F");
  REQUIRE(ww.state.lts.decaytree.size() == 2);
  CHECK(ww.state.lts.final_state_parton_masses.empty());
  for (const auto &branch : ww.state.lts.decaytree) {
    CHECK(branch.p.mass == Approx(*w_mass));
    CHECK(branch.p.width == Approx(*w_width));
  }
}

// Reject external and nested pole requests that conflict with generated model inputs
TEST_CASE("Generated processes reject incompatible PDG parameters during configuration", "[parton-selector][MG5][model]") {
  ModelParamRestoreGuard restore;
  gra::MODELPARAM = "TUNE0";
  for (const std::string phase : {"F", "C"}) {
    for (const std::string request : {"@PDG[11]{M:0}", "@PDG[24]{M:1}", "@PDG[24]{W:0}", "@PDG[11]{mass:0}"}) {
      CAPTURE(phase, request);
      PartonProposalProbe process;
      REQUIRE_THROWS_AS(process.ConfigureProcess("yy[WW]<" + phase + ">",
          "W+ > {e+ ve} W- > {mu- vm~}", "yy", phase, request), std::invalid_argument);
    }
    PartonProposalProbe valid;
    REQUIRE_NOTHROW(valid.ConfigureProcess("yy[WW]<" + phase + ">",
        "W+ > {e+ ve} W- > {mu- vm~}", "yy", phase));
    const auto model = gra::CreatePhotonMG5Process("MG5_YY_WW")->Particles();
    for (const auto *leaf : valid.Leaves()) {
      REQUIRE(leaf->p.mass == Approx(model.at(std::abs(leaf->p.pdg)).mass));
    }
  }
  PartonProposalProbe parton;
  REQUIRE_THROWS_AS(parton.ConfigureProcess("IPp[Zj]<F>", "Z > {mu+ mu-} b", "IPp", "F",
                                           "@PDG[13]{M:0}"), std::invalid_argument);
}
