// Process help derived from amplitude definitions and model cards
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "Graniitti/Process/MProcessTable.h"

#include <algorithm>
#include <cmath>
#include <filesystem>
#include <set>
#include <tuple>

#include "Graniitti/Nuclear/MBeam.h"
#include "Graniitti/Regge/MReggeInit.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {
namespace {

// Identify nonzero production or decay couplings in the selected basis
bool ActiveCoupling(const nlohmann::json& value, double minimum) {
  if (value.is_array()) {
    return std::any_of(value.begin(), value.end(), [minimum](const auto& child) {
      return ActiveCoupling(child, minimum);
    });
  }
  if (!value.is_object()) { return false; }
  if (value.contains("BR")) { return value.at("BR").get<double>() > minimum; }
  if (value.contains("g")) {
    return value.at("g").at(0).is_null() || std::abs(value.at("g").at(0).get<double>()) > minimum;
  }
  const std::string basis = value.value("basis", std::string{});
  for (const std::string key : {"g_tensor", "g_ls", "helicity", "alpha_ls"}) {
    if (!value.contains(key)) { continue; }
    if ((key == "helicity" || key == "alpha_ls") && basis != key) { continue; }
    for (const auto& row : value.at(key)) {
      const auto& magnitude = row.is_array() ? row.at(2) : row;
      if (magnitude.is_null() || std::abs(magnitude.get<double>()) > minimum) { return true; }
    }
  }
  return std::any_of(value.begin(), value.end(),
                     [minimum](const auto& child) { return ActiveCoupling(child, minimum); });
}

// Compute a flavour summary while retaining the registered decay topology
std::string StatePattern(const AmplitudeTopology& topology, const MPDG& pdg) {
  std::string out;
  for (const auto& node : topology) {
    const int code = node.allowed_pdgs.front();
    const int abs_code = std::abs(code);
    std::string name;
    if (node.allowed_pdgs.size() > 1) {
      name = "j";
    } else if (abs_code >= 1 && abs_code <= 6) {
      name = code > 0 ? "q" : "qbar";
    } else if (abs_code == 11 || abs_code == 13 || abs_code == 15) {
      name = code > 0 ? "l-" : "l+";
    } else if (abs_code == 12 || abs_code == 14 || abs_code == 16) {
      name = code > 0 ? "nu" : "nubar";
    } else {
      name = pdg.FindByPDG(code).name;
    }
    if (!node.daughters.empty()) { name += " > {" + StatePattern(node.daughters, pdg) + "}"; }
    out += (out.empty() ? "" : " ") + name;
  }
  if (out == "l+ l-" || out == "l- l+") { return "l+l-"; }
  if (out == "q qbar" || out == "qbar q") { return "qqbar"; }
  return out;
}

// Join distinct amplitude patterns without duplicating generated flavour channels
std::string StateSummary(const MProc& process, const MPDG& pdg) {
  if (process.ISTATE == "X") { return process.DESCRIPTION.process; }
  std::set<std::string> patterns;
  std::map<std::string, std::set<std::string>> topologies;
  for (const auto& state : process.Processes()) {
    if (!state.topology.empty()) {
      topologies[StatePattern(state.topology, pdg)].insert(state.final_state_syntax);
    } else {
      for (const auto& pattern : aux::SplitStr2Str(state.final_state_syntax)) { patterns.insert(pattern); }
    }
  }
  for (const auto& [pattern, states] : topologies) {
    if (patterns.contains(pattern)) { continue; }
    patterns.insert(states.size() > 1 ? pattern : *states.begin());
  }
  std::string out;
  for (const auto& pattern : patterns) { out += (out.empty() ? "" : ", ") + pattern; }
  return out;
}

}  // namespace

// Read particle names, production channels and decays once for the help tables
MProcessTable::MProcessTable(const std::string& modelparam) {
  pdg.ReadParticleData(MPDG::DataFile(), modelparam);
  general              = nlohmann::json::parse(aux::GetInputData(ResolveModelDataFile(modelparam, "GENERAL.json")));
  decays               = nlohmann::json::parse(aux::GetInputData(ResolveModelDataFile(modelparam, "DECAYS.json")));
  const auto numerics  = nlohmann::json::parse(aux::GetInputData(ResolveModelDataFile(modelparam, "NUMERICS.json")));
  coupling_min         = numerics.at("NUMERICS_GLOBAL").at("coupling_min").get<double>();
  const auto directory = std::filesystem::path(ResolveModelDataFile(modelparam, "GENERAL.json")).parent_path() / "RES";
  for (const auto& entry : std::filesystem::directory_iterator(directory)) {
    if (entry.path().extension() != ".json") { continue; }
    const auto card = nlohmann::json::parse(aux::GetInputData(entry.path().string()));
    resonances.emplace(entry.path().stem().string(), card.at("PARAM_RES"));
  }
}

// Construct a physical decay tree from a signed PDG channel key
std::vector<MDecayBranch> MProcessTable::Tree(const std::string& key) const {
  std::vector<MDecayBranch> tree;
  for (const int code : nlohmann::json::parse(key).get<std::vector<int>>()) {
    MDecayBranch branch;
    branch.p = pdg.FindByPDG(code);
    tree.push_back(std::move(branch));
  }
  return tree;
}

// Expand unstable vector daughters using their largest configured branching ratios
bool MProcessTable::Expand(std::vector<MDecayBranch>& tree) const {
  bool expanded = false;
  for (auto& branch : tree) {
    const auto channels = decays.find(std::to_string(branch.p.pdg));
    if (branch.p.spinX2 != 2 || channels == decays.end()) { continue; }
    double br = 0.0;
    for (const auto& [key, decay] : channels->items()) {
      if (key.front() != '[' || !decay.contains("BR")) { continue; }
      const double next = decay.at("BR").get<double>();
      if (next > br && next > coupling_min) {
        branch.legs = Tree(key);
        br          = next;
      }
    }
    expanded = expanded || !branch.legs.empty();
  }
  return expanded;
}

// Compute steering syntax using the same particle names as the input parser
std::string MProcessTable::Syntax(const std::vector<MDecayBranch>& tree) const {
  std::string out;
  for (const auto& branch : tree) {
    if (!out.empty()) { out += " "; }
    out += branch.p.name;
    if (!branch.legs.empty()) { out += " > {" + Syntax(branch.legs) + "}"; }
  }
  return out;
}

// Compute final-state examples from process definitions and active model cards
std::vector<std::string> MProcessTable::Examples(const MProc& process, nuclear::CollisionType collision) const {
  std::vector<std::string> out;
  const auto& info = process.Info();
  if (!process.SupportsCollision(collision, general, pdg)) { return out; }
  if (info.example_decay == RootDecayMode::None) { return out; }
  const auto append = [&](const std::vector<MDecayBranch>& tree) {
    if (process.MatchProcess(tree).has_value()) { out.push_back(Syntax(tree)); }
  };
  for (const auto& state : process.Processes()) {
    if (!state.topology.empty()) { out.push_back(state.final_state_syntax); }
  }
  const bool pomeron = info.model != ReggeProductionModel::None;
  if (pomeron) {
    const auto& con =
        general.at(info.model == ReggeProductionModel::TP ? "PARAM_TENSORPOM" : "PARAM_REGGE").at("PARAM_CON").at(process.ISTATE);
    if (info.continuum && info.resonance == ResonanceType::None) {
      for (const auto& [key, channels] : con.items()) {
        if (key == "[*]" || key.front() != '[' || channels.empty()) { continue; }
        auto tree = Tree(key);
        append(tree);
        if (Expand(tree)) { append(tree); }
      }
    }
    if (info.resonance == ResonanceType::Required) {
      for (const auto& [name, resonance] : resonances) {
        if (resonance.at("PDG").get<int>() == 0) { continue; }
        const auto& models = resonance.at("MODELS");
        if (!models.contains(process.ISTATE) || !ActiveCoupling(models.at(process.ISTATE), coupling_min)) { continue; }
        if (collision != nuclear::CollisionType::PP) {
          MParticle particle;
          particle.pdg = resonance.at("PDG");
          particle.spinX2 = resonance.at("spinX2");
          particle.P = resonance.at("P");
          particle.C = resonance.at("C");
          bool supported = true;
          for (const auto& [key, channel] : models.at(process.ISTATE).items()) {
            if (key.front() != '[' || !ActiveCoupling(channel, coupling_min)) { continue; }
            const auto exchange = nlohmann::json::parse(key).get<std::array<int, 2>>();
            supported = supported && PhotoProductionError(info.model, collision, particle, exchange, pdg, general).empty();
          }
          if (!supported) { continue; }
        }
        const auto channels = decays.find(std::to_string(resonance.at("PDG").get<int>()));
        if (channels == decays.end()) { continue; }
        std::array<std::string, 2> preferred;
        std::array<double, 2> branching = {0.0, 0.0};
        const std::string suffix = " @RES{" + name + ":1}";
        // Prefer the largest configured branching ratio for direct and cascaded examples
        const auto select = [&](const auto& tree, std::size_t i, double br) {
          if (!process.MatchProcess(tree).has_value()) { return; }
          if (preferred[i].empty() || br > branching[i]) {
            preferred[i] = Syntax(tree) + suffix;
            branching[i] = br;
          }
        };
        for (const auto& [key, decay] : channels->items()) {
          if (key.front() != '[' || !ActiveCoupling(decay, coupling_min)) { continue; }
          if (info.continuum && !con.contains(key)) { continue; }
          auto tree = Tree(key);
          const double br = decay.value("BR", 0.0);
          select(tree, 0, br);
          if (!info.continuum && Expand(tree)) { select(tree, 1, br); }
        }
        for (const auto& example : preferred) { if (!example.empty()) { out.push_back(example); } }
      }
    }
  }
  if (!pomeron || info.resonance == ResonanceType::Optional) {
    LORENTZSCALAR lts;
    const int     root        = process.RootResonancePDG(lts);
    const auto    root_decays = decays.find(std::to_string(root));
    if (root != 0 && root_decays != decays.end()) {
      // Choose the largest configured branching ratio for a fixed resonance
      std::string preferred;
      double      br = 0.0;
      for (const auto& [key, decay] : root_decays->items()) {
        if (key.front() != '[' || !ActiveCoupling(decay, coupling_min)) { continue; }
        const auto tree = Tree(key);
        if (!process.MatchProcess(tree).has_value()) { continue; }
        const double next = decay.value("BR", 0.0);
        if (preferred.empty() || next > br) {
          preferred = Syntax(tree);
          br        = next;
        }
      }
      if (!preferred.empty()) { out.push_back(preferred); }
    } else if (root == PDG::PDG_monopolium) {
      // Use a diphoton kinematic example for isolated monopolium decays
      const auto gamma = pdg.FindByPDG(PDG::PDG_gamma).name;
      out.push_back(gamma + " " + gamma);
    } else if (root == 0) {
      for (const auto& [code, particle] : pdg.PDG_table) {
        if (code <= 0 || particle.name == "j") { continue; }
        const int conjugate = pdg.PDG_table.count(-code) ? -code : code;
        append(Tree("[" + std::to_string(code) + "," + std::to_string(conjugate) + "]"));
      }
    }
  }
  std::set<std::string> seen;
  std::erase_if(out, [&seen](const auto& value) { return !seen.insert(value).second; });
  return out;
}

// Compute compact overview rows in process registration order
std::vector<std::vector<std::string>> MProcessTable::Rows(const MSubProc& subprocess) const {
  auto processes = subprocess.CreateAllProcesses();
  std::sort(processes.begin(), processes.end(),
            [](const auto& a, const auto& b) { return a->DESCRIPTION.display_order < b->DESCRIPTION.display_order; });
  std::vector<std::vector<std::string>> rows;
  std::string istate;
  for (const auto& process : processes) {
    const auto& info = process->Info();
    const std::string prefix = process->ISTATE + "[" + process->CHANNEL + "]<";
    for (const auto& [command, descriptor] : subprocess.ProcessRegistry) {
      if (command.compare(0, prefix.size(), prefix) != 0) { continue; }
      const bool shared = subprocess.ProcessExist(prefix + "F>") && subprocess.ProcessExist(prefix + "C>");
      if (shared && command == prefix + "C>") { continue; }
      const std::string shown_command = shared ? prefix + "F|C>" : command;
      if (!istate.empty() && istate != process->ISTATE) { rows.emplace_back(); }
      istate = process->ISTATE;

      std::vector<std::pair<std::string, std::vector<std::string>>> groups;
      for (const auto collision : {nuclear::CollisionType::PP, nuclear::CollisionType::EE, nuclear::CollisionType::EP,
                                    nuclear::CollisionType::PA, nuclear::CollisionType::EA, nuclear::CollisionType::AA}) {
        if (!process->SupportsCollision(collision, general, pdg)) { continue; }
        auto examples = Examples(*process, collision);
        if (examples.empty() && info.example_decay != RootDecayMode::None) { continue; }
        // Prefer f mesons and compact syntax within each supported beam class
        std::stable_sort(examples.begin(), examples.end(), [](const auto& a, const auto& b) {
          return std::tuple{a.find("@RES{f") == std::string::npos, a.size()} <
                 std::tuple{b.find("@RES{f") == std::string::npos, b.size()};
        });
        const bool isolated = info.example_decay == RootDecayMode::Isolated;
        for (auto& example : examples) { example = (isolated ? "&> " : "-> ") + example; }
        std::vector<std::string> shown = {examples.empty() ? "-" : examples.front()};
        const auto nested = std::find_if(examples.begin(), examples.end(),
                                         [](const auto& example) { return example.find(" > {") != std::string::npos; });
        if (nested != examples.end() && *nested != shown.front()) { shown.push_back(*nested); }
        const auto group = std::find_if(groups.begin(), groups.end(), [&](const auto& entry) { return entry.second == shown; });
        if (group == groups.end()) {
          groups.emplace_back(nuclear::CollisionName(collision), std::move(shown));
        } else {
          group->first += ", " + nuclear::CollisionName(collision);
        }
      }
      for (const auto& [beams, examples] : groups) {
        rows.push_back({shown_command, beams, descriptor.model, descriptor.channels, StateSummary(*process, pdg), examples.front()});
        for (std::size_t i = 1; i < examples.size(); ++i) { rows.push_back({"", "", "", "", "", examples[i]}); }
      }
    }
  }
  return rows;
}

}  // namespace gra
