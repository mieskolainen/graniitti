// Process configuration and worker-local mutable state operations
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

// C++
#include <algorithm>
#include <cmath>
#include <cstddef>
#include <functional>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

// Own
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Process/MProcessState.h"

namespace gra {

namespace mg5helas {

// Sample one finite hard color-flow weight, treating one eligible assignment deterministically
std::optional<std::size_t> SelectHardColorFlow(const std::vector<HardColorFlow> &flow, MRandom &random,
                                               const std::optional<std::size_t> channel) {
  std::size_t eligible = 0;
  std::size_t unique   = 0;
  double      total    = 0.0;
  for (std::size_t i = 0; i < flow.size(); ++i) {
    if (channel.has_value() && flow[i].channel != channel) { continue; }
    ++eligible;
    unique              = i;
    const double weight = flow[i].Weight();
    if (std::isfinite(weight) && weight > 0.0) { total += weight; }
  }
  if (eligible == 1) { return unique; }
  if (!std::isfinite(total) || !(total > 0.0)) { return std::nullopt; }

  const double target = random.U(0.0, total);
  double       sum    = 0.0;
  std::size_t  last   = 0;
  for (std::size_t i = 0; i < flow.size(); ++i) {
    if (channel.has_value() && flow[i].channel != channel) { continue; }
    const double weight = flow[i].Weight();
    if (!std::isfinite(weight) || !(weight > 0.0)) { continue; }
    last = i;
    sum += weight;
    if (target < sum) { return i; }
  }
  return last;
}

}  // namespace mg5helas

namespace {

// Compute true when one selector entry matches a candidate PDG id
bool SelectorEntryMatches(const FIDPDGCUT &cut, const std::size_t index,
                          const int candidate) {
  // The generic jet selector acts on resolved quarks, antiquarks and gluons
  if (cut.pdg[index] == PDG::PDG_hard_jet) {
    return candidate == PDG::PDG_gluon || (std::abs(candidate) >= 1 && std::abs(candidate) <= 5);
  }
  if (cut.pdg_abs[index]) {
    return std::abs(static_cast<long long>(candidate)) == cut.pdg[index];
  }
  return candidate == cut.pdg[index];
}

// Compute true when one candidate PDG id is present in a selector
bool SelectorContains(const FIDPDGCUT &cut, const int candidate) {
  for (std::size_t i = 0; i < cut.pdg.size(); ++i) {
    if (SelectorEntryMatches(cut, i, candidate)) {
      return true;
    }
  }
  return false;
}

// Collect stable leaves from one decay branch recursively
void CollectStableLeaves(const MDecayBranch &branch,
                         std::vector<const MDecayBranch *> &leaves) {
  if (branch.legs.empty()) {
    leaves.push_back(&branch);
    return;
  }
  for (const auto &leg : branch.legs) {
    CollectStableLeaves(leg, leaves);
  }
}

// Collect stable leaves from the complete central decay tree
std::vector<const MDecayBranch *>
CollectStableLeaves(const std::vector<MDecayBranch> &tree) {
  std::vector<const MDecayBranch *> leaves;
  for (const auto &branch : tree) {
    CollectStableLeaves(branch, leaves);
  }
  return leaves;
}

// Compute stable leaves whose PDG id is present in one selector
std::vector<const MDecayBranch *>
SelectLeaves(const std::vector<const MDecayBranch *> &leaves,
             const FIDPDGCUT &cut) {
  std::vector<const MDecayBranch *> selected;
  for (const auto *leaf : leaves) {
    if (SelectorContains(cut, leaf->p.pdg)) {
      selected.push_back(leaf);
    }
  }
  return selected;
}

// Compute the four-momentum sum of selected stable leaves
M4Vec LeafSystemMomentum(const std::vector<const MDecayBranch *> &leaves) {
  M4Vec system;
  for (const auto *leaf : leaves) {
    system += leaf->p4;
  }
  return system;
}

// Compute true when selected leaves match every selector entry exactly once
bool SelectorMultisetMatches(
    const FIDPDGCUT &cut, const std::vector<const MDecayBranch *> &selected) {
  if (selected.size() != cut.pdg.size()) {
    return false;
  }

  std::vector<bool> used(selected.size(), false);
  std::function<bool(std::size_t)> match = [&](const std::size_t selector) {
    if (selector == cut.pdg.size()) {
      return true;
    }
    for (std::size_t leaf = 0; leaf < selected.size(); ++leaf) {
      if (used[leaf] ||
          !SelectorEntryMatches(cut, selector, selected[leaf]->p.pdg)) {
        continue;
      }
      used[leaf] = true;
      if (match(selector + 1)) {
        return true;
      }
      used[leaf] = false;
    }
    return false;
  };
  return match(0);
}

// Compute true when one four-momentum passes a selected observable cut
bool PassSelectedObservables(const FIDPDGCUT &cut, const M4Vec &system) {
  return cut.M.Contains(system.M()) && cut.Rap.Contains(system.Rap()) &&
         cut.Eta.Contains(system.Eta()) && cut.Pt.Contains(system.Pt()) &&
         cut.Et.Contains(system.Et());
}

// Compute true when one PDG-selected fiducial cut passes
bool PassSelectedCut(const FIDPDGCUT &cut,
                     const std::vector<const MDecayBranch *> &leaves) {
  if (cut.pdg.size() != cut.pdg_abs.size()) {
    return false;
  }
  if (cut.pdg.empty()) {
    return true;
  }

  const auto selected = SelectLeaves(leaves, cut);
  if (selected.empty()) {
    return false;
  }
  if (cut.pdg.size() == 1) {
    return std::all_of(selected.begin(), selected.end(),
                       [&cut](const MDecayBranch *leaf) {
                         return PassSelectedObservables(cut, leaf->p4);
                       });
  }
  return SelectorMultisetMatches(cut, selected) &&
         PassSelectedObservables(cut, LeafSystemMomentum(selected));
}

// Compute true when one stable central leaf passes the common particle cuts
bool PassCentralLeaf(const FIDCUT &cuts, const MDecayBranch &branch) {
  if (branch.legs.empty()) {
    return (!cuts.particle_pt_active ||
            (branch.p4.Pt() >= cuts.pt_min && branch.p4.Pt() <= cuts.pt_max)) &&
           (!cuts.particle_Et_active ||
            (branch.p4.Et() >= cuts.Et_min && branch.p4.Et() <= cuts.Et_max)) &&
           (!cuts.particle_eta_active || (branch.p4.Eta() >= cuts.eta_min &&
                                          branch.p4.Eta() <= cuts.eta_max)) &&
           (!cuts.particle_rap_active || (branch.p4.Rap() >= cuts.rap_min &&
                                          branch.p4.Rap() <= cuts.rap_max));
  }
  return std::all_of(
      branch.legs.begin(), branch.legs.end(),
      [&cuts](const MDecayBranch &leg) { return PassCentralLeaf(cuts, leg); });
}

// Compute true when one stable leaf avoids all selected veto domains
bool PassVetoLeaf(const MDecayBranch &branch, const bool source_forward,
                  const std::vector<VETODOMAIN> &domains) {
  if (!branch.legs.empty()) {
    return std::all_of(branch.legs.begin(), branch.legs.end(),
                       [source_forward, &domains](const MDecayBranch &leg) {
                         return PassVetoLeaf(leg, source_forward, domains);
                       });
  }

  for (const auto &domain : domains) {
    const bool source_selected =
        source_forward ? domain.source_forward : domain.source_central;
    bool charge_selected = domain.charge == VetoCharge::Any;
    if (domain.charge == VetoCharge::Charged) {
      charge_selected = branch.p.chargeX3 != 0;
    } else if (domain.charge == VetoCharge::Neutral) {
      charge_selected = branch.p.chargeX3 == 0;
    }
    if (source_selected && charge_selected && branch.p4.Pt() >= domain.pt_min &&
        branch.p4.Pt() <= domain.pt_max && branch.p4.Eta() >= domain.eta_min &&
        branch.p4.Eta() <= domain.eta_max) {
      return false;
    }
  }
  return true;
}

} // namespace

// Compute true when one value is inside this active inclusive range
bool FIDCUTRANGE::Contains(const double value) const {
  return !active || (value >= min && value <= max);
}

// Compute this PDG selector in steering-card syntax
std::string FIDPDGCUT::SelectorString() const {
  if (pdg.size() != pdg_abs.size()) {
    throw std::invalid_argument(
        "FIDPDGCUT::SelectorString: malformed PDG selector state");
  }
  std::ostringstream output;
  output << "[";
  for (std::size_t i = 0; i < pdg.size(); ++i) {
    if (i != 0) {
      output << ", ";
    }
    const std::string pdg_name = pdg[i] == PDG::PDG_hard_jet ? "j" : std::to_string(pdg[i]);
    output << (pdg_abs[i] ? "ABS(" + pdg_name + ")" : pdg_name);
  }
  output << "]";
  return output.str();
}

// Compute true when the forward event kinematics pass all configured cuts
bool FIDCUT::PassForward(const LORENTZSCALAR &lts) const {
  if (!forward_t1.Contains(std::abs(lts.t1)) || !forward_t2.Contains(std::abs(lts.t2)) ||
      !PassForwardDeltaPhi(lts) || !PassForwardXi(lts)) {
    return false;
  }
  if (forward_t_active && !(std::abs(lts.t1) >= forward_t_min &&
                            std::abs(lts.t1) <= forward_t_max &&
                            std::abs(lts.t2) >= forward_t_min &&
                            std::abs(lts.t2) <= forward_t_max)) {
    return false;
  }
  if (forward_M_active && (lts.excite1 || lts.excite2) &&
      lts.pfinal.size() <= 2) {
    return false;
  }
  if (forward_M_active &&
      ((lts.excite1 && !(lts.pfinal[1].M() >= forward_M_min &&
                         lts.pfinal[1].M() <= forward_M_max)) ||
       (lts.excite2 && !(lts.pfinal[2].M() >= forward_M_min &&
                         lts.pfinal[2].M() <= forward_M_max)))) {
    return false;
  }
  return true;
}

// Compute true when the two forward legs pass the azimuthal-separation cut
bool FIDCUT::PassForwardDeltaPhi(const LORENTZSCALAR &lts) const {
  if (!forward_dPhi_active) {
    return true;
  }
  if (lts.pfinal.size() <= 2) {
    return false;
  }
  const double delta_phi =
      math::Rad2Deg(lts.pfinal[1].DeltaPhiAbs(lts.pfinal[2]));
  return delta_phi >= forward_dPhi_min && delta_phi <= forward_dPhi_max;
}

// Compute true when the available forward legs pass the momentum-loss cut
bool FIDCUT::PassForwardXi(const LORENTZSCALAR &lts) const {
  if (!forward_xi_active) {
    return true;
  }
  return (lts.has_xi1 || lts.has_xi2) &&
         (!lts.has_xi1 ||
          (lts.xi1 >= forward_xi_min && lts.xi1 <= forward_xi_max)) &&
         (!lts.has_xi2 ||
          (lts.xi2 >= forward_xi_min && lts.xi2 <= forward_xi_max));
}

// Compute true when the central system passes all configured cuts
bool FIDCUT::PassCentralSystem(const LORENTZSCALAR &lts) const {
  return (!system_M_active ||
          (std::isfinite(lts.m2) && lts.m2 >= 0.0 &&
           math::msqrt(lts.m2) >= M_min && math::msqrt(lts.m2) <= M_max)) &&
         (!system_Rap_active || (lts.Y >= Y_min && lts.Y <= Y_max)) &&
         (!system_Pt_active || (lts.Pt >= Pt_min && lts.Pt <= Pt_max));
}

// Compute true when every stable central particle passes the common cuts
bool FIDCUT::PassCentralParticles(const std::vector<MDecayBranch> &tree) const {
  return std::all_of(tree.begin(), tree.end(),
                     [this](const MDecayBranch &branch) {
                       return PassCentralLeaf(*this, branch);
                     });
}

// Compute true when all PDG-selected particle and system cuts pass
bool FIDCUT::PassSelectedParticles(
    const std::vector<MDecayBranch> &tree) const {
  if (pdg_cuts.empty()) {
    return true;
  }
  const auto leaves = CollectStableLeaves(tree);
  return std::all_of(
      pdg_cuts.begin(), pdg_cuts.end(),
      [&leaves](const FIDPDGCUT &cut) { return PassSelectedCut(cut, leaves); });
}

// Compute true when no stable forward or central particle enters a veto domain
bool VETOCUT::Pass(const LORENTZSCALAR &lts) const {
  return Pass(lts, lts.decaytree);
}

// Test veto domains with an explicit central momentum tree
bool VETOCUT::Pass(const LORENTZSCALAR &lts,
                   const std::vector<MDecayBranch> &central) const {
  if (!active) {
    return true;
  }
  if (!PassVetoLeaf(lts.decayforward1, true, cuts) ||
      !PassVetoLeaf(lts.decayforward2, true, cuts)) {
    return false;
  }
  return std::all_of(central.begin(), central.end(),
                     [this](const MDecayBranch &branch) {
                       return PassVetoLeaf(branch, false, cuts);
                     });
}

} // namespace gra
