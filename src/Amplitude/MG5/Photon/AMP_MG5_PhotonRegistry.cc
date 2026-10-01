// Generated photon MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]
// [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181]

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_PhotonRegistry.h"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_ProcessRegistry.h"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <utility>

#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Helicity.h"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"

#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ll.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_uubarg.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_uubargg.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_jj.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_ww.h"
#include "Graniitti/Amplitude/MG5/Photon/AMP_MG5_yy_zjj.h"

namespace gra {
namespace {

// Compare an ordered final state without allocating a temporary vector
bool Matches(const std::vector<int> &actual,
             std::initializer_list<int> expected) {
  return actual.size() == expected.size() &&
         std::equal(actual.begin(), actual.end(), expected.begin());
}

// Compute true when one complete registry process matches these stable PDGs
bool ProcessMatches(const amplitude::Process &process,
                    const std::string &process_name,
                    std::initializer_list<int> stable_pdgs) {
  const auto expected =
      amplitude::FindProcess("PHOTON", process_name);
  return expected.has_value() && process == *expected &&
         Matches(process.stable_pdgs, stable_pdgs);
}

// Collect mutable stable leaves in decay-tree order
void CollectStableLeaves(MDecayBranch &branch,
                         std::vector<MDecayBranch *> &leaves) {
  if (branch.legs.empty()) {
    leaves.push_back(&branch);
    return;
  }
  for (auto &leg : branch.legs) { CollectStableLeaves(leg, leaves); }
}

// Convert one symbolic photon flow to stable-leaf event tags
bool PhotonColorCandidate(const mg5helas::ExternalColorFlow &external,
                          std::vector<MColorFlow> &candidate) {
  candidate.clear();
  if (external.size() < 2) {
    return false;
  }
  if (external[0].color != 0 || external[0].anticolor != 0 ||
      external[1].color != 0 || external[1].anticolor != 0) {
    return false;
  }

  std::map<int, int> tag_map;
  const auto convert_tag = [&tag_map](int tag) {
    if (tag == 0) { return 0; }
    const int sign = tag > 0 ? 1 : -1;
    const auto [entry, inserted] = tag_map.try_emplace(std::abs(tag), 501 + static_cast<int>(tag_map.size()));
    return sign * entry->second;
  };

  candidate.reserve(external.size() - 2);
  for (std::size_t leg = 2; leg < external.size(); ++leg) {
    candidate.push_back({convert_tag(external[leg].color),
                         convert_tag(external[leg].anticolor)});
  }

  std::map<int, std::array<std::size_t, 2>> tag_counts;
  for (const auto &leg : candidate) {
    if (leg.flow1 != 0) {
      ++tag_counts[std::abs(leg.flow1)][leg.flow1 < 0 ? 1 : 0];
    }
    if (leg.flow2 != 0) {
      ++tag_counts[std::abs(leg.flow2)][leg.flow2 > 0 ? 1 : 0];
    }
  }
  for (const auto &entry : tag_counts) {
    if (entry.second[0] != 1 || entry.second[1] != 1) {
      candidate.clear();
      return false;
    }
  }
  return true;
}

// Build structured helicity components from one projected color block
std::vector<mg5helas::HelicityComponent> BuildComponents(
    const std::vector<std::complex<double>> &projected,
    std::size_t row_offset, std::size_t row_count,
    const std::vector<int> &helicities, std::size_t helicity_count,
    std::size_t external_count, double normalization) {
  std::vector<mg5helas::HelicityComponent> components;
  components.reserve(row_count * helicity_count);
  for (std::size_t color = 0; color < row_count; ++color) {
    for (std::size_t ihel = 0; ihel < helicity_count; ++ihel) {
      const std::size_t helicity_offset = ihel * external_count;
      std::vector<int> outgoing(
          helicities.begin() + helicity_offset + 2,
          helicities.begin() + helicity_offset + external_count);
      components.push_back({
          {helicities[helicity_offset], helicities[helicity_offset + 1]},
          std::move(outgoing), color,
          normalization *
              projected[(row_offset + color) * helicity_count + ihel]});
    }
  }
  return components;
}

// Represent one concrete generated matrix element in one photon process
template <class MatrixElement>
class GeneratedPhotonProcess final : public PhotonMG5Process {
 public:
  // Store immutable color, helicity and shower-flow process_data
  GeneratedPhotonProcess(
      std::string name, std::vector<int> final_pdgs,
      std::vector<int> final_color_representations,
      std::size_t rank, int process_denominator,
      int final_symmetry_factor, int initial_state_denominator,
      double default_alpha_s, double default_alpha_qed, int alpha_s_power,
      int alpha_qed_power,
      std::vector<int> helicities,
      std::vector<std::complex<double>> combined_projectors,
      std::vector<std::vector<MColorFlow>> flow_candidates)
      : PhotonMG5Process(amplitude::Processes("PHOTON", name), default_alpha_qed),
        name_(std::move(name)),
        final_pdgs_(std::move(final_pdgs)),
        final_color_representations_(std::move(final_color_representations)),
        rank_(rank),
        process_denominator_(process_denominator),
        final_symmetry_factor_(final_symmetry_factor),
        initial_state_denominator_(initial_state_denominator),
        default_alpha_s_(default_alpha_s),
        default_alpha_qed_(default_alpha_qed),
        alpha_s_power_(alpha_s_power),
        alpha_qed_power_(alpha_qed_power),
        helicities_(std::move(helicities)),
        combined_projectors_(std::move(combined_projectors)),
        flow_candidates_(std::move(flow_candidates)) {
    const std::size_t external_count = final_pdgs_.size() + 2;
    const std::size_t combined_rows = rank_ * (1 + flow_candidates_.size());
    if (rank_ == 0 || process_denominator_ <= 0 ||
        final_symmetry_factor_ <= 0 || initial_state_denominator_ <= 0 ||
        final_symmetry_factor_ * initial_state_denominator_ !=
            process_denominator_ ||
        std::abs(
            mg5helas::FinalStateSymmetryFactor(final_pdgs_) -
            static_cast<double>(final_symmetry_factor_)) > 0.5 ||
        !std::isfinite(default_alpha_s_) || default_alpha_s_ < 0.0 ||
        !std::isfinite(default_alpha_qed_) || !(default_alpha_qed_ > 0.0) ||
        alpha_s_power_ < 0 || alpha_qed_power_ < 0 ||
        final_color_representations_.size() != final_pdgs_.size() ||
        helicities_.size() != MatrixElement::nhelicity * external_count ||
        combined_projectors_.size() != combined_rows * MatrixElement::ncolor ||
        !gra::AllFinite(combined_projectors_)) {
      throw std::invalid_argument(
          "GeneratedPhotonProcess: inconsistent generated data for " + name_);
    }
    projected_amplitudes_.resize(combined_rows * MatrixElement::nhelicity);
  }

  // Compute the one generated matrix element owned by this exact process
  std::size_t SubprocessCount() const override { return 1; }

  // Initialize all model parameters from one complete SLHA card
  void InitParameters(SLHAReader card) override { matrix_element_.InitParameters(std::move(card)); }

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override { return matrix_element_.Particles(); }

  // Compute the evaluated UFO electromagnetic coupling
  double DefaultAlphaQED() const override { return matrix_element_.AlphaQED(); }

  // Compute the Born matrix-element power of alpha_s
  int AlphaSPower() const noexcept override { return alpha_s_power_; }

  // Compute the Born matrix-element power of alpha_QED
  int AlphaQEDPower() const noexcept override { return alpha_qed_power_; }

  // Evaluate exact and shower-partitioned photon amplitudes in one HELAS pass
  mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) override {
    mg5helas::MatrixElementEvaluation result;
    amplitudes_.clear();
    flow_amplitudes_.clear();
    lts.epa_hard.Clear();
    lts.hard_color_flows.clear();
    const double effective_alpha_s =
        alpha_s <= 0.0 ? default_alpha_s_ : alpha_s;
    const std::size_t combined_rows = rank_ * (1 + flow_candidates_.size());
    std::fill(projected_amplitudes_.begin(), projected_amplitudes_.end(), 0.0);
    M4Vec k1;
    M4Vec k2;
    bool epa_hard = false;
    mg5helas::EPAHardFrame epa_frame;
    if (coherent_epa && mg5helas::HasEPAHardBeamState(lts)) {
      result.status = matrix_element_.CalcColorProjectedEPAHardHelicity(
          lts, effective_alpha_s, combined_projectors_.data(),
          static_cast<int>(combined_rows), projected_amplitudes_.data(),
          &epa_frame);
      k1 = epa_frame.incoming[0];
      k2 = epa_frame.incoming[1];
      epa_hard = result.Valid();
    } else {
      result.status = matrix_element_.CalcColorProjectedHelicity(
          lts, effective_alpha_s, combined_projectors_.data(),
          static_cast<int>(combined_rows), projected_amplitudes_.data(), &k1,
          &k2);
    }
    if (!result.Valid()) {
      lts.hamp.clear();
      return result;
    }

    const std::size_t external_count = final_pdgs_.size() + 2;
    const double normalization = std::sqrt(
        mg5helas::AppliedFinalStateSymmetryFactor(lts) /
        static_cast<double>(process_denominator_));
    const auto exact_components = BuildComponents(
        projected_amplitudes_, 0, rank_, helicities_, MatrixElement::nhelicity,
        external_count, normalization);
    auto amplitudes = mg5::PhotonAmplitudes(
        lts, exact_components, k1, k2, coherent_epa,
        epa_hard ? &epa_frame : nullptr);

    std::vector<std::vector<std::complex<double>>> flow_amplitudes(
        flow_candidates_.size());
    std::vector<std::vector<mg5helas::HelicityComponent>> hard_flow_components;
    hard_flow_components.reserve(flow_candidates_.size());
    for (std::size_t flow = 0; flow < flow_candidates_.size(); ++flow) {
      auto flow_components = BuildComponents(
          projected_amplitudes_, rank_ * (flow + 1), rank_,
          helicities_, MatrixElement::nhelicity, external_count, normalization);
      flow_amplitudes[flow] = mg5::PhotonAmplitudes(
          lts, flow_components, k1, k2, coherent_epa,
          epa_hard ? &epa_frame : nullptr);
      if (epa_hard) {
        hard_flow_components.push_back(std::move(flow_components));
      }
    }
    result.amp2 = gra::SquaredNorm(amplitudes);
    if (amplitudes.empty() ||
        !std::all_of(flow_amplitudes.cbegin(), flow_amplitudes.cend(),
                     [](const auto &flow) {
                       return !flow.empty();
                     })) {
      lts.hamp.clear();
      result.amp2 = 0.0;
      result.status = mg5helas::EvaluationStatus::AmplitudeFailure;
      return result;
    }

    amplitudes_ = std::move(amplitudes);
    flow_amplitudes_ = std::move(flow_amplitudes);
    for (const auto &flow : indices(flow_amplitudes_)) {
      mg5helas::ExternalColorFlow external(2);
      for (const auto &leg : flow_candidates_[flow]) { external.push_back({leg.flow1, leg.flow2}); }
      lts.hard_color_flows.push_back({flow_amplitudes_[flow], std::move(external)});
    }
    lts.hamp = amplitudes_;
    // Screening applies the incoming photon spin average after coherent summation
    gra::Scale(lts.hamp,
               std::sqrt(static_cast<double>(initial_state_denominator_)));
    if (alpha_s_power_ > 0) {
      lts.id1 = PDG::PDG_gamma;
      lts.id2 = PDG::PDG_gamma;
      lts.alphaQCD = effective_alpha_s;
    }
    if (epa_hard) {
      lts.epa_hard.frame         = epa_frame;
      lts.epa_hard.amplitude     = exact_components;
      lts.epa_hard.normalization = std::sqrt(static_cast<double>(initial_state_denominator_));
      for (const auto &flow : indices(hard_flow_components)) {
        mg5helas::ExternalColorFlow external(2);
        for (const auto &leg : flow_candidates_[flow]) { external.push_back({leg.flow1, leg.flow2}); }
        lts.epa_hard.color.push_back({std::move(hard_flow_components[flow]), std::move(external)});
      }
    }
    return result;
  }

  // Sample one generated shower-flow candidate using the evaluated amplitudes
  bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) override {
    if (flow_candidates_.empty()) { return ClearHardColorFlow(lts); }
    return SampleHardColorFlow(lts, random, lts.hard_color_flows, flow_candidates_.size());
  }

  // Compute whether this exact process has colored stable final states
  bool HasFinalStateColor() const override {
    return std::any_of(
        final_color_representations_.begin(),
        final_color_representations_.end(),
        [](int representation) { return representation != 1; });
  }

 private:
  MatrixElement matrix_element_;
  std::string name_;
  std::vector<int> final_pdgs_;
  std::vector<int> final_color_representations_;
  std::size_t rank_;
  int process_denominator_;
  int final_symmetry_factor_;
  int initial_state_denominator_;
  double default_alpha_s_;
  double default_alpha_qed_;
  int alpha_s_power_;
  int alpha_qed_power_;
  std::vector<int> helicities_;
  std::vector<std::complex<double>> combined_projectors_;
  std::vector<std::vector<MColorFlow>> flow_candidates_;
  std::vector<std::complex<double>> projected_amplitudes_;
  std::vector<std::complex<double>> amplitudes_;
  std::vector<std::vector<std::complex<double>>> flow_amplitudes_;
};

}  // namespace

// Assign one shower-flow candidate to stable final-state leaves
bool AssignHardColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate) {
  if (DecayTreeHasColoredIntermediate(lts.decaytree)) {
    return false;
  }
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : lts.decaytree) { CollectStableLeaves(branch, leaves); }
  if (candidate.size() != leaves.size()) { return false; }

  for (auto &branch : lts.decaytree) {
    ClearDecayBranchColorFlow(branch);
  }
  for (std::size_t i = 0; i < leaves.size(); ++i) {
    leaves[i]->p.color_flow = candidate[i];
  }
  return true;
}

// Assign the unique color-singlet flow of one quark-antiquark pair
bool AssignPhotonSingletQuarkPairColorFlow(LORENTZSCALAR &lts) {
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : lts.decaytree) { CollectStableLeaves(branch, leaves); }
  if (leaves.size() != 2) { return false; }

  const int first = leaves[0]->p.pdg;
  const int second = leaves[1]->p.pdg;
  if (first == 0 || std::abs(first) > 6 || first != -second) {
    return false;
  }

  std::vector<MColorFlow> candidate;
  candidate.reserve(2);
  for (const auto *leaf : leaves) {
    candidate.push_back(
        leaf->p.pdg > 0 ? MColorFlow{501, 0} : MColorFlow{0, 501});
  }
  return AssignHardColorFlowCandidate(lts, candidate);
}

// Clear color tags from all stable photon-process final states
bool ClearHardColorFlow(LORENTZSCALAR &lts) {
  for (auto &branch : lts.decaytree) {
    ClearDecayBranchColorFlow(branch);
  }
  return true;
}

// Sample one prepared hard process color flow
bool SampleHardColorFlow(
    LORENTZSCALAR &lts, MRandom &random,
    const std::vector<mg5helas::HardColorFlow> &flows,
    std::size_t expected_flows) {
  if (flows.empty() || flows.size() != expected_flows) {
    return false;
  }
  std::vector<double> weights(flows.size(), 0.0);
  for (std::size_t flow = 0; flow < flows.size(); ++flow) {
    weights[flow] = SquaredNorm(flows[flow].amplitudes);
    if (!std::isfinite(weights[flow]) || weights[flow] <= 0.0) {
      weights[flow] = 0.0;
    }
  }
  const double total = Sum(weights);
  if (!(std::isfinite(total) && total > 0.0)) {
    return false;
  }
  const double target = random.U(0.0, total);
  double cumulative = 0.0;
  std::size_t selected = weights.size() - 1;
  for (std::size_t flow = 0; flow < weights.size(); ++flow) {
    cumulative += weights[flow];
    if (target < cumulative) {
      selected = flow;
      break;
    }
  }
  std::vector<MColorFlow> candidate;
  return PhotonColorCandidate(flows[selected].external, candidate) &&
         AssignHardColorFlowCandidate(lts, candidate);
}

// Compute true when a generated matrix element matches this exact process
bool HasPhotonMG5Process(const amplitude::Process &process) {
  if (ProcessMatches(process, "yy_ll", {-11, 11})) { return true; }
  if (ProcessMatches(process, "yy_uubarg", {2, -2, 21})) { return true; }
  if (ProcessMatches(process, "yy_uubargg", {2, -2, 21, 21})) { return true; }
  return false;
}

// Construct the generated matrix element for this exact process
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const amplitude::Process &process) {
  if (ProcessMatches(process, "yy_ll", {-11, 11})) {
    return std::make_unique<GeneratedPhotonProcess<AMP_MG5_yy_ll>>(
        "yy_ll", std::vector<int>{-11, 11},
        std::vector<int>{1, 1},
        1, 4,
        1,
        4,
        0.11799999999999999,
        0.0075467711139788835,
        0,
        2,
        std::vector<int>{
      -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, -1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1
  },
        std::vector<std::complex<double>>{
      std::complex<double>(1, 0.0)
  },
        std::vector<std::vector<MColorFlow>>{});
  }
  if (ProcessMatches(process, "yy_uubarg", {2, -2, 21})) {
    return std::make_unique<GeneratedPhotonProcess<AMP_MG5_yy_uubarg>>(
        "yy_uubarg", std::vector<int>{2, -2, 21},
        std::vector<int>{3, -3, 8},
        1, 4,
        1,
        4,
        0.11799999999999999,
        0.0075467711139788835,
        1,
        2,
        std::vector<int>{
      -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, -1, -1, -1, -1, 1, -1, -1, 1, -1, 1, -1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1, -1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, -1, -1, 1, 1, 1, -1, 1, -1, -1, 1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1
  },
        std::vector<std::complex<double>>{
      std::complex<double>(2, 0.0),
      std::complex<double>(2, 0.0)
  },
        std::vector<std::vector<MColorFlow>>{
      {MColorFlow{501, 0}, MColorFlow{0, 502}, MColorFlow{502, 501}}
  });
  }
  if (ProcessMatches(process, "yy_uubargg", {2, -2, 21, 21})) {
    return std::make_unique<GeneratedPhotonProcess<AMP_MG5_yy_uubargg>>(
        "yy_uubargg", std::vector<int>{2, -2, 21, 21},
        std::vector<int>{3, -3, 8, 8},
        2, 8,
        2,
        4,
        0.11799999999999999,
        0.0075467711139788835,
        2,
        2,
        std::vector<int>{
      -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, -1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, -1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, -1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, -1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, -1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, -1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, -1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, -1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, -1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, -1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1, -1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, -1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, -1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, -1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, -1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, -1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1, -1, 1, 1, 1, 1, 1, 1
  },
        std::vector<std::complex<double>>{
      std::complex<double>(2.3094010767585029, 0.0),
      std::complex<double>(-0.28867513459481287, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(2.2912878474779199, 0.0),
      std::complex<double>(2.3094010767585029, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(-0.28867513459481287, 0.0),
      std::complex<double>(0.0, 0.0),
      std::complex<double>(2.2912878474779199, 0.0)
  },
        std::vector<std::vector<MColorFlow>>{
      {MColorFlow{501, 0}, MColorFlow{0, 503}, MColorFlow{502, 501}, MColorFlow{503, 502}},
      {MColorFlow{501, 0}, MColorFlow{0, 502}, MColorFlow{502, 503}, MColorFlow{503, 501}}
  });
  }
  return nullptr;
}

// Compute generated photon families and their public process channels
std::vector<MG5ProcessInfo> PhotonMG5ProcessInfos() {
  return {
      MG5ProcessInfo{"MG5_YY_JJ", "jj", "MG5cards/Photon/MG5_YY_JJ/param_card.dat"},
      MG5ProcessInfo{"MG5_YY_WW", "WW", "MG5cards/Photon/MG5_YY_WW/param_card.dat"},
      MG5ProcessInfo{"MG5_YY_ZJJ", "Zjj", "MG5cards/Photon/MG5_YY_ZJJ/param_card.dat"}
  };
}

// Construct one generated photon subprocess family
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const std::string &process_family) {
  if (process_family == "MG5_YY_JJ") {
    return std::make_unique<AMP_MG5_yy_jj>();
  }
  if (process_family == "MG5_YY_WW") {
    return std::make_unique<AMP_MG5_yy_ww>();
  }
  if (process_family == "MG5_YY_ZJJ") {
    return std::make_unique<AMP_MG5_yy_zjj>();
  }
  return nullptr;
}

}  // namespace gra
