#!/usr/bin/env python3
#
# Generate MG5 gamma-gamma metadata, amplitudes and runtime registries
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from pathlib import Path
from typing import Any

from core.io.serialize import load_json_file

from . import mg5_color, output_layout, process_registry
from .process_registry import parse_process_syntax, stable_tokens, topology_signature
from .registry_support import (
    combined_projectors,
    cpp_complex_matrix,
    cpp_float,
    cpp_flow_candidates,
    cpp_int_vector,
    exact_int,
    finite_float,
    int_vector,
    validate_color_data,
)


# Format one flattened generated helicity table
def cpp_helicities(table: list[list[int]]) -> str:
    values = [str(int(value)) for row in table for value in row]
    return "{\n      " + ", ".join(values) + "\n  }"


# Load and validate one checked-in photon process_data sidecar
def load_process_data(base_dir: Path, entry: dict[str, Any]) -> dict[str, Any]:
    path = base_dir / output_layout.color_path(entry["projection"], entry["name"])
    if not path.is_file():
        raise RuntimeError(f"Missing generated photon process_data: {path}")
    process_data = load_json_file(path)
    validate_color_data(process_data, entry, [22, 22], 1, "photon")
    required = {
        "has_decay_chain",
        "process_denominator",
        "final_symmetry_factor",
        "initial_state_denominator",
        "default_alpha_s",
        "default_alpha_qed",
        "coupling_orders",
        "alpha_s_power",
        "alpha_qed_power",
        "color_metric_residual",
        "helicities",
    }
    missing = required - set(process_data)
    if missing:
        raise RuntimeError(f"Photon process_data for {entry['name']} lacks {sorted(missing)}")
    parsed = parse_process_syntax(entry["process"])
    if len(stable_tokens(parsed)) != len(process_data["final_pdgs"]):
        raise RuntimeError(f"Invalid photon stable final-state size for {entry['name']}")
    if (
        type(process_data["has_decay_chain"]) is not bool
        or process_data["has_decay_chain"] != parsed.has_decay_chain
    ):
        raise RuntimeError(f"Invalid photon decay-chain flag for {entry['name']}")
    denominator = exact_int(process_data, "process_denominator", "photon", positive=True)
    symmetry = exact_int(process_data, "final_symmetry_factor", "photon", positive=True)
    initial = exact_int(process_data, "initial_state_denominator", "photon", positive=True)
    alpha_s_power = exact_int(process_data, "alpha_s_power", "photon")
    alpha_qed_power = exact_int(process_data, "alpha_qed_power", "photon")
    orders = process_data["coupling_orders"]
    if not isinstance(orders, list) or not orders or any(
        not isinstance(row, dict) or any(type(value) is not int or value < 0 for value in row.values())
        for row in orders
    ):
        raise RuntimeError(f"Invalid photon coupling orders for {entry['name']}")
    for name, power in (("QCD", alpha_s_power), ("QED", alpha_qed_power)):
        if power != max(row.get(name, 0) for row in orders):
            raise RuntimeError(f"Invalid photon {name} power for {entry['name']}")
    if (
        symmetry != mg5_color.final_state_symmetry_factor(process_data["final_pdgs"])
        or symmetry * initial != denominator
    ):
        raise RuntimeError(f"Invalid photon denominator split for {entry['name']}")
    finite_float(
        process_data["default_alpha_s"],
        f"photon default alpha_s for {entry['name']}",
        nonnegative=True,
    )
    finite_float(
        process_data["default_alpha_qed"],
        f"photon default alpha_QED for {entry['name']}",
        nonnegative=True,
    )
    finite_float(
        process_data["color_metric_residual"],
        f"photon color-metric residual for {entry['name']}",
        nonnegative=True,
    )
    external_count = 2 + len(process_data["final_pdgs"])
    helicities = process_data["helicities"]
    if type(helicities) is not list or not helicities:
        raise RuntimeError(f"Invalid photon helicity table for {entry['name']}")
    for row in helicities:
        int_vector(
            row,
            f"photon helicity row for {entry['name']}",
            external_count,
        )
    return process_data


# Validate that every ordered topology and stable final state has one process
def validate_signatures(entries: list[tuple[dict[str, Any], dict[str, Any]]]) -> None:
    owners: dict[tuple[str, tuple[int, ...]], str] = {}
    for entry, process_data in entries:
        signature = (
            topology_signature(parse_process_syntax(entry["process"])),
            tuple(int(pdg) for pdg in process_data["final_pdgs"]),
        )
        if signature in owners:
            raise RuntimeError(
                "Duplicate photon topology and stable-final-state signature "
                f"{signature}: {owners[signature]} and {entry['name']}"
            )
        owners[signature] = entry["name"]


# Generate the photon process base class
def registry_header() -> str:
    return (
        """// Generated photon MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PHOTONREGISTRY_H
#define AMP_MG5_PHOTONREGISTRY_H

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "@HELICITY_INCLUDE@"
#include "@SUBPROCESS_SUM_INCLUDE@"
#include "@PROCESS_INCLUDE@"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {

namespace flux {
struct EPASectorWeights;
}

// Clear color tags from all stable photon-process final states
bool ClearHardColorFlow(LORENTZSCALAR &lts);

// Common runtime interface for every generated gamma-gamma amplitude
class PhotonMG5Process : public MG5Process {
 public:
  // Store the exact processes represented by this generated runtime
  PhotonMG5Process(std::vector<amplitude::Process> processes, double alpha_qed)
      : MG5Process(std::move(processes)), alpha_qed_(alpha_qed) {}

  // Destroy one generated photon runtime through its base class
  virtual ~PhotonMG5Process() = default;

  // Evaluate the generated photon helicity amplitudes
  virtual mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) = 0;

  // Compute the Born matrix-element power of alpha_s
  virtual int AlphaSPower() const noexcept = 0;

  // Compute the Born matrix-element power of alpha_QED
  virtual int AlphaQEDPower() const noexcept = 0;

  // Compute the reference electromagnetic coupling
  virtual double DefaultAlphaQED() const { return alpha_qed_; }

  // Sample one generated final-state color flow
  virtual bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) = 0;

  // Compute whether the selected event channel has colored stable final states
  virtual bool HasFinalStateColor() const = 0;

 private:
  double alpha_qed_ = 0.0;
};

// Assign one shower-flow candidate to stable final-state leaves
bool AssignHardColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate);

// Assign the unique color-singlet flow of one quark-antiquark pair
bool AssignPhotonSingletQuarkPairColorFlow(LORENTZSCALAR &lts);

// Sample one prepared generated photon color flow
bool SampleHardColorFlow(
    LORENTZSCALAR &lts, MRandom &random,
    const std::vector<mg5helas::HardColorFlow> &flows,
    std::size_t expected_flows);

// Compute true when a generated photon matrix element matches this exact process
bool HasPhotonMG5Process(const amplitude::Process &process);

// Construct the generated photon matrix element for this exact process
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const amplitude::Process &process);

// Compute generated photon families and their public process channels
std::vector<MG5ProcessInfo> PhotonMG5ProcessInfos();

// Construct one generated photon subprocess family
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const std::string &process_family);

}  // namespace gra

#endif
""".replace(
            "@HELICITY_INCLUDE@",
            output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Helicity.h"),
        )
        .replace(
            "@SUBPROCESS_SUM_INCLUDE@",
            output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_SubprocessSum.h"),
        )
        .replace(
            "@PROCESS_INCLUDE@",
            output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Process.h"),
        )
    )


# Generate one matrix element construction block for a photon process
def matrix_element_block(entry: dict[str, Any], process_data: dict[str, Any]) -> str:
    name = entry["name"]
    class_name = f"AMP_MG5_{name}"
    signature = cpp_int_vector(process_data["final_pdgs"])
    combined = combined_projectors(process_data)
    flows = cpp_flow_candidates(process_data["flow_candidates"])
    return f"""  if (ProcessMatches(process, "{name}", {signature})) {{
    return std::make_unique<GeneratedPhotonProcess<{class_name}>>(
        "{name}", std::vector<int>{signature},
        std::vector<int>{cpp_int_vector(process_data["final_color_representations"])},
        {int(process_data["rank"])}, {int(process_data["process_denominator"])},
        {int(process_data["final_symmetry_factor"])},
        {int(process_data["initial_state_denominator"])},
        {cpp_float(process_data["default_alpha_s"])},
        {cpp_float(process_data["default_alpha_qed"])},
        {int(process_data["alpha_s_power"])},
        {int(process_data["alpha_qed_power"])},
        std::vector<int>{cpp_helicities(process_data["helicities"])},
        std::vector<std::complex<double>>{cpp_complex_matrix(combined)},
        std::vector<std::vector<MColorFlow>>{flows});
  }}
"""


# Generate one exact process match block
def process_match_block(entry: dict[str, Any], process_data: dict[str, Any]) -> str:
    signature = cpp_int_vector(process_data["final_pdgs"])
    return f'  if (ProcessMatches(process, "{entry["name"]}", {signature})) {{ return true; }}\n'


# Generate the concrete registry implementation for all photon entries
def registry_source(
    entries: list[tuple[dict[str, Any], dict[str, Any]]],
    families: list[dict[str, Any]],
) -> str:
    standalone_includes = "\n".join(
        f'#include "{output_layout.include_path(output_layout.PHOTON, f"AMP_MG5_{entry['name']}.h")}"'
        for entry, _ in entries
    )
    family_includes = "\n".join(
        f'#include "{output_layout.include_path(output_layout.PHOTON, f"{process_registry.family_wrapper_name(family)}.h")}"'
        for family in families
    )
    includes = "\n".join(block for block in (standalone_includes, family_includes) if block)
    matches = "".join(process_match_block(entry, process_data) for entry, process_data in entries)
    matrix_elements = "".join(
        matrix_element_block(entry, process_data) for entry, process_data in entries
    )
    infos = ",\n".join(
        f'      MG5ProcessInfo{{"{family["name"]}", "{family["channel"]}", '
        f'"{output_layout.source_path(output_layout.parameter_card_path(family["projection"], family["name"]))}"}}'
        for family in families
    )
    factories = "".join(
        f'  if (process_family == "{family["name"]}") {{\n'
        f"    return std::make_unique<{process_registry.family_wrapper_name(family)}>();\n"
        "  }\n"
        for family in families
    )
    return f"""// Generated photon MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.
//
// [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]
// [REFERENCE: Budnev et al., Phys. Rept. 15 (1975) 181]

#include "{output_layout.include_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.h")}"
#include "{output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h")}"

#include <algorithm>
#include <array>
#include <cmath>
#include <initializer_list>
#include <map>
#include <stdexcept>
#include <utility>

#include "{output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Helicity.h")}"
#include "Graniitti/Math/MMath.h"
#include "Graniitti/Math/MMatrix.h"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Photon/MFlux.h"
#include "Graniitti/Photon/MQED.h"

{includes}

namespace gra {{
namespace {{

// Compare an ordered final state without allocating a temporary vector
bool Matches(const std::vector<int> &actual,
             std::initializer_list<int> expected) {{
  return actual.size() == expected.size() &&
         std::equal(actual.begin(), actual.end(), expected.begin());
}}

// Compute true when one complete registry process matches these stable PDGs
bool ProcessMatches(const amplitude::Process &process,
                    const std::string &process_name,
                    std::initializer_list<int> stable_pdgs) {{
  const auto expected =
      amplitude::FindProcess("PHOTON", process_name);
  return expected.has_value() && process == *expected &&
         Matches(process.stable_pdgs, stable_pdgs);
}}

// Collect mutable stable leaves in decay-tree order
void CollectStableLeaves(MDecayBranch &branch,
                         std::vector<MDecayBranch *> &leaves) {{
  if (branch.legs.empty()) {{
    leaves.push_back(&branch);
    return;
  }}
  for (auto &leg : branch.legs) {{ CollectStableLeaves(leg, leaves); }}
}}

// Convert one symbolic photon flow to stable-leaf event tags
bool PhotonColorCandidate(const mg5helas::ExternalColorFlow &external,
                          std::vector<MColorFlow> &candidate) {{
  candidate.clear();
  if (external.size() < 2) {{
    return false;
  }}
  if (external[0].color != 0 || external[0].anticolor != 0 ||
      external[1].color != 0 || external[1].anticolor != 0) {{
    return false;
  }}

  std::map<int, int> tag_map;
  const auto convert_tag = [&tag_map](int tag) {{
    if (tag == 0) {{ return 0; }}
    const int sign = tag > 0 ? 1 : -1;
    const auto [entry, inserted] = tag_map.try_emplace(std::abs(tag), 501 + static_cast<int>(tag_map.size()));
    return sign * entry->second;
  }};

  candidate.reserve(external.size() - 2);
  for (std::size_t leg = 2; leg < external.size(); ++leg) {{
    candidate.push_back({{convert_tag(external[leg].color),
                         convert_tag(external[leg].anticolor)}});
  }}

  std::map<int, std::array<std::size_t, 2>> tag_counts;
  for (const auto &leg : candidate) {{
    if (leg.flow1 != 0) {{
      ++tag_counts[std::abs(leg.flow1)][leg.flow1 < 0 ? 1 : 0];
    }}
    if (leg.flow2 != 0) {{
      ++tag_counts[std::abs(leg.flow2)][leg.flow2 > 0 ? 1 : 0];
    }}
  }}
  for (const auto &entry : tag_counts) {{
    if (entry.second[0] != 1 || entry.second[1] != 1) {{
      candidate.clear();
      return false;
    }}
  }}
  return true;
}}

// Build structured helicity components from one projected color block
std::vector<mg5helas::HelicityComponent> BuildComponents(
    const std::vector<std::complex<double>> &projected,
    std::size_t row_offset, std::size_t row_count,
    const std::vector<int> &helicities, std::size_t helicity_count,
    std::size_t external_count, double normalization) {{
  std::vector<mg5helas::HelicityComponent> components;
  components.reserve(row_count * helicity_count);
  for (std::size_t color = 0; color < row_count; ++color) {{
    for (std::size_t ihel = 0; ihel < helicity_count; ++ihel) {{
      const std::size_t helicity_offset = ihel * external_count;
      std::vector<int> outgoing(
          helicities.begin() + helicity_offset + 2,
          helicities.begin() + helicity_offset + external_count);
      components.push_back({{
          {{helicities[helicity_offset], helicities[helicity_offset + 1]}},
          std::move(outgoing), color,
          normalization *
              projected[(row_offset + color) * helicity_count + ihel]}});
    }}
  }}
  return components;
}}

// Represent one concrete generated matrix element in one photon process
template <class MatrixElement>
class GeneratedPhotonProcess final : public PhotonMG5Process {{
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
        flow_candidates_(std::move(flow_candidates)) {{
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
        !gra::AllFinite(combined_projectors_)) {{
      throw std::invalid_argument(
          "GeneratedPhotonProcess: inconsistent generated data for " + name_);
    }}
    projected_amplitudes_.resize(combined_rows * MatrixElement::nhelicity);
  }}

  // Compute the one generated matrix element owned by this exact process
  std::size_t SubprocessCount() const override {{ return 1; }}

  // Initialize all model parameters from one complete SLHA card
  void InitParameters(SLHAReader card) override {{ matrix_element_.InitParameters(std::move(card)); }}

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override {{ return matrix_element_.Particles(); }}

  // Compute the evaluated UFO electromagnetic coupling
  double DefaultAlphaQED() const override {{ return matrix_element_.AlphaQED(); }}

  // Compute the Born matrix-element power of alpha_s
  int AlphaSPower() const noexcept override {{ return alpha_s_power_; }}

  // Compute the Born matrix-element power of alpha_QED
  int AlphaQEDPower() const noexcept override {{ return alpha_qed_power_; }}

  // Evaluate exact and shower-partitioned photon amplitudes in one HELAS pass
  mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) override {{
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
    if (coherent_epa && mg5helas::HasEPAHardBeamState(lts)) {{
      result.status = matrix_element_.CalcColorProjectedEPAHardHelicity(
          lts, effective_alpha_s, combined_projectors_.data(),
          static_cast<int>(combined_rows), projected_amplitudes_.data(),
          &epa_frame);
      k1 = epa_frame.incoming[0];
      k2 = epa_frame.incoming[1];
      epa_hard = result.Valid();
    }} else {{
      result.status = matrix_element_.CalcColorProjectedHelicity(
          lts, effective_alpha_s, combined_projectors_.data(),
          static_cast<int>(combined_rows), projected_amplitudes_.data(), &k1,
          &k2);
    }}
    if (!result.Valid()) {{
      lts.hamp.clear();
      return result;
    }}

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
    for (std::size_t flow = 0; flow < flow_candidates_.size(); ++flow) {{
      auto flow_components = BuildComponents(
          projected_amplitudes_, rank_ * (flow + 1), rank_,
          helicities_, MatrixElement::nhelicity, external_count, normalization);
      flow_amplitudes[flow] = mg5::PhotonAmplitudes(
          lts, flow_components, k1, k2, coherent_epa,
          epa_hard ? &epa_frame : nullptr);
      if (epa_hard) {{
        hard_flow_components.push_back(std::move(flow_components));
      }}
    }}
    result.amp2 = gra::SquaredNorm(amplitudes);
    if (amplitudes.empty() ||
        !std::all_of(flow_amplitudes.cbegin(), flow_amplitudes.cend(),
                     [](const auto &flow) {{
                       return !flow.empty();
                     }})) {{
      lts.hamp.clear();
      result.amp2 = 0.0;
      result.status = mg5helas::EvaluationStatus::AmplitudeFailure;
      return result;
    }}

    amplitudes_ = std::move(amplitudes);
    flow_amplitudes_ = std::move(flow_amplitudes);
    for (const auto &flow : indices(flow_amplitudes_)) {{
      mg5helas::ExternalColorFlow external(2);
      for (const auto &leg : flow_candidates_[flow]) {{ external.push_back({{leg.flow1, leg.flow2}}); }}
      lts.hard_color_flows.push_back({{flow_amplitudes_[flow], std::move(external)}});
    }}
    lts.hamp = amplitudes_;
    // Screening applies the incoming photon spin average after coherent summation
    gra::Scale(lts.hamp,
               std::sqrt(static_cast<double>(initial_state_denominator_)));
    if (alpha_s_power_ > 0) {{
      lts.id1 = PDG::PDG_gamma;
      lts.id2 = PDG::PDG_gamma;
      lts.alphaQCD = effective_alpha_s;
    }}
    if (epa_hard) {{
      lts.epa_hard.frame         = epa_frame;
      lts.epa_hard.amplitude     = exact_components;
      lts.epa_hard.normalization = std::sqrt(static_cast<double>(initial_state_denominator_));
      for (const auto &flow : indices(hard_flow_components)) {{
        mg5helas::ExternalColorFlow external(2);
        for (const auto &leg : flow_candidates_[flow]) {{ external.push_back({{leg.flow1, leg.flow2}}); }}
        lts.epa_hard.color.push_back({{std::move(hard_flow_components[flow]), std::move(external)}});
      }}
    }}
    return result;
  }}

  // Sample one generated shower-flow candidate using the evaluated amplitudes
  bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) override {{
    if (flow_candidates_.empty()) {{ return ClearHardColorFlow(lts); }}
    return SampleHardColorFlow(lts, random, lts.hard_color_flows, flow_candidates_.size());
  }}

  // Compute whether this exact process has colored stable final states
  bool HasFinalStateColor() const override {{
    return std::any_of(
        final_color_representations_.begin(),
        final_color_representations_.end(),
        [](int representation) {{ return representation != 1; }});
  }}

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
}};

}}  // namespace

// Assign one shower-flow candidate to stable final-state leaves
bool AssignHardColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate) {{
  if (DecayTreeHasColoredIntermediate(lts.decaytree)) {{
    return false;
  }}
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : lts.decaytree) {{ CollectStableLeaves(branch, leaves); }}
  if (candidate.size() != leaves.size()) {{ return false; }}

  for (auto &branch : lts.decaytree) {{
    ClearDecayBranchColorFlow(branch);
  }}
  for (std::size_t i = 0; i < leaves.size(); ++i) {{
    leaves[i]->p.color_flow = candidate[i];
  }}
  return true;
}}

// Assign the unique color-singlet flow of one quark-antiquark pair
bool AssignPhotonSingletQuarkPairColorFlow(LORENTZSCALAR &lts) {{
  std::vector<MDecayBranch *> leaves;
  for (auto &branch : lts.decaytree) {{ CollectStableLeaves(branch, leaves); }}
  if (leaves.size() != 2) {{ return false; }}

  const int first = leaves[0]->p.pdg;
  const int second = leaves[1]->p.pdg;
  if (first == 0 || std::abs(first) > 6 || first != -second) {{
    return false;
  }}

  std::vector<MColorFlow> candidate;
  candidate.reserve(2);
  for (const auto *leaf : leaves) {{
    candidate.push_back(
        leaf->p.pdg > 0 ? MColorFlow{{501, 0}} : MColorFlow{{0, 501}});
  }}
  return AssignHardColorFlowCandidate(lts, candidate);
}}

// Clear color tags from all stable photon-process final states
bool ClearHardColorFlow(LORENTZSCALAR &lts) {{
  for (auto &branch : lts.decaytree) {{
    ClearDecayBranchColorFlow(branch);
  }}
  return true;
}}

// Sample one prepared hard process color flow
bool SampleHardColorFlow(
    LORENTZSCALAR &lts, MRandom &random,
    const std::vector<mg5helas::HardColorFlow> &flows,
    std::size_t expected_flows) {{
  if (flows.empty() || flows.size() != expected_flows) {{
    return false;
  }}
  std::vector<double> weights(flows.size(), 0.0);
  for (std::size_t flow = 0; flow < flows.size(); ++flow) {{
    weights[flow] = SquaredNorm(flows[flow].amplitudes);
    if (!std::isfinite(weights[flow]) || weights[flow] <= 0.0) {{
      weights[flow] = 0.0;
    }}
  }}
  const double total = Sum(weights);
  if (!(std::isfinite(total) && total > 0.0)) {{
    return false;
  }}
  const double target = random.U(0.0, total);
  double cumulative = 0.0;
  std::size_t selected = weights.size() - 1;
  for (std::size_t flow = 0; flow < weights.size(); ++flow) {{
    cumulative += weights[flow];
    if (target < cumulative) {{
      selected = flow;
      break;
    }}
  }}
  std::vector<MColorFlow> candidate;
  return PhotonColorCandidate(flows[selected].external, candidate) &&
         AssignHardColorFlowCandidate(lts, candidate);
}}

// Compute true when a generated matrix element matches this exact process
bool HasPhotonMG5Process(const amplitude::Process &process) {{
{matches}  return false;
}}

// Construct the generated matrix element for this exact process
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const amplitude::Process &process) {{
{matrix_elements}  return nullptr;
}}

// Compute generated photon families and their public process channels
std::vector<MG5ProcessInfo> PhotonMG5ProcessInfos() {{
  return {{
{infos}
  }};
}}

// Construct one generated photon subprocess family
std::unique_ptr<PhotonMG5Process> CreatePhotonMG5Process(
    const std::string &process_family) {{
{factories}  return nullptr;
}}

}}  // namespace gra
"""


# Load all photon process and family data and generate one registry
def generate_registry(base_dir: Path, manifest: dict[str, Any]) -> tuple[str, str]:
    entries = [
        (entry, load_process_data(base_dir, entry))
        for entry in manifest["processes"]
        if entry["projection"] == process_registry.PHOTON_PROJECTION
    ]
    families = photon_families(manifest)
    validate_signatures(entries)
    return registry_header(), registry_source(entries, families)


# Compute generated families evaluated with incoming photon helicities
def photon_families(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    families: list[dict[str, Any]] = []
    for family in manifest.get("families", []):
        if family["projection"] != process_registry.PHOTON_PROJECTION:
            continue
        process_registry.family_process_family(family)
        families.append(family)
    return families


# Generate one complete photon-family amplitude declaration
def family_amplitude_header(family: dict[str, Any]) -> str:
    family_name = family["name"]
    amplitude = process_registry.family_wrapper_name(family)
    channel = family["channel"]
    guard = f"{amplitude.upper()}_H"
    family_directory = output_layout.family_cpp_directory(family["projection"], family_name)
    registry_include = output_layout.include_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.h")
    subprocess_include = output_layout.include_path(
        output_layout.RUNTIME, "AMP_MG5_SubprocessSum.h"
    )
    process_base_include = output_layout.include_path(family_directory, "ProcessBase.h")
    return f"""// Generated MG5 amplitude for gamma-gamma {channel} production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef {guard}
#define {guard}

#include <complex>
#include <cstddef>
#include <vector>

#include "{registry_include}"
#include "{subprocess_include}"
#include "{process_base_include}"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Sampling/MRandom.h"

namespace gra {{

// Sum every generated subprocess for gamma-gamma {channel} production
class {amplitude} : public PhotonMG5Process {{
 public:
  // Construct all generated subprocesses
  {amplitude}();

  // Destroy all generated subprocesses
  ~{amplitude}() override = default;

  {amplitude}(const {amplitude} &) = delete;
  {amplitude} &operator=(const {amplitude} &) = delete;

  // Evaluate the summed matrix element squared
  mg5helas::MatrixElementEvaluation Evaluate(
      LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) override;

  // Compute the Born matrix-element power of alpha_s
  int AlphaSPower() const noexcept override;

  // Compute the Born matrix-element power of alpha_QED
  int AlphaQEDPower() const noexcept override;

  // Sample one generated final-state color flow
  bool SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) override;

  // Compute whether the selected event channel has colored stable final states
  bool HasFinalStateColor() const override;

  // Compute the number of generated subprocesses
  std::size_t SubprocessCount() const override;

  // Initialize all model parameters before sampling
  void InitParameters(SLHAReader card) override {{ subprocess_sum.InitParameters(std::move(card)); }}

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override {{ return subprocess_sum.Particles(); }}

  // Compute the evaluated UFO electromagnetic coupling
  double DefaultAlphaQED() const override {{ return subprocess_sum.AlphaQED(); }}

 private:
  mg5::SubprocessSum<{family_name}::ProcessBase> subprocess_sum;
  bool has_final_state_color_ = false;
}};

}}  // namespace gra

#endif
"""


# Generate one complete photon-family amplitude implementation
def family_amplitude_source(family: dict[str, Any]) -> str:
    family_name = family["name"]
    amplitude = process_registry.family_wrapper_name(family)
    channel = family["channel"]
    parameter_card = output_layout.source_path(
        output_layout.parameter_card_path(family["projection"], family_name)
    )
    family_directory = output_layout.family_cpp_directory(family["projection"], family_name)
    amplitude_include = output_layout.include_path(output_layout.PHOTON, f"{amplitude}.h")
    processes_include = output_layout.include_path(family_directory, "Processes.h")
    return f"""// Generated MG5 amplitude for gamma-gamma {channel} production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <utility>
#include <vector>

#include "{amplitude_include}"
#include "{processes_include}"
#include "Graniitti/Particle/MPDG.h"
#include "Graniitti/Tech/MAux.h"

namespace gra {{

// Construct all generated subprocesses
{amplitude}::{amplitude}()
    : PhotonMG5Process(amplitude::Processes("{family_name}"), 0.0) {{
  std::vector<mg5::Subprocess<{family_name}::ProcessBase>> subprocesses;
  {family_name}::BuildSubprocesses(
      subprocesses,
      aux::ResolveProjectPath("{parameter_card}"));
  mg5::ValidateSubprocessInitialStates(
      subprocesses, mg5::PhotonInitialStates(), "{amplitude}");
  subprocess_sum.SetSubprocesses(std::move(subprocesses));
}}

// Compute the Born matrix-element power of alpha_s
int {amplitude}::AlphaSPower() const noexcept {{ return subprocess_sum.AlphaSPower(); }}

// Compute the Born matrix-element power of alpha_QED
int {amplitude}::AlphaQEDPower() const noexcept {{ return subprocess_sum.AlphaQEDPower(); }}

// Evaluate the summed matrix element squared
mg5helas::MatrixElementEvaluation {amplitude}::Evaluate(
    LORENTZSCALAR &lts, double alpha_s, bool coherent_epa) {{
  has_final_state_color_ = false;
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {{
    return {{mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}};
  }}
  const bool epa_hard = coherent_epa &&
                        mg5helas::HasEPAHardBeamState(lts);
  const auto status = subprocess_sum.Prepare(lts, alpha_s, epa_hard);
  if (!mg5helas::EvaluationSucceeded(status)) {{
    return {{status, 0.0}};
  }}
  const auto *channel = subprocess_sum.MatchingChannel(
      lts, {{PDG::PDG_gamma, PDG::PDG_gamma}});
  if (channel == nullptr) {{
    return {{mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}};
  }}
  has_final_state_color_ = mg5::ChannelHasFinalStateColor(*channel);
  double amp2 = 0.0;
  const auto amplitude_status =
      subprocess_sum.CalcPreparedPhotonHelicityAmp2(lts, coherent_epa, amp2);
  if (!mg5helas::EvaluationSucceeded(amplitude_status)) {{
    has_final_state_color_ = false;
  }}
  return {{amplitude_status, amp2}};
}}

// Sample one generated final-state color flow
bool {amplitude}::SampleColorFlow(LORENTZSCALAR &lts, MRandom &random) {{
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {{
    return false;
  }}
  if (!has_final_state_color_) {{
    return ClearHardColorFlow(lts);
  }}
  const auto *channel = subprocess_sum.MatchingChannel(
      lts, {{PDG::PDG_gamma, PDG::PDG_gamma}});
  if (channel == nullptr) {{
    return false;
  }}
  return SampleHardColorFlow(
      lts, random, lts.hard_color_flows,
      channel->external_color_flows.size());
}}

// Compute whether the selected event channel has colored stable final states
bool {amplitude}::HasFinalStateColor() const {{
  return has_final_state_color_;
}}

// Compute the number of generated subprocesses
std::size_t {amplitude}::SubprocessCount() const {{
  return subprocess_sum.SubprocessCount();
}}

}}  // namespace gra
"""
