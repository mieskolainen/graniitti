#!/usr/bin/env python3
#
# Generate the registry for MG5 Durham QCD hard processes
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from pathlib import Path
from typing import Any

from core.io.serialize import load_json_file

from . import output_layout, process_registry
from .process_registry import parse_process_syntax, stable_tokens, topology_signature
from .registry_support import (
    combined_projectors,
    cpp_complex_matrix,
    cpp_flow_candidates,
    cpp_int_vector,
    validate_color_data,
)


# Load and validate one checked-in Durham process data file
def load_process_data(base_dir: Path, entry: dict[str, Any]) -> dict[str, Any]:
    path = base_dir / output_layout.color_path(entry["projection"], entry["name"])
    if not path.is_file():
        raise RuntimeError(f"Missing generated Durham process_data: {path}")
    process_data = load_json_file(path)
    validate_color_data(process_data, entry, [21, 21], 1, "Durham")
    parsed = parse_process_syntax(entry["process"])
    if len(stable_tokens(parsed)) != len(process_data["final_pdgs"]):
        raise RuntimeError(f"Invalid Durham stable final-state size for {entry['name']}")
    return process_data


# Validate every generated Durham topology and stable final state
def validate_signatures(entries: list[tuple[dict[str, Any], dict[str, Any]]]) -> None:
    owners: dict[tuple[str, tuple[int, ...]], str] = {}
    for entry, process_data in entries:
        signature = (
            topology_signature(parse_process_syntax(entry["process"])),
            tuple(int(pdg) for pdg in process_data["final_pdgs"]),
        )
        if signature in owners:
            raise RuntimeError(
                "Duplicate Durham topology and final-state signature "
                f"{signature}: {owners[signature]} and {entry['name']}"
            )
        owners[signature] = entry["name"]


# Generate the Durham process base class
def registry_header() -> str:
    return """// Generated Durham MadGraph process registry
//
// (c) 2017-2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_DURHAMREGISTRY_H
#define AMP_MG5_DURHAMREGISTRY_H

#include <complex>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "@HELICITY_INCLUDE@"
#include "@PROCESS_INCLUDE@"
#include "Graniitti/Kinematics/MKinematics.h"
#include "Graniitti/Particle/MParticle.h"

namespace gra {

// Amplitudes from one generated Durham hard-process evaluation
struct DurhamMG5Evaluation {
  mg5helas::EvaluationStatus status = mg5helas::EvaluationStatus::Success;
  std::vector<std::complex<double>> projected;
  std::vector<std::vector<std::complex<double>>> flow_projected;

  // Compute true when event kinematics and amplitudes were evaluated
  bool Valid() const { return mg5helas::EvaluationSucceeded(status); }
};

// Base class for one generated finite-Nc Durham process
class DurhamMG5Process : public MG5Process {
 public:
  // Store the exact process represented by this matrix element
  explicit DurhamMG5Process(
      std::vector<amplitude::Process> processes)
      : MG5Process(std::move(processes)) {}

  // Destroy one generated matrix element through its base class
  virtual ~DurhamMG5Process() = default;

  // Compute the stable MG2GRA process name
  virtual const std::string &Name() const = 0;
  // Compute final-state PDG codes in MadGraph momentum order
  virtual const std::vector<int> &FinalPDGs() const = 0;
  // Compute final-state SU(3) representations in MadGraph momentum order
  virtual const std::vector<int> &FinalColorRepresentations() const = 0;
  // Expose the event-dependent decay structure of the registered process
  using amplitude::ProcessRegistry::DecayStructureFor;
  // Compute the generated matrix-element decay structure
  DecayStructure DecayStructureFor() const {
    return Processes().at(0).decay_structure;
  }
  // Compute the number of MadGraph color-basis tensors
  virtual std::size_t ColorCount() const = 0;
  // Compute the rank of the incoming-singlet restricted color space
  virtual std::size_t ColorRank() const = 0;
  // Compute the complete MadGraph helicity count
  virtual std::size_t HelicityCount() const = 0;
  // Compute exact finite-Nc projector rows in MadGraph color-basis order
  virtual const std::vector<std::complex<double>> &ExactProjectors() const = 0;
  // Compute leading-color final-state candidates used only for shower tags
  virtual const std::vector<std::vector<MColorFlow>> &FlowCandidates() const = 0;
  // Evaluate exact and shower-partitioned projected helicity amplitudes
  virtual DurhamMG5Evaluation Evaluate(LORENTZSCALAR &lts, double alpha_s,
                                       M4Vec *hard_k1 = nullptr,
                                       M4Vec *hard_k2 = nullptr) = 0;
};

// Assign one generated shower-flow candidate to ordered stable decay leaves
bool AssignDurhamColorFlowCandidate(
    LORENTZSCALAR &lts, const std::vector<MColorFlow> &candidate);

// Compute true when an exact generated Durham matrix element matches this process
bool HasDurhamMG5Process(const amplitude::Process &process);

// Construct the generated Durham matrix element for this exact process
std::unique_ptr<DurhamMG5Process> CreateDurhamMG5Process(
    const amplitude::Process &process);

}  // namespace gra

#endif
""".replace(
        "@HELICITY_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Helicity.h"),
    ).replace(
        "@PROCESS_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Process.h"),
    )


# Generate one matrix element construction block for a Durham process
def matrix_element_block(entry: dict[str, Any], process_data: dict[str, Any]) -> str:
    name = entry["name"]
    class_name = f"AMP_MG5_{name}"
    signature = cpp_int_vector(process_data["final_pdgs"])
    if entry.get("type", "tree") == "loop":
        return f"""  if (ProcessMatches(process, "{name}", {signature})) {{
    return std::make_unique<{class_name}>();
  }}
"""
    exact = process_data["projectors"]
    combined = combined_projectors(process_data)
    flows = cpp_flow_candidates(process_data["flow_candidates"])
    return f"""  if (ProcessMatches(process, "{name}", {signature})) {{
    return std::make_unique<GeneratedDurhamProcess<{class_name}>>(
        "{name}", std::vector<int>{signature},
        std::vector<int>{cpp_int_vector(process_data["final_color_representations"])},
        {int(process_data["rank"])},
        std::vector<std::complex<double>>{cpp_complex_matrix(exact)},
        std::vector<std::complex<double>>{cpp_complex_matrix(combined)},
        std::vector<std::vector<MColorFlow>>{flows});
  }}
"""


# Generate one process match block without constructing a matrix element
def process_match_block(entry: dict[str, Any], process_data: dict[str, Any]) -> str:
    signature = cpp_int_vector(process_data["final_pdgs"])
    return f'  if (ProcessMatches(process, "{entry["name"]}", {signature})) {{ return true; }}\n'


# Generate the concrete registry implementation for all Durham entries
def registry_source(entries: list[tuple[dict[str, Any], dict[str, Any]]]) -> str:
    includes = "\n".join(
        f'#include "{output_layout.include_path(output_layout.DURHAM, f"AMP_MG5_{entry['name']}.h")}"'
        for entry, _ in entries
    )
    matches = "".join(process_match_block(entry, process_data) for entry, process_data in entries)
    matrix_elements = "".join(
        matrix_element_block(entry, process_data) for entry, process_data in entries
    )
    return f"""// Generated Durham MadGraph process registry
//
// [REFERENCE: Alwall et al., JHEP 07 (2014) 079, arXiv:1405.0301]
// [REFERENCE: Sjodahl, JHEP 09 (2009) 087, arXiv:0906.1121]

#include "{output_layout.include_path(output_layout.DURHAM, "AMP_MG5_DurhamRegistry.h")}"
#include "{output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h")}"

#include <algorithm>
#include <initializer_list>
#include <stdexcept>
#include <utility>

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
      amplitude::FindProcess("DURHAM", process_name);
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

// Represent one concrete generated matrix element in one Durham process
template <class MatrixElement>
class GeneratedDurhamProcess final : public DurhamMG5Process {{
 public:
  // Store immutable exact projectors and their leading-color shower partition
  GeneratedDurhamProcess(
      std::string name, std::vector<int> final_pdgs,
      std::vector<int> final_color_representations, std::size_t rank,
      std::vector<std::complex<double>> exact_projectors,
      std::vector<std::complex<double>> combined_projectors,
      std::vector<std::vector<MColorFlow>> flow_candidates)
      : DurhamMG5Process(
            amplitude::Processes("DURHAM", name)),
        name_(std::move(name)),
        final_pdgs_(std::move(final_pdgs)),
        final_color_representations_(std::move(final_color_representations)),
        rank_(rank),
        exact_projectors_(std::move(exact_projectors)),
        combined_projectors_(std::move(combined_projectors)),
        flow_candidates_(std::move(flow_candidates)) {{
    const std::size_t exact_size = rank_ * MatrixElement::ncolor;
    const std::size_t combined_rows =
        rank_ * (1 + flow_candidates_.size());
    if (rank_ == 0 ||
        final_color_representations_.size() != final_pdgs_.size() ||
        exact_projectors_.size() != exact_size ||
        combined_projectors_.size() != combined_rows * MatrixElement::ncolor ||
        MatrixElement::nhelicity % 4 != 0 ||
        !gra::AllFinite(exact_projectors_) || !gra::AllFinite(combined_projectors_)) {{
      throw std::invalid_argument(
          "GeneratedDurhamProcess: inconsistent generated data for " + name_);
    }}
    combined_amplitudes_.resize(combined_rows * MatrixElement::nhelicity);
  }}

  // Compute the stable MG2GRA process name
  const std::string &Name() const override {{ return name_; }}

  // Compute final-state PDG codes in MadGraph momentum order
  const std::vector<int> &FinalPDGs() const override {{ return final_pdgs_; }}

  // Compute final-state SU(3) representations in momentum order
  const std::vector<int> &FinalColorRepresentations() const override {{
    return final_color_representations_;
  }}

  // Compute the number of MadGraph color-basis tensors
  std::size_t ColorCount() const override {{ return MatrixElement::ncolor; }}

  // Compute the finite-Nc restricted color rank
  std::size_t ColorRank() const override {{ return rank_; }}

  // Compute the complete MadGraph helicity count
  std::size_t HelicityCount() const override {{ return MatrixElement::nhelicity; }}

  // Compute the one generated matrix element owned by this exact process
  std::size_t SubprocessCount() const override {{ return 1; }}

  // Initialize all model parameters before sampling
  void InitParameters(SLHAReader card) override {{ matrix_element_.InitParameters(std::move(card)); }}

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override {{ return matrix_element_.Particles(); }}

  // Compute exact finite-Nc projector rows
  const std::vector<std::complex<double>> &ExactProjectors() const override {{
    return exact_projectors_;
  }}

  // Compute shower-compatible leading-color final-state candidates
  const std::vector<std::vector<MColorFlow>> &FlowCandidates() const override {{
    return flow_candidates_;
  }}

  // Evaluate all exact and shower-partitioned projections in one HELAS pass
  DurhamMG5Evaluation Evaluate(LORENTZSCALAR &lts, double alpha_s,
                               M4Vec *hard_k1, M4Vec *hard_k2) override {{
    DurhamMG5Evaluation result;
    const std::size_t combined_rows =
        rank_ * (1 + flow_candidates_.size());
    std::fill(combined_amplitudes_.begin(), combined_amplitudes_.end(), 0.0);
    result.status = matrix_element_.CalcColorProjectedHelicity(
        lts, alpha_s, combined_projectors_.data(),
        static_cast<int>(combined_rows), combined_amplitudes_.data(), hard_k1,
        hard_k2);
    if (!result.Valid()) {{
      return result;
    }}

    const std::size_t exact_size = rank_ * MatrixElement::nhelicity;
    result.projected.assign(
        combined_amplitudes_.begin(), combined_amplitudes_.begin() + exact_size);
    result.flow_projected.resize(flow_candidates_.size());
    for (std::size_t flow = 0; flow < flow_candidates_.size(); ++flow) {{
      const auto first =
          combined_amplitudes_.begin() + exact_size * (flow + 1);
      result.flow_projected[flow].resize(exact_size);
      std::copy(first, first + exact_size, result.flow_projected[flow].begin());
    }}
    return result;
  }}

 private:
  MatrixElement matrix_element_;
  std::string name_;
  std::vector<int> final_pdgs_;
  std::vector<int> final_color_representations_;
  std::size_t rank_;
  std::vector<std::complex<double>> exact_projectors_;
  std::vector<std::complex<double>> combined_projectors_;
  std::vector<std::vector<MColorFlow>> flow_candidates_;
  std::vector<std::complex<double>> combined_amplitudes_;
}};

}}  // namespace

// Assign one shower-flow candidate to stable final-state leaves
bool AssignDurhamColorFlowCandidate(
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

// Compute true when an exact generated matrix element matches this process
bool HasDurhamMG5Process(const amplitude::Process &process) {{
{matches}  return false;
}}

// Construct the generated matrix element for this exact process
std::unique_ptr<DurhamMG5Process> CreateDurhamMG5Process(
    const amplitude::Process &process) {{
{matrix_elements}  return nullptr;
}}

}}  // namespace gra
"""


# Load all Durham process data and generate the registry
def generate_registry(base_dir: Path, manifest: dict[str, Any]) -> tuple[str, str]:
    entries = [
        (entry, load_process_data(base_dir, entry))
        for entry in manifest["processes"]
        if entry["projection"] == process_registry.DURHAM_PROJECTION
    ]
    validate_signatures(entries)
    return registry_header(), registry_source(entries)
