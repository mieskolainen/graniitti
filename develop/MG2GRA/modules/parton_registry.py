#!/usr/bin/env python3
#
# Generate the MadGraph parton runtime registry
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

from typing import Any

from . import output_layout, process_registry


# Compute generated families which implement the common parton event
def parton_families(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    families: list[dict[str, Any]] = []
    for family in manifest.get("families", []):
        if family["projection"] != process_registry.PARTON_PROJECTION:
            continue
        process_registry.family_process_family(family)
        families.append(family)
    return families


# Generate one complete parton-family amplitude declaration
def amplitude_header(family: dict[str, Any]) -> str:
    family_name = family["name"]
    amplitude = process_registry.family_wrapper_name(family)
    channel = family["channel"]
    guard = f"{amplitude.upper()}_H"
    family_directory = output_layout.family_cpp_directory(
        family["projection"], family_name
    )
    registry_include = output_layout.include_path(
        output_layout.PARTON, "AMP_MG5_PartonRegistry.h"
    )
    subprocess_include = output_layout.include_path(
        output_layout.RUNTIME, "AMP_MG5_SubprocessSum.h"
    )
    process_base_include = output_layout.include_path(
        family_directory, "ProcessBase.h"
    )
    return f"""// Generated MG5 amplitude for parton {channel} production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef {guard}
#define {guard}

#include <cstddef>
#include <vector>

#include "{registry_include}"
#include "{subprocess_include}"
#include "{process_base_include}"
#include "Graniitti/Kinematics/MKinematics.h"

namespace gra {{

// Sum every generated subprocess for parton {channel} production
class {amplitude} : public PartonMG5Process {{
 public:
  // Construct all generated subprocesses
  {amplitude}();

  // Destroy all generated subprocesses
  ~{amplitude}() override = default;

  {amplitude}(const {amplitude} &) = delete;
  {amplitude} &operator=(const {amplitude} &) = delete;

  // Prepare event-local momenta and final-state channels
  mg5helas::EvaluationStatus Prepare(LORENTZSCALAR &lts,
                                     double alpha_s) override;

  // Evaluate one prepared event and return all event-local results
  PartonMG5Evaluation EvaluatePrepared(LORENTZSCALAR &lts,
                                     double alpha_s) override;

  // Compute the number of generated subprocesses
  std::size_t SubprocessCount() const override;

  // Initialize all model parameters before sampling
  void InitParameters(SLHAReader card) override {{ subprocess_sum.InitParameters(std::move(card)); }}

  // Compute the evaluated pole parameters of the generated model
  mg5::ParticleMap Particles() const override {{ return subprocess_sum.Particles(); }}

 private:
  mg5::SubprocessSum<{family_name}::ProcessBase> subprocess_sum;
}};

}}  // namespace gra

#endif
"""


# Generate one complete parton-family amplitude implementation
def amplitude_source(family: dict[str, Any]) -> str:
    family_name = family["name"]
    amplitude = process_registry.family_wrapper_name(family)
    channel = family["channel"]
    parameter_card = output_layout.source_path(
        output_layout.parameter_card_path(family["projection"], family_name)
    )
    family_directory = output_layout.family_cpp_directory(
        family["projection"], family_name
    )
    amplitude_include = output_layout.include_path(
        output_layout.PARTON, f"{amplitude}.h"
    )
    processes_include = output_layout.include_path(
        family_directory, "Processes.h"
    )
    return f"""// Generated MG5 amplitude for parton {channel} production
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include <utility>
#include <vector>

#include "{amplitude_include}"
#include "{processes_include}"
#include "Graniitti/Tech/MAux.h"

namespace gra {{

// Construct all generated subprocesses
{amplitude}::{amplitude}()
    : PartonMG5Process(amplitude::Processes("{family_name}")) {{
  std::vector<mg5::Subprocess<{family_name}::ProcessBase>> subprocesses;
  {family_name}::BuildSubprocesses(
      subprocesses,
      aux::ResolveProjectPath("{parameter_card}"));
  mg5::ValidateSubprocessInitialStates(
      subprocesses, mg5::MasslessQCDInitialStates(), "{amplitude}");
  subprocess_sum.SetSubprocesses(std::move(subprocesses));
}}

// Prepare event-local momenta and final-state channels
mg5helas::EvaluationStatus {amplitude}::Prepare(LORENTZSCALAR &lts,
                                                double alpha_s) {{
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {{
    return mg5helas::EvaluationStatus::AmplitudeFailure;
  }}
  return subprocess_sum.Prepare(lts, alpha_s);
}}

// Evaluate one prepared event and return all event-local results
PartonMG5Evaluation {amplitude}::EvaluatePrepared(LORENTZSCALAR &lts,
                                                double alpha_s) {{
  PartonMG5Evaluation result;
  if (!subprocess_sum.RequireGeneratedTopology(*this, lts)) {{
    result.status = mg5helas::EvaluationStatus::AmplitudeFailure;
    return result;
  }}
  if (!subprocess_sum.PreparedStateMatches(lts, alpha_s)) {{
    result.status = Prepare(lts, alpha_s);
  }}
  if (!result.Valid()) {{
    return result;
  }}
  result.status = subprocess_sum.CalcPreparedHelicityAmp2(lts, result.amp2);
  result.contributing_subprocesses =
      subprocess_sum.ContributingSubprocesses();
  result.color_flows = subprocess_sum.LastHardColorFlows();
  return result;
}}

// Compute the number of generated subprocesses
std::size_t {amplitude}::SubprocessCount() const {{
  return subprocess_sum.SubprocessCount();
}}

}}  // namespace gra
"""


# Generate the parton process base class
def registry_header() -> str:
    return """// Generated parton MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PARTONREGISTRY_H
#define AMP_MG5_PARTONREGISTRY_H

#include <complex>
#include <cstddef>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "@HELICITY_INCLUDE@"
#include "@PROCESS_INCLUDE@"
#include "@PROCESS_REGISTRY_INCLUDE@"
#include "Graniitti/Kinematics/MKinematics.h"

namespace gra {

// Own all values produced by one generated parton-process evaluation
struct PartonMG5Evaluation {
  mg5helas::EvaluationStatus status = mg5helas::EvaluationStatus::Success;
  double                     amp2   = 0.0;
  std::size_t                contributing_subprocesses = 0;
  std::vector<mg5helas::HardColorFlow>           color_flows;

  // Compute true when event kinematics and amplitudes were evaluated
  bool Valid() const { return mg5helas::EvaluationSucceeded(status); }
};

// Base class for one generated parton subprocess family
class PartonMG5Process : public MG5Process {
 public:
  // Store the exact process family represented by this matrix element
  explicit PartonMG5Process(std::vector<amplitude::Process> processes)
      : MG5Process(std::move(processes)) {}

  // Destroy one generated parton process family through its base class
  virtual ~PartonMG5Process() = default;

  // Prepare event momenta and final-state channels
  virtual mg5helas::EvaluationStatus Prepare(LORENTZSCALAR &lts, double alpha_s) = 0;

  // Evaluate one prepared event and return all event-local results by value
  virtual PartonMG5Evaluation EvaluatePrepared(LORENTZSCALAR &lts, double alpha_s) = 0;

  // Compute the number of generated subprocess channels
  virtual std::size_t SubprocessCount() const = 0;
};

// Compute generated parton families and their public process channels
std::vector<MG5ProcessInfo> PartonMG5ProcessInfos();

// Construct one generated parton process family
std::unique_ptr<PartonMG5Process> CreatePartonMG5Process(
    const std::string &process_family);

}  // namespace gra

#endif
""".replace(
        "@HELICITY_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Helicity.h"),
    ).replace(
        "@PROCESS_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Process.h"),
    ).replace(
        "@PROCESS_REGISTRY_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h"),
    )


# Generate one manifest driven parton process factory implementation
def registry_source(families: list[dict[str, Any]]) -> str:
    includes = "\n".join(
        f'#include "{output_layout.include_path(output_layout.PARTON, f"{process_registry.family_wrapper_name(family)}.h")}"'
        for family in families
    )
    if includes:
        includes = f"\n{includes}\n"
    process_infos = ",\n".join(
        f'      MG5ProcessInfo{{"{family["name"]}", "{family["channel"]}", '
        f'"{output_layout.source_path(output_layout.parameter_card_path(family["projection"], family["name"]))}"}}'
        for family in families
    )
    matrix_elements = "".join(
        f'  if (process_family == "{family["name"]}") {{\n'
        f"    return std::make_unique<{process_registry.family_wrapper_name(family)}>();\n"
        "  }\n"
        for family in families
    )
    return f"""// Generated parton MadGraph process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "{output_layout.include_path(output_layout.PARTON, "AMP_MG5_PartonRegistry.h")}"

#include <memory>
#include <string>
#include <vector>
{includes}
namespace gra {{

// Compute generated parton families and their public process channels
std::vector<MG5ProcessInfo> PartonMG5ProcessInfos() {{
  return {{
{process_infos}
  }};
}}

// Construct one generated parton process family
std::unique_ptr<PartonMG5Process> CreatePartonMG5Process(
    const std::string &process_family) {{
{matrix_elements}  return nullptr;
}}

}}  // namespace gra
"""


# Generate the parton process base class and process family construction
def generate_registry(manifest: dict[str, Any]) -> tuple[str, str]:
    families = parton_families(manifest)
    return registry_header(), registry_source(families)
