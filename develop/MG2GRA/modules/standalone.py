#!/usr/bin/env python3
#
# Convert MG5 standalone C++ amplitudes for GRANIITTI
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re

from . import output_layout
from .cpp_support import narrow_std_namespace as qualify_std_namespace
from .cpp_support import normalize_whitespace

MESSAGE_ORIG = "// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch"
MESSAGE_NEW = (
    "// Visit launchpad.net/madgraph5 and amcatnlo.web.cern.ch\n"
    "// @@@@ MadGraph to GRANIITTI conversion done @@@@"
)


def fail(message: str) -> None:
    raise RuntimeError(message)


def replace_or_fail(text: str, old: str, new: str, count: int = -1) -> str:
    occurrences = text.count(old)
    if occurrences == 0:
        fail(f"Could not find replacement anchor: {old}")

    if count != -1 and occurrences != count:
        fail(f"Expected {count} replacement anchor(s), found {occurrences}: {old}")

    if count == -1:
        return text.replace(old, new)
    return text.replace(old, new, count)


def sub_or_fail(pattern: str, repl: str, text: str, *, count: int = 0, flags: int = 0) -> str:
    updated, nmatch = re.subn(pattern, repl, text, count=count, flags=flags)
    if nmatch == 0:
        fail(f"Could not match pattern: {pattern}")
    return updated


def extract_or_fail(pattern: str, text: str, *, flags: int = 0) -> re.Match[str]:
    match = re.search(pattern, text, flags)
    if match is None:
        fail(f"Could not extract pattern: {pattern}")
    return match


def indent_block(text: str, spaces: int) -> str:
    prefix = " " * spaces
    return "\n".join(prefix + line if line.strip() else line for line in text.splitlines())


# Replace broad generated std imports with explicit qualifications
def narrow_std_namespace(text: str) -> str:
    return normalize_whitespace(qualify_std_namespace(text))


def extract_color_structure(source: str) -> dict[str, str | int] | None:
    ncolor_match = re.search(r"const int ncolor = (\d+);", source)
    if ncolor_match is None:
        return None

    nhelicity = int(extract_or_fail(r"const int ncomb = (\d+);", source).group(1))
    matrix_name = extract_or_fail(r"double CPPProcess::(matrix_[^(]+)\(\)", source).group(1)
    ngraphs = extract_or_fail(r"const int ngraphs = (\d+);", source).group(1)

    denom = extract_or_fail(
        r"static const double denom\[ncolor\] = (\{.*?\});", source, flags=re.S
    ).group(1)
    cf = extract_or_fail(
        r"static const double cf\[ncolor\]\[ncolor\] = (\{\{.*?\}\});", source, flags=re.S
    ).group(1)
    color_flows = extract_or_fail(
        r"// Calculate color flows\n(.*?)\n\s*// Sum and square the color flows", source, flags=re.S
    ).group(1)
    color_flows = color_flows.strip("\n")

    return {
        "matrix_name": matrix_name,
        "ncolor": int(ncolor_match.group(1)),
        "nhelicity": nhelicity,
        "ngraphs": int(ngraphs),
        "denom": denom,
        "cf": cf,
        "color_flows": color_flows,
    }


def transform_header(
    process: str,
    raw_header: str,
    color_meta: dict[str, str | int] | None,
    projection: str,
    model: str,
) -> str:
    class_name = f"AMP_MG5_{process}"
    photon_process = projection == "photon"
    parameter_card = output_layout.source_path(
        output_layout.parameter_card_path(projection, process)
    )
    header = normalize_whitespace(raw_header)
    parameter_match = extract_or_fail(r'#include "(Parameters_[^"]+)\.h"', raw_header)
    parameter_class = parameter_match.group(1)
    parameter_header = f"{parameter_class}.h"
    model_directory = output_layout.model_cpp_directory(model)
    runtime_directory = output_layout.RUNTIME
    header = replace_or_fail(header, MESSAGE_ORIG, MESSAGE_NEW, 1)
    header = replace_or_fail(header, "CPPProcess", class_name)
    header = replace_or_fail(header, "virtual ", "")
    support_includes = (
        f'#include "{output_layout.include_path(model_directory, parameter_header)}"\n'
        f'#include "{output_layout.include_path(runtime_directory, "AMP_MG5_Helicity.h")}"\n'
        f'#include "{output_layout.include_path(runtime_directory, "AMP_MG5_ProcessRegistry.h")}"\n'
        f'#include "{output_layout.include_path(runtime_directory, "AMP_MG5_Kinematics.h")}"\n'
    )
    replacement_includes = (
        "#include <algorithm>\n"
        "#include <array>\n"
        "#include <cmath>\n"
        "#include <complex>\n"
        "#include <stdexcept>\n"
        "#include <vector>\n\n" + support_includes + '#include "Graniitti/Tech/MAux.h"\n'
        '#include "Graniitti/Particle/MForm.h"\n'
        '#include "Graniitti/Kinematics/MKinematics.h"\n'
        '#include "Graniitti/Math/MMath.h"'
    )

    header = sub_or_fail(
        rf'#include <complex>\s*\n#include <vector>\s*\n\s*#include "{re.escape(parameter_header)}"',
        replacement_includes,
        header,
        flags=re.S,
    )

    header = sub_or_fail(
        rf"{class_name}\(\)\s*\{{\s*\}}",
        (
            f'{class_name}() {{ initProc(gra::aux::ResolveProjectPath("{parameter_card}")); }}\n'
            "  // Keep event-local MG5 work arrays isolated between process instances\n"
            f"  {class_name}(const {class_name} &) = delete;\n"
            f"  {class_name} &operator=(const {class_name} &) = delete;"
        ),
        header,
        count=1,
    )
    pointer_declarations = (f"{parameter_class} * pars;", f"{parameter_class} *pars;")
    pointer_matches = sum(header.count(declaration) for declaration in pointer_declarations)
    if pointer_matches != 1:
        fail(f"Expected one generated parameter pointer, found {pointer_matches}")
    for declaration in pointer_declarations:
        header = header.replace(declaration, f"{parameter_class} pars;  // GRANIITTI")
    calc_signature = (
        "gra::mg5helas::MatrixElementEvaluation Evaluate(gra::LORENTZSCALAR &lts, double alphas"
    )
    calc_signature += ", bool coherent_epa" if photon_process else ""
    calc_signature += ");"
    header = replace_or_fail(header, "void sigmaKin();", calc_signature, 1)

    if color_meta is not None:
        header = sub_or_fail(
            rf"(class {class_name}\s*\{{\s*public:\s*)",
            (
                r"\1\n"
                "  using ColorFlowVector          = std::vector<std::complex<double>>;\n"
                "  using ColorFlowHelicityMatrix  = std::vector<ColorFlowVector>;\n"
            ),
            header,
            count=1,
            flags=re.S,
        )
        header = replace_or_fail(
            header,
            calc_signature,
            (
                calc_signature + "\n"
                "  void   CalcColorFlowHelicity(gra::LORENTZSCALAR &lts, double alphas,\n"
                "                               ColorFlowHelicityMatrix &jamp_matrix);\n"
                "  // Contract arbitrary color projectors while the generated amplitudes are live\n"
                "  gra::mg5helas::EvaluationStatus CalcColorProjectedHelicity(\n"
                "      gra::LORENTZSCALAR &lts, double alphas,\n"
                "      const std::complex<double> *color_projectors, int projector_count,\n"
                "      std::complex<double> *projected, gra::M4Vec *hard_k1 = nullptr,\n"
                "      gra::M4Vec *hard_k2 = nullptr);"
                + (
                    "\n  // Prepare the transfer-independent EPA hard tensor and frame\n"
                    "  gra::mg5helas::EvaluationStatus CalcColorProjectedEPAHardHelicity(\n"
                    "      gra::LORENTZSCALAR &lts, double alphas,\n"
                    "      const std::complex<double> *color_projectors, int projector_count,\n"
                    "      std::complex<double> *projected,\n"
                    "      gra::mg5helas::EPAHardFrame *frame);"
                    if photon_process
                    else ""
                )
            ),
            1,
        )

    header = sub_or_fail(
        r"  // Constants for array limits\s+"
        r"static const int ninitial\s*=\s*2;\s+"
        r"static const int nexternal\s*=\s*",
        "  // Constants for array limits\n"
        "  static const int ninitial   = 2;\n"
        "  static const int nexternal  = ",
        header,
        count=1,
    )
    header = sub_or_fail(
        r"\s*static const int nprocesses\s*=\s*1;",
        "\n  static const int nprocesses = 1;",
        header,
        count=1,
    )

    if color_meta is not None:
        nexternal = extract_or_fail(r"static const int nexternal = (\d+);", raw_header).group(1)
        header = sub_or_fail(
            rf"  static const int nexternal\s*=\s*{nexternal};\s*"
            r"static const int nprocesses\s*=\s*1;",
            "  static const int nexternal  = "
            + nexternal
            + ";\n"
            + f"  static const int ncolor     = {color_meta['ncolor']};\n"
            + f"  static constexpr int nhelicity = {color_meta['nhelicity']};\n"
            + "  static const int nprocesses = 1;\n\n"
            + "  std::vector<double>              ColorDenominators() const;\n"
            + "  std::vector<std::vector<double>> ColorMetric() const;",
            header,
            count=1,
        )

    header = sub_or_fail(
        r"\s*private:\s+// Private functions to calculate the matrix element for all subprocesses\s+",
        (
            "\n private:\n"
            "  // Private functions to calculate the matrix element for all subprocesses\n"
            "  // Prepare one physical on-shell HELAS phase-space point\n"
            "  bool setup_kinematics(gra::LORENTZSCALAR &lts"
            + (
                ", gra::mg5helas::EPAHardFrame *frame = nullptr);\n"
                if photon_process
                else ");\n"
            )
            + (
                "  // Contract prepared HELAS amplitudes into arbitrary color projectors\n"
                "  gra::mg5helas::EvaluationStatus CalcColorProjectedPrepared(\n"
                "      const std::complex<double> *color_projectors, int projector_count,\n"
                "      std::complex<double> *projected);\n"
                if photon_process and color_meta is not None
                else ""
            )
            + (
                "  void calculate_color_flows(std::complex<double> jamp[ncolor]) const;\n"
                if color_meta is not None
                else ""
            )
        ),
        header,
        count=1,
    )

    header = sub_or_fail(
        rf"class {class_name}\s*\{{",
        (f"class {class_name}\n    : public gra::amplitude::MG5ProcessRegistry_{process} {{"),
        header,
        count=1,
    )

    header = re.sub(
        r"void initProc\((?:std::)?string param_card_name\);",
        "void initProc(std::string param_card_name);\n"
        "    // Initialize every mass and coupling through the generated model\n"
        "    void InitParameters(SLHAReader slha);\n"
        "    // Compute the evaluated particle parameters of the generated model\n"
        "    gra::mg5::ParticleMap Particles() const { return pars.Particles(); }\n"
        "    // Compute the evaluated UFO electromagnetic coupling\n"
        "    double AlphaQED() const { return pars.AlphaQED(); }",
        header,
        count=1,
    )
    header = header.replace("    void calculate_wavefunctions", "  void calculate_wavefunctions")
    header = header.replace(
        "    static const int nwavefuncs = ", "  static const int nwavefuncs = "
    )
    header = header.replace("    std::complex<double> w", "  std::complex<double> w")
    header = header.replace(
        "    static const int namplitudes = ", "  static const int namplitudes = "
    )
    header = header.replace("    std::complex<double> amp", "  std::complex<double> amp")
    header = header.replace("    double matrix_", "  double matrix_")
    header = header.replace("    double matrix_element", "  double matrix_element")
    header = header.replace("    double * jamp2", "  double *jamp2")
    header = sub_or_fail(
        r"  double \*jamp2\[nprocesses\];\s*",
        "  std::array<std::array<double, ncolor>, nprocesses> jamp2 = {};\n",
        header,
        count=1,
    )
    header = header.replace(
        f"    {parameter_class} pars;  // GRANIITTI",
        f"  {parameter_class} pars;  // GRANIITTI",
    )
    header = header.replace("    vector<double> mME;", "  vector<double> mME;")
    header = sub_or_fail(
        r"\s*vector\s*<\s*double\s*\*\s*>\s*p;\s*// Initial particle ids",
        "\n  vector<double *> p;\n"
        "  vector<array<double, 4>> momentum_buffer;\n"
        "  vector<gra::M4Vec> final_buffer;\n"
        "  // Initial particle ids",
        header,
        count=1,
    )
    header = sub_or_fail(r"\s*int id1,\s*id2;", "\n  int id1, id2;", header, count=1)

    # Compact indentation from the original MadGraph formatting
    header = header.replace("\n{\n", "\n{\n")
    header = header.replace("\n  public:\n", "\n public:\n")

    return narrow_std_namespace(header)


# Generate the physical HELAS momenta for an initialized decay topology
def make_setup_kinematics(class_name: str, decay_chain: bool, photon_process: bool) -> str:
    final_state_setup = (
        """  const auto leaves = gra::mg5::StableDecayLeaves(lts.decaytree);
  if (leaves.size() != nexternal - ninitial) { return false; }
  final_buffer.reserve(leaves.size());
  for (const auto *leaf : leaves) {
    final_buffer.push_back(leaf->p4);
  }"""
        if decay_chain
        else """  final_buffer.reserve(lts.decaytree.size());
  for (std::size_t i = 0; i < lts.decaytree.size(); ++i) {
    final_buffer.push_back(lts.decaytree[i].p4);
  }"""
    )
    signature = (
        f"bool {class_name}::setup_kinematics(\n"
        "    gra::LORENTZSCALAR &lts, gra::mg5helas::EPAHardFrame *frame)"
        if photon_process
        else f"bool {class_name}::setup_kinematics(gra::LORENTZSCALAR &lts)"
    )
    incoming_setup = (
        """  gra::M4Vec p1_;
  gra::M4Vec p2_;
  if (frame != nullptr) {
    if (!gra::mg5helas::PrepareEPAHardFrame(lts, final_buffer, *frame)) { return false; }
    p1_ = frame->incoming[0];
    p2_ = frame->incoming[1];
  } else if (!gra::mg5helas::PrepareOnShellKinematics(
                 lts, final_buffer, p1_, p2_)) {
    return false;
  }"""
        if photon_process
        else """  gra::M4Vec p1_;
  gra::M4Vec p2_;
  if (!gra::mg5helas::PrepareOnShellKinematics(
          lts, final_buffer, p1_, p2_)) { return false; }"""
    )
    return f"""
// Prepare one physical on-shell HELAS phase-space point
{signature} {{
  // *** MADGRAPH CONVENTION IS [E,px,py,pz] ! ***

  // Reuse event buffers because this function runs at every screening node
  final_buffer.clear();

{final_state_setup}
  if (!gra::mg5::OnShellFinal(final_buffer, mME)) {{ return false; }}

{incoming_setup}

  p.clear();
  momentum_buffer.clear();
  momentum_buffer.reserve(ninitial + final_buffer.size());

  momentum_buffer.push_back({{p1_.E(), p1_.Px(), p1_.Py(), p1_.Pz()}});
  momentum_buffer.push_back({{p2_.E(), p2_.Px(), p2_.Py(), p2_.Pz()}});

  for (const auto &p4 : final_buffer) {{
    momentum_buffer.push_back({{p4.E(), p4.Px(), p4.Py(), p4.Pz()}});
  }}

  for (auto &mom : momentum_buffer) {{ p.push_back(mom.data()); }}
  return true;
}}
"""


# Generate complex amplitudes in the MG5 color basis
def make_calc_color_flow_helicity(
    class_name: str, ncomb: int, helicity_block: str, ncolor: int, incoming_pdg: int,
    alpha_zero: bool = True,
) -> str:
    alpha_update = (
        "  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }\n"
        if incoming_pdg == 22 and alpha_zero
        else ""
    )
    return f"""
// Evaluate every generated color flow in the fixed helicity basis
void {class_name}::CalcColorFlowHelicity(gra::LORENTZSCALAR &lts, double alphas,
                                         ColorFlowHelicityMatrix &jamp_matrix) {{
  jamp_matrix.clear();
  if (!std::isfinite(alphas) || alphas < 0.0) {{ return; }}
  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
{alpha_update}

  const int ncomb = {ncomb};
  jamp_matrix.assign(ncomb, ColorFlowVector(ncolor, 0.0));
  if (!setup_kinematics(lts)) {{ return; }}

  {helicity_block.replace("static const int helicities", "const int helicities")}

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) {{ perm[i] = i; }}

  for (int ihel = 0; ihel < ncomb; ++ihel) {{
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);

    for (int icolor = 0; icolor < ncolor; ++icolor) {{
      jamp_matrix[ihel][icolor] = jamp[icolor];
    }}
  }}
}}
"""


# Generate complex amplitudes contracted with fixed color projectors
def make_calc_color_projected_helicity(
    class_name: str, ncomb: int, helicity_block: str, incoming_pdg: int,
    alpha_zero: bool = True,
) -> str:
    alpha_update = (
        "  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }\n"
        if incoming_pdg == 22 and alpha_zero
        else ""
    )
    return f"""
// Contract generated color flows directly into arbitrary hard-process projectors
gra::mg5helas::EvaluationStatus {class_name}::CalcColorProjectedHelicity(
    gra::LORENTZSCALAR &lts, double alphas,
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected, gra::M4Vec *hard_k1,
    gra::M4Vec *hard_k2) {{
  const std::size_t output_size =
      projector_count > 0
          ? static_cast<std::size_t>(projector_count) * nhelicity
          : 0;
  if (projected != nullptr && output_size > 0) {{
    std::fill(projected, projected + output_size,
              std::complex<double>(0.0));
  }}
  if (hard_k1 != nullptr) {{ *hard_k1 = gra::M4Vec(); }}
  if (hard_k2 != nullptr) {{ *hard_k2 = gra::M4Vec(); }}
  if (color_projectors == nullptr || projected == nullptr ||
      projector_count <= 0) {{
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }}
  if (!std::isfinite(alphas) || alphas < 0.0) {{
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }}

  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
{alpha_update}
  if (!setup_kinematics(lts)) {{
    return gra::mg5helas::EvaluationStatus::KinematicsFailure;
  }}

  const int ncomb = nhelicity;
  {helicity_block.replace("static const int helicities", "const int helicities")}

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) {{ perm[i] = i; }}

  for (int ihel = 0; ihel < nhelicity; ++ihel) {{
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    const std::span<const std::complex<double>> jamp_view(jamp, ncolor);

    for (int projector = 0; projector < projector_count; ++projector) {{
      const std::span<const std::complex<double>> projector_view(
          color_projectors + projector * ncolor, ncolor);
      const std::complex<double> value =
          gra::BilinearProduct(projector_view, jamp_view);
      projected[static_cast<std::size_t>(projector) * nhelicity + ihel] = value;
    }}
  }}
  if (hard_k1 != nullptr) {{
    *hard_k1 = gra::M4Vec(p[0][1], p[0][2], p[0][3], p[0][0]);
  }}
  if (hard_k2 != nullptr) {{
    *hard_k2 = gra::M4Vec(p[1][1], p[1][2], p[1][3], p[1][0]);
  }}
  return gra::mg5helas::EvaluationStatus::Success;
}}
"""


# Generate the common projector contraction for prepared photon kinematics
def make_calc_color_projected_prepared(class_name: str, helicity_block: str) -> str:
    return f"""
// Contract prepared HELAS amplitudes into arbitrary color projectors
gra::mg5helas::EvaluationStatus {class_name}::CalcColorProjectedPrepared(
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected) {{
  const int ncomb = nhelicity;
  {helicity_block.replace("static const int helicities", "const int helicities")}

  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) {{ perm[i] = i; }}

  for (int ihel = 0; ihel < nhelicity; ++ihel) {{
    calculate_wavefunctions(perm, helicities[ihel]);

    std::complex<double> jamp[ncolor];
    calculate_color_flows(jamp);
    const std::span<const std::complex<double>> jamp_view(jamp, ncolor);

    for (int projector = 0; projector < projector_count; ++projector) {{
      const std::span<const std::complex<double>> projector_view(
          color_projectors + projector * ncolor, ncolor);
      const std::complex<double> value = gra::BilinearProduct(projector_view, jamp_view);
      projected[static_cast<std::size_t>(projector) * nhelicity + ihel] = value;
    }}
  }}
  return gra::mg5helas::EvaluationStatus::Success;
}}
"""


# Generate the transfer-independent EPA hard-tensor evaluator
def make_calc_color_projected_epa_hard(class_name: str, alpha_zero: bool = True) -> str:
    alpha_update = "  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }" if alpha_zero else ""
    return f"""
// Prepare the transfer-independent EPA hard tensor and its reusable frame
gra::mg5helas::EvaluationStatus
{class_name}::CalcColorProjectedEPAHardHelicity(
    gra::LORENTZSCALAR &lts, double alphas,
    const std::complex<double> *color_projectors, int projector_count,
    std::complex<double> *projected, gra::mg5helas::EPAHardFrame *frame) {{
  const std::size_t output_size =
      projector_count > 0
          ? static_cast<std::size_t>(projector_count) * nhelicity
          : 0;
  if (projected != nullptr && output_size > 0) {{
    std::fill(projected, projected + output_size,
              std::complex<double>(0.0));
  }}
  if (frame != nullptr) {{ *frame = {{}}; }}
  if (frame == nullptr ||
      color_projectors == nullptr || projected == nullptr ||
      projector_count <= 0) {{
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }}
  if (!std::isfinite(alphas) || alphas < 0.0) {{
    return gra::mg5helas::EvaluationStatus::AmplitudeFailure;
  }}

  pars.setDependentParameters(alphas);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
{alpha_update}
  if (!setup_kinematics(lts, frame)) {{
    return gra::mg5helas::EvaluationStatus::KinematicsFailure;
  }}
  return CalcColorProjectedPrepared(
      color_projectors, projector_count, projected);
}}
"""


def make_matrix_tail(class_name: str, color_meta: dict[str, str | int]) -> str:
    matrix_name = color_meta["matrix_name"]
    ngraphs = color_meta["ngraphs"]
    denom = color_meta["denom"]
    cf = color_meta["cf"]
    color_flows = indent_block(str(color_meta["color_flows"]), 2)

    return f"""
double {class_name}::{matrix_name}() {{
  int i, j;
  // Local variables
  const int            ngraphs = {ngraphs};
  std::complex<double> ztemp;
  std::complex<double> jamp[ncolor];
  const auto           denom   = ColorDenominators();
  const auto           cf      = ColorMetric();

  calculate_color_flows(jamp);

  // Sum and square the color flows to get the matrix element
  double matrix = 0;
  for (i = 0; i < ncolor; i++) {{
    ztemp = 0.;
    for (j = 0; j < ncolor; j++) ztemp = ztemp + cf[i][j] * jamp[j];
    matrix = matrix + real(ztemp * conj(jamp[i])) / denom[i];
  }}

  // Store the leading color flows for choice of color
  for (i = 0; i < ncolor; i++) jamp2[0][i] += real(jamp[i] * conj(jamp[i]));

  return matrix;
}}

void {class_name}::calculate_color_flows(std::complex<double> jamp[ncolor]) const {{
{color_flows}
}}

std::vector<double> {class_name}::ColorDenominators() const {{ return {denom}; }}

std::vector<std::vector<double>> {class_name}::ColorMetric() const {{
  return {cf};
}}
"""


# Convert one standalone amplitude and its physical helicity evaluation paths
def transform_source(
    process: str,
    raw_source: str,
    color_meta: dict[str, str | int] | None,
    projection: str,
    model: str,
    alpha_zero: bool = True,
) -> str:
    class_name = f"AMP_MG5_{process}"
    photon_process = projection == "photon"
    incoming_pdg = 22 if photon_process else 21
    source = normalize_whitespace(raw_source)
    if photon_process:
        source = (
            f"// Generated MadGraph photon amplitude for {process}\n"
            "//\n"
            "// (c) 2026 Mikael Mieskolainen\n"
            "// Licensed under the MIT License <http://opensource.org/licenses/MIT>.\n\n" + source
        )
    helas_match = extract_or_fail(r'#include "(HelAmps_([^"]+))\.h"', raw_source)
    helas_header = f"{helas_match.group(1)}.h"
    process_directory = output_layout.projection_cpp_directory(projection)
    model_directory = output_layout.model_cpp_directory(model)
    parameter_class = f"Parameters_{helas_match.group(2)}"
    decay_chain = "// *   Decay:" in raw_source

    source = replace_or_fail(source, MESSAGE_ORIG, MESSAGE_NEW, 1)
    source = replace_or_fail(
        source,
        '#include "CPPProcess.h"',
        f"#include <cmath>\n#include <span>\n\n"
        f'#include "{output_layout.include_path(process_directory, f"{class_name}.h")}"',
        1,
    )
    source = replace_or_fail(
        source,
        f'#include "{helas_header}"',
        f'#include "{output_layout.include_path(model_directory, helas_header)}"',
        1,
    )
    source = replace_or_fail(source, "CPPProcess", class_name)

    source = replace_or_fail(
        source,
        f"pars = {parameter_class}::getInstance();",
        f"pars = {parameter_class}();  // GRANIITTI",
        1,
    )
    source = replace_or_fail(source, "pars->", "pars.")
    source = sub_or_fail(
        rf"void {class_name}::initProc\((?:std::)?string param_card_name\)\s*\{{",
        f"void {class_name}::initProc(std::string param_card_name) {{\n"
        "  InitParameters(SLHAReader(param_card_name));\n}\n\n"
        "// Initialize model parameters before constructing external wavefunctions\n"
        f"void {class_name}::InitParameters(SLHAReader slha)\n{{",
        source,
        count=1,
    )
    source = sub_or_fail(r"[ \t]*SLHAReader slha\(param_card_name\);[^\S\n]*\n", "", source, count=1)
    source = source.replace("  pars.setIndependentParameters(slha);",
                            "  pars.setIndependentParameters(slha);\n"
                            "  gra::mg5::ValidateModel(slha, pars.Particles());", 1)

    source = sub_or_fail(
        r"  static bool firsttime = true;\s*\n"
        r"  if \(firsttime\)\s*\n"
        r"  \{\s*\n"
        r"    pars\.printDependentParameters\(\);\s*\n"
        r"    pars\.printDependentCouplings\(\);\s*\n"
        r"    firsttime = false;\s*\n"
        r"  \}\s*\n",
        "",
        source,
        count=1,
    )
    source = sub_or_fail(
        r"  static bool goodhel\[ncomb\]\s*=\s*\{ncomb \* false\};",
        "  bool goodhel[ncomb] = {};",
        source,
        count=1,
    )
    source = replace_or_fail(
        source,
        "if (tsum != 0. && !goodhel[ihel])",
        "if (std::fpclassify(tsum) != FP_ZERO && !goodhel[ihel])",
        1,
    )
    source = sub_or_fail(
        r"  static int ntry\s*=\s*0,\s*sum_hel\s*=\s*0,\s*ngood\s*=\s*0;",
        "  int ntry = 0, sum_hel = 0, ngood = 0;",
        source,
        count=1,
    )
    source = sub_or_fail(
        r"  static int igood\[ncomb\];",
        "  int igood[ncomb + 1] = {};",
        source,
        count=1,
    )
    source = sub_or_fail(
        r"  static int jhel;",
        "  int jhel = 0;",
        source,
        count=1,
    )
    source = replace_or_fail(source, "static const int helicities", "const int helicities", 1)
    source = replace_or_fail(source, "ntry = ntry + 1;", "ntry = 1;  // GRANIITTI", 1)
    source = replace_or_fail(
        source, "setDependentParameters();", "setDependentParameters(alphas);", 1
    )

    source = replace_or_fail(
        source,
        "pars.printIndependentParameters();",
        "// pars.printIndependentParameters();  // GRANIITTI",
        1,
    )
    source = replace_or_fail(
        source,
        "pars.printIndependentCouplings();",
        "// pars.printIndependentCouplings();  // GRANIITTI",
        1,
    )

    source_signature = (
        f"gra::mg5helas::MatrixElementEvaluation {class_name}::Evaluate("
        "gra::LORENTZSCALAR &lts, double alphas"
    )
    source_signature += ", bool coherent_epa" if photon_process else ""
    source_signature += ")"
    source = sub_or_fail(
        rf"void {class_name}::sigmaKin\(\)",
        source_signature,
        source,
        count=1,
    )
    evaluation_guard = (
        r"\1\n  lts.hamp.clear();"
        r"\n  if (!std::isfinite(alphas) || alphas < 0.0) {"
        r" return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}; }"
    )
    if decay_chain:
        evaluation_guard += (
            r"\n  if (!gra::mg5::UsesFullDecayChainMode(lts)) {"
            r" return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0}; }"
        )
    source = sub_or_fail(
        rf"({re.escape(source_signature)}\s*\{{)",
        evaluation_guard,
        source,
        count=1,
    )

    if photon_process:
        source = sub_or_fail(
            r"  pars\.setDependentCouplings\(\);",
            (
                "  pars.setIndependentCouplings();\n"
                "  pars.setDependentCouplings();"
                + ("\n  if (gra::mg5helas::AlphaQEDAtZero(lts)) { pars.setAlphaQEDZero(); }" if alpha_zero else "")
            ),
            source,
            count=1,
        )

    source = replace_or_fail(
        source,
        "//--------------------------------------------------------------------------\n"
        "// Evaluate |M|^2, part independent of incoming flavour.\n",
        make_setup_kinematics(class_name, decay_chain, photon_process)
        + "\n//--------------------------------------------------------------------------\n"
        + "// Evaluate |M|^2, part independent of incoming flavour.\n",
        1,
    )

    source = sub_or_fail(
        r"  jamp2\[0\] = new double\[\d+\];",
        "  // Reinitialization must clear the fixed color-flow buffer\n  jamp2[0].fill(0.0);",
        source,
        count=1,
    )
    source = replace_or_fail(
        source,
        "  // Set external particle masses for this matrix element\n  mME.push_back",
        "  // Reinitialization must replace the external mass table\n"
        "  mME.clear();\n"
        "  // Set external particle masses for this matrix element\n"
        "  mME.push_back",
        1,
    )

    ncomb = int(extract_or_fail(r"const int ncomb = (\d+);", source).group(1))
    source = sub_or_fail(
        rf"  const int ncomb\s*=\s*{ncomb};",
        f"  const int ncomb = {ncomb}; \n"
        "  if (!setup_kinematics(lts)) {\n"
        "    return {gra::mg5helas::EvaluationStatus::KinematicsFailure, 0.0};\n"
        "  }",
        source,
        count=1,
    )

    source = replace_or_fail(
        source,
        "  }\n\n  if (sum_hel == 0 || ntry < 10)\n",
        "  }\n\n  goto SKIPLABEL;  // GRANIITTI: Skip this block\n  if (sum_hel == 0 || ntry < 10)\n",
        1,
    )

    tail_pattern = (
        r"  for \(int i = 0; i < nprocesses; i\+\+ \)\s*\n"
        r"\s*matrix_element\[i\]\s*/=\s*denominators\[i\];\s*\n"
        r"\s*\}\n"
    )
    if photon_process:
        if color_meta is None:
            fail(f"Photon process {process} has no generated color metric")
        amplitude_fill = (
            "  std::vector<gra::mg5helas::HelicityComponent> components;\n"
            "  components.reserve(ncomb * ncolor);\n"
            "  const auto color_denominators = ColorDenominators();\n"
            "  const auto color_metric = ColorMetric();\n"
            "  for (int ihel = 0; ihel < ncomb; ++ihel) {\n"
            "    calculate_wavefunctions(perm, helicities[ihel]);\n"
            "    std::complex<double> jamp[ncolor];\n"
            "    calculate_color_flows(jamp);\n"
            "    // Contract the JAMP color basis with the full finite-Nc metric\n"
            "    // [REFERENCE: Lifson and Mattelaer, Eur. Phys. J. C 82, 1144 (2022), arXiv:2210.07267]\n"
            "    const std::vector<std::complex<double>> color_amplitudes =\n"
            "        gra::mg5helas::ColorMetricAmplitudes(\n"
            "            std::vector<std::complex<double>>(jamp, jamp + ncolor),\n"
            "            color_denominators, color_metric);\n"
            "    if (color_amplitudes.size() != ncolor) {\n"
            "      lts.hamp.clear();\n"
            "      return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};\n"
            "    }\n"
            "    std::vector<int> outgoing;\n"
            "    for (int leg = ninitial; leg < nexternal; ++leg) { outgoing.push_back(helicities[ihel][leg]); }\n"
            "    for (int color = 0; color < ncolor; ++color) {\n"
            "      components.push_back({{helicities[ihel][0], helicities[ihel][1]},\n"
            "                            outgoing, static_cast<std::size_t>(color),\n"
            "                            color_amplitudes[color]});\n"
            "    }\n"
            "  }\n\n"
            "  lts.hamp.clear();\n"
            "  if (coherent_epa) {\n"
            "    const gra::M4Vec k1(p[0][1], p[0][2], p[0][3], p[0][0]);\n"
            "    const gra::M4Vec k2(p[1][1], p[1][2], p[1][3], p[1][0]);\n"
            "    lts.hamp = gra::mg5helas::ContractEPAPhotonSources(lts, components, k1, k2);\n"
            "  } else {\n"
            "    for (const auto &component : components) { lts.hamp.push_back(component.value); }\n"
            "  }\n"
        )
    else:
        if color_meta is None:
            fail(f"Colored process {process} has no generated color metric")
        amplitude_fill = (
            "  lts.hamp.clear();\n"
            "  lts.hamp.reserve(ncomb * ncolor);\n"
            "  const auto color_denominators = ColorDenominators();\n"
            "  const auto color_metric = ColorMetric();\n"
            "  for (int ihel = 0; ihel < ncomb; ++ihel) {\n"
            "    calculate_wavefunctions(perm, helicities[ihel]);\n"
            "    std::complex<double> jamp[ncolor];\n"
            "    calculate_color_flows(jamp);\n"
            "    // Contract the JAMP color basis with the full finite-Nc metric\n"
            "    // [REFERENCE: Lifson and Mattelaer, Eur. Phys. J. C 82, 1144 (2022), arXiv:2210.07267]\n"
            "    const auto color_amplitudes = gra::mg5helas::ColorMetricAmplitudes(\n"
            "        std::vector<std::complex<double>>(jamp, jamp + ncolor),\n"
            "        color_denominators, color_metric);\n"
            "    if (color_amplitudes.size() != ncolor) {\n"
            "      lts.hamp.clear();\n"
            "      return {gra::mg5helas::EvaluationStatus::AmplitudeFailure, 0.0};\n"
            "    }\n"
            "    for (const auto amplitude : color_amplitudes) {\n"
            "      lts.hamp.push_back(amplitude);\n"
            "    }\n"
            "  }\n"
        )
    amplitude_norm = (
        "  const double final_symmetry_factor =\n"
        "      gra::mg5helas::AppliedFinalStateSymmetryFactor(lts);\n"
        "  const double normalization = std::sqrt(\n"
        "      final_symmetry_factor / static_cast<double>(denominators[0]));\n"
        "  gra::Scale(lts.hamp, normalization);\n"
        "  const double amp2 = gra::SquaredNorm(lts.hamp);\n"
        "  // Screening applies the incoming photon spin average after coherent summation\n"
        "  gra::Scale(lts.hamp, 2.0);\n"
        if photon_process
        else (
            "  const double amp2 = gra::SquaredNorm(lts.hamp) / denominators[0];  "
            "// spin average matrix element squared\n"
        )
    )

    tail_replacement = (
        "  for (int i = 0; i < nprocesses; i++ )\n"
        "    matrix_element[i] /= denominators[i]; \n\n"
        "SKIPLABEL:\n\n"
        "  // @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
        "  // Define permutation\n"
        "  for (int i = 0; i < nexternal; ++i) { perm[i] = i; }\n\n"
        "  // Loop over helicity combinations in the fixed GRANIITTI basis\n"
        + amplitude_fill
        + "\n"
        "  // Total amplitude squared over all helicity combinations individually\n"
        + amplitude_norm
        + "\n"
        + "  return {gra::mg5helas::EvaluationStatus::Success, amp2};\n"
        "                // @@@@@@@@@@@@@@@@@@@@@@@@@@ GRANIITTI @@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n"
        "}\n"
    )
    source = sub_or_fail(tail_pattern, tail_replacement, source, count=1)

    if color_meta is not None:
        helicity_block = extract_or_fail(
            r"((?:static )?const int helicities\[ncomb\]\[nexternal\] = \{.*?\};)",
            source,
            flags=re.S,
        ).group(1)

        projection_helpers = make_calc_color_flow_helicity(
            class_name, ncomb, helicity_block, int(color_meta["ncolor"]), incoming_pdg, alpha_zero
        ) + make_calc_color_projected_helicity(class_name, ncomb, helicity_block, incoming_pdg, alpha_zero)
        if photon_process:
            projection_helpers += make_calc_color_projected_epa_hard(
                class_name, alpha_zero
            ) + make_calc_color_projected_prepared(
                class_name, helicity_block
            )

        source = replace_or_fail(
            source,
            "//--------------------------------------------------------------------------\n"
            "// Evaluate |M|^2, including incoming flavour dependence.\n",
            projection_helpers
            + "\n//--------------------------------------------------------------------------\n"
            + "// Evaluate |M|^2, including incoming flavour dependence.\n",
            1,
        )

        matrix_anchor = f"double {class_name}::{color_meta['matrix_name']}()"
        start = source.find(matrix_anchor)
        if start == -1:
            fail(f"Could not locate matrix function {matrix_anchor}")
        source = source[:start] + make_matrix_tail(class_name, color_meta)

    return narrow_std_namespace(source)
