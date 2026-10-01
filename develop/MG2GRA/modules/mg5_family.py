#!/usr/bin/env python3
#
# Transform multi-subprocess MG5 exports into isolated GRANIITTI families
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re
import sys
from contextlib import chdir
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

from core.io.files import ensure_dir

from . import output_layout
from .cpp_support import cpp_code_only, function_block
from .mg5_color import coupling_orders
from .mg5_support import (
    narrow_std_namespace,
    transform_helas_header,
    transform_helas_source,
    transform_parameters_header,
    transform_parameters_source,
)
from .process_registry import (
    FamilyGenerationSetup,
    MG5Source,
    ParsedProcessSyntax,
    ProcessNode,
    family_generation_setup_data,
    family_generation_setup_from_data,
    mg5_source_data,
    mg5_source_from_data,
    parse_process_clause,
    parse_process_syntax,
    stable_tokens,
    standard_particle_selectors,
)


# Replace one generated anchor only when its multiplicity is exact
def replace_or_fail(text: str, old: str, new: str, count: int = 1) -> str:
    occurrences = text.count(old)
    if occurrences != count:
        raise RuntimeError(
            f"Expected {count} generated family anchor(s), found {occurrences}: {old}"
        )
    return text.replace(old, new, count)


# Replace one generated regular-expression anchor only when it is present
def sub_or_fail(pattern: str, repl: Any, text: str, *, count: int = 1) -> str:
    updated, matches = re.subn(pattern, repl, text, count=count)
    if matches != count:
        raise RuntimeError(
            f"Expected {count} generated family pattern(s), found {matches}: {pattern}"
        )
    return updated


@dataclass(frozen=True)
class ChannelTopologyNode:
    """One exact numeric generated production or decay node"""

    pdg: int
    daughters: tuple[ChannelTopologyNode, ...] = ()


@dataclass(frozen=True)
class Channel:
    """One exact incoming flavour and generated final-state topology"""

    initial: tuple[int, int]
    final: tuple[int, ...]
    topology: tuple[ChannelTopologyNode, ...]
    external_color_representations: tuple[int, ...] = ()
    external_color_flows: tuple[tuple[tuple[int, int], ...], ...] = ()


@dataclass(frozen=True)
class ChannelCrossing:
    """One exact channel and its generated beam orientation"""

    channel: Channel
    mirrored: bool


@dataclass(frozen=True)
class ExternalColorStructure:
    """Exact external signature, representations, and raw JAMP color rows"""

    signature: tuple[int, ...]
    representations: tuple[int, ...]
    flows: tuple[tuple[tuple[int, int], ...], ...]


@dataclass(frozen=True)
class SubprocessColorStructure:
    """MG Python color structure for one standalone C++ subprocess export"""

    generated_name: str
    canonical: tuple[ExternalColorStructure, ...]
    mirrored: tuple[ExternalColorStructure, ...]
    orders: tuple[tuple[tuple[str, int], ...], ...] = ()


@dataclass(frozen=True)
class SubprocessChannels:
    """Complete exact channels for one transformed subprocess"""

    generated_name: str
    channels: tuple[Channel, ...]


@dataclass(frozen=True)
class FamilyChannels:
    """Complete persistent channels for one generated family"""

    family: str
    generation_setup: FamilyGenerationSetup
    mg5_source: MG5Source
    subprocesses: tuple[SubprocessChannels, ...]


@dataclass(frozen=True)
class FamilySubprocess:
    """One transformed subprocess and its exact channels"""

    class_name: str
    header: str
    source: str
    channels: tuple[Channel, ...]


# Reuse exact MG5 particle tokens from the process registry schema
PDG_NAMES = {
    token: selector[0]
    for token, selector in standard_particle_selectors().items()
    if len(selector) == 1
}


# Wrap a transformed header declaration in one process-family namespace
def wrap_header_namespace(text: str, family: str, declaration: str) -> str:
    start = text.find(declaration)
    if start < 0:
        raise RuntimeError(f"Could not find family declaration {declaration}")
    endif = text.rfind("#endif")
    if endif < 0:
        raise RuntimeError("Could not find generated header guard terminator")
    text = text[:start] + f"namespace {family} {{\n\n" + text[start:]
    endif += len(f"namespace {family} {{\n\n")
    return text[:endif] + f"}}  // namespace {family}\n\n" + text[endif:]


# Wrap transformed source definitions in one process-family namespace
def wrap_source_namespace(text: str, family: str, declaration: str) -> str:
    start = text.find(declaration)
    if start < 0:
        raise RuntimeError(f"Could not find family definition {declaration}")
    return (
        text[:start]
        + f"namespace {family} {{\n\n"
        + text[start:]
        + (f"\n}}  // namespace {family}\n")
    )


# Replace a generated include guard with a repository-unique family guard
def replace_include_guard(text: str, guard: str) -> str:
    text, first = re.subn(r"^#ifndef\s+\S+", f"#ifndef {guard}", text, count=1, flags=re.M)
    text, second = re.subn(r"^#define\s+\S+", f"#define {guard}", text, count=1, flags=re.M)
    if first != 1 or second != 1:
        raise RuntimeError("Could not replace generated include guard")
    return text


# Transform and isolate one family's generated HELAS declarations
def transform_family_helas_header(
    raw_header: str, family: str, helas_header: str, model_suffix: str
) -> str:
    header = transform_helas_header(raw_header)
    guard = f"GRANIITTI_AMPLITUDE_{family}_{helas_header.replace('.', '_').upper()}"
    header = replace_include_guard(header, guard)
    return wrap_header_namespace(header, family, f"namespace MG5_{model_suffix}")


# Transform and isolate one family's generated HELAS definitions
def transform_family_helas_source(
    raw_source: str,
    family: str,
    helas_header: str,
    model_suffix: str,
    include_directory: str,
) -> str:
    source = transform_helas_source(raw_source, helas_header, include_directory)
    return wrap_source_namespace(source, family, f"namespace MG5_{model_suffix}")


# Transform and isolate one family's generated parameter declarations
def transform_family_parameters_header(
    raw_header: str,
    family: str,
    parameter_class: str,
    parameter_header: str,
    alpha_zero: bool,
) -> str:
    header = transform_parameters_header(raw_header, parameter_class, alpha_zero)
    guard = f"GRANIITTI_AMPLITUDE_{family}_{parameter_header.replace('.', '_').upper()}"
    header = replace_include_guard(header, guard)
    return wrap_header_namespace(header, family, f"class {parameter_class}")


# Transform and isolate one family's generated parameter definitions
def transform_family_parameters_source(
    raw_source: str,
    family: str,
    parameter_class: str,
    parameter_header: str,
    charge: str | None,
    charge_square: str | None,
    inverse_alpha: float,
    include_directory: str,
) -> str:
    source = transform_parameters_source(
        raw_source,
        parameter_class,
        parameter_header,
        charge,
        charge_square,
        inverse_alpha,
        include_directory,
    )
    return wrap_source_namespace(
        source, family, f"void {parameter_class}::setIndependentParameters"
    )


# Compute each generated production comment with its following decay comments
def process_comment_blocks(raw_source: str) -> list[tuple[str, tuple[str, ...]]]:
    blocks: list[tuple[str, tuple[str, ...]]] = []
    production: str | None = None
    decays: list[str] = []
    for line in raw_source.splitlines():
        process_match = re.match(r"^// Process:\s*(.+)$", line)
        if process_match is not None:
            if production is not None:
                blocks.append((production, tuple(decays)))
            production = process_match.group(1)
            decays = []
            continue
        decay_match = re.match(r"^// \*\s+Decay:\s*(.+)$", line)
        if decay_match is not None and production is not None:
            decays.append(decay_match.group(1))
    if production is not None:
        blocks.append((production, tuple(decays)))
    if not blocks:
        raise RuntimeError("No generated production process comments found")
    return blocks


# Compute generated initial particles and parsed final topology in MG5 order
def production_processes(
    raw_source: str,
) -> list[tuple[tuple[str, str], ParsedProcessSyntax]]:
    processes: list[tuple[tuple[str, str], ParsedProcessSyntax]] = []
    for production, decays in process_comment_blocks(raw_source):
        initial, _ = parse_process_clause(production)
        if len(initial) != 2:
            raise RuntimeError("Generated process does not have two incoming particles")
        parsed = parse_process_syntax(", ".join((production, *decays)))
        processes.append(((initial[0], initial[1]), parsed))
    return processes


# Convert one generated branch to an exact numeric channel topology node
def channel_topology_node(node: ProcessNode, pdg_names: dict[str, int]) -> ChannelTopologyNode:
    try:
        pdg = pdg_names[node.particle]
    except KeyError as error:
        raise RuntimeError(f"Unsupported generated particle name {error.args[0]}") from error
    return ChannelTopologyNode(
        pdg,
        tuple(channel_topology_node(daughter, pdg_names) for daughter in node.daughters),
    )


# Convert one parsed process to an exact numeric channel topology
def channel_topology(
    parsed: ParsedProcessSyntax, pdg_names: dict[str, int]
) -> tuple[ChannelTopologyNode, ...]:
    return tuple(channel_topology_node(node, pdg_names) for node in parsed.production)


# Convert generated particle names into exact subprocess channels
def channel_crossings(
    raw_source: str,
    particle_pdgs: dict[str, int] | None = None,
) -> tuple[ChannelCrossing, ...]:
    channels: list[ChannelCrossing] = []
    pdg_names = {**PDG_NAMES, **(particle_pdgs or {})}
    mirror = "Mirror initial state momenta" in raw_source
    for initial_names, parsed in production_processes(raw_source):
        try:
            initial = (pdg_names[initial_names[0]], pdg_names[initial_names[1]])
            stable_final = tuple(pdg_names[name] for name in stable_tokens(parsed))
        except KeyError as error:
            raise RuntimeError(f"Unsupported generated particle name {error.args[0]}") from error
        topology = channel_topology(parsed, pdg_names)
        channel = Channel(initial, stable_final, topology)
        crossing = ChannelCrossing(channel, False)
        if all(existing.channel != channel for existing in channels):
            channels.append(crossing)
        mirrored = Channel((initial[1], initial[0]), stable_final, topology)
        mirror_crossing = ChannelCrossing(mirrored, True)
        if mirror and all(existing.channel != mirrored for existing in channels):
            channels.append(mirror_crossing)
    return tuple(channels)


# Compute source derived exact channels without external color structure
def channels(
    raw_source: str,
    particle_pdgs: dict[str, int] | None = None,
) -> tuple[Channel, ...]:
    return tuple(crossing.channel for crossing in channel_crossings(raw_source, particle_pdgs))


# Compute physical external legs in generated momentum order
def external_legs(process: Any) -> list[Any]:
    legs = list(process.get_legs_with_decays())
    if len(legs) < 2:
        raise RuntimeError("Generated process has fewer than two external legs")
    return legs


# Compute one exact external PDG signature in generated momentum order
def external_signature(process: Any) -> tuple[int, ...]:
    return tuple(int(leg.get("id")) for leg in external_legs(process))


# Compute signed external color representations in physical leg order
def external_color_representations(process: Any, legs: list[Any]) -> tuple[int, ...]:
    model = process.get("model")
    return tuple(int(model.get_particle(int(leg.get("id"))).get_color()) for leg in legs)


# Convert one MG color decomposition into compact physical external rows
def external_color_flows(
    matrix_element: Any,
    process: Any,
    legs: list[Any],
) -> tuple[tuple[tuple[int, int], ...], ...]:
    basis = matrix_element.get("color_basis")
    if not basis:
        return (tuple((0, 0) for _ in legs),)
    representations = {
        int(leg.get("number")): representation
        for leg, representation in zip(
            legs,
            external_color_representations(process, legs),
            strict=True,
        )
    }
    decomposition = basis.color_flow_decomposition(representations, 2)
    flows = tuple(
        tuple(tuple((1 if tag >= 0 else -1) * (abs(int(tag)) % 500)
                    for tag in flow[int(leg.get("number"))]) for leg in legs)
        for flow in decomposition
    )
    if len(flows) != len(basis):
        raise RuntimeError("MadGraph color decomposition and matrix element JAMP counts disagree")
    return flows


# Validate grouped exact processes share the exported color representation layout
def validate_grouped_color_layout(
    processes: list[Any],
    template_numbers: tuple[int, ...],
    template_representations: tuple[int, ...],
) -> None:
    for process in processes:
        legs = external_legs(process)
        numbers = tuple(int(leg.get("number")) for leg in legs)
        representations = external_color_representations(process, legs)
        if numbers != template_numbers or representations != template_representations:
            raise RuntimeError(
                "Grouped MadGraph processes have incompatible external color layouts"
            )


# Encode the standalone C++ subprocess directory name from one MG process object
def exported_subprocess_name(matrix_element: Any) -> str:
    process = matrix_element.get("processes")[0]
    process_string = process.base_string().replace(" ", "")
    process_string = process_string.replace(">", "_")
    process_string = process_string.replace("+", "p")
    process_string = process_string.replace("-", "m")
    process_string = process_string.replace("~", "x")
    process_string = process_string.replace("/", "_no_")
    process_string = process_string.replace("$", "_nos_")
    process_string = process_string.replace("|", "_or_")
    model_name = process.get("model").get("name").replace("-", "_")
    model_name = model_name.replace("+", "_plus_")
    return f"P{int(process.get('id'))}_Sigma_{model_name}_{process_string}"


# Extract exact canonical and mirror color structures from one MG matrix element
def matrix_element_color_structure(matrix_element: Any) -> SubprocessColorStructure:
    canonical = list(matrix_element.get("processes"))
    if not canonical:
        raise RuntimeError("MadGraph matrix element has no exact processes")
    template_process = canonical[0]
    canonical_legs = external_legs(template_process)
    canonical_numbers = tuple(int(leg.get("number")) for leg in canonical_legs)
    canonical_representations = external_color_representations(template_process, canonical_legs)
    validate_grouped_color_layout(canonical, canonical_numbers, canonical_representations)
    canonical_flows = external_color_flows(matrix_element, template_process, canonical_legs)
    canonical_structures = tuple(
        ExternalColorStructure(
            external_signature(process),
            canonical_representations,
            canonical_flows,
        )
        for process in canonical
    )

    mirrored = list(matrix_element.get_mirror_processes())
    mirror_structures: tuple[ExternalColorStructure, ...] = ()
    if mirrored:
        mirror_legs = list(canonical_legs)
        mirror_legs[0:2] = [mirror_legs[1], mirror_legs[0]]
        mirror_numbers = tuple(int(leg.get("number")) for leg in mirror_legs)
        mirror_representations = external_color_representations(template_process, mirror_legs)
        validate_grouped_color_layout(mirrored, mirror_numbers, mirror_representations)
        mirror_flows = external_color_flows(matrix_element, template_process, mirror_legs)
        mirror_structures = tuple(
            ExternalColorStructure(
                external_signature(process),
                mirror_representations,
                mirror_flows,
            )
            for process in mirrored
        )
    return SubprocessColorStructure(
        exported_subprocess_name(matrix_element),
        canonical_structures,
        mirror_structures,
        tuple(tuple(row.items()) for row in coupling_orders(matrix_element)),
    )


# Load the MadGraph command and HELAS APIs from one installation
def load_madgraph_family_api(mg5_root: Path) -> tuple[Any, Any]:
    root = str(mg5_root.resolve())
    if root not in sys.path:
        sys.path.insert(0, root)
    from madgraph.core import helas_objects
    from madgraph.interface.master_interface import MasterCmd

    return MasterCmd, helas_objects


# Generate exact family color structures directly from the MadGraph Python API
def generate_family_color_structure(
    mg5_root: Path,
    work_dir: Path,
    model_import: str,
    definitions: list[str],
    processes: list[str],
    complex_mass_scheme: bool = False,
) -> tuple[SubprocessColorStructure, ...]:
    if not processes:
        raise RuntimeError("MadGraph family color generation setup requires a process")
    ensure_dir(work_dir)
    with chdir(work_dir):
        MasterCmd, helas_objects = load_madgraph_family_api(mg5_root)
        command = MasterCmd()
        command.no_notification()
        command.exec_cmd(
            "set complex_mass_scheme True --allow_qed"
            if complex_mass_scheme
            else "set complex_mass_scheme False",
            printcmd=False,
            precmd=True,
            postcmd=True,
        )
        command.exec_cmd(f"import model {model_import}", printcmd=False, precmd=True, postcmd=True)
        for definition in definitions:
            command.exec_cmd(definition, printcmd=False, precmd=True, postcmd=True)
        command.exec_cmd(f"generate {processes[0]}", printcmd=False, precmd=True, postcmd=True)
        for process in processes[1:]:
            command.exec_cmd(f"add process {process}", printcmd=False, precmd=True, postcmd=True)
        matrix_elements = helas_objects.HelasMultiProcess.generate_matrix_elements(
            command._curr_amps
        )
    color_structures = tuple(
        matrix_element_color_structure(matrix_element) for matrix_element in matrix_elements
    )
    names = [subprocess.generated_name for subprocess in color_structures]
    if len(names) != len(set(names)):
        raise RuntimeError("MadGraph generated duplicate standalone subprocess names")
    return color_structures


# Compute exact signature-indexed color rows with duplicate checks
def indexed_external_color_structure(
    rows: tuple[ExternalColorStructure, ...],
) -> dict[tuple[int, ...], ExternalColorStructure]:
    indexed: dict[tuple[int, ...], ExternalColorStructure] = {}
    for row in rows:
        if row.signature in indexed:
            raise RuntimeError("Duplicate exact MadGraph external color signature")
        indexed[row.signature] = row
    return indexed


# Compute subprocess color structures indexed by export name
def indexed_subprocess_color_structure(
    rows: tuple[SubprocessColorStructure, ...],
) -> dict[str, SubprocessColorStructure]:
    indexed: dict[str, SubprocessColorStructure] = {}
    for row in rows:
        if row.generated_name in indexed:
            raise RuntimeError("Duplicate MadGraph standalone subprocess color structure")
        indexed[row.generated_name] = row
    return indexed


# Attach authoritative external color rows to source-derived topology channels
def complete_channels(
    raw_source: str,
    color_structure: SubprocessColorStructure,
    particle_pdgs: dict[str, int] | None = None,
    source_class_name: str = "CPPProcess",
) -> tuple[Channel, ...]:
    sources = channel_crossings(raw_source, particle_pdgs)
    canonical = indexed_external_color_structure(color_structure.canonical)
    mirrored = indexed_external_color_structure(color_structure.mirrored)
    expected_canonical = {
        crossing.channel.initial + crossing.channel.final
        for crossing in sources
        if not crossing.mirrored
    }
    expected_mirrored = {
        crossing.channel.initial + crossing.channel.final
        for crossing in sources
        if crossing.mirrored
    }
    if expected_canonical != set(canonical):
        raise RuntimeError(
            f"Source and MadGraph canonical signatures disagree for "
            f"{color_structure.generated_name}"
        )
    if expected_mirrored != set(mirrored):
        raise RuntimeError(
            f"Source and MadGraph mirror signatures disagree for {color_structure.generated_name}"
        )

    channels = []
    for crossing in sources:
        signature = crossing.channel.initial + crossing.channel.final
        color = (mirrored if crossing.mirrored else canonical)[signature]
        channels.append(
            replace(
                crossing.channel,
                external_color_representations=color.representations,
                external_color_flows=color.flows,
            )
        )
    ncolor, _, _, _ = color_data(raw_source, source_class_name)
    if any(len(channel.external_color_flows) != ncolor for channel in channels):
        raise RuntimeError(
            f"MadGraph API and standalone matrix element JAMP counts disagree for "
            f"{color_structure.generated_name}"
        )
    return tuple(channels)


# Compute the generated number of outgoing stable external particles
def generated_stable_final_count(raw_header: str) -> int:
    values: dict[str, int] = {}
    for name in ("ninitial", "nexternal"):
        match = re.search(
            rf"\bstatic\s+const\s+int\s+{name}\s*=\s*(\d+)\s*;",
            raw_header,
        )
        if match is None:
            raise RuntimeError(f"Could not find generated {name} constant")
        values[name] = int(match.group(1))
    if values["nexternal"] < values["ninitial"]:
        raise RuntimeError("Generated external particle count is smaller than initial count")
    return values["nexternal"] - values["ninitial"]


# Make generated sigmaKin scratch state local to one call
def localize_process_scratch(source: str) -> str:
    source, firsttime_blocks = re.subn(
        r"  static bool firsttime = true;\s*"
        r"if \(firsttime\)\s*\{\s*"
        r"pars\.printDependentParameters\(\);\s*"
        r"pars\.printDependentCouplings\(\);\s*"
        r"firsttime = false;\s*\}",
        "",
        source,
    )
    if firsttime_blocks > 1:
        raise RuntimeError("Expected at most one generated first-time print block")
    scratch_replacements = (
        ("static bool goodhel[ncomb] = {ncomb * false};", "bool goodhel[ncomb] = {};"),
        (
            "static int ntry = 0, sum_hel = 0, ngood = 0;",
            "int ntry = 0, sum_hel = 0, ngood = 0;",
        ),
        ("static int igood[ncomb];", "int igood[ncomb + 1] = {};"),
        ("static int jhel;", "int jhel = 0;"),
        ("static const int helicities", "const int helicities"),
    )
    if any(old in source for old, _ in scratch_replacements):
        for old, new in scratch_replacements:
            count = source.count(old)
            if count != 1:
                raise RuntimeError(f"Expected one generated family anchor, found {count}: {old}")
            source = source.replace(old, new, 1)

    print_replacements = (
        ("pars.printIndependentParameters();", "// pars.printIndependentParameters();"),
        ("pars.printIndependentCouplings();", "// pars.printIndependentCouplings();"),
        ("pars.printDependentParameters();", "// pars.printDependentParameters();"),
        ("pars.printDependentCouplings();", "// pars.printDependentCouplings();"),
    )
    for old, new in print_replacements:
        count = source.count(old)
        if count > 1:
            raise RuntimeError(
                f"Expected at most one generated family anchor, found {count}: {old}"
            )
        if count == 1:
            source = source.replace(old, new, 1)
    return source


# Transform one generated family subprocess header
def transform_subprocess_header(
    raw_header: str,
    family: str,
    class_name: str,
    parameter_header: str,
    include_directory: str,
) -> str:
    header = raw_header.replace("\r\n", "\n")
    cpp_process_count = header.count("CPPProcess")
    if cpp_process_count == 0:
        raise RuntimeError(f"Could not find generated class anchors for {class_name}")
    header = header.replace("CPPProcess", class_name)
    guard = f"GRANIITTI_AMPLITUDE_{family}_{class_name}_H"
    header = replace_include_guard(header, guard)
    header = replace_or_fail(
        header,
        f'#include "{parameter_header}"',
        f'#include "{output_layout.include_path(include_directory, parameter_header)}"\n'
        f'#include "{output_layout.include_path(include_directory, "ProcessBase.h")}"',
    )
    header = replace_or_fail(
        header, f"class {class_name}", f"class {class_name} : public ProcessBase"
    )
    header = sub_or_fail(
        rf"{re.escape(class_name)}\(\)\s*\{{\s*\}}",
        (f"{class_name}() = default;\n    ~{class_name}() override = default;"),
        header,
    )
    header = replace_or_fail(
        header, "virtual void sigmaKin();", "virtual void sigmaKin() override;"
    )
    header = replace_or_fail(
        header, "virtual double sigmaHat();", "virtual double sigmaHat() override;"
    )
    header = replace_or_fail(
        header,
        "const vector<double> & getMasses() const {return mME;}",
        "const vector<double> & getMasses() const override {return mME;}",
    )
    header = replace_or_fail(
        header,
        "void setMomenta(vector < double * > & momenta){p = momenta;}",
        "void setMomenta(vector < double * > & momenta) override {p = momenta;}",
    )
    header = replace_or_fail(
        header,
        "void setInitial(int inid1, int inid2){id1 = inid1; id2 = inid2;}",
        "void setInitial(int inid1, int inid2) override {id1 = inid1; id2 = inid2;}\n"
        "    void setAlphaS(double in) override {alphaS = in;}",
    )
    header = re.sub(
        r"void initProc\((?:std::)?string param_card_name\);",
        "void initProc(std::string param_card_name);\n"
        "    // Initialize every mass and coupling through the generated model\n"
        "    void InitParameters(SLHAReader slha) override;\n"
        "    // Compute evaluated model particle parameters\n"
        "    gra::mg5::ParticleMap Particles() const override { return pars.Particles(); }\n"
        "    // Compute the evaluated UFO electromagnetic coupling\n"
        "    double AlphaQED() const override { return pars.AlphaQED(); }",
        header,
        count=1,
    )
    header = sub_or_fail(
        r"Parameters_\w+\s*\*\s*pars\s*;",
        lambda match: match.group(0).replace("*", "").replace(";", ";  // GRANIITTI"),
        header,
    )
    header = sub_or_fail(
        r"double\s*\*\s*jamp2\[nprocesses\]\s*;",
        "std::vector<std::vector<double>> jamp2 =\n"
        "        std::vector<std::vector<double>>(nprocesses);",
        header,
    )
    header = replace_or_fail(
        header,
        "    // Initial particle ids",
        "    // Event-dependent strong coupling\n"
        "    double alphaS = 0.118;\n\n"
        "    // Initial particle ids",
    )
    header = wrap_header_namespace(header, family, f"class {class_name}")
    return narrow_std_namespace(header)


# Transform one generated family subprocess source
def transform_subprocess_source(
    raw_source: str,
    family: str,
    class_name: str,
    helas_header: str,
    parameter_class: str,
    model_suffix: str,
    alpha_zero: bool,
    include_directory: str,
) -> str:
    source = raw_source.replace("\r\n", "\n")
    cpp_process_count = source.count("CPPProcess")
    if cpp_process_count == 0:
        raise RuntimeError(f"Could not find generated source anchors for {class_name}")
    source = source.replace("CPPProcess", class_name)
    source = replace_or_fail(
        source,
        f'#include "{class_name}.h"',
        f"#include <cmath>\n\n"
        f'#include "{output_layout.include_path(include_directory, f"{class_name}.h")}"',
    )
    source = replace_or_fail(
        source,
        f'#include "{helas_header}"',
        f'#include "{output_layout.include_path(include_directory, helas_header)}"\n'
        '#include "Graniitti/Particle/MForm.h"',
    )
    source = sub_or_fail(
        rf"pars\s*=\s*{re.escape(parameter_class)}::getInstance\(\)\s*;",
        f"pars = {parameter_class}();",
        source,
    )
    pars_pointer_count = source.count("pars->")
    if pars_pointer_count == 0:
        raise RuntimeError(f"Could not find generated parameter access for {class_name}")
    source = source.replace("pars->", "pars.")
    source = sub_or_fail(
        rf"void {class_name}::initProc\((?:std::)?string param_card_name\)\s*\{{",
        f"void {class_name}::initProc(std::string param_card_name) {{\n"
        "  InitParameters(SLHAReader(param_card_name));\n}\n\n"
        "// Initialize model parameters before constructing external wavefunctions\n"
        f"void {class_name}::InitParameters(SLHAReader slha)\n{{",
        source,
        count=1,
    )
    source = sub_or_fail(r"[ \t]*SLHAReader slha\(param_card_name\);[^\S\n]*\n", "", source)
    source = source.replace("  pars.setIndependentParameters(slha);",
                            "  pars.setIndependentParameters(slha);\n"
                            "  gra::mg5::ValidateModel(slha, pars.Particles());", 1)

    if "if (tsum != 0. && !goodhel[ihel])" in source:
        source = replace_or_fail(
            source,
            "if (tsum != 0. && !goodhel[ihel])",
            "if (std::fpclassify(tsum) != FP_ZERO && !goodhel[ihel])",
        )
    if source.count("pars.setDependentParameters();") != 1:
        raise RuntimeError(f"Could not replace event-dependent parameters for {class_name}")
    source = source.replace(
        "pars.setDependentParameters();", "pars.setDependentParameters(alphaS);", 1
    )
    if alpha_zero:
        source = replace_or_fail(
            source,
            "pars.setDependentCouplings();",
            "pars.setIndependentCouplings();\n  pars.setDependentCouplings();\n  if (alphaQEDZero()) { pars.setAlphaQEDZero(); }",
        )
    if "jamp2" in source:
        if source.count("jamp2[0] = new double") != 1:
            raise RuntimeError(f"Could not replace the color-flow buffer for {class_name}")
        source = sub_or_fail(
            r"jamp2\[0\] = new double\[(\d+)\];",
            r"jamp2[0].assign(\1, 0.0);",
            source,
        )
    mass_anchor = "  // Set external particle masses for this matrix element\n  mME.push_back"
    if "mME.push_back" in source:
        if source.count(mass_anchor) != 1:
            raise RuntimeError(f"Could not replace the external mass table for {class_name}")
        source = source.replace(
            mass_anchor,
            "  // Reinitialization must replace the external mass table\n"
            "  mME.clear();\n"
            "  // Set external particle masses for this matrix element\n"
            "  mME.push_back",
            1,
        )
    std_namespace_count = source.count("using namespace std;")
    if std_namespace_count > 1:
        raise RuntimeError("Generated family source has duplicate std namespace imports")
    source = source.replace("using namespace std;", "", 1)
    source = replace_or_fail(source, f"using namespace MG5_{model_suffix};", "")
    source = localize_process_scratch(source)
    marker = "//==========================================================================\n// Class member functions"
    if marker not in source:
        raise RuntimeError(f"Could not find class implementation marker for {class_name}")
    source = replace_or_fail(
        source, marker, f"namespace {family} {{\n\nusing namespace MG5_{model_suffix};\n\n{marker}"
    )
    source += f"\n}}  // namespace {family}\n"
    return narrow_std_namespace(source)


# Transform one matrix_elements generated subprocess pair
def transform_subprocess(
    raw_header: str,
    raw_source: str,
    family: str,
    generated_name: str,
    parameter_header: str,
    helas_header: str,
    parameter_class: str,
    model_suffix: str,
    alpha_zero: bool,
    include_directory: str,
    color_structure: SubprocessColorStructure | None = None,
    particle_pdgs: dict[str, int] | None = None,
) -> FamilySubprocess:
    class_name = f"{family}_{generated_name}"
    source_channels = channels(raw_source, particle_pdgs)
    expected_stable_finals = generated_stable_final_count(raw_header)
    for channel in source_channels:
        if len(channel.final) != expected_stable_finals:
            raise RuntimeError(
                f"Stable final count {len(channel.final)} does not match generated "
                f"external count {expected_stable_finals} for {class_name}"
            )
    if color_structure is None:
        raise RuntimeError(f"Missing MadGraph color structure for {class_name}")
    if color_structure.generated_name != generated_name:
        raise RuntimeError(
            f"MadGraph color structure name {color_structure.generated_name} does not "
            f"match generated subprocess {generated_name}"
        )
    exact_channels = complete_channels(raw_source, color_structure, particle_pdgs)
    validate_selected_multiplicities(raw_source, "CPPProcess", exact_channels)
    header = transform_subprocess_header(
        raw_header, family, class_name, parameter_header, include_directory
    )
    for name, order in (("AlphaSPower", "QCD"), ("AlphaQEDPower", "QED")):
        power = max((dict(row).get(order, 0) for row in color_structure.orders), default=0)
        header = header.replace("public:", f"public:\n    // Compute the maximum generated amplitude order\n    int {name}() const override {{ return {power}; }}", 1)
    return FamilySubprocess(
        class_name=class_name,
        header=header,
        source=transform_subprocess_source(
            raw_source,
            family,
            class_name,
            helas_header,
            parameter_class,
            model_suffix,
            alpha_zero,
            include_directory,
        ),
        channels=exact_channels,
    )


# Generate the common matrix element base class for one family
def process_base_header(family: str) -> str:
    guard = f"GRANIITTI_AMPLITUDE_{family}_PROCESSBASE_H"
    return f"""#ifndef {guard}
#define {guard}

#include <string>
#include <vector>

#include "@HELICITY_INCLUDE@"
#include "Graniitti/Amplitude/MG5/Runtime/AMP_MG5_Model.h"

namespace {family} {{

// Common base class for one generated MG5 subprocess
class ProcessBase {{
 public:
  // Construct one generated subprocess matrix element
  ProcessBase() = default;
  ProcessBase(const ProcessBase &) = delete;
  ProcessBase &operator=(const ProcessBase &) = delete;
  ProcessBase(ProcessBase &&) = delete;
  ProcessBase &operator=(ProcessBase &&) = delete;

  // Destroy one generated subprocess instance
  virtual ~ProcessBase() {{}}

  // Initialize this subprocess from the shared parameter card
  virtual void initProc(std::string param_card_name) = 0;

  // Evaluate the kinematic matrix element for the current momenta
  virtual void sigmaKin() = 0;

  // Compute the flavour-filtered matrix element
  virtual double sigmaHat() = 0;

  // Compute spin-color averaged complex helicity amplitudes
  virtual std::vector<gra::mg5helas::HelicityComponent> helicityAmplitudes() = 0;

  // Initialize masses and couplings together from one model card
  virtual void InitParameters(SLHAReader slha) = 0;

  // Compute evaluated model particle parameters
  virtual gra::mg5::ParticleMap Particles() const = 0;

  // Compute the evaluated UFO electromagnetic coupling
  virtual double AlphaQED() const = 0;

  // Compute the maximum generated amplitude orders
  virtual int AlphaSPower() const = 0;
  virtual int AlphaQEDPower() const = 0;

  // Compute the current external mass table
  virtual const std::vector<double> &getMasses() const = 0;

  // Set the external momenta in MadGraph convention
  virtual void setMomenta(std::vector<double *> &momenta) = 0;

  // Set the current incoming PDG flavours
  virtual void setInitial(int inid1, int inid2) = 0;

  // Set alpha_s for event-dependent couplings
  virtual void setAlphaS(double in) = 0;

  // Select the fixed Thomson-limit QED coupling
  void setAlphaQEDZero(bool value) noexcept {{ alpha_qed_zero_ = value; }}

 protected:
  // Compute whether the Thomson-limit QED coupling is selected
  bool alphaQEDZero() const noexcept {{ return alpha_qed_zero_; }}

 private:
  bool alpha_qed_zero_ = true;
}};

}}  // namespace {family}

#endif
""".replace(
        "@HELICITY_INCLUDE@",
        output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_Helicity.h"),
    )


# Remove generated trailing whitespace without changing substantive formatting
def clean_cpp_block(text: str) -> str:
    return "\n".join(line.rstrip() for line in text.splitlines())


# Describe one generated matrix evaluation and its external-leg permutation
@dataclass(frozen=True)
class ProcessEvaluation:
    index: int
    permutation: tuple[int, ...]
    matrix_function: str


# Extract one generated integer constant from a subprocess header
def generated_int_constant(header: str, name: str, class_name: str) -> int:
    values = {int(value) for value in re.findall(rf"\b{re.escape(name)}\s*=\s*(\d+)\b", header)}
    if len(values) != 1:
        raise RuntimeError(f"Expected one {name} value for {class_name}")
    return values.pop()


# Validate one generated external-leg permutation
def validate_permutation(permutation: tuple[int, ...], ninitial: int, class_name: str) -> None:
    nexternal = len(permutation)
    if sorted(permutation) != list(range(nexternal)):
        raise RuntimeError(f"Invalid generated permutation in {class_name}: {permutation}")
    if set(permutation[:ninitial]) != set(range(ninitial)):
        raise RuntimeError(f"Generated permutation crosses initial and final legs in {class_name}")


# Validate that sigmaHat channeles exactly the generated subprocess indices
def validate_sigma_hat_indices(text: str, class_name: str, nprocesses: int) -> None:
    block = function_block(text, f"double {class_name}::sigmaHat()")[2]
    raw_indices = re.findall(r"matrix_element\s*\[\s*([^\]]+)\s*\]", cpp_code_only(block))
    if not raw_indices or any(not re.fullmatch(r"\d+", index) for index in raw_indices):
        raise RuntimeError(f"Unsupported sigmaHat matrix index in {class_name}")
    indices = {int(index) for index in raw_indices}
    expected = set(range(nprocesses))
    invalid = indices - expected
    if invalid:
        raise RuntimeError(f"Out-of-range sigmaHat matrix index in {class_name}: {sorted(invalid)}")
    missing = expected - indices
    if missing:
        raise RuntimeError(f"Missing sigmaHat matrix index in {class_name}: {sorted(missing)}")


# Decode one exhaustive generated sigmaKin subprocess evaluation batch
def process_evaluations(text: str, header: str, class_name: str) -> tuple[ProcessEvaluation, ...]:
    ninitial = generated_int_constant(header, "ninitial", class_name)
    nexternal = generated_int_constant(header, "nexternal", class_name)
    nprocesses = generated_int_constant(header, "nprocesses", class_name)
    if ninitial != 2 or nexternal <= ninitial or nprocesses <= 0:
        raise RuntimeError(f"Unsupported generated process dimensions in {class_name}")

    sigma_kin = cpp_code_only(function_block(text, f"void {class_name}::sigmaKin()")[2])
    statement = re.compile(
        r"(?P<assign>\bperm\s*\[\s*(?P<assign_leg>\d+)\s*\]\s*=\s*"
        r"(?P<assign_value>\d+)\s*;)"
        r"|(?P<swap>\b(?:std::)?swap\s*\(\s*perm\s*\[\s*(?P<swap_left>\d+)\s*\]"
        r"\s*,\s*perm\s*\[\s*(?P<swap_right>\d+)\s*\]\s*\)\s*;)"
        r"|(?P<calculate>\bcalculate_wavefunctions\s*\(\s*perm\s*,)"
        r"|(?P<matrix>\bt\s*\[\s*(?P<matrix_index>\d+)\s*\]\s*=\s*"
        r"(?P<matrix_rhs>[^;]+);)"
    )
    permutation = list(range(nexternal))
    wavefunction_permutation = None
    evaluations: dict[int, ProcessEvaluation] = {}
    for match in statement.finditer(sigma_kin):
        if match.group("assign") is not None:
            leg = int(match.group("assign_leg"))
            value = int(match.group("assign_value"))
            if leg >= nexternal or value >= nexternal:
                raise RuntimeError(f"Out-of-range generated permutation index in {class_name}")
            permutation[leg] = value
        elif match.group("swap") is not None:
            left = int(match.group("swap_left"))
            right = int(match.group("swap_right"))
            if left >= nexternal or right >= nexternal:
                raise RuntimeError(f"Out-of-range generated permutation index in {class_name}")
            permutation[left], permutation[right] = permutation[right], permutation[left]
        elif match.group("calculate") is not None:
            wavefunction_permutation = tuple(permutation)
            validate_permutation(wavefunction_permutation, ninitial, class_name)
        else:
            index = int(match.group("matrix_index"))
            if index >= nprocesses:
                raise RuntimeError(f"Out-of-range sigmaKin matrix index in {class_name}: {index}")
            if index in evaluations:
                raise RuntimeError(f"Duplicate sigmaKin matrix index in {class_name}: {index}")
            if wavefunction_permutation is None:
                raise RuntimeError(f"Matrix evaluation without wavefunctions in {class_name}")
            function_match = re.fullmatch(
                r"\s*([A-Za-z_]\w*)\s*\(\s*\)\s*", match.group("matrix_rhs")
            )
            if function_match is None:
                raise RuntimeError(f"Unsupported sigmaKin matrix expression in {class_name}")
            matrix_function = function_match.group(1)
            signature = (
                rf"\bdouble\s+{re.escape(class_name)}::{re.escape(matrix_function)}\s*\(\s*\)"
            )
            if re.search(signature, text) is None:
                raise RuntimeError(f"Matrix function {matrix_function} not found for {class_name}")
            evaluations[index] = ProcessEvaluation(index, wavefunction_permutation, matrix_function)
            if len(evaluations) == nprocesses:
                break

    missing = set(range(nprocesses)) - set(evaluations)
    if missing:
        raise RuntimeError(f"Missing sigmaKin matrix index in {class_name}: {sorted(missing)}")
    validate_sigma_hat_indices(text, class_name, nprocesses)
    return tuple(evaluations[index] for index in range(nprocesses))


# Identify direct or beam-swapped evaluations of the same matrix element
def legacy_evaluation_structure(evaluations: tuple[ProcessEvaluation, ...]) -> bool:
    identity = tuple(range(len(evaluations[0].permutation)))
    if len(evaluations) == 1:
        return evaluations[0].permutation == identity
    mirror = (1, 0, *identity[2:])
    return (
        len(evaluations) == 2
        and evaluations[0].permutation == identity
        and evaluations[1].permutation == mirror
        and evaluations[0].matrix_function == evaluations[1].matrix_function
    )


# Extract generated color-flow data from one subprocess source
def color_data(
    text: str, class_name: str, matrix_function: str | None = None
) -> tuple[int, str, str, str]:
    function_pattern = (
        re.escape(matrix_function) if matrix_function is not None else r"matrix_[^(\s]+"
    )
    match = re.search(
        rf"double\s+{re.escape(class_name)}::{function_pattern}\s*\(\s*\)\s*\{{", text
    )
    if match is None:
        detail = f" {matrix_function}" if matrix_function is not None else ""
        raise RuntimeError(f"Matrix function{detail} not found for {class_name}")
    block = function_block(text, match.group(0).split("{")[0].strip())[2]
    ncolor_match = re.search(r"const int ncolor\s*=\s*(\d+)\s*;", block)
    denom_match = re.search(
        r"const double denom\s*\[\s*ncolor\s*\]\s*=\s*(\{.*?\})\s*;", block, re.S
    )
    factors_match = re.search(
        r"const double cf\s*\[\s*ncolor\s*\]\s*\[\s*ncolor\s*\]\s*=\s*"
        r"(\{\{.*?\}\})\s*;",
        block,
        re.S,
    )
    flows_match = re.search(
        r"// Calculate color flows\s*(.*?)\s*// Sum and square the color flows", block, re.S
    )
    if any(value is None for value in (ncolor_match, denom_match, factors_match, flows_match)):
        raise RuntimeError(f"Incomplete color-flow definition in {match.group(0).split('::')[1]}")
    ncolor = int(ncolor_match.group(1))
    if ncolor <= 0:
        raise RuntimeError(f"Invalid color-flow dimension in {class_name}")
    return ncolor, denom_match.group(1), factors_match.group(1), flows_match.group(1).strip()


# Build a flavour selector equivalent to the generated sigmaHat channel selection
def selected_process_definition(text: str, class_name: str) -> str:
    old = function_block(text, f"double {class_name}::sigmaHat()")[2]
    body = old[old.index("{") :]
    body = re.sub(r"return matrix_element\[(\d+)\](?: \* [0-9.]+)?;", r"return \1;", body)
    body = body.replace("return 0.;", "return -1;")
    return clean_cpp_block(
        f"// Compute the generated subprocess selected by the incoming flavours\n"
        f"int {class_name}::selectedProcess() const\n{body}"
    )


# Build the multiplicity accompanying the selected generated subprocess
def selected_multiplicity_definition(text: str, class_name: str) -> str:
    old = function_block(text, f"double {class_name}::sigmaHat()")[2]
    body = old[old.index("{") :]
    body = re.sub(r"return matrix_element\[\d+\] \* ([0-9.]+);", r"return \1;", body)
    body = re.sub(r"return matrix_element\[\d+\];", "return 1.0;", body)
    body = body.replace("return 0.;", "return 0.0;")
    return clean_cpp_block(
        f"// Compute the identical-flavour multiplicity of the selected subprocess\n"
        f"double {class_name}::selectedProcessMultiplicity() const\n{body}"
    )


# Compute the exact discrete multiplicity selected for each incoming pair
def selected_process_multiplicities(text: str, class_name: str) -> dict[tuple[int, int], int]:
    block = function_block(text, f"double {class_name}::sigmaHat()")[2]
    rows = re.findall(
        r"(?:if|else\s+if)\s*\(\s*id1\s*==\s*(-?\d+)\s*&&\s*"
        r"id2\s*==\s*(-?\d+)\s*\)\s*\{.*?"
        r"return\s+matrix_element\[\d+\](?:\s*\*\s*([1-9]\d*))?\s*;",
        block,
        re.S,
    )
    if len(rows) != block.count("return matrix_element["):
        raise RuntimeError(f"Could not decode every generated multiplicity for {class_name}")
    multiplicities: dict[tuple[int, int], int] = {}
    for first, second, factor in rows:
        initial = (int(first), int(second))
        if initial in multiplicities:
            raise RuntimeError(f"Duplicate generated incoming multiplicity for {class_name}")
        multiplicities[initial] = int(factor) if factor else 1
    return multiplicities


# Verify MG5 sigmaHat grouping against the exact generated channel count
def validate_selected_multiplicities(
    text: str, class_name: str, exact_channels: tuple[Channel, ...]
) -> None:
    expected: dict[tuple[int, int], int] = {}
    for channel in exact_channels:
        expected[channel.initial] = expected.get(channel.initial, 0) + 1
    generated = selected_process_multiplicities(text, class_name)
    if generated != expected:
        raise RuntimeError(
            f"Generated sigmaHat multiplicities and exact channels disagree for {class_name}: "
            f"{generated} != {expected}"
        )


# Indent generated color-flow expressions inside one implementation scope
def indented_color_flows(flows: str, spaces: int) -> str:
    return "\n".join(" " * spaces + line.strip() for line in flows.splitlines())


# Build one selected color-metric projection case
def color_projection_case(
    text: str, class_name: str, evaluation: ProcessEvaluation, leading: bool
) -> str:
    ncolor, denom, factors, flows = color_data(text, class_name, evaluation.matrix_function)
    indented = indented_color_flows(flows, 6)
    if leading:
        result = "return std::vector<std::complex<double>>(jamp, jamp + ncolor);"
    else:
        result = f"""const std::vector<std::complex<double>> color_jamp(jamp, jamp + ncolor);
      const std::vector<double> denominators = {denom};
      const std::vector<std::vector<double>> color_factors = {factors};
      return gra::mg5helas::ColorMetricAmplitudes(color_jamp, denominators, color_factors);"""
    return f"""    case {evaluation.index}: {{
      constexpr int ncolor = {ncolor};
      std::complex<double> jamp[ncolor];
{indented}
      {result}
    }}"""


# Build the selected generated color projection implementation
def selected_color_definition(
    text: str, class_name: str, evaluations: tuple[ProcessEvaluation, ...], leading: bool
) -> str:
    description = (
        "// Compute raw MG5 leading-color flow amplitudes"
        if leading
        else "// Compute orthogonal complex components of the generated MG5 color sum"
    )
    method = "leadingColorAmplitudes" if leading else "colorAmplitudes"
    cases = "\n".join(
        color_projection_case(text, class_name, evaluation, leading) for evaluation in evaluations
    )
    return f"""{description}
std::vector<std::complex<double>> {class_name}::{method}() const
{{
  switch (selectedProcess()) {{
{cases}
    default:
      return {{}};
  }}
}}"""


# Build the generated color-metric projection implementation
def color_definition(
    text: str, class_name: str, evaluations: tuple[ProcessEvaluation, ...] | None = None
) -> str:
    if evaluations is not None:
        return selected_color_definition(text, class_name, evaluations, leading=False)
    ncolor, denom, factors, flows = color_data(text, class_name)
    indented = indented_color_flows(flows, 2)
    return f"""// Compute orthogonal complex components of the generated MG5 color sum
std::vector<std::complex<double>> {class_name}::colorAmplitudes() const
{{
  constexpr int ncolor = {ncolor};
  std::complex<double> jamp[ncolor];
{indented}
  const std::vector<std::complex<double>> color_jamp(jamp, jamp + ncolor);
  const std::vector<double> denominators = {denom};
  const std::vector<std::vector<double>> color_factors = {factors};
  return gra::mg5helas::ColorMetricAmplitudes(color_jamp, denominators, color_factors);
}}"""


# Build the generated leading-color flow amplitudes used for event sampling
def leading_color_definition(
    text: str, class_name: str, evaluations: tuple[ProcessEvaluation, ...] | None = None
) -> str:
    if evaluations is not None:
        return selected_color_definition(text, class_name, evaluations, leading=True)
    ncolor, _, _, flows = color_data(text, class_name)
    indented = indented_color_flows(flows, 2)
    return f"""// Compute raw MG5 leading-color flow amplitudes
std::vector<std::complex<double>> {class_name}::leadingColorAmplitudes() const
{{
  constexpr int ncolor = {ncolor};
  std::complex<double> jamp[ncolor];
{indented}
  return std::vector<std::complex<double>>(jamp, jamp + ncolor);
}}"""


# Extract the common generated helicity table data
def helicity_data(text: str, header: str, class_name: str) -> tuple[int, str, str, int]:
    ncomb_match = re.search(r"const int ncomb\s*=\s*(\d+)\s*;", text)
    helicities_match = re.search(
        r"const int helicities\s*\[\s*ncomb\s*\]\s*\[\s*nexternal\s*\]\s*=\s*"
        r"(\{.*?\})\s*;",
        text,
        re.S,
    )
    denominators_match = re.search(
        r"const int denominators\s*\[\s*nprocesses\s*\]\s*=\s*(\{.*?\})\s*;",
        text,
        re.S,
    )
    if any(value is None for value in (ncomb_match, helicities_match, denominators_match)):
        raise RuntimeError(f"Incomplete helicity definition in {class_name}")
    nprocesses = generated_int_constant(header, "nprocesses", class_name)
    denominator_entries = [
        entry.strip() for entry in denominators_match.group(1)[1:-1].split(",") if entry.strip()
    ]
    if len(denominator_entries) != nprocesses:
        raise RuntimeError(f"Invalid denominator count in {class_name}")
    return (
        int(ncomb_match.group(1)),
        helicities_match.group(1),
        denominators_match.group(1),
        nprocesses,
    )


# Build the existing one-process or initial-mirror helicity implementation
def legacy_helicity_definition(
    ncomb: int,
    helicities: str,
    denominators: str,
    class_name: str,
    alpha_zero: bool,
) -> str:
    alpha_update = "  if (alphaQEDZero()) { pars.setAlphaQEDZero(); }\n" if alpha_zero else ""
    return f"""// Compute fixed-basis complex helicity and orthogonal color components
std::vector<gra::mg5helas::HelicityComponent>
{class_name}::helicityAmplitudes()
{{
  const int selected = selectedProcess();
  if (selected < 0) {{ return {{}}; }}
  pars.setDependentParameters(alphaS);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
{alpha_update.rstrip()}

  const int ncomb = {ncomb};
  const int helicities[ncomb][nexternal] = {helicities};
  const int denominators[nprocesses] = {denominators};
  int perm[nexternal];
  for (int i = 0; i < nexternal; ++i) {{ perm[i] = i; }}
  if (selected == 1) {{ std::swap(perm[0], perm[1]); }}

  // Keep the complete MadGraph denominator until the stable final state is known
  // The subprocess_sum restores its final-state symmetry factor once
  const double average = std::sqrt(selectedProcessMultiplicity() /
                                   static_cast<double>(denominators[selected]));
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (int ihel = 0; ihel < ncomb; ++ihel) {{
    calculate_wavefunctions(perm, helicities[ihel]);
    const auto colors = colorAmplitudes();
    const auto leading_flows = leadingColorAmplitudes();
    if (colors.empty() || leading_flows.empty()) {{ return {{}}; }}
    std::array<int, 2> physical_helicities = {{0, 0}};
    physical_helicities[perm[0]] = helicities[ihel][0];
    physical_helicities[perm[1]] = helicities[ihel][1];
    std::vector<int> outgoing;
    for (int leg = ninitial; leg < nexternal; ++leg) {{
      outgoing.push_back(helicities[ihel][leg]);
    }}
    for (std::size_t color = 0; color < colors.size(); ++color) {{
      gra::mg5helas::HelicityComponent component{{physical_helicities, outgoing, color,
                                                   average * colors[color], {{}}}};
      // Store raw flows once because the orthogonal color components already span this helicity
      if (color == 0) {{
        component.flow_values.reserve(leading_flows.size());
        for (const auto flow : leading_flows) {{
          component.flow_values.push_back(average * flow);
        }}
      }}
      components.push_back(std::move(component));
    }}
  }}
  return components;
}}"""


# Format all selected external-leg permutations as one C++ table
def permutation_table(evaluations: tuple[ProcessEvaluation, ...]) -> str:
    rows = ("{" + ", ".join(str(index) for index in row.permutation) + "}" for row in evaluations)
    return "{" + ", ".join(rows) + "}"


# Build the general grouped-subprocess helicity implementation
def grouped_helicity_definition(
    ncomb: int,
    helicities: str,
    denominators: str,
    class_name: str,
    alpha_zero: bool,
    evaluations: tuple[ProcessEvaluation, ...],
) -> str:
    alpha_update = "  if (alphaQEDZero()) { pars.setAlphaQEDZero(); }\n" if alpha_zero else ""
    permutations = permutation_table(evaluations)
    return f"""// Compute fixed-basis complex helicity and orthogonal color components
std::vector<gra::mg5helas::HelicityComponent>
{class_name}::helicityAmplitudes()
{{
  const int selected = selectedProcess();
  if (selected < 0 || selected >= nprocesses) {{ return {{}}; }}
  pars.setDependentParameters(alphaS);
  pars.setIndependentCouplings();
  pars.setDependentCouplings();
{alpha_update.rstrip()}

  const int ncomb = {ncomb};
  const int helicities[ncomb][nexternal] = {helicities};
  const int denominators[nprocesses] = {denominators};
  const int permutations[nprocesses][nexternal] = {permutations};
  const int *perm = permutations[selected];

  // Keep the complete MadGraph denominator until the stable final state is known
  // The subprocess_sum restores its final-state symmetry factor once
  const double average = std::sqrt(selectedProcessMultiplicity() /
                                   static_cast<double>(denominators[selected]));
  std::vector<gra::mg5helas::HelicityComponent> components;
  for (int ihel = 0; ihel < ncomb; ++ihel) {{
    calculate_wavefunctions(perm, helicities[ihel]);
    const auto colors = colorAmplitudes();
    const auto leading_flows = leadingColorAmplitudes();
    if (colors.empty() || leading_flows.empty()) {{ return {{}}; }}
    std::array<int, 2> physical_helicities = {{0, 0}};
    std::vector<int> outgoing(nexternal - ninitial);
    for (int leg = 0; leg < nexternal; ++leg) {{
      const int physical_leg = perm[leg];
      if (physical_leg < ninitial) {{
        physical_helicities[physical_leg] = helicities[ihel][leg];
      }} else {{
        outgoing[physical_leg - ninitial] = helicities[ihel][leg];
      }}
    }}
    for (std::size_t color = 0; color < colors.size(); ++color) {{
      gra::mg5helas::HelicityComponent component{{physical_helicities, outgoing, color,
                                                   average * colors[color], {{}}}};
      // Store raw flows once because the orthogonal color components already span this helicity
      if (color == 0) {{
        component.flow_values.reserve(leading_flows.size());
        for (const auto flow : leading_flows) {{
          component.flow_values.push_back(average * flow);
        }}
      }}
      components.push_back(std::move(component));
    }}
  }}
  return components;
}}"""


# Build the fixed-basis complex helicity implementation
def helicity_definition(
    text: str,
    header: str,
    class_name: str,
    alpha_zero: bool = True,
    evaluations: tuple[ProcessEvaluation, ...] | None = None,
) -> str:
    ncomb, helicities, denominators, nprocesses = helicity_data(text, header, class_name)
    decoded = evaluations
    if decoded is None and f"void {class_name}::sigmaKin()" in text:
        decoded = process_evaluations(text, header, class_name)
    if decoded is not None:
        if len(decoded) != nprocesses:
            raise RuntimeError(f"Generated process count mismatch in {class_name}")
        if legacy_evaluation_structure(decoded):
            return legacy_helicity_definition(
                ncomb, helicities, denominators, class_name, alpha_zero
            )
        return grouped_helicity_definition(
            ncomb, helicities, denominators, class_name, alpha_zero, decoded
        )
    if nprocesses > 2 or (nprocesses == 2 and "Mirror initial state momenta" not in text):
        raise RuntimeError(f"Unsupported generated permutation structure in {class_name}")
    return legacy_helicity_definition(ncomb, helicities, denominators, class_name, alpha_zero)


# Add declarations and implementations to one generated subprocess
def update_subprocess(process: FamilySubprocess, alpha_zero: bool = True) -> FamilySubprocess:
    class_name = process.class_name
    header = process.header
    source = process.source

    if "helicityAmplitudes() override" not in header:
        header = header.replace(
            "    // Info on the subprocess.",
            "    // Compute fixed-basis complex helicity and color components\n"
            "    std::vector<gra::mg5helas::HelicityComponent> helicityAmplitudes() override;\n\n"
            "    // Info on the subprocess.",
            1,
        )
    if "selectedProcess() const" not in header:
        header = header.replace(
            "    // Private functions to calculate the matrix element for all subprocesses",
            "    // Compute the generated subprocess selected by the incoming flavours\n"
            "    int selectedProcess() const;\n\n"
            "    // Compute the identical-flavour multiplicity of the selected subprocess\n"
            "    double selectedProcessMultiplicity() const;\n\n"
            "    // Compute orthogonal complex components of the generated color sum\n"
            "    std::vector<std::complex<double>> colorAmplitudes() const;\n\n"
            "    // Compute raw MG5 leading-color flow amplitudes\n"
            "    std::vector<std::complex<double>> leadingColorAmplitudes() const;\n\n"
            "    // Private functions to calculate the matrix element for all subprocesses",
            1,
        )
    if "leadingColorAmplitudes() const" not in header:
        header = header.replace(
            "    // Private functions to calculate the matrix element for all subprocesses",
            "    // Compute raw MG5 leading-color flow amplitudes\n"
            "    std::vector<std::complex<double>> leadingColorAmplitudes() const;\n\n"
            "    // Private functions to calculate the matrix element for all subprocesses",
            1,
        )

    for comment in (
        "// Compute the generated subprocess selected by the incoming flavours\n",
        "// Compute the identical-flavour multiplicity of the selected subprocess\n",
        "// Compute orthogonal complex components of the generated MG5 color sum\n",
        "// Compute raw MG5 leading-color flow amplitudes\n",
        "// Compute fixed-basis complex helicity and orthogonal color components\n",
    ):
        source = source.replace(comment, "")

    for signature in (
        f"int {class_name}::selectedProcess() const",
        f"double {class_name}::selectedProcessMultiplicity() const",
        f"std::vector<std::complex<double>> {class_name}::colorAmplitudes() const",
        f"std::vector<std::complex<double>> {class_name}::leadingColorAmplitudes() const",
        f"std::vector<gra::mg5helas::HelicityComponent>\n{class_name}::helicityAmplitudes()",
    ):
        if signature in source:
            source = source.replace(function_block(source, signature)[2], "", 1)

    evaluations = process_evaluations(source, header, class_name)
    selected_evaluations = None if legacy_evaluation_structure(evaluations) else evaluations
    generated = "\n\n".join(
        (
            selected_process_definition(source, class_name),
            selected_multiplicity_definition(source, class_name),
            color_definition(source, class_name, selected_evaluations),
            leading_color_definition(source, class_name, selected_evaluations),
            helicity_definition(source, header, class_name, alpha_zero, evaluations),
        )
    )
    marker = "//==========================================================================\n// Private class member functions"
    source = source.replace(marker, generated + "\n\n" + marker, 1)
    # Keep repeated regeneration text-idempotent after replacing old generated blocks
    source = re.sub(
        r"\n{3,}(?=// Compute the generated subprocess selected by the incoming flavours)",
        "\n\n",
        source,
        count=1,
    )
    return FamilySubprocess(
        class_name=class_name,
        header=header,
        source=source,
        channels=process.channels,
    )


# Format one exact channel topology node in C++
def format_topology_node(node: ChannelTopologyNode) -> str:
    daughters = ", ".join(format_topology_node(daughter) for daughter in node.daughters)
    return f"AmplitudeTopologyNode{{{{{node.pdg}}}, {{{daughters}}}}}"


# Format one generated subprocess channel in C++
def format_channel(channel: Channel) -> str:
    if not channel.external_color_representations or not channel.external_color_flows:
        raise RuntimeError("Cannot format a subprocess channel without external color flows")
    final_values = ", ".join(str(value) for value in channel.final)
    topology = ", ".join(format_topology_node(node) for node in channel.topology)
    representations = ", ".join(str(value) for value in channel.external_color_representations)
    color_flows = ", ".join(
        "{" + ", ".join(f"ColorFlowLeg{{{color}, {anticolor}}}" for color, anticolor in flow) + "}"
        for flow in channel.external_color_flows
    )
    return (
        f"Channel{{{{{channel.initial[0]}, {channel.initial[1]}}}, "
        f"{{{final_values}}}, {{{topology}}}, {{{representations}}}, "
        f"{{{color_flows}}}}}"
    )


# Append exact stable leaves from one channel topology node
def append_channel_stable_pdgs(node: ChannelTopologyNode, stable_pdgs: list[int]) -> None:
    if node.pdg == 0 or abs(node.pdg) == 89:
        raise RuntimeError("Generated channel topology contains an unresolved particle")
    if not node.daughters:
        stable_pdgs.append(node.pdg)
        return
    for daughter in node.daughters:
        append_channel_stable_pdgs(daughter, stable_pdgs)


# Validate exact channel topology across all generated subprocess channels
def validate_channel_structure(channel_rows: list[tuple[Channel, ...]]) -> None:
    identities: set[
        tuple[
            tuple[int, int],
            tuple[int, ...],
            tuple[ChannelTopologyNode, ...],
        ]
    ] = set()
    for channels in channel_rows:
        if not channels:
            raise RuntimeError("Generated subprocess has no exact channels")
        for channel in channels:
            if (
                0 in channel.initial
                or any(abs(pdg) == 89 for pdg in channel.initial)
                or not channel.final
                or not channel.topology
            ):
                raise RuntimeError("Generated subprocess channel is not exact")
            topology_final: list[int] = []
            for node in channel.topology:
                append_channel_stable_pdgs(node, topology_final)
            if tuple(topology_final) != channel.final:
                raise RuntimeError("Generated channel topology and stable final state disagree")
            identity = (channel.initial, channel.final, channel.topology)
            if identity in identities:
                raise RuntimeError("Duplicate exact generated subprocess channel")
            identities.add(identity)


# Reject one physical color slot that disagrees with its signed representation
def validate_external_color_leg(representation: int, leg: tuple[int, int]) -> None:
    color, anticolor = leg
    if representation == 1 and (color != 0 or anticolor != 0):
        raise RuntimeError("Generated singlet leg has external color tags")
    if representation == 3 and (color <= 0 or anticolor != 0):
        raise RuntimeError("Generated triplet leg has invalid external color tags")
    if representation == -3 and (color != 0 or anticolor <= 0):
        raise RuntimeError("Generated antitriplet leg has invalid external color tags")
    if representation == 8 and (color <= 0 or anticolor <= 0 or color == anticolor):
        raise RuntimeError("Generated octet leg has invalid external color tags")
    if representation == 6 and not (color > 0 and anticolor < 0):
        raise RuntimeError("Generated sextet leg has invalid external color tags")
    if representation == -6 and not (color < 0 and anticolor > 0):
        raise RuntimeError("Generated antisextet leg has invalid external color tags")


# Validate signed representations and crossed color-line balance for one channel
def validate_external_color_structure(channel: Channel) -> None:
    external_count = 2 + len(channel.final)
    if len(channel.external_color_representations) != external_count:
        raise RuntimeError("Generated external color representation count disagrees")
    if any(
        representation not in (1, 3, -3, 6, -6, 8)
        for representation in channel.external_color_representations
    ):
        raise RuntimeError("Generated external color representation is unsupported")
    if not channel.external_color_flows:
        raise RuntimeError("Generated subprocess channel has no external color flows")
    for flow in channel.external_color_flows:
        if len(flow) != external_count:
            raise RuntimeError("Generated external color flow leg count disagrees")
        for representation, leg in zip(channel.external_color_representations, flow, strict=True):
            validate_external_color_leg(representation, leg)
        crossed = list(flow)
        crossed[0:2] = [
            (crossed[0][1], crossed[0][0]),
            (crossed[1][1], crossed[1][0]),
        ]
        colors = [color for color, _ in crossed if color > 0] + [-anti for _, anti in crossed if anti < 0]
        anticolors = [anti for _, anti in crossed if anti > 0] + [-color for color, _ in crossed if color < 0]
        if any(tag >= 500 for tag in (*colors, *anticolors)):
            raise RuntimeError("Generated external color tag is outside local MG range")
        if sorted(colors) != sorted(anticolors) or len(colors) != len(set(colors)):
            raise RuntimeError("Generated external color flow is not line balanced")


# Validate complete exact channels across generated subprocesses
def validate_channels(channel_rows: list[tuple[Channel, ...]]) -> None:
    validate_channel_structure(channel_rows)
    for channels in channel_rows:
        for channel in channels:
            validate_external_color_structure(channel)


# Serialize one exact topology node into the family channel schema
def topology_node_data(node: ChannelTopologyNode) -> dict[str, Any]:
    return {
        "pdg": node.pdg,
        "daughters": [topology_node_data(daughter) for daughter in node.daughters],
    }


# Decode one exact topology node from the family channel schema
def topology_node_from_data(data: Any) -> ChannelTopologyNode:
    if not isinstance(data, dict) or set(data) != {"pdg", "daughters"}:
        raise RuntimeError("Invalid generated channel topology data")
    if type(data["pdg"]) is not int or not isinstance(data["daughters"], list):
        raise RuntimeError("Invalid generated channel topology values")
    return ChannelTopologyNode(
        data["pdg"],
        tuple(topology_node_from_data(child) for child in data["daughters"]),
    )


# Decode one strict integer array from persistent family channels
def integer_tuple(value: Any, name: str) -> tuple[int, ...]:
    if not isinstance(value, list) or any(type(item) is not int for item in value):
        raise RuntimeError(f"Invalid {name} in generated family channels")
    return tuple(value)


# Serialize one complete exact subprocess channel
def channel_to_data(channel: Channel) -> dict[str, Any]:
    return {
        "initial": list(channel.initial),
        "final": list(channel.final),
        "topology": [topology_node_data(node) for node in channel.topology],
        "external_color_representations": list(channel.external_color_representations),
        "external_color_flows": [
            [value for leg in flow for value in leg] for flow in channel.external_color_flows
        ],
    }


# Decode one complete exact subprocess channel
def channel_from_data(data: Any) -> Channel:
    required = {
        "initial",
        "final",
        "topology",
        "external_color_representations",
        "external_color_flows",
    }
    if not isinstance(data, dict) or set(data) != required:
        raise RuntimeError("Invalid exact channel in generated family channels")
    initial = integer_tuple(data["initial"], "initial state")
    if len(initial) != 2:
        raise RuntimeError("Generated family channels require two incoming particles")
    if not isinstance(data["topology"], list):
        raise RuntimeError("Invalid topology in generated family channels")
    if not isinstance(data["external_color_flows"], list):
        raise RuntimeError("Invalid color flows in generated family channels")
    flow_rows = []
    for row in data["external_color_flows"]:
        flat = integer_tuple(row, "external color flow")
        if len(flat) % 2 != 0:
            raise RuntimeError("Generated external color flow has an odd value count")
        flow_rows.append(tuple(zip(flat[0::2], flat[1::2], strict=True)))
    return Channel(
        (initial[0], initial[1]),
        integer_tuple(data["final"], "final state"),
        tuple(topology_node_from_data(node) for node in data["topology"]),
        integer_tuple(
            data["external_color_representations"],
            "external color representations",
        ),
        tuple(flow_rows),
    )


# Serialize complete deterministic family channels
def family_channels_data(
    family_channels: FamilyChannels,
) -> dict[str, Any]:
    return {
        "schema": 1,
        "family": family_channels.family,
        "generation_setup": family_generation_setup_data(family_channels.generation_setup),
        "mg5_source": mg5_source_data(family_channels.mg5_source),
        "subprocesses": [
            {
                "generated_name": subprocess.generated_name,
                "channels": [channel_to_data(channel) for channel in subprocess.channels],
            }
            for subprocess in family_channels.subprocesses
        ],
    }


# Decode and validate complete persistent family channels
def family_channels_from_data(data: Any) -> FamilyChannels:
    if (
        not isinstance(data, dict)
        or set(data) != {"schema", "family", "generation_setup", "mg5_source", "subprocesses"}
        or type(data["schema"]) is not int
        or data["schema"] != 1
        or not isinstance(data["family"], str)
        or not data["family"]
        or not isinstance(data["subprocesses"], list)
    ):
        raise RuntimeError("Invalid generated family channel schema")
    subprocesses = []
    for subprocess in data["subprocesses"]:
        if (
            not isinstance(subprocess, dict)
            or set(subprocess) != {"generated_name", "channels"}
            or not isinstance(subprocess["generated_name"], str)
            or not subprocess["generated_name"]
            or not isinstance(subprocess["channels"], list)
        ):
            raise RuntimeError("Invalid generated family subprocess channels")
        subprocesses.append(
            SubprocessChannels(
                subprocess["generated_name"],
                tuple(channel_from_data(channel) for channel in subprocess["channels"]),
            )
        )
    names = [subprocess.generated_name for subprocess in subprocesses]
    if len(names) != len(set(names)):
        raise RuntimeError("Duplicate generated family subprocess channels")
    family_channels = FamilyChannels(
        data["family"],
        family_generation_setup_from_data(
            data["generation_setup"], "Generated family channel generation setup"
        ),
        mg5_source_from_data(data["mg5_source"], "Generated family channel source mg5_source"),
        tuple(subprocesses),
    )
    validate_channels([subprocess.channels for subprocess in family_channels.subprocesses])
    return family_channels


# Build complete persistent channels in transformed subprocess order
def make_family_channels(
    family: str,
    generation_setup: FamilyGenerationSetup,
    mg5_source: MG5Source,
    subprocesses: list[tuple[str, FamilySubprocess]],
) -> FamilyChannels:
    family_channels = FamilyChannels(
        family,
        generation_setup,
        mg5_source,
        tuple(
            SubprocessChannels(generated_name, subprocess.channels)
            for generated_name, subprocess in subprocesses
        ),
    )
    validate_channels([subprocess.channels for subprocess in family_channels.subprocesses])
    return family_channels


# Generate the converter-owned C++ subprocess channels
def subprocess_channels_block(channel_rows: list[tuple[Channel, ...]]) -> str:
    validate_channels(channel_rows)
    rows = ",\n".join(
        "      {" + ", ".join(format_channel(channel) for channel in channels) + "}"
        for channels in channel_rows
    )
    return f"""inline std::vector<std::vector<gra::mg5::Channel>> SubprocessChannels() {{
  using gra::AmplitudeTopologyNode;
  using gra::mg5helas::ColorFlowLeg;
  using gra::mg5::Channel;
  return {{
{rows}
  }};
}}"""


# Generate subprocess construction and exact channels
def processes_header(
    family: str,
    subprocesses: list[FamilySubprocess],
    include_directory: str,
) -> str:
    guard = f"GRANIITTI_AMPLITUDE_{family}_PROCESSES_H"
    includes = "\n".join(
        f'#include "{output_layout.include_path(include_directory, f"{process.class_name}.h")}"'
        for process in subprocesses
    )
    builders = "\n".join(
        (
            "  {\n"
            f"    auto proc = std::make_unique<{process.class_name}>();\n"
            "    proc->initProc(param_card);\n"
            "    out.push_back(std::move(proc));\n"
            "  }"
        )
        for process in subprocesses
    )
    subprocess_sum_include = f'#include "{output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_SubprocessSum.h")}"\n'
    process_base_include = output_layout.include_path(include_directory, "ProcessBase.h")
    channels_block = subprocess_channels_block([process.channels for process in subprocesses])
    channel_code = f"""

// Compute generated subprocess channels in BuildProcesses order
{channels_block}
"""

    channel_code += f"""
// Build generated subprocesses with their exact channels
inline void BuildSubprocesses(std::vector<gra::mg5::Subprocess<ProcessBase>> &out,
                              const std::string &param_card) {{
  std::vector<std::unique_ptr<ProcessBase>> matrix_elements;
  BuildProcesses(matrix_elements, param_card);
  auto channels = SubprocessChannels();
  if (matrix_elements.size() != channels.size()) {{
    throw std::invalid_argument(
        "{family}::BuildSubprocesses: subprocess and channel counts disagree");
  }}
  out.reserve(out.size() + matrix_elements.size());
  for (std::size_t i = 0; i < matrix_elements.size(); ++i) {{
    out.push_back({{std::move(matrix_elements[i]), std::move(channels[i])}});
  }}
}}
"""
    return f"""#ifndef {guard}
#define {guard}

#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

{subprocess_sum_include}#include "{process_base_include}"
{includes}

namespace {family} {{

// Build all generated MG5 subprocesses for this process family
inline void BuildProcesses(std::vector<std::unique_ptr<ProcessBase>> &out,
                           const std::string &param_card) {{
{builders}
}}
{channel_code}
}}  // namespace {family}

#endif
"""
