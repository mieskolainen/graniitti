#!/usr/bin/env python3
#
# Generate the common MG5 amplitude process registry
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import hashlib
import json
import math
import re
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from core.io.serialize import load_json_file

from . import output_layout

STANDALONE_PROCESS_FAMILIES = frozenset({"DURHAM", "PHOTON"})
PROCESS_FAMILY_PREFIX = "MG5_"
DURHAM_PROJECTION = output_layout.DURHAM_PROJECTION
PHOTON_PROJECTION = output_layout.PHOTON_PROJECTION
PARTON_PROJECTION = output_layout.PARTON_PROJECTION
PROJECTIONS = frozenset({DURHAM_PROJECTION, PHOTON_PROJECTION, PARTON_PROJECTION})
STANDALONE_PROJECTIONS = frozenset({DURHAM_PROJECTION, PHOTON_PROJECTION})
FAMILY_PROJECTIONS = frozenset({PHOTON_PROJECTION, PARTON_PROJECTION})
FAMILY_CHANNELS_SCHEMA = 1


@dataclass(frozen=True)
class ProcessNode:
    """One occurrence in an ordered MG5 production and decay tree"""

    particle: str
    daughters: tuple[ProcessNode, ...] = ()


@dataclass(frozen=True)
class ParsedProcessSyntax:
    """Ordered occurrence-aware production tree from one MG5 process"""

    incoming: tuple[str, ...]
    production: tuple[ProcessNode, ...]

    @property
    def has_decay_chain(self) -> bool:
        """Return whether at least one produced occurrence has daughters"""

        return any(node.daughters for node in self.production)


# Store the model QED replacement which changes transformed family parameters
@dataclass(frozen=True)
class AlphaQEDGenerationSetup:
    """Normalized QED replacement inputs for one generated family"""

    charge: str
    charge_square: str
    inverse_alpha_zero: float


# Store every manifest input which changes one generated family
@dataclass(frozen=True)
class FamilyGenerationSetup:
    """Normalized MadGraph generation and family transformation inputs"""

    model: str
    model_import: str
    complex_mass_scheme: bool
    definitions: tuple[str, ...]
    processes: tuple[str, ...]
    particle_pdgs: tuple[tuple[str, int], ...]
    alpha_qed: AlphaQEDGenerationSetup | None
    mass_overrides: tuple[tuple[int, float], ...]


# Store the MG5 release and installed generated family source digest
@dataclass(frozen=True)
class MG5Source:
    """MG5 version and source digest for one generated family"""

    mg5_version: str
    source_digest: str


# Decode one strict ordered string array from a generation setup
def generation_string_tuple(value: Any, context: str, *, allow_empty: bool) -> tuple[str, ...]:
    if (
        not isinstance(value, list)
        or (not allow_empty and not value)
        or any(not isinstance(item, str) or not item for item in value)
        or len(value) != len(set(value))
    ):
        raise RuntimeError(f"{context} is invalid")
    return tuple(value)


# Decode one optional strict QED replacement setup
def alpha_qed_generation_setup_from_data(data: Any, context: str) -> AlphaQEDGenerationSetup | None:
    if data is None:
        return None
    required = {"charge", "charge_square", "inverse_alpha_zero"}
    if not isinstance(data, dict) or set(data) != required:
        raise RuntimeError(f"{context} alpha_qed is invalid")
    charge = data["charge"]
    charge_square = data["charge_square"]
    inverse = data["inverse_alpha_zero"]
    if (
        not isinstance(charge, str)
        or not charge
        or not isinstance(charge_square, str)
        or not charge_square
        or not isinstance(inverse, (int, float))
        or isinstance(inverse, bool)
        or not math.isfinite(inverse)
        or inverse <= 0.0
    ):
        raise RuntimeError(f"{context} alpha_qed is invalid")
    return AlphaQEDGenerationSetup(charge, charge_square, float(inverse))


# Decode one complete strict family generation setup
def family_generation_setup_from_data(
    data: Any, context: str = "Family generation setup"
) -> FamilyGenerationSetup:
    required = {
        "model",
        "model_import",
        "complex_mass_scheme",
        "definitions",
        "processes",
        "particle_pdgs",
        "alpha_qed",
        "mass_overrides",
    }
    if not isinstance(data, dict) or set(data) != required:
        raise RuntimeError(f"{context} is invalid")
    model = data["model"]
    model_import = data["model_import"]
    complex_mass_scheme = data["complex_mass_scheme"]
    particle_pdgs = data["particle_pdgs"]
    mass_overrides = data["mass_overrides"]
    if (
        not isinstance(model, str)
        or not model
        or not isinstance(model_import, str)
        or not model_import
        or not isinstance(complex_mass_scheme, bool)
        or not isinstance(particle_pdgs, dict)
        or any(
            not isinstance(particle, str) or not particle or type(pdg) is not int or pdg == 0
            for particle, pdg in particle_pdgs.items()
        )
        or not isinstance(mass_overrides, dict)
        or any(
            not isinstance(pdg, str)
            or re.fullmatch(r"-?[1-9][0-9]*", pdg) is None
            or not isinstance(mass, (int, float))
            or isinstance(mass, bool)
            or not math.isfinite(mass)
            or mass < 0.0
            for pdg, mass in mass_overrides.items()
        )
    ):
        raise RuntimeError(f"{context} is invalid")
    return FamilyGenerationSetup(
        model=model,
        model_import=model_import,
        complex_mass_scheme=complex_mass_scheme,
        definitions=generation_string_tuple(
            data["definitions"], f"{context} definitions", allow_empty=True
        ),
        processes=generation_string_tuple(
            data["processes"], f"{context} processes", allow_empty=False
        ),
        particle_pdgs=tuple(sorted(particle_pdgs.items())),
        alpha_qed=alpha_qed_generation_setup_from_data(data["alpha_qed"], context),
        mass_overrides=tuple(
            sorted((int(pdg), float(mass)) for pdg, mass in mass_overrides.items())
        ),
    )


# Normalize generation inputs from one validated family and model
def family_generation_setup_from_manifest(
    family: dict[str, Any], model: dict[str, Any]
) -> FamilyGenerationSetup:
    alpha = model.get("alpha_qed")
    if isinstance(alpha, dict):
        alpha = {
            "charge": alpha.get("charge"),
            "charge_square": alpha.get("charge_square"),
            "inverse_alpha_zero": alpha.get("inverse_alpha_zero", 137.03599908),
        }
    return family_generation_setup_from_data(
        {
            "model": family["model"],
            "model_import": model.get("import"),
            "complex_mass_scheme": model.get("complex_mass_scheme", False),
            "definitions": family.get("definitions"),
            "processes": family.get("processes"),
            "particle_pdgs": model.get("particle_pdgs", {}),
            "alpha_qed": alpha,
            "mass_overrides": family.get("mass_overrides", {}),
        },
        f"Family {family.get('name', '<unknown>')} generation setup",
    )


# Serialize one normalized family generation setup
def family_generation_setup_data(
    setup: FamilyGenerationSetup,
) -> dict[str, Any]:
    alpha_qed = None
    if setup.alpha_qed is not None:
        alpha_qed = {
            "charge": setup.alpha_qed.charge,
            "charge_square": setup.alpha_qed.charge_square,
            "inverse_alpha_zero": setup.alpha_qed.inverse_alpha_zero,
        }
    return {
        "model": setup.model,
        "model_import": setup.model_import,
        "complex_mass_scheme": setup.complex_mass_scheme,
        "definitions": list(setup.definitions),
        "processes": list(setup.processes),
        "particle_pdgs": dict(setup.particle_pdgs),
        "alpha_qed": alpha_qed,
        "mass_overrides": {str(pdg): mass for pdg, mass in setup.mass_overrides},
    }


# Decode one strict generated MG5 source
def mg5_source_from_data(data: Any, context: str = "Family MG5 source") -> MG5Source:
    if not isinstance(data, dict) or set(data) != {
        "mg5_version",
        "source_digest",
    }:
        raise RuntimeError(f"{context} is invalid")
    version = data["mg5_version"]
    digest = data["source_digest"]
    if (
        not isinstance(version, str)
        or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9.+_-]*", version) is None
        or not isinstance(digest, str)
        or re.fullmatch(r"[0-9a-f]{64}", digest) is None
    ):
        raise RuntimeError(f"{context} is invalid")
    return MG5Source(version, digest)


# Serialize one generated MG5 source
def mg5_source_data(
    source: MG5Source,
) -> dict[str, str]:
    return {
        "mg5_version": source.mg5_version,
        "source_digest": source.source_digest,
    }


# Compute true for one MG5 coupling-order or diagram-selection token
def is_mg5_modifier(token: str) -> bool:
    return (
        "=" in token
        or token.startswith(("/", "$", "@"))
        or token in {"WEIGHTED", "WEIGHTED<=", "WEIGHTED>="}
    )


# Extract particle tokens before any MG5 coupling-order modifiers
def particle_tokens(text: str) -> tuple[str, ...]:
    particles: list[str] = []
    for token in text.strip().split():
        if is_mg5_modifier(token):
            break
        particle = token[1:-1] if token.startswith("(") and token.endswith(")") else token
        particles.append(particle)
    if not particles or any(not particle for particle in particles):
        raise RuntimeError(f"Invalid MG5 particle list: {text}")
    return tuple(particles)


# Split comma-separated MG5 clauses without cutting nested decay groups
def split_top_level_clauses(text: str) -> tuple[str, ...]:
    clauses: list[str] = []
    start = 0
    depth = 0
    for index, character in enumerate(text):
        if character == "(":
            depth += 1
        elif character == ")":
            depth -= 1
            if depth < 0:
                raise RuntimeError(f"Unmatched ')' in MG5 process syntax: {text}")
        elif character == "," and depth == 0:
            clauses.append(text[start:index].strip())
            start = index + 1
    if depth != 0:
        raise RuntimeError(f"Unmatched '(' in MG5 process syntax: {text}")
    clauses.append(text[start:].strip())
    if any(not clause for clause in clauses):
        raise RuntimeError(f"Invalid MG5 process syntax: {text}")
    return tuple(clauses)


# Remove balanced parentheses enclosing one complete decay group
def strip_outer_group(text: str) -> str:
    group = text.strip()
    while group.startswith("("):
        depth = 0
        closing = -1
        for index, character in enumerate(group):
            if character == "(":
                depth += 1
            elif character == ")":
                depth -= 1
                if depth == 0:
                    closing = index
                    break
        if closing != len(group) - 1:
            raise RuntimeError(f"Malformed parenthesized MG5 decay group: {text}")
        group = group[1:-1].strip()
        if not group:
            raise RuntimeError(f"Empty parenthesized MG5 decay group: {text}")
    return group


# Parse one production or decay clause into incoming and outgoing particles
def parse_process_clause(clause: str) -> tuple[tuple[str, ...], tuple[str, ...]]:
    sides = re.split(r"\s+>\s+", clause, maxsplit=1)
    if len(sides) != 2:
        raise RuntimeError(f"Unsupported MG5 process clause: {clause}")
    incoming, outgoing = sides
    return particle_tokens(incoming), particle_tokens(outgoing)


# Find the leftmost matching unexpanded occurrence at the shallowest depth
def shallowest_decay_path(nodes: tuple[ProcessNode, ...], particle: str) -> tuple[int, ...] | None:
    level = [((index,), node) for index, node in enumerate(nodes)]
    while level:
        for path, node in level:
            if node.particle == particle and not node.daughters:
                return path
        level = [
            ((*path, index), daughter)
            for path, node in level
            for index, daughter in enumerate(node.daughters)
        ]
    return None


# Replace one occurrence selected by its ordered tree path
def replace_decay_path(
    nodes: tuple[ProcessNode, ...], path: tuple[int, ...], decay: ProcessNode
) -> tuple[ProcessNode, ...]:
    updated = list(nodes)
    index = path[0]
    if len(path) == 1:
        updated[index] = decay
    else:
        node = updated[index]
        updated[index] = ProcessNode(
            node.particle, replace_decay_path(node.daughters, path[1:], decay)
        )
    return tuple(updated)


# Attach one decay to the shallowest matching unexpanded occurrence
def attach_decay(
    nodes: tuple[ProcessNode, ...], decay: ProcessNode
) -> tuple[tuple[ProcessNode, ...], bool]:
    path = shallowest_decay_path(nodes, decay.particle)
    if path is None:
        return nodes, False
    return replace_decay_path(nodes, path, decay), True


# Parse one optionally parenthesized decay chain into a scoped branch
def parse_decay_group(group: str, process_syntax: str) -> ProcessNode:
    clauses = split_top_level_clauses(strip_outer_group(group))
    incoming, outgoing = parse_process_clause(clauses[0])
    if len(incoming) != 1:
        raise RuntimeError(f"MG5 decay clause has multiple mothers: {process_syntax}")
    root = ProcessNode(incoming[0], tuple(ProcessNode(particle) for particle in outgoing))
    for clause in clauses[1:]:
        nested = parse_decay_group(clause, process_syntax)
        attached, found = attach_decay((root,), nested)
        if not found:
            raise RuntimeError(
                f"MG5 decay mother {nested.particle} is absent from its scoped branch: "
                f"{process_syntax}"
            )
        root = attached[0]
    return root


# Parse one MG5 production process and its optional decay-chain groups
def parse_process_syntax(process_syntax: str) -> ParsedProcessSyntax:
    clauses = split_top_level_clauses(process_syntax.strip())
    incoming, outgoing = parse_process_clause(clauses[0])
    if not incoming:
        raise RuntimeError(f"Invalid MG5 incoming state: {process_syntax}")
    production = tuple(ProcessNode(particle) for particle in outgoing)
    for clause in clauses[1:]:
        decay = parse_decay_group(clause, process_syntax)
        production, attached = attach_decay(production, decay)
        if not attached:
            raise RuntimeError(
                f"MG5 decay mother {decay.particle} is absent or already assigned: {process_syntax}"
            )
    return ParsedProcessSyntax(incoming=incoming, production=production)


# Translate one MG5 particle token to the GRANIITTI process spelling
def graniitti_particle(token: str) -> str:
    if token == "a":
        return "gamma"
    if token == "ta-":
        return "tau-"
    if token == "ta+":
        return "tau+"
    electroweak = re.fullmatch(r"([wzh])([+-]?)", token)
    if electroweak is not None:
        return electroweak.group(1).upper() + electroweak.group(2)
    return token


# Render one nested GRANIITTI decay branch
def render_branch(node: ProcessNode) -> str:
    particle = graniitti_particle(node.particle)
    if not node.daughters:
        return particle
    daughters = " ".join(render_branch(daughter) for daughter in node.daughters)
    return f"{particle} > {{{daughters}}}"


# Render valid nested GRANIITTI final-state syntax
def final_state_syntax(parsed: ParsedProcessSyntax) -> str:
    return " ".join(render_branch(node) for node in parsed.production)


# Collect stable MG5 particle tokens from one occurrence in momentum order
def stable_branch_tokens(node: ProcessNode) -> tuple[str, ...]:
    if not node.daughters:
        return (node.particle,)
    leaves: list[str] = []
    for daughter in node.daughters:
        leaves.extend(stable_branch_tokens(daughter))
    return tuple(leaves)


# Flatten one MG5 production and decay tree to ordered stable tokens
def stable_tokens(parsed: ParsedProcessSyntax) -> tuple[str, ...]:
    leaves: list[str] = []
    for node in parsed.production:
        leaves.extend(stable_branch_tokens(node))
    return tuple(leaves)


# Format a stable final state and compact a pure multi-gluon state
def stable_final_state(parsed: ParsedProcessSyntax) -> str:
    leaves = tuple(graniitti_particle(token) for token in stable_tokens(parsed))
    if leaves and all(token == "g" for token in leaves):
        return "g" * len(leaves)
    return " ".join(leaves)


# Render one branch of the deterministic occurrence-aware topology signature
def topology_branch(node: ProcessNode) -> str:
    if not node.daughters:
        return node.particle
    daughters = ",".join(topology_branch(daughter) for daughter in node.daughters)
    return f"{node.particle}({daughters})"


# Compute the ordered production and decay-tree signature
def topology_signature(parsed: ParsedProcessSyntax) -> str:
    return ",".join(topology_branch(node) for node in parsed.production)


# Compute canonical Standard Model particle selectors used by MG5 syntax
def standard_particle_selectors() -> dict[str, tuple[int, ...]]:
    return {
        "d": (1,),
        "d~": (-1,),
        "u": (2,),
        "u~": (-2,),
        "s": (3,),
        "s~": (-3,),
        "c": (4,),
        "c~": (-4,),
        "b": (5,),
        "b~": (-5,),
        "t": (6,),
        "t~": (-6,),
        "e-": (11,),
        "e+": (-11,),
        "ve": (12,),
        "ve~": (-12,),
        "mu-": (13,),
        "mu+": (-13,),
        "vm": (14,),
        "vm~": (-14,),
        "ta-": (15,),
        "ta+": (-15,),
        "vt": (16,),
        "vt~": (-16,),
        "g": (21,),
        "a": (22,),
        "z": (23,),
        "w+": (24,),
        "w-": (-24,),
        "h": (25,),
        "j": (-89, -5, -4, -3, -2, -1, 1, 2, 3, 4, 5, 21, 89),
    }


# Normalize one model-supplied exact particle PDG selector
def normalize_particle_pdg(value: Any, context: str) -> tuple[int, ...]:
    if isinstance(value, bool) or not isinstance(value, int) or value == 0:
        raise RuntimeError(f"{context} must be one nonzero integer PDG code")
    return (value,)


# Merge optional model-specific particle PDG overrides with the SM map
def model_particle_selectors(
    manifest: dict[str, Any], model_name: str
) -> dict[str, tuple[int, ...]]:
    selectors = standard_particle_selectors()
    model = manifest.get("models", {}).get(model_name, {})
    overrides = model.get("particle_pdgs", {})
    if not isinstance(overrides, dict):
        raise RuntimeError(f"Model {model_name} particle_pdgs must be an object")
    for token, value in overrides.items():
        if not isinstance(token, str) or not token:
            raise RuntimeError(f"Model {model_name} has an invalid particle token")
        selectors[token] = normalize_particle_pdg(value, f"Model {model_name} particle {token}")
    return selectors


# Parse family-level MG5 define aliases without resolving their dependencies
def parse_family_aliases(definitions: list[str]) -> dict[str, tuple[str, ...]]:
    aliases: dict[str, tuple[str, ...]] = {}
    for definition in definitions:
        match = re.fullmatch(r"\s*define\s+(\S+)\s*=\s*(.*?)\s*", definition)
        if match is None or not match.group(2):
            raise RuntimeError(f"Unsupported MG5 family definition: {definition}")
        alias = match.group(1)
        particles = tuple(match.group(2).split())
        if alias in aliases:
            raise RuntimeError(f"Duplicate MG5 family alias: {alias}")
        aliases[alias] = particles
    return aliases


# Resolve one particle or recursively defined alias to an explicit PDG set
def resolve_particle_selector(
    token: str,
    selectors: dict[str, tuple[int, ...]],
    aliases: dict[str, tuple[str, ...]],
    resolving: tuple[str, ...] = (),
) -> tuple[int, ...]:
    if token in resolving:
        chain = " -> ".join((*resolving, token))
        raise RuntimeError(f"Cyclic MG5 particle alias: {chain}")
    if token in aliases:
        allowed: set[int] = set()
        for particle in aliases[token]:
            allowed.update(
                resolve_particle_selector(particle, selectors, aliases, (*resolving, token))
            )
    elif token in selectors:
        allowed = set(selectors[token])
    else:
        raise RuntimeError(f"Unknown MG5 particle token {token}; add a model particle_pdgs entry")
    if token == "j":
        allowed.update((-89, 89))
    return tuple(sorted(allowed))


# Convert one parsed branch to an ordered numeric topology node
def topology_node(
    node: ProcessNode,
    selectors: dict[str, tuple[int, ...]],
    aliases: dict[str, tuple[str, ...]],
    stable_leaf_pdgs: list[int] | None,
    leaf_index: int,
) -> tuple[dict[str, Any], int]:
    if not node.daughters and stable_leaf_pdgs is not None:
        if leaf_index >= len(stable_leaf_pdgs):
            raise RuntimeError("Stable final-state PDG list is shorter than the topology")
        stable_pdg = stable_leaf_pdgs[leaf_index]
        if node.particle in selectors or node.particle in aliases:
            allowed_pdgs = resolve_particle_selector(node.particle, selectors, aliases)
            if stable_pdg not in allowed_pdgs:
                raise RuntimeError(
                    f"Stable PDG {stable_pdg} disagrees with particle token {node.particle}"
                )
        allowed_pdgs = (stable_pdg,)
        return {
            "allowed_pdgs": list(allowed_pdgs),
            "daughters": [],
        }, leaf_index + 1

    allowed_pdgs = resolve_particle_selector(node.particle, selectors, aliases)
    daughters: list[dict[str, Any]] = []
    next_index = leaf_index
    for daughter in node.daughters:
        converted, next_index = topology_node(
            daughter, selectors, aliases, stable_leaf_pdgs, next_index
        )
        daughters.append(converted)
    return {
        "allowed_pdgs": list(allowed_pdgs),
        "daughters": daughters,
    }, next_index


# Build the complete ordered numeric topology for one parsed MG5 process
def numeric_topology(
    parsed: ParsedProcessSyntax,
    selectors: dict[str, tuple[int, ...]],
    aliases: dict[str, tuple[str, ...]] | None = None,
    stable_leaf_pdgs: list[int] | None = None,
) -> list[dict[str, Any]]:
    topology: list[dict[str, Any]] = []
    next_index = 0
    for node in parsed.production:
        converted, next_index = topology_node(
            node, selectors, aliases or {}, stable_leaf_pdgs, next_index
        )
        topology.append(converted)
    if stable_leaf_pdgs is not None and next_index != len(stable_leaf_pdgs):
        raise RuntimeError("Stable final-state PDG list is longer than the topology")
    return topology


# Append the stable selector projection of one numeric topology node
def append_stable_pattern(node: dict[str, Any], stable_pdgs: list[int]) -> bool:
    allowed_pdgs = node["allowed_pdgs"]
    has_selector = len(allowed_pdgs) > 1
    daughters = node["daughters"]
    if not daughters:
        stable_pdgs.append(allowed_pdgs[0] if len(allowed_pdgs) == 1 else 0)
        return has_selector
    for daughter in daughters:
        has_selector = append_stable_pattern(daughter, stable_pdgs) or has_selector
    return has_selector


# Compute stable leaf selectors and whether any topology node is non-exact
def topology_stable_pattern(topology: list[dict[str, Any]]) -> tuple[list[int], bool]:
    stable_pdgs: list[int] = []
    has_selector = False
    for node in topology:
        has_selector = append_stable_pattern(node, stable_pdgs) or has_selector
    return stable_pdgs, has_selector


# Compute one strict projection selected by a manifest entry
def manifest_projection(entry: dict[str, Any], allowed: frozenset[str], context: str) -> str:
    projection = entry.get("projection")
    if not isinstance(projection, str) or projection not in allowed:
        choices = ", ".join(sorted(allowed))
        raise RuntimeError(f"{context} requires projection in {{{choices}}}")
    return projection


# Compute the projection selected by one standalone process
def standalone_projection(entry: dict[str, Any]) -> str:
    return manifest_projection(
        entry,
        STANDALONE_PROJECTIONS,
        f"Standalone process {entry.get('name', '<unknown>')}",
    )


# Compute the projection selected by one generated family
def family_projection(family: dict[str, Any]) -> str:
    return manifest_projection(
        family,
        FAMILY_PROJECTIONS,
        f"Family {family.get('name', '<unknown>')}",
    )


# Derive one standalone beam PDG from its projection and parsed incoming state
def standalone_incoming_pdg(entry: dict[str, Any]) -> int:
    projection = standalone_projection(entry)
    expected = {
        DURHAM_PROJECTION: (("g", "g"), 21),
        PHOTON_PROJECTION: (("a", "a"), 22),
    }[projection]
    process = entry.get("process")
    if not isinstance(process, str) or not process:
        raise RuntimeError(
            f"Standalone process {entry.get('name', '<unknown>')} has invalid MG5 syntax"
        )
    incoming = parse_process_syntax(process).incoming
    if incoming != expected[0]:
        state = " ".join(expected[0])
        raise RuntimeError(
            f"Standalone process {entry.get('name', '<unknown>')} projection "
            f"{projection} requires incoming state {state}"
        )
    return expected[1]


# Load generated standalone process data containing stable final-state PDGs
def load_standalone_process_data(base_dir: Path, entry: dict[str, Any]) -> dict[str, Any]:
    projection = standalone_projection(entry)
    data_path = base_dir / output_layout.color_path(projection, entry["name"])
    if not data_path.is_file():
        raise RuntimeError(f"Missing generated MG5 process data: {data_path}")
    data = load_json_file(data_path)
    if not isinstance(data, dict):
        raise RuntimeError(f"Invalid generated MG5 process data for {entry['name']}")
    if data.get("process") != entry["process"]:
        raise RuntimeError(f"Stale generated MG5 process data for {entry['name']}")
    final_pdgs = data.get("final_pdgs")
    if (
        not isinstance(final_pdgs, list)
        or not final_pdgs
        or any(type(pdg) is not int or pdg == 0 or abs(pdg) == 89 for pdg in final_pdgs)
    ):
        raise RuntimeError(f"Invalid stable final-state PDGs for {entry['name']}")
    if "has_decay_chain" in data and not isinstance(data["has_decay_chain"], bool):
        raise RuntimeError(f"Invalid decay-chain flag for {entry['name']}")
    return data


# Convert one exact family channel branch to a numeric process node
def family_channels_topology_node(data: Any, context: str) -> dict[str, Any]:
    if not isinstance(data, dict) or set(data) != {"pdg", "daughters"}:
        raise RuntimeError(f"{context} must contain pdg and daughters")
    pdg = data["pdg"]
    daughters = data["daughters"]
    if type(pdg) is not int or pdg == 0 or abs(pdg) == 89:
        raise RuntimeError(f"{context} requires one exact nonzero PDG id")
    if not isinstance(daughters, list):
        raise RuntimeError(f"{context} daughters must be an ordered list")
    return {
        "allowed_pdgs": [pdg],
        "daughters": [
            family_channels_topology_node(daughter, f"{context} daughter {index}")
            for index, daughter in enumerate(daughters)
        ],
    }


# Load exact generated family topologies from persistent subprocess channels
def load_family_channel_topologies(
    base_dir: Path, family: dict[str, Any], model: dict[str, Any]
) -> list[list[dict[str, Any]]]:
    family_name = family["name"]
    path = base_dir / output_layout.channels_path(family_projection(family), family_name)
    if not path.is_file():
        raise RuntimeError(f"Missing generated family channels: {path}")
    family_channels = load_json_file(path)
    if (
        not isinstance(family_channels, dict)
        or set(family_channels)
        != {"schema", "family", "generation_setup", "mg5_source", "subprocesses"}
        or type(family_channels["schema"]) is not int
        or family_channels["schema"] != FAMILY_CHANNELS_SCHEMA
        or family_channels["family"] != family_name
        or not isinstance(family_channels["subprocesses"], list)
        or not family_channels["subprocesses"]
    ):
        raise RuntimeError(f"Invalid generated family channels for {family_name}")
    generation_setup = family_generation_setup_from_data(
        family_channels["generation_setup"],
        f"Generated family {family_name} generation setup",
    )
    expected_generation = family_generation_setup_from_manifest(family, model)
    if generation_setup != expected_generation:
        raise RuntimeError(f"Stale generated family generation setup for {family_name}")
    mg5_source_from_data(
        family_channels["mg5_source"], f"Generated family {family_name} MG5 source"
    )

    generated: list[list[dict[str, Any]]] = []
    channel_fields = {
        "initial",
        "final",
        "topology",
        "external_color_representations",
        "external_color_flows",
    }
    for subprocess_index, subprocess in enumerate(family_channels["subprocesses"]):
        context = f"Generated family {family_name} subprocess {subprocess_index}"
        if (
            not isinstance(subprocess, dict)
            or set(subprocess) != {"generated_name", "channels"}
            or not isinstance(subprocess["generated_name"], str)
            or not subprocess["generated_name"]
            or not isinstance(subprocess["channels"], list)
            or not subprocess["channels"]
        ):
            raise RuntimeError(f"{context} is invalid")
        for channel_index, channel in enumerate(subprocess["channels"]):
            channel_context = f"{context} channel {channel_index}"
            if not isinstance(channel, dict) or set(channel) != channel_fields:
                raise RuntimeError(f"{channel_context} has invalid fields")
            initial = channel["initial"]
            final = channel["final"]
            topology_data = channel["topology"]
            representations = channel["external_color_representations"]
            flows = channel["external_color_flows"]
            if (
                not isinstance(initial, list)
                or len(initial) != 2
                or any(type(pdg) is not int or pdg == 0 or abs(pdg) == 89 for pdg in initial)
                or not isinstance(final, list)
                or not final
                or any(type(pdg) is not int or pdg == 0 or abs(pdg) == 89 for pdg in final)
                or not isinstance(topology_data, list)
                or not topology_data
            ):
                raise RuntimeError(f"{channel_context} has invalid exact particles")
            topology = [
                family_channels_topology_node(node, f"{channel_context} topology node {index}")
                for index, node in enumerate(topology_data)
            ]
            topology_final, has_selector = topology_stable_pattern(topology)
            if has_selector or topology_final != final:
                raise RuntimeError(f"{channel_context} topology and stable final state disagree")
            external_count = len(final) + 2
            if (
                not isinstance(representations, list)
                or len(representations) != external_count
                or any(
                    type(value) is not int or value not in {1, 3, -3, 8}
                    for value in representations
                )
                or not isinstance(flows, list)
                or not flows
                or any(
                    not isinstance(flow, list)
                    or len(flow) != 2 * external_count
                    or any(type(value) is not int for value in flow)
                    for flow in flows
                )
            ):
                raise RuntimeError(f"{channel_context} has invalid external color structure")
            if topology not in generated:
                generated.append(topology)
    return generated


# Assign every exact generated channel to one narrowest advertised source pattern
def assign_channel_topologies(
    family_name: str,
    entries: list[dict[str, Any]],
    channel_topologies: list[list[dict[str, Any]]],
) -> None:
    for generated in channel_topologies:
        candidates = [
            index
            for index, data in enumerate(entries)
            if topology_is_subset(generated, data["topology"])
        ]
        narrowest = [
            index
            for index in candidates
            if not any(
                topology_is_strict_subset(entries[other]["topology"], entries[index]["topology"])
                for other in candidates
                if other != index
            )
        ]
        if len(narrowest) != 1:
            raise RuntimeError(
                f"Generated family {family_name} topology maps to {len(narrowest)} source patterns"
            )
        entries[narrowest[0]]["channel_topologies"].append(generated)
    for data in entries:
        if not data["channel_topologies"]:
            raise RuntimeError(
                f"Generated family {family_name} source {data['process_syntax']} "
                "has no exact channels"
            )


# Compute the standalone process family
def standalone_process_family(entry: dict[str, Any]) -> str:
    projection = standalone_projection(entry)
    standalone_incoming_pdg(entry)
    return {
        DURHAM_PROJECTION: "DURHAM",
        PHOTON_PROJECTION: "PHOTON",
    }[projection]


# Compute one generated process family in the reserved MG5 namespace
def family_process_family(family: dict[str, Any]) -> str:
    name = family.get("name")
    if not isinstance(name, str) or not name:
        raise RuntimeError("Generated process family is invalid")
    if name in STANDALONE_PROCESS_FAMILIES:
        raise RuntimeError(f"Generated process family {name} is reserved")
    if not name.startswith(PROCESS_FAMILY_PREFIX):
        raise RuntimeError(
            f"Generated process family {name!r} must use the {PROCESS_FAMILY_PREFIX} namespace"
        )
    suffix = name.removeprefix(PROCESS_FAMILY_PREFIX)
    if re.fullmatch(r"[A-Za-z][A-Za-z0-9_]*", suffix) is None:
        raise RuntimeError(f"Generated process family {name!r} has an invalid suffix")
    return name


# Derive one public wrapper class from its generated MG5 family name
def family_wrapper_name(family: dict[str, Any]) -> str:
    name = family_process_family(family)
    suffix = name.removeprefix(PROCESS_FAMILY_PREFIX).lower()
    return f"AMP_MG5_{suffix}"


# Build one JSON-ready process data
def process_data(
    process_family: str,
    process_name: str,
    process_syntax: str,
    parsed: ParsedProcessSyntax,
    topology: list[dict[str, Any]],
    allows_isolated_resonance: bool,
    channel_topologies: list[list[dict[str, Any]]] | None = None,
) -> dict[str, Any]:
    stable_pdgs, has_selector = topology_stable_pattern(topology)
    exact_topologies = channel_topologies
    if exact_topologies is None:
        exact_topologies = [] if has_selector else [topology]
    return {
        "process_family": process_family,
        "process_name": process_name,
        "stable_final_state": stable_final_state(parsed),
        "final_state_syntax": final_state_syntax(parsed),
        "topology_signature": topology_signature(parsed),
        "process_syntax": process_syntax,
        "topology": topology,
        "channel_topologies": exact_topologies,
        "stable_pdgs": stable_pdgs,
        "decay_structure": {
            "type": "full",
            "allows_isolated_resonance": allows_isolated_resonance,
        },
        "matrix_element_form": "generated",
        "topology_mode": "flavour_set" if has_selector else "exact",
    }


# Convert one MG5 token to a stable process-name component
def process_name_token(token: str) -> str:
    suffix = "bar" if token.endswith("~") else ""
    core = token.removesuffix("~").replace("+", "p").replace("-", "m")
    return re.sub(r"[^a-zA-Z0-9]", "", core).lower() + suffix


# Compute the stable base name for one generated MG5 family
def family_base_name(family: dict[str, Any]) -> str:
    name = family_process_family(family).removeprefix(PROCESS_FAMILY_PREFIX)
    normalized = re.sub(r"[^a-zA-Z0-9]+", "_", name).strip("_").lower()
    if not normalized:
        raise RuntimeError(f"Cannot derive a process name for family {family['name']}")
    return normalized


# Compute explicit names or derive membership-independent full-topology names
def family_process_names(
    family: dict[str, Any], parsed_rows: list[ParsedProcessSyntax]
) -> list[str]:
    explicit = family.get("process_names")
    if explicit is not None:
        if (
            not isinstance(explicit, list)
            or len(explicit) != len(parsed_rows)
            or len(explicit) != len(set(explicit))
            or any(
                not isinstance(name, str) or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", name) is None
                for name in explicit
            )
        ):
            raise RuntimeError(f"Family {family['name']} has invalid process_names")
        return list(explicit)

    base = family_base_name(family)
    names: list[str] = []
    for parsed in parsed_rows:
        suffix = "".join(process_name_token(token) for token in stable_tokens(parsed))
        signature = topology_signature(parsed)
        digest = hashlib.sha256(signature.encode("ascii")).hexdigest()[:12]
        prefix = f"{base}_{suffix}" if suffix else base
        names.append(f"{prefix}_topology_{digest}")
    if len(names) != len(set(names)):
        raise RuntimeError(f"Cannot derive unique process names for family {family['name']}")
    return names


# Derive every standalone and family process from the converter inputs
def build_registry(base_dir: Path, manifest: dict[str, Any]) -> dict[str, Any]:
    entries: list[dict[str, Any]] = []
    for entry in manifest["processes"]:
        source = str(entry["process"])
        parsed = parse_process_syntax(source)
        standalone_data = load_standalone_process_data(base_dir, entry)
        has_decay_chain = standalone_data.get("has_decay_chain")
        if has_decay_chain is not None and has_decay_chain != parsed.has_decay_chain:
            raise RuntimeError(f"Decay-chain process data mismatch for {entry['name']}")
        pdgs = list(standalone_data["final_pdgs"])
        if len(pdgs) != len(stable_tokens(parsed)):
            raise RuntimeError(f"Stable final-state size mismatch for {entry['name']}")
        selectors = model_particle_selectors(manifest, entry["model"])
        topology = numeric_topology(parsed, selectors, stable_leaf_pdgs=pdgs)
        entries.append(
            process_data(
                standalone_process_family(entry),
                entry["name"],
                source,
                parsed,
                topology,
                allows_isolated_resonance=not parsed.has_decay_chain,
            )
        )

    for family in manifest.get("families", []):
        family_projection(family)
        process_family = family_process_family(family)
        sources = [str(source) for source in family["processes"]]
        parsed_rows = [parse_process_syntax(source) for source in sources]
        names = family_process_names(family, parsed_rows)
        selectors = model_particle_selectors(manifest, family["model"])
        aliases = parse_family_aliases(family.get("definitions", []))
        family_datas = [
            process_data(
                process_family,
                name,
                source,
                parsed,
                numeric_topology(parsed, selectors, aliases),
                allows_isolated_resonance=False,
                channel_topologies=[],
            )
            for name, source, parsed in zip(names, sources, parsed_rows, strict=True)
        ]
        assign_channel_topologies(
            process_family,
            family_datas,
            load_family_channel_topologies(base_dir, family, manifest["models"][family["model"]]),
        )
        entries.extend(family_datas)
    validate_registry(entries)
    return {"version": 1, "processes": entries}


# Validate one numeric topology node and return its stable selector projection
def validate_topology_node(node: Any, context: str) -> tuple[list[int], bool]:
    if not isinstance(node, dict) or set(node) != {"allowed_pdgs", "daughters"}:
        raise RuntimeError(f"{context} must contain allowed_pdgs and daughters")
    allowed_pdgs = node["allowed_pdgs"]
    if (
        not isinstance(allowed_pdgs, list)
        or not allowed_pdgs
        or any(
            isinstance(pdg, bool) or not isinstance(pdg, int) or pdg == 0 for pdg in allowed_pdgs
        )
        or allowed_pdgs != sorted(set(allowed_pdgs))
    ):
        raise RuntimeError(f"{context} selectors must be nonempty, sorted, unique and nonzero")
    daughters = node["daughters"]
    if not isinstance(daughters, list):
        raise RuntimeError(f"{context} daughters must be an ordered list")
    has_selector = len(allowed_pdgs) > 1
    if not daughters:
        return [allowed_pdgs[0] if len(allowed_pdgs) == 1 else 0], has_selector
    stable_pdgs: list[int] = []
    for index, daughter in enumerate(daughters):
        stable, daughter_selector = validate_topology_node(daughter, f"{context} daughter {index}")
        stable_pdgs.extend(stable)
        has_selector = has_selector or daughter_selector
    return stable_pdgs, has_selector


# Compute whether two numeric topology nodes accept a common branch
def topology_nodes_overlap(first: dict[str, Any], second: dict[str, Any]) -> bool:
    if len(first["daughters"]) != len(second["daughters"]):
        return False
    if not set(first["allowed_pdgs"]).intersection(second["allowed_pdgs"]):
        return False
    return all(
        topology_nodes_overlap(first_daughter, second_daughter)
        for first_daughter, second_daughter in zip(
            first["daughters"], second["daughters"], strict=True
        )
    )


# Compute whether two ordered numeric topologies accept a common decay tree
def topologies_overlap(first: list[dict[str, Any]], second: list[dict[str, Any]]) -> bool:
    return len(first) == len(second) and all(
        topology_nodes_overlap(first_node, second_node)
        for first_node, second_node in zip(first, second, strict=True)
    )


# Compute whether every branch accepted by the first node is accepted by the second
def topology_node_is_subset(first: dict[str, Any], second: dict[str, Any]) -> bool:
    if len(first["daughters"]) != len(second["daughters"]):
        return False
    if not set(first["allowed_pdgs"]).issubset(second["allowed_pdgs"]):
        return False
    return all(
        topology_node_is_subset(first_daughter, second_daughter)
        for first_daughter, second_daughter in zip(
            first["daughters"], second["daughters"], strict=True
        )
    )


# Compute whether every branch in one ordered topology is accepted by another
def topology_is_subset(first: list[dict[str, Any]], second: list[dict[str, Any]]) -> bool:
    return len(first) == len(second) and all(
        topology_node_is_subset(first_node, second_node)
        for first_node, second_node in zip(first, second, strict=True)
    )


# Compute whether one ordered numeric topology is a strict subset of another
def topology_is_strict_subset(first: list[dict[str, Any]], second: list[dict[str, Any]]) -> bool:
    return first != second and topology_is_subset(first, second)


# Validate process keys, topologies and process family intersections
def validate_registry(entries: list[dict[str, Any]]) -> None:
    names: set[tuple[str, str]] = set()
    for data in entries:
        key = (data["process_family"], data["process_name"])
        if key in names:
            raise RuntimeError(f"Duplicate MG5 process key: {key}")
        names.add(key)
        if not all(
            data[field]
            for field in (
                "process_family",
                "process_name",
                "stable_final_state",
                "final_state_syntax",
                "topology_signature",
                "process_syntax",
            )
        ):
            raise RuntimeError(f"Incomplete MG5 process: {key}")
        topology = data.get("topology")
        if not isinstance(topology, list) or not topology:
            raise RuntimeError(f"MG5 process {key} requires a nonempty topology")
        expected_stable: list[int] = []
        has_selector = False
        for index, node in enumerate(topology):
            stable, node_selector = validate_topology_node(
                node, f"MG5 process {key} topology node {index}"
            )
            expected_stable.extend(stable)
            has_selector = has_selector or node_selector
        if data.get("stable_pdgs") != expected_stable:
            raise RuntimeError(f"MG5 process {key} stable_pdgs disagree with its topology")
        expected_topology_mode = "flavour_set" if has_selector else "exact"
        if (
            data.get("matrix_element_form") != "generated"
            or data.get("topology_mode") != expected_topology_mode
        ):
            raise RuntimeError(f"MG5 process {key} topology mode disagrees with its selectors")
        channel_topologies = data.get("channel_topologies")
        if not isinstance(channel_topologies, list):
            raise RuntimeError(f"MG5 process {key} has invalid generated topologies")
        for exact_index, exact in enumerate(channel_topologies):
            if not isinstance(exact, list) or not exact:
                raise RuntimeError(f"MG5 process {key} has an empty generated topology")
            exact_has_selector = False
            for node_index, node in enumerate(exact):
                _, node_selector = validate_topology_node(
                    node,
                    f"MG5 process {key} generated topology {exact_index} node {node_index}",
                )
                exact_has_selector = exact_has_selector or node_selector
            if exact_has_selector:
                raise RuntimeError(
                    f"MG5 process {key} generated topology {exact_index} is not exact"
                )
            if not topology_is_subset(exact, topology):
                raise RuntimeError(
                    f"MG5 process {key} generated topology {exact_index} is outside "
                    "its selector pattern"
                )
            if exact in channel_topologies[:exact_index]:
                raise RuntimeError(f"MG5 process {key} has a duplicate generated topology")
        decay_structure = data.get("decay_structure")
        if (
            not isinstance(decay_structure, dict)
            or set(decay_structure) != {"type", "allows_isolated_resonance"}
            or decay_structure["type"] != "full"
            or not isinstance(decay_structure["allows_isolated_resonance"], bool)
        ):
            raise RuntimeError(f"MG5 process {key} must declare a complete decay amplitude")

    for first_index, first in enumerate(entries):
        for second in entries[first_index + 1 :]:
            if first["process_family"] != second["process_family"]:
                continue
            first_topology = first["topology"]
            second_topology = second["topology"]
            if not topologies_overlap(first_topology, second_topology):
                continue
            if first_topology == second_topology:
                signature = (first["process_family"], first["topology_signature"])
                raise RuntimeError(f"Duplicate MG5 process topology: {signature}")
            if not topology_is_strict_subset(
                first_topology, second_topology
            ) and not topology_is_strict_subset(second_topology, first_topology):
                raise RuntimeError(
                    "Overlapping incomparable MG5 process topologies: "
                    f"{first['process_name']} and {second['process_name']}"
                )


# Serialize the generated JSON with compact entries and PDG arrays
def registry_json(registry: dict[str, Any]) -> str:
    entries = [
        "    " + json.dumps(data, ensure_ascii=True, separators=(", ", ": "))
        for data in registry["processes"]
    ]
    return (
        "{\n"
        f'  "version": {int(registry["version"])},\n'
        '  "processes": [\n' + ",\n".join(entries) + "\n  ]\n}\n"
    )


# Escape one generated C++ string literal
def cpp_string(value: str) -> str:
    return json.dumps(value, ensure_ascii=True)


# Format one generated C++ integer vector
def cpp_pdgs(values: list[int]) -> str:
    return "{" + ", ".join(str(int(value)) for value in values) + "}"


# Format one generated C++ boolean literal
def cpp_bool(value: bool) -> str:
    return "true" if value else "false"


# Format one generated nested C++ topology node initializer
def cpp_topology_node(node: dict[str, Any]) -> str:
    daughters = ", ".join(cpp_topology_node(child) for child in node["daughters"])
    return (
        "AmplitudeTopologyNode{"
        f"std::vector<int>{cpp_pdgs(node['allowed_pdgs'])}, "
        f"std::vector<AmplitudeTopologyNode>{{{daughters}}}"
        "}"
    )


# Format one generated ordered C++ topology initializer
def cpp_topology(topology: list[dict[str, Any]]) -> str:
    nodes = ", ".join(cpp_topology_node(node) for node in topology)
    return f"AmplitudeTopology{{{nodes}}}"


# Format one generated C++ vector of exact ordered topologies
def cpp_topologies(topologies: list[list[dict[str, Any]]]) -> str:
    values = ", ".join(cpp_topology(topology) for topology in topologies)
    return f"std::vector<AmplitudeTopology>{{{values}}}"


# Compute a C++ identifier fragment for one generated registry key
def cpp_identifier(value: str) -> str:
    identifier = re.sub(r"[^a-zA-Z0-9_]", "_", value)
    if not identifier or identifier[0].isdigit():
        identifier = "_" + identifier
    return identifier


# Generate one immutable registry base backed by converter process lookups
def cpp_registry_class(
    class_name: str, process_family: str, process_name: str | None = None
) -> str:
    arguments = cpp_string(process_family)
    if process_name is not None:
        arguments += ", " + cpp_string(process_name)
    return f"""// Own converter-generated processes for {process_family}
class {class_name}
    : public ProcessRegistry {{
 public:
  // Construct the immutable generated process view
  {class_name}()
      : ProcessRegistry(
            ::gra::amplitude::Processes({arguments})) {{}}
}};
"""


# Generate exact standalone registry bases used by raw amplitudes
def cpp_registry_classes(registry: dict[str, Any]) -> str:
    entries = registry["processes"]
    classes = list(
        cpp_registry_class(
            f"MG5ProcessRegistry_{cpp_identifier(str(data['process_name']))}",
            str(data["process_family"]),
            str(data["process_name"]),
        )
        for data in entries
        if data["process_family"] in STANDALONE_PROCESS_FAMILIES
    )
    return "\n".join(classes).rstrip()


# Format one generated C++ process initializer
def cpp_process(data: dict[str, Any]) -> str:
    allows_isolated_resonance = cpp_bool(data["decay_structure"]["allows_isolated_resonance"])
    matrix_element_form = {"generated": "Generated"}[data["matrix_element_form"]]
    topology_mode = {"exact": "Exact", "flavour_set": "FlavourSet"}[data["topology_mode"]]
    fields = [
        cpp_string(data["process_family"]),
        cpp_string(data["process_name"]),
        cpp_string(data["stable_final_state"]),
        cpp_string(data["final_state_syntax"]),
        cpp_string(data["topology_signature"]),
        cpp_string(data["process_syntax"]),
        cpp_topology(data["topology"]),
        cpp_topologies(data["channel_topologies"]),
        f"std::vector<int>{cpp_pdgs(data['stable_pdgs'])}",
        "DecayStructure{"
        "DecayType::Full, "
        f"{allows_isolated_resonance}"
        "}",
        f"MatrixElementForm::{matrix_element_form}",
        f"TopologyMode::{topology_mode}",
    ]
    return "      Process{" + ",\n       ".join(fields) + "}"


# Generate the value-returning common registry declaration
def registry_header(registry: dict[str, Any]) -> str:
    registries = cpp_registry_classes(registry)
    return f"""// Generated MadGraph amplitude process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#ifndef AMP_MG5_PROCESSREGISTRY_H
#define AMP_MG5_PROCESSREGISTRY_H

#include <optional>
#include <string>
#include <vector>

#include "Graniitti/Process/MAmpMatch.h"

namespace gra::amplitude {{

// Compute all generated MG5 amplitude processes by value
std::vector<Process> AllProcesses();

// Compute processes in one generated amplitude family by value
std::vector<Process> Processes(
    const std::string &process_family);

// Compute one exact generated process as a vector
std::vector<Process> Processes(
    const std::string &process_family, const std::string &process_name);

// Find one generated process by its family and process names
std::optional<Process> FindProcess(
    const std::string &process_family, const std::string &process_name);

// Find one generated process by typed topology and ordered stable PDGs
std::optional<Process> FindProcess(
    const std::string &process_family, const AmplitudeTopology &topology,
    const std::vector<int> &stable_pdgs);

// Compute the parameter card owned by one exact generated process
std::optional<std::string> ParameterCard(
    const std::string &process_family, const std::string &process_name);

{registries}

}}  // namespace gra::amplitude

#endif
"""


# Generate the immutable-by-construction common registry implementation
def registry_source(registry: dict[str, Any], manifest: dict[str, Any]) -> str:
    entries = ",\n".join(cpp_process(data) for data in registry["processes"])
    card_owners = [
        (
            standalone_process_family(entry),
            entry["name"],
            output_layout.parameter_card_path(entry["projection"], entry["name"]),
        )
        for entry in manifest["processes"]
    ]
    for family in manifest.get("families", []):
        process_family = family_process_family(family)
        card = output_layout.parameter_card_path(family["projection"], family["name"])
        card_owners.extend(
            (process_family, process["process_name"], card)
            for process in registry["processes"]
            if process["process_family"] == process_family
        )
    parameter_cards = "".join(
        f"  if (process_family == {cpp_string(process_family)} &&\n"
        f"      process_name == {cpp_string(process_name)}) {{\n"
        f"    return {cpp_string(output_layout.source_path(card))};\n"
        "  }\n"
        for process_family, process_name, card in card_owners
    )
    return f"""// Generated MadGraph amplitude process registry
//
// (c) 2026 Mikael Mieskolainen
// Licensed under the MIT License <http://opensource.org/licenses/MIT>.

#include "{output_layout.include_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h")}"

#include <utility>

namespace gra::amplitude {{

// Compute all generated MG5 amplitude processes by value
std::vector<Process> AllProcesses() {{
  return {{
{entries}
  }};
}}

// Compute processes in one generated amplitude family by value
std::vector<Process> Processes(
    const std::string &process_family) {{
  std::vector<Process> selected;
  for (auto process : AllProcesses()) {{
    if (process.process_family == process_family) {{
      selected.push_back(std::move(process));
    }}
  }}
  return selected;
}}

// Compute one exact generated process as a vector
std::vector<Process> Processes(
    const std::string &process_family, const std::string &process_name) {{
  const auto process =
      FindProcess(process_family, process_name);
  if (!process.has_value()) {{ return {{}}; }}
  return {{*process}};
}}

// Find one generated process by its family and process names
std::optional<Process> FindProcess(
    const std::string &process_family, const std::string &process_name) {{
  for (auto process : AllProcesses()) {{
    if (process.process_family == process_family &&
        process.process_name == process_name) {{
      return process;
    }}
  }}
  return std::nullopt;
}}

// Compute the parameter card owned by one exact generated process
std::optional<std::string> ParameterCard(
    const std::string &process_family, const std::string &process_name) {{
{parameter_cards}  return std::nullopt;
}}

// Find one generated process by typed topology and ordered stable PDGs
std::optional<Process> FindProcess(
    const std::string &process_family, const AmplitudeTopology &topology,
    const std::vector<int> &stable_pdgs) {{
  for (auto process : AllProcesses()) {{
    if (process.process_family == process_family &&
        process.topology == topology &&
        process.stable_pdgs == stable_pdgs) {{
      return process;
    }}
  }}
  return std::nullopt;
}}

}}  // namespace gra::amplitude
"""


# Generate JSON and C++ views from the same converter-fed process values
def generate_registry(base_dir: Path, manifest: dict[str, Any]) -> tuple[str, str, str]:
    registry = build_registry(base_dir, manifest)
    return registry_json(registry), registry_header(registry), registry_source(registry, manifest)
