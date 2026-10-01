#!/usr/bin/env python3
#
# Fully regenerate standalone MG5 amplitudes and shared model support
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import argparse
import copy
import errno
import hashlib
import json
import math
import os
import re
import shutil
import stat
import subprocess
import tempfile
from pathlib import Path
from typing import Any

from core.io.files import ensure_dir
from core.io.serialize import load_json_file
from modules import (
    durham_registry,
    madloop,
    mg5_color,
    mg5_family,
    output_layout,
    parton_registry,
    photon_registry,
    process_registry,
    standalone,
)
from modules.mg5_support import (
    ModelFiles,
    add_model_particles,
    discover_process_model,
    supports_alpha_qed_zero,
    transform_helas_header,
    transform_helas_source,
    transform_parameters_header,
    transform_parameters_source,
    transform_slha_header,
    transform_slha_source,
    validate_dependency_union,
)

BASE_DIR = Path(__file__).resolve().parent
ROOT = BASE_DIR.parents[1]
DEFAULT_MANIFEST = BASE_DIR / "processes.json"
CARDS_ROOT = ROOT / output_layout.CARDS_ROOT
CONVERTER_STATE_PATH = CARDS_ROOT / output_layout.converter_path()
CONVERTER_STATE_SCHEMA = 1
CPP_IDENTIFIER = re.compile(r"[A-Za-z_][A-Za-z0-9_]*")
MODEL_IMPORT = re.compile(r"[A-Za-z0-9_./+~-]+")
SUPPORTED_MG5_VERSIONS = {"2.9.27"}
PARTICLE_NAME = re.compile(r"[A-Za-z][A-Za-z0-9_+~-]*")
MG5_DEFINITION = re.compile(
    r"define[ \t]+[A-Za-z][A-Za-z0-9_+~-]*[ \t]*=[ \t]*"
    r"[A-Za-z0-9_+~-]+(?:[ \t]+[A-Za-z0-9_+~-]+)*"
)
FIXED_AMPLITUDE_SUFFIXES = frozenset(
    {
        "DurhamRegistry",
        "Helicity",
        "PartonRegistry",
        "PhotonRegistry",
        "Process",
        "SubprocessSum",
        "Utils",
        "ZProcessUtils",
    }
)
RESERVED_OUTPUT_NAMES = FIXED_AMPLITUDE_SUFFIXES | {"family_channels"}
FAMILY_WRAPPER_INCLUDE = re.compile(
    r'/Amplitude/MG5/(?:Photon|Parton)/(MG5_[A-Za-z0-9_]+)/(?:ProcessBase|Processes)\.h"'
)


# Compute one converter-owned generated header below the selected repository root
def cpp_header_path(directory: str, filename: str, root: Path | None = None) -> Path:
    return (ROOT if root is None else root) / output_layout.header_path(directory, filename)


# Compute one converter-owned generated source below the selected repository root
def cpp_source_path(directory: str, filename: str, root: Path | None = None) -> Path:
    return (ROOT if root is None else root) / output_layout.cpp_source_path(directory, filename)


# Compute the physical C++ directory of one standalone process
def standalone_cpp_directory(entry: dict[str, Any]) -> str:
    return output_layout.projection_cpp_directory(entry["projection"])


# Compute the physical C++ directory of one process family
def family_cpp_directory(family: dict[str, Any]) -> str:
    return output_layout.family_cpp_directory(family["projection"], family["name"])


# Compute true when one JSON value is compact enough for a single line
def is_compact_json_value(value: Any) -> bool:
    if isinstance(value, list):
        return all(not isinstance(item, (dict, list)) for item in value)
    if isinstance(value, dict):
        return all(not isinstance(item, (dict, list)) for item in value.values())
    return True


# Format JSON with compact scalar arrays and entries
def format_json_value(value: Any, level: int = 0) -> str:
    if not isinstance(value, (dict, list)) or is_compact_json_value(value):
        return json.dumps(value, ensure_ascii=False, separators=(", ", ": "))

    indentation = "  " * (level + 1)
    if isinstance(value, list):
        items = [indentation + format_json_value(item, level + 1) for item in value]
        return "[\n" + ",\n".join(items) + "\n" + "  " * level + "]"

    items = [
        indentation + json.dumps(key) + ": " + format_json_value(item, level + 1)
        for key, item in value.items()
    ]
    return "{\n" + ",\n".join(items) + "\n" + "  " * level + "}"


# Serialize the process manifest without expanding compact arrays
def dump_manifest(manifest: dict[str, Any]) -> str:
    return format_json_value(manifest) + "\n"


# Compute the parameter card owned by one manifest process or family
def parameter_card_path(entry: dict[str, Any]) -> Path:
    return CARDS_ROOT / output_layout.parameter_card_path(entry["projection"], entry["name"])


# Compute the finite color data owned by one standalone process
def standalone_color_path(entry: dict[str, Any]) -> Path:
    return CARDS_ROOT / output_layout.color_path(entry["projection"], entry["name"])


# Compute the persistent exact channel path for one process family
def family_channels_path(family: dict[str, Any]) -> Path:
    return CARDS_ROOT / output_layout.channels_path(family["projection"], family["name"])


# Load persistent exact channels and verify their generation setup
def load_family_channels(
    family: dict[str, Any], model: dict[str, Any]
) -> mg5_family.FamilyChannels:
    family_name = family["name"]
    path = family_channels_path(family)
    if not path.is_file():
        raise RuntimeError(f"Missing generated family channels for {family_name}")
    family_channels = mg5_family.family_channels_from_data(load_json_file(path))
    if family_channels.family != family_name:
        raise RuntimeError(f"Generated family channel name mismatch for {family_name}")
    expected_generation = process_registry.family_generation_setup_from_manifest(family, model)
    if family_channels.generation_setup != expected_generation:
        raise RuntimeError(f"Stale generated family generation setup for {family_name}")
    return family_channels


# Require one JSON object at a manifest boundary
def require_manifest_object(value: Any, context: str) -> dict[str, Any]:
    if not isinstance(value, dict):
        raise RuntimeError(f"{context} must be an object")
    return value


# Reject misspelled or unsupported fields at one manifest object boundary
def require_manifest_keys(
    value: dict[str, Any],
    required: set[str],
    optional: set[str],
    context: str,
) -> None:
    missing = required - set(value)
    unknown = set(value) - required - optional
    if missing:
        raise RuntimeError(f"{context} lacks fields {sorted(missing)}")
    if unknown:
        raise RuntimeError(f"{context} has unknown fields {sorted(unknown)}")


# Require one present manifest field
def require_manifest_field(data: dict[str, Any], field: str, context: str) -> Any:
    if field not in data:
        raise RuntimeError(f"{context} requires {field}")
    return data[field]


# Require one safe C++ and path identifier
def require_manifest_identifier(value: Any, context: str) -> str:
    if not isinstance(value, str) or CPP_IDENTIFIER.fullmatch(value) is None:
        raise RuntimeError(f"{context} must be a safe nonempty identifier")
    return value


# Require one nonempty single-line manifest string
def require_manifest_string(value: Any, context: str) -> str:
    if (
        not isinstance(value, str)
        or not value
        or value.strip() != value
        or any(character in value for character in "\r\n\0")
    ):
        raise RuntimeError(f"{context} must be a nonempty single-line string")
    return value


# Require one JSON list containing strings
def require_manifest_string_list(value: Any, context: str, *, allow_empty: bool) -> list[str]:
    if not isinstance(value, list) or (not allow_empty and not value):
        raise RuntimeError(f"{context} must be a {'nonempty ' if not allow_empty else ''}list")
    return [
        require_manifest_string(item, f"{context} entry {index}")
        for index, item in enumerate(value)
    ]


# Validate one model manifest object
def validate_manifest_model(model_name: str, value: Any) -> dict[str, Any]:
    model = require_manifest_object(value, f"Model {model_name}")
    require_manifest_keys(
        model,
        {"import"},
        {
            "alpha_qed",
            "complex_mass_scheme",
            "particle_pdgs",
        },
        f"Model {model_name}",
    )
    model_import = require_manifest_string(
        require_manifest_field(model, "import", f"Model {model_name}"),
        f"Model {model_name} import",
    )
    if MODEL_IMPORT.fullmatch(model_import) is None:
        raise RuntimeError(f"Model {model_name} has an unsafe import")
    complex_mass_scheme = model.get("complex_mass_scheme", False)
    if not isinstance(complex_mass_scheme, bool):
        raise RuntimeError(f"Model {model_name} has a non-boolean complex_mass_scheme")

    particle_pdgs = require_manifest_object(
        model.get("particle_pdgs", {}), f"Model {model_name} particle_pdgs"
    )
    for particle, pdg in particle_pdgs.items():
        if PARTICLE_NAME.fullmatch(particle) is None or type(pdg) is not int or pdg == 0:
            raise RuntimeError(f"Model {model_name} has invalid particle_pdgs")

    if "alpha_qed" not in model:
        return model
    alpha = require_manifest_object(model["alpha_qed"], f"Model {model_name} alpha_qed")
    require_manifest_keys(
        alpha,
        {"charge", "charge_square"},
        {"inverse_alpha_zero"},
        f"Model {model_name} alpha_qed",
    )
    for field in ("charge", "charge_square"):
        require_manifest_identifier(
            require_manifest_field(alpha, field, f"Model {model_name} alpha_qed"),
            f"Model {model_name} alpha_qed {field}",
        )
    inverse = alpha.get("inverse_alpha_zero", 137.03599908)
    if (
        not isinstance(inverse, (int, float))
        or isinstance(inverse, bool)
        or not math.isfinite(inverse)
        or inverse <= 0.0
    ):
        raise RuntimeError(f"Model {model_name} has invalid inverse_alpha_zero")
    return model


# Validate explicit SLHA mass replacements keyed by signed PDG id
def validate_mass_overrides(value: Any, context: str) -> dict[str, Any]:
    masses = require_manifest_object(value, context)
    for pdg, mass in masses.items():
        if not isinstance(pdg, str) or re.fullmatch(r"-?[1-9][0-9]*", pdg) is None:
            raise RuntimeError(f"{context} has an invalid PDG id")
        if (
            not isinstance(mass, (int, float))
            or isinstance(mass, bool)
            or not math.isfinite(mass)
            or mass < 0.0
        ):
            raise RuntimeError(f"{context} has an invalid mass")
    return masses


# Validate one MG5 process line without restricting its generated topology
def validate_manifest_process_syntax(value: Any, context: str) -> str:
    process = require_manifest_string(value, context)
    try:
        process_registry.parse_process_syntax(process)
    except RuntimeError as error:
        raise RuntimeError(f"{context} has invalid MG5 process syntax") from error
    return process


# Validate one standalone process manifest object
def validate_manifest_process(value: Any, models: dict[str, Any]) -> dict[str, Any]:
    entry = require_manifest_object(value, "Standalone process")
    require_manifest_keys(
        entry,
        {"model", "name", "process", "projection"},
        {"kernel", "type"},
        "Standalone process",
    )
    name = require_manifest_identifier(
        require_manifest_field(entry, "name", "Standalone process"),
        "Standalone process name",
    )
    model = require_manifest_identifier(
        require_manifest_field(entry, "model", f"Process {name}"),
        f"Process {name} model",
    )
    if model not in models:
        raise RuntimeError(f"Process {name} references unknown model {model}")
    process = validate_manifest_process_syntax(
        require_manifest_field(entry, "process", f"Process {name}"),
        f"Process {name} process",
    )
    projection = require_manifest_string(
        require_manifest_field(entry, "projection", f"Process {name}"),
        f"Process {name} projection",
    )
    if projection not in process_registry.STANDALONE_PROJECTIONS:
        raise RuntimeError(f"Process {name} has invalid standalone projection {projection}")
    parsed = process_registry.parse_process_syntax(process)
    expected_incoming = (
        ("g", "g") if projection == process_registry.DURHAM_PROJECTION else ("a", "a")
    )
    if parsed.incoming != expected_incoming:
        raise RuntimeError(
            f"Process {name} projection {projection} requires incoming "
            f"{' '.join(expected_incoming)}"
        )
    if "kernel" in entry:
        require_manifest_identifier(entry["kernel"], f"Process {name} kernel")
        if projection != process_registry.PHOTON_PROJECTION:
            raise RuntimeError(f"Process {name} kernel requires photon projection")
    process_type = entry.get("type", "tree")
    if process_type not in {"tree", "loop"}:
        raise RuntimeError(f"Process {name} has invalid type {process_type}")
    if process_type == "loop":
        if projection != process_registry.DURHAM_PROJECTION:
            raise RuntimeError(f"Loop process {name} currently requires Durham projection")
        if "kernel" in entry:
            raise RuntimeError(f"Loop process {name} cannot select a tree kernel")
        if "noborn" not in process.lower():
            raise RuntimeError(f"Loop process {name} must use an explicit noborn selector")
    return entry


# Compute the explicit MG2GRA matrix-element implementation type
def process_type(entry: dict[str, Any]) -> str:
    return str(entry.get("type", "tree"))


# Compute only C++ standalone tree processes
def tree_processes(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    return [entry for entry in manifest["processes"] if process_type(entry) == "tree"]


# Compute only generated MadLoop processes
def loop_processes(manifest: dict[str, Any]) -> list[dict[str, Any]]:
    return [entry for entry in manifest["processes"] if process_type(entry) == "loop"]


# Validate one generated family manifest object
def validate_manifest_family(value: Any, models: dict[str, Any]) -> dict[str, Any]:
    family = require_manifest_object(value, "Generated family")
    require_manifest_keys(
        family,
        {"channel", "definitions", "model", "name", "processes", "projection"},
        {
            "mass_overrides",
            "process_names",
        },
        "Generated family",
    )
    name = require_manifest_identifier(
        require_manifest_field(family, "name", "Generated family"),
        "Generated family name",
    )
    if name in process_registry.STANDALONE_PROCESS_FAMILIES:
        raise RuntimeError(f"Family name {name} is a reserved standalone process family")
    if not name.startswith(process_registry.PROCESS_FAMILY_PREFIX):
        raise RuntimeError(
            f"Family name {name} must use the {process_registry.PROCESS_FAMILY_PREFIX} namespace"
        )
    model = require_manifest_identifier(
        require_manifest_field(family, "model", f"Family {name}"),
        f"Family {name} model",
    )
    if model not in models:
        raise RuntimeError(f"Family {name} references unknown model {model}")

    definitions = require_manifest_string_list(
        require_manifest_field(family, "definitions", f"Family {name}"),
        f"Family {name} definitions",
        allow_empty=True,
    )
    if any(MG5_DEFINITION.fullmatch(definition) is None for definition in definitions):
        raise RuntimeError(f"Family {name} has invalid MG5 definitions")
    if len(definitions) != len(set(definitions)):
        raise RuntimeError(f"Family {name} has duplicate MG5 definitions")

    processes = require_manifest_string_list(
        require_manifest_field(family, "processes", f"Family {name}"),
        f"Family {name} processes",
        allow_empty=False,
    )
    for index, process in enumerate(processes):
        validate_manifest_process_syntax(process, f"Family {name} process {index}")
    if len(processes) != len(set(processes)):
        raise RuntimeError(f"Family {name} has duplicate generated processes")
    if "process_names" in family:
        process_names = require_manifest_string_list(
            family["process_names"], f"Family {name} process_names", allow_empty=False
        )
        if len(process_names) != len(processes):
            raise RuntimeError(f"Family {name} process_names must match processes")
        for index, process_name in enumerate(process_names):
            require_manifest_identifier(process_name, f"Family {name} process_names entry {index}")
        if len(process_names) != len(set(process_names)):
            raise RuntimeError(f"Family {name} has duplicate process_names")

    projection = require_manifest_string(
        require_manifest_field(family, "projection", f"Family {name}"),
        f"Family {name} projection",
    )
    if projection not in process_registry.FAMILY_PROJECTIONS:
        raise RuntimeError(f"Family {name} has invalid family projection {projection}")
    expected_incoming = (
        ("a", "a") if projection == process_registry.PHOTON_PROJECTION else ("p", "p")
    )
    for index, process in enumerate(processes):
        if process_registry.parse_process_syntax(process).incoming != expected_incoming:
            raise RuntimeError(
                f"Family {name} process {index} projection {projection} requires "
                f"incoming {' '.join(expected_incoming)}"
            )
    if "mass_overrides" in family:
        validate_mass_overrides(family["mass_overrides"], f"Family {name} mass_overrides")
    require_manifest_identifier(
        require_manifest_field(family, "channel", f"Family {name}"),
        f"Family {name} channel",
    )
    process_registry.family_wrapper_name(family)
    return family


# Strictly validate one in-memory declarative process registry
def validate_manifest(value: Any) -> dict[str, Any]:
    manifest = require_manifest_object(value, "MG5 registry")
    require_manifest_keys(
        manifest,
        {"models", "processes"},
        {"families"},
        "MG5 registry",
    )
    models = require_manifest_object(
        require_manifest_field(manifest, "models", "MG5 registry"),
        "MG5 registry models",
    )
    if not models:
        raise RuntimeError("MG5 registry models must be nonempty")
    for model_name, model in models.items():
        require_manifest_identifier(model_name, "Model name")
        validate_manifest_model(model_name, model)

    processes_value = require_manifest_field(manifest, "processes", "MG5 registry")
    if not isinstance(processes_value, list):
        raise RuntimeError("MG5 registry processes must be a list")
    processes = [validate_manifest_process(entry, models) for entry in processes_value]
    process_names = [entry["name"] for entry in processes]
    if len(process_names) != len(set(process_names)):
        raise RuntimeError("MG5 registry contains duplicate process names")
    for name in process_names:
        if name in RESERVED_OUTPUT_NAMES:
            raise RuntimeError(f"Standalone process name {name} is reserved by a fixed MG5 output")

    families_value = manifest.get("families", [])
    if not isinstance(families_value, list):
        raise RuntimeError("MG5 registry families must be a list")
    families = [validate_manifest_family(family, models) for family in families_value]
    family_names = [family["name"] for family in families]
    if len(family_names) != len(set(family_names)):
        raise RuntimeError("MG5 registry contains duplicate family names")
    amplitudes = [process_registry.family_wrapper_name(family) for family in families]
    if len(amplitudes) != len(set(amplitudes)):
        raise RuntimeError("MG5 registry contains duplicate amplitude classes")
    channels = [(family["projection"], family["channel"]) for family in families]
    if len(channels) != len(set(channels)):
        raise RuntimeError("MG5 registry contains duplicate public family channels")
    return manifest


# Load and strictly validate the declarative process registry
def load_manifest(path: Path) -> dict[str, Any]:
    return validate_manifest(load_json_file(path))


# Compute one unused sibling path carrying the requested suffix
def unused_sibling(path: Path, suffix: str) -> Path:
    candidate = path.with_name(path.name + suffix)
    index = 1
    while candidate.exists():
        candidate = path.with_name(path.name + f"{suffix}.{index}")
        index += 1
    return candidate


# Preserve one existing file before replacing it with regenerated content
def backup_existing(path: Path) -> Path | None:
    if not path.exists():
        return None
    candidate = path.with_name(path.name + "._old")
    index = 1
    while candidate.exists():
        candidate = path.with_name(path.name + f"._old.{index}")
        index += 1
    path.rename(candidate)
    return candidate


# Move one file with a copy fallback when source and destination filesystems differ
def move_file(source: Path, destination: Path) -> None:
    try:
        source.replace(destination)
    except OSError as error:
        if error.errno != errno.EXDEV:
            raise
        shutil.copy2(source, destination)
        source.unlink()


# Record generated installations and restore them together after a failure
class InstallTransaction:
    # Construct an empty replacement journal
    def __init__(self, journal_root: Path | None = None) -> None:
        self._changes: list[tuple[Path, Path | None]] = []
        self._paths: set[Path] = set()
        self._retained: dict[Path, Path] = {}
        journal_root = journal_root or ROOT / "tmp"
        ensure_dir(journal_root)
        self._journal = Path(tempfile.mkdtemp(prefix="MG2GRA_transaction_", dir=journal_root))

    # Move one active output into the private transaction journal
    def _backup(self, path: Path, label: str) -> Path | None:
        if not path.exists():
            return None
        backup = self._journal / f"{len(self._changes):06d}_{label}_{path.name}"
        move_file(path, backup)
        return backup

    # Preserve one destination once before its first transaction write
    def preserve(self, path: Path) -> None:
        path = path.resolve()
        if path in self._paths:
            return
        self._paths.add(path)
        self._changes.append((path, self._backup(path, "old")))

    # Preserve one file while leaving a working copy for in-place transforms
    def preserve_in_place(self, path: Path) -> None:
        path = path.resolve()
        if path in self._paths:
            return
        self._paths.add(path)
        backup = self._backup(path, "old")
        self._changes.append((path, backup))
        if backup is not None:
            shutil.copy2(backup, path)

    # Retain every replaced output beside its active path after successful validation
    def commit(self) -> None:
        for path, backup in self._changes:
            if backup is None:
                continue
            retained = unused_sibling(path, "._old")
            try:
                move_file(backup, retained)
            finally:
                if retained.exists() and not backup.exists():
                    self._retained[path] = retained
        self._journal.rmdir()
        self._changes.clear()
        self._paths.clear()
        self._retained.clear()

    # Restore every replaced destination without deleting failed outputs
    def rollback(self) -> None:
        for path, backup in reversed(self._changes):
            if path.exists():
                failed = self._journal / f"failed_{len(self._changes):06d}_{path.name}"
                move_file(path, unused_sibling(failed, ""))
            retained = self._retained.get(path, backup)
            if retained is not None and retained.exists():
                move_file(retained, path)
        self._changes.clear()
        self._paths.clear()
        self._retained.clear()


# Preserve one active generated output during registry removal
def retire_output(path: Path, transaction: InstallTransaction) -> None:
    if path.is_file() and "._old" not in path.name:
        transaction.preserve(path)


# Preserve every active file owned by one generated directory
def retire_directory_outputs(path: Path, transaction: InstallTransaction) -> None:
    if not path.is_dir():
        return
    for output in sorted(path.iterdir()):
        retire_output(output, transaction)


# Preserve every active file below one generated directory tree
def retire_tree_outputs(path: Path, transaction: InstallTransaction) -> None:
    if not path.is_dir():
        return
    for output in sorted(item for item in path.rglob("*") if item.is_file()):
        retire_output(output, transaction)


# Preserve every obsolete active top-level MG5 C++ output as one retained path
def retire_legacy_cpp_layout(transaction: InstallTransaction) -> None:
    for relative_root in (output_layout.CPP_INCLUDE_ROOT, output_layout.CPP_SOURCE_ROOT):
        root = ROOT / relative_root
        if not root.is_dir():
            continue
        for path in sorted(root.iterdir()):
            if path.name in {
                output_layout.RUNTIME,
                output_layout.DURHAM,
                output_layout.PHOTON,
                output_layout.PARTON,
            }:
                continue
            if "._old" not in path.name:
                transaction.preserve(path)


# Find generated family wrappers and every family source included by each file
def installed_family_wrappers(directory: Path, suffix: str) -> dict[Path, frozenset[str]]:
    wrappers: dict[Path, frozenset[str]] = {}
    for path in sorted(directory.glob(f"AMP_MG5_*.{suffix}")):
        if "._old" in path.name:
            continue
        families = frozenset(FAMILY_WRAPPER_INCLUDE.findall(path.read_text()))
        if families:
            wrappers[path] = families
    return wrappers


# Compute the exact wrapper stem and family ownership declared by the manifest
def expected_family_wrappers(manifest: dict[str, Any]) -> dict[str, str]:
    return {
        process_registry.family_wrapper_name(family): family["name"]
        for family in manifest.get("families", [])
    }


# Validate exact generated family wrapper stems and header-source pairing
def validate_family_wrappers(
    manifest: dict[str, Any], include_root: Path, source_root: Path
) -> None:
    expected = expected_family_wrappers(manifest)
    source_files = installed_family_wrappers(source_root, "cc")
    header_files = installed_family_wrappers(include_root, "h")
    if any(len(owners) != 1 for owners in (*source_files.values(), *header_files.values())):
        raise RuntimeError("One MG5 family amplitude references multiple families")
    source_wrappers = {path.stem: next(iter(owners)) for path, owners in source_files.items()}
    header_wrappers = {path.stem: next(iter(owners)) for path, owners in header_files.items()}
    if source_wrappers != expected or header_wrappers != expected:
        raise RuntimeError("Orphaned or missing MG5 family amplitudes")


# Retire only outputs exclusively owned by removed registry entries
def retire_removed_outputs(
    original: dict[str, Any],
    current: dict[str, Any],
    transaction: InstallTransaction,
) -> None:
    current_processes = {entry["name"] for entry in current["processes"]}
    removed_processes = [
        entry for entry in original["processes"] if entry["name"] not in current_processes
    ]
    for entry in removed_processes:
        name = entry["name"]
        directory = standalone_cpp_directory(entry)
        for path in (
            cpp_header_path(directory, f"AMP_MG5_{name}.h"),
            cpp_source_path(directory, f"AMP_MG5_{name}.cc"),
        ):
            retire_output(path, transaction)
        retire_tree_outputs(
            CARDS_ROOT / output_layout.cards_directory(entry["projection"], name), transaction
        )

    current_families = {family["name"] for family in current.get("families", [])}
    for family in original.get("families", []):
        if family["name"] in current_families:
            continue
        family_name = family["name"]
        directory = family_cpp_directory(family)
        for path in (
            ROOT / output_layout.CPP_INCLUDE_ROOT / directory,
            ROOT / output_layout.CPP_SOURCE_ROOT / directory,
        ):
            if path.is_dir():
                transaction.preserve(path)
        amplitude = process_registry.family_wrapper_name(family)
        projection_directory = output_layout.projection_cpp_directory(family["projection"])
        retire_output(cpp_header_path(projection_directory, f"{amplitude}.h"), transaction)
        retire_output(cpp_source_path(projection_directory, f"{amplitude}.cc"), transaction)
        retire_directory_outputs(
            CARDS_ROOT / output_layout.cards_directory(family["projection"], family_name),
            transaction,
        )


# Retire installed outputs which are no longer present in the active manifest
def retire_orphaned_outputs(manifest: dict[str, Any], transaction: InstallTransaction) -> None:
    for projection in process_registry.STANDALONE_PROJECTIONS:
        directory = output_layout.projection_cpp_directory(projection)
        expected = {
            entry["name"] for entry in manifest["processes"] if entry["projection"] == projection
        }
        include_root = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
        source_root = ROOT / output_layout.CPP_SOURCE_ROOT / directory
        for source in sorted(source_root.glob("AMP_MG5_*.cc")):
            if (
                "._old" in source.name
                or "MadGraph to GRANIITTI conversion done" not in source.read_text()
            ):
                continue
            name = source.stem.removeprefix("AMP_MG5_")
            if name not in expected:
                retire_output(source, transaction)
                retire_output(include_root / f"AMP_MG5_{name}.h", transaction)

    expected_wrappers = expected_family_wrappers(manifest)
    for projection in process_registry.FAMILY_PROJECTIONS:
        directory = output_layout.projection_cpp_directory(projection)
        include_root = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
        source_root = ROOT / output_layout.CPP_SOURCE_ROOT / directory
        family_names = {
            family["name"]
            for family in manifest.get("families", [])
            if family["projection"] == projection
        }
        for root in (include_root, source_root):
            for family_path in sorted(root.glob("MG5_*")):
                if (
                    family_path.is_dir()
                    and "._old" not in family_path.name
                    and family_path.name not in family_names
                ):
                    transaction.preserve(family_path)
        wrappers = {
            **installed_family_wrappers(source_root, "cc"),
            **installed_family_wrappers(include_root, "h"),
        }
        for path, owners in wrappers.items():
            expected_family = expected_wrappers.get(path.stem)
            if expected_family is None or owners != frozenset({expected_family}):
                retire_output(path, transaction)

    active_owners = {
        projection: {
            entry["name"]
            for entry in (*manifest["processes"], *manifest.get("families", []))
            if entry["projection"] == projection
        }
        for projection in process_registry.PROJECTIONS
    }
    for projection in process_registry.PROJECTIONS:
        projection_root = CARDS_ROOT / output_layout.projection_directory(projection)
        for owner in sorted(projection_root.iterdir() if projection_root.is_dir() else ()):
            if (
                owner.is_dir()
                and "._old" not in owner.name
                and owner.name not in active_owners[projection]
            ):
                retire_tree_outputs(owner, transaction)

    runtime_root = CARDS_ROOT / output_layout.RUNTIME
    expected_runtime = {output_layout.PROCESS_REGISTRY, output_layout.CONVERTER_STATE}
    for path in sorted(runtime_root.iterdir() if runtime_root.is_dir() else ()):
        if path.is_file() and "._old" not in path.name and path.name not in expected_runtime:
            retire_output(path, transaction)

    # Retire the replaced flat layout during the next complete regeneration
    for path in sorted(CARDS_ROOT.glob("param_card_*.dat")):
        retire_output(path, transaction)
    retire_output(CARDS_ROOT / output_layout.PROCESS_REGISTRY, transaction)
    legacy_registry = CARDS_ROOT / "registry"
    if legacy_registry.is_dir():
        transaction.preserve(legacy_registry)


# Atomically replace one generated text destination while retaining its file mode
def write_atomic(
    path: Path,
    content: str,
    failed_root: Path | None = None,
    output_mode: int | None = None,
) -> None:
    if output_mode is None:
        output_mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "w") as output:
            os.fchmod(output.fileno(), output_mode)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        move_file(temporary, path)
    except Exception:
        if temporary.exists():
            failed = (failed_root or path.parent) / f"failed_write_{path.name}"
            move_file(temporary, unused_sibling(failed, ""))
        raise


# Atomically replace one generated binary destination
def write_atomic_bytes(
    path: Path,
    content: bytes,
    failed_root: Path | None = None,
    output_mode: int | None = None,
) -> None:
    if output_mode is None:
        output_mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else 0o644
    descriptor, temporary_name = tempfile.mkstemp(prefix=f".{path.name}.", dir=path.parent)
    temporary = Path(temporary_name)
    try:
        with os.fdopen(descriptor, "wb") as output:
            os.fchmod(output.fileno(), output_mode)
            output.write(content)
            output.flush()
            os.fsync(output.fileno())
        move_file(temporary, path)
    except Exception:
        if temporary.exists():
            failed = (failed_root or path.parent) / f"failed_write_{path.name}"
            move_file(temporary, unused_sibling(failed, ""))
        raise


# Install text without touching an already identical output
def install_text(path: Path, content: str, transaction: InstallTransaction | None = None) -> bool:
    content = standalone.normalize_whitespace(content)
    if path.exists() and path.read_text() == content:
        return False
    ensure_dir(path.parent)
    output_mode = stat.S_IMODE(path.stat().st_mode) if path.exists() else None
    if transaction is None:
        backup_existing(path)
    else:
        transaction.preserve(path)
    write_atomic(
        path,
        content,
        transaction._journal if transaction is not None else None,
        output_mode,
    )
    return True


# Install normalized generated text without touching an already identical output
def install_file(
    source: Path, destination: Path, transaction: InstallTransaction | None = None
) -> bool:
    content = standalone.normalize_whitespace(source.read_text())
    if destination.exists() and destination.read_text() == content:
        return False
    ensure_dir(destination.parent)
    output_mode = stat.S_IMODE(destination.stat().st_mode) if destination.exists() else None
    if transaction is None:
        backup_existing(destination)
    else:
        transaction.preserve(destination)
    write_atomic(
        destination,
        content,
        transaction._journal if transaction is not None else None,
        output_mode,
    )
    return True


# Install one exact generated binary file without touching identical output
def install_binary(
    source: Path, destination: Path, transaction: InstallTransaction | None = None
) -> bool:
    content = source.read_bytes()
    if destination.exists() and destination.read_bytes() == content:
        return False
    ensure_dir(destination.parent)
    output_mode = stat.S_IMODE(destination.stat().st_mode) if destination.exists() else None
    if transaction is None:
        backup_existing(destination)
    else:
        transaction.preserve(destination)
    write_atomic_bytes(
        destination,
        content,
        transaction._journal if transaction is not None else None,
        output_mode,
    )
    return True


# Install a generated default card without overwriting a checked in process card
def install_default_card(
    source: Path, destination: Path, transaction: InstallTransaction | None = None
) -> bool:
    if destination.exists():
        return False
    return install_file(source, destination, transaction)


# Compute requested masses parsed structurally from the SLHA MASS block
def slha_masses(content: str, pdgs: set[int]) -> dict[int, float]:
    masses: dict[int, float] = {}
    in_mass = False
    for line in standalone.normalize_whitespace(content).splitlines():
        stripped = line.split("#", 1)[0].strip()
        if re.match(r"(?i)^block\s+mass(?:\s|$)", stripped):
            in_mass = True
            continue
        if in_mass and re.match(r"(?i)^(?:block|decay)\s", stripped):
            in_mass = False
        if not in_mass or not stripped:
            continue
        fields = stripped.split()
        if len(fields) < 2 or re.fullmatch(r"-?[0-9]+", fields[0]) is None:
            continue
        pdg = int(fields[0])
        if pdg not in pdgs:
            continue
        if pdg in masses:
            raise RuntimeError(f"SLHA MASS block contains duplicate PDG id {pdg}")
        try:
            mass = float(fields[1].replace("d", "e").replace("D", "E"))
        except ValueError as error:
            raise RuntimeError(f"SLHA MASS entry {pdg} has an invalid value") from error
        if not math.isfinite(mass):
            raise RuntimeError(f"SLHA MASS entry {pdg} has a non-finite value")
        masses[pdg] = mass
    missing = pdgs - set(masses)
    if missing:
        values = ", ".join(str(pdg) for pdg in sorted(missing))
        raise RuntimeError(f"SLHA MASS block is missing PDG ids {values}")
    return masses


# Apply explicit family mass replacements to one generated SLHA card
def set_slha_masses(content: str, overrides: dict[str, Any]) -> str:
    values = {int(pdg): float(mass) for pdg, mass in overrides.items()}
    slha_masses(content, set(values))
    lines = standalone.normalize_whitespace(content).splitlines()
    in_mass = False
    replaced: set[int] = set()
    for index, line in enumerate(lines):
        stripped = line.split("#", 1)[0].strip()
        if re.match(r"(?i)^block\s+mass(?:\s|$)", stripped):
            in_mass = True
            continue
        if in_mass and re.match(r"(?i)^(?:block|decay)\s", stripped):
            in_mass = False
        if not in_mass:
            continue
        match = re.match(r"^(\s*)(-?[0-9]+)(\s+)(\S+)(.*)$", line)
        if match is None:
            continue
        pdg = int(match.group(2))
        if pdg not in values:
            continue
        lines[index] = (
            f"{match.group(1)}{match.group(2)}{match.group(3)}{values[pdg]:.16e}{match.group(5)}"
        )
        replaced.add(pdg)
    if replaced != set(values):
        raise RuntimeError("Could not apply every SLHA MASS override")
    return "\n".join(lines) + "\n"


# Load the MASS entries owned by the previous installed family generation
def installed_family_mass_overrides(
    family: dict[str, Any], root: Path | None = None
) -> dict[str, Any]:
    root = ROOT if root is None else root
    family_name = family["name"]
    path = (
        root
        / output_layout.CARDS_ROOT
        / output_layout.channels_path(family["projection"], family_name)
    )
    if not path.is_file():
        return {}
    try:
        data = load_json_file(path)
    except json.JSONDecodeError as error:
        raise RuntimeError(f"Invalid generated family channels for {family_name}") from error
    if not isinstance(data, dict) or data.get("family") != family_name:
        raise RuntimeError(f"Generated family channel name mismatch for {family_name}")
    generation_setup = data.get("generation_setup")
    if not isinstance(generation_setup, dict):
        raise RuntimeError(f"Invalid generated family generation setup for {family_name}")
    return validate_mass_overrides(
        generation_setup.get("mass_overrides", {}),
        f"Previous family {family_name} mass_overrides",
    )


# Install one generated family card with its automatic scattering mass scheme
def install_family_card(
    source: Path,
    destination: Path,
    family: dict[str, Any],
    transaction: InstallTransaction | None = None,
    previous_overrides: dict[str, Any] | None = None,
) -> bool:
    overrides = family.get("mass_overrides", {})
    if destination.exists():
        previous = previous_overrides or {}
        removed = {int(pdg) for pdg in previous if pdg not in overrides}
        if not overrides and not removed:
            return False
        content = destination.read_text()
        replacements = {
            str(pdg): mass for pdg, mass in slha_masses(source.read_text(), removed).items()
        }
        replacements.update(overrides)
    else:
        content = source.read_text()
        replacements = overrides
    if replacements:
        content = set_slha_masses(content, replacements)
    return install_text(destination, content, transaction)


# Locate the requested MadGraph command-line frontend
def resolve_mg5(explicit: str | None) -> Path:
    candidate = explicit or os.environ.get("MG5AMC") or shutil.which("mg5_aMC")
    if candidate is None:
        raise RuntimeError("Pass --mg5 /path/to/bin/mg5_aMC or set MG5AMC")
    path = Path(candidate).expanduser().resolve()
    if not path.is_file():
        raise RuntimeError(f"MadGraph executable does not exist: {path}")
    return path


# Read the exact particle pole symbols from the selected restricted UFO model
def load_model_particles(mg5_root: Path, model: dict) -> list[dict]:
    mg5_root = mg5_root.resolve()
    work_dir = Path(tempfile.mkdtemp(prefix="MG5_model_", dir=ROOT / "tmp"))
    with mg5_family.chdir(work_dir):
        command_class, _ = mg5_family.load_madgraph_family_api(mg5_root)
        command = command_class()
        command.no_notification()
        command.exec_cmd("set complex_mass_scheme " +
                         ("True --allow_qed" if model.get("complex_mass_scheme", False) else "False"),
                         printcmd=False)
        command.exec_cmd(f"import model {model['import']}", printcmd=False)
        particles = command._curr_model.get("particles")
        model["particle_pdgs"] = {
            name: sign * int(particle.get("pdg_code"))
            for particle in particles
            for name, sign in ((particle.get("name"), 1), (particle.get("antiname"), -1))
            if sign == 1 or not particle.get("self_antipart")
        }
        return [{"pdg": abs(int(particle.get("pdg_code"))),
                 "mass": particle.get("mass"), "width": particle.get("width"),
                 "signed_mass": particle.is_fermion() and particle.get("self_antipart")}
                for particle in particles]


# Group registered processes by their UFO model
def model_processes(manifest: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    grouped = {name: [] for name in manifest["models"]}
    for process in manifest["processes"]:
        grouped[process["model"]].append(process)
    return grouped


# Build one MadGraph batch file containing individual and aggregate exports
def build_mg5_commands(
    manifest: dict[str, Any], work_dir: Path
) -> tuple[str, dict[str, Path], dict[str, Path], dict[str, Path]]:
    lines: list[str] = []
    process_exports: dict[str, Path] = {}
    support_exports: dict[str, Path] = {}
    family_exports: dict[str, Path] = {}
    grouped = model_processes(manifest)
    for model_name, model in manifest["models"].items():
        entries = grouped[model_name]
        families = [
            family for family in manifest.get("families", []) if family["model"] == model_name
        ]
        if not entries and not families:
            continue
        if model.get("complex_mass_scheme", False):
            lines.append("set complex_mass_scheme True --allow_qed")
        else:
            lines.append("set complex_mass_scheme False")
        lines.append(f"import model {model['import']}")
        tree_entries = [entry for entry in entries if process_type(entry) == "tree"]
        loop_entries = [entry for entry in entries if process_type(entry) == "loop"]
        for entry in tree_entries:
            output = work_dir / f"process_{entry['name']}"
            process_exports[entry["name"]] = output
            lines.extend(
                (
                    f"generate {entry['process']}",
                    f"output standalone_cpp {output}",
                )
            )
        if tree_entries:
            support = work_dir / f"support_{model_name}"
            support_exports[model_name] = support
            lines.append(f"generate {tree_entries[0]['process']}")
            lines.extend(f"add process {entry['process']}" for entry in tree_entries[1:])
            lines.append(f"output standalone_cpp {support}")
        for entry in loop_entries:
            output = work_dir / f"process_{entry['name']}"
            process_exports[entry["name"]] = output
            lines.extend(
                (
                    f"generate {entry['process']}",
                    f"output {output}",
                )
            )
        for family in families:
            lines.extend(family.get("definitions", []))
            family_processes = family["processes"]
            lines.append(f"generate {family_processes[0]}")
            lines.extend(f"add process {process}" for process in family_processes[1:])
            output = work_dir / f"family_{family['name']}"
            family_exports[family["name"]] = output
            lines.append(f"output standalone_cpp {output}")
    return (
        "\n".join(lines) + "\n",
        process_exports,
        support_exports,
        family_exports,
    )


# Compute the sole subprocess source directory of an individual export
def individual_subprocess(export: Path) -> Path:
    sources = sorted(export.glob("SubProcesses/*/CPPProcess.cc"))
    if len(sources) != 1:
        raise RuntimeError(f"Expected one subprocess in {export}, found {len(sources)}")
    return sources[0].parent


# Reject stale generated sources and return their MG5 release
def validate_generator_stamp(path: Path) -> str:
    header = path.read_text()[:1000]
    match = re.search(r"MadGraph5_aMC@NLO v\.\s*([^,\n]+),\s*(\d{4})-", header)
    if match is None:
        raise RuntimeError(f"No MadGraph version stamp in {path}")
    version = match.group(1).strip()
    if version not in SUPPORTED_MG5_VERSIONS:
        supported = ", ".join(sorted(SUPPORTED_MG5_VERSIONS))
        raise RuntimeError(
            f"Unsupported MadGraph version {version} in {path}; supported versions: {supported}"
        )
    return version


# Compute every installed file owned by one generated family
def family_installed_source_paths(family: dict[str, Any], root: Path) -> tuple[Path, ...]:
    family_name = family["name"]
    directory = family_cpp_directory(family)
    directories = (
        root / output_layout.CPP_INCLUDE_ROOT / directory,
        root / output_layout.CPP_SOURCE_ROOT / directory,
    )
    paths = sorted(
        (
            path
            for directory in directories
            if directory.is_dir()
            for path in directory.iterdir()
            if path.is_file() and path.suffix in {".h", ".cc"} and "._old" not in path.name
        ),
        key=lambda path: path.relative_to(root).as_posix(),
    )
    if not paths:
        raise RuntimeError(f"Family {family_name} has no installed generated sources")
    return tuple(paths)


# Hash exact relative paths and bytes for one installed generated family
def family_installed_source_digest(family: dict[str, Any], root: Path) -> str:
    digest = hashlib.sha256()
    for path in family_installed_source_paths(family, root):
        relative = path.relative_to(root).as_posix().encode("utf-8")
        content = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


# Hash one exact file while retaining its repository-relative ownership
def managed_file_hashes(paths: list[Path], root: Path) -> dict[str, str]:
    hashes: dict[str, str] = {}
    for path in sorted(paths, key=lambda item: item.relative_to(root).as_posix()):
        if not path.is_file():
            raise RuntimeError(f"Missing generated MG5 output {path}")
        relative = path.relative_to(root).as_posix()
        hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return hashes


# Hash every converter source which can change installed generated C++
def converter_source_digest() -> str:
    paths = [BASE_DIR / "regenerate.py", BASE_DIR / "run.sh"]
    paths.extend(sorted((BASE_DIR / "modules").glob("*.py")))
    digest = hashlib.sha256()
    for path in paths:
        relative = path.relative_to(BASE_DIR).as_posix().encode("utf-8")
        content = path.read_bytes()
        digest.update(len(relative).to_bytes(8, "big"))
        digest.update(relative)
        digest.update(len(content).to_bytes(8, "big"))
        digest.update(content)
    return digest.hexdigest()


# Discover installed global support ownership from standalone process pairs
def installed_model_support(manifest: dict[str, Any], root: Path) -> dict[str, ModelFiles]:
    support: dict[str, ModelFiles] = {}
    suffix_owners: dict[str, str] = {}
    for entry in tree_processes(manifest):
        name = entry["name"]
        directory = standalone_cpp_directory(entry)
        header_path = cpp_header_path(directory, f"AMP_MG5_{name}.h", root)
        source_path = cpp_source_path(directory, f"AMP_MG5_{name}.cc", root)
        if not header_path.is_file() or not source_path.is_file():
            raise RuntimeError(f"Incomplete installed standalone process {name}")
        files = discover_process_model(header_path.read_text(), source_path.read_text())
        previous = support.setdefault(entry["model"], files)
        if previous != files:
            raise RuntimeError(f"Model support mismatch for {entry['model']}")
        owner = suffix_owners.setdefault(files.model_suffix, entry["model"])
        if owner != entry["model"]:
            raise RuntimeError(
                f"Standalone models {owner} and {entry['model']} share installed "
                f"support suffix {files.model_suffix}"
            )
    return support


# Compute every generated process sidecar which participates in runtime selection
def converter_sidecar_paths(manifest: dict[str, Any], root: Path) -> list[Path]:
    cards_root = root / output_layout.CARDS_ROOT
    paths = [
        cards_root / output_layout.color_path(entry["projection"], entry["name"])
        for entry in manifest["processes"]
    ]
    paths.extend(
        cards_root / output_layout.channels_path(family["projection"], family["name"])
        for family in manifest.get("families", [])
    )
    paths.extend(
        cards_root
        / output_layout.cards_directory(entry["projection"], entry["name"])
        / "source"
        / madloop.SOURCE_DATA
        for entry in loop_processes(manifest)
    )
    return paths


# Compute every generated photon registry and family wrapper output
def photon_registry_paths(manifest: dict[str, Any], root: Path) -> list[Path]:
    paths = [
        cpp_header_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.h", root),
        cpp_source_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.cc", root),
    ]
    for family in photon_registry.photon_families(manifest):
        amplitude = process_registry.family_wrapper_name(family)
        paths.extend(
            [
                cpp_header_path(output_layout.PHOTON, f"{amplitude}.h", root),
                cpp_source_path(output_layout.PHOTON, f"{amplitude}.cc", root),
            ]
        )
    return paths


# Build exact converter and managed standalone source provenance
def converter_state_data(manifest: dict[str, Any], root: Path) -> dict[str, Any]:
    model_files = installed_model_support(manifest, root)
    standalones = []
    for entry in manifest["processes"]:
        name = entry["name"]
        directory = standalone_cpp_directory(entry)
        paths = [
            cpp_header_path(directory, f"AMP_MG5_{name}.h", root),
            cpp_source_path(directory, f"AMP_MG5_{name}.cc", root),
        ]
        standalones.append(
            {"name": name, "model": entry["model"], "files": managed_file_hashes(paths, root)}
        )

    models = []
    for model_name, files in model_files.items():
        directory = output_layout.model_cpp_directory(model_name)
        paths = [
            cpp_header_path(directory, files.helas_header, root),
            cpp_source_path(directory, files.helas_header.replace(".h", ".cc"), root),
            cpp_header_path(directory, files.parameter_header, root),
            cpp_source_path(directory, files.parameter_header.replace(".h", ".cc"), root),
        ]
        models.append(
            {
                "model": model_name,
                "model_suffix": files.model_suffix,
                "files": managed_file_hashes(paths, root),
            }
        )

    slha = managed_file_hashes(
        [cpp_header_path(output_layout.RUNTIME, "read_slha.h", root),
         cpp_source_path(output_layout.RUNTIME, "read_slha.cc", root)], root)
    return {
        "schema": CONVERTER_STATE_SCHEMA,
        "converter_source_digest": converter_source_digest(),
        "manifest_digest": hashlib.sha256(dump_manifest(manifest).encode("utf-8")).hexdigest(),
        "standalones": standalones,
        "models": models,
        "read_slha": slha,
        "sidecars": managed_file_hashes(converter_sidecar_paths(manifest, root), root),
        "photon": managed_file_hashes(photon_registry_paths(manifest, root), root),
    }


# Install exact converter provenance after all managed sources are current
def install_converter_state(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    install_text(
        CONVERTER_STATE_PATH,
        dump_manifest(converter_state_data(manifest, ROOT)),
        transaction,
    )


# Validate converter freshness and every managed standalone support byte
def validate_converter_state(manifest: dict[str, Any]) -> None:
    if not CONVERTER_STATE_PATH.is_file():
        raise RuntimeError("Missing generated MG5 converter state")
    try:
        installed = load_json_file(CONVERTER_STATE_PATH)
    except json.JSONDecodeError as error:
        raise RuntimeError("Invalid generated MG5 converter state") from error
    expected = converter_state_data(manifest, ROOT)
    if installed != expected:
        raise RuntimeError("Stale generated MG5 converter state or managed support")

    expected_support = {
        ROOT / relative for model in expected["models"] for relative in model["files"]
    }
    installed_support: set[Path] = set()
    for root in (
        ROOT / output_layout.CPP_INCLUDE_ROOT / output_layout.RUNTIME / "Models",
        ROOT / output_layout.CPP_SOURCE_ROOT / output_layout.RUNTIME / "Models",
    ):
        for pattern in ("HelAmps_*.h", "HelAmps_*.cc", "Parameters_*.h", "Parameters_*.cc"):
            installed_support.update(
                path for path in root.rglob(pattern) if "._old" not in path.name
            )
    if installed_support != expected_support:
        raise RuntimeError("Orphaned or missing global MG5 model support")


# Derive the common MG5 release and installed source digest for one family
def family_mg5_source(family: dict[str, Any], root: Path) -> process_registry.MG5Source:
    family_name = family["name"]
    source_directory = root / output_layout.CPP_SOURCE_ROOT / family_cpp_directory(family)
    sources = [
        path
        for path in family_installed_source_paths(family, root)
        if path.parent == source_directory and path.suffix == ".cc"
    ]
    versions = {validate_generator_stamp(path) for path in sources}
    if len(versions) != 1:
        raise RuntimeError(f"Family {family_name} has inconsistent MG5 version stamps")
    return process_registry.MG5Source(
        mg5_version=versions.pop(),
        source_digest=family_installed_source_digest(family, root),
    )


# Validate the stored MG5 source against the installed generated files
def validate_mg5_source(
    family: dict[str, Any], family_channels: mg5_family.FamilyChannels, root: Path
) -> None:
    expected = family_mg5_source(family, root)
    if family_channels.mg5_source != expected:
        raise RuntimeError(f"Stale MG5 source for {family_channels.family}")


# Discover the model support file set in one aggregate export
def aggregate_model_files(export: Path) -> ModelFiles:
    headers = sorted((export / "src").glob("Parameters_*.h"))
    helas = sorted((export / "src").glob("HelAmps_*.h"))
    if len(headers) != 1 or len(helas) != 1:
        raise RuntimeError(f"Expected one model support set in {export}")
    suffix = headers[0].stem.removeprefix("Parameters_")
    if helas[0].stem != f"HelAmps_{suffix}":
        raise RuntimeError(f"Inconsistent aggregate model support in {export}")
    return ModelFiles(
        parameter_class=headers[0].stem,
        parameter_header=headers[0].name,
        helas_header=helas[0].name,
        model_suffix=suffix,
    )


# Reject logical standalone models which would install the same global C++ names
def validate_support_owners(support_exports: dict[str, Path]) -> dict[str, ModelFiles]:
    files_by_model: dict[str, ModelFiles] = {}
    owners: dict[str, str] = {}
    for model_name, export in support_exports.items():
        files = aggregate_model_files(export)
        previous = owners.setdefault(files.model_suffix, model_name)
        if previous != model_name:
            raise RuntimeError(
                f"Standalone models {previous} and {model_name} both generate "
                f"global support suffix {files.model_suffix}"
            )
        files_by_model[model_name] = files
    return files_by_model


# Convert the common SLHA reader from one standalone or family export
def slha_support(export: Path) -> tuple[str, str]:
    source_dir = export / "src"
    return (
        transform_slha_header((source_dir / "read_slha.h").read_text()),
        transform_slha_source((source_dir / "read_slha.cc").read_text()),
    )


# Verify all active exports provide one compatible common SLHA reader
def validate_slha_exports(
    support_exports: dict[str, Path], family_exports: dict[str, Path]
) -> Path | None:
    exports = [*support_exports.items(), *family_exports.items()]
    installed: tuple[str, str] | None = None
    selected: Path | None = None
    for owner, export in exports:
        candidate = slha_support(export)
        if installed is not None and candidate != installed:
            raise RuntimeError(f"Generated output {owner} has incompatible read_slha support")
        if installed is None:
            installed = candidate
            selected = export
    return selected


# Compute global model support paths produced by active standalone models
def model_support_paths(support_exports: dict[str, Path]) -> set[Path]:
    paths: set[Path] = set()
    for model_name, files in validate_support_owners(support_exports).items():
        directory = output_layout.model_cpp_directory(model_name)
        paths.update(
            {
                cpp_header_path(directory, files.helas_header),
                cpp_source_path(directory, files.helas_header.replace(".h", ".cc")),
                cpp_header_path(directory, files.parameter_header),
                cpp_source_path(directory, files.parameter_header.replace(".h", ".cc")),
            }
        )
    return paths


# Retire global model support no longer owned by a standalone logical model
def retire_orphaned_model_support(
    support_exports: dict[str, Path], transaction: InstallTransaction
) -> None:
    active = model_support_paths(support_exports)
    roots = (
        ROOT / output_layout.CPP_INCLUDE_ROOT / output_layout.RUNTIME / "Models",
        ROOT / output_layout.CPP_SOURCE_ROOT / output_layout.RUNTIME / "Models",
    )
    for root in roots:
        for pattern in ("HelAmps_*.h", "HelAmps_*.cc", "Parameters_*.h", "Parameters_*.cc"):
            for path in sorted(root.rglob(pattern)):
                if "._old" not in path.name and path not in active:
                    retire_output(path, transaction)


# Transform and install one aggregate model support set
def install_model_support(
    model_name: str,
    model: dict[str, Any],
    export: Path,
    process_sources: list[str],
    particles: list[dict],
    transaction: InstallTransaction | None = None,
) -> tuple[ModelFiles, str, str]:
    files = aggregate_model_files(export)
    source_dir = export / "src"
    raw_helas_header = (source_dir / files.helas_header).read_text()
    raw_helas_source = (source_dir / files.helas_header.replace(".h", ".cc")).read_text()
    raw_parameter_header = (source_dir / files.parameter_header).read_text()
    raw_parameter_source = (source_dir / files.parameter_header.replace(".h", ".cc")).read_text()

    validate_dependency_union(
        process_sources,
        raw_helas_header,
        raw_parameter_header,
        files.model_suffix,
    )

    alpha = model.get("alpha_qed")
    charge = alpha.get("charge") if alpha else None
    charge_square = alpha.get("charge_square") if alpha else None
    inverse_alpha = float(alpha.get("inverse_alpha_zero", 137.03599908)) if alpha else 0.0
    alpha_zero = alpha is not None and supports_alpha_qed_zero(
        raw_parameter_source, files.parameter_class, charge, charge_square
    )
    helas_header = transform_helas_header(raw_helas_header)
    directory = output_layout.model_cpp_directory(model_name)
    helas_source = transform_helas_source(raw_helas_source, files.helas_header, directory)
    parameter_header = add_model_particles(transform_parameters_header(
        raw_parameter_header, files.parameter_class, alpha_zero=alpha_zero
    ), particles, raw_parameter_source, charge)
    parameter_source = transform_parameters_source(
        raw_parameter_source,
        files.parameter_class,
        files.parameter_header,
        charge if alpha_zero else None,
        charge_square if alpha_zero else None,
        inverse_alpha,
        directory,
    )

    install_text(cpp_header_path(directory, files.helas_header), helas_header, transaction)
    install_text(
        cpp_source_path(directory, files.helas_header.replace(".h", ".cc")),
        helas_source,
        transaction,
    )
    install_text(cpp_header_path(directory, files.parameter_header), parameter_header, transaction)
    install_text(
        cpp_source_path(directory, files.parameter_header.replace(".h", ".cc")),
        parameter_source,
        transaction,
    )
    return files, helas_header, parameter_header


# Transform and install the generated common SLHA reader
def install_slha_reader(
    export: Path | None, transaction: InstallTransaction | None = None
) -> tuple[str, str]:
    header, source = slha_support(export) if export is not None else (transform_slha_header(), transform_slha_source())
    install_text(cpp_header_path(output_layout.RUNTIME, "read_slha.h"), header, transaction)
    install_text(cpp_source_path(output_layout.RUNTIME, "read_slha.cc"), source, transaction)
    return header, source


# Transform and install standalone exports directly from the temporary MG5 run
def install_standalone_processes(
    manifest: dict[str, Any],
    exports: dict[str, Path],
    transaction: InstallTransaction | None = None,
) -> dict[str, list[str]]:
    sources_by_model: dict[str, list[str]] = {name: [] for name in manifest["models"]}
    for entry in tree_processes(manifest):
        subprocess_dir = individual_subprocess(exports[entry["name"]])
        validate_generator_stamp(subprocess_dir / "CPPProcess.cc")
        raw_header = (subprocess_dir / "CPPProcess.h").read_text()
        raw_source = (subprocess_dir / "CPPProcess.cc").read_text()
        color = standalone.extract_color_structure(raw_source)
        files = discover_process_model(raw_header, raw_source)
        alpha = manifest["models"][entry["model"]].get("alpha_qed", {})
        alpha_zero = supports_alpha_qed_zero(
            (exports[entry["name"]] / "src" / files.parameter_header.replace(".h", ".cc")).read_text(),
            files.parameter_class, alpha.get("charge"), alpha.get("charge_square"),
        )
        header = standalone.transform_header(
            entry["name"], raw_header, color, entry["projection"], entry["model"]
        )
        source = standalone.transform_source(
            entry["name"], raw_source, color, entry["projection"], entry["model"], alpha_zero
        )
        directory = standalone_cpp_directory(entry)
        install_text(
            cpp_header_path(directory, f"AMP_MG5_{entry['name']}.h"),
            header,
            transaction,
        )
        install_text(
            cpp_source_path(directory, f"AMP_MG5_{entry['name']}.cc"),
            source,
            transaction,
        )
        install_file(
            exports[entry["name"]] / "Cards/param_card.dat",
            parameter_card_path(entry),
            transaction,
        )
        sources_by_model[entry["model"]].append(raw_source)
    return sources_by_model


# Install one converted MadLoop runtime source bundle and Durham wrapper
def install_loop_processes(
    manifest: dict[str, Any],
    mg5: Path,
    exports: dict[str, Path],
    work_dir: Path,
    transaction: InstallTransaction | None = None,
) -> None:
    mg5_root = mg5.parent.parent
    for entry in loop_processes(manifest):
        export = exports[entry["name"]]
        loop_model = madloop.loop_model(mg5_root, export, entry, manifest["models"][entry["model"]])
        metadata = madloop.inspect_output(export, entry, loop_model)
        staging = work_dir / f"bundle_{entry['name']}"
        madloop.create_bundle(export, staging, mg5_root, metadata)
        owner = CARDS_ROOT / output_layout.cards_directory(entry["projection"], entry["name"])
        source = owner / "source"
        for filename in madloop.SOURCE_ARCHIVES.values():
            install_binary(staging / filename, source / filename, transaction)
        install_file(staging / madloop.SOURCE_DATA, source / madloop.SOURCE_DATA, transaction)
        install_binary(
            staging / madloop.SOURCE_LICENSE,
            source / madloop.SOURCE_LICENSE,
            transaction,
        )
        install_file(export / "Cards" / "param_card.dat", parameter_card_path(entry), transaction)
        process_data = madloop.durham_color_data(entry, metadata)
        install_text(
            standalone_color_path(entry),
            json.dumps(process_data, indent=2) + "\n",
            transaction,
        )
        directory = standalone_cpp_directory(entry)
        install_text(
            cpp_header_path(directory, f"AMP_MG5_{entry['name']}.h"),
            madloop.durham_header(entry, metadata),
            transaction,
        )
        install_text(
            cpp_source_path(directory, f"AMP_MG5_{entry['name']}.cc"),
            madloop.durham_source(entry, metadata),
            transaction,
        )


# Generate exact finite-Nc Durham process data
def install_durham_process_data(
    manifest: dict[str, Any],
    mg5: Path,
    exports: dict[str, Path],
    work_dir: Path,
    transaction: InstallTransaction | None = None,
) -> None:
    mg5_root = mg5.parent.parent
    if not (mg5_root / "madgraph").is_dir():
        raise RuntimeError(f"Could not locate the MadGraph Python package under {mg5_root}")

    for entry in tree_processes(manifest):
        if entry["projection"] != process_registry.DURHAM_PROJECTION:
            continue
        model_import = manifest["models"][entry["model"]]["import"]
        process_data = mg5_color.generate_durham_data(
            mg5_root,
            model_import,
            entry["process"],
            work_dir / "mg5_color_durham",
        )
        subprocess_dir = individual_subprocess(exports[entry["name"]])
        raw_source = (subprocess_dir / "CPPProcess.cc").read_text()
        color = standalone.extract_color_structure(raw_source)
        if color is None or int(color["ncolor"]) != int(process_data["ncolor"]):
            raise RuntimeError(f"MadGraph color-basis mismatch for Durham process {entry['name']}")
        install_text(
            standalone_color_path(entry),
            json.dumps(process_data, indent=2) + "\n",
            transaction,
        )


# Generate exact finite-Nc photon process data
def install_photon_process_data(
    manifest: dict[str, Any],
    mg5: Path,
    exports: dict[str, Path],
    work_dir: Path,
    transaction: InstallTransaction | None = None,
) -> None:
    mg5_root = mg5.parent.parent
    if not (mg5_root / "madgraph").is_dir():
        raise RuntimeError(f"Could not locate the MadGraph Python package under {mg5_root}")

    for entry in tree_processes(manifest):
        if entry["projection"] != process_registry.PHOTON_PROJECTION:
            continue
        subprocess_dir = individual_subprocess(exports[entry["name"]])
        raw_source = (subprocess_dir / "CPPProcess.cc").read_text()
        color = standalone.extract_color_structure(raw_source)
        if color is None:
            raise RuntimeError(
                f"MadGraph generated no color metric for photon process {entry['name']}"
            )
        model_import = manifest["models"][entry["model"]]["import"]
        process_data = mg5_color.generate_photon_data(
            mg5_root,
            model_import,
            entry["process"],
            raw_source,
            (exports[entry["name"]] / "Cards/param_card.dat").read_text(),
            color,
            work_dir / "photon_metadata",
            manifest["models"][entry["model"]].get("alpha_qed", {}).get("charge"),
            manifest["models"][entry["model"]].get("complex_mass_scheme", False),
        )
        if int(color["ncolor"]) != int(process_data["ncolor"]):
            raise RuntimeError(f"MadGraph color-basis mismatch for photon process {entry['name']}")
        install_text(
            standalone_color_path(entry),
            json.dumps(process_data, indent=2) + "\n",
            transaction,
        )


# Generate and install the Durham process registry
def install_durham_registry(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    header, source = durham_registry.generate_registry(CARDS_ROOT, manifest)
    install_text(
        cpp_header_path(output_layout.DURHAM, "AMP_MG5_DurhamRegistry.h"), header, transaction
    )
    install_text(
        cpp_source_path(output_layout.DURHAM, "AMP_MG5_DurhamRegistry.cc"), source, transaction
    )


# Generate and install all photon process amplitudes and registries
def install_photon_registry(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    for legacy in (
        ROOT / output_layout.CPP_INCLUDE_ROOT / "AMP_MG5_PhotonFamilyRegistry.h",
        ROOT / output_layout.CPP_SOURCE_ROOT / "AMP_MG5_PhotonFamilyRegistry.cc",
    ):
        if transaction is None:
            backup_existing(legacy)
        else:
            retire_output(legacy, transaction)
    header, source = photon_registry.generate_registry(CARDS_ROOT, manifest)
    install_text(
        cpp_header_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.h"), header, transaction
    )
    install_text(
        cpp_source_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.cc"), source, transaction
    )
    for family in photon_registry.photon_families(manifest):
        amplitude = process_registry.family_wrapper_name(family)
        install_text(
            cpp_header_path(output_layout.PHOTON, f"{amplitude}.h"),
            photon_registry.family_amplitude_header(family),
            transaction,
        )
        install_text(
            cpp_source_path(output_layout.PHOTON, f"{amplitude}.cc"),
            photon_registry.family_amplitude_source(family),
            transaction,
        )


# Generate and install the hard process registry
def install_parton_registry(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    for family in parton_registry.parton_families(manifest):
        amplitude = process_registry.family_wrapper_name(family)
        install_text(
            cpp_header_path(output_layout.PARTON, f"{amplitude}.h"),
            parton_registry.amplitude_header(family),
            transaction,
        )
        install_text(
            cpp_source_path(output_layout.PARTON, f"{amplitude}.cc"),
            parton_registry.amplitude_source(family),
            transaction,
        )
    header, source = parton_registry.generate_registry(manifest)
    install_text(
        cpp_header_path(output_layout.PARTON, "AMP_MG5_PartonRegistry.h"), header, transaction
    )
    install_text(
        cpp_source_path(output_layout.PARTON, "AMP_MG5_PartonRegistry.cc"), source, transaction
    )


# Generate and install the common amplitude process registry
def install_process_registry(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    registry, header, source = process_registry.generate_registry(CARDS_ROOT, manifest)
    install_text(CARDS_ROOT / output_layout.registry_path(), registry, transaction)
    install_text(
        cpp_header_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h"), header, transaction
    )
    install_text(
        cpp_source_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.cc"), source, transaction
    )


# Preserve generated family files which are no longer present in a fresh export
def backup_stale_family_files(
    family: dict[str, Any],
    desired_headers: set[str],
    desired_sources: set[str],
    transaction: InstallTransaction | None = None,
) -> None:
    directory = family_cpp_directory(family)
    include_dir = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
    source_dir = ROOT / output_layout.CPP_SOURCE_ROOT / directory
    for path in include_dir.glob("*"):
        if path.is_file() and "._old" not in path.name and path.name not in desired_headers:
            if transaction is None:
                backup_existing(path)
            else:
                transaction.preserve(path)
    for path in source_dir.glob("*"):
        if path.is_file() and "._old" not in path.name and path.name not in desired_sources:
            if transaction is None:
                backup_existing(path)
            else:
                transaction.preserve(path)


# Transform and install one complete multi-subprocess MadGraph family
def install_family(
    family: dict[str, Any],
    model: dict[str, Any],
    export: Path,
    color_structure: tuple[mg5_family.SubprocessColorStructure, ...],
    particles: list[dict],
    transaction: InstallTransaction | None = None,
) -> None:
    family_name = family["name"]
    previous_mass_overrides = installed_family_mass_overrides(family)
    files = aggregate_model_files(export)
    source_dir = export / "src"
    raw_helas_header = (source_dir / files.helas_header).read_text()
    raw_helas_source = (source_dir / files.helas_header.replace(".h", ".cc")).read_text()
    raw_parameter_header = (source_dir / files.parameter_header).read_text()
    raw_parameter_source = (source_dir / files.parameter_header.replace(".h", ".cc")).read_text()

    raw_processes: list[tuple[str, str, str]] = []
    for raw_source_path in sorted(export.glob("SubProcesses/*/CPPProcess.cc")):
        raw_header_path = raw_source_path.with_name("CPPProcess.h")
        validate_generator_stamp(raw_source_path)
        raw_processes.append(
            (
                raw_source_path.parent.name,
                raw_header_path.read_text(),
                raw_source_path.read_text(),
            )
        )
    if not raw_processes:
        raise RuntimeError(f"No generated subprocesses found for {family_name}")
    validate_dependency_union(
        [source for _, _, source in raw_processes],
        raw_helas_header,
        raw_parameter_header,
        files.model_suffix,
    )

    alpha = model.get("alpha_qed")
    charge = alpha.get("charge") if alpha else None
    charge_square = alpha.get("charge_square") if alpha else None
    inverse_alpha = float(alpha.get("inverse_alpha_zero", 137.03599908)) if alpha else 0.0
    alpha_zero = alpha is not None and supports_alpha_qed_zero(
        raw_parameter_source, files.parameter_class, charge, charge_square
    )
    indexed_colors = mg5_family.indexed_subprocess_color_structure(color_structure)
    directory = family_cpp_directory(family)
    raw_names = {generated_name for generated_name, _, _ in raw_processes}
    if raw_names != set(indexed_colors):
        raise RuntimeError(
            f"MadGraph API and standalone subprocess exports disagree for {family_name}"
        )
    transformed_rows = []
    for generated_name, raw_header, raw_source in raw_processes:
        process = mg5_family.transform_subprocess(
            raw_header,
            raw_source,
            family_name,
            generated_name,
            files.parameter_header,
            files.helas_header,
            files.parameter_class,
            files.model_suffix,
            alpha_zero,
            directory,
            indexed_colors[generated_name],
            model.get("particle_pdgs"),
        )
        transformed_rows.append(
            (
                generated_name,
                mg5_family.update_subprocess(process, alpha_zero),
            )
        )
    transformed = [process for _, process in transformed_rows]
    desired_headers = {
        files.helas_header,
        files.parameter_header,
        "ProcessBase.h",
        "Processes.h",
        *(f"{process.class_name}.h" for process in transformed),
    }
    desired_sources = {
        files.helas_header.replace(".h", ".cc"),
        files.parameter_header.replace(".h", ".cc"),
        *(f"{process.class_name}.cc" for process in transformed),
    }
    backup_stale_family_files(family, desired_headers, desired_sources, transaction)

    include_dir = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
    source_out = ROOT / output_layout.CPP_SOURCE_ROOT / directory
    install_text(
        include_dir / files.helas_header,
        mg5_family.transform_family_helas_header(
            raw_helas_header, family_name, files.helas_header, files.model_suffix
        ),
        transaction,
    )
    install_text(
        source_out / files.helas_header.replace(".h", ".cc"),
        mg5_family.transform_family_helas_source(
            raw_helas_source,
            family_name,
            files.helas_header,
            files.model_suffix,
            directory,
        ),
        transaction,
    )
    install_text(
        include_dir / files.parameter_header,
        add_model_particles(mg5_family.transform_family_parameters_header(
            raw_parameter_header,
            family_name,
            files.parameter_class,
            files.parameter_header,
            alpha_zero,
        ), particles, raw_parameter_source, charge),
        transaction,
    )
    install_text(
        source_out / files.parameter_header.replace(".h", ".cc"),
        mg5_family.transform_family_parameters_source(
            raw_parameter_source,
            family_name,
            files.parameter_class,
            files.parameter_header,
            charge if alpha_zero else None,
            charge_square if alpha_zero else None,
            inverse_alpha,
            directory,
        ),
        transaction,
    )
    install_text(
        include_dir / "ProcessBase.h", mg5_family.process_base_header(family_name), transaction
    )
    install_text(
        include_dir / "Processes.h",
        mg5_family.processes_header(family_name, transformed, directory),
        transaction,
    )
    for process in transformed:
        install_text(include_dir / f"{process.class_name}.h", process.header, transaction)
        install_text(source_out / f"{process.class_name}.cc", process.source, transaction)

    generation_setup = process_registry.family_generation_setup_from_manifest(family, model)
    mg5_source = family_mg5_source(family, ROOT)
    family_channels = mg5_family.make_family_channels(
        family_name, generation_setup, mg5_source, transformed_rows
    )
    install_text(
        family_channels_path(family),
        dump_manifest(mg5_family.family_channels_data(family_channels)),
        transaction,
    )

    install_family_card(
        export / "Cards/param_card.dat",
        parameter_card_path(family),
        family,
        transaction,
        previous_mass_overrides,
    )


# Install registries from the generated standalone and family metadata
def install_registries(
    manifest: dict[str, Any], transaction: InstallTransaction | None = None
) -> None:
    install_durham_registry(manifest, transaction)
    install_photon_registry(manifest, transaction)
    install_parton_registry(manifest, transaction)
    install_process_registry(manifest, transaction)
    for family in manifest.get("families", []):
        path = family_channels_path(family)
        data = load_json_file(path)
        data["mg5_source"] = process_registry.mg5_source_data(family_mg5_source(family, ROOT))
        install_text(path, dump_manifest(data), transaction)
    install_converter_state(manifest, transaction)


# Remove C++ comments before structural source validation
def strip_cpp_comments(source: str) -> str:
    return re.sub(r"//[^\n]*|/\*.*?\*/", "", source, flags=re.DOTALL)


# Compute one comment-free C++ class declaration and body
def cpp_class_definition(source: str, class_name: str) -> tuple[str, str] | None:
    declaration = re.search(
        rf"\bclass\s+{re.escape(class_name)}\b(?P<declaration>[^;{{]*){{",
        source,
        flags=re.DOTALL,
    )
    if declaration is None:
        return None
    depth = 1
    index = declaration.end()
    while index < len(source) and depth > 0:
        if source[index] == "{":
            depth += 1
        elif source[index] == "}":
            depth -= 1
        index += 1
    if depth != 0:
        return None
    return declaration.group("declaration"), source[declaration.end() : index - 1]


# Validate one advertised family amplitude
def validate_family_amplitude(family: dict[str, Any], root: Path) -> None:
    amplitude = process_registry.family_wrapper_name(family)
    directory = output_layout.projection_cpp_directory(family["projection"])
    header = cpp_header_path(directory, f"{amplitude}.h", root)
    source = cpp_source_path(directory, f"{amplitude}.cc", root)
    if not header.is_file() or not source.is_file():
        raise RuntimeError(f"Family {family['name']} amplitude {amplitude} is not installed")
    projection = process_registry.family_projection(family)
    if projection == process_registry.PARTON_PROJECTION:
        expected_header = parton_registry.amplitude_header(family)
        expected_source = parton_registry.amplitude_source(family)
    elif projection == process_registry.PHOTON_PROJECTION:
        expected_header = photon_registry.family_amplitude_header(family)
        expected_source = photon_registry.family_amplitude_source(family)
    else:
        raise RuntimeError(f"Family {family['name']} has no generated amplitude type")
    if header.read_text() != standalone.normalize_whitespace(expected_header):
        raise RuntimeError(f"Family {family['name']} amplitude header is stale")
    if source.read_text() != standalone.normalize_whitespace(expected_source):
        raise RuntimeError(f"Family {family['name']} amplitude source is stale")
    header_text = strip_cpp_comments(header.read_text())
    definition = cpp_class_definition(header_text, amplitude)
    hard_family = projection == process_registry.PARTON_PROJECTION
    hard_inheritance = re.compile(r"(?:^|[:,])\s*public\s+(?:gra::)?PartonMG5Process\s*(?:,|$)")
    if hard_family and (definition is None or hard_inheritance.search(definition[0]) is None):
        raise RuntimeError(f"Family {family['name']} amplitude does not inherit PartonMG5Process")
    photon_inheritance = re.compile(r"(?:^|[:,])\s*public\s+(?:gra::)?PhotonMG5Process\s*(?:,|$)")
    if not hard_family and (definition is None or photon_inheritance.search(definition[0]) is None):
        raise RuntimeError(f"Family {family['name']} amplitude does not inherit PhotonMG5Process")
    subprocess_sum = re.compile(
        rf"\b(?:gra::)?mg5::SubprocessSum\s*<\s*"
        rf"{re.escape(family['name'])}::ProcessBase\s*>\s+"
        r"[A-Za-z_][A-Za-z0-9_]*\s*;"
    )
    if subprocess_sum.search(definition[1]) is None:
        raise RuntimeError(f"Family {family['name']} amplitude has no subprocess sum")
    source_text = strip_cpp_comments(source.read_text())
    channel_builder = re.compile(rf"\b{re.escape(family['name'])}::BuildSubprocesses\s*\(")
    if channel_builder.search(source_text) is None:
        raise RuntimeError(
            f"Family {family['name']} amplitude does not build generated subprocesses"
        )


# Validate exact subprocess channels against installed sources
def validate_family_channels(family: dict[str, Any], model: dict[str, Any], root: Path) -> None:
    family_name = family["name"]
    directory = family_cpp_directory(family)
    persistent = load_family_channels(family, model)
    validate_mg5_source(family, persistent, root)
    header_path = cpp_header_path(directory, "Processes.h", root)
    header = header_path.read_text()
    class_names = re.findall(r"make_unique<([^>]+)>", header)
    if not class_names:
        raise RuntimeError(f"Family {family_name} has no generated process builders")
    expected_classes = [
        f"{family_name}_{subprocess.generated_name}" for subprocess in persistent.subprocesses
    ]
    if class_names != expected_classes:
        raise RuntimeError(f"Stale generated subprocess order for {family_name}")
    for class_name, subprocess_channels in zip(class_names, persistent.subprocesses, strict=True):
        source_path = cpp_source_path(directory, f"{class_name}.cc", root)
        if not source_path.is_file():
            raise RuntimeError(f"Family {family_name} is missing source {class_name}.cc")
        source_channels = mg5_family.channels(source_path.read_text(), model.get("particle_pdgs"))
        source_structure = [
            (channel.initial, channel.final, channel.topology) for channel in source_channels
        ]
        persistent_structure = [
            (channel.initial, channel.final, channel.topology)
            for channel in subprocess_channels.channels
        ]
        if source_structure != persistent_structure:
            raise RuntimeError(f"Source and persistent channel topology disagree for {class_name}")
    expected = mg5_family.subprocess_channels_block(
        [subprocess.channels for subprocess in persistent.subprocesses]
    )
    marker = "inline std::vector<std::vector<gra::mg5::Channel>> SubprocessChannels()"
    if marker not in header:
        raise RuntimeError(f"Family {family_name} has no generated subprocess channels")
    start = header.index(marker)
    builder_marker = "// Build generated subprocesses with their exact channels"
    end_marker = (
        f"\n\n{builder_marker}"
        if builder_marker in header[start:]
        else f"\n\n}}  // namespace {family_name}"
    )
    end = header.index(end_marker, start)
    if header[start:end] != expected:
        raise RuntimeError(f"Stale generated subprocess channels for {family_name}")
    if builder_marker not in header:
        raise RuntimeError(f"Stale generated subprocess construction for {family_name}")


# Validate complete model-support pairs and dependency coverage for one family
def validate_family_support(
    family_name: str, include_dir: Path, source_dir: Path, sources: list[Path]
) -> None:
    parameter_headers = sorted(include_dir.glob("Parameters_*.h"))
    parameter_sources = sorted(source_dir.glob("Parameters_*.cc"))
    helas_headers = sorted(include_dir.glob("HelAmps_*.h"))
    helas_sources = sorted(source_dir.glob("HelAmps_*.cc"))
    if not (
        len(parameter_headers)
        == len(parameter_sources)
        == len(helas_headers)
        == len(helas_sources)
        == 1
    ):
        raise RuntimeError(f"Incomplete model support for MG5 family {family_name}")
    suffix = parameter_headers[0].stem.removeprefix("Parameters_")
    if (
        parameter_sources[0].stem != f"Parameters_{suffix}"
        or helas_headers[0].stem != f"HelAmps_{suffix}"
        or helas_sources[0].stem != f"HelAmps_{suffix}"
    ):
        raise RuntimeError(f"Mismatched model support for MG5 family {family_name}")
    validate_dependency_union(
        [source.read_text() for source in sources],
        helas_headers[0].read_text(),
        parameter_headers[0].read_text(),
        suffix,
    )


# Validate the deterministic common SLHA reader and its source pair
def validate_slha_reader(manifest: dict[str, Any]) -> None:
    header_path = cpp_header_path(output_layout.RUNTIME, "read_slha.h")
    source_path = cpp_source_path(output_layout.RUNTIME, "read_slha.cc")
    expected_header = transform_slha_header()
    expected_source = transform_slha_source()
    if not header_path.is_file() or header_path.read_text() != expected_header:
        raise RuntimeError("Missing or stale common MG5 read_slha header")
    if not source_path.is_file() or source_path.read_text() != expected_source:
        raise RuntimeError("Missing or stale common MG5 read_slha source")


# Validate the exact active file ownership of converter-managed trees
def validate_managed_outputs(manifest: dict[str, Any]) -> None:
    active_directories = {
        output_layout.RUNTIME,
        output_layout.DURHAM,
        output_layout.PHOTON,
        output_layout.PARTON,
    }
    for relative_root in (output_layout.CPP_INCLUDE_ROOT, output_layout.CPP_SOURCE_ROOT):
        root = ROOT / relative_root
        active = {path.name for path in root.iterdir() if "._old" not in path.name}
        if active != active_directories:
            raise RuntimeError("MG5 C++ roots contain misplaced active outputs")

    for projection in process_registry.STANDALONE_PROJECTIONS:
        directory = output_layout.projection_cpp_directory(projection)
        include_root = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
        source_root = ROOT / output_layout.CPP_SOURCE_ROOT / directory
        expected = {
            entry["name"] for entry in manifest["processes"] if entry["projection"] == projection
        }
        installed_sources = {
            path.stem.removeprefix("AMP_MG5_")
            for path in source_root.glob("AMP_MG5_*.cc")
            if "MadGraph to GRANIITTI conversion done" in path.read_text()
        }
        installed_headers = {
            path.stem.removeprefix("AMP_MG5_")
            for path in include_root.glob("AMP_MG5_*.h")
            if "MadGraph to GRANIITTI conversion done" in path.read_text()
        }
        if installed_sources != expected or installed_headers != expected:
            raise RuntimeError("Orphaned or missing standalone MG5 amplitudes")

    expected_owners = {
        projection: {
            entry["name"]
            for entry in (*manifest["processes"], *manifest.get("families", []))
            if entry["projection"] == projection
        }
        for projection in process_registry.PROJECTIONS
    }
    for projection, expected in expected_owners.items():
        projection_root = CARDS_ROOT / output_layout.projection_directory(projection)
        installed = {
            path.name
            for path in projection_root.iterdir()
            if path.is_dir()
            and "._old" not in path.name
            and any(item.is_file() and "._old" not in item.name for item in path.iterdir())
        }
        if installed != expected:
            raise RuntimeError(f"Orphaned or missing MG5 {projection} card owners")
    for entry in manifest["processes"]:
        owner = parameter_card_path(entry).parent
        installed = {path.name for path in owner.iterdir() if "._old" not in path.name}
        expected = {output_layout.PARAMETER_CARD, output_layout.COLOR_DATA}
        if process_type(entry) == "loop":
            expected.add("source")
        if installed != expected:
            raise RuntimeError("Orphaned or missing standalone MG5 card data")
        if process_type(entry) == "loop":
            source = owner / "source"
            installed_source = {path.name for path in source.iterdir() if "._old" not in path.name}
            expected_source = {
                *madloop.SOURCE_ARCHIVES.values(),
                madloop.SOURCE_DATA,
                madloop.SOURCE_LICENSE,
            }
            if installed_source != expected_source:
                raise RuntimeError("Orphaned or missing MadLoop source data")
    for family in manifest.get("families", []):
        owner = parameter_card_path(family).parent
        installed = {path.name for path in owner.iterdir() if "._old" not in path.name}
        if installed != {output_layout.PARAMETER_CARD, output_layout.CHANNEL_DATA}:
            raise RuntimeError("Orphaned or missing MG5 family card data")

    source_wrappers: dict[Path, frozenset[str]] = {}
    header_wrappers: dict[Path, frozenset[str]] = {}
    installed_families: set[str] = set()
    for projection in process_registry.FAMILY_PROJECTIONS:
        directory = output_layout.projection_cpp_directory(projection)
        include_root = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
        source_root = ROOT / output_layout.CPP_SOURCE_ROOT / directory
        for root in (include_root, source_root):
            installed_families.update(
                path.name for path in root.glob("MG5_*")
                if path.is_dir() and "._old" not in path.name
            )
        source_wrappers.update(installed_family_wrappers(source_root, "cc"))
        header_wrappers.update(installed_family_wrappers(include_root, "h"))
    family_names = {family["name"] for family in manifest.get("families", [])}
    if installed_families != family_names:
        raise RuntimeError("Orphaned or missing MG5 family source directories")
    expected = expected_family_wrappers(manifest)
    installed_sources = {path.stem: next(iter(owners)) for path, owners in source_wrappers.items()}
    installed_headers = {path.stem: next(iter(owners)) for path, owners in header_wrappers.items()}
    if installed_sources != expected or installed_headers != expected:
        raise RuntimeError("Orphaned or missing MG5 family amplitudes")


# Validate installed amplitudes, support dependencies, cards, and generator ages
def validate_installed(manifest: dict[str, Any]) -> None:
    validate_slha_reader(manifest)
    grouped_sources: dict[str, list[str]] = {name: [] for name in manifest["models"]}
    model_files: dict[str, ModelFiles] = {}
    for entry in manifest["processes"]:
        directory = standalone_cpp_directory(entry)
        header_path = cpp_header_path(directory, f"AMP_MG5_{entry['name']}.h")
        source_path = cpp_source_path(directory, f"AMP_MG5_{entry['name']}.cc")
        card_path = parameter_card_path(entry)
        if not header_path.is_file() or not source_path.is_file() or not card_path.is_file():
            raise RuntimeError(f"Incomplete installed standalone process {entry['name']}")
        validate_generator_stamp(source_path)
        header = header_path.read_text()
        source = source_path.read_text()
        if "MadGraph to GRANIITTI conversion done" not in header or (
            "MadGraph to GRANIITTI conversion done" not in source
        ):
            raise RuntimeError(f"Unconverted standalone process {entry['name']}")
        if process_type(entry) == "loop":
            source_root = card_path.parent / "source"
            data = madloop.load_source_data(source_root)
            metadata = madloop.MadLoopMetadata(**data["metadata"])
            if (
                metadata.name != entry["name"]
                or metadata.process != entry["process"]
                or metadata.model != entry["model"]
            ):
                raise RuntimeError(f"Stale MadLoop source data for {entry['name']}")
            process_data = durham_registry.load_process_data(CARDS_ROOT, entry)
            if int(process_data["ncolor"]) != metadata.ncolor:
                raise RuntimeError(f"Stale Durham loop color data for {entry['name']}")
            continue
        files = discover_process_model(header, source)
        previous = model_files.setdefault(entry["model"], files)
        if previous != files:
            raise RuntimeError(f"Model support mismatch for {entry['model']}")
        grouped_sources[entry["model"]].append(source)
        ncolor_match = re.search(r"static const int ncolor\s*=\s*(\d+);", header)
        if ncolor_match is None:
            raise RuntimeError(f"Missing generated color dimension for {entry['name']}")
        ncolor = int(ncolor_match.group(1))

        if entry["projection"] == process_registry.DURHAM_PROJECTION:
            process_data = durham_registry.load_process_data(CARDS_ROOT, entry)
            if int(process_data["ncolor"]) != ncolor:
                raise RuntimeError(f"Stale Durham color data for {entry['name']}")
        if entry["projection"] == process_registry.PHOTON_PROJECTION:
            process_data = photon_registry.load_process_data(CARDS_ROOT, entry)
            if int(process_data["ncolor"]) != ncolor:
                raise RuntimeError(f"Stale photon color data for {entry['name']}")

    expected_registry_header, expected_registry_source = durham_registry.generate_registry(
        CARDS_ROOT, manifest
    )
    expected_registry_header = standalone.normalize_whitespace(expected_registry_header)
    expected_registry_source = standalone.normalize_whitespace(expected_registry_source)
    registry_header = cpp_header_path(output_layout.DURHAM, "AMP_MG5_DurhamRegistry.h")
    registry_source = cpp_source_path(output_layout.DURHAM, "AMP_MG5_DurhamRegistry.cc")
    if not registry_header.is_file() or registry_header.read_text() != expected_registry_header:
        raise RuntimeError("Stale generated Durham registry header")
    if not registry_source.is_file() or registry_source.read_text() != expected_registry_source:
        raise RuntimeError("Stale generated Durham registry source")

    expected_photon_header, expected_photon_source = photon_registry.generate_registry(
        CARDS_ROOT, manifest
    )
    expected_photon_header = standalone.normalize_whitespace(expected_photon_header)
    expected_photon_source = standalone.normalize_whitespace(expected_photon_source)
    photon_header = cpp_header_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.h")
    photon_source = cpp_source_path(output_layout.PHOTON, "AMP_MG5_PhotonRegistry.cc")
    if not photon_header.is_file() or photon_header.read_text() != expected_photon_header:
        raise RuntimeError("Stale generated photon registry header")
    if not photon_source.is_file() or photon_source.read_text() != expected_photon_source:
        raise RuntimeError("Stale generated photon registry source")
    legacy_photon_registry = (
        ROOT / output_layout.CPP_INCLUDE_ROOT / "AMP_MG5_PhotonFamilyRegistry.h",
        ROOT / output_layout.CPP_SOURCE_ROOT / "AMP_MG5_PhotonFamilyRegistry.cc",
    )
    if any(path.is_file() for path in legacy_photon_registry):
        raise RuntimeError("Orphaned MG5 photon family registry")

    expected_process_json, expected_process_header, expected_process_source = (
        process_registry.generate_registry(CARDS_ROOT, manifest)
    )
    expected_process_json = standalone.normalize_whitespace(expected_process_json)
    expected_process_header = standalone.normalize_whitespace(expected_process_header)
    expected_process_source = standalone.normalize_whitespace(expected_process_source)
    process_json = CARDS_ROOT / output_layout.registry_path()
    process_header = cpp_header_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.h")
    process_source = cpp_source_path(output_layout.RUNTIME, "AMP_MG5_ProcessRegistry.cc")
    if not process_json.is_file() or process_json.read_text() != expected_process_json:
        raise RuntimeError("Stale generated MG5 process registry JSON")
    if not process_header.is_file() or process_header.read_text() != expected_process_header:
        raise RuntimeError("Stale generated MG5 process registry header")
    if not process_source.is_file() or process_source.read_text() != expected_process_source:
        raise RuntimeError("Stale generated MG5 process registry source")

    expected_hard_header, expected_hard_source = parton_registry.generate_registry(manifest)
    expected_hard_header = standalone.normalize_whitespace(expected_hard_header)
    expected_hard_source = standalone.normalize_whitespace(expected_hard_source)
    hard_header = cpp_header_path(output_layout.PARTON, "AMP_MG5_PartonRegistry.h")
    hard_source = cpp_source_path(output_layout.PARTON, "AMP_MG5_PartonRegistry.cc")
    if not hard_header.is_file() or hard_header.read_text() != expected_hard_header:
        raise RuntimeError("Stale generated MG5 hard process registry header")
    if not hard_source.is_file() or hard_source.read_text() != expected_hard_source:
        raise RuntimeError("Stale generated MG5 hard process registry source")

    for model_name, sources in grouped_sources.items():
        if not sources:
            continue
        files = model_files[model_name]
        directory = output_layout.model_cpp_directory(model_name)
        helas_header_path = cpp_header_path(directory, files.helas_header)
        helas_source_path = cpp_source_path(directory, files.helas_header.replace(".h", ".cc"))
        parameters_header_path = cpp_header_path(directory, files.parameter_header)
        parameters_source_path = cpp_source_path(
            directory, files.parameter_header.replace(".h", ".cc")
        )
        if not all(
            path.is_file()
            for path in (
                helas_header_path,
                helas_source_path,
                parameters_header_path,
                parameters_source_path,
            )
        ):
            raise RuntimeError(f"Incomplete installed model support for {model_name}")
        helas_header = helas_header_path.read_text()
        parameters_header = parameters_header_path.read_text()
        validate_dependency_union(sources, helas_header, parameters_header, files.model_suffix)

    for family in manifest.get("families", []):
        family_name = family["name"]
        directory = family_cpp_directory(family)
        include_dir = ROOT / output_layout.CPP_INCLUDE_ROOT / directory
        source_dir = ROOT / output_layout.CPP_SOURCE_ROOT / directory
        headers = sorted(include_dir.glob(f"{family_name}_P*.h"))
        sources = sorted(source_dir.glob(f"{family_name}_P*.cc"))
        if not headers or {header.stem for header in headers} != {
            source.stem for source in sources
        }:
            raise RuntimeError(f"Incomplete installed MG5 family {family_name}")
        if (
            not (include_dir / "ProcessBase.h").is_file()
            or not (include_dir / "Processes.h").is_file()
        ):
            raise RuntimeError(
                f"Missing generated family base or subprocess definitions for {family_name}"
            )
        for source in sources:
            validate_generator_stamp(source)
        validate_family_support(family_name, include_dir, source_dir, sources)
        card_path = parameter_card_path(family)
        if not card_path.is_file():
            raise RuntimeError(f"Missing parameter card for {family_name}")
        overrides = family.get("mass_overrides", {})
        if overrides:
            card_masses = slha_masses(
                card_path.read_text(),
                {int(pdg) for pdg in overrides},
            )
            for pdg, expected in overrides.items():
                if not math.isclose(
                    card_masses[int(pdg)], float(expected), rel_tol=1e-12, abs_tol=1e-15
                ):
                    raise RuntimeError(f"Stale MASS override for {family_name} PDG {pdg}")
        validate_family_channels(family, manifest["models"][family["model"]], ROOT)
        validate_family_amplitude(family, ROOT)
    validate_managed_outputs(manifest)
    validate_converter_state(manifest)


# Add one optional UFO model requested by a registry command
def add_registry_model(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    if args.model not in manifest["models"]:
        if args.model_import is None:
            raise RuntimeError("A new model requires --model-import")
        model: dict[str, Any] = {"import": args.model_import}
        if args.complex_mass_scheme:
            model["complex_mass_scheme"] = True
        if args.alpha_charge is not None:
            model["alpha_qed"] = {
                "charge": args.alpha_charge,
                "charge_square": args.alpha_charge_square,
                "inverse_alpha_zero": args.inverse_alpha_zero,
            }
        manifest["models"][args.model] = model
    elif args.complex_mass_scheme and not manifest["models"][args.model].get(
        "complex_mass_scheme", False
    ):
        raise RuntimeError(
            "--complex-mass-scheme requires a new model name or an existing complex-mass model"
        )


# Add one standalone process without hand-editing the registry
def add_standalone_registry_entry(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    if args.mg5_process is None:
        raise RuntimeError("--add requires --process")
    if any(entry["name"] == args.add for entry in manifest["processes"]):
        raise RuntimeError(f"Process {args.add} already exists")
    if args.projection not in process_registry.STANDALONE_PROJECTIONS:
        raise RuntimeError("--add requires --projection durham or photon")
    add_registry_model(args, manifest)
    entry = {
        "name": args.add,
        "model": args.model,
        "process": args.mg5_process,
        "projection": args.projection,
    }
    if getattr(args, "process_type", "tree") == "loop":
        entry["type"] = "loop"
    manifest["processes"].append(entry)


# Add one multi-subprocess family without hand-editing the registry
def add_family_registry_entry(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    if not args.family_process:
        raise RuntimeError("--add-family requires at least one --family-process")
    if any(entry["name"] == args.add_family for entry in manifest["families"]):
        raise RuntimeError(f"Family {args.add_family} already exists")
    add_registry_model(args, manifest)
    if args.projection not in process_registry.FAMILY_PROJECTIONS:
        raise RuntimeError("--add-family requires --projection photon or parton")
    if args.channel is None:
        raise RuntimeError("--add-family requires --channel")

    family: dict[str, Any] = {
        "name": args.add_family,
        "model": args.model,
        "definitions": args.definition,
        "processes": args.family_process,
        "projection": args.projection,
        "channel": args.channel,
    }
    if args.mass_override:
        masses = {pdg: mass for pdg, mass in args.mass_override}
        if len(masses) != len(args.mass_override):
            raise RuntimeError("--mass-override contains duplicate PDG ids")
        family["mass_overrides"] = masses
    manifest["families"].append(family)


# Remove one requested registry entry
def remove_registry_entry(args: argparse.Namespace, manifest: dict[str, Any]) -> None:
    if args.remove is not None:
        matches = [entry for entry in manifest["processes"] if entry["name"] == args.remove]
        if len(matches) != 1:
            raise RuntimeError(f"Unknown standalone process {args.remove}")
        manifest["processes"].remove(matches[0])
        return
    if args.remove_family is None:
        return
    matches = [
        family for family in manifest.get("families", []) if family["name"] == args.remove_family
    ]
    if len(matches) != 1:
        raise RuntimeError(f"Unknown generated family {args.remove_family}")
    removed = matches[0]
    manifest["families"].remove(removed)


# Add one requested registry entry only after validating a deep candidate
def add_registry_entry(args: argparse.Namespace, manifest: dict[str, Any]) -> dict[str, Any]:
    if args.add_family is None and (
        args.channel is not None or args.complex_mass_scheme or args.mass_override
    ):
        raise RuntimeError(
            "--channel, --complex-mass-scheme and --mass-override require --add-family"
        )
    changes = [
        args.add is not None,
        args.add_family is not None,
        args.remove is not None,
        args.remove_family is not None,
    ]
    if sum(changes) > 1:
        raise RuntimeError(
            "--add, --add-family, --remove and --remove-family are mutually exclusive"
        )
    if not any(changes):
        if args.projection is not None:
            raise RuntimeError("--projection requires --add or --add-family")
        return manifest
    candidate = copy.deepcopy(manifest)
    if args.add is not None:
        add_standalone_registry_entry(args, candidate)
    elif args.add_family is not None:
        add_family_registry_entry(args, candidate)
    else:
        remove_registry_entry(args, candidate)
    return validate_manifest(candidate)


# Parse one explicit PDG mass replacement
def parse_mass_override(value: str) -> tuple[str, float]:
    fields = value.split("=", 1)
    if len(fields) != 2:
        raise argparse.ArgumentTypeError("mass override must use PDG=MASS")
    pdg, mass_text = fields
    try:
        mass = float(mass_text)
        validate_mass_overrides({pdg: mass}, "Mass override")
    except (RuntimeError, ValueError) as error:
        raise argparse.ArgumentTypeError(str(error)) from error
    return pdg, mass


# Parse the regeneration command line
def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, default=DEFAULT_MANIFEST)
    parser.add_argument("--mg5", help="Path to bin/mg5_aMC")
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--add", metavar="NAME")
    parser.add_argument("--process", dest="mg5_process")
    parser.add_argument("--add-family", metavar="NAME")
    parser.add_argument("--remove", metavar="NAME")
    parser.add_argument("--remove-family", metavar="NAME")
    parser.add_argument("--family-process", action="append", default=[])
    parser.add_argument("--definition", action="append", default=[])
    parser.add_argument("--projection", choices=sorted(process_registry.PROJECTIONS))
    parser.add_argument("--type", dest="process_type", choices=("tree", "loop"), default="tree")
    parser.add_argument("--channel")
    parser.add_argument("--model", default="sm")
    parser.add_argument("--model-import")
    parser.add_argument("--complex-mass-scheme", action="store_true")
    parser.add_argument("--alpha-charge")
    parser.add_argument("--alpha-charge-square")
    parser.add_argument("--inverse-alpha-zero", type=float, default=137.03599908)
    parser.add_argument(
        "--mass-override",
        action="append",
        type=parse_mass_override,
        default=[],
        metavar="PDG=MASS",
    )
    return parser.parse_args()


# Execute full regeneration or validate the installed generated outputs
def main() -> None:
    args = parse_args()
    args.manifest = args.manifest.resolve()
    original_manifest = load_manifest(args.manifest)
    manifest = add_registry_entry(args, original_manifest)
    registry_change = manifest != original_manifest
    if args.validate_only:
        if registry_change:
            raise RuntimeError("Registry changes require full regeneration")
        validate_installed(manifest)
        print(
            f"Validated {len(manifest['processes'])} standalone processes and "
            f"{len(manifest.get('families', []))} subprocess families"
        )
        return
    transaction = InstallTransaction()
    try:
        mg5 = resolve_mg5(args.mg5)
        work_dir = Path(tempfile.mkdtemp(prefix="MG2GRA_", dir=ROOT / "tmp"))
        commands, process_exports, support_exports, family_exports = build_mg5_commands(
            manifest, work_dir
        )
        command_file = work_dir / "generate.mg5"
        command_file.write_text(commands)
        configuration = work_dir / "mg5_configuration"
        ensure_dir(configuration, exist_ok=False)
        (configuration / "mg5_configuration.txt").write_text(
            "nb_core = 4\n"
            "ninja = None\n"
            "collier = None\n"
            "output_dependencies = internal\n"
            "automatic_html_opening = False\n"
            "notification_center = False\n"
        )
        mg5_env = os.environ.copy()
        mg5_env["PYTHONHASHSEED"] = "0"
        mg5_env["MADGRAPH_BASE"] = str(configuration)
        subprocess.run([str(mg5), str(command_file)], cwd=work_dir, check=True, env=mg5_env)

        validate_support_owners(support_exports)
        slha_export = validate_slha_exports(support_exports, family_exports)
        retire_legacy_cpp_layout(transaction)
        retire_removed_outputs(original_manifest, manifest, transaction)
        retire_orphaned_outputs(manifest, transaction)
        retire_orphaned_model_support(support_exports, transaction)
        sources_by_model = install_standalone_processes(manifest, process_exports, transaction)
        install_loop_processes(manifest, mg5, process_exports, work_dir, transaction)
        install_durham_process_data(manifest, mg5, process_exports, work_dir, transaction)
        install_photon_process_data(manifest, mg5, process_exports, work_dir, transaction)
        model_particles = {name: load_model_particles(mg5.parent.parent, model)
                           for name, model in manifest["models"].items()}
        for model_name, export in support_exports.items():
            model = manifest["models"][model_name]
            install_model_support(
                model_name, model, export, sources_by_model[model_name], model_particles[model_name], transaction
            )
        install_slha_reader(slha_export, transaction)
        family_by_name = {family["name"]: family for family in manifest.get("families", [])}
        for family_name, export in family_exports.items():
            family = family_by_name[family_name]
            model = manifest["models"][family["model"]]
            color_structure = mg5_family.generate_family_color_structure(
                mg5.parent.parent,
                work_dir,
                model["import"],
                family.get("definitions", []),
                family["processes"],
                model.get("complex_mass_scheme", False),
            )
            install_family(family, model, export, color_structure, model_particles[family["model"]], transaction)
        install_registries(manifest, transaction)
        install_text(args.manifest, dump_manifest(manifest), transaction)
        validate_installed(manifest)
        transaction.commit()
        print(
            f"Regenerated {len(manifest['processes'])} standalone processes and "
            f"{len(family_exports)} subprocess families under {work_dir}"
        )
    except Exception:
        transaction.rollback()
        raise


if __name__ == "__main__":
    main()
