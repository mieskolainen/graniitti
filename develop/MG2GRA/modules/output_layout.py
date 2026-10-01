#!/usr/bin/env python3
#
# Define the generated MG5 C++ output layout
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re

RUNTIME = "Runtime"
DURHAM = "Durham"
PHOTON = "Photon"
PARTON = "Parton"

INCLUDE_ROOT = "Graniitti/Amplitude/MG5"
CPP_INCLUDE_ROOT = f"include/{INCLUDE_ROOT}"
CPP_SOURCE_ROOT = "src/Amplitude/MG5"
CARDS_ROOT = "MG5cards"

DURHAM_PROJECTION = "durham"
PHOTON_PROJECTION = "photon"
PARTON_PROJECTION = "parton"
PROJECTIONS = frozenset({DURHAM_PROJECTION, PHOTON_PROJECTION, PARTON_PROJECTION})

PARAMETER_CARD = "param_card.dat"
COLOR_DATA = "color.json"
CHANNEL_DATA = "channels.json"
PROCESS_REGISTRY = "process_registry.json"
CONVERTER_STATE = "converter.json"

_PROJECTION_DIRECTORIES = {
    DURHAM_PROJECTION: DURHAM,
    PHOTON_PROJECTION: PHOTON,
    PARTON_PROJECTION: PARTON,
}

_CPP_ROOTS = frozenset({RUNTIME, DURHAM, PHOTON, PARTON})


# Validate one relative generated C++ directory
def validate_cpp_directory(directory: str) -> str:
    if not isinstance(directory, str) or not directory:
        raise RuntimeError(f"Invalid MG5 C++ directory {directory!r}")
    parts = directory.split("/")
    if parts[0] not in _CPP_ROOTS or any(
        not part or part in {".", ".."} or validate_owner(part) != part
        for part in parts
    ):
        raise RuntimeError(f"Invalid MG5 C++ directory {directory!r}")
    return directory


# Validate one generated C++ filename
def validate_cpp_filename(filename: str, suffix: str) -> str:
    if (
        not isinstance(filename, str)
        or re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*\." + re.escape(suffix), filename) is None
    ):
        raise RuntimeError(f"Invalid MG5 C++ filename {filename!r}")
    return filename


# Compute one public generated MG5 include path
def include_path(directory: str, filename: str) -> str:
    return (
        f"{INCLUDE_ROOT}/{validate_cpp_directory(directory)}/"
        f"{validate_cpp_filename(filename, 'h')}"
    )


# Compute one repository relative generated MG5 header path
def header_path(directory: str, filename: str) -> str:
    return (
        f"{CPP_INCLUDE_ROOT}/{validate_cpp_directory(directory)}/"
        f"{validate_cpp_filename(filename, 'h')}"
    )


# Compute one repository relative generated MG5 source path
def cpp_source_path(directory: str, filename: str) -> str:
    return (
        f"{CPP_SOURCE_ROOT}/{validate_cpp_directory(directory)}/"
        f"{validate_cpp_filename(filename, 'cc')}"
    )


# Compute the generated C++ directory selected by one physical projection
def projection_cpp_directory(projection: str) -> str:
    return projection_directory(projection)


# Compute one generated family C++ directory
def family_cpp_directory(projection: str, family: str) -> str:
    if projection not in {PHOTON_PROJECTION, PARTON_PROJECTION}:
        raise RuntimeError(f"Projection {projection!r} has no generated family sources")
    return f"{projection_cpp_directory(projection)}/{validate_owner(family)}"


# Compute one common model support C++ directory
def model_cpp_directory(model: str) -> str:
    return f"{RUNTIME}/Models/{validate_owner(model)}"


# Compute the physical card directory selected by one projection
def projection_directory(projection: str) -> str:
    try:
        return _PROJECTION_DIRECTORIES[projection]
    except (KeyError, TypeError) as error:
        raise RuntimeError(f"Invalid MG5 projection {projection!r}") from error


# Validate one manifest process or family identifier used as a directory
def validate_owner(owner: str) -> str:
    if not isinstance(owner, str) or re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_]*", owner) is None:
        raise RuntimeError(f"Invalid MG5 card owner {owner!r}")
    return owner


# Compute one process or family card directory relative to MG5cards
def cards_directory(projection: str, owner: str) -> str:
    return f"{projection_directory(projection)}/{validate_owner(owner)}"


# Compute one process or family parameter card relative to MG5cards
def parameter_card_path(projection: str, owner: str) -> str:
    return f"{cards_directory(projection, owner)}/{PARAMETER_CARD}"


# Compute one standalone color sidecar relative to MG5cards
def color_path(projection: str, owner: str) -> str:
    if not isinstance(projection, str) or projection not in {
        DURHAM_PROJECTION,
        PHOTON_PROJECTION,
    }:
        raise RuntimeError(f"Projection {projection!r} has no standalone color data")
    return f"{cards_directory(projection, owner)}/{COLOR_DATA}"


# Compute one family channel sidecar relative to MG5cards
def channels_path(projection: str, owner: str) -> str:
    if not isinstance(projection, str) or projection not in {
        PHOTON_PROJECTION,
        PARTON_PROJECTION,
    }:
        raise RuntimeError(f"Projection {projection!r} has no family channel data")
    return f"{cards_directory(projection, owner)}/{CHANNEL_DATA}"


# Compute one generated runtime data path relative to MG5cards
def runtime_path(filename: str) -> str:
    if not isinstance(filename, str) or filename not in {PROCESS_REGISTRY, CONVERTER_STATE}:
        raise RuntimeError(f"Invalid MG5 runtime data file {filename!r}")
    return f"{RUNTIME}/{filename}"


# Compute the generated common process registry relative to MG5cards
def registry_path() -> str:
    return runtime_path(PROCESS_REGISTRY)


# Compute the generated converter provenance relative to MG5cards
def converter_path() -> str:
    return runtime_path(CONVERTER_STATE)


# Prefix one path relative to MG5cards for generated runtime source
def source_path(relative_path: str) -> str:
    if not isinstance(relative_path, str) or not relative_path:
        raise RuntimeError(f"Invalid MG5 card path {relative_path!r}")
    parts = relative_path.split("/")
    if parts[0] not in {RUNTIME, DURHAM, PHOTON, PARTON} or any(
        not part or part in {".", ".."} for part in parts
    ):
        raise RuntimeError(f"Invalid MG5 card path {relative_path!r}")
    return f"{CARDS_ROOT}/{relative_path}"
