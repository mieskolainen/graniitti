#!/usr/bin/env python3
#
# Shared sidecar validation and C++ formatting for MG5 registries
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import math
import re
from typing import Any


# Compute one exact integer field and reject bool or floating substitutions
def exact_int(data: dict[str, Any], field: str, label: str, positive: bool = False) -> int:
    value = data[field]
    if type(value) is not int or (positive and value <= 0):
        raise RuntimeError(f"Invalid {label} {field}")
    return value


# Compute one finite numeric value and reject bool or encoded strings
def finite_float(
    value: Any, label: str, nonnegative: bool = False
) -> float:
    if type(value) not in (int, float):
        raise RuntimeError(f"Invalid {label}")
    numeric = float(value)
    if not math.isfinite(numeric) or (nonnegative and numeric < 0.0):
        raise RuntimeError(f"Invalid {label}")
    return numeric


# Compute one exact integer vector
def int_vector(value: Any, label: str, size: int | None = None) -> list[int]:
    if type(value) is not list or (size is not None and len(value) != size):
        raise RuntimeError(f"Invalid {label}")
    if any(type(item) is not int for item in value):
        raise RuntimeError(f"Invalid {label}")
    return value


# Validate the shared finite-Nc color sidecar schema
def validate_color_data(
    process_data: dict[str, Any],
    entry: dict[str, Any],
    incoming_pdgs: list[int],
    version: int,
    label: str,
) -> None:
    if type(process_data) is not dict:
        raise RuntimeError(f"Invalid {label} process_data for {entry['name']}")
    required = {
        "version",
        "process",
        "incoming_pdgs",
        "final_pdgs",
        "final_color_representations",
        "ncolor",
        "rank",
        "basis_sha256",
        "factorization_residual",
        "projectors",
        "flow_candidates",
        "flow_weights",
    }
    missing = required - set(process_data)
    if missing:
        raise RuntimeError(
            f"{label} process_data for {entry['name']} lacks {sorted(missing)}"
        )
    if type(process_data["version"]) is not int or process_data["version"] != version:
        raise RuntimeError(f"Stale {label} process_data for {entry['name']}")
    if (
        type(process_data["process"]) is not str
        or process_data["process"] != entry["process"]
    ):
        raise RuntimeError(f"Stale {label} process_data for {entry['name']}")

    incoming = int_vector(
        process_data["incoming_pdgs"],
        f"{label} incoming state for {entry['name']}",
        len(incoming_pdgs),
    )
    if incoming != incoming_pdgs:
        raise RuntimeError(f"Invalid {label} incoming state for {entry['name']}")
    final_pdgs = int_vector(
        process_data["final_pdgs"], f"{label} final state for {entry['name']}"
    )
    if not final_pdgs:
        raise RuntimeError(f"Invalid {label} final state for {entry['name']}")
    representations = int_vector(
        process_data["final_color_representations"],
        f"{label} final-state colors for {entry['name']}",
        len(final_pdgs),
    )
    supported_representations = {1, 3, -3, 6, -6, 8}
    if any(value not in supported_representations for value in representations):
        raise RuntimeError(
            f"Unsupported {label} final-state color representation for {entry['name']}"
        )

    ncolor = exact_int(process_data, "ncolor", f"{label} process_data", positive=True)
    rank = exact_int(process_data, "rank", f"{label} process_data", positive=True)
    if rank > ncolor:
        raise RuntimeError(f"Invalid {label} projector rank for {entry['name']}")
    digest = process_data["basis_sha256"]
    if type(digest) is not str or re.fullmatch(r"[0-9a-f]{64}", digest) is None:
        raise RuntimeError(f"Invalid {label} color-basis digest for {entry['name']}")
    finite_float(
        process_data["factorization_residual"],
        f"{label} factorization residual for {entry['name']}",
        nonnegative=True,
    )

    projectors = process_data["projectors"]
    if type(projectors) is not list or len(projectors) != rank:
        raise RuntimeError(f"Invalid {label} projector rank for {entry['name']}")
    for row in projectors:
        if type(row) is not list or len(row) != ncolor:
            raise RuntimeError(f"Invalid {label} projector width for {entry['name']}")
        for value in row:
            if type(value) is not list or len(value) != 2:
                raise RuntimeError(f"Invalid {label} projector value for {entry['name']}")
            finite_float(value[0], f"{label} projector value for {entry['name']}")
            finite_float(value[1], f"{label} projector value for {entry['name']}")

    candidates = process_data["flow_candidates"]
    weights = process_data["flow_weights"]
    if type(candidates) is not list or type(weights) is not list:
        raise RuntimeError(f"Invalid {label} shower-flow partition for {entry['name']}")
    if len(candidates) != len(weights):
        raise RuntimeError(f"Invalid {label} shower-flow partition for {entry['name']}")
    has_final_color = any(value != 1 for value in representations)
    if has_final_color == (not candidates):
        raise RuntimeError(
            f"Invalid {label} shower-flow availability for {entry['name']}"
        )
    for candidate in candidates:
        if type(candidate) is not list or len(candidate) != len(final_pdgs):
            raise RuntimeError(f"Invalid {label} shower-flow size for {entry['name']}")
        tag_counts: dict[int, list[int]] = {}
        for representation, flow in zip(representations, candidate, strict=True):
            pair = int_vector(flow, f"{label} shower-flow value for {entry['name']}", 2)
            color, anticolor = pair
            valid_slots = (
                (representation == 1 and color == 0 and anticolor == 0)
                or (representation == 3 and color > 0 and anticolor == 0)
                or (representation == -3 and color == 0 and anticolor > 0)
                or (representation == 6 and color > 0 and anticolor < 0)
                or (representation == -6 and color < 0 and anticolor > 0)
                or (
                    representation == 8
                    and color > 0
                    and anticolor > 0
                    and color != anticolor
                )
            )
            if not valid_slots:
                raise RuntimeError(
                    f"Invalid {label} shower-flow representation for {entry['name']}"
                )
            if color:
                tag_counts.setdefault(abs(color), [0, 0])[int(color < 0)] += 1
            if anticolor:
                tag_counts.setdefault(abs(anticolor), [0, 0])[int(anticolor > 0)] += 1
        if any(counts != [1, 1] for counts in tag_counts.values()):
            raise RuntimeError(
                f"Unbalanced {label} shower-flow tags for {entry['name']}"
            )
    for row in weights:
        if type(row) is not list or len(row) != ncolor:
            raise RuntimeError(f"Invalid {label} shower-flow weights for {entry['name']}")
        for weight in row:
            finite_float(
                weight,
                f"{label} shower-flow weight for {entry['name']}",
                nonnegative=True,
            )
    if weights:
        for column in range(ncolor):
            total = math.fsum(float(row[column]) for row in weights)
            if not math.isclose(total, 1.0, rel_tol=0.0, abs_tol=1.0e-12):
                raise RuntimeError(
                    f"Non-unit {label} shower-flow partition for {entry['name']}"
                )


# Format one floating-point literal reproducibly for generated C++
def cpp_float(value: float) -> str:
    numeric = finite_float(value, "registry floating-point literal")
    if not numeric:
        return "0.0"
    return format(numeric, ".17g")


# Format one serialized complex value for generated C++
def cpp_complex(value: list[float]) -> str:
    if type(value) is not list or len(value) != 2:
        raise RuntimeError("Invalid registry complex value")
    return f"std::complex<double>({cpp_float(value[0])}, {cpp_float(value[1])})"


# Format one flat complex matrix for generated C++
def cpp_complex_matrix(matrix: list[list[list[float]]]) -> str:
    values = [cpp_complex(value) for row in matrix for value in row]
    return "{\n      " + ",\n      ".join(values) + "\n  }"


# Format one integer vector for generated C++
def cpp_int_vector(values: list[int]) -> str:
    int_vector(values, "registry integer vector")
    return "{" + ", ".join(str(value) for value in values) + "}"


# Format all generated shower color-flow candidates
def cpp_flow_candidates(candidates: list[list[list[int]]]) -> str:
    if not candidates:
        return "{}"
    formatted = []
    for candidate in candidates:
        entries = ", ".join(
            f"MColorFlow{{{flow[0]}, {flow[1]}}}" for flow in candidate
        )
        formatted.append("{" + entries + "}")
    return "{\n      " + ",\n      ".join(formatted) + "\n  }"


# Expand exact projector rows into exact and leading-color candidate rows
def combined_projectors(process_data: dict[str, Any]) -> list[list[list[float]]]:
    exact = process_data["projectors"]
    weights = process_data["flow_weights"]
    combined = list(exact)
    for flow in weights:
        for projector in exact:
            combined.append(
                [
                    [
                        float(value[0]) * float(flow[column]),
                        float(value[1]) * float(flow[column]),
                    ]
                    for column, value in enumerate(projector)
                ]
            )
    return combined
