#!/usr/bin/env python3
#
# Generate scale-specific couplings from converted MG5 model support
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import re

from .cpp_support import compact_cpp_whitespace, cpp_function_statements, replace_cpp_identifiers


# Extract generated scalar assignments in their dependency order
def assignments(source: str, parameter_class: str, method: str) -> list[tuple[str, str]]:
    signature = f"void {parameter_class}::{method}()"
    if method == "setIndependentParameters":
        signature = f"void {parameter_class}::{method}(SLHAReader"
    if signature not in source:
        return []
    rows = []
    for statement in cpp_function_statements(source, signature):
        match = re.fullmatch(r"\s*([A-Za-z_]\w*)\s*=(?!=)\s*(.*?)\s*;\s*", statement, flags=re.S)
        if match is not None:
            rows.append(match.groups())
    return rows


# Generate the alpha(0) assignments required by the selected model support
def alpha_zero_definition(
    raw_source: str,
    parameter_class: str,
    charge: str | None,
    charge_square: str | None,
    inverse_alpha: float,
) -> str | None:
    if charge is None:
        return None
    replacements = {charge: f"{charge}_NEW"}
    if charge_square is not None:
        replacements[charge_square] = f"{charge_square}_NEW"
    parameters = []
    for method in ("setIndependentParameters", "setDependentParameters"):
        for target, expression in assignments(raw_source, parameter_class, method):
            if target in replacements:
                continue
            replacement, used = replace_cpp_identifiers(expression, replacements)
            if used:
                replacements[target] = f"{target}_NEW"
                parameters.append((target, compact_cpp_whitespace(replacement)))
    selected: list[tuple[str, str]] = []
    for target, expression in (
        assignments(raw_source, parameter_class, "setIndependentCouplings")
        + assignments(raw_source, parameter_class, "setDependentCouplings")
    ):
        replacement, used = replace_cpp_identifiers(expression, replacements)
        if not used:
            continue
        selected.append((target, compact_cpp_whitespace(replacement)))

    if not selected:
        raise RuntimeError(
            f"No generated coupling depends on configured electromagnetic charge {charge}"
        )

    needed = set(re.findall(r"\b\w+_NEW\b", " ".join(expression for _, expression in selected)))
    retained = []
    for target, expression in reversed(parameters):
        if f"{target}_NEW" in needed:
            retained.append((target, expression))
            needed.update(re.findall(r"\b\w+_NEW\b", expression))
    parameters = list(reversed(retained))

    lines = [
        f"void {parameter_class}::setAlphaQEDZero() {{",
        f"  const double {charge}_NEW = 2.0 * std::sqrt(1.0 / {inverse_alpha:.11f}) * std::sqrt(M_PI);",
    ]
    if charge_square is not None and f"{charge_square}_NEW" in needed:
        lines.append(f"  const double {charge_square}_NEW = {charge}_NEW * {charge}_NEW;")
    lines.append("")
    lines.extend(f"  const auto {target}_NEW = {expression};" for target, expression in parameters)
    lines.extend(f"  {target} = {expression};" for target, expression in selected)
    lines.append("}")
    return "\n".join(lines)
