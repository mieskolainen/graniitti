# Shared phase and output helpers for repository tools
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from __future__ import annotations

import json
import math
from collections.abc import Iterable, Sequence
from typing import Any


# Compute the equivalent coupling phase in the C++ card interval [-pi, pi)
def phase(phase: float) -> float:
    if not math.isfinite(phase):
        raise ValueError("coupling phase must be finite")
    value = math.remainder(phase, 2.0 * math.pi)
    return value - 2.0 * math.pi if value >= math.pi else value


# Compute an aligned plain-text table with optional blank lines between groups
def table(
    headers: Sequence[object],
    rows: Iterable[Sequence[object]],
    *,
    right_align: set[int] | None = None,
    group_by: int | None = None,
) -> str:
    body = [[str(value) for value in row] for row in rows]
    if not body:
        return "(no rows)"
    head = [str(value) for value in headers]
    if any(len(row) != len(head) for row in body):
        raise ValueError("table rows and headers must have equal column counts")
    widths = [max(len(head[i]), *(len(row[i]) for row in body)) for i in range(len(head))]
    numeric = right_align or set()

    # Align one row according to its column roles
    def aligned(row: Sequence[str]) -> str:
        cells = []
        for i, value in enumerate(row):
            cells.append(value.rjust(widths[i]) if i in numeric else value.ljust(widths[i]))
        return "  ".join(cells).rstrip()

    separator = "  ".join("-" * width for width in widths)
    lines = [aligned(head), separator]
    previous = None
    for row in body:
        current = row[group_by] if group_by is not None else None
        if previous is not None and current != previous:
            lines.append("")
        lines.append(aligned(row))
        previous = current
    return "\n".join(lines)


# Compute true for values that fit on one compact JSON line
def _scalar(value: Any) -> bool:
    return value is None or isinstance(value, (bool, int, float, str))


# Serialize JSON while keeping scalar vectors and numeric matrix rows compact
def dumps(value: Any, indent: int = 0) -> str:
    pad = " " * indent
    if _scalar(value):
        return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"))
    if isinstance(value, list):
        if not value:
            return "[]"
        if all(_scalar(item) for item in value):
            return "[" + ",".join(dumps(item) for item in value) + "]"
        if all(isinstance(item, list) and all(_scalar(x) for x in item) for item in value):
            rows = ["[" + ",".join(dumps(x) for x in item) + "]" for item in value]
            row_pad = " " * (indent + 1)
            return "[" + (",\n" + row_pad).join(rows) + "]"
        rows = [dumps(item, indent + 2) for item in value]
        return "[\n" + ",\n".join(" " * (indent + 2) + row for row in rows) + "\n" + pad + "]"
    if isinstance(value, dict):
        rows = []
        for key, item in value.items():
            label = json.dumps(str(key), ensure_ascii=False)
            rows.append(" " * (indent + 2) + label + ": " + dumps(item, indent + 2))
        return "{\n" + ",\n".join(rows) + "\n" + pad + "}"
    raise TypeError(f"cannot serialize value of type {type(value).__name__}")


# Print a titled table
def print_table(title, headers, rows, **options):
    print(title)
    print(table(headers, rows, **options))
