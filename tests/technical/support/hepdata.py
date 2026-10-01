# Read original HEPData fields for independent measurement comparisons in tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
from pathlib import Path

import numpy as np


# Extract quoted values and absolute down and up errors directly from JSON fields
def original(path, group=0, rows=None):
    table = json.loads(Path(path).read_bytes())
    selected = table["values"] if rows is None else [table["values"][index] for index in rows]
    points = [next(point for point in row["y"] if point["group"] == group) for row in selected]
    values = np.asarray([float(point["value"]) if point["value"] not in (None, "", "-") else np.nan
                         for point in points])
    errors = {}
    for index, point in enumerate(points):
        for error in point["errors"]:
            pair = [error["symerror"]] * 2 if "symerror" in error else [error["asymerror"][side] for side in ("minus", "plus")]
            absolute = [abs(float(str(value).rstrip("%"))) * (abs(values[index]) / 100 if str(value).endswith("%") else 1)
                        for value in pair]
            errors.setdefault(error.get("label", "error"), np.zeros((len(points), 2)))[index] = absolute
    return selected, values, errors


# Compare all numeric values and published intervals after an optional row selection
def assert_original(data, path, group=0, rows=None, density=False):
    selected, values, errors = original(path, group, rows)
    if density:
        widths = np.asarray([float(row["x"][0]["high"]) - float(row["x"][0]["low"]) for row in selected])
        values = values / widths
        errors = {name: pair / widths[:, None] for name, pair in errors.items()}
    if len(selected) > 1:
        order = np.argsort([float(row["x"][0].get("value", row["x"][0].get("low", 0))) for row in selected])
    else:
        order = np.arange(len(selected))
    valid = data.get("valid", np.ones(len(data["y"]), dtype=bool))
    np.testing.assert_allclose(data["y"][valid], values[order][np.isfinite(values[order])])
    if all("low" in row["x"][0] and "high" in row["x"][0] for row in selected):
        edges = np.asarray([[row["x"][0][side] for side in ("low", "high")] for row in selected], dtype=float)[order]
        np.testing.assert_allclose(np.column_stack((data["bins"][:-1], data["bins"][1:]))[valid], edges[np.isfinite(values[order])])
    if len(selected) > 1 and all("value" in row["x"][0] for row in selected):
        np.testing.assert_allclose(data["x"][valid], np.asarray([row["x"][0]["value"] for row in selected], dtype=float)[order])
    return values, errors
