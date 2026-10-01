# Shared reader for official UPC HEPData JSON table responses
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import doublet2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, residual_uncertainty, uncertainty_source

from icepack._common.hepdata import (
    error_pairs,
    fill_bins,
    normalized,
    number,
    read_covariance,
    read_table,
    resolve_group,
    select_rows,
    table_bins,
    table_values,
)


# Convert one HEPData label into a concise uncertainty source name
def source_name(value):
    output = re.sub(r"[^a-z0-9]+", "_", str(value).lower()).strip("_")
    if not output:
        raise ValueError("ReadHEPData_UPC: empty uncertainty label")
    return output


# Read published bins, point coordinates or the standard fiducial cross section bin
def read_bins(table, point_axis=False, reference=None):
    if len(table["values"]) == 1:
        coordinates = table["values"][0].get("x", [])
        if (
            len(coordinates) == 1
            and "low" not in coordinates[0]
            and "high" not in coordinates[0]
            and normalized(coordinates[0].get("value", "")) == "crosssection"
        ):
            return (
                np.asarray([0.5], dtype=float),
                np.asarray([[0.0, 1.0]], dtype=float),
                np.asarray([0.0, 1.0], dtype=float),
            )

    unbinned = all("low" not in row["x"][0] and "high" not in row["x"][0] for row in table["values"])
    centers, edges = table_bins(table, reference=reference, points=point_axis or unbinned)
    return centers, edges, doublet2linear(edges)


# Classify one named HEPData uncertainty without inventing covariance
def source_type(label):
    key = normalized(label)
    if "stat" in key and "syst" not in key:
        return "statistical", "uncorrelated"
    if "uncorr" in key:
        return "combined" if key == "uncorr" else "systematic", "uncorrelated"
    if "corr" in key or "lumi" in key:
        return "systematic", "collective"
    return "systematic", "uncorrelated"


# Build source aware uncertainties from the selected dependent group
def read_uncertainties(values, scope):
    labels = [str(error["label"]) for error in values[0].get("errors", [])]
    if len(labels) != len(set(labels)):
        raise ValueError("ReadHEPData_UPC: duplicate uncertainty labels")
    for value in values[1:]:
        current = [str(error["label"]) for error in value.get("errors", [])]
        if current != labels:
            raise ValueError("ReadHEPData_UPC: uncertainty labels change between bins")

    normalized_labels = [normalized(label) for label in labels]
    total_index = normalized_labels.index("total") if "total" in normalized_labels else None
    statistical_indices = [index for index, label in enumerate(labels) if source_type(label)[0] == "statistical"]
    residual = None
    if total_index is not None and statistical_indices:
        total_pairs = np.abs(error_pairs(values, total_index))
        known_pairs = [np.abs(error_pairs(values, index)) for index in statistical_indices]
        residual = (
            residual_uncertainty(
                total_pairs[:, 0], [component[:, 0] for component in known_pairs]
            ),
            residual_uncertainty(
                total_pairs[:, 1], [component[:, 1] for component in known_pairs]
            ),
        )

    sources = []
    for index, label in enumerate(labels):
        if index == total_index and residual is not None:
            sources.append(
                uncertainty_source(
                    "systematic_from_total",
                    residual[0],
                    down=residual[1],
                    category="systematic",
                    correlation="uncorrelated",
                    scope=f"{scope}:systematic_from_total",
                    provenance="Quadrature residual of official HEPData total and statistical fields",
                )
            )
            continue
        pairs = np.abs(error_pairs(values, index))
        up, down = pairs.T
        category, correlation = source_type(label)
        sources.append(
            uncertainty_source(
                source_name(label),
                up,
                down=down,
                category=category,
                correlation=correlation,
                scope=f"{scope}:{source_name(label)}",
                provenance=(f"Official HEPData {label} field. " + (
                    "Same-sign common-shift approximation to the correlated marginal" if correlation == "collective"
                    else "Diagonal covariance where bin correlations are not supplied")),
            )
        )
    return sources


# Read one selected official UPC HEPData dependent variable
@cache
def _read(filename, file_filter, rows, point_axis, covariance, uncertainties):
    original = read_table(filename, dimensions=1)
    table = select_rows(original, rows)
    indices = np.arange(len(table["values"])) if rows is None else np.asarray(rows)
    if len(table["values"]) > 1:
        coordinates = [row["x"][0] for row in table["values"]]
        key = "value" if point_axis or all("low" not in item for item in coordinates) else "low"
        order = np.argsort([number(item[key]) for item in coordinates])
        indices = indices[order]
        table = {**table, "values": [table["values"][index] for index in order]}
    group = resolve_group(table, file_filter, headers=True)
    values = table_values(table, group)
    x, binedges, bins = read_bins(table, point_axis=point_axis, reference=original)
    y = np.asarray([number(value["value"]) for value in values], dtype=float)
    scope = f"{Path(filename).parent.name}:{Path(filename).stem}:group{group}"
    ds = {
        "x": x,
        "binedges": binedges,
        "bins": bins,
        "binwidth": doublet2binwidth(binedges),
        "xlim": np.asarray([bins[0], bins[-1]], dtype=float),
        "y": y,
        "group": group,
        "headers": table["headers"],
        "qualifiers": table["qualifiers"],
    }
    sources = read_uncertainties(values, scope) if uncertainties is None else uncertainties(table, values, scope)
    finalize_uncertainties(ds, sources)
    if covariance is not None:
        marginal = {key: value for key, value in ds.items() if key.startswith("y_err_") and key != "y_err_syst"}
        matrix = read_covariance(str(Path(filename).parent / covariance), tuple(indices))
        finalize_uncertainties(ds, [uncertainty_source(
            "published", matrix, category="combined", correlation="covariance",
            provenance=f"Original HEPData covariance table {covariance}",
        )])
        ds.update(marginal)
    return fill_bins(ds)


# Read one UPC table with an optional dependent variable filter
def read(filename, file_filter=None, rebin_factor=None, hist=None, uncertainties=None, **_kwargs):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_UPC: rebinning is not supported")
    rows = None if hist is None or "rows" not in hist else tuple(hist["rows"])
    point_axis = False if hist is None else bool(hist.get("point_axis", False))
    covariance = None if hist is None else hist.get("covariance")
    return _read(filename, file_filter, rows, point_axis, covariance, uncertainties)
