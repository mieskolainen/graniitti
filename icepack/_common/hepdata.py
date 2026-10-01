# Shared HEPData table decoding and measured bin transformations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import json
import re

import numpy as np
from core.io.cache import cache
from core.stats.hist import (
    bins2binwidth,
    center2edgebins,
    doublet2linear,
    edge2centerbins,
    rebin_histogram_groups,
    rebin_partition,
)
from core.stats.uncertainty import (
    covariance_square_root,
    finalize_uncertainties,
    linear_source,
    split_uncertainty,
    uncertainty_source,
)
from termcolor import cprint


# Read original HEPData JSON without interpreting steering references
def read_table(filename, dimensions=None):
    with open(filename, encoding="utf-8") as stream:
        table = json.load(stream)
    if not isinstance(table, dict) or not {"headers", "qualifiers", "values", "x_count"}.issubset(table):
        raise ValueError(f"Invalid HEPData table JSON: {filename}")
    if dimensions is not None and (int(table["x_count"]) != dimensions or not table["values"]):
        raise ValueError(f"HEPData requires {dimensions} nonempty independent variables: {filename}")
    return table


# Convert one HEPData numeric field into a finite float
def number(value):
    if isinstance(value, str) and value.strip().endswith("%"):
        raise ValueError("HEPData: relative percent errors are not supported")
    output = float(value)
    if not np.isfinite(output):
        raise ValueError("HEPData: nonfinite numeric field")
    return output


# Normalize one HEPData label for exact source and group matching
def normalized(value):
    return re.sub(r"[^a-z0-9]+", "", str(value).lower())


# Resolve an optional HEPData dependent variable qualifier
def resolve_group(table, file_filter, *, headers=False):
    groups = sorted({int(item["group"]) for row in table["values"] for item in row["y"]})
    if file_filter is None:
        if len(groups) != 1:
            raise ValueError("HEPData: file_filter is required for a multigroup table")
        return groups[0]

    selected = normalized(file_filter)
    if not selected:
        raise ValueError("HEPData: empty dependent variable filter")
    matches = set()
    for entries in table["qualifiers"].values():
        for entry in entries:
            group = int(entry["group"])
            if group in groups and int(entry.get("colspan", 1)) == 1 and normalized(entry["value"]) == selected:
                matches.add(group)
    if not matches and headers:
        for group, header in zip(groups, table["headers"][1:], strict=True):
            if selected in normalized(header["name"]):
                matches.add(group)
    if len(matches) != 1:
        raise ValueError(f"HEPData: filter {file_filter!r} matched groups {sorted(matches)}")
    return matches.pop()


# Decode signed absolute or percentage errors without discarding asymmetric variations
def error_pairs(values, label):
    pairs = []
    for value in values:
        if isinstance(label, str):
            matches = [error for error in value.get("errors", []) if error.get("label", "error") == label]
            if len(matches) != 1:
                raise ValueError(f"HEPData requires one error labelled {label!r}")
            error = matches[0]
        else:
            error = value["errors"][label]
        if "symerror" in error:
            raw = [error["symerror"]] * 2
        elif "asymerror" in error and {"plus", "minus"}.issubset(error["asymerror"]):
            raw = [error["asymerror"][side] for side in ("plus", "minus")]
        else:
            raise ValueError("HEPData: uncertainty requires symerror or complete asymerror")
        pair = []
        for item in raw:
            relative = isinstance(item, str) and item.strip().endswith("%")
            if relative and "value" not in value:
                raise ValueError("HEPData: relative errors require a central value")
            pair.append(number(item.strip().rstrip("%")) * number(value["value"]) / 100.0
                        if relative else number(item))
        if "symerror" in error:
            pair[1] = -pair[1]
        pairs.append(pair)
    result = np.asarray(pairs, dtype=float)
    if not np.isfinite(result).all():
        raise ValueError(f"HEPData has nonfinite errors for {label!r}")
    return result


# Compute absolute magnitudes with the common HEPData error decoder
def error_pair(error):
    return tuple(np.abs(error_pairs([{"errors": [error]}], 0)[0]))


# Select one dependent variable by its group or an exact published qualifier
def table_values(table, group=None):
    groups = {int(value["group"]) for row in table["values"] for value in row["y"]}
    if isinstance(group, str):
        matches = {int(entry["group"]) for entries in table["qualifiers"].values() for entry in entries
                   if entry["value"] == group and int(entry.get("colspan", 1)) == 1}
        if len(matches) != 1:
            raise ValueError(f"HEPData qualifier {group!r} does not select one dependent variable")
        group = matches.pop()
    if group is None:
        if len(groups) != 1:
            raise ValueError("HEPData table requires a dependent variable selection")
        group = groups.pop()
    values = []
    for row in table["values"]:
        selected = [value for value in row["y"] if int(value["group"]) == group]
        if len(selected) != 1:
            raise ValueError(f"HEPData row requires one dependent value for group {group}")
        values.append(selected[0])
    return values


# Select explicitly published rows without changing the source table
def select_rows(table, rows):
    if rows is None:
        return table
    if any(isinstance(index, (bool, np.bool_)) or not isinstance(index, (int, np.integer)) for index in rows):
        raise ValueError("HEPData: rows must contain integer table row indices")
    selected = [int(index) for index in rows]
    if not selected or len(selected) != len(set(selected)):
        raise ValueError("HEPData: rows must contain unique table row indices")
    if min(selected) < 0 or max(selected) >= len(table["values"]):
        raise ValueError("HEPData: table row index is out of range")
    output = copy.deepcopy(table)
    output["values"] = [output["values"][index] for index in selected]
    return output


# Decode published intervals and infer point intervals before any row selection
def table_bins(table, axis=0, *, reference=None, points=False):
    coordinates = [row["x"][axis] for row in table["values"]]
    if not coordinates or any(("low" in point) != ("high" in point) for point in coordinates):
        raise ValueError("HEPData axis is empty or has incomplete bin bounds")
    bounded = ["low" in point and "high" in point for point in coordinates]
    if all(bounded) and not points:
        edges = np.asarray([[point["low"], point["high"]] for point in coordinates], dtype=float)
        centers = np.asarray([float(point["value"]) if "value" in point else np.mean(pair)
                              for point, pair in zip(coordinates, edges, strict=True)])
    elif points or not any(bounded):
        centers = np.asarray([point["value"] for point in coordinates], dtype=float)
        full = centers if reference is None else np.sort([number(row["x"][axis]["value"]) for row in reference["values"]])
        if len(full) == 1 and not points:
            raise ValueError("A single HEPData point requires explicit bin bounds")
        bins = center2edgebins(full)
        indices = np.searchsorted(full, centers)
        if np.any(indices >= len(full)) or not np.array_equal(full[indices], centers):
            raise ValueError("HEPData selected points are missing from the original axis")
        edges = np.column_stack((bins[indices], bins[indices + 1]))
    else:
        raise ValueError("HEPData axis mixes points and bounded bins")
    if not np.isfinite(centers).all() or not np.isfinite(edges).all() or np.any(edges[:, 1] <= edges[:, 0]):
        raise ValueError("HEPData has invalid coordinates or bin widths")
    if np.any(centers < edges[:, 0]) or np.any(centers > edges[:, 1]):
        raise ValueError("HEPData coordinates lie outside their bin bounds")
    return centers, edges


# Embed measured intervals and all their uncertainties into a histogram with masked gaps
def fill_bins(ds):
    bins = doublet2linear(ds["binedges"])
    indices = np.searchsorted(bins, np.asarray(ds["binedges"])[:, 1]) - 1
    valid = np.zeros(len(bins) - 1, dtype=bool)
    valid[indices] = ds.get("valid", True)
    if len(valid) != len(ds["y"]):
        transform = np.zeros((len(valid), len(ds["y"])))
        transform[indices, np.arange(len(ds["y"]))] = 1.0
        marginal = {key: transform @ value for key, value in ds.items() if key.startswith("y_err_")}
        centers = edge2centerbins(bins)
        centers[indices] = ds["x"]
        ds["x"], ds["y"] = centers, transform @ ds["y"]
        finalize_uncertainties(ds, [linear_source(source, transform) for source in ds["uncertainties"]])
        ds.update(marginal)
    ds.update(bins=bins, binedges=np.column_stack((bins[:-1], bins[1:])),
              binwidth=bins2binwidth(bins), valid=valid, xlim=bins[[0, -1]])
    return ds


# Separate a common luminosity scale using its original cross section and quoted uncertainty
def split_luminosity(ds, table, label, source):
    value = table_values(table, group=0)[0]
    fraction = np.abs(error_pairs([value], label)[0]) / number(value["value"])
    common = [uncertainty_source(
        "luminosity", ds["y"] * fraction[0], down=ds["y"] * fraction[1], category="systematic",
        correlation="collective", scope=f"{table['doi'].rsplit('/', 1)[0]}:luminosity", effect="multiplicative",
        provenance=f"Original HEPData {table['doi']}, {label} divided by its cross section",
    )]
    return finalize_uncertainties(ds, split_uncertainty(ds["uncertainties"], common, source=source))


# Rebin measured intervals and their uncertainty sources without crossing gaps
def rebin_common(ds: dict, rebin_factor: int = None) -> dict:
    old_edges = np.asarray(ds["bins"], dtype=float)
    valid = np.asarray(ds.get("valid", np.ones(len(ds["y"]), dtype=bool)))
    groups = rebin_partition(valid, rebin_factor)
    new_edges, new_content, new_error = rebin_histogram_groups(
        old_edges, ds["y"], ds["y_err"], groups, differential=True
    )
    centers = edge2centerbins(new_edges)
    single = np.diff(groups) == 1
    centers[single] = np.asarray(ds["x"])[groups[:-1][single]]
    ds["x"] = centers
    ds["bins"] = new_edges
    ds["binwidth"] = bins2binwidth(new_edges)
    ds["y"] = new_content
    ds["y_err"] = new_error
    if "binedges" in ds:
        ds["binedges"] = np.column_stack((new_edges[:-1], new_edges[1:]))
    if "valid" in ds:
        ds["valid"] = valid[groups[:-1]]
    if "uncertainties" in ds:
        transform = np.zeros((len(new_content), len(old_edges) - 1), dtype=float)
        old_width = bins2binwidth(old_edges)
        for row, (start, stop) in enumerate(zip(groups[:-1], groups[1:], strict=True)):
            transform[row, start:stop] = old_width[start:stop] / ds["binwidth"][row]
        finalize_uncertainties(ds, [linear_source(source, transform, full=True) for source in ds["uncertainties"]])
    else:
        ds["y_err_stat"] = np.zeros_like(new_content)
        ds["y_err_syst"] = np.zeros_like(new_content)
    return ds


# Print one common reader rebin message
def print_rebin(reader: str, rebin_factor: int) -> None:
    cprint(f"{reader}: Rebin histogram with factor // {rebin_factor}", "red")


# Read a complete HEPData covariance table in its published bin order
@cache
def read_covariance(filename, indices=None):
    rows = read_table(filename)["values"]
    bins = np.asarray([[int(item["value"]) - 1 for item in row["x"]] for row in rows])
    if bins.ndim != 2 or bins.shape[1] != 2 or np.any(bins < 0):
        raise ValueError(f"Invalid covariance bin indices: {filename}")
    size = int(np.max(bins)) + 1
    if len(rows) != size * size or len(np.unique(bins, axis=0)) != len(rows):
        raise ValueError(f"Incomplete or duplicate covariance entries: {filename}")
    matrix = np.full((size, size), np.nan)
    matrix[bins[:, 0], bins[:, 1]] = [float(row["y"][0]["value"]) for row in rows]
    if not np.isfinite(matrix).all() or not np.allclose(matrix, matrix.T):
        raise ValueError(f"Nonfinite or asymmetric covariance: {filename}")
    covariance_square_root(matrix)
    selected = np.arange(size) if indices is None else np.asarray(indices)
    return matrix[np.ix_(selected, selected)]


# Read signed nuisance responses shared between categories in each source bin
def signed_sources(table, values, label, scope):
    sources = []
    for index, (row, value) in enumerate(zip(table["values"], values, strict=True)):
        errors = [error for error in value.get("errors", []) if error["label"] == label]
        if not errors:
            continue
        if len(errors) != 1 or "asymerror" not in errors[0]:
            raise ValueError(f"HEPData: {label} requires one signed asymmetric error")
        plus, minus = error_pairs([value], label)[0]
        coordinate = row["x"][0]
        bounds = tuple(number(coordinate[key]) for key in ("low", "high"))
        up, down = np.zeros(len(values)), np.zeros(len(values))
        up[index], down[index] = abs(plus), abs(minus)
        source = uncertainty_source(
            label, up, down=down, category="systematic", correlation="collective",
            scope=f"{scope}:{bounds}", provenance=f"Signed HEPData {label} responses at the same source bin",
        )
        source["shift"][index] = 0.5 * (plus - minus)
        sources.append(source)
    return sources
