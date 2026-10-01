# Integrated and differential theory reference reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.io.cache import cache
from core.io.serialize import load_json_file
from core.stats.hist import bins2binwidth, edge2centerbins
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source


# Convert a published theory prediction into a histogram with separately quoted uncertainties
def histogram(row, key):
    if row["type"] not in {"integrated", "differential"}:
        raise ValueError(f"Reference {key!r}: unknown type {row['type']!r}")
    if row["result"] == "measurement":
        raise ValueError("Measured datapoints require an original experimental data reader")
    if row["result"] != "theory":
        raise ValueError(f"Reference {key!r}: {row['result']} is not a measured value or prediction")
    measurement = row["measurement"]
    differential = row["type"] == "differential"
    bins = np.asarray(row["bins"] if differential else [0.0, 1.0], dtype=float)
    if bins.ndim != 1 or len(bins) < 2 or not np.isfinite(bins).all() or np.any(np.diff(bins) <= 0.0):
        raise ValueError(f"Reference {key!r}: invalid bin edges")
    shape = (len(bins) - 1,) if differential else ()
    values = np.asarray(measurement["value"], dtype=float)
    if values.shape != shape or not np.isfinite(values).all() or np.any(values < 0.0):
        raise ValueError(f"Reference {key!r}: invalid measurement values")
    values = np.atleast_1d(values)
    valid = np.asarray(row.get("valid", np.ones(len(values), dtype=bool)), dtype=bool)
    if valid.shape != values.shape:
        raise ValueError(f"Reference {key!r}: invalid bin mask")
    units = {None: "1", 0: "b", 3: "mb", 6: "ub", 9: "nb", 12: "pb", 15: "fb"}
    if row["unit"] not in units:
        raise ValueError(f"Reference {key!r}: invalid cross section unit")
    source = row.get("source") or row.get("bibtex") or key
    data = {
        "x": edge2centerbins(bins),
        "bins": bins,
        "binwidth": bins2binwidth(bins),
        "xlim": np.array([bins[0], bins[-1]]),
        "y": values,
        "valid": valid,
        "unit": units[row["unit"]],
        "source": source,
        "channel": row.get("channel"),
        "cuts": row.get("cuts"),
    }
    sources = []
    for name, category, correlation in (
        ("stat", "statistical", "uncorrelated"),
        ("syst", "systematic", "collective"),
        ("lumi", "systematic", "collective"),
        ("total", "combined", "uncorrelated"),
        ("extrap", "systematic", "collective"),
    ):
        if measurement.get(name) is None:
            continue
        errors = np.asarray(measurement[name], dtype=float)
        if errors.shape != (2, *shape) or not np.isfinite(errors).all() or np.any(errors < 0.0):
            raise ValueError(f"Reference {key!r}: invalid {name} uncertainties")
        sources.append(uncertainty_source(
            name, np.atleast_1d(errors[1]), down=np.atleast_1d(errors[0]), category=category,
            correlation=row.get("correlations", {}).get(name, correlation),
            provenance=f"Quoted {name} uncertainty in {source}",
        ))
    if not sources:
        data["reference"] = f"{source}: central values without quoted uncertainties"
    return finalize_uncertainties(data, sources)


# Read one reference ID using the common integrated or differential format
@cache
def read(filename, file_filter=None, rebin_factor=None):
    if rebin_factor is not None:
        raise ValueError("Reference reader: rebinning is not supported")
    if not file_filter:
        raise ValueError("Reference reader: file_filter is required")
    reference = load_json_file(filename)
    if file_filter not in reference:
        raise ValueError(f"Reference reader: unknown reference ID {file_filter!r}")
    return histogram(reference[file_filter], file_filter)
