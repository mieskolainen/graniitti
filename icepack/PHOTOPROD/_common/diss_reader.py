# H1 proton dissociation HEPData spectra with published photon fluxes and covariances
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth
from core.stats.uncertainty import finalize_uncertainties, offset_source, uncertainty_source

from icepack._common.hepdata import error_pairs, read_table

from . import hera_reader


# Select the published rho proton dissociation rows using their target mass interval
def rho_rows(table):
    rows = [row for row in table["values"] if "low" in row["x"][1]]
    if not rows or any(
        not np.allclose([float(row["x"][1]["low"]), float(row["x"][1]["high"])], [1.0, 10.0]) for row in rows
    ):
        raise ValueError("H1 dissociation: expected the published 1 < M_Y < 10 GeV rows")
    return rows


# Construct the rho statistical covariance from the complete published correlation table
def rho_stat(rows, filename, bin_index=2):
    indices = {int(row["x"][bin_index]["value"]): i for i, row in enumerate(rows)}
    stat = np.abs(error_pairs([row["y"][0] for row in rows], 0)[:, 0])
    if filename is None:
        return uncertainty_source(
            "statistical", stat, category="statistical", correlation="uncorrelated",
            provenance="Official H1 statistical errors, no statistical correlation table supplied",
        )
    correlation = np.full((len(rows), len(rows)), np.nan)
    for row in read_table(filename)["values"]:
        first, second = (int(item["value"]) for item in row["x"])
        if first in indices and second in indices:
            i, j = indices[first], indices[second]
            correlation[i, j] = correlation[j, i] = float(row["y"][0]["value"])
    if not np.isfinite(correlation).all() or not np.allclose(np.diag(correlation), 1.0):
        raise ValueError("H1 dissociation: incomplete statistical correlations")
    return uncertainty_source(
        "statistical",
        correlation * np.outer(stat, stat),
        category="statistical",
        correlation="covariance",
        provenance="Official H1 statistical correlation table",
    )


# Retain every signed correlated rho systematic variation
def rho_sources(rows, filename, bin_index=2):
    sources = [rho_stat(rows, filename, bin_index)]
    errors = [row["y"][0]["errors"] for row in rows]
    labels = [error["label"] for error in errors[0]]
    if any([error["label"] for error in entries] != labels for entries in errors):
        raise ValueError("H1 dissociation: inconsistent systematic variations")
    for i, label in enumerate(labels[1:], start=1):
        up, down = error_pairs([row["y"][0] for row in rows], i).T
        sources.append(offset_source(label, up, down, f"Official H1 {label} offset variation"))
    return sources


# Read the original gamma p spectra and their published uncertainties
# [REFERENCE: H1 arXiv:1304.5162, Tables 4-6 and H1 arXiv:2005.14471, Tables 7, 10]
@cache
def _read(filename, file_filter, covariance, flux):
    if file_filter not in {None, "W", "pd_t", "pd_W"}:
        raise ValueError(f"H1 dissociation: unknown selection {file_filter!r}")
    if file_filter is None:
        return hera_reader.read(filename, hist={"flux": flux})
    path = Path(filename)
    table = read_table(path)
    rho = file_filter.startswith("pd_")
    rows = rho_rows(table) if rho else table["values"]
    x, binedges, bins = hera_reader.read_bins({"values": rows})
    widths = bins2binwidth(bins)
    values = [row["y"][0] for row in rows]
    dataset = {
        "x": x,
        "binedges": binedges,
        "bins": bins,
        "binwidth": widths,
        "xlim": np.asarray([bins[0], bins[-1]]),
        "y": np.asarray([float(value["value"]) for value in values]),
    }
    dataset["mc_scale"] = 1.0 / hera_reader.photon_flux(table, rows)
    if rho:
        if covariance is None:
            raise ValueError("H1 dissociation requires its HEPData statistical correlation table")
        sources = rho_sources(rows, path.parent / covariance)
    else:
        sources = hera_reader.read_uncertainties(values, str(path))
    finalize_uncertainties(dataset, sources)
    return dataset


# Cache only the table selection while accepting the standard analysis reader arguments
def read(filename, file_filter=None, rebin_factor=None, hist=None, **_kwargs):
    if rebin_factor is not None:
        raise ValueError("H1 dissociation: rebinning is not supported")
    hist = {} if hist is None else hist
    return _read(filename, file_filter, hist.get("covariance"), hist.get("flux"))
