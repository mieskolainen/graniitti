# LHCb HEPData 1277076 exclusive charmonium reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
from core.io.cache import cache
from core.stats.hist import doublet2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, norm, uncertainty_source

from icepack._common.hepdata import error_pairs, read_table, table_bins, table_values

STATE_COLUMNS = {
    "J/psi": "J/PSI",
    "psi(2S)": "PSI(2S)",
}

REGIONS = {"fiducial"}


# Compute the published integrated cross sections and uncertainties by meson
@cache
def read_integrated(filename):
    table = read_table(filename)
    values = table_values(table)
    output = {}
    for row, value in zip(table["values"], values, strict=True):
        reaction = row["x"][0]["value"]
        state = next(state for state, name in STATE_COLUMNS.items() if f" {name} " in reaction)
        stat = float(np.mean(np.abs(error_pairs([value], "stat"))))
        syst = float(np.mean(np.abs(error_pairs([value], "sys"))))
        output[state] = {"value": float(value["value"]), "stat": stat, "syst": syst,
                         "total": float(norm(np.asarray([stat, syst]))), "unit": "pb"}
    return output


# Read the fiducial charmonium rapidity distribution with the original correlated totals
@cache
def _read_rapidity(filename, state, region, rebin_factor):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_1277076: rebinning is not supported")
    if state not in STATE_COLUMNS:
        raise ValueError(f"ReadHEPData_1277076: unsupported state {state!r}")
    table = read_table(filename)
    values = table_values(table, f"P P--> {STATE_COLUMNS[state]} <MU+ MU-> P")
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    statistical = error_pairs(values, "error")
    sources = [uncertainty_source(
        "counting_statistical", statistical[:, 0], down=statistical[:, 1],
        category="statistical", correlation="uncorrelated",
        provenance="HEPData Table 2 uncorrelated statistical errors",
    )]
    for category, label in (
        ("statistical", f"total correlated statistical uncertainty ({STATE_COLUMNS[state]})"),
        ("systematic", "total correlated systematic uncertainty"),
    ):
        pair = error_pairs(values, f"sys,{label}")
        sources.append(uncertainty_source(
            f"correlated_{category}", pair[:, 0], down=pair[:, 1], category=category, correlation="collective",
            scope=f"LHCb_1277076:{state}:{category}", effect="multiplicative",
            provenance=f"HEPData Table 2 {label}, no cross-state covariance supplied",
        ))
    return finalize_uncertainties({
        "state": state, "region": region, "x": x, "binedges": edges, "bins": bins,
        "binwidth": doublet2binwidth(edges), "xlim": bins[[0, -1]],
        "y": np.asarray([value["value"] for value in values], dtype=float),
        "y_err_uncorrelated": np.mean(np.abs(statistical), axis=1),
    }, sources)


# Resolve the state explicitly or from one iceplot dataset set
def _resolve_state(state, dataset):
    if state is not None:
        return state
    if dataset is None or "state" not in dataset:
        raise ValueError("ReadHEPData_1277076: dataset state is required")
    return str(dataset["state"])


# Resolve and validate one explicitly selected measurement region
def _resolve_region(region, dataset):
    selected = region
    if selected is None:
        if dataset is None or "region" not in dataset:
            raise ValueError("ReadHEPData_1277076: dataset region is required")
        selected = dataset["region"]
    selected = str(selected)
    if selected not in REGIONS:
        raise ValueError(f"ReadHEPData_1277076: unsupported region {selected!r}")
    return selected


# Read one explicitly selected fiducial measurement
def read(filename, rebin_factor=None, state=None, region=None, dataset=None):
    return _read_rapidity(
        filename,
        _resolve_state(state, dataset),
        _resolve_region(region, dataset),
        rebin_factor,
    )
