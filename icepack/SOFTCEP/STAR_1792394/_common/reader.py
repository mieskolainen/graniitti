# HEPData 1792394 reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from html import unescape
from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import (
    error_pairs,
    fill_bins,
    print_rebin,
    read_table,
    rebin_common,
    select_rows,
    table_bins,
    table_values,
)


# Read the Figure 13 proton cuts from the original HEPData description
def read_cuts(filename):
    path = Path(__file__).resolve().parents[4] / "HEPData/CEP/HEPData-ins1792394-v1-json" / filename
    description = unescape(read_table(path)["description"].split("\n")[0])
    bounds = re.findall(r"([<>])\s*([0-9.]+)", description)
    if len(bounds) != 2 or bounds[0][0] != "<":
        raise ValueError(f"STAR Figure 13 requires azimuth and transverse-momentum cuts: {path}")
    return {"dPhi": float(bounds[0][1]), "dPt": float(bounds[1][1]), "above": bounds[1][0] == ">"}


# Select one fiducial cross section from the original multichannel HEPData table
def read_integrated(filename, reaction, rows):
    table = read_table(filename)
    if rows is None or len(rows) != 1:
        raise ValueError("STAR Table 2 requires one reaction and one azimuthal interval")
    values = table_values(select_rows(table, rows), reaction)
    sources = []
    for name, label, category in (("stat", "stat.", "statistical"), ("syst", "syst.", "systematic")):
        pair = error_pairs(values, label)
        sources.append(uncertainty_source(
            name, pair[:, 0], down=pair[:, 1], category=category, correlation="uncorrelated",
            provenance=f"Original HEPData Table 2 {label} errors",
        ))
    return finalize_uncertainties({
        "x": np.array([0.5]), "bins": np.array([0.0, 1.0]), "binwidth": np.array([1.0]),
        "xlim": np.array([0.0, 1.0]), "y": np.array([float(values[0]["value"])]), "unit": "nb",
        "source": table["doi"],
    }, sources)


# Read every published bin, retaining gaps as unmeasured intervals
@cache
def read(filename: str, rebin_factor: int = None, file_filter=None, hist=None) -> dict:
    """Read STAR HEPData 1792394 JSON content"""
    if Path(filename).name == "Table2.json":
        if rebin_factor is not None:
            raise ValueError("STAR integrated cross sections cannot be rebinned")
        return read_integrated(filename, file_filter, None if hist is None else hist.get("rows"))
    table = read_table(filename)
    values = table_values(table)
    centers, edges = table_bins(table)
    yerr_stat = error_pairs(values, "stat.")
    yerr_syst = error_pairs(values, "syst. (luminosity)")
    yerr_sexp = error_pairs(values, "syst. (experimental)")
    ds = {"x": centers, "binedges": edges, "y": np.asarray([float(value["value"]) for value in values])}
    table_scope = Path(filename).stem
    # [REFERENCE: HEPData 1792394 differential-table correlation description]
    finalize_uncertainties(
        ds,
        [
            uncertainty_source(
                "statistical",
                yerr_stat[:, 0],
                down=yerr_stat[:, 1],
                category="statistical",
                correlation="uncorrelated",
                provenance="HEPData statistical columns",
            ),
            uncertainty_source(
                "luminosity",
                yerr_syst[:, 0],
                down=yerr_syst[:, 1],
                category="systematic",
                correlation="collective",
                scope="STAR_1792394:luminosity",
                effect="multiplicative",
                provenance="HEPData states luminosity is fully correlated over all points",
            ),
            uncertainty_source(
                "experimental",
                yerr_sexp[:, 0],
                down=yerr_sexp[:, 1],
                category="systematic",
                correlation="collective",
                scope=f"STAR_1792394:{table_scope}:experimental",
                provenance="HEPData states the experimental systematic is a collective bin variation",
            ),
        ],
    )

    ds = fill_bins(ds)
    if rebin_factor is not None:
        print_rebin(__name__, rebin_factor)
        ds = rebin_common(ds=ds, rebin_factor=rebin_factor)

    return ds
