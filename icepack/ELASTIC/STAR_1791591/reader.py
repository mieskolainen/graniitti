# HEPData 1791591 reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.stats.hist import bins2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import (
    error_pairs,
    print_rebin,
    read_table,
    rebin_common,
    split_luminosity,
    table_bins,
    table_values,
)


# Read the STAR elastic differential cross section
def read(filename: str, rebin_factor: int = None) -> dict:
    """Read STAR elastic HEPData 1791591 JSON content"""
    table = read_table(filename)
    values = table_values(table)
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    binwidth = bins2binwidth(bins)
    xlim = np.array([bins[0], bins[-1]])
    y = np.asarray([value["value"] for value in values], dtype=float)
    statistical = error_pairs(values, "Statistical uncertainties")
    systematic = error_pairs(values, "Full systematic uncertainties")
    ds = {
        "x": x,
        "bins": bins,
        "binwidth": binwidth,
        "xlim": xlim,
        "y": y,
    }
    sources = [
        uncertainty_source(
            "statistical", statistical[:, 0], down=statistical[:, 1], category="statistical",
            correlation="uncorrelated", provenance="HEPData statistical columns",
        ),
        uncertainty_source(
            "systematic", systematic[:, 0], down=systematic[:, 1], category="systematic",
            correlation="uncorrelated", provenance="HEPData full systematic columns",
        ),
    ]
    # [REFERENCE: STAR, arXiv:2003.12136, Sec. VI and Table I, luminosity is a common cross-section scale]
    finalize_uncertainties(ds, sources)
    split_luminosity(ds, read_table(Path(filename).with_name("Cross-Sections.json")), "SysUnc luminosity", "systematic")

    if rebin_factor is not None:
        print_rebin(__name__, rebin_factor)
        ds = rebin_common(ds=ds, rebin_factor=rebin_factor)

    return ds
