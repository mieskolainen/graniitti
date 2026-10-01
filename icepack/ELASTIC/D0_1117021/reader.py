# D0 elastic HEPData 1117021 reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, print_rebin, read_table, rebin_common, table_bins, table_values


# Read published D0 point errors and correlated luminosity uncertainty
@cache
def read(filename: str, rebin_factor: int = None) -> dict:
    """Read D0 elastic HEPData 1117021 JSON content"""
    table = read_table(filename)
    values = table_values(table)
    x, edges = table_bins(table)
    bins = doublet2linear(edges)
    binwidth = bins2binwidth(bins)
    xlim = np.array([bins[0], bins[-1]])
    y = np.asarray([value["value"] for value in values], dtype=float)
    y_err_data = error_pairs(values, "error")
    y_err_lumi = error_pairs(values, "sys,luminosity uncertainty")
    ds = {
        "x": x,
        "bins": bins,
        "binwidth": binwidth,
        "xlim": xlim,
        "y": y,
    }
    # [REFERENCE: HEPData 1117021 Table 1 description and luminosity columns]
    finalize_uncertainties(
        ds,
        [
            uncertainty_source(
                "combined_data",
                y_err_data[:, 0],
                down=y_err_data[:, 1],
                category="combined",
                correlation="uncorrelated",
                provenance="HEPData combines statistical and non-luminosity systematic errors",
            ),
            uncertainty_source(
                "luminosity",
                y_err_lumi[:, 0],
                down=y_err_lumi[:, 1],
                category="systematic",
                correlation="collective",
                scope="D0_1117021:luminosity",
                effect="multiplicative",
                provenance="HEPData luminosity percentage columns",
            ),
        ],
    )

    if rebin_factor is not None:
        print_rebin(__name__, rebin_factor)
        ds = rebin_common(ds=ds, rebin_factor=rebin_factor)

    return ds
