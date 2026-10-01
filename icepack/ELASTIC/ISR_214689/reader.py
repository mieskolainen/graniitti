# ISR elastic proton proton and antiproton proton HEPData reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth, doublet2linear
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, print_rebin, read_table, rebin_common, table_bins, table_values


# Read one ISR small angle elastic differential cross section table
@cache
def read(filename: str, rebin_factor: int | None = None) -> dict:
    table = read_table(filename)
    values = table_values(table)
    x, edges = table_bins(table)
    x, edges = 1.0e-3 * x, 1.0e-3 * edges
    bins = doublet2linear(edges)
    y = np.asarray([value["value"] for value in values], dtype=float)
    errors = error_pairs(values, "error")
    dataset = {
        "x": x,
        "bins": bins,
        "binwidth": bins2binwidth(bins),
        "xlim": np.array([bins[0], bins[-1]]),
        "y": y,
    }
    # [REFERENCE: HEPData 214689 Tables 2 to 7 error columns]
    finalize_uncertainties(
        dataset,
        [
            uncertainty_source(
                "published",
                errors[:, 0],
                down=errors[:, 1],
                category="combined",
                correlation="uncorrelated",
                provenance="HEPData error columns without a source decomposition",
            ),
        ],
    )

    if rebin_factor is not None:
        print_rebin(__name__, rebin_factor)
        dataset = rebin_common(ds=dataset, rebin_factor=rebin_factor)
    return dataset
