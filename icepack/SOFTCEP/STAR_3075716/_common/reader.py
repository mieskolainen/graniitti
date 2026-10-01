# Original STAR 510 GeV HEPData spectra and experimental uncertainties
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
from core.io.cache import cache
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, fill_bins, read_table, rebin_common, table_bins, table_values


# Read published bins and preserve unmeasured intervals as masked gaps
@cache
def read(filename, rebin_factor=None):
    table = read_table(filename)
    values = table_values(table)
    centers, edges = table_bins(table)
    y = np.asarray([float(value["value"]) for value in values])
    sources = []
    for name, label, category in (("statistical", "Stat. error", "statistical"),
                                  ("experimental", "Syst. error", "systematic")):
        errors = error_pairs(values, label)
        sources.append(uncertainty_source(
            name, errors[:, 0], down=errors[:, 1], category=category, correlation="uncorrelated",
            provenance="Original HEPData errors, bin correlations unpublished and assumed uncorrelated",
        ))
    data = fill_bins(finalize_uncertainties({
        "x": centers, "y": y, "binedges": edges, "source": table["doi"],
        "unit": "pb" if "[pb" in table["headers"][1]["name"] else "nb",
    }, sources))
    return data if rebin_factor is None else rebin_common(data, rebin_factor)
