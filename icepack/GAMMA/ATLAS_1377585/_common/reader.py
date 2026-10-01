# ATLAS HEPData 1377585 shared exclusive dilepton reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, fill_bins, read_table, table_bins, table_values


# Read one measured dilepton fiducial cross section from the original HEPData blocks
def read_integrated(filename, channel):
    table = read_table(filename)
    values = [value for value in table_values(table, channel) if value["value"] != "-"]
    if len(values) != 1:
        raise ValueError("ATLAS Table 1 requires one measured fiducial cross section")
    sources = []
    for name, label, category in (("stat", "stat", "statistical"), ("syst", "sys", "systematic")):
        pair = error_pairs(values, label)
        sources.append(uncertainty_source(
            name, pair[:, 0], down=pair[:, 1], category=category, correlation="uncorrelated",
            provenance=f"Original HEPData Table 1 {label} errors",
        ))
    # The HEPData record abstract specifies pb for the integrated SIG values
    return finalize_uncertainties({
        "x": np.array([0.5]), "bins": np.array([0.0, 1.0]), "binwidth": np.array([1.0]),
        "xlim": np.array([0.0, 1.0]), "y": np.array([float(values[0]["value"])]), "unit": "pb",
        "source": table["doi"],
    }, sources)


# Read one unfolded ATLAS acoplanarity event count distribution
@cache
def read(filename, rebin_factor=None, file_filter=None):
    if rebin_factor is not None:
        raise ValueError("ReadHEPData_1377585: rebinning is not supported")

    if Path(filename).name == "Table1.json":
        return read_integrated(filename, file_filter)
    table = read_table(filename)
    if table["headers"][0]["name"] != "ACO" or table["headers"][1]["name"] != "N":
        raise ValueError("ReadHEPData_1377585: expected acoplanarity event counts")
    values = table_values(table)
    x, edges = table_bins(table)
    stat, syst = error_pairs(values, "stat"), error_pairs(values, "sys")
    ds = {"x": x, "binedges": edges, "y": np.asarray([float(value["value"]) for value in values])}
    # [REFERENCE: HEPData 1377585 Tables 5 and 6]
    finalize_uncertainties(
        ds,
        [
            uncertainty_source(
                "statistical",
                stat[:, 0],
                down=stat[:, 1],
                category="statistical",
                correlation="uncorrelated",
                provenance="HEPData statistical column",
            ),
            uncertainty_source(
                "systematic",
                syst[:, 0],
                down=syst[:, 1],
                category="systematic",
                correlation="uncorrelated",
                provenance="HEPData supplies no covariance matrix",
            ),
        ],
    )
    return fill_bins(ds)
