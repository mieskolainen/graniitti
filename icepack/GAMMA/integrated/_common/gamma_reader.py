# Original HEPData integrated dilepton cross section reader
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re

import numpy as np
from core.io.cache import cache
from core.stats.uncertainty import finalize_uncertainties, uncertainty_source

from icepack._common.hepdata import error_pairs, read_table


# Read exactly one integrated measurement without copying or supplementing its errors
@cache
def read(filename, rebin_factor=None):
    if rebin_factor is not None:
        raise ValueError("Integrated HEPData measurements cannot be rebinned")
    table = read_table(filename)
    if len(table["values"]) != 1 or len(table["values"][0]["y"]) != 1:
        raise ValueError("Integrated HEPData reader requires one measured value")
    row = table["values"][0]["y"][0]
    value = float(row["value"])
    sources = []
    for index, entry in enumerate(row["errors"]):
        label = entry.get("label", "combined")
        name = "lumi" if "lumi" in label.lower() else "stat" if label.lower() == "stat" else "syst"
        category = "statistical" if name == "stat" else "combined" if label == "combined" else "systematic"
        up, down = error_pairs([row], index)[0]
        sources.append(uncertainty_source(
            name, [abs(up)], down=[abs(down)], category=category, correlation="uncorrelated",
            effect="multiplicative" if name == "lumi" else "additive", provenance=f"Original HEPData {label} field",
        ))
    unit = re.search(r"\[([^]]+)\]", table["headers"][-1]["name"])[1].lower()
    return finalize_uncertainties({
        "x": np.array([0.5]), "bins": np.array([0.0, 1.0]), "binwidth": np.array([1.0]),
        "xlim": np.array([0.0, 1.0]), "y": np.array([value]), "unit": unit, "source": table["doi"],
    }, sources)
