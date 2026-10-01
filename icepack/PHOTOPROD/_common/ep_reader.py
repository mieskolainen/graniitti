# Integrated HERA positron-proton vector meson cross sections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import re
from pathlib import Path

import numpy as np
from core.stats.uncertainty import finalize_uncertainties

from icepack._common.hepdata import number, read_table, resolve_group, select_rows, table_values
from icepack.PHOTOPROD._common.hera_reader import read_uncertainties


# Read one published ep cross section without a photon flux conversion
def read(filename, file_filter=None, hist=None, rebin_factor=None, **_kwargs):
    if rebin_factor is not None:
        raise ValueError("Integrated HERA cross sections cannot be rebinned")
    table = select_rows(read_table(filename, dimensions=1), (hist or {}).get("rows"))
    group = resolve_group(table, file_filter, headers=True)
    values = table_values(table, group)
    if len(values) != 1:
        raise ValueError("Select exactly one integrated HERA measurement")
    headers = [header["name"] for header in table["headers"][table["x_count"]:]
               for _ in range(header.get("colspan", 1))]
    unit = re.search(r"\[(PB|NB)\]$", headers[group])
    if unit is None:
        raise ValueError("Integrated HERA measurement requires a cross section unit")
    return finalize_uncertainties({
        "x": np.array([0.5]), "bins": np.array([0.0, 1.0]), "binwidth": np.array([1.0]),
        "xlim": np.array([0.0, 1.0]), "y": np.array([number(values[0]["value"])]),
        "unit": unit[1].lower(), "source": table["doi"], "group": group,
    }, read_uncertainties(values, f"{Path(filename).parent.name}:group{group}"))
