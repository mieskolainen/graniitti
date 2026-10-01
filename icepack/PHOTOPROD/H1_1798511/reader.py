# H1 elastic pion pair photoproduction from original HEPData
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
from core.io.cache import cache
from core.stats.hist import bins2binwidth
from core.stats.uncertainty import finalize_uncertainties

from icepack._common.hepdata import read_table, select_rows
from icepack.PHOTOPROD._common.diss_reader import rho_sources
from icepack.PHOTOPROD._common.hera_reader import photon_flux, read_bins


# Select the elastic mass spectrum and propagate its published signed variations
@cache
def read(filename, file_filter=None, rebin_factor=None, hist=None):
    if rebin_factor is not None:
        raise ValueError("ReadH1_1798511: rebinning is not supported")
    if file_filter not in {None, "elastic"}:
        raise ValueError(f"ReadH1_1798511: unsupported channel {file_filter!r}")
    hist = {} if hist is None else hist
    path = Path(filename)
    table = read_table(path)
    headers = [item["name"] for item in table["headers"][:table["x_count"]]]
    mass = next(i for i, name in enumerate(headers) if name.startswith(r"$m_{\pi\pi}$"))
    target = next(i for i, name in enumerate(headers) if name.startswith("$m_Y$"))
    global_bin = next(i for i, name in enumerate(headers) if name.startswith("globalBinNumber"))
    rows = [row for row in select_rows(table, hist.get("rows"))["values"] if "value" in row["x"][target]]
    x, binedges, bins = read_bins({"values": [{"x": [row["x"][mass]]} for row in rows]})
    mc_scale = 1.0 / photon_flux(table, rows)
    for i, name in enumerate(headers):
        if name.startswith("$t$"):
            mc_scale /= np.asarray([float(row["x"][i]["high"]) - float(row["x"][i]["low"]) for row in rows])
    covariance = hist.get("covariance")
    if covariance is None and not any(name.startswith("$t$") for name in headers):
        raise ValueError("H1 mass spectrum requires its HEPData statistical correlation table")
    dataset = {
        "x": x, "binedges": binedges, "bins": bins, "binwidth": bins2binwidth(bins),
        "xlim": bins[[0, -1]], "y": np.asarray([float(row["y"][0]["value"]) for row in rows]),
        "valid": np.ones(len(rows), dtype=bool), "mc_scale": mc_scale,
    }
    return finalize_uncertainties(dataset, rho_sources(rows, None if covariance is None else path.parent / covariance, global_bin))
