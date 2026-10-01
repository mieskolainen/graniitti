# ALICE rho rapidity folding of original HEPData values and uncertainties
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy

import numpy as np
from core.stats.hist import bins2binwidth, edge2centerbins
from core.stats.uncertainty import finalize_uncertainties, linear_source

from icepack.UPC._common import upc_reader as reader


# Fold the independent absolute rapidity bins and propagate the HEPData errors
def read(filename, file_filter=None, rebin_factor=None, **_kwargs):
    ds = copy.deepcopy(reader.read(filename, file_filter, rebin_factor, hist={"rows": [2, 3, 4]}))
    ds["bins"][0] = ds["binedges"][0, 0] = 0.0
    ds["binwidth"] = bins2binwidth(ds["bins"])
    ds["x"] = edge2centerbins(ds["bins"])
    ds["xlim"] = ds["bins"][[0, -1]]
    # Folding both rapidity signs gives d sigma / d |y| = 2 d sigma / dy
    ds["y"] *= 2.0
    return finalize_uncertainties(ds, [linear_source(source, 2.0 * np.eye(len(ds["y"])))
                                      for source in ds["uncertainties"]])
