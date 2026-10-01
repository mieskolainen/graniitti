# ATLAS photon scattering HEPData reader with a shared luminosity uncertainty
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from pathlib import Path

from icepack._common.hepdata import read_table, split_luminosity
from icepack.UPC._common import upc_reader


# Separate the Table 11 luminosity term already included in differential total errors
# [REFERENCE: ATLAS, arXiv:2008.05355, Secs. 7 and 8, HEPData 10.17182/hepdata.95747.v1/t11]
def read(filename, **kwargs):
    ds = copy.deepcopy(upc_reader.read(filename, **kwargs))
    if "normalis" in Path(filename).stem.lower():
        return ds
    integrated = any(source["name"] == "lumi" for source in ds["uncertainties"])
    if integrated:
        ds["uncertainties"] = [source for source in ds["uncertainties"] if source["name"] != "lumi"]
    return split_luminosity(ds, read_table(Path(filename).with_name("Table11.json")), "lumi",
                            None if integrated else "systematic_from_total")
