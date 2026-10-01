# TOTEM elastic HEPData reader with the original common luminosity uncertainty
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
from pathlib import Path

from icepack._common.hepdata import read_table, rebin_common, split_luminosity
from icepack.ELASTIC._common import elastic_reader


# Table 4 gives the optical-point luminosity error as its third systematic component
# [REFERENCE: HEPData 10.17182/hepdata.66456.v1/t4, systematic source order in the table description]
def read(filename, rebin_factor=None):
    ds = copy.deepcopy(elastic_reader.read(filename))
    split_luminosity(ds, read_table(Path(filename).with_name("Table4.json")), "sys_3", "systematic")
    return ds if rebin_factor is None else rebin_common(ds, rebin_factor)
