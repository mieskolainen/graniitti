# ALICE coherent J/psi reader with neutron migration anticorrelations
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
from icepack._common.hepdata import signed_sources
from icepack.UPC._common import upc_reader


# Correlate migration shifts between neutron classes at the same rapidity
def uncertainties(table, values, scope):
    sources = upc_reader.read_uncertainties(values, scope)
    if "anti-correlated" not in table["description"]:
        return sources
    sources = [source for source in sources if source["name"] != "migr"]
    return sources + signed_sources(table, values, "migr", f"{scope.split(':')[0]}:migr")


# Read the original table and its published neutron migration convention
def read(filename, **kwargs):
    return upc_reader.read(filename, uncertainties=uncertainties, **kwargs)
