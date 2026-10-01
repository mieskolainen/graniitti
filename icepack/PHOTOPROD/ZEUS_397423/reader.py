# HERA photoproduction comparison definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

from icepack.PHOTOPROD._common import hera_reader
from icepack.PHOTOPROD._common.flux import rho_flux

from .cuts import cut_param


# Preserve measured cross sections and unfold the experimental VMD photon flux
def read(filename, **kwargs):
    data = dict(hera_reader.read(filename, **kwargs))
    data['mc_scale'] = 1.0 / rho_flux(str(Path(__file__).with_name('gencard.json')), cut_param)
    return data
