# CMS 13 TeV tensor-region pion-pair event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.cuts import fiducial_cuts

cut_param = {}


# Apply the published CMS tensor region pion pair fiducial selection
def cut_func(event):
    return fiducial_cuts(event, (1.80, 2.20))
