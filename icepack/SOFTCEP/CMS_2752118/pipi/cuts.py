# CMS 13 TeV pion-pair event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.cuts import fiducial_cuts

cut_param = {}


# Apply the published CMS pion pair mass and fiducial selection
def cut_func(event):
    return fiducial_cuts(event, (0.30, 2.50))
