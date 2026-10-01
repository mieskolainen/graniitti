# CMS 13 TeV pion-pair 0.35 to 0.65 GeV event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.cuts import cut_func as dphi_acceptance
from .._common.cuts import fiducial_cuts

cut_param = {}


# Apply the missing Delta phi bin veto and selected pion pair mass window
def cut_func(event):
    return bool(fiducial_cuts(event, (0.35, 0.65)) and dphi_acceptance(event))
