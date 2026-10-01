# ZEUS elastic J/psi photoproduction selection at 90 to 110 GeV
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common import exclusive_jpsi

cut_param = {"W": [90.0, 110.0], "Q2_MAX": 1.0, "ABS_T_MAX": 1.8}


# Apply the published ZEUS elastic phase space
def cut_func(event):
    return exclusive_jpsi.accepted(event)
