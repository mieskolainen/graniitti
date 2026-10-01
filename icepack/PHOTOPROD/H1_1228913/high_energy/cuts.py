# H1 high-energy elastic J/psi photoproduction selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import exclusive_jpsi

cut_param = {"W": [40.0, 110.0], "Q2_MAX": 2.5, "ABS_T_MAX": 1.2}


# Apply the published H1 high-energy elastic phase space
def cut_func(event):
    return exclusive_jpsi.accepted(event)
