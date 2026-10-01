# H1 low-energy elastic J/psi photoproduction selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import exclusive_jpsi

cut_param = {"W": [25.0, 80.0], "Q2_MAX": 2.5, "ABS_T_MAX": 1.2}


# Apply the published H1 low-energy elastic phase space
def cut_func(event):
    return exclusive_jpsi.accepted(event)
