# ZEUS J/psi helicity angular distribution selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import exclusive_jpsi

# [REFERENCE: ZEUS, arXiv:hep-ex/0201043, Section 8]
cut_param = {"W": [30.0, 170.0], "Q2_MAX": 1.0, "ABS_T_MAX": 1.0}


# Apply the phase space of the acceptance corrected angular distributions
def cut_func(event):
    return exclusive_jpsi.accepted(event)
