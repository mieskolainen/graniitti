# ZEUS unresolved Upsilon positron-proton selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common import exclusive_jpsi

cut_param = {"W": [80.0, 160.0], "Q2_MAX": 1.0}


# Select the published elastic ep region for all three Upsilon dimuon states
# [REFERENCE: ZEUS, arXiv:hep-ex/9807020, Table 1]
def cut_func(event):
    values = exclusive_jpsi.kinematics(event)
    low, high = event.cut_param["W"]
    return low < values["w"] < high and 0 <= values["q2"] < event.cut_param["Q2_MAX"]
