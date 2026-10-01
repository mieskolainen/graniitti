# H1 low energy proton dissociation momentum transfer selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import dissociative

# [REFERENCE: H1 arXiv:1304.5162, Table 1]
cut_param = {"W": [25.0, 80.0], "Q2_MAX": 2.5, "ABS_T_MAX": 5.0, "MY_MAX": 10.0, "M": [3.0, 3.2]}


# Apply the smaller transfer range used in the low energy differential measurement
def cut_func(event):
    return dissociative.accepted(event)
