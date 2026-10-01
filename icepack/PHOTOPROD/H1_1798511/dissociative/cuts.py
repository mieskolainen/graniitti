# H1 proton dissociation fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import dissociative

cut_param = {"W": [20.0, 80.0], "Q2_MAX": 2.5, "ABS_T_MAX": 1.5, "MY_MAX": 10.0, "M": [0.2792, 1.53]}


# Apply the published H1 proton dissociation cuts
def cut_func(event):
    return dissociative.accepted(event)
