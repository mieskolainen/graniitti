# H1 proton dissociation fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from ..._common import dissociative

cut_param = {"W": [40.0, 110.0], "Q2_MAX": 2.5, "ABS_T_MAX": 8.0, "MY_MAX": 10.0, "M": [3.0, 3.2]}


# Apply the published H1 proton dissociation cuts
def cut_func(event):
    return dissociative.accepted(event)
