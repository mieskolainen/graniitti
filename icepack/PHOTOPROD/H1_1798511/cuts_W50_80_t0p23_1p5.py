# H1 elastic dipion mass spectrum fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from functools import partial

from icepack.PHOTOPROD._common.exclusive_jpsi import accepted

# Reuse cut-independent kinematics between the H1 fiducial selections
cut_param = {}
cut_func = partial(accepted, cuts={'W': [50, 80], 'Q2_MAX': 2.5, 'M': [0.5, 1.2], 'ABS_T_MAX': 1.5, 'ABS_T_MIN': 0.23})
