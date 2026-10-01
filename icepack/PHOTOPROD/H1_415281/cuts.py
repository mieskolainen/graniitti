# HERA photoproduction comparison definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from icepack.PHOTOPROD._common.exclusive_jpsi import accepted as cut_func  # noqa: F401

cut_param = {'W': [40.0, 80.0], 'Q2_MAX': 0.5, 'M': [0.52, 1.17], 'ABS_T_MAX': 0.5}
