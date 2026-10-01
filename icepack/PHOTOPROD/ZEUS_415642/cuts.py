# ZEUS elastic phi differential photoproduction selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.exclusive_jpsi import accepted as cut_func  # noqa: F401

# [REFERENCE: arXiv:hep-ex/9601009, Sections 2 and 6.3]
cut_param = {"W": [60.0, 80.0], "Q2_MAX": 4.0, "ABS_T_MIN": 0.1, "ABS_T_MAX": 0.5}
