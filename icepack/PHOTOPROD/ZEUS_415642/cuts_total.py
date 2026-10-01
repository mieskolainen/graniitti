# ZEUS elastic phi cross section extrapolated to zero momentum transfer
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.exclusive_jpsi import accepted as cut_func  # noqa: F401
from .cuts import cut_param

cut_param = {**cut_param, "ABS_T_MIN": 0.0}
