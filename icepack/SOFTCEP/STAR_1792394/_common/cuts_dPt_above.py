# STAR 200 GeV proton transverse-momentum difference above the published boundary
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .reader import read_cuts
from .selections import dpt_cut as cut_func

__all__ = ["cut_func", "cut_param"]

cut_param = read_cuts("Figure13(right).json")
