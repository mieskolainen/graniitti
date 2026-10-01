# ALICE pPb exclusive dimuon upper rapidity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from . import cuts as common_cuts

cut_param = dict(common_cuts.cut_param, rap=[3.25, 4.0])
cut_func = common_cuts.cut_func
