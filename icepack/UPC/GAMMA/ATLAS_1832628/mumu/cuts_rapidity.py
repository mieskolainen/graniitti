# ATLAS PbPb dimuon rapidity spectrum selection at low mass
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from . import cuts as common_cuts

cut_param = dict(common_cuts.cut_param, mass_max=20.0)
cut_func = common_cuts.fiducial_cut
