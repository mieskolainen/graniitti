# ATLAS PbPb dimuon mass spectrum selection at central rapidity
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from . import cuts as common_cuts

cut_param = dict(common_cuts.cut_param, abs_rap_max=0.8)


# Apply the fiducial cuts and the published central rapidity interval
def cut_func(event):
    system = common_cuts.dimuon(event)
    return (
        system is not None
        and common_cuts.fiducial_cut(event)
        and abs(system.rapidity) < event.cut_param["abs_rap_max"]
    )
