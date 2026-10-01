# ALICE pPb exclusive dimuon fiducial selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.analysis.forward import neutron_class

cut_param = {"rap": [2.5, 4.0], "pt_max": 3.0, "mass": [1.0, 2.5]}


# Apply the published dimuon rapidity and transverse momentum selection
def cut_func(event):
    system = obs.proj_central_system(event)
    lower, upper = event.cut_param["rap"]
    mass_min, mass_max = event.cut_param["mass"]
    return (
        neutron_class(event) in ("0n0n", "Xn0n")
        and lower < system.rapidity < upper
        and mass_min < system.m < mass_max
        and system.pt < event.cut_param["pt_max"]
    )
