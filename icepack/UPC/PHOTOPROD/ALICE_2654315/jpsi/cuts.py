# ALICE pPb J/psi rapidity and transverse momentum selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.analysis.forward import neutron_class

cut_param = {"rap": [2.5, 4.0], "pt_max": 3.0}


# Apply the published J/psi rapidity and transverse momentum selection
def cut_func(event):
    rapidity = obs.proj_1D_Rap(event)
    lower, upper = event.cut_param["rap"]
    return neutron_class(event) in ("0n0n", "Xn0n") and lower < rapidity < upper and obs.proj_1D_Pt(event) < event.cut_param["pt_max"]
