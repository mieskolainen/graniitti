# STAR coherent rho rapidity and neutron selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.analysis.forward import Acceptance

cut_param = {"abs_rap_max": 1.0, "neutrons": "XnXn"}


# Select the physical forward neutron multiplicities and central meson rapidity
def cut_func(event):
    if abs(obs.proj_1D_Rap(event)) >= event.cut_param["abs_rap_max"]:
        return False
    counts = [len(Acceptance(pid=(2112,), beam=beam).particles(event)) for beam in (1, 2)]
    return all(n == 1 for n in counts) if event.cut_param["neutrons"] == "1n1n" else all(n > 0 for n in counts)
