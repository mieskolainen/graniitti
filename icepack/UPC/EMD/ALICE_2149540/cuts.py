# Inclusive EMD cascade event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis.forward import Acceptance

cut_param = {"emd_side": 1}


# Select the neutron-emitting detector side using final particles
def cut_func(event):
    side = event.cut_param.get("emd_side", 1)
    return bool(Acceptance(pid=(2112,), side=side).particles(event))
