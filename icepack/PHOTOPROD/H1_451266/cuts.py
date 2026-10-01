# H1 quasi-elastic psi(2S) positron-proton selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

from .._common import dissociative
from .._common.exclusive_jpsi import unique_momentum

cut_param = {"W": [40.0, 160.0], "Q2_MAX": 4.0, "Z_MIN": 0.95}


# Select elastic and proton-dissociative events with the published invariant elasticity
# [REFERENCE: H1, arXiv:hep-ex/9711012, Sections 2 and 6]
def cut_func(event):
    values = dissociative.kinematics(event)
    proton = unique_momentum(event, 2212, True)
    photon = unique_momentum(event, 11, True) - unique_momentum(event, 11, False)
    z = proton.dot4(obs.proj_central_system(event)) / proton.dot4(photon)
    low, high = event.cut_param["W"]
    return low < values["w"] < high and 0 <= values["q2"] < event.cut_param["Q2_MAX"] and z > event.cut_param["Z_MIN"]
