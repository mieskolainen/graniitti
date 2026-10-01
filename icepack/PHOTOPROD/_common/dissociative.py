# HERA proton dissociation reconstruction and fiducial cuts
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs
from core.io.cache import cache

from .exclusive_jpsi import unique_momentum


# Reconstruct the full target system, including every forward fragmentation daughter
@cache
def kinematics(event):
    photon = unique_momentum(event, 11, True) - unique_momentum(event, 11, False)
    proton = unique_momentum(event, 2212, True)
    central = obs.proj_central_system(event)
    target = proton + photon - central
    w2 = (photon + proton).m2
    return {
        "q2": -photon.m2,
        "w": max(w2, 0.0) ** 0.5,
        "abs_t": abs((photon - central).m2),
        "my2": target.m2,
        "mp2": proton.m2,
        "mass": central.m,
    }


# Apply the published proton dissociation phase space
def accepted(event):
    values = kinematics(event)
    w_min, w_max = event.cut_param["W"]
    m_min, m_max = event.cut_param["M"]
    return (
        w_min < values["w"] < w_max
        and 0.0 <= values["q2"] < event.cut_param["Q2_MAX"]
        and values["abs_t"] < event.cut_param["ABS_T_MAX"]
        and values["mp2"] < values["my2"] < min(
            event.cut_param.get("MY_MAX", float("inf")) ** 2,
            event.cut_param.get("MY2_W2_MAX", float("inf")) * values["w"] ** 2,
        )
        and m_min < values["mass"] < m_max
    )


# Compute the absolute momentum transfer at the proton vertex
def abs_t(event):
    return kinematics(event)["abs_t"]


# Compute the photon-proton centre-of-mass energy
def w(event):
    return kinematics(event)["w"]
