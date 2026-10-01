# STAR 510 GeV mass and proton azimuth selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs


# Select one published mass interval and proton azimuth region
def parameters(mass=None, phi=None):
    bounds = [[0.0, 1.0], [1.0, 1.5], [1.5, None]][mass] if mass is not None else [0.0, None]
    return {"M": bounds,
            "phi": phi, "dPhi": [[0.0, 90.0], [90.0, 180.0]][phi] if phi is not None else None}


# Select the central pair and tagged proton azimuth in their published intervals
def cut_func(event):
    mass = obs.proj_1D_M(event)
    phi = obs.proj_1D_dPhi_pp(event)
    low, high = event.cut_param["M"]
    side = event.cut_param["phi"]
    passed = side is None or (phi < event.cut_param["dPhi"][1] if side == 0
                              else phi > event.cut_param["dPhi"][0])
    return low <= mass and (high is None or mass < high) and passed


cut_param = parameters()
