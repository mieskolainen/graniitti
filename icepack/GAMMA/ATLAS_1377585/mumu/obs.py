# ATLAS 7 TeV exclusive dimuon histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
#
# [REFERENCE: ATLAS, Phys. Lett. B 749 (2015) 242, arXiv:1506.07098]

import math

from core.analysis import obs
from core.io.cache import cache


# Compute the two stable central muons ordered by transverse momentum
@cache
def selected_muons(event):
    records = [
        record
        for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 13
    ]
    if len(records) != 2:
        raise ValueError(f"selected_muons: expected two final muons, found {len(records)}")
    return sorted(
        (record["p4"] for record in records), key=lambda momentum: momentum.pt, reverse=True
    )


# Compute ATLAS acoplanarity
@cache
def proj_acoplanarity(event):
    first, second = selected_muons(event)
    value = 1.0 - first.abs_delta_phi(second) / math.pi
    return min(1.0, max(0.0, value))


obs_acoplanarity = {
    "tag": "acoplanarity",
    "xlim": None,
    "ylim": None,
    "xlabel": r"$1-|\Delta\phi_{\mu\mu}|/\pi$",
    "ylabel": r"$d\sigma/dA_\phi$",
    "density_ylabel": r"$\frac{1}{\sigma}\,d\sigma/dA_\phi$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"ATLAS semi-exclusive dimuon acoplanarity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_acoplanarity,
}
