# Common final-state acceptance for the hard diffraction comparison
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.
from pathlib import Path

import json5
from core.analysis import obs

from ..z_with_pythia import cuts as base

with (Path(__file__).resolve().parents[3] / "modeldata/TUNE0/GENERAL.json").open() as source:
    model = json5.load(source)["PARAM_HARDPOMERON"]
cut_param = dict(base.cut_param)


# Require the same final-state proton xi and t support for both generators
def cut_func(event):
    if not base.cut_func(event):
        return False
    records = obs.proj_event_particles(event)
    beams = [r["p4"] for r in records if r["status"] == 4 and r["pid"] == 2212]
    protons = [r["p4"] for r in records if r["is_final"] and r["pid"] == 2212]
    for beam in beams:
        side = 1 if beam.z > 0 else -1
        candidates = [p for p in protons if side * p.z > 0]
        if not candidates:
            continue
        p = max(candidates, key=lambda v: v.t + side * v.z)
        xi = 1 - (p.t + side * p.z) / (beam.t + side * beam.z)
        t = (beam - p).m2
        if model["xi_range"][0] <= xi <= model["xi_range"][1] and model["t_range"][0] <= t <= model["t_range"][1]:
            return True
    return False
