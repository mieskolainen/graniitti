# Common observables for soft central exclusive process comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import copy
import importlib

import numpy as np
from core.analysis import obs, pdg
from core.io.cache import cache

FOUR_BODY = importlib.import_module("core.analysis.observables.four_body")


# Compute the stable pair attached to the direct central production vertex
@cache
def central_pair(event):
    if len(event.pid) != 2 or event.pid[0] != -event.pid[1]:
        raise ValueError(
            f"Two body continuum observables require a particle antiparticle pair, found {event.pid}"
        )

    records = obs.proj_central_particle_records(event)
    pair = []
    for pid_value in event.pid:
        candidates = [record["p4"] for record in records if record["pid"] == pid_value]
        if len(candidates) != 1:
            raise ValueError(
                "Two body continuum observables require exactly one central particle "
                f"with PDG {pid_value}, found {len(candidates)}"
            )
        pair.append(candidates[0])
    return pair


# Transform the direct central pair to the Collins Soper frame
@cache
def collins_soper_pair(event):
    pair = central_pair(event)
    central_system = sum(pair)
    beams = [
        record["p4"]
        for record in obs.proj_event_particles(event)
        if record["status"] == pdg.INITIAL_STATE
    ]
    if len(beams) != 2:
        raise ValueError(
            f"Two body continuum observables require two beam particles, found {len(beams)}"
        )

    beam_1, beam_2, central_pair_in_frame = obs.LorentFramePrepare(
        pbeam1=beams[0], pbeam2=beams[1], particles=pair, X=central_system
    )
    return obs.LorentzFrame(
        pb1boost=beam_1,
        pb2boost=beam_2,
        pfboost=central_pair_in_frame,
        frametype="CS",
    )


# Compute the positive particle Collins Soper polar angle
def collins_soper_costheta(event):
    return collins_soper_pair(event)[0].costheta


# Compute the positive particle Collins Soper azimuth in the full angular range
def collins_soper_phi(event):
    return np.mod(obs.rad2deg(collins_soper_pair(event)[0].phi), 360.0)


obs_M = {
    "tag": "M",
    "xlim": (0.0, 6.0),
    "ylim": None,
    "xlabel": r"$M_X$",
    "ylabel": r"$d\sigma/dM_X$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Central system invariant mass",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 6.0, 61),
    "density": False,
    "func": obs.proj_1D_M,
}

obs_Rap = {
    "tag": "Rap",
    "xlim": (-3.0, 3.0),
    "ylim": None,
    "xlabel": r"$Y_X$",
    "ylabel": r"$d\sigma/dY_X$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Central system rapidity",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-3.0, 3.0, 49),
    "density": False,
    "func": obs.proj_1D_Rap,
}

obs_Pt = {
    "tag": "Pt",
    "xlim": (0.0, 2.0),
    "ylim": None,
    "xlabel": r"$p_{T,X}$",
    "ylabel": r"$d\sigma/dp_{T,X}$",
    "units": {"x": r"GeV", "y": r"pb"},
    "label": r"Central system transverse momentum",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 2.0, 41),
    "density": False,
    "func": obs.proj_1D_Pt,
}

obs_t1 = {
    "tag": "t1",
    "xlim": (0.0, 2.0),
    "ylim": None,
    "xlabel": r"$|t_1|$",
    "ylabel": r"$d\sigma/d|t_1|$",
    "units": {"x": r"GeV$^{2}$", "y": r"pb"},
    "label": r"Positive side proton momentum transfer",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 2.0, 41),
    "density": False,
    "func": obs.proj_1D_t1,
}

obs_t2 = {
    "tag": "t2",
    "xlim": (0.0, 2.0),
    "ylim": None,
    "xlabel": r"$|t_2|$",
    "ylabel": r"$d\sigma/d|t_2|$",
    "units": {"x": r"GeV$^{2}$", "y": r"pb"},
    "label": r"Negative side proton momentum transfer",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 2.0, 41),
    "density": False,
    "func": obs.proj_1D_t2,
}

obs_Abs_t1t2 = {
    "tag": "Abs_t1t2",
    "xlim": (0.0, 2.0),
    "ylim": None,
    "xlabel": r"$|t_1+t_2|$",
    "ylabel": r"$d\sigma/d|t_1+t_2|$",
    "units": {"x": r"GeV$^{2}$", "y": r"pb"},
    "label": r"Absolute sum of proton momentum transfers",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 2.0, 41),
    "density": False,
    "func": obs.proj_1D_Abs_t1t2,
}

obs_costheta_CS = {
    "tag": "costheta_CS",
    "xlim": (-1.0, 1.0),
    "ylim": None,
    "xlabel": r"$\cos\theta_{\mathrm{CS}}$",
    "ylabel": r"$d\sigma/d\cos\theta_{\mathrm{CS}}$",
    "units": {"x": r"unit", "y": r"pb"},
    "label": r"Collins Soper polar angle",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(-1.0, 1.0, 41),
    "density": False,
    "func": collins_soper_costheta,
}

obs_phi_CS = {
    "tag": "phi_CS",
    "xlim": (0.0, 360.0),
    "ylim": None,
    "xlabel": r"$\phi_{\mathrm{CS}}$",
    "ylabel": r"$d\sigma/d\phi_{\mathrm{CS}}$",
    "units": {"x": r"deg", "y": r"pb"},
    "label": r"Collins Soper azimuthal angle",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 360.0, 49),
    "density": False,
    "func": collins_soper_phi,
}

obs_dPhi_pp = {
    "tag": "dPhi_pp",
    "xlim": (0.0, 180.0),
    "ylim": None,
    "xlabel": r"$\Delta\phi_{pp}$",
    "ylabel": r"$d\sigma/d\Delta\phi_{pp}$",
    "units": {"x": r"deg", "y": r"pb"},
    "label": r"Outgoing proton azimuthal separation",
    "ylim_ratio": (0.0, 2.0),
    "ytick_ratio_step": 0.5,
    "bins": np.linspace(0.0, 180.0, 37),
    "density": False,
    "func": obs.proj_1D_dPhi_pp,
}

# Reuse the standard four body definitions for topology specific comparisons
obs_4body_M_A = copy.deepcopy(FOUR_BODY.obs_4body_M_A)
obs_4body_M_B = copy.deepcopy(FOUR_BODY.obs_4body_M_B)
obs_4body_cos1 = copy.deepcopy(FOUR_BODY.obs_4body_cos1)
obs_4body_cos2 = copy.deepcopy(FOUR_BODY.obs_4body_cos2)
obs_4body_phi12 = copy.deepcopy(FOUR_BODY.obs_4body_phi12)
