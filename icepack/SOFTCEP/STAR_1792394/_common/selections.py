# Shared STAR 200 GeV central-production event selections
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs


# Compute independent STAR selection parameters
def cut_parameters():
    return {}


# Accept one event without an additional analysis-level selection
def accept_all(event):
    return True


# Select a central system at or below one invariant-mass boundary
def mass_at_most(event, upper):
    return obs.proj_1D_M(event) <= upper


# Select a central system inside one open-closed invariant-mass interval
def mass_interval(event, lower, upper):
    mass = obs.proj_1D_M(event)
    return lower < mass <= upper


# Select a central system above one invariant-mass boundary
def mass_above(event, lower):
    return obs.proj_1D_M(event) > lower


# Select events below one outgoing-proton azimuthal-separation boundary
def azimuth_below(event, upper):
    return obs.proj_1D_dPhi_pp(event) < upper


# Select events above one outgoing-proton azimuthal-separation boundary
def azimuth_above(event, lower):
    return obs.proj_1D_dPhi_pp(event) > lower


# Select the proton azimuth and transverse vector difference [REFERENCE: arXiv:2004.11078, Figure 13]
def dpt_cut(event):
    if not azimuth_below(event, event.cut_param["dPhi"]):
        return False
    _, protons = obs.proj_init_final_protons(event)
    dpt = (protons[0] - protons[1]).pt
    return dpt > event.cut_param["dPt"] if event.cut_param["above"] else dpt < event.cut_param["dPt"]
