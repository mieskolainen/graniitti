# CMS 7 TeV dipion selection at stable particle level
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {"PION_PT_MIN": 0.2, "PION_ABS_RAP_MAX": 2.0}


# Require exactly two opposite charge pions in the published fiducial region
def cut_func(event):
    # [REFERENCE: CMS arXiv:1706.08310, Sec. 4.3]
    pions = [
        record for record in obs.proj_event_particles(event)
        if record["is_final"] and abs(record["pid"]) == 211
        and record["pt"] > event.cut_param["PION_PT_MIN"]
        and abs(record["p4"].rapidity) < event.cut_param["PION_ABS_RAP_MAX"]
    ]
    return len(pions) == 2 and pions[0]["pid"] + pions[1]["pid"] == 0
