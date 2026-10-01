# Common CMS UPC charged, neutral and ZDC exclusivity selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math

from core.analysis import obs
from core.analysis.forward import Acceptance

cut_param = {
    "charged_pt_min": 0.3,
    "charged_abs_eta_max": 2.4,
    "neutral_et_min": 1.0,
    "neutral_abs_eta_max": 5.2,
    "neutral_deta_max": 0.15,
    "ecal_barrel_abs_eta_max": 1.479,
    "electron_dphi_max": (0.7, 0.4),
    "photon_dphi_max": 0.15,
    "zdc_abs_eta_min": 8.3,
    "zdc_energy_max": 7000.0,
}

CHARGED_ABS_PDGS = frozenset((11, 13, 15, 211, 321, 2212, 3112, 3222, 3312, 3334))


# Exclude ECAL photons in the candidate windows while retaining all HCAL activity
# [REFERENCE: arXiv:2412.15413, Section 5]
def in_ecal_window(record, selected, param):
    if record["pid"] != 22:
        return False
    for candidate in selected:
        endcap = abs(candidate["eta"]) >= param["ecal_barrel_abs_eta_max"]
        dphi = param["electron_dphi_max"][endcap] if abs(candidate["pid"]) == 11 else param["photon_dphi_max"]
        if (abs(record["eta"] - candidate["eta"]) < param["neutral_deta_max"]
                and record["p4"].abs_delta_phi(candidate["p4"]) < dphi):
            return True
    return False


# Apply final-particle exclusivity and require a quiet ZDC on at least one side
# [REFERENCE: arXiv:2412.15413, Section 5, Tables 1, 3 and 6]
def exclusivity_cut(event, selected):
    param = event.cut_param
    ids = {record["id"] for record in selected}
    for record in obs.proj_event_particles(event):
        if not record["is_final"] or record["id"] in ids:
            continue
        pid, eta = abs(record["pid"]), abs(record["eta"])
        charged = pid in CHARGED_ABS_PDGS or (pid >= 1000000000 and (pid // 10000) % 1000 > 0)
        if charged:
            if record["pt"] > param["charged_pt_min"] and eta < param["charged_abs_eta_max"]:
                return False
        elif pid not in (12, 14, 16) and eta < param["neutral_abs_eta_max"]:
            et = record["p4"].e / math.cosh(eta)
            if et > param["neutral_et_min"] and not in_ecal_window(record, selected, param):
                return False
    return any(
        Acceptance(pid=(2112, -2112), side=side, abs_eta=(param["zdc_abs_eta_min"], math.inf)).energy_sum(event)
        < param["zdc_energy_max"] for side in (-1, 1)
    )
