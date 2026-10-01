# ALICE forward neutron multiplicity observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis.forward import Acceptance


# Count final neutrons on the selected detector side after particle cuts
def proj_neutron_multiplicity(event):
    side = getattr(event, "cut_param", {}).get("emd_side", 1)
    multiplicity = len(Acceptance(pid=(2112,), side=side).particles(event))
    if multiplicity < 1:
        raise ValueError("Single EMD projection requires a nonzero neutron multiplicity")
    return float(multiplicity)


obs_neutron_multiplicity = {
    "tag": "neutron_multiplicity",
    "xlim": None,
    "ylim": None,
    "xlabel": r"Forward neutron multiplicity $i$",
    "ylabel": r"$1/N_{\mathrm{evt}}\,dN_{\mathrm{evt}}/di$",
    "density_ylabel": r"$1/N_{\mathrm{evt}}\,dN_{\mathrm{evt}}/di$",
    "units": {"x": r"unit", "y": r"unit"},
    "label": r"ALICE single EMD",
    "ylim_ratio": (0.0, 2.5),
    "ytick_ratio_step": 0.5,
    "bins": None,
    "density": False,
    "func": proj_neutron_multiplicity,
}
