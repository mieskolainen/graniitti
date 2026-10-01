# STAR coherent rho integrated cross section observable
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

# Compute a single bin for the fiducial cross section
def proj_cross_section(event):
    return 0.5


obs_cross_section = {
    "tag": "cross_section", "xlim": None, "ylim": None, "bins": None, "density": False,
    "xlabel": "Fiducial cross section", "ylabel": r"$\sigma$",
    "units": {"x": "unit", "y": "mb"}, "label": "STAR coherent rho",
    "ylim_ratio": (0.0, 2.0), "ytick_ratio_step": 0.5, "func": proj_cross_section,
}
