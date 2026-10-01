# TOTEM 8 TeV elastic CNI histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.obs import obs_mandelstam_t as _obs_mandelstam_t

obs_mandelstam_t = {**_obs_mandelstam_t, "xscale": "log", "xlim": (0.0006, 0.2000)}
