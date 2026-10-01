# HERA photoproduction comparison definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from icepack.PHOTOPROD._common.angles import obs_costheta_HX, obs_phi_HX

obs_costheta_HX = {**obs_costheta_HX, 'density': False, 'ylabel': r'$d\sigma/d\cos\theta_h$', 'units': {'x':'', 'y':'pb'}}
obs_phi_HX = {**obs_phi_HX, 'density': False, 'ylabel': r'$d\sigma/d\phi_h$', 'units': {'x':'rad', 'y':'pb'}}
