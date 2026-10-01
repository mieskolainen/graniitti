# STAR 200 GeV proton-antiproton histogram definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.observables import (
    absolute_momentum_transfer_sum,
    collins_soper_azimuthal_angle,
    collins_soper_polar_angle,
    invariant_mass,
    proton_azimuthal_separation,
    rapidity,
)

obs_M = invariant_mass()
obs_dPhi_pp = proton_azimuthal_separation()
obs_Abs_t1t2 = absolute_momentum_transfer_sum()
obs_Rap = rapidity()
obs_costheta_CS = collins_soper_polar_angle()
obs_phi_CS = collins_soper_azimuthal_angle()
