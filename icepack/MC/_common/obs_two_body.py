# Phase space closure observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .observables import observable, production

globals().update(production(0.0, 4.0))
obs_costheta_CS = observable("costheta_CS", -1.0, 1.0)
obs_Rap = observable("Rap", -0.9, 0.9)
