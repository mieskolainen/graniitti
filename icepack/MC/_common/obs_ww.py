# Phase space closure observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .observables import flavour_pair, production

globals().update(production(150.0, 300.0))
globals().update(flavour_pair(((-13, 14), (11, -12)), (60.0, 100.0), (60.0, 100.0)))
