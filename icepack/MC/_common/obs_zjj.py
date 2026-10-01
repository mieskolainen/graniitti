# Phase space closure observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .observables import flavour_pair, production

globals().update(production(120.0, 300.0))
globals().update(flavour_pair(((-13, 13), (2, -2)), (70.0, 110.0), (0.0, 220.0)))
