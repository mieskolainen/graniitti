# Phase space closure observables
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .observables import production

globals().update(production(500.0, 600.0))
