# STAR 510 GeV published mass and azimuth region
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .cuts import cut_func as cut_func
from .cuts import parameters

cut_param = parameters(mass=0, phi=None)
