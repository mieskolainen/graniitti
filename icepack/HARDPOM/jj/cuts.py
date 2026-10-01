# Select the two central partons in hard diffractive dijet events
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis import obs

cut_param = {}


# Require two selected central partons within the generator fiducial phase space
def cut_func(event):
    return len(obs.proj_central_particles(event)) == 2
