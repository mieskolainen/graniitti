# Shared STAR 200 GeV above-90-degree event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .selections import azimuth_above, cut_parameters

cut_param = cut_parameters()


# Select events with outgoing-proton azimuthal separation above 90 degrees
def cut_func(event):
    return azimuth_above(event, 90.0)
