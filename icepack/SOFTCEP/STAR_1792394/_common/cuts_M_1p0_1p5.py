# Shared STAR 200 GeV 1 to 1.5 GeV event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .selections import cut_parameters, mass_interval

cut_param = cut_parameters()


# Select central systems with invariant mass above 1 and at or below 1.5 GeV
def cut_func(event):
    return mass_interval(event, 1.0, 1.5)
