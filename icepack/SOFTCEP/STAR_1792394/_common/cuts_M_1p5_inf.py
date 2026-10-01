# Shared STAR 200 GeV above-1.5-GeV event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .selections import cut_parameters, mass_above

cut_param = cut_parameters()


# Select central systems with invariant mass above 1.5 GeV
def cut_func(event):
    return mass_above(event, 1.5)
