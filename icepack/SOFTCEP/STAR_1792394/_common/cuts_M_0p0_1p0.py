# Shared STAR 200 GeV 0 to 1 GeV event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .selections import cut_parameters, mass_at_most

cut_param = cut_parameters()


# Select central systems with invariant mass at or below 1 GeV
def cut_func(event):
    return mass_at_most(event, 1.0)
