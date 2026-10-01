# Neutron selection from final HepMC3 beam remnants
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from core.analysis.forward import neutron_class

from . import cuts as common

cut_param = common.cut_param


# Apply the fiducial cuts and the measured neutron class
def cut_func(event):
    return common.cut_func(event) and neutron_class(event) in ('Xn0n',)
