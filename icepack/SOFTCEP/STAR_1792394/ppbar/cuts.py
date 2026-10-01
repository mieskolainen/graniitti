# STAR 200 GeV proton-antiproton event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from .._common.selections import accept_all, cut_parameters

cut_param = cut_parameters()


# Accept events passing the generator-level STAR fiducial selection
def cut_func(event):
    return accept_all(event)
