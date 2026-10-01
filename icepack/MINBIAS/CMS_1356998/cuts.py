# CMS 7 TeV diffraction event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

cut_param = {}


# Accept events for projection into the analysis-specific diffraction regions
def cut_func(event):
    return True
