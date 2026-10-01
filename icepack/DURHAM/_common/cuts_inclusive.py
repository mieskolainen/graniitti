# Inclusive event selection for Durham QCD icepack comparisons
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

cut_param = {}


# Accept events inside the generator level fiducial phase space
def cut_func(_event):
    return True
