# Shared elastic event selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

cut_param = {}


# Accept events within the generator level elastic fiducial phase space
def cut_func(event):
    return True
