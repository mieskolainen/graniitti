# Inclusive diffraction validation selection
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

cut_param = {}


# Accept every event generated inside the steering-card fiducial phase space
def cut_func(event):
    return True
