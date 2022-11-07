# None event selection applied on HepMC3 files
# 
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np

import pdg
from icevec import *

# ------------------------------------
# Define cut parameters here (used in projectors in iceobs)
CENTRAL_ETA = 999
CENTRAL_PT  = 0

# ------------------------------------

def cut_func(event):
    return True
