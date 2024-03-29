# Default minimal event selection applied on HepMC3 files
# 
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np

import pdg
from icevec import *

# ------------------------------------
# Define cut parameters here (used in projectors in iceobs)
CENTRAL_ETA = 6
CENTRAL_PT  = 0

# ------------------------------------

def cut_func(event):
    
    plist = filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() in event.pid \
        and np.abs(p.momentum().eta()) < CENTRAL_ETA, event.evt.particles())
    
    p = [hepmc2vec4(p.momentum()) for p in plist]
    
    if len(p) == len(event.pid):
        return True
    else:
        return False
