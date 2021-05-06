# Default minimal event selection applied on HepMC3 files
# 
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np

from iceobs import *


def cut_func(event, central_eta=6):
    plist = filter(lambda p: p.status == FINAL_STATE and p.pid in event.pid \
        and np.abs(p.momentum.eta()) < central_eta, event.evt.particles)
    
    p = [hepmc2vec4(p.momentum) for p in plist]
    
    if len(p) == len(event.pid):

    	# Sum the system momentum
    	#X = reduce(lambda a,b: a + b, p)
    	
        return True
    else:
        return False
