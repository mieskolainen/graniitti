# MC observables constructions
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from functools import lru_cache
from cachetools import cached
import copy
import numpy as np
import pyhepmc_ng as hepmc3
import numba
import re
import time

from functools import reduce
from tqdm import tqdm

from icecache import icecache
from icevec import *


# PDG codes (p,pi+,K+)
PDG_PROTON  = 2212
PDG_PI      = 211
PDG_K       = 321


# PDG (HepMC) status codes
INITIAL_STATE = 4
FINAL_STATE   = 1



# ========================================================================
# Helper functions

@numba.njit
def unitvec(v):
    """ Return unit vector """
    return v / np.linalg.norm(v)

@numba.njit
def rad2deg(x):
    return x / (2*np.pi) * 360

@numba.njit
def deg2rad(x):
    return x / 360 * (2*np.pi)


# ========================================================================

## Common projectors
#
#
@icecache
def proj_central_particles(event, beam_eta=6):
    return filter(lambda p: p.status == FINAL_STATE and p.pid in event.pid \
        and np.abs(p.momentum.eta()) < beam_eta, event.evt.particles)

@icecache
def sum_4_momenta(particles):
    return reduce(lambda a,b: a.momentum + b.momentum, particles)

@icecache
def proj_1D_CS(event, beam_eta=6):
    
    # Central system sum 4-momentum
    plist     = filter(lambda p: p.status == FINAL_STATE and p.pid in event.pid \
        and np.abs(p.momentum.eta()) < beam_eta, event.evt.particles)
    particles = [hepmc2vec4(p.momentum) for p in plist]
    X         = reduce(lambda a,b: a + b, particles)
    
    # Beams
    beam  = [hepmc2vec4(p.momentum) for p in filter(lambda p: p.status == INITIAL_STATE, event.evt.particles)]
    
    # Do the frame transform
    pb1boost,pb2boost,pfboost = LorentFramePrepare(pbeam1=beam[0], pbeam2=beam[1], particles=particles, X=X)
    pfout = LorentzFrame(pb1boost=pb1boost, pb2boost=pb2boost, pfboost=pfboost, frametype='CS')

    return pfout


## Observables
#
#
def proj_1D_M(event):
    """
    Central system invariant mass (GeV)
    """
    plist = proj_central_particles(event)
    X     = sum_4_momenta(plist)

    return X.m()


def proj_1D_Rap(event):
    """
    Central system rapidity
    """
    plist = proj_central_particles(event)
    X     = sum_4_momenta(plist)
    
    return X.rap()


def proj_1D_Abs_t1t2(event, beam_eta=6):
    """
    |t1 + t2| (GeV^2)
    """
    forward  = filter(lambda p: p.pid == PDG_PROTON and p.momentum.eta() >  beam_eta, event.evt.particles)
    backward = filter(lambda p: p.pid == PDG_PROTON and p.momentum.eta() < -beam_eta, event.evt.particles)

    # Mandelstam invariants
    t1 = reduce(lambda a,b: (a.momentum - b.momentum).m2(), forward)
    t2 = reduce(lambda a,b: (a.momentum - b.momentum).m2(), backward)

    return np.abs(t1 + t2)


def proj_1D_costheta_CS(event):
    """
    Cos(theta) angle of pi+(K+) in the Collins-Soper rest frame [-1,1]
    """
    pfout = proj_1D_CS(event)
    
    return pfout[0].costheta


def proj_1D_phi_CS(event, in_deg=True):
    """
    Cos(theta) angle of pi+(K+, ...) in the Collins-Soper rest frame [-1,1]
    """
    pfout = proj_1D_CS(event)
    
    # Phi angle in the new rest frame
    phi = pfout[0].phi

    return rad2deg(phi) if in_deg else phi


def proj_1D_dPhi_pp(event, beam_eta=6):
    """
    Forward proton pair deltaphi in the lab frame (in deg [0,180])
    """
    
    plist = filter(lambda p: p.status == FINAL_STATE and p.pid == PDG_PROTON \
        and np.abs(p.momentum.eta()) > beam_eta, event.evt.particles)
    beam  = [hepmc2vec4(p.momentum) for p in plist]
    
    return rad2deg(beam[0].abs_delta_phi(beam[1]))


# ------------------------------------------------------------------------


def LorentFramePrepare(pbeam1, pbeam2, particles, X):
    """
    Lorentz Transform prepatation function for LorentzFrame()
    """

    # Boost each particle to the system rest frame
    for i in range(len(particles)):
        particles[i].boost(b=X, sign=-1)

    # Boost initial state protons
    pbeam1.boost(b=X, sign=-1)
    pbeam2.boost(b=X, sign=-1)

    return pbeam1, pbeam2, particles


def LorentzFrame(pb1boost, pb2boost, pfboost, frametype='CS', direction=1):
    """
    Rotated Lorentz frame transformation
    
    "CS" : Collins-Soper
    "AH" : Anti-Helicity (Anti-CS)
    "HE" : Helicity
    "PG" : Pseudo-Gottfried-Jackson
    
    Args:
        pb1boost,pb2boost : Boosted beam 4-momenta
        pfboost           : Boosted 4-momenta
        frametype         : Lorentz frame type
    
    Returns:
        rotated 4-momenta
    """

    pb1boost3 = pb1boost.p3
    pb2boost3 = pb2boost.p3

    # Frame rotation x-y-z-axes

    # @@ NON-ROTATED FRAME AXIS DEFINITION @@
    if   (frametype == "CM"):
        zaxis = np.array([0, 0, 1])
        yaxis = np.array([0, 1, 0])

    # @@ COLLINS-SOPER FRAME POLARIZATION AXIS DEFINITION @@
    elif (frametype == "CS" or frametype == "G1"):
        zaxis = unitvec(unitvec(pb1boost3) - unitvec(pb2boost3))

    # @@ ANTI-HELICITY FRAME POLARIZATION AXIS DEFINITION @@
    elif (frametype == "AH" or frametype == "G2"):
        zaxis = unitvec(unitvec(pb1boost3) + unitvec(pb2boost3))

    # @@ HELICITY FRAME POLARIZATION AXIS DEFINITION @@
    elif (frametype == "HX" or frametype == "G3"):
        zaxis = unitvec(-(pb1boost3 + pb2boost3))
    
    # @@ PSEUDO-GOTTFRIED-JACKSON AXIS DEFINITION: [1] or [2] @@
    elif (frametype == "PG" or frametype == "G4"):
        if   (direction == -1):
            zaxis = unitvec(pb1boost3)
        elif (direction ==  1):
            zaxis = unitvec(pb2boost3)
        else:
            raise Exception(f'LorentzFrame: Invalid direction {direction}')
    else:
        raise Exception(f'LorentzFrame: Unknown frame <{frametype}>')

    # y-axis
    if (frametype != "CM"):
        yaxis = unitvec(np.cross(unitvec(pb1boost3), unitvec(pb2boost3)))

    # x-axis
    xaxis = unitvec(np.cross(yaxis, zaxis))  # x = y [cross product] z

    # Create SO(3) rotation matrix for the new coordinate axes
    R = np.array([xaxis, yaxis, zaxis])      # Axes as rows

    # Rotate all vectors
    for k in range(len(pfboost)):
        pfboost[k].rotateSO3(R)

    return pfboost

