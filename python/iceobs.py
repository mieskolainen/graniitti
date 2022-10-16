# MC observables constructions
#
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from functools import lru_cache
from cachetools import cached
import copy
import numpy as np
import numba
import re
import time
from pyHepMC3 import HepMC3 as hepmc3

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
DECAYED_STATE = 2
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
def proj_central_particles(event, central_eta=6):
    return filter(lambda p: p.status() == FINAL_STATE and p.pid() in event.pid \
        and np.abs(p.momentum().eta()) < central_eta, event.evt.particles())

@icecache
def sum_4_momenta(particles):
    return reduce(lambda a,b: a.momentum() + b.momentum(), particles)

@icecache
def proj_1D_CS(event, central_eta=6):
    
    # Central system sum 4-momentum
    plist     = filter(lambda p: p.status() == FINAL_STATE and p.pid() in event.pid \
        and np.abs(p.momentum().eta()) < central_eta, event.evt.particles())
    particles = [hepmc2vec4(p.momentum()) for p in plist]
    X         = reduce(lambda a,b: a + b, particles)
    
    # Beams
    beam  = [hepmc2vec4(p.momentum()) for p in filter(lambda p: p.status() == INITIAL_STATE, event.evt.particles())]
    
    # Do the frame transform
    pb1boost,pb2boost,pfboost = LorentFramePrepare(pbeam1=beam[0], pbeam2=beam[1], particles=particles, X=X)
    pfout = LorentzFrame(pb1boost=pb1boost, pb2boost=pb2boost, pfboost=pfboost, frametype='CS')

    return pfout

#@icecache (icecache does not work for this)
def proj_1D_4body_helicity(event, central_eta=6, target_m=0.775):
    """
    Sequential 2->2 body projector, e.g. resonance -> rho > {pi+ pi-} rho > {pi+ pi-}

    Returns:
        pi+ 4-momentum in the helicity frame of the intermediate mother, respectively.
    """
    
    # Pions
    pi_pos = filter(lambda p: p.status() == FINAL_STATE and p.pid() == PDG_PI \
        and np.abs(p.momentum().eta()) < central_eta, event.evt.particles())
    pi_pos_p4 = [hepmc2vec4(p.momentum()) for p in pi_pos]
    
    # Anti-pions
    pi_neg = filter(lambda p: p.status() == FINAL_STATE and p.pid() == -PDG_PI \
        and np.abs(p.momentum().eta()) < central_eta, event.evt.particles())    
    pi_neg_p4 = [hepmc2vec4(p.momentum()) for p in pi_neg]

    # Find the closest combination (two combinations)
    Pe   = [ [[0,0],    # A- C1
              [1,1]],   # B- C1
             [[0,1],    # A- C2
              [1,0]] ]  # B- C2
    loss = np.zeros(len(Pe))
    
    A,B  = 0,1
    for n in range(len(Pe)):

        i,j,k,l = Pe[n][A][0],Pe[n][A][1], Pe[n][B][0],Pe[n][B][1]

        mA      = (pi_pos_p4[i] + pi_neg_p4[j]).m
        mB      = (pi_pos_p4[k] + pi_neg_p4[l]).m
        loss[n] = (target_m - mA)**2 + (target_m - mB)**2

    # Pick the solution
    n = np.argmin(loss)

    i,j,k,l   = Pe[n][A][0],Pe[n][A][1], Pe[n][B][0],Pe[n][B][1]
    particles = [ [pi_pos_p4[i], pi_neg_p4[j]],
                  [pi_pos_p4[k], pi_neg_p4[l]] ]

    # Beams
    beam = [hepmc2vec4(p.momentum()) for p in filter(lambda p: p.status() == INITIAL_STATE, event.evt.particles())]
    
    # Over both intermediate systems
    pions = []
    for i in range(2):
        
        # Intermediate system X_i
        X_i = particles[i][0] + particles[i][1]
        
        # Boost & rotate to the intermediate mother system X_i helicity frame
        pb1boost,pb2boost,pfboost = LorentFramePrepare(pbeam1=beam[0], pbeam2=beam[1], particles=particles[i], X=X_i)
        pfout = LorentzFrame(pb1boost=pb1boost, pb2boost=pb2boost, pfboost=pfboost, frametype="HX")

        # Get the first (positive) daughter
        pions.append(pfout[0])

    return pions


## Observables
#
#

def proj_1D_4body_cos1(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        angle cos(theta_{A1}) [-1,1] in the mother A helicity frame
    """
    p = proj_1D_4body_helicity(event=event)
    return p[0].costheta

def proj_1D_4body_cos2(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        angle cos(theta_{B1}) [-1,1] in the mother B helicity frame
    """
    p = proj_1D_4body_helicity(event=event)
    return p[1].costheta

def proj_1D_4body_phi12(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        sum angle (phi_{A1} + phi_{B1})
    
    Notes:
        phi_{A1} and phi_{A2} defined in the mother A and B helicity frames, respectively
    """
    p = proj_1D_4body_helicity(event=event)
    angle = p[0].phi + p[1].phi

    return rad2deg(np.arctan2(np.sin(angle), np.cos(angle))) # Wrap to [-pi,pi]
    #return rad2deg(angle) # raw [-2pi, 2pi] (gives ~ sin/cos modulated 'triangle')


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


def proj_1D_Abs_t1t2(event, central_eta=6):
    """
    |t1 + t2| (GeV^2)
    """
    forward  = filter(lambda p: p.pid() == PDG_PROTON and p.momentum().eta() >  central_eta, event.evt.particles())
    backward = filter(lambda p: p.pid() == PDG_PROTON and p.momentum().eta() < -central_eta, event.evt.particles())

    # Mandelstam invariants
    t1 = reduce(lambda a,b: (a.momentum() - b.momentum()).m2(), forward)
    t2 = reduce(lambda a,b: (a.momentum() - b.momentum()).m2(), backward)

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


def proj_1D_dPhi_pp(event, central_eta=6):
    """
    Forward proton pair deltaphi in the lab frame (in deg [0,180])
    """
    
    plist = filter(lambda p: p.status() == FINAL_STATE and p.pid() == PDG_PROTON \
        and np.abs(p.momentum().eta()) > central_eta, event.evt.particles())
    beam  = [hepmc2vec4(p.momentum()) for p in plist]
    
    return rad2deg(beam[0].abs_delta_phi(beam[1]))


# ------------------------------------------------------------------------


def LorentFramePrepare(pbeam1, pbeam2, particles, X):
    """
    Lorentz Transform preparation function for LorentzFrame()
    """

    # Boost each particle to the system X rest frame
    for i in range(len(particles)):
        particles[i].boost(b=X, sign=-1)

    # Boost initial state protons
    pbeam1.boost(b=X, sign=-1)
    pbeam2.boost(b=X, sign=-1)

    return pbeam1, pbeam2, particles


def LorentzFrame(pb1boost, pb2boost, pfboost, frametype='CS', direction=1):
    """
    Rotated Lorentz frame transformation
    
    "CM" : No-rotation, just boost
    "CS" : Collins-Soper
    "AH" : Anti-Helicity (Anti-CS)
    "HE" : Helicity
    "PG" : Pseudo-Gottfried-Jackson
    
    For more description, see MKinematics.h (LorentzFramePrepare function)
    
    N.B. For the helicity frame, this function is compatible with symmetric beam energies
    (LHC proton-proton type)
    
    Args:
        pb1boost,pb2boost : Boosted beam 4-momenta
        pfboost           : Boosted 4-momenta
        frametype         : Lorentz frame type
        direction         : Forward/backward direction choice for PG-frame

    Returns:
        rotated 4-momenta
    """

    pb1boost3 = pb1boost.p3
    pb2boost3 = pb2boost.p3

    # Frame rotation x-y-z-axes

    # @@ NON-ROTATED FRAME AXIS DEFINITION @@
    if   (frametype == "CM"):
        zaxis = np.array([0.0, 0.0, 1.0],dtype=float)
        yaxis = np.array([0.0, 1.0, 0.0],dtype=float)
    
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

