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

import pdg

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
# N.B: icecache functions needs to be called as follows:
# e.g. proj_central_particles(event), but not proj_central_particles(event=event)
#
#
@icecache
def sum_4_momenta(plist):
    particles = [hepmc2vec4(p.momentum()) for p in plist]
    return reduce(lambda a,b: a + b, particles)

@icecache
def proj_central_particles(event):
    return filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() in event.pid \
        and np.abs(p.momentum().eta()) < event.cut_module.CENTRAL_ETA, event.evt.particles())

@icecache
def proj_1D_CS(event):
    
    # Central system sum 4-momentum
    plist     = filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() in event.pid \
        and np.abs(p.momentum().eta()) < event.cut_module.CENTRAL_ETA, event.evt.particles())
    
    particles = [hepmc2vec4(p.momentum()) for p in plist]
    X         = reduce(lambda a,b: a + b, particles)
    
    # Beams
    beam  = [hepmc2vec4(p.momentum()) for p in filter(lambda p: p.status() == pdg.INITIAL_STATE, event.evt.particles())]
    
    # Do the frame transform
    pb1boost,pb2boost,pfboost = LorentFramePrepare(pbeam1=beam[0], pbeam2=beam[1], particles=particles, X=X)
    pfout = LorentzFrame(pb1boost=pb1boost, pb2boost=pb2boost, pfboost=pfboost, frametype='CS')

    return pfout

@icecache
def proj_1D_4body_helicity(event, target_m=0.775):
    """
    Sequential 2->2 body projector, e.g. resonance -> rho > {pi+ pi-} rho > {pi+ pi-}
    
    Returns:
        pions:      4-momentum per system in mother helicity frames
        system:     intermediate system 4-momentum in the lab frame
        particles:  matched particles in the lab frame (two pairs)
    """
    
    # Pions
    pi_pos = filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() == pdg.PDG_PI \
        and np.abs(p.momentum().eta()) < event.cut_module.CENTRAL_ETA, event.evt.particles())
    pi_pos_p4 = [hepmc2vec4(p.momentum()) for p in pi_pos]

    # Anti-pions
    pi_neg = filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() == -pdg.PDG_PI \
        and np.abs(p.momentum().eta()) < event.cut_module.CENTRAL_ETA, event.evt.particles())    
    pi_neg_p4 = [hepmc2vec4(p.momentum()) for p in pi_neg]

    # Beams
    beam = [hepmc2vec4(p.momentum()) for p in filter(lambda p: p.status() == pdg.INITIAL_STATE, event.evt.particles())]
    
    # -----------------------------------------------------
    # Find the closest combination (two combinations)
    Pe   = [ [[0,0],    # Combination 0: system A
              [1,1]],   # Combination 0: system B
             [[0,1],    # Combination 1: system A
              [1,0]] ]  # Combination 1: system B

    loss = np.zeros(len(Pe))
    
    A,B  = 0,1
    for n in range(len(Pe)):

        ind_A, ind_B = Pe[n][A], Pe[n][B]

        mA      = (pi_pos_p4[ind_A[0]] + pi_neg_p4[ind_A[1]]).m
        mB      = (pi_pos_p4[ind_B[0]] + pi_neg_p4[ind_B[1]]).m
        loss[n] = (target_m - mA)**2 + (target_m - mB)**2

    # Pick the better solution
    n = np.argmin(loss)

    ind_A, ind_B = Pe[n][A], Pe[n][B]

    particles = [ [pi_pos_p4[ind_A[0]], pi_neg_p4[ind_A[1]]],
                  [pi_pos_p4[ind_B[0]], pi_neg_p4[ind_B[1]]] ]

    # -----------------------------------------------------
    X = particles[0][0] + particles[0][1] + particles[1][0] + particles[1][1]

    # Over both intermediate systems
    pions_in_HX   = []
    system_in_lab = []
    for i in range(len(particles)):
        
        # Intermediate system Xi in the lab
        Xi_in_lab = particles[i][0] + particles[i][1]
        system_in_lab.append(Xi_in_lab)

        # Intermediate system Xi in the central system resonance X rest frame
        Xi_in_X_frame = copy.deepcopy(Xi_in_lab)
        Xi_in_X_frame.boost(b=X, sign=-1)

        # Transform to the Helicity Frame
        pfout = HXFrame(p=particles[i], X=Xi_in_X_frame)

        # Get the pions
        pions_in_HX.append(pfout)
            
    return pions_in_HX, system_in_lab, particles


def HXFrame(p, X):
    """
    From lab to the "Helicity frame"
    Quantization z-axis as the direction spanned by X in the frame of reference.
    
    Be careful with the frame definitions in multibody cascaded decays, i.e.
    then the intermediate X typically is defined in _its_ own mother frame.
    
    Args:
        p = Set of 4-momentum in the lab to be transformed
            (N.B. sum over [rotated] p is used to define the boost to their rest frame)
        X = Helicity direction 4-momentum
            (N.B. this vector defines only the helicity direction but not the boost here)
    
    This function is identical to HXFrame() in MKinematics.h, which is used
    in the generator (Jacob-Wick 2-body [sequential] amplitudes).
    """

    # ZYZ-sequence rotation angles as defined by the system X direction, note the minus
    Z_angle = -X.phi
    Y_angle = -X.theta

    def rotate(particle):
        particle.rotateZ(Z_angle)
        particle.rotateY(Y_angle)
        particle.rotateZ(np.pi)  # Reflection

    # Rotate particles
    pout = copy.deepcopy(p)
    for i in range(len(pout)):
        rotate(pout[i])

    # Boost direction defined as a sum over the rotated particles
    p_boost = vec4()
    for i in range(len(pout)):
        p_boost += pout[i]

    # Boost each particle
    for i in range(len(pout)):
        pout[i].boost(b=p_boost, sign=-1)

    return pout


## Observables
#
#

def proj_1D_4body_DeltaM_AB(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}
    
    Returns:
        invariant mass difference between intermediate mothers A - B
    """
    _, system_in_lab, _ = proj_1D_4body_helicity(event)
    return system_in_lab[0].m - system_in_lab[1].m


def proj_1D_4body_Deltacos_12(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        angle difference between A1 and B2
    """
    pions_in_HX, _, _ = proj_1D_4body_helicity(event)
    return pions_in_HX[0][0].costheta - pions_in_HX[1][0].costheta


def proj_1D_4body_cos1_dotprod(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}
    
    X-check routine
    
    Gives identical answer to proj_1D_4body_cos1(),
    but using an implicit rotation (i.e. boosts + dot product)
    instead of an explicit frame transform made in HXFrame().
    
    See also: https://en.wikipedia.org/wiki/Wigner_rotation
    
    Returns:
        angle cos(theta_{A1}) [-1,1] in the mother A helicity frame
    """
    _, _, p = proj_1D_4body_helicity(event)
    
    A,B,A1,A2,B1,B2 = 0,1,0,1,0,1
    
    X = p[A][A1] + p[A][A2] + p[B][B1] + p[B][B2]
    
    # Boost final states to the central system X-rest frame
    S_in_X = copy.deepcopy(p[A])
    for i in range(len(S_in_X)):
        S_in_X[i].boost(b=X, sign=-1)
    
    V_in_X = S_in_X[A1] + S_in_X[A2]
    
    # Boost V daughters to the V-rest frame
    V_in_lab = p[A][A1] + p[A][A2] 
    S_in_V = copy.deepcopy(p[A])
    for i in range(len(S_in_V)):
        S_in_V[i].boost(b=V_in_lab, sign=-1)
    
    # Take cos(theta) angle between S in the V-frame and V in the X-frame
    # (V ~ A and S ~ A1)
    costheta = np.cos(S_in_V[A1].angle(V_in_X))
    
    return costheta


def proj_1D_4body_cos1(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        angle cos(theta_{A1}) [-1,1] in the mother A helicity frame
    """
    pions_in_HX, _, _ = proj_1D_4body_helicity(event)
    return pions_in_HX[0][0].costheta


def proj_1D_4body_cos2(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}

    Returns:
        angle cos(theta_{B1}) [-1,1] in the mother B helicity frame
    """
    pions_in_HX, _, _ = proj_1D_4body_helicity(event)
    return pions_in_HX[1][0].costheta


def proj_1D_4body_phi12(event):
    """
    X -> A > {A1 + A2} B > {B1 + B2}
    
    Returns:
        sum angle (phi_{A1} + phi_{B1})
    
    Notes:
        phi_{A1} and phi_{A2} defined in the mother A and B helicity frames, respectively
    """
    pions_in_HX, _, _ = proj_1D_4body_helicity(event)
    angle = pions_in_HX[0][0].phi + pions_in_HX[1][0].phi

    return rad2deg(np.arctan2(np.sin(angle), np.cos(angle))) # Wrap to [-pi,pi]
    #return rad2deg(angle) # raw [-2pi, 2pi] (gives ~ sin/cos modulated 'triangle')


def proj_1D_M(event):
    """
    Central system invariant mass (GeV)
    """
    plist = proj_central_particles(event)
    X     = sum_4_momenta(plist)

    return X.m


def proj_1D_Rap(event):
    """
    Central system rapidity
    """
    plist = proj_central_particles(event)
    X     = sum_4_momenta(plist)
    
    return X.rapidity


def proj_1D_Abs_t1t2(event):
    """
    |t1 + t2| (GeV^2)
    """
    forward  = filter(lambda p: p.pid() == pdg.PDG_PROTON and
        p.momentum().eta() >  event.cut_module.CENTRAL_ETA, event.evt.particles())
    backward = filter(lambda p: p.pid() == pdg.PDG_PROTON and
        p.momentum().eta() < -event.cut_module.CENTRAL_ETA, event.evt.particles())

    # Mandelstam invariants
    t1 = reduce(lambda a,b: (a.momentum() - b.momentum()).m2(), forward)
    t2 = reduce(lambda a,b: (a.momentum() - b.momentum()).m2(), backward)

    return np.abs(t1 + t2)


def proj_1D_dPhi_pp(event):
    """
    Forward proton pair deltaphi in the lab frame (in deg [0,180])
    """
    
    plist = filter(lambda p: p.status() == pdg.FINAL_STATE and p.pid() == pdg.PDG_PROTON \
        and np.abs(p.momentum().eta()) > event.cut_module.CENTRAL_ETA, event.evt.particles())
    beam  = [hepmc2vec4(p.momentum()) for p in plist]
    
    return rad2deg(beam[0].abs_delta_phi(beam[1]))


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


# ------------------------------------------------------------------------


def LorentFramePrepare(pbeam1, pbeam2, particles, X):
    """
    Lorentz Transform preparation function for LorentzFrame()
    """

    pfboost = copy.deepcopy(particles)
    pbeam1b = copy.deepcopy(pbeam1) 
    pbeam2b = copy.deepcopy(pbeam2)

    # Boost each particle to the system X rest frame
    for i in range(len(pfboost)):
        pfboost[i].boost(b=X, sign=-1)

    # Boost initial state protons
    pbeam1b.boost(b=X, sign=-1)
    pbeam2b.boost(b=X, sign=-1)

    return pbeam1b, pbeam2b, pfboost


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
    (LHC proton-proton type). Use HXFrame() for generic helicity frame transforms
    necessary e.g. in cascaded 2-body decays.
    
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

