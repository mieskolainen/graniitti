# MC observables constructed fast
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
from functools import reduce
from tqdm import tqdm
import time

from importlib import import_module


import iceplot
from icecache import icecache
from icevec import *


# PDG codes
PDG_PROTON  = 2212
PDG_PI      = 211
PDG_K       = 321


# PDG (HepMC) status codes
INITIAL_STATE = 4
FINAL_STATE   = 1


@numba.njit
def rad2deg(x):
    return x / (2*np.pi) * 360

@numba.njit
def deg2rad(x):
    return x / 360 * (2*np.pi)


# Object wrapper (for icecache)
class ice_event:
    def __init__(self, evt, pid):
        self.evt = evt
        self.pid = pid


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


@numba.njit
def unitvec(v):
    """ Return unit vector """
    return v / np.linalg.norm(v)


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


def read_hepmc3(hepmc3file, all_obs, pid, cuts=None, maxevents=None, chunk_range=[0, int(1E12)]):
    """
    Read in a single HepMC3 file and compute observables
    """
    mc_data = {}
    for OBS in all_obs.keys():
        mc_data[OBS] = []

    mc_weights       = []
    mc_xsec_pb       = 0
    mc_xsec_pb_err   = 0


    mc_weights_tot   = []

    # --------------------------------------------------------------------
    def identity_cut(event):
        return True

    if cuts is not None:
        MySelection = import_module(cuts)
        cut_func = MySelection.cut_func
    else:
        cut_func = identity_cut
    # --------------------------------------------------------------------

    # Open HepMC3 file
    t = time.time()

    with hepmc3.open(hepmc3file, "r") as f:

        # Loop over events
        for k, evt in enumerate(f):

            event = ice_event(evt=evt, pid=pid)

            # This wont work with hepmc3-ng (python bindings) currently
            #xsection = evt.cross_section
            
            # ------------------------------------------------------------
            # Read it from GRANIITTI header            
            if k == 0:
                string = repr(evt.run_info)
                tmp = re.search("{(.*?)}", string).group(1).split(",")
                mc_xsec_pb     = float(tmp[0].split(":")[1])
                mc_xsec_err_pb = float(tmp[1].split(":")[1])
            # ------------------------------------------------------------

            # -------------------------------
            # This after the header reader
            if k < chunk_range[0]:
                continue
            if k > chunk_range[1]:
                break
            # -------------------------------

            # Event weight for each event
            mc_weights_tot.append(evt.weights[0])

            # -------------------------------
            # Apply additional fiducial cuts here ...
            if not cut_func(event):
                continue;
            # -------------------------------

            # Construct all observables
            for OBS in all_obs.keys():
                try:
                    mc_data[OBS].append( all_obs[OBS]['func'](event=event))
                except:
                    continue;

            # Event weight for passed events
            mc_weights.append(evt.weights[0])

    mc_weights_tot = np.array(mc_weights_tot)
    mc_weights     = np.array(mc_weights)

    # Count the user cuts fiducial acceptance
    acceptance     = np.sum(mc_weights) / np.sum(mc_weights_tot)

    elapsed        = time.time() - t
    print(f'{hepmc3file}')
    print(f'Before cuts: mc_xs_pb: {mc_xsec_pb:.3g} +- {mc_xsec_err_pb:.3g} | events read = {len(mc_weights_tot)} | wsum = {np.sum(mc_weights_tot):0.3E} | chunk = {chunk_range}')
    
    if acceptance < 1.0:
        mc_xsec_pb     *= acceptance
        mc_xsec_err_pb *= acceptance
    
    print(f'After  cuts: mc_xs_pb: {mc_xsec_pb:.3g} +- {mc_xsec_err_pb:.3g} | events read = {len(mc_weights)} | wsum = {np.sum(mc_weights):0.3E} | cutset <{cuts}> acceptance = {acceptance:0.3E}')
    
    return {'data': mc_data, 'weights': mc_weights, 'xsection_pb': mc_xsec_pb, 'xsection_pb_err' : mc_xsec_pb, 'acceptance': acceptance}


def read_hepdata(dataset, all_obs, cdir=None):
    """
    Read in HepData files
    """

    datapath = f"{cdir}/{dataset['DATAPATH']}"
    hepdata  = {}

    # Over different observables    
    for k in range(len(dataset["HIST"])):
            
        filename = datapath + dataset["HIST"][k]['FILE']
        OBS      = dataset["HIST"][k]['OBS']
        
        ### Read out HEPDATA
        if   dataset['TYPE'] == 'HEPDATA':
            hepdata[OBS]          = ReadHEPData(filename)
            hepdata[OBS]['scale'] = dataset["HIST"][k]['SCALE']
        
        ### Raw scalar values in ascii file
        elif dataset['TYPE'] == 'RAW_SCALAR':
            
            with open(filename) as f:
                content = f.readlines()
            content     = [float(x.strip()) for x in content]

            x           = np.array(content)
            bins        = np.linspace(dataset["HIST"][k]['XMIN'], dataset["HIST"][k]['XMAX'], dataset["HIST"][k]['NBINS']+1)

            counts, errs, bins, cbins = iceplot.hist(x=x, bins=bins, density=all_obs[OBS]['density'])

            hepdata[OBS] = {}
            hepdata[OBS]['y']        = counts
            hepdata[OBS]['y_err']    = errs
            hepdata[OBS]['bins']     = bins
            hepdata[OBS]['x']        = cbins
            hepdata[OBS]['binwidth'] = bins[1:] - bins[0:-1]
            hepdata[OBS]['scale']    = dataset["HIST"][k]['SCALE']
        
        else:
            raise Exception(__name__ + f'.read_hepdata: Unknown dataset "TYPE" = {dataset["TYPE"]}')
        
        ### Set the same binning for the observables as HepDATA has
        all_obs[OBS]['bins']     = hepdata[OBS]['bins']
        all_obs[OBS]['binwidth'] = hepdata[OBS]['binwidth']
        all_obs[OBS]['xlim']     = (np.min(all_obs[OBS]['bins']), np.max(all_obs[OBS]['bins']))
    
    
    ### Remove observables not found in the steering card
    poplist = []
    for OBS in all_obs.keys():
        if OBS not in hepdata.keys():
            poplist.append(OBS)

    for OBS in poplist:
        all_obs.pop(OBS)

    return hepdata, all_obs


def doublet2linear(binedges):
    """
    Convert doublet bins to linear array assuming
    bins are continous (no holes)
    """
    bins = []
    bins.append(binedges[0][0])
    for i in range(len(binedges)):
        bins.append(binedges[i][1])

    return np.array(bins)


def doublet2binwidth(binedges):
    """
    Convert doublet bins to binwidth
    """
    binwidth = np.zeros(len(binedges))
    for i in range(len(binedges)):
        binwidth[i] = binedges[i][1] - binedges[i][0]

    return binwidth


def orthogonal_error(data):
    """
    Total orthogonal error combined from different sources
    """
    total_error = None
    for key in data.keys():
        if total_error is None: # Not yet initialized
            total_error = np.zeros(len(data[key]))

        total_error += (data[key][:,1] - data[key][:,0])**2

    return np.sqrt(total_error)


def ReadHEPData(filename):
    """
    Read in HEPData .csv file
    """

    xval        = 0.0
    xval_low    = 0.0
    xval_high   = 0.0
    yval        = 0.0

    yerr_stat_p = 0.0
    yerr_stat_m = 0.0

    yerr_syst_p = 0.0
    yerr_syst_m = 0.0

    yerr_sexp_p = 0.0
    yerr_sexp_m = 0.0

    x   = []
    y   = []
    err = []

    file   = open(filename, 'r')
    Lines  = file.readlines()
    found_header = False

    dataset = {
        'x'          : [],
        'binedges'   : [],
        'y'          : [],
        'y_err'      : [],
        'y_err_stat' : [],
        'y_err_syst' : [],
        'y_err_sexp' : []
    }

    # Read line-by-line
    i = 0
    for line in Lines:

        thisline = line#.strip()
        result   = thisline.split(',')

        if len(result) == 10:

            if not found_header:
                found_header = True
                continue

            xval        = float(result[0])
            xval_low    = float(result[1])
            xval_high   = float(result[2])
            yval        = float(result[3])
            yerr_stat_p = float(result[4])
            yerr_stat_m = float(result[5])
            yerr_syst_p = float(result[6])
            yerr_syst_m = float(result[7])
            yerr_sexp_p = float(result[8])
            yerr_sexp_m = float(result[9])

        else:
            continue

        # Push values
        dataset['x'].append(xval)
        dataset['binedges'].append([xval_low, xval_high])
        dataset['y'].append(yval)
        dataset['y_err_stat'].append([yerr_stat_m, yerr_stat_p])
        dataset['y_err_syst'].append([yerr_syst_m, yerr_syst_p])
        dataset['y_err_sexp'].append([yerr_sexp_m, yerr_sexp_p])


    if (len(dataset['x']) == 0):
        raise Exception(__name__ + f".ReadHEPData: Found 0 row entries from the file {filename}!")

    # Change datatypes
    for key in dataset.keys():
        dataset[key] = np.array(dataset[key])

    dataset['bins']     = doublet2linear(dataset['binedges'])
    dataset['binwidth'] = doublet2binwidth(dataset['binedges'])
    dataset['y_err']    = orthogonal_error({'stat':dataset['y_err_stat'], 'syst':dataset['y_err_syst'], 'sexp':dataset['y_err_sexp']})

    return dataset
