# HepMC3 and HEPData input/output functions
#
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from importlib import import_module

from iceobs import *

from pyHepMC3 import HepMC3 as hepmc3


# Object wrapper (for icecache)
class ice_event:
    def __init__(self, evt, pid):
        self.evt = evt
        self.pid = pid

def generate_partitions(totalsize, chunksize):
    """
    Generate partition indices
    """

    # Special case
    if totalsize <= chunksize:
        yield (0, totalsize-1)
        return # Needs return here

    # General case
    N = int(np.ceil(totalsize/chunksize))
    remainder = totalsize - chunksize * N
    a = 0
    for i in range(N):
        b = a + chunksize + (i < remainder)
        yield (a, b - 1)
        a = b


def generate_chunks(lst, n):
    """Yield successive n-sized chunks from the list
    """
    for i in range(0, len(lst), n):
        yield lst[i:i + n]


def read_hepmc3(hepmc3file, all_obs, pid, cuts=None, maxevents=None, chunk_range=[0, int(1E12)], xsmode='header'):
    """
    Read in a single HepMC3 file and compute observables

    Args:
    
        xsmode = 'header' or 'sample'
    """
    mc_data = {}
    for OBS in all_obs.keys():
        mc_data[OBS] = []

    mc_weights       = []
    mc_xsec_pb       = 0
    mc_xsec_pb_err   = 0

    mc_weights_tot   = []
    mc_trials        = 0

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
    t    = time.time()
    file = hepmc3.ReaderAscii(hepmc3file)
    k    = -1

    while not file.failed():
        k  += 1
        evt = hepmc3.GenEvent()
        file.read_event(evt)

        if file.failed(): # end-of-file
            break

        # -------------------------------
        # This after the header reader
        if k < chunk_range[0]:
            continue
        if k > chunk_range[1]:
            break
        # -------------------------------

        # Event weight for each event
        w = float(evt.weights()[0])
        mc_weights_tot.append(w)

        # Get GenCrossSection info attributes
        mc_xsec_pb     = float(evt.cross_section().xsec())
        mc_xsec_pb_err = float(evt.cross_section().xsec_err())
        mc_trials      = float(evt.cross_section().get_attempted_events())

        event = ice_event(evt=evt, pid=pid)

        # -------------------------------
        # Apply additional fiducial cuts here ...
        if not cut_func(event):
            continue;
        # -------------------------------

        # Construct all observables
        for OBS in all_obs.keys():
            try:
                value = all_obs[OBS]['func'](event=event)
                mc_data[OBS].append(value)
            except:
                continue;

        # Event weight for passed events
        mc_weights.append(w)

    mc_weights_tot = np.array(mc_weights_tot)
    mc_weights     = np.array(mc_weights)

    # --------------------------------------------------------------------
    # Compute cross section from the sample weights (works with weighted events)
    if xsmode == 'sample':
        SCALE          = 1E12
        N              = mc_trials
        W              = np.sum(mc_weights_tot)
        W2             = np.sum(mc_weights_tot**2)

        mc_xsec_pb     = SCALE * W/N
        mc_xsec_pb_err = SCALE * np.sqrt((W2/N - (W/N)**2)/N);
    # --------------------------------------------------------------------
    
    # Count the python user cuts fiducial acceptance
    acceptance     = np.sum(mc_weights) / np.sum(mc_weights_tot)

    elapsed        = time.time() - t
    print(f'{hepmc3file}')
    print(f'Before cuts: mc_xs_pb: {mc_xsec_pb:.3g} +- {mc_xsec_pb_err:.3g} | events read = {len(mc_weights_tot)} | wsum = {np.sum(mc_weights_tot):0.3E} | chunk = {chunk_range}')
    
    if acceptance < 1.0:
        mc_xsec_pb     *= acceptance
        mc_xsec_pb_err *= acceptance
    
    print(f'After  cuts: mc_xs_pb: {mc_xsec_pb:.3g} +- {mc_xsec_pb_err:.3g} | events read = {len(mc_weights)} | wsum = {np.sum(mc_weights):0.3E} | cuts <{cuts}> acceptance = {acceptance:0.3E}')
    
    return {'data': mc_data, 'weights': mc_weights, 'xsection_pb': mc_xsec_pb, 'xsection_pb_err' : mc_xsec_pb_err, 'acceptance': acceptance}


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
        
        ### Read out HEPData
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
