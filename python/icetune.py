# IO & MC tuning steering functions
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import numpy as np
import pyjson5 as json5
import json as json
import os
import subprocess
import re

import pathlib
from tqdm import tqdm
import socket
from importlib import import_module
from pprint import pprint
import sys


# Command line arguments
from glob import glob
from braceexpand import braceexpand

import matplotlib
matplotlib.use('Agg') # Important for multithreaded applications
from matplotlib import pyplot as plt


import iceio
import iceplot


def setpaths(cdir, libdir, PYTHON_VERSION="3.8"):

    try:
        LD_LIBRARY_PATH = os.environ['LD_LIBRARY_PATH']
    except:
        LD_LIBRARY_PATH = ''

    try:
        PYTHONPATH = os.environ['PYTHONPATH']
    except:
        PYTHONPATH = ''
    
    LD_LIBRARY_PATH = f"{libdir}/HEPMC3/lib:{libdir}/HEPMC3/lib64:{libdir}/LHAPDF/lib:{libdir}/LHAPDF/lib64:{LD_LIBRARY_PATH}"
    PYTHONPATH      = f"{cdir}/python:{libdir}/HEPMC3/lib/python{PYTHON_VERSION}/site-packages:{PYTHONPATH}"
    
    os.environ['LD_LIBRARY_PATH'] = LD_LIBRARY_PATH
    os.environ['PYTHONPATH']      = PYTHONPATH
    
    sys.path.append(PYTHONPATH) # This is crucial with Ray!

    return LD_LIBRARY_PATH, PYTHONPATH


def get_observables(module_name):
    """ 
    Read out content from module under 'obs'
    """
    MyConfig = import_module('config_obs.' + module_name)

    all_obs = {}
    for key in MyConfig.__dict__.keys():
        if 'obs_' in key:
            obj = MyConfig.__dict__[key]
            all_obs[obj['tag']] = obj

    return all_obs


def get_ip():
    """
    Return IP address
    """
    s = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    try:
        # doesn't even have to be reachable
        s.connect(('10.255.255.255', 1))
        IP = s.getsockname()[0]
    except Exception:
        IP = '127.0.0.1'
    finally:
        s.close()
    return IP


def compute(thread_id, datacards, obs_module, mc_steer, compare_steer, cdir=None, best_cost=None):
    """
    Compute MC sample and compare with HEPData over all datacards input
    """
    
    print(__name__ + f'.compute: Running with thread_id = {thread_id} ...')

    def recompute(inputfile, output):
        cmd = f"{cdir}/bin/gr -h 0 -i {inputfile} -o {output} -n 0 -l {'true' if mc_steer['POMLOOP'] else 'false'}"
        #result = os.system(cmd)
        result = subprocess.check_output(cmd, shell=True, text=True)
    
    if cdir is None:
        cdir = os.getcwd()
        print(__name__ + f'.compute: cdir : {cdir}')

    ### Loop over all datasets of the card
    datasets    = []
    datacards = list(braceexpand(datacards.replace(" ", "")))

    for k in range(len(datacards)):

        # Do wildcard expansion
        #print(path)
        this_ = glob(cdir + '/datacard/' + datacards[k])
        #print(this_)
        
        for j in range(len(this_)):
            with open(this_[j], 'r') as file:
                json_data = json5.load(file)
            datasets.append(json_data)
    
    #pprint(datasets)
    chi2 = []
    ndf  = []

    if len(datasets) == 0:
        raise Exception(__name__ + f".compute: Error: did not find any input datasets -- check your input")


    ### Loop over different datasets
    plot_index = 0

    for i in range(len(datasets)):

        if datasets[i]['ACTIVE']:
            
            # ========================================================================
            ### Read input and initialize
            
            inputfile = cdir + datasets[i]["GENCARD"][0] # Take the first one here only

            with open(inputfile, 'r') as file:
                output = json5.load(file)['GENERALPARAM']['OUTPUT']

            output   = output + f"_POMLOOP_{mc_steer['POMLOOP']}"
            
            # ------------------------------------------------------------------------
            # Integration grid reset requested
            if mc_steer['XSMODE'] == 'reset':
                output = f'{output}_{thread_id}'
            # ------------------------------------------------------------------------

            gridfile = cdir + '/vgrid/' + output + '.vgrid'

            # If not yet initialized
            if not os.path.exists(gridfile):
                print(f'lossfunc: Running integration grid initialization for {inputfile}')
                recompute(inputfile=inputfile, output=output)

            # It exists, check it
            else:
                with open(gridfile, 'r') as file:
                    griddata = json5.load(file)

                if griddata['POMLOOP'] != mc_steer['POMLOOP']:
                    print('vgrid file with non-matching POMLOOP -- re-computing')
                    recompute(inputfile=inputfile, output=output)

            # ========================================================================
            ### Generate events
            
            if thread_id == 0:
                tunename = f'TUNE0'
            else:
                tunename = f'TUNE-icetune-{thread_id}'
            hepmc3output = f'{output}_{thread_id}'

            # Generation command
            cmd         = f"{cdir}/bin/gr -h 0 -w true -m {tunename} -d {gridfile} -i {inputfile} -o {hepmc3output } -n {mc_steer['NEVENTS']} -l {'true' if mc_steer['POMLOOP'] else 'false'}"
            #print(__name__ + f".compute: Generating events: {cmd}")
            result      = subprocess.check_output(cmd, shell=True, text=True)
            
            outputfile  = f'{cdir}/output/' + hepmc3output
            hepmc3file  = outputfile + '.hepmc3'
            pid         = datasets[i]['PID']
            cuts        = 'config_cuts.' + datasets[i]['CUTS']
            
            # ========================================================================
            ### Get observables
            
            all_obs = get_observables(obs_module)
            if len(all_obs) == 0: raise Exception(__name__ + '.compute: Observables not found!')

            # Read in data observables
            hepdata, all_obs = iceio.read_hepdata(dataset=datasets[i], all_obs=all_obs, cdir=cdir)

            # Read in MC observables
            xsmode  = 'sample' if mc_steer['XSMODE'] == 'sample' else 'header'
            mcdata           = iceio.read_hepmc3(hepmc3file=hepmc3file, all_obs=all_obs, pid=pid, cuts=cuts, xsmode=xsmode)     
            
            
            # ========================================================================
            ### Histogramming
            
            # Only one MC source here
            mc   = iceplot.histmc(mcdata=mcdata, all_obs=all_obs, density=compare_steer['density'], scale=mc_steer['kfactor'], color=iceplot.colors(0), label='GRANIITTI')
            
            # Only one data source here
            data = iceplot.histhepdata(hepdata=hepdata, all_obs=all_obs, density=compare_steer['density'], label='Data')


            # ========================================================================
            ### Change y-axis labels (for plots)
            
            if compare_steer['density']:
                all_obs = iceplot.change2density_labels(all_obs=all_obs)
            

            # ========================================================================
            ### Cost per histogram

            ### Loop over observables
            for OBS in all_obs.keys():

                ### ------------------------------------------------------------------
                ### Chi2 for this observable
                chi2_  = iceplot.chi2_cost(h_mc=mc[OBS]['hdata'], h_data=data[OBS]['hdata'])
                ndf_   = len(mc[OBS]['hdata'].counts) # Count here only the number of bins (not parameters)
                
                chi2.append(chi2_)
                ndf.append(ndf_)
                # --------------------------------------------------------------------
                
                # ========================================================================
                ### Plotting function

                if compare_steer['plot'] and best_cost is not None: # Only if plot output asked

                    if chi2[plot_index] < best_cost[plot_index]:

                        hist_objs = [data[OBS], mc[OBS]] # Create a list, data first!

                        for yscale in ['linear', 'log']:
                            fig, ax = iceplot.superplot(hist_objs, observable=all_obs[OBS], yscale=yscale)

                            # Create path and save
                            ax[0].set_title(f'$\\chi^2$ / ndf = {chi2_:0.1f} / {ndf_:0.0f} = {chi2_/ndf_:0.1f}')
                            fullpath = f'{cdir}/figs/icetune--{output}'
                            pathlib.Path(fullpath).mkdir(parents=True, exist_ok=True)
                            fig.savefig(f"{fullpath}/hplot__{all_obs[OBS]['tag']}_{yscale}.pdf", bbox_inches='tight')
                            plt.close()            

                plot_index += 1

            # ========================================================================
            # ** Remove temporary outputfile **
            for file in ['hepmc3']:
                subprocess.check_output(f'rm {outputfile}.{file}', shell=True, text=True)
            
    # ------------------------------------------------------------
    # ** Remove temporary tune folder **
    if thread_id != 0:
        tunefolder = f'{cdir}/modeldata/TUNE-icetune-{thread_id}'
        subprocess.check_output(f'rm {tunefolder} -f -r', shell=True, text=True)
    # ------------------------------------------------------------

    return np.array(chi2), np.array(ndf)


def update_all(thread_id, param, cdir=None):
    """
    Update parameter json files
    """

    if cdir is None:
        cdir = os.getcwd()
        print(__name__ + f'.update_all: cdir : {cdir}')

    # --------------------------------------------------------------------
    # Create tune folder
    cmd = f'mkdir {cdir}/modeldata/TUNE-icetune-{thread_id} -p'
    result = subprocess.check_output(cmd, shell=True, text=True)

    # Copy content from default tune
    cmd = f'cp {cdir}/modeldata/TUNE0/* {cdir}/modeldata/TUNE-icetune-{thread_id}/'
    result = subprocess.check_output(cmd, shell=True, text=True)
    # --------------------------------------------------------------------

    path = f'{cdir}/modeldata/TUNE-icetune-{thread_id}'

    # General parameters
    orig_gen_values = update_general_param(path=path, param=param)

    # Resonance parameters
    orig_res_values = update_resonance_param(path=path, param=param)

    # Branching parameters, TBD!
    #orig_br_values = update_branching_param(thread_id=thread_id, param=param, cdir=cdir)

    return {**orig_gen_values, **orig_res_values}


def find_str_between(s, start, end):
    return s[s.find(start)+len(start):s.rfind(end)]


def substring_after(s, delim):
    return s.partition(delim)[2]


def substring_before(s, delim):
    return s.split(delim)[0]


def string_to_index(s, delim):
    """ Returns tuple with indices """
    return tuple([int(item) for item in s.split(delim)])


def update_block(parname, blockname, new_value, json_data):
    """
    Updates parameters in json
    
    Returns:
        old_value from json
    """
    # Vector parameter: parname[index]
    if '[' in parname and ']' in parname:

        index   = string_to_index(s=find_str_between(s=parname, start='[', end=']'), delim=',')
        parname = substring_before(s=parname, delim='[')

        if parname in json_data[blockname]:

            if   len(index) == 1: # Vector
                old_value = json_data[blockname][parname][index[0]]
                json_data[blockname][parname][index[0]]           = new_value

                return old_value

            elif len(index) == 2: # Matrix
                old_value = json_data[blockname][parname][index[0]][index[1]]
                json_data[blockname][parname][index[0]][index[1]] = new_value

                return old_value
            else:
                raise Exception(__name__ + f'.update_block: Too many [] indices with {parname}')
        else:
            raise Exception(__name__ + f'.update_block: Unknown parameter {parname}')

    # Scalar parameter: parname
    else:
        if parname in json_data[blockname]:
            old_value = json_data[blockname][parname]
            json_data[blockname][parname] = new_value

            return old_value
        else:
            raise Exception(__name__ + f'.update_block: Unknown parameter {parname}')


def update_general_param(path, param):
    """
    Update general parameters.
    
    Tune parameter string format is: CLASS|parameter[index]
    """

    filename = f'{path}/GENERAL.json'

    # 1. READ
    with open(filename, 'r') as file:
        json_data = json5.load(file)

    orig_values = {}

    # 2. UPDATE PARAMETERS
    for par in param.keys():

        # Extract NAME out
        CLASS   = substring_before(s=par, delim='|')

        ### Resonance or branching parameters not treated here
        if CLASS == 'RES' or CLASS == 'BR': continue

        parname   = substring_after(s=par, delim='|')
        orig_values[par] = update_block(parname=parname, blockname=f'PARAM_{CLASS}', new_value=param[par], json_data=json_data)

    # 3. WRITE
    with open(filename, 'w') as file:
        json.dump(json_data, file, indent=2)

    return orig_values


def update_resonance_param(path, param):
    """
    Update resonance parameters.
    
    Tune parameter string format is: RES|resonance:parameter[index]
    """
    RES = set()

    # Extract unique resonance names (use a set)
    for par in param.keys():
        if "RES|" in par:
            RES.add(find_str_between(s=par, start='|', end=':'))

    orig_values = {}

    # Loop over resonances found associated with the free parameters
    for res in RES:

        # ----------------------------------------------------------
        # Open up file
        filename = f'{path}/RES_{res}.json'

        # 1. READ
        with open(filename, 'r') as file:
            json_data = json5.load(file)
        # ----------------------------------------------------------

        # Loop over all tune parameters
        for par in [j for j in param.keys() if f'RES|{res}' in j]:

            # Update the parameter
            parname = substring_after(s=par, delim=':')
            orig_values[par] = update_block(parname=parname, blockname='PARAM_RES', new_value=param[par], json_data=json_data)

        # 3. WRITE
        with open(filename, 'w') as file:
            json.dump(json_data, file, indent=2)

    return orig_values


def update_branching_param(path, param):

    """
    # TBD
    
    filename = f'{path}/BRANCHING.json'
    
    # 1. READ
    with open(filename, 'r') as file:
        json_data = json5.load(file)
    
    # 2. UPDATE PARAMETERS (if any)
    #
    # ...
        
    # 3. WRITE
    with open(filename, 'w') as file:
        json.dump(json_data, file, indent=2)
    
    """
    return

