# GRANIITTI comparison with data
#
# Run with: pytest ./tests/{filename} -s
# 
# (c) 2017-2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import subprocess
import pytest
import os
import pyjson5 as json5

from test_helpers import *


#@pytest.mark.parametrize("POMLOOP", [False, True])
def test_cepdata(POMLOOP, NEVENTS, MODELPARAM):
    
    print(f'POMLOOP={POMLOOP} NEVENTS={NEVENTS} MODELPARAM={MODELPARAM}')

    cdir      = os.getcwd()
    
    # Filenames under /fitcard/
    filenames = ['STAR_1792394_pipi_less_90', 'STAR_1792394_KK_less_90', 'STAR_1792394_ppbar_less_90',
                 'STAR_1792394_pipi_more_90', 'STAR_1792394_KK_more_90', 'STAR_1792394_ppbar_more_90',
                 'STAR_1792394_pipi', 'STAR_1792394_KK', 'STAR_1792394_ppbar']
    

    # Apply fast screening (absorption) factor at \sqrt{s} = 200 GeV
    # (applied to plots)
    S2_factor = 0.3 if not POMLOOP else 1.0
    
    
    ### Loop over all datasets
    for name in filenames:
            
        datasetfile = f'dataset_{name}'
        with open(f'{cdir}/fitcard/{datasetfile}.json', 'r') as file:
            json_data = json5.load(file)
        
        pid       = json_data['PID']
        
        inputcard = f"{cdir}/{json_data['GENCARD'][0]}"
        output    = f'{name}_POMLOOP_{POMLOOP}'
        gridfile  = f'{cdir}/vgrid/{output}.vgrid'       
        
        # ----------------------------------------------------------------
        
        ### Compute MC grids
        cmd = f"./bin/gr -i {inputcard} -m {MODELPARAM} -l {'true' if POMLOOP else 'false'} -h 0 -n 0 -o {output}"
        execute(cmd)

        ### Generate events with pre-computed MC grids
        cmd = f"./bin/gr -i {inputcard} -m {MODELPARAM} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w true -d {gridfile} -o {output}"
        execute(cmd)
        
        ### Produce comparison figures
        cmd = f"python ./python/iceshot --pid '[{pid}]' --hepmc3 {output} --hepdata {datasetfile} --mcscale {S2_factor} --chi2 --unit nb --mclabel 'GRANIITTI'"
        execute(cmd, expect="[iceshot: done]")
