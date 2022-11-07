# GRANIITTI comparisons for cascaded X > A > {A1 A2} B > {B1 B2} spin correlations
#
# Run with: pytest tests/{filename} -s
# 
# (c) 2017-2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import subprocess
import pytest
import copy
import os
import pyjson5 as json5

from test_helpers import *


#@pytest.mark.parametrize("POMLOOP", [False, True])
def test_2to4_spin(POMLOOP, NEVENTS, GENERATE=True, ANALYZE=True):
    
    cdir      = os.getcwd()    
    inputcard = f"{cdir}/tests/LHC_TEVATRON_RHIC/full_solid_angle.json"

    PROCESS  = 'PP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}'
    RESOS    = ['f0_1710', 'f0_1710_neg_P']
    PID      = '[211,-211, 211,-211]'
    WEIGHTED = 'false'
    
    if GENERATE:

        # Pure phase space
        for RES in RESOS:
            cmd = f"./bin/gr -p '{PROCESS} @RES{{{RES}:1}} @SPINDEC:false' -o '{RES}@SPINDEC:false' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w {WEIGHTED}"
            print(cmd)
            os.system(cmd)
        
        # Spin correlations
        for RES in RESOS:
            if 'f_2' in RES:
                EXTRA = '@R[{RES}]{{JZ0:0.0, JZ1:0.0, JZ2:1.0}} @FRAME:CM'
            else:
                EXTRA = ''
            cmd = f"./bin/gr -p '{PROCESS} @RES{{{RES}:1}} {EXTRA}' -o '{RES}' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w {WEIGHTED}"
            print(cmd)
            os.system(cmd)

    if ANALYZE:

        # ----------------------------
        # Generate steering txt
        RESOS_txt = ''
        for i in range(len(RESOS)):
            RESOS_txt += f"'{RESOS[i]}@SPINDEC:false' "
        for i in range(len(RESOS)):
            RESOS_txt += f"'{RESOS[i]}' "

        PID_txt = copy.deepcopy(PID)
        for i in range(len(RESOS)*2 - 1):
            PID_txt += f', {PID}'
        PID_txt = f"'{PID_txt}'"
        # ----------------------------

        cmd  = f"python python/iceshot --obs 4body --cuts null"
        cmd += f" --title '{PROCESS}' --density --unit nb"
        cmd += f" --pid {PID_txt} --hepmc3 {RESOS_txt}"
        
        print(cmd)
        os.system(cmd)

if __name__ == "__main__":
    test_2to4_spin(POMLOOP=False, NEVENTS=int(1E4), GENERATE=True, ANALYZE=True)
