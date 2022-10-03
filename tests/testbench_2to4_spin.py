# GRANIITTI comparisons
#
# Run with: pytest ./tests/{filename} -s
# 
# (c) 2017-2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import subprocess
import pytest
import os
import pyjson5 as json5

from test_helpers import *


#@pytest.mark.parametrize("POMLOOP", [False, True])
def test_2to4_spin(POMLOOP, NEVENTS, GENERATE, ANALYZE):
    
    cdir      = os.getcwd()    
    inputcard = f"{cdir}/tests/LHC_TEVATRON_RHIC/ATLAS17_4pi_0.json"

    RESOS = ['f0_1500', 'f2_1525']
    PID   = '211,-211, 211,-211'

    if GENERATE:

        # Pure phase space
        for RES in RESOS:
            cmd = f"./bin/gr -p 'PP[RES]<F> -> rho(770)0 > {{pi+ pi-}} rho(770)0 > {{pi+ pi-}} @RES{{{RES}:1}} @SPINDEC:false' -o '{RES}_@SPINDEC:false' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w true"
            print(cmd)
            os.system(cmd)

        for RES in RESOS:
            if 'f_2' in RES:
                EXTRA = '@R[{RES}]{{JZ0:0.0, JZ1:0.0, JZ2:1.0}} @FRAME:CM'
            else:
                EXTRA = ''
            cmd = f"./bin/gr -p 'PP[RES]<F> -> rho(770)0 > {{pi+ pi-}} rho(770)0 > {{pi+ pi-}} @RES{{{RES}:1}} {EXTRA}' -o '{RES}' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w true"
            print(cmd)
            os.system(cmd)

    if ANALYZE:

        cmd = f"python python/iceshot --obs '4body' --title 'PP[RES]<F> -> rho(770)0 > {{pi+ pi-}} rho(770)0 > {{pi+ pi-}}' --density --pid '[[{PID}], [{PID}], [{PID}], [{PID}]]' --hepmc3 '{RESOS[0]}_@SPINDEC:false' '{RESOS[1]}_@SPINDEC:false' '{RESOS[0]}' '{RESOS[1]}' --unit nb"
        print(cmd)
        os.system(cmd)

if __name__ == "__main__":
    test_2to4_spin(POMLOOP=False, NEVENTS=int(1E5), GENERATE=False, ANALYZE=True)
