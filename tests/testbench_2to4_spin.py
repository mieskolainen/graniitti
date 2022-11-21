# GRANIITTI comparisons for cascaded X > A > {A1 A2} B > {B1 B2} spin correlations
#
# Run with: pytest tests/{filename} -s -rP --NEVENTS 10000 --GENERATE 1 --ANALYZE 1
# 
# (c) 2017-2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import subprocess
import pytest
import copy
import os
import pyjson5 as json5

from test_helpers import *


def test_2to4_spin(POMLOOP, NEVENTS, GENERATE, ANALYZE):
    
    cdir      = os.getcwd()    
    inputcard = f"{cdir}/tests/LHC_TEVATRON_RHIC/full_solid_angle.json"

    PROCESS   = 'PP[RES]<F> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}'
    # Increase via @syntax the f0_1710 mass so we are not too near the 2-body threshold of 2 x rho
    # (otherwise phase-space squeezing non-flatness effect in the angular cos(theta),phi distributions)
    MASS      = ['@R[f0_1710]{M:2.7,W:0.15}', '@R[f0_1710_neg_P]{M:2.7,W:0.15}']
    RESOS     = ['f0_1710', 'f0_1710_neg_P']
    OUTPUT    = ['f0_2700', 'f0_2700_neg_P']
    PID       = '[211,-211, 211,-211]'
    WEIGHTED  = 'false'
    
    if GENERATE:

        # Pure phase space
        for i in range(len(RESOS)):
            cmd = f"./bin/gr -p '{PROCESS} {MASS[i]} @RES{{{RESOS[i]}:1}} @SPINDEC:false' -o '{OUTPUT[i]}@SPINDEC:false' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w {WEIGHTED}"
            print(cmd)
            execute(cmd)
        
        # Spin correlations
        for i in range(len(RESOS)):
            cmd = f"./bin/gr -p '{PROCESS} {MASS[i]} @RES{{{RESOS[i]}:1}}' -o '{OUTPUT[i]}' -i {inputcard} -l {'true' if POMLOOP else 'false'} -h 0 -n {NEVENTS} -w {WEIGHTED}"
            print(cmd)
            execute(cmd)


    if ANALYZE:

        # ----------------------------
        # Generate HepMC3 input txt
        HEPMC3_txt = ''
        for i in range(len(OUTPUT)):
            HEPMC3_txt += f"'{OUTPUT[i]}@SPINDEC:false' "
        for i in range(len(OUTPUT)):
            HEPMC3_txt += f"'{OUTPUT[i]}' "

        # Create PID txt
        PID_txt = copy.deepcopy(PID)
        for i in range(len(RESOS)*2 - 1):
            PID_txt += f', {PID}'
        PID_txt = f"'{PID_txt}'"
        # ----------------------------
        
        cmd  = f"python python/iceshot --obs 4body --cuts none"
        cmd += f" --title '{PROCESS}' --density --unit nb"
        cmd += f" --pid {PID_txt} --hepmc3 {HEPMC3_txt}"

        print(cmd)
        execute(cmd, expect="[iceshot: done]")
