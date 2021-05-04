# Global tests (slow process tests)
#
# 
# Run with: pytest ./tests/{filename} -s
# 
# (c) 2017-2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import subprocess
import pytest

from test_helpers import *


# ------------------------------------------------------------------------

def test_PP():
    """
    Test soft processes
    """
    proclist = ['PP[RESTENSOR]', 'PP[RES+CONTENSOR]']

    cmd = []
    for proc in proclist:
        for fs in ['pi+ pi-', 'K+ K-']:
            cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p '{proc}<F> -> {fs} @RES{{rho_770:1, f0_980:1, f2_1270:1}}' -o test_PP_{proc}")
    execute(cmd)


def test_gg():
    """
    Test Durham QCD processes
    """
    cmd = []

    proclist = ['gg[CON]<C> -> g g', 'gg[FLUX]<C> -> g g'];
    for proc in proclist:
        cmd.append(f"./bin/gr -i ./tests/processes/gg2gg.json -w true -l false -h 0 -n 100 -p '{proc}' -o test_gg")
    
    proclist = ['gg[CON]<F> -> pi+ pi-', 'gg[CON]<F> -> pi+ pi-']
    for proc in proclist:
        cmd.append(f"./bin/gr -i ./tests/processes/gg2chic0.json -w true -l false -h 0 -n 100 -p '{proc}' -o test_chic0")
    
    proclist = ['gg[chic(0)]<F> &> pi+ pi-', 'gg[chic(0)]<F> &> K+ K-']
    for proc in proclist:
        cmd.append(f"./bin/gr -i ./tests/processes/gg2MMbar.json -w true -l false -h 0 -n 100 -p '{proc}' -o test_MMbar")
    
    execute(cmd)


def test_pomloop():
    cmd = []
    for screening in ['false', 'true']:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l {screening} -h 0 -n 0 -p 'PP[CON]<C> -> pi+ pi-'")
    execute(cmd)


def test_lhapdf():
    """
    Test setting lhapdf
    """
    cmd = []
    for pdf in ['CT10nlo', 'MMHT2014lo68cl']:
        cmd.append(f"./bin/gr -i ./tests/processes/gg2gg.json -w true -l false -h 0 -n 0 -p 'gg[CON]<C> -> g g' -o gg2gg -f hepmc3 -q {pdf}")
    execute(cmd)


