# Global tests (fast process tests)
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
# INPUT-OUTPUT

def test_input_output():
    """
    Test input-output procedures
    """

    cmd = []
    cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -o CON")
    cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[RES+CON]<F> -> pi+ pi-' -o RES_CON")
    cmd.append(f"./bin/gr -i ./input/test.json -d ./vgrid/CON.vgrid -w true -l false -h 0 -n 1000 -p 'PP[CON]<F> -> pi+ pi-' -o CON -f hepmc3")
    cmd.append(f"./bin/gr -i ./input/test.json -d ./vgrid/RES_CON.vgrid -w true -l false -h 0 -n 1000 -p 'PP[RES+CON]<F> -> pi+ pi-' -o RES_CON -f hepmc3")
    execute(cmd)
    
    cmd = "python ./python/iceshot --hepmc3 CON RES_CON --maxevents 1000 --pid '[[211,-211], [211,-211]]'"
    execute([cmd], expect='[iceshot: done]')


# ------------------------------------------------------------------------
# GENERAL


def test_vgrid():
    """
    Test reading pre-computed MC grid
    """
    cmd = []

    command = f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -o vgrid_test"
    cmd.append(f"{command}")
    cmd.append(f"{command} -d ./vgrid/vgrid_test.vgrid")
    execute(cmd)


def test_output():
    """
    Test output name
    """
    cmd = []
    for output in ['output_A', 'output_B']:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -o {output} -o output_test_{output}")
    execute(cmd)


def test_format():
    """
    Test output formats
    """
    cmd = []
    for f in ['hepmc2', 'hepmc3', 'hepevt']:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -f {f} -o format_test_{f}")
    execute(cmd)


def test_cores():
    """
    Test cores
    """
    cmd = []
    for cores in [0, 16, 128, 1024, 4096]:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -c {cores} -o cores_test_{cores}")
    execute(cmd)


def test_nevents():
    """
    Test number of events
    """
    cmd = []
    for n in [0, 100]:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -p 'PP[CON]<F> -> pi+ pi-' -n {n}  -o nevents_test_{n}")
    execute(cmd)


def test_weighted():
    """
    Test (un)weighted event generation
    """
    cmd = []
    for weighted in ['false', 'true']:
        cmd.append(f"./bin/gr -i ./input/test.json -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi-' -w {weighted} -o weighted_test_{weighted}")
    execute(cmd)


def test_modelparam():
    """
    Test changing modeltune
    """
    cmd = []
    for modelparam in ['TUNE0', 'TUNE0']:
        cmd.append(f"./bin/gr -i ./input/test.json -w false -m {modelparam} -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -o modelparam_test_{modelparam}")
    execute(cmd)


# ---------------------------------------------------------
# PROCESS

def test_process():
    cmd = []
    for proc in ['PP[RES+CON]<F> -> pi+ pi-', 'PP[CON]<C> -> K+ K-', \
                 'PP[CON]<C> -> rho(770)0 > {pi+ pi-} rho(770)0 > {pi+ pi-}', \
                 'PP[RES+CON]<F> -> pi+ pi- @RES{f0_500:0,rho_770:1,f0_980:1,f2_1270:1} @R[f0_980]{M:0.98,W:0.065}']:
        
        cmd.append(f"./bin/gr -i ./input/test.json -w true -h 0 -n 0 -p '{proc}' -o process_test")
    execute(cmd)


def test_energy():
    cmd = []
    for energy in ['7000', '13000']:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -h 0 -n 0 -p 'PP[CON]<C> -> pi+ pi-' -e {energy} -o energy_test_{energy}")
    execute(cmd)


def test_nstars():
    """
    Test forward excitation (N*)
    """
    cmd = []
    for excite in [0,1,2]:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 0 -p 'PP[CON]<F> -> pi+ pi-' -s {excite} -o excite_test_{excite}")
    execute(cmd)


def test_hist():
    """
    Test fast embedded histograms
    """
    cmd = []
    for screening in ['true', 'false']:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -h 0 -n 0 -p 'PP[CON]<C> -> pi+ pi-' -l {screening} -o hist_test_{screening}")
    execute(cmd)

def test_rndseed():
    """
    Test random seeding
    """
    cmd = []
    for seed in [123, 456]:
        cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<C> -> pi+ pi-' -r {seed} -o rndseed_test_{seed}")
    execute(cmd)


# ------------------------------------------------------------------------
# On-the-flight parameter syntax


def test_otf():

    cmd = []

    # Tree sampling
    clist  = ['@FLATAMP:1', '@FLATMASS2:true', '@OFFSHELL:7']
    
    for OF in clist: cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi- {OF}' -o test1")
    execute(cmd)

    # Mass and width
    clist  = ['@PDG[992]{M:350.0, W:5.0}']
    clist += ['@R[f0_980]{M:0.990, W:0.065}']

    for OF in clist: cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi- {OF}' -o test2")
    execute(cmd)

    # Active resonances
    clist  = ['@RES{f0_980:1, f2_1270:1}']
    
    for OF in clist: cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi- {OF}' -o test3")
    execute(cmd)

    # Spin parameters
    clist  = ['@SPINGEN:true', '@SPINDEC:true']
    clist += ['@FRAME:CS', '@FRAME:HX', '@FRAME:CM']
    clist += ['@R[f2_1270]{JZ0:0.5, JZ1:0.0, JZ2:0.5}']
    clist += ['@JMAX:4']
    
    for OF in clist: cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi- {OF}' -o test4")
    execute(cmd)

    # Tensor Pomeron parameters
    clist  = ['@R[f0_980]{g0:1.0, g1:0.2}']
    clist += ['@R[f0_1270]{g0:1.0, g1:0.1, g2:0.4, g3:0.3, g4:0.2, g5:0.1, g6:0.0}']

    for OF in clist: cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p 'PP[CON]<F> -> pi+ pi- {OF}' -o test5")
    execute(cmd)


# ------------------------------------------------------------------------
# PROCESSES

def test_PP():
    """
    Test soft processes
    """
    proclist = ['OP[RES]', 'PP[CON-]', 'PP[CON]', 'PP[RES+CON]', 'PP[RESHEL]', 'PP[RES]']

    cmd = []
    for proc in proclist:
        for fs in ['pi+ pi-', 'K+ K-']:
            cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p '{proc}<F> -> {fs}' -o test_{proc}")
    execute(cmd)


def test_yP():
    """
    # Test gamma-Pomeron processes
    """
    proclist = ['yP[RESTENSOR]', 'yP[RES]']
    RES = '@RES{rho_770:1}'

    cmd = []
    for proc in proclist:
        for fs in ['pi+ pi-']:
            cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p '{proc}<F> -> {fs} {RES}' -o test_yP")
    execute(cmd)


def test_yy_DZ():
    """
    Test gamma-gamma processes
    """
    proclist = ['yy_DZ[CON]<P>', 'yy_DZ[FLUX]<P>', 'yy_LUX[CON]<P>']
    
    cmd = []
    for proc in proclist:
        for fs in ['e+ e-', 'mu+ mu-', "tau+ tau-"]:
            cmd.append(f"./bin/gr -i ./tests/processes/yy2ee.json -w true -l false -h 0 -n 100 -p '{proc} -> {fs}' -o test_yy_DZ")
    execute(cmd)


def test_yy():
    """
    Test gamma-gamma processes
    """
    proclist = ['yy[CON]', 'yy[FLUX]', 'yy[QED]', 'yy[RES]']
    
    RES = '@RES{rho_770:1}'

    # Higgs
    execute([f"./bin/gr -i ./tests/processes/yy2Higgs.json -w true -l false -h 0 -n 100 -p 'yy[Higgs]<F> &> b b~' -o test_Higgs_bbar"])
    
    # Others
    cmd = []
    for proc in proclist:
        for fs in ['e+ e-', 'mu+ mu-']:
            cmd.append(f"./bin/gr -i ./input/test.json -w true -l false -h 0 -n 100 -p '{proc}<F> -> {fs} {RES}' -o test_yy")
    execute(cmd)

