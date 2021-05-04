# Unit test MC integration & grids
#
# Run with: pytest ./tests/{filename} -s
#
# (c) 2017-2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import pytest
import os
import numpy as np
import json
import matplotlib.pyplot as plt
import pathlib

from test_helpers import *



@pytest.mark.parametrize("FS", ['pi+ pi-', 'K+ K-', 'pi+ pi- pi+ pi-'])
def test_multicore(FS):
    """
    Test multithreaded sampling consistency of different amplitudes.
    Should yield the same results, if no multithreading related problems.
    """
    print(__name__ + ': Generating events ...')
    
    sigma = {'C': [], 'F': []}

    for CORES in [1, 16, 64]:
        for PS in sigma.keys():
            
            cmd = f'./bin/gr -i ./tests/run_compare_PS/card.json -p "PP[CON]<{PS}> -> {FS}" -h 0 -n 0 -c {CORES} -o PS_{PS}_CORES_{CORES}'
            execute(cmd)
            
            ## Read integration results
            with open(f'./vgrid/PS_{PS}_CORES_{CORES}.vgrid') as f:
                data = json.load(f)
                sigma[PS].append(data['STAT']['sigma'])

    # To Numpy
    for PS in sigma.keys():
        sigma[PS] = np.array(sigma[PS])

    # ** Self-comparison **
    for PS1 in sigma.keys():
        for PS2 in sigma.keys():
            for i in range(len(sigma[PS2])):
                assert sigma[PS1] == pytest.approx(sigma[PS2][i], rel=0.05)


def test_vgrid():

    ## Execute MC generator
    print(__name__ + ': Generating events ...')

    cmd = f"make -j4 && ./bin/gr -i ./input/test.json -p 'PP[CON]<F> -> pi+ pi-' -w true -l false -h 0 -n 0"
    execute(cmd)

    ## Read integration output
    with open('./vgrid/test.vgrid') as f:
      data = json.load(f)
    
    print(data.keys())
    print(data['VEGAS'].keys())

    pathlib.Path('./testfigs/').mkdir(parents=True, exist_ok=True)

    keys = ['xmat', 'f2mat', 'fmat', 'region', 'rvec', 'xcache']

    for key in keys:
        fig,ax = plt.subplots()
        plt.title(f'{key}')
        plt.plot(data['VEGAS'][key])
        fig.savefig(f"./testfigs/testbench_vgrid_{key}.pdf", bbox_inches='tight')

    assert True
