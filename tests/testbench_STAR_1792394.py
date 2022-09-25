# Python steering for STAR \\sqrt{s} = 0.2 TeV pp -> p + X + p
# 
# Run all simulations and produce comparison plots with HEPData
# for Low-x and Diffraction 2022.
#
# Run e.g. with "pytest filename.py --POMLOOP False --NEVENTS 1000"
# 
# (c) 2017-2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import argparse

from test_helpers import *


def datasets():
    """
    Labelling and filenames
    """
    sets,labels,pids = {},{},{}
    
    # pipi
    sets['pipi']    = ['pipi', 'pipi_0--1', 'pipi_1--1-5', 'pipi_1-5--inf', 'pipi_less_90', 'pipi_more_90']
    labels['pipi']  = ['$\\pi^+\\pi^-$',
                       '$\\pi^+\\pi^-, M_{{X}} < 1$ GeV',
                       '$\\pi^+\\pi^-, 1 < M_{{X}} < 1.5$ GeV',
                       '$\\pi^+\\pi^-, M_{{X}} > 1.5$ GeV',
                       '$\\pi^+\\pi^-, \\Delta \\phi_{{pp}} < 90^{\\circ}$',
                       '$\\pi^+\\pi^-, \\Delta \\phi_{{pp}} > 90^{\\circ}$']
    pids['pipi']    = f'[[211,-211]]'

    # KK
    sets['KK']      = ['KK', 'KK_less_90', 'KK_more_90']
    labels['KK']    = ['$K^+K^-$',
                       '$K^+K^-, \\Delta \\phi_{{pp}} < 90^{\\circ}$',
                       '$K^+K^-, \\Delta \\phi_{{pp}} > 90^{\\circ}$']
    pids['KK']      = f'[[321,-321]]'

    # ppbar
    sets['ppbar']   = ['ppbar', 'ppbar_less_90', 'ppbar_more_90']
    labels['ppbar'] = ['$p\\bar{p}$',
                       '$p\\bar{p}, \\Delta \\phi_{{pp}} < 90^{\\circ}$',
                       '$p\\bar{p}, \\Delta \\phi_{{pp}} > 90^{\\circ}$']
    pids['ppbar']   = f'[[2212,-2212]]'

    return sets,labels,pids

def test_generate_analyze(POMLOOP, NEVENTS, GENERATE=True, ANALYZE=True):

    # Event generation parameters
    WEIGHTED = 'true' # (lowercase string)! Weighted events is faster (we do not need unweighted events here)
    
    # Get datasets
    HEPDATA  = 'STAR_1792394'
    SETS, LABELS, PIDS = datasets()

    # Loop over
    if GENERATE:
        for channel in SETS.keys():
            for i in range(len(SETS[channel])):
                cmd = f"./bin/gr -i gencard/{HEPDATA}_{SETS[channel][i]}.json -w {WEIGHTED} -l {'true' if POMLOOP else 'false'} -n {NEVENTS}"

                execute(cmd)
    
    if ANALYZE:
        for channel in SETS.keys():
            for i in range(len(SETS[channel])):
                cmd  = f"python python/iceshot --hepmc3 {HEPDATA}_{SETS[channel][i]} --hepdata dataset_{HEPDATA}_{SETS[channel][i]} --pid '{PIDS[channel]}' "
                cmd += f" --mclabel 'GRANIITTI MC'"
                cmd += f" --datalabel 'STAR ({LABELS[channel][i]}), $\\sqrt{{s}}=0.2$ TeV'"
                cmd += f" --chi2"

                execute(cmd, expect="[iceshot: done]")

if __name__ == "__main__":
    test_generate_analyze(POMLOOP=True, NEVENTS=1000000)
