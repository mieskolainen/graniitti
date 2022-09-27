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
    sets,labels,pids,names = {},{},{},{}
    

    # pipi
    sets['pipi']    = ['pipi', 'pipi_0--1', 'pipi_1--1-5', 'pipi_1-5--inf', 'pipi_less_90', 'pipi_more_90']
    labels['pipi']  = ['$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts)',
                       '$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts & $M_{{X}} < 1$ GeV)',
                       '$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts & $1 < M_{{X}} < 1.5$ GeV)',
                       '$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts & $M_{{X}} > 1.5$ GeV)',
                       '$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts & $\\Delta \\phi_{{pp}} < 90^{\\circ})$',
                       '$pp \\rightarrow p + \\pi^+\\pi^- + p$ (fid. cuts & $\\Delta \\phi_{{pp}} > 90^{\\circ})$']
    pids['pipi']    = f'[[211,-211], [211,-211]]'
    names['pipi']   = f'pi+ pi-'


    # KK
    sets['KK']      = ['KK', 'KK_less_90', 'KK_more_90']
    labels['KK']    = ['$pp \\rightarrow p + K^+K^- + p$ (fid. cuts)',
                       '$pp \\rightarrow p + K^+K^- + p$ (fid. cuts & $\\Delta \\phi_{{pp}} < 90^{\\circ})$',
                       '$pp \\rightarrow p + K^+K^- + p$ (fid. cuts & $\\Delta \\phi_{{pp}} > 90^{\\circ})$']
    pids['KK']      = f'[[321,-321], [321,-321]]'
    names['KK']     = f'K+ K-'


    # ppbar
    sets['ppbar']   = ['ppbar', 'ppbar_less_90', 'ppbar_more_90']
    labels['ppbar'] = ['$pp \\rightarrow p + p\\bar{p} + p$ (fid. cuts)',
                       '$pp \\rightarrow p + p\\bar{p} + p$ (fid. cuts & $\\Delta \\phi_{{pp}} < 90^{\\circ})$',
                       '$pp \\rightarrow p + p\\bar{p} + p$ (fid. cuts & $\\Delta \\phi_{{pp}} > 90^{\\circ})$']
    pids['ppbar']   = f'[[2212,-2212], [2212,-2212]]'
    names['ppbar']  = f'p+ p-'


    return sets,labels,pids,names


def test_generate_analyze(POMLOOP, NEVENTS, GENERATE=True, ANALYZE=True):

    # Event generation parameters
    WEIGHTED = 'true' # (lowercase string)! Weighted events is faster (we do not need unweighted events here)
    
    # Get datasets
    HEPDATA  = 'STAR_1792394'
    SETS, LABELS, PIDS, NAMES = datasets()

    # Loop over
    if GENERATE:
        for channel in SETS.keys():
            for i in range(len(SETS[channel])):

                # Process as steered by the input card
                cmd = f"./bin/gr -i gencard/{HEPDATA}_{SETS[channel][i]}.json -w {WEIGHTED} -l {'true' if POMLOOP else 'false'} -n {NEVENTS}"
                
                print(cmd)
                execute(cmd)

                # Continuum only process
                PROCESS = f'PP[CON]<F> -> {NAMES[channel]}'
                OUTPUT  = f'{HEPDATA}_{SETS[channel][i]}_CON' 
                cmd     = f"./bin/gr -p '{PROCESS}' -o {OUTPUT} -i gencard/{HEPDATA}_{SETS[channel][i]}.json -w {WEIGHTED} -l {'true' if POMLOOP else 'false'} -n {NEVENTS}"
                
                print(cmd)
                execute(cmd)

    if ANALYZE:
        for channel in SETS.keys():
            for i in range(len(SETS[channel])):

                cmd  = f"python python/iceshot --hepmc3 {HEPDATA}_{SETS[channel][i]} {HEPDATA}_{SETS[channel][i]}_CON --hepdata dataset_{HEPDATA}_{SETS[channel][i]} --pid '{PIDS[channel]}'"
                cmd += f" --mclabel 'GRANIITTI [RES+CON]' 'GRANIITTI [CON]'"
                cmd += f" --datalabel 'STAR, $\\sqrt{{s}}=0.2$ TeV'"
                cmd += f" --title '{LABELS[channel][i]}'"
                #cmd += f" --maxevents 10000"
                
                print(cmd)
                execute(cmd, expect="[iceshot: done]")

if __name__ == "__main__":
    test_generate_analyze(POMLOOP=True, NEVENTS=int(1E5), GENERATE=True, ANALYZE=True)
