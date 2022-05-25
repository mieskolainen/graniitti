# Loop over fiducial measurements and construct MC & data comparison tex array,
# together with asserted tests comparing cross sections
# 
# 
# Run with: pytest ./tests/{filename} -s
# 
# m.mieskolainen@imperial.ac.uk, 2021


import pyjson5 as json5
import os
import subprocess
import pytest

from pytest import approx
import numpy as np


BARNS = {3: 'mb', 6: '$\\mu\\$b', 9: 'nb', 12: 'pb', 15: 'fb'}


def print_tex_header():

    print('\n\n\n')
    print('\\begin{table}')
    print('\\begin{center}')
    print('\\begin{tabular}{|lllr|rrc|}')
    print('\\hline')
    print(' & & & \\tiny{MEASUREMENT} & & \\tiny{GRANIITTI} & \\\\')
    print('\\hline')
    print('$\\sqrt{s}$ TeV & channel & \\tiny{PHASE SPACE CUTS} & value $\\pm$ stat $\\pm$ syst & $\\sigma_{S^2}$ & $\\sigma_0$ & $\\langle S^2 \\rangle$ \\\\')
    print('\\hline')


def print_tex_footer():
    print('\\hline')
    print('\\end{tabular}')
    print('\\caption{Measurements versus GRANIITTI.}')
    print('\\label{table:xstable}')
    print('\\end{center}')
    print('\\end{table}')
    print('')


#@pytest.mark.parametrize("POMLOOP", [False, True])
def test_exloop(POMLOOP, datacard='/HEPData/LHC_TEVATRON_RHIC.json', ASSERT_ON=False):
    
    print(f'POMLOOP={POMLOOP} datacard={datacard} ASSERT_ON={ASSERT_ON}')

    print_tex_header()
    
    # ====================================================================
    # Set screening true / false

    POMLOOP_values = ['false', 'true'] if POMLOOP else ['false']

    # ====================================================================
    # Read in JSON data
    
    cdir    = os.getcwd()
    SYS_FOLDERS = ['/gencard/', '/tests/LHC_TEVATRON_RHIC', '/tests/processes']
    
    with open(f'{cdir}/{datacard}', 'r') as file:
        j = json5.load(file)
    
    # Loop over JSON data
    for key in j.keys():
        
        CITE = j[key]['BIBTEX']
        if CITE is None: CITE = ''
        
        ### Run MC
        for folder in SYS_FOLDERS: # Check different folders

            inputcard = cdir + f'/{folder}' + '/' + key + '.json'
            if os.path.isfile(inputcard):

                for P in POMLOOP_values:
                    cmd    = f'{cdir}/bin/gr -i {inputcard} -l {P} -h 0 -n 0 -o {key}_POMLOOP_{P}'
                    result = subprocess.check_output(cmd, shell=True, text=True)
            else:
                continue
        
        # MC x-sections here
        XS          = {}
        XS['true']  = 0
        XS['false'] = 0
        
        ### Collect MC output values
        for P in POMLOOP_values:
            with open(f'{cdir}/vgrid/{key}_POMLOOP_{P}.vgrid', 'r') as file:
                vgrid = json5.load(file)

                # Read-out integrated cross-section
                XS[P] = float(vgrid['STAT']['sigma']) * (10**j[key]['UNIT']) 
        
        # Integrated screening (absorption) factor computed
        if POMLOOP:
            S2 = XS['true'] / XS['false']
        else:
            S2 = 0
        
        ### Print output
        m            = j[key]['MEASUREMENT']
        MEASUREMENT  = f'${m[0]} \\pm {m[1]} \\pm {m[2]}$' if m[0] is not None else ''
        CUTS_et_CITE = j[key]['CUTS'] + ' \\cite{' + CITE + '}' if m[0] is not None else j[key]['CUTS'] 
        MC_XS        = f"{XS['true']:0.3g} & {XS['false']:0.3g} {BARNS[j[key]['UNIT']]}"
        
        print(f"{j[key]['SQRTS']:5s} & {j[key]['CHANNEL']:25s} & {CUTS_et_CITE:60s} & {MEASUREMENT:30s} & {MC_XS:20s} & {S2:0.3g} \\\\")
        
        
        ### ** Autocomparison of the measurement and screened MC result **
        if m[0] is not None and POMLOOP and ASSERT_ON:
            error = np.sqrt(m[1]**2 + m[2]**2) # Measurement error stat & syst in quadrature
            abs_tolerance = error * j[key]['ASSERT_NSIGMA']

            assert XS['true'] == approx(m[0], abs=abs_tolerance)
    
    # ====================================================================
    
    print_tex_footer()

