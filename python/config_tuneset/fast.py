# MC tuning steering file (FAST TEST)
#
# Be careful with parameter names, i.e. no extra white space etc.
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
from ray import tune


# Event generation setup
mc_steer = {
  'NEVENTS' : 30000,     # At least ~ 30k to 50k events per parameter space trial
  'POMLOOP' : False,     # Screening on/off (set False for fast test trials, True for precision tuning)
  'XSMODE'  : 'sample',  # Cross section normalization mode:
                         # -- 'reset'   compute full cross section initialization (vgrid file) for every parameter trial
                         # -- 'vgrid'   use once pre-computed value from .vgrid file (use this if parameters under tune do not impact the integrated xs)
                         # -- 'sample'  (default) compute new cross-section from the event sample weights, otherwise as 'vgrid' (use this if parameters impact the integrated xs)
  
  'kfactor' : 1.0        # MC cross-section scale factor
}


# ** Data input collection **
# Use glob wildcard / braces here

SET = 0


if   SET == 0:
  datacards = '{ \
    dataset_STAR_1792394_pipi_1--1.5.json, \
    dataset_STAR_1792394_KK.json, \
    dataset_STAR_1792394_KK_less_90.json, \
    dataset_STAR_1792394_pipi.json, \
    dataset_STAR_1792394_pipi_less_90.json }'

  if mc_steer['POMLOOP'] == False:
    mc_steer['kfactor'] = 0.3 # Simple fast approximation

elif SET == 1:
  datacards = '{ \
    dataset_STAR_1792394_KK.json, \
    dataset_STAR_1792394_KK_less_90.json}'
  
  if mc_steer['POMLOOP'] == False:
    mc_steer['kfactor'] = 0.3

elif SET == 2:
  datacards = '{dataset_STAR_1792394_pipi*, dataset_STAR_1792394_KK*}'
  
  if mc_steer['POMLOOP'] == False:
    mc_steer['kfactor'] = 0.3

elif SET == 3:
  datacards = '{dataset_ALICE*, dataset_ATLAS*}'

  if mc_steer['POMLOOP'] == False:
    mc_steer['kfactor'] = 0.15

elif SET == 4:
  datacards = '{ \
    dataset_ALICE*, dataset_ATLAS*, \
    dataset_STAR_1792394_pipi_1--1.5.json, \
    dataset_STAR_1792394_KK.json, \
    dataset_STAR_1792394_KK_less_90.json, \
    dataset_STAR_1792394_pipi.json, \
    dataset_STAR_1792394_pipi_less_90.json }'

else:
  raise Exception(__name__ + 'Unknown dataset chosen !')


# ------------------------------------------------------------------------



# General parameters under GENERAL.json
#
# SYNTAX EXAMPLE: 'CLASS|a0[2]', where [index] are used only with vector (matrix) parameters
#
config = {
  
  ### Pomeron trajectory
  #'REGGE|a0[0]':        tune.uniform(1.05, 1.15),
  'REGGE|ap[0]':        tune.uniform(0.05, 0.30),

  ### Offshell meson form factors
  #'REGGE|offshellFF':   "EXP",
  #'REGGE|b_EXP':        tune.uniform(0.1, 3.0),

  'REGGE|offshellFF':   "POW",
  'REGGE|b_POW':        tune.uniform(0.1, 1.5),

  #'REGGE|offshellFF':  "OREAR",
  #'REGGE|a_OREAR':     tune.uniform(0.1, 3.0),
  #'REGGE|b_OREAR':     tune.uniform(0.1, 3.0)
}


# Resonance parameters under RES_resname.json(s)
# 
# SYNTAX EXAMPLE: 'RES|resname:param[1,1]', where [] are used only with vector (matrix) parameters


# Set here the resonances which are subject to the fit (process cards can contain also other resonances)
free_resonances = ['f0_500', 'f0_980', 'f2_1270', 'f0_1500', 'f2_1525', 'f0_1710']
#free_resonances = ['f0_1500', 'f2_1525', 'f0_1710']
#free_resonances = []


## Masses and widths
"""
config['RES|f0_500:M'] = tune.uniform(0.3, 0.7)
config['RES|f0_500:W'] = tune.uniform(0.3, 0.6)
"""

"""
config['RES|f0_980:M'] = tune.uniform(0.950, 1.0)
config['RES|f0_980:W'] = tune.uniform(0.05, 0.20)
"""

# Production coupling amplitudes
for res in free_resonances:
  #if res not in ['rho_770']:
  config[f'RES|{res}:g_A']   = tune.uniform(0.0001, 1.0)

# Production coupling phases
for res in free_resonances:
  #if res in ['f0_980', 'f2_1270']:
  config[f'RES|{res}:g_phi'] = tune.uniform(-3.14159, 3.14159) # Continuum phase
  #config[f'RES|{res}:g_phi'] = tune.choice(np.array([-3.14159, -1.57, 0.0, 1.57])) # Quantized phase
  #config[f'RES|{res}:g_phi'] = tune.choice(np.linspace(-np.pi, np.pi, 17)[:-1]) # Quantized phase


raytune_setup = {
   'config'      : config,
   'metric'      : "chi2",
   'mode'        : "min",
   'num_samples' : 5000
}

