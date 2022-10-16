# Histogram definitions for the observables
# 
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import iceobs

obs_4body_cos1 = {

'tag'     : '4body_cos1',

# Axis limits
'xlim'    : (-1, 1),
'ylim'    : None,
'xlabel'  : r'$\cos \theta_{1}$',
'ylabel'  : r'$d\sigma/d\cos \theta_{1}$',
'units'   : {'x': r'unit', 'y': r'pb'},
'label'   : r'Angle',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-1, 1, 40),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_4body_cos1
}


obs_4body_cos2 = {

'tag'     : '4body_cos2',

# Axis limits
'xlim'    : (-1, 1),
'ylim'    : None,
'xlabel'  : r'$\cos \theta_{2}$',
'ylabel'  : r'$d\sigma/d\cos \theta_{2}$',
'units'   : {'x': r'unit', 'y': r'pb'},
'label'   : r'Angle',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-1, 1, 40),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_4body_cos2
}


obs_4body_phi12 = {

'tag'     : '4body_phi12',

# Axis limits
'xlim'    : (-180, 180),
'ylim'    : None,
'xlabel'  : r'$[\phi_{1} + \phi_{2}]_{atan2}$',
'ylabel'  : r'$d\sigma/d[\phi_{1} + \phi_{2}]_{atan2}$',
'units'   : {'x': r'deg', 'y': r'pb'},
'label'   : r'Sum angle',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-180, 180, 40),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_4body_phi12
}

