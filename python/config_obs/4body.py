# Histogram definitions for the observables
# 
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import iceobs


obs_4body_DeltaM_AB = {

'tag'     : '4body_DeltaM_AB',

# Axis limits
'xlim'    : (-1.0, 1.0),
'ylim'    : None,
'xlabel'  : r'$\Delta M_{AB}$',
'ylabel'  : r'$d\sigma/d \Delta M_{AB}$',
'units'   : {'x': r'GeV', 'y': r'pb'},
'label'   : r'Invariant mass difference',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-1.0, 1.0, 50),
'density' : False,
 
# Function to calculate
'func'    : iceobs.proj_1D_4body_DeltaM_AB
}


obs_4body_Deltacos_12 = {

'tag'     : '4body_Deltacos_12',

# Axis limits
'xlim'    : (-2, 2),
'ylim'    : None,
'xlabel'  : r'$\Delta \cos \theta_{12}$',
'ylabel'  : r'$d\sigma/d\Delta \cos \theta_{12}$',
'units'   : {'x': r'unit', 'y': r'pb'},
'label'   : r'Angle',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-2, 2, 40),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_4body_Deltacos_12
}

obs_M = {

'tag'     : 'M',

# Axis limits
'xlim'    : (1.0, 3.0),
'ylim'    : None,
'xlabel'  : r'$M$',
'ylabel'  : r'$d\sigma/dM$',
'units'   : {'x': r'GeV', 'y': r'pb'},
'label'   : r'Invariant mass',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(1, 3, 50),
'density' : False,
 
# Function to calculate
'func'    : iceobs.proj_1D_M
}

obs_4body_cos1_dotprod = {

'tag'     : '4body_cos1_dotprod',

# Axis limits
'xlim'    : (-1, 1),
'ylim'    : None,
'xlabel'  : r'$\cos \theta_{1}$ [boost + dot product]',
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
'func'    : iceobs.proj_1D_4body_cos1_dotprod
}


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
