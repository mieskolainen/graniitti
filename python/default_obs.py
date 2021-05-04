# Histogram definitions for the observables
# 
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import iceobs


obs_M = {

'tag'     : 'M',

# Axis limits
'xlim'    : (0, 3.0),
'ylim'    : None,
'xlabel'  : r'$M$',
'ylabel'  : r'$d\sigma/dM$',
'units'   : {'x': r'GeV', 'y': r'pb'},
'label'   : r'Invariant mass',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(0, 3, 50),
'density' : False,
 
# Function to calculate
'func'    : iceobs.proj_1D_M
}


obs_Abs_t1t2 = {

'tag'     : 'Abs_t1t2',

# Axis limits
'xlim'    : (0, 2),
'ylim'    : None,
'xlabel'  : r'$|t_1+t_2|$',
'ylabel'  : r'$d\sigma/d|t_1+t_2|$',
'units'   : {'x': r'GeV$^{2}$', 'y': r'pb'},
'label'   : r'Sum of momentum transfers squared',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(0, 2, 50),
'density' : False,
 
# Function to calculate
'func'    : iceobs.proj_1D_Abs_t1t2
}


obs_Rap = {

'tag'     : 'Rap',

# Axis limits
'xlim'    : (-4,4),
'ylim'    : None,
'xlabel'  : r'$Y$',
'ylabel'  : r'$d\sigma/dY$',
'units'   : {'x': r'unit', 'y': r'pb'},
'label'   : r'Rapidity',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(-3, 3, 50),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_Rap
}


obs_dPhi_pp = {

'tag'     : 'dPhi_pp',

# Axis limits
'xlim'    : (0, 180),
'ylim'    : None,
'xlabel'  : r'$\Delta\phi_{pp}$',
'ylabel'  : r'$d\sigma/d\Delta\phi_{pp}$',
'units'   : {'x': r'deg', 'y': r'pb'},
'label'   : r'Transverse angle separation',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(0, 180, 50),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_dPhi_pp
}
"""

"""
obs_costheta_CS = {

'tag'     : 'costheta_CS',

# Axis limits
'xlim'    : (-1, 1),
'ylim'    : None,
'xlabel'  : r'cos$(\theta)_{CS}$',
'ylabel'  : r'$d\sigma/dcos(\theta)_{CS}$',
'units'   : {'x': r'unit', 'y': r'pb'},
'label'   : r'Angle in Collins-Soper frame',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(-1, 1, 50),
'density' : False,
 
# Function to calculate
'func'    : iceobs.proj_1D_costheta_CS
}


obs_phi_CS = {

'tag'     : 'phi_CS',

# Axis limits
'xlim'    : (0, 180),
'ylim'    : None,
'xlabel'  : r'$\phi_{CS}$',
'ylabel'  : r'$d\sigma/d\phi$',
'units'   : {'x': r'deg', 'y': r'pb'},
'label'   : r'Angle in Collins-Soper frame',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.7, 1.3),

# Default histogramming
'bins'    : np.linspace(0, 180, 50),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_phi_CS
}




