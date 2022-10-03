# Histogram definitions for the observables
# 
# (c) 2022 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
import iceobs


obs_phi12 = {

'tag'     : 'phi12',

# Axis limits
'xlim'    : (-360, 360),
'ylim'    : None,
'xlabel'  : r'$\phi_{1} - \phi_{2}$',
'ylabel'  : r'$d\sigma/(d\phi_{1} - d\phi_{2})$',
'units'   : {'x': r'deg', 'y': r'pb'},
'label'   : r'Sum angle',
'figsize' : (4, 3.75),

# Ratio
'ylim_ratio' : (0.0, 2.0),
'ytick_ratio_step': 0.5,

# Default histogramming
'bins'    : np.linspace(-360, 360, 40),
'density' : False,

# Function to calculate
'func'    : iceobs.proj_1D_phi12
}


