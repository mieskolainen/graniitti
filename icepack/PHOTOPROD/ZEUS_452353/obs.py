# HERA photoproduction comparison definitions
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from icepack.PHOTOPROD._common.exclusive_jpsi import abs_t
from icepack.PHOTOPROD.H1_1798511.obs import obs_M  # noqa: F401

obs_t = {'tag':'t', 'xlim':None, 'ylim':None, 'xlabel':r'$|t|$', 'ylabel':r'$d\sigma/d|t|$',
         'units':{'x':r'GeV$^2$', 'y':'pb'}, 'label':'Proton momentum transfer',
         'ylim_ratio':(0.0,2.0), 'ytick_ratio_step':0.5, 'bins':None, 'density':False, 'func':abs_t}
