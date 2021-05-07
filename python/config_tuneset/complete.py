# MC tuning steering file
#
# Be careful with parameter names, i.e. no extra white space etc.
#
# (c) 2021 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.


import numpy as np
from ray import tune

from config_tuneset.default import *

mc_steer['FRESHGRID'] = True
