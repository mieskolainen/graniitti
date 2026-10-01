# Installed tools for simulation, inference and data analysis
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from importlib.resources import files
from pathlib import Path

__version__ = "0.1.0"
__RELEASE__ = "beta"
__AUTHOR__  = "Mikael Mieskolainen (mikael.mieskolainen@cern.ch)"

# Locate an installed configuration, template or execution helper
def resource(name: str) -> Path:
    return Path(str(files("core").joinpath(name)))
