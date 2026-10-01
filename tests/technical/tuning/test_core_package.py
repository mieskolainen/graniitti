# Python source package and resource tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import os
import shutil
import subprocess
import sys

from core import resource


# Exercise resources and discovery outside the simulator working directory
def test_package_outside_simulator(tmp_path):
    script = '''import sys
from core import resource
from core.stats.hist import hist
from core.tune.drivers.registry import driver_entries
from core.tune.parameters.space import normalize_param_space
assert set(driver_entries()) >= {"GRANIITTI", "PANDORA"}
assert not {"torch", "ray", "pyHepMC3"}.intersection(sys.modules)
assert not any(name.endswith(".driver") for name in sys.modules)
assert "submit" not in sys.modules
assert resource("tune/settings/ampfit.json").is_file()
assert normalize_param_space({"x": {"type": "float", "lower": 0, "upper": 1}})["x"]["upper"] == 1
assert resource("templates/iceweb.html").is_file()
'''
    subprocess.run([sys.executable, "-c", script], cwd=tmp_path, check=True)
    subprocess.run([sys.executable, "-m", "core.icetune", "--help"], cwd=tmp_path, check=True, capture_output=True)


# Discover supplied drivers and resources without installed Python packages
def test_source_package_without_site_packages(tmp_path):
    shutil.copytree(resource(""), tmp_path / "core", ignore=shutil.ignore_patterns("__pycache__"))
    script = '''from importlib.metadata import distributions
from core import resource
from core.tune.drivers.registry import driver_entries
assert not any(dist.metadata["Name"] == "core" for dist in distributions())
entries = driver_entries()
assert entries["GRANIITTI"].value == "core.tune.drivers.graniitti.driver:GraniittiDriver"
assert entries["PANDORA"].value == "core.tune.drivers.pandora.driver:PandoraDriver"
assert resource("tune/settings/ampfit.json").is_file()
'''
    subprocess.run(
        [sys.executable, "-S", "-c", script], cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(tmp_path)}, check=True)
