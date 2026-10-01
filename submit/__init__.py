# Repository resources for simulation campaign submission
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

CAMPAIGN_DIR = Path(__file__).resolve().parent
SHELL_DIR = CAMPAIGN_DIR / "shell"


# Resolve the JSON tune card selected by a repository campaign
def campaign_source(name, catalog=None):
    import yaml

    path = Path(catalog or CAMPAIGN_DIR / "campaigns.yml").resolve()
    campaign = yaml.safe_load(path.read_text())["campaigns"][name]
    return str((path.parent / campaign["tunesetup"]).resolve())
