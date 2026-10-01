# Shared real HEPData input for tuning initialization tests
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import shutil
from pathlib import Path

import pyjson5 as json5
import pytest
from core.io import jsonref


# Copy the analysis bundle while retaining the unmodified published data
@pytest.fixture
def data_card(tmp_path):
    root = Path(__file__).resolve().parents[3]
    bundle = tmp_path / 'icepack/PHOTOPROD/LHCb_2825384'
    shutil.copytree(root / 'icepack/PHOTOPROD/LHCb_2825384', bundle)
    shutil.copytree(root / 'icepack/_common', tmp_path / 'icepack/_common')
    source = jsonref.JsonReader(json5.load)
    source.read(root / 'icepack/PHOTOPROD/LHCb_2825384/jpsi/dataset.json')
    for original in source.documents:
        copied = tmp_path / original.relative_to(root)
        copied.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, copied)
    path = bundle / 'jpsi/dataset.json'
    dataset = json5.loads(path.read_text())
    dataset['datapath'] = str(root / dataset['datapath'])
    path.write_text(json.dumps(dataset))
    return path
