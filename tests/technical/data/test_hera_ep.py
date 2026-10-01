# Validate integrated HERA ep measurements and their published uncertainties
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

from pathlib import Path

import numpy as np
import pytest
from core.io import readers, steering
from core.plot import plot

from icepack.PHOTOPROD._common import ep_reader
from tests.technical.support.hepdata import original

ROOT = Path(__file__).resolve().parents[3]


# Preserve the original ep column and asymmetric errors through the plotting API
@pytest.mark.parametrize(("card", "record", "table", "group"), [
    ("H1_451266", 451266, 2, 1), ("ZEUS_473522", 473522, 1, 0),
])
def test_ep_measurement(card, record, table, group):
    dataset, path = steering.load_dataset(str(ROOT / "icepack/PHOTOPROD" / card / "dataset.json"), cdir=ROOT)
    entry = dataset["sets"][0]
    observables, _ = steering.load_observables(entry["obs"], dataset_path=path, cdir=ROOT)
    data, observables = readers.read_hepdata(
        dataset=entry, datapath=dataset["datapath"], datatype=dataset["type"], all_obs=observables,
        cdir=str(ROOT), reader=dataset["reader"], dataset_path=path,
    )
    source = ROOT / f"HEPData/PHOTOPROD/HEPData-ins{record}-v1-json/Table{table}.json"
    _, values, errors = original(source, group=group)
    result = data["cross_section"]
    assert result["y"] == pytest.approx(values)
    total = np.sqrt(sum(pair**2 for pair in errors.values()))
    assert result["y_err_down"] == pytest.approx(total[:, 0])
    assert result["y_err_up"] == pytest.approx(total[:, 1])
    weights = values * result["scale"] * 1e12
    mc = plot.histmc(
        mcdata={"data": {"cross_section": result["x"]}, "weights": weights, "xsection_pb": weights.sum()},
        obs=observables,
    )["cross_section"]["hdata"]
    assert mc.counts_scaled == pytest.approx(weights)
    with pytest.raises(ValueError, match="file_filter"):
        ep_reader.read(str(source))


# Apply invariant elasticity before and after proton fragmentation in either beam direction
@pytest.mark.parametrize("direction", [-1, 1])
@pytest.mark.parametrize("fragment", [False, True])
@pytest.mark.parametrize("target_mass", [2.0, 12.0])
def test_psi_elasticity(direction, fragment, target_mass):
    from icepack.PHOTOPROD.H1_451266 import cuts
    from tests.technical.data.test_hera_dissociation import photo_event

    event, _ = photo_event(direction, fragment, target_mass)
    event.cut_param = cuts.cut_param
    assert cuts.cut_func(event) == (target_mass < 10.0)
