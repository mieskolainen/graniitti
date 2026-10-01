# Direct MP polarization through native amplitudes and the icetune driver
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import math
import shutil
import subprocess
from itertools import product
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pyjson5
import pytest
import torch
from core import resource
from core.io.serialize import load_json_file, write_json_file
from core.tune.drivers.graniitti.ampfit.amplitude import AmplitudeBank, input_card
from core.tune.drivers.graniitti.driver import GraniittiDriver
from core.tune.drivers.graniitti.tunesetup import card
from core.tune.drivers.graniitti.tunesetup.domains import angular_row_name
from core.tune.optimizers.ampfit.config import load_settings
from core.tune.parameters.tools import raw_phase_key, spherical_angle_key

from tests.technical.integration.test_ampfit import dataset_card

ROOT = Path(__file__).resolve().parents[3]


# Compare differentiated spin coordinates with native amplitudes on the same VEGAS events
@pytest.mark.parametrize(
    ("frame", "screened", "resonance", "basis", "derivative"),
    [(frame, screened, resonance, "auto_min_L", False)
     for frame, screened, resonance in product(("CM", "CS", "HX"), (False, True),
                                               ("f2_1270", "rho3_1690", "f2_1270_yy"))]
    + [("CS", screened, "f2_1270", basis, True)
       for screened, basis in product((False, True),
                                      ("auto_min_S", "auto_equal_ls", "auto_equal_helicity", "g_ls"))]
    + [("HX", screened, "rho3_1690", "helicity", True) for screened in (False, True)],
)
def test_polarization_reweighting(tmp_path, frame, screened, resonance, basis, derivative):
    driver = GraniittiDriver()
    dataset = dataset_card(tmp_path, "CMS_2752118/pipi_0p7")
    gencard = load_json_file(tmp_path / "gencard.json", loader=pyjson5.load)
    gencard["SCATTERING"]["RES"] = [resonance]
    write_json_file(tmp_path / "gencard.json", gencard)
    prefix = f"RES|{resonance}:MP:"
    pole = load_json_file(ROOT / f"modeldata/TUNE0/RES/{resonance}.json", loader=pyjson5.load)
    source = "TUNE0"
    magnitude, coupling_phase = prefix + "g[0]", prefix + "g[1]"
    block = driver._resonance_block(pole, "MP")
    if basis in {"g_ls", "helicity"}:
        source = str(tmp_path / "source_tune")
        shutil.copytree(ROOT / "modeldata/TUNE0", source)
        rows = [row.copy() for row in driver._resonance_block(pole, "XP")[basis]]
        if len(rows) > 1:
            scale = max(row[2] for row in rows)
            for index, row in enumerate(rows):
                row[2:4] = [scale / (index + 1), 0.23 * index] if index < 2 else [0.0, 0.0]
        block.update(basis=basis, **{basis: rows})
        block.pop("g")
        write_json_file(Path(source) / f"RES/{resonance}.json", pole)
        magnitude = prefix + angular_row_name(basis, rows[0]) + "@MAG"
        coupling_phase = prefix + angular_row_name(basis, rows[0]) + "@PHASE"
        coupling = rows[0][2]
    else:
        coupling = block["g"][0]
    space = {
        prefix + "phi": {"type": "uniform", "lower": -math.pi, "upper": math.pi},
    }
    if coupling is None:
        space[coupling_phase] = {"type": "uniform", "lower": -math.pi, "upper": math.pi}
    else:
        space[magnitude] = {"type": "uniform", "lower": 0.5 * coupling, "upper": 1.5 * coupling}
    angle = spherical_angle_key(prefix + "polarization.a_Jz", 0)
    spin_phase = raw_phase_key(prefix + "polarization.phi_Jz[2]")
    if frame != "CM":
        space[angle] = {"type": "uniform", "lower": -math.pi, "upper": math.pi}
        space[spin_phase] = {"type": "uniform", "lower": -math.pi, "upper": math.pi}
    auxiliary = {"REGGE|MP_FRAME": frame, "REGGE|DERIVATIVE_FACTOR.MP": derivative,
                 prefix + "basis": basis, prefix + "polarization.mode": "a_Jz"}
    # A pure m=0 state preserves azimuthal symmetry in fixed CM axes
    if frame == "CM":
        auxiliary[prefix + "polarization.a_Jz"] = [[0, 1.0, 0.0]]
    tunesetup = SimpleNamespace(
        param_space=space, aux_param_space=auxiliary,
        datacards=[card(dataset, nevents=8, loopscreen=screened, xsmode="sample",
                        swap_process="MP[RES+CON]<F>", integrator="VEGAS")],
    )
    controls = load_settings(resource("tune/settings/ampfit.json"))["bank"]
    steering = dict(tune_default=source, tunesetup_name=tmp_path.name, data_covariance_mode="diagonal", ampfit=controls)
    initial = driver.initialize(run_name=str(tmp_path / "run"), tunesetup=tunesetup, mc_steer=steering,
                                obs_module="default", cdir=str(ROOT), init_force=True,
                                max_t=600, pickle_dump=False, processes=1, rngseed=9182)
    directory = driver.amplitude_directory(cdir=str(ROOT), run_name=str(tmp_path / "run")) / "0/nominal"
    bank = AmplitudeBank(driver=driver, directory=directory, obs=driver.obs[0], pid=driver.pid[0],
                         cuts=driver.cuts[0], controls=controls)
    bank.validate_source(bank.amplitude(initial), controls["closure_rtol"], directory / "closure.json")
    varied = {**initial, prefix + "phi": 0.37}
    if coupling is None:
        varied[coupling_phase] = 0.23
    else:
        varied[magnitude] = 1.13 * initial[magnitude]
    if frame != "CM":
        varied[angle] = 0.47
        varied[spin_phase] = 0.29
    tune = tmp_path / "direct_tune"
    driver.create_steering_card(param_space=varied | auxiliary, tunename=str(tune), cdir=str(ROOT), tune_default=source)
    run = driver._run_cards(index=0, datacard=tunesetup.datacards[0], mc_steer=steering, tunename="TUNE0", cdir=str(ROOT))[0]
    direct = tmp_path / "direct.json"
    write_json_file(direct, input_card(driver, run, tune))
    output = tmp_path / "direct.bin"
    result = subprocess.run(
        [str(ROOT / "bin/ampfit"), "--input", str(direct), "--events", str(directory / "events.hepmc3"),
         "--output", str(output), "--batch-events", str(controls["batch_events"]),
         "--closure-rtol", str(controls["closure_rtol"])],
        cwd=ROOT, capture_output=True, text=True, timeout=120, check=False,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    expected = np.fromfile(output, dtype=np.complex128).reshape(bank.amplitudes[0].shape)
    np.testing.assert_allclose(bank.coherent(bank.coefficients(varied)).detach().numpy(), expected, rtol=1e-9, atol=1e-12)
    names = sorted(space)
    theta = torch.tensor([varied[name] for name in names], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(
        lambda value: bank.intensity(bank.coefficients(varied | dict(zip(names, value, strict=True)))) / bank.reference,
        (theta,),
    )
