# Tests for direct MP polarization coordinates in icetune
#
# (c) 2026 Mikael Mieskolainen
# Licensed under the MIT License <http://opensource.org/licenses/MIT>.

import json
import math
import shutil
from pathlib import Path

import numpy as np
import pytest
from core.tune.drivers.graniitti import driver as graniitti_driver
from core.tune.drivers.graniitti.tunesetup.domains import mp_spin_sectors
from core.tune.parameters import tools

json5 = pytest.importorskip("pyjson5")
ROOT = Path(__file__).resolve().parents[3]


# Prepare direct production cards for spin coordinate tests
def _prepare_cdir(tmp_path):
    cdir = tmp_path / "graniitti"
    tune = cdir / "modeldata" / "TUNE0"
    tune.parent.mkdir(parents=True)
    shutil.copytree(ROOT / "modeldata" / "TUNE0", tune)
    for name in ("f2_1270", "f4_2300", "f6_2510"):
        path = tune / "RES" / f"{name}.json"
        card = json5.loads(path.read_text())
        block = _block(card)
        block["basis"] = "auto_min_L"
        block["Lambda"] = 1.0
        block["polarization"]["mode"] = "a_Jz"
        j = card["PARAM_RES"]["spinX2"] // 2
        block["polarization"]["a_Jz"] = [[m, float(m == 0), 0.0] for m in range(-j, 1)]
        path.write_text(json.dumps(card))
    return cdir


# Select the single production channel used by these physical resonance tests
def _block(card):
    return next(value for key, value in card["PARAM_RES"]["MODELS"]["MP"].items() if key.startswith("["))


# Read the complex compact coordinates written by the actual driver
def _amplitudes(cdir, tune, name):
    path = cdir / "modeldata" / tune / "RES" / f"{name}.json"
    block = _block(json.loads(path.read_text()))
    return {int(m): mag * np.exp(1j * phase) for m, mag, phase in block["polarization"]["a_Jz"]}


# Check signed coherent states and the physical beam exchange sectors in both geometries
@pytest.mark.parametrize("frame", ["CM", "CS", "HX"])
@pytest.mark.parametrize("geometry", ["projective", "sphere"])
@pytest.mark.parametrize("name,j", [("f2_1270", 2), ("f4_2300", 4), ("f6_2510", 6), ("eta", 0), ("chi_c1", 1), ("rho3_1690", 3), ("eta2_1645", 2)])
def test_direct_spin_coordinates(tmp_path, frame, geometry, name, j):
    cdir = _prepare_cdir(tmp_path)
    driver = graniitti_driver.GraniittiDriver()
    labels = mp_spin_sectors(j, frame)
    base = f"RES|{name}:MP:polarization.a_Jz"
    key = tools.ajzp_angle_key if geometry == "projective" else tools.spherical_angle_key
    angles = [0.6] * (len(labels) - 1)
    if angles:
        angles[-1] = -0.4 if geometry == "projective" else 2.1
    params = {key(base, k): value for k, value in enumerate(angles)}
    params.update({"REGGE|MP_FRAME": frame, f"RES|{name}:MP:polarization.mode": "a_Jz"})
    driver.create_steering_card(param_space=params, tunename="DIRECT", cdir=str(cdir))
    actual = _amplitudes(cdir, "DIRECT", name)
    vector = tools.projective_vector_from_angles(angles) if geometry == "projective" else tools.spherical_vector_from_angles(angles)
    expected = dict(zip(labels, vector, strict=True))
    for m, value in actual.items():
        assert value == pytest.approx(expected.get(abs(m), 0.0) / math.sqrt(1.0 if m == 0 else 2.0))
    assert sum((1 if m == 0 else 2) * abs(value)**2 for m, value in actual.items()) == pytest.approx(1.0)


# Preserve relative phases by extracting and replaying the actual card coordinates
@pytest.mark.parametrize("frame", ["CS", "HX"])
def test_coherent_spin_phase_roundtrip(tmp_path, frame):
    cdir = _prepare_cdir(tmp_path)
    driver = graniitti_driver.GraniittiDriver()
    path = cdir / "modeldata" / "TUNE0" / "RES" / "f4_2300.json"
    card = json.loads(path.read_text())
    rows = [[-4, 0.3, 0.7], [-3, 0.0, 0.0], [-2, 0.4, -0.2], [-1, 0.0, 0.0], [0, math.sqrt(0.5), 0.3]]
    _block(card)["polarization"]["a_Jz"] = rows
    path.write_text(json.dumps(card))
    base = "RES|f4_2300:MP:polarization"
    params = {tools.spherical_angle_key(base + ".a_Jz", k): 0.0 for k in range(2)}
    params.update({tools.raw_phase_key(f"{base}.phi_Jz[{m}]"): 0.0 for m in (2, 4)})
    params["REGGE|MP_FRAME"] = frame
    original = driver.create_steering_card(param_space=params, tunename="PROBE", cdir=str(cdir))
    driver.create_steering_card(param_space=original, tunename="REPLAY", cdir=str(cdir))
    expected = {m: mag * np.exp(1j * phase) for m, mag, phase in rows}
    actual = _amplitudes(cdir, "REPLAY", "f4_2300")
    for m in expected:
        assert actual[m] == pytest.approx(expected[m], abs=1e-12)


# Check derivatives of the same coherent amplitude map used in ampfit and card writing
def test_coherent_spin_autograd():
    import torch

    driver = graniitti_driver.GraniittiDriver()
    name = "f4_2300"
    card = json5.loads((ROOT / "modeldata" / "TUNE0" / "RES" / f"{name}.json").read_text())
    base = f"RES|{name}:MP:polarization"

    # Keep the relative phase derivative active at a vanishing amplitude coordinate
    def evaluate(theta):
        param = {tools.spherical_angle_key(base + ".a_Jz", k): theta[k] for k in range(2)}
        param.update({f"{base}.phi_Jz[{m}]": theta[k + 2] for k, m in enumerate((2, 4))})
        sectors = driver._ajz_sectors(res=name, j=4, geometry="sphere", param=param, json_data=card)
        return torch.view_as_real(torch.stack(list(sectors.values())))

    theta = torch.tensor([1.1, 0.0, 0.4, -0.3], dtype=torch.float64, requires_grad=True)
    assert torch.autograd.gradcheck(evaluate, (theta,))
    assert torch.autograd.gradgradcheck(evaluate, (theta,))


# Reject incomplete directions and production prescriptions with no target polarization
@pytest.mark.parametrize("name", ["rho_770", "f1_1420", "rho3_1690"])
def test_spin_coords_require_coherent_mode(tmp_path, name):
    cdir = _prepare_cdir(tmp_path)
    driver = graniitti_driver.GraniittiDriver()
    path = cdir / "modeldata" / "TUNE0" / "RES" / f"{name}.json"
    source = json5.loads(path.read_text())
    block = _block(source)
    block.update(basis="auto_min_L", g=[1.0, 0.0], Lambda=1.0)
    block["polarization"]["mode"] = "none"
    path.write_text(json.dumps(source))
    params = {tools.ajzp_angle_key(f"RES|{name}:MP:polarization.a_Jz", 0): 0.6}
    with pytest.raises(ValueError, match="polarization.mode a_Jz"):
        driver.create_steering_card(param_space=params, tunename="INVALID", cdir=str(cdir))


# An incomplete normalized direction has no unambiguous physical state
def test_incomplete_spin_direction(tmp_path):
    cdir = _prepare_cdir(tmp_path)
    driver = graniitti_driver.GraniittiDriver()
    params = {tools.ajzp_angle_key("RES|f4_2300:MP:polarization.a_Jz", 0): 0.6}
    with pytest.raises(Exception, match="Incomplete direction"):
        driver.create_steering_card(param_space=params, tunename="INCOMPLETE", cdir=str(cdir))
